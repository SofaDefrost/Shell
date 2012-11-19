#include <sofa/helper/fixed_array.h>
#include <sofa/gpu/cuda/CudaCommon.h>
#include <sofa/gpu/cuda/CudaMath.h>
//#include <cuda.h>

#if defined(__cplusplus) && CUDA_VERSION < 2000
namespace sofa
{
namespace gpu
{
namespace cuda
{
#endif

extern "C"
{
void Test2DAdapterCuda3f_computeTriangleNormal(unsigned int size, const void* x, void* tri);
void Test2DAdapterCuda3f_functionalGeom(unsigned int size, const void* x, void* tri);
void Test2DAdapterCuda3f_reduceStep(unsigned int size, void* x, void* pt, void* indices);
void Test2DAdapterCuda3f_restoreUnchanged(unsigned int size, void* x, void* pt, void* indices);
void Test2DAdapterCuda3f_smooth(unsigned int size, void* x, void* tri, void* pt, void* indices);
void Test2DAdapterCuda3f_testAcceptable(unsigned int size, void* x, const void* tri, void* pt, void* indices, float tolerance);
///// parallel
void Test2DAdapterCuda3f_prepareGradients(unsigned int size, const void* x, void* tri);
void Test2DAdapterCuda3f_smoothParallel(unsigned int size, void* x, const void* tri, void* pt);
void Test2DAdapterCuda3f_reduceStepP(unsigned int size, void* x, void* pt);
void Test2DAdapterCuda3f_testAcceptableP(unsigned int size, void* x, const void* tri, void* pt, float tolerance);
void Test2DAdapterCuda3f_restoreUnchangedP(unsigned int size, void* x, void* pt);
//#ifdef SOFA_GPU_CUDA_DOUBLE
//void Test2DAdapterCuda3d_computeTriangleNormal(const void* x, const void* n);
//#endif
}// extern "C"

typedef unsigned int Index;

// NOTE: must be equivalent to the Test2DAdapterData::TriangleData
template <class Real>
struct TriangleData {
    Index nodes[3];
    CudaVec3<Real> normal;
    Real functional;
    CudaVec3<Real> gradient[3];
};

// NOTE: must be equivalent to the Test2DAdapterData::PointData
template <class Real>
struct PointData {
    bool bBoundary;

    unsigned int nNeighboursPt;
    const Index *neighboursPt;

    unsigned int nNeighboursTri;
    const Index *neighboursTri;

    bool bAccepted; /// New position has been accepted in current step.
    CudaVec3<Real> oldpos;
    Real oldworst;
    Real newworst;

    Index mintri;
    CudaVec3<Real> grad;
};

//////////////////////
// GPU-side methods //
//////////////////////

template<class Real>
__device__ CudaVec3<Real> computeTriangleNormal(const CudaVec3<Real>* x,
    const Index nodes[3])
{
    CudaVec3<Real> A, B;
    A = x[ nodes[1] ] - x[ nodes[0] ];
    B = x[ nodes[2] ] - x[ nodes[0] ];

    CudaVec3<Real> normal = CudaVec3<Real>::make(0.0, 0.0, 0.0);

    Real An = invnorm(A), Bn = invnorm(B);
    if (An > 1e-20 && Bn > 1e-20) {
        A = A*An;
        B = B*Bn;
        normal = cross(A, B);
        normal = normal * invnorm(normal);
    }

    return normal;
}

template<class Real>
__device__ Real getMinFunc(Index v, const TriangleData<Real>* tri,
    const PointData<Real>* pt)
{
    unsigned int nElem = pt[v].nNeighboursTri;
    Real value = 1.0;
    // TODO: do some unrolling?
    for (Index it=0; it<nElem; it++) {
        if (value > tri[ pt[v].neighboursTri[it] ].functional) {
            value = tri[ pt[v].neighboursTri[it] ].functional;
        }
    }

    return value;
}

template<class Real>
__device__ Real functionalGeom(const Index t,
    const CudaVec3<Real>* x, const TriangleData<Real>* tri)
{
    // TODO: move outside, pass nodes as arguments

    // TODO: we can precompute these, is it worth it?
    CudaVec3<Real> ab = x[ tri[t].nodes[1] ] - x[ tri[t].nodes[0] ];
    CudaVec3<Real> ca = x[ tri[t].nodes[0] ] - x[ tri[t].nodes[2] ];
    CudaVec3<Real> cb = x[ tri[t].nodes[1] ] - x[ tri[t].nodes[2] ];

    // Normalizing factor so that the value is 1 in maximum
    // TODO: does compiler precompute this? NOTE: is float
    Real m = 2 * sqrt(3.0f);

    m *= norm(cross(ca,cb)); // || CA × CB ||
    m /= norm2(ca) + norm2(ab) + norm2(cb);

    // Is triangle inverted?
    CudaVec3<Real> nnew = computeTriangleNormal<Real>(x, tri[t].nodes);
    if (dot(nnew, tri[t].normal) < 0.0) {
        m *= -1.0;
    }

    return m;
}

template<class Real>
__device__ __inline__ Index translateIndexInTriangle(Index index,
    const TriangleData<Real> &tri)
{
    if (tri.nodes[0] == index) return 0;
    if (tri.nodes[1] == index) return 1;
    return 2;
}

//////////////////////
// Kernels          //
//////////////////////

template<class Real>
__global__ void Test2DAdapterCuda3t_computeTriangleNormal_kernel(unsigned int size,
    const CudaVec3<Real>* x, TriangleData<Real>* tri)
{
    int index = umul24(blockIdx.x,BSIZE)+threadIdx.x;
    if (index < size) {
        tri[index].normal = computeTriangleNormal(x, tri[index].nodes);
    }
}

template<class Real>
__global__ void Test2DAdapterCuda3t_functionalGeom_kernel(unsigned int size,
    const CudaVec3<Real>* x, TriangleData<Real>* tri)
{
    int index = umul24(blockIdx.x,BSIZE)+threadIdx.x;
    if (index < size) {
        tri[index].functional = functionalGeom(index, x, tri);
    }
}


// Laplacian smoothing
template<class Real>
__global__ void Test2DAdapterCuda3t_smoothLaplacian_kernel(unsigned int size,
    CudaVec3<Real>* x, const TriangleData<Real>* tri, PointData<Real>* pt,
    const Index *indices)
{
    typedef CudaVec3<Real> Vec3;

    int index = umul24(blockIdx.x,BSIZE)+threadIdx.x;
    if (index < size) {

        Index v = indices[index];

        pt[v].oldpos = x[v];
        pt[v].oldworst = getMinFunc(v, tri, pt);
        pt[v].bAccepted = false;

        // Compute centroid of polygon from 1-ring around the vertex
        Vec3 xnew = Vec3::make(0.0, 0.0, 0.0);
        for (Index ie=0; ie<pt[v].nNeighboursPt; ie++) {
            xnew += x[ pt[v].neighboursPt[ie] ];
        }
        x[v] = xnew / Real(pt[v].nNeighboursPt);
    }
}

// Search for maximum of the functional
template<class Real>
__global__ void Test2DAdapterCuda3t_smoothOptimize_kernel(unsigned int size,
    CudaVec3<Real>* x, const TriangleData<Real>* tri, PointData<Real>* pt,
    const Index *indices)
{
    typedef CudaVec3<Real> Vec3;

    int index = umul24(blockIdx.x,BSIZE)+threadIdx.x;
    if (index < size) {

        Index v = indices[index];
        Vec3 xold = x[v];

        pt[v].oldpos = x[v];
        pt[v].oldworst = getMinFunc(v, tri, pt);
        pt[v].bAccepted = false;

        unsigned int nElem = pt[v].nNeighboursTri;

        // Compute gradients
        // TODO: do it once for all elements!
        Vec3 grad[10]; // TODO: Vec3 grad[nElem];
        if (nElem > 10) nElem = 10; // XXX
        Real delta = 1e-5;

        // NOTE: Constrained to 2D!
        // TODO: can we use shared memory here?
        // -- X
        x[v].x += delta;
        for (Index it=0; it<nElem; it++) {
            Real m = functionalGeom<Real>(pt[v].neighboursTri[it], x, tri);
            grad[it].x = (m - tri[ pt[v].neighboursTri[it] ].functional)/delta;
        }
        // -- Y
        x[v].x = xold.x;
        x[v].y += delta;
        for (Index it=0; it<nElem; it++) {
            Real m = functionalGeom<Real>(pt[v].neighboursTri[it], x, tri);
            grad[it].y = (m - tri[ pt[v].neighboursTri[it] ].functional)/delta;
        }

        // Find smallest functional with non-zero gradient
        Index imin = 0;
        Real fmin = 1.0;
        for (Index it=0; it<nElem; it++) {
            if ((tri[ pt[v].neighboursTri[it] ].functional < fmin) &&
                (norm2(grad[it]) > 1e-15)) {
                fmin = tri[ pt[v].neighboursTri[it] ].functional;
                imin = it;
            }
        }

        Vec3 step = grad[imin];
        // Find out step size
        Real gamma = 0.05;
        //gamma *= step.norm();
        step = step * invnorm(step);

        x[v] = xold + gamma*step;
    }
}

template<class Real>
__global__ void Test2DAdapterCuda3t_reduceStep_kernel(unsigned int size,
    CudaVec3<Real>* x, const PointData<Real>* pt, const Index *indices)
{
    typedef CudaVec3<Real> Vec3;

    int index = umul24(blockIdx.x,BSIZE)+threadIdx.x;
    if (index < size) {

        Index v = indices[index];
        if (!pt[v].bAccepted) {
            // The correct step size is best found empiricaly
            x[v] = pt[v].oldpos + (x[v] - pt[v].oldpos) * Real(2.0/3.0);
            //x[v] = (x[v] + pt[v].oldpos)/2.0;
        }
    }
}

template<class Real>
__global__ void Test2DAdapterCuda3t_restoreUnchanged_kernel(unsigned int size,
    CudaVec3<Real>* x, const PointData<Real>* pt, const Index *indices)
{
    int index = umul24(blockIdx.x,BSIZE)+threadIdx.x;
    if (index < size) {
        Index v = indices[index];
        if (!pt[v].bAccepted) {
            x[v] = pt[v].oldpos;
        }
    }
}

template<class Real>
__global__ void Test2DAdapterCuda3t_testAcceptable_kernel(unsigned int size,
    CudaVec3<Real>* x, const TriangleData<Real>* tri, PointData<Real>* pt, const Index *indices, float tolerance)
{
    int index = umul24(blockIdx.x,BSIZE)+threadIdx.x;
    if (index < size) {

        //// This check is not worth the effort
        //if ((xold - x[v]).norm2() < 1e-8) {
        //    // No change in position
        //    //std::cout << "No change in position for " << v << "\n";
        //    break;
        //}

        Index v = indices[index];
        //if (!pt[v].bAccepted) { // TODO

        // We accept any change that doesn't decrease worst metric for the
        // triangle set.
        Real newworst = getMinFunc(v, tri, pt);
        if (newworst >= (pt[v].oldworst + tolerance)) {
            pt[v].bAccepted = true;
        }
        pt[v].newworst = newworst;

        //}
    }
}

///// parallel

template<class Real>
__global__ void Test2DAdapterCuda3t_prepareGradients_kernel(unsigned int size,
    CudaVec3<Real>* x, TriangleData<Real>* tri)
{
    typedef CudaVec3<Real> Vec3;

    int index = umul24(blockIdx.x,BSIZE)+threadIdx.x;
    if (index < size) {
        // TODO: handle boundary vertices

        for (int i=0; i<3; i++) {
            // For each corner

            Real delta = 1e-5;

            Index v = tri[index].nodes[i];
            Vec3 xold = x[v];

            // NOTE: Constrained to 2D!
            // TODO: can we use shared memory here?
            // -- X
            x[v].x += delta;
            Real m = functionalGeom<Real>(index, x, tri);
            tri[index].gradient[i].x = (m - tri[index].functional)/delta;
            // -- Y
            x[v].x = xold.x;
            x[v].y += delta;
            m = functionalGeom<Real>(index, x, tri);
            tri[index].gradient[i].y = (m - tri[index].functional)/delta;

            x[v].y = xold.y;
        }
    }
}

template<class Real>
__global__ void Test2DAdapterCuda3t_smoothParallel_kernel(unsigned int size,
    CudaVec3<Real>* x, const TriangleData<Real>* tri, PointData<Real>* pt)
{
    typedef CudaVec3<Real> Vec3;

    int index = umul24(blockIdx.x,BSIZE)+threadIdx.x;
    if (index < size) {

        pt[index].oldpos = x[index];
        pt[index].oldworst = getMinFunc(index, tri, pt);
        pt[index].bAccepted = false;

        unsigned int nElem = pt[index].nNeighboursTri;

        // Find smallest functional with non-zero gradient
        Index tmin = Index(-1);
        Real fmin = 1.0;
        Vec3 step = Vec3::make(0.0, 0.0, 0.0);
        if (!pt[index].bBoundary) {
        for (Index it=0; it<nElem; it++) {
            Index t = pt[index].neighboursTri[it];
            Vec3 tg = tri[t].gradient[ translateIndexInTriangle(index, tri[t]) ];
            if ((tri[t].functional < fmin) && (norm2(tg) > 1e-15)) {
                fmin = tri[t].functional;
                tmin = t;
                step = tg;
            }
        }
        }
        pt[index].mintri = tmin;
        pt[index].grad = step;

        // Sync -- we need mintri to be available for all points
        __syncthreads();

        // Consult neighbourhood and make an estimate
        Vec3 ns = Vec3::make(0.0, 0.0, 0.0);
        if (!pt[index].bBoundary) {

        for (Index it=0; it<nElem; it++) {
            Index t = pt[index].neighboursTri[it];

            Index otherp[2];
            Index othert[2];
            for (int v=0, i=0; v<3; v++) {
                if (tri[t].nodes[v] == index) continue;
                otherp[i] = tri[t].nodes[v];
                othert[i] = pt[ otherp[i] ].mintri;
                i++;
            }

            Vec3 tmp = step;
            if (tmin == othert[0]) tmp = tmp / Real(2.0);
            if (tmin == othert[1]) tmp = tmp / Real(2.0);

            if (tmin != othert[0] && tmin != othert[1]) {
                tmp += (pt[ otherp[0] ].grad + pt[ otherp[1] ].grad)/Real(2.0);
            } else if (tmin != othert[0]) {
                tmp += pt[ otherp[0] ].grad;
            } else if (tmin != othert[1]) {
                tmp += pt[ otherp[1] ].grad;
            } 

            ns += tmp;
        }
        ns = ns / Real(nElem);

        }

        Real gamma = 0.005;

        x[index] += gamma*ns;
    }
}

template<class Real>
__global__ void Test2DAdapterCuda3t_reduceStepP_kernel(unsigned int size,
    CudaVec3<Real>* x, const PointData<Real>* pt)
{
    typedef CudaVec3<Real> Vec3;

    int index = umul24(blockIdx.x,BSIZE)+threadIdx.x;
    if (index < size) {

        if (!pt[index].bAccepted) {
            // The correct step size is best found empiricaly
            x[index] = pt[index].oldpos
                + (x[index] - pt[index].oldpos) * Real(2.0/3.0);
            //x[index] = (x[v] + pt[index].oldpos)/2.0;
        }
    }
}

template<class Real>
__global__ void Test2DAdapterCuda3t_testAcceptableP_kernel(unsigned int size,
    CudaVec3<Real>* x, const TriangleData<Real>* tri, PointData<Real>* pt,
    float tolerance)
{
    int index = umul24(blockIdx.x,BSIZE)+threadIdx.x;
    if (index < size) {

        //// This check is not worth the effort
        //if ((xold - x[index]).norm2() < 1e-8) {
        //    // No change in position
        //    //std::cout << "No change in position for " << index << "\n";
        //    break;
        //}

        //if (!pt[index].bAccepted) { // TODO

        // We accept any change that doesn't decrease worst metric for the
        // triangle set.
        Real newworst = getMinFunc(index, tri, pt);
        if (newworst >= (pt[index].oldworst + tolerance)) {
            pt[index].bAccepted = true;
        }
        pt[index].newworst = newworst;

        //}
    }
}

template<class Real>
__global__ void Test2DAdapterCuda3t_restoreUnchangedP_kernel(unsigned int size,
    CudaVec3<Real>* x, const PointData<Real>* pt)
{
    int index = umul24(blockIdx.x,BSIZE)+threadIdx.x;
    if (index < size) {
        if (!pt[index].bAccepted) {
            x[index] = pt[index].oldpos;
        }
    }
}



//////////////////////
// CPU-side methods //
//////////////////////


void Test2DAdapterCuda3f_computeTriangleNormal(unsigned int size, const void* x, void* tri)
{
    dim3 threads(BSIZE,1);
    dim3 grid((size+BSIZE-1)/BSIZE,1);
    Test2DAdapterCuda3t_computeTriangleNormal_kernel<float><<< grid, threads >>>(size, (const CudaVec3<float>*)x, (TriangleData<float>*)tri);
    mycudaDebugError("Test2DAdapterCuda3t_computeTriangleNormal_kernel<float>");
}

void Test2DAdapterCuda3f_functionalGeom(unsigned int size, const void* x, void* tri)
{
    dim3 threads(BSIZE,1);
    dim3 grid((size+BSIZE-1)/BSIZE,1);
    Test2DAdapterCuda3t_functionalGeom_kernel<float><<< grid, threads >>>(size, (const CudaVec3<float>*)x, (TriangleData<float>*)tri);
    mycudaDebugError("Test2DAdapterCuda3t_functionalGeom_kernel<float>");
}

void Test2DAdapterCuda3f_reduceStep(unsigned int size, void* x, void* pt, void* indices)
{
    dim3 threads(BSIZE,1);
    dim3 grid((size+BSIZE-1)/BSIZE,1);
    Test2DAdapterCuda3t_reduceStep_kernel<float><<< grid, threads >>>(size, (CudaVec3<float>*)x, (const PointData<float>*)pt, (const Index*) indices);
    mycudaDebugError("Test2DAdapterCuda3t_reduceStep_kernel<float>");
}

void Test2DAdapterCuda3f_restoreUnchanged(unsigned int size, void* x, void* pt, void* indices)
{
    dim3 threads(BSIZE,1);
    dim3 grid((size+BSIZE-1)/BSIZE,1);
    Test2DAdapterCuda3t_restoreUnchanged_kernel<float><<< grid, threads >>>(size, (CudaVec3<float>*)x, (const PointData<float>*)pt, (const Index*) indices);
    mycudaDebugError("Test2DAdapterCuda3t_restoreUnchanged_kernel<float>");
}

void Test2DAdapterCuda3f_smooth(unsigned int size, void* x, void* tri, void* pt, void* indices)
{
    dim3 threads(BSIZE,1);
    dim3 grid((size+BSIZE-1)/BSIZE,1);
    //Test2DAdapterCuda3t_smoothLaplacian_kernel<float><<< grid, threads >>>(
    //    size, (CudaVec3<float>*)x, (const TriangleData<float>*)tri,
    //    (PointData<float>*)pt, (const Index*) indices);
    //mycudaDebugError("Test2DAdapterCuda3t_smoothLaplacian_kernel<float>");
    Test2DAdapterCuda3t_smoothOptimize_kernel<float><<< grid, threads >>>(
        size, (CudaVec3<float>*)x, (const TriangleData<float>*)tri,
        (PointData<float>*)pt, (const Index*) indices);
    mycudaDebugError("Test2DAdapterCuda3t_smoothOptimize_kernel<float>");
}

void Test2DAdapterCuda3f_testAcceptable(unsigned int size, void* x, const void* tri, void* pt, void* indices, float tolerance)
{
    dim3 threads(BSIZE,1);
    dim3 grid((size+BSIZE-1)/BSIZE,1);
    Test2DAdapterCuda3t_testAcceptable_kernel<float><<< grid, threads >>>(size, (CudaVec3<float>*)x, (const TriangleData<float>*)tri, (PointData<float>*)pt, (const Index*) indices, tolerance);
    mycudaDebugError("Test2DAdapterCuda3t_testAcceptable_kernel<float>");
}


/// Specific to parallel version

void Test2DAdapterCuda3f_prepareGradients(unsigned int size, const void* x,
    void* tri)
{
    dim3 threads(BSIZE,1);
    dim3 grid((size+BSIZE-1)/BSIZE,1);
    Test2DAdapterCuda3t_prepareGradients_kernel<float>
        <<< grid, threads >>>(
            size, (CudaVec3<float>*)x, (TriangleData<float>*)tri);
    mycudaDebugError("Test2DAdapterCuda3t_prepareGradients_kernel<float>");
}

void Test2DAdapterCuda3f_smoothParallel(unsigned int size, void* x,
    const void* tri, void* pt)
{
    dim3 threads(BSIZE,1);
    dim3 grid((size+BSIZE-1)/BSIZE,1);
    Test2DAdapterCuda3t_smoothParallel_kernel<float>
        <<< grid, threads >>>(
            size, (CudaVec3<float>*)x, (const TriangleData<float>*)tri,
            (PointData<float>*)pt);
    mycudaDebugError("Test2DAdapterCuda3t_smoothParallel_kernel<float>");
}

void Test2DAdapterCuda3f_reduceStepP(unsigned int size, void* x, void* pt)
{
    dim3 threads(BSIZE,1);
    dim3 grid((size+BSIZE-1)/BSIZE,1);
    Test2DAdapterCuda3t_reduceStepP_kernel<float>
        <<< grid, threads >>>(
            size, (CudaVec3<float>*)x, (const PointData<float>*)pt);
    mycudaDebugError("Test2DAdapterCuda3t_reduceStepP_kernel<float>");
}

void Test2DAdapterCuda3f_testAcceptableP(unsigned int size, void* x,
    const void* tri, void* pt, float tolerance)
{
    dim3 threads(BSIZE,1);
    dim3 grid((size+BSIZE-1)/BSIZE,1);
    Test2DAdapterCuda3t_testAcceptableP_kernel<float>
        <<< grid, threads >>>(
            size, (CudaVec3<float>*)x, (const TriangleData<float>*)tri,
            (PointData<float>*)pt, tolerance);
    mycudaDebugError("Test2DAdapterCuda3t_testAcceptableP_kernel<float>");
}

void Test2DAdapterCuda3f_restoreUnchangedP(unsigned int size, void* x,
    void* pt)
{
    dim3 threads(BSIZE,1);
    dim3 grid((size+BSIZE-1)/BSIZE,1);
    Test2DAdapterCuda3t_restoreUnchangedP_kernel<float>
        <<< grid, threads >>>(
            size, (CudaVec3<float>*)x, (const PointData<float>*)pt);
    mycudaDebugError("Test2DAdapterCuda3t_restoreUnchangedP_kernel<float>");
}


#if defined(__cplusplus) && CUDA_VERSION < 2000
} // namespace cuda
} // namespace gpu
} // namespace sofa
#endif
