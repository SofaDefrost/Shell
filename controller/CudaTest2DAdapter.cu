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
//#ifdef SOFA_GPU_CUDA_DOUBLE
//void Test2DAdapterCuda3d_computeTriangleNormal(const void* x, const void* n);
//#endif
}// extern "C"

typedef unsigned int Index;

// NOTE: should be equivalent to the Test2DAdapterData::TriangleData
template <class Real>
struct TriangleData {
    Index nodes[3];
    CudaVec3<Real> normal;
    Real functional;
};

template <class Real>
struct PointData {
    unsigned int nNeighboursPt;
    const Index *neighboursPt;

    unsigned int nNeighboursTri;
    const Index *neighboursTri;

    bool bAccepted; /// New position has been accepted in current step.
    CudaVec3<Real> oldpos;
    Real oldworst;
    Real newworst;
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
        // TODO: we can precompute these, is it worth it?
        CudaVec3<Real> ab = x[ tri[index].nodes[1] ] - x[ tri[index].nodes[0] ];
        CudaVec3<Real> ca = x[ tri[index].nodes[0] ] - x[ tri[index].nodes[2] ];
        CudaVec3<Real> cb = x[ tri[index].nodes[1] ] - x[ tri[index].nodes[2] ];

        // Normalizing factor so that the value is 1 in maximum
        Real m = 2 * sqrt(3.0f); // TODO: does compiler precompute this? NOTE: is float
        m *= norm(cross(ca,cb)); // || CA × CB ||
        m /= norm2(ca) + norm2(ab) + norm2(cb);

        // Is triangle inverted?
        CudaVec3<Real> nnew = computeTriangleNormal<Real>(x, tri[index].nodes);
        if (dot(nnew, tri[index].normal) < 0.0) {
            m *= -1.0;
        }

        tri[index].functional = m;
    }
}


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

        // Compute centroid of polygon from 1-ring around the vertex
        Vec3 xnew = Vec3::make(0.0, 0.0, 0.0);
        for (Index ie=0; ie<pt[v].nNeighboursPt; ie++) {
            xnew += x[ pt[v].neighboursPt[ie] ];
        }
        x[v] = xnew / Real(pt[v].nNeighboursPt);

        pt[v].bAccepted = false;
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
            x[v] = (x[v] + pt[v].oldpos) * Real(2.0/3.0);
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
    Test2DAdapterCuda3t_smoothLaplacian_kernel<float><<< grid, threads >>>(size, (CudaVec3<float>*)x, (const TriangleData<float>*)tri, (PointData<float>*)pt, (const Index*) indices);
    mycudaDebugError("Test2DAdapterCuda3t_smoothLaplacian_kernel<float>");
}

void Test2DAdapterCuda3f_testAcceptable(unsigned int size, void* x, const void* tri, void* pt, void* indices, float tolerance)
{
    dim3 threads(BSIZE,1);
    dim3 grid((size+BSIZE-1)/BSIZE,1);
    Test2DAdapterCuda3t_testAcceptable_kernel<float><<< grid, threads >>>(size, (CudaVec3<float>*)x, (const TriangleData<float>*)tri, (PointData<float>*)pt, (const Index*) indices, tolerance);
    mycudaDebugError("Test2DAdapterCuda3t_testAcceptable_kernel<float>");
}



#if defined(__cplusplus) && CUDA_VERSION < 2000
} // namespace cuda
} // namespace gpu
} // namespace sofa
#endif
