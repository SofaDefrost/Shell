#ifndef CUDA_TEST2DADAPTER_INL
#define CUDA_TEST2DADAPTER_INL

#include <sofa/gpu/cuda/mycuda.h>
#include <sofa/helper/system/thread/debug.h>
#include "CudaTest2DAdapter.h"
#include "Test2DAdapter.inl"

extern "C"
{
void Test2DAdapterCuda3f_computeTriangleNormal(unsigned int size, const void* x, const void* tri);
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


namespace sofa
{

namespace component
{

namespace controller
{

template<>
void Test2DAdapter< gpu::cuda::CudaVec3fTypes >::onEndAnimationStep(const double /*dt*/)
{
    using namespace sofa::gpu::cuda;

    std::cout << "GPU step\n";

    if ((m_container == NULL) || (m_state == NULL))
        return;

    Index nTriangles = m_container->getNbTriangles();
    if (nTriangles == 0)
        return;

    Data<VecCoord>* datax = m_state->write(sofa::core::VecCoordId::position());
    VecCoord& x = *datax->beginEdit();

    //////////////////////////
    // Initialization
    if (data.colours.size() == 0) {
        // Do the initialization
        // TODO: move this elsewhere

        // Colour the graph of vertices/edges to get the independent sets.
        colourGraph();
        // TODO: use a better container for the sets?

        // Initialize triangle data
        data.triangles.resize(nTriangles);
        for (Index i=0; i<nTriangles; i++) {
            data.triangles[i].nodes = m_container->getTriangle(i);
        }

        // Initialize point data
        Index nPoints = x.size();
        data.points.resize(nPoints);
        data.pointsHost.resize(nPoints);
        for (Index v=0; v<nPoints; v++) {
            // List of point neighbours
            EdgesAroundVertex N1e = m_container->getEdgesAroundVertex(v);
            data.pointsHost[v].neighboursPt.resize(N1e.size()); 

            for (Index ie=0; ie<N1e.size(); ie++) {
                Edge e = m_container->getEdge(N1e[ie]); 
                for (int n=0; n<2; n++) {
                    if (e[n] != v) {
                        data.pointsHost[v].neighboursPt[ie] = e[n];
                    }
                }
            }

            data.points[v].nNeighboursPt = N1e.size();
            data.points[v].neighboursPt = reinterpret_cast<const Index*>
                (data.pointsHost[v].neighboursPt.deviceRead());

            // List of triangle neighbours
            TrianglesAroundVertex N1 = m_container->getTrianglesAroundVertex(v);
            data.pointsHost[v].neighboursTri.resize(N1.size()); 

            for (Index it=0; it<N1.size(); it++) {
                data.pointsHost[v].neighboursTri[it] = N1[it];
            }

            data.points[v].bBoundary = pointInfo.getValue()[v].isBoundary() || pointInfo.getValue()[v].isFixed();
            data.points[v].nNeighboursTri = N1.size();
            data.points[v].neighboursTri= reinterpret_cast<const Index*>
                (data.pointsHost[v].neighboursTri.deviceRead());
        }
    }
    // End of Initialization
    //////////////////////////

    stepCounter++;

    // Array sizes
    //unsigned int szNormals = nTriangles * 3 * sizeof(Real);

    // Compute initial metrics and normals
    // TODO: kernel
    //Real *normals_gpu;
    //mycudaMalloc((void**)&normals_gpu, szNormals);

    //Real *normals_buf;
    //normals_buf = (Real*) malloc(szNormals);
    //mycudaMemcpyDeviceToHost(normals_buf, normals_gpu, szNormals);

    //vector<Vec3> normals(nTriangles);
    //for (unsigned int i=0; i<szNormals; i++) {
    //    normals[i/3][i%3] = normals_buf[i];
    //}

    sofa::helper::system::thread::ctime_t start, stop;
    sofa::helper::system::thread::CTime timer;
    start = timer.getTime();

    //smoothLinear();
    smoothParallel();

    stop = timer.getTime();
    std::cout << "---------- GPU time = " << stop-start << "\n";


    data.triangles.hostRead();

    vector<Real> &functionals = *m_functionals.beginEdit();
    functionals.resize(nTriangles);
    for (Index i=0; i<nTriangles; i++) {
        functionals[i] = data.triangles[i].functional;
    }
    m_functionals.endEdit();

    // Write metrics to file
    std::ofstream of("/tmp/metrics.csv", std::ios::app);
    of << "geom," << stepCounter;
    for (Index i=0; i < m_functionals.getValue().size(); i++) {
        of << "," << m_functionals.getValue()[i];
    }
    of << "\n";
    of.close();

#if 0
    //ngamma = 0;
    //sumgamma = maxgamma = 0.0;
    //mingamma = 1.0;

    //Real maxdelta=0.0;
    //unsigned int moved=0;

    // Evaluate improvement
    Real sum=0.0, sum2=0.0, min = 1.0;
    for (Index i=0; i < nTriangles; i++) {
        if (m_functionals[i] < min) {
            min = m_functionals[i];
        }
        sum += m_functionals[i];
        sum2 += m_functionals[i] * m_functionals[i];
    }
    sum /= nTriangles;
    sum2 = helper::rsqrt(sum2/nTriangles);

    std::cout << stepCounter << "] moved " << moved << " points, max delta=" << helper::rsqrt(maxdelta)
        << " gamma min/avg/max: " << mingamma << "/" << sumgamma/ngamma
        << "/" << maxgamma

        << " Quality min/avg/RMS: " << min << "/" << sum << "/" << sum2

        << "\n";

#endif

    datax->endEdit();
}


template<>
void Test2DAdapter< gpu::cuda::CudaVec3fTypes >::smoothLinear()
{
    Index nTriangles = m_container->getNbTriangles();
    if (nTriangles == 0)
        return;

    Data<VecCoord>* datax = m_state->write(sofa::core::VecCoordId::position());
    VecCoord& x = *datax->beginEdit();

    // Persistent on gpu, we don't need to use this on host
    void *triangles_gpu = data.triangles.deviceWrite();
    void *points_gpu = data.points.deviceWrite();

    // Compute initial normals
    Test2DAdapterCuda3f_computeTriangleNormal(nTriangles, x.deviceRead(),
        triangles_gpu);

    // Compute initial values of the functional
    Test2DAdapterCuda3f_functionalGeom(nTriangles, x.deviceRead(),
        triangles_gpu);

    //data.triangles.hostRead();
    //std::cout << "m: ";
    //for (unsigned int i=0; i<nTriangles; i++) {
    //    std::cout << data.triangles[i].functional << " ";
    //}
    //std::cout << "\n";

    for (unsigned int c=0; c<data.colours.size(); c++) {
        unsigned int nPoints = data.colours[c].size();
        // Compute new position
        Test2DAdapterCuda3f_smooth(nPoints, x.deviceWrite(), triangles_gpu,
            points_gpu, data.colours[c].deviceWrite());
        // Re-evaluate functional
        Test2DAdapterCuda3f_functionalGeom(nTriangles, x.deviceRead(),
            triangles_gpu);
        // Test for acceptance
        Test2DAdapterCuda3f_testAcceptable(nPoints, x.deviceWrite(),
            triangles_gpu, points_gpu, data.colours[c].deviceWrite(),
            m_sigma.getValue());

        for (int i=0; i<10; i++) {
            // Try smaller step
            Test2DAdapterCuda3f_reduceStep(nPoints, x.deviceWrite(),
                points_gpu, data.colours[c].deviceWrite());
            // Re-evaluate functional
            Test2DAdapterCuda3f_functionalGeom(nTriangles, x.deviceRead(),
                triangles_gpu);
            // Test for acceptance
            Test2DAdapterCuda3f_testAcceptable(nPoints, x.deviceWrite(),
                triangles_gpu, points_gpu, data.colours[c].deviceWrite(),
                m_sigma.getValue());
        }

        // Restore positions of unchanged nodes
        Test2DAdapterCuda3f_restoreUnchanged(nPoints, x.deviceWrite(),
            points_gpu, data.colours[c].deviceWrite());

        // Re-evaluate functional
        Test2DAdapterCuda3f_functionalGeom(nTriangles, x.deviceRead(),
            triangles_gpu);

        // TODO: can we track count of accepted solutions to stop sooner (than
        // after 10 loops)? Or, does it even make sense? We expect large sets
        // and then chance that all points will move is small imho. What about
        // including points that cannot move at all?
    }

    //data.points.hostRead();
    //for (unsigned int i=0; i<x.size(); i++) {
    //    std::cout << i << ": " << data.points[i].oldworst << "/" <<
    //        data.points[i].newworst << " " << data.points[i].bAccepted << "\n";
    //}


    // Copy results
    x.hostRead();
    datax->endEdit();
}

template<>
void Test2DAdapter< gpu::cuda::CudaVec3fTypes >::smoothParallel()
{
    Index nTriangles = m_container->getNbTriangles();
    Index nPoints = m_container->getNbPoints();
    if (nTriangles == 0)
        return;

    Data<VecCoord>* datax = m_state->write(sofa::core::VecCoordId::position());
    VecCoord& x = *datax->beginEdit();

    // Persistent on gpu, we don't need to use this on host
    void *triangles_gpu = data.triangles.deviceWrite();
    void *points_gpu = data.points.deviceWrite();

    // Compute initial normals
    Test2DAdapterCuda3f_computeTriangleNormal(nTriangles, x.deviceRead(),
        triangles_gpu);

    // Compute initial values of the functional
    Test2DAdapterCuda3f_functionalGeom(nTriangles, x.deviceRead(),
        triangles_gpu);

    // Prepare gradients
    Test2DAdapterCuda3f_prepareGradients(nTriangles, x.deviceRead(),
        triangles_gpu);

    // Per vertex operations: gradient, assumed step
    Test2DAdapterCuda3f_smoothParallel(nPoints, x.deviceWrite(),
        triangles_gpu, points_gpu);

    // Re-evaluate functional
    Test2DAdapterCuda3f_functionalGeom(nTriangles, x.deviceRead(),
        triangles_gpu);

    // Test for acceptance
    Test2DAdapterCuda3f_testAcceptableP(nPoints, x.deviceWrite(),
        triangles_gpu, points_gpu, m_sigma.getValue());

    data.points.hostRead();
    for (unsigned int i=0; i<x.size(); i++) {
        std::cout << i << ": " << data.points[i].oldworst << "/" <<
            data.points[i].newworst << " " << data.points[i].bAccepted << "\n";
    }

    for (int i=0; i<10; i++) {
        // Try smaller step
        // TODO: do we want something "smarter"?
        Test2DAdapterCuda3f_reduceStepP(nPoints, x.deviceWrite(), points_gpu);
        // Re-evaluate functional
        Test2DAdapterCuda3f_functionalGeom(nTriangles, x.deviceRead(),
            triangles_gpu);
        // Test for acceptance
        Test2DAdapterCuda3f_testAcceptableP(nPoints, x.deviceWrite(),
            triangles_gpu, points_gpu, m_sigma.getValue());
    }

    // Restore positions of unchanged nodes
    // TODO: we can't use the same method as for linear. By reseting the
    // position we may change the functional for neighbouring regions. Also we
    // should somehow keep the position inside the new convex polygon formed by
    // moving the N1-ring instead of simply restoring the original position.
    Test2DAdapterCuda3f_restoreUnchangedP(nPoints, x.deviceWrite(),
        points_gpu);

    // Re-evaluate functional
    Test2DAdapterCuda3f_functionalGeom(nTriangles, x.deviceRead(),
        triangles_gpu);

    // Copy results
    x.hostRead();
    datax->endEdit();
}

template<>
void Test2DAdapter< gpu::cuda::CudaVec3fTypes >::colourGraph()
{
    data.colours.clear();

    int ncolours = 0;
    helper::vector<int> c(m_container->getNbPoints(), -1);

    for (Index v=0; (int)v<m_container->getNbPoints(); v++) {

        if (pointInfo.getValue()[v].isBoundary() || pointInfo.getValue()[v].isFixed()) continue; // Skip boundary vertices

        c[v] = 0;

        EdgesAroundVertex N1e = m_container->getEdgesAroundVertex(v);
        Index ie = 0;
        while (ie < N1e.size()) {
            Edge e = m_container->getEdge(N1e[ie]); 
            Index other = (e[0] == v) ? e[1] : e[0];
            if (c[v] == c[other]) {
                c[v]++;
                ie = 0;
                continue;
            }
            ie++;
        }

        if (c[v] >= ncolours) {
            ncolours = c[v]+1;
            data.colours.resize(ncolours);
        }
        data.colours[ c[v] ].push_back(v);
    }

    std::cout << "GC: " << data.colours.size() << " colours\n";
}

} // namespace controller

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_INL
