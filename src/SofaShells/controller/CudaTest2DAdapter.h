#ifndef CUDA_TEST2DADAPTER_H
#define CUDA_TEST2DADAPTER_H

#include <sofa/gpu/cuda/CudaTypes.h>
#include <SofaShells/controller/Test2DAdapter.h>

#include <sofa/helper/fixed_array.h>
#include <sofa/helper/vector.h>


namespace sofa
{

namespace gpu
{

namespace cuda
{

template<class DataTypes>
class CudaKernelsTest2DAdapter;

} // namespace cuda

} // namespace gpu

namespace component
{

namespace controller
{

/**
 * @brief Internal data for Test2DAdapter specialized for Cuda version.
 */
template<>
class Test2DAdapterData< gpu::cuda::CudaVec3fTypes >
{
public:
    typedef gpu::cuda::CudaVec3fTypes DataTypes;
    typedef Test2DAdapter<DataTypes> Main;
    typedef DataTypes::Coord Coord;
    typedef DataTypes::Deriv Deriv;
    typedef Coord::value_type Real;

    typedef sofa::component::topology::TriangleSetTopologyContainer::TriangleID     Index;

    // Graph colours (of independent sets)
    helper::vector< sofa::gpu::cuda::CudaVector<Index> > colours;

    struct TriangleData {
        helper::fixed_array<Index,3> nodes;
        Coord normal;
        Real functional;
        helper::fixed_array<Deriv,3> gradient;
    };

    struct PointData {

        bool bBoundary;   /// Is the pont on a border?

        unsigned int nNeighboursPt;     /// Number of points in N1-ring.
        const Index *neighboursPt;      /// List of points in N1-ring.

        unsigned int nNeighboursTri;    /// Number of triangles in N1-ring.
        const Index *neighboursTri;     /// List of triangles in N1-ring.

        bool bAccepted; /// New position has been accepted in current step.
        Coord oldpos;   /// Temporary holder for original point position.
        Real oldworst;
        Real newworst;

        Index mintri;           // Triangle from N1-ring with smallest functional
        Deriv grad;    // Gradient for mintri
    };

    struct PointDataHost {
        sofa::gpu::cuda::CudaVector<Index> neighboursPt;
        sofa::gpu::cuda::CudaVector<Index> neighboursTri;
    };

    sofa::gpu::cuda::CudaVector<TriangleData> triangles;
    sofa::gpu::cuda::CudaVector<PointData> points;
    sofa::gpu::cuda::CudaVector<PointDataHost> pointsHost;

};

template<>
class Test2DAdapterData< gpu::cuda::CudaVec3fTypes >;

template<>
void Test2DAdapter< gpu::cuda::CudaVec3fTypes >::onEndAnimationStep(const double dt);

template<>
void Test2DAdapter< gpu::cuda::CudaVec3fTypes >::smoothLinear();

template<>
void Test2DAdapter< gpu::cuda::CudaVec3fTypes >::smoothParallel();

template<>
void Test2DAdapter< gpu::cuda::CudaVec3fTypes >::colourGraph();


} // namespace controller

} // namespace component

} // namespace sofa

#endif // CUDA_TEST2DADAPTER_H
