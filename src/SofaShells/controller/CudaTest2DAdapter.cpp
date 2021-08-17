
#include <SofaShells/config.h>
#include <sofa/gpu/cuda/CudaTypes.h>
#include <SofaShells/controller/CudaTest2DAdapter.inl>
#include <sofa/core/ObjectFactory.h>


namespace sofa
{

namespace component
{

namespace controller
{

template class Test2DAdapter<gpu::cuda::CudaVec3fTypes>;
//#ifdef SOFA_GPU_CUDA_DOUBLE
//template class Test2DAdapter<gpu::cuda::CudaVec3dTypes>;
//#endif // SOFA_GPU_CUDA_DOUBLE

} // namespace controller

} // namespace component

namespace gpu
{

namespace cuda
{

SOFA_DECL_CLASS(CudaTest2DAdapter)

// Register in the Factory
int CudaTest2DAdapterClass = core::RegisterObject("Adaptive mesh improvement component for 2D triangular meshes (for testing) -- the CUDA version")
.add< component::controller::Test2DAdapter<gpu::cuda::CudaVec3fTypes> >()
//#ifdef SOFA_GPU_CUDA_DOUBLE
//.add< component::controller::Test2DAdapter<gpu::cuda::CudaVec3dTypes> >()
//#endif // SOFA_GPU_CUDA_DOUBLE
;

} // namespace cuda

} // namespace gpu

} // namespace sofa
