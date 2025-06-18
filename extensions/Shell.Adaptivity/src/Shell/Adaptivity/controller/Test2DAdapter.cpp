#include <Shell/Adaptivity/config.h>
#include <Shell/Adaptivity/controller/Test2DAdapter.inl>
#include <sofa/core/ObjectFactory.h>


namespace sofa
{

namespace component
{

namespace controller
{

using namespace sofa::defaulttype;

SOFA_DECL_CLASS(Test2DAdapter)

// Register in the Factory
int Test2DAdapterClass = core::RegisterObject("Adaptive mesh improvement component for 2D triangular meshes (for testing)")
#ifdef SOFA_FLOAT
.add< Test2DAdapter<defaulttype::Vec3fTypes> >(true) // default template
#else
.add< Test2DAdapter<defaulttype::Vec3dTypes> >(true) // default template
# ifndef SOFA_DOUBLE
.add< Test2DAdapter<defaulttype::Vec3fTypes> >()
# endif
#endif
;

#ifndef SOFA_FLOAT
template class SHELL_ADAPTIVITY_API Test2DAdapter<defaulttype::Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SHELL_ADAPTIVITY_API Test2DAdapter<defaulttype::Vec3fTypes>;
#endif //SOFA_DOUBLE

} // namespace controller

} // namespace component

} // namespace sofa
