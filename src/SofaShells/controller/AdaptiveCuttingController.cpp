#include <SofaShells/config.h>
#include <SofaShells/cutting/AdaptiveCuttingController.inl>
#include <sofa/core/ObjectFactory.h>


namespace sofa
{

namespace component
{

namespace controller
{

using namespace sofa::defaulttype;

SOFA_DECL_CLASS(AdaptiveCuttingController)

// Register in the Factory
int AdaptiveCuttingControllerClass = core::RegisterObject(
    "Controller that handles the cutting method based on mesh adaptivity.")
#ifdef SOFA_FLOAT
.add< AdaptiveCuttingController<type::Vec3fTypes> >(true) // default template
#else
.add< AdaptiveCuttingController<type::Vec3dTypes> >(true) // default template
# ifndef SOFA_DOUBLE
.add< AdaptiveCuttingController<type::Vec3fTypes> >()
# endif
#endif
;

#ifndef SOFA_FLOAT
template class SOFA_SHELLS_API AdaptiveCuttingController<type::Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_SHELLS_API AdaptiveCuttingController<type::Vec3fTypes>;
#endif //SOFA_DOUBLE

} // namespace controller

} // namespace component

} // namespace sofa
