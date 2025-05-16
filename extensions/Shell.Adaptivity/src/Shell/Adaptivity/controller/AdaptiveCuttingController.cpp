#include <Shell/config.h>
#include <Shell/Adaptivity/controller/AdaptiveCuttingController.inl>
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
.add< AdaptiveCuttingController<defaulttype3fTypes> >(true) // default template
#else
.add< AdaptiveCuttingController<defaulttype3dTypes> >(true) // default template
# ifndef SOFA_DOUBLE
.add< AdaptiveCuttingController<defaulttype3fTypes> >()
# endif
#endif
;

#ifndef SOFA_FLOAT
template class SHELL_ADAPTIVITY_API AdaptiveCuttingController<defaulttype3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SHELL_ADAPTIVITY_API AdaptiveCuttingController<defaulttype3fTypes>;
#endif //SOFA_DOUBLE

} // namespace controller

} // namespace component

} // namespace sofa
