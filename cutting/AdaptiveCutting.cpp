#include "../initPluginShells.h"

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>
//#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/helper/Factory.inl>
#include <sofa/gui/OperationFactory.h>

//#include "AdaptiveCutting.inl"
#include "AdaptiveCutting.h"

namespace sofa
{

namespace component
{

// ----------  Settings  -----------------------------------------------------

namespace configurationsetting
{

SOFA_DECL_CLASS(AdaptiveCuttingSetting)
int AdaptiveCuttingSettingClass = core::RegisterObject("Configuration for adaptive cutting operation.")
    .add< AdaptiveCuttingSetting >()
    //.addAlias("AdaptiveCutting")
;

} // namespace configurationsetting

// ----------  Performer  ----------------------------------------------------

namespace collision
{

#ifndef SOFA_DOUBLE
helper::Creator<
    InteractionPerformer::InteractionPerformerFactory,
    AdaptiveCuttingPerformer<defaulttype::Vec3fTypes> >
        AdaptiveCuttingPerformerVec3fClass("AdaptiveCutting", true);
#endif
#ifndef SOFA_FLOAT
helper::Creator<
    InteractionPerformer::InteractionPerformerFactory,
    AdaptiveCuttingPerformer<defaulttype::Vec3dTypes> >
        AdaptiveCuttingPerformerVec3dClass("AdaptiveCutting", true);
#endif

#ifndef SOFA_DOUBLE
template class SOFA_SHELLS_API  AdaptiveCuttingPerformer<defaulttype::Vec3fTypes>;
//template class SOFA_SHELLS_API  AdaptiveCuttingPerformer<defaulttype::Rigid3fTypes>;
#endif
#ifndef SOFA_FLOAT
template class SOFA_SHELLS_API  AdaptiveCuttingPerformer<defaulttype::Vec3dTypes>;
//template class SOFA_SHELLS_API  AdaptiveCuttingPerformer<defaulttype::Rigid3dTypes>;
#endif

} // namespace collision

} // namespace component


// ----------  Operation  ----------------------------------------------------

namespace gui
{

int AdaptiveCuttingOperationReg = RegisterOperation("AdaptiveCutting")
    .add< AdaptiveCuttingOperation >();


} // namespace gui

} // namespace sofa
