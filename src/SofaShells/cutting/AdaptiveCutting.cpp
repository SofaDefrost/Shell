#include "../initPluginShells.h"

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/helper/Factory.inl>
#include <sofa/gui/PickHandler.h>
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

int AdaptiveCuttingOperationReg = sofa::gui::RegisterOperation("AdaptiveCutting")
    .add< AdaptiveCuttingOperation >();

void AdaptiveCuttingOperation::start()
{
    // Add cutting point
    CuttingAdapter *ca = getAdapter();
    if (ca) {
        ca->addCuttingPoint();
    }
    Operation::start(); // TODO: do we need this?
}


void AdaptiveCuttingOperation::endOperation()
{
    // Free the tracked point
    CuttingAdapter *ca = getAdapter();
    if (ca) {
        ca->freeTrackedPoint();
    }
    this->end(); // TODO: do we need this?
}

void AdaptiveCuttingOperation::wait()
{
    // Update the position in the adaptivity component
    if (!pickHandle) return;
    sofa::component::collision::BodyPicked *picked = pickHandle->getLastPicked();
    if (!picked) return;

    CuttingAdapter *ca = getAdapter();
    if (ca) {
        ca->setTrackedPoint(*picked);
    }
}

component::controller::CuttingAdapter* AdaptiveCuttingOperation::getAdapter()
{
    if (!pickHandle) return NULL;

    sofa::component::collision::BodyPicked *picked = pickHandle->getLastPicked();
    if (!picked) return NULL;

    component::controller::CuttingAdapter *ca = NULL;
    if (picked->body) {
        if (!picked->body->getContext()) std::cout << "no context!\n";
        picked->body->getContext()->get(ca);
    } else if (picked->mstate) {
        if (!picked->mstate->getContext()) std::cout << "no context!2\n";
        picked->mstate->getContext()->get(ca);
    }

    return ca;
}

} // namespace gui

} // namespace sofa
