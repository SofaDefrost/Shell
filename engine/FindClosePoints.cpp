#define SOFA_COMPONENT_ENGINE_FINDCLOSEPOINTS_CPP
#include "FindClosePoints.inl"
#include <sofa/core/ObjectFactory.h>
#include "../initPluginShells.h"

namespace sofa
{

namespace component
{

namespace engine
{

int FindClosePointsClass = core::RegisterObject(
    "Finds points whose distance is smaller than defined threshold")

.add< FindClosePoints<defaulttype::Vec3Types> >(true) // default template
.add< FindClosePoints<defaulttype::Vec1Types> >()
.add< FindClosePoints<defaulttype::Vec2Types> >()
.add< FindClosePoints<defaulttype::Rigid2Types> >()
.add< FindClosePoints<defaulttype::Rigid3Types> >()
;

template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec1Types>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec2Types>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec3Types>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Rigid2Types>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Rigid3Types>;


} // namespace constraint

} // namespace component

} // namespace sofa
