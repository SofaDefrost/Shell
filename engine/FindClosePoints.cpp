#define SOFA_COMPONENT_ENGINE_FINDCLOSEPOINTS_CPP
#include <engine/FindClosePoints.inl>
#include <sofa/core/ObjectFactory.h>
#include "../initPluginShells.h"

namespace sofa
{

namespace component
{

namespace engine
{

SOFA_DECL_CLASS(FindClosePoints)

int FindClosePointsClass = core::RegisterObject(
    "Finds points whose distance is smaller than defined threshold")

#ifdef SOFA_FLOAT
.add< FindClosePoints<defaulttype::Vec3fTypes> >(true) // default template
#else
.add< FindClosePoints<defaulttype::Vec3dTypes> >(true) // default template
# ifndef SOFA_DOUBLE
.add< FindClosePoints<defaulttype::Vec3fTypes> >()
# endif
#endif

#ifndef SOFA_FLOAT
.add< FindClosePoints<defaulttype::Vec1dTypes> >()
.add< FindClosePoints<defaulttype::Vec2dTypes> >()
.add< FindClosePoints<defaulttype::Rigid2dTypes> >()
.add< FindClosePoints<defaulttype::Rigid3dTypes> >()
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
.add< FindClosePoints<defaulttype::Vec1fTypes> >()
.add< FindClosePoints<defaulttype::Vec2fTypes> >()
.add< FindClosePoints<defaulttype::Rigid2fTypes> >()
.add< FindClosePoints<defaulttype::Rigid3fTypes> >()
#endif //SOFA_DOUBLE
;

#ifndef SOFA_FLOAT
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec1dTypes>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec2dTypes>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec3dTypes>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Rigid2dTypes>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Rigid3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec1fTypes>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec2fTypes>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec3fTypes>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Rigid2fTypes>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Rigid3fTypes>;
#endif //SOFA_DOUBLE


} // namespace constraint

} // namespace component

} // namespace sofa
