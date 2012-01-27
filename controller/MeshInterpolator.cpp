#define SOFA_COMPONENT_CONTROLLER_MESHINTERPOLATION_CPP
#include "MeshInterpolator.inl"
#include <sofa/core/ObjectFactory.h>
#include "../initPluginShells.h"

namespace sofa
{

namespace component
{

namespace controller
{

using namespace sofa::defaulttype;

SOFA_DECL_CLASS(MeshInterpolator)

// Register in the Factory
int MeshInterpolatorClass = core::RegisterObject("Performs linear interpolation between two meshes")
#ifdef SOFA_FLOAT
.add< MeshInterpolator<defaulttype::Vec3fTypes> >(true) // default template
#else
.add< MeshInterpolator<defaulttype::Vec3dTypes> >(true) // default template
# ifndef SOFA_DOUBLE
.add< MeshInterpolator<defaulttype::Vec3fTypes> >()
# endif
#endif

#ifndef SOFA_FLOAT
.add< MeshInterpolator<defaulttype::Vec1dTypes> >()
.add< MeshInterpolator<defaulttype::Vec2dTypes> >()
.add< MeshInterpolator<defaulttype::Rigid2dTypes> >()
.add< MeshInterpolator<defaulttype::Rigid3dTypes> >()
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
.add< MeshInterpolator<defaulttype::Vec1fTypes> >()
.add< MeshInterpolator<defaulttype::Vec2fTypes> >()
.add< MeshInterpolator<defaulttype::Rigid2fTypes> >()
.add< MeshInterpolator<defaulttype::Rigid3fTypes> >()
#endif //SOFA_DOUBLE
;

#ifndef SOFA_FLOAT
template class SOFA_SHELLS_API MeshInterpolator<defaulttype::Vec1dTypes>;
template class SOFA_SHELLS_API MeshInterpolator<defaulttype::Vec2dTypes>;
template class SOFA_SHELLS_API MeshInterpolator<defaulttype::Vec3dTypes>;
template class SOFA_SHELLS_API MeshInterpolator<defaulttype::Rigid2dTypes>;
template class SOFA_SHELLS_API MeshInterpolator<defaulttype::Rigid3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_SHELLS_API MeshInterpolator<defaulttype::Vec1fTypes>;
template class SOFA_SHELLS_API MeshInterpolator<defaulttype::Vec2fTypes>;
template class SOFA_SHELLS_API MeshInterpolator<defaulttype::Vec3fTypes>;
template class SOFA_SHELLS_API MeshInterpolator<defaulttype::Rigid2fTypes>;
template class SOFA_SHELLS_API MeshInterpolator<defaulttype::Rigid3fTypes>;
#endif //SOFA_DOUBLE

} // namespace controller

} // namespace component

} // namespace sofa
