#define SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_CPP
#include "JoinMeshPoints.inl"
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace engine
{

SOFA_DECL_CLASS(JoinMeshPoints)

int JoinMeshPointsClass = core::RegisterObject(
    "Join two or more points of the mesh together while updating the topology"
    " accordingly")

#ifdef SOFA_FLOAT
.add< JoinMeshPoints<defaulttype::Vec3fTypes> >(true) // default template
#else
.add< JoinMeshPoints<defaulttype::Vec3dTypes> >(true) // default template
# ifndef SOFA_DOUBLE
.add< JoinMeshPoints<defaulttype::Vec3fTypes> >()
# endif
#endif

#ifndef SOFA_FLOAT
.add< JoinMeshPoints<defaulttype::Vec1dTypes> >()
.add< JoinMeshPoints<defaulttype::Vec2dTypes> >()
.add< JoinMeshPoints<defaulttype::Rigid2dTypes> >()
.add< JoinMeshPoints<defaulttype::Rigid3dTypes> >()
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
.add< JoinMeshPoints<defaulttype::Vec1fTypes> >()
.add< JoinMeshPoints<defaulttype::Vec2fTypes> >()
.add< JoinMeshPoints<defaulttype::Rigid2fTypes> >()
.add< JoinMeshPoints<defaulttype::Rigid3fTypes> >()
#endif //SOFA_DOUBLE
;

#ifndef SOFA_FLOAT
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec1dTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec2dTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec3dTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Rigid2dTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Rigid3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec1fTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec2fTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec3fTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Rigid2fTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Rigid3fTypes>;
#endif //SOFA_DOUBLE


} // namespace constraint

} // namespace component

} // namespace sofa
