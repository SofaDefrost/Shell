#define SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_CPP
#include "JoinMeshPoints.inl"
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace engine
{

int JoinMeshPointsClass = core::RegisterObject(
    "Join two or more points of the mesh together while updating the topology"
    " accordingly")

.add< JoinMeshPoints<defaulttype::Vec3Types> >(true) // default template

.add< JoinMeshPoints<defaulttype::Vec1Types> >()
.add< JoinMeshPoints<defaulttype::Vec2Types> >()
.add< JoinMeshPoints<defaulttype::Rigid2Types> >()
.add< JoinMeshPoints<defaulttype::Rigid3Types> >()
;

template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec1Types>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec2Types>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec3Types>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Rigid2Types>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Rigid3Types>;

} // namespace constraint

} // namespace component

} // namespace sofa
