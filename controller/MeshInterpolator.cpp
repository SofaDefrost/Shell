#define SOFA_COMPONENT_CONTROLLER_MESHINTERPOLATION_CPP
#include "MeshInterpolator.inl"
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace controller
{

using namespace sofa::defaulttype;

// Register in the Factory
int MeshInterpolatorClass = core::RegisterObject("Performs linear interpolation between two meshes")
.add< MeshInterpolator<defaulttype::Rigid3Types> >(true) // default template
;

template class SOFA_SHELLS_API MeshInterpolator<defaulttype::Rigid3Types>;

} // namespace controller

} // namespace component

} // namespace sofa
