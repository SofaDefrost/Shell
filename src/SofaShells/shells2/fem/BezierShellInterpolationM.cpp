#include "../../initPluginShells.h"
#include "BezierShellInterpolationM.inl"
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>


namespace sofa
{

namespace component
{

namespace fem
{

using namespace sofa::defaulttype;

// Register in the Factory
int BezierShellInterpolationMClass = core::RegisterObject("Bezier Shell Interpolation supporting Mechanical Mapping")
.add< BezierShellInterpolationM<Rigid3Types, Vec3Types> >()
;

template class SOFA_SHELLS_API BezierShellInterpolationM<Rigid3Types, Vec3Types>;

} // namespace fem

} // namespace component

} // namespace sofa

