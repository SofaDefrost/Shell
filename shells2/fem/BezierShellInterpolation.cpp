#include "../../initPluginShells.h"
#include "BezierShellInterpolation.inl"
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
int BezierShellInterpolationClass = core::RegisterObject("Bezier Shell Interpolation")
.add< BezierShellInterpolation<Rigid3Types> >()
;

template class SOFA_SHELLS_API BezierShellInterpolation<Rigid3Types>;

} // namespace fem

} // namespace component

} // namespace sofa

