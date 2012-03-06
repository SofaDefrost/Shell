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


SOFA_DECL_CLASS(BezierShellInterpolation)

// Register in the Factory
int BezierShellInterpolationClass = core::RegisterObject("Bezier Shell Interpolation")
#ifndef SOFA_FLOAT
.add< BezierShellInterpolation<Rigid3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
.add< BezierShellInterpolation<Rigid3fTypes> >()
#endif
;

#ifndef SOFA_FLOAT
template class SOFA_SHELLS_API BezierShellInterpolation<Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_SHELLS_API BezierShellInterpolation<Rigid3fTypes>;
#endif

} // namespace fem

} // namespace component

} // namespace sofa

