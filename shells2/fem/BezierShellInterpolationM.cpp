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


SOFA_DECL_CLASS(BezierShellInterpolationM)

// Register in the Factory
int BezierShellInterpolationMClass = core::RegisterObject("Bezier Shell Interpolation supporting Mechanical Mapping")
#ifndef SOFA_FLOAT
.add< BezierShellInterpolationM<Rigid3dTypes, Vec3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
.add< BezierShellInterpolationM<Rigid3fTypes, Vec3fTypes> >()
#endif
#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
//.add< BezierShellInterpolationM< Rigid3dTypes, Vec3fTypes > >()
//.add< BezierShellInterpolationM< Rigid3fTypes, Vec3dTypes > >()
#endif
#endif
;

#ifndef SOFA_FLOAT
template class SOFA_SHELLS_API BezierShellInterpolationM<Rigid3dTypes, Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_SHELLS_API BezierShellInterpolationM<Rigid3fTypes, Vec3fTypes>;
#endif

//#ifndef SOFA_FLOAT
//#ifndef SOFA_DOUBLE
//template class SOFA_SHELLS_API BezierShellInterpolationM< Rigid3dTypes, Vec3fTypes >;
//template class SOFA_SHELLS_API BezierShellInterpolationM< Rigid3fTypes, Vec3dTypes >;
//#endif
//#endif


} // namespace fem

} // namespace component

} // namespace sofa

