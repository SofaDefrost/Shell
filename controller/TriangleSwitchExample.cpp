#define SOFA_COMPONENT_CONTROLLER_TRIANGLESWITCHEXAMPLE_CPP

#include "../initPluginShells.h"
#include "TriangleSwitchExample.inl"
#include <sofa/core/ObjectFactory.h>


namespace sofa
{

namespace component
{

namespace controller
{

using namespace sofa::defaulttype;

SOFA_DECL_CLASS(TriangleSwitchExample)

// Register in the Factory
int TriangleSwitchExampleClass = core::RegisterObject("Repeatedly switches topology of two neighbouring triangles in predefined interval")
#ifdef SOFA_FLOAT
//.add< TriangleSwitchExample<defaulttype::Vec3fTypes> >(true) // default template
.add< TriangleSwitchExample<defaulttype::Rigid3fTypes> >(true) // default template
#else
//.add< TriangleSwitchExample<defaulttype::Vec3dTypes> >(true) // default template
.add< TriangleSwitchExample<defaulttype::Rigid3dTypes> >(true) // default template
# ifndef SOFA_DOUBLE
//.add< TriangleSwitchExample<defaulttype::Vec3fTypes> >()
.add< TriangleSwitchExample<defaulttype::Rigid3fTypes> >()
# endif
#endif

#ifndef SOFA_FLOAT
//.add< TriangleSwitchExample<defaulttype::Vec1dTypes> >()
//.add< TriangleSwitchExample<defaulttype::Vec2dTypes> >()
//.add< TriangleSwitchExample<defaulttype::Vec3dTypes> >()
//.add< TriangleSwitchExample<defaulttype::Rigid2dTypes> >()
//.add< TriangleSwitchExample<defaulttype::Rigid3dTypes> >()
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
//.add< TriangleSwitchExample<defaulttype::Vec1fTypes> >()
//.add< TriangleSwitchExample<defaulttype::Vec2fTypes> >()
//.add< TriangleSwitchExample<defaulttype::Vec3fTypes> >()
//.add< TriangleSwitchExample<defaulttype::Rigid2fTypes> >()
//.add< TriangleSwitchExample<defaulttype::Rigid3fTypes> >()
#endif //SOFA_DOUBLE
;

#ifndef SOFA_FLOAT
//template class SOFA_SHELLS_API TriangleSwitchExample<defaulttype::Vec1dTypes>;
//template class SOFA_SHELLS_API TriangleSwitchExample<defaulttype::Vec2dTypes>;
//template class SOFA_SHELLS_API TriangleSwitchExample<defaulttype::Vec3dTypes>;
//template class SOFA_SHELLS_API TriangleSwitchExample<defaulttype::Rigid2dTypes>;
template class SOFA_SHELLS_API TriangleSwitchExample<defaulttype::Rigid3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
//template class SOFA_SHELLS_API TriangleSwitchExample<defaulttype::Vec1fTypes>;
//template class SOFA_SHELLS_API TriangleSwitchExample<defaulttype::Vec2fTypes>;
//template class SOFA_SHELLS_API TriangleSwitchExample<defaulttype::Vec3fTypes>;
//template class SOFA_SHELLS_API TriangleSwitchExample<defaulttype::Rigid2fTypes>;
template class SOFA_SHELLS_API TriangleSwitchExample<defaulttype::Rigid3fTypes>;
#endif //SOFA_DOUBLE

} // namespace controller

} // namespace component

} // namespace sofa
