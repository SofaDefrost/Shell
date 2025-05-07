#define SOFA_COMPONENT_CONTROLLER_TRIANGLESWITCHEXAMPLE_CPP

#include <Shell/config.h>
#include <Shell/controller/TriangleSwitchExample.inl>
#include <sofa/core/ObjectFactory.h>


namespace sofa
{

namespace component
{

namespace controller
{

using namespace sofa::defaulttype;


// Register in the Factory
int TriangleSwitchExampleClass = core::RegisterObject("Repeatedly switches topology of two neighbouring triangles in predefined interval")
.add< TriangleSwitchExample<defaulttype::Rigid3Types> >(true) // default template
;

template class SOFA_SHELL_API TriangleSwitchExample<defaulttype::Rigid3Types>;

} // namespace controller

} // namespace component

} // namespace sofa
