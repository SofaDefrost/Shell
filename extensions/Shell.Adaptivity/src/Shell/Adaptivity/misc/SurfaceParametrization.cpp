//
// Class for parametrizing 3D surface in 2D domain.
//

#include <Shell/Adaptivity/config.h>
#include <Shell/Adaptivity/misc/SurfaceParametrization.inl>

namespace sofa
{

#ifndef SOFA_FLOAT
template class SHELL_ADAPTIVITY_API SurfaceParametrization<double>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SHELL_ADAPTIVITY_API SurfaceParametrization<float>;
#endif //SOFA_DOUBLE

}
