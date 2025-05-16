//
// Class for parametrizing 3D surface in 2D domain.
//

#include <SofaShells/config.h>
#include <SofaShells/misc/SurfaceParametrization.inl>

namespace sofa
{

#ifndef SOFA_FLOAT
template class SOFA_SHELLS_API SurfaceParametrization<double>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_SHELLS_API SurfaceParametrization<float>;
#endif //SOFA_DOUBLE

}
