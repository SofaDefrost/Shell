//
// Class for parametrizing 3D surface in 2D domain.
//

#include "../initPluginShells.h"
#include "SurfaceParametrization.inl"
//#include <sofa/core/ObjectFactory.h>

namespace sofa
{

#ifndef SOFA_FLOAT
template class SOFA_SHELLS_API SurfaceParametrization<double>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_SHELLS_API SurfaceParametrization<float>;
#endif //SOFA_DOUBLE

}
