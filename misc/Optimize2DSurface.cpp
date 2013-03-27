//
// Class for optimizing 2D function serving as a core for smoothing
// triangular networks.
//

#include "../initPluginShells.h"
#include "Optimize2DSurface.inl"

namespace sofa
{

#ifndef SOFA_FLOAT
//template class SOFA_SHELLS_API Optimize2DSurface<double>;
template class SOFA_SHELLS_API Optimize2DSurface<defaulttype::Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
//template class SOFA_SHELLS_API Optimize2DSurface<float>;
template class SOFA_SHELLS_API Optimize2DSurface<defaulttype::Vec3fTypes>;
#endif //SOFA_DOUBLE

}
