//
// Class for optimizing 2D function serving as a core for smoothing
// triangular networks.
//

#include <Shell/Adaptivity/config.h>
#include <Shell/Adaptivity/misc/Optimize2DSurface.inl>

namespace sofa
{

#ifndef SOFA_FLOAT
//template class SHELL_ADAPTIVITY_API Optimize2DSurface<double>;
template class SHELL_ADAPTIVITY_API Optimize2DSurface<defaulttype::Vec3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
//template class SHELL_ADAPTIVITY_API Optimize2DSurface<float>;
template class SHELL_ADAPTIVITY_API Optimize2DSurface<defaulttype3fTypes>;
#endif //SOFA_DOUBLE

}
