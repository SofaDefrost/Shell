#ifndef SOFA_CORE_OBJECTMODEL_MESHCHANGEDEVENT_H
#define SOFA_CORE_OBJECTMODEL_MESHCHANGEDEVENT_H

#include <sofa/core/objectmodel/Event.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/defaulttype/Quat.h>
#include <sofa/helper/system/config.h>

namespace sofa
{

namespace core
{

namespace objectmodel
{

using namespace sofa::defaulttype;

// Sent by MeshInterpolator when the mesh changes
class MeshChangedEvent : public sofa::core::objectmodel::Event
{
public:

    MeshChangedEvent(double _alpha) : alpha(_alpha)
    {}

    ~MeshChangedEvent() {}

    double getAlpha() const { return alpha; }

private:

    double alpha; // interpolation parameter
};

} // namespace objectmodel

} // namespace core

} // namespace sofa

#endif
