#ifndef SOFA_CORE_OBJECTMODEL_MESHCHANGEDEVENT_H
#define SOFA_CORE_OBJECTMODEL_MESHCHANGEDEVENT_H

#include <sofa/core/objectmodel/Event.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/Quat.h>
#include <SofaShells/config.h>

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

    SOFA_EVENT_H( MeshChangedEvent )

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
