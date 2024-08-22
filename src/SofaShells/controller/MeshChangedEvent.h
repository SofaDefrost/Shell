#pragma once

#include <sofa/core/objectmodel/Event.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/Quat.h>
#include <SofaShells/config.h>

namespace sofa::core::objectmodel
{

using namespace sofa::type;

// Sent by MeshInterpolator when the mesh changes
class SOFA_SHELLS_API MeshChangedEvent : public sofa::core::objectmodel::Event
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
} // namespace
