#ifndef SOFA_COMPONENT_FORCEFIELD_MESHINTERPOLATION_INL
#define SOFA_COMPONENT_FORCEFIELD_MESHINTERPOLATION_INL

#include "MeshInterpolator.h"
#include "MeshChangedEvent.h"


namespace sofa
{

namespace component
{

namespace controller
{

template<class DataTypes>
MeshInterpolator<DataTypes>::MeshInterpolator()
: f_startTime(initData(&f_startTime, (Real)0, "startTime", "Time when the interpolation starts"))
, f_nbSteps(initData(&f_nbSteps, (unsigned int)10, "nbSteps", "Perform a transition every nbStep steps"))
, f_increment(initData(&f_increment, (Real)0.05, "increment", "How quickly converge to the final positions"))
, f_startPosition(initData(&f_startPosition, "startPosition","Starting positions of the nodes"))
, f_endPosition(initData(&f_endPosition, "endPosition","Final positions of the nodes"))
, f_position(initData(&f_position, "position","Interpolated positions of the nodes"))
, stepCounter(0)
{
}

template<class DataTypes>
MeshInterpolator<DataTypes>::~MeshInterpolator()
{
}

template<class DataTypes>
void MeshInterpolator<DataTypes>::init()
{
    reinit();
}

template<class DataTypes>
void MeshInterpolator<DataTypes>::reinit()
{
    *this->f_listening.beginEdit() = true;
    this->f_listening.endEdit();

    unsigned int lenStart = f_startPosition.getValue().size();
    unsigned int lenEnd = f_endPosition.getValue().size();

    // Check that the number of nodes is the same
    if (lenStart != lenEnd) {
        serr << "Number of start end end nodes has to be the same" << sendl;
        if (lenStart > lenEnd) {
            VecCoord &pts = *f_startPosition.beginEdit();
            pts.resize(lenEnd);
            f_startPosition.endEdit();
        } else {
            VecCoord &pts = *f_endPosition.beginEdit();
            pts.resize(lenStart);
            f_endPosition.endEdit();
        }
    }

    // Check startTime
    if (f_startTime.getValue() < 0.0) {
        serr << "startTime has to be greater then or equal to 0" << sendl;
        *f_startTime.beginEdit() = 0.0;
        f_startTime.endEdit();
    }

    // Check nbSteps
    if (f_nbSteps.getValue() == 0) {
        serr << "nbSteps has to be nonzero" << sendl;
        *f_nbSteps.beginEdit() = 1;
        f_nbSteps.endEdit();
    }

    if ((f_increment.getValue() <= 0.0) || (f_increment.getValue() > 1.0)) {
        serr << "Increment has to be geater than 0 and "
            "smaller than or equall to 1" << sendl;
        *f_increment.beginEdit() = 0.05;
        f_increment.endEdit();
    }

    // XXX: does it make sense to reinit also internal state?
    stepCounter = 0;
    alpha = 0;

    // Start with starting point
    *f_position.beginEdit() = f_startPosition.getValue();
    f_position.endEdit();

}


template<class DataTypes>
void MeshInterpolator<DataTypes>::onEndAnimationStep(const double /*dt*/)
{
    //std::cout << "alpha=" << alpha <<
    //    " time=" << getContext()->getTime() <<
    //    " step=" << stepCounter << "\n";

    if (alpha >= 1.0)
        return; // Nothing more to do

    if (getContext()->getTime() < f_startTime.getValue())
        return; // Not yet ...

    stepCounter++;

    if (stepCounter < f_nbSteps.getValue())
        return; // Stay still for a few steps

    stepCounter = 0;

    // Increase the linear factor
    alpha += f_increment.getValue();

    // Update positions
    interpolate();

    // Send signal about change
    core::objectmodel::MeshChangedEvent mcEvent(alpha);
    this->getContext()->propagateEvent(sofa::core::ExecParams::defaultInstance(), &mcEvent);
}

template<class DataTypes>
void MeshInterpolator<DataTypes>::interpolate()
{
    const VecCoord &startPt = f_startPosition.getValue();
    const VecCoord &endPt = f_endPosition.getValue();

    VecCoord &pt = *f_position.beginEdit();

    for (unsigned int i=0; i<startPt.size(); i++) {
        pt[i] = startPt[i] * (1.0-alpha) + endPt[i] * alpha;
    }
}

#if 0

I do not know how to switch between these two implementations cleanly

template<class DataTypes>
void MeshInterpolator<DataTypes>::interpolate()
{
    const VecCoord &startPt = f_startPosition.getValue();
    const VecCoord &endPt = f_endPosition.getValue();

    VecCoord &pt = *f_position.beginEdit();

    for (unsigned int i=0; i<startPt.size(); i++) {

        pt[i].getCenter() =
            startPt[i].getCenter() * (1-alpha) +
            endPt[i].getCenter() * alpha;

        pt[i].getOrientation().slerp(
            startPt[i].getOrientation(),
            endPt[i].getOrientation(),
            alpha, false);
    }
}
#endif

} // namespace controller

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_MESHINTERPOLATION_INL
