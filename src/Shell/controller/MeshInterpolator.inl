/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once

#include <Shell/controller/MeshInterpolator.h>
#include <Shell/controller/MeshChangedEvent.h>


namespace shell::controller
{

template<class DataTypes>
MeshInterpolator<DataTypes>::MeshInterpolator()
: f_startTime(initData(&f_startTime, (Real)0, "startTime", "Time when the interpolation starts"))
, f_nbSteps(initData(&f_nbSteps, (unsigned int)10, "nbSteps", "Perform a transition every nbStep steps"))
, f_increment(initData(&f_increment, (Real)0.05, "increment", "How quickly converge to the final positions"))
, f_startPosition(initData(&f_startPosition, "startPosition","Starting positions of the nodes"))
, f_startNormals(initData(&f_startNormals, "startNormals","Starting normals of the nodes"))
, f_endPosition(initData(&f_endPosition, "endPosition","Final positions of the nodes"))
, f_endNormals(initData(&f_endNormals, "endNormals","Final normals of the nodes"))
, f_position(initData(&f_position, "position","Interpolated positions of the nodes"))
, f_normals(initData(&f_normals, "normals","Interpolated normals of the nodes"))
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
        msg_warning() << "Number of start and end nodes has to be the same.";
        if (lenStart > lenEnd) {
            VecCoord &pts = *f_startPosition.beginEdit();
            pts.resize(lenEnd);
            f_startPosition.endEdit();
            lenStart = lenEnd;
        } else {
            VecCoord &pts = *f_endPosition.beginEdit();
            pts.resize(lenStart);
            f_endPosition.endEdit();
        }
    }

    unsigned int lenStartN = f_startNormals.getValue().size();
    unsigned int lenEndN = f_endNormals.getValue().size();

    // Check that the number of nodes is the same
    if (((lenStartN != 0) || (lenEndN != 0)) &&
        ((lenStartN != lenEndN) || (lenStart != lenStartN))) {
        msg_warning() << "Number of start normals, end normals and positions has to be the same.";
        f_startNormals.beginEdit()->clear();
        f_startNormals.endEdit();
        f_endNormals.beginEdit()->clear();
        f_endNormals.endEdit();
    }

    // Check startTime
    if (f_startTime.getValue() < 0.0) {
         msg_warning() << "startTime has to be greater then or equal to 0.";
        *f_startTime.beginEdit() = 0.0;
        f_startTime.endEdit();
    }

    // Check nbSteps
    if (f_nbSteps.getValue() == 0) {
         msg_warning() << "nbSteps has to be nonzero.";
        *f_nbSteps.beginEdit() = 1;
        f_nbSteps.endEdit();
    }

    if ((f_increment.getValue() <= 0.0) || (f_increment.getValue() > 1.0)) {
         msg_warning() << "Increment has to be geater than 0 and "
            "smaller than or equall to 1.";
        *f_increment.beginEdit() = 0.05;
        f_increment.endEdit();
    }

    // XXX: does it make sense to reinit also internal state?
    stepCounter = 0;
    alpha = 0;

    // Start with starting point
    *f_position.beginEdit() = f_startPosition.getValue();
    f_position.endEdit();

    // Start with starting normals
    *f_normals.beginEdit() = f_startNormals.getValue();
    f_normals.endEdit();

}


template<class DataTypes>
void MeshInterpolator<DataTypes>::onEndAnimationStep(const double /*dt*/)
{
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
    shell::objectmodel::MeshChangedEvent mcEvent(alpha);
    this->getContext()->propagateEvent(sofa::core::execparams::defaultInstance(), &mcEvent);
}

template<class DataTypes>
void MeshInterpolator<DataTypes>::interpolate()
{
    const VecCoord &startPt = f_startPosition.getValue();
    const VecCoord &endPt = f_endPosition.getValue();

    const sofa::type::vector<Vec3> &startNorm = f_startNormals.getValue();
    const sofa::type::vector<Vec3> &endNorm = f_endNormals.getValue();

    VecCoord &pt = *f_position.beginEdit();
    sofa::type::vector<Vec3> &norm = *f_normals.beginEdit();

    for (unsigned int i=0; i<startPt.size(); i++) {

        pt[i].getCenter() =
            startPt[i].getCenter() * (1-alpha) +
            endPt[i].getCenter() * alpha;

        pt[i].getOrientation().slerp(
            startPt[i].getOrientation(),
            endPt[i].getOrientation(),
            alpha, false);
    }

    if (startNorm.size() > 0) {
        for (unsigned int i=0; i<startNorm.size(); i++) {
            norm[i] = startNorm[i] * (1.0-alpha) + endNorm[i] * alpha;
        }
    }

    f_position.endEdit();
    f_normals.endEdit();
}

} // namespace
