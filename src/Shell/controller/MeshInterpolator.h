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

#include <sofa/component/controller/Controller.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/Vec.h>

#include <Shell/config.h>

namespace shell::controller
{

template<class DataTypes>
class MeshInterpolator : public sofa::component::controller::Controller
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(MeshInterpolator,DataTypes), sofa::component::controller::Controller);

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef sofa::type::Vec<3,Real> Vec3;

protected:

    MeshInterpolator();

    virtual ~MeshInterpolator();

public:

    void init() override;
    void reinit() override;

    sofa::Data<Real>              f_startTime;
    sofa::Data<unsigned int>      f_nbSteps;
    sofa::Data<Real>              f_increment;

    sofa::Data<VecCoord>                  f_startPosition;
    sofa::Data<sofa::type::vector<Vec3> > f_startNormals;
    sofa::Data<VecCoord>                  f_endPosition;
    sofa::Data<sofa::type::vector<Vec3> > f_endNormals;

    sofa::Data<VecCoord>                  f_position;
    sofa::Data<sofa::type::vector<Vec3> > f_normals;

    void onEndAnimationStep(const double dt) override;

    Real getInterpolationVar() { return alpha; }

private:

    unsigned int    stepCounter;
    Real            alpha;

    void interpolate();
};

} // namespace
