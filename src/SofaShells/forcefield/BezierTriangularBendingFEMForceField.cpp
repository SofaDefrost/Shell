/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#define SOFA_COMPONENT_FORCEFIELD_BEZIER_TRIANGULAR_BENDING_FEM_FORCEFIELD_CPP
#include <SofaShells/forcefield/BezierTriangularBendingFEMForceField.inl>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/defaulttype/RigidTypes.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/ObjectFactory.h>
#include <SofaShells/config.h>


namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;

// Register in the Factory
int BezierTriangularBendingFEMForceFieldClass = core::RegisterObject("Triangular finite elements with bending based on Bézier triangle formulation")
.add< BezierTriangularBendingFEMForceField<Rigid3Types> >()
;

template class SOFA_SHELLS_API BezierTriangularBendingFEMForceField<Rigid3Types>;

} // namespace forcefield

} // namespace component

} // namespace sofa
