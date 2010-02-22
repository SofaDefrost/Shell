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
#define SOFA_COMPONENT_FORCEFIELD_TRIANGULAR_BENDING_FEM_FORCEFIELD_CPP
#include "TriangularBendingFEMForceField.inl"
#include <sofa/core/componentmodel/behavior/ForceField.inl>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/core/ObjectFactory.h>


namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;



SOFA_DECL_CLASS(TriangularBendingFEMForceField)


// Register in the Factory
int TriangularBendingFEMForceFieldClass = core::RegisterObject("Triangular finite elements with bending")
#ifndef SOFA_FLOAT
.add< TriangularBendingFEMForceField<Rigid3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
        // WARNING: THE TEMPLATE BELOW HAS BEEN CHANGED ONLY BECAUSE MACROS SOFA_FLOAT AND SOFA_DOUBLE DO NOT WORK
//.add< TriangularBendingFEMForceField<Rigid3fTypes> >()
#endif
;

#ifndef SOFA_FLOAT
template class SOFA_COMPONENT_FORCEFIELD_API TriangularBendingFEMForceField<Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
    // WARNING: THE TEMPLATE BELOW HAS BEEN CHANGED ONLY BECAUSE MACROS SOFA_FLOAT AND SOFA_DOUBLE DO NOT WORK
//template class SOFA_COMPONENT_FORCEFIELD_API TriangularBendingFEMForceField<Rigid3fTypes>;
#endif


} // namespace forcefield

} // namespace component

} // namespace sofa
