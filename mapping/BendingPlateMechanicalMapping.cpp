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
#include "BendingPlateMechanicalMapping.inl"

#include <sofa/core/componentmodel/behavior/MechanicalMapping.inl>
#include <sofa/core/componentmodel/behavior/MappedModel.h>
#include <sofa/core/componentmodel/behavior/MechanicalState.h>

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::defaulttype;
using namespace core;
using namespace core::componentmodel::behavior;

SOFA_DECL_CLASS(BendingPlateMechanicalMapping)

int BendingPlateMechanicalMappingClass = core::RegisterObject("Mechanical mapping between triangular shell and a cloud of points")
#ifndef SOFA_FLOAT
.add< BendingPlateMechanicalMapping< MechanicalMapping< MechanicalState<Rigid3dTypes>, MechanicalState<Vec3dTypes> > > >()
.add< BendingPlateMechanicalMapping< Mapping< State<Rigid3dTypes>, MappedModel<Vec3dTypes> > > >()

#endif
#ifndef SOFA_DOUBLE
.add< BendingPlateMechanicalMapping< MechanicalMapping< MechanicalState<Rigid3fTypes>, MechanicalState<Vec3fTypes> > > >()
.add< BendingPlateMechanicalMapping< Mapping< State<Rigid3fTypes>, MappedModel<Vec3fTypes> > > >()
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
//.add< BendingPlateMechanicalMapping< MechanicalMapping< MechanicalState<Rigid3fTypes>, MechanicalState<Vec3dTypes> > > >()
//.add< BendingPlateMechanicalMapping< MechanicalMapping< MechanicalState<Rigid3dTypes>, MechanicalState<Vec3fTypes> > > >()
//.add< BendingPlateMechanicalMapping< Mapping< State<Rigid3fTypes>, MappedModel<Vec3dTypes> > > >()
//.add< BendingPlateMechanicalMapping< Mapping< State<Rigid3dTypes>, MappedModel<Vec3fTypes> > > >()
#endif
#endif
;


#ifndef SOFA_FLOAT
template class BendingPlateMechanicalMapping< MechanicalMapping< MechanicalState<Rigid3dTypes>, MechanicalState<Vec3dTypes> > >;
template class BendingPlateMechanicalMapping< Mapping< State<Rigid3dTypes>, MappedModel<Vec3dTypes> > >;
#endif
#ifndef SOFA_DOUBLE
template class BendingPlateMechanicalMapping< MechanicalMapping< MechanicalState<Rigid3fTypes>, MechanicalState<Vec3fTypes> > >;
template class BendingPlateMechanicalMapping< Mapping< State<Rigid3fTypes>, MappedModel<Vec3fTypes> > >;
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
//template class BendingPlateMechanicalMapping< MechanicalMapping< MechanicalState<Rigid3dTypes>, MechanicalState<Vec3fTypes> > >;
//template class BendingPlateMechanicalMapping< MechanicalMapping< MechanicalState<Rigid3fTypes>, MechanicalState<Vec3dTypes> > >;
//template class BendingPlateMechanicalMapping< Mapping< State<Rigid3dTypes>, MappedModel<Vec3fTypes> > >;
//template class BendingPlateMechanicalMapping< Mapping< State<Rigid3fTypes>, MappedModel<Vec3dTypes> > >;
#endif
#endif

// Mech -> Mech

// Mech -> Mapped

// Mech -> ExtMapped

} // namespace mapping

} // namespace component

} // namespace sofa
