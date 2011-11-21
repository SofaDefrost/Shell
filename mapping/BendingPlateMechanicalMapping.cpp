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
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>
#include "initPluginShells.h"


namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::defaulttype;
using namespace core;
using namespace core::behavior;

SOFA_DECL_CLASS(BendingPlateMechanicalMapping)

int BendingPlateMechanicalMappingClass = core::RegisterObject("Mechanical mapping between triangular shell (TriangularBendingFEMForceField) and a cloud of points")
#ifndef SOFA_FLOAT
.add< BendingPlateMechanicalMapping< Rigid3dTypes, Vec3dTypes > >()
#endif
#ifndef SOFA_DOUBLE
.add< BendingPlateMechanicalMapping< Rigid3fTypes, Vec3fTypes > >()
#endif
#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
.add< BendingPlateMechanicalMapping< Rigid3dTypes, Vec3fTypes > >()
.add< BendingPlateMechanicalMapping< Rigid3fTypes, Vec3dTypes > >()
#endif
#endif
;


#ifndef SOFA_FLOAT
template class SOFA_SHELLS_API BendingPlateMechanicalMapping< Rigid3dTypes, Vec3dTypes >;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_SHELLS_API BendingPlateMechanicalMapping< Rigid3fTypes, Vec3fTypes >;
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_SHELLS_API BendingPlateMechanicalMapping< Rigid3dTypes, Vec3fTypes >;
template class SOFA_SHELLS_API BendingPlateMechanicalMapping< Rigid3fTypes, Vec3dTypes >;
#endif
#endif


} // namespace mapping

} // namespace component

} // namespace sofa
