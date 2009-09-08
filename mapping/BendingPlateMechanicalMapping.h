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
#ifndef SOFA_COMPONENT_MAPPING_BENDINGPLATEMAPPING_INL
#define SOFA_COMPONENT_MAPPING_BENDINGPLATEMAPPING_INL

#include <sofa/core/componentmodel/behavior/MechanicalMapping.h>
#include <sofa/core/componentmodel/behavior/MechanicalState.h>
#include <sofa/helper/vector.h>

#include "../forcefield/TriangularBendingFEMForceField.h"
#include <sofa/component/topology/TriangleSubdivisionTopologicalMapping.h>

#include <sofa/defaulttype/VecTypes.h>

#include <sofa/helper/system/thread/CTime.h>


namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::defaulttype;
using namespace sofa::component::forcefield;
using namespace sofa::component::topology;
using namespace sofa::helper::system::thread;

template <class BasicMapping>
class BendingPlateMechanicalMapping : public BasicMapping, public virtual core::objectmodel::BaseObject
{
public:
    typedef BasicMapping Inherit;
    typedef typename Inherit::In In;
    typedef typename Inherit::Out Out;
    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;   
    typedef typename std::map<unsigned int, OutDeriv>::const_iterator OutConstraintIterator;        

    typedef typename In::VecCoord InVecCoord;
    typedef typename In::VecDeriv InVecDeriv;
    typedef typename In::Coord InCoord;
    typedef typename In::Deriv InDeriv;
    typedef typename In::DataTypes InDataTypes;

    typedef core::componentmodel::topology::BaseMeshTopology::Edge	Edge;
    typedef core::componentmodel::topology::BaseMeshTopology::SeqEdges	SeqEdges;
    typedef core::componentmodel::topology::BaseMeshTopology::Triangle	Triangle;
    typedef core::componentmodel::topology::BaseMeshTopology::SeqTriangles SeqTriangles;

    typedef typename Out::Real Real;
    typedef Vec<3, Real> Vec3;

    typedef typename TriangularBendingFEMForceField<InDataTypes>::TriangleInformation TriangleInformation;


    
    BendingPlateMechanicalMapping(In* from, Out* to);

    void init();
    void reinit();
    
    virtual ~BendingPlateMechanicalMapping();
    
    void apply( typename Out::VecCoord& out, const typename In::VecCoord& in );    
    void applyJ( typename Out::VecDeriv& out, const typename In::VecDeriv& in );    
    void applyJT( typename In::VecDeriv& out, const typename Out::VecDeriv& in );
    void applyJT( typename In::VecConst& out, const typename Out::VecConst& in );

protected:
	core::componentmodel::topology::BaseMeshTopology* inputTopo;
	core::componentmodel::topology::BaseMeshTopology* outputTopo;

        // Pointer on the forcefield associated with the in topology
        TriangularBendingFEMForceField<InDataTypes>* triangularBendingForcefield;

        // Computes the barycentric coordinates of a vertex within a triangle
        void computeBaryCoefs(Vec3 &baryCoefs, const Vec3 &p, const Vec3 &a, const Vec3 &b, const Vec3 &c);

        Real FindClosestPoints(const Vec3& point1, sofa::helper::vector<unsigned int>& listClosestVertices);
        Real FindClosestEdges(const Vec3& point, sofa::helper::vector<unsigned int>& listClosestEdges);
        Real FindClosestTriangles(const Vec3& point, sofa::helper::vector<unsigned int>& listClosestTriangles);

        // Contains the list of base triangles a vertex belongs to
        sofa::helper::vector< sofa::helper::vector<int> > listBaseTriangles;
        // Contains the barycentric coordinates of the same vertex within all base triangles
        sofa::helper::vector< sofa::helper::vector<Vec3> > barycentricCoordinates;
        // Coefficients ci for each triangle
        sofa::helper::vector< Vec<9, Real> > listCoeffs;
        // List of non-null indices within the displacemente vector u
        int* nonNullIndices;

        vector<Quat> previousOrientation;

};

} // namespace mapping

} // namespace component

} // namespace sofa

#endif
