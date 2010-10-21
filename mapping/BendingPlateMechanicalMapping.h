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

#include <sofa/core/behavior/MechanicalMapping.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/vector.h>
#include <sofa/helper/gl/GLSLShader.h>
//#include <sofa/core/VisualModel.h>

#include <sofa/component/topology/TriangleSetTopologyContainer.h>

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
using namespace core::topology;

template <class BasicMapping>
class BendingPlateMechanicalMapping : public BasicMapping, public virtual core::objectmodel::BaseObject
// public core::VisualModel,
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(BendingPlateMechanicalMapping,BasicMapping), BasicMapping);
    typedef BasicMapping Inherit;
    typedef typename Inherit::In In;
    typedef typename Inherit::Out Out;
    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename Out::DataTypes OutDataTypes;
//    typedef typename std::map<unsigned int, OutDeriv>::const_iterator OutConstraintIterator;

    typedef typename In::VecCoord InVecCoord;
    typedef typename In::VecDeriv InVecDeriv;
    typedef typename In::Coord InCoord;
    typedef typename In::Deriv InDeriv;
    typedef typename In::DataTypes InDataTypes;

    typedef BaseMeshTopology::Edge	Edge;
    typedef BaseMeshTopology::SeqEdges	SeqEdges;
    typedef BaseMeshTopology::Triangle	Triangle;
    typedef BaseMeshTopology::SeqTriangles SeqTriangles;

    typedef typename Out::Real Real;
    typedef Vec<3, Real> Vec3;

    typedef typename TriangularBendingFEMForceField<InDataTypes>::TriangleInformation TriangleInformation;


    
    BendingPlateMechanicalMapping(In* from, Out* to);

    void init();
    void reinit();
    virtual void draw();
//    virtual void drawVisual();

//    Vec3 direction;
    
    virtual ~BendingPlateMechanicalMapping();
    
    void apply( typename Out::VecCoord& out, const typename In::VecCoord& in );    
    void applyJ( typename Out::VecDeriv& out, const typename In::VecDeriv& in );    
    void applyJT( typename In::VecDeriv& out, const typename Out::VecDeriv& in );
    void applyJT( typename In::MatrixDeriv& out, const typename Out::MatrixDeriv& in );

protected:

        helper::gl::GLSLShader shader;

        BaseMeshTopology* inputTopo;
	BaseMeshTopology* outputTopo;
        
        Data<bool> measureError;
        Data<std::string> nameTargetTopology;

        TriangleSetTopologyContainer* topologyTarget;
        OutVecCoord verticesTarget;
        SeqTriangles trianglesTarget;
        
        helper::vector<Vec3> colourMapping;
        helper::vector<Vec3> coloursPerVertex;
        helper::vector<Real> vectorErrorCoarse;
        helper::vector<Real> vectorErrorTarget;

        // Pointer on the forcefield associated with the in topology
        TriangularBendingFEMForceField<InDataTypes>* triangularBendingForcefield;

        // Pointer on the topological mapping to retrieve the list of edges
        TriangleSubdivisionTopologicalMapping* triangleSubdivisionTopologicalMapping;

        void HSL2RGB(Vec3 &rgb, Real h, Real sl, Real l);
        void MeasureError();
        Real DistanceHausdorff(BaseMeshTopology *topo1, BaseMeshTopology *topo2, helper::vector<Real> &vectorError);
        void ComputeNormals(helper::vector<Vec3> &normals);
        void FindTriangleInNormalDirection(const InVecCoord& highResVertices, const SeqTriangles highRestriangles, const helper::vector<Vec3> &normals);

        // Computes the barycentric coordinates of a vertex within a triangle
        void computeBaryCoefs(Vec3 &baryCoefs, const Vec3 &p, const Vec3 &a, const Vec3 &b, const Vec3 &c);

        Real FindClosestPoints(sofa::helper::vector<unsigned int>& listClosestVertices, const Vec3& point, const OutVecCoord &inVertices);
        Real FindClosestEdges(sofa::helper::vector<unsigned int>& listClosestEdges, const Vec3& point, const OutVecCoord &inVertices, const SeqEdges &inEdges);
        Real FindClosestTriangles(sofa::helper::vector<unsigned int>& listClosestEdges, const Vec3& point, const OutVecCoord &inVertices, const SeqTriangles &inTriangles);

        // Contains the list of base triangles a vertex belongs to
        sofa::helper::vector< sofa::helper::vector<int> > listBaseTriangles;
        // Contains the barycentric coordinates of the same vertex within all base triangles
        sofa::helper::vector< sofa::helper::vector<Vec3> > barycentricCoordinates;
};

} // namespace mapping

} // namespace component

} // namespace sofa

#endif
