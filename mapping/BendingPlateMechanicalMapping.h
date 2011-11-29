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


#include <sofa/core/Mapping.h>



#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/vector.h>


#include <sofa/helper/gl/GLSLShader.h>

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


template <class TIn, class TOut>
class BendingPlateMechanicalMapping : public core::Mapping<TIn, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(BendingPlateMechanicalMapping,TIn,TOut), SOFA_TEMPLATE2(core::Mapping,TIn,TOut));
    typedef core::Mapping<TIn, TOut> Inherit;
    typedef TIn In;
    typedef TOut Out;

    typedef typename In::VecCoord               InVecCoord;
    typedef typename In::VecDeriv               InVecDeriv;
    //typedef typename In::Coord                  InCoord;
    //typedef typename In::Deriv                  InDeriv;
    typedef typename In::MatrixDeriv            InMatrixDeriv;

    typedef typename Out::VecCoord              OutVecCoord;
    typedef typename Out::VecDeriv              OutVecDeriv;
    //typedef typename Out::Coord                 OutCoord;
    //typedef typename Out::Deriv                 OutDeriv;
    typedef typename Out::MatrixDeriv           OutMatrixDeriv;
    typedef typename Out::Real                  Real;

    typedef Vec<3, Real> Vec3;


    //typedef BaseMeshTopology::Edge              Edge;
    typedef BaseMeshTopology::SeqEdges          SeqEdges;
    typedef BaseMeshTopology::Triangle          Triangle;
    typedef BaseMeshTopology::SeqTriangles      SeqTriangles;



    typedef typename TriangularBendingFEMForceField<In>::TriangleInformation TriangleInformation;


    BendingPlateMechanicalMapping(core::State<In>* from, core::State<Out>* to)
    : Inherit(from, to)
    , inputTopo(NULL)
    , outputTopo(NULL)
    , measureError(initData(&measureError, false, "measureError","Error with high resolution mesh"))
    , targetTopology(initLink("targetTopology","Targeted high resolution topology"))
    {
    }

    virtual ~BendingPlateMechanicalMapping()
    {
    }

    void init();
    void reinit();
    virtual void draw(const core::visual::VisualParams* vparams);


    void apply(const core::MechanicalParams *mparams, Data<OutVecCoord>& out, const Data<InVecCoord>& in);
    void applyJ(const core::MechanicalParams *mparams, Data<OutVecDeriv>& out, const Data<InVecDeriv>& in);
    void applyJT(const core::MechanicalParams *mparams, Data<InVecDeriv>& out, const Data<OutVecDeriv>& in);
    void applyJT(const core::ConstraintParams *cparams, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in);

protected:

    BendingPlateMechanicalMapping()
    : Inherit()
    , inputTopo(NULL)
    , outputTopo(NULL)
    , measureError(initData(&measureError, false, "measureError","Error with high resolution mesh"))
    , targetTopology(initLink("targetTopology","Targeted high resolution topology"))
    {
    }

        helper::gl::GLSLShader shader;

        BaseMeshTopology* inputTopo;
        BaseMeshTopology* outputTopo;

        Data<bool> measureError;
        SingleLink<BendingPlateMechanicalMapping<TIn, TOut>,
            sofa::core::topology::BaseMeshTopology,
            BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> targetTopology;

        TriangleSetTopologyContainer* topologyTarget;
        OutVecCoord verticesTarget;
        SeqTriangles trianglesTarget;

        helper::vector<Vec3> colourMapping;
        helper::vector<Vec3> coloursPerVertex;
        helper::vector<Real> vectorErrorCoarse;
        helper::vector<Real> vectorErrorTarget;

        // Pointer on the forcefield associated with the in topology
        TriangularBendingFEMForceField<In>* triangularBendingForcefield;

        // Pointer on the topological mapping to retrieve the list of edges
        // XXX: The edges are no longer there!!!
        //TriangleSubdivisionTopologicalMapping* triangleSubdivisionTopologicalMapping;

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
