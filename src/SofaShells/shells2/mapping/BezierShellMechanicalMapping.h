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
#ifndef SOFA_COMPONENT_MAPPING_BEZIERSHELLMECHANICALMAPPING_H
#define SOFA_COMPONENT_MAPPING_BEZIERSHELLMECHANICALMAPPING_H


#include <sofa/core/Mapping.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/type/vector.h>

#include <sofa/gl/GLSLShader.h>

#include <sofa/linearalgebra/CompressedRowSparseMatrix.h>
#include <sofa/component/topology/container/dynamic/TriangleSetTopologyContainer.h>
#include <sofa/simulation/AnimateBeginEvent.h>

#include <sofa/defaulttype/VecTypes.h>

#include <SofaShells/shells2/fem/BezierShellInterpolationM.h>


namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::type;
using namespace sofa::core::topology;
using namespace core::topology;
using namespace sofa::core::behavior;

template <class TIn, class TOut>
class BezierShellMechanicalMapping : public core::Mapping<TIn, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(BezierShellMechanicalMapping,TIn,TOut), SOFA_TEMPLATE2(core::Mapping,TIn,TOut));
    typedef core::Mapping<TIn, TOut> Inherit;
    typedef TIn In;
    typedef TOut Out;

    typedef typename In::VecCoord               InVecCoord;
    typedef typename In::VecDeriv               InVecDeriv;
    typedef typename In::Coord                  InCoord;
    typedef typename In::Deriv                  InDeriv;
    typedef typename In::MatrixDeriv            InMatrixDeriv;
    typedef typename In::Real                   InReal;

    typedef typename Out::VecCoord              OutVecCoord;
    typedef typename Out::VecDeriv              OutVecDeriv;
    typedef typename Out::Coord                 OutCoord;
    typedef typename Out::Deriv                 OutDeriv;
    typedef typename Out::MatrixDeriv           OutMatrixDeriv;
    typedef typename Out::Real                  OutReal;

    typedef InReal Real;

    typedef Vec<3, Real> Vec3;
    typedef Mat<3, 3, Real> Mat33;

    typedef sofa::type::Quat<Real> Quat;

    typedef typename sofa::component::fem::BezierShellInterpolationM<TIn,TOut>::ShapeFunctions ShapeFunctions;
    typedef typename sofa::component::fem::BezierShellInterpolationM<TIn,TOut>::VecShapeFunctions VecShapeFunctions;

    typedef sofa::Index Index;
    typedef BaseMeshTopology::SeqEdges          SeqEdges;
    typedef BaseMeshTopology::Triangle          Triangle;
    typedef BaseMeshTopology::SeqTriangles      SeqTriangles;

    enum { NIn = sofa::defaulttype::DataTypeInfo<InDeriv>::Size };
    enum { NOut = sofa::defaulttype::DataTypeInfo<OutDeriv>::Size };
    typedef type::Mat<NOut, NIn, Real> MBloc;
    typedef sofa::linearalgebra::CompressedRowSparseMatrix<MBloc> MatrixType;

    BezierShellMechanicalMapping(core::State<In>* from, core::State<Out>* to)
    : Inherit(from, to)
    , inputTopo(NULL)
    , outputTopo(NULL)
    , bsInterpolation(initLink("bsInterpolation","Attached BezierShellInterpolationM object"))
    , measureError(initData(&measureError, false, "measureError","Error with high resolution mesh"))
    , measureStress(initData(&measureStress, false, "measureStress","Tell forcefield to measure stress values at mapped points"))
    , targetTopology(initLink("targetTopology","Targeted high resolution topology"))
    , matrixJ()
    , updateJ(false)
    {
    }

    virtual ~BezierShellMechanicalMapping()
    {
    }

    void init() override;
    void reinit() override;
    //virtual void draw(const core::visual::VisualParams* vparams);


    void apply(const core::MechanicalParams *mparams, Data<OutVecCoord>& out, const Data<InVecCoord>& in) override;
    void applyJ(const core::MechanicalParams *mparams, Data<OutVecDeriv>& out, const Data<InVecDeriv>& in) override;
    void applyJT(const core::MechanicalParams *mparams, Data<InVecDeriv>& out, const Data<OutVecDeriv>& in) override;
    void applyJT(const core::ConstraintParams *cparams, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in) override;


#if 0
    /// For checkJacobian and to hide some deprecation warnings
    const sofa::linearalgebra::BaseMatrix* getJ() { return getJ(NULL); }
    void applyJ(Data<OutVecDeriv>& out, const Data<InVecDeriv>& in)
    { applyJ(NULL, out, in); }
    void applyJ( OutVecDeriv& out, const InVecDeriv& in)
    {
        Data<OutVecDeriv> dout;
        Data<InVecDeriv> din;
        *din.beginEdit() = in;
        din.endEdit();
        dout.beginEdit()->resize(out.size());
        dout.endEdit();
        applyJ(NULL, dout, din);
        out = dout.getValue();
    }
    void applyJT(InVecDeriv& out, const OutVecDeriv& in)
    {
        Data<OutVecDeriv> din;
        Data<InVecDeriv> dout;
        *din.beginEdit() = in;
        din.endEdit();
        dout.beginEdit()->resize(out.size());
        dout.endEdit();
        applyJT(NULL, dout, din);
        out = dout.getValue();
    }
    ///
#endif


    void handleEvent(sofa::core::objectmodel::Event *event) override {
        if (dynamic_cast<simulation::AnimateBeginEvent*>(event))
        {
            //std::cout << "begin\n";
            // We have to update the matrix at every step
            updateJ = true;
        }
    }

protected:

    BezierShellMechanicalMapping()
    : Inherit()
    , inputTopo(NULL)
    , outputTopo(NULL)
    , bsInterpolation(initLink("bsInterpolation","Attached BezierShellInterpolationM object"))
    , measureError(initData(&measureError, false, "measureError","Error with high resolution mesh"))
    , measureStress(initData(&measureStress, false, "measureStress","Tell forcefield to measure stress values at mapped points"))
    , targetTopology(initLink("targetTopology","Targeted high resolution topology"))
    , matrixJ()
    , updateJ(false)
    {
    }

    typedef struct {

        // Nodes of the Bézier triangle
        sofa::type::fixed_array<Vec3, 10> bezierNodesV;

        // Contains the list of poinst connected to this triangle
        sofa::type::vector<int> attachedPoints;

    } TriangleInformation;

        gl::GLSLShader shader;

        BaseMeshTopology* inputTopo;
        BaseMeshTopology* outputTopo;

        SingleLink<BezierShellMechanicalMapping<TIn, TOut>,
            sofa::component::fem::BezierShellInterpolationM<TIn,TOut>,
            BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> bsInterpolation;

        Data<bool> measureError;
        Data<bool> measureStress;
        SingleLink<BezierShellMechanicalMapping<TIn, TOut>,
            sofa::core::topology::BaseMeshTopology,
            BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> targetTopology;

        topology::container::dynamic::TriangleSetTopologyContainer* topologyTarget;

        type::vector<Vec3> colourMapping;
        type::vector<Vec3> coloursPerVertex;
        type::vector<Real> vectorErrorCoarse;
        type::vector<Real> vectorErrorTarget;

        type::vector<TriangleInformation> triangleInfo;
        type::vector<Vec3> projBaryCoords;    // Barycentric coordinates
        VecShapeFunctions projN;                // Precomputed shape functions
        type::vector<Index> projElements;

        std::unique_ptr<MatrixType> matrixJ;
        bool updateJ;

        // Pointer on the topological mapping to retrieve the list of edges
        // XXX: The edges are no longer there!!!
        //TriangleSubdivisionTopologicalMapping* triangleSubdivisionTopologicalMapping;

        void HSL2RGB(Vec3 &rgb, Real h, Real sl, Real l);
        void MeasureError();
        Real DistanceHausdorff(BaseMeshTopology *topo1, BaseMeshTopology *topo2, type::vector<Real> &vectorError);
        void ComputeNormals(type::vector<Vec3> &normals);
        void FindTriangleInNormalDirection(const InVecCoord& highResVertices, const SeqTriangles highRestriangles, const type::vector<Vec3> &normals);

        // Contains the barycentric coordinates of the point within a triangle
        sofa::type::vector<Vec3> barycentricCoordinates;
};

} // namespace mapping

} // namespace component

} // namespace sofa

#endif // #ifndef SOFA_COMPONENT_MAPPING_BEZIERSHELLMECHANICALMAPPING_H
