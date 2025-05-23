//
// Interpolation of Bézier triangles
//
// Author: Tomáš Golembiovský
//
// Copyright:
//
#ifndef SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_H
#define SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_H

#include <Shell/config.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/component/mapping/linear/Mesh2PointTopologicalMapping.h>

#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/topology/TopologyData.h>
#include <sofa/simulation/AnimateBeginEvent.h>

#include <sofa/type/vector.h>
#include <sofa/type/Vec.h>

#include <sofa/core/objectmodel/BaseObject.h>


namespace sofa
{

namespace component
{

namespace fem
{

template<class DataTypes>
class BezierShellInterpolation : public virtual sofa::core::objectmodel::BaseObject
{
    public:

        SOFA_CLASS(SOFA_TEMPLATE(BezierShellInterpolation,DataTypes), sofa::core::objectmodel::BaseObject);

        typedef typename DataTypes::VecCoord VecCoord;
        typedef typename DataTypes::VecDeriv VecDeriv;
        typedef typename DataTypes::VecReal VecReal;
        typedef typename DataTypes::MatrixDeriv MatrixDeriv;
        typedef typename DataTypes::Coord Coord;
        typedef typename DataTypes::Deriv Deriv;
        typedef typename Coord::value_type Real;

        typedef sofa::Index Index;
        typedef sofa::core::topology::BaseMeshTopology::TriangleID ElementID;
        typedef type::vector<Index> VecIndex;

        typedef typename  sofa::defaulttype::SolidTypes<Real>::Transform Transform;
        typedef typename  sofa::defaulttype::SolidTypes<Real>::SpatialVector SpatialVector;

        typedef sofa::type::Vec<2, Real> Vec2;
        typedef sofa::type::Vec<3, Real> Vec3;
        typedef sofa::type::Mat<3,3,Real> Mat33;
        typedef type::vector<Vec3> VecVec3;

        // These vector is related to Bézier nodes
        typedef type::vector<sofa::type::Vec<3, double> > VecVec3d;

        typedef sofa::type::Vec<10,Real>     ShapeFunctions;
        typedef sofa::type::Vec<10,Index>    BTri;
        typedef type::vector<ShapeFunctions>      VecShapeFunctions;
        typedef type::vector<BTri>                VecBTri;

        class PointInformation
        {
            public:
                // How are the internal bezier nodes attached to the corners of the triangle
                // (valid only for internal edge nodes)
                Vec3 segment;

                PointInformation() {}

                /// Output stream
                inline friend std::ostream& operator<< ( std::ostream& os, const PointInformation& /*pi*/ ) { return os; }
                /// Input stream
                inline friend std::istream& operator>> ( std::istream& in, PointInformation& /*pi*/ ) { return in; }
        };

        class PointInfoHandler : public core::topology::TopologyDataHandler<core::topology::BaseMeshTopology::Point, sofa::type::vector<PointInformation> >
        {
            typedef core::topology::TopologyDataHandler<core::topology::BaseMeshTopology::Point, sofa::type::vector<PointInformation> > Inherited;

            public:
                PointInfoHandler(BezierShellInterpolation<DataTypes> *_bsi, core::topology::PointData<sofa::type::vector<PointInformation> >* _data) : Inherited(_data), bsi(_bsi) {}

                //void applyCreateFunction(
                //    unsigned int pointIndex,
                //    PointInformation &pInfo,
                //    const sofa::type::vector< unsigned int > &ancestors,
                //    const sofa::type::vector< double > &coeffs)
                //{ applyCreateFunction(pointIndex, pInfo, topology::BaseMeshTopology::InvalidID, ancestors, coeffs); }

                //void applyCreateFunction(
                //    unsigned int pointIndex,
                //    PointInformation &pInfo,
                //    const topology::Point &elem,
                //    const sofa::type::vector< unsigned int > &ancestors,
                //    const sofa::type::vector< double > &coeffs);

                void swap( unsigned int i1, unsigned int i2 );

            protected:
                BezierShellInterpolation<DataTypes> *bsi;
        };

        class TriangleInformation
        {
            public:
                BTri btri;

                TriangleInformation() {}

                /// Output stream
                inline friend std::ostream& operator<< ( std::ostream& os, const TriangleInformation& /*ti*/ ) { return os; }
                /// Input stream
                inline friend std::istream& operator>> ( std::istream& in, TriangleInformation& /*ti*/ ) { return in; }
        };

        class TriangleInfoHandler : public core::topology::TopologyDataHandler<core::topology::BaseMeshTopology::Triangle, sofa::type::vector<TriangleInformation> >
        {
            typedef core::topology::TopologyDataHandler<core::topology::BaseMeshTopology::Triangle, sofa::type::vector<TriangleInformation> > Inherited;

            public:
                TriangleInfoHandler(BezierShellInterpolation<DataTypes> *_bsi, core::topology::TriangleData<sofa::type::vector<TriangleInformation> >* _data) : Inherited(_data), bsi(_bsi) {}

                void applyCreateFunction(
                    unsigned int triIndex,
                    TriangleInformation &tInfo,
                    const sofa::type::vector< unsigned int > &ancestors,
                    const sofa::type::vector< double > &coeffs)
                {
                    applyCreateFunction(triIndex, tInfo,
                    core::topology::BaseMeshTopology::Triangle(core::topology::BaseMeshTopology::InvalidID, core::topology::BaseMeshTopology::InvalidID, core::topology::BaseMeshTopology::InvalidID),
                    ancestors, coeffs);
                }

                void applyCreateFunction(
                    unsigned int triIndex,
                    TriangleInformation &tInfo,
                    const core::topology::BaseMeshTopology::Triangle &elem,
                    const sofa::type::vector< unsigned int > &ancestors,
                    const sofa::type::vector< double > &coeffs);

            protected:
                BezierShellInterpolation<DataTypes> *bsi;
        };

        Data< sofa::type::vector<sofa::type::vector< unsigned int > > > f_interpolationIndices; //output interpolation indices
        Data< sofa::type::vector<sofa::type::vector< Real > > > f_interpolationValues;          //output interpolation values

        BezierShellInterpolation();
        ~BezierShellInterpolation()
        {
            if(pointHandler) delete pointHandler;
            if(triHandler) delete triHandler;
        }

        /**
         * @brief SceneGraph callback initialization method.
         */
        void init() override;

        /**
         * @brief SceneGraph callback backward initialization method.
         */
        void bwdInit() override;

        void reinit() override { init(); bwdInit(); }

        /**
         * @brief SceneGraph callback backward reset method.
         */
        //void reset() { bwdInit(); }

        /**
         * @brief SceneGraph callback to handle event
         * Update the positions of Bézier points
         */
        void handleEvent(core::objectmodel::Event *event) override
        {
            if (dynamic_cast< sofa::simulation::AnimateBeginEvent *>(event))
            {
                this->updateBezierPoints();
            }
        }

        void draw(const core::visual::VisualParams* vparams) override;

        //
        // Accessors
        //

        sofa::core::topology::BaseMeshTopology* getInputTopology()
        {
            return this->inputTopology;
        }

        const BTri& getBezierTriangle(unsigned int t) const
        {
            return this->triInfo.getValue()[t].btri;
        }

        //
        // Interface
        //

        void getBezierNodes(ElementID elemID, sofa::type::fixed_array<Vec3,10>& bn)
        {
            const BTri& btri = getBezierTriangle(elemID);
            const VecVec3d& x = mStateNodes->read(sofa::core::vec_id::read_access::position)->getValue();
            for (int i=0; i<10; i++) {
                bn[i] = x[ btri[i] ];
            }
        }

        void getDOFtoLocalTransform(sofa::core::topology::Triangle tri,
            Transform DOF0_H_local0, Transform DOF1_H_local1, Transform DOF2_H_local2);

        void computeShapeFunctions(const Vec3& baryCoord, ShapeFunctions &N);

        //
        // Simple interpolation of a point:
        //
        void interpolateOnBTriangle(Index triID, const VecVec3d& nodes, const ShapeFunctions& N,
            Vec3& point);

        void interpolateOnBTriangle(Index triID, const VecVec3d& nodes, const Vec3& baryCoord,
            Vec3& point)
        {
            ShapeFunctions N;
            computeShapeFunctions(baryCoord, N);
            interpolateOnBTriangle(triID, nodes, N, point);
        }

        void interpolateOnBTriangle(Index triID, const ShapeFunctions& N, Vec3& point) {
            interpolateOnBTriangle(triID, mStateNodes->read(sofa::core::vec_id::read_access::position)->getValue(), N, point);
        }

        void interpolateOnBTriangle(Index triID, const Vec3& baryCoord, Vec3& point) {
            interpolateOnBTriangle(triID, mStateNodes->read(sofa::core::vec_id::read_access::position)->getValue(), baryCoord, point);
        }

        //
        // Interpolation + normal
        //
        void interpolateOnBTriangle(Index triID, const VecVec3d& nodes, const Vec3& baryCoord,
            Vec3& point, Vec3& normal, Vec3& t0, Vec3 &t1);
        void interpolateOnBTriangle(Index triID, const Vec3& baryCoord,
            Vec3& point, Vec3& normal, Vec3& t0, Vec3 &t1) {
            interpolateOnBTriangle(triID, mStateNodes->read(sofa::core::vec_id::read_access::position)->getValue(), baryCoord,
                point, normal, t0, t1);
        }

        // interpolation + normal + second derivatives
        //void interpolateOnBTriangle(Index triID, const VecVec3d& nodes, const Vec3& baryCoord,
        //    Vec3& point, Vec3& normal, Vec3& t0, Vec3 &t1,
        //    Vec3& D2t0, Vec3& D2t01, Vec3& D2t1);

    protected:
        // pointer on mechanical state of the simulation
        sofa::core::behavior::MechanicalState<DataTypes> *mState;
        // pointer to the topology
        sofa::core::topology::BaseMeshTopology* inputTopology;
        sofa::component::mapping::linear::Mesh2PointTopologicalMapping *bezierM2P;

        // Mapping that creates bezier points and MO to store positions of the points
        //sofa::component::topology::Mesh2PointTopologicalMapping::SPtr bezierM2P;
        // Mechanical state holding the Bézier points and their velocities
        sofa::component::statecontainer::MechanicalObject<sofa::defaulttype::Vec3dTypes>::SPtr mStateNodes;

        Data< VecVec3 > inputNormals;

        core::topology::PointData< sofa::type::vector<PointInformation> > pointInfo;
        PointInfoHandler* pointHandler;

        core::topology::TriangleData< sofa::type::vector<TriangleInformation> > triInfo;
        TriangleInfoHandler* triHandler;

        // init process=> computes the position of the bezier point given the positions and the normals
        void computeBezierPointsUsingNormals(const Index& inputTri, VecVec3d& x, const VecVec3& normals);
        void updateBezierPoints();
        void updateBezierPoints(Index triIndex);

        /////// projection of points
        //Data< VecVec3>  pointsToProject; // input : position of the points to project of the surface
        //Data< sofa::type::vector<sofa::type::vector< unsigned int > > > f_interpolationIndices; //output interpolation indices
        //Data< sofa::type::vector<sofa::type::vector< Real > > > f_interpolationValues;          //output interpolation values

        //type::vector< Index > triangleProjectionBuf;  // buffer of the triangle on which was the previous projection of the point
        //VecDeriv pointBaryCoordBuf;       // buffer of barycoordinates to store the projection of point


        //struct SegOfTriInfo {
        //    Index edgeA; //edge that correspond to Alpha<0 (between nodes 1 and 2 of the triangle)
        //    Index edgeB; //edge that correspond to Beta<0 (between nodes 0 and 2 of the triangle)
        //    Index edgeG; //edge that correspond to Gamma<0 (between nodes 0 and 1 of the triangle)
        //};

        //type::vector<SegOfTriInfo> SegOfTriInfoVector;

        //// store for the points inof on the border (especially the id of the borderEdge)
        //struct BorderPointInfo{
        //    Index IdPoint;
        //    type::vector< Index>  borderEdge;
        //};

        //type::vector< BorderPointInfo > PointsOnBorderInfo;

        void initTriangle(Index triIndex, TriangleInformation &tInfo);

        const Vec3& getSegment(Index point)
        {
            return this->pointInfo.getValue()[point].segment;
        }


};

} // namespace fem

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_H
