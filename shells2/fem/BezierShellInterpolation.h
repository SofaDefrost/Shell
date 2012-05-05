//
// Interpolation of Bézier triangles
//
// Author: Tomáš Golembiovský
//
// Copyright:
//
#ifndef SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_H
#define SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_H

#include "../../initPluginShells.h"
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/component/container/MechanicalObject.h>
#include <sofa/component/topology/Mesh2PointTopologicalMapping.h>
//#include <sofa/core/behavior/Mass.h>
//#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/component/topology/TopologyData.inl>
#include <sofa/simulation/common/AnimateBeginEvent.h>

#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
//#include <sofa/defaulttype/Mat.h>

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
        typedef typename DataTypes::Coord Coord;
        typedef typename DataTypes::Deriv Deriv;
        typedef typename Coord::value_type Real;

        typedef sofa::core::topology::BaseMeshTopology::index_type Index;
        typedef sofa::core::topology::BaseMeshTopology::TriangleID ElementID;
        //typedef sofa::helper::vector<BaseMeshTopology::EdgeID> VecElementID;
        //typedef sofa::helper::vector<BaseMeshTopology::Edge> VecEdges;
        typedef helper::vector<Index> VecIndex;

        typedef typename  sofa::defaulttype::SolidTypes<Real>::Transform Transform;
        typedef typename  sofa::defaulttype::SolidTypes<Real>::SpatialVector SpatialVector;

        typedef sofa::defaulttype::Vec<2, Real> Vec2;
        typedef sofa::defaulttype::Vec<3, Real> Vec3;
        typedef sofa::defaulttype::Mat<3,3,Real> Mat33;
        typedef helper::vector<Vec3> VecVec3;

        // These vector is related to Bézier nodes
        typedef helper::vector<sofa::defaulttype::Vec<3, double> > VecVec3d;

        //for newton iterations
        //typedef Mat<2,2,Real> Mat22;
        //bezier segment
        //typedef sofa::defaulttype::Vec<4,Index> BSeg;
        //typedef helper::vector<BSeg> VecBSeg;
        //bezier triangle
        typedef sofa::defaulttype::Vec<10,Index> BTri;
        typedef helper::vector<BTri> VecBTri;

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

        class PointInfoHandler : public topology::TopologyDataHandler<topology::Point, sofa::helper::vector<PointInformation> >
        {
            typedef topology::TopologyDataHandler<topology::Point, sofa::helper::vector<PointInformation> > Inherited;

            public:
                PointInfoHandler(BezierShellInterpolation<DataTypes> *_bsi, topology::PointData<sofa::helper::vector<PointInformation> >* _data) : Inherited(_data), bsi(_bsi) {}

                //void applyCreateFunction(
                //    unsigned int pointIndex,
                //    PointInformation &pInfo,
                //    const sofa::helper::vector< unsigned int > &ancestors,
                //    const sofa::helper::vector< double > &coeffs)
                //{ applyCreateFunction(pointIndex, pInfo, topology::BaseMeshTopology::InvalidID, ancestors, coeffs); }

                //void applyCreateFunction(
                //    unsigned int pointIndex,
                //    PointInformation &pInfo,
                //    const topology::Point &elem,
                //    const sofa::helper::vector< unsigned int > &ancestors,
                //    const sofa::helper::vector< double > &coeffs);

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

        class TriangleInfoHandler : public topology::TopologyDataHandler<topology::Triangle, sofa::helper::vector<TriangleInformation> >
        {
            typedef topology::TopologyDataHandler<topology::Triangle, sofa::helper::vector<TriangleInformation> > Inherited;

            public:
                TriangleInfoHandler(BezierShellInterpolation<DataTypes> *_bsi, topology::TriangleData<sofa::helper::vector<TriangleInformation> >* _data) : Inherited(_data), bsi(_bsi) {}

                void applyCreateFunction(
                    unsigned int triIndex,
                    TriangleInformation &tInfo,
                    const sofa::helper::vector< unsigned int > &ancestors,
                    const sofa::helper::vector< double > &coeffs)
                {
                    applyCreateFunction(triIndex, tInfo,
                    topology::Triangle(topology::BaseMeshTopology::InvalidID, topology::BaseMeshTopology::InvalidID, topology::BaseMeshTopology::InvalidID),
                    ancestors, coeffs);
                }

                void applyCreateFunction(
                    unsigned int triIndex,
                    TriangleInformation &tInfo,
                    const topology::Triangle &elem,
                    const sofa::helper::vector< unsigned int > &ancestors,
                    const sofa::helper::vector< double > &coeffs);

            protected:
                BezierShellInterpolation<DataTypes> *bsi;
        };

        Data< sofa::helper::vector<sofa::helper::vector< unsigned int > > > f_interpolationIndices; //output interpolation indices
        Data< sofa::helper::vector<sofa::helper::vector< Real > > > f_interpolationValues;          //output interpolation values

        BezierShellInterpolation();
        ~BezierShellInterpolation()
        {
            if(pointHandler) delete pointHandler;
            if(triHandler) delete triHandler;
        }

        virtual std::string getTemplateName() const
        {
            return templateName(this);
        }

        static std::string templateName(const BezierShellInterpolation<DataTypes>* = NULL)
        {
            return DataTypes::Name();
        }

        //// Pre-construction check method called by ObjectFactory.
        //// Check that DataTypes matches the MechanicalState.
        //template<class T>
        //static bool canCreate(T* obj, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg)
        //{
        //    if (dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes>*>(context->getMechanicalState()) == NULL)
        //    {
        //        return false;
        //    }
        //    return BaseObject::canCreate(obj, context, arg);
        //}

        /**
         * @brief SceneGraph callback initialization method.
         */
        void init();

        /**
         * @brief SceneGraph callback backward initialization method.
         */
        void bwdInit();

        void reinit() { init(); bwdInit(); }

        /**
         * @brief SceneGraph callback backward reset method.
         */
        //void reset() { bwdInit(); }

        /**
         * @brief SceneGraph callback to handle event
         * Update the positions of Bézier points
         */
        void handleEvent(core::objectmodel::Event *event)
        {
            if (dynamic_cast< sofa::simulation::AnimateBeginEvent *>(event))
            {
                this->updateBezierPoints();
            }
        }

        void draw(const core::visual::VisualParams* vparams);

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

        void getBezierNodes(ElementID elemID, sofa::helper::fixed_array<Vec3,10>& bn)
        {
            const BTri& btri = getBezierTriangle(elemID);
            const VecVec3d& x = *mStateNodes->getX();
            for (int i=0; i<10; i++) {
                bn[i] = x[ btri[i] ];
            }
        }

        void getDOFtoLocalTransform(sofa::core::topology::Triangle tri,
            Transform DOF0_H_local0, Transform DOF1_H_local1, Transform DOF2_H_local2);

        // simple interpolation of a point:
        void interpolateOnBTriangle(Index triID, const VecVec3d& nodes, const Vec3& baryCoord,
            Vec3& point);
        void interpolateOnBTriangle(Index triID, const Vec3& baryCoord, Vec3& point) {
            interpolateOnBTriangle(triID, *mStateNodes->getX(), baryCoord, point);
        }

        // interpolation + normal
        void interpolateOnBTriangle(Index triID, const VecVec3d& nodes, const Vec3& baryCoord,
            Vec3& point, Vec3& normal, Vec3& t0, Vec3 &t1);

        // interpolation + normal + second derivatives
        //void interpolateOnBTriangle(Index triID, const VecVec3d& nodes, const Vec3& baryCoord,
        //    Vec3& point, Vec3& normal, Vec3& t0, Vec3 &t1,
        //    Vec3& D2t0, Vec3& D2t01, Vec3& D2t1);

        void applyOnBTriangle(VecVec3 projBaryCoords, VecIndex projElements, VecVec3& out);
        void applyJOnBTriangle(VecVec3 projBaryCoords, VecIndex projElements, const VecDeriv& in, VecVec3& out);
        void applyJTOnBTriangle(VecVec3 projBaryCoords, VecIndex projElements, const VecVec3& in, VecDeriv& out);

    protected:
        // pointer on mechanical state holding the Bézier points and to
        // mechanical state of the simulation
        sofa::core::behavior::MechanicalState<DataTypes> *mState;
        sofa::core::behavior::MechanicalState<sofa::defaulttype::Vec3dTypes> *mStateNodes;
        // pointer to the topology
        sofa::core::topology::BaseMeshTopology* inputTopology;
        sofa::component::topology::Mesh2PointTopologicalMapping *bezierM2P;

        // Mapping that creates bezier points and MO to store positions of the points
        //sofa::component::topology::Mesh2PointTopologicalMapping::SPtr bezierM2P;
        //sofa::component::container::MechanicalObject<sofa::defaulttype::Vec3dTypes>::SPtr bezierState;

        Data< VecVec3 > inputNormals;

        topology::PointData< sofa::helper::vector<PointInformation> > pointInfo;
        PointInfoHandler* pointHandler;

        topology::TriangleData< sofa::helper::vector<TriangleInformation> > triInfo;
        TriangleInfoHandler* triHandler;

        // init process=> computes the position of the bezier point given the positions and the normals
        void computeBezierPointsUsingNormals(const Index& inputTri, VecVec3d& x, const VecVec3& normals);
        void updateBezierPoints();
        void updateBezierPoints(Index triIndex);

        /////// projection of points
        //Data< VecVec3>  pointsToProject; // input : position of the points to project of the surface
        //Data< sofa::helper::vector<sofa::helper::vector< unsigned int > > > f_interpolationIndices; //output interpolation indices
        //Data< sofa::helper::vector<sofa::helper::vector< Real > > > f_interpolationValues;          //output interpolation values

        //helper::vector< Index > triangleProjectionBuf;  // buffer of the triangle on which was the previous projection of the point
        //VecDeriv pointBaryCoordBuf;       // buffer of barycoordinates to store the projection of point


        //struct SegOfTriInfo {
        //    Index edgeA; //edge that correspond to Alpha<0 (between nodes 1 and 2 of the triangle)
        //    Index edgeB; //edge that correspond to Beta<0 (between nodes 0 and 2 of the triangle)
        //    Index edgeG; //edge that correspond to Gamma<0 (between nodes 0 and 1 of the triangle)
        //};

        //helper::vector<SegOfTriInfo> SegOfTriInfoVector;

        //// store for the points inof on the border (especially the id of the borderEdge)
        //struct BorderPointInfo{
        //    Index IdPoint;
        //    helper::vector< Index>  borderEdge;
        //};

        //helper::vector< BorderPointInfo > PointsOnBorderInfo;

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
