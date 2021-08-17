//
// Class for parametrizing 3D surface in 2D domain.
//
// WARNING: Right now we assume surface of genus 1 (i.e. single boundary).
//          The method is very naive and may not work for all topologies.
//

#ifndef SURFACEPARAMETRIZATION_H
#define SURFACEPARAMETRIZATION_H

#include <sofa/core/topology/TopologyHandler.h>
#include <sofa/core/visual/VisualParams.h>

#include <SofaBaseTopology/TriangleSetTopologyContainer.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/Vec.h>

#include <sofa/type/vector.h>

namespace sofa
{

/**
 * @brief Parametrization of 3D surface in 2D domain.
 *
 * Parametrization of 3D surface in 2D domain. Serves also as a base class for
 * more specialized projections.
 *
 * @tparam Real Real type to use.
 */
template <class Real>
class SurfaceParametrization
{
    public:

        typedef sofa::type::Vec<2, Real> Vec2;
        typedef sofa::type::Vec<3, Real> Vec3;

        typedef sofa::type::Mat<2,2, Real> Mat22;
        typedef sofa::type::Mat<3,3, Real> Mat33;

        typedef sofa::component::topology::TriangleSetTopologyContainer::Edge                   Edge;
        typedef sofa::component::topology::TriangleSetTopologyContainer::TriangleID             Index;
        typedef sofa::component::topology::TriangleSetTopologyContainer::Triangle               Triangle;
        typedef sofa::component::topology::TriangleSetTopologyContainer::TrianglesAroundVertex  TrianglesAroundVertex;
        typedef sofa::component::topology::TriangleSetTopologyContainer::TrianglesAroundEdge    TrianglesAroundEdge;
        typedef sofa::component::topology::TriangleSetTopologyContainer::EdgesInTriangle        EdgesInTriangle;
        typedef sofa::component::topology::TriangleSetTopologyContainer::SeqEdges               SeqEdges;
        typedef sofa::component::topology::TriangleSetTopologyContainer::SeqTriangles           SeqTriangles;

        typedef sofa::type::vector<Vec2>      VecVec2;
        typedef sofa::type::vector<Vec3>      VecVec3;
        typedef sofa::type::vector<Index>     VecIndex;
        typedef sofa::type::vector<Mat22>     VecMat22;

        enum { InvalidID = sofa::core::topology::Topology::InvalidID };

        SurfaceParametrization() : m_topology(NULL), m_storedId(InvalidID) {}

        /**
         * @brief Initialize parametrization.
         *
         * @param _topology Topology defining the surface.
         * @param x         Position of points of the surface.
         */
        void init(
            sofa::component::topology::TriangleSetTopologyContainer *_topology,
            const VecVec3 &x);

        void draw(const core::visual::VisualParams* vparams);

        /**
         * @name Hadling of topological changes (also handles point relocation).
         * @{
         */
        void pointAdd(unsigned int pointIndex, const sofa::core::topology::Point &elem,
            const sofa::type::vector< unsigned int > &ancestors,
            const sofa::type::vector< double > &coeffs);
        void pointRemove(unsigned int pointIndex);
        void pointSwap(unsigned int i1, unsigned int i2);

        /**  @} */

        /**
         * @brief Store point parameters for later use.
         *
         * Store point parameters so they can later be restored. Only one
         * point can be stored at a time.
         *
         * @param pt Index of a point to store.
         */
        void storePoint(Index pt) {
            m_storedId = pt;
            m_storedPos = m_points[pt];
        }

        /**
         * @brief Restore previously stored point parameters.
         */
        void restorePoint() {
            if (m_storedId == InvalidID) return;
            m_points[m_storedId] = m_storedPos;
            m_storedId = InvalidID;
        }

        /**
         * @brief Returns position of a point in 3D.
         *
         * @param position  Position in parameter space.
         * @param tId       Triangle in which the point is.
         * @param x         Position of surface points in real space. Note
         *                  that it should be same as used in init().
         *
         * @return Coordinates of the point in real space.
         */
        Vec3 getPointPosition(Vec2 position, Index tId, const VecVec3 &x) const;

        /**
         * @brief Positions of points in parameter space.
         */
        const VecVec2& getPositions() { return m_points; }

        /**
         * @brief Moves point in parameter spac.
         *
         * @param p         Point to move.
         * @param newPos    Target position.
         * @param targetTri Triangle in which the new position lies.
         */
        void movePoint(Index p, const Vec2 &newPos, Index targetTri);

        /**
         * @name Distances and areas.
         * @{ */

        // NOTE: Big note! The following is (in a strict point of view) invalid
        // because it assumes that the metric tensor varies linearly over the
        // edge/triangle. For our barycentric mapping however the resulting
        // metric space is piecewise-linear (for original metric) or it is
        // constant over each triangle (for updated metric). But since
        // piecewise-linear space would be difficult to integrate we use linear
        // interpolation because it is more versatile than constant metric.

        /**
         * @brief Compute length of an edge.
         *
         * @param e Edge.
         */
        Real dist(const Edge &e) const {
            return helper::rsqrt(dist2(e));
        }

        /**
         * @brief Compute square of the length of an edge.
         *
         * @param e Edge.
         */
        Real dist2(const Edge &e) const {
            Vec2 v = m_points[ e[1] ] - m_points[ e[0] ];
            Mat22 M = (m_metrics[ e[1] ] + m_metrics[ e[0] ])/2.0;
            return v*(M*v);
        }

        /**
         * @brief Compute area of the triangle.
         *
         * @param t Triangle.
         */
        Real area(const Triangle &t) const {

            Real area = helper::rabs(cross(
                    m_points[t[1]] - m_points[t[0]],
                    m_points[t[2]] - m_points[t[0]]))/2.0;

            // === Reddy. Introduction to the Finite Element Method. 1993
            // n = 1, d = 1 
            //const double GX[] = { 1.0/3.0 };
            //const double GY[] = { 1.0/3.0 };
            //const double GW[] = { 1.0 };
            // n = 3, d = 2
            const double GX[] = { 0.5,     0.5,     0.0 };
            const double GY[] = { 0.5,     0.0,     0.5 };
            const double GW[] = { 1.0/3.0, 1.0/3.0, 1.0/3.0 };
            // n = 4, d = 3
            //const double GX[] = { 1.0/3.0,    0.6,       0.2,       0.2 };
            //const double GY[] = { 1.0/3.0,    0.2,       0.6,       0.2 };
            //const double GW[] = { -27.0/48.0, 25.0/48.0, 25.0/48.0, 25.0/48.0 };
            // n = 7, d = 5
            //const double GX[] = { 1.0/3.0, 0.797426985353, 0.101286507323, 0.101286507323, 0.059715871789, 0.470142064105, 0.470142064105 };
            //const double GY[] = { 1.0/3.0, 0.101286507323, 0.797426985353, 0.101286507323, 0.470142064105, 0.059715871789, 0.470142064105 };
            //const double GW[] = { 0.225,   0.125939180544, 0.125939180544, 0.125939180544, 0.132394152788, 0.132394152788, 0.132394152788 };
            // === http://www.electromagnetics.biz/integration.htm
            // NOTE: The follwing assume 2A*I ... the weights sum to 0.5
            // n = 6, d = 4 
            //const double GX[] = { 0.8168476,  0.09157621, 0.09157621, 0.1081030, 0.4459485, 0.4459485 };
            //const double GY[] = { 0.09157621, 0.8168476,  0.09157621, 0.4459485, 0.1081030, 0.4459485 };
            //const double GW[] = { 0.05497587, 0.05497587, 0.05497587, 0.1116908, 0.1116908, 0.1116908 };
            // n = 12, d = 6
            //const double GX[] = { 0.5014265,  0.2492867,  0.2492867,  0.8738220,  0.06308901, 0.06308901, 0.6365025,  0.6365025,  0.05314505, 0.05314505, 0.3103525,  0.3103525  };
            //const double GY[] = { 0.2492867,  0.5014265,  0.2492867,  0.06308901, 0.8738220,  0.06308901, 0.05314505, 0.3103525,  0.6365025,  0.3103525,  0.6365025,  0.05314505 };
            //const double GW[] = { 0.05839314, 0.05839314, 0.05839314, 0.02542245, 0.02542245, 0.02542245, 0.04142554, 0.04142554, 0.04142554, 0.04142554, 0.04142554, 0.04142554 };
            // n = 27, d = 11
            //const double GX[] = { 0.9352701,   0.03236495,  0.03236495,  0.7612982,  0.1193509,  0.1193509,  0.06922210,   0.5346110, 0.5346110, 0.5933802, 0.2033099, 0.2033099, 0.2020614, 0.3989693, 0.3989693, 0.05017814, 0.05017814, 0.5932012, 0.5932012, 0.3566206, 0.3566206, 0.02102202, 0.02102202, 0.8074890, 0.8074890, 0.1714890, 0.1714890 };
            //const double GY[] = { 0.03236495,  0.9352701,   0.03236495,  0.1193509,  0.7612982,  0.1193509,  0.5346110,    0.06922210, 0.5346110, 0.2033099, 0.5933802, 0.2033099, 0.3989693, 0.2020614, 0.3989693, 0.5932012, 0.3566206, 0.05017814, 0.3566206, 0.05017814, 0.5932012, 0.8074890, 0.1714890, 0.02102202, 0.1714890, 0.02102202, 0.8074890 };
            //const double GW[] = { 0.006829866, 0.006829866, 0.006829866, 0.01809227, 0.01809227, 0.01809227, 0.0004635032, 0.0004635032, 0.0004635032, 0.02966149, 0.02966149, 0.02966149, 0.03857477, 0.03857477, 0.03857477, 0.02616856, 0.02616856, 0.02616856, 0.02616856, 0.02616856, 0.02616856, 0.01035383, 0.01035383, 0.01035383, 0.01035383, 0.01035383, 0.01035383 };

            const int N=3;
            Real I = 0.0;
            //std::cout <<
            //    determinant(M1b) << " " <<
            //    determinant(M2b) << " " <<
            //    determinant(M3b) << "\n  :: ";
            for (int i = 0; i < N; i++) {
                Mat22 M = m_metrics[t[0]]*GX[i] +
                    m_metrics[t[1]]*GY[i] +
                    m_metrics[t[2]]*(1.0-GX[i]-GY[i]);
                I += GW[i] * sqrt(determinant(M));
                //std::cout << determinant(M) << " ";
            }
            //std::cout << "\n  :: " << I << "\n";
            return area*I;
        }

        /**
         * @brief Return metric tensor of a point inside the triangle.
         *
         * @param tId   Triangle ID.
         * @param bary  Barycentric coordinates of the point.
         * @param M     Computed metric tensor.
         */
        virtual void getMetricTensor(Index tId, const Vec3 &bary, Mat22 &M) const;

        /**  @} */

    private:

        enum ElementState {
            FREE,
            BOUNDARY,
            FIXED
        };

        /// Topology defining the surface.
        sofa::component::topology::TriangleSetTopologyContainer *m_topology;

        /// Positions in the parameter space.
        VecVec2 m_points;

        /// Metric tensors (first fundamental form).
        VecMat22 m_metrics;

        // Data for initialization.
        type::vector<bool> ptDone, triDone, triBoundary;

        Index   m_storedId;
        Vec2    m_storedPos;
        //Mat22 m_storedMetric;

        virtual void initMetricTensors();
        void projectPoint(const Index pId, Vec2 &pt, const VecVec3 &x) const;
        void computeFrame(Mat33 &frame, const Triangle &t, const VecVec3 &x) const;
        ElementState getEdgeState(Index edgeId) const;
        ElementState getTriangleState(Index triId) const;
        void getAngle(const Vec2 &u, const Vec2 &v, Real &alpha, Real &calpha) const;
};

}

#endif // #ifndef SURFACEPARAMETRIZATION_H
