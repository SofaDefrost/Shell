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

#include <sofa/component/topology/TriangleSetTopologyContainer.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/Vec.h>

#include <sofa/helper/vector.h>

namespace sofa
{

template <class Real>
class SurfaceParametrization
{
    public:

        typedef sofa::defaulttype::Vec<2, Real> Vec2;
        typedef sofa::defaulttype::Vec<3, Real> Vec3;

        typedef sofa::defaulttype::Mat<2,2, Real> Mat22;
        typedef sofa::defaulttype::Mat<3,3, Real> Mat33;

        typedef sofa::component::topology::TriangleSetTopologyContainer::Edge                   Edge;
        typedef sofa::component::topology::TriangleSetTopologyContainer::TriangleID             Index;
        typedef sofa::component::topology::TriangleSetTopologyContainer::Triangle               Triangle;
        typedef sofa::component::topology::TriangleSetTopologyContainer::TrianglesAroundVertex  TrianglesAroundVertex;
        typedef sofa::component::topology::TriangleSetTopologyContainer::TrianglesAroundEdge    TrianglesAroundEdge;
        typedef sofa::component::topology::TriangleSetTopologyContainer::EdgesInTriangle        EdgesInTriangle;
        typedef sofa::component::topology::TriangleSetTopologyContainer::SeqEdges               SeqEdges;
        typedef sofa::component::topology::TriangleSetTopologyContainer::SeqTriangles           SeqTriangles;

        typedef sofa::helper::vector<Vec2>      VecVec2;
        typedef sofa::helper::vector<Vec3>      VecVec3;
        typedef sofa::helper::vector<Index>     VecIndex;

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
            const sofa::helper::vector< unsigned int > &ancestors,
            const sofa::helper::vector< double > &coeffs);
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

        const VecVec2& getPositions() { return m_points; }

        /**
         * @brief Moves point in parameter spac.
         *
         * @param p         Point to move.
         * @param newPos    Target position.
         * @param targetTri Triangle in which the new position lies.
         */
        void movePoint(Index p, const Vec2 &newPos, Index targetTri);

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

        // Data for initialization.
        helper::vector<bool> ptDone, triDone, triBoundary;

        Index   m_storedId;
        Vec2    m_storedPos;
        //Mat22 m_storedMetric;

        void projectPoint(const Index pId, Vec2 &pt, const VecVec3 &x) const;
        void computeFrame(Mat33 &frame, const Triangle &t, const VecVec3 &x) const;
        ElementState getEdgeState(Index edgeId) const;
        ElementState getTriangleState(Index triId) const;
        void getAngle(const Vec2 &u, const Vec2 &v, Real &alpha, Real &calpha) const;
};

}

#endif // #ifndef SURFACEPARAMETRIZATION_H
