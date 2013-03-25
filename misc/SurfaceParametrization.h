//
// Class for parametrizing 3D surface in 2D domain.
//
// WARNING: Right now we assume surface of genus 1 (i.e. single boundary).
//          The method is very naive and may not work for all topologies.
//

#ifndef SURFACEPARAMETRIZATION_H
#define SURFACEPARAMETRIZATION_H

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

        SurfaceParametrization() : topology(NULL) {}

        void init(
            sofa::component::topology::TriangleSetTopologyContainer *_topology,
            const VecVec3 &x);

        void draw(const core::visual::VisualParams* vparams);

    private:

        enum ElementState {
            FREE,
            BOUNDARY,
            FIXED
        };

        sofa::component::topology::TriangleSetTopologyContainer *topology;

        VecVec2 points;
        helper::vector<bool> ptDone, triDone, triBoundary; // Data for initialization.


        void projectPoint(const Index pId, Vec2 &pt, const VecVec3 &x) const;
        void computeFrame(Mat33 &frame, const Triangle &t, const VecVec3 &x) const;
        ElementState getEdgeState(Index edgeId) const;
        ElementState getTriangleState(Index triId) const;
        void getAngle(const Vec2 &u, const Vec2 &v, Real &alpha, Real &calpha) const;
};

}

#endif // #ifndef SURFACEPARAMETRIZATION_H
