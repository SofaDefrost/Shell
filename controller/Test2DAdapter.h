#ifndef SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_H
#define SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_H

#include <sofa/component/component.h>
#include <sofa/component/controller/Controller.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/Vec.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/component/topology/TriangleSetTopologyContainer.h>

#include <sofa/helper/vector.h>

namespace sofa
{

namespace component
{

namespace controller
{

/**
 *
 * @brief Component for adaptivity/smoothing of 2D triangular meshes.
 *
 * References:
 *
 * [CTS98] Canann, S. A.; Tristano, J. R. & Staten, M. L. An Approach to
 *         Combined Laplacian and Optimization-Based Smoothing for Triangular,
 *         Quadrilateral, and Quad-Dominant Meshes International Meshing
 *         Roundtable, 1998, 479-494.
 *
 * [VL99] Y. Vasilevskii, K. Lipnikov, An adaptive algorithm for quasioptimal
 *        mesh generation, Computational mathematics and mathematical physics
 *        39 (9) (1999) 1468–1486.
 */
template<class DataTypes>
class Test2DAdapter : public Controller
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(Test2DAdapter,DataTypes),Controller);

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;

    //typedef sofa::defaulttype::Vec<2, Real> Vec2;
    typedef sofa::defaulttype::Vec<3, Real> Vec3;
    //typedef sofa::defaulttype::Mat<3,3,Real> Mat33;
    //typedef helper::vector<Vec2> VecVec2;
    typedef helper::vector<Vec3> VecVec3;


    typedef sofa::component::topology::EdgeSetTopologyContainer::Edge               Edge;
    typedef sofa::component::topology::EdgeSetTopologyContainer::EdgesAroundVertex  EdgesAroundVertex;
    typedef sofa::component::topology::TriangleSetTopologyContainer::TriangleID     Index;
    typedef sofa::component::topology::TriangleSetTopologyContainer::Triangle       Triangle;
    typedef sofa::component::topology::TriangleSetTopologyContainer::TrianglesAroundVertex  TrianglesAroundVertex;
    typedef sofa::component::topology::TriangleSetTopologyContainer::TrianglesAroundEdge    TrianglesAroundEdge;
    typedef sofa::helper::vector<Index> VecIndex;

protected:

    Test2DAdapter();

    virtual ~Test2DAdapter();

public:

    virtual void init();
    virtual void reinit();

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const Test2DAdapter<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    void onEndAnimationStep(const double dt);
    void onKeyPressedEvent(core::objectmodel::KeypressedEvent *key);

    /**
     * @brief Distortion metric for a triangle from [CTS98] (but probably due
     * to somebody else)
     *
     * @param t         Triangle to compute the metric for.
     * @param x         Vertices.
     * @param normal    Surface normal for the triangle (to check for inversion).
     */
    Real metricGeom(const Triangle &t, const VecVec3 &x, const Vec3 &normal) {
        // TODO: we can precompute these
        Vec3 ab = x[ t[1] ] - x[ t[0] ];
        Vec3 ca = x[ t[0] ] - x[ t[2] ];
        Vec3 cb = x[ t[1] ] - x[ t[2] ];

        // Normalizing factor so that the matric is 1 in maximum
        Real m = 2 * sqrt(3);
        m *= ca.cross(cb).norm(); // || CA × CB ||
        m /= ca.norm2() + ab.norm2() + cb.norm2();

        // Is triangle inverted?
        Vec3 nt;
        computeTriangleNormal(t, x, nt);
        if (dot(nt, normal) < 0) {
            m *= -1.0;
        }

        return m;
    }

    /**
     * @brief Distortion metric for a triangle from [VL99] extended for
     * handling of inverted triangles.
     *
     * @param t         Triangle to compute the metric for.
     * @param x         Vertices.
     * @param normal    Surface normal for the triangle (to check for inversion).
     */
    Real metricGeom2(const Triangle &t, const VecVec3 &x, const Vec3 &normal) {
        // TODO: we can precompute these
        Vec3 ab = x[ t[1] ] - x[ t[0] ];
        Vec3 ca = x[ t[0] ] - x[ t[2] ];
        Vec3 cb = x[ t[1] ] - x[ t[2] ];

        Real p = ca.norm() + ab.norm() + cb.norm(); // perimeter
        Real ip = 1.0/p; // inverse
        // Normalizing factor so that the matric is 1 in maximum
        // TODO: verify this, the maximum seems to be at 2
        //Real m = 12 * sqrt(3);
        Real m = 6 * helper::rsqrt(3);
        m *= ca.cross(cb).norm(); // || CA × CB ||
        m /= helper::rsqrt(p);
        Real f;
        if (p < ip) {
            f = (p * (2-p));
        } else {
            f = (ip * (2-ip));
        }

        m *= ip * ip * ip;

        // Is triangle inverted?
        Vec3 nt;
        computeTriangleNormal(t, x, nt);
        if (dot(nt, normal) < 0) {
            m *= -1.0;
        }

        return m;
    }

private:

    unsigned int stepCounter;
    sofa::component::topology::TriangleSetTopologyContainer*  m_container;
    sofa::core::behavior::MechanicalState<DataTypes>* m_state;
    helper::vector<bool> m_boundary; /// marks whether node lies on the boundary

    vector<Real> m_metrics;

    Real sumgamma, mingamma, maxgamma;
    int ngamma;

    /**
     * @brief Constrained Laplacian smoothing
     *
     * @param v         Vertex to move
     * @param x         Current positions
     * @param metrics   Current metrice values for elements
     * @param normals   Original normals (to check for inversion)
     */
    bool smoothLaplacian(Index v, VecVec3 &x, vector<Real>metrics, vector<Vec3> normals);

    /**
     * @brief Optimization based smoothing
     *
     * @param v         Vertex to move
     * @param x         Current positions
     * @param metrics   Current metrice values for elements
     * @param normals   Original normals (to check for inversion)
     */
    bool smoothOptimize(Index v, VecVec3 &x, vector<Real>metrics, vector<Vec3> normals);

    /**
     * @brief Detect which nodes lie on the boundary.
     */
    void detectBoundary();

    /**
     * @brief Compute triangle normal.
     *
     * @param a         First vertex.
     * @param b         Second vertex.
     * @param c         Third vertex.
     * @param normal    The computed normal.
     *
     */
    // TODO: isn't this in geometry algorithms or CGAL?
    void computeTriangleNormal(const Triangle &t, const VecVec3 &x, Vec3 &normal);
};


} // namespace controller

} // namespace component

} // namespace sofa

#endif //SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_H
