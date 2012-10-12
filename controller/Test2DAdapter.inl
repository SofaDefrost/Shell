#ifndef SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_INL
#define SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_INL

#include "Test2DAdapter.h"

namespace sofa
{

namespace component
{

namespace controller
{

template<class DataTypes>
Test2DAdapter<DataTypes>::Test2DAdapter()
//: f_interval(initData(&f_interval, (unsigned int)10, "interval", "Switch triangles every number of steps"))
//, stepCounter(0)
{
}

template<class DataTypes>
Test2DAdapter<DataTypes>::~Test2DAdapter()
{
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::init()
{
    m_state = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes>*> (this->getContext()->getMechanicalState());
    if (!m_state)
        serr << "Unable to find MechanicalState" << sendl;

    this->getContext()->get(m_container);
    if (m_container == NULL)
        serr << "Unable to find triangular topology" << sendl;

    reinit();
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::reinit()
{
    *this->f_listening.beginEdit() = true;
    this->f_listening.endEdit();
    
    detectBoundary();

    //stepCounter = 0;
}


template<class DataTypes>
void Test2DAdapter<DataTypes>::onEndAnimationStep(const double /*dt*/)
{
    if ((m_container == NULL) || (m_state == NULL))
        return;

    //stepCounter++;

    //if (stepCounter < f_interval.getValue())
    //    return; // Stay still for a few steps

    //stepCounter = 0;

    Data<VecVec3>* datax = m_state->write(sofa::core::VecCoordId::position());
    VecVec3& x = *datax->beginEdit();

    Index nTriangles = m_container->getNbTriangles();
    if (nTriangles == 0)
        return;
    
    vector<Real> metrics(nTriangles);
    vector<Vec3> normals(nTriangles);

    // Compute initial metrics and normals
    for (Index i=0; i < nTriangles; i++) {
        Triangle t = m_container->getTriangle(i);
        computeTriangleNormal(t, x, normals[i]);
        metrics[i] = metricGeom(t, x, normals[i]);
    }
    std::cout << "m: " << metrics << "\n";

    for (Index i=0; i<x.size(); i++) {
        if (m_boundary[i]) {
            std::cout << "skipping boundary " << i << "\n";
            continue;
        }

        Vec3 xold = x[i];
        if (!smoothLaplacian(i, x, metrics, normals)) {
            x[i] = xold;
        }
    }

    datax->endEdit();
}

template<class DataTypes>
bool Test2DAdapter<DataTypes>::smoothLaplacian(Index v, VecVec3 &x, vector<Real>metrics, vector<Vec3> normals)
{
    Vec3 xold = x[v];

    // Compute new position
    EdgesAroundVertex N1 = m_container->getEdgesAroundVertex(v);

    // Compute centroid of polygon from 1-ring around the vertex
    Vec3 xnew(0,0,0);
    for (Index ie=0; ie<N1.size(); ie++) {
        Edge e = m_container->getEdge(N1[ie]); 
        for (int n=0; n<2; n++) {
            if (e[n] != v) {
                xnew += x[ e[n] ];
            }
        }
    }
    x[v] = xnew / N1.size();

    // Check if this improves the mesh
    //
    // Note: To track element inversion we either need a normal computed
    // from vertex normals, or assume the triangle was originaly not
    // inverted. Now we do the latter.
    //
    // We accept any change that doesn't decreas worst metric for the
    // triangle set.

    bool bAccepted = false;
    for (int iter=10; iter>0 && !bAccepted; iter--) {

        if ((xold - x[v]).norm2() < 1e-12) {
            // No change in position
            std::cout << "No change in position for " << v << "\n";
            break;
        }

        Real oldworst = 1.0, newworst = 1.0;
        for (Index it=0; it<N1.size(); it++) {
            Real newmetric = metricGeom(m_container->getTriangle(it), x,
                normals[v]);

            if (metrics[it] < oldworst) {
                oldworst = metrics[it];
            }

            if (newmetric < newworst) {
                newworst = newmetric;
            }
        }
        std::cout << "cmp: " << newworst << " vs. " << oldworst << "\n";
        if (newworst < oldworst) {
            std::cout << "   --rejected " << xold << " -> " << x[v] << "\n";
            x[v] = (x[v] + xold)/2.0;
        } else {
            std::cout << "   --accepted: " << xold << " -> " << x[v] << "\n";
            bAccepted = true;
        }
    }
    return bAccepted;
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::detectBoundary()
{
    if ((m_container == NULL) || (m_container->getNbTriangles() == 0))
        return;

    // TODO: maybe use TriangleSetTopologyContainer::getPointsOnBorder()
    m_boundary.resize(m_container->getNbPoints());
    for (Index i=0; i<(Index)m_container->getNbEdges(); i++) {
        TrianglesAroundEdge tlist = m_container->getTrianglesAroundEdge(i);
        if (tlist.size() != 2) {
            Edge e = m_container->getEdge(i);
            m_boundary[e[0]] = true;
            m_boundary[e[1]] = true;
        }
    }
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::computeTriangleNormal(const Triangle &t, const VecVec3 &x, Vec3 &normal)
{
    Vec3 A, B;
    A = x[ t[1] ] - x[ t[0] ];
    B = x[ t[2] ] - x[ t[0] ];

    Real An = A.norm(), Bn = B.norm();
    if (An < 1e-20 || Bn < 1e-20) {
        serr << "Found degenerated triangle: "
            << x[ t[0] ] << " / " << x[ t[1] ] << " / " << x[ t[2] ]
            << " :: " << An << ", " << Bn << sendl;

        normal = Vec3(0,0,0);
        return;
    }

    A /= An;
    B /= Bn;
    normal = cross(A, B);
    normal.normalize();
}

} // namespace controller

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_INL
