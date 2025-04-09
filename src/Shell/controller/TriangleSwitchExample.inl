#ifndef SOFA_COMPONENT_CONTROLLER_TRIANGLESWITCHEXAMPLE_INL
#define SOFA_COMPONENT_CONTROLLER_TRIANGLESWITCHEXAMPLE_INL

#include <SofaShells/controller/TriangleSwitchExample.h>

namespace sofa
{

namespace component
{

namespace controller
{

template<class DataTypes>
TriangleSwitchExample<DataTypes>::TriangleSwitchExample()
: f_interval(initData(&f_interval, (unsigned int)10, "interval", "Switch triangles every number of steps"))
, stepCounter(0)
{
}

template<class DataTypes>
TriangleSwitchExample<DataTypes>::~TriangleSwitchExample()
{
}

template<class DataTypes>
void TriangleSwitchExample<DataTypes>::init()
{
    m_state = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes>*> (this->getContext()->getMechanicalState());
    if (!m_state)
        msg_error() << "Unable to find MechanicalState.";

    this->getContext()->get(m_container);
    if (m_container == NULL)
        msg_error() << "Unable to find triangular topology.";

    this->getContext()->get(m_modifier);
    if (m_modifier == NULL)
        msg_error() << "Unable to find TriangleSetTopologyModifier.";

    reinit();
}

template<class DataTypes>
void TriangleSwitchExample<DataTypes>::reinit()
{
    *this->f_listening.beginEdit() = true;
    this->f_listening.endEdit();

    // Check interval
    if (f_interval.getValue() == 0) {
        msg_warning() << "interval has to be nonzero.";
        *f_interval.beginEdit() = 1;
        f_interval.endEdit();
    }

    stepCounter = 0;
}


template<class DataTypes>
void TriangleSwitchExample<DataTypes>::onEndAnimationStep(const double /*dt*/)
{
    if ((m_container == NULL) || (m_modifier == NULL) || (m_state == NULL))
        return;

    stepCounter++;

    if (stepCounter < f_interval.getValue())
        return; // Stay still for a few steps

    stepCounter = 0;


    SeqTriangles add_list;
    VecIndex remove_list;

    for (Index i=0; i < m_container->getNbTriangles()-1; i+=2) {
        Triangle t1, t2, t1n, t2n;
        t1 = m_container->getTriangle(i);
        t2 = m_container->getTriangle(i+1);

        remove_list.push_back(i);
        remove_list.push_back(i+1);

        // Find common vertices (the shared edge) and the other two
        VecIndex common, other;

        int v1, v2;
        for (v1=0; v1 < 3; v1++) {
            for (v2=0; v2 < 3; v2++) {
                if (t1[v1] == t2[v2]) break;
            }
            if (v2 < 3) {
                // found
                common.push_back(t1[v1]);
            } else {
                // not found
                other.push_back(t1[v1]);
            }
        }

        for (v2=0; v2 < 3; v2++) {
            for (v1=0; v1 < 3; v1++) {
                if (t1[v1] == t2[v2]) break;
            }
            if (v1 >= 3) {
                // not found
                other.push_back(t2[v2]);
                break;
            }
        }

        if ((common.size() != 2) || (other.size() != 2)) {
            msg_error() << "Invalid triangle pair " << i << ", " << i+1;
            continue;
        }

        //std::cout << "common= " << common << " other= " << other << std::endl;

        // Construct new triangles
        t1n[0] = t2n[0] = other[0];
        t1n[1] = t2n[1] = other[1];
        t1n[2] = common[0];
        t2n[2] = common[1];

        // Make sure the vertex order of the triangle is correct

        Vec3 nt1, nt1n;
        computeTriangleNormal(t1, nt1);
        computeTriangleNormal(t1n, nt1n);
        if (dot(nt1, nt1n) > 0) {
            // Invert the vertex order
            Index tmp = t1n[1];
            t1n[1] = t1n[2];
            t1n[2] = tmp;
        }

        Vec3 nt2, nt2n;
        computeTriangleNormal(t2, nt2);
        computeTriangleNormal(t2n, nt2n);
        if (dot(nt2, nt2n) > 0) {
            // Invert the vertex order
            Index tmp = t2n[1];
            t2n[1] = t2n[2];
            t2n[2] = tmp;
        }

        //  Add triangles to the list
        add_list.push_back(t1n);
        add_list.push_back(t2n);
    }

    // Add new and remove old triangles
    m_modifier->addTriangles(add_list);
    m_modifier->removeTriangles(remove_list, true, true);
    m_modifier->notifyEndingEvent();
}

// Computes triangle normal
template<class DataTypes>
void TriangleSwitchExample<DataTypes>::computeTriangleNormal(const Triangle &t, Vec3 &normal)
{
    const VecCoord& x = m_state->read(sofa::core::ConstVecCoordId::restPosition())->getValue();

        //std::cout << "tri: " << t << " n=" << x.size() << std::endl;

    Vec3 A, B;
    A = x[t[1]].getCenter() - x[t[0]].getCenter();
    B = x[t[2]].getCenter() - x[t[0]].getCenter();

    Real An = A.norm(), Bn = B.norm();
    if (An < 1e-20 || Bn < 1e-20) {
        msg_error() << "Found degenerated triangle " << t <<
            " -- "  << x[t[0]].getCenter() << " / " << x[t[1]].getCenter() <<
            " / " << x[t[2]].getCenter() << " :: " << An << ", " << Bn;

        normal = Vec3(0,0,0);
        return;
    }

    A.normalize();
    B.normalize();
    normal = cross(A, B);
    normal.normalize();
}

} // namespace controller

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONTROLLER_TRIANGLESWITCHEXAMPLE_INL
