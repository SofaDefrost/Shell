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
#ifndef SOFA_COMPONENT_FORCEFIELD_TRIANGULAR_BENDING_FEM_FORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_TRIANGULAR_BENDING_FEM_FORCEFIELD_INL

#include <SofaShells/forcefield/TriangularBendingFEMForceField.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/gl/template.h>
#include <sofa/gl/gl.h>
#include <sofa/helper/system/thread/debug.h>
#include <fstream> // for reading the file
#include <iostream> //for debugging
#include <vector>
#include <algorithm>
#include <sofa/defaulttype/VecTypes.h>
#include <assert.h>
#include <map>
#include <utility>
#include <sofa/core/topology/TopologyData.inl>
#include <sofa/component/topology/container/dynamic/TriangleSetTopologyContainer.h>

#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/AnimateEndEvent.h>

#include <SofaShells/controller/MeshChangedEvent.h>

#ifdef _WIN32
#include <windows.h>
#endif

namespace sofa
{
	namespace component
	{
		namespace forcefield
		{
			using namespace sofa::type;
			using namespace	sofa::component::topology;



// --------------------------------------------------------------------------------------
// ---  Topology Creation/Destruction functions
// --------------------------------------------------------------------------------------
template<class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::TRQSTriangleHandler::applyCreateFunction(unsigned int triangleIndex, TriangleInformation &, const Triangle &t, const sofa::type::vector<unsigned int> &, const sofa::type::vector<double> &)
{
    if (ff)
    {
        ff->initTriangleOnce(triangleIndex, t[0], t[1], t[2]);
        ff->initTriangle(triangleIndex);
        ff->computeMaterialStiffness(triangleIndex);
    }
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
TriangularBendingFEMForceField<DataTypes>::TriangularBendingFEMForceField()
: d_poisson(initData(&d_poisson,(Real)0.45,"poissonRatio","Poisson ratio in Hooke's law"))
, d_young(initData(&d_young,(Real)3000.,"youngModulus","Young modulus in Hooke's law"))
, d_bending(initData(&d_bending,false,"bending","Adds bending"))
, d_thickness(initData(&d_thickness,(Real)0.1,"thickness","Thickness of the plates"))
, d_membraneRatio(initData(&d_membraneRatio,(Real)1.0,"membraneRatio","In plane forces ratio"))
, d_bendingRatio(initData(&d_bendingRatio,(Real)1.0,"bendingRatio","Bending forces ratio"))
, d_refineMesh(initData(&d_refineMesh, false, "refineMesh","Hierarchical refinement of the mesh"))
, d_iterations(initData(&d_iterations,(int)0,"iterations","Iterations for refinement"))
, l_targetTopology(initLink("targetTopology","Targeted high resolution topology"))
, l_restShape(initLink("restShape","MeshInterpolator component for variable rest shape"))
, m_mapTopology(false)
, l_topologyMapper(initLink("topologyMapper","Component supplying different topology for the rest shape"))
, m_exportFilename(initData(&m_exportFilename, "exportFilename", "file name to export coefficients into"))
, d_exportEveryNbSteps(initData(&d_exportEveryNbSteps, (unsigned int)0, "exportEveryNumberOfSteps", "export file only at specified number of steps (0=disable)"))
, d_exportAtBegin(initData(&d_exportAtBegin, false, "exportAtBegin", "export file at the initialization"))
, d_exportAtEnd(initData(&d_exportAtEnd, false, "exportAtEnd", "export file when the simulation is finished"))
, m_stepCounter(0)
, triangleInfo(initData(&triangleInfo, "triangleInfo", "Internal triangle data"))
{
    m_triangleHandler = new TRQSTriangleHandler(this, &triangleInfo);
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes> void TriangularBendingFEMForceField<DataTypes>::handleTopologyChange()
{
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
TriangularBendingFEMForceField<DataTypes>::~TriangularBendingFEMForceField()
{
    if(m_triangleHandler) delete m_triangleHandler;
}

// --------------------------------------------------------------------------------------
// --- Initialization stage
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::init()
{
    this->d_componentState.setValue(core::objectmodel::ComponentState::Valid);
    this->Inherited::init();

    _topology = this->getContext()->getMeshTopology();

    if (_topology->getNbTriangles()==0)
    {
        msg_error() << "TriangularBendingFEMForceField: object must have a Triangular Set Topology.";
        this->d_componentState.setValue(core::objectmodel::ComponentState::Invalid);
        return;
    }

    // Create specific handler for TriangleData
    triangleInfo.createTopologyHandler(_topology);

    reinit();

    if (d_refineMesh.getValue())
    {
        sofa::core::topology::BaseMeshTopology* _topologyTarget = l_targetTopology.get();

        if (_topologyTarget)
        {
            MechanicalState<defaulttype::Vec3Types>* mStateTarget = dynamic_cast<MechanicalState<defaulttype::Vec3Types>*> (_topologyTarget->getContext()->getMechanicalState());
            if (mStateTarget)
            {
                m_targetTriangles = _topologyTarget->getTriangles();
                m_targetVertices = mStateTarget->read(sofa::core::ConstVecCoordId::position())->getValue();
            }
            else
            {
                msg_warning() << "No mechanical state for target high resolution topology" ;
                return;
            }
        }
        else
        {
            msg_warning() << "No target high resolution mesh found";
            return;
        }

        // Run procedure for shell remeshing
        refineCoarseMeshToTarget();
    }
}


template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::refineCoarseMeshToTarget(void)
{
    msg_info() << "Refining a mesh of " << _topology->getNbTriangles() << " triangles towards a target surface of " << m_targetTriangles.size() << " triangles.";

    // List of vertices
    const VecCoord& x = this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();
    // List of triangles
    const SeqTriangles triangles = _topology->getTriangles();

    // Creates new mesh
    sofa::type::vector<Vec3> subVertices;
    SeqTriangles subTriangles;

    // Initialises list of subvertices and triangles
    for (unsigned int i=0; i<x.size(); i++)
    {
        subVertices.push_back(x[i].getCenter());
    }
    for (unsigned int t=0; t<triangles.size(); t++)
    {
        subTriangles.push_back(triangles[t]);
    }

    // Adjusts position of each subvertex to get closer to actual surface before iterating again
    for (unsigned int i=0; i<subVertices.size(); i++)
    {
        movePoint(subVertices[i]);
    }


    // Refines mesh
    for (int n=0; n<d_iterations.getValue(); n++)
    {
        // Subdivides each triangle into 4 smaller ones
        subTriangles.clear();
        for (unsigned int t=0; t<triangles.size(); t++)
        {
            Vec3 a = subVertices[(int)triangles[t][0]];
            Vec3 b = subVertices[(int)triangles[t][1]];
            Vec3 c = subVertices[(int)triangles[t][2]];

            subdivide(a, b, c, subVertices, subTriangles);
        }

        // Adjusts position of each subvertex to get closer to actual surface before iterating again
        for (unsigned int i=0; i<subVertices.size(); i++)
        {
            movePoint(subVertices[i]);
        }
    }


    msg_info() << "Number of vertices of the resulting mesh = " << subVertices.size();
    msg_info() << "Number of shells of the resulting mesh   = " << subTriangles.size();

    // Writes in Gmsh format
    std::ofstream myfile;
    myfile.open ("mesh_refined.obj");
    for (unsigned int vertex=0; vertex<subVertices.size(); vertex++)
    {
        myfile << "v " << subVertices[vertex] << "\n";
    }
    for (unsigned int element=0; element<subTriangles.size(); element++)
    {
        myfile << "f " << subTriangles[element][0]+1 << " " << subTriangles[element][1]+1 << " " << subTriangles[element][2]+1 << "\n";
    }
    myfile.close();
    msg_info() << "Mesh written in mesh_refined.obj";
}

// --------------------------------------------------------------------------------------
// Subdivides each triangle into 4 by taking the middle of each edge
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::subdivide(const Vec3& a, const Vec3& b, const Vec3& c, sofa::type::vector<Vec3> &subVertices, SeqTriangles &subTriangles)
{
    // Global coordinates
    Vec3 mAB, mAC, mBC;
    mAB = (a+b)/2;
    mAC = (a+c)/2;
    mBC = (b+c)/2;

    // Adds vertex if we deal with a new point
    int indexAB, indexAC, indexBC;
    addVertexAndFindIndex(subVertices, mAB, indexAB);
    addVertexAndFindIndex(subVertices, mAC, indexAC);
    addVertexAndFindIndex(subVertices, mBC, indexBC);

    // Finds index of the 3 original vertices
    int indexA, indexB, indexC;
    addVertexAndFindIndex(subVertices, a, indexA);
    addVertexAndFindIndex(subVertices, b, indexB);
    addVertexAndFindIndex(subVertices, c, indexC);

    // Adds the 4 subdivided triangles to the list
    subTriangles.push_back(Triangle(indexA, indexAB, indexAC));
    subTriangles.push_back(Triangle(indexAB, indexB, indexBC));
    subTriangles.push_back(Triangle(indexAC, indexBC, indexC));
    subTriangles.push_back(Triangle(indexBC, indexAC, indexAB));
}


// --------------------------------------------------------------------------------------
// Adds a vertex if it is not already in the list
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addVertexAndFindIndex(sofa::type::vector<Vec3> &subVertices, const Vec3 &vertex, int &index)
{
    bool alreadyHere = false;

    for (unsigned int v=0; v<subVertices.size(); v++)
    {
        if ( (subVertices[v]-vertex).norm() < 0.0000001)
        {
            alreadyHere = true;
            index = v;
        }
    }
    if (alreadyHere == false)
    {
        subVertices.push_back(vertex);
        index = (int)subVertices.size()-1;
    }
}


template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::movePoint(Vec3& pointToMove)
{
    sofa::type::vector<Vec3> listClosestPoints;
    findClosestGravityPoints(pointToMove, listClosestPoints);
    pointToMove = (listClosestPoints[0]+listClosestPoints[1]+listClosestPoints[2])/3;
}


// --------------------------------------------------------------------------------------
// Finds the list of the 3 closest gravity points of targeted surface
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::findClosestGravityPoints(const Vec3& point, sofa::type::vector<Vec3>& listClosestPoints)
{
    std::multimap<Real, Vec3> closestTrianglesData;

    for (unsigned int t=0; t<m_targetTriangles.size(); t++)
    {
        Vec3 pointTriangle1 = m_targetVertices[ m_targetTriangles[t][0] ];
        Vec3 pointTriangle2 = m_targetVertices[ m_targetTriangles[t][1] ];
        Vec3 pointTriangle3 = m_targetVertices[ m_targetTriangles[t][2] ];

        Vec3 G = (pointTriangle1+pointTriangle2+pointTriangle3)/3;

        // Distance between the point and current triangle
        Real distance = (G-point).norm2();

        // Stores distances (automatically sorted)
        std::make_pair<Real,Vec3 >(1.0,Vec3());
        closestTrianglesData.insert( std::make_pair(distance,G));
    }

    // Returns the 3 closest points
    int count = 0;
    typename std::multimap<Real,Vec3>::iterator it;
    for (it = closestTrianglesData.begin(); it!=closestTrianglesData.end(); it++)
    {
        if (count < 3)
        {
            listClosestPoints.push_back((*it).second);
        }
        count++;
    }

}


// --------------------------------------------------------------------------------------
// --- Re-initialization (called when we change a parameter through the GUI)
// --------------------------------------------------------------------------------------
template <class DataTypes>void TriangularBendingFEMForceField<DataTypes>::reinit()
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

    if (l_topologyMapper.get() != nullptr) {

        sofa::component::engine::JoinMeshPoints<DataTypes>* jmp = l_topologyMapper.get();
        if (jmp->f_output_triangles.getValue().size() == 0)
        {
            msg_warning() << "Mapped topology must be triangular. No triangles found." ;
        } else {
            m_mapTopology = true;
        }
    }

    if (l_restShape.get() != nullptr) {
        // Listen for MeshChangedEvent
        *this->f_listening.beginEdit() = true;
        this->f_listening.endEdit();

        // Check if there is same number of nodes
        if (!m_mapTopology) {
            if (l_restShape.get()->f_position.getValue().size() !=
                this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue().size()) {
                msg_warning() << "Different number of nodes in rest shape and mechanical state." ;
            }
        } else if (l_restShape.get()->f_position.getValue().size() !=
            l_topologyMapper.get()->f_input_position.getValue().size()) {
            msg_warning() << "Different number of nodes in rest shape and (original) mapped topology." ;
        }
    }

    /// Prepare to store info in the triangle array
    triangleInf.resize(_topology->getNbTriangles());

    for (sofa::Index i=0; i<_topology->getNbTriangles(); ++i)
    {
        m_triangleHandler->applyCreateFunction(i, triangleInf[i],  _topology->getTriangle(i),  (const sofa::type::vector< unsigned int > )0, (const sofa::type::vector< double >)0);
    }

    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
double TriangularBendingFEMForceField<DataTypes>::getPotentialEnergy(const VecCoord& /*x*/) const
{
    msg_warning() << "TriangularBendingFEMForceField::getPotentialEnergy is not implemented.";
    return 0;
}


// --------------------------------------------------------------------------------------
// Computes the quaternion that embodies the rotation from triangle to world
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeRotation(Quat& Qframe, const VecCoord &x, const Index &a, const Index &b, const Index &c)
{
    // First vector on first edge
    // Second vector in the plane of the two first edges
    // Third vector orthogonal to first and second

    Vec3 edgex = x[b].getCenter() - x[a].getCenter();
    edgex.normalize();

    Vec3 edgey = x[c].getCenter() - x[a].getCenter();
    edgey.normalize();

    Vec3 edgez;
    edgez = cross(edgex, edgey);
    edgez.normalize();

    edgey = cross(edgez, edgex);
    edgey.normalize();

    Transformation R;
    R[0][0] = edgex[0];
    R[0][1] = edgex[1];
    R[0][2] = edgex[2];
    R[1][0] = edgey[0];
    R[1][1] = edgey[1];
    R[1][2] = edgey[2];
    R[2][0] = edgez[0];
    R[2][1] = edgez[1];
    R[2][2] = edgez[2];

    Qframe.fromMatrix(R);
}


// --------------------------------------------------------------------------------------
// --- Initialization of a triangle, is called *only* once
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::initTriangleOnce(const int i, const Index&a, const Index&b, const Index&c)
{

    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[i];

    // Store indices of each vertex
    tinfo->a = a;
    tinfo->b = b;
    tinfo->c = c;

    Index a0=a, b0=b, c0=c;
    if (m_mapTopology) {
        sofa::component::engine::JoinMeshPoints<DataTypes>* jmp = l_topologyMapper.get();

        // Get indices in original topology
        a0 = jmp->getSrcNodeFromTri(i, a0);
        b0 = jmp->getSrcNodeFromTri(i, b0);
        c0 = jmp->getSrcNodeFromTri(i, c0);
    }

    tinfo->a0 = a0;
    tinfo->b0 = b0;
    tinfo->c0 = c0;

    triangleInfo.endEdit();
}

// --------------------------------------------------------------------------------------
// --- Initialization of a triangle. Can be called more than once if there is a
// --- changing rest shape.
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::initTriangle(const int i)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[i];

    const Index a = tinfo->a;
    const Index b = tinfo->b;
    const Index c = tinfo->c;

    const Index a0 = tinfo->a0;
    const Index b0 = tinfo->b0;
    const Index c0 = tinfo->c0;

    // Gets vertices of rest and initial positions respectively
    const VecCoord& x0 = (l_restShape.get() != nullptr)
        // if having changing rest shape take it
        ? l_restShape.get()->f_position.getValue()
        : (m_mapTopology
            // if rest shape is fixed but we have mapped topology use it
            ? l_topologyMapper.get()->f_input_position.getValue()
            // otherwise just take rest shape in mechanical state
            : this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue()
          );
    const VecCoord& x = this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();

    // Rotation from triangle to world at rest and initial positions (respectively)
    Quat Qframe0, Qframe;
    computeRotation(Qframe0, x0, a0, b0, c0 );
    computeRotation(Qframe, x, a, b, c );
    tinfo->Qframe = Qframe;

    // The positions of each point is expressed into the local frame at rest position
    tinfo->restLocalPositions[0] = Qframe0.rotate(x0[b0].getCenter() - x0[a0].getCenter());
    tinfo->restLocalPositions[1] = Qframe0.rotate(x0[c0].getCenter() - x0[a0].getCenter());

    if (d_bending.getValue())
    {
        // Computes inverse of C for initial position (in case of the latter is different than the rest_position)
        tinfo->localB = Qframe.rotate(x[b].getCenter()-x[a].getCenter());
        tinfo->localC = Qframe.rotate(x[c].getCenter()-x[a].getCenter());
        computeStrainDisplacementMatrixBending(tinfo, tinfo->localB, tinfo->localC);

        // Computes triangles' surface
        StrainDisplacement J;
        computeStrainDisplacementMatrix(J, i, tinfo->localB, tinfo->localC);

        // Local rest orientations (Evaluates the difference between the rest position and the flat position to allow the use of a deformed rest shape)
        tinfo->restLocalOrientations[0] = qDiffZ(x0[a0].getOrientation(), Qframe0);
        tinfo->restLocalOrientations[1] = qDiffZ(x0[b0].getOrientation(), Qframe0);
        tinfo->restLocalOrientations[2] = qDiffZ(x0[c0].getOrientation(), Qframe0);

        // Creates a vector u_rest matching the difference between rest and flat positions
        tinfo->u_rest.clear();
        tinfo->u_rest[1] = tinfo->restLocalOrientations[0].toEulerVector()[0];
        tinfo->u_rest[2] = tinfo->restLocalOrientations[0].toEulerVector()[1];
        tinfo->u_rest[4] = tinfo->restLocalOrientations[1].toEulerVector()[0];
        tinfo->u_rest[5] = tinfo->restLocalOrientations[1].toEulerVector()[1];
        tinfo->u_rest[7] = tinfo->restLocalOrientations[2].toEulerVector()[0];
        tinfo->u_rest[8] = tinfo->restLocalOrientations[2].toEulerVector()[1];

        // Computes vector displacement between initial position and rest positions (actual displacements that define the amount of stress within the structure)
        DisplacementBending Disp_bending;
        computeDisplacementBending(Disp_bending, x, i);

    }

    triangleInfo.endEdit();
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::applyStiffness(VecDeriv& v, const VecDeriv& dx, const Index elementIndex, const double kFactor)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Get the indices of the 3 vertices for the current triangle
    const Index& a = tinfo->a;
    const Index& b = tinfo->b;
    const Index& c = tinfo->c;

    // Computes displacements
    Displacement Disp;
    Vec3 x_a, x_b, x_c;
    x_a = tinfo->Qframe.rotate(getVCenter(dx[a]));
    Disp[0] = x_a[0];
    Disp[1] = x_a[1];

    x_b = tinfo->Qframe.rotate(getVCenter(dx[b]));
    Disp[2] = x_b[0];
    Disp[3] = x_b[1];

    x_c = tinfo->Qframe.rotate(getVCenter(dx[c]));
    Disp[4] = x_c[0];
    Disp[5] = x_c[1];

    // Compute dF
    Displacement dF;
    dF = tinfo->stiffnessMatrix * Disp;

    // Transfer into global frame
    getVCenter(v[a]) += tinfo->Qframe.inverseRotate(Vec3(-dF[0], -dF[1], 0)) * kFactor;
    getVCenter(v[b]) += tinfo->Qframe.inverseRotate(Vec3(-dF[2], -dF[3], 0)) * kFactor;
    getVCenter(v[c]) += tinfo->Qframe.inverseRotate(Vec3(-dF[4], -dF[5], 0)) * kFactor;

    // If bending is requested
    if (d_bending.getValue())
    {
        // Bending displacements
        DisplacementBending Disp_bending;
        Vec3 u;
        u = tinfo->Qframe.rotate(getVOrientation(dx[a]));
        Disp_bending[0] = x_a[2];
        Disp_bending[1] = u[0];
        Disp_bending[2] = u[1];

        u = tinfo->Qframe.rotate(getVOrientation(dx[b]));
        Disp_bending[3] = x_b[2];
        Disp_bending[4] = u[0];
        Disp_bending[5] = u[1];

        u = tinfo->Qframe.rotate(getVOrientation(dx[c]));
        Disp_bending[6] = x_c[2];
        Disp_bending[7] = u[0];
        Disp_bending[8] = u[1];

        // Compute dF
        DisplacementBending dF_bending;
        dF_bending = tinfo->stiffnessMatrixBending * Disp_bending;

        // Go back into global frame
        Vec3 fa1, fa2, fb1, fb2, fc1, fc2;
        fa1 = tinfo->Qframe.inverseRotate(Vec3(0.0, 0.0, dF_bending[0]));
        fa2 = tinfo->Qframe.inverseRotate(Vec3(dF_bending[1], dF_bending[2], 0.0));

        fb1 = tinfo->Qframe.inverseRotate(Vec3(0.0, 0.0, dF_bending[3]));
        fb2 = tinfo->Qframe.inverseRotate(Vec3(dF_bending[4], dF_bending[5], 0.0));

        fc1 = tinfo->Qframe.inverseRotate(Vec3(0.0, 0.0, dF_bending[6]));
        fc2 = tinfo->Qframe.inverseRotate(Vec3(dF_bending[7], dF_bending[8], 0.0));

        v[a] += Deriv(-fa1, -fa2) * kFactor;
        v[b] += Deriv(-fb1, -fb2) * kFactor;
        v[c] += Deriv(-fc1, -fc2) * kFactor;
    }


    triangleInfo.endEdit();
}

// -------------------------------------------------------------------------------------------------------------
// --- Compute displacement vector D as the difference between current current position 'p' and initial position
// --- expressed in the co-rotational frame of reference
// -------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeDisplacement(Displacement &Disp, const VecCoord &x, const Index elementIndex)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    Index a = tinfo->a;
    Index b = tinfo->b;
    Index c = tinfo->c;

    // Compute local postions of vertices b and c in the co-rotational frame (a is always (0,0,0))
    tinfo->localB = tinfo->Qframe.rotate(x[b].getCenter()-x[a].getCenter());
    tinfo->localC = tinfo->Qframe.rotate(x[c].getCenter()-x[a].getCenter());

    // In-plane local displacements
    Vec3 uAB, uAC;
    uAB = tinfo->localB - tinfo->restLocalPositions[0];
    uAC = tinfo->localC - tinfo->restLocalPositions[1];

    Disp[0] = 0;
    Disp[1] = 0;
    Disp[2] = uAB[0];
    Disp[3] = 0;
    Disp[4] = uAC[0];
    Disp[5] = uAC[1];

    triangleInfo.endEdit();
}


// -------------------------------------------------------------------------------------------------------------
// --- Compute bending displacement vector D as the difference between current current position 'p' and initial position
// --- expressed in the co-rotational frame of reference
// -------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeDisplacementBending(DisplacementBending &Disp, const VecCoord &x, const Index elementIndex)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    Index a = tinfo->a;
    Index b = tinfo->b;
    Index c = tinfo->c;

    Quat dQA = qDiffZ(x[a].getOrientation(), tinfo->Qframe);     // Rotation of axis Z from Qframe to axis Z from x[a].getOrientation()
    Quat dQB = qDiffZ(x[b].getOrientation(), tinfo->Qframe);     // Rotation of axis Z from Qframe to axis Z from x[b].getOrientation()
    Quat dQC = qDiffZ(x[c].getOrientation(), tinfo->Qframe);     // Rotation of axis Z from Qframe to axis Z from x[c].getOrientation()

    // Difference between the current and the rest orientations yields displacement (into the triangle's frame)
    dQA = qDiff(tinfo->restLocalOrientations[0].inverse(), dQA.inverse());
    dQB = qDiff(tinfo->restLocalOrientations[1].inverse(), dQB.inverse());
    dQC = qDiff(tinfo->restLocalOrientations[2].inverse(), dQC.inverse());

    // Takes the Euler vector to get the rotation's axis
    Vec3 rA, rB, rC;
    rA = dQA.toEulerVector();
    rB = dQB.toEulerVector();
    rC = dQC.toEulerVector();

    // Writes the computed displacements
    Disp[0] = 0;          // z displacement in A
    Disp[1] = rA[0];      // x rotation in A
    Disp[2] = rA[1];      // y rotation in A

    Disp[3]  = 0;         // z displacement in B
    Disp[4] = rB[0];     // x rotation in B
    Disp[5] = rB[1];     // y rotation in B

    Disp[6] = 0;         // z displacement in C
    Disp[7] = rC[0];     // x rotation in C
    Disp[8] = rC[1];     // y rotation in C

    // Stores the vector u of displacements (used by the mechanical mapping for rendering)
    tinfo->u = Disp;

    triangleInfo.endEdit();
}

// ------------------------------------------------------------------------------------------------------------
// --- Compute the strain-displacement matrix where (a, b, c) are the local coordinates of the 3 nodes of a triangle
// ------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeStrainDisplacementMatrix(StrainDisplacement &J, const Index elementIndex, const Vec3& b, const Vec3& c)
{
    Real determinant;
    determinant = b[0] * c[1];

    Real x13 = -c[0] / determinant; // since a=(0,0)
    Real x21 = b[0] / determinant; // since a=(0,0)
    Real x32 = (c[0]-b[0]) / determinant;
    Real y12 = 0;       // since a=(0,0) and b[1] = 0
    Real y23 = -c[1] / determinant; // since a=(0,0) and b[1] = 0
    Real y31 = c[1] / determinant; // since a=(0,0)

    J[0][0] = y23;
    J[0][1] = 0;
    J[0][2] = x32;

    J[1][0] = 0;
    J[1][1] = x32;
    J[1][2] = y23;

    J[2][0] = y31;
    J[2][1] = 0;
    J[2][2] = x13;

    J[3][0] = 0;
    J[3][1] = x13;
    J[3][2] = y31;

    J[4][0] = y12;
    J[4][1] = 0;
    J[4][2] = x21;

    J[5][0] = 0;
    J[5][1] = x21;
    J[5][2] = y12;

    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    triangleInf[elementIndex].area = 0.5*determinant;
    triangleInfo.endEdit();

}


// ------------------------------------------------------------------------------------------------------------
// --- Compute the bending strain-displacement matrix where (a, b, c) are the coordinates of the 3 nodes of a triangle
// ------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeStrainDisplacementMatrixBending(TriangleInformation *tinfo, const Vec3& b, const Vec3& c)
{
    // Computation of the inverse of matrix C. Its inverse gives the coefficients c1, c2, ..., c9 of the deflection function given by:
    // Uz = c1 + c2*x+ c3*y + c4*x^2 + c5*x*y + c6*y^2 + c7*x^3 + c8*(x*y^2 + x^2*y) + c9*y^3
    // Source: Tocher's deflection function presented by Przemieniecki
    // Corrected deflection to get a symmetrical deformation: Uz = c1 + c2*x+ c3*y + c4*x^2 + c5*x*y + c6*y^2 + c7*x^3 + c8*x*y^2 + c9*y^3

    Mat<9, 9, Real> C;
    C.clear();

    // Corrected
    C(0,0) = 1;
    C(1,2) = 1;
    C(2,1) = -1;
    C(3,0) = 1;         C(3,1) = b[0];          C(3,3) = b[0]*b[0];         C(3,6) = b[0]*b[0]*b[0];
    C(4,2) = 1;         C(4,4) = b[0];
    C(5,1) = -1;        C(5,3) = -2*b[0];       C(5,6) = -3*b[0]*b[0];
    C(6,0) = 1;         C(6,1) = c[0];          C(6,2) = c[1];              C(6,3) = c[0]*c[0];             C(6,4) = c[0]*c[1];         C(6,5) = c[1]*c[1];
    C(6,6) = c[0]*c[0]*c[0];                    C(6,7) = c[0]*c[1]*c[1];    C(6,8) = c[1]*c[1]*c[1];
    C(7,2) = 1;         C(7,4) = c[0];          C(7,5) = 2*c[1];            C(7,7) = 2*c[0]*c[1];           C(7,8) = 3*c[1]*c[1];
    C(8,1) = -1;        C(8,3) = -2*c[0];       C(8,4) = -c[1];             C(8,6) = -3*c[0]*c[0];          C(8,7) = -c[1]*c[1];

    // Inverse of C
    Mat<9, 9, Real> invC;
    invC.invert(C);
    tinfo->invC = invC;

    // Calculation of the 3 Gauss points taken in the centre of each edge
    Vec3 gaussPoint1 = b*0.5;       // a is (0, 0, 0)
    Vec3 gaussPoint2 = (b+c)*0.5;
    Vec3 gaussPoint3 = c*0.5;       // a is (0, 0, 0)

    // Retrieves the strain tensor used in flat-plate theory at each Gauss point
    Mat<3, 9, Real> D1, D2, D3;
    tensorFlatPlate(D1, gaussPoint1);
    tensorFlatPlate(D2, gaussPoint2);
    tensorFlatPlate(D3, gaussPoint3);

    // Compute strain-displacement matrix
    tinfo->strainDisplacementMatrix1 = D1 * invC;
    tinfo->strainDisplacementMatrix2 = D2 * invC;
    tinfo->strainDisplacementMatrix3 = D3 * invC;
}


// --------------------------------------------------------------------------------------------------------
// --- Computes the strain tensor used in flat-plate theory in a given point
// --------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::tensorFlatPlate(Mat<3, 9, Real>& D, const Vec3 &P)
{
#ifdef DEBUG_TRIANGLEFEM
    msg_info() << "TriangleBendingFEMForceField::tensorFlatPlate";
#endif

    // Flat-plat theory gives:
    // e = D * c with
    //        [ 0  0  0  2  0  0  6x   2y   0  ]
    // D = -z | 0  0  0  0  0  2  0    2x   6y |
    //        [ 0  0  0  0  2  0  0  4(x+y) 0  ]
    // where e is the strain vector and c the coefficient vector of the deflection function
    //
    // CORRECTED:
    //        [ 0  0  0  2  0  0  6x  0   0  ]
    // D = -z | 0  0  0  0  0  2  0   2x  6y |
    //        [ 0  0  0  0  2  0  0   4y  0  ]

    // Corrected
    D.clear();
    D(0,3) = 2;         D(0,6) = 6*P[0];
    D(1,5) = 2;         D(1,7) = 2*P[0];        D(1,8) = 6*P[1];
    D(2,4) = 2;         D(2,7) = 4*P[1];
}


// ----------------------------------------------------------------------------------------------------------------------
// --- Compute the stiffness matrix K = J * M * Jt where J is the strain-displacement matrix and M the material matrix
// ----------------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeStiffnessMatrix(StiffnessMatrix &K, const StrainDisplacement &J, const MaterialStiffness &M)
{
    Mat<3,6,Real> Jt;
    Jt.transpose(J);

    K = J * M * Jt;
}


// ----------------------------------------------------------------------------------------------------------------------
// --- Compute the stiffness matrix for bending K = J * M * Jt where J is the strain-displacement matrix and M the material matrix
// ----------------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeStiffnessMatrixBending(StiffnessMatrixBending &K, TriangleInformation *tinfo)
{
    Mat<9, 3, Real> J1t, J2t, J3t;
    J1t.transpose(tinfo->strainDisplacementMatrix1);
    J2t.transpose(tinfo->strainDisplacementMatrix2);
    J3t.transpose(tinfo->strainDisplacementMatrix3);

    // K = (surface/3)*(t^3)*(J1t*material*J1 + J2t*material*J2 + J3t*material*J3)
    K = J1t * tinfo->materialMatrix * tinfo->strainDisplacementMatrix1 +
        J2t * tinfo->materialMatrix * tinfo->strainDisplacementMatrix2 +
        J3t * tinfo->materialMatrix * tinfo->strainDisplacementMatrix3;
}

// --------------------------------------------------------------------------------------
// ---  Compute material stiffness
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeMaterialStiffness(const int i)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

    TriangleInformation *tinfo = &triangleInf[i];

    tinfo->materialMatrix[0][0] = 1;
    tinfo->materialMatrix[0][1] = d_poisson.getValue();
    tinfo->materialMatrix[0][2] = 0;
    tinfo->materialMatrix[1][0] = d_poisson.getValue();
    tinfo->materialMatrix[1][1] = 1;
    tinfo->materialMatrix[1][2] = 0;
    tinfo->materialMatrix[2][0] = 0;
    tinfo->materialMatrix[2][1] = 0;
    tinfo->materialMatrix[2][2] = 0.5f * (1 - d_poisson.getValue());

    tinfo->materialMatrix *= (d_young.getValue() / (12 *  (1 - d_poisson.getValue() * d_poisson.getValue())));

    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---  Compute force F = J * material * Jt * u
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeForce(Displacement &F, const Displacement& D, const Index elementIndex)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Compute strain-displacement matrix J
    StrainDisplacement J;
    computeStrainDisplacementMatrix(J, elementIndex, tinfo->localB, tinfo->localC);
    tinfo->strainDisplacementMatrix = J;

    // Compute stiffness matrix K = J*material*Jt
    StiffnessMatrix K;
    computeStiffnessMatrix(K, J, tinfo->materialMatrix);
    tinfo->stiffnessMatrix = K;

    // Compute forces
    F = K * D;

    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---  Compute force F = Jt * material * J * u
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeForceBending(DisplacementBending &F_bending, const DisplacementBending& D_bending, const Index elementIndex)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Compute strain-displacement matrix J
    computeStrainDisplacementMatrixBending(tinfo, tinfo->localB, tinfo->localC);

    // Compute stiffness matrix K = Jt * material * J
    StiffnessMatrixBending K_bending;
    computeStiffnessMatrixBending(K_bending, tinfo);
    Real t = d_thickness.getValue();
    K_bending *= t*t*t*(triangleInf[elementIndex].area)/3;
    tinfo->stiffnessMatrixBending = K_bending;

    // Compute forces
    F_bending = K_bending * D_bending;

    triangleInfo.endEdit();
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::accumulateForce(VecDeriv &f, const VecCoord &x, const Index elementIndex)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Get the indices of the 3 vertices for the current triangle
    const Index& a = tinfo->a;
    const Index& b = tinfo->b;
    const Index& c = tinfo->c;

    // Compute the quaternion that embodies the rotation between the triangle and world frames (co-rotational method)
    Quat Qframe;
    computeRotation(Qframe, x, a, b, c);
    tinfo->Qframe = Qframe;

    // Compute in-plane displacement into the triangle's frame
    Displacement D;
    computeDisplacement(D, x, elementIndex);

    // Compute in-plane forces on this element (in the co-rotational space)
    Displacement F;
    computeForce(F, D, elementIndex);

    // Transform forces back into global reference frame
    getVCenter(f[a]) -= tinfo->Qframe.inverseRotate(Vec3(F[0], F[1], 0));
    getVCenter(f[b]) -= tinfo->Qframe.inverseRotate(Vec3(F[2], F[3], 0));
    getVCenter(f[c]) -= tinfo->Qframe.inverseRotate(Vec3(F[4], F[5], 0));

    if (d_bending.getValue())
    {
        // Compute bending displacement for bending into the triangle's frame
        DisplacementBending D_bending;
        computeDisplacementBending(D_bending, x, elementIndex);

        // Compute bending forces on this element (in the co-rotational space)
        DisplacementBending F_bending;
        computeForceBending(F_bending, D_bending, elementIndex);

        // Transform forces back into global reference frame
        Vec3 fa1 = tinfo->Qframe.inverseRotate(Vec3(0.0, 0.0, F_bending[0]));
        Vec3 fa2 = tinfo->Qframe.inverseRotate(Vec3(F_bending[1], F_bending[2], 0.0));

        Vec3 fb1 = tinfo->Qframe.inverseRotate(Vec3(0.0, 0.0, F_bending[3]));
        Vec3 fb2 = tinfo->Qframe.inverseRotate(Vec3(F_bending[4], F_bending[5], 0.0));

        Vec3 fc1 = tinfo->Qframe.inverseRotate(Vec3(0.0, 0.0, F_bending[6]));
        Vec3 fc2 = tinfo->Qframe.inverseRotate(Vec3(F_bending[7], F_bending[8], 0.0));

        f[a] += Deriv(-fa1, -fa2);
        f[b] += Deriv(-fb1, -fb2);
        f[c] += Deriv(-fc1, -fc2);
    }

    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& dataF, const DataVecCoord& dataX, const DataVecDeriv& /*dataV*/ )
{
    VecDeriv& f        = *(dataF.beginEdit());
    const VecCoord& p  =   dataX.getValue()  ;

    int nbTriangles=_topology->getNbTriangles();
    f.resize(p.size());

    for (int i=0; i<nbTriangles; i++)
    {
        accumulateForce(f, p, i);
    }

    dataF.endEdit();
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& datadF, const DataVecDeriv& datadX )
{
    VecDeriv& df        = *(datadF.beginEdit());
    const VecDeriv& dp  =   datadX.getValue()  ;

    double kFactor = mparams->kFactor();

    int nbTriangles=_topology->getNbTriangles();
    df.resize(dp.size());

    for (int i=0; i<nbTriangles; i++)
    {
        applyStiffness(df, dp, i, kFactor);
    }

    datadF.endEdit();
}


template<class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::convertStiffnessMatrixToGlobalSpace(StiffnessMatrixGlobalSpace &K_gs, TriangleInformation *tinfo)
{
    // Stiffness matrix of current triangle
    const StiffnessMatrix &K = tinfo->stiffnessMatrix;

    // Firstly, add all degrees of freedom (we add the unused translation in z)
    StiffnessMatrixGlobalSpace K_18x18;
    K_18x18.clear();
    unsigned int ig, jg;

    // Copy the stiffness matrix by block 2x2 into 18x18 matrix (the new index of each bloc into global matrix is a combination of 0, 6 and 12 in indices)
    for (unsigned int bx=0; bx<3; bx++)
    {
        // Global row index
        ig = 6*bx;

        for (unsigned int by=0; by<3; by++)
        {
            // Global column index
            jg = 6*by;

            // Iterates over the indices of the bloc 2x2
            for (unsigned int i=0; i<2; i++)
            {
                for (unsigned int j=0; j<2; j++)
                {
                    K_18x18[ig+i][jg+j] = K[2*bx+i][2*by+j];
                }
            }
        }
    }


    if (d_bending.getValue())
    {
        // Stiffness matrix in bending of current triangle
        const StiffnessMatrixBending &K_bending = tinfo->stiffnessMatrixBending;

        // Copy the stiffness matrix by block 3x3 into global matrix (the new index of each bloc into global matrix is a combination of 2, 8 and 15 in indices)
        for (unsigned int bx=0; bx<3; bx++)
        {
            // Global row index
            ig = 6*bx+2;

            for (unsigned int by=0; by<3; by++)
            {
                // Global column index
                jg = 6*by+2;

                // Iterates over the indices of the bloc 3x3
                for (unsigned int i=0; i<3; i++)
                {
                    for (unsigned int j=0; j<3; j++)
                    {
                        K_18x18[ig+i][jg+j] += K_bending[3*bx+i][3*by+j];
                    }
                }

            }
        }

    }

    // Extend rotation matrix and its transpose
    Transformation R, Rt;
    tinfo->Qframe.toMatrix(R);
    Rt.transpose(R);

    StiffnessMatrixGlobalSpace R18x18, Rt18x18;

    for(unsigned int i=0;i<3;++i)
    {
        for(unsigned int j=0;j<3;++j)
        {
            R18x18[i][j] = R18x18[i+3][j+3] = R18x18[i+6][j+6] = R18x18[i+9][j+9] = R18x18[i+12][j+12] = R18x18[i+15][j+15] = R[i][j];
            Rt18x18[i][j] = Rt18x18[i+3][j+3] = Rt18x18[i+6][j+6] = Rt18x18[i+9][j+9] = Rt18x18[i+12][j+12] = Rt18x18[i+15][j+15] = Rt[i][j];
        }
    }

    // Then we put the stifness matrix into the global frame
    K_gs = Rt18x18 * K_18x18 * R18x18;

}

#define ASSEMBLED_K

#ifdef ASSEMBLED_K

template<class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
    StiffnessMatrixGlobalSpace K_gs;

    // Build Matrix Block for this ForceField
    unsigned int i, j ,n1, n2, row, column, ROW, COLUMN;
    Index node1, node2;

    sofa::core::behavior::MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

    double kFactor = mparams->kFactor();

    for(sofa::Index t=0 ; t != _topology->getNbTriangles() ; ++t)
    {
            TriangleInformation *tinfo = &triangleInf[t];
            const Triangle triangle = _topology->getTriangle(t);

            convertStiffnessMatrixToGlobalSpace(K_gs, tinfo);

            // find index of node 1
            for (n1=0; n1<3; n1++)
            {
                    node1 = triangle[n1];

                    for(i=0; i<6; i++)
                    {
                            ROW = r.offset+6*node1+i;
                            row = 6*n1+i;
                            // find index of node 2
                            for (n2=0; n2<3; n2++)
                            {
                                    node2 = triangle[n2];

                                    for (j=0; j<6; j++)
                                    {
                                            COLUMN = r.offset+6*node2+j;
                                            column = 6*n2+j;
                                            r.matrix->add(ROW, COLUMN, - K_gs[row][column] * kFactor);
                                    }
                            }
                    }
            }
    }

    triangleInfo.endEdit();
}


#else

template<class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addKToMatrix(sofa::linearalgebra::BaseMatrix *mat, SReal /*k*/, unsigned int &offset)
{
    VecCoord X = *this->mstate->getX();
    VecDeriv df, dx;

    dx.resize(X.size());
    df.resize(X.size());

    for (unsigned int i=0; i<X.size(); i++)
    {
        for (unsigned int j=0; j<6; j++)
        {
            dx.clear();
            df.clear();
            dx.resize(X.size());
            df.resize(X.size());
            dx[i][j] = 1;
            addDForce(df, dx);

            for (unsigned int k=0; k<X.size(); k++)
            {
                for (unsigned int l=0; l<6; l++)
                {
                    mat->add(6*k+l, 6*i+j, -df[k][l]);
                }
            }

        }
    }
}

#endif


template<class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addBToMatrix(sofa::linearalgebra::BaseMatrix * /*mat*/, double /*bFact*/, unsigned int &/*offset*/)
{
}


template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowForceFields())
        return;

    const auto stateLifeCycle = vparams->drawTool()->makeStateLifeCycle();
    vparams->drawTool()->disableLighting();

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0, true);

    std::vector<sofa::type::RGBAColor> colorVector;
    std::vector<sofa::type::Vec3> vertices;

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();

    const SeqTriangles triangles = _topology->getTriangles();
    for (const Triangle& tri: triangles)
    {
        colorVector.push_back(sofa::type::RGBAColor::green());
        vertices.push_back(sofa::type::Vec3(DataTypes::getCPos(x[tri[0]])));
        colorVector.push_back(sofa::type::RGBAColor(0, 0.5, 0.5, 1));
        vertices.push_back(sofa::type::Vec3(DataTypes::getCPos(x[tri[1]])));
        colorVector.push_back(sofa::type::RGBAColor(0, 0, 1, 1));
        vertices.push_back(sofa::type::Vec3(DataTypes::getCPos(x[tri[2]])));
    }
    vparams->drawTool()->drawTriangles(vertices, colorVector);

    if (m_mapTopology){
        // Draw lines to visualize the mapping between nodes
        std::vector<Vec<3,double> > points;

        // Gets vertices of rest and initial positions respectively
        const VecCoord& x0 = (l_restShape.get() != nullptr)
            // if having changing rest shape take it
            ? l_restShape.get()->f_position.getValue()
            : (m_mapTopology
                // if rest shape is fixed but we have mapped topology use it
                ? l_topologyMapper.get()->f_input_position.getValue()
                // otherwise just take rest shape in mechanical state
                : this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue()
              );
        const VecCoord& x = this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();

        int nbTriangles=_topology->getNbTriangles();

        for (sofa::Index i=0; i<nbTriangles; i++) {
            const TriangleInformation &tinfo = triangleInfo.getValue()[i];
            if ((x[tinfo.a].getCenter() - x0[tinfo.a0].getCenter()).norm() > 1e-8) {
                points.push_back(x[tinfo.a].getCenter());
                points.push_back(x0[tinfo.a0].getCenter());
            }
            if ((x[tinfo.b].getCenter() - x0[tinfo.b0].getCenter()).norm() > 1e-8) {
                points.push_back(x[tinfo.b].getCenter());
                points.push_back(x0[tinfo.b0].getCenter());
            }
            if ((x[tinfo.c].getCenter() - x0[tinfo.c0].getCenter()).norm() > 1e-8) {
                points.push_back(x[tinfo.c].getCenter());
                points.push_back(x0[tinfo.c0].getCenter());
            }
        }

        if (points.size() > 0)
            vparams->drawTool()->drawLines(points, 1.0f, type::RGBAColor(1.0, 0.0, 1.0, 1.0));
    }
}

template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::handleEvent(sofa::core::objectmodel::Event *event)
{
    if (dynamic_cast<sofa::core::objectmodel::MeshChangedEvent*>(event))
    {
        // Update of the rest shape
        // NOTE: the number of triangles should be the same in all topologies
        unsigned int nbTriangles = _topology->getNbTriangles();
        for (unsigned int i=0; i<nbTriangles; i++) {
            initTriangle(i);
        }
    }
}

} // namespace forcefield

} // namespace component

} // namespace sofa


#endif
