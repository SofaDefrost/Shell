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

#include <SofaShells/forcefield/TriangularShellForceField.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/ColorMap.h>
#include <sofa/type/RGBAColor.h>

#include <sofa/core/topology/TopologyData.inl>
#include <sofa/gl/template.h>
#include <sofa/helper/rmath.h>
#include <sofa/gl/gl.h>
#include <sofa/gl/template.h>
#include <sofa/helper/system/thread/debug.h>
#include <fstream> // for reading the file
#include <iostream> //for debugging
#include <vector>
#include <algorithm>
#include <sofa/defaulttype/VecTypes.h>
#include <assert.h>
#include <map>
#include <utility>

#include <sofa/helper/ScopedAdvancedTimer.h>

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
using namespace	sofa::core::topology;

// --------------------------------------------------------------------------------------
// ---  Topology Creation/Destruction functions
// --------------------------------------------------------------------------------------
template< class DataTypes>
void TriangularShellForceField<DataTypes>::TRQSTriangleHandler::applyCreateFunction(unsigned int triangleIndex, TriangleInformation &, const Triangle &t, const sofa::type::vector<unsigned int> &, const sofa::type::vector<double> &, const VecCoord& x0)
{
    if (ff)
        ff->initTriangle(triangleIndex, t[0], t[1], t[2], x0);
}


// --------------------------------------------------------------------------------------
// --- Constructor
// --------------------------------------------------------------------------------------
template <class DataTypes>
TriangularShellForceField<DataTypes>::TriangularShellForceField()
    : d_poisson(initData(&d_poisson,(Real)0.45,"poissonRatio","Poisson ratio in Hooke's law"))
    , d_young(initData(&d_young,(Real)3000.,"youngModulus","Young modulus in Hooke's law"))
    , d_thickness(initData(&d_thickness,(Real)0.1,"thickness","Thickness of the plates"))
    , d_membraneElement(initData(&d_membraneElement, "membraneElement", "The membrane element to use"))
    , d_bendingElement(initData(&d_bendingElement, "bendingElement", "The bending plate element to use"))
    , d_corotated(initData(&d_corotated, true, "corotated", "Compute forces in corotational frame"))
    , d_measure(initData(&d_measure, "measure", "Compute the strain or stress"))
    , d_measuredValues(initData(&d_measuredValues, "measuredValues", "Measured values for stress or strain"))
    , d_isShellveryThin(initData(&d_isShellveryThin, false, "isShellveryThin", "This bool is to adapt "
                                                                               "computation in case we are using verry tiny(thickness) shell element"))
    , d_use_rest_position(initData(&d_use_rest_position, true, "use_rest_position", "Use the rest position inteat of using postion to update the restposition"))
    , triangleInfo(initData(&triangleInfo, "triangleInfo", "Internal triangle data"))
    , d_showArrow(initData(&d_showArrow, false, "showArrow", "Show the orientation of the triangle as arrow"))
    , d_arrowRadius(initData(&d_arrowRadius, (Real)0.1, "arrow_radius", "the arrow radius"))
    , d_showTriangle(initData(&d_showTriangle, false, "showTriangle", "Show the triangle"))
    , d_showMeasuredValue(initData(&d_showMeasuredValue, false, "showMeasuredValue", "Show the measured value"))
{
    d_membraneElement.beginEdit()->setNames(7,
                                            "None",     // No membrane element
                                            "CST",      // Constant strain triangle
                                            // ANDES templates
                                            "ALL-3I",   // Allman 88 element integrated by 3-point interior rule
                                            "ALL-3M",   // Allman 88 element integrated by 3-midpoint rule
                                            "ALL-LS",   // Allman 88 element, least-square strain fit
                                            "LST-Ret",  // Retrofitted LST with α_b=1⁄4
                                            "ANDES-OPT" // Optimal ANDES element
                                            );
    d_membraneElement.beginEdit()->setSelectedItem("ANDES-OPT");
    d_membraneElement.endEdit();

    d_bendingElement.beginEdit()->setNames(2,
                                           "None",     // No bending element
                                           "DKT"       // Discrete Kirchhoff Triangle
                                           );
    d_bendingElement.beginEdit()->setSelectedItem("DKT");
    d_bendingElement.endEdit();

    d_measure.beginEdit()->setNames(3,
                                    "None",                 // Draw nothing
                                    "Strain (norm)",        // L_2 norm of strain in x and y directions
                                    "Von Mises stress"      // Von Mises stress criterion
                                    );
    d_measure.beginEdit()->setSelectedItem("None");
    d_measure.endEdit();

    triangleHandler = new TRQSTriangleHandler(this, &triangleInfo);

    d_drawColorMap = helper::ColorMap(256, "Blue to Red");
}


// --------------------------------------------------------------------------------------
// --- Destructor
// --------------------------------------------------------------------------------------
template <class DataTypes>
TriangularShellForceField<DataTypes>::~TriangularShellForceField()
{
    if(triangleHandler) delete triangleHandler;
}

// --------------------------------------------------------------------------------------
// --- Initialization stage
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularShellForceField<DataTypes>::init()
{
    this->Inherited::init();

    _topology = this->getContext()->getMeshTopology();

    if (_topology->getNbTriangles()==0)
    {
        msg_warning() << "TriangularShellForceField: object must have a Triangular Set Topology.";
        return;
    }

    // Create specific handler for TriangleData
    triangleInfo.createTopologyHandler(_topology);

    reinit();
}


// --------------------------------------------------------------------------------------
// --- Re-initialization (called when we change a parameter through the GUI)
// --------------------------------------------------------------------------------------
template <class DataTypes> void TriangularShellForceField<DataTypes>::reinit()
{
    type::vector<TriangleInformation>& ti = *(triangleInfo.beginEdit());

    // Prepare material matrices
    computeMaterialStiffness();

    // Decode the selected elements to use
    if (d_membraneElement.getValue().getSelectedItem() == "None") {
        csMembrane = NULL;
        for (unsigned int t=0; t<ti.size(); ++t) {
            for (unsigned int i=0; i<ti[t].measure.size(); ++i) {
                ti[t].measure[i].B.clear();
            }
        }
    } else if (d_membraneElement.getValue().getSelectedItem() == "CST") {
        csMembrane = &TriangularShellForceField<DataTypes>::computeStiffnessMatrixCST;
    } else if (d_membraneElement.getValue().getSelectedItem() == "ALL-3I") {
        csMembrane = &TriangularShellForceField<DataTypes>::computeStiffnessMatrixAll3I;
    } else if (d_membraneElement.getValue().getSelectedItem() == "ALL-3M") {
        csMembrane = &TriangularShellForceField<DataTypes>::computeStiffnessMatrixAll3M;
    } else if (d_membraneElement.getValue().getSelectedItem() == "ALL-LS") {
        csMembrane = &TriangularShellForceField<DataTypes>::computeStiffnessMatrixAllLS;
    } else if (d_membraneElement.getValue().getSelectedItem() == "LST-Ret") {
        csMembrane = &TriangularShellForceField<DataTypes>::computeStiffnessMatrixLSTRet;
    } else if (d_membraneElement.getValue().getSelectedItem() == "ANDES-OPT") {
        csMembrane = &TriangularShellForceField<DataTypes>::computeStiffnessMatrixAndesOpt;
    } else {
        msg_warning() << "Invalid membrane element '" << d_membraneElement.getValue().getSelectedItem() << "'" ;
        return;
    }

    if (d_bendingElement.getValue().getSelectedItem() == "None") {
        csBending = NULL;
        for (unsigned int t=0; t<ti.size(); ++t) {
            for (unsigned int i=0; i<ti[t].measure.size(); ++i) {
                ti[t].measure[i].Bb.clear();
            }
        }
    } else if (d_bendingElement.getValue().getSelectedItem() == "DKT") {
        csBending = &TriangularShellForceField<DataTypes>::computeStiffnessMatrixDKT;
    } else {
        msg_warning() << "Invalid bending plate element '" << d_bendingElement.getValue().getSelectedItem() << "'" ;
        return;
    }

    // What to compute?
    if (d_measure.getValue().getSelectedItem() == "None") {
        bMeasureStrain = false;  bMeasureStress = false;
    } else if (d_measure.getValue().getSelectedItem() == "Strain (norm)") {
        bMeasureStrain = true;  bMeasureStress = false;
    } else if (d_measure.getValue().getSelectedItem() == "Von Mises stress") {
        bMeasureStrain = false;  bMeasureStress = true;
    } else {
        msg_warning() << "Invalid value for measure'" << d_measure.getValue().getSelectedItem() << "'" ;
        return;
    }

    if (bMeasureStrain || bMeasureStress)
    {
        d_measuredValues.beginEdit()->resize(_topology->getNbPoints());
        d_measuredValues.endEdit();
    }

    /// Prepare to store info in the triangle array
    ti.resize(_topology->getNbTriangles());

    const auto use_rest_position = d_use_rest_position.getValue();

    if (use_rest_position)
    {
        const VecCoord& x0 =this->mstate->read(sofa::core::ConstVecCoordId::restPosition())->getValue();
        for (sofa::Index i=0; i<_topology->getNbTriangles(); ++i)
            triangleHandler->applyCreateFunction(i, ti[i], _topology->getTriangle(i),
                                                 (const sofa::type::vector< unsigned int >)0, (const sofa::type::vector< double >)0, x0);
    }
    else
    {
        const VecCoord& x0 =this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();
        for (sofa::Index i=0; i<_topology->getNbTriangles(); ++i)
            triangleHandler->applyCreateFunction(i, ti[i], _topology->getTriangle(i),
                                                 (const sofa::type::vector< unsigned int >)0, (const sofa::type::vector< double >)0, x0);
    }
    triangleInfo.endEdit();
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularShellForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& dataF, const DataVecCoord& dataX, const DataVecDeriv& /*dataV*/ )
{
    VecDeriv& f        = *(dataF.beginEdit());
    const VecCoord& p  =   dataX.getValue()  ;

    //dmsg_info() << "--addForce" ;
    //    sofa::helper::system::thread::ctime_t start, stop;
    //    sofa::helper::system::thread::CTime timer;
    //
    //    start = timer.getTime();

    int nbTriangles=_topology->getNbTriangles();
    f.resize(p.size());

    for (int i=0; i<nbTriangles; i++)
    {
        accumulateForce(f, p, i);
    }

    dataF.endEdit();

    //    stop = timer.getTime();
    //    dmsg_info() << "---------- time addForce = " << stop-start ;
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularShellForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& datadF, const DataVecDeriv& datadX )
{
    VecDeriv& df        = *(datadF.beginEdit());
    const VecDeriv& dp  =   datadX.getValue()  ;

    double kFactor = mparams->kFactor();

    //dmsg_info() << "--addDForce" ;

    //    sofa::helper::system::thread::ctime_t start, stop;
    //    sofa::helper::system::thread::CTime timer;
    //
    //    start = timer.getTime();

    int nbTriangles=_topology->getNbTriangles();
    df.resize(dp.size());

    for (int i=0; i<nbTriangles; i++)
    {
        applyStiffness(df, dp, i, kFactor);
    }

    datadF.endEdit();

    //    stop = timer.getTime();
    //    dmsg_info() << "time addDForce = " << stop-start ;
}

//#define PRINT

template<class DataTypes>
void TriangularShellForceField<DataTypes>::addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
    StiffnessMatrixFull K_gs;

    // Build Matrix Block for this ForceField
    unsigned int i, j ,n1, n2, row, column, ROW, COLUMN;
    Index node1, node2;

    sofa::core::behavior::MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

    double kFactor = mparams->kFactor();

#ifdef PRINT
    r.matrix->clear();
    dmsg_info() << "Global matrix (" << r.matrix->rowSize() << "x" << r.matrix->colSize() << ")" <<
        " kFactor=" << kFactor ;
    for (unsigned int i=0; i<r.matrix->rowSize(); i++)
    {
        for (unsigned int j=0; j<r.matrix->colSize(); j++)
        {
            dmsg_info() << r.matrix->element(i,j) << ",";
            //dmsg_info() << -r.matrix->element(i,j) / kFactor << ",";
        }
        dmsg_info() ;
    }
#endif
    // XXX: Matrix not necessarily empty!

    for(sofa::Index t=0 ; t != _topology->getNbTriangles() ; ++t)
    {
        const TriangleInformation &tinfo = triangleInf[t];
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
                        //r.matrix->add(ROW, COLUMN, K_gs[row][column]);
                    }
                }
            }
        }
    }

#ifdef PRINT
    dmsg_info() << "Global matrix (" << r.matrix->rowSize() << "x" << r.matrix->colSize() << ")" <<
        " kFactor=" << kFactor ;
    for (unsigned int i=0; i<r.matrix->rowSize(); i++)
    {
        for (unsigned int j=0; j<r.matrix->colSize(); j++)
        {
            dmsg_info() << r.matrix->element(i,j) << ",";
            //dmsg_info() << -r.matrix->element(i,j) / kFactor << ",";
        }
        dmsg_info() ;
    }
#endif

    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// --- Store the initial position of the nodes
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularShellForceField<DataTypes>::initTriangle(const int i, const Index&a, const Index&b, const Index&c, const VecCoord& x0)
{
    type::vector<TriangleInformation>& ti = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &ti[i];

    // Store indices of each vertex
    tinfo->a = a;
    tinfo->b = b;
    tinfo->c = c;

    tinfo->measure.resize(3);
    tinfo->measure[0].id = a; tinfo->measure[0].point = Vec3(0,0,0);
    tinfo->measure[1].id = b; tinfo->measure[1].point = Vec3(1,0,0);
    tinfo->measure[2].id = c; tinfo->measure[2].point = Vec3(0,1,0);

    // Rotation from global to local frame
    Transformation R0;

    // Rest positions expressed in the local frame
    tinfo->restPositions[0] = x0[a].getCenter();
    tinfo->restPositions[1] = x0[b].getCenter();
    tinfo->restPositions[2] = x0[c].getCenter();

    if (d_corotated.getValue()) {
        // Center the element
        Vec3 center = (tinfo->restPositions[0] + tinfo->restPositions[1] +
                       tinfo->restPositions[2])/3;

        tinfo->restPositions[0] -= center;
        tinfo->restPositions[1] -= center;
        tinfo->restPositions[2] -= center;
    }

    computeRotation(R0, tinfo->restPositions);
    tinfo->R = R0;
    tinfo->Rt.transpose(R0);
#ifdef CRQUAT
    tinfo->Q.fromMatrix(tinfo->R);
#endif

    tinfo->restPositions[0] = R0 * tinfo->restPositions[0];
    tinfo->restPositions[1] = R0 * tinfo->restPositions[1];
    tinfo->restPositions[2] = R0 * tinfo->restPositions[2];

// Rest orientations -- inverted (!)
#ifdef CRQUAT
    tinfo->restOrientationsInv[0] = (tinfo->Q * x0[a].getOrientation()).inverse();
    tinfo->restOrientationsInv[1] = (tinfo->Q * x0[b].getOrientation()).inverse();
    tinfo->restOrientationsInv[2] = (tinfo->Q * x0[c].getOrientation()).inverse();
#else
    x0[a].getOrientation().toMatrix(tinfo->restOrientationsInv[0]);
    x0[b].getOrientation().toMatrix(tinfo->restOrientationsInv[1]);
    x0[c].getOrientation().toMatrix(tinfo->restOrientationsInv[2]);

    tinfo->restOrientationsInv[0].transpose( tinfo->R * tinfo->restOrientationsInv[0] );
    tinfo->restOrientationsInv[1].transpose( tinfo->R * tinfo->restOrientationsInv[1] );
    tinfo->restOrientationsInv[2].transpose( tinfo->R * tinfo->restOrientationsInv[2] );
#endif

    // Do some precomputations
    // - directional vectors
    tinfo->d[0] = tinfo->restPositions[0] - tinfo->restPositions[1];
    tinfo->d[1] = tinfo->restPositions[1] - tinfo->restPositions[2];
    tinfo->d[2] = tinfo->restPositions[2] - tinfo->restPositions[0];

    // - squared lengths
    tinfo->l2[0] = tinfo->d[0][0]*tinfo->d[0][0] + tinfo->d[0][1]*tinfo->d[0][1];
    tinfo->l2[1] = tinfo->d[1][0]*tinfo->d[1][0] + tinfo->d[1][1]*tinfo->d[1][1];
    tinfo->l2[2] = tinfo->d[2][0]*tinfo->d[2][0] + tinfo->d[2][1]*tinfo->d[2][1];

    // - triangle area
    tinfo->area = helper::rabs(tinfo->d[2][0]*(-tinfo->d[0][1]) - (-tinfo->d[0][0])*tinfo->d[2][1])/2;

    // Compute stiffness matrix for membrane element
    computeStiffnessMatrixMembrane(tinfo->stiffnessMatrixMembrane, *tinfo);
    //dmsg_info() << "Km^e=" << tinfo->stiffnessMatrixMembrane ;

    // Compute stiffness matrix for bending plate elemnt
    computeStiffnessMatrixBending(tinfo->stiffnessMatrixBending, *tinfo);
    //dmsg_info() << "Kb^e=" << tinfo->stiffnessMatrixBending ;


    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// Computes the rotation from global frame to local triangle frame
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularShellForceField<DataTypes>::computeRotation(Transformation& R, const VecCoord &x, const Index &a, const Index &b, const Index &c)
{

    if (!d_corotated.getValue()) {
        // Return identity matrix
        R = Transformation(Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,1));
        return;
    }
    
    // First vector on first edge
    // Second vector in the plane of the two first edges
    // Third vector orthogonal to first and second

    Vec3 edgex = x[b].getCenter() - x[a].getCenter();

    Vec3 edgey = x[c].getCenter() - x[a].getCenter();

    Vec3 edgez = cross(edgex, edgey);

    edgey = cross(edgez, edgex);

    edgex.normalize();
    edgey.normalize();
    edgez.normalize();

    R = Transformation(edgex, edgey, edgez);
    //Qframe.fromMatrix(Transformation(edgex, edgey, edgez));
}

template <class DataTypes>
void TriangularShellForceField<DataTypes>::computeRotation(Transformation& R, const type::fixed_array<Vec3, 3> &x)
{
    if (!d_corotated.getValue()) {
        // Return identity matrix
        R = Transformation(Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,1));
        return;
    }
    
    // First vector on first edge
    // Second vector in the plane of the two first edges
    // Third vector orthogonal to first and second

    Vec3 edge_x = x[1] - x[0];
    Vec3 edge_y = x[2] - x[0];

    Vec3 edge_z = cross(edge_x, edge_y);

    edge_y = cross(edge_z, edge_x);

    edge_x.normalize();
    edge_y.normalize();
    edge_z.normalize();

    R = Transformation(edge_x, edge_y, edge_z);
    //Qframe.fromMatrix(Transformation(edge_x, edge_y, edge_z));
}


// --------------------------------------------------------------------------------------
// ---  Compute material stiffness
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularShellForceField<DataTypes>::computeMaterialStiffness()
{
    Real E = d_young.getValue(),
        nu = d_poisson.getValue(),
        t = d_thickness.getValue();

    materialMatrix[0][0] = 1.0;
    materialMatrix[0][1] = nu;
    materialMatrix[0][2] = 0;
    materialMatrix[1][0] = nu;
    materialMatrix[1][1] = 1.0;
    materialMatrix[1][2] = 0;
    materialMatrix[2][0] = 0;
    materialMatrix[2][1] = 0;
    materialMatrix[2][2] = (1 - nu)/2.0;

    materialMatrix *= E / (1.0 - nu * nu);

    materialMatrixMembrane = materialMatrix;
    materialMatrixBending = materialMatrix;

    // Integrate through the shell thickness
    materialMatrixMembrane *= t;
    if (d_isShellveryThin.getValue())
        materialMatrixBending *= t*t;
    else
        materialMatrixBending *= t*t*t / 12;
}


// -------------------------------------------------------------------------------------------------------------
// --- Compute displacement vector D as the difference between current position and initial position
// -------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularShellForceField<DataTypes>::computeDisplacement(Displacement &Dm, Displacement &Db, const VecCoord &x, const Index elementIndex)
{
    type::vector<TriangleInformation>& ti = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &ti[elementIndex];

    Index a = tinfo->a;
    Index b = tinfo->b;
    Index c = tinfo->c;

    // Compute local (in-plane) positions
    tinfo->deformedPositions[0] = x[a].getCenter();
    tinfo->deformedPositions[1] = x[b].getCenter();
    tinfo->deformedPositions[2] = x[c].getCenter();

    if (d_corotated.getValue()) {
        Vec3 center = (tinfo->deformedPositions[0] + tinfo->deformedPositions[1] +
                       tinfo->deformedPositions[2])/3;

        tinfo->deformedPositions[0] -= center;
        tinfo->deformedPositions[1] -= center;
        tinfo->deformedPositions[2] -= center;
    }

    // Compute rotation to local (in-plane) frame
    computeRotation(tinfo->R, tinfo->deformedPositions);
    tinfo->Rt.transpose(tinfo->R);
#ifdef CRQUAT
    tinfo->Q.fromMatrix(tinfo->R);
#endif

    // Compute local (in-plane) positions
    tinfo->deformedPositions[0] = tinfo->R * tinfo->deformedPositions[0];
    tinfo->deformedPositions[1] = tinfo->R * tinfo->deformedPositions[1];
    tinfo->deformedPositions[2] = tinfo->R * tinfo->deformedPositions[2];

    // Displacements
    Vec3 uA = tinfo->deformedPositions[0] - tinfo->restPositions[0];
    Vec3 uB = tinfo->deformedPositions[1] - tinfo->restPositions[1];
    Vec3 uC = tinfo->deformedPositions[2] - tinfo->restPositions[2];

// Rotations
#ifdef CRQUAT
    Quat qA = (tinfo->Q * x[a].getOrientation()) * tinfo->restOrientationsInv[0];
    Quat qB = (tinfo->Q * x[b].getOrientation()) * tinfo->restOrientationsInv[1];
    Quat qC = (tinfo->Q * x[c].getOrientation()) * tinfo->restOrientationsInv[2];
#else
    Transformation tmpA, tmpB, tmpC;
    x[a].getOrientation().toMatrix(tmpA);
    x[b].getOrientation().toMatrix(tmpB);
    x[c].getOrientation().toMatrix(tmpC);

    Quat qA; qA.fromMatrix( tinfo->R * tmpA * tinfo->restOrientationsInv[0] );
    Quat qB; qB.fromMatrix( tinfo->R * tmpB * tinfo->restOrientationsInv[1] );
    Quat qC; qC.fromMatrix( tinfo->R * tmpC * tinfo->restOrientationsInv[2] );
    // TODO: can we do this without the quaternions?
#endif
    Vec3 rA = qA.toEulerVector();
    Vec3 rB = qB.toEulerVector();
    Vec3 rC = qC.toEulerVector();
    //dmsg_info() << "Θ: " << rA << " | " << rB << " | " << rC ;

    // Membrane
    Dm[0] = uA[0];
    Dm[1] = uA[1];
    Dm[2] = rA[2];

    Dm[3] = uB[0];
    Dm[4] = uB[1];
    Dm[5] = rB[2];

    Dm[6] = uC[0];
    Dm[7] = uC[1];
    Dm[8] = rC[2];

    // Bending plate
    Db[0] = uA[2];
    Db[1] = rA[0];
    Db[2] = rA[1];

    Db[3] = uB[2];
    Db[4] = rB[0];
    Db[5] = rB[1];

    Db[6] = uC[2];
    Db[7] = rC[0];
    Db[8] = rC[1];

    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularShellForceField<DataTypes>::accumulateForce(VecDeriv &f, const VecCoord &x, const Index elementIndex)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Get the indices of the 3 vertices for the current triangle
    const Index& a = tinfo->a;
    const Index& b = tinfo->b;
    const Index& c = tinfo->c;

    // Compute in-plane displacements
    Displacement Dm, Db;
    computeDisplacement(Dm, Db, x, elementIndex);

    // Compute the membrane and bending plate forces on this element
    Displacement Fm, Fb;
    computeForce(Fm, Dm, Fb, Db, elementIndex);

    if (this->f_printLog.getValue()) {
        dmsg_info() << "E: " << elementIndex << "\tu: " << Dm << "\tf: " << Fm << "\n";
        dmsg_info() << "E: " << elementIndex << "\tuB: " << Db << "\tfB: " << Fb << "\n";
        dmsg_info() << "   xg [ " << a << "/" << b << "/" << c << " - "
                    << x[a] << ", " << x[b] << ", " << x[c] << "\n";
        dmsg_info() << "   xl [ " << tinfo->deformedPositions[0] << ", " << tinfo->deformedPositions[1] << ", " << tinfo->deformedPositions[2] << "\n";
        dmsg_info() << "   fg: " <<
            tinfo->Rt * Vec3(Fm[0], Fm[1], Fb[0]) << " " << tinfo->R * Vec3(Fb[1], Fb[2], Fm[2]) << " | " <<
            tinfo->Rt * Vec3(Fm[3], Fm[4], Fb[3]) << " " << tinfo->R * Vec3(Fb[4], Fb[5], Fm[5]) << " | " <<
            tinfo->Rt * Vec3(Fm[6], Fm[7], Fb[6]) << " " << tinfo->R * Vec3(Fb[7], Fb[8], Fm[8]) ;

    }

    m_minStress = std::numeric_limits<Real>::max();
    m_maxStress = std::numeric_limits<Real>::lowest();
    // Compute the measure (stress/strain)
    if (bMeasureStrain) {
        type::vector<Real> &values = *d_measuredValues.beginEdit();
        for (unsigned int i=0; i< tinfo->measure.size(); i++) {
            Vec3 strain = tinfo->measure[i].B * Dm + tinfo->measure[i].Bb * Db;
            // Norm from strain in x and y
            // NOTE: Shear strain is not included
            values[ tinfo->measure[i].id ] = helper::rsqrt(
                strain[0] * strain[0] + strain[1] * strain[1]);
            if (values[ tinfo->measure[i].id ] < m_minStress)
                m_minStress = values[ tinfo->measure[i].id ];
            if (values[ tinfo->measure[i].id ] > m_maxStress)
                m_maxStress = values[ tinfo->measure[i].id ];
        }
        d_measuredValues.endEdit();
    } else if (bMeasureStress) {
        type::vector<Real> &values = *d_measuredValues.beginEdit();
        for (unsigned int i=0; i< tinfo->measure.size(); i++) {
            Vec3 stress = materialMatrix * tinfo->measure[i].B * Dm
                          + materialMatrix * tinfo->measure[i].Bb * Db;
            // Von Mises stress criterion (plane stress)
            values[ tinfo->measure[i].id ] = helper::rsqrt(
                stress[0] * stress[0] - stress[0] * stress[1]
                + stress[1] * stress[1] + 3 * stress[2] * stress[2]);
            if (values[ tinfo->measure[i].id ] < m_minStress)
                m_minStress = values[ tinfo->measure[i].id ];
            if (values[ tinfo->measure[i].id ] > m_maxStress)
                m_maxStress = values[ tinfo->measure[i].id ];
        }
        d_measuredValues.endEdit();
    }

    // Transform forces back into global frame
    getVCenter(f[a]) -= tinfo->Rt * Vec3(Fm[0], Fm[1], Fb[0]);
    getVCenter(f[b]) -= tinfo->Rt * Vec3(Fm[3], Fm[4], Fb[3]);
    getVCenter(f[c]) -= tinfo->Rt * Vec3(Fm[6], Fm[7], Fb[6]);

    getVOrientation(f[a]) -= tinfo->Rt * Vec3(Fb[1], Fb[2], Fm[2]);
    getVOrientation(f[b]) -= tinfo->Rt * Vec3(Fb[4], Fb[5], Fm[5]);
    getVOrientation(f[c]) -= tinfo->Rt * Vec3(Fb[7], Fb[8], Fm[8]);

    triangleInfo.endEdit();
}


// ----------------------------------------------------------------------------------------------------------------------
// --- Compute the stiffness matrix for membrane element
// ----------------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularShellForceField<DataTypes>::computeStiffnessMatrixMembrane(StiffnessMatrix &K, TriangleInformation &tinfo)
{
    if (csMembrane == NULL)
        return;

    (this->*csMembrane)(K, tinfo);

    if (this->f_printLog.getValue())
        dmsg_info() << "Km = " << K ;
}


// ----------------------------------------------------------------------------------------------------------------------
// --- Compute the stiffness matrix for bending plate element
// ----------------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularShellForceField<DataTypes>::computeStiffnessMatrixBending(StiffnessMatrix &K, TriangleInformation &tinfo)
{
    if (csBending == NULL)
        return;

    (this->*csBending)(K, tinfo);

    if (this->f_printLog.getValue())
        dmsg_info() << "Kb = " << K ;
}

// --------------------------------------------------------------------------------------
// ---  Compute force F = K * u
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularShellForceField<DataTypes>::computeForce(Displacement &Fm, const Displacement& Dm, Displacement &Fb, const Displacement& Db,const Index elementIndex)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation &tinfo = triangleInf[elementIndex];

    // Compute forces
    Fm = tinfo.stiffnessMatrixMembrane * Dm;
    Fb = tinfo.stiffnessMatrixBending * Db;

    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularShellForceField<DataTypes>::applyStiffness(VecDeriv& v, const VecDeriv& dx, const Index elementIndex, const double kFactor)
{
    type::vector<TriangleInformation>& ti = *(triangleInfo.beginEdit());
    TriangleInformation &tinfo = ti[elementIndex];

    // Get the indices of the 3 vertices for the current triangle
    const Index& a = tinfo.a;
    const Index& b = tinfo.b;
    const Index& c = tinfo.c;

    // Computes displacements
    Displacement Dm, Db;
    Vec3 x_a, x_b, x_c;
    Vec3 r_a, r_b, r_c;

    x_a = tinfo.R * getVCenter(dx[a]);
    r_a = tinfo.R * getVOrientation(dx[a]);

    x_b = tinfo.R * getVCenter(dx[b]);
    r_b = tinfo.R * getVOrientation(dx[b]);

    x_c = tinfo.R * getVCenter(dx[c]);
    r_c = tinfo.R * getVOrientation(dx[c]);

    Dm[0] = x_a[0];
    Dm[1] = x_a[1];
    Dm[2] = r_a[2];

    Dm[3] = x_b[0];
    Dm[4] = x_b[1];
    Dm[5] = r_b[2];

    Dm[6] = x_c[0];
    Dm[7] = x_c[1];
    Dm[8] = r_c[2];

    Db[0] = x_a[2];
    Db[1] = r_a[0];
    Db[2] = r_a[1];

    Db[3] = x_b[2];
    Db[4] = r_b[0];
    Db[5] = r_b[1];

    Db[6] = x_c[2];
    Db[7] = r_c[0];
    Db[8] = r_c[1];

    // Compute dF
    Displacement dFm, dFb;
    dFm = tinfo.stiffnessMatrixMembrane * Dm;
    dFb = tinfo.stiffnessMatrixBending * Db;

    if (this->f_printLog.getValue()) {
        dmsg_info() << "E: " << elementIndex << "\tdu: " << Dm << "\tdf: " << dFm << "\n";
        dmsg_info() << "E: " << elementIndex << "\tduB: " << Db << "\tdfB: " << dFb << "\n";
        dmsg_info() << "   dfg: " <<
            tinfo.Rt * Vec3(dFm[0], dFm[1], dFb[0]) * kFactor << " " << tinfo.Rt * Vec3(dFb[1], dFb[2], dFm[2]) * kFactor << " | " <<
            tinfo.Rt * Vec3(dFm[3], dFm[4], dFb[3]) * kFactor << " " << tinfo.Rt * Vec3(dFb[4], dFb[5], dFm[5]) * kFactor << " | " <<
            tinfo.Rt * Vec3(dFm[6], dFm[7], dFb[6]) * kFactor << " " << tinfo.Rt * Vec3(dFb[7], dFb[8], dFm[8]) * kFactor ;
    }

    // Transform into global frame
    getVCenter(v[a]) -= tinfo.Rt * Vec3(dFm[0], dFm[1], dFb[0]) * kFactor;
    getVCenter(v[b]) -= tinfo.Rt * Vec3(dFm[3], dFm[4], dFb[3]) * kFactor;
    getVCenter(v[c]) -= tinfo.Rt * Vec3(dFm[6], dFm[7], dFb[6]) * kFactor;

    getVOrientation(v[a]) -= tinfo.Rt * Vec3(dFb[1], dFb[2], dFm[2]) * kFactor;
    getVOrientation(v[b]) -= tinfo.Rt * Vec3(dFb[4], dFb[5], dFm[5]) * kFactor;
    getVOrientation(v[c]) -= tinfo.Rt * Vec3(dFb[7], dFb[8], dFm[8]) * kFactor;

    triangleInfo.endEdit();
}


template<class DataTypes>
void TriangularShellForceField<DataTypes>::convertStiffnessMatrixToGlobalSpace(StiffnessMatrixFull &K_gs, const TriangleInformation &tinfo)
{
    // Firstly, add all degrees of freedom (we add the unused translation in z)
    StiffnessMatrixFull K_18x18;
    K_18x18.clear();
    unsigned int ig, jg;


    // Copy the stiffness matrix into 18x18 matrix (the new index of each block in global matrix is a combination of 0, 6 and 12 in indices)
    const StiffnessMatrix &K = tinfo.stiffnessMatrixMembrane;
    for (unsigned int bx=0; bx<3; bx++)
    {
        // Global row index
        ig = 6*bx;

        for (unsigned int by=0; by<3; by++)
        {
            // Global column index
            jg = 6*by;

            // linear X
            K_18x18[ig+0][jg+0] = K[3*bx+0][3*by+0]; // linear X
            K_18x18[ig+0][jg+1] = K[3*bx+0][3*by+1]; // linear Y
            K_18x18[ig+0][jg+5] = K[3*bx+0][3*by+2]; // angular Z

            // linear Y
            K_18x18[ig+1][jg+0] = K[3*bx+1][3*by+0]; // linear X
            K_18x18[ig+1][jg+1] = K[3*bx+1][3*by+1]; // linear Y
            K_18x18[ig+1][jg+5] = K[3*bx+1][3*by+2]; // angular Z

            // angular Z
            K_18x18[ig+5][jg+0] = K[3*bx+2][3*by+0]; // linear X
            K_18x18[ig+5][jg+1] = K[3*bx+2][3*by+1]; // linear Y
            K_18x18[ig+5][jg+5] = K[3*bx+2][3*by+2]; // angular Z
        }
    }


    // Copy the stiffness matrix by block 3x3 into global matrix (the new index of each bloc into global matrix is a combination of 2, 8 and 15 in indices)
    const StiffnessMatrix &K_bending = tinfo.stiffnessMatrixBending;
    for (unsigned int bx=0; bx<3; bx++)
    {
        // Global row index
        ig = 6*bx+2;

        for (unsigned int by=0; by<3; by++)
        {
            // Global column index
            jg = 6*by+2;

            // Iterates over the indices of the 3x3 block
            for (unsigned int i=0; i<3; i++)
            {
                for (unsigned int j=0; j<3; j++)
                {
                    K_18x18[ig+i][jg+j] += K_bending[3*bx+i][3*by+j];
                }
            }

        }
    }


    // Extend rotation matrix and its transpose
    StiffnessMatrixFull R18x18, Rt18x18;

    for(unsigned int i=0;i<3;++i)
    {
        for(unsigned int j=0;j<3;++j)
        {
            R18x18[i][j] = R18x18[i+3][j+3] = R18x18[i+6][j+6] = R18x18[i+9][j+9] = R18x18[i+12][j+12] = R18x18[i+15][j+15] = tinfo.R[i][j];
            Rt18x18[i][j] = Rt18x18[i+3][j+3] = Rt18x18[i+6][j+6] = Rt18x18[i+9][j+9] = Rt18x18[i+12][j+12] = Rt18x18[i+15][j+15] = tinfo.Rt[i][j];
        }
    }

    // Then we put the stifness matrix into the global frame
    K_gs = Rt18x18 * K_18x18 * R18x18;
    //dmsg_info() << "R=" << R18x18 << " -- " << Rt18x18 ;
    //dmsg_info() << "K_gs=" << K_gs ;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Membrane elements

template <class DataTypes>
void TriangularShellForceField<DataTypes>::computeStiffnessMatrixCST(StiffnessMatrix &K, TriangleInformation &tinfo)
{
    // NOTE: we could use the ANDES template below for CST too, but this is a little faster
    Mat<9,3, Real> L;
    Mat<3,9, Real> Lt;
    L.clear();

    L[0][0] = tinfo.d[1][1];
    L[0][1] = 0;
    L[0][2] = -tinfo.d[1][0];

    L[1][0] = 0;
    L[1][1] = -tinfo.d[1][0];
    L[1][2] = tinfo.d[1][1];

    L[3][0] = tinfo.d[2][1];
    L[3][1] = 0;
    L[3][2] = -tinfo.d[2][0];

    L[4][0] = 0;
    L[4][1] = -tinfo.d[2][0];
    L[4][2] = tinfo.d[2][1];

    L[6][0] = tinfo.d[0][1];
    L[6][1] = 0;
    L[6][2] = -tinfo.d[0][0];

    L[7][0] = 0;
    L[7][1] = -tinfo.d[0][0];
    L[7][2] = tinfo.d[0][1];

    // Now should be:
    //  L *= h/2;
    //  K = 1/(h*Area) * L * Em * L^T;
    // But we may simplify this a little and we also have 'h' already in Em

    Lt.transpose(L);

    K = L * materialMatrixMembrane * Lt / (4*tinfo.area);

    // Compute strain-displacement matrix
    for (unsigned int i=0; i< tinfo.measure.size(); i++) {
        tinfo.measure[i].B = Lt / (2*tinfo.area);
    }

}

// The ANDES template for membrane element
// See for example: C. A. Felippa, A study of optimal membrane triangles with drilling freedoms, 2003
template <class DataTypes>
void TriangularShellForceField<DataTypes>::andesTemplate(StiffnessMatrix &K, const TriangleInformation &tinfo,
                                                         const Real alpha, const AndesBeta &beta)
{
    Real h = d_thickness.getValue();
    Real A4 = 4*tinfo.area;

    // Force-lumping matrix
    Mat<9,3, Real> L;
    Mat<3,9, Real> Lt;

    // y23,    0,      x32,
    L[0][0] = tinfo.d[1][1];
    L[0][1] = 0;
    L[0][2] = -tinfo.d[1][0];

    // 0,      x32,    y23,
    L[1][0] = 0;
    L[1][1] = -tinfo.d[1][0];
    L[1][2] = tinfo.d[1][1];

    // alpha/6*y23*(y13 - y21), alpha/6*x32*(x31 - x12), alpha/3*(x31*y13 - x12*y21),
    L[2][0] = alpha/6.0*tinfo.d[1][1]*(tinfo.d[0][1] - tinfo.d[2][1]);
    L[2][1] = alpha/6.0*(-tinfo.d[1][0])*(tinfo.d[2][0] - tinfo.d[0][0]);
    L[2][2] = alpha/3.0*(tinfo.d[0][0]*tinfo.d[0][1] - tinfo.d[2][0]*tinfo.d[2][1]),

        // y31,    0,      x13,
        L[3][0] = tinfo.d[2][1];
    L[3][1] = 0;
    L[3][2] = -tinfo.d[2][0];

    // 0,      x13,    y31,
    L[4][0] = 0;
    L[4][1] = -tinfo.d[2][0];
    L[4][2] = tinfo.d[2][1];

    // alpha/6*y31*(y21 - y32), alpha/6*x13*(x12 - x23), alpha/3*(x12*y21 - x23*y32),
    L[5][0] = alpha/6*tinfo.d[2][1]*(tinfo.d[1][1] - tinfo.d[0][1]);
    L[5][1] = alpha/6*(-tinfo.d[2][0])*(tinfo.d[0][0] - tinfo.d[1][0]);
    L[5][2] = alpha/3*(tinfo.d[1][0]*tinfo.d[1][1] - tinfo.d[0][0]*tinfo.d[0][1]);

    // y12,    0,      x21,
    L[6][0] = tinfo.d[0][1];
    L[6][1] = 0;
    L[6][2] = -tinfo.d[0][0];

    // 0,      x21,    y12,
    L[7][0] = 0;
    L[7][1] = -tinfo.d[0][0];
    L[7][2] = tinfo.d[0][1];

    // alpha/6*y12*(y32 - y13), alpha/6*x21*(x23 - x31), alpha/3*(x23*y32 - x31*y13)
    L[8][0] = alpha/6*tinfo.d[0][1]*(tinfo.d[2][1] - tinfo.d[1][1]);
    L[8][1] = alpha/6*(-tinfo.d[0][0])*(tinfo.d[1][0] - tinfo.d[2][0]);
    L[8][2] = alpha/3*(tinfo.d[2][0]*tinfo.d[2][1] - tinfo.d[1][0]*tinfo.d[1][1]);

    // Now should be:
    //  L *= h/2;
    //  Kb = 1/(h*Area) * L * Em * L^T;
    // But we may simplify this a little and we also have 'h' already in Em

    Lt.transpose(L);

    // Add base stiffness matrix
    K = L * materialMatrixMembrane * Lt / A4;

    // Transformation matrix from DOFs to hierarchical drilling freedoms
    Mat<3,9, Real> T(
        Vec9(tinfo.d[1][0], tinfo.d[1][1], -A4, tinfo.d[2][0], tinfo.d[2][1],   0, tinfo.d[0][0], tinfo.d[0][1],   0),
        Vec9(tinfo.d[1][0], tinfo.d[1][1],   0, tinfo.d[2][0], tinfo.d[2][1], -A4, tinfo.d[0][0], tinfo.d[0][1],   0),
        Vec9(tinfo.d[1][0], tinfo.d[1][1],   0, tinfo.d[2][0], tinfo.d[2][1],   0, tinfo.d[0][0], tinfo.d[0][1], -A4)
        );
    T /= -A4;

    Mat<9,3, Real> Tt;
    Tt.transpose(T);

    Mat<3,3, Real> Q1(
        Vec3(beta[1], beta[2], beta[3])/tinfo.l2[0],
        Vec3(beta[4], beta[5], beta[6])/tinfo.l2[1],
        Vec3(beta[7], beta[8], beta[9])/tinfo.l2[2]);
    Q1 *= 2.0*tinfo.area/3.0;

    Mat<3,3, Real> Q2(
        Vec3(beta[9], beta[7], beta[8])/tinfo.l2[0],
        Vec3(beta[3], beta[1], beta[2])/tinfo.l2[1],
        Vec3(beta[6], beta[4], beta[5])/tinfo.l2[2]);
    Q2 *= 2.0*tinfo.area/3.0;

    Mat<3,3, Real> Q3(
        Vec3(beta[5], beta[6], beta[4])/tinfo.l2[0],
        Vec3(beta[8], beta[9], beta[7])/tinfo.l2[1],
        Vec3(beta[2], beta[3], beta[1])/tinfo.l2[2]);
    Q3 *= 2.0*tinfo.area/3.0;

    Mat<3,3, Real> Q4 = (Q1 + Q2)/2.0, Q4t;
    Mat<3,3, Real> Q5 = (Q2 + Q3)/2.0, Q5t;
    Mat<3,3, Real> Q6 = (Q3 + Q1)/2.0, Q6t;
    Q4t.transpose(Q4);
    Q5t.transpose(Q5);
    Q6t.transpose(Q6);

    Mat<3,3, Real> Te, Tet, Enat;
    Te[0][0] = -tinfo.d[1][1]*tinfo.d[2][1]*tinfo.l2[0];
    Te[0][1] = -tinfo.d[2][1]*tinfo.d[0][1]*tinfo.l2[1];
    Te[0][2] = -tinfo.d[0][1]*tinfo.d[1][1]*tinfo.l2[2];

    Te[1][0] = -tinfo.d[1][0]*tinfo.d[2][0]*tinfo.l2[0];
    Te[1][1] = -tinfo.d[2][0]*tinfo.d[0][0]*tinfo.l2[1];
    Te[1][2] = -tinfo.d[1][0]*tinfo.d[1][0]*tinfo.l2[2];

    Te[2][0] = (tinfo.d[1][1]*tinfo.d[2][0] + tinfo.d[1][0]*tinfo.d[2][1])*tinfo.l2[0];
    Te[2][1] = (tinfo.d[2][1]*tinfo.d[1][0] + tinfo.d[2][0]*tinfo.d[0][1])*tinfo.l2[1];
    Te[2][2] = (tinfo.d[0][1]*tinfo.d[1][0] + tinfo.d[0][0]*tinfo.d[1][1])*tinfo.l2[2];

    Te /= tinfo.area*A4;

    Tet.transpose(Te);
    Enat = Tet * materialMatrixMembrane * Te;

    // Add higher order stiffness matrix
    K += beta[0] * Real(3.0/4.0) * h * tinfo.area * Tt *
         // Higher order stiffness in terms of hierarchical rotations
         (Q4t * Enat * Q4 + Q5t * Enat * Q5 + Q6t * Enat * Q6)
         * T;
}

template <class DataTypes>
void TriangularShellForceField<DataTypes>::computeStiffnessMatrixAll3I(StiffnessMatrix &K, TriangleInformation &tinfo)
{
    return andesTemplate(K, tinfo, 1.0, AndesBeta(4.0/9.0, 1.0/12.0, 5.0/12.0, 1.0/2.0, 0.0, 1.0/3.0, -1.0/3.0, -1.0/12.0, -1.0/2.0, -5.0/12.0));
}

template <class DataTypes>
void TriangularShellForceField<DataTypes>::computeStiffnessMatrixAll3M(StiffnessMatrix &K, TriangleInformation &tinfo)
{
    return andesTemplate(K, tinfo, 1.0, AndesBeta(4.0/9.0, 1.0/4.0, 5.0/4.0, 3.0/2.0, 0.0, 1.0, -1.0, -1.0/4.0, -3.0/2.0, -5.0/4.0));
}

template <class DataTypes>
void TriangularShellForceField<DataTypes>::computeStiffnessMatrixAllLS(StiffnessMatrix &K, TriangleInformation &tinfo)
{
    return andesTemplate(K, tinfo, 1.0, AndesBeta(4.0/9.0, 3.0/20.0, 3.0/4.0, 9.0/10.0, 0.0, 3.0/5.0, -3.0/5.0, -3.0/20.0, -9.0/10.0, -3.0/4.0));
}

template <class DataTypes>
void TriangularShellForceField<DataTypes>::computeStiffnessMatrixLSTRet(StiffnessMatrix &K, TriangleInformation &tinfo)
{
    return andesTemplate(K, tinfo, 4.0/3.0, AndesBeta(1.0/2.0, 2.0/3.0, -2.0/3.0, 0.0, 0.0, -4.0/3.0, 4.0/3.0, -2.0/3.0, 0.0, 2.0/3.0));
}

// Optimal ANDES membrane element
// See: C. A. Felippa, A study of optimal membrane triangles with drilling freedoms, 2003
template <class DataTypes>
void TriangularShellForceField<DataTypes>::computeStiffnessMatrixAndesOpt(StiffnessMatrix &K, TriangleInformation &tinfo)
{
    Real beta0 = helper::rmax(0.5 - 2.0*d_poisson.getValue()*d_poisson.getValue(), 0.01);
    return andesTemplate(K, tinfo, 3.0/2.0, AndesBeta(beta0, 1.0, 2.0, 1.0, 0.0, 1.0, -1.0, -1.0, -1.0, -2.0));
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//  Bending plate elements

template <class DataTypes>
void TriangularShellForceField<DataTypes>::computeStiffnessMatrixDKT(StiffnessMatrix &K, TriangleInformation &tinfo)
{
    // Weights and abscissa for 6-point Gaussian quadrature of a triangle
    // (source: http://www.electromagnetics.biz/integration.htm)
    Vec<6, Real> gx(0.8168476,  0.09157621, 0.09157621, 0.1081030, 0.4459485, 0.4459485);
    Vec<6, Real> gy(0.09157621, 0.8168476,  0.09157621, 0.4459485, 0.1081030, 0.4459485);
    Vec<6, Real> gw(0.05497587, 0.05497587, 0.05497587, 0.1116908, 0.1116908, 0.1116908);

    // Integrage over triangle area
    K.clear();
    //dmsg_info() << "B=\n";
    for (int i=0; i<6; i++) {
        StrainDisplacement B;
        Mat<9,3, Real> Bt;
        dktSD(B, tinfo, gx[i], gy[i]);  // Compute strain-displacement matrix
        //dmsg_info() << "  " << B ;
        Bt.transpose(B);
        K += gw[i] * Bt * materialMatrixBending * B;
    }
    // Compute strain-displacement matrix
    for (unsigned int i=0; i< tinfo.measure.size(); i++) {
        dktSD(tinfo.measure[i].Bb, tinfo,
              tinfo.measure[i].point[0],
              tinfo.measure[i].point[1]);
    }

    //dmsg_info() ;
    K *= 2*tinfo.area;
}

template <class DataTypes>
void TriangularShellForceField<DataTypes>::dktSD(StrainDisplacement &B, const TriangleInformation &tinfo, const Real xi, const Real eta)
{
    Real P4 = -6*tinfo.d[1][0]/tinfo.l2[1];
    Real P5 = -6*tinfo.d[2][0]/tinfo.l2[2];
    Real P6 = -6*tinfo.d[0][0]/tinfo.l2[0];

    Real t4 = -6*tinfo.d[1][1]/tinfo.l2[1];
    Real t5 = -6*tinfo.d[2][1]/tinfo.l2[2];
    Real t6 = -6*tinfo.d[0][1]/tinfo.l2[0];

    Real q4 = 3*tinfo.d[1][0]*tinfo.d[1][1]/tinfo.l2[1];
    Real q5 = 3*tinfo.d[2][0]*tinfo.d[2][1]/tinfo.l2[2];
    Real q6 = 3*tinfo.d[0][0]*tinfo.d[0][1]/tinfo.l2[0];

    Real r4 = 3*tinfo.d[1][1]*tinfo.d[1][1]/tinfo.l2[1];
    Real r5 = 3*tinfo.d[2][1]*tinfo.d[2][1]/tinfo.l2[2];
    Real r6 = 3*tinfo.d[0][1]*tinfo.d[0][1]/tinfo.l2[0];


    Real tmp = 1.0 - 2.0*xi;
    Vec9 Hx_dxi(
        P6*tmp + eta * (P5 - P6),
        q6*tmp - eta * (q5 + q6),
        -4.0 + 6.0*(xi + eta) + r6*tmp - eta*(r5 + r6),
        -P6*tmp + eta*(P4 + P6),
        q6*tmp - eta*(q6 - q4),
        -2.0 + 6.0*xi + r6*tmp + eta*(r4 - r6),
        - eta*(P5 + P4),
        eta*(q4 - q5),
        - eta*(r5 - r4) 
        );

    //tmp = 1.0 - 2.0*xi;
    Vec9 Hy_dxi(
        t6*tmp + eta*(t5 - t6),
        1.0 + r6*tmp - eta*(r5 + r6),
        -q6*tmp + eta*(q5 + q6),
        -t6*tmp + eta*(t4 + t6),
        -1.0 + r6*tmp + eta*(r4 - r6),
        -q6*tmp - eta*(q4 - q6),
        -eta*(t4 + t5),
        eta*(r4 - r5),
        -eta*(q4 - q5)
        );

    tmp = 1.0 - 2.0*eta;
    Vec9 Hx_deta(
        -P5*tmp - xi*(P6 - P5),
        q5*tmp - xi*(q5 + q6),
        -4.0 + 6.0*(xi + eta) + r5*tmp - xi*(r5 + r6),
        xi*(P4 + P6),
        xi*(q4 - q6),
        -xi*(r6 - r4),
        P5*tmp - xi*(P4 + P5),
        q5*tmp + xi*(q4 - q5),
        -2.0 + 6.0*eta + r5*tmp + xi*(r4 - r5)
        );

    // tmp = 1.0 - 2.0*eta;
    Vec9 Hy_deta(
        -t5 * tmp - xi*(t6 - t5),
        1 + r5*tmp - xi*(r5 + r6),
        -q5 * tmp + xi*(q5 + q6),
        xi * (t4 + t6),
        xi * (r4 - r6),
        -xi * (q4 - q6),
        t5 * tmp - xi*(t4 + t5),
        -1 + r5*tmp + xi*(r4 - r5),
        -q5*tmp - xi*(q4 - q5)
        );

    B = StrainDisplacement(
        Hx_dxi*tinfo.d[2][1] + Hx_deta*tinfo.d[0][1],
        -Hy_dxi*tinfo.d[2][0] - Hy_deta*tinfo.d[0][0],
        -Hx_dxi*tinfo.d[2][0] - Hx_deta*tinfo.d[0][0] + Hy_dxi*tinfo.d[2][1] + Hy_deta*tinfo.d[0][1]
        );

    B /= 2*tinfo.area;
}

template <class DataTypes>
void TriangularShellForceField<DataTypes>::draw(const core::visual::VisualParams* vparams){
    SCOPED_TIMER("draw_TriangularShellForceField");

    if (!vparams->displayFlags().getShowForceFields()) {
        return;
    }

    // a flag indicating if there are infomation to draw
    bool p_computeDrawInfo;
    p_computeDrawInfo = d_showArrow.getValue() ||
                        d_showTriangle.getValue() ||
                        d_showMeasuredValue.getValue();

    if (!p_computeDrawInfo) {
        return;
    }
    

    const auto stateLifeCycle = vparams->drawTool()->makeStateLifeCycle();

    if (vparams->displayFlags().getShowWireFrame())
        vparams->drawTool()->setPolygonMode(0, true);

    vparams->drawTool()->disableLighting();

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
    const type::vector<TriangleInformation>& triangleInf = triangleInfo.getValue();
    const auto& triangles = _topology->getTriangles();
    const Size nbTriangles = triangles.size();

   // Draw arrows from using shell triangle info (tinfo.R) in the center of the triangle
    if (d_showArrow.getValue())
    {
        const VecCoord& positions =this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();
        const Real radius = d_arrowRadius.getValue();
        for (unsigned int i=0; i< nbTriangles; i++)
        {
          const TriangleInformation &tinfo = triangleInf[i];
          Vec3 vx = tinfo.R * Vec3(1, 0, 0);
          Vec3 vy = tinfo.R * Vec3(0, 1, 0);
          Vec3 vz = tinfo.R * Vec3(0, 0, 1);
     
          auto triangle = _topology->getTriangle(i);
          // get triangle i indices
          auto a = triangle[0];
          auto b = triangle[1];
          auto c = triangle[2];
     
          // compute the center of the triangle
          Vec3 center = (positions[a].getCenter() + positions[b].getCenter() + positions[c].getCenter())/3;
          vparams->drawTool()->drawArrow(center, center + vx, radius, type::RGBAColor(1.0, 0.0, 0.0, 1.0));
          vparams->drawTool()->drawArrow(center, center + vy, radius, type::RGBAColor(0.0, 1.0, 0.0, 1.0));
          vparams->drawTool()->drawArrow(center, center + vz, radius, type::RGBAColor(0.0, 0.0, 1.0, 1.0));
        }
        triangleInfo.endEdit();
    }

    // Draw triangles using shell triangle info
    if (d_showTriangle.getValue())
    {
        std::vector<sofa::type::Vec3> vertices;
        for (Size i = 0; i < nbTriangles; ++i)
        {
            const Triangle& tri = triangles[i];
            Index a = tri[0];
            Index b = tri[1];
            Index c = tri[2];

            auto vertex_vector = x[a];
            auto vertex_position = sofa::type::Vec3(vertex_vector[0], vertex_vector[1], vertex_vector[2]);
            vertices.push_back(sofa::type::Vec3(vertex_position));
            vertex_vector = x[b];
            vertex_position = sofa::type::Vec3(vertex_vector[0], vertex_vector[1], vertex_vector[2]);
            vertices.push_back(sofa::type::Vec3(vertex_position));
            vertex_vector = x[c];
            vertex_position = sofa::type::Vec3(vertex_vector[0], vertex_vector[1], vertex_vector[2]);
            vertices.push_back(sofa::type::Vec3(vertex_position));        
            }
        vparams->drawTool()->drawTriangles(vertices, sofa::type::RGBAColor(1, 0, 0, 0.9));
        vertices.clear();
    }

    // Visualize the stress on the shell
    if (d_showMeasuredValue.getValue())
    {        
        std::vector<sofa::type::Vec3> vertices;
        std::vector<sofa::type::RGBAColor> colorVector;

        auto evalColor = d_drawColorMap.getEvaluator(m_minStress, m_maxStress);
        const auto& measuredValues = sofa::helper::getReadAccessor(d_measuredValues);

        if(measuredValues.size()){
            for (Size i = 0; i < nbTriangles; ++i)
            {
                const Triangle& tri = triangles[i];
                Index a = tri[0];
                Index b = tri[1];
                Index c = tri[2];

                auto vertex_vector = x[a];
                auto vertex_position = sofa::type::Vec3(vertex_vector[0], vertex_vector[1], vertex_vector[2]);
                auto measured_value = measuredValues[ triangleInf[i].measure[0].id ];
                colorVector.push_back(evalColor(measured_value));
                vertices.push_back(vertex_position);

                vertex_vector = x[b];
                vertex_position = sofa::type::Vec3(vertex_vector[0], vertex_vector[1], vertex_vector[2]);
                measured_value = measuredValues[ triangleInf[i].measure[1].id ];
                colorVector.push_back(evalColor(measured_value));
                vertices.push_back(vertex_position);

                vertex_vector = x[c];
                vertex_position = sofa::type::Vec3(vertex_vector[0], vertex_vector[1], vertex_vector[2]);
                measured_value = measuredValues[ triangleInf[i].measure[2].id ]; 
                colorVector.push_back(evalColor(measured_value));
                vertices.push_back(vertex_position);

            }
            vparams->drawTool()->drawTriangles(vertices, colorVector);
            vertices.clear();
            colorVector.clear();
        }

    }


  }


} // namespace forcefield

} // namespace component

} // namespace sofa


#endif
