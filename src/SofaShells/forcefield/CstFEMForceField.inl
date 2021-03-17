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
#ifndef SOFA_COMPONENT_FORCEFIELD_CST_FEM_FORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_CST_FEM_FORCEFIELD_INL

#include <SofaShells/forcefield/CstFEMForceField.h>
#include <SofaBaseTopology/TopologyData.inl>
#include <sofa/helper/rmath.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>

#ifdef _WIN32
#include <windows.h>
#endif


namespace sofa
{
	namespace component
	{
		namespace forcefield
		{
			using namespace sofa::defaulttype;
			using namespace	sofa::component::topology;

// ----------------------------------------------------------------------------
// ---  Topology Creation/Destruction functions
// ----------------------------------------------------------------------------
template< class DataTypes>
void CstFEMForceField<DataTypes>::TRQSTriangleHandler::applyCreateFunction(unsigned int triangleIndex, TriangleInformation &, const Triangle &t, const sofa::helper::vector<unsigned int> &, const sofa::helper::vector<double> &)
{
    if (ff)
    {
        ff->initTriangle(triangleIndex, t[0], t[1], t[2]);
    }
}


// ----------------------------------------------------------------------------
// --- Constructor
// ----------------------------------------------------------------------------
template <class DataTypes>
CstFEMForceField<DataTypes>::CstFEMForceField()
: f_poisson(initData(&f_poisson,(Real)0.45,"poissonRatio","Poisson ratio in Hooke's law"))
, f_young(initData(&f_young,(Real)3000.,"youngModulus","Young modulus in Hooke's law"))
, f_thickness(initData(&f_thickness,(Real)0.1,"thickness","Thickness of the plates"))
, f_corotated(initData(&f_corotated, true, "corotated", "Compute forces in corotational frame"))
//, f_measure(initData(&f_measure, "measure", "Compute the strain or stress"))
//, f_measuredValues(initData(&f_measuredValues, "measuredValues", "Measured values for stress or strain"))
, f_stiffnessFactor(initData(&f_stiffnessFactor, Real(1.0), "stiffnessFactor", "stiffness factor between 0 and 1 to reduce the weight of the forcefield"))
, triangleInfo(initData(&triangleInfo, "triangleInfo", "Internal triangle data"))

{
    //f_measure.beginEdit()->setNames(3,
    //    "None",                 // Draw nothing
    //    "Strain (norm)",        // L_2 norm of strain in x and y directions
    //    "Von Mises stress"      // Von Mises stress criterion
    //    );
    //f_measure.beginEdit()->setSelectedItem("None");
    //f_measure.endEdit();

    triangleHandler = new TRQSTriangleHandler(this, &triangleInfo);
}


// ----------------------------------------------------------------------------
// --- Destructor
// ----------------------------------------------------------------------------
template <class DataTypes>
CstFEMForceField<DataTypes>::~CstFEMForceField()
{
    if(triangleHandler) delete triangleHandler;
}

// ----------------------------------------------------------------------------
// --- Initialization stage
// ----------------------------------------------------------------------------
template <class DataTypes>
void CstFEMForceField<DataTypes>::init()
{
    this->Inherited::init();

    _topology = this->getContext()->getMeshTopology();

    if (_topology->getNbTriangles()==0)
    {
            serr << "CstFEMForceField: object must have a Triangular Set Topology."<<sendl;
            return;
    }

    // Create specific handler for TriangleData
    triangleInfo.createTopologicalEngine(_topology, triangleHandler);
    triangleInfo.registerTopologicalData();

    reinit();
}


// ----------------------------------------------------------------------------
// --- Re-initialization (called when we change a parameter through the GUI)
// ----------------------------------------------------------------------------
template <class DataTypes> void CstFEMForceField<DataTypes>::reinit()
{
    helper::vector<TriangleInformation>& ti = *(triangleInfo.beginEdit());

    // Prepare material matrices
    computeMaterialStiffness();

    //// What to compute?
    //if (f_measure.getValue().getSelectedItem() == "None") {
    //    bMeasureStrain = false;  bMeasureStress = false;
    //} else if (f_measure.getValue().getSelectedItem() == "Strain (norm)") {
    //    bMeasureStrain = true;  bMeasureStress = false;
    //} else if (f_measure.getValue().getSelectedItem() == "Von Mises stress") {
    //    bMeasureStrain = false;  bMeasureStress = true;
    //} else {
    //    serr << "Invalid value for measure'" << f_measure.getValue().getSelectedItem() << "'" << sendl;
    //    return;
    //}

    //if (bMeasureStrain || bMeasureStress)
    //{
    //    f_measuredValues.beginEdit()->resize(_topology->getNbPoints());
    //    f_measuredValues.endEdit();
    //}

    /// Prepare to store info in the triangle array
    ti.resize(_topology->getNbTriangles());

    for (sofa::Index i=0; i<_topology->getNbTriangles(); ++i)
    {

        triangleHandler->applyCreateFunction(i, ti[i], _topology->getTriangle(i),
            (const sofa::helper::vector< unsigned int >)0, (const sofa::helper::vector< double >)0);
    }

    triangleInfo.endEdit();    
}


template <class DataTypes>
void CstFEMForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& dataF, const DataVecCoord& dataX, const DataVecDeriv& /*dataV*/ )
{    
    VecDeriv& f        = *(dataF.beginEdit());
    const VecCoord& p  =   dataX.getValue()  ;

    //std::cout << "--addForce" << std::endl;

    int nbTriangles=_topology->getNbTriangles();
    f.resize(p.size());

    for (int i=0; i<nbTriangles; i++)
    {
        accumulateForce(f, p, i);
    }

    dataF.endEdit();
}


template <class DataTypes>
void CstFEMForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& datadF, const DataVecDeriv& datadX )
{    
    VecDeriv& df        = *(datadF.beginEdit());
    const VecDeriv& dp  =   datadX.getValue()  ;

    Real kFactor = (Real)sofa::core::mechanicalparams::kFactor(mparams);

    //std::cout << "--addDForce" << std::endl;

    int nbTriangles=_topology->getNbTriangles();
    df.resize(dp.size());

    for (int i=0; i<nbTriangles; i++)
    {
        applyStiffness(df, dp, i, kFactor);
    }

    datadF.endEdit();
}

// ----------------------------------------------------------------------------
// --- Store the initial position of the nodes
// ----------------------------------------------------------------------------
template <class DataTypes>
void CstFEMForceField<DataTypes>::initTriangle(const int i, const Index&a, const Index&b, const Index&c)
{
    helper::vector<TriangleInformation>& ti = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &ti[i];

    // Store indices of each vertex
    tinfo->a = a;
    tinfo->b = b;
    tinfo->c = c;

    // Gets vertices of rest positions
    const VecCoord& x0 = this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();

    // Rotation from global to local frame
    Transformation R0;

    // Rest positions expressed in the local frame
    tinfo->restPositions[0] = x0[a];
    tinfo->restPositions[1] = x0[b];
    tinfo->restPositions[2] = x0[c];

    computeRotation(R0, tinfo->restPositions);
    tinfo->R = R0;
    tinfo->Rt.transpose(R0);

    tinfo->restPositions[0] = R0 * tinfo->restPositions[0];
    tinfo->restPositions[1] = R0 * tinfo->restPositions[1];
    tinfo->restPositions[2] = R0 * tinfo->restPositions[2];

    // Do some precomputations
    // - directional vectors
    tinfo->d[0] = tinfo->restPositions[0] - tinfo->restPositions[1];
    tinfo->d[1] = tinfo->restPositions[1] - tinfo->restPositions[2];
    tinfo->d[2] = tinfo->restPositions[2] - tinfo->restPositions[0];

    // - triangle area
    tinfo->area = helper::rabs(tinfo->d[2][0]*(-tinfo->d[0][1]) - (-tinfo->d[0][0])*tinfo->d[2][1])/2.0;

    // Compute stiffness matrix
    computeStiffnessMatrix(tinfo->stiffnessMatrix, *tinfo);
    //std::cout << "Km^e=" << tinfo->stiffnessMatrix << std::endl;

    triangleInfo.endEdit();
}


template <class DataTypes>
void CstFEMForceField<DataTypes>::computeRotation(Transformation& R, const helper::fixed_array<Vec3, 3> &x)
{
    if (!f_corotated.getValue()) {
        // Return identity matrix
        R = Transformation(Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,1));
        return;
    }
    
    // First vector on first edge
    // Second vector in the plane of the two first edges
    // Third vector orthogonal to first and second

    Vec3 edgex = x[1] - x[0];
    Vec3 edgey = x[2] - x[0];

    Vec3 edgez = cross(edgex, edgey);

    edgey = cross(edgez, edgex);

    edgex.normalize();
    edgey.normalize();
    edgez.normalize();

    R = Transformation(edgex, edgey, edgez);
}


// ----------------------------------------------------------------------------
// ---  Compute material stiffness
// ----------------------------------------------------------------------------
template <class DataTypes>
void CstFEMForceField<DataTypes>::computeMaterialStiffness()
{
    Real E = f_young.getValue(),
         nu = f_poisson.getValue(),
         t = f_thickness.getValue();

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

    // Integrate through the shell thickness
    materialMatrixMembrane *= t;
}

// ----------------------------------------------------------------------------
// --- Compute the stiffness matrix for membrane element
// ----------------------------------------------------------------------------
template <class DataTypes>
void CstFEMForceField<DataTypes>::computeStiffnessMatrix(StiffnessMatrix &K, TriangleInformation &tinfo)
{
    StrainDisplacement L;
    StrainDisplacementT Lt;
    L.clear();

    L[0][0] = tinfo.d[1][1];
    L[0][1] = 0;
    L[0][2] = tinfo.d[2][1];
    L[0][3] = 0;
    L[0][4] = tinfo.d[0][1];
    L[0][5] = 0;

    L[1][0] = 0;
    L[1][1] = -tinfo.d[1][0];
    L[1][2] = 0;
    L[1][3] = -tinfo.d[2][0];
    L[1][4] = 0;
    L[1][5] = -tinfo.d[0][0];

    L[2][0] = -tinfo.d[1][0];
    L[2][1] = tinfo.d[1][1];
    L[2][2] = -tinfo.d[2][0];
    L[2][3] = tinfo.d[2][1];
    L[2][4] = -tinfo.d[0][0];
    L[2][5] = tinfo.d[0][1];

    // Now should be:
    //  L *= h/2;
    //  K = 1/(h*Area) * L * E * L^T;
    // But we may simplify this a little and we also have 'h' already in E

    Lt.transpose(L);

    K = Lt * materialMatrixMembrane * L / (4.0*tinfo.area);

    // Compute strain-displacement matrix
    //for (unsigned int i=0; i< tinfo.measure.size(); i++) {
    //    tinfo.measure[i].B = Lt / (2*tinfo.area);
    //}


    if (this->f_printLog.getValue())
    {
        Displacement u = Vec<6,Real>(1, -5, 1, -5, 1, -5);
        Vec<6,Real> f = K*u;
        if (helper::rabs(f.sum()) > 1e-12) {
            sout << "Area = " << tinfo.area << sendl;
            sout << "a / b / c = " <<
                tinfo.restPositions[0] << " / " <<
                tinfo.restPositions[1] << " / " <<
                tinfo.restPositions[2] << sendl;
            sout << "d = " << tinfo.d << sendl;
            sout << "Km = " << K << sendl;
            sout << "-- Disp test Km (u=" << u << ")" <<
                " : " << f << " ... should be zero" << sendl;
        }
    }
}



// ----------------------------------------------------------------------------
// --- Compute displacement vector D as the difference between current current
// --- position and initial position.
// ----------------------------------------------------------------------------
template <class DataTypes>
void CstFEMForceField<DataTypes>::computeDisplacement(Displacement &D, const VecCoord &x, const Index elementIndex)
{
    helper::vector<TriangleInformation>& ti = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &ti[elementIndex];

    Index a = tinfo->a;
    Index b = tinfo->b;
    Index c = tinfo->c;

    // Compute local (in-plane) postions
    tinfo->deformedPositions[0] = x[a];
    tinfo->deformedPositions[1] = x[b];
    tinfo->deformedPositions[2] = x[c];

    // Compute rotation to local (in-plane) frame
    computeRotation(tinfo->R, tinfo->deformedPositions);
    tinfo->Rt.transpose(tinfo->R);

    // Compute local (in-plane) postions
    tinfo->deformedPositions[0] = tinfo->R * tinfo->deformedPositions[0];
    tinfo->deformedPositions[1] = tinfo->R * tinfo->deformedPositions[1];
    tinfo->deformedPositions[2] = tinfo->R * tinfo->deformedPositions[2];

    // Displacements
    Vec3 uA = tinfo->deformedPositions[0] - tinfo->restPositions[0];
    Vec3 uB = tinfo->deformedPositions[1] - tinfo->restPositions[1];
    Vec3 uC = tinfo->deformedPositions[2] - tinfo->restPositions[2];

    // Membrane
    D[0] = uA[0];
    D[1] = uA[1];

    D[2] = uB[0];
    D[3] = uB[1];

    D[4] = uC[0];
    D[5] = uC[1];

    triangleInfo.endEdit();
}


template <class DataTypes>
void CstFEMForceField<DataTypes>::accumulateForce(VecDeriv &f, const VecCoord &x, const Index elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Get the indices of the 3 vertices for the current triangle
    const Index& a = tinfo->a;
    const Index& b = tinfo->b;
    const Index& c = tinfo->c;

    // Compute in-plane displacements
    Displacement D;
    computeDisplacement(D, x, elementIndex);

    // Compute the membrane and bending plate forces on this element
    Displacement F;
    computeForce(F, D, elementIndex);

    // Compute the measure (stress/strain)
    //if (bMeasureStrain) {
    //    helper::vector<Real> &values = *f_measuredValues.beginEdit();
    //    for (unsigned int i=0; i< tinfo->measure.size(); i++) {
    //        Vec3 strain = tinfo->measure[i].B * Dm + tinfo->measure[i].Bb * Db;
    //        // Norm from strain in x and y
    //        // NOTE: Shear strain is not included
    //        values[ tinfo->measure[i].id ] = helper::rsqrt(
    //            strain[0] * strain[0] + strain[1] * strain[1]);
    //    }
    //    f_measuredValues.endEdit();
    //} else if (bMeasureStress) {
    //    helper::vector<Real> &values = *f_measuredValues.beginEdit();
    //    for (unsigned int i=0; i< tinfo->measure.size(); i++) {
    //        Vec3 stress = materialMatrix * tinfo->measure[i].B * Dm
    //            + materialMatrix * tinfo->measure[i].Bb * Db;
    //        // Von Mises stress criterion (plane stress)
    //        values[ tinfo->measure[i].id ] = helper::rsqrt(
    //              stress[0] * stress[0] - stress[0] * stress[1]
    //            + stress[1] * stress[1] + 3 * stress[2] * stress[2]);
    //    }
    //    f_measuredValues.endEdit();
    //}

    // Transform forces back into global frame
    f[a] -= Deriv(tinfo->Rt * Vec3(F[0], F[1], 0))*f_stiffnessFactor.getValue();
    f[b] -= Deriv(tinfo->Rt * Vec3(F[2], F[3], 0))*f_stiffnessFactor.getValue();
    f[c] -= Deriv(tinfo->Rt * Vec3(F[4], F[5], 0))*f_stiffnessFactor.getValue();

    triangleInfo.endEdit();
}


// ----------------------------------------------------------------------------
// ---  Compute force F = K * u
// ----------------------------------------------------------------------------
template <class DataTypes>
void CstFEMForceField<DataTypes>::computeForce(Displacement &F, const Displacement& D, const Index elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation &tinfo = triangleInf[elementIndex];

    // Compute forces
    F = tinfo.stiffnessMatrix * D;

    triangleInfo.endEdit();
}


// ----------------------------------------------------------------------------
// ---
// ----------------------------------------------------------------------------
template <class DataTypes>
void CstFEMForceField<DataTypes>::applyStiffness(VecDeriv& v, const VecDeriv& dx, const Index elementIndex, const double kFactor)
{
    helper::vector<TriangleInformation>& ti = *(triangleInfo.beginEdit());
    TriangleInformation &tinfo = ti[elementIndex];

    // Computes displacement
    Displacement D;
    Vec3 x_a, x_b, x_c;

    x_a = tinfo.R * dx[tinfo.a];
    x_b = tinfo.R * dx[tinfo.b];
    x_c = tinfo.R * dx[tinfo.c];

    D[0] = x_a[0];
    D[1] = x_a[1];

    D[2] = x_b[0];
    D[3] = x_b[1];

    D[4] = x_c[0];
    D[5] = x_c[1];

    // Compute dF
    Displacement dF;
    dF = tinfo.stiffnessMatrix * D;

    // Transform into global frame
    v[tinfo.a] -= Deriv(tinfo.Rt * Vec3(dF[0], dF[1], 0)) * kFactor * f_stiffnessFactor.getValue();
    v[tinfo.b] -= Deriv(tinfo.Rt * Vec3(dF[2], dF[3], 0)) * kFactor * f_stiffnessFactor.getValue();
    v[tinfo.c] -= Deriv(tinfo.Rt * Vec3(dF[4], dF[5], 0)) * kFactor * f_stiffnessFactor.getValue();

    triangleInfo.endEdit();
}

//#define PRINT

template<class DataTypes>
void CstFEMForceField<DataTypes>::addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
    StiffnessMatrixFull Kg;

    // Build Matrix Block for this ForceField
    unsigned int i, j ,n1, n2, row, column, ROW, COLUMN;
    Index node1, node2;

    sofa::core::behavior::MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

    Real kFactor = (Real)sofa::core::mechanicalparams::kFactor(mparams);

#ifdef PRINT
    r.matrix->clear();
    std::cout << "Initial global matrix (" << r.matrix->rowSize() << "x" << r.matrix->colSize() << ")" <<
        " kFactor=" << kFactor << std::endl;
    for (unsigned int i=0; i<r.matrix->rowSize(); i++)
    {
        for (unsigned int j=0; j<r.matrix->colSize(); j++)
        {
            std::cout << r.matrix->element(i,j) << ",";
            //std::cout << -r.matrix->element(i,j) / kFactor << ",";
        }
        std::cout << std::endl;
    }
#endif

    for(sofa::Index t=0 ; t != _topology->getNbTriangles() ; ++t)
    {
        const TriangleInformation &tinfo = triangleInf[t];
        const Triangle triangle = _topology->getTriangle(t);

        convertStiffnessMatrixToGlobalSpace(Kg, tinfo);

        // First node
        for (n1=0; n1<3; n1++)
        {
            node1 = triangle[n1];

            for(i=0; i<3; i++)
            {
                ROW = r.offset+3*node1+i;
                row = 3*n1+i;

                // Second node
                for (n2=0; n2<3; n2++)
                {
                    node2 = triangle[n2];

                    for (j=0; j<3; j++)
                    {
                        COLUMN = r.offset+3*node2+j;
                        column = 3*n2+j;

                        r.matrix->add(ROW, COLUMN, - Kg[row][column] * kFactor);
                    }
                }
            }
        }
    }

#ifdef PRINT
    std::cout << "Global matrix (" << r.matrix->rowSize() << "x" << r.matrix->colSize() << ")" <<
        " kFactor=" << kFactor << std::endl;
    for (unsigned int i=0; i<r.matrix->rowSize(); i++)
    {
        for (unsigned int j=0; j<r.matrix->colSize(); j++)
        {
            std::cout << r.matrix->element(i,j) << ",";
            //std::cout << -r.matrix->element(i,j) / kFactor << ",";
        }
        std::cout << std::endl;
    }
#endif

    triangleInfo.endEdit();
}



template<class DataTypes>
void CstFEMForceField<DataTypes>::convertStiffnessMatrixToGlobalSpace(StiffnessMatrixFull &Kg, const TriangleInformation &tinfo)
{
    StiffnessMatrixFull K1;
    unsigned int ig, jg;


    // Copy the stiffness matrix
    const StiffnessMatrix &K = tinfo.stiffnessMatrix;
    for (unsigned int bx=0; bx<3; bx++)
    {
        // Global row index
        ig = 3*bx;

        for (unsigned int by=0; by<3; by++)
        {
            // Global column index
            jg = 3*by;

            // X
            K1[ig+0][jg+0] = K[2*bx+0][2*by+0]; // X
            K1[ig+0][jg+1] = K[2*bx+0][2*by+1]; // Y

            // Y
            K1[ig+1][jg+0] = K[2*bx+1][2*by+0]; // X
            K1[ig+1][jg+1] = K[2*bx+1][2*by+1]; // Y
        }
    }

    // Extend rotation matrix and its transpose
    StiffnessMatrixFull R, Rt;

    for(unsigned int i=0;i<3;++i)
    {
        for(unsigned int j=0;j<3;++j)
        {
            R[i][j] = R[i+3][j+3] = R[i+6][j+6] = tinfo.R[i][j];
            Rt[i][j] = Rt[i+3][j+3] = Rt[i+6][j+6] = tinfo.Rt[i][j];
        }
    }

    // Transform stifness matrix into the global frame
    Kg = Rt * K1 * R;
    //std::cout << "Rorig=" << tinfo.R << " -- " << tinfo.Rt << std::endl;
    //std::cout << "K=" << K << std::endl;
    //std::cout << "K1=" << K1 << std::endl;
    //std::cout << "R=" << R << " -- " << Rt << std::endl;
    //std::cout << "Kg=" << Kg << std::endl;
}

} // namespace forcefield

} // namespace component

} // namespace sofa


#endif
