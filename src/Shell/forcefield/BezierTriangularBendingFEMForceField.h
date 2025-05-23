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
#ifndef SOFA_COMPONENT_FORCEFIELD_BEZIER_TRIANGULAR_BENDING_FEM_FORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_BEZIER_TRIANGULAR_BENDING_FEM_FORCEFIELD_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif


#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/Data.h>

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/topology/TopologyData.h>

#include <Shell/controller/MeshInterpolator.h>
#include <Shell/engine/JoinMeshPoints.h>


// Uncomment the following to use quaternions instead of matrices for
// rotations. Quaternions are slightly faster but numericaly much, much *less*
// stable. I don't recommend that!
//#define CRQUAT


namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::type;
using sofa::type::vector;
using namespace sofa::core::topology;
using namespace sofa::core::behavior;

/// This class can be overridden if needed for additionnal storage within template specializations.
template<class DataTypes>
class BezierTriangularBendingFEMForceFieldInternalData
{
public:
};


template<class DataTypes>
class BezierTriangularBendingFEMForceField : public core::behavior::ForceField<DataTypes>
{
    public:
        SOFA_CLASS(SOFA_TEMPLATE(BezierTriangularBendingFEMForceField,DataTypes), SOFA_TEMPLATE(core::behavior::ForceField,DataTypes));

        typedef core::behavior::ForceField<DataTypes>       Inherited;
        typedef typename DataTypes::VecCoord                VecCoord;
        typedef typename DataTypes::VecDeriv                VecDeriv;

        typedef typename DataTypes::Coord                   Coord;
        typedef typename DataTypes::Deriv                   Deriv;
        typedef typename Coord::value_type                  Real;

        typedef Vec<3,Real> Vec3;
        typedef Vec<2,Real> Vec2;

        typedef Mat<3,3,Real> Mat33;

        typedef sofa::type::Quat<Real> Quat;

        typedef Data<VecCoord>                              DataVecCoord;
        typedef Data<VecDeriv>                              DataVecDeriv;

        typedef sofa::defaulttype::Vec3Types::VecCoord VecCoordHigh;

        typedef sofa::Index Index;
        typedef sofa::core::topology::BaseMeshTopology::Triangle Triangle;
        typedef sofa::core::topology::BaseMeshTopology::SeqTriangles SeqTriangles;

        class TriangleInformation;

    protected:

        // Displacement vector for in-plane forces:
        //  [ U1x, U1y, dT1z, U2x, U2y, dT2z, U3x, U3y, dT3z ]
        typedef Vec<9, Real> Displacement;
        // Displacement vector for bending forces:
        //  [ U1z, dT1x, dT1y, U2z, dT2x, dT2y, U3z, dT3x, dT3y ]
        typedef Vec<9, Real> DisplacementBending;
        typedef Mat<3, 3, Real> MaterialStiffness;              ///< matrix of material stiffness
        typedef Mat<3, 9, Real> StrainDisplacement;             ///< strain-displacement matrix for in-plane forces
        typedef Mat<3, 9, Real> StrainDisplacementBending;
        typedef Mat<3, 3, Real > Transformation;                ///< matrix for rigid transformations like rotations
        typedef Mat<9, 9, Real> StiffnessMatrix;
        typedef Mat<9, 9, Real> StiffnessMatrixBending;
        typedef Mat<18, 18, Real> StiffnessMatrixGlobalSpace;

        sofa::core::topology::BaseMeshTopology* _topology;
        sofa::core::topology::BaseMeshTopology* _topologyTarget;

        ////////////////////////// Inherited attributes ////////////////////////////
        /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
        /// Bring inherited attributes and function in the current lookup context.
        /// otherwise any access to the base::attribute would require
        /// the "this->" approach.
        using  ForceField<DataTypes>::d_componentState ;
        ////////////////////////////////////////////////////////////////////////////

public:

        class TriangleInformation
        {
            public:

                type::fixed_array <Vec3, 3> restLocalPositions;
#ifdef CRQUAT
                type::fixed_array <Quat, 3> restLocalOrientationsInv;
#else
                type::fixed_array <Transformation, 3> restLocalOrientationsInv;
#endif

                // Indices of each vertex
                Index a, b, c;

                // Indices in rest shape topology. Normaly is the same as a, b,
                // c, but is different if topologyMapper is set.
                Index a0, b0, c0;

                // Corotational frame
                Vec3 frameCenter;
                Transformation frameOrientation;    // frame orientation
                Transformation frameOrientationInv; // it's inverse (transposition)
#ifdef CRQUAT
                Quat frameOrientationQ;             // representation as quaternion
#endif

                // Matrix of interpolation functions
                Mat<3,3> interpol;

                // Segments in rest position used to keep Bézier points rigidly fixed
                Vec3 P0_P1_inFrame0;
                Vec3 P0_P2_inFrame0;
                Vec3 P1_P2_inFrame1;
                Vec3 P1_P0_inFrame1;
                Vec3 P2_P0_inFrame2;
                Vec3 P2_P1_inFrame2;

                // Nodes of the Bezier triangle
                type::fixed_array<Vec3, 10> bezierNodes;  // ... in global frame
                type::fixed_array<Vec3, 10> pts;          // ... in local frame
                // ... of the rest shape
                type::fixed_array<Vec3, 10> bezierNodes0; // ... in global frame

                // the strain-displacement matrices at each Gauss point
                StrainDisplacement strainDisplacementMatrix1;
                StrainDisplacement strainDisplacementMatrix2;
                StrainDisplacement strainDisplacementMatrix3;
                StrainDisplacement strainDisplacementMatrix4;

                // the strain-displacement matrices at each Gauss point
                StrainDisplacementBending strainDisplacementMatrixB1;
                StrainDisplacementBending strainDisplacementMatrixB2;
                StrainDisplacementBending strainDisplacementMatrixB3;
                StrainDisplacementBending strainDisplacementMatrixB4;

                // Stiffness matrix K = J * M * Jt
                StiffnessMatrix stiffnessMatrix;

                // Stiffness matrix for bending K = Jt * M * J
                StiffnessMatrixBending stiffnessMatrixBending;

                // Surface Area * 2
                Real area2;

                TriangleInformation() { }

                /// Output stream
                inline friend std::ostream& operator<< ( std::ostream& os, const TriangleInformation& /*ti*/ )
                {
                    return os;
                }

                /// Input stream
                inline friend std::istream& operator>> ( std::istream& in, TriangleInformation& /*ti*/ )
                {
                    return in;
                }
        };

        class TRQSTriangleHandler : public TopologyDataHandler<Triangle, sofa::type::vector<TriangleInformation> >
        {
            public:
                TRQSTriangleHandler(BezierTriangularBendingFEMForceField<DataTypes>* _ff, TriangleData<sofa::type::vector<TriangleInformation> >* _data)
                    : TopologyDataHandler<Triangle, sofa::type::vector<TriangleInformation> >(_data)
                    , ff(_ff)
                {
                }

                void applyCreateFunction(unsigned int triangleIndex, TriangleInformation& ,
                    const Triangle & t,
                    const sofa::type::vector< unsigned int > &,
                    const sofa::type::vector< double > &);

            protected:
                BezierTriangularBendingFEMForceField<DataTypes>* ff;
        };

        BezierTriangularBendingFEMForceField();

        virtual ~BezierTriangularBendingFEMForceField();
        void init() override;
        void reinit() override;
        void addForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& dataF, const DataVecCoord& dataX, const DataVecDeriv& /*dataV*/ ) override ;
        void addDForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& datadF, const DataVecDeriv& datadX ) override ;
        void addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix) override;
        void addBToMatrix(sofa::linearalgebra::BaseMatrix * /*mat*/, double /*bFact*/, unsigned int &/*offset*/) override;
        void handleTopologyChange() override;

        SReal getPotentialEnergy(const sofa::core::MechanicalParams* /*mparams*/, const DataVecCoord& /*x*/) const override { return 0; }

        void draw(const core::visual::VisualParams* vparams) override;

        sofa::core::topology::BaseMeshTopology* getTopology() {return _topology;}

        Data<Real> f_poisson;
        Data<Real> f_young;
        Data <Real> f_thickness;
        Data< type::vector<Vec3> > normals;

        // Allow transition between rest shapes
        SingleLink<BezierTriangularBendingFEMForceField<DataTypes>,
            shell::controller::MeshInterpolator<DataTypes>,
            BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> restShape;

        // Indirect rest shape indexing (e.g. for "joining" two meshes)
        bool mapTopology;
        SingleLink<BezierTriangularBendingFEMForceField<DataTypes>,
            shell::engine::JoinMeshPoints<DataTypes>,
            BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> topologyMapper;


        static void computeEdgeBezierPoints(const Index& a, const Index& b, const Index& c,
            const VecCoord& x, const type::vector<Vec3>& norms,
            type::fixed_array<Vec3,10> &bezierPoints);

        void handleEvent(sofa::core::objectmodel::Event *event) override;

        const TriangleInformation& getTriangleInfo(Index t) {
            return triangleInfo.getValue()[t];
        }

protected :

        TriangleData< sofa::type::vector<TriangleInformation> > triangleInfo;
        TRQSTriangleHandler* triangleHandler;

        /// Material stiffness matrices for plane stress and bending
        MaterialStiffness materialMatrix;
        MaterialStiffness materialMatrixBending;


        void initTriangleOnce(const int i, const Index&a, const Index&b, const Index&c);
        void initTriangle(const int i);

        void computeLocalTriangle(const VecCoord &x, const Index elementIndex);

        void computeDisplacements( Displacement &Disp, DisplacementBending &BDisp, const VecCoord &x, TriangleInformation *tinfo);
        void computeStrainDisplacementMatrixMembrane(TriangleInformation &tinfo);
        void computeStrainDisplacementMatrixBending(TriangleInformation &tinfo);
        void computeStiffnessMatrixMembrane(StiffnessMatrix &K, const TriangleInformation &tinfo);
        void computeStiffnessMatrixBending(StiffnessMatrixBending &K, const TriangleInformation &tinfo);
        void computeForceMembrane(Displacement &F, const Displacement& D, const Index elementIndex);
        void computeForceBending(DisplacementBending &F, const DisplacementBending& D, const Index elementIndex);

        // Strain-displacement matrices
        void matrixSDM(StrainDisplacement &J, const Vec3 &GP, const TriangleInformation& tinfo);
        void matrixSDB(StrainDisplacementBending &J, const Vec3 &GP, const TriangleInformation& tinfo);

        /// f += Kx where K is the stiffness matrix and x a displacement
        virtual void applyStiffness(VecDeriv& f, const VecDeriv& dx, const Index elementIndex, const double kFactor);
        virtual void computeMaterialMatrix();

        void computePosBezierPoint(const TriangleInformation *tinfo, const Index& a, const Index& b, const Index& c, const VecCoord& x, sofa::type::fixed_array<Vec3,10> &X_bezierPoints);
        void bezierFunctions(const Vec2& baryCoord, sofa::type::fixed_array<Real,10> &f_bezier);
        void bezierDerivateFunctions(const Vec2& baryCoord, sofa::type::fixed_array<Real,10> &df_dx_bezier, sofa::type::fixed_array<Real,10> &df_dy_bezier);
        void interpolateRefFrame(TriangleInformation *tinfo, const Vec2& baryCoord, const Index& a, const Index& b, const Index& c, const VecCoord& x, sofa::type::fixed_array<Vec3,10>& X_bezierPoints );


        void accumulateForce(VecDeriv& f, const VecCoord & p, const Index elementIndex);

        void convertStiffnessMatrixToGlobalSpace(StiffnessMatrixGlobalSpace &K_gs, TriangleInformation *tinfo);
};


#if defined(WIN32) && !defined(SOFA_COMPONENT_FORCEFIELD_BEZIER_TRIANGULAR_BENDING_FEM_FORCEFIELD_CPP)
#pragma warning(disable : 4231)
#ifndef SOFA_FLOAT
extern template class BezierTriangularBendingFEMForceField<defaulttype::Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class BezierTriangularBendingFEMForceField<defaulttype::Rigid3fTypes>;
#endif
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // #ifndef SOFA_COMPONENT_FORCEFIELD_BEZIER_TRIANGULAR_BENDING_FEM_FORCEFIELD_H
