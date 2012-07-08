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
#ifndef SOFA_COMPONENT_FORCEFIELD_BEZIERSHELLFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_BEZIERSHELLFORCEFIELD_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/Data.h>

#include <sofa/component/component.h>

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/component/topology/TopologyData.h>

#include "../../controller/MeshInterpolator.h"
#include "../../engine/JoinMeshPoints.h"
#include "../fem/BezierShellInterpolation.h"


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

using namespace sofa::defaulttype;
using sofa::helper::vector;
using namespace sofa::component::topology;
using namespace sofa::core::behavior;

/// This class can be overridden if needed for additionnal storage within template specializations.
template<class DataTypes>
class BezierShellForceFieldInternalData
{
public:
};


template<class DataTypes>
class BezierShellForceField : public core::behavior::ForceField<DataTypes>
{
    public:
        SOFA_CLASS(SOFA_TEMPLATE(BezierShellForceField,DataTypes), SOFA_TEMPLATE(core::behavior::ForceField,DataTypes));

        typedef core::behavior::ForceField<DataTypes>       Inherited;
        typedef typename DataTypes::VecCoord                VecCoord;
        typedef typename DataTypes::VecDeriv                VecDeriv;

        typedef typename DataTypes::Coord                   Coord;
        typedef typename DataTypes::Deriv                   Deriv;
        typedef typename Coord::value_type                  Real;

        typedef Vec<3,Real> Vec3;
        typedef Vec<2,Real> Vec2;

        typedef Mat<2,2,Real> Mat22;
        typedef Mat<3,3,Real> Mat33;

        typedef helper::Quater<Real> Quat;

        typedef Data<VecCoord>                              DataVecCoord;
        typedef Data<VecDeriv>                              DataVecDeriv;

        typedef Vec3Types::VecCoord VecCoordHigh;

        typedef sofa::core::topology::BaseMeshTopology::index_type Index;
        typedef sofa::core::topology::BaseMeshTopology::Triangle Triangle;
        typedef sofa::core::topology::BaseMeshTopology::SeqTriangles SeqTriangles;

        typedef typename sofa::defaulttype::SolidTypes<Real>::Transform Transform;
        typedef typename sofa::defaulttype::SolidTypes<Real>::SpatialVector SpatialVector;

    protected:

        // Displacement vectorshfor in-plane forces:
        //  [ U1x, U1y, dT1z, U2x, U2y, dT2z, U3x, U3y, dT3z ]
        typedef Vec<9, Real> Displacement;
        // Displacement vector for bending forces:
        //  [ U1z, dT1x, dT1y, U2z, dT2x, dT2y, U3z, dT3x, dT3y ]
        typedef Vec<9, Real> DisplacementBending;
        typedef Mat<3, 3, Real> MaterialStiffness;              ///< matrix of material stiffness
        typedef Mat<3, 9, Real> StrainDisplacement;             ///< strain-displacement matrix for in-plane forces
        typedef Mat<3, 9, Real> StrainDisplacementBending;
        typedef Mat33  Transformation;                          ///< matrix for rigid transformations like rotations
        typedef Mat<9, 9, Real> StiffnessMatrix;
        typedef Mat<9, 9, Real> StiffnessMatrixBending;
        typedef Mat<18, 18, Real> StiffnessMatrixGlobalSpace;
        typedef Mat<4, 6, Real> GradDisplacement;


        sofa::core::topology::BaseMeshTopology* _topology;
        sofa::core::topology::BaseMeshTopology* _topologyTarget;

        // Nodes of the Bezier triangles

public:

        // Data for Gaussian quadrature
        static const int Gn = 6; // Number of Gauss points

        class TriangleInformation
        {
            public:

                helper::fixed_array <Vec3, 3> restLocalPositions;
#ifdef CRQUAT
                helper::fixed_array <Quat, 3> restLocalOrientationsInv;
#else
                helper::fixed_array <Transformation, 3> restLocalOrientationsInv;
#endif

                // Index of this element
                Index elementID;
                // Indices of each vertex
                Index a, b, c;

                // Corotational frame
                Vec3 frameCenter;
                Transformation frameOrientation;    // frame orientation    *
                Transformation frameOrientationInv; // it's inverse (transposition)
#ifdef CRQUAT
                Quat frameOrientationQ;             // representation as quaternion
#endif

                // Matrix for computing displacement gradient
                GradDisplacement gradU;

                // Matrix of interpolation functions
                // NOTE: we might need to always use double here, with
                // floats the matrix makes the strain-displacement and
                // stiffness matrices unusable due to lack of precision.
                Mat33 interpol;

                // Nodes of the Bezier triangle
                helper::fixed_array<Vec3, 10> pts;          // ... in local frame

                // the strain-displacement matrices at each Gauss point
                StrainDisplacement strainDisplacementMatrix1;
                StrainDisplacement strainDisplacementMatrix2;
                StrainDisplacement strainDisplacementMatrix3;
                StrainDisplacement strainDisplacementMatrix4;

                StrainDisplacement strainDisplacementMatrix[Gn];

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

        class TriangleHandler : public topology::TopologyDataHandler<topology::Triangle, sofa::helper::vector<TriangleInformation> >
        {
            typedef topology::TopologyDataHandler<topology::Triangle, sofa::helper::vector<TriangleInformation> > Inherited;
            public:
                TriangleHandler(BezierShellForceField<DataTypes>* _ff, TriangleData<sofa::helper::vector<TriangleInformation> >* _data) : Inherited(_data), ff(_ff) {}

                void applyCreateFunction(unsigned int triangleIndex, TriangleInformation& ,
                    const Triangle & t,
                    const sofa::helper::vector< unsigned int > &,
                    const sofa::helper::vector< double > &);

                void swap(unsigned int i1, unsigned int i2);

            protected:
                BezierShellForceField<DataTypes>* ff;
        };

        BezierShellForceField();

        virtual ~BezierShellForceField();
        virtual void init();
        virtual void reinit();
        virtual void addForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& dataF, const DataVecCoord& dataX, const DataVecDeriv& /*dataV*/ ) ;
        virtual void addDForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& datadF, const DataVecDeriv& datadX ) ;
        virtual void addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix);
        virtual void addBToMatrix(sofa::defaulttype::BaseMatrix * /*mat*/, double /*bFact*/, unsigned int &/*offset*/);
        //virtual double getPotentialEnergy(const VecCoord& x) const;
        virtual void handleTopologyChange();

        virtual void draw(const core::visual::VisualParams* vparams);

        sofa::core::topology::BaseMeshTopology* getTopology() {return _topology;}

        Data<Real> f_poisson;
        Data<Real> f_young;
        Data <Real> f_thickness;
        Data<unsigned int> f_polarMaxIters;
        Data<Real> f_polarMinTheta;

        // Allow transition between rest shapes
        SingleLink<BezierShellForceField<DataTypes>,
            sofa::component::controller::MeshInterpolator<DataTypes>,
            BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> restShape;

        // Indirect rest shape indexing (e.g. for "joining" two meshes)
        bool mapTopology;
        SingleLink<BezierShellForceField<DataTypes>,
            sofa::component::engine::JoinMeshPoints<DataTypes>,
            BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> topologyMapper;

        // Bezier shell interpolation
        SingleLink<BezierShellForceField<DataTypes>,
            sofa::component::fem::BezierShellInterpolation<DataTypes>,
            BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> bsInterpolation;



        static void computeEdgeBezierPoints(const Index& a, const Index& b, const Index& c,
            const VecCoord& x, const helper::vector<Vec3>& norms,
            helper::fixed_array<Vec3,10> &bezierPoints);

        void handleEvent(sofa::core::objectmodel::Event *event);

        const TriangleInformation& getTriangleInfo(Index t) {
            return triangleInfo.getValue()[t];
        }

protected :

        TriangleData< sofa::helper::vector<TriangleInformation> > triangleInfo;
        TriangleHandler* triangleHandler;

        /// Material stiffness matrices for plane stress and bending
        MaterialStiffness materialMatrix;
        MaterialStiffness materialMatrixBending;

        Real polarMinSinTheta;

        //unsigned int pditers;

        void initTriangleOnce(const int i, const Index&a, const Index&b, const Index&c);
        void initTriangle(const int i);

        void computeLocalTriangle(const Index elementIndex, bool bFast);

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

        // inPlane Gradient
        void computeInPlaneDisplacementGradient(GradDisplacement& gradU, const Vec3& GP, const TriangleInformation &tinfo);

        // Use polar decomposition to fix inPlane rotation
        void fixFramePolar(const Displacement &Disp, Mat22 &R, TriangleInformation &tinfo);

        /// f += Kx where K is the stiffness matrix and x a displacement
        virtual void applyStiffness(VecDeriv& f, const VecDeriv& dx, const Index elementIndex, const double kFactor);
        virtual void computeMaterialMatrix();

        //void bezierFunctions(const Vec2& baryCoord, sofa::helper::fixed_array<Real,10> &f_bezier);
        //void bezierDerivateFunctions(const Vec2& baryCoord, sofa::helper::fixed_array<Real,10> &df_dx_bezier, sofa::helper::fixed_array<Real,10> &df_dy_bezier);
        void interpolateRefFrame(TriangleInformation *tinfo, const Vec2& baryCoord);


        void accumulateForce(VecDeriv& f, const VecCoord & p, const Index elementIndex);

        void convertStiffnessMatrixToGlobalSpace(StiffnessMatrixGlobalSpace &K_gs, TriangleInformation *tinfo);
};


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // #ifndef SOFA_COMPONENT_FORCEFIELD_BEZIERSHELLFORCEFIELD_H
