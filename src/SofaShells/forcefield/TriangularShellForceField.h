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
#ifndef SOFA_COMPONENT_FORCEFIELD_TRIANGULAR_BENDING_FEM_FORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_TRIANGULAR_BENDING_FEM_FORCEFIELD_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif


#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/helper/OptionsGroup.h>

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/topology/TopologyData.h>

#include <sofa/helper/ColorMap.h>

// Uncomment the following to use quaternions instead of matrices for
// rotations. Quaternions are slightly faster but numericaly quite unstable
// (and in my oppinion unusable). I don't recommend that!
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

/// This class can be overridden if needed for additional storage within template specializations.
template<class DataTypes>
class TriangularShellForceFieldInternalData
{
public:
};


template<class DataTypes>
class TriangularShellForceField : public core::behavior::ForceField<DataTypes>
{
    public:
        SOFA_CLASS(SOFA_TEMPLATE(TriangularShellForceField,DataTypes), SOFA_TEMPLATE(core::behavior::ForceField,DataTypes));

        typedef core::behavior::ForceField<DataTypes>       Inherited;
        typedef typename DataTypes::VecCoord                VecCoord;
        typedef typename DataTypes::VecDeriv                VecDeriv;
        //typedef typename DataTypes::VecReal                 VecReal;

        typedef typename DataTypes::Coord                   Coord;
        typedef typename DataTypes::Deriv                   Deriv;
        typedef typename Coord::value_type                  Real;

        typedef Vec<2,Real> Vec2;
        typedef Vec<3,Real> Vec3;
        typedef Vec<9,Real> Vec9;

        typedef sofa::type::Quat<Real> Quat;

        typedef Data<VecCoord>                              DataVecCoord;
        typedef Data<VecDeriv>                              DataVecDeriv;

        typedef sofa::defaulttype::Vec3Types::VecCoord VecCoordHigh;

        typedef sofa::Index Index;
        typedef sofa::core::topology::BaseMeshTopology::Triangle Triangle;
        typedef sofa::core::topology::BaseMeshTopology::SeqTriangles SeqTriangles;

        class TriangleInformation;

        // Flags for the visualization
        Data<bool> d_showArrow;
        Data<Real> d_arrowRadius;
        Data<bool> d_showTriangle;
        Data<bool> d_showMeasuredValue;

    protected:

        typedef Vec<9, Real> Displacement;                      // the displacement vector
        typedef Mat<3, 3, Real> MaterialStiffness;              // the matrix of material stiffness
        typedef Mat<3, 9, Real> StrainDisplacement;             // the strain-displacement matrix
        typedef Mat<3, 3, Real> Transformation;                 // matrix for rigid transformations like rotations
        typedef Mat<9, 9, Real> StiffnessMatrix;                // element stiffness matrix
        typedef Mat<18, 18, Real> StiffnessMatrixFull;          // stiffness matrix for shell (= bending plate + membrane)

        typedef void (TriangularShellForceField<DataTypes>::*compstiff)(StiffnessMatrix &K, TriangleInformation &tinfo);
        typedef type::fixed_array<Real, 10> AndesBeta;

        sofa::core::topology::BaseMeshTopology* _topology;

public:

        class TriangleInformation
        {
            public:

                // Indices of each vertex
                Index a, b, c;

                // Rest position in local (in-plane) coordinates
                type::fixed_array <Vec3, 3> restPositions;
#ifdef CRQUAT
                type::fixed_array <Quat, 3> restOrientationsInv;
#else
                type::fixed_array <Transformation, 3> restOrientationsInv;
#endif

                // Deformed position in local (in-plane) coordinates
                type::fixed_array <Vec3, 3> deformedPositions;

                // Frame rotation as matrix and quaternion
                Transformation R, Rt;
#ifdef CRQUAT
                Quat Q;
#endif

                // The strain-displacement matrices at Gauss points
                StrainDisplacement strainDisplacementMatrixMembrane[4];
                StrainDisplacement strainDisplacementMatrixBending[4];

                // Stiffness matrix
                StiffnessMatrix stiffnessMatrixMembrane;
                StiffnessMatrix stiffnessMatrixBending;

                // Measure stress or strain
                struct MeasurePoint {
                    Vec3 point;             // Barycentric coordinates
                    StrainDisplacement B;   // Strain-displacement Matrix
                    StrainDisplacement Bb;  // Strain-displacement Matrix bending
                    Index id;               // Index into the result array
                };
                type::vector<MeasurePoint> measure;


                // The following are in rest shape
                // - element area
                Real area;
                // - directional vectors: 1-2, 2-3, 3-1
                type::fixed_array <Vec2, 3> d;
                // - squared lengths of 'd'
                type::fixed_array <Real, 3> l2;

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
                TRQSTriangleHandler(TriangularShellForceField<DataTypes>* _ff, TriangleData<sofa::type::vector<TriangleInformation> >* _data) : TopologyDataHandler<Triangle, sofa::type::vector<TriangleInformation> >(_data), ff(_ff) {}

                void applyCreateFunction(unsigned int triangleIndex, TriangleInformation& ,
                    const Triangle & t,
                    const sofa::type::vector< unsigned int > &,
                    const sofa::type::vector< double > &, const VecCoord& x0);

            protected:
                TriangularShellForceField<DataTypes>* ff;
        };

        TriangularShellForceField();

        virtual ~TriangularShellForceField();
        void init() override;
        void reinit() override;
        void addForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& dataF, const DataVecCoord& dataX, const DataVecDeriv& /*dataV*/ ) override ;
        void addDForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& datadF, const DataVecDeriv& datadX ) override ;
        void addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix) override;
        void draw(const core::visual::VisualParams* vparams) override;

        SReal getPotentialEnergy(const sofa::core::MechanicalParams* /*mparams*/, const DataVecCoord& /*x*/) const override { return 0; }

        sofa::core::topology::BaseMeshTopology* getTopology() {return _topology;}

        Data<Real> d_poisson;
        Data<Real> d_young;
        Data <Real> d_thickness;
        Data <sofa::helper::OptionsGroup> d_membraneElement;
        Data <sofa::helper::OptionsGroup> d_bendingElement;
        Data<bool> d_corotated;
        Data<sofa::helper::OptionsGroup> d_measure;
        Data<type::vector<Real> > d_measuredValues;
        Data<bool> d_isShellveryThin;
        Data<bool> d_use_rest_position;

        TRQSTriangleHandler* triangleHandler;

private:
        sofa::helper::ColorMap* p_drawColorMap; ///< colormap to display the gradiant of stress if @sa showStressValue is set to true
        Real m_minStress = 0; ///< min stress computed for @sa showStressValue
        Real m_maxStress = 0; ///< max stress computed for @sa showStressValue

protected :

        // Selected elements
        compstiff csMembrane;
        compstiff csBending;

        /// Material stiffness matrix
        MaterialStiffness materialMatrix, materialMatrixMembrane, materialMatrixBending;
        TriangleData< sofa::type::vector<TriangleInformation> > triangleInfo;

        // What to measure
        bool bMeasureStrain;
        bool bMeasureStress;

        void initTriangle(const int i, const Index&a, const Index&b, const Index&c, const VecCoord& x0);

        void computeRotation(Transformation& R, const VecCoord &x, const Index &a, const Index &b, const Index &c);
        void computeRotation(Transformation& R, const type::fixed_array<Vec3, 3> &x);
        void computeMaterialStiffness();

        void computeDisplacement(Displacement &Dm, Displacement &Db, const VecCoord &x, const Index elementIndex);
        void accumulateForce(VecDeriv& f, const VecCoord & p, const Index elementIndex);
        void computeStiffnessMatrixMembrane(StiffnessMatrix &K, TriangleInformation &tinfo);
        void computeStiffnessMatrixBending(StiffnessMatrix &K, TriangleInformation &tinfo);
        void computeForce(Displacement &Fm, const Displacement& Dm, Displacement &Fb, const Displacement& Db,const Index elementIndex);
        virtual void applyStiffness(VecDeriv& f, const VecDeriv& dx, const Index elementIndex, const double kFactor);

        void convertStiffnessMatrixToGlobalSpace(StiffnessMatrixFull &K_gs, const TriangleInformation &tinfo);

        // Membrane Elements
        void computeStiffnessMatrixCST(StiffnessMatrix &K, TriangleInformation &tinfo);
        void computeStiffnessMatrixAll3I(StiffnessMatrix &K, TriangleInformation &tinfo);
        void computeStiffnessMatrixAll3M(StiffnessMatrix &K, TriangleInformation &tinfo);
        void computeStiffnessMatrixAllLS(StiffnessMatrix &K, TriangleInformation &tinfo);
        void computeStiffnessMatrixLSTRet(StiffnessMatrix &K, TriangleInformation &tinfo);
        void computeStiffnessMatrixAndesOpt(StiffnessMatrix &K, TriangleInformation &tinfo);

        // Bending plate elements
        void computeStiffnessMatrixDKT(StiffnessMatrix &K, TriangleInformation &tinfo);

        // Helper functions for the elements
        void andesTemplate(StiffnessMatrix &K, const TriangleInformation &tinfo, const Real alpha, const AndesBeta &beta);
        void dktSD(StrainDisplacement &B, const TriangleInformation &tinfo, const Real xi, const Real eta);

};


#if defined(WIN32) && !defined(SOFA_COMPONENT_FORCEFIELD_TRIANGULAR_BENDING_FEM_FORCEFIELD_CPP)
#pragma warning(disable : 4231)
#ifndef SOFA_FLOAT
extern template class TriangularShellForceField<defaulttype::Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class TriangularShellForceField<defaulttype::Rigid3fTypes>;
#endif
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif
