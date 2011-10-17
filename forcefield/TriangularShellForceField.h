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

#include <sofa/component/component.h>

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/component/topology/TriangleData.h>
#include <sofa/core/objectmodel/ObjectRef.h>



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

        typedef Data<VecCoord>                              DataVecCoord;
        typedef Data<VecDeriv>                              DataVecDeriv;

        typedef Vec3Types::VecCoord VecCoordHigh;

        typedef sofa::core::topology::BaseMeshTopology::index_type Index;
        typedef sofa::core::topology::BaseMeshTopology::Triangle Triangle;
        typedef sofa::core::topology::BaseMeshTopology::SeqTriangles SeqTriangles;

        class TriangleInformation;

    protected:

        typedef Vec<9, Real> Displacement;                      // the displacement vector
        typedef Mat<3, 3, Real> MaterialStiffness;              // the matrix of material stiffness
        typedef Mat<3, 9, Real> StrainDisplacement;             // the strain-displacement matrix
        typedef Mat<3, 3, Real> Transformation;                 // matrix for rigid transformations like rotations
        typedef Mat<9, 9, Real> StiffnessMatrix;                // element stiffness matrix
        typedef Mat<18, 18, Real> StiffnessMatrixFull;          // stiffness matrix for shell (= bending plate + membrane)

        typedef void (TriangularShellForceField<DataTypes>::*compstiff)(StiffnessMatrix &K, const TriangleInformation &tinfo);
        typedef helper::fixed_array<Real, 10> AndesBeta;

        sofa::core::topology::BaseMeshTopology* _topology;

public:

        class TriangleInformation
        {
            public:

                // Indices of each vertex
                Index a, b, c;

                // Rest position in local (in-plane) coordinates
                helper::fixed_array <Vec3, 3> restPositions;
                helper::fixed_array <Quat, 3> restOrientations;

                // Deformed position in local (in-plane) coordinates
                helper::fixed_array <Vec3, 3> deformedPositions;

                // Transformation rotation
                Transformation R, Rt;
                //Quat Qframe;

                // The strain-displacement matrices at Gauss points
                StrainDisplacement strainDisplacementMatrixMembrane[4];
                StrainDisplacement strainDisplacementMatrixBending[4];

                // Stiffness matrix
                StiffnessMatrix stiffnessMatrixMembrane;
                StiffnessMatrix stiffnessMatrixBending;

                // The following are in rest shape
                // - element area
                Real area;
                // - directional vectors: 1-2, 2-3, 3-1
                helper::fixed_array <Vec2, 3> d;
                // - squared lengths of 'd'
                helper::fixed_array <Real, 3> l2;

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

        TriangularShellForceField();

        virtual ~TriangularShellForceField();
        virtual void init();
        virtual void reinit();
        virtual void addForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& dataF, const DataVecCoord& dataX, const DataVecDeriv& /*dataV*/ ) ;
        virtual void addDForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& datadF, const DataVecDeriv& datadX ) ;
        virtual void addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix);
        virtual void addBToMatrix(sofa::defaulttype::BaseMatrix * /*mat*/, double /*bFact*/, unsigned int &/*offset*/);
        ////virtual void handleTopologyChange();

        virtual void draw();

        sofa::core::topology::BaseMeshTopology* getTopology() {return _topology;}
        //TriangleData<TriangleInformation>& getTriangleInfo() {return triangleInfo;}

        Data<Real> f_poisson;
        Data<Real> f_young;
        Data <Real> f_thickness;
        Data <sofa::helper::OptionsGroup> f_membraneElement;
        Data <sofa::helper::OptionsGroup> f_bendingElement;
        Data<bool> f_corotated;

protected :

        // Selected elements
        compstiff csMembrane;
        compstiff csBending;

        /// Material stiffness matrix
        MaterialStiffness materialMatrixMembrane, materialMatrixBending;
        TriangleData<TriangleInformation> triangleInfo;

        void initTriangle(const int i, const Index&a, const Index&b, const Index&c);

        static void TriangleCreationFunction (int , void* , TriangleInformation &, const Triangle& , const sofa::helper::vector< unsigned int > &, const sofa::helper::vector< double >&);

        void computeRotation(Transformation& R, const VecCoord &x, const Index &a, const Index &b, const Index &c);
        void computeRotation(Transformation& R, const helper::fixed_array<Vec3, 3> &x);
        void computeMaterialStiffness();

        void computeDisplacement(Displacement &Dm, Displacement &Db, const VecCoord &x, const Index elementIndex);
        void accumulateForce(VecDeriv& f, const VecCoord & p, const Index elementIndex);
        void computeStiffnessMatrixMembrane(StiffnessMatrix &K, TriangleInformation &tinfo);
        void computeStiffnessMatrixBending(StiffnessMatrix &K, const TriangleInformation &tinfo);
        void computeForce(Displacement &Fm, const Displacement& Dm, Displacement &Fb, const Displacement& Db,const Index elementIndex);
        virtual void applyStiffness(VecDeriv& f, const VecDeriv& dx, const Index elementIndex, const double kFactor);

        void convertStiffnessMatrixToGlobalSpace(StiffnessMatrixFull &K_gs, const TriangleInformation &tinfo);

        // Membrane Elements
        void computeStiffnessMatrixCST(StiffnessMatrix &K, const TriangleInformation &tinfo);
        void computeStiffnessMatrixAll3I(StiffnessMatrix &K, const TriangleInformation &tinfo);
        void computeStiffnessMatrixAll3M(StiffnessMatrix &K, const TriangleInformation &tinfo);
        void computeStiffnessMatrixAllLS(StiffnessMatrix &K, const TriangleInformation &tinfo);
        void computeStiffnessMatrixLSTRet(StiffnessMatrix &K, const TriangleInformation &tinfo);
        void computeStiffnessMatrixAndesOpt(StiffnessMatrix &K, const TriangleInformation &tinfo);

        // Bending plate elements
        void computeStiffnessMatrixDKT(StiffnessMatrix &K, const TriangleInformation &tinfo);

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
