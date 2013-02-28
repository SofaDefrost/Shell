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
#ifndef SOFA_COMPONENT_FORCEFIELD_CST_FEM_FORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_CST_FEM_FORCEFIELD_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/Data.h>

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/component/topology/TopologyData.h>


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
class CstFEMForceFieldInternalData
{
public:
};


template<class DataTypes>
class CstFEMForceField : public core::behavior::ForceField<DataTypes>
{
    public:
        SOFA_CLASS(SOFA_TEMPLATE(CstFEMForceField,DataTypes), SOFA_TEMPLATE(core::behavior::ForceField,DataTypes));

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

        typedef helper::Quater<Real> Quat;

        typedef Data<VecCoord>                              DataVecCoord;
        typedef Data<VecDeriv>                              DataVecDeriv;

        typedef Vec3Types::VecCoord VecCoordHigh;

        typedef sofa::core::topology::BaseMeshTopology::index_type Index;
        typedef sofa::core::topology::BaseMeshTopology::Triangle Triangle;
        typedef sofa::core::topology::BaseMeshTopology::SeqTriangles SeqTriangles;

        class TriangleInformation;

    protected:

        typedef Vec<6, Real> Displacement;                      // the displacement vector
        typedef Mat<3, 3, Real> MaterialStiffness;              // the matrix of material stiffness
        typedef Mat<3, 6, Real> StrainDisplacement;             // the strain-displacement matrix
        typedef Mat<6, 3, Real> StrainDisplacementT;            // the strain-displacement matrix transposed
        typedef Mat<3, 3, Real> Transformation;                 // matrix for rigid transformations like rotations
        typedef Mat<6, 6, Real> StiffnessMatrix;                // element stiffness matrix
        typedef Mat<9, 9, Real> StiffnessMatrixFull;            // element stiffness matrix in 3D

        sofa::core::topology::BaseMeshTopology* _topology;

public:

        class TriangleInformation
        {
            public:

                // Indices of each vertex
                Index a, b, c;

                // Rest position in local (in-plane) coordinates
                helper::fixed_array <Vec3, 3> restPositions;

                // Deformed position in local (in-plane) coordinates
                helper::fixed_array <Vec3, 3> deformedPositions;

                // Frame rotation as matrix and quaternion
                Transformation R, Rt;

                // Stiffness matrix
                StiffnessMatrix stiffnessMatrix;

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

        class TRQSTriangleHandler : public TopologyDataHandler<Triangle, sofa::helper::vector<TriangleInformation> >
        {
            public:
                TRQSTriangleHandler(CstFEMForceField<DataTypes>* _ff, TriangleData<sofa::helper::vector<TriangleInformation> >* _data) : TopologyDataHandler<Triangle, sofa::helper::vector<TriangleInformation> >(_data), ff(_ff) {}

                void applyCreateFunction(unsigned int triangleIndex, TriangleInformation& ,
                    const Triangle & t,
                    const sofa::helper::vector< unsigned int > &,
                    const sofa::helper::vector< double > &);

            protected:
                CstFEMForceField<DataTypes>* ff;
        };

        CstFEMForceField();

        virtual ~CstFEMForceField();
        virtual void init();
        virtual void reinit();
        virtual void addForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& dataF, const DataVecCoord& dataX, const DataVecDeriv& /*dataV*/ ) ;
        virtual void addDForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& datadF, const DataVecDeriv& datadX ) ;
        virtual void addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix);

        sofa::core::topology::BaseMeshTopology* getTopology() {return _topology;}

        Data<Real> f_poisson;
        Data<Real> f_young;
        Data<Real> f_thickness;
        Data<bool> f_corotated;
        Data<int> f_startDTAppl;
        Data<int> f_numDTAppl;

        Real applyFactor;
        int actualStep;
        //Data<sofa::helper::OptionsGroup> f_measure;
        //Data< helper::vector<Real> > f_measuredValues;

        TRQSTriangleHandler* triangleHandler;

protected :

        /// Material stiffness matrix
        MaterialStiffness materialMatrix, materialMatrixMembrane;
        TriangleData< sofa::helper::vector<TriangleInformation> > triangleInfo;

        // What to measure
        //bool bMeasureStrain;
        //bool bMeasureStress;

        void initTriangle(const int i, const Index&a, const Index&b, const Index&c);

        //void computeRotation(Transformation& R, const VecCoord &x, const Index &a, const Index &b, const Index &c);
        void computeRotation(Transformation& R, const helper::fixed_array<Vec3, 3> &x);
        void computeMaterialStiffness();
        void computeStiffnessMatrix(StiffnessMatrix &K, TriangleInformation &tinfo);

        void computeDisplacement(Displacement &D, const VecCoord &x, const Index elementIndex);
        void accumulateForce(VecDeriv& f, const VecCoord & p, const Index elementIndex);
        void computeForce(Displacement &F, const Displacement& D, const Index elementIndex);
        virtual void applyStiffness(VecDeriv& f, const VecDeriv& dx, const Index elementIndex, const double kFactor);

        void convertStiffnessMatrixToGlobalSpace(StiffnessMatrixFull &K_gs, const TriangleInformation &tinfo);

        void handleEvent(core::objectmodel::Event *event);
};


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif
