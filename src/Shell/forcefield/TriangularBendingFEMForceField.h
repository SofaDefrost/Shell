/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/Data.h>

#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/topology/TopologyData.h>

#include <Shell/controller/MeshInterpolator.h>
#include <Shell/engine/JoinMeshPoints.h>


namespace shell::forcefield
{

using namespace sofa::type;
using sofa::type::vector;
using namespace sofa::core::topology;
using namespace sofa::core::behavior;

/// This class can be overridden if needed for additionnal storage within template specializations.
template<class DataTypes>
class TriangularBendingFEMForceFieldInternalData
{
public:
};

template<class DataTypes>
class TriangularBendingFEMForceField : public sofa::core::behavior::ForceField<DataTypes>
{
    public:
        SOFA_CLASS(SOFA_TEMPLATE(TriangularBendingFEMForceField,DataTypes), SOFA_TEMPLATE(sofa::core::behavior::ForceField,DataTypes));

        typedef sofa::core::behavior::ForceField<DataTypes>       Inherited;
        typedef typename DataTypes::VecCoord                VecCoord;
        typedef typename DataTypes::VecDeriv                VecDeriv;

        typedef typename DataTypes::Coord                   Coord;
        typedef typename DataTypes::Deriv                   Deriv;
        typedef typename Coord::value_type                  Real;

        typedef Vec<3,Real> Vec3;
        typedef Vec<2,Real> Vec2;
        typedef sofa::type::Quat<Real> Quat;

        typedef sofa::Data<VecCoord>                              DataVecCoord;
        typedef sofa::Data<VecDeriv>                              DataVecDeriv;

        typedef sofa::defaulttype::Vec3Types::VecCoord VecCoordHigh;

        typedef sofa::Index Index;
        typedef sofa::core::topology::BaseMeshTopology::Triangle Triangle;
        typedef sofa::core::topology::BaseMeshTopology::SeqTriangles SeqTriangles;

    protected:

        typedef Vec<6, Real> Displacement;                      ///< the displacement vector
        typedef Vec<9, Real> DisplacementBending;               ///< the displacement vector for bending
        typedef Mat<3, 3, Real> MaterialStiffness;              ///< the matrix of material stiffness
        typedef Mat<6, 3, Real> StrainDisplacement;             ///< the strain-displacement matrix
        typedef Mat<3, 9, Real> StrainDisplacementBending;
        typedef Mat<3, 3, Real > Transformation;                ///< matrix for rigid transformations like rotations
        typedef Mat<6, 6, Real> StiffnessMatrix;
        typedef Mat<9, 9, Real> StiffnessMatrixBending;
        typedef Mat<18, 18, Real> StiffnessMatrixGlobalSpace;

        sofa::core::topology::BaseMeshTopology* _topology;

public:

        class TriangleInformation
        {
            public:

                // position of B, C vertices in local (in-plane) coordinates (A is naturaly [0,0])
                sofa::type::fixed_array <Vec3, 2> restLocalPositions;
                sofa::type::fixed_array <Quat, 3> restLocalOrientations;

                /// material stiffness matrices of each tetrahedron
                MaterialStiffness materialMatrix;
                // the strain-displacement matrices vector
                StrainDisplacement strainDisplacementMatrix;
                // the strain-displacement matrices vector (at Gauss points)
                StrainDisplacementBending strainDisplacementMatrix1;
                StrainDisplacementBending strainDisplacementMatrix2;
                StrainDisplacementBending strainDisplacementMatrix3;
                // Indices of each vertex
                Index a, b, c;
                // Indices in rest shape topology. Normaly is the same as a, b,
                // c, but is different if topologyMapper is set.
                Index a0, b0, c0;
                // Local coordinates
                Vec3 localB, localC;
                // Transformation rotation;
                Quat Qframe;
                // Stiffness matrix K = J * M * Jt
                StiffnessMatrix stiffnessMatrix;
                // Stiffness matrix for bending K = Jt * M * J
                StiffnessMatrixBending stiffnessMatrixBending;

                // Surface
                Real area;
                // Variables needed for drawing the shell
                Vec<9, Real> u; // displacement vector
                Mat<9, 9, Real> invC; // inverse of C (used in bending mode only)
                Vec<9, Real> coefficients; // coefficients Ci computed from tinfo->invC * (tinfo->u + tinfo->u_flat)
                Vec <9, Real> u_rest; // difference between the initial position and the flate position to allow the use of an initial deformed shape

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
                TRQSTriangleHandler(TriangularBendingFEMForceField<DataTypes>* _ff, TriangleData<sofa::type::vector<TriangleInformation> >* _data) : TopologyDataHandler<Triangle, sofa::type::vector<TriangleInformation> >(_data), ff(_ff) {}

                void applyCreateFunction(unsigned int triangleIndex, TriangleInformation& ,
                    const Triangle & t,
                    const sofa::type::vector< unsigned int > &,
                    const sofa::type::vector< double > &);

            protected:
                TriangularBendingFEMForceField<DataTypes>* ff;
        };

        TriangularBendingFEMForceField();

        virtual ~TriangularBendingFEMForceField();
        void init() override;
        void reinit() override;
        void addForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& dataF, const DataVecCoord& dataX, const DataVecDeriv& /*dataV*/ ) override ;
        void addDForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& datadF, const DataVecDeriv& datadX ) override ;
        void addKToMatrix(const sofa::core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix) override;
        void addBToMatrix(sofa::linearalgebra::BaseMatrix * /*mat*/, double /*bFact*/, unsigned int &/*offset*/) override;
        double getPotentialEnergy(const VecCoord& x) const;
        void handleTopologyChange() override;

        SReal getPotentialEnergy(const sofa::core::MechanicalParams* /*mparams*/, const DataVecCoord& /*x*/) const override { return 0; }

        void draw(const sofa::core::visual::VisualParams* vparams) override;

        sofa::core::topology::BaseMeshTopology* getTopology() {return _topology;}
        TriangleData< sofa::type::vector<TriangleInformation> >& getTriangleInfo() {return triangleInfo;}

        sofa::Data<Real> d_poisson;
        sofa::Data<Real> d_young;
        sofa::Data<bool> d_bending;
        sofa::Data <Real> d_thickness;
        sofa::Data <Real> d_membraneRatio;
        sofa::Data <Real> d_bendingRatio;
        sofa::Data<bool> d_refineMesh;
        sofa::Data<int> d_iterations;
        sofa::SingleLink<TriangularBendingFEMForceField<DataTypes>,
            sofa::core::topology::BaseMeshTopology,
            sofa::BaseLink::FLAG_STOREPATH|sofa::BaseLink::FLAG_STRONGLINK> l_targetTopology;
        VecCoordHigh m_targetVertices;
        SeqTriangles m_targetTriangles;

        // Allow transition between rest shapes
        sofa::SingleLink<TriangularBendingFEMForceField<DataTypes>,
            shell::controller::MeshInterpolator<DataTypes>,
            sofa::BaseLink::FLAG_STOREPATH|sofa::BaseLink::FLAG_STRONGLINK> l_restShape;

        // Indirect rest shape indexing (e.g. for "joining" two meshes)
        bool m_mapTopology;
        sofa::SingleLink<TriangularBendingFEMForceField<DataTypes>,
            shell::engine::JoinMeshPoints<DataTypes>,
            sofa::BaseLink::FLAG_STOREPATH|sofa::BaseLink::FLAG_STRONGLINK> l_topologyMapper;


        sofa::core::objectmodel::DataFileName m_exportFilename;
        sofa::Data<unsigned int> d_exportEveryNbSteps;
        sofa::Data<bool> d_exportAtBegin;
        sofa::Data<bool> d_exportAtEnd;
        unsigned int m_stepCounter;

        TRQSTriangleHandler* m_triangleHandler;

protected :

        TriangleData< sofa::type::vector<TriangleInformation> > triangleInfo;

        void computeDisplacement(Displacement &Disp, const VecCoord &x, const Index elementIndex);
        void computeDisplacementBending(DisplacementBending &Disp, const VecCoord &x, const Index elementIndex);
        void computeStrainDisplacementMatrix(StrainDisplacement &J, const Index elementIndex, const Vec3& b, const Vec3& c);
        void computeStrainDisplacementMatrixBending(TriangleInformation *tinfo, const Vec3& b, const Vec3& c);
        void tensorFlatPlate(Mat<3, 9, Real>& D, const Vec3 &P);
        void computeStiffnessMatrix(StiffnessMatrix &K, const StrainDisplacement &J, const MaterialStiffness &M);
        void computeStiffnessMatrixBending(StiffnessMatrixBending &K, TriangleInformation *tinfo);
        void computeForce(Displacement &F, const Displacement& D, const Index elementIndex);
        void computeForceBending(DisplacementBending &F, const DisplacementBending& D, const Index elementIndex);

        static void TRQSTriangleCreationFunction (unsigned int , void* , TriangleInformation &, const Triangle& , const sofa::type::vector< unsigned int > &, const sofa::type::vector< double >&);

        /// f += Kx where K is the stiffness matrix and x a displacement
        virtual void applyStiffness(VecDeriv& f, const VecDeriv& dx, const Index elementIndex, const double kFactor);
        virtual void computeMaterialStiffness(const int i);

        void initTriangleOnce(const int i, const Index&a, const Index&b, const Index&c);
        void initTriangle(const int i);
        void computeRotation(Quat &Qframe, const VecCoord &p, const Index &a, const Index &b, const Index &c);
        void accumulateForce(VecDeriv& f, const VecCoord & p, const Index elementIndex);

        void convertStiffnessMatrixToGlobalSpace(StiffnessMatrixGlobalSpace &K_gs, TriangleInformation *tinfo);

        void refineCoarseMeshToTarget(void);
        void subdivide(const Vec3& a, const Vec3& b, const Vec3& c, sofa::type::vector<Vec3> &subVertices, SeqTriangles &subTriangles);
        void addVertexAndFindIndex(sofa::type::vector<Vec3> &subVertices, const Vec3 &vertex, int &index);
        void movePoint(Vec3& pointToMove);
        void findClosestGravityPoints(const Vec3& point, sofa::type::vector<Vec3>& listClosestPoints);

        void handleEvent(sofa::core::objectmodel::Event *event) override;


        Quat qDiff(Quat a, const Quat& b)
        {
            // If the axes are not oriented in the same direction, flip the axis and angle of a to get the same convention than b
            if (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]<0)
            {
                a[0] = -a[0];
                a[1] = -a[1];
                a[2] = -a[2];
                a[3] = -a[3];
            }
            Quat q = b.inverse() * a;
            return q;
        }

        Quat qDiffZ(const Quat& vertex, const Quat& Qframe)
        {
            // dQ is the quaternion that embodies the rotation between the z axis of the vertex and the z axis of the local triangle's frame (in local space)
            Quat dQ;

            // u = z axis of the triangle's frame
            Vec3d u(0,0,1);

            // v = z axis of the vertex's frame is expressed into world space
            Vec3d v = vertex.rotate(Vec3d(0.0, 0.0, 1.0));
            // v0 = v expressed into local triangle's frame
            Vec3d v0 = Qframe.rotate(v);

            // Axis of rotation between the 2 vectors u and v lies into the plan of the 2 vectors
            Vec3d axis = cross(u, v0);
            // Shortest angle between the 2 vectors
            double angle = acos(dot(u, v0));

            // Quaternion associated to this axis and this angle
            if (fabs(angle)>1e-6)
            {
                dQ.axisToQuat(axis,angle);
            }
            else
            {
                dQ = Quat(0,0,0,1);
            }

            return dQ;
        }
};

} // namespace
