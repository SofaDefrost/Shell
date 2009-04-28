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
#ifndef SOFA_COMPONENT_FORCEFIELD_TRIANGULARBENDINGFEMFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_TRIANGULARBENDINGFEMFORCEFIELD_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/core/componentmodel/behavior/ForceField.h>
#include <sofa/core/componentmodel/topology/BaseMeshTopology.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/component/topology/TriangleData.h>
#include <sofa/component/topology/EdgeData.h>
#include <sofa/component/topology/PointData.h>
#include <newmat/newmat.h>
#include <newmat/newmatap.h>



namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;
using sofa::helper::vector;
using namespace sofa::component::topology;


template<class DataTypes>
class TriangularBendingFEMForceField : public core::componentmodel::behavior::ForceField<DataTypes>, public virtual core::objectmodel::BaseObject
{
public:
  typedef core::componentmodel::behavior::ForceField<DataTypes> Inherited;
	typedef typename DataTypes::VecCoord VecCoord;
	typedef typename DataTypes::VecDeriv VecDeriv;
	typedef typename DataTypes::VecReal VecReal;
	typedef VecCoord Vector;
	typedef typename DataTypes::Coord    Coord   ;
	typedef typename DataTypes::Deriv    Deriv   ;
	typedef typename Coord::value_type   Real    ;
	typedef Vec<3,Real> Vec3;

	typedef sofa::core::componentmodel::topology::BaseMeshTopology::index_type Index;
	typedef sofa::core::componentmodel::topology::BaseMeshTopology::Triangle Element;
	typedef sofa::core::componentmodel::topology::BaseMeshTopology::SeqTriangles VecElement;

protected:

	typedef Vec<15, Real> Displacement;                                      ///< the displacement vector
	typedef Mat<3, 3, Real> MaterialStiffness;				///< the matrix of material stiffness
	typedef sofa::helper::vector<MaterialStiffness> VecMaterialStiffness;   ///< a vector of material stiffness matrices
	typedef Mat<6, 3, Real> StrainDisplacement;				///< the strain-displacement matrix
	typedef sofa::helper::vector<StrainDisplacement> VecStrainDisplacement;	///< a vector of strain-displacement matrices
	typedef Mat<3, 3, Real > Transformation;				///< matrix for rigid transformations like rotations
        typedef helper::fixed_array <Vec3, 3> RenderingTriangle;                ///> contains the 3 summets of a triangle
        typedef sofa::helper::vector<RenderingTriangle> ListTriangles;          ///> vector of triangles

	class TriangleInformation
        {
            public:

                helper::fixed_array <Vec3, 2> restLocalPositions;
                helper::fixed_array <Quat, 3> initialOrientations;

                /// material stiffness matrices of each tetrahedron
                MaterialStiffness materialMatrix;
                // the strain-displacement matrices vector
                StrainDisplacement strainDisplacementMatrix;
                // bending strain-displacement matrices at each Gauss point
                Mat<3, 9, Real> b1;
                Mat<3, 9, Real> b2;
                Mat<3, 9, Real> b3;
                // Transformation rotation;
                Quat Qframe0;
                Quat Qframe;
                // strain vector
                Vec3 strain;
                // strain caused by bending at each Gauss point
                Vec3 bendingStrain1;
                Vec3 bendingStrain2;
                Vec3 bendingStrain3;
                // stress vector
                Vec3 stress;
                // stress caused by bending at each Gauss point
                Vec3 bendingStress1;
                Vec3 bendingStress2;
                Vec3 bendingStress3;
                
                Real thirdSurface;

                TriangleInformation() { }

                // variables needed for drawing the shell
                Vec<9, Real> u; // displacement vector
                Mat<9, 9, Real> invC; // inverse of C (used in bending mode only)
                
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

	TriangleData<TriangleInformation> triangleInfo;

	sofa::core::componentmodel::topology::BaseMeshTopology* _topology;

public:

    TriangularBendingFEMForceField();

	virtual ~TriangularBendingFEMForceField();
	virtual void init();
	virtual void reinit();
	virtual void addForce (VecDeriv& f, const VecCoord& x, const VecDeriv& v);
	virtual void addDForce (VecDeriv& df, const VecDeriv& dx);
	virtual double getPotentialEnergy(const VecCoord& x);
	virtual void handleTopologyChange();
	virtual void draw();

	Data<Real> f_poisson;
	Data<Real> f_young;
        Data<bool> f_bending;
        Data <Real> f_thickness;
        Data<int> subdivisions;


protected :

	void computeDisplacement(Displacement &Disp, Index elementIndex, const VecCoord &p);
        void computeDisplacementBending(Displacement &Disp, Index elementIndex, const VecCoord &p);
	void computeStrainDisplacement( StrainDisplacement &J, Vec3 a, Vec3 b, Vec3 c );
        void computeStrainDisplacementBending(const Index elementIndex, Vec3& /*a*/, Vec3& b, Vec3& c );
        void tensorFlatPlate(Mat<3, 9, Real>& D, Vec3 &P);
	void computeStrain(Vec<3,Real> &strain, const StrainDisplacement &J, const Displacement &D);
        void computeStrainBending(const Index& elementIndex, const Displacement &D);
	void computeStress(Vec<3,Real> &stress, MaterialStiffness &K, Vec<3,Real> &strain);
        void computeStressBending(const Index& elementIndex);
	void computeForce(Displacement &F, Index elementIndex, const VecCoord &p);

	static void TRQSTriangleCreationFunction (int , void* , TriangleInformation &, const Triangle& , const sofa::helper::vector< unsigned int > &, const sofa::helper::vector< double >&);

	/// f += Kx where K is the stiffness matrix and x a displacement
	virtual void applyStiffness( VecDeriv& f, Real h, const VecDeriv& dx );
	virtual void computeMaterialStiffness(int i, Index& a, Index& b, Index& c);

	////////////// large displacements method
	//sofa::helper::vector< helper::fixed_array <Coord, 3> > _rotatedInitialElements;   ///< The initials positions in its frame
	//sofa::helper::vector< Transformation > _rotations;
	void initLarge(int i, Index&a, Index&b, Index&c);
	void computeRotation(Quat &Qframe, const VecCoord &p, const Index &a, const Index &b, const Index &c);
	void accumulateForce(VecDeriv& f, const VecCoord & p, Index elementIndex);

        void subdivide(const ListTriangles listTriangles, ListTriangles& newListTriangles);
        void computeDeflection(ListTriangles &listTriangles, const Vec3 &a0, const Quat &Qframe, const Mat<9, 9, Real> &invC, const Vec <9, Real> &u);
        void renderTriangles(const ListTriangles& listTriangles);
};


#if defined(WIN32) && !defined(SOFA_COMPONENT_FORCEFIELD_TRIANGULARBENDINGFEMFORCEFIELD_CPP)
#pragma warning(disable : 4231)
#ifndef SOFA_FLOAT
extern template class SOFA_COMPONENT_FORCEFIELD_API TriangularBendingFEMForceField<defaulttype::Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_COMPONENT_FORCEFIELD_API TriangularBendingFEMForceField<defaulttype::Rigid3fTypes>;
#endif
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif
