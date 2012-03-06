//
// Interpolation of Bézier triangles
//
// Author: Tomáš Golembiovský
//
// Copyright:
//
//
#ifndef SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_H
#define SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_H

#include "../../initPluginShells.h"
#include <sofa/core/behavior/MechanicalState.h>
//#include <sofa/core/behavior/Mass.h>
//#include <sofa/core/objectmodel/Data.h>
//#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
//#include <sofa/defaulttype/Mat.h>

#include <sofa/core/objectmodel/BaseObject.h>

namespace sofa
{

namespace component
{

namespace fem
{

template<class DataTypes>
class BezierShellInterpolation : public virtual sofa::core::objectmodel::BaseObject
{
    public:

        SOFA_CLASS(SOFA_TEMPLATE(BezierShellInterpolation,DataTypes), sofa::core::objectmodel::BaseObject);

        typedef typename DataTypes::VecCoord VecCoord;
        typedef typename DataTypes::VecDeriv VecDeriv;
        typedef typename DataTypes::VecReal VecReal;
        typedef typename DataTypes::Coord Coord;
        typedef typename DataTypes::Deriv Deriv;
        typedef typename Coord::value_type Real;

        typedef unsigned int Index;
        typedef sofa::core::topology::BaseMeshTopology::TriangleID ElementID;
        //typedef sofa::helper::vector<BaseMeshTopology::EdgeID> VecElementID;
        //typedef sofa::helper::vector<BaseMeshTopology::Edge> VecEdges;
        typedef helper::vector<unsigned int> VecIndex;

        typedef typename  sofa::defaulttype::SolidTypes<Real>::Transform Transform;
        typedef typename  sofa::defaulttype::SolidTypes<Real>::SpatialVector SpatialVector;

        typedef sofa::defaulttype::Vec<2, Real> Vec2;
        typedef sofa::defaulttype::Vec<3, Real> Vec3;

        enum PosType {
            RestPos,
            FreePos,
            CurrentPos
        };


        BezierShellInterpolation()
            : _topology(NULL),
            _mstate(NULL) { }

        ~BezierShellInterpolation(){}

        //void init();
        //void bwdInit();
        //void reinit(){init(); bwdInit(); }
        //void reset(){bwdInit(); }

        // Compute bezier points
        void computeBezierierPoints(VecCoord const x, PosType type);

        void Interpolate(ElementID tri, const Vec2 &baryCoord, PosType type, Vec3 &xOut);
        void Interpolate(ElementID tri, const Vec2 &baryCoord, PosType type, Vec3 &xOut, Vec3 &normal);
        //void Interpolate(ElementID tri, const Vec2 &baryCoord, PosType type, Rigid xIn, Rigid xOut)
        //ComputeBezierierApplyJ(VecDeriv dx, type);
        //ComputeBezierierApplyJT(VecDeriv dx, type);

    protected:
        // pointer to the topology
        sofa::core::topology::BaseMeshTopology* _topology;

        // pointer on mechanical state
        sofa::core::behavior::MechanicalState<DataTypes> *_mstate;


};

} // namespace fem

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_H
