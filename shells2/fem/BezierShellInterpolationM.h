//
// Interpolation of Bézier triangles (with mechanical mapping)
//
// Author: Tomáš Golembiovský
//
// Copyright:
//
#ifndef SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATIONM_H
#define SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATIONM_H

#include "../../initPluginShells.h"
#include "BezierShellInterpolation.h"

#include <sofa/core/behavior/MechanicalState.h>
#include <MechanicalObject.h>
#include <Mesh2PointTopologicalMapping.h>
//#include <sofa/core/behavior/Mass.h>
//#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <TopologyData.h>
#include <sofa/simulation/common/AnimateBeginEvent.h>

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

template <class TIn, class TOut>
class BezierShellInterpolationM : public BezierShellInterpolation<TIn>
{
    public:

        SOFA_CLASS(SOFA_TEMPLATE2(BezierShellInterpolationM,TIn,TOut), SOFA_TEMPLATE(BezierShellInterpolation,TIn));

        typedef TIn In;
        typedef TOut Out;

        typedef typename In::VecCoord               InVecCoord;
        typedef typename In::VecDeriv               InVecDeriv;
        typedef typename In::Coord                  InCoord;
        typedef typename In::Deriv                  InDeriv;
        typedef typename In::MatrixDeriv            InMatrixDeriv;
        typedef typename In::Real                   InReal;

        typedef typename Out::VecCoord              OutVecCoord;
        typedef typename Out::VecDeriv              OutVecDeriv;
        typedef typename Out::Coord                 OutCoord;
        typedef typename Out::Deriv                 OutDeriv;
        typedef typename Out::MatrixDeriv           OutMatrixDeriv;
        typedef typename Out::Real                  OutReal;

        typedef InReal Real;

        typedef typename Inherit1::Index Index;
        typedef typename Inherit1::VecIndex VecIndex;

        typedef typename Inherit1::Vec2 Vec2;
        typedef typename Inherit1::Vec3 Vec3;
        typedef typename Inherit1::Mat33 Mat33;
        typedef typename Inherit1::VecVec3 VecVec3;

        typedef typename Inherit1::VecVec3d VecVec3d;

        typedef typename Inherit1::ShapeFunctions ShapeFunctions;
        typedef typename Inherit1::BTri BTri;
        typedef typename Inherit1::VecShapeFunctions VecShapeFunctions;
        typedef typename Inherit1::VecBTri VecBTri;

        BezierShellInterpolationM() {}
        virtual ~BezierShellInterpolationM() {}

        virtual std::string getTemplateName() const
        {
            return templateName(this);
        }

        static std::string templateName(const BezierShellInterpolationM<TIn, TOut>* = NULL)
        {
            return TIn::Name() + std::string(",") + TOut::Name();
        }

        void applyOnBTriangle(VecShapeFunctions projShapeFunctions, VecIndex projElements, helper::WriteAccessor< Data<VecVec3> > &out);
        void applyJOnBTriangle(VecShapeFunctions projShapeFunctions, VecIndex projElements, const InVecDeriv& in, helper::WriteAccessor< Data<OutVecDeriv> > &out);
        void applyJTOnBTriangle(VecShapeFunctions projShapeFunctions, VecIndex projElements, const OutVecDeriv& in, helper::WriteAccessor< Data<InVecDeriv> > &out);
        void applyJTOnBTriangle(VecShapeFunctions projN, VecIndex projElements,
            const OutMatrixDeriv& in, InMatrixDeriv &out);

    protected:

        void applyJTCore(const InVecCoord &xSim, const VecVec3d &x,
            const Index &triId, const ShapeFunctions &N, const OutCoord &force,
            Vec3 &f1, Vec3 &f2, Vec3 &f3, Vec3 &f1r, Vec3 &f2r, Vec3 &f3r);


};

} // namespace fem

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATIONM_H
