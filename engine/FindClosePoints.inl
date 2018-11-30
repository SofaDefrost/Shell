#ifndef SOFA_COMPONENT_ENGINE_FINDCLOSEPOINTS_INL
#define SOFA_COMPONENT_ENGINE_FINDCLOSEPOINTS_INL

#include "FindClosePoints.h"

namespace sofa
{

namespace component
{

namespace engine
{

//using namespace sofa::helper;
using namespace sofa::defaulttype;
using namespace core::objectmodel;

template <class DataTypes>
FindClosePoints<DataTypes>::FindClosePoints()
: f_input_threshold(initData(&f_input_threshold, (Real)1e-5, "threshold","Threshold"))
, f_input_position(initData(&f_input_position,"position","Vertices"))
, f_output_closePoints(initData(&f_output_closePoints,"closePoints","Pairs of indices of points considered close"))
{
}

template <class DataTypes>
FindClosePoints<DataTypes>::~FindClosePoints()
{
}

template <class DataTypes>
void FindClosePoints<DataTypes>::init()
{
    addInput(&f_input_threshold);
    addInput(&f_input_position);
    addOutput(&f_output_closePoints);

    setDirtyValue();
}

template <class DataTypes>
void FindClosePoints<DataTypes>::reinit()
{
    update();
}

template <class DataTypes>
void FindClosePoints<DataTypes>::doUpdate()
{

    Real threshold = f_input_threshold.getValue();
    const VecCoord& points = f_input_position.getValue();

    helper::vector< helper::fixed_array<Index,2> >& list = *f_output_closePoints.beginEdit();
    list.clear();

    if (points.size() > 1) {
        for (Index i=0; i<points.size()-1; i++) {
            for (Index j=i+1; j<points.size(); j++) {
                if ((points[i] - points[j]).norm() <= threshold) {
                    list.push_back(helper::fixed_array<Index,2>(i,j));
                }
            }
        }
    }


    f_output_closePoints.endEdit();
}

} // namespace engine

} // namespace component

} // namespace sofa

#endif // #ifndef SOFA_COMPONENT_ENGINE_FINDCLOSEPOINTS_INL
