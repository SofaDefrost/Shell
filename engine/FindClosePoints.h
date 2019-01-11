#ifndef SOFA_COMPONENT_ENGINE_FINDCLOSEPOINTS_H
#define SOFA_COMPONENT_ENGINE_FINDCLOSEPOINTS_H

#include <sofa/defaulttype/Vec.h>
#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
//#include <sofa/core/loader/MeshLoader.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
//#include <sofa/component/component.h>

namespace sofa
{

namespace component
{

namespace engine
{

/**
 * @brief Finds points whose distance is smaller than defined threshold.
 *
 * @tparam DataTypes Associated data type.
 */
template <class DataTypes>
class FindClosePoints : public sofa::core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(FindClosePoints,DataTypes),core::DataEngine);

    typedef typename DataTypes::VecCoord    VecCoord;
    typedef typename DataTypes::Coord       Coord;
    typedef typename Coord::value_type      Real;

    typedef unsigned int Index;

protected:
    FindClosePoints();

    virtual ~FindClosePoints();

public:

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const FindClosePoints<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    void init();
    void reinit();
    void doUpdate();

    // Data
    Data<Real>      f_input_threshold;
    Data<VecCoord>  f_input_position;

    Data< helper::vector< helper::fixed_array<Index,2> > > f_output_closePoints;

};

#if defined(WIN32) && !defined(SOFA_COMPONENT_ENGINE_FINDCLOSEPOINTS_CPP)
#pragma warning(disable : 4231)
#ifndef SOFA_FLOAT
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec1Types>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec2Types>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec3Types>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Rigid2Types>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Rigid3Types>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec1Types>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec2Types>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Vec3Types>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Rigid2Types>;
template class SOFA_SHELLS_API FindClosePoints<defaulttype::Rigid3Types>;
#endif //SOFA_DOUBLE


#endif

} // namespace engine

} // namespace component

} // namespace sofa
#endif // #ifndef SOFA_COMPONENT_ENGINE_FINDCLOSEPOINTS_H
