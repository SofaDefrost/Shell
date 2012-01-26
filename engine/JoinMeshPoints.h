#ifndef SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_H
#define SOFA_COMPONENT_ENGINE_JOINMESHSPOINTS_H

#include <sofa/defaulttype/Vec.h>
#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/loader/MeshLoader.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/component/component.h>

namespace sofa
{

namespace component
{

namespace engine
{

template <class DataTypes>
class JoinMeshPoints : public sofa::core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(JoinMeshPoints,DataTypes),core::DataEngine);

    typedef typename DataTypes::VecCoord VecCoord;
    typedef helper::vector<unsigned int> VecIndex;

protected:
    JoinMeshPoints();

    ~JoinMeshPoints();

public:

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const JoinMeshPoints<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    void init();
    void reinit();
    void update();

    // Data

    Data< helper::vector< helper::fixed_array <unsigned int,2> > > f_input_joinPoints;

    Data<VecCoord> f_input_positions;
    Data< helper::vector< helper::fixed_array <unsigned int,2> > > f_input_edges;
    Data< helper::vector< helper::fixed_array <unsigned int,3> > > f_input_triangles;
    Data< helper::vector< helper::fixed_array <unsigned int,4> > > f_input_quads;
    Data< helper::vector< helper::fixed_array<unsigned int,4> > > f_input_tetrahedra;
    Data< helper::vector< helper::fixed_array<unsigned int,8> > > f_input_hexahedra;

    Data<VecCoord> f_output_positions;
    Data< helper::vector< helper::fixed_array <unsigned int,2> > > f_output_edges;
    Data< helper::vector< helper::fixed_array <unsigned int,3> > > f_output_triangles;
    Data< helper::vector< helper::fixed_array <unsigned int,4> > > f_output_quads;
    Data< helper::vector< helper::fixed_array<unsigned int,4> > > f_output_tetrahedra;
    Data< helper::vector< helper::fixed_array<unsigned int,8> > > f_output_hexahedra;

private:

    template<unsigned int N> void createElements(
        std::map<unsigned int, unsigned int> mapInIn,
        helper::vector<unsigned int> mapInOut,
        const Data< helper::vector< helper::fixed_array<unsigned int,N> > > &inElements,
        Data< helper::vector< helper::fixed_array<unsigned int,N> > > &outElements);

};

#if defined(WIN32) && !defined(SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_CPP)
#pragma warning(disable : 4231)
#ifndef SOFA_FLOAT
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec1dTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec2dTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec3dTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Rigid2dTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Rigid3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec1fTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec2fTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Vec3fTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Rigid2fTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Rigid3fTypes>;
#endif //SOFA_DOUBLE


#endif

} // namespace engine

} // namespace component

} // namespace sofa
#endif // #ifndef SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_H
