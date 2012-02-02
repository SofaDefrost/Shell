#ifndef SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_H
#define SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_H

#include <sofa/defaulttype/Vec.h>
#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/component/component.h>
#include "../initPluginShells.h"

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
    typedef typename DataTypes::Coord    Coord;
    typedef typename Coord::value_type   Real;

    typedef defaulttype::Vec<3,Real> Vec3;
    typedef unsigned int Index;

protected:
    JoinMeshPoints();

    virtual ~JoinMeshPoints();

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

    // TODO: add methods to find out the inverse mappings

    // Return index of a node in input topology based on the triangle and node
    // ID in output topology.
    Index getSrcNodeFromTri(Index triangleId, Index nodeId) {
        if (triangleId < f_output_triangles.getValue().size()) {
            if (nodeId < f_output_position.getValue().size()) {
                const std::map<Index, Index> &n2n = imapTriNode2Node[triangleId];
                std::map<Index, Index>::const_iterator iter = n2n.find(nodeId);
                if (iter != n2n.end()) {
                    return iter->second;
                } else {
                    // Not a joined point
                    return imapNode2Node[nodeId];
                }
            }
        }
        // Invalid indices
        serr << "Requested invalid triangle " << triangleId << " and node " << nodeId << sendl;
        return 0;
    }


    // Data

    Data< helper::vector< helper::fixed_array<Index,2> > > f_input_joinPoints;

    Data<VecCoord> f_input_position;
    Data< helper::vector<Vec3> > f_input_normals;
    Data< helper::vector< helper::fixed_array<Index,2> > > f_input_edges;
    Data< helper::vector< helper::fixed_array<Index,3> > > f_input_triangles;
    Data< helper::vector< helper::fixed_array<Index,4> > > f_input_quads;
    Data< helper::vector< helper::fixed_array<Index,4> > > f_input_tetrahedra;
    Data< helper::vector< helper::fixed_array<Index,8> > > f_input_hexahedra;

    Data<VecCoord> f_output_position;
    Data< helper::vector<Vec3> > f_output_normals;
    Data< helper::vector< helper::fixed_array<Index,2> > > f_output_edges;
    Data< helper::vector< helper::fixed_array<Index,3> > > f_output_triangles;
    Data< helper::vector< helper::fixed_array<Index,4> > > f_output_quads;
    Data< helper::vector< helper::fixed_array<Index,4> > > f_output_tetrahedra;
    Data< helper::vector< helper::fixed_array<Index,8> > > f_output_hexahedra;
    Data<VecCoord> f_output_mergedPosition;
    Data< helper::vector<Vec3> > f_output_mergedNormals;

private:

    template<unsigned int N> void createElements(
        std::map<Index, Index> mapInIn,
        helper::vector<Index> mapInOut,
        const Data< helper::vector< helper::fixed_array<Index,N> > > &inElements,
        Data< helper::vector< helper::fixed_array<Index,N> > > &outElements);

    // Inverse mappings: from output to input topology
    // ... node -> node (only for nodes that are not joined and thus have one to one mapping)
    std::map<Index, Index> imapNode2Node;
    // ... triangle + node -> node
    std::map<Index, std::map<Index, Index> > imapTriNode2Node;
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
