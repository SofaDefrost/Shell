#ifndef SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_H
#define SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_H

#include <sofa/type/Vec.h>
#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/type/Vec.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <SofaShells/config.h>

namespace sofa
{

namespace component
{

namespace engine
{

/**
 * @brief Joins two or more points of the mesh together and updates the
 * topology accordingly.
 *
 * @tparam DataTypes Associated data type.
 */
template <class DataTypes>
class JoinMeshPoints : public sofa::core::DataEngine
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(JoinMeshPoints,DataTypes),core::DataEngine);

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::Coord    Coord;
    typedef typename Coord::value_type   Real;

    typedef type::Vec<3,Real> Vec3;
    typedef unsigned int Index;

protected:
    JoinMeshPoints();

    virtual ~JoinMeshPoints();

public:

    std::string getTemplateName() const override
    {
        return templateName(this);
    }

    static std::string templateName(const JoinMeshPoints<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    void init() override;
    void reinit() override;
    void doUpdate() override;

    // TODO: add methods to find out the inverse mappings

    /**
     * @brief Get index of a node in input topology based on the triangle and
     * node ID in output topology.
     *
     * @param triangleId Triangle from output topology.
     * @param nodeId     Node in output topology.
     *
     * @return Node from input topology.
     */
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

    Data< type::vector< type::fixed_array<Index,2> > > f_input_joinPoints;

    Data<VecCoord> f_input_position;
    Data< type::vector<Vec3> > f_input_normals;
    Data< type::vector< type::fixed_array<Index,2> > > f_input_edges;
    Data< type::vector< type::fixed_array<Index,3> > > f_input_triangles;
    Data< type::vector< type::fixed_array<Index,4> > > f_input_quads;
    Data< type::vector< type::fixed_array<Index,4> > > f_input_tetrahedra;
    Data< type::vector< type::fixed_array<Index,8> > > f_input_hexahedra;

    Data<VecCoord> f_output_position;
    Data< type::vector<Vec3> > f_output_normals;
    Data< type::vector< type::fixed_array<Index,2> > > f_output_edges;
    Data< type::vector< type::fixed_array<Index,3> > > f_output_triangles;
    Data< type::vector< type::fixed_array<Index,4> > > f_output_quads;
    Data< type::vector< type::fixed_array<Index,4> > > f_output_tetrahedra;
    Data< type::vector< type::fixed_array<Index,8> > > f_output_hexahedra;
    Data<VecCoord> f_output_mergedPosition;
    Data< type::vector<Vec3> > f_output_mergedNormals;

private:

    template<unsigned int N> void createElements(
        std::map<Index, Index> mapInIn,
        type::vector<Index> mapInOut,
        const Data< type::vector< type::fixed_array<Index,N> > > &inElements,
        Data< type::vector< type::fixed_array<Index,N> > > &outElements);

    // Inverse mappings: from output to input topology
    // ... node -> node (only for nodes that are not joined and thus have one to one mapping)
    std::map<Index, Index> imapNode2Node;
    // ... triangle + node -> node
    std::map<Index, std::map<Index, Index> > imapTriNode2Node;
};

#if defined(WIN32) && !defined(SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_CPP)
#pragma warning(disable : 4231)
#ifndef SOFA_FLOAT
template class SOFA_SHELLS_API JoinMeshPoints<type::Vec1dTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<type::Vec2dTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<type::Vec3dTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Rigid2dTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Rigid3dTypes>;
#endif //SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_SHELLS_API JoinMeshPoints<type::Vec1fTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<type::Vec2fTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<type::Vec3fTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Rigid2fTypes>;
template class SOFA_SHELLS_API JoinMeshPoints<defaulttype::Rigid3fTypes>;
#endif //SOFA_DOUBLE


#endif

} // namespace engine

} // namespace component

} // namespace sofa
#endif // #ifndef SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_H
