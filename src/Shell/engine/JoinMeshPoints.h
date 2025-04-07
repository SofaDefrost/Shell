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

#include <sofa/type/Vec.h>
#include <sofa/core/DataEngine.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/type/Vec.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <Shell/config.h>

namespace shell::engine
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
    SOFA_CLASS(SOFA_TEMPLATE(JoinMeshPoints,DataTypes), sofa::core::DataEngine);

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::Coord    Coord;
    typedef typename Coord::value_type   Real;

    typedef sofa::type::Vec<3,Real> Vec3;
    typedef unsigned int Index;

protected:
    JoinMeshPoints();

    virtual ~JoinMeshPoints();

public:

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
        msg_error() << "Requested invalid triangle " << triangleId << " and node " << nodeId;
        return 0;
    }


    // Data

    sofa::Data< sofa::type::vector< sofa::type::fixed_array<Index,2> > > f_input_joinPoints;

    sofa::Data<VecCoord> f_input_position;
    sofa::Data< sofa::type::vector<Vec3> > f_input_normals;
    sofa::Data< sofa::type::vector< sofa::type::fixed_array<Index,2> > > f_input_edges;
    sofa::Data< sofa::type::vector< sofa::type::fixed_array<Index,3> > > f_input_triangles;
    sofa::Data< sofa::type::vector< sofa::type::fixed_array<Index,4> > > f_input_quads;
    sofa::Data< sofa::type::vector< sofa::type::fixed_array<Index,4> > > f_input_tetrahedra;
    sofa::Data< sofa::type::vector< sofa::type::fixed_array<Index,8> > > f_input_hexahedra;

    sofa::Data<VecCoord> f_output_position;
    sofa::Data< sofa::type::vector<Vec3> > f_output_normals;
    sofa::Data< sofa::type::vector< sofa::type::fixed_array<Index,2> > > f_output_edges;
    sofa::Data< sofa::type::vector< sofa::type::fixed_array<Index,3> > > f_output_triangles;
    sofa::Data< sofa::type::vector< sofa::type::fixed_array<Index,4> > > f_output_quads;
    sofa::Data< sofa::type::vector< sofa::type::fixed_array<Index,4> > > f_output_tetrahedra;
    sofa::Data< sofa::type::vector< sofa::type::fixed_array<Index,8> > > f_output_hexahedra;
    sofa::Data<VecCoord> f_output_mergedPosition;
    sofa::Data< sofa::type::vector<Vec3> > f_output_mergedNormals;

private:

    template<unsigned int N> void createElements(
        std::map<Index, Index> mapInIn,
        sofa::type::vector<Index> mapInOut,
        const sofa::Data< sofa::type::vector< sofa::type::fixed_array<Index,N> > > &inElements,
        sofa::Data< sofa::type::vector< sofa::type::fixed_array<Index,N> > > &outElements);

    // Inverse mappings: from output to input topology
    // ... node -> node (only for nodes that are not joined and thus have one to one mapping)
    std::map<Index, Index> imapNode2Node;
    // ... triangle + node -> node
    std::map<Index, std::map<Index, Index> > imapTriNode2Node;
};

} // namespace
