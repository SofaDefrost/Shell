#ifndef SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_INL
#define SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_INL

#include "JoinMeshPoints.h"

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
JoinMeshPoints<DataTypes>::JoinMeshPoints()
: f_input_joinPoints(initData(&f_input_joinPoints,"joinPoints","Which points to join"))
  , f_input_position(initData(&f_input_position,"position","Input Vertices"))
  , f_input_normals(initData(&f_input_normals,"normals","Input Normals"))
  , f_input_edges(initData(&f_input_edges,"edges","Input Edges"))
  , f_input_triangles(initData(&f_input_triangles,"triangles","Input Triangles"))
  , f_input_quads(initData(&f_input_quads,"quads","Input Quads"))
  , f_input_tetrahedra(initData(&f_input_tetrahedra,"tetrahedra","Input Tetrahedra"))
  , f_input_hexahedra(initData(&f_input_hexahedra,"hexahedra","Input Hexahedra"))

  , f_output_position(initData(&f_output_position,"joinedPosition","Output Vertices of the joined mesh"))
  , f_output_normals(initData(&f_output_normals,"joinedNormals","Output Normals of the joined mesh"))
  , f_output_edges(initData(&f_output_edges,"joinedEdges","Output Edges of the joined mesh"))
  , f_output_triangles(initData(&f_output_triangles,"joinedTriangles","Output Triangles of the joined mesh"))
  , f_output_quads(initData(&f_output_quads,"joinedQuads","Output Quads of the joined mesh"))
  , f_output_tetrahedra(initData(&f_output_tetrahedra,"joinTetrahedra","Output Tetrahedra of the joined mesh"))
  , f_output_hexahedra(initData(&f_output_hexahedra,"joinHexahedra","Output Hexahedra of the joined mesh"))
  , f_output_mergedPosition(initData(&f_output_mergedPosition,"mergedPosition","Positions of the merged Vertices of the input topology"))
  , f_output_mergedNormals(initData(&f_output_mergedNormals,"mergedNormals","Normals of the merged Vertices of the input topology"))
{
}

template <class DataTypes>
JoinMeshPoints<DataTypes>::~JoinMeshPoints()
{
}

template <class DataTypes>
void JoinMeshPoints<DataTypes>::init()
{
    addInput(&f_input_joinPoints);

    addInput(&f_input_position);
    addInput(&f_input_normals);
    addInput(&f_input_edges);
    addInput(&f_input_triangles);
    addInput(&f_input_quads);
    addInput(&f_input_tetrahedra);
    addInput(&f_input_hexahedra);

    addOutput(&f_output_position);
    addOutput(&f_output_normals);
    addOutput(&f_output_edges);
    addOutput(&f_output_triangles);
    addOutput(&f_output_quads);
    addOutput(&f_output_tetrahedra);
    addOutput(&f_output_hexahedra);

    addOutput(&f_output_mergedPosition);
    addOutput(&f_output_mergedNormals);

    setDirtyValue();
}

template <class DataTypes>
void JoinMeshPoints<DataTypes>::reinit()
{
    update();
}

template <class DataTypes>
void JoinMeshPoints<DataTypes>::doUpdate()
{

    const helper::vector< helper::fixed_array <Index,2> >& inJP = f_input_joinPoints.getValue();
	const VecCoord& inPt = f_input_position.getValue();
	const helper::vector<Vec3>& inNorm = f_input_normals.getValue();

	VecCoord& outPt = *f_output_position.beginEdit();
	helper::vector<Vec3>& outNorm = *f_output_normals.beginEdit();
	VecCoord& outMPt = *f_output_mergedPosition.beginEdit();
	helper::vector<Vec3>& outMNorm = *f_output_mergedNormals.beginEdit();
    outPt.clear();
    outNorm.clear();
    outMPt.resize(inPt.size());
    outMNorm.resize(inPt.size());

    bool bHasNormals = false;
    if (inPt.size() == inNorm.size()) {
        bHasNormals = true;
    } else if (inNorm.size() != 0) {
        serr << "Normal count does not match node count! Ignoring normals." << sendl;
    }

    // Map analogy to the value of f_input_joinPoints. It maps index of a node
    // from input array to index into output array.
    std::map<Index, Index> mapInIn;

#if 0
    for (Index i=0; i<inJP.size(); i++) {
        if (inJP[i][0] == inJP[i][1]) continue;
        Index end = inJP[i][1];
        // Traverse the map to find the last index
        do {
            std::map<Index, Index>::const_iterator iter =
                mapInIn.find(end);
            if (iter != mapInIn.end()) {
                end = iter->second;
            } else {
                break;
            }
        } while (1);
        mapInIn[inJP[i][0]] = end;
    }

    // Order the values so they are allways mapped from larger to smaller.
    // I.e: if i maps to j then i > j.
    for (std::map<Index, Index>::iterator iter = mapInIn.begin();
        iter != mapInIn.end(); iter++) {

        if (iter->first < iter->second) {
            // Invert the mapping
            mapInIn[iter->second] = iter->first;

            std::map<Index, Index>::iterator tmp = iter;
            iter++; // This is safe since iter != mapInIn.end()
            mapInIn.erase(tmp);
            iter--;
        }
    }
#endif

    // Add items so that if i maps to j then i > j
    for (Index i=0; i<inJP.size(); i++) {
        Index from = inJP[i][0];
        Index to = inJP[i][1];

        if (from == to) continue;

        if (from < to) {
            from = inJP[i][1];
            to = inJP[i][0];
        }

        // If there already is 'from' in the map we might need to reassign it
        std::map<Index, Index>::iterator iter =
            mapInIn.find(from);
        if (iter != mapInIn.end()) {
            if (iter->second == to) {
                continue;
            } else if (iter->second < to) {
                // Reassign the current item and keep the old one untact
                from = to;
                to = iter->second;
            } else {
                // Reassign previous item
                Index newFrom = iter->second;
                Index newTo = to;
                mapInIn.erase(iter);
                mapInIn[newFrom] = newTo;
                // ... and add the new item
            }
        }
        mapInIn[from] = to;
    }

    // Traverse the map to find the smallest indices
    for(std::map<Index, Index>::iterator iter = mapInIn.begin();
        iter != mapInIn.end(); iter++) {
        Index end = iter->second;
        do {
            std::map<Index, Index>::const_iterator iter =
                mapInIn.find(end);
            if (iter != mapInIn.end()) {
                end = iter->second;
            } else {
                break;
            }
        } while (1);
        iter->second = end;
    }

    // Create output list
    helper::vector<Index> mapInOut;
    mapInOut.resize(inPt.size());

    for (Index i=0, newId=0; i<inPt.size(); i++) {
        std::map<Index, Index>::const_iterator iter =
            mapInIn.find(i);
        if (iter == mapInIn.end()) {
            outPt.push_back(inPt[i]);
            mapInOut[i] = newId;
            imapNode2Node[newId] = i;
            newId++;
            // Take the position
            outMPt[i] = inPt[i];
            outMNorm[i] = inNorm[i];
        } else {
            mapInOut[i] = mapInOut[iter->second];
            //imapNode2Node[ mapInOut[iter->second] ] = i;
            // Take the position of the point we join with
            outMPt[i] = inPt[iter->second];
            outMNorm[i] = inNorm[iter->second];
        }
    }

    // Compute new normals by averaging the normals of joined nodes.
    if (bHasNormals) {
        outNorm.resize(outPt.size());
        for (Index i=0; i<inPt.size(); i++) {
            // Normalization is necessary for even contribution
            Vec3 n = inNorm[i];
            n.normalize();
            outNorm[mapInOut[i]] += n;
        }

        // We already have the sum, now just divide.
        for (Index i=0; i<outNorm.size(); i++) {
            outNorm[i].normalize();
        }
    }

    f_output_position.endEdit();
    f_output_normals.endEdit();
    f_output_mergedPosition.endEdit();
    f_output_mergedNormals.endEdit();

    // Create output elements with new indices
    createElements<2>(mapInIn, mapInOut, f_input_edges, f_output_edges);
    createElements<3>(mapInIn, mapInOut, f_input_triangles, f_output_triangles);
    createElements<4>(mapInIn, mapInOut, f_input_quads, f_output_quads);
    createElements<4>(mapInIn, mapInOut, f_input_tetrahedra, f_output_tetrahedra);
    createElements<8>(mapInIn, mapInOut, f_input_hexahedra, f_output_hexahedra);
}

template <class DataTypes>
template<unsigned int N>
void JoinMeshPoints<DataTypes>::createElements(
        std::map<Index, Index> mapInIn,
        helper::vector<Index> mapInOut,
        const Data< helper::vector< helper::fixed_array<Index,N> > > &inElements,
        Data< helper::vector< helper::fixed_array<Index,N> > > &outElements)
{
	const helper::vector< helper::fixed_array <Index, N> >& inEle = inElements.getValue();

	helper::vector< helper::fixed_array <Index, N> >& outEle = *outElements.beginEdit();
    outEle.resize(inEle.size());

    for (Index i= 0; i<inEle.size(); i++) {
        helper::fixed_array <Index, N>& out = outEle[i];
        for (Index j=0; j<N; j++) {
            Index id = inEle[i][j];
            std::map<Index, Index>::const_iterator iter =
                mapInIn.find(id);

            bool bJoined = false;
            // In -> In mapping (joining nodes)
            if (iter != mapInIn.end()) {
                id = iter->second;
                bJoined = true;
            }

            // In -> Out mapping (reindexing nodes)
            if (id < f_input_position.getValue().size()) {
                out[j] = mapInOut[id];
            } else {
                // Invalid node ID, what now?
                serr << "Invalid node ID in elements! " <<
                    id << " is not a valid node ID!" << sendl; 
                out[j] = 0;
            }

            if (N == 3) {
                if (bJoined)
                    imapTriNode2Node[i][out[j]] = inEle[i][j];
            }

        }
    }

    outElements.endEdit();
}

} // namespace engine

} // namespace component

} // namespace sofa

#endif // #ifndef SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_INL
