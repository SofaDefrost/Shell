#ifndef SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_INL
#define SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_INL

#include <engine/JoinMeshPoints.h>
//#include <sofa/core/visual/VisualParams.h>
//#include <sofa/helper/gl/template.h>

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
  , f_input_positions(initData(&f_input_positions,"position","Input Vertices"))
  , f_input_edges(initData(&f_input_edges,"edges","Input Edges"))
  , f_input_triangles(initData(&f_input_triangles,"triangles","Input Triangles"))
  , f_input_quads(initData(&f_input_quads,"quads","Input Quads"))
  , f_input_tetrahedra(initData(&f_input_tetrahedra,"tetrahedra","Input Tetrahedra"))
  , f_input_hexahedra(initData(&f_input_hexahedra,"hexahedra","Input Hexahedra"))

  , f_output_positions(initData(&f_output_positions,"joinedPosition","Output Vertices of the joined mesh"))
  , f_output_edges(initData(&f_output_edges,"joinedEdges","Output Edges of the joined mesh"))
  , f_output_triangles(initData(&f_output_triangles,"joinedTriangles","Output Triangles of the joined mesh"))
  , f_output_quads(initData(&f_output_quads,"joinedQuads","Output Quads of the joined mesh"))
  , f_output_tetrahedra(initData(&f_output_tetrahedra,"joinTetrahedra","Output Tetrahedra of the joined mesh"))
  , f_output_hexahedra(initData(&f_output_hexahedra,"joinHexahedra","Output Hexahedra of the joined mesh"))
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
    addInput(&f_input_positions);
    addInput(&f_input_edges);
    addInput(&f_input_triangles);
    addInput(&f_input_quads);
    addInput(&f_input_tetrahedra);
    addInput(&f_input_hexahedra);

    addOutput(&f_output_positions);
    addOutput(&f_output_edges);
    addOutput(&f_output_triangles);
    addOutput(&f_output_quads);
    addOutput(&f_output_tetrahedra);
    addOutput(&f_output_hexahedra);

    setDirtyValue();
}

template <class DataTypes>
void JoinMeshPoints<DataTypes>::reinit()
{
    update();
}

template <class DataTypes>
void JoinMeshPoints<DataTypes>::update()
{
    cleanDirty();

    const helper::vector< helper::fixed_array <unsigned int,2> >& inJP = f_input_joinPoints.getValue();
	const VecCoord& inPt = f_input_positions.getValue();

	VecCoord& outPt = *f_output_positions.beginEdit();
    outPt.clear();


    // Map analogy to the value of f_input_joinPoints. It maps index of a node
    // from input array to index into output array.
    std::map<unsigned int, unsigned int> mapInIn;

    for (unsigned int i=0; i<inJP.size(); i++) {
        if (inJP[i][0] == inJP[i][1]) continue;
        unsigned int end = inJP[i][1];
        // Traverse the map to find the real index
        do {
            std::map<unsigned int, unsigned int>::const_iterator iter =
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
    for (std::map<unsigned int, unsigned int>::iterator iter = mapInIn.begin();
        iter != mapInIn.end(); iter++) {

        if (iter->first < iter->second) {
            // Invert the mapping
            mapInIn[iter->second] = iter->first;

            std::map<unsigned int, unsigned int>::iterator tmp = iter;
            iter++; // This is safe since iter != mapInIn.end()
            mapInIn.erase(tmp);
            iter--;
        }
    }

    // Create output list
    helper::vector<unsigned int> mapInOut;
    mapInOut.resize(inPt.size());

    for (unsigned int i=0, newId=0; i<inPt.size(); i++) {
        std::map<unsigned int, unsigned int>::const_iterator iter =
            mapInIn.find(i);
        if (iter == mapInIn.end()) {
            outPt.push_back(inPt[i]);
            mapInOut[i] = newId;
            newId++;
        }
    }

    f_output_positions.endEdit();

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
        std::map<unsigned int, unsigned int> mapInIn,
        helper::vector<unsigned int> mapInOut,
        const Data< helper::vector< helper::fixed_array<unsigned int,N> > > &inElements,
        Data< helper::vector< helper::fixed_array<unsigned int,N> > > &outElements)
{
	const helper::vector< helper::fixed_array <unsigned int, N> >& inEle = inElements.getValue();

	helper::vector< helper::fixed_array <unsigned int, N> >& outEle = *outElements.beginEdit();
    outEle.resize(inEle.size());

    for (unsigned int i= 0; i<inEle.size(); i++) {
        helper::fixed_array <unsigned int, N>& out = outEle[i];
        for (unsigned int j=0; j<N; j++) {
            unsigned int id = inEle[i][j];
            std::map<unsigned int, unsigned int>::const_iterator iter =
                mapInIn.find(id);

            // In -> In mapping (joining nodes)
            if (iter != mapInIn.end()) {
                id = iter->second;
            }

            // In -> Out mapping (reindexing nodes)
            if (id < f_input_positions.getValue().size()) {
                out[j] = mapInOut[id];
            } else {
                // Invalid node ID, what now?
                serr << "Invalid node ID in elements! " <<
                    id << " is not a valid node ID!" << sendl; 
                out[j] = 0;
            }
        }
    }

    outElements.endEdit();
}

} // namespace engine

} // namespace component

} // namespace sofa

#endif // #ifndef SOFA_COMPONENT_ENGINE_JOINMESHPOINTS_INL
