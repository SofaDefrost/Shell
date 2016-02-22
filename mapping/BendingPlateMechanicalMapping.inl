/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_MAPPING_MESH2POINTMAPPING_INL
#define SOFA_COMPONENT_MAPPING_MESH2POINTMAPPING_INL

#include "BendingPlateMechanicalMapping.h"
#include <SofaBaseTopology/TriangleSetTopologyContainer.h>
#include <SofaBaseCollision/MinProximityIntersection.h>
#include <sofa/core/visual/VisualParams.h>

#include <SofaBoundaryCondition/ConstantForceField.h>

namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::component::collision;


template <class TIn, class TOut>
void BendingPlateMechanicalMapping<TIn, TOut>::init()
{
//    std::cout << "BendingPlateMechanicalMapping::init()" << std::endl;

    // Retrieves topology
    inputTopo = this->fromModel->getContext()->getMeshTopology();
    outputTopo = this->toModel->getContext()->getMeshTopology();

    if (inputTopo && outputTopo && inputTopo->getNbTriangles() > 0)
    {
        const OutVecCoord &outVertices = this->toModel->read(sofa::core::ConstVecCoordId::position())->getValue();

        listBaseTriangles.clear();
        barycentricCoordinates.clear();
        listBaseTriangles.resize(outVertices.size());
        barycentricCoordinates.resize(outVertices.size());

        // Retrieves 'in' vertices and triangles
        const InVecCoord &inVerticesRigid = this->fromModel->read(sofa::core::ConstVecCoordId::position())->getValue();
        // Conversion to Vec3Types to be able to call same methods used by Hausdorff distance
        OutVecCoord inVertices;
        for (unsigned int i=0; i<inVerticesRigid.size(); i++)
        {
            inVertices.push_back(inVerticesRigid[i].getCenter());
        }
        const SeqEdges &inEdges = inputTopo->getEdges();
        const SeqTriangles &inTriangles = inputTopo->getTriangles();

        // Iterates over 'out' vertices
        Real minimumDistanceVertices, minimumDistanceEdges, minimumDistanceTriangles, minimumDistance;
        int caseToProcess, triangleID;
        Vec3 vertexBaryCoord;
        for (unsigned int i=0; i<outVertices.size(); i++)
        {
            // Iterates over 'in' vertices
            sofa::helper::vector<unsigned int> listClosestVertices;
            minimumDistanceVertices = FindClosestPoints(listClosestVertices, outVertices[i], inVertices);

            // Iterates over 'in' edges
            sofa::helper::vector<unsigned int> listClosestEdges;
            minimumDistanceEdges = FindClosestEdges(listClosestEdges, outVertices[i], inVertices, inEdges);

            // Iterates over 'in' triangles
            sofa::helper::vector<unsigned int> listClosestTriangles;
            minimumDistanceTriangles = FindClosestTriangles(listClosestTriangles, outVertices[i], inVertices, inTriangles);

            // Finds out which type of primitive is the closest
            minimumDistance = std::min(minimumDistanceVertices, std::min(minimumDistanceEdges, minimumDistanceTriangles));


            // Adds the list of triangles attached to the primitives found
            caseToProcess = 0;
            if ( minimumDistance == minimumDistanceVertices )
                caseToProcess = 1;
            if ( minimumDistance == minimumDistanceEdges )
                caseToProcess = 2;
            if ( minimumDistance == minimumDistanceTriangles )
                caseToProcess = 3;

            switch(caseToProcess)
            {
                // If it is a vertex, consider the triangles attached to it
                case 1 :
                    for (unsigned int j=0; j<listClosestVertices.size(); j++)
                    {
                        BaseMeshTopology::TrianglesAroundVertex trianglesAroundVertex = inputTopo->getTrianglesAroundVertex( listClosestVertices[j] );
                        for (unsigned int t=0; t<trianglesAroundVertex.size(); t++)
                        {
                            triangleID = trianglesAroundVertex[t];
                            listBaseTriangles[i].push_back(triangleID);

                            // Computes barycentric coordinates within each triangles
                            Vec3 v1 = inVertices[ inTriangles[triangleID][0] ];
                            Vec3 v2 = inVertices[ inTriangles[triangleID][1] ];
                            Vec3 v3 = inVertices[ inTriangles[triangleID][2] ];

//                            // Computes triangle's normal
//                            Vec3 M = (Vec3) (v2-v1).cross(v3-v1);
//                            double norm2_M = M*(M);
//                            Vec3 N =  M/norm2_M;
//
//                            // Computes projection of the vertex onto this triangle
//                            const Vec3 AP = outVertices[i]-v1;
//                            Vec3 proj = outVertices[i] - N*(AP*N);

                            computeBaryCoefs(vertexBaryCoord, outVertices[i], v1, v2, v3);

                            // Adds the barycentric coordinates to the list
                            barycentricCoordinates[i].push_back(vertexBaryCoord);

                        }
                    }
                    break;


                // If it is an edge, consider the triangles attached to it
                case 2 :
                    for (unsigned int j=0; j<listClosestEdges.size(); j++)
                    {
                        BaseMeshTopology::TrianglesAroundEdge trianglesAroundEdge = inputTopo->getTrianglesAroundEdge( listClosestEdges[j] );
                        for (unsigned int t=0; t<trianglesAroundEdge.size(); t++)
                        {
                            triangleID = trianglesAroundEdge[t];
                            listBaseTriangles[i].push_back(triangleID);

                            // Computes barycentric coordinates within each triangles
                            Vec3 v1 = inVertices[ inTriangles[triangleID][0] ];
                            Vec3 v2 = inVertices[ inTriangles[triangleID][1] ];
                            Vec3 v3 = inVertices[ inTriangles[triangleID][2] ];

//                            // Computes triangle's normal
//                            Vec3 M = (Vec3) (v2-v1).cross(v3-v1);
//                            double norm2_M = M*(M);
//                            Vec3 N =  M/norm2_M;
//
//                            // Computes projection of the vertex onto this triangle
//                            const Vec3 AP = outVertices[i]-v1;
//                            Vec3 proj = outVertices[i] - N*(AP*N);

                            computeBaryCoefs(vertexBaryCoord, outVertices[i], v1, v2, v3);

                            // Adds the barycentric coordinates to the list
                            barycentricCoordinates[i].push_back(vertexBaryCoord);
                         }
                    }
                    break;


                // If it is a triangle, consider the list of triangles
                case 3 :
                    for (unsigned int j=0; j<listClosestTriangles.size(); j++)
                    {
                        listBaseTriangles[i].push_back(listClosestTriangles[j]);

                        // Computes barycentric coordinates within each triangles
                        Vec3 v1 = inVertices[ inTriangles[listClosestTriangles[j]][0] ];
                        Vec3 v2 = inVertices[ inTriangles[listClosestTriangles[j]][1] ];
                        Vec3 v3 = inVertices[ inTriangles[listClosestTriangles[j]][2] ];

//                        // Computes triangle's normal
//                        Vec3 M = (Vec3) (v2-v1).cross(v3-v1);
//                        double norm2_M = M*(M);
//                        Vec3 N =  M/norm2_M;
//
//                        // Computes projection of the vertex onto this triangle
//                        const Vec3 AP = outVertices[i]-v1;
//                        Vec3 proj = outVertices[i] - N*(AP*N);

                        computeBaryCoefs(vertexBaryCoord, outVertices[i], v1, v2, v3);

                        // Adds the barycentric coordinates to the list
                        barycentricCoordinates[i].push_back(vertexBaryCoord);
                    }
                    break;


                default :
                    serr << "BendingPlateMechanicalMapping init(): No closest primitive has been found for vertex " << i << sendl;
                    return;
            }
        }

    }
    else
    {
        serr << "BendingPlateMechanicalMapping requires an input triangular topology" << sendl;
        return;
    }


    // Retrieves Forcefield to compute deflection at runtime
    triangularBendingForcefield = NULL;
    this->getContext()->get(triangularBendingForcefield);
    if (!triangularBendingForcefield)
    {
        serr << "WARNING(BendingPlateMechanicalMapping): triangularBendingForcefield was not found" << sendl;
        return;
    }

#if 0
    // Retrieves topological mapping to retrieve list of edges (to render in wireframe mode)
    triangleSubdivisionTopologicalMapping = NULL;
//    this->getContext()->get(triangleSubdivisionTopologicalMapping, nameHighTopology.getValue(), sofa::core::objectmodel::BaseContext::SearchRoot);
    this->getContext()->get(triangleSubdivisionTopologicalMapping);
    if (!triangleSubdivisionTopologicalMapping)
    {
        serr << "WARNING(BendingPlateMechanicalMapping): triangleSubdivisionTopologicalMapping was not found" << sendl;
        return;
    }
#endif

    // Call of apply() and applyJ()
    this->Inherit::init();

    // Set each colour of each vertex to default
    const OutVecCoord &outVertices = this->toModel->read(sofa::core::ConstVecCoordId::position())->getValue();
    for (unsigned int i=0; i<outVertices.size(); i++)
    {
        coloursPerVertex.push_back(Vec3(0.56, 0.14, 0.6));    // purple
    }

    // If we want to measure the error between the two meshes using Hausdorff distance
    if (measureError.getValue())
    {
        // List of colours to create a colour map
        Vec3 colour;
        Real incr = (float)2/3/240; // (2/3) is chosen stop the gradient to blue
        for (int i=0; i<240; i++)
        {
            HSL2RGB(colour, (float)2/3-i*incr, 0.8, 0.5);
            colourMapping.push_back(colour);
        }

        // Retrieves high resolution mesh topology
        if (targetTopology.get() == NULL) {
            serr << "Missing target topology" << sendl;
        } else {
            // Computes two-sided Hausdorff distance
            MeasureError();
        }

        // Overwrites colour for each vertex based on the error and colour map
        Real maximum = 0;
        // Normalises the error
        for (unsigned int i=0; i<vectorErrorCoarse.size(); i++)
        {
            if (fabs(vectorErrorCoarse[i])>maximum)
            {
                maximum = fabs(vectorErrorCoarse[i]);
            }
        }
        Real correctedError;
        for (unsigned int i=0; i<vectorErrorCoarse.size(); i++)
        {
            correctedError = fabs(vectorErrorCoarse[i])*5;
            if (correctedError > maximum)
                correctedError = maximum;
            coloursPerVertex[i] = colourMapping[ (int)((correctedError/maximum)*239) ];
        }

        // Initialises shader
        shader.InitShaders("shaders/errorMap.vert", "shaders/errorMap.frag");
    }

}


template <class TIn, class TOut>
void BendingPlateMechanicalMapping<TIn, TOut>::reinit()
{
    std::cout << "BendingPlateMechanicalMapping<TIn, TOut>::reinit()" << std::endl;
    init();
}


// Given H,S,L in range of 0-1
// Returns a RGB colour in range of 0-255
// http://www.geekymonkey.com/Programming/CSharp/RGB2HSL_HSL2RGB.htm
template <class TIn, class TOut>
void BendingPlateMechanicalMapping<TIn, TOut>::HSL2RGB(Vec3 &rgb, Real h, Real sl, Real l)
{
    Real v;
    Real r,g,b;

    r = l;   // default to gray
    g = l;
    b = l;
    v = (l <= 0.5) ? (l * (1.0 + sl)) : (l + sl - l * sl);
    if (v > 0)
    {
          Real m;
          Real sv;
          int sextant;
          Real fract, vsf, mid1, mid2;

          m = l + l - v;
          sv = (v - m ) / v;
          h *= 6.0;
          sextant = (int)h;
          fract = h - sextant;
          vsf = v * sv * fract;
          mid1 = m + vsf;
          mid2 = v - vsf;
          switch (sextant)
          {
                case 0:
                      r = v;
                      g = mid1;
                      b = m;
                      break;
                case 1:
                      r = mid2;
                      g = v;
                      b = m;
                      break;
                case 2:
                      r = m;
                      g = v;
                      b = mid1;
                      break;
                case 3:
                      r = m;
                      g = mid2;
                      b = v;
                      break;
                case 4:
                      r = mid1;
                      g = m;
                      b = v;
                      break;
                case 5:
                      r = v;
                      g = m;
                      b = mid2;
                      break;
          }
    }

    rgb[0] = r;
    rgb[1] = g;
    rgb[2] = b;
}


template <class TIn, class TOut>
void BendingPlateMechanicalMapping<TIn, TOut>::MeasureError()
{
    Real distance1;
    std::cout << "Computing Hausdorff distance high res->coarse" << std::endl;
    distance1 = DistanceHausdorff(targetTopology.get(), outputTopo, vectorErrorTarget);
    std::cout << "Hausdorff distance between high res mesh and coarse mesh = " << distance1 << std::endl;

    Real average = 0;
    for (unsigned int i=0; i<vectorErrorTarget.size(); i++)
    {
        average += vectorErrorTarget[i];
    }
    std::cout << "Mean Hausdorff distance = " << average/vectorErrorTarget.size() << std::endl;



    Real distance2;
    std::cout << "Computing Hausdorff distance coarse->high res" << std::endl;
    distance2 = DistanceHausdorff(outputTopo, targetTopology.get(), vectorErrorCoarse);
    std::cout << "Hausdorff distance between coarse mesh and high res mesh = " << distance2 << std::endl;

    average = 0;
    for (unsigned int i=0; i<vectorErrorCoarse.size(); i++)
    {
        average += vectorErrorCoarse[i];
    }
    std::cout << "Mean Hausdorff distance = " << average/vectorErrorCoarse.size() << std::endl;

}

template <class TIn, class TOut>
typename BendingPlateMechanicalMapping<TIn, TOut>::Real BendingPlateMechanicalMapping<TIn, TOut>::DistanceHausdorff(BaseMeshTopology *topo1, BaseMeshTopology *topo2, helper::vector<Real> &vectorError)
{
    // Mesh 1
    MechanicalState<Out>* mState1 = dynamic_cast<MechanicalState<Out>*> (topo1->getContext()->getMechanicalState());
    const OutVecCoord &vertices1 = mState1->read(sofa::core::ConstVecCoordId::position())->getValue();

    // Mesh 2
    MechanicalState<Out>* mState2 = dynamic_cast<MechanicalState<Out>*> (topo2->getContext()->getMechanicalState());
    const OutVecCoord &vertices2 = mState2->read(sofa::core::ConstVecCoordId::position())->getValue();
    const SeqEdges edges2 = topo2->getEdges();
    const SeqTriangles triangles2 = topo2->getTriangles();

    // Dummy list (lists of primitives are useless here)
    sofa::helper::vector<unsigned int> listDummy;

    // Iterates over 'mesh1' vertices
    Real minimumDistanceVertices, minimumDistanceEdges, minimumDistanceTriangles, minimumDistance;
    Real HausdorffDistance = -1;
    for (unsigned int i=0; i<vertices1.size(); i++)
    {
        // Iterates over 'mesh2' vertices
        minimumDistanceVertices = FindClosestPoints(listDummy, vertices1[i], vertices2);

        // Iterates over 'mesh2' edges
        minimumDistanceEdges = FindClosestEdges(listDummy, vertices1[i], vertices2, edges2);

        // Iterates over 'mesh2' triangles
        minimumDistanceTriangles = FindClosestTriangles(listDummy, vertices1[i], vertices2, triangles2);

        // Finds out which type of primitive is the closest
        minimumDistance = std::min(minimumDistanceVertices, std::min(minimumDistanceEdges, minimumDistanceTriangles));
        // And stores the distance between the vertex of mesh1 and mesh2
        vectorError.push_back(minimumDistance);

        // The maximum distance is the Hausdorff distance between mesh1 and mesh2
        if (minimumDistance > HausdorffDistance)
        {
            HausdorffDistance = minimumDistance;
        }
    }

    return HausdorffDistance;
}


// --------------------------------------------------------------------------------------
// Finds the list of the closest points to a point
// --------------------------------------------------------------------------------------
template <class TIn, class TOut>
typename BendingPlateMechanicalMapping<TIn, TOut>::Real BendingPlateMechanicalMapping<TIn, TOut>::FindClosestPoints(sofa::helper::vector<unsigned int>& listClosestVertices, const Vec3& point, const OutVecCoord &inVertices)
{
    Real minimumDistance = 10e12;
    for (unsigned int v=0; v<inVertices.size(); v++)
    {
        Real distance = (inVertices[v] - point).norm2();

        Real threshold = 1e-12;
        if (distance < minimumDistance)
        {
            // We deal with a new minimum, so clear the previous list
            listClosestVertices.clear();
            // Adds the new minimum
            listClosestVertices.push_back(v);
            // Updates the minimum's value
            minimumDistance = distance;
        }
        else if ( distance-minimumDistance < threshold )
        {
            // Adds the new minimum
            listClosestVertices.push_back(v);
        }
    }

    return minimumDistance;
}


// --------------------------------------------------------------------------------------
// Finds the list of the closest edges to a point
// --------------------------------------------------------------------------------------
template <class TIn, class TOut>
typename BendingPlateMechanicalMapping<TIn, TOut>::Real BendingPlateMechanicalMapping<TIn, TOut>::FindClosestEdges(sofa::helper::vector<unsigned int>& listClosestEdges, const Vec3& point, const OutVecCoord &inVertices, const SeqEdges &inEdges)
{
    Real minimumDistance = 10e12;
    for (unsigned int e=0; e<inEdges.size(); e++)
    {
        Vec3 pointEdge1 = inVertices[ inEdges[e][0] ];
        Vec3 pointEdge2 = inVertices[ inEdges[e][1] ];

        const Vec3 AB = pointEdge2-pointEdge1;
        const Vec3 AP = point-pointEdge1;

        double A;
        double b;
        A = AB*AB;
        b = AP*AB;

        double alpha = b/A;

        // If the point is on the edge
        if (alpha >= 0 && alpha <= 1)
        {
            Vec3 P, Q, PQ;
            P = point;
            Q = pointEdge1 + AB * alpha;
            PQ = Q-P;

            Real distance = PQ.norm2();
            Real threshold = 1e-12;
            if (distance < minimumDistance)
            {
                // We deal with a new minimum, so clear the previous list
                listClosestEdges.clear();
                // Adds the new minimum
                listClosestEdges.push_back(e);
                // Updates the minimum's value
                minimumDistance = distance;
            }
            else if ( distance-minimumDistance < threshold )
            {
                // Adds the new minimum
                listClosestEdges.push_back(e);
            }
        }

    }

    return minimumDistance;
}


// --------------------------------------------------------------------------------------
// Finds the list of the closest triangles to a point
// --------------------------------------------------------------------------------------
template <class TIn, class TOut>
typename BendingPlateMechanicalMapping<TIn, TOut>::Real BendingPlateMechanicalMapping<TIn, TOut>::FindClosestTriangles(sofa::helper::vector<unsigned int>& listClosestTriangles, const Vec3& point, const OutVecCoord &inVertices, const SeqTriangles &inTriangles)
{
    Real minimumDistance = 10e12;
    for (unsigned int t=0; t<inTriangles.size(); t++)
    {
        Vec3 pointTriangle1 = inVertices[ inTriangles[t][0] ];
        Vec3 pointTriangle2 = inVertices[ inTriangles[t][1] ];
        Vec3 pointTriangle3 = inVertices[ inTriangles[t][2] ];

        const Vector3 AB = pointTriangle2-pointTriangle1;
        const Vector3 AC = pointTriangle3-pointTriangle1;
        const Vector3 AP = point-pointTriangle1;
        Matrix2 A;
        Vector2 b;

        // We want to find alpha,beta so that:
        // AQ = AB*alpha+AC*beta
        // PQ.AB = 0 and PQ.AC = 0
        // (AQ-AP).AB = 0 and (AQ-AP).AC = 0
        // AQ.AB = AP.AB and AQ.AC = AP.AC
        //
        // (AB*alpha+AC*beta).AB = AP.AB and
        // (AB*alpha+AC*beta).AC = AP.AC
        //
        // AB.AB*alpha + AC.AB*beta = AP.AB and
        // AB.AC*alpha + AC.AC*beta = AP.AC
        //
        // A . [alpha beta] = b
        A[0][0] = AB*AB;
        A[1][1] = AC*AC;
        A[0][1] = A[1][0] = AB*AC;
        b[0] = AP*AB;
        b[1] = AP*AC;
        const double det = determinant(A);

        double alpha = (b[0]*A[1][1] - b[1]*A[0][1])/det;
        double beta  = (b[1]*A[0][0] - b[0]*A[1][0])/det;

        // If point is on one of the edge, returns
        Real threshold = 1e-12;
        if (alpha >= 0 && beta >= 0 && alpha + beta <= 1 )
        {
            const Vector3 PQ = AB * alpha + AC * beta - AP;

            Real distance = PQ.norm2();
            if (distance < minimumDistance)
            {
                // We deal with a new minimum, so clear the previous list
                listClosestTriangles.clear();
                // Adds the new minimum
                listClosestTriangles.push_back(t);
                // Updates the minimum's value
                minimumDistance = distance;
            }
            else if ( distance-minimumDistance < threshold )
            {
                // Adds the new minimum
                listClosestTriangles.push_back(t);
            }
        }
    }

    return minimumDistance;
}


// --------------------------------------------------------------------------------------
// Barycentric coefficients of point p in triangle whose vertices are a, b and c
// --------------------------------------------------------------------------------------
template <class TIn, class TOut>
void BendingPlateMechanicalMapping<TIn, TOut>::computeBaryCoefs(Vec3 &baryCoefs, const Vec3 &p, const Vec3 &a, const Vec3 &b, const Vec3 &c)
{
    const double ZERO = 1e-20;

    Vec3 M = (Vec3) (b-a).cross(c-a);
    double norm2_M = M*(M);

    double coef_a, coef_b, coef_c;

    //if(norm2_M==0.0) // triangle (a,b,c) is flat
    if(norm2_M < ZERO) // triangle (a,b,c) is flat
    {
        coef_a = (double) (1.0/3.0);
        coef_b = (double) (1.0/3.0);
        coef_c = (double) (1.0 - (coef_a + coef_b));
    }
    else
    {
        Vec3 N =  M/norm2_M;

        coef_a = N*((b-p).cross(c-p));
        coef_b = N*((c-p).cross(a-p));
        coef_c = (double) (1.0 - (coef_a + coef_b)); //N*((a-p).cross(b-p));
    }

    baryCoefs[0] = coef_a;
    baryCoefs[1] = coef_b;
    baryCoefs[2] = coef_c;
}


// Updates positions of the visual mesh from mechanical vertices
template <class TIn, class TOut>
//void BendingPlateMechanicalMapping<TIn, TOut>::apply( typename Out::VecCoord& out, const typename In::VecCoord& in )
void BendingPlateMechanicalMapping<TIn, TOut>::apply(const core::MechanicalParams * /*mparams*/, Data<OutVecCoord>& dOut, const Data<InVecCoord>& dIn)
{

    helper::WriteAccessor< Data<OutVecCoord> > out = dOut;
    helper::ReadAccessor< Data<InVecCoord> > in = dIn;


//    std::cout << "---------------- Apply ----------------------------" << std::endl;

//    sofa::helper::system::thread::ctime_t start, stop;
//    sofa::helper::system::thread::CTime timer;
//
//    start = timer.getTime();

    if (!inputTopo || !outputTopo)
    {
        serr << "BendingPlateMechanicalMapping apply() was called before init()" << sendl;
        return;
    }
    if (inputTopo->getNbTriangles() <= 0)
    {
        serr << "BendingPlateMechanicalMapping apply() requires an input triangular topology" << sendl;
        return;
    }

    if (!triangularBendingForcefield)
    {
        serr << "No TriangularBendingForcefield has been found" << sendl;
        this->getContext()->get(triangularBendingForcefield);
        return;
    }
    else
    {
        helper::vector<TriangleInformation>& triangleInf = *(triangularBendingForcefield->getTriangleInfo().beginEdit());
        TriangleInformation *tinfo = NULL;

        // List of in triangles
        const SeqTriangles& inTriangles = inputTopo->getTriangles();

        // Computes the coefficients ci for each triangle
        for (unsigned int t=0; t<inTriangles.size();t++)
        {
            tinfo = &triangleInf[t];
            tinfo->coefficients = tinfo->invC * (tinfo->u + tinfo->u_rest);
        }

        Vec3 a, b, c, baryCoord, vertexLocal;
        Real z;
        for (unsigned int i=0; i<out.size(); i++)
        {
            // Gets the first triangle that the vertex belongs to
            Triangle triangle = inTriangles[ listBaseTriangles[i][0] ];

            // Gets its 3 vertices
            a = in[ triangle[0] ].getCenter();
            b = in[ triangle[1] ].getCenter();
            c = in[ triangle[2] ].getCenter();

            baryCoord = barycentricCoordinates[i][0];
            out[i] = a*baryCoord[0] + b*baryCoord[1] + c*baryCoord[2];

            Vec3 Uz(0.0, 0.0, 0.0);
            Vec3 w(0, 0, 0);
            for (unsigned int t=0; t<listBaseTriangles[i].size();t++)
            {
                triangle = triangularBendingForcefield->getTopology()->getTriangle(listBaseTriangles[i][t]);
                tinfo = &triangleInf[listBaseTriangles[i][t]];

                // Local coordinates needed to compute deflection
                a = in[ triangle[0] ].getCenter();
                vertexLocal = tinfo->Qframe.rotate(out[i]-a);

                // Adds deflection
                z = tinfo->coefficients[0] + tinfo->coefficients[1]*vertexLocal[0] + tinfo->coefficients[2]*vertexLocal[1] + tinfo->coefficients[3]*vertexLocal[0]*vertexLocal[0] + tinfo->coefficients[4]*vertexLocal[0]*vertexLocal[1] + tinfo->coefficients[5]*vertexLocal[1]*vertexLocal[1] + tinfo->coefficients[6]*vertexLocal[0]*vertexLocal[0]*vertexLocal[0] + tinfo->coefficients[7]*vertexLocal[0]*vertexLocal[1]*vertexLocal[1] + tinfo->coefficients[8]*vertexLocal[1]*vertexLocal[1]*vertexLocal[1];

                w += tinfo->Qframe.inverseRotate(Vec3(0, 0, z));
            }

            out[i] += w/listBaseTriangles[i].size();
        }

        triangularBendingForcefield->getTriangleInfo().endEdit();
    }

//    stop = timer.getTime();
//    std::cout << "time apply = " << stop-start << std::endl;
}


// Updates velocities of the visual mesh from mechanical vertices
template <class TIn, class TOut>
void BendingPlateMechanicalMapping<TIn, TOut>::applyJ(const core::MechanicalParams* /*mparams*/, Data<OutVecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    helper::WriteAccessor< Data<OutVecDeriv> > out = dOut;
    helper::ReadAccessor< Data<InVecDeriv> > in = dIn;

//    std::cout << "---------------- ApplyJ ----------------------------" << std::endl;

//    sofa::helper::system::thread::ctime_t start, stop;
//    sofa::helper::system::thread::CTime timer;
//
//    start = timer.getTime();

    if (!inputTopo || !outputTopo)
    {
        serr << "BendingPlateMechanicalMapping applyJ() was called before init()" << sendl;
        return;
    }
    if (inputTopo->getNbTriangles() <= 0)
    {
        serr << "BendingPlateMechanicalMapping applyJ() requires an input triangular topology" << sendl;
        return;
    }

    if (!triangularBendingForcefield)
    {
        serr << "No TriangularBendingForcefield has been found" << sendl;
        this->getContext()->get(triangularBendingForcefield);
        return;
    }
    else
    {
        helper::vector<TriangleInformation>& triangleInf = *(triangularBendingForcefield->getTriangleInfo().beginEdit());
        TriangleInformation *tinfo = NULL;

         // List of 'in' positions
        const InVecCoord &inVertices = this->fromModel->read(sofa::core::ConstVecCoordId::position())->getValue();

        // List of 'out' positions
        const OutVecCoord &outVertices = this->toModel->read(sofa::core::ConstVecCoordId::position())->getValue();

        // List of in triangles
        const SeqTriangles& inTriangles = inputTopo->getTriangles();

        // Computes the coefficients ci for each triangle
        Vec3 va_a, va_b, va_c, va_a_local, va_b_local, va_c_local;
        Vec <9, Real> v_u;
        for (unsigned int t=0; t<inTriangles.size();t++)
        {
            Triangle triangle = triangularBendingForcefield->getTopology()->getTriangle(t);
            tinfo = &triangleInf[t];

            // Gets the angular velocities of each vertex
            va_a = getVOrientation(in[ triangle[0] ]);
            va_b = getVOrientation(in[ triangle[1] ]);
            va_c = getVOrientation(in[ triangle[2] ]);
            // In local frame
            va_a_local = tinfo->Qframe.rotate(va_a);
            va_b_local = tinfo->Qframe.rotate(va_b);
            va_c_local = tinfo->Qframe.rotate(va_c);
            // Fills in du/dt
            v_u.clear();
            v_u[1] = va_a_local[0];   v_u[2] = va_a_local[1];
            v_u[4] = va_b_local[0];   v_u[5] = va_b_local[1];
            v_u[7] = va_c_local[0];   v_u[8] = va_c_local[1];

            tinfo->coefficients = tinfo->invC * v_u;
        }


        // Iterates over out vertices to update coordinates
        Vec3 v_a, v_b, v_c, baryCoord, a, vertexLocal;
        Real v_z;
        for (unsigned int i=0; i<out.size(); i++)
        {
            // Gets the first triangle that the vertex belongs to
            Triangle triangle = inTriangles[ listBaseTriangles[i][0] ];

            // Gets the linear velocities of each vertex
            v_a = getVCenter(in[ triangle[0] ]);
            v_b = getVCenter(in[ triangle[1] ]);
            v_c = getVCenter(in[ triangle[2] ]);

            baryCoord = barycentricCoordinates[i][0];
            out[i] = v_a*baryCoord[0] + v_b*baryCoord[1] + v_c*baryCoord[2];

            Vec3 Uz(0.0, 0.0, 0.0);
            Real w = 0;
            for (unsigned int t=0; t<listBaseTriangles[i].size();t++)
            {
                triangle = triangularBendingForcefield->getTopology()->getTriangle(listBaseTriangles[i][t]);
                tinfo = &triangleInf[listBaseTriangles[i][t]];

                // Local coordinates needed to compute deflection
                a = inVertices[ triangle[0] ].getCenter();
                vertexLocal = tinfo->Qframe.rotate(outVertices[i]-a);

                // Adds deflection velocity
                v_z = tinfo->coefficients[0] + tinfo->coefficients[1]*vertexLocal[0] + tinfo->coefficients[2]*vertexLocal[1] + tinfo->coefficients[3]*vertexLocal[0]*vertexLocal[0] + tinfo->coefficients[4]*vertexLocal[0]*vertexLocal[1] + tinfo->coefficients[5]*vertexLocal[1]*vertexLocal[1] + tinfo->coefficients[6]*vertexLocal[0]*vertexLocal[0]*vertexLocal[0] + tinfo->coefficients[7]*vertexLocal[0]*vertexLocal[1]*vertexLocal[1] + tinfo->coefficients[8]*vertexLocal[1]*vertexLocal[1]*vertexLocal[1];
                w += v_z;

            }

            // Position in Qframe
            Vec3 out0 = tinfo->Qframe.rotate(out[i]);
            // Computed deflection w in Qframe
            out0[2] += w/listBaseTriangles[i].size();
            out[i] = tinfo->Qframe.inverseRotate(out0);
        }

        triangularBendingForcefield->getTriangleInfo().endEdit();

    }

//    stop = timer.getTime();
//    std::cout << "time applyJ = " << stop-start << std::endl;
}


// Updates positions of the mechanical vertices from visual    f(n-1) = JT * fn
template <class TIn, class TOut>
void BendingPlateMechanicalMapping<TIn, TOut>::applyJT(const core::MechanicalParams * /*mparams*/, Data<InVecDeriv>& dOut, const Data<OutVecDeriv>& dIn)
{
    helper::WriteAccessor< Data<InVecDeriv> > out = dOut;
    helper::ReadAccessor< Data<OutVecDeriv> > in = dIn;

//    std::cout << "---------------- ApplyJT ----------------------------" << std::endl;

//    sofa::helper::system::thread::ctime_t start, stop;
//    sofa::helper::system::thread::CTime timer;
//
//    start = timer.getTime();

    if (!inputTopo || !outputTopo)
    {
        serr << "BendingPlateMechanicalMapping applyJT() was called before init()" << sendl;
        return;
    }
    if (inputTopo->getNbTriangles() <= 0)
    {
        serr << "BendingPlateMechanicalMapping applyJT() requires an input triangular topology" << sendl;
        return;
    }

    if (!triangularBendingForcefield)
    {
        serr << "No TriangularBendingForcefield has been found" << sendl;
        this->getContext()->get(triangularBendingForcefield);
        return;
    }
    else
    {
        helper::vector<TriangleInformation>& triangleInf = *(triangularBendingForcefield->getTriangleInfo().beginEdit());
        TriangleInformation *tinfo = NULL;

        // List of in triangles
        const SeqTriangles& inTriangles = inputTopo->getTriangles();

        sofa::helper::vector< Mat<9, 9, Real> > vecTransposedInvC;
        for (unsigned int t=0; t<inTriangles.size();t++)
        {
            tinfo = &triangleInf[t];
            vecTransposedInvC.push_back(tinfo->invC.transposed());
        }

        // List of 'in' positions
        const OutVecCoord &inVertices = this->toModel->read(sofa::core::ConstVecCoordId::position())->getValue();

        // List of 'out' positions (mechanical points)
        const InVecCoord &outVertices = this->fromModel->read(sofa::core::ConstVecCoordId::position())->getValue();

        // Iterates over in vertices
        Triangle triangle;
        Vec3 baryCoord, inLocal, a, vertexLocal, torqueA, torqueB, torqueC;
        Real Fz;
        Vec<9, Real> polynomDeflection, a_u;
        for (unsigned int i=0; i<in.size(); i++)
        {
            triangle = inTriangles[ listBaseTriangles[i][0] ];
            tinfo = &triangleInf[listBaseTriangles[i][0]];

            // Linear acceleration
            baryCoord = barycentricCoordinates[i][0];
            getVCenter(out[ triangle[0] ]) += in[i] * baryCoord[0];
            getVCenter(out[ triangle[1] ]) += in[i] * baryCoord[1];
            getVCenter(out[ triangle[2] ]) += in[i] * baryCoord[2];

            // Iterates over triangles
            for (unsigned int t=0; t<listBaseTriangles[i].size();t++)
            {
                triangle = triangularBendingForcefield->getTopology()->getTriangle(listBaseTriangles[i][t]);
                tinfo = &triangleInf[listBaseTriangles[i][t]];

                // Applied force into local frame
                inLocal = tinfo->Qframe.rotate(in[i]);
                Fz = inLocal[2];

                if (Fz != 0)
                {
                    // Local coordinates needed to compute deflection
                    a = outVertices[ triangle[0] ].getCenter();
                    vertexLocal = tinfo->Qframe.rotate(inVertices[i]-a); // WARNING: SHOULD NOT WE NEED TO PROJECT THE INVERTICES INTO THE TRIANGLE'S PLAN FIRST?

                    // Uz = c1 + c2*x+ c3*y + c4*x^2 + c5*x*y + c6*y^2 + c7*x^3 + c8*x*y^2 + c9*y^3
                    polynomDeflection[0] = 1;
                    polynomDeflection[1] = vertexLocal[0];
                    polynomDeflection[2] = vertexLocal[1];
                    polynomDeflection[3] = vertexLocal[0]*vertexLocal[0];
                    polynomDeflection[4] = vertexLocal[0]*vertexLocal[1];
                    polynomDeflection[5] = vertexLocal[1]*vertexLocal[1];
                    polynomDeflection[6] = vertexLocal[0]*vertexLocal[0]*vertexLocal[0];
                    polynomDeflection[7] = vertexLocal[0]*vertexLocal[1]*vertexLocal[1];
                    polynomDeflection[8] = vertexLocal[1]*vertexLocal[1]*vertexLocal[1];

                    polynomDeflection *= Fz;

                    // Moments at each point
                    a_u = (vecTransposedInvC[ listBaseTriangles[i][t] ] * polynomDeflection) / listBaseTriangles[i].size();

                    // Moments into global frame
                    torqueA = tinfo->Qframe.inverseRotate(Vec3(a_u[1], a_u[2], 0));
                    torqueB = tinfo->Qframe.inverseRotate(Vec3(a_u[4], a_u[5], 0));
                    torqueC = tinfo->Qframe.inverseRotate(Vec3(a_u[7], a_u[8], 0));

                    getVOrientation(out[ triangle[0] ]) += torqueA;
                    getVOrientation(out[ triangle[1] ]) += torqueB;
                    getVOrientation(out[ triangle[2] ]) += torqueC;
                }
            }
        }

    }


//    stop = timer.getTime();
//    std::cout << "time applyJT = " << stop-start << std::endl;
}


template <class TIn, class TOut>
void BendingPlateMechanicalMapping<TIn, TOut>::applyJT(const core::ConstraintParams * /*cparams*/, Data<InMatrixDeriv>& /*dOut*/, const Data<OutMatrixDeriv>& /*dIn*/)
{

}


template <class TIn, class TOut>
void BendingPlateMechanicalMapping<TIn, TOut>::draw(const core::visual::VisualParams* vparams)
{
    const OutVecCoord &outVertices = this->toModel->read(sofa::core::ConstVecCoordId::position())->getValue();

    if(vparams->displayFlags().getShowVisualModels())
    {
        glDisable(GL_LIGHTING);

        if (measureError.getValue())
        {
            shader.TurnOn();


            const SeqTriangles &outTriangles = outputTopo->getTriangles();
            unsigned int index;

            glEnable(GL_DEPTH_TEST);
            glPolygonOffset(1.0, 1.0);

            if(vparams->displayFlags().getShowWireFrame())
            {
                glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                glEnable(GL_POLYGON_OFFSET_LINE);
                glColor4f(0.0, 0.0, 1.0, 1.0);
                glBegin(GL_TRIANGLES);
                for (unsigned int i=0; i<outTriangles.size(); i++)
                {
                    index = outTriangles[i][0];
                    glVertex3f(outVertices[index][0], outVertices[index][1], outVertices[index][2]);

                    index = outTriangles[i][1];
                    glVertex3f(outVertices[index][0], outVertices[index][1], outVertices[index][2]);

                    index = outTriangles[i][2];
                    glVertex3f(outVertices[index][0], outVertices[index][1], outVertices[index][2]);
                }
                glEnd();
                glDisable(GL_POLYGON_OFFSET_LINE);
            }
            else
            {
                glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                glEnable(GL_POLYGON_OFFSET_FILL);
                //            glColor4f(0.0, 0.0, 0.0, 0.0);
                glBegin(GL_TRIANGLES);
                for (unsigned int i=0; i<outTriangles.size(); i++)
                {
                    index = outTriangles[i][0];
                    glColor4f(coloursPerVertex[index][0], coloursPerVertex[index][1], coloursPerVertex[index][2], 1.0);
                    glVertex3f(outVertices[index][0], outVertices[index][1], outVertices[index][2]);

                    index = outTriangles[i][1];
                    glColor4f(coloursPerVertex[index][0], coloursPerVertex[index][1], coloursPerVertex[index][2], 1.0);
                    glVertex3f(outVertices[index][0], outVertices[index][1], outVertices[index][2]);

                    index = outTriangles[i][2];
                    glColor4f(coloursPerVertex[index][0], coloursPerVertex[index][1], coloursPerVertex[index][2], 1.0);
                    glVertex3f(outVertices[index][0], outVertices[index][1], outVertices[index][2]);
                }
                glEnd();
                glDisable(GL_POLYGON_OFFSET_FILL);
            }

            shader.TurnOff();
        }

#if 0
        if (!triangleSubdivisionTopologicalMapping) return;

        const SeqEdges &outEdges = triangleSubdivisionTopologicalMapping->getSubEdges();
//        const SeqEdges &outEdges = outputTopo->getEdges();

        // Render shells' contours (subdivision of edges)
        glColor4f(1.0, 1.0, 1.0, 1.0);
//        glColor4f(0.0, 0.0, 0.0, 1.0);
        glLineWidth(0.5);
        glBegin(GL_LINES);
        for (unsigned int i=0; i<outEdges.size(); i++)
        {
            index = outEdges[i][0];
            glVertex3f(outVertices[index][0], outVertices[index][1], outVertices[index][2]);

            index = outEdges[i][1];
            glVertex3f(outVertices[index][0], outVertices[index][1], outVertices[index][2]);
        }
        glEnd();
#endif

    }
}


} // namespace mapping

} // namespace component

} // namespace sofa

#endif
