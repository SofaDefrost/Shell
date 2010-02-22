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
#include <sofa/component/topology/TriangleSetTopologyContainer.h>
#include <sofa/component/collision/MinProximityIntersection.h>

#include <sofa/component/forcefield/ConstantForceField.h>

namespace sofa
{

namespace component
{

namespace mapping
{

using namespace	sofa::component::collision;

template <class BaseMapping>
BendingPlateMechanicalMapping<BaseMapping>::BendingPlateMechanicalMapping(In* from, Out* to)
: Inherit(from, to)
, inputTopo(NULL)
, outputTopo(NULL)
, measureError(initData(&measureError, false, "measureError","Error with high resolution mesh"))
, targetVertices(initData(&targetVertices, "targetVertices","Vertices of the target mesh"))
, targetTriangles(initData(&targetTriangles, "targetTriangles","Triangles of the target mesh"))
{
}

template <class BaseMapping>
BendingPlateMechanicalMapping<BaseMapping>::~BendingPlateMechanicalMapping()
{
}

template <class BaseMapping>
void BendingPlateMechanicalMapping<BaseMapping>::init()
{
    // Retrieves topology
    inputTopo = this->fromModel->getContext()->getMeshTopology();
    outputTopo = this->toModel->getContext()->getMeshTopology();

    if (inputTopo && outputTopo && inputTopo->getNbTriangles() > 0)
    {
        OutVecCoord &outVertices = *this->toModel->getX();
        listBaseTriangles.clear();
        barycentricCoordinates.clear();
        listBaseTriangles.resize(outVertices.size());
        barycentricCoordinates.resize(outVertices.size());

        // Retrieves 'in' vertices and triangles
        InVecCoord &inVertices = *this->fromModel->getX();
        const SeqTriangles &inTriangles = inputTopo->getTriangles();

        // Iterates over 'out' vertices
        Real minimumDistanceVertices, minimumDistanceEdges, minimumDistanceTriangles, minimumDistance;
        int caseToProcess, triangleID;
        Vec3 vertexBaryCoord;
        for (unsigned int i=0; i<outVertices.size(); i++)
        {
            // Iterates over 'in' vertices
            sofa::helper::vector<unsigned int> listClosestVertices;
            minimumDistanceVertices = FindClosestPoints(outVertices[i], listClosestVertices);
            
            // Iterates over 'in' edges
            sofa::helper::vector<unsigned int> listClosestEdges;
            minimumDistanceEdges = FindClosestEdges(outVertices[i], listClosestEdges);

            // Iterates over 'in' triangles
            sofa::helper::vector<unsigned int> listClosestTriangles;
            minimumDistanceTriangles = FindClosestTriangles(outVertices[i], listClosestTriangles);

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
                        TrianglesAroundVertex trianglesAroundVertex = inputTopo->getTrianglesAroundVertex( listClosestVertices[j] );
                        for (unsigned int t=0; t<trianglesAroundVertex.size(); t++)
                        {
                            triangleID = trianglesAroundVertex[t];
                            listBaseTriangles[i].push_back(triangleID);

                            // Computes barycentric coordinates within each triangles
                            Vec3 v1 = inVertices[ inTriangles[triangleID][0] ].getCenter();
                            Vec3 v2 = inVertices[ inTriangles[triangleID][1] ].getCenter();
                            Vec3 v3 = inVertices[ inTriangles[triangleID][2] ].getCenter();

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
                        TrianglesAroundEdge trianglesAroundEdge = inputTopo->getTrianglesAroundEdge( listClosestEdges[j] );
                        for (unsigned int t=0; t<trianglesAroundEdge.size(); t++)
                        {
                            triangleID = trianglesAroundEdge[t];
                            listBaseTriangles[i].push_back(triangleID);

                            // Computes barycentric coordinates within each triangles
                            Vec3 v1 = inVertices[ inTriangles[triangleID][0] ].getCenter();
                            Vec3 v2 = inVertices[ inTriangles[triangleID][1] ].getCenter();
                            Vec3 v3 = inVertices[ inTriangles[triangleID][2] ].getCenter();

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


                // If it is a vertex, consider the list of triangles
                case 3 :
                    for (unsigned int j=0; j<listClosestTriangles.size(); j++)
                    {
                        listBaseTriangles[i].push_back(listClosestTriangles[j]);

                        // Computes barycentric coordinates within each triangles
                        Vec3 v1 = inVertices[ inTriangles[listClosestTriangles[j]][0] ].getCenter();
                        Vec3 v2 = inVertices[ inTriangles[listClosestTriangles[j]][1] ].getCenter();
                        Vec3 v3 = inVertices[ inTriangles[listClosestTriangles[j]][2] ].getCenter();

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
                    serr << "BendingPlateMechanicalMapping init(): No closest primitive has been found" << sendl;
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
        return;

    // Retrieves topological mapping to retrieve list of edges (for rendering of wireframe)
    triangleSubdivisionTopologicalMapping = NULL;
    this->getContext()->get(triangleSubdivisionTopologicalMapping);
    if (!triangleSubdivisionTopologicalMapping)
        return;

    // Call of apply() and applyJ()
    this->Inherit::init();


    // Retrieves high resolution mesh topology and measures the error
    if (measureError.getValue())
    {
        MeasureError();
    }
    else
    {
        OutVecCoord &outVertices = *this->toModel->getX();
        for (unsigned int i=0; i<outVertices.size(); i++)
        {
            vectorError.push_back(0);
            normals.push_back(Vec3(0,0,0));
        }
    }


}


template <class BaseMapping>
void BendingPlateMechanicalMapping<BaseMapping>::reinit()
{
    std::cout << "BendingPlateMechanicalMapping<BaseMapping>::reinit()" << std::endl;
    init();
}


template <class BaseMapping>
void BendingPlateMechanicalMapping<BaseMapping>::MeasureError()
{
    // Retrieves vertices of high resolution mesh
    const InVecCoord& highResVertices = targetVertices.getValue();
    // Retrieves triangles of high resolution mesh
    const SeqTriangles highRestriangles = targetTriangles.getValue();

    // Computes normal for each subvertex
//    helper::vector<Vec3> normals;
    ComputeNormals(normals);
    
    // Finds the triangle on high res mesh towards normal direction for each point and measure the distance to it
    FindTriangleInNormalDirection(highResVertices, highRestriangles, normals);
}


template <class BaseMapping>
void BendingPlateMechanicalMapping<BaseMapping>::ComputeNormals(helper::vector<Vec3> &normals)
{
    // Retrieves vertices and triangles of subdivided coarse mesh
    OutVecCoord &outVertices = *this->toModel->getX();
    const SeqTriangles &outTriangles = outputTopo->getTriangles();

    /**
     * Computes normal at each subvertex
     */
    for (unsigned int i=0; i<outVertices.size(); i++)
    {
        normals.push_back(Vec3(0, 0, 0));
    }

    for (unsigned int t=0; t<outTriangles.size(); t++)
    {
        Vec3 a = outVertices[ outTriangles[t][0] ];
        Vec3 b = outVertices[ outTriangles[t][1] ];
        Vec3 c = outVertices[ outTriangles[t][2] ];

        Vec3 z = cross(b-a, c-a);
        z.normalize();

        normals[ outTriangles[t][0] ] += z;
        normals[ outTriangles[t][1] ] += z;
        normals[ outTriangles[t][2] ] += z;
    }

    for (unsigned int i=0; i<normals.size(); i++)
    {
        normals[i].normalize();
    }
}


// ---------------------------------------------------------------------------------------------
// Finds the intersection between the normal of subvertex and targeted surface
//
// http://en.wikipedia.org/wiki/Line-plane_intersection
// A triangle with its 3 vertices P1, P2 and P3 defines a plan. A subvertex P and its normal
// defines a line. The intersection can be parametered with t (distance along the line from P)
// and (u,v) the coordinates within the triangles (u along side P0P1 and v along P0P2).
//
// ---------------------------------------------------------------------------------------------
template <class BaseMapping>
void BendingPlateMechanicalMapping<BaseMapping>::FindTriangleInNormalDirection(const InVecCoord& highResVertices, const SeqTriangles highRestriangles, const helper::vector<Vec3> &normals)
{
    OutVecCoord &outVertices = *this->toModel->getX();
    
    Real minimumDistance, distance, error;
    Vec3 point, normal;

    // Iterates over each subvertex
    for (unsigned int i=0; i<outVertices.size(); i++)
    {
        point = outVertices[i];
        normal = normals[i];
        
        minimumDistance = 10e12;
        for (unsigned int t=0; t<highRestriangles.size(); t++)
        {
            Vec3 P0 = highResVertices[ highRestriangles[t][0] ].getCenter();
            Vec3 P1 = highResVertices[ highRestriangles[t][1] ].getCenter();
            Vec3 P2 = highResVertices[ highRestriangles[t][2] ].getCenter();

            Mat<3, 3, Real> M, invM;
            Vec3 P0P1 = P1-P0;
            Vec3 P0P2 = P2-P0;
            M[0][0] = -normal[0];   M[0][1] = P0P1[0];   M[0][2] = P0P2[0];
            M[1][0] = -normal[1];   M[1][1] = P0P1[1];   M[1][2] = P0P2[1];
            M[2][0] = -normal[2];   M[2][1] = P0P1[2];   M[2][2] = P0P2[2];

            Vec3 P0point = point-P0;

            // Intersection containts (t, u, v)
            Vec3 intersection;
            invM.invert(M);
            intersection = invM*P0point;

            // If intersection is within triangle
            if (intersection[1] >= 0 && intersection[2] >= 0 && intersection[1] + intersection[2] <= 1)
            {
                // Distance from the point
                distance = P0point.norm2();

                // We test the distance to only keep the closest one
                if (distance < minimumDistance)
                {
                    // We set the new minimum
                    minimumDistance = distance;

                    error = intersection[0];
                }
            }
        }

        // Stores the error
        vectorError.push_back(error);
    }
}


// --------------------------------------------------------------------------------------
// Finds the list of the closest points to a point
// --------------------------------------------------------------------------------------
template <class BaseMapping>
typename BaseMapping::Out::Real BendingPlateMechanicalMapping<BaseMapping>::FindClosestPoints(const Vec3& point1, sofa::helper::vector<unsigned int>& listClosestVertices)
{
    const InVecCoord &inVertices = *this->fromModel->getX();
    Real minimumDistance = 10e12;

    for (unsigned int v=0; v<inVertices.size(); v++)
    {
        Vec3 point2 = inVertices[v].getCenter();
        Real distance = (point2 - point1).norm2();

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
template <class BaseMapping>
typename BaseMapping::Out::Real BendingPlateMechanicalMapping<BaseMapping>::FindClosestEdges(const Vec3& point, sofa::helper::vector<unsigned int>& listClosestEdges)
{
    const InVecCoord &inVertices = *this->fromModel->getX();
    const SeqEdges &inEdges = inputTopo->getEdges();
    Real minimumDistance = 10e12;

    for (unsigned int e=0; e<inEdges.size(); e++)
    {
        Vec3 pointEdge1 = inVertices[ inEdges[e][0] ].getCenter();
        Vec3 pointEdge2 = inVertices[ inEdges[e][1] ].getCenter();

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
template <class BaseMapping>
typename BaseMapping::Out::Real BendingPlateMechanicalMapping<BaseMapping>::FindClosestTriangles(const Vec3& point, sofa::helper::vector<unsigned int>& listClosestTriangles)
{
    const InVecCoord &inVertices = *this->fromModel->getX();
    const SeqTriangles &inTriangles = inputTopo->getTriangles();
    Real minimumDistance = 10e12;

    for (unsigned int t=0; t<inTriangles.size(); t++)
    {
        Vec3 pointTriangle1 = inVertices[ inTriangles[t][0] ].getCenter();
        Vec3 pointTriangle2 = inVertices[ inTriangles[t][1] ].getCenter();
        Vec3 pointTriangle3 = inVertices[ inTriangles[t][2] ].getCenter();

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
// Barycentric coefficients of point p in triangle whose vertices are indexed by (a, b, c)
// --------------------------------------------------------------------------------------
template <class BaseMapping>
void BendingPlateMechanicalMapping<BaseMapping>::computeBaryCoefs(Vec3 &baryCoefs, const Vec3 &p, const Vec3 &a, const Vec3 &b, const Vec3 &c)
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
template <class BaseMapping>
void BendingPlateMechanicalMapping<BaseMapping>::apply( typename Out::VecCoord& out, const typename In::VecCoord& in )
{
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
//            std::cout << "vertex " << i << std::endl;

            // Gets the first triangle that the vertex belongs to
            Triangle triangle = inTriangles[ listBaseTriangles[i][0] ];

            // Gets its 3 vertices
            a = in[ triangle[0] ].getCenter();
            b = in[ triangle[1] ].getCenter();
            c = in[ triangle[2] ].getCenter();

            baryCoord = barycentricCoordinates[i][0];
//            std::cout << "baryCoord for vertex " << i << ":    " << baryCoord << std::endl;
            out[i] = a*baryCoord[0] + b*baryCoord[1] + c*baryCoord[2];

            Vec3 Uz(0.0, 0.0, 0.0);
            Vec3 w(0, 0, 0);
            for (unsigned int t=0; t<listBaseTriangles[i].size();t++)
            {
//                std::cout << "triangle " << t << std::endl;
//                std::cout << "tinfo->u " << tinfo->u << std::endl;
//                std::cout << "tinfo->u_rest " << tinfo->u_rest << std::endl;

                triangle = triangularBendingForcefield->getTopology()->getTriangle(listBaseTriangles[i][t]);
                tinfo = &triangleInf[listBaseTriangles[i][t]];

                // Local coordinates needed to compute deflection
                a = in[ triangle[0] ].getCenter();
                vertexLocal = tinfo->Qframe.rotate(out[i]-a);

                // Adds deflection
                z = tinfo->coefficients[0] + tinfo->coefficients[1]*vertexLocal[0] + tinfo->coefficients[2]*vertexLocal[1] + tinfo->coefficients[3]*vertexLocal[0]*vertexLocal[0] + tinfo->coefficients[4]*vertexLocal[0]*vertexLocal[1] + tinfo->coefficients[5]*vertexLocal[1]*vertexLocal[1] + tinfo->coefficients[6]*vertexLocal[0]*vertexLocal[0]*vertexLocal[0] + tinfo->coefficients[7]*vertexLocal[0]*vertexLocal[1]*vertexLocal[1] + tinfo->coefficients[8]*vertexLocal[1]*vertexLocal[1]*vertexLocal[1];
//                std::cout << "z = " << z << std::endl;

//                Mat<3, 3, Real > R;
//                tinfo->Qframe.toMatrix(R);
//                std::cout << "R = " << R << std::endl;

                w += tinfo->Qframe.inverseRotate(Vec3(0, 0, z));
            }

            out[i] += w/listBaseTriangles[i].size();

//            // Position in Qframe
//            out0 = tinfo->Qframe.rotate(out[i]);
//            // Computed deflection w in Qframe
//            out0[2] += w/listBaseTriangles[i].size();
//            out[i] = tinfo->Qframe.inverseRotate(out0);
        }

        triangularBendingForcefield->getTriangleInfo().endEdit();
    }

//    stop = timer.getTime();
//    std::cout << "time apply = " << stop-start << std::endl;
}


// Updates velocities of the visual mesh from mechanical vertices
template <class BaseMapping>
void BendingPlateMechanicalMapping<BaseMapping>::applyJ( typename Out::VecDeriv& out, const typename In::VecDeriv& in )
{
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
        const InVecCoord &inVertices = *this->fromModel->getX();

        // List of 'out' positions
        const OutVecCoord &outVertices = *this->toModel->getX();

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
            va_a = in[ triangle[0] ].getVOrientation();
            va_b = in[ triangle[1] ].getVOrientation();
            va_c = in[ triangle[2] ].getVOrientation();
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
//        Vec <9, Real> coeff;
        Real v_z;
        for (unsigned int i=0; i<out.size(); i++)
        {
            // Gets the first triangle that the vertex belongs to
            Triangle triangle = inTriangles[ listBaseTriangles[i][0] ];

            // Gets the linear velocities of each vertex
            v_a = in[ triangle[0] ].getVCenter();
            v_b = in[ triangle[1] ].getVCenter();
            v_c = in[ triangle[2] ].getVCenter();

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
template <class BaseMapping>
void BendingPlateMechanicalMapping<BaseMapping>::applyJT( typename In::VecDeriv& out, const typename Out::VecDeriv& in )
{
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
        const OutVecCoord &inVertices = *this->toModel->getX();

        // List of 'out' positions (mechanical points)
        const InVecCoord &outVertices = *this->fromModel->getX();

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
            out[ triangle[0] ].getVCenter() += in[i] * baryCoord[0];
            out[ triangle[1] ].getVCenter() += in[i] * baryCoord[1];
            out[ triangle[2] ].getVCenter() += in[i] * baryCoord[2];

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

                    out[ triangle[0] ].getVOrientation() += torqueA;
                    out[ triangle[1] ].getVOrientation() += torqueB;
                    out[ triangle[2] ].getVOrientation() += torqueC;
                }
            }
        }

    }


//    stop = timer.getTime();
//    std::cout << "time applyJT = " << stop-start << std::endl;
}


template <class BaseMapping>
void BendingPlateMechanicalMapping<BaseMapping>::applyJT( typename In::VecConst& /*out*/, const typename Out::VecConst& /*in*/ )
{

}


template <class BaseMapping>
void BendingPlateMechanicalMapping<BaseMapping>::drawVisual()
{
    OutVecCoord &outVertices = *this->toModel->getX();

    Real maximum = 0;
    if (measureError.getValue())
    {
        // Normalises the error
        for (unsigned int i=0; i<vectorError.size(); i++)
        {
            if (fabs(vectorError[i])>maximum)
            {
                maximum = fabs(vectorError[i]);
            }
        }
        for (unsigned int i=0; i<vectorError.size(); i++)
        {
            vectorError[i] = vectorError[i]/maximum;
        }
    }

    
    if(this->getContext()->getShowBehaviorModels())
    {
        glDisable(GL_LIGHTING);

        const SeqEdges &outEdges = triangleSubdivisionTopologicalMapping->getSubEdges();
//        const SeqEdges &outEdges = outputTopo->getEdges();
        const SeqTriangles &outTriangles = outputTopo->getTriangles();
        unsigned int index;

        glEnable(GL_DEPTH_TEST);
        glPolygonOffset(1.0, 1.0);

        if(this->getContext()->getShowWireFrame())
        {
            glEnable(GL_POLYGON_OFFSET_LINE);
            glColor4f(0.0, 0.0, 1.0, 1.0);
            glBegin(GL_LINES);
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
            glEnable(GL_POLYGON_OFFSET_FILL);
            glColor4f(0.0, 0.0, 0.0, 0.0);
            glBegin(GL_TRIANGLES);
            for (unsigned int i=0; i<outTriangles.size(); i++)
            {
                index = outTriangles[i][0];
                glTexCoord1f((float)vectorError[index]);
                glVertex3f(outVertices[index][0], outVertices[index][1], outVertices[index][2]);

                index = outTriangles[i][1];
                glTexCoord1f((float)vectorError[index]);
                glVertex3f(outVertices[index][0], outVertices[index][1], outVertices[index][2]);

                index = outTriangles[i][2];
                glTexCoord1f((float)vectorError[index]);
                glVertex3f(outVertices[index][0], outVertices[index][1], outVertices[index][2]);
            }
            glEnd();
            glDisable(GL_POLYGON_OFFSET_FILL);
        }

        // Render shells' contours (subdivision of edges)
        glColor4f(1.0, 1.0, 1.0, 1.0);
        glBegin(GL_LINES);
        for (unsigned int i=0; i<outEdges.size(); i++)
        {
            index = outEdges[i][0];
            glVertex3f(outVertices[index][0], outVertices[index][1], outVertices[index][2]);

            index = outEdges[i][1];
            glVertex3f(outVertices[index][0], outVertices[index][1], outVertices[index][2]);
        }
        glEnd();

    }
    else
    {
        glBegin(GL_LINES);
        for (unsigned int i=0; i<outVertices.size(); i++)
        {
            float alpha = vectorError[i];

            float scale = 3.0;
            alpha = scale * alpha;

            if (alpha<0.0)
            {
                glColor4f(-alpha, 0.0, 1.0+alpha, 1.0);
            }
            else
            {
                glColor4f(0.0, alpha, 1.0-alpha, 1.0);
            }

            helper::gl::glVertexT(outVertices[i]);
            helper::gl::glVertexT(outVertices[i]+normals[i]*vectorError[i]*maximum);
        }
        glEnd();
    }
}


} // namespace mapping

} // namespace component

} // namespace sofa

#endif
