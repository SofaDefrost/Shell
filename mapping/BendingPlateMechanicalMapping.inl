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
{
}

template <class BaseMapping>
BendingPlateMechanicalMapping<BaseMapping>::~BendingPlateMechanicalMapping()
{
}

template <class BaseMapping>
void BendingPlateMechanicalMapping<BaseMapping>::init()
{
//    std::cout << "BendingPlateMechanicalMapping<BaseMapping>::init()" << std::endl;

    this->Inherit::init();

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
        listCoeffs.resize(inputTopo->getNbTriangles());

        // Iterates over 'out' vertices
        for (unsigned int i=0; i<outVertices.size(); i++)
        {
            // Iterates over 'in' vertices
            Real minimumDistanceVertices;
            sofa::helper::vector<unsigned int> listClosestVertices;
            minimumDistanceVertices = FindClosestPoints(outVertices[i], listClosestVertices);
            
            // Iterates over 'in' edges
            Real minimumDistanceEdges;
            sofa::helper::vector<unsigned int> listClosestEdges;
            minimumDistanceEdges = FindClosestEdges(outVertices[i], listClosestEdges);

            // Iterates over 'in' triangles
            Real minimumDistanceTriangles;
            sofa::helper::vector<unsigned int> listClosestTriangles;
            minimumDistanceTriangles = FindClosestTriangles(outVertices[i], listClosestTriangles);

            // Finds out which type of primitive is the closest
            Real minimumDistance = std::min(minimumDistanceVertices, std::min(minimumDistanceEdges, minimumDistanceTriangles));

            // Retrieves 'in' vertices and triangles
            InVecCoord &inVertices = *this->fromModel->getX();
            const SeqTriangles &inTriangles = inputTopo->getTriangles();

            // Adds the list of triangles attached to the found primitives
            int caseToProcess = 0;
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
                            int triangleID = trianglesAroundVertex[t];
                            listBaseTriangles[i].push_back(triangleID);

                            // Computes barycentric coordinates within each triangles
                            Vec3 v1 = inVertices[ inTriangles[triangleID][0] ].getCenter();
                            Vec3 v2 = inVertices[ inTriangles[triangleID][1] ].getCenter();
                            Vec3 v3 = inVertices[ inTriangles[triangleID][2] ].getCenter();

                            Vec3 vertexBaryCoord;
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
                            int triangleID = trianglesAroundEdge[t];
                            listBaseTriangles[i].push_back(triangleID);

                            // Computes barycentric coordinates within each triangles
                            Vec3 v1 = inVertices[ inTriangles[triangleID][0] ].getCenter();
                            Vec3 v2 = inVertices[ inTriangles[triangleID][1] ].getCenter();
                            Vec3 v3 = inVertices[ inTriangles[triangleID][2] ].getCenter();

                            Vec3 vertexBaryCoord;
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

                        Vec3 vertexBaryCoord;
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

    // List of non-null indices within the displacemente vector u
    nonNullIndices = new int[6];
    nonNullIndices[0] = 1;  nonNullIndices[1] = 2;
    nonNullIndices[2] = 4;  nonNullIndices[3] = 5;
    nonNullIndices[4] = 7;  nonNullIndices[5] = 8;

    // Stores initial orientation for each vertex
    InVecCoord& x = *this->fromModel->getX();
    previousOrientation.resize(x.size());
    for (unsigned int i=0; i<previousOrientation.size(); i++)
    {
        previousOrientation[i] = x[i].getOrientation();
    }

    // Retrieves Forcefield to compute deflection at runtime
    triangularBendingForcefield = NULL;
    this->getContext()->get(triangularBendingForcefield);

//    ConstantForceField<InDataTypes>* ptr;
//    this->getContext()->get(ptr, 2);
//    std::cout << "ConstantForceField init() = " << ptr << std::endl;

//    std::cout << "triangularBendingForcefield init = " << triangularBendingForcefield << std::endl;

    if (!triangularBendingForcefield)
        return;

}


template <class BaseMapping>
void BendingPlateMechanicalMapping<BaseMapping>::reinit()
{
    std::cout << "BendingPlateMechanicalMapping<BaseMapping>::reinit()" << std::endl;
    init();
}

// --------------------------------------------------------------------------------------
// Finds the list of the closest points to a point
// --------------------------------------------------------------------------------------
template <class BaseMapping>
typename BaseMapping::Out::Real BendingPlateMechanicalMapping<BaseMapping>::FindClosestPoints(const Vec3& point1, sofa::helper::vector<unsigned int>& listClosestVertices)
{
    InVecCoord &inVertices = *this->fromModel->getX();
    Real minimumDistance = 10e12;

    for (unsigned int v=0; v<inVertices.size(); v++)
    {
        Vec3 point2 = inVertices[v].getCenter();
        Real distance = (point2 - point1).norm2();

        Real threshold = 0.00000001;
        if ( distance < minimumDistance - threshold )
        {
            // We deal with a new minimum, so clear the previous list
            listClosestVertices.clear();
            // Adds the new minimum
            listClosestVertices.push_back(v);
            // Updates the minimum's value
            minimumDistance = distance;
        }
        else if ( fabs(distance-minimumDistance) < threshold )
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
    InVecCoord &inVertices = *this->fromModel->getX();
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
            Real threshold = 0.000001;
            if (distance < minimumDistance - threshold)
            {
                // We deal with a new minimum, so clear the previous list
                listClosestEdges.clear();
                // Adds the new minimum
                listClosestEdges.push_back(e);
                // Updates the minimum's value
                minimumDistance = distance;
            }
            else if ( fabs(distance-minimumDistance) < threshold )
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
    InVecCoord &inVertices = *this->fromModel->getX();
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
        Real threshold = 0.000001;
        if (alpha + threshold >= 0 && beta + threshold >= 0 && alpha + beta - threshold <= 1 )
        {
            const Vector3 PQ = AB * alpha + AC * beta - AP;

            Real distance = PQ.norm2();
            if (distance < minimumDistance - threshold)
            {
                // We deal with a new minimum, so clear the previous list
                listClosestTriangles.clear();
                // Adds the new minimum
                listClosestTriangles.push_back(t);
                // Updates the minimum's value
                minimumDistance = distance;
            }
            else if ( fabs(distance-minimumDistance) < threshold )
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
        const helper::vector<TriangleInformation>& triangleInf = *(triangularBendingForcefield->getTriangleInfo().beginEdit());

        const TriangleInformation *tinfo = NULL;

        // List of in triangles
        const SeqTriangles& inTriangles = inputTopo->getTriangles();

        // Computes the coefficients ci for each triangle
        for (unsigned int t=0; t<inTriangles.size();t++)
        {
            tinfo = &triangleInf[t];
//            listCoeffs[t] = tinfo->invC.multiplyBySparseVector(tinfo->u + tinfo->u_flat, nonNullIndices, 6);
            listCoeffs[t] = tinfo->invC * (tinfo->u + tinfo->u_flat);
        }

        Vec3 a, b, c, baryCoord, vertexLocal, out0;
        Vec<9, Real> coeffs;
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
            Real w = 0;
            for (unsigned int t=0; t<listBaseTriangles[i].size();t++)
            {
                triangle = triangularBendingForcefield->getTopology()->getTriangle(listBaseTriangles[i][t]);
                tinfo = &triangleInf[listBaseTriangles[i][t]];

                // Local coordinates needed to compute deflection
                a = in[ triangle[0] ].getCenter();
                vertexLocal = tinfo->Qframe.inverseRotate(out[i]-a);

                coeffs = listCoeffs[listBaseTriangles[i][t]];

                // Adds deflection
                z = coeffs[0] + coeffs[1]*vertexLocal[0] + coeffs[2]*vertexLocal[1] + coeffs[3]*vertexLocal[0]*vertexLocal[0] + coeffs[4]*vertexLocal[0]*vertexLocal[1] + coeffs[5]*vertexLocal[1]*vertexLocal[1] + coeffs[6]*vertexLocal[0]*vertexLocal[0]*vertexLocal[0] + coeffs[7]*vertexLocal[0]*vertexLocal[1]*vertexLocal[1] + coeffs[8]*vertexLocal[1]*vertexLocal[1]*vertexLocal[1];
                w += z;

            }

            // Position in Qframe
            out0 = tinfo->Qframe.inverseRotate(out[i]);
            // Computed deflection w in Qframe
            out0[2] += w/listBaseTriangles[i].size();
            out[i] = tinfo->Qframe.rotate(out0);
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
            va_a_local = tinfo->Qframe.inverseRotate(va_a);
            va_b_local = tinfo->Qframe.inverseRotate(va_b);
            va_c_local = tinfo->Qframe.inverseRotate(va_c);
            // Fills in du/dt
            v_u.clear();
            v_u[1] = va_a_local[0];   v_u[2] = va_a_local[1];
            v_u[4] = va_b_local[0];   v_u[5] = va_b_local[1];
            v_u[7] = va_c_local[0];   v_u[8] = va_c_local[1];

//            if (t == 0)
//            {
//                std::cout << "v_u = " << v_u << std::endl;
//            }

    //        listCoeffs[t] = tinfo->invC.multiplyBySparseVector(v_u, nonNullIndices, 6);
            listCoeffs[t] = tinfo->invC * v_u;
        }


        // Iterates over out vertices to update coordinates
        Vec3 v_a, v_b, v_c, baryCoord, a, vertexLocal;
        Vec <9, Real> coeff;
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
                vertexLocal = tinfo->Qframe.inverseRotate(outVertices[i]-a);

                // Retrieves coefficients ci
                coeff = listCoeffs[listBaseTriangles[i][t]];

                // Adds deflection velocity
                v_z = coeff[0] + coeff[1]*vertexLocal[0] + coeff[2]*vertexLocal[1] + coeff[3]*vertexLocal[0]*vertexLocal[0] + coeff[4]*vertexLocal[0]*vertexLocal[1] + coeff[5]*vertexLocal[1]*vertexLocal[1] + coeff[6]*vertexLocal[0]*vertexLocal[0]*vertexLocal[0] + coeff[7]*vertexLocal[0]*vertexLocal[1]*vertexLocal[1] + coeff[8]*vertexLocal[1]*vertexLocal[1]*vertexLocal[1];
                w += v_z;

            }

            // Position in Qframe
            Vec3 out0 = tinfo->Qframe.inverseRotate(out[i]);
            // Computed deflection w in Qframe
            out0[2] += w/listBaseTriangles[i].size();
            out[i] = tinfo->Qframe.rotate(out0);
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
                inLocal = tinfo->Qframe.inverseRotate(in[i]);
                Fz = inLocal[2];

                if (Fz != 0)
                {
                    // Local coordinates needed to compute deflection
                    a = outVertices[ triangle[0] ].getCenter();
                    vertexLocal = tinfo->Qframe.inverseRotate(inVertices[i]-a);

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
                    torqueA = tinfo->Qframe.rotate(Vec3(a_u[1], a_u[2], 0));
                    torqueB = tinfo->Qframe.rotate(Vec3(a_u[4], a_u[5], 0));
                    torqueC = tinfo->Qframe.rotate(Vec3(a_u[7], a_u[8], 0));

                    out[ triangle[0] ].getVOrientation() += torqueA;
                    out[ triangle[1] ].getVOrientation() += torqueB;
                    out[ triangle[2] ].getVOrientation() += torqueC;
                }
            }
        }


        // HACK TO PREVENT ROTATION AROUND Z AXIS
    //    InVecCoord& x = *this->fromModel->getX();
    //    for (unsigned int i=0; i<out.size(); i++)
    //    {
    //        // Current orientations
    //        Quat Q = x[i].getOrientation();
    //
    //        // Previous orientations
    //        Quat Q_prev = previousOrientation[i];
    //
    //        Vec3 edgez, edgex_prev, edgex, edgey;
    //        Mat<3, 3, Real > R;
    //
    //        edgez = Q.rotate(Vec3(0.0, 0.0, 1.0));
    //        edgex_prev = Q_prev.rotate(Vec3(1.0, 0.0, 0.0));
    //        edgey = cross(edgez, edgex_prev);
    //        edgex = cross(edgey, edgez);
    //        R[0][0] = edgex[0];    R[0][1] = edgex[1];    R[0][2] = edgex[2];
    //        R[1][0] = edgey[0];    R[1][1] = edgey[1];    R[1][2] = edgey[2];
    //        R[2][0] = edgez[0];    R[2][1] = edgez[1];    R[2][2] = edgez[2];
    //
    //        Quat newOrientation;
    //        newOrientation.fromMatrix(R.transposed());
    //        x[i].getOrientation() = newOrientation;
    //
    //        // Stores orientations for next iteration
    //        previousOrientation[i] = newOrientation;
    //    }

    //    std::cout << "------------------------------------" << std::endl;
    }


//    stop = timer.getTime();
//    std::cout << "time applyJT = " << stop-start << std::endl;
}


template <class BaseMapping>
void BendingPlateMechanicalMapping<BaseMapping>::applyJT( typename In::VecConst& /*out*/, const typename Out::VecConst& /*in*/ )
{

}


} // namespace mapping

} // namespace component

} // namespace sofa

#endif
