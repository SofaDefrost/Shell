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

#include "BezierTriangleMechanicalMapping.h"
#include <sofa/component/topology/TriangleSetTopologyContainer.h>
#include <sofa/component/collision/MinProximityIntersection.h>

#include <sofa/component/forcefield/ConstantForceField.h>

namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::component::collision;


template <class TIn, class TOut>
void BezierTriangleMechanicalMapping<TIn, TOut>::init()
{
//    std::cout << "BezierTriangleMechanicalMapping::init()" << std::endl;

    // Retrieves topology
    inputTopo = this->fromModel->getContext()->getMeshTopology();
    outputTopo = this->toModel->getContext()->getMeshTopology();

    if (inputTopo && outputTopo && inputTopo->getNbTriangles() > 0)
    {
        const OutVecCoord &outVertices = *this->toModel->getX();

        listBaseTriangles.clear();
        barycentricCoordinates.clear();
        listBaseTriangles.resize(outVertices.size());
        barycentricCoordinates.resize(outVertices.size());

        // Retrieves 'in' vertices and triangles
        const InVecCoord &inVerticesRigid = *this->fromModel->getX();
        const InVecCoord &inVerticesRigid0 = *this->fromModel->getX0();

        // Conversion to Vec3Types to be able to call same methods used by Hausdorff distance
        OutVecCoord inVertices;
        for (unsigned int i=0; i<inVerticesRigid.size(); i++)
        {
            inVertices.push_back(inVerticesRigid[i].getCenter());
        }
        const SeqEdges &inEdges = inputTopo->getEdges();
        const SeqTriangles &inTriangles = inputTopo->getTriangles();

        // Iterates over 'in' triangles
        triangleInfo.resize(inTriangles.size());
        for (unsigned int t=0; t<inTriangles.size(); t++) {

            InCoord x0[3] = {
                inVerticesRigid0[ inTriangles[t][0] ],
                inVerticesRigid0[ inTriangles[t][1] ],
                inVerticesRigid0[ inTriangles[t][2] ] };

            TriangleInformation &tinfo = triangleInfo[t];

            // get the segment positions in the reference frames of the rest-shape
            tinfo.P0_P1 = x0[0].getOrientation().inverseRotate(
                x0[1].getCenter() - x0[0].getCenter());
            tinfo.P0_P2 = x0[0].getOrientation().inverseRotate(
                x0[2].getCenter() - x0[0].getCenter());

            tinfo.P1_P2 = x0[1].getOrientation().inverseRotate(
                x0[2].getCenter() - x0[1].getCenter());
            tinfo.P1_P0 = x0[1].getOrientation().inverseRotate(
                x0[0].getCenter() - x0[1].getCenter());

            tinfo.P2_P0 = x0[2].getOrientation().inverseRotate(
                x0[0].getCenter() - x0[2].getCenter());
            tinfo.P2_P1 = x0[2].getOrientation().inverseRotate(
                x0[1].getCenter() - x0[2].getCenter());
        }

        // Iterates over 'out' vertices
        for (unsigned int i=0; i<outVertices.size(); i++)
        {
            unsigned int closestVertex, closestEdge, closestTriangle;
            Real minVertex, minEdge, minTriangle;
            int triangleID;
            Vec3 vertexBaryCoord;

            // Iterates over 'in' vertices
            minVertex = FindClosestPoint(closestVertex, outVertices[i], inVertices);

            // Iterates over 'in' edges
            minEdge = FindClosestEdge(closestEdge, outVertices[i], inVertices, inEdges);

            // Iterates over 'in' triangles
            minTriangle = FindClosestTriangle(closestTriangle, outVertices[i], inVertices, inTriangles);

            if ((minVertex <= minEdge) && (minVertex <= minTriangle))
            {
                // If it is a vertex, consider one of the triangles attached to it
                TrianglesAroundVertex trianglesAroundVertex = inputTopo->getTrianglesAroundVertex(closestVertex);
                triangleID = trianglesAroundVertex[0];
                listBaseTriangles[i] = triangleID;

                // Computes barycentric coordinates within the triangle
                computeBaryCoefs(vertexBaryCoord, outVertices[i],
                    inVertices[ inTriangles[triangleID][0] ],
                    inVertices[ inTriangles[triangleID][1] ],
                    inVertices[ inTriangles[triangleID][2] ]);

                // Adds the barycentric coordinates to the list
                barycentricCoordinates[i] = vertexBaryCoord;

            }
            else if ((minEdge <= minTriangle) && (minEdge <= minVertex))
            {
                // If it is an edge, consider one of the triangles attached to it
                TrianglesAroundEdge trianglesAroundEdge = inputTopo->getTrianglesAroundEdge(closestEdge);
                triangleID = trianglesAroundEdge[0];
                listBaseTriangles[i] = triangleID;

                // Computes barycentric coordinates within the triangle
                computeBaryCoefs(vertexBaryCoord, outVertices[i],
                    inVertices[ inTriangles[triangleID][0] ],
                    inVertices[ inTriangles[triangleID][1] ],
                    inVertices[ inTriangles[triangleID][2] ]);

                // Adds the barycentric coordinates to the list
                barycentricCoordinates[i] = vertexBaryCoord;

            }
            else
            {
                // If it is a triangle, consider it
                listBaseTriangles[i] = closestTriangle;

                // Computes barycentric coordinates within each triangles
                computeBaryCoefs(vertexBaryCoord, outVertices[i],
                    inVertices[ inTriangles[closestTriangle][0] ],
                    inVertices[ inTriangles[closestTriangle][1] ],
                    inVertices[ inTriangles[closestTriangle][2] ]);

                // Adds the barycentric coordinates to the list
                barycentricCoordinates[i] = vertexBaryCoord;
            }
        }

    }
    else
    {
        serr << "BezierTriangleMechanicalMapping requires an input triangular topology" << sendl;
        return;
    }


    // Retrieves Forcefield to compute deflection at runtime
    //triangularBendingForcefield = NULL;
    //this->getContext()->get(triangularBendingForcefield);
    //if (!triangularBendingForcefield)
    //{
    //    serr << "WARNING(BezierTriangleMechanicalMapping): triangularBendingForcefield was not found" << sendl;
    //    return;
    //}

    // Retrieves topological mapping to retrieve list of edges (to render in wireframe mode)
    triangleSubdivisionTopologicalMapping = NULL;
//    this->getContext()->get(triangleSubdivisionTopologicalMapping, nameHighTopology.getValue(), sofa::core::objectmodel::BaseContext::SearchRoot);
    this->getContext()->get(triangleSubdivisionTopologicalMapping);
    if (!triangleSubdivisionTopologicalMapping)
    {
        serr << "WARNING(BezierTriangleMechanicalMapping): triangleSubdivisionTopologicalMapping was not found" << sendl;
        return;
    }

    // Call of apply() and applyJ()
    this->Inherit::init();

    // Set each colour of each vertex to default
    const OutVecCoord &outVertices = *this->toModel->getX();
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
        const core::objectmodel::ObjectRef& refTopo = nameTargetTopology.getValue();
        topologyTarget = refTopo.getObject<TriangleSetTopologyContainer>(this->getContext());
        if (topologyTarget == NULL)
        {
            serr << "Target mesh " << nameTargetTopology.getValue() << " was not found" << sendl;
            return;
        }

        // Computes two-sided Hausdorff distance
        MeasureError();

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
        //shader.InitShaders("applications/plugins/shells/shaders/errorMap.vert", "applications/plugins/shells/shaders/errorMap.frag");

    }

}


template <class TIn, class TOut>
void BezierTriangleMechanicalMapping<TIn, TOut>::reinit()
{
    std::cout << "BezierTriangleMechanicalMapping<TIn, TOut>::reinit()" << std::endl;
    init();
}


// Given H,S,L in range of 0-1
// Returns a RGB colour in range of 0-255
// http://www.geekymonkey.com/Programming/CSharp/RGB2HSL_HSL2RGB.htm
template <class TIn, class TOut>
void BezierTriangleMechanicalMapping<TIn, TOut>::HSL2RGB(Vec3 &rgb, Real h, Real sl, Real l)
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
void BezierTriangleMechanicalMapping<TIn, TOut>::MeasureError()
{
    Real distance1;
    std::cout << "Computing Hausdorff distance high res->coarse" << std::endl;
    distance1 = DistanceHausdorff(topologyTarget, outputTopo, vectorErrorTarget);
    std::cout << "Hausdorff distance between high res mesh and coarse mesh = " << distance1 << std::endl;

    Real average = 0;
    for (unsigned int i=0; i<vectorErrorTarget.size(); i++)
    {
        average += vectorErrorTarget[i];
    }
    std::cout << "Mean Hausdorff distance = " << average/vectorErrorTarget.size() << std::endl;



    Real distance2;
    std::cout << "Computing Hausdorff distance coarse->high res" << std::endl;
    distance2 = DistanceHausdorff(outputTopo, topologyTarget, vectorErrorCoarse);
    std::cout << "Hausdorff distance between coarse mesh and high res mesh = " << distance2 << std::endl;

    average = 0;
    for (unsigned int i=0; i<vectorErrorCoarse.size(); i++)
    {
        average += vectorErrorCoarse[i];
    }
    std::cout << "Mean Hausdorff distance = " << average/vectorErrorCoarse.size() << std::endl;

}

template <class TIn, class TOut>
typename BezierTriangleMechanicalMapping<TIn, TOut>::Real BezierTriangleMechanicalMapping<TIn, TOut>::DistanceHausdorff(BaseMeshTopology *topo1, BaseMeshTopology *topo2, helper::vector<Real> &vectorError)
{
    // Mesh 1
    MechanicalState<Out>* mState1 = dynamic_cast<MechanicalState<Out>*> (topo1->getContext()->getMechanicalState());
    const OutVecCoord &vertices1 = *mState1->getX();

    // Mesh 2
    MechanicalState<Out>* mState2 = dynamic_cast<MechanicalState<Out>*> (topo2->getContext()->getMechanicalState());
    const OutVecCoord &vertices2 = *mState2->getX();
    const SeqEdges edges2 = topo2->getEdges();
    const SeqTriangles triangles2 = topo2->getTriangles();

    // The primitive is useless here
    unsigned int dummy;

    // Iterates over 'mesh1' vertices
    Real minVertex, minEdge, minTriangle, minDistance;
    Real HausdorffDistance = -1;
    for (unsigned int i=0; i<vertices1.size(); i++)
    {
        // Iterates over 'mesh2' vertices
        minVertex = FindClosestPoint(dummy, vertices1[i], vertices2);

        // Iterates over 'mesh2' edges
        minEdge = FindClosestEdge(dummy, vertices1[i], vertices2, edges2);

        // Iterates over 'mesh2' triangles
        minTriangle = FindClosestTriangle(dummy, vertices1[i], vertices2, triangles2);

        // Finds out which type of primitive is the closest
        minDistance = std::min(minVertex, std::min(minEdge, minTriangle));

        // And stores the distance between the vertex of mesh1 and mesh2
        vectorError.push_back(minDistance);

        // The maximum distance is the Hausdorff distance between mesh1 and mesh2
        if (minDistance > HausdorffDistance)
        {
            HausdorffDistance = minDistance;
        }
    }

    return HausdorffDistance;
}


// --------------------------------------------------------------------------------------
// Finds the closest point to a specified point
// --------------------------------------------------------------------------------------
template <class TIn, class TOut>
typename BezierTriangleMechanicalMapping<TIn, TOut>::Real BezierTriangleMechanicalMapping<TIn, TOut>::FindClosestPoint(unsigned int& closestVertex, const Vec3& point, const OutVecCoord &inVertices)
{
    Real minimumDistance = 10e12;
    for (unsigned int v=0; v<inVertices.size(); v++)
    {
        Real distance = (inVertices[v] - point).norm2();

        if (distance < minimumDistance)
        {
            // Store the new closest vertex
            closestVertex = v;

            // Updates the minimum's value
            minimumDistance = distance;
        }
    }

    return minimumDistance;
}


// --------------------------------------------------------------------------------------
// Finds the list of the closest edges to a point
// --------------------------------------------------------------------------------------
template <class TIn, class TOut>
typename BezierTriangleMechanicalMapping<TIn, TOut>::Real BezierTriangleMechanicalMapping<TIn, TOut>::FindClosestEdge(unsigned int& closestEdge, const Vec3& point, const OutVecCoord &inVertices, const SeqEdges &inEdges)
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
            if (distance < minimumDistance)
            {
                // Store the new closest edge
                closestEdge = e;

                // Updates the minimum's value
                minimumDistance = distance;
            }
        }
    }

    return minimumDistance;
}


// --------------------------------------------------------------------------------------
// Finds the list of the closest triangles to a point
// --------------------------------------------------------------------------------------
template <class TIn, class TOut>
typename BezierTriangleMechanicalMapping<TIn, TOut>::Real BezierTriangleMechanicalMapping<TIn, TOut>::FindClosestTriangle(unsigned int& closestTriangle, const Vec3& point, const OutVecCoord &inVertices, const SeqTriangles &inTriangles)
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
        if (alpha >= 0 && beta >= 0 && alpha + beta <= 1 )
        {
            const Vector3 PQ = AB * alpha + AC * beta - AP;

            Real distance = PQ.norm2();
            if (distance < minimumDistance)
            {
                // Store the new closest triangle
                closestTriangle = t;

                // Updates the minimum's value
                minimumDistance = distance;
            }
        }
    }

    return minimumDistance;
}


// --------------------------------------------------------------------------------------
// Barycentric coefficients of point p in triangle whose vertices are a, b and c
// --------------------------------------------------------------------------------------
template <class TIn, class TOut>
void BezierTriangleMechanicalMapping<TIn, TOut>::computeBaryCoefs(Vec3 &baryCoefs, const Vec3 &p, const Vec3 &a, const Vec3 &b, const Vec3 &c)
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
//void BezierTriangleMechanicalMapping<TIn, TOut>::apply( typename Out::VecCoord& out, const typename In::VecCoord& in )
void BezierTriangleMechanicalMapping<TIn, TOut>::apply(const core::MechanicalParams * /*mparams*/, Data<OutVecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    helper::WriteAccessor< Data<OutVecCoord> > out = dOut;
    helper::ReadAccessor< Data<InVecCoord> > in = dIn;


    //std::cout << "---------------- Apply ----------------------------" << std::endl;

    //sofa::helper::system::thread::ctime_t start, stop;
    //sofa::helper::system::thread::CTime timer;
    //
    //start = timer.getTime();

    if (!inputTopo || !outputTopo)
    {
        serr << "BezierTriangleMechanicalMapping apply() was called before init()" << sendl;
        return;
    }
    if (inputTopo->getNbTriangles() <= 0)
    {
        serr << "BezierTriangleMechanicalMapping apply() requires an input triangular topology" << sendl;
        return;
    }


    // List of in triangles
    const SeqTriangles& inTriangles = inputTopo->getTriangles();

    // Compute nodes of the Bézier triangle for each input triangle
    triangleInfo.resize(inTriangles.size());
    for (unsigned int t=0; t<inTriangles.size();t++)
    {
        //Vec3 p, cp;

        TriangleInformation &tinfo = triangleInfo[t];
        Triangle triangle = inTriangles[t];

        // Corners
        tinfo.bezierNodes[0] = in[ triangle[0] ].getCenter();
        tinfo.bezierNodes[1] = in[ triangle[1] ].getCenter();
        tinfo.bezierNodes[2] = in[ triangle[2] ].getCenter();

        Quaternion q[3] = {
            in[ triangle[0] ].getOrientation(),
            in[ triangle[1] ].getOrientation(),
            in[ triangle[2] ].getOrientation() };

#define BN(i, p, seg) do { \
    tinfo.bezierNodes[(i)] = tinfo.bezierNodes[(p)] + \
        q[(p)].rotate(tinfo.seg/3.0); \
} while (0)

        BN(3, 0, P0_P1);
        BN(4, 0, P0_P2);
        BN(5, 1, P1_P2);
        BN(6, 1, P1_P0);
        BN(7, 2, P2_P0);
        BN(8, 2, P2_P1);

#undef BN

        // Center
        tinfo.bezierNodes[9] =
            (tinfo.bezierNodes[0] + q[0].rotate( (tinfo.P0_P1 + tinfo.P0_P2)/3.0 ))/3.0 +
            (tinfo.bezierNodes[1] + q[1].rotate( (tinfo.P1_P0 + tinfo.P1_P2)/3.0 ))/3.0 +
            (tinfo.bezierNodes[2] + q[2].rotate( (tinfo.P2_P0 + tinfo.P2_P1)/3.0 ))/3.0;
    }

    for (unsigned int i=0; i<out.size(); i++)
    {
        // Gets the first triangle that the vertex belongs to
        Triangle triangle = inTriangles[ listBaseTriangles[i] ];

        Vec3 bc = barycentricCoordinates[i];

        sofa::helper::fixed_array<Vec3,10> &bn = triangleInfo[ listBaseTriangles[i] ].bezierNodes;
        // TODO: precompute the coefficients
        out[i] = bn[0] * bc[0]*bc[0]*bc[0] +
            bn[1] * bc[1]*bc[1]*bc[1] +
            bn[2] * bc[2]*bc[2]*bc[2] +
            bn[3] * 3*bc[0]*bc[0]*bc[1] +
            bn[4] * 3*bc[0]*bc[0]*bc[2] +
            bn[5] * 3*bc[1]*bc[1]*bc[2] +
            bn[6] * 3*bc[0]*bc[1]*bc[1] +
            bn[7] * 3*bc[0]*bc[2]*bc[2] +
            bn[8] * 3*bc[1]*bc[2]*bc[2] +
            bn[9] * 6*bc[0]*bc[1]*bc[2];
    }

    //{
    //    std::cout << "| In[" << in.size() << "] : ";
    //    for (unsigned int i=0; i<in.size(); i++) {
    //        if (i != 0) std::cout << ", ";
    //        std::cout << in[i];
    //    }
    //    std::cout << std::endl;

    //    helper::ReadAccessor< Data<OutVecCoord> > rout = dOut;
    //    std::cout << "| Out[" << rout.size() << "]: ";
    //    for (unsigned int i=0; i<rout.size(); i++) {
    //        if (i != 0) std::cout << ", ";
    //        std::cout << rout[i];
    //    }
    //    std::cout << std::endl;
    //}


    //stop = timer.getTime();
    //std::cout << "time apply = " << stop-start << std::endl;
}


// Updates velocities of the visual mesh from mechanical vertices
template <class TIn, class TOut>
void BezierTriangleMechanicalMapping<TIn, TOut>::applyJ(const core::MechanicalParams* /*mparams*/, Data<OutVecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    helper::WriteAccessor< Data<OutVecDeriv> > out = dOut;
    helper::ReadAccessor< Data<InVecDeriv> > in = dIn;

    //std::cout << "---------------- ApplyJ ----------------------------" << std::endl;

    //sofa::helper::system::thread::ctime_t start, stop;
    //sofa::helper::system::thread::CTime timer;
    //
    //start = timer.getTime();

    if (!inputTopo || !outputTopo)
    {
        serr << "BezierTriangleMechanicalMapping applyJ() was called before init()" << sendl;
        return;
    }
    if (inputTopo->getNbTriangles() <= 0)
    {
        serr << "BezierTriangleMechanicalMapping applyJ() requires an input triangular topology" << sendl;
        return;
    }

    // List of in triangles
    const SeqTriangles& inTriangles = inputTopo->getTriangles();
    //const helper::vector<Rigid>& inVertices = *this->fromModel->getX();
    const InVecCoord& inVertices = *this->fromModel->getX();

    // Compute nodes of the Bézier triangle for each input triangle
    for (unsigned int t=0; t<inTriangles.size();t++)
    {
        Triangle triangle = inTriangles[t];
        TriangleInformation &tinfo = triangleInfo[t];

        // Node rotations
        Quaternion q[3] = {
            inVertices[ triangle[0] ].getOrientation(),
            inVertices[ triangle[1] ].getOrientation(),
            inVertices[ triangle[2] ].getOrientation() };

        // Combined orientation quaternion of x+v
        Quaternion xvq[3] = {
            inVertices[ triangle[0] ].getOrientation()
                + Quat::createQuaterFromEuler(in[ triangle[0] ].getVOrientation()),
            inVertices[ triangle[1] ].getOrientation()
                + Quat::createQuaterFromEuler(in[ triangle[1] ].getVOrientation()),
            inVertices[ triangle[2] ].getOrientation()
                + Quat::createQuaterFromEuler(in[ triangle[2] ].getVOrientation()) };

        // Velocities in corner nodes
        tinfo.bezierNodesV[0] = in[ triangle[0] ].getVCenter();
        tinfo.bezierNodesV[1] = in[ triangle[1] ].getVCenter();
        tinfo.bezierNodesV[2] = in[ triangle[2] ].getVCenter();

#define BN(i, p, seg) do { \
    tinfo.bezierNodesV[(i)] = tinfo.bezierNodesV[(p)] \
        + xvq[(p)].rotate(tinfo.seg) \
        - q[(p)].rotate(tinfo.seg); \
} while (0)


        BN(3, 0, P0_P1);
        BN(4, 0, P0_P2);
        BN(5, 1, P1_P2);
        BN(6, 1, P1_P0);
        BN(7, 2, P2_P0);
        BN(8, 2, P2_P1);

#undef BN

        // Center
        tinfo.bezierNodesV[9] = (tinfo.bezierNodesV[0]
            + xvq[0].rotate((tinfo.P0_P1 + tinfo.P0_P2)/3.0) 
            - q[0].rotate((tinfo.P0_P1 + tinfo.P0_P2)/3.0))/3.0
        + (tinfo.bezierNodesV[1]
            + xvq[1].rotate((tinfo.P1_P2 + tinfo.P1_P0)/3.0) 
            - q[1].rotate((tinfo.P1_P2 + tinfo.P1_P0)/3.0))/3.0
        + (tinfo.bezierNodesV[2]
            + xvq[0].rotate((tinfo.P2_P0 + tinfo.P2_P1)/3.0) 
            - q[0].rotate((tinfo.P2_P0 + tinfo.P2_P1)/3.0))/3.0;
    }

    for (unsigned int i=0; i<out.size(); i++)
    {
        // Gets the first triangle that the vertex belongs to
        Triangle triangle = inTriangles[ listBaseTriangles[i] ];
        TriangleInformation &tinfo = triangleInfo[ listBaseTriangles[i] ];

        Vec3 bc = barycentricCoordinates[i];

        out[i] =
            tinfo.bezierNodesV[0] * bc[0]*bc[0]*bc[0] +
            tinfo.bezierNodesV[1] * bc[1]*bc[1]*bc[1] +
            tinfo.bezierNodesV[2] * bc[2]*bc[2]*bc[2] +
            tinfo.bezierNodesV[3] * 3*bc[0]*bc[0]*bc[1] +
            tinfo.bezierNodesV[4] * 3*bc[0]*bc[0]*bc[2] +
            tinfo.bezierNodesV[5] * 3*bc[1]*bc[1]*bc[2] +
            tinfo.bezierNodesV[6] * 3*bc[0]*bc[1]*bc[1] +
            tinfo.bezierNodesV[7] * 3*bc[0]*bc[2]*bc[2] +
            tinfo.bezierNodesV[8] * 3*bc[1]*bc[2]*bc[2] +
            tinfo.bezierNodesV[9] * 6*bc[0]*bc[1]*bc[2];
    }

    //{
    //    std::cout << "In[" << in.size() << "] : ";
    //    for (unsigned int i=0; i<in.size(); i++) {
    //        if (i != 0) std::cout << ", ";
    //        std::cout << in[i];
    //    }
    //    std::cout << std::endl;

    //    helper::ReadAccessor< Data<OutVecCoord> > rout = dOut;
    //    std::cout << "Out[" << rout.size() << "]: ";
    //    for (unsigned int i=0; i<rout.size(); i++) {
    //        if (i != 0) std::cout << ", ";
    //        std::cout << out[i];
    //    }
    //    std::cout << std::endl;
    //}

    // Approximate the velocity as two-point difference

    //Real epsilon = 1e-6;

    //Data<OutVecCoord>
    //    dOutVertices(*this->toModel->getX()),
    //    dOutVertices2(*this->toModel->getX());

    //Data<InVecCoord>
    //    dInVertices(*this->fromModel->getX()),
    //    dInVertices2(*this->fromModel->getX());

    //helper::WriteAccessor< Data<InVecCoord> > iv = dInVertices;
    //helper::WriteAccessor< Data<InVecCoord> > iv2 = dInVertices2;

    //for (unsigned int i=0; i<iv.size(); i++) {
    //    iv[i] += -in[i]*epsilon;
    //    iv2[i] += in[i]*epsilon;
    //}

    //apply(NULL, dOutVertices, dInVertices);
    //apply(NULL, dOutVertices2, dInVertices2);

    //helper::ReadAccessor< Data<OutVecCoord> > ov = dOutVertices;
    //helper::ReadAccessor< Data<OutVecCoord> > ov2 = dOutVertices2;

    //std::cout << "Dif[" << ov.size() << "]: ";
    //for (unsigned int i=0; i<out.size(); i++) {
    //    out[i] = (ov2[i] - ov[i])/(2.0*epsilon);
    //    if (i != 0) std::cout << ", ";
    //    std::cout << out[i];
    //}
    //std::cout << std::endl;

    //stop = timer.getTime();
    //std::cout << "time applyJ = " << stop-start << std::endl;
}


// Updates positions of the mechanical vertices from visual    f(n-1) = JT * fn
template <class TIn, class TOut>
void BezierTriangleMechanicalMapping<TIn, TOut>::applyJT(const core::MechanicalParams * /*mparams*/, Data<InVecDeriv>& dOut, const Data<OutVecDeriv>& dIn)
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
        serr << "BezierTriangleMechanicalMapping applyJT() was called before init()" << sendl;
        return;
    }
    if (inputTopo->getNbTriangles() <= 0)
    {
        serr << "BezierTriangleMechanicalMapping applyJT() requires an input triangular topology" << sendl;
        return;
    }

    // List of in triangles
    const SeqTriangles& inTriangles = inputTopo->getTriangles();

    Vec3 f1, f2, f3;    // resulting linear forces on corner nodes 
    Vec3 f1r, f2r, f3r; // resulting torques
    Vec3 fn;

    for (unsigned int i=0; i<in.size(); i++)
    {
        if (in[i] == Vec3(0,0,0)) continue;

        // The corresponding triangle
        int t = listBaseTriangles[i];

        Triangle triangle = inTriangles[t];
        TriangleInformation &tinfo = triangleInfo[t];

        sofa::helper::fixed_array<Vec3,10> &bn = tinfo.bezierNodesV;
        Vec3 bc = barycentricCoordinates[t];

        // Compute the influence on the corner nodes
        f1 = in[i] * bc[0]*bc[0]*bc[0];
        f2 = in[i] * bc[1]*bc[1]*bc[1];
        f3 = in[i] * bc[2]*bc[2]*bc[2];
        // TODO: what about f1r, ...?

        // Now the influence through other nodes

        fn = in[i] * 3*bc[0]*bc[0]*bc[1];
        if (fn != Vec3(0,0,0))
        {
            f1 += fn;
            f1r += cross((bn[3]-bn[0]), fn);
        }

        fn = in[i] * 3*bc[0]*bc[0]*bc[2];
        if (fn != Vec3(0,0,0))
        {
            f1 += fn;
            f1r += cross((bn[4]-bn[0]), fn);
        }

        fn = in[i] * 3*bc[1]*bc[1]*bc[2];
        if (fn != Vec3(0,0,0))
        {
            f2 += fn;
            f2r += cross((bn[5]-bn[1]), fn);
        }

        fn = in[i] * 3*bc[0]*bc[1]*bc[1];
        if (fn != Vec3(0,0,0))
        {
            f2 += fn;
            f2r += cross((bn[6]-bn[1]), fn);
        }

        fn = in[i] * 3*bc[0]*bc[2]*bc[2];
        if (fn != Vec3(0,0,0))
        {
            f3 += fn;
            f3r += cross((bn[7]-bn[2]), fn);
        }

        fn = in[i] * 3*bc[1]*bc[2]*bc[2];
        if (fn != Vec3(0,0,0))
        {
            f3 += fn;
            f3r += cross((bn[8]-bn[2]), fn);
        }

        fn = in[i] * 6*bc[0]*bc[1]*bc[2];
        if (fn != Vec3(0,0,0))
        {
            f1 += fn;
            f2 += fn;
            f3 += fn;
            f1r += cross((bn[9]-bn[0]), fn);
            f2r += cross((bn[9]-bn[1]), fn);
            f3r += cross((bn[9]-bn[2]), fn);
        }

        getVCenter(out[ triangle[0] ]) += f1;
        getVCenter(out[ triangle[1] ]) += f2;
        getVCenter(out[ triangle[2] ]) += f3;

        getVOrientation(out[ triangle[0] ]) += f1r;
        getVOrientation(out[ triangle[1] ]) += f2r;
        getVOrientation(out[ triangle[2] ]) += f3r;
    }

//    stop = timer.getTime();
//    std::cout << "time applyJT = " << stop-start << std::endl;
}


template <class TIn, class TOut>
void BezierTriangleMechanicalMapping<TIn, TOut>::applyJT(const core::ConstraintParams * /*cparams*/, Data<InMatrixDeriv>& /*dOut*/, const Data<OutMatrixDeriv>& /*dIn*/)
{

}


template <class TIn, class TOut>
void BezierTriangleMechanicalMapping<TIn, TOut>::draw()
{
    const OutVecCoord &outVertices = *this->toModel->getX();

    //bool merror = measureError.getValue();
    //if (merror)
    //    shader.TurnOn();

    if(this->getContext()->getShowVisualModels())
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

#if 0
        // Render nodes of the Bézier triangles
        glPointSize(8);
        glDisable(GL_LIGHTING);
        glBegin(GL_POINTS);
        for (unsigned int i=0; i<inputTopo->getTriangles().size(); i++)
        {
            sofa::helper::fixed_array<Vec3,10> &bn = triangleInfo[i].bezierNodes;
            for (int j=0; j<10; j++)
            {
                glColor4f(0.0, 0.5, 0.3, 1.0);
                glVertex3f(bn[j][0], bn[j][1], bn[j][2]);
            }
        } 
        glEnd();
        glPointSize(1);
        // TODO: visualise the mesh
#endif

    }

    //if (merror)
    //    shader.TurnOff();
}


} // namespace mapping

} // namespace component

} // namespace sofa

#endif
