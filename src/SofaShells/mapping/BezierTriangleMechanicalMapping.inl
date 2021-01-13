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
using namespace sofa::helper;

// Returns the skew-symetric matrix for computing a cross-product with the 
// vector @x
template <typename Real>
inline void crossMatrix(const Vec<3, Real>& x, Mat<3,3, Real>& m)
{
    m[0][0] = 0;
    m[0][1] = -x[2];
    m[0][2] = x[1];

    m[1][0] = x[2];
    m[1][1] = 0;
    m[1][2] = -x[0];

    m[2][0] = -x[1];
    m[2][1] = x[0];
    m[2][2] = 0;
}



template <class TIn, class TOut>
void BezierTriangleMechanicalMapping<TIn, TOut>::init()
{
//    std::cout << "BezierTriangleMechanicalMapping::init()" << std::endl;

    *this->f_listening.beginEdit() = true;
    this->f_listening.endEdit();

    if (this->fromModel == NULL)
    {
        serr << "Missing input Mechanical state!" << sendl;
        return;
    }

    if (this->toModel == NULL)
    {
        serr << "Missing output Mechanical state!" << sendl;
        return;
    }

    // Retrieves topology
    inputTopo = this->fromModel->getContext()->getMeshTopology();
    outputTopo = this->toModel->getContext()->getMeshTopology();

    if (!inputTopo || (inputTopo->getNbTriangles() <= 0))
    {
        serr << "BezierTriangleMechanicalMapping requires an input triangular topology" << sendl;
        return;
    }

    if (!outputTopo || (outputTopo->getNbTriangles() <= 0))
    {
        serr << "BezierTriangleMechanicalMapping requires an output triangular topology" << sendl;
        return;
    }

    // Retrieve associated Force Field (if available)
    this->getContext()->get(bezierForcefield);
    if (bezierForcefield)
    {
        sout << "TriangularBendingForcefield was found" << sendl;
        if (normals.getValue().size() != 0) {
            serr << "Ignoring normals, using configuration from force field" <<
                sendl;
        }
    }
    else
    {
        sout << "TriangularBendingForcefield was NOT found" << sendl;

        if (normals.getValue().size() == 0) {
            serr << "No normals defined, assuming flat triangles" << sendl;
        } else if (normals.getValue().size() != this->fromModel->read(sofa::core::ConstVecCoordId::position())->getValue().size()) {
            serr << "Normals count doesn't correspond with nodes count" << sendl;
            return;
        }
    }

    const OutVecCoord &outVertices = this->toModel->read(sofa::core::ConstVecCoordId::position())->getValue();

    barycentricCoordinates.clear();
    barycentricCoordinates.resize(outVertices.size());

    // Retrieves 'in' vertices and triangles
    const InVecCoord &inVerticesRigid = this->fromModel->read(sofa::core::ConstVecCoordId::position())->getValue();
    const InVecCoord &inVerticesRigid0 = this->fromModel->read(sofa::core::ConstVecCoordId::position())->getValue();

    // Conversion to Vec3Types to be able to call same methods used by Hausdorff distance
    OutVecCoord inVertices;
    for (unsigned int i=0; i<inVerticesRigid.size(); i++)
    {
        inVertices.push_back(inVerticesRigid[i].getCenter());
    }
    const SeqEdges &inEdges = inputTopo->getEdges();
    const SeqTriangles &inTriangles = inputTopo->getTriangles();

    const helper::vector<Vec3>& norms = normals.getValue();

    // Iterates over 'in' triangles
    triangleInfo.resize(inTriangles.size());
    for (unsigned int t=0; t<inTriangles.size(); t++) {
        Index a = inTriangles[t][0];
        Index b = inTriangles[t][1];
        Index c = inTriangles[t][2];

        TriangleInformation &tinfo = triangleInfo[t];

        if (!bezierForcefield)
        {

            // Compute initial bezier points at the edges 
            BezierFF::computeEdgeBezierPoints(
                a, b, c, inVerticesRigid0, norms, tinfo.bezierNodes);

            // Get the segments' position in the reference frames of the rest-shape
            tinfo.P0_P1 = inVerticesRigid0[a].getOrientation().inverseRotate( tinfo.bezierNodes[3] );
            tinfo.P0_P2 = inVerticesRigid0[a].getOrientation().inverseRotate( tinfo.bezierNodes[4] );

            tinfo.P1_P2 = inVerticesRigid0[b].getOrientation().inverseRotate( tinfo.bezierNodes[5] );
            tinfo.P1_P0 = inVerticesRigid0[b].getOrientation().inverseRotate( tinfo.bezierNodes[6] );

            tinfo.P2_P1 = inVerticesRigid0[c].getOrientation().inverseRotate( tinfo.bezierNodes[7] );
            tinfo.P2_P0 = inVerticesRigid0[c].getOrientation().inverseRotate( tinfo.bezierNodes[8] );

        }

        tinfo.attachedPoints.clear();
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

        int which = 2; /* 0 vertex, 1 edge, 2 triangle */

        if ((minVertex <= minEdge) && (minVertex <= minTriangle))
        {
            which = 0;
        }
        else if ((minEdge <= minTriangle) && (minEdge <= minVertex))
        {
            which = 1;
        }

        if (which == 0)
        {
            // If it is a vertex, consider one of the triangles attached to it
            BaseMeshTopology::TrianglesAroundVertex trianglesAroundVertex = inputTopo->getTrianglesAroundVertex(closestVertex);
            if (trianglesAroundVertex.size() <= 0)
            {
                serr << "No triangles attached to vertex " << closestVertex << sendl;
                which = (minEdge <= minTriangle) ? 1 : 2;
            }
            else
            {
                triangleID = trianglesAroundVertex[0];
                triangleInfo[triangleID].attachedPoints.push_back(i);

                // Computes barycentric coordinates within the triangle
                computeBaryCoefs(vertexBaryCoord, outVertices[i],
                    inVertices[ inTriangles[triangleID][0] ],
                    inVertices[ inTriangles[triangleID][1] ],
                    inVertices[ inTriangles[triangleID][2] ]);

                // Adds the barycentric coordinates to the list
                barycentricCoordinates[i] = vertexBaryCoord;
            }
        }

        if (which == 1)
        {
            // If it is an edge, consider one of the triangles attached to it
            BaseMeshTopology::TrianglesAroundEdge trianglesAroundEdge = inputTopo->getTrianglesAroundEdge(closestEdge);
            if (trianglesAroundEdge.size() <= 0)
            {
                serr << "No triangles attached to edge " << closestEdge << sendl;
                which = 3;
            }
            else
            {
                triangleID = trianglesAroundEdge[0];
                triangleInfo[triangleID].attachedPoints.push_back(i);

                // Computes barycentric coordinates within the triangle
                computeBaryCoefs(vertexBaryCoord, outVertices[i],
                    inVertices[ inTriangles[triangleID][0] ],
                    inVertices[ inTriangles[triangleID][1] ],
                    inVertices[ inTriangles[triangleID][2] ]);

                // Adds the barycentric coordinates to the list
                barycentricCoordinates[i] = vertexBaryCoord;
            }
        }

        if (which == 2)
        {
            // If it is a triangle, consider it
            triangleInfo[closestTriangle].attachedPoints.push_back(i);

            // Computes barycentric coordinates within each triangles
            computeBaryCoefs(vertexBaryCoord, outVertices[i],
                inVertices[ inTriangles[closestTriangle][0] ],
                inVertices[ inTriangles[closestTriangle][1] ],
                inVertices[ inTriangles[closestTriangle][2] ]);

            // Adds the barycentric coordinates to the list
            barycentricCoordinates[i] = vertexBaryCoord;
        }
    }


#if 0
    // Retrieves topological mapping to get list of edges  (for contour rendering)
    triangleSubdivisionTopologicalMapping = NULL;
    //    this->getContext()->get(triangleSubdivisionTopologicalMapping, nameHighTopology.getValue(), sofa::core::objectmodel::BaseContext::SearchRoot);
    this->getContext()->get(triangleSubdivisionTopologicalMapping);
    if (!triangleSubdivisionTopologicalMapping)
    {
        // This is not fatal
        serr << "triangleSubdivisionTopologicalMapping was not found" << sendl;
    }
#endif

    // Call of apply() and applyJ()
    this->Inherit::init();

    // Set each colour of each vertex to default
    for (unsigned int i=0; i<outVertices.size(); i++)
    {
        coloursPerVertex.push_back(Vec3(0.56, 0.14, 0.6));    // purple
    }

    // If we want to measure the error between the two meshes using Hausdorff distance
    if (measureError.getValue())
    {
        // List of colours to create a colour map
        Vec3 colour;
        Real incr = 2.0/3.0/240.0; // (2/3) is chosen stop the gradient to blue
        for (int i=0; i<240; i++)
        {
            HSL2RGB(colour, 2.0/3.0-i*incr, 0.8, 0.5);
            colourMapping.push_back(colour);
        }

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
    }
}


template <class TIn, class TOut>
void BezierTriangleMechanicalMapping<TIn, TOut>::reinit()
{
    sout << "reinit()" << sendl;
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
    sout << "Computing Hausdorff distance high res->coarse" << sendl;
    distance1 = DistanceHausdorff(targetTopology.get(), outputTopo, vectorErrorTarget);
    sout << "Hausdorff distance between high res mesh and coarse mesh = " << distance1 << sendl;

    Real average = 0;
    for (unsigned int i=0; i<vectorErrorTarget.size(); i++)
    {
        average += vectorErrorTarget[i];
    }
    sout << "Mean Hausdorff distance = " << average/vectorErrorTarget.size() << sendl;



    Real distance2;
    sout << "Computing Hausdorff distance coarse->high res" << sendl;
    distance2 = DistanceHausdorff(outputTopo, targetTopology.get(), vectorErrorCoarse);
    sout << "Hausdorff distance between coarse mesh and high res mesh = " << distance2 << sendl;

    average = 0;
    for (unsigned int i=0; i<vectorErrorCoarse.size(); i++)
    {
        average += vectorErrorCoarse[i];
    }
    sout << "Mean Hausdorff distance = " << average/vectorErrorCoarse.size() << sendl;

}

template <class TIn, class TOut>
typename BezierTriangleMechanicalMapping<TIn, TOut>::Real BezierTriangleMechanicalMapping<TIn, TOut>::DistanceHausdorff(BaseMeshTopology *topo1, BaseMeshTopology *topo2, helper::vector<Real> &vectorError)
{
    // Mesh 1
    MechanicalState<Out>* mState1 = dynamic_cast<MechanicalState<Out>*> (topo1->getContext()->getMechanicalState());
    const OutVecCoord &vertices1 = mState1->read(sofa::core::ConstVecCoordId::position())->getValue();

    // Mesh 2
    MechanicalState<Out>* mState2 = dynamic_cast<MechanicalState<Out>*> (topo2->getContext()->getMechanicalState());
    const OutVecCoord &vertices2 = mState2->read(sofa::core::ConstVecCoordId::position())->getValue();
    const SeqEdges edges2 = topo2->getEdges();
    const SeqTriangles triangles2 = topo2->getTriangles();

    // The primitive is useless here
    unsigned int dummy;

    // Iterates over 'in' vertices
    Real minVertex, minEdge, minTriangle, minDistance;
    Real HausdorffDistance = -1;
    for (unsigned int i=0; i<vertices1.size(); i++)
    {
        // Iterates over 'out' vertices
        minVertex = FindClosestPoint(dummy, vertices1[i], vertices2);

        // Iterates over 'out' edges
        minEdge = FindClosestEdge(dummy, vertices1[i], vertices2, edges2);

        // Iterates over 'out' triangles
        minTriangle = FindClosestTriangle(dummy, vertices1[i], vertices2, triangles2);

        // Finds out which type of primitive is the closest
        minDistance = std::min(minVertex, std::min(minEdge, minTriangle));

        // And stores the distance for the vertex
        vectorError.push_back(minDistance);

        // The maximum distance is the Hausdorff distance
        if (minDistance > HausdorffDistance)
        {
            HausdorffDistance = minDistance;
        }
    }

    return HausdorffDistance;
}


// --------------------------------------------------------------------------------------
// Finds the closest point to a point
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
// Finds the closest edge to a point
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
// Finds the closest triangle to a point
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

        Vec3 bary;
        computeBaryCoefs(bary, point,
            pointTriangle1, pointTriangle2, pointTriangle3, false);
        if ((bary[0] < 0.0) || (bary[1] < 0.0) || (bary[2] < 0.0) ||
            (rabs(1.0 - (bary[0] + bary[1] + bary[2])) > 1e-10)) {
            // Point projected onto the plane of the triangle lies outside
            // of the triangle. Some vertex or edge will be more
            // appropriate.
            continue;
        }

        Vector3 N = cross(AB, AC);
        //Real distance = N*point - N*pointTriangle1;
        Real distance = N*(point - pointTriangle1);
        distance = distance*distance / N.norm2();

        if (distance < minimumDistance)
        {
            // Store the new closest triangle
            closestTriangle = t;

            // Updates the minimum's value
            minimumDistance = distance;
        }
    }

    return minimumDistance;
}


// --------------------------------------------------------------------------------------
// Barycentric coefficients of point p in triangle whose vertices are a, b and c.
// If bConstraint is true constraint the coordinates to lie inside the triangle.
// --------------------------------------------------------------------------------------
template <class TIn, class TOut>
void BezierTriangleMechanicalMapping<TIn, TOut>::computeBaryCoefs(Vec3 &baryCoefs,
    const Vec3 &p, const Vec3 &a, const Vec3 &b, const Vec3 &c, bool bConstraint)
{
    const double ZERO = 1e-20;

    Vec3 M = (Vec3) (b-a).cross(c-a);
    double norm2_M = M*(M);

    double coef_a, coef_b, coef_c;

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
        if (bConstraint)
        {
            // Do some magic to constraint the coordinates inside the triangle
            // the requirements are:
            //    coef_a, coef_b, coef_c ≥ 0
            //    coef_a + coef_b + coef_c = 1
            if (coef_a < 0.0) coef_a = 0.0;
            if (coef_b < 0.0) coef_b = 0.0;
            coef_c = 1.0 - (coef_a + coef_b);
            if (coef_c < 0.0)
            {
                // We have to be carefull so as not to overshoot some other
                // coefficient
                if (coef_a < -coef_c/2.0) {
                    coef_c += coef_a;
                    coef_b += coef_c;
                    coef_a = 0.0;
                } else if (coef_b < -coef_c/2.0) {
                    coef_c += coef_b;
                    coef_a += coef_c;
                    coef_b = 0.0;
                } else {
                    coef_a += coef_c/2.0;
                    coef_b += coef_c/2.0;
                }

                coef_c = 0.0;
            }
        }
        else
        {
            coef_c = N*((a-p).cross(b-p));
        }
    }

    baryCoefs[0] = coef_a;
    baryCoefs[1] = coef_b;
    baryCoefs[2] = coef_c;
}

// Updates positions of the visual mesh from mechanical vertices
template <class TIn, class TOut>
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
        serr << "apply() was called before init()" << sendl;
        return;
    }
    if (inputTopo->getNbTriangles() <= 0)
    {
        serr << "apply() requires an input triangular topology" << sendl;
        return;
    }

    // List of in triangles
    const SeqTriangles& inTriangles = inputTopo->getTriangles();

    for (unsigned int t=0; t<inTriangles.size();t++)
    {

        TriangleInformation &tinfo = triangleInfo[t];
        Triangle triangle = inTriangles[t];

        if (bezierForcefield)
        {
            // Use data from the force field
            const typename BezierFF::TriangleInformation fftinfo = 
                bezierForcefield->getTriangleInfo(t);
            tinfo.bezierNodes = fftinfo.bezierNodes; 

            tinfo.P0_P1 = fftinfo.P0_P1_inFrame0;
            tinfo.P0_P2 = fftinfo.P0_P2_inFrame0;
            tinfo.P1_P2 = fftinfo.P1_P2_inFrame1;
            tinfo.P1_P0 = fftinfo.P1_P0_inFrame1;
            tinfo.P2_P0 = fftinfo.P2_P0_inFrame2;
            tinfo.P2_P1 = fftinfo.P2_P1_inFrame2;

        }
        else
        {
            // Compute nodes of the Bézier triangle

            // Corners
            tinfo.bezierNodes[0] = in[ triangle[0] ].getCenter();
            tinfo.bezierNodes[1] = in[ triangle[1] ].getCenter();
            tinfo.bezierNodes[2] = in[ triangle[2] ].getCenter();

#ifdef ROTQ
            Quaternion q[3] = {
                in[ triangle[0] ].getOrientation(),
                in[ triangle[1] ].getOrientation(),
                in[ triangle[2] ].getOrientation() };

#define BN(i, p, seg) do { \
    tinfo.bezierNodes[(i)] = tinfo.bezierNodes[(p)] + \
        q[(p)].rotate(tinfo.seg); \
} while (0)

#else
            // Rotation matrices at corner nodes
            Mat33 R[3];
            in[ triangle[0] ].getOrientation().toMatrix(R[0]);
            in[ triangle[1] ].getOrientation().toMatrix(R[1]);
            in[ triangle[2] ].getOrientation().toMatrix(R[2]);

#define BN(i, p, seg) do { \
    tinfo.bezierNodes[(i)] = tinfo.bezierNodes[(p)] + \
        R[(p)] * tinfo.seg; \
} while (0)

#endif
            BN(3, 0, P0_P1);
            BN(4, 0, P0_P2);
            BN(5, 1, P1_P2);
            BN(6, 1, P1_P0);
            BN(7, 2, P2_P0);
            BN(8, 2, P2_P1);

#undef BN

#ifdef ROTQ
            // Center
            tinfo.bezierNodes[9] =
                (tinfo.bezierNodes[0] + q[0].rotate( tinfo.P0_P1 + tinfo.P0_P2 ))/3.0 +
                (tinfo.bezierNodes[1] + q[1].rotate( tinfo.P1_P0 + tinfo.P1_P2 ))/3.0 +
                (tinfo.bezierNodes[2] + q[2].rotate( tinfo.P2_P0 + tinfo.P2_P1 ))/3.0;
#else
            // Center
            tinfo.bezierNodes[9] =
                (tinfo.bezierNodes[0] + R[0]*( tinfo.P0_P1 + tinfo.P0_P2 ))/3.0 +
                (tinfo.bezierNodes[1] + R[1]*( tinfo.P1_P0 + tinfo.P1_P2 ))/3.0 +
                (tinfo.bezierNodes[2] + R[2]*( tinfo.P2_P0 + tinfo.P2_P1 ))/3.0;
#endif

        }

        // Go through the attached points
        for (unsigned int i=0; i<tinfo.attachedPoints.size(); i++)
        {
            Index pt = tinfo.attachedPoints[i] ;
            Vec3 bc = barycentricCoordinates[pt];

            // TODO: precompute the coefficients
            out[pt] = tinfo.bezierNodes[0] * bc[0]*bc[0]*bc[0] +
                tinfo.bezierNodes[1] * bc[1]*bc[1]*bc[1] +
                tinfo.bezierNodes[2] * bc[2]*bc[2]*bc[2] +
                tinfo.bezierNodes[3] * 3*bc[0]*bc[0]*bc[1] +
                tinfo.bezierNodes[4] * 3*bc[0]*bc[0]*bc[2] +
                tinfo.bezierNodes[5] * 3*bc[1]*bc[1]*bc[2] +
                tinfo.bezierNodes[6] * 3*bc[0]*bc[1]*bc[1] +
                tinfo.bezierNodes[7] * 3*bc[0]*bc[2]*bc[2] +
                tinfo.bezierNodes[8] * 3*bc[1]*bc[2]*bc[2] +
                tinfo.bezierNodes[9] * 6*bc[0]*bc[1]*bc[2];
        }

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

//    sofa::helper::system::thread::ctime_t start, stop;
//    sofa::helper::system::thread::CTime timer;
//
//    start = timer.getTime();

    if (!inputTopo || !outputTopo)
    {
        serr << "applyJ() was called before init()" << sendl;
        return;
    }
    if (inputTopo->getNbTriangles() <= 0)
    {
        serr << "applyJ() requires an input triangular topology" << sendl;
        return;
    }

    // List of in triangles
    const SeqTriangles& inTriangles = inputTopo->getTriangles();
    //const helper::vector<Rigid>& inVertices = *this->fromModel->getX();
    const InVecCoord& inVertices = this->fromModel->read(sofa::core::ConstVecCoordId::position())->getValue();

    // Compute nodes of the Bézier triangle for each input triangle
    for (unsigned int t=0; t<inTriangles.size();t++)
    {
        Triangle triangle = inTriangles[t];
        TriangleInformation &tinfo = triangleInfo[t];

        // Velocities in corner nodes
        tinfo.bezierNodesV[0] = in[ triangle[0] ].getVCenter();
        tinfo.bezierNodesV[1] = in[ triangle[1] ].getVCenter();
        tinfo.bezierNodesV[2] = in[ triangle[2] ].getVCenter();

        // Angular velocities in cross-product matrix
        Mat33 Omega0, Omega1, Omega2;
        crossMatrix<Real>(in[ triangle[0] ].getVOrientation(), Omega0);
        crossMatrix<Real>(in[ triangle[1] ].getVOrientation(), Omega1);
        crossMatrix<Real>(in[ triangle[2] ].getVOrientation(), Omega2);

        Mat33 dR0, dR1, dR2;

        // Rotation matrices at corner nodes
        inVertices[ triangle[0] ].getOrientation().toMatrix(dR0);
        inVertices[ triangle[1] ].getOrientation().toMatrix(dR1);
        inVertices[ triangle[2] ].getOrientation().toMatrix(dR2);

        // Derivatives of the rotation matrix
        dR0 = Omega0*dR0;
        dR1 = Omega1*dR1;
        dR2 = Omega2*dR2;

        // Velocities at other nodes
        tinfo.bezierNodesV[3] = tinfo.bezierNodesV[0] + dR0*tinfo.P0_P1;
        tinfo.bezierNodesV[4] = tinfo.bezierNodesV[0] + dR0*tinfo.P0_P2;

        tinfo.bezierNodesV[5] = tinfo.bezierNodesV[1] + dR1*tinfo.P1_P2;
        tinfo.bezierNodesV[6] = tinfo.bezierNodesV[1] + dR1*tinfo.P1_P0;

        tinfo.bezierNodesV[7] = tinfo.bezierNodesV[2] + dR2*tinfo.P2_P0;
        tinfo.bezierNodesV[8] = tinfo.bezierNodesV[2] + dR2*tinfo.P2_P1;

        tinfo.bezierNodesV[9] = (
            (tinfo.bezierNodesV[0] + dR0*(tinfo.P0_P1 + tinfo.P0_P2)) +
            (tinfo.bezierNodesV[1] + dR1*(tinfo.P1_P0 + tinfo.P1_P2)) +
            (tinfo.bezierNodesV[2] + dR2*(tinfo.P2_P0 + tinfo.P2_P1))
        )/3.0;

        for (unsigned int i=0; i<tinfo.attachedPoints.size(); i++)
        {

            Index pt = tinfo.attachedPoints[i];
            Vec3 bc = barycentricCoordinates[pt];

            out[pt] =
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
    }

    // The following code compares the result with results obtained using
    // getJ() because checkJacobian sucks (at this point in time).
#ifdef CHECK_J
    const sofa::defaulttype::BaseMatrix* J = getJ(NULL);
    if (J != NULL) {
        Real* in_alloc = NULL;
        Real* out_alloc = NULL;

        // Prepare in vector
        in_alloc = new Real[in.size()*NIn];
        for (unsigned int i=0;i<in.size();++i)
            for (int j=0;j<NIn;++j)
                in_alloc[i*NIn+j] = (Real)in[i][j];

        // Multiply
        out_alloc = new Real[out.size()*NOut];
        J->opMulV(out_alloc, in_alloc);

        // Compare results
        Real amax = 0; Index maxi=0;
        //std::cout << "Delta with getJ():";
        Real dif;
        for (unsigned int i=0;i<out.size();++i)
            for (int j=0;j<NOut;++j) {
                dif = out_alloc[i*NOut+j] - out[i][j];
                //std::cout << " " << dif;
                //out[i][j] = out_alloc[i*NOut+j];
                //if (rabs(dif[k]) < 1e-5) dif[k] = 0;
                if (rabs(dif) > amax) { amax = rabs(dif); maxi = i; }
            }
        //std::cout << "\n";
        if (amax > 1e-9)
            std::cout << "check J: amax=" << amax << " i=" << maxi << " phi=" <<
                barycentricCoordinates[maxi][0] << "/" <<
                barycentricCoordinates[maxi][1] << "/" <<
                barycentricCoordinates[maxi][2] << "\n";

        // Cleanup
        delete[] in_alloc;
        delete[] out_alloc;
    }
#endif

    // Dump input and output vectors
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

    //Real amax = 0; Index maxi=0;
    ////std::cout << "Dif[" << ov.size() << "]: ";
    //for (unsigned int i=0; i<out.size(); i++) {
    //    OutCoord tmp = (ov2[i] - ov[i])/(2.0*epsilon);
    //    OutCoord dif = tmp - out[i];
    //    out[i] = tmp;
    //    for (int k=0; k<3; k++) { //if (rabs(dif[k]) < 1e-5) dif[k] = 0;
    //        if (rabs(dif[k]) > amax) { amax = rabs(dif[k]); maxi = i; } }
    //    //if (i != 0) std::cout << ", ";
    //    //std::cout << dif;
    //    //std::cout << out[i];
    //}
    ////std::cout << std::endl;
    //if (amax > 1e-9)
    //std::cout << "amax=" << amax << " i=" << maxi << " phi=" <<
    //    barycentricCoordinates[maxi][0] << "/" <<
    //    barycentricCoordinates[maxi][1] << "/" <<
    //    barycentricCoordinates[maxi][2] << "\n";

//    stop = timer.getTime();
//    std::cout << "time applyJ = " << stop-start << std::endl;
}

template <class TIn, class TOut>
const BaseMatrix* BezierTriangleMechanicalMapping<TIn, TOut>::getJ(const core::MechanicalParams * /*mparams*/)
{
    //std::cout << "---------------- getJ ----------------------------" << std::endl;

    if (matrixJ.get() == NULL || updateJ)
    {
        if (!inputTopo || !outputTopo)
        {
            serr << "getJ() was called before init()" << sendl;
            return NULL;
        }
        if (inputTopo->getNbTriangles() <= 0)
        {
            serr << "getJ() requires an input triangular topology" << sendl;
            return NULL;
        }

        const OutVecCoord& out = this->toModel->read(sofa::core::ConstVecCoordId::position())->getValue();
        const InVecCoord& in = this->fromModel->read(sofa::core::ConstVecCoordId::position())->getValue();

        // Initialize the matrix
        if (matrixJ.get() == 0 ||
            matrixJ->rowBSize() != out.size() ||
            matrixJ->colBSize() != in.size())
        {
            matrixJ.reset(new MatrixType(
                    out.size() * NOut,
                    in.size() * NIn));
        }
        else
        {
            matrixJ->clear();
        }

        Mat33 I(Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1));

        const SeqTriangles& inTriangles = inputTopo->getTriangles();
        const InVecCoord& inVertices = this->fromModel->read(sofa::core::ConstVecCoordId::position())->getValue();

        // Go through all input triangles
        for (unsigned int t=0; t<inTriangles.size();t++)
        {
            Triangle triangle = inTriangles[t];
            TriangleInformation &tinfo = triangleInfo[t];

            // Rotation matrices at corner nodes
            Mat33 R[3];
            inVertices[ triangle[0] ].getOrientation().toMatrix(R[0]);
            inVertices[ triangle[1] ].getOrientation().toMatrix(R[1]);
            inVertices[ triangle[2] ].getOrientation().toMatrix(R[2]);

            // Cross matrices for rotated control nodes
            Mat33 Ap1[3]; // nodes 3, 5, 7
            crossMatrix(R[0]*tinfo.P0_P1, Ap1[0]);
            crossMatrix(R[1]*tinfo.P1_P2, Ap1[1]);
            crossMatrix(R[2]*tinfo.P2_P0, Ap1[2]);

            Mat33 Ap2[3]; // nodes 4, 6, 8
            crossMatrix(R[0]*tinfo.P0_P2, Ap2[0]);
            crossMatrix(R[1]*tinfo.P1_P0, Ap2[1]);
            crossMatrix(R[2]*tinfo.P2_P1, Ap2[2]);

            // Transpose the matrices to change the order of arguments in the
            // cross product
            for (int i=0; i<3; i++) {
                Ap1[i].transpose();
                Ap2[i].transpose();
            }

            // Go through all attached nodes
            for (unsigned int i=0; i<tinfo.attachedPoints.size(); i++)
            {
                Index pt = tinfo.attachedPoints[i];
                Vec3 bc = barycentricCoordinates[pt];

                // Go through the three nodes of a trinagle and consider their
                // respective influences
                for (int k=0; k<3; k++)
                {
                    MBloc& block = *matrixJ->wbloc(pt, triangle[k], true);
                    Mat33 trans, ang;
                    block.getsub(0, 0, trans);  // Translational DOFS
                    block.getsub(0, 3, ang);    // Angular DOFS

                    // Corner node
                    trans += I * bc[k]*bc[k]*bc[k];

                    //// 3 / 5 / 7
                    int l = (k + 1) % 3;
                    trans += I * 3*bc[k]*bc[k]*bc[l];
                    ang += Ap1[k] * 3*bc[k]*bc[k]*bc[l];

                    // 4 / 6 / 8
                    l = (k + 2) % 3;
                    trans += I * 3*bc[k]*bc[k]*bc[l];
                    ang += Ap2[k] * 3*bc[k]*bc[k]*bc[l];

                    // Central node
                    trans += I * 2*bc[0]*bc[1]*bc[2]; // <-- 6/3 = 2
                    ang += (Ap1[k] + Ap2[k]) * 2*bc[0]*bc[1]*bc[2]; // <-- 6/3 = 2

                    block.setsub(0, 0, trans);
                    block.setsub(0, 3, ang);
                }
            }
        }

    } // if (matrixJ.get() == NULL || updateJ)

    return matrixJ.get();
}

// Updates positions of the mechanical vertices from visual    f(n-1) = JT * fn
template <class TIn, class TOut>
void BezierTriangleMechanicalMapping<TIn, TOut>::applyJT(const core::MechanicalParams * /*mparams*/, Data<InVecDeriv>& dOut, const Data<OutVecDeriv>& dIn)
{
    helper::WriteAccessor< Data<InVecDeriv> > out = dOut;
    helper::ReadAccessor< Data<OutVecDeriv> > in = dIn;

    //std::cout << "---------------- ApplyJT ----------------------------" << std::endl;

//    sofa::helper::system::thread::ctime_t start, stop;
//    sofa::helper::system::thread::CTime timer;
//
//    start = timer.getTime();

    if (!inputTopo || !outputTopo)
    {
        serr << "applyJT() was called before init()" << sendl;
        return;
    }
    if (inputTopo->getNbTriangles() <= 0)
    {
        serr << "applyJT() requires an input triangular topology" << sendl;
        return;
    }

#ifdef CHECK_J
    const sofa::defaulttype::BaseMatrix* J = getJ(NULL);
    Real* in_alloc = NULL;
    Real* out_alloc = NULL;
    if (J != NULL) {
        // Prepare in vector
        in_alloc = new Real[in.size()*NOut];
        for (unsigned int i=0;i<in.size();++i)
            for (int j=0;j<NOut;++j)
                in_alloc[i*NOut+j] = (Real)in[i][j];

        // Prepare out vector
        out_alloc = new Real[out.size()*NIn];
        for (unsigned int i=0;i<out.size();++i)
            for (int j=0;j<NIn;++j)
                out_alloc[i*NIn+j] = (Real)out[i][j];
    }
#endif

    // List of in triangles
    const SeqTriangles& inTriangles = inputTopo->getTriangles();
    const InVecCoord& inVertices = this->fromModel->read(sofa::core::ConstVecCoordId::position())->getValue();

    // Compute nodes of the Bézier triangle for each input triangle
    for (unsigned int t=0; t<inTriangles.size();t++)
    {
        Triangle triangle = inTriangles[t];
        TriangleInformation &tinfo = triangleInfo[t];

        // Rotation matrices at corner nodes
        Mat33 R[3];
        inVertices[ triangle[0] ].getOrientation().toMatrix(R[0]);
        inVertices[ triangle[1] ].getOrientation().toMatrix(R[1]);
        inVertices[ triangle[2] ].getOrientation().toMatrix(R[2]);

        for (unsigned int i=0; i<tinfo.attachedPoints.size(); i++)
        {
            Vec3 f1, f2, f3;    // resulting linear forces on corner nodes 
            Vec3 f1r, f2r, f3r; // resulting torques
            Vec3 fn;

            Index pt = tinfo.attachedPoints[i];
            Vec3 bc = barycentricCoordinates[pt];

            if (in[pt] == Vec3(0,0,0)) continue;

            // Compute the influence on the corner nodes
            f1 = in[pt] * (bc[0]*bc[0]*bc[0]);
            f2 = in[pt] * (bc[1]*bc[1]*bc[1]);
            f3 = in[pt] * (bc[2]*bc[2]*bc[2]);

            // Now the influence through other nodes

            sofa::helper::fixed_array<Vec3,10> &bn = tinfo.bezierNodes;

            fn = in[pt] * (3*bc[0]*bc[0]*bc[1]);
            if (fn != Vec3(0,0,0))
            {
                f1 += fn;
                f1r += cross((bn[3]-bn[0]), fn);
            }

            fn = in[pt] * (3*bc[0]*bc[0]*bc[2]);
            if (fn != Vec3(0,0,0))
            {
                f1 += fn;
                f1r += cross((bn[4]-bn[0]), fn);
            }

            fn = in[pt] * (3*bc[1]*bc[1]*bc[2]);
            if (fn != Vec3(0,0,0))
            {
                f2 += fn;
                f2r += cross((bn[5]-bn[1]), fn);
            }

            fn = in[pt] * (3*bc[0]*bc[1]*bc[1]);
            if (fn != Vec3(0,0,0))
            {
                f2 += fn;
                f2r += cross((bn[6]-bn[1]), fn);
            }

            fn = in[pt] * (3*bc[0]*bc[2]*bc[2]);
            if (fn != Vec3(0,0,0))
            {
                f3 += fn;
                f3r += cross((bn[7]-bn[2]), fn);
            }

            fn = in[pt] * (3*bc[1]*bc[2]*bc[2]);
            if (fn != Vec3(0,0,0))
            {
                f3 += fn;
                f3r += cross((bn[8]-bn[2]), fn);
            }

            fn = in[pt] * (2*bc[0]*bc[1]*bc[2]);
            if (fn != Vec3(0,0,0))
            {
                f1 += fn;
                f2 += fn;
                f3 += fn;
                f1r += cross(R[0]*(tinfo.P0_P1 + tinfo.P0_P2), fn);
                f2r += cross(R[1]*(tinfo.P1_P2 + tinfo.P1_P0), fn);
                f3r += cross(R[2]*(tinfo.P2_P0 + tinfo.P2_P1), fn);
            }

            getVCenter(out[ triangle[0] ]) += f1;
            getVCenter(out[ triangle[1] ]) += f2;
            getVCenter(out[ triangle[2] ]) += f3;

            getVOrientation(out[ triangle[0] ]) += f1r;
            getVOrientation(out[ triangle[1] ]) += f2r;
            getVOrientation(out[ triangle[2] ]) += f3r;
        }
    }

    // The following code compares the result with results obtained using
    // getJ() because checkJacobian sucks (at this point in time).
#ifdef CHECK_J
    if (J != NULL) {
        J->opPMulTV(out_alloc, in_alloc);

        // Compare results
        Real amax = 0; Index maxi=0;
        //std::cout << "Delta with getJT():";
        Real dif;
        for (unsigned int i=0;i<out.size();++i)
            for (int j=0;j<NIn;++j) {
                dif = out_alloc[i*NIn+j] - out[i][j];
                //std::cout << " " << dif;
                out[i][j] = out_alloc[i*NIn+j];
                if (rabs(dif) > amax) { amax = rabs(dif); maxi = i; }
            }
        //std::cout << "\n";
        if (amax > 1e-9)
            std::cout << "check JT: amax=" << amax << " i=" << maxi << " phi=" <<
                barycentricCoordinates[maxi][0] << "/" <<
                barycentricCoordinates[maxi][1] << "/" <<
                barycentricCoordinates[maxi][2] << 
                //" val= " << out_alloc[maxi*NIn+0] << " " <<
                //out_alloc[maxi*NIn+1] << " " <<
                //out_alloc[maxi*NIn+2] << " " <<
                //out_alloc[maxi*NIn+3] << " " <<
                //out_alloc[maxi*NIn+4] << " " <<
                //out_alloc[maxi*NIn+5] << " / " <<
                //out[maxi] <<
                "\n";

        // Cleanup
        delete[] in_alloc;
        delete[] out_alloc;
    }
#endif

//    stop = timer.getTime();
//    std::cout << "time applyJT = " << stop-start << std::endl;
}


template <class TIn, class TOut>
void BezierTriangleMechanicalMapping<TIn, TOut>::applyJT(const core::ConstraintParams * /*cparams*/, Data<InMatrixDeriv>& /*dOut*/, const Data<OutMatrixDeriv>& /*dIn*/)
{
    //serr << "applyJT(const core::ConstraintParams*, Data<InMatrixDeriv>&, const Data<OutMatrixDeriv>&) NOT implemented" << sendl;
}


template <class TIn, class TOut>
void BezierTriangleMechanicalMapping<TIn, TOut>::draw(const core::visual::VisualParams* vparams)
{
    if (!inputTopo || !outputTopo)
    {
        serr << "draw() was called before init()" << sendl;
        return;
    }
    if (inputTopo->getNbTriangles() <= 0)
    {
        serr << "draw() requires an input triangular topology" << sendl;
        return;
    }
    if (outputTopo->getNbTriangles() <= 0)
    {
        serr << "draw() requires an output triangular topology" << sendl;
        return;
    }

    const OutVecCoord &outVertices = this->toModel->read(sofa::core::ConstVecCoordId::position())->getValue();

    if(vparams->displayFlags().getShowVisualModels())
    {
        glDisable(GL_LIGHTING);

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
                glVertex3f(outVertices[index][0],
                    outVertices[index][1],
                    outVertices[index][2]);

                index = outTriangles[i][1];
                glVertex3f(outVertices[index][0],
                    outVertices[index][1],
                    outVertices[index][2]);

                index = outTriangles[i][2];
                glVertex3f(outVertices[index][0],
                    outVertices[index][1],
                    outVertices[index][2]);
            }
            glEnd();
            glDisable(GL_POLYGON_OFFSET_LINE);
        }
        else
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glEnable(GL_POLYGON_OFFSET_FILL);
            glBegin(GL_TRIANGLES);
            for (unsigned int i=0; i<outTriangles.size(); i++)
            {
                index = outTriangles[i][0];
                glColor4f(coloursPerVertex[index][0],
                    coloursPerVertex[index][1],
                    coloursPerVertex[index][2],
                    1.0);
                glVertex3f(outVertices[index][0],
                    outVertices[index][1],
                    outVertices[index][2]);

                index = outTriangles[i][1];
                glColor4f(coloursPerVertex[index][0],
                    coloursPerVertex[index][1],
                    coloursPerVertex[index][2],
                    1.0);
                glVertex3f(outVertices[index][0],
                    outVertices[index][1],
                    outVertices[index][2]);

                index = outTriangles[i][2];
                glColor4f(coloursPerVertex[index][0],
                    coloursPerVertex[index][1],
                    coloursPerVertex[index][2],
                    1.0);
                glVertex3f(outVertices[index][0],
                    outVertices[index][1],
                    outVertices[index][2]);
            }
            glEnd();
            glDisable(GL_POLYGON_OFFSET_FILL);
        }

#if 0
        // Render shells' contours (subdivision of edges)
        if (triangleSubdivisionTopologicalMapping)
        {
            const SeqEdges &outEdges = triangleSubdivisionTopologicalMapping->getSubEdges();
            //const SeqEdges &outEdges = outputTopo->getEdges();
        glColor4f(1.0, 1.0, 1.0, 1.0);
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
        }
#endif
    }

    if(vparams->displayFlags().getShowMechanicalMappings())
    {
        // Render nodes of the Bézier triangles
        glPointSize(8);
        glDisable(GL_LIGHTING);
        glBegin(GL_POINTS);
        //for (unsigned int i=0; i<inputTopo->getTriangles().size(); i++)
        unsigned int i=3;
        {
            sofa::helper::fixed_array<Vec3,10> &bn = triangleInfo[i].bezierNodes;
            for (int j=0; j<10; j++)
            {
                //glColor4f(0.0, 0.5, 0.3, 1.0);
                glColor4f(0.5, 1.0, 0.5, 1.0);
                glVertex3f(bn[j][0], bn[j][1], bn[j][2]);
            }
        } 
        glEnd();
        glPointSize(1);

    }
        // TODO: visualise the mesh
}


} // namespace mapping

} // namespace component

} // namespace sofa

#endif
