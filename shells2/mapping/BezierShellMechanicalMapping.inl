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
#ifndef SOFA_COMPONENT_MAPPING_BEZIERSHELLMECHANICALMAPPING_INL
#define SOFA_COMPONENT_MAPPING_BEZIERSHELLMECHANICALMAPPING_INL

#include "BezierShellMechanicalMapping.h"
#include <sofa/component/topology/TriangleSetTopologyContainer.h>
#include <sofa/component/collision/MinProximityIntersection.h>
#include <sofa/simulation/common/Simulation.h>
//#include <sofa/helper/system/thread/CTime.h>

#include <sofa/component/forcefield/ConstantForceField.h>

#include "../forcefield/BezierShellForceField.h"
#include "../../misc/PointProjection.h"

// We have own code to check the getJ() because checkJacobian sucks (at this
// point in time).
//#define CHECK_J

namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::component::collision;
using namespace sofa::helper;


template <class TIn, class TOut>
void BezierShellMechanicalMapping<TIn, TOut>::init()
{
//    std::cout << "BezierShellMechanicalMapping::init()" << std::endl;

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
        serr << "BezierShellMechanicalMapping requires an input triangular topology" << sendl;
        return;
    }

    if (!outputTopo || (outputTopo->getNbTriangles() <= 0))
    {
        serr << "BezierShellMechanicalMapping requires an output triangular topology" << sendl;
        return;
    }

    const OutVecCoord &outVertices = *this->toModel->getX();

    barycentricCoordinates.clear();
    barycentricCoordinates.resize(outVertices.size());

    // Retrieves 'in' vertices and triangles
    const InVecCoord &inVerticesRigid = *this->fromModel->getX();

    // Conversion to Vec3Types to be able to call same methods used by Hausdorff distance
    OutVecCoord inVertices;
    for (unsigned int i=0; i<inVerticesRigid.size(); i++)
    {
        inVertices.push_back(inVerticesRigid[i].getCenter());
    }
    //const SeqEdges &inEdges = inputTopo->getEdges();
    const SeqTriangles &inTriangles = inputTopo->getTriangles();

    // Iterates over 'in' triangles
    triangleInfo.resize(inTriangles.size());
    projBaryCoords.clear();
    projN.clear();
    projElements.clear();
    for (unsigned int t=0; t<inTriangles.size(); t++) {
        TriangleInformation &tinfo = triangleInfo[t];
        tinfo.attachedPoints.clear();
    }

    PointProjection<Real> proj(*dynamic_cast<TriangleSetTopologyContainer*>(inputTopo));

    // Iterates over 'out' vertices
    for (unsigned int i=0; i<outVertices.size(); i++)
    {
        Index triangleID = 0;
        Vec3 vertexBaryCoord;

        proj.ProjectPoint(vertexBaryCoord, triangleID, outVertices[i], inVertices);

        // Mark attached point
        triangleInfo[triangleID].attachedPoints.push_back(i);

        // Add the barycentric coordinates to the list
        barycentricCoordinates[i] = vertexBaryCoord;

        projBaryCoords.push_back(vertexBaryCoord);
        ShapeFunctions N;
        bsInterpolation->computeShapeFunctions(projBaryCoords.back(), N);
        projN.push_back(N);
        projElements.push_back(triangleID);
    }

    if (measureStress.getValue())
    {
        forcefield::BezierShellForceField<TIn> *ff;
        this->getContext()->get(ff);

        if (ff)
        {
            ff->stressAtPoints(projBaryCoords, projElements);
        }
        else
        {
            serr << "Unable to find force field component! Ignoring 'measureStress' option." << sendl;
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
void BezierShellMechanicalMapping<TIn, TOut>::reinit()
{
    sout << "reinit()" << sendl;
    init();
}


// Given H,S,L in range of 0-1
// Returns a RGB colour in range of 0-255
// http://www.geekymonkey.com/Programming/CSharp/RGB2HSL_HSL2RGB.htm
template <class TIn, class TOut>
void BezierShellMechanicalMapping<TIn, TOut>::HSL2RGB(Vec3 &rgb, Real h, Real sl, Real l)
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
void BezierShellMechanicalMapping<TIn, TOut>::MeasureError()
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
typename BezierShellMechanicalMapping<TIn, TOut>::Real BezierShellMechanicalMapping<TIn, TOut>::DistanceHausdorff(BaseMeshTopology *topo1, BaseMeshTopology *topo2, helper::vector<Real> &vectorError)
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

    PointProjection<Real> proj(*dynamic_cast<TriangleSetTopologyContainer*>(inputTopo));

    // Iterates over 'in' vertices
    Real minVertex, minEdge, minTriangle, minDistance;
    Real HausdorffDistance = -1;
    for (unsigned int i=0; i<vertices1.size(); i++)
    {
        // Iterates over 'out' vertices
        minVertex = proj.FindClosestPoint(dummy, vertices1[i], vertices2);

        // Iterates over 'out' edges
        minEdge = proj.FindClosestEdge(dummy, vertices1[i], vertices2, edges2);

        // Iterates over 'out' triangles
        minTriangle = proj.FindClosestTriangle(dummy, vertices1[i], vertices2, triangles2);

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

// Updates positions of the visual mesh from mechanical vertices
template <class TIn, class TOut>
void BezierShellMechanicalMapping<TIn, TOut>::apply(const core::MechanicalParams * /*mparams*/, Data<OutVecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    helper::WriteAccessor< Data<OutVecCoord> > out = dOut;
    helper::ReadAccessor< Data<InVecCoord> > in = dIn;


    //std::cout << "---------------- Apply ----------------------------" << std::endl;

    //sofa::helper::system::thread::ctime_t start, stop;
    //sofa::helper::system::thread::CTime timer;
    //
    //start = timer.getTime();

    bsInterpolation->applyOnBTriangle(projN, projElements, out);

    //stop = timer.getTime();
    //std::cout << "time apply = " << stop-start << std::endl;
}


// Updates velocities of the visual mesh from mechanical vertices
template <class TIn, class TOut>
void BezierShellMechanicalMapping<TIn, TOut>::applyJ(const core::MechanicalParams* /*mparams*/, Data<OutVecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    helper::WriteAccessor< Data<OutVecDeriv> > out = dOut;
    helper::ReadAccessor< Data<InVecDeriv> > in = dIn;

    //std::cout << "---------------- ApplyJ ----------------------------" << std::endl;

    //sofa::helper::system::thread::ctime_t start, stop;
    //sofa::helper::system::thread::CTime timer;

    //start = timer.getTime();

    bsInterpolation->applyJOnBTriangle(projN, projElements, dIn.getValue(), out);

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

    // Dump input and output vectors {{{
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
    // }}}

    //stop = timer.getTime();
    //std::cout << "time applyJ = " << stop-start << std::endl;
}

#if 0 // TODO {{{
template <class TIn, class TOut>
const BaseMatrix* BezierShellMechanicalMapping<TIn, TOut>::getJ(const core::MechanicalParams * /*mparams*/)
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

        const OutVecCoord& out = *this->toModel->getX();
        const InVecCoord& in = *this->fromModel->getX();

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
        const InVecCoord& inVertices = *this->fromModel->getX();

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
#endif // }}}

// Updates positions of the mechanical vertices from visual    f(n-1) = JT * fn
template <class TIn, class TOut>
void BezierShellMechanicalMapping<TIn, TOut>::applyJT(const core::MechanicalParams * /*mparams*/, Data<InVecDeriv>& dOut, const Data<OutVecDeriv>& dIn)
{
    helper::WriteAccessor< Data<InVecDeriv> > out = dOut;
    helper::ReadAccessor< Data<OutVecDeriv> > in = dIn;

    //std::cout << "---------------- ApplyJT ----------------------------" << std::endl;

    //sofa::helper::system::thread::ctime_t start, stop;
    //sofa::helper::system::thread::CTime timer;

    //start = timer.getTime();

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

    bsInterpolation->applyJTOnBTriangle(projN, projElements,
        in.ref(), out);

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

    //stop = timer.getTime();
    //std::cout << "time applyJT = " << stop-start << std::endl;
}

template <class TIn, class TOut>
void BezierShellMechanicalMapping<TIn, TOut>::applyJT(const ConstraintParams* cparams, Data<InMatrixDeriv>& dOut, const Data<OutMatrixDeriv>& dIn)
{
    //std::cout << "---------------- ApplyJT (constraints) --------------" << std::endl;

    helper::WriteAccessor< Data<InMatrixDeriv> > out = dOut;
    helper::ReadAccessor< Data<OutMatrixDeriv> > in = dIn;

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

    bsInterpolation->applyJTOnBTriangle(projN, projElements, dIn.getValue(cparams), *dOut.beginEdit(cparams));
    dOut.endEdit(cparams);
}


#if 0 // TODO {{{
template <class TIn, class TOut>
void BezierShellMechanicalMapping<TIn, TOut>::draw(const core::visual::VisualParams* vparams)
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

    const OutVecCoord &outVertices = *this->toModel->getX();

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
#endif // }}}


} // namespace mapping

} // namespace component

} // namespace sofa

#endif // #ifndef SOFA_COMPONENT_MAPPING_BEZIERSHELLMECHANICALMAPPING_INL
