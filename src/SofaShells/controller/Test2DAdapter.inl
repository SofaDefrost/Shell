//
// Component for dynamic remeshing.
//
// Few notes for all the brave adventurers.
// * There is a thin line between operations performed in rest shape and in
//   deformed shape. Specificaly the element geometry is analysed in rest
//   shape, but point tracking for cutting support is (mostly) done in deformed
//   shape.
// * The point tracking is not ideal because the optimizion is slow (read:
//   makes small updates). Idealy we would want to run several optimization
//   steps per simulation step. Unfortunately the work related to relocation of
//   points (updating mechanical state, propagation of topololgical changes,
//   ...) is really big so we cannot do that.
//
// TODO:
// - Bugs
//    * triangles get inverted
//    * attaching to fixed points and then relocating them
//
// - reattach points in N1-ring during cutting to avoid degeneration
// - refactor into several components
//
// Edge cases that are not handled:
//   - cut is too short (doesn't span two edges)
//   - cutting near the edge (handling of boundary/fixed nodes)
//
// TODO parametrization
//   + updates on point relocation
//       * done by propagation of topological changes
//   + factor out optimization related stuff
//   + modify optimization process to use parametrized 2D surface 
//   - check if cutting still works
//   - fix parametrization
//   - compute metric tensors
//   - add Bézier surfaces
//   - fix edge swapping

#ifndef SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_INL
#define SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_INL

#include "../initPluginShells.h"

#include <map>
#include <float.h>

#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/core/visual/VisualParams.h>  
#include <sofa/helper/rmath.h>
#include <sofa/helper/SimpleTimer.h>
#include <SofaBaseTopology/TopologyData.inl>

#include <SofaMeshCollision/TriangleModel.h>

#include "../misc/PointProjection.h"

#include "Test2DAdapter.h"

#define OTHER(x, a, b) ((x == a) ? b : a)

// Return non-zero if triangle with points (a,b,c) is defined in
// counter-clockwise direction. 
// NOTE: Constrained to 2D!
#define CCW(a,b,c) (\
    cross(Vec2(b[0]-a[0], b[1]-a[1]), \
        Vec2(c[0]-a[0], c[1]-a[1])) > 1e-15)

namespace sofa
{

namespace component
{

namespace controller
{

//sofa::helper::TSimpleTimer<1,10> mytimer;

template<class DataTypes>
void Test2DAdapter<DataTypes>::PointInfoHandler:: applyCreateFunction(
    unsigned int pointIndex,
    PointInformation &/*pInfo*/,
    const sofa::core::topology::Point &point,
    const sofa::helper::vector< unsigned int > &ancestors,
    const sofa::helper::vector< double > &coeffs)
{
    adapter->m_toUpdate[pointIndex] = true;
    adapter->m_surf.pointAdd(pointIndex, point, ancestors, coeffs);
}



template<class DataTypes>
void Test2DAdapter<DataTypes>::PointInfoHandler::applyDestroyFunction(unsigned int pointIndex, PointInformation& /*pInfo*/)
{
    //std::cout << "pt " << __FUNCTION__ << pointIndex << std::endl;
    adapter->m_toUpdate.erase(pointIndex);
    adapter->m_surf.pointRemove(pointIndex);
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::PointInfoHandler::swap(unsigned int i1,
    unsigned int i2)
{
    //std::cout << "pt " << __FUNCTION__ << " " << i1 << " " << i2 << std::endl;
    Inherited::swap(i1, i2);

    // m_pointId may change
    if (adapter->m_pointId == i1) {
        adapter->m_pointId = i2;
    } else if (adapter->m_pointId == i2) {
        adapter->m_pointId = i1;
    }

    // Update indices in m_toUpdate
    if (adapter->m_toUpdate[i1] && !adapter->m_toUpdate[i2]) {
        adapter->m_toUpdate.erase(i1);
        adapter->m_toUpdate[i2] = true;
    } else if (adapter->m_toUpdate[i2] && !adapter->m_toUpdate[i1]) {
        adapter->m_toUpdate.erase(i2);
        adapter->m_toUpdate[i1] = true;
    }

    adapter->m_surf.pointSwap(i1, i2);
}




template<class DataTypes>
void Test2DAdapter<DataTypes>::TriangleInfoHandler::applyCreateFunction(
    unsigned int triangleIndex, TriangleInformation &tInfo,
    const sofa::core::topology::Triangle& elem,
    const sofa::helper::vector< unsigned int > &/*ancestors*/,
    const sofa::helper::vector< double > &/*coeffs*/)
{
    //std::cout << "tri " << __FUNCTION__ << triangleIndex << " (" << elem << ") [" << adapter->m_container->getTriangle(triangleIndex) << "]" << std::endl;

    // Check if vertices are on the boundary
    // NOTE: We cannot do any complicated checking here because the topology
    //       may be inconsistent.
    const Triangle& t = adapter->m_container->getTriangle(triangleIndex);
    adapter->m_toUpdate[t[0]] = true;
    adapter->m_toUpdate[t[1]] = true;
    adapter->m_toUpdate[t[2]] = true;

    // Compute initial normal
    adapter->computeTriangleNormal(elem, adapter->m_state->read(sofa::core::VecCoordId::restPosition())->getValue(), tInfo.normal);
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::TriangleInfoHandler::applyDestroyFunction(unsigned int triangleIndex, TriangleInformation &/*tInfo*/)
{
    //std::cout << "tri " << __FUNCTION__ << triangleIndex << std::endl;
    // Check if vertices are on the boundary
    // NOTE: We cannot do any complicated checking here because the topology
    //       may be inconsistent.
    const Triangle& t = adapter->m_container->getTriangle(triangleIndex);
    adapter->m_toUpdate[t[0]] = true;
    adapter->m_toUpdate[t[1]] = true;
    adapter->m_toUpdate[t[2]] = true;
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::TriangleInfoHandler::swap(unsigned int i1,
    unsigned int i2)
{
    //std::cout << "tri " << __FUNCTION__ << " " << i1 << " " << i2 << std::endl;
    Inherited::swap(i1, i2);

    if (adapter->m_pointTriId == i1) {
        adapter->m_pointTriId = i2;
    } else if (adapter->m_pointTriId == i2) {
        adapter->m_pointTriId = i1;
    }
}

#if 0
template<class DataTypes>
void Test2DAdapter<DataTypes>::TriangleInfoHandler::addOnMovedPosition(const sofa::helper::vector<unsigned int> &indexList,
    const sofa::helper::vector< topology::Triangle > & elems)
{
    std::cout << "tri " << __FUNCTION__ << " " << indexList << std::endl;
    Inherited::addOnMovedPosition(indexList, elems);
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::TriangleInfoHandler::removeOnMovedPosition(const sofa::helper::vector<unsigned int> &indices)
{
    std::cout << "tri " << __FUNCTION__ << " " << indices << std::endl;
    Inherited::removeOnMovedPosition(indices);
}


#endif








template<class DataTypes>
Test2DAdapter<DataTypes>::Test2DAdapter()
: m_sigma(initData(&m_sigma, (Real)0.01, "sigma", "Minimal increase in functional to accept the change"))
, m_functionals(initData(&m_functionals, "functionals", "Current values of the functional for each triangle"))
//, m_mappedState(initLink("mappedState", "Points to project onto the topology."))
, m_projectedPoints(initData(&m_projectedPoints, "projectedPoints", "Points to project onto the topology."))
, m_interpolationIndices(initData(&m_interpolationIndices, "interpolationIndices",
        "Interpolation indices for projected points."))
, m_interpolationValues(initData(&m_interpolationValues, "interpolationValues",
        "Interpolation values for projected points."))
, stepCounter(0)
, m_precision(1e-8)
, m_pointId(InvalidID)
, m_pointTriId(InvalidID)
, m_opt(this, m_surf)
, pointInfo(initData(&pointInfo, "pointInfo", "Internal point data"))
, triInfo(initData(&triInfo, "triInfo", "Internal triangle data"))
{
    pointHandler = new PointInfoHandler(this, &pointInfo);
    triHandler = new TriangleInfoHandler(this, &triInfo);
}

template<class DataTypes>
Test2DAdapter<DataTypes>::~Test2DAdapter()
{
    if(pointHandler) delete pointHandler;
    if(triHandler) delete triHandler;
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::init()
{
    m_state = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes>*> (this->getContext()->getMechanicalState());
    if (!m_state) {
        serr << "Unable to find MechanicalState" << sendl;
        return;
    }

    this->getContext()->get(m_container);
    if (m_container == NULL) {
        serr << "Unable to find triangular topology" << sendl;
        return;
    }

    this->getContext()->get(m_modifier);
    if (m_modifier == NULL) {
        serr << "Unable to find TriangleSetTopologyModifier" << sendl;
        return;
    }

    this->getContext()->get(m_algoGeom);
    if (m_algoGeom == NULL) {
        serr << "Unable to find TriangleSetGeometryAlgorithms" << sendl;
        return;
    }

    this->getContext()->get(m_algoTopo);
    if (m_algoTopo == NULL) {
        serr << "Unable to find TriangleSetTopologyAlgorithms" << sendl;
        return;
    }

    pointInfo.createTopologicalEngine(m_container, pointHandler);
    pointInfo.registerTopologicalData();

    triInfo.createTopologicalEngine(m_container, triHandler);
    triInfo.registerTopologicalData();


    reinit();
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::reinit()
{
    *this->f_listening.beginEdit() = true;
    this->f_listening.endEdit();

    if ((m_sigma.getValue() < (Real)0.0) || (m_sigma.getValue() > (Real)1.0)) {
        serr << "The value of sigma must be between 0 and 1." << sendl;
        *m_sigma.beginEdit() = 0.01;
        m_sigma.endEdit();
    }
    
    m_functionals.beginEdit()->resize(m_container->getNbTriangles(), Real(0.0));
    m_functionals.endEdit();

    helper::vector<PointInformation>& pts = *pointInfo.beginEdit();
    pts.resize(m_container->getNbPoints());
    //for (int i=0; i<m_container->getNbPoints(); i++)
    //{
    //    pointHandler->applyCreateFunction(
    //        i,
    //        pts[i],
    //        m_container->getPoint(i),
    //        (const sofa::helper::vector< unsigned int > )0,
    //        (const sofa::helper::vector< double >)0);
    //}
    pointInfo.endEdit();


    helper::vector<TriangleInformation>& tris = *triInfo.beginEdit();
    tris.resize(m_container->getNbTriangles());
    for (int i=0; i<m_container->getNbTriangles(); i++)
    {
        triHandler->applyCreateFunction(
            i,
            tris[i],
            m_container->getTriangle(i),
            (const sofa::helper::vector< unsigned int > )0,
            (const sofa::helper::vector< double >)0);
    }
    triInfo.endEdit();

    recheckBoundary();

    projectionInit();

    stepCounter = 0;

    m_surf.init(m_container, m_state->read(sofa::core::VecCoordId::restPosition())->getValue());
}


template<class DataTypes>
void Test2DAdapter<DataTypes>::onEndAnimationStep(const double /*dt*/)
{
    //std::cout << "CPU step\n";

    if ((m_container == NULL) || (m_state == NULL))
        return;

    stepCounter++;

    // Update boundary vertices
    recheckBoundary();

    //if (stepCounter < f_interval.getValue())
    //    return; // Stay still for a few steps

    //stepCounter = 0;

    Data<VecCoord>* datax = m_state->write(sofa::core::VecCoordId::restPosition());
    //VecCoord& x = *datax->beginEdit();
    // WARNING: Notice that we're working on a copy that is NOT updated by
    //          external changes!
    VecCoord x = datax->getValue();
    //
    const VecCoord& x0 = datax->getValue();

    Index nTriangles = m_container->getNbTriangles();
    if (nTriangles == 0)
        return;

    // Update projection of tracked point in rest shape
    if (m_pointId != InvalidID) {
        Triangle tri = m_container->getTriangle(m_pointTriId);
        helper::vector< double > bary = m_algoGeom->compute3PointsBarycoefs(
            m_point, tri[0], tri[1], tri[2], false);
        m_pointRest =
            x0[ tri[0] ] * bary[0] +
            x0[ tri[1] ] * bary[1] +
            x0[ tri[2] ] * bary[2];
    }


    //sofa::helper::system::thread::ctime_t start, stop;
    //sofa::helper::system::thread::CTime timer;
    //start = timer.getTime();

    //mytimer.start("Init");
    sofa::helper::vector<Real> &functionals = *m_functionals.beginEdit();

    functionals.resize(nTriangles);
    m_opt.initValues(functionals, m_container);

    //std::cout << "m: " << functionals << "\n";
    //mytimer.stop();

    ngamma = 0;
    sumgamma = maxgamma = 0.0;
    mingamma = 1.0;

    Real maxdelta=0.0;
    unsigned int moved=0;
    //mytimer.step("Optimize");
    for (Index i=0; i<x.size(); i++) {
        if (pointInfo.getValue()[i].isFixed()) {
            //std::cout << "skipping fixed node " << i << "\n";
            continue;
        }

        //mytimer.start("Optimize");
        Vec2 newPos;
        Index tId=InvalidID;
        if (m_opt.smooth(i, newPos, tId,
                functionals, m_sigma.getValue(), m_precision)) {
            if (tId == InvalidID) {
                std::cout << "BUG!\n";
            }
            // Move the point
            Vec3 xold = x[i];
            Vec3 xnew = m_surf.getPointPosition(newPos, tId, x0);
            //std::cout << "    moving " << xold << " -- " << xnew << "\n";
            //mytimer.step("Relocate");
            relocatePoint(i, xnew);

            // Update boundary vertices
            //mytimer.step("Recheck boundary");
            recheckBoundary();

            moved++;
            Real delta = (x[i] - xold).norm2();
            if (delta > maxdelta) {
                maxdelta = delta;
            }

            // Update projection of tracked point in rest shape
            if (m_pointId == i) {
                Triangle tri = m_container->getTriangle(m_pointTriId);
                helper::vector< double > bary = m_algoGeom->compute3PointsBarycoefs(
                    m_point, tri[0], tri[1], tri[2], false);
                m_pointRest =
                    x0[ tri[0] ] * bary[0] +
                    x0[ tri[1] ] * bary[1] +
                    x0[ tri[2] ] * bary[2];
            }
        }
        //mytimer.stop();
    }
    //mytimer.step("Eval");

    //stop = timer.getTime();
    //std::cout << "---------- CPU time = " << stop-start << "\n";

    // Evaluate improvement
    Real sum=0.0, sum2=0.0, min = DBL_MAX;
    Index minTriID = InvalidID;
    for (Index i=0; i < nTriangles; i++) {
        if (functionals[i] < min) {
            min = functionals[i];
            minTriID = i;
        }
        sum += functionals[i];
        sum2 += functionals[i] * functionals[i];
    }
    sum /= nTriangles;
    sum2 = helper::rsqrt(sum2/nTriangles);

    // Try swapping edge for the worst triangle
    //mytimer.step("Swap");
    //swapEdge(minTriID);
    //NOTE: we do some work twice, this can be optimized.
    //for (Index i=0; i<m_container->getNumberOfTriangles(); i++) {
    //    swapEdge(i);
    //}
    //mytimer.stop();


    //std::cout << stepCounter << "] moved " << moved << " points, max delta=" << helper::rsqrt(maxdelta)
    //    << " gamma min/avg/max: " << mingamma << "/" << sumgamma/ngamma
    //    << "/" << maxgamma

    //    << " Quality min/avg/RMS: " << min << "/" << sum << "/" << sum2

    //    << "\n";

    datax->endEdit();

    m_functionals.endEdit();

    //// Write metrics to file
    //std::ofstream of("/tmp/metrics.csv", std::ios::app);
    //of << "geom," << stepCounter;
    //for (Index i=0; i < m_functionals.getValue().size(); i++) {
    //    of << "," << m_functionals.getValue()[i];
    //}
    //of << "\n";
    //of.close();
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::onKeyPressedEvent(core::objectmodel::KeypressedEvent *key)
{
    if (key->getKey() == 'M') { // Ctrl+M
        std::cout << "Writing metrics to file /tmp/metrics.csv\n";
        // Write metrics to file
        std::ofstream of("/tmp/metrics.csv", std::ios::app);
        of << "geom," << stepCounter;
        for (Index i=0; i < m_functionals.getValue().size(); i++) {
            of << "," << m_functionals.getValue()[i];
        }
        of << "\n";
        of.close();
    }
}

#if 0 // TODO
template<class DataTypes>
void Test2DAdapter<DataTypes>::swapEdge(Index triID)
{
    Index swapTri = InvalidID;
    Real swapVal1, swapVal2, swapMin = DBL_MIN;
    Triangle swapT1, swapT2;
    // Keep old value
    Real oldValue = m_functionals.getValue()[triID];

    const VecCoord& x0 = m_state->read(
        sofa::core::ConstVecCoordId::restPosition())->getValue();
    const VecCoord& x = m_state->read(
        sofa::core::ConstVecCoordId::position())->getValue();
    const Triangle &t1 = m_container->getTriangle(triID);

    const EdgesInTriangle &elist = m_container->getEdgesInTriangle(triID);

    int orientation0 = CCW(x0[t1[0]], x0[t1[1]], x0[t1[2]]);
    int orientation = CCW(x[t1[0]], x[t1[1]], x[t1[2]]);

    for (int ie=0; ie < 3; ie++) {

        const TrianglesAroundEdge &tlist =
            m_container->getTrianglesAroundEdge(elist[ie]);
        if (tlist.size() != 2)
            continue; // No other triangle or non-manifold edge

        bool bProtected = false;
        for (unsigned int i=0; i<m_protectedEdges.size(); i++) {
            if (m_protectedEdges[i] == elist[ie]) {
                bProtected = true;
                break;
            }
        }

        if (bProtected)
            continue; // Do not touch!

        // TODO: is the pair t1-t2 convex?
        //       -- this check may not be necessary, swapping should lead to
        //          inverted triangle
        //       -- then again, pair may be convex in rest shape but
        //          concave in deformed shape

        Index triID2 = OTHER(triID, tlist[0], tlist[1]);
        //Real oldValue2 = m_functionals.getValue()[triID2];
        const Edge &e = m_container->getEdge(elist[ie]);
        const Triangle &t2 = m_container->getTriangle(triID2);

        if (orientation0 != CCW(x0[t2[0]], x0[t2[1]], x0[t2[2]])) {
            std::cout << "Unable to swap edge for triangles with different "
                "orientation in rest shape.\n";
            continue;
        } else if (orientation != CCW(x[t2[0]], x[t2[1]], x[t2[2]])) {
            std::cout << "Unable to swap edge for triangles with different "
                "orientation.\n";
            continue;
        }

        // Detect indices of what
        Index ne1 = InvalidID, ne2 = InvalidID; // Points not on the shared edge
        int swap1 = -1, swap2 = -1;
        for (int i=0; i<3; i++) {
            if ((t1[i] != e[0]) && (t1[i] != e[1])) {
                ne1 = t1[i];
            } else if (t1[i] == e[0]) {
                swap1 = i;
            }
            if ((t2[i] != e[0]) && (t2[i] != e[1])) {
                ne2 = t2[i];
            } else if (t2[i] == e[1]) {
                swap2 = i;
            }
        }


        // Swap edge shared between the triangles
        Triangle nt1 = t1, nt2 = t2;
        nt1[swap1] = ne2;
        nt2[swap2] = ne1;

        // Check if this helps
        Real newValue = funcTriangle(nt1, x0,
            triInfo.getValue()[triID].normal);
        Real newValue2 = funcTriangle(nt2, x0,
            triInfo.getValue()[triID2].normal);

        // Is there improvement?
        // NOTE: We assume that oldValue < oldValue2.
        if ((newValue > oldValue) && (newValue2 > oldValue) &&
            (newValue > swapMin) && (newValue2 > swapMin)) {
            // Yes, this helps!
            swapTri = triID2;
            swapMin = (newValue < newValue2) ? newValue : newValue2;
            swapVal1 = newValue;
            swapVal2 = newValue2;
            swapT1 = nt1;
            swapT2 = nt2;
        }
    }

    // Perform the topological operation
    if (swapTri != InvalidID) {
        //std::cout << "Swapping\n";

        if ((triID == m_pointTriId) || (swapTri == m_pointTriId)) {
            std::cout << "Cannot swap tracked triangle!\n";
            // TODO: we need to project tracked position to find out which
            // of the triangles contains it.
        }

        //if (CCW(x0[swapT1[0]], x0[swapT1[1]], x0[swapT1[2]]) !=
        //    CCW(x0[swapT2[0]], x0[swapT2[1]], x0[swapT2[2]])) {
        //    std::cout << "WHOOPS! Inverting triangle! (rest)\n";
        //} else if (CCW(x[swapT1[0]], x[swapT1[1]], x[swapT1[2]]) !=
        //    CCW(x[swapT2[0]], x[swapT2[1]], x[swapT2[2]])) {
        //    std::cout << "WHOOPS! Inverting triangle!\n";
        //}

        sofa::helper::vector< unsigned int > del;
        del.push_back(triID);
        del.push_back(swapTri);
        m_modifier->removeTriangles(del, true, false);

        sofa::helper::vector<Triangle> add;
        add.push_back(swapT1);
        add.push_back(swapT2);
        m_modifier->addTriangles(add);
    }

    // TODO: projection
}
#endif

template<class DataTypes>
void Test2DAdapter<DataTypes>::computeTriangleNormal(const Triangle &t, const VecCoord &x, Vec3 &normal) const
{
    Vec3 A, B;
    A = x[ t[1] ] - x[ t[0] ];
    B = x[ t[2] ] - x[ t[0] ];

    Real An = A.norm(), Bn = B.norm();
    if (An < 1e-20 || Bn < 1e-20) {
        serr << "Found degenerated triangle: "
            << x[ t[0] ] << " / " << x[ t[1] ] << " / " << x[ t[2] ]
            << " :: " << An << ", " << Bn << sendl;

        normal = Vec3(0,0,0);
        return;
    }

    A /= An;
    B /= Bn;
    normal = cross(A, B);
    normal.normalize();
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::recheckBoundary()
{
    helper::vector<PointInformation>& pts = *pointInfo.beginEdit();
    for (std::map<Index,bool>::const_iterator i=m_toUpdate.begin();
        i != m_toUpdate.end(); i++) {
        pts[i->first].type = detectNodeType(i->first, pts[i->first].boundary);
    }
    pointInfo.endEdit();
    m_toUpdate.clear();
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::relocatePoint(Index pt, Coord target,
    Index hint, bool bInRest)
{
    if (m_modifier == NULL ||
        m_algoGeom == NULL ||
        m_container == NULL ||
        m_state == NULL)
        return;

    //mytimer.step("Relocate");

    Index tId = hint;

    const VecCoord& x = m_state->read(
        bInRest
        ? sofa::core::ConstVecCoordId::restPosition()
        : sofa::core::ConstVecCoordId::position()
        )->getValue();

    if ((x[pt] - target).norm() < m_precision) return; // Nothing to do

    if (tId == InvalidID) {
        tId = m_algoGeom->getTriangleInDirection(pt, target - x[pt]);
    }

    if (tId == InvalidID) {
        // Try once more and more carefully
        TrianglesAroundVertex N1 = m_container->getTrianglesAroundVertex(pt);
        for (Index it=0; it<N1.size(); it++) {
            Index t;
            if (m_algoGeom->isPointInsideTriangle(N1[it], false, target, t, bInRest)) {
                Triangle tri = m_container->getTriangle(N1[it]);
                helper::vector<double> bary =
                    m_algoGeom->compute3PointsBarycoefs(
                        target, tri[0], tri[1], tri[2], bInRest);
                bool bGood = true;
                for (int bc=0; bc<3; bc++) {
                    if (bary[bc] < -1e-15) {
                        bGood = false;
                        break;
                    }
                }

                if (!bGood) continue;

                tId = N1[it];
                break;
            }
        }
    }

    if (tId == InvalidID) {
        // We screwed up something. Probably a triangle got inverted accidentaly.
        serr << "Unexpected triangle Id -1! Cannot move point "
            << pt << ", marking point as fixed." << sendl;

        (*pointInfo.beginEdit())[pt].forceFixed = true;
        pointInfo.endEdit();

        //EdgesAroundVertex N1e = m_container->getEdgesAroundVertex(pt);
        //for (Index ip=0; ip<N1e.size(); ip++) {
        //    Edge e = m_container->getEdge(N1e[ip]);
        //    std::cout << m_container->getPointDataArray().getValue()[ e[0] ] << " " << m_container->getPointDataArray().getValue()[ e[1] ] << "\n";
        //}

        return;
    }

    Triangle tri = m_container->getTriangle(tId);

    helper::vector <unsigned int> move_ids;
    helper::vector< helper::vector< unsigned int > > move_ancestors;
    helper::vector< helper::vector< double > > move_coefs;

    move_ids.push_back(pt);

    move_ancestors.resize(1);
    move_ancestors.back().push_back(tri[0]);
    move_ancestors.back().push_back(tri[1]);
    move_ancestors.back().push_back(tri[2]);

    move_coefs.push_back(
        m_algoGeom->compute3PointsBarycoefs(target,
            tri[0], tri[1], tri[2], bInRest));

    // Do the real work
    m_modifier->movePointsProcess(move_ids, move_ancestors, move_coefs);
    m_modifier->notifyEndingEvent();
    m_modifier->propagateTopologicalChanges();

    //mytimer.step("Projection");
    projectionUpdate(pt);

    // Check
    //const VecCoord& xnew = m_state->read(
    //    sofa::core::ConstVecCoordId::restPosition())->getValue();
    //std::cout << "requested " << target << "\tgot " << xnew[pt] << "\tdelta"
    //    << (target - xnew[pt]).norm() << "\n";
}

template<class DataTypes>
typename Test2DAdapter<DataTypes>::PointInformation::NodeType Test2DAdapter<DataTypes>::detectNodeType(Index pt, Vec3 &boundaryDirection)
{
    if (m_container == NULL)
        return PointInformation::NORMAL;

    const VecCoord& x0 = m_state->read(sofa::core::ConstVecCoordId::restPosition())->getValue();

    VecVec3 dirlist; // Directions of boundary edges

    //std::cout << "checking " << pt << "\n";
    EdgesAroundVertex N1e = m_container->getEdgesAroundVertex(pt);
    for (Index ie=0; ie<N1e.size(); ie++) {
        TrianglesAroundEdge tlist = m_container->getTrianglesAroundEdge(N1e[ie]);
        unsigned int count = tlist.size();
        //std::cout << "  " << count << " (" << tlist << ")\n";
        if (count < 1 || count > 2) {
            // What a strange topology! Better not touch ...
            return PointInformation::FIXED;
        } else if (count == 1) {
            Edge e = m_container->getEdge(N1e[ie]);
            Vec3 dir;
            if (e[0] == pt) {
                dir = x0[e[1]] - x0[e[0]];
            } else {
                dir = x0[e[0]] - x0[e[1]];
            }
            dir.normalize();
            dirlist.push_back(dir);
        }
    }

    if (dirlist.size() == 0) {
        // No boundary edges
        return PointInformation::NORMAL;
    } else if (dirlist.size() != 2) {
        // What a strange topology! Better not touch ...
        return PointInformation::FIXED;
    }

    if ((1.0 + dirlist[0] * dirlist[1]) > 1e-15) {
        // Not on a line
        return PointInformation::FIXED;
    }

    boundaryDirection = dirlist[0];

    return PointInformation::BOUNDARY;
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::projectionInit()
{
    if (!m_container) return;

    const VecCoord& x0= m_state->read(
        sofa::core::ConstVecCoordId::restPosition())->getValue();

    const VecCoord& xProj = m_projectedPoints.getValue();
    unsigned int nVertices = xProj.size();

    sofa::helper::vector<sofa::helper::vector< unsigned int > > &indices =
        *m_interpolationIndices.beginEdit();
    sofa::helper::vector<sofa::helper::vector< Real > > &values =
        *m_interpolationValues.beginEdit();

    indices.resize(nVertices);
    values.resize(nVertices);

    helper::vector<TriangleInformation>& tris = *triInfo.beginEdit();
    // Clear list of previously attached points.
    for (Index i=0; i<m_container->getNumberOfTriangles(); i++) {
        tris[i].attachedPoints.clear();
    }

    PointProjection<Real> proj(*m_container);

    for (Index i=0; i<nVertices; i++)
    {
        Index triangleID;
        Vec3 vertexBaryCoord;

        proj.ProjectPoint(vertexBaryCoord, triangleID, xProj[i], x0);
        if (triangleID == InvalidID) {
            serr << "Failed to project point " << i << "!" << sendl;
            break;
        }

        // Mark attached point
        tris[triangleID].attachedPoints.push_back(i);

        // Add node indices to the list
        Triangle t = m_container->getTriangle(triangleID);
        indices[i].clear();
        indices[i].push_back(t[0]);
        indices[i].push_back(t[1]);
        indices[i].push_back(t[2]);

        // Add the barycentric coordinates to the list
        values[i].clear();
        values[i].push_back(vertexBaryCoord[0]);
        values[i].push_back(vertexBaryCoord[1]);
        values[i].push_back(vertexBaryCoord[2]);
    }

    m_interpolationIndices.endEdit();
    m_interpolationValues.endEdit();
    triInfo.endEdit();
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::projectionUpdate(Index pt)
{
    if (!m_container) return;

    PointProjection<Real> proj(*m_container);

    const VecCoord& x0= m_state->read(
        sofa::core::ConstVecCoordId::restPosition())->getValue();

    const VecCoord& xProj = m_projectedPoints.getValue();

    sofa::helper::vector<sofa::helper::vector< unsigned int > > &indices =
        *m_interpolationIndices.beginEdit();
    sofa::helper::vector<sofa::helper::vector< Real > > &values =
        *m_interpolationValues.beginEdit();

    helper::vector<TriangleInformation>& tris = *triInfo.beginEdit();

    TrianglesAroundVertex N1 = m_container->getTrianglesAroundVertex(pt);

    // List of newly attached points for each triangle
    // We keep it separated to avoid double work.
    std::map<Index, helper::vector<Index> > newAttached;

    // Recompute all barycentric coordinates
    for (unsigned int it=0; it<N1.size(); it++) {
        Index t = N1[it];
        Triangle tri = m_container->getTriangle(t);

        sofa::helper::vector<Index> oldAttached = tris[t].attachedPoints;
        tris[t].attachedPoints.clear();

        for (unsigned int ip=0; ip<oldAttached.size(); ip++) {
            Index ptAttached = oldAttached[ip];
            Vec3 newBary;

            proj.ComputeBaryCoords(newBary, xProj[ptAttached],
                x0[ tri[0] ], x0[ tri[1] ], x0[ tri[2] ], false);

            int neg = 0;
            if (newBary[0] < -1e-20) neg++;
            if (newBary[1] < -1e-20) neg++;
            if (newBary[2] < -1e-20) neg++;

            if (neg == 0) {
                // Point still inside the triangle
                tris[t].attachedPoints.push_back(ptAttached);
                values[ptAttached][0] = newBary[0];
                values[ptAttached][1] = newBary[1];
                values[ptAttached][2] = newBary[2];
                continue;
            //} else if (neg == 1) {
            //    // Only one negative coordinates means the point lies in one of
            //    // the neighbouring triangles.
            //    // TODO
            //    continue;
            }

            // Retry projection but only inside the N1-ring
            Index newTri;
            proj.ProjectPoint(newBary, newTri, xProj[ptAttached], x0, N1);
            if (newTri == InvalidID) {
                serr << "Failed to project point " << ptAttached << "!" << sendl;
                continue;
            }

            // Reassign point to new triangle
            newAttached[newTri].push_back(ptAttached);

            // Update points
            Triangle newTri2 = m_container->getTriangle(newTri);
            indices[ptAttached][0] = newTri2[0];
            indices[ptAttached][1] = newTri2[1];
            indices[ptAttached][2] = newTri2[2];

            // Update barycentric coordinates
            values[ptAttached][0] = newBary[0];
            values[ptAttached][1] = newBary[1];
            values[ptAttached][2] = newBary[2];
        }
    }

    // Add newly attached points to the lists
    for (std::map<Index, helper::vector<Index> >::const_iterator i =
        newAttached.begin();
        i != newAttached.end(); i++) {
        for (unsigned int j=0; j < i->second.size(); j++) {
            tris[i->first].attachedPoints.push_back(i->second[j]);
        }
    }

    m_interpolationIndices.endEdit();
    m_interpolationValues.endEdit();
    triInfo.endEdit();
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if ((!vparams->displayFlags().getShowBehaviorModels()))
        return;

    if (!m_state) return;
    const VecCoord& x = m_state->read(sofa::core::ConstVecCoordId::position())->getValue();

    helper::vector<defaulttype::Vector3> boundary;
    helper::vector<defaulttype::Vector3> fixed;
    const helper::vector<PointInformation> &pts = pointInfo.getValue();
    if (pts.size() != x.size()) return;
    for (Index i=0; i < x.size(); i++) {
        if (pts[i].isFixed()) {
            fixed.push_back(x[i]);
        } else if (pts[i].isBoundary()) {
            boundary.push_back(x[i]);
        }
    }
    vparams->drawTool()->drawPoints(boundary, 10,
        sofa::defaulttype::Vec<4,float>(0.5, 0.5, 1.0, 1.0));
    vparams->drawTool()->drawPoints(fixed, 10,
        sofa::defaulttype::Vec<4,float>(0.8, 0.0, 0.8, 1.0));

    if (m_pointId != InvalidID) {
        // Draw tracked position
        helper::vector<defaulttype::Vector3> vv;
        vv.push_back(
            defaulttype::Vector3(m_point[0], m_point[1], m_point[2]));
        vparams->drawTool()->drawPoints(vv, 6,
            sofa::defaulttype::Vec<4,float>(1.0, 1.0, 1.0, 1.0));
        // Draw tracked position (projected in rest shape)
        //vv.clear();
        //vv.push_back(m_pointRest);
        //vparams->drawTool()->drawPoints(vv, 4,
        //    sofa::defaulttype::Vec<4,float>(1.0, 1.0, 0.0, 1.0));
        // Draw attached point
        vv.clear();
        vv.push_back(x[m_pointId]);
        vparams->drawTool()->drawPoints(vv, 6,
            /*cutting()
            ? sofa::defaulttype::Vec<4,float>(1.0, 1.0, 0.0, 1.0)
            :*/ sofa::defaulttype::Vec<4,float>(1.0, 1.0, 1.0, 1.0));
    }


    if (m_protectedEdges.size() > 0) {
        helper::vector<defaulttype::Vector3> points;
        for (VecIndex::const_iterator i=m_protectedEdges.begin();
            i != m_protectedEdges.end(); i++) {
            const Edge &e = m_container->getEdge(*i);
            points.push_back(x[ e[0] ]);
            points.push_back(x[ e[1] ]);
        }
        vparams->drawTool()->drawLines(points, 4,
            sofa::defaulttype::Vec<4,float>(0.0, 1.0, 0.0, 1.0));
    }

    m_surf.draw(vparams);
}

} // namespace controller

} // namespace component

} // namespace sofa

#undef OTHER
#undef CCW

#endif // SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_INL
