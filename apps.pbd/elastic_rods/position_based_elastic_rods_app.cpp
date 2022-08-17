//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "apps.pbd/common/demo_base.h"
#include "apps.pbd/common/mini_gl.h"
#include "apps.pbd/elastic_rods/position_based_elastic_rods_model.h"
#include "apps.pbd/elastic_rods/position_based_elastic_rods_tsc.h"
#include "vox.base/logging.h"
#include "vox.base/time_manager.h"
#include "vox.base/timing.h"
#include "vox.pbd/simulation.h"

// Enable memory leak detection
#if defined(_DEBUG) && !defined(EIGEN_ALIGN)
#define new DEBUG_NEW
#endif

using namespace vox;
using namespace Eigen;
using namespace std;

void timeStep();
void buildModel();
void createRod();
void render();
void reset();

void TW_CALL setDistanceStiffness(const void *value, void *clientData);
void TW_CALL getDistanceStiffness(void *value, void *clientData);

void TW_CALL setRestDarbouxX(const void *value, void *clientData);
void TW_CALL setRestDarbouxY(const void *value, void *clientData);
void TW_CALL setRestDarbouxZ(const void *value, void *clientData);

void TW_CALL getRestDarbouxX(void *value, void *clientData);
void TW_CALL getRestDarbouxY(void *value, void *clientData);
void TW_CALL getRestDarbouxZ(void *value, void *clientData);

void TW_CALL getBendingAndTwistingStiffnessX(void *value, void *clientData);
void TW_CALL getBendingAndTwistingStiffnessY(void *value, void *clientData);
void TW_CALL getBendingAndTwistingStiffnessZ(void *value, void *clientData);

void TW_CALL setBendingAndTwistingStiffnessX(const void *value, void *clientData);
void TW_CALL setBendingAndTwistingStiffnessY(const void *value, void *clientData);
void TW_CALL setBendingAndTwistingStiffnessZ(const void *value, void *clientData);

DemoBase *base;
PositionBasedElasticRodsTSC sim;
Real distanceStiffness = 1.0;

const int numberOfPoints = 32;

// main
int main(int argc, char **argv) {
    REPORT_MEMORY_LEAKS

    base = new DemoBase();
    base->init(argc, argv, "Elastic rod demo");

    auto *model = new PositionBasedElasticRodsModel();
    model->init();
    Simulation::getCurrent()->setModel(model);
    auto *tsc = new PositionBasedElasticRodsTSC();
    tsc->init();
    delete Simulation::getCurrent()->getTimeStep();
    Simulation::getCurrent()->setTimeStep(tsc);

    buildModel();

    base->createParameterGUI();

    // OpenGL
    MiniGL::setClientIdleFunc(timeStep);
    MiniGL::addKeyFunc('r', reset);
    MiniGL::setClientSceneFunc(render);
    MiniGL::setViewport(40.0f, 0.1f, 500.0f, Vector3r(5.0, 5.0, 10.0), Vector3r(5.0, 0.0, 0.0));

    TwAddVarCB(MiniGL::getTweakBar(), "DistanceStiffness", TW_TYPE_REAL, setDistanceStiffness, getDistanceStiffness,
               model, " label='Distance constraint stiffness'  min=0.0 step=0.1 precision=4 group='BendTwist' ");

    TwAddVarCB(MiniGL::getTweakBar(), "RestDarbouxX", TW_TYPE_REAL, setRestDarbouxX, getRestDarbouxX, model,
               " label='Rest Darboux X'  min=-1.0 max = 1.0 step=0.01 precision=2 group=BendTwist ");
    TwAddVarCB(MiniGL::getTweakBar(), "RestDarbouxY", TW_TYPE_REAL, setRestDarbouxY, getRestDarbouxY, model,
               " label='Rest Darboux Y'  min=-1.0 max = 1.0 step=0.01 precision=2 group=BendTwist ");
    TwAddVarCB(MiniGL::getTweakBar(), "RestDarbouxZ", TW_TYPE_REAL, setRestDarbouxZ, getRestDarbouxZ, model,
               " label='Rest Darboux Z'  min=-1.0 max = 1.0 step=0.01 precision=2 group=BendTwist ");

    TwAddVarCB(MiniGL::getTweakBar(), "BendingAndTwistingStiffnessX", TW_TYPE_REAL, setBendingAndTwistingStiffnessX,
               getBendingAndTwistingStiffnessX, model,
               " label='Bending X stiffness'  min=0.01 max = 1.0 step=0.01 precision=2 group=BendTwist ");
    TwAddVarCB(MiniGL::getTweakBar(), "BendingAndTwistingStiffnessY", TW_TYPE_REAL, setBendingAndTwistingStiffnessY,
               getBendingAndTwistingStiffnessY, model,
               " label='Bending Y stiffness'  min=0.01 max = 1.0 step=0.01 precision=2 group=BendTwist ");
    TwAddVarCB(MiniGL::getTweakBar(), "BendingAndTwistingStiffnessZ", TW_TYPE_REAL, setBendingAndTwistingStiffnessZ,
               getBendingAndTwistingStiffnessZ, model,
               " label='Twisting stiffness'  min=0.01 max = 1.0 step=0.01 precision=2 group=BendTwist ");

    MiniGL::mainLoop();

    Timing::printAverageTimes();
    Timing::printTimeSums();

    delete Simulation::getCurrent();
    delete base;
    delete model;

    return 0;
}

void reset() {
    Timing::printAverageTimes();
    Timing::reset();

    Simulation::getCurrent()->reset();
    base->getSelectedParticles().clear();
}

void buildModel() {
    TimeManager::getCurrent()->setTimeStepSize(0.002f);
    auto *model = (PositionBasedElasticRodsModel *)Simulation::getCurrent()->getModel();
    model->setBendingAndTwistingStiffness(Vector3r(0.5, 0.5, 0.5));

    sim.setDamping(0.001f);

    createRod();
}

void timeStep() {
    const Real pauseAt = base->getValue<Real>(DemoBase::PAUSE_AT);
    if ((pauseAt > 0.0) && (pauseAt < TimeManager::getCurrent()->getTime())) base->setValue(DemoBase::PAUSE, true);

    if (base->getValue<bool>(DemoBase::PAUSE)) return;

    // Simulation code
    SimulationModel *model = Simulation::getCurrent()->getModel();
    const auto numSteps = base->getValue<unsigned int>(DemoBase::NUM_STEPS_PER_RENDER);
    for (unsigned int i = 0; i < numSteps; i++) {
        START_TIMING("SimStep")
        Simulation::getCurrent()->getTimeStep()->step(*model);
        STOP_TIMING_AVG
    }
}

void render() {
    base->render();

    // Draw sim model
    auto *model = (PositionBasedElasticRodsModel *)Simulation::getCurrent()->getModel();
    ParticleData &pd = model->getParticles();
    ParticleData &ghostParticles = model->getGhostParticles();
    SimulationModel::ConstraintVector &constraints = model->getConstraints();

    float pointColor[4] = {0.1f, 0.2f, 0.6f, 1};
    float ghostPointColor[4] = {0.1f, 0.1f, 0.1f, 0.5f};
    float edgeColor[4] = {0.0f, 0.6f, 0.2f, 1};

    for (unsigned int i = 0; i < numberOfPoints; i++) {
        MiniGL::drawSphere(pd.getPosition(i), 0.07f, pointColor);
    }

    for (unsigned int i = 0; i < numberOfPoints - 1; i++) {
        MiniGL::drawSphere(ghostParticles.getPosition(i), 0.07f, ghostPointColor);
        MiniGL::drawVector(pd.getPosition(i), pd.getPosition(i + 1), 0.2f, edgeColor);
    }

    MiniGL::drawTime(TimeManager::getCurrent()->getTime());
}

/** Create the elastic rod model
 */
void createRod() {
    auto *model = (PositionBasedElasticRodsModel *)Simulation::getCurrent()->getModel();
    ParticleData &particles = model->getParticles();
    ParticleData &ghostParticles = model->getGhostParticles();
    SimulationModel::ConstraintVector &constraints = model->getConstraints();

    // centreline points
    for (unsigned int i = 0; i < numberOfPoints; i++) {
        particles.addVertex(Vector3r(static_cast<Real>(0.25 * i), 0.0, 0.0));
    }

    // edge ghost points
    for (unsigned int i = 0; i < numberOfPoints - 1; i++) {
        ghostParticles.addVertex(
                Vector3r(static_cast<Real>(0.25 * i) + static_cast<Real>(0.125), static_cast<Real>(0.25), 0.0));
    }

    // lock two first particles and first ghost point
    particles.setMass(0, 0.0f);
    particles.setMass(1, 0.0f);
    ghostParticles.setMass(0, 0.0f);

    for (unsigned int i = 0; i < numberOfPoints - 1; i++) {
        model->addDistanceConstraint(i, i + 1, distanceStiffness);
        model->addPerpendiculaBisectorConstraint(i, i + 1, i);
        model->addGhostPointEdgeDistanceConstraint(i, i + 1, i);

        if (i < numberOfPoints - 2) {
            //  Single rod element:
            //      D   E		//ghost points
            //		|	|
            //  --A---B---C--	// rod points
            int pA = i;
            int pB = i + 1;
            int pC = i + 2;
            int pD = i;
            int pE = i + 1;
            model->addDarbouxVectorConstraint(pA, pB, pC, pD, pE);
        }
    }
}

void TW_CALL setRestDarbouxX(const void *value, void *clientData) {
    ((PositionBasedElasticRodsModel *)clientData)->getRestDarbouxVector()[0] = *(const Real *)(value);
}

void TW_CALL setRestDarbouxY(const void *value, void *clientData) {
    ((PositionBasedElasticRodsModel *)clientData)->getRestDarbouxVector()[1] = *(const Real *)(value);
}

void TW_CALL setRestDarbouxZ(const void *value, void *clientData) {
    ((PositionBasedElasticRodsModel *)clientData)->getRestDarbouxVector()[2] = *(const Real *)(value);
}

void TW_CALL getRestDarbouxX(void *value, void *clientData) {
    *(Real *)(value) = ((PositionBasedElasticRodsModel *)clientData)->getRestDarbouxVector()[0];
}

void TW_CALL getRestDarbouxY(void *value, void *clientData) {
    *(Real *)(value) = ((PositionBasedElasticRodsModel *)clientData)->getRestDarbouxVector()[1];
}

void TW_CALL getRestDarbouxZ(void *value, void *clientData) {
    *(Real *)(value) = ((PositionBasedElasticRodsModel *)clientData)->getRestDarbouxVector()[2];
}

void TW_CALL setBendingAndTwistingStiffnessX(const void *value, void *clientData) {
    ((PositionBasedElasticRodsModel *)clientData)->getBendingAndTwistingStiffness()[0] = *(const Real *)(value);
}

void TW_CALL setBendingAndTwistingStiffnessY(const void *value, void *clientData) {
    ((PositionBasedElasticRodsModel *)clientData)->getBendingAndTwistingStiffness()[1] = *(const Real *)(value);
}

void TW_CALL setBendingAndTwistingStiffnessZ(const void *value, void *clientData) {
    ((PositionBasedElasticRodsModel *)clientData)->getBendingAndTwistingStiffness()[2] = *(const Real *)(value);
}

void TW_CALL getBendingAndTwistingStiffnessX(void *value, void *clientData) {
    *(Real *)(value) = ((PositionBasedElasticRodsModel *)clientData)->getBendingAndTwistingStiffness()[0];
}

void TW_CALL getBendingAndTwistingStiffnessY(void *value, void *clientData) {
    *(Real *)(value) = ((PositionBasedElasticRodsModel *)clientData)->getBendingAndTwistingStiffness()[1];
}

void TW_CALL getBendingAndTwistingStiffnessZ(void *value, void *clientData) {
    *(Real *)(value) = ((PositionBasedElasticRodsModel *)clientData)->getBendingAndTwistingStiffness()[2];
}

void TW_CALL setDistanceStiffness(const void *value, void *clientData) {
    distanceStiffness = *(const Real *)(value);
    ((SimulationModel *)clientData)
            ->setConstraintValue<DistanceConstraint, Real, &DistanceConstraint::m_stiffness>(distanceStiffness);
    ((SimulationModel *)clientData)
            ->setConstraintValue<DistanceConstraint_XPBD, Real, &DistanceConstraint_XPBD::m_stiffness>(
                    distanceStiffness);
}

void TW_CALL getDistanceStiffness(void *value, void *clientData) { *(Real *)(value) = distanceStiffness; }
