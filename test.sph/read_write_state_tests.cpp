//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <gtest/gtest.h>
#include <omp.h>
#include "vox.base/logging.h"
#include "vox.editor/simulator_base.h"
#include "vox.base/file_system.h"
#include "vox.sph/simulation.h"
#include "vox.base/time_manager.h"
#include "vox.sph/viscosity/viscosity_weiler2018.h"
#include "vox.sph/vorticity/micropolar_model_bender2017.h"

// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
#define new DEBUG_NEW
#endif
#endif

using namespace vox;
using namespace std;

SimulatorBase* base = nullptr;
Real eps = static_cast<Real>(1.0e-9);

struct FluidModelData {
    std::string m_id;

    std::vector<Real> m_masses;
    std::vector<Vector3r> m_a;
    std::vector<Vector3r> m_v0;
    std::vector<Vector3r> m_x0;
    std::vector<Vector3r> m_x;
    std::vector<Vector3r> m_v;
    std::vector<Real> m_density;
    std::vector<unsigned int> m_particleId;
    std::vector<unsigned int> m_objectId;
    std::vector<unsigned int> m_numNeighbors;
    std::vector<ParticleState> m_particleState;
    Real m_V;
    int m_surfaceTensionMethod;
    int m_viscosityMethod;
    int m_vorticityMethod;
    int m_dragMethod;
    int m_elasticityMethod;

    explicit FluidModelData(const FluidModel* fm) {
        Simulation* sim = Simulation::getCurrent();
        m_id = fm->getId();
        m_V = fm->getVolume(0);

        auto numParticles = fm->numActiveParticles();
        m_masses.resize(numParticles);
        m_a.resize(numParticles);
        m_v0.resize(numParticles);
        m_x0.resize(numParticles);
        m_x.resize(numParticles);
        m_v.resize(numParticles);
        m_density.resize(numParticles);
        m_particleId.resize(numParticles);
        m_objectId.resize(numParticles);
        m_particleState.resize(numParticles);
        m_numNeighbors.resize(numParticles);
        for (auto i = 0u; i < numParticles; i++) {
            m_masses[i] = fm->getMass(i);
            m_a[i] = fm->getAcceleration(i);
            m_v0[i] = fm->getVelocity0(i);
            m_x0[i] = fm->getPosition0(i);
            m_x[i] = fm->getPosition(i);
            m_v[i] = fm->getVelocity(i);
            m_density[i] = fm->getDensity(i);
            m_particleId[i] = fm->getParticleId(i);
            m_objectId[i] = fm->getObjectId(i);
            m_particleState[i] = fm->getParticleState(i);
            m_numNeighbors[i] = sim->numberOfNeighbors(0, 0, i);
        }

        m_surfaceTensionMethod = fm->getSurfaceTensionMethod();
        m_viscosityMethod = fm->getViscosityMethod();
        m_vorticityMethod = fm->getVorticityMethod();
        m_dragMethod = fm->getDragMethod();
        m_elasticityMethod = fm->getElasticityMethod();
    }
};

struct BoundaryModelData {
    bool m_sorted;
    unsigned int m_pointSetIndex;
    std::vector<Vector3r> m_x0;
    std::vector<Vector3r> m_x;
    std::vector<Vector3r> m_v;
    std::vector<Real> m_V;

    explicit BoundaryModelData(const BoundaryModel_Akinci2012* bm) {
        Simulation* sim = Simulation::getCurrent();
        m_sorted = bm->isSorted();
        m_pointSetIndex = bm->getPointSetIndex();

        auto numParticles = bm->numberOfParticles();
        m_x0.resize(numParticles);
        m_x.resize(numParticles);
        m_v.resize(numParticles);
        m_V.resize(numParticles);
        for (auto i = 0u; i < numParticles; i++) {
            m_x0[i] = bm->getPosition0(i);
            m_x[i] = bm->getPosition(i);
            m_v[i] = bm->getVelocity(i);
            m_V[i] = bm->getVolume(i);
        }
    }
};

struct NonpressureData {
    std::vector<Vector3r> m_omega;
    std::vector<Vector3r> m_vDiff;

    explicit NonpressureData(FluidModel* model) {
        auto numParticles = model->numActiveParticles();
        auto* micro = static_cast<MicropolarModel_Bender2017*>(model->getVorticityBase());
        auto* visco = static_cast<Viscosity_Weiler2018*>(model->getViscosityBase());
        m_omega.resize(numParticles);
        m_vDiff.resize(numParticles);
        for (auto i = 0u; i < numParticles; i++) {
            m_omega[i] = micro->getAngularVelocity(i);
            m_vDiff[i] = visco->getVDiff(i);
        }
    }
};

bool cmp(const Vector3r& x1, const Vector3r& x2) {
#ifdef USE_DOUBLE
    return (x1 - x2).norm() < eps;
#else
    return x1 == x2;
#endif
}

void compareStateFluidData(const FluidModelData* f1, const FluidModelData* f2) {
    auto numParticles = f1->m_x.size();
    auto numParticles2 = f2->m_x.size();
    ASSERT_TRUE(numParticles == numParticles2);

    bool chkPos = true;
    bool chkVel = true;
    bool chkId = true;
    bool chkObjId = true;
    for (auto i = 0; i < numParticles; i++) {
        if (!cmp(f1->m_x[i], f2->m_x[i])) chkPos = false;
        if (!cmp(f1->m_v[i], f2->m_v[i])) chkVel = false;
        if (f1->m_particleId[i] != f2->m_particleId[i]) chkId = false;
        if (f1->m_objectId[i] != f2->m_objectId[i]) chkObjId = false;
    }

    EXPECT_EQ(f1->m_id, f2->m_id);
    EXPECT_EQ(f1->m_surfaceTensionMethod, f2->m_surfaceTensionMethod);
    EXPECT_EQ(f1->m_viscosityMethod, f2->m_viscosityMethod);
    EXPECT_EQ(f1->m_vorticityMethod, f2->m_vorticityMethod);
    EXPECT_EQ(f1->m_dragMethod, f2->m_dragMethod);
    EXPECT_EQ(f1->m_elasticityMethod, f2->m_elasticityMethod);
    ASSERT_TRUE(chkPos);
    ASSERT_TRUE(chkVel);
    ASSERT_TRUE(chkId);
    ASSERT_TRUE(chkObjId);
}

void compareStateBoundaryData(const BoundaryModelData* b1, const BoundaryModelData* b2) {
    auto numParticles = b1->m_x.size();
    auto numParticles2 = b2->m_x.size();
    ASSERT_TRUE(numParticles == numParticles2);

    bool chkPos = true;
    bool chkVel = true;
    bool chkPos0 = true;
    bool chkV = true;
    for (auto i = 0; i < numParticles; i++) {
        if (!cmp(b1->m_x[i], b2->m_x[i])) chkPos = false;
        if (!cmp(b1->m_v[i], b2->m_v[i])) chkVel = false;
        if (!cmp(b1->m_x0[i], b2->m_x0[i])) chkPos0 = false;
        if (b1->m_V[i] != b2->m_V[i]) chkV = false;
    }

    ASSERT_TRUE(b1->m_sorted == b2->m_sorted);
    ASSERT_TRUE(b1->m_pointSetIndex == b2->m_pointSetIndex);
    ASSERT_TRUE(chkPos);
    ASSERT_TRUE(chkVel);
    ASSERT_TRUE(chkPos0);
    ASSERT_TRUE(chkV);
}

void compareStateMicropolarData(const NonpressureData* b1, const NonpressureData* b2) {
    auto numParticles = b1->m_omega.size();
    auto numParticles2 = b2->m_omega.size();
    ASSERT_TRUE(numParticles == numParticles2);

    bool chkOmega = true;
    bool chkVDiff = true;
    for (auto i = 0; i < numParticles; i++) {
        if (!cmp(b1->m_omega[i], b2->m_omega[i])) chkOmega = false;
        if (!cmp(b1->m_vDiff[i], b2->m_vDiff[i])) chkVDiff = false;
    }

    ASSERT_TRUE(chkOmega);
    ASSERT_TRUE(chkVDiff);
}

void compareAdditionalData(const FluidModelData* f1, const FluidModelData* f2) {
    auto numParticles = f1->m_x.size();
    auto numParticles2 = f2->m_x.size();
    ASSERT_TRUE(numParticles == numParticles2);

    bool chkPos0 = true;
    bool chkVel0 = true;
    bool chkDensity = true;
    bool chkNumNeighbors = true;
    for (auto i = 0; i < numParticles; i++) {
        if (!cmp(f1->m_x0[i], f2->m_x0[i])) chkPos0 = false;
        if (!cmp(f1->m_v0[i], f2->m_v0[i])) chkVel0 = false;
        if (f1->m_numNeighbors[i] != f2->m_numNeighbors[i]) chkNumNeighbors = false;
        if (abs(f1->m_density[i] - f2->m_density[i]) > 1.e-3) chkDensity = false;
    }

    ASSERT_TRUE(f1->m_id == f2->m_id);
    ASSERT_TRUE(f1->m_surfaceTensionMethod == f2->m_surfaceTensionMethod);
    ASSERT_TRUE(f1->m_viscosityMethod == f2->m_viscosityMethod);
    ASSERT_TRUE(f1->m_vorticityMethod == f2->m_vorticityMethod);
    ASSERT_TRUE(f1->m_dragMethod == f2->m_dragMethod);
    ASSERT_TRUE(f1->m_elasticityMethod == f2->m_elasticityMethod);
    ASSERT_TRUE(chkPos0);
    ASSERT_TRUE(chkVel0);
    ASSERT_TRUE(chkNumNeighbors);
    ASSERT_TRUE(chkDensity);
}

TEST(ReadWriteStateFileTests, read_write_state) {
    REPORT_MEMORY_LEAKS

    // no parallelism
    omp_set_num_threads(1);

    base = new SimulatorBase();

    // init parameters for simulator call
    std::string exePath = vox::utility::FileSystem::getProgramPath();
    std::string sceneFile =
            vox::utility::FileSystem::normalizePath(exePath + "/../data/Scenes/ReadWriteStateTest.json");

    int argc = 2;
    char** argv = new char*[2];
    argv[0] = new char[exePath.length() + 1];
    strcpy(argv[0], exePath.c_str());

    argv[1] = new char[sceneFile.length() + 1];
    strcpy(argv[1], sceneFile.c_str());

    base->init(argc, argv, "SPlisHSPlasH");

    base->setUseGUI(false);
    base->setUseParticleCaching(false);

    base->initSimulation();

    // set parameters
    Simulation* sim = Simulation::getCurrent();
    sim->setValue<int>(Simulation::CFL_METHOD, 0);
    sim->setValue<bool>(Simulation::ENABLE_Z_SORT, false);

    TimeManager::getCurrent()->setTimeStepSize(static_cast<Real>(0.01));

    // perform simulation for 10 steps and save state
    base->getBoundarySimulator()->initBoundaryData();
    for (unsigned int i = 0; i < 5; i++) base->timeStepNoGUI();

    base->setValue<Real>(SimulatorBase::PAUSE_AT, 2.342);
    base->setValue<Real>(SimulatorBase::STOP_AT, 5.2334);
    base->setValue<unsigned int>(SimulatorBase::NUM_STEPS_PER_RENDER, 9);
    base->setValue<Real>(SimulatorBase::DATA_EXPORT_FPS, 21.0);
    base->setValue<std::string>(SimulatorBase::PARTICLE_EXPORT_ATTRIBUTES, "velocity;density;factor");
    sim->setValue<Real>(Simulation::CFL_FACTOR, 0.77);
    sim->setValue<Real>(Simulation::CFL_MAX_TIMESTEPSIZE, 0.011);

    std::string stateFile = vox::utility::FileSystem::normalizePath(exePath + "/ReadWriteStateTest.bin");
    base->saveState(stateFile);

    // store data
    const Real W_zero = sim->W_zero();
    const bool doPause = base->getValue<bool>(SimulatorBase::PAUSE);
    const Real pauseAt = base->getValue<Real>(SimulatorBase::PAUSE_AT);
    const Real stopAt = base->getValue<Real>(SimulatorBase::STOP_AT);
    const auto numStepsPerRender = base->getValue<unsigned int>(SimulatorBase::NUM_STEPS_PER_RENDER);
    const int renderWalls = base->getValue<int>(SimulatorBase::RENDER_WALLS);
    const bool stateExport = base->getValue<bool>(SimulatorBase::STATE_EXPORT);
    const Real framesPerSecond = base->getValue<Real>(SimulatorBase::DATA_EXPORT_FPS);
    const Real framesPerSecondState = base->getValue<Real>(SimulatorBase::STATE_EXPORT_FPS);
    const auto particleAttributes = base->getValue<std::string>(SimulatorBase::PARTICLE_EXPORT_ATTRIBUTES);

    const bool sim2D = sim->getValue<bool>(Simulation::SIM_2D);
    const bool enableZSort = sim->getValue<bool>(Simulation::ENABLE_Z_SORT);
    const Real particleRadius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
    const Vector3r gravitation(sim->getVecValue<Real>(Simulation::GRAVITATION));
    const int cflMethod = sim->getValue<int>(Simulation::CFL_METHOD);
    const Real cflFactor = sim->getValue<Real>(Simulation::CFL_FACTOR);
    const Real cflMinTimeStepSize = sim->getValue<Real>(Simulation::CFL_MIN_TIMESTEPSIZE);
    const Real cflMaxTimeStepSize = sim->getValue<Real>(Simulation::CFL_MAX_TIMESTEPSIZE);
    const int kernel = sim->getValue<int>(Simulation::KERNEL_METHOD);
    const int gradKernel = sim->getValue<int>(Simulation::GRAD_KERNEL_METHOD);
    const int simulationMethod = sim->getValue<int>(Simulation::SIMULATION_METHOD);
    const int boundaryHandlingMethod = sim->getValue<int>(Simulation::BOUNDARY_HANDLING_METHOD);

    // copy data
    sim->performNeighborhoodSearch();
    sim->getTimeStep()->computeDensities(0);
    auto* fluidCopy1 = new FluidModelData(sim->getFluidModel(0));
    auto* boundaryCopy1 = new BoundaryModelData(static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(0)));
    auto* micropolarCopy1 = new NonpressureData(sim->getFluidModel(0));
    base->cleanup();

    delete base;
    base = nullptr;

    // second run
    base = new SimulatorBase();
    base->init(argc, argv, "SPlisHSPlasH");

    base->setUseGUI(false);

    base->initSimulation();
    sim = Simulation::getCurrent();
    sim->setValue<int>(Simulation::CFL_METHOD, 0);
    sim->setValue<bool>(Simulation::ENABLE_Z_SORT, false);

    base->getBoundarySimulator()->initBoundaryData();
    base->loadState(stateFile);
    sim->getNeighborhoodSearch()->set_radius(sim->getSupportRadius());

    sim->performNeighborhoodSearch();
    sim->getTimeStep()->computeDensities(0);

    ASSERT_TRUE(abs(TimeManager::getCurrent()->getTimeStepSize() - static_cast<Real>(0.01)) < eps);

    auto* fluidCopy2 = new FluidModelData(sim->getFluidModel(0));
    auto* boundaryCopy2 = new BoundaryModelData(static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(0)));
    auto* micropolarCopy2 = new NonpressureData(sim->getFluidModel(0));

    ASSERT_TRUE(W_zero == sim->W_zero());

    // compare write/read of state
    {
        compareStateFluidData(fluidCopy1, fluidCopy2);
        compareStateBoundaryData(boundaryCopy1, boundaryCopy2);
        compareStateMicropolarData(micropolarCopy1, micropolarCopy2);
    }

    // compare additional data
    { compareAdditionalData(fluidCopy1, fluidCopy2); }

    // parameters
    {
        ASSERT_TRUE(doPause == base->getValue<bool>(SimulatorBase::PAUSE));
        ASSERT_TRUE(pauseAt == base->getValue<Real>(SimulatorBase::PAUSE_AT));
        ASSERT_TRUE(stopAt == base->getValue<Real>(SimulatorBase::STOP_AT));
        ASSERT_TRUE(numStepsPerRender == base->getValue<unsigned int>(SimulatorBase::NUM_STEPS_PER_RENDER));
        ASSERT_TRUE(renderWalls == base->getValue<int>(SimulatorBase::RENDER_WALLS));
        ASSERT_TRUE(stateExport == base->getValue<bool>(SimulatorBase::STATE_EXPORT));
        ASSERT_TRUE(framesPerSecond == base->getValue<Real>(SimulatorBase::DATA_EXPORT_FPS));
        ASSERT_TRUE(framesPerSecondState == base->getValue<Real>(SimulatorBase::STATE_EXPORT_FPS));
        ASSERT_TRUE(particleAttributes == base->getValue<std::string>(SimulatorBase::PARTICLE_EXPORT_ATTRIBUTES));

        ASSERT_TRUE(sim2D == sim->getValue<bool>(Simulation::SIM_2D));
        ASSERT_TRUE(enableZSort == sim->getValue<bool>(Simulation::ENABLE_Z_SORT));
        ASSERT_TRUE(particleRadius == sim->getValue<Real>(Simulation::PARTICLE_RADIUS));
        const Vector3r gravitation2(sim->getVecValue<Real>(Simulation::GRAVITATION));
        ASSERT_TRUE(gravitation == gravitation2);
        ASSERT_TRUE(cflMethod == sim->getValue<int>(Simulation::CFL_METHOD));
        ASSERT_TRUE(cflFactor == sim->getValue<Real>(Simulation::CFL_FACTOR));
        ASSERT_TRUE(cflMinTimeStepSize == sim->getValue<Real>(Simulation::CFL_MIN_TIMESTEPSIZE));
        ASSERT_TRUE(cflMaxTimeStepSize == sim->getValue<Real>(Simulation::CFL_MAX_TIMESTEPSIZE));
        ASSERT_TRUE(kernel == sim->getValue<int>(Simulation::KERNEL_METHOD));
        ASSERT_TRUE(gradKernel == sim->getValue<int>(Simulation::GRAD_KERNEL_METHOD));
        ASSERT_TRUE(simulationMethod == sim->getValue<int>(Simulation::SIMULATION_METHOD));
        ASSERT_TRUE(boundaryHandlingMethod == sim->getValue<int>(Simulation::BOUNDARY_HANDLING_METHOD));
    }

    base->cleanup();
    delete base;
    delete fluidCopy1;
    delete fluidCopy2;
    delete[] argv[0];
    delete[] argv[1];
    delete[] argv;
}
