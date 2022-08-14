//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>

#include "py.pbd//bind_pointer_vector.h"
#include "py.pbd/common.h"
#include "vox.base/discrete_grid/triangle_mesh_distance.h"
#include "vox.pbd/cubic_sdf_collision_detection.h"
#include "vox.pbd/simulation.h"

namespace py = pybind11;

vox::CubicSDFCollisionDetection::GridPtr generateSDF(const std::vector<Vector3r>& vertices,
                                                     const std::vector<unsigned int>& faces,
                                                     const Eigen::Matrix<unsigned int, 3, 1>& resolution) {
    const unsigned int nFaces = faces.size() / 3;
#ifdef USE_DOUBLE
    vox::TriangleMesh sdfMesh((vertices[0]).data(), faces.data(), vertices.size(), nFaces);
#else
    // if type is float, copy vector to double vector
    std::vector<double> doubleVec;
    doubleVec.resize(3 * vertices.size());
    for (unsigned int i = 0; i < vertices.size(); i++)
        for (unsigned int j = 0; j < 3; j++) doubleVec[3 * i + j] = vertices[i][j];
    vox::TriangleMesh sdfMesh(&doubleVec[0], faces.data(), vertices.size(), nFaces);
#endif
    vox::TriangleMeshDistance md(sdfMesh);
    Eigen::AlignedBox3d domain;
    for (auto const& x : sdfMesh.vertices()) {
        domain.extend(x);
    }
    domain.max() += 0.1 * Eigen::Vector3d::Ones();
    domain.min() -= 0.1 * Eigen::Vector3d::Ones();

    std::cout << "Set SDF resolution: " << resolution[0] << ", " << resolution[1] << ", " << resolution[2] << std::endl;
    auto sdf = std::make_shared<vox::CubicSDFCollisionDetection::Grid>(
            domain, std::array<unsigned int, 3>({resolution[0], resolution[1], resolution[2]}));
    auto func = vox::DiscreteGrid::ContinuousFunction{};
    func = [&md](Eigen::Vector3d const& xi) { return md.signed_distance(xi).distance; };
    std::cout << "Generate SDF\n";
    sdf->addFunction(func, true);
    return sdf;
}

void SimulationModelModule(const py::module& m_sub) {
    py::class_<vox::TriangleModel>(m_sub, "TriangleModel")
            .def(py::init<>())
            .def("getParticleMesh", (const vox::TriangleModel::ParticleMesh& (vox::TriangleModel::*)()
                                             const)(&vox::TriangleModel::getParticleMesh))
            .def("cleanupModel", &vox::TriangleModel::cleanupModel)
            .def("getIndexOffset", &vox::TriangleModel::getIndexOffset)
            .def("initMesh", &vox::TriangleModel::initMesh)
            .def("updateMeshNormals", &vox::TriangleModel::updateMeshNormals)
            .def("getRestitutionCoeff", &vox::TriangleModel::getRestitutionCoeff)
            .def("setRestitutionCoeff", &vox::TriangleModel::setRestitutionCoeff)
            .def("getFrictionCoeff", &vox::TriangleModel::getFrictionCoeff)
            .def("setFrictionCoeff", &vox::TriangleModel::setFrictionCoeff);

    py::class_<vox::TetModel>(m_sub, "TetModel")
            .def(py::init<>())
            .def("getInitialX", &vox::TetModel::getInitialX)
            .def("setInitialX", &vox::TetModel::setInitialX)
            .def("getInitialR", &vox::TetModel::getInitialR)
            .def("setInitialR", &vox::TetModel::setInitialR)
            .def("getInitialScale", &vox::TetModel::getInitialScale)
            .def("setInitialScale", &vox::TetModel::setInitialScale)

            .def("getSurfaceMesh", &vox::TetModel::getSurfaceMesh)
            .def("getVisVertices", &vox::TetModel::getVisVertices)
            .def("getVisMesh", &vox::TetModel::getVisMesh)
            .def("getParticleMesh",
                 (const vox::TetModel::ParticleMesh& (vox::TetModel::*)() const)(&vox::TetModel::getParticleMesh))
            .def("cleanupModel", &vox::TetModel::cleanupModel)
            .def("getIndexOffset", &vox::TetModel::getIndexOffset)
            .def("initMesh", &vox::TetModel::initMesh)
            .def("updateMeshNormals", &vox::TetModel::updateMeshNormals)
            .def("attachVisMesh", &vox::TetModel::attachVisMesh)
            .def("updateVisMesh", &vox::TetModel::updateVisMesh)
            .def("getRestitutionCoeff", &vox::TetModel::getRestitutionCoeff)
            .def("setRestitutionCoeff", &vox::TetModel::setRestitutionCoeff)
            .def("getFrictionCoeff", &vox::TetModel::getFrictionCoeff)
            .def("setFrictionCoeff", &vox::TetModel::setFrictionCoeff);

    py::bind_pointer_vector<std::vector<vox::TriangleModel*>>(m_sub, "VecTriangleModels");
    py::bind_pointer_vector<std::vector<vox::TetModel*>>(m_sub, "VecTetModels");
    py::bind_pointer_vector<std::vector<vox::RigidBody*>>(m_sub, "VecRigidBodies");
    py::bind_vector<std::vector<vox::Constraint*>>(m_sub, "VecConstraints");

    py::class_<vox::SimulationModel, vox::ParameterObject>(m_sub, "SimulationModel")
            .def(py::init<>())
            .def("init", &vox::SimulationModel::init)
            .def("reset", &vox::SimulationModel::reset)
            .def("cleanup", &vox::SimulationModel::cleanup)
            .def("resetContacts", &vox::SimulationModel::resetContacts)
            .def("updateConstraints", &vox::SimulationModel::updateConstraints)
            .def("initConstraintGroups", &vox::SimulationModel::initConstraintGroups)
            .def(
                    "addTriangleModel",
                    [](vox::SimulationModel& model, std::vector<Vector3r>& points, std::vector<unsigned int>& indices,
                       const vox::TriangleModel::ParticleMesh::UVIndices& uvIndices,
                       const vox::TriangleModel::ParticleMesh::UVs& uvs, const bool testMesh) {
                        auto& triModels = model.getTriangleModels();
                        int i = triModels.size();
                        model.addTriangleModel(points.size(), indices.size() / 3, points.data(), indices.data(),
                                               uvIndices, uvs);
                        if (testMesh) {
                            vox::ParticleData& pd = model.getParticles();
                            const unsigned int nVert = triModels[i]->getParticleMesh().numVertices();
                            unsigned int offset = triModels[i]->getIndexOffset();
                            vox::Simulation* sim = vox::Simulation::getCurrent();
                            auto* cd = dynamic_cast<vox::CubicSDFCollisionDetection*>(
                                    sim->getTimeStep()->getCollisionDetection());
                            if (cd != nullptr)
                                cd->addCollisionObjectWithoutGeometry(
                                        i, vox::CollisionDetection::CollisionObject::TriangleModelCollisionObjectType,
                                        &pd.getPosition(offset), nVert, true);
                        }
                        return triModels[i];
                    },
                    py::arg("points"), py::arg("indices"),
                    py::arg("uvIndices") = vox::TriangleModel::ParticleMesh::UVIndices(),
                    py::arg("uvs") = vox::TriangleModel::ParticleMesh::UVs(), py::arg("testMesh") = false,
                    py::return_value_policy::reference)
            .def(
                    "addRegularTriangleModel",
                    [](vox::SimulationModel& model, const int width, const int height, const Vector3r& translation,
                       const Matrix3r& rotation, const Vector2r& scale, const bool testMesh) {
                        auto& triModels = model.getTriangleModels();
                        int i = triModels.size();
                        model.addRegularTriangleModel(width, height, translation, rotation, scale);
                        if (testMesh) {
                            vox::ParticleData& pd = model.getParticles();
                            const unsigned int nVert = triModels[i]->getParticleMesh().numVertices();
                            unsigned int offset = triModels[i]->getIndexOffset();
                            vox::Simulation* sim = vox::Simulation::getCurrent();
                            auto* cd = dynamic_cast<vox::CubicSDFCollisionDetection*>(
                                    sim->getTimeStep()->getCollisionDetection());
                            if (cd != nullptr)
                                cd->addCollisionObjectWithoutGeometry(
                                        i, vox::CollisionDetection::CollisionObject::TriangleModelCollisionObjectType,
                                        &pd.getPosition(offset), nVert, true);
                        }
                        return triModels[i];
                    },
                    py::arg("width"), py::arg("height"), py::arg("translation") = Vector3r::Zero(),
                    py::arg("rotation") = Matrix3r::Identity(), py::arg("scale") = Vector2r::Ones(),
                    py::arg("testMesh") = false, py::return_value_policy::reference)
            .def(
                    "addTetModel",
                    [](vox::SimulationModel& model, std::vector<Vector3r>& points, std::vector<unsigned int>& indices,
                       const bool testMesh, bool generateCollisionObject,
                       const Eigen::Matrix<unsigned int, 3, 1>& resolution) {
                        auto& tetModels = model.getTetModels();
                        int i = tetModels.size();
                        model.addTetModel(points.size(), indices.size() / 4, points.data(), indices.data());

                        vox::ParticleData& pd = model.getParticles();
                        vox::TetModel* tetModel = tetModels[i];
                        const unsigned int nVert = tetModel->getParticleMesh().numVertices();
                        unsigned int offset = tetModel->getIndexOffset();
                        vox::Simulation* sim = vox::Simulation::getCurrent();
                        if (generateCollisionObject) {
                            auto& surfaceMesh = tetModel->getSurfaceMesh();
                            vox::CubicSDFCollisionDetection::GridPtr sdf =
                                    generateSDF(points, surfaceMesh.getFaces(), resolution);
                            if (sdf != nullptr) {
                                auto* cd = dynamic_cast<vox::CubicSDFCollisionDetection*>(
                                        sim->getTimeStep()->getCollisionDetection());
                                if (cd != nullptr) {
                                    auto index = cd->getCollisionObjects().size();
                                    cd->addCubicSDFCollisionObject(
                                            i, vox::CollisionDetection::CollisionObject::TetModelCollisionObjectType,
                                            &pd.getPosition(offset), nVert, sdf, Vector3r::Ones(), testMesh, false);

                                    const unsigned int modelIndex = cd->getCollisionObjects()[index]->m_bodyIndex;
                                    vox::TetModel* tm = tetModels[modelIndex];
                                    const unsigned int offset = tm->getIndexOffset();
                                    const vox::utility::IndexedTetMesh& mesh = tm->getParticleMesh();

                                    ((vox::DistanceFieldCollisionDetection::DistanceFieldCollisionObject*)
                                             cd->getCollisionObjects()[index])
                                            ->initTetBVH(&pd.getPosition(offset), mesh.numVertices(),
                                                         mesh.getTets().data(), mesh.numTets(), cd->getTolerance());
                                }
                            }
                        } else if (testMesh) {
                            auto* cd = dynamic_cast<vox::CubicSDFCollisionDetection*>(
                                    sim->getTimeStep()->getCollisionDetection());
                            if (cd != nullptr) {
                                auto index = cd->getCollisionObjects().size();
                                cd->addCollisionObjectWithoutGeometry(
                                        i, vox::CollisionDetection::CollisionObject::TetModelCollisionObjectType,
                                        &pd.getPosition(offset), nVert, true);

                                const unsigned int modelIndex = cd->getCollisionObjects()[index]->m_bodyIndex;
                                vox::TetModel* tm = tetModels[modelIndex];
                                const unsigned int offset = tm->getIndexOffset();
                                const vox::utility::IndexedTetMesh& mesh = tm->getParticleMesh();

                                ((vox::DistanceFieldCollisionDetection::DistanceFieldCollisionObject*)
                                         cd->getCollisionObjects()[index])
                                        ->initTetBVH(&pd.getPosition(offset), mesh.numVertices(), mesh.getTets().data(),
                                                     mesh.numTets(), cd->getTolerance());
                            }
                        }
                        return tetModel;
                    },
                    py::arg("points"), py::arg("indices"), py::arg("testMesh") = false,
                    py::arg("generateCollisionObject") = false,
                    py::arg("resolution") = Eigen::Matrix<unsigned int, 3, 1>(30, 30, 30),
                    py::return_value_policy::reference)
            .def(
                    "addRegularTetModel",
                    [](vox::SimulationModel& model, const int width, const int height, const int depth,
                       const Vector3r& translation, const Matrix3r& rotation, const Vector3r& scale,
                       const bool testMesh) {
                        auto& tetModels = model.getTetModels();
                        int i = tetModels.size();
                        model.addRegularTetModel(width, height, depth, translation, rotation, scale);
                        if (testMesh) {
                            vox::ParticleData& pd = model.getParticles();
                            const unsigned int nVert = tetModels[i]->getParticleMesh().numVertices();
                            unsigned int offset = tetModels[i]->getIndexOffset();
                            vox::Simulation* sim = vox::Simulation::getCurrent();
                            auto* cd = dynamic_cast<vox::CubicSDFCollisionDetection*>(
                                    sim->getTimeStep()->getCollisionDetection());
                            if (cd != nullptr)
                                cd->addCollisionObjectWithoutGeometry(
                                        i, vox::CollisionDetection::CollisionObject::TetModelCollisionObjectType,
                                        &pd.getPosition(offset), nVert, true);
                        }
                        return tetModels[i];
                    },
                    py::arg("width"), py::arg("height"), py::arg("depth"), py::arg("translation") = Vector3r::Zero(),
                    py::arg("rotation") = Matrix3r::Identity(), py::arg("scale") = Vector3r::Ones(),
                    py::arg("testMesh") = false, py::return_value_policy::reference)
            .def("addLineModel",
                 [](vox::SimulationModel& model, const unsigned int nPoints, const unsigned int nQuaternions,
                    std::vector<Vector3r>& points, std::vector<Quaternionr>& quaternions,
                    std::vector<unsigned int>& indices, std::vector<unsigned int>& indicesQuaternions) {
                     model.addLineModel(nPoints, nQuaternions, points.data(), quaternions.data(), indices.data(),
                                        indicesQuaternions.data());
                 })
            .def("addBallJoint", &vox::SimulationModel::addBallJoint)
            .def("addBallOnLineJoint", &vox::SimulationModel::addBallOnLineJoint)
            .def("addHingeJoint", &vox::SimulationModel::addHingeJoint)
            .def("addTargetAngleMotorHingeJoint", &vox::SimulationModel::addTargetAngleMotorHingeJoint)
            .def("addTargetVelocityMotorHingeJoint", &vox::SimulationModel::addTargetVelocityMotorHingeJoint)
            .def("addUniversalJoint", &vox::SimulationModel::addUniversalJoint)
            .def("addSliderJoint", &vox::SimulationModel::addSliderJoint)
            .def("addTargetPositionMotorSliderJoint", &vox::SimulationModel::addTargetPositionMotorSliderJoint)
            .def("addTargetVelocityMotorSliderJoint", &vox::SimulationModel::addTargetVelocityMotorSliderJoint)
            .def("addRigidBodyParticleBallJoint", &vox::SimulationModel::addRigidBodyParticleBallJoint)
            .def("addRigidBodySpring", &vox::SimulationModel::addRigidBodySpring)
            .def("addDistanceJoint", &vox::SimulationModel::addDistanceJoint)
            .def("addDamperJoint", &vox::SimulationModel::addDamperJoint)
            .def("addRigidBodyContactConstraint", &vox::SimulationModel::addRigidBodyContactConstraint)
            .def("addParticleRigidBodyContactConstraint", &vox::SimulationModel::addParticleRigidBodyContactConstraint)
            .def("addParticleSolidContactConstraint", &vox::SimulationModel::addParticleSolidContactConstraint)
            .def("addDistanceConstraint", &vox::SimulationModel::addDistanceConstraint)
            .def("addDistanceConstraint_XPBD", &vox::SimulationModel::addDistanceConstraint_XPBD)
            .def("addDihedralConstraint", &vox::SimulationModel::addDihedralConstraint)
            .def("addIsometricBendingConstraint", &vox::SimulationModel::addIsometricBendingConstraint)
            .def("addIsometricBendingConstraint_XPBD", &vox::SimulationModel::addIsometricBendingConstraint_XPBD)
            .def("addFEMTriangleConstraint", &vox::SimulationModel::addFEMTriangleConstraint)
            .def("addStrainTriangleConstraint", &vox::SimulationModel::addStrainTriangleConstraint)
            .def("addVolumeConstraint", &vox::SimulationModel::addVolumeConstraint)
            .def("addVolumeConstraint_XPBD", &vox::SimulationModel::addVolumeConstraint_XPBD)
            .def("addFEMTetConstraint", &vox::SimulationModel::addFEMTetConstraint)
            .def("addStrainTetConstraint", &vox::SimulationModel::addStrainTetConstraint)
            .def("addShapeMatchingConstraint",
                 [](vox::SimulationModel& model, const unsigned int numberOfParticles,
                    const std::vector<unsigned int>& particleIndices, const std::vector<unsigned int>& numClusters,
                    const Real stiffness) {
                     model.addShapeMatchingConstraint(numberOfParticles, particleIndices.data(), numClusters.data(),
                                                      stiffness);
                 })

            .def("addStretchShearConstraint", &vox::SimulationModel::addStretchShearConstraint)
            .def("addBendTwistConstraint", &vox::SimulationModel::addBendTwistConstraint)
            .def("addStretchBendingTwistingConstraint", &vox::SimulationModel::addStretchBendingTwistingConstraint)
            .def("addDirectPositionBasedSolverForStiffRodsConstraint",
                 &vox::SimulationModel::addDirectPositionBasedSolverForStiffRodsConstraint)

            .def("getParticles", &vox::SimulationModel::getParticles, py::return_value_policy::reference)
            .def("getRigidBodies", &vox::SimulationModel::getRigidBodies, py::return_value_policy::reference)
            .def("getTriangleModels", &vox::SimulationModel::getTriangleModels, py::return_value_policy::reference)
            .def("getTetModels", &vox::SimulationModel::getTetModels, py::return_value_policy::reference)
            .def("getLineModels", &vox::SimulationModel::getLineModels, py::return_value_policy::reference)
            .def("getConstraints", &vox::SimulationModel::getConstraints, py::return_value_policy::reference)
            .def("getOrientations", &vox::SimulationModel::getOrientations, py::return_value_policy::reference)
            .def("getRigidBodyContactConstraints", &vox::SimulationModel::getRigidBodyContactConstraints,
                 py::return_value_policy::reference)
            .def("getParticleRigidBodyContactConstraints",
                 &vox::SimulationModel::getParticleRigidBodyContactConstraints, py::return_value_policy::reference)
            .def("getParticleSolidContactConstraints", &vox::SimulationModel::getParticleSolidContactConstraints,
                 py::return_value_policy::reference)
            .def("getConstraintGroups", &vox::SimulationModel::getConstraintGroups, py::return_value_policy::reference)
            .def("resetContacts", &vox::SimulationModel::resetContacts)

            .def("addClothConstraints", &vox::SimulationModel::addClothConstraints)
            .def("addBendingConstraints", &vox::SimulationModel::addBendingConstraints)
            .def("addSolidConstraints", &vox::SimulationModel::addSolidConstraints)

            .def("getContactStiffnessRigidBody", &vox::SimulationModel::getContactStiffnessRigidBody)
            .def("setContactStiffnessRigidBody", &vox::SimulationModel::setContactStiffnessRigidBody)
            .def("getContactStiffnessParticleRigidBody", &vox::SimulationModel::getContactStiffnessParticleRigidBody)
            .def("setContactStiffnessParticleRigidBody", &vox::SimulationModel::setContactStiffnessParticleRigidBody)

            .def(
                    "addRigidBody",
                    [](vox::SimulationModel& model, const Real density, const vox::VertexData& vertices,
                       const vox::utility::IndexedFaceMesh& mesh, const Vector3r& translation, const Matrix3r& rotation,
                       const Vector3r& scale, const bool testMesh, const vox::CubicSDFCollisionDetection::GridPtr& sdf) {
                        vox::SimulationModel::RigidBodyVector& rbs = model.getRigidBodies();
                        auto* rb = new vox::RigidBody();
                        rb->initBody(density, translation, Quaternionr(rotation), vertices, mesh, scale);
                        rbs.push_back(rb);
                        if (sdf != nullptr) {
                            vox::Simulation* sim = vox::Simulation::getCurrent();
                            auto* cd = dynamic_cast<vox::CubicSDFCollisionDetection*>(
                                    sim->getTimeStep()->getCollisionDetection());
                            if (cd != nullptr) {
                                const std::vector<Vector3r>& vertices =
                                        rb->getGeometry().getVertexDataLocal().getVertices();
                                const auto nVert = static_cast<unsigned int>(vertices.size());
                                cd->addCubicSDFCollisionObject(
                                        rbs.size() - 1,
                                        vox::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType,
                                        vertices.data(), nVert, sdf, scale, testMesh, false);
                            }
                        }
                        return rb;
                    },
                    py::arg("density"), py::arg("vertices"), py::arg("mesh"), py::arg("translation") = Vector3r::Zero(),
                    py::arg("rotation") = Matrix3r::Identity(), py::arg("scale") = Vector3r::Ones(),
                    py::arg("testMesh") = false, py::arg("sdf"), py::return_value_policy::reference)
            .def(
                    "addRigidBody",
                    [](vox::SimulationModel& model, const Real density, const vox::VertexData& vertices,
                       const vox::utility::IndexedFaceMesh& mesh, const Vector3r& translation, const Matrix3r& rotation,
                       const Vector3r& scale, const bool testMesh, const bool generateCollisionObject,
                       const Eigen::Matrix<unsigned int, 3, 1>& resolution) {
                        vox::Simulation* sim = vox::Simulation::getCurrent();
                        vox::SimulationModel::RigidBodyVector& rbs = model.getRigidBodies();
                        auto* rb = new vox::RigidBody();
                        auto i = rbs.size();
                        rb->initBody(density, translation, Quaternionr(rotation), vertices, mesh, scale);
                        rbs.push_back(rb);

                        if (generateCollisionObject) {
                            vox::CubicSDFCollisionDetection::GridPtr sdf =
                                    generateSDF(vertices.getVertices(), mesh.getFaces(), resolution);
                            if (sdf != nullptr) {
                                auto* cd = dynamic_cast<vox::CubicSDFCollisionDetection*>(
                                        sim->getTimeStep()->getCollisionDetection());
                                if (cd != nullptr) {
                                    const std::vector<Vector3r>& vertices =
                                            rb->getGeometry().getVertexDataLocal().getVertices();
                                    const auto nVert = static_cast<unsigned int>(vertices.size());
                                    cd->addCubicSDFCollisionObject(
                                            rbs.size() - 1,
                                            vox::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType,
                                            vertices.data(), nVert, sdf, scale, testMesh, false);
                                }
                            }
                        } else if (testMesh) {
                            auto* cd = dynamic_cast<vox::CubicSDFCollisionDetection*>(
                                    sim->getTimeStep()->getCollisionDetection());
                            if (cd != nullptr) {
                                auto index = cd->getCollisionObjects().size();
                                const std::vector<Vector3r>& vertices =
                                        rbs[i]->getGeometry().getVertexDataLocal().getVertices();
                                const auto nVert = static_cast<unsigned int>(vertices.size());
                                cd->addCollisionObjectWithoutGeometry(
                                        i, vox::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType,
                                        vertices.data(), nVert, true);
                            }
                        }
                        return rb;
                    },
                    py::arg("density"), py::arg("vertices"), py::arg("mesh"), py::arg("translation") = Vector3r::Zero(),
                    py::arg("rotation") = Matrix3r::Identity(), py::arg("scale") = Vector3r::Ones(),
                    py::arg("testMesh") = false, py::arg("generateCollisionObject") = false,
                    py::arg("resolution") = Eigen::Matrix<unsigned int, 3, 1>(30, 30, 30),
                    py::return_value_policy::reference);
}