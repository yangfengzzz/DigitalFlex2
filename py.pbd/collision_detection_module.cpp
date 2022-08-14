//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include <utility>

#include "py.pbd/common.h"
#include "vox.base/geometry/triangle_mesh_distance.h"
#include "vox.base/mesh/triangle_mesh.h"
#include "vox.pbd/collision_detection.h"
#include "vox.pbd/cubic_sdf_collision_detection.h"

namespace py = pybind11;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

void CollisionDetectionModule(const py::module& m_sub) {
    py::class_<vox::CollisionDetection::CollisionObject>(m_sub, "CollisionObject")
            .def_readonly_static("RigidBodyCollisionObjectType",
                                 &vox::CollisionDetection::CollisionObject::RigidBodyCollisionObjectType)
            .def_readonly_static("TriangleModelCollisionObjectType",
                                 &vox::CollisionDetection::CollisionObject::TriangleModelCollisionObjectType)
            .def_readonly_static("TetModelCollisionObjectType",
                                 &vox::CollisionDetection::CollisionObject::TetModelCollisionObjectType)
            .def_readwrite("aabb", &vox::CollisionDetection::CollisionObject::m_aabb)
            .def_readwrite("bodyIndex", &vox::CollisionDetection::CollisionObject::m_bodyIndex)
            .def_readwrite("bodyType", &vox::CollisionDetection::CollisionObject::m_bodyType);

    py::class_<vox::CollisionDetection::CollisionObjectWithoutGeometry, vox::CollisionDetection::CollisionObject>(
            m_sub, "CollisionObjectWithoutGeometry")
            .def(py::init<>())
            .def_readwrite_static("TYPE_ID", &vox::CollisionDetection::CollisionObjectWithoutGeometry::TYPE_ID)
            .def("getTypeId", &vox::CollisionDetection::CollisionObjectWithoutGeometry::getTypeId);

    py::class_<vox::CollisionDetection>(m_sub, "CollisionDetection")
            .def_readonly_static("RigidBodyContactType", &vox::CollisionDetection::RigidBodyContactType)
            .def_readonly_static("ParticleContactType", &vox::CollisionDetection::ParticleContactType)
            .def_readonly_static("ParticleRigidBodyContactType", &vox::CollisionDetection::ParticleRigidBodyContactType)
            .def_readonly_static("ParticleSolidContactType", &vox::CollisionDetection::ParticleSolidContactType)

            .def("cleanup", &vox::CollisionDetection::cleanup)
            .def("getTolerance", &vox::CollisionDetection::getTolerance)
            .def("setTolerance", &vox::CollisionDetection::setTolerance)
            .def("addRigidBodyContact", &vox::CollisionDetection::addRigidBodyContact)
            .def("addParticleRigidBodyContact", &vox::CollisionDetection::addParticleRigidBodyContact)
            .def("addParticleSolidContact", &vox::CollisionDetection::addParticleSolidContact)
            .def("addCollisionObject", &vox::CollisionDetection::addCollisionObject)
            .def("getCollisionObjects", &vox::CollisionDetection::getCollisionObjects)
            .def("collisionDetection", &vox::CollisionDetection::collisionDetection)
            //.def("setContactCallback", &vox::CollisionDetection::setContactCallback)
            //.def("setSolidContactCallback", &vox::CollisionDetection::setSolidContactCallback)
            .def("updateAABBs", &vox::CollisionDetection::updateAABBs)
            .def("updateAABB", overload_cast_<vox::SimulationModel&, vox::CollisionDetection::CollisionObject*>()(
                                       &vox::CollisionDetection::updateAABB));

    py::class_<vox::DistanceFieldCollisionDetection, vox::CollisionDetection>(m_sub, "DistanceFieldCollisionDetection")
            .def(py::init<>())
            .def("addCollisionBox",
                 [](vox::DistanceFieldCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::ParticleData& pd, const unsigned int offset, const unsigned int numVertices,
                    const Vector3r& box, const bool testMesh, const bool invertSDF) {
                     cd.addCollisionBox(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, box, testMesh,
                                        invertSDF);
                 })
            .def("addCollisionBox",
                 [](vox::DistanceFieldCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::VertexData& pd, const unsigned int offset, const unsigned int numVertices,
                    const Vector3r& box, const bool testMesh, const bool invertSDF) {
                     cd.addCollisionBox(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, box, testMesh,
                                        invertSDF);
                 })

            .def("addCollisionSphere",
                 [](vox::DistanceFieldCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::ParticleData& pd, const unsigned int offset, const unsigned int numVertices,
                    const Real radius, const bool testMesh, const bool invertSDF) {
                     cd.addCollisionSphere(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, radius, testMesh,
                                           invertSDF);
                 })
            .def("addCollisionSphere",
                 [](vox::DistanceFieldCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::VertexData& pd, const unsigned int offset, const unsigned int numVertices,
                    const Real radius, const bool testMesh, const bool invertSDF) {
                     cd.addCollisionSphere(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, radius, testMesh,
                                           invertSDF);
                 })

            .def("addCollisionTorus",
                 [](vox::DistanceFieldCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::ParticleData& pd, const unsigned int offset, const unsigned int numVertices,
                    const Vector2r& radii, const bool testMesh, const bool invertSDF) {
                     cd.addCollisionTorus(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, radii, testMesh,
                                          invertSDF);
                 })
            .def("addCollisionTorus",
                 [](vox::DistanceFieldCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::VertexData& pd, const unsigned int offset, const unsigned int numVertices,
                    const Vector2r& radii, const bool testMesh, const bool invertSDF) {
                     cd.addCollisionTorus(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, radii, testMesh,
                                          invertSDF);
                 })

            .def("addCollisionCylinder",
                 [](vox::DistanceFieldCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::ParticleData& pd, const unsigned int offset, const unsigned int numVertices,
                    const Vector2r& dim, const bool testMesh, const bool invertSDF) {
                     cd.addCollisionCylinder(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, dim, testMesh,
                                             invertSDF);
                 })
            .def("addCollisionCylinder",
                 [](vox::DistanceFieldCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::VertexData& pd, const unsigned int offset, const unsigned int numVertices,
                    const Vector2r& dim, const bool testMesh, const bool invertSDF) {
                     cd.addCollisionCylinder(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, dim, testMesh,
                                             invertSDF);
                 })

            .def("addCollisionHollowSphere",
                 [](vox::DistanceFieldCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::ParticleData& pd, const unsigned int offset, const unsigned int numVertices,
                    const Real radius, const Real thickness, const bool testMesh, const bool invertSDF) {
                     cd.addCollisionHollowSphere(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, radius,
                                                 thickness, testMesh, invertSDF);
                 })
            .def("addCollisionHollowSphere",
                 [](vox::DistanceFieldCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::VertexData& pd, const unsigned int offset, const unsigned int numVertices,
                    const Real radius, const Real thickness, const bool testMesh, const bool invertSDF) {
                     cd.addCollisionHollowSphere(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, radius,
                                                 thickness, testMesh, invertSDF);
                 })

            .def("addCollisionHollowBox",
                 [](vox::DistanceFieldCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::ParticleData& pd, const unsigned int offset, const unsigned int numVertices,
                    const Vector3r& box, const Real thickness, const bool testMesh, const bool invertSDF) {
                     cd.addCollisionHollowBox(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, box, thickness,
                                              testMesh, invertSDF);
                 })
            .def("addCollisionHollowBox",
                 [](vox::DistanceFieldCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::VertexData& pd, const unsigned int offset, const unsigned int numVertices,
                    const Vector3r& box, const Real thickness, const bool testMesh, const bool invertSDF) {
                     cd.addCollisionHollowBox(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, box, thickness,
                                              testMesh, invertSDF);
                 })

            .def("addCollisionObjectWithoutGeometry",
                 [](vox::DistanceFieldCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::ParticleData& pd, const unsigned int offset, const unsigned int numVertices,
                    const bool testMesh) {
                     cd.addCollisionObjectWithoutGeometry(bodyIndex, bodyType, &pd.getPosition(offset), numVertices,
                                                          testMesh);
                 })
            .def("addCollisionObjectWithoutGeometry",
                 [](vox::DistanceFieldCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::VertexData& vd, const unsigned int offset, const unsigned int numVertices,
                    const bool testMesh) {
                     cd.addCollisionObjectWithoutGeometry(bodyIndex, bodyType, &vd.getPosition(offset), numVertices,
                                                          testMesh);
                 })

            .def("isDistanceFieldCollisionObject",
                 &vox::DistanceFieldCollisionDetection::isDistanceFieldCollisionObject);

    py::class_<vox::CubicSDFCollisionDetection::Grid, std::shared_ptr<vox::CubicSDFCollisionDetection::Grid>>(
            m_sub, "CubicSDFCollisionDetectionGridPtr");

    py::class_<vox::CubicSDFCollisionDetection, vox::DistanceFieldCollisionDetection>(m_sub,
                                                                                      "CubicSDFCollisionDetection")
            .def(py::init<>())
            .def("addCubicSDFCollisionObject",
                 [](vox::CubicSDFCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::ParticleData& pd, const unsigned int offset, const unsigned int numVertices,
                    const std::string& sdfFile, const Vector3r& scale, const bool testMesh, const bool invertSDF) {
                     cd.addCubicSDFCollisionObject(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, sdfFile,
                                                   scale, testMesh, invertSDF);
                 })
            .def("addCubicSDFCollisionObject",
                 [](vox::CubicSDFCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::VertexData& pd, const unsigned int offset, const unsigned int numVertices,
                    const std::string& sdfFile, const Vector3r& scale, const bool testMesh, const bool invertSDF) {
                     cd.addCubicSDFCollisionObject(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, sdfFile,
                                                   scale, testMesh, invertSDF);
                 })
            .def("addCubicSDFCollisionObject",
                 [](vox::CubicSDFCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::ParticleData& pd, const unsigned int offset, const unsigned int numVertices,
                    vox::CubicSDFCollisionDetection::GridPtr sdf, const Vector3r& scale, const bool testMesh,
                    const bool invertSDF) {
                     cd.addCubicSDFCollisionObject(bodyIndex, bodyType, &pd.getPosition(offset), numVertices,
                                                   std::move(sdf), scale, testMesh, invertSDF);
                 })
            .def("addCubicSDFCollisionObject",
                 [](vox::CubicSDFCollisionDetection& cd, const unsigned int bodyIndex, const unsigned int bodyType,
                    const vox::VertexData& pd, const unsigned int offset, const unsigned int numVertices,
                    const vox::CubicSDFCollisionDetection::GridPtr& sdf, const Vector3r& scale, const bool testMesh,
                    const bool invertSDF) {
                     cd.addCubicSDFCollisionObject(bodyIndex, bodyType, &pd.getPosition(offset), numVertices, sdf,
                                                   scale, testMesh, invertSDF);
                 })
            .def_static("generateSDF",
                        [](const vox::VertexData& vd, const vox::utility::IndexedFaceMesh& mesh,
                           const Eigen::Matrix<unsigned int, 3, 1>& resolution)
                                -> vox::CubicSDFCollisionDetection::GridPtr {
                            const std::vector<unsigned int>& faces = mesh.getFaces();
                            const unsigned int nFaces = mesh.numFaces();

#ifdef USE_DOUBLE
                            vox::TriangleMesh sdfMesh(&vd.getPosition(0)[0], faces.data(), vd.size(), nFaces);
#else
                // if type is float, copy vector to double vector
                std::vector<double> doubleVec;
                doubleVec.resize(3 * vd.size());
                for (unsigned int i = 0; i < vd.size(); i++)
                    for (unsigned int j = 0; j < 3; j++)
                        doubleVec[3 * i + j] = vd.getPosition(i)[j];
                vox::TriangleMesh sdfMesh(&doubleVec[0], faces.data(), vd.size(), nFaces);
#endif
                            vox::TriangleMeshDistance md(sdfMesh);
                            Eigen::AlignedBox3d domain;
                            for (auto const& x : sdfMesh.vertices()) {
                                domain.extend(x);
                            }
                            domain.max() += 0.1 * Eigen::Vector3d::Ones();
                            domain.min() -= 0.1 * Eigen::Vector3d::Ones();

                            std::cout << "Set SDF resolution: " << resolution[0] << ", " << resolution[1] << ", "
                                      << resolution[2] << std::endl;
                            // vox::CubicSDFCollisionDetection::Grid *sdf = new
                            // vox::CubicSDFCollisionDetection::Grid(domain, std::array<unsigned int, 3>({
                            // resolution[0], resolution[1], resolution[2] }));
                            auto sdf = std::make_shared<vox::CubicSDFCollisionDetection::Grid>(
                                    domain, std::array<unsigned int, 3>({resolution[0], resolution[1], resolution[2]}));
                            auto func = vox::DiscreteGrid::ContinuousFunction{};
                            func = [&md](Eigen::Vector3d const& xi) { return md.signed_distance(xi).distance; };
                            std::cout << "Generate SDF\n";
                            sdf->addFunction(func, true);
                            return sdf;
                        });
}