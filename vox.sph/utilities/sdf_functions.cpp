//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/utilities/sdf_functions.h"

#include "vox.base/discrete_grid/triangle_mesh_distance.h"
#include "vox.base/mesh/triangle_mesh.h"
#include "vox.base/timing.h"
#include "vox.base/obj_loader.h"

using namespace Eigen;
using namespace std;
namespace vox::utility {

AlignedBox3r SDFFunctions::computeBoundingBox(const unsigned int numVertices, const Vector3r *vertices) {
    AlignedBox3r box;

    // compute bounding box
    box.min() = vertices[0];
    box.max() = box.min();
    box.setEmpty();
    for (unsigned int i = 1; i < numVertices; ++i) {
        const Vector3r &p = vertices[i];
        box.extend(p);
    }
    return box;
}

double SDFFunctions::distance(CubicLagrangeDiscreteGrid *sdf,
                              const Vector3r &x,
                              const Real thickness,
                              Vector3r &normal,
                              Vector3r &nextSurfacePoint) {
    Eigen::Vector3d n;
    const double dist = sdf->interpolate(0, x.template cast<double>(), &n);
    if (dist == std::numeric_limits<double>::max()) return dist;
    normal.normalize();
    normal = n.template cast<Real>();

    nextSurfacePoint = (x - dist * normal);

    return dist - thickness;
}

double SDFFunctions::distance(CubicLagrangeDiscreteGrid *sdf, const Vector3r &x, const Real thickness) {
    const double dist = sdf->interpolate(0, x.template cast<double>());
    if (dist == std::numeric_limits<double>::max()) return dist;
    return dist - thickness;
}

CubicLagrangeDiscreteGrid *SDFFunctions::generateSDF(const unsigned int numVertices,
                                                     const Vector3r *vertices,
                                                     const unsigned int numFaces,
                                                     const unsigned int *faces,
                                                     const AlignedBox3r &bbox,
                                                     const std::array<unsigned int, 3> &resolution,
                                                     const bool invert) {
    START_TIMING("SDF Generation")
    //////////////////////////////////////////////////////////////////////////
    // Generate distance field of object using Discregrid
    //////////////////////////////////////////////////////////////////////////
#ifdef USE_DOUBLE
    TriangleMesh sdfMesh(&vertices[0][0], faces, numVertices, numFaces);
#else
    // if type is float, copy vector to double vector
    std::vector<double> doubleVec;
    doubleVec.resize(3 * numVertices);
    for (unsigned int i = 0; i < numVertices; i++)
        for (unsigned int j = 0; j < 3; j++) doubleVec[3 * i + j] = vertices[i][j];
    TriangleMesh sdfMesh(&doubleVec[0], faces, numVertices, numFaces);
#endif

    TriangleMeshDistance md(sdfMesh);
    Eigen::AlignedBox3d domain;
    domain.extend(bbox.min().cast<double>());
    domain.extend(bbox.max().cast<double>());
    domain.max() += 1.0e-3 * domain.diagonal().norm() * Eigen::Vector3d::Ones();
    domain.min() -= 1.0e-3 * domain.diagonal().norm() * Eigen::Vector3d::Ones();

    auto *distanceField = new CubicLagrangeDiscreteGrid(domain, resolution);
    auto func = DiscreteGrid::ContinuousFunction{};
    Real factor = 1.0;
    if (invert) factor = -1.0;
    func = [&md, &factor](Eigen::Vector3d const &xi) { return factor * md.signed_distance(xi).distance; };

    distanceField->addFunction(func, false);
    STOP_TIMING_PRINT

    return distanceField;
}
}  // namespace vox::utility
