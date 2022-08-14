//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.sph/common.h"
#include "vox.sph/fluid_model.h"
#include "vox.sph/viscosity/viscosity_base.h"

namespace vox {
/** \brief This class implements the standard method for viscosity descibed e.g. by
 * Ihmsen et al. [IOS+14].\n\n
 * The method evaluates the term \f$\nu \nabla^2 \mathbf{v}\f$ and uses an approximation
 * of the kernel Laplacian to improve the stability. This approximation is given in [IOS+14].
 *
 * References:
 * - [IOS+14] Markus Ihmsen, Jens Orthmann, Barbara Solenthaler, Andreas Kolb, and Matthias Teschner. SPH Fluids in
 * Computer Graphics. In Sylvain Lefebvre and Michela Spagnuolo, editors, Eurographics 2014 - State of the Art Reports.
 * The Eurographics Association, 2014. URL: http://dx.doi.org/10.2312/egst.20141034
 */
class Viscosity_Standard : public ViscosityBase {
protected:
    Real m_boundaryViscosity;

    void initParameters() override;

public:
    static int VISCOSITY_COEFFICIENT_BOUNDARY;

    explicit Viscosity_Standard(FluidModel* model);
    ~Viscosity_Standard() override;

    void step() override;
    void reset() override;

    static NonPressureForceBase* creator(FluidModel* model) { return new Viscosity_Standard(model); }
};
}  // namespace vox
