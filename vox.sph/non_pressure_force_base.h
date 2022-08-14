//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.base/reflect/parameter_object.h"
#include "vox.sph/common.h"
#include "vox.sph/fluid_model.h"

namespace vox {
/** \brief Base class for all non-pressure force methods.
 */
class NonPressureForceBase : public ParameterObject {
protected:
    FluidModel *m_model;

public:
    explicit NonPressureForceBase(FluidModel *model);
    NonPressureForceBase(const NonPressureForceBase &) = delete;
    NonPressureForceBase &operator=(const NonPressureForceBase &) = delete;
    ~NonPressureForceBase() override;

    virtual void step() = 0;
    virtual void reset(){};

    virtual void performNeighborhoodSearchSort(){};
    virtual void emittedParticles(const unsigned int startIndex){};

    virtual void saveState(BinaryFileWriter &binWriter){};
    virtual void loadState(BinaryFileReader &binReader){};

    FluidModel *getModel() { return m_model; }

    virtual void init();
    /** This function is called after the simulation scene is loaded and all
     * parameters are initialized. While reading a scene file several parameters
     * can change. The deferred init function should initialize all values which
     * depend on these parameters.
     */
    virtual void deferredInit(){};
};
}  // namespace vox