//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include <utility>
#include <vector>

#include "vox.base/binary_file_reader_writer.h"
#include "vox.base/common.h"
#include "vox.base/reflect/parameter_object.h"
#include "vox.base/sph_kernels.h"
#include "vox.sph/rigid_body_object.h"

namespace vox {
class TimeStep;
class ViscosityBase;
class SurfaceTensionBase;
class VorticityBase;
class DragBase;
class ElasticityBase;
class EmitterSystem;

enum FieldType { Scalar = 0, Vector3, Vector6, Matrix3, Matrix6, UInt };
struct FieldDescription {
    std::string name;
    FieldType type;
    // getFct(particleIndex)
    std::function<void *(const unsigned int)> getFct;
    bool storeData;

    FieldDescription(std::string n,
                     const FieldType &t,
                     std::function<void *(const unsigned int)> fct,
                     const bool s = false)
        : name(std::move(n)), type(t), getFct(std::move(fct)), storeData(s) {}
};

enum class ParticleState { Active = 0, AnimatedByEmitter, Fixed };

/** \brief The fluid model stores the particle and simulation information
 */
class FluidModel : public ParameterObject {
public:
    static int NUM_PARTICLES;
    static int NUM_REUSED_PARTICLES;
    static int DENSITY0;

    static int DRAG_METHOD;
    static int SURFACE_TENSION_METHOD;
    static int VISCOSITY_METHOD;
    static int VORTICITY_METHOD;
    static int ELASTICITY_METHOD;

    FluidModel();
    FluidModel(const FluidModel &) = delete;
    FluidModel &operator=(const FluidModel &) = delete;
    ~FluidModel() override;

    void init();
    /** This function is called after the simulation scene is loaded and all
     * parameters are initialized. While reading a scene file several parameters
     * can change. The deferred init function should initialize all values which
     * depend on these parameters.
     */
    void deferredInit();

    [[nodiscard]] std::string getId() const { return m_id; }

protected:
    std::string m_id;
    EmitterSystem *m_emitterSystem;

    // Mass
    // If the mass is zero, the particle is static
    std::vector<Real> m_masses;
    std::vector<Vector3r> m_a;
    std::vector<Vector3r> m_v0;
    std::vector<Vector3r> m_x0;
    std::vector<Vector3r> m_x;
    std::vector<Vector3r> m_v;
    std::vector<Real> m_density;
    std::vector<unsigned int> m_particleId;
    std::vector<unsigned int> m_objectId;
    std::vector<unsigned int> m_objectId0;
    std::vector<ParticleState> m_particleState;
    Real m_V{};

    unsigned int m_surfaceTensionMethod;
    SurfaceTensionBase *m_surfaceTension;
    unsigned int m_viscosityMethod;
    ViscosityBase *m_viscosity;
    unsigned int m_vorticityMethod;
    VorticityBase *m_vorticity;
    unsigned int m_dragMethod;
    DragBase *m_drag;
    unsigned int m_elasticityMethod;
    ElasticityBase *m_elasticity;
    std::vector<FieldDescription> m_fields;

    std::function<void()> m_dragMethodChanged;
    std::function<void()> m_surfaceTensionMethodChanged;
    std::function<void()> m_viscosityMethodChanged;
    std::function<void()> m_vorticityMethodChanged;
    std::function<void()> m_elasticityMethodChanged;

    Real m_density0;
    unsigned int m_pointSetIndex;

    unsigned int m_numActiveParticles{};
    unsigned int m_numActiveParticles0{};

    void initParameters() override;

    void initMasses();

    /** Resize the arrays containing the particle data.
     */
    virtual void resizeFluidParticles(unsigned int newSize);

    /** Release the arrays containing the particle data.
     */
    virtual void releaseFluidParticles();

public:
    [[nodiscard]] FORCE_INLINE Real getDensity0() const { return m_density0; }
    void setDensity0(Real v);

    [[nodiscard]] unsigned int getPointSetIndex() const { return m_pointSetIndex; }

    void addField(const FieldDescription &field);
    const std::vector<FieldDescription> &getFields() { return m_fields; }
    const FieldDescription &getField(const unsigned int i) { return m_fields[i]; }
    const FieldDescription &getField(const std::string &name);
    unsigned int numberOfFields() { return static_cast<unsigned int>(m_fields.size()); }
    void removeFieldByName(const std::string &fieldName);

    void setNumActiveParticles(unsigned int num);
    [[nodiscard]] unsigned int numberOfParticles() const { return static_cast<unsigned int>(m_x.size()); }

    EmitterSystem *getEmitterSystem() { return m_emitterSystem; }

    virtual void reset();

    void performNeighborhoodSearchSort();

    void initModel(const std::string &id,
                   unsigned int nFluidParticles,
                   Vector3r *fluidParticles,
                   Vector3r *fluidVelocities,
                   const unsigned int *fluidObjectIds,
                   unsigned int nMaxEmitterParticles);

    [[nodiscard]] unsigned int numParticles() const { return static_cast<unsigned int>(m_masses.size()); }
    [[nodiscard]] unsigned int numActiveParticles() const;

    [[nodiscard]] unsigned int getNumActiveParticles0() const { return m_numActiveParticles0; }
    void setNumActiveParticles0(unsigned int val) { m_numActiveParticles0 = val; }

    void emittedParticles(unsigned int startIndex);

    [[nodiscard]] unsigned int getSurfaceTensionMethod() const { return m_surfaceTensionMethod; }
    void setSurfaceTensionMethod(const std::string &val);
    void setSurfaceTensionMethod(unsigned int val);
    [[nodiscard]] unsigned int getViscosityMethod() const { return m_viscosityMethod; }
    void setViscosityMethod(const std::string &val);
    void setViscosityMethod(unsigned int val);
    [[nodiscard]] unsigned int getVorticityMethod() const { return m_vorticityMethod; }
    void setVorticityMethod(const std::string &val);
    void setVorticityMethod(unsigned int val);
    [[nodiscard]] unsigned int getDragMethod() const { return m_dragMethod; }
    void setDragMethod(const std::string &val);
    void setDragMethod(unsigned int val);
    [[nodiscard]] unsigned int getElasticityMethod() const { return m_elasticityMethod; }
    void setElasticityMethod(const std::string &val);
    void setElasticityMethod(unsigned int val);

    SurfaceTensionBase *getSurfaceTensionBase() { return m_surfaceTension; }
    ViscosityBase *getViscosityBase() { return m_viscosity; }
    VorticityBase *getVorticityBase() { return m_vorticity; }
    DragBase *getDragBase() { return m_drag; }
    ElasticityBase *getElasticityBase() { return m_elasticity; }

    void setDragMethodChangedCallback(std::function<void()> const &callBackFct);
    void setSurfaceMethodChangedCallback(std::function<void()> const &callBackFct);
    void setViscosityMethodChangedCallback(std::function<void()> const &callBackFct);
    void setVorticityMethodChangedCallback(std::function<void()> const &callBackFct);
    void setElasticityMethodChangedCallback(std::function<void()> const &callBackFct);

    void computeSurfaceTension();
    void computeViscosity();
    void computeVorticity();
    void computeDragForce();
    void computeElasticity();

    void saveState(BinaryFileWriter &binWriter);
    void loadState(BinaryFileReader &binReader);

    FORCE_INLINE Vector3r &getPosition0(const unsigned int i) { return m_x0[i]; }

    [[nodiscard]] FORCE_INLINE const Vector3r &getPosition0(const unsigned int i) const { return m_x0[i]; }

    FORCE_INLINE void setPosition0(const unsigned int i, const Vector3r &pos) { m_x0[i] = pos; }

    FORCE_INLINE Vector3r &getPosition(const unsigned int i) { return m_x[i]; }

    [[nodiscard]] FORCE_INLINE const Vector3r &getPosition(const unsigned int i) const { return m_x[i]; }

    FORCE_INLINE void setPosition(const unsigned int i, const Vector3r &pos) { m_x[i] = pos; }

    FORCE_INLINE Vector3r &getVelocity(const unsigned int i) { return m_v[i]; }

    [[nodiscard]] FORCE_INLINE const Vector3r &getVelocity(const unsigned int i) const { return m_v[i]; }

    FORCE_INLINE void setVelocity(const unsigned int i, const Vector3r &vel) { m_v[i] = vel; }

    FORCE_INLINE Vector3r &getVelocity0(const unsigned int i) { return m_v0[i]; }

    [[nodiscard]] FORCE_INLINE const Vector3r &getVelocity0(const unsigned int i) const { return m_v0[i]; }

    FORCE_INLINE void setVelocity0(const unsigned int i, const Vector3r &vel) { m_v0[i] = vel; }

    FORCE_INLINE Vector3r &getAcceleration(const unsigned int i) { return m_a[i]; }

    [[nodiscard]] FORCE_INLINE const Vector3r &getAcceleration(const unsigned int i) const { return m_a[i]; }

    FORCE_INLINE void setAcceleration(const unsigned int i, const Vector3r &accel) { m_a[i] = accel; }

    [[nodiscard]] FORCE_INLINE Real getMass(const unsigned int i) const { return m_masses[i]; }

    FORCE_INLINE Real &getMass(const unsigned int i) { return m_masses[i]; }

    FORCE_INLINE void setMass(const unsigned int i, const Real mass) { m_masses[i] = mass; }

    [[nodiscard]] FORCE_INLINE const Real &getDensity(const unsigned int i) const { return m_density[i]; }

    FORCE_INLINE Real &getDensity(const unsigned int i) { return m_density[i]; }

    FORCE_INLINE void setDensity(const unsigned int i, const Real &val) { m_density[i] = val; }

    FORCE_INLINE unsigned int &getParticleId(const unsigned int i) { return m_particleId[i]; }

    [[nodiscard]] FORCE_INLINE const unsigned int &getParticleId(const unsigned int i) const { return m_particleId[i]; }

    FORCE_INLINE unsigned int &getObjectId(const unsigned int i) { return m_objectId[i]; }

    [[nodiscard]] FORCE_INLINE const unsigned int &getObjectId(const unsigned int i) const { return m_objectId[i]; }

    FORCE_INLINE void setObjectId(const unsigned int i, const unsigned int val) { m_objectId[i] = val; }

    [[nodiscard]] FORCE_INLINE const ParticleState &getParticleState(const unsigned int i) const {
        return m_particleState[i];
    }

    FORCE_INLINE ParticleState &getParticleState(const unsigned int i) { return m_particleState[i]; }

    FORCE_INLINE void setParticleState(const unsigned int i, const ParticleState &val) { m_particleState[i] = val; }

    [[nodiscard]] FORCE_INLINE Real getVolume(const unsigned int i) const { return m_V; }

    FORCE_INLINE Real &getVolume(const unsigned int i) { return m_V; }
};
}  // namespace vox
