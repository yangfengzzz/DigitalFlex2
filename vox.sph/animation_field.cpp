//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/animation_field.h"

#include <tinyexpr.h>

#include <iostream>
#include <utility>

#include "vox.base/logging.h"
#include "vox.base/time_manager.h"
#include "vox.sph/simulation.h"
#include "vox.sph/time_step.h"

using namespace vox;

AnimationField::AnimationField(std::string particleFieldName,
                               Vector3r pos,
                               Matrix3r rotation,
                               Vector3r scale,
                               const std::string expression[3],
                               const unsigned int type)
    : m_particleFieldName(std::move(particleFieldName)),
      m_x(std::move(pos)),
      m_rotation(std::move(rotation)),
      m_scale(std::move(scale)),
      m_type(type),
      m_startTime(0),
      m_endTime(REAL_MAX) {
    for (int i = 0; i < 3; i++) m_expression[i] = expression[i];
}

AnimationField::~AnimationField() = default;

void AnimationField::reset() {}

double getTime() { return TimeManager::getCurrent()->getTime(); }

void AnimationField::step() {
    Simulation *sim = Simulation::getCurrent();
    TimeManager *tm = TimeManager::getCurrent();
    const Real t = tm->getTime();
    const Real dt = tm->getTimeStepSize();

    if (t >= m_startTime && t <= m_endTime) {
        // animate particles
        const unsigned int nModels = sim->numberOfFluidModels();
        for (unsigned int m = 0; m < nModels; m++) {
            FluidModel *fm = sim->getFluidModel(m);
            const unsigned int numParticles = fm->numActiveParticles();

            // find animated field
            const FieldDescription *particleField = nullptr;
            for (unsigned int j = 0; j < fm->numberOfFields(); j++) {
                const FieldDescription &field = fm->getField(j);
                if (field.name == m_particleFieldName) {
                    particleField = &field;
                    break;
                }
            }

            if (particleField == nullptr) continue;

#pragma omp parallel for schedule(static) default(shared)
            for (int i = 0; i < (int)numParticles; i++) {
                const Vector3r &xi = fm->getPosition(i);
                const Vector3r &vi = fm->getVelocity(i);

                const Eigen::Vector3d xi_double = xi.cast<double>();
                const Eigen::Vector3d vi_double = vi.cast<double>();

                if (inShape(m_type, xi, m_x, m_rotation, m_scale)) {
                    Eigen::Map<Vector3r> value((Real *)particleField->getFct(i));
                    const Eigen::Vector3d value_double = Vector3r(value).cast<double>();

                    const auto t_double = static_cast<double>(t);
                    const auto dt_double = static_cast<double>(dt);

                    te_variable vars[] = {
                            {"t", &t_double},
                            {"dt", &dt_double},
                            {"x", &xi_double[0]},
                            {"y", &xi_double[1]},
                            {"z", &xi_double[2]},
                            {"vx", &vi_double[0]},
                            {"vy", &vi_double[1]},
                            {"vz", &vi_double[2]},
                            {"valuex", &value_double[0]},
                            {"valuey", &value_double[1]},
                            {"valuez", &value_double[2]},
                    };
                    const int numVars = 11;
                    int err;

                    //////////////////////////////////////////////////////////////////////////
                    // v_x
                    //////////////////////////////////////////////////////////////////////////
                    if (!m_expression[0].empty()) {
                        te_expr *expr_vx = te_compile(m_expression[0].c_str(), vars, numVars, &err);
                        if (expr_vx) value[0] = static_cast<Real>(te_eval(expr_vx));
                        te_free(expr_vx);

                        if (err != 0) LOGE("Animation field: expression for x is wrong.")
                    }

                    //////////////////////////////////////////////////////////////////////////
                    // v_y
                    //////////////////////////////////////////////////////////////////////////
                    if (!m_expression[1].empty()) {
                        te_expr *expr_vy = te_compile(m_expression[1].c_str(), vars, numVars, &err);
                        if (expr_vy) value[1] = static_cast<Real>(te_eval(expr_vy));
                        te_free(expr_vy);

                        if (err != 0) LOGE("Animation field: expression for y is wrong.")
                    }

                    //////////////////////////////////////////////////////////////////////////
                    // v_z
                    //////////////////////////////////////////////////////////////////////////
                    if (!m_expression[2].empty()) {
                        te_expr *expr_vz = te_compile(m_expression[2].c_str(), vars, numVars, &err);
                        if (expr_vz) value[2] = static_cast<Real>(te_eval(expr_vz));
                        te_free(expr_vz);

                        if (err != 0) LOGE("Animation field: expression for z is wrong.")
                    }
                }
            }
        }
    }
}
