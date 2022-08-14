//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include <pybind11/pybind11.h>

#include <iostream>

#include "py.sph/common.h"
#include "vox.editor/exporter/exporter_base.h"

namespace py = pybind11;

class pyExporterBase : public vox::ExporterBase {
    using vox::ExporterBase::ExporterBase;

    void init(const std::string& outputPath) override { PYBIND11_OVERRIDE(void, vox::ExporterBase, init, outputPath); };
    void reset() override { PYBIND11_OVERRIDE(void, vox::ExporterBase, reset, ); };
    void step(const unsigned int frame) override { PYBIND11_OVERRIDE_PURE(void, vox::ExporterBase, step, frame); };
    void setActive(const bool active) override { PYBIND11_OVERRIDE(void, vox::ExporterBase, setActive, active); }
};

void ExporterModule(const py::module& m_sub) {
    // ---------------------------------------
    // Exporter Base
    // ---------------------------------------
    py::class_<vox::ExporterBase, pyExporterBase>(m_sub, "ExporterBase")
            .def(py::init<vox::SimulatorBase*>())
            .def("step", &vox::ExporterBase::step)
            .def("reset", &vox::ExporterBase::reset)
            .def("setActive", &vox::ExporterBase::setActive)
            .def("getActive", &vox::ExporterBase::getActive)
            .def("init", &vox::ExporterBase::init);
}
