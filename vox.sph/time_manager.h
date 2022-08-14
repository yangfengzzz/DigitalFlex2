//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#pragma once

#include "vox.sph/binary_file_reader_writer.h"
#include "vox.sph/common.h"

namespace vox {
/** \brief Class to manage the current simulation time and the time step size.
 * This class is a singleton.
 */
class TimeManager {
private:
    Real time;
    static TimeManager *current;
    Real h;

public:
    TimeManager();
    ~TimeManager();

    // Singleton
    static TimeManager *getCurrent();
    static void setCurrent(TimeManager *tm);
    static bool hasCurrent();

    Real getTime() const;
    void setTime(Real t);
    Real getTimeStepSize() const;
    void setTimeStepSize(Real tss);

    void saveState(BinaryFileWriter &binWriter) const;
    void loadState(BinaryFileReader &binReader);
};
}  // namespace vox