#pragma once

#include <goofit/GlobalCudaDefines.h>

#include <string>
#include <vector>

namespace GooFit {

/// A CPU version of the SmartVector.
///
/// Warning: Changing this vector's size
/// will invalidate the data - while the GPU
/// Version will remain valid (but static)
/// between syncs. This could be corrected,
/// but you should not be touching this
/// in between sync and "device" usage, so it
/// is safe.
template <typename T>
class SmartVector : public std::vector<T> {
    std::string name;

  public:
    SmartVector(std::string name)
        : std::vector<T>()
        , name(name) {}
    /// Sync data to GPU (no op on CPU)
    template <typename Q>
    void sync(Q &device_pointer) {
        device_pointer = this->data();
        GOOFIT_DEBUG("Syncing {}", name);
    }

    /// Only sync changed values
    template <typename Q>
    void smart_sync(Q &device_pointer) {
        sync(device_pointer);
        GOOFIT_DEBUG("Smart syncing {}", name);
    };

    // If this is in a global or static variable,
    // you should clear the device
    // information before the program ends. (no op on CPU)
    void clear_device() { GOOFIT_DEBUG("Clearing {}", name); }
};

} // namespace GooFit
