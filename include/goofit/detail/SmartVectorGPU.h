#pragma once

#include <goofit/GlobalCudaDefines.h>

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include <string>
#include <vector>

namespace GooFit {

/// A GPU version of SmartVector
///
/// This version works by holding two
/// vectors besides the built in one;
/// The `local_copy` stores a copy of "last synced"
/// data and `remote_copy` is on the GPU. This allows
/// `smart_sync` to only copy changed values
/// (verify this is truly faster!)

template <typename T>
class SmartVector : public std::vector<T> {
  protected:
    thrust::host_vector<T> local_copy;
    thrust::device_vector<T> device_copy;
    std::string name;

  public:
    SmartVector(std::string name)
        : std::vector<T>()
        , name(name) {}

    /// Sync data to GPU (no op on CPU)
    template <typename Q>
    void sync(Q &device_pointer) {
        GOOFIT_DEBUG("Syncing {}", name);
        local_copy.resize(this->size());
        std::copy(this->begin(), this->end(), local_copy.begin());
        device_copy = local_copy;

        void *thrust_data_ptr = thrust::raw_pointer_cast(device_copy.data());
        MEMCPY_TO_SYMBOL(device_pointer, &thrust_data_ptr, sizeof(void *), 0, cudaMemcpyHostToDevice);
    }

    /// Only sync changed values
    template <typename Q>
    void smart_sync(Q &device_pointer) {
        if(local_copy.size() != this->size()) {
            sync(device_pointer);
        } else {
            GOOFIT_DEBUG("Smart syncing {}", name);
            for(size_t i = 0; i < this->size(); i++)
                if(local_copy[i] != this->operator[](i))
                    device_copy[i] = local_copy[i] = this->operator[](i);
        }
    }

    // If this is in a global or static variable,
    // you should clear the device
    // information before the program ends. (no op on CPU)
    void clear_device() {
        GOOFIT_DEBUG("Clearing {}", name);
        device_copy.clear();
        device_copy.shrink_to_fit();
    }
};

} // namespace GooFit
