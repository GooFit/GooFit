#if !defined(FAKE_HOST_DEFINES_H)
#define FAKE_HOST_DEFINES_H 1
#include <string>

#define __device__
#define __host__
#define __shared__

enum {
    cudaErrorAddressOfConstant,
    cudaErrorCudartUnloading,
    cudaErrorECCUncorrectable,
    cudaErrorInitializationError,
    cudaErrorInsufficientDriver,
    cudaErrorInvalidChannelDescriptor,
    cudaErrorInvalidConfiguration,
    cudaErrorInvalidDevice,
    cudaErrorInvalidDeviceFunction,
    cudaErrorInvalidDevicePointer,
    cudaErrorInvalidFilterSetting,
    cudaErrorInvalidHostPointer,
    cudaErrorInvalidMemcpyDirection,
    cudaErrorInvalidNormSetting,
    cudaErrorInvalidPitchValue,
    cudaErrorInvalidResourceHandle,
    cudaErrorInvalidSymbol,
    cudaErrorInvalidTexture,
    cudaErrorInvalidTextureBinding,
    cudaErrorInvalidValue,
    cudaErrorLaunchFailure,
    cudaErrorLaunchOutOfResources,
    cudaErrorLaunchTimeout,
    cudaErrorMapBufferObjectFailed,
    cudaErrorMemoryAllocation,
    cudaErrorMemoryValueTooLarge,
    cudaErrorMissingConfiguration,
    cudaErrorMixedDeviceExecution,
    cudaErrorNoDevice,
    cudaErrorNotReady,
    cudaErrorNotYetImplemented,
    cudaErrorPriorLaunchFailure,
    cudaErrorSetOnActiveProcess,
    cudaErrorStartupFailure,
    cudaErrorSynchronizationError,
    cudaErrorTextureFetchFailed,
    cudaErrorTextureNotBound,
    cudaErrorUnknown,
    cudaErrorUnmapBufferObjectFailed,
    cudaMemcpyDeviceToDevice,
    cudaSuccess
};

enum cudaMemcpyKind { cudaMemcpyDeviceToHost, cudaMemcpyHostToDevice };

typedef int cudaError_t;
struct cudaDeviceProp {
    int major;
    int minor;
    int maxGridSize[3];
    int maxThreadsPerBlock;
    int maxThreadsPerMultiProcessor;
    int multiProcessorCount;
    int regsPerBlock;
    int sharedMemPerBlock;
    int warpSize;
};

const int cudaErrorApiFailureBase = 0;

char* cudaGetErrorString(cudaError_t);

cudaError_t cudaMalloc(void**, std::size_t);
cudaError_t cudaGetDeviceProperties(cudaDeviceProp*, std::size_t);
cudaError_t cudaGetDevice(int*);
cudaError_t cudaThreadSynchronize();
cudaError_t cudaMemcpy(void*, const void*, size_t, enum cudaMemcpyKind);
cudaError_t cudaFree(void*);
void cudaMemGetInfo(std::size_t*, std::size_t*);

#endif
