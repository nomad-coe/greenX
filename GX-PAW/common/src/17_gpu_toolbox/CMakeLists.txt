# Obtain files under the folder
if (USEGPU)
    include_directories(common/src/17_gpu_toolbox)
    set(SRC_COMMON "dev_spec.cpp
                    gpu_fft.cpp
                    gpu_linalg.cpp
                    m_gpu_detect.F90
                    m_gpu_toolbox.F90
                    timing.cpp")
    if(NVIDIAGPU)
        set(SRC_VENDOR "dev_spec_cuda.cpp
                        gpu_fft_cuda.cpp
                        gpu_fft.cu
                        gpu_linalg_cuda.cpp
                        m_initcuda.F90
                        m_nvtx.F90
                        timing_cuda.cpp")
    elseif(AMDGPU)
        set(SRC_VENDOR "dev_spec_cuda.cpp
                        gpu_fft_hip.cpp
                        gpu_linalg_hip.cpp
                        timing_hip.cpp")
    else()
        message(FATAL_ERROR "Only AMD and NVIDIA are supported")
    endif()
endif()
set(ABINIT_COMMON_SRC ${ABINIT_COMMON_SRC} ${ABINIT_COMMON_SRC_FOLDER} PARENT_SCOPE)
