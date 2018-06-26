cmake_minimum_required(VERSION 2.8)
find_package(CUDA)

# Pass options to NVCC
set(CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    -O3 -arch sm_30
    )

# The source files (no 'main')
set(src
   geom_pbc.cu
   geom_pbc.h
   randgpu.cu
   randgpu.h
   global.h
   )

# The executable (has 'main')
cuda_add_executable(action action.cu ${src})
