cmake_minimum_required(VERSION 2.8)
find_package(CUDA)

# Pass options to NVCC
set(CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    -O3 -lm -arch sm_30
    )

# The source files (no 'main')
set(src
   geom_pbc.cu
   geom_pbc.h
   global.h
   )

# The executable (has 'main')
cuda_add_executable(test-complex test.cu ${src})
target_link_libraries(test-complex m)
