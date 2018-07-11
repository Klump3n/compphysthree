cmake_minimum_required(VERSION 2.8)
find_package(CUDA)

# user defined options
set(PHI_0 "0.0" CACHE STRING "Start-Konfiguration (Realteil) (0.0: random)")
set(HFIELD "0.0" CACHE STRING "Magnetfeld (Realteil)")

# Pass options to NVCC
set(CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    -O3 -arch sm_30 -DPHI_0=${PHI_0} -DHFILED=${HFIELD}
    )

# The source files (no 'main')
set(src
   geom_pbc.cu
   geom_pbc.h
   spin.cu
   spin.h
   metropolis.cu
   metropolis.h
   global.h
   randgpu.cu
   randgpu.h
   added_stuff.cu
   added_stuff.h
   added_stuff_gpu.cu
   added_stuff_gpu.h
   )

#The executable (has 'main')
cuda_add_executable(check check.cu ${src})
cuda_add_executable(run run.cu ${src})
