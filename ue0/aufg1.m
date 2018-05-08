mat_1 = rand(1024);
mat_2 = rand(1024);

mat_res = mat_1 * mat_2;

gpu_mat_1 = gpuArray(mat_1);
gpu_mat_2 = gpuArray(mat_2);

gpu_mat_res = gpu_mat_1 * gpu_mat_2;
gpu_cpu_mat_res = gather(gpu_mat_res);

sum(sum(mat_res - gpu_cpu_mat_res))
