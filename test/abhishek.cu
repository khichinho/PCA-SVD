#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cmath>
 
// CUDA kernel. Each thread takes care of one element of c
__global__ void vecAdd(double *a, double *b, double *c, int n)
{
    // Get our global thread ID
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    // printf("this is in GPU\n");
 
    // Make sure we do not go out of bounds
    if (id < n){
        c[id] = a[id] + b[id];
        printf("%lf\n",c[id]);
    }
        
}

__global__ void sum(double* a, double* b, double* c){
    int id = blockIdx.x*blockDim.x + threadIdx.x;
    c[id] = a[id]+b[id];
}
__global__ void start(double* a, double* b,int n){
    int id = blockIdx.x*blockDim.x + threadIdx.x;
    // c[id] = a[id]+b[id];
    a[id] = id;
    b[id] = n-id;
}

 
int main( int argc, char* argv[] )
{
    // Size of vectors
    int n = 1000;
    // int m = 4;
    
    // Host input vectors
    double *h_a;
    double *h_b;

    // give_input(&h_a, n, m);
    //Host output vector
    double *h_c;
 
    // Device input vectors
    double *d_a;
    double *d_b;
    //Device output vector
    double *d_c;
 
    // Size, in bytes, of each vector
    size_t bytes = (n*n)*sizeof(double);
    
 
    // Allocate memory for each vector on host
    h_a = (double*)malloc(bytes);
    h_b = (double*)malloc(bytes);

    h_c = (double*)malloc(bytes);
 
    // Allocate memory for each vector on GPU
    cudaMalloc((void **)&d_a, bytes);
    cudaMalloc((void **)&d_b, bytes);
    cudaMalloc((void **)&d_c, bytes);
 
    int i;
    // Initialize vectors on host
    for( i = 0; i < n; i++ ) {
        h_a[i] = i;
		h_b[i] = n-i;
    }
    // start<<<n,n>>>(d_a,d_b,n);
	// for(i=0;i<n/100;i++){
	// 	printf("h_a[ %d ]-> %lf\n",i,h_a[i*100]);
	// 	printf("h_b[ %d ]-> %lf\n",i,h_b[i*100]);
	// }
 
    // Copy host vectors to device
    // cudaMemcpy( d_a, h_a, bytes, cudaMemcpyHostToDevice);
    // cudaMemcpy( d_b, h_b, bytes, cudaMemcpyHostToDevice);
 
    // int blockSize, gridSize;
 
    // Number of threads in each thread block
    // blockSize = 1024;
 
    // Number of thread blocks in grid
    // gridSize = (int)ceil((float)n/blockSize);
 
    // Execute the kernel
    // vecAdd<<<gridSize, blockSize>>>(d_a, d_b, d_c, n);
    // sum<<<n, n>>> (d_a, d_b, d_c);
    // Copy array back to host
    // cudaMemcpy( h_c, d_c, bytes, cudaMemcpyDeviceToHost );
 
    // Sum up vector c and print result divided by n, this should equal 1 within error
	// double sum = 0;
    for(i=0; i<n; i++){
        // sum += h_c[i];
        for(int j=0; j<n; j++){
            h_c[i*n+j] = h_a[i*n+j] + h_b[i*n+j];
            // printf("final hc: %f\n", h_c[i]);
        }
        
    }
        
    // printf("final result: %f\n", sum/n);
 
    // Release device memory
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);
 
    // Release host memory
    free(h_a);
    free(h_b);
    free(h_c);
 
    return 0;
}