#include <malloc.h>
#include <omp.h>


// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{
    // Matrix 1D t array
    printf("M = \n");
    float D_mat[M][N];
    float D_transpose[N][M];
    for(int i=0; i<M; i++)
    {
        for(int j=0; j<N; j++)
        {
            D_mat[i][j] = D[i*N + j];
            printf("%f ", D_mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");


    // Matrix transpose
    printf("M_t = \n");
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<M; j++)
        {
            D_transpose[i][j] = D_mat[j][i];
            printf("%f ", D_transpose[i][j]);
        }
        printf("\n");
    }
    printf("\n");


    // Matrix multiplication
    printf("M_t * M = \n");
    float sum=0;
    float Mt_M[N][N];
    for(int i = 0; i<N ; i++)
    {
        for(int j = 0; j<M; j++)
        {
            for(int k = 0; k<M; k++)
            {
                sum+=D_transpose[i][k] * D_mat[k][j];
            }
            Mt_M[i][j] = sum;
            sum=0;
        }
    }
    // print
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            printf("%f ", Mt_M[i][j]);
        }
        printf("\n");
    }


    // Gram Schmidt begins here:
    // for(int i = 0; i=n ; i++)
    // {
    //     /* code */
    // }




}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    
}
