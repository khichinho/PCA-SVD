#include <malloc.h>
#include <omp.h>

#include <math.h>

int max_iterations = 30000;

void Matrix_print(float *A, int M, int N);
void Matrix_mult(float *A, float *B, float *result, int M, int N, int K);

void Gram_Schmidt_process(float *A, float *eigen_values, float *eigen_vectors, int N);
void Gram_Schmidt_iteration(float *A, float *Q, float *R, int N);

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{
    // Assign Matrix D and Matrix D_transpose
    float M_mat[M][N];
    float M_transpose[N][M];
    for(int i=0; i<M; i++)
    {
        for(int j=0; j<N; j++)
        {
            M_mat[i][j] = D[i*N + j];
            M_transpose[j][i] = D[i*N + j];
        }
    }
    printf("M = \n");
    Matrix_print((float *)M_mat, M, N);
    printf("\nM_t = \n");
    Matrix_print((float *)M_transpose, N, M);


    // Matrix multiplication
    printf("\nM_t * M = \n");
    float Mt_M[N][N];
    Matrix_mult((float *)M_transpose, (float *)M_mat, (float *)Mt_M, N, M, N);
    Matrix_print((float *)Mt_M, N, N);


    // Gram Schmidt begins here:
    printf("\nGram Schmidt = \n");
    float eigen_values[N][N];
    float eigen_vectors[N][N];

    Gram_Schmidt_process((float *)Mt_M, (float *)eigen_values, (float *)eigen_vectors, N);


    printf("\nEigen Values = \n");
    for(int i = 0; i < N; i++)
    {
        eigen_values[i][i] =  fabs(eigen_values[i][i]);
        // absolute conversion of eigen values
        printf("%f ", eigen_values[i][i]);
    }
    printf("\n\nEigen Vectors = \n");
    Matrix_print((float *)eigen_vectors, N, N);



    // arrange Eigen Values/Vectors in descending order
    for (int i = 0; i < N; ++i) 
    {
        for (int j = i + 1; j < N; ++j)
        {
            if (eigen_values[i][i] < eigen_values[j][j])
            {
                float a =  eigen_values[i][i];
                eigen_values[i][i] = eigen_values[j][j];
                eigen_values[j][i] = a;

                for(int k = 0; k < N; k++)
                {
                    a =  eigen_vectors[k][i];
                    eigen_vectors[k][i] = eigen_vectors[k][j];
                    eigen_vectors[k][j] = a;
                }
            }
        }
    }

    // sigma matrix
    float sigma[N][N];
    // sigma inverse matrix
    float sigma_inverse[N][N];
    // V matrix
    float V[N][N];
    // V Transpose matrix
    float V_transpose[N][N];


    for(int i=0; i<N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            sigma[i][j] = 0.0f;
            sigma_inverse[i][j] = 0.0f;
            V[i][j] = eigen_vectors[i][j];
            V_transpose[j][i] = eigen_vectors[i][j];
        }
        sigma[i][i] = (sqrt(eigen_values[i][i]));
        sigma_inverse[i][i] = 1/(sqrt(eigen_values[i][i]));
    }

    printf("\nSigma = \n");
    Matrix_print((float *)sigma, N, N);
    printf("\nSigma_inverse = \n");
    Matrix_print((float *)sigma_inverse, N, N);
    printf("\nV = \n");
    Matrix_print((float *)V, N, N);
    printf("\nV_t = \n");
    Matrix_print((float *)V_transpose, N, N);

    // U matrix
    float U_mat[M][N];
    // U_mult matrix
    float U_mult[M][N];

    Matrix_mult((float *)M_mat, (float *)V, (float *)U_mult, M, N, N);
    Matrix_mult((float *)U_mult, (float *)sigma_inverse, (float *)U_mat, M, N, N);

    // printf("\nU_mat = \n");
    // Matrix_print((float *)U_mat, M, N);
    // printf("\nD = \n");
    // Matrix_print((float *)M_mat, M, N);

    // Return Matrices U
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            *(*(U)+ i*N + j) = U_mat[i][j];
        }
    }

    // Return Matrices SIGMA, V_T
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            *(SIGMA[0] + i*N + j) = sigma[i][j];
            *(V_T[0] + i*N + j) = V_transpose[i][j];
        }
    }
}

void Matrix_print(float *A, int M, int N)
{
    // print
    for(int i=0; i<M; i++)
    {
        for(int j=0; j<N; j++)
        {
            printf("%f ", *( A + i*N + j));
        }
        printf("\n");
    }
}

void Matrix_mult(float *A, float *B, float *result, int M, int N, int K)
{
    // Initializing elements of matrix mult to 0.
    for(int i = 0; i < M; ++i)
    {
        for(int j = 0; j < K; ++j)
        {
            *(result + i*K + j)=0;
        }
    }

    // Multiplying matrix a and b and storing in array mult.
    for(int i = 0; i < M; ++i)
    {
        for(int j = 0; j < K; ++j)
        {
            for(int k = 0; k < N; ++k)
            {
                *(result + i*K + j) += (*(A + i*N + k)) * (*(B + k*K + j));
            }
        }
    }

}

void Gram_Schmidt_iteration(float *A, float *Q, float *R, int N)
{
    for(int i=0; i < N; i++)
    {
        if(i==0)
        {
            // calculate norm
            float norm = 0;
            for(int k = 0; k < N; k++)
            {
                norm += (*(A + k*N +i))*(*(A + k*N +i));
            }
            norm = sqrt(norm);

            // Return Q
            if(norm == 0){
                for(int k = 0; k < N; k++)
                {
                    (*(Q + k*N +i)) = norm;
                }
            }
            else{
                for(int k = 0; k < N; k++)
                {
                    (*(Q + k*N +i)) = (*(A + k*N +i))/norm;
                }
            }

            // Return R11
            (*(R + i*N +i)) = norm;
        }
        else{
            // Qi = Ai
            for(int k=0; k < N; k++)
            {
                *(Q + k*N +i) = *(A + k*N +i);
            }

            // Calculate projections with previous columns
            for(int j = i-1; j > -1; j=j-1){
                float proj = 0;
                for(int k = 0; k < N; k++)
                {
                    proj += (*(Q + k*N +j)) * (*(A + k*N +i));
                }

                // Return previous R
                (*(R + j*N +i)) = proj;

                for(int k=0; k < N; k++)
                {
                    (*(Q + k*N +i)) -= (proj * (*(Q + k*N +j)));
                }
            }

            // calculate norm
            float norm = 0;
            for(int k = 0; k < N; k++)
            {
                norm += (*(Q + k*N +i)) * (*(Q + k*N +i));
            }
            norm = sqrt(norm);

            // Return self R
            (*(R + i*N +i)) = norm;

            // Return Q
            if(norm == 0){
                for(int k = 0; k < N; k++)
                {
                    (*(Q + k*N +i)) = norm;
                }
            }
            else{
                for(int k = 0; k < N; k++)
                {
                    (*(Q + k*N +i)) = (*(Q + k*N +i))/norm;
                }
            }
        }

    }
}

void Gram_Schmidt_process(float *A, float *eigen_values, float *eigen_vectors, int N)
{
    float old_error;
    float new_error;

    float D_old_ARR[N*N];
    float* D_old = D_old_ARR;

    for(int j = 0; j < 10; j++)
    {
        for(int k = 0; k < N; k++)
        {
            *(D_old + j*N + k) = *(A + j*N + k);
        }
    }

    float D_new_ARR[N*N];
    float* D_new = D_new_ARR;

    float Q_ARR[N*N];
    float* Q = Q_ARR;

    float R_ARR[N*N];
    float* R = R_ARR;

    float E_ARR[N*N];
    float* E = E_ARR;

    float E_copy_ARR[N*N];
    float* E_copy = E_copy_ARR;

    // initialize arrays to stop seg faults in future during multiplication
    for(int j = 0; j < N; j++)
    {
        for(int k = 0; k < N; k++)
        {
            *(Q + j*N + k) = 0;
            *(R + j*N + k) = 0;
            *(E + j*N + k) = 0;
        }
    }

    for(int j = 0; j < N; j++)
    {
        E_ARR[j*N + j] = 1.0f;
    }

    for(int i = 0; i < max_iterations; i++)
    {
        printf("############################## iteration = %d ##############################\n", i+1);
        Gram_Schmidt_iteration(D_old, Q, R, N);

        Matrix_mult(R, Q, D_new, N, N, N);

        for(int j = 0; j < N; j++)
        {
            for(int k = 0; k < N; k++)
            {
                new_error += fabs(D_new[j*N +k] - D_old[j*N +k]);
                D_old[j*N +k] = D_new[j*N +k];

                E_copy[j*N +k] = E[j*N +k];
            }
        }

        Matrix_mult(E_copy, Q, E, N, N, N);

        printf("\nD = \n");
        Matrix_print(D_new, N, N);
        printf("\nE = \n");
        Matrix_print(E, N, N);

        printf("\nerror = %f\n", (fabs(new_error-old_error)) );
        if((fabs(new_error-old_error) <= 0) || (i == max_iterations - 1) )
        // if(i == max_iterations - 1)
        {
            for(int j = 0; j < N; j++)
            {
                for(int k = 0; k < N; k++)
                {
                    *(eigen_vectors + j*N + k) = E[j*N + k];
                }
                *(eigen_values + j*N + j) = D_new[j*N + j];
            }
            break;
        }

        old_error = new_error;
        new_error = 0;
    }
    printf("###########################################################################\n");
}










// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    
}