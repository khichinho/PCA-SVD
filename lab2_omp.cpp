#include <malloc.h>
#include <omp.h>

#include <math.h>
#include <vector>
using namespace std;

int max_iterations = 30000;

void Matrix_print(vector<vector<float>> &M);
void Matrix_transpose(vector<vector<float>> &M, vector<vector<float>> &M_t);
void Matrix_mult(vector<vector<float>> &A, vector<vector<float>> &B, vector<vector<float>> &result);

void Gram_Schmidt_process(vector<vector<float>> &A, vector<float> &eigen_values, vector<vector<float>> &eigen_vectors);
void Gram_Schmidt_iteration(vector<vector<float>> &A, vector<vector<float>> &Q, vector<vector<float>> &R);

static float D_HAT_ARR[10000000];

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{
    // // Assign Matrix D
    // printf("M = \n");
    vector<vector<float>> M_mat;
    M_mat.assign(M,vector<float>(N));
    for(int i=0; i<M; i++)
    {
        for(int j=0; j<N; j++)
        {
            M_mat[i][j] = D[i*N + j];
        }
    }
    // Matrix_print(M_mat);

    // // Matrix transpose
    // printf("\nM_t = \n");
    static vector<vector<float>> M_transpose;
    M_transpose.assign(N,vector<float>(M));
    Matrix_transpose(M_mat, M_transpose);
    // Matrix_print(M_transpose);



    // // Matrix multiplication
    // printf("\nM_t * M = \n");
    static vector<vector<float>> Mt_M;
    Mt_M.assign(N,vector<float>(N));
    Matrix_mult(M_transpose, M_mat, Mt_M);
    // Matrix_print(Mt_M);


    // // Gram Schmidt begins here:
    // printf("\nGram Schmidt = \n");
    static vector<float> eigen_values;
    eigen_values.assign(N,0);
    static vector<vector<float>> eigen_vectors;
    eigen_vectors.assign(N,vector<float>(N,0));

    Gram_Schmidt_process(Mt_M, eigen_values, eigen_vectors);

    // printf("\nEigen Values = \n");
    // for(int i = 0; i < N; i++)
    // {
    //     printf("%f ", eigen_values[i]);
    // }
    // printf("\n\nEigen Vectors = \n");
    // Matrix_print(eigen_vectors);

    // absolute conversion of eigen values
    for (int i = 0; i < eigen_values.size(); ++i) 
    {
        for (int j = i + 1; j < eigen_values.size(); ++j)
        {
            eigen_values[i] =  fabs(eigen_values[i]);
        }
    }


    // // arrange Eigen Values/Vectors in descending order
    for (int i = 0; i < eigen_values.size(); ++i) 
    {
        for (int j = i + 1; j < eigen_values.size(); ++j)
        {
            if (eigen_values[i] < eigen_values[j] )
            {
                float a =  eigen_values[i];
                eigen_values[i] = eigen_values[j];
                eigen_values[j] = a;

                for(int k = 0; k < N; k++)
                {
                    a =  eigen_vectors[k][i];
                    eigen_vectors[k][i] = eigen_vectors[k][j];
                    eigen_vectors[k][j] = a;
                }
            }
        }
    }

    // // sigma matrix
    static vector<vector<float>> sigma;
    sigma.assign(N,vector<float>(M,0));
    // // sigma inverse matrix
    static vector<vector<float>> sigma_inverse;
    sigma_inverse.assign(N,vector<float>(M,0));

    // // V matrix
    static vector<vector<float>> V;
    V.assign(N,vector<float>(N,0));
    // // V Transpose matrix
    static vector<vector<float>> V_transpose;
    V_transpose.assign(N,vector<float>(N,0));

    for(int i=0; i<N; i++)
    {
        sigma[i][i] = (sqrt(eigen_values[i]));
        sigma_inverse[i][i] = 1/sigma[i][i];
        for (int j = 0; j < N; j++)
        {
            V[i][j] = eigen_vectors[i][j];
            V_transpose[j][i] = eigen_vectors[i][j];
        }
    }

    // printf("\nSigma = \n");
    // Matrix_print(sigma);

    // printf("\nSigma_inverse = \n");
    // Matrix_print(sigma_inverse);

    // printf("\nV = \n");
    // Matrix_print(V);


    // Return Matrices SIGMA, U
    for(int i=0; i<N; i++)
    {
        // TODO: should have returned sigma[i][i]
        *(SIGMA[0] + i) = sigma[i][i];
        for(int j=0; j<N; j++)
        {
            // should have returned U_mat[i][j]
            *( *(U) + i*N + j) = V.at(i).at(j);
        }
    }

    // printf("\nV_t = \n");
    // Matrix_print(V_transpose);

    // U matrix
    static vector<vector<float>> U_mat;
    U_mat.assign(M,vector<float>(N,0));

    // U_mult matrix
    static vector<vector<float>> U_mult;
    U_mult.assign(M,vector<float>(M,0));

    Matrix_mult(M_mat, V, U_mult);
    Matrix_mult(U_mult, sigma_inverse, U_mat);

    // Return Matrices V_T
    for(int i=0; i<M; i++)
    {
        for(int j=0; j<M; j++)
        {
            // should have returned V_transpose[i][j]
            *( *(V_T) + i*M + j) = U_mat[j][i];
        }
    }

    // printf("\nU_mat = \n");
    // Matrix_print(U_mat);
    // printf("\nD = \n");
    // Matrix_print(M_mat);

    // printf("\nSuccess !\n");

    // clear vectors
    M_mat.clear(); 
    M_transpose.clear();
    Mt_M.clear(); eigen_values.clear(); eigen_vectors.clear(); sigma.clear(); sigma_inverse.clear(); V.clear(); V_transpose.clear(); U_mat.clear(); U_mult.clear();


}

void Matrix_print(vector<vector<float>> &M)
{
    // print
    for(int i=0; i<M.size(); i++)
    {
        for(int j=0; j<M[0].size(); j++)
        {
            printf("%f ", M[i][j]);
        }
        printf("\n");
    }
}

void Matrix_transpose(vector<vector<float>> &M, vector<vector<float>> &M_t)
{
    for(int i=0; i<M[0].size(); i++)
    {
        for(int j=0; j<M.size(); j++)
        {
            M_t[i][j] = M[j][i];
        }
    }
}

void Matrix_mult(vector<vector<float>> &A, vector<vector<float>> &B, vector<vector<float>> &result)
{
    for(int i = 0; i < A.size(); ++i)
    {
        for(int j = 0; j < B[0].size(); ++j)
        {
            result[i][j]=0;
        }
    }

    for(int i = 0; i<A.size() ; i++)
    {
        for(int j = 0; j<B[0].size(); j++)
        {
            for(int k = 0; k<A[0].size(); k++)
            {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void Gram_Schmidt_iteration(vector<vector<float>> &A, vector<vector<float>> &Q, vector<vector<float>> &R)
{
    for(int i=0; i < A[0].size(); i++)
    {
        if(i==0)
        {
            // calculate norm
            float norm = 0;
            for(int k = 0; k < A.size(); k++)
            {
                norm += A[k][i]*A[k][i];
            }
            norm = sqrt(norm);

            // Return Q
            for(int k = 0; k < A.size(); k++)
            {
                Q[k][i] = A[k][i]/norm;
            }

            // Return R11
            R[i][i] = norm;
        }

        else{
            // Qi = Ai
            for(int k=0; k < A.size(); k++)
            {
                Q[k][i] = A[k][i];
            }

            // Calculate projections with previous columns
            for(int j = i-1; j > -1; j=j-1){
                float proj = 0;
                for(int k = 0; k < A.size(); k++)
                {
                    proj += Q[k][j] * A[k][i];
                }

                // Return previous R
                R[j][i] = proj;

                for(int k=0; k < A.size(); k++)
                {
                    Q[k][i] -= (proj * Q[k][j]);
                }
            }

            // calculate norm
            float norm = 0;
            for(int k = 0; k < A.size(); k++)
            {
                norm += Q[k][i]*Q[k][i];
            }
            norm = sqrt(norm);

            // Return self R
            R[i][i] = norm;

            // Return Q
            for(int k = 0; k < A.size(); k++)
            {
                Q[k][i] = Q[k][i]/norm;
            }
        }
    }
}

void Gram_Schmidt_process(vector<vector<float>> &A, vector<float> &eigen_values, vector<vector<float>> &eigen_vectors){

    float old_error;
    float new_error;

    vector<vector<float>> D_old;
    D_old.assign(A.size(),vector<float>(A[0].size()));

    for(int j = 0; j < A.size(); j++)
    {
        for(int k = 0; k < A[0].size(); k++)
        {
            D_old[j][k] = A[j][k];
        }
    }

    vector<vector<float>> D_new;
    D_new.assign(A.size(),vector<float>(A[0].size()));

    vector<vector<float>> Q;
    Q.assign(A.size(),vector<float>(A[0].size()));

    vector<vector<float>> R;
    R.assign(A[0].size(),vector<float>(A[0].size()));

    vector<vector<float>> E;
    E.assign(A[0].size(),vector<float>(A[0].size(),0));
    for(int j = 0; j < A.size(); j++)
    {
        E[j][j] = 1.0f;
    }

    vector<vector<float>> E_copy;
    E_copy.assign(A[0].size(),vector<float>(A[0].size(),0));

    for(int i = 0; i < max_iterations; i++)
    {
        // printf("\niteration= %d", i);
        Gram_Schmidt_iteration(D_old, Q, R);

        Matrix_mult(R, Q, D_new);
        // printf("\nD_new = \n");
        // Matrix_print(D_new);

        for(int j = 0; j < A.size(); j++)
        {
            for(int k = 0; k < A[0].size(); k++)
            {
                new_error += (D_new[j][k] - D_old[j][k])*(D_new[j][k] - D_old[j][k]);
                D_old[j][k] = D_new[j][k];

                E_copy[j][k] = E[j][k];
            }
        }

        Matrix_mult(E_copy,Q,E);

        // printf("\nE = \n");
        // Matrix_print(E);

        // printf("\nerror = %f\n", fabs(new_error-old_error));
        if(fabs(new_error-old_error) < 0.000001)
        // if(i == max_iterations - 1)
        {
            for(int j = 0; j < A.size(); j++)
            {
                for(int k = 0; k < A[0].size(); k++)
                {
                    eigen_vectors[j][k] = E[j][k];
                }
                eigen_values[j] = D_new[j][j];
            }
            break;
        }

        old_error = new_error;

    }

    D_old.clear(); D_new.clear(); Q.clear(); R.clear(); E.clear(); E_copy.clear();

}










// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{

    float my_retention = 0;
    float eigen_value_sum = 0;

    // printf("\nSigma = \n");
    for(int i=0; i<N; i++)
    {
        // printf("%f ", SIGMA[i]);
        eigen_value_sum += SIGMA[i];
    }

    // D mat
    // printf("\n\nD = \n");
    static vector<vector<float>> D_mat;
    D_mat.assign(M,vector<float>(N));
    for(int i=0; i<M; i++)
    {
        for(int j=0; j<N; j++)
        {
            D_mat[i][j] = D[i*N + j];
        }
    }
    // Matrix_print(D_mat);

    for(int i=0; i<N; i++){
        my_retention += SIGMA[i]/eigen_value_sum;
        if(my_retention > ((float)retention)/100)
        {
            *K = i+1;
            break;
        }
    }
    // printf("\nK = %d\n", *K);

    // U matrix
    static vector<vector<float>> W_mat;
    W_mat.assign(N,vector<float>(*K,0));

    for(int i=0; i<N; i++)
    {
        for(int j=0; j < (*K); j++)
        {
            W_mat[i][j] = D[i*N + j];
        }
    }

    // static float* D_HAT_ARR = (float*)malloc(sizeof(float) * N*(*K));

    // printf("\nD_HAT = \n");
    for(int i=0; i<M; ++i)
    {
        for(int j=0; j < (*K); ++j)
        {
            D_HAT_ARR[i*N + j] = 0;
            // *(D_HAT_ARR + i*N + j) = 0;
            for(int k = 0; k<N; ++k)
            {
                // *(D_HAT_ARR + i*N + j) += D_mat[i][k] * W_mat[k][j];
                D_HAT_ARR[i*N + j] += D_mat[i][k] * W_mat[k][j];
            }
            // printf("%f ", *(D_HAT_ARR + i*N + j));
            // printf("%f ", D_HAT_ARR[i*N + j]);
        }
            // printf("\n");
    }

    *D_HAT = D_HAT_ARR;

    // clear vectors
    D_mat.clear(); W_mat.clear();
}