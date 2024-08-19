#include <vector>
#include <algorithm>
#include <math.h>
#include <iterator>

#define NDIM 3


inline double dotProduct3(double v1[3], double v2[3])
{
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

inline double dotProduct2(double v1[2], double v2[2])
{
    return v1[0]*v2[0] + v1[1]*v2[1];
}


template <int N>
inline void matrixTranspose( double mat[N][N] )
{
    for (int i=0; i<N; ++i)
        for (int j=i+1; j<N; ++j)
        {
            double temp = mat[i][j];
            mat[i][j] = mat[j][i];
            mat[j][i] = temp;
        }
}


/* Computes the matrix multiplication of two matrices. Matrice sizes are NxN & NxN */
template <int N>
inline void matrixMultiplication(double M1[N][N],
                                 double M2[N][N],
                                 double result[N][N])
{
    for (int i=0; i<N; ++i)
        for (int j=0; j<N; ++j)
        {
            result[i][j] = 0.;
            for (int k=0; k<N; ++k)
                result[i][j] += M1[i][k] * M2[k][j];
        }
}
template <int N>
inline void matrixMultiplication(double M1[N][N],
                                 double *M2,
                                 double *result )
{
    for (int i=0; i<N; ++i)
    {
        result[i] = 0.;
        for (int j=0; j<N; ++j)
            result[i] += M1[i][j] * M2[j];
    }
}



/* Computes the minor determinant coresponding to entry (row,column) of a 3x3 matrix. */
inline double minorDeterminant(double matrix[][3],
                             size_t const row,
                             size_t const column)
{
    size_t i1=(row+1)%3,
    i2=(row+2)%3,
    j1=(column+1)%3,
    j2=(column+2)%3;
    return matrix[j1][i1]*matrix[j2][i2] - matrix[j1][i2]*matrix[j2][i1];
}


/* Computes the inverse of a 2x2 matrix. */
void matrixInverse_2x2(double matrix[][2],
                   double result[][2])
{
    double det = matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0];
    double det_inv = 0.0 ;
    double temp = fabs(matrix[0][0]) + fabs(matrix[0][1]) + fabs(matrix[1][0]) + fabs(matrix[1][1]);
    if( fabs(det) > temp*1e-15 )
        det_inv = 1.0/det;
    result[0][0] = matrix[1][1] *det_inv ;
    result[0][1] = -matrix[0][1] *det_inv ;
    result[1][0] = -matrix[1][0] *det_inv ;
    result[1][1] = matrix[0][0] *det_inv ;
}


/* Computes the inverse of a 3x3 matrix. */
void matrixInverse_3x3(double matrix[][3],
                       double result[][3])
{
    double det =  matrix[0][0] * minorDeterminant( matrix, 0, 0)
                + matrix[0][1] * minorDeterminant( matrix, 1, 0)
                + matrix[0][2] * minorDeterminant( matrix, 2, 0);
    double det_inv = 0.0 ;
    double temp = fabs(matrix[0][0]) + fabs(matrix[0][1]) + fabs(matrix[0][2]) + fabs(matrix[1][0]) + fabs(matrix[1][1]) + fabs(matrix[1][2]) + fabs(matrix[2][0]) + fabs(matrix[2][1]) + fabs(matrix[2][2]);
    if( fabs(det) > temp*1e-15 )
        det_inv = 1.0/det;
    for (size_t row=0; row<3; ++row)
        for (size_t column=0; column<3; ++column)
            result[row][column] = minorDeterminant( matrix, row, column) * det_inv;
}

