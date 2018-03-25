#include <stdio.h>
#include <iostream>
#include <math.h>
#include "matmul.h"

using namespace std;

matmul::matmul(bool inp){

    flag = inp;
}

Matrix matmul::naive_multiplication(Matrix m1, Matrix m2){
    
    Matrix product(m1.rows, m2.cols, 0);

    for (int i = 0; i < m1.rows; i++)
    {
        for (int j = 0; j < m1.cols; j++)
        {
            for (int k = 0; k < m1.rows; k++)
            {
                product.mat[i][j] += m1.mat[i][k]*m2.mat[k][j];
            }
        }
    }
    
    return product;
}
