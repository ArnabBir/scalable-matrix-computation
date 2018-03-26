#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include "matrix.h"

using namespace std;

int nextPowerOf2(int n)
{
    return pow(2, int(ceil(log2(n))));
}

//------------------------------------------------------------------
//           Constructor for the matrix class             
//-----------------------------------------------------------------
Matrix::Matrix(int a=0, int b=0, double val)
{
    rows = a;
    cols = b;
    int rows_ = nextPowerOf2(rows);
    int cols_ = nextPowerOf2(cols);
    
    vector<vector <double> > temp_mat; 
    
    for(int i = 0; i < rows_; ++i)
    {
        vector<double> v(cols_, val);
        temp_mat.push_back(v);
    }
    
    mat = temp_mat;
}

//------------------------------------------------------------------
//           Function to randomize the values (0, 1)             
//-----------------------------------------------------------------

void Matrix::randomize()
{

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j) mat[i][j] = ((double) rand()) / RAND_MAX;
    }
}

//------------------------------------------------------------------
//           Function to read a matrix             
//-----------------------------------------------------------------
void Matrix::read_matrix()
{
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j<cols; j++)
        {
            cin>>mat[i][j];
        }
    }
}

//------------------------------------------------------------------
//           Function to display a matrix             
//-----------------------------------------------------------------
void Matrix::display_matrix()
{
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j<cols; j++)
        {
            //cout<<mat[i][j]<<"\t";
            printf("% 05.4lf  ", mat[i][j]);
        }
        cout<<endl;
    }
    cout<<endl;
}

//------------------------------------------------------------------
//           Function to perform deep copy of matrix             
//-----------------------------------------------------------------
Matrix Matrix::copy()
{
    int i,j;
    Matrix clone(rows, cols);
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j) clone.mat[i][j] = this->mat[i][j];
    }
    return clone;
}

//------------------------------------------------------------------
//           Function to find transpose of a matrix             
//-----------------------------------------------------------------
Matrix Matrix::transpose()
{
    Matrix result(cols, rows);
    int i, j;
    for(i = 0; i < rows; i++)
    {
        for(j = 0; j < cols; j++)
        {
            result.mat[j][i]  = this->mat[i][j] ;
        }
    }
    return result;
}
//------------------------------------------------------------------
//           Function to concatenate a matrix horizontally            
//-----------------------------------------------------------------
Matrix Matrix::horzcat(Matrix b)
{
    Matrix result(rows, this->cols + b.cols);
    int i, j;
    for(i = 0; i < this->rows; i++)
    {
        for(j = 0; j < this->cols; j++)
        {
            result.mat[i][j] = this->mat[i][j];            
        }
    }
    for(i = 0; i < this->rows; i++)
    {
        for(j = 0; j < b.cols; j++)
        {
            result.mat[i][cols+j]  = b.mat[i][j] ;
        }
    }
    return result;
}

//------------------------------------------------------------------
//           Function to find the first non-zero element              
//-----------------------------------------------------------------
void Matrix::update_leading_0s(int *leading_0s, Matrix a)
{
    int i, j; 
    for(i = 0; i < a.rows; ++i)
    {
        leading_0s[i] = 0;
        for(j=0; (fabs(a.mat[i][j]) < 0.00001) && (j < a.cols) ;++j) leading_0s[i]++;
    }
}


//-------------------------------------------------------------------
// Function to rearrange arrange A such that pivot positions have 1               
//-------------------------------------------------------------------
void Matrix::pivot_rearrange(int *leading_0s, Matrix & a)
{
    int l, remrow, i, k, lastrow, large;
    double *rowtemn = new double[a.cols];
    lastrow = a.rows-1;

    for(l = 0; l < a.rows; ++l)
    {
        large = leading_0s[0];
        for(i = 0; i < a.rows; ++i)
        {
            if( large <= leading_0s[i])
            {
                large=leading_0s[i];
                remrow=i;
            }
        }

    leading_0s[remrow] = leading_0s[lastrow];
    leading_0s[lastrow] = -1;

    for(k = 0; k < a.cols; ++k)   rowtemn[k] = a.mat[lastrow][k];
    for(k = 0; k < a.cols; ++k)   a.mat[lastrow][k] = a.mat[remrow][k];
    for(k = 0; k < a.cols; ++k)   a.mat[remrow][k] = rowtemn[k];
 
    lastrow--;
    }
}

//---------------------------------------------------------------------------
//           Function definition to scale a A                              
//---------------------------------------------------------------------------
void Matrix::scale_A(int *leading_0s, Matrix  & a)
{
    int i,j;
    double divisor;
    for(i = 0; i < a.rows; ++i)
    {  
        if(leading_0s[i] == a.cols) continue;
        divisor = a.mat[i][leading_0s[i]];
        for(j = leading_0s[i]; j < a.cols; ++j)   a.mat[i][j] = a.mat[i][j]/ divisor;
    }

}


//---------------------------------------------------------------------------
//           Function definition to convert to row reduced form                               
//---------------------------------------------------------------------------
int* Matrix::row_reduced()
{
    int i, next_row = 1, grp, p, r, j, *leading_0s, t, m, rank;
    leading_0s = new int[rows];
    update_leading_0s(leading_0s, * this);
    pivot_rearrange(leading_0s, * this);

    if(fabs(mat[0][0]) < 0.00001)
    {
        cout << "Not a valid matrix as pivot element is 0" << endl;
        return NULL; 
    }

    update_leading_0s(leading_0s, * this);
    scale_A(leading_0s, * this);
    
   //for(i = 0; i < a.rows; ++i) cout<<leading_0s[i]<<endl;

    while(next_row == 1)
    {
        
        grp = 0;
        
        for( i = 0; i < rows; ++i)
        {
            p = 0;
            while(leading_0s[i+p] == leading_0s[i+p+1] && (i+p+1) < rows)
            {
                grp++;
                p++;
            }

            if(grp != 0)
            {
                while(grp != 0)
                {
                    for(j = 0; j < cols; ++j) mat[i+grp][j] = mat[i+grp][j] - mat[i][j];
                    grp--;
                }
                break;
            }
        }   

        update_leading_0s(leading_0s, * this);
        pivot_rearrange(leading_0s, * this);
        update_leading_0s(leading_0s, * this);
        scale_A(leading_0s, * this);

        next_row = 0;

        for(r = 0; r < rows ; ++r)
        {
            if(leading_0s[r] == leading_0s[r+1] && r+1 < rows)
            {
                if(leading_0s[r] != cols) next_row = 1;
            }
        }

    }
    return leading_0s;
}

//------------------------------------------------------------------
//           Function to find rank of a matrix             
//-----------------------------------------------------------------
int Matrix::rank()
{
    Matrix a = this->copy();
    int *leading_0s = row_reduced();
    int rank = 0;
    for (int i = 0; i < a.rows; ++i)
    {
        if (leading_0s[i] != a.cols)  ++rank;
    }
    return rank;
}

//------------------------------------------------------------------
//           Function to find the determinant              
//-----------------------------------------------------------------
double Matrix::determinant()
{

    if(rows != cols) { cout<<" Not a square matrix !!"; return 0;}

    int j,p,q;
    double det =0;

    if(rows == 2){ return (mat[0][0]*mat[1][1]) - (mat[0][1]*mat[1][0]); }
    Matrix b(rows-1, rows-1);
    for(j = 0; j < cols; j++)
    {
        int r = 0, s = 0;
        for(p = 0; p < rows; p++)
        {
            for(q = 0; q < cols; q++)
            {
                if(p !=0 && q != j)
                {
                    b.mat[r][s] = mat[p][q];
                    s++;
                    if(s > cols-2)
                        {
                            r++;
                            s = 0;
                        }
                }
            }
        }
        det += (j%2 ? -1:1) * mat[0][j] * b.determinant();
    }
    return det;
}

//------------------------------------------------------------------
//           Function to truncate extremely small float values to 0             
//------------------------------------------------------------------
Matrix Matrix::readjust()
{
    for(int i = 0; i < rows; i++)
    {
        for(int j = 0; j<cols; j++)
        {
            if(fabs(mat[i][j]) < 1e-5) mat[i][j] = 0;
        }
    }
    return *this;
}

//------------------------------------------------------------------
//           Function to concatenate a matrix vertically           
//-----------------------------------------------------------------

Matrix Matrix::vertcat(Matrix b)
{
    Matrix result(this->rows + b.rows, this->cols);
    int i, j;
    for(i = 0; i < this->rows; i++)
    {
        for(j = 0; j < this->cols; j++)
        {
            result.mat[i][j] = this->mat[i][j];            
        }
    }
    for(i = 0; i < b.rows; i++)
    {
        for(j = 0; j < b.cols; j++)
        {
            result.mat[i+this->rows][j]  = b.mat[i][j] ;
        }
    }
    return result;
}

//------------------------------------------------------------------
//           Function to naive multiplication of 2 matrices             
//-----------------------------------------------------------------

Matrix Matrix::multiply(Matrix other)
{
    if(cols != other.rows) { cout<< " Invalid dimensions !"; return *this;}

    Matrix product(this-> rows, other.cols);
    int i, j, k;

    for(i = 0; i < this-> rows; i++)
    {
        for(j = 0; j < other.cols; j++)
        {
            for(k = 0; k < this-> cols; k++)
            {
                product.mat[i][j] += this->mat[i][k] * other.mat[k][j];
            }
        }
    }
    return product;
}

//------------------------------------------------------------------
//           Multiply matrices using strassen's algorithm          
//------------------------------------------------------------------

void add_matrix(vector<vector<double> > A, vector<vector<double> > B, vector<vector<double> > & resultant, int n) {
    
    int i, j;
    for(i=0; i< n; ++i) {
        for(j=0; j<n; ++j)  resultant[i][j] = A[i][j] + B[i][j];
    }

}

void substract_matrix(vector<vector<double> > A, vector<vector<double> > B, vector<vector<double> > & resultant, int n) {
    
    int i, j;
    for(i=0; i< n; ++i) {
        for(j=0; j<n; ++j)  resultant[i][j] = A[i][j] - B[i][j];
    }

}

void strassen(vector<vector<double> > A, vector<vector<double> > B, vector<vector<double> > & C, int n) {

    if(n== 1) {
        C[0][0] = A[0][0] * B[0][0];
    }
    else{
        int divide_  =(n/2),i,j;
        vector<vector<double> > A11(n, vector<double>(n));
        vector<vector<double> > A12(n, vector<double>(n));
        vector<vector<double> > A21(n, vector<double>(n));
        vector<vector<double> > A22(n, vector<double>(n));
        vector<vector<double> > B11(n, vector<double>(n));
        vector<vector<double> > B12(n, vector<double>(n));
        vector<vector<double> > B21(n, vector<double>(n));
        vector<vector<double> > B22(n, vector<double>(n));
        vector<vector<double> > C11(n, vector<double>(n));
        vector<vector<double> > C12(n, vector<double>(n));
        vector<vector<double> > C21(n, vector<double>(n));
        vector<vector<double> > C22(n, vector<double>(n));
        vector<vector<double> > P1(n, vector<double>(n));
        vector<vector<double> > P2(n, vector<double>(n));
        vector<vector<double> > P3(n, vector<double>(n));
        vector<vector<double> > P4(n, vector<double>(n));
        vector<vector<double> > P5(n, vector<double>(n));
        vector<vector<double> > P6(n, vector<double>(n));
        vector<vector<double> > P7(n, vector<double>(n));
        vector<vector<double> > AResultant(n, vector<double>(n));
        vector<vector<double> > BResultant(n, vector<double>(n));
        
        for (i = 0; i < divide_; i++)        {
            for (j = 0; j < divide_; j++) {
                A11[i][j] = A[i][j];
                A12[i][j] = A[i][j + divide_];
                A21[i][j] = A[i + divide_][j];
                A22[i][j] = A[i + divide_][j + divide_];
                
                B11[i][j] = B[i][j];
                B12[i][j] = B[i][j + divide_];
                B21[i][j] = B[i + divide_][j];
                B22[i][j] = B[i + divide_][j + divide_];
            }
        }
        
        add_matrix(A11, A22, AResultant, divide_);   // A11 + A22
        add_matrix(B11, B22, BResultant, divide_);   // B11 + B22
        strassen(AResultant, BResultant, P1, divide_);  // P1 = (A11+A22) * (B11+B22)
        
        add_matrix(A21, A22, AResultant, divide_);   // A21 + A22
        strassen(AResultant, B11, P2, divide_); // P2 = (A21+A22) * (B11)
        
        substract_matrix(B12, B22, BResultant, divide_); // B12 - B22
        strassen(A11, BResultant, P3, divide_); // P3 = (A11) * (B12 - B22)
        
        substract_matrix(B21, B11, BResultant, divide_); // B21 - B11
        strassen(A22, BResultant, P4, divide_); // P4 = (A22) * (B21 - B11)
        
        add_matrix(A11, A12, AResultant, divide_); // A11 + A12
        strassen(AResultant, B22, P5, divide_); // P5 = (A11+A12) * (B22)
        
        substract_matrix(A21, A11, AResultant, divide_); // A21 - A11
        add_matrix(B11, B12, BResultant, divide_); // B11 + B12
        strassen(AResultant, BResultant, P6, divide_); // p6 = (A21-A11) * (B11+B12)
        
        substract_matrix(A12, A22, AResultant, divide_); // A12 - A22
        add_matrix(B21, B22, BResultant, divide_); // B21 + B22
        strassen(AResultant, BResultant, P7, divide_); // P7 = (A12-A22) * (B21+B22)
        
        
        add_matrix(P3, P5, C12, divide_); // C12 = P3 + P5
        add_matrix(P2, P4, C21, divide_); // C21 = P2 + P4
        
        add_matrix(P1, P4, AResultant, divide_); // P1 + P4
        add_matrix(AResultant, P7, BResultant, divide_); // P1 + P4 + P7
        substract_matrix(BResultant, P5, C11, divide_); // C11 = P1 + P4 - P5 + P7
        
        add_matrix(P1, P3, AResultant, divide_); // P1 + P3
        add_matrix(AResultant, P6, BResultant, divide_); // P1 + P3 + P6
        substract_matrix(BResultant, P2, C22, divide_); // C22 = P1 + P3 - P2 + P6

        for (i = 0; i < divide_ ; i++)   {
            for (j = 0 ; j < divide_ ; j++)  {
                C[i][j] = C11[i][j];
                C[i][j + divide_] = C12[i][j];
                C[i + divide_][j] = C21[i][j];
                C[i + divide_][j + divide_] = C22[i][j];
            }
        }
    }
}

Matrix Matrix::strassen_multiply(Matrix m2){

    Matrix m(rows, cols, 0.0);
    int rows_ = nextPowerOf2(rows);

    int dim = rows, count_rows = rows, cont = 0;
    
    if(dim > 1) {
         while(dim>=2) {
             dim/=2;
             cont++;
        }
        
        dim = count_rows;
        if(dim != (pow(2.0,cont))) {
            count_rows = pow(2.0,cont+1);
    }}

    strassen(mat, m2.mat, m.mat, count_rows);
    
    return m; 
}

// Matrix Matrix::strassen_multiply(Matrix m2){

//     Matrix m(rows, cols, 0.0);
//     int rows_ = nextPowerOf2(rows);

//     vector<vector<double> > A(rows_, vector<double>(rows_));
//     vector<vector<double> > B(rows_, vector<double>(rows_));
//     vector<vector<double> > C(rows_, vector<double>(rows_));

//     for(int i = 0; i < rows; ++i){
//         for(int j = 0; j < cols; ++j){
//             A[i][j] = mat[i][j];
//             B[i][j] = m2.mat[i][j];
//         }
//     }

//     int dim = rows, count_rows = rows, cont = 0;
    
//     if(dim > 1) {
//         while(dim>=2) {
//             dim/=2;
//             cont++;
//         }
        
//         dim = count_rows;
//         if(dim != (pow(2.0,cont))) {
//             count_rows = pow(2.0,cont+1);
//             for(int i=0; i<count_rows; i++)  {
//                 for(int j=0; j<count_rows; j++)  {
//                     if((i>=dim) || (j>=dim))  {
//                         A[i][j] = 0.0;
//                         B[i][j] = 0.0;
//     }}}}}

//     strassen(A, B, C, count_rows);

//     for(int i = 0; i < rows; ++i){
//         for(int j = 0; j < rows; ++j){
//             m.mat[i][j] = C[i][j];
//         }
//     }

//     return m; 
// }

//------------------------------------------------------------------
//           Function to inpmelent Gaussian Elimination              
//------------------------------------------------------------------

void Matrix::gaussian_elimination(){

    row_reduced();
}


//------------------------------------------------------------------
//           Function to find the naive inverse              
//------------------------------------------------------------------

Matrix Matrix::inverse()
{
    double det = this->determinant();

    if( fabs(det) < 0.01)
    {
        cout << "Matrix is not invertible!" << endl;
    }

    int i,j, q, p, sign, r, s;
    double cofdet;
    Matrix inv( rows, cols);
    Matrix cof( rows - 1, cols - 1);

    for(i = 0; i < rows; i++)
    {
        for(j = 0; j < cols; j++)
        {
            sign = ((i+j)%2 ? -1 : 1);
            r = 0, s = 0;
            for(p = 0; p < rows; p++)
            {
                for(q = 0; q < cols; q++)
                {
                    if(p != i && q != j)
                    {
                        cof.mat[r][s] = mat[p][q];
                        s++;
                        if(s > cols-2)
                            {
                                r++;
                                s = 0;
                            }
                    }
                }
            }
            cofdet = cof.determinant();
            inv.mat[i][j] = (fabs(cofdet) < 0.1 ? 0 : sign) * cofdet / det;
        }
    }
    return inv.transpose();
}



