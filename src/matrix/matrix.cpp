#include <stdio.h>
#include <iostream>
#include <math.h>
#include <vector>
#include "matrix.h"

using namespace std;

//------------------------------------------------------------------
//           Constructor for the matrix class             
//-----------------------------------------------------------------
Matrix::Matrix(int a=0, int b=0, double val)
{
    rows = a;
    cols = b;
    
    vector<vector <double> > temp_mat; 
    
    for(int i = 0; i < rows; ++i)
    {
        vector<double> v(cols);
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
            printf("%.4lf  ", mat[i][j]);
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
void Matrix::pivot_rearrange(int *leading_0s, Matrix a)
{
    int l, remrow, i, k, lastrow, large;
    double *rowtemn = new double[a.cols];
    lastrow = rows-1;

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
void Matrix::scale_A(int *leading_0s, Matrix a)
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
int* Matrix::row_reduced(Matrix a)
{
    int i, next_row = 1, grp, p, r, j, *leading_0s, t, m, rank;
    leading_0s = new int[a.rows];

    update_leading_0s(leading_0s, a);
    pivot_rearrange(leading_0s, a);

    if(fabs(a.mat[0][0]) < 0.00001)
    {
        cout << "Not a valid matrix as pivot element is 0" << endl;
        return NULL; 
    }

    update_leading_0s(leading_0s, a);
    scale_A(leading_0s, a);

    while(next_row == 1)
    {
        grp = 0;
        for( i = 0; i < rows; ++i)
        {
            p = 0;
            while(leading_0s[i+p] == leading_0s[i+p+1] && (i+p+1) < a.rows)
            {
                grp++;
                p++;
            }

            if(grp != 0)
            {
                while(grp != 0)
                {
                    for(j = 0; j < a.cols; ++j) a.mat[i+grp][j] = a.mat[i+grp][j] - a.mat[i][j];
                    grp--;
                }
                break;
            }
        }   

        update_leading_0s(leading_0s, a);
        pivot_rearrange(leading_0s, a);
        update_leading_0s(leading_0s, a);
        scale_A(leading_0s, a);

        next_row = 0;

        for(r = 0; r < a.rows ; ++r)
        {
            if(leading_0s[r] == leading_0s[r+1] && r+1 < a.rows)
            {
                if(leading_0s[r] != a.cols) next_row = 1;
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
    int *leading_0s = row_reduced(a);
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

vector<vector<double> > add_matrix(vector<vector<double> > A, vector<vector<double> > B, int n) {
    
    int i, j;
    vector<vector<double> > resultant(n, vector<double>(n));
    for(i=0; i< n; i++) {
        for(j=0; j<n; j++)  resultant[i][j] = A[i][j] + B[i][j];
    }

    return resultant;
}

vector<vector<double> > substract_matrix(vector<vector<double> > A, vector<vector<double> > B, int n) {
    
    int i, j;
    vector<vector<double> > resultant(n, vector<double>(n));
    for(i=0; i< n; i++) {
        for(j=0; j<n; j++)  resultant[i][j] = A[i][j] - B[i][j];
    }

    return resultant;
}

int nextPowerOf2(int n)
{
    return pow(2, int(ceil(log2(n))));
}

vector<vector<double> > strassen(vector<vector<double> > A, vector<vector<double> > B, vector<vector<double> > C, int n) {

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
        
        AResultant = add_matrix(A11, A22, divide_);   // a11 + a22
        BResultant = add_matrix(B11, B22, divide_);   // b11 + b22
        P1 = strassen(AResultant, BResultant, P1, divide_);  // p1 = (a11+a22) * (b11+b22)
        
        AResultant = add_matrix(A21, A22, divide_);   // a21 + a22
        P2 = strassen(AResultant, B11, P2, divide_); // p2 = (a21+a22) * (b11)
        
        BResultant = substract_matrix(B12, B22, divide_); // b12 - b22
        P3 = strassen(A11, BResultant, P3, divide_); // p3 = (a11) * (b12 - b22)
        
        BResultant = substract_matrix(B21, B11, divide_); // b21 - b11
        P4 = strassen(A22, BResultant, P4, divide_); // p4 = (a22) * (b21 - b11)
        
        AResultant = add_matrix(A11, A12, divide_); // a11 + a12
        P5 = strassen(AResultant, B22, P5, divide_); // p5 = (a11+a12) * (b22)
        
        AResultant = substract_matrix(A21, A11, divide_); // a21 - a11
        BResultant = add_matrix(B11, B12, divide_); // b11 + b12
        P6 = strassen(AResultant, BResultant, P6, divide_); // p6 = (a21-a11) * (b11+b12)
        
        AResultant = substract_matrix(A12, A22, divide_); // a12 - a22
        BResultant = add_matrix(B21, B22, divide_); // b21 + b22
        P7 = strassen(AResultant, BResultant, P7, divide_); // p7 = (a12-a22) * (b21+b22)
        
        
        C12 = add_matrix(P3, P5, divide_); // c12 = p3 + p5
        C21 = add_matrix(P2, P4, divide_); // c21 = p2 + p4
        
        AResultant = add_matrix(P1, P4, divide_); // p1 + p4
        BResultant = add_matrix(AResultant, P7, divide_); // p1 + p4 + p7
        C11 = substract_matrix(BResultant, P5, divide_); // c11 = p1 + p4 - p5 + p7
        
        AResultant = add_matrix(P1, P3, divide_); // p1 + p3
        BResultant = add_matrix(AResultant, P6, divide_); // p1 + p3 + p6
        C22 = substract_matrix(BResultant, P2, divide_); // c22 = p1 + p3 - p2 + p6

        for (i = 0; i < divide_ ; i++)   {
            for (j = 0 ; j < divide_ ; j++)  {
                C[i][j] = C11[i][j];
                C[i][j + divide_] = C12[i][j];
                C[i + divide_][j] = C21[i][j];
                C[i + divide_][j + divide_] = C22[i][j];
            }
        }
    }

    return C;
}

Matrix Matrix::strassen_multiply(Matrix m2){

    Matrix m(rows, cols, 0.0);
    int rows_ = nextPowerOf2(rows);

    vector<vector<double> > A(rows_, vector<double>(rows_));
    vector<vector<double> > B(rows_, vector<double>(rows_));
    vector<vector<double> > C(rows_, vector<double>(rows_));

    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < cols; ++j){
            A[i][j] = mat[i][j];
            B[i][j] = m2.mat[i][j];
        }
    }

    int dim = rows, count_rows = rows, cont = 0;
    
    if(dim > 1) {
        while(dim>=2) {
            dim/=2;
            cont++;
        }
        
        dim = count_rows;
        if(dim != (pow(2.0,cont))) {
            count_rows = pow(2.0,cont+1);
            for(int i=0; i<count_rows; i++)  {
                for(int j=0; j<count_rows; j++)  {
                    if((i>=dim) || (j>=dim))  {
                        A[i][j] = 0.0;
                        B[i][j] = 0.0;
    }}}}}

    C = strassen(A, B, C, count_rows);

    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < rows; ++j){
            m.mat[i][j] = C[i][j];
        }
    }

    return m; 
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



