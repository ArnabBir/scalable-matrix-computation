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
//           Function to create Indentity Matrix             
//-----------------------------------------------------------------

void Matrix::identity(){

    for (int i = 0; i < rows; ++i){

        mat[i][i] = 1.0;
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
//           Function to perform copy constructor of a matrix             
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
    else if(n <= 32){
        for(int i = 0; i < n; ++i){
            for(int j = 0; j < n; ++j){
                for(int k = 0; k < n; ++k){
                    C[i][j] += A[i][k]*B[k][j];
                }
            }
        }
    }
    else{
        int divide_  =(n/2),i,j;
        vector<vector<double> > A11(divide_, vector<double>(divide_));
        vector<vector<double> > A12(divide_, vector<double>(divide_));
        vector<vector<double> > A21(divide_, vector<double>(divide_));
        vector<vector<double> > A22(divide_, vector<double>(divide_));

        vector<vector<double> > B11(divide_, vector<double>(divide_));
        vector<vector<double> > B12(divide_, vector<double>(divide_));
        vector<vector<double> > B21(divide_, vector<double>(divide_));
        vector<vector<double> > B22(divide_, vector<double>(divide_));
        
        vector<vector<double> > C11(divide_, vector<double>(divide_));
        vector<vector<double> > C12(divide_, vector<double>(divide_));
        vector<vector<double> > C21(divide_, vector<double>(divide_));
        vector<vector<double> > C22(divide_, vector<double>(divide_));
        
        vector<vector<double> > P1(divide_, vector<double>(divide_));
        vector<vector<double> > P2(divide_, vector<double>(divide_));
        vector<vector<double> > P3(divide_, vector<double>(divide_));
        vector<vector<double> > P4(divide_, vector<double>(divide_));
        vector<vector<double> > P5(divide_, vector<double>(divide_));
        vector<vector<double> > P6(divide_, vector<double>(divide_));
        vector<vector<double> > P7(divide_, vector<double>(divide_));
        
        vector<vector<double> > AResultant(divide_, vector<double>(divide_));
        vector<vector<double> > BResultant(divide_, vector<double>(divide_));
        
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


//----------------------------------------------------------------------------------
//           Multiply matrices using strassen's algorithm and naive multiplication         
//---------------------------------------------------------------------------------

void inv_strassen(Matrix A, Matrix & C, int n) {

    if(n== 1) {
         C.mat[0][0] = 1.0 / A.mat[0][0];
    }

    //cout<<n<<endl;
    // if(n <= 2){

    //     C = A.inverse();
    // }
    else{


        int divide_  =(n/2),i,j;
        Matrix A11(divide_, divide_, 0.0);
        Matrix A12(divide_, divide_, 0.0);
        Matrix A21(divide_, divide_, 0.0);
        Matrix A22(divide_, divide_, 0.0);

        Matrix C11(divide_, divide_, 0.0);
        Matrix C12(divide_, divide_, 0.0);
        Matrix C21(divide_, divide_, 0.0);
        Matrix C22(divide_, divide_, 0.0);

        Matrix P1(divide_, divide_, 0.0);
        Matrix P2(divide_, divide_, 0.0);
        Matrix P3(divide_, divide_, 0.0);
        Matrix P4(divide_, divide_, 0.0);
        Matrix P5(divide_, divide_, 0.0);
        Matrix P6(divide_, divide_, 0.0);
        Matrix P7(divide_, divide_, 0.0);
        
        for (i = 0; i < divide_; i++)        {
            for (j = 0; j < divide_; j++) {
                A11.mat[i][j] = A.mat[i][j];
                A12.mat[i][j] = A.mat[i][j + divide_];
                A21.mat[i][j] = A.mat[i + divide_][j];
                A22.mat[i][j] = A.mat[i + divide_][j + divide_];

            }
        }
        
        cout<<divide_<<endl;
        inv_strassen(A11, P1, divide_);   // P1 = inv(A11)
        //P1 = A11.strassen_inverse(); // TEST
        P2 = A21.multiply(P1);  // P2 = A21 * PI
        P3 = P1.multiply(A12);  // P3 = P1 * A12
        P4 = A21.multiply(P3);  // P4 = A21 * P3
        substract_matrix(P4.mat, A22.mat, P5.mat, divide_);  // P5 = P4 - A22
        inv_strassen(P5, P6, divide_);  // P6 = inv(P5)
        //P6 = P5.strassen_inverse();
        C12 = P3.multiply(P6);   // C12 = P3 * P6
        C21 = P6.multiply(P2);   // C21 = P4 * P2
        P7 = P3.multiply(C21);   // P7 = P3 * C21
        substract_matrix(P1.mat, P7.mat, C11.mat, divide_);  // C11 = P1 - P7        
        substract_matrix(Matrix(divide_, divide_, 0.0).mat, P6.mat, C22.mat, divide_);  // C22 = - P6

        // C22 = P6;
        // C12 = P3.multiply(P6);
        // C21 = P6.multiply(P2);

        for (i = 0; i < divide_ ; i++)   {
            for (j = 0 ; j < divide_ ; j++)  {
                C.mat[i][j] = C11.mat[i][j];
                C.mat[i][j + divide_] = C12.mat[i][j];
                C.mat[i + divide_][j] = C21.mat[i][j];
                C.mat[i + divide_][j + divide_] = C22.mat[i][j];
            }
        }
    }
}

Matrix Matrix::strassen_inverse(){

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

    inv_strassen( * this, m, count_rows);
    
    return m; 
}

//----------------------------------------------------------------------------------
//           Multiply matrices using strassen's algorithm and strassen's multiplication         
//---------------------------------------------------------------------------------

void inv_strassen_mul_strassen(Matrix A, Matrix & C, int n) {

    if(n== 1) {
         C.mat[0][0] = 1.0 / A.mat[0][0];
    }

    //cout<<n<<endl;
    // if(n <= 2){

    //     C = A.inverse();
    // }
    else{


        int divide_  =(n/2),i,j;
        Matrix A11(divide_, divide_, 0.0);
        Matrix A12(divide_, divide_, 0.0);
        Matrix A21(divide_, divide_, 0.0);
        Matrix A22(divide_, divide_, 0.0);

        Matrix C11(divide_, divide_, 0.0);
        Matrix C12(divide_, divide_, 0.0);
        Matrix C21(divide_, divide_, 0.0);
        Matrix C22(divide_, divide_, 0.0);

        Matrix P1(divide_, divide_, 0.0);
        Matrix P2(divide_, divide_, 0.0);
        Matrix P3(divide_, divide_, 0.0);
        Matrix P4(divide_, divide_, 0.0);
        Matrix P5(divide_, divide_, 0.0);
        Matrix P6(divide_, divide_, 0.0);
        Matrix P7(divide_, divide_, 0.0);
        
        for (i = 0; i < divide_; i++)        {
            for (j = 0; j < divide_; j++) {
                A11.mat[i][j] = A.mat[i][j];
                A12.mat[i][j] = A.mat[i][j + divide_];
                A21.mat[i][j] = A.mat[i + divide_][j];
                A22.mat[i][j] = A.mat[i + divide_][j + divide_];

            }
        }
        
        
        inv_strassen_mul_strassen(A11, P1, divide_);   // P1 = inv(A11)
        //P1 = A11.strassen_inverse(); // TEST
        
        P2 = A21.strassen_multiply(P1);  // P2 = A21 * PI
        P3 = P1.strassen_multiply(A12);  // P3 = P1 * A12
        P4 = A21.strassen_multiply(P3);  // P4 = A21 * P3

        //strassen(A21.mat, P1.mat, P2.mat, divide_);
        //strassen(P1.mat, A12.mat, P3.mat, divide_);
        //strassen(A21.mat, P3.mat, P4.mat, divide_);
        
        //P2 = A21.multiply(P1);  // P2 = A21 * PI
        //P3 = P1.multiply(A12);  // P3 = P1 * A12
        //P4 = A21.multiply(P3);  // P4 = A21 * P3
        
        substract_matrix(P4.mat, A22.mat, P5.mat, divide_);  // P5 = P4 - A22
        inv_strassen_mul_strassen(P5, P6, divide_);  // P6 = inv(P5)
        //P6 = P5.strassen_inverse();
        C12 = P3.strassen_multiply(P6);   // C12 = P3 * P6
        C21 = P6.strassen_multiply(P2);   // C21 = P4 * P2
        P7 = P3.strassen_multiply(C21);   // P7 = P3 * C21

        //strassen(P3.mat, P6.mat, C12.mat, divide_);
        //strassen(P6.mat, P2.mat, C21.mat, divide_);
        //strassen(P3.mat, C21.mat, P7.mat, divide_);

        //C12 = P3.multiply(P6);   // C12 = P3 * P6
        //C21 = P6.multiply(P2);   // C21 = P4 * P2
        //P7 = P3.multiply(C21);   // P7 = P3 * C21

        substract_matrix(P1.mat, P7.mat, C11.mat, divide_);  // C11 = P1 - P7        
        substract_matrix(Matrix(divide_, divide_, 0.0).mat, P6.mat, C22.mat, divide_);  // C22 = - P6

        // C22 = P6;
        // C12 = P3.multiply(P6);
        // C21 = P6.multiply(P2);

        for (i = 0; i < divide_ ; i++)   {
            for (j = 0 ; j < divide_ ; j++)  {
                C.mat[i][j] = C11.mat[i][j];
                C.mat[i][j + divide_] = C12.mat[i][j];
                C.mat[i + divide_][j] = C21.mat[i][j];
                C.mat[i + divide_][j + divide_] = C22.mat[i][j];
            }
        }
    }
}

Matrix Matrix::strassen_inverse_strassen_multiplication(){

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

    inv_strassen_mul_strassen( * this, m, count_rows);
    
    return m; 
}

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

    if( fabs(det) < 0.00001)
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
            inv.mat[i][j] = (fabs(cofdet) < 0.000001 ? 0 : sign) * cofdet / det;
        }
    }
    return inv.transpose();
}

//------------------------------------------------------------------
//           Function to find the Gauss Joardan Elimination              
//------------------------------------------------------------------

void Matrix::gauss_joardan_elimination(){

    double d;
     
    for (int i = 0; i < rows; i++)
    {
        

        for (int j = 0; j < rows; j++)
            if (j != i){
                
                //cout<<mat[i][j]<<" "<<j<<endl;
                d = mat[j][i] / mat[i][i];
                //cout<<d<<endl;
                for (int k = 0; k < cols; k++){

                    mat[j][k] -= mat[i][k] * d;

                }
            }

    }
        
    for (int i = 0; i < rows; i++){

        d = mat[i][i];
        for (int j = 0; j < cols; j++)
            mat[i][j] = mat[i][j] / d;
    }

}

//------------------------------------------------------------------
//           Function to find the Gauss Joardan Inverse              
//------------------------------------------------------------------


Matrix Matrix::gauss_joardan_inverse(){

    Matrix identity(rows, cols, 0.0), inverse(rows, cols, 0.0);
    
    for(int i = 0; i < rows; ++i){
        identity.mat[i][i] = 1.0;
    }

    Matrix joardan_matrix = horzcat(identity);

    joardan_matrix.gauss_joardan_elimination();

    for(int i = 0; i < rows; ++i){

        for(int j = 0; j < cols; ++j){
            inverse.mat[i][j] = joardan_matrix.mat[i][j+cols];
        }
    }

    return inverse;

}

//------------------------------------------------------------------
//           Function for LU decomposition              
//------------------------------------------------------------------

void Matrix::LU(Matrix & L, Matrix & U){

    for (int i = 0; i < rows; i++) {

        for (int j = 0; j < rows; j++){

            if (j < i)
                L.mat[j][i] = 0;
            else{

                L.mat[j][i] = mat[j][i];
                for (int k = 0; k < i; k++){

                    L.mat[j][i] = L.mat[j][i] - L.mat[j][k] * L.mat[k][i];
                }
            }
        }

        for (int j = 0; j < rows; j++){

            if (j < i)
                U.mat[i][j] = 0;
            else if (j == i)
                U.mat[i][j] = 1;
            else{

                U.mat[i][j] = mat[i][j] / L.mat[i][i];
                for (int k = 0; k < i; k++){
                    
                    U.mat[i][j] = U.mat[i][j] - ((L.mat[i][k] * U.mat[k][j]) / L.mat[i][i]);
                }
            }
        }
    }
}

void Matrix::LU_(Matrix & L, Matrix & U){

    double sum;

    for(int k=1;k<= rows;k++){

        U.mat[k-1][k-1]=1;
        for(int i=k;i<=rows;i++){

            sum=0;
            for(int p=1; p<=k-1; p++)
                sum += L.mat[i-1][p-1] * U.mat[p-1][k-1];
            L.mat[i-1][k-1] = mat[i-1][k-1]-sum;
        }

        for(int j=k+1;j<=rows;j++){
            
            sum = 0;
            for(int p=1;p<=k-1;p++)
                sum += L.mat[k-1][p-1] * U.mat[p-1][j-1];
            U.mat[k-1][j-1] = (mat[k-1][j-1] - sum) / L.mat[k-1][k-1];
        }
    }

}

  void Matrix::compute_minor(Matrix & mat_, int d) {

 
    for (int i = 0; i < d; i++)  mat[i][i] = 1.0;
    for (int i = d; i < rows; i++) {
      for (int j = d; j < cols; j++)  {
        mat[i][j] = mat_.mat[i][j];
      }
    }
  }

void Matrix::extract_column(Vector & v, int c) {
  
  if (rows != v.size) {
    std::cerr << "[Matrix::extract_column]: Matrix and Vector sizes don't match\n";
    return;
  }
 
  for (int i = 0; i < rows; i++)
    v(i) = mat[i][c];
}

void vmadd(const Vector & a, const Vector & b, double s, Vector & c){
 
  for (int i = 0; i < c.size; i++)
    c(i) = a(i) + s * b(i);
}

void compute_householder_factor(Matrix & mat_, const Vector & v)
{
 
  int n = v.size;

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      mat_.mat[i][j] = -2 *  v(i) * v(j);

  for (int i = 0; i < n; i++)
    mat_.mat[i][i] += 1;  
}
 
void Matrix::householder(Matrix & R, Matrix & Q)  {
 
  int m = rows;
  int n = cols;
 
  // array of factor Q1, Q2, ... Qm
  vector<Matrix> qv(m);
 
  // temp array
  Matrix z = * this;
  Matrix z1(m, n, 0.0);
 
  for (int k = 0; k < n && k < m - 1; k++) {
 
    Vector e(m), x(m);
    double a;
 
    // compute minor
    //cout<<"hello"<<endl;
    //z.display_matrix();
    cout<<"----------"<<endl;
    z1.compute_minor(z, k);
    //z1.display_matrix();
    
    // extract k-th column into x
    z1.extract_column(x, k);
    
    a = x.norm();
    if (mat[k][k] > 0) a = -a;
    
    
    for (int i = 0; i < e.size; i++)
      e(i) = (i == k) ? 1 : 0;
 
    // e = x + a*e
    vmadd(x, e, a, e);
    
    // e = e / ||e||
    e.rescale_unit();
    
    // qv[k] = I - 2 *e*e^Tq   
    qv[k] = Matrix(rows, cols, 0.0);
    compute_householder_factor(qv[k], e);
    //qv[k].display_matrix();
    // z = qv[k] * z1
    z = qv[k].multiply(z1);
    z.display_matrix();
  }
 
  Q = qv[0];
 
  // after this loop, we will obtain Q (up to a transpose operation)
  for (int i = 1; i < n && i < m - 1; i++) {
 
    z1 = qv[i].multiply(Q);
    Q = z1;
 
  }
 
  R = Q.multiply(* this);
  Q = Q.transpose();
}
