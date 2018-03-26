#include <iostream>
#include <vector>

using namespace std;

class Matrix{
public:

    int rows;
    int cols;
    vector<vector<double> > mat;

    Matrix(int, int, double = 0);
    void randomize(void);
    void read_matrix(void);
    void display_matrix(void);
    Matrix copy();
    Matrix transpose();
    int* row_reduced();
    void gaussian_elimination();
    int rank(void);
    double determinant();
    Matrix inverse();
    Matrix gauss_joardan_inverse();
    Matrix multiply(Matrix);
    Matrix readjust();
    Matrix horzcat(Matrix);
    Matrix vertcat(Matrix);
    void scale_A(int *, Matrix &);
    void pivot_rearrange(int *, Matrix &);
    void update_leading_0s(int *, Matrix);
    Matrix strassen_multiply(Matrix);
    void gauss_joardan_elimination();
    void LU(Matrix & L, Matrix & U);


}; 