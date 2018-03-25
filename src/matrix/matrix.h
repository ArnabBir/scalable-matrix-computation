class Matrix{
public:

    int rows;
    int cols;
    double** mat;

    Matrix(int, int, double = 0);
    void randomize(void);
    void read_matrix(void);
    void display_matrix(void);
    Matrix copy();
    Matrix transpose();
    int* row_reduced(Matrix);
    int rank(void);
    double determinant();
    Matrix inverse();
    Matrix multiply(Matrix);
    Matrix readjust();
    Matrix horzcat(Matrix);
    Matrix vertcat(Matrix);
    void scale_A(int *, Matrix);
    void pivot_rearrange(int *, Matrix);
    void update_leading_0s(int *, Matrix);
    Matrix strassen_multiply(Matrix);

}; 