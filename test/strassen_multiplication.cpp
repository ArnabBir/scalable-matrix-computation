#include <stdlib.h>
#include <time.h>
#include <limits>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "../src/matrix/matrix.cpp"
#include "../src/matmul/matmul.cpp"

using namespace std;

int main(){
    
    Matrix  m1(5, 5, 0);
    m1.randomize();
    cout<<"Matrix 1 : "<<endl;
    m1.display_matrix();

    Matrix  m2(5, 5, 0);
    m2.randomize();
    cout<<"Matrix 2 : "<<endl;
    m2.display_matrix();

    matmul multiplier = matmul(true);
    
    Matrix m = multiplier.naive_multiplication(m1, m2);
    cout<<"Product Matrix : "<<endl;
    m.display_matrix();

    Matrix m_ = m1.strassen_multiply(m2);
    cout<<"Strassen Product Matrix : "<<endl;
    m_.display_matrix();

    //double a[5][5] = {{1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5},{1, 2, 3, 4, 5},{1, 2, 3, 4, 5}};
    //double b[5][5] = {{1, 1, 1, 1, 1}, {1, 1, 1, 1, 1}, {1, 1, 1, 1, 1},{1, 1, 1, 1, 1},{1, 1, 1, 1, 1}};
    
    //double x[5][5];
    //strassen(a, b, x, 5);
    //for(int i = 0; i < 5; ++i){  for(int j = 0; j < 5; ++j){ cout<<x[i][j]<<'\t';}   cout<<endl;}
    return 0;
}

