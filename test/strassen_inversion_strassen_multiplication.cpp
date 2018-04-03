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
    
    Matrix  m(16, 16, 0);
    m.randomize();
    cout<<"Matrix 1 : "<<endl;
    m.display_matrix();

    // Matrix m1 = m.inverse();
    // cout<<"Inverse Matrix : "<<endl;
    // m1.display_matrix();
    
    Matrix m2 = m.strassen_inverse_strassen_multiplication();
    //Matrix m2 = m.strassen_inverse();
    cout<<"Inverse Matrix : "<<endl;
    m2.display_matrix();

   // m1.multiply(m).display_matrix();

    m2.multiply(m).display_matrix();
    return 0;
}

