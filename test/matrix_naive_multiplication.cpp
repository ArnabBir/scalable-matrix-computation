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
    
    Matrix  m1(10, 10, 0);
    m1.randomize();
    cout<<"Matrix 1 : "<<endl;
    m1.display_matrix();

    Matrix  m2(10, 10, 0);
    m2.randomize();
    cout<<"Matrix 2 : "<<endl;
    m2.display_matrix();

    matmul multiplier = matmul(true);
    
    Matrix m = multiplier.naive_multiplication(m1, m2);
    cout<<"Product Matrix : "<<endl;
    m.display_matrix();

    return 0;
}

