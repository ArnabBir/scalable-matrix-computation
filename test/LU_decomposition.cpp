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
    
    Matrix  m(10, 10, 0.0), L(10, 10, 0.0), U(10, 10, 0.0);
    m.randomize();
    cout<<"Matrix 1 : "<<endl;
    m.display_matrix();

    m.LU_(L, U);
    cout<<"Decomposed L matrix : "<<endl;
    L.display_matrix();

    cout<<"Decomposed U matrix : "<<endl;
    U.display_matrix();

    Matrix product = L.multiply(U);

    cout<<"Product matrix : "<<endl;
    product.display_matrix();

    return 0;
}

