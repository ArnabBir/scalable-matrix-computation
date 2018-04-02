/************** QR Decomposition for solving linear equations ***********/
#include<iostream>
#include<cstdlib>
#include "../src/matrix/matrix.cpp"

using namespace std;

int main(){
    
    Matrix  m(10, 10, 0.0), L(10, 10, 0.0), U(10, 10, 0.0);
    m.randomize();
    cout<<"Matrix 1 : "<<endl;
    m.display_matrix();

    m.householder(L, U);
    cout<<"Decomposed R matrix : "<<endl;
    L.display_matrix();

    cout<<"Decomposed Q matrix : "<<endl;
    U.display_matrix();

    Matrix product = L.multiply(U);

    cout<<"Product matrix : "<<endl;
    product.display_matrix();

    return 0;
}