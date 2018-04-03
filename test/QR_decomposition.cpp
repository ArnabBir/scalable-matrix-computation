/************** QR Decomposition for solving linear equations ***********/
#include<iostream>
#include<cstdlib>
#include "../src/matrix/matrix.cpp"

using namespace std;

int main(){
    
    Matrix  m(10, 10, 0.0), R(10, 10, 0.0), Q(10, 10, 0.0);
    m.randomize();
    cout<<"Matrix 1 : "<<endl;
    m.display_matrix();

    m.QR(R, Q);
    cout<<"Decomposed R matrix : "<<endl;
    R.display_matrix();

    cout<<"Decomposed Q matrix : "<<endl;
    Q.display_matrix();

    Matrix product = Q.multiply(R);

    cout<<"Product matrix : "<<endl;
    product.display_matrix();

    return 0;
}