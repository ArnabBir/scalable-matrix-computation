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
    
    Matrix  m(10, 10, 0);
    m.randomize();
    cout<<"Matrix 1 : "<<endl;
    m.display_matrix();
    
    m = m.inverse();
    cout<<"Inverse Matrix : "<<endl;
    m.display_matrix();

    return 0;
}

