#include <stdlib.h>
#include <time.h>
#include <limits>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include<unistd.h>
#include "../src/matrix/matrix.cpp"

using namespace std;

int main(){
    
    Matrix  m1(50, 50, 0);
    m1.randomize();
    cout<<"Matrix 1 : "<<endl;
    m1.display_matrix();
    Matrix  m2(10, 10, 0);
    int * matrix = m1.row_reduced();
    cout<<"Matrix after pivoting : "<<endl;
    m1.display_matrix();;

    return 0;
}

