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
    
    Matrix  m1(5, 6, 0);
    m1.randomize();
    cout<<"Matrix 1 : "<<endl;
    m1.display_matrix();
    Matrix  m2(10, 10, 0);
    m1.gaussian_elimination();
    cout<<"Matrix after Gaussian Elimination : "<<endl;
    m1.display_matrix();;

    return 0;
}

