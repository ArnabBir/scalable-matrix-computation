#include <stdlib.h>
#include <time.h>
#include <limits>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include "../src/matrix/matrix.cpp"

using namespace std;

int main(){
    
    Matrix  m(10, 10, 0);
    m.display_matrix();

    m.randomize();
    m.display_matrix();

    m = m.transpose();
    m.display_matrix();

    return 0;
}

