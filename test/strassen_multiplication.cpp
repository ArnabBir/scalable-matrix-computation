#include <stdlib.h>
#include <time.h>
#include <limits>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include<unistd.h>
#include "../src/matrix/matrix.cpp"
#include "../src/matmul/matmul.cpp"

using namespace std;

int main(){
    
    Matrix  m1(100, 100, 0);
    m1.randomize();
    cout<<"Matrix 1 : "<<endl;
    //m1.display_matrix();

    Matrix  m2(100, 100, 0);
    m2.randomize();
    cout<<"Matrix 2 : "<<endl;
    //m2.display_matrix();

    //matmul multiplier = matmul(true);
    int start_naive_multiplication = clock();
    Matrix m = m1.multiply(m2);
    //sleep(1);
    int stop_naive_multiplication = clock();
    
    //sleep(1);
    cout<<"Product Matrix : "<<endl;
    //m.display_matrix();
    cout<<start_naive_multiplication<<" "<<stop_naive_multiplication<<endl;
    cout<<"Time = "<<(stop_naive_multiplication - start_naive_multiplication)/double(CLOCKS_PER_SEC)*1000<<endl;
    
    //sleep(1);

    int start_strassen_multiply = clock();
    Matrix m_ = m1.strassen_multiply(m2);
    int stop_strassen_multiply = clock();
    
    cout<<"Strassen Product Matrix : "<<endl;
    //m_.display_matrix();
    cout<<start_strassen_multiply<<" "<<stop_strassen_multiply<<endl;
    cout<<"Time = "<<(stop_strassen_multiply - start_strassen_multiply)/double(CLOCKS_PER_SEC)*1000<<endl;
    
    //double a[5][5] = {{1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5},{1, 2, 3, 4, 5},{1, 2, 3, 4, 5}};
    //double b[5][5] = {{1, 1, 1, 1, 1}, {1, 1, 1, 1, 1}, {1, 1, 1, 1, 1},{1, 1, 1, 1, 1},{1, 1, 1, 1, 1}};
    
    //double x[5][5];
    //strassen(a, b, x, 5);
    //for(int i = 0; i < 5; ++i){  for(int j = 0; j < 5; ++j){ cout<<x[i][j]<<'\t';}   cout<<endl;}
    return 0;
}

