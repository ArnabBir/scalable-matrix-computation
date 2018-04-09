#include <stdlib.h>
#include <time.h>
#include <limits>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include<unistd.h>
#include <sys/time.h> 
#include <sys/resource.h>
#include "../src/matrix/matrix.cpp"
#include "../src/matmul/matmul.cpp"

using namespace std;

int main(){
    
    Matrix  m1(2048, 2048, 0);
    m1.randomize();
    //cout<<"Matrix 1 : "<<endl;
    //m1.display_matrix();

    Matrix  m2(2048, 2048, 0);
    m2.randomize();
    //cout<<"Matrix 2 : "<<endl;
    //m2.display_matrix();

    //matmul multiplier = matmul(true);
    // int start_naive_multiplication = clock();
    // Matrix m = m1.multiply(m2);
    // int stop_naive_multiplication = clock();
    
    //sleep(1);
    // cout<<"Product Matrix : "<<endl;
    // //m.display_matrix();
    // cout<<start_naive_multiplication<<" "<<stop_naive_multiplication<<endl;
    // cout<<"Time = "<<(stop_naive_multiplication - start_naive_multiplication)/double(CLOCKS_PER_SEC)*1000<<endl;
    
    struct rusage usage;

    int start_strassen_multiply = clock();
    Matrix m_ = m1.strassen_multiply(m2);
    int stop_strassen_multiply = clock();
    
    getrusage(RUSAGE_SELF, &usage);
    cout<<"resident set size = "<<usage.ru_maxrss<<endl;
    cout<<"user time = "<<(usage.ru_utime.tv_sec*1000000.0 + usage.ru_utime.tv_usec)/double(CLOCKS_PER_SEC)*1000<<endl;
    cout<<"sys time = "<<(usage.ru_stime.tv_sec  *1000000.0 + usage.ru_stime.tv_usec)/double(CLOCKS_PER_SEC)*1000<<endl;    
    //cout<<"Strassen Product Matrix : "<<endl;
    //m_.display_matrix();
    //cout<<start_strassen_multiply<<" "<<stop_strassen_multiply<<endl;
    cout<<"Time = "<<(stop_strassen_multiply - start_strassen_multiply)/double(CLOCKS_PER_SEC)*1000<<endl;

    return 0;
}

