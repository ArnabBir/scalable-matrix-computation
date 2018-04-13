#include <stdlib.h>
#include <time.h>
#include <limits>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <sys/time.h> 
#include <sys/resource.h>
#include "../src/matrix/matrix.cpp"
#include "../src/matmul/matmul.cpp"

using namespace std;

int main(){
    
    Matrix  m(2048, 2048, 0);
    m.randomize();
    //cout<<"Matrix 1 : "<<endl;
    //m.display_matrix();

    struct rusage usage;

    int start_naive_multiplication = clock();
    m = m.inverse_LU();
    int stop_naive_multiplication = clock();   
    
    getrusage(RUSAGE_SELF, &usage);
    
    cout<<"resident set size = "<<usage.ru_maxrss<<endl;
    cout<<"user time = "<<(usage.ru_utime.tv_sec*1000000.0 + usage.ru_utime.tv_usec)/double(CLOCKS_PER_SEC)*1000<<endl;
    cout<<"sys time = "<<(usage.ru_stime.tv_sec  *1000000.0 + usage.ru_stime.tv_usec)/double(CLOCKS_PER_SEC)*1000<<endl;
    cout<<"Time = "<<(stop_naive_multiplication - start_naive_multiplication)/double(CLOCKS_PER_SEC)*1000<<endl;    
    
    return 0;
}
