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

// struct rusage {
//                struct timeval ru_utime; /* user CPU time used */
//                struct timeval ru_stime; /* system CPU time used */
//                long   ru_maxrss;        /* maximum resident set size */
//                long   ru_ixrss;         /* integral shared memory size */
//                long   ru_idrss;         /* integral unshared data size */
//                long   ru_isrss;         /* integral unshared stack size */
//                long   ru_minflt;        /* page reclaims (soft page faults) */
//                long   ru_majflt;        /* page faults (hard page faults) */
//                long   ru_nswap;         /* swaps */
//                long   ru_inblock;       /* block input operations */
//                long   ru_oublock;       /* block output operations */
//                long   ru_msgsnd;        /* IPC messages sent */
//                long   ru_msgrcv;        /* IPC messages received */
//                long   ru_nsignals;      /* signals received */
//                long   ru_nvcsw;         /* voluntary context switches */
//                long   ru_nivcsw;        /* involuntary context switches */
//            };

int main(){

    Matrix  m1(4096, 4096, 0);
    m1.randomize();
    //cout<<"Matrix 1 : "<<endl;
    //m1.display_matrix();

    Matrix  m2(4096, 4096, 0);
    m2.randomize();
    //cout<<"Matrix 2 : "<<endl;
    //m2.display_matrix();

    struct rusage usage;

    int start_naive_multiplication = clock();
    Matrix m = m1.multiply(m2);
    int stop_naive_multiplication = clock();
    
    getrusage(RUSAGE_SELF, &usage);
    cout<<"resident set size = "<<usage.ru_maxrss<<endl;
    cout<<"user time = "<<(usage.ru_utime.tv_sec*1000000.0 + usage.ru_utime.tv_usec)/double(CLOCKS_PER_SEC)*1000<<endl;
    cout<<"sys time = "<<(usage.ru_stime.tv_sec  *1000000.0 + usage.ru_stime.tv_usec)/double(CLOCKS_PER_SEC)*1000<<endl;
    //cout<<"Product Matrix : "<<endl;
    //m.display_matrix();
    
    //cout<<start_naive_multiplication<<" "<<stop_naive_multiplication<<endl;
    cout<<"Time = "<<(stop_naive_multiplication - start_naive_multiplication)/double(CLOCKS_PER_SEC)*1000<<endl;

    return 0;
}

