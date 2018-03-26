/************** LU Decomposition for solving linear equations ***********/
#include<iostream>
#include<cstdlib>
#include "../src/matrix/matrix.cpp"

int main()
{
    int n,i,k,j,p;
    //float m.mat[10][10],l.mat[10][10]={0},u.mat[10][10]={0},sum,b[10],z[10]={0},x[10]={0};
    double sum;
    Matrix  m(10, 10, 0.0), l(10, 10, 0.0), u(10, 10, 0.0);
    m.randomize();
    cout<<"Matrix 1 : "<<endl;
    m.display_matrix();

    n = 10;

    //********** LU decomposition *****//
    for(k=1;k<=n;k++)
    {
        u.mat[k-1][k-1]=1;
        for(i=k;i<=n;i++)
        {
            sum=0;
            for(p=1;p<=k-1;p++)
                sum+=l.mat[i-1][p-1]*u.mat[p-1][k-1];
            l.mat[i-1][k-1]=m.mat[i-1][k-1]-sum;
        }

        for(j=k+1;j<=n;j++)
        {
            sum=0;
            for(p=1;p<=k-1;p++)
                sum+=l.mat[k-1][p-1]*u.mat[p-1][j-1];
            u.mat[k-1][j-1]=(m.mat[k-1][j-1]-sum)/l.mat[k-1][k-1];
        }
    }
    //******** Displaying LU matrix**********//
    cout<<" L ="<<endl;
    l.display_matrix();
    cout<<endl;
    
    cout<<" U ="<<endl;
    u.display_matrix();
    cout<<endl;

    l.multiply(u).display_matrix();


    return 0;
}