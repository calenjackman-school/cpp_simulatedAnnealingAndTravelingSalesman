//
//  Rand.cpp
//  StatTests
//
//  Created by Calen Jackman on 2/18/19.
//  Copyright 2019 Calen Jackman. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <iomanip>

using namespace std;
random_device rd1;
mt19937 gen1(rd1());

// Primes:
// 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199

double *MTUni(int N,int a,int b){
    uniform_real_distribution<double>uni(a,b);
    double *x=new double[N];
    for(int i=0;i<N;i++){
        x[i]=uni(gen1);
    }
    return x;
}

int *MTInt(int N,int a,int b){
    uniform_int_distribution<int>uni(a,b);
    int *x=new int[N];
    for(int i=0;i<N;i++){
        x[i]=uni(gen1);
    }
    return x;
}

int *UniqueSeq(int N,int a,int b){
    if(N>(b-a+1)){
        cout<<"N too large for given constraints. N re-declared as N=(a-b+1).\n";
        N=(b-a+1);
    }
    uniform_int_distribution<int>uni(a,b);
    int *x=new int[N];
    for(int i=0;i<N;i++){
        x[i]=uni(gen1);
    }
    int flag=0;
    do{
        flag=0;
        for(int k=1;k<N;k++){
            for(int l=0;l<k;l++){
                if (x[k]==x[l]){
                    x[k]=uni(gen1);
                    flag++;
                }
            }
        }
    }while(flag>0);
    return x;
}

double Halt(int i, int b){
    double r;
    double f=1;
    r=0;
    do{
        f=f/b;
        r=r+f*(i%b);
        i=i/b;
    } while(i>0);
    return r;
}

double simp(int n){
    int i;
    double a,b,h,sum,x[n],y[n];
    a=0;
    b=1;
    h=(b-a)/n;
    sum=0;
    for(i=0; i<=n; i++){
        x[i]=a+i*h;
        y[i]=exp(x[i]);
    }
    for(i=0; i<=n; i++){
        if (i==0){
            sum+=y[i];
        }
        if (i==n){
            sum+=y[i];
        }
        if (i%2==0){
            sum+=2*y[i];
        }
        else {
            sum+=4*y[i];
        }
    }
    sum*=h/3;
    return sum;
}

double getMean(double data[], int size) {
    double sum=0.0;
    for(int i=0; i<size; i++){
        sum+=data[i];
    }
    return sum/size;
}

double getVariance(double data[], int size) {
    double mean=getMean(data, size);
    double SQD=0;
    for(int i=0; i<size; i++){
        SQD+=pow((data[i]-mean),2);
    }
    return SQD/(size-1);
}

double getStdDev(double data[], int size) {
    return sqrt(getVariance(data, size));
}

double getrealVariance(double data[], double mean, int size) {
    double SQD=0;
    for(int i=0; i<size; i++){
        SQD+=pow((data[i]-mean),2);
    }
    return SQD/(size-1);
}

double getrealStdDev(double data[], double mean, int size) {
    return sqrt(getrealVariance(data, mean, size));
}

double BoxMX(double a, double b){
    double rsq=-2*log(a);
    double theta=2*M_PI*b;
    return sqrt(rsq)*cos(theta);
}

double BoxMY(double a, double b){
    double rsq=-2*log(a);
    double theta=2*M_PI*b;
    return sqrt(rsq)*sin(theta);
}

double BSMval(double u){
    double y, r, x=0;
    double a_0=2.50662823884,
    a_1=-18.61500062529,
    a_2=41.39119773534,
    a_3=-25.44106049637;
    double b_0=-8.47351093090,
    b_1=23.08336743743,
    b_2=-21.06224101826,
    b_3=3.13082909833;
    double c_0=0.3374754822726147,
    c_1=0.9761690190917186,
    c_2=0.1607979714918209,
    c_3=0.0276438810333863,
    c_4=0.0038405729373609,
    c_5=0.0003951896511919,
    c_6=0.0000321767881768,
    c_7=0.0000002888167364,
    c_8=0.0000003960315187;
    y=u-0.5;
    if (fabs(y)<0.42){
        r=pow(y, 2);
        x=y*(((a_3*r+a_2)*r+a_1)*r+a_0)/((((b_3*r+b_2)*r+b_1)*r+b_0)*r+1);
    }else{
        r=u;
        if (y>0){
            r=1-u;
        }
        r=log(-log(r));
        x=c_0+r*(c_1+r*(c_2+r*(c_3+r*(c_4+r*(c_5+r*(c_6+r*(c_7+r*c_8)))))));
        if (y<0){
            x=-x;
        }
    }
    return x;
}

double ADStat(double array[], int size){
    double sum=0;
    for (int i=1; i <= size; i++){
        sum+=((2*i)-1)*(log(array[i-1])+log(1-array[size-i]));
    }
    double A=-size-pow(size, -1)*sum;
    return A;
}

double NormProb(double x){
    return 0.5*erfc(-x/sqrt(2));
}

double shift(double a, double u){
    double sh=a+u;
    if(sh>1){
        sh-=1;
    }
    return sh;
}

double maxim(double a, double b){
    double result;
    if (a <= b){
        result=b;
    }else{
        result=a;
    }
    return result;
}

bool comp(int a, int b){
    return (a<b);
}

int RANDUX(int length, int Mod){
    int x[length];
    for(int i=0; i<length; i++){
        x[i+1]=(65539*x[i])%Mod;
    }
    return *x;
}

double RANDUU(int length, int Mod, int x[]){
    double u[length];
    for(int i=0; i<length; i++){
        u[i]=(double)x[i]/(double)Mod;
    }
    return *u;
}

double RandArr(int length, double rand){
    double x[length];
    for (int i=0; i<length; i++){
        x[i]=rand;
    }
    return *x;
}

double ShHalt(int length, float x[]){
    double sx[length];
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(0.0, 1.0);
    for (int i=0; i<length; i++){
        sx[i]=x[i]+dis(gen);
        if (sx[i]>1){
            sx[i]-=1;
        }
    }
    return *sx;
}

double regressM(double x[], double meanX, double y[], double meanY, int size){
    double sum1, sum2;
    sum1=sum2=0;
    for (int i=0; i<size; i++){
        sum1+=(x[i]-meanX)*(y[i]-meanY);
        sum2+=pow((x[i]-meanX), 2);
    }
    return sum1/sum2;
}

double regressB(double meanX, double meanY, int slope){
    return (meanY-meanX*slope);
}

double cov(double x[], double meanX, double y[], double meanY, int size){
    double sum1;
    sum1=0;
    for (int i=0; i<size; i++){
        sum1+=(x[i]-meanX)*(y[i]-meanY);
    }
    return sum1/((double)size-1);
}
