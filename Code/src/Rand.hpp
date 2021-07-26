//
//  Rand.hpp
//  RandGen
//
//  Created by Calen Jackman on 2/18/19.
//  Copyright 2019 Calen Jackman. All rights reserved.
//

double Halt(int i, int b);
double simp(int n);
double getMean(double data[], int size);
double getVariance(double data[], int size);
double getStdDev(double data[], int size);
double getrealVariance(double data[], double mean, int size);
double getrealStdDev(double data[], double mean, int size);
int RANDUX(int length, int Mod);
double RANDUU(int length, int Mod, int x[]);
double RandArr(int length, double rand);
double ShHalt(int length, double x[]);
int factorial(int n);
double BoxMX(double a, double b);
double BoxMY(double a, double b);
double BSMval(double u);
double ADStat(double array[], int size);
double NormProb(double x);
double shift(double a, double u);
double maxim(double a, double b);
bool comp(int a, int b);
double cov(double x[], double meanX, double y[], double meanY, int size);
double regressM(double x[], double meanX, double y[], double meanY, int size);
double regressB(double meanX, double meanY, int slope);
double correl(double cov, double sigX, double sigY);
int *UniqueSeq(int N,int a,int b);
int *MTInt(int N,int a,int b);
double *MTUni(int N,int a,int b);
