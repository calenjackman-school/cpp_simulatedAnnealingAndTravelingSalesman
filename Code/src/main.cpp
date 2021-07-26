//
//  main.cpp
//  TermPaper
//
//  Created by Calen Jackman on 4/16/19.
//  Copyright 2019 Calen Jackman. All rights reserved.
//

#include<iostream>
#include<fstream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<random>
#include<iomanip>
#include<algorithm>
#include<time.h>
#include"Rand.hpp"
using namespace std;
random_device rd;
mt19937 gen(rd());

//Holder Table -abs(sin(x)*cos(y)*exp(abs(1-(sqrt(pow(x,2)+pow(y,2))/M_PI))))
//Bukin N.6 100*sqrt(abs(y-0.01*pow(x,2)))+0.01*abs(x+10)
//sphere function pow(x,2)+pow(y,2)
//booth func pow((x+2*y-7),2)+pow((2*x+y-5),2)
//3 hum camel 2*pow(x,2)-1.05*pow(x,4)+(1.0/6.0)*pow(x,6)+x*y+pow(y,2)
//6 hum camel (4-2.1*pow(x,2)+(1.0/3.0)*pow(x,4))*pow(x,2)+x*y+(-4+4*pow(y,2))*pow(y,2)
//Ackley func -20*exp(-0.2*sqrt(0.5*(x*x+y*y)))-exp(0.5*cos(2*M_PI*x+cos(2*M_PI*y)))+exp(1)+20
//Easom func -cos(x)*cos(y)*exp(-pow((x-M_PI),2)-pow((y-M_PI),2))

double f(double x,double y){
    return pow((x+2*y-7),2)+pow((2*x+y-5),2);
}

double *TwoGrad(double in[2],double beta){
    double *gradf=new double[2];
    gradf[0]=(f(in[0]+beta,in[1])-f(in[0]-beta,in[1]))/(2*beta);
    gradf[1]=(f(in[0],in[1]+beta)-f(in[0],in[1]-beta))/(2*beta);
    return gradf;
}

double *GradDescentConst(double start[2]){
    double *Coord=new double[3];
    double *gvals={nullptr};
    int counter=0;
    double in[2]={0};
    in[0]=start[0];
    in[1]=start[1];
    clock_t t=clock();
    do{
        gvals=TwoGrad(in,0.05);
        //cout<<in[0]<<" "<<in[1]<<" "<<gvals[0]<<" "<<gvals[1]<<endl;
        in[0]-=0.05*gvals[0];
        in[1]-=0.05*gvals[1];
        counter++;
    }while(abs(gvals[0])>0.0001||abs(gvals[1])>0.0001);
    Coord[2]=(double)(clock()-t)/(CLOCKS_PER_SEC/1000);
    Coord[0]=in[0];
    Coord[1]=in[1];
    Coord[3]=counter;
    return Coord;
}

double *GradDescentMC(double start[2]){
    uniform_real_distribution<double>alpha(0,0.1);
    uniform_real_distribution<double>beta(0,0.1);
    double *Coord=new double[4];
    double *gvals={nullptr};
    int counter=0;
    double in[2]={0};
    in[0]=start[0];
    in[1]=start[1];
    double alpha1[200]={0},alpha2[200]={0},beta1[200]={0};
    for(int i=0;i<200;i++){
        alpha1[i]=alpha(gen);
        alpha2[i]=alpha(gen);
        beta1[i]=beta(gen);
    }
    clock_t t=clock();
    do{
        double ax=alpha1[counter];
        double ay=alpha2[counter];
        gvals=TwoGrad(in,beta1[counter]);
        //cout<<in[0]<<" "<<in[1]<<" "<<gvals[0]<<" "<<gvals[1]<<endl;
        in[0]-=ax*gvals[0];
        in[1]-=ay*gvals[1];
        counter++;
    }while(abs(gvals[0])>0.0001||abs(gvals[1])>0.0001);
    Coord[2]=(double)(clock()-t)/(CLOCKS_PER_SEC/1000);
    Coord[0]=in[0];
    Coord[1]=in[1];
    Coord[3]=counter;
    return Coord;
}

double *GradDescentQMC(double start[2]){
    double *Coord=new double[4];
    double *gvals={nullptr};
    int counter=0;
    double in[2]={0};
    in[0]=start[0];
    in[1]=start[1];
    double Halt1[200]={0},Halt2[200]={0},Halt3[200]={0};
    for(int i=0;i<200;i++){
        Halt1[i]=0.1*Halt(i+1,2);
        Halt2[i]=0.1*Halt(i+1,3);;
        Halt3[i]=0.1*Halt(i+1,5);;
    }
    clock_t t=clock();
    do{
        double ax=Halt1[counter];
        double ay=Halt1[counter];
        gvals=TwoGrad(in,Halt1[counter]);
        //cout<<in[0]<<" "<<in[1]<<" "<<gvals[0]<<" "<<gvals[1]<<endl;
        in[0]-=ax*gvals[0];
        in[1]-=ay*gvals[1];
        counter++;
    }while(abs(gvals[0])>0.0001||abs(gvals[1])>0.0001);
    Coord[2]=(double)(clock()-t)/(CLOCKS_PER_SEC/1000);
    Coord[0]=in[0];
    Coord[1]=in[1];
    Coord[3]=counter;
    return Coord;
}

double *GradDescentRQMC(double start[2]){
    uniform_real_distribution<double>rand(0,1);
    double *Coord=new double[4];
    double *gvals={nullptr};
    int counter=0;
    double in[2]={0};
    in[0]=start[0];
    in[1]=start[1];
    double seed1=rand(gen);
    double seed2=rand(gen);
    double seed3=rand(gen);
    double Halt1[200]={0},Halt2[200]={0},Halt3[200]={0};
    for(int i=0;i<200;i++){
        Halt1[i]=0.1*shift(Halt(i+1,2), seed1);
        Halt2[i]=0.1*shift(Halt(i+1,2), seed2);
        Halt3[i]=0.1*shift(Halt(i+1,2), seed3);
    }
    clock_t t=clock();
    do{
        double ax=Halt1[counter];
        double ay=Halt2[counter];
        gvals=TwoGrad(in,Halt3[counter]);
        //cout<<in[0]<<" "<<in[1]<<" "<<gvals[0]<<" "<<gvals[1]<<endl;
        in[0]-=ax*gvals[0];
        in[1]-=ay*gvals[1];
        counter++;
    }while(abs(gvals[0])>0.0001||abs(gvals[1])>0.0001);
    Coord[2]=(double)(clock()-t)/(CLOCKS_PER_SEC/1000);
    Coord[0]=in[0];
    Coord[1]=in[1];
    Coord[3]=counter;
    return Coord;
}

/*double *DoubleTrial(double in[2]){
    uniform_real_distribution<double>alpha(0,0.1);
    double *Coord=new double[3];
    double *gvals={nullptr};
    clock_t t=clock();
    do{
        //double ax=alpha(gen);
        //double ay=alpha(gen);
        //cout<<ax<<" "<<ay<<endl;
        gvals=TwoGrad(in);
        cout<<in[0]<<" "<<in[1]<<" "<<gvals[0]<<" "<<gvals[1]<<endl;
        in[0]-=0.05*gvals[0];
        in[1]-=0.05*gvals[1];
    }while(abs(gvals[0])>0.0001||abs(gvals[1])>0.0001);
    Coord[0]=in[0];
    Coord[1]=in[1];
    Coord[2]=(double)(clock()-t)/(CLOCKS_PER_SEC/1000);
    return Coord;
}*/

double distcalc(double coord1[3], double coord2[3]){
    return sqrt(pow((coord2[1]-coord1[1]),2)+pow((coord2[2]-coord1[2]),2));
}

#define cit 16

double SimAnn(){
    uniform_real_distribution<double>cond(0,1);
    uniform_int_distribution<int>swi(0,cit-1);
    double Halt1[700]={0},Halt2[700]={0};
    double seed1=cond(gen),seed2=cond(gen);
    for(int i=0;i<700;i++){
        Halt1[i]=shift(Halt(i+1,2), seed1);
        Halt2[i]=shift(Halt(i+1,2), seed2);
    }
    int *initialroll={nullptr};
    double points[cit][3];
    initialroll=UniqueSeq(cit,0,cit-1);
    for(int i=0;i<cit;i++){
        points[i][0]=initialroll[i];
    }
    double candidate[cit][3];
    double Temp=5,cool=0.95;
    int divisor=sqrt(cit);
    for(int i=0;i<cit;i++){
        points[i][1]=(int)points[i][0]%divisor;
        points[i][2]=floor(points[i][0]/(double) divisor);
    }
    int counter=0;
    int accept=0;
    double Eprev=0,Efinal=0;
    do{
        for(int i=0;i<cit;i++){
            for(int j=0;j<3;j++){
                candidate[i][j]=points[i][j];
            }
        }
        int pick=Halt1[counter];
        if(pick==0){
            pick=1;
        }
        double hold[3]={0};
        for(int i=0;i<3;i++){
            hold[i]=candidate[pick][i];
        }
        for(int i=0;i<3;i++){
            candidate[pick][i]=candidate[pick-1][i];
            candidate[pick-1][i]=hold[i];
        }
        double Ep=0,Ec=0;
        for(int i=0;i<cit-1;i++){
            Ep+=distcalc(points[i],points[i+1]);
            Ec+=distcalc(candidate[i],candidate[i+1]);
        }
        if(Ec<=Ep){
            for(int i=0;i<cit;i++){
                for(int j=0;j<3;j++){
                    points[i][j]=candidate[i][j];
                }
            }
            Eprev=Ep;
            Ep=Ec;
            accept=0;
        }else{
            double test=Halt2[counter];
            double prob=exp(-(Ec-Ep)/Temp);
            if(test<prob){
                for(int i=0;i<cit;i++){
                    for(int j=0;j<3;j++){
                        points[i][j]=candidate[i][j];
                    }
                }
                Eprev=Ep;
                Ep=Ec;
                accept=0;
            }else{
                accept++;
            }
        }
        Temp*=cool;
        Efinal=Ep;
        counter++;
    }while((accept<10)&&(counter<700));
    return Efinal;
}

int main(){
    uniform_real_distribution<double>uni1(-10,10);
    uniform_real_distribution<double>uni2(-2,2);
    /*double start[2]={0};
    start[0]=uni1(gen);
    start[1]=uni1(gen);
    double *FD1=GradDescentConst(start);
    double *FD2=GradDescentMC(start);
    double *FD3=GradDescentQMC(start);
    double *FD4=GradDescentRQMC(start);
    cout<<"Cons "<<FD1[0]<<" "<<FD1[1]<<" "<<FD1[2]<<" "<<FD1[3]<<endl;
    cout<<"MC "<<FD2[0]<<" "<<FD2[1]<<" "<<FD2[2]<<" "<<FD2[3]<<endl;
    cout<<"QMC "<<FD3[0]<<" "<<FD3[1]<<" "<<FD3[2]<<" "<<FD3[3]<<endl;
    cout<<"RQMC "<<FD4[0]<<" "<<FD4[1]<<" "<<FD4[2]<<" "<<FD4[3]<<endl;*/
    double GDResults[8][50]={0};
    for(int i=0;i<50;i++){
        double start[2]={0};
        start[0]=uni1(gen);
        start[1]=uni1(gen);
        //cout<<"With starting points of: x:"<<start[0]<<" y:"<<start[1]<<endl;
        double *FD1=GradDescentConst(start);
        double *FD2=GradDescentMC(start);
        double *FD3=GradDescentQMC(start);
        double *FD4=GradDescentRQMC(start);
        /*cout<<"Cons "<<FD1[0]<<" "<<FD1[1]<<" "<<FD1[2]<<" "<<FD1[3]<<endl;
        cout<<"MC "<<FD2[0]<<" "<<FD2[1]<<" "<<FD2[2]<<" "<<FD2[3]<<endl;
        cout<<"QMC "<<FD3[0]<<" "<<FD3[1]<<" "<<FD3[2]<<" "<<FD3[3]<<endl;
        cout<<"RQMC "<<FD4[0]<<" "<<FD4[1]<<" "<<FD4[2]<<" "<<FD4[3]<<endl;
        cout<<"\n\n";*/
        GDResults[0][i]=FD1[2];
        GDResults[1][i]=FD1[3];
        GDResults[2][i]=FD2[2];
        GDResults[3][i]=FD2[3];
        GDResults[4][i]=FD3[2];
        GDResults[5][i]=FD3[3];
        GDResults[6][i]=FD4[2];
        GDResults[7][i]=FD4[3];
    }
    ofstream File2;
    File2.open("/Users/ctlj13/Documents/School/MonteCarlo/TermPaper/GDResults.txt", ios::trunc);
    File2<<"Const "<<getMean(GDResults[0],50)<<" "<<getMean(GDResults[1],50)<<" "<<getVariance(GDResults[0],50)<<" "<<getVariance(GDResults[1],50)<<endl;
    File2<<"MC "<<getMean(GDResults[2],50)<<" "<<getMean(GDResults[3],50)<<" "<<getVariance(GDResults[2],50)<<" "<<getVariance(GDResults[3],50)<<endl;
    File2<<"QMC "<<getMean(GDResults[4],50)<<" "<<getMean(GDResults[5],50)<<" "<<getVariance(GDResults[4],50)<<" "<<getVariance(GDResults[5],50)<<endl;
    File2<<"RQMC "<<getMean(GDResults[6],50)<<" "<<getMean(GDResults[7],50)<<" "<<getVariance(GDResults[6],50)<<" "<<getVariance(GDResults[7],50)<<endl;
    /*ofstream File;
    File.open("/Users/ctlj13/Documents/School/MonteCarlo/TermPaper/RNGResults.txt", ios::trunc);
    ofstream File1;
    File1.open("/Users/ctlj13/Documents/School/MonteCarlo/TermPaper/RNGRawResults.txt", ios::trunc);
    double AnnErr[100]={0};
    //testArr[0]=10;
    for(int i=1;i<30;i++){
        testArr[i]=testArr[i-1]+2;
    }
    for(int i=0;i<100;i++){
        AnnErr[i]=SimAnn()-15;
    }
    for (int i=0;i<100;i++){
        File1<<setprecision(3)<<AnnErr[i]<<endl;
    }
    File<<getMean(AnnErr, 100)<<" "<<getVariance(AnnErr, 100)<<endl;
    cout<<getMean(AnnErr, 100)<<" "<<getVariance(AnnErr, 100)<<endl;*/
    return 0;
}
