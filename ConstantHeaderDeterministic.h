#ifndef CONSTANTHEADERRK_H
#define CONSTANTHEADERRK_H

#include <math.h>


double Time = 1.0*pow(10,5);
double Vol = 2.0*pow(10,-14);
double NA = 6.02*pow(10,23);
double Cal = NA*Vol*pow(10,-6);
double RT = 0;
double N;double IRna;double I;
double IKKa;
double IKKi;
double G; double B;

////########### Contants for N %%%%%%%%%%%%%%%                                                                     
double kNin = 5.4;double Ntot = 1*Cal;double KN = 0.029*Cal;double kIin = 0.018;
///// ########### Contants for IRNA %%%%%%%%%%%%%%%                                                                
double kt = 1.03/Cal;double ym = 0.017;
///########### Contants for I %%%%%%%%%%%%%%%                                                                      
double ktl = 0.24;double a = 1.05/Cal;double KI = 0.035*Cal;
////########### Contants for IKKa %%%%%%%%%%%%%%                                                                  
double ka = 0.24;double IKKtot = 2*Cal;double ki = 0.18;
////########## Contants for IKKi %%%%%%%%%%%%%%%                                                                   
double kp = 0.036;double kA20 = 0.0018;double A20 = 0.0026;


double IKKn = IKKtot-IKKa-IKKi;double Ncy = Ntot-N; double TNF = 0.5; double nu; double Ac;
int sw; double Nmax; double Nmin; double Tmax; double Tmin;
int swt;
double cdec; double kBp; double kBm; double kGp; double kGm; 
double kTon;double kToff;double kNon;double kNoff; double PNon; double PTon; double PNoff; double PToff; int ONp; int OTp;

double ka1;double ka2;double ka3;double ka4;
double ja1;double ja2;double ja3;double ja4;
double la1;double la2;double la3;double la4;
double ma1;double ma2;double ma3;double ma4;
double na1;double na2;double na3;double na4;
double d1,d2,d3,d4,d5,d6,d7,d8,d9,d10;

double dt = 0.02;double DT = 1.0;double click = 0.0;double it = 0.0;
double t1 = 0;double t2 = 0;double t3 = 0;double t4 = 0;

double Nsta; double IRnasta; double Ista; double IKKasta; double IKKista; 

double inject = 0;double injectT = 0;
double RandT;
double RT0;
double injc;
#endif 
