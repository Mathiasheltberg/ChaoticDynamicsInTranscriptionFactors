#ifndef CONSTANTHEADER_H
#define CONSTANTHEADER_H
#include <math.h>

////////// Parameters to change //////////////
double Vol = 1.0*pow(10,-14); double nu; double frac; double Ac;
///////////// Simulation Constants /////////////
double Time = 100*pow(10,3); double NA = 6.02*pow(10,23);double Cal = NA*Vol*pow(10,-6);
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
////########## Contants for IKKi %%%%%%%%%%%%%%%                                                                                                                                   
double Hobj = 0;
//////////// Definition of variables ///////////////
long double RT = 0; double RT0 = 0;
int N;int IRna;int I;int IKKa;int IKKi;int IKKn;int Ncy;
double TNF;
int Pp[5];
////////// Script - parameters /////////////
double R [10]; double DT = 2.01;
int Upd [10]; 
int click = 0;int click2 = 0;double dt = 2.0;
int NPeak, IRnaPeak, IPeak, IKKaPeak, IKKiPeak, TimePeak, Nmax;
#endif 
