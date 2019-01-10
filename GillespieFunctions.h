#include "ConstantHeaderGillespie.h"
#include <math.h> 
#include <stdlib.h> 
#include <time.h>
#include <fstream>
#include <ctype.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

double GilUpdTime(int Part[], double s6, double randT) {
  Ncy = Ntot-Part[0];   IKKn = IKKtot-Part[3]-Part[4];

  R[0] = kNin*Ncy*KI/(KI+Part[2]);      R[1] = kIin*Part[2]*Part[0]/(KN+Part[0]);
  R[2] = kt*pow(Part[0],2);     R[3] = ym*Part[1];
  R[4] = ktl*Part[1];        R[5] = a*Part[3]*Ncy*Part[2]/(KI+Part[2]);
  R[6] = ka*s6*IKKn;     R[7] = ki*Part[3];
  R[8] = ki*Part[3];         R[9] = kp*Part[4]*kA20/(kA20+A20*s6);
  
  double Tot = R[0] + R[1] +R[2] + R[3] +R[4] + R[5] +R[6] + R[7] +R[8] + R[9];
  if (randT == 1) {         randT = 0.999999999;          }
  if (randT == 0) {         randT = 0.000000001;          }
  double t1 = log(randT)/Tot;
  return t1;
}

int * GilUpdReact(int Part[], int U[], double s6, double A) {
  static int P[2];
  Ncy = Ntot-Part[0];   IKKn = IKKtot-Part[3]-Part[4];
  R[0] = kNin*Ncy*KI/(KI+Part[2]);      R[1] = kIin*Part[2]*Part[0]/(KN+Part[0]);
  R[2] = kt*pow(Part[0],2);     R[3] = ym*Part[1];
  R[4] = ktl*Part[1];        R[5] = a*Part[3]*Ncy*Part[2]/(KI+Part[2]);
  R[6] = ka*s6*IKKn;     R[7] = ki*Part[3];
  R[8] = ki*Part[3];         R[9] = kp*Part[4]*kA20/(kA20+A20*s6);

  double Tot = R[0] + R[1] +R[2] + R[3] +R[4] + R[5] +R[6] + R[7] +R[8] + R[9];
  if (A == 0) {     A = 0.00000001;       }
  if (A == 1) {     A = 0.9999999;        }
  
  double r1 = 0;      int ci = 0;  int x = 0;
  while (A > r1){
    r1 += R[ci]/Tot;
    if (A > r1){
      ci++;
    }
  }
  x = floor((double)(ci)/2.0 + 0.001);
  if (ci < 10){
    Part[x] += Upd[ci];
  }
  else {            std::cout << "Wrong" << std::endl;    }
  
  P[0] = x; P[1] = Part[x]; 
  return P;
}

double SpikePulse(double trt, double trt1, double T, double Ampc){
  double TNFa = T + (trt - trt1)*log(0.5)/85.0*T;
  if (TNFa < 0) {            TNFa = 0;  }
  return TNFa;
}

int Nover(int Nprev[], int L){
  double c3 = 0;
  for (unsigned int i = 0; i < L-1; i++){
    if (Nprev[i] > Nprev[i+1]){
      c3++;
    }
  }
  return c3;
}

int UpDo(int P1, double mu1, int no, int s){
  int s2 = 0;
  if (P1 > mu1 && no > 20){
    s2 = 1;
  }
  else if (P1 < mu1 && no < 9){
    s2 = -1;
  }
  else{
    s2 = s;
  }
  return s2;
}

