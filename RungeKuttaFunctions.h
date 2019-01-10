#include "ConstantHeaderDeterministic.h"
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

double f1 (double v1, double v3){
  double r1;
  r1 = kNin*(Ntot-v1)*KI/(KI+v3)-kIin*v3*v1/(KN+v1);
  return r1;
}
double f2 (double v1, double v2){
  double r2;
  r2 = kt*pow(v1,2)-ym*v2;
  return r2;
}
double f3 (double v1, double v2, double v3, double v4){
  double r3;
  r3 = ktl*v2-a*v4*(Ntot-v1)*v3/(KI+v3);
  return r3;
}
double f4 (double v4, double v5, double v6){
  double r4;
  r4 = ka*v6*(IKKtot-v5-v4)-ki*v4;
  return r4;
}
double f5 (double v4, double v5, double v6){
  double r5;
  r5 = ki*v4 - kp*v5*kA20/(kA20+A20*v6);
  return r5;
}

