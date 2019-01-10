#include "ConstantHeaderGillespie.h"
#include "GillespieFunctions.h"
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

int main() {
  srand (time(0));
  ///////////////////// Set FileName ///////////////////// 
  Upd[0] = 1; Upd[2] = 1; Upd[4] = 1; Upd[6] = 1; Upd[8] = 1; Upd[1] = -1; Upd[3] = -1; Upd[5] = -1; Upd[7] = -1; Upd[9] = -1;
  
  double V1,V2,V3,V4,V5;  double P0a,P0b,P0c,P0d,P0e,P0f,P1a,P1b,P1c,P1d,P1e,P1f;
  ///////////////////////  RunScript  ///////////////////// 


  for (int vi = 0; vi < 1; vi++){
	nu = 50.0;
	Ac = vi*200.0;
	
	double gamma1,gamma2,delta1,delta2,Gam1,Gam2,Del1,Del2,h1,h2,K1,K2,h3,K3;
	double P1 = 5.0;  double P2 = 5.0; double m1 = 10.0; double m2 = 10.0; double m3 = 10.0; double P3 = 5.0;
	double ran1,ran2,uran1,uran2;
	double dN,dIRna,dI,dIKKa,dIKKi,dm1,dm2,dP1,dP2,dm3,dP3;

	gamma1 = 1.0;  gamma2 = 1.0;
	delta1 = 0.03;  delta2 = 0.03;
	Gam1 = 1.0;  Gam2 = 1.0;
	Del1 = 0.001;  Del2 = 0.001;
	h1 = 2.0;  h2 = 4.0; h3 = 3.0;
	K1 = 1.0/6.023;  K2 = 4.5/6.023; K3 = 2.0/6.023;

	double rf,drf1,drf2,wf,dwf,noirf3,noirf1,noirf2,noiwf; double con = 1.0;

	double nois1 = 0.30;
	double nois2 = 5.0;
	rf = 0.1; wf = 0;

	ostringstream trajname;
	trajname << "DataFiles/Distributions_V2_"; trajname << vi; trajname << "_4.txt";
	std::ofstream trajfile (trajname.str().c_str());
	double nVol = 2.0;
  	Vol = nVol*pow(10,-14); Cal = NA*Vol*pow(10,-6);
	Ntot = 1*Cal; KN = 0.029*Cal; kt = 1.03/Cal; a = 1.05/Cal;KI = 0.035*Cal; IKKtot = 2*Cal;

	N= 1500; P0b= 4982; P0c= 35444; P0d= 2285; P0e = 19073; P0f = 16.192;
	N= 1500; P1b= 3058; P1c= 33733; P1d= 2387; P1e = 18820; P1f = 11.58;
	V1 = P0b - P1b;     V2 = P0c - P1c;     V3 = P0d - P1d;     V4 = P0e - P1e;     V5 = P0f - P1f;
	int Place = 30;
	N = 1500.0;     IRna = P1b + V1/30.0*Place;     I = P1c + V2/30.0*Place;
	IKKa = P1d + V3/30.0*Place; IKKi = P1e + V4/30.0*Place;      RT = P1f + V5/30.0*Place;

	nu = 1/nu;      Ac = Ac/10000.0;

	double IPoin = I;	int NPoin = N;	double TPoin = RT;

	int MaxN = 0;	double Tmax = 500000.0;

	//////////// Helping objects ///////////////////////////
	int c = 0; click = 0; int cmax = 400;	int *p;
	Pp[0] = N; Pp[1] = IRna; 	Pp[2] = I; Pp[3] = IKKa; 	Pp[4] = IKKi;
	int swmu = 0; int I0 = 0;	int N0 = 0;
	/////////////////////////////////////////////////////////////////////////////////////////////////
	//                                      Start actual run                                       //
	/////////////////////////////////////////////////////////////////////////////////////////////////
	cout << vi << endl;
	while (RT < Tmax) {

	  ////////////////////////////// Updates //////////////////////////////////////////////////////
	  dt = - GilUpdTime(Pp,TNF,double(rand()/double(RAND_MAX)));
	  RT = RT + dt; // TimeStep                  //
	  p = GilUpdReact(Pp , Upd , TNF,double(rand()/double(RAND_MAX))); // Reaction               //
	  Pp[*(p+0)] = *(p+1); // Update                                                             //
	  /////////////////////////////////////////////////////////////////////////////////////////////
	  
	  TNF = 0.5 + Ac*sin(2*3.141592*nu*RT);

	  dm1   = gamma1*pow(Pp[0]/(nVol*1000*6.023),h1)/(pow(Pp[0]/(nVol*1000*6.023),h1) + pow(K1,h1)) - delta1*m1;
	  dP1   = Gam1*m1 - Del1*P1;
	  dm2   = gamma2*pow(Pp[0]/(nVol*1000*6.023),h2)/(pow(Pp[0]/(nVol*1000*6.023),h2) + pow(K2,h2)) - delta2*m2;
	  dP2   = Gam2*m2 - Del2*P2;
	  dm3   = gamma2*pow(Pp[0]/(nVol*1000*6.023),h3)/(pow(Pp[0]/(nVol*1000*6.023),h3) + pow(K3,h3)) - delta2*m3;
	  dP3   = Gam2*m3 - Del2*P3;


	  m1   = m1   + dt*(dm1);
	  m2   = m2   + dt*(dm2);
	  m3   = m3   + dt*(dm3);
	  P1   = P1   + dt*(dP1);
	  P2   = P2   + dt*(dP2);
          P3   = P3   + dt*(dP3);
	  /////////////////// PeakFinding ///////////////////////////////

	  if (RT > DT*click){
	    click++;
	    trajfile << RT << "\t" << TNF << "\t" << Pp[0] << "\t" << P1 << "\t" << P2 << "\t" << P3 << "\n";
	  }
	  I0 = Pp[2];
	  N0 = Pp[0];

	}

  }
}

