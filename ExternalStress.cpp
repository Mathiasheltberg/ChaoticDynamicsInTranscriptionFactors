#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <ctype.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <random>
using namespace std;

int main(){
  
  int L = 20000; int M = 500; int N = 2;  double dt = 0.01; double RT;
  std::random_device rd2;    std::mt19937 genuni(rd2());    std::uniform_real_distribution<> dUni2(0.0 , 1.0);  


  double gamma1 [N];  double delta1 [N];  
  double h [N];  double K [N];  
  double Gamma2 [N];  double Delta2 [N];
  double P1 = 0.0;  double P2 = 0.0; double P0 = 0.000001;

  gamma1[0] = 1.0;  gamma1[1] = 5.0;
  delta1[0] = 0.03;  delta1[1] = 0.03;
  Gamma2[0] = 1.0;  Gamma2[1] = 5.0; 
  Delta2[0] = 0.001;  Delta2[1] = 0.001; 
  //////////// Specific genes ////////////
  h[0] = 2.0;  h[1] = 4.0; 
  K[0] = 1.0/6.023;  K[1] = 4.5/6.023;

  ostringstream filenameL;  filenameL << "SecondData/LiveFileChaos.txt"; 
  std::ofstream TriFileL (filenameL.str().c_str());
  ostringstream filename1;  filename1 << "SecondData/SaveDynamicsDeterChaos1.txt";  
  std::ofstream TriFile1 (filename1.str().c_str());
  ostringstream filename2;  filename2 << "SecondData/SaveDynamicsDeterChaos2.txt";  
  std::ofstream TriFile2 (filename2.str().c_str());

  for (int i = 0; i < M; i++){
    ifstream readFile ("Datafiles/DeterTraj_Chaos" + to_string(i) + ".txt");
    ////////// Read in NFKB
    double NFKB [L];        double a1,a2,a3;    int c = 0;
    while (readFile >> a1 >> a2 >> a3){
      NFKB[c] = a1/(2000*6.023);
      c++;
      if (c == L){	break;      }      
    }

    ////////// Initialize
    double NFRec;  double mRNA [N];     double Pro [N]; 
    int click = 0;
    Pro[0] = 0.75*5000;      mRNA[0] = 0.1*500;    Pro[1] = 0.75*1000;      mRNA[1] = 0.01*500;
    int Live = 1; int dea1 = 0; int dea2 = 0;    double D1 = 3000.0;  double D2 = 3000.0;

    RT = 0;
    while (Live == 1 && RT < L){
      RT = RT + dt;
      for (unsigned int u2 = 0; u2 < N; u2++){
	mRNA[u2] = mRNA[u2] + dt*(gamma1[u2]*pow(NFRec,h[u2])/(pow(NFRec,h[u2]) + pow(K[u2],h[u2])) - delta1[u2]*mRNA[u2]);
	Pro[u2] = Pro[u2] + dt*(Gamma2[u2]*mRNA[u2] - Delta2[u2]*Pro[u2]);	  
      }
      
      P1 = P0*pow(D1,4)/(pow(D1,4) + pow(Pro[0],4));
      P2 = P0*pow(D2,4)/(pow(D2,4) + pow(Pro[1],4));
      if (dUni2(genuni) < P1){
	Live = 0;
	TriFileL << RT << "\t" << 1 << "\n";
      }
      if (dUni2(genuni) < P2){
	Live = 0;
	TriFileL << RT << "\t" << 2 << "\n";       
      }
   
      if (RT > click){
	click++;
	NFRec = NFKB[click];
	TriFile1 << Pro[0] << "\t";	
	TriFile2 << Pro[1] << "\t";		
      }
    } // End of time loopq
    cout << Live << endl;
  TriFile1 << "\n"; TriFile2 << "\n";
  if (Live == 1){
    TriFileL << RT << "\t" << 0 << "\n";
  }
  }
}
