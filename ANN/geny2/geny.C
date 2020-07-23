#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "TBHDVCS.h"
#include "TBHDVCS.cxx"
#include "TFormFactors.cxx"
#include "TFormFactors.h"

using namespace std;

TBHDVCS *bhdvcs = new TBHDVCS;

TFormFactors *ff = new TFormFactors;

//_________________________________________________________________________
Double_t TotalUUXS(Double_t *angle, Double_t *par)
{
	Double_t phi = angle[0];

	// Set QQ, xB, t and k and calculate 4-vector products
	bhdvcs ->SetKinematics( par[0], par[1], par[2], par[3] );
	bhdvcs ->Set4VectorsPhiDep(phi);
	bhdvcs ->Set4VectorProducts(phi);

	Double_t AUUI, BUUI, CUUI;

	Double_t xsbhuu	 = bhdvcs ->GetBHUUxs(phi, par[4], par[5]);
	Double_t xsiuu	 = bhdvcs ->GetIUUxs(phi, par[4], par[5], par[6], par[7], par[8], AUUI, BUUI, CUUI);

	Double_t tot_sigma_uu = xsbhuu + xsiuu + par[9]; // Constant added to account for DVCS contribution

	return tot_sigma_uu;
}
//_________________________________________________________________________
void geny(){

 Double_t k; Double_t PhiX;
 Int_t seed = 400;
 //Use Hall A physics as inputs for the generated kinematics
 const Int_t NumOfSets = 20;
 Double_t QQ[NumOfSets] = {1.82, 1.933, 1.964, 1.986, 1.999, 2.218, 2.318, 2.348, 2.36, 2.375, 2.012, 2.054, 2.074, 2.084, 2.091, 2.161, 2.19, 2.194, 2.191, 2.193};
 Double_t xB[NumOfSets] = {0.343, 0.368, 0.375, 0.379, 0.381, 0.345, 0.363, 0.368, 0.371, 0.373, 0.378, 0.392, 0.398, 0.4, 0.401, 0.336, 0.342, 0.343, 0.342, 0.342};
 Double_t t[NumOfSets] = {-0.172, -0.232, -0.278, -0.323, -0.371, -0.176, -0.232, -0.279, -0.325, -0.372, -0.192, -0.233, -0.279, -0.324, -0.371, -0.171, -0.231, -0.278, -0.324, -0.371};

 TF1 *fl = new TF1("fl",TotalUUXS,0,360,10);

   ofstream myfile;
   myfile.open ("DVCS_cross.csv");

   for (Int_t m=0; m<6; m++) {
   k=2.75+m;
   for (Int_t n=0; n<NumOfSets; n++) {
    //generate the data from random seed
   TRandom *r   = new TRandom3();
   TRandom *s   = new TRandom3();
   r->SetSeed(seed+n);
   s->SetSeed(seed-n);
   //kinematic variation based on Hall A physics
   Double_t d_QQ =  s->Gaus(QQ[n],0.1*QQ[n]); //using 10% variation
   Double_t d_xB =  s->Gaus(xB[n],0.1*xB[n]); //using 10% variation
   Double_t d_t  =  s->Gaus(t[n],0.1*t[n]); //using 10% variation
	 cout<< "----- Kinematics -----"<<endl;
	 cout<< "k: "<<k<<" QQ: "<<d_QQ<<" xB: "<<d_xB<<" t: "<<d_t<<endl;

   fl->FixParameter(0, d_QQ); //QQ
   fl->FixParameter(1, d_xB); //xB
   fl->FixParameter(2, d_t); //t
   fl->FixParameter(3, k); //k

   fl->FixParameter(4, ff->ffF1(t[n])); //F1
   fl->FixParameter(5, ff->ffF2(t[n])); //F2
   //Using our standard model as a starting point
   fl->FixParameter(6, ff->ffF1(t[n])); //ReH
   fl->FixParameter(7, ff->ffF2(t[n])); //ReE
   fl->FixParameter(8, ff->ffGA(t[n])); //ReHtilde
   fl->FixParameter(9, 0.014863); //DVCSXS

   PhiX=0; //Generate cross section from kinematics
   for (Int_t i=0; i<36; i++) {
      PhiX = i*10+r->Gaus(5,3);// choose the phi
      Double_t F = r->Gaus(fl->Eval(PhiX),0.1*fl->Eval(PhiX));// find cross section with 10% variation in output
      Double_t f1 = r->Uniform(F,0.4);// Generate Uniform variation
      Double_t f2 = 0.1*F + r->Gaus(F*.01,0.001); // A simulated physical error (with absolute and relative contributions)
      //if(F>0){myfile <<k<<" "<<d_QQ<<","<<d_xB<<","<<d_t<<","<<PhiX<<","<<F<<","<<TMath::Abs(f2)<<endl;}
      if(F>0)myfile<<k<<" "<<d_QQ<<","<<d_xB<<","<<d_t<<","<<ff->ffF1(t[n])<<","<<ff->ffF2(t[n])<<","<<ff->ffGA(t[n])<<endl;
			if(F>0) cout<<"phi: "<<PhiX<<" F: "<<F<<"  f1: "<<f1<<" f2: "<<f2<<" ReH: "<<fl->GetParameter(6)<<" ReE: "<< fl->GetParameter(7)<<"ReHtilde: "<<fl->GetParameter(8)<<endl;
   }

	 Double_t flQQ = fl->GetParameter(0);
	 cout << "Func QQ: "<<flQQ<<endl;
	 fl->ReleaseParameter(0);

 }
}
}
