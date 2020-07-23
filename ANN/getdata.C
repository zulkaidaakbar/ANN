#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include "TROOT.h"
#include "TStyle.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TH1I.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
//#include "defineHistos.h"
#include "TH1F.h"



using namespace std; 

//
Double_t TProduct(TLorentzVector v1, TLorentzVector v2);

void SetKinematics(Double_t _QQ, Double_t _x, Double_t _t, Double_t _k);

void Set4VectorsPhiDep(Double_t phi);

void Set4VectorProducts(Double_t phi);

Double_t GetBHUUxs(Double_t F1, Double_t F2);

Double_t GetIUUxs(Double_t phi, Double_t F1, Double_t F2, Double_t ReH, Double_t ReE, Double_t ReHtilde);



void getdata()
{

//variable
double temp[38844];
float k[83];
float QQ[83];
float xb[83];
float t[83];
double phi;
double phi_rad;
double F;
double errF;

//
Double_t ALP_INV = 137.0359998; // 1 / Electromagnetic Fine Structure Constant
  Double_t PI = 3.1415926535;
  Double_t RAD = PI / 180.;
	Double_t M = 0.938272; //Mass of the proton in GeV
  Double_t GeV2nb = .389379*1000000; // Conversion from GeV to NanoBarn
  // Elastic FF
  // Double_t F1;  // Dirac FF - helicity conserving (non spin flip)
  // Double_t F2;  // Pauli FF - helicity non-conserving (spin flip)

	Double_t QQ, x, t, k;// creo q no hace falta
	Double_t y, e, xi, tmin, kpr, gg, q, qp, po, pmag, cth, theta, sth, sthl, cthl, cthpr, sthpr, M2, tau;

  // 4-momentum vectors
  TLorentzVector K, KP, Q, QP, D, p, P;
  // 4 - vector products independent of phi
  Double_t kkp, kq, kp, kpp;
  // 4 - vector products dependent of phi
  Double_t kd, kpd, kP, kpP, kqp, kpqp, dd, Pq, Pqp, qd, qpd;

	//Double_t KK_T, KQP_T, KKP_T, KXQP_T, KD_T, DD_T;
  Double_t kk_T, kqp_T, kkp_T, kd_T, dd_T;

  Double_t s;     // Mandelstam variable s which is the center of mass energy
  Double_t Gamma; // Factor in front of the cross section
  Double_t jcob;  //Defurne's Jacobian

  Double_t AUUBH, BUUBH; // Coefficients of the BH unpolarized structure function FUU_BH
  Double_t AUUI, BUUI, CUUI; // Coefficients of the BHDVCS interference unpolarized structure function FUU_I
  Double_t con_AUUBH, con_BUUBH, con_AUUI, con_BUUI, con_CUUI;  // Coefficients times the conversion to nb and the jacobian
  Double_t bhAUU, bhBUU; // Auu and Buu term of the BH cross section
  Double_t iAUU, iBUU, iCUU; // Terms of the interference containing AUUI, BUUI and CUUI
  Double_t xbhUU; // Unpolarized BH cross section
  Double_t xIUU; // Unpolarized interference cross section

  
//

char name[100];
char name1[100];
char name2[100];

double pi = 3.141592;
int nbin = 126;

TFile f ("testdata.root", "RECREATE", "Histograms from ntuples" ); 

//get the data from text file  
ifstream inputfile;
sprintf(name,"data1.txt"); //copy from dvcc_cross_fixed_t.csv
inputfile.open(name,ios::in);
cout<<" Get the data from "<<name<<endl<<endl;
  
for (int i = 0; i < 38844 ; i++) 
 {
  inputfile >> temp[i];
 } 

//define the histograms - filling the histogram - write the histograms
gStyle->SetOptFit(1);
TH1D *hist_F[83];
TF1 *fcos = new TF1("fcos","[0]*([1]*cos(x)+1)");
for (int i=0; i<83; i++)
{
 k[i]= temp[i*468 + 1];
 QQ[i]= temp[i*468 + 2];
 xb[i]= temp[i*468 + 3];
 t[i]= temp[i*468 + 4];
 sprintf(name1,"kinematics: k=%f, QQ=%f, xb=%f, t=%f",k[i],QQ[i],xb[i],t[i]);
 cerr<<name1<<endl; 
   
 hist_F[i] = new TH1D(name1,name1, nbin,0,2*pi);
 hist_F[i]->Sumw2();

 for (int j=0; j<36; j++)
  {
   phi = temp[i*468+j*13+5];
   F = temp[i*468+j*13+6];
   errF = temp[i*468+j*13+7];

   phi_rad = (phi/180)*pi;
   double bin_space = /*2*phi/nbin*/0.049867;
   int bin_phi = (phi_rad-0)/bin_space + 1;
   hist_F[i]->SetBinContent(bin_phi, F);
   hist_F[i]->SetBinError(bin_phi, errF);   
  }
 hist_F[i]->Fit(fcos);
 hist_F[i]->Write();
}

//save the canvas for each histogram
TCanvas *cn[83];
for (int i=0; i<20; i++)
{
 sprintf(name1,"c%i",i);
 sprintf(name2,"c%i_.png",i);
 cn[i] = new TCanvas(name1,name1);
 cn[i]->cd();
 //hist_F[i]->SetStats(0);
 hist_F[i]->SetOption("p");
hist_F[i]->GetXaxis()->SetTitle("#phi_{x}");
hist_F[i]->GetYaxis()->SetTitle("F");
 hist_F[i]->Draw("e");
 cn[i]->SaveAs(name2);
}



//plotting the histograms
/*
float xmin = 0.0;
float xmax = 360.0;
float ymin = 0.00; 
float ymax = 13;
int divx = 5;
int divy = 4;
bool ylog =0;

char xaxis[100]=" #phi_{x}";
char yaxis[100]="F ";

TLatex xTitle,yTitle;
TGaxis *yax,*xax;
Float_t small = 1e-5;
  
gROOT->SetStyle("Plain");

gStyle->SetOptTitle(0);
gStyle->SetOptStat(0000000);
gStyle->SetPadBorderMode(0);
gStyle->SetFrameBorderMode(0);
gStyle->SetLabelSize(.1,"Y");
gStyle->SetLabelSize(.1,"X");
//gStyle->SetErrorY(0);
//gROOT->SetBatch();

TCanvas *tC;
TPad *c;
char NPad[100];
char NCanvas[100];
sprintf(NCanvas,"F_dvcs");
sprintf(NPad,"c_cross");
tC = new TCanvas(NCanvas,"DVCS",30,10,1000,800);

c = new TPad(NPad,"pad", .1, .1, .95, .95);
c->Draw();
c->Divide(divx,divy,small,small);

for (int e=1; e<20; e++) 
{ 
   c->cd(e);
    
    gPad->SetFillStyle(0);
    gPad->SetFrameFillStyle(0);
    gPad->SetTopMargin(small);
    gPad->SetBottomMargin(small);
    gPad->SetRightMargin(small);
    gPad->SetLeftMargin(small);
    if(ylog==1){
      gPad->SetLogy();
      }
    
    hist_F[e]->SetMaximum(ymax);
    hist_F[e]->SetMinimum(ymin);
    hist_F[e]->SetAxisRange(xmin,xmax);
    hist_F[e]->SetLineColor(kRed);
    hist_F[e]->SetMarkerColor(2);
    hist_F[e]->SetMarkerStyle(20);
    hist_F[e]->SetMarkerSize(0.8);
    hist_F[e]->Draw("P ");
    


    /////////////////////////////////
    // Titles for individual plots //
    /////////////////////////////////
    TLatex plotTitle;
    char title[200], title2[200];
    float title_k= temp[e*468 + 1];
    float title_QQ= temp[e*468 + 2];
    float title_xb= temp[e*468 + 3];
    float title_t= temp[e*468 + 4];
    
    sprintf(title,"k=%f, QQ=%f, xb=%f, t=%f",title_k,title_QQ,title_xb,title_t);
    plotTitle.SetTextSize(0.092);
    plotTitle.DrawLatex(50,12,title);

  
  }

//////////
  // Axes //
  //////////
  tC->cd();

  for(int i=0;i<divy;i++) {
    if(ylog==1){
      yax = new TGaxis(.1,(((.95-.1)/divy)*i)+.1,.1,(((.95-.1)/divy)*(i+1))+.1,ymin,ymax,3,"+G");
      }
    else{
      yax = new TGaxis(.1,(((.95-.1)/divy)*i)+.1,.1,(((.95-.1)/divy)*(i+1))+.1,ymin,ymax,3,"+");
      }
    yax->SetLabelSize(.025);
    yax->SetLabelOffset(-.015);
    yax->Draw();
    }

  for (int i=0;i<divx;i++) {

    xax = new TGaxis((((.95-.1)/divx)*i)+.1,.1,(((.95-.1)/divx)*(i+1))+.1,.1,xmin,xmax,5,"-");
    xax->SetLabelSize(.025);
    xax->SetLabelOffset(-.03);
    xax->Draw();
  }

  ///////////////////
  // Axes labeling //
  /////////////////// 
  xTitle.SetTextSize(0.037);
  xTitle.DrawLatex(0.5,0.02,xaxis); 

  yTitle.SetTextSize(0.035);
  yTitle.SetTextAngle(90);
  yTitle.DrawLatex(0.04,0.45,yaxis); 

  
tC->Print("plot_test.png");


*/





}    
