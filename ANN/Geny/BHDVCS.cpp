//Calculation of the Bethe-Heitler Deeply Virtual Compton Scattering Cross Section 
//Interference with all Beam and Target Polarization Configurations
//Referenced paper arXiv : 1903.05742
//Calculated by Brandon Kriesten: btk8bh@virginia.edu
//Work is done as part of the paper cited above
//Used for extraction of Observables from Deeply Virtual Compton Scattering Experiments

#include <cmath>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include "four_vector.h"
#include "form_factor.h"

using namespace std;

int main() {

  //Opens Unpolarized Cross Section Data file
  ofstream myfile;
  myfile.open("XSXBHDVCSUU.dat");

  //Opens Polarized Cross Section Data File
  ofstream myfile2;
  myfile2.open("XSXBHDVCSLL.dat");
 
  ofstream myfile3;
  myfile3.open("XSXBHDVCSUL.dat");

  ofstream myfile4;
  myfile4.open("XSXBHDVCSLU.dat");

  //Sets kinematic variables used
  double Q2, k0, xbj, t, M;
  double nu, xi, gamma, y;
  double k0p;
  double costl, sintl, sintlp, costlp;
  double p0p;
  double t0;
  double eps;
  double cost, sint, costp, sintp;
  double q0p;
  double tau;
  double alpha;
  double s;
  long double Gamma;
  double conversion;

  //Sets the value of pi
  const double pi = 3.14159265358979323846;

  //Initialized Kinematic Variables
  //Q^2 value
  Q2 = 1.82;

  //Electron Beam Energy
  k0 = 5.75;

  //Bjorken x
  xbj = 0.34;

  //t value as the momentum transfer squared
  t = -0.172;

  //Mass of the proton in GeV
  M = .938;

  //Electromagnetic Fine Structure Constant
  alpha = 0.00729927007;
  
  //Conversion from GeV to NanoBarn
  conversion = .389379*1000000;

  //Secondary Kinematic Variables
  //Energy of the virtual photon
  nu = Q2/(2*M*xbj);

  //Skewness parameter set by xbj, t, and Q^2
  xi =xbj *((1+(t/(2*Q2)))/(2-xbj+((xbj*t)/Q2)));

  //gamma variable ratio of virtuality to energy of virtual photon
  gamma = sqrt(Q2)/nu;

  //fractional energy of virtual photon
  y = sqrt(Q2)/(gamma*k0);

  //Final Lepton energy
  k0p = k0*(1-y);

  //Minimum t value
  t0 = -(4*xi*xi*M*M)/(1-(xi*xi));
  
  //Lepton Angle Kinematics of initial lepton
  costl = -(1/(sqrt(1+gamma*gamma)))*(1 + (y*gamma*gamma/2));
  sintl = (gamma/(sqrt(1+gamma*gamma)))*sqrt(1-y-(y*y*gamma*gamma/4));

  //Lepton Angle Kinematics of final lepton
  sintlp = sintl/(1-y);
  costlp = (costl + y*sqrt(1+gamma*gamma))/(1-y);
  
  //Final proton Energy
  p0p = M - (t/(2*M));

  //ratio of longitudinal to transverse virtual photon flux
  eps = (1-y-0.25*y*y*gamma*gamma)/(1-y+.5*y*y+.25*y*y*gamma*gamma);

  //Angular Kinematics of outgoing photon
  cost = - (1/(sqrt(1+gamma*gamma)))*(1 + (.5*gamma*gamma)*((1+(t/Q2))/(1+((xbj*t)/(Q2)))));
  sint = sqrt(1-cost*cost);
  
  //outgoing photon energy
  q0p = (sqrt(Q2)/gamma)*(1+((xbj*t)/Q2));
  
  //ratio of momentum transfer to proton mass
  tau = -t/(4*M*M);

  //Creates arrays of 4-vector kinematics used in Bethe Heitler Cross Section
  double p[] = {M,0,0,0};
  double k[] = {k0,k0*sintl,0,k0*costl};
  double kp[] = {k0p,k0p*sintlp,0,k0p*costlp};
  double q[] = {k[0]-kp[0],k[1]-kp[1],k[2]-kp[2],k[3]-kp[3]};

  //Initializes doubles
  double plp, qq, kk ,kkp, kq, pk, pkp;

  //Uses product function from four_vector.h which takes a 4-vector product
  //using the Minkowski Metric
  plp = product(p,p);  //(pp)
  qq = product(q,q);  //(qq)
  kk = product(k,k);  //(kk)
  kkp = product(k,kp);  //(kk')
  kq = product(k,q);  //(kq)
  pk = product(k,p);  //(pk)
  pkp = product(kp,p);  //(pk')
  
  //Sets the Mandelstam variable s which is the center of mass energy
  s = product(k,k) + 2*product(k,p) + product(p,p);
  
  //The Gamma factor in front of the cross section
  Gamma = (alpha*alpha*alpha)/(16*pi*pi*(s-M*M)*(s-M*M)*sqrt(1+gamma*gamma)*xbj);

  //Looping over phi
  for (int i=0; i<=360; i++){
    
    //Sets the phi in radians from degrees
    double phi = i *.0174532951;

    //Initializes the 4-vectors that depend on phi
    double qp[] = {q0p,q0p*sint*cos(phi),q0p*sint*sin(phi),q0p*cost};
    double d[] = {q[0]-qp[0], q[1]-qp[1],q[2]-qp[2],q[3]-qp[3]};
    double pp[] = {p[0]+d[0],p[1]+d[1],p[2]+d[2],p[3]+d[3]};
    double P[] = {.5*(p[0]+pp[0]),.5*(pp[1]),.5*(pp[2]),.5*(pp[3])};

    //Initializes doubles used that depend on phi
    double kd,kpd,kP,kpP,kqp,kpqp,dd;
    double Pq, Pqp, qd, qpd;
    double kkT, kqpT, kkpT, ddT, kdT;
    double kpqpT, qpdT, kPT, kpPT, qpPT, kpdT;
    double kplus, kpplus, kminus,kpminus, qplus,qpplus;
    double qminus,qpminus,Pplus,Pminus,dplus,dminus;
    double Dplus,Dminus;
    double AUUBHDVCS, BUUBHDVCS, CUUBHDVCS;
    double ALLBHDVCS, BLLBHDVCS, CLLBHDVCS;
    double AULBHDVCS, BULBHDVCS, CULBHDVCS;
    double ALUBHDVCS, BLUBHDVCS, CLUBHDVCS;

    //4-vector products denoted in the paper by the parentheses
    kd = product(k,d);  //(kΔ)
    kpd = product(kp,d);  //(k'Δ)
    kP = product(k,P);  //(kP)
    kpP = product(kp,P);  //(k'P)
    kqp = product(k,qp);  //(kq')
    kpqp = product(kp,qp);  //(k'q')
    dd = product(d,d);  //(ΔΔ)
    Pq = product(P,q);  //(Pq)
    Pqp = product(P,qp);  //(Pq')
    qd = product(q,d);  //(qΔ)
    qpd = product(qp,d);  //(q'Δ)

    //Transverse vector products
    kkT = tproduct(k,k);  //(kk)T
    kqpT = tproduct(k,qp);
    kkpT = tproduct(k,kp);
    ddT = tproduct(d,d);
    kdT = tproduct(k,d);
    kpqpT = tproduct(kp,qp);
    qpdT = tproduct(qp,d);
    kPT = tproduct(k,P);
    kpPT = tproduct(kp,P);
    qpPT = tproduct(qp,P);
    kpdT = tproduct(kp,d);

    //Light cone variables expressed as A^{+-} = 1/sqrt(2)(A^{0} +- A^{3})
    kplus =( 1/sqrt(2))*(k[0]+k[3]);
    kpplus = ( 1/sqrt(2))*(kp[0]+kp[3]);
    kminus = ( 1/sqrt(2))*(k[0]-k[3]);
    kpminus = ( 1/sqrt(2))*(kp[0]-kp[3]);
    qplus = ( 1/sqrt(2))*(q[0]+q[3]);
    qpplus = ( 1/sqrt(2))*(qp[0]+qp[3]);
    qminus = ( 1/sqrt(2))*(q[0]-q[3]);
    qpminus = ( 1/sqrt(2))*(qp[0]-qp[3]);
    Pplus = ( 1/sqrt(2))*(P[0]+P[3]);
    Pminus = ( 1/sqrt(2))*(P[0]-P[3]);
    dplus = ( 1/sqrt(2))*(d[0]+d[3]);
    dminus = ( 1/sqrt(2))*(d[0]-d[3]);

    //Expressions used that appear in coefficient calculations
    Dplus = (1/(2*kpqp))-(1/(2*kqp));
    Dminus = (1/(2*kpqp))+(1/(2*kqp));

    //Calculates the Unpolarized Coefficients in front of the Elastic Form Factors and
    //Compton Form Factors
    //No conversion factor to nano barn, nor the Gamma/t factor is included
    AUUBHDVCS = -16*Dplus*((kqpT - 2*kkT-2*kqp)*kpP+(2*kpqp-2*kkpT-kpqpT)*kP+kpqp*kPT+kqp*kpPT-2*kkp*kPT)*cos(phi) - 16*Dminus*((2*kkp-kpqpT-kkpT)*Pqp + 2*kkp*qpPT - kpqp*kPT-kqp*kpPT)*cos(phi);
    BUUBHDVCS = -8*xi*Dplus*((kqpT - 2*kkT - 2*kqp)*kpd+(2*kpqp-2*kkpT-kpqpT)*kd+kpqp*kdT+kqp*kpdT-2*kkp*kdT)*cos(phi)-8*xi*Dminus*((2*kkp-kpqpT-kkpT)*qpd+2*kkp*qpdT-kpqp*kdT-kqp*kpdT)*cos(phi);
    CUUBHDVCS = -8*Dplus*((2*kkp*kdT - kpqp*kdT-kqp*kpdT)+2*xi*(2*kkp*kPT - kpqp*kPT-kqp*kpPT))*cos(phi)-8*Dminus*((kkp*qpdT-kpqp*kdT-kqp*kpdT)+2*xi*(kkp*qpPT-kpqp*kPT-kqp*kpPT))*cos(phi);

    //Calculates the Unpolarized Beam Polarized Target Coefficients in front of the 
    //Elastic Form Factors and Compton Form Factors                                         
    //No conversion factor to nano barn, nor the Gamma/t factor is included                  
    AULBHDVCS = 16*Dplus*(kpP*(2*kkT-kqpT+2*kqp)+kP*(2*kkpT-kpqpT+2*kpqp)+2*kkp*kPT-kpqp*kPT-kqp*kpPT)*sin(phi)+16*Dminus*(Pqp*(kkpT+kpqpT-2*kkp)-(2*kkp*qpPT-kpqp*kPT-kqp*kpPT))*sin(phi);
    BULBHDVCS = 8*xi*Dplus*(kpd*(2*kkT - kqpT + 2*kqp)+kd*(2*kkpT-kpqpT+2*kpqp)+2*kkp*kpdT-kpqp*kdT-kqp*kpdT)*sin(phi)+8*xi*Dminus*(qpd*(kkpT+kpqpT-2*kkp)-(2*kkp*qpdT-kpqp*kdT-kqp*kpdT))*sin(phi);
    CULBHDVCS = 4*Dplus*(2*(2*kkp*kdT-kpqp*kdT-kqp*kpdT)+4*xi*(2*kkp*kPT -kpqp*kPT-kqp*kpPT))*sin(phi)+4*Dminus*(-2*(kkp*qpdT-kpqp*kdT-kqp*kpdT)-4*xi*(kkp*qpPT-kpqp*kPT-kqp*kpPT))*sin(phi);

    //Calculates the Polarized Beam Unpolarized Target Coefficients in front of the          
    //Elastic Form Factors and Compton Form Factors                                          
    //No conversion factor to nano barn, nor the Gamma/t factor is included                  
    ALUBHDVCS = 16*Dplus*(2*(k[1]*Pplus*kp[1]*kminus - k[1]*Pplus*kpminus*k[1]+k[1]*Pminus*kpplus*k[1]-k[1]*Pminus*kp[1]*kplus+k[1]*P[1]*kpminus*kplus-k[1]*P[1]*kpplus*kminus)+kp[1]*Pplus*qpminus*k[1]-kp[1]*Pplus*qp[1]*kminus+kp[1]*Pminus*qp[1]*kplus-kp[1]*Pminus*qpplus*k[1]+kp[1]*P[1]*qpplus*kminus-kp[1]*P[1]*qpminus*kplus+k[1]*Pplus*qpminus*kp[1]-k[1]*Pplus*qp[1]*kpminus+k[1]*Pminus*qp[1]*kpplus-k[1]*Pminus*qpplus*kp[1]+k[1]*P[1]*qpplus*kpminus-k[1]*P[1]*qpminus*kpplus+2*(qpminus*Pplus-qpplus*Pminus)*kkp)*sin(phi)+16*Dminus*(2*(kminus*kpplus-kplus*kpminus)*Pqp+kpminus*kplus*qp[1]*P[1]+kpplus*k[1]*qpminus*P[1]+kp[1]*kminus*qpplus*P[1]-kpplus*kminus*qp[1]*P[1]-kp[1]*kplus*qpminus*P[1]-kpminus*k[1]*qpplus*P[1]+kpminus*kplus*qp[2]*P[2]-kpplus*kminus*qp[2]*P[2])*sin(phi);
    BLUBHDVCS = 8*xi*Dplus*(2*(k[1]*dplus*kp[1]*kminus-k[1]*dplus*kpminus*k[1]+k[1]*dminus*kpplus*k[1]-k[1]*dminus*kp[1]*kplus+k[1]*d[1]*kpminus*kplus-k[1]*d[1]*kpplus*kminus)+kp[1]*dplus*qpminus*k[1] - kp[1]*dplus*qp[1]*kminus+kp[1]*dminus*qp[1]*kplus-kp[1]*dminus*qpplus*k[1]+kp[1]*d[1]*qpplus*kminus-kp[1]*d[1]*qpminus*kplus+k[1]*dplus*qpminus*kp[1]-k[1]*dplus*qp[1]*kpminus+k[1]*dminus*qp[1]*kpplus-k[1]*dminus*qpplus*kp[1]+k[1]*d[1]*qpplus*kpminus-k[1]*d[1]*qpminus*kpplus+2*(qpminus*dplus-qpplus*dminus)*kkp)*sin(phi)+8*xi*Dminus*(2*(kminus*kpplus-kplus*kpminus)*qpd+kpminus*kplus*qp[1]*d[1]+kpplus*k[1]*qpminus*d[1]+kp[1]*kminus*qpplus*d[1]-kpplus*kminus*qp[1]*d[1]-kp[1]*kplus*qpminus*d[1]-kpminus*k[1]*qpplus*d[1]+kpminus*kplus*qp[2]*d[2]-kpplus*kminus*qp[2]*d[2])*sin(phi);
    CLUBHDVCS = -8*Dplus*(2*(kp[1]*kpminus*kplus*d[1]-kp[1]*kpplus*kminus*d[1])+kp[1]*qpminus*kplus*d[1]-kp[1]*qpplus*kminus*d[1]+k[1]*qpminus*kpplus*d[1]-k[1]*qpplus*kpminus*d[1])*sin(phi)-8*Dminus*(-kpminus*k[1]*qpplus*d[1]+kpminus*kplus*qp[1]*d[1]+kpplus*k[1]*qpminus*d[1]-kpplus*kminus*qp[1]*d[1]+kp[1]*kminus*qpplus*d[1]-kp[1]*kplus*qpminus*d[1]-qp[2]*d[2]*(kpplus*kminus-kpminus*kplus))*sin(phi);

    //Calculates the Longitudinally Polarized Coefficients in front of the EFFs
    //No conversion factor to nano barn, nor the Gamma/t factor is included
    ALLBHDVCS = 16*Dplus*(2*kp[1]*(kp[1]*kminus-kpminus*k[1])*Pplus+2*kp[1]*(kpplus*k[1]-kp[1]*kplus)*Pminus+2*kp[1]*(kpminus*kplus-kpplus*kminus)*P[1]+kp[1]*(qpminus*k[1]-qp[1]*kminus)*Pplus+kp[1]*(qp[1]*kplus-qpplus*k[1])*Pminus+kp[1]*(qpplus*kminus-qpminus*kplus)*P[1]+k[1]*(qpminus*kp[1]-qp[1]*kpminus)*Pplus+k[1]*(qp[1]*kpplus-qpplus*kp[1])*Pminus+k[1]*(qpplus*kpminus-qpminus*kpplus)*P[1]-2*kkp*(qpplus*Pminus-qpminus*Pplus))*cos(phi)+16*Dminus*(2*Pqp*(kpplus*kminus-kpminus*kplus)+P[2]*qp[2]*(kpplus*kminus-kpminus*kplus)+P[1]*(kp[1]*kplus*qpminus-kp[1]*kminus*qpplus+kpminus*k[1]*qpplus-kpminus*kplus*qp[1]+kpplus*kminus*qp[1]-kpplus*k[1]*qpminus))*cos(phi);
    BLLBHDVCS = 8*xi*Dplus*(2*kp[1]*(kp[1]*kminus-kpminus*k[1])*dplus+2*kp[1]*(kpplus*k[1]-kp[1]*kplus)*dminus+2*kp[1]*(kpminus*kplus-kpplus*kminus)*d[1]+kp[1]*(qpminus*k[1]-qp[1]*kminus)*dplus+kp[1]*(qp[1]*kplus-qpplus*k[1])*dminus+kp[1]*d[1]*(qpplus*kminus-qpminus*kplus)+k[1]*(qpminus*kp[1]-qp[1]*kpminus)*dplus+k[1]*(qp[1]*kpplus-qpplus*kp[1])*dminus+k[1]*(qpplus*kpminus-qpminus*kpplus)*d[1]-2*kkp*(qpplus*dminus-qpminus*dplus))*cos(phi)+8*xi*Dminus*(2*qpd*(kpplus*kminus-kpminus*kplus)+d[2]*qp[2]*(kpplus*kminus-kpminus*kplus)+d[1]*(kp[1]*kplus*qpminus-kp[1]*kminus*qpplus+kpminus*k[1]*qpplus-kpminus*kplus*qp[1]+kpplus*kminus*qp[1]-kpplus*k[1]*qpminus))*cos(phi);
    CLLBHDVCS = 16*Dplus*(2*(k[1]*kminus*kpplus*d[1]-k[1]*kpminus*kplus*d[1])+kp[1]*qpplus*kminus*d[1]-kp[1]*qpminus*kplus*d[1]+k[1]*qpplus*kpminus*d[1]-k[1]*qpminus*kpplus*d[1])*cos(phi)-16*Dminus*(-d[1]*(kpminus*kplus*qp[1]-kpminus*k[1]*qpplus+kpplus*k[1]*qpminus-kpplus*kminus*qp[1]+kp[1]*kminus*qpplus-kp[1]*kplus*qpminus)+qp[1]*kpplus*kminus*d[2]-qp[2]*kpminus*kplus*d[2])*cos(phi);

    //Converted Unpolarized Coefficients with the Gamma factor and in nano-barn
    double con_AUUBHDVCS = (Gamma/(Q2*-t))*AUUBHDVCS*conversion;
    double con_BUUBHDVCS = (Gamma/(Q2*-t))*BUUBHDVCS*conversion;
    double con_CUUBHDVCS = (Gamma/(Q2*-t))*CUUBHDVCS*conversion;

    //Converted Longitudinally Polarized Coefficients with the Gamma Factor and in nano-barn
    double con_ALLBHDVCS = (Gamma/(Q2*-t))*ALLBHDVCS*conversion;
    double con_BLLBHDVCS = (Gamma/(Q2*-t))*BLLBHDVCS*conversion;
    double con_CLLBHDVCS = (Gamma/(Q2*-t))*CLLBHDVCS*conversion;

    //Converted Longitudinally Polarized Beam Unpolarized Target Coefficients with 
    //the Gamma Factor and in nano-barn  
    double con_ALUBHDVCS = (Gamma/(Q2*-t))*ALUBHDVCS*conversion;
    double con_BLUBHDVCS = (Gamma/(Q2*-t))*BLUBHDVCS*conversion;
    double con_CLUBHDVCS = (Gamma/(Q2*-t))*CLUBHDVCS*conversion;

    //Converted Longitudinally Polarized Target Unpolarized Beam Coefficients with 
    //the Gamma Factor and in nano-barn 
    double con_AULBHDVCS = (Gamma/(Q2*-t))*AULBHDVCS*conversion;
    double con_BULBHDVCS = (Gamma/(Q2*-t))*BULBHDVCS*conversion;
    double con_CULBHDVCS = (Gamma/(Q2*-t))*CULBHDVCS*conversion;

    //Unpolarized Coefficients multiplied by the Form Factors calculated in form_factor.cpp
    //We use the Galster Form Factors as approximations to the Compton Form Factors
    double bhdvcsAUU = con_AUUBHDVCS*(ffF1(t)*ffF1(t) - tau*ffF2(t)*ffF2(t));
    double bhdvcsBUU = con_BUUBHDVCS*(ffGM(t)*ffGM(t));
    double bhdvcsCUU = con_CUUBHDVCS*(ffGM(t)*ffGA(t));

    //Polarized Coefficients multiplied by the Form Factors calculated in form_factor.cpp
    //Using the Galster Form Factor Model
    double bhdvcsALL = con_ALLBHDVCS*(ffF1(t)*ffF1(t) - tau*ffF2(t)*ffF2(t));
    double bhdvcsBLL = con_BLLBHDVCS*(ffGM(t)*ffGM(t));
    double bhdvcsCLL = con_CLLBHDVCS*(ffGM(t)*ffGA(t));

    //Unpolarized Beam Polarized Target Coefficients multiplied by the Form Factors 
    //calculated in form_factor.cpp  
    //We use the Galster Form Factors as approximations                                      
    double bhdvcsAUL = con_AULBHDVCS*(ffF1(t)*ffGA(t) - xi*ffF1(t)*ffGP(t)-tau*ffF2(t)*ffGP(t));
    double bhdvcsBUL = con_BULBHDVCS*(ffGM(t)*ffGA(t));
    double bhdvcsCUL = con_CULBHDVCS*(ffGM(t)*ffGM(t));

    //Polarized Beam Unpolarized Target Coefficients multiplied by the Form Factors 
    //calculated in form_factor.cpp    
    //Using the Galster Form Factor Model                                                    
    double bhdvcsALU = con_ALUBHDVCS*(ffF1(t)*ffGA(t) - xi*ffF1(t)*ffGP(t)-tau*ffF2(t)*ffGP(t));
    double bhdvcsBLU = con_BLUBHDVCS*(ffGM(t)*ffGA(t));
    double bhdvcsCLU = con_CLUBHDVCS*(ffGM(t)*ffGM(t));

    //Calculation of the Total Unpolarized Bethe-Heitler Cross Section
    double XSXUUBHDVCS = bhdvcsAUU + bhdvcsBUU + bhdvcsCUU;
    double XSXLLBHDVCS = bhdvcsALL + bhdvcsBLL + bhdvcsCLL;
    double XSXULBHDVCS = bhdvcsAUL + bhdvcsBUL + bhdvcsCUL;
    double XSXLUBHDVCS = bhdvcsALU + bhdvcsBLU + bhdvcsCLU;

    //Writing Unpolarized Data to the file XSXBHDVCSUU.dat
    //First Column: phi
    //Second Column: Total Cross Section
    //Third Column: A Coefficient with elastic form factors
    //Fourth Column: B Coefficient with elastic form factors
    //Fifth Column: C Coefficient with elastic form factors
    myfile << i << " " <<  XSXUUBHDVCS << " " << bhdvcsAUU << " " << bhdvcsBUU << " " << bhdvcsCUU<< endl;
    
    //Writing Polarized Data to the file XSXBHDVCSLL.dat
    //First Column: phi
    //Second Column: Total Polarized Cross Section
    //Third Column: A Coefficient with elastic form factors
    //Fourth Column: B Coefficient with elastic form factors
    //Fifth Column: C Coefficient with elastic form factors
    myfile2 << i << " " << XSXLLBHDVCS << " " << bhdvcsALL << " " << bhdvcsBLL << " " << bhdvcsCLL <<  endl;

    //Writing Unpolarized Data to the file XSXBHDVCSUL.dat                                  
    //First Column: phi                                                                      
    //Second Column: Total Cross Section                                                     
    //Third Column: A Coefficient with elastic form factors                                  
    //Fourth Column: B Coefficient with elastic form factors
    //Fifth Column: C Coefficient with elastic form factors                                 
    myfile3 << i << " " <<  XSXULBHDVCS << " " << bhdvcsAUL << " " << bhdvcsBUL << " " << bhdvcsCUL <<  endl;

    //Writing Polarized Data to the file XSXBHDVCSLU.dat                                     
    //First Column: phi                                                                      
    //Second Column: Total Polarized Cross Section                                           
    //Third Column: A Coefficient with elastic form factors                                  
    //Fourth Column: B Coefficient with elastic form factors
    //Fifth Column: C Coefficient with elastic form factors                                 
    myfile4 << i << " " << XSXLUBHDVCS << " " << bhdvcsALU << " " << bhdvcsBLU << " " << bhdvcsCLU <<  endl;
  }

  //Closing Files
  myfile.close();
  myfile2.close();
  myfile3.close();
  myfile4.close();

  return 0;
}
