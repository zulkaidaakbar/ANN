/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
//  Calculation of the elastic form factors using the Galster parametrization  //
//                                                                             //
//  Calculates GM, GE, F1, F2, GA, and GP and from their t dependence          //
//                                                                             //
//  Written as part of the work referencing arXiv: 1903.05742                  //
//                                                                             //
//  Written by: Brandon Kriesten                                               //
//                                                                             //
//  Email: btk8bh@virginia.edu                                                 //
//                                                                             //
///////////////////////////////////////////////////////////////////////////////// 

#include <cmath>
#include "form_factor.h"

double ffGE(double t) {
  double piece = t/.710649;
  double shape = (1 + piece)*(1+piece);
  double GE = 1/shape;
  return GE;
}

double ffGM(double t) {
  double shape = ffGE(t);
  double GM0 = 2.792847337;
  return GM0*shape;
}

double ffF2(double t) { 
  double f2 = (ffGM(t) - ffGE(t))/(1-(t/(4*.938*.938)));
  return f2;
}

double ffF1(double t) { 
  double f1 = ffGM(t)- ffF2(t); 
  return f1;
}

double ffGA(double t) {
  double ga = 1.2695;
  double ma = 1.026;
  double part = t/(ma*ma);
  double dif = (1-part)*(1-part);
  double GA = ga/dif;
  return GA;
}

double ffGP(double t) {
  return 8.13;
}
