///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Calculates four vector products and the transverse four vector products  //
//                                                                           //
//  References work: 1903.05742                                              //
//                                                                           //
//  Written by: Brandon Kriesten                                             //
//                                                                           //
//  Email: btk8bh@virginia.edu                                               //
//                                                                           //
/////////////////////////////////////////////////////////////////////////////// 

#include "four_vector.h"
#include <iostream>

double product(double arr[], double arr2[]) {

  double sum;

  sum = arr[0]*arr2[0] - arr[1]*arr2[1] - arr[2]*arr2[2] - arr[3]*arr2[3];

  return sum;
}

double tproduct(double arr[], double arr2[]) {

  double tprod;
  
  tprod = arr[1]*arr2[1] + arr[2]*arr2[2];

  return tprod;
}
