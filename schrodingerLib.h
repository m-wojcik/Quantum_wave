#ifndef _SCHRODINGERLIB_H_
#define _SCHRODINGERLIB_H_

#include <fftw3.h>
#include <complex>

using namespace std;

void FourierX(complex<double> *Fx, complex<double> *Fp, double dx, int NPoints);
void FourierP(complex<double> *Fp, complex<double> *Fx, double dp, int NPoints);

#endif 
