#include "schrodingerLib.h"
#include <complex>

#define REAL 0
#define IMAG 1
#define PI 3.14159265359

using namespace std;

void FourierX(complex<double> *Fx, complex<double> *Fp, double dx, int NPoints){
	fftw_complex signal[NPoints];
	fftw_complex result[NPoints];
	
	for (int i = 0; i < NPoints; ++i) {
		signal[i][REAL] = Fx[i].real(); 
		signal[i][IMAG] = Fx[i].imag();       
	}
	fftw_plan plan = fftw_plan_dft_1d(NPoints, signal, result, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
	
	for (int i = 0; i < NPoints; ++i) {
		Fp[i] = (result[i][REAL] + result[i][IMAG] * 1i)*dx/sqrt(2*PI);
	}
}

	
	
void FourierP(complex<double> *Fp, complex<double> *Fx, double dp, int NPoints){
	fftw_complex signal[NPoints];
	fftw_complex result[NPoints];
	
	for (int i = 0; i < NPoints; ++i) {
		signal[i][REAL] = Fp[i].real(); 
		signal[i][IMAG] = Fp[i].imag();       
	}
	fftw_plan plan = fftw_plan_dft_1d(NPoints, signal, result, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
	
	for (int i = 0; i < NPoints; ++i) {
		Fx[i] = (result[i][REAL] + result[i][IMAG] * 1i)*dp/sqrt(2*PI);
	}
	
}
