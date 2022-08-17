#include <iostream>
#include <fstream>
#include <math.h>
#include <fftw3.h>
#include <complex>
#include <string>
#include "schrodingerLib.h"

using namespace std;

#define REAL 0
#define IMAG 1
#define PI 3.14159265359
#define SQRT_2PI 2.50662827463

extern "C" {
extern int dstev_(char*,int*,double*,double*,double*,int*,double*,int*);
}

double V(double x){
	return x*x/2;
}


double F(double x, double mu, double sigma){
	return exp(-(x-mu)*(x-mu)/(2*sigma)) * sin(3*sin(5*x)*cos(2*x)*x)*x; 
	//if(abs(x)<1) return 0.5;
	//else return 0;
}


// DO PRZENIESIENIA DO BIBLIOTEKI


/*
complex<double> *f_RK4(int NPoints, double dx, double *VPot, complex<double> *Function){
	complex<double> *Function_tmp = new complex<double>[NPoints];
	
	Function_tmp[0] = 1i*((Function[1]+Function[NPoints]-(double)2*Function[0])/(2*dx*dx) - VPot[0]*Function[0]);
	for(int i=1; i<NPoints-1; i++){
		Function_tmp[i] = 1i*((Function[i+1]+Function[i-1]-(double)2*Function[i])/(2*dx*dx) - VPot[i]*Function[i]);
	}
	Function_tmp[0] = 1i*((Function[0]+Function[NPoints-1]-(double)2*Function[NPoints])/(2*dx*dx) - VPot[NPoints]*Function[NPoints]);
	return Function_tmp;
}

complex<double> *RK4(int NPoints, double dx, double *VPot, complex<double> *Function){
	complex<double> *K1 = new complex<double>[NPoints];
	complex<double> *K2 = new complex<double>[NPoints];
	complex<double> *K3 = new complex<double>[NPoints];
	complex<double> *K4 = new complex<double>[NPoints];
	
	// K1
	K1[0] = 1i*((Function[1]+Function[NPoints]-(double)2*Function[0])/(2*dx*dx) - VPot[0]*Function[0]);
	for(int i=1; i<NPoints-1; i++){
		K1[i] = 1i*((Function[i+1]+Function[i-1]-(double)2*Function[i])/(2*dx*dx) - VPot[i]*Function[i]);
	}
	K1[0] = 1i*((Function[0]+Function[NPoints-1]-(double)2*Function[NPoints])/(2*dx*dx) - VPot[NPoints]*Function[NPoints]);
	
	// K2
	for(int i=1; i<NPoints-1; i++){
		K2[i] = 1i*((Function[i+1]+Function[i-1]-(double)2*Function[i])/(2*dx*dx) - VPot[i]*Function[i]);
	}
	Function_tmp[0] = 1i*((Function[0]+Function[NPoints-1]-(double)2*Function[NPoints])/(2*dx*dx) - VPot[NPoints]*Function[NPoints]);
	return Function_tmp;
}
*/





int main(int argc, char** argv){

// LATTICE PARAMETERS
double LMax = 4;
int NPoints = 256;

// CREATING SPATIAL LATTICE
double dx = (double)2*LMax/(NPoints-1);
double *VecX = new double[NPoints];
for(int i=0; i<NPoints; i++){
	VecX[i] = (double)dx*i - LMax;
}

// CREATING MOMENTUM LATTICE
double dp = (double) 2*PI/(NPoints*dx); 
double *VecP = new double[NPoints];
for(int i=0; i<NPoints/2; i++){
	VecP[i] = (double) i*dp;
}
for(int i=NPoints/2; i<NPoints; i++){
	VecP[i] = (double) dp*(-NPoints) + (double) i*dp;
}

// CREATING POTENTIAL VECTOR
double *Vpot = new double[NPoints];
for(int i=0; i<NPoints; i++){
	Vpot[i] = V(VecX[i]);
}

// CREATING HAMILTONIAN MATRIX
double *D = new double[NPoints];
double *E = new double[NPoints-1];

for(int i=0; i<NPoints; i++) D[i] = 1/(dx*dx) + Vpot[i];
for(int i=0; i<NPoints-1; i++) E[i] = -1/(2*dx*dx);

// LAPACK ROUTINE
char Vchar = 'V';
int LDZ = NPoints;
double *Z = new double[NPoints * LDZ];
double *work = new double[2*NPoints - 2];
int info;

dstev_(&Vchar, &NPoints, D, E, Z, &LDZ, work, &info);

if (info!=0){
cout << "Error: dstev returned error code " << info << endl;
return -1;
}

// NORMALIZING AND SETTING SIGN OF EIGENVECTORS
//sign is set in such a way that the value of the eigenfunction f(x) is positive as x->0^{+}.
int xZeroIndex = (int) ceil((double)NPoints/2);
complex<double> EigVecX[NPoints][NPoints] = {0};
double sign;
for (int i=0;i<NPoints;i++){
	for (int j=0;j<NPoints;j++){
		sign = 1.0;
		if(Z[j*NPoints+xZeroIndex] < 0) sign = -1.0;
		EigVecX[j][i] = sign * Z[j*NPoints+i]/sqrt(dx);
	}
}                             

// FOURIER TRANSFORM USING FFTW3 LIBRARY
complex<double> EigVecP[NPoints][NPoints] = {0};

for(int j=0; j<NPoints; j++){
	 FourierX(EigVecX[j], EigVecP[j], dx, NPoints);
}

// EXPORTING SPATIAL EIGENVECTORS TO FILE
FILE *fp;
fp = fopen("eigX.txt","w");
for (int i=0;i<NPoints;i++){
	fprintf(fp, "%f ", VecX[i]);
	for(int j=0; j<NPoints; j++){
		fprintf(fp, "%f ", norm(EigVecX[j][i]));
	}
	fprintf(fp, "\n");
}
fclose(fp);

// EXPORTING MOMENTUM EIGENVECTORS TO FILE
FILE *fp2;
fp2 = fopen("eigP.txt","w");
for (int i=0;i<NPoints;i++){
	fprintf(fp, "%f ", VecP[i]);
	for(int j=0; j<NPoints; j++){
		fprintf(fp, "%f ", norm(EigVecP[j][i]));
	}
	fprintf(fp, "\n");
}
fclose(fp2);


//CHECKING SOMETHING
double suma = 0;
int index = 1;
for(int i=0; i<NPoints; i++){
	suma += norm(EigVecX[index][i])*dx;
}
//cout<<"suma: "<<suma<<endl;

double sumaP = 0;
int indexP = 0;
for(int i=0; i<NPoints; i++){
	sumaP += norm(EigVecP[indexP][i])*dp;
}
//cout<<"sumaP: "<<sumaP<<endl;

// CHECKING PARSEVAL'S THEOREM 
int ind = 0;
double sumP = 0;
double sumX = 0;
for (int j=0;j<NPoints;j++){
	sumX += norm(EigVecX[ind][j])*dx;
	sumP += norm(EigVecP[ind][j])*dp; 
}
//printf("|psi(x)|^2(n=%d) = %f\n", ind, sumX);
//printf("|psi(f)|^2(n=%d) = %f\n", ind, sumP);

for (int j=0;j<NPoints;j++){
	//cout<<j<<": "<<VecX[j]<<"\t"<<VecP[j]<<endl;
}

//CREATING FUNCTION VECTOR
double *Function = new double[NPoints];
double FunctionNorm = 0;
for(int i=0; i<NPoints; i++){
	Function[i] = F(VecX[i], 0, 1);
	FunctionNorm += norm(Function[i]); 
}

double sumaF = 0;
for(int i=0; i<NPoints; i++){
	Function[i] /= sqrt(FunctionNorm*dx);
	sumaF += norm(Function[i])*dx;
}
printf("sumaF: %f\n", sumaF);

//DECOMPOSITION INTO EIGENVECTORS
complex<double> *coeffs = new complex<double>[NPoints];

for(int i=0; i<NPoints; i++){
	for(int j=0; j<NPoints; j++){
		coeffs[i] += Function[j]*conj(EigVecX[i][j])*dx;
	}
}

double coeffsum=0;
for(int i=0; i<NPoints; i++){
	
	coeffsum += norm(coeffs[i]);
	//printf("a_%d = %f\n", i, norm(coeffs[i]));
}
printf("suma wspolczynnikow: %.19f\n", coeffsum);

// EXPORTING ORIGINAL FUNCTION TO FILE
FILE *fp3;
fp3 = fopen("funcOrig.txt","w");
for (int i=0;i<NPoints;i++){
	fprintf(fp3, "%f %f\n", VecX[i], norm(Function[i]));
}
fclose(fp3);

// EXPORTING DECOMPOSED FUNCTION TO FILE
FILE *fp4;
complex<double> node = 0;
fp4 = fopen("funcDeco.txt","w");
for (int i=0;i<NPoints;i++){
	node = 0;
	for (int j=0;j<NPoints;j++){
		node += coeffs[j]*EigVecX[j][i];
	}
	fprintf(fp4, "%f %f\n", VecX[i], norm(node));
}
fclose(fp4);


// CREATING ANIMATION
FILE *fp5;

double t = 0;
double dt = 0.01;
int ok;
for(int k=0; k<(int)2*PI/dt; k++){
	char buffer [50];
	ok = sprintf(buffer, "animation/funcAnim%05d.txt", k);
	fp5 = fopen(buffer, "w");
	for (int i=0;i<NPoints;i++){
		node = 0;
		for (int j=0;j<NPoints;j++){
			node += coeffs[j]*exp(-1i*t*((double)j+(1/2)))*EigVecX[j][i];
		}
		fprintf(fp5, "%f %f\n", VecX[i], norm(node));
	}
	fclose(fp5);
	t += dt;
}

// EVOLUTION WITH RK4
complex<double> *Function_RK4 = new complex<double>[NPoints];
for(int i=0; i<NPoints; i++){
	Function_RK4[i] = Function[i];
}



// FREEING ALLOCATED MEMORY
delete [] VecX;
delete [] VecP;
delete [] Vpot;

delete [] work;
delete [] D;
delete [] E;
delete [] Z;

return 0;
}

// COMPILATION:
// g++ schrodinger.cpp -L/usr/local/lib/ -llapack -lblas -lfftw3

// Do zrobienia na 01.06.2022:
// wtorek 10 - spotkania						
// ogarnij gnuplota								W MIARE
// przynies dokumenty							NIEOK
// prawidlowa siatka P 							CHYBA OK
// dobra kolejnosc P							OK
// dobre unormowanie wektorow w P Tw Parsevala	OK
// numerical recipes - strona 520				OK
// EigVecP 										OK
// rutyna dstev() 								OK

// Do zrobienia na 14.06.2022:
// DFT jako funkcja FourierX(N, f(x), dx) (analogicznie FourierP)												OK
// Podefiniowac stałe takie jak PI, sqrt(2*PI)																	OK
// zrobic zewnetrzna biblioteke z funkcjami jak FourierX, FourierP												OK
// ogarnij jak podac funkcji jedna kolumne macierzy																OK
// dla dowolnej funkcji, np psi(x) = exp(-(x-a)^2/(2*sigma)), znajdz rozklad na funkcje wlasne operatora H		OK
// (poprzez calkowanie z poszczegolnymi funkcjami wlasnymi)
// sprawdz czy znalezione wspolczynniki zsumowane po modulokwadratowaniu dają jeden (im blizej tym lepiej)		NIEOK
// wyrysuj razem funkcje analityczną i funkcję przybliżoną po rozkladzie na funkcje wlasne						OK
// zrob ewolucje w czasie: Psi(x,t) = SUM a_i*exp(-i*E_i/hbar*t)*psi_i(x)										CHYBA OK
// zrob film z ewolucji (gnuplot albo python)																	NIEOK
// OPCJONALNIE
// zrob te samą ewolucje, ale korzystając po prostu z rownania schrodingera, metoda Runge-Kutta					NIE

// ogarnij skrypt do animacji																					OK
// ogarnij normalizacje funkcji rozkladanej																		OK
