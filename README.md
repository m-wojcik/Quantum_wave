# Quantum_wave
A program for various operations related to the quantum wavefunction in a given potential. It operates on a discrete lattice in a finite spatial domain.

# Requirements
1. LAPACK
2. BLAS
3. Gfortran
4. FFTW3
5. Gnuplot
6. FFMPEG
7. ImageMagick

LAPACK, BLAS, and Gfortran can be downloaded from: https://netlib.org/lapack/

FFTW3 can be downloaded from: http://fftw.org/download.html

Gnuplot and FFMPEG can easily be downloaded via command prompt.

ImageMagick can be downloaded from: https://imagemagick.org/script/download.php
 
# Manual
1. Copy files from the repository to a chosen location. 
2. Modify contents of `schrodinger.cpp` as you wish. In particular you might want to modify the function `double F(...)` in the beggining of the file. This function provides the initial shape for the wavefunction. The function `double V(double x)` defines the potential in a quantum system, and by default its a simple harmonic potential: $V(x) = \dfrac{x^2}{2}$.
3. Open command prompt and move to a directory with the copied files.
4. Type: `bash schrodinger.sh FILENAME` with your chosen filename, without any file extension.

This creates five text files:
- `EigX.txt` with eigenfunctions of the Schrodinger equation for a given potential.
- `EigP.txt`the same eigenfunctions, but in momentum domain (Fourier transform is done using FFTW3).
- `funcOrig.txt` original  numerical representation of the initial wavefunction.
- `funcDeco.txt` initial wavefunction decomposed using the eigenfunctions. It should be similar to the original wavefunction.
- `function.txt` one-liner text file containing the formula for the initial wavefunction from the latest simulation.

In addition, it creates (if not already present) two directories: `movies` and `gifs`, and creates a `.MP4` file and a `.gif` of the evolution in time of the given wavefunction, and puts them in their respectful folders. The `.MP4` file contains a metadata comment in which the formulat for the initial wavefunction is preserved.
