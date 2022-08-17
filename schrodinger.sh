#!/bin/bash
if [[ $# -eq 0 ]] ; then
    echo 'Too few arguments'
    exit 1
fi

if [ ! -d "./animation" ] 
then
	echo "asdf"
    mkdir ./animation
fi

if [ ! -d "./gifs" ] 
then
    mkdir ./gifs
fi

if [ ! -d "./movies" ] 
then
    mkdir ./movies
fi


rm ./animation/*
g++ -o schrodinger schrodinger.cpp schrodingerLib.cpp -L/usr/local/lib/ -llapack -lblas -lgfortran -lfftw3
./schrodinger
cd ./animation
numfiles=$(ls | wc -l)
xrange=$(cat *.txt | awk '{if ($2 > max) max=$2}END{print max}')
gnuplot -e "do for [ii=0:$numfiles-1] {set output sprintf('funcAnim%05d.png',ii); set terminal pngcairo; filename=sprintf('funcAnim%05d.txt',ii); plot [ ] [0:$xrange] filename u 1:2 w l s u notitle}"
convert -delay 2 -loop 0 *.png $1.gif
ffmpeg -r 50 -i funcAnim%05d.png -y -an -pix_fmt yuv420p $1.mp4
mv $1.gif ../gifs/
mv $1.mp4 ../movies/
cd ..
rm -r ./animation
