#!/bin/bash
if [[ $# -eq 0 ]] ; then
    echo 'Too few arguments'
    exit 1
fi

if [ ! -d "./animation" ] 
then
    mkdir ./animation
else
	rm ./animation/*
fi

if [ ! -d "./gifs" ] 
then
    mkdir ./gifs
fi

if [ ! -d "./movies" ] 
then
    mkdir ./movies
fi

start_time=$(date +%s)
echo "Running simulation"
g++ -o schrodinger schrodinger.cpp schrodingerLib.cpp -L/usr/local/lib/ -llapack -lblas -lgfortran -lfftw3
./schrodinger
initialfunction=$(cat schrodinger.cpp | grep "double F(")
cd ./animation

echo "Creating frames"
numfiles=$(ls | wc -l)
xrange=$(cat *.txt | awk '{if ($2 > max) max=$2}END{print max}')
gnuplot -e "do for [ii=0:$numfiles-1] {set output sprintf('funcAnim%05d.png',ii); set terminal pngcairo; filename=sprintf('funcAnim%05d.txt',ii); plot [ ] [0:$xrange] filename u 1:2 w l s u notitle}"

echo "Creating gif"
convert -delay 2 -loop 0 *.png $1.gif

echo "Creating MP4 movie"
ffmpeg -loglevel warning -r 50 -i funcAnim%05d.png -y -an -pix_fmt yuv420p -metadata description="$initialfunction" $1.mp4 
mv $1.gif ../gifs/
mv $1.mp4 ../movies/
cd ..
echo "$initialfunction" > function.txt
rm -r ./animation
end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
echo "Done in $elapsed seconds."
