#!/bin/bash
source /opt/intel/compilers_and_libraries_2016/mac/bin/compilervars.sh intel64 
MODE=$1
INDIR=$2
mkdir -p dats pngs
if [ "$MODE" -eq "0" ]
then
	ITER=0
	for IN in ${INDIR}density.*
	do
		ITER=`echo "scale=2; $ITER+0.1" | bc | sed 's/^\./0./'`
		echo $ITER ps
		OUT1=particles${IN#${INDIR}density}
		OUT2=${OUT1%dat}png
		echo "./NParticles < $IN > dats/$OUT1"
		./NParticles < $IN > dats/$OUT1
		echo "gnuplot -e \"DATA='dats/$OUT1'; OUT='pngs/$OUT2'; STEP='$ITER ps'\" NParticles.plt"
		gnuplot -e "DATA='dats/$OUT1'; OUT='pngs/$OUT2'; STEP='$ITER ps'" NParticles.plt
		echo
	done
fi

if [ "$MODE" -eq "1" ]
then
	ITER=0
	for IN in dats/particles.*.dat
	do
		ITER=`echo "scale=2; $ITER+0.1" | bc | sed 's/^\./0./'`
		echo $ITER ps
		OUT1=$IN
		OUT2=${IN%dat}png
		echo "gnuplot -e \"DATA='$OUT1'; OUT='pngs/$OUT2'; STEP='$ITER ps'\" NParticles.plt"
		gnuplot -e "DATA='$OUT1'; OUT='$OUT2'; STEP='$ITER ps'" NParticles.plt
		echo
	done
fi

ffmpeg -f image2 -start_number 0001 -r 10 -i dats/particles.%04d.png -vcodec libx264 -preset ultrafast -crf 18 -pix_fmt yuv420p NParticles.mp4