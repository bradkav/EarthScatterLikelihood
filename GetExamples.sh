SIGMA=1e-34
SIGMASTR=-34
#LAT=46
LAT=-37
HEM=S

./calcContour 0.1 $SIGMA 1 0 $LAT test 0.7 > results/example_mx100_sig${SIGMASTR}_${HEM}_E.txt
#./calcContour 0.1 $SIGMA 1 0 $LAT test 0.8 > results/example_mx100_sig${SIGMASTR}_${HEM}_B.txt
#./calcContour 0.1 $SIGMA 1 0 $LAT test 1.0 > results/example_mx100_sig${SIGMASTR}_${HEM}_C.txt
#./calcContour 0.1 $SIGMA 1 0 $LAT test 2.0 > results/example_mx100_sig${SIGMASTR}_${HEM}_D.txt