SIGMA=3.1622777e-33
SIGMASTR=-32.5
#LAT=46
LAT=-37
HEM=S

./calcContour 0.1 $SIGMA 1 0 $LAT test 0.3162 > results/example_mx100_sig${SIGMASTR}_${HEM}_A.txt
./calcContour 0.1 $SIGMA 1 0 $LAT test 1.0 > results/example_mx100_sig${SIGMASTR}_${HEM}_B.txt
./calcContour 0.1 $SIGMA 1 0 $LAT test 3.16 > results/example_mx100_sig${SIGMASTR}_${HEM}_C.txt
#./calcContour 0.1 $SIGMA 1 0 $LAT test 2.0 > results/example_mx100_sig${SIGMASTR}_${HEM}_D.txt
