#!/bin/bash

LOG_DIR=./LOG
DRIVERS_DIR=../DRIVERS/DOUBLE
MATRIX_FILE=../DATA/laplacian221.mtx

if [ ! -d $LOG_DIR ]; then \
    echo "mkdir $LOG_DIR"; \
    mkdir $LOG_DIR; \
fi



###### the Lanczos algorithm (w/ partial reorthogonalization) for eigenvalue computations

echo "$DRIVERS_DIR/laneig_driver $MATRIX_FILE -part=SA -nev=7 > $LOG_DIR/a0.log"
$DRIVERS_DIR/laneig_driver $MATRIX_FILE -part=SA -nev=7 > $LOG_DIR/a0.log
# this computes the 7 smallest eigenvalues
# all less than 0.5

echo "$DRIVERS_DIR/laneig_driver $MATRIX_FILE -part=SA -nev=17 > $LOG_DIR/a1.log"
$DRIVERS_DIR/laneig_driver $MATRIX_FILE -part=SA -nev=17 > $LOG_DIR/a1.log
# this computes the 17 smallest eigenvalues
# the first 7 less than 0.5, the next 9 less than 1.0, the last one greater than 1.0



###### preview the polynomial filter

echo "$DRIVERS_DIR/polyfilt_driver 10 20 0.0805287,0.5,1.0,7.91947 > $LOG_DIR/d20.dat"
$DRIVERS_DIR/polyfilt_driver 10 20 0.0805287,0.5,1.0,7.91947 > $LOG_DIR/d20.dat

echo "cd $LOG_DIR"
cd $LOG_DIR

# plot the data in "d20.dat" just generated
echo "gnuplot ../d20.gnuplot"
gnuplot ../d20.gnuplot

# view the generated eps file, assuming your machine has ghostview installed
# note that the command is "ghostview" instead of "gv" in some machines
echo "gv d20.eps &"
gv d20.eps &

# assumes that we can go back by "cd .."
echo "cd .."
cd ..



###### the filtered Lanczos algorithm (w/ partial reorthogonalization) for eigenvalue computations

echo "$DRIVERS_DIR/filtlan_driver $MATRIX_FILE -basedeg=10 -polydeg=10 -bound1=0.5 > $LOG_DIR/b0.log"
$DRIVERS_DIR/filtlan_driver $MATRIX_FILE -basedeg=10 -polydeg=10 -bound1=0.5 > $LOG_DIR/b0.log
# this computes the eigenvalues less than 1.0 (with polynomial degree 10)
# in total there are 7 eigenvalues

echo "$DRIVERS_DIR/filtlan_driver $MATRIX_FILE -basedeg=10 -polydeg=12 -bound1=1.0 > $LOG_DIR/b1.log"
$DRIVERS_DIR/filtlan_driver $MATRIX_FILE -basedeg=10 -polydeg=12 -bound1=1.0 > $LOG_DIR/b1.log
# this computes the eigenvalues less than 1.0 (with polynomial degree 12)
# in total there are 16 eigenvalues

echo "$DRIVERS_DIR/filtlan_driver $MATRIX_FILE -basedeg=10 -polydeg=20 -bound0=0.5 -bound1=1.0 > $LOG_DIR/b2.log"
$DRIVERS_DIR/filtlan_driver $MATRIX_FILE -basedeg=10 -polydeg=20 -bound0=0.5 -bound1=1.0 > $LOG_DIR/b2.log
# this computes the eigenvalues in [0.5,1.0] (with polynomial degree 20)
# in total there are 9 eigenvalues



###### estimation of number of eigenvalues

echo "$DRIVERS_DIR/numeig_driver $MATRIX_FILE 2:2:20 -basedeg=10 -polydeg=200 -bound1=0.5 > $LOG_DIR/c0.log"
$DRIVERS_DIR/numeig_driver $MATRIX_FILE 2:2:20 -basedeg=10 -polydeg=200 -bound1=0.5 > $LOG_DIR/c0.log
# this shows the result with 20 random trials
# in a test the estimated number of eigenvalues below 0.5 is 7.78 (the actual number of eigenvalues below 0.5 is 7)

echo "$DRIVERS_DIR/numeig_driver $MATRIX_FILE 2:2:20 -basedeg=10 -polydeg=200 -bound1=1.0 > $LOG_DIR/c1.log"
$DRIVERS_DIR/numeig_driver $MATRIX_FILE 2:2:20 -basedeg=10 -polydeg=200 -bound1=1.0 > $LOG_DIR/c1.log
# this shows that the result 20 random trials
# in a test the estimated number of eigenvalues below 1.0 is 17.20 (the actual number of eigenvalues below 0.5 is 16)

echo "$DRIVERS_DIR/numeig_driver $MATRIX_FILE 2:2:20 -basedeg=10 -polydeg=200 -bound0=0.5 -bound1=1.0 > $LOG_DIR/c2.log"
$DRIVERS_DIR/numeig_driver $MATRIX_FILE 2:2:20 -basedeg=10 -polydeg=200 -bound0=0.5 -bound1=1.0 > $LOG_DIR/c2.log
# this shows the result with 20 random trials
# in a test the estimated number of eigenvalues between 0.5 and 1.0 is 9.38 (the actual number of eigenvalues in this range is 9)
