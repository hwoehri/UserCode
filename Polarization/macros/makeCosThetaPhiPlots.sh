#!/bin/bash

thisFrame=0
while [ $thisFrame -lt 2 ]; do

    thisPT=1
    while [ $thisPT -lt 13 ]; do 

	thisRap=1
	while [ $thisRap -lt 6 ]; do 

	    root -b -q "plotRec.C+(${thisFrame}, ${thisPT}, ${thisRap})"
	    thisRap=$(( $thisRap + 1 ))
	done
	thisPT=$(( $thisPT + 1 ))
    done
    thisFrame=$(( $thisFrame + 1 ))
done
