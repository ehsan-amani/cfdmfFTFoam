#!/bin/sh

# run e.g by      
# ./run.sh 1e-06 0:5     
# parameters are $1=timeFolder $2=xRange

gnuplot logPlot.gplt
gnuplot logPlot2.gplt
gnuplot logPlot3.gplt

