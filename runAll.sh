#!/bin/bash
g++ -std=c++11 main.cpp -o simplex
for VARIABLE in {1..14}
do
    echo "LP$VARIABLE"
	./simplex lp$VARIABLE
    echo "===================================================="
done