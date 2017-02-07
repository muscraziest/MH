#!/bin/bash

g++ -o p5 practica5.cpp libgtest.a libarff.a -std=c++11

cd ..

mv FUENTES/p5 BIN

cd FUENTES
