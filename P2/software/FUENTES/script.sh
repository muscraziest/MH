#!/bin/bash

g++ -o p2 practica2.cpp libgtest.a libarff.a -std=c++11

cd ..

mv FUENTES/p2 BIN

cd FUENTES
