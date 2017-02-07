#!/bin/bash

g++ -o p1 practica1.cpp libgtest.a libarff.a -std=c++11

cd ..

mv FUENTES/p1 BIN

cd FUENTES
