#!/bin/bash

g++ -o p3 practica3.cpp libgtest.a libarff.a -std=c++11

cd ..

mv FUENTES/p3 BIN

cd FUENTES
