#!/bin/bash

g++ -O2 -Wno-deprecated $1.cc -o $1 -pthread -I$ROOTSYS/include -L$ROOTSYS/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lGui -lMinuit -pthread -lm -ldl -rdynamic
