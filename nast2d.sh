#! /bin/bash
#
g++ -c -Wall main.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ -c -Wall boundary.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ -c -Wall extras.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ -c -Wall init.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ -c -Wall surface.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ -c -Wall uvp.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ -c -Wall visual.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ main.o boundary.o extras.o init.o surface.o uvp.o visual.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm *.o
mv a.out ~/bincpp/nast2d
#
echo "Normal end of execution."
