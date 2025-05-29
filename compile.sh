SRC=src
BIN=bin
DAT=data
TARGET=main

#NO espacios alrededor del =, para asignar valor a variables en shell
FC=gfortran
FFLAGS="-c -J$BIN -I$BIN -Wall -Wextra -std=f2008"

for file in parameters arrays functions statistics measurements main
do
  $FC $FFLAGS $SRC/$file.f90 -o $BIN/$file.o
done

$FC $BIN/*.o -o $BIN/$TARGET
$BIN/$TARGET
