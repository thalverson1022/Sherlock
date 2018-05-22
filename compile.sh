#!/bin/bash 


OUT="sherlock"
MODS="NDSearch_MODS-2"
SOURCE="NDSearch_2-0"
FLAGS=""

echo gfortran -o ${OUT}.exe ${MODS}.f90 ${SOURCE}.f90 $FLAGS
gfortran -o ${OUT}.exe ${MODS}.f90 ${SOURCE}.f90  $FLAGS
