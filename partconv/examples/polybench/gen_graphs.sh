#/bin/bash

#this script takes the name of a bench without the extension (like jacobi-1d-imper) as argument
#and generates two graphs : NAME.gra with the original schedule, and NAME.tile.gra tiled with polly

if ! [ -e polybench.o ]
then
   clang -c -I. -flto -w polybench.c
fi

echo "compiling executable"
clang -O0 -c -I. -flto -w -DSMALL_DATASET "${1}.c"
llvm-link "${1}.o" polybench.o -o ${1}.ll
echo "instrumentation"
opt -load libddg-instr.so -mem2reg -indvars \
  -instrument-ddg -instrument-indvars ${1}.ll -o ${1}.instr.ll
llc ${1}.instr.ll  -o ${1}.instr.s
clang -lddg-rt -lm ${1}.instr.s -o ${1}.bin
echo "generating trace"
export LD_LIBRARY_PATH=~/llvm-3.5.2/lib/ #change this for your installation
./${1}.bin > ${1}.out 2> ${1}.log
echo "printing graph"
dynamic-graph -noPhi ${1}.instr.ll
mv dynamic-graph.txt ${1}.gra

