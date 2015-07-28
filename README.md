# TeLaPo
Tetrahedral Lattice based Polymers

TeLaPo is a command line operated program that creates polymer conformations on the fly and analyzes them. The program can generate chains ensembles of different types: non-self-avoiding and self-avoiding walks, where there exist different hierarchies of self-avoidance. 
All options are displayed when using the flags "-h", "--help", or when the input was incorrect. 

The program is divided into many files which are compiled to the main program using the MAKEFILE. If you type make telapo a executable program named telapo will be compiled. The MAKEFILE contains two lines for compiling depending on the system your running on. Comment one line out with a hash. You need the flag -lm if running on linux. For MacOS this flag needs to be removed.
