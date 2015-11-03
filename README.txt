### miRPANDA: miRNA extension to PANDA: Passing Attributes between Networks for Data Assimilation
### the code is based on PANDA version 2 C++
### the original PANDA C++ code is available at: http://sourceforge.net/projects/panda-net/
### PANDA was published as "Passing messages between biological networks to refine predicted interactions"
### by Glass K, Huttenhower C, Quackenbush J, Yuan GC. PLoS One. 2013 May 31l8(5):e64983
### doi: 10.1371/journal.pone.0064832
### PMID: 23741402


# miRPANDA can be compiled by:
g++ miRPANDA.c -O3 -o miRPANDA

# run miRPANDA assuming all regulators can form complexes (estimate "responsibility" for all regulators)
./miRPANDA -e ToyExpressionData.txt -m ToyMotifData.txt -o ToyOutput

# run miRPANDA with a set of regulators that cannot form complexes ("miRlist.txt"). This option does not calculate the "responsibility" between these regulators and any other type of regulator.
./miRPANDA -e ToyExpressionData.txt -m ToyMotifData.txt -u miRlist.txt -o ToyOutput
