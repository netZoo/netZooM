#PUMA

or **P**ANDA **U**sing **M**icroRNA **A**ssociations, is an extension of PANDA (**P**assing **A**ttributes between **N**etworks for **D**ata **A**ssimilation. PUMA can reconstruct gene regulatory networks using both transcription factors and microRNAs as regulators of mRNA expression levels. This is the PUMA C++ code, which is based on PANDA version 2. The original PANDA C++ code is available at: http://sourceforge.net/projects/panda-net/
PANDA was published as "Passing messages between biological networks to refine predicted interactions" by Glass K, Huttenhower C, Quackenbush J, Yuan GC. PLoS One. 2013 May 31l8(5):e64983, doi: 10.1371/journal.pone.0064832, PMID: 23741402

PUMA can be compiled using:
'''
g++ PUMA.c -O3 -o PUMA
'''

Run PUMA assuming all regulators can form complexes (estimate "responsibility" for all regulators)
'''
./PUMA -e ToyExpressionData.txt -m ToyMotifData.txt -o ToyOutput
'''

Run PUMA with a set of regulators that cannot form complexes ("miRlist.txt"). This option does not calculate the "responsibility" between these regulators and any other type of regulator.
'''
./PUMA -e ToyExpressionData.txt -m ToyMotifData.txt -u miRlist.txt -o ToyOutput
'''
