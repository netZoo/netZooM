#PUMA
##PANDA Using MicroRNA Associations
PUMA, or **P**ANDA **U**sing **M**icroRNA **A**ssociations, is an extension of the gene regulatory network reconstruction algorithm PANDA, which was published in "Passing messages between biological networks to refine predicted interactions" by Glass K, Huttenhower C, Quackenbush J, Yuan GC. *PLoS One*. 2013 May 31l8(5):e64983, doi: 10.1371/journal.pone.0064832, PMID: 23741402

PUMA can reconstruct gene regulatory networks using both transcription factors and microRNAs as regulators of mRNA expression levels. This is the PUMA C++ code, which is based on PANDA version 2. The original PANDA C++ code is available at: http://sourceforge.net/projects/panda-net/


PUMA can be compiled using:
```
g++ PUMA.c -O3 -o PUMA
```

To run PUMA with the assumption that all regulators can form complexes (estimate *responsibility* for all regulators, *eg* a gene regulatory prior with transcription factors only):
```
./PUMA -e ToyExpressionData.txt -m ToyMotifData.txt -o ToyOutput
```

To tell PUMA to discriminate between regulators that can, and regulators that cannot cannot form complexes (for example a list of microRNAs in `miRlist.txt`), run:
```
./PUMA -e ToyExpressionData.txt -m ToyMotifData.txt -u miRlist.txt -o ToyOutput
```
