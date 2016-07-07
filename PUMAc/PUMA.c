/*

Changed:
#define MAXTFS 1500    # was 1000
#define MAXMIRS 1500   # was 1000

*OBS! this means marieke updated/changed this line

Version 1 Modifications (May 2013):
1) added "randomseed" variable which allows the user to specific the random number generator seed
("srand(randomseed)").  This is useful when doing paired randomizations (e.g. if one wants the gene labels to
be swapped the same way for two different sets of input data).
2) added in a second method of randweights (randweight=2 set by specifying -r 3 at the command prompt).  This
generates a random-weight value for any covariate weight that was initially greater than zero, but leaves
zero weights unchanged.
3) Removed criteria that a protein in the PPI must be a member of the regulatory prior.  Instead add in any
"new" proteins into the regulatory prior assuming no known regulatory interactions.
4) Removed criteria that the TF/motif in the regulatory prior must also be a gene in the expression data. 
This allows regulators to take names that aren't gene names (e.g. a regulator could be TAL1::GATA1, but the
genes are in RefSeq annotation).  One limitation is that correlation in expresion levels between "TFs" and
genes is no longer calculated.  This information, however, was never used by PANDA, so removing the
calculation had the added benefit of freeing up memory.
5) added in another verboseoutput option.  Now setting "-v 2" at the command prompt will cause PANDA to print
out additional files recording the initial and final protein-interaction and co-regulatory networks.  Changed
the behavior of the code such that an initial regulatory network is only printed out when using this option.
6) Increased the number of Regulators allowed by Program to 1000.
7) Modified function that reads in regulatory and co-regulating information such that it can handle
multiple instances.  Single instances in the input files will result in the initial value being set equal to
the value in the "weight" column.  If an interaction is multiply defined, the initial value will be set
equal to the sum of values associated with these instances in the "weight" column.  Undefined instances
are given a default value of 0.

Version 2 Modifications (July 2014):
1) added in "LeaveOutSample" to leave out a single sample from the network reconstruction.
2) fixed missing string termination that could cause the terminal window to become bold.
3) added in the "JackKnife" option to designate number of samples to use in a jack-knife network.
4) modified the ReadInExpression function to allow users to add additional rows to their expression file (likely
header rows), so long as the first character in these rows is a hashtag (#).
5) Changed length of "TF" in regulation struct to be 64 characters, in preparation to longer miRNA names.
6) Created default value for outtag so that the -o command-line option is now optional instead of required.
7) Defined MAXGENES, MAXTFS, MAXCONDITIONS, BUFSIZE, and MAXPATHLENGTH to allow easier manipulation of these
values should they need to be altered.
8) Added some additional outputs that tell the user what the code is doing.
9) Put in catch when normalizing the prior for the case where a TF and its potential target gene both have no
targets/inputs (variance of 0 in the prior).
10) removed index, indegree and outdegree parts of the genes and regulation structures as they were not being used.
11) Restructured code so that there are much fewer global variables and most are declared locally.
12) Modified ReadInExpression function to sum over multiple entries of a gene in the expression file.
13) Added in program exits triggers for if input files have more than MAXGENES MAXTFS or MAXCONDITIONS.

Potential Future improvements:
* change the Genes.corr and Regulation.P from symmetrix matrices to vectors to save on memory space (especially
for the former)
* remove the exp and stderr portions of the Genes and Regulation structs and make them local vectors (currently
they are only used to normalize the initial PPI and corr matrices and nothing else).
* parralize the code by adding options for multi-threading for-loops using the openmp library, shoud enhance speed.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <signal.h>
#include <string.h>

#define CR 13
#define LF 10
#define BOLD "\033[1m"
#define NORMAL "\033[0m"

#define MAXGENES 20000
#define MAXTFS 1500
#define MAXMIRS 1500 // *OBS!
#define MAXCONDITIONS 500
#define BUFSIZE 10000
#define MAXPATHLENGTH 500

FILE *fid;
FILE *statfid;
char temp[BUFSIZE];
char temp1[32];
char * temp2;

char covariate_file[MAXPATHLENGTH];
char output_file[MAXPATHLENGTH];
char outtag[MAXPATHLENGTH];
char pname[MAXPATHLENGTH];

// global variables 
float alpha;
int killstep;
int randweights;
int verboseoutput;
int weightedpearson;
int randomseed;
int NumGenes;
int NumTFs;
int NumTFsm;	// *OBS! number of TFs
int NumMiRs;	// *OBS! number of miRs
int NumConditions;
int NumInteractions, NumUniqueInteractions;
int noExp;
int LeaveOutSample;
int JackKnife;
int nomirlist;	// *OBS!

// structs to hold information about the networks

typedef struct{
char name[64];			// name of regulator
float W[MAXGENES];		// Current Weight of edge
float M[MAXGENES];		// Prior knowledge from chip-chip/motif data
float P[MAXTFS];		// PPI data
float T[MAXGENES][2];		// Temporary storage of estimates
double exp;
double stdev;
} REGULATION;
REGULATION Regulation[MAXTFS];

typedef struct{		// *OBS! start
char name[64];		
} MIR;
MIR MiR[MAXMIRS];	// *OBS! end

typedef struct{
char name[32];
char temp[32];
float expression[MAXCONDITIONS];
double exp;			// expected value of expression across conditions
double stdev;			// standard deviation of expression across conditions
float corr[MAXGENES];
} GENES;
GENES Genes[MAXGENES];

// functions
int Initialize(REGULATION *reg, GENES *gen);
int ReadInExpressionData(GENES *gen, char filename[]);
int ReadInPriorData(REGULATION *reg, GENES *gen, char filename[]);
int ReadInMiRlist(MIR *MiR, char filename[]); 		// *OBS! function to read in miR list
int CompareMiRwithReg(MIR *MiR, REGULATION *reg); 	// *OBS! function to compare miRs with all regulators
int ReadInInteractionData(REGULATION *reg, GENES *gen, char filename[]);
int NormalizePriorData(REGULATION *reg, GENES *gen);
int Correlation(REGULATION *reg, GENES *gen);
int IdentityPPI(REGULATION *reg, GENES *gen);
int IdentityCorrelation(GENES *gen);

float LearnNetwork(REGULATION *reg, GENES *gen, int step);
float UpdateCorrelation(REGULATION *reg, GENES *gen, int step);
float UpdatePPI(REGULATION *reg, GENES *gen, MIR *MiR, int step); // *OBS!

int PrintStats(REGULATION *reg, GENES *gen, char filename[]);
int PrintCoReg(REGULATION *reg, GENES *gen, char filename[]);
int PrintPPI(REGULATION *reg, GENES *gen, char filename[]);

double CDF(double value);
double inverseCDF(double value);
int RandPerm(GENES *gen);
int sort(const void *x, const void *y){if(*(const float*)y < *(const float*)x) return -1; return *(const float*)y > *(const float*)x;}
void SignalHandler(int signum);
void useage();

int main(int argc, char *argv[])
{
	// set program defaults for global variables
	strcpy(outtag, "PANDA_prediction");
	NumGenes=6000;
	NumTFs=1000;
	NumTFsm=1000; 	// *OBS!
	NumMiRs=1000; 	// *OBS!
	NumConditions=500;
	alpha=0.1;
	noExp=1;
	LeaveOutSample=0;
	JackKnife=0;
	randweights=0;
	verboseoutput=0;
	weightedpearson=0;
	randomseed=0;
	nomirlist=1; 	// *OBS! default is run PANDA without list of miRs

	// set local variables
	int maxstep=1000;
	int outputstep=0;
	int noPPI=1;
	int noMotif=1;
	int randlabels=0;
	float hamming;

	// variables to store file names
	char interaction_file[MAXPATHLENGTH];
	char expression_file[MAXPATHLENGTH];
	char motif_file[MAXPATHLENGTH];
	char mirlist_file[MAXPATHLENGTH]; 	// *OBS!

	extern char *optarg;
	int errflg=2;
	int s;
	strcpy(pname, argv[0]); 

	while((s = getopt(argc, argv, "a:e:m:o:p:k:l:j:n:r:s:u:v:w:")) != -1) 	// *OBS! added u
	switch(s)
	{
		case 'e':	// file name of expression data 
			strcpy(expression_file, optarg); errflg--;
			break;
		case 'a':	// set value for alpha (default 0.1)
			alpha=atof(optarg);
			break;
		case 'm':	// pair of motif edges 
			strcpy(motif_file, optarg); errflg--; noMotif=0;
			break;
		case 'p':	//pair file of PPI
			strcpy(interaction_file, optarg);  noPPI=0;
			break;
		case 'o':	// tag for output files
			strcpy(outtag, optarg);
			break;
		case 'k':
			maxstep=atoi(optarg);
			break;
		case 'n':
			outputstep=atoi(optarg);
			break;
		case 'r':
			randlabels=atoi(optarg);
			break;
		case 's':
			randomseed=atoi(optarg);
			break;
		case 'l':
			LeaveOutSample=atoi(optarg);
			break;
		case 'j':
			JackKnife=atoi(optarg);
			break;
		case 'u': 	// file name of miR list		// *OBS!
			strcpy(mirlist_file, optarg);  nomirlist=0; 	// *OBS!
			break;			 			// *OBS!
		case 'w':
			weightedpearson=1;
			strcpy(covariate_file, optarg);
			break;
		case 'v':
			verboseoutput=atoi(optarg);
			break;
		default:
			errflg++;
	}
	
	if (errflg)
	{
		useage();
		exit (2);
	}

	signal(SIGINT, SignalHandler);
	if(randomseed==0) randomseed=time(NULL);

	fprintf(stderr, "Reading in data.  Note that the code only allocates space for up to %u Regulators, %u Target Genes, and %u expression conditions.  If you have more than this you may need to alter the code to avoid a segfault.\n", MAXTFS, MAXGENES, MAXCONDITIONS);
	// Step (0): Initialize Values
	Initialize((REGULATION *) &Regulation, (GENES *) &Genes);
	ReadInExpressionData((GENES*) &Genes, expression_file);
	if(randlabels==1) RandPerm((GENES*) &Genes);
	if(randlabels==2) randweights=1;
	if(randlabels==3) randweights=2;
	ReadInPriorData((REGULATION*) &Regulation, (GENES*) &Genes, motif_file);	

	// Report what program thinks is going on
	fprintf(stderr, "Running PANDA using the following parameters:\n");
	fprintf(stderr, "alpha=%f\n", alpha);
	fprintf(stderr, "Data-Types Being Used include:\n");
	if(noMotif==0) fprintf(stderr, "Regulation Data\n");
	if(nomirlist==0) fprintf(stderr, "List of miRs\n"); 	// *OBS!
	if(noPPI==0) fprintf(stderr, "Protein Interaction Data\n");
	if(noExp==0) fprintf(stderr, "Expression Data\n");
	
	// Read in data-types
	if(noPPI==0) ReadInInteractionData((REGULATION*) &Regulation, (GENES*) &Genes, interaction_file);
	else IdentityPPI((REGULATION*) &Regulation, (GENES *) &Genes);
	if(noExp==0) Correlation((REGULATION *) &Regulation, (GENES *) &Genes);
	else IdentityCorrelation((GENES *) &Genes);

	if(nomirlist==0) { 	// *OBS! start
		ReadInMiRlist((MIR*) &MiR, mirlist_file);	
		CompareMiRwithReg((MIR*) &MiR, (REGULATION*) &Regulation);
	}
	else NumMiRs = 0; 	// *OBS! end

	NumTFsm=NumTFs-NumMiRs; // *OBS! NumTFsm = number of TFs (number of regulators - number of miRs)

	if(verboseoutput==2)
	{
		sprintf(output_file, "%s_InitialNetwork.pairs", outtag);
		PrintStats((REGULATION *) & Regulation, (GENES *) & Genes, output_file);
	}
	NormalizePriorData((REGULATION*) &Regulation, (GENES *) &Genes);
	
	fprintf(stderr, "\nLearning Network!\n");
	// Learn Network!
	hamming=NumTFs*NumGenes;
	killstep=0;
	while(hamming>1e-3 && killstep<=maxstep)
	{
		hamming=0;
		// Step (1): Learn Network
		hamming=LearnNetwork((REGULATION *) &Regulation, (GENES *) &Genes, killstep);
		
		// Step (2): Update Correlation
		UpdateCorrelation((REGULATION *) &Regulation, (GENES *) &Genes, killstep);
		UpdatePPI((REGULATION *) &Regulation, (GENES *) &Genes, (MIR *) &MiR, killstep); // *OBS!
		
		if(outputstep>0 && killstep % outputstep == 0)
		{
			sprintf(output_file, "%s_Step%u.stats", outtag, killstep);
			PrintStats((REGULATION *) & Regulation, (GENES *) & Genes, output_file);
		}
		fprintf(stderr, "Step#%u, hamming=%f\n", killstep, hamming);
		killstep++;
	}
	sprintf(output_file, "%s_FinalNetwork.pairs", outtag);
	PrintStats((REGULATION *) & Regulation, (GENES *) & Genes, output_file);
	
	if(verboseoutput==2)
	{
		sprintf(output_file, "%s_FinalCoReg.pairs", outtag);
		PrintCoReg((REGULATION *) & Regulation, (GENES *) & Genes, output_file);
		sprintf(output_file, "%s_FinalPPI.pairs", outtag);
		PrintPPI((REGULATION *) & Regulation, (GENES *) & Genes, output_file);
	}
	return 0;
}


float LearnNetwork(REGULATION *reg, GENES *gen, int step)
{
	// local variables
	float ExpUpdate, PPIUpdate;
	float temphamming;
	double ExpMean, ExpStd, PPIMean, PPIStd;
	int cnt, cnt1, cnt2, d;
	double A, B, C;

	ExpMean=0;
	ExpStd=0;
	for(cnt1=0; cnt1<NumTFs; cnt1++)
	{
		for(cnt2=0; cnt2<NumGenes; cnt2++)
		{
			A=0;
			B=0;
			C=0;
			d=0;
			for(cnt=0; cnt<NumGenes; cnt++)
			{
				A+=gen[cnt2].corr[cnt]*reg[cnt1].W[cnt];
				B+=gen[cnt2].corr[cnt]*gen[cnt2].corr[cnt];
				C+=reg[cnt1].W[cnt]*reg[cnt1].W[cnt];
				d++;
			}
			A=A/sqrt(B+C-fabs(A));
			
			ExpMean+=A;
			ExpStd+=A*A;
			reg[cnt1].T[cnt2][0]=A;
		}
	}
	ExpMean=ExpMean/(NumTFs*NumGenes);
	ExpStd=sqrt(ExpStd/(NumGenes*NumTFs)-ExpMean*ExpMean);
	if(verboseoutput==1) fprintf(stderr, "%f,%f;", ExpMean, ExpStd);

	PPIMean=0;
	PPIStd=0;
	for(cnt2=0; cnt2<NumGenes; cnt2++)
	{
		for(cnt1=0; cnt1<NumTFs; cnt1++)
		{
			A=0;
			B=0;
			C=0;
			d=0;
			for(cnt=0; cnt<NumTFs; cnt++)
			{
				A+=reg[cnt].P[cnt1]*reg[cnt].W[cnt2];
				B+=reg[cnt].P[cnt1]*reg[cnt].P[cnt1];
				C+=reg[cnt].W[cnt2]*reg[cnt].W[cnt2];

			}
			A=A/sqrt(B+C-fabs(A));

			PPIMean+=A;
			PPIStd+=A*A;
			reg[cnt1].T[cnt2][1]=A;
		}
	}
	PPIMean=PPIMean/(NumTFs*NumGenes);
	PPIStd=sqrt(PPIStd/(NumGenes*NumTFs)-PPIMean*PPIMean);
	if(verboseoutput==1) fprintf(stderr, "%f,%f;", PPIMean, PPIStd);

	ExpUpdate=0.5;
	PPIUpdate=0.5;

	temphamming=0;
	for(cnt1=0; cnt1<NumTFs; cnt1++)
	{
		for(cnt=0; cnt<NumGenes; cnt++)
		{
			reg[cnt1].T[cnt][0]=ExpUpdate*reg[cnt1].T[cnt][0]+PPIUpdate*reg[cnt1].T[cnt][1];
			temphamming+=fabs(reg[cnt1].W[cnt]-reg[cnt1].T[cnt][0]);
			reg[cnt1].W[cnt]=(1-alpha)*reg[cnt1].W[cnt]+alpha*reg[cnt1].T[cnt][0];
		}
	}
	return temphamming/(NumTFs*NumGenes);
}

float UpdateCorrelation(REGULATION *reg, GENES *gen, int step)
{
	// local variables
	int cnt, cnt1, cnt2, c;
	double A, B, C;
	double CorrMean, CorrStd;

	for(cnt=0; cnt<NumGenes; cnt++) {gen[cnt].exp=0; gen[cnt].stdev=0;}

	CorrMean=0;
	CorrStd=0;
	c=0;
	for(cnt=0; cnt<NumGenes; cnt++)
	{
		for(cnt2=cnt+1; cnt2<NumGenes; cnt2++)
		{
			A=0;
			B=0;
			C=0;
			for(cnt1=0; cnt1<NumTFs; cnt1++)
			{
				A+=reg[cnt1].W[cnt2]*reg[cnt1].W[cnt];
				B+=reg[cnt1].W[cnt]*reg[cnt1].W[cnt];
				C+=reg[cnt1].W[cnt2]*reg[cnt1].W[cnt2];
			}
			A=A/sqrt(B+C-fabs(A));
			gen[cnt].corr[cnt2]=A;
			
			CorrMean+=A;
			CorrStd+=A*A;
			c++;

			gen[cnt].exp+=A;
			gen[cnt2].exp+=A;
			gen[cnt].stdev+=A*A;
			gen[cnt2].stdev+=A*A;
		}
	}

	CorrMean=CorrMean/c;
	CorrStd=sqrt(CorrStd/c-CorrMean*CorrMean);
	// if(verboseoutput==1) fprintf(stderr, "%f,%f;", CorrMean, CorrStd);

	c=0;
	for(cnt=0; cnt<NumGenes; cnt++)
	{
		gen[cnt].exp=gen[cnt].exp/(NumGenes-1);
		gen[cnt].stdev=gen[cnt].stdev/(NumGenes-1)-gen[cnt].exp*gen[cnt].exp;
		gen[cnt].exp=(NumGenes)*(sqrt(gen[cnt].stdev))*exp(2*alpha*((float) step));
		gen[cnt].corr[cnt]=(1-alpha)*gen[cnt].corr[cnt]+alpha*gen[cnt].exp;
		for(cnt2=cnt+1; cnt2<NumGenes; cnt2++)
		{
			c++;
			gen[cnt].corr[cnt2]=(1-alpha)*gen[cnt2].corr[cnt]+alpha*(gen[cnt].corr[cnt2]);
			gen[cnt2].corr[cnt]=gen[cnt].corr[cnt2];
		}
	}
	return 0;
}

float UpdatePPI(REGULATION *reg, GENES *gen, MIR *mir, int step) 	// *OBS!
{
	// local variables
	double PPIMean, PPIStd;
	int cnt, cnt1, cnt2, c, d; 	// *OBS! d added
	double A, B, C;

	for(cnt1=0; cnt1<NumTFs; cnt1++) {reg[cnt1].exp=0; reg[cnt1].stdev=0;}
	
	PPIMean=0;
	PPIStd=0;
	c=0;
	for(cnt1=0; cnt1<NumTFs; cnt1++)
	{
		for(cnt=cnt1+1; cnt<NumTFs; cnt++)
		{
			A=0;
			B=0;
			C=0;
			for(cnt2=0; cnt2<NumGenes; cnt2++)
			{
				A+=reg[cnt1].W[cnt2]*reg[cnt].W[cnt2];
				B+=reg[cnt1].W[cnt2]*reg[cnt1].W[cnt2];
				C+=reg[cnt].W[cnt2]*reg[cnt].W[cnt2];
			}
			A=A/sqrt(B+C-fabs(A));
			
			PPIMean+=A;
			PPIStd+=A*A;
			c++;

			// *OBS! i've commented out the next line, i think cnt==cnt1 does not happen
			// if(cnt==cnt1) reg[cnt].P[cnt]=(1-alpha)*reg[cnt].P[cnt]+A*alpha;
		// *OBS! start
			for(d=0; d<NumMiRs; d++)
			{
				if(strcmp(reg[cnt].name,mir[d].name)==0) break;
				if(strcmp(reg[cnt1].name,mir[d].name)==0) break;
			}
			if (d==NumMiRs) reg[cnt1].P[cnt]=A; // *OBS! this was // else reg[cnt1].P[cnt]=A;
		// *OBS! end
			reg[cnt1].exp+=A;
			reg[cnt].exp+=A;
			reg[cnt1].stdev+=A*A;
			reg[cnt].stdev+=A*A;
		}
	}
	PPIMean=PPIMean/c;
	PPIStd=sqrt(PPIStd/c-PPIMean*PPIMean);
	// if(verboseoutput==1) fprintf(stderr,"%f,%f;", PPIMean, PPIStd);

	for(cnt1=0; cnt1<NumTFs; cnt1++)
	{
		reg[cnt1].exp=reg[cnt1].exp/(NumTFs-1);
		reg[cnt1].stdev=reg[cnt1].stdev/(NumTFs-1)-reg[cnt1].exp*reg[cnt1].exp;
		reg[cnt1].exp=(NumTFs)*(sqrt(reg[cnt1].stdev))*exp(2*alpha*((float) step));
		reg[cnt1].P[cnt1]=(1-alpha)*reg[cnt1].P[cnt1]+alpha*reg[cnt1].exp;
		for(cnt=cnt1+1; cnt<NumTFs; cnt++)
		{
			reg[cnt].P[cnt1]=(1-alpha)*reg[cnt].P[cnt1]+alpha*reg[cnt1].P[cnt];
			reg[cnt1].P[cnt]=reg[cnt].P[cnt1];
		}
	}
	return 0;
}

int Initialize(REGULATION *reg, GENES *gen)
{
	int cnt,cnt1,cnt2;

	for(cnt1=0; cnt1<NumTFs; cnt1++)
	{
		for(cnt2=0; cnt2<NumGenes; cnt2++)
		{
			reg[cnt1].M[cnt2]=0;
			reg[cnt1].W[cnt2]=0;
			reg[cnt1].T[cnt2][0]=0;
			reg[cnt1].T[cnt2][1]=0;
		}
		for(cnt=0; cnt<NumTFs; cnt++) reg[cnt1].P[cnt]=0;
	}
	return 0;
}

int ReadInExpressionData(GENES *gen, char filename[])
{
	int cnt, cnt2;
	int gcnt, c;

	if((fid=fopen(filename, "r"))==NULL)
	{
		printf("ERROR OPENING EXPRESSION DATA FILE\n");
		exit(1);
	}

	fprintf(stderr, "Reading In Expression Data!\n");

	cnt=0;
	cnt2=0;
	NumGenes=0;
	while(!feof(fid))
	{
		if(fgets(temp, BUFSIZE , fid) != NULL)
		{
			if(strncmp(temp,"#",1)>0)
			{
				temp2 = strtok (temp,"\t\n");
				cnt2++;
				gcnt=NumGenes;
				for(c=0; c<NumGenes; c++)
				{
					if(strcmp(temp2,gen[c].name)==0)
					{
						gcnt=c;
						break;
					}
				}
				if(gcnt==NumGenes)
				{
					cnt=0;
					strcpy(gen[gcnt].name, temp2);
					temp2 = strtok (NULL,"\t\n");
	  				while (temp2 != NULL)
  					{
						gen[gcnt].expression[cnt]=atof(temp2);
						cnt++;
						temp2 = strtok (NULL,"\t\n");
					}
					if(gcnt==0)
					{
						NumConditions=cnt;
						if(NumConditions>MAXCONDITIONS)
						{
							printf("TOO MANY CONDITIONS DETECTED IN EXPRESSION FILE. CODE UPDATE NECESSARY TO AVOID SEG-FAULT.\n");
							exit(1);
						}
					}
					NumGenes++;
					if(NumGenes>MAXGENES)
					{
						printf("TOO MANY GENES DETECTED. CODE UPDATE NECESSARY TO AVOID SEG-FAULT.\n");
						exit(1);
					}

				}
				else
				{
					temp2 = strtok (NULL,"\t\n");
  					while (temp2 != NULL)
  					{
						gen[gcnt].expression[cnt]+=atof(temp2);
						temp2 = strtok (NULL,"\t\n");
					}
				}
			}
		}
	}
	NumGenes=NumGenes;
	fclose(fid);
	fprintf(stderr, "Num Unique Genes in Expression File:%u (%u total entries).\n", NumGenes, cnt2);

	if(NumConditions>2) noExp=0;
	else
	{
		NumConditions=0;
		fprintf(stderr, "Note that an insufficient number of conditions have been found to utilize Expression Data.  Using input file as a gene list.\n");
	}
	if(NumGenes>6000)
	{
		fprintf(stderr, "Data includes information for more than 6000 genes, allocating more space\n");
		Initialize((REGULATION *) &Regulation, (GENES *) &Genes);
	}
	return 0;
}

// *OBS! start
int ReadInMiRlist(MIR *mir, char filename[]) // reads in list of miRs
{

	int cnt, c;
	NumMiRs=0;

	if((fid=fopen(filename, "r"))==NULL)
	{
		printf("ERROR OPENING MIRNA LIST FILE\n");
		exit(1);
	}

	fprintf(stderr, "Reading in the miRNA list\n");

	while(fscanf(fid, "%s", temp)==1) //!feof(fid))
	{
		cnt=NumMiRs;

		for(c=0; c<NumMiRs; c++)
		{
			if(strcmp(temp,mir[c].name)==0)
			{
				cnt=c;
				break;
			}
		}

		if(cnt==NumMiRs)
		{
			strcpy(mir[cnt].name, temp);
			NumMiRs++;
		}

	}
	fclose(fid);
	return 0;
}

int CompareMiRwithReg(MIR *mir, REGULATION *reg) // compares list of miRs with regulators
{

	int c, d;
	for(c=0; c<NumMiRs; c++)
	{
		int d=0;
		for(d==0; d<NumTFs; d++)
		{
			if(strcmp(mir[c].name,reg[d].name)==0)
			{
			break;
			}
		}
		if(d==NumTFs)
		{
			fprintf(stderr, "miRNA '%s' found in mirlist but not Regulation Data.  Please verify your input files.\n", mir[c].name);

			exit(1); // *OBS! this should be changed to: add miR to Regulation Data without edges
		}
	}
}
// *OBS! end

int ReadInPriorData(REGULATION *reg, GENES *gen, char filename[])
{
	int cnt, cnt2, c;
	float P;

	NumTFs=0;
	NumInteractions=0;
	NumUniqueInteractions=0;

	if((fid=fopen(filename, "r"))==NULL)
	{
		printf("ERROR OPENING TF-GENE FILE\n");
		exit(1);
	}

	fprintf(stderr, "Reading in Motif Data!\n");

	while(fscanf(fid, "%s\t%s\t%f", temp,temp1,&P)==3) //!feof(fid))
	{
		cnt=NumTFs;
		for(c=0; c<NumTFs; c++)
		{
			if(strcmp(temp,reg[c].name)==0)
			{
				cnt=c;
				break;
			}
		}

		cnt2=NumGenes;
		for(c=0; c<NumGenes; c++)
		{
			if(strcmp(temp1,gen[c].name)==0)
			{
				cnt2=c;
				break;
			}
		}
		if(cnt2==NumGenes)
		{
			if(strcmp(temp1, "none")==0) cnt2=-1;
			else
			{
				fprintf(stderr, "WARNING: Gene '%s' found in Regulation Data but not Gene Expression Data.  It will be added with assumed equal expression across all conditions.  If this is an error, exit the program and verify your input files.\n", temp1);
				strcpy(gen[cnt2].name, temp1);
				for(c=0; c<NumConditions; c++) gen[cnt2].expression[c]=0;
				for(c=0; c<200; c++)
				{
					reg[c].M[cnt2]=0;
					reg[c].W[cnt2]=0;
					reg[c].T[cnt2][0]=0;
					reg[c].T[cnt2][1]=0;
				}
				NumGenes++;
				if(NumGenes>MAXGENES)
				{
					printf("TOO MANY GENES DETECTED. CODE UPDATE NECESSARY TO AVOID SEG-FAULT.\n");
					exit(1);
				}
			}
		}
		
		if(cnt==NumTFs)
		{
			strcpy(reg[cnt].name, temp);
			NumTFs++;
			if(NumTFs>MAXTFS)
			{
				printf("TOO MANY REGULATORS DETECTED. CODE UPDATE NECESSARY TO AVOID SEG-FAULT.\n");
				exit(1);
			}

		}

		if(cnt2>=0)
		{
			NumInteractions++;
			if(reg[cnt].M[cnt2]==0) NumUniqueInteractions++;
			reg[cnt].M[cnt2]=reg[cnt].M[cnt2]+P;
			reg[cnt].W[cnt2]=reg[cnt].W[cnt2]+P;
		}
	}
	fclose(fid);
	return 0;
}

int NormalizePriorData(REGULATION *reg, GENES *gen)
{
	int cnt, cnt1, c;
	double A, B;
	double PriorMean, PriorStd;

	int LocConditions=NumConditions;
	if(LeaveOutSample) LocConditions--;
	else if(JackKnife) LocConditions=JackKnife;
	fprintf(stderr, "\nNetwork Data Stats:\n");
//	fprintf(stderr, "NumRegulators:%u, NumGenes:%u, NumConditions:%u (%u used for Network), NumRegulatoryInteractions:%u (%u unique)\n", NumTFs, NumGenes, NumConditions, LocConditions, NumInteractions, NumUniqueInteractions);
	fprintf(stderr, "NumReg:%u, NumTFs:%u, NumMiRs:%u, NumGenes:%u, NumConditions:%u (%u used for Network), NumRegulatoryInteractions:%u (%u unique)\n", NumTFs, NumTFsm, NumMiRs, NumGenes, NumConditions, LocConditions, NumInteractions, NumUniqueInteractions); // *OBS!

	for(cnt1=0; cnt1<NumTFs; cnt1++)
	{
		reg[cnt1].exp=0;
		reg[cnt1].stdev=0;
		for(cnt=0; cnt<NumGenes; cnt++)
		{
			reg[cnt1].exp+=reg[cnt1].W[cnt];
			reg[cnt1].stdev+=reg[cnt1].W[cnt]*reg[cnt1].W[cnt];
		}
		reg[cnt1].exp=reg[cnt1].exp/NumGenes;
		reg[cnt1].stdev=sqrt(reg[cnt1].stdev/NumGenes-reg[cnt1].exp*reg[cnt1].exp);
	}
	
	A=0;
	B=0;
	c=0;
	for(cnt=0; cnt<NumGenes; cnt++)
	{
		gen[cnt].exp=0;
		gen[cnt].stdev=0;
		for(cnt1=0; cnt1<NumTFs; cnt1++)
		{
			A+=reg[cnt1].W[cnt];
			B+=reg[cnt1].W[cnt]*reg[cnt1].W[cnt];
			c++;
			gen[cnt].exp+=reg[cnt1].W[cnt];
			gen[cnt].stdev+=reg[cnt1].W[cnt]*reg[cnt1].W[cnt];
		}
		gen[cnt].exp=gen[cnt].exp/NumTFs;
		gen[cnt].stdev=sqrt(gen[cnt].stdev/NumTFs-gen[cnt].exp*gen[cnt].exp);
	}
	A=A/c;
	B=sqrt(B/c-A*A);
	PriorMean=A;
	PriorStd=B;
	
	for(cnt1=0; cnt1<NumTFs; cnt1++)
	{
		for(cnt=0; cnt<NumGenes; cnt++)
		{
			if(gen[cnt].stdev==0)
			{
				if(reg[cnt1].stdev==0) reg[cnt1].W[cnt]=2*(reg[cnt1].W[cnt]-PriorMean)/(sqrt(2)*PriorStd);
				else reg[cnt1].W[cnt]=(reg[cnt1].W[cnt]-reg[cnt1].exp)/(sqrt(2)*reg[cnt1].stdev)+(reg[cnt1].W[cnt]-PriorMean)/(sqrt(2)*PriorStd);

			}
			else if(reg[cnt1].stdev==0) reg[cnt1].W[cnt]=(reg[cnt1].W[cnt]-PriorMean)/(sqrt(2)*PriorStd)+(reg[cnt1].W[cnt]-gen[cnt].exp)/(sqrt(2)*gen[cnt].stdev);
			else reg[cnt1].W[cnt]=(reg[cnt1].W[cnt]-reg[cnt1].exp)/(sqrt(2)*reg[cnt1].stdev)+(reg[cnt1].W[cnt]-gen[cnt].exp)/(sqrt(2)*gen[cnt].stdev);
			// reg[cnt1].M[cnt]=reg[cnt1].W[cnt];
		}
	}
	return 0;
}

int ReadInInteractionData(REGULATION *reg, GENES *gen, char filename[])
{
	int cnt, cnt1, cnt2, c;
	double A,B;
	float P;

	if((fid=fopen(filename, "r"))==NULL)
	{
		printf("ERROR OPENING PROTEIN INTERACTION DATA FILE\n");
		exit(1);
	}

	fprintf(stderr, "Reading in Protein Interation Data!\n");

	while(fscanf(fid, "%s\t%s\t%f", temp,temp1,&P)==3)
	{
		cnt1=NumTFs;
		for(c=0; c<NumTFs; c++)
		{
			if(strcmp(temp,reg[c].name)==0)
			{
				cnt1=c;
				break;
			}
		}
		if(cnt1==NumTFs)
		{
			NumTFs++;
			if(NumTFs>MAXTFS)
			{
				printf("TOO MANY REGULATORS DETECTED. CODE UPDATE NECESSARY TO AVOID SEG-FAULT.\n");
				exit(1);
			}
			fprintf(stderr, "WARNING: Protein %s found in Interaction Data but not Regulation Data.  It will be added into Regulation matrix with no known targets.  If this is an error, exit the program and verify your input files.\n", temp);
			strcpy(reg[cnt1].name, temp);
			for(cnt=0; cnt<NumGenes; cnt++)
			{
				reg[cnt1].M[cnt]=0;
				reg[cnt1].W[cnt]=0;
				reg[cnt1].T[cnt][0]=0;
				reg[cnt1].T[cnt][1]=0;
			}
			for(cnt=0; cnt<NumTFs; cnt++) reg[cnt1].P[cnt]=0;
		}

		cnt2=NumTFs;
		for(c=0; c<NumTFs; c++)
		{
			if(strcmp(temp1,reg[c].name)==0)
			{
				cnt2=c;
				break;
			}
		}
		if(cnt2==NumTFs)
		{
			NumTFs++;
			if(NumTFs>MAXTFS)
			{
				printf("TOO MANY REGULATORS DETECTED. CODE UPDATE NECESSARY TO AVOID SEG-FAULT.\n");
				exit(1);
			}
			fprintf(stderr, "WARNING: Protein %s found in Interaction Data but not Regulation Data.  It will be added into Regulation matrix with no known targets.  If this is an error, exit the program and verify your input files.\n", temp1);
			strcpy(reg[cnt2].name, temp1);
			for(cnt=0; cnt<NumGenes; cnt++)
			{
				reg[cnt2].M[cnt]=0;
				reg[cnt2].W[cnt]=0;
				reg[cnt2].T[cnt][0]=0;
				reg[cnt2].T[cnt][1]=0;
			}
			for(cnt=0; cnt<NumTFs; cnt++) reg[cnt2].P[cnt]=0;
		}
		
		reg[cnt1].P[cnt2]=reg[cnt1].P[cnt2]+P;
		reg[cnt2].P[cnt1]=reg[cnt2].P[cnt1]+P;
	}
	fclose(fid);
	for(cnt=0; cnt<NumTFs; cnt++) reg[cnt].P[cnt]=1;

	if(verboseoutput==2)
	{
		sprintf(output_file, "%s_InitialPPI.pairs", outtag);
		fid=fopen(output_file, "w");
		for(cnt=0; cnt<NumTFs; cnt++)
		{
			for(cnt1=cnt+1; cnt1<NumTFs; cnt1++)
			{
				fprintf(fid, "%s\t%s\t%f\n", reg[cnt1].name, reg[cnt].name,reg[cnt1].P[cnt]);
			}
		}
		fclose(fid);
	}

	A=0;
	B=0;
	c=0;
	for(cnt=0; cnt<NumTFs; cnt++)
	{
		reg[cnt].exp=0;
		reg[cnt].stdev=0;
		for(cnt1=0; cnt1<NumTFs; cnt1++)
		{
			A+=reg[cnt].P[cnt1];
			B+=reg[cnt].P[cnt1]*reg[cnt].P[cnt1];
			c++;
			reg[cnt].exp+=reg[cnt].P[cnt1];
			reg[cnt].stdev+=reg[cnt].P[cnt1]*reg[cnt].P[cnt1];
		}
		reg[cnt].exp=reg[cnt].exp/NumTFs;
		reg[cnt].stdev=sqrt(reg[cnt].stdev/NumTFs-reg[cnt].exp*reg[cnt].exp);
	}
	A=A/c;
	B=sqrt(B/c-A*A);
		
	for(cnt1=0; cnt1<NumTFs; cnt1++)
	{
		for(cnt=cnt1; cnt<NumTFs; cnt++)
		{
			reg[cnt].P[cnt1]=(reg[cnt].P[cnt1]-reg[cnt].exp)/(sqrt(2)*reg[cnt].stdev)+(reg[cnt].P[cnt1]-reg[cnt1].exp)/(sqrt(2)*reg[cnt1].stdev);
			reg[cnt1].P[cnt]=reg[cnt].P[cnt1];
		}
	}
	return 0;
}

int Correlation(REGULATION *reg, GENES *gen)
{
	int cnt, cnt1, cnt2, c, v;
	double A, B, C, F;
	float P;
	float covariateweight[MAXCONDITIONS];

	// Intitalize Weights Vector
	for(cnt=0; cnt<NumConditions; cnt++) covariateweight[cnt]=1;
	if(LeaveOutSample>0)
	{
		if(LeaveOutSample>NumConditions)
		{
			fprintf(stderr, "WARNING: Sample chosen to exclude does not exist. Will use all samples to build network.\n");
			LeaveOutSample=0;
		}
		else covariateweight[LeaveOutSample-1]=0;
	}
	if(JackKnife>=NumConditions)
	{
		fprintf(stderr, "WARNING: Number of Samples Chosen to Keep in Jack-knife is greater than or equal to the total number of Samples. No samples will be removed.\n");
		JackKnife=NumConditions;
	}
	else if(JackKnife>0)
	{
		if(JackKnife<3)
		{
			fprintf(stderr, "WARNING: Number of Samples Chosen to Keep in Jack-knife is too few to calculate a correlation. Defaulting to using all samples.\n");
			JackKnife=NumConditions;
		}
		else
		{
			for(cnt=JackKnife; cnt<NumConditions; cnt++) covariateweight[cnt]=0;
			randweights=1;
			fprintf(stderr, "Choosing %u Random Samples!\n", JackKnife);
			int RandConditions[MAXCONDITIONS];
			float tempweight[MAXCONDITIONS];
			srand(randomseed);
			for(c=0; c<NumConditions; ++c)
			{
				v= rand() % (c+1);
				RandConditions[c]=RandConditions[v];
				RandConditions[v]=c;
			}
			for(c=0; c<NumConditions; c++) tempweight[c]=covariateweight[RandConditions[c]];
			for(c=0; c<NumConditions; c++) covariateweight[c]=tempweight[c];
		}
	}

	if(weightedpearson==1)
	{
		if((fid=fopen(covariate_file, "r"))==NULL)
		{
			printf("ERROR OPENING COVARIARE WEIGHT FILE, DEFAULTING TO UNWEIGHTED PEARSON\n");
		}

		cnt=0;
		while(fscanf(fid, "%f", &P)==1)
		{
			covariateweight[cnt]=P;
			cnt++;
		}
		fclose(fid);
		
		if(randweights==1)
		{
			fprintf(stderr, "Randomizing Condition Labels!\n");
			int RandConditions[MAXCONDITIONS];
			float tempweight[MAXCONDITIONS];
			srand(randomseed);
			for(c=0; c<NumConditions; ++c)
			{
				v= rand() % (c+1);
				RandConditions[c]=RandConditions[v];
				RandConditions[v]=c;
			}
			for(c=0; c<NumConditions; c++) tempweight[c]=covariateweight[RandConditions[c]];
			for(c=0; c<NumConditions; c++) covariateweight[c]=tempweight[c];
		}
		
		if(cnt!=NumConditions) fprintf(stderr, "WARNING!!! NUMBER OF ROWS IN COVARIATE FILE DO NOT MATCH NUMBER OF NUMERIC COLUMNS IN EXPRESSION FILE. YOU MAY WANT TO KILL PROGRAM AND CHECK INPUT FILES\n");
	}
	else if(randweights==1)
	{
		fprintf(stderr, "Note that the random condition-labelling options only work if you invoke the weighted-pearson option and suppply a covariate file.\n");
	}

	if(randweights==2)
	{
		fprintf(stderr, "Generating Random Weights for Conditions!\n");
		A=0;
		srand(randomseed);
		for(c=0; c<NumConditions; c++)
		{
			F=((double) (rand() % 100000+1))/100000;
			if(covariateweight[c]>0) {covariateweight[c]=F; A+=F;}
		}
		sprintf(output_file, "%s_RandWeights.pairs", outtag);
		fid=fopen(output_file, "w");
		for(c=0; c<NumConditions; c++)
		{
			// covariateweight[c]=covariateweight[c]/A;
			fprintf(fid, "%g\n", covariateweight[c]);
		}
		fclose(fid);
	}

	for(cnt2=0; cnt2<NumGenes; cnt2++)
	{
		gen[cnt2].exp=0;
		gen[cnt2].stdev=0;
		c=0;
		C=0;
		for(cnt=0; cnt<NumConditions; cnt++)
		{
			if(gen[cnt2].expression[cnt]==gen[cnt2].expression[cnt]) // make sure expression data exists (not nan)
			{
				gen[cnt2].exp+=gen[cnt2].expression[cnt]*covariateweight[cnt];
				gen[cnt2].stdev+=covariateweight[cnt]*pow(gen[cnt2].expression[cnt],2);
				C+=covariateweight[cnt];
				c++;
			}
		}
		gen[cnt2].exp=gen[cnt2].exp/C;
		gen[cnt2].stdev=gen[cnt2].stdev/C-pow(gen[cnt2].exp,2);
		if(gen[cnt2].stdev==0) fprintf(stderr, "No variation of expression for %s, will give default correlation of zero to all its interactions.\n", gen[cnt2].name);
	}

	for(cnt1=0; cnt1<NumGenes; cnt1++)
	{
		for(cnt2=cnt1; cnt2<NumGenes; cnt2++)
		{
			if(cnt2==cnt1) gen[cnt1].corr[cnt1]=1;
			else
			{
				A=0;
				B=0;
				C=0;
				c=0;
				for(cnt=0; cnt<NumConditions; cnt++)
				{
					if((gen[cnt2].expression[cnt]==gen[cnt2].expression[cnt] && gen[cnt1].expression[cnt]==gen[cnt1].expression[cnt]))
					{
						A+=covariateweight[cnt]*(gen[cnt1].expression[cnt]-gen[cnt1].exp)*(gen[cnt2].expression[cnt]-gen[cnt2].exp);
						B+=covariateweight[cnt]*pow(gen[cnt1].expression[cnt]-gen[cnt1].exp,2);
						C+=covariateweight[cnt]*pow(gen[cnt2].expression[cnt]-gen[cnt2].exp,2);
						c++;
					}
				}
				if(c>2 && gen[cnt1].stdev>0 && gen[cnt2].stdev>0) gen[cnt1].corr[cnt2]=A/(sqrt(B)*sqrt(C));
				else gen[cnt1].corr[cnt2]=0;
				gen[cnt2].corr[cnt1]=gen[cnt1].corr[cnt2];
			}
		}
	}

	if(verboseoutput==2)
	{
		sprintf(output_file, "%s_InitialCoReg.pairs", outtag);
		fid=fopen(output_file, "w");
		for(cnt=0; cnt<NumGenes; cnt++)
		{
			for(cnt1=cnt+1; cnt1<NumGenes; cnt1++)
			{
				fprintf(fid, "%s\t%s\t%f\n", gen[cnt].name, gen[cnt1].name,gen[cnt1].corr[cnt]);
			}
		}
		fclose(fid);
	}

	// Recast Co-regulation in Z-score space
	A=0;
	B=0;
	cnt=0;
	for(cnt1=0; cnt1<NumGenes; cnt1++)
	{
		gen[cnt1].exp=0;
		gen[cnt1].stdev=0;
		for(cnt2=0; cnt2<NumGenes; cnt2++)
		{
			A+=gen[cnt1].corr[cnt2];
			B+=gen[cnt1].corr[cnt2]*gen[cnt1].corr[cnt2];
			cnt++;
			gen[cnt1].exp+=gen[cnt1].corr[cnt2];
			gen[cnt1].stdev+=gen[cnt1].corr[cnt2]*gen[cnt1].corr[cnt2];
		}
		gen[cnt1].exp=gen[cnt1].exp/NumGenes;
		gen[cnt1].stdev=sqrt(gen[cnt1].stdev/NumGenes-gen[cnt1].exp*gen[cnt1].exp);
	}
	A=A/cnt;
	B=sqrt(B/cnt-A*A);

	for(cnt1=0; cnt1<NumGenes; cnt1++)
	{
		for(cnt2=cnt1; cnt2<NumGenes; cnt2++)
		{
			gen[cnt1].corr[cnt2]=(gen[cnt1].corr[cnt2]-gen[cnt1].exp)/(sqrt(2)*gen[cnt1].stdev)+(gen[cnt1].corr[cnt2]-gen[cnt2].exp)/(sqrt(2)*gen[cnt2].stdev);
			gen[cnt2].corr[cnt1]=gen[cnt1].corr[cnt2];
		}
	}
	return 0;
}

int IdentityPPI(REGULATION *reg, GENES *gen)
{
	int cnt, cnt1;

	for(cnt=0; cnt<NumTFs; cnt++) reg[cnt].P[cnt]=1;

	if(verboseoutput==2)
	{
		sprintf(output_file, "%s_InitialPPI.pairs", outtag);
		fid=fopen(output_file, "w");
		for(cnt=0; cnt<NumTFs; cnt++)
		{
			for(cnt1=cnt+1; cnt1<NumTFs; cnt1++)
			{
				fprintf(fid, "%s\t%s\t%f\n", reg[cnt1].name, reg[cnt].name,reg[cnt1].P[cnt]);
			}
		}
		fclose(fid);
	}

	for(cnt=0; cnt<NumTFs; cnt++)
	{
		reg[cnt].exp=0;
		reg[cnt].stdev=0;
		for(cnt1=0; cnt1<NumTFs; cnt1++)
		{
			reg[cnt].exp+=reg[cnt].P[cnt1];
			reg[cnt].stdev+=reg[cnt].P[cnt1]*reg[cnt].P[cnt1];
		}
		reg[cnt].exp=reg[cnt].exp/NumTFs;
		reg[cnt].stdev=sqrt(reg[cnt].stdev/NumTFs-reg[cnt].exp*reg[cnt].exp);
	}
	
	for(cnt1=0; cnt1<NumTFs; cnt1++)
	{
		for(cnt=0; cnt<NumTFs; cnt++)
		{
			reg[cnt].P[cnt1]=(reg[cnt].P[cnt1]-reg[cnt].exp)/(sqrt(2)*reg[cnt].stdev)+(reg[cnt].P[cnt1]-reg[cnt1].exp)/(sqrt(2)*reg[cnt1].stdev);
		}
	}
	return 0;
}
int IdentityCorrelation(GENES *gen)
{
	int cnt, cnt1, cnt2;
	double A, B;

	for(cnt=0; cnt<NumGenes; cnt++) gen[cnt].corr[cnt]=1;

	if(verboseoutput==2)
	{
		sprintf(output_file, "%s_InitialCoReg.pairs", outtag);
		fid=fopen(output_file, "w");
		for(cnt=0; cnt<NumGenes; cnt++)
		{
			for(cnt1=cnt+1; cnt1<NumGenes; cnt1++)
			{
				fprintf(fid, "%s\t%s\t%f\n", gen[cnt].name, gen[cnt1].name,gen[cnt1].corr[cnt]);
			}
		}
		fclose(fid);
	}

	A=0;
	B=0;
	cnt=0;
	for(cnt1=0; cnt1<NumGenes; cnt1++)
	{
		gen[cnt1].exp=0;
		gen[cnt1].stdev=0;
		for(cnt2=0; cnt2<NumGenes; cnt2++)
		{
			A+=gen[cnt1].corr[cnt2];
			B+=gen[cnt1].corr[cnt2]*gen[cnt1].corr[cnt2];
			cnt++;
			gen[cnt1].exp+=gen[cnt1].corr[cnt2];
			gen[cnt1].stdev+=gen[cnt1].corr[cnt2]*gen[cnt1].corr[cnt2];
		}
		gen[cnt1].exp=gen[cnt1].exp/NumGenes;
		gen[cnt1].stdev=sqrt(gen[cnt1].stdev/NumGenes-gen[cnt1].exp*gen[cnt1].exp);
	}
	A=A/cnt;
	B=sqrt(B/cnt-A*A);

	for(cnt1=0; cnt1<NumGenes; cnt1++)
	{
		for(cnt2=0; cnt2<NumGenes; cnt2++)
		{
			gen[cnt1].corr[cnt2]=(gen[cnt1].corr[cnt2]-gen[cnt1].exp)/(sqrt(2)*gen[cnt1].stdev)+(gen[cnt1].corr[cnt2]-gen[cnt2].exp)/(sqrt(2)*gen[cnt2].stdev);
		}
	}
	return 0;
}


int PrintStats(REGULATION *reg, GENES *gen, char filename[])
{
	int cnt1, cnt2;

	statfid=fopen(filename, "w");
	for(cnt1=0; cnt1<NumTFs; cnt1++)
	{
		for(cnt2=0; cnt2<NumGenes; cnt2++)
		{
			fprintf(statfid,"%s\t%s\t%f\t%f\n",reg[cnt1].name,gen[cnt2].name,reg[cnt1].M[cnt2],reg[cnt1].W[cnt2]);
		}
	}
	fclose(statfid);
	return 0;
}

int PrintPPI(REGULATION *reg, GENES *gen, char filename[])
{
	int cnt, cnt1;
	statfid=fopen(filename, "w");
	for(cnt=0; cnt<NumTFs; cnt++)
	{
		for(cnt1=cnt+1; cnt1<NumTFs; cnt1++)
		{
			fprintf(statfid, "%s\t%s\t%f\n", reg[cnt1].name, reg[cnt].name,reg[cnt1].P[cnt]);
		}
	}
	fclose(statfid);
	return 0;
}

int PrintCoReg(REGULATION *reg, GENES *gen, char filename[])
{
	int cnt,cnt1;

	statfid=fopen(filename, "w");
	for(cnt=0; cnt<NumGenes; cnt++)
	{
		for(cnt1=cnt+1; cnt1<NumGenes; cnt1++)
		{
			fprintf(statfid, "%s\t%s\t%f\n", gen[cnt].name, gen[cnt1].name,gen[cnt1].corr[cnt]);
		}
	}
	fclose(statfid);
	return 0;
}

double CDF(double value)
{
	// CDF of a Z-Score using approximation from Abramowitz and Stegun 26.2.17

	// CDF variables
	float b0;
	float b1=0.319381530;
	float b2=-0.356563782;
	float b3=1.781477937;
	float b4=-1.821255978;
	float b5=1.330274429;

	if(value>=0)
	{
		b0=1/(1+0.2316419*value);
		value=1-0.39894228*exp(-0.5*pow(value,2))*b0*(b0*(b0*(b0*(b0*b5+b4)+b3)+b2)+b1);
	}
	else
	{
		b0=1/(1-0.2316419*value);
		value=0.39894228*exp(-0.5*pow(value,2))*b0*(b0*(b0*(b0*(b0*b5+b4)+b3)+b2)+b1);
	}
	return value;
}

double inverseCDF(double value)
{
	double c1[6] = {-3.969683028665376e+01,2.209460984245205e+02,-2.759285104469687e+02,1.383577518672690e+02,-3.066479806614716e+01,2.506628277459239e+00};
	double c2[5] = {-5.447609879822406e+01,1.615858368580409e+02,-1.556989798598866e+02,6.680131188771972e+01,-1.328068155288572e+01};
	double c3[6] = {-7.784894002430293e-03,-3.223964580411365e-01,-2.400758277161838e+00,-2.549732539343734e+00,4.374664141464968e+00,2.938163982698783e+00};
	double c4[4] = {7.784695709041462e-03,3.224671290700398e-01,2.445134137142996e+00,3.754408661907416e+00};
	double p_low=0.02425;
	double p_high=1-p_low;
	double val1, val2;
	if(value<=0) value=-40;
	else if(value>=1) value=40;
	else if(value < p_low)
	{
		val1=sqrt(-2*log(value));
		value=(((((c3[0]*val1+c3[1])*val1+c3[2])*val1+c3[3])*val1+c3[4])*val1+c3[5])/((((c4[0]*val1+c4[1])*val1+c4[2])*val1+c4[3])*val1+1);
	}
	else if(value > p_high)
	{
		val1=sqrt(-2*log(1-value));
		value=-1*(((((c3[0]*val1+c3[1])*val1+c3[2])*val1+c3[3])*val1+c3[4])*val1+c3[5])/((((c4[0]*val1+c4[1])*val1+c4[2])*val1+c4[3])*val1+1);
	}
	else
	{
		val2=value-0.5;
		val1=val2*val2;
		value=(((((c1[0]*val1+c1[1])*val1+c1[2])*val1+c1[3])*val1+c1[4])*val1+c1[5])*val2/(((((c2[0]*val1+c2[1])*val1+c2[2])*val1+c2[3])*val1+c2[4])*val1+1);
	}
	return value;
}

int RandPerm(GENES *gen)
{
	int c,v;

	fprintf(stderr, "Randomizing Gene Labels!\n");
	int RandSignature[MAXGENES];
	srand(randomseed);
	
	for(c=0; c<NumGenes; ++c)
	{
		v= rand() % (c+1);
		RandSignature[c]=RandSignature[v];
		RandSignature[v]=c;
	}

	for(c=0; c<NumGenes; c++) strcpy(gen[c].temp, gen[RandSignature[c]].name);
	for(c=0; c<NumGenes; c++) strcpy(gen[c].name, gen[c].temp);
	return 0;
}

void useage()
{
	printf ("%s Useage %s\n", BOLD, pname);
	printf ("%s\t-e (required) file of expression values (can alternately be a list of gene names)%s\n", BOLD, NORMAL);
	printf ("%s\t-m (required) pair file of motif edges%s\n", BOLD, NORMAL);
	printf ("%s\t-p (optional) pair file of PPI edges%s\n", BOLD, NORMAL);
	printf ("%s\t-o (optional) tag for output files%s\n", BOLD, NORMAL);
	printf ("%s\t-a (optional) value to be used for update variable, alpha (default=0.1)%s\n", BOLD, NORMAL);
	printf ("%s Additional options (see README): %s\n", BOLD, NORMAL);
	printf ("%s\t-k (optional) kill the program after it has run k steps (default=1000)%s\n", BOLD, NORMAL);
	printf ("%s\t-n (optional) output a \"stats\" file every n steps (default, no stats file)%s\n", BOLD, NORMAL);
	printf ("%s\t-w (optional) file with list of covariate weights%s\n", BOLD, NORMAL);
	printf ("%s\t-l (optional) leave out the lth sample when building the network%s\n", BOLD, NORMAL);
	printf ("%s\t-j (optional) retain only j samples when building the network%s\n", BOLD, NORMAL);
	printf ("%s\t-r (optional) randomization options%s\n", BOLD, NORMAL);
	printf ("%s\t-s (optional) value to seed the random number generator (defaults to system time)%s\n", BOLD, NORMAL);
	printf ("%s\t-v (optional) verbose output options%s\n", BOLD, NORMAL);
}

void SignalHandler(int signum)
{
	sprintf(output_file, "%s_Step%u.stats", outtag, killstep);
	PrintStats((REGULATION *) & Regulation, (GENES *) & Genes, output_file);
	fprintf(stderr, "Caught signal %d.  Current values printed to %s\n", signum, output_file);
	exit(signum);
}
