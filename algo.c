#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<time.h>
#include<string.h>
#include"structAdj.h"
#include"chromosome.h"
#include"fitness.h"

#define NUM_CHROMOSOMES 1000
#define CROSS_PROBABILITY 0.7
#define MUTATE_PROBABILITY 0.06
#define NUM_GENERATIONS 10000

clock_t start,end;

int getChromaticNumChromosome(Chromosome chromosome){
	if(chromosome.numConflicts==0){
		int max=0;

		for(int i=0;i<chromosome.seqLength;i++){
			if(max<chromosome.sequence[i])
				max=chromosome.sequence[i];
		}

		return max;
	}
	return -1;
}

void displayHavingConflict(Chromosome chromosomes[],int numChromosomes,int targetConflict){
	for(int i=0;i<numChromosomes;i++){
		if(chromosomes[i].numConflicts==targetConflict){
			for(int j=0;j<chromosomes[i].seqLength;j++){
				printf("%d ",chromosomes[i].sequence[j]);
			}

			printf("\n");
		}
	}
}

/*
	Driver Code:
		Give input from the command line. 
		After the object file name write the graph's name with extension.

		Example:
			Compile: gcc -o opt selection.c findFitness.c getAdjList.c -lm
			Run: ./opt myciel3.col
*/

//Dummy driver code
int main(int argc,char *argv[]){
	srand((unsigned)time(NULL));

	if(argc<2){	//If argument count is less than 2 then
		printf("Please Provide the file name\n");

		exit(1);
	}
	
	//If the git repo is correctly cloned the graphs must be within the same directory.
	char filePath[100]="GCP_DATASET/";
	strcat(filePath,argv[1]);	//Hence the relative file is "GCP_DATASET/"+argv[1]
	
	FILE *file;
	file=fopen(filePath,"r");	//Open the file in read mode
	
	if(file==NULL){
		printf("Can not open file\n");

		exit(1);
	}
	
	//File has successfully been opened
	start=clock();
	int knownChromaticNum,numVertices,numEdges;

	fscanf(file,"%d %d %d",&knownChromaticNum,&numVertices,&numEdges);
	int edges[numEdges][2];
	buildEdgeMatrix(file,edges,&numEdges,numVertices);
	
	//displayEdges(edges,numEdges);

	Chromosome chromosomes[NUM_CHROMOSOMES];

	getRandomChromosomes(chromosomes,NUM_CHROMOSOMES,numVertices,knownChromaticNum);

	findConflictsAndFitnesses(edges,numEdges,chromosomes,NUM_CHROMOSOMES);
	//displayChromosomes(chromosomes,NUM_CHROMOSOMES);

	int minConflictIndex=0;
	Chromosome matingPool[NUM_CHROMOSOMES];
	
	//displayChromosomes(chromosomes,NUM_CHROMOSOMES);
	//printf("Generation,Min Conflict,Max Conflict,Avg Conflict\n");

	int i=1;
	for(i=1;i<=NUM_GENERATIONS;i++){
		minConflictIndex=0;

		for(int j=0;j<NUM_CHROMOSOMES;j++){
			if(chromosomes[j].numConflicts<chromosomes[minConflictIndex].numConflicts)
				minConflictIndex=j;
		}

		if(chromosomes[minConflictIndex].numConflicts==0)
			break;

		selectChromosomes(chromosomes,matingPool,NUM_CHROMOSOMES);
		
		crossChromosomes(matingPool,NUM_CHROMOSOMES,CROSS_PROBABILITY);
		
		mutateChromosomes(matingPool,NUM_CHROMOSOMES,MUTATE_PROBABILITY,knownChromaticNum);

		findConflictsAndFitnesses(edges,numEdges,matingPool,NUM_CHROMOSOMES);

		for(int j=0;j<NUM_CHROMOSOMES;j++){
			copyChromosome(&matingPool[j],&chromosomes[j]);
		}
	}

	end=clock();
	
	printf("Chomatic_Num Iteration_Num Time_Taken\n");
	printf("%d %d %lf\n",getChromaticNumChromosome(chromosomes[minConflictIndex]),i,(double)(end-start)/CLOCKS_PER_SEC);
	//displayHavingConflict(chromosomes,NUM_CHROMOSOMES,0);

	return 0;
}
