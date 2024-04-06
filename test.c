#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include"chromosome.h"

void main(){
	srand(time(0));

	int chromosome1[]={2,3,6,1,8,9,5,4,7,0};
	int chromosome2[]={7,1,3,9,5,6,8,2,0,4};
	int numGenes=10;

	printf("Before Crossover\n");

	for(int i=0;i<numGenes;i++)
		printf("%d ",chromosome1[i]);
	printf("\n");
	
	for(int i=0;i<numGenes;i++)
		printf("%d ",chromosome2[i]);
	printf("\n");

	crossover(chromosome1,chromosome2,numGenes);

	printf("After Crossover\n");

	for(int i=0;i<numGenes;i++)
		printf("%d ",chromosome1[i]);
	printf("\n");

	for(int i=0;i<numGenes;i++)
		printf("%d ",chromosome2[i]);
	printf("\n");
}

//	int indices[] = {0,1,2,3,4,5}
//	Chromosome chormosomes[] = {{3},{7},{2},{6},{10},{9}}
//	sort indices based on chromosomes[indices[i]].fitness
//	{2,0,3,1,5,4}
//	for i:=0 to n do acess chromosome[indices[i]]
