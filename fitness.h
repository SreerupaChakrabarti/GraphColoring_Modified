#ifndef __FITNESS__
	#define __FITNESS__

	int findConflict(int edges[][2],int numEdges,int color[],int colorLength){
		int conflict=0;

		for(int i=0;i<numEdges;i++){
			if(color[edges[i][0]]==color[edges[i][1]]){
				conflict++;
			}
		}

		return conflict;
	}

	void findConflicts(int edges[][2],int numEdges,Chromosome chromosomes[],int numChromosomes){
		for(int i=0;i<numChromosomes;i++){
			chromosomes[i].numConflicts=findConflict(edges,numEdges,chromosomes[i].sequence,chromosomes[i].seqLength);
		}
	}

	void findFitnesses(Chromosome chromosomes[],int numChromosomes){
		int maxConflicts=chromosomes[0].numConflicts;
		
		for(int i=1;i<numChromosomes;i++){
			if(chromosomes[i].numConflicts > maxConflicts)
				maxConflicts = chromosomes[i].numConflicts;
		}
		
		//printf("%d\n",maxFitness);
		
		for(int i=0;i<numChromosomes;i++){
			chromosomes[i].fitness = maxConflicts-chromosomes[i].numConflicts;
		}
		
		return ;
	}

	void findConflictsAndFitnesses(int edges[][2],int numEdges,Chromosome chromosomes[],int numChromosomes){
		findConflicts(edges,numEdges,chromosomes,numChromosomes);
		findFitnesses(chromosomes,numChromosomes);
		return ;
	}
#endif
