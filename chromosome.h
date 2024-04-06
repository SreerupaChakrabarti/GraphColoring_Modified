#ifndef __CHROMOSOME__
    #define __CHROMOSOME__

    #include <stdlib.h>
    #include <stdio.h>

        //Structure for each chromosome
        typedef struct Chromosome {
            int seqLength;  // Length of the sequence
            int *sequence;  // An array of length seqLength. Defines a coloration of the graph
            int numConflicts;   // number of conflicting edges for the given coloring
            int fitness;    // fitness based on the number of conflicting edges
            int rank;       // rank of the chromosome
            double rankProbability; // rank probability of the chromosome
        } Chromosome;

        void swap(int *num1, int *num2) {
            int t = *num1;
            *num1 = *num2;
            *num2 = t;
        }

        /*
            copyChromosome():
                Input: The source and destination chromosomes.
                Output: A deep copy of the source chromosome to the destination chromosome
        */
        void copyChromosome(Chromosome *srcChromosome, Chromosome *destChromosome) {
            destChromosome->numConflicts = srcChromosome->numConflicts;
            destChromosome->fitness = srcChromosome->fitness;

            destChromosome->seqLength = srcChromosome->seqLength;

            destChromosome->sequence = (int *)calloc(destChromosome->seqLength, sizeof(int));
            for (int i = 0; i < srcChromosome->seqLength; i++) {
                destChromosome->sequence[i] = srcChromosome->sequence[i];
            }
        }

        // Initiates each chromosome in the given array
        void getRandomChromosomes(Chromosome chromosomes[], int numChromosomes, int numVertices, int highestColor) {
            for (int i = 0; i < numChromosomes; i++) {
                chromosomes[i].seqLength = numVertices;
                chromosomes[i].sequence = (int *)calloc(numVertices, sizeof(int));

                for (int j = 0; j < numVertices; j++) {
                    chromosomes[i].sequence[j] = rand() % highestColor + 1;
                }
            }
        }

        void displayChromosomes(Chromosome chromosomes[], int numChromosomes) {
            printf("Chromosome\tConflicts\tFitness\n");
            for (int i = 0; i < numChromosomes; i++) {
                for (int j = 0; j < chromosomes[i].seqLength; j++) {
                    printf("%d ", chromosomes[i].sequence[j]);
                }

                printf("\t%d\t%d\n", chromosomes[i].numConflicts, chromosomes[i].fitness);
            }
        }

    /*
        selectChromosomes():
            Input: The fitnesses of chromosomes and the number of chromosomes
            Output: A selection/ mating pool of chromosomes based on their fitness values
        [Selection is done using Rank Selection method]
    */
    // Rank selection function
    void selectChromosomes(Chromosome chromosomes[], Chromosome matingPool[], int numChromosomes) {

        for (int i = 0; i < numChromosomes - 1; i++) { //chromosmes are sorted based on their fitness value
            for (int j = 0; j < numChromosomes - i - 1; j++) {
                if (chromosomes[j].fitness > chromosomes[j + 1].fitness) {
                    Chromosome temp = chromosomes[j];
                    chromosomes[j] = chromosomes[j + 1];
                    chromosomes[j + 1] = temp;
                }
            }
        }

        //assigning the ranks
        for (int i = 0; i < numChromosomes; i++) {
            chromosomes[i].rank = numChromosomes - i;
        }

        //assigning the ranks
        for (int i = 0; i < numChromosomes; i++) {
            chromosomes[i].rank = numChromosomes - i;
        }

        double totalRank = (numChromosomes * (numChromosomes + 1)) / 2.0;  //calculating the sum of ranks

        for (int i = 0; i < numChromosomes; i++) {
            chromosomes[i].rankProbability = (double)chromosomes[i].rank / totalRank;  //calculatng rank probability
        }

        double cumulativeProb = 0.0;

        for (int i = 0; i < numChromosomes; i++) {
            cumulativeProb = cumulativeProb + chromosomes[i].rankProbability;

            double randNum = (double)rand() / RAND_MAX;

            //select chromosomes if random number is less than or equal to cumulative probability
            if (randNum <= cumulativeProb) {
                copyChromosome(&chromosomes[i], &matingPool[i]);
            }
        }
        return;
    }


    /*
            crossover():
                Input: Two chromosomes and their length i.e., number of genes
                Output: A two point crossover between the given chromosomes to produce two childrens
        */
        void crossover(Chromosome chromosome1,Chromosome chromosome2,int numGenes){
            int point1=0;
            int point2=0;

            while(point1>=point2){
                point1=rand()%numGenes;
                point2=rand()%numGenes;
            }
            
            for(int i=point1;i<=point2;i++){
                swap(&chromosome1.sequence[i],&chromosome2.sequence[i]);
            }
            
            return ;
        }

        /*
            crossChromosomes():
                Input: A list of chromosomes, the length of the list i.e., the number of chromosomes and crossover probability
                Output: Crossover between all the chromosomes based on the chrossover probability
        */
        void crossChromosomes(Chromosome chromosomes[],int numChromosomes,double probability){
            int index1,index2;
            int numCrossover=0;

            for(int i=0;i<numChromosomes/2;i++){
                double random=1.0*rand()/RAND_MAX;
                
                //printf("Rand number: %lf\n",random);
                
                index1=0;
                index2=0;

                if(random<=probability){
                    while(index1==index2){
                        index1=rand()%numChromosomes;
                        index2=rand()%numChromosomes;
                    }

                    crossover(chromosomes[index1],chromosomes[index2],chromosomes[index1].seqLength);
                    numCrossover++;

                }
            }

            return ;
        }

    // /*
    //     crossover():
    //         Input: Two chromosomes and their length i.e., number of genes
    //         Output: Order based crossover between the given chromosomes to produce two children
    // */
    // void crossover(Chromosome chromosome1, Chromosome chromosome2, Chromosome *child1, Chromosome *child2) {
    //     int *visited1 = (int *)calloc(chromosome1.seqLength, sizeof(int));
    //     int *visited2 = (int *)calloc(chromosome2.seqLength, sizeof(int));

    //     int point1 = rand() % chromosome1.seqLength;
    //     int point2 = rand() % chromosome2.seqLength;
    //     int temp = point2;

    //     // Mark visited genes
    //     for (int i = point1; i != temp; i = (i + 1) % chromosome1.seqLength) {
    //         child1->sequence[i] = chromosome1.sequence[i];
    //         visited1[child1->sequence[i]] = 1;
    //     }
    //     child1->sequence[temp] = chromosome2.sequence[temp];
    //     visited1[child1->sequence[temp]] = 1;

    //     for (int i = point2; i != point1; i = (i + 1) % chromosome2.seqLength) {
    //         child2->sequence[i] = chromosome2.sequence[i];
    //         visited2[child2->sequence[i]] = 1;
    //     }
    //     child2->sequence[point1] = chromosome1.sequence[point1];
    //     visited2[child2->sequence[point1]] = 1;

    //     // Fill remaining positions in children
    //     int index1 = (point1 + 1) % chromosome1.seqLength;
    //     int index2 = (point2 + 1) % chromosome2.seqLength;

    //     for (int i = (point2 + 1) % chromosome2.seqLength; i != point2; i = (i + 1) % chromosome2.seqLength) {
    //         if (!visited1[chromosome2.sequence[i]]) {
    //             child1->sequence[index1] = chromosome2.sequence[i];
    //             index1 = (index1 + 1) % chromosome1.seqLength;
    //             visited1[child1->sequence[index1]] = 1;
    //         }
    //     }

    //     for (int i = (point1 + 1) % chromosome1.seqLength; i != point1; i = (i + 1) % chromosome1.seqLength) {
    //         if (!visited2[chromosome1.sequence[i]]) {
    //             child2->sequence[index2] = chromosome1.sequence[i];
    //             index2 = (index2 + 1) % chromosome2.seqLength;
    //             visited2[child2->sequence[index2]] = 1;
    //         }
    //     }

    //     free(visited1);
    //     free(visited2);
    // }

    // /*
    //     crossChromosomes():
    //         Input: A list of chromosomes, the length of that list i.e., the number of chromosomes, and crossover probability
    //         Output: Order based crossover between all the chromosomes based on the crossover probability
    // */
    // void crossChromosomes(Chromosome chromosomes[], int numChromosomes, double probability) {
    //     for (int i = 0; i < numChromosomes / 2; i++) {
    //         double random = (double)rand() / RAND_MAX;

    //         if (random <= probability) {
    //             int index1 = rand() % numChromosomes;
    //             int index2 = rand() % numChromosomes;

    //             Chromosome child1, child2;
    //             crossover(chromosomes[index1], chromosomes[index2], &child1, &child2);
    //             copyChromosome(&child1, &chromosomes[index1]);
    //             copyChromosome(&child2, &chromosomes[index2]);
    //         }
    //     }
    // }

    /*
        mutate():
            Input: A chromosome, length of that chromosome, and the highest order gene
            Output: Reverses the sequence between two randomly selected points to generate a new mutated sequence
    */

   void mutate(Chromosome chromosome,int numGenes,int highestGene){
		int randIndex=rand()%numGenes;
		
		int randGene=chromosome.sequence[randIndex];
		while(randGene==chromosome.sequence[randIndex]){
			randGene=rand()%highestGene+1;
		}
		
		chromosome.sequence[randIndex]=randGene;

		return ;
	}
	
	/*
		mutateChromosomes():
			Input: A list of chromosomes, length of that list, mutation probability, known chromatic number
			Output: Mutation of the chromosomes based on the mutation probability
	*/
	void mutateChromosomes(Chromosome chromosomes[],int numChromosomes,double probability,int chromaticNum){
		int chromosome;
		int numMutation=0;

		for(int i=0;i<numChromosomes;i++){
			double random=1.0*rand()/RAND_MAX;
			
			chromosome=0;

			if(random<=probability){
				chromosome=rand()%numChromosomes;

				mutate(chromosomes[chromosome],chromosomes[chromosome].seqLength,chromaticNum);
				numMutation++;
			}
		}

		return ;
	}
    // void mutate(Chromosome *chromosome, int numGenes, int highestGene) {
    //     int point1=rand() % numGenes;
    //     int point2=rand() % numGenes;  //generates two random points in the sequence chromosome

    //     int start=point1 < point2 ? point1 : point2;
    //     int end=point1 < point2 ? point2 : point1;  //sets the start and end point for the swapping to take place

    //     while (start<end) {
    //         swap(&chromosome->sequence[start], &chromosome->sequence[end]);
    //         start++;
    //         end--;
    //     }
    // }

    // /*
    //     mutateChromosomes():
    //         Input: A list of chromosomes, length of that list, mutation probability, and known chromatic number
    //         Output: Mutation of the chromosomes based on the mutation probability
    // */

    // void mutateChromosomes(Chromosome chromosomes[], int numChromosomes, double probability, int chromaticNum) {
    //     int numMutation=0;
    //     for (int i= 0; i< numChromosomes; i++) {
    //         double random=1.0*rand()/RAND_MAX; 

    //         if (random<=probability) {
    //             mutate(&chromosomes[i], chromosomes[i].seqLength, chromaticNum); //function call to mutate the sequence
    //             numMutation++;
    //         }
    //     }
    // }

#endif 

