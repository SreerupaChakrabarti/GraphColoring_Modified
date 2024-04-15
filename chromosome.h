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

            int indices[numChromosomes];  
            
            for (int i = 0; i < numChromosomes; i++) {
                indices[i] = i;
            }

            for (int i = 0; i < numChromosomes - 1; i++) {
                for (int j = 0; j < numChromosomes - i - 1; j++) {
                    if (chromosomes[indices[j]].fitness > chromosomes[indices[j + 1]].fitness) {
                        int temp = indices[j];
                        indices[j] = indices[j + 1];
                        indices[j + 1] = temp;
                    }
                }
            }

            // for (int i = 0; i < numChromosomes; i++) {
            //     printf("%d\n",chromosomes[indices[i]].fitness);
            // }

            //assigning the ranks
            for (int i = 0; i < numChromosomes; i++) {
                chromosomes[indices[i]].rank = numChromosomes - i;
            }

            double totalRank = (numChromosomes * (numChromosomes + 1)) / 2.0;  //calculating the sum of ranks

            for (int i = 0; i < numChromosomes; i++) {
                chromosomes[indices[i]].rankProbability = (double)chromosomes[i].rank / totalRank;  //calculatng rank probability
            }


            double cumulativeProb = 0.0;

            for (int i = 0; i < numChromosomes; i++) {
                cumulativeProb = cumulativeProb + chromosomes[indices[i]].rankProbability;

                double randNum = (double)rand() / RAND_MAX;

                //select chromosomes if random number is less than or equal to cumulative probability
                if (randNum <= cumulativeProb) {
                    copyChromosome(&chromosomes[indices[i]], &matingPool[i]);
                }
            }
            return;
        }

	 /*
            crossover():
                Input: Two chromosomes and their length i.e., number of genes
                Output: Crossover based on ordered crossover
        */
/*	void crossover(int sequence1[],int sequence2[],int numGenes){
		//Choose the random point
		int point1 = rand()%numGenes;
		
		//Choice variable determines which half to swap
		int choice = rand()%2;

		int startIndex = (choice==0)?0:point1;
		int endIndex = (choice==0)?point1:(numGenes-1);

		for(int i=startIndex;i<=endIndex;i++){
			swap(&sequence1[i],&sequence2[i]);
		}
	}
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

/*	void crossChromosomes(Chromosome chromosomes[], int numChromosomes, double probability) {
            
            int numCrossover = 0;
            for (int i = 0; i < numChromosomes - 1; i++) {
                for (int j = i + 1; j < numChromosomes; j++) {
                    double random = 1.0*rand()/RAND_MAX;
                    
                    if (random <= probability) {
                        crossover(chromosomes[i].sequence, chromosomes[j].sequence, chromosomes[i].seqLength);
                        numCrossover++;
                    }
                }
            }
	    return;
        }
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

				crossover(chromosomes[index1].sequence,chromosomes[index2].sequence,chromosomes[index1].seqLength);
				numCrossover++;

			}
		}

		return ;
	}

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

