#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int conflictingVertices[5]; // Global array to track conflicting vertices

typedef struct Chromosome {
    int seqLength;  // Length of the sequence
    int *sequence;  // An array of length seqLength. Defines a coloration of the graph
    int numConflicts;   // Number of conflicting edges for the given coloring
} Chromosome;

void findConflictingVertex(const Chromosome *chromosome, int numVertices, int edges[][2], int numEdges) {
    // Reset conflictingVertices array
    for (int i = 0; i < numVertices; i++) {
        conflictingVertices[i] = 0;
    }

    for (int i = 0; i < numEdges; i++) {
        int vertex1 = edges[i][0];
        int vertex2 = edges[i][1];

        if (chromosome->sequence[vertex1] == chromosome->sequence[vertex2]) {
            conflictingVertices[vertex1] = 1;
            conflictingVertices[vertex2] = 1;
        }
    }

    printf("Conflicting Vertices:\n");
    int conflictFound = 0;
    for (int i = 0; i < numVertices; i++) {
        if (conflictingVertices[i]) {
            printf("%d ", i);
            conflictFound = 1;
        }
    }
    if (!conflictFound) {
        printf("No conflicting vertices found.");
    }
    printf("\n");
}

// Mutation function for conflicting vertices
void mutateConflictingVertices(int numVertices, int edges[][2], int numEdges, const Chromosome *chromosome) {
    int numConflicting = 0;
    for (int i = 0; i < numVertices; i++) {
        if (conflictingVertices[i]) {
            numConflicting++;
        }
    }

    if (numConflicting > 0) {
        // Choose a random conflicting vertex
        int randomIndex = rand() % numConflicting;
        int count = 0;
        int conflictingVertex = -1;
        for (int i = 0; i < numVertices; i++) {
            if (conflictingVertices[i]) {
                if (count == randomIndex) {
                    conflictingVertex = i;
                    break;
                }
                count++;
            }
        }

        if (conflictingVertex != -1) {
            printf("Selected Conflicting Vertex: %d\n", conflictingVertex);

            // Find adjacent vertices and their colors
            int *adjacentColors = (int *)malloc(numVertices * sizeof(int)); // Allocate memory for adjacent colors
            if (adjacentColors == NULL) {
                printf("Memory allocation failed.\n");
                return;
            }

            int adjacentCount = 0;
            for (int i = 0; i < numEdges; i++) {
                if (edges[i][0] == conflictingVertex) {
                    adjacentColors[adjacentCount++] = chromosome->sequence[edges[i][1]];
                } else if (edges[i][1] == conflictingVertex) {
                    adjacentColors[adjacentCount++] = chromosome->sequence[edges[i][0]];
                }
            }

            // Display adjacent vertices' colors
            printf("Adjacent Colors Array:\n");
            for (int i = 0; i < adjacentCount; i++) {
                printf("%d ", adjacentColors[i]);
            }
            printf("\n");

            // Find the available color not present in adjacent colors
            int maxColorUsed = 0;
            for (int i = 0; i < chromosome->seqLength; i++) {
                if (chromosome->sequence[i] > maxColorUsed) {
                    maxColorUsed = chromosome->sequence[i];
                }
            }

            int availableColor = -1;
            for (int i = 1; i <= maxColorUsed; i++) {
                int colorFound = 0;
                for (int j = 0; j < adjacentCount; j++) {
                    if (adjacentColors[j] == i) {
                        colorFound = 1;
                        break;
                    }
                }
                if (!colorFound && chromosome->sequence[conflictingVertex] != i) {
                    availableColor = i;
                    break;
                }
            }

            if (availableColor != -1) {
                printf("Available Color for Conflicting Vertex %d: %d\n", conflictingVertex, availableColor);
            } else {
                printf("No available color found for mutation.\n");
            }

            free(adjacentColors); // Free dynamically allocated memory
        }
    } else {
        printf("No conflicting vertices found for mutation.\n");
    }
}

int main() {
    srand((unsigned)time(NULL));

    // Define the edges of the graph
    int numVertices = 5;
    int numEdges = 6;  // Updated to include all edges
    int edges[][2] = {{0, 1}, {1, 2}, {1, 3}, {2, 3}, {1, 4}, {0, 2}};  // Added all edges including conflicting ones

    // Create a dummy chromosome with a known sequence that includes conflicts
    Chromosome dummyChromosome;
    dummyChromosome.seqLength = numVertices;
    dummyChromosome.sequence = (int *)malloc(numVertices * sizeof(int));
    // Manually assign colors that create conflicts
    dummyChromosome.sequence[0] = 1;  // Vertex 0 is color 1
    dummyChromosome.sequence[1] = 3;  // Vertex 1 is color 3 
    dummyChromosome.sequence[2] = 1;  // Vertex 2 is color 1 
    dummyChromosome.sequence[3] = 3;  // Vertex 3 is color 3
    dummyChromosome.sequence[4] = 2;  // Vertex 4 is color 2 

    // Display the original chromosome
    printf("Original Chromosome:\n");
    for (int i = 0; i < numVertices; i++) {
        printf("%d ", dummyChromosome.sequence[i]);
    }
    printf("\n");

    // Find conflicting vertices in the dummy chromosome
    findConflictingVertex(&dummyChromosome, numVertices, edges, numEdges);

    // Perform mutation on conflicting vertices (outside of findConflictingVertex)
    mutateConflictingVertices(numVertices, edges, numEdges, &dummyChromosome);

    free(dummyChromosome.sequence);

    return 0;
}
