#include <iostream>
#include <algorithm>
#include <fstream>

#include <mpi.h>
#define INF 101
#define ROOT_PROC 0

using namespace std;

void printMatrix(int *matrix, int numberOfVert) {
    for (int i = 0; i < numberOfVert; i++) {
        for (int j = 0; j < numberOfVert; j++) {
            cout << matrix[i*numberOfVert + j] << " ";
        }
        cout << endl;
    }
}

int main(int argc, char** argv) {
    srand(time(NULL));
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int numberOfVert;
    ifstream file("matrix");
    file >> numberOfVert;
    MPI_Request request;
    MPI_Status status;

    if(rank == ROOT_PROC) {
        // numberOfVert = 5;
        int *matrix;


 
        if (size >= numberOfVert + 1)
            size = numberOfVert + 1;
        matrix = (int *)malloc(sizeof(int) * (numberOfVert * numberOfVert));
        for (int i = 0; i < numberOfVert; i++) {
            for (int j = 0; j < numberOfVert; j++) {
                file >> matrix[i*numberOfVert + j];
            }
        }

        file.close();

        cout << "Root process. Number of vertex sent to other processes" << endl;
        // MPI_Bcast(&numberOfVert, 1, MPI_INT, ROOT_PROC, MPI_COMM_WORLD);

        cout << "Root process. Print first matrix:" << endl;
        printMatrix(matrix, numberOfVert);
        int workSize = size - 1;
        int *rowCounts = (int*)malloc(sizeof(int) * workSize);

        for(int i = 0; i < workSize - 1; i++) {
            rowCounts[i] = numberOfVert / workSize;
        }
        if(numberOfVert % workSize == 0) {
            rowCounts[workSize-1] = numberOfVert / workSize;
        }
        else {
            rowCounts[workSize - 1] = numberOfVert - ((numberOfVert / workSize) * (workSize - 1));
        }

        cout << "Root process. Row counts sent to other processes" << endl;
        for(int i = 1; i < size; i++) {
            MPI_Send(&rowCounts[i-1], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        cout << "Root process. Main loop start" << endl;
        for(int k = 0; k < numberOfVert; k++) {
            for(int p = 1; p < size; p++) {
                MPI_Isend(&matrix[k*numberOfVert], numberOfVert, MPI_INT, p, 0, MPI_COMM_WORLD, &request);
            }
            int num = 0;
            for(int i = 0; i < workSize; i++) { 
                for(int j = 0; j < rowCounts[i]; j++) { 
                    MPI_Isend(&matrix[num*numberOfVert], numberOfVert, MPI_INT, i+1, 0, MPI_COMM_WORLD, &request);
                    num++;
                }
            }
            num = 0; 
            for(int i = 0; i < workSize; i++) { 
                for(int j = 0; j < rowCounts[i]; j++) {
                    MPI_Recv(&matrix[num*numberOfVert], numberOfVert, MPI_INT, i+1, 0, MPI_COMM_WORLD, &status);
                    num++;
                }
            }
        }

        cout << "Root process. Print final matrix:" << endl;
        printMatrix(matrix, numberOfVert);
        free(matrix);
    }
    else {if (rank <= numberOfVert){
        file.close();
        MPI_Bcast(&numberOfVert, 1, MPI_INT, ROOT_PROC, MPI_COMM_WORLD);
        int rowCount;
        MPI_Recv(&rowCount, 1, MPI_INT, ROOT_PROC, 0, MPI_COMM_WORLD, &status);
        cout << rank << " : " << rowCount << endl;
        int *kRow = (int *)malloc(sizeof(int) * numberOfVert);
        int *rows = (int *)malloc(sizeof(int) * (rowCount * numberOfVert));
        for(int k = 0; k < numberOfVert; k++) {
            MPI_Recv(kRow, numberOfVert, MPI_INT, ROOT_PROC, 0, MPI_COMM_WORLD, &status);
            for(int i = 0; i < rowCount; i++) {
                MPI_Recv(&rows[i*numberOfVert], numberOfVert, MPI_INT, ROOT_PROC, 0, MPI_COMM_WORLD, &status);
            }
            for (int i = 0; i < rowCount; i++) {
                for (int j = 0; j < numberOfVert; j++) {
                    rows[i*numberOfVert + j] = min(
                            rows[i*numberOfVert + j],
                            rows[i*numberOfVert + k] + kRow[j]
                    );
                }
            }
            for(int i = 0; i < rowCount; i++) {
                MPI_Send(&rows[i*numberOfVert], numberOfVert, MPI_INT, ROOT_PROC, 0, MPI_COMM_WORLD);
            }
        }
        free(kRow);
        free(rows);
    }
}
    MPI_Finalize();
    return 0;
}
