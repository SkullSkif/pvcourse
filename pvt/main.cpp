#include <iostream>
#include <algorithm>
#include <fstream>

#include <mpi.h>

//Максимальное значение веса = 100
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
    double time, timeResult;
    MPI_Request request;
    MPI_Status status;

    if(rank == ROOT_PROC) {
        numberOfVert = 5;
        int *matrix;
        time = MPI_Wtime();
        //Главный процесс считывает матрицу из файла и раздает ее остальным
        ifstream file("matrix");

        file >> numberOfVert;

        //Матрица смежности с весами ребер графа(101 - ребра нет, 0 ребро в себя)
        //Матрица будет лежать в строку, чтобы проще было передавать процессам
        matrix = (int *)malloc(sizeof(int) * (numberOfVert * numberOfVert));

        //Считываем матрицу весов ребер
        for (int i = 0; i < numberOfVert; i++) {
            for (int j = 0; j < numberOfVert; j++) {
                file >> matrix[i*numberOfVert + j];
            }
        }

        file.close();

        cout << "Root process. Number of vertex sent to other processes" << endl;
        //Раздать процессам количество вершин
        MPI_Bcast(&numberOfVert, 1, MPI_INT, ROOT_PROC, MPI_COMM_WORLD);

        cout << "Root process. Print first matrix:" << endl;
        printMatrix(matrix, numberOfVert);

        //Раздать каждому процессу нужно количество строк матрицы

        //Вычислить по сколько строк отдать каждому процессу
        int workSize = size - 1;//(руту не отдаем)
        int *rowCounts = (int*)malloc(sizeof(int) * workSize);

        for(int i = 0; i < workSize - 1; i++) {
            rowCounts[i] = numberOfVert / workSize;
        }
        //Если поровну делится то последнему столько же
        if(numberOfVert % workSize == 0) {
            rowCounts[workSize-1] = numberOfVert / workSize;
        }
            //Если не поровну, то в последний скидываем остатки
        else {
            rowCounts[workSize - 1] = numberOfVert - ((numberOfVert / workSize) * (workSize - 1));
        }

        cout << "Root process. Row counts sent to other processes" << endl;
        //Сказать каждому процессу сколько ему ждать строк
        for(int i = 1; i < size; i++) {
            MPI_Send(&rowCounts[i-1], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        cout << "Root process. Main loop start" << endl;
        //Основной цикл алгоритма Флойда
        for(int k = 0; k < numberOfVert; k++) {
            //Раздаем k-ую строку всем процессам(кроме рута)
            for(int p = 1; p < size; p++) {
                MPI_Isend(&matrix[k*numberOfVert], numberOfVert, MPI_INT, p, 0, MPI_COMM_WORLD, &request);
            }

            //Отдать каждому процессу его строки матрицы
            int num = 0; //Номер строки
            for(int i = 0; i < workSize; i++) { //Бежим по массиву с кол-вом строк
                for(int j = 0; j < rowCounts[i]; j++) { //Отдаем нужное количество строк
                    MPI_Isend(&matrix[num*numberOfVert], numberOfVert, MPI_INT, i+1, 0, MPI_COMM_WORLD, &request);
                    num++;
                }
            }

            //Получаем от процессов строки, измененные на этой итерации
            num = 0; //Номер строки
            for(int i = 0; i < workSize; i++) { //Бежим по массиву с кол-вом строк
                for(int j = 0; j < rowCounts[i]; j++) { //Принимаем нужное количество строк
                    MPI_Recv(&matrix[num*numberOfVert], numberOfVert, MPI_INT, i+1, 0, MPI_COMM_WORLD, &status);
                    num++;
                }
            }
        }

        cout << "Root process. Print final matrix:" << endl;
        printMatrix(matrix, numberOfVert);

        timeResult = MPI_Wtime() - time;
        cout << "Root process. Time spent: " << timeResult << endl;

        free(matrix);
    }
    else {
        //Получить количество вершин
        MPI_Bcast(&numberOfVert, 1, MPI_INT, ROOT_PROC, MPI_COMM_WORLD);

        //Получить количество строк, которые нужно будет изменить
        int rowCount;
        MPI_Recv(&rowCount, 1, MPI_INT, ROOT_PROC, 0, MPI_COMM_WORLD, &status);
        cout << rank << " : " << rowCount << endl;

        //Выделить память
        int *kRow = (int *)malloc(sizeof(int) * numberOfVert);
        int *rows = (int *)malloc(sizeof(int) * (rowCount * numberOfVert));

        //Основной цикл алгоритма Флойда
        for(int k = 0; k < numberOfVert; k++) {
            //Получить k-ю строку
            MPI_Recv(kRow, numberOfVert, MPI_INT, ROOT_PROC, 0, MPI_COMM_WORLD, &status);

            //Получить свои строки из матрицы
            for(int i = 0; i < rowCount; i++) {
                MPI_Recv(&rows[i*numberOfVert], numberOfVert, MPI_INT, ROOT_PROC, 0, MPI_COMM_WORLD, &status);
            }

            //Изменить нужные строки в матрице(параллельная работа)
            for (int i = 0; i < rowCount; i++) {
                for (int j = 0; j < numberOfVert; j++) {
                    rows[i*numberOfVert + j] = min(
                            rows[i*numberOfVert + j],
                            rows[i*numberOfVert + k] + kRow[j]
                    );
                }
            }

            //Отправить изменения руту
            for(int i = 0; i < rowCount; i++) {
                MPI_Send(&rows[i*numberOfVert], numberOfVert, MPI_INT, ROOT_PROC, 0, MPI_COMM_WORLD);
            }
        }

        free(kRow);
        free(rows);
    }

    MPI_Finalize();
    return 0;
}