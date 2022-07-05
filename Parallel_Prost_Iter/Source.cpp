#include <mpi.h>
#include <iostream>
#include <cmath>

using namespace std;

double eps = 0.00001;   // pri zadannoy tochnosti, idelno tochniy otvet s zadannoy matrix poluchaetsya za 11 iterasiy na 3 proccesah
static const int n = 3;
double masA[n][n] = { { 0.0, -0.1, -0.1 },
                      { -0.2, 0.0, -0.1 },
                      { -0.2, -0.2, 0.0 },
                      };

double masB[n] = { 1.2, 1.3, 1.4 };

long double x[n], x0[n], max1[n];
int counter = 0, counter1[1] = { 0 };
double rbufA[n * n];
double rbuf[n * n];
double rbufB[n * n];
double rbufC[n * n];
double rbufT[n * n];


int main() {
    MPI_Init(NULL, NULL);
    double starttime, endtime,all_time[n];
    starttime = MPI_Wtime();

    int rank, size;

    //MPI_Status status;
    //MPI_Request request;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int i = 0; i < n; i++) {
        x0[i] = masB[i];
    }

    do
    {
        MPI_Bcast(masA, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(masB, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (rbufC != NULL)
        {
            for (size_t i = 0; i < n; i++)
            {
                x0[i] = rbufC[i];
            }
        }
       
        x[rank] = 0.0;
        for (int j = 0; j < n; j++) {
            x[rank] += masA[rank][j] * x0[j];
        };

        x[rank] += masB[rank];

        max1[rank] = 0.0;

        if (max1[rank] < fabs(x[rank] - x0[rank])) {
            max1[rank] = fabs(x[rank] - x0[rank]);  
        }
        x0[rank] = x[rank];

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(&x0[rank], 1, MPI_DOUBLE, rbufC, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(rbufC, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&max1[rank], 1, MPI_DOUBLE, rbufB, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        counter++;
        counter1[0] = 0;
        if (rank == 0)
        {
            for (int i = 0; i < n; i++) {
                if (rbufB[i] < eps)
                {
                    counter1[0]++;
                }
            }
        }
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(counter1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } while (counter1[0] != n);

    MPI_Gather(&x[rank], 1, MPI_DOUBLE, rbuf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        cout << endl << "Kol iter: " << counter << endl << endl;
        for (size_t i = 0; i < n; i++)
        {
            cout << "x" << i + 1 << "=" << rbuf[i] << " " << endl;
        }
    }

    endtime = MPI_Wtime();
    all_time[rank] = endtime - starttime;
    MPI_Gather(&all_time[rank], 1, MPI_DOUBLE, rbufT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank==0)
    {
        rbufT[0] += rbufT[1] + rbufT[2];
        printf("vipolnenie zanyalo %f seconds\n", rbufT[0]);
    }
    MPI_Finalize();
    return 1;
}