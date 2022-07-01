#include <mpi.h>
#include <iostream>
#include <cmath>

using namespace std;

double eps = 0.001;
static const int n = 4;
double masA[n][n] = { { 0.05, -0.06, -0.12, 0.14 },
                      { 0.04, -0.12, 0.68, 0.11 },
                      { 0.34, 0.08, -0.06, 0.44 },
                      { 0.11, 0.12, -0.03, -0.8 } };

double masB[n] = { -2.17, 1.4, -2.1, -0.8 };

double x[n] = { 0.0, 0.0, 0.0, 0.0 }, x0[n] = { 0.0, 0.0, 0.0, 0.0 }, max1[n] = { 0.0, 0.0, 0.0, 0.0 };
int counter = 0, counter1[1] = { 0 };
double rbufA[n];
double rbuf[n];
double rbufB[n];


int main() {
    MPI_Init(NULL, NULL);
    double starttime, endtime;
    starttime = MPI_Wtime();

    int rank, size;

    //MPI_Status status;
    //MPI_Request request;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  

    for (int i = 0; i < n; i++) {
        x0[i] = masB[i];
    }

    MPI_Bcast(masA, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(masB, n, MPI_DOUBLE, rbufA, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    do
    {
        x[rank] = 0.0;
        for (int j = 0; j < n; j++) {
            x[rank] += masA[rank][j] * x0[j];
        }


        x[rank] += rbufA[rank];
       
        max1[rank] = 0.0;


        if (max1[rank] < fabs(x[rank] - x0[rank])) {
            max1[rank] = fabs(x[rank] - x0[rank]);
        }
        x0[rank] = x[rank];

        counter++;
       // cout << "rank = " << rank << "max1 = " << max1[rank] << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(&max1[0], n, MPI_DOUBLE, rbufB, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        counter1[0] = 0;
        if (rank == 0)
        {
            cout << "rank = " << rank << " max1 1 = " << max1[0] << endl;
            cout << "rank = " << rank << " max1 2 = " << max1[1] << endl;
            cout << "rank = " << rank << " max1 3 = " << max1[2] << endl;
            cout << "rank = " << rank << " max1 4 = " << max1[3] << endl;
            for (int i = 0; i < 4; i++) {
                if (rbufB[i] <= eps)
                {
                    counter1[0]++;
                    cout << "rank = "<< rank << " counter = "<<counter1[0] << endl;
                    cout << "rank = " << rank << " max1 = " << rbufB[i] << endl;
                    fflush(stdout);
                }
            }            
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(counter1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    } while (counter1[0] != n);

    MPI_Gather(&x[0], n, MPI_DOUBLE, rbuf, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        cout << endl << "Kol iter: " << counter << endl << endl;
        for (size_t i = 0; i < n; i++)
        {
            cout << "x" << i + 1 << "=" << rbuf[i] << " " << endl;
        }

    }



    endtime = MPI_Wtime();
    printf("vipolnenie zanyalo %f seconds\n", endtime - starttime);
    MPI_Finalize();
    return 1;
}