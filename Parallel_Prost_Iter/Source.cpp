#include <mpi.h>
#include <iostream>
#include <cmath>

using namespace std;

double eps = 0.001;
static const int n = 3;
double masA[n][n] = { { 0.0, -0.1, -0.1 },
                      { -0.2, 0.0, -0.1 },
                      { -0.2, -0.2, 0.0 },
                      };

double masB[n] = { 1.2, 1.3, 1.4 };

long double x[n] = { 0.0, 0.0, 0.0}, x0[n] = { 0.0, 0.0, 0.0 }, max1[n] = { 0.0, 0.0, 0.0 };
int counter = 0, counter1[1] = { 0 };
double rbufA[n*n];
double rbuf[n*n];
double rbufB[n*n];


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
    MPI_Scatter(masB, 1, MPI_DOUBLE, rbufA, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    do
    {
        x[rank] = 0.0;
        for (int j = 0; j < n; j++) {
            x[rank] += masA[rank][j] * x0[j];
        };
        cout << "rank = " << rank << "  do + bufa x  = " << x[rank] << endl;
        cout << "rank = " << rank << "  do + bufa x0  = " << x0[rank] << endl;
        cout << "rank = " << rank << "  do + bufa x - x0  = " << x[rank] - x0[rank] << endl;
        cout << "rank = " << rank << "  do + bufa |x - x0|  = " << fabs(x[rank] - x0[rank]) << endl;

        x[rank] += rbufA[0];

        //cout << "rank = " << rank << "  posle + bufa x  = " << x[rank] << endl;
        //cout << "rank = " << rank << "  posle + bufa x0  = " << x0[rank] << endl;
        //cout << "rank = " << rank << "  posle + bufa x - x0  = " << x[rank] - x0[rank] << endl;
        cout << "rank = " << rank << "  posle + bufa pered obnuleniem max1 |x - x0|  = " << fabs(x[rank] - x0[rank]) << endl;
        max1[rank] = 0.0;


        if (max1[rank] <= fabs(x[rank] - x0[rank])) {
            max1[rank] = fabs(x[rank] - x0[rank]);
            cout << "rank = " << rank << " max1  = " << max1[rank] << endl;
        }
        x0[rank] = x[rank];
        /*     cout << "rank = " << rank << " x  = " << x[rank] << endl;
               cout << "rank = " << rank << " x0  = " << x0[rank] << endl;
               cout << "rank = " << rank << " x - x0  = " << x[rank] - x0[rank] << endl;
               cout << "rank = " << rank << " |x - x0|  = " << fabs(x[rank] - x0[rank]) << endl;*/
            // cout << "rank = " << rank << " max1 1 = " << max1[0] << endl;
        counter++;
        // cout << "rank = " << rank << "max1 = " << max1[rank] << endl;
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Gather(&max1[rank], 1, MPI_DOUBLE, rbufB, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        counter1[0] = 0;
        if (rank == 0)
        {
            cout << " counter Iter = " << counter << endl;
            //cout << "rank = " << rank << " max1 1 = " << rbufB[0] << endl;
            //cout << "rank = " << rank << " max1 2 = " << rbufB[1] << endl;
            //cout << "rank = " << rank << " max1 3 = " << rbufB[2] << endl;
            //cout << "rank = " << rank << " max1 4 = " << rbufB[3] << endl;
            for (int i = 0; i < n; i++) {
                if (rbufB[i] < eps)
                {
                    counter1[0]++;
                    cout << "rank = " << rank << " counter = " << counter1[0] << endl;
                    //cout << "rank = " << rank << " rbufB = " << rbufB[i] << endl;
                    fflush(stdout);
                }
            }
        }
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
    printf("vipolnenie zanyalo %f seconds\n", endtime - starttime);
    MPI_Finalize();
    return 1;
}