#include <mpi.h>
#include <iostream>
#include <cmath>

using namespace std;

double eps = 0.00001;   // pri zadannoy tochnosti, idelno tochniy otvet s zadannoy matrix poluchaetsya za 11 iterasiy na 3 proccesah
static const int n = 3;
double masA[n][n] = { { 0.0, -0.1, -0.1 },      //uzge privedennaya k nughnomu vidu matrix
                      { -0.2, 0.0, -0.1 },
                      { -0.2, -0.2, 0.0 },
                      };

double masB[n] = { 1.2, 1.3, 1.4 };   // stolbec svobodnih chlenov

long double x[n], x0[n], max1[n];
int counter = 0, counter1[1] = { 0 };
double rbufA[n * n];   // buffers neobhodimie dlya otpravki dannih mej rankov
double rbuf[n * n];
double rbufB[n * n];
double rbufC[n * n];
double rbufT[n * n];


int main() {
    MPI_Init(NULL, NULL);                        //nachalo
    double starttime, endtime,all_time[n];
    starttime = MPI_Wtime();    //uznaem vremya nachala

    int rank, size;

    //MPI_Status status;
    //MPI_Request request;
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // poluchaem kolichestvo proccesov
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // poluchaem rank

    for (int i = 0; i < n; i++) {    // zapolnyaem massive x0 svobodnimi chlenami
        x0[i] = masB[i];
    }
    MPI_Bcast(masA, n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);    // otpravlyaem vsu matrix na vse proccessi
    MPI_Bcast(masB, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);          // otpravlyaem ves vector svobodnih chlenov 

    do
    {

        if (rbufC != NULL)                      // tut zapolnyaetcya x0 na 2+ iter
        {
            for (size_t i = 0; i < n; i++)
            {
                x0[i] = rbufC[i];
            }
        }
       
        x[rank] = 0.0;                                 // x obnulyaetsya i vichislyaetsya zanovo na kajdoi iter na kajdom ranke svoy
        for (int j = 0; j < n; j++) {
            x[rank] += masA[rank][j] * x0[j];
        };

        x[rank] += masB[rank];

        max1[rank] = 0.0;                              //obnul max i posle emu prisvaevaetsya znachenie modulya raznosti

        if (max1[rank] < fabs(x[rank] - x0[rank])) {
            max1[rank] = fabs(x[rank] - x0[rank]);  
        }
        x0[rank] = x[rank];               

        MPI_Barrier(MPI_COMM_WORLD);   // programma jdet chtobi vse proccesi dowli do etogo momenta i tolko potom razrewaet idti dalwe
        MPI_Gather(&x0[rank], 1, MPI_DOUBLE, rbufC, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  //sobiraet vse znacheniya x0 so vseh rankov
        MPI_Bcast(rbufC, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);   // otpravlyaet poluchennie znacheniya x0 na vse ranki
        MPI_Gather(&max1[rank], 1, MPI_DOUBLE, rbufB, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // sobir vse max so vseh rankov
        counter++;  // schetchik iteratciy
        counter1[0] = 0;
        if (rank == 0)   // proverka chto vse  x naydenni s zadannoy tochnost, esli hot' odin ne tochen, to provoditsya ewe 1 iter
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
        MPI_Bcast(counter1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  // kak vse poschitanno otpravlyaet kol tochnih znach na vse ranki
    } while (counter1[0] != n);  //esli kol tochnih ravno kol iskomih, then cicle ends.

    MPI_Gather(&x[rank], 1, MPI_DOUBLE, rbuf, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  // sobiraet finalniye x na 0 ranke

    if (rank == 0)  // vivodit kol iter i otvet v console
    {
        cout << endl << "Kol iter: " << counter << endl << endl;
        for (size_t i = 0; i < n; i++)
        {
            cout << "x" << i + 1 << "=" << rbuf[i] << " " << endl;
        }
    }

    endtime = MPI_Wtime();  // polushem vremya zaverweniya
    all_time[rank] = endtime - starttime; // uznayem vremya raboti dlya kajdogo proccessa
    MPI_Gather(&all_time[rank], 1, MPI_DOUBLE, rbufT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // sobiraem ih na 0 ranke 
    if (rank==0)
    {
        rbufT[0] += rbufT[1] + rbufT[2];  //skladivaem i vivodim v console obwee vremya raboti
        printf("vipolnenie zanyalo %f seconds\n", rbufT[0]);
    }
    MPI_Finalize(); // The end.
    return 1;
}