//#include "mainwindow.h"
//#include <QApplication>
#include <math.h>
#include <string>

int fileRW(FILE* file, char o ){

}

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  MainWindow w;
  w.show();
  //--defination-->
    const long double Pi = 3.14159265358979323846264338327950288419716939937510;
    const long double F = 96485.3329;
    const long double R = 8.3144621;


    int Nmax = 4000;

    int I = 0;
    int J = 0;
    int N = 0;

    long double Psif[8100];
    long double PsiTotal[8100];
    long double Rqd[8100];
    long double Pota[8100];
    long double Dtau, Csi, CsiP;
    long double Tau(0), Tauf(0);
    long double u;
    long double Sum;
    long double Ct1, Ct2, Ct3, Ct4, Ct5;
    long double Lambda;
    long double Alpha = 0.5;
    long double PsiProv = 0;
    long double CsiPrime;
    long double Theta = 0; // Zero pour lest tests
    long double Gamma = 0; // Zero pour les tests
    long double Rau, CsiM;
    long double PsiP, PsiM;



    long double FRT;
    long double FSC0;
    long double DFRT;


    long double Ep1, Ep2;
    long double Ip1, Ip2;

    struct inputData{
      double Ei, Ef, E0;
      double S, C0, T, D;
      double Ru, Cd, ks;
      double v;
    }inputData;



 /*-----文件读取部分------*/

    FILE* file=NULL;

    file=fopen("input.txt","r");
    if(!feof(file) )
      {
        fseek(file,0L,SEEK_SET);
        //--------指针返回文件头
        fscanf(file,"lf",&inputData.Ei);
        fscanf(file,"lf",&inputData.Ef);
        fscanf(file,"lf",&inputData.E0);
        fscanf(file,"lf",&inputData.S);  inpuData.S/=10000;
        //-----------------------
        fscanf(file,"lf",&inputData.C0); inpuData.C0*=1000.0;
        fscanf(file,"lf",&inputData.T);
        fscanf(file,"lf",&inputData.D);  inputData.D/=10000;
        //------------------------
        fscanf(file,"lf",&inputData.Ru);
        fscanf(file,"lf",&inputData.Cd);
        fscanf(file,"lf",&inputData.ks); inputData.ks/=100;
        //-------------------
        fscanf(file,"lf",&inputData.v);

        //-------------关闭文件
//        if(feof(file)){fclose(file);}
      }
    fclose(file);


 //
    /*/-------------------------------------------
     * 本想改为用C++做流的输入输出，先用C语言实现
    ifstream file("input.txt", ios::in);

    if(!file.is_open())
      {cout << "Can not opening file";exit(1);}
    while(!file.eof()){
        inputData.
      }
    */

 /*-----文件读取部分结束------*/

/*----参数定义---*/
    FRT = F / R / inputData.T;
    FSC0 = F * inpuData.S * inpuData.C0;
    DFRT = inputData.D * FRT;

    //parameters are defined for oxidation
    u = -(inpuData.Ei - inpuData.E0) * FRT;
    CsiM = (inputData.Ef - inpuData.E0) * FRT;
    Lambda = inputData.ks / sqrt (v * DFRT);
    Rau = FRT * inputData.Ru * FSC0 * sqrt (v *DFRT);
    Gamma = inputData.Cd * sqrt (v) / FSC0 /sqrt (DFRT);


    //------------------------------------------
    printf (" ************************************************************* \n Ohmic Drop - Slow Charge Transfer - Version 3.2b 09/15/2016 \n " );
    printf ("                     CNRS / Univ. Rennes 1  \n +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n\n " );

    printf (" Ei = %f ; Ef = %f ; S = %e ;  C0 = %e  \n", inputData.Ei, inputData.Ef, inputData.S, inputData.C0);
    printf (" T = %f ; D = %e ; Ru = %e ;  Cd = %e  \n", inputData.T, inputData.D, inputData.Ru, inputData.Cd);
    printf (" ks = %f ; v = %f \n", inputData.ks, inputData.v);
    //-------------------------------------------

    Psif[0] = 0;
    Psif[1] = 0;
    Dtau = (CsiM + u) / Nmax; // calculation of the time step.
    Tau = Dtau;
    Csi = -u;
    Theta = Rau * Gamma;
    Ct3 = Theta/Dtau;

    //-----------------------
    printf (" Lambda = %LF ; Rau = %Lf ; Gamma = %Lf ; Dtau = %Lf \n", Lambda, Rau, Gamma, Dtau);
    //-----------------------


    //Rgd项
    for (I = 1 ; I <= 2.0 * Nmax + 1; I++)
      {
        Rqd[I] = (powl ((I + 1.0) , 1.5)) - (2.0 * powl (I, 1.5)) + (powl ((I - 1) , 1.5));
      }



    /*----------------数据扫描计算--------------------*/
 //Forward Scan
    for (I = 1 ; I <= Nmax ; I++)
      {
        Ct1 = expl (-Alpha * Csi) / Lambda;
        Ct2 = 4.0 / 3.0 * sqrtl (Dtau / Pi) * (1.0 + expl (-Csi));
        Pota[I] = Csi;

        Sum = 0;

        for (N = 1 ; N < I ; N++)

          {
            Sum = Sum + (Psif[N] * (1.0 + Ct3) - Psif[N - 1] * Ct3) * Rqd[I - N]; //P335末项的求和项
          }

        N = 0;

        Psif[I] = Psif[I - 1];

        do

          {
            Psif[I] = (Psif[I] + 100.0 * PsiProv) /101.0; // Il faut calmer la convergence //*********使Psif偏移约%1，记作PsiProv，然后用此Psif回代入方程，求解此条件下的Psif
            PsiProv = Psif[I];
            CsiPrime = Csi - Rau * Psif[I] - Theta * (1.0 - expl (-Tau/Theta));   //方程(13)
            Ct1 = expl (-Alpha * CsiPrime) / Lambda;
            Ct2 = 4.0 / 3.0 * sqrtl (Dtau / Pi) * (1.0 + expl (-CsiPrime)); //方程（11）的变形项

            Psif[I] = (1.0 - Sum * Ct2 + (Ct1 * Ct3 + Ct2 * Ct3) * Psif[I - 1] ) / (Ct1 + Ct2 + Ct1 * Ct3 + Ct2 * Ct3);

            N++; //???????????????为什么要在这里进行N++

          }
        while ( fabs ( Psif[I] - PsiProv ) > 0.0001 ); // 控制精度为delt<0.0001


        PsiTotal[I] = Psif[I] + Gamma * (1.0 - expl (-Tau / Theta));  //方程（12）

        if ((I % 200) == 0)

          {

            printf (" I= %d Csi= %Lf Psif= %Lf PsiTotal= %Lf N= %d \n", I, Csi, Psif[I], PsiTotal[I], N);

          }

        Csi = Csi + Dtau;
        Tau = Tau + Dtau;
      }

       Tauf = Tau - Dtau;




       // Reverse Scan

       for (I = Nmax + 1 ; I <= 2.0 * Nmax ; I++)
         {
           Ct1 = expl (-Alpha * Csi) / Lambda;
           Ct2 = 4.0 / 3.0 * sqrtl (Dtau / Pi) * (1.0 + expl (-Csi));

           Pota[I] = Csi;


           Sum = 0;
           for (N = 1 ; N < I ; N++)

             {
               Sum = Sum + (Psif[N] * (1.0 + Ct3) - Psif[N - 1] * Ct3) * Rqd[I - N];
             }

           N = 0;

           Psif[I] = Psif[I - 1]; // We start with psi calculated at N-1. This is the guest value.

           /*和Forward Scan的区别在这里*/
           if (Theta > 0.0001 ) // Do not need to consider the double layer if Theta is small

             {
               Ct4 = Theta * (2.0 * expl ((Tauf-Tau)/Theta) - expl (-Tau / Theta) - 1.0);
             }

           else

             {
               Ct4 = 0.0;

             }

           Ct5 = Gamma * (2.0 * expl ( (Tauf - Tau) / Theta) - expl (-Tau / Theta) - 1.0);
           /*------------------------------*/

           do

             {
               Psif[I] = (Psif[I] + 100.0 * PsiProv) /101.0; // This is required to dump the convergence.
               PsiProv = Psif[I];
               CsiPrime = Csi - Rau * Psif[I] -  Ct4;
               Ct1 = expl (-Alpha * CsiPrime) / Lambda;
               Ct2 = 4.0 / 3.0 * sqrtl (Dtau / Pi) * (1.0 + expl (-CsiPrime));

               Psif[I] = (1.0 - Sum * Ct2 + (Ct1 * Ct3 + Ct2 * Ct3) * Psif[I - 1] ) / (Ct1 + Ct2 + Ct1 * Ct3 + Ct2 * Ct3);

               N++;

             }

           while ( fabs ( Psif[I] - PsiProv ) > 0.0001 );


           PsiTotal[I] = Psif[I] + Ct5;


           if ((I % 200) == 0)
             {
               printf (" I= %d Csi= %Lf Psif= %Lf PsiTotal= %Lf N= %d \n", I, Csi, Psif[I], PsiTotal[I], N);
             }

           Csi = Csi - Dtau;
           Tau = Tau + Dtau;
         }
 /*----------------数据扫描计算部分结束--------------------*/
 /*----------------计算各种参数--------------------------*/
// Calculations of the maximum and minimum values, Csip, Psip, etc...

       PsiP = -100;
       CsiP = -500;
       for (N = 200 ; N <  Nmax; N++)
         {
           if ( PsiTotal[N] > PsiP )
             {
               PsiP = ( PsiTotal[N] + PsiTotal[N - 1] ) / 2.0;
               CsiP = Pota[N] + Dtau / 2.0;
               I = N;
             }
         }

       PsiM = 100;
       CsiM = 500;

       for (N = Nmax + 200 ; N < 2 * Nmax; N++)
         {
           if ( PsiTotal[N] < PsiM )
             {
               PsiM = ( PsiTotal[N] + PsiTotal[N - 1] ) / 2.0;
               CsiM = Pota[N] + Dtau/2.0;
               J = N;
             }
         }

       Ep1 = inpuData.E0 + CsiP / FRT;
       Ep2 = inpuData.E0 + CsiM / FRT;
       Ip1 = PsiP * FSC0 * sqrt(v * DFRT);
       Ip2 = PsiM * FSC0 * sqrt(v * DFRT);

       printf ("\n Calculations of peak potentials and peak currents :\n");
       printf (" IMax = %d ; CsiP = %Lf ; PsiP = %Lf \n", I, CsiP, PsiP);
       printf (" IMin = %d ; CsiM = %Lf ; PsiM = %Lf \n", J, CsiM, PsiM);
       printf (" Ep1 = %Lf V ; Ip1 = %Le A \n", Ep1, Ip1);
       printf (" Ep2 = %Lf V ; Ip2 = %Le A \n", Ep2, Ip2);

/*-----------------参数计算部分结束--------------*/
/*----------------写入文件--------------------*/
       file = open("results.txt","w");

       for (N = 1 ; N < 2 * Nmax - 1; N++)

         {
           if ((N % 20) == 0)

             {
               fprintf (fichier," %Lf ; %Le \n", inpuData.E0 + Pota[N]/FRT, PsiTotal[N] * FSC0 * sqrt(v * DFRT));
             }

         }

       fclose(fichier);
    return a.exec();
  }
