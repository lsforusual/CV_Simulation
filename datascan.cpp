#include<inputdata.h>
#include<outputdata.h>
#include<math.h>
#include<string>

using namespace::std;

outputData dataScan(inputData idata)
{
  const long double Pi = 3.14159265358979323846264338327950288419716939937510;
  const long double F = 96485.3329;
  const long double R = 8.3144621;
  int I,N;
  int Nmax = 4000;

  long double Psif[8100];
  long double PsiTotal[8100];
  long double Rqd[8100];
  long double Pota[8100];
  long double Dtau, Csi;
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




  long double FRT;
  long double FSC0;
  long double DFRT;

  outputData odata;

  /*----参数定义---*/
      FRT = F / R / idata.T;
      FSC0 = F * idata.S * idata.C0;
      DFRT = idata.D * FRT;

      //parameters are defined for oxidation
      u = -(idata.Ei - idata.E0) * FRT;
      CsiM = (idata.Ef - idata.E0) * FRT;
      Lambda = idata.ks / sqrt (idata.v * DFRT);
      Rau = FRT * idata.Ru * FSC0 * sqrt (idata.v *DFRT);
      Gamma = idata.Cd * sqrt (idata.v) / FSC0 /sqrt (DFRT);


      Psif[0] = 0;
      Psif[1] = 0;
      Dtau = (CsiM + u) / Nmax; // calculation of the time step.
      Tau = Dtau;
      Csi = -u;
      Theta = Rau * Gamma;
      Ct3 = Theta/Dtau;

      //-----------------------
      //printf (" Lambda = %LF ; Rau = %Lf ; Gamma = %Lf ; Dtau = %Lf \n", Lambda, Rau, Gamma, Dtau);
      //-----------------------
      odata.Lambda=Lambda;
      odata.Rau=Rau;
      odata.Gamma=Gamma;
      odata.Dtau=Dtau;

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
          Psif[I] = (Psif[I] + 100.0 * PsiProv) /101.0; //*********使Psif偏移约%1，记作PsiProv，然后用此Psif回代入方程，求解此条件下的Psif
          PsiProv = Psif[I];
          CsiPrime = Csi - Rau * Psif[I] - Theta * (1.0 - expl (-Tau/Theta));   //方程(13)
          Ct1 = expl (-Alpha * CsiPrime) / Lambda;
          Ct2 = 4.0 / 3.0 * sqrtl (Dtau / Pi) * (1.0 + expl (-CsiPrime)); //方程（11）的变形项

          Psif[I] = (1.0 - Sum * Ct2 + (Ct1 * Ct3 + Ct2 * Ct3) * Psif[I - 1] ) / (Ct1 + Ct2 + Ct1 * Ct3 + Ct2 * Ct3);

          N++; //------?为什么要在这里进行N++

        }
      while ( fabs ( Psif[I] - PsiProv ) > 0.0001 ); // 控制精度为delt<0.0001


      PsiTotal[I] = Psif[I] + Gamma * (1.0 - expl (-Tau / Theta));  //方程（12）

//      if ((I % 200) == 0)

//        {

//          printf (" I= %d Csi= %Lf Psif= %Lf PsiTotal= %Lf N= %d \n", I, Csi, Psif[I], PsiTotal[I], N);

//        }

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


//         if ((I % 200) == 0)
//           {
//             //printf (" I= %d Csi= %Lf Psif= %Lf PsiTotal= %Lf N= %d \n", I, Csi, Psif[I], PsiTotal[I], N);
//             str << "I = "        << I           <<";";
//             str << "Csi = "      << Csi         <<";";
//             str << "PsiTotal = " << PsiTotal[I] <<";";
//             str << "N = "        << N           <<";";
//           }

         Csi = Csi - Dtau;
         Tau = Tau + Dtau;
       }
/*----------------数据扫描计算部分结束--------------------*/


     /*----------------计算各种参数--------------------------*/
    // Calculations of the maximum and minimum values, Csip, Psip, etc...

           long double PsiP = -100;
           long double CsiP = -500;
           int IMax=0;
           int IMin=0;
           for (int N = 200 ; N <  Nmax; N++)
             {
               if ( PsiTotal[N] > PsiP )
                 {
                   PsiP = ( PsiTotal[N] + PsiTotal[N - 1] ) / 2.0;
                   CsiP = Pota[N] + Dtau / 2.0;
                   IMax = N;
                 }
             }

            long double PsiM = 100;
            CsiM = 500;

           for (int N = Nmax + 200 ; N < 2 * Nmax; N++)
             {
               if ( PsiTotal[N] < PsiM )
                 {
                   PsiM = ( PsiTotal[N] + PsiTotal[N - 1] ) / 2.0;
                   CsiM = Pota[N] + Dtau/2.0;
                   IMin = N;
                 }
             }
           //-----------------------------------------

           long double Ep1 = idata.E0 + CsiP / FRT;
           long double Ep2 = idata.E0 + CsiM / FRT;
           long double Ip1 = PsiP * FSC0 * sqrt(idata.v * DFRT);
           long double Ip2 = PsiM * FSC0 * sqrt(idata.v * DFRT);


   //        printf ("\n Calculations of peak potentials and peak currents :\n");
   //        printf (" IMax = %d ; CsiP = %Lf ; PsiP = %Lf \n", I, CsiP, PsiP);
   //        printf (" IMin = %d ; CsiM = %Lf ; PsiM = %Lf \n", J, CsiM, PsiM);
   //        printf (" Ep1 = %Lf V ; Ip1 = %Le A \n", Ep1, Ip1);
   //        printf (" Ep2 = %Lf V ; Ip2 = %Le A \n", Ep2, Ip2);

           //----------------------------
           odata.IMax=IMax;
           odata.IMin=IMin;
           odata.Ip1=Ip1;
           odata.Ip2=Ip2;
           odata.Ep1=Ep1;
           odata.Ep2=Ep2;
           odata.PsiM=PsiM;
           odata.PsiP=PsiP;
           odata.CsiM=CsiM;
           odata.CsiP=CsiP;
           odata.DEp=Ep1-Ep2;


           /*-----------------参数计算部分结束--------------*/
           /*-----------------输出数据点-------------------*/

           for (int i=0,N = 1; (N < 2 * Nmax - 1); N++)
             {

               if ((N % 20) == 0)
                 {
                   odata.plotdata[0][i] = idata.E0 + Pota[N]/FRT;
                   odata.plotdata[1][i] = PsiTotal[N] * FSC0 * sqrt(idata.v * DFRT);
                   i++;
                 }

             }

     return odata;
}
