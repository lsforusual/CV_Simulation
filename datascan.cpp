#include<inputdata.h>
#include<outputdata.h>
#include<windows.h>
#include<cmath>
#include<string>
#include<vector>
#include<algorithm>
#include<iostream>
#include<fstream>

using namespace::std;


outputData dataScan(inputData idata)
{

	//Count the time the programme spend
	int long startTime = GetTickCount();

	//Dimentionless

	outputData odata;

	// sepcify the const
	const double Pi = 3.14159265358979323846264338327950288419716939937510;
	const double F = 96485.3329;
	const double R = 8.3144621;
	const double T = 298.15;
	const double FRT = F/(R*T);
	//Specify simulation parameters
	double theta_i = 20;//FRT * (idata.Ei - idata.E0);//!!!!!!!!!!!
	double theta_v = -20;//FRT * (idata.Ef - idata.E0);//!!!!!!!!!!!!!!!!
	bool flag = theta_i > theta_v;  //add a flag in case of theta_i < theta_v;
	//  double deltaX = 2E-4;
	double deltaTheta = 0.02;
	double Ru = idata.Ru;
	double Cd = idata.Cd;
	double D = idata.D;
	double A = idata.S;
	double sigma = idata.v * (FRT) * (A/(Pi*D));
	double C0 = idata.C0;


	//Calculate other parameters
	double deltaT = deltaTheta / sigma;
	double maxT = 2 * std::abs(theta_v - theta_i) / sigma;
	double maxX = 6.0 * sqrt(maxT);
	//	int iterX = (int) (maxX / deltaX) ; //number of spacesteps
	int iterT = (int) (maxT / deltaT);  //number of timesteps

	//Expanding grid
	double omega = 1.01;
	double h = 1E-4;

	std::vector<double> X;
	X.push_back(0.0);

	while(X.back() < maxX)
	{
		X.push_back(X.back() + h);
		h *= omega;
	}
	int iterX = X.size();

	//Calculate Thomas coefficients
	std::vector<double> alpha(iterX,0);
	std::vector<double> beta(iterX,1);
	std::vector<double> gamma(iterX,0);
	for(int i=1; i<iterX-1; i++)
	{
		double deltaX_p = X[i+1] - X[i];
		double deltaX_m = X[i] - X[i-1];

		alpha[i] = -2 * deltaT / (deltaX_m * deltaX_m + deltaX_p * deltaX_m);
		gamma[i] = -2 * deltaT / (deltaX_p * deltaX_p + deltaX_p * deltaX_m);
		beta[i] = 1 - alpha[i] -gamma[i];
	}

	////Create containers
	std::vector<double> g_mod(iterX, 0);
	std::vector<double> C(iterX, 1.0); //concentration profile

	//Modify gamma coefficients
	g_mod[0] = 0;
	for(int i=1; i<iterX-1; i++)
		g_mod[i] = gamma[i] / (beta[i] - g_mod[i-1] * alpha[i]);


	//BEGIN SIMULATION
	std::vector<double> flux(iterT,0);

	//Calculate only C[A,0] C[A,1]
	double Theta = theta_i;
	double a = 0.5;
	double K0 =  idata.ks * sqrt(A/Pi) / D;

	for(int k=0; k<iterT; k++)
	  {
		flag ? (k<iterT/2 ? Theta -= deltaTheta : Theta += deltaTheta):
			(k<iterT/2 ? Theta += deltaTheta : Theta -= deltaTheta);

		//Forward sweep -- creat modified deltas
		double BV_Theta = exp(-a*Theta);

		alpha[0] = 0;
		beta[0] = 1 + h * BV_Theta * K0 * (1+exp(Theta));
		gamma[0] = -1;
		// C[0] = 1.0 / (1.0 + exp(-Theta)); //Create delta[0]
		C[0] = h * BV_Theta * K0 *exp(Theta);
		C[0] /= beta[0]; //Create modified delta[0]

		for(int i=1; i<iterX; i++)
			C[i] = ( C[i] - C[i-1] * alpha[i] ) / ( beta[i] - g_mod[i-1] * alpha[i] );

		//Back Substitution
		C[iterX-1] = 1.0;
		for(int i=iterX-2; i>=0; i--)
			C[i] = C[i] - g_mod[i] * C[i+1];

		//Output current
		flux[k] = -(C[1] - C[0] ) / h ;
	}

//	alpha.resize(iterT);
//	beta.resize(iterT);
//	gamma.resize(iterT);
//	g_mod.resize(iterT);

//	//Recalculate Thomas coefficient for totalI
//	alpha.assign(iterT,0);
//	beta.assign(iterT, deltaT-Ru*Cd);
//	gamma.assign(iterT,Ru*Cd);
//	g_mod.assign(iterT,0);

//	g_mod[0] = 0;
//	for(int i=1; i<iterT-1; i++)
//		g_mod[i] = gamma[i] / (beta[i] - g_mod[i-1] * alpha[i]);

//	//Create Thomas coefficient delta WITH dimension
//	double factor1 = D*Pi/A * Cd * sigma /FRT * deltaT;//idata.v * idata.Cd * deltaT
//	double factor2 = sqrt(Pi*A)*F*D*C0*deltaT;
//	std::vector<double> totalI(iterT,0);
//	for(int t=0; t<iterT; t++)
//		totalI[t] = factor1+factor2*flux[t];

//	//Calculate the modified delta
//	flag ? totalI[0] =  factor1/deltaT : totalI[0] =  -factor1/deltaT ; //idata.v * idata.Cd
//	for(int i=1; i<iterT; i++)
//		totalI[i] = (totalI[i]-totalI[i-1]*alpha[i]) / -( beta[i] - g_mod[i-1] * alpha[i] );

//	//Caculate the totalI
//	flag ? totalI[iterT-1] = -D*Pi/A * Cd * sigma / FRT : totalI[iterT-1] = D*Pi/A * Cd * sigma / FRT;
//	for(int t=iterT-2; t>=0; t--)
//	    totalI[t] = totalI[t] - g_mod[t]*totalI[t+1]; //it cannot work correctly for some reason
//	//for(int i=1; i<iterT-1; i++)
//	//	totalI[i]=totalI[i]-beta[i]*totalI[i-1]/gamma[i];

	std::vector<double> totalI(iterT,0);
	for(int i=0;i<iterT;i++)
	  totalI[i] = flux[i]*sqrt(A*Pi)*F*D;
	//Open file to output CV
	Theta = theta_i;
	std::ofstream CV("CV_output.txt");
	for(int i=0; i<iterT; i++)
	{
	    flag ? (i<iterT/2 ? Theta -= deltaTheta : Theta += deltaTheta):
		   (i<iterT/2 ? Theta += deltaTheta : Theta -= deltaTheta);
		CV << (Theta/FRT + idata.E0) << "\t" << totalI[i] << "\n";
	}
	CV.close();

	Theta = theta_i;


	if(iterT<4000)
	{
		for(int i=0; i<iterT; i++)
		{
			i<iterT/2 ? Theta -= deltaTheta:Theta += deltaTheta;

			odata.plotdata[0][i] = Theta / FRT + idata.E0 ;
			odata.plotdata[1][i] = totalI[i];
		}
		for(int i=iterT; i<4000; i++)
		{
			odata.plotdata[0][i] = 0 ;
			odata.plotdata[1][i] = 0;
		}
	}else{
		for(int i=0; i<iterT/2; i++)
		{
			i<iterT/2 ? Theta -= deltaTheta:Theta += deltaTheta;

			odata.plotdata[0][i] = Theta / FRT - idata.E0 ;
			odata.plotdata[1][i] = totalI[i];
		}
		for(int i=iterT/2; i<4000; i++)
		{
			odata.plotdata[0][i] = 0 ;
			odata.plotdata[1][i] = 0;
		}
	}


	//END SIMULATION

	int long timeElapsed = GetTickCount() - startTime;
	double B4 = timeElapsed / (1000.0 * iterX *iterT);

	std::cout << B4 << std::endl;
	std::vector<double>::iterator Ip1=std::max_element(totalI.begin(),totalI.end());
	odata.Ip1 = *Ip1;
	odata.Ep1 = (theta_i-(std::distance(totalI.begin(),Ip1)*deltaTheta))/FRT + idata.E0;
	std::vector<double>::iterator Ip2=std::min_element(totalI.begin(),totalI.end());
	odata.Ip2 = *Ip2;
	odata.Ep2 = (theta_v+(std::distance(totalI.begin(),Ip2)*deltaTheta)-(theta_i-theta_v))/FRT + idata.E0;

	odata.DEp = odata.Ep1 - odata.Ep2;

	return odata;
}

/*
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

   //----参数定义---
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

     // -----------------------
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



  //----------------数据扫描计算--------------------
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

         //和Forward Scan的区别在这里
         if (Theta > 0.0001 ) // Do not need to consider the double layer if Theta is small

           {
             Ct4 = Theta * (2.0 * expl ((Tauf-Tau)/Theta) - expl (-Tau / Theta) - 1.0);
           }

         else

           {
             Ct4 = 0.0;

           }

         Ct5 = Gamma * (2.0 * expl ( (Tauf - Tau) / Theta) - expl (-Tau / Theta) - 1.0);
         //-----------------------------

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

		
         //if ((I % 200) == 0)
         //  {
         //    //printf (" I= %d Csi= %Lf Psif= %Lf PsiTotal= %Lf N= %d \n", I, Csi, Psif[I], PsiTotal[I], N);
         //    str << "I = "        << I           <<";";
         //    str << "Csi = "      << Csi         <<";";
         //    str << "PsiTotal = " << PsiTotal[I] <<";";
         //    str << "N = "        << N           <<";";
         //  }

         Csi = Csi - Dtau;
         Tau = Tau + Dtau;
       }
	  
     //---------------数据扫描计算部分结束--------------------


     //----------------计算各种参数--------------------------
     //Calculations of the maximum and minimum values, Csip, Psip, etc...

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


     //      printf ("\n Calculations of peak potentials and peak currents :\n");
     //      printf (" IMax = %d ; CsiP = %Lf ; PsiP = %Lf \n", IMax, CsiP, PsiP);
		   //printf (" IMin = %d ; CsiM = %Lf ; PsiM = %Lf \n", IMin, CsiM, PsiM);
     //      printf (" Ep1 = %Lf V ; Ip1 = %Le A \n", Ep1, Ip1);
     //      printf (" Ep2 = %Lf V ; Ip2 = %Le A \n", Ep2, Ip2);

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


           //-----------------参数计算部分结束--------------//
           //-----------------输出数据点-------------------//

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

*/
