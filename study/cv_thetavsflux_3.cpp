#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
//on Linux
//#include <time.h>
#include <windows.h>

int main()
{
  //Count the time the programme spend
  //        int long startTime = GetTickCount();

  // sepcify the const
  const double Pi = 3.14159265358979323846264338327950288419716939937510;
  const double F = 96485.3329;
  const double R = 8.3144621;
  const double T = 298.15;
  //Specify simulation parameters
  double theta_i = 20.0;
  double theta_v = -20.0;
  double sigma = 100.0;
  double deltaX = 2E-4;
  double deltaTheta = 0.02;
  double Ru = 424.0;
  double Cd = 6.7e-8;
  double D = 1.0;
  double A = 0.0102;

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

  //Create containers
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
  double K0 = 1.0;

  for(int k=0; k<iterT; k++)
    {
      if( k<iterT/2 )
        Theta -= deltaTheta;
      else
        Theta += deltaTheta;

      //Forward sweep -- creat modified deltas
      double BV_Theta = exp(-a*Theta);

      alpha[0] = 0;
      beta[0] = 1 + h * BV_Theta * K0 * (1+exp(Theta));
      gamma[0] = -1;
      // C[0] = 1.0 / (1.0 + exp(-Theta)); //Create delta[0]
      C[0] = h * BV_Theta * K0 *exp(Theta);
      C[0] /= beta[0]; //Create modified delta[0]

      for(int i=1; i<iterX+1; i++)
        C[i] = ( C[i] - C[i-1] * alpha[i] ) / ( beta[i] - g_mod[i-1] * alpha[i] );

      //Back Substitution
      C[iterX-1] = 1.0;
      for(int i=iterX-2; i>=0; i--)
        C[i] = C[i] - g_mod[i] * C[i+1];

      //Output current
      flux[k] = -(C[1] - C[0] ) / h ;
    }

  alpha.resize(iterT);
  beta.resize(iterT);
  gamma.resize(iterT);
  g_mod.resize(iterT);

  //Recalculate Thomas coefficient for totalI
  for(int i=0; i<iterT; i++)
    {
      alpha[i] = 0;
      beta[i] = deltaT-Ru*Cd;
      gamma[i] = Ru*Cd;
    }

  g_mod[0] = 0;
  for(int i=1; i<iterT-1; i++)
    g_mod[i] = gamma[i] / (beta[i] - g_mod[i-1] * alpha[i]);

  //Create Thomas coefficient delta
  double factor1 = R*T/F * D*Pi/A * Cd * sigma * deltaT;
  double factor2 = sqrt(Pi*A)*F*D*C[iterX-1]*deltaT;
  std::vector<double> totalI(iterT,0);
  for(int t=0; t<iterT; t++)
    totalI[t] = factor1+factor2*flux[t];

  //Calculate the modified delta
  totalI[0] =  R*T/F * D*Pi/A * Cd * sigma;
  for(int i=1; i<iterT+1; i++)
    totalI[i] = (totalI[i]-totalI[i-1]*alpha[i]) / ( beta[i] - g_mod[i-1] * alpha[i] );

  //Caculate the totalI
  totalI[iterT-1] = R*T/F * D*Pi/A * Cd * sigma;
  for(int t=iterT-2; t>=0; t--)
    totalI[t] = totalI[t] - g_mod[t]*totalI[t+1];


  //Open file to output CV
  Theta = theta_i;
  std::ofstream CV("CV_output.txt");
  for(int i=0; i<iterT; i++)
    {
      if( i<iterT/2 )
        Theta -= deltaTheta;
      else
        Theta += deltaTheta;
      CV << Theta << "\t" << totalI[i] << "\n";
    }
  //END SIMULATION
  CV.close();

  //        int long timeElapsed = GetTickCount() - startTime;
  //        double B4 = timeElapsed / (1000 * iterX *iterT);

  //        std::cout << B4 << std::endl;

  exit(0);

}
