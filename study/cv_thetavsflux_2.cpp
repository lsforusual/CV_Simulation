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

  //Specifiy simulation parameters
  double theta_i = 20.0;
  double theta_v = -20.0;
  double sigma = 100.0;
  double deltaX = 2E-4;
  double deltaTheta = 0.02;

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
  std::vector<double> beta(iterX,0);
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

  //Open file to output CV
  std::ofstream CV("CV_output.txt");

  //BEGIN SIMULATION
  double Theta = theta_i;

  for(int k=0; k<iterT; k++)
    {
      if( k<iterT/2 )
        Theta += deltaTheta;
      else
        Theta -= deltaTheta;

      //Forward sweep -- creat modified deltas
      C[0] = 1.0 / (1.0 + exp(-Theta));
      for(int i=1; i<iterX+1; i++)
        C[i] = ( C[i] - C[i-1] * alpha[i] ) / ( beta[i] - g_mod[i-1] * alpha[i] );

      //Back Substitution
      C[iterX-1] = 1.0;
      for(int i=iterX-2; i>=0; i--)
        C[i] = C[i] - g_mod[i] * C[i+1];

      //Output current
      double flux = -( -C[2] + 4 * C[1] - 3 * C[0] ) / (2*deltaX);
      CV << Theta << "\t" << flux << "\n";
    }


  //END SIMULATION
  CV.close();

//        int long timeElapsed = GetTickCount() - startTime;
//        double B4 = timeElapsed / (1000 * iterX *iterT);

//        std::cout << B4 << std::endl;

  exit(0);

}
