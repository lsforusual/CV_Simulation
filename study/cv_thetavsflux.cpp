#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
//on Linux
//#include <time.h>
/*
//on Windows
#include <windows.h>
*/ 

int main()
{
	//Count the time the programme spend
//	long int startTime = clock();

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
	int iterX = (int) (maxX / deltaX) ; //number of spacesteps
	int iterT = (int) (maxT / deltaT);  //number of timesteps

	//Calculate Thomas coefficients
	double lambda = deltaT / (deltaX * deltaX);
	double alpha = -lambda;
	double beta = 2 * lambda + 1.0;
	double gamma = -lambda;

	//Create containers
	std::vector<double> g_mod(iterX, 0);
	std::vector<double> C(iterT, 1.0); //concentration profile

	//Modify gama coefficients
	g_mod[0] = 0;
	for(int i=1; i<iterX-1; i++)
		g_mod[i] = gamma / (beta - g_mod[i-1] * alpha);

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
			 C[i] = ( C[i] - C[i-1] * alpha ) / ( beta - g_mod[i-1] * alpha );

		 //Back Substitution
		 C[iterX-1] = 1.0;
		 for(int i=iterX-1; i>=0; i--)
			 C[i] = C[i] - g_mod[i] * C[i+1];

		 //Output current
		 double flux = -( -C[2] + 4 * C[1] - 3 * C[0] ) / (2.0 * deltaX);
		 CV << Theta << "\t" << flux << "\n";
	 }

	//END SIMULATION
	 CV.close();

//	long int timeElapsed = clock() - startTime;
//	double B4 = timeElapsed / (1000 * iterX *iterT);
	
//	std::cout << B4 << std::endl;

}
