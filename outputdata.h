#ifndef OUTPUTDATA_H
#define OUTPUTDATA_H
#include<string>
using namespace::std;
class outputData
{
public:
  string getData();
  string format_GetData();
  string getTitle();
  string getPlotData();

    long double Lambda;
    long double Rau;
    long double Gamma;
    long double Dtau;

    long double IMax;
    long double IMin;

    long double CsiP;
    long double CsiM;

    long double PsiP;
    long double PsiM;

    long double Ep1;
    long double Ip1;

    long double Ep2;
    long double Ip2;

    long double DEp;

    long double ks;
	long double S;

    long double plotdata[2][4000];



};
#endif // OUTPUTDATA_H
