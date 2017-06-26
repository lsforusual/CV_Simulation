#ifndef INPUTDATA_H
#define INPUTDATA_H
#include<string>
using namespace::std;
class inputData
{
public:
  int setData(string str);
  string getData();
  int getRows(char *filename);

  double Ei, Ef,E0;
  double S, C0, T, D;
  double Ru, Cd, ks;
  double v;

  double realIp1;
  double realdE;

};

#endif // INPUTDATA_H
