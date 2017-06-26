#include<iostream>
#include<vector>
#include<string>
#include<inputdata.h>
#include<outputdata.h>
#include<fstream>
#include<stdio.h>
#include<sstream>
#include<iomanip>
using namespace::std;



int inputData::setData(string str="0 0 0 0 0 0 0 0 0 0 0 0 0")
{
  std::stringstream ss(str);

  int i=0;

  while(i <= 12)
    {
      switch (i % 13)
        {
        case 0:
          ss >> this->Ei;
          cout << "Ei =" << Ei <<";";
          break;
        case 1:
          ss >> this->Ef;
          cout << "Ef = " << Ef << ";" ;
          break;
        case 2:
          ss >> this->E0;
          cout << "E0 = " << E0 << ";" ;
          break;
        case 3:
          ss >> this->S;
          this->S/=10000;
          cout << "S = "  << this->S  << ";" ;
          break;
        case 4:
          ss >> this->C0;
          this->C0*=1000.0;
          cout << "C0 = " << this->C0 << ";" ;
          cout << "\n";
          break;
        case 5:
          ss >> this->T;
          cout << "T = "  << this->T  << ";" ;
          break;
        case 6:
          ss >> this->D;
          this->D/=10000;
          cout << "D = "  << this->D  << ";" ;
          break;
        case 7:
          ss >> this->Ru;
          cout << "Ru = " << this->Ru << ";" ;
          break;
        case 8:
          ss >> this->Cd;
          cout << "Cd = " << this->Cd << ";" ;
          cout << "\n";
          break;
        case 9:
          ss >> this->ks;
          this->ks/=100;
          cout << "ks = " << this->ks << ";" ;
          break;
        case 10:
          ss >> this->v;
          cout << "v = "  << this->v  << ";" ;
          cout << endl;
          break;
        case 11:
          ss >> this->realdE;
          cout << "realdE = "  << this->realdE  << ";" ;
          cout << endl;
          break;        
		case 12:
          ss >> this->realIp1;
          cout << "realIp1 = "  << this->realIp1  << ";" ;
          cout << endl;
          break;
        default:
          cout << "Can not split the data correctly, please check the dataformat!!" <<endl;
          return 0;
        }

      i++;
    }

  if(str =="0 0 0 0 0 0 0 0 0 0 0")
    {
      cout << "idata.data is set as default 0. " <<endl;
	  return 1;
    }


}
string inputData::getData()
{
  std::stringstream str;
  str << "Ei = " << Ei << ";" ;
  str << "Ef = " << Ef << ";" ;
  str << "S = "  << S  << ";" ;
  str << "C0 = " << C0 << ";" ;
  str << "\n";
  //---------------------------------
  str << "T = "  << T  << ";" ;
  str << "D = "  << D  << ";" ;
  str << "Ru = " << Ru << ";" ;
  str << "Cd = " << Cd << ";" ;
  str << "\n";
  //---------------------------------
  str << "ks = " << ks << ";" ;
  str << "v = "  << v  << ";" ;
  str << "realdE="<<realdE<< ";";
  str << "realIp1="<<realIp1<< ";";

  str << "\n";

  return str.str();
}

string outputData::getData()
{
  std::stringstream str;
  str << "-------------------------OutputData of Ohmic Drop-------------------------" <<"\n";
  str << "Lambda = " << Lambda << ", " ;
  str << "Rau = "    << Rau    << ", " ;
  str << "Gamma = "  << Gamma  << ", " ;
  str << "Dtau = "   << Dtau   << ", " ;
  str << "\n";
  //-----------------------
  str << "IMax = " << IMax << ", " ;
  str << "IMin = " << IMin << ", " ;
  str << "\n";
  str << "CsiP = " << CsiP << ", " ;
  str << "CsiM = " << CsiM << ", " ;
  str << "\n";
  str << "PsiP = " << PsiP << ", " ;
  str << "PsiM = " << PsiM << ", " ;
  str << "\n";
  //--------------
  str << "Ep1 = " << Ep1 << ", " ;
  str << "Ip1 = " << Ip1 << ", " ;
  str << "\n";
  str << "Ep2 = " << Ep2 << ", " ;
  str << "Ip2 = " << Ip2 << ", " ;
  str << "\n";
  DEp=Ep1-Ep2;
  str << "DEp = " << DEp << ", " ;
  str << "ks= "   << ks  << ", " ;
  str << "S =" << S <<", ";
  str << "\n";
  //-----------------
  str << "---------------------end of OutputData of Ohmic Drop----------------------"<< endl;
  return str.str();
}
string outputData::format_GetData()
{

  std::stringstream str;

  str.setf(ios::scientific);
  str.precision(9);

  str << Lambda <<"," ;
  str << Rau <<",";
  str << Gamma << ",";
  str << Dtau <<",";
  str << IMax << ",";
  str << IMin << ",";
  str << CsiP << ",";
  str << CsiM << ",";
  str << PsiP << ",";
  str << PsiM << ",";
  str << Ep1 << ",";
  str << Ep2 << ",";
  str << Ip1 << ",";
  str << Ip2 << ",";
  str << DEp << ",";
  str << ks <<  ",";
  str << S	<< ",";

  return str.str();
}




string outputData::getPlotData()
{
  std::stringstream str;
  str.setf(ios::scientific);
  str.precision(9);
  for(int i=0;i<400;i++)
    if((this->plotdata[0][i] !=0) && this->plotdata[0][i] !=0 )
      str << this->plotdata[0][i] << "," <<this->plotdata[1][i] <<"\n";
  return str.str();
}


string outputData::getTitle()
{
  std::stringstream str;
  str << "Lambda," << "Rau," << "Gamma," << "Dtau," ;
  str << "IMax," << "IMin," << "CsiP," << "CsiM," << "PsiP," << "PsiM,";
  str << "Ep1," << "Ep2," << "Ip1," << "Ip2," << "DEp," ;
  str << "ks," << "S";
  return str.str();
}

