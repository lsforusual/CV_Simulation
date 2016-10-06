#include<iostream>
#include<vector>
#include<string>
#include<inputdata.h>
#include<outputdata.h>
#include<fstream>


using namespace::std;

extern int rows;

int getRows(char* filename)
{
  int Rows=0;
  ifstream file(filename);
  if(file.is_open())
    for(char str[256];!file.eof();Rows++)
      file.getline(str,256);
  else
    {
      cout << "Can not opening file !!!\n";
      getchar();
      exit(1);
    }
  file.close();
  return Rows;

}

void writeDataPlot(char* filename, outputData odata)
{
  //std::stringstream str;
  //str << odata.getPlotData();
  ofstream file(filename);
  file << odata.getPlotData() <<endl;
  file.close();
}

void dataRW(char* filename,const char* o,inputData* idata)
{
  int i=0;
  if(*o=='r')
    {
      ifstream file(filename);

      if(!file.is_open())
        {cout << "Can not opening file !!!\n";getchar();exit(1);}
      while(file.peek() != EOF)
        {
          char str[256];

          file.getline(str,256);
          idata[i].setData(str);
          //return idata;
          i++;
        }
      file.close();

    }
}
  //---------------------
void dataRW(char* filename, const char* o, outputData *odata)
  {
      fstream file;

      if (*o=='w') {
          file.open(filename,ios::in);
          if(!file.is_open())
            {
              file.close();
              file.open(filename,ios::out);
              file << (*odata).getTitle() <<endl;
            }
          else
            {

              file.close();
              cout << "File Exist!! Will save data in append mode!!" <<endl;
              file.open(filename,ios::app);
            }
        }else if (*o=='a'){
          file.open(filename,ios::app);
        }


      if(file.is_open())
        {
          //file << odata.getData() <<endl;
          int i=0;
          while(i<rows)
            {
              file << odata[i++].format_GetData() <<endl;
            }
        }
      file.close();


}
