//#include "mainwindow.h"
//#include <QApplication>
#include <cmath>
#include <string>
#include <inputdata.h>
#include <outputdata.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <windows.h>


using namespace::std;

int rows=0;

extern outputData dataScan(inputData idata);
extern void dataRW(char* filename, const char* o,outputData* odata);
extern void dataRW(char* filename, const char* o,inputData* idata);
extern int getRows(char* filename);
extern void writeDataPlot(char* filename, outputData odata);
extern int fitFunc(inputData idata, outputData& odata);



int main(int argc, char *argv[])
{

	//  MainWindow w;
	//  w.show();

	char filename[_MAX_FNAME];
	char file_path[_MAX_PATH];
	char dir[_MAX_DIR];
	char drive[_MAX_DRIVE];


	char* inputFilename;
	char* outputFilename;

	inputFilename="input";

	_splitpath(argv[0],drive,dir,filename,NULL);
	_makepath(file_path,drive,dir,inputFilename,".txt");

	cout << file_path << endl;
	//getchar();

	rows = getRows(file_path);

	inputData *idata  = new inputData[rows];
	outputData *odata = new outputData[rows];

	dataRW(file_path,"r",idata);

	int fit=0;

 	//for(int i=0;i<rows;i++)

	//if(fit == 0) {cout << "can not find the exactly ks"<<endl; system("pause");
	//exit(1);}

	for(int i=0;i<rows;i++)
	{
		char outputFilename[20];

		odata[i] = dataScan(idata[i]);
		sprintf(outputFilename,"Origin_output.csv");
		dataRW(outputFilename,"w",&odata[i]);

		fit = fitFunc(idata[i],odata[i]);

		sprintf(outputFilename,"Fit_output.csv");
		//_makepath(file_path,drive,dir,outputFilename,".csv");
		dataRW(outputFilename,"w",&(odata[i]));

 		sprintf(outputFilename,"output_data_%d_.csv",i);
		//_makepath(file_path,drive,dir,outputFilename,".csv");
		writeDataPlot(outputFilename,odata[i]);
		cout << i <<endl;
	}



	system("pause");
		exit(0);
	//--defination-->
}
