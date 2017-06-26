#include <outputdata.h>
#include <inputdata.h>
#include <cmath>
#include <iostream>

extern outputData dataScan(inputData idata);

long double fitS(inputData idata, outputData odata)
{
	int i=0;
	double step=0.1;
	odata = dataScan(idata);
	long double _s=idata.S;
	long double last[3] ;
	long double realIp1=idata.realIp1;

	//fit for S
	idata.v = 0.1;

	//flag
	//odata = dataScan(idata);
	last[0] = 0.0;

	idata.S *=(1.0 + step);
	odata = dataScan(idata);
	last[1] = odata.Ip1;
	idata.S = _s;

	idata.S *=(1.0 - step);
	odata = dataScan(idata);
	last[2] = odata.Ip1;
	idata.S = _s;

	odata = dataScan(idata);

	while (abs(realIp1-odata.Ip1) > 0.000001 && i<=1000)
	  {
	    if (abs(realIp1-last[1]) < abs(realIp1-last[2]))
	      {
		last[0] = abs(realIp1-odata.Ip1);
		idata.S *=(1.0 + step);
		odata = dataScan(idata);
		if(last[0] < abs(realIp1-odata.Ip1))
		  {
		    idata.S /=(1.0 + step);
		    step /=10;
		  }else if (last[0] == abs(realIp1-odata.Ip1)){
		    break;
		  }

	      }else{
		last[0] = abs(realIp1-odata.Ip1);
		idata.S /=(1.0 + step);
		odata = dataScan(idata);
		if(last[0] < abs(realIp1-odata.Ip1))
		  {
		    idata.S *=(1.0 + step);
		    step /=10;
		  }else if (last[0] == abs(realIp1-odata.Ip1)){
		    break;
		  }
	      }
	    cout << "No_" << i << "interotor the err of dIp1 is :" << abs(realIp1-odata.Ip1) << " the S is " << idata.S*10000 << endl;
	    i++;
	  }
	return 	idata.S;
;
}

long double fitKs(inputData idata, outputData odata)
{
	int i=0;
	double step=0.1;
	odata = dataScan(idata);
	long double realdE=idata.realdE;
	long double _k=idata.ks;
	
	long double last[3] ;

	//idata.S = fitS(idata,odata);

	//fit for Ep
	last[0] = 0.0;

	idata.ks *=(1.0 + step);
	odata = dataScan(idata);
	last[1] = odata.DEp;
	idata.ks = _k;

	idata.ks *=(1.0 - step);
	odata = dataScan(idata);
	last[2] = odata.DEp;
	idata.ks = _k;

	odata = dataScan(idata);

	while (abs(realdE-odata.DEp) > 0.000001 && i<=1000)
	{
		if (abs(realdE-last[1]) < abs(realdE-last[2]))
		{
			last[0] = abs(realdE-odata.DEp);
			idata.ks *=(1.0 + step);
			odata = dataScan(idata);
			if(last[0] < abs(realdE-odata.DEp))
			{
				idata.ks /=(1.0 + step);
				step /=10;
			}else if (last[0] == abs(realdE-odata.DEp)){ 
				break;
			}

		}else{
			last[0] = abs(realdE-odata.DEp);
			idata.ks /=(1.0 + step);
			odata = dataScan(idata);
			if(last[0] < abs(realdE-odata.DEp))
			{
				idata.ks *=(1.0 + step);
				step /=10;
			}else if (last[0] == abs(realdE-odata.DEp)){ 
			    cout << "No_" << i << "interotor the err of delta_dE is :" << abs(realdE-odata.DEp) << " the ks is " << idata.ks*100.0 << endl;
				break;
			}
		};

		cout << "No_" << i << "interotor the err of dE is :" << abs(realdE-odata.DEp) << " the ks is " << idata.ks*100.0 << endl;
		i++;
	}

	return idata.ks;

}

int fitFunc(inputData idata, outputData& odata)
{
	idata.S = fitS(idata,odata);
	idata.ks= fitKs(idata,odata);
	odata = dataScan(idata);
	odata.S = idata.S * 10000.0;
	odata.ks = idata.ks * 100.0;
	return 1;
}
