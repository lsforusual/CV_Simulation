#include <outputdata.h>
#include <inputdata.h>
#include <cmath>
#include <iostream>

extern outputData dataScan(inputData idata);

int fitFunc(inputData idata, outputData& odata){
  int i=0;
  double step=0.1;
  odata = dataScan(idata);
  long double realdE=idata.realdE;
  long double k=idata.ks;
  long double last_dep = 0.0 ;
  while(abs(realdE-odata.DEp)>0.0005 && i <=1000 )
    {
      last_dep=abs(realdE-odata.DEp);

      idata.ks /= (1.0 + step) ;
      odata = dataScan(idata);

      if (last_dep < abs(realdE-odata.DEp))
        {

          idata.ks *= (1.0 + step);
          step /= 10;
        }

      i++;
      if(last_dep == abs(realdE-odata.DEp)) break;
      cout << "No " << i << " interator the err is : " << abs(realdE-odata.DEp) << "   the ks is: " <<  idata.ks*100 <<endl;

    }

  if(i>1000)
    {
      idata.ks =k;
      while(abs(realdE-odata.DEp)>0.0005 && i <= 2000)
      {
        last_dep=abs(realdE-odata.DEp);

        idata.ks *= (1.0 + step) ;
        odata = dataScan(idata);

        if(last_dep == abs(realdE-odata.DEp)) break;

        if (last_dep < abs(realdE-odata.DEp))
          {

            idata.ks *= (1.0 + step);
            step /= 10;
          }

        i++;

        cout << i << " : " << abs(realdE-odata.DEp) <<endl;
      }
    }
  if(i>2000) return 0;

  odata.ks = idata.ks*100.0;
  return i;
}
