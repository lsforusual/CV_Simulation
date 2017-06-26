#include <stdlib.h>
#include <stdio.h>
#include <SDL/SDL.h>
#include <math.h>

void pause();
void setPixel(SDL_Surface *surface, int x, int y, Uint32 pixel);

int main(int argc, char *argv[])

{
   int Nmax = 4000;
   int I = 0;
   int J = 0;
   int N = 0;
   long double Psif[8100];
   long double PsiTotal[8100];
   long double Rqd[8100];
   long double Pota[8100];
   long double Dtau;
   long double Csi;
   long double Tau = 0;
   long double Tauf = 0;
   long double u;
   long double Pi = 3.14159265358979323846264338327950288419716939937510;
   long double Sum;
   long double Ct1;
   long double Ct2;
   long double Ct3;
   long double Ct4;
   long double Ct5;
   long double Lambda;
   long double Alpha = 0.5;
   long double PsiProv = 0;
   long double CsiPrime;
   long double Theta = 0; // Zero pour lest tests
   long double Gamma = 0; // Zero pour les tests
   long double Rau;
   long double CsiP;
   long double PsiP;
   long double CsiM;
   long double PsiM;
   double Ei;
   double Ef;
   double E0;
   double S;
   double C0;
   double T;
   double D;
   double Ru;
   double Cd;
   double ks;
   double v;
   long double FRT;
   long double FSC0;
   long double DFRT;
   long double F = 96485.3329;
   long double R = 8.3144621;
   long double Ep1;
   long double Ep2;
   long double Ip1;
   long double Ip2;




   FILE* fichier = NULL;

    // initialisation, read data from input.txt


    fichier = fopen("input.txt", "r");

     if (fichier != NULL)
        {
        // On peut lire dans le fichier

            int ret = fscanf(fichier,"%lf", &Ei);
            ret = fscanf(fichier,"%lf", &Ef);
            ret = fscanf(fichier,"%lf", &E0);
            ret = fscanf(fichier,"%lf", &S);
            S = S /10000; //Conversion cm2 en m2
            ret = fscanf(fichier,"%lf", &C0);
            C0 = C0 *1000.0; // Conversion mol.L-1 en mol.m-3
            ret = fscanf(fichier,"%lf", &T);
            ret = fscanf(fichier,"%lf", &D);
            D = D / 10000; // Conversion cm2.s-1 en m2.s-1
            ret = fscanf(fichier,"%lf", &Ru);
            ret = fscanf(fichier,"%lf", &Cd);
            ret = fscanf(fichier,"%lf", &ks);
            ks = ks / 100; // Conversion cm.s-1 en m.s-1
            ret = fscanf(fichier,"%lf", &v);

            if ( ret == 0 )
            {
            // test on input data
                printf ("Problem Input Data, ret = %d \n", ret);
            }
        }
        else

            {
        // test on file
            printf("Impossible to open file input.txt \n Invalid Input File ");
            return EXIT_SUCCESS;

            }

    fclose(fichier);


    FRT = F / R / T;
    FSC0 = F * S * C0;
    DFRT = D * FRT;

    //parameters are defined for oxidation
    u = -(Ei - E0) * FRT;
    CsiM = (Ef - E0) * FRT;
    Lambda = ks / sqrt (v * DFRT);
    Rau = FRT * Ru * FSC0 * sqrt (v *DFRT);
    Gamma = Cd * sqrt (v) / FSC0 /sqrt (DFRT);

    printf (" ************************************************************* \n Ohmic Drop - Slow Charge Transfer - Version 3.2b 09/15/2016 \n " );
    printf ("                     CNRS / Univ. Rennes 1  \n +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ \n\n " );

    printf (" Ei = %f ; Ef = %f ; S = %e ;  C0 = %e  \n", Ei, Ef, S, C0);
    printf (" T = %f ; D = %e ; Ru = %e ;  Cd = %e  \n", T, D, Ru, Cd);
    printf (" ks = %f ; v = %f \n", ks, v);


    Psif[0] = 0;
    Psif[1] = 0;
    Dtau = (CsiM + u) / Nmax; // calculation of the time step.
    Tau = Dtau;
    Csi = -u;
    Theta = Rau * Gamma;
    Ct3 = Theta/Dtau;

    printf (" Lambda = %LF ; Rau = %Lf ; Gamma = %Lf ; Dtau = %Lf \n", Lambda, Rau, Gamma, Dtau);

    for (I = 1 ; I <= 2.0 * Nmax + 1; I++)

        {
            Rqd[I] = (powl ((I + 1.0) , 1.5)) - (2.0 * powl (I, 1.5)) + (powl ((I - 1) , 1.5));
        }

 // Forward Scan


    for (I = 1 ; I <= Nmax ; I++)
        {
            Ct1 = expl (-Alpha * Csi) / Lambda;
            Ct2 = 4.0 / 3.0 * sqrtl (Dtau / Pi) * (1.0 + expl (-Csi));
            Pota[I] = Csi;

            Sum = 0;

            for (N = 1 ; N < I ; N++)

                {
                    Sum = Sum + (Psif[N] * (1.0 + Ct3) - Psif[N - 1] * Ct3) * Rqd[I - N]; //P335末项的求和项
                }

            N = 0;

            Psif[I] = Psif[I - 1];

            do

                {
                    Psif[I] = (Psif[I] + 100.0 * PsiProv) /101.0; // Il faut calmer la convergence //*********使Psif偏移约%1，记作PsiProv，然后用此Psif回代入方程，求解此条件下的Psif
                    PsiProv = Psif[I];
                    CsiPrime = Csi - Rau * Psif[I] - Theta * (1.0 - expl (-Tau/Theta));   //方程(13)
                    Ct1 = expl (-Alpha * CsiPrime) / Lambda;
                    Ct2 = 4.0 / 3.0 * sqrtl (Dtau / Pi) * (1.0 + expl (-CsiPrime)); //方程（11）的变形项

                    Psif[I] = (1.0 - Sum * Ct2 + (Ct1 * Ct3 + Ct2 * Ct3) * Psif[I - 1] ) / (Ct1 + Ct2 + Ct1 * Ct3 + Ct2 * Ct3);

                    N++; //???????????????为什么要在这里进行N++

                }
            while ( fabs ( Psif[I] - PsiProv ) > 0.0001 ); // 控制精度为delt<0.0001


            PsiTotal[I] = Psif[I] + Gamma * (1.0 - expl (-Tau / Theta));  //方程（12）

            if ((I % 200) == 0)

            {

                printf (" I= %d Csi= %Lf Psif= %Lf PsiTotal= %Lf N= %d \n", I, Csi, Psif[I], PsiTotal[I], N);

            }

            Csi = Csi + Dtau;
            Tau = Tau + Dtau;
        }

    Tauf = Tau - Dtau;


// Reverse Scan

    for (I = Nmax + 1 ; I <= 2.0 * Nmax ; I++)
    {
        Ct1 = expl (-Alpha * Csi) / Lambda;
        Ct2 = 4.0 / 3.0 * sqrtl (Dtau / Pi) * (1.0 + expl (-Csi));

        Pota[I] = Csi;


        Sum = 0;
        for (N = 1 ; N < I ; N++)

            {
                Sum = Sum + (Psif[N] * (1.0 + Ct3) - Psif[N - 1] * Ct3) * Rqd[I - N];
            }

        N = 0;

        Psif[I] = Psif[I - 1]; // We start with psi calculated at N-1. This is the guest value.

        if (Theta > 0.0001 ) // Do not need to consider the double layer if Theta is small

            {
                Ct4 = Theta * (2.0 * expl ((Tauf-Tau)/Theta) - expl (-Tau / Theta) - 1.0);
            }

        else

            {
                Ct4 = 0.0;

            }

        Ct5 = Gamma * (2.0 * expl ( (Tauf - Tau) / Theta) - expl (-Tau / Theta) - 1.0);

                do

                {
                    Psif[I] = (Psif[I] + 100.0 * PsiProv) /101.0; // This is required to dump the convergence.
                    PsiProv = Psif[I];
                    CsiPrime = Csi - Rau * Psif[I] -  Ct4;
                    Ct1 = expl (-Alpha * CsiPrime) / Lambda;
                    Ct2 = 4.0 / 3.0 * sqrtl (Dtau / Pi) * (1.0 + expl (-CsiPrime));

                    Psif[I] = (1.0 - Sum * Ct2 + (Ct1 * Ct3 + Ct2 * Ct3) * Psif[I - 1] ) / (Ct1 + Ct2 + Ct1 * Ct3 + Ct2 * Ct3);

                    N++;

                   }

                while ( fabs ( Psif[I] - PsiProv ) > 0.0001 );


        PsiTotal[I] = Psif[I] + Ct5;


        if ((I % 200) == 0)
        {
             printf (" I= %d Csi= %Lf Psif= %Lf PsiTotal= %Lf N= %d \n", I, Csi, Psif[I], PsiTotal[I], N);
        }

        Csi = Csi - Dtau;
        Tau = Tau + Dtau;
    }



// Calculations of the maximum and minimum values, Csip, Psip, etc...

    PsiP = -100;
    CsiP = -500;

    for (N = 200 ; N <  Nmax; N++)
        {
            if ( PsiTotal[N] > PsiP )
                {
                    PsiP = ( PsiTotal[N] + PsiTotal[N - 1] ) / 2.0;
                    CsiP = Pota[N] + Dtau / 2.0;
                    I = N;
                }
        }

    PsiM = 100;
    CsiM = 500;

        for (N = Nmax + 200 ; N < 2 * Nmax; N++)
        {
            if ( PsiTotal[N] < PsiM )
                {
                    PsiM = ( PsiTotal[N] + PsiTotal[N - 1] ) / 2.0;
                    CsiM = Pota[N] + Dtau/2.0;
                    J = N;
                }
        }

    Ep1 = E0 + CsiP / FRT;
    Ep2 = E0 + CsiM / FRT;
    Ip1 = PsiP * FSC0 * sqrt(v * DFRT);
    Ip2 = PsiM * FSC0 * sqrt(v * DFRT);

    printf ("\n Calculations of peak potentials and peak currents :\n");
    printf (" IMax = %d ; CsiP = %Lf ; PsiP = %Lf \n", I, CsiP, PsiP);
    printf (" IMin = %d ; CsiM = %Lf ; PsiM = %Lf \n", J, CsiM, PsiM);
    printf (" Ep1 = %Lf V ; Ip1 = %Le A \n", Ep1, Ip1);
    printf (" Ep2 = %Lf V ; Ip2 = %Le A \n", Ep2, Ip2);


// Results are written in a file. results.txt


    fichier = fopen("results.txt", "w");

    for (N = 1 ; N < 2 * Nmax - 1; N++)

        {
            if ((N % 20) == 0)

                {
                    fprintf (fichier," %Lf ; %Le \n", E0 + Pota[N]/FRT, PsiTotal[N] * FSC0 * sqrt(v * DFRT));
                }

        }

    fclose(fichier);




// The following is to draw the graphics windows with the two voltamograms. This part depends on SDL and could be changed for another library.
    SDL_Surface *ecran = NULL, *rectangle = NULL;

    SDL_Init(SDL_INIT_VIDEO);

    ecran = SDL_SetVideoMode(1024, 720, 32, SDL_HWSURFACE);

    // Allocation de la surface

    rectangle = SDL_CreateRGBSurface(SDL_HWSURFACE, 220, 180, 32, 0, 0, 0, 0);
    SDL_WM_SetCaption("Voltamogram - Total current in Blue - Apparent current (minus capacitive) in Yellow", NULL);
    SDL_FillRect(ecran, NULL, SDL_MapRGB(ecran->format, 157, 188, 188));

    SDL_LockSurface(ecran); // On bloque la surface

    for (N = 1 ; N < 2 * Nmax - 1 ; N++)

    {
        I = ceil ((Pota[N] - Pota[1]) / (Pota[Nmax] - Pota[1]) * 1020 ) ;
        J = ceil (360 - Psif[N] * 500);

       if (I > 0 && I < 1024 && J > 0 && J < 720 )

        {
            setPixel(ecran, I, J, SDL_MapRGB(ecran->format, 255, 255, 0)); // We draw a red pixel
        }

     }

    for (N = 1 ; N < 2 * Nmax - 1; N++)

    {
        I = ceil ((Pota[N] - Pota[1]) / (Pota[Nmax] - Pota[1]) * 1020 ) ;
        J = ceil (360 - PsiTotal[N] * 500);

        if (I > 0 && I < 1024 && J > 0 && J < 720 )

         {
            setPixel(ecran, I, J, SDL_MapRGB(ecran->format, 0, 0, 255)); // We darw a blue pixel
         }

    }

    SDL_UnlockSurface(ecran); // Unlock the surface
    SDL_Flip(ecran); // Update of the screen.

    pause();

    SDL_FreeSurface(rectangle); // We free the surface.
    SDL_Quit();



    return EXIT_SUCCESS;

}



void pause()

{
    int continuer = 1;
    SDL_Event event;

    while (continuer)
    {
        SDL_WaitEvent(&event);
        switch(event.type)

        {
            case SDL_QUIT:
                continuer = 0;

        }
    }
}


void setPixel(SDL_Surface *surface, int x, int y, Uint32 pixel)
{
    int bpp = surface->format->BytesPerPixel;

    Uint8 *p = (Uint8 *)surface->pixels + y * surface->pitch + x * bpp;

    switch(bpp) {
    case 1:
        *p = pixel;
        break;

    case 2:
        *(Uint16 *)p = pixel;
        break;

    case 3:
        if(SDL_BYTEORDER == SDL_BIG_ENDIAN) {
            p[0] = (pixel >> 16) & 0xff;
            p[1] = (pixel >> 8) & 0xff;
            p[2] = pixel & 0xff;
        } else {
            p[0] = pixel & 0xff;
            p[1] = (pixel >> 8) & 0xff;
            p[2] = (pixel >> 16) & 0xff;
        }
        break;

    case 4:
        *(Uint32 *)p = pixel;
        break;
    }
}
