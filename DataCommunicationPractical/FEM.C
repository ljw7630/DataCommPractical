
/*******************************************************************************************************************

TO CALCULATE THE FIELD SCATTERED FROM AN UNDULATING SURFACE WHICH HAS BEEN EXCITED
BY AN INFINITE LINE SOURCE. SEE PAGE 680 - BALANIS.

***********************************************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <complex>
using namespace std;

#define EXP	    2.718281828
#define PI          3.14159265358979323846
#define CONJ -1
#define NO_CONJ 1
#define TRANS 1

#define complex std::complex<double>

#define SWAP(a,b) tempr=complex(a);(a)=complex(b);(b)=tempr

typedef int boolean;
#define TRUE 1
#define FALSE 0

complex    Exp(double);
complex    Exp(complex);

complex    G(double, boolean);
complex    H02(double);
complex    H02_FF(double);

complex    cplx_log(complex);

complex	   Ei_Rad(double);

double	   x(int);
double	   y(int,double*);

double     R_source_p(int);
double     R_source_obs(int n);
double	   R_p_q(int,int);
double     R_surf_obs(int m, int n);

const int Groupsize=13;

const double  Epsilon_0=8.854e-12,Mu_0=PI*(4.0e-7), 
	Gross_Step=10.0,c=(1.0/sqrt(Mu_0*Epsilon_0)),
	f=970e6,Lambda=c/f, Omega=2.0*PI*f, Beta_0=Omega*(sqrt(Mu_0*Epsilon_0)), 
	Xsource = 0.0, Ysource=442.0;
double  X[385],Y[385];
const double Seg_Length=Lambda/4.0; 

const complex j=complex(0.0,1.0);

int fem()
{
	FILE *fp; 
	ofstream coutput1;
	ofstream coutput2;
	ifstream cinput;

	int m, n, Start;

	complex Zself, Hself, Sigma, Sum, Const;

	const int Terr_Length=700;
	const int No_Grps=(int)((((double)(Terr_Length)))/(Groupsize*Seg_Length));

	Zself=complex(Seg_Length, Seg_Length*((-2.0/PI)*(log((1.781*Beta_0*Seg_Length)/4.0)-1.0)));
	Hself=complex(1.0,((-2.0/PI)*(log((1.781*Beta_0*Seg_Length)/4.0)-1.0)));

	vector<complex> J(No_Grps+1), Field(No_Grps+1);

	/*****************************************************************************************

	SETUP

	****************************************************************************************/

	Start=90;

	Sum=complex(0.0,0.0);

	for(n=1;n<=(Groupsize/2-1);n++)
	{
		Sum+=abs(Exp(Beta_0*n*Seg_Length)*H02((Groupsize/2-n)*Seg_Length));
	}

	Const=((2.0*(Sum)+abs(Hself)));
	cout << Const << endl;

	fp=fopen("X.04","r");
	for (n=0;n<384;n++)
	{
		fscanf(fp, "%lf %lf", &X[n],&Y[n]);
	}
	fclose(fp);
	/***********************************************************************************************************************************

	CALCULATIONS
	***********************************************************************************************************************************/

	coutput1.open("E.dat");
	coutput2.open("J.dat");

	for (n=0;n<=No_Grps;n++)
	{
		Sigma=complex(0.0,0.0);
		Field[n]=Ei_Rad(R_source_obs(n));

		for(m=Start;m<n;m++)
		{
			Sigma+=J[m]*Const*H02(R_p_q(m,n));
			Field[n]-=J[m]*Const*H02(R_surf_obs(m,n));
		}

		J[n]=(Ei_Rad(R_source_p(n))-(Sigma))/(Hself);

		coutput2 << x(n) << "  " << abs(J[n]) << endl;
		coutput1 << x(n) << "  " << 20*log10(abs(Field[n])/sqrt(R_source_obs(n+1))) << endl;
	}

	coutput1.close();
	coutput2.close();
	
} /* Main*/    

complex Exp(double d)
{
	return(complex(cos(d),-sin(d)));
}

complex Exp(complex Z)
{
	return(exp(Z.imag())*Exp(Z.real()));
}

complex Ei_Rad(double Arg)
{
	return(complex(-j0(Beta_0*Arg),y0(Beta_0*Arg)));
}

complex H02(double Arg)
{
	return(complex(j0(Beta_0*Arg),-y0(Beta_0*Arg)));
}

complex H02_FF(double Arg)
{
	return(sqrt(2.0/(PI*Beta_0*Arg))*Exp((Beta_0*Arg)-(PI/4.0)));
}

complex cplx_log(complex Z)
{
	return(complex(log(abs(Z)),atan2(Z.imag(),Z.real())));
}

double x(int a)
{
	return((double)a*Seg_Length*Groupsize);
}

double y(int a,double Y[])
{
	double Temp, Prop, s;
	int Index; 
	Temp=((double)a*Seg_Length*Groupsize)/Gross_Step;
	Index=(int)Temp;
	Prop=Temp-(double)Index;

	s=Y[Index]+(Prop*(Y[Index+1]-Y[Index]));

	return(s);
}

double R_source_p(int p)
{
	return(sqrt(((Xsource-x(p))*(Xsource-x(p)))+
		((Ysource-y(p,Y))*(Ysource-y(p,Y)))));
}

double R_source_obs(int p)
{
	return(sqrt(((Xsource-x(p))*(Xsource-x(p)))+
		((Ysource-y(p,Y)-2.4)*(Ysource-y(p,Y)-2.4))));
}

double R_p_q(int p, int q)
{
	return(sqrt(((x(q)-x(p))*(x(q)-x(p)))
		+((y(q,Y)-y(p,Y))*(y(q,Y)-y(p,Y)))));
}

double R_surf_obs(int p, int q)
{
	return(sqrt(((x(p)-x(q))*(x(p)-
		x(q)))+(((y(q,Y)+2.4)-y(p,Y))*((y(q,Y)+2.4)-
		y(p,Y)))));
}

int main()
{
	return fem();
}

//Copright Eamonn O Nuallain, Department of Computer Science, Trinity College Dublin 2002/3
