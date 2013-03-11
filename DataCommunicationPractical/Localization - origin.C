
/*******************************************************************************************************************

TO CALCULATE THE FIELD SCATTERED FROM AN UNDULATING SURFACE WHICH HAS BEEN EXCITED
BY AN INFINITE LINE SOURCE. SEE PAGE 680 -BALANIS.

***********************************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
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

double	   Mod(complex);
double     R_source_p(double, double, int);
double     R_source_obs(double, double, double, int n);
double	   R_p_q(int,int);
double     R_surf_obs(double, int m, int n);

const int Groupsize=80;

const double  Epsilon_0=8.854e-12,Mu_0=PI*(4.0e-7), 
	Gross_Step=50.0,c=(1.0/sqrt(Mu_0*Epsilon_0)),
	f=435e6,Lambda=c/f, Omega=2.0*PI*f, Beta_0=Omega*(sqrt(Mu_0*Epsilon_0)),
	Cconst=Omega*Mu_0/4.0;
double  X[385],Y[385],Mx[1110], My[1110], Seg_Length=Lambda/4.0; 

const complex j=complex(0.0,1.0);

void localization ()
{
	FILE *fp1; 

	ofstream coutput1;
	ofstream coutput2;
	ofstream coutput3;
	ifstream cinput1;
	ifstream cinput2;

	int Terr_Length, m, n, z, w, X_Index, Y_Index, Tuple_No, No_Grps, Start, Grp_No, Group_No_Receiver_1, Group_No_Receiver_2, Group_No_Source, No_Residual_Surfaces=7, Window_Size=2;

	complex Zself, Hself, Sigma, Sum, Const;

	Terr_Length=7950;
	No_Grps=(int)((((double)(Terr_Length)))/(Groupsize*Seg_Length));

	Zself=complex(Seg_Length, Seg_Length*((-2.0/PI)*(log((1.781*Beta_0*Seg_Length)/4.0)-1.0)));
	Hself=complex(1.0,((-2.0/PI)*(log((1.781*Beta_0*Seg_Length)/4.0)-1.0)));

	vector<complex> J(No_Grps), Field(No_Grps), Primary_Source_Field(No_Grps), Secondary_Source_Field_1(No_Grps), Secondary_Source_Field_2(No_Grps) ;

	double Observed_Primary_Field_Difference_1_2, Observed_Primary_Field_Difference_2_1,Total_SS, Minimum[3], Temp;
	vector<double> Residual(No_Grps), Difference_Secondary_Fields_1_2(No_Grps);

	double Xsource = 0.0, Ysource=16.4, 

		Xreceiver_A=1500.0, Yreceiver_A=7.4, Mreceiver_A=-100.573, 
		Xreceiver_B=5500.0, Yreceiver_B=46.4, Mreceiver_B=-107.0, 

		Xreceiver_C=5250.0, Yreceiver_C=41.0,Mreceiver_C=-112.97,
		Xreceiver_D=2500.0, Yreceiver_D=16.4, Mreceiver_D=-97.2,

		Xreceiver_E=4500.0, Yreceiver_E=34.4, Mreceiver_E=-120.252,
		Xreceiver_F=7000.0, Yreceiver_F=48.4, Mreceiver_F=-122.224,

		Xreceiver_G=2000.0, Yreceiver_G=7.4, Mreceiver_G=-102.51,
		Xreceiver_H=6000.0, Yreceiver_H=53.5,Mreceiver_H=-122.41,

		Xreceiver_J=6500.0, Yreceiver_J=42.4, Mreceiver_J=-140.62,
		Xreceiver_K=5250.0, Yreceiver_K=41.0,Mreceiver_K=-112.97,

		Xreceiver_L=3000.0, Yreceiver_L=26.4, Mreceiver_L=-99.06,
		Xreceiver_M=5750.0, Yreceiver_M=58.6,Mreceiver_M=-109.4,

		Xreceiver_Err1=3800.0, Yreceiver_Err1=36.4, Mreceiver_Err1=-115.01,
		Xreceiver_Err2=6300.0, Yreceiver_Err2=40.8,Mreceiver_Err2=-132.54,

		Obs;

	vector<vector<double> > Receiver_Tuple_Array(No_Residual_Surfaces, vector<double> (7));

	// Receiver Pairs are: A,D  B,E  C,F

	Receiver_Tuple_Array[0][1]=Xreceiver_A;
	Receiver_Tuple_Array[0][2]=Yreceiver_A;
	Receiver_Tuple_Array[0][3]=Mreceiver_A;
	Receiver_Tuple_Array[0][4]=Xreceiver_B;
	Receiver_Tuple_Array[0][5]=Yreceiver_B;
	Receiver_Tuple_Array[0][6]=Mreceiver_B;

	Receiver_Tuple_Array[1][1]=Xreceiver_C;
	Receiver_Tuple_Array[1][2]=Yreceiver_C;
	Receiver_Tuple_Array[1][3]=Mreceiver_C;
	Receiver_Tuple_Array[1][4]=Xreceiver_D;
	Receiver_Tuple_Array[1][5]=Yreceiver_D;
	Receiver_Tuple_Array[1][6]=Mreceiver_D;

	Receiver_Tuple_Array[2][1]=Xreceiver_E;
	Receiver_Tuple_Array[2][2]=Yreceiver_E;
	Receiver_Tuple_Array[2][3]=Mreceiver_E;
	Receiver_Tuple_Array[2][4]=Xreceiver_F;
	Receiver_Tuple_Array[2][5]=Yreceiver_F;
	Receiver_Tuple_Array[2][6]=Mreceiver_F;

	Receiver_Tuple_Array[3][1]=Xreceiver_G;
	Receiver_Tuple_Array[3][2]=Yreceiver_G;
	Receiver_Tuple_Array[3][3]=Mreceiver_G;
	Receiver_Tuple_Array[3][4]=Xreceiver_H;
	Receiver_Tuple_Array[3][5]=Yreceiver_H;
	Receiver_Tuple_Array[3][6]=Mreceiver_H;

	Receiver_Tuple_Array[4][1]=Xreceiver_J;
	Receiver_Tuple_Array[4][2]=Yreceiver_J;
	Receiver_Tuple_Array[4][3]=Mreceiver_J;
	Receiver_Tuple_Array[4][4]=Xreceiver_K;
	Receiver_Tuple_Array[4][5]=Yreceiver_K;
	Receiver_Tuple_Array[4][6]=Mreceiver_K;

	Receiver_Tuple_Array[5][1]=Xreceiver_L;
	Receiver_Tuple_Array[5][2]=Yreceiver_L;
	Receiver_Tuple_Array[5][3]=Mreceiver_L;
	Receiver_Tuple_Array[5][4]=Xreceiver_M;
	Receiver_Tuple_Array[5][5]=Yreceiver_M;
	Receiver_Tuple_Array[5][6]=Mreceiver_M;

	Receiver_Tuple_Array[6][1]=Xreceiver_Err1;
	Receiver_Tuple_Array[6][2]=Yreceiver_Err1;
	Receiver_Tuple_Array[6][3]=Mreceiver_Err1;
	Receiver_Tuple_Array[6][4]=Xreceiver_Err2;
	Receiver_Tuple_Array[6][5]=Yreceiver_Err2;
	Receiver_Tuple_Array[6][6]=Mreceiver_Err2;

	double Xreceiver1, Yreceiver1, Mreceiver1, Xreceiver2,Yreceiver2, Mreceiver2;

	vector< vector< vector<double> > > Residual_Signal_Strength_Surfaces(No_Residual_Surfaces, vector<vector<double> >(No_Grps, vector<double> (50)) );
	vector< vector<double> > Averaged_Sum_of_Residual_Surfaces(No_Grps, vector<double> (50)), Sum_of_Residual_Surfaces(No_Grps, vector<double> (50));

	/*****************************************************************************************

	SETUP

	****************************************************************************************/

	Sum=complex(0.0,0.0);

	for(n=1;n<=(Groupsize/2-1);n++)
	{
		Sum+=abs(Exp(Beta_0*n*Seg_Length)*H02((Groupsize/2-n)*Seg_Length));
	}

	Const=((2.0*(Sum)+abs(Hself)));

	fp1=fopen("hadsund.dhm","r");
	for (n=0;n<160;n++)
	{
		fscanf(fp1, "%lf %lf", &X[n],&Y[n]);
	}
	fclose(fp1);

	Start=0;

	// For Each Tuple (Receiver Pair
	for(Tuple_No=0;(Tuple_No<=No_Residual_Surfaces-1)&&(Tuple_No !=6);Tuple_No++)
	{
		Xreceiver1=Receiver_Tuple_Array[Tuple_No][1];
		Yreceiver1=Receiver_Tuple_Array[Tuple_No][2];
		Mreceiver1=Receiver_Tuple_Array[Tuple_No][3];

		Xreceiver2=Receiver_Tuple_Array[Tuple_No][4];
		Yreceiver2=Receiver_Tuple_Array[Tuple_No][5];
		Mreceiver2=Receiver_Tuple_Array[Tuple_No][6];

		// For Each Observation Point (Y-value)
		Obs=0.0;
		for(Y_Index=0; Y_Index<=49; Y_Index++)
		{
			Obs+=1.0;

			/***********************************************************************************************************************************************

			CALCULATIONS

			**************************************************************************************************************************************/

			Group_No_Receiver_1=(int)(Xreceiver1/((double)Groupsize*Seg_Length));

			Group_No_Receiver_2=(int)(Xreceiver2/((double)Groupsize*Seg_Length));

			/********************************************************************************************************************************************

			RECEIVER 1

			**************************************************************************************************************************************/

			/*****************************
			FWD Receiver 1
			******************************/

			for (n=Group_No_Receiver_1;n<=No_Grps-1;n++)
			{
				Sigma=complex(0.0,0.0);
				Secondary_Source_Field_1[n]=Ei_Rad(R_source_obs(Obs, Xreceiver1, Yreceiver1, n));

				for(m=Group_No_Receiver_1+Start;m<n;m++)
				{
					Sigma+=J[m]*Const*H02(R_p_q(m,n));
					Secondary_Source_Field_1[n]-=J[m]*Const*H02(R_surf_obs(Obs, m,n));
				}
				//Sigma=complex(0.0,0.0);
				J[n]=(Ei_Rad(R_source_p(Xreceiver1, Yreceiver1, n))-(Sigma))/(Hself);
			}

			/***************************
			BK Receiver 1
			****************************/

			for (n=Group_No_Receiver_1-Start;n>=0;n--)
			{
				Sigma=complex(0.0,0.0);
				Secondary_Source_Field_1[n]=Ei_Rad(R_source_obs(Obs, Xreceiver1, Yreceiver1,n));

				for(m=Group_No_Receiver_1;m>n;m--)
				{
					Sigma+=J[m]*Const*H02(R_p_q(m,n));
					Secondary_Source_Field_1[n]-=J[m]*Const*H02(R_surf_obs(Obs, m,n));
				}
				//Sigma=complex(0.0,0.0);
				J[n]=(Ei_Rad(R_source_p(Xreceiver1, Yreceiver1, n))-(Sigma))/(Hself);
			}

			coutput1.open("E_Receiver_1.dat");

			for(n=0;n<=No_Grps-1;n++)
			{
				coutput1 << x(n) << "  " << 20*log10(abs(Secondary_Source_Field_1[n])/*-sqrt(R_source_obs(Obs,Xreceiver1, Yreceiver1,n+1))*/)<< endl;
			}

			coutput1.close();


			/********************************************************************************************************************************************

			RECEIVER 2

			**************************************************************************************************************************************/


			/*****************************
			FWD Receiver 2
			******************************/

			for (n=Group_No_Receiver_2;n<=No_Grps-1;n++)
			{
				Sigma=complex(0.0,0.0);
				Secondary_Source_Field_2[n]=Ei_Rad(R_source_obs(Obs, Xreceiver2, Yreceiver2, n));

				for(m=Group_No_Receiver_2+Start;m<n;m++)
				{
					Sigma+=J[m]*Const*H02(R_p_q(m,n));
					Secondary_Source_Field_2[n]-=J[m]*Const*H02(R_surf_obs(Obs,m,n));
				}
				//Sigma=complex(0.0,0.0);
				J[n]=(Ei_Rad(R_source_p(Xreceiver2, Yreceiver2, n))-(Sigma))/(Hself);
			}

			/***************************
			BK Receiver 2
			****************************/

			for (n=Group_No_Receiver_2;n>=0;n--)
			{
				Sigma=complex(0.0,0.0);
				Secondary_Source_Field_2[n]=Ei_Rad(R_source_obs(Obs,Xreceiver2, Yreceiver2,n));

				for(m=Group_No_Receiver_2-Start;m>n;m--)
				{
					Sigma+=J[m]*Const*H02(R_p_q(m,n));
					Secondary_Source_Field_2[n]-=J[m]*Const*H02(R_surf_obs(Obs,m,n));
				}
				//Sigma=complex(0.0,0.0);
				J[n]=(Ei_Rad(R_source_p(Xreceiver2, Yreceiver2, n))-(Sigma))/(Hself);
			}

			coutput1.open("E_Receiver_2.dat");

			for(n=0;n<=No_Grps-1;n++)
			{

				coutput1 << x(n) << "  " << 20*log10(abs(Secondary_Source_Field_2[n])/*-sqrt(R_source_obs(Obs,Xreceiver2, Yreceiver2,n+1))*/) << endl;
			}

			coutput1.close();

			/***********************************************************************************************************************************

			CALCULATIONS

			***********************************************************************************************************************************/

			Observed_Primary_Field_Difference_1_2 = abs(abs(Mreceiver1)-abs(Mreceiver2));

			cout << "Observed Primary Field Difference  " << Observed_Primary_Field_Difference_1_2 << endl;

			for (n=0;n<=No_Grps-1;n++)
			{
				Difference_Secondary_Fields_1_2[n]
					= abs(abs(20*log10(abs(Secondary_Source_Field_1[n])))-abs(20*log10(abs(Secondary_Source_Field_2[n]))));

				Residual[n]=abs((Difference_Secondary_Fields_1_2[n])-(Observed_Primary_Field_Difference_1_2));
			}

			coutput3.open("Residual.dat");

			for(n=0;n<=No_Grps-1;n++)
			{
				coutput3 << x(n) << "  " << Residual[n] << endl;
			}

			coutput3.close();

			// Input Residual to  Residual SS_Surface

			for(X_Index=0;X_Index<=No_Grps-1;X_Index++)
			{
				Residual_Signal_Strength_Surfaces[Tuple_No][X_Index][Y_Index]=Residual[X_Index];
			}

			coutput3.open("Residual.dat");

			for(n=0;n<=No_Grps-1;n++)
			{
				coutput3 << x(n) << "  " << Residual[n] << endl;
			}

			coutput3.close();

		} // End Y_Index (Obs)

	} // End Tuple_No (Surface_No)

	// Perform Sum of SS Surfaces
	for(Tuple_No=0;Tuple_No<=No_Residual_Surfaces-1;Tuple_No++)
	{
		for(X_Index=0;X_Index<=No_Grps-1;X_Index++)
		{
			for(Y_Index=0;Y_Index<=49;Y_Index++)
			{
				Sum_of_Residual_Surfaces[X_Index][Y_Index]+=(Residual_Signal_Strength_Surfaces[Tuple_No][X_Index][Y_Index]);
			}
		}
	}

	// Generate Moving Average of Sum of Residual SS Surface
	Averaged_Sum_of_Residual_Surfaces[0][0]=Sum_of_Residual_Surfaces[0][0];

	// X Domains
	for(Y_Index=0;Y_Index<=49;Y_Index++)
	{
		for(X_Index=1;X_Index<=No_Grps-1;X_Index++)
		{
			Averaged_Sum_of_Residual_Surfaces[X_Index][Y_Index]
				=(Sum_of_Residual_Surfaces[X_Index][Y_Index]+Sum_of_Residual_Surfaces[X_Index-1][Y_Index])/2.0;
		}
	}

	// Y Domains
	for(X_Index=0;X_Index<=No_Grps-1;X_Index++)
	{
		for(Y_Index=1;Y_Index<=49;Y_Index++)
		{
			Averaged_Sum_of_Residual_Surfaces[X_Index][Y_Index]
				=(Sum_of_Residual_Surfaces[X_Index][Y_Index]+Sum_of_Residual_Surfaces[X_Index][Y_Index-1])/2.0;
		}
	}

	coutput1.open("Averaged_Sum_of_Residual_Surface.dat");
	for(X_Index=0;X_Index<=No_Grps-1;X_Index++)
	{
		for(Y_Index=0;Y_Index<=49;Y_Index++)
		{
			coutput1 << x(X_Index) << "  " << Y_Index+1 << "  " << Averaged_Sum_of_Residual_Surfaces[X_Index][Y_Index] << endl;
		}
		coutput1 << endl;
	}
	coutput1.close();

	// Obtain minimum of sum of surfaces
	coutput1.open("Location_of_Minimum.dat");

	Minimum[0]=0;
	Minimum[1]=0;
	Minimum[2]=100;

	for(X_Index=0;X_Index<=No_Grps-1;X_Index++)
	{
		for(Y_Index=0;Y_Index<=49;Y_Index++)
		{
			if (Averaged_Sum_of_Residual_Surfaces[X_Index][Y_Index]< Minimum[2])
			{
				Minimum[0]=X_Index;
				Minimum[1]=Y_Index;
				Minimum[2]=Averaged_Sum_of_Residual_Surfaces[X_Index][Y_Index];
			}

			coutput1 << x(X_Index) << "  " << Y_Index+1 << "  " <<  Averaged_Sum_of_Residual_Surfaces[X_Index][Y_Index] << endl;

			//coutput1 << x(Minimum[0]) << "  " << Minimum[1]+1 << "  " << Minimum[2] << endl;
		}
	}
	coutput1.close();

	// OVERALL RESULT
	cout << "X_Ordinate =" << x(Minimum[0]) << endl;
	cout << "Y_Ordinate =" << Minimum[1]+1 << endl;
	cout << "Minimum Value =" << Minimum[2] << endl;
	cout << "Error =" << sqrt( ((Xsource-x(Minimum[0]))*(Xsource-x(Minimum[0])))-((Ysource-Minimum[1]+1)*(Ysource-Minimum[1]+1))  ) << endl;

} /* Main*/    

double Mod(complex Z)
{
	return(sqrt((Z.real()*Z.real())+(Z.imag()*Z.imag())));
}

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

double R_source_p(double xsrc, double ysrc, int p)
{
	return(sqrt(((xsrc-x(p))*(xsrc-x(p)))+
		((ysrc-y(p,Y))*(ysrc-y(p,Y)))));
}

double R_source_obs(double obs, double xsrc, double ysrc, int p)
{
	return(sqrt(((xsrc-x(p))*(xsrc-x(p)))+
		((ysrc-y(p,Y)-obs)*(ysrc-y(p,Y)-obs))));
}

double R_p_q(int p, int q)
{
	return(sqrt(((x(q)-x(p))*(x(q)-x(p)))
		+((y(q,Y)-y(p,Y))*(y(q,Y)-y(p,Y)))));
}

double R_surf_obs(double obs, int p, int q)
{
	return(sqrt(((x(p)-x(q))*(x(p)-x(q)))+
		(((y(q,Y)+obs)-y(p,Y))*((y(q,Y)+obs)-y(p,Y)))));
}

int main()
{
	localization();
	system("pause");
	return 0;
}

//Copright Eamonn O Nuallain, Department of Computer Science, Trinity College Dublin 2002/3
