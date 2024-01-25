#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include "./azar.cc"


// Poronguita
using namespace std;

double elec_crity = 4;

int num = 100;

double perchar = 0.5;
int ite = 270;
int bc = 0;
int Time;
int lostch;
int numberdis;
int cutoff = 80;
double aproxval = 2.*elec_crity; //VALUE FOR THE APROXIMATION FOR INTEGER CHARGE
double resta;
vector < double > pot, charge, elec, tempcharge;
int gaussiana(int);
double Potential(int, vector < double >);
double Electricy(int, vector < double >);
double Green(int, int );
double GreenN(int, int);
double GreenP(int, int);
double Energy(vector < double >);



int main()
{
	randomize();
	ofstream datos("./output/datos.txt");
	ofstream inf("./output/pos_cant.txt");
	datos << num << " " << ite << endl;
	datos.close();


   	//INITIALIZING THE GRID AT 0 POTENCIAL
	for(double i = 0; i < num ; i++)
	{
		int a = 0;
		pot.push_back(a);
    }

	//COPY OF THE GRID
	charge = pot;
	elec = pot;
	//DATA WITH THE CHARGE
	ofstream inicialcarga("./output/initcarg_1D.txt");
	//DATA WITH THE POTENTIAL
	ofstream inicialpotencial("./output/initpot_1D.txt");
	//DATA WITH THE ELECTRIC FIELD
	ofstream inicialelectric("./output/initelec_1D.txt");
	//DATA WITH THE VARIATION OF MOBILE CHARGES BETWEEN DISCHARGES ( I'LL TAKE ONLY ~ 100 STEPS)
	ofstream QM("./output/corriente.txt", ios::app);

   //INITIALIZING THE DISTRIBUTION OF THE INITIAL CHARGE DISTRIBUTION

//RANDOM WITH GAUSSIAN DISTRIBUTION
	charge[gaussiana(num)] = iazar(1,7);

//USING A POST TRANSIENT CONDITION OF THE CHARGES

/*   int run = 0;
	ifstream transiente("fin.txt");
	while(!transiente.eof())
	{
		transiente >> charge[run];
		run++;
		//cout << run << endl;
	}
*/

   	// CALCULATING THE POTENTIAL AND THE ELECTRIC FIELD
	for(int i = 0; i < int(pot.size()); i++)
    {
		pot[i]  = Potential(i,charge);
		elec[i] = sqrt(Electricy(i,charge)*Electricy(i, charge));
	}

   //SAVING THE INITIAL CONDITION TO BOTH FIELDS
	for(int i = 0; i < int(charge.size()); i++)
	{
		inicialelectric  << scientific << setw(15) << setfill(' ') << elec[i] << " " << endl;
		inicialpotencial << scientific << setw(15) << setfill(' ') <<  pot[i] << " " << endl;
		inicialcarga << scientific << setw(15) << setfill(' ') << charge[i] << " " << endl;
	}
	inicialelectric.close();
	inicialpotencial.close();
	inicialcarga.close();

	ofstream enerTot("./output/enerTOT.txt", ios::app);
	enerTot <<  scientific << setw(15) << setfill(' ') << Energy(charge)  << endl;

   	// INITIALIZING THE PROCESS
 	// int porc = cutoff*num/100;
	for(int l = 1; l <= ite; l++)
	{
		int aux1 = gaussiana(num);
		int val = iazar(1,10);
		inf << aux1 << " " << val << endl;
		charge[ aux1 ] = charge[aux1] + val;
		double ini = Energy(charge);
		Time = 0;
		lostch = 0;
		int reset = 1;

	//RECALCULATING THE POTENTIAL AND ELECTRIC FIELD AFTER ADDING THE CHARGE
	for(int i = 0; i < int(pot.size()); i++)
	{
		pot[i]  = Potential(i,charge);
		elec[i] = sqrt(Electricy(i,charge)*Electricy(i, charge));
	     //cout << scientific << setw(10) << setfill(' ') << elec[i] << " " << charge[i] << " " << "parte1" << endl;
	}
	tempcharge = charge;// BACKUP OF THE CHARGE GRID

	int contador = 0;
	int mobileQ = 0;
	////AVALANCHE/////
	while(reset > 0)
	{
		for(int j = 1; j < num; j++)
		{
			bool cut = round(elec[j]/aproxval);
			int tempval = cut*round(perchar*tempcharge[j]);
			tempcharge[j] = tempcharge[j] - tempval;
			tempcharge[j-1] = tempcharge[j-1] + tempval;
			Time = Time + cut;
			contador = contador + tempval;
			mobileQ = mobileQ + tempval;
		}
	    //SECCION PARA LOS ARCHIVOS QUE TENDRAN LA VARIACION DE CARGA Y CAMPO ELECTRICO PARA TODO INDICE N  -- ACA
		//EL CALCULO ES PARA INESTABILIDADES
		stringstream c;
		c << l;
		string caele = "./output/CampoElectrico/elec_"+c.str()+".txt";

		string car = "./output/Carga/carga_"+c.str()+".txt";
	    //DATA WITH THE MOVING CHARGE
		ofstream carga(car.c_str(), ios::app);

	    //DATA WITH THE VARIATION OF POTENTIAL
		ofstream potencial("./output/Potencial/pot_"+c.str()+".txt");

	     //DATA WITH THE VARIATION ELECTRIC FIELD
		ofstream electrico(caele.c_str(), ios::app);

		 //ACA TERMINA EL CALCULO DE INESTABILIDADES

	     //RECALCULATING THE POTENTIAL AND ELECTRIC FIELD FOR RELAXING
		for(int i = 0; i < int(pot.size()); i++)
		{
			elec[i] = sqrt(Electricy(i, tempcharge)*Electricy(i, tempcharge));
			pot[i]  = Potential(i,tempcharge);
			electrico << scientific << setw(15) << setfill(' ') << elec[i] << endl;
			potencial  << scientific << setw(15) << setfill(' ') << pot[i] << endl;
			carga << scientific << setw(15) << setfill(' ') << tempcharge[i]  << endl;
		}
		reset = contador;
		contador = 0;
		electrico.close();
		carga.close();
		}

	QM << mobileQ << endl;
	charge = tempcharge;
	//WRITING THE LOST CHARGE WHEN REACH THE EARTH
	lostch = charge[0];
	ofstream qper("./output/qperdida.txt", ios::app);
	qper <<  scientific << setw(15) << setfill(' ') <<  lostch << endl;
	charge[0] = 0;
	qper.close();

	//TIME DURATION OF THE CHARGE
	ofstream dur("./output/duracion.txt", ios::app);
	dur << Time << endl;
	dur.close();

	//RECALCULATING THE POTENTIAL AND ELECTRIC FIELD AFTER  THE DISCHARGE AND SAVING THE FINAL CHARGE, POTENTIAL, AND ELECTRIC FIELD
	/*	for(int k = 0; k < int(pot.size()); k++)
		{
			pot[k]  = Potential(k,charge);
			elec[k] = sqrt(Electricy(k,charge)*Electricy(k, charge));
			electrico << scientific << setw(15) << setfill(' ') << elec[k]  << endl;
			potencial  << scientific << setw(15) << setfill(' ') << pot[k] <<  endl;
			carga << charge[k]  << endl;
		}
	*/
	double fin = Energy(charge);
	resta = ini - fin;
	ofstream ener("./output/deltaenergia.txt", ios::app);
	ener <<  scientific << setw(15) << setfill(' ') << resta  << endl;
	enerTot <<  scientific << setw(15) << setfill(' ') << fin  << endl;
	ener.close();
	}

	return 0;
}



int gaussiana(int bla)
{
	int acept = 0;
	int azar1;
	while(acept == 0)
	{
        azar1 = iazar(bla-4,bla-1);
		double azar2 = dazar(0,1);
		double gauss = exp(-pow((azar1-bla),2)/4.0);
		if(azar2 < gauss)
			acept = 1;
	}

	return azar1;
}


double Potential(int y, vector < double > prueba)
{
	int tempq = 0;
	double value = 0.;
	for(int i = 0; i < int(prueba.size()); i++ )
	{
		tempq = tempq + prueba[i];
		value  = value + tempq*Green(y,i);
		tempq = 0;
	}
	return value;
}


double Electricy(int y, vector < double > prueba)
{
	int tempq = 0;
	double value = 0;
	for(int i = 0; i < int(prueba.size()); i++)
		{
			tempq = tempq + prueba[i];
			value = value + tempq*((i-y)*GreenN(i,y)*GreenN(i,y)*GreenN(i,y) - (i+y)*GreenP(i,y)*GreenP(i,y)*GreenP(i,y));
			tempq = 0;
		}
	return value;
}


double GreenN(int y,int y0)
{
	if(y0 == y)
	{
		return 0;
	}
	else
	{
		return (1/(sqrt((y-y0)*(y-y0))));
	}
}


double GreenP(int y,int y0)
{
	if(y0 == y)
	{
		return 0;
	}
	else
	{
		return (1/(sqrt((y+y0)*(y+y0))));
	}
}


double Energy(vector < double > prueba)
{
	double value = 0;
	for(int j = 0; j < int(prueba.size()); j++)
	{
		value = value + prueba[j]*Potential(j,prueba);
	}
	return value;
}


double Green(int y, int y0)
{
	if((y0 == y))
	{
		return 0;
	}
	else
	{
		return  1/sqrt((y-y0)*(y-y0)) - 1/sqrt((y+y0)*(y+y0));
	}
}


/*
double ElectricField(vector<double> prueba, int x, int y)
{
	int tempq = 0;
	double value = 0.;

	for(int i = 0; i < int(prueba.size()); i++)
	{
		for(int j = 0; j < int(prueba.size()); j++)
		{
			tempq = tempq + prueba[i][j];
			value = value + tempq*( (1*(i-x)*Green(y, x,i,j)*Green(y,x,i,j)*Green(y,x, i,j)  )  + (1*(j-y)*Green(y, x,i,j)*Green(y,x,i,j)*Green(y,x, i,j)  )  );
			// el  1 es porque en le ecuación hay un a^2 qu es el radio de la esfera-
			tempq = 0;
		}
	}

	return value;
}
*/



