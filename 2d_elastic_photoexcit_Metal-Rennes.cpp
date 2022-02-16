#define _CRT_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <conio.h> 
#include <string.h>

#include <vector>
#include <random>

#include "alglibmisc.h"
#include "gsl/gsl_sort_double.h"
using namespace alglib;
using namespace std;

#define grafic 1
//#undef grafic

#define MHL 1
#undef MHL

#define PHOTO
//#undef PHOTO

constexpr auto n_lat = 50;
constexpr auto n_part = 9900;
constexpr auto n_equ = 39600;
constexpr auto n_max_vec = 6;

constexpr auto radius = 0.2;	// radius trebuie sa ramana 1, e doar ce e la LS
constexpr auto L = 0.6;
constexpr auto rmic = radius;// 0.20;
constexpr auto rmare = 1.1 * radius; // 0.22;

constexpr auto m = 1.0;
constexpr auto Kf_mic_mic = 5.0;
constexpr auto Kf_poly = 1.0;
constexpr auto mu = 0.04;

constexpr auto n_steps = 301;
constexpr auto T_LIM_DWN = 100.0;
constexpr auto T_LIM_UP = 350.0;
constexpr auto delta_T = (T_LIM_UP - T_LIM_DWN) / (n_steps - 1);
constexpr auto T_EXCITATION = 500.0;
constexpr auto N_MAX_STEPS = 100000;

constexpr auto H = 1100.0;		//1100;
constexpr auto S = 5.5;			//7;
constexpr auto E = 400.0;		//200;
constexpr auto ka = 2000.0;		//700;
constexpr auto tau = 100.0;		//50;

constexpr auto CoefTerm = 0.01; //% din diferenta de temperaturi ce se schimba per pas
constexpr auto CoefTermExt = 0.01;

typedef struct
{
public: double x, y, z, raza, theta, k;
}sReadData;

typedef struct
{
public: int vecin, tip_vecin;
}sPosCoef;

sReadData Medium[n_part];
sPosCoef Position_Coef[n_part][n_max_vec];
int   neighbours[n_part];
double sol[n_equ], sol_old[n_equ];
double T[n_part];
double probabilitateHL[n_part], probabilitateLH[n_part], pres[n_part];

double depth = 0.0;

constexpr char fis_particule[500] = "E:\\Stoleriu\\C\\special\\3d\\generare\\2022\\Elastic\\50x50_RektHex_L06_LS.dat"; // HS: r=1.1 L=2

constexpr char fis_solutiiMHL[500] = "E:\\Stoleriu\\C\\special\\3d\\res\\2022\\elastic\\TiOX\\50x50_RektHex_Sol_MHL";
constexpr char fis_volumeMHL[500] = "E:\\Stoleriu\\C\\special\\3d\\res\\2022\\elastic\\TiOX\\50x50_RektHex_Sol_MHL.dat";
constexpr char fis_volumePHOTO[500] = "E:\\Stoleriu\\C\\special\\3d\\res\\2022\\elastic\\TiOX\\50x50_RektHex_Sol_PHOTO0.5_TExcit500_Exo00.dat";

char file[200]      = "E:\\Stoleriu\\C\\special\\3d\\res\\2022\\elastic\\TiOX\\50x50_RektHex_PHOTOViz";
char fileHisto[200] = "E:\\Stoleriu\\C\\special\\3d\\res\\2022\\elastic\\TiOX\\50x50_RektHex_PHOTOHisto";

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

int initializare(void);
void alglib_function_neighbours(void);
int TemperaturiExchange(void);
int Temperaturi(void);
double Suprafata(bool save);
double Finvers(double x);
void doTheHisto(int step);
int Funct_Dopri(double time, double *input, double *deriv);
int Dopri5(double x, double xend, double eps, double hmax, double has, double *sol);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
	FILE *fvol;
	fvol = fopen(fis_volumeMHL, "w");
	fclose(fvol);
	int i, j;

	initializare();

	//////////////////////////////////////////////////////////////////////////
	// PAS DE ADAPTARE LA PRESIUNEA EXTERIOARA
	//////////////////////////////////////////////////////////////////////////

	double timp;
	double t_init = 0.0;
	double step_t = 0.1;
	double eps = 1.0e-5;
	double diference;
	double arie_ini = 0.0;
	double arie = 0.0;
	{
		timp = t_init;


		for (i = 0; i < n_part; i++)
		{
			T[i] = T_LIM_DWN;
		}



		for (i = 0; i < n_part; i++) //Conditii initiale
		{
			sol[4 * i + 0] = Medium[i].x;
			sol[4 * i + 1] = 0.0;
			sol[4 * i + 2] = Medium[i].y;
			sol[4 * i + 3] = 0.0;
			//Medium[i].raza = 1.1 * radius;	//1.1 inseamna HS, 1.0 inseamna LS
		}

		diference = radius;
		j = 0;

		while (((diference) > 1.0e-3) || (j < 10))
		{
			j++;

			for (i = 0; i < n_part; i++)
			{
				sol_old[4 * i + 0] = sol[4 * i + 0];
				sol_old[4 * i + 2] = sol[4 * i + 2];
			}

			Dopri5(timp, timp + step_t, eps, step_t, step_t / 4.0, &sol[0]);

			diference = fabs(sol[0] - sol_old[0]);
			for (i = 0; i < n_part; i++)
			{
				Medium[i].x = sol[4 * i + 0];
				Medium[i].y = sol[4 * i + 2];
				if (fabs(sol[4 * i + 0] - sol_old[4 * i + 0]) > diference) diference = fabs(sol[4 * i + 0] - sol_old[4 * i + 0]);
				if (fabs(sol[4 * i + 2] - sol_old[4 * i + 2]) > diference) diference = fabs(sol[4 * i + 2] - sol_old[4 * i + 2]);
			}
			timp += step_t;
			printf("S  T  E  P :   %d , \t\t diference: %7.5lg \t\t in: %7.3lf ms\n", j, diference, 0.0);
		}

		arie_ini = Suprafata(false);
	}

#ifdef MHL
	//////////////////////////////////////////////////////////////////////////
	// MHL - TO - UP
	//////////////////////////////////////////////////////////////////////////

	int n_H = 0;
	int n_L = n_part;
	int n_L_vechi = n_part;

	timp = t_init;

	for (i = 0; i < n_part; i++)
	{
		T[i] = T_LIM_DWN;
	}

	for (i = 0; i < n_part; i++) //Conditii initiale
	{
		sol[4 * i + 0] = Medium[i].x;
		sol[4 * i + 1] = 0.0;
		sol[4 * i + 2] = Medium[i].y;
		sol[4 * i + 3] = 0.0;
		//Medium[i].raza = 1.1 * radius;
	}

	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_real_distribution<double> rand_dis(0.0, 1.0);  // use with rand_dis(gen)

	int contor_pasi_la_o_temperatura = 0;

	for (int temperature_step = 0; temperature_step < n_steps; temperature_step++)
	{
		for (int i = 0; i < n_part; i++)
		{
			T[i] = T_LIM_DWN + delta_T * temperature_step;
		}

		contor_pasi_la_o_temperatura = 0;
		while ((contor_pasi_la_o_temperatura < 100) && (n_L > 0.025 * n_part))
		{
			contor_pasi_la_o_temperatura++;

			Dopri5(timp, timp + step_t, eps, step_t, step_t / 4.0, &sol[0]);

			for (int i = 0; i < n_part; i++)
			{
				Medium[i].x = sol[4 * i + 0];
				Medium[i].y = sol[4 * i + 2];
			}
			timp += step_t;

			Temperaturi();

			for (int i = 0; i < n_part; i++)
			{
				if ((Medium[i].raza > 1.05 * radius) && (probabilitateHL[i] > rand_dis(gen)))
				{
					Medium[i].raza = rmic;
					n_H--; n_L++;
				}
				else
				{
					if ((Medium[i].raza < 1.05 * radius) && (probabilitateLH[i] > rand_dis(gen)))
					{
						Medium[i].raza = rmare;
						n_H++; n_L--;
					}
				}
			}
		}

		if (fabs((double)(n_L - n_L_vechi)) > ((double)n_part / 100))
		{
#ifdef grafic
			{
				int p, i, j, v1, v2;
				int count = 0;
				int *count_switched;
				double d1;

				count_switched = (int *)calloc(n_part, sizeof(int));

				for (p = 0; p < n_part; p++)
				{
					count_switched[p] = 0;

					for (i = 0; i < neighbours[p]; i++)
					{
						if (Medium[Position_Coef[p][i].vecin].raza > 1.05)
						{
							count_switched[p]++;
						}
					}

					for (i = 0; i < neighbours[p] - 1; i++)
					{
						for (j = i + 1; j < neighbours[p]; j++)
						{
							v1 = Position_Coef[p][i].vecin;
							v2 = Position_Coef[p][j].vecin;

							d1 = sqrt((Medium[v1].x - Medium[v2].x) * (Medium[v1].x - Medium[v2].x) + (Medium[v1].y - Medium[v2].y) * (Medium[v1].y - Medium[v2].y) + (Medium[v1].z - Medium[v2].z) * (Medium[v1].z - Medium[v2].z));

							if ((d1 < 5.0)) //ca sa nu luam doi vecini ambii de pe Ox sau ambii de pe Oy
							{
								count++;
							}
						}
					}
				}


				char fis_save_vis[500];
				sprintf(fis_save_vis, "%s_ucd_%06d.inp", file, (int)timp);

				FILE *fpout;
				fpout = fopen(fis_save_vis, "w");

				fprintf(fpout, "%d %d 1 0 0\n", n_part, count);
				printf("SAVING UCD %d %d\n", n_part, count);

				for (i = 0; i < n_part; i++)
				{
					fprintf(fpout, "%d %f %f %f\n", i + 1, Medium[i].x, Medium[i].y, Medium[i].z);
				}

				count = 0;
				for (p = 0; p < n_part; p++)
				{
					for (i = 0; i < neighbours[p] - 1; i++)
					{
						for (j = i + 1; j < neighbours[p]; j++)
						{
							v1 = Position_Coef[p][i].vecin;
							v2 = Position_Coef[p][j].vecin;

							d1 = sqrt((Medium[v1].x - Medium[v2].x) * (Medium[v1].x - Medium[v2].x) + (Medium[v1].y - Medium[v2].y) * (Medium[v1].y - Medium[v2].y) + (Medium[v1].z - Medium[v2].z) * (Medium[v1].z - Medium[v2].z));

							if ((d1 < 5.0))  //ca sa nu luam doi vecini ambii de pe Ox sau ambii de pe Oy
							{
								count++;
								fprintf(fpout, "%d 1 tri  %d  %d  %d \n", count, p + 1, v1 + 1, v2 + 1);
							}

						}
					}
				}

				fprintf(fpout, "5 1 1 1 1 1\n");
				fprintf(fpout, "raza, nm\n");
				fprintf(fpout, "phase, au\n");
				fprintf(fpout, "PLH, au\n");
				fprintf(fpout, "PHL, au\n");
				fprintf(fpout, "pressure, au\n");

				for (i = 0; i < n_part; i++)
				{
					fprintf(fpout, "%d %lf %lf %f %f %f\n", i + 1, Medium[i].raza, Medium[i].k,
						((probabilitateLH[i] < 1.0) ? (probabilitateLH[i]) : (1.0)),
						((probabilitateHL[i] < 1.0) ? (probabilitateHL[i]) : (1.0)),
						pres[i]);
				}

				fclose(fpout);

				free(count_switched);
			}

#endif
			n_L_vechi = n_L;
		}

		printf("Timp %lf \t\t Temp %lf \t\t HS %d \n", timp, T[0], n_H);

		fvol = fopen(fis_volumeMHL, "a");
		// ************************* SALVARI VOL********************************
		fprintf(fvol, "%lf   %lf   %lf\n", timp, T[0], (double)n_H / n_part);
		// ************************ END SALVARI *****************************
		fclose(fvol);

	}


	//////////////////////////////////////////////////////////////////////////
	// MHL - TO - DOWN
	//////////////////////////////////////////////////////////////////////////

		//n_H = n_part;
		//n_L = 0;
		//n_L_vechi = 0;

		//timp = t_init;

	for (i = 0; i < n_part; i++)
	{
		T[i] = T_LIM_UP;
	}

	for (i = 0; i < n_part; i++) //Conditii initiale
	{
		sol[4 * i + 0] = Medium[i].x;
		sol[4 * i + 1] = 0.0;
		sol[4 * i + 2] = Medium[i].y;
		sol[4 * i + 3] = 0.0;
		//Medium[i].raza = 1.1 * radius;
	}

	for (int temperature_step = 0; temperature_step < n_steps; temperature_step++)
	{
		for (int i = 0; i < n_part; i++)
		{
			T[i] = T_LIM_UP - delta_T * temperature_step;
		}

		contor_pasi_la_o_temperatura = 0;
		while ((contor_pasi_la_o_temperatura < 100) && (n_H > 0.025 * n_part))
		{
			contor_pasi_la_o_temperatura++;

			Dopri5(timp, timp + step_t, eps, step_t, step_t / 4.0, &sol[0]);

			for (int i = 0; i < n_part; i++)
			{
				Medium[i].x = sol[4 * i + 0];
				Medium[i].y = sol[4 * i + 2];
			}
			timp += step_t;

			Temperaturi();

			for (int i = 0; i < n_part; i++)
			{
				if ((Medium[i].raza > 1.05 * radius) && (probabilitateHL[i] > rand_dis(gen)))
				{
					Medium[i].raza = rmic;
					n_H--; n_L++;
				}
				else
				{
					if ((Medium[i].raza < 1.05 * radius) && (probabilitateLH[i] > rand_dis(gen)))
					{
						Medium[i].raza = rmare;
						n_H++; n_L--;
					}
				}
			}
		}

		if (fabs((double)(n_L - n_L_vechi)) > ((double)n_part / 100))
		{
#ifdef grafic
			{
				int p, i, j, v1, v2;
				int count = 0;
				int *count_switched;
				double d1;

				count_switched = (int *)calloc(n_part, sizeof(int));

				for (p = 0; p < n_part; p++)
				{
					count_switched[p] = 0;

					for (i = 0; i < neighbours[p]; i++)
					{
						if (Medium[Position_Coef[p][i].vecin].raza > 1.05)
						{
							count_switched[p]++;
						}
					}

					for (i = 0; i < neighbours[p] - 1; i++)
					{
						for (j = i + 1; j < neighbours[p]; j++)
						{
							v1 = Position_Coef[p][i].vecin;
							v2 = Position_Coef[p][j].vecin;

							d1 = sqrt((Medium[v1].x - Medium[v2].x) * (Medium[v1].x - Medium[v2].x) + (Medium[v1].y - Medium[v2].y) * (Medium[v1].y - Medium[v2].y) + (Medium[v1].z - Medium[v2].z) * (Medium[v1].z - Medium[v2].z));

							if ((d1 < 5.0)) //ca sa nu luam doi vecini ambii de pe Ox sau ambii de pe Oy
							{
								count++;
							}
						}
					}
				}


				char fis_save_vis[500];
				sprintf(fis_save_vis, "%s_ucd_%06d.inp", file, (int)timp);

				FILE *fpout;
				fpout = fopen(fis_save_vis, "w");

				fprintf(fpout, "%d %d 1 0 0\n", n_part, count);
				printf("SAVING UCD %d %d\n", n_part, count);

				for (i = 0; i < n_part; i++)
				{
					fprintf(fpout, "%d %f %f %f\n", i + 1, Medium[i].x, Medium[i].y, Medium[i].z);
				}

				count = 0;
				for (p = 0; p < n_part; p++)
				{
					for (i = 0; i < neighbours[p] - 1; i++)
					{
						for (j = i + 1; j < neighbours[p]; j++)
						{
							v1 = Position_Coef[p][i].vecin;
							v2 = Position_Coef[p][j].vecin;

							d1 = sqrt((Medium[v1].x - Medium[v2].x) * (Medium[v1].x - Medium[v2].x) + (Medium[v1].y - Medium[v2].y) * (Medium[v1].y - Medium[v2].y) + (Medium[v1].z - Medium[v2].z) * (Medium[v1].z - Medium[v2].z));

							if ((d1 < 5.0))  //ca sa nu luam doi vecini ambii de pe Ox sau ambii de pe Oy
							{
								count++;
								fprintf(fpout, "%d 1 tri  %d  %d  %d \n", count, p + 1, v1 + 1, v2 + 1);
							}

						}
					}
				}

				fprintf(fpout, "5 1 1 1 1 1\n");
				fprintf(fpout, "raza, nm\n");
				fprintf(fpout, "phase, au\n");
				fprintf(fpout, "PLH, au\n");
				fprintf(fpout, "PHL, au\n");
				fprintf(fpout, "pressure, au\n");

				for (i = 0; i < n_part; i++)
				{
					fprintf(fpout, "%d %lf %lf %f %f %f\n", i + 1, Medium[i].raza, Medium[i].k,
						((probabilitateLH[i] < 1.0) ? (probabilitateLH[i]) : (1.0)),
						((probabilitateHL[i] < 1.0) ? (probabilitateHL[i]) : (1.0)),
						pres[i]);
				}

				fclose(fpout);

				free(count_switched);
			}

#endif	
			n_L_vechi = n_L;
		}

		printf("Timp %lf \t\t Temp %lf \t\t HS %d \n", timp, T[0], n_H);


		fvol = fopen(fis_volumeMHL, "a");
		// ************************* SALVARI VOL********************************
		fprintf(fvol, "%lf   %lf   %lf\n", timp, T[0], (double)n_H / n_part);
		// ************************ END SALVARI *****************************
		fclose(fvol);
		//}
	}
#endif

#ifdef PHOTO
	//////////////////////////////////////////////////////////////////////////
	// PHOTOEXCITATION
	//////////////////////////////////////////////////////////////////////////
	fvol = fopen(fis_volumePHOTO, "w");

	int n_H = 0;
	int n_L = n_part;
	int n_L_vechi = n_part;

	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_real_distribution<double> rand_dis(0.0, 1.0);  // use with rand_dis(gen)

	timp = t_init;

	for (i = 0; i < n_part; i++)
	{
		if ( (0.05*Medium[i].x/depth) < (Finvers(rand_dis(gen)) - 0.25*0.05) )
		{
			T[i] = T_EXCITATION;
			Medium[i].raza = rmare;
			n_H++; n_L--;
		}
		else
		{
			T[i] = T_LIM_DWN;
			Medium[i].raza = rmic;
	 	}
	}

	for (i = 0; i < n_part; i++) //Conditii initiale
	{
		sol[4 * i + 0] = Medium[i].x;
		sol[4 * i + 1] = 0.0;
		sol[4 * i + 2] = Medium[i].y;
		sol[4 * i + 3] = 0.0;
		//Medium[i].raza = 1.1 * radius;
	}

	int contor_pasi = 0;
#ifdef grafic
	int v1, v2;
	int count = 0;
	int *count_switched;
	double d1;
#endif

	while ((contor_pasi < N_MAX_STEPS))
	{
		Dopri5(timp, timp + step_t, eps, step_t, step_t / 4.0, &sol[0]);

		for (int i = 0; i < n_part; i++)
		{
			Medium[i].x = sol[4 * i + 0];
			Medium[i].y = sol[4 * i + 2];
		}
		timp += step_t;

		TemperaturiExchange();
		Temperaturi();

		for (int i = 0; i < n_part; i++)
		{
			if ((Medium[i].raza > 1.05 * radius) && (T[i] <= T_LIM_UP) && (probabilitateHL[i] > rand_dis(gen)))
			{
				Medium[i].raza = rmic;
				//T[i] += 10.0;			// exotermic la H-to-L
				for (int j = 0; j < neighbours[i]; j++)
				{
					T[Position_Coef[i][j].vecin] += 10.0;
				}
				n_H--; n_L++;
			}
			else
			{
				if ((Medium[i].raza < 1.05 * radius) && (T[i] >= T_LIM_UP) /*(probabilitateLH[i] > rand_dis(gen))*/)		//neconditionat
				{
					Medium[i].raza = rmare;
					n_H++; n_L--;
				}
			}
		}

		arie = Suprafata(false);

		if (((contor_pasi <= 1000) && !(contor_pasi % 25)) || ((contor_pasi > 1000) && !(contor_pasi % 100)))
		{
#ifdef grafic
			{
				count = 0;
				count_switched = (int *)calloc(n_part, sizeof(int));

				for (int p = 0; p < n_part; p++)
				{
					count_switched[p] = 0;

					for (int i = 0; i < neighbours[p]; i++)
					{
						if (Medium[Position_Coef[p][i].vecin].raza > 1.05)
						{
							count_switched[p]++;
						}
					}

					for (int i = 0; i < neighbours[p] - 1; i++)
					{
						for (int j = i + 1; j < neighbours[p]; j++)
						{
							v1 = Position_Coef[p][i].vecin;
							v2 = Position_Coef[p][j].vecin;

							d1 = sqrt((Medium[v1].x - Medium[v2].x) * (Medium[v1].x - Medium[v2].x) + (Medium[v1].y - Medium[v2].y) * (Medium[v1].y - Medium[v2].y) + (Medium[v1].z - Medium[v2].z) * (Medium[v1].z - Medium[v2].z));

							if ((d1 < L + 3.0 * radius)) //ca sa nu luam doi vecini ambii de pe Ox sau ambii de pe Oy
							{
								count++;
							}
						}
					}
				}


				char fis_save_vis[500];
				sprintf(fis_save_vis, "%s_ucd_%06d.inp", file, (int)timp);

				FILE *fpout;
				fpout = fopen(fis_save_vis, "w");

				fprintf(fpout, "%d %d 1 0 0\n", n_part, count);
				printf("SAVING UCD %d %d\n", n_part, count);

				for (int i = 0; i < n_part; i++)
				{
					fprintf(fpout, "%d %f %f %f\n", i + 1, Medium[i].x, Medium[i].y, Medium[i].z);
				}

				count = 0;
				for (int p = 0; p < n_part; p++)
				{
					for (int i = 0; i < neighbours[p] - 1; i++)
					{
						for (int j = i + 1; j < neighbours[p]; j++)
						{
							v1 = Position_Coef[p][i].vecin;
							v2 = Position_Coef[p][j].vecin;

							d1 = sqrt((Medium[v1].x - Medium[v2].x) * (Medium[v1].x - Medium[v2].x) + (Medium[v1].y - Medium[v2].y) * (Medium[v1].y - Medium[v2].y) + (Medium[v1].z - Medium[v2].z) * (Medium[v1].z - Medium[v2].z));

							if ((d1 < L + 3.0 * radius))  //ca sa nu luam doi vecini ambii de pe Ox sau ambii de pe Oy
							{
								count++;
								fprintf(fpout, "%d 1 tri  %d  %d  %d \n", count, p + 1, v1 + 1, v2 + 1);
							}

						}
					}
				}

				fprintf(fpout, "6 1 1 1 1 1 1\n");
				fprintf(fpout, "raza, nm\n");
				fprintf(fpout, "phase, au\n");
				fprintf(fpout, "PLH, au\n");
				fprintf(fpout, "PHL, au\n");
				fprintf(fpout, "Temperature, K\n");
				fprintf(fpout, "pressure, au\n");

				for (int i = 0; i < n_part; i++)
				{
					fprintf(fpout, "%d %lf %lf %lf %lf %lf %lf\n", i + 1, Medium[i].raza, Medium[i].k,
						((probabilitateLH[i] < 1.0) ? (probabilitateLH[i]) : (1.0)),
						((probabilitateHL[i] < 1.0) ? (probabilitateHL[i]) : (1.0)),
						T[i],
						pres[i]);
				}

				fclose(fpout);

				free(count_switched);
			}

#endif
			doTheHisto((int)timp);
			printf("Timp %5.2lf \t Temp %5.2lf \t HS %d \t Surf  %6.4lf \n", timp, T[0], n_H, arie);
		}

		// ************************* SALVARI VOL********************************
		fprintf(fvol, "%lf   %lf   %lf   %lf\n", timp, T[0], (double)n_H / n_part, arie);
		// ************************ END SALVARI *****************************

		contor_pasi++;
	}
	fclose(fvol);

#endif

	printf("\n\n D  O  N  E : \n");

	return 0;
}

//*************************************************************************

int initializare(void)
{
	FILE *fp;
	long i;

	/// READ Medium
	fp = fopen(fis_particule, "r");
	for (i = 0; i < n_part; i++)
	{
		fscanf(fp, "%lG %lG %lG %lG %lG %lG \n", &Medium[i].x, &Medium[i].z, &Medium[i].y, &Medium[i].raza, &Medium[i].theta, &Medium[i].k);
		if (Medium[i].x > depth)
			depth = Medium[i].x;
	}
	fclose(fp);

	///// COMPUTE NEIGHBOURS
	alglib_function_neighbours();

	///// PRINT NEIGHBOURS
// 	for (i = 0; i < n_part; i++)
// 	{
// 		printf("%d: ", i);
// 		for (j = 0; j < neighbours[i]; j++)
// 		{
// 			printf(" %d ", Position_Coef[i][j].vecin);
// 		}
// 		printf("\n");
// 	}

	printf("Citit %d particule  \n", n_part);

	return(0);
}

//*************************************************************************

void alglib_function_neighbours(void)
{
	int	neighbours_max = 0, neighbours_med = 0;
	int i, j, local_index;
	double distance;
	real_2d_array a;

	a.setlength(n_part, 3);
	for (i = 0; i < n_part; i++)
	{
		a(i, 0) = Medium[i].x;
		a(i, 1) = Medium[i].y;
		a(i, 2) = Medium[i].z;
	}

	integer_1d_array tags;
	tags.setlength(n_part);
	for (int i = 0; i < n_part; i++)
	{
		tags(i) = i;
	}

	ae_int_t nx = 3;
	ae_int_t ny = 0;
	ae_int_t normtype = 2;

	kdtree kdt;
	//kdtreebuild(a, nx, ny, normtype, kdt);
	kdtreebuildtagged(a, tags, nx, ny, normtype, kdt);

	real_1d_array x;
	x.setlength(3);
	real_2d_array r = "[[]]";

	integer_1d_array indexes;

	for (i = 0; i < n_part; i++)
	{
		x(0) = Medium[i].x;
		x(1) = Medium[i].y;
		x(2) = Medium[i].z;

		ae_int_t k;
		//k = kdtreequeryknn(kdt, x, 2, false);
		distance = 1.1 * (2.0 * radius + L);
		k = kdtreequeryrnn(kdt, x, distance, false);

		neighbours[i] = (int)k;

		// 		if (neighbours[i] + 1 > n_max_vec - 1)
		// 		{
		// 			printf("\n\n PREA MULTI VECINI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n\n");
		// 			break;
		// 		}
		// 		else
		{
			kdtreequeryresultstags(kdt, indexes);

			//SORT FOR TYPE
			std::vector<int> myvector(neighbours[i]);
			for (j = 0; j < neighbours[i]; j++)
				myvector[j] = (int)indexes(j);

			std::sort(myvector.begin(), myvector.end());

			for (j = 0; j < neighbours[i]; j++)
			{
				local_index = myvector[j];// (int)indexes(j);

				Position_Coef[i][j].vecin = local_index;  //<<<---- asta e vecin
				Position_Coef[i][j].tip_vecin = 0;
				if (local_index < (i - 1))
					Position_Coef[i][j].tip_vecin = 1;
				if (local_index == (i - 1))
					Position_Coef[i][j].tip_vecin = 2;
				if (local_index == (i + 1))
					Position_Coef[i][j].tip_vecin = 3;
				if (local_index > (i + 1))
					Position_Coef[i][j].tip_vecin = 4;
			}
			//n_vec++;
		}

		//printf("particula: %d  cu  %d  vecini \n", i, neighbours[i]);
		neighbours_med += neighbours[i];
		if (neighbours[i] > neighbours_max)	neighbours_max = neighbours[i];
	}
	printf("Numar maxim de vecini: %d    Numar mediu de vecini: %f \n", neighbours_max, (double)neighbours_med / n_part);
	//getchar();
}

//*************************************************************************

int TemperaturiExchange(void)
{
	int j, v;
	double Q;
	for (int i = 0; i < n_part; i++)
	{
// 		if (!(i % 59))	// latura din stanga ia temp substrat
// 		{
// 			T[i] += (T_EXCITATION - T[i]) * CoefTermExt;
// 		}
// 		else
// 		{
			if (neighbours[i] < 5)
			{
				T[i] -= (T[i] - T_LIM_DWN) * CoefTermExt;
				//T[i] = T_LIM_DWN;
			}
			else
			{
				for (j = 0; j < neighbours[i]; j++)
				{
					v = Position_Coef[i][j].vecin;
					Q = (T[i] - T[v]) * CoefTerm;
					T[i] -= Q;
					T[v] += Q;
				}
			}
// 		}
	}
	return(0);
}

//**************************************************************************

int Temperaturi(void)
{
	double Fe, radical;
	int v, j;

	for (int i = 0; i < n_part; i++)
	{
		Fe = 0.0;

		for (j = 0; j < neighbours[i]; j++)
		{
			v = Position_Coef[i][j].vecin;
			radical = sqrt((Medium[v].x - Medium[i].x) * (Medium[v].x - Medium[i].x) + (Medium[v].y - Medium[i].y) * (Medium[v].y - Medium[i].y) + (Medium[v].z - Medium[i].z) * (Medium[v].z - Medium[i].z));

			Fe -= Kf_mic_mic * (radical - Medium[v].raza - Medium[i].raza - L);
		}

		pres[i] = Fe;
		probabilitateLH[i] = 1.0 / tau * exp(-(E + ka * Fe) / T[i]) * exp(-(H - T[i] * S) / T[i]);
		probabilitateHL[i] = 1.0 / tau * exp(-(E - ka * Fe) / T[i]);
	}
	return(0);
}

//**************************************************************************

double Suprafata(bool save)
{
	int i;
	double a, b, c, p;												//laturile triunghiului
	//double x0 = Medium[1785].x, y0 = Medium[1785].y;	//30x30
	double x0 = Medium[4875].x, y0 = Medium[4875].y;	//50x50
	double x1, y1, x2, y2;

	double pante_inainte_ordonare[6 * n_lat - 3];
	double index_margini[6 * n_lat - 3];
	double index_sort[6 * n_lat - 3];

	int contor = 0;
	double Surf;

	/////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////// calcul arie ////
	/////////////////////////////////////////////////////////////////////

	Surf = 0.0;
	contor = 0;
	for (i = 0; i < n_part; i++)
	{
		if (neighbours[i] < 4.9)
		{

			x1 = Medium[i].x;
			y1 = Medium[i].y;
			index_margini[contor] = i;
			index_sort[contor] = contor;
			pante_inainte_ordonare[contor] = atan2((y1 - y0), (x1 - x0));

			contor++;
		}
	}

	//imsl_d_sort((6 * n_lat - 6), pante_inainte_ordonare, IMSL_PERMUTATION_USER, indecsi_sort, IMSL_RETURN_USER, pante_dupa_ordonare, 0);
	gsl_sort2(pante_inainte_ordonare, 1, index_sort, 1, (6 * n_lat - 3));

	for (i = 0; i < contor - 1; i++)
	{
		x1 = Medium[(long)index_margini[(long)index_sort[i]]].x;
		y1 = Medium[(long)index_margini[(long)index_sort[i]]].y;
		x2 = Medium[(long)index_margini[(long)index_sort[i + 1]]].x;
		y2 = Medium[(long)index_margini[(long)index_sort[i + 1]]].y;
		a = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
		b = sqrt((x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2));
		c = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
		p = (a + b + c) / 2.0;
		Surf += sqrt(p * (p - a) * (p - b) * (p - c));
	}
	x1 = Medium[(long)index_margini[(long)index_sort[contor - 1]]].x;
	y1 = Medium[(long)index_margini[(long)index_sort[contor - 1]]].y;
	x2 = Medium[(long)index_margini[(long)index_sort[0]]].x;
	y2 = Medium[(long)index_margini[(long)index_sort[0]]].y;
	a = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
	b = sqrt((x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2));
	c = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
	p = (a + b + c) / 2.0;
	Surf += sqrt(p * (p - a) * (p - b) * (p - c));

	if (save)
	{
		FILE *fp;
		fp = fopen("E:\\Stoleriu\\C\\special\\3d\\res\\2021\\Elastic\\run1.dat", "w");
		for (i = 0; i < contor - 1; i++)
		{
			x1 = Medium[(long)index_margini[(long)index_sort[i]]].x;
			y1 = Medium[(long)index_margini[(long)index_sort[i]]].y;
			x2 = Medium[(long)index_margini[(long)index_sort[i + 1]]].x;
			y2 = Medium[(long)index_margini[(long)index_sort[i + 1]]].y;
			a = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
			b = sqrt((x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2));
			c = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
			p = (a + b + c) / 2.0;
			fprintf(fp, "%ld %ld %lf %lf %lf %lf %lf %lf %lf %lf\n", (long)index_margini[(long)index_sort[i]], (long)index_margini[(long)index_sort[i + 1]], x1, y1, x2, y2, a, b, c, sqrt(p * (p - a) * (p - b) * (p - c)));
		}
		x1 = Medium[(long)index_margini[(long)index_sort[contor - 1]]].x;
		y1 = Medium[(long)index_margini[(long)index_sort[contor - 1]]].y;
		x2 = Medium[(long)index_margini[(long)index_sort[0]]].x;
		y2 = Medium[(long)index_margini[(long)index_sort[0]]].y;
		a = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
		b = sqrt((x0 - x2) * (x0 - x2) + (y0 - y2) * (y0 - y2));
		c = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
		p = (a + b + c) / 2.0;
		fprintf(fp, "%ld %ld %lf %lf %lf %lf %lf %lf %lf %lf\n", (long)index_margini[(long)index_sort[i]], (long)index_margini[(long)index_sort[i + 1]], x1, y1, x2, y2, a, b, c, sqrt(p * (p - a) * (p - b) * (p - c)));

		fclose(fp);
	}

	return Surf;
}

//**************************************************************************

double Finvers(double x)
{
	return( log(1.0 + (x * 2.0 / depth) * 6.3890560989306502272) / 2.6230812603996638992); //log(1.0 + x * (exp(depth) - 1.0)) scalat la depth = 3.0;
}

//**************************************************************************

void doTheHisto(int step)
{
	int binSize = 2, indexBin;
	int numBins = (int)(depth / binSize) + 1;

	vector<double> histo(numBins, 0.0);
	vector<int> countBins(numBins, 0);

	char fis_save_histo[500];
	sprintf(fis_save_histo, "%s_%d.dat", fileHisto, step);

	FILE *fpout;
	fpout = fopen(fis_save_histo, "w");

	for (int i=0; i < n_part; i++)
	{
		indexBin = (int)(fabs(Medium[i].x) / binSize);
		histo[indexBin] += T[i];
		countBins[indexBin]++;
	}

	for (int i=0; i<numBins; i++)
	{
		if (countBins[i])
			histo[i] /= countBins[i];
		fprintf(fpout, "%d  %lf\n", i * binSize, histo[i]);
	}

	fclose(fpout);
}

//**************************************************************************

int Funct_Dopri(double time, double *input, double *deriv)
{
	long i, j, k, kk;
	double radical, Fe, Fex, Fey;
	double alungirea;
	double distanta_normala = L;

	for (i = 0; i < n_part; i++)
	{
		k = 4 * i;
		Fex = 0.0;
		Fey = 0.0;

		for (j = 0; j < neighbours[i]; j++)
		{
			kk = 4 * Position_Coef[i][j].vecin;
			radical = sqrt((input[kk + 0] - input[k + 0]) * (input[kk + 0] - input[k + 0]) + (input[kk + 2] - input[k + 2]) * (input[kk + 2] - input[k + 2]));

			alungirea = (radical - Medium[Position_Coef[i][j].vecin].raza - Medium[i].raza - distanta_normala);

			Fe = Kf_mic_mic * alungirea;

			Fex += Fe * (input[kk + 0] - input[k + 0]) / radical;
			Fey += Fe * (input[kk + 2] - input[k + 2]) / radical;
		}

		deriv[k] = input[k + 1];
		deriv[k + 1] = (Fex - mu * input[k + 1]) / m;
		deriv[k + 2] = input[k + 3];
		deriv[k + 3] = (Fey - mu * input[k + 3]) / m;
	}
	return (0);
}

//*************************************************************************

int Dopri5(double x, double xend, double eps, double hmax, double has, double *sol)

{
	const double uround = 1.0e-7;
	int i;
	int n = n_equ, nmax, reject, nrejct, naccpt, nfcn, nstep;
	double posneg, auxmax, xph, er, auxmax1, auxmax2, denom, hnew, auxmin, aux, fac;

	double verif;

	//#if defined(GPU_aux) && defined(GPU)

	//#else
	double *k1, *k2, *k3, *k4, *k5, *ygrec1;
	k1 = (double *)malloc(n_equ * sizeof(double));
	k2 = (double *)malloc(n_equ * sizeof(double));
	k3 = (double *)malloc(n_equ * sizeof(double));
	k4 = (double *)malloc(n_equ * sizeof(double));
	k5 = (double *)malloc(n_equ * sizeof(double));
	ygrec1 = (double *)malloc(n_equ * sizeof(double));
	//#endif

#ifndef GPU_aux
	for (i = 0; i < n; i++)
	{
		ygrec1[i] = k1[i] = k2[i] = k3[i] = k4[i] = k5[i] = 0.0;
	}
#else
	Dopri5_00_init0 << <Blk, Thrd >> > (&k1[0], &k2[0], &k3[0], &k4[0], &k5[0], &ygrec1[0]);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Dopri5_00_init0 launch failed: %s\n", cudaGetErrorString(cudaStatus));
		return(cudaStatus);
	}
#endif

	nmax = 3000;
	//posneg = sign(1.0, xend-x);
	((xend - x) >= 0.0) ? (posneg = 1.0) : (posneg = -1.0);

	//***************************
	//   initializare
	//***************************

	hmax = fabs(hmax);
	if (1.0e-4 > fabs(has))  auxmax = 1.0e-4;
	else auxmax = fabs(has);
	if (auxmax > hmax)       has = hmax;
	else has = auxmax;
	if (posneg > 0.0)        has = fabs(has);
	else has = -fabs(has);
	//eps    = Amax1r(eps, 7.0*uround);
	if (eps >= 7.0 * uround) (eps = 7.0 * uround);
	reject = 0;
	naccpt = 0;
	nrejct = 0;
	nfcn = 1;
	nstep = 0;

#ifndef GPU
	Funct_Dopri(x, &sol[0], &k1[0]);
#else
	Funct_DopriCUDA << <BlkD, ThrdD >> > (&k1[0], &sol[0], &Medium_k[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0], Ext_force_x, Ext_force_y);
	cudaDeviceSynchronize();
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Funct_DopriCUDA launch failed: %s\n", cudaGetErrorString(cudaStatus));
		return(cudaStatus);
	}
#endif

	//****************************
	//   pas de integrare de baza
	//****************************

	while ((x - xend) * posneg + uround <= 0.0)
	{
		verif = x + 1.0e-1 * has - x;
		if ((nstep > nmax) || (verif == 0.0))
		{
			//AfxMessageBox("OUT!!!");

			return (1);
		}

		if ((x + has - xend) * posneg > 0.0)
			has = xend - x;
		nstep++;

		//****************************
		//   primele 6 etape
		//****************************

#ifndef GPU_DOPRI

		/////////////////////////////////////////////////////////////////////////////
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			ygrec1[i] = sol[i] + has * 2.0e-1 * k1[i];
#else
		Dopri5_01_add1 << <Blk, Thrd >> > (&ygrec1[0], has, &sol[0], &k1[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Dopri5_01_add1 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU
		Funct_Dopri(x + 2.0e-1 * has, &ygrec1[0], &k2[0]);
#else
		Funct_DopriCUDA << <BlkD, ThrdD >> > (&k2[0], &ygrec1[0], &Medium_k[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0], Ext_force_x, Ext_force_y);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Funct_DopriCUDA launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			ygrec1[i] = sol[i] + has * ((3.0e0 / 4.0e1) * k1[i] + (9.0e0 / 4.0e1) * k2[i]);
#else
		Dopri5_02_add2 << <Blk, Thrd >> > (&ygrec1[0], has, &sol[0], &k1[0], &k2[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Dopri5_02_add2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU
		Funct_Dopri(x + 3.0e-1 * has, &ygrec1[0], &k3[0]);
#else
		Funct_DopriCUDA << <BlkD, ThrdD >> > (&k3[0], &ygrec1[0], &Medium_k[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0], Ext_force_x, Ext_force_y);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Funct_DopriCUDA launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			ygrec1[i] = sol[i] + has * ((4.4e1 / 4.5e1) * k1[i] - (5.6e1 / 1.5e1) * k2[i] + (3.2e1 / 9.0e0) * k3[i]);
#else
		Dopri5_03_add3 << <Blk, Thrd >> > (&ygrec1[0], has, &sol[0], &k1[0], &k2[0], &k3[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Dopri5_03_add3 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU
		Funct_Dopri(x + 8.0e-1 * has, &ygrec1[0], &k4[0]);
#else
		Funct_DopriCUDA << <BlkD, ThrdD >> > (&k4[0], &ygrec1[0], &Medium_k[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0], Ext_force_x, Ext_force_y);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Funct_DopriCUDA launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			ygrec1[i] = sol[i] + has * ((1.9372e4 / 6.561e3) * k1[i] - (2.5360e4 / 2.187e3) * k2[i] + (6.4448e4 / 6.561e3) * k3[i] - (2.1200e2 / 7.290e2) * k4[i]);
#else
		Dopri5_04_add4 << <Blk, Thrd >> > (&ygrec1[0], has, &sol[0], &k1[0], &k2[0], &k3[0], &k4[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Dopri5_04_add4 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU
		Funct_Dopri(x + (8.0e0 / 9.0e0) * has, &ygrec1[0], &k5[0]);
#else
		Funct_DopriCUDA << <BlkD, ThrdD >> > (&k5[0], &ygrec1[0], &Medium_k[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0], Ext_force_x, Ext_force_y);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Funct_DopriCUDA launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			ygrec1[i] = sol[i] + has * ((9.0170e3 / 3.1680e3) * k1[i] - (3.5500e2 / 3.3000e1) * k2[i] + (4.6732e4 / 5.2470e3) * k3[i] + (4.9000e1 / 1.7600e2) * k4[i] - (5.1030e3 / 1.8656e4) * k5[i]);
#else
		Dopri5_05_add5 << <Blk, Thrd >> > (&ygrec1[0], has, &sol[0], &k1[0], &k2[0], &k3[0], &k4[0], &k5[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Dopri5_05_add5 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
		xph = x + has;
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU
		Funct_Dopri(xph, &ygrec1[0], &k2[0]);
#else
		Funct_DopriCUDA << <BlkD, ThrdD >> > (&k2[0], &ygrec1[0], &Medium_k[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0], Ext_force_x, Ext_force_y);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Funct_DopriCUDA launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			ygrec1[i] = sol[i] + has * ((3.500e1 / 3.840e2) * k1[i] + (1.100e1 / 8.400e1) * k2[i] + (5.000e2 / 1.113e3) * k3[i] + (1.250e2 / 1.920e2) * k4[i] - (2.187e3 / 6.784e3) * k5[i]);
#else
		Dopri5_06_add6 << <Blk, Thrd >> > (&ygrec1[0], has, &sol[0], &k1[0], &k2[0], &k3[0], &k4[0], &k5[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Dopri5_02_add2 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
		//*******************************
		//   calculul sumei intermediare
		//   pentru economie de memorie
		//*******************************
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			k2[i] = (7.1000e1 / 5.7600e4) * k1[i] + (2.2000e1 / 5.2500e2) * k2[i] - (7.1000e1 / 1.6695e4) * k3[i] + (7.1000e1 / 1.9200e3) * k4[i] - (1.7263e4 / 3.3920e5) * k5[i];
#else
		Dopri5_07_add7 << <Blk, Thrd >> > (&k2[0], &k1[0], &k3[0], &k4[0], &k5[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Dopri5_07_add7 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
		//****************************
		//   ultima etapa
		//****************************
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU
		Funct_Dopri(xph, &ygrec1[0], &k3[0]);
#else
		Funct_DopriCUDA << <BlkD, ThrdD >> > (&k3[0], &ygrec1[0], &Medium_k[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0], Ext_force_x, Ext_force_y);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Funct_DopriCUDA launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#ifndef GPU_aux
		for (i = 0; i < n; i++)
			k4[i] = (k2[i] - (1.0e0 / 4.0e1) * k3[i]) * has;
#else
		Dopri5_08_add8 << <Blk, Thrd >> > (&k4[0], has, &k2[0], &k3[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Dopri5_08_add8 launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////
#else
		Dopri5_XX_FullStep << <BlkD, ThrdD >> > (&ygrec1[0], has, &sol[0], &k1[0], &k2[0], &k3[0], &k4[0], &k5[0], &Medium_raza[0], &Position_Coef_vecin[0], &neighbours[0]);
		cudaDeviceSynchronize();
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Dopri5_XX_FullStep launch failed: %s\n", cudaGetErrorString(cudaStatus));
			return(cudaStatus);
		}
		xph = x + has;
#endif
		nfcn += 6;

		//****************************
		//   estimarea erorii
		//****************************

		er = 0.0e0;
		for (i = 0; i < n; i++)
		{
			//auxmax1 = Amax1r(1.0e-5, fabs(ygrec1[i]));
			//auxmax2 = Amax1r(fabs(sol[i]), 2.0e0 * uround / eps);
			//denom   = Amax1r(auxmax1, auxmax2);
			//er     += (k4[i] / denom) * (k4[i] / denom);

			(1.0e-5 >= fabs(ygrec1[i])) ? (auxmax1 = 1.0e-5) : (auxmax1 = fabs(ygrec1[i]));
			(fabs(sol[i]) >= 2.0e0 * uround / eps) ? (auxmax2 = fabs(sol[i])) : (auxmax2 = 2.0e0 * uround / eps);
			(auxmax1 >= auxmax2) ? (denom = auxmax1) : (denom = auxmax2);
			er += (k4[i] / denom) * (k4[i] / denom);
		};

		er = sqrt(er / n);

		//***************************************
		//   calculul variabilei hnew
		//   se cere 0.2d0 <= hnew / has <= 10.0d0
		//***************************************

		aux = exp(log(er / eps) / 5.0e0);
		//auxmin = Amin1r(5.0e0, aux / 9.0e-1);
		//fac    = Amax1r(1.0e-1, auxmin);
		(5.0e0 <= aux / 9.0e-1) ? (auxmin = 5.0e0) : (auxmin = aux / 9.0e-1);
		(1.0e-1 >= auxmin) ? (fac = 1.0e-1) : (fac = auxmin);
		hnew = has / fac;
		if (er <= eps)
		{
			//****************************
			//   pasul este acceptat
			//****************************

			naccpt++;
			for (i = 0; i < n; i++)
			{
				k1[i] = k3[i];
				sol[i] = ygrec1[i];
			};

			x = xph;
			if (fabs(hnew) > hmax)
				hnew = posneg * hmax;
			//auxmin = Amin1r(fabs(hnew), fabs(has));
			(fabs(hnew) <= fabs(has)) ? (auxmin = fabs(hnew)) : (auxmin = fabs(has));
			if (reject == 1)
				hnew = posneg * auxmin;
			reject = 0;
		};

		//****************************
		//   pasul nu este acceptat
		//****************************

		if (er > eps)
		{
			reject = 1;
			if (naccpt >= 1)
				nrejct++;
		};

		has = hnew;
	};

	free(k1); free(k2); free(k3); free(k4); free(k5); free(ygrec1);

	return(0);
}

//**************************************************************************
