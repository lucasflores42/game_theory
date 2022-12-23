#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "pointers.h"
#include "regular_lattices.h"
#include "mc.h"

/********************************************************************
***                          C.I.                                 ***
********************************************************************/	//gcc dilema.c -lm -lgsl -lgslcblas ---- ./a.out r gama delta rede L ruido
#define dens_aleat
#define dens_inic_c  0.5
#define dens_inic_p  0.
#define dens_inic_pc 0.
#define dens_inic_dc 0.

//#define dens_quad
//#define dens_quad2
//#define dens_listras

/********************************************************************
***                          Plot                                 ***
********************************************************************/
//#define snapshot_gnuplot    
//#define snapshot_gif
//#define snapshot_hexa
//#define snapshot_kagome
#define fps 0000.00000100005 // fps^-1
//#define densidade_arquivo
#define densidade_terminal	
//#define histograma_plot

/********************************************************************
***                        Topologia                              ***
********************************************************************/
//#define viz_well_mixed
#define viz_espacial

//#define topo_aleat
#define topo_fixa

/********************************************************************
***                          Jogo                                 ***
********************************************************************/
//#define PRISONERS_DILEMMA
#define PGG_FOCAL
//#define PGG


//------------------------------------------------------------------------------------

enum tipo_rede {UNIDIMENSIONAL,QUADRADA, CUBICA, QUADRIDIMENSIONAL, HEXAGONAL, KAGOME,
					TRIANGULAR, MOORE};
enum tipo_rede REDE_ATUAL;
enum tipo_estrategia {DEFECTOR, COOPERATOR, PUNISHER};


double RUIDO;
double INV_RUIDO;
#define prob_mobil 0.
#define EPSILON (1e-8)
const int tmax=100000;//5000;
double measure_time = 1;

long long int N;
int L;
long int L2;//= L*L; 
long int L3;// = L*L*L;
long long int L4;// = L*L*L*L; 
int G; // =4 (rede hexagonal) =5 (rede quadrada / kagome) =7 (rede triangular / quadrada 3D) =9 (rede moore / quadrada 4D)

int ND=0;
int NC=0;
int NP=0;
int NPC=0;
int NDC=0;

int nd=0;
int nc=0;
int np=0;
int ndc=0;
int npc=0;

double investimento_total;
double commitment_total;
//------------------------------------------------------------------------------------
unsigned long rseed;
const gsl_rng_type * T;
gsl_rng * rand_vec;

//#define DEBUG

void set_gsl_rng(void)
{
#ifdef DEBUG
	rseed=423;
#else
	rseed=time(NULL);
#endif
  
	gsl_rng_env_setup();
	T    = gsl_rng_default;
	rand_vec = gsl_rng_alloc (T);
	gsl_rng_set (rand_vec, rseed);

  return;
}
/********************************************************************
***                  distribuição dos estados                     ***
********************************************************************/		
void calculo_ci_estado(int state[N], int **viz, double *investimento)
{
	int n;

	for(n=0; n < N; n++)
	{	

		#ifdef dens_aleat
		double temp = gsl_rng_uniform(rand_vec);
		
		if(temp < dens_inic_c ) state[n] = COOPERATOR;
		else 
		{
			if(temp < dens_inic_c + dens_inic_p) state[n] = PUNISHER;
			else 
			{
				if(temp < dens_inic_c + dens_inic_p + dens_inic_pc ) state[n] = 3;
				else
				{	
					if(temp < dens_inic_c + dens_inic_p + dens_inic_pc + dens_inic_dc) state[n] = 4;
					else state[n] = DEFECTOR;
				}
			}			
		}
		#endif
	
		#ifdef dens_quad
		int i,j;
		j=n%L;
		i=n/L;
		state[n] = DEFECTOR;

		if(j>1*L/5 && j<2*L/5 && i>1*L/5 && i<2*L/5){state[j+i*L] = COOPERATOR;investimento[j+i*L] = 0.;}
		
		if(j>3*L/5 && j<4*L/5 && i>1*L/5 && i<2*L/5)
		{
			double temp = gsl_rng_uniform(rand_vec);
			if(temp < 0.5) {investimento[j+i*L] = 2.;state[j+i*L] = COOPERATOR;}
			else{investimento[j+i*L] = 2.0;state[j+i*L] = COOPERATOR;}
		}
		if(j>1*L/5 && j<2*L/5 && i>3*L/5 && i<4*L/5)
		{
			double temp = gsl_rng_uniform(rand_vec);
			if(temp < 1./3) {state[j+i*L] = COOPERATOR;investimento[j+i*L] = 0.5;}
			else if(temp < 2./3) {state[j+i*L] = COOPERATOR;investimento[j+i*L] = 1.0;}			
			else{state[j+i*L] = COOPERATOR;investimento[j+i*L] = 2.0;}
		}
		if(j>3*L/5 && j<4*L/5 && i>3*L/5 && i<4*L/5){state[j+i*L] = COOPERATOR;investimento[j+i*L] = 0.;}
		#endif

		#ifdef dens_quad2
		int i,j;
		j=n%L;
		i=n/L;
		//if(j<L/2){state[j+i*L] = COOPERATOR; investimento[j+i*L] = 0.5;}
		//else{state[j+i*L] = COOPERATOR; investimento[j+i*L] = 2.0;}
		state[j+i*L] = COOPERATOR; investimento[j+i*L] = 2.0;
		if(j==L/3 || j==2*L/3 || i==L/3 || i==2*L/3){state[j+i*L] = COOPERATOR;investimento[j+i*L] = 0.;}
		if(j>L/3 && j<2*L/3 && i>L/3 && i<2*L/3){state[j+i*L] = COOPERATOR;investimento[j+i*L] = 0.5;}
		#endif

		#ifdef dens_listras
		int i,j;
		j=n%L;
		i=n/L;
		if(i<L/3){state[j+i*L] = COOPERATOR;investimento[j+i*L] = .5;}
		if(i>=L/3 && i<2*L/3){state[j+i*L] = COOPERATOR;investimento[j+i*L] = 1.0;}	
		if(i>=2*L/3){state[j+i*L] = COOPERATOR;investimento[j+i*L] = 2.0;}
		#endif
		
	//state[n] = 0;
	}
	
	return;
	
}
/********************************************************************
***                  distribuição da topologia                    ***
********************************************************************/
void calculo_ci_topologia(int topologia[N])
{
	int i;
	for(i=0;i<N;i++)
	{
		#ifdef topo_aleat
		double temp = gsl_rng_uniform(rand_vec);
		if( temp < 1/3. ) topologia[n] = QUADRADA; //ou KAGOME
		else
		{
			if( temp < 2/3.) topologia[n] = TRIANGULAR; // ou CUBICA
			else 
				topologia[n] = MOORE;
		}
		#endif

		#ifdef topo_fixa		
		topologia[i] = G;
		#endif
		//printf("%d\n",topologia[i]);
	}
}

//#define M_PI
/********************************************************************
***                  distribuição da contribuiçao                 ***
********************************************************************/
void calculo_ci_investimento(double  *investimento, int state[N], double gama)
{	

	int i;
	for(i=0;i<N;i++)
	{
		double ga1 = gsl_rng_uniform(rand_vec);
		double ga2 = gsl_rng_uniform(rand_vec);
		double gaussian = 4 + 1*sqrt(-2*log(ga1))*cos(2*M_PI*ga2); // media 3 desvio 1
		int k = 1 + (int) 3*gsl_rng_uniform(rand_vec);
		double kk;
		if(k==1){kk=1.0;}
		if(k==2){kk=2.0;}
		if(k==3){kk=0.5;}
		if(k==4){kk=1.5;}
		if(k==5){kk=1.0;}
		if(k==6){kk=1.5;}		

		if(state[i] == 0){investimento[i] = 0;}
		if(state[i] == 1){investimento[i] = 1 ;}
		if(state[i] == 2){investimento[i] = 1 ;}		
		
		/*int jj,j;
		j=i%L;
		jj=i/L;
		if(jj<L/3){investimento[j+jj*L] = 1;}
		if(jj>=L/3 && jj<2*L/3){investimento[j+jj*L] = k;}	
		if(jj>=2*L/3){investimento[j+jj*L] = 0;}*/
		
		//printf(" %d %d %lf\n",i,state[i],investimento[i]);
	}
}
/********************************************************************
***                  distribuição do comprometimento              ***
********************************************************************/
void calculo_ci_commitment(double  *commitment, int state[N])
{	
	int i;
	for(i=0;i<N;i++)
	{
		double a = gsl_rng_uniform(rand_vec);

		if(a <= 0.5){commitment[i] = 1 ;}
		else{commitment[i] = 0 ;}

		//if(state[i] == 0){investimento[i] = 0;}
		//if(state[i] == 1){investimento[i] = 1 ;}
		//if(state[i] == 2){investimento[i] = 2 ;}		

	}
}
/********************************************************************
***           		          pool do grupo                       ***
********************************************************************/
void calculo_pool( int state[N], int sitio2, int sitio, int vizinho, int **viz, int topologia[N], double *investimento, double *commitment)
{

	int i;
	investimento_total = 0.;
	commitment_total = 0.;

	#ifdef viz_espacial

		for(i=0;i<topologia[sitio2];i++)
		{
			investimento_total += investimento[viz[sitio2][i]];
			commitment_total += commitment[viz[sitio2][i]];
		}
	#endif
	#ifdef viz_well_mixed

		for(i=2;i<topologia[sitio2];i++) //i=1 se nao conta o invadido e i=2 se sim
		{
			int k;
			k = (int) N*gsl_rng_uniform(rand_vec);
			investimento_total += state[k]*investimento[k];

		}
		//conta sitio invadido
		//comentar if/else se i=1
		if(sitio2 != viz[sitio][G-1])
		{
			int k;
			k = (int) N*gsl_rng_uniform(rand_vec);
			investimento_total += state[k]*investimento[k];				

		}
		else 
		{
			investimento_total += state[vizinho]*investimento[vizinho];
			//printf("%d %d\n",sitio2,viz[sitio][G-1]);
		}
		//conta o proprio sitio
		investimento_total += state[sitio]*investimento[sitio];

	#endif

	return;
}
/********************************************************************
***                  numero de cada estrategia                    ***
********************************************************************/
void calculo_numero( int state[N], int sitio2, int sitio, int **viz, int topologia[N])
{
	int i=0;
	
	nd=0;
	nc=0;
	np=0;
	ndc=0;
	npc=0;

	#ifdef viz_espacial
	for(i=0;i<topologia[sitio2];i++) // vizinho 0 e o proprio sitio
	{
		switch(state[viz[sitio2][i]])
		{
			case DEFECTOR:   ++nd; break;
			case COOPERATOR: ++nc; break;
			case PUNISHER:   ++np; break;
			case 3: ++npc; break;
			case 4: ++ndc; break;
			default: 
				fprintf(stderr,"ERRO calculo vizinhos\n");
				fflush(stderr);
		}
	}
	#endif

	#ifdef viz_well_mixed
	for(i=2;i<topologia[sitio2];i++) //i=1 se nao conta o invadido
	{
		int k;
		k = (int) N*gsl_rng_uniform(rand_vec);
		
		if(k != sitio2)
		{
			if(state[k]==0) 
			{
			 	nd = nd+1; 		
		  	}
			if(state[k]==2) 
			{
			 	np = np+1; 	
		  	}
			if(state[k]==1) 
			{
			 	nc = nc+1; 	
		  	}
		}
	}

	//conta o proprio sitio
	switch(state[viz[sitio][0]])
	{
		case DEFECTOR:   ++nd; break;
		case COOPERATOR: ++nc; break;
		case PUNISHER:   ++np; break;
		default:
			fprintf(stderr,"ERRO contagem estrategia proprio sitio\n");
			fflush(stderr);	
	}
	//conta o invasdido sitio
	int aaa = abs(state[sitio] - 1);
	if(aaa==0) 
	{
	 	nd = nd+1; 		
  	}
	if(aaa==1) 
	{
	 	nc = nc+1; 	
  	}
  	
	#endif
	//printf("nc=%d nd=%d\n",state[sitio],aaa);
	return;
}

/********************************************************************
***                          Payoff                               ***
********************************************************************/

void calculo_payoff ( double *payoff, double r, double gama, double delta, int x, int y, int topologia[N], int state[N], double *investimento, int **viz )
{

	/*//puniçao c/ propina
	double q = gsl_rng_uniform(rand_vec);
	double p = gsl_rng_uniform(rand_vec);

	payoff[COOPERATOR] = (1./G)*r*c* (np+nc) - c; //C
	payoff[PUNISHER]   = (1./G)*r*c* (np+nc) - c + p*q*delta*nd - (1 - p*q)*gama*nd; //P
	payoff[DEFECTOR]   = (1./G)*r*c* (np+nc) - (1 - p*q)*gama*np - p*q*delta*np; //D*/


	/*//Ising
	payoff[COOPERATOR] = 2*1*(nc-nd + 0.0); 
	payoff[DEFECTOR]   = 2*(-1)*(nc-nd + 0.0); //W = exp(-payoff/T)*/
	
	//puniçao 
	#if defined(PGG) || defined(PGG_FOCAL)

	double pool = (r/topologia[x])*investimento_total;

	payoff[COOPERATOR] = pool - investimento[y];
	payoff[PUNISHER]   = pool - investimento[y]; // - gama*nd;
	payoff[DEFECTOR]   = pool - investimento[y]; // - delta*np;

	//corrupçao com propina
	/*double pool = (r/topologia[x])*(nc+np);

	payoff[COOPERATOR] = pool - (r/topologia[x])*gama*np + delta*gama*np/nc;
	payoff[DEFECTOR]   = pool - (r/topologia[x])*gama*np;
	
	if(nc==0)
	{
	payoff[PUNISHER]   = pool - (r/topologia[x])*gama*np + (1)*gama;
	}
	else
	{
	payoff[PUNISHER]   = pool - (r/topologia[x])*gama*np + (1-delta)*gama;
	}*/
	
	

	#endif
	
	
	#ifdef PRISONERS_DILEMMA
	double b=3.6;
	double g=0.1;	
	switch(state[x])
	{		

		case COOPERATOR:
			//payoff[COOPERATOR] = nc;   // <- com autointeracao   //RR*(nc - 1) + SS*nd; // sem autointeracao
			payoff[COOPERATOR] = (b-g)*(nc-1) - g*nd;
			break;
		case DEFECTOR: 
			//payoff[DEFECTOR]   = TT*nc;// <- com autointeracao  //TT*nc + PP*(nd-1);   // sem autointeracao
			payoff[DEFECTOR]   = b*nc ;
			break;
	}
	#endif
	
	/*//recompensa
	payoff[COOPERATOR] = (1./G)*r*c* (np+nc) - c + delta*np;
	payoff[PUNISHER]   = (1./G)*r*c* (np+nc) - c - gama*nc;
	payoff[DEFECTOR]   = (1./G)*r*c* (np+nc);*/



	return;
}

/********************************************************************
***                         Update rule                           ***
********************************************************************/
void update_rule( int x, int vizinho, int state[N], double Px, double Py ,int topologia[N], double *investimento, int t)
{

	double Wxy = 1.0/(1.0 + exp(-(Py-Px)*INV_RUIDO));//Wxy = 1.0/(1.0 + exp(-(Py-Px)/RUIDO));
	//double Wxy = exp(-Px/r); //Ising

	double l = gsl_rng_uniform(rand_vec);

	//printf("Px=%lf Py=%lf delta=%lf W=%lf \n",Px,Py,Py-Px,Wxy);

	//troca de estado
	if(Wxy > l)
	//if(Py > Px)
	{   
		//if(state[x]==1) //comentar e ver todas transiçoes que acontecem
		//{histograma[(int) ((invX-invY)/0.5)] += 1;}
   
		state[x] = state[vizinho];
		topologia[x] = topologia[vizinho];	
		investimento[x] = investimento[vizinho]; // tirar pra fixa investimento com local da rede
		//state[x] = fabs(state[x]-1); //Ising
	}
}


/********************************************************************
***                         MCS                              ***
********************************************************************/
void calculo_mcs(double *payoff, int state[N], int **viz, double r, double gama, double delta, int topologia[N], double *investimento, double *commitment, int t)
{

	int n;
	int x,y,sitio,vizinho;
	double Px, Py;
	int sitio2, vizinho2;	

	for(n=0; n<N; n++)
	{					
		Px = 0;
		Py = 0;

		x = (int) N*gsl_rng_uniform(rand_vec);
   		sitio=x;

		#ifdef viz_espacial
		y = 1 + (int) (topologia[x]-1)*gsl_rng_uniform(rand_vec); 
		vizinho = viz[x][y];
		#endif

		#ifdef viz_well_mixed
		y = (int) N*gsl_rng_uniform(rand_vec);
		vizinho = y;
		#endif
			
		if (investimento[vizinho] != investimento[x] || state[vizinho] != state[x])
		{
		#if defined(PGG_FOCAL)

			#ifdef viz_espacial
			sitio2 = sitio;
			vizinho2 = vizinho;			
			#endif		
			#ifdef viz_well_mixed
			sitio2 = viz[sitio][G-1];
			vizinho2 = viz[vizinho][G-1];
			#endif

			calculo_numero(state,sitio2,sitio,viz,topologia);
			calculo_pool(state,sitio2,sitio,vizinho,viz,topologia,investimento,commitment);
			calculo_payoff(payoff,r,gama,delta,sitio,sitio,topologia,state,investimento,viz);
			Px = payoff[state[sitio]];

			calculo_numero(state,vizinho2,vizinho,viz,topologia);
			calculo_pool(state,vizinho2,vizinho,sitio,viz,topologia,investimento,commitment);
			calculo_payoff(payoff,r,gama,delta,vizinho,vizinho,topologia,state,investimento,viz);
			Py = payoff[state[vizinho]];

		#endif	

		#ifdef PGG

			for(y=0;y<topologia[sitio];y++)
			{
				sitio2 = viz[sitio][y];
				calculo_numero(state,sitio2,sitio,viz,topologia);
				calculo_pool(state,sitio2,sitio,vizinho,viz,topologia,investimento,commitment);
				calculo_payoff(payoff,r,gama,delta,sitio2,sitio,topologia,state,investimento,viz);
				Px += payoff[state[sitio]];
				//printf("%d %d %d %d %lf \n", sitio, state[sitio],sitio2, state[sitio2],Px);
			}

			for(y=0;y<topologia[vizinho];y++)
			{			
				vizinho2 = viz[vizinho][y];
				calculo_numero(state,vizinho2,vizinho,viz,topologia);
				calculo_pool(state,vizinho2,vizinho,sitio,viz,topologia,investimento,commitment	);
				calculo_payoff(payoff,r,gama,delta,vizinho2,vizinho,topologia,state,investimento,viz);
				Py += payoff[state[vizinho]];
				//printf("%d %d %d %d %lf \n", vizinho, state[vizinho],vizinho2, state[vizinho2],Py);

			}
			
		#endif
		
		#if defined(PRISONERS_DILEMMA)

			#ifdef viz_espacial
			sitio2 = sitio;
			vizinho2 = vizinho;			
			#endif		
			#ifdef viz_well_mixed
			sitio2 = viz[sitio][G-1];
			vizinho2 = viz[vizinho][G-1];
			#endif

			calculo_numero(state,sitio2,sitio,viz,topologia);
			calculo_pool(state,sitio2,sitio,vizinho,viz,topologia,investimento);
			calculo_payoff(payoff,r,gama,delta,sitio,topologia,state,investimento,viz);
			Px = payoff[state[sitio]];

			calculo_numero(state,vizinho2,vizinho,viz,topologia);
			calculo_pool(state,vizinho2,vizinho,sitio,viz,topologia,investimento);
			calculo_payoff(payoff,r,gama,delta,vizinho,topologia,state,investimento,viz);
			Py = payoff[state[vizinho]];

		#endif
		
			update_rule(x,vizinho,state,Px,Py,topologia,investimento,t);
		}//if (state[vizinho] != state[x])
  	}
}
/********************************************************************
***                          Mobilidade                           ***
********************************************************************/				
void calculo_mobilidade( int state[N], int **viz)
{
	
	double g = gsl_rng_uniform(rand_vec);
	if(g<=prob_mobil)
	{		
		int y = 1 + (int) 4*gsl_rng_uniform(rand_vec);
		int x = (int) N*gsl_rng_uniform(rand_vec);
		int e = state[x];
		state[x] = state[viz[x][y]];
		state[viz[x][y]] = e;
	}
}
/********************************************************************
***		                  Hoshen - Kopelman                       ***
********************************************************************/				
void calculo_percolacao( int state[N], int **viz, int label[N])
{
	int largest_label = 0;
	int n,i;

	for(n=0;n<N;n++)
	{
		label[n] = state[n]*(n+1);
	}
	
	for(n=0;n<N;n++)
	{
		if(state[n]==1)
		{
			if(state[viz[n][1]] == 0 && state[viz[n][4]] != 0)
			{
				int a = label[n];
				label[n] = label[viz[n][4]];
			
				for(i=0;i<N;i++)
				{
					if(label[i]==a)
					{label[i]=label[n];}
				}
			}
			else if(state[viz[n][1]] != 0 && state[viz[n][4]] == 0)
			{
				int a = label[n];
				label[n] = label[viz[n][1]];
			
				for(i=0;i<N;i++)
				{
					if(label[i]==a)
					{label[i]=label[n];}
				}
			}
			else if(state[viz[n][1]] != 0 && state[viz[n][4]] != 0)
			{
				int a = label[viz[n][1]];
				int b = label[n];
				label[viz[n][1]] = label[viz[n][4]];
				label[n] = label[viz[n][4]];
			
				for(i=0;i<N;i++)
				{
					if(label[i]==a || label[i]==b)
					{label[i]=label[n];}
				}
			}
		}
	}
}
/********************************************************************
***		                  contribuiçao por cluster                ***
********************************************************************/				
void calculo_tamanho_cluster( int state[N], int **viz, int label[N], double *investimento, int t)
{
	int n,m,flag;
	double tamanho,inv_medio, desvio, desvio_acum;
	long double inve, inve2,inve22;
	desvio_acum = 0.;

	for(m=0;m<N;m++)
	{	
		if(label[m]!=0)
		{
			tamanho = 0.;
			inve = 0.;
			inve2 = 0.;
			desvio = 0.;
			//inv_medio = 0.;

			for(n=0;n<N;n++)
			{
				if(label[n]==label[m])
				{
					tamanho += 1.;
					inve += investimento[n];
					inve2 += investimento[n]*investimento[n];
					inve22 = inve*inve;
					//desvio += sqrt(inve2/tamanho-inve22/(tamanho*tamanho));
					desvio += inve2/tamanho-inve22/(tamanho*tamanho); // variancia
				}
			}
			inv_medio = inve/tamanho;
			desvio_acum += desvio/tamanho; // soma desvio de cada cluster por passo
			//printf("%d %d %lf %lf\n", t, label[m], inv_medio, desvio);
		}
	}
	printf("%d %lf\n", t, desvio_acum);
	//printf("%d %.26Lf %.26Lf %lf %Lf\n",t, inve22/(tamanho*tamanho),inve2/tamanho,tamanho, inve2/tamanho-inve22/(tamanho*tamanho));
}

/********************************************************************
***                           Histograma                          ***
********************************************************************/				
void calculo_histograma( int t, double *investimento, int histograma[N], int state[N])
{
	int i,j;
	
	for(j=0;j<10;j++)
	{
		histograma[j] = 0;
	}
	for(i=0;i<N;i++)
	{
		if(investimento[i]==0.5 && state[i]==1){histograma[0] += 1;}
		else if(investimento[i]==1.0 && state[i]==1){histograma[1] += 1;}
		else if(investimento[i]==1.5 && state[i]==1){histograma[2] += 1;}
		else if(investimento[i]==2.0 && state[i]==1){histograma[3] += 1;}
		else if(investimento[i]<=1.0 && state[i]==1){histograma[4] += 1;}
		else if(investimento[i]<=1.2 && state[i]==1){histograma[5] += 1;}
		else if(investimento[i]<=1.4 && state[i]==1){histograma[6] += 1;}
		else if(investimento[i]<=1.6 && state[i]==1){histograma[7] += 1;}
		else if(investimento[i]<=1.8 && state[i]==1){histograma[8] += 1;}	
		else if(investimento[i]<=2.0 && state[i]==1){histograma[9] += 1;}											
	}

printf("%d %lf %lf %lf %lf\n",t, (double)histograma[0]/N, (double)histograma[1]/N, (double)histograma[2]/N, (double)histograma[3]/N);
}
/********************************************************************
***                          Densidades                           ***
********************************************************************/
void calculo_densidades(int state[N], double *investimento, int t, FILE *fp)
{
	int k;
	long double inv=0.;
	long double inv2=0.;


	ND=0;
	NC=0;
	NP=0;
	NPC=0;
	NDC=0;
	
	for(k=0;k<N;k++)
	{
		switch (state[k])
		{
			case COOPERATOR: ++NC; break;
			case PUNISHER:   ++NP; break;
			case DEFECTOR:   ++ND; break;
			case 3: ++NPC; break;
			case 4: ++NDC; break;
			default: 
				fprintf(stderr,"ERRO - tipo de estrategia\n");
				fflush(stderr);
		}	
		if(state[k] == 1)
		{
			inv += investimento[k];
			inv2 += investimento[k]*investimento[k];			
		}
	}	

	long double inv22 = inv*inv;
	double desvio = (double)sqrt((inv2-(inv*inv)/(NC))/NC);
 	double inve = (double)inv/(NC);
	
	if(NC == 0)
	{
		desvio = 0.;
		inve = 0.;
	}

	#ifdef densidade_terminal
	printf("%d %lf %lf %lf %lf %lf\n", t, (double)NC/(N), (double)ND/(N), (double)NP/(N), inve ,desvio);
	#endif

	#ifdef densidade_arquivo
	fprintf(fp,"%d %lf %lf %lf %lf %lf\n", t, (double)NC/(N), (double)ND/(N), (double)NP/(N), inve ,desvio);	
	#endif

}

/********************************************************************
***                        Snapshots                              ***
********************************************************************/
void snap_gnuplot(int state[N], int topologia[N], double *investimento, int label[N], int t)	//    ./a.out | gnuplot
{
	// plota de cima pra baixo, esquerda pra direita
	// gnuplot inverte baixo e cima
	int i,j;
	
	printf("set title \"MCS = %d\" \n",t);
	printf("set autoscale keepfix\n");
	printf("set palette model RGB\n");
	printf("unset border\n");
	printf("unset colorbox\n");	
	printf("unset xtics\n");
	printf("unset ytics\n");
	printf("set palette defined ( 0 \"dark-red\", 0.5 \"light-red\",  1  \"#0000B3\", 2 \"#000057\")\n");
	//printf("set palette defined ( 0 \"#a6611a\", 0.5 \"#dfc27d\",  1  \"#80cdc1\", 2 \"#018571\")\n");
	printf("set cbrange[0:2]\n");
	printf("set xrange[0:%d]\n",L);
	printf("set yrange[0:%d]\n",L);
	printf("set size square\n");
			
	printf("plot \"-\" matrix with image\n");
			

	for(i=0;i<L;i++)
	{
		/*for(j=0;j<L;j++)
		{
			printf("%d ",state[j+i*L]);
		
		}*/
		for(j=0;j<L;j++)
		{
			printf("%d ",state[j+i*L]);
		}
		printf("\n");
	}
	printf("\n"); 
	printf("e\n");    printf("pause(%lf)\n",fps);
	
}
void snap_gif(int state[N], int topologia[N], int t)
{

	int i,j;
	
	for(i=0;i<L;i++)
	{
		for(j=0;j<L;j++)
		{
			printf("%d ",state[j+i*L]);
		
		}
		/*for(j=0;j<L;j++)
		{
		printf("%d ",topologia[j+i*L]);
		}*/
		printf("\n");
	}
	printf("\n"); 

}
void snap_kagome(int state[N], int topologia[N], int t)
{

	int n,i,j; //n = i + 2j N
	double x,y;

	printf("set title \"MCS = %d\" \n",t);
	printf("set object 1 rect from graph 0, graph 0 to graph 1, graph 1 back\n");
	printf("set object 1 rect fc rgb \"black\" fillstyle solid 1.0 \n");
	printf("set view map\n");
	printf("set size ratio .9\n");
	//printf("set xrange[30:50]\n");
	//printf("set yrange[30:60]\n");
	printf("set xrange[00:100]\n");
	printf("set yrange[0:40]\n");
	printf("set palette defined ( 0 \"red\", 1 \"blue\",  2 \"green\")\n");
	printf("set cbrange[0:2]\n");
	printf("unset key\n");
	printf("splot \"-\" using 1:2:3 with points pointtype 7 pointsize 1. palette linewidth 10\n");

	for(n=0;n<N;n++)
	{
	//int L = 100;
	int d = 1;
	double h = d*pow(3/4,1/2);
	int n_celula = (int) (n/3);
	int linha_celula = (int) (n_celula/L);
	int coluna_celula = (n_celula%(L));	
		if(n % 3 == 0)
		{
			x = 2*d*coluna_celula + linha_celula*d;
			y = 0 + 2*h*linha_celula;
		}
		if(n % 3 == 1)
		{
			x = d + 2*d*coluna_celula  + linha_celula*d;
			y = 0 + 2*h*linha_celula;
		}
		if(n % 3 == 2)
		{
			x = (d/2.) + 2*d*coluna_celula  + linha_celula*d;
			y = h + 2*h*linha_celula;
		}
		printf("%lf %lf %d\n",x,y,state[n]);
		//printf("%lf %lf %d %lf %lf\n",x,y,a,h,2*h* (int) (n_celula/L));
	}
	printf("\n"); 
	printf("e\n");    printf("pause(%lf)\n",fps);
}
void snap_hexa(int state[N], int topologia[N], int t)
{

	int n,i,j; //n = i + 2j N
	double x,y,d=1,h=70.;

	printf("set title \"MCS = %d\" \n",t);
	printf("set object 1 rect from graph 0, graph 0 to graph 1, graph 1 back\n");
	printf("set object 1 rect fc rgb \"black\" fillstyle solid 1.0 \n");
	printf("set view map\n");
	printf("set size ratio .9\n");
	//printf("set xrange[00:100]\n");
	//printf("set yrange[0:2000]\n");
	printf("set palette defined ( 0 \"red\", 1 \"blue\",  2 \"green\")\n");
	printf("set cbrange[0:2]\n");
	printf("unset key\n");
	printf("splot \"-\" using 1:2:3 with points pointtype 7 pointsize 1. palette linewidth 10\n");

	for(n=0;n<N;n++)
	{
		i = n % (2*L);
		j = n / (2*L);

		x = d*(i+j);
		y = j*(L+h) - h*(i%2);

		printf("%lf %lf %d\n",x,y,state[n]);

	}
	printf("\n"); 
	printf("e\n");    printf("pause(%lf)\n",fps);
}

/********************************************************************
***                     Initialize Network                        ***
********************************************************************/
void initialize_network(char REDE_ATUAL_NAME[], char filename[], double r, double gama, double delta, unsigned long rseed)
{	    
       if (!strcmp(REDE_ATUAL_NAME, "unidimensional") || !strcmp(REDE_ATUAL_NAME, "UNIDIMENSIONAL"))
	{
		N = L; G = 3;
		REDE_ATUAL = UNIDIMENSIONAL;
		//fprintf(stderr,"%s escolhida!\n",REDE_ATUAL_NAME); fflush(stderr);
		sprintf(filename,"unidimensional_dados_r%f_g%f_d%f_seed%ld.txt",r,gama,delta,rseed);
	} 
    	if (!strcmp(REDE_ATUAL_NAME, "quadrada") || !strcmp(REDE_ATUAL_NAME, "QUADRADA"))
	{
		N = L2; G = 5;
		REDE_ATUAL = QUADRADA;
		//fprintf(stderr,"%s escolhida!\n",REDE_ATUAL_NAME); fflush(stderr);
		sprintf(filename,"quadrada_dados_r%f_g%f_d%f_seed%ld.txt",r,gama,delta,rseed);
	}
	if (!strcmp(REDE_ATUAL_NAME, "cubica") || !strcmp(REDE_ATUAL_NAME, "CUBICA"))
	{
		N = L3; G = 7;
		REDE_ATUAL = CUBICA;
		//fprintf(stderr,"%s escolhida!\n",REDE_ATUAL_NAME); fflush(stderr);
		sprintf(filename,"cubica_dados_r%f_g%f_d%f_seed%ld.txt",r,gama,delta,rseed);
	} 
	if (!strcmp(REDE_ATUAL_NAME, "quadridimensional") || !strcmp(REDE_ATUAL_NAME, "QUADRIDIMENSIONAL"))
	{
		N = L4;	G = 9;
		REDE_ATUAL=QUADRIDIMENSIONAL;
		//fprintf(stderr,"%s escolhida!\n",REDE_ATUAL_NAME); fflush(stderr);
		sprintf(filename,"quadridim_dados_r%f_g%f_d%f_seed%ld.txt",r,gama,delta,rseed);
	}
	if (!strcmp(REDE_ATUAL_NAME, "hexagonal") || !strcmp(REDE_ATUAL_NAME, "HEXAGONAL"))
	{
		N = 2*L2; G = 4;
		REDE_ATUAL=HEXAGONAL;
		//fprintf(stderr,"%s escolhida!\n",REDE_ATUAL_NAME); fflush(stderr);
		sprintf(filename,"hexagonal_dados_r%f_g%f_d%f_seed%ld.txt",r,gama,delta,rseed);
	} 
	if (!strcmp(REDE_ATUAL_NAME, "kagome")  || !strcmp(REDE_ATUAL_NAME, "KAGOME"))
	{
		N = 3*L2; G = 5;
		REDE_ATUAL=KAGOME;
		//fprintf(stderr,"%s escolhida!\n",REDE_ATUAL_NAME); fflush(stderr);
		sprintf(filename,"kagome_dados_r%f_g%f_d%f_seed%ld.txt",r,gama,delta,rseed);
	}
	if (!strcmp(REDE_ATUAL_NAME, "triangular")  || !strcmp(REDE_ATUAL_NAME, "TRIANGULAR"))
	{
		N = L2; G = 7;
		REDE_ATUAL=TRIANGULAR;
		//fprintf(stderr,"%s escolhida!\n",REDE_ATUAL_NAME); fflush(stderr);
		sprintf(filename,"triangular_dados_r%f_g%f_d%f_seed%ld.txt",r,gama,delta,rseed);
	}	
	if (!strcmp(REDE_ATUAL_NAME, "moore")  || !strcmp(REDE_ATUAL_NAME, "MOORE"))
	{
		N = L2; G = 9;
		REDE_ATUAL=MOORE;
		//fprintf(stderr,"%s escolhida!\n",REDE_ATUAL_NAME); fflush(stderr);
		sprintf(filename,"moore_dados_r%f_g%f_d%f_seed%ld.txt",r,gama,delta,rseed);
	}
	/*if (!strcmp(REDE_ATUAL_NAME, "random")  || !strcmp(REDE_ATUAL_NAME, "RANDOM"))
	{
		N = L2; G = 5;
		REDE_ATUAL=RANDOM;
		//fprintf(stderr,"%s escolhida!\n",REDE_ATUAL_NAME); fflush(stderr);
		sprintf(filename,"random_dados_r%f_g%f_d%f_seed%ld.txt",r,gama,delta,rseed);
	}*/	
	//fprintf(stderr,"%s escolhida!\n G=%d\nN=%lld\n",REDE_ATUAL_NAME,G,N); fflush(stderr);
	return;
}
/********************************************************************
***                          Main                                 ***
********************************************************************/
int main(int argc, char *argv[])
{	
	//gera arquivo
	char filename[200];
	char REDE_ATUAL_NAME[50];	
	FILE *fp;		
  	int t=0;	
	int **viz;
	double payoff[5], gama, delta,  r;

	if (argc!=7)
	{
		printf("O programa precisa de 6 argumentos:\n1) r\n2) gama\n3) delta\n4) tipo de rede\n5) L\n6) RUIDO");
		exit(1);
	}
	r=atof(argv[1]);
	gama=atof(argv[2]);
	delta=atof(argv[3]);
	sprintf(REDE_ATUAL_NAME,"%s",argv[4]);
	L=atoi(argv[5]);
	RUIDO=atof(argv[6]); 
    INV_RUIDO = 1.0/RUIDO;

   	L2 = L*L; 
	L3 = L2*L;
	L4 = L3*L; 
      
	set_gsl_rng(); // algo da gsl   
	initialize_network(REDE_ATUAL_NAME,filename,r,gama,delta,rseed);  
	viz = create_2d_int_pointer_h(N,G);	

	int state[N],topologia[N], histograma[N], label[N];
	double investimento[N], commitment[N];

	fp=fopen(filename,"w"); //abre .dat
//------------------------------------------------------------------------------------
	switch (REDE_ATUAL)
	{
		case UNIDIMENSIONAL: unidimensional_lattice(viz,L);
			break;
		case QUADRADA: square_lattice(viz,L);
			break;
		case CUBICA: cubic_lattice(viz,L);
			break;
		case QUADRIDIMENSIONAL: hypercubic_lattice(viz,L); 
			break;
		case HEXAGONAL: honeycomb_lattice(viz,L);
			break;
		case KAGOME: kagome_lattice(viz,L);
			break;
		case TRIANGULAR: triangular_lattice(viz,L);
			break;
		case MOORE: moore_lattice(viz,L);
			break;
		/*case RANDOM: rede_random(viz,L);
			break;*/
		default:
			fprintf(stderr,"ERRO inicializacao rede!\n");
			fflush(stderr);
	} 

//-------------------C.I.-------------------------------------------------------------

	calculo_ci_estado(state,viz,investimento);
	calculo_ci_topologia(topologia);
	#ifdef dens_aleat
	calculo_ci_investimento(investimento,state,gama);
	#endif
	
	calculo_densidades(state,investimento,t,fp);

	#ifdef histograma_plot
	calculo_histograma(t,investimento,histograma,state);
	#endif
//------------------M.C.S.------------------------------------------------------------

	for(t=1; t < tmax; t++)
	{	

		#ifdef snapshot_gnuplot	
		snap_gnuplot(state, topologia,investimento,label,t);
		#endif 
		#ifdef snaphot_gif	
		snap_gif(state, topologia,t);
		#endif 
		#ifdef snapshot_hexa	
		snap_hexa(state, topologia,t);
		#endif 
		#ifdef snapshot_kagome
		snap_kagome(state, topologia,t);
		#endif 


		/*if(t==1 || t==40 || t==90 || t==120 || t==300){
			printf("set terminal postscript eps enhanced color\n");
			printf("set output 'snapshotC1_MCS%d.eps'\n",t);
			snap_gnuplot(state, topologia,investimento,label,t);
		}*/

//------------------------------------------------------------------------------------

		calculo_mcs(payoff,state,viz,r,gama,delta,topologia,investimento,commitment,t);
		/*if(t>70000)
		{
		calculo_percolacao(state,viz,label);	
		calculo_tamanho_cluster(state,viz,label,investimento,t);
		}*/
		//calculo_mobilidade(state,viz);
		
		//printf("%d %d %d %d %d %d %d %d %d %d %d %d \n",t, histograma[-5], histograma[-4], histograma[-3], histograma[-2], histograma[-1], histograma[0], histograma[1], histograma[2], histograma[3], histograma[4], histograma[5]);

//------------------------------------------------------------------------------------

		if (t >= measure_time)
		{
			calculo_densidades(state,investimento,t,fp);
			
			#ifdef histograma_plot
			calculo_histograma(t,investimento,histograma,state);
			#endif

			measure_time = rint(1.1*t); //printf("%d %lf\n",t, measure_time);
		}    
  		if((NC/(N)==1) || (NP/(N)==1) || (NPC/(N)==1) || (ND/(N)==1) )
		{
			do
			{
				if ( t >= measure_time)
				{
					calculo_densidades(state,investimento,t,fp);
					measure_time = rint(1.1*measure_time);
				}
				t=t+1;
			}
			while(t<tmax);
		 break;
		}
	
	}//MCS

fclose(fp);

return 0;
free_2d_int_pointer(viz,N,G);

} //main
