 //author : Zarrouk Rim
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include<time.h>
#include <limits.h>
#include <stdlib.h>
#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) < (b) ? a : b)
#define SWARM_SIZE 500
#define MAX_ITERATION 500
#define maxswarmPSO2 10
#define maxiterPSO2 10

#define nbPatient    12
#define nbTech      2
int MaxbestFitness = 300;
int Ll[nbPatient][nbTech+1];
int deplacement [nbPatient+1][nbPatient+1] = {
{0, 20, 20,20, 20, 20, 20, 20, 20, 20, 20, 20 ,20},
{20,0, 28, 25,26, 20, 14, 23, 31, 23, 20, 22, 24 },
{20,28, 0, 27, 39, 37, 25, 25, 23, 27, 40, 26, 36},
{20,25, 27, 0, 33, 34, 22, 12, 20, 32, 37, 14, 33},
{20,26, 39, 33, 0, 24, 23, 34, 42, 38, 39, 33, 39},
{20,20, 37, 34, 24, 0, 24, 32, 40, 33, 18, 31, 34},
{20,14, 25, 22, 23, 24, 0, 20, 28, 26, 27, 19, 27},
{20,23, 25, 12, 34, 32, 20, 0, 8, 30, 34, 9, 31},
{20,31, 23, 20, 42, 40, 28, 8, 0, 38, 42, 17, 39},
{20,23, 27, 32, 38, 33, 26, 30, 38, 0, 29, 30, 14},
{20, 40, 37, 39, 18, 27, 34, 42, 29, 0, 34, 31},
{20,22, 26, 14, 33, 31, 19, 9, 17, 30, 34, 0, 31},
{20,24, 36, 33, 39, 34, 27, 31, 39, 14, 31, 31, 0}
};

double C1 , C2;
double W_UPPERBOUND =1.2;
double W_LOWERBOUND = 0.2;
double w , wPSO2;
double r1,r2;
double fitness_avg ;
int fitness_min ;
int tab_fitness[SWARM_SIZE] ;
int gBestScheduling [2 * nbPatient + 2*nbTech][3];
int newLocation[2*nbPatient];

typedef struct {
  int PreviousFitness;
  int bestFitness;
  int PreviousLocation[2*nbPatient];
  int bestLocation[2*nbPatient];
  double PreviousVelocity[2*nbPatient];
  int Scheduling[nbPatient + 2*nbTech][3];
  int bestIteration;
  int num_Particule;
} Particule;

 Particule bestParticule;
 Particule swarm[SWARM_SIZE];
 int vecteur_initial[nbPatient] ;
int iteris;
/****************************************************************************** Step_zero **/
void  Step_zero(int i, Particule*  p,   int* tab_deplacement[2] [nbTech ], int* position[2*nbPatient],  double* vitessee[2*nbPatient],   int* tab_Scheduling [nbPatient + 2*nbTech][3]) {
  int solution[2*nbPatient];
   double vitesse[2*nbPatient];
 int tab_depl[2] [nbTech ] ;
//printf("\nappel ApproachByRandom\n");
            ApproachByRandom(&solution);
/*
printf("\nappel ApproachByRandomVitesse\n");
            //ApproachByRandomVitesse(&vitesse);
            */
  //       printf("\nappel Scheduling\n");
         Scheduling(&tab_depl,solution , &tab_Scheduling);
    //                printf("\n  fin     Scheduling\n");
    for (int x= 0; x < nbTech; x++) {
			tab_deplacement [0][x]  = tab_depl[0][x] ;
			tab_deplacement [1][x]  = tab_depl[1][x] ;
		}
   for (int x= 0; x < 2*nbPatient; x++) {
			position [x]  = solution[x] ;
			vitessee [x]  =  (int) vitesse[x] ;

		}
    }// end step_zero

/****************************************************************************** initializeSwarm **/
void initializeSwarm(){
    int i,j;
    int position[2*nbPatient];
    double vitessee[2*nbPatient];
      int tab_deplacement[2] [nbTech ];
    int tab_Scheduling [nbPatient + 2*nbTech][3];
      bestParticule.bestFitness = 999999;
    for (i = 0; i < SWARM_SIZE; i++) {
        Particule p ;
        p.num_Particule=i;

 //      printf("\nappel Step_zero\n");
        Step_zero(i,&p, &tab_deplacement, &position , &vitessee , tab_Scheduling);

        int y;
    for (y = 0; y < 2*nbPatient; y++) {
            vitessee[y] =1;//rand()% 2;
         }

       evaluation (&tab_deplacement,position,vitessee , 0, &p, tab_Scheduling, i ); // remplir les élt de la particule
 swarm[i]=p;

		}
    printf("\n FIN initializeSwarm\n");
}// end initializeSwarm
/**************************************************************************** ApproachByRandom **/
int ApproachByRandom(int* position_init[2*nbPatient]) {
int i,j;

		for (i = 0; i < nbPatient; i++) {
			if  (i<nbPatient/2) position_init[i] = 1;
            else position_init[i] = 2;
		}
int occup[nbPatient];
			for (i = 0; i < nbPatient; i++) {
			      occup[i] = 0;
			// printf("  %d  \n" , occup[i] );
			}

int patient =0;
	for (i = 0; i < nbPatient ; i++) {
//printf("  %d  \n" , i);
   do{
   patient = rand()% nbPatient  +1;
   }
    while (occup[patient -1] ==1  );
     //  printf("  %d  \n" , patient-1);
       occup[patient-1] = 1;
       position_init[i+nbPatient] = patient ;

}

/**affichage**/
/*
		printf("Solution au niveau ApproachByRandom\n");
		 for (i = 0; i < 2; i++) {
                for (j = 0; j < nbPatient; j++) {
                 printf(" %d  ",  position_init[i][j]);
            }
		 }
            printf("\n");
*/
}//end ApproachByRandom

/************************************************************************ ApproachByRandomVitesse **/

void ApproachByRandomVitesse(double*  vit[2*nbPatient]) {
    int y;
    for (y = 0; y <2* nbPatient; y++) {
            vit[y] =rand()% 2;
         }
/**affichage**/

    printf("Vitesse au niveau  ApproachByRandomVitesse\n");

                for (int j = 0; j < 2*nbPatient; j++) {
                 printf(" %f  ",  vit[j]);

		 }

}//end ApproachByRandomVitesse

/***************************************************** iterate **/
void iterate(){
    int i,y,h , j  ;
 //srand(time(NULL)); // initialisation de rand
//  BestMachineParticle.MbestFitness = 9999999 ;
while (iteris < MAX_ITERATION  ) {
         printf("ttttttttttttttttt %d ----> %d\n",iteris , bestParticule.bestFitness);
//if (iteris == 0)
//{
    	for (i = 0; i < SWARM_SIZE; i++) {
                 // 	w = W_UPPERBOUND - iteris * (W_UPPERBOUND - W_LOWERBOUND) / MAX_ITERATION;
            //      C1=0.9 ; C2=0.6;
		 w = W_UPPERBOUND -(W_UPPERBOUND - W_LOWERBOUND) * (iteris-1) / (MAX_ITERATION-1);
	//	 printf("%f  ", w);
       C1 = 0.5 * sqrt(w+1);
  //       printf("%f  ", C1);
        C2 = MIN(4.2 , 2*(W_LOWERBOUND+1)) -0.5*sqrt(w+1)- 0.000001 ;
    //     printf("%f  ", C2);

              UpdateLocationVelocity(iteris,&swarm[i]  , i);

			}
			iteris++;
}
//}
     printf("\n Le Best est trouve a l'iteration %d :\n", bestParticule.bestIteration);
		printf(")\t Fitness %d \n",bestParticule.bestFitness);
		printf("bestParticule est %d \n", bestParticule.num_Particule);
	printf("\n\nPlaning table : \n    Staff|  Pat |  Ttravel | Tarrived | Tleave \n");
		for (y = 0; y < nbPatient + 2*nbTech; y++) {
		    printf("\n");
			for (h = 0; h < 5; h++) {
					printf(" %d         ",bestParticule.Scheduling[y][h]);
            }
		}

}// end iterate
void mappingCorrection(int fuzzyLocation [2*nbPatient])
{
        int ratte = 0;
        int occup[nbPatient];
        int patient =0, tech =0;

        int c =0 ;
        int occup_tec [nbTech];
           for ( int t =0 ; t<nbTech; t++)
        { occup_tec [t] =0;
        }
         /**** supprimer les doublons de techniciens**********/

            for (int i = 0; i < nbPatient ; i++) {
 //printf(" %d  ",  fuzzyLocation[i]);
            if (fuzzyLocation [i]  == 0)
             {  //   printf(" je suis if  0  \n");
                 fuzzyLocation [i]  = 999;
            }
            else
            {if (fuzzyLocation [i] >nbTech)
                       {fuzzyLocation [i]  = 999;
                          //   printf(" je suis if > nbTech \n");
                          }

              else {
                   //     printf(" je suis else  \n");
                for ( int t =0 ; t<nbTech; t++)
                { // c =0 ;
                 if (fuzzyLocation [i] == t +1 )
                 { if (  occup_tec[t] <  nbPatient/2)
                    {
                        occup_tec[t] =  occup_tec[t] +  1;
                        fuzzyLocation [i]  = t+1;
                    }
                        else  {fuzzyLocation [i]  = 999;}
                 }

            }
              }
            }
             //  printf("%d occup_tec[t-1] \n",   occup_tec[t]);
        }

      /**** Correction des doublons de techniciens**********/
        int infir = 0;
        for ( int  i=0 ; i< nbPatient; i++)
        {  if   (fuzzyLocation [i]  == 999)
            { // printf(" je suis uzzyLocation [i]    technixc  %d\n" , fuzzyLocation [i]);
                  do{
                 //       printf("%d ",  occup_tec[tech -1]);
                      infir = rand()% nbTech  ;
                 //     printf("%d  rand ", patient);
                    }  while (occup_tec[infir] >= nbPatient/nbTech  );
           //   printf("  %d  \n" , infir);
               occup_tec[infir] ++;
               fuzzyLocation [i]  =  infir +1;
            }
        }

         /**** supprimer les doublons de patients**********/
for (int i = 0; i < nbPatient; i++) {
                          occup[i] = 0;
                    }

     for (int i = nbPatient; i < nbPatient*2 ; i++) {
  // printf(" %d  \n",  fuzzyLocation[i]);
                    if (abs(fuzzyLocation [i]) < 1)
                    {
                        fuzzyLocation [i] = fuzzyLocation [i]*10;
          //   printf(" je suis if < 1 \n");
                    }
           if (abs( fuzzyLocation [i] ) == 0)
             {
                 fuzzyLocation [i]  = 999;
            }
            else if ( abs( fuzzyLocation [i] )  > nbPatient)
                { fuzzyLocation [i]  = 999;
                }
          else
          {
                for ( int t =1 ; t<=nbPatient; t++)
                {
                if (abs( fuzzyLocation [i] ) == t  )
                 {  //     printf(" je suis == t+1   %d  \n"  , occup [t-1]);
                     if (occup [t-1]== 0 )
                        {occup [t-1] = 1 ;
                            fuzzyLocation [i]  = t ;
                            // printf(" mabadaltach  %d   \n",  fuzzyLocation [i] );
                       }
                     else  {fuzzyLocation [i]  = 999;
                    //   printf(" badalt  %d   \n  *************** \n",  fuzzyLocation [i] );
                    }
                 }

                 }
          }
     }//for

        /**** correction les doublons de patients**********/
        for (int i = nbPatient; i <= nbPatient*2 ; i++) {
              //  printf("  je suis iiiiiiiiiiiiii %d   je suis fuzzy loc %d  \n", i , fuzzyLocation [i]  );
            if   (fuzzyLocation [i]  == 999)
            {
                  do{
                        patient = rand()% nbPatient ;
                      }    while (occup[patient ] ==1  );

               occup[patient] = 1;
               fuzzyLocation [i]  = patient+1 ;
               //    printf(" %d    %d \n",  occup[patient],    fuzzyLocation [i] );
            }
             }
// printf("  lahnaaaaaaaaaaaaaaaaaaaaaaaaaaa \n" );

             for (int i = 0; i < nbPatient*2 ; i++) {
                   newLocation[i] =  fuzzyLocation [i];
             }

/*
        printf("\n ******************Solution au niveau mappingCorrection  \n");

                        for (int j = 0; j < nbPatient *2; j++) {
                         printf(" %d ",  newLocation[j]);
                    }

                    printf("********************\n");*/
}

/*********************************************************************************** Evaluation **/
void evaluation(int* tab_deplacement [2][nbTech], int solution[2*nbPatient], double vitesse[2*nbPatient] , int iteration, Particule* p ,   int tab_Scheduling [nbPatient + 2*nbTech][3], int num_p) {
int fitness =0 , x ; // la fonction Z

	//	p->bestFitness = 999999;
	//	p->PreviousFitness = 999999;
        int coordonne [nbPatient*2] ;
        double coordonne_vitesse [nbPatient*2] ;
        int y =nbPatient;

         y =0;
        for (x= 0; x < nbTech; x++) {
             fitness =  MAX(fitness ,  tab_deplacement[1] [x])  ;
        }
    //    	printf( "Fitnesssssss %d  particule  %d    iteration %d \n ", fitness,   p->num_Particule , iteris);

  p->PreviousFitness = fitness ; //999999;

//printf( "\n Fitnesssssss  p->PreviousFitness  %f   \n",p->PreviousFitness);
 for (x= 0; x < nbPatient*2; x++) {
		p->PreviousLocation[x]=   solution[x];
    		p->PreviousVelocity[x] = vitesse[x];
		}

 for (x= 0; x < nbPatient + 2*nbTech ; x++) {
     for (y = 0; y< 3; y++)
     {
          p->Scheduling [x][y] = tab_Scheduling [x][y] ;
     }
 }
		     for (int i = 0; i < nbPatient+2*nbTech; i++){
                for (int j = 0; j < 3; j++){
           p->Scheduling[i][j] = tab_Scheduling [i][j];
                }//end for j
            }//end for i
if (iteration==0) {
              p->bestFitness = fitness ; //999999;
		    bestParticule.bestIteration=0;
			bestParticule.num_Particule=0;
			bestParticule.PreviousFitness = p->PreviousFitness;
   for (x= 0; x < nbPatient*2; x++) {
            bestParticule.bestLocation[x] = 0;
			bestParticule.PreviousLocation[x]=  0;
			bestParticule.PreviousVelocity[x] =  0.0;
		}}
	              if (fitness < p->PreviousFitness)
            {
                  p->bestFitness = fitness ; //999999;
		          p->bestIteration = iteration;
		        p->  num_Particule=num_p;
		           for (int i = 0; i < 2*nbPatient; i++){
                        p->bestLocation[i] =   p->PreviousLocation[i];
               // printf(" %f" ,  p.bestLocation[i]) ;
        }
         for (int i = 0; i < nbPatient+2*nbTech; i++){
                for (int j = 0; j < 3; j++){
                    p->Scheduling[i][j] = tab_Scheduling[i][j] ;
                }//end for j
            }//end for i
		}
}

/*********************************************************************************** Scheduling **/
void Scheduling(int *tab_deplacement [2][nbTech],   int solution[2*nbPatient],      int* tab_Scheduling [nbPatient + 2*nbTech][3]) {
        int y,k;
        int tabTechnicien [nbPatient/nbTech];

			for (y = 0; y < nbPatient + 2*nbTech ; y++) {
			tab_Scheduling[y][0] = 0; //Technicien
			tab_Scheduling[y][1] = 0; //Patient
            tab_Scheduling[y][2] = 0; //temp
			}
	int tech =1;
	int deb =0 ;
	int fin =(nbPatient/nbTech -1)+2;
	int indimm =0;

         // on remplit  dans tab_Scheduling
		for (y = deb; y <= fin; y++) {
			tab_Scheduling[y][0] = tech; //Technicien
			tab_Scheduling[y][1] = 0; //Patient
            tab_Scheduling[y][2] = 0; //temp
			if (y == fin && y <= (nbPatient -(nbPatient/nbTech)+1))
            {   deb = fin;
			   fin = fin+ (nbPatient/nbTech) +2;
			     tech =tech+1;
			}
		}

int indice =0;

for (int t=1; t<=nbTech; t++)
{  tab_Scheduling[indice] [1] =0 ;
    indice ++;
  //  printf(" %d",t);
    for (int i =0; i<nbPatient   ;i++)
    {

        if (solution[i] == t)
        {
            tab_Scheduling[indice] [1] =solution[i+nbPatient] ;
            //       printf(" \n solution[i+nbPatient] %d  ", solution[i+nbPatient] );
            indice ++;
        }
    }
    tab_Scheduling[indice] [1] =0 ;
    indice ++;
}

/**********Calcule de temps de deplacement*************/
 tech =1;
 deb =0 ;
 fin =nbPatient/nbTech +1;
int depart = 0, arrive =0, somme_dep=0;
		for (y = deb; y <= fin; y++) {
            //    printf("y  if  %d\n",  tech);
        arrive = tab_Scheduling[y][1];
		tab_Scheduling[y][2]  =  deplacement [arrive] [depart]	; //temps de deplacement
        depart =  arrive ;
            if (y == fin && y <= (nbPatient -(nbPatient/nbTech)+2))
            {
                   deb = fin;
			     fin = fin+ (nbPatient/nbTech)+2;
			     tech =tech+1;
			}
		}

int mk , ff,  pppp =0;
    for (y = 0; y <nbTech; y++) {
          //  printf("tech num %d\n", y+1);
            somme_dep = 0;
             mk =0; ff=0 ;
            tab_deplacement [0] [y]  = y+1;
            tab_deplacement [1] [y] = 0;
	for (int i=0 ; i<(nbPatient + 2*nbTech); i++){
	       if(tab_Scheduling[i][0] == y+1)
	       {
	            mk =  tab_Scheduling[i][2]   ;
	  //     printf("\n je suis mk   %d\n", mk);
         ff = tab_deplacement [1] [y] ;
      //  printf("\n je suis ff   %d\n",  tab_deplacement [1] [y] );
         tab_deplacement [1] [y]  =  mk +ff    ;

	       }
	}
	//printf("\n je suis tab_deplacement   %d\n",  tab_deplacement [1] [y] );
	 pppp =  tab_deplacement [1] [y] ;
	tab_deplacement [1] [y]  = pppp   + 20*(nbPatient/nbTech);
	//	printf("\n je suis tab_deplacement ****  %d\n",  tab_deplacement [1] [y] );
	}

}//end Scheduling

/*************** ***************************** UpdateLocationVelocity **/
void UpdateLocationVelocity (int t,Particule* pp , int ind){

    double newVelocity [nbPatient*2];
    int fuzzyLocation [nbPatient*2];
    int tab_Scheduling [nbPatient + 2*nbTech][3];
    int i,j,y,x,bb;
    Particule p;
    p.bestFitness = pp->bestFitness;
    p.PreviousFitness= pp->PreviousFitness;
    p.num_Particule= pp->num_Particule;
    p.bestIteration = pp->bestIteration;
    for (int j = 0; j < nbPatient *2; j++) {
        p.PreviousVelocity[j] = pp->PreviousVelocity[j];
        p.PreviousLocation[j] = pp->PreviousLocation[j];
        p.bestLocation[j]= pp->bestLocation[j];
    }
//update velocity

   for (x = 0; x < nbPatient*2; x++) {
     r1 = (double)rand() / (double)RAND_MAX;
     r2 = (double)rand() / (double)RAND_MAX;
       //   double c= ;
    //   printf(" \n %f      %f    %f    %f    %f",  r1, r2, C1, C2, w );
   //   printf(" \n %f  ",  w );
         p.PreviousVelocity [x] = (w *  (double) p.PreviousVelocity[x]
                + r1  * C1 *(double) (p.bestLocation [x] - p.PreviousLocation[x])) +r2 * C2*(double)(bestParticule.bestLocation[x]- p.PreviousLocation[x]);
        //   printf(" \n %f     ",  p.PreviousVelocity[x] );
  }

 // update position -. fuzzyLocation [nbPatient][2]
  for (x = 0; x < nbPatient*2; x++) {
//        printf(" \n %f     %f     ", p.PreviousLocation[x], p.PreviousVelocity[x]);
   fuzzyLocation[x] = abs( p.PreviousLocation[x] + p.PreviousVelocity [x] );
  }

    mappingCorrection(fuzzyLocation); //return newLocation[nbPatient][2]

  for (x = 0; x < nbPatient*2; x++) {
   p.PreviousLocation[x]= newLocation[x];
 //  printf("%d  \t",(int)newLocation[x] );
  }

  //(int *tab_deplacement [2][nbTech], int solution[2][nbPatient], int vitesse[2][nbPatient] , int i, Particule *p
int tab_deplacement [2][nbTech];

  Scheduling(tab_deplacement, newLocation, &tab_Scheduling) ;
 //printf("\n fin Scheduling  \n");
    evaluation (&tab_deplacement,newLocation,p.PreviousVelocity , t, &p, &tab_Scheduling, ind );

           if ( p.bestFitness  <  bestParticule.bestFitness) {
		    bestParticule.bestIteration=iteris;
			bestParticule.num_Particule=ind;
			bestParticule.PreviousFitness = p.PreviousFitness;
			bestParticule.bestFitness = p.bestFitness  ;
        for (x= 0; x < nbPatient*2; x++) {
            bestParticule.bestLocation[x] = p.bestLocation[x] ;
			bestParticule.PreviousLocation[x]=  p.PreviousLocation[x];
			bestParticule.PreviousVelocity[x] =  p.PreviousVelocity[x] ;
           for (int i = 0; i < nbPatient+2*nbTech; i++){
                for (int j = 0; j < 3; j++){
                    bestParticule.Scheduling[i][j] = p.Scheduling[i][j] ;
                }//end for j
            }//end for i
 }
		}
        pp->bestFitness = p.bestFitness;
        pp->PreviousFitness= p.PreviousFitness;
        pp->num_Particule=  p.num_Particule;
        pp->bestIteration =  p.bestIteration;
        for (int j = 0; j < nbPatient *2; j++) {
        pp->PreviousVelocity[j] =  p.PreviousVelocity[j];
        pp->PreviousLocation[j] =  p.PreviousLocation[j];
        pp->bestLocation[j]=  p.bestLocation[j];
}
} //end UpdateLocationVelocity

/**************************************************************************** main **/
int main(){
   //printf("appel initializeSwarm\n");
    initializeSwarm();
    //printf("fin appel initializeSwarm\n");
    iterate();
   // getch();

}

