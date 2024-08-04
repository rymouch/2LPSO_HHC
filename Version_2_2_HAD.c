 //author : Zarrouk Rim
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include<time.h>
#include <limits.h>
#include <stdlib.h>
#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) < (b) ? a : b)
#define SWARM_SIZE  1000
#define MAX_ITERATION 600
#  include <windows.h>
#  define psleep(sec) Sleep ((sec) * 1000)
#define nbP    30
#define nbPatient    12
#define nbTech      2
int MaxbestFitness = 300;
int Ll[nbPatient][nbTech+1];
int deplacement [nbPatient+1][nbPatient+1] = {
{0   , 20  , 20   ,20   , 20   , 20    , 20    , 20    , 20     , 20     , 20   , 20      ,20    },
{20 ,0     , 28   , 25  ,26    , 20    , 14    , 23    , 31     , 23     , 20   , 22       , 24   },
{20 ,28   , 0,    27    , 39    , 37    , 25    , 25    , 23     , 27    , 40   , 26       , 36  },
{20 ,25   , 27  , 0     , 33    , 34    , 22    , 12     , 20     , 32    , 37   , 14      , 33  },
{20 ,26   , 39  , 33    , 0     , 24    , 23    , 34     , 42     , 38    , 39    , 33      , 39 },
{20 ,20   , 37  , 34    , 24   , 0       , 24    , 32     , 40    , 33    , 18    , 31     , 34  },
{20 ,14   , 25  , 22    , 23    , 24    , 0      , 20    , 28     , 26    , 27     , 19     , 27 },
{20 ,23   , 25  , 12    , 34     , 32    , 20    , 0     , 8       , 30    , 34     , 9       , 31 },
{20 ,31   , 23  , 20    , 42     , 40     , 28    , 8    , 0       , 38    , 42     , 17      , 39},
{20 ,23   , 27  , 32     , 38    , 33     , 26    , 30    , 38    , 0      , 29    , 30      , 14},
{20 , 20    ,40  , 37  , 39    , 18     , 27    , 34      , 42    , 29    ,0      , 34     , 31  },
{20 ,22   , 26  , 14    , 33     , 31    , 19     , 9       , 17    , 30     , 34    , 0      , 31},
{20 ,24   , 36   , 33   , 39     , 34    , 27     , 31     , 39    , 14     , 31     , 31    , 0 }
};
#define nbPatient_rq 1
 int new_deplacement [nbPatient+1+nbPatient_rq][nbPatient+1 +nbPatient_rq]  = {
{0   , 20  , 20   ,20   , 20   , 20    , 20    , 20    , 20     , 20     , 20   , 20      ,20      , 20},
{20 ,0     , 28   , 25  ,26    , 20    , 14    , 23    , 31     , 23     , 20   , 22       , 24     , 10 },
{20 ,28   , 0,    27    , 39    , 37    , 25    , 25    , 23     , 27    , 40   , 26       , 36     , 20},
{20 ,25   , 27  , 0     , 33    , 34    , 22    , 12     , 20     , 32    , 37   , 14      , 33     , 15},
{20 ,26   , 39  , 33    , 0     , 24    , 23    , 34     , 42     , 38    , 39    , 33      , 39    , 10},
{20 ,20   , 37  , 34    , 24   , 0       , 24    , 32     , 40    , 33    , 18    , 31     , 34     , 30},
{20 ,14   , 25  , 22    , 23    , 24    , 0      , 20    , 28     , 26    , 27     , 19     , 27     , 15},
{20 ,23   , 25  , 12    , 34     , 32    , 20    , 0     , 8       , 30    , 34     , 9       , 31     , 39},
{20 ,31   , 23  , 20    , 42     , 40     , 28    , 8    , 0       , 38    , 42     , 17      , 39     , 40 },
{20 ,23   , 27  , 32     , 38    , 33     , 26    , 30    , 38    , 0      , 29    , 30      , 14     , 20},
{20 , 20    ,40  , 37  , 39    , 18     , 27    , 34      , 42    , 29    ,0      , 34     , 31     , 19},
{20 ,22   , 26  , 14    , 33     , 31    , 19     , 9       , 17    , 30     , 34    , 0      , 31     , 10},
{20 ,24   , 36   , 33   , 39     , 34    , 27     , 31     , 39    , 14     , 31     , 31    , 0      ,  17},
{20 ,10   , 20   , 15    , 10    , 30    , 15     , 39     , 40    , 20     , 19     , 10    , 17    , 0 }
};
 int occup[nbPatient+1];
int occup_tec [nbTech];
double C1 , C2;
double W_UPPERBOUND =1.2;
double W_LOWERBOUND = 0.2;
double w , wPSO2;
double r1,r2;
int nbPatient_annul = 0 ;
int nbPatientLocal = nbPatient;
int gBestScheduling [nbPatient +1 + 2*nbTech][5];
int newLocation[2*(nbPatient+1)];
  int patient_annule =0 ;
   int patient_new =0 ;
typedef struct {
  int PreviousFitness;
  int bestFitness;
  int PreviousLocation[2*nbP ];
  int bestLocation[2*nbP];
  double PreviousVelocity[2*nbP];
  int Scheduling[nbP + 2*nbTech][5];
  int bestIteration;
  int num_Particule;
} Particule;
int annul_p[nbP];
int annul_t[nbTech];
 Particule bestParticule;
Particule secondbestParticule;
Particule ThirdbestParticule;
 Particule swarm[SWARM_SIZE];
Particule Secondswarm[SWARM_SIZE];
 int vecteur_initial[nbP] ;
int iteris=0;
int rest_patient = 0;
   int pt_depart[nbTech][5];
      int newNbPatient = nbPatient  ;
      int solution1[2*(nbPatient+1)] ;
/****************************************************************************** initializeSwarm **/
void initializeSwarm(Particule* pbest , int choix ){
    int i,j;
if (choix == 1)
     nbPatientLocal  = nbPatient- nbPatient_annul ;
 if (choix == 2)
    nbPatientLocal  = nbPatient + nbPatient_rq ;

   //  printf("%d ", nbPatientLocal) ;
    int position[2*nbPatientLocal];
    double vitessee[2*nbPatientLocal];
    int tab_Scheduling [nbPatientLocal + 2*nbTech][5];
      pbest->bestFitness = 999999;

    for (i = 0; i < SWARM_SIZE; i++) {
        Particule p ;
         swarm[i]=p;
        p.bestFitness = 99999;
        p.num_Particule=i;
 //      printf("\nappel Step_zero\n");
      //  Step_zero(i,&p, &tab_deplacement, &position , &vitessee , &tab_Scheduling);
//printf("\nappel ApproachByRandom\n");
            ApproachByRandom(&position);

//printf("\nappel ApproachByRandomVitesse\n");
//            ApproachByRandomVitesse(&vitesse);

 //     printf("\nappel Scheduling\n");
         Scheduling(position , &tab_Scheduling , choix);
      //   printf("\n  fin     Scheduling\n");
   // printf(" \n **********  Step_zero tab scheduling \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<=4 ; j++){
              p.Scheduling [i] [j] = tab_Scheduling[i][j];
            }
		 }
 //  printf("\nfin appel Step_zero\n");
    for (int y = 0; y < 2*nbPatientLocal; y++) {
           p.PreviousVelocity [y]= 1.0;//rand()% 2;
          p.PreviousLocation [y]= position[y];
         }

int     fitness= 0 ;
        for (int x= 0; x < nbPatientLocal + 2*nbTech; x++) {
            fitness = MAX(fitness, p.Scheduling [x] [3]);
               //printf( "Fitnesssssss   %d \n ", fitness);
        }
//   printf( "evaluation  iteration 0  Fitnesssssss   %d \n ", fitness);

  	            if (p.bestFitness >=   fitness)
            {     //printf( ":p  ***    p->bestFitness  %d      p->PreviousFitness,  %d   p->num_Particule  %d     iteration %d \n ", p->bestFitness  , p->PreviousFitness,   p->num_Particule , iteris);

                  p.bestFitness = fitness ; //999999;
		          p.bestIteration = 0;
                 p.num_Particule= i;
		           for (int ii= 0; ii < 2*nbPatientLocal; ii++){
                        p.bestLocation[ii] =   p.PreviousLocation[ii];
               // printf(" %f" ,  p.bestLocation[i]) ;
                for (int ii= 0; ii< nbPatientLocal+2*nbTech; ii++)
                {
                for (int j = 0; j < 4; j++){
                    p.Scheduling[ii][j] =tab_Scheduling[ii][j] ;
                }//end for j
            }
        }

            }

            pbest->bestIteration=0;
			pbest->num_Particule=0;
			pbest->PreviousFitness = p.PreviousFitness;
   for (int x= 0; x < nbPatientLocal*2; x++)
    {
            pbest->bestLocation[x] = 0;
			pbest->PreviousLocation[x]=  0;
			pbest->PreviousVelocity[x] =  0.0;
    }

           for (int ii= 0; ii< nbPatientLocal+2*nbTech; ii++){
                for (int j = 0; j < 4; j++){
                    pbest->Scheduling[ii][j] =0;
                }//end for j
            }//end for


 swarm[i]=p;

/*
		 for (int j = 0; j < nbPatientLocal *2; j++) {
                 printf(" \n %lf     %f",    p.bestLocation[j]  , p.PreviousLocation[j]);
            }
     printf("\n *************initializeSwarm  \n");
    for (int j = 0; j < nbPatientLocal *2; j++) {
                 printf(" \n %lf     %f",    p.bestLocation[j]  , p.PreviousLocation[j]);
            }
printf("\n *************  \n");
*/
		}
   printf("\n FIN initializeSwarm\n");
}// end initializeSwarm
/**************************************************************************** ApproachByRandom **/
int ApproachByRandom(int* position_init[2*(nbPatientLocal)]) {
int i,j;
		for (i = 0; i < nbPatientLocal; i++) {
			if  (i<nbPatientLocal/2) position_init[i] = 1;
            else position_init[i] = 2;
		}
//int occup[nbPatientLocal];
			for (i = 0; i < nbPatientLocal; i++) {
			      occup[i] = 0;
			// printf("  %d  \n" , occup[i] );
			}


int patient =0;
	for (i = 0; i < nbPatientLocal ; i++) {
//printf("  %d  \n" , i);
   do{
   patient = rand()% nbPatientLocal  +1;
   }
    while (occup[patient -1] ==1  );
     //  printf("  %d  \n" , patient-1);
       occup[patient-1] = 1;
       position_init[i+nbPatientLocal] = patient ;
    /* printf(" 11111 ApproachByRandom  %d  \n" ,  position_init[i]);
     for (int u = 0; u < nbPatientLocal ; u++) {
        printf("  %d  \n" , occup[u] );
     }*/
}

/**affichage**/
/*
		printf("Solution au niveau ApproachByRandom\n");
		 for (i = 0; i < 2; i++) {
                for (j = 0; j < nbPatientLocal; j++) {
                 printf(" %d  ",  position_init[i][j]);
            }
		 }
            printf("\n");
*/
}//end ApproachByRandom

/***************************************************** iterate **/
void iterate(int choix){
    int i,y,h , j  ;
 //srand(time(NULL)); // initialisation de rand
//  BestMachineParticle.MbestFitness = 9999999 ;
iteris = 0;
while (iteris < MAX_ITERATION  ) {
//printf("ttttttttttttttttt %d ----> %d\n",iteris , bestParticule.bestFitness);
printf(".");
 if (choix == 0)
 {
      for (i=0;i<nbPatientLocal;i++)
    {
        annul_p[i]=0;
    }
    for (i=0;i<nbTech;i++)
    {
        annul_t[i]=0;
    }

 }
 /*else{
        printf("\n dans itere  Patient déja traité ");
    for (y = 0; y < nbPatient ; y++) {
    printf( "\n %d  " ,   annul_p[y]);
}
printf("\n charge tech ");
for (y = 0; y < nbTech ; y++) {
    printf( "\n %d  " ,   annul_t[y]);
}
 }*/

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

    for (h=0;h<nbPatientLocal;h++)
    {
         occup[h] = annul_p[h];
    }
 for (h=0;h<nbTech;h++)
        {
            occup_tec [h] =annul_t[h];
        }
  //        if (choix == 2)  printf("\n  debut  UpdateLocationVelocity") ;

                  if (choix == 2)
                  {
                      UpdateLocationVelocity(iteris,&Secondswarm[i]  , i  , choix);
                      }

                  else  {
                        UpdateLocationVelocity(iteris,&swarm[i]  , i  , choix);
                  }

			}
			iteris++;
}
//}
  if (choix == 2)
        {
                printf("\n Le Best est trouve a l'iteration %d :\n", ThirdbestParticule.bestIteration);
		printf(")\t Fitness %d \n",ThirdbestParticule.bestFitness);

		printf("bestParticule est %d \n", ThirdbestParticule.num_Particule);

		printf("\n\nPlaning table : \n    Staff|  Pat |  Ttravel | Tarrived | Tleave \n");
		for (y = 0; y < newNbPatient + 2*nbTech; y++) {
		    printf("\n");
			for (h = 0; h < 5; h++) {
					printf(" %d         ",ThirdbestParticule.Scheduling[y][h]);
            }
		}
    }
else{
     printf("\n Le Best est trouve a l'iteration %d :\n", bestParticule.bestIteration);
		printf(")\t Fitness %d \n",bestParticule.bestFitness);

		printf("bestParticule est %d \n", bestParticule.num_Particule);

		printf("\n\nPlaning table : \n    Staff|  Pat |  Ttravel | Tarrived | Tleave \n");
		for (y = 0; y < nbPatientLocal + 2*nbTech; y++) {
		    printf("\n");
			for (h = 0; h < 5; h++) {
					printf(" %d         ",bestParticule.Scheduling[y][h]);
            }
		}
}
}// end iterate
/***************************************************** mappingCorrection **/
void mappingCorrection(int fuzzyLocation [2*(nbPatientLocal)])
{ //  printf("\n debut mapping correction \n") ;
      // printf("\n ");
            /*            for (int j = 0; j < nbPatient *2; j++) {
                        printf(" %d  ",  fuzzyLocation[j]);
                    }*/
                    int moitier = 0 ;
        if (nbPatientLocal%2 !=0) moitier = (nbPatientLocal/2)+1 ;
                   //  printf(" %d  ",  nbPatient_annul);
        int ratte = 0;
     //   int nbPatientLocal  = nbPatient - nbPatient_annul ;
        int patient =0, tech =0;
        int c =0 ;

         /**** supprimer les doublons de techniciens**********/

            for (int i = 0; i < nbPatientLocal ; i++) {
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
                 { if (  occup_tec[t] <= moitier)
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
 /*      printf("\n ");
                        for (int j = 0; j < nbPatientLocal *2; j++) {
                        printf(" %d  ",  fuzzyLocation[j]);
                    }
      /**** Correction des doublons de techniciens**********/
        int infir = 0;
        for ( int  i=0 ; i< nbPatientLocal; i++)
        {  if   (fuzzyLocation [i]  == 999)
            { // printf(" je suis uzzyLocation [i]    technixc  %d\n" , fuzzyLocation [i]);
                  do{
                 //       printf("%d ",  occup_tec[tech -1]);
                      infir = rand()% nbTech ;
                   // printf("%d  rand ", infir);
                    }  while (occup_tec[infir] >= nbPatientLocal/nbTech +1 );
             //printf("  %d  \n" , infir);
               occup_tec[infir] ++;
               fuzzyLocation [i]  =  infir +1;
            }
        }
//printf("11111");
         /**** supprimer les doublons des patients**********/

     for (int i = nbPatientLocal; i < nbPatientLocal*2 ; i++) {
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
            else if ( abs( fuzzyLocation [i] )  > nbPatientLocal)
                { fuzzyLocation [i]  = 999;
                }
          else
          {
                for ( int t =1 ; t<=nbPatientLocal; t++)
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

       /*     if (fuzzyLocation [i] == patient_annule)
                       {fuzzyLocation [i]  = 999;

                          //   printf(" je suis if > nbTech \n");
                          } */
     }//for
 //   printf("9999911111");
  /*      printf("\n ");
                        for (int j = 0; j < nbPatientLocal *2; j++) {
                        printf(" %d  ",  fuzzyLocation[j]);
                    }
 printf("   \n  occup teck \n");
                        for (int j = 0; j < nbPatient ; j++) {
                        printf(" %d   \n",  occup[j]);
                    }
 printf("\n ");
        /**** correction les doublons de patients**********/

        for (int i = nbPatientLocal;  i <=nbPatientLocal*2 ; i++) {
              //  printf("  je suis iiiiiiiiiiiiii %d   je suis fuzzy loc %d  \n", i , fuzzyLocation [i]  );
            if   (fuzzyLocation [i]  == 999)
            {
                  do{
                        patient = rand()% nbPatientLocal ;
                        //printf("%d  rand ", patient);
                      }while (occup[patient ] ==1  );

               occup[patient] = 1;
               fuzzyLocation [i]  = patient+1 ;
          //        printf(" %d    %d \n",  occup[patient],    fuzzyLocation [i] );
            }
             }
// printf("  lahnaaaaaaaaaaaaaaaaaaaaaaaaaaa \n" );

             for (int i = 0; i < nbPatientLocal*2 ; i++) {
                   newLocation[i] =  fuzzyLocation [i];
             }
   // printf("fin mapping cor");

   /*     printf("\n ******************Solution au niveau mappingCorrection  \n");

                        for (int j = 0; j < nbPatientLocal *2; j++) {
                         printf(" %d ",  newLocation[j]);
                    }

                    printf("********************\n");
                    */
}



/***************************************************** mappingCorrection **/
void mappingCorrection2(int fuzzyLocation [2*newNbPatient])
{ //  printf("\n debut mapping correction \n") ;
      // printf("\n ");

                    int moitier = 0 ;
        if (nbPatientLocal%2 !=0) moitier = (nbPatientLocal/2)+1 ;
        else moitier = nbPatientLocal/2;
                   //  printf(" %d  ",  nbPatient_annul);
        int ratte = 0;
     //   int nbPatientLocal  = nbPatient - nbPatient_annul ;
        int patient =0, tech =0;
        int c =0 ;

 for (int y = 0; y < nbPatientLocal *2+ 2; y++) {
     solution1[y] =0;
 }
           /**** Solution 1 preparing*********/
for (int h=0 ; h<nbTech ; h++)
{
    for (int y = 0; y < nbPatient+ 2*nbTech; y++) {
         if ( secondbestParticule.Scheduling[y] [0] == h+1 && secondbestParticule.Scheduling[y] [1] !=0)
        {  // printf("\n   bestParticule.Scheduling[y] [3]   %d  "   , bestParticule.Scheduling[y] [3]);
            if(secondbestParticule.Scheduling[y] [4] <=110 )
                {

                      //  printf(" bestParticule %d", bestParticule.Scheduling[y] [1]) ;
                        solution1[y+nbPatientLocal-1] = secondbestParticule.Scheduling[y] [1] ;
                        solution1[y-1] = h +1 ;
                       // printf("\n solution %d    %d ", solution [i]  , solution[i+nbPatientLocal] ) ;

            }
        }
    }
}
/*
 printf("\n ******************Solution   solution1 avant \n");

                        for (int j = 0; j < nbPatientLocal *2; j++) {
                         printf(" %d ",  solution1[j]);
                    }

                    printf("********************\n");

         /**** supprimer les doublons de techniciens**********/

            for (int i = 0; i < newNbPatient; i++) {
 //printf(" %d  ",  fuzzyLocation[i]);
            if (fuzzyLocation [i]  == 0)
             {  //   printf(" je suis if  0  \n");
                 fuzzyLocation [i]  = 999;
            }
            else
            {
                if (fuzzyLocation [i] > nbTech)
                       {
                           fuzzyLocation [i]  = 999;
                          //   printf(" je suis if > nbTech \n");
                        }

              else {
                   //     printf(" je suis else  \n");
                for ( int t =0 ; t<nbTech; t++)
                { // c =0 ;
                 if (fuzzyLocation [i] == t +1 )
                 { if (occup_tec[t] <= moitier)
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
   /*          printf("\n  newNbPatient %d", newNbPatient);
      printf("\n hhhhhhhhhhhhhhhhhh************* ");
                       for (int j = 0; j < newNbPatient *2; j++) {
                        printf(" %d  ",  fuzzyLocation[j]);
                    }
         printf("\n  ");*/
         int y = newNbPatient ;
//         printf( "\n y /2 = %d \n ", y/2) ;
    //      printf("Length of the array fuzzyLocation is: %d\n", (int)( sizeof(fuzzyLocation) / sizeof(fuzzyLocation[0]) ));
      /**** Correction des doublons de techniciens**********/
        int infir = 0;
 /*       printf("\n charge tech  avans correction");
for (y = 0; y < nbTech ; y++) {
    printf( "\n %d  " ,   occup_tec[y]);
}*/
/*printf("\n charge patient  avans correction");
for (y = 0; y < newNbPatient; y++) {
    printf( "\n %d  " ,   occup[y]);
}*/
        for ( int  i=0 ; i< newNbPatient; i++)
        {
            if   (fuzzyLocation [i]  == 999)
            { // printf(" je suis uzzyLocation [i]    technixc  %d\n" , fuzzyLocation [i]);
                  do{
                 //       printf("%d ",  occup_tec[tech -1]);
                      infir = rand()% nbTech ;
                      // printf(" rand %d ", infir);
                    }  while (occup_tec[infir] >= moitier);
           //   printf("\n");
           //  printf("  ******************   %d   /////////////////////////// \n" , infir);
               occup_tec[infir] ++;
               fuzzyLocation [i]  =  infir +1;

/*
printf("\n charge tech ");
for (y = 0; y < nbTech ; y++) {
    printf( "\n %d  " ,   occup_tec[y]);
}*/
            }
        }
//printf("\n 11111");
         /**** supprimer les doublons des patients**********/

     for (int i = newNbPatient; i < newNbPatient*2 ; i++) {
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
            else if ( abs( fuzzyLocation [i] )  > nbPatientLocal)
                { fuzzyLocation [i]  = 999;
                }
          else
          {
                for ( int t =1 ; t<=nbPatientLocal; t++)
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

       /*     if (fuzzyLocation [i] == patient_annule)
                       {fuzzyLocation [i]  = 999;

                          //   printf(" je suis if > nbTech \n");
                          } */
     }//for
//   printf("\n 9999911111");
  //     printf("\n ");
    /*  printf("\n lllllllllllllllllllllllllllllllll************* ");
                       for (int j = 0; j < newNbPatient *2; j++) {
                        printf(" %d  ",  fuzzyLocation[j]);
                    }
         printf("\n  ");
/*
 printf("   \n  occup patient  \n");
                        for (int j = 0; j < nbPatientLocal; j++) {
                        printf(" %d   \n",  occup[j]);
                    }
 printf("\n ");
        /**** correction les doublons de patients**********/
//printf("y =%d \n", y);
        for (int i = y-1;  i < y*2 ; i++) {
      //         printf("  je suis iiiiiiiiiiiiii %d   je suis fuzzy loc %d  \n", i , fuzzyLocation [i]  );
            if   (fuzzyLocation [i]  == 999)
            {
                  do{
                        patient = rand()% (nbPatientLocal) ;
                       // printf("rand patient %d  ", patient);
                      }while (occup[patient ] ==1  );
                // printf("\n   lahna %d \n" );
               occup[patient] = 1;
               fuzzyLocation [i]  = patient+1 ;
              //  printf(" %d    %d \n",  occup[patient],    fuzzyLocation [i] );

             /*printf("\n Position PSO  correction  : ");
              for (int j = 0; j < y *2; j++) {
                 printf(" %d",  fuzzyLocation[j]);
            }*/
            }
             }
/*printf("  lahnaaaaaaaaaaaaaaaa njib b9iyet el vecteur\n" );
   printf("\n ******************Solution   fuzzyLocation\n");

                        for (int j = 0; j < newNbPatient*2; j++) {
                         printf(" %d ",  fuzzyLocation[j]);
                    }

                    printf("********************\n");
*/
 int q=0;
   for (int j =0; j< nbTech; j++)
        {
 for (int i = 0; i < newNbPatient; i++) {

            if (fuzzyLocation [i] ==  j+1)
            {
                for (int  g= 0; g< nbPatientLocal ; g++) {

                if (solution1[g] ==0)
                {
                       solution1[g] =  fuzzyLocation [i];
                       solution1[g+nbPatientLocal] =  fuzzyLocation [i+newNbPatient];
                        g = nbPatientLocal;
                }
            }
        }
    }
        }
 /*    printf("\n ******************Solution   solution1 aprés  \n");

                        for (int j = 0; j < nbPatientLocal *2; j++) {
                         printf(" %d ",  solution1[j]);
                    }
printf("********************\n");
*/
        for (int i = 0; i < nbPatientLocal*2 ; i++) {
                   newLocation[i] =  0;
             }
                 //   printf("********************\n");
             for (int i = 0; i < newNbPatient*2 ; i++) {
                   newLocation[i] =  fuzzyLocation [i];
             }
   // printf("fin mapping cor");

   /*   printf("\n ******************Solution au niveau mappingCorrection  \n");

                        for (int j = 0; j < nbPatientLocal *2; j++) {
                         printf(" %d ",  newLocation[j]);
                    }

                    printf("********************\n");
*/
}



/*********************************************************************************** Scheduling **/
void Scheduling( int solution[2*newNbPatient],      int* tab_Scheduling [nbPatientLocal+ 2*nbTech][5] ,int choix) {
        int y,k,g;
           int moitier =0;
        if (nbPatientLocal%2 !=0) moitier = (nbPatientLocal/2)+1 ;
        else moitier = nbPatientLocal/2;
        int tabTechnicien [moitier];
//       int tab_Scheduling [nbPatientLocal + 2*nbTech][3];
/*printf(" \n solution vecteur tab scheduling \n");
    for (int i =0; i<2*nbPatientLocal   ;i++)
    {printf("%d    \n", solution[i] ); }
*/
           //     printf("nbPatient local %d \n ", nbPatientLocal);
			for (y = 0; y < nbPatientLocal + 2*nbTech ; y++) {
			tab_Scheduling[y][0] = 0; //Technicien
			tab_Scheduling[y][1] = 0; //Patient
            tab_Scheduling[y][2] = 0; //
            tab_Scheduling[y][3] = 0; //
			tab_Scheduling[y][4] = 0; //
			}

 /*printf(" \n first tab scheduling \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<3 ; j++){
                printf(" %d",tab_Scheduling[i][j]);
            }printf("\n");
		 }printf("\n");
*/
int indice =0;
//printf("new nb patient = %d" , newNbPatient);
for (int t=1; t<=nbTech; t++)
{

  //  printf(" %d",t);
    for (int i =0; i<nbPatientLocal   ;i++)
    {

        if (solution[i] == t)
        { tab_Scheduling[indice] [0] =solution1[i] ;
            tab_Scheduling[indice] [1] =solution1[i+nbPatientLocal] ;
            //       printf(" \n solution[i+nbPatientLocal] %d  ", solution[i+nbPatientLocal] );
            indice ++;
        }
    }

tab_Scheduling[indice] [1] =0 ;
    indice ++;
}
/**********Calcule de temps de deplacement*************/
int  tech =1;
int  deb =0 ;
int fin =nbPatientLocal/nbTech +1;
int depart = 0, arrive =0, somme_dep=0;
		for (y = deb; y <= fin; y++) {
            //    printf("y  if  %d\n",  tech);
        arrive = tab_Scheduling[y][1];
        if (choix == 2)
        {
            tab_Scheduling[y][2]  =  new_deplacement [arrive] [depart]	;
        }
        else
        {
            tab_Scheduling[y][2]  =  deplacement [arrive] [depart]	; //temps de deplacement
        }

        depart =  arrive ;
      /*  somme_dep = somme_dep + tab_Scheduling[y][2] ;
         if (y == fin )
            {   tab_deplacement [0] [tech-1] = tech;
                 tab_deplacement [1] [tech-1] = somme_dep;
            }
            */
        if (y == fin && y <= (nbPatientLocal -(moitier)+2))
            {
                   deb = fin;
			     fin = fin+ (moitier)+2;
			     tech =tech+1;
			}
		}
  /*  printf(" \n third tab scheduling \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<5 ; j++){
                printf(" %d",tab_Scheduling[i][j]);
            }printf("\n");
	}printf("\n");

/	/**********Calcule de temps d'arrivé  + tmp depart*************/
int o =0,   oo=0;
		 for (int h = 0; h < nbTech; h++) {
	for (int y = 0; y < nbPatientLocal + 2*nbTech; y++) {
                   if  (  tab_Scheduling  [y][2] ==0){
                    tab_Scheduling  [y][3] =  0;
                     tab_Scheduling  [y][4] =  0;
                }
                 else{
                       o =   tab_Scheduling [y][2] ;
                       oo =  tab_Scheduling  [y-1][4] ;
                        tab_Scheduling  [y][3] =  o+ oo ;
                        tab_Scheduling  [y][4] =   o+ oo +20;
                 }
		}

}

// correction de tpm d'arrivé neud final
	 deb =0 ;
	 fin =(moitier -1)+2;
		for (y = deb; y <= fin; y++) {
			if (y == fin && y <= (nbPatientLocal -(moitier)+1))
            {   deb = fin;
			     fin = fin+ (moitier) +2;
			      		tab_Scheduling[y][4] = 0; //
			}
           if (y == nbPatientLocal + 2* nbTech -1)
            {
			      		tab_Scheduling[y][4] = 0; //
			}
		}
/**affichage**

    printf(" \n **********tab scheduling \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<5; j++){
                printf(" %d",tab_Scheduling[i][j]);
            }printf("\n");
		 }printf("\n");

 	/* printf("\n *****tab deplacement \n");
		 for (int i=0 ; i<nbTech; i++){
                printf(" %d  %d  \n",tab_deplacement[0][i], tab_deplacement[1][i]);
		 }printf("\n");

/***********************************/
}//end Scheduling

/*************** ***************************** UpdateLocationVelocity **/
void UpdateLocationVelocity (int t, Particule* p , int ind , int choix){

    if (choix == 1)
     newNbPatient  = nbPatientLocal - nbPatient_annul  +1;
 if (choix == 2)
    newNbPatient  = nbPatientLocal -  nbPatient_annul  +1  ;

   // printf("newNbPatient %d " , newNbPatient);
    double newVelocity [newNbPatient*2];
    int fuzzyLocation [newNbPatient*2];
    int tab_Scheduling [nbPatient + 2*nbTech][5];
    int i,j,y,x,bb;
 for (x = 0; x < newNbPatient*2; x++) {
//        printf(" \n %f     %f     ", p.PreviousLocation[x], p.PreviousVelocity[x]);
   fuzzyLocation[x] = 0;
  }

   for (x = 0; x < newNbPatient*2; x++) {
     r1 = (double)rand() / (double)RAND_MAX;
     r2 = (double)rand() / (double)RAND_MAX;
       //   double c= ;
    //   printf(" \n %f      %f    %f    %f    %f",  r1, r2, C1, C2, w );
   //   printf(" \n %f  ",  w );
         p->PreviousVelocity [x] = (w *  (double)  p->PreviousVelocity[x]
                + r1  * C1 *(double) ( p->bestLocation [x] -  p->PreviousLocation[x])) +r2 * C2*(double)(bestParticule.bestLocation[x]-  p->PreviousLocation[x]);
        //   printf(" \n %f     ",  p.PreviousVelocity[x] );
  }


 // update position -. fuzzyLocation [nbPatientLocal][2]
  for (x = 0; x < newNbPatient*2; x++) {
//        printf(" \n %f     %f     ", p.PreviousLocation[x], p.PreviousVelocity[x]);
   fuzzyLocation[x] = abs(  p->PreviousLocation[x] +  p->PreviousVelocity [x] );
  }
/*
if (choix == 2)
{
    printf("\n Position PSO  au niveau update  \n");

                for (int j = 0; j < newNbPatient *2; j++) {
                 printf(" %d",  fuzzyLocation[j]);
            }
            printf("\n");
}
*/
if (choix== 2)
{for (int h=0;h<nbPatientLocal;h++)
    {
         occup[h] = annul_p[h];
    }


        for ( int t =0 ; t<nbTech; t++)
                    {
                        occup_tec[t] = annul_t[t]  ;
                    }

//  printf("\n  mappingCorrection  \n");
    mappingCorrection2(fuzzyLocation); //return newLocation[nbPatientLocal][2]
}
else
{
    if (choix == 1)
    {
         mappingCorrection3(fuzzyLocation);
    }
    else {
         //  printf("\n  mappingCorrection  \n");
         mappingCorrection(fuzzyLocation); //return newLocation[nbPatientLocal][2]
    }

}

  for (x = 0; x < newNbPatient*2; x++) {
    p->PreviousLocation[x]= newLocation[x];
 //printf("%d  \t",(int)newLocation[x] );
  }
 // printf("\n fin mappingCorrection  \n");
/*if (choix == 2)
 {
  printf("  pp->bestLocation[j]   newLocationn[j]  p.PreviousVelocity[j]");
    for (int j = 0; j < nbPatientLocal *2; j++) {
                 printf(" \n %d             %d                       %f   ",   p->bestLocation[j], newLocation[j], p->PreviousVelocity[j]);
            }
printf("\n *************  \n ");
 }*/

     Scheduling2(newLocation, &tab_Scheduling,2) ;

 /* printf(" \n  tab_Scheduling up  \n");
		 for (int i=0 ; i<(newNbPatient + 2*nbTech); i++){
            for (int j=0 ; j<4 ; j++){
                printf(" %d", tab_Scheduling[i][j]);
            }printf("\n");
		 }printf("\n");*/
// printf("\n fin Scheduling  \n");
	for (int y = 0; y < newNbPatient + 2*nbTech; y++) {
            for(int k=0; k<5; k++)
		{
		      p->Scheduling  [y][k] = tab_Scheduling[y][k] ;
		}

}
 //printf("\n fin Scheduling  \n");
  int fitnessValue =0 ; // =  evaluation (&p,  ind );
                 int r = 0;
        for (x= 0; x < newNbPatient + 2*nbTech; x++) {
             r =   p->Scheduling [x] [3];
           if( fitnessValue < r)
            {
             fitnessValue = p->Scheduling [x] [3];
            }
             //  fitnessValue = MAX(fitnessValue, p->Scheduling [x] [3]);
        }
if (choix == 2)
{
    if(iteris == 0) {
   // printf("iteris 0 particule %d",  ind) ;

        p->bestFitness = fitnessValue ;
     if (fitnessValue <= ThirdbestParticule.bestFitness)
        {
            ThirdbestParticule.num_Particule=p->num_Particule;
            ThirdbestParticule.bestFitness = fitnessValue;
            ThirdbestParticule.bestIteration=0;

             for (int i = 0; i < newNbPatient+2*nbTech; i++){
                for (int j = 0; j < 5; j++){
                    ThirdbestParticule.Scheduling[i][j] = p->Scheduling[i][j] ;
                }//end for j
            }//end for i
        }
        //printf("\n fin fitnessValue %d  ******* p->bestFitness  %d ********* bestParticule.bestFitness %d\n" , fitnessValue , p->bestFitness,  bestParticule.bestFitness );
}
     //   printf("\n fin fitnessValue %d  ******* p->bestFitness  %d ********* bestParticule.bestFitness %d\n" , fitnessValue , p->bestFitness,  bestParticule.bestFitness );

          // printf("\n fin fitnessValue %d   bestParticule.bestFitness %d\n" , fitnessValue ,  bestParticule.bestFitness);

// printf( "newNbPatient  %d  rest %d"  , newNbPatient, rest_patient);
     if ( p->bestFitness >= fitnessValue){
        p->bestIteration = t;
        p->bestFitness = fitnessValue;
        for (i = 0; i < newNbPatient+2*nbTech; i++){
                p->bestLocation[i]=  p->PreviousLocation[i];
        }

     }
        if (fitnessValue <= ThirdbestParticule.bestFitness) {
            ThirdbestParticule.num_Particule=p->num_Particule;
            ThirdbestParticule.bestFitness = fitnessValue;
            ThirdbestParticule.bestIteration=t;
        for (i = 0; i < newNbPatient+2*nbTech; i++){
                ThirdbestParticule.bestLocation[i] =  p->bestLocation[i];
        }
              for (int i = 0; i < newNbPatient+2*nbTech; i++){
                for (int j = 0; j < 5; j++){
                    ThirdbestParticule.Scheduling[i][j] = p->Scheduling[i][j] ;
                }//end for j
            }//end for i
        }
/* printf(" \n  ThirdbestParticule.Scheduling[i][j]  \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<4 ; j++){
                printf(" %d", ThirdbestParticule.Scheduling[i][j]);
            }printf("\n");
		 }printf("\n");
        // end if
*/
}
else{
if(iteris == 0) {
   // printf("iteris 0 particule %d",  ind) ;

        p->bestFitness = fitnessValue ;
     if (fitnessValue <= bestParticule.bestFitness)
        {
            bestParticule.num_Particule=p->num_Particule;
            bestParticule.bestFitness = fitnessValue;
            bestParticule.bestIteration=0;
            }
        //printf("\n fin fitnessValue %d  ******* p->bestFitness  %d ********* bestParticule.bestFitness %d\n" , fitnessValue , p->bestFitness,  bestParticule.bestFitness );
}
 //printf( "newNbPatient  %d"  , newNbPatient);
     if ( p->bestFitness >= fitnessValue){
        p->bestIteration = t;
        p->bestFitness = fitnessValue;
        for (i = 0; i < newNbPatient+2*nbTech; i++){
                p->bestLocation[i]=  p->PreviousLocation[i];
        }

     }
        if (fitnessValue <= bestParticule.bestFitness) {
            bestParticule.num_Particule=p->num_Particule;
            bestParticule.bestFitness = fitnessValue;
            bestParticule.bestIteration=t;
        for (i = 0; i < newNbPatient+2*nbTech; i++){
                bestParticule.bestLocation[i] =  p->bestLocation[i];
        }
              for (int i = 0; i < newNbPatient+2*nbTech; i++){
                for (int j = 0; j < 5; j++){
                    bestParticule.Scheduling[i][j] = p->Scheduling[i][j] ;
                }//end for j
            }//end for i
        }
 /* printf(" \n  bestParticule.Scheduling[i][j]  \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<4 ; j++){
                printf(" %d", bestParticule.Scheduling[i][j]);
            }printf("\n");
		 }printf("\n");
        // end if
        */
}
} //end UpdateLocationVelocity
void  fn_new_request ()
{
      nbPatientLocal = nbPatient +1;
      printf("\n nbPatient = %d", nbPatient) ;
        printf("\n nbPatientLocal= %d", nbPatientLocal) ;
        printf("\n ");
    int y,k, h, S;
    //     int nbPatientLocal  = nbPatient- nbPatient_annul ;
// on suppose que le tmps d'annulation = 120 min
  //  patient_new=  7; //rand()% nbPatientLocal  +1;
  int tmp_request =150;
/*       printf("\n charge tech ");
for (int y = 0; y < nbTech ; y++) {
    printf( "\n %d  " ,   annul_t[y]);
}*/
int g ;
int solution[2*nbPatientLocal] ;

for (int t =0 ; t <=  nbPatientLocal* 2; t++)
{solution[t] =0 ;
     solution1[t] =0;
}
k =0;
   int val ;
   int i=0 ;
   int f=0 ;

for (h=0 ; h<nbTech ; h++)
{
    for (y = 0; y < nbPatient+ 2*nbTech; y++) {
         if ( bestParticule.Scheduling[y] [0] == h+1)
        {  // printf("\n   bestParticule.Scheduling[y] [3]   %d  "   , bestParticule.Scheduling[y] [3]);
         val = bestParticule.Scheduling[y] [3]  ;

             if(bestParticule.Scheduling[y] [3] > tmp_request )
                {
                   if ( bestParticule.Scheduling[y] [1] !=0)
                    {
                      //  printf(" bestParticule %d", bestParticule.Scheduling[y] [1]) ;
                        solution[i+nbPatientLocal] = bestParticule.Scheduling[y] [1] ;
                        solution [i] = h +1 ;
                       // printf("\n solution %d    %d ", solution [i]  , solution[i+nbPatientLocal] ) ;
                         i++;
                    }

            }
            else{
                for (g=0; g<5;g++)
                    { //patient_annule ++ ;
                        pt_depart[h][g] =  bestParticule.Scheduling [y][g]  ;
                         //      printf(" %d   ", pt_depart[k][g]);
                    }

        }

                        }
    }

}
secondbestParticule =bestParticule ;
           /**** Solution 1 preparing*********/
for (int h=0 ; h<nbTech ; h++)
{
    for (int y = 0; y < nbPatient+ 2*nbTech; y++) {
         if ( secondbestParticule.Scheduling[y] [0] == h+1 && secondbestParticule.Scheduling[y] [1] !=0)
        {  // printf("\n   bestParticule.Scheduling[y] [3]   %d  "   , bestParticule.Scheduling[y] [3]);
            if(secondbestParticule.Scheduling[y] [4] <=110 )
                {

                      //  printf(" bestParticule %d", bestParticule.Scheduling[y] [1]) ;
                        solution1[y+nbPatientLocal-1] = secondbestParticule.Scheduling[y] [1] ;
                        solution1[y-1] = h +1 ;
                       // printf("\n solution %d    %d ", solution [i]  , solution[i+nbPatientLocal] ) ;

            }
        }
    }
}
/**********************************/
 nbPatient_annul  =1 ;
annul_t[0] = 0;
annul_t[1] = 0;
  	for (y = 0; y < nbPatientLocal + 2*nbTech; y++) {
		    if (bestParticule.Scheduling[y][3] <=  tmp_request &&  bestParticule.Scheduling[y][3] !=0)// &&  bestParticule.Scheduling[y][1] != 0 )
		    {     nbPatient_annul ++;
         //  printf ( "\n bestParticule.Scheduling[y][3] =  %d ", bestParticule.Scheduling[y][3]);
		             annul_p [bestParticule.Scheduling[y][1] -1] =1  ;
		             annul_t  [bestParticule.Scheduling[y][0] -1 ] ++ ;
		    }
		}
/*		        printf("\n dans fct request Patient déja traité ");
    for (y = 0; y < nbPatient ; y++) {
    printf( "\n %d  " ,   annul_p[y]);
}
printf("\n charge tech ");
for (y = 0; y < nbTech ; y++) {
    printf( "\n %d  " ,   annul_t[y]);
}*/
printf("\n point de depard \n");

for (y = 0; y < nbTech; y++) {
                    for (g=0; g<5;g++)
                    {
                        printf(" %d   ", pt_depart[y][g]);
                    }printf("\n");
                }

printf("\n point à traiter \n");

for (int t =0 ; t < nbPatientLocal* 2; t++)
{
     printf("  %d   ", solution[t]);
     if (solution[t] !=0)   rest_patient ++ ;
}
printf("\n point déja traité \n");
for (int t =0 ; t <  nbPatientLocal* 2; t++)
{
     printf("  %d   ", solution1[t]);
}
int off ;
int tab_Scheduling [nbPatientLocal+ 2*nbTech][5] ;
printf("\n taper 1 si  on va continuer avec le mm planing; taper 2 si on va réaffecter les patients \n");
 // scanf("%d", &off);
 //if (off == 1)
 //{
   solution[i]= 1;
   solution[i+nbPatientLocal] = nbPatientLocal;
   printf("\n point à traiter \n");
for (int t =0 ; t < nbPatientLocal* 2; t++)
{
     printf("  %d   ", solution[t]);
}
   Scheduling3 (solution, tab_Scheduling, 2) ;

		printf("\n\nPlaning table : \n    Staff|  Pat |  Ttravel | Tarrived | Tleave \n");
		for (int y = 0; y < nbPatientLocal + 2*nbTech; y++) {
		    printf("\n");
			for (int h = 0; h < 5; h++) {
					printf(" %d         ",tab_Scheduling[y][h]);
            }
		}


 //}
 //else{
   /*     printf("\n point à traiter \n");
for (int t =0 ; t < nbPatientLocal* 2; t++)
{
     printf("  %d   ", solution[t]);
}*/
    iteris = 0;
    //Particule ParticuleRequest;
   // nbPatientLocal = nbPatientLocal - patient_annule;
  //  printf("deb init \n") ;
     initializeSwarm(&ThirdbestParticule , 2);
  //  printf("fin init \n") ;
    iterate(2);
 //}
}

/**************************** Fonction pour les threads des patients. ******************************/
 void  fn_annull ()
{   nbPatientLocal = nbPatient -1;
      printf("\n nbPatient = %d", nbPatient) ;
        printf("\n nbPatientLocal= %d", nbPatientLocal) ;
        printf("\n ");
    int y,k, h, S;
    //     int nbPatientLocal  = nbPatient- nbPatient_annul ;
// on suppose que le tmps d'annulation = 120 min
    patient_annule =  6; //rand()% nbPatientLocal  +1;
  int tmp_change =100;
  /*printf("\nPatient annulé \n");
  scanf("%d", &patient_annule);
  printf("tmp d'annulation ");
  scanf("%d", &tmp_change);
     annul_p[patient_annule-1] =1;
    occup[patient_annule-1] =1;
 printf("\n %d  \n", tmp_change);*/
int g ;
int solution[2*(nbPatientLocal)] ;
k =0;
   int val ;
   int i=0 ;
for (h=0 ; h<nbTech ; h++)
    {   	for (y = 0; y < nbPatientLocal+ 2*nbTech; y++) {
        if ( bestParticule.Scheduling[y] [0] == h+1)
        {   //printf("\n   bestParticule.Scheduling[y] [3]   %d  "   , bestParticule.Scheduling[y] [3]);
         val = bestParticule.Scheduling[y] [4]  ;
             if(val < tmp_change )
                {
                    for (g=0; g<5;g++)
                    {
                        pt_depart[h][g] =  bestParticule.Scheduling [y][g]  ;
                         //      printf(" %d   ", pt_depart[k][g]);
                    }
                }//printf("\n");
            else
            {

               solution [i] =  bestParticule.Scheduling[y] [1] ;
               i++;
            }
        }

                        }
    }
 //printf("\n %d  \n", tmp_change);
/*
for (y = 0; y < nbTech; y++) {
                    for (g=0; g<5;g++)
                    {
                                    printf(" %d   ", pt_depart[y][g]);
                    }printf("\n");
                }
*/
  nbPatient_annul  =1 ;

  	for (y = 0; y < nbPatientLocal + 2*nbTech; y++) {
		    if (bestParticule.Scheduling[y][4] <  tmp_change)// &&  bestParticule.Scheduling[y][1] != 0 )
		    {     nbPatient_annul ++;
		             annul_p [bestParticule.Scheduling[y][1] -1] =1  ;
		             annul_t  [bestParticule.Scheduling[y][0] -1 ] ++ ;
		    }
		}
/*for (y = 0; y < nbPatient ; y++) {
    printf( "\n %d  " ,   annul_p[y]);
}
printf("\n tech ");
for (y = 0; y < nbTech ; y++) {
    printf( "\n %d  " ,   annul_t[y]);
}*/
int off ;
int tab_Scheduling [nbPatientLocal - nbPatient_annul + 2*nbTech][5] ;
printf("taper 1 si  on va continuer avec le mm planing; taper 2 si on va réafficter les patients restants \n");
//  scanf("%d", &off);
 //if (off == 1)
 //{
     Scheduling (solution, tab_Scheduling , 1) ;

		printf("\n\nPlaning table : \n    Staff|  Pat |  Ttravel | Tarrived | Tleave \n");
		for (int y = 0; y < nbPatientLocal + 2*nbTech; y++) {
		    printf("\n");
			for (int h = 0; h < 5; h++) {
					printf(" %d         ",tab_Scheduling[y][h]);
            }
		}


 //}
 //else{
    iteris = 0;
    initializeSwarm(&bestParticule , 1);
    iterate(1);
 //}
}

/****************************************************** Scheduling 2 **/
void Scheduling2( int solution[2*nbPatientLocal],      int* tab_Scheduling [newNbPatient+ 2*nbTech][5] , int choix) {
        int y,k,g;
       int moitier =0;
        if (nbPatientLocal%2 !=0) moitier = (nbPatientLocal/2)+1 ;
        else moitier = nbPatientLocal/2;
        int tabTechnicien [moitier];
//       int tab_Scheduling [nbPatientLocal + 2*nbTech][3];
/*printf(" \n  vecteur tab scheduling solution \n");
    for (int i =0; i<2*nbPatientLocal   ;i++)
    {printf("%d    \n", solution[i] ); }
*/
            //    printf("nbPatient local %d \n ", nbPatientLocal);
			for (y = 0; y < nbPatientLocal + 2*nbTech ; y++) {
			tab_Scheduling[y][0] = 0; //Technicien
			tab_Scheduling[y][1] = 0; //Patient
            tab_Scheduling[y][2] = 0; //
            tab_Scheduling[y][3] = 0; //
			tab_Scheduling[y][4] = 0; //
			}

	int tech =1;
	int deb =0 ;
	int fin =(moitier -1)+2;
	int indimm =0;
	//  printf("dans ordonn  \n" );
         // on remplit  dans tab_Scheduling

int indice =0;
for (int t=1; t<=nbTech; t++)
{
    indice ++;
  //  printf(" %d",t);
    for (int i =0; i<newNbPatient ;i++)
    {

        if (solution[i] == t)
        {
            tab_Scheduling[indice] [0] =solution[i] ;
            tab_Scheduling[indice] [1] =solution[i+newNbPatient] ;
            //       printf(" \n solution[i+nbPatientLocal] %d  ", solution[i+nbPatientLocal] );
            indice ++;
        }
    }
  tab_Scheduling[indice] [0] =t ;
  tab_Scheduling[indice] [1] =0 ;
    indice ++;
}
/*
 printf(" \n first tab scheduling \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<5 ; j++){
                printf(" %d",tab_Scheduling[i][j]);
            }printf("\n");
		 }printf("\n");

/*******************correction point de départ***/
//if (choix !=2)
//{
      for (int i =0 ; i< nbTech ; i++)
    { y =0;
             if ( i==0 && y==0)
            {
          //   tab_Scheduling[y][0] =pt_depart[i][0] ; //Patient
           tab_Scheduling[y][0] =pt_depart[i][0] ; //Patient
            tab_Scheduling[y][1] =pt_depart[i][1] ; //Patient
            tab_Scheduling[y][2] = pt_depart[i][2] ; //temp
            tab_Scheduling[y][3] = pt_depart[i][3] ; //temp
             tab_Scheduling[y][4] = pt_depart[i][4] ; //temp
			}

    for (y = 1; y < nbPatientLocal+ 2*nbTech; y++) {
  			if ( tab_Scheduling[y][0]  == i+2)
            {
       //      tab_Scheduling[y][0] =pt_depart[i+1][0] ; //Patient
        tab_Scheduling[y-1][0] =pt_depart[i+1][0] ; //Patient
            tab_Scheduling[y-1][1] =pt_depart[i+1][1] ; //Patient
            tab_Scheduling[y-1][2] = pt_depart[i+1][2] ; //temp
            tab_Scheduling[y-1][3] = pt_depart[i+1][3] ; //temp
             tab_Scheduling[y-1][4] = pt_depart[i+1][4] ; //temp
              y = nbPatientLocal+ 2*nbTech;
			}
		}
    }
//}
/*
printf(" \n second tab scheduling \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<5 ; j++){
                printf(" %d",tab_Scheduling[i][j]);
            }printf("\n");
	}printf("\n");

/**********Calcul de temps de deplacement*************/
 tech =1;
 deb =0 ;
 fin =nbPatientLocal/nbTech +1;
 int depart = 0;
 if (choix == 2)
 {
     depart = pt_depart[1][0];
 }

int  arrive =0, somme_dep=0;
		for (y = deb; y <= fin; y++) {
            //    printf("y  if  %d\n",  tech);
        arrive = tab_Scheduling[y][1];
        if (choix == 2)
        {
            tab_Scheduling[y][2]  =  new_deplacement [arrive] [depart]	;
        }
        else{
		tab_Scheduling[y][2]  =  deplacement [arrive] [depart]	; //temps de deplacement
        }
        depart =  arrive ;

        if (y == fin && y <= (nbPatientLocal -(moitier)+2))
            {
                   deb = fin;
			     fin = fin+ (moitier)+2;
			     tech =tech+1;
			}
		}
/*
printf(" \n third tab scheduling \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<5 ; j++){
                printf(" %d",tab_Scheduling[i][j]);
            }printf("\n");
	}printf("\n");

	/**********Calcule de temps d'arrivé  + tmp depart*************/
int o =0,   oo=0;
		 for (int h = 0; h < nbTech; h++) {
	for (int y = 0; y < nbPatientLocal + 2*nbTech; y++) {
            if  (  tab_Scheduling  [y][0] ==h+1 &&  tab_Scheduling  [y+1][0] ==h+1  ){
                     o =   tab_Scheduling [y+1][2] ;
                       oo =  tab_Scheduling  [y][4] ;
                               tab_Scheduling  [y+1][3] =  o+ oo ;
                                tab_Scheduling  [y+1][4] =   o+ oo +20;
                }

		}

}
/******************
printf(" \n fourth tab scheduling \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<5 ; j++){
                printf(" %d",tab_Scheduling[i][j]);
            }printf("\n");
	}printf("\n");

// correction de tpm d'arrivé neud final
	 deb =0 ;
	 fin =(nbPatientLocal/nbTech -1)+2;
		for (y = deb; y <= fin; y++) {
			if (y == fin && y <= (nbPatientLocal -(nbPatientLocal/nbTech)+1))
            {   deb = fin;
			     fin = fin+ (nbPatientLocal/nbTech) +2;
			      		tab_Scheduling[y][4] = 0; //
			}
           if (y == nbPatientLocal + 2* nbTech -1)
            {
			      		tab_Scheduling[y][4] = 0; //
			}
		}
		/*****************/
/**affichage**

  printf(" \n **********tab scheduling \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<5; j++){
                printf(" %d",tab_Scheduling[i][j]);
            }printf("\n");
		 }printf("\n");

 	/* printf("\n *****tab deplacement \n");
		 for (int i=0 ; i<nbTech; i++){
                printf(" %d  %d  \n",tab_deplacement[0][i], tab_deplacement[1][i]);
		 }printf("\n");

/***********************************/
}//end Scheduling

/**************************************************************************** main **/
int main(){
//readBenchmark() ;
    int i ;
    for (int k=0 ; k< nbTech;k++)
 {   pt_depart[k][0] =  k+1;
      for (i=1; i<5;i++)
                    {
                        pt_depart[k][i] =  0;
                    }
 }
    initializeSwarm(&bestParticule , 0);
    iterate(0);
    //  fn_new_request () ;
    fn_annull();
   // getch();
}
// changer tt tab schuduling à 4 colonne dés l'initialisation
void Scheduling3( int solution[2*newNbPatient],      int* tab_Scheduling [newNbPatient+ 2*nbTech][5] , int choix) {
        int y,k,g;
        int tabTechnicien [newNbPatient/nbTech];
    //       printf("\n nbPatient = %d", nbPatient) ;
    //     printf("\n newNbPatient = %d", newNbPatient) ;
    //    printf("\n nbPatientLocal= %d", nbPatientLocal) ;
//       int tab_Scheduling [nbPatientLocal + 2*nbTech][3];
/*printf(" \n solution vecteur tab scheduling \n");
    for (int i =0; i<2*nbPatientLocal   ;i++)
    {printf("%d    \n", solution[i] ); }
*/
            //    printf("nbPatient local %d \n ", nbPatientLocal);
			for (y = 0; y < newNbPatient + 2*nbTech ; y++) {
			tab_Scheduling[y][0] = 0; //Technicien
			tab_Scheduling[y][1] = 0; //Patient
            tab_Scheduling[y][2] = 0; //
            tab_Scheduling[y][3] = 0; //
			tab_Scheduling[y][4] = 0; //
			}

	int tech =1;
	int deb =0 ;
	int fin =(newNbPatient/nbTech -1)+2;
	int indimm =0;
	//  printf("dans ordonn  \n" );
         // on remplit  dans tab_Scheduling
 /*printf(" \n first tab scheduling \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<3 ; j++){
                printf(" %d",tab_Scheduling[i][j]);
            }printf("\n");
		 }printf("\n");
*/
int indice =0;
for (int t=1; t<=nbTech; t++)
{
   indice ++;
  //  printf(" %d",t);
    for (int i =0; i<newNbPatient   ;i++)
    {

        if (solution[i] == t)
        {
            tab_Scheduling[indice] [0] =solution[i] ;
            tab_Scheduling[indice] [1] =solution[i+nbPatientLocal] ;
            //       printf(" \n solution[i+nbPatientLocal] %d  ", solution[i+nbPatientLocal] );
            indice ++;
        }
    }
  tab_Scheduling[indice] [0] =t ;
  tab_Scheduling[indice] [1] =0 ;
    indice ++;
}
/* printf(" \n first tab scheduling \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<5 ; j++){
                printf(" %d",tab_Scheduling[i][j]);
            }printf("\n");
		 }printf("\n");

/*******************correction point de départ***/

  for (int i =0 ; i< nbTech ; i++)
    { y =0;
             if ( i==0 && y==0)
            {
          //   tab_Scheduling[y][0] =pt_depart[i][0] ; //Patient
           tab_Scheduling[y][0] =pt_depart[i][0] ; //Patient
            tab_Scheduling[y][1] =pt_depart[i][1] ; //Patient
            tab_Scheduling[y][2] = pt_depart[i][2] ; //temp
            tab_Scheduling[y][3] = pt_depart[i][3] ; //temp
             tab_Scheduling[y][4] = pt_depart[i][4] ; //temp
			}

    for (y = 1; y < newNbPatient+ 2*nbTech; y++) {
  			if ( tab_Scheduling[y][0]  == i+2)
            {
       //      tab_Scheduling[y][0] =pt_depart[i+1][0] ; //Patient
        tab_Scheduling[y-1][0] =pt_depart[i+1][0] ; //Patient
            tab_Scheduling[y-1][1] =pt_depart[i+1][1] ; //Patient
            tab_Scheduling[y-1][2] = pt_depart[i+1][2] ; //temp
            tab_Scheduling[y-1][3] = pt_depart[i+1][3] ; //temp
             tab_Scheduling[y-1][4] = pt_depart[i+1][4] ; //temp
              y = nbPatientLocal+ 2*nbTech;
			}
		}
    }
/*printf(" \n second tab scheduling \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<5 ; j++){
                printf(" %d",tab_Scheduling[i][j]);
            }printf("\n");
	}printf("\n");

/**********Calcul de temps de deplacement*************/
 tech =1;
 deb =0 ;
 fin =nbPatientLocal/nbTech +1;
int depart = 0, arrive =0, somme_dep=0;
		for (y = deb; y <= fin; y++) {
            //    printf("y  if  %d\n",  tech);
        arrive = tab_Scheduling[y][1];
        if (choix == 2)
        {
            tab_Scheduling[y][2]  =  new_deplacement [arrive] [depart]	;
        }
        else{
		tab_Scheduling[y][2]  =  deplacement [arrive] [depart]	; //temps de deplacement
        }
        depart =  arrive ;

        if (y == fin && y <= (newNbPatient -(newNbPatient/nbTech)+2))
            {
                   deb = fin;
			     fin = fin+ (newNbPatient/nbTech)+2;
			     tech =tech+1;
			}
		}

/*printf(" \n third tab scheduling \n");
		 for (int i=0 ; i<(nbPatientLocal + 2*nbTech); i++){
            for (int j=0 ; j<5 ; j++){
                printf(" %d",tab_Scheduling[i][j]);
            }printf("\n");
	}printf("\n");

	/**********Calcule de temps d'arrivé  + tmp depart*************/
int o =0,   oo=0;
		 for (int h = 0; h < nbTech; h++) {
	for (int y = 0; y < newNbPatient + 2*nbTech; y++) {
            if  (  tab_Scheduling  [y][0] ==h+1 &&  tab_Scheduling  [y+1][0] ==h+1  ){

                     o =   tab_Scheduling [y+1][2] ;
                       oo =  tab_Scheduling  [y][4] ;
                               tab_Scheduling  [y+1][3] =  o+ oo ;
                                tab_Scheduling  [y+1][4] =   o+ oo +20;
                }

		}

}

}//end Scheduling
/***************************************************** mappingCorrection **/
void mappingCorrection3(int fuzzyLocation [2*(nbPatientLocal)])
{ //  printf("\n debut mapping correction \n") ;
      // printf("\n ");
            /*            for (int j = 0; j < nbPatient *2; j++) {
                        printf(" %d  ",  fuzzyLocation[j]);
                    }*/
        int moitier =0;
        if (nbPatientLocal%2 !=0) moitier = (nbPatientLocal/2)+1 ;
        else moitier = nbPatientLocal/2;
                   //  printf(" %d  ",  nbPatient_annul);
        int ratte = 0;
     //   int nbPatientLocal  = nbPatient - nbPatient_annul ;
        int patient =0, tech =0;
        int c =0 ;

         /**** supprimer les doublons de techniciens**********/

            for (int i = 0; i < nbPatientLocal ; i++) {
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
                 { if (  occup_tec[t] <= moitier)
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
 /*      printf("\n ");
                        for (int j = 0; j < nbPatientLocal *2; j++) {
                        printf(" %d  ",  fuzzyLocation[j]);
                    }
      /**** Correction des doublons de techniciens**********/
        int infir = 0;
        for ( int  i=0 ; i< nbPatientLocal-1; i++)
        {  if   (fuzzyLocation [i]  == 999)
            { // printf(" je suis uzzyLocation [i]    technixc  %d\n" , fuzzyLocation [i]);
                  do{
                 //       printf("%d ",  occup_tec[tech -1]);
                      infir = rand()% nbTech ;
                   // printf("%d  rand ", infir);
                    }  while (occup_tec[infir] >= nbPatientLocal/nbTech +1 );
             //printf("  %d  \n" , infir);
               occup_tec[infir] ++;
               fuzzyLocation [i]  =  infir +1;
            }
        }
//printf("11111");
         /**** supprimer les doublons des patients**********/

     for (int i = nbPatientLocal; i < nbPatientLocal*2 ; i++) {
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
            else if ( abs( fuzzyLocation [i] )  > nbPatientLocal)
                { fuzzyLocation [i]  = 999;
                }
          else
          {
                for ( int t =1 ; t<=nbPatientLocal; t++)
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

       /*     if (fuzzyLocation [i] == patient_annule)
                       {fuzzyLocation [i]  = 999;

                          //   printf(" je suis if > nbTech \n");
                          } */
     }//for
 //   printf("9999911111");
  /*      printf("\n ");
                        for (int j = 0; j < nbPatientLocal *2; j++) {
                        printf(" %d  ",  fuzzyLocation[j]);
                    }
 printf("   \n  occup teck \n");
                        for (int j = 0; j < nbPatient ; j++) {
                        printf(" %d   \n",  occup[j]);
                    }
 printf("\n ");
        /**** correction les doublons de patients**********/

        for (int i = nbPatientLocal;  i <=nbPatientLocal*2 ; i++) {
              //  printf("  je suis iiiiiiiiiiiiii %d   je suis fuzzy loc %d  \n", i , fuzzyLocation [i]  );
            if   (fuzzyLocation [i]  == 999)
            {
                  do{
                        patient = rand()% nbPatientLocal ;
                        //printf("%d  rand ", patient);
                      }while (occup[patient ] ==1  );

               occup[patient] = 1;
               fuzzyLocation [i]  = patient+1 ;
          //        printf(" %d    %d \n",  occup[patient],    fuzzyLocation [i] );
            }
             }
// printf("  lahnaaaaaaaaaaaaaaaaaaaaaaaaaaa \n" );

             for (int i = 0; i < nbPatientLocal*2 ; i++) {
                   newLocation[i] =  fuzzyLocation [i];
             }
   // printf("fin mapping cor");

   /*     printf("\n ******************Solution au niveau mappingCorrection  \n");

                        for (int j = 0; j < nbPatientLocal *2; j++) {
                         printf(" %d ",  newLocation[j]);
                    }

                    printf("********************\n");*/
}


