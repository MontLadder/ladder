#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <gmp.h>
#include <inttypes.h>
#include "intel_measurement_stuff.c"

int size;
gmp_randstate_t r;
mpz_t A, K, X, M, p, tempo, Xgmp, Xmont;
int i, nb_limbs, taille;
const mp_limb_t *p_limbs ;
mp_limb_t *A_limbs, *Z_limbs, *ZZA_limbs, *ZZB_limbs, *ZZC_limbs, *ZA_limbs, *C_limbs,*M_limbs, *X_limbs, 
*C1_limbs, *C2_limbs, *aux_limbs, *temp_limbs, *aux2_limbs, *tempo_limbs, *q_limbs;
mp_bitcnt_t bit;


void from_limbs_to_mpz_t(mpz_t rop, mp_limb_t *op, int op__nb_limbs){
	
	mpz_set_ui(rop, 0);
	
	if (op__nb_limbs == 0)
		return;
	
	int i, nb_dec;
	
	nb_dec = 8*sizeof(mp_limb_t);
    // nb_dec = 64
	
	mpz_set_ui(rop, op[op__nb_limbs - 1]);
	
	for(i=(op__nb_limbs - 2); i >=0; i--){
		mpz_mul_2exp (rop, rop, nb_dec);
		mpz_add_ui (rop, rop, op[i]);
	}
}

void mont_ladder()
{
    X_limbs[0] = 1;
    for (i = taille-1; i>=0; i--)
    {
        bit = mpz_tstbit(K,i);
        if (bit == 1)
        {
            mpn_mul_n(ZA_limbs,X_limbs,A_limbs,nb_limbs);
            mpn_tdiv_qr(q_limbs, X_limbs, 0, ZA_limbs, (nb_limbs*2), p_limbs, nb_limbs);
            mpn_sqr(ZA_limbs,A_limbs,nb_limbs);
            mpn_tdiv_qr(q_limbs, A_limbs, 0, ZA_limbs, (nb_limbs*2), p_limbs, nb_limbs); 
        }
        else
        {
            mpn_mul_n(ZA_limbs,X_limbs,A_limbs,nb_limbs);
            mpn_tdiv_qr(q_limbs, A_limbs, 0, ZA_limbs, (nb_limbs*2), p_limbs, nb_limbs);
            mpn_sqr(ZA_limbs,X_limbs,nb_limbs);
            mpn_tdiv_qr(q_limbs, X_limbs, 0, ZA_limbs, (nb_limbs*2), p_limbs, nb_limbs); 
        }
   }
}

void semi_ladder()
{
    X_limbs[0] = 1;
    // Calcul de C1=MA
  //  gmp_printf("limb M %Nx\n",M_limbs,nb_limbs);
  //  gmp_printf("limb A %Nx\n",A_limbs,nb_limbs);
    mpn_mul_n(ZA_limbs,M_limbs,A_limbs,nb_limbs);
    mpn_tdiv_qr(q_limbs, C1_limbs, 0, ZA_limbs, (nb_limbs*2), p_limbs, nb_limbs);
    // Calcul de C2=1-(C1.A+M) = p-(C1A+M)+1
    //C1A+M
    mpn_mul_n(ZA_limbs,C1_limbs,A_limbs,nb_limbs);
    mpn_add(aux2_limbs,ZA_limbs,2*nb_limbs,M_limbs,nb_limbs);
    mpn_tdiv_qr(q_limbs, C2_limbs, 0, aux2_limbs, (nb_limbs*2+1), p_limbs, nb_limbs);
    // p -(C1A+M)
    mpn_sub_n(C2_limbs,p_limbs,C2_limbs,nb_limbs);
    // p - (C1A+M)+1
    mpn_add_n(C2_limbs,C2_limbs,X_limbs,nb_limbs);
    //from_limbs_to_mpz_t(X,C1_limbs,nb_limbs);
    //gmp_printf("C1: %Zd\n",X);
    //from_limbs_to_mpz_t(X,C2_limbs,nb_limbs);
    //gmp_printf("C2: %Zd\n",X);
    for (i = taille-1; i>=0; i--)
    {
        bit = mpz_tstbit(K,i);
        //printf(" i = %d\n",i);
        if (bit == 1)
        {
            // z <- a^2
            mpn_sqr(ZA_limbs,A_limbs,nb_limbs);
            mpn_tdiv_qr(q_limbs, Z_limbs, 0, ZA_limbs, (nb_limbs*2), p_limbs, nb_limbs);
            /*from_limbs_to_mpz_t(X,Z_limbs,nb_limbs);
            gmp_printf("A^2: %Zd\n",X);
            gmp_printf("limb C1_limbs %Nx\n",C1_limbs,nb_limbs);*/

            // calcul de c_2xa 
            mpn_mul_n(ZA_limbs,X_limbs,A_limbs,nb_limbs);
            //mpn_tdiv_qr(q_limbs, aux_limbs, 0, ZA_limbs, (nb_limbs*2), p_limbs, nb_limbs);
            mpn_mul(ZZA_limbs,ZA_limbs,2*nb_limbs, C2_limbs,nb_limbs);
            //mpn_tdiv_qr(q_limbs, ZA_limbs, 0, aux2_limbs, (nb_limbs*2), p_limbs, nb_limbs);
            /*from_limbs_to_mpz_t(X,ZZA_limbs,3*nb_limbs);
            gmp_printf("C2XA: %Zd\n",X);           
            gmp_printf("limb C1_limbs %Nx\n",C1_limbs,nb_limbs);*/

            // calcul de c_1(x^2+z)
            mpn_sqr(aux_limbs,X_limbs,nb_limbs);
            /*from_limbs_to_mpz_t(X,aux_limbs,2*nb_limbs);
            gmp_printf("X^2: %Zd\n",X);           
            gmp_printf("limb C1_limbs %Nx\n",C1_limbs,nb_limbs);*/
            //mpn_tdiv_qr(q_limbs, temp_limbs, 0, aux_limbs, (nb_limbs*2), p_limbs, nb_limbs);
            //mpn_zero(aux2_limbs,2*nb_limbs+1);
            mpn_add(aux2_limbs,aux_limbs,2*nb_limbs,Z_limbs,nb_limbs);
            /*from_limbs_to_mpz_t(X,aux2_limbs,2*nb_limbs+1);
            gmp_printf("X^2+Z: %Zd\n",X);           
            gmp_printf("limb C1_limbs %Nx\n",C1_limbs,nb_limbs);*/
            //mpn_zero(ZZB_limbs,3*nb_limbs+1);    
            /*gmp_printf("limb aux2_limbs %Nx\n",aux2_limbs,2*nb_limbs+1);
            gmp_printf("limb C1_limbs %Nx\n",C1_limbs,nb_limbs);*/
            mpn_mul(ZZB_limbs,aux2_limbs,2*nb_limbs+1,C1_limbs,nb_limbs);
            /*from_limbs_to_mpz_t(X,ZZB_limbs,3*nb_limbs+1);
            gmp_printf("C1(X^2+Z): %Zd\n",X);           
            gmp_printf("limb C1_limbs %Nx\n",C1_limbs,nb_limbs);
            gmp_printf("limb aux2 %Nx\n",aux2_limbs,2*nb_limbs);*/
       
            // calcul de c_1(x^2+z)+c_2xa
            //mpn_zero(ZZC_limbs,3*nb_limbs+2);
            //gmp_printf("0 limb C1_limbs %Nx\n",C1_limbs,nb_limbs);
            mpn_add(ZZC_limbs,ZZB_limbs,3*nb_limbs+1,ZZA_limbs,(nb_limbs*3));
            //gmp_printf(" 1 limb C1_limbs %Nx\n",C1_limbs,nb_limbs); 
            //mpn_zero(X_limbs,nb_limbs);
            mpn_tdiv_qr(q_limbs, X_limbs, 0, ZZC_limbs, (2+nb_limbs*3), p_limbs, nb_limbs);
            //gmp_printf(" 2 limb C1_limbs %Nx\n",C1_limbs,nb_limbs);
            //from_limbs_to_mpz_t(X,X_limbs,nb_limbs);
           // gmp_printf("X: %Zd\n",X);
            //gmp_printf("limb C1_limbs %Nx\n",C1_limbs,nb_limbs);

            // a <- z
            mpn_copyi(A_limbs, Z_limbs, nb_limbs);
        }
        else
        {
            // z <- x^2
            mpn_sqr(ZA_limbs,X_limbs,nb_limbs);
            mpn_tdiv_qr(q_limbs, Z_limbs, 0, ZA_limbs, (nb_limbs*2), p_limbs, nb_limbs);
            // calcul de c_2xa 
            mpn_mul_n(ZA_limbs,X_limbs,A_limbs,nb_limbs);
            //mpn_tdiv_qr(q_limbs, aux_limbs, 0, ZA_limbs, (nb_limbs*2), p_limbs, nb_limbs);
            mpn_mul(ZZA_limbs,ZA_limbs,2*nb_limbs, C2_limbs,nb_limbs);
            //mpn_tdiv_qr(q_limbs, ZA_limbs, 0, aux2_limbs, (nb_limbs*2), p_limbs, nb_limbs);
            // calcul de c_1(a^2+z)
            mpn_sqr(aux_limbs,A_limbs,nb_limbs);
//            mpn_tdiv_qr(q_limbs, temp_limbs, 0, aux_limbs, (nb_limbs*2), p_limbs, nb_limbs);
            mpn_add(aux2_limbs,aux_limbs,2*nb_limbs,Z_limbs,nb_limbs);
            mpn_mul(ZZB_limbs,aux2_limbs,2*nb_limbs+1,C1_limbs,nb_limbs);
           // gmp_printf("limb aux2 %Nx\n",aux2_limbs,2*nb_limbs);
        // calcul de c_1(x^2+z)+c_2xa
            mpn_add(ZZC_limbs,ZZB_limbs,3*nb_limbs+1,ZZA_limbs,(nb_limbs*3));
            //mpn_zero(A_limbs,nb_limbs);
            mpn_tdiv_qr(q_limbs, A_limbs, 0, ZZC_limbs, (2+nb_limbs*3), p_limbs, nb_limbs);
            // x <- z
            mpn_copyi(X_limbs, Z_limbs, nb_limbs);
        }

   }
}

int main(int argc, char *argv[])
{

    unsigned long long timermin, timermax, t1, t2, diff_t;
    uint64_t cycles1[NTEST*NSAMPLES]={0}, cycles2[NTEST*NSAMPLES]={0}, cycles3[NTEST*NSAMPLES]={0};;
    unsigned long long meanTimer1_min=0, meanTimer2_min=0, meanTimer3_min=0;
    unsigned long long meanTimer1_max=0, meanTimer2_max=0, meanTimer3_max=0;
    unsigned long long *statTimer1, *statTimer2, *statTimer3;
    
    unsigned long long int START, STOP, START1,STOP1;
    uint64_t mini = (uint64_t)-1L, mini1 = (uint64_t)-1L;
    unsigned long long int timer=0, timer1=0;

    size = atoi(argv[1]);
    mpz_inits(A,K,X,M,p,tempo,Xgmp,Xmont,NULL);
	unsigned long seed = time(NULL);
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);


    srand(time(NULL));
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ////////////////////////// cache memory heating ////////////////////
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    mpz_urandomb(p,r,size);
    mpz_setbit(p,size-1);
    if (mpz_tstbit(p,0) == 0)
       mpz_setbit(p,0);
    //mpz_nextprime(p,p);
    nb_limbs = mpz_size(p);
    p_limbs = mpz_limbs_read(p); 
    X_limbs = (mp_limb_t *)calloc(nb_limbs, sizeof(mp_limb_t));
    Z_limbs = (mp_limb_t *)calloc(nb_limbs, sizeof(mp_limb_t));
    ZA_limbs = (mp_limb_t *)calloc(2*nb_limbs, sizeof(mp_limb_t));
    q_limbs = (mp_limb_t *)calloc(2*nb_limbs+4, sizeof(mp_limb_t));
    C_limbs = (mp_limb_t *)calloc(nb_limbs, sizeof(mp_limb_t));
    C1_limbs =  (mp_limb_t *)calloc(nb_limbs, sizeof(mp_limb_t));
    C2_limbs =  (mp_limb_t *)calloc(nb_limbs, sizeof(mp_limb_t));
    temp_limbs =  (mp_limb_t *)calloc(nb_limbs+1, sizeof(mp_limb_t));
    tempo_limbs =  (mp_limb_t *)calloc(nb_limbs, sizeof(mp_limb_t));
    aux_limbs =  (mp_limb_t *)calloc(2*nb_limbs, sizeof(mp_limb_t));
    aux2_limbs =  (mp_limb_t *)calloc(2*nb_limbs+1, sizeof(mp_limb_t));
    ZZA_limbs = (mp_limb_t *)calloc(3*nb_limbs, sizeof(mp_limb_t));
    ZZB_limbs = (mp_limb_t *)calloc(3*nb_limbs+1, sizeof(mp_limb_t));
    ZZC_limbs = (mp_limb_t *)calloc(3*nb_limbs+2, sizeof(mp_limb_t));
     
    int cptsemi = 0;
    int cptmont = 0;
    for(int i=0;i<NTEST;i++)
    {
        mpz_urandomm(A,r,p);
        mpz_urandomm(K,r,p);
        mpz_urandomm(M,r,p);
      /*  gmp_printf("A = %Zd\n",A);
        gmp_printf("K = %Zd\n",K);
        gmp_printf("M = %Zd\n",M); */
        taille = mpz_sizeinbase(K,2); 
        A_limbs = mpz_limbs_modify(A,nb_limbs);
        mpn_copyi(tempo_limbs,A_limbs,nb_limbs);
        //gmp_printf("A : %Nx\n",A_limbs,nb_limbs);
        M_limbs = mpz_limbs_modify(M,nb_limbs);
        mpz_powm_sec(Xgmp,A,K,p);
        mpn_zero(X_limbs,nb_limbs);
        semi_ladder();
        from_limbs_to_mpz_t(X,X_limbs,nb_limbs);
        mpz_limbs_finish(A,nb_limbs);
        mpn_copyi(A_limbs,tempo_limbs,nb_limbs);
        mpn_zero(X_limbs,nb_limbs);
        mont_ladder();
        from_limbs_to_mpz_t(Xmont,X_limbs,nb_limbs);
        if ((mpz_cmp(X,Xgmp) == 0)) 
            cptsemi += 1;
        if ((mpz_cmp(Xmont,Xgmp)==0))
            cptmont += 1;
    }
    printf("OK %d/%d %d/%d\n",cptsemi,NTEST,cptmont,NTEST);
    printf("\t  /*********************/\n");
    printf("\t / Timings !!!!!!!!!!!!/\n");
   printf("\t/*********************/\n\n");
   for(int i=0;i<NSAMPLES;i++){
	// Génération d'un jeu de paramètres aléatoires
    mpz_urandomm(A,r,p);
    mpz_urandomm(K,r,p);
    mpz_urandomm(M,r,p);
    taille = mpz_sizeinbase(K,2); 
    A_limbs = mpz_limbs_modify(A,nb_limbs);
    M_limbs = mpz_limbs_modify(M,nb_limbs);
    mpn_copyi(tempo_limbs,A_limbs,nb_limbs);
	timermin = (unsigned long long int)0x1<<63;
	timermax = 0;
	for(int j=0;j<NTEST;j++){
  		t1 = cpucyclesStart();
        mpz_powm_sec(X,A,K,p);
		t2 = cpucyclesStop();
		if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2+diff_t+1;
			}
		else
			diff_t = t2-t1;
			
		if(timermin > diff_t) 
			timermin = diff_t;
		else if(timermax < diff_t) 
			timermax = diff_t;
		cycles1[i*NTEST+j]=diff_t;
		}
	meanTimer1_min += timermin;
	meanTimer1_max += timermax; 
	timermin = (unsigned long long int)0x1<<63;
	timermax = 0;
	for(int j=0;j<NTEST;j++){
        mpn_copyi(A_limbs,tempo_limbs,nb_limbs);
        mpn_zero(X_limbs,nb_limbs);
  		t1 = cpucyclesStart();
        semi_ladder();
		t2 = cpucyclesStop();
		if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2+diff_t+1;
			}
		else
			diff_t = t2-t1;
			
		if(timermin > diff_t) 
			timermin = diff_t;
		else if(timermax < diff_t) 
			timermax = diff_t;
		cycles2[i*NTEST+j]=diff_t;
		}
	meanTimer2_min += timermin;
	meanTimer2_max += timermax; 
	timermin = (unsigned long long int)0x1<<63;
	timermax = 0;
	for(int j=0;j<NTEST;j++){
        mpn_copyi(A_limbs,tempo_limbs,nb_limbs);
        mpn_zero(X_limbs,nb_limbs);
  		t1 = cpucyclesStart();
        mont_ladder();
		t2 = cpucyclesStop();
		if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2+diff_t+1;
			}
		else
			diff_t = t2-t1;
			
		if(timermin > diff_t) 
			timermin = diff_t;
		else if(timermax < diff_t) 
			timermax = diff_t;
		cycles3[i*NTEST+j]=diff_t;
		}
	meanTimer3_min += timermin;
	meanTimer3_max += timermax;
  	} 
statTimer1 = quartiles(cycles1, NTEST*NSAMPLES);        
statTimer2 = quartiles(cycles2, NTEST*NSAMPLES);        
statTimer3 = quartiles(cycles3, NTEST*NSAMPLES);        
printf("               |    Q1   |    Q2   |    Q3   |\n");
printf(" ---------------------------------------------\n");
//printf("| gmp          | %7lld | %7lld | %7lld | \n", meanTimer1_min/NSAMPLES, meanTimer1_max/NSAMPLES, statTimer1[0], statTimer1[1], statTimer1[2]);
printf("| semi         | %7lld | %7lld | %7lld | \n",  statTimer2[0], statTimer2[1], statTimer2[2]);
printf("| mont         | %7lld | %7lld | %7lld | \n",  statTimer3[0], statTimer3[1], statTimer3[2]);
printf(" ----------------------------------------------\n");
printf("   semi/mont                %.3f\n",((float)statTimer2[1])/statTimer3[1]);
printf("\t  /**************************/\n");
printf("\t / Instructions !!!!!!!!!!!!/\n");
printf("\t/**************************/\n\n");
   for(int i=0;i<NSAMPLES;i++){
    // Génération d'un jeu de paramètres aléatoires
    mpz_urandomm(A,r,p);
    mpz_urandomm(K,r,p);
    mpz_urandomm(M,r,p);
    taille = mpz_sizeinbase(K,2); 
    A_limbs = mpz_limbs_modify(A,nb_limbs);
    M_limbs = mpz_limbs_modify(M,nb_limbs);
    mpn_copyi(tempo_limbs,A_limbs,nb_limbs);
    mini = (uint64_t)-1L, mini1 = (uint64_t)-1L;
    for(int j=0;j<NTEST;j++){
        mpn_copyi(A_limbs,tempo_limbs,nb_limbs);
        mpn_zero(X_limbs,nb_limbs);
        semi_ladder();
	}
    for(int j=0;j<NTEST;j++){
        mpn_copyi(A_limbs,tempo_limbs,nb_limbs);
        mpn_zero(X_limbs,nb_limbs);
  	START = rdpmc_instructions();
        semi_ladder();
	STOP = rdpmc_instructions();
	if(mini>STOP-START) mini = STOP-START;
	}
    timer += mini;
    for(int j=0;j<NTEST;j++){
        mpn_copyi(A_limbs,tempo_limbs,nb_limbs);
        mpn_zero(X_limbs,nb_limbs);
        mont_ladder();
        }
    for(int j=0;j<NTEST;j++){
        mpn_copyi(A_limbs,tempo_limbs,nb_limbs);
        mpn_zero(X_limbs,nb_limbs);
  	START1 = rdpmc_instructions();
        mont_ladder();
	STOP1 = rdpmc_instructions();
	if(mini1>STOP1-START1) mini1 = STOP1-START1;
        }
    timer1 += mini1;    
}
printf("IC semi ladder        : %llu\n",timer/NSAMPLES);	
printf("IC Montgomery ladder  : %llu\n",timer1/NSAMPLES);
}
