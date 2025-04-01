/* -------------------------------------------------------------

   file: solvde4.c

   Richard E. Zeebe    (zeebe@soest.hawaii.edu)
   Dieter Wolf-Gladrow (Dieter.Wolf-Gladrow@awi.de)

   updates:


     11/17/08 # of output files reduced
        12/02 boron isotope runs
              (Zeebe et al., Paleoceanography 2003)
     06.06.01 solvde4.c
              pH scale -> total scale
              revise all reactions and rate constants
              (cf. Zeebe and Wolf-Gladrow, 2001)
     09.09.96 k-1/k'-1 = 1.022 k-4/k'-4 = 1.011
     	      (O'Leary, pers. comm.)
         6/96 Erratum (Mar. Chem.) K_1 slightly changed (0.0846834)
              and dicosw (0.0007379) (dwg).
     07.05.96 km4 = kp4 * Kwater * km1s / kp1s
	          small error in Kwater corrected (REZ)
     04.05.96 convert all other equilibrium coeff. to 'free' scale
     30.04.96 introduce loop for runs for different ph values
     25.03.96 start including 13C Isotopes.
         1/96 borates included (DWG)
     13.12.95 fdco3: D_CO3-- = 0.955e-9 at 25 degrees.
     12.12.95 error in fvisc corrected
     04.07.95 k_CO2 (T,S)
     31.03.95 modify K_1, K_2 from pH_SWS to pH_NBS
              My version, including Dieter's updates from
              31.03.95, 4.07.95, and 04.05.96
     23.03.95 Richard Zeebe
              (email: rzeebe@awi-bremerhaven.de)
     10.08.94 small error in OH corrected
     13.04.94 include oxygen (forams)
     19.07.93 debugged !!!
     15.06.93 solvde3.c, Dieter Wolf-Gladrow
              (email: wolf@awi-bremerhaven.de)

   purpose: solve a system of ordinary differential equations:
            diffusion-reaction system

   remarks:
     1. units: SI except length in [mu = 1.e-6 m] and
               concentrations in [mumol/kg]
     2. uptake = 4 * PI * RADIUS**2 * D * dc/dr

     3. IMAX-arrays -> M-arrays; r1,... remove?

     4.       order of equations:
              1. CO2    2. HCO3-   3. CO3--   4. H+     5. OH-
	     (can vary: 6. 13CO2	7. H13CO3- 8. 13CO3-- 9. BOH3  10. BOH4
	      11. O2)
   C:

        K&R or ANSI: ANSI

     1. C is case-sensitive
     2. C: call by value; FORTRAN: call by reference
     3. convention: "define-variables" uppercase

Notation:

     paper                C-code

\def\kthon{k'_{+1}}     % kp1s
\def\konth{k'_{-1}}     % km1s
\def\konfo{k_{+4}}      % kp4
\def\kfoon{k_{-4}}      % km4
\def\konfi{k_{+5}}      % kp5h
\def\kfion{k_{-5}}      % km5h
			% kp5oh
			% km5oh
\def\ksion{k_{+6}}      % kp6
\def\konsi{k_{-6}}      % km6
\def\kbf{k_{-7}}        % km7 bor forward
\def\kbb{k_{+7}}        % kp7 bor backward


--------------------------------------------------- */

#include <stdio.h>
#include <math.h>
#include "nrutil.c"

#define SQ(x) ((x)*(x))
#define KU(x) ((x)*(x)*(x))

/* -----   debug options   ----- */

#define UDEB05
#define UDEB07

/* -------- run options --------------- */

#define INITA	       /* init arrays        */
#define UPTAKE
#define REACTION
#define UANASOL    1   /* analytical solution */

/*---------  TOTAL or FREE ph scale -------  */

#define UFREEPH

/*---------  OLD RATE CONSTANTS (SOLVDE3)-------  */

#define URC3

/* -------- switch reactions on or off ----- */

/* <----- comment !!
#define NOCO2
#define NOHCO3
#define NOCO3
#define NOHPLUS
#define NOOH		*/


#define UCARTEST
/*
#define NOK5
*/

/* ------------------------------------------------------------------------ */


#define PI 3.1415926535897932384626433832795	   /*             Gaudi Pi  */
#define MOLH2O 55.56e6   /*   	         H2O is 55 molar = 55e6 [mumol/kg]  */
#define EACTK 1.5e4      /* = 15 kcal/mol activation energy for rate coeff. */


/* -----  Forams, Aggregates, Coccos, Diatoms                         ----- */
#define PYTHON    /* FORAMW,FORAMOUD,FORAMOUL,FORAMGSD,FORAMGSL
			           FORAMRES,FORAMPHS,FORAMCLC
			           FORAMBORD FORAMBORL
                 PYTHON  -  parameters set by python                        */
#define UAGG
#define UCOCCOW
#define UDIATOMW



#define UPRINT

#define UFLSYMUPT	/* write data of C uptake in file 	*/

/*-------    BORON, OXYGEN AND CALCIUM: on/off		--------*/

/* OXYGEN  can only be included if BORON is on!			*/
/* CALCIUM can only be included if BORON and OXYGEN is on!	*/


#define BORON		/* BORTBULK, after SALINITY	*/
#define OXYGEN		/* O2BULK,   see FORAMS		*/
#define CALCIUM


#ifdef BORON
#define BORONRC4
#endif

#define F13_CO3		   /*  CO3UPT for calc.*/
#define UF13_HCO3      /* HCO3UPT for calc.*/

#define UBORISTP
#define B10B11		/* 10B and 11B are variables, not 11B and total B */
#ifdef BORISTP
#define BSTAND 4.0014	/* (Hemming & Hanson,92), SRM 951 NBS */
#define D11BSEA 39.5    /* delta 11B value of sea water */
#define EPSB ((1.0194-1.)*1.e3)	/* Kakihana et al., 1977	*/
#endif

#define UC13ISTP		/* include 13C Isotopes	        */
#define UCISTP		/* off: total C and 13C		*/
			/* on :     12C and 13C		*/

#ifdef C13ISTP
#ifdef SYMTCIN
 #define DC13DIC 1.5	/* delta C13 of total CO2		*/
#else
 #define DC13DIC 2.0	/* delta C13 of total CO2		*/
#endif
#define RSTAND  0.01124 /* (13C/12C)_standard, O'Leary (1980) 	*/
#define D13CRES (-21.9) /* delta value of respired CO2		*/
#define D13CPHS (-22.)  /* delta value of photosynth. prod.	*/
			/* (NO MIMECO2SYM !!)			*/
#define UJASPER
#define CEPSP
#ifdef JASPER
 #define CO2EPSP 10.0    /* umol/kg  CO2 concentration for partitioning
			   between Jaspers formula for CO2 fractionation
			   (epsp = A - B/co2) and linear decrease
			   until co2 = 0.			*/
 #define AEPSP	27.
 #define BEPSP	130.
#endif
#ifdef CEPSP
 #define EPSPCO2   18.	/* fractionation associated with CO2  upt.*/
 #define EPSPHCO3  18.	/* fractionation associated with HCO3 upt.*/
#endif

#endif

#define UCBNS		/* pH and dic or alk as arguments of main*/

#define UCBNSLOOP	/* run the program for several ph values*/
                    /* remark: CBNS + CBNSLOOP not possible!*/

#ifdef  CBNSLOOP

#define READ 		/* read pH and DIC from file		*/
#define ULIDA		/* DARK for 12 hours light/dark		*/
#ifdef  READ
#define LDAT
 #ifdef DDAT
  #define NDAT    13	/* DARK: 13 (cont. dark) x (dark/light)*/
 #endif
 #ifdef LDAT
  #define NDAT    14	  /* LIGHT: 14 				    */
  #define SYMTCUPPH89 2.8 /* C uptake of symbionts at pH 8.9 nmol/h */
 #endif
#endif

#ifndef READ		/* increase pH linear	 		*/
#define PHMIN 	7.7
#define PHMAX	8.6
#define ICLMAX	4	/* Number of pH values			*/
#endif
#endif


/* Added 25/01/2017 by Branson and Holland! */
#ifdef PYTHON                /* Parameters defined by Python */
#define RADIUS 2.500000e+02      /* radius of foram [mu] diatoms */
  /*  UPTAKE at the shell :           */
#define CO3UPT   8.333333e-13  /* .75/3.25 direct calcification */
#define CO2UPT  -5.555556e-13 /* 3+2[mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT 0.000000e+00
#define O2UPT (-CO2UPT)
#define CAUPT   CO3UPT
#define OHUPT 0.0
#define HUPT  (0.0e-20)
#define BOH4UPT 0.0   /* (1.0e-9/3600.*6.8e-5) */
#define REDF  1.0   /* Redfield ratio O2=R*CO2 foram resp.*/
        /* see: O2UPT at the shell  */
#define PHBULK  8.063000e+00    /* 8.2 /8.16      */
#define DICBULK 4.035000e+03       /* 2167.      */
#define UALKBULK 4.671000e+03    /* 2723.      */

/*  UPTAKE of symbionts : 16.2 = 18 (gross) - 1.8 (symb. resp.)*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT 0.000000e+00  /* CO2  uptake by symbionts */
#define SYMHCO3UPT  0.000000e+00  /* HCO3 uptake by symbionts */
#define SYMHUPT SYMHCO3UPT    /* H    uptake by symbionts */

#define MIMECO2SYM/*----   2. MIME: set total carbon  uptake     ---*/
#define SYMTCUPT 2.777778e-12
#define VMAX 3.333333e-12    /* Michaelis-Menten max upt. for CO2 */
#define KS 5.         /* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT 3.0   /* increase of uptake per iteration */
#define SYMRAD (RADIUS+5.000000e+02)   /* 500/200 outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.000000e+00    /* Redfield ratio symbiont photosynth.  */
#define O2BULK 210.  /* [mumol/kg] Joergensen, 1985 */

#endif


#ifdef FORAMW
#define RADIUS 200.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/
	/*	UPTAKE at the shell : 					*/
#define CO3UPT   (0.0e-9/3600.)  /* .75/3.25 direct calcification	*/
#define CO2UPT  (-3.0e-9/3600.)	/* -3 [mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT (0.0e-9/3600.)  /* 2 */
#define O2UPT (-CO2UPT)
#define CAUPT   (0.0e-9/3600.)
#define OHUPT 0.0
#define HUPT  (0.0e-20)
#define BOH4UPT	0.0		/* (1.0e-9/3600.*6.8e-5) */
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK  8.2		/* 8.25				*/
#define DICBULK 2200. 		/* 2167.			*/
#define UALKBULK 2400.	 	/* 2723.			*/

/*	UPTAKE of symbionts : 16.2 = 18 (gross) - 1.8 (symb. resp.)*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.)  /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake     ---*/
#define SYMTCUPT (12.e-9/3600.) /* 12 */
#define VMAX (12.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 5.      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Joergensen, 1985 */

#endif


#ifdef FORAMBORD		/* HCO3- uptake */
#define RADIUS 250.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/
	/*	UPTAKE at the shell : 					*/
#define CO3UPT   (1.0e-9/3600.)  /* .75/3.25 direct calcification	*/
#define CO2UPT  (-3.0e-9/3600.)	/* 1+2 [mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT (2.0e-9/3600.)
#define O2UPT (-CO2UPT)
#define CAUPT   CO3UPT
#define OHUPT 0.0
#define HUPT  (0.0e-20)
#define BOH4UPT	0.0		/* (1.0e-9/3600.*6.8e-5) */
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK  7.623	 	 /* 7.9,8.16,8.5 / 8.2	     	 */
#define DICBULK 2054. 		/* 2167.			*/
#define UALKBULK 2214.	 	/* 2723.			*/

/*	UPTAKE of symbionts : 16.2 = 18 (gross) - 1.8 (symb. resp.)*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.)  /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake     ---*/
#define SYMTCUPT (10.e-9/3600.)
#define VMAX (12.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 5.      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Joergensen, 1985 */

#endif

#ifdef FORAMBORD2		/* CO32- uptake */
#define RADIUS 250.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/
	/*	UPTAKE at the shell : 					*/
#define CO3UPT   (1.0e-9/3600.)  /* .75/3.25 direct calcification	*/
#define CO2UPT  (-3.0e-9/3600.)	/* 1+2 [mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT (0.0e-9/3600.)
#define O2UPT (-CO2UPT)
#define CAUPT   CO3UPT
#define OHUPT 0.0
#define HUPT  (0.0e-20)
#define BOH4UPT	0.0		/* (1.0e-9/3600.*6.8e-5) */
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK  8.063	  /* 7.9,8.16,8.5 / 8.2	  	  	  */
#define DICBULK 4035. 	    /* 2167.			*/
#define UALKBULK 4671.	 	/* 2723.			*/

/*	UPTAKE of symbionts : 16.2 = 18 (gross) - 1.8 (symb. resp.)*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.)  /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake     ---*/
#define SYMTCUPT (10.e-9/3600.)
#define VMAX (12.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 5.      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Joergensen, 1985 */

#endif

#ifdef FORAMBORD3		/* SIZE */
#define RADIUS 125.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/
	/*	UPTAKE at the shell : 					*/
#define FSZ (RADIUS/250.)     /* Size factor */
#define CO3UPT   (0.0e-9/3600.)*FSZ  /* .75/3.25 direct calcification	*/
#define CO2UPT  (-3.0e-9/3600.)*FSZ	/* 1+2 [mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT (2.0e-9/3600.)*FSZ
#define O2UPT (-CO2UPT)
#define CAUPT   CO3UPT
#define OHUPT 0.0
#define HUPT  (0.0e-20)
#define BOH4UPT	0.0		/* (1.0e-9/3600.*6.8e-5) */
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK  8.16		/* 7.9,8.16,8.5 / 8.2			*/
#define DICBULK 2000. 		/* 2167.			*/
#define UALKBULK 2400.	 	/* 2723.			*/

/*	UPTAKE of symbionts : 16.2 = 18 (gross) - 1.8 (symb. resp.)*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.)  /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake     ---*/
#define SYMTCUPT (12.e-9/3600.)
#define VMAX (12.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 5.      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Joergensen, 1985 */

#endif



#ifdef FORAMBORL1		/* HCO3- uptake */
#define RADIUS 250.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/
	/*	UPTAKE at the shell : 					*/
#define CO3UPT   (0.0e-9/3600.)  /* .75/3.25 direct calcification	*/
#define CO2UPT  (-5.0e-9/3600.)	   /*3+2[mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT (6.0e-9/3600.)
#define O2UPT (-CO2UPT)
#define CAUPT   CO3UPT
#define OHUPT 0.0
#define HUPT  (0.0e-20)
#define BOH4UPT	0.0		/* (1.0e-9/3600.*6.8e-5) */
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK  7.623	  /* 8.2 /8.16	  	  	  */
#define DICBULK 2054. 		/* 2167.			*/
#define UALKBULK 2214.	 	/* 2723.			*/

/*	UPTAKE of symbionts : 16.2 = 18 (gross) - 1.8 (symb. resp.)*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.)  /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define MIMECO2SYM/*----   2. MIME: set total carbon  uptake     ---*/
#define SYMTCUPT (10.e-9/3600.)
#define VMAX (12.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 5.      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* 500/200 outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Joergensen, 1985 */

#endif


#ifdef FORAMBORL2                /* CO32- uptake */
#define RADIUS 250.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/
	/*	UPTAKE at the shell : 					*/
#define CO3UPT   (3.0e-9/3600.)  /* .75/3.25 direct calcification	*/
#define CO2UPT  (-2.0e-9/3600.)	/* 3+2[mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT (0.0e-9/3600.)
#define O2UPT (-CO2UPT)
#define CAUPT   CO3UPT
#define OHUPT 0.0
#define HUPT  (0.0e-20)
#define BOH4UPT	0.0		/* (1.0e-9/3600.*6.8e-5) */
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK  8.063	 	 /* 8.2 /8.16	 	 	 */
#define DICBULK 4035. 	    /* 2167.			*/
#define UALKBULK 4671.	 	/* 2723.			*/

/*	UPTAKE of symbionts : 16.2 = 18 (gross) - 1.8 (symb. resp.)*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.)  /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define MIMECO2SYM/*----   2. MIME: set total carbon  uptake     ---*/
#define SYMTCUPT (10.0e-9/3600.)
#define VMAX (12.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 5.      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* 500/200 outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Joergensen, 1985 */

#endif

#ifdef FORAMBORL		/* HCO3- uptake */
#define RADIUS 125.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/
	/*	UPTAKE at the shell :            		*/
#define FSZ KU(RADIUS/250.)     /* Size factor */
#define CO3UPT   (0.0e-9/3600.)*FSZ  /* .75/3.25 direct calcification	*/
#define CO2UPT  (-5.0e-9/3600.)*FSZ	/* 3+2[mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT (6.0e-9/3600.)*FSZ
#define O2UPT (-CO2UPT)
#define CAUPT   CO3UPT
#define OHUPT 0.0
#define HUPT  (0.0e-20)
#define BOH4UPT	0.0		/* (1.0e-9/3600.*6.8e-5) */
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK  8.16		/* 8.2 /8.16			*/
#define DICBULK 2000. 		/* 2167.			*/
#define UALKBULK 2400.	 	/* 2723.			*/

/*	UPTAKE of symbionts : 16.2 = 18 (gross) - 1.8 (symb. resp.)*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.)  /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define MIMECO2SYM/*----   2. MIME: set total carbon  uptake     ---*/
#define SYMTCUPT (10.e-9/3600.)*FSZ
#define VMAX (12.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 5.      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	1.0      	/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+200.0*RADIUS/250.)   /* 500/200 outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Joergensen, 1985 */

#endif


#ifdef FORAMOUD		/* Orbulina Universa, Dark */
	/*					*/
        /* realistic values :			*/
	/* 2. respiration of foram		*/
	/*   -2.0 nmol/h 			*/
	/* 3. calcification			*/
	/*    1. nmol/h			*/

#define RADIUS 300.0   		/* radius of foram [mu] diatoms     */
				/* Bulk radius, see BULKFAC	    */
				/* Num. of mesh points, see #define M */

/*	UPTAKE at the shell : 	(0.0e-9/3600.)			*/
#define CO3UPT  (1.0e-9/3600.)  /*  1. direct calcification	*/
#define CO2UPT  (-2.1e-9/3600.)	/* -2.1[mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT 0.0
#define CAUPT (1.e-9/3600.)     /* (100.7176e-9/3600.) */
#define OHUPT 0.0
#define HUPT  (-HCO3UPT)	/* 0.0 				*/
#define BOH4UPT	(1.0e-9/3600.*6.8e-5)
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 	8.2	   	   /* >= 7.9 !!! 7.7, 8.15, 8.55	*/
#define UDICBULK 2015. 		/* 2010.			*/
#define UALKBULK 2400.	 	/* 2723.			*/

/*	UPTAKE of symbionts : 					*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(-0.0e-9/3600.)  /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake   -----*/
#define SYMTCUPT (0.0e-9/3600.)
#define VMAX (10.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 1.5      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 200.  /* [mumol/kg] Bijma et al., 1992, p. 252 */

#endif

#ifdef FORAMOUD2	/* Orbulina Universa, Dark (Steffi) */
	/*					*/
        /* realistic values :			*/
	/* 2. respiration of foram		*/
	/*   -2.0 nmol/h 			*/
	/* 3. calcification			*/
	/*    1. nmol/h			*/

#define RADIUS 242.0   		/* radius of foram [mu] diatoms     */
				/* Bulk radius, see BULKFAC	    */
				/* Num. of mesh points, see #define M */

/*	UPTAKE at the shell : 	(0.0e-9/3600.)			*/
#define CO3UPT  (1.0e-9/3600.)  /*  direct calcification	*/
#define CO2UPT  (-1.96e-9/3600.)/* [mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT 0.0
#define CAUPT CO3UPT
#define OHUPT 0.0
#define HUPT  (-HCO3UPT)	/* 0.0 				*/
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 	8.2		/* >= 7.9 !!!			*/
#define UDICBULK 2015. 		/* 2010.			*/
#define UALKBULK 2400.	 	/* 2723.			*/

/*	UPTAKE of symbionts : 					*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(-1.39e-9/3600.)  /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake   -----*/
#define SYMTCUPT (0.0e-9/3600.)
#define VMAX (10.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 1.5      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 200.  /* [mumol/kg] Bijma et al., 1992, p. 252 */

#endif

#ifdef FORAMOUD3	/* Orbulina Universa, Dark */
	/*					*/
        /* realistic values :			*/
	/* 2. respiration of foram		*/
	/*   -2.0 nmol/h 			*/
	/* 3. calcification			*/
	/*    1. nmol/h				*/
#define D13CRES (-14.537) 	/* delta value of respired CO2      */
#define RADIUS 300.0   		/* radius of foram [mu] diatoms     */
				/* Bulk radius, see BULKFAC	    */
				/* Num. of mesh points, see #define M */

/*	UPTAKE at the shell : 	(0.0e-9/3600.)			*/
#define CO3UPT  (0.0e-9/3600.)  /*  1. direct calcification	*/
#define CO2UPT  (-3.1e-9/3600.)	/* -2.1[mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT (2.*1.0e-9/3600.)
#define CAUPT (100.7176e-9/3600.)
#define OHUPT 0.0
#define HUPT  0.0	/* (-HCO3UPT) 				*/
#define BOH4UPT	(1.0e-9/3600.*6.8e-5)
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 	8.2	   	   /* >= 7.9 !!!	   	   */
#define UDICBULK 2400. 		/* 2010.			*/
#define UALKBULK 2819.	 	/* 2723.			*/

/*	UPTAKE of symbionts : 					*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(-0.0e-9/3600.)  /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake   -----*/
#define SYMTCUPT (0.0e-9/3600.)
#define VMAX (10.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 1.5      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 200.  /* [mumol/kg] Bijma et al., 1992, p. 252 */

#endif


#ifdef FORAMOUL	/* Orbulina Universa, Light*/
	/*					*/
        /* realistic values :			*/
	/* 1. photos. of symbionts		*/
	/*    10. nmol/h 			*/
	/* 2. respiration of foram		*/
	/*   -2.5 nmol/h 			*/
	/* 3. calcification			*/
	/*    2.2 nmol/h			*/

#define RADIUS 300.0   		/* radius of foram [mu] diatoms     */
				/* Bulk radius, see BULKFAC	    */
				/* Num. of mesh points, see #define M */

/*	UPTAKE at the shell : 	(0.0e-9/3600.)			*/
#define CO2UPT  (-2.1e-9/3600.)	/* [mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT (0.0e-9/3600.)
#define OHUPT 0.0
#define HUPT  0.0		/* 0.0 				*/
#define BOH4UPT	0.0
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 	8.2	   	   /* >= 7.	   	   	   */
#define UDICBULK 2015. 		/* 2010.			*/
#define UALKBULK 2400.	 	/* 2400.			*/

/*	UPTAKE of symbionts : 					*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#ifdef  SYMCO2IN
 #define SYMCO2UPT	(SYMCO2IN*1.e-9/3600.)
 #define SYMHCO3UPT	(SYMHCO3IN*1.e-9/3600.)
 #define CO3UPT  	(CO3IN*1.e-9/3600.)
#else
 #define SYMCO2UPT  (2.0e-9/3600.)  /* CO2  uptake by symbionts */
 #define SYMHCO3UPT (5.0e-9/3600.)  /* HCO3 uptake by symbionts */
 #define CO3UPT     (3.0e-9/3600.)  /*  direct calcification	*/
#endif
#define CAUPT CO3UPT
#define SYMHUPT	SYMHCO3UPT	    /* H    uptake by symbionts */

#define MIMECO2SYM/*----   2. MIME: set total carbon  uptake   -----*/
#ifdef  SYMTCIN
 #define SYMTCUPT	(SYMTCIN*1.e-9/3600.)
#else
 #define SYMTCUPT (7.2e-9/3600.)
#endif
#define VMAX (16.2e-9*.768/3600.)/* Michaelis-Menten max upt. for CO2 */
#define KS 5.0      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 200.  /* [mumol/kg] Bijma et al., 1992, p. 252 */

#endif

#ifdef FORAMOUL2	/* Orbulina Universa, Light (Steffi's data) */
	/*					*/
        /* realistic values :			*/
	/* 1. photos. of symbionts		*/
	/*    10. nmol/h 			*/
	/* 2. respiration of foram		*/
	/*   -2.5 nmol/h 			*/
	/* 3. calcification			*/
	/*    2.2 nmol/h			*/

#define RADIUS 277.0   		/* radius of foram [mu] diatoms     */
				/* Bulk radius, see BULKFAC	    */
				/* Num. of mesh points, see #define M */

/*	UPTAKE at the shell : 	(0.0e-9/3600.)			*/
#define CO2UPT  (-3.78e-9/3600.)	/* [mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT (0.0e-9/3600.)
#define OHUPT 0.0
#define HUPT  0.0		/* 0.0 				*/
#define BOH4UPT	0.0
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 	8.2		/* >= 7.			*/
#define UDICBULK 2400. 		/* 2010.			*/
#define UALKBULK 2819.	 	/* 2750.			*/

/*	UPTAKE of symbionts : 					*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#ifdef  SYMCO2IN
 #define SYMCO2UPT	(SYMCO2IN*1.e-9/3600.)
 #define SYMHCO3UPT	(SYMHCO3IN*1.e-9/3600.)
 #define CO3UPT  	(CO3IN*1.e-9/3600.)
#else
 #define SYMCO2UPT  (0.0e-9/3600.)  /* CO2  uptake by symbionts */
 #define SYMHCO3UPT (0.0e-9/3600.)  /* HCO3 uptake by symbionts */
 #define CO3UPT     (3.0e-9/3600.)  /*  direct calcification	*/
#endif
#define CAUPT CO3UPT
#define SYMHUPT	SYMHCO3UPT	    /* H    uptake by symbionts */

#define MIMECO2SYM/*----   2. MIME: set total carbon  uptake   -----*/
#ifdef  SYMTCIN
 #define SYMTCUPT	(SYMTCIN*1.e-9/3600.)
#else
 #define SYMTCUPT (9.61e-9/3600.)
#endif
#define VMAX (12.44e-9/3600.)   /* Michaelis-Menten max upt. for CO2 */
#define KS 5.0      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	2.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+100.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 200.  /* [mumol/kg] Bijma et al., 1992, p. 252 */

#endif


#ifdef FORAMGSD	/* Globigerinoides Sacculifer, Joergensen 85 */
     	/*		DARK:			*/
        /* 					*/
	/* 1. photos. of symbionts		*/
	/*    0. nmol/h 			*/
	/* 2. respiration of foram		*/
	/*   -2.7  nmol/h 			*/
	/* 3. calcification			*/
	/*    0.4 nmol/h 				*/

#define RADIUS 200.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/

/*	UPTAKE at the shell : 					*/
#define CO3UPT   (0.0e-9/3600.)  /* .4  direct calcification	*/
#define CO2UPT  (-1.3e-9/3600.)	/* -0.9[mol CO2 / s / foram] respiration */
				/* -1.3 with -0.4 from HCO3- upt for calc*/
#define TCO2UPT CO2UPT
#define HCO3UPT (0.4e-9/3600.)
#define CAUPT CO3UPT
#define OHUPT 0.0
#define HUPT  0.0
#define BOH4UPT	0.0
#define REDF	1.0		/* Redfield ratio foram O2=R*CO2 resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 8.25		/* 8.25				*/
#define DICBULK 2167. 		/* 2167.			*/

/*	UPTAKE of symbionts : 					*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(-1.8e-9/3600.) /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake -------*/
#define SYMTCUPT (0.0e-9/3600.)
#define VMAX (16.2e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 1.5      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Jorgensen, 1985 */

#endif



#ifdef FORAMGSL1	/* Globigerinoides Sacculifer, Joergensen 85 */
     	/*		LIGHT:			*/
        /* 					*/
	/* 1.  gross photos. of symbionts	*/
	/*     18   nmol/h 			*/
	/*     (net photos. system: 15 nmol/h)	*/
	/* 2a. respiration of foram		*/
	/*    -1.2  nmol/h 			*/
	/* 2b. respiration of symbionts		*/
	/*    -1.8  nmol/h 			*/
	/* 3.  calcification			*/
	/*     3.25 nmol/h 			*/

#define RADIUS 200.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/

/*	UPTAKE at the shell : 					*/
#define CO3UPT   (3.25e-9/3600.)  /* 3.25 direct calcification	*/
#define CO2UPT  (-1.2e-9/3600.)/* -1.2 [mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT (0.0e-9/3600.)  /* (0.0e-9/3600.) */
#define O2UPT (-CO2UPT)
#define CAUPT   CO3UPT
#define OHUPT 0.0
#define HUPT  0.0
#define BOH4UPT	0.0 		/* (0.0e-9/3600.*6.8e-5) */
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 8.25		/* 8.25				*/
#define DICBULK 2167. 		/* 2167.			*/

/*	UPTAKE of symbionts : 16.2 = 18 (gross) - 1.8 (symb. resp.)*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.)  /* CO2  uptake by symbionts (2.0)*/
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts (14.2)*/
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define MIMECO2SYM/*----   2. MIME: set total carbon  uptake  ------*/
#define SYMTCUPT (16.2e-9*.768/3600.)/* Redfield CO2:O2 106:138 = 0.768	 */
#define VMAX (16.2e-9*.768/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 5.      		/* [mumol/l] Michaelis-Menten half sat.  */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution*/
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	 */
#define O2BULK 210.  /* [mumol/kg] Joergensen, 1985 */

#endif

#ifdef FORAMGSL	/* Globigerinoides Sacculifer, Joergensen 85 */
     	/*		LIGHT:			*/
        /* 					*/
	/* 1.  gross photos. of symbionts	*/
	/*     18   nmol/h 			*/
	/*     (net photos. system: 15 nmol/h)	*/
	/* 2a. respiration of foram		*/
	/*    -1.2  nmol/h 			*/
	/* 2b. respiration of symbionts		*/
	/*    -1.8  nmol/h 			*/
	/* 3.  calcification			*/
	/*     3.25 nmol/h 			*/

#define RADIUS 200.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/

/*	UPTAKE at the shell : 					*/
#define CO3UPT   (0.0e-9/3600.)  /* 0 direct calcification	*/
#define CO2UPT  (-4.45e-9/3600.) /* (1.2+3.25) */
				 /* -4.45 with -3.25 from HCO3- upt for calc*/
#define TCO2UPT CO2UPT
#define HCO3UPT (6.5e-9/3600.)
#define O2UPT (-CO2UPT)
#define CAUPT   CO3UPT
#define OHUPT 0.0
#define HUPT  0.0
#define BOH4UPT	0.0 		/* (0.0e-9/3600.*6.8e-5) */
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 8.25		/* 8.25				*/
#define DICBULK 2167. 		/* 2167.			*/

/*	UPTAKE of symbionts : 16.2 = 18 (gross) - 1.8 (symb. resp.)*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.)  /* CO2  uptake by symbionts (2.0)*/
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts (14.2)*/
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define MIMECO2SYM/*----   2. MIME: set total carbon  uptake  ------*/
#define SYMTCUPT (16.2e-9*.768/3600.)/* Redfield CO2:O2 106:138 = 0.768	 */
#define VMAX (16.2e-9*.768/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 5.      		/* [mumol/l] Michaelis-Menten half sat.  */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution*/
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	 */
#define O2BULK 210.  /* [mumol/kg] Joergensen, 1985 */

#endif


#ifdef FORAMRES /* Respiration experiment for paper 	*/
 	/*  respiration of foram+symbionts -3.0  nmol/h */


#define RADIUS 200.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/

/*	UPTAKE at the shell : 					*/
#define CO3UPT   (0.0e-9/3600.)  /* .75/3.25 direct calcification	*/
#define CO2UPT  (-1.3e-9/3600.)	/* -1.0[mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT 0.0
#define CAUPT CO3UPT
#define OHUPT 0.0
#define HUPT  0.0
#define BOH4UPT	0.0
#define REDF	1.0		/* Redfield ratio foram O2=R*CO2 resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 8.25		/* 8.25				*/
#define DICBULK 2167. 		/* 2167.			*/

/*	UPTAKE of symbionts : 					*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(-1.7e-9/3600.) /* -1.7 CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake -------*/
#define SYMTCUPT (0.0e-9/3600.)
#define VMAX (16.2e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 1.5      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Jorgensen, 1985 */

#endif

#ifdef FORAMRES2 /* Respiration experiment for diss. 	*/
 	/*  respiration of foram -2.0  nmol/h */
#define DISS
#define RADIUS 200.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/

/*	UPTAKE at the shell : 					*/
#define CO3UPT   (0.0e-9/3600.)  /* .75/3.25 direct calcification	*/
#define CO2UPT  (-2.0e-9/3600.)	/* -1.0[mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT 0.0
#define CAUPT CO3UPT
#define OHUPT 0.0
#define HUPT  0.0
#define BOH4UPT	0.0
#define REDF	1.0		/* Redfield ratio foram O2=R*CO2 resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 8.2		/* 8.25				*/
#define DICBULK 1820. 		/* 2167.			*/

/*	UPTAKE of symbionts : 					*/
#define USYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.) /* -1.7 CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake -------*/
#define SYMTCUPT (0.0e-9/3600.)
#define VMAX (16.2e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 1.5      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Jorgensen, 1985 */

#endif


#ifdef FORAMPHS /* Photosynthesis experiment for paper 	*/
 	/*  photosynthesis of symbionts 12.5  nmol/h */


#define RADIUS 200.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/

/*	UPTAKE at the shell : 					*/
#define CO3UPT   (0.0e-9/3600.)  /* .75/3.25 direct calcification	*/
#define CO2UPT  (0.0e-9/3600.)	/* -1.0[mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT 0.0
#define CAUPT CO3UPT
#define OHUPT 0.0
#define HUPT  0.0
#define BOH4UPT	0.0
#define REDF	1.0		/* Redfield ratio foram O2=R*CO2 resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 8.25		/* 8.25				*/
#define DICBULK 2167. 		/* 2167.			*/

/*	UPTAKE of symbionts : 					*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.) /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define MIMECO2SYM/*----   2. MIME: set total carbon  uptake -------*/
#define SYMTCUPT (16.2e-9/3600.)
#define VMAX (16.2e-9*.768/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 5.      		/* [mumol/l] Michaelis-Menten half sat.  */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Jorgensen, 1985 */

#endif

#ifdef FORAMCLC /* Calcification experiment for paper 	*/
 	/*  calcfication of forams 3.25  nmol/h */


#define RADIUS 200.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/

/*	UPTAKE at the shell : 					*/
#define CO3UPT   (3.25e-9/3600.)  /* 3.25 direct calcification	*/
#define CO2UPT   (0.0e-9/3600.) /* -1.0[mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT 0.0
#define CAUPT CO3UPT
#define OHUPT 0.0
#define HUPT  0.0
#define BOH4UPT  0.0
#define REDF	1.0		/* Redfield ratio foram O2=R*CO2 resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 8.25		/* 8.25				*/
#define DICBULK 2167. 		/* 2167.			*/

/*	UPTAKE of symbionts : 					*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.) /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake -------*/
#define SYMTCUPT (0.0e-9/3600.)
#define VMAX (10.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 1.5      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Jorgensen, 1985 */

#endif

#ifdef FORAMCLC2 /* Calcification experiment for paper 	*/
 		/*  calcfication of forams 6.5 HCO3- nmol/h HCO3- */


#define RADIUS 200.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/

/*	UPTAKE at the shell : 					*/
#define CO3UPT   (0.0e-9/3600.)  /* 0.0 direct calcification	*/
#define CO2UPT   (-3.25e-9/3600.) /* -3.25[mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT (6.5e-9/3600.)
#define CAUPT CO3UPT
#define OHUPT 0.0
#define HUPT  0.0
#define BOH4UPT	0.0
#define REDF	1.0		/* Redfield ratio foram O2=R*CO2 resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 8.25		/* 8.25				*/
#define DICBULK 2167. 		/* 2167.			*/

/*	UPTAKE of symbionts : 					*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.) /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake -------*/
#define SYMTCUPT (0.0e-9/3600.)
#define VMAX (10.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 1.5      		/* [mumol/l] Michaelis-Menten half sat. */
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Jorgensen, 1985 */

#endif


#ifdef FORAMSYM1 /* Symbionts in halo around shell (500/50 um) 	*/

#define RADIUS 300.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/

/*	UPTAKE at the shell : 					*/
#define CO3UPT   (1.5e-9/3600.)  /* 3. direct calcification	*/
#define CO2UPT   (-2.1e-9/3600.)  /*  [mol CO2 / s / foram] respiration  */
#define TCO2UPT CO2UPT
#define HCO3UPT (0.0e-9/3600.)
#define CAUPT CO3UPT
#define OHUPT 0.0
#define HUPT  HCO3UPT
#define BOH4UPT	0.0
#define REDF	1.0		/* Redfield ratio foram O2=R*CO2 resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 	8.1		/* >= 7.			*/
#define UDICBULK 2010. 		/* 2010.			*/
#define UALKBULK 2400.	 	/* 2750.			*/


/*	UPTAKE of symbionts : 					*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.) /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define MIMECO2SYM/*----   2. MIME: set total carbon  uptake -------*/
#define SYMTCUPT (7.0e-9/3600.)
#define VMAX (10.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 5.      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+100.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Jorgensen, 1985 */

#endif

#ifdef FORAMSYM2 /* Symbiont halo 100um from shell Exp. for paper */

#define RADIUS 300.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/

/*	UPTAKE at the shell : 					*/
#define CO3UPT   (1.5e-9/3600.)  /* .75/3.25 direct calcification	*/
#define CO2UPT   (-2.1e-9/3600.)  /*  [mol CO2 / s / foram] respiration  */
#define TCO2UPT CO2UPT
#define HCO3UPT (0.0e-9/3600.)
#define CAUPT CO3UPT
#define OHUPT 0.0
#define HUPT  HCO3UPT
#define BOH4UPT	0.0
#define REDF	1.0		/* Redfield ratio foram O2=R*CO2 resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 	8.1		/* >= 7.			*/
#define UDICBULK 2010. 		/* 2010.			*/
#define UALKBULK 2400.	 	/* 2750.			*/

/*	UPTAKE of symbionts : 					*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.) /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define MIMECO2SYM/*----   2. MIME: set total carbon  uptake -------*/
#define SYMTCUPT (7.0e-9/3600.)
#define VMAX (10.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 5.      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRADMIN (RADIUS+100.)   /* inner radius for symbiont distrib.*/
#define SYMRAD (SYMRADMIN+500.0)   /* outer radius for symbiont distrib.*/
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Jorgensen, 1985 */

#endif

#ifdef FORAMR /* Vary radius and P,C,R 	*/

#define RADIUS 200.0   	/* 50.0 200.0 300.0 radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/

/*	UPTAKE at the shell : 					*/
#define CO3UPT   (2.0e-9/3600.)  /* 0.5 2.0 3.0 direct calcification	*/
#define CO2UPT   (-1.5e-9/3600.)  /* 0.5 1.5 2.1 [mol CO2 / s / foram] respiration  */
#define TCO2UPT CO2UPT
#define HCO3UPT (0.0e-9/3600.)
#define CAUPT CO3UPT
#define OHUPT 0.0
#define HUPT  HCO3UPT
#define BOH4UPT	0.0
#define REDF	1.0		/* Redfield ratio foram O2=R*CO2 resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 	8.1		/* >= 7.			*/
#define UDICBULK 2010. 		/* 2010.			*/
#define UALKBULK 2400.	 	/* 2750.			*/


/*	UPTAKE of symbionts : 					*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.) /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define MIMECO2SYM/*----   2. MIME: set total carbon  uptake -------*/
#define SYMTCUPT (4.0e-9/3600.) /* 1.0 4.0 7.0			*/
#define VMAX (10.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 5.      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+300.0)   /* 50.0 300.0 500.0 */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Jorgensen, 1985 */

#endif


#ifdef AGG
#define AGGRADIUS 1300.   	/* radius of aggregate  */
#define DIFFAGG 1.0		/* reduce diffusion in aggregate */
#define RADIUS 50.		/* inner boundary	*/
#define AGGVOL (4.*PI*KU(AGGRADIUS)/3.) /* agg. volume (mum3)*/
				/* uptake volume, see agguptvol */
#define CO3CRIT 60.		/* mu mol kg-1 */
#define KDISS (501./24./3600./100.)	/* diss. constant (frac per s) */
#define NU 3.			/* reaction order */

#define REHUX 3.75		/* radius Ehux cell (mum) */
#define CEHUX 7.e-13		/* carbon Ehux cell (mol C) */
#define TIME (20.*24.*3600)	/* current time (s)*/
#define TAU  (10.*24.*3600)	/* degradation time const. (s)*/

	/*	UPTAKE at the shell : 					*/
#define CO3UPT   (0.0e-9/3600.)  /* .75/3.25 direct calcification	*/
#define CO2UPT  (-0.0e-9/3600.)  /* [mol CO2 / s / ] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT (0.0e-16)
#define O2UPT (-0.0)
#define CAUPT   CO3UPT
#define OHUPT 0.0
#define HUPT  (0.0e-20)
#define BOH4UPT	(0.0e-9/3600.)
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 7.8094		/* 7.7299 7.8094 7.9887 */
#define DICBULK 2200. 		/* 2167. 2200 2323.3 2181.2	*/
#define UALKBULK 2400.	 	/* 2723.			*/

/*	UPTAKE of symbionts : 16.2 = 18 (gross) - 1.8 (symb. resp.)*/
#define SYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define USYMCO2UPT	(-0.0e-9/3600.)  /* CO2 release, symco2upt */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define USYMCO3UPT	(-0.0e-9/3600.)/* CO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake     ---*/
#define SYMTCUPT (3.e-9/3600.)
#define VMAX (10.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 1.5      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (AGGRADIUS)      /* outer radius for symbiont distribution */
#define USYMO2UPT 0.0
#define REDAGG (138./106.)      /* Redfield ratio Corg oxidation	*/
#define O2BULK 35.  /* [mumol/kg] Joergensen, 1985 */

#endif


#ifdef AGG1
#define RADIUS 45.   		/* radius of aggregate */
				/* Bulk radius, see BULKFAC	*/
#define CO3CRIT 60.		/* mu mol kg-1 */
#define KDISSK 5.		/* per day */
#define NU 4.5

	/*	UPTAKE at the shell : 					*/
#define CO3UPT   (-0.0e-9/3600.)  /* .75/3.25 direct calcification	*/
#define CO2UPT  (-0.0083e-9/3600.)/* [mol CO2 / s / foram] respiration */
#define TCO2UPT CO2UPT
#define HCO3UPT (0.0e-16)
#define O2UPT (-CO2UPT)
#define CAUPT   CO3UPT
#define OHUPT 0.0
#define HUPT  (0.0e-20)
#define BOH4UPT	(0.0e-9/3600.)
#define REDF	1.0		/* Redfield ratio O2=R*CO2 foram resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 7.8844449	/* 8.25 7.9716 7.8601 7.765 7.8844449	*/
#define DICBULK 2181.2		/* 2167. 2200 2323.3		*/
#define UALKBULK 2400.	 	/* 2723.			*/

/*	UPTAKE of symbionts : 16.2 = 18 (gross) - 1.8 (symb. resp.)*/
#define USYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define SYMCO2UPT	(0.0e-9/3600.)  /* CO2  uptake by symbionts */
#define SYMHCO3UPT	(0.0e-9/3600.)	/* HCO3 uptake by symbionts */
#define SYMHUPT	SYMHCO3UPT		/* H    uptake by symbionts */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake     ---*/
#define SYMTCUPT (3.e-9/3600.)
#define VMAX (10.0e-9/3600.)    /* Michaelis-Menten max upt. for CO2 */
#define KS 1.5      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRAD (RADIUS+500.0)   /* outer radius for symbiont distribution */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Joergensen, 1985 */

#endif



#ifdef DIATOM1	/* model internal of algal cell (chloroplast) */
#define CLPL
#define DRAIN
#define VMAXDRAIN (3.0e-16)
#define KSDRAIN    5.0
#define DC13DIC 0.0	    /* delta C13 of total CO2		*/
#define EPSPCO2UPT (-27.0)  /* delta value of respired CO2	*/
#define RADIUS 0.1   	    /* 0.1 inner boundary in the center of the cell */
			    /* Bulk radius, see BULKFAC		*/
			    /* note BULKFAC = ?, M = ? !!!!!!!! */
#define CELLRADIUS 10.	    /* radius of cell			*/
#define MEMBTHK	1.	    /* thickness of cell membrane	*/
#define PERMCO2 1.e2	    /* 1.e2 permeability of cell membrane CO2 (?m/s)*/
#define PERMBOH3 1.e2	    /* permeability of cell membrane for BOH3	*/
#define PERMION 1.e1	    /* permeability of cell membrane for ions	*/
#define REDDCO2 (PERMCO2*MEMBTHK*2.)	    /* (PERMCO2*MEMBTHK*2.)	*/
/*	UPTAKE at the cell wall : 					*/
#define CO3UPT   (0.0e-9/3600.)  /* 	*/
#define CO2UPT  (0.0e-16)	/*  */
#define TCO2UPT CO2UPT
#define HCO3UPT 0.0
#define CAUPT CO3UPT
#define OHUPT 0.0
#define HUPT  0.0
#define REDF	1.0		/* Redfield ratio foram O2=R*CO2 resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 8.2 		/* 8.25				*/
#define DICBULK 2000. 		/* 2167.			*/

/*	UPTAKE within the cell : 				*/
#define SYMBIONTS/*----     1. set CO2 uptake           --------*/
#define SYMCO2UPT	(5.0e-16)       /* CO2  uptake in chloroplast */
#define SYMHCO3UPT	(0.0e-16)	/* HCO3 uptake  */
#define SYMHUPT SYMHCO3UPT		/* H    uptake  */
#define SYMHUPT 0.0			/* H    uptake  */

#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake -------*/
#define SYMTCUPT (5.0e-16) /* 1.0 4.0 7.0			*/
#define VMAX (10.e-16)    /* Michaelis-Menten max upt. for CO2 */
#define KS 1.5      		/* [mumol/l] Michaelis-Menten half sat. */
#define DVDIT	3.0		/* increase of uptake per iteration	*/
#define SYMRADMIN 4.	        /* inner bound. of chloroplast (mu m)*/
#define SYMRAD 6.		/* outer bound. of chloroplast (mu m) */
#define SYMO2UPT (-(SYMCO2UPT+SYMHCO3UPT))+
#define REDS 1.0		/* Redfield ratio symbiont photosynth.	*/
#define O2BULK 210.  /* [mumol/kg] Jorgensen, 1985 */
#endif


#ifdef DIATOMW	/* Ulf */
#define DC13DIC 0.0	/* delta C13 of total CO2		*/
#define EPSPCO2UPT (-4.0)  /* delta value of respired CO2		*/
#define RADIUS 10.0   		/* radius of foram [mu] diatoms */
				/* Bulk radius, see BULKFAC	*/

/*	UPTAKE at the shell : 					*/
#define CO3UPT   (0.0e-9/3600.)  /* .4/3.25 direct calcification	*/
#define CO2UPT  (5.0e-16)	/*  */
#define TCO2UPT CO2UPT
#define HCO3UPT 0.0
#define CAUPT CO3UPT
#define OHUPT 0.0
#define HUPT  0.0
#define REDF	1.0		/* Redfield ratio foram O2=R*CO2 resp.*/
				/* see: O2UPT at the shell	*/
#define PHBULK 8.2 		/* 8.25				*/
#define DICBULK 2000. 		/* 2167.			*/

/*	UPTAKE of symbionts : 					*/
#define USYMBIONTS/*----     1. set CO2 and HCO3- uptake     --------*/
#define UMIMECO2SYM/*----   2. MIME: set total carbon  uptake -------*/
#define O2BULK 210.  /* [mumol/kg] Jorgensen, 1985 */

#endif


#ifdef DIATOM3
#define RADIUS 20.0   /* radius of algae [mu] diatoms */
#define CO2UPT (0.e-15)  /* def=1.e-15; 4.e-15 [mol/s*C/cell] for r=20 mu */
#define TCO2UPT CO2UPT   /* 5e-17[mol/s * C/cell] */
#define HCO3UPT (1.e-16)
#define CO3UPT 0.0
#define CAUPT 0.0
#define OHUPT 0.0
#define HUPT (1.e-16)
#define PHBULK 8.2       /* def = 8.2; CO2UPT=1.e-16 for pH=9.0 */
#define DICBULK 2167.
#endif




#ifdef COCCOW
#define RADIUS 4.0   /* radius of algae [mu] coccos */
#define PICO 20.0    /* uptake [pg C/12h/cell]  */
#define CO2UPT (PICO/12.e12/12./3600.)  /* [mol C/s/cell] */
#define HCO3UPT 0.0
#define UCO3UPT CO2UPT        /* Alk:C_org = 1:1 */
#define UCO3UPT 0.0           /* Alk:C_org = 0:1 */
#define CO3UPT (2.0*CO2UPT)   /* Alk:C_org = 2:1 */
#define TCO2UPT (CO2UPT+HCO3UPT+CO3UPT)
#define CAUPT   (CO3UPT+0.5*HCO3UPT)
#define OHUPT   0.0
#define HUPT    0.0
#define PHBULK 8.5
#define DICBULK 2167.
#endif



#ifdef CLPL
#define BULKFAC 1000   			/* 1000 */
#else
#define BULKFAC 10   			/* def. = 10 set outer radius */
#endif
#ifdef AGG
#define RBULK AGGRADIUS
#else
#define RBULK (RADIUS*(double)BULKFAC)  /* "bulk radius" [mu] */
#endif
#define CO2ATM 355.0                    /* [ppmv] */
#define CABULK (10.3*1000.)             /* [mumol/kg]
						Stumm and Morgan 1981 */

/* -----   Order of equations   ----- */


#ifdef C13ISTP
#ifndef CISTP
#define EQCO2   1	/* this is total [CO2]  = [12CO2]  + [13CO2]  */
#define EQHCO3  2	/* this is total [HCO3] = [H12CO3] + [H13CO3] */
#define EQCO3   3	/* this is total [CO3]  = [12CO3]  + [13CO3]  */
#define EQHP    4
#define EQOH    5
#define EQCCO2  6	/* 13CO2  */
#define EQHCCO3 7	/* 13HCO3 */
#define EQCCO3  8	/* 13CO3  */
#define EQBOH3  9
#define EQBOH4 10
#define EQO2   11
#define EQCA   12
#endif
#ifdef CISTP
#define EQCO2   1	/* this is [12CO2]  */
#define EQHCO3  2	/* this is [H12CO3] */
#define EQCO3   3	/* this is [12CO3]  */
#define EQHP    4
#define EQOH    5
#define EQCCO2  6	/* 13CO2  */
#define EQHCCO3 7	/* 13HCO3 */
#define EQCCO3  8	/* 13CO3  */
#define EQBOH3  9
#define EQBOH4 10
#define EQO2   11
#define EQCA   12
#endif
#else
#ifndef BORISTP
#define EQCO2  1
#define EQHCO3 2
#define EQCO3  3
#define EQHP   4
#define EQOH   5
#define EQBOH3 6
#define EQBOH4 7
#define EQO2   8
#define EQCA   9
#endif
#ifdef BORISTP
#define EQCO2  1
#define EQHCO3 2
#define EQCO3  3
#define EQHP   4
#define EQOH   5
#define EQBOH3 6		/* 10B+11B or 10B alone   */
#define EQBOH4 7
#define EQBBOH3   8		/* 11B 	(80% of total B3) */
#define EQBBOH4   9
#define EQO2    11
#define EQCA    12
#endif
#endif

/* -----   parameter ----- */

#define RGASC 1.987                /* [cal/K/mol] gas constant */
#define RGASJ 8.3145               /* [J/K/mol] gas constant */



/* Added 25/01/2017 by Branson and Holland! */
#if defined PYTHON
#define SALINITY 3.330000e+01
#define TEMP 2.200000e+01
#endif

#ifdef  AGG
#define SALINITY 35.0
#define TEMP 4. 			   /* 13.34  */
#endif

#ifdef  FORAMW
#define SALINITY 35.
#define TEMP 25.
#endif

#if defined (FORAMBORD2) || defined (FORAMBORL2)
#define SALINITY 33.3
#define TEMP 22.
#endif

#if defined (FORAMOUD) || defined (FORAMOUL) /* O. Universa	*/
#define SALINITY 38.	     	    /* (Catalina Island)*/
#define TEMP 28.
#endif

#if defined (FORAMGSD) || defined (FORAMGSL) /* G. Sacculifer		*/
#define SALINITY 40.7		  	     /* Joergensen 85		*/
#define TEMP 24.5
#endif

#if defined (FORAMRES) || defined (FORAMPHS) || defined (FORAMCLC)
#ifndef DISS
#define SALINITY 40.7		  	     /* experiments for paper	*/
#define TEMP 24.5
#endif
#ifdef DISS
#define SALINITY 35.		  	     /* experiments for diss	*/
#define TEMP 20.
#endif
#endif


#if defined (FORAMSYM1) || defined (FORAMSYM2)
#define SALINITY 33.5
#define TEMP 20.
#endif

#ifdef FORAMR
#define SALINITY 33.5
#define TEMP 20.
#endif


#ifdef  DIATOMW
#define SALINITY 32.0
#define TEMP 15.0
#endif

#ifdef  COCCOW
#define SALINITY 35.0
#define TEMP 15.0
#endif

/* Added 25/01/2017 by Branson and Holland! */
#ifdef PYTHON
#define BORTBULK 2.055086e+03 /* [mumol/kg], DOE94  */
#else
#define BORTBULK (432.*(SALINITY/35.)) /* [mumol/kg], DOE94	 */
#endif

#define TKELVIN 273.15
#define EPS0 8.8542e-12
#define KBOLTZ 1.38e-23
#define ECHARGE 1.602e-19
#define TK (TEMP + TKELVIN)



#ifdef C13ISTP

#ifndef BORON
#ifndef OXYGEN
#ifndef CALCIUM
#define  N2 8	/* number of second order equations */
#endif
#endif
#endif


#ifdef BORON
#ifndef OXYGEN
#ifndef CALCIUM
#define  N2 10	/* number of second order equations */
#endif
#endif
#endif

#ifdef BORON
#ifdef OXYGEN
#ifndef CALCIUM
#define  N2 11	/* number of second order equations */
#endif
#endif
#endif

#ifdef BORON
#ifdef OXYGEN
#ifdef CALCIUM
#define  N2 12	/* number of second order equations */
#endif
#endif
#endif

#else		/* else C13ISTP			    */

#ifndef BORON
#ifndef OXYGEN
#ifndef CALCIUM
#define  N2 5	/* number of second order equations */
#endif
#endif
#endif


#ifdef BORON
#ifndef OXYGEN
#ifndef CALCIUM
#ifdef BORISTP
#define  N2 9	/* number of second order equations */
#else
#define  N2 7	/* number of second order equations */
#endif
#endif
#endif
#endif

#ifdef BORON
#ifdef OXYGEN
#ifndef CALCIUM
#define  N2 8	/* number of second order equations */
#endif
#endif
#endif

#ifdef BORON
#ifdef OXYGEN
#ifdef CALCIUM
#define  N2 9	/* number of second order equations */
#endif
#endif
#endif


#endif		/* end C13ISTP			    */


#define NE (2*N2)        /* number of first order equations 	*/
#define NLB N2           /* number of left boundary conditions 	*/
#define NRB (NE-NLB)     /* number of right boundary conditions */
#define NB NLB           /* number of left boundary conditions 	*/
#define NR (NE-NB)       /* number of right boundary conditions */
#define M (BULKFAC*100)  /* number of mesh points (BLKFC = 10 or 20)*/
			 /* set M = 10*100 or 20*50		*/
#ifdef CLPL
#define M (int)(BULKFAC*4)
#endif
#define NSJ (2*NE+1)
#define NCJ (NE-NB+1)
#define NCK (M+1)

#define IMAX M       /* number of points (discretization) */
#define IMAX1 (IMAX-1)
#define IMAX2 (IMAX-2)

#define CONV 5.0e-6  /* convergence criterion def. = 5.0e-6 */
#define SLOWC 1.0    /* fraction of correction def. = 1.0 */
#define ITMAX 30     /* def. = 30 max. number of iterations */

                   /*  Diatom-Michaelis-Menten for CO2 */
/* #define MIMECO2DIA  */
#define VMAXDIA 0.2    /* 0.2 [mol/kg/mu] Michaelis-Menten */
#define KSDIA 2.0      /* [mumol/kg] Michaelis-Menten half sat. */

/* ===================== global (begin) ==================== */

FILE *fppara,*fpr,*fpanasol,
     *fpco2,*fpco3,*fph,*fphco3,*fpoh,
     *fpdco2,*fpdco3,*fpdh,*fpdhco3,*fpdoh,
     *fpscale,*fpequi,*fpks,
     *fpdc
#ifdef C13ISTP
    ,*fpcco2,*fphcco3,*fpcco3,*fpdcco2,*fpdhcco3,*fpdcco3,*fpdc13s,
     *fpdcc
#endif
#ifdef OXYGEN
    ,*fpo2,*fpdo2
#endif
#ifdef BORON
    ,*fpboh3,*fpboh4,*fpdboh3,*fpdboh4
#ifdef BORISTP
    ,*fpbboh3,*fpbboh4,*fpdbboh3,*fpdbboh4
#endif
#endif
#ifdef CALCIUM
    ,*fpca,*fpdca
#endif

#ifdef FLSYMUPT
    ,*fpsyup
#ifdef C13ISTP
    ,*fpsyup13
#endif
#endif

#ifdef CBNSLOOP
    ,*fpdc13slp
#ifdef READ
    ,*fpread
#endif
#endif
#ifdef CLPL
    ,*fpdifcofm
#endif
     ;

int debug02=0,ir,nsymrad
#if defined (FORAMSYM2) || defined (CLPL)
     ,nsymradmin
#endif
#ifdef CLPL
     ,ncmembmin,ncmembmax,dumk1,dumk2,calldifeq
#endif
#ifdef CBNSLOOP
    ,icl,n
#endif
     ;

double
       aux1,aux2,aux3,dt,h2obulk,surface,
       co2bulk,co3bulk,hco3bulk,dicbulk,dr,phbulk,hbulk,hco3bulk,ohbulk,
       co2flux,hflux,hco3flux,co3flux,ohflux, /* fluxes (left boundary) */
       co2fluxs,hfluxs,hco3fluxs,co3fluxs,ohfluxs, /* fluxes (save) */
       data,dco2,dh,doh,dhco3,dco3,      /* diffusion coefficients */
       h,hh,alkbulk,alkbulkd,k1d,k2d,scratch1,scratch2,
       km1s,
       k12,k21,k13,k31,k23,k32,kp4,km4,kp5h,km5h,km6,kp6,
       kp1s,Kh2co3,kh2co3,kw,ohminus,dummy,tmp,tmp1,tmp2,tmp3,tmp4,
       alk[M+1],ata[M+1],ef[M+1],
       co2[M+1],hplus[M+1],hco3[M+1],co3[M+1],oh[M+1],h2o[M+1],
       r[M+1]
#ifdef C13ISTP
      ,cco2[M+1],cco2bulk,cco2conv1,dcco2,cco2flux,cco2fluxs
      ,hcco3[M+1],hcco3bulk,hcco3conv1,dhcco3,hcco3flux,hcco3fluxs
      ,cco3[M+1],cco3bulk,cco3conv1,dcco3,cco3flux,cco3fluxs
      ,kp1scc,km1scc,kp4cc,km4cc,kp5cc,km5cc,tmp1cc,tmp2cc,tmp4cc
      ,alpha,alpha5,alphac,alphahc
      ,eps1,eps2,eps3,eps4,eps5,eps6	/* fractionation  between
				           different carbon species.
			 	      (see initeps() for explanation) */
      ,d13co2bulk,d13hco3bulk,d13co3bulk
      ,d13co2,d13hco3
      ,cco2upt,hcco3upt,cco3upt,symcco2upt,symhcco3upt
#ifdef DIATOMW
      ,d13co2upt
#endif
#endif

#ifdef CBNS
      ,phbulkarg,dicbulkarg,alkbulkarg
#endif
#ifdef CBNSLOOP
      ,phbulkv
 #ifdef READ
      ,phin[NDAT],dicin[NDAT],alkin[NDAT],d,symtcup,co3uptldat
 #endif
#endif
#ifdef OXYGEN
      ,o2[M+1],o2bulk,o2conv1,do2,o2flux,o2fluxs
#endif
#ifdef BORON
      ,boh3[M+1],boh4[M+1],boh3bulk,boh4bulk,dboh3,dboh4,
       boh3flux,boh4flux,boh3fluxs,boh4fluxs,
       kp7,km7,Kbor,kbdum
#ifdef BORONRC4
       ,bortmp = 88
#endif
#ifdef BORISTP
      ,bboh3[M+1],bboh4[M+1],bboh3bulk,bboh4bulk,dbboh3,dbboh4,
       bboh3flux,bboh4flux,bboh3fluxs,bboh4fluxs,
       kp7bb,km7bb,
       d11boh3bulk,d11boh4bulk,btmp,alphab,alphabp,d11boh4upt,
       K10bor,K11bor,tmpb
#endif
#endif
#ifdef CALCIUM
      ,ca[M+1],cabulk,dca,caflux,cafluxs
#endif
#ifdef CLPL
      ,difcofm[N2+1][M+1],fdrain,dumdr
#endif
#ifdef AGG
      ,c0,co2upt,symco2upt,symo2upt,agguptvol
#endif
      ;


#ifdef MIMECO2SYM
      double dummyd,vmaxco2=VMAX,vmaxit
#ifdef C13ISTP
      ,d13co2phyt,d13co2phytkm1,d13hco3phyt,d13hco3phytkm1,dumf,dumfkm1
      ,z,zkm1,dz_dx
#endif
      ;
#endif
      int co2negflag = 0;


/* ===================== global (end)   ==================== */




/* =========================================================
   =========================================================

               ''chemical'' functions (begin)

   =========================================================
   ========================================================= */




double chlorinity(double s)
{
/* =============================================================
   Stumm and Morgan [1981, p. 567]
   check value: chlorinity(34.32445) = 19.0
   ============================================================= */
   aux1 = s / 1.80655;
   return(aux1);
}   /* --- end of chlorinity --- */

double scl(double cl)
{
/* ============= salinity -> chlorinity =======================
   Stumm and Morgan [1981, p. 567]
   check value: chlorinity(34.32445) = 19.0
   ============================================================= */
   aux1 = 1.80655 * cl;
   return(aux1);
}   /* --- end of scl --- */

double ionicst(double cl)
{
/* =============================================================
   ionic strength
   Stumm and Morgan [1981, p. 204]
   ============================================================= */
   aux1 = 0.00147 + 0.03592 * cl + 0.000068 * cl * cl;
   return(aux1);
}   /* --- end of ionicst --- */

double fkdico(double tk,double z1, double z2,
              double d1, double d2,double ab)
{
/* ===========================================================

   Debye (1942)
     diffusion controlled reaction rates
   units: [m**3/mol/s] * 1.e3 -> [l/mol/s]

    DICO1 = E * E / (EPS0*epsrh2o) / KB     % 2.6248e-06 m K
    DICO2 = 4. * PI * NA * DICO1            % 1.9865e+19 m K/mol

   =========================================================== */

   double expfac,arg;

#define DICO1 2.6248e-06
#define DICO2 1.9865e+19

   arg = DICO1 * z1 * z2 / ab / tk;

#ifdef PRINT
   printf("%e   arg \n",arg);
#endif /* PRINT */

   expfac = 0.0;
   if(arg > (-20.)) expfac = exp(arg);
#ifdef PRINT
   printf("%e   expfac \n",expfac);
#endif /* PRINT */

   /* 1.e-12 [mu**2/s] -> [m**2/s] */

   arg = 1.0e-9 * DICO2 * z1 * z2 / tk * (d1 + d2) / (expfac - 1.0);

#ifdef PRINT
   printf("%e  %e  d1 d2 \n",d1,d2);
   printf("%e   fkdico \n",arg);
#endif /* PRINT */

/*   exit(1); */

   return(arg);
}

double fk21(double tk)
{
/* ===========================================================

   Eigen and Hammes (1963)
     k21 = 9.e6 [1/s] at 25 C

   Eact = 15 kcal/mol (same as for kp1s assumed!!!)

   =========================================================== */
   aux1 = EACTK / RGASC * (1. / (25. + TKELVIN) - 1. / tk);
   aux2 = 9.e6 * exp(aux1);

   return(aux2);
}

double fk31(double tk)
{
/* ===========================================================

   a) not known !!! -> 0.0
   b) k31' = 0.043 [1/s] Eigen et al. (1961)

   Eact = 15 kcal/mol (same as for kp1s assumed!!!)

   =========================================================== */
   aux1 = EACTK / RGASC * (1. / (25. + TKELVIN) - 1. / tk);
   aux2 = 0.043 * exp(aux1);

   return(aux2);
}

double fk13(double tk)
{
/* ===========================================================

   k_{-1}' in new nomenclature

   a) not known !!! -> 0.0
   b) k13' = 5.6e4 [l/mol/s] Eigen et al. (1961)

   Eact = 15 kcal/mol (same as for kp1s assumed!!!)

   =========================================================== */

   aux1 = EACTK / RGASC * (1. / (25. + TKELVIN) - 1. / tk);
   aux2 = 5.6e4 * exp(aux1);

   return(aux2);
}

double fk32(double tk)
{
/* ===========================================================

   Gavis and Ferguson (1975, p. 214)
     k32 = 0.03 [1/s]   (k_4 in their notation)

   Eact = 15 kcal/mol (same as for kp1s assumed!!!)

   =========================================================== */

   aux1 = EACTK / RGASC * (1. / (25. + TKELVIN) - 1. / tk);
   aux2 = 0.03 * exp(aux1);

   return(aux2);
}

double fkp4(double tk)
{
   double tk0,kp4s,a4,e4;

/* ===========================================================

   Johnson (1982), cf. Zeebe and Wolf-Gladrow, 2001

   k+4 = A4 * exp(-E4/RT)

   A4 = 4.70e7 kg/mol/s, E4 = 23.2e3 J/mol

   =========================================================== */

   a4   = 4.7e7*1.e-6; /* *1e-6 -> per mumol */
   e4   = 23.2e3;

   aux2 = a4*exp(-e4/RGASJ/tk);

#ifdef RC3

   /* Johnson 1982 gives kp4*kw = 13.4e-11 at T=25, S=33.77 	*/
   /* 		    and  kp4*kw = 15.2e-11 at T=25, S=37.06 	*/
   /* 		    =>   kp4*kw = 14.1e-11 at T=25, S=35.00 	*/
   /* on NBS scale. 	kw = 5.9299e-14  Hansson scale
   		    	kw = 4.6411e-14  Free   scale		*/
   /*  we assume NBS = Free					*/

   kp4s = 14.1e-11/4.6411e-14/1.0e6; /* [kg/mumol/s] (3038 kg/mol/s) */
   tk0  = 25. + TKELVIN;
   aux1 = EACTK / RGASC * (1. / tk0 - 1. / tk);
   aux2 = kp4s * exp(aux1);
#endif

   return(aux2);
}   /* --- end of fkp4 --- */



double fkp5h(double tk)
{

   double ab,arg;

/* ===========================================================

   Eigen (1964), cf. Zeebe and Wolf-Gladrow, 2001

   kp5h = 5.e10; kg/mol/s

   =========================================================== */

   arg = 5.e10;


#ifdef RC3
/* ===========================================================

   diffusion controlled -> O(1e10 l/mol/s)
   compare Knoche (1980)

   =========================================================== */

   ab = 4.0e-10;    /* 0.4 nm */

   arg = fkdico(tk,-1.0,1.0,dco3,dh,ab);

   arg = 1.e10;      /* Knoche (1980)  */
#endif

   return(arg);

}

double fkp6(double tk)
{
/* ===========================================================

   Eigen (1964), cf. Zeebe and Wolf-Gladrow, 2001
     kp6 = 1.4e-3 mol/kg/s

   =========================================================== */

   aux2 = 1.4e-3*1e6; /* *1e6 -> mu mol/kg/s */


#ifdef RC3
/* ===========================================================

   Atkins (1990,p. 796)
     kp6 = 2.4e-5 [1/s] * MOLH2O [mumol/kg]

   =========================================================== */

   aux1 = EACTK / RGASC * (1. / (25. + TKELVIN) - 1. / tk);
   aux2 = 2.4e-5*MOLH2O * exp(aux1);

#endif

   return(aux2);
}


#ifdef BORONRC4
double fkp7(double tk)
{
/* ===========================================================
   Rate coefficient of boric acid-borate system.


   Waton et al. (1984), cf. Zeebe and Wolf-Gladrow, 2001

		  kp7
   B(OH)3 + OH-   -->    B(OH)4-
		  <--
		  km7

   Kbor/Kw = [H+][B(OH)4-]/[B(OH)3]/Kw = kp7 / km7


   k+7 = A7 * exp(-E7/RT)

   A7 = 4.58e10 kg/mol/s
   E7 = 20.8e3   J/mol

   =========================================================== */
   double arg, a7, e7;

   a7 = 4.58e10;
   e7 = 20.8e3;

   arg = a7 * exp(-e7/RGASJ/tk); /* kg/mol/s */

   return(arg);
}
#endif


#ifdef BORONRC3
double fkm7(double tk)
{
/* ===========================================================
   Rate coefficient of boron.

		  kp7
   H2O + B(OH)3   -->   H+  +  B(OH)4-
		  <--
		  km7

   Kbor = [H+][B(OH)4-]/[B(OH)3] = kp7 / km7

   diffusion controlled -> km7 ~ O(1e10 l/mol/s)
   (pers. communication Millero)

   =========================================================== */

   double arg;


   arg = 1.e10;      /* [kg/mol/s]  */

#define UBORKTEST
#ifdef BORKTEST
   arg = 1.e7;      /* [kg/mol/s]  */
#endif

   return(arg);
}
#endif



double fkp1s(double s, double tk)
{

/* ===========================================================

   Johnson (1982)

     ln k_+1 (T,S) = 1246.98 - 6.19e4/tk - 183.0*log(tk)

   =========================================================== */

   aux2 = exp(1246.98 - 6.19e4/tk - 183.0*log(tk));


#ifdef RC3
/* ===========================================================

   Dieter's fit to Johnson (1982) ->
     k_CO2 (T,S) = (-4.7228e+08 * S + 6.8046e+10) * exp(-8.3467e+03/tk)

   =========================================================== */

   aux2 = (-4.7228e+08 * s + 6.8046e+10) * exp(-8.3467e+03/tk);

#endif

   return(aux2);

}   /* --- end of kp1s --- */

double fkh2co3(double tk)
{
   double ekh2co3,tk0,kh2co30;
/* ===========================================================

   Stumm and Morgan (1981, p.211)
     kh2co3 = 10 to 20 [1/s] at 20-25 C
     Eact = 16 kcal/mol

   Knoche (1980)
     kh2co3 = 28 \pm 3 [1/s] at 25 C

   =========================================================== */

#define KH2CKN80

#ifdef KH2CSM81
   kh2co30 = 10;
#endif

#ifdef KH2CKN80
   kh2co30 = 28;
#endif

   ekh2co3 = 1.6e4;
   tk0 = 25. + TKELVIN;
   aux1 = ekh2co3 / RGASC * (1. / tk0 - 1. / tk);
   aux2 = kh2co30 * exp(aux1);
   return(aux2);
}   /* --- end of kh2co3 --- */

double khenry(double s, double tk)
{
/* ==============================================================
   solubility of CO2 in seawater (Henry's law)
        [CO_2]_water = khenry       *    pCO_2
   units:[mol/kgiter]   [mol/kgiter/atm]   [atm]
   s [psu] salinity; tk [K] absolute temperature
   Stumm and Morgan [1981, p. 179 and 204]
   check value: -log K_H (s=34.32445 (cl=19.),t=25 C) = 1.53
   ============================================================== */
   double cl,is;
   cl = chlorinity(s);
   is = ionicst(cl);
   aux1 = pow(10.0,(2385.73 / tk - 14.0184 + 0.0152642 * tk
             - is * (0.28596 - 6.167e-4 * tk)));
   return(aux1);
}   /* --- end of khenry --- */


double K_S(double s, double tk)
{
/* ===========================================================

   % K_S  Dickson and Goyet (1994, Chapter 5, p. 13)
     neccessary for TOTAL2FREE()

     Equilibrium constant for HSO4- = H+ + SO4--

     K_S  = [H+]free [SO4--] / [HSO4-]
     pH-scale: free scale !!!

   =========================================================== */

   double iom0,laux1;

   iom0 = 19.924*s/(1000-1.005*s);   /* ionic strength */
   laux1 = -4276.1/tk + 141.328 -23.093*log(tk)
          +(-13856/tk + 324.57 - 47.986 * log(tk))*sqrt(iom0)
           +(35474/tk - 771.54 + 114.723 * log(tk))*iom0
             -2698/tk * sqrt(iom0)*iom0 + 1776/tk * iom0 * iom0
           + log(1-0.001005*s);
   return(exp(laux1));
}    /* --- end of K_S --- */



double TOTAL2FREE(double s, double tk)
{
/* ===========================================================
  convert from pH_Hansson ('total`) to pH_NBS ('free`):
  pH_Hansson = pH_NBS - log(1+S_T/K_S(s,tk))
  =========================================================== */

   double laux1,S_T;

   S_T = 0.14/96.062/1.80655*s;   /* [mol/kg soln] total sulfate
                                     Dickson and Goyet (1994) Ch.5 p.11 */
   laux1 = 1.0+S_T/K_S(s,tk);

   return(laux1);
}


double Kwater(double s, double tk)
{
/* ==============================================================

   Kwater = [H+] [OH-] : units: (mol/kg)**2
   ion product of water as a function of salinity (s [psu]) and
   absolute temperature (t [K])

	Millero (1995)(in Dickson and Goyet (1994, Chapter 5, p.18))
	$K_w$ in mol/kg-soln.
	pH-scale: pH$_{Hansson}$ ('total` scale).

   ============================================================== */

   aux1 = -13847.26 / tk + 148.96502 - 23.6521 * log(tk)
        + (118.67 / tk - 5.977 + 1.0495 * log(tk)) * sqrt(s)
        - 0.01615 * s;
#ifdef FREEPH
   aux1 -= log(TOTAL2FREE(s,tk));
#endif
   aux2  = exp(aux1);
   return(aux2);

}   /* --- end of Kwater --- */





double K_1(double s, double tk)
{
/* ===========================================================
   first acidity constant:
   [H^+] [HCO_3^-] / [H_2CO_3] = K_1

   (Roy et al., 1993 in Dickson and Goyet, 1994, Chapter 5, p. 14)
   pH-scale: 'total'

   =========================================================== */

   aux1 = 2.83655 - 2307.1266 / tk - 1.5529413 * log(tk)
        - (0.20760841 + 4.0484 / tk) * sqrt(s)
        + 0.0846834 * s - 0.00654208 * s * sqrt(s)
        + log(1 - 0.001005 * s);
#ifdef FREEPH
   aux1 -= log(TOTAL2FREE(s,tk));
#endif
   aux2 = 1.0e6*exp(aux1);       /* 10^6  <-> mol -> mumol */
   return(aux2);

}   /* --- end of K_1 --- */

double K_2(double s, double tk)
{
/* ===========================================================
   second acidity constant:
   [H^+] [CO_3^--] / [HCO_3^-] = K_2

   (Roy et al., 1993 in Dickson and Goyet, 1994, Chapter 5, p. 15)
   pH-scale: 'total'
   =========================================================== */


   aux1 =  -9.226508 - 3351.6106 / tk - 0.2005743 * log(tk)
        - (0.106901773 + 23.9722 / tk) * sqrt(s)
        + 0.1130822 * s - 0.00846934 * s * sqrt(s)
        + log(1 - 0.001005 * s);
#ifdef FREEPH
   aux1 -= log(TOTAL2FREE(s,tk));
#endif
   aux2 = 1.0e6*exp(aux1);       /* 10^6  <-> mol -> mumol */
   return(aux2);

}   /* --- end of K_2 --- */

#ifdef BORON
double fKbor(double s, double tk)
{
/* ===========================================================
   Boric acid constant:

   Kbor = [H+][B(OH)4-]/[B(OH)3]

   (Dickson, 1990 in Dickson and Goyet, 1994, Chapter 5, p. 14)
   pH-scale: 'total'


   =========================================================== */
   aux1  = (-8966.90-2890.53*sqrt(s)-77.942*s+1.728*pow(s,1.5)
            -0.0996*s*s);
   aux2  =  +148.0248+137.1942*sqrt(s)+1.62142*s;
   aux2 += (-24.4344-25.085*sqrt(s)-0.2474*s)*log(tk);

   aux3  = aux1 / tk + aux2 + 0.053105*sqrt(s)*tk;
#ifdef FREEPH
   aux3 -= log(TOTAL2FREE(s,tk));
#endif
   return(exp(aux3));	/* [mol/kg]  */

}   /* --- end of Kbor --- */
#endif

double fKh2co3(double tk)
{

/* ===========================================================

   Eigen and Hammes (1963)
     K_{H2CO3} = 2e-4 mol/l

   =========================================================== */
   return(2.e-4);
}   /* --- end of fKh2co3 --- */


/* =========================================================
   =========================================================

               ''chemical'' functions (end)

   =========================================================
   ========================================================= */



/* =========================================================
   =========================================================

               diffusion coefficients (begin)

   =========================================================
   ========================================================= */

double dicosw(double tk)
{
   double tc;
/* ============= diffusion coefficient correction for seawater ====
                 --                    --             -  -
   Li and Gregory, 1974, eq. (8) or Boudreau, 1994, p.52

   dicosw = 0.9508 - 0.0007379 * tc       (linear fit by dwg)

   =========================================================== */

   tc = tk - TKELVIN;
   return(0.9508 - 0.0007389 * tc);
}




double fdco2(double tk)
{
   double aco2,eaco2,laux1;
/* ============= CO2 diffusion coefficient ===================
   fresh water:
   D = A exp(-E_a/R T)
  $A = 5019\cdot 10^{-9}$ m$^2$/s, $E_a = 19.51$ kJ/mol activation energy;
  $R = 8.3143$ J/K/mol gas constant; $T$ absolute temperature.
  J"ahne et al. (1987)
   =========================================================== */

    aco2 = 5019.e-9 * 1.e12;   /* [mu**2/s] */
   eaco2 = 19.51e3; /* [J/mol] */
   laux1 = aco2 * exp(-eaco2/RGASJ/tk) * dicosw(tk);
   return(laux1);
}   /* --- end of fdco2 --- */

#ifdef C13ISTP
double fdcco2(double tk)
{
   double acco2,eacco2,laux1;
/* ============= 13CO2 diffusion coefficient ===================

   D(13CO2) = D(12CO2) * 0.9993, O'Leary (1984).

     =========================================================== */

    acco2 = 5019.e-9 * 1.e12;   /* [mu**2/s] */
   eacco2 = 19.51e3; /* [J/mol] */
   laux1 = acco2 * exp(-eacco2/RGASJ/tk) * dicosw(tk);
#define UD13EQD12
#ifdef D13EQD12
   return(laux1);
#else
   return(laux1*0.9993);
#endif
}   /* --- end of fdcco2 --- */
#endif

double fvisc(double tk)
{
   double d0,tc,mu0,tk0,laux1,laux2;
/* ============= dynamic viscosity of fresh water ========================
   Siedler and Peters (1986) :
     \log \frac{\mu (TC)}{\mu (20^\circ C)} =
          \frac{a (20-TC) - b (TC-20)^2}{TC+c}

   where $\mu (20^\circ C) = 1.002\cdot 10^{-3}$ N s m${-2}$,
   $a = 1.1709$, $b = 1.827 \cdot 10^{-3}$, $c = 89.93$ and TC is
   in $[^\circ C]$.
   =========================================================== */
   tc = tk - TKELVIN;
   mu0 = 1.002e-3;
   laux2 = (1.1709*(20.-tc) - 1.827e-3*(tc-20.)*(tc-20.))/(tc+89.93);
   laux1 = pow(10.,laux2)*mu0;
   return(laux1);
}   /* --- end of fvisc --- */

double fdh(double tk)
{
   double tk0,d0,laux1;
/* ============= H+ diffusion coefficient ===================
   D = k T [k] / f(T)  Stokes-Einstein relation [Atkins, 1990, p.765]
   D (25 C) = 9.31e-9 [m**2/s] [Li and Gregory, 1974]
   =========================================================== */
   tk0 = 25. + TKELVIN;
   d0 = 9.31e-9 * 1.e12;   /* [mu**2/s] */
   laux1 = d0 * tk/tk0 * fvisc(tk0)/fvisc(tk) * dicosw(tk);
   return(laux1);
}   /* --- end of fdh --- */


double fdoh(double tk)
{
   double tk0,d0,laux1;
/* ============= OH- diffusion coefficient ===================
   D = k T [k] / f(T)  Stokes-Einstein relation [Atkins, 1990, p.765]
   D (25 C) = 5.27e-9 [m**2/s] [Li and Gregory, 1974]
   =========================================================== */
   tk0 = 25. + TKELVIN;
   d0 = 5.27e-9 * 1.e12;   /* [mu**2/s] */
   laux1 = d0 * tk/tk0 * fvisc(tk0)/fvisc(tk) * dicosw(tk);
   return(laux1);
}   /* --- end of fdoh --- */

#ifdef BORON
double fdboh3(double tk)
{
   double tk0,dboh3,laux1;
/* ============= BOH3 diffusion coefficient ===================
   D = k T [k] / f(T)  Stokes-Einstein relation [Atkins, 1990, p.765]
   D (25 C) = 1.11e-9 [m**2/s] [Mackin, 1986] (seawater, ph 6-8)
   =========================================================== */
   tk0 = 25. + TKELVIN;
   dboh3 = 1.11e-9 * 1.e12;   /* [mu**2/s] */
   laux1 = dboh3 * tk/tk0 * fvisc(tk0)/fvisc(tk) ;
   return(laux1);
}   /* --- end of fdboh3 --- */

double fdboh4(double tk)
{
   double tk0,dboh4,laux1;
/* ============= BOH4 diffusion coefficient ===================
   D = k T [k] / f(T)  Stokes-Einstein relation [Atkins, 1990, p.765]
   D (25 C) = 0.97e-9 [m**2/s] [Boudreau and Canfield, 1993] (seawater)
   =========================================================== */
   tk0 = 25. + TKELVIN;
   dboh4 = 0.97e-9 * 1.e12;   /* [mu**2/s] */
   laux1 = dboh4 * tk/tk0 * fvisc(tk0)/fvisc(tk) ;
   return(laux1);
}   /* --- end of fdboh4 --- */
#endif

#ifdef OXYGEN

double fdo2(double tk)
{
   double tk0,do20,laux1;
/* ============= O2 diffusion coefficient ===================
   D = k T [k] / f(T)  Stokes-Einstein relation [Atkins, 1990, p.765]
   (D (21 C) = 2.33e-9 [m**2/s] [Lax, 1967, p. 1-1443])
   D (25 C) = 2.26e-9 [m**2/s] [Ramsing und Gunderson, 1994]
   =========================================================== */
   tk0 = 25. + TKELVIN;
   do20 = 2.26e-9 * 1.e12;   /* [mu**2/s] */
   /*laux1 = do20 * tk/tk0 * dicosw(tk);*/

   laux1 = do20 * tk/tk0 * fvisc(tk0)/fvisc(tk) ;

   return(laux1);
}   /* --- end of fdo2 --- */

#endif



#ifdef CALCIUM
double fdca(double tk)
{
   double tk0,d0,laux1;
/* ============= Ca++ diffusion coefficient ===================
   D = k T [k] / f(T)  Stokes-Einstein relation [Atkins, 1990, p.765]
   D (25 C) = 0.793e-9 [m**2/s] [Li and Gregory, 1974]
   =========================================================== */
   tk0 = 25. + TKELVIN;
   d0 = 0.793e-9 * 1.e12;   /* [mu**2/s] */
   laux1 = d0 * tk/tk0 * fvisc(tk0)/fvisc(tk) * dicosw(tk);
   return(laux1);
}   /* --- end of fdca --- */
#endif

double fdhco3(double tk)
{
   double tk0,d0,laux1;
/* ============= HCO3- diffusion coefficient ===================
   D = k T [k] / f(T)  Stokes-Einstein relation [Atkins, 1990, p.765]
   D (25 C) = 1.18e-9 [m**2/s] [Li and Gregory, 1974]
   =========================================================== */
   tk0 = 25. + TKELVIN;
   d0 = 1.18e-9 * 1.e12;   /* [mu**2/s] */
   laux1 = d0 * tk/tk0 * fvisc(tk0)/fvisc(tk) * dicosw(tk);

#define UDTEST
#ifdef DTEST
   printf("\n------------ DIFFC. TEST --------------\n");
   printf("%e %e %e   d0 tk tk0 \n",d0,tk,tk0);
   printf("%e %e   fvisc(tk) fvisc(tk0) \n",fvisc(tk),fvisc(tk0));
   printf("%e %e   dicosw(tk) dhco3\n\n",dicosw(tk),laux1);

#endif

   return(laux1);
}   /* --- end of fdhco3 --- */

#ifdef C13ISTP
double fdhcco3(double tk)
{
   double tk0,d0,laux1;
/* ============= H13CO3- diffusion coefficient ===================
   Assumed the same as for H12CO3-
   =========================================================== */
   tk0 = 25. + TKELVIN;
   d0 = 1.18e-9 * 1.e12;   /* [mu**2/s] */
   laux1 = d0 * tk/tk0 * fvisc(tk0)/fvisc(tk) * dicosw(tk);
   return(laux1);
}   /* --- end of fdhcco3 --- */
#endif

double fdco3(double tk)
{
   double tk0,d0,laux1;
/* ============= CO3-- diffusion coefficient ===================
   D = k T [k] / f(T)  Stokes-Einstein relation [Atkins, 1990, p.765]
   D (25 C) = 0.955e-9 [m**2/s] [Li and Gregory, 1974]
   =========================================================== */
   tk0 = 25. + TKELVIN;
   d0 = 0.955e-9 * 1.e12;   /* [mu**2/s] */
   fprintf(fppara,"%e   fvisc(tk0) \n",fvisc(tk0));
   fprintf(fppara,"%e   fvisc(tk)  \n",fvisc(tk));
   laux1 = d0 * tk/tk0 * fvisc(tk0)/fvisc(tk) * dicosw(tk);
   return(laux1);
}   /* --- end of fdco3 --- */

#ifdef C13ISTP
double fdcco3(double tk)
{
   double tk0,d0,laux1;
/* ============= 13CO3-- diffusion coefficient ===================
   Assumed the same as for 12CO3-
   =========================================================== */
   tk0 = 25. + TKELVIN;
   d0 = 0.955e-9 * 1.e12;   /* [mu**2/s] */
   laux1 = d0 * tk/tk0 * fvisc(tk0)/fvisc(tk) * dicosw(tk);
   return(laux1);
}   /* --- end of fdcco3 --- */
#endif

/* =========================================================
   =========================================================

               diffusion coefficients (end)

   =========================================================
   ========================================================= */




/* =========================================================
   =========================================================

               init (begin)

   =========================================================
   ========================================================= */



#ifdef C13ISTP
void initeps()
{
   /* ==== init eps-values for fractionation between carbon species === */

   /* -- Determine fractonation between C species.    ---*/
   /* -- Using formula given by Mook (1986).          ---*/

   /*--- 1. eps(CO2(g) - HCO3-)    --*/

   eps1 = -9483./TK + 23.89;

   /*--- 2. eps(CO2(aq) - CO2(g))  --*/

   eps2 = -373./TK + 0.19;

   /*--- 3. eps(CO2(aq) - HCO3-)   --*/

   eps3 = -9866./TK + 24.12;

  /*--- 4. eps(CO3-- - HCO3-)     --*/

   eps4 = -867./TK + 2.52;

  /*--- 5. eps(CaCO3 - HCO3-)     --*/

   eps5 = -4232./TK + 15.10;

  /*--- 6. eps(CaCO3 - CO3--), using 4 and 5 --*/

   eps6 = -3341./TK + 12.54;
}
#endif

void initc()
{
/* =============== init concentrations ======================= */

   double a0,a1,a2,dhds;
   /* ------   pH -> [H^+] and [OH^-]   ----- */

   kw = Kwater(SALINITY,TK);
   phbulk = PHBULK;
#ifdef CBNS
   phbulk = phbulkarg;   /* pHbulk is arg of main */
#endif
#ifdef CBNSLOOP
 #ifndef READ
   phbulk = phbulkv;
   printf("\n\n                   CBNSLOOP                 \n");
   printf("CBNSLOOP        %d ICLMAX %d icl       CBNSLOOP\n\n",ICLMAX,icl);
 #endif
 #ifdef READ
   phbulk = phin[icl-1];
   printf("\n\n                   CBNSLOOP                 \n");
   printf("CBNSLOOP        %d NDAT   %d icl       CBNSLOOP\n\n",NDAT,icl);
 #endif
#endif
   hbulk  = pow(10.0,-phbulk);   /* NBS-pH */
   ohbulk = kw / hbulk;
   /* ---- [mol/kg]  ->  [mumol/kg] ----- */
   hbulk  = hbulk  * 1.0e6;
   ohbulk = ohbulk * 1.0e6;
   kw = hbulk * ohbulk;

   printf("%f  salinity [psu] (input)    \n",SALINITY);
   printf("%f  temperature [C] (input)   \n",TEMP);
   printf("%f  pH-bulk (input)           \n",phbulk);
#ifdef PRINT
   printf("%e  K_water    [(mumol/kg)**2] \n",kw);
   printf("%e  [H+]-bulk  [mumol/kg]      \n",hbulk);
   printf("%e  [OH-]-bulk [mumol/kg]      \n",ohbulk);
#endif /* PRINT */

#ifdef CALCIUM
   cabulk = CABULK;
#ifdef PRINT
   printf("%f  Ca++ bulk (input) \n",cabulk); */
#endif /* PRINT */
#endif
#ifdef OXYGEN
   o2bulk = O2BULK;
#ifdef PRINT
   printf("%f  O2-bulk (input)           \n",o2bulk);
#endif /* PRINT */
#endif

#ifdef BORON
  kbdum = fKbor(SALINITY,TK)*1.e6; /* [mol/kg] -> [mumol/kg] */
  boh4bulk = BORTBULK / (1. + hbulk/kbdum);
  boh3bulk = BORTBULK - boh4bulk;

#ifdef BORISTP

  d11boh4bulk = ( D11BSEA*BORTBULK-EPSB*boh3bulk )
  		/ ( (1.+EPSB/1000.)*boh3bulk + boh4bulk );
  d11boh3bulk = d11boh4bulk*(1.+EPSB/1000.) + EPSB;


  btmp 		= (d11boh3bulk/1000. + 1.)*BSTAND;
  bboh3bulk 	= boh3bulk*btmp/(1.+btmp);
  btmp 		= (d11boh4bulk/1000. + 1.)*BSTAND;
  bboh4bulk 	= boh4bulk*btmp/(1.+btmp);

#ifndef B10B11
  printf("\n%e d11b3\n",d11boh3bulk);
  printf("%e d11b4\n",d11boh4bulk);
  printf("%e 10b3\n",boh3bulk-bboh3bulk);
  printf("%e 11b3\n",bboh3bulk);
  printf("%e b3\n",boh3bulk);
  printf("%e 10b4\n",boh4bulk-bboh4bulk);
  printf("%e 11b4\n",bboh4bulk);
  printf("%e b4\n\n",boh4bulk);

  printf("\n%e epsb\n\n",( (bboh3bulk/(boh3bulk-bboh3bulk))/(bboh4bulk/(boh4bulk-bboh4bulk)) - 1. )*1000.);
#endif

#ifdef B10B11

  /* boh3 and boh4 is 10B(OH)3 and 10B(OH)4-	*/

  boh3bulk = boh3bulk-bboh3bulk;
  boh4bulk = boh4bulk-bboh4bulk;

  printf("\n%e d11b3\n",d11boh3bulk);
  printf("%e d11b4\n",d11boh4bulk);
  printf("%e 10b3\n",boh3bulk);
  printf("%e 11b3\n",bboh3bulk);
  printf("%e b3\n",boh3bulk+bboh3bulk);
  printf("%e 10b4\n",boh4bulk);
  printf("%e 11b4\n",bboh4bulk);
  printf("%e b4\n\n",boh4bulk+bboh4bulk);

  printf("\n%e epsb\n\n",( (bboh3bulk/(boh3bulk))/(bboh4bulk/(boh4bulk)) - 1. )*1000.);
#endif


#endif
#endif

   /* ------   ALK, CO2, pH -> HCO3^-, CO3^--   --------- */

  k1d = K_1(SALINITY,TK);
  k2d = K_2(SALINITY,TK);
  a0 = 1.0 / (1.0 + k1d / hbulk + k1d * k2d / hbulk / hbulk);
  a1 = 1.0 / (hbulk / k1d + 1.0 + k2d / hbulk);
  a2 = 1.0 / (hbulk * hbulk / k1d / k2d + hbulk / k2d + 1.0);


#ifdef DICBULK

  dicbulk  = DICBULK;
#ifdef CBNS
  dicbulk  = dicbulkarg;	/* DIC is arg of main	*/
#endif
#ifdef CBNSLOOP
 #ifndef READ
  dicbulk  = DICBULK;		/* leave DIC const.	*/
 #endif
 #ifdef READ
  dicbulk  = dicin[icl-1];
 #endif
#endif
  co2bulk  = a0 * dicbulk;
  hco3bulk = a1 * dicbulk;
  co3bulk  = a2 * dicbulk;

#endif /* DICBULK */

#ifdef ALKBULK

  alkbulk  = ALKBULK;
#ifdef CBNS
  alkbulk  = alkbulkarg;	/* ALK is arg of main	*/
#endif
#ifdef CBNSLOOP
 #ifndef READ
  alkbulk  = ALKBULK;		/* leave  ALK const.	*/
 #endif
 #ifdef READ
  alkbulk  = alkin[icl-1];
 #endif
#endif

#ifdef BORON
  co2bulk = (alkbulk-kw/hbulk+hbulk-kbdum*BORTBULK/(kbdum+hbulk)) /
            (k1d/hbulk+2.*k1d*k2d/hbulk/hbulk);
#else
  co2bulk = (alkbulk-kw/hbulk+hbulk) /
            (k1d/hbulk+2.*k1d*k2d/hbulk/hbulk);
#endif
  dicbulk  = co2bulk / a0;
  hco3bulk = a1 * dicbulk;
  co3bulk  = a2 * dicbulk;

#endif /* ALKBULK */


#ifdef C13ISTP
   /*--- calculate d13C values of carbon species ---------*/
   /*    co2bulk is 12CO2+13CO2				  */
   /*    From mass balance for bulk concentrations:	  */
   /*	 d13C_DIC * [DIC] =       d13C_CO2  *  [CO2]
				+ d13C_HCO3 * [HCO3]
				+ d13C_CO3  *  [CO3]	  */
   /*    and
    	 d13C_CO2 = d13C_HCO3*(1+eps3/1000)   + eps3
    	 d13C_CO3 = d13C_HCO3*(1+eps4/1000)   + eps4

   	 follows:					  */


   d13hco3bulk  = DC13DIC*dicbulk - (eps3*co2bulk + eps4*co3bulk);
   d13hco3bulk /= ((1.+eps3/1000.)*co2bulk
		 + 		   hco3bulk
		 + (1.+eps4/1000.)*co3bulk);
   d13co2bulk   = d13hco3bulk*(1.+eps3/1000.) + eps3;
   d13co3bulk   = d13hco3bulk*(1.+eps4/1000.) + eps4;





  /*-----------  calculate [13CO2]_bulk from: ---------------*/
  /*
  d13CO2_bulk = (([13CO2_b.]/[12CO2_b.] / [13C_st.]/[12C_st.])-1)*1,000

  where

  [13C_st.]/[12C_st.] = R_st = 0.01124, O'Leary (1980).

  Note that [12CO2_b.] = [CO2_b.] - [13CO2_b.]		---*/


  cco2bulk   = co2bulk*RSTAND*(1+d13co2bulk/1000.);
  cco2bulk  /= (1 + RSTAND*(1+d13co2bulk/1000.));
  hcco3bulk  = hco3bulk*RSTAND*(1+d13hco3bulk/1000.);
  hcco3bulk /= (1 + RSTAND*(1+d13hco3bulk/1000.));
  cco3bulk   = co3bulk*RSTAND*(1+d13co3bulk/1000.);
  cco3bulk  /= (1 + RSTAND*(1+d13co3bulk/1000.));


   tmp  = d13hco3bulk*hco3bulk+d13co2bulk*co2bulk+d13co3bulk*co3bulk;
   tmp /= dicbulk;

#ifdef CISTP	/* NOW: BULK CONCENTRATIONS ARE SET TO 12C CONC. */

   co2bulk  = co2bulk  -  cco2bulk;
   hco3bulk = hco3bulk - hcco3bulk;
   co3bulk  = co3bulk  -  cco3bulk;

   tmp  = d13hco3bulk*(hco3bulk+hcco3bulk)
   	 +d13co2bulk*(co2bulk+cco2bulk)
   	 +d13co3bulk*(co3bulk+cco3bulk);
   tmp /= dicbulk;

   /*printf("%e\n",(cco2bulk/co2bulk/RSTAND-1.)*1000.);
   printf("%e\n",(hcco3bulk/hco3bulk/RSTAND-1.)*1000.);
   printf("%e\n",(cco3bulk/co3bulk/RSTAND-1.)*1000.);*/
#endif

#ifdef PRINT
   printf("\n ========== D13C values  =========== \n");
   printf("DC13DIC     (0/oo) %f \n",tmp);
   printf("D13CO2bulk  (0/oo) %f \n",d13co2bulk);
   printf("D13HCO3bulk (0/oo) %f \n",d13hco3bulk);
   printf("D13CO3bulk  (0/oo) %f \n\n",d13co3bulk);
   printf("D13CRES     (0/oo) %f \n",D13CRES);
   printf("D13CPHS     (0/oo) %f \n",D13CPHS);
#endif


#endif

  h2obulk = MOLH2O;

   /* ------ calculate d[H+]/d[CO2] -------------
        compare eq. (11) Wolf-Gladrow (1992)        */

   dhds = (k1d/hbulk + 2.*k1d*k2d/hbulk)
        / (1.+kw/hbulk/hbulk
          +co2bulk*(k1d/hbulk/hbulk+4.*k1d*k2d/hbulk/hbulk/hbulk));
#ifdef PRINT
   printf("%f  dhds = d[H+]/d[CO2]       \n",dhds);
   printf("%f  CO2-atm [ppmv] (input)    \n",CO2ATM);
   printf("%f  CO2-bulk [mumol/kg]        \n",co2bulk);
   printf("%f  h2obulk  [mumol/kg]        \n",h2obulk);
   printf("%e  k1d first acidity constant [mumol/kg]   \n",k1d);
   printf("%e  k2d second acidity constant [mumol/kg]  \n",k2d);
#ifdef BORON
   printf("%e  boric acid constant [mumol/kg] \n",kbdum);
#endif
   printf("%f  alpha-0 * 100 (CO2)   \n",(a0 * 100.0));
   printf("%f  alpha-1 * 100 (HCO3-) \n",(a1 * 100.0));
   printf("%f  alpha-2 * 100 (CO3--) \n",(a2 * 100.0));
#endif /* PRINT */
#ifdef ALKBULK
   printf("%f   alkbulk    [mueq /kg] \n",alkbulk);
#endif
   printf("%f   dicbulk    [mumol/kg] \n",dicbulk);
   printf("%f   co2bulk    [mumol/kg] \n",co2bulk);
   printf("%f  hco3bulk    [mumol/kg] \n",hco3bulk);
   printf("%f   co3bulk    [mumol/kg] \n",co3bulk);
#ifdef BORON
   printf("%f  total boron [mumol/kg] \n",BORTBULK);
   printf("%f  boh3bulk    [mumol/kg] \n",boh3bulk);
   printf("%f  boh4bulk    [mumol/kg] \n",boh4bulk);
#endif
#ifdef PRINT
#ifdef C13ISTP
   printf("%f 13co2bulk    [mumol/kg] \n",cco2bulk);
   printf("%f 13hco3bulk   [mumol/kg] \n",hcco3bulk);
   printf("%f 13co3bulk    [mumol/kg] \n",cco3bulk);
#endif
#endif /* PRINT */

   fprintf(fppara,"=====================");
   fprintf(fppara,"============= initc start =============== \n");
   fprintf(fppara,"salinity    [psu] (input)       %f \n",SALINITY);
   fprintf(fppara,"temperature [C] (input)         %f \n",TEMP);
   fprintf(fppara,"K_water     [(mumol/kg)**2]      %e \n",kw);
   fprintf(fppara,"CO2-atm     [ppmv] (input)      %f \n",CO2ATM);
   fprintf(fppara,"first acidity constant  [mumol/kg]  %e \n",k1d);
   fprintf(fppara,"second acidity constant [mumol/kg]  %e \n",k2d);
#ifdef BORON
   fprintf(fppara,"boric acid constant     [mumol/kg]  %e \n",kbdum);
#endif
   fprintf(fppara,"alpha-0 * 100 (CO2)       %f \n",(a0 * 100.0));
   fprintf(fppara,"alpha-1 * 100 (HCO3-)     %f \n",(a1 * 100.0));
   fprintf(fppara,"alpha-2 * 100 (CO3--)     %f \n",(a2 * 100.0));
   fprintf(fppara,"============= bulk concentration =============== \n");
   fprintf(fppara,"pH-bulk (input)            %f \n",phbulk);
   fprintf(fppara,"[H+]-bulk  [mumol/kg]       %e \n",hbulk);
   fprintf(fppara,"[OH-]-bulk [mumol/kg]       %e \n",ohbulk);
#ifdef CALCIUM
   fprintf(fppara,"Ca++ bulk (input)    %f \n",CABULK);
#endif
   fprintf(fppara,"CO2-bulk    [mumol/kg]   %f \n",co2bulk);
   fprintf(fppara,"DIC-bulk    [mumol/kg]   %f \n",dicbulk);
   fprintf(fppara,"HCO3- bulk  [mumol/kg]   %f \n",hco3bulk);
   fprintf(fppara,"CO3-- bulk  [mumol/kg]   %f \n",co3bulk);
   fprintf(fppara,"H2O bulk    [mumol/kg]   %f \n",h2obulk);
#ifdef C13ISTP
   fprintf(fppara,"eps1        (0/oo) %f \n",eps1);
   fprintf(fppara,"eps2        (0/oo) %f \n",eps2);
   fprintf(fppara,"eps3        (0/oo) %f \n",eps3);
   fprintf(fppara,"eps4        (0/oo) %f \n",eps4);
   fprintf(fppara,"eps5        (0/oo) %f \n",eps5);
   fprintf(fppara,"eps6        (0/oo) %f \n",eps6);
   fprintf(fppara,"DC13DIC     (0/oo) %f \n",DC13DIC);
   fprintf(fppara,"D13CO2bulk  (0/oo) %f \n",d13co2bulk);
   fprintf(fppara,"D13HCO3bulk (0/oo) %f \n",d13hco3bulk);
   fprintf(fppara,"D13CO3bulk  (0/oo) %f \n",d13co3bulk);
   fprintf(fppara,"13CO2-bulk  [mumol/kg]   %f \n",cco2bulk);
   fprintf(fppara,"13HCO3-bulk [mumol/kg]   %f \n",hcco3bulk);
   fprintf(fppara,"13CO3-bulk  [mumol/kg]   %f \n",cco3bulk);
#endif
#ifdef BORON
   fprintf(fppara,"total boron [mumol/kg]   %f \n",BORTBULK);
   fprintf(fppara,"boh3bulk    [mumol/kg]   %f \n",boh3bulk);
   fprintf(fppara,"boh4bulk    [mumol/kg]   %f \n",boh4bulk);
#endif
   fprintf(fppara,"dhds = d[H+]/d[CO2]    %f \n",dhds);
   alkbulkd = hco3bulk + 2. * co3bulk;
#ifdef CISTP
   alkbulkd = hco3bulk+hcco3bulk + 2. * (co3bulk+cco3bulk);
#endif
   fprintf(fppara,"alkbulk1 (CA)       [mueq/kg]       %f \n",alkbulkd);
   alkbulkd = alkbulkd + ohbulk - hbulk;
   fprintf(fppara,"alkbulk2 (CA+WA)    [mueq/kg]       %f \n",alkbulkd);
#ifdef BORON
   alkbulkd = alkbulkd + boh4bulk;
   fprintf(fppara,"alkbulk3 (CA+WA+BA) [mueq/kg]       %f \n",alkbulkd);
#endif
   fprintf(fppara,"=====================");
   fprintf(fppara,"============= initc end =============== \n");

}   /* ----- end of initc ----- */


void initk()
{

/* ============= init kinetic and diffusion coefficients ============= */

/* ------------ diffusion coefficients ---------- */

#ifdef AGG
#ifdef OXYGEN
   do2   = DIFFAGG*fdo2(TK);
#endif
#ifdef BORON
   dboh3 = DIFFAGG*fdboh3(TK);
   dboh4 = DIFFAGG*fdboh4(TK);
#endif
   dco2  = DIFFAGG*fdco2(TK);
   dhco3 = DIFFAGG*fdhco3(TK);
   dco3  = DIFFAGG*fdco3(TK);
   dh    = DIFFAGG*fdh(TK);
   doh   = DIFFAGG*fdoh(TK);

#else
#ifdef OXYGEN
   do2   = fdo2(TK);
#endif

#ifdef C13ISTP
   dcco2  = fdcco2(TK);
   dhcco3 = fdhcco3(TK);
   dcco3  = fdcco3(TK);
#endif

#ifdef BORON
   dboh3 = fdboh3(TK);
   dboh4 = fdboh4(TK);
#ifdef BORISTP
   dbboh3 = fdboh3(TK);
   dbboh4 = fdboh4(TK);
#endif
#endif

   dco2  = fdco2(TK);
   dhco3 = fdhco3(TK);
   dco3  = fdco3(TK);
   dh    = fdh(TK);
   doh   = fdoh(TK);
#ifdef CALCIUM
   dca   = fdca(TK);
   data  = dca;
#endif

#endif /* AGG */

#ifdef REDDH
   dh    = dco2;
   fprintf(fppara,"diffusion coeff. of H+ reduced \n");
#endif

#ifdef BORON
   Kbor = fKbor(SALINITY,TK);
#endif



   kh2co3 = fkh2co3(TK);

/*       Kh2co3 = 2.0e-4 * 1.e6 ;
        [mumol/kg] Stumm and Morgan [1981, p. 210]     */

   Kh2co3 = fKh2co3(TK) * 1.e6;   /* e6: mol/kg -> mumol/kg */

   kp1s = fkp1s(SALINITY,TK);

   km1s = kh2co3/Kh2co3;

   aux1 = kp1s/km1s;

   fprintf(fppara,"------------------------- \n");

   fprintf(fppara,"kp1s                   %e \n",kp1s);

   fprintf(fppara,"km1s (before consist)  %e \n",km1s);

   fprintf(fppara,"kp1s/km1s              %e \n",aux1);

   aux2 = K_1(SALINITY,TK);

   fprintf(fppara,"K1                     %e \n",aux2);

   km1s = kp1s / K_1(SALINITY,TK);   /* consist */

   fprintf(fppara,"km1s (after consist)   %e \n",km1s);

   fprintf(fppara,"------------------------- \n");

   kp4 = fkp4(TK);


   /* before 5/96
   km4 = kp4 * co2bulk * ohbulk / hco3bulk;     	     consist */

   /* after 5/96 */
   km4 = kp4 * Kwater(SALINITY,TK) * km1s / kp1s * 1.e12; /* consist */

   fprintf(fppara,"kp4(T)                    %e \n",kp4);

   fprintf(fppara,"km4(T) (consist value)    %e \n",km4);

   fprintf(fppara,"------------------------- \n");


#define CONKH2CO3

#ifdef CONKH2CO3

   kh2co3 = ((kp1s + kp4 * ohbulk) * co2bulk / hco3bulk - km4)
            * Kh2co3 / hbulk;

   fprintf(fppara," --- make kh2co3 consistent --- \n");
   fprintf(fppara,"kh2co3          [1/s]         %f \n",kh2co3);
   fprintf(fppara,"------------------------- \n");

#endif


/* ---          k21  k12
        H2CO3   ->   <-   H+ + HCO3-
        k21 = 9.e6 [1/s]
        k21/k12 = Kh2co3                 --- */

   k21 = fk21(TK);    /* [1/s] Eigen and Hammes, 1963 */
   k12 = k21/Kh2co3;


   k32 = fk32(TK);

   k13 = fk13(TK);

   k31 = fk31(TK);                     /* [1/s] */
   k13 = k31 / K_1(SALINITY,TK);    /* [kg/mumol/s] */

/* ---          km5h  kp5h
        HCO3-   ->   <-    H+ + CO3--
   values not known:
   km5h/kp5h = K_2                               kp1s
   assume: reaction is fast compared to CO2 + H2O  ->  H+ + HCO3- or
                                                       H2CO3  --- */

   kp5h = 1.e-6 * fkp5h(TK);              /* [kg/mumol/s] */
   km5h = K_2(SALINITY,TK) * kp5h;     /* [1/s] */

/* ---      kp6  km6
       H2O  ->   <-  H+ + OH-
       kp6/km6 = Kwater          --- */

   kp6 = fkp6(TK);  /* [mumol/kg/s] */

   km6 = kp6/Kwater(SALINITY,TK)/1.e12;
       /* -> K_W o.k.; e12 -> mol^2 -> mumol^2 */
       /* km6: [kg/mumol/s] */


#ifdef BORONRC4

   kp7 = fkp7(TK)*1.e-6;	   /* [kg/mol/s] -> [kg/mumol/s] */
   km7 = kp7*Kwater(SALINITY,TK)/fKbor(SALINITY,TK)*1.e6; /* [1/s] */

#ifdef BORISTP

#ifndef B10B11

   /*alphab = (1.+EPSB/1000.)*(boh3bulk-bboh3bulk)*boh4bulk
   	    / ( (boh4bulk-bboh4bulk)*boh3bulk );*/

   alphab  = (1.+EPSB/1000.);
   alphabp = alphab*(1.-bboh3bulk/boh3bulk)/(1.-bboh4bulk/boh4bulk);
   printf("%e alpha(BOR)\n",alphab);
   printf("%e alpha'(BOR)\n",alphabp);


   /* NOTE: kp7 is for total B (80% 11B, 20% 10B) 	*/
   /*       kp7bb is for 11B, so kp7 is 		*/
   /*	    a LITTLE bit higher than kp7bb	 	*/

   km7bb = km7/1.030;
   kp7bb = kp7 * km7bb / km7 / alphabp;

#endif

#ifdef B10B11

   alphab  = (1.+EPSB/1000.);
   alphabp = alphab*(1.-bboh3bulk/(boh3bulk+bboh3bulk))
            /(1.-bboh4bulk/(boh4bulk+bboh4bulk));
   printf("%e alpha(BOR)\n",alphab);
   printf("%e alpha'(BOR)\n",alphabp);

   /* NOTE: kp7,km7 is for total B (80% 11B, 20% 10B) 	*/
   /*       Kbor   = kp7/km7 * Kwater		*/
   /*	    11Kbor = Kbor/alphabp	 	*/
   /*	    10Kbor = Kbor*alphab/alphabp	*/

   K11bor = Kwater(SALINITY,TK)*kp7/km7/alphabp;
   K10bor = Kwater(SALINITY,TK)*kp7*alphab/km7/alphabp;

   /* kp7 and km7 are now used as 10kp7 and 10km7! 	*/
   /* let km7 unchanged, set kp7 = km7*K10bor		*/

   kp7   = km7*K10bor/Kwater(SALINITY,TK);
   km7bb = km7/1.030;
   kp7bb = kp7 * km7bb / km7 / alphab;

#endif
   printf("\n%e km7\n",km7);
   printf("%e km7bb\n",km7bb);
   printf("%e kp7\n",kp7);
   printf("%e kp7bb\n",kp7bb);

   printf("%e kp7/kp7bb\n",kp7/kp7bb);
   printf("%e km7/km7bb\n",km7/km7bb);

#endif
#endif



#ifdef BORONRC3

   km7 = fkm7(TK) * 1.e-6;	   /* [kg/mol/s] -> [kg/mumol/s] */
   kp7 = km7 * fKbor(SALINITY,TK) * 1.e6; /* [1/s] */

#ifdef BORISTP

#ifndef B10B11

   /*alphab = (1.+EPSB/1000.)*(boh3bulk-bboh3bulk)*boh4bulk
   	    / ( (boh4bulk-bboh4bulk)*boh3bulk );*/

   alphab  = (1.+EPSB/1000.);
   alphabp = alphab*(1.-bboh3bulk/boh3bulk)/(1.-bboh4bulk/boh4bulk);
   printf("%e alpha(BOR)\n",alphab);
   printf("%e alpha'(BOR)\n",alphabp);


   /* NOTE: kp7 is for total B (80% 11B, 20% 10B) 	*/
   /*       kp7bb is for 11B, so kp7 is 		*/
   /*	    a LITTLE bit higher than kp7bb	 	*/

   km7bb = km7/1.030;
   kp7bb = kp7 * km7bb / km7 / alphabp;

#endif

#ifdef B10B11

   alphab  = (1.+EPSB/1000.);
   alphabp = alphab*(1.-bboh3bulk/(boh3bulk+bboh3bulk))
            /(1.-bboh4bulk/(boh4bulk+bboh4bulk));
   printf("%e alpha(BOR)\n",alphab);
   printf("%e alpha'(BOR)\n",alphabp);

   /* NOTE: kp7,km7 is for total B (80% 11B, 20% 10B) 	*/
   /*       Kbor   = kp7/km7 			*/
   /*	    11Kbor = Kbor/alphabp	 	*/
   /*	    10Kbor = Kbor*alphab/alphabp	*/

   K11bor = kp7/km7/alphabp;
   K10bor = kp7*alphab/km7/alphabp;

   /* kp7 and km7 are now used as 10kp7 and 10km7! 	*/
   /* let km7 unchanged, set kp7 = km7*K10bor		*/

   kp7   = km7*K10bor;
   km7bb = km7/1.030;
   kp7bb = kp7 * km7bb / km7 / alphab;

#endif
   printf("\n%e km7\n",km7);
   printf("%e km7bb\n",km7bb);
   printf("%e kp7\n",kp7);
   printf("%e kp7bb\n",kp7bb);

   printf("%e kp7/kp7bb\n",kp7/kp7bb);
   printf("%e km7/km7bb\n",km7/km7bb);

#endif

   printf("%e %e K11bor K10bor\n",K11bor,K10bor);

#endif

#ifdef C13ISTP
   /* Calculate reaction rates for C13.				*/
   /* Determine fractionation factors from eps-values by Mook   */
   /*

	alpha   = [13CO2]/[12CO2] / [H13CO3]/[H12CO3]

	alpha5  = [13CO3]/[12CO3] / [H13CO3]/[H12CO3]

        alphac  = [13CaCO3]/[12CaCO3] / [13CO3]/[12CO3]

        alphahc = [13CaCO3]/[12CaCO3] / [13HCO3]/[12HCO3]  	*/


   alpha   =  eps3/1000. + 1.;	/* eps(CO2(aq) - HCO3-) 	*/
   alpha5  =  eps4/1000. + 1.;	/* eps(CO3--   - HCO3-)		*/
   alphac  =  eps6/1000. + 1.;	/* eps(CaCO3   - CO3--)		*/
   alphahc =  eps5/1000. + 1.;	/* eps(CaCO3   - HCO3-)		*/

   /* Reactions 1: CO2 + H2O -> H+ + HCO3- and
	        4: CO2 + OH- -> HCO3-
      fractionate equally in eq. (equilibrium fractionation is
      independent of the way of the reaction).


      ==>    alpha = K1 / K1'  = kp1s/km1s / kp1scc/km1scc

	     alpha = K1 / K1'  = kp4 /km4  / kp4cc /km4cc 	*/

   /* Note that alpha is given for 13C/12C by Mook '86.		*/
   /* In this program, however, the reaction constants must be  */
   /* given for 13C/Cges. Thus, alpha is slightly modified.	*/
   /* From

	alpha  = 13co2/12co2 / h13co3/h12co3

	alpha' = 13co2/(12co2+13co2) / h13co3/(h12co3+h12co3)   */

   /* follows:							*/

#ifndef CISTP
   /* For 12C and 13C alpha's are not modified	*/
   alpha  *= (1.-cco2bulk/co2bulk)/(1.-hcco3bulk/hco3bulk);
   alpha5 *= (1.-cco3bulk/co3bulk)/(1.-hcco3bulk/hco3bulk);
#endif
#ifdef CISTP
   /* However, reaction rates must be slightly modified !!	*/
   /* k's are calc. for total C, not for 12C			*/
   /* Thus, leave k+'s unchanged and calc. 12k's from		*/
   /* 12C equilibrium (co2 == 12co2 ,...).			*/


   km1s = kp1s*co2bulk/hbulk/hco3bulk;
   km4  = kp4*ohbulk*co2bulk/hco3bulk;
   km5h = kp5h*hbulk*co3bulk/hco3bulk;

#endif

   /* And now .... the reaction rates for 13C !			*/


   /* Marlier and O'Leary (1984)
   kp1scc  =         kp1s / 1.0069;
   km1scc  = alpha * km1s / 1.0069; */

   /* O'Leary (pers. comm. 1996) k-1/k'-1 = 1.022 */
   kp1scc  =         kp1s / alpha / 1.022;
   km1scc  = 	     km1s / 1.022;


   /* Siegenthaler and Muennich, 1981
   kp4cc   =        .973*kp4;
   km4cc   =  alpha*.973*km4;	*/

   /* O'Leary (pers. comm. 1996) k-4/k'-4 = 1.011 */
   kp4cc  =          kp4 / alpha / 1.011;
   km4cc  = 	     km4 / 1.011;


   /* before 17.06.96 						*/
   /* Both, kpcc and kmcc are reduced and increased proportional*/
   /* to kp and km ( kpcc = kp * (1+d) and  kmcc = km * (1-d) ).*/
   /* ("Symmetric change of kp' and km' ")

   kp1scc  = (kp1s*2.)/(1.+alpha);
   km1scc  = (km1s*2.)*alpha/(1.+alpha);

   kp4cc   = (kp4*2.)/(1.+alpha);
   km4cc   = (km4*2.)*alpha/(1.+alpha);				 */

   kp5cc   = (kp5h*2.)/(1.+alpha5);
   km5cc   = (km5h*2.)*alpha5/(1.+alpha5);

   /* printf("%e kp5h/kp5cc\n",kp5h/kp5cc);
   printf("%e km5h/km5cc\n",km5h/km5cc);
   printf("%e K/K'\n",(kp5h/km5h)/(kp5cc/km5cc));

   kp5cc   =         .999*kp5h;
   km5cc   =  alpha5*.999*km5h;

   printf("%e kp5h/kp5cc\n",kp5h/kp5cc);
   printf("%e km5h/km5cc\n",km5h/km5cc);
   printf("%e K/K'\n",(kp5h/km5h)/(kp5cc/km5cc));*/

#endif

#ifdef PRINT

   printf("============= initk start =============== \n");
   printf("%e  kp1s   [1/s]        \n",kp1s);
   printf("%e  km1s   [kg/mumol/s] \n",km1s);
   printf("%e  kp4    [kg/mumol/s] \n",kp4);
   printf("%e  km4    [1/s]        \n",km4);
   printf("%e  kp5h   [kg/mumol/s] \n",kp5h);
   printf("%e  km5h   [1/s]        \n",km5h);
#ifdef C13ISTP
   printf("-------------- 13C ---------------------- \n");
   printf("%e  kp1scc [1/s]        \n",kp1scc);
   printf("%e  km1scc [kg/mumol/s] \n",km1scc);
   printf("%e  kp4cc  [kg/mumol/s] \n",kp4cc);
   printf("%e  km4cc  [1/s]        \n",km4cc);
   printf("%e  kp5cc  [kg/mumol/s] \n",kp5cc);
   printf("%e  km5cc  [1/s]        \n",km5cc);
   printf("----------------------------------------- \n");
#endif
#ifdef BORONRC4
   printf("%e  kp7 [kg/mumol/s] \n",kp7);
   printf("%e  km7 [1/s] \n",km7);
#endif
#ifdef BORONRC3
   printf("%e  kp7 [1/s] \n",kp7);
   printf("%e  km7 [kg/mumol/s] \n",km7);
#endif
   printf("%e  kh2co3 [1/s] \n",kh2co3);
   printf("%e  Kh2co3 [mumol/kg] \n",Kh2co3);

   printf("--- ratio of reaction (1) and (4) ----\n");
   printf("%e  kp1s      [CO2]_b \n",kp1s*co2bulk);
   printf("%e  km1s[H+][HCO3-]_b \n",km1s*hbulk*hco3bulk);
   printf("%e  kp4[OH-]  [CO2]_b \n",kp4*ohbulk*co2bulk);
   printf("%e  km4     [HCO3-]_b \n",km4*hco3bulk);
   printf("%12.3f  fraction (1) \n",kp1s*co2bulk/(kp1s*co2bulk+km4*hco3bulk));
   printf("%12.3f  fraction (4) \n",km4*hco3bulk/(kp1s*co2bulk+km4*hco3bulk));


   printf("============= initk end   =============== \n");

#endif /*PRINT*/

#ifdef MIMECO2SYM
   /* if defined SYMTCUPT < VMAX 		*/
   /*   set vmaxco2 = SYMTCUPT.			*/
   /* (else vmaxco2 = VMAX, see initializ. above)*/

   if(SYMTCUPT < VMAX)
		vmaxco2 = SYMTCUPT;
#ifdef LDAT
   /* calculate symtcup as a function of [CO3--].
      data from:
      symtcup(ph=8.15,co3=206) = 7.2 Spero et al., 91
      symtcup(ph=8.9 ,co3=590) = 2.8 Rink unp. data */

   			/* C uptake linear function of CO3 */
   tmp = (7.2-SYMTCUPPH89)/(590.-206);
   symtcup = ( ((7.2+tmp*206.) - tmp*co3bulk)*1.e-9/3600.);

#define SYMUPTF_CO2 	/* C uptake linear function of CO2 */
#ifdef  SYMUPTF_CO2
   tmp = (7.2-SYMTCUPPH89)/(13.4-1.2);
   symtcup = ( ((2.8-tmp*1.2) + tmp*co2bulk)*1.e-9/3600.);
   if(symtcup > VMAX)
   		symtcup = VMAX;

#endif

   if(symtcup <= VMAX)
		vmaxco2 = symtcup;
#endif /* LDAT */

#endif /* MIMECO2SYM */

#ifdef LDAT
   /* calculate co3uptldat as a function of [CO3--]. */
   if(co3bulk <= 200.){
   	tmp = (3.-1.)/(200.-50.);
   	co3uptldat = ( ((1.-tmp*50.) + tmp*co3bulk)*1.e-9/3600.);
   }
   if(co3bulk > 200.)
   	co3uptldat = 3.e-9/3600.;

   printf("%e co3uptldat\n",co3uptldat*3600.);
#endif


#ifdef C13ISTP
   /* Calculate C13 uptake from d13F_j.	(F_1 = F_resp.)		*/
   /*					(F_2 = F_phot.)		*/

   /* We have

      d13F_j = ( (13F_j/12F_j)/R_st.) - 1) * 1000

      from wich follows: (note that F_j = 13F_j + 12F_j)     */

   /* --------- 13CO2 UPTAKE AT THE SHELL (respiration) -----*/

   tmp      = RSTAND*(D13CRES/1000. + 1.);
   cco2upt  = CO2UPT*tmp/(1.+tmp);

#ifdef DIATOMW
   d13co2upt = d13co2bulk + EPSPCO2UPT;
   printf("d13co2upt %e \n",d13co2upt);
   tmp      = RSTAND*(d13co2upt/1000. + 1.);
   cco2upt  = CO2UPT*tmp/(1.+tmp);
#endif


#ifdef F13_CO3

   hcco3upt = 0.0;


   /* --------- 13CO3 UPTAKE AT THE SHELL (calcification) ---*/

   /* 	The d13C_shell is the final goal of the model (see
	difeq(). Here, an initial guess of the 13CO3--
	flux at the shell is given (equilibrium fractionation
	between HCO3- (bulk) and CaCO3).		     */


   tmp      = RSTAND*((d13hco3bulk+eps5)/1000. + 1.);
   cco3upt  = CO3UPT*tmp/(1.+tmp);	/* first guess !!!!  */
#ifdef LDAT
   cco3upt  = co3uptldat*tmp/(1.+tmp);	/* first guess !!!!  */
#endif


   /* --------- 13C UPTAKE OF SYMBIONTS (photosynth.)   -----*/

#ifdef SYMBIONTS
   tmp         = RSTAND*(D13CPHS/1000. + 1.);
   symcco2upt  = SYMCO2UPT*tmp/(1.+tmp);

   tmp         = RSTAND*(D13CPHS/1000. + 1.);
   symhcco3upt = SYMHCO3UPT*tmp/(1.+tmp);
#endif

#ifdef PRINT
   if(CO2UPT != 0.0){
   	tmp = (cco2upt/(CO2UPT-cco2upt)/RSTAND - 1.)*1000;
   	printf("%e d13F_res\n",tmp);
   }
   if(CO3UPT != 0.0){
   	tmp = (cco3upt/(CO3UPT-cco3upt)/RSTAND - 1.)*1000;
#ifdef LDAT
   	tmp = (cco3upt/(co3uptldat-cco3upt)/RSTAND - 1.)*1000;
#endif
   	printf("%e d13hco3bulk            \n",d13hco3bulk);
   	printf("%e eps(CaCO3 - HCO3-)     \n",eps5);
   	printf("%e d13F_calc (init. guess)\n",tmp);
   }
   if(SYMCO2UPT != 0.0){
   	tmp = (symcco2upt/(SYMCO2UPT-symcco2upt)/RSTAND - 1.)*1000;
   	printf("%e d13F_phs(CO2)\n",tmp);
   }
   if(SYMHCO3UPT != 0.0){
   	tmp = (symhcco3upt/(SYMHCO3UPT-symhcco3upt)/RSTAND - 1.)*1000;
   	printf("%e d13F_phs(HCO3)\n",tmp);
   }
#endif /* PRINT   */
#endif /* F13_CO3 */


#ifdef F13_HCO3

   cco3upt = 0.0;


   /* --------- H13CO3 UPTAKE AT THE SHELL (calcification) ---*/

   /* 	The d13C_shell is the final goal of the model (see
	difeq(). Here, an initial guess of the H13CO3-
	flux at the shell is given (equilibrium fractionation
	between HCO3- (bulk) and CaCO3).		     */


   tmp       = RSTAND*((d13hco3bulk+eps5)/1000. + 1.);
   hcco3upt  = HCO3UPT*tmp/(1.+tmp);	/* first guess !!!!  */
#ifdef LDAT
   hcco3upt  = hco3uptldat*tmp/(1.+tmp);	/* first guess !!!!  */
#endif


   /* --------- 13C UPTAKE OF SYMBIONTS (photosynth.)   -----*/

#ifdef SYMBIONTS
   tmp         = RSTAND*(D13CPHS/1000. + 1.);
   symcco2upt  = SYMCO2UPT*tmp/(1.+tmp);

   tmp         = RSTAND*(D13CPHS/1000. + 1.);
   symhcco3upt = SYMHCO3UPT*tmp/(1.+tmp);
#endif

#ifdef PRINT
   if(CO2UPT != 0.0){
   	tmp = (cco2upt/(CO2UPT-cco2upt)/RSTAND - 1.)*1000;
   	printf("%e d13F_res\n",tmp);
   }
   if(HCO3UPT != 0.0){
   	tmp = (hcco3upt/(HCO3UPT-hcco3upt)/RSTAND - 1.)*1000;
#ifdef LDAT
   	tmp = (hcco3upt/(hco3uptldat-hcco3upt)/RSTAND - 1.)*1000;
#endif
   	printf("%e d13hco3bulk            \n",d13hco3bulk);
   	printf("%e eps(CaCO3 - HCO3-)     \n",eps5);
   	printf("%e d13F_calc (init. guess)\n",tmp);
   }
   if(SYMCO2UPT != 0.0){
   	tmp = (symcco2upt/(SYMCO2UPT-symcco2upt)/RSTAND - 1.)*1000;
   	printf("%e d13F_phs(CO2)\n",tmp);
   }
   if(SYMHCO3UPT != 0.0){
   	tmp = (symhcco3upt/(SYMHCO3UPT-symhcco3upt)/RSTAND - 1.)*1000;
   	printf("%e d13F_phs(HCO3)\n",tmp);
   }
#endif /* PRINT   */
#endif /* F13_HCO3 */
#endif /* C13ISTP */


#ifdef PRINT

#ifdef FORAMW
   printf("==== WORK FORAM: UPTAKE (INPUT)   ======= \n");
   printf("%e  SYMCO2UPT  [mol/h] \n",SYMCO2UPT*3600.);
   printf("%e  SYMHCO3UPT [mol/h] \n",SYMHCO3UPT*3600.);
   printf("%e  SYMHUPT    [mol/h] \n",SYMHUPT*3600.);
   printf("%e  CO3UPT     [mol/h] \n",CO3UPT*3600.);
   printf("%e  CO2UPT     [mol/h] \n",CO2UPT*3600.);
#ifdef C13ISTP
   printf("-------------- 13C ---------------------- \n");
   printf("%e  SYMCCO2UPT  [mol/h] \n",symcco2upt*3600.);
   printf("%e  SYMHCCO3UPT [mol/h] \n",symhcco3upt*3600.);
   printf("%e  SYMHUPT     [mol/h] \n",SYMHUPT*3600.);
   printf("%e  CCO3UPT     [mol/h] \n",cco3upt*3600.);
   printf("%e  CCO2UPT     [mol/h] \n",cco2upt*3600.);
#endif
#ifdef OXYGEN
   printf("%e  O2UPT      [mol/h] \n",-REDF*CO2UPT*3600.);
   printf("%e  REDF               \n",REDF);
   printf("%e  SYMO2UPT   [mol/h] \n",SYMO2UPT*3600.);
   printf("%e  REDF               \n",REDF);
#endif
#ifdef CALCIUM
   printf("%e  CAUPT      [mol/h] \n",CAUPT*3600.);
#endif
   printf("======  UPTAKE (INPUT) end  ============== \n");
#endif

#if defined (FORAMOUD) || defined (FORAMOUL)
   printf("==== O. UNIVERSA: UPTAKE (INPUT)   ======= \n");
   printf("%e  SYMCO2UPT  [mol/h] \n",SYMCO2UPT*3600.);
   printf("%e  SYMHCO3UPT [mol/h] \n",SYMHCO3UPT*3600.);
   printf("%e  SYMHUPT    [mol/h] \n",SYMHUPT*3600.);
   printf("%e  CO3UPT     [mol/h] \n",CO3UPT*3600.);
   printf("%e  CO2UPT     [mol/h] \n",CO2UPT*3600.);
#ifdef C13ISTP
   printf("-------------- 13C ---------------------- \n");
   printf("%e  SYMCCO2UPT  [mol/h] \n",symcco2upt*3600.);
   printf("%e  SYMHCCO3UPT [mol/h] \n",symhcco3upt*3600.);
   printf("%e  SYMHUPT     [mol/h] \n",SYMHUPT*3600.);
   printf("%e  CCO3UPT     [mol/h] \n",cco3upt*3600.);
   printf("%e  CCO2UPT     [mol/h] \n",cco2upt*3600.);
#endif
#ifdef OXYGEN
   printf("%e  O2UPT      [mol/h] \n",-REDF*CO2UPT*3600.);
   printf("%e  REDF               \n",REDF);
   printf("%e  SYMO2UPT   [mol/h] \n",SYMO2UPT*3600.);
   printf("%e  REDF               \n",REDF);
#endif
#ifdef CALCIUM
   printf("%e  CAUPT      [mol/h] \n",CAUPT*3600.);
#endif
   printf("======  UPTAKE (INPUT) end  ============== \n");
#endif


#if defined (FORAMGSD) || defined (FORAMGSL) /* G. Sacculifer		*/
   printf("==== G. SACCULIFER: UPTAKE (INPUT)   ====== \n");
   printf("%e  SYMCO2UPT  [mol/h] \n",SYMCO2UPT*3600.);
   printf("%e  SYMHCO3UPT [mol/h] \n",SYMHCO3UPT*3600.);
   printf("%e  SYMHUPT    [mol/h] \n",SYMHUPT*3600.);
   printf("%e  CO3UPT     [mol/h] \n",CO3UPT*3600.);
   printf("%e  CO2UPT     [mol/h] \n",CO2UPT*3600.);
#ifdef C13ISTP
   printf("-------------- 13C ---------------------- \n");
   printf("%e  SYMCCO2UPT  [mol/h] \n",symcco2upt*3600.);
   printf("%e  SYMHCCO3UPT [mol/h] \n",symhcco3upt*3600.);
   printf("%e  SYMHUPT     [mol/h] \n",SYMHUPT*3600.);
   printf("%e  CCO3UPT     [mol/h] \n",cco3upt*3600.);
   printf("%e  CCO2UPT     [mol/h] \n",cco2upt*3600.);
#endif
#ifdef OXYGEN
   printf("%e  O2UPT      [mol/h] \n",-REDF*CO2UPT*3600.);
   printf("%e  REDF               \n",REDF);
   printf("%e  SYMO2UPT   [mol/h] \n",SYMO2UPT*3600.);
   printf("%e  REDS               \n",REDS);
#endif
#ifdef CALCIUM
   printf("%e  CAUPT      [mol/h] \n",CAUPT*3600.);
#endif
   printf("======  UPTAKE (INPUT) end  ============== \n");
#endif

#if defined (FORAMRES) || defined (FORAMPHS) || defined (FORAMCLC)
   printf("==== RES,PHS,CLC-EXP.: UPTAKE (INPUT)   ====== \n");
   printf("%e  SYMCO2UPT  [mol/h] \n",SYMCO2UPT*3600.);
   printf("%e  SYMHCO3UPT [mol/h] \n",SYMHCO3UPT*3600.);
   printf("%e  SYMHUPT    [mol/h] \n",SYMHUPT*3600.);
   printf("%e  CO3UPT     [mol/h] \n",CO3UPT*3600.);
   printf("%e  CO2UPT     [mol/h] \n",CO2UPT*3600.);
#ifdef C13ISTP
   printf("-------------- 13C ---------------------- \n");
   printf("%e  SYMCCO2UPT  [mol/h] \n",symcco2upt*3600.);
   printf("%e  SYMHCCO3UPT [mol/h] \n",symhcco3upt*3600.);
   printf("%e  SYMHUPT     [mol/h] \n",SYMHUPT*3600.);
   printf("%e  CCO3UPT     [mol/h] \n",cco3upt*3600.);
   printf("%e  CCO2UPT     [mol/h] \n",cco2upt*3600.);
#endif
#ifdef OXYGEN
   printf("%e  O2UPT      [mol/h] \n",-REDF*CO2UPT*3600.);
   printf("%e  REDF               \n",REDF);
   printf("%e  SYMO2UPT   [mol/h] \n",SYMO2UPT*3600.);
   printf("%e  REDS               \n",REDS);
#endif
#ifdef CALCIUM
   printf("%e  CAUPT      [mol/h] \n",CAUPT*3600.);
#endif
   printf("======  UPTAKE (INPUT) end  ============== \n");
#endif

#endif /* PRINT */

#ifdef FORAMW
   fprintf(fppara,"==== WORK FORAM: UPTAKE (INPUT)   ======= \n");
   fprintf(fppara,"SYMCO2UPT  %e  [mol/h] \n",SYMCO2UPT*3600.);
   fprintf(fppara,"SYMHCO3UPT %e  [mol/h] \n",SYMHCO3UPT*3600.);
   fprintf(fppara,"SYMHUPT    %e  [mol/h] \n",SYMHUPT*3600.);
   fprintf(fppara,"CO3UPT     %e  [mol/h] \n",CO3UPT*3600.);
   fprintf(fppara,"CO2UPT     %e  [mol/h] \n",CO2UPT*3600.);
#ifdef C13ISTP
   fprintf(fppara,"-------------- 13C ---------------------- \n");
   fprintf(fppara,"SYMCCO2UPT  %e  [mol/h] \n",symcco2upt*3600.);
   fprintf(fppara,"SYMHCCO3UPT %e  [mol/h] \n",symhcco3upt*3600.);
   fprintf(fppara,"SYMHUPT     %e  [mol/h] \n",SYMHUPT*3600.);
   fprintf(fppara,"CCO3UPT     %e  [mol/h] \n",cco3upt*3600.);
   fprintf(fppara,"CCO2UPT     %e  [mol/h] \n",cco2upt*3600.);
#endif
#ifdef OXYGEN
   fprintf(fppara,"O2UPT      %e  [mol/h] \n",-REDF*CO2UPT*3600.);
   fprintf(fppara,"REDF       %e          \n",REDF);
   fprintf(fppara,"SYMO2UPT   %e  [mol/h] \n",SYMO2UPT*3600.);
   fprintf(fppara,"REDF       %e          \n",REDF);
#endif
#ifdef CALCIUM
   fprintf(fppara,"CAUPT      %e  [mol/h] \n",CAUPT*3600.);
#endif
   fprintf(fppara,"======  UPTAKE (INPUT) end  ============== \n");
#endif


#if defined (FORAMBORD) || defined (FORAMBORL)
   fprintf(fppara,"==== BORON FORAM: UPTAKE (INPUT)   ======= \n");
   fprintf(fppara,"SYMCO2UPT  %e  [mol/h] \n",SYMCO2UPT*3600.);
   fprintf(fppara,"SYMHCO3UPT %e  [mol/h] \n",SYMHCO3UPT*3600.);
   fprintf(fppara,"SYMHUPT    %e  [mol/h] \n",SYMHUPT*3600.);
   fprintf(fppara,"CO3UPT     %e  [mol/h] \n",CO3UPT*3600.);
   fprintf(fppara,"CO2UPT     %e  [mol/h] \n",CO2UPT*3600.);
#ifdef C13ISTP
   fprintf(fppara,"-------------- 13C ---------------------- \n");
   fprintf(fppara,"SYMCCO2UPT  %e  [mol/h] \n",symcco2upt*3600.);
   fprintf(fppara,"SYMHCCO3UPT %e  [mol/h] \n",symhcco3upt*3600.);
   fprintf(fppara,"SYMHUPT     %e  [mol/h] \n",SYMHUPT*3600.);
   fprintf(fppara,"CCO3UPT     %e  [mol/h] \n",cco3upt*3600.);
   fprintf(fppara,"CCO2UPT     %e  [mol/h] \n",cco2upt*3600.);
#endif
#ifdef OXYGEN
   fprintf(fppara,"O2UPT      %e  [mol/h] \n",-REDF*CO2UPT*3600.);
   fprintf(fppara,"REDF       %e          \n",REDF);
   fprintf(fppara,"SYMO2UPT   %e  [mol/h] \n",SYMO2UPT*3600.);
   fprintf(fppara,"REDF       %e          \n",REDF);
#endif
#ifdef CALCIUM
   fprintf(fppara,"CAUPT      %e  [mol/h] \n",CAUPT*3600.);
#endif
   fprintf(fppara,"======  UPTAKE (INPUT) end  ============== \n");
#endif

#if defined (FORAMOUD) || defined (FORAMOUL)
   fprintf(fppara,"==== O. UNIVERSA: UPTAKE (INPUT)   ======= \n");
   fprintf(fppara,"SYMCO2UPT  %e  [mol/h] \n",SYMCO2UPT*3600.);
   fprintf(fppara,"SYMHCO3UPT %e  [mol/h] \n",SYMHCO3UPT*3600.);
   fprintf(fppara,"SYMHUPT    %e  [mol/h] \n",SYMHUPT*3600.);
   fprintf(fppara,"CO3UPT     %e  [mol/h] \n",CO3UPT*3600.);
   fprintf(fppara,"CO2UPT     %e  [mol/h] \n",CO2UPT*3600.);
#ifdef C13ISTP
   fprintf(fppara,"-------------- 13C ---------------------- \n");
   fprintf(fppara,"SYMCCO2UPT  %e  [mol/h] \n",symcco2upt*3600.);
   fprintf(fppara,"SYMHCCO3UPT %e  [mol/h] \n",symhcco3upt*3600.);
   fprintf(fppara,"SYMHUPT     %e  [mol/h] \n",SYMHUPT*3600.);
   fprintf(fppara,"CCO3UPT     %e  [mol/h] \n",cco3upt*3600.);
   fprintf(fppara,"CCO2UPT     %e  [mol/h] \n",cco2upt*3600.);
#endif
#ifdef OXYGEN
   fprintf(fppara,"O2UPT      %e  [mol/h] \n",-REDF*CO2UPT*3600.);
   fprintf(fppara,"REDF       %e          \n",REDF);
   fprintf(fppara,"SYMO2UPT   %e  [mol/h] \n",SYMO2UPT*3600.);
   fprintf(fppara,"REDF       %e          \n",REDF);
#endif
#ifdef CALCIUM
   fprintf(fppara,"CAUPT      %e  [mol/h] \n",CAUPT*3600.);
#endif
   fprintf(fppara,"======  UPTAKE (INPUT) end  ============== \n");
#endif


#if defined (FORAMGSD) || defined (FORAMGSL) /* G. Sacculifer		*/
   fprintf(fppara,"==== G. SACCULIFER: UPTAKE (INPUT)   ====== \n");
   fprintf(fppara,"SYMCO2UPT  %e  [mol/h] \n",SYMCO2UPT*3600.);
   fprintf(fppara,"SYMHCO3UPT %e  [mol/h] \n",SYMHCO3UPT*3600.);
   fprintf(fppara,"SYMHUPT    %e  [mol/h] \n",SYMHUPT*3600.);
   fprintf(fppara,"CO3UPT     %e  [mol/h] \n",CO3UPT*3600.);
   fprintf(fppara,"CO2UPT     %e  [mol/h] \n",CO2UPT*3600.);
#ifdef C13ISTP
   fprintf(fppara,"-------------- 13C ---------------------- \n");
   fprintf(fppara,"SYMCCO2UPT  %e  [mol/h] \n",symcco2upt*3600.);
   fprintf(fppara,"SYMHCCO3UPT %e  [mol/h] \n",symhcco3upt*3600.);
   fprintf(fppara,"SYMHUPT     %e  [mol/h] \n",SYMHUPT*3600.);
   fprintf(fppara,"CCO3UPT     %e  [mol/h] \n",cco3upt*3600.);
   fprintf(fppara,"CCO2UPT     %e  [mol/h] \n",cco2upt*3600.);
#endif
#ifdef OXYGEN
   fprintf(fppara,"O2UPT      %e  [mol/h] \n",-REDF*CO2UPT*3600.);
   fprintf(fppara,"REDF       %e          \n",REDF);
   fprintf(fppara,"SYMO2UPT   %e  [mol/h] \n",SYMO2UPT*3600.);
   fprintf(fppara,"REDS       %e          \n",REDS);
#endif
#ifdef CALCIUM
   fprintf(fppara,"CAUPT      %e  [mol/h] \n",CAUPT*3600.);
#endif
   fprintf(fppara,"======  UPTAKE (INPUT) end  ============== \n");
#endif

#if defined (FORAMRES) || defined (FORAMPHS) || defined (FORAMCLC)
   fprintf(fppara,"==== RES,PHS,CLC-EXP.: UPTAKE (INPUT)   ====== \n");
   fprintf(fppara,"SYMCO2UPT  %e  [mol/h] \n",SYMCO2UPT*3600.);
   fprintf(fppara,"SYMHCO3UPT %e  [mol/h] \n",SYMHCO3UPT*3600.);
   fprintf(fppara,"SYMHUPT    %e  [mol/h] \n",SYMHUPT*3600.);
   fprintf(fppara,"CO3UPT     %e  [mol/h] \n",CO3UPT*3600.);
   fprintf(fppara,"CO2UPT     %e  [mol/h] \n",CO2UPT*3600.);
#ifdef C13ISTP
   fprintf(fppara,"-------------- 13C ---------------------- \n");
   fprintf(fppara,"SYMCCO2UPT  %e  [mol/h] \n",symcco2upt*3600.);
   fprintf(fppara,"SYMHCCO3UPT %e  [mol/h] \n",symhcco3upt*3600.);
   fprintf(fppara,"SYMHUPT     %e  [mol/h] \n",SYMHUPT*3600.);
   fprintf(fppara,"CCO3UPT     %e  [mol/h] \n",cco3upt*3600.);
   fprintf(fppara,"CCO2UPT     %e  [mol/h] \n",cco2upt*3600.);
#endif
#ifdef OXYGEN
   fprintf(fppara,"O2UPT      %e  [mol/h] \n",-REDF*CO2UPT*3600.);
   fprintf(fppara,"REDF       %e          \n",REDF);
   fprintf(fppara,"SYMO2UPT   %e  [mol/h] \n",SYMO2UPT*3600.);
   fprintf(fppara,"REDS       %e          \n",REDS);
#endif
#ifdef CALCIUM
   fprintf(fppara,"CAUPT      %e  [mol/h] \n",CAUPT*3600.);
#endif
   fprintf(fppara,"======  UPTAKE (INPUT) end  ============== \n");
#endif




   fprintf(fppara,"==================================== \n");

#ifdef COCCOW
   fprintf(fppara,"============= COCCOW =============== \n");
#endif
#ifdef COCCO1
   fprintf(fppara,"============= COCCO1 =============== \n");
#endif
#ifdef COCCO2
   fprintf(fppara,"============= COCCO2 =============== \n");
#endif

#ifdef DIATOMW
   fprintf(fppara,"============= DIATOMW =============== \n");
#endif

#ifdef DIATOM1
   fprintf(fppara,"============= DIATOM1 =============== \n");
#endif
#ifdef DIATOM2
   fprintf(fppara,"============= DIATOM2 =============== \n");
#endif
#ifdef DIATOM3
   fprintf(fppara,"============= DIATOM3 =============== \n");
#endif
#ifdef DIATOM4
   fprintf(fppara,"============= DIATOM4 =============== \n");
#endif
#ifdef DIATOM5
   fprintf(fppara,"============= DIATOM5 =============== \n");
#endif

#ifdef FORAM1
   fprintf(fppara,"============= FORAM1 =============== \n");
#endif
#ifdef FORAM2
   fprintf(fppara,"============= FORAM2 =============== \n");
#endif

   fprintf(fppara,"========================================= ");
   fprintf(fppara,"============= initk start =============== \n");
   fprintf(fppara,"======== rate coefficients ============== \n");
   fprintf(fppara,"k12    [kg/mumol/s]              %e \n",k12);
   fprintf(fppara,"k21    [1/s]                     %e \n",k21);
   fprintf(fppara,"k13    [kg/mumol/s]              %e \n",k13);
   fprintf(fppara,"k31    [1/s]                     %e \n",k31);
   fprintf(fppara,"kp4    [kg/mumol/s]              %e \n",kp4);
   fprintf(fppara,"km4    [1/s]                     %e \n",km4);
   fprintf(fppara,"kp5h   [kg/mumol/s]              %e \n",kp5h);
   fprintf(fppara,"km5h   [1/s]                     %e \n",km5h);
   fprintf(fppara,"kp6    [mumol/kg/s]              %e \n",kp6);
   fprintf(fppara,"km6    [kg/mumol/s]              %e \n",km6);
#ifdef BORONRC4
   fprintf(fppara,"kp7    [kg/mumol/s]              %e \n",kp7);
   fprintf(fppara,"km7    [1/s]                     %e \n",km7);
#endif
#ifdef BORONRC3
   fprintf(fppara,"kp7    [1/s]                     %e \n",kp7);
   fprintf(fppara,"km7    [kg/mumol/s]              %e \n",km7);
#endif
#ifdef C13ISTP
   fprintf(fppara,"----------------  C  -------------- \n");
   fprintf(fppara,"kp1s   [1/s]                     %e \n",kp1s);
   fprintf(fppara,"km1s   [kg/mumol/s]              %e \n",km1s);
   fprintf(fppara,"kp4    [kg/mumol/s]              %e \n",kp4);
   fprintf(fppara,"km4    [1/s]                     %e \n",km4);
   fprintf(fppara,"kp5h   [kg/mumol/s]              %e \n",kp5h);
   fprintf(fppara,"km5h   [1/s]                     %e \n",km5h);
   fprintf(fppara,"----------------  13C  ------------ \n");
   fprintf(fppara,"kp1scc [1/s]                     %e \n",kp1scc);
   fprintf(fppara,"km1scc [kg/mumol/s]              %e \n",km1scc);
   fprintf(fppara,"kp4cc  [kg/mumol/s]              %e \n",kp4cc);
   fprintf(fppara,"km4cc  [1/s]                     %e \n",km4cc);
   fprintf(fppara,"kp5cc  [kg/mumol/s]              %e \n",kp5cc);
   fprintf(fppara,"km5cc  [1/s]                     %e \n",km5cc);
   fprintf(fppara,"----------------------------------- \n");
#endif
   fprintf(fppara,"kh2co3 [1/s]                     %e \n",kh2co3);
   fprintf(fppara,"======== equi. coefficients ============= \n");
   fprintf(fppara,"Kh2co3 [mumol/kg]               %e \n",Kh2co3);
   fprintf(fppara,"======== diff. coefficients ============== \n");
   fprintf(fppara,"diff. coeff. CO2   [mu**2/s]   %e \n",dco2);
   fprintf(fppara,"diff. coeff. CO2   [m**2/s]    %e \n",dco2/1.0e12);
   fprintf(fppara,"diff. coeff. HCO3- [m**2/s]    %e \n",dhco3/1.0e12);
   fprintf(fppara,"diff. coeff. CO3-- [m**2/s]    %e \n",dco3/1.0e12);
   fprintf(fppara,"diff. coeff. H+    [m**2/s]    %e \n",dh/1.0e12);
   fprintf(fppara,"diff. coeff. OH-   [m**2/s]    %e \n",doh/1.0e12);
   fprintf(fppara,"dicosw(T) (correction)         %e \n",dicosw(TK));
   /* Diffusion coefficients in seperate file */
   fprintf(fpdc,"%e\n",dco2/1.0e12);
   fprintf(fpdc,"%e\n",dhco3/1.0e12);
   fprintf(fpdc,"%e\n",dco3/1.0e12);
   fprintf(fpdc,"%e\n",dh/1.0e12);
   fprintf(fpdc,"%e\n",doh/1.0e12);
#ifdef C13ISTP
   fprintf(fpdcc,"%e\n",dcco2/1.0e12);
   fprintf(fpdcc,"%e\n",dhcco3/1.0e12);
   fprintf(fpdcc,"%e\n",dcco3/1.0e12);
#endif
#ifdef OXYGEN
   fprintf(fpdc,"%e\n",do2/1.0e12);
#endif
#ifdef CALCIUM
   fprintf(fpdc,"%e\n",dca/1.0e12);
#endif

#ifdef CALCIUM
   fprintf(fppara,"diff. coeff. Ca++  [m**2/s]    %e \n",dca/1.0e12);
   fprintf(fppara,"diff. coeff. ATA   [m**2/s]    %e \n",data/1.0e12);
#endif
#ifdef OXYGEN
   fprintf(fppara,"diff. coeff. O2    [m**2/s]    %e \n",do2/1.0e12);
#endif
#ifdef BORON
   fprintf(fppara,"diff. coeff. BOH3  [m**2/s]    %e \n",dboh3/1.0e12);
   fprintf(fppara,"diff. coeff. BOH4  [m**2/s]    %e \n",dboh4/1.0e12);
#endif
#ifdef C13ISTP
   fprintf(fppara,"diff. coeff. 13CO2 [m**2/s]    %e \n",dcco2/1.0e12);
   fprintf(fppara,"diff. coeff. 13HCO3[m**2/s]    %e \n",dhcco3/1.0e12);
   fprintf(fppara,"diff. coeff. 13CO3 [m**2/s]    %e \n",dcco3/1.0e12);
#endif


   fprintf(fppara,"=====================");
   fprintf(fppara,"============= initk end =============== \n");

}   /* --- end of initk --- */


#ifdef CLPL
void initdm()
{

/* ============= init matrix of diffusion coefficients ============= */


int jvar,k,sumkmem;

 /* initialize matrix of diffusion coefficients */
 /* first index NVAR, second k (grid)		*/
 /* grid: k=1 L.B., 2<=k<=M interior, k=M+1 R.B.*/

 for(jvar=0; jvar<=N2; jvar++)
   	difcofm[jvar][0] = 0.0;	/* not used */

 for(k=1; k<=M; k++){
   difcofm[0][k] = 0.0;		/* not used */
#ifdef OXYGEN
   difcofm[EQO2][k] = do2;
#endif

#ifdef C13ISTP
   difcofm[EQCCO2][k] = dcco2;
   difcofm[EQHCCO3][k] = dhcco3;
   difcofm[EQCCO3][k] = dcco3;
#endif

#ifdef BORON
   difcofm[EQBOH3][k] = dboh3;
   difcofm[EQBOH4][k] = dboh4;
#endif

   difcofm[EQCO2][k] = dco2;
   difcofm[EQHCO3][k] = dhco3;
   difcofm[EQCO3][k] = dco3;
   difcofm[EQHP][k] = dh;
   difcofm[EQOH][k] = doh;

#ifdef CALCIUM
   difcofm[EQCA] = dca;
#endif

    if(k == 2){
	printf("\n %e dco2\n",difcofm[EQCO2][k]*1.e-12);
	printf(" %e dhco3\n",difcofm[EQHCO3][k]*1.e-12);
	printf(" %e dco3\n",difcofm[EQCO3][k]*1.e-12);
	printf(" %e dh\n",difcofm[EQHP][k]*1.e-12);
	printf(" %e doh\n",difcofm[EQOH][k]*1.e-12);
	printf(" %e dboh3\n",difcofm[EQBOH3][k]*1.e-12);
	printf(" %e dboh4\n",difcofm[EQBOH4][k]*1.e-12);
    }

#define LINDIF
#ifdef LINDIF
   /* linear increase of diff.-coefficient towards cell membrane */
   sumkmem = (ncmembmax-ncmembmin);
   dumk1    = ncmembmin - (ncmembmax-ncmembmin);
   if(k >= dumk1 && k < ncmembmin){
#ifdef OXYGEN
   difcofm[EQO2][k]    = do2 - (k-dumk1)*(do2-PERMO2*MEMBTHK*2.)/(sumkmem);
#endif

#ifdef C13ISTP
   difcofm[EQCCO2][k]  = dcco2  - (k-dumk1)*(dcco2-PERMCO2*MEMBTHK*2.)/(sumkmem);
   difcofm[EQHCCO3][k] = dhcco3 - (k-dumk1)*(dhcco3-PERMION*MEMBTHK*2.)/(sumkmem);
   difcofm[EQCCO3][k]  = dcco3  - (k-dumk1)*(dcco3-PERMION*MEMBTHK*2.)/(sumkmem);
#endif


#ifdef BORON
   difcofm[EQBOH3][k] = dboh3  - (k-dumk1)*(dboh3-PERMBOH3*MEMBTHK*2.)/(sumkmem);
   difcofm[EQBOH4][k] = dboh4  - (k-dumk1)*(dboh4-PERMION*MEMBTHK*2.)/(sumkmem);
#endif
   difcofm[EQCO2][k]  = dco2  - (k-dumk1)*(dco2-REDDCO2)/(sumkmem);
   difcofm[EQHCO3][k] = dhco3 - (k-dumk1)*(dhco3-PERMION*MEMBTHK*2.)/(sumkmem);
   difcofm[EQCO3][k]  = dco3  - (k-dumk1)*(dco3-PERMION*MEMBTHK*2.)/(sumkmem);
   difcofm[EQHP][k]   = dh    - (k-dumk1)*(dh-PERMION*MEMBTHK*2.)/(sumkmem);
   difcofm[EQOH][k]   = doh   - (k-dumk1)*(doh-PERMION*MEMBTHK*2.)/(sumkmem);


#ifdef CALCIUM
   difcofm[EQCA][k]   = dca   - (k-dumk1)*(dca-PERMION*MEMBTHK*2.)/(sumkmem);
#endif
   }   /* end if */
#endif /* LINDIF */



   /* reduce diffusion coefficients within cell membrane */
   if(k >= ncmembmin && k < ncmembmax){
#ifdef OXYGEN
   difcofm[EQO2][k]    = PERMO2*MEMBTHK*2.;
#endif

#ifdef C13ISTP
   difcofm[EQCCO2][k]  = PERMCO2*MEMBTHK*2.;
   difcofm[EQHCCO3][k] = PERMION*MEMBTHK*2.;
   difcofm[EQCCO3][k]  = PERMION*MEMBTHK*2.;
#endif


#ifdef BORON
   difcofm[EQBOH3][k] = PERMBOH3*MEMBTHK*2.;
   difcofm[EQBOH4][k] = PERMION*MEMBTHK*2.;
#endif

   difcofm[EQCO2][k]  = REDDCO2;

   /* difcofm[EQCO2][k]  = PERMCO2*MEMBTHK*2.; */
   difcofm[EQHCO3][k] = PERMION*MEMBTHK*2.;
   difcofm[EQCO3][k]  = PERMION*MEMBTHK*2.;
   difcofm[EQHP][k]   = PERMION*MEMBTHK*2.;
   difcofm[EQOH][k]   = PERMION*MEMBTHK*2.;

    if(k == ncmembmin){
	printf("\n %e dco2\n",difcofm[EQCO2][k]*1.e-12);
	printf(" %e dhco3\n",difcofm[EQHCO3][k]*1.e-12);
	printf(" %e dco3\n",difcofm[EQCO3][k]*1.e-12);
	printf(" %e dh\n",difcofm[EQHP][k]*1.e-12);
	printf(" %e doh\n",difcofm[EQOH][k]*1.e-12);
	printf(" %e dboh3\n",difcofm[EQBOH3][k]*1.e-12);
	printf(" %e dboh4\n",difcofm[EQBOH4][k]*1.e-12);
    }
#ifdef CALCIUM
   difcofm[EQCA][k]   = PERMION*MEMBTHK*2.;
#endif
   } /* end  if  */

#ifdef LINDIF
   /* linear decrease of diff-coefficient to cell membrane */
   dumk2    = ncmembmax + (ncmembmax-ncmembmin);
   if(k >= ncmembmax  && k < dumk2){

#ifdef OXYGEN
   difcofm[EQO2][k]    = PERMO2*MEMBTHK*2.
   + (k-ncmembmax)*(dco2-PERMO2*MEMBTHK*2.)/(sumkmem);
#endif

#ifdef C13ISTP
   difcofm[EQCCO2][k]  = PERMCO2*MEMBTHK*2.
	+ (k-ncmembmax)*(dcco2-PERMION*MEMBTHK*2.)/(sumkmem);
   difcofm[EQHCCO3][k] = PERMION*MEMBTHK*2.
	+ (k-ncmembmax)*(dhcco3-PERMION*MEMBTHK*2.)/(sumkmem);
   difcofm[EQCCO3][k]  = PERMION*MEMBTHK*2.
	+ (k-ncmembmax)*(dcco3-PERMION*MEMBTHK*2.)/(sumkmem);
#endif


#ifdef BORON
   difcofm[EQBOH3][k] = PERMBOH3*MEMBTHK*2.
	+ (k-ncmembmax)*(dboh3-PERMBOH3*MEMBTHK*2.)/(sumkmem);
   difcofm[EQBOH4][k] = PERMION*MEMBTHK*2.
	+ (k-ncmembmax)*(dboh4-PERMION*MEMBTHK*2.)/(sumkmem);
#endif
   difcofm[EQCO2][k]  = REDDCO2
	+ (k-ncmembmax)*(dco2-REDDCO2)/(sumkmem);
   difcofm[EQHCO3][k] = PERMION*MEMBTHK*2.
	+ (k-ncmembmax)*(dhco3-PERMION*MEMBTHK*2.)/(sumkmem);
   difcofm[EQCO3][k]  = PERMION*MEMBTHK*2.
	+ (k-ncmembmax)*(dco3-PERMION*MEMBTHK*2.)/(sumkmem);
   difcofm[EQHP][k]   = PERMION*MEMBTHK*2.
	+ (k-ncmembmax)*(dh-PERMION*MEMBTHK*2.)/(sumkmem);
   difcofm[EQOH][k]   = PERMION*MEMBTHK*2.
	+ (k-ncmembmax)*(doh-PERMION*MEMBTHK*2.)/(sumkmem);
#ifdef CALCIUM
   difcofm[EQCA][k]   = PERMION*MEMBTHK*2.
	+ (k-ncmembmax)*(dca-PERMION*MEMBTHK*2.)/(sumkmem);
#endif
   }   /* end  if  */
#endif /* LINDIF */

   }   /* end  for */


   for(jvar=1;jvar<=EQOH;jvar++)
   	printf("%e difcofm[%d][memb]\n",difcofm[jvar][ncmembmin]*1.e-12,jvar);

}	/* end initdm() */
#endif  /* CLPL 	*/

#ifdef INITA

void inita()
{
/* =============== init all arrays ======================= */

   /* ---------- concentrations ---------- */

   for(ir = 0; ir <= M; ir++) {
      co2[ir]   = co2bulk;
      hco3[ir]  = hco3bulk;
      co3[ir]   = co3bulk;
      hplus[ir] = hbulk;
      oh[ir]    = ohbulk;
      h2o[ir]   = h2obulk;
#ifdef C13ISTP
      cco2[ir]    = cco2bulk;
      hcco3[ir]   = hcco3bulk;
      cco3[ir]    = cco3bulk;
#endif
#ifdef OXYGEN
      o2[ir]    = o2bulk;
#endif
#ifdef BORON
      boh3[ir]  = boh3bulk;
      boh4[ir]  = boh4bulk;
#ifdef BORISTP
      bboh3[ir]  = bboh3bulk;
      bboh4[ir]  = bboh4bulk;
#endif
#endif

#ifdef CALCIUM
      ca[ir]    = cabulk;
#endif
   }

#ifdef ALKALI
     for(ir = 0; ir < IMAX; ir++)
       ata[ir] = -(2. * co3[ir] + hco3[ir] + oh[ir] - hplus[ir]);
#ifdef CISTP
     for(ir = 0; ir < IMAX; ir++)
       ata[ir] = -(2. * (co3[ir]+cco3[ir])
       		 +(hco3[ir]+hcco3[ir])
       		 + oh[ir] - hplus[ir]);
#endif
#endif

#if defined (ALKALI) && defined (BORON)
     for(ir = 0; ir < IMAX; ir++)
       ata[ir] = -(2. * co3[ir] + hco3[ir] + boh4[ir] + oh[ir] - hplus[ir]);
#ifdef CISTP
     for(ir = 0; ir < IMAX; ir++)
       ata[ir] = -(2. * (co3[ir]+cco3[ir])
       		 +(hco3[ir]+hcco3[ir]) + boh4[ir]
       		 + oh[ir] - hplus[ir]);
#endif
#endif

   surface = 4. * PI * RADIUS * RADIUS;

  /* ----------------------- units:
   #define CO2UPT (5.e-17)  [mol/s * C/cell]
   #define DCO2 (1.6e-9*1.0e12)     [mu**2/s] diffusion coefficient
   units: co2flux [mumol/kg/mu] =
    CO2UPT [mol/s] * 1.e6 (mol -> mumol)
      / surface [mu**2] / DCO2 [mu**2/s] * 1.0e15 ( mu**3 -> l)
  */

   dummy = HUPT - HCO3UPT - 2.*CO3UPT - OHUPT + 2.*CAUPT;   /* neutral? */

#ifdef AGG
   c0          = pow((RBULK/REHUX),(1./0.44))*CEHUX; 	   /* mol C */
#ifndef SYMCO2UPT
   symco2upt   = -c0*exp(-TIME/TAU)/TAU;
#endif
#ifndef SYMO2UPT
   symo2upt    = -REDAGG*symco2upt;
#endif
#ifndef CO2UPT
   c0 	    = pow((RBULK/REHUX),(1./0.44))*CEHUX; /* mol C */
   co2upt   = -c0*exp(-TIME/TAU)/TAU;
   /*printf("\n %e c0\n",c0);
     printf("\n %e co2upt\n",co2upt);*/
   co2flux  = co2upt  / surface / dco2  * 1.e21;
#endif
#else
   co2flux  = CO2UPT  / surface / dco2  * 1.e21;
#endif
   hflux    = HUPT    / surface / dh    * 1.e21;
   hco3flux = HCO3UPT / surface / dhco3 * 1.e21;
   co3flux  = CO3UPT  / surface / dco3  * 1.e21;

#ifdef LDAT
   co3flux  = co3uptldat  / surface / dco3  * 1.e21;
#endif
   ohflux   = OHUPT   / surface / doh   * 1.e21;

#ifdef CISTP
   /* cco3upt was set in initk() as first guess */
   co2flux  = (CO2UPT - cco2upt)  / surface / dco2  * 1.e21;
   hco3flux = (HCO3UPT-hcco3upt)  / surface / dhco3 * 1.e21;
   co3flux  = (CO3UPT - cco3upt)  / surface / dco3  * 1.e21;
#ifdef LDAT
   co3flux  = (co3uptldat - cco3upt)  / surface / dco3  * 1.e21;
#endif
#endif

   co2fluxs  = co2flux;
   hfluxs    = hflux;
   hco3fluxs = hco3flux;
   co3fluxs  = co3flux;
   ohfluxs   = ohflux;
#ifdef C13ISTP
   cco2flux  = cco2upt  / surface / dcco2  * 1.e21;
   hcco3flux = hcco3upt / surface / dhcco3 * 1.e21;
   cco3flux  = cco3upt  / surface / dcco3  * 1.e21;

   cco2fluxs  = cco2flux;
   hcco3fluxs = hcco3flux;
   cco3fluxs  = cco3flux;
#endif

#ifdef CALCIUM
   caflux   = CAUPT   / surface / dca   * 1.e21;
   cafluxs  = caflux;
#endif
#ifdef OXYGEN
 #ifdef O2UPT
   o2flux =       O2UPT / surface / do2  * 1.e21;
 #endif
 #ifdef REDFYEA
   o2flux = -REDF*CO2UPT/ surface / do2  * 1.e21;
			  /* REDF: Redfield ratio O2=R*CO2 of foram resp. */
			  /* this is O2UPT at the shell		    */
 #endif
 #ifdef AGG
   o2flux =     -co2upt / surface / do2  * 1.e21;
   /*printf("%e co2flux\n",co2flux);
     printf("%e  o2flux\n",o2flux);*/
 #ifdef O2UPT
   o2flux =       O2UPT / surface / do2  * 1.e21;
 #endif
 #endif
   o2fluxs = o2flux;
#endif
#ifdef BORON
   boh3flux  = 0.0;
   boh4flux  = BOH4UPT / surface / dboh4  * 1.e21;
#ifdef BORISTP
#ifndef B10B11
   d11boh4upt = d11boh4bulk;
   bboh3flux  = 0.0;
   bboh4flux  = boh4flux*(1./(1.+1./(BSTAND*(1.+d11boh4upt/1000.))));
#endif
#ifdef B10B11
   bboh3flux  = 0.0;
   boh4flux   = BOH4UPT / surface / dboh4  * 1.e21
   		/(BSTAND*(1.+d11boh4upt/1000.));
   bboh4flux  = BOH4UPT / surface / dboh4  * 1.e21
   		*(1./(1.+1./(BSTAND*(1.+d11boh4upt/1000.))));
#endif
   bboh3fluxs = bboh3flux;
   bboh4fluxs = bboh4flux;
#endif
   boh3fluxs = boh3flux;
   boh4fluxs = boh4flux;
#endif



   fprintf(fppara,"=====================");
   fprintf(fppara,"============= inita start =============== \n");
   fprintf(fppara,"TKELVIN (input)          %f \n",TKELVIN);
   fprintf(fppara,"temperature [K] (input)  %f \n",TK);
   fprintf(fppara,"cell radius [mu] (input) %f \n",RADIUS);
   fprintf(fppara,"bulk radius [mu] (input) %f \n",RBULK);
#ifdef AGG
   fprintf(fppara,"CO2   uptake   [mol/s] (input) %e \n",co2upt);
#else
   fprintf(fppara,"TCO2  uptake   [mol/s] (input) %e \n",TCO2UPT);
   fprintf(fppara,"CO2   uptake   [mol/s] (input) %e \n",CO2UPT);
#endif
   fprintf(fppara,"HCO3- uptake   [mol/s] (input) %e \n",HCO3UPT);
   fprintf(fppara,"CO3-- uptake   [mol/s] (input) %e \n",CO3UPT);
#ifdef C13ISTP
   fprintf(fppara,"13CO2   uptake [mol/s] (input) %e \n",cco2upt);
   fprintf(fppara,"13HCO3- uptake [mol/s] (input) %e \n",hcco3upt);
   fprintf(fppara,"13CO3-- uptake [mol/s] (input) %e \n",cco3upt);
#endif
#ifdef CALCIUM
   fprintf(fppara,"CA++  uptake [mol/s] (input) %e \n",CAUPT);
#endif
   fprintf(fppara,"OH-   uptake [mol/s] (input) %e \n",OHUPT);
   fprintf(fppara,"H+    uptake [mol/s] (input) %e \n",HUPT);
#ifdef CALCIUM
   fprintf(fppara,"CA++ bulk [mol/kg] (input)         %e \n",CABULK);
#endif
   fprintf(fppara,"neutral?     [mol/s/cell] (input) %e \n",dummy);
   fprintf(fppara,"M grid points (input)     %d \n",M);
   fprintf(fppara,"co2flux  [mol/kg/mu]      %e \n",co2flux);
#ifdef OXYGEN
   fprintf(fppara,"o2flux   [mol/kg/mu]      %e \n",o2flux);
#endif
#ifdef BORON
   fprintf(fppara,"boh3flux [mol/kg/mu]      %e \n",boh3flux);
   fprintf(fppara,"boh4flux [mol/kg/mu]      %e \n",boh4flux);
#endif
   fprintf(fppara,"hflux    [mol/kg/mu]      %e \n",hflux);
   fprintf(fppara,"hco3flux [mol/kg/mu]      %e \n",hco3flux);
   fprintf(fppara,"co3flux  [mol/kg/mu]      %e \n",co3flux);
   fprintf(fppara,"ohflux   [mol/kg/mu]      %e \n",ohflux);
#ifdef C13ISTP
   fprintf(fppara,"cco2flux  [mol/kg/mu]      %e \n",cco2flux);
   fprintf(fppara,"hcco3flux [mol/kg/mu]      %e \n",hcco3flux);
   fprintf(fppara,"cco3flux  [mol/kg/mu]      %e \n",cco3flux);
#endif
#ifdef CALCIUM
   fprintf(fppara,"caflux    [mol/kg/mu]      %e \n",caflux);
#endif
   fprintf(fppara,"surface area [mu**2]    %e \n",surface);
#ifdef OXYGEN
   fprintf(fppara,"REDF: Redf. foram resp. %e \n",REDF);
#endif
   fprintf(fppara,"=====================");
   fprintf(fppara,"============= inita end =============== \n");

}   /* ----- end of inita ----- */
#endif


/* =========================================================
   =========================================================

               init (end)

   =========================================================
   ========================================================= */




/* =========================================================
   =========================================================

             Numerical Recipes   (begin)

   =========================================================
   ========================================================= */


void solvde(itmax,conv,slowc,scalv,indexv,ne,nb,m,y,c,s)
       /* ----- 6/93 dwg Numerical Recipes: float -> double ----- */
int itmax,ne,nb,m;
double conv,slowc,scalv[],**y,***c,**s;
int indexv[];
{
	int ic1,ic2,ic3,ic4,it,j,j1,j2,j3,j4,j5,j6,j7,j8,j9;
	int jc1,jcf,jv,k,k1,k2,km,kp,nvars,*kmax,*ivector();
	double err,errj,fac,vmax,vz,*ermax,*dvector(),x;
	void pinvs(),difeq(),red(),bksub(),nrerror(),free_dvector(),
		free_ivector();
	kmax=ivector(1,ne);
	ermax=dvector(1,ne);
	k1=1;
	k2=m;
	nvars=ne*m;
	j1=1;
	j2=nb;
	j3=nb+1;
	j4=ne;
	j5=j4+j1;
	j6=j4+j2;
	j7=j4+j3;
	j8=j4+j4;
	j9=j8+j1;
	ic1=1;
	ic2=ne-nb;
	ic3=ic2+1;
	ic4=ne;
	jc1=1;
	jcf=ic3;
	for (it=1;it<=itmax;it++) {

/* -----   store data: -> difeq   ----- */

          for(j=1; j <= M; j++) {
            co2[j] = y[EQCO2][j];
           hco3[j] = y[EQHCO3][j];
#ifdef EQCO3
            co3[j] = y[EQCO3][j];
          hplus[j] = y[EQHP][j];
             oh[j] = y[EQOH][j];
#endif
#ifdef C13ISTP
           cco2[j] = y[EQCCO2][j];
          hcco3[j] = y[EQHCCO3][j];
           cco3[j] = y[EQCCO3][j];
#endif
#ifdef OXYGEN
             o2[j] = y[EQO2][j];
#endif
#ifdef BORON
           boh3[j] = y[EQBOH3][j];
           boh4[j] = y[EQBOH4][j];
#ifdef BORISTP
           bboh3[j] = y[EQBBOH3][j];
           bboh4[j] = y[EQBBOH4][j];
#endif
#endif
#ifdef CALCIUM
             ca[j] = y[EQCA][j];
#endif
          }

		co2negflag = 0;
		for(j=1;j<=M;j++){
		 	if(y[EQCO2][j] < 0.0) co2negflag = 1;
		}
		if(co2negflag == 1)
			printf("\n ! too bad - CO2 is negative !\n");
#ifdef MIMECO2SYM
		/* set vmaxit: linear increase with step of iteration (it)*/
		/* to avoid negative values of co2. Neg. values occur	*/
		/* if the initial uptake is too large.			*/
		/* Vmax(it) = a * it					*/
		/* slope: a = dV/d(it)					*/

		if(it <= (int)(vmaxco2*1.e9*3600./DVDIT))
			vmaxit = DVDIT*(double)(it*1.e-9/3600.);
		if(it >  (int)(vmaxco2*1.e9*3600./DVDIT))
			vmaxit = vmaxco2;

		printf("\n-----  before iteration ------\n");
		printf("%d co2negflag\n",co2negflag);
		printf("%d it\n",it);
		printf("%e vmaxit\n",vmaxit*3600.);
#endif

		k=k1;
		difeq(k,k1,k2,j9,ic3,ic4,indexv,ne,s,y);
		pinvs(ic3,ic4,j5,j9,jc1,k1,c,s);
#if defined (CLPL) && defined (DRAIN)
		calldifeq = 1;
		fdrain = 0.0;
		for (k=k1+1;k<=k2;k++) {
			difeq(k,k1,k2,j9,ic1,ic4,indexv,ne,s,y);
		}
		calldifeq = 2;
		printf("calldifeq %d %e \n",calldifeq,fdrain);
		for (k=k1+1;k<=k2;k++) {
			kp=k-1;
			difeq(k,k1,k2,j9,ic1,ic4,indexv,ne,s,y);
			red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp,c,s);
			pinvs(ic1,ic4,j3,j9,jc1,k,c,s);
		}
#else
		for (k=k1+1;k<=k2;k++) {
			kp=k-1;
			difeq(k,k1,k2,j9,ic1,ic4,indexv,ne,s,y);
			red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp,c,s);
			pinvs(ic1,ic4,j3,j9,jc1,k,c,s);
		}
#endif
		k=k2+1;
		difeq(k,k1,k2,j9,ic1,ic2,indexv,ne,s,y);
		red(ic1,ic2,j5,j6,j7,j8,j9,ic3,jc1,jcf,k2,c,s);
		pinvs(ic1,ic2,j7,j9,jcf,k2+1,c,s);
		bksub(ne,nb,jcf,k1,k2,c);
		err=0.0;
		for (j=1;j<=ne;j++) {
			jv=indexv[j];
			errj=vmax=0.0;
			km=0;
			for (k=k1;k<=k2;k++) {
				vz=fabs(c[j][1][k]);
				if (vz > vmax) {
					 vmax=vz;
					 km=k;
				}
				errj += vz;
			}
			err += errj/scalv[jv];
			ermax[j]=c[j][1][km]/scalv[jv];
			kmax[j]=km;
		}
		err /= nvars;
		fac=(err > slowc ? slowc/err : 1.0);
              /* set new values */

		for (jv=1;jv<=ne;jv++) {
			j=indexv[jv];
			for (k=k1;k<=k2;k++)
				y[j][k] -= fac*c[jv][1][k];
		}

#ifdef NOHPLUS
          for (k=k1;k<=k2;k++) y[4][k] = hplus[k];  /* diagnostic */
#endif

		printf("\n%8s %9s %9s\n","Iter.","Error","FAC");
		printf("%6d %12.6f %11.6f\n",it,err,fac);
#ifdef PRINT
		printf("%8s %8s %14s\n","Var.","Kmax","Max. Error");
		for (j=1;j<=ne;j++)
			printf("%6d %9d %14.6f \n",indexv[j],kmax[j],ermax[j]);
#endif

#ifdef MIMECO2SYM
		if (err < conv && vmaxit >= vmaxco2) {
			free_dvector(ermax,1,ne);
			free_ivector(kmax,1,ne);
			return;
		}
#else
		if (err < conv) {
			free_dvector(ermax,1,ne);
			free_ivector(kmax,1,ne);
			return;
		}
#endif


	}

/*        debug   (only if too many iterations in SOLVDE)   */


   fpdco2   = fopen("dco2.sv4","w");
   fpdhco3  = fopen("dhco3.sv4","w");
   fpdco3   = fopen("dco3.sv4","w");
   fpdh     = fopen("dh.sv4","w");
   fpdoh    = fopen("doh.sv4","w");
#ifdef C13ISTP
   fpdcco2   = fopen("dcco2.sv4","w");
   fpdhcco3  = fopen("dhcco3.sv4","w");
   fpdcco3   = fopen("dcco3.sv4","w");
#endif
#ifdef OXYGEN
   fpdo2    = fopen("do2.sv4","w");
#endif
#ifdef BORON
   fpdboh3  = fopen("dboh3.sv4","w");
   fpdboh4  = fopen("dboh4.sv4","w");
#ifdef BORISTP
   fpdbboh3  = fopen("dbboh3.sv4","w");
   fpdbboh4  = fopen("dbboh4.sv4","w");
#endif
#endif
#ifdef CALCIUM
   fpdca    = fopen("dca.sv4","w");
#endif


      for(j=1;j<=M;j++) {
        fprintf(fpr,"%f\n",r[j]);
        fprintf(fpco2,"%e\n", y[EQCO2][j]);
        fprintf(fphco3,"%e\n",y[EQHCO3][j]);
        fprintf(fpco3,"%e\n", y[EQCO3][j]);
        fprintf(fph,"%e\n",   y[EQHP][j]);
        fprintf(fpoh,"%e\n",  y[EQOH][j]);
#ifdef C13ISTP
        fprintf(fpcco2,"%e\n", y[EQCCO2][j]);
        fprintf(fphcco3,"%e\n",y[EQHCCO3][j]);
        fprintf(fpcco3,"%e\n", y[EQCCO3][j]);
#endif
#ifdef OXYGEN
        fprintf(fpo2,"%e\n",  y[EQO2][j]);
#endif
#ifdef BORON
        fprintf(fpboh3,"%e\n",  y[EQBOH3][j]);
        fprintf(fpboh4,"%e\n",  y[EQBOH4][j]);
#ifdef BORISTP
        fprintf(fpbboh3,"%e\n",  y[EQBBOH3][j]);
        fprintf(fpbboh4,"%e\n",  y[EQBBOH4][j]);
#endif
#endif
#ifdef CALCIUM
        fprintf(fpca,"%e\n",  y[EQCA][j]);
#endif

        fprintf(fpdco2,"%e\n", y[N2+EQCO2][j]);
        fprintf(fpdhco3,"%e\n",y[N2+EQHCO3][j]);
        fprintf(fpdco3,"%e\n", y[N2+EQCO3][j]);
        fprintf(fpdh,"%e\n",   y[N2+EQHP][j]);
        fprintf(fpdoh,"%e\n",  y[N2+EQOH][j]);
#ifdef C13ISTP
 	    fprintf(fpdcco2,"%e\n", y[N2+EQCCO2][j]);
        fprintf(fpdhcco3,"%e\n",y[N2+EQHCCO3][j]);
        fprintf(fpdcco3,"%e\n", y[N2+EQCCO3][j]);
#endif
#ifdef OXYGEN
        fprintf(fpdo2,"%e\n",  y[N2+EQO2][j]);
#endif
#ifdef BORON
        fprintf(fpdboh3,"%e\n",  y[N2+EQBOH3][j]);
        fprintf(fpdboh4,"%e\n",  y[N2+EQBOH4][j]);
#ifdef BORISTP
        fprintf(fpdbboh3,"%e\n",  y[N2+EQBBOH3][j]);
        fprintf(fpdbboh4,"%e\n",  y[N2+EQBBOH4][j]);
#endif
#endif
#ifdef CALCIUM
        fprintf(fpdca,"%e\n",  y[N2+EQCA][j]);
#endif
      } /* for */

      fclose(fpr);

      fclose(fpco2);
      fclose(fphco3);
      fclose(fpco3);
      fclose(fph);
      fclose(fpoh);
#ifdef C13ISTP
      fclose(fpcco2);
      fclose(fphcco3);
      fclose(fpcco3);
#endif
#ifdef OXYGEN
      fclose(fpo2);
#endif
#ifdef BORON
      fclose(fpboh3);
      fclose(fpboh4);
#ifdef BORISTP
      fclose(fpbboh3);
      fclose(fpbboh4);
#endif
#endif
#ifdef CALCIUM
      fclose(fpca);
#endif

      fclose(fpdco2);
      fclose(fpdhco3);
      fclose(fpdco3);
      fclose(fpdh);
      fclose(fpdoh);
#ifdef C13ISTP
      fclose(fpdcco2);
      fclose(fpdhcco3);
      fclose(fpdcco3);
#endif
#ifdef OXYGEN
      fclose(fpdo2);
#endif
#ifdef BORON
      fclose(fpdboh3);
      fclose(fpdboh4);
#ifdef BORISTP
      fclose(fpdbboh3);
      fclose(fpdbboh4);
#endif
#endif
#ifdef CALCIUM
      fclose(fpdca);
#endif
	nrerror("Too many iterations in SOLVDE");
}

void bksub(ne,nb,jf,k1,k2,c)
int ne,nb,jf,k1,k2;
double ***c;
{
	int nbf,im,kp,k,j,i;
	double xx;

	nbf=ne-nb;
	im=1;
	for (k=k2;k>=k1;k--) {
		if (k == k1) im=nbf+1;
		kp=k+1;
		for (j=1;j<=nbf;j++) {
			xx=c[j][jf][kp];
			for (i=im;i<=ne;i++)
				c[i][jf][k] -= c[i][j][k]*xx;
		}
	}
	for (k=k1;k<=k2;k++) {
		kp=k+1;
		for (i=1;i<=nb;i++) c[i][1][k]=c[i+nbf][jf][k];
		for (i=1;i<=nbf;i++) c[i+nb][1][k]=c[i][jf][kp];
	}
}


void pinvs(ie1,ie2,je1,jsf,jc1,k,c,s)
int ie1,ie2,je1,jsf,jc1,k;
double ***c,**s;
{
	int js1,jpiv,jp,je2,jcoff,j,irow,ipiv,id,icoff,i,*indxr,*ivector();
	double pivinv,piv,dum,big,*pscl,*dvector();
	void nrerror(),free_dvector(),free_ivector();

	indxr=ivector(ie1,ie2);
	pscl=dvector(ie1,ie2);
	je2=je1+ie2-ie1;
	js1=je2+1;
	for (i=ie1;i<=ie2;i++) {
		big=0.0;
		for (j=je1;j<=je2;j++)
			if (fabs(s[i][j]) > big) big=fabs(s[i][j]);
		if (big == 0.0) nrerror("Singular matrix - row all 0, in PINVS");
		pscl[i]=1.0/big;
		indxr[i]=0;
	}
	for (id=ie1;id<=ie2;id++) {
		piv=0.0;
		for (i=ie1;i<=ie2;i++) {
			if (indxr[i] == 0) {
				big=0.0;
				for (j=je1;j<=je2;j++) {
					if (fabs(s[i][j]) > big) {
						jp=j;
						big=fabs(s[i][j]);
					}
				}
				if (big*pscl[i] > piv) {
					ipiv=i;
					jpiv=jp;
					piv=big*pscl[i];
				}
			}
		}
		if (s[ipiv][jpiv] == 0.0) nrerror("Singular matrix in routine PINVS");
		indxr[ipiv]=jpiv;
		pivinv=1.0/s[ipiv][jpiv];
		for (j=je1;j<=jsf;j++) s[ipiv][j] *= pivinv;
		s[ipiv][jpiv]=1.0;
		for (i=ie1;i<=ie2;i++) {
			if (indxr[i] != jpiv) {
				if (s[i][jpiv]) {
					dum=s[i][jpiv];
					for (j=je1;j<=jsf;j++)
						s[i][j] -= dum*s[ipiv][j];
					s[i][jpiv]=0.0;
				}
			}
		}
	}
	jcoff=jc1-js1;
	icoff=ie1-je1;
	for (i=ie1;i<=ie2;i++) {
		irow=indxr[i]+icoff;
		for (j=js1;j<=jsf;j++) c[irow][j+jcoff][k]=s[i][j];
	}
	free_dvector(pscl,ie1,ie2);
	free_ivector(indxr,ie1,ie2);
}

void red(iz1,iz2,jz1,jz2,jm1,jm2,jmf,ic1,jc1,jcf,kc,c,s)
double ***c,**s;
int iz1,iz2,jz1,jz2,jm1,jm2,jmf,ic1,jc1,jcf,kc;
{
	int loff,l,j,ic,i;
	double vx;

	loff=jc1-jm1;
	ic=ic1;
	for (j=jz1;j<=jz2;j++) {
		for (l=jm1;l<=jm2;l++) {
			vx=c[ic][l+loff][kc];
			for (i=iz1;i<=iz2;i++) s[i][l] -= s[i][j]*vx;
		}
		vx=c[ic][jcf][kc];
		for (i=iz1;i<=iz2;i++) s[i][jmf] -= s[i][j]*vx;
		ic += 1;
	}
}

/* =========================================================
   =========================================================

             Numerical Recipes   (end)

   =========================================================
   ========================================================= */



void difeq(k,k1,k2,jsf,is1,isf,indexv,ne,s,y)
   int k,k1,k2,jsf,is1,isf,indexv[],ne;
   double **s,**y;
/* --- returns matrix s for solvde */
{

   int a,b;


#ifdef DEB07
   printf("--- difeq: in --- \n");
#endif

#ifdef NOK5NO
     /* ----- diagnostic ----- */
      for (a = 0; a <= M; a++)
         co3[a] = dicbulk / (hplus[a] * hplus[a] / k1d / k2d
                           + hplus[a] / k2d + 1.0);
#endif

/*
                --- pure diffusion (g = 0) ---


                --- diffusion-reaction ---

       spherical symmetric:

       c_j'' + 2/r c_j' = - reaction_j / D_c_j = g_j / D_c_j

       y_j = c_j                              j = 1,...,N2

       y_j'      = y_{j+N2}                   j = 1,...,N2
       y_{j+N2}' = g_j - 2/r y_{j+N2}         j = 1,...,N2
       y_{j+N2}' + reaction + 2/r y_{j+N2}         j = 1,...,N2
       j = 1,...,N2
       E_j,k    = y_j,k    - y_j,k-1    - h/2 [y_j+N2,k  + y_j+N2,k-1]  = 0
       E_j+N2,k = y_j+N2,k - y_j+N2,k-1 - h/2 [g_j (y_k) + g_j (y_k-1)]
                  + h [y_k / r_k + y_k-1 / r_k] = 0

       S_j,n    = d E_j,k / d y_n,k-1     j = 1,..., NE; n = 1, ..., NE
       S_j,n+NE = d E_j,k / d y_n,k       j = 1,..., NE; n = 1, ..., NE

       S_1,1 = -1  \
       S_2,2 = -1   \  S_j,n = - delta_j,n      1 <= j,n <= N2
       S_1,2 =  0   /
       S_2,1 =  0  /

       S_1,3 = -h/2  \
       S_2,4 = -h/2   \  S_j,n+N2 = - h/2 * delta_j,n   1 <= j,n <= N2
       S_1,4 = 0      /
       S_2,3 = 0     /

       S_1,5 =  1  \
       S_2,6 =  1   \  S_j,n+2*N2 = delta_j,n      1 <= j,n <= N2
       S_1,6 =  0   /
       S_2,5 =  0  /

       S_1,7 = -h/2  \
       S_2,8 = -h/2   \  S_j,n+3*N2 = - h/2 * delta_j,n   1 <= j,n <= N2
       S_1,8 = 0      /
       S_2,7 = 0     /

       S_3,1 = -h/2 dg_1/dy_1,k-1  \
       S_4,2 = -h/2 dg_2/dy_2,k-1   \  S_j+N2,n = - h/2 dg_j/dy_n,k-1
       S_3,2 = -h/2 dg_1/dy_2,k-1   /        1 <= j,n <= N2
       S_4,1 = -h/2 dg_2/dy_1,k-1  /

       S_3,3 = -1  \
       S_4,4 = -1   \  S_j,n+N2 = - h/2 * delta_j,n   1 <= j,n <= N2
       S_4,4 =  0   /
       S_3,3 =  0  /

       S_3,5 = -h/2 dg_1/dy_1,k \
       S_4,6 = -h/2 dg_2/dy_2,k  \  S_j,n+2*N2 = -h/2 dg_j/dy_n,k
       S_3,6 = -h/2 dg_1/dy_2,k  /        1 <= j,n <= N2
       S_4,5 = -h/2 dg_2/dy_1,k /

       S_3,7 = 1  \
       S_4,8 = 1   \  S_j,n+3*N2 = delta_j,n   1 <= j,n <= N2
       S_3,8 = 0   /
       S_4,7 = 0  /

       -------------------------------------------

               boundary conditions

       NRB  = number of right boundary conditions
       NLB  = number of left boundary conditions

       left (x_1):

       E_LB_j = E_1
              = (E_1,1, E_2,1, E_3,1,        E_4,1)
              = (0,     0,     y_1,1 - C1X1, y_2,1 - C2X1)

       S-LB_j,n+NE = d E_j,1 / d y_n,1       j=NRB+1,..., NE; n = 1,...,NE

       S-LB_3,5 = +1
       S-LB_4,5 =  0

       S-LB_3,6 =  0
       S-LB_4,6 = +1

       S-LB_3,7 =  0
       S-LB_4,7 =  0

       S-LB_3,8 =  0
       S-LB_4,8 =  0


       right (x_2):

       E_RB = E_M+1
            = (E_1,M+1,        E_2,M+1,        0, 0)
            = (y_1,M+1 - C1X2, y_2,M+1 - C2X2, 0, 0)


       S-RB_j,n = d E_j,M+1 / d y_n,M+1       j=1,..., NRB; n = 1,...,NE

       S-RB_1,1 = +1
       S-RB_2,1 =  0

       S-RB_1,2 =  0
       S-RB_2,2 = +1

       S-RB_1,3 =  0
       S-RB_2,3 =  0

       S-RB_1,4 =  0
       S-RB_2,4 =  0


*/

   if(k == k1) {

   /* --- boundary conditions at left boundary --- */


      for(a=1; a <= N2; a++)
      for(b=1; b <= N2; b++) {
         s[NRB+a][NE+indexv[b]]    = 0.0;
         s[NRB+a][NE+indexv[N2+b]] = 0.0;
      }

      for(a=1; a <= N2; a++) s[NRB+a][NE+indexv[N2+a]] = 1.0;



#ifdef MIMECO2DIA

               /* Michaelis-Menten at the shell (diatoms) */

   printf("----- Michaelis-Menten at the shell (diatoms) ----- \n");
   printf("%e   y[EQCO2][1] \n",y[EQCO2][1]);
   printf("%e   VMAXDIA \n",VMAXDIA);
   printf("%e   KSDIA \n",KSDIA);
   printf("%e   surface \n",surface);
   printf("%e   dco2 \n",dco2);
   dummy = VMAXDIA*y[EQCO2][1]/(KSDIA+y[EQCO2][1]);/*/surface/dco2;*/
   printf("%e   dummy \n",dummy);
   s[NRB+1][jsf] = y[N2+EQCO2][1] -  dummy;
   printf("%e  s[NRB+1][jsf]  \n",s[NRB+1][jsf]);
   s[NRB+1][NE+indexv[1]] =
   VMAXDIA * KSDIA / (KSDIA + y[EQCO2][1]) / (KSDIA + y[EQCO2][1]);
   printf("%e  s[NRB+1][NE+indexv[1]]  \n",s[NRB+1][NE+indexv[1]]);
#else
      s[NRB+1][jsf] = y[N2+EQCO2][1]  -  co2flux;
#endif
      s[NRB+2][jsf] = y[N2+EQHCO3][1] - hco3flux;
      s[NRB+3][jsf] = y[N2+EQCO3][1]  -  co3flux;
      s[NRB+4][jsf] = y[N2+EQHP][1]   -    hflux;
      s[NRB+5][jsf] = y[N2+EQOH][1]   -   ohflux;
#ifdef C13ISTP





#ifdef F13_CO3
      s[NRB+EQCCO2][jsf]  = y[N2+EQCCO2][1]  -  cco2flux;
      s[NRB+EQHCCO3][jsf] = y[N2+EQHCCO3][1] - hcco3flux;

      /* Here is the left boundary condition for the		*/
      /* 13CO3-- flux which gives the final d13C of the shell!!	*/
      /* The ratio of 13Fcalc and 12Fcalc is given by:		*/

      /* 13Fc/12Fc = alphac * 13CO3(R1) / 12CO3(R1)		*/

      /* where alphac is the eq. fractionation factor between 	*/
      /* CO3-- and CaCO3. It is calcul. by a combination of eps-*/
      /* values for HCO3-CO3 and HCO3-CaCO3, see Mook (1986).   */
      /* Solve for 13Fc, note that Fc = 12Fc + 13Fc.		*/
      /*

      => cco3upt = CO3UPT/( (CO3(R1)-13CO3(R1))/alphac/13CO3(R1) + 1 )
								*/
#ifdef CISTP
      /*
      => cco3upt = CO3UPT/( CO3(R1)/alphac/13CO3(R1) + 1 )
								*/
#endif


#ifdef PRINT
      if(CO3UPT != 0.0){
        tmp = (cco3upt/(CO3UPT-cco3upt)/RSTAND - 1.)*1000.;
        printf("\n%e d13F_calc\n",tmp);
      }
#endif

      dummy = (y[EQCO3][1]-y[EQCCO3][1])/alphac/y[EQCCO3][1] + 1.;
#ifdef CISTP
      dummy = (y[EQCO3][1]             )/alphac/y[EQCCO3][1] + 1.;
#endif
      cco3upt  = CO3UPT / dummy;
#ifdef LDAT
      cco3upt  = co3uptldat / dummy;
#endif
      cco3flux = cco3upt/surface/dcco3*1.e21;

#ifdef PRINT
      if(CO3UPT != 0.0){
	tmp = (cco3upt/(CO3UPT-cco3upt)/RSTAND - 1.)*1000.;
      	printf("%e d13F_calc\n",tmp);
#ifndef CISTP
	printf("%e resid\n",cco3upt/(CO3UPT-cco3upt)
		- alphac*y[EQCCO3][1]/(y[EQCO3][1]-y[EQCCO3][1]));
#endif
#ifdef CISTP
	printf("%e resid\n",cco3upt/(CO3UPT-cco3upt)
		- alphac*y[EQCCO3][1]/(y[EQCO3][1]));
#endif
      }
#endif /* PRINT */

      s[NRB+EQCCO3][jsf]   = y[N2+EQCCO3][1]  -  cco3flux;
      s[NRB+EQCCO3][NE+indexv[EQCO3]] =
	   CO3UPT*1.e21/surface/dcco3/SQ(dummy)/alphac/y[EQCCO3][1];
      s[NRB+EQCCO3][NE+indexv[EQCCO3]] =
	  -CO3UPT*1.e21*y[EQCO3][1]/surface/dcco3;
      s[NRB+EQCCO3][NE+indexv[EQCCO3]] /=
	   SQ(dummy)*alphac*SQ(y[EQCCO3][1]);
#ifdef LDAT
      s[NRB+EQCCO3][jsf]   = y[N2+EQCCO3][1]  -  cco3flux;
      s[NRB+EQCCO3][NE+indexv[EQCO3]] =
	   co3uptldat*1.e21/surface/dcco3/SQ(dummy)/alphac/y[EQCCO3][1];
      s[NRB+EQCCO3][NE+indexv[EQCCO3]] =
	  -co3uptldat*1.e21*y[EQCO3][1]/surface/dcco3;
      s[NRB+EQCCO3][NE+indexv[EQCCO3]] /=
	   SQ(dummy)*alphac*SQ(y[EQCCO3][1]);
#endif
#ifdef CISTP
      /* set new left boundary cond. for 12C */
      s[NRB+3][jsf] = y[N2+EQCO3][1]  -
      		(CO3UPT/surface/dco3*1.e21  - cco3flux);
#ifdef LDAT
      /* set new left boundary cond. for 12C */
      s[NRB+3][jsf] = y[N2+EQCO3][1]  -
      		(co3uptldat/surface/dco3*1.e21  - cco3flux);
#endif
#endif

#else /* F13_CO3 */
      s[NRB+EQCCO3][jsf]   = y[N2+EQCCO3][1]  -  cco3flux;
      if(CO3UPT != 0.0){
     	tmp1 = (cco3upt/(CO3UPT-cco3upt)/RSTAND - 1.)*1000;
     	tmp2 = (cco3[1]/(co3[1]-cco3[1])/RSTAND - 1.)*1000;
     	printf("%e d13F_calc\n",tmp1);
     	printf("%e d13CO3(R1)\n\n",tmp2);
      }
#ifdef CISTP
      /* set new left boundary cond. for 12C */
      s[NRB+3][jsf] = y[N2+EQCO3][1]  -
      		(CO3UPT/surface/dco3*1.e21  - cco3flux);
#ifdef LDAT
      /* set new left boundary cond. for 12C */
      s[NRB+3][jsf] = y[N2+EQCO3][1]  -
      		(co3uptldat/surface/dco3*1.e21  - cco3flux);
#endif
#endif
#endif /* F13_CO3 */



#ifdef F13_HCO3
      s[NRB+EQCCO2][jsf]  = y[N2+EQCCO2][1]  - cco2flux;
      s[NRB+EQCCO3][jsf]  = y[N2+EQCCO3][1]  - cco3flux;

      /* Here is the left boundary condition for the		*/
      /* 13HCO3- flux which gives the final d13C of the shell!!	*/
      /* The ratio of 13Fcalc and 12Fcalc is given by:		*/

      /* 13Fc/12Fc = alpha * 13HCO3(R1) / 12HCO3(R1)		*/

      /* where alpha is the eq. fractionation factor between 	*/
      /* HCO3- and CaCO3, Mook (1986).   			*/
      /* Solve for 13Fc, note that Fc = 12Fc + 13Fc.		*/
      /*

      => hcco3upt = HCO3UPT/( (HCO3(R1)-H13CO3(R1))/alphahc/H13CO3(R1) + 1 )
								*/
#ifdef CISTP
      /*
      => hcco3upt = HCO3UPT/( HCO3(R1)/alphahc/H13CO3(R1) + 1 )
								*/
#endif

#ifdef PRINT
      if(HCO3UPT != 0.0){
        tmp = (hcco3upt/(HCO3UPT-hcco3upt)/RSTAND - 1.)*1000.;
        printf("\n%e d13F_calc\n",tmp);
      }
#endif

      dummy = (y[EQHCO3][1]-y[EQHCCO3][1])/alphahc/y[EQHCCO3][1] + 1.;
#ifdef CISTP
      dummy = (y[EQHCO3][1]             )/alphahc/y[EQHCCO3][1] + 1.;
#endif
      hcco3upt  = HCO3UPT / dummy;
#ifdef LDAT
      hcco3upt  = hco3uptldat / dummy;
#endif
      hcco3flux = hcco3upt/surface/dhcco3*1.e21;

#ifdef PRINT
      if(HCO3UPT != 0.0){
	tmp = (hcco3upt/(HCO3UPT-hcco3upt)/RSTAND - 1.)*1000.;
      	printf("%e d13F_calc\n",tmp);
#ifndef CISTP
	printf("%e resid\n",hcco3upt/(HCO3UPT-hcco3upt)
		- alphahc*y[EQHCCO3][1]/(y[EQHCO3][1]-y[EQHCCO3][1]));
#endif
#ifdef CISTP
	printf("%e resid\n",hcco3upt/(HCO3UPT-hcco3upt)
		- alphahc*y[EQHCCO3][1]/(y[EQHCO3][1]));
#endif
      }
#endif /* PRINT */

      s[NRB+EQHCCO3][jsf]   = y[N2+EQHCCO3][1]  -  hcco3flux;
      s[NRB+EQHCCO3][NE+indexv[EQHCO3]] =
	   HCO3UPT*1.e21/surface/dhcco3/SQ(dummy)/alphahc/y[EQHCCO3][1];
      s[NRB+EQHCCO3][NE+indexv[EQHCCO3]] =
	  -HCO3UPT*1.e21*y[EQHCO3][1]/surface/dhcco3;
      s[NRB+EQHCCO3][NE+indexv[EQHCCO3]] /=
	   SQ(dummy)*alphahc*SQ(y[EQHCCO3][1]);
#ifdef LDAT
      s[NRB+EQHCCO3][jsf]   = y[N2+EQHCCO3][1]  -  hcco3flux;
      s[NRB+EQHCCO3][NE+indexv[EQHCO3]] =
	   hco3uptldat*1.e21/surface/dhcco3/SQ(dummy)/alphahc/y[EQHCCO3][1];
      s[NRB+EQHCCO3][NE+indexv[EQHCCO3]] =
	  -hco3uptldat*1.e21*y[EQHCO3][1]/surface/dhcco3;
      s[NRB+EQHCCO3][NE+indexv[EQHCCO3]] /=
	   SQ(dummy)*alphahc*SQ(y[EQHCCO3][1]);
#endif
#ifdef CISTP
      /* set new left boundary cond. for 12C */
      s[NRB+EQHCO3][jsf] = y[N2+EQHCO3][1]  -
      		(HCO3UPT/surface/dhco3*1.e21  - hcco3flux);
#ifdef LDAT
      /* set new left boundary cond. for 12C */
      s[NRB+EQHCO3][jsf] = y[N2+EQHCO3][1]  -
      		(hco3uptldat/surface/dhco3*1.e21  - hcco3flux);
#endif
#endif

#else /* F13_HCO3 */
      s[NRB+EQHCCO3][jsf]   = y[N2+EQHCCO3][1]  -  hcco3flux;
      if(HCO3UPT != 0.0){
     	tmp1 = (hcco3upt/(HCO3UPT-hcco3upt)/RSTAND - 1.)*1000;
     	tmp2 = (hcco3[1]/(hco3[1]-hcco3[1])/RSTAND - 1.)*1000;
     	printf("%e d13F_calc\n",tmp1);
     	printf("%e d13HCO3(R1)\n\n",tmp2);
      }
#ifdef CISTP
      /* set new left boundary cond. for 12C */
      s[NRB+EQHCO3][jsf] = y[N2+EQHCO3][1]  -
      		(HCO3UPT/surface/dhco3*1.e21  - hcco3flux);
#ifdef LDAT
      /* set new left boundary cond. for 12C */
      s[NRB+EQHCO3][jsf] = y[N2+EQHCO3][1]  -
      		(hco3uptldat/surface/dhco3*1.e21  - hcco3flux);
#endif
#endif
#endif /* F13_HCO3 */

#endif /* C13ISTP */


#ifdef BORON
      s[NRB+EQBOH3][jsf] = y[N2+EQBOH3][1]   -   boh3flux;
      s[NRB+EQBOH4][jsf] = y[N2+EQBOH4][1]   -   boh4flux;
#ifdef BORISTP
      s[NRB+EQBBOH3][jsf] = y[N2+EQBBOH3][1]   -   bboh3flux;
      s[NRB+EQBBOH4][jsf] = y[N2+EQBBOH4][1]   -   bboh4flux;
#endif
#endif
#ifdef OXYGEN
      s[NRB+EQO2][jsf] = y[N2+EQO2][1]   -   o2flux;
#endif
#ifdef CALCIUM
      s[NRB+EQCA][jsf] = y[N2+EQCA][1]   -   caflux;
#endif



   } else if(k > k2) {

   /* --- boundary conditions at right boundary --- */



      for(a=1; a <= N2; a++)
      for(b=1; b <= N2; b++) {
         s[a][NE+indexv[b]]    = 0.0;
         s[a][NE+indexv[N2+b]] = 0.0;
      }

      for(a=1; a <= N2; a++)  s[a][NE+indexv[a]] = 1.0; /* prior 28.11.93 */

      s[1][jsf] = y[EQCO2][M]  -  co2bulk;
      s[2][jsf] = y[EQHCO3][M] - hco3bulk;
      s[3][jsf] = y[EQCO3][M]  -  co3bulk;
      s[4][jsf] = y[EQHP][M]   -    hbulk;
      s[5][jsf] = y[EQOH][M]   -   ohbulk;
#ifdef C13ISTP
      s[EQCCO2][jsf]  = y[EQCCO2][M]   -  cco2bulk;
      s[EQHCCO3][jsf] = y[EQHCCO3][M]  - hcco3bulk;
      s[EQCCO3][jsf]  = y[EQCCO3][M]   -  cco3bulk;
#endif
#ifdef BORON
      s[EQBOH3][jsf] = y[EQBOH3][M]   -   boh3bulk;
      s[EQBOH4][jsf] = y[EQBOH4][M]   -   boh4bulk;
#ifdef BORISTP
      s[EQBBOH3][jsf] = y[EQBBOH3][M]   -   bboh3bulk;
      s[EQBBOH4][jsf] = y[EQBBOH4][M]   -   bboh4bulk;
#endif
#endif
#ifdef OXYGEN
      s[EQO2][jsf] = y[EQO2][M]   -   o2bulk;
#endif
#ifdef CALCIUM
      s[EQCA][jsf] = y[EQCA][M]   -   cabulk;
#endif

   } else {



   /* --- interior point --- */

/* first index: equation;
  second index: dependent variable */

/*
       S_1,1 = -1  \
       S_2,2 = -1   \  S_j,n = - delta_j,n      1 <= j,n <= N2
       S_1,2 =  0   /
       S_2,1 =  0  /

       S_1,3 = -h/2  \
       S_2,4 = -h/2   \  S_j,n+N2 = - h/2 * delta_j,n   1 <= j,n <= N2
       S_1,4 = 0      /
       S_2,3 = 0     /

       S_1,5 =  1  \
       S_2,6 =  1   \  S_j,n+2*N2 = delta_j,n      1 <= j,n <= N2
       S_1,6 =  0   /
       S_2,5 =  0  /

       S_1,7 = -h/2  \
       S_2,8 = -h/2   \  S_j,n+3*N2 = - h/2 * delta_j,n   1 <= j,n <= N2
       S_1,8 = 0      /
       S_2,7 = 0     /

       S_3,1 = -h/2 dg_1/dy_1,k-1  \
       S_4,2 = -h/2 dg_2/dy_2,k-1   \  S_j+N2,n = - h/2 dg_j/dy_n,k-1
       S_3,2 = -h/2 dg_1/dy_2,k-1   /        1 <= j,n <= N2
       S_4,1 = -h/2 dg_2/dy_1,k-1  /

       S_3,3 = -1  \
       S_4,4 = -1   \  S_j,n+N2 = - h/2 * delta_j,n   1 <= j,n <= N2
       S_4,4 =  0   /
       S_3,3 =  0  /

       S_3,5 = -h/2 dg_1/dy_1,k \
       S_4,6 = -h/2 dg_2/dy_2,k  \  S_j,n+2*N2 = -h/2 dg_j/dy_n,k
       S_3,6 = -h/2 dg_1/dy_2,k  /        1 <= j,n <= N2
       S_4,5 = -h/2 dg_2/dy_1,k /

       S_3,7 = 1  \
       S_4,8 = 1   \  S_j,n+3*N2 = delta_j,n   1 <= j,n <= N2
       S_3,8 = 0   /
       S_4,7 = 0  /

*/

      for(a=1; a <= N2; a++)
      for(b=1; b <= N2; b++) {
        s[   a][   indexv[   b]] = 0.0;
        s[   a][   indexv[N2+b]] = 0.0;
        s[   a][NE+indexv[   b]] = 0.0;
        s[   a][NE+indexv[N2+b]] = 0.0;
        s[N2+a][   indexv[N2+b]] = 0.0;
        s[N2+a][NE+indexv[N2+b]] = 0.0;
      }
#ifdef CLPL
      for(a=1; a <= N2; a++) {
        s[   a][   indexv[   a]] = -1.0;
        s[   a][   indexv[N2+a]] = -hh;
        s[   a][NE+indexv[   a]] =  1.0;
        s[   a][NE+indexv[N2+a]] = -hh;
        s[N2+a][   indexv[N2+a]] = -1.0 +
        (r[k]*r[k]*difcofm[a][k]-r[k-1]*r[k-1]*difcofm[a][k-1])
        /2./difcofm[a][k-1]/r[k-1]/r[k-1];
        s[N2+a][NE+indexv[N2+a]] =  1.0 +
        (r[k]*r[k]*difcofm[a][k]-r[k-1]*r[k-1]*difcofm[a][k-1])/
        2./difcofm[a][k]/r[k]/r[k];
      }
#else
      for(a=1; a <= N2; a++) {
        s[   a][   indexv[   a]] = -1.0;
        s[   a][   indexv[N2+a]] = -hh;
        s[   a][NE+indexv[   a]] =  1.0;
        s[   a][NE+indexv[N2+a]] = -hh;
        s[N2+a][   indexv[N2+a]] = -1.0 + h / r[k-1];
        s[N2+a][NE+indexv[N2+a]] =  1.0 + h / r[k];
      }
#endif

      /* --- all to zero --- */

      for(a=1; a <= N2; a++)
      for(b=1; b <= N2; b++) {
        s[N2+a][   indexv[b]] = 0.0;
        s[N2+a][NE+indexv[b]] = 0.0;
      }

#ifdef REACTION

#ifdef CLPL
     if(k < dumk1 || k >= dumk2){
#endif

#ifndef NOCO2

  /* ===== CO2 ================= */

  /* -----   dCO2 / dCO2   ----- */

  s[N2+EQCO2][indexv[1]]    =  -hh * (kp1s + kp4 * oh[k-1]) / dco2;
  s[N2+EQCO2][NE+indexv[1]] =  -hh * (kp1s + kp4 * oh[k]  ) / dco2;

  /* -----   dCO2 / dHCO3   ----- */

  /* kh2co3/Kh2co3 = km1s */

  s[N2+EQCO2][indexv[2]] =
     hh * (km1s * hplus[k-1] + km4) / dco2;
  s[N2+EQCO2][NE+indexv[2]] =
     hh * (km1s * hplus[k]   + km4) / dco2;

  /* -----   dCO2 / dCO3   ----- */

  s[N2+EQCO2][indexv[3]]    = 0.0;
  s[N2+EQCO2][NE+indexv[3]] = 0.0;

  /* -----   dCO2 / dH   ----- */

  s[N2+EQCO2][indexv[4]]    = hh * km1s * hco3[k-1] / dco2;
  s[N2+EQCO2][NE+indexv[4]] = hh * km1s * hco3[k]   / dco2;

  /* -----   dCO2 / dOH   ----- */

  s[N2+EQCO2][indexv[5]]    = -hh * kp4 * co2[k-1] / dco2;
  s[N2+EQCO2][NE+indexv[5]] = -hh * kp4 * co2[k]   / dco2;

#endif

#ifndef NOHCO3


  /* ===== HCO3 ================= */

  /* -----   dHCO3 / dCO2   ----- */

  s[N2+EQHCO3][indexv[1]]    =
     hh * (kp1s + kp4 * oh[k-1]) / dhco3;
  s[N2+EQHCO3][NE+indexv[1]] =
     hh * (kp1s + kp4 * oh[k])   / dhco3;

  /* -----   dHCO3 / dHCO3   ----- */

  s[N2+EQHCO3][indexv[2]] =
     - hh * (km1s*hplus[k-1] + km4 + km5h) / dhco3;
  s[N2+EQHCO3][NE+indexv[2]] =
     - hh * (km1s*hplus[k]   + km4 + km5h) / dhco3;

  /* -----   dHCO3 / dCO3   ----- */

  s[N2+EQHCO3][indexv[3]]    = hh * kp5h * hplus[k-1] / dhco3;
  s[N2+EQHCO3][NE+indexv[3]] = hh * kp5h * hplus[k]   / dhco3;

  /* -----   dHCO3 / dH     ----- */

  s[N2+EQHCO3][indexv[4]]    =
   hh * (kp5h * co3[k-1] - km1s*hco3[k-1]) / dhco3;
  s[N2+EQHCO3][NE+indexv[4]] =
   hh * (kp5h * co3[k]   - km1s*hco3[k])   / dhco3;

  /* -----   dHCO3 / dOH    ----- */

  s[N2+EQHCO3][indexv[5]]    =
      hh * kp4 * co2[k-1] / dhco3;
  s[N2+EQHCO3][NE+indexv[5]] =
      hh * kp4 * co2[k]   / dhco3;

#endif

#ifndef NOCO3

  /* ===== CO3 ================= */


  /* -----   dCO3 / dCO2   ----- */

  s[N2+EQCO3][indexv[1]]    = 0.0;
  s[N2+EQCO3][NE+indexv[1]] = 0.0;

  /* -----   dCO3 / dHCO3   ----- */

  s[N2+EQCO3][indexv[2]]    = hh * km5h / dco3;
  s[N2+EQCO3][NE+indexv[2]] = hh * km5h / dco3;

  /* -----   dCO3 / dCO3   ----- */

  s[N2+EQCO3][indexv[3]]    = - hh * kp5h * hplus[k-1] / dco3;
  s[N2+EQCO3][NE+indexv[3]] = - hh * kp5h * hplus[k]   / dco3;

  /* -----   dCO3 / dH   ----- */

  s[N2+EQCO3][indexv[4]]    = - hh * kp5h * co3[k-1] / dco3;
  s[N2+EQCO3][NE+indexv[4]] = - hh * kp5h * co3[k]   / dco3;

  /* -----   dCO3 / dOH   ----- */

  s[N2+EQCO3][indexv[5]]    = 0.0;
  s[N2+EQCO3][NE+indexv[5]] = 0.0;

#ifdef CARTEST
if(k == 20){
	printf("%e hco3[k]\n",hco3[k]);
	printf("%e hco3bulk\n",hco3bulk);
	printf("%e hplus[k]\n",hplus[k]);
	printf("%e hbulk\n",hbulk);
	printf("%e km5h*hco3[k]\n",km5h*hco3[k]);
	printf("%e kp5h*hplus[k]*co3[k]\n",kp5h*hplus[k]*co3[k]);
	printf("%e kp6\n",kp6);
	printf("%e km6*hplus[k]*oh[k]\n",km6*hplus[k]*oh[k]);
}
#endif

#endif

#ifndef NOHPLUS

  /* ===== H ================= */

#ifndef CARTEST

  /* -----   dH / dCO2   ----- */

  s[N2+EQHP][indexv[1]]    = hh * kp1s / dh;
  s[N2+EQHP][NE+indexv[1]] = hh * kp1s / dh;

  /* -----   dH / dHCO3   ----- */

  s[N2+EQHP][indexv[2]] =
    hh * (km5h - km1s*hplus[k-1]) / dh;
  s[N2+EQHP][NE+indexv[2]] =
    hh * (km5h - km1s*hplus[k])   / dh;

#endif

  /* -----   dH / dCO3   ----- */

  s[N2+EQHP][indexv[3]]    = - hh * kp5h * hplus[k-1] / dh;
  s[N2+EQHP][NE+indexv[3]] = - hh * kp5h * hplus[k]   / dh;


  /* -----   dH / dH   ----- */

  s[N2+EQHP][indexv[4]] =
    - hh * (
#ifndef CARTEST
    km1s * hco3[k-1]
#endif
    + kp5h * co3[k-1] + km6 * oh[k-1]
#ifdef BORONRC3
    + km7 * boh4[k-1]
#ifdef BORISTP
#ifdef B10B11
    + km7bb * bboh4[k-1]
#endif
#endif
#endif
#ifdef CISTP
    + km1scc * hcco3[k-1] + kp5cc * cco3[k-1]
#endif
   ) / dh;

  s[N2+EQHP][NE+indexv[4]] =
    - hh * (
#ifndef CARTEST
    km1s * hco3[k]
#endif
         + kp5h * co3[k] + km6 * oh[k]
#ifdef BORONRC3
    + km7 * boh4[k]
#ifdef BORISTP
#ifdef B10B11
    + km7bb * bboh4[k]
#endif
#endif
#endif
#ifdef CISTP
    + km1scc * hcco3[k] + kp5cc * cco3[k]
#endif
    ) / dh;

  /* -----   dH / dOH   ----- */

  s[N2+EQHP][indexv[5]]    = - hh * km6 * hplus[k-1] / dh;
  s[N2+EQHP][NE+indexv[5]] = - hh * km6 * hplus[k]   / dh;


#ifdef BORONRC3
  /* -----   dH / dB(OH)3   ----- */

  s[N2+EQHP][indexv[EQBOH3]] =
    hh * kp7 / dh;
  s[N2+EQHP][NE+indexv[EQBOH3]] =
    hh * kp7 / dh;

  /* -----   dH / dB(OH)4   ----- */

  s[N2+EQHP][indexv[EQBOH4]] =
    -hh * km7*hplus[k-1] / dh;
  s[N2+EQHP][NE+indexv[EQBOH4]] =
    -hh * km7*hplus[k]   / dh;
#ifdef BORISTP
#ifdef B10B11
  /* -----   dH / d11B(OH)3   ----- */

  s[N2+EQHP][indexv[EQBBOH3]] =
    hh * kp7bb / dh;
  s[N2+EQHP][NE+indexv[EQBBOH3]] =
    hh * kp7bb / dh;

  /* -----   dH / d11B(OH)4   ----- */

  s[N2+EQHP][indexv[EQBBOH4]] =
    -hh * km7bb*hplus[k-1] / dh;
  s[N2+EQHP][NE+indexv[EQBBOH4]] =
    -hh * km7bb*hplus[k]   / dh;
#endif
#endif
#endif


#ifdef CISTP

  /* -----   dH / d13CO2   ----- */

  s[N2+EQHP][indexv[EQCCO2]]    = hh * kp1scc / dh;
  s[N2+EQHP][NE+indexv[EQCCO2]] = hh * kp1scc / dh;

  /* -----   dH / dH13CO3   ----- */

  s[N2+EQHP][indexv[EQHCCO3]] =
    hh * (km5cc - km1scc*hplus[k-1]) / dh;
  s[N2+EQHP][NE+indexv[EQHCCO3]] =
    hh * (km5cc - km1scc*hplus[k])   / dh;

  /* -----   dH / d13CO3   ----- */

  s[N2+EQHP][indexv[EQCCO3]]    = - hh * kp5cc * hplus[k-1] / dh;
  s[N2+EQHP][NE+indexv[EQCCO3]] = - hh * kp5cc * hplus[k]   / dh;
#endif

#endif

#ifndef NOOH

  /* ===== OH ================= */

#ifndef CARTEST

  /* -----   dOH / dCO2   ----- */

  s[N2+EQOH][indexv[1]]    = - hh * kp4 * oh[k-1] / doh;
  s[N2+EQOH][NE+indexv[1]] = - hh * kp4 * oh[k]   / doh;

  /* -----   dOH / dHCO3   ----- */

  s[N2+EQOH][indexv[2]]    =  hh * km4 / doh;
  s[N2+EQOH][NE+indexv[2]] =  hh * km4 / doh;

#endif

  /* -----   dOH / dCO3   ----- */

  s[N2+EQOH][indexv[3]]    = 0.0;
  s[N2+EQOH][NE+indexv[3]] = 0.0;

  /* -----   dOH / dH   ----- */

  s[N2+EQOH][indexv[4]]    = - hh * km6 * oh[k-1] / doh;
  s[N2+EQOH][NE+indexv[4]] = - hh * km6 * oh[k]   / doh;

  /* -----   dOH / dOH   ----- */

  s[N2+EQOH][indexv[5]] =
     - hh * (
#ifndef CARTEST
     kp4 * co2[k-1]
#endif
#ifdef BORONRC4
    + kp7 * boh3[k-1]
#ifdef BORISTP
#ifdef B10B11
    + kp7bb * bboh3[k-1]
#endif
#endif
#endif
#ifdef CISTP
   + kp4cc * cco2[k-1]
#endif
   + km6 * hplus[k-1]) / doh;



  s[N2+EQOH][NE+indexv[5]] =
     - hh * (
#ifndef CARTEST
     kp4 * co2[k]
#endif
#ifdef BORONRC4
    + kp7 * boh3[k]
#ifdef BORISTP
#ifdef B10B11
    + kp7bb * bboh3[k]
#endif
#endif
#endif
#ifdef CISTP
   + kp4cc * cco2[k]
#endif
     + km6 * hplus[k]  ) / doh;



#ifdef BORONRC4

  /* -----   dOH / dB(OH)3   ----- */

  s[N2+EQOH][indexv[EQBOH3]]    =  - hh * kp7 * oh[k-1] / doh;
  s[N2+EQOH][NE+indexv[EQBOH3]] =  - hh * kp7 * oh[k]   / doh;

  /* -----   dOH / dB(OH)4   ----- */

  s[N2+EQOH][indexv[EQBOH4]]    =  hh * km7 / doh;
  s[N2+EQOH][NE+indexv[EQBOH4]] =  hh * km7 / doh;

  /* if(k == 2){
	printf("\n %e %e\n",km7,kp7*oh[k]);
	bortmp = 0;
               } */

#ifdef BORISTP
#ifdef B10B11

  /* -----   dOH / d11B(OH)3   ----- */

  s[N2+EQOH][indexv[EQBBOH3]]    =  - hh * kp7bb * oh[k-1] / doh;
  s[N2+EQOH][NE+indexv[EQBBOH3]] =  - hh * kp7bb * oh[k]   / doh;

  /* -----   dOH / d11B(OH)4   ----- */

  s[N2+EQOH][indexv[EQBBOH4]]    =  hh * km7bb / doh;
  s[N2+EQOH][NE+indexv[EQBBOH4]] =  hh * km7bb / doh;

#endif
#endif

#endif

#ifdef CISTP
  /* -----   dOH / d13CO2   ----- */

  s[N2+EQOH][indexv[EQCCO2]]    = - hh * kp4cc * oh[k-1] / doh;
  s[N2+EQOH][NE+indexv[EQCCO2]] = - hh * kp4cc * oh[k]   / doh;

  /* -----   dOH / dH13CO3   ----- */

  s[N2+EQOH][indexv[EQHCCO3]]    =  hh * km4cc / doh;
  s[N2+EQOH][NE+indexv[EQHCCO3]] =  hh * km4cc / doh;

  /* -----   dOH / dCO3   ----- */

  s[N2+EQOH][indexv[EQCCO3]]    = 0.0;
  s[N2+EQOH][NE+indexv[EQCCO3]] = 0.0;
#endif

#endif /* OH */



#ifdef C13ISTP

  /* ===== 13CO2 ================= */

  /* -----   d13CO2 / d13CO2   ----- */

  s[N2+EQCCO2][indexv[EQCCO2]]    =  -hh*(kp1scc + kp4cc * oh[k-1])/dcco2;
  s[N2+EQCCO2][NE+indexv[EQCCO2]] =  -hh*(kp1scc + kp4cc * oh[k]  )/dcco2;

  /* -----   d13CO2 / d13HCO3   ----- */

  /* kh2co3/Kh2co3 = km1s */

  s[N2+EQCCO2][indexv[EQHCCO3]] =
     hh * (km1scc * hplus[k-1] + km4cc) / dcco2;
  s[N2+EQCCO2][NE+indexv[EQHCCO3]] =
     hh * (km1scc * hplus[k]   + km4cc) / dcco2;

  /* -----   d13CO2 / d13CO3   ----- */

  s[N2+EQCCO2][indexv[EQCCO3]]    = 0.0;
  s[N2+EQCCO2][NE+indexv[EQCCO3]] = 0.0;

  /* -----   d13CO2 / dH   ----- */

  s[N2+EQCCO2][indexv[4]]    = hh * km1scc * hcco3[k-1] / dcco2;
  s[N2+EQCCO2][NE+indexv[4]] = hh * km1scc * hcco3[k]   / dcco2;

  /* -----   d13CO2 / dOH   ----- */

  s[N2+EQCCO2][indexv[5]]    = -hh * kp4cc * cco2[k-1] / dcco2;
  s[N2+EQCCO2][NE+indexv[5]] = -hh * kp4cc * cco2[k]   / dcco2;





  /* ===== H13CO3 ================= */

  /* -----   dH13CO3 / d13CO2   ----- */

  s[N2+EQHCCO3][indexv[EQCCO2]]    =
     hh * (kp1scc + kp4cc * oh[k-1]) / dhcco3;
  s[N2+EQHCCO3][NE+indexv[EQCCO2]] =
     hh * (kp1scc + kp4cc * oh[k])   / dhcco3;

  /* -----   dH13CO3 / dH13CO3   ----- */

  s[N2+EQHCCO3][indexv[EQHCCO3]] =
     - hh * (km1scc*hplus[k-1] + km4cc + km5cc) / dhcco3;
  s[N2+EQHCCO3][NE+indexv[EQHCCO3]] =
     - hh * (km1scc*hplus[k]   + km4cc + km5cc) / dhcco3;

  /* -----   dH13CO3 / d13CO3   ----- */

  s[N2+EQHCCO3][indexv[EQCCO3]]    = hh * kp5cc * hplus[k-1] / dhcco3;
  s[N2+EQHCCO3][NE+indexv[EQCCO3]] = hh * kp5cc * hplus[k]   / dhcco3;

  /* -----   dH13CO3 / dH     ----- */

  s[N2+EQHCCO3][indexv[4]]    =
   hh * (kp5cc * cco3[k-1] - km1scc*hcco3[k-1]) / dhcco3;
  s[N2+EQHCCO3][NE+indexv[4]] =
   hh * (kp5cc * cco3[k]   - km1scc*hcco3[k])   / dhcco3;

  /* -----   dH13CO3 / dOH    ----- */

  s[N2+EQHCCO3][indexv[5]]    =
      hh * kp4cc * cco2[k-1] / dhcco3;
  s[N2+EQHCCO3][NE+indexv[5]] =
      hh * kp4cc * cco2[k]   / dhcco3;





  /* ===== 13CO3 ================= */


  /* -----   d13CO3 / d13CO2   ----- */

  s[N2+EQCCO3][indexv[EQCCO2]]    = 0.0;
  s[N2+EQCCO3][NE+indexv[EQCCO2]] = 0.0;

  /* -----   d13CO3 / dH13CO3   ----- */

  s[N2+EQCCO3][indexv[EQHCCO3]]    = hh * km5cc / dcco3;
  s[N2+EQCCO3][NE+indexv[EQHCCO3]] = hh * km5cc / dcco3;

  /* -----   d13CO3 / d13CO3   ----- */

  s[N2+EQCCO3][indexv[EQCCO3]]    = - hh * kp5cc * hplus[k-1] / dcco3;
  s[N2+EQCCO3][NE+indexv[EQCCO3]] = - hh * kp5cc * hplus[k]   / dcco3;

  /* -----   d13CO3 / dH   ----- */

  s[N2+EQCCO3][indexv[4]]    = - hh * kp5cc * cco3[k-1] / dcco3;
  s[N2+EQCCO3][NE+indexv[4]] = - hh * kp5cc * cco3[k]   / dcco3;

  /* -----   d13CO3 / dOH   ----- */

  s[N2+EQCCO3][indexv[5]]    = 0.0;
  s[N2+EQCCO3][NE+indexv[5]] = 0.0;

#endif


#ifdef BORONRC4

  /* ===== B(OH)3 ================= */


  /* -----   dB(OH)3 / dOH   ----- */

  s[N2+EQBOH3][indexv[EQOH]]    =  - hh * kp7 * boh3[k-1] / dboh3;
  s[N2+EQBOH3][NE+indexv[EQOH]] =  - hh * kp7 * boh3[k]   / dboh3;

  /* -----   dB(OH)3 / dB(OH)3   ----- */

  s[N2+EQBOH3][indexv[EQBOH3]]    =  - hh * kp7 * oh[k-1] / dboh3;
  s[N2+EQBOH3][NE+indexv[EQBOH3]] =  - hh * kp7 * oh[k]   / dboh3;

  /* -----   dB(OH)3 / dB(OH)4   ----- */

  s[N2+EQBOH3][indexv[EQBOH4]]    =  hh * km7 / dboh3;
  s[N2+EQBOH3][NE+indexv[EQBOH4]] =  hh * km7 / dboh3;



  /* ===== B(OH)4 ================= */

  /* -----   dB(OH)4 / dOH   ----- */

  s[N2+EQBOH4][indexv[EQOH]]    =   hh * kp7 * boh3[k-1] / dboh4;
  s[N2+EQBOH4][NE+indexv[EQOH]] =   hh * kp7 * boh3[k]   / dboh4;

  /* -----   dB(OH)4 / dB(OH)3   ----- */

  s[N2+EQBOH4][indexv[EQBOH3]]    =  hh * kp7 * oh[k-1] / dboh4;
  s[N2+EQBOH4][NE+indexv[EQBOH3]] =  hh * kp7 * oh[k]   / dboh4;

  /* -----   dB(OH)4 / dB(OH)4   ----- */

  s[N2+EQBOH4][indexv[EQBOH4]]    = - hh * km7 / dboh4;
  s[N2+EQBOH4][NE+indexv[EQBOH4]] = - hh * km7 / dboh4;


#ifdef BORISTP

  /* ===== 11B(OH)3 ================= */

  /* -----   d11B(OH)3 / dOH   ----- */

  s[N2+EQBBOH3][indexv[EQOH]]    =  - hh * kp7bb * bboh3[k-1] / dbboh3;
  s[N2+EQBBOH3][NE+indexv[EQOH]] =  - hh * kp7bb * bboh3[k]   / dbboh3;


  /* -----   d11B(OH)3 / d11B(OH)3   ----- */

  s[N2+EQBBOH3][indexv[EQBBOH3]]    =  - hh * kp7bb * oh[k-1] / dbboh3;
  s[N2+EQBBOH3][NE+indexv[EQBBOH3]] =  - hh * kp7bb * oh[k]   / dbboh3;


  /* -----   d11B(OH)3 / d11B(OH)4   ----- */

  s[N2+EQBBOH3][indexv[EQBBOH4]]    =  hh * km7bb / dbboh3;
  s[N2+EQBBOH3][NE+indexv[EQBBOH4]] =  hh * km7bb / dbboh3;



  /* ===== 11B(OH)4 ================= */


  /* -----   d11B(OH)4 / dOH   ----- */

  s[N2+EQBBOH4][indexv[EQOH]]    =   hh * kp7bb * bboh3[k-1] / dbboh4;
  s[N2+EQBBOH4][NE+indexv[EQOH]] =   hh * kp7bb * bboh3[k]   / dbboh4;

  /* -----   d11B(OH)4 / d11B(OH)3   ----- */

  s[N2+EQBBOH4][indexv[EQBBOH3]]    =  hh * kp7bb * oh[k-1] / dbboh4;
  s[N2+EQBBOH4][NE+indexv[EQBBOH3]] =  hh * kp7bb * oh[k]   / dbboh4;

  /* -----   d11B(OH)4 / d11B(OH)4   ----- */

  s[N2+EQBBOH4][indexv[EQBBOH4]]    = - hh * km7bb / dbboh4;
  s[N2+EQBBOH4][NE+indexv[EQBBOH4]] = - hh * km7bb / dbboh4;


#endif
#endif


#ifdef BORONRC3

  /* ===== B(OH)3 ================= */


  /* -----   dB(OH)3 / dH   ----- */

  s[N2+EQBOH3][indexv[EQHP]]    =  hh * km7 * boh4[k-1] / dboh3;
  s[N2+EQBOH3][NE+indexv[EQHP]] =  hh * km7 * boh4[k]   / dboh3;

  /* -----   dB(OH)3 / dB(OH)3   ----- */

  s[N2+EQBOH3][indexv[EQBOH3]]    =  - hh * kp7 / dboh3;
  s[N2+EQBOH3][NE+indexv[EQBOH3]] =  - hh * kp7 / dboh3;

  /* -----   dB(OH)3 / dB(OH)4   ----- */

  s[N2+EQBOH3][indexv[EQBOH4]]    =  hh * km7 * hplus[k-1] / dboh3;
  s[N2+EQBOH3][NE+indexv[EQBOH4]] =  hh * km7 * hplus[k]   / dboh3;

  /* ===== B(OH)4 ================= */


  /* -----   dB(OH)4 / dH   ----- */

  s[N2+EQBOH4][indexv[EQHP]]    =  - hh * km7 * boh4[k-1] / dboh4;
  s[N2+EQBOH4][NE+indexv[EQHP]] =  - hh * km7 * boh4[k]   / dboh4;

  /* -----   dB(OH)4 / dB(OH)3   ----- */

  s[N2+EQBOH4][indexv[EQBOH3]]    =  hh * kp7 / dboh4;
  s[N2+EQBOH4][NE+indexv[EQBOH3]] =  hh * kp7 / dboh4;

  /* -----   dB(OH)4 / dB(OH)4   ----- */

  s[N2+EQBOH4][indexv[EQBOH4]]    = - hh * km7 * hplus[k-1] / dboh4;
  s[N2+EQBOH4][NE+indexv[EQBOH4]] = - hh * km7 * hplus[k]   / dboh4;


#ifdef BORISTP

  /* ===== 11B(OH)3 ================= */


  /* -----   d11B(OH)3 / dH   ----- */

  s[N2+EQBBOH3][indexv[EQHP]]    =  hh * km7bb * bboh4[k-1] / dbboh3;
  s[N2+EQBBOH3][NE+indexv[EQHP]] =  hh * km7bb * bboh4[k]   / dbboh3;

  /* -----   d11B(OH)3 / d11B(OH)3   ----- */

  s[N2+EQBBOH3][indexv[EQBBOH3]]    =  - hh * kp7bb / dbboh3;
  s[N2+EQBBOH3][NE+indexv[EQBBOH3]] =  - hh * kp7bb / dbboh3;

  /* -----   d11B(OH)3 / d11B(OH)4   ----- */

  s[N2+EQBBOH3][indexv[EQBBOH4]]    =  hh * km7bb * hplus[k-1] / dbboh3;
  s[N2+EQBBOH3][NE+indexv[EQBBOH4]] =  hh * km7bb * hplus[k]   / dbboh3;



  /* ===== 11B(OH)4 ================= */


  /* -----   d11B(OH)4 / dH   ----- */

  s[N2+EQBBOH4][indexv[EQHP]]    =  - hh * km7bb * bboh4[k-1] / dbboh4;
  s[N2+EQBBOH4][NE+indexv[EQHP]] =  - hh * km7bb * bboh4[k]   / dbboh4;

  /* -----   d11B(OH)4 / d11B(OH)3   ----- */

  s[N2+EQBBOH4][indexv[EQBBOH3]]    =  hh * kp7bb / dbboh4;
  s[N2+EQBBOH4][NE+indexv[EQBBOH3]] =  hh * kp7bb / dbboh4;

  /* -----   d11B(OH)4 / d11B(OH)4   ----- */

  s[N2+EQBBOH4][indexv[EQBBOH4]]    = - hh * km7bb * hplus[k-1] / dbboh4;
  s[N2+EQBBOH4][NE+indexv[EQBBOH4]] = - hh * km7bb * hplus[k]   / dbboh4;


#endif
#endif

#ifdef CLPL
     }
#endif

/* ----- end of REACTION ----- */

#endif

/* --- these are the discretizations at interior points ---

       j = 1,...,N2
       E_j,k    = y_j,k    - y_j,k-1    - h/2 [y_j+N2,k  + y_j+N2,k-1]  = 0
       E_j+N2,k = y_j+N2,k - y_j+N2,k-1 - h/2 [g_j (y_k) + g_j (y_k-1)] = 0

*/

        /* -----       this is diffusion        ----- */

#ifdef CLPL
      for(a=1; a <= N2; a++) {
        s[a][jsf] = y[a][k] - y[a][k-1]
            - hh * (y[N2+a][k] + y[N2+a][k-1]);
        s[N2+a][jsf] = y[N2+a][k]        - y[N2+a][k-1]
            + 0.5*(y[N2+a][k]/r[k]/r[k]/difcofm[a][k] +
                   y[N2+a][k-1]/r[k-1]/r[k-1]/difcofm[a][k-1])
                 *(r[k]*r[k]*difcofm[a][k] - r[k-1]*r[k-1]*difcofm[a][k-1]);
      }
#else
      for(a=1; a <= N2; a++) {
        s[a][jsf] = y[a][k] - y[a][k-1]
            - hh * (y[N2+a][k] + y[N2+a][k-1]);
        s[N2+a][jsf] = y[N2+a][k]        - y[N2+a][k-1]
                + h * (y[N2+a][k] / r[k] + y[N2+a][k-1] / r[k-1]);
      }
#endif


#ifdef SYMBIONTS

 /* -----   CO2 uptake by symbionts   ----- */

 /* units: SYMCO2UPT[mol/s]
         ------------------------- -> 1.e15 (1/mu**3 -> 1/l) 1.e6 (mol-> mumol)
         dco2[mu**2/s] r**2 [mu**2]

             ->   mumol
                  -----
                  l mu			   */

#if defined (FORAMSYM2) || defined (CLPL)
  if(k >= nsymradmin && k <= nsymrad){
#else
  if(k <= nsymrad){
#endif

#ifdef MIMECO2SYM
 	  a = 1;	/* CO2 */

	  /* right hand side		*/

	  tmp1  = vmaxit*co2[k]*1.e21;
	  tmp1 /= (dco2*4.*PI*r[k]*r[k]*(double)nsymrad*(KS+co2[k]));
	  s[N2+a][jsf] -= tmp1;

	  /*  derivatives dCO2/dCO2 	*/

	  dummyd  = vmaxit*KS*1.e21;
	  dummyd /= (dco2*2.*4.*PI*r[k]*r[k]*(double)nsymrad*SQ(KS+co2[k-1]));
	  s[N2+a][indexv[1]]    -= dummyd;
	  dummyd  = vmaxit*KS*1.e21;
	  dummyd /= (dco2*2.*4.*PI*r[k]*r[k]*(double)nsymrad*SQ(KS+co2[k]));
	  s[N2+a][NE+indexv[1]] -= dummyd;

	  a = 2;	/* HCO3- */

	  /* right hand side		*/

	  tmp2  = vmaxit*1.e21;
	  tmp2 /= (dhco3*4.*PI*r[k]*r[k]*(double)nsymrad);
	  tmp2 *= (1. - co2[k]/(KS+co2[k]));
	  s[N2+a][jsf] -= tmp2;

	  /*  derivatives dHCO3/dCO2	*/

	  dummyd  = (-1.)*vmaxit*KS*1.e21;
	  dummyd /= (dhco3*2.*4.*PI*r[k]*r[k]*(double)nsymrad*SQ(KS+co2[k-1]));
	  s[N2+a][indexv[1]]    -= dummyd;
	  dummyd  = (-1.)*vmaxit*KS*1.e21;
	  dummyd /= (dhco3*2.*4.*PI*r[k]*r[k]*(double)nsymrad*SQ(KS+co2[k]));
	  s[N2+a][NE+indexv[1]] -= dummyd;

#define HSYM

#ifdef HSYM
	  a = 4;	/* H+ */

	  /* right hand side		*/

	  tmp4  = vmaxit*1.e21;
	  tmp4 /= (dh*4.*PI*r[k]*r[k]*(double)nsymrad);
	  tmp4 *= (1. - co2[k]/(KS+co2[k]));
	  s[N2+a][jsf] -= tmp4;

	  /*  derivatives dH/dCO2 	*/

	  dummyd  = (-1.)*vmaxit*KS*1.e21;
	  dummyd /= (dh*2.*4.*PI*r[k]*r[k]*(double)nsymrad*SQ(KS+co2[k-1]));
	  s[N2+a][indexv[1]]    -= dummyd;
	  dummyd  = (-1.)*vmaxit*KS*1.e21;
	  dummyd /= (dh*2.*4.*PI*r[k]*r[k]*(double)nsymrad*SQ(KS+co2[k]));
	  s[N2+a][NE+indexv[1]] -= dummyd;

#endif

#ifdef OHSYM
	  a = 5;	/* OH- */

	  /* right hand side		*/

	  tmp4  = vmaxit*1.e21;
	  tmp4 /= (doh*4.*PI*r[k]*r[k]*(double)nsymrad);
	  tmp4 *= (1. - co2[k]/(KS+co2[k]));
	  s[N2+a][jsf] -= (-1.)*tmp4;

	  /*  derivatives dOH/dCO2 	*/

	  dummyd  = (-1.)*vmaxit*KS*1.e21;
	  dummyd /= (doh*2.*4.*PI*r[k]*r[k]*(double)nsymrad*SQ(KS+co2[k-1]));
	  s[N2+a][indexv[1]]    -= (-1.)*dummyd;
	  dummyd  = (-1.)*vmaxit*KS*1.e21;
	  dummyd /= (doh*2.*4.*PI*r[k]*r[k]*(double)nsymrad*SQ(KS+co2[k]));
	  s[N2+a][NE+indexv[1]] -= (-1.)*dummyd;

#endif


#ifdef C13ISTP

    a = EQCCO2;

    /*--------------------------------------------------*/
    /*							*/
    /* 		 	13CO2 			  	*/
    /* 						  	*/
    /*--------------------------------------------------*/

    	/* Total CO2 uptake at k for MiMe is calculated by
    	   tmp1. Calculate 13CO2 uptake according to:

    	   13F        =   F_ges *  R * z_i / (1 + R * z_i)

    	   R = RSTAND

    	   with  z_i for CO2:

    1. JASPER

    	   z_1  =   1/R cco2/(co2-cco2)
    	   	  - ( AEPSP - BEPS/co2) / 1000.  (1)
    	   	  				for co2 > CO2EPSP


	   z_2  =   1/R cco2/(co2-cco2)
    	   	  - ( AEPSP - BEPS/CO2EPS) * co2 / CO2EPSP / 1000.
    	   					 (2)
    	   					 for co2 < CO2EPSP

    	(formula (1) from Jasper et al. 94)
    	((2) linear for co2 goes to zero)

    2. CONSTANT EPSP

    	   z_1  =   1/R cco2/(co2-cco2)
    	   	  - EPSPCO2 / 1000.
    	   	  				for all co2    */

#ifdef JASPER
    if(co2[k] >= CO2EPSP){
      zkm1 =   cco2[k-1]/(co2[k-1]-cco2[k-1])/RSTAND
             - (AEPSP - BEPSP / co2[k-1]) / 1000.;
      z    =   cco2[k]/(co2[k]-cco2[k])/RSTAND
             - (AEPSP - BEPSP / co2[k]) / 1000.;
    }

    if(co2[k] <  CO2EPSP){
      zkm1 =   cco2[k-1]/(co2[k-1]-cco2[k-1])/RSTAND
             - (AEPSP - BEPSP / CO2EPSP)*co2[k-1]/CO2EPSP / 1000.;
      z    =   cco2[k]/(co2[k]-cco2[k])/RSTAND
             - (AEPSP - BEPSP / CO2EPSP)*co2[k]/CO2EPSP / 1000.;
    }
#endif
#ifdef CEPSP
      zkm1 =   cco2[k-1]/(co2[k-1]-cco2[k-1])/RSTAND
             - EPSPCO2 / 1000.;
      z    =   cco2[k]/(co2[k]-cco2[k])/RSTAND
             - EPSPCO2 / 1000.;
#endif


    /*--------------------------------------------------*/
    /*                                                  */
    /* This is the 13CO2 - uptake                       */
    /*                                                  */
    /*------------------------------------------------- */

    tmp1cc	  = tmp1*dco2*RSTAND*z/(1.+RSTAND*z)/dcco2;
    s[N2+a][jsf] -= tmp1cc;

#ifdef CISTP	/* correct 12C UPT
    s[N2+EQCO2][jsf] += tmp1cc;	 missing   !!
    --> skip, results are identical        !! */
#endif

    /*--------------------------------------------------*/
    /*                                                  */
    /*  DERIVATIVES                                     */
    /*                                                  */
    /*------------------------------------------------- */
    /*
    	d13F / dx_j =   F_ges * R *  (1 + R * z_i)^2
    		      * dz_i / dx_j
    		      					*/

    /* derivatives d13CO2/dCO2 */

#ifdef JASPER
    if(co2[k] >= CO2EPSP){
      dz_dx   = - cco2[k-1]/SQ(co2[k-1]-cco2[k-1])/RSTAND
    		- BEPSP / SQ(co2[k-1]) / 1000.;

      s[N2+a][indexv[1]]    -=    tmp1*RSTAND/SQ(1.+RSTAND*zkm1)*dco2
      				* dz_dx / dcco2;

      dz_dx   = - cco2[k]/SQ(co2[k]-cco2[k])/RSTAND
    		- BEPSP / SQ(co2[k]) / 1000.;

      s[N2+a][NE+indexv[1]] -=    tmp1*RSTAND/SQ(1.+RSTAND*z)*dco2
      				* dz_dx / dcco2;
     }

    if(co2[k] <  CO2EPSP){
      dz_dx   = - cco2[k-1]/SQ(co2[k-1]-cco2[k-1])/RSTAND
    		- (AEPSP - BEPSP / CO2EPSP) / CO2EPSP / 1000.;

      s[N2+a][indexv[1]]    -=    tmp1*RSTAND/SQ(1.+RSTAND*zkm1)*dco2
      				* dz_dx / dcco2;

      dz_dx   = - cco2[k]/SQ(co2[k]-cco2[k])/RSTAND
    		- (AEPSP - BEPSP / CO2EPSP) / CO2EPSP / 1000.;

      s[N2+a][NE+indexv[1]] -=    tmp1*RSTAND/SQ(1.+RSTAND*z)*dco2
      				* dz_dx / dcco2;
     }
#endif
#ifdef CEPSP
      dz_dx   = - cco2[k-1]/SQ(co2[k-1]-cco2[k-1])/RSTAND;

      s[N2+a][indexv[1]]    -=    tmp1*RSTAND/SQ(1.+RSTAND*zkm1)*dco2
      				* dz_dx / dcco2;

      dz_dx   = - cco2[k]/SQ(co2[k]-cco2[k])/RSTAND;

      s[N2+a][NE+indexv[1]] -=    tmp1*RSTAND/SQ(1.+RSTAND*z)*dco2
      				* dz_dx / dcco2;
#endif


    /* derivatives d13CO2 / d13CO2 */

      dz_dx = co2[k-1]/SQ(co2[k-1]-cco2[k-1])/RSTAND;

      s[N2+a][indexv[EQCCO2]]    -=   tmp1*RSTAND/SQ(1.+RSTAND*zkm1)*dco2
      				* dz_dx / dcco2;

      dz_dx = co2[k]/SQ(co2[k]-cco2[k])/RSTAND;

      s[N2+a][NE+indexv[EQCCO2]] -= tmp1*RSTAND/SQ(1.+RSTAND*z)*dco2
      				* dz_dx / dcco2;


  a = EQHCCO3;

    /*--------------------------------------------------*/
    /*                                                  */
    /* 		 	H13CO3                                  */
    /*                                                  */
    /*--------------------------------------------------*/



#ifdef JASPER
      zkm1 =   hcco3[k-1]/(hco3[k-1]-hcco3[k-1])/RSTAND
             - (AEPSP + eps3) / 1000.;
      z    =   hcco3[k]/(hco3[k]-hcco3[k])/RSTAND
             - (AEPSP + eps3) / 1000.;
#endif

#ifdef CEPSP
      zkm1 =   hcco3[k-1]/(hco3[k-1]-hcco3[k-1])/RSTAND
             - EPSPHCO3 / 1000.;
      z    =   hcco3[k]/(hco3[k]-hcco3[k])/RSTAND
             - EPSPHCO3 / 1000.;
#endif

    /*--------------------------------------------------*/
    /*                                                  */
    /*      This is the H13CO3 - uptake                 */
    /*                                                  */
    /*------------------------------------------------- */


    tmp2cc	  = tmp2*dhco3*RSTAND*z/(1.+RSTAND*z)/dhcco3;
    s[N2+a][jsf] -= tmp2cc;

    /*--------------------------------------------------*/
    /*                                                  */
    /*      DERIVATIVES                                 */
    /*                                                  */
    /*------------------------------------------------- */
    /*
    	d13F / dx_j =   F_ges * R *  (1 + R * z_i)^2
    		      * dz_i / dx_j
    		      					*/

    /* derivatives dH13CO3 / dHCO3 */



    dz_dx = - hcco3[k-1]/SQ(hco3[k-1]-hcco3[k-1])/RSTAND;

    s[N2+a][indexv[EQHCO3]]    -=   tmp2*RSTAND/SQ(1.+RSTAND*zkm1)*dhco3
    				  * dz_dx / dhcco3;

    dz_dx = - hcco3[k]/SQ(hco3[k]-hcco3[k])/RSTAND;

    s[N2+a][NE+indexv[EQHCO3]] -=   tmp2*RSTAND/SQ(1.+RSTAND*z   )*dhco3
    				  * dz_dx / dhcco3;


    /* derivatives dH13CO3 / dH13CO3 */



    dz_dx =  hco3[k-1]/SQ(hco3[k-1]-hcco3[k-1])/RSTAND;

    s[N2+a][indexv[EQHCCO3]]    -= tmp2*RSTAND/SQ(1.+RSTAND*zkm1)*dhco3
    				  * dz_dx / dhcco3;


    dz_dx =  hco3[k]/SQ(hco3[k]-hcco3[k])/RSTAND;

    s[N2+a][NE+indexv[EQHCCO3]] -= tmp2*RSTAND/SQ(1.+RSTAND*z   )*dhco3
    				  * dz_dx / dhcco3;

       if(k == 113){
    	printf("\n----- MIME 13C -----\n\n",r[k]);
    	printf("%e r\n",r[k]);
    	printf("%e co2\n",co2[k]);
    	printf("%e d13co2phyt \n",(tmp1cc/(tmp1-tmp1cc)/RSTAND -1.)*1000.);
    	printf("%e d13hco3phyt\n",(tmp2cc/(tmp2-tmp2cc)/RSTAND -1.)*1000.);
    	/*printf("%e d13co2\n",d13co2);
    	printf("%e cco2\n",cco2[k]);
        printf("%e hco3\n",hco3[k]);
    	printf("%e d13hco3\n",d13hco3);
    	printf("%e hcco3\n",hcco3[k]);*/
     }

#endif /* C13ISTP */


#else  /* NO MIMECO2SYM */

#ifdef CLPL
 	a = 1;	/* CO2 */
	/*
	tmp1 	      = 1.e21*SYMCO2UPT/dco2/4./PI/r[k]/r[k]
			/(double)(nsymrad-nsymradmin);
	tmp1 	      = 1.e21*SYMCO2UPT*3./dco2/4./PI/h/h
			/( KU((double)(nsymrad))-KU((double)(nsymradmin)) );
	tmp1 	      = 1.e21*SYMCO2UPT*4.*r[k]/dco2/4./PI/h/h/h
			/( (double)(nsymrad)*KU((double)(nsymrad))
			  -(double)(nsymradmin)*KU((double)(nsymradmin)) );
	*/
	tmp1 	      = 1.e21*SYMCO2UPT*5.*r[k]*r[k]/dco2/4./PI/h/h/h/h
			/( SQ((double)(nsymrad))*KU((double)(nsymrad))
			  -SQ((double)(nsymradmin))*KU((double)(nsymradmin)) );
	s[N2+a][jsf] -= tmp1;
#endif /* CLPL  */

#ifndef CLPL
#ifdef FORAMSYM2
        a = 1;	/* CO2 */

	tmp1 	      = 1.e21*SYMCO2UPT/dco2/4./PI/r[k]/r[k]/(double)(nsymrad-nsymradmin);
	s[N2+a][jsf] -= tmp1;
#else
 #ifdef AGG
        a = 1;	/* CO2 */
	tmp1 	      = 1.e21*symco2upt*h/dco2/agguptvol;
	s[N2+a][jsf] -= tmp1;
 #else
        a = 1;	/* CO2 */
	tmp1 	      = 1.e21*SYMCO2UPT/dco2/4./PI/r[k]/r[k]/(double)nsymrad;
	s[N2+a][jsf] -= tmp1;
 #endif
#endif
#endif
	a = 2; 	/* HCO3- */
	tmp2	      = 1.e21*SYMHCO3UPT/dhco3/4./PI/r[k]/r[k]/(double)nsymrad;
	s[N2+a][jsf] -= tmp2;

	a = 4;  /* H+    */
	tmp4	      = 1.e21*SYMHUPT/dh/4./PI/r[k]/r[k]/(double)nsymrad;
	s[N2+a][jsf] -= tmp4;

#ifdef C13ISTP
    a = EQCCO2;		/* 13CO2 */
    tmp1cc  	  = 1.e21*symcco2upt/dcco2/4./PI/r[k]/r[k]/(double)nsymrad;
    s[N2+a][jsf] -= tmp1cc;

    a = EQHCCO3;	/* 13HCO3- */
    tmp2cc  	  = 1.e21*symhcco3upt/dhcco3/4./PI/r[k]/r[k]/(double)nsymrad;
    s[N2+a][jsf] -= tmp2cc;


    /* H+ uptake for H13CO3 is included in total H+ uptake */
    /* a = 4;  		 H+
    tmp4cc	  = 1.e21*symhcco3upt/dh/4./PI/r[k]/r[k]/(double)nsymrad;
    s[N2+a][jsf] -= tmp4cc; */

#endif

#endif  /* MIMECO2SYM */


#ifdef FLSYMUPT
	fpsyup    = fopen("syup.sv4","a");
	fprintf(fpsyup,"%d %e %e %e %e\n",k,tmp1,tmp2,tmp4,co2[k]);
	fclose(fpsyup);
#ifdef C13ISTP
	fpsyup13    = fopen("syup13.sv4","a");
	fprintf(fpsyup13,"%d %e %e %e %e\n",k,tmp1cc,tmp2cc,0.,cco2[k]);
	fclose(fpsyup13);
#endif
#endif


#ifdef OXYGEN
#ifdef AGG
        a = EQO2;
        s[N2+a][jsf] -= 1.e21*symo2upt*h/do2/agguptvol;

	/* dissolution: release of CO32- */
#ifndef SYMCO3UPT

	a = EQCO3;
	if(co3[k] < CO3CRIT){

	tmp3	      = -pow((1.-.5*(co3[k-1]+co3[k])/CO3CRIT),NU);
	tmp3	     *=  1.e21*h*KDISS*c0/dco3/agguptvol;
	if(k == 500){
	printf("%e r\n",r[k]);
	printf("%e tmp3*dco3\n",tmp3*dco3);
	}
	if(k == 2){
	printf("%e r\n",r[k]);
	printf("%e tmp3*dco3\n",tmp3*dco3);
	}
	s[N2+a][jsf] -= tmp3;


  	/* -----   dCO3 / dCO3   -----  	*/

	dummy	      =  .5*pow((1.-.5*(co3[k-1]+co3[k])/CO3CRIT),(NU-1.))*NU/CO3CRIT;
	dummy	     *=   1.e21*h*KDISS*c0/dco3/agguptvol;
	s[N2+a][indexv[EQCO3]] -= dummy;

	dummy	      =  .5*pow((1.-.5*(co3[k-1]+co3[k])/CO3CRIT),(NU-1.))*NU/CO3CRIT;
	dummy	     *=   1.e21*h*KDISS*c0/dco3/agguptvol;
	s[N2+a][NE+indexv[EQCO3]] -= dummy;

	}
#else
	a = EQCO3;
	tmp2	      = 1.e21*SYMCO3UPT/dco3/4./PI/r[k]/r[k]/(double)nsymrad;
	s[N2+a][jsf] -= tmp2;
#endif

#else
        a = EQO2;
        s[N2+a][jsf] -= 1.e21*SYMO2UPT/do2/4./PI/r[k]/r[k]/(double)nsymrad;
#endif
#endif
	} /* end if(k < nsymrad) */


#ifdef CLPL
	/* Michaelis -Menten uptake for CO2 in the interior of the
	   cell (except chloroplast). The total uptake is stored in fdrain. */

#ifdef DRAIN

    if(k < nsymradmin){
 	  a = 1;	/* CO2 */

	  /* right hand side		*/

	  tmp1 	      = 1.e21*VMAXDRAIN*co2[k]*5.*r[k]*r[k]/dco2/4./PI/h/h/h/h
			/( SQ((double)(nsymradmin))*KU((double)(nsymradmin)) )
			/(KSDRAIN+co2[k]) ;
	  s[N2+a][jsf] -= tmp1;
	  fdrain += dco2*4.*PI*r[k]*r[k]*tmp1/1.e21;

	  /*printf("calldifeq %d %e \n",calldifeq,fdrain);*/


	  /*  derivatives dCO2/dCO2 	*/


	  dumdr  = 1.e21*VMAXDRAIN*KSDRAIN*5.*r[k]*r[k]/dco2/2./4./PI/h/h/h/h
			/( SQ((double)(nsymradmin))*KU((double)(nsymradmin)) )
			/SQ(KSDRAIN+co2[k]) ;
	  s[N2+a][NE+indexv[1]] -= dumdr;
    }
   if(k > nsymrad && k < dumk1){
 	  a = 1;	/* CO2 */

	  /* right hand side		*/

	  tmp1 	= 1.e21*VMAXDRAIN*co2[k]*5.*r[k]*r[k]/dco2/4./PI/h/h/h/h
			/( SQ((double)(dumk1))*KU((double)(dumk1))
			  -SQ((double)(nsymrad))*KU((double)(nsymrad)) )
			/(KSDRAIN+co2[k]) ;
	  s[N2+a][jsf] -= tmp1;
	  fdrain += dco2*4.*PI*r[k]*r[k]*tmp1/1.e21;

	  /*  derivatives dCO2/dCO2 	*/

	  dumdr  = 1.e21*VMAXDRAIN*KSDRAIN*5.*r[k]*r[k]/dco2/2./4./PI/h/h/h/h
			/( SQ((double)(dumk1))*KU((double)(dumk1))
			  -SQ((double)(nsymrad))*KU((double)(nsymrad)) )
			/SQ(KSDRAIN+co2[k]) ;
	  s[N2+a][NE+indexv[1]] -= dumdr;
      }

   if(calldifeq == 2){
     if(k >= nsymradmin && k <= nsymrad){

 	a = 1;	/* CO2 */

	tmp1 	      = 1.e21*fdrain*5.*r[k]*r[k]/dco2/4./PI/h/h/h/h
			/( SQ((double)(nsymrad))*KU((double)(nsymrad))
			  -SQ((double)(nsymradmin))*KU((double)(nsymradmin)) );
	s[N2+a][jsf] += tmp1;
      }
    }
#endif /* DRAIN */
#endif /* CLPL  */

#endif /* SYMBIONTS */


#ifdef REACTION

        /* -----       this is reaction        ----- */
        /* --- add up everything -> consistent    --- */

#ifdef CLPL
     if(k < dumk1 || k >= dumk2){
#endif

#ifndef NOCO2

     /* -----                CO2                ----- */


 s[N2+EQCO2][jsf] += (
     - hh * (kp1s + kp4 * oh[k-1]) * co2[k-1] / dco2
     - hh * (kp1s + kp4 * oh[k]  ) * co2[k]   / dco2
     + hh * (km1s * hplus[k-1] + km4)* hco3[k-1] / dco2
     + hh * (km1s * hplus[k]   + km4)* hco3[k]   / dco2
     );

#endif

#ifndef NOHCO3


     /* -----                HCO3                ----- */


    s[N2+EQHCO3][jsf] += (
    hh / dhco3 * (
     + kp1s * (co2[k-1] + co2[k])
     - km1s * (hplus[k-1] * hco3[k-1] + hplus[k] * hco3[k])
     - km4 * (hco3[k-1] + hco3[k])
     + kp4 * (co2[k-1] * oh[k-1] + co2[k] * oh[k])
     - km5h * (hco3[k-1] + hco3[k])
     + kp5h * (hplus[k-1] * co3[k-1] + hplus[k] * co3[k])
    ) );

#endif

#ifndef NOCO3

     /* -----                CO3                ----- */

   s[N2+EQCO3][jsf] += (
    hh / dco3 * (
   + km5h * (hco3[k-1] + hco3[k])
   - kp5h * (hplus[k-1] * co3[k-1] + hplus[k] * co3[k])
    ) );


#endif

#ifndef NOHPLUS

     /* -----                H                ----- */

#ifdef CARTEST
    s[N2+EQHP][jsf] += (
    hh / dh * (
    + (km5h ) * hco3[k-1]
    + (km5h )   * hco3[k]
    - kp5h * (hplus[k-1] * co3[k-1] + hplus[k] * co3[k])
    + 2. * kp6
    - km6 * (hplus[k-1] * oh[k-1] + hplus[k] * oh[k])
    ) );
#else
    s[N2+EQHP][jsf] += (
    hh / dh * (
    + (km5h - km1s*hplus[k-1]) * hco3[k-1]
    + (km5h - km1s*hplus[k])   * hco3[k]
    + kp1s * (co2[k-1] + co2[k])
    - kp5h * (hplus[k-1] * co3[k-1] + hplus[k] * co3[k])
    + 2. * kp6
    - km6 * (hplus[k-1] * oh[k-1] + hplus[k] * oh[k])
#ifdef BORONRC3
    + kp7 * (boh3[k-1] + boh3[k])
    - km7 * (hplus[k-1] * boh4[k-1] + hplus[k] * boh4[k])
#ifdef BORISTP
#ifdef B10B11
    + kp7bb * (bboh3[k-1] + bboh3[k])
    - km7bb * (hplus[k-1] * bboh4[k-1] + hplus[k] * bboh4[k])
#endif
#endif
#endif
#ifdef CISTP
    + (km5cc - km1scc*hplus[k-1]) * hcco3[k-1]
    + (km5cc - km1scc*hplus[k])   * hcco3[k]
    + kp1scc * (cco2[k-1] + cco2[k])
    - kp5cc * (hplus[k-1] * cco3[k-1] + hplus[k] * cco3[k])
#endif
    ) );
#endif

#endif

#ifndef NOOH
     /* -----              OH                ----- */


     s[N2+EQOH][jsf] += (
       hh / doh * (
#ifndef CARTEST
     + km4 * (hco3[k-1] + hco3[k])
     - kp4 * (co2[k-1] * oh[k-1] + co2[k] * oh[k])  /* 8/94 corrected */
#endif
#ifdef BORONRC4
    + km7 * (boh4[k-1] + boh4[k])
    - kp7 * (oh[k-1] * boh3[k-1] + oh[k] * boh3[k])
#ifdef BORISTP
#ifdef B10B11
    + km7bb * (bboh4[k-1] + bboh4[k])
    - kp7bb * (oh[k-1] * bboh3[k-1] + oh[k] * bboh3[k])
#endif
#endif
#endif
#ifdef CISTP
     + km4cc * (hcco3[k-1] + hcco3[k])
     - kp4cc * (cco2[k-1] * oh[k-1] + cco2[k] * oh[k])
#endif
     + 2.* kp6
     - km6 * (hplus[k-1] * oh[k-1] + hplus[k] * oh[k])
     ) );


#endif

#ifdef C13ISTP

     /* -----              13CO2                ----- */


 s[N2+EQCCO2][jsf] += (
     - hh * (kp1scc + kp4cc * oh[k-1]) * cco2[k-1] / dcco2
     - hh * (kp1scc + kp4cc * oh[k]  ) * cco2[k]   / dcco2
     + hh * (km1scc * hplus[k-1] + km4cc)* hcco3[k-1] / dcco2
     + hh * (km1scc * hplus[k]   + km4cc)* hcco3[k]   / dcco2
     );




     /* -----              13HCO3                ----- */


    s[N2+EQHCCO3][jsf] += (
    hh / dhcco3 * (
     + kp1scc * (cco2[k-1] + cco2[k])
     - km1scc * (hplus[k-1] * hcco3[k-1] + hplus[k] * hcco3[k])
     - km4cc * (hcco3[k-1] + hcco3[k])
     + kp4cc * (cco2[k-1] * oh[k-1] + cco2[k] * oh[k])
     - km5cc * (hcco3[k-1] + hcco3[k])
     + kp5cc * (hplus[k-1] * cco3[k-1] + hplus[k] * cco3[k])
    ) );



     /* -----               13CO3                ----- */

   s[N2+EQCCO3][jsf] += (
    hh / dcco3 * (
   + km5cc * (hcco3[k-1] + hcco3[k])
   - kp5cc * (hplus[k-1] * cco3[k-1] + hplus[k] * cco3[k])
    ) );


#endif



#ifdef BORONRC4

     /* -----              B(OH)3               ----- */

     s[N2+EQBOH3][jsf] += (
       hh / dboh3 * (
     - kp7 * (oh[k-1] * boh3[k-1] + oh[k] * boh3[k])
     + km7 * (boh4[k-1] + boh4[k]) ) );

     /* -----              B(OH)4               ----- */

     s[N2+EQBOH4][jsf] += (
       hh / dboh4 * (
     - km7 * (boh4[k-1] + boh4[k])
     + kp7 * (oh[k-1] * boh3[k-1] + oh[k] * boh3[k]) ) );

#ifdef BORISTP

    /* -----              11B(OH)3               ----- */

     s[N2+EQBBOH3][jsf] += (
       hh / dbboh3 * (
     - kp7bb * (oh[k-1] * bboh3[k-1] + oh[k] * bboh3[k])
     + km7bb * (bboh4[k-1] + bboh4[k]) ) );


     /* -----              11B(OH)4               ----- */

     s[N2+EQBBOH4][jsf] += (
       hh / dbboh4 * (
     - km7bb * (bboh4[k-1] + bboh4[k])
     + kp7bb * (oh[k-1] * bboh3[k-1] + oh[k] * bboh3[k]) ) );


#endif

#endif

#ifdef BORONRC3

     /* -----              B(OH)3               ----- */

     s[N2+EQBOH3][jsf] += (
       hh / dboh3 * (
       km7 * (hplus[k-1] * boh4[k-1] + hplus[k] * boh4[k])
     - kp7 * (boh3[k-1] + boh3[k]) ) );

     /* -----              B(OH)4               ----- */

     s[N2+EQBOH4][jsf] += (
       hh / dboh4 * (
       kp7 * (boh3[k-1] + boh3[k])
     - km7 * (hplus[k-1] * boh4[k-1] + hplus[k] * boh4[k]) ) );

#ifdef BORISTP

     /* -----              11B(OH)3               ----- */

     s[N2+EQBBOH3][jsf] += (
       hh / dbboh3 * (
       km7bb * (hplus[k-1] * bboh4[k-1] + hplus[k] * bboh4[k])
     - kp7bb * (bboh3[k-1] + bboh3[k]) ) );

     /* -----              11B(OH)4               ----- */

     s[N2+EQBBOH4][jsf] += (
       hh / dbboh4 * (
       kp7bb * (bboh3[k-1] + bboh3[k])
     - km7bb * (hplus[k-1] * bboh4[k-1] + hplus[k] * bboh4[k]) ) );

#endif

#endif


#ifdef CLPL
     }
#endif

   /* --- end of REACTION ---- */
#endif

   }

}



void precision()
{
   double prearr[20],pre1;

   int prei;

   FILE *fppre;

   fppre = fopen("pre.sv4","w");

   fprintf(fppre," precision test: double \n");

   for(prei=0; prei < 20; prei++)
      prearr[prei] = 1. + pow(0.1,(double)prei);
   for(prei=0; prei < 20; prei++) {
      pre1 = (prearr[prei] - 1.) * pow(10.,(double)prei);
      fprintf(fppre,"%d   %e   %e   \n",prei,prearr[prei],pre1);
      }

   fclose(fppre);
}
/* =================    main (begin)    ======================= */

#ifdef CBNS
void main(int argc,char *argv[])
#else
int main()
#endif
{

   int i,indexv[NE+1],j,k;
   double scalv[NE+1],x, ***c, **s, **y, **dmatrix(),
          a1,b1,a2,b2,a3,b3,err[NE+1];
#ifdef C13ISTP
   double acc1,bcc1,acc2,bcc2,acc3,bcc3;
#endif
   double cr,cinfty,ca,ak;

#ifdef CBNS
#ifdef DICBULK
   if(argc < 2){
	printf("please give arguments: 1. pH    2. DIC\n");
	return;
   }
   dicbulkarg = atof(argv[2]);
#endif
#ifdef ALKBULK
   if(argc < 2){
	printf("please give arguments: 1. pH    2. ALK\n");
	return;
   }
   alkbulkarg = atof(argv[2]);
#endif
   phbulkarg  = atof(argv[1]);
#endif

   y = dmatrix(1,NE,1,M);
   s = dmatrix(1,NE,1,NSJ);
   c = (double ***)malloc((unsigned) NE*sizeof(double **))-1;

   setbuf(stdout,NULL);

   /* ---------------------- open files --------------------- */

   fppara     = fopen("par.sv4","w");

   setbuf(fppara,NULL);

   fpks       = fopen("ks.sv4","w");

   fpr        = fopen("r.sv4","w");

   fpco2   = fopen("co2.sv4","w");
   fphco3  = fopen("hco3.sv4","w");
   fpco3   = fopen("co3.sv4","w");
   fph     = fopen("h.sv4","w");
   fpoh    = fopen("oh.sv4","w");

   fpdc    = fopen("dc.sv4","w");

#ifdef C13ISTP
   fpcco2   = fopen("cco2.sv4","w");
   fphcco3  = fopen("hcco3.sv4","w");
   fpcco3   = fopen("cco3.sv4","w");
   fpdc13s  = fopen("dc13s.sv4","w");
   fpdcc    = fopen("dcc.sv4","w");

#ifdef CBNSLOOP
   fpdc13slp = fopen("dc13slp.sv4","w");
 #ifdef READ
  #ifdef DDAT
   #ifndef LIDA
   	fpread    = fopen("alkcd.dat","r");
   	if (fpread == 0){
   		printf("file alkcd.dat not available.\n");
    		exit(1);
   	}
   #endif
   #ifdef LIDA
   	fpread    = fopen("alkcl.dat","r");
   	if (fpread == 0){
   		printf("file alkcl.dat not available.\n");
    		exit(1);
   	}
   #endif
  #endif
  #ifdef LDAT
   fpread    = fopen("alkcl.dat","r");
   if (fpread == 0){
   	printf("file alkcl.dat not available.\n");
    	exit(1);
   }
  #endif
 #endif
#endif

#endif /* C13ISTP */

#ifdef OXYGEN
   fpo2    = fopen("o2.sv4","w");
#endif
#ifdef BORON
   fpboh3  = fopen("boh3.sv4","w");
   fpboh4  = fopen("boh4.sv4","w");
#ifdef BORISTP
   fpbboh3   = fopen("bboh3.sv4","w");
   fpbboh4   = fopen("bboh4.sv4","w");
#endif
#endif
#ifdef CALCIUM
   fpca    = fopen("ca.sv4","w");
#endif



   /* ------------------------------------------------------- */

   fprintf(fppara,"==================================== \n");
   fprintf(fppara,"=========== par.solvde4 ============ \n");
   fprintf(fppara,"==================================== \n");


#ifdef CBNSLOOP
 #ifndef READ
   for(icl=1;icl<=ICLMAX;icl++){
	phbulkv  = PHMIN  + (icl-1)*(PHMAX-PHMIN)/(ICLMAX-1);
 #endif
 #ifdef READ
  for(icl=1;icl<=NDAT;icl++){
    n  = fscanf(fpread,"%lf %lf %lf\n",&d,&phin[icl-1],&alkin[icl-1]);
    n += fscanf(fpread,"%lf %lf %lf %lf %lf\n",&dicin[icl-1],&d,&d,&d,&d);
    if(n!=8) printf("read error in alkc*.dat\n");
  }
  for(icl=1;icl<=NDAT;icl++){
 #endif
#endif


#ifdef C13ISTP
   initeps(); /* initialize eps-values for fractionation of 13C */
#endif

   initc();    /* initialize concentrations */

   initk();    /* initialize rate coefficients k */

#ifdef INITA
   dt = 0.005;    /* time step */
   inita();
#endif

   /* precision(); */

/* write rate coefficients k (-> Matlab) */

   fprintf(fpks,"%e \n",kp1s);
   fprintf(fpks,"%e \n",km1s);
   fprintf(fpks,"%e \n",kp4);
   fprintf(fpks,"%e \n",km4);
   fprintf(fpks,"%e \n",kp5h);
   fprintf(fpks,"%e \n",km5h);
   fprintf(fpks,"%e \n",kp6);
   fprintf(fpks,"%e \n",km6);

   fprintf(fpks,"%e \n",dco2);
   fprintf(fpks,"%e \n",dhco3);
   fprintf(fpks,"%e \n",dco3);
   fprintf(fpks,"%e \n",dh);
   fprintf(fpks,"%e \n",doh);

   fclose(fpks);

   for(i=1; i<=NE; i++) c[i] = dmatrix(1,NCJ,1,NCK);

   fprintf(fppara,"--- solvde4.c  --- \n");

#ifdef REACTION
   fprintf(fppara,"--- diffusion and reaction  --- \n");
#endif

#ifdef NOCO2
   fprintf(fppara,"--- NOCO2 defined --- \n");
#endif

#ifdef NOHCO3
   fprintf(fppara,"--- NOHCO3 defined --- \n");
#endif

#ifdef NOCO3
   fprintf(fppara,"--- NOCO3 defined --- \n");
#endif

#ifdef NOHPLUS
   fprintf(fppara,"--- NOHPLUS defined --- \n");
#endif

#ifdef NOOH
   fprintf(fppara,"--- NOOH defined --- \n");
#endif

   fprintf(fppara,"EQCO2   %d \n",EQCO2);
   fprintf(fppara,"EQHCO3  %d \n",EQHCO3);
   fprintf(fppara,"EQCO3   %d \n",EQCO3);
   fprintf(fppara,"EQHP    %d \n",EQHP);
   fprintf(fppara,"EQOH    %d \n",EQOH);
#ifdef C13ISTP
   fprintf(fppara,"EQCCO2   %d \n",EQCCO2);
   fprintf(fppara,"EQHCCO3  %d \n",EQHCCO3);
   fprintf(fppara,"EQCCO3   %d \n",EQCCO3);
#endif
#ifdef BORON
   fprintf(fppara,"EQBOH3  %d \n",EQBOH3);
   fprintf(fppara,"EQBOH4  %d \n",EQBOH4);
#endif
#ifdef OXYGEN
   fprintf(fppara,"EQO2    %d \n",EQO2);
#endif
#ifdef CALCIUM
   fprintf(fppara,"EQCA    %d \n",EQCA);
#endif

#ifdef PRINT
   printf("EQCO2    %d  \n",EQCO2);
   printf("EQHCO3   %d  \n",EQHCO3);
   printf("EQCO3    %d  \n",EQCO3);
   printf("EQHP     %d  \n",EQHP);
   printf("EQOH     %d  \n",EQOH);
#ifdef C13ISTP
   printf("EQCCO2   %d \n",EQCCO2);
   printf("EQHCCO3  %d \n",EQHCCO3);
   printf("EQCCO3   %d \n",EQCCO3);
#endif
#ifdef BORON
   printf("EQBOH3   %d  \n",EQBOH3);
   printf("EQBOH4   %d  \n",EQBOH4);
#endif
#ifdef OXYGEN
   printf("EQO2     %d  \n",EQO2);
#endif
#ifdef CALCIUM
   printf("EQCA     %d  \n",EQCA);
#endif

#endif /* PRINT */

   fprintf(fppara,"(N2) number of equations of 2. order %d \n",N2);
   fprintf(fppara,"(NE) number of equations of 1. order %d \n",NE);
   fprintf(fppara,"(NLB) number of left b.c.            %d \n",NLB);
   fprintf(fppara,"(NRB) number of right b.c.           %d \n",NRB);

   fprintf(fppara,"dimensions of y-array  %d  %d \n",NE,M);
   fprintf(fppara,"dimensions of s-array  %d  %d \n",NE,NSJ);
   fprintf(fppara,"dimensions of c-array  %d  %d  %d \n",NE,NCJ,NCK);



   indexv[0] = -1;     /* not used */

#ifdef UPTAKE
   for(i = 1; i <= N2; i++) {
      indexv[   i] = N2+i;
      indexv[N2+i] = i;
   }
#else
   for(i = 1; i <= N2; i++) {
      indexv[   i] = i;
      indexv[N2+i] = N2+i;
   }
#endif

#ifdef PRINT
      for(i=1;i <= NE; i++)
         printf("%d  %d  indexv[i]\n",indexv[i],i);
#endif

   /* --- initial guess --- */


   /* --- test: pure diffusion --- */

   /*
        c'' + 2/r c' = 0                  -> c(r) = a/r + b
        dc/dr (r1) = alpha = - a / r1**2  -> a = -alpha * r1**2
        c(r2) = beta = a/r2 + b           -> b = beta + alpha * r1**2 / r2
   */


      a1 = -co2flux * RADIUS * RADIUS;
      b1 =  co2bulk - a1 / RBULK;

      a2 = -hco3flux * RADIUS * RADIUS;
      b2 =  hco3bulk - a2 / RBULK;

      a3 = -co3flux * RADIUS * RADIUS;
      b3 =  co3bulk - a3 / RBULK;

#ifdef C13ISTP
      acc1 = -cco2flux * RADIUS * RADIUS;
      bcc1 = cco2bulk - acc1 / RBULK;

      acc2 = -hcco3flux * RADIUS * RADIUS;
      bcc2 = hcco3bulk - acc2 / RBULK;

      acc3 = -cco3flux * RADIUS * RADIUS;
      bcc3 =  cco3bulk - acc3 / RBULK;
#endif

      h = (RBULK - RADIUS)/(double)(M-1);
      hh = 0.5 * h;

      fprintf(fppara,"h [mu] grid spacing       %e \n",h);

#ifdef DEB05
   printf("--- before initial guess loop --- \n");
   printf("%e  %e  a1,b1 \n",a1,b1);
   printf("%e  %e  a2,b2 \n",a2,b2);
#endif

   r[0] = 0.0;    /* not used */
   for(k=1; k<=M; k++) {
     x = RADIUS + h*(double)(k-1);
     r[k] = x;
        y[EQCO2][k]  =  a1/x + b1;
     y[N2+EQCO2][k]  = -a1/x/x;
        y[EQHCO3][k] =  a2/x + b2;
     y[N2+EQHCO3][k] = -a2/x/x;
        y[EQCO3][k]  =  co3bulk;
     y[N2+EQCO3][k]  =  0.0;
        y[EQHP][k]   =   hbulk;
        y[EQOH][k]   =  ohbulk;
     y[N2+EQHP][k]   =  0.0;
     y[N2+EQOH][k]   =  0.0;
#ifdef C13ISTP
        y[EQCCO2][k]   =  acc1/x + bcc1;
     y[N2+EQCCO2][k]   = -acc1/x/x;
        y[EQHCCO3][k]  =  acc2/x + bcc2;
     y[N2+EQHCCO3][k]  = -acc2/x/x;
        y[EQCCO3][k]   =  cco3bulk;
     y[N2+EQCCO3][k]   =  0.0;
#endif
#ifdef BORON
        y[EQBOH3][k]   = boh3bulk;
     y[N2+EQBOH3][k]   = 0.0;
        y[EQBOH4][k]   = boh4bulk;
     y[N2+EQBOH4][k]   = 0.0;
#ifdef BORISTP
        y[EQBBOH3][k]   = bboh3bulk;
     y[N2+EQBBOH3][k]   = 0.0;
        y[EQBBOH4][k]   = bboh4bulk;
     y[N2+EQBBOH4][k]   = 0.0;
#endif
#endif
#ifdef OXYGEN
        y[EQO2][k]   = o2bulk;
     y[N2+EQO2][k]   = 0.0;
#endif
#ifdef CALCIUM
        y[EQCA][k]   = cabulk;
     y[N2+EQCA][k]   = 0.0;
#endif
   }

   fprintf(fppara,"------------------------- \n");
   fprintf(fppara,"%e   r[1]                 \n",r[1]);
   fprintf(fppara,"%e   r[M]                 \n",r[M]);
   fprintf(fppara,"------------------------- \n");

#ifdef SYMBIONTS
#if defined (FORAMSYM2) || defined (CLPL)
     for(k=1; k<=M; k++) {
		if(r[k] > SYMRADMIN){
			nsymradmin = k - 1;
			k = M + 2;
       	}
     }
     for(k=1; k<=M; k++) {
		if(r[k] > SYMRAD){
			nsymrad = k - 1;
			k = M + 2;
       	}
     }
     printf("%d  %e   nsymradmin,r[nsymradmin] \n",nsymradmin,r[nsymradmin]);
     printf("%d  %e   nsymrad,r[nsymrad] \n",nsymrad,r[nsymrad]);
     /* write rmin and rmax of symbiont halo in file */
     fprintf(fpdc,"%e\n",SYMRADMIN);
     fprintf(fpdc,"%e\n",r[nsymrad]);
#ifdef CLPL
     /* determine ncmembmin and max */
     for(k=1; k<=M; k++) {
		if(r[k] > CELLRADIUS){
			ncmembmin = k - 1;
			k = M + 2;
       	}
     }
     for(k=1; k<=M; k++) {
		if(r[k] > (CELLRADIUS+MEMBTHK)){
			ncmembmax = k - 1;
			k = M + 2;
       	}
     }
     printf("%d  %e   ncmembmin,r[ncmembmin] \n",ncmembmin,r[ncmembmin]);
     printf("%d  %e   ncmembmax,r[ncmembmax] \n",ncmembmax,r[ncmembmax]);

     initdm();	/* init matrix of diffusion coefficients */

   fpdifcofm = fopen("difcofm.sv4","w");
   for(k=1;k<=M;k++){
   	for(j=1;j<=N2;j++)
   		fprintf(fpdifcofm,"%e ",difcofm[j][k]*1.e-12);
	fprintf(fpdifcofm,"\n");
   }
   fclose(fpdifcofm);

#endif
#else
#ifdef AGG
     agguptvol = (4.*PI*(KU(AGGRADIUS-h)-KU(RADIUS))/3.); /* upt. volume (mum3)*/
     for(k=1; k<=M; k++) {
        if(r[k] > (SYMRAD-h)) {
          nsymrad = k - 1;
#else
     for(k=1; k<=M; k++) {
        if(r[k] > SYMRAD) {
          nsymrad = k - 1;
#endif
#ifdef PRINT
          printf("%d  %e   nsymrad,r[nsymrad] \n",nsymrad,r[nsymrad]);
#endif /* PRINT */
	  /* write rmin and rmax of symbiont halo in file */
	  fprintf(fpdc,"%e\n",RADIUS);
          fprintf(fpdc,"%e\n",r[nsymrad]);
          k = M + 2;
        }
      }
#endif
#endif


#ifdef DEB05
   printf("--- after initial guess loop --- \n");
#endif

   /* --- scalv: typical values --- */

#ifdef DEB05
   printf("--- before scalv loop --- \n");
#endif

   scalv[0]    = 0.0;     /*   not used   */

   for(i=1; i <= N2; i++)  {
     scalv[i]    = fabs(y[i][M]);
     scalv[N2+i] = fabs(y[i][M] / (RBULK - RADIUS) );
#ifdef PRINT
     printf("%e  %e  %d  scalv[i] scalv[N2+i] \n",
            scalv[i],scalv[N2+i],i);
#endif
   }

#ifdef DEB05
   printf("--- after scalv loop --- \n");
#endif

   /*   exit(1);  */

   printf("--- before solvde --- \n");

   /* -----   test: set some k's to zero   ----- */

#ifdef SCR
   km6 = 0.0;
   kp6 = 0.0;
   printf("--- km6 = kp6 = 0 before solvde --- \n");
#endif

      for(i=1;i <= NE; i++) {
      for(j=1;j <= NSJ;j++) {
         s[i][j] = 0.0;
      }}

   solvde(ITMAX,CONV,SLOWC,scalv,indexv,NE,NB,M,y,c,s);

   printf("--- after solvde --- \n");

   printf("\n %e co2 at the shell\n\n",y[EQCO2][1]);

#ifdef C13ISTP
#ifdef F13_CO3

      if(CO3UPT != 0.0){
     	tmp1 = (cco3upt/(CO3UPT-cco3upt)/RSTAND - 1.)*1000;
#ifdef LDAT
     	tmp1 = (cco3upt/(co3uptldat-cco3upt)/RSTAND - 1.)*1000;
#endif
      	tmp2 = (y[EQCCO3][1]/(y[EQCO3][1]-y[EQCCO3][1])/RSTAND - 1.)*1000.;
#ifdef CISTP
      	tmp2 = (y[EQCCO3][1]/(y[EQCO3][1]             )/RSTAND - 1.)*1000.;
#endif
      	printf("%e d13C of the shell\n",tmp1);
      	printf("%e d13CO3(R1)       \n",tmp2);
	fprintf(fpdc13s,"%e\n",tmp1);
	fprintf(fpdc13s,"%e\n",TK);
	fprintf(fpdc13s,"%e\n",CO3UPT);
	fprintf(fpdc13s,"%e\n",CO2UPT);
#ifdef SYMBIONTS
	fprintf(fpdc13s,"%e\n",SYMCO2UPT);
#else
	fprintf(fpdc13s,"%e\n",0.0);
#endif
	fprintf(fpdc13s,"%e\n",kp4);
      }
      if(CO3UPT == 0.0){
	tmp1 = d13hco3bulk*(1.+eps5/1000.) + eps5;
      	tmp2 = (y[EQCCO3][1]/(y[EQCO3][1]-y[EQCCO3][1])/RSTAND - 1.)*1000;
#ifdef CISTP
      	tmp2 = (y[EQCCO3][1]/(y[EQCO3][1]             )/RSTAND - 1.)*1000.;
#endif
      	printf("%e d13C of the shell\n",tmp1);
      	printf("%e d13CO3(R1)       \n",tmp2);
	fprintf(fpdc13s,"%e\n",tmp1);
	fprintf(fpdc13s,"%e\n",TK);
	fprintf(fpdc13s,"%e\n",CO3UPT);
	fprintf(fpdc13s,"%e\n",CO2UPT);
	fprintf(fpdc13s,"%e\n",kp4);
      }
#endif

#ifdef F13_HCO3

      if(HCO3UPT != 0.0){
     	tmp1 = (hcco3upt/(HCO3UPT-hcco3upt)/RSTAND - 1.)*1000;
#ifdef LDAT
     	tmp1 = (hcco3upt/(hco3uptldat-hcco3upt)/RSTAND - 1.)*1000;
#endif
      	tmp2 = (y[EQHCCO3][1]/(y[EQHCO3][1]-y[EQHCCO3][1])/RSTAND - 1.)*1000.;
#ifdef CISTP
      	tmp2 = (y[EQHCCO3][1]/(y[EQHCO3][1]             )/RSTAND - 1.)*1000.;
#endif
      	printf("%e d13C of the shell\n",tmp1);
      	printf("%e d13HCO3(R1)       \n",tmp2);
	fprintf(fpdc13s,"%e\n",tmp1);
	fprintf(fpdc13s,"%e\n",TK);
	fprintf(fpdc13s,"%e\n",CO3UPT);
	fprintf(fpdc13s,"%e\n",CO2UPT);
#ifdef SYMBIONTS
	fprintf(fpdc13s,"%e\n",SYMCO2UPT);
#else
	fprintf(fpdc13s,"%e\n",0.0);
#endif
	fprintf(fpdc13s,"%e\n",kp4);
      }
      if(HCO3UPT == 0.0){
	tmp1 = d13hco3bulk*(1.+eps5/1000.) + eps5;
      	tmp2 = (y[EQHCCO3][1]/(y[EQHCO3][1]-y[EQHCCO3][1])/RSTAND - 1.)*1000;
#ifdef CISTP
      	tmp2 = (y[EQHCCO3][1]/(y[EQHCO3][1]             )/RSTAND - 1.)*1000.;
#endif
      	printf("%e d13C of the shell\n",tmp1);
      	printf("%e d13HCO3(R1)       \n",tmp2);
	fprintf(fpdc13s,"%e\n",tmp1);
	fprintf(fpdc13s,"%e\n",TK);
	fprintf(fpdc13s,"%e\n",CO3UPT);
	fprintf(fpdc13s,"%e\n",CO2UPT);
	fprintf(fpdc13s,"%e\n",kp4);
      }
#endif


#endif


#ifdef BORISTP
	/* d11B(BOH)4 at the shell */

  tmpb = ( (y[EQBBOH4][1]/y[EQBOH4][1])/BSTAND - 1.)*1000;
  printf("\n%e d11boh4[1] at the shell\n",tmpb);
  tmpb = ( D11BSEA*BORTBULK-EPSB*(y[EQBBOH3][1]+y[EQBOH3][1]) )
  		/ ( (1.+EPSB/1000.)*(y[EQBBOH3][1]+y[EQBOH3][1]) + (y[EQBBOH4][1]+y[EQBOH4][1]) );
  printf("%e d11boh4(pH[1]) at the shell\n\n",tmpb);
#endif


#ifdef CBNSLOOP
	fprintf(fpdc13slp,"%e %e %e %e ",TEMP,SALINITY,phbulk,dicbulk);
        fprintf(fpdc13slp,"%e %e %e ",co2bulk,hco3bulk,co3bulk);
        fprintf(fpdc13slp,"%e\n",tmp1);
      }		/* end of Carbon-system-loop */
#endif
      for(j=1;j<=M;j++) {
        fprintf(fpr,"%f\n",r[j]);
        fprintf(fpco2,"%e\n", y[EQCO2][j]);
        fprintf(fphco3,"%e\n",y[EQHCO3][j]);
        fprintf(fpco3,"%e\n", y[EQCO3][j]);
        fprintf(fph,"%e\n",   y[EQHP][j]);
        fprintf(fpoh,"%e\n",  y[EQOH][j]);
#ifdef C13ISTP
  	fprintf(fpcco2,"%e\n", y[EQCCO2][j]);
        fprintf(fphcco3,"%e\n",y[EQHCCO3][j]);
        fprintf(fpcco3,"%e\n", y[EQCCO3][j]);
#endif
#ifdef OXYGEN
        fprintf(fpo2,"%e\n", y[EQO2][j]);
#endif
#ifdef BORON
        fprintf(fpboh3,"%e\n", y[EQBOH3][j]);
        fprintf(fpboh4,"%e\n", y[EQBOH4][j]);
#ifdef BORISTP
        fprintf(fpbboh3,"%e\n", y[EQBBOH3][j]);
        fprintf(fpbboh4,"%e\n", y[EQBBOH4][j]);
#endif
#endif
#ifdef CALCIUM
        fprintf(fpca,"%e\n", y[EQCA][j]);
#endif
      }


#ifdef ANASOL

/* ===================================================

   analytical solution:
   (c-cinfty)/(ca-cinfty) = a/r * exp( (a-r)/ak);
   ak = sqrt(D/k');
   k' = kp4*oh+k32;
   Qa = 4*pi*a*D*(1+a/ak)*(cinfty-ca)

   =================================================== */

   cinfty = y[EQCO2][M];
   ak = sqrt(dco2/(kp4*y[EQOH][M]+kp1s));

   fprintf(fppara,"------------------------- \n");
   fprintf(fppara,"%e  ak[mu]                \n",ak);
   fprintf(fppara,"------------------------- \n");

   printf("%e k32    [1/s]       \n",k32);
   printf("%e kp1s   [1/s]       \n",kp1s);

/* e6 mol -> mumol; e15 1/mu**3 -> 1/l */
   ca = cinfty - CO2UPT*1.e21/(4.*PI*RADIUS*dco2*(1.+RADIUS/ak));

   cr = 4.*PI*RADIUS*dco2*(1.+RADIUS/ak);
   printf("%e scr1    [mu**3/s]       \n",cr);

   cr = CO2UPT*1.e21/(4.*PI*RADIUS*dco2*(1.+RADIUS/ak));
   printf("%e scr2    [mumol/kg]  \n",cr);

   printf("%e CO2UPT  [mol/s]  \n",CO2UPT);
   printf("%e PI      []       \n",PI);
   printf("%e RADIUS  []       \n",RADIUS);
   printf("%e dco2    []       \n",dco2);

   printf("%e cinfty [mumol/kg] \n",cinfty);
   printf("%e ak     [mu]      \n",ak);
   printf("%e ca     [mumol/kg] \n",ca);

   fpanasol = fopen("anasol.sv4","w");

   for(j=1;j<=M;j++) {
     cr = cinfty + (ca - cinfty)*RADIUS/r[j]
         *exp( (RADIUS - r[j])/ak );
     fprintf(fpanasol,"%e \n",cr);
   }

   fclose(fpanasol);

#endif

#define FLUXTEST

#ifdef FLUXTEST

   fprintf(fppara,"---   dc/dr at the surface (solution)   --- \n");

   co2flux = (y[EQCO2][2]  - y[EQCO2][1])  / (r[2] - r[1]);
     hflux = (y[EQHP][2]   - y[EQHP][1] )  / (r[2] - r[1]);
  hco3flux = (y[EQHCO3][2] - y[EQHCO3][1]) / (r[2] - r[1]);
   co3flux = (y[EQCO3][2]  - y[EQCO3][1])  / (r[2] - r[1]);
    ohflux = (y[EQOH][2]   - y[EQOH][1])   / (r[2] - r[1]);

#ifdef CARTEST
    printf("%e\n",y[EQCO3][2]);
    printf("%e\n",y[EQCO3][1]);
    printf("%e\n",r[2]);
    printf("%e\n",r[1]);
#endif

   fprintf(fppara,"co2flux  [mol/kg/mu]       %e \n",co2flux);
   fprintf(fppara,"hflux    [mol/kg/mu]       %e \n",hflux);
   fprintf(fppara,"hco3flux [mol/kg/mu]       %e \n",hco3flux);
   fprintf(fppara,"co3flux  [mol/kg/mu]       %e \n",co3flux);
   fprintf(fppara,"ohflux   [mol/kg/mu]       %e \n",ohflux);

   if(fabs(co2fluxs) > 0.001) {
     aux1 = 100.*(co2flux - co2fluxs)/co2fluxs;
     fprintf(fppara,"d(co2flux)             [percent]    %e \n",aux1);

     aux2 = 100.*(hflux - hfluxs)/co2fluxs;
     fprintf(fppara,"d(hflux)/d(co2flux)    [percent]    %e \n",aux2);

     aux2 = 100.*(hco3flux - hco3fluxs)/co2fluxs;
     fprintf(fppara,"d(hco3flux)/d(co2flux) [percent]    %e \n",aux2);

     aux2 = 100.*(co3flux - co3fluxs)/co2fluxs;
     fprintf(fppara,"d(co3flux)/d(co2flux)  [percent]    %e \n",aux2);

     aux2 = 100.*(ohflux - ohfluxs)/co2fluxs;
     fprintf(fppara,"d(ohflux)/d(co2flux)   [percent]    %e \n",aux2);
   }

   if(fabs(hfluxs) > 0.001) {
     aux1 = 100.*(hflux - hfluxs)/hfluxs;
     fprintf(fppara,"d(hflux) [percent]          %e \n",aux1);
   }

   if(fabs(hco3fluxs) > 0.001) {
     aux1 = 100.*(hco3flux - hco3fluxs)/hco3fluxs;
     fprintf(fppara,"d(hco3flux) [percent]       %e \n",aux1);
   }

   if(fabs(co3fluxs) > 0.001) {
     aux1 = 100.*(co3flux - co3fluxs)/co3fluxs;
     fprintf(fppara,"d(co3flux) [percent]        %e \n",aux1);
   }
   if(fabs(ohfluxs) > 0.001) {
     aux1 = 100.*(ohflux - ohfluxs)/ohfluxs;
     fprintf(fppara,"d(ohflux) [percent]         %e \n",aux1);
   }


  fprintf(fppara,"---  uptake --- \n");

   co2flux  =  co2flux * surface * dco2  / 1.e21;
     hflux  =    hflux * surface * dh    / 1.e21;
  hco3flux  = hco3flux * surface * dhco3 / 1.e21;
   co3flux  =  co3flux * surface * dco3  / 1.e21;
    ohflux  =   ohflux * surface * doh   / 1.e21;

   fprintf(fppara,"co2upt        %e \n",co2flux);
   fprintf(fppara,"hupt          %e \n",hflux);
   fprintf(fppara,"hco3upt       %e \n",hco3flux);
   fprintf(fppara,"co3upt        %e \n",co3flux);
   fprintf(fppara,"ohupt         %e \n",ohflux);

   dummy = hflux - hco3flux - 2.*co3flux - ohflux;   /* neutral? */

   fprintf(fppara,"neutral?      %e \n",dummy);


   fprintf(fppara,"---   dc/dr bulk (solution)   --- \n");

   co2flux = (y[EQCO2][M]  - y[EQCO2][M-1])  / (r[2] - r[1]);
     hflux = (y[EQHP][M]   - y[EQHP][M-1] )  / (r[2] - r[1]);
  hco3flux = (y[EQHCO3][M] - y[EQHCO3][M-1]) / (r[2] - r[1]);
   co3flux = (y[EQCO3][M]  - y[EQCO3][M-1])  / (r[2] - r[1]);
    ohflux = (y[EQOH][M]   - y[EQOH][M-1])   / (r[2] - r[1]);

   fprintf(fppara,"co2flux        %e \n",co2flux);
   fprintf(fppara,"hflux          %e \n",hflux);
   fprintf(fppara,"hco3flux       %e \n",hco3flux);
   fprintf(fppara,"co3flux        %e \n",co3flux);
   fprintf(fppara,"ohflux         %e \n",ohflux);

  fprintf(fppara,"---  uptake --- \n");

   surface = surface / RADIUS / RADIUS * RBULK * RBULK;

   co2flux  =  co2flux * surface * dco2  / 1.e21;
     hflux  =    hflux * surface * dh    / 1.e21;
  hco3flux  = hco3flux * surface * dhco3 / 1.e21;
   co3flux  =  co3flux * surface * dco3  / 1.e21;
    ohflux  =   ohflux * surface * doh   / 1.e21;

   fprintf(fppara,"co2upt        %e \n",co2flux);
   fprintf(fppara,"hupt          %e \n",hflux);
   fprintf(fppara,"hco3upt       %e \n",hco3flux);
   fprintf(fppara,"co3upt        %e \n",co3flux);
   fprintf(fppara,"ohupt         %e \n",ohflux);

   dummy = hflux - hco3flux - 2.*co3flux - ohflux;   /* neutral? */

   fprintf(fppara,"neutral?      %e \n",dummy);

#endif






      fclose(fpr);
      fclose(fpco2);
      fclose(fphco3);
      fclose(fpco3);
      fclose(fph);
      fclose(fpoh);

      fclose(fpdc);

#ifdef C13ISTP
      fclose(fpcco2);
      fclose(fphcco3);
      fclose(fpcco3);
      fclose(fpdc13s);
      fclose(fpdcc);

#ifdef CBNSLOOP
      fclose(fpdc13slp);
 #ifdef READ
      fclose(fpread);
 #endif
#endif
#endif
#ifdef OXYGEN
      fclose(fpo2);
#endif
#ifdef BORON
      fclose(fpboh3);
      fclose(fpboh4);
#ifdef BORISTP
      fclose(fpbboh3);
      fclose(fpbboh4);
      fclose(fpdbboh3);
      fclose(fpdbboh4);
#endif
#endif
#ifdef CALCIUM
      fclose(fpca);
#endif
      fclose(fppara);


/* ====================    end of main   =============================== */

}

