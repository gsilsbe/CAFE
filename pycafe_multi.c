#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <dirent.h> 

float * cafe(float adg_443, 
	     float aph_443, 
             float bbp_443,
             float bbp_s,
	     float bbwater,
    	     float chl,
	     float lat,
	     float mld,
	     float PAR,
	     int yday)
{

/*------------------------------------------------------------------------------
	Step 1: Declare all variables and dependent functions                       */

int w; /* Wavelength step */ 
int t; /* Time step       */
int z; /* Depth step      */

/* Inherent Optical Properties*/

float aphi[31];   /* Phytoplankton absorption coefficient [m-1] */
float a[31];      /* Total absorption coefficient [m-1]         */ 
float bb[31];     /* total backscattering coefficient [m-1]     */
float bbw[31];    /* Backscattering of pure water [m-1]         */

/* Wavelength [nm] */
int wv[] = { 400, 410, 420, 430, 440, 450, 460, 470, 480, 490, 500, 510, 520, 
             530, 540, 550, 560, 570, 580, 590, 600, 610, 620, 630, 640, 650, 
             660, 670, 680, 690, 700 };

/* Pure water absorption coefficient [m-1] */
float aw[] = { 0.00663, 0.00473, 0.00454, 0.00495, 0.00635, 0.00922, 0.00979, 
               0.0106, 0.0127, 0.015, 0.0204, 0.0325, 0.0409, 0.0434, 0.0474, 
               0.0565, 0.0619, 0.0695, 0.0896, 0.1351, 0.2224, 0.2644, 0.2755, 
               0.2916, 0.3108, 0.34, 0.41, 0.439, 0.465, 0.516, 0.624 };

/* Spectral shape of aphi (Bricaud et al. 1998) */
float A_Bricaud[] = { 0.0241, 0.0287, 0.0328, 0.0359, 0.0378, 0.0350, 0.0328, 
                      0.0309, 0.0281, 0.0254, 0.0210, 0.0162, 0.0126, 0.0103, 
                      0.0085, 0.0070, 0.0057, 0.0050, 0.0051, 0.0054, 0.0052, 
                      0.0055, 0.0061, 0.0066, 0.0071, 0.0078, 0.0108, 0.0174, 
                      0.0161, 0.0069, 0.0025 };

/* Spectral shape of aphi (Bricaud et al. 1998) */
float E_Bricaud[] = { 0.6877, 0.6834, 0.6664, 0.6478, 0.6266, 0.5993, 0.5961, 
                      0.5970, 0.5890, 0.6074, 0.6529, 0.7212, 0.7939, 0.8500, 
                      0.9036, 0.9312, 0.9345, 0.9298, 0.8933, 0.8589, 0.8410, 
                      0.8548, 0.8704, 0.8638, 0.8524, 0.8155, 0.8233, 0.8138, 
                      0.8284, 0.9255, 1.0286 };

/* Assumed solar PAR spectrum */
float PAR_spectrum[] = { 0.00227, 0.00218, 0.00239, 0.00189, 0.00297, 0.00348, 
                         0.00345, 0.00344, 0.00373, 0.00377, 0.00362, 0.00364, 
                         0.00360, 0.00367, 0.00354, 0.00368, 0.00354, 0.00357, 
                         0.00363, 0.00332, 0.00358, 0.00357, 0.00359, 0.00340, 
                         0.00350, 0.00332, 0.00342, 0.00347, 0.00342, 0.00290, 
                         0.00314 };

/* Spectral Shape of bbw Normalized to bbw(400 nm) */
float bbw_spectrum[] = { 1.0000, 0.8991, 0.8107, 0.7330, 0.6645, 0.6039, 
                         0.5500, 0.5021, 0.4591, 0.4210, 0.3866, 0.3557, 
			 0.3278, 0.3026, 0.2798, 0.2591, 0.2403, 0.2232, 
			 0.2075, 0.1933, 0.1802, 0.1682, 0.1572, 0.1471, 
			 0.1378, 0.1292, 0.1213, 0.1139, 0.1071, 0.1009, 
			 0.0950};

float AP_sat = 0.0; 	        /* Photons absorbed by phyto [PAR * ap/a]      */
float decl;			/* Solar declination                           */
float m0;                      	/* Kd coefficient (Lee et al. 2005)            */ 
float m1 = 4.18;		/* Kd coefficient (Lee et al. 2005)            */ 
float m2 = 0.52;		/* Kd coefficient (Lee et al. 2005)            */ 
float m3 = 10.8;		/* Kd coefficient (Lee et al. 2005)            */ 
float DL;			/* Dalength (hours)                            */
float solzen;			/* Solar zenith angle                          */
float kd[31];			/* Downwelling attenuation coefficient [m-1]   */		
float zeu;			/* Euphotic Depth [m]                          */		
float kdpar;                   	/* Downwelling attenuation coefficient of PAR  */
float tseq[51];			/* Time sequence                               */
float zseq[51];               	/* Depth sequence                              */
float delz;                    	/* Magnitude of depth sequence [m]             */
float PAR_noon[31];            	/* Maximum Daily PAR                           */
float AP_tzw[51][51][31];    	/* Absorbed photons(time, depth, wv)           */
float AP_tz[51][51] = {0.0}; 	/* Absorbed photons(time, depth)               */
float AP_t[51] = {0.0};       	/* Absorbed photons(depth)                     */
float AP = {0.0};              	/* Absorbed photons                            */
float Eu;                      	/* Fraction of Upwelling Irradiance            */
float IML;			/* Median Irradiance in the mixed layer        */
float Ek[51] = {0.0};         	/* Photoacclimation parameter                  */
float Eg_mld;                  	/* Growth irradiance in the mixed layer        */
float Eg;                      	/* Growth irradiance                           */
float phimax[51];                  /* photon yield of net carbon fixation    */
float phirange[2] = {0.018, 0.030}; /* bounds for phimax                      */
float Ekrange[2] = {12.96, 0.864};  /* EK correspoding to phimax              */
float E_tzw[51][51][31] = {0.0};  /* Irradiance(time, depth, wv)            */
float E_tz[51][51] = {0.0};       /* Irradiance(time, depth)                */
float mean_aphi = 0.0;              /* Spectrally average aphi                */
float numerator = 0.0;		    /* For Spectral Correction Factor         */
float denominator = 0.0;	    /* For Spectral Correction Factor         */
float KPUR[51] = {0.0};            /* Spectrally scaled EK                   */
float NPP_tz[51][51] = {0.0};     /* Net Primary Production(time, depth)    */
float NPP_t[51] = {0.0};           /* Net Primary Production(depth}          */
float NPP = 0.0;                    /* Net Primary Production [mol C m-2 day-1]*/
float NPP_mld = 0.0;                /* NPP in MLD             [mol C m-2 day-1]*/
float GPP_tz[51][51] = {0.0};     /* Gross Primary Production (time, depth)    */
float GPP_t[51] = {0.0};           /* Gross Primary Production (depth)          */
float GPP = 0.0;                    /* Gross Primary Production [mol photons m-2 day-1]*/
float GPP_mld = 0.0;                /* Gross Primary Production in the mld [mol photons m-2 day-1]*/
float AP_mld = 0.0;                /* Absorbed Photons in the mld [mol photons m-2 day-1]*/


/*------------------------------------------------------------------------------  
	Step 2: Derive IOPs at 10 nm increments from 400 to 700 nm                 
 
	Comments:
	IOP data from GIOP-DC (Werdell et al. 2013)
	GIOP-DC assumes slope of adg = 0.018 [m-1]
	GIOP-DC assumes spectral shape of phyto absorption coefficient is a function 
	of Chl (Bricuad et al. 1998)                                                */ 

for (w=0; w<31; w++){
  bbw[w]  = bbwater * bbw_spectrum[w];  
  aphi[w] = aph_443 * A_Bricaud[w] * pow(chl, E_Bricaud[w]) / 
            (0.03711 * pow(chl, 0.61479));
  a[w]    = aw[w] + aphi[w] + adg_443 * exp(-0.018 * (wv[w] - 443.0));
  bb[w]   = bbw[w] + bbp_443 * powf(443.0 / wv[w], bbp_s);
}

/* -----------------------------------------------------------------------------  
  Step 3: Calculate Absorbed Energy                                           */ 

for (w=0; w<30; w++){
  AP_sat += 5.0 * (PAR_spectrum[w] * aphi[w]/a[w] + 
                                  PAR_spectrum[w+1] * aphi[w+1]/a[w+1]);
}

AP_sat *= PAR * 0.95;

/* -----------------------------------------------------------------------------  
	Step 4: Derive Kd following Lee et al 2005, Eq. 11                     */
 
lat *= M_PI/180.0;
decl = 23.5 * cos (2.0 * M_PI * (yday - 172.0) / 365.0) * M_PI / 180.0;

DL = -1.0 * tan(lat) * tan(decl);

if (DL > 1.0)  {DL = 1.0;}   /* Check for daylength less than 0 hours */
if (DL < -1.0) {DL = -1.0;}  /* Check for daylength greater than 24 hours */  
    
DL = acos(DL) / M_PI; /* Daylength in days */

solzen = 90.0 - asin (sin(lat) * sin(decl) - 
         cos(lat) * cos(decl) *  cos(M_PI)) * 180.0 / M_PI;
m0     = sqrt((1.0 + 0.005 * solzen) * (1.0 + 0.005 * solzen)); 

for (w=0; w<31; w++){
  kd[w] =  m0 * a[w] + m1 * (1 - m2 * exp(-m3 * a[w])) * bb[w];
}

/*------------------------------------------------------------------------------  
	Step 5: Construct water column irradiance through time (t) and depth (z)    */

kdpar = 0.0665 + 0.874 * kd[9] - 0.00121 / kd[9];
zeu   = -1.0 * log (0.1 /(PAR * 0.95)) / kdpar;

for (t=0; t<51; t++){
  tseq[t] = (float)t / 50.0;
  zseq[t] = (float)t / 50.0 * ceil(zeu);
}

delz = zseq[1] - zseq[0];

for (w=0; w<31; w++){
  PAR_noon[w] = M_PI / 2.0 * PAR * 0.95  * PAR_spectrum[w];
}
  
for (t=0; t<51; t++){
  for (z=0; z<51; z++){
    for (w=0; w<31; w++){
      E_tzw[t][z][w]  = PAR_noon[w] * sin(M_PI*tseq[t]) * 
                        exp(-1.0 * kd[w] * zseq[z]);
      AP_tzw[t][z][w] = E_tzw[t][z][w] * aphi[w]; 
    }
  }
}

/* Integrate E and AP through wavelength                                      */
for (t=0; t<51; t++){
  for (z=0; z<51; z++){  
    for (w=0; w<30; w++){
      E_tz[t][z]  += 5.0 * (E_tzw[t][z][w] + E_tzw[t][z][w+1]);
      AP_tz[t][z] += 5.0 * (AP_tzw[t][z][w] + AP_tzw[t][z][w+1]);
    }
  }
}
  
for (t=0; t<51; t++){                      /* Integrate AP through depth     */
  for (z=0; z<50; z++){  
    AP_t[t] += delz * 0.5 * (AP_tz[t][z] + AP_tz[t][z+1]);
  }
}

for (t=0; t<50; t++){                      /* Integrate AP through time      */
	AP += 0.01 * (AP_t[t] + AP_t[t+1]);
}

/* Derive Upwelling Irradiance, absorbed energy is from Section 2 */
Eu = AP_sat / AP;  


/* -----------------------------------------------------------------------------  
	Step 6: CALCULATE EK through depth	                              */  

float foo;
float bar;

foo = PAR * 0.95 / (DL* 24.0);

IML   = foo * exp(-0.5 * kdpar * mld);  

for (z=0; z<51; z++){
  Ek[z] = 19.0 * exp(0.038 * powf(foo, 0.45) / kdpar);
}

if (mld < zeu){
  
  Eg_mld = (PAR / DL)  * exp(-1.0 * kdpar * mld);
  
  for (z=0; z<51; z++){
    Ek[z] *= (1.0 + exp(-0.15 * foo)) / (1.0 + exp(-3.0 * IML));
	if (zseq[z] > mld){
		bar = Ek[z];
      		Eg = (PAR / DL)  * exp(-1.0 * kdpar * zseq[z]);  
      		Ek[z] = 10.0 + (bar - 10.0) / (Eg_mld - 0.1) * (Eg - 0.1);
	}
    if (Ek[z] < 10.0) {Ek[z] = 10.0;}       
  }
}

for (z=0; z<51; z++){
  Ek[z] *= 0.0864;                    /* 0.0864 Convert to mol photons/m2/day */
}


/* -----------------------------------------------------------------------------  
	Step 7: Calculate phimax as a function of Ek 			      */    

for (z=0; z<51; z++){
  
  phimax[z] = phirange[1] + (Ek[z] - Ekrange[1]) * (-0.000992);
  
  if (phimax[z] > phirange[1]) {phimax[z] = phirange[1];}
  if (phimax[z] < phirange[0]) {phimax[z] = phirange[0];}
}

/*------------------------------------------------------------------------------  
	Step 8: Derive Spectral Correction Factor and apply to Ek  	      */

for (w=0; w<31; w++){
  mean_aphi += aphi[w];
}
mean_aphi /=  31.0;

for (z=0; z<51; z++){
	numerator = 0.0;
	denominator = 0.0;
  for (w=0; w<30; w++){   
    numerator   += 5.0 * (AP_tzw[25][z][w] + AP_tzw[25][z][w+1]);
    denominator += 5.0 * (E_tzw[25][z][w] + E_tzw[25][z][w+1]) * mean_aphi;
  }  
  KPUR[z] = Ek[z] / (numerator / denominator / 1.3);
}


/* -----------------------------------------------------------------------------  
	Step 9: Scale phytoplankton beneath the MLD to Ek 
           Need to rescale E and absorbed energy beneath the MLD */

float aphi_fact[51] = {1.0};

for (z=0; z<51; z++){
  if (zseq[z] > mld){
    aphi_fact[z] = 1.0 + Ek[0]/Ek[z] * 0.15;
  } else {
	  aphi_fact[z] = 1.0; 
	}
}

/* -----------------------------------------------------------------------------  
 Step 10: Reconstruct water column irradiance     			*/

/* Modify Irradiance and absorbed energy to account for upwelled irradiance   */

for (t=0; t<51; t++){
  for (w=0; w<31; w++){
    E_tzw[t][0][w]  = PAR_noon[w] * sin(M_PI*tseq[t]);
    AP_tzw[t][0][w] = E_tzw[t][0][w] * aphi[w] * aphi_fact[0];
  } 
  for (z=1; z<51; z++){
    for (w=0; w<31; w++){
      a[w]            = aw[w] + aphi[w] * aphi_fact[z] + adg_443 * 
                        exp(-0.018 * (wv[w] - 443.0));
      kd[w]           = m0 * a[w] + m1 * (1.0 - m2 * exp(-m3 * a[w])) * bb[w];
      E_tzw[t][z][w]  = E_tzw[t][z-1][w] * exp(-1.0 * kd[w] * delz);
      AP_tzw[t][z][w] = E_tzw[t][z][w] * aphi[w] * aphi_fact[z];
     }
  }
}

/* Integrate E and AP through wavelength and multiply by Eu                   */
for (t=0; t<51; t++){
  for (z=0; z<51; z++){
    E_tz[t][z] = 0.0;
    AP_tz[t][z] = 0.0;
    for (w=0; w<30; w++){
      E_tz[t][z] += 5.0 * (E_tzw[t][z][w] + E_tzw[t][z][w+1]);
      AP_tz[t][z] += 5.0 * (AP_tzw[t][z][w] + AP_tzw[t][z][w+1]);
    }
  E_tz[t][z]  *= Eu; 
  AP_tz[t][z] *= Eu;
  }
}

/* ----------------------------------------------------------------------------- 
 Step 11: CALCULATE NPP, GPP and absorbed photons                              */ 

            
for (t=0; t<51; t++){
  for (z=0; z<51; z++){
    GPP_tz[t][z] = AP_tz[t][z] * tanh(KPUR[z] / E_tz[t][z]); 
    NPP_tz[t][z] = GPP_tz[t][z] * phimax[z];  
  }
}    

/* Integrate through depth */
for (t=0; t<51; t++){
  AP_t[t] = 0.0;
  for (z=0; z<50; z++){  
    AP_t[t] += delz * 0.5 * (AP_tz[t][z] + AP_tz[t][z+1]);
    GPP_t[t] += delz * 0.5 * (GPP_tz[t][z] + GPP_tz[t][z+1]);
    NPP_t[t] += delz * 0.5 * (NPP_tz[t][z] + NPP_tz[t][z+1]);
    
  }
}

/* Integrate through time */
AP = 0.0;
for (t=0; t<50; t++){  
  fnps += 0.01 * (fnps_t[t] + fnps_t[t+1]);
  GPP += 0.01 * (GPP_t[t] + GPP_t[t+1]);
  NPP += 0.01 * (NPP_t[t] + NPP_t[t+1]);
  AP += 0.01 * (AP_t[t] + AP_t[t+1]);
}

/* Integrate over MLD */
if (zseq[50] < mld){
	NPP_mld = NPP;
	GPP_mld = GPP;
	AP_mld = AP;
}else{
	for (t=0; t<51; t++){
		fnps_t[t] = 0.0;	
		AP_t[t] = 0.0;	
		GPP_t[t] = 0.0;
		NPP_t[t] = 0.0;		
		z = 1;
  		while (zseq[z] < mld) {
			AP_t[t] += delz * 0.5 * (AP_tz[t][z] + AP_tz[t][z+1]);
			GPP_t[t] += delz * 0.5* (GPP_tz[t][z] + GPP_tz[t][z+1]);
    			NPP_t[t] += delz * 0.5 * (NPP_tz[t][z] + NPP_tz[t][z+1]);
			z++;
  		}
	}
	for (t=0; t<50; t++){  
		AP_mld += 0.01 * (AP_t[t] + AP_t[t+1]);
		GPP_mld += 0.01 * (GPP_t[t] + GPP_t[t+1]);
		NPP_mld += 0.01 * (NPP_t[t] + NPP_t[t+1]);
	}
}

/* ----------------------------------------------------------------------------- 
 Step 12: Return Data                                                         */ 

float *results = (float *)malloc(sizeof(float)*6);

results[0] = NPP;
results[1] = GPP;
results[2] = AP;
results[3] = NPP_mld;
results[4] = GPP_mld;
results[5] = AP_mld;

return results;

free(results);

}
