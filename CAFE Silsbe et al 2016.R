# This code is complete but inefficient. 
# Please note that single pixel data is passed to  the CAFE function 
# (i.e chl is not a global array, rather it is the chl at a single point). 
# Therefore a simple wrapper outside of this code is used to incrementally 
# pass data to the function.



###########################################################################################
#                                    List of Functions

#CAFE: Main code

#trap.t, trap.w, trap.z: Trapeziod integration functions (time, wavelength, depth)

#backscattering of pure water functions:
#   betaws_ZHH2009, density_sw, dlnasw_ds
#   IsoComp, PMH, RInw, 



###########################################################################################

CAFE <- function(
          aph_443, # Phytoplankton absorption coefficient [m-1]
          adg_443, # Absorption coefficient of detritus and gelbstoff [m-1]
          bbp_443, # Particulate backscaterring coefficient [m-1]
          bbp_s,   # Slope of backscattering coefficient [dimensionless]
          PAR,     # Daily integrated PAR [mol photons m-2 day-1]
          chl,     # Cholorphyll a concentration [mg m-3]
          sst,     # Sea surface temperature [Deg C]
          sal,     # Salinity of water [PSU]
          mld,     # Mixed layer Depth [m]
          lat,     # Latitude, nothern hemisphere is postive, decimal degrees     
          yd){     # Day of year

  #Notes on Input
  # 1) If Mixed layer depth is not known, can pass deep value (>250 m) to assume
  #    vertically homogenous water column       
  # 2) IOP data from GIOP-DC (Werdell et al. 2013)
  #    - GIOP-DC assumes slope of adg = 0.018 [m-1]
  #    - GIOP-DC assumes spectral shape of phyto absorption coefficient 
  #      is a function of Chl (Bricuad et al. 1998)

  
  #########################################################################################
  # Step 1: Derive IOPs at 10 nm increments from 400 to 700 nm
  
  #Variables
  #wv: wavelength [nm]
  #aw: Absorption of Pure Water (Pope and Fry 1997) 
  #aphi: Phytoplantkon absorption coefficient
  #A.Bricaud, E.Bricaud: Spectral shape of phyto. abs. coeff (Bricaud et al. 1998)
  #bbw:  Pure water backscattering coefficients derived from Zhang et al. (2009) 
  
  #Comments
  #  GIOP-DC assumes slope of adg = 0.018 [m-1]
  #  GIOP-DC assumes spectral shape of aphi is a function of Chl (Bricuad et al. 1998)
  #########################################################################################
  
  wv   <- seq(400, 700, by=10) 
  
  aw   <- c(0.00663, 0.00473, 0.00454, 0.00495, 0.00635, 0.00922, 0.00979, 0.0106, 0.0127, 
            0.015, 0.0204, 0.0325, 0.0409, 0.0434, 0.0474, 0.0565, 0.0619, 0.0695, 0.0896,
            0.1351, 0.2224,0.2644, 0.2755, 0.2916, 0.3108, 0.34, 0.41, 0.439, 0.465, 0.516,
            0.624)  
  
  A.Bricaud <- c(0.0241, 0.0287, 0.0328, 0.0359, 0.0378, 0.0350, 0.0328, 0.0309, 0.0281,
                 0.0254, 0.0210, 0.0162, 0.0126, 0.0103, 0.0085, 0.0070, 0.0057, 0.0050, 
                 0.0051, 0.0054, 0.0052, 0.0055, 0.0061, 0.0066, 0.0071, 0.0078, 0.0108, 
                 0.0174, 0.0161, 0.0069, 0.0025)
  
  E.Bricaud <- c(0.6877, 0.6834, 0.6664, 0.6478, 0.6266, 0.5993, 0.5961, 0.5970, 0.5890, 
                 0.6074, 0.6529, 0.7212, 0.7939, 0.8500, 0.9036, 0.9312, 0.9345, 0.9298, 
                 0.8933, 0.8589, 0.8410, 0.8548, 0.8704, 0.8638, 0.8524, 0.8155, 0.8233, 
                 0.8138, 0.8284, 0.9255, 1.0286)

  if (is.finite(aph_443)){
    aphi <- aph_443 *  A.Bricaud * chl^(E.Bricaud) / (0.03711 * chl^(0.61479))
  }else{ 
    aphi <- A.Bricaud * chl^(E.Bricaud)
  }
  

  a    <- aw + aphi + adg_443 * exp( -0.018 * (wv - 443))
  
  bbw <- rep(NA, 31)
  for (w in 1:31){
    bbw[w] <- betasw_ZHH2009(wv[w], sal, sst, delta= 0.039)/2} 
  
  bb   <- bbw + bbp_443 * (443 / wv)^(bbp_s)
  
  #########################################################################################
  # Step 2: Calculate the fraction of photons absorbed by phytoplankton
  
  # Assume %5 is reflected/upwelled from surface
  
  # PARfraction derived from ASTMG173 reference spectrum 
  # http://rredc.nrel.gov/solar/spectra/am1.5/astmg173/astmg173.html
  #########################################################################################
  
  PARfraction <- c(0.00227, 0.00218, 0.00239, 0.00189, 0.00297, 0.00348, 0.00345, 0.00344,
                   0.00373, 0.00377, 0.00362, 0.00364, 0.00360, 0.00367, 0.00354, 0.00368, 
                   0.00354, 0.00357, 0.00363, 0.00332, 0.00358, 0.00357, 0.00359, 0.00340, 
                   0.00350, 0.00332, 0.00342, 0.00347, 0.00342, 0.00290, 0.00314)
  
  absorbed_photons <-  PAR * 0.95 * trap.w(PARfraction * aphi/a)
  
  #########################################################################################
  # Step 3: Derive Kd following Lee et al 2005, Eq. 11

  #dec: declination of sun
  #DL: Daylength calculated following Westberry et al 2008 
  #solzen: Solar zenith angle
  #m0, m1, m2, m3: Coefficients from Lee et al. 2005
  #########################################################################################
  
  lat <- lat * pi/180   
  dec <- 23.5 * cos(2 * pi * (yd - 172)/365) *pi/180 
  
  DL           <- -1 * tan(lat) * tan(dec)
  DL[DL>1]     <- 1                         #Check for daylengths less than 0 hours
  DL[DL<(-1)]  <- -1                        #Check for daylengths greater than 24 hours
  DL           <- acos(DL) / pi             #Daylength in days
  
  solzen <- 90 - asin (sin(lat) * sin(dec) - cos(lat) * cos(dec) *  cos(pi)) *180/pi
  m0     <- abs(1+0.005*solzen)
  m1     <- 4.18
  m2     <- 0.52
  m3     <- 10.8
  kd     <-  m0 * a + m1 * (1 - m2 * exp(-m3 * a)) * bb
  
  #########################################################################################
  # Step 4: Construct water column irradiance through time (t) and depth (z)
  
  #tseq:   Divides diurnal period into 101 increments [dimensionless]
  #kdpar:  Attenuation coefficient of downwelling PAR [m-1]
  #        Derived from Kd(490) following Morel et al. 2007. 
  #zeu:    Euphotic depth [m]. Taken as 0.1 mol photons m-2 day-1 isolume
  #zseq:   Divides euphotic depth into 101 increments [m]
  #delz:   Incremental change  in depth [m]
  #E.tz:   Array of irradiance through time and depth [mol photons m-2 d-1]
  #AP.tz:  Array of absorbed photons through time and depth [mol photons m-3 d-1]
  #########################################################################################
 
  tseq   <- seq(0, 1, length.out=101)
  kdpar  <- 0.0665 + 0.874 * kd[10]  - 0.00121 / kd[10] 
  zeu    <- -1 * log (0.1 /(PAR * 0.95)) / kdpar
  zseq   <- seq(0, ceiling(zeu), length.out=101)
  delz   <- zseq[2] - zseq[1]
  
  #Setup arrays
  E.tz   <- array(0, c(101, 101)) #Irradiance
  AP.tz  <- array(0, c(101, 101)) #Absorbed photons

  #Max surface PAR at solar noon
  PAR.noon <- pi * PAR * 0.95 / 2 * PARfraction
  
  for (t in 1:101){
    for (z in 1:101){
      E.tz[t,z]   <- trap.w(PAR.noon * sin(pi * tseq[t]) * exp(-kd * zseq[z]))
      AP.tz[t,z]  <- trap.w(PAR.noon * sin(pi * tseq[t]) * exp(-kd * zseq[z]) * aphi)
    }
  }
  
  #Integrate AP through time and depth 
  AP.z  <- apply(AP.tz, 2, trap.t) 
  AP    <- trap.z(AP.z, delz, 101)
 
  
  # Try C implementation
  AP.z_C = rep(0, 101)
  for (z in 1:101){
    for (t in 1:100){
      AP.z_C[z] = 0.005*(AP.tz[t+1][z] + AP.tz[t][z])
    }
  }
  cbind(AP.z, AP.z_C)
  
  
   
  #Derive Upwelling Irradiance, absorbed energy is from Section 2
  Eu    <- absorbed_photons / AP  
  
  #Modify Irradiance and absorbed energy to account for upwelled irradiance
  E.tz  <-  E.tz * Eu
  AP.tz <-  AP.tz * Eu

  
  #########################################################################################
  # Step 5: CALCULATE EK through depth
  
  #Surface Ek is calculated following: 
  #I extended this algorithm to propagate Ek beneath the MLD                    
  
  #IML: Median Mixed Layer Irradiance [mol photons m-2 hour-1]
  
  #########################################################################################

  IML   <- (PAR * 0.95 / (DL* 24)) * exp(-0.5 * kdpar * mld)  
  Ek    <- rep(19 * exp(0.038 * (PAR * 0.95 / (DL * 24)) ^ 0.45 /kdpar), 101) 
  
  if (mld < zeu){
    
    Ek <- Ek * (1 + exp(-0.15 * (PAR * 0.95 / (DL * 24)))) / (1 + exp(-3 * IML))
    
    #Find indices of zseq deeper than MLD
    deep     <- which(zseq > mld)
    Eg       <- (PAR / DL)  * exp(-1*kdpar*zseq)  
    Eg.mld   <- (PAR / DL)  * exp(-1*kdpar*mld)  
    Ek[deep] <- 10 + (Ek[deep] - 10) / (Eg.mld - 0.1) * (Eg[deep] -0.1)
  }
  
  Ek[Ek < 10] <- 10          #Ensure Ek is no smaller than 10 umol m-2 s-1
  Ek          <- Ek * 0.0864 #Convert to mol photons/m2/day

  #########################################################################################
  # Step 6: CALCULATE SCF (Spectral Correction Factor)
  # KPUR is spectrally scaled EK
  #########################################################################################
  SCF    <- rep(NA, 101)
  
  for (z in 1:101){
    E.tzw    <- PAR.noon * sin(pi * tseq[50]) * exp(-kd * zseq[z]) 
    AQ.tzw   <- E.tzw * aphi      
    SCF[z]   <- trap.w(AQ.tzw) / (trap.w(E.tzw) * mean(aphi)) / 1.3
  }
  
  KPUR <- Ek / SCF  

  #########################################################################################
  # Step 7: Tie PHIMax to Ek
  #########################################################################################
  
  phirange                     <- c(0.018, 0.030)
  ekrange                      <- c(150*86400/1e6, 10*86400/1e6)
  slope                        <- (phirange[2] - phirange[1]) / (ekrange[2] - ekrange[1])
  phimax                       <- phirange[2] + (Ek - ekrange[2]) * slope
  phimax[phimax < phirange[1]] <- phirange[1]
  phimax[phimax > phirange[2]] <- phirange[2]
  
  #########################################################################################
  # Step 8: Deep Chlorophyll Maxima is scaled to Ek
  #########################################################################################
  
  if (mld < zeu){
    
    aphi.fact       <- rep(1, 101) 
    aphi.fact[deep] <- 1 + Ek[1]/Ek[deep] * 0.15
    
    #Recalculate irradiance and absorbed energy over time depth array
    aphi.z <- matrix(rep(aphi,101), nrow=101, byrow=T) * aphi.fact
    a.z    <- matrix(rep(aw, 101), nrow=101, byrow=T) +
      matrix(rep(adg_443 * exp(-0.018 * (wv-443)), 101), nrow=101, byrow=T) +
      aphi.z
    mat    <- m1 *(1- m2*exp(-m3*a.z))
    kd.z   <- m0*a.z + sweep(mat, 2, bb, '*')
    
    for (t in 1:101){
      E.wv       <- PAR.noon * sin(pi * tseq[t])
      AP.tz[t,1] <- trap.w(E.wv * aphi.z[1,])
      
      for (z in 2:101){
        E.wv       <- E.wv * exp(-kd.z[z,] * (zseq[2] - zseq[1]))
        AP.tz[t,z] <- trap.w(E.wv * aphi.z[z,])
      }
    }
    
    #Modify absorbed energy to account for upwelled irradiance
    AP.tz <-  AP.tz * Eu
  }else{
    aphi.z <- matrix(rep(aphi,101), nrow=101, byrow=T)
  }


  
  #########################################################################################
  # Step 9: Final Calculations
  
  #Calculate carbon specific growth Rate
  # Phytoplankton carbon biomass derived from bbp(470 nm) following
  #Graff, J.R., T.K. Westberry, A.J. Milligan, M.B. Brown, G. Dall'Olmo, 
  #V. van Dongen-Vogels, K.M. Reifel, and M.J. Behrenfeld. 2015. Analytical phytoplankton
  #carbon measuremetns spanning diverse ecosystems. Deeo-Sea Research I. 102: 16-25.
  
  
  #########################################################################################

  #Net Primary Production
  GPP.tz <-  t(t(AP.tz * tanh(t(KPUR/t(E.tz))))) * 12000
  GPP.z  <- apply(GPP.tz, 2, trap.t)
  GPP    <- trap.z(GPP.z, delz, 101)
    
  #Net Primary Production
  NPP.tz <-  t(phimax * t(AP.tz * tanh(t(KPUR/t(E.tz))))) * 12000
  NPP.z  <- apply(NPP.tz, 2, trap.t)
  NPP    <- trap.z(NPP.z, delz, 101)
  
  #Absorbed photons
  AP.z   <- apply(AP.tz, 2, trap.t) 
  AP     <- trap.z(AP.z, delz, 101)
  
  
  #Fraction of absorbed photons dissipated
  NPQ      <- AP.tz * phimax * 12000
  NPQ.z    <- apply(NPQ, 2, trap.t)
  NPQ      <- 1 - NPP/trap.z(NPQ.z, delz, 101)
  
  #Calculate NPP, AQ, and NPQ and mu in the mixed layer

  
  if (mld > max(zseq)){   #MLD is deeper than euphotic depth
    NPP.mld <- NPP
    AP.mld  <- AP     
    NPQ.mld <- NPQ
  }else{                 #MLD is shallower than euphotic depth (i.e. vertical structure)
    
    ind     <- which.min(abs(mld-zseq))
    AP.mld  <- trap.z(AP.z, delz, zseq[ind])
    NPP.mld <- trap.z(NPP.z, delz, zseq[ind])
  
  }
  

  
  bbp.470 <- bbp_443*(443/470)^(bbp_s)
  Cphyto  <- 12128 * bbp.470 + 0.59
  mu.mld  <- NPP.mld / (Cphyto * mld)
  mu      <- NPP / (Cphyto * zeu)
  
  #####################################################
  # Return Data
  #####################################################

  return(list(AP=AP,                       # photons absorbed by phytoplankton [mol photons m-2 d-1] 
              AP.mld=AP.mld,               # photons absorbed by phytoplankton in the Mixed layer [mol photons m-2 d-1]
              AP.sat = absorbed_photons,   # photons absorbed by phytoplankton (Satellite Derived) [mol photons m-2 d-1] (equivalent to AP when mld > zeu, otherwise lower due to DCM)
              NPP=NPP,                     # Net Primary Production [mg C m-2 d-1]
              NPP.mld=NPP.mld,             # Net Primary Production in the mixed layer [mg C m-2 d-1]
              GPP=GPP,                     # Gross Primary Production [mol photons C m-2 d-1]
              GPP.mld=GPP.mld,             # Gross Primary Production in the mixed layer [mol photons m-2 d-1]
              mu=mu,                       # Effective Growth Rate through euphotic depth [d-1]
              mu.mld = mu.mld))            # Growth rate in the mixed layer [m-1]          
          
              
}

#Trapezoid Integration Functions

trap.t <- function(f){return(0.5 * sum(0.01*(f[2:101] + f[1:100])))}

trap.w <- function(f){return(0.5 * sum(10 * (f[2:31] + f[1:30])))}

trap.z <- function(f, delz, zlim) {return(0.5 * sum(delz * (f[2:zlim] + f[1:(zlim-1)])))}


betasw_ZHH2009 <- function(lambda, S, Tc, delta= 0.039){
  
  #function [theta,betasw,bsw,beta90sw]= betasw_ZHH2009(lambda,S,Tc,delta)
  # Scatteirng by pure seawater: Effect of salinity
  # Xiaodong Zhang, Lianbo Hu, and Ming-Xia He, Optics Express, 2009, accepted
  # lambda (nm): wavelength
  # Tc: temperauter in degree Celsius, must be a scalar
  # S: salinity, must be scalar
  # delta: depolarization ratio, if not provided, default = 0.039 will be
  # used.
  # betasw: volume scattering at angles defined by theta. Its size is [x y],
  # where x is the number of angles (x = length(theta)) and y is the number
  # of wavelengths in lambda (y = length(lambda))
  # beta90sw: volume scattering at 90 degree. Its size is [1 y]
  # bw: total scattering coefficient. Its size is [1 y]
  # for backscattering coefficients, divide total scattering by 2
  #
  # Xiaodong Zhang, March 10, 2009
  
  # values of the constants
  Na  = 6.0221417930e23 ;   #  Avogadro's constant
  Kbz = 1.3806503e-23 ;     #  Boltzmann constant
  Tk  = Tc + 273.15 ;       #  Absolute tempearture
  M0  = 18e-3;              #  Molecular weigth of water in kg/mol
  
  theta= seq(0,180, by=1)
  
  rad = theta * pi/180; # angle in radian as a colum variable
  
  # nsw: absolute refractive index of seawater
  # dnds: partial derivative of seawater refractive index w.r.t. salinity
  nsw  = RInw(lambda,Tc,S)$nsw
  dnds = RInw(lambda,Tc,S)$dnswds
  
  # isothermal compressibility is from Lepple & Millero (1971,Deep
  # Sea-Research), pages 10-11
  # The error ~ ?0.004e-6 bar-1
  IsoCom = IsoComp(Tc,S)
  
  density_sw = density_sw(Tc, S);
  
  # water activity data of seawater is from Millero and Leung (1976,American
  # Journal of Science,276,1035-1077). Table 19 was reproduced using
  # Eq.(14,22,23,88,107) then were fitted to polynominal equation.
  # dlnawds is partial derivative of natural logarithm of water activity
  # w.r.t.salinity
  dlnawds = dlnasw_ds(Tc, S);
  
  # density derivative of refractive index from PMH model
  DFRI = PMH(nsw);  ## PMH model
  
  # volume scattering at 90 degree due to the density fluctuation
  beta_df = pi*pi/2*((lambda*1e-9)^(-4)) * Kbz * Tk * IsoCom * DFRI^2 * (6 + 6*delta) / (6-7*delta)
  # volume scattering at 90 degree due to the concentration fluctuation
  flu_con = S * M0 * dnds^2 / density_sw / (-1*dlnawds)/Na
  beta_cf = 2*pi*pi*((lambda*1e-9)^(-4)) * nsw^2 * (flu_con) * (6+6*delta) / (6-7*delta)
  # total volume scattering at 90 degree
  beta90sw = beta_df+beta_cf
  bsw      = 8 * pi/3 * beta90sw * (2+delta) / (1+delta)
  return(bsw)
  
}


RInw <- function(lambda,Tc,S){
  # refractive index of air is from Ciddor (1996,Applied Optics)
  n_air = 1.0 + (5792105.0 / (238.0185 - 1 / (lambda / 1e3)^2) + 167917.0 / (57.362 - 1 / (lambda/1e3)^2)) /1e8
  
  # refractive index of seawater is from Quan and Fry (1994, Applied Optics)
  n0 = 1.31405; n1 = 1.779e-4 ; n2 = -1.05e-6 ; n3 = 1.6e-8 ; n4 = -2.02e-6 ;
  n5 = 15.868; n6 = 0.01155;  n7 = -0.00423;  n8 = -4382 ; n9 = 1.1455e6;
  
  nsw = n0 + (n1 + n2 * Tc + n3 * Tc^2) * S + n4 * Tc^2 + (n5+n6*S+n7*Tc)/lambda+n8/lambda^2+n9/lambda^3 # pure seawater
  nsw = nsw * n_air
  dnswds = (n1+n2*Tc+n3*Tc^2+n6/lambda)*n_air
  return(list(nsw=nsw, dnswds=dnswds))
}

IsoComp <- function(Tc, S){
  # pure water secant bulk Millero (1980, Deep-sea Research)
  kw = 19652.21 + 148.4206 * Tc - 2.327105 * Tc^2 + 1.360477e-2 * Tc^3 - 5.155288e-5 * Tc^4
  Btw_cal = 1/kw
  
  # isothermal compressibility from Kell sound measurement in pure water
  # Btw = (50.88630+0.717582*Tc+0.7819867e-3*Tc.^2+31.62214e-6*Tc.^3-0.1323594e-6*Tc.^4+0.634575e-9*Tc.^5)./(1+21.65928e-3*Tc)*1e-6;
  
  # seawater secant bulk
  a0 = 54.6746 - 0.603459 * Tc + 1.09987e-2 * Tc^2 - 6.167e-5 * Tc^3
  b0 = 7.944e-2 + 1.6483e-2 * Tc - 5.3009e-4 * Tc^2
  
  Ks = kw + a0*S + b0*S^1.5
  
  # calculate seawater isothermal compressibility from the secant bulk
  return(1/Ks*1e-5) # unit is pa
}

density_sw <- function(Tc, S){
  
  # density of water and seawater,unit is Kg/m3, from UNESCO,38,1981
  a0 = 8.24493e-1;  a1 = -4.0899e-3; a2 = 7.6438e-5; a3 = -8.2467e-7; a4 = 5.3875e-9;
  a5 = -5.72466e-3; a6 = 1.0227e-4;  a7 = -1.6546e-6; a8 = 4.8314e-4;
  b0 = 999.842594; b1 = 6.793952e-2; b2 = -9.09529e-3; b3 = 1.001685e-4;
  b4 = -1.120083e-6; b5 = 6.536332e-9;
  
  # density for pure water 
  density_w = b0+b1*Tc+b2*Tc^2+b3*Tc^3+b4*Tc^4+b5*Tc^5
  # density for pure seawater
  return(density_w +((a0+a1*Tc+a2*Tc^2+a3*Tc^3+a4*Tc^4)*S+(a5+a6*Tc+a7*Tc^2)*S^1.5+a8*S^2))
}

dlnasw_ds <- function(Tc, S){
  # water activity data of seawater is from Millero and Leung (1976,American
  # Journal of Science,276,1035-1077). Table 19 was reproduced using
  # Eqs.(14,22,23,88,107) then were fitted to polynominal equation.
  # dlnawds is partial derivative of natural logarithm of water activity
  # w.r.t.salinity
  # lnaw = (-1.64555e-6-1.34779e-7*Tc+1.85392e-9*Tc.^2-1.40702e-11*Tc.^3)+......
  #            (-5.58651e-4+2.40452e-7*Tc-3.12165e-9*Tc.^2+2.40808e-11*Tc.^3).*S+......
  #            (1.79613e-5-9.9422e-8*Tc+2.08919e-9*Tc.^2-1.39872e-11*Tc.^3).*S.^1.5+......
  #            (-2.31065e-6-1.37674e-9*Tc-1.93316e-11*Tc.^2).*S.^2;
  
  return((-5.58651e-4 + 2.40452e-7 * Tc -3.12165e-9 * Tc^2 + 2.40808e-11 * Tc^3) +
           1.5*(1.79613e-5 -9.9422e-8 *Tc +2.08919e-9 * Tc^2 - 1.39872e-11 * Tc^3)* S^0.5 +
           2*(-2.31065e-6-1.37674e-9*Tc-1.93316e-11*Tc^2)*S)
}

# density derivative of refractive index from PMH model
PMH <- function(n_wat){
  n_wat2 = n_wat^2
  return( (n_wat2-1) * (1+2/3*(n_wat2+2 )* (n_wat/3 - 1/3/n_wat)^2) )
}

#########################################################################################
#References
#########################################################################################

#Behrenfeld, M.J.et al. 2016. Revaluating ocean warming impacts on global phytoplankton. 
#Nature Climate Change. 6: 323-330.

#Bricaud et al. 1998. Variation in light absorption by suspended particles with chlorophyll
#concentration in oceanic (case 1) waters: analysis and implications for bio-optical models
#J.Geophys. Res. 103: 31033-31044.

#Lee, Z et al. 2005. A model for the diffuse attenuation coefficient of 
#downwelling irradiance. J. Geophys. Res. 110: C02016. 

# Morel, A. et. al 2007. Examining the consitency of products derived from various ocean 
# color sensors in open ocean (Case 1) water in the perspective of a multi-sensor approach. 
# Remote Sensing of Environment. 111: 69-88