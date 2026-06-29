*-------------------------------------------------------------------------*
* Copyright 2013. Wageningen University, Plant Production Systems group,  *
* P.O. Box 430, 6700 AK Wageningen, The Netherlands.                      *
* You may not use this work except in compliance with the Licence.        *
* You may obtain a copy of the Licence at:                                *
*                                                                         *
* http://models.pps.wur.nl/content/licence-agreement                      *
*                                                                         *
* Unless required by applicable law or agreed to in writing, software     *
* distributed under the Licence is distributed on an "AS IS" basis,       *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
*-------------------------------------------------------------------------*

DEFINE_CALL GLA   (INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
                   INPUT,           OUTPUT)
DEFINE_CALL SUBEAI(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,  OUTPUT)
DEFINE_CALL TOTASS(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
                   OUTPUT,OUTPUT,OUTPUT)
DEFINE_CALL SUBGRT(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,...
                   INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,  OUTPUT)
DEFINE_CALL SUBFR (INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,  OUTPUT)


TITLE Crop growth for water-limited production (SUCROS2)
*     Spring wheat, Version September 1997 (SUCROS2_97 V1.0)       

* 2.1 INTRODUCTION

***********************************************************************
*    This model for potential crop growth (example: spring wheat)     * 
*    is described in:                                                 *
*    J. Goudriaan & H.H. van Laar, 1994. Modelling Potential Crop     *
*       Growth Processes. Textbook with Exercises. Kluwer Academic    *
*       Publishers, Dordrecht, The Netherlands, 238 pp.               *

*    Earlier versions of the model are described in:                  *
*    H. van Keulen, F.W.T. Penning de Vries & E.M. Drees, 1982.       *
*       A summary model for crop growth. In: Simulation of plant      *
*       growth and crop production. Eds F.W.T. Penning de Vries &     *
*       H.H. van Laar, Simulation Monographs, Pudoc, Wageningen,      *
*       pp. 87-97.                                                    *
*    Spitters, C.J.T., H. van Keulen & D.W.G. van Kraalingen, 1989.   *
*       A simple and universal crop growth model: SUCROS87. In:       *
*       Simulation and systems management in crop protection.         *
*       Eds R. Rabbinge, S.A. Ward & H.H. van Laar, Simulation        *
*       Monographs, Pudoc, Wageningen, pp.147-181.                    *   
*    The CSMP version of this model is described in:                  *
*    H. van Keulen, J. Goudriaan, L. Stroosnijder, E.A. Lantinga &    *
*    H.H. van Laar, 1992.                                             *
*       Crop growth model for water-limited production (SUCROS2).     *
*       In: Eds H.H. van Laar, J. Goudriaan & H. van Keulen,          *
*       Simulation of crop growth for potential and water-limited     *
*       production situations (as applied to wheat).                  *
*       Simulation Reports CABO-TT, December 1992, 78 pp.             * 
*       Department of Theoretical Production Ecology, Wageningen      *
*       Agricultural University, and DLO-Centre for Agrobiological    *
*       Research, Wageningen.                                         *
***********************************************************************

* 2.2 Initial conditions

INITIAL
INCON ZERO       = 0.
PARAMETER DOYEM  = 90.
INCON WLVI       = 0.5;  WSTI = 0.3; WRTI = 0.8
INCON WLVDI      = 0.;   WSOI = 0. ; ILAI = 0.012
INCON IDVS       = 0.;   IEAI = 0.

INCON ZRTI       = 5.
      WL1I       = WCLI1 * TKL1
      WL2I       = WCLI2 * TKL2
      WL3I       = WCLI3 * TKL3
      WL4I       = WCLI4 * TKL4
PARAM WCLI1      = 0.2;  WCLI2 = 0.2;  WCLI3 = 0.2;  WCLI4 = 0.2
PARAM TKL1       = 200.; TKL2  = 400.; TKL3  = 600.; TKL4  = 800.
      TKLT       = TKL1 + TKL2 + TKL3 + TKL4
      WCUMI      = WL1I + WL2I + WL3I + WL4I

INCON IDSLR      = 1.
PARAM MDRATE     = 50.

      ZRTM       = MIN(ZRTMC, ZRTMS, TKLT)
PARAM ZRTMS      = 1200.
PARAM ZRTMC      = 1200.


* Initialization of TNASS: total CO2 equivalents initially available
      TNASSI     = (WLVI*CFLV + WSTI*CFST + WRTI*CFRT) * 44./12.

* 2.3 Crop development

DYNAMIC
      DVS        = INTGRL(IDVS, DVR)
      DVR        = INSW(DVS-1., AFGEN(DVRVT, DAVTMP),...
                                AFGEN(DVRRT, DAVTMP)) * EMERG
      EMERG      = INSW(TIME-DOYEM, 0., 1.) 
FUNCTION DVRVT   = -10.,0., 0.,0., 30.,0.027
FUNCTION DVRRT   = -10.,0., 0.,0., 30.,0.031


* 2.4 Leaf CO2 assimilation

      AMAX       = AMX * AMDVS * AMTMP * EMERG
      AMDVS      = AFGEN(AMDVST, DVS)
      AMTMP      = AFGEN(AMTMPT, DDTMP)
PARAMETER AMX    = 1.11E-3
* 1.11 milligram CO2/m2/s  = 40 kg CO2/ha/h
FUNCTION AMDVST  = 0.0,1.0, 1.0,1.0, 2.0,0.5, 2.5,0.0
FUNCTION AMTMPT  = -10.,0., 0.,0., 10.,1., 25.,1., 35.,0., 50.,0.


* 2.5 Daily gross CO2 assimilation

     CALL TOTASS(DOY,LAT,DTR,SCP,AMAX,EFF,KDF,TAI,  ...
                 DAYL,DTGA,DS0)
PARAMETER EFF    = 12.5E-6
* 12.5 microgram CO2/J  = 0.45 kg (CO2/ha/h)/(J/m2/s) 
PARAMETER KDF    = 0.60
PARAMETER SCP    = 0.20
*PARAMETER LATT   = 52.


* 2.6 Carbohydrate production

      GPHOT      = DTGA  * PCEW * 30./44.


* 2.7 Maintenance respiration

      MAINT      = MAINTS * TEFF * MNDVS * EMERG
      MAINTS     = MAINLV*WLVG + MAINST*WST + MAINRT*WRT + MAINSO*WSO
      MNDVS      = WLVG / NOTNUL(WLV)
      TEFF       = Q10**((DAVTMP-TREF)/10.)
PARAMETER Q10    = 2.;    TREF   = 25.
PARAMETER MAINLV = 0.03;  MAINST = 0.015
PARAMETER MAINRT = 0.015; MAINSO = 0.01


* 2.8 Dry matter partitioning

      FSHP       = AFGEN(FSHTB, DVS)
      FSH        = (FSHP * CPEW) / (1. + (CPEW-1.) * FSHP)
      FRT        = 1. - FSH
FUNCTION FSHTB   = 0.00,0.50, 0.10,0.50, 0.20,0.60, 0.35,0.78,...
                   0.40,0.83, 0.50,0.87, 0.60,0.90, 0.70,0.93,...
                   0.80,0.95, 0.90,0.97, 1.00,0.98, 1.10,0.99,...
                   1.20,1.00, 2.50,1.00

      FLV        = AFGEN(FLVTB, DVS) 
      FST        = AFGEN(FSTTB, DVS) 
      FSO        = AFGEN(FSOTB, DVS)
FUNCTION FLVTB   = 0.00,0.65, 0.10,0.65, 0.25,0.70, 0.50,0.50,...
                   0.70,0.15, 0.95,0.00, 2.50,0.00
FUNCTION FSTTB   = 0.00,0.35, 0.10,0.35, 0.25,0.30, 0.50,0.50,...
                   0.70,0.85, 0.95,1.00, 1.05,0.00, 2.50,0.00
FUNCTION FSOTB   = 0.00,0.00, 0.95,0.00, 1.05,1.00, 2.50,1.00

      ERRSH      = ABS(FLV+FST+FSO - 1.)
FINISH ERRSH     > 1.E-6

* 2.9 Growth of plant organs and translocation

      ASRQ       = FSH * (ASRQLV*FLV + ASRQST*FST + ASRQSO*FSO) + ...
                          ASRQRT*FRT
      
      TRANSL     = INSW(DVS-1., 0., WST * DVR * FRTRL)

      GTW        = (GPHOT - MAINT + CONVL*TRANSL*CFST*30./12.) / ASRQ
      GRT        = FRT * GTW 
      GLV        = FLV * FSH * GTW
      GST        = FST * FSH * GTW - TRANSL
      GSO        = FSO * FSH * GTW

* The following values are calculated without 
* the costs ofnitrate reduction:
PARAMETER ASRQRT = 1.444; ASRQLV = 1.463
PARAMETER ASRQST = 1.513; ASRQSO = 1.415
PARAMETER FRTRL  = 0.20 ; CONVL  = 0.947


* 2.10 Leaf and ear development

      TAI        = 0.5 * EAI + LAI
      LAI        = INTGRL(ILAI, RLAI)
      RLAI       = GLAI - DLAI
      CALL GLA(TIME, DOYEM, DTEFF, DVS,           ...
               RGRL, DELT, SLA  , LAI, GLV,      GLAI)
PARAMETER RGRL   = 0.009
PARAMETER SLA    = 0.022

      EAI        = INTGRL(IEAI, REAI)
      CALL SUBEAI(DELT,DVS,EAR,TADRW,RDRDV,EAI,  REAI)
PARAMETER EAR    = 0.63E-3

      DLAI       = LAI * RDR 
      RDR        = MAX(RDRDV, RDRSH)
      RDRDV      = INSW(DVS-1.0, 0., DVR/(MAX(0.1, 2.-DVS))*FRDR)
      RDRSH      = LIMIT(0., 0.03, 0.03 * (LAI-LAICR) / LAICR)
PARAMETER LAICR  = 4.0 ; FRDR = 1.

      DLV        = WLVG * DLAI/NOTNUL(LAI)


* 2.11 Dry matter production

      WRT        = INTGRL(WRTI,  GRT)
      WLVG       = INTGRL(WLVI, RWLVG)
      RWLVG      = GLV - DLV
      WLVD       = INTGRL(WLVDI, DLV)
      WST        = INTGRL(WSTI,  GST)
      WSO        = INTGRL(WSOI,  GSO)

      WLV        = WLVG  + WLVD
      TADRW      = WLV   + WST  + WSO
      TDRW       = TADRW + WRT

      HI         = WSO / NOTNUL(TADRW)


* 2.12 Weather data

WEATHER WTRDIR='C:\SYS\WEATHER\';CNTR='NLD';ISTN=1;IYEAR=1990
* Reading weather data from the weather file:
* RDD    Daily global radiation              J/m2/d
* TMMN   Daily minimum temperature           degree C
* TMMX   Daily maximum temperature           degree C
* VP     Vapour pressure                     kPa
* WN     Wind speed                          m/s
* RAIN   Precipitation                       mm/d
* LAT    Latitude of the site                degree
* DOY    Daynumber of year = TIME            d


      DTR        = RDD
      
      DAVTMP     = 0.5 * (TMMX + TMMN)
      DDTMP      = TMMX - 0.25 * (TMMX-TMMN)
      DTEFF      = MAX(0., DAVTMP-TBASE)
PARAMETER TBASE  = 0.

      AVP        = VP     
      WDS        = WN
      RRAIN      = RAIN
      
      TRAIN      = INTGRL(ZERO, RRAIN)

* 2.13 Penman-Monteith combination equation

      PENMAN     = EVAPR + EVAPD
      EVAPR      = (1./LHVAP) * (SLOPE/(SLOPE+PSYCH)) * NRAD
      SLOPE      = 4158.6 * SVP / (DAVTMP + 239.)**2
      SVP        = 0.611 * EXP(17.4 * DAVTMP / (DAVTMP + 239.))
PARAMETER LHVAP  = 2.4E6
PARAMETER PSYCH  = 0.067

      NRAD       = (1.-ALB)* DTR - RLWN
      ALB        = ALBS*EXP(-0.5*LAI) + 0.25*(1.-EXP(-0.5*LAI))
      ALBS       = 0.25 * (1.-0.5*WCL1/WCST1)
      
      RLWN       = BBRAD * FVAP * FCLEAR * 86400.
      BBRAD      = BOLTZM * (DAVTMP+273.)**4
PARAMETER BOLTZM = 5.668E-8
      FVAP       = 0.56-0.079*SQRT(AVP*10.)
      FCLEAR     = 0.1+0.9*CLEAR
      CLEAR      = LIMIT(0., 1., ((DTR/DS0)-A)/B)
PARAMETER A      = 0.25; B=0.45

      WDF        = 2.63 * (1.0 + 0.54 * WDS) * PSYCH
      DRYP       = (SVP-AVP) * WDF
      EVAPD      = DRYP/(SLOPE+PSYCH)

      TPENM      = INTGRL(ZERO, PENMAN)
      TEVAPR     = INTGRL(ZERO, EVAPR)
      TEVAPD     = INTGRL(ZERO, EVAPD)


* 2.14 The soil water balance

      AINTC      = MIN(RRAIN, INTC*LAI)
PARAMETER INTC   = 0.25

      RNOFF      = MAX(0., 0.15*(RRAIN-AINTC-10.),                 ...
                           RRAIN-AINTC-(WCST1*TKL1-WL1)/(2.*DELT))
      
      WLFL1      = RRAIN - AINTC - RNOFF
      WLFL2      = MAX(0., MIN(WL1-WCFC1*TKL1,                    ...
                           WCST2*TKL2-WL2)/(2.*DELT))
      WLFL3      = MAX(0., MIN(WL2-WCFC2*TKL2,                    ...
                           WCST3*TKL3-WL3)/(2.*DELT))
      WLFL4      = MAX(0., MIN(WL3-WCFC3*TKL3,                    ...
                           WCST4*TKL4-WL4)/(2.*DELT))
      WLFL5      = MAX(0., MIN((WL4-WCFC4*TKL4)/(2.*DELT),MDRATE))

PARAMETER WCFC1  = 0.23; WCFC2 = 0.23; WCFC3 = 0.23; WCFC4 = 0.23
PARAMETER WCST1  = 0.40; WCST2 = 0.40; WCST3 = 0.40; WCST4 = 0.40

      DRAIN      = WLFL5
      
      RWL1       = WLFL1-WLFL2-EVSW1-TRWL1
      RWL2       = WLFL2-WLFL3-EVSW2-TRWL2
      RWL3       = WLFL3-WLFL4-EVSW3-TRWL3
      RWL4       = WLFL4-WLFL5-EVSW4-TRWL4

      WL1        = INTGRL(WL1I,RWL1)
      WL2        = INTGRL(WL2I,RWL2)
      WL3        = INTGRL(WL3I,RWL3)
      WL4        = INTGRL(WL4I,RWL4)
      WCL1       = WL1/TKL1
      WCL2       = WL2/TKL2
      WCL3       = WL3/TKL3
      WCL4       = WL4/TKL4
      RWCL1      = (WCL1-WCWP1)/(WCFC1-WCWP1)
      RWCL2      = (WCL2-WCWP2)/(WCFC2-WCWP2)
      RWCL3      = (WCL3-WCWP3)/(WCFC3-WCWP3)
      RWCL4      = (WCL4-WCWP4)/(WCFC4-WCWP4)
      
      TDRAIN     = INTGRL(ZERO, DRAIN)
      TAINTC     = INTGRL(ZERO, AINTC)
      TRNOFF     = INTGRL(ZERO, RNOFF)

      WCUM       = WL1+WL2+WL3+WL4
      CHECK      = TRAIN+WCUMI-TAINTC-TRNOFF-TDRAIN-WCUM-...
                   TATRAN-TAEVAP


* 2.15 Rooted depth

      ZRT        = INTGRL(ZRTI, EZRT)
      EZRT       = EZRTC * WSERT * AMTMP
* Temperature effect included, same as for AMAX (AMTMP)

     CALL SUBGRT (ZRT,ZRTM,DVS,TKL1,TKL2,TKL3,TKL4,WCL1,WCL2, ...
                  WCL3,WCL4,WCWP1,WCWP2,WCWP3,WCWP4,  WSERT)
PARAMETER EZRTC  = 12.



* 2.16 Transpiration

      PTRANS     = (1. - EXP(-0.5*LAI)) * EVAPR + EVAPD * ...
                    MIN(2.0, LAI) - 0.5 * AINTC

FUNCTION EDPTFT  = -.50,0., -.05,0., 0.,.15, .15,.6, ...
                     .3,.8, .5,1., 2.,1.
      ERLB       = ZRT1*AFGEN(EDPTFT, RWCL1) +...
                   ZRT2*AFGEN(EDPTFT, RWCL2) +...
                   ZRT3*AFGEN(EDPTFT, RWCL3) +...
                   ZRT4*AFGEN(EDPTFT, RWCL4)
      TRRM       = PTRANS/NOTNUL(ERLB)
      TRWL1      = TRRM*WSE1*ZRT1*AFGEN(EDPTFT, RWCL1)
      TRWL2      = TRRM*WSE2*ZRT2*AFGEN(EDPTFT, RWCL2)
      TRWL3      = TRRM*WSE3*ZRT3*AFGEN(EDPTFT, RWCL3)
      TRWL4      = TRRM*WSE4*ZRT4*AFGEN(EDPTFT, RWCL4)
      ZRT1       = LIMIT(0.,TKL1,ZRT)
      ZRT2       = LIMIT(0.,TKL2,ZRT-TKL1)
      ZRT3       = LIMIT(0.,TKL3,ZRT-TKL1-TKL2)
      ZRT4       = LIMIT(0.,TKL4,ZRT-TKL1-TKL2-TKL3)

      ATRANS     = TRWL1+TRWL2+TRWL3+TRWL4

      TPTRAN     = INTGRL(ZERO, PTRANS)
      TATRAN     = INTGRL(ZERO, ATRANS)


* 2.17 Evaporation

      PEVAP      = EXP(-0.5*LAI) * (EVAPR + EVAPD)
      
      DSLR       = INTGRL(IDSLR, RDSLR)
      RDSLR      = INSW(RRAIN-0.5, 1., -(DSLR-1.)/DELT)
      
      AEVAP      = INSW(RRAIN-0.5, EVSD, EVSH)
      EVSH       = MIN(PEVAP, (WL1-WCAD1*TKL1)/DELT+WLFL1)
PARAMETER WCAD1  = 0.025; WCAD2 = 0.025; WCAD3 = 0.025; WCAD4 = 0.025

      EVSD       = MIN(PEVAP, 0.6*PEVAP*(SQRT(DSLR+1.)-...
                   SQRT(DSLR))+WLFL1)

      FEVL1      = MAX(WL1-WCAD1*TKL1, 0.1)*EXP(-EES*(0.5*TKL1))
      FEVL2      = MAX(WL2-WCAD2*TKL2, 0.1)*EXP(-EES*(TKL1+...
                   (0.5*TKL2)))
      FEVL3      = MAX(WL3-WCAD3*TKL3, 0.1)*EXP(-EES*(TKL1+TKL2+...
                   (0.5*TKL3)))
      FEVL4      = MAX(WL4-WCAD4*TKL4, 0.1)*EXP(-EES*(TKL1+TKL2+TKL3+...
                   (0.5*TKL4)))
PARAMETER EES    = 0.002

      FEVLT      = FEVL1+FEVL2+FEVL3+FEVL4
      EVSW1      = AEVAP*(FEVL1/FEVLT)
      EVSW2      = AEVAP*(FEVL2/FEVLT)
      EVSW3      = AEVAP*(FEVL3/FEVLT)
      EVSW4      = AEVAP*(FEVL4/FEVLT)

      TPEVAP     = INTGRL(ZERO, PEVAP)
      TAEVAP     = INTGRL(ZERO, AEVAP)


* 2.18 Effects of water stress

      P          = TRANSC/(TRANSC+PTRANS)
PARAMETER TRANSC = 9.

      CALL SUBFR (WCL1,WCFC1,P,WCWP1,WCWET1,WCST1, WSE1)
      CALL SUBFR (WCL2,WCFC2,P,WCWP2,WCWET2,WCST2, WSE2)
      CALL SUBFR (WCL3,WCFC3,P,WCWP3,WCWET3,WCST3, WSE3)
      CALL SUBFR (WCL4,WCFC4,P,WCWP4,WCWET4,WCST4, WSE4)

PARAMETER WCWET1 = 0.35; WCWET2 = 0.35; WCWET3 = 0.35; WCWET4 = 0.35
PARAMETER WCWP1  = 0.075;WCWP2  = 0.075;WCWP3  = 0.075; WCWP4 = 0.075

      PCEW       = ATRANS/NOTNUL(PTRANS)
      CPEW       = MIN(1., 0.5+ATRANS/NOTNUL(PTRANS))


* 2.19 Water use efficiency

      TDTGA      = INTGRL(ZERO, DTGA)
      TAR        = TATRAN*1.E3/NOTNUL(TDTGA)
*      TRC        = TATRAN*1.E4/(NOTNUL(TDRW)
      CROPF      = (PTRANS+PEVAP)/PENMAN


*  2.20 Carbon balance check
      
      CHKIN      = WLV * CFLV + WST * CFST + ...
                   WRT * CFRT + WSO * CFSO
      CHKFL      = TNASS * (12./44.) 
      TNASS      = INTGRL(TNASSI, RTNASS)
      RTNASS     = ((GPHOT - MAINT)*44./30.) -  ...
                                 (GRT*CO2RT + GLV*CO2LV +  ...
                         (GST+TRANSL)*CO2ST + GSO*CO2SO +  ...
                         (1.-CONVL)* TRANSL*CFST*44./12.)
      CO2RT      = 44./12. * (ASRQRT*12./30. - CFRT)
      CO2LV      = 44./12. * (ASRQLV*12./30. - CFLV)
      CO2ST      = 44./12. * (ASRQST*12./30. - CFST)
      CO2SO      = 44./12. * (ASRQSO*12./30. - CFSO)

      CHKDIF     = (CHKIN-CHKFL)/NOTNUL(CHKIN)

PARAM CFLV=0.459; CFST=0.494; CFRT=0.467; CFSO=0.471


* 2.21 Run control

*      DAY        = 1. + AMOD(TIME-1., 365.)

FINISH DVS       > 2.
TIMER STTIME     = 80.; FINTIM = 300.; DELT = 1.; PRDEL = 5.
TRANSLATION_GENERAL DRIVER='EUDRIV'

PRINT DOY,DVS,DAYL,TDRW,TADRW,WLVG,WLVD,WLV,WST,...
      WSO,WRT,LAI,EAI,HI,WL1,WL2,WL3,WL4,...
      TPENM,TEVAPR,TEVAPD,TRAIN,TAINTC,TRNOFF,TDRAIN,...
      TPTRAN,TATRAN,TPEVAP,TAEVAP,CHECK,TAR,EVAPR,...
      EVAPD,CROPF,WSE1,WSE2,WSE3,WSE4,PCEW,CPEW,ZRT,ERLB,WCUM,...
      RWCL1,RWCL2,RWCL3,RWCL4,CHKDIF, P

END
STOP


* 2.22 SUBROUTINES

* ---------------------------------------------------------------------*
*  SUBROUTINE GLA                                                      *
*  Purpose: This subroutine computes daily increase of leaf area index *
*           (m2 leaf/ m2 ground/ d)                                    *
*  Version: September 1997 (SUCROS97 V1.0)                             *
*                                                                      *
*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)     *
*  name   type meaning                                    units  class *
*  ----   ---- -------                                    -----  ----- *
*  TIME    R4  Time in simulation                            d      T  *
*  DOYEM   R4  Day number of crop emergence                  d      T  *
*  DTEFF   R4  Daily effective temperature                  oC      I  *
*  DVS     R4  Development stage of the crop                 -      I  *
*  RGRL    R4  Relative leaf growth rate                 oC-1 d-1   I  *
*  DELT    R4  Time step of integration                     d       T  *
*  SLA     R4  Specific leaf area                          m2 g-1   I  *
*  LAI     R4  Leaf area index                             m2 m-2   I  *
*  GLV     R4  Growth rate of the leaves                  g m-2 d-1 I  *
*  GLAI    R4  Growth rate of leaf area index            m2 m-2 d-1 O  *
* ---------------------------------------------------------------------*

      SUBROUTINE GLA (TIME,DOYEM,DTEFF,DVS,RGRL,DELT,SLA,LAI,GLV,
     $                GLAI)
      IMPLICIT REAL (A-Z)

*-----growth during maturation stage
      GLAI  = SLA * GLV

*-----growth during juvenile stage
      IF ((DVS.LT.0.3).AND.(LAI.LT.0.75)) THEN
      GLAI  = (LAI * (EXP(RGRL*DTEFF*DELT)-1.))/DELT
      ENDIF

*-----growth before seedling emergence
      IF (TIME.LE.DOYEM) GLAI = 0.

      RETURN
      END


*----------------------------------------------------------------------*
*  SUBROUTINE SUBEAI                                                   *
*  Purpose: This subroutine calculates ear area index                  *
*  Version: September 1997 (SUCROS97 V1.0)                             *
*                                                                      *
*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)     *
*  name   type meaning                                    units  class *
*  ----   ---- -------                                    -----  ----- *
*  DELT    R4  Time step of integration                     d       T  *
*  DVS     R4  Development stage of the crop                -       I  *
*  EAR     R4  Ear area/weight ratio                       g m-2    I  *
*  TADRW   R4  Total above-ground dry weight               g m-2    I  *
*  RDR     R4  Relative death rate                          d-1     I  *
*  EAI     R4  Ear area index                             m2 m-2    I  *
*  REAI    R4  Growth rate ear area index               m2 m-2 d-1  O  *
*----------------------------------------------------------------------*

      SUBROUTINE SUBEAI(DELT,DVS,EAR,TADRW,RDRDV,EAI, REAI)
      IMPLICIT REAL(A-Z)

      IF (DVS.LT.0.8) REAI = 0.
      IF (DVS.GE.0.8 .AND. EAI.EQ.0.) THEN
                      REAI = (EAR * TADRW)/DELT
      ELSE
                      REAI = 0.
      ENDIF
      IF (DVS.GE.1.3) REAI = -RDRDV * EAI
      RETURN
      END


*----------------------------------------------------------------------*
*  SUBROUTINE ASTRO                                                    *
*  Purpose: This subroutine calculates astronomic daylength,           *
*           diurnal radiation characteristics such as the daily        *
*           integral of sine of solar elevation and solar constant.    *
*  Version: September 1997 (SUCROS97 V1.0)                             *
*                                                                      *
*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)     *
*  name   type meaning                                    units  class *
*  ----   ---- -------                                    -----  ----- *
*  DOY     R4  Daynumber (Jan 1st = 1)                       -      T  *
*  LAT     R4  Latitude of the site                       degrees   I  *
*  SC      R4  Solar constant                             J m-2 s-1 O  *
*  DS0     R4  Daily extraterrestrial radiation           J m-2 d-1 O  *
*  SINLD   R4  Seasonal offset of sine of solar height       -      O  *
*  COSLD   R4  Amplitude of sine of solar height             -      O  *
*  DAYL    R4  Astronomic daylength (base = 0 degrees)       h      O  *
*  DSINB   R4  Daily total of sine of solar height           s      O  *
*  DSINBE  R4  Daily total of effective solar height         s      O  *
*                                                                      *
*  FATAL ERROR CHECKS (execution terminated, message)                  *
*  condition: LAT > 67, LAT < -67                                      *
*----------------------------------------------------------------------*

      SUBROUTINE ASTRO (DOY, LAT,
     &                  SC , DS0, SINLD, COSLD, DAYL, DSINB, DSINBE)
      IMPLICIT REAL (A-Z)

*-----PI and conversion factor from degrees to radians
      PI    = 3.141592654
      RAD   = PI/180.

*-----check on input range of parameters
      IF (LAT.GT.67.)  STOP 'ERROR IN ASTRO: LAT> 67'
      IF (LAT.LT.-67.) STOP 'ERROR IN ASTRO: LAT>-67'

*-----declination of the sun as function of daynumber (DOY)
      DEC   = -ASIN (SIN (23.45*RAD)*COS (2.*PI*(DOY+10.)/365.))

*-----SINLD, COSLD and AOB are intermediate variables

      SINLD = SIN (RAD*LAT)*SIN (DEC)
      COSLD = COS (RAD*LAT)*COS (DEC)
      AOB   = SINLD/COSLD

*-----daylength (DAYL) 
      DAYL   = 12.0*(1.+2.*ASIN (AOB)/PI)

      DSINB  = 3600.*(DAYL*SINLD+24.*COSLD*SQRT (1.-AOB*AOB)/PI)
      DSINBE = 3600.*(DAYL*(SINLD+0.4*(SINLD*SINLD+COSLD*COSLD*0.5))+
     &         12.0*COSLD*(2.0+3.0*0.4*SINLD)*SQRT (1.-AOB*AOB)/PI)

*-----solar constant (SC) and daily extraterrestrial radiation (DS0) 
      SC  = 1370.*(1.+0.033*COS (2.*PI*DOY/365.))
      DS0 = SC*DSINB

      RETURN
      END


*----------------------------------------------------------------------*
*  SUBROUTINE TOTASS                                                   *
*  Purpose: This subroutine calculates daily total gross               *
*           assimilation (DTGA) by performing a Gaussian integration   *
*           over time. At three different times of the day,            *
*           radiation is computed and used to determine assimilation   *
*           whereafter integration takes place.                        *
*  Version: September 1997 (SUCROS97 V1.0)                             *
*                                                                      *
*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)     *
*  name   type meaning                                    units  class *
*  ----   ---- -------                                    -----  ----- *
*  DOY     R4  Day number (January 1 = 1)                     -     T  *
*  LAT     R4  Latitude of the site                       degrees   I  *
*  DTR     R4  Daily total of global radiation            J/m2/d    I  *
*  SCP     R4  Scattering coefficient of leaves for visible            *
*              radiation (PAR)                                -     I  *
*  AMAX    R4  Assimilation rate at light saturation       g CO2/   I  *
*                                                        m2 leaf/s     *
*  EFF     R4  Initial light conversion factor            g CO2/J   I  *
*  KDF     R4  Extinction coefficient diffuse flux leaves    -      I  *
*  LAI     R4  Leaf area index as used for photosynthesis  m2/m2    I  *
*              Note: This can involve stem, flower or                  *
*                    ear area index!!                                  *
*  DAYL    R4  Astronomic daylength (base = 0 degrees)       h      O  *
*  DTGA    R4  Daily total gross assimilation            g CO2/m2/d O  *
*  DS0     R4  Daily extraterrestrial radiation            J/m2/s   O  *
*                                                                      *
*                                                                      *
*  SUBROUTINES called : ASTRO, ASSIM                                   *
*----------------------------------------------------------------------*

      SUBROUTINE TOTASS (DOY,  LAT , DTR, SCP, AMAX, EFF, KDF, LAI,
     &                   DAYL, DTGA, DS0)
      IMPLICIT REAL(A-Z)
      REAL XGAUSS(3), WGAUSS(3)
      INTEGER I1, IGAUSS

      DATA IGAUSS /3/
      DATA XGAUSS /0.112702, 0.500000, 0.887298/
      DATA WGAUSS /0.277778, 0.444444, 0.277778/

      PI   = 3.141592654

      CALL ASTRO(DOY,LAT,SC,DS0,SINLD,COSLD,DAYL,DSINB,DSINBE)

*-----assimilation set to zero and three different times of the day (HOUR)
      DTGAS = 0.
      
      DO 10 I1=1,IGAUSS

*--------at the specified HOUR, radiation is computed and used to compute
*        assimilation
         HOUR = 12.0+DAYL*0.5*XGAUSS(I1)
 
*--------sine of solar elevation
         SINB  = MAX (0., SINLD+COSLD*COS (2.*PI*(HOUR+12.)/24.))

*--------diffuse light fraction (FRDF) from atmospheric 
*        transmission (ATMTR)
         PAR   = 0.5*DTR*SINB*(1.+0.4*SINB)/DSINBE
         ATMTR = PAR/(0.5*SC*SINB)

         IF (ATMTR.LE.0.22) THEN
            FRDF = 1.
         ELSE IF (ATMTR.GT.0.22 .AND. ATMTR.LE.0.35) THEN
            FRDF = 1.-6.4*(ATMTR-0.22)**2
         ELSE
            FRDF = 1.47-1.66*ATMTR
         END IF

         FRDF  = MAX (FRDF, 0.15+0.85*(1.-EXP (-0.1/SINB)))

*--------diffuse PAR (PARDF) and direct PAR (PARDR)
         PARDF = PAR * FRDF
         PARDR = PAR - PARDF

         CALL ASSIM (SCP,AMAX,EFF,KDF,LAI,SINB,PARDR,PARDF,FGROS)

*--------integration of assimilation rate to a daily total (DTGA) 
         DTGAS = DTGAS+FGROS*WGAUSS(I1)

10    CONTINUE

      DTGA = DTGAS * DAYL * 3600.

      RETURN
      END


*----------------------------------------------------------------------*
*  SUBROUTINE ASSIM                                                    *
*  Purpose: This subroutine performs a Gaussian integration over       *
*           depth of canopy by selecting five different LAI's and      *
*           computing assimilation at these LAI levels. The            *
*           integrated variable is FGROS.                              *
*  Version: September 1997 (SUCROS97 V1.0)                             *
*                                                                      *
*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)     *
*  name   type meaning                                    units  class *
*  ----   ---- -------                                    -----  ----- *
*  SCP     R4  Scattering coefficient of leaves for visible            *
*              radiation (PAR)                              -       I  *
*  AMAX    R4  Assimilation rate at light saturation       g CO2/   I  *
*                                                        m2 leaf/s     *
*  EFF     R4  Initial light conversion factor            g CO2/J   I  *
*  KDF     R4  Extinction coefficient diffuse flux leaves    -      I  *
*  LAI     R4  Leaf area index as used for photosynthesis  m2/m2    I  *
*              Note: This can involve stem, flower or                  *
*                    ear area index!!                                  *
*  SINB    R4  Sine of solar height                          -      I  *
*  PARDR   R4  Instantaneous flux of direct radiation (PAR) W/m2    I  *
*  PARDF   R4  Instantaneous flux of diffuse radiation(PAR) W/m2    I  *
*  FGROS   R4  Instantaneous assimilation rate of          g CO2/   O  *
*              whole canopy                              m2 soil/s     *
*                                                                      *
*----------------------------------------------------------------------*

      SUBROUTINE ASSIM (SCP, AMAX, EFF, KDF, LAI, SINB, PARDR, PARDF,
     &                  FGROS)
      IMPLICIT REAL(A-Z)
      REAL XGAUSS(5), WGAUSS(5)
      INTEGER I1, I2, IGAUSS

*-----Gauss weights for five point Gauss
      DATA IGAUSS /5/
      DATA XGAUSS /0.0469101,0.2307534,0.5      ,0.7692465,0.9530899/
      DATA WGAUSS /0.1184635,0.2393144,0.2844444,0.2393144,0.1184635/

*-----reflection of horizontal and spherical leaf angle distribution
      SQV  = SQRT(1.-SCP)
      REFH = (1.-SQV)/(1.+SQV)
      REFS = REFH*2./(1.+2.*SINB)

*-----extinction coefficient for direct radiation and total direct flux
      CLUSTF = KDF / (0.8*SQV) 
      KBL    = (0.5/SINB) * CLUSTF
      KDRT   = KBL * SQV

*-----selection of depth of canopy, canopy assimilation is set to zero
      FGROS = 0.
      
      DO 10 I1=1,IGAUSS
         LAIC = LAI * XGAUSS(I1)

*--------absorbed fluxes per unit leaf area: diffuse flux, total direct
*        flux, direct component of direct flux.
         VISDF = (1.-REFH)*PARDF*KDF  *EXP (-KDF  *LAIC)
         VIST  = (1.-REFS)*PARDR*KDRT *EXP (-KDRT *LAIC)
         VISD  = (1.-SCP) *PARDR*KBL  *EXP (-KBL  *LAIC)

*--------absorbed flux (J/M2 leaf/s) for shaded leaves and assimilation of
*        shaded leaves
         VISSHD = VISDF + VIST - VISD
         IF (AMAX.GT.0.) THEN
            FGRSH  = AMAX * (1.-EXP(-VISSHD*EFF/AMAX))
         ELSE
            FGRSH = 0.
         END IF
 
*--------direct flux absorbed by leaves perpendicular on direct beam and
*        assimilation of sunlit leaf area

         VISPP  = (1.-SCP) * PARDR / SINB
         FGRSUN = 0.
         DO 20 I2=1,IGAUSS
            VISSUN = VISSHD + VISPP * XGAUSS(I2)
            IF (AMAX.GT.0.) THEN
               FGRS = AMAX * (1.-EXP(-VISSUN*EFF/AMAX))
            ELSE
               FGRS = 0.
            END IF
            FGRSUN = FGRSUN + FGRS * WGAUSS(I2)
20       CONTINUE

*--------fraction sunlit leaf area (FSLLA) and local assimilation 
*        rate (FGL)
         FSLLA = CLUSTF * EXP(-KBL*LAIC)
         FGL   = FSLLA  * FGRSUN + (1.-FSLLA) * FGRSH

*--------integration of local assimilation rate to canopy 
*        assimilation (FGROS)
         FGROS = FGROS + FGL * WGAUSS(I1)

10    CONTINUE
      FGROS = FGROS * LAI

      RETURN
      END

*----------------------------------------------------------------------*
*  Subroutine SUBGRT                                                   *
*  Purpose: To decide whether root extension growth continues or       *
*           ceases (value either 0 or 1)                               *
*  Version: September 1997 (SUCROS97 V1.0)                             *
*                                                                      *
*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)     *
*  name   type meaning                                    units  class *
*  ----   ---- -------                                    -----  ----- *
*  ZRT     R4  Rooted depth                                 mm     I   *
*  ZRTM    R4  Maximum value for rooted depth               mm     I   *
*  DVS     R4  Development stage of the crop                 -     I   *
*  TKL1-4  R4  Thickness of the soil layers                 mm     I   *
*  WCL1-4  R4  Volumetric water content in soil layers    cm3/cm3  I   *
*  WCWP1-4 R4  Volumetric water content at wilting point  cm3/cm3  I   *
*  WSERT   R4  Variable to calculate root extension          -     O   *
*----------------------------------------------------------------------*

      SUBROUTINE SUBGRT(ZRT,ZRTM,DVS,TKL1,TKL2,TKL3,TKL4,WCL1,WCL2, 
     $                  WCL3,WCL4,WCWP1,WCWP2,WCWP3,WCWP4, WSERT)
      IMPLICIT REAL(A-Z)
      SAVE

      WSERT = 1.
      IF (ZRT.LT.TKL1 .AND. WCL1.LT.WCWP1) WSERT = 0.
      IF (ZRT.GT.TKL1 .AND. ZRT.LT.(TKL1+TKL2) .AND. 
     $    WCL2.LT.WCWP2) WSERT = 0.
      IF (ZRT.GT.(TKL1+TKL2) .AND. ZRT.LT.(TKL1+TKL2+TKL3) .AND. 
     $    WCL3.LT.WCWP3)  WSERT = 0.
      IF (ZRT.GT.(TKL1+TKL2+TKL3) .AND. ZRT.LT.(TKL1+TKL2+TKL3+TKL4) 
     $    .AND. WCL4.LT.WCWP4)  WSERT = 0.

      IF (DVS.GE.1.)    WSERT=0.
      IF (ZRT.GT.ZRTM)  WSERT=0.
      RETURN
      END


*----------------------------------------------------------------------*
*  Subroutine SUBFR                                                    *
*  Purpose: To compute factors accounting for water stress effect on   *
*           water uptake                                               *
*           Revision 24 Nov 96 by Paul Kiepe                           *
*           To avoid division by zero, if P=0.                         *
*  Version: September 1997 (SUCROS97 V1.0)                             *
*                                                                      *
*  FORMAL PARAMETERS:  (I=input,O=output,C=control,IN=init,T=time)     *
*  name   type meaning                                    units  class *
*  ----   ---- -------                                    -----  ----- *
*  WCL    R4   Volumetric water content in soil layers    cm3/cm3   I  *
*  WCFC   R4   Volumetric water content at field capacity cm3/cm3   I  *
*  P      R4   Soil water depletion factor                   -      I  *
*  WCWP   R4   Volumetric water content at wilting point  cm3/cm3   I  *
*  WCWET  R4   Volumetric water content where                          *
*              water logging begins                       cm3/cm3   I  *
*  WCST   R4   Volumetric water contentat saturation      cm3/cm3   I  *
*  WSE    R4   Factor accounting for effect of uptake                  *
*              availability of soil water                    -      O  *
*----------------------------------------------------------------------*

      SUBROUTINE SUBFR(WCL,WCFC,P,WCWP,WCWET,WCST, WSE)
      IMPLICIT REAL (A-Z)
      SAVE

      WCCR = WCWP + (1.-P) * (WCFC - WCWP)

      IF (WCL.GT.WCWET) THEN 
          FR = (WCST-WCL)/(WCST-WCWET)
      ELSE IF (WCL.LE.WCWET .AND. WCL.GT.WCCR) THEN
          FR = 1.
      ELSE IF (WCL.LE.WCCR .AND. WCL.GT.WCWP) THEN
          FR = (WCL-WCWP)/(WCCR-WCWP)
      ELSE
          FR = 0.
      ENDIF

      WSE = MIN(1.,MAX(0., FR))

      RETURN
      END
