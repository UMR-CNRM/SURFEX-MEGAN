!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!    ###############################
     SUBROUTINE COUPLING_MEGAN_n(MGN, CHI, GK, PEK, KK,     &
                     KYEAR, KMONTH, KDAY, PTIME, OTR_ML, &
                     KSLTYP, PPFT, PEF,                  &
                     PTEMP, PIACAN, PLEAFT, PRN_SUNLIT, PRN_SHADE, &
                     PWIND, PPRES, PQV, PSFTS, T24, P_PPFD, PPFD24, ICOUNT, &
                     PGAMMA_SM, PGAMMA_AGE, PACTIVITY_FACTOR, PCO2, PGAMMA_CO2) 
!    ###############################
!!
!!***  *BVOCEM*
!! 
!!    PURPOSE
!!    -------
!!    Calculate the biogenic emission fluxes upon the MEGAN code
!!    http://lar.wsu.edu/megan/
!!
!!    METHOD
!!    ------
!!
!!
!!    AUTHOR
!!    ------
!!    P. Tulet (LACy)
!!    
!!    MODIFICATIONS
!!    -------------
!!    Original: 25/10/2014
!!    Modified: 06/07/2017, J. Pianezze, adaptation for SurfEx v8.0
!!    Modified: 06/07/2018, P. Tulet, correction for T leaf
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!
USE MODD_MEGAN_n, ONLY : MEGAN_t
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_ISBA_n, ONLY: ISBA_PE_t, ISBA_K_t
USE MODD_SFX_GRID_n, ONLY: GRID_t
!
USE MODD_CSTS, ONLY : XAVOGADRO
!
#ifdef MNH_MEGAN
USE MODD_MEGAN
USE MODI_JULIAN
USE MODI_EMPROC
USE MODI_MGN2MECH
#endif
!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
!
IMPLICIT NONE
!
TYPE(MEGAN_t), INTENT(INOUT) :: MGN
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(GRID_t), INTENT(INOUT) :: GK
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
TYPE(ISBA_K_t), INTENT(INOUT) :: KK
!
!*       0.1  declaration of arguments
!
INTEGER,             INTENT(IN)     :: KYEAR     ! I current year (UTC)
INTEGER,             INTENT(IN)     :: KMONTH    ! I current month (UTC)
INTEGER,             INTENT(IN)     :: KDAY      ! I current day (UTC)
REAL,                INTENT(IN)     :: PTIME     ! I current time since midnight (UTC, s)
INTEGER                             :: ICOUNT    ! day counter 
LOGICAL,             INTENT(IN)     :: OTR_ML    ! new radiation for leaves temperatures
!
REAL, DIMENSION(:),  INTENT(IN)     :: PTEMP     ! I Air temperature (K)
REAL, DIMENSION(:),INTENT(IN)     :: PIACAN    ! I PAR (umol/m2.s)
REAL, DIMENSION(:),  INTENT(IN)     :: PLEAFT    ! I Leaf temperature (K)
REAL, DIMENSION(:),  INTENT(IN)     :: PRN_SUNLIT! I Leaf RN
REAL, DIMENSION(:),  INTENT(IN)     :: PRN_SHADE ! I Leaf RN
REAL, DIMENSION(:),  INTENT(IN)     :: PWIND
REAL, DIMENSION(:),  INTENT(IN)     :: PPRES     ! I Atmospheric pressure (Pa)
REAL, DIMENSION(:),  INTENT(IN)     :: PQV       ! I Air humidity (kg/kg)
REAL, DIMENSION(:),  INTENT(IN)     :: PCO2       ! CO2 concentration (ppm)
REAL, DIMENSION(:,:),INTENT(IN)     :: PPFT, PEF
INTEGER, DIMENSION(:),  INTENT(IN)     :: KSLTYP
REAL, DIMENSION(:), INTENT(IN) :: T24            ! average temperature over the past 24h (K)
REAL, DIMENSION(:), INTENT(IN) :: PPFD24         ! average PAR over the past 24h         (umol/m2.s)
REAL, DIMENSION(:,:), INTENT(INOUT)  :: PSFTS    ! O Scalar flux in molecules/m2/s
REAL, DIMENSION(:), INTENT(OUT)  :: P_PPFD, PGAMMA_SM, PGAMMA_AGE, PACTIVITY_FACTOR, PGAMMA_CO2
#ifdef MNH_MEGAN
!*   0.1 Declaration of local variables
!
INTEGER, PARAMETER :: NROWS = 1
INTEGER            ::   ITIME              ! Time of the day HHMMSS
INTEGER            ::   IDATE              ! Date YYYYDDD
INTEGER            ::   IDAY               !  julian day
REAL               ::   ZHOUR, ZMIN, ZSEC  ! conversion ptime to itime format
REAL, DIMENSION(SIZE(PTEMP)) :: ZLAIC, ZWWILT, ZLAIP    ! Current monthly LAI
REAL, DIMENSION(SIZE(PTEMP)) :: ZPFD     ! Calculated PAR (umol/m2.s)
REAL, DIMENSION(SIZE(PTEMP)) :: ZLSUT    ! Leaf on sun temperature (K)
REAL, DIMENSION(SIZE(PTEMP)) :: ZLSHT    ! Leaf on shade temperature (K)
REAL, DIMENSION(SIZE(PTEMP)) :: ZRN
REAL, DIMENSION(SIZE(PTEMP)) :: ZCFNO    ! NO correction factor
REAL, DIMENSION(SIZE(PTEMP)) :: ZCFNOG   ! NO correction factor for grass
REAL, DIMENSION(N_MGN_SPC,SIZE(PTEMP)) :: ZCFSPEC  ! Output emission buffer
REAL, DIMENSION(MGN%NVARS3D,SIZE(PTEMP)) :: ZFLUX   ! Output emission megan flux
REAL, DIMENSION(SIZE(PTEMP)) :: ZD_TEMP  ! Daily temperature (K)
REAL, DIMENSION(SIZE(PTEMP)) :: ZD_PPFD  ! Daily PAR (umol/m2.s)


REAL, DIMENSION(SIZE(PTEMP)) :: ZALPHAP, ZALPHAT, ZPPFD, ZGAMMA_SM, ZGAMMA_AGE, ZACTIVITY_FACTOR, ZGAMMA_CO2 !!


!
REAL :: ZDI      ! Drought Index (0 normal, -2  moderate drought, -3 severe drought, -4 extreme drought)
REAL :: ZREC_ADJ ! Rain adjustment factor
REAL :: MM_ISO   ! Masse Molaire de l'isoprene en g/mol
REAL :: MM_MYRC  ! MASSE Molaire du myrcene en g/mol
REAL :: NB_AVO   ! nombre d'AVOGADRO  en mol ⁻¹
!
INTEGER,DIMENSION(SIZE(PTEMP)) ::   ISLTYP  !Soil category (function of silt, clay and sand))
INTEGER :: JSV, JSM
!
!
! Input parameters
MM_ISO = 68.117 
MM_MYRC = 136.23 
ZHOUR = FLOAT(INT(PTIME/3600.))
ZMIN  = FLOAT(INT((PTIME - ZHOUR*3600) / 60.))
ZSEC  = FLOAT(INT(PTIME - ZHOUR*3600. - ZMIN * 60.))
ITIME = INT(ZHOUR)*10000 + INT(ZMIN)*100 + ZSEC
IDAY  = JULIAN(KYEAR, KMONTH, KDAY)
IDATE = KYEAR*1000 + IDAY
!
ZLAIC(:) = MIN(MAX(0.001,PEK%XLAI(:)),8.)
ZLAIP(:) = MIN(MAX(0.001,PEK%XLAIp(:)),8.)
ZWWILT(:) = KK%XWWILT(:,5)
!
ZDI      = MGN%XDROUGHT
ZREC_ADJ = MGN%XMODPREC  
!
!PPFD24 AND T24
IF (ICOUNT == 0) THEN
    ZD_TEMP  = SPREAD(MGN%XDAILYTEMP, DIM = 1, NCOPIES = SIZE(PTEMP))
    ZD_PPFD  = SPREAD(MGN%XDAILYPAR, DIM = 1, NCOPIES = SIZE(PTEMP))
ELSE
    ZD_TEMP = T24
    ZD_PPFD = PPFD24
ENDIF
!
ZCFNO   = 0.
ZCFNOG  = 0.
ZCFSPEC = 0.
ZALPHAP = 0.
ZALPHAT = 0.
ZPPFD = 0.
!
ZPFD(:) = 0.
! Compute PPFD from incoming solar radiation
ZPFD(:) = PIACAN(:)*0.5*4.6
P_PPFD(:) = ZPFD(:)
!
! UPG*PT en attendat un calcul propre. Temperature des feuilles à l'ombre egale a la
! température de l'air. La temparature des feuilles au soleil egale a la valeur
! max entre la temperature de l'air et la temperaure radiative.
ZLSUT(:) = MAX(PLEAFT(:),PTEMP(:))
ZLSHT(:) = PTEMP(:)
!UPG*PT

!
! MEGAN : calcul des facteurs d'ajustement et de perte dans la canopée.
! ZCFSPEC: classe de sorties MEGAN (voir SPC_NOCONVER.EXT)
! 1: ISOP isoprene
! 2: MYRC myrcene
! 3: SABI sabinene
! 4: LIMO limonene 
! 5: A_3CAR carene_3
! 6: OCIM ocimene_t_b 
! 7: BPIN pinene_b
! 8: APIN pinene_a
! 9: OMTP A_2met_styrene + cymene_p + cymene_o + phellandrene_a + thujene_a + terpinene_a 
!          + terpinene_g + terpinolene + phellandrene_b + camphene + bornene + fenchene_a
!          + ocimene_al + ....
! 10: FARN
! 11: BCAR
! 12: OSQT
! 13: MBO
! 14: MEOH
! 15: ACTO
! 16: CO
! 17: NO
! 18: BIDER
! 19: STRESS
! 20: OTHER
!
!PRINT*, ZLAIC(:)

CALL EMPROC(ITIME, IDATE, ZD_PPFD, ZD_TEMP, ZDI, ZREC_ADJ, &
            GK%XLAT, GK%XLON, ZLAIP, ZLAIC, PTEMP,         &
            ZPFD,  PWIND, PPRES, PQV,  KSLTYP,             &
            PEK%XWG(:,5), PEK%XTG(:,1), PPFT,              &
            CHI%LSOILNOX, ZCFNO, ZCFNOG, ZCFSPEC, ZPPFD, ZWWILT, &
            ZGAMMA_SM, ZGAMMA_AGE, ZACTIVITY_FACTOR, PCO2, ZGAMMA_CO2)
            
!output the activation factors            
PGAMMA_SM(:)        =  ZGAMMA_SM(:)
PGAMMA_AGE(:)       =  ZGAMMA_AGE(:)
PACTIVITY_FACTOR(:) =  ZACTIVITY_FACTOR(:)
PGAMMA_CO2(:)       =  ZGAMMA_CO2(:)

!  
!
! MEGAN : calcul des flux d'émission 
! Dans cette partie du programme les sorties des 20 catégories obtenues à l'issu de la partie
!EMPROC sont multipliées par les valeurs des facteurs d'émissions correspondants, puis converties
!en 150 espèces, et associées en différentes catégories chimiques en fonction du schéma de chimie
!atmosphérique choisi parmi RADM2, RACM, SAPRCII, SAPRC99, CBMZ, SAPRC99X,
!SAPRC99Q, CB05, CB6, SOAX .
!
CALL MGN2MECH(IDATE, GK%XLAT, PEF, PPFT, ZCFNO, ZCFNOG, ZCFSPEC, &
              MGN%NSPMH_MAP, MGN%NMECH_MAP, MGN%XCONV_FAC,       &
              MGN%LCONVERSION, ZFLUX)

!
! Conversion ZFLUX from MEGAN mole/m2/s into molec/m2/s  
!ZFLUX(:,:) = ZFLUX(:,:) * XAVOGADRO
!
! Case of the same species between megan and mesonh
!DO JSV=1, SIZE(CHI%SVI%CSV)
 ! DO JSM=1, MGN%NVARS3D
  ! IF (TRIM(CHI%SVI%CSV(JSV)) == TRIM(MGN%CVNAME3D(JSM))) THEN
   !  print*, "TRUE"
    ! PSFTS(:,JSV) = PSFTS(:,JSV) + ZFLUX(JSM,:)
   !ELSE
    ! print*, 'FALSE'
   !END IF
  !END DO
!END DO
!



! Case of special treatment : ReLACS 1, 2, 3 scheme or CACM scheme
! Megan conversion is upon SOAX species
IF ( TRIM(MGN%CMECHANISM)=="RELACS" ) THEN
    PSFTS(:,MGN%NBIO ) = PSFTS(:,MGN%NBIO ) + (ZFLUX(MGN%NISOPRENE,:)* MM_ISO * 0.001) !mole/m2/s into Kg/m2/s
ENDIF
!
IF ( TRIM(MGN%CMECHANISM)=="RELACS2") THEN
  PSFTS(:,MGN%NORA1) = PSFTS(:,MGN%NORA1) + ZFLUX(MGN%NHCOOH,:) 
  PSFTS(:,MGN%NORA2) = PSFTS(:,MGN%NORA2) + ZFLUX(MGN%NCCO_OH,:) 
  PSFTS(:,MGN%NACID) = PSFTS(:,MGN%NACID) + ZFLUX(MGN%NRCO_OH,:)   
END IF
!
IF ( TRIM(MGN%CMECHANISM)=="CACM" ) THEN
  PSFTS(:,MGN%NACID) = PSFTS(:,MGN%NACID) + ZFLUX(MGN%NHCOOH,:) + ZFLUX(MGN%NCCO_OH,:) + ZFLUX(MGN%NRCO_OH,:)
ENDIF

IF ( TRIM(MGN%CMECHANISM)=="CACM".OR.TRIM(MGN%CMECHANISM)=="RELACS2" ) THEN
   PSFTS(:,MGN%NISOP) = PSFTS(:,MGN%NISOP) + ZFLUX(MGN%NISOPRENE,:) 
   PSFTS(:,MGN%NBIOH) = PSFTS(:,MGN%NBIOH) + 0.75*ZFLUX(MGN%NTRP1,:)
   PSFTS(:,MGN%NBIOL) = PSFTS(:,MGN%NBIOL) + 0.25*ZFLUX(MGN%NTRP1,:) 
   PSFTS(:,MGN%NKETL) = PSFTS(:,MGN%NKETL) + ZFLUX(MGN%NACET,:) + ZFLUX(MGN%NMEK,:) 
   PSFTS(:,MGN%NARAL) = PSFTS(:,MGN%NARAL) + ZFLUX(MGN%NBALD,:)
   PSFTS(:,MGN%NETHE) = PSFTS(:,MGN%NETHE) + ZFLUX(MGN%NETHENE,:)
   PSFTS(:,MGN%NALKL) = PSFTS(:,MGN%NALKL) + ZFLUX(MGN%NALK4,:)
   PSFTS(:,MGN%NALKM) = PSFTS(:,MGN%NALKM) + 0.5*ZFLUX(MGN%NALK5,:)
   PSFTS(:,MGN%NALKH) = PSFTS(:,MGN%NALKH) + 0.5*ZFLUX(MGN%NALK5,:)
   PSFTS(:,MGN%NAROH) = PSFTS(:,MGN%NAROH) + 0.5*ZFLUX(MGN%NARO1,:)
   PSFTS(:,MGN%NAROL) = PSFTS(:,MGN%NAROL) + 0.5*ZFLUX(MGN%NARO1,:)
   PSFTS(:,MGN%NAROO) = PSFTS(:,MGN%NAROO) + ZFLUX(MGN%NARO2,:)
   PSFTS(:,MGN%NOLEL) = PSFTS(:,MGN%NOLEL) + 0.5*ZFLUX(MGN%NOLE1,:)
   PSFTS(:,MGN%NOLEH) = PSFTS(:,MGN%NOLEH) + 0.5*ZFLUX(MGN%NOLE1,:)
END IF
!print*, PSFTS(:,MGN%NBIO)
!
!
#endif
END SUBROUTINE COUPLING_MEGAN_n
