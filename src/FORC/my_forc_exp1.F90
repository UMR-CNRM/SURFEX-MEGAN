SUBROUTINE MY_FORC_EXP1(HEXPER,KNI,KNPTS,              &
                   KYEAR,KMONTH,KDAY,PTIME,                        &
                   PLON, PLAT, PZS, PZREF, PUREF,                  &
                   PTA, PQA, PPS, PWINDSPEED, PWINDDIR,            &
                   PDIR_SW, PSCA_SW, PLW, PRAIN, PSNOW, PCO2       )

!----------------------------
!!
!!    PURPOSE
!!    -------
!!   This subroutine allows the user to build atm. forcing of his(her) run.
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson and P. Lemoigne                 Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     06/06
!!
!
!----------------------------------------------------------------------------
!      
!*    0.     Declaration of dummy arguments
!            ------------------------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

!----------------------------------------------------------------------------

CHARACTER(LEN=12), INTENT(IN) :: HEXPER    ! experiment name
INTEGER, INTENT(IN)          :: KNI       ! number of grid cells
INTEGER, INTENT(IN)          :: KNPTS     ! number of forcing instants
INTEGER, INTENT(OUT)         :: KYEAR     ! year  of simulation begining
INTEGER, INTENT(OUT)         :: KMONTH    ! month of simulation begining
INTEGER, INTENT(OUT)         :: KDAY      ! day   of simulation begining
REAL,    INTENT(OUT)         :: PTIME     ! time  of simulation begining (s)
REAL*4, DIMENSION(KNPTS,KNI), INTENT(OUT) :: PCO2      ! CO2 concentration (kg/m3) 
REAL*4, DIMENSION(KNPTS,KNI), INTENT(OUT) :: PDIR_SW   ! Solar direct   radiation (W/m2)
REAL*4, DIMENSION(KNPTS,KNI), INTENT(OUT) :: PSCA_SW   ! Solar diffused radiation (W/m2)
REAL*4, DIMENSION(KNPTS,KNI), INTENT(OUT) :: PLW       ! Longwave radiation (W/m2)
REAL*4, DIMENSION(KNPTS,KNI), INTENT(OUT) :: PWINDSPEED! Wind speed (m/s)
REAL*4, DIMENSION(KNPTS,KNI), INTENT(OUT) :: PWINDDIR  ! Wind dir. (deg. from N, clockwise)
REAL*4, DIMENSION(KNPTS,KNI), INTENT(OUT) :: PRAIN     ! rain rate (kg/m2/s)
REAL*4, DIMENSION(KNPTS,KNI), INTENT(OUT) :: PSNOW     ! snow rate (kg/m2/s)
REAL*4, DIMENSION(KNPTS,KNI), INTENT(OUT) :: PTA       ! temperature (K)
REAL*4, DIMENSION(KNPTS,KNI), INTENT(OUT) :: PQA       ! humidity (kg/kg)
REAL*4, DIMENSION(KNPTS,KNI), INTENT(OUT) :: PPS       ! pressure (Pa)
REAL*4, DIMENSION(KNI),       INTENT(OUT) :: PZREF     ! height of temperature forcing (m)
REAL*4, DIMENSION(KNI),       INTENT(OUT) :: PUREF     ! height of wind forcing (m)
REAL*4, DIMENSION(KNI),       INTENT(OUT) :: PZS       ! orography (m)
REAL, DIMENSION(KNI),       INTENT(OUT) :: PLON      ! longitude (degrees)
REAL, DIMENSION(KNI),       INTENT(OUT) :: PLAT      ! latitude  (degrees)
!
!*    1.     Declaration of user local variables
!            -----------------------------------
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Input file:
!
CHARACTER(LEN=*), PARAMETER       :: YFILE_FORCIN_01 = '../DATA/exp1/Q'  !Humidité
CHARACTER(LEN=*), PARAMETER       :: YFILE_FORCIN_02 = '../DATA/exp1/SSRD'  !rayonnement solaire incident
CHARACTER(LEN=*), PARAMETER       :: YFILE_FORCIN_03 = '../DATA/exp1/T'     !température
CHARACTER(LEN=*), PARAMETER       :: YFILE_FORCIN_04 = '../DATA/exp1/U'     
CHARACTER(LEN=*), PARAMETER       :: YFILE_FORCIN_05 = '../DATA/exp1/V'
CHARACTER(LEN=*), PARAMETER       :: YFILE_FORCIN_06 = '../DATA/exp1/CRR'   !précipitation convective
CHARACTER(LEN=*), PARAMETER       :: YFILE_FORCIN_07 = '../DATA/exp1/SP'    !pression
CHARACTER(LEN=*), PARAMETER       :: YFILE_FORCIN_08 = '../DATA/exp1/STRD'  !rayonnement solaire LW
CHARACTER(LEN=*), PARAMETER       :: YFILE_FORCIN_09 = '../DATA/exp1/LSRR'  !précipitation large scale
CHARACTER(LEN=*), PARAMETER       :: YFILE_FORCIN_10 = '../DATA/exp1/CSFR'  !neige convection
CHARACTER(LEN=*), PARAMETER       :: YFILE_FORCIN_11 = '../DATA/exp1/LSSFR' !neige large scale
CHARACTER(LEN=*), PARAMETER       :: YFILE_FORCIN_12 = '../DATA/exp1/Z'     !orographie
!       
!----------------------------------------------------------------------------
! Declare variables to translate forcing from obs to surfex atmospheric forcing
!
!
REAL, DIMENSION(KNPTS,KNI) :: ZTA, ZQA, ZUA, ZVA, ZSP, ZRG, ZPRECIP_CONV, ZPRECIP_LS, ZSNOW_CONV, ZSNOW_LS, ZLG
REAL, DIMENSION(KNI) :: ZS
				
!
REAL, DIMENSION(KNI) :: ZLON, ZLAT



INTEGER :: I, T, J, N! loop counters


REAL(KIND=JPRB) :: DUMMY
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
!*    2.     Initialization of date (UTC)
!            ------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_MY_FORC_EXP1:MY_FORC_EXP1',0,ZHOOK_HANDLE)

KDAY    = 01          ! starting day 

KMONTH  = 09          ! starting month

KYEAR   = 2018        ! starting year

PTIME   =    0.       ! starting time (s)
!
!-----------------------------------------------------------------------------

!      3.    grid definition
!            ---------------
!
PRINT*, 'YFILE_FORCIN=', YFILE_FORCIN_01
OPEN(UNIT=11, FILE=YFILE_FORCIN_01, FORM='FORMATTED', STATUS='OLD')
DO I =1, KNI
   READ(11,*)ZLAT(I),ZLON(I) ! N A TOUTES LES LONGITUTES ET LES LATITUDES DES FICHIERS DE FORCAGE
ENDDO

PLAT(:) = ZLAT(:)
PLON(:) = ZLON(:)

CLOSE(UNIT=11)
!
!----------------------------------------------------------------------------
!      
!        4.    orography definition
!               --------------------

PRINT*, 'YFILE_FORCIN_12=', YFILE_FORCIN_12
OPEN(UNIT=22, FILE=YFILE_FORCIN_12, FORM='FORMATTED', STATUS='OLD')

DO I=1,KNI
   READ(22,*)DUMMY,DUMMY,ZS(I)
ENDDO

CLOSE(UNIT=22)

!      
!-----------------------------------------------------------------------------
!      
!      5.    Forcing height
!            --------------
!
PZREF(:)   = 2.
PUREF(:)   = 10.
!
!----------------------------------------------------------------------------


!*      6.   Initialization of forcing variables
!            -----------------------------------
PRINT*, 'YFILE_FORCIN_1=', YFILE_FORCIN_01
OPEN(UNIT=11, FILE=YFILE_FORCIN_01, FORM='FORMATTED', STATUS='OLD')
PRINT*, 'YFILE_FORCIN_2=', YFILE_FORCIN_02
OPEN(UNIT=12, FILE=YFILE_FORCIN_02, FORM='FORMATTED', STATUS='OLD')
PRINT*, 'YFILE_FORCIN-3=', YFILE_FORCIN_03
OPEN(UNIT=13, FILE=YFILE_FORCIN_03, FORM='FORMATTED', STATUS='OLD')
PRINT*, 'YFILE_FORCIN_4=', YFILE_FORCIN_04
OPEN(UNIT=14, FILE=YFILE_FORCIN_04, FORM='FORMATTED', STATUS='OLD')
PRINT*, 'YFILE_FORCIN_5=', YFILE_FORCIN_05
OPEN(UNIT=15, FILE=YFILE_FORCIN_05, FORM='FORMATTED', STATUS='OLD')
PRINT*, 'YFILE_FORCIN_6=', YFILE_FORCIN_06
OPEN(UNIT=16, FILE=YFILE_FORCIN_06, FORM='FORMATTED', STATUS='OLD')
PRINT*, 'YFILE_FORCIN_7=', YFILE_FORCIN_07
OPEN(UNIT=17, FILE=YFILE_FORCIN_07, FORM='FORMATTED', STATUS='OLD')
PRINT*, 'YFILE_FORCIN_8=', YFILE_FORCIN_08
OPEN(UNIT=18, FILE=YFILE_FORCIN_08, FORM='FORMATTED', STATUS='OLD')
PRINT*, 'YFILE_FORCIN_9=', YFILE_FORCIN_09
OPEN(UNIT=19, FILE=YFILE_FORCIN_09, FORM='FORMATTED', STATUS='OLD')
PRINT*, 'YFILE_FORCIN_10=', YFILE_FORCIN_10
OPEN(UNIT=20, FILE=YFILE_FORCIN_10, FORM='FORMATTED', STATUS='OLD')
PRINT*, 'YFILE_FORCIN_11=', YFILE_FORCIN_11
OPEN(UNIT=21, FILE=YFILE_FORCIN_11, FORM='FORMATTED', STATUS='OLD')


DO T=1,KNPTS
   DO I=1,KNI
      READ(11,*)DUMMY,DUMMY,ZQA(T,I)
      READ(12,*)DUMMY,DUMMY,ZRG(T,I)   
      READ(13,*)DUMMY,DUMMY,ZTA(T,I)
      READ(14,*)DUMMY,DUMMY,ZUA(T,I)    
      READ(15,*)DUMMY,DUMMY,ZVA(T,I)
      READ(16,*)DUMMY,DUMMY,ZPRECIP_CONV(T,I)
      READ(17,*)DUMMY,DUMMY,ZSP(T,I)
      READ(18,*)DUMMY,DUMMY,ZLG(T,I)
      READ(19,*)DUMMY,DUMMY,ZPRECIP_LS(T,I)
      READ(20,*)DUMMY,DUMMY,ZSNOW_CONV(T,I)
      READ(21,*)DUMMY,DUMMY,ZSNOW_LS(T,I)
   ENDDO
ENDDO

DO T=1,KNPTS
   DO I=1,KNI
      IF (ZQA(T,I) .LT. 0._JPRB) ZQA(T,I) = 0.000001
   ENDDO
ENDDO
          
CLOSE(UNIT=11)
CLOSE(UNIT=12)
CLOSE(UNIT=13)
CLOSE(UNIT=14)
CLOSE(UNIT=15)
CLOSE(UNIT=16)
CLOSE(UNIT=17)
CLOSE(UNIT=18)
CLOSE(UNIT=19)
CLOSE(UNIT=20)
CLOSE(UNIT=21)

!
!---------------------------------------------------------------------------------------------------------------------------------
!
!        6. Fills Surfex forcing variables
!           ------------------------------
!
PZS(:) = ZS(:)
PCO2(:,:)    = 0.000620   ! (kg/m3, equivalent to 350 ppm) 
PDIR_SW(:,:) = ZRG(:,:)
PSCA_SW(:,:) = 0.
PWINDSPEED(:,:) = SQRT(ZVA(:,:)**2+ZUA(:,:)**2)
PWINDDIR  (:,:) = 57.29578*(ATAN(ZUA(:,:),ZVA(:,:)))+180.
PRAIN(:,:) =  ZPRECIP_CONV(:,:)+ZPRECIP_LS(:,:)
PSNOW = ZSNOW_CONV(:,:)+ZSNOW_LS(:,:)
PLW(:,:) = ZLG(:,:)
PTA(:,:) = ZTA(:,:)
PPS(:,:) = ZSP(:,:)
PQA(:,:) = ZQA(:,:)
print*,minval(PTA(:,:)),maxval(PTA(:,:))
print*,minval(PQA(:,:)),maxval(PQA(:,:))
!
IF (LHOOK) CALL DR_HOOK('MODI_MY_FORC_EXP1:MY_FORC_EXP1',1,ZHOOK_HANDLE)

!----------------------------------------------------------------------------
END SUBROUTINE MY_FORC_EXP1
!==============================================================================
