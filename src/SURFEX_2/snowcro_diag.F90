!     #########

SUBROUTINE SNOWCRO_DIAG(HSNOWMETAMO, &
                        PSNOWDZ, PSNOWSWE, PSNOWRHO, PSNOWGRAN1, PSNOWGRAN2, PSNOWAGE, &
                        PSNOWHIST, PSNOWTEMP, PSNOWLIQ, PDIRCOSZW, PSNOWDEND, PSNOWSPHER, &
                        PSNOWSIZE, PSNOWSSA, PSNOWTYPEMEPRA, PSNOWRAM, PSNOWSHEAR, &
                        PSNOWDEPTH_1DAYS, PSNOWDEPTH_3DAYS, PSNOWDEPTH_5DAYS, &
                        PSNOWDEPTH_7DAYS, PSNOWSWE_1DAYS, PSNOWSWE_3DAYS, PSNOWSWE_5DAYS,&
                        PSNOWSWE_7DAYS, PSNOWRAM_SONDE, PSNOW_WETTHICKNESS, PSNOW_REFROZENTHICKNESS)

! Diagnostics of Crocus snowpack model
! Author: M. Lafaysse, Meteo-France, October 2015

USE MODD_SURF_PAR,      ONLY : XUNDEF

USE MODD_CSTS,ONLY : XRHOLI, XRHOLW

USE MODD_SNOW_PAR,ONLY : ICRIS_DEND1D, ICRIS_NONDEND1D, &
                         IFR, IFR_LB, ILB, ILB_FIN, ILB_ANG, IROUL, IFIN, IFIN_AR, IFIN_ANG, &
                         IPL, IPL_GOB, IGOB, IGEL, IGOB_FON, IRON_ANG, XX, XD1, XD2, XD3

IMPLICIT NONE

CHARACTER(3), INTENT(IN)         :: HSNOWMETAMO ! metamorphism option
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ
REAL, DIMENSION(:,:), INTENT(IN) :: PSNOWSWE
REAL, DIMENSION(:,:), INTENT(IN) :: PSNOWRHO
REAL, DIMENSION(:,:), INTENT(IN) :: PSNOWGRAN1
REAL, DIMENSION(:,:), INTENT(IN) :: PSNOWGRAN2
REAL, DIMENSION(:,:), INTENT(IN) :: PSNOWAGE
REAL, DIMENSION(:,:), INTENT(IN) :: PSNOWHIST
REAL, DIMENSION(:,:), INTENT(IN) :: PSNOWTEMP
REAL, DIMENSION(:,:), INTENT(IN) :: PSNOWLIQ
REAL, DIMENSION(:),   INTENT(IN) :: PDIRCOSZW !cosine of slope
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWDEND
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWSPHER
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWSIZE
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWSSA
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWTYPEMEPRA
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWRAM
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWSHEAR
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWDEPTH_1DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWDEPTH_3DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWDEPTH_5DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWDEPTH_7DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWSWE_1DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWSWE_3DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWSWE_5DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWSWE_7DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWRAM_SONDE
REAL, DIMENSION(:), INTENT(OUT) :: PSNOW_WETTHICKNESS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOW_REFROZENTHICKNESS

REAL :: ZG1, ZG2, ZRFIN, ZRDEN, ZRFGF

REAL :: ZDIAM

LOGICAL,DIMENSION(SIZE(PSNOWSWE,1)) :: GRAM, GWET, GREFROZEN

INTEGER :: ICLASS_DEND, ICLASS_SPHER, ICLASS_SIZE, ICLASS_HIST
INTEGER :: ICLASS
LOGICAL :: LTHERM

INTEGER :: JJ,JST

! PRINT*,ASSOCIATED(PSNOWDEPTH_1DAYS),ASSOCIATED(PSNOWDEPTH_3DAYS)
! PRINT*,ASSOCIATED(PSNOWDEPTH_5DAYS),ASSOCIATED(PSNOWDEPTH_7DAYS)
! PRINT*,"in snowcrodiag"
! PRINT*,ALLOCATED(PSNOWDEPTH_1DAYS),ALLOCATED(PSNOWDEPTH_3DAYS)
! PRINT*,ALLOCATED(PSNOWDEPTH_5DAYS),ALLOCATED(PSNOWDEPTH_7DAYS)
! PRINT*,SIZE(PSNOWDEPTH_1DAYS),SIZE(PSNOWDEPTH_3DAYS)
! PRINT*,SIZE(PSNOWDEPTH_5DAYS),SIZE(PSNOWDEPTH_7DAYS)

! Initializations

PSNOWDEND      = XUNDEF
PSNOWSPHER     = XUNDEF
PSNOWSIZE      = XUNDEF
PSNOWSSA       = XUNDEF
PSNOWTYPEMEPRA = XUNDEF
PSNOWRAM       = XUNDEF
PSNOWSHEAR     = XUNDEF

PSNOWDEPTH_1DAYS = 0.
PSNOWDEPTH_3DAYS = 0.
PSNOWDEPTH_5DAYS = 0.
PSNOWDEPTH_7DAYS = 0.
PSNOWSWE_1DAYS = 0.
PSNOWSWE_3DAYS = 0.
PSNOWSWE_5DAYS = 0.
PSNOWSWE_7DAYS = 0.
PSNOWRAM_SONDE = 0.
PSNOW_WETTHICKNESS = 0.
PSNOW_REFROZENTHICKNESS = 0.
GRAM = .TRUE.
GWET = .TRUE.
GREFROZEN = .TRUE.


DO JST=1,SIZE(PSNOWSWE,2)
  DO JJ=1,SIZE(PSNOWSWE,1)
    
    IF (PSNOWSWE(JJ,JST)>0) THEN
    
      ! In this routine Crocus diagnostics are perpendicular to the slope.
      ! The projection is done in diag_misc_isban             
    
      ZG1 = PSNOWGRAN1(JJ,JST)/99.
      ZRFIN = 0.17*PSNOWRHO(JJ,JST)-31
    
      IF (PSNOWGRAN1(JJ,JST)>=0) THEN
        !Non dendritic case

        !Dendricity,sphericty and grain size
        PSNOWSIZE(JJ,JST)  = PSNOWGRAN2(JJ,JST)
        PSNOWDEND(JJ,JST)  = 0
        PSNOWSPHER(JJ,JST) = PSNOWGRAN1(JJ,JST) / XX

        !Optical diameter for SSA diagnostic
        ZDIAM = PSNOWSIZE(JJ,JST) * PSNOWSPHER(JJ,JST) + &
        MAX( 0.0004, 0.5*PSNOWSIZE(JJ,JST) ) * ( 1.-PSNOWSPHER(JJ,JST) )

        !10 classes of sphericity 0:[0,0.05[, 1:[0.05,0.15[, ..., 9:[0.85,1.0]
        !###########Strange way of defining sphericity classes -> Check with very old versions
        !###########ICLASS_SPHER = MIN(NINT(10 * PSNOWSPHER(JJ,JST)),9) Not exactly the same
        ICLASS_SPHER = MIN(INT(10 * PSNOWSPHER(JJ,JST) + 0.05),9)


        !6 classes of historical variable {0,1,...,5}
        !##########Why is PSNOWHIST stored as float64. It always takes only 6 integer values
        ICLASS_HIST = NINT(PSNOWHIST(JJ,JST))

        !3 classes of grain size in mm 0:[0,0.55[, 1:[0.55,1.05[, 2:[1.05, +inf[
        !#########Strange +0.05
        IF (PSNOWSIZE(JJ,JST) < 0.00055) THEN
          ICLASS_SIZE = 0
        ELSEIF (PSNOWSIZE(JJ,JST) < 0.00105) THEN
          ICLASS_SIZE = 1
        ELSE
          ICLASS_SIZE = 2
        ENDIF

        !Overall 10x3x6 classes from 1 to 180 (included)
        ICLASS = 1 + ICLASS_SPHER + ICLASS_SIZE*10 + ICLASS_HIST*30

        !Snow type obtained in table ICRIS_NONDEND1D
        PSNOWTYPEMEPRA(JJ,JST) = ICRIS_NONDEND1D(ICLASS)


        ! Ram resistance (non dendritic case)
        ! PSNOWLIQ tel massique en m
        ! PSNOWLIQ*XHROLW/PSNOWDZ tel volumique en kg m-3
        ! PSNOWLIQ/PSNOWDZ : rapport sans unité
        ! Seuil à 0.5% soit 0.005
        LTHERM=((PSNOWTEMP(JJ,JST)<272.96).OR.(PSNOWLIQ(JJ,JST)/PSNOWDZ(JJ,JST)<=0.005))

        SELECT CASE (NINT(PSNOWTYPEMEPRA(JJ,JST)))
        CASE (IFIN)
          PSNOWRAM(JJ,JST)=MAX(3.,ZRFIN)
        CASE (IFIN_ANG)
          IF (PSNOWRHO(JJ,JST)<200) THEN
            PSNOWRAM(JJ,JST)=ZRFIN*PSNOWSPHER(JJ,JST)+&
            (1- PSNOWSPHER(JJ,JST))*(ZRFIN*(0.8-PSNOWSIZE(JJ,JST))+2*PSNOWSIZE(JJ,JST))
          ELSE
            PSNOWRAM(JJ,JST)=2
          ENDIF
        CASE (IFIN_AR,IGEL,IGOB_FON,IRON_ANG)
      
          IF (LTHERM) THEN
            PSNOWRAM(JJ,JST)=MAX(10.,0.103*PSNOWRHO(JJ,JST)-19.666)
          ELSE
            IF (PSNOWRHO(JJ,JST)<250) THEN
              PSNOWRAM(JJ,JST)=1
            ELSE
              PSNOWRAM(JJ,JST)=MAX(2.,0.16*PSNOWRHO(JJ,JST)-54)
            ENDIF
          END IF
      
        CASE (IPL,IPL_GOB)
          IF (PSNOWSIZE(JJ,JST)>0.8) THEN
            PSNOWRAM(JJ,JST)=MAX(3.,ZRFIN)*(0.8-PSNOWSIZE(JJ,JST))+2*PSNOWSIZE(JJ,JST)
          ELSE
            PSNOWRAM(JJ,JST)=2
          ENDIF
        CASE DEFAULT
        END SELECT

      ELSE
        !Dendritic case

        !Dendricity,sphericty and grain size
        PSNOWSIZE(JJ,JST)  =  XUNDEF  !Grain size not defined for dendritic snow
        PSNOWDEND(JJ,JST)  = -PSNOWGRAN1(JJ,JST) / XX
        PSNOWSPHER(JJ,JST) =  PSNOWGRAN2(JJ,JST) / XX

        !Optical diameter for SSA diagnostic
        ZDIAM = PSNOWDEND(JJ,JST) * XD1 + (1 - PSNOWDEND(JJ,JST)) * &
        (PSNOWSPHER(JJ,JST) * XD2 + (1 - PSNOWSPHER(JJ,JST)) * XD3)
        !ZDIAM =  -PSNOWGRAN1(JJ,JST)*XD1/XX + (1.+PSNOWGRAN1(JJ,JST)/XX) * &
        !      ( PSNOWGRAN2(JJ,JST)*XD2/XX + (1.-PSNOWGRAN2(JJ,JST)/XX) * XD3 )
        ZDIAM = ZDIAM/10000.

        !10 classes of dendricity 0:[0,0.1[, ..., 9:[0.9,1.0[ (value 1.0 does not exist)
        ICLASS_DEND = INT(10 * PSNOWDEND(JJ,JST))

        !10 classes of sphericity 0:[0,0.05[, 1:[0.05,0.15[, ..., 9:[0.85,1.0]
        !###########Strange way of defining sphericity classes -> Check with very old versions
        !###########ICLASS_SPHER = MIN(NINT(10 * PSNOWSPHER(JJ,JST)),9) Not exactly the same
        ICLASS_SPHER = MIN(INT(10 * PSNOWSPHER(JJ,JST) + 0.05),9)

        !Overall 10x10 classes from 1 to 100 (included)
        ICLASS = 1 + ICLASS_DEND + ICLASS_SPHER*10

        !Snow type obtained in table ICRIS_DEND1D
        PSNOWTYPEMEPRA(JJ,JST) = ICRIS_DEND1D(ICLASS)

      ENDIF

      ! All cases
      ! Compute depth and SWE of recent snow
      IF(PSNOWAGE(JJ,JST)<=1)THEN
        PSNOWDEPTH_1DAYS(JJ) = PSNOWDEPTH_1DAYS(JJ) + PSNOWDZ (JJ,JST)
        PSNOWSWE_1DAYS  (JJ) = PSNOWSWE_1DAYS  (JJ) + PSNOWSWE(JJ,JST)
      ENDIF

      IF(PSNOWAGE(JJ,JST)<=3)THEN
        PSNOWDEPTH_3DAYS(JJ) = PSNOWDEPTH_3DAYS(JJ) + PSNOWDZ (JJ,JST)    
        PSNOWSWE_3DAYS  (JJ) = PSNOWSWE_3DAYS  (JJ) + PSNOWSWE(JJ,JST)
      ENDIF
    
      IF(PSNOWAGE(JJ,JST)<=5)THEN
        PSNOWDEPTH_5DAYS(JJ) = PSNOWDEPTH_5DAYS(JJ) + PSNOWDZ (JJ,JST)    
        PSNOWSWE_5DAYS  (JJ) = PSNOWSWE_5DAYS  (JJ) + PSNOWSWE(JJ,JST)
      ENDIF
  
      IF(PSNOWAGE(JJ,JST)<=7)THEN
        PSNOWDEPTH_7DAYS(JJ) = PSNOWDEPTH_7DAYS(JJ) + PSNOWDZ (JJ,JST)    
        PSNOWSWE_7DAYS  (JJ) = PSNOWSWE_7DAYS  (JJ) + PSNOWSWE(JJ,JST)
      END IF
    
      ! Ram sonde penetration
      IF ((GRAM(JJ)).AND.(PSNOWRAM(JJ,JST)<=2.)) THEN
        PSNOWRAM_SONDE(JJ)=PSNOWRAM_SONDE(JJ)+PSNOWDZ(JJ,JST)
      ELSE
        GRAM(JJ)=.FALSE.
      ENDIF

      ! Depth of wet snow
      IF ((GWET(JJ)).AND.(PSNOWLIQ(JJ,JST)>0.)) THEN
        PSNOW_WETTHICKNESS(JJ)=PSNOW_WETTHICKNESS(JJ)+PSNOWDZ(JJ,JST)
      ELSE
        GWET(JJ)=.FALSE.
      ENDIF
      ! Depth of refrozen snow
      IF ((GREFROZEN(JJ)).AND.(PSNOWHIST(JJ,JST)>=2).AND.(PSNOWTEMP(JJ,JST)<273.15)) THEN
        PSNOW_REFROZENTHICKNESS(JJ)=PSNOW_REFROZENTHICKNESS(JJ)+PSNOWDZ(JJ,JST)
      ELSE
        GREFROZEN(JJ)=.FALSE.
      ENDIF    
    
      ! Specific surface area
      IF ( HSNOWMETAMO=='B92' ) THEN
        PSNOWSSA(JJ,JST) = 6. / (XRHOLI*ZDIAM)
      ELSE
        PSNOWSSA(JJ,JST) = 6. / (XRHOLI*PSNOWGRAN1(JJ,JST))
      END IF

    END IF
  END DO
END DO
  
END SUBROUTINE SNOWCRO_DIAG
