
































MODULE saprc99_mosaic_4bin_vbs2_Integrator

 USE saprc99_mosaic_4bin_vbs2_Parameters
 USE saprc99_mosaic_4bin_vbs2_Precision
 USE saprc99_mosaic_4bin_vbs2_JacobianSP

  IMPLICIT NONE
 


  INTEGER, PARAMETER :: ifun=1, ijac=2, istp=3, iacc=4, &
    irej=5, idec=6, isol=7, isng=8, itexit=1, ihexit=2
    

  
  CHARACTER(LEN=50), PARAMETER, DIMENSION(-8:1) :: IERR_NAMES = (/ &
    'Matrix is repeatedly singular                     ', & 
    'Step size too small                               ', & 
    'No of steps exceeds maximum bound                 ', & 
    'Improper tolerance values                         ', & 
    'FacMin/FacMax/FacRej must be positive             ', & 
    'Hmin/Hmax/Hstart must be positive                 ', & 
    'Selected Rosenbrock method not implemented        ', & 
    'Improper value for maximal no of steps            ', & 
    '                                                  ', & 
    'Success                                           ' /) 

CONTAINS

SUBROUTINE  saprc99_mosaic_4bin_vbs2_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE saprc99_mosaic_4bin_vbs2_Parameters

   IMPLICIT NONE
   REAL(kind=dp), INTENT(INOUT), DIMENSION(NFIX) :: FIX
   REAL(kind=dp), INTENT(INOUT), DIMENSION(NVAR) :: VAR
   REAL(kind=dp), INTENT(INOUT) :: IRR_WRK(NREACT)
   REAL(kind=dp), INTENT(IN), DIMENSION(NSPEC) :: ATOL, RTOL
   REAL(kind=dp), INTENT(IN), DIMENSION(NREACT) :: RCONST
   REAL(kind=dp), INTENT(IN) :: TIN  
   REAL(kind=dp), INTENT(IN) :: TOUT 
   
   INTEGER,  INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
   REAL(kind=dp), INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
   INTEGER,  INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
   REAL(kind=dp), INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
   INTEGER,  INTENT(OUT), OPTIONAL :: IERR_U

   REAL(kind=dp)  :: STEPMIN


   INTEGER :: N_stp, N_acc, N_rej, N_sng
   SAVE N_stp, N_acc, N_rej, N_sng
   INTEGER :: i, IERR
   REAL(kind=dp) :: RCNTRL(20), RSTATUS(20)
   INTEGER :: ICNTRL(20), ISTATUS(20)


   ICNTRL(:)  = 0
   RCNTRL(:)  = 0.0_dp
   ISTATUS(:) = 0
   RSTATUS(:) = 0.0_dp

   
   
   IF (PRESENT(ICNTRL_U)) THEN
     WHERE(ICNTRL_U(:) > 0) ICNTRL(:) = ICNTRL_U(:)
   END IF
   IF (PRESENT(RCNTRL_U)) THEN
     WHERE(RCNTRL_U(:) > 0) RCNTRL(:) = RCNTRL_U(:)
   END IF

   CALL saprc99_mosaic_4bin_vbs2_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  saprc99_mosaic_4bin_vbs2_INTEGRATE


SUBROUTINE  saprc99_mosaic_4bin_vbs2_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE saprc99_mosaic_4bin_vbs2_Parameters

  IMPLICIT NONE


   REAL(kind=dp), INTENT(INOUT) :: Y(NVAR)
   REAL(kind=dp), INTENT(INOUT) :: IRR_WRK(NREACT)
   REAL(kind=dp), INTENT(IN), DIMENSION(NFIX) :: FIX
   REAL(kind=dp), INTENT(IN), DIMENSION(NREACT) :: RCONST
   REAL(kind=dp), INTENT(IN)   :: Tstart,Tend
   REAL(kind=dp), INTENT(IN)   :: AbsTol(NVAR),RelTol(NVAR)
   INTEGER, INTENT(IN)    :: ICNTRL(20)
   REAL(kind=dp), INTENT(IN)   :: RCNTRL(20)
   INTEGER, INTENT(INOUT) :: ISTATUS(20)
   REAL(kind=dp), INTENT(INOUT) :: RSTATUS(20)
   INTEGER, INTENT(OUT)   :: IERR

   INTEGER, PARAMETER :: Smax = 6
   INTEGER  :: Method, ros_S
   REAL(kind=dp), DIMENSION(Smax) :: ros_M, ros_E, ros_Alpha, ros_Gamma
   REAL(kind=dp), DIMENSION(Smax*(Smax-1)/2) :: ros_A, ros_C
   REAL(kind=dp) :: ros_ELO
   LOGICAL, DIMENSION(Smax) :: ros_NewF
   CHARACTER(LEN=12) :: ros_Name


  INTEGER :: Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng



   REAL(kind=dp) :: Roundoff, FacMin, FacMax, FacRej, FacSafe
   REAL(kind=dp) :: Hmin, Hmax, Hstart, Hexit
   REAL(kind=dp) :: Texit
   INTEGER :: i, UplimTol, Max_no_steps
   LOGICAL :: Autonomous, VectorTol

   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp


   Nfun = ISTATUS(ifun)
   Njac = ISTATUS(ijac)
   Nstp = ISTATUS(istp)
   Nacc = ISTATUS(iacc)
   Nrej = ISTATUS(irej)
   Ndec = ISTATUS(idec)
   Nsol = ISTATUS(isol)
   Nsng = ISTATUS(isng)


   Autonomous = .NOT.(ICNTRL(1) == 0)







   IF (ICNTRL(2) == 0) THEN
      VectorTol = .TRUE.
         UplimTol  = NVAR
   ELSE
      VectorTol = .FALSE.
         UplimTol  = 1
   END IF


   IF (ICNTRL(3) == 0) THEN
      Method = 4
   ELSEIF ( (ICNTRL(3) >= 1).AND.(ICNTRL(3) <= 5) ) THEN
      Method = ICNTRL(3)
   ELSE
      PRINT * , 'User-selected Rosenbrock method: ICNTRL(3)=', Method
      CALL saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = saprc99_mosaic_4bin_vbs2_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL saprc99_mosaic_4bin_vbs2_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL saprc99_mosaic_4bin_vbs2_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL saprc99_mosaic_4bin_vbs2_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL saprc99_mosaic_4bin_vbs2_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL saprc99_mosaic_4bin_vbs2_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL saprc99_mosaic_4bin_vbs2_ros_Integrator(Y,Tstart,Tend,Texit,      &
        AbsTol, RelTol,                          &

        ros_S, ros_M, ros_E, ros_A, ros_C,       &
        ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, &

        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart, Hexit,     &
        FacMin, FacMax, FacRej, FacSafe,         &

        IRR_WRK,IERR,                            &

         Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng,&

         RCONST, FIX &
)



   ISTATUS(ifun) = Nfun
   ISTATUS(ijac) = Njac
   ISTATUS(istp) = Nstp
   ISTATUS(iacc) = Nacc
   ISTATUS(irej) = Nrej
   ISTATUS(idec) = Ndec
   ISTATUS(isol) = Nsol
   ISTATUS(isng) = Nsng

   RSTATUS(itexit) = Texit
   RSTATUS(ihexit) = Hexit


CONTAINS 



 SUBROUTINE  saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(Code,T,H,IERR)



   USE saprc99_mosaic_4bin_vbs2_Precision

   REAL(kind=dp), INTENT(IN) :: T, H
   INTEGER, INTENT(IN)  :: Code
   INTEGER, INTENT(OUT) :: IERR

   IERR = Code
   PRINT * , &
     'Forced exit from Rosenbrock due to the following error:'
   IF ((Code>=-8).AND.(Code<=-1)) THEN
     PRINT *, IERR_NAMES(Code)
   ELSE
     PRINT *, 'Unknown Error code: ', Code
   ENDIF

   PRINT *, "T=", T, "and H=", H

 END SUBROUTINE  saprc99_mosaic_4bin_vbs2_ros_ErrorMsg


 SUBROUTINE  saprc99_mosaic_4bin_vbs2_ros_Integrator (Y, Tstart, Tend, T,     &
        AbsTol, RelTol,                          &

        ros_S, ros_M, ros_E, ros_A, ros_C,       &
        ros_Alpha, ros_Gamma, ros_ELO, ros_NewF, &

        Autonomous, VectorTol, Max_no_steps,     &
        Roundoff, Hmin, Hmax, Hstart, Hexit,     &
        FacMin, FacMax, FacRej, FacSafe,         &

        IRR_WRK,IERR,                            &

        Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng, &

        RCONST, FIX )






  IMPLICIT NONE


   REAL(kind=dp), INTENT(INOUT) :: Y(NVAR)

   REAL(kind=dp), INTENT(INOUT) :: IRR_WRK(NREACT)

   REAL(kind=dp), INTENT(IN) :: Tstart,Tend

   REAL(kind=dp), INTENT(OUT) ::  T

   REAL(kind=dp), INTENT(IN) ::  AbsTol(NVAR), RelTol(NVAR)

   INTEGER, INTENT(IN) ::  ros_S
   REAL(kind=dp), INTENT(IN) :: ros_M(ros_S), ros_E(ros_S),  &
       ros_Alpha(ros_S), ros_A(ros_S*(ros_S-1)/2), &
       ros_Gamma(ros_S), ros_C(ros_S*(ros_S-1)/2), ros_ELO
   LOGICAL, INTENT(IN) :: ros_NewF(ros_S)

   LOGICAL, INTENT(IN) :: Autonomous, VectorTol
   REAL(kind=dp), INTENT(IN) :: Hstart, Hmin, Hmax
   INTEGER, INTENT(IN) :: Max_no_steps
   REAL(kind=dp), INTENT(IN) :: Roundoff, FacMin, FacMax, FacRej, FacSafe

   REAL(kind=dp), INTENT(OUT) :: Hexit

   INTEGER, INTENT(OUT) :: IERR

   REAL(kind=dp), INTENT(IN), DIMENSION(NFIX) :: FIX

   REAL(kind=dp), INTENT(IN), DIMENSION(NREACT) :: RCONST


  INTEGER, INTENT(INOUT)  :: Nfun,Njac,Nstp,Nacc,Nrej,Ndec,Nsol,Nsng


   REAL(kind=dp) :: Ynew(NVAR), Fcn0(NVAR), Fcn(NVAR)
   REAL(kind=dp) :: K(NVAR*ros_S), dFdT(NVAR)
   REAL(kind=dp) :: Jac0(LU_NONZERO), Ghimj(LU_NONZERO)
   REAL(kind=dp) :: H, Hnew, HC, HG, Fac, Tau
   REAL(kind=dp) :: Err, Yerr(NVAR)
   INTEGER :: Pivot(NVAR), Direction, ioffset, j, istage
   LOGICAL :: RejectLastH, RejectMoreH, Singular

   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp, ONE  = 1.0_dp
   REAL(kind=dp), PARAMETER :: DeltaMin = 1.0E-5_dp







   T = Tstart
   Hexit = 0.0_dp
   H = MIN(Hstart,Hmax)
   IF (ABS(H) <= 10.0_dp*Roundoff) H = DeltaMin

   IF (Tend  >=  Tstart) THEN
     Direction = +1
   ELSE
     Direction = -1
   END IF

   RejectLastH=.FALSE.
   RejectMoreH=.FALSE.



TimeLoop: DO WHILE ( (Direction > 0).AND.((T-Tend)+Roundoff <= ZERO) &
       .OR. (Direction < 0).AND.((Tend-T)+Roundoff <= ZERO) )

   IF ( Nstp > Max_no_steps ) THEN  
      CALL saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL saprc99_mosaic_4bin_vbs2_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL saprc99_mosaic_4bin_vbs2_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL saprc99_mosaic_4bin_vbs2_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL saprc99_mosaic_4bin_vbs2_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL saprc99_mosaic_4bin_vbs2_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL saprc99_mosaic_4bin_vbs2_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL saprc99_mosaic_4bin_vbs2_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL saprc99_mosaic_4bin_vbs2_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL saprc99_mosaic_4bin_vbs2_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL saprc99_mosaic_4bin_vbs2_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL saprc99_mosaic_4bin_vbs2_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL saprc99_mosaic_4bin_vbs2_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL saprc99_mosaic_4bin_vbs2_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL saprc99_mosaic_4bin_vbs2_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL saprc99_mosaic_4bin_vbs2_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL saprc99_mosaic_4bin_vbs2_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL saprc99_mosaic_4bin_vbs2_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL saprc99_mosaic_4bin_vbs2_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = saprc99_mosaic_4bin_vbs2_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL saprc99_mosaic_4bin_vbs2_WCOPY(NVAR,Ynew,1,Y,1)
      T = T + Direction*H
      Hnew = MAX(Hmin,MIN(Hnew,Hmax))
      IF (RejectLastH) THEN  
         Hnew = MIN(Hnew,H)
      END IF
      RejectLastH = .FALSE.
      RejectMoreH = .FALSE.
      H = Hnew
      EXIT UntilAccepted 
   ELSE           
      IF (RejectMoreH) THEN
         Hnew = H*FacRej
      END IF
      RejectMoreH = RejectLastH
      RejectLastH = .TRUE.
      H = Hnew
      IF (Nacc >= 1) THEN
         Nrej = Nrej+1
      END IF
   END IF 

   END DO UntilAccepted

   END DO TimeLoop


   IERR = 1  

  END SUBROUTINE  saprc99_mosaic_4bin_vbs2_ros_Integrator



  REAL(kind=dp) FUNCTION  saprc99_mosaic_4bin_vbs2_ros_ErrorNorm ( Y, Ynew, Yerr, &
               AbsTol, RelTol, VectorTol )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: Y(NVAR), Ynew(NVAR), &
          Yerr(NVAR), AbsTol(NVAR), RelTol(NVAR)
   LOGICAL, INTENT(IN) ::  VectorTol

   REAL(kind=dp) :: Err, Scale, Ymax
   INTEGER  :: i
   REAL(kind=dp), PARAMETER :: ZERO = 0.0_dp

   Err = ZERO
   DO i=1,NVAR
     Ymax = MAX(ABS(Y(i)),ABS(Ynew(i)))
     IF (VectorTol) THEN
       Scale = AbsTol(i)+RelTol(i)*Ymax
     ELSE
       Scale = AbsTol(1)+RelTol(1)*Ymax
     END IF
     Err = Err+(Yerr(i)/Scale)**2
   END DO
   Err  = SQRT(Err/NVAR)

    saprc99_mosaic_4bin_vbs2_ros_ErrorNorm = Err

  END FUNCTION  saprc99_mosaic_4bin_vbs2_ros_ErrorNorm



  SUBROUTINE saprc99_mosaic_4bin_vbs2_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL saprc99_mosaic_4bin_vbs2_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL saprc99_mosaic_4bin_vbs2_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL saprc99_mosaic_4bin_vbs2_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  saprc99_mosaic_4bin_vbs2_ros_FunTimeDeriv



  SUBROUTINE  saprc99_mosaic_4bin_vbs2_ros_PrepareMatrix ( H, Direction, gam, &
             Jac0, Ghimj, Pivot, Singular, Ndec,  Nsng  )








   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) ::  Jac0(LU_NONZERO)
   REAL(kind=dp), INTENT(IN) ::  gam
   INTEGER, INTENT(IN) ::  Direction

   REAL(kind=dp), INTENT(OUT) :: Ghimj(LU_NONZERO)
   LOGICAL, INTENT(OUT) ::  Singular
   INTEGER, INTENT(OUT) ::  Pivot(NVAR)

   REAL(kind=dp), INTENT(INOUT) :: H   
   INTEGER, INTENT(INOUT) ::  Ndec, Nsng

   INTEGER  :: i, ising, Nconsecutive
   REAL(kind=dp) :: ghinv
   REAL(kind=dp), PARAMETER :: ONE  = 1.0_dp, HALF = 0.5_dp

   Nconsecutive = 0
   Singular = .TRUE.

   DO WHILE (Singular)


     CALL saprc99_mosaic_4bin_vbs2_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL saprc99_mosaic_4bin_vbs2_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL saprc99_mosaic_4bin_vbs2_ros_Decomp( Ghimj, Pivot, ising, Ndec )
     IF (ising == 0) THEN

        Singular = .FALSE.
     ELSE 

        Nsng = Nsng+1
        Nconsecutive = Nconsecutive+1
        Singular = .TRUE.
        PRINT*,'Warning: LU Decomposition returned ising = ',ising
        IF (Nconsecutive <= 5) THEN 
           H = H*HALF
        ELSE  
           RETURN
        END IF  
      END IF    

   END DO 

  END SUBROUTINE  saprc99_mosaic_4bin_vbs2_ros_PrepareMatrix



  SUBROUTINE  saprc99_mosaic_4bin_vbs2_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 



CALL decomp_saprc99_mosaic_4bin_vbs2 ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  saprc99_mosaic_4bin_vbs2_ros_Decomp



  SUBROUTINE  saprc99_mosaic_4bin_vbs2_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL saprc99_mosaic_4bin_vbs2_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  saprc99_mosaic_4bin_vbs2_ros_Solve




  SUBROUTINE  saprc99_mosaic_4bin_vbs2_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
            ros_Gamma,ros_NewF,ros_ELO,ros_Name)




  IMPLICIT NONE

   INTEGER, PARAMETER :: S=2
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name

    REAL(kind=dp) :: g

    g = 1.0_dp + 1.0_dp/SQRT(2.0_dp)


    ros_Name = 'ROS-2'

    ros_S = S








    ros_A(1) = (1.0_dp)/g
    ros_C(1) = (-2.0_dp)/g


    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.

    ros_M(1)= (3.0_dp)/(2.0_dp*g)
    ros_M(2)= (1.0_dp)/(2.0_dp*g)

    ros_E(1) = 1.0_dp/(2.0_dp*g)
    ros_E(2) = 1.0_dp/(2.0_dp*g)


    ros_ELO = 2.0_dp

    ros_Alpha(1) = 0.0_dp
    ros_Alpha(2) = 1.0_dp

    ros_Gamma(1) = g
    ros_Gamma(2) =-g

 END SUBROUTINE  saprc99_mosaic_4bin_vbs2_Ros2



  SUBROUTINE  saprc99_mosaic_4bin_vbs2_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
           ros_Gamma,ros_NewF,ros_ELO,ros_Name)




  IMPLICIT NONE

   INTEGER, PARAMETER :: S=3
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name


   ros_Name = 'ROS-3'

   ros_S = S








   ros_A(1)= 1.0_dp
   ros_A(2)= 1.0_dp
   ros_A(3)= 0.0_dp

   ros_C(1) = -0.10156171083877702091975600115545E+01_dp
   ros_C(2) =  0.40759956452537699824805835358067E+01_dp
   ros_C(3) =  0.92076794298330791242156818474003E+01_dp


   ros_NewF(1) = .TRUE.
   ros_NewF(2) = .TRUE.
   ros_NewF(3) = .FALSE.

   ros_M(1) =  0.1E+01_dp
   ros_M(2) =  0.61697947043828245592553615689730E+01_dp
   ros_M(3) = -0.42772256543218573326238373806514E+00_dp

   ros_E(1) =  0.5E+00_dp
   ros_E(2) = -0.29079558716805469821718236208017E+01_dp
   ros_E(3) =  0.22354069897811569627360909276199E+00_dp


   ros_ELO = 3.0_dp

   ros_Alpha(1)= 0.0E+00_dp
   ros_Alpha(2)= 0.43586652150845899941601945119356E+00_dp
   ros_Alpha(3)= 0.43586652150845899941601945119356E+00_dp

   ros_Gamma(1)= 0.43586652150845899941601945119356E+00_dp
   ros_Gamma(2)= 0.24291996454816804366592249683314E+00_dp
   ros_Gamma(3)= 0.21851380027664058511513169485832E+01_dp

  END SUBROUTINE  saprc99_mosaic_4bin_vbs2_Ros3





  SUBROUTINE  saprc99_mosaic_4bin_vbs2_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
           ros_Gamma,ros_NewF,ros_ELO,ros_Name)










  IMPLICIT NONE

   INTEGER, PARAMETER :: S=4
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(4), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(6), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(4), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name

   REAL(kind=dp) :: g



   ros_Name = 'ROS-4'

   ros_S = S








   ros_A(1) = 0.2000000000000000E+01_dp
   ros_A(2) = 0.1867943637803922E+01_dp
   ros_A(3) = 0.2344449711399156E+00_dp
   ros_A(4) = ros_A(2)
   ros_A(5) = ros_A(3)
   ros_A(6) = 0.0_dp

   ros_C(1) =-0.7137615036412310E+01_dp
   ros_C(2) = 0.2580708087951457E+01_dp
   ros_C(3) = 0.6515950076447975E+00_dp
   ros_C(4) =-0.2137148994382534E+01_dp
   ros_C(5) =-0.3214669691237626E+00_dp
   ros_C(6) =-0.6949742501781779E+00_dp


   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .TRUE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .FALSE.

   ros_M(1) = 0.2255570073418735E+01_dp
   ros_M(2) = 0.2870493262186792E+00_dp
   ros_M(3) = 0.4353179431840180E+00_dp
   ros_M(4) = 0.1093502252409163E+01_dp

   ros_E(1) =-0.2815431932141155E+00_dp
   ros_E(2) =-0.7276199124938920E-01_dp
   ros_E(3) =-0.1082196201495311E+00_dp
   ros_E(4) =-0.1093502252409163E+01_dp


   ros_ELO  = 4.0_dp

   ros_Alpha(1) = 0.0_dp
   ros_Alpha(2) = 0.1145640000000000E+01_dp
   ros_Alpha(3) = 0.6552168638155900E+00_dp
   ros_Alpha(4) = ros_Alpha(3)

   ros_Gamma(1) = 0.5728200000000000E+00_dp
   ros_Gamma(2) =-0.1769193891319233E+01_dp
   ros_Gamma(3) = 0.7592633437920482E+00_dp
   ros_Gamma(4) =-0.1049021087100450E+00_dp

  END SUBROUTINE  saprc99_mosaic_4bin_vbs2_Ros4


  SUBROUTINE  saprc99_mosaic_4bin_vbs2_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
            ros_Gamma,ros_NewF,ros_ELO,ros_Name)




  IMPLICIT NONE

   INTEGER, PARAMETER :: S=4
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name

   REAL(kind=dp) :: g


   ros_Name = 'RODAS-3'

   ros_S = S








   ros_A(1) = 0.0E+00_dp
   ros_A(2) = 2.0E+00_dp
   ros_A(3) = 0.0E+00_dp
   ros_A(4) = 2.0E+00_dp
   ros_A(5) = 0.0E+00_dp
   ros_A(6) = 1.0E+00_dp

   ros_C(1) = 4.0E+00_dp
   ros_C(2) = 1.0E+00_dp
   ros_C(3) =-1.0E+00_dp
   ros_C(4) = 1.0E+00_dp
   ros_C(5) =-1.0E+00_dp
   ros_C(6) =-(8.0E+00_dp/3.0E+00_dp)



   ros_NewF(1)  = .TRUE.
   ros_NewF(2)  = .FALSE.
   ros_NewF(3)  = .TRUE.
   ros_NewF(4)  = .TRUE.

   ros_M(1) = 2.0E+00_dp
   ros_M(2) = 0.0E+00_dp
   ros_M(3) = 1.0E+00_dp
   ros_M(4) = 1.0E+00_dp

   ros_E(1) = 0.0E+00_dp
   ros_E(2) = 0.0E+00_dp
   ros_E(3) = 0.0E+00_dp
   ros_E(4) = 1.0E+00_dp


   ros_ELO  = 3.0E+00_dp

   ros_Alpha(1) = 0.0E+00_dp
   ros_Alpha(2) = 0.0E+00_dp
   ros_Alpha(3) = 1.0E+00_dp
   ros_Alpha(4) = 1.0E+00_dp

   ros_Gamma(1) = 0.5E+00_dp
   ros_Gamma(2) = 1.5E+00_dp
   ros_Gamma(3) = 0.0E+00_dp
   ros_Gamma(4) = 0.0E+00_dp

  END SUBROUTINE  saprc99_mosaic_4bin_vbs2_Rodas3


  SUBROUTINE  saprc99_mosaic_4bin_vbs2_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
             ros_Gamma,ros_NewF,ros_ELO,ros_Name)









  IMPLICIT NONE

   INTEGER, PARAMETER :: S=6
   INTEGER, INTENT(OUT) ::  ros_S
   REAL(kind=dp), DIMENSION(S), INTENT(OUT) :: ros_M,ros_E,ros_Alpha,ros_Gamma
   REAL(kind=dp), DIMENSION(S*(S-1)/2), INTENT(OUT) :: ros_A, ros_C
   REAL(kind=dp), INTENT(OUT) :: ros_ELO
   LOGICAL, DIMENSION(S), INTENT(OUT) :: ros_NewF
   CHARACTER(LEN=12), INTENT(OUT) :: ros_Name

    REAL(kind=dp) :: g


    ros_Name = 'RODAS-4'

    ros_S = 6


    ros_Alpha(1) = 0.000_dp
    ros_Alpha(2) = 0.386_dp
    ros_Alpha(3) = 0.210_dp
    ros_Alpha(4) = 0.630_dp
    ros_Alpha(5) = 1.000_dp
    ros_Alpha(6) = 1.000_dp


    ros_Gamma(1) = 0.2500000000000000E+00_dp
    ros_Gamma(2) =-0.1043000000000000E+00_dp
    ros_Gamma(3) = 0.1035000000000000E+00_dp
    ros_Gamma(4) =-0.3620000000000023E-01_dp
    ros_Gamma(5) = 0.0_dp
    ros_Gamma(6) = 0.0_dp







    ros_A(1) = 0.1544000000000000E+01_dp
    ros_A(2) = 0.9466785280815826E+00_dp
    ros_A(3) = 0.2557011698983284E+00_dp
    ros_A(4) = 0.3314825187068521E+01_dp
    ros_A(5) = 0.2896124015972201E+01_dp
    ros_A(6) = 0.9986419139977817E+00_dp
    ros_A(7) = 0.1221224509226641E+01_dp
    ros_A(8) = 0.6019134481288629E+01_dp
    ros_A(9) = 0.1253708332932087E+02_dp
    ros_A(10) =-0.6878860361058950E+00_dp
    ros_A(11) = ros_A(7)
    ros_A(12) = ros_A(8)
    ros_A(13) = ros_A(9)
    ros_A(14) = ros_A(10)
    ros_A(15) = 1.0E+00_dp

    ros_C(1) =-0.5668800000000000E+01_dp
    ros_C(2) =-0.2430093356833875E+01_dp
    ros_C(3) =-0.2063599157091915E+00_dp
    ros_C(4) =-0.1073529058151375E+00_dp
    ros_C(5) =-0.9594562251023355E+01_dp
    ros_C(6) =-0.2047028614809616E+02_dp
    ros_C(7) = 0.7496443313967647E+01_dp
    ros_C(8) =-0.1024680431464352E+02_dp
    ros_C(9) =-0.3399990352819905E+02_dp
    ros_C(10) = 0.1170890893206160E+02_dp
    ros_C(11) = 0.8083246795921522E+01_dp
    ros_C(12) =-0.7981132988064893E+01_dp
    ros_C(13) =-0.3152159432874371E+02_dp
    ros_C(14) = 0.1631930543123136E+02_dp
    ros_C(15) =-0.6058818238834054E+01_dp


    ros_M(1) = ros_A(7)
    ros_M(2) = ros_A(8)
    ros_M(3) = ros_A(9)
    ros_M(4) = ros_A(10)
    ros_M(5) = 1.0E+00_dp
    ros_M(6) = 1.0E+00_dp


    ros_E(1) = 0.0E+00_dp
    ros_E(2) = 0.0E+00_dp
    ros_E(3) = 0.0E+00_dp
    ros_E(4) = 0.0E+00_dp
    ros_E(5) = 0.0E+00_dp
    ros_E(6) = 1.0E+00_dp



    ros_NewF(1) = .TRUE.
    ros_NewF(2) = .TRUE.
    ros_NewF(3) = .TRUE.
    ros_NewF(4) = .TRUE.
    ros_NewF(5) = .TRUE.
    ros_NewF(6) = .TRUE.



    ros_ELO = 4.0_dp

  END SUBROUTINE  saprc99_mosaic_4bin_vbs2_Rodas4




END SUBROUTINE  saprc99_mosaic_4bin_vbs2_Rosenbrock




SUBROUTINE  saprc99_mosaic_4bin_vbs2_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE saprc99_mosaic_4bin_vbs2_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL saprc99_mosaic_4bin_vbs2_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  saprc99_mosaic_4bin_vbs2_FunTemplate



SUBROUTINE  saprc99_mosaic_4bin_vbs2_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE saprc99_mosaic_4bin_vbs2_Parameters
 
 USE saprc99_mosaic_4bin_vbs2_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL saprc99_mosaic_4bin_vbs2_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  saprc99_mosaic_4bin_vbs2_JacTemplate

















SUBROUTINE saprc99_mosaic_4bin_vbs2_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(92)
  A(2) = RCT(2)*V(84)*F(2)
  A(3) = RCT(3)*V(84)*V(88)
  A(4) = RCT(4)*V(84)*V(100)*F(2)
  A(5) = RCT(5)*V(84)*V(92)
  A(6) = RCT(6)*V(84)*V(92)
  A(7) = RCT(7)*V(88)*V(100)
  A(8) = RCT(8)*V(88)*V(92)
  A(9) = RCT(9)*V(93)*V(100)
  A(10) = RCT(10)*V(100)*V(100)*F(2)
  A(11) = RCT(11)*V(92)*V(93)
  A(12) = RCT(12)*V(42)
  A(13) = RCT(13)*V(42)*F(1)
  A(14) = RCT(14)*V(92)*V(93)
  A(15) = RCT(15)*V(93)
  A(16) = RCT(16)*V(93)
  A(17) = RCT(17)*V(88)
  A(18) = RCT(18)*V(88)
  A(19) = RCT(19)*V(26)*F(1)
  A(20) = RCT(20)*V(26)*F(2)
  A(21) = RCT(21)*V(98)*V(100)
  A(22) = RCT(22)*V(43)
  A(23) = RCT(23)*V(43)
  A(24) = RCT(24)*V(43)*V(98)
  A(25) = RCT(25)*V(92)*V(98)
  A(26) = RCT(26)*V(93)*V(98)
  A(27) = RCT(27)*V(66)*V(98)
  A(28) = RCT(28)*V(66)
  A(29) = RCT(29)*V(65)*V(98)
  A(30) = RCT(30)*V(88)*V(98)
  A(31) = RCT(31)*V(97)*V(100)
  A(32) = RCT(32)*V(92)*V(97)
  A(33) = RCT(33)*V(51)
  A(34) = RCT(34)*V(51)
  A(35) = RCT(35)*V(51)*V(98)
  A(36) = RCT(36)*V(88)*V(97)
  A(37) = RCT(37)*V(97)*V(97)
  A(38) = RCT(38)*V(97)*V(97)*F(1)
  A(39) = RCT(39)*V(93)*V(97)
  A(40) = RCT(40)*V(93)*V(93)
  A(41) = RCT(41)*V(38)
  A(42) = RCT(42)*V(38)*V(98)
  A(43) = RCT(43)*V(97)*V(98)
  A(44) = RCT(44)*V(32)*V(98)
  A(45) = RCT(45)*V(98)*F(2)
  A(46) = RCT(46)*V(94)*V(100)
  A(47) = RCT(47)*V(94)*V(97)
  A(48) = RCT(48)*V(93)*V(94)
  A(49) = RCT(49)*V(94)*V(94)
  A(50) = RCT(50)*V(94)*V(94)
  A(51) = RCT(51)*V(90)*V(100)
  A(52) = RCT(52)*V(90)*V(97)
  A(53) = RCT(53)*V(90)*V(93)
  A(54) = RCT(54)*V(90)*V(94)
  A(55) = RCT(55)*V(90)*V(90)
  A(56) = RCT(56)*V(73)*V(100)
  A(57) = RCT(57)*V(73)*V(97)
  A(58) = RCT(58)*V(73)*V(93)
  A(59) = RCT(59)*V(73)*V(94)
  A(60) = RCT(60)*V(73)*V(90)
  A(61) = RCT(61)*V(73)*V(73)
  A(62) = RCT(62)*V(91)*V(100)
  A(63) = RCT(63)*V(91)*V(97)
  A(64) = RCT(64)*V(91)*V(94)
  A(65) = RCT(65)*V(91)*V(93)
  A(66) = RCT(66)*V(90)*V(91)
  A(67) = RCT(67)*V(73)*V(91)
  A(68) = RCT(68)*V(91)*V(91)
  A(69) = RCT(69)*V(89)*V(92)
  A(70) = RCT(70)*V(34)
  A(71) = RCT(71)*V(89)*V(100)
  A(72) = RCT(72)*V(89)*V(97)
  A(73) = RCT(73)*V(89)*V(93)
  A(74) = RCT(74)*V(89)*V(94)
  A(75) = RCT(75)*V(89)*V(90)
  A(76) = RCT(76)*V(73)*V(89)
  A(77) = RCT(77)*V(89)*V(91)
  A(78) = RCT(78)*V(89)*V(89)
  A(79) = RCT(79)*V(92)*V(99)
  A(80) = RCT(80)*V(35)
  A(81) = RCT(81)*V(99)*V(100)
  A(82) = RCT(82)*V(97)*V(99)
  A(83) = RCT(83)*V(93)*V(99)
  A(84) = RCT(84)*V(94)*V(99)
  A(85) = RCT(85)*V(90)*V(99)
  A(86) = RCT(86)*V(73)*V(99)
  A(87) = RCT(87)*V(91)*V(99)
  A(88) = RCT(88)*V(89)*V(99)
  A(89) = RCT(89)*V(99)*V(99)
  A(90) = RCT(90)*V(92)*V(95)
  A(91) = RCT(91)*V(36)
  A(92) = RCT(92)*V(95)*V(100)
  A(93) = RCT(93)*V(95)*V(97)
  A(94) = RCT(94)*V(93)*V(95)
  A(95) = RCT(95)*V(94)*V(95)
  A(96) = RCT(96)*V(90)*V(95)
  A(97) = RCT(97)*V(73)*V(95)
  A(98) = RCT(98)*V(91)*V(95)
  A(99) = RCT(99)*V(89)*V(95)
  A(100) = RCT(100)*V(95)*V(99)
  A(101) = RCT(101)*V(95)*V(95)
  A(102) = RCT(102)*V(92)*V(96)
  A(103) = RCT(103)*V(37)
  A(104) = RCT(104)*V(96)*V(100)
  A(105) = RCT(105)*V(96)*V(97)
  A(106) = RCT(106)*V(93)*V(96)
  A(107) = RCT(107)*V(94)*V(96)
  A(108) = RCT(108)*V(90)*V(96)
  A(109) = RCT(109)*V(73)*V(96)
  A(110) = RCT(110)*V(91)*V(96)
  A(111) = RCT(111)*V(89)*V(96)
  A(112) = RCT(112)*V(96)*V(99)
  A(113) = RCT(113)*V(95)*V(96)
  A(114) = RCT(114)*V(96)*V(96)
  A(115) = RCT(115)*V(45)*V(92)
  A(116) = RCT(116)*V(45)
  A(117) = RCT(117)*V(70)*V(92)
  A(118) = RCT(118)*V(70)*V(97)
  A(119) = RCT(119)*V(70)
  A(120) = RCT(120)*V(49)*V(92)
  A(121) = RCT(121)*V(49)*V(97)
  A(122) = RCT(122)*V(49)
  A(123) = RCT(123)*V(82)
  A(124) = RCT(124)*V(82)
  A(125) = RCT(125)*V(82)*V(98)
  A(126) = RCT(126)*V(82)*V(97)
  A(127) = RCT(127)*V(48)
  A(128) = RCT(128)*V(48)*V(100)
  A(129) = RCT(129)*V(82)*V(93)
  A(130) = RCT(130)*V(79)*V(98)
  A(131) = RCT(131)*V(79)
  A(132) = RCT(132)*V(79)*V(93)
  A(133) = RCT(133)*V(85)*V(98)
  A(134) = RCT(134)*V(85)
  A(135) = RCT(135)*V(85)*V(93)
  A(136) = RCT(136)*V(68)*V(98)
  A(137) = RCT(137)*V(68)
  A(138) = RCT(138)*V(86)*V(98)
  A(139) = RCT(139)*V(86)
  A(140) = RCT(140)*V(52)*V(98)
  A(141) = RCT(141)*V(41)*V(98)
  A(142) = RCT(142)*V(47)*V(98)
  A(143) = RCT(143)*V(47)
  A(144) = RCT(144)*V(60)*V(98)
  A(145) = RCT(145)*V(60)
  A(146) = RCT(146)*V(77)
  A(147) = RCT(147)*V(77)
  A(148) = RCT(148)*V(77)*V(98)
  A(149) = RCT(149)*V(77)*V(93)
  A(150) = RCT(150)*V(64)
  A(151) = RCT(151)*V(64)*V(98)
  A(152) = RCT(152)*V(64)*V(93)
  A(153) = RCT(153)*V(40)
  A(154) = RCT(154)*V(63)*V(98)
  A(155) = RCT(155)*V(63)*V(93)
  A(156) = RCT(156)*V(57)*V(98)
  A(157) = RCT(157)*V(57)*V(93)
  A(158) = RCT(158)*V(61)*V(93)
  A(159) = RCT(159)*V(62)*V(98)
  A(160) = RCT(160)*V(62)
  A(161) = RCT(161)*V(62)*V(93)
  A(162) = RCT(162)*V(74)*V(98)
  A(163) = RCT(163)*V(74)*V(88)
  A(164) = RCT(164)*V(74)*V(93)
  A(165) = RCT(165)*V(74)*V(84)
  A(166) = RCT(166)*V(74)
  A(167) = RCT(167)*V(81)*V(98)
  A(168) = RCT(168)*V(81)*V(88)
  A(169) = RCT(169)*V(81)*V(84)
  A(170) = RCT(170)*V(81)
  A(171) = RCT(171)*V(78)*V(98)
  A(172) = RCT(172)*V(78)*V(88)
  A(173) = RCT(173)*V(78)*V(93)
  A(174) = RCT(174)*V(78)
  A(175) = RCT(175)*V(87)*V(98)
  A(176) = RCT(176)*V(87)
  A(177) = RCT(177)*V(83)*V(98)
  A(178) = RCT(178)*V(83)
  A(179) = RCT(179)*V(58)*V(98)
  A(180) = RCT(180)*V(58)*V(88)
  A(181) = RCT(181)*V(55)*V(98)
  A(182) = RCT(182)*V(55)
  A(183) = RCT(183)*V(56)*V(98)
  A(184) = RCT(184)*V(56)
  A(185) = RCT(185)*V(31)*V(98)
  A(186) = RCT(186)*V(67)*V(98)
  A(187) = RCT(187)*V(67)*V(88)
  A(188) = RCT(188)*V(67)*V(93)
  A(189) = RCT(189)*V(67)*V(84)
  A(190) = RCT(190)*V(71)*V(98)
  A(191) = RCT(191)*V(71)*V(88)
  A(192) = RCT(192)*V(71)*V(93)
  A(193) = RCT(193)*V(71)*V(84)
  A(194) = RCT(194)*V(76)*V(98)
  A(195) = RCT(195)*V(76)*V(88)
  A(196) = RCT(196)*V(76)*V(93)
  A(197) = RCT(197)*V(76)*V(84)
  A(198) = RCT(198)*V(75)*V(98)
  A(199) = RCT(199)*V(75)*V(88)
  A(200) = RCT(200)*V(75)*V(93)
  A(201) = RCT(201)*V(75)*V(84)
  A(202) = RCT(202)*V(33)*V(98)
  A(203) = RCT(203)*V(39)*V(98)
  A(204) = RCT(204)*V(59)*V(98)
  A(205) = RCT(205)*V(44)*V(98)
  A(206) = RCT(206)*V(53)*V(98)
  A(207) = RCT(207)*V(46)*V(98)
  A(208) = RCT(208)*V(54)*V(98)
  A(209) = RCT(209)*V(50)*V(98)
  A(210) = RCT(210)*V(72)*V(98)
  A(211) = RCT(211)*V(72)*V(88)
  A(212) = RCT(212)*V(72)*V(93)
  A(213) = RCT(213)*V(72)*V(84)
  A(214) = RCT(214)*V(80)*V(98)
  A(215) = RCT(215)*V(80)*V(88)
  A(216) = RCT(216)*V(80)*V(93)
  A(217) = RCT(217)*V(80)*V(84)
  A(218) = RCT(218)*V(59)*V(88)
  A(219) = RCT(219)*V(69)*V(98)
  A(220) = RCT(220)*V(69)*V(88)
  A(221) = RCT(221)*V(69)*V(93)
  A(222) = RCT(222)*V(69)*V(84)
  A(223) = RCT(223)*V(32)
  A(224) = RCT(224)*V(97)
  A(225) = RCT(225)*V(32)
  A(226) = RCT(226)*V(1)
  A(227) = RCT(227)*V(66)
  A(228) = RCT(228)*V(38)
  A(229) = RCT(229)*V(2)
  A(230) = RCT(230)*V(53)*V(98)
  A(231) = RCT(231)*V(46)*V(98)
  A(232) = RCT(232)*V(72)*V(98)
  A(233) = RCT(233)*V(80)*V(98)
  A(234) = RCT(234)*V(54)*V(98)
  A(235) = RCT(235)*V(50)*V(98)
  A(236) = RCT(236)*V(71)*V(98)
  A(237) = RCT(237)*V(76)*V(98)
  A(238) = RCT(238)*V(75)*V(98)
  A(239) = RCT(239)*V(3)*V(98)
  A(240) = RCT(240)*V(4)*V(98)
  A(241) = RCT(241)*V(28)*V(98)
  A(242) = RCT(242)*V(22)*V(98)
  A(243) = RCT(243)*V(5)*V(98)
  A(244) = RCT(244)*V(6)*V(98)
  A(245) = RCT(245)*V(30)*V(98)
  A(246) = RCT(246)*V(24)*V(98)
  A(247) = RCT(247)*V(27)*V(98)
  A(248) = RCT(248)*V(23)*V(98)
  A(249) = RCT(249)*V(29)*V(98)
  A(250) = RCT(250)*V(25)*V(98)


  Vdot(1) = A(44)+A(223)-A(226)
  Vdot(2) = 0.5*A(218)+0.135*A(220)-A(229)
  Vdot(3) = 0
  Vdot(4) = 0
  Vdot(5) = 0
  Vdot(6) = 0
  Vdot(7) = A(128)+0.333*A(163)+0.351*A(168)+0.1*A(172)+0.37*A(187)+0.204*A(191)+0.103*A(195)+0.103*A(199)+0.297*A(204)&
              &+0.185*A(211)+0.073*A(215)+0.185*A(220)
  Vdot(8) = 0.25*A(72)+A(74)+A(75)+A(77)+0.05*A(211)+0.129*A(215)+0.17*A(220)
  Vdot(9) = 0.25*A(82)+A(84)+A(85)+A(87)+0.25*A(93)+A(95)+A(96)+A(98)+0.25*A(105)+A(107)+A(108)+2*A(110)+0.372*A(172)&
              &+0.15*A(191)+0.189*A(195)+0.189*A(199)+0.119*A(211)+0.247*A(215)
  Vdot(10) = 0.75*A(72)
  Vdot(11) = 0.75*A(82)+0.75*A(93)+0.75*A(105)
  Vdot(12) = 2*A(120)+A(221)
  Vdot(13) = 6*A(120)+7*A(160)+0.048*A(219)+0.07*A(220)+2.693*A(221)+0.55*A(222)
  Vdot(14) = A(46)+A(48)+A(51)+A(53)+A(56)+A(58)+A(62)+A(65)
  Vdot(15) = A(47)+A(49)+A(50)+A(52)+A(54)+A(55)+A(57)+A(59)+A(60)+A(61)+A(63)+A(64)+A(66)+A(67)+A(68)
  Vdot(16) = A(230)+A(231)+A(232)+A(233)
  Vdot(17) = A(234)+A(235)
  Vdot(18) = A(236)
  Vdot(19) = A(237)+A(238)
  Vdot(20) = A(21)+A(24)+A(25)+A(26)+A(27)+A(29)+A(30)+A(43)+A(44)+A(45)+A(125)+A(130)+A(133)+A(136)+A(138)+A(140)&
               &+A(141)+A(142)+A(144)+A(148)+A(151)+A(154)+A(156)+A(159)+A(162)+A(167)+A(171)+A(175)+A(177)+A(179)+A(181)&
               &+A(183)+A(185)+A(186)+A(190)+A(194)+A(198)+A(202)+A(203)+A(204)+A(205)+A(206)+A(207)+A(208)+A(209)+A(210)&
               &+A(214)+A(219)
  Vdot(21) = A(239)+A(240)+A(241)+A(242)+A(243)+A(244)+A(245)+A(246)+A(247)+A(248)+A(249)+A(250)
  Vdot(22) = 0.5*A(247)+A(248)
  Vdot(23) = -A(248)
  Vdot(24) = 0.5*A(249)+A(250)
  Vdot(25) = -A(250)
  Vdot(26) = A(18)-A(19)-A(20)
  Vdot(27) = -A(247)
  Vdot(28) = A(247)
  Vdot(29) = -A(249)
  Vdot(30) = A(249)
  Vdot(31) = -A(185)
  Vdot(32) = -A(44)-A(223)-A(225)
  Vdot(33) = -A(202)
  Vdot(34) = A(69)-A(70)
  Vdot(35) = A(79)-A(80)
  Vdot(36) = A(90)-A(91)
  Vdot(37) = A(102)-A(103)
  Vdot(38) = A(37)+A(38)-A(41)-A(42)-A(228)
  Vdot(39) = -A(203)
  Vdot(40) = -A(153)+0.031*A(195)+0.031*A(199)+0.087*A(209)
  Vdot(41) = -A(141)
  Vdot(42) = A(11)-A(12)-A(13)
  Vdot(43) = A(21)-A(22)-A(23)-A(24)
  Vdot(44) = -A(205)
  Vdot(45) = -A(115)-A(116)+0.236*A(205)
  Vdot(46) = -A(207)
  Vdot(47) = A(47)-A(142)-A(143)
  Vdot(48) = A(126)-A(127)-A(128)
  Vdot(49) = -A(120)-A(121)-A(122)+A(158)
  Vdot(50) = -A(209)
  Vdot(51) = A(32)-A(33)-A(34)-A(35)
  Vdot(52) = A(49)+0.25*A(54)+0.25*A(64)-A(140)
  Vdot(53) = -A(206)
  Vdot(54) = -A(208)
  Vdot(55) = -A(181)-A(182)+0.108*A(208)+0.099*A(209)
  Vdot(56) = -A(183)-A(184)+0.051*A(208)+0.093*A(209)
  Vdot(57) = -A(156)-A(157)+0.207*A(208)+0.187*A(209)
  Vdot(58) = -A(179)-A(180)+0.491*A(208)+0.561*A(209)
  Vdot(59) = -A(204)-A(218)
  Vdot(60) = A(52)+A(63)-A(144)-A(145)
  Vdot(61) = A(117)+A(121)+A(122)-A(158)
  Vdot(62) = -A(159)-A(160)-A(161)+0.059*A(208)+0.05*A(209)+0.061*A(214)+0.042*A(215)+0.015*A(216)
  Vdot(63) = A(118)+A(119)-A(154)-A(155)+0.017*A(208)
  Vdot(64) = -A(150)-A(151)-A(152)+0.23*A(156)+0.084*A(162)+0.9*A(163)+0.3*A(167)+0.95*A(168)+0.174*A(171)+0.742*A(172)&
               &+0.008*A(173)+0.5*A(182)+0.5*A(184)+0.119*A(208)+0.287*A(209)
  Vdot(65) = -A(29)+A(123)+A(124)+A(125)+A(129)+A(131)+0.034*A(133)+A(134)+2*A(146)+A(147)+1.26*A(148)+1.26*A(149)&
               &+A(150)+A(151)+A(152)+0.416*A(162)+0.45*A(163)+0.5*A(164)+0.67*A(166)+0.475*A(168)+0.7*A(170)+0.336*A(171)&
               &+0.498*A(172)+0.572*A(173)+1.233*A(174)+A(179)+1.5*A(180)+A(182)+A(184)+0.5*A(187)+0.491*A(189)+0.275*A(191)&
               &+0.157*A(195)+0.157*A(199)+0.393*A(204)+0.002*A(206)+0.345*A(211)+0.265*A(215)+0.012*A(217)+1.5*A(218)+0.51&
               &*A(220)
  Vdot(66) = 2*A(13)+A(25)-A(27)-A(28)+0.2*A(39)+A(129)+A(132)+A(135)+A(149)+A(152)+A(155)+A(157)+A(158)+A(161)+0.5&
               &*A(164)+0.15*A(173)-A(227)
  Vdot(67) = -A(186)-A(187)-A(188)-A(189)
  Vdot(68) = A(116)-A(136)-A(137)+0.006*A(177)+0.02*A(178)+0.13*A(195)+0.13*A(199)+0.704*A(203)+0.024*A(205)+0.452&
               &*A(206)+0.072*A(207)+0.005*A(210)+0.001*A(211)+0.024*A(212)+0.127*A(214)+0.045*A(215)+0.102*A(216)
  Vdot(69) = -A(219)-A(220)-A(221)-A(222)
  Vdot(70) = A(92)+A(94)+A(99)+A(100)+2*A(101)+A(113)-A(117)-A(118)-A(119)+0.24*A(154)+A(155)+0.24*A(156)+A(157)
  Vdot(71) = -A(190)-A(191)-A(192)-A(193)
  Vdot(72) = -A(210)-A(211)-A(212)-A(213)
  Vdot(73) = -A(56)-A(57)-A(58)-A(59)-A(60)-A(67)-A(76)-A(86)+A(92)+A(94)-A(97)+A(99)+A(100)+2*A(101)-A(109)+A(113)&
               &+A(136)+0.616*A(138)+0.675*A(167)+0.515*A(176)+0.596*A(177)+0.152*A(178)+A(181)+A(182)+A(183)+A(184)+0.079&
               &*A(190)+0.126*A(191)+0.187*A(192)+0.24*A(193)+0.5*A(194)+0.729*A(195)+0.75*A(196)+0.5*A(198)+0.729*A(199)&
               &+0.75*A(200)+0.559*A(205)+0.936*A(206)+0.948*A(207)+0.205*A(210)+0.488*A(212)+0.001*A(214)+0.137*A(215)&
               &+0.711*A(216)
  Vdot(74) = -A(162)-A(163)-A(164)-A(165)-A(166)+0.23*A(190)+0.39*A(191)+0.025*A(214)+0.026*A(215)+0.012*A(217)
  Vdot(75) = -A(198)-A(199)-A(200)-A(201)
  Vdot(76) = -A(194)-A(195)-A(196)-A(197)
  Vdot(77) = -A(146)-A(147)-A(148)-A(149)+0.23*A(154)+0.15*A(171)+0.023*A(172)+A(180)+0.5*A(182)+0.5*A(184)+0.009*A(189)&
               &+0.001*A(195)+0.001*A(199)+0.607*A(204)+0.118*A(208)+0.097*A(209)
  Vdot(78) = -A(171)-A(172)-A(173)-A(174)+0.357*A(190)+0.936*A(192)+0.025*A(214)
  Vdot(79) = A(81)+A(83)+A(88)+2*A(89)+A(100)+A(112)-A(130)-A(131)-A(132)+0.034*A(133)+A(134)+0.482*A(138)+A(139)+0.96&
               &*A(141)+0.129*A(171)+0.047*A(172)+0.467*A(174)+0.084*A(175)+0.246*A(176)+0.439*A(177)+0.431*A(178)+0.195&
               &*A(186)+0.25*A(189)+A(202)+0.445*A(205)+0.455*A(206)+0.099*A(207)+0.294*A(210)+0.154*A(211)+0.009*A(212)&
               &+0.732*A(214)+0.456*A(215)+0.507*A(216)+0.984*A(219)+0.5*A(220)
  Vdot(80) = -A(214)-A(215)-A(216)-A(217)
  Vdot(81) = -A(167)-A(168)-A(169)-A(170)+0.32*A(190)+0.16*A(191)+0.019*A(215)+0.048*A(216)
  Vdot(82) = A(46)+A(48)+A(49)+2*A(50)+0.75*A(54)+0.75*A(64)+A(74)+A(84)+A(95)+A(104)+A(106)+A(107)+A(111)+A(112)+A(113)&
               &+2*A(114)-A(123)-A(124)-A(125)-A(126)+A(127)-A(129)+A(136)+0.115*A(138)+A(140)+0.081*A(141)+0.35*A(142)&
               &+A(143)+A(147)+0.084*A(162)+0.2*A(163)+0.67*A(166)+0.3*A(167)+0.1*A(168)+0.055*A(171)+0.125*A(172)+0.227&
               &*A(173)+0.3*A(174)+0.213*A(175)+0.506*A(176)+0.01*A(177)+0.134*A(178)+1.61*A(186)+A(187)+0.191*A(189)+0.624&
               &*A(190)+0.592*A(191)+0.24*A(193)+0.276*A(194)+0.235*A(195)+0.276*A(198)+0.235*A(199)+0.096*A(204)+0.026&
               &*A(205)+0.024*A(206)+0.026*A(207)+0.732*A(210)+0.5*A(211)+0.244*A(214)+0.269*A(215)+0.079*A(216)+0.984&
               &*A(219)+0.5*A(220)
  Vdot(83) = A(62)+A(115)+0.572*A(173)-0.69*A(177)-A(178)+0.276*A(196)+0.276*A(200)+0.511*A(212)+0.321*A(216)
  Vdot(84) = A(1)-A(2)-A(3)-A(4)-A(5)-A(6)+A(16)+A(17)+A(20)-A(165)-A(169)-A(189)-A(193)-A(197)-A(201)-A(213)-A(217)&
               &-A(222)
  Vdot(85) = -A(133)-A(134)-A(135)+0.37*A(138)+A(144)+A(145)+A(165)+0.675*A(167)+0.45*A(169)+0.013*A(171)+0.218*A(173)&
               &+0.558*A(175)+0.71*A(176)+0.213*A(177)+0.147*A(178)+A(179)+A(181)+A(183)+A(188)+0.474*A(194)+0.205*A(195)&
               &+0.474*A(196)+0.147*A(197)+0.474*A(198)+0.205*A(199)+0.474*A(200)+0.147*A(201)+0.261*A(203)+0.122*A(205)&
               &+0.244*A(206)+0.204*A(207)+0.497*A(210)+0.363*A(211)+0.037*A(212)+0.45*A(213)+0.511*A(214)+0.305*A(215)&
               &+0.151*A(216)+0.069*A(217)+0.45*A(222)
  Vdot(86) = 0.5*A(64)+A(65)+0.5*A(66)+A(68)-A(138)-A(139)+0.416*A(162)+0.55*A(169)+0.15*A(171)+0.21*A(172)+0.233*A(174)&
               &+0.115*A(175)+0.177*A(177)+0.243*A(178)+0.332*A(205)+0.11*A(206)+0.089*A(207)+0.437*A(213)+0.072*A(214)&
               &+0.026*A(215)+0.001*A(216)+0.659*A(217)+0.55*A(222)
  Vdot(87) = 0.5*A(64)+0.5*A(66)+A(68)+A(77)+A(87)+A(98)+0.7*A(170)+0.332*A(171)-0.671*A(175)-A(176)+0.048*A(177)+0.435&
               &*A(178)+0.1*A(191)+0.75*A(193)+0.276*A(194)+0.276*A(195)+0.853*A(197)+0.276*A(198)+0.276*A(199)+0.853*A(201)&
               &+0.125*A(206)+0.417*A(207)+0.055*A(208)+0.119*A(210)+0.215*A(211)+0.113*A(213)+0.043*A(215)+0.259*A(217)
  Vdot(88) = A(2)-A(3)-A(7)-A(8)-A(17)-A(18)-A(30)-A(36)+0.25*A(72)+0.25*A(82)+0.25*A(93)+0.25*A(105)-A(163)-A(168)&
               &-A(172)-A(180)-A(187)-A(191)-A(195)-A(199)-A(211)-A(215)-A(218)-A(220)
  Vdot(89) = -A(69)+A(70)-A(71)-A(72)-A(73)-A(74)-A(75)-A(77)-2*A(78)-A(88)-A(99)+A(104)+A(106)+A(112)+A(113)+2*A(114)&
               &+A(130)+A(132)+A(136)+A(137)+0.492*A(138)+A(139)+A(150)+A(151)+A(152)+2*A(153)+0.67*A(166)+0.675*A(167)&
               &+0.467*A(174)+0.029*A(175)+0.667*A(176)+A(181)+0.5*A(182)+A(183)+0.5*A(184)+0.123*A(195)+0.123*A(199)+0.011&
               &*A(206)+0.137*A(215)
  Vdot(90) = -A(51)-A(52)-A(53)-A(54)-2*A(55)-A(66)-A(75)+A(81)+A(83)-A(85)+A(88)+2*A(89)-A(96)+A(100)-A(108)+A(112)&
               &+0.034*A(133)+A(134)+0.37*A(138)+A(139)+0.05*A(141)+0.34*A(144)+0.76*A(154)+0.76*A(156)+0.5*A(162)+0.1&
               &*A(163)+0.5*A(164)+0.33*A(166)+0.3*A(167)+0.05*A(168)+0.67*A(171)+0.048*A(172)+0.799*A(173)+0.473*A(175)&
               &+0.96*A(176)+0.376*A(177)+0.564*A(178)+A(179)+A(182)+A(184)+A(186)+A(188)+0.2*A(189)+0.907*A(190)+0.066&
               &*A(191)+0.749*A(192)+0.75*A(194)+0.031*A(195)+0.276*A(196)+0.75*A(198)+0.031*A(199)+0.276*A(200)+A(202)&
               &+0.965*A(203)+0.1*A(204)+0.695*A(205)+0.835*A(206)+0.653*A(207)+0.765*A(208)+0.804*A(209)+0.91*A(210)+0.022&
               &*A(211)+0.824*A(212)+0.918*A(214)+0.033*A(215)+0.442*A(216)+0.012*A(217)+0.984*A(219)+0.949*A(221)
  Vdot(91) = -A(62)-A(63)-A(64)-A(65)-A(66)-2*A(68)-A(77)-A(87)-A(98)-A(110)+0.001*A(133)+0.042*A(138)+0.025*A(167)&
               &+0.041*A(171)+0.051*A(173)+0.07*A(175)+0.04*A(176)+0.173*A(177)+0.095*A(178)+0.093*A(190)+0.008*A(191)+0.064&
               &*A(192)+0.01*A(193)+0.25*A(194)+0.18*A(195)+0.25*A(196)+0.25*A(198)+0.18*A(199)+0.25*A(200)+0.035*A(203)&
               &+0.07*A(205)+0.143*A(206)+0.347*A(207)+0.011*A(208)+0.009*A(209)+0.09*A(210)+0.001*A(211)+0.176*A(212)+0.082&
               &*A(214)+0.002*A(215)+0.136*A(216)+0.001*A(217)+0.016*A(219)+0.051*A(221)
  Vdot(92) = -A(1)+A(4)-A(5)-A(6)+A(7)-A(8)+2*A(9)+2*A(10)-A(11)+A(12)+A(16)+A(23)+A(24)-A(25)+A(26)+A(28)+A(31)-A(32)&
               &+A(33)+0.61*A(34)+A(35)+0.8*A(39)+2*A(40)+A(46)+A(48)+A(51)+A(53)+A(56)+A(58)+A(65)-A(69)+A(70)+A(71)+A(73)&
               &-A(79)+A(80)+A(81)+A(83)-A(90)+A(91)+A(92)+A(94)-A(102)+A(103)+A(104)+A(106)-A(115)-A(117)-A(120)+A(128)&
               &+0.338*A(177)+A(178)+0.187*A(192)+0.474*A(196)+0.474*A(200)+0.391*A(216)
  Vdot(93) = A(6)+A(8)-A(9)-A(11)+A(12)-A(14)-A(15)-A(16)-A(26)+A(27)+0.39*A(34)-A(39)-2*A(40)-A(48)-A(53)-A(58)-A(65)&
               &-A(73)-A(83)-A(94)-A(106)-A(129)-A(132)-A(135)-A(149)-A(152)-A(155)-A(157)-A(158)-A(161)-A(164)-A(173)&
               &-A(188)-A(192)-A(196)-A(200)-A(212)-A(216)-A(221)
  Vdot(94) = -A(46)-A(47)-A(48)-2*A(49)-2*A(50)-A(54)-A(64)+A(71)+A(73)-A(74)+2*A(78)-A(84)+A(88)-A(95)+A(99)-A(107)&
               &+A(111)+A(116)+A(131)+A(137)+0.65*A(142)+0.3*A(170)+A(185)+0.3*A(189)+0.25*A(193)+0.011*A(206)+0.076*A(211)&
               &+0.197*A(215)+0.03*A(216)+0.26*A(220)
  Vdot(95) = -A(90)+A(91)-A(92)-A(93)-A(94)-A(95)-A(96)-A(98)-A(99)-A(100)-2*A(101)-A(113)+A(159)+A(161)
  Vdot(96) = -A(102)+A(103)-A(104)-A(105)-A(106)-A(107)-A(108)-A(110)-A(111)-A(112)-A(113)-2*A(114)+0.5*A(162)+0.5&
               &*A(164)+0.33*A(166)+0.3*A(170)+0.289*A(171)+0.15*A(173)+0.192*A(191)+0.24*A(193)
  Vdot(97) = A(23)+A(26)+A(29)+A(30)-A(31)-A(32)+A(33)+0.61*A(34)-A(36)-2*A(37)-2*A(38)-A(39)+A(42)-A(43)+A(44)+A(45)&
               &+A(46)-A(47)+A(48)+2*A(50)+A(51)-A(52)+A(53)+A(54)+A(55)-A(63)+A(64)+A(65)+A(66)+A(68)-A(72)-A(82)-A(93)&
               &-A(105)-A(118)-A(121)+2*A(123)+A(125)-A(126)+A(127)+A(128)+A(129)+A(131)+A(134)+A(140)+0.95*A(141)+A(143)&
               &+A(145)+2*A(146)+0.63*A(148)+0.63*A(149)+A(150)+0.008*A(163)+0.34*A(166)+0.064*A(168)+0.4*A(172)+1.233&
               &*A(174)+0.379*A(175)+0.113*A(177)+0.341*A(178)+1.5*A(180)+0.5*A(182)+0.5*A(184)+0.12*A(187)+0.5*A(189)+0.033&
               &*A(195)+0.033*A(199)+0.297*A(204)+0.224*A(208)+0.187*A(209)+0.056*A(211)+0.003*A(215)+0.013*A(217)+1.5&
               &*A(218)+0.06*A(220)-A(224)
  Vdot(98) = 2*A(19)-A(21)+A(22)-A(24)-A(25)-A(26)-A(27)+A(28)-A(29)-A(30)+A(31)+0.39*A(34)-A(35)+A(36)+0.8*A(39)+2&
               &*A(41)-A(42)-A(43)-A(44)-A(45)-A(125)-A(130)-A(133)-A(136)-A(138)-A(140)-A(141)-0.65*A(142)+A(143)-0.34&
               &*A(144)+A(145)-A(148)-A(151)-A(154)-A(156)-A(159)-A(162)+0.208*A(163)+0.33*A(166)-A(167)+0.164*A(168)-A(171)&
               &+0.285*A(172)-A(175)-A(177)-A(179)+0.5*A(180)-A(181)-A(183)-A(185)-A(186)+0.12*A(187)-A(190)+0.266*A(191)&
               &-A(194)+0.567*A(195)-A(198)+0.567*A(199)-A(202)-A(203)-0.397*A(204)-A(205)-A(206)-A(207)-A(208)-A(209)&
               &-A(210)+0.155*A(211)-A(214)+0.378*A(215)+0.5*A(218)-A(219)+0.32*A(220)-A(239)-A(241)-A(243)-A(245)-A(247)&
               &-A(249)
  Vdot(99) = -A(79)+A(80)-A(81)-A(82)-A(83)-A(84)-A(85)-A(87)-A(88)-2*A(89)-A(100)-A(112)+0.965*A(133)+A(135)+0.096&
               &*A(138)+0.37*A(148)+0.37*A(149)+0.1*A(163)+0.05*A(168)+0.048*A(172)+0.3*A(174)+0.049*A(175)+0.333*A(176)&
               &+0.201*A(195)+0.201*A(199)+0.006*A(215)
  Vdot(100) = A(1)-A(4)+A(5)-A(7)-A(9)-2*A(10)+A(14)+A(15)-A(21)+A(22)-A(31)-A(46)-A(51)-A(56)-A(62)-A(71)-A(81)-A(92)&
                &-A(104)-A(128)
      
END SUBROUTINE saprc99_mosaic_4bin_vbs2_Fun
















SUBROUTINE saprc99_mosaic_4bin_vbs2_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(92)
  IRR(2) = RCT(2)*V(84)*F(2)
  IRR(3) = RCT(3)*V(84)*V(88)
  IRR(4) = RCT(4)*V(84)*V(100)*F(2)
  IRR(5) = RCT(5)*V(84)*V(92)
  IRR(6) = RCT(6)*V(84)*V(92)
  IRR(7) = RCT(7)*V(88)*V(100)
  IRR(8) = RCT(8)*V(88)*V(92)
  IRR(9) = RCT(9)*V(93)*V(100)
  IRR(10) = RCT(10)*V(100)*V(100)*F(2)
  IRR(11) = RCT(11)*V(92)*V(93)
  IRR(12) = RCT(12)*V(42)
  IRR(13) = RCT(13)*V(42)*F(1)
  IRR(14) = RCT(14)*V(92)*V(93)
  IRR(15) = RCT(15)*V(93)
  IRR(16) = RCT(16)*V(93)
  IRR(17) = RCT(17)*V(88)
  IRR(18) = RCT(18)*V(88)
  IRR(19) = RCT(19)*V(26)*F(1)
  IRR(20) = RCT(20)*V(26)*F(2)
  IRR(21) = RCT(21)*V(98)*V(100)
  IRR(22) = RCT(22)*V(43)
  IRR(23) = RCT(23)*V(43)
  IRR(24) = RCT(24)*V(43)*V(98)
  IRR(25) = RCT(25)*V(92)*V(98)
  IRR(26) = RCT(26)*V(93)*V(98)
  IRR(27) = RCT(27)*V(66)*V(98)
  IRR(28) = RCT(28)*V(66)
  IRR(29) = RCT(29)*V(65)*V(98)
  IRR(30) = RCT(30)*V(88)*V(98)
  IRR(31) = RCT(31)*V(97)*V(100)
  IRR(32) = RCT(32)*V(92)*V(97)
  IRR(33) = RCT(33)*V(51)
  IRR(34) = RCT(34)*V(51)
  IRR(35) = RCT(35)*V(51)*V(98)
  IRR(36) = RCT(36)*V(88)*V(97)
  IRR(37) = RCT(37)*V(97)*V(97)
  IRR(38) = RCT(38)*V(97)*V(97)*F(1)
  IRR(39) = RCT(39)*V(93)*V(97)
  IRR(40) = RCT(40)*V(93)*V(93)
  IRR(41) = RCT(41)*V(38)
  IRR(42) = RCT(42)*V(38)*V(98)
  IRR(43) = RCT(43)*V(97)*V(98)
  IRR(44) = RCT(44)*V(32)*V(98)
  IRR(45) = RCT(45)*V(98)*F(2)
  IRR(46) = RCT(46)*V(94)*V(100)
  IRR(47) = RCT(47)*V(94)*V(97)
  IRR(48) = RCT(48)*V(93)*V(94)
  IRR(49) = RCT(49)*V(94)*V(94)
  IRR(50) = RCT(50)*V(94)*V(94)
  IRR(51) = RCT(51)*V(90)*V(100)
  IRR(52) = RCT(52)*V(90)*V(97)
  IRR(53) = RCT(53)*V(90)*V(93)
  IRR(54) = RCT(54)*V(90)*V(94)
  IRR(55) = RCT(55)*V(90)*V(90)
  IRR(56) = RCT(56)*V(73)*V(100)
  IRR(57) = RCT(57)*V(73)*V(97)
  IRR(58) = RCT(58)*V(73)*V(93)
  IRR(59) = RCT(59)*V(73)*V(94)
  IRR(60) = RCT(60)*V(73)*V(90)
  IRR(61) = RCT(61)*V(73)*V(73)
  IRR(62) = RCT(62)*V(91)*V(100)
  IRR(63) = RCT(63)*V(91)*V(97)
  IRR(64) = RCT(64)*V(91)*V(94)
  IRR(65) = RCT(65)*V(91)*V(93)
  IRR(66) = RCT(66)*V(90)*V(91)
  IRR(67) = RCT(67)*V(73)*V(91)
  IRR(68) = RCT(68)*V(91)*V(91)
  IRR(69) = RCT(69)*V(89)*V(92)
  IRR(70) = RCT(70)*V(34)
  IRR(71) = RCT(71)*V(89)*V(100)
  IRR(72) = RCT(72)*V(89)*V(97)
  IRR(73) = RCT(73)*V(89)*V(93)
  IRR(74) = RCT(74)*V(89)*V(94)
  IRR(75) = RCT(75)*V(89)*V(90)
  IRR(76) = RCT(76)*V(73)*V(89)
  IRR(77) = RCT(77)*V(89)*V(91)
  IRR(78) = RCT(78)*V(89)*V(89)
  IRR(79) = RCT(79)*V(92)*V(99)
  IRR(80) = RCT(80)*V(35)
  IRR(81) = RCT(81)*V(99)*V(100)
  IRR(82) = RCT(82)*V(97)*V(99)
  IRR(83) = RCT(83)*V(93)*V(99)
  IRR(84) = RCT(84)*V(94)*V(99)
  IRR(85) = RCT(85)*V(90)*V(99)
  IRR(86) = RCT(86)*V(73)*V(99)
  IRR(87) = RCT(87)*V(91)*V(99)
  IRR(88) = RCT(88)*V(89)*V(99)
  IRR(89) = RCT(89)*V(99)*V(99)
  IRR(90) = RCT(90)*V(92)*V(95)
  IRR(91) = RCT(91)*V(36)
  IRR(92) = RCT(92)*V(95)*V(100)
  IRR(93) = RCT(93)*V(95)*V(97)
  IRR(94) = RCT(94)*V(93)*V(95)
  IRR(95) = RCT(95)*V(94)*V(95)
  IRR(96) = RCT(96)*V(90)*V(95)
  IRR(97) = RCT(97)*V(73)*V(95)
  IRR(98) = RCT(98)*V(91)*V(95)
  IRR(99) = RCT(99)*V(89)*V(95)
  IRR(100) = RCT(100)*V(95)*V(99)
  IRR(101) = RCT(101)*V(95)*V(95)
  IRR(102) = RCT(102)*V(92)*V(96)
  IRR(103) = RCT(103)*V(37)
  IRR(104) = RCT(104)*V(96)*V(100)
  IRR(105) = RCT(105)*V(96)*V(97)
  IRR(106) = RCT(106)*V(93)*V(96)
  IRR(107) = RCT(107)*V(94)*V(96)
  IRR(108) = RCT(108)*V(90)*V(96)
  IRR(109) = RCT(109)*V(73)*V(96)
  IRR(110) = RCT(110)*V(91)*V(96)
  IRR(111) = RCT(111)*V(89)*V(96)
  IRR(112) = RCT(112)*V(96)*V(99)
  IRR(113) = RCT(113)*V(95)*V(96)
  IRR(114) = RCT(114)*V(96)*V(96)
  IRR(115) = RCT(115)*V(45)*V(92)
  IRR(116) = RCT(116)*V(45)
  IRR(117) = RCT(117)*V(70)*V(92)
  IRR(118) = RCT(118)*V(70)*V(97)
  IRR(119) = RCT(119)*V(70)
  IRR(120) = RCT(120)*V(49)*V(92)
  IRR(121) = RCT(121)*V(49)*V(97)
  IRR(122) = RCT(122)*V(49)
  IRR(123) = RCT(123)*V(82)
  IRR(124) = RCT(124)*V(82)
  IRR(125) = RCT(125)*V(82)*V(98)
  IRR(126) = RCT(126)*V(82)*V(97)
  IRR(127) = RCT(127)*V(48)
  IRR(128) = RCT(128)*V(48)*V(100)
  IRR(129) = RCT(129)*V(82)*V(93)
  IRR(130) = RCT(130)*V(79)*V(98)
  IRR(131) = RCT(131)*V(79)
  IRR(132) = RCT(132)*V(79)*V(93)
  IRR(133) = RCT(133)*V(85)*V(98)
  IRR(134) = RCT(134)*V(85)
  IRR(135) = RCT(135)*V(85)*V(93)
  IRR(136) = RCT(136)*V(68)*V(98)
  IRR(137) = RCT(137)*V(68)
  IRR(138) = RCT(138)*V(86)*V(98)
  IRR(139) = RCT(139)*V(86)
  IRR(140) = RCT(140)*V(52)*V(98)
  IRR(141) = RCT(141)*V(41)*V(98)
  IRR(142) = RCT(142)*V(47)*V(98)
  IRR(143) = RCT(143)*V(47)
  IRR(144) = RCT(144)*V(60)*V(98)
  IRR(145) = RCT(145)*V(60)
  IRR(146) = RCT(146)*V(77)
  IRR(147) = RCT(147)*V(77)
  IRR(148) = RCT(148)*V(77)*V(98)
  IRR(149) = RCT(149)*V(77)*V(93)
  IRR(150) = RCT(150)*V(64)
  IRR(151) = RCT(151)*V(64)*V(98)
  IRR(152) = RCT(152)*V(64)*V(93)
  IRR(153) = RCT(153)*V(40)
  IRR(154) = RCT(154)*V(63)*V(98)
  IRR(155) = RCT(155)*V(63)*V(93)
  IRR(156) = RCT(156)*V(57)*V(98)
  IRR(157) = RCT(157)*V(57)*V(93)
  IRR(158) = RCT(158)*V(61)*V(93)
  IRR(159) = RCT(159)*V(62)*V(98)
  IRR(160) = RCT(160)*V(62)
  IRR(161) = RCT(161)*V(62)*V(93)
  IRR(162) = RCT(162)*V(74)*V(98)
  IRR(163) = RCT(163)*V(74)*V(88)
  IRR(164) = RCT(164)*V(74)*V(93)
  IRR(165) = RCT(165)*V(74)*V(84)
  IRR(166) = RCT(166)*V(74)
  IRR(167) = RCT(167)*V(81)*V(98)
  IRR(168) = RCT(168)*V(81)*V(88)
  IRR(169) = RCT(169)*V(81)*V(84)
  IRR(170) = RCT(170)*V(81)
  IRR(171) = RCT(171)*V(78)*V(98)
  IRR(172) = RCT(172)*V(78)*V(88)
  IRR(173) = RCT(173)*V(78)*V(93)
  IRR(174) = RCT(174)*V(78)
  IRR(175) = RCT(175)*V(87)*V(98)
  IRR(176) = RCT(176)*V(87)
  IRR(177) = RCT(177)*V(83)*V(98)
  IRR(178) = RCT(178)*V(83)
  IRR(179) = RCT(179)*V(58)*V(98)
  IRR(180) = RCT(180)*V(58)*V(88)
  IRR(181) = RCT(181)*V(55)*V(98)
  IRR(182) = RCT(182)*V(55)
  IRR(183) = RCT(183)*V(56)*V(98)
  IRR(184) = RCT(184)*V(56)
  IRR(185) = RCT(185)*V(31)*V(98)
  IRR(186) = RCT(186)*V(67)*V(98)
  IRR(187) = RCT(187)*V(67)*V(88)
  IRR(188) = RCT(188)*V(67)*V(93)
  IRR(189) = RCT(189)*V(67)*V(84)
  IRR(190) = RCT(190)*V(71)*V(98)
  IRR(191) = RCT(191)*V(71)*V(88)
  IRR(192) = RCT(192)*V(71)*V(93)
  IRR(193) = RCT(193)*V(71)*V(84)
  IRR(194) = RCT(194)*V(76)*V(98)
  IRR(195) = RCT(195)*V(76)*V(88)
  IRR(196) = RCT(196)*V(76)*V(93)
  IRR(197) = RCT(197)*V(76)*V(84)
  IRR(198) = RCT(198)*V(75)*V(98)
  IRR(199) = RCT(199)*V(75)*V(88)
  IRR(200) = RCT(200)*V(75)*V(93)
  IRR(201) = RCT(201)*V(75)*V(84)
  IRR(202) = RCT(202)*V(33)*V(98)
  IRR(203) = RCT(203)*V(39)*V(98)
  IRR(204) = RCT(204)*V(59)*V(98)
  IRR(205) = RCT(205)*V(44)*V(98)
  IRR(206) = RCT(206)*V(53)*V(98)
  IRR(207) = RCT(207)*V(46)*V(98)
  IRR(208) = RCT(208)*V(54)*V(98)
  IRR(209) = RCT(209)*V(50)*V(98)
  IRR(210) = RCT(210)*V(72)*V(98)
  IRR(211) = RCT(211)*V(72)*V(88)
  IRR(212) = RCT(212)*V(72)*V(93)
  IRR(213) = RCT(213)*V(72)*V(84)
  IRR(214) = RCT(214)*V(80)*V(98)
  IRR(215) = RCT(215)*V(80)*V(88)
  IRR(216) = RCT(216)*V(80)*V(93)
  IRR(217) = RCT(217)*V(80)*V(84)
  IRR(218) = RCT(218)*V(59)*V(88)
  IRR(219) = RCT(219)*V(69)*V(98)
  IRR(220) = RCT(220)*V(69)*V(88)
  IRR(221) = RCT(221)*V(69)*V(93)
  IRR(222) = RCT(222)*V(69)*V(84)
  IRR(223) = RCT(223)*V(32)
  IRR(224) = RCT(224)*V(97)
  IRR(225) = RCT(225)*V(32)
  IRR(226) = RCT(226)*V(1)
  IRR(227) = RCT(227)*V(66)
  IRR(228) = RCT(228)*V(38)
  IRR(229) = RCT(229)*V(2)
  IRR(230) = RCT(230)*V(53)*V(98)
  IRR(231) = RCT(231)*V(46)*V(98)
  IRR(232) = RCT(232)*V(72)*V(98)
  IRR(233) = RCT(233)*V(80)*V(98)
  IRR(234) = RCT(234)*V(54)*V(98)
  IRR(235) = RCT(235)*V(50)*V(98)
  IRR(236) = RCT(236)*V(71)*V(98)
  IRR(237) = RCT(237)*V(76)*V(98)
  IRR(238) = RCT(238)*V(75)*V(98)
  IRR(239) = RCT(239)*V(3)*V(98)
  IRR(240) = RCT(240)*V(4)*V(98)
  IRR(241) = RCT(241)*V(28)*V(98)
  IRR(242) = RCT(242)*V(22)*V(98)
  IRR(243) = RCT(243)*V(5)*V(98)
  IRR(244) = RCT(244)*V(6)*V(98)
  IRR(245) = RCT(245)*V(30)*V(98)
  IRR(246) = RCT(246)*V(24)*V(98)
  IRR(247) = RCT(247)*V(27)*V(98)
  IRR(248) = RCT(248)*V(23)*V(98)
  IRR(249) = RCT(249)*V(29)*V(98)
  IRR(250) = RCT(250)*V(25)*V(98)
      
END SUBROUTINE saprc99_mosaic_4bin_vbs2_IRRFun
















SUBROUTINE saprc99_mosaic_4bin_vbs2_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(443)


  B(1) = RCT(1)

  B(2) = RCT(2)*F(2)

  B(4) = RCT(3)*V(88)

  B(5) = RCT(3)*V(84)

  B(6) = RCT(4)*V(100)*F(2)

  B(7) = RCT(4)*V(84)*F(2)

  B(9) = RCT(5)*V(92)

  B(10) = RCT(5)*V(84)

  B(11) = RCT(6)*V(92)

  B(12) = RCT(6)*V(84)

  B(13) = RCT(7)*V(100)

  B(14) = RCT(7)*V(88)

  B(15) = RCT(8)*V(92)

  B(16) = RCT(8)*V(88)

  B(17) = RCT(9)*V(100)

  B(18) = RCT(9)*V(93)

  B(19) = RCT(10)*2*V(100)*F(2)

  B(21) = RCT(11)*V(93)

  B(22) = RCT(11)*V(92)

  B(23) = RCT(12)

  B(24) = RCT(13)*F(1)

  B(26) = RCT(14)*V(93)

  B(27) = RCT(14)*V(92)

  B(28) = RCT(15)

  B(29) = RCT(16)

  B(30) = RCT(17)

  B(31) = RCT(18)

  B(32) = RCT(19)*F(1)

  B(34) = RCT(20)*F(2)

  B(36) = RCT(21)*V(100)

  B(37) = RCT(21)*V(98)

  B(38) = RCT(22)

  B(39) = RCT(23)

  B(40) = RCT(24)*V(98)

  B(41) = RCT(24)*V(43)

  B(42) = RCT(25)*V(98)

  B(43) = RCT(25)*V(92)

  B(44) = RCT(26)*V(98)

  B(45) = RCT(26)*V(93)

  B(46) = RCT(27)*V(98)

  B(47) = RCT(27)*V(66)

  B(48) = RCT(28)

  B(49) = RCT(29)*V(98)

  B(50) = RCT(29)*V(65)

  B(51) = RCT(30)*V(98)

  B(52) = RCT(30)*V(88)

  B(53) = RCT(31)*V(100)

  B(54) = RCT(31)*V(97)

  B(55) = RCT(32)*V(97)

  B(56) = RCT(32)*V(92)

  B(57) = RCT(33)

  B(58) = RCT(34)

  B(59) = RCT(35)*V(98)

  B(60) = RCT(35)*V(51)

  B(61) = RCT(36)*V(97)

  B(62) = RCT(36)*V(88)

  B(63) = RCT(37)*2*V(97)

  B(64) = RCT(38)*2*V(97)*F(1)

  B(66) = RCT(39)*V(97)

  B(67) = RCT(39)*V(93)

  B(68) = RCT(40)*2*V(93)

  B(69) = RCT(41)

  B(70) = RCT(42)*V(98)

  B(71) = RCT(42)*V(38)

  B(72) = RCT(43)*V(98)

  B(73) = RCT(43)*V(97)

  B(74) = RCT(44)*V(98)

  B(75) = RCT(44)*V(32)

  B(76) = RCT(45)*F(2)

  B(78) = RCT(46)*V(100)

  B(79) = RCT(46)*V(94)

  B(80) = RCT(47)*V(97)

  B(81) = RCT(47)*V(94)

  B(82) = RCT(48)*V(94)

  B(83) = RCT(48)*V(93)

  B(84) = RCT(49)*2*V(94)

  B(85) = RCT(50)*2*V(94)

  B(86) = RCT(51)*V(100)

  B(87) = RCT(51)*V(90)

  B(88) = RCT(52)*V(97)

  B(89) = RCT(52)*V(90)

  B(90) = RCT(53)*V(93)

  B(91) = RCT(53)*V(90)

  B(92) = RCT(54)*V(94)

  B(93) = RCT(54)*V(90)

  B(94) = RCT(55)*2*V(90)

  B(95) = RCT(56)*V(100)

  B(96) = RCT(56)*V(73)

  B(97) = RCT(57)*V(97)

  B(98) = RCT(57)*V(73)

  B(99) = RCT(58)*V(93)

  B(100) = RCT(58)*V(73)

  B(101) = RCT(59)*V(94)

  B(102) = RCT(59)*V(73)

  B(103) = RCT(60)*V(90)

  B(104) = RCT(60)*V(73)

  B(105) = RCT(61)*2*V(73)

  B(106) = RCT(62)*V(100)

  B(107) = RCT(62)*V(91)

  B(108) = RCT(63)*V(97)

  B(109) = RCT(63)*V(91)

  B(110) = RCT(64)*V(94)

  B(111) = RCT(64)*V(91)

  B(112) = RCT(65)*V(93)

  B(113) = RCT(65)*V(91)

  B(114) = RCT(66)*V(91)

  B(115) = RCT(66)*V(90)

  B(116) = RCT(67)*V(91)

  B(117) = RCT(67)*V(73)

  B(118) = RCT(68)*2*V(91)

  B(119) = RCT(69)*V(92)

  B(120) = RCT(69)*V(89)

  B(121) = RCT(70)

  B(122) = RCT(71)*V(100)

  B(123) = RCT(71)*V(89)

  B(124) = RCT(72)*V(97)

  B(125) = RCT(72)*V(89)

  B(126) = RCT(73)*V(93)

  B(127) = RCT(73)*V(89)

  B(128) = RCT(74)*V(94)

  B(129) = RCT(74)*V(89)

  B(130) = RCT(75)*V(90)

  B(131) = RCT(75)*V(89)

  B(132) = RCT(76)*V(89)

  B(133) = RCT(76)*V(73)

  B(134) = RCT(77)*V(91)

  B(135) = RCT(77)*V(89)

  B(136) = RCT(78)*2*V(89)

  B(137) = RCT(79)*V(99)

  B(138) = RCT(79)*V(92)

  B(139) = RCT(80)

  B(140) = RCT(81)*V(100)

  B(141) = RCT(81)*V(99)

  B(142) = RCT(82)*V(99)

  B(143) = RCT(82)*V(97)

  B(144) = RCT(83)*V(99)

  B(145) = RCT(83)*V(93)

  B(146) = RCT(84)*V(99)

  B(147) = RCT(84)*V(94)

  B(148) = RCT(85)*V(99)

  B(149) = RCT(85)*V(90)

  B(150) = RCT(86)*V(99)

  B(151) = RCT(86)*V(73)

  B(152) = RCT(87)*V(99)

  B(153) = RCT(87)*V(91)

  B(154) = RCT(88)*V(99)

  B(155) = RCT(88)*V(89)

  B(156) = RCT(89)*2*V(99)

  B(157) = RCT(90)*V(95)

  B(158) = RCT(90)*V(92)

  B(159) = RCT(91)

  B(160) = RCT(92)*V(100)

  B(161) = RCT(92)*V(95)

  B(162) = RCT(93)*V(97)

  B(163) = RCT(93)*V(95)

  B(164) = RCT(94)*V(95)

  B(165) = RCT(94)*V(93)

  B(166) = RCT(95)*V(95)

  B(167) = RCT(95)*V(94)

  B(168) = RCT(96)*V(95)

  B(169) = RCT(96)*V(90)

  B(170) = RCT(97)*V(95)

  B(171) = RCT(97)*V(73)

  B(172) = RCT(98)*V(95)

  B(173) = RCT(98)*V(91)

  B(174) = RCT(99)*V(95)

  B(175) = RCT(99)*V(89)

  B(176) = RCT(100)*V(99)

  B(177) = RCT(100)*V(95)

  B(178) = RCT(101)*2*V(95)

  B(179) = RCT(102)*V(96)

  B(180) = RCT(102)*V(92)

  B(181) = RCT(103)

  B(182) = RCT(104)*V(100)

  B(183) = RCT(104)*V(96)

  B(184) = RCT(105)*V(97)

  B(185) = RCT(105)*V(96)

  B(186) = RCT(106)*V(96)

  B(187) = RCT(106)*V(93)

  B(188) = RCT(107)*V(96)

  B(189) = RCT(107)*V(94)

  B(190) = RCT(108)*V(96)

  B(191) = RCT(108)*V(90)

  B(192) = RCT(109)*V(96)

  B(193) = RCT(109)*V(73)

  B(194) = RCT(110)*V(96)

  B(195) = RCT(110)*V(91)

  B(196) = RCT(111)*V(96)

  B(197) = RCT(111)*V(89)

  B(198) = RCT(112)*V(99)

  B(199) = RCT(112)*V(96)

  B(200) = RCT(113)*V(96)

  B(201) = RCT(113)*V(95)

  B(202) = RCT(114)*2*V(96)

  B(203) = RCT(115)*V(92)

  B(204) = RCT(115)*V(45)

  B(205) = RCT(116)

  B(206) = RCT(117)*V(92)

  B(207) = RCT(117)*V(70)

  B(208) = RCT(118)*V(97)

  B(209) = RCT(118)*V(70)

  B(210) = RCT(119)

  B(211) = RCT(120)*V(92)

  B(212) = RCT(120)*V(49)

  B(213) = RCT(121)*V(97)

  B(214) = RCT(121)*V(49)

  B(215) = RCT(122)

  B(216) = RCT(123)

  B(217) = RCT(124)

  B(218) = RCT(125)*V(98)

  B(219) = RCT(125)*V(82)

  B(220) = RCT(126)*V(97)

  B(221) = RCT(126)*V(82)

  B(222) = RCT(127)

  B(223) = RCT(128)*V(100)

  B(224) = RCT(128)*V(48)

  B(225) = RCT(129)*V(93)

  B(226) = RCT(129)*V(82)

  B(227) = RCT(130)*V(98)

  B(228) = RCT(130)*V(79)

  B(229) = RCT(131)

  B(230) = RCT(132)*V(93)

  B(231) = RCT(132)*V(79)

  B(232) = RCT(133)*V(98)

  B(233) = RCT(133)*V(85)

  B(234) = RCT(134)

  B(235) = RCT(135)*V(93)

  B(236) = RCT(135)*V(85)

  B(237) = RCT(136)*V(98)

  B(238) = RCT(136)*V(68)

  B(239) = RCT(137)

  B(240) = RCT(138)*V(98)

  B(241) = RCT(138)*V(86)

  B(242) = RCT(139)

  B(243) = RCT(140)*V(98)

  B(244) = RCT(140)*V(52)

  B(245) = RCT(141)*V(98)

  B(246) = RCT(141)*V(41)

  B(247) = RCT(142)*V(98)

  B(248) = RCT(142)*V(47)

  B(249) = RCT(143)

  B(250) = RCT(144)*V(98)

  B(251) = RCT(144)*V(60)

  B(252) = RCT(145)

  B(253) = RCT(146)

  B(254) = RCT(147)

  B(255) = RCT(148)*V(98)

  B(256) = RCT(148)*V(77)

  B(257) = RCT(149)*V(93)

  B(258) = RCT(149)*V(77)

  B(259) = RCT(150)

  B(260) = RCT(151)*V(98)

  B(261) = RCT(151)*V(64)

  B(262) = RCT(152)*V(93)

  B(263) = RCT(152)*V(64)

  B(264) = RCT(153)

  B(265) = RCT(154)*V(98)

  B(266) = RCT(154)*V(63)

  B(267) = RCT(155)*V(93)

  B(268) = RCT(155)*V(63)

  B(269) = RCT(156)*V(98)

  B(270) = RCT(156)*V(57)

  B(271) = RCT(157)*V(93)

  B(272) = RCT(157)*V(57)

  B(273) = RCT(158)*V(93)

  B(274) = RCT(158)*V(61)

  B(275) = RCT(159)*V(98)

  B(276) = RCT(159)*V(62)

  B(277) = RCT(160)

  B(278) = RCT(161)*V(93)

  B(279) = RCT(161)*V(62)

  B(280) = RCT(162)*V(98)

  B(281) = RCT(162)*V(74)

  B(282) = RCT(163)*V(88)

  B(283) = RCT(163)*V(74)

  B(284) = RCT(164)*V(93)

  B(285) = RCT(164)*V(74)

  B(286) = RCT(165)*V(84)

  B(287) = RCT(165)*V(74)

  B(288) = RCT(166)

  B(289) = RCT(167)*V(98)

  B(290) = RCT(167)*V(81)

  B(291) = RCT(168)*V(88)

  B(292) = RCT(168)*V(81)

  B(293) = RCT(169)*V(84)

  B(294) = RCT(169)*V(81)

  B(295) = RCT(170)

  B(296) = RCT(171)*V(98)

  B(297) = RCT(171)*V(78)

  B(298) = RCT(172)*V(88)

  B(299) = RCT(172)*V(78)

  B(300) = RCT(173)*V(93)

  B(301) = RCT(173)*V(78)

  B(302) = RCT(174)

  B(303) = RCT(175)*V(98)

  B(304) = RCT(175)*V(87)

  B(305) = RCT(176)

  B(306) = RCT(177)*V(98)

  B(307) = RCT(177)*V(83)

  B(308) = RCT(178)

  B(309) = RCT(179)*V(98)

  B(310) = RCT(179)*V(58)

  B(311) = RCT(180)*V(88)

  B(312) = RCT(180)*V(58)

  B(313) = RCT(181)*V(98)

  B(314) = RCT(181)*V(55)

  B(315) = RCT(182)

  B(316) = RCT(183)*V(98)

  B(317) = RCT(183)*V(56)

  B(318) = RCT(184)

  B(319) = RCT(185)*V(98)

  B(320) = RCT(185)*V(31)

  B(321) = RCT(186)*V(98)

  B(322) = RCT(186)*V(67)

  B(323) = RCT(187)*V(88)

  B(324) = RCT(187)*V(67)

  B(325) = RCT(188)*V(93)

  B(326) = RCT(188)*V(67)

  B(327) = RCT(189)*V(84)

  B(328) = RCT(189)*V(67)

  B(329) = RCT(190)*V(98)

  B(330) = RCT(190)*V(71)

  B(331) = RCT(191)*V(88)

  B(332) = RCT(191)*V(71)

  B(333) = RCT(192)*V(93)

  B(334) = RCT(192)*V(71)

  B(335) = RCT(193)*V(84)

  B(336) = RCT(193)*V(71)

  B(337) = RCT(194)*V(98)

  B(338) = RCT(194)*V(76)

  B(339) = RCT(195)*V(88)

  B(340) = RCT(195)*V(76)

  B(341) = RCT(196)*V(93)

  B(342) = RCT(196)*V(76)

  B(343) = RCT(197)*V(84)

  B(344) = RCT(197)*V(76)

  B(345) = RCT(198)*V(98)

  B(346) = RCT(198)*V(75)

  B(347) = RCT(199)*V(88)

  B(348) = RCT(199)*V(75)

  B(349) = RCT(200)*V(93)

  B(350) = RCT(200)*V(75)

  B(351) = RCT(201)*V(84)

  B(352) = RCT(201)*V(75)

  B(353) = RCT(202)*V(98)

  B(354) = RCT(202)*V(33)

  B(355) = RCT(203)*V(98)

  B(356) = RCT(203)*V(39)

  B(357) = RCT(204)*V(98)

  B(358) = RCT(204)*V(59)

  B(359) = RCT(205)*V(98)

  B(360) = RCT(205)*V(44)

  B(361) = RCT(206)*V(98)

  B(362) = RCT(206)*V(53)

  B(363) = RCT(207)*V(98)

  B(364) = RCT(207)*V(46)

  B(365) = RCT(208)*V(98)

  B(366) = RCT(208)*V(54)

  B(367) = RCT(209)*V(98)

  B(368) = RCT(209)*V(50)

  B(369) = RCT(210)*V(98)

  B(370) = RCT(210)*V(72)

  B(371) = RCT(211)*V(88)

  B(372) = RCT(211)*V(72)

  B(373) = RCT(212)*V(93)

  B(374) = RCT(212)*V(72)

  B(375) = RCT(213)*V(84)

  B(376) = RCT(213)*V(72)

  B(377) = RCT(214)*V(98)

  B(378) = RCT(214)*V(80)

  B(379) = RCT(215)*V(88)

  B(380) = RCT(215)*V(80)

  B(381) = RCT(216)*V(93)

  B(382) = RCT(216)*V(80)

  B(383) = RCT(217)*V(84)

  B(384) = RCT(217)*V(80)

  B(385) = RCT(218)*V(88)

  B(386) = RCT(218)*V(59)

  B(387) = RCT(219)*V(98)

  B(388) = RCT(219)*V(69)

  B(389) = RCT(220)*V(88)

  B(390) = RCT(220)*V(69)

  B(391) = RCT(221)*V(93)

  B(392) = RCT(221)*V(69)

  B(393) = RCT(222)*V(84)

  B(394) = RCT(222)*V(69)

  B(395) = RCT(223)

  B(396) = RCT(224)

  B(397) = RCT(225)

  B(398) = RCT(226)

  B(399) = RCT(227)

  B(400) = RCT(228)

  B(401) = RCT(229)

  B(402) = RCT(230)*V(98)

  B(403) = RCT(230)*V(53)

  B(404) = RCT(231)*V(98)

  B(405) = RCT(231)*V(46)

  B(406) = RCT(232)*V(98)

  B(407) = RCT(232)*V(72)

  B(408) = RCT(233)*V(98)

  B(409) = RCT(233)*V(80)

  B(410) = RCT(234)*V(98)

  B(411) = RCT(234)*V(54)

  B(412) = RCT(235)*V(98)

  B(413) = RCT(235)*V(50)

  B(414) = RCT(236)*V(98)

  B(415) = RCT(236)*V(71)

  B(416) = RCT(237)*V(98)

  B(417) = RCT(237)*V(76)

  B(418) = RCT(238)*V(98)

  B(419) = RCT(238)*V(75)

  B(420) = RCT(239)*V(98)

  B(421) = RCT(239)*V(3)

  B(422) = RCT(240)*V(98)

  B(423) = RCT(240)*V(4)

  B(424) = RCT(241)*V(98)

  B(425) = RCT(241)*V(28)

  B(426) = RCT(242)*V(98)

  B(427) = RCT(242)*V(22)

  B(428) = RCT(243)*V(98)

  B(429) = RCT(243)*V(5)

  B(430) = RCT(244)*V(98)

  B(431) = RCT(244)*V(6)

  B(432) = RCT(245)*V(98)

  B(433) = RCT(245)*V(30)

  B(434) = RCT(246)*V(98)

  B(435) = RCT(246)*V(24)

  B(436) = RCT(247)*V(98)

  B(437) = RCT(247)*V(27)

  B(438) = RCT(248)*V(98)

  B(439) = RCT(248)*V(23)

  B(440) = RCT(249)*V(98)

  B(441) = RCT(249)*V(29)

  B(442) = RCT(250)*V(98)

  B(443) = RCT(250)*V(25)



  JVS(1) = -B(398)

  JVS(2) = B(74)+B(395)

  JVS(3) = B(75)

  JVS(4) = -B(401)

  JVS(5) = 0.5*B(385)

  JVS(6) = 0.135*B(389)

  JVS(7) = 0.5*B(386)+0.135*B(390)

  JVS(8) = 0

  JVS(9) = 0

  JVS(10) = 0

  JVS(11) = 0

  JVS(12) = 0

  JVS(13) = B(223)

  JVS(14) = 0.297*B(357)

  JVS(15) = 0.37*B(323)

  JVS(16) = 0.185*B(389)

  JVS(17) = 0.204*B(331)

  JVS(18) = 0.185*B(371)

  JVS(19) = 0.333*B(282)

  JVS(20) = 0.103*B(347)

  JVS(21) = 0.103*B(339)

  JVS(22) = 0.1*B(298)

  JVS(23) = 0.073*B(379)

  JVS(24) = 0.351*B(291)

  JVS(25) = 0.333*B(283)+0.351*B(292)+0.1*B(299)+0.37*B(324)+0.204*B(332)+0.103*B(340)+0.103*B(348)+0.185*B(372)+0.073&
              &*B(380)+0.185*B(390)

  JVS(26) = 0.297*B(358)

  JVS(27) = B(224)

  JVS(28) = 0

  JVS(29) = 0.17*B(389)

  JVS(30) = 0.05*B(371)

  JVS(31) = 0.129*B(379)

  JVS(32) = 0.05*B(372)+0.129*B(380)+0.17*B(390)

  JVS(33) = 0.25*B(124)+B(128)+B(130)+B(134)

  JVS(34) = B(131)

  JVS(35) = B(135)

  JVS(36) = B(129)

  JVS(37) = 0.25*B(125)

  JVS(38) = 0

  JVS(39) = 0.15*B(331)

  JVS(40) = 0.119*B(371)

  JVS(41) = 0.189*B(347)

  JVS(42) = 0.189*B(339)

  JVS(43) = 0.372*B(298)

  JVS(44) = 0.247*B(379)

  JVS(45) = 0.372*B(299)+0.15*B(332)+0.189*B(340)+0.189*B(348)+0.119*B(372)+0.247*B(380)

  JVS(46) = B(148)+B(168)+B(190)

  JVS(47) = B(152)+B(172)+2*B(194)

  JVS(48) = B(146)+B(166)+B(188)

  JVS(49) = 0.25*B(162)+B(167)+B(169)+B(173)

  JVS(50) = 0.25*B(184)+B(189)+B(191)+2*B(195)

  JVS(51) = 0.25*B(142)+0.25*B(163)+0.25*B(185)

  JVS(52) = 0.25*B(143)+B(147)+B(149)+B(153)

  JVS(53) = 0

  JVS(54) = 0.75*B(124)

  JVS(55) = 0.75*B(125)

  JVS(56) = 0

  JVS(57) = 0.75*B(162)

  JVS(58) = 0.75*B(184)

  JVS(59) = 0.75*B(142)+0.75*B(163)+0.75*B(185)

  JVS(60) = 0.75*B(143)

  JVS(61) = 0

  JVS(62) = 2*B(211)

  JVS(63) = B(391)

  JVS(64) = 2*B(212)

  JVS(65) = B(392)

  JVS(66) = 0

  JVS(67) = 6*B(211)

  JVS(68) = 7*B(277)

  JVS(69) = 0.048*B(387)+0.07*B(389)+2.693*B(391)+0.55*B(393)

  JVS(70) = 0.55*B(394)

  JVS(71) = 0.07*B(390)

  JVS(72) = 6*B(212)

  JVS(73) = 2.693*B(392)

  JVS(74) = 0.048*B(388)

  JVS(75) = 0

  JVS(76) = B(95)+B(99)

  JVS(77) = B(86)+B(90)

  JVS(78) = B(106)+B(112)

  JVS(79) = B(82)+B(91)+B(100)+B(113)

  JVS(80) = B(78)+B(83)

  JVS(81) = B(79)+B(87)+B(96)+B(107)

  JVS(82) = 0

  JVS(83) = B(97)+B(101)+B(103)+B(105)+B(116)

  JVS(84) = B(88)+B(92)+B(94)+B(104)+B(114)

  JVS(85) = B(108)+B(110)+B(115)+B(117)+B(118)

  JVS(86) = B(80)+B(84)+B(85)+B(93)+B(102)+B(111)

  JVS(87) = B(81)+B(89)+B(98)+B(109)

  JVS(88) = 0

  JVS(89) = B(404)

  JVS(90) = B(402)

  JVS(91) = B(406)

  JVS(92) = B(408)

  JVS(93) = B(403)+B(405)+B(407)+B(409)

  JVS(94) = 0

  JVS(95) = B(412)

  JVS(96) = B(410)

  JVS(97) = B(411)+B(413)

  JVS(98) = 0

  JVS(99) = B(414)

  JVS(100) = B(415)

  JVS(101) = 0

  JVS(102) = B(418)

  JVS(103) = B(416)

  JVS(104) = B(417)+B(419)

  JVS(105) = 0

  JVS(106) = B(319)

  JVS(107) = B(74)

  JVS(108) = B(353)

  JVS(109) = B(355)

  JVS(110) = B(245)

  JVS(111) = B(40)

  JVS(112) = B(359)

  JVS(113) = B(363)

  JVS(114) = B(247)

  JVS(115) = B(367)

  JVS(116) = B(243)

  JVS(117) = B(361)

  JVS(118) = B(365)

  JVS(119) = B(313)

  JVS(120) = B(316)

  JVS(121) = B(269)

  JVS(122) = B(309)

  JVS(123) = B(357)

  JVS(124) = B(250)

  JVS(125) = B(275)

  JVS(126) = B(265)

  JVS(127) = B(260)

  JVS(128) = B(49)

  JVS(129) = B(46)

  JVS(130) = B(321)

  JVS(131) = B(237)

  JVS(132) = B(387)

  JVS(133) = B(329)

  JVS(134) = B(369)

  JVS(135) = B(280)

  JVS(136) = B(345)

  JVS(137) = B(337)

  JVS(138) = B(255)

  JVS(139) = B(296)

  JVS(140) = B(227)

  JVS(141) = B(377)

  JVS(142) = B(289)

  JVS(143) = B(218)

  JVS(144) = B(306)

  JVS(145) = B(232)

  JVS(146) = B(240)

  JVS(147) = B(303)

  JVS(148) = B(51)

  JVS(149) = B(42)

  JVS(150) = B(44)

  JVS(151) = B(72)

  JVS(152) = B(36)+B(41)+B(43)+B(45)+B(47)+B(50)+B(52)+B(73)+B(75)+B(76)+B(219)+B(228)+B(233)+B(238)+B(241)+B(244)&
               &+B(246)+B(248)+B(251)+B(256)+B(261)+B(266)+B(270)+B(276)+B(281)+B(290)+B(297)+B(304)+B(307)+B(310)+B(314)&
               &+B(317)+B(320)+B(322)+B(330)+B(338)+B(346)+B(354)+B(356)+B(358)+B(360)+B(362)+B(364)+B(366)+B(368)+B(370)&
               &+B(378)+B(388)

  JVS(153) = B(37)

  JVS(154) = B(420)

  JVS(155) = B(422)

  JVS(156) = B(428)

  JVS(157) = B(430)

  JVS(158) = 0

  JVS(159) = B(426)

  JVS(160) = B(438)

  JVS(161) = B(434)

  JVS(162) = B(442)

  JVS(163) = B(436)

  JVS(164) = B(424)

  JVS(165) = B(440)

  JVS(166) = B(432)

  JVS(167) = B(421)+B(423)+B(425)+B(427)+B(429)+B(431)+B(433)+B(435)+B(437)+B(439)+B(441)+B(443)

  JVS(168) = 0

  JVS(169) = B(438)

  JVS(170) = 0.5*B(436)

  JVS(171) = 0.5*B(437)+B(439)

  JVS(172) = -B(438)

  JVS(173) = -B(439)

  JVS(174) = 0

  JVS(175) = B(442)

  JVS(176) = 0.5*B(440)

  JVS(177) = 0.5*B(441)+B(443)

  JVS(178) = -B(442)

  JVS(179) = -B(443)

  JVS(180) = -B(32)-B(34)

  JVS(181) = B(31)

  JVS(182) = -B(436)

  JVS(183) = -B(437)

  JVS(184) = B(436)

  JVS(185) = 0

  JVS(186) = B(437)

  JVS(187) = -B(440)

  JVS(188) = -B(441)

  JVS(189) = B(440)

  JVS(190) = 0

  JVS(191) = B(441)

  JVS(192) = -B(319)

  JVS(193) = -B(320)

  JVS(194) = -B(74)-B(395)-B(397)

  JVS(195) = -B(75)

  JVS(196) = -B(353)

  JVS(197) = -B(354)

  JVS(198) = -B(121)

  JVS(199) = B(119)

  JVS(200) = B(120)

  JVS(201) = -B(139)

  JVS(202) = B(137)

  JVS(203) = B(138)

  JVS(204) = -B(159)

  JVS(205) = B(157)

  JVS(206) = B(158)

  JVS(207) = -B(181)

  JVS(208) = B(179)

  JVS(209) = B(180)

  JVS(210) = -B(69)-B(70)-B(400)

  JVS(211) = B(63)+B(64)

  JVS(212) = -B(71)

  JVS(213) = -B(355)

  JVS(214) = -B(356)

  JVS(215) = -B(264)

  JVS(216) = 0.087*B(367)

  JVS(217) = 0.031*B(347)

  JVS(218) = 0.031*B(339)

  JVS(219) = 0.031*B(340)+0.031*B(348)

  JVS(220) = 0.087*B(368)

  JVS(221) = -B(245)

  JVS(222) = -B(246)

  JVS(223) = -B(23)-B(24)

  JVS(224) = B(21)

  JVS(225) = B(22)

  JVS(226) = -B(38)-B(39)-B(40)

  JVS(227) = B(36)-B(41)

  JVS(228) = B(37)

  JVS(229) = -B(359)

  JVS(230) = -B(360)

  JVS(231) = 0.236*B(359)

  JVS(232) = -B(203)-B(205)

  JVS(233) = -B(204)

  JVS(234) = 0.236*B(360)

  JVS(235) = -B(363)

  JVS(236) = -B(364)

  JVS(237) = -B(247)-B(249)

  JVS(238) = B(80)

  JVS(239) = B(81)

  JVS(240) = -B(248)

  JVS(241) = -B(222)-B(223)

  JVS(242) = B(220)

  JVS(243) = B(221)

  JVS(244) = -B(224)

  JVS(245) = -B(211)-B(213)-B(215)

  JVS(246) = B(273)

  JVS(247) = -B(212)

  JVS(248) = B(274)

  JVS(249) = -B(214)

  JVS(250) = -B(367)

  JVS(251) = -B(368)

  JVS(252) = -B(57)-B(58)-B(59)

  JVS(253) = B(55)

  JVS(254) = B(56)

  JVS(255) = -B(60)

  JVS(256) = -B(243)

  JVS(257) = 0.25*B(92)

  JVS(258) = 0.25*B(110)

  JVS(259) = B(84)+0.25*B(93)+0.25*B(111)

  JVS(260) = -B(244)

  JVS(261) = -B(361)

  JVS(262) = -B(362)

  JVS(263) = -B(365)

  JVS(264) = -B(366)

  JVS(265) = 0.099*B(367)

  JVS(266) = 0.108*B(365)

  JVS(267) = -B(313)-B(315)

  JVS(268) = -B(314)+0.108*B(366)+0.099*B(368)

  JVS(269) = 0.093*B(367)

  JVS(270) = 0.051*B(365)

  JVS(271) = -B(316)-B(318)

  JVS(272) = -B(317)+0.051*B(366)+0.093*B(368)

  JVS(273) = 0.187*B(367)

  JVS(274) = 0.207*B(365)

  JVS(275) = -B(269)-B(271)

  JVS(276) = -B(272)

  JVS(277) = -B(270)+0.207*B(366)+0.187*B(368)

  JVS(278) = 0.561*B(367)

  JVS(279) = 0.491*B(365)

  JVS(280) = -B(309)-B(311)

  JVS(281) = -B(312)

  JVS(282) = -B(310)+0.491*B(366)+0.561*B(368)

  JVS(283) = -B(357)-B(385)

  JVS(284) = -B(386)

  JVS(285) = -B(358)

  JVS(286) = -B(250)-B(252)

  JVS(287) = B(88)

  JVS(288) = B(108)

  JVS(289) = B(89)+B(109)

  JVS(290) = -B(251)

  JVS(291) = B(213)+B(215)

  JVS(292) = -B(273)

  JVS(293) = B(206)

  JVS(294) = B(207)

  JVS(295) = -B(274)

  JVS(296) = B(214)

  JVS(297) = 0.05*B(367)

  JVS(298) = 0.059*B(365)

  JVS(299) = -B(275)-B(277)-B(278)

  JVS(300) = 0.061*B(377)+0.042*B(379)+0.015*B(381)

  JVS(301) = 0.042*B(380)

  JVS(302) = -B(279)+0.015*B(382)

  JVS(303) = -B(276)+0.059*B(366)+0.05*B(368)+0.061*B(378)

  JVS(304) = 0.017*B(365)

  JVS(305) = -B(265)-B(267)

  JVS(306) = B(208)+B(210)

  JVS(307) = -B(268)

  JVS(308) = B(209)

  JVS(309) = -B(266)+0.017*B(366)

  JVS(310) = 0.287*B(367)

  JVS(311) = 0.119*B(365)

  JVS(312) = 0.5*B(315)

  JVS(313) = 0.5*B(318)

  JVS(314) = 0.23*B(269)

  JVS(315) = -B(259)-B(260)-B(262)

  JVS(316) = 0.084*B(280)+0.9*B(282)

  JVS(317) = 0.174*B(296)+0.742*B(298)+0.008*B(300)

  JVS(318) = 0.3*B(289)+0.95*B(291)

  JVS(319) = 0.9*B(283)+0.95*B(292)+0.742*B(299)

  JVS(320) = -B(263)+0.008*B(301)

  JVS(321) = -B(261)+0.23*B(270)+0.084*B(281)+0.3*B(290)+0.174*B(297)+0.119*B(366)+0.287*B(368)

  JVS(322) = 0.002*B(361)

  JVS(323) = B(315)

  JVS(324) = B(318)

  JVS(325) = B(309)+1.5*B(311)

  JVS(326) = 0.393*B(357)+1.5*B(385)

  JVS(327) = B(259)+B(260)+B(262)

  JVS(328) = -B(49)

  JVS(329) = 0.5*B(323)+0.491*B(327)

  JVS(330) = 0.51*B(389)

  JVS(331) = 0.275*B(331)

  JVS(332) = 0.345*B(371)

  JVS(333) = 0.416*B(280)+0.45*B(282)+0.5*B(284)+0.67*B(288)

  JVS(334) = 0.157*B(347)

  JVS(335) = 0.157*B(339)

  JVS(336) = 2*B(253)+B(254)+1.26*B(255)+1.26*B(257)

  JVS(337) = 0.336*B(296)+0.498*B(298)+0.572*B(300)+1.233*B(302)

  JVS(338) = B(229)

  JVS(339) = 0.265*B(379)+0.012*B(383)

  JVS(340) = 0.475*B(291)+0.7*B(295)

  JVS(341) = B(216)+B(217)+B(218)+B(225)

  JVS(342) = 0.491*B(328)+0.012*B(384)

  JVS(343) = 0.034*B(232)+B(234)

  JVS(344) = 0.45*B(283)+0.475*B(292)+0.498*B(299)+1.5*B(312)+0.5*B(324)+0.275*B(332)+0.157*B(340)+0.157*B(348)+0.345&
               &*B(372)+0.265*B(380)+1.5*B(386)+0.51*B(390)

  JVS(345) = B(226)+1.26*B(258)+B(263)+0.5*B(285)+0.572*B(301)

  JVS(346) = -B(50)+B(219)+0.034*B(233)+1.26*B(256)+B(261)+0.416*B(281)+0.336*B(297)+B(310)+0.393*B(358)+0.002*B(362)

  JVS(347) = 2*B(24)

  JVS(348) = B(271)

  JVS(349) = B(273)

  JVS(350) = B(278)

  JVS(351) = B(267)

  JVS(352) = B(262)

  JVS(353) = -B(46)-B(48)-B(399)

  JVS(354) = 0

  JVS(355) = 0.5*B(284)

  JVS(356) = B(257)

  JVS(357) = 0.15*B(300)

  JVS(358) = B(230)

  JVS(359) = 0

  JVS(360) = 0

  JVS(361) = B(225)

  JVS(362) = B(235)

  JVS(363) = 0

  JVS(364) = B(42)

  JVS(365) = 0.2*B(66)+B(226)+B(231)+B(236)+B(258)+B(263)+B(268)+B(272)+B(274)+B(279)+0.5*B(285)+0.15*B(301)

  JVS(366) = 0.2*B(67)

  JVS(367) = B(43)-B(47)

  JVS(368) = -B(321)-B(323)-B(325)-B(327)

  JVS(369) = -B(328)

  JVS(370) = -B(324)

  JVS(371) = -B(326)

  JVS(372) = -B(322)

  JVS(373) = 0.704*B(355)

  JVS(374) = 0.024*B(359)

  JVS(375) = B(205)

  JVS(376) = 0.072*B(363)

  JVS(377) = 0.452*B(361)

  JVS(378) = -B(237)-B(239)

  JVS(379) = 0.005*B(369)+0.001*B(371)+0.024*B(373)

  JVS(380) = 0.13*B(347)

  JVS(381) = 0.13*B(339)

  JVS(382) = 0.127*B(377)+0.045*B(379)+0.102*B(381)

  JVS(383) = 0.006*B(306)+0.02*B(308)

  JVS(384) = 0.13*B(340)+0.13*B(348)+0.001*B(372)+0.045*B(380)

  JVS(385) = 0

  JVS(386) = 0.024*B(374)+0.102*B(382)

  JVS(387) = -B(238)+0.006*B(307)+0.704*B(356)+0.024*B(360)+0.452*B(362)+0.072*B(364)+0.005*B(370)+0.127*B(378)

  JVS(388) = -B(387)-B(389)-B(391)-B(393)

  JVS(389) = -B(394)

  JVS(390) = -B(390)

  JVS(391) = -B(392)

  JVS(392) = -B(388)

  JVS(393) = 0.24*B(269)+B(271)

  JVS(394) = 0.24*B(265)+B(267)

  JVS(395) = -B(206)-B(208)-B(210)

  JVS(396) = B(174)

  JVS(397) = -B(207)

  JVS(398) = B(164)+B(268)+B(272)

  JVS(399) = B(160)+B(165)+B(175)+B(176)+2*B(178)+B(200)

  JVS(400) = B(201)

  JVS(401) = -B(209)

  JVS(402) = 0.24*B(266)+0.24*B(270)

  JVS(403) = B(177)

  JVS(404) = B(161)

  JVS(405) = -B(329)-B(331)-B(333)-B(335)

  JVS(406) = -B(336)

  JVS(407) = -B(332)

  JVS(408) = -B(334)

  JVS(409) = -B(330)

  JVS(410) = -B(369)-B(371)-B(373)-B(375)

  JVS(411) = -B(376)

  JVS(412) = -B(372)

  JVS(413) = -B(374)

  JVS(414) = -B(370)

  JVS(415) = 0.559*B(359)

  JVS(416) = 0.948*B(363)

  JVS(417) = 0.936*B(361)

  JVS(418) = B(313)+B(315)

  JVS(419) = B(316)+B(318)

  JVS(420) = B(237)

  JVS(421) = 0.079*B(329)+0.126*B(331)+0.187*B(333)+0.24*B(335)

  JVS(422) = 0.205*B(369)+0.488*B(373)

  JVS(423) = -B(95)-B(97)-B(99)-B(101)-B(103)-B(116)-B(132)-B(150)-B(170)-B(192)

  JVS(424) = 0.5*B(345)+0.729*B(347)+0.75*B(349)

  JVS(425) = 0.5*B(337)+0.729*B(339)+0.75*B(341)

  JVS(426) = 0.001*B(377)+0.137*B(379)+0.711*B(381)

  JVS(427) = 0.675*B(289)

  JVS(428) = 0.596*B(306)+0.152*B(308)

  JVS(429) = 0.24*B(336)

  JVS(430) = 0.616*B(240)

  JVS(431) = 0.515*B(305)

  JVS(432) = 0.126*B(332)+0.729*B(340)+0.729*B(348)+0.137*B(380)

  JVS(433) = -B(133)+B(174)

  JVS(434) = -B(104)

  JVS(435) = -B(117)

  JVS(436) = 0

  JVS(437) = -B(100)+B(164)+0.187*B(334)+0.75*B(342)+0.75*B(350)+0.488*B(374)+0.711*B(382)

  JVS(438) = -B(102)

  JVS(439) = B(160)+B(165)-B(171)+B(175)+B(176)+2*B(178)+B(200)

  JVS(440) = -B(193)+B(201)

  JVS(441) = -B(98)

  JVS(442) = B(238)+0.616*B(241)+0.675*B(290)+0.596*B(307)+B(314)+B(317)+0.079*B(330)+0.5*B(338)+0.5*B(346)+0.559*B(360)&
               &+0.936*B(362)+0.948*B(364)+0.205*B(370)+0.001*B(378)

  JVS(443) = -B(151)+B(177)

  JVS(444) = -B(96)+B(161)

  JVS(445) = 0.23*B(329)+0.39*B(331)

  JVS(446) = -B(280)-B(282)-B(284)-B(286)-B(288)

  JVS(447) = 0.025*B(377)+0.026*B(379)+0.012*B(383)

  JVS(448) = -B(287)+0.012*B(384)

  JVS(449) = -B(283)+0.39*B(332)+0.026*B(380)

  JVS(450) = -B(285)

  JVS(451) = -B(281)+0.23*B(330)+0.025*B(378)

  JVS(452) = -B(345)-B(347)-B(349)-B(351)

  JVS(453) = -B(352)

  JVS(454) = -B(348)

  JVS(455) = -B(350)

  JVS(456) = -B(346)

  JVS(457) = -B(337)-B(339)-B(341)-B(343)

  JVS(458) = -B(344)

  JVS(459) = -B(340)

  JVS(460) = -B(342)

  JVS(461) = -B(338)

  JVS(462) = 0.097*B(367)

  JVS(463) = 0.118*B(365)

  JVS(464) = 0.5*B(315)

  JVS(465) = 0.5*B(318)

  JVS(466) = B(311)

  JVS(467) = 0.607*B(357)

  JVS(468) = 0.23*B(265)

  JVS(469) = 0.009*B(327)

  JVS(470) = 0

  JVS(471) = 0.001*B(347)

  JVS(472) = 0.001*B(339)

  JVS(473) = -B(253)-B(254)-B(255)-B(257)

  JVS(474) = 0.15*B(296)+0.023*B(298)

  JVS(475) = 0.009*B(328)

  JVS(476) = 0.023*B(299)+B(312)+0.001*B(340)+0.001*B(348)

  JVS(477) = 0

  JVS(478) = 0

  JVS(479) = -B(258)

  JVS(480) = 0

  JVS(481) = 0

  JVS(482) = 0

  JVS(483) = -B(256)+0.23*B(266)+0.15*B(297)+0.607*B(358)+0.118*B(366)+0.097*B(368)

  JVS(484) = 0

  JVS(485) = 0

  JVS(486) = 0.357*B(329)+0.936*B(333)

  JVS(487) = -B(296)-B(298)-B(300)-B(302)

  JVS(488) = 0.025*B(377)

  JVS(489) = 0

  JVS(490) = -B(299)

  JVS(491) = -B(301)+0.936*B(334)

  JVS(492) = -B(297)+0.357*B(330)+0.025*B(378)

  JVS(493) = B(353)

  JVS(494) = 0.96*B(245)

  JVS(495) = 0.445*B(359)

  JVS(496) = 0.099*B(363)

  JVS(497) = 0.455*B(361)

  JVS(498) = 0.195*B(321)+0.25*B(327)

  JVS(499) = 0.984*B(387)+0.5*B(389)

  JVS(500) = 0.294*B(369)+0.154*B(371)+0.009*B(373)

  JVS(501) = 0.129*B(296)+0.047*B(298)+0.467*B(302)

  JVS(502) = -B(227)-B(229)-B(230)

  JVS(503) = 0.732*B(377)+0.456*B(379)+0.507*B(381)

  JVS(504) = 0.439*B(306)+0.431*B(308)

  JVS(505) = 0.25*B(328)

  JVS(506) = 0.034*B(232)+B(234)

  JVS(507) = 0.482*B(240)+B(242)

  JVS(508) = 0.084*B(303)+0.246*B(305)

  JVS(509) = 0.047*B(299)+0.154*B(372)+0.456*B(380)+0.5*B(390)

  JVS(510) = B(154)

  JVS(511) = B(144)-B(231)+0.009*B(374)+0.507*B(382)

  JVS(512) = B(176)

  JVS(513) = B(198)

  JVS(514) = -B(228)+0.034*B(233)+0.482*B(241)+0.96*B(246)+0.129*B(297)+0.084*B(304)+0.439*B(307)+0.195*B(322)+B(354)&
               &+0.445*B(360)+0.455*B(362)+0.099*B(364)+0.294*B(370)+0.732*B(378)+0.984*B(388)

  JVS(515) = B(140)+B(145)+B(155)+2*B(156)+B(177)+B(199)

  JVS(516) = B(141)

  JVS(517) = -B(377)-B(379)-B(381)-B(383)

  JVS(518) = -B(384)

  JVS(519) = -B(380)

  JVS(520) = -B(382)

  JVS(521) = -B(378)

  JVS(522) = 0.32*B(329)+0.16*B(331)

  JVS(523) = 0.019*B(379)+0.048*B(381)

  JVS(524) = -B(289)-B(291)-B(293)-B(295)

  JVS(525) = -B(294)

  JVS(526) = -B(292)+0.16*B(332)+0.019*B(380)

  JVS(527) = 0.048*B(382)

  JVS(528) = -B(290)+0.32*B(330)

  JVS(529) = 0.081*B(245)

  JVS(530) = 0.026*B(359)

  JVS(531) = 0.026*B(363)

  JVS(532) = 0.35*B(247)+B(249)

  JVS(533) = B(222)

  JVS(534) = B(243)

  JVS(535) = 0.024*B(361)

  JVS(536) = 0.096*B(357)

  JVS(537) = 1.61*B(321)+B(323)+0.191*B(327)

  JVS(538) = B(237)

  JVS(539) = 0.984*B(387)+0.5*B(389)

  JVS(540) = 0.624*B(329)+0.592*B(331)+0.24*B(335)

  JVS(541) = 0.732*B(369)+0.5*B(371)

  JVS(542) = 0.084*B(280)+0.2*B(282)+0.67*B(288)

  JVS(543) = 0.276*B(345)+0.235*B(347)

  JVS(544) = 0.276*B(337)+0.235*B(339)

  JVS(545) = B(254)

  JVS(546) = 0.055*B(296)+0.125*B(298)+0.227*B(300)+0.3*B(302)

  JVS(547) = 0.244*B(377)+0.269*B(379)+0.079*B(381)

  JVS(548) = 0.3*B(289)+0.1*B(291)

  JVS(549) = -B(216)-B(217)-B(218)-B(220)-B(225)

  JVS(550) = 0.01*B(306)+0.134*B(308)

  JVS(551) = 0.191*B(328)+0.24*B(336)

  JVS(552) = 0.115*B(240)

  JVS(553) = 0.213*B(303)+0.506*B(305)

  JVS(554) = 0.2*B(283)+0.1*B(292)+0.125*B(299)+B(324)+0.592*B(332)+0.235*B(340)+0.235*B(348)+0.5*B(372)+0.269*B(380)&
               &+0.5*B(390)

  JVS(555) = B(128)+B(196)

  JVS(556) = 0.75*B(92)

  JVS(557) = 0.75*B(110)

  JVS(558) = 0

  JVS(559) = B(82)+B(186)-B(226)+0.227*B(301)+0.079*B(382)

  JVS(560) = B(78)+B(83)+B(84)+2*B(85)+0.75*B(93)+0.75*B(111)+B(129)+B(146)+B(166)+B(188)

  JVS(561) = B(167)+B(200)

  JVS(562) = B(182)+B(187)+B(189)+B(197)+B(198)+B(201)+2*B(202)

  JVS(563) = -B(221)

  JVS(564) = -B(219)+B(238)+0.115*B(241)+B(244)+0.081*B(246)+0.35*B(248)+0.084*B(281)+0.3*B(290)+0.055*B(297)+0.213&
               &*B(304)+0.01*B(307)+1.61*B(322)+0.624*B(330)+0.276*B(338)+0.276*B(346)+0.096*B(358)+0.026*B(360)+0.024&
               &*B(362)+0.026*B(364)+0.732*B(370)+0.244*B(378)+0.984*B(388)

  JVS(565) = B(147)+B(199)

  JVS(566) = B(79)+B(183)

  JVS(567) = B(203)

  JVS(568) = 0.511*B(373)

  JVS(569) = 0.276*B(349)

  JVS(570) = 0.276*B(341)

  JVS(571) = 0.572*B(300)

  JVS(572) = 0.321*B(381)

  JVS(573) = -0.69*B(306)-B(308)

  JVS(574) = 0

  JVS(575) = 0

  JVS(576) = B(106)

  JVS(577) = B(204)

  JVS(578) = 0.572*B(301)+0.276*B(342)+0.276*B(350)+0.511*B(374)+0.321*B(382)

  JVS(579) = -0.69*B(307)

  JVS(580) = B(107)

  JVS(581) = B(34)

  JVS(582) = -B(327)

  JVS(583) = -B(393)

  JVS(584) = -B(335)

  JVS(585) = -B(375)

  JVS(586) = -B(286)

  JVS(587) = -B(351)

  JVS(588) = -B(343)

  JVS(589) = -B(383)

  JVS(590) = -B(293)

  JVS(591) = -B(2)-B(4)-B(6)-B(9)-B(11)-B(287)-B(294)-B(328)-B(336)-B(344)-B(352)-B(376)-B(384)-B(394)

  JVS(592) = -B(5)+B(30)

  JVS(593) = B(1)-B(10)-B(12)

  JVS(594) = B(29)

  JVS(595) = 0

  JVS(596) = -B(7)

  JVS(597) = 0.261*B(355)

  JVS(598) = 0.122*B(359)

  JVS(599) = 0.204*B(363)

  JVS(600) = 0.244*B(361)

  JVS(601) = B(313)

  JVS(602) = B(316)

  JVS(603) = B(309)

  JVS(604) = B(250)+B(252)

  JVS(605) = B(325)

  JVS(606) = 0.45*B(393)

  JVS(607) = 0.497*B(369)+0.363*B(371)+0.037*B(373)+0.45*B(375)

  JVS(608) = B(286)

  JVS(609) = 0.474*B(345)+0.205*B(347)+0.474*B(349)+0.147*B(351)

  JVS(610) = 0.474*B(337)+0.205*B(339)+0.474*B(341)+0.147*B(343)

  JVS(611) = 0.013*B(296)+0.218*B(300)

  JVS(612) = 0.511*B(377)+0.305*B(379)+0.151*B(381)+0.069*B(383)

  JVS(613) = 0.675*B(289)+0.45*B(293)

  JVS(614) = 0.213*B(306)+0.147*B(308)

  JVS(615) = B(287)+0.45*B(294)+0.147*B(344)+0.147*B(352)+0.45*B(376)+0.069*B(384)+0.45*B(394)

  JVS(616) = -B(232)-B(234)-B(235)

  JVS(617) = 0.37*B(240)

  JVS(618) = 0.558*B(303)+0.71*B(305)

  JVS(619) = 0.205*B(340)+0.205*B(348)+0.363*B(372)+0.305*B(380)

  JVS(620) = 0

  JVS(621) = 0

  JVS(622) = 0

  JVS(623) = -B(236)+0.218*B(301)+B(326)+0.474*B(342)+0.474*B(350)+0.037*B(374)+0.151*B(382)

  JVS(624) = 0

  JVS(625) = -B(233)+0.37*B(241)+B(251)+0.675*B(290)+0.013*B(297)+0.558*B(304)+0.213*B(307)+B(310)+B(314)+B(317)+0.474&
               &*B(338)+0.474*B(346)+0.261*B(356)+0.122*B(360)+0.244*B(362)+0.204*B(364)+0.497*B(370)+0.511*B(378)

  JVS(626) = 0

  JVS(627) = 0.332*B(359)

  JVS(628) = 0.089*B(363)

  JVS(629) = 0.11*B(361)

  JVS(630) = 0.55*B(393)

  JVS(631) = 0.437*B(375)

  JVS(632) = 0.416*B(280)

  JVS(633) = 0.15*B(296)+0.21*B(298)+0.233*B(302)

  JVS(634) = 0.072*B(377)+0.026*B(379)+0.001*B(381)+0.659*B(383)

  JVS(635) = 0.55*B(293)

  JVS(636) = 0.177*B(306)+0.243*B(308)

  JVS(637) = 0.55*B(294)+0.437*B(376)+0.659*B(384)+0.55*B(394)

  JVS(638) = -B(240)-B(242)

  JVS(639) = 0.115*B(303)

  JVS(640) = 0.21*B(299)+0.026*B(380)

  JVS(641) = 0.5*B(114)

  JVS(642) = 0.5*B(110)+B(112)+0.5*B(115)+B(118)

  JVS(643) = 0

  JVS(644) = B(113)+0.001*B(382)

  JVS(645) = 0.5*B(111)

  JVS(646) = -B(241)+0.416*B(281)+0.15*B(297)+0.115*B(304)+0.177*B(307)+0.332*B(360)+0.11*B(362)+0.089*B(364)+0.072&
               &*B(378)

  JVS(647) = 0

  JVS(648) = 0.417*B(363)

  JVS(649) = 0.125*B(361)

  JVS(650) = 0.055*B(365)

  JVS(651) = 0.1*B(331)+0.75*B(335)

  JVS(652) = 0.119*B(369)+0.215*B(371)+0.113*B(375)

  JVS(653) = 0.276*B(345)+0.276*B(347)+0.853*B(351)

  JVS(654) = 0.276*B(337)+0.276*B(339)+0.853*B(343)

  JVS(655) = 0.332*B(296)

  JVS(656) = 0.043*B(379)+0.259*B(383)

  JVS(657) = 0.7*B(295)

  JVS(658) = 0.048*B(306)+0.435*B(308)

  JVS(659) = 0.75*B(336)+0.853*B(344)+0.853*B(352)+0.113*B(376)+0.259*B(384)

  JVS(660) = -0.671*B(303)-B(305)

  JVS(661) = 0.1*B(332)+0.276*B(340)+0.276*B(348)+0.215*B(372)+0.043*B(380)

  JVS(662) = B(134)

  JVS(663) = 0.5*B(114)

  JVS(664) = 0.5*B(110)+0.5*B(115)+B(118)+B(135)+B(152)+B(172)

  JVS(665) = 0

  JVS(666) = 0

  JVS(667) = 0.5*B(111)

  JVS(668) = B(173)

  JVS(669) = 0.332*B(297)-0.671*B(304)+0.048*B(307)+0.276*B(338)+0.276*B(346)+0.125*B(362)+0.417*B(364)+0.055*B(366)&
               &+0.119*B(370)

  JVS(670) = B(153)

  JVS(671) = 0

  JVS(672) = -B(311)

  JVS(673) = -B(385)

  JVS(674) = -B(323)

  JVS(675) = -B(389)

  JVS(676) = -B(331)

  JVS(677) = -B(371)

  JVS(678) = -B(282)

  JVS(679) = -B(347)

  JVS(680) = -B(339)

  JVS(681) = -B(298)

  JVS(682) = -B(379)

  JVS(683) = -B(291)

  JVS(684) = B(2)-B(4)

  JVS(685) = -B(5)-B(13)-B(15)-B(30)-B(31)-B(51)-B(61)-B(283)-B(292)-B(299)-B(312)-B(324)-B(332)-B(340)-B(348)-B(372)&
               &-B(380)-B(386)-B(390)

  JVS(686) = 0.25*B(124)

  JVS(687) = -B(16)

  JVS(688) = 0

  JVS(689) = 0.25*B(162)

  JVS(690) = 0.25*B(184)

  JVS(691) = -B(62)+0.25*B(125)+0.25*B(142)+0.25*B(163)+0.25*B(185)

  JVS(692) = -B(52)

  JVS(693) = 0.25*B(143)

  JVS(694) = -B(14)

  JVS(695) = B(121)

  JVS(696) = 2*B(264)

  JVS(697) = 0

  JVS(698) = 0.011*B(361)

  JVS(699) = B(313)+0.5*B(315)

  JVS(700) = B(316)+0.5*B(318)

  JVS(701) = B(259)+B(260)+B(262)

  JVS(702) = B(237)+B(239)

  JVS(703) = 0

  JVS(704) = 0.67*B(288)

  JVS(705) = 0.123*B(347)

  JVS(706) = 0.123*B(339)

  JVS(707) = 0.467*B(302)

  JVS(708) = B(227)+B(230)

  JVS(709) = 0.137*B(379)

  JVS(710) = 0.675*B(289)

  JVS(711) = 0

  JVS(712) = 0

  JVS(713) = 0

  JVS(714) = 0.492*B(240)+B(242)

  JVS(715) = 0.029*B(303)+0.667*B(305)

  JVS(716) = 0.123*B(340)+0.123*B(348)+0.137*B(380)

  JVS(717) = -B(119)-B(122)-B(124)-B(126)-B(128)-B(130)-B(134)-2*B(136)-B(154)-B(174)

  JVS(718) = -B(131)

  JVS(719) = -B(135)

  JVS(720) = -B(120)

  JVS(721) = -B(127)+B(186)+B(231)+B(263)

  JVS(722) = -B(129)

  JVS(723) = -B(175)+B(200)

  JVS(724) = B(182)+B(187)+B(198)+B(201)+2*B(202)

  JVS(725) = -B(125)

  JVS(726) = B(228)+B(238)+0.492*B(241)+B(261)+0.675*B(290)+0.029*B(304)+B(314)+B(317)+0.011*B(362)

  JVS(727) = -B(155)+B(199)

  JVS(728) = -B(123)+B(183)

  JVS(729) = B(353)

  JVS(730) = 0.965*B(355)

  JVS(731) = 0.05*B(245)

  JVS(732) = 0.695*B(359)

  JVS(733) = 0.653*B(363)

  JVS(734) = 0.804*B(367)

  JVS(735) = 0.835*B(361)

  JVS(736) = 0.765*B(365)

  JVS(737) = B(315)

  JVS(738) = B(318)

  JVS(739) = 0.76*B(269)

  JVS(740) = B(309)

  JVS(741) = 0.1*B(357)

  JVS(742) = 0.34*B(250)

  JVS(743) = 0.76*B(265)

  JVS(744) = B(321)+B(325)+0.2*B(327)

  JVS(745) = 0.984*B(387)+0.949*B(391)

  JVS(746) = 0

  JVS(747) = 0.907*B(329)+0.066*B(331)+0.749*B(333)

  JVS(748) = 0.91*B(369)+0.022*B(371)+0.824*B(373)

  JVS(749) = 0.5*B(280)+0.1*B(282)+0.5*B(284)+0.33*B(288)

  JVS(750) = 0.75*B(345)+0.031*B(347)+0.276*B(349)

  JVS(751) = 0.75*B(337)+0.031*B(339)+0.276*B(341)

  JVS(752) = 0.67*B(296)+0.048*B(298)+0.799*B(300)

  JVS(753) = 0.918*B(377)+0.033*B(379)+0.442*B(381)+0.012*B(383)

  JVS(754) = 0.3*B(289)+0.05*B(291)

  JVS(755) = 0.376*B(306)+0.564*B(308)

  JVS(756) = 0.2*B(328)+0.012*B(384)

  JVS(757) = 0.034*B(232)+B(234)

  JVS(758) = 0.37*B(240)+B(242)

  JVS(759) = 0.473*B(303)+0.96*B(305)

  JVS(760) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.066*B(332)+0.031*B(340)+0.031*B(348)+0.022*B(372)+0.033*B(380)

  JVS(761) = -B(130)+B(154)

  JVS(762) = -B(86)-B(88)-B(90)-B(92)-2*B(94)-B(114)-B(131)-B(148)-B(168)-B(190)

  JVS(763) = -B(115)

  JVS(764) = 0

  JVS(765) = -B(91)+B(144)+0.5*B(285)+0.799*B(301)+B(326)+0.749*B(334)+0.276*B(342)+0.276*B(350)+0.824*B(374)+0.442&
               &*B(382)+0.949*B(392)

  JVS(766) = -B(93)

  JVS(767) = -B(169)+B(176)

  JVS(768) = -B(191)+B(198)

  JVS(769) = -B(89)

  JVS(770) = 0.034*B(233)+0.37*B(241)+0.05*B(246)+0.34*B(251)+0.76*B(266)+0.76*B(270)+0.5*B(281)+0.3*B(290)+0.67*B(297)&
               &+0.473*B(304)+0.376*B(307)+B(310)+B(322)+0.907*B(330)+0.75*B(338)+0.75*B(346)+B(354)+0.965*B(356)+0.1*B(358)&
               &+0.695*B(360)+0.835*B(362)+0.653*B(364)+0.765*B(366)+0.804*B(368)+0.91*B(370)+0.918*B(378)+0.984*B(388)

  JVS(771) = B(140)+B(145)-B(149)+B(155)+2*B(156)+B(177)+B(199)

  JVS(772) = -B(87)+B(141)

  JVS(773) = 0.035*B(355)

  JVS(774) = 0.07*B(359)

  JVS(775) = 0.347*B(363)

  JVS(776) = 0.009*B(367)

  JVS(777) = 0.143*B(361)

  JVS(778) = 0.011*B(365)

  JVS(779) = 0.016*B(387)+0.051*B(391)

  JVS(780) = 0.093*B(329)+0.008*B(331)+0.064*B(333)+0.01*B(335)

  JVS(781) = 0.09*B(369)+0.001*B(371)+0.176*B(373)

  JVS(782) = 0.25*B(345)+0.18*B(347)+0.25*B(349)

  JVS(783) = 0.25*B(337)+0.18*B(339)+0.25*B(341)

  JVS(784) = 0.041*B(296)+0.051*B(300)

  JVS(785) = 0.082*B(377)+0.002*B(379)+0.136*B(381)+0.001*B(383)

  JVS(786) = 0.025*B(289)

  JVS(787) = 0.173*B(306)+0.095*B(308)

  JVS(788) = 0.01*B(336)+0.001*B(384)

  JVS(789) = 0.001*B(232)

  JVS(790) = 0.042*B(240)

  JVS(791) = 0.07*B(303)+0.04*B(305)

  JVS(792) = 0.008*B(332)+0.18*B(340)+0.18*B(348)+0.001*B(372)+0.002*B(380)

  JVS(793) = -B(134)

  JVS(794) = -B(114)

  JVS(795) = -B(106)-B(108)-B(110)-B(112)-B(115)-2*B(118)-B(135)-B(152)-B(172)-B(194)

  JVS(796) = 0

  JVS(797) = -B(113)+0.051*B(301)+0.064*B(334)+0.25*B(342)+0.25*B(350)+0.176*B(374)+0.136*B(382)+0.051*B(392)

  JVS(798) = -B(111)

  JVS(799) = -B(173)

  JVS(800) = -B(195)

  JVS(801) = -B(109)

  JVS(802) = 0.001*B(233)+0.042*B(241)+0.025*B(290)+0.041*B(297)+0.07*B(304)+0.173*B(307)+0.093*B(330)+0.25*B(338)+0.25&
               &*B(346)+0.035*B(356)+0.07*B(360)+0.143*B(362)+0.347*B(364)+0.011*B(366)+0.009*B(368)+0.09*B(370)+0.082&
               &*B(378)+0.016*B(388)

  JVS(803) = -B(153)

  JVS(804) = -B(107)

  JVS(805) = B(121)

  JVS(806) = B(139)

  JVS(807) = B(159)

  JVS(808) = B(181)

  JVS(809) = B(23)

  JVS(810) = B(39)+B(40)

  JVS(811) = -B(203)

  JVS(812) = B(223)

  JVS(813) = -B(211)

  JVS(814) = B(57)+0.61*B(58)+B(59)

  JVS(815) = 0

  JVS(816) = B(48)

  JVS(817) = -B(206)

  JVS(818) = 0.187*B(333)

  JVS(819) = B(95)+B(99)

  JVS(820) = 0

  JVS(821) = 0.474*B(349)

  JVS(822) = 0.474*B(341)

  JVS(823) = 0

  JVS(824) = 0

  JVS(825) = 0

  JVS(826) = 0.391*B(381)

  JVS(827) = 0

  JVS(828) = 0

  JVS(829) = 0.338*B(306)+B(308)

  JVS(830) = B(6)-B(9)-B(11)

  JVS(831) = 0

  JVS(832) = 0

  JVS(833) = 0

  JVS(834) = B(13)-B(15)

  JVS(835) = -B(119)+B(122)+B(126)

  JVS(836) = B(86)+B(90)

  JVS(837) = B(112)

  JVS(838) = -B(1)-B(10)-B(12)-B(16)-B(21)-B(42)-B(55)-B(120)-B(137)-B(157)-B(179)-B(204)-B(207)-B(212)

  JVS(839) = 2*B(17)-B(22)+B(29)+B(44)+0.8*B(66)+2*B(68)+B(82)+B(91)+B(100)+B(113)+B(127)+B(144)+B(164)+B(186)+0.187&
               &*B(334)+0.474*B(342)+0.474*B(350)+0.391*B(382)

  JVS(840) = B(78)+B(83)

  JVS(841) = -B(158)+B(160)+B(165)

  JVS(842) = -B(180)+B(182)+B(187)

  JVS(843) = B(53)-B(56)+0.8*B(67)

  JVS(844) = B(41)-B(43)+B(45)+B(60)+0.338*B(307)

  JVS(845) = -B(138)+B(140)+B(145)

  JVS(846) = B(7)+B(14)+2*B(18)+2*B(19)+B(54)+B(79)+B(87)+B(96)+B(123)+B(141)+B(161)+B(183)+B(224)

  JVS(847) = B(23)

  JVS(848) = 0.39*B(58)

  JVS(849) = -B(271)

  JVS(850) = -B(273)

  JVS(851) = -B(278)

  JVS(852) = -B(267)

  JVS(853) = -B(262)

  JVS(854) = B(46)

  JVS(855) = -B(325)

  JVS(856) = -B(391)

  JVS(857) = 0

  JVS(858) = -B(333)

  JVS(859) = -B(373)

  JVS(860) = -B(99)

  JVS(861) = -B(284)

  JVS(862) = -B(349)

  JVS(863) = -B(341)

  JVS(864) = -B(257)

  JVS(865) = -B(300)

  JVS(866) = -B(230)

  JVS(867) = -B(381)

  JVS(868) = 0

  JVS(869) = -B(225)

  JVS(870) = 0

  JVS(871) = B(11)

  JVS(872) = -B(235)

  JVS(873) = 0

  JVS(874) = 0

  JVS(875) = B(15)

  JVS(876) = -B(126)

  JVS(877) = -B(90)

  JVS(878) = -B(112)

  JVS(879) = B(12)+B(16)-B(21)-B(26)

  JVS(880) = -B(17)-B(22)-B(27)-B(28)-B(29)-B(44)-B(66)-2*B(68)-B(82)-B(91)-B(100)-B(113)-B(127)-B(144)-B(164)-B(186)&
               &-B(226)-B(231)-B(236)-B(258)-B(263)-B(268)-B(272)-B(274)-B(279)-B(285)-B(301)-B(326)-B(334)-B(342)-B(350)&
               &-B(374)-B(382)-B(392)

  JVS(881) = -B(83)

  JVS(882) = -B(165)

  JVS(883) = -B(187)

  JVS(884) = -B(67)

  JVS(885) = -B(45)+B(47)

  JVS(886) = -B(145)

  JVS(887) = -B(18)

  JVS(888) = B(319)

  JVS(889) = B(205)

  JVS(890) = 0.65*B(247)

  JVS(891) = 0.011*B(361)

  JVS(892) = 0.3*B(327)

  JVS(893) = B(239)

  JVS(894) = 0.26*B(389)

  JVS(895) = 0.25*B(335)

  JVS(896) = 0.076*B(371)

  JVS(897) = 0

  JVS(898) = 0

  JVS(899) = B(229)

  JVS(900) = 0.197*B(379)+0.03*B(381)

  JVS(901) = 0.3*B(295)

  JVS(902) = 0

  JVS(903) = 0.3*B(328)+0.25*B(336)

  JVS(904) = 0

  JVS(905) = 0

  JVS(906) = 0

  JVS(907) = 0.076*B(372)+0.197*B(380)+0.26*B(390)

  JVS(908) = B(122)+B(126)-B(128)+2*B(136)+B(154)+B(174)+B(196)

  JVS(909) = -B(92)

  JVS(910) = -B(110)

  JVS(911) = 0

  JVS(912) = -B(82)+B(127)+0.03*B(382)

  JVS(913) = -B(78)-B(80)-B(83)-2*B(84)-2*B(85)-B(93)-B(111)-B(129)-B(146)-B(166)-B(188)

  JVS(914) = -B(167)+B(175)

  JVS(915) = -B(189)+B(197)

  JVS(916) = -B(81)

  JVS(917) = 0.65*B(248)+B(320)+0.011*B(362)

  JVS(918) = -B(147)+B(155)

  JVS(919) = -B(79)+B(123)

  JVS(920) = B(159)

  JVS(921) = B(275)+B(278)

  JVS(922) = 0

  JVS(923) = 0

  JVS(924) = 0

  JVS(925) = -B(174)

  JVS(926) = -B(168)

  JVS(927) = -B(172)

  JVS(928) = -B(157)

  JVS(929) = -B(164)+B(279)

  JVS(930) = -B(166)

  JVS(931) = -B(158)-B(160)-B(162)-B(165)-B(167)-B(169)-B(173)-B(175)-B(176)-2*B(178)-B(200)

  JVS(932) = -B(201)

  JVS(933) = -B(163)

  JVS(934) = B(276)

  JVS(935) = -B(177)

  JVS(936) = -B(161)

  JVS(937) = B(181)

  JVS(938) = 0.192*B(331)+0.24*B(335)

  JVS(939) = 0.5*B(280)+0.5*B(284)+0.33*B(288)

  JVS(940) = 0.289*B(296)+0.15*B(300)

  JVS(941) = 0

  JVS(942) = 0.3*B(295)

  JVS(943) = 0.24*B(336)

  JVS(944) = 0.192*B(332)

  JVS(945) = -B(196)

  JVS(946) = -B(190)

  JVS(947) = -B(194)

  JVS(948) = -B(179)

  JVS(949) = -B(186)+0.5*B(285)+0.15*B(301)

  JVS(950) = -B(188)

  JVS(951) = -B(200)

  JVS(952) = -B(180)-B(182)-B(184)-B(187)-B(189)-B(191)-B(195)-B(197)-B(198)-B(201)-2*B(202)

  JVS(953) = -B(185)

  JVS(954) = 0.5*B(281)+0.289*B(297)

  JVS(955) = -B(199)

  JVS(956) = -B(183)

  JVS(957) = B(74)

  JVS(958) = B(70)

  JVS(959) = 0.95*B(245)

  JVS(960) = B(39)

  JVS(961) = B(249)

  JVS(962) = B(222)+B(223)

  JVS(963) = -B(213)

  JVS(964) = 0.187*B(367)

  JVS(965) = B(57)+0.61*B(58)

  JVS(966) = B(243)

  JVS(967) = 0.224*B(365)

  JVS(968) = 0.5*B(315)

  JVS(969) = 0.5*B(318)

  JVS(970) = 1.5*B(311)

  JVS(971) = 0.297*B(357)+1.5*B(385)

  JVS(972) = B(252)

  JVS(973) = 0

  JVS(974) = B(259)

  JVS(975) = B(49)

  JVS(976) = 0.12*B(323)+0.5*B(327)

  JVS(977) = 0.06*B(389)

  JVS(978) = -B(208)

  JVS(979) = 0

  JVS(980) = 0.056*B(371)

  JVS(981) = 0.008*B(282)+0.34*B(288)

  JVS(982) = 0.033*B(347)

  JVS(983) = 0.033*B(339)

  JVS(984) = 2*B(253)+0.63*B(255)+0.63*B(257)

  JVS(985) = 0.4*B(298)+1.233*B(302)

  JVS(986) = B(229)

  JVS(987) = 0.003*B(379)+0.013*B(383)

  JVS(988) = 0.064*B(291)

  JVS(989) = 2*B(216)+B(218)-B(220)+B(225)

  JVS(990) = 0.113*B(306)+0.341*B(308)

  JVS(991) = 0.5*B(328)+0.013*B(384)

  JVS(992) = B(234)

  JVS(993) = 0

  JVS(994) = 0.379*B(303)

  JVS(995) = B(51)-B(61)+0.008*B(283)+0.064*B(292)+0.4*B(299)+1.5*B(312)+0.12*B(324)+0.033*B(340)+0.033*B(348)+0.056&
               &*B(372)+0.003*B(380)+1.5*B(386)+0.06*B(390)

  JVS(996) = -B(124)

  JVS(997) = B(86)-B(88)+B(90)+B(92)+B(94)+B(114)

  JVS(998) = -B(108)+B(110)+B(112)+B(115)+B(118)

  JVS(999) = -B(55)

  JVS(1000) = B(44)-B(66)+B(82)+B(91)+B(113)+B(226)+0.63*B(258)

  JVS(1001) = B(78)-B(80)+B(83)+2*B(85)+B(93)+B(111)

  JVS(1002) = -B(162)

  JVS(1003) = -B(184)

  JVS(1004) = -B(53)-B(56)-B(62)-2*B(63)-2*B(64)-B(67)-B(72)-B(81)-B(89)-B(109)-B(125)-B(142)-B(163)-B(185)-B(209)&
                &-B(214)-B(221)-B(396)

  JVS(1005) = B(45)+B(50)+B(52)+B(71)-B(73)+B(75)+B(76)+B(219)+B(244)+0.95*B(246)+0.63*B(256)+0.379*B(304)+0.113*B(307)&
                &+0.297*B(358)+0.224*B(366)+0.187*B(368)

  JVS(1006) = -B(143)

  JVS(1007) = -B(54)+B(79)+B(87)+B(224)

  JVS(1008) = -B(420)

  JVS(1009) = -B(428)

  JVS(1010) = 2*B(32)

  JVS(1011) = -B(436)

  JVS(1012) = -B(424)

  JVS(1013) = -B(440)

  JVS(1014) = -B(432)

  JVS(1015) = -B(319)

  JVS(1016) = -B(74)

  JVS(1017) = -B(353)

  JVS(1018) = 2*B(69)-B(70)

  JVS(1019) = -B(355)

  JVS(1020) = -B(245)

  JVS(1021) = B(38)-B(40)

  JVS(1022) = -B(359)

  JVS(1023) = -B(363)

  JVS(1024) = -0.65*B(247)+B(249)

  JVS(1025) = -B(367)

  JVS(1026) = 0.39*B(58)-B(59)

  JVS(1027) = -B(243)

  JVS(1028) = -B(361)

  JVS(1029) = -B(365)

  JVS(1030) = -B(313)

  JVS(1031) = -B(316)

  JVS(1032) = -B(269)

  JVS(1033) = -B(309)+0.5*B(311)

  JVS(1034) = -0.397*B(357)+0.5*B(385)

  JVS(1035) = -0.34*B(250)+B(252)

  JVS(1036) = -B(275)

  JVS(1037) = -B(265)

  JVS(1038) = -B(260)

  JVS(1039) = -B(49)

  JVS(1040) = -B(46)+B(48)

  JVS(1041) = -B(321)+0.12*B(323)

  JVS(1042) = -B(237)

  JVS(1043) = -B(387)+0.32*B(389)

  JVS(1044) = 0

  JVS(1045) = -B(329)+0.266*B(331)

  JVS(1046) = -B(369)+0.155*B(371)

  JVS(1047) = -B(280)+0.208*B(282)+0.33*B(288)

  JVS(1048) = -B(345)+0.567*B(347)

  JVS(1049) = -B(337)+0.567*B(339)

  JVS(1050) = -B(255)

  JVS(1051) = -B(296)+0.285*B(298)

  JVS(1052) = -B(227)

  JVS(1053) = -B(377)+0.378*B(379)

  JVS(1054) = -B(289)+0.164*B(291)

  JVS(1055) = -B(218)

  JVS(1056) = -B(306)

  JVS(1057) = 0

  JVS(1058) = -B(232)

  JVS(1059) = -B(240)

  JVS(1060) = -B(303)

  JVS(1061) = -B(51)+B(61)+0.208*B(283)+0.164*B(292)+0.285*B(299)+0.5*B(312)+0.12*B(324)+0.266*B(332)+0.567*B(340)+0.567&
                &*B(348)+0.155*B(372)+0.378*B(380)+0.5*B(386)+0.32*B(390)

  JVS(1062) = 0

  JVS(1063) = 0

  JVS(1064) = 0

  JVS(1065) = -B(42)

  JVS(1066) = -B(44)+0.8*B(66)

  JVS(1067) = 0

  JVS(1068) = 0

  JVS(1069) = 0

  JVS(1070) = B(53)+B(62)+0.8*B(67)-B(72)

  JVS(1071) = -B(36)-B(41)-B(43)-B(45)-B(47)-B(50)-B(52)-B(60)-B(71)-B(73)-B(75)-B(76)-B(219)-B(228)-B(233)-B(238)&
                &-B(241)-B(244)-B(246)-0.65*B(248)-0.34*B(251)-B(256)-B(261)-B(266)-B(270)-B(276)-B(281)-B(290)-B(297)&
                &-B(304)-B(307)-B(310)-B(314)-B(317)-B(320)-B(322)-B(330)-B(338)-B(346)-B(354)-B(356)-0.397*B(358)-B(360)&
                &-B(362)-B(364)-B(366)-B(368)-B(370)-B(378)-B(388)-B(421)-B(425)-B(429)-B(433)-B(437)-B(441)

  JVS(1072) = 0

  JVS(1073) = -B(37)+B(54)

  JVS(1074) = B(139)

  JVS(1075) = 0.1*B(282)

  JVS(1076) = 0.201*B(347)

  JVS(1077) = 0.201*B(339)

  JVS(1078) = 0.37*B(255)+0.37*B(257)

  JVS(1079) = 0.048*B(298)+0.3*B(302)

  JVS(1080) = 0.006*B(379)

  JVS(1081) = 0.05*B(291)

  JVS(1082) = 0

  JVS(1083) = 0.965*B(232)+B(235)

  JVS(1084) = 0.096*B(240)

  JVS(1085) = 0.049*B(303)+0.333*B(305)

  JVS(1086) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.201*B(340)+0.201*B(348)+0.006*B(380)

  JVS(1087) = -B(154)

  JVS(1088) = -B(148)

  JVS(1089) = -B(152)

  JVS(1090) = -B(137)

  JVS(1091) = -B(144)+B(236)+0.37*B(258)

  JVS(1092) = -B(146)

  JVS(1093) = -B(176)

  JVS(1094) = -B(198)

  JVS(1095) = -B(142)

  JVS(1096) = 0.965*B(233)+0.096*B(241)+0.37*B(256)+0.049*B(304)

  JVS(1097) = -B(138)-B(140)-B(143)-B(145)-B(147)-B(149)-B(153)-B(155)-2*B(156)-B(177)-B(199)

  JVS(1098) = -B(141)

  JVS(1099) = B(38)

  JVS(1100) = -B(223)

  JVS(1101) = -B(95)

  JVS(1102) = 0

  JVS(1103) = 0

  JVS(1104) = 0

  JVS(1105) = 0

  JVS(1106) = 0

  JVS(1107) = 0

  JVS(1108) = -B(6)+B(9)

  JVS(1109) = 0

  JVS(1110) = 0

  JVS(1111) = -B(13)

  JVS(1112) = -B(122)

  JVS(1113) = -B(86)

  JVS(1114) = -B(106)

  JVS(1115) = B(1)+B(10)+B(26)

  JVS(1116) = -B(17)+B(27)+B(28)

  JVS(1117) = -B(78)

  JVS(1118) = -B(160)

  JVS(1119) = -B(182)

  JVS(1120) = -B(53)

  JVS(1121) = -B(36)

  JVS(1122) = -B(140)

  JVS(1123) = -B(7)-B(14)-B(18)-2*B(19)-B(37)-B(54)-B(79)-B(87)-B(96)-B(107)-B(123)-B(141)-B(161)-B(183)-B(224)
      
END SUBROUTINE saprc99_mosaic_4bin_vbs2_Jac_SP














SUBROUTINE saprc99_mosaic_4bin_vbs2_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(1123), W(100), a
      INTEGER  :: k, kk, j, jj

      a = 0. 
      IER = 0
      DO k=1,NVAR

        
        IF ( ABS(JVS(LU_DIAG(k))) < TINY(a) ) THEN
            IER = k
            RETURN
        END IF
        DO kk = LU_CROW(k), LU_CROW(k+1)-1
              W( LU_ICOL(kk) ) = JVS(kk)
        END DO
        DO kk = LU_CROW(k), LU_DIAG(k)-1
            j = LU_ICOL(kk)
            a = -W(j) / JVS( LU_DIAG(j) )
            W(j) = -a
            DO jj = LU_DIAG(j)+1, LU_CROW(j+1)-1
               W( LU_ICOL(jj) ) = W( LU_ICOL(jj) ) + a*JVS(jj)
            END DO
         END DO
         DO kk = LU_CROW(k), LU_CROW(k+1)-1
            JVS(kk) = W( LU_ICOL(kk) )
         END DO
      END DO
      
END SUBROUTINE saprc99_mosaic_4bin_vbs2_KppDecomp



SUBROUTINE saprc99_mosaic_4bin_vbs2_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(1123), W(100), a
      INTEGER  :: k, kk, j, jj

      IER = 0
      DO k=1,NVAR
        IF ( JVS( LU_DIAG(k) ) .EQ. 0. ) THEN
            IER = k
            RETURN
        END IF
        DO kk = LU_CROW(k), LU_CROW(k+1)-1
              W( LU_ICOL(kk) ) = JVS(kk)
        END DO
        DO kk = LU_CROW(k), LU_DIAG(k)-1
            j = LU_ICOL(kk)
            a = -W(j) / JVS( LU_DIAG(j) )
            W(j) = -a
            DO jj = LU_DIAG(j)+1, LU_CROW(j+1)-1
               W( LU_ICOL(jj) ) = W( LU_ICOL(jj) ) + a*JVS(jj)
            END DO
         END DO
         DO kk = LU_CROW(k), LU_CROW(k+1)-1
            JVS(kk) = W( LU_ICOL(kk) )
         END DO
      END DO
      
END SUBROUTINE saprc99_mosaic_4bin_vbs2_KppDecompCmplx


SUBROUTINE saprc99_mosaic_4bin_vbs2_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(1123), X(100), sum

      DO i=1,NVAR
         DO j = LU_CROW(i), LU_DIAG(i)-1 
             X(i) = X(i) - JVS(j)*X(LU_ICOL(j));
         END DO  
      END DO

      DO i=NVAR,1,-1
        sum = X(i);
        DO j = LU_DIAG(i)+1, LU_CROW(i+1)-1
          sum = sum - JVS(j)*X(LU_ICOL(j));
        END DO
        X(i) = sum/JVS(LU_DIAG(i));
      END DO
      
END SUBROUTINE saprc99_mosaic_4bin_vbs2_KppSolveIndirect


SUBROUTINE saprc99_mosaic_4bin_vbs2_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(1123), X(100), sum

      DO i=1,NVAR
         DO j = LU_CROW(i), LU_DIAG(i)-1 
             X(i) = X(i) - JVS(j)*X(LU_ICOL(j));
         END DO  
      END DO

      DO i=NVAR,1,-1
        sum = X(i);
        DO j = LU_DIAG(i)+1, LU_CROW(i+1)-1
          sum = sum - JVS(j)*X(LU_ICOL(j));
        END DO
        X(i) = sum/JVS(LU_DIAG(i));
      END DO
      
END SUBROUTINE saprc99_mosaic_4bin_vbs2_KppSolveCmplx













SUBROUTINE saprc99_mosaic_4bin_vbs2_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(21) = X(21)-JVS(154)*X(3)-JVS(155)*X(4)-JVS(156)*X(5)-JVS(157)*X(6)
  X(28) = X(28)-JVS(184)*X(27)
  X(30) = X(30)-JVS(189)*X(29)
  X(45) = X(45)-JVS(231)*X(44)
  X(55) = X(55)-JVS(265)*X(50)-JVS(266)*X(54)
  X(56) = X(56)-JVS(269)*X(50)-JVS(270)*X(54)
  X(57) = X(57)-JVS(273)*X(50)-JVS(274)*X(54)
  X(58) = X(58)-JVS(278)*X(50)-JVS(279)*X(54)
  X(61) = X(61)-JVS(291)*X(49)
  X(62) = X(62)-JVS(297)*X(50)-JVS(298)*X(54)
  X(63) = X(63)-JVS(304)*X(54)
  X(64) = X(64)-JVS(310)*X(50)-JVS(311)*X(54)-JVS(312)*X(55)-JVS(313)*X(56)-JVS(314)*X(57)
  X(65) = X(65)-JVS(322)*X(53)-JVS(323)*X(55)-JVS(324)*X(56)-JVS(325)*X(58)-JVS(326)*X(59)-JVS(327)*X(64)
  X(66) = X(66)-JVS(347)*X(42)-JVS(348)*X(57)-JVS(349)*X(61)-JVS(350)*X(62)-JVS(351)*X(63)-JVS(352)*X(64)
  X(68) = X(68)-JVS(373)*X(39)-JVS(374)*X(44)-JVS(375)*X(45)-JVS(376)*X(46)-JVS(377)*X(53)
  X(70) = X(70)-JVS(393)*X(57)-JVS(394)*X(63)
  X(73) = X(73)-JVS(415)*X(44)-JVS(416)*X(46)-JVS(417)*X(53)-JVS(418)*X(55)-JVS(419)*X(56)-JVS(420)*X(68)-JVS(421)*X(71)&
            &-JVS(422)*X(72)
  X(74) = X(74)-JVS(445)*X(71)
  X(77) = X(77)-JVS(462)*X(50)-JVS(463)*X(54)-JVS(464)*X(55)-JVS(465)*X(56)-JVS(466)*X(58)-JVS(467)*X(59)-JVS(468)*X(63)&
            &-JVS(469)*X(67)-JVS(470)*X(70)-JVS(471)*X(75)-JVS(472)*X(76)
  X(78) = X(78)-JVS(486)*X(71)
  X(79) = X(79)-JVS(493)*X(33)-JVS(494)*X(41)-JVS(495)*X(44)-JVS(496)*X(46)-JVS(497)*X(53)-JVS(498)*X(67)-JVS(499)*X(69)&
            &-JVS(500)*X(72)-JVS(501)*X(78)
  X(81) = X(81)-JVS(522)*X(71)-JVS(523)*X(80)
  X(82) = X(82)-JVS(529)*X(41)-JVS(530)*X(44)-JVS(531)*X(46)-JVS(532)*X(47)-JVS(533)*X(48)-JVS(534)*X(52)-JVS(535)*X(53)&
            &-JVS(536)*X(59)-JVS(537)*X(67)-JVS(538)*X(68)-JVS(539)*X(69)-JVS(540)*X(71)-JVS(541)*X(72)-JVS(542)*X(74)&
            &-JVS(543)*X(75)-JVS(544)*X(76)-JVS(545)*X(77)-JVS(546)*X(78)-JVS(547)*X(80)-JVS(548)*X(81)
  X(83) = X(83)-JVS(567)*X(45)-JVS(568)*X(72)-JVS(569)*X(75)-JVS(570)*X(76)-JVS(571)*X(78)-JVS(572)*X(80)
  X(84) = X(84)-JVS(581)*X(26)-JVS(582)*X(67)-JVS(583)*X(69)-JVS(584)*X(71)-JVS(585)*X(72)-JVS(586)*X(74)-JVS(587)*X(75)&
            &-JVS(588)*X(76)-JVS(589)*X(80)-JVS(590)*X(81)
  X(85) = X(85)-JVS(597)*X(39)-JVS(598)*X(44)-JVS(599)*X(46)-JVS(600)*X(53)-JVS(601)*X(55)-JVS(602)*X(56)-JVS(603)*X(58)&
            &-JVS(604)*X(60)-JVS(605)*X(67)-JVS(606)*X(69)-JVS(607)*X(72)-JVS(608)*X(74)-JVS(609)*X(75)-JVS(610)*X(76)&
            &-JVS(611)*X(78)-JVS(612)*X(80)-JVS(613)*X(81)-JVS(614)*X(83)-JVS(615)*X(84)
  X(86) = X(86)-JVS(627)*X(44)-JVS(628)*X(46)-JVS(629)*X(53)-JVS(630)*X(69)-JVS(631)*X(72)-JVS(632)*X(74)-JVS(633)*X(78)&
            &-JVS(634)*X(80)-JVS(635)*X(81)-JVS(636)*X(83)-JVS(637)*X(84)
  X(87) = X(87)-JVS(648)*X(46)-JVS(649)*X(53)-JVS(650)*X(54)-JVS(651)*X(71)-JVS(652)*X(72)-JVS(653)*X(75)-JVS(654)*X(76)&
            &-JVS(655)*X(78)-JVS(656)*X(80)-JVS(657)*X(81)-JVS(658)*X(83)-JVS(659)*X(84)
  X(88) = X(88)-JVS(672)*X(58)-JVS(673)*X(59)-JVS(674)*X(67)-JVS(675)*X(69)-JVS(676)*X(71)-JVS(677)*X(72)-JVS(678)*X(74)&
            &-JVS(679)*X(75)-JVS(680)*X(76)-JVS(681)*X(78)-JVS(682)*X(80)-JVS(683)*X(81)-JVS(684)*X(84)
  X(89) = X(89)-JVS(695)*X(34)-JVS(696)*X(40)-JVS(697)*X(50)-JVS(698)*X(53)-JVS(699)*X(55)-JVS(700)*X(56)-JVS(701)*X(64)&
            &-JVS(702)*X(68)-JVS(703)*X(72)-JVS(704)*X(74)-JVS(705)*X(75)-JVS(706)*X(76)-JVS(707)*X(78)-JVS(708)*X(79)&
            &-JVS(709)*X(80)-JVS(710)*X(81)-JVS(711)*X(83)-JVS(712)*X(84)-JVS(713)*X(85)-JVS(714)*X(86)-JVS(715)*X(87)&
            &-JVS(716)*X(88)
  X(90) = X(90)-JVS(729)*X(33)-JVS(730)*X(39)-JVS(731)*X(41)-JVS(732)*X(44)-JVS(733)*X(46)-JVS(734)*X(50)-JVS(735)*X(53)&
            &-JVS(736)*X(54)-JVS(737)*X(55)-JVS(738)*X(56)-JVS(739)*X(57)-JVS(740)*X(58)-JVS(741)*X(59)-JVS(742)*X(60)&
            &-JVS(743)*X(63)-JVS(744)*X(67)-JVS(745)*X(69)-JVS(746)*X(70)-JVS(747)*X(71)-JVS(748)*X(72)-JVS(749)*X(74)&
            &-JVS(750)*X(75)-JVS(751)*X(76)-JVS(752)*X(78)-JVS(753)*X(80)-JVS(754)*X(81)-JVS(755)*X(83)-JVS(756)*X(84)&
            &-JVS(757)*X(85)-JVS(758)*X(86)-JVS(759)*X(87)-JVS(760)*X(88)-JVS(761)*X(89)
  X(91) = X(91)-JVS(773)*X(39)-JVS(774)*X(44)-JVS(775)*X(46)-JVS(776)*X(50)-JVS(777)*X(53)-JVS(778)*X(54)-JVS(779)*X(69)&
            &-JVS(780)*X(71)-JVS(781)*X(72)-JVS(782)*X(75)-JVS(783)*X(76)-JVS(784)*X(78)-JVS(785)*X(80)-JVS(786)*X(81)&
            &-JVS(787)*X(83)-JVS(788)*X(84)-JVS(789)*X(85)-JVS(790)*X(86)-JVS(791)*X(87)-JVS(792)*X(88)-JVS(793)*X(89)&
            &-JVS(794)*X(90)
  X(92) = X(92)-JVS(805)*X(34)-JVS(806)*X(35)-JVS(807)*X(36)-JVS(808)*X(37)-JVS(809)*X(42)-JVS(810)*X(43)-JVS(811)*X(45)&
            &-JVS(812)*X(48)-JVS(813)*X(49)-JVS(814)*X(51)-JVS(815)*X(61)-JVS(816)*X(66)-JVS(817)*X(70)-JVS(818)*X(71)&
            &-JVS(819)*X(73)-JVS(820)*X(74)-JVS(821)*X(75)-JVS(822)*X(76)-JVS(823)*X(77)-JVS(824)*X(78)-JVS(825)*X(79)&
            &-JVS(826)*X(80)-JVS(827)*X(81)-JVS(828)*X(82)-JVS(829)*X(83)-JVS(830)*X(84)-JVS(831)*X(85)-JVS(832)*X(86)&
            &-JVS(833)*X(87)-JVS(834)*X(88)-JVS(835)*X(89)-JVS(836)*X(90)-JVS(837)*X(91)
  X(93) = X(93)-JVS(847)*X(42)-JVS(848)*X(51)-JVS(849)*X(57)-JVS(850)*X(61)-JVS(851)*X(62)-JVS(852)*X(63)-JVS(853)*X(64)&
            &-JVS(854)*X(66)-JVS(855)*X(67)-JVS(856)*X(69)-JVS(857)*X(70)-JVS(858)*X(71)-JVS(859)*X(72)-JVS(860)*X(73)&
            &-JVS(861)*X(74)-JVS(862)*X(75)-JVS(863)*X(76)-JVS(864)*X(77)-JVS(865)*X(78)-JVS(866)*X(79)-JVS(867)*X(80)&
            &-JVS(868)*X(81)-JVS(869)*X(82)-JVS(870)*X(83)-JVS(871)*X(84)-JVS(872)*X(85)-JVS(873)*X(86)-JVS(874)*X(87)&
            &-JVS(875)*X(88)-JVS(876)*X(89)-JVS(877)*X(90)-JVS(878)*X(91)-JVS(879)*X(92)
  X(94) = X(94)-JVS(888)*X(31)-JVS(889)*X(45)-JVS(890)*X(47)-JVS(891)*X(53)-JVS(892)*X(67)-JVS(893)*X(68)-JVS(894)*X(69)&
            &-JVS(895)*X(71)-JVS(896)*X(72)-JVS(897)*X(75)-JVS(898)*X(76)-JVS(899)*X(79)-JVS(900)*X(80)-JVS(901)*X(81)&
            &-JVS(902)*X(83)-JVS(903)*X(84)-JVS(904)*X(85)-JVS(905)*X(86)-JVS(906)*X(87)-JVS(907)*X(88)-JVS(908)*X(89)&
            &-JVS(909)*X(90)-JVS(910)*X(91)-JVS(911)*X(92)-JVS(912)*X(93)
  X(95) = X(95)-JVS(920)*X(36)-JVS(921)*X(62)-JVS(922)*X(80)-JVS(923)*X(84)-JVS(924)*X(88)-JVS(925)*X(89)-JVS(926)*X(90)&
            &-JVS(927)*X(91)-JVS(928)*X(92)-JVS(929)*X(93)-JVS(930)*X(94)
  X(96) = X(96)-JVS(937)*X(37)-JVS(938)*X(71)-JVS(939)*X(74)-JVS(940)*X(78)-JVS(941)*X(80)-JVS(942)*X(81)-JVS(943)*X(84)&
            &-JVS(944)*X(88)-JVS(945)*X(89)-JVS(946)*X(90)-JVS(947)*X(91)-JVS(948)*X(92)-JVS(949)*X(93)-JVS(950)*X(94)&
            &-JVS(951)*X(95)
  X(97) = X(97)-JVS(957)*X(32)-JVS(958)*X(38)-JVS(959)*X(41)-JVS(960)*X(43)-JVS(961)*X(47)-JVS(962)*X(48)-JVS(963)*X(49)&
            &-JVS(964)*X(50)-JVS(965)*X(51)-JVS(966)*X(52)-JVS(967)*X(54)-JVS(968)*X(55)-JVS(969)*X(56)-JVS(970)*X(58)&
            &-JVS(971)*X(59)-JVS(972)*X(60)-JVS(973)*X(61)-JVS(974)*X(64)-JVS(975)*X(65)-JVS(976)*X(67)-JVS(977)*X(69)&
            &-JVS(978)*X(70)-JVS(979)*X(71)-JVS(980)*X(72)-JVS(981)*X(74)-JVS(982)*X(75)-JVS(983)*X(76)-JVS(984)*X(77)&
            &-JVS(985)*X(78)-JVS(986)*X(79)-JVS(987)*X(80)-JVS(988)*X(81)-JVS(989)*X(82)-JVS(990)*X(83)-JVS(991)*X(84)&
            &-JVS(992)*X(85)-JVS(993)*X(86)-JVS(994)*X(87)-JVS(995)*X(88)-JVS(996)*X(89)-JVS(997)*X(90)-JVS(998)*X(91)&
            &-JVS(999)*X(92)-JVS(1000)*X(93)-JVS(1001)*X(94)-JVS(1002)*X(95)-JVS(1003)*X(96)
  X(98) = X(98)-JVS(1008)*X(3)-JVS(1009)*X(5)-JVS(1010)*X(26)-JVS(1011)*X(27)-JVS(1012)*X(28)-JVS(1013)*X(29)-JVS(1014)&
            &*X(30)-JVS(1015)*X(31)-JVS(1016)*X(32)-JVS(1017)*X(33)-JVS(1018)*X(38)-JVS(1019)*X(39)-JVS(1020)*X(41)&
            &-JVS(1021)*X(43)-JVS(1022)*X(44)-JVS(1023)*X(46)-JVS(1024)*X(47)-JVS(1025)*X(50)-JVS(1026)*X(51)-JVS(1027)&
            &*X(52)-JVS(1028)*X(53)-JVS(1029)*X(54)-JVS(1030)*X(55)-JVS(1031)*X(56)-JVS(1032)*X(57)-JVS(1033)*X(58)&
            &-JVS(1034)*X(59)-JVS(1035)*X(60)-JVS(1036)*X(62)-JVS(1037)*X(63)-JVS(1038)*X(64)-JVS(1039)*X(65)-JVS(1040)&
            &*X(66)-JVS(1041)*X(67)-JVS(1042)*X(68)-JVS(1043)*X(69)-JVS(1044)*X(70)-JVS(1045)*X(71)-JVS(1046)*X(72)&
            &-JVS(1047)*X(74)-JVS(1048)*X(75)-JVS(1049)*X(76)-JVS(1050)*X(77)-JVS(1051)*X(78)-JVS(1052)*X(79)-JVS(1053)&
            &*X(80)-JVS(1054)*X(81)-JVS(1055)*X(82)-JVS(1056)*X(83)-JVS(1057)*X(84)-JVS(1058)*X(85)-JVS(1059)*X(86)&
            &-JVS(1060)*X(87)-JVS(1061)*X(88)-JVS(1062)*X(89)-JVS(1063)*X(90)-JVS(1064)*X(91)-JVS(1065)*X(92)-JVS(1066)&
            &*X(93)-JVS(1067)*X(94)-JVS(1068)*X(95)-JVS(1069)*X(96)-JVS(1070)*X(97)
  X(99) = X(99)-JVS(1074)*X(35)-JVS(1075)*X(74)-JVS(1076)*X(75)-JVS(1077)*X(76)-JVS(1078)*X(77)-JVS(1079)*X(78)&
            &-JVS(1080)*X(80)-JVS(1081)*X(81)-JVS(1082)*X(84)-JVS(1083)*X(85)-JVS(1084)*X(86)-JVS(1085)*X(87)-JVS(1086)&
            &*X(88)-JVS(1087)*X(89)-JVS(1088)*X(90)-JVS(1089)*X(91)-JVS(1090)*X(92)-JVS(1091)*X(93)-JVS(1092)*X(94)&
            &-JVS(1093)*X(95)-JVS(1094)*X(96)-JVS(1095)*X(97)-JVS(1096)*X(98)
  X(100) = X(100)-JVS(1099)*X(43)-JVS(1100)*X(48)-JVS(1101)*X(73)-JVS(1102)*X(75)-JVS(1103)*X(76)-JVS(1104)*X(80)&
             &-JVS(1105)*X(81)-JVS(1106)*X(82)-JVS(1107)*X(83)-JVS(1108)*X(84)-JVS(1109)*X(86)-JVS(1110)*X(87)-JVS(1111)&
             &*X(88)-JVS(1112)*X(89)-JVS(1113)*X(90)-JVS(1114)*X(91)-JVS(1115)*X(92)-JVS(1116)*X(93)-JVS(1117)*X(94)&
             &-JVS(1118)*X(95)-JVS(1119)*X(96)-JVS(1120)*X(97)-JVS(1121)*X(98)-JVS(1122)*X(99)
  X(100) = X(100)/JVS(1123)
  X(99) = (X(99)-JVS(1098)*X(100))/(JVS(1097))
  X(98) = (X(98)-JVS(1072)*X(99)-JVS(1073)*X(100))/(JVS(1071))
  X(97) = (X(97)-JVS(1005)*X(98)-JVS(1006)*X(99)-JVS(1007)*X(100))/(JVS(1004))
  X(96) = (X(96)-JVS(953)*X(97)-JVS(954)*X(98)-JVS(955)*X(99)-JVS(956)*X(100))/(JVS(952))
  X(95) = (X(95)-JVS(932)*X(96)-JVS(933)*X(97)-JVS(934)*X(98)-JVS(935)*X(99)-JVS(936)*X(100))/(JVS(931))
  X(94) = (X(94)-JVS(914)*X(95)-JVS(915)*X(96)-JVS(916)*X(97)-JVS(917)*X(98)-JVS(918)*X(99)-JVS(919)*X(100))/(JVS(913))
  X(93) = (X(93)-JVS(881)*X(94)-JVS(882)*X(95)-JVS(883)*X(96)-JVS(884)*X(97)-JVS(885)*X(98)-JVS(886)*X(99)-JVS(887)&
            &*X(100))/(JVS(880))
  X(92) = (X(92)-JVS(839)*X(93)-JVS(840)*X(94)-JVS(841)*X(95)-JVS(842)*X(96)-JVS(843)*X(97)-JVS(844)*X(98)-JVS(845)&
            &*X(99)-JVS(846)*X(100))/(JVS(838))
  X(91) = (X(91)-JVS(796)*X(92)-JVS(797)*X(93)-JVS(798)*X(94)-JVS(799)*X(95)-JVS(800)*X(96)-JVS(801)*X(97)-JVS(802)&
            &*X(98)-JVS(803)*X(99)-JVS(804)*X(100))/(JVS(795))
  X(90) = (X(90)-JVS(763)*X(91)-JVS(764)*X(92)-JVS(765)*X(93)-JVS(766)*X(94)-JVS(767)*X(95)-JVS(768)*X(96)-JVS(769)&
            &*X(97)-JVS(770)*X(98)-JVS(771)*X(99)-JVS(772)*X(100))/(JVS(762))
  X(89) = (X(89)-JVS(718)*X(90)-JVS(719)*X(91)-JVS(720)*X(92)-JVS(721)*X(93)-JVS(722)*X(94)-JVS(723)*X(95)-JVS(724)&
            &*X(96)-JVS(725)*X(97)-JVS(726)*X(98)-JVS(727)*X(99)-JVS(728)*X(100))/(JVS(717))
  X(88) = (X(88)-JVS(686)*X(89)-JVS(687)*X(92)-JVS(688)*X(93)-JVS(689)*X(95)-JVS(690)*X(96)-JVS(691)*X(97)-JVS(692)&
            &*X(98)-JVS(693)*X(99)-JVS(694)*X(100))/(JVS(685))
  X(87) = (X(87)-JVS(661)*X(88)-JVS(662)*X(89)-JVS(663)*X(90)-JVS(664)*X(91)-JVS(665)*X(92)-JVS(666)*X(93)-JVS(667)&
            &*X(94)-JVS(668)*X(95)-JVS(669)*X(98)-JVS(670)*X(99)-JVS(671)*X(100))/(JVS(660))
  X(86) = (X(86)-JVS(639)*X(87)-JVS(640)*X(88)-JVS(641)*X(90)-JVS(642)*X(91)-JVS(643)*X(92)-JVS(644)*X(93)-JVS(645)&
            &*X(94)-JVS(646)*X(98)-JVS(647)*X(100))/(JVS(638))
  X(85) = (X(85)-JVS(617)*X(86)-JVS(618)*X(87)-JVS(619)*X(88)-JVS(620)*X(90)-JVS(621)*X(91)-JVS(622)*X(92)-JVS(623)&
            &*X(93)-JVS(624)*X(97)-JVS(625)*X(98)-JVS(626)*X(100))/(JVS(616))
  X(84) = (X(84)-JVS(592)*X(88)-JVS(593)*X(92)-JVS(594)*X(93)-JVS(595)*X(98)-JVS(596)*X(100))/(JVS(591))
  X(83) = (X(83)-JVS(574)*X(84)-JVS(575)*X(88)-JVS(576)*X(91)-JVS(577)*X(92)-JVS(578)*X(93)-JVS(579)*X(98)-JVS(580)&
            &*X(100))/(JVS(573))
  X(82) = (X(82)-JVS(550)*X(83)-JVS(551)*X(84)-JVS(552)*X(86)-JVS(553)*X(87)-JVS(554)*X(88)-JVS(555)*X(89)-JVS(556)&
            &*X(90)-JVS(557)*X(91)-JVS(558)*X(92)-JVS(559)*X(93)-JVS(560)*X(94)-JVS(561)*X(95)-JVS(562)*X(96)-JVS(563)*X(97)&
            &-JVS(564)*X(98)-JVS(565)*X(99)-JVS(566)*X(100))/(JVS(549))
  X(81) = (X(81)-JVS(525)*X(84)-JVS(526)*X(88)-JVS(527)*X(93)-JVS(528)*X(98))/(JVS(524))
  X(80) = (X(80)-JVS(518)*X(84)-JVS(519)*X(88)-JVS(520)*X(93)-JVS(521)*X(98))/(JVS(517))
  X(79) = (X(79)-JVS(503)*X(80)-JVS(504)*X(83)-JVS(505)*X(84)-JVS(506)*X(85)-JVS(507)*X(86)-JVS(508)*X(87)-JVS(509)&
            &*X(88)-JVS(510)*X(89)-JVS(511)*X(93)-JVS(512)*X(95)-JVS(513)*X(96)-JVS(514)*X(98)-JVS(515)*X(99)-JVS(516)&
            &*X(100))/(JVS(502))
  X(78) = (X(78)-JVS(488)*X(80)-JVS(489)*X(84)-JVS(490)*X(88)-JVS(491)*X(93)-JVS(492)*X(98))/(JVS(487))
  X(77) = (X(77)-JVS(474)*X(78)-JVS(475)*X(84)-JVS(476)*X(88)-JVS(477)*X(89)-JVS(478)*X(92)-JVS(479)*X(93)-JVS(480)&
            &*X(95)-JVS(481)*X(96)-JVS(482)*X(97)-JVS(483)*X(98)-JVS(484)*X(99)-JVS(485)*X(100))/(JVS(473))
  X(76) = (X(76)-JVS(458)*X(84)-JVS(459)*X(88)-JVS(460)*X(93)-JVS(461)*X(98))/(JVS(457))
  X(75) = (X(75)-JVS(453)*X(84)-JVS(454)*X(88)-JVS(455)*X(93)-JVS(456)*X(98))/(JVS(452))
  X(74) = (X(74)-JVS(447)*X(80)-JVS(448)*X(84)-JVS(449)*X(88)-JVS(450)*X(93)-JVS(451)*X(98))/(JVS(446))
  X(73) = (X(73)-JVS(424)*X(75)-JVS(425)*X(76)-JVS(426)*X(80)-JVS(427)*X(81)-JVS(428)*X(83)-JVS(429)*X(84)-JVS(430)&
            &*X(86)-JVS(431)*X(87)-JVS(432)*X(88)-JVS(433)*X(89)-JVS(434)*X(90)-JVS(435)*X(91)-JVS(436)*X(92)-JVS(437)*X(93)&
            &-JVS(438)*X(94)-JVS(439)*X(95)-JVS(440)*X(96)-JVS(441)*X(97)-JVS(442)*X(98)-JVS(443)*X(99)-JVS(444)*X(100))&
            &/(JVS(423))
  X(72) = (X(72)-JVS(411)*X(84)-JVS(412)*X(88)-JVS(413)*X(93)-JVS(414)*X(98))/(JVS(410))
  X(71) = (X(71)-JVS(406)*X(84)-JVS(407)*X(88)-JVS(408)*X(93)-JVS(409)*X(98))/(JVS(405))
  X(70) = (X(70)-JVS(396)*X(89)-JVS(397)*X(92)-JVS(398)*X(93)-JVS(399)*X(95)-JVS(400)*X(96)-JVS(401)*X(97)-JVS(402)&
            &*X(98)-JVS(403)*X(99)-JVS(404)*X(100))/(JVS(395))
  X(69) = (X(69)-JVS(389)*X(84)-JVS(390)*X(88)-JVS(391)*X(93)-JVS(392)*X(98))/(JVS(388))
  X(68) = (X(68)-JVS(379)*X(72)-JVS(380)*X(75)-JVS(381)*X(76)-JVS(382)*X(80)-JVS(383)*X(83)-JVS(384)*X(88)-JVS(385)&
            &*X(92)-JVS(386)*X(93)-JVS(387)*X(98))/(JVS(378))
  X(67) = (X(67)-JVS(369)*X(84)-JVS(370)*X(88)-JVS(371)*X(93)-JVS(372)*X(98))/(JVS(368))
  X(66) = (X(66)-JVS(354)*X(70)-JVS(355)*X(74)-JVS(356)*X(77)-JVS(357)*X(78)-JVS(358)*X(79)-JVS(359)*X(80)-JVS(360)&
            &*X(81)-JVS(361)*X(82)-JVS(362)*X(85)-JVS(363)*X(88)-JVS(364)*X(92)-JVS(365)*X(93)-JVS(366)*X(97)-JVS(367)&
            &*X(98))/(JVS(353))
  X(65) = (X(65)-JVS(329)*X(67)-JVS(330)*X(69)-JVS(331)*X(71)-JVS(332)*X(72)-JVS(333)*X(74)-JVS(334)*X(75)-JVS(335)&
            &*X(76)-JVS(336)*X(77)-JVS(337)*X(78)-JVS(338)*X(79)-JVS(339)*X(80)-JVS(340)*X(81)-JVS(341)*X(82)-JVS(342)*X(84)&
            &-JVS(343)*X(85)-JVS(344)*X(88)-JVS(345)*X(93)-JVS(346)*X(98))/(JVS(328))
  X(64) = (X(64)-JVS(316)*X(74)-JVS(317)*X(78)-JVS(318)*X(81)-JVS(319)*X(88)-JVS(320)*X(93)-JVS(321)*X(98))/(JVS(315))
  X(63) = (X(63)-JVS(306)*X(70)-JVS(307)*X(93)-JVS(308)*X(97)-JVS(309)*X(98))/(JVS(305))
  X(62) = (X(62)-JVS(300)*X(80)-JVS(301)*X(88)-JVS(302)*X(93)-JVS(303)*X(98))/(JVS(299))
  X(61) = (X(61)-JVS(293)*X(70)-JVS(294)*X(92)-JVS(295)*X(93)-JVS(296)*X(97))/(JVS(292))
  X(60) = (X(60)-JVS(287)*X(90)-JVS(288)*X(91)-JVS(289)*X(97)-JVS(290)*X(98))/(JVS(286))
  X(59) = (X(59)-JVS(284)*X(88)-JVS(285)*X(98))/(JVS(283))
  X(58) = (X(58)-JVS(281)*X(88)-JVS(282)*X(98))/(JVS(280))
  X(57) = (X(57)-JVS(276)*X(93)-JVS(277)*X(98))/(JVS(275))
  X(56) = (X(56)-JVS(272)*X(98))/(JVS(271))
  X(55) = (X(55)-JVS(268)*X(98))/(JVS(267))
  X(54) = (X(54)-JVS(264)*X(98))/(JVS(263))
  X(53) = (X(53)-JVS(262)*X(98))/(JVS(261))
  X(52) = (X(52)-JVS(257)*X(90)-JVS(258)*X(91)-JVS(259)*X(94)-JVS(260)*X(98))/(JVS(256))
  X(51) = (X(51)-JVS(253)*X(92)-JVS(254)*X(97)-JVS(255)*X(98))/(JVS(252))
  X(50) = (X(50)-JVS(251)*X(98))/(JVS(250))
  X(49) = (X(49)-JVS(246)*X(61)-JVS(247)*X(92)-JVS(248)*X(93)-JVS(249)*X(97))/(JVS(245))
  X(48) = (X(48)-JVS(242)*X(82)-JVS(243)*X(97)-JVS(244)*X(100))/(JVS(241))
  X(47) = (X(47)-JVS(238)*X(94)-JVS(239)*X(97)-JVS(240)*X(98))/(JVS(237))
  X(46) = (X(46)-JVS(236)*X(98))/(JVS(235))
  X(45) = (X(45)-JVS(233)*X(92)-JVS(234)*X(98))/(JVS(232))
  X(44) = (X(44)-JVS(230)*X(98))/(JVS(229))
  X(43) = (X(43)-JVS(227)*X(98)-JVS(228)*X(100))/(JVS(226))
  X(42) = (X(42)-JVS(224)*X(92)-JVS(225)*X(93))/(JVS(223))
  X(41) = (X(41)-JVS(222)*X(98))/(JVS(221))
  X(40) = (X(40)-JVS(216)*X(50)-JVS(217)*X(75)-JVS(218)*X(76)-JVS(219)*X(88)-JVS(220)*X(98))/(JVS(215))
  X(39) = (X(39)-JVS(214)*X(98))/(JVS(213))
  X(38) = (X(38)-JVS(211)*X(97)-JVS(212)*X(98))/(JVS(210))
  X(37) = (X(37)-JVS(208)*X(92)-JVS(209)*X(96))/(JVS(207))
  X(36) = (X(36)-JVS(205)*X(92)-JVS(206)*X(95))/(JVS(204))
  X(35) = (X(35)-JVS(202)*X(92)-JVS(203)*X(99))/(JVS(201))
  X(34) = (X(34)-JVS(199)*X(89)-JVS(200)*X(92))/(JVS(198))
  X(33) = (X(33)-JVS(197)*X(98))/(JVS(196))
  X(32) = (X(32)-JVS(195)*X(98))/(JVS(194))
  X(31) = (X(31)-JVS(193)*X(98))/(JVS(192))
  X(30) = (X(30)-JVS(191)*X(98))/(JVS(190))
  X(29) = (X(29)-JVS(188)*X(98))/(JVS(187))
  X(28) = (X(28)-JVS(186)*X(98))/(JVS(185))
  X(27) = (X(27)-JVS(183)*X(98))/(JVS(182))
  X(26) = (X(26)-JVS(181)*X(88))/(JVS(180))
  X(25) = (X(25)-JVS(179)*X(98))/(JVS(178))
  X(24) = (X(24)-JVS(175)*X(25)-JVS(176)*X(29)-JVS(177)*X(98))/(JVS(174))
  X(23) = (X(23)-JVS(173)*X(98))/(JVS(172))
  X(22) = (X(22)-JVS(169)*X(23)-JVS(170)*X(27)-JVS(171)*X(98))/(JVS(168))
  X(21) = (X(21)-JVS(159)*X(22)-JVS(160)*X(23)-JVS(161)*X(24)-JVS(162)*X(25)-JVS(163)*X(27)-JVS(164)*X(28)-JVS(165)&
            &*X(29)-JVS(166)*X(30)-JVS(167)*X(98))/(JVS(158))
  X(20) = (X(20)-JVS(106)*X(31)-JVS(107)*X(32)-JVS(108)*X(33)-JVS(109)*X(39)-JVS(110)*X(41)-JVS(111)*X(43)-JVS(112)&
            &*X(44)-JVS(113)*X(46)-JVS(114)*X(47)-JVS(115)*X(50)-JVS(116)*X(52)-JVS(117)*X(53)-JVS(118)*X(54)-JVS(119)*X(55)&
            &-JVS(120)*X(56)-JVS(121)*X(57)-JVS(122)*X(58)-JVS(123)*X(59)-JVS(124)*X(60)-JVS(125)*X(62)-JVS(126)*X(63)&
            &-JVS(127)*X(64)-JVS(128)*X(65)-JVS(129)*X(66)-JVS(130)*X(67)-JVS(131)*X(68)-JVS(132)*X(69)-JVS(133)*X(71)&
            &-JVS(134)*X(72)-JVS(135)*X(74)-JVS(136)*X(75)-JVS(137)*X(76)-JVS(138)*X(77)-JVS(139)*X(78)-JVS(140)*X(79)&
            &-JVS(141)*X(80)-JVS(142)*X(81)-JVS(143)*X(82)-JVS(144)*X(83)-JVS(145)*X(85)-JVS(146)*X(86)-JVS(147)*X(87)&
            &-JVS(148)*X(88)-JVS(149)*X(92)-JVS(150)*X(93)-JVS(151)*X(97)-JVS(152)*X(98)-JVS(153)*X(100))/(JVS(105))
  X(19) = (X(19)-JVS(102)*X(75)-JVS(103)*X(76)-JVS(104)*X(98))/(JVS(101))
  X(18) = (X(18)-JVS(99)*X(71)-JVS(100)*X(98))/(JVS(98))
  X(17) = (X(17)-JVS(95)*X(50)-JVS(96)*X(54)-JVS(97)*X(98))/(JVS(94))
  X(16) = (X(16)-JVS(89)*X(46)-JVS(90)*X(53)-JVS(91)*X(72)-JVS(92)*X(80)-JVS(93)*X(98))/(JVS(88))
  X(15) = (X(15)-JVS(83)*X(73)-JVS(84)*X(90)-JVS(85)*X(91)-JVS(86)*X(94)-JVS(87)*X(97))/(JVS(82))
  X(14) = (X(14)-JVS(76)*X(73)-JVS(77)*X(90)-JVS(78)*X(91)-JVS(79)*X(93)-JVS(80)*X(94)-JVS(81)*X(100))/(JVS(75))
  X(13) = (X(13)-JVS(67)*X(49)-JVS(68)*X(62)-JVS(69)*X(69)-JVS(70)*X(84)-JVS(71)*X(88)-JVS(72)*X(92)-JVS(73)*X(93)&
            &-JVS(74)*X(98))/(JVS(66))
  X(12) = (X(12)-JVS(62)*X(49)-JVS(63)*X(69)-JVS(64)*X(92)-JVS(65)*X(93))/(JVS(61))
  X(11) = (X(11)-JVS(57)*X(95)-JVS(58)*X(96)-JVS(59)*X(97)-JVS(60)*X(99))/(JVS(56))
  X(10) = (X(10)-JVS(54)*X(89)-JVS(55)*X(97))/(JVS(53))
  X(9) = (X(9)-JVS(39)*X(71)-JVS(40)*X(72)-JVS(41)*X(75)-JVS(42)*X(76)-JVS(43)*X(78)-JVS(44)*X(80)-JVS(45)*X(88)-JVS(46)&
           &*X(90)-JVS(47)*X(91)-JVS(48)*X(94)-JVS(49)*X(95)-JVS(50)*X(96)-JVS(51)*X(97)-JVS(52)*X(99))/(JVS(38))
  X(8) = (X(8)-JVS(29)*X(69)-JVS(30)*X(72)-JVS(31)*X(80)-JVS(32)*X(88)-JVS(33)*X(89)-JVS(34)*X(90)-JVS(35)*X(91)-JVS(36)&
           &*X(94)-JVS(37)*X(97))/(JVS(28))
  X(7) = (X(7)-JVS(13)*X(48)-JVS(14)*X(59)-JVS(15)*X(67)-JVS(16)*X(69)-JVS(17)*X(71)-JVS(18)*X(72)-JVS(19)*X(74)-JVS(20)&
           &*X(75)-JVS(21)*X(76)-JVS(22)*X(78)-JVS(23)*X(80)-JVS(24)*X(81)-JVS(25)*X(88)-JVS(26)*X(98)-JVS(27)*X(100))&
           &/(JVS(12))
  X(6) = X(6)/JVS(11)
  X(5) = X(5)/JVS(10)
  X(4) = X(4)/JVS(9)
  X(3) = X(3)/JVS(8)
  X(2) = (X(2)-JVS(5)*X(59)-JVS(6)*X(69)-JVS(7)*X(88))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(32)-JVS(3)*X(98))/(JVS(1))
      
END SUBROUTINE saprc99_mosaic_4bin_vbs2_KppSolve
























      SUBROUTINE saprc99_mosaic_4bin_vbs2_WCOPY(N,X,incX,Y,incY)








      
      INTEGER i,incX,incY,M,MP1,N
      REAL(kind=dp) X(N),Y(N)

      IF (N.LE.0) RETURN

      M = MOD(N,8)
      IF( M .NE. 0 ) THEN
        DO i = 1,M
          Y(i) = X(i)
        END DO
        IF( N .LT. 8 ) RETURN
      END IF    
      MP1 = M+1
      DO i = MP1,N,8
        Y(i) = X(i)
        Y(i + 1) = X(i + 1)
        Y(i + 2) = X(i + 2)
        Y(i + 3) = X(i + 3)
        Y(i + 4) = X(i + 4)
        Y(i + 5) = X(i + 5)
        Y(i + 6) = X(i + 6)
        Y(i + 7) = X(i + 7)
      END DO

      END SUBROUTINE saprc99_mosaic_4bin_vbs2_WCOPY



      SUBROUTINE saprc99_mosaic_4bin_vbs2_WAXPY(N,Alpha,X,incX,Y,incY)









      INTEGER i,incX,incY,M,MP1,N
      REAL(kind=dp) X(N),Y(N),Alpha
      REAL(kind=dp) ZERO
      PARAMETER( ZERO = 0.0_dp )

      IF (Alpha .EQ. ZERO) RETURN
      IF (N .LE. 0) RETURN

      M = MOD(N,4)
      IF( M .NE. 0 ) THEN
        DO i = 1,M
          Y(i) = Y(i) + Alpha*X(i)
        END DO
        IF( N .LT. 4 ) RETURN
      END IF
      MP1 = M + 1
      DO i = MP1,N,4
        Y(i) = Y(i) + Alpha*X(i)
        Y(i + 1) = Y(i + 1) + Alpha*X(i + 1)
        Y(i + 2) = Y(i + 2) + Alpha*X(i + 2)
        Y(i + 3) = Y(i + 3) + Alpha*X(i + 3)
      END DO
      
      END SUBROUTINE saprc99_mosaic_4bin_vbs2_WAXPY




      SUBROUTINE saprc99_mosaic_4bin_vbs2_WSCAL(N,Alpha,X,incX)









      INTEGER i,incX,M,MP1,N
      REAL(kind=dp) X(N),Alpha
      REAL(kind=dp) ZERO, ONE
      PARAMETER( ZERO = 0.0_dp ) 
      PARAMETER( ONE  = 1.0_dp )

      IF (Alpha .EQ. ONE) RETURN
      IF (N .LE. 0) RETURN

      M = MOD(N,5)
      IF( M .NE. 0 ) THEN
        IF (Alpha .EQ. (-ONE)) THEN
          DO i = 1,M
            X(i) = -X(i)
          END DO
        ELSEIF (Alpha .EQ. ZERO) THEN
          DO i = 1,M
            X(i) = ZERO
          END DO
        ELSE
          DO i = 1,M
            X(i) = Alpha*X(i)
          END DO
        END IF
        IF( N .LT. 5 ) RETURN
      END IF
      MP1 = M + 1
      IF (Alpha .EQ. (-ONE)) THEN
        DO i = MP1,N,5
          X(i)     = -X(i)
          X(i + 1) = -X(i + 1)
          X(i + 2) = -X(i + 2)
          X(i + 3) = -X(i + 3)
          X(i + 4) = -X(i + 4)
        END DO
      ELSEIF (Alpha .EQ. ZERO) THEN
        DO i = MP1,N,5
          X(i)     = ZERO
          X(i + 1) = ZERO
          X(i + 2) = ZERO
          X(i + 3) = ZERO
          X(i + 4) = ZERO
        END DO
      ELSE
        DO i = MP1,N,5
          X(i)     = Alpha*X(i)
          X(i + 1) = Alpha*X(i + 1)
          X(i + 2) = Alpha*X(i + 2)
          X(i + 3) = Alpha*X(i + 3)
          X(i + 4) = Alpha*X(i + 4)
        END DO
      END IF

      END SUBROUTINE saprc99_mosaic_4bin_vbs2_WSCAL


      REAL(kind=dp) FUNCTION saprc99_mosaic_4bin_vbs2_WLAMCH( C )








      CHARACTER C
      INTEGER   i
      REAL(kind=dp)  ONE, HALF, Eps, Sum
      PARAMETER (ONE  = 1.0_dp)
      PARAMETER (HALF = 0.5_dp)
      LOGICAL   First
      SAVE     First, Eps
      DATA     First /.TRUE./
      
      IF (First) THEN
        First = .FALSE.
        Eps = HALF**(16)
        DO i = 17, 80
          Eps = Eps*HALF
          CALL saprc99_mosaic_4bin_vbs2_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      saprc99_mosaic_4bin_vbs2_WLAMCH = Eps

      END FUNCTION saprc99_mosaic_4bin_vbs2_WLAMCH
     
      SUBROUTINE saprc99_mosaic_4bin_vbs2_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE saprc99_mosaic_4bin_vbs2_WLAMCH_ADD




      SUBROUTINE saprc99_mosaic_4bin_vbs2_SET2ZERO(N,Y)




      
      INTEGER ::  i,M,MP1,N
      REAL(kind=dp) ::  Y(N)
      REAL(kind=dp), PARAMETER :: ZERO = 0.0d0

      IF (N.LE.0) RETURN

      M = MOD(N,8)
      IF( M .NE. 0 ) THEN
        DO i = 1,M
          Y(i) = ZERO
        END DO
        IF( N .LT. 8 ) RETURN
      END IF    
      MP1 = M+1
      DO i = MP1,N,8
        Y(i)     = ZERO
        Y(i + 1) = ZERO
        Y(i + 2) = ZERO
        Y(i + 3) = ZERO
        Y(i + 4) = ZERO
        Y(i + 5) = ZERO
        Y(i + 6) = ZERO
        Y(i + 7) = ZERO
      END DO

      END SUBROUTINE saprc99_mosaic_4bin_vbs2_SET2ZERO



      REAL(kind=dp) FUNCTION saprc99_mosaic_4bin_vbs2_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      saprc99_mosaic_4bin_vbs2_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        saprc99_mosaic_4bin_vbs2_WDOT = saprc99_mosaic_4bin_vbs2_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         saprc99_mosaic_4bin_vbs2_WDOT = saprc99_mosaic_4bin_vbs2_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          saprc99_mosaic_4bin_vbs2_WDOT = saprc99_mosaic_4bin_vbs2_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        saprc99_mosaic_4bin_vbs2_WDOT = saprc99_mosaic_4bin_vbs2_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION saprc99_mosaic_4bin_vbs2_WDOT                                          




   SUBROUTINE decomp_saprc99_mosaic_4bin_vbs2( JVS, IER )
   
     IMPLICIT NONE
   
      INTEGER  :: IER
      REAL(kind=dp) :: JVS(LU_NONZERO), W(NVAR), a
   
   
  a = 0._dp
  ier = 0 
   
  IF ( ABS(  JVS( 1 )) < TINY(a) ) THEN
         IER = 1                                       
         RETURN
  END IF
   W( 1 ) = JVS( 1 )
   W( 32 ) = JVS( 2 )
   W( 98 ) = JVS( 3 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 32 )
  JVS( 3) = W( 98 )
  IF ( ABS(  JVS( 4 )) < TINY(a) ) THEN
         IER = 2                                       
         RETURN
  END IF
   W( 2 ) = JVS( 4 )
   W( 59 ) = JVS( 5 )
   W( 69 ) = JVS( 6 )
   W( 88 ) = JVS( 7 )
  JVS( 4) = W( 2 )
  JVS( 5) = W( 59 )
  JVS( 6) = W( 69 )
  JVS( 7) = W( 88 )
  IF ( ABS(  JVS( 8 )) < TINY(a) ) THEN
         IER = 3                                       
         RETURN
  END IF
   W( 3 ) = JVS( 8 )
  JVS( 8) = W( 3 )
  IF ( ABS(  JVS( 9 )) < TINY(a) ) THEN
         IER = 4                                       
         RETURN
  END IF
   W( 4 ) = JVS( 9 )
  JVS( 9) = W( 4 )
  IF ( ABS(  JVS( 10 )) < TINY(a) ) THEN
         IER = 5                                       
         RETURN
  END IF
   W( 5 ) = JVS( 10 )
  JVS( 10) = W( 5 )
  IF ( ABS(  JVS( 11 )) < TINY(a) ) THEN
         IER = 6                                       
         RETURN
  END IF
   W( 6 ) = JVS( 11 )
  JVS( 11) = W( 6 )
  IF ( ABS(  JVS( 12 )) < TINY(a) ) THEN
         IER = 7                                       
         RETURN
  END IF
   W( 7 ) = JVS( 12 )
   W( 48 ) = JVS( 13 )
   W( 59 ) = JVS( 14 )
   W( 67 ) = JVS( 15 )
   W( 69 ) = JVS( 16 )
   W( 71 ) = JVS( 17 )
   W( 72 ) = JVS( 18 )
   W( 74 ) = JVS( 19 )
   W( 75 ) = JVS( 20 )
   W( 76 ) = JVS( 21 )
   W( 78 ) = JVS( 22 )
   W( 80 ) = JVS( 23 )
   W( 81 ) = JVS( 24 )
   W( 88 ) = JVS( 25 )
   W( 98 ) = JVS( 26 )
   W( 100 ) = JVS( 27 )
  JVS( 12) = W( 7 )
  JVS( 13) = W( 48 )
  JVS( 14) = W( 59 )
  JVS( 15) = W( 67 )
  JVS( 16) = W( 69 )
  JVS( 17) = W( 71 )
  JVS( 18) = W( 72 )
  JVS( 19) = W( 74 )
  JVS( 20) = W( 75 )
  JVS( 21) = W( 76 )
  JVS( 22) = W( 78 )
  JVS( 23) = W( 80 )
  JVS( 24) = W( 81 )
  JVS( 25) = W( 88 )
  JVS( 26) = W( 98 )
  JVS( 27) = W( 100 )
  IF ( ABS(  JVS( 28 )) < TINY(a) ) THEN
         IER = 8                                       
         RETURN
  END IF
   W( 8 ) = JVS( 28 )
   W( 69 ) = JVS( 29 )
   W( 72 ) = JVS( 30 )
   W( 80 ) = JVS( 31 )
   W( 88 ) = JVS( 32 )
   W( 89 ) = JVS( 33 )
   W( 90 ) = JVS( 34 )
   W( 91 ) = JVS( 35 )
   W( 94 ) = JVS( 36 )
   W( 97 ) = JVS( 37 )
  JVS( 28) = W( 8 )
  JVS( 29) = W( 69 )
  JVS( 30) = W( 72 )
  JVS( 31) = W( 80 )
  JVS( 32) = W( 88 )
  JVS( 33) = W( 89 )
  JVS( 34) = W( 90 )
  JVS( 35) = W( 91 )
  JVS( 36) = W( 94 )
  JVS( 37) = W( 97 )
  IF ( ABS(  JVS( 38 )) < TINY(a) ) THEN
         IER = 9                                       
         RETURN
  END IF
   W( 9 ) = JVS( 38 )
   W( 71 ) = JVS( 39 )
   W( 72 ) = JVS( 40 )
   W( 75 ) = JVS( 41 )
   W( 76 ) = JVS( 42 )
   W( 78 ) = JVS( 43 )
   W( 80 ) = JVS( 44 )
   W( 88 ) = JVS( 45 )
   W( 90 ) = JVS( 46 )
   W( 91 ) = JVS( 47 )
   W( 94 ) = JVS( 48 )
   W( 95 ) = JVS( 49 )
   W( 96 ) = JVS( 50 )
   W( 97 ) = JVS( 51 )
   W( 99 ) = JVS( 52 )
  JVS( 38) = W( 9 )
  JVS( 39) = W( 71 )
  JVS( 40) = W( 72 )
  JVS( 41) = W( 75 )
  JVS( 42) = W( 76 )
  JVS( 43) = W( 78 )
  JVS( 44) = W( 80 )
  JVS( 45) = W( 88 )
  JVS( 46) = W( 90 )
  JVS( 47) = W( 91 )
  JVS( 48) = W( 94 )
  JVS( 49) = W( 95 )
  JVS( 50) = W( 96 )
  JVS( 51) = W( 97 )
  JVS( 52) = W( 99 )
  IF ( ABS(  JVS( 53 )) < TINY(a) ) THEN
         IER = 10                                      
         RETURN
  END IF
   W( 10 ) = JVS( 53 )
   W( 89 ) = JVS( 54 )
   W( 97 ) = JVS( 55 )
  JVS( 53) = W( 10 )
  JVS( 54) = W( 89 )
  JVS( 55) = W( 97 )
  IF ( ABS(  JVS( 56 )) < TINY(a) ) THEN
         IER = 11                                      
         RETURN
  END IF
   W( 11 ) = JVS( 56 )
   W( 95 ) = JVS( 57 )
   W( 96 ) = JVS( 58 )
   W( 97 ) = JVS( 59 )
   W( 99 ) = JVS( 60 )
  JVS( 56) = W( 11 )
  JVS( 57) = W( 95 )
  JVS( 58) = W( 96 )
  JVS( 59) = W( 97 )
  JVS( 60) = W( 99 )
  IF ( ABS(  JVS( 61 )) < TINY(a) ) THEN
         IER = 12                                      
         RETURN
  END IF
   W( 12 ) = JVS( 61 )
   W( 49 ) = JVS( 62 )
   W( 69 ) = JVS( 63 )
   W( 92 ) = JVS( 64 )
   W( 93 ) = JVS( 65 )
  JVS( 61) = W( 12 )
  JVS( 62) = W( 49 )
  JVS( 63) = W( 69 )
  JVS( 64) = W( 92 )
  JVS( 65) = W( 93 )
  IF ( ABS(  JVS( 66 )) < TINY(a) ) THEN
         IER = 13                                      
         RETURN
  END IF
   W( 13 ) = JVS( 66 )
   W( 49 ) = JVS( 67 )
   W( 62 ) = JVS( 68 )
   W( 69 ) = JVS( 69 )
   W( 84 ) = JVS( 70 )
   W( 88 ) = JVS( 71 )
   W( 92 ) = JVS( 72 )
   W( 93 ) = JVS( 73 )
   W( 98 ) = JVS( 74 )
  JVS( 66) = W( 13 )
  JVS( 67) = W( 49 )
  JVS( 68) = W( 62 )
  JVS( 69) = W( 69 )
  JVS( 70) = W( 84 )
  JVS( 71) = W( 88 )
  JVS( 72) = W( 92 )
  JVS( 73) = W( 93 )
  JVS( 74) = W( 98 )
  IF ( ABS(  JVS( 75 )) < TINY(a) ) THEN
         IER = 14                                      
         RETURN
  END IF
   W( 14 ) = JVS( 75 )
   W( 73 ) = JVS( 76 )
   W( 90 ) = JVS( 77 )
   W( 91 ) = JVS( 78 )
   W( 93 ) = JVS( 79 )
   W( 94 ) = JVS( 80 )
   W( 100 ) = JVS( 81 )
  JVS( 75) = W( 14 )
  JVS( 76) = W( 73 )
  JVS( 77) = W( 90 )
  JVS( 78) = W( 91 )
  JVS( 79) = W( 93 )
  JVS( 80) = W( 94 )
  JVS( 81) = W( 100 )
  IF ( ABS(  JVS( 82 )) < TINY(a) ) THEN
         IER = 15                                      
         RETURN
  END IF
   W( 15 ) = JVS( 82 )
   W( 73 ) = JVS( 83 )
   W( 90 ) = JVS( 84 )
   W( 91 ) = JVS( 85 )
   W( 94 ) = JVS( 86 )
   W( 97 ) = JVS( 87 )
  JVS( 82) = W( 15 )
  JVS( 83) = W( 73 )
  JVS( 84) = W( 90 )
  JVS( 85) = W( 91 )
  JVS( 86) = W( 94 )
  JVS( 87) = W( 97 )
  IF ( ABS(  JVS( 88 )) < TINY(a) ) THEN
         IER = 16                                      
         RETURN
  END IF
   W( 16 ) = JVS( 88 )
   W( 46 ) = JVS( 89 )
   W( 53 ) = JVS( 90 )
   W( 72 ) = JVS( 91 )
   W( 80 ) = JVS( 92 )
   W( 98 ) = JVS( 93 )
  JVS( 88) = W( 16 )
  JVS( 89) = W( 46 )
  JVS( 90) = W( 53 )
  JVS( 91) = W( 72 )
  JVS( 92) = W( 80 )
  JVS( 93) = W( 98 )
  IF ( ABS(  JVS( 94 )) < TINY(a) ) THEN
         IER = 17                                      
         RETURN
  END IF
   W( 17 ) = JVS( 94 )
   W( 50 ) = JVS( 95 )
   W( 54 ) = JVS( 96 )
   W( 98 ) = JVS( 97 )
  JVS( 94) = W( 17 )
  JVS( 95) = W( 50 )
  JVS( 96) = W( 54 )
  JVS( 97) = W( 98 )
  IF ( ABS(  JVS( 98 )) < TINY(a) ) THEN
         IER = 18                                      
         RETURN
  END IF
   W( 18 ) = JVS( 98 )
   W( 71 ) = JVS( 99 )
   W( 98 ) = JVS( 100 )
  JVS( 98) = W( 18 )
  JVS( 99) = W( 71 )
  JVS( 100) = W( 98 )
  IF ( ABS(  JVS( 101 )) < TINY(a) ) THEN
         IER = 19                                      
         RETURN
  END IF
   W( 19 ) = JVS( 101 )
   W( 75 ) = JVS( 102 )
   W( 76 ) = JVS( 103 )
   W( 98 ) = JVS( 104 )
  JVS( 101) = W( 19 )
  JVS( 102) = W( 75 )
  JVS( 103) = W( 76 )
  JVS( 104) = W( 98 )
  IF ( ABS(  JVS( 105 )) < TINY(a) ) THEN
         IER = 20                                      
         RETURN
  END IF
   W( 20 ) = JVS( 105 )
   W( 31 ) = JVS( 106 )
   W( 32 ) = JVS( 107 )
   W( 33 ) = JVS( 108 )
   W( 39 ) = JVS( 109 )
   W( 41 ) = JVS( 110 )
   W( 43 ) = JVS( 111 )
   W( 44 ) = JVS( 112 )
   W( 46 ) = JVS( 113 )
   W( 47 ) = JVS( 114 )
   W( 50 ) = JVS( 115 )
   W( 52 ) = JVS( 116 )
   W( 53 ) = JVS( 117 )
   W( 54 ) = JVS( 118 )
   W( 55 ) = JVS( 119 )
   W( 56 ) = JVS( 120 )
   W( 57 ) = JVS( 121 )
   W( 58 ) = JVS( 122 )
   W( 59 ) = JVS( 123 )
   W( 60 ) = JVS( 124 )
   W( 62 ) = JVS( 125 )
   W( 63 ) = JVS( 126 )
   W( 64 ) = JVS( 127 )
   W( 65 ) = JVS( 128 )
   W( 66 ) = JVS( 129 )
   W( 67 ) = JVS( 130 )
   W( 68 ) = JVS( 131 )
   W( 69 ) = JVS( 132 )
   W( 71 ) = JVS( 133 )
   W( 72 ) = JVS( 134 )
   W( 74 ) = JVS( 135 )
   W( 75 ) = JVS( 136 )
   W( 76 ) = JVS( 137 )
   W( 77 ) = JVS( 138 )
   W( 78 ) = JVS( 139 )
   W( 79 ) = JVS( 140 )
   W( 80 ) = JVS( 141 )
   W( 81 ) = JVS( 142 )
   W( 82 ) = JVS( 143 )
   W( 83 ) = JVS( 144 )
   W( 85 ) = JVS( 145 )
   W( 86 ) = JVS( 146 )
   W( 87 ) = JVS( 147 )
   W( 88 ) = JVS( 148 )
   W( 92 ) = JVS( 149 )
   W( 93 ) = JVS( 150 )
   W( 97 ) = JVS( 151 )
   W( 98 ) = JVS( 152 )
   W( 100 ) = JVS( 153 )
  JVS( 105) = W( 20 )
  JVS( 106) = W( 31 )
  JVS( 107) = W( 32 )
  JVS( 108) = W( 33 )
  JVS( 109) = W( 39 )
  JVS( 110) = W( 41 )
  JVS( 111) = W( 43 )
  JVS( 112) = W( 44 )
  JVS( 113) = W( 46 )
  JVS( 114) = W( 47 )
  JVS( 115) = W( 50 )
  JVS( 116) = W( 52 )
  JVS( 117) = W( 53 )
  JVS( 118) = W( 54 )
  JVS( 119) = W( 55 )
  JVS( 120) = W( 56 )
  JVS( 121) = W( 57 )
  JVS( 122) = W( 58 )
  JVS( 123) = W( 59 )
  JVS( 124) = W( 60 )
  JVS( 125) = W( 62 )
  JVS( 126) = W( 63 )
  JVS( 127) = W( 64 )
  JVS( 128) = W( 65 )
  JVS( 129) = W( 66 )
  JVS( 130) = W( 67 )
  JVS( 131) = W( 68 )
  JVS( 132) = W( 69 )
  JVS( 133) = W( 71 )
  JVS( 134) = W( 72 )
  JVS( 135) = W( 74 )
  JVS( 136) = W( 75 )
  JVS( 137) = W( 76 )
  JVS( 138) = W( 77 )
  JVS( 139) = W( 78 )
  JVS( 140) = W( 79 )
  JVS( 141) = W( 80 )
  JVS( 142) = W( 81 )
  JVS( 143) = W( 82 )
  JVS( 144) = W( 83 )
  JVS( 145) = W( 85 )
  JVS( 146) = W( 86 )
  JVS( 147) = W( 87 )
  JVS( 148) = W( 88 )
  JVS( 149) = W( 92 )
  JVS( 150) = W( 93 )
  JVS( 151) = W( 97 )
  JVS( 152) = W( 98 )
  JVS( 153) = W( 100 )
  IF ( ABS(  JVS( 158 )) < TINY(a) ) THEN
         IER = 21                                      
         RETURN
  END IF
   W( 3 ) = JVS( 154 )
   W( 4 ) = JVS( 155 )
   W( 5 ) = JVS( 156 )
   W( 6 ) = JVS( 157 )
   W( 21 ) = JVS( 158 )
   W( 22 ) = JVS( 159 )
   W( 23 ) = JVS( 160 )
   W( 24 ) = JVS( 161 )
   W( 25 ) = JVS( 162 )
   W( 27 ) = JVS( 163 )
   W( 28 ) = JVS( 164 )
   W( 29 ) = JVS( 165 )
   W( 30 ) = JVS( 166 )
   W( 98 ) = JVS( 167 )
  a = -W( 3 ) / JVS(            8  )
  W( 3 ) = -a
  a = -W( 4 ) / JVS(            9  )
  W( 4 ) = -a
  a = -W( 5 ) / JVS(           10  )
  W( 5 ) = -a
  a = -W( 6 ) / JVS(           11  )
  W( 6 ) = -a
  JVS( 154) = W( 3 )
  JVS( 155) = W( 4 )
  JVS( 156) = W( 5 )
  JVS( 157) = W( 6 )
  JVS( 158) = W( 21 )
  JVS( 159) = W( 22 )
  JVS( 160) = W( 23 )
  JVS( 161) = W( 24 )
  JVS( 162) = W( 25 )
  JVS( 163) = W( 27 )
  JVS( 164) = W( 28 )
  JVS( 165) = W( 29 )
  JVS( 166) = W( 30 )
  JVS( 167) = W( 98 )
  IF ( ABS(  JVS( 168 )) < TINY(a) ) THEN
         IER = 22                                      
         RETURN
  END IF
   W( 22 ) = JVS( 168 )
   W( 23 ) = JVS( 169 )
   W( 27 ) = JVS( 170 )
   W( 98 ) = JVS( 171 )
  JVS( 168) = W( 22 )
  JVS( 169) = W( 23 )
  JVS( 170) = W( 27 )
  JVS( 171) = W( 98 )
  IF ( ABS(  JVS( 172 )) < TINY(a) ) THEN
         IER = 23                                      
         RETURN
  END IF
   W( 23 ) = JVS( 172 )
   W( 98 ) = JVS( 173 )
  JVS( 172) = W( 23 )
  JVS( 173) = W( 98 )
  IF ( ABS(  JVS( 174 )) < TINY(a) ) THEN
         IER = 24                                      
         RETURN
  END IF
   W( 24 ) = JVS( 174 )
   W( 25 ) = JVS( 175 )
   W( 29 ) = JVS( 176 )
   W( 98 ) = JVS( 177 )
  JVS( 174) = W( 24 )
  JVS( 175) = W( 25 )
  JVS( 176) = W( 29 )
  JVS( 177) = W( 98 )
  IF ( ABS(  JVS( 178 )) < TINY(a) ) THEN
         IER = 25                                      
         RETURN
  END IF
   W( 25 ) = JVS( 178 )
   W( 98 ) = JVS( 179 )
  JVS( 178) = W( 25 )
  JVS( 179) = W( 98 )
  IF ( ABS(  JVS( 180 )) < TINY(a) ) THEN
         IER = 26                                      
         RETURN
  END IF
   W( 26 ) = JVS( 180 )
   W( 88 ) = JVS( 181 )
  JVS( 180) = W( 26 )
  JVS( 181) = W( 88 )
  IF ( ABS(  JVS( 182 )) < TINY(a) ) THEN
         IER = 27                                      
         RETURN
  END IF
   W( 27 ) = JVS( 182 )
   W( 98 ) = JVS( 183 )
  JVS( 182) = W( 27 )
  JVS( 183) = W( 98 )
  IF ( ABS(  JVS( 185 )) < TINY(a) ) THEN
         IER = 28                                      
         RETURN
  END IF
   W( 27 ) = JVS( 184 )
   W( 28 ) = JVS( 185 )
   W( 98 ) = JVS( 186 )
  a = -W( 27 ) / JVS(          182  )
  W( 27 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 183 )
  JVS( 184) = W( 27 )
  JVS( 185) = W( 28 )
  JVS( 186) = W( 98 )
  IF ( ABS(  JVS( 187 )) < TINY(a) ) THEN
         IER = 29                                      
         RETURN
  END IF
   W( 29 ) = JVS( 187 )
   W( 98 ) = JVS( 188 )
  JVS( 187) = W( 29 )
  JVS( 188) = W( 98 )
  IF ( ABS(  JVS( 190 )) < TINY(a) ) THEN
         IER = 30                                      
         RETURN
  END IF
   W( 29 ) = JVS( 189 )
   W( 30 ) = JVS( 190 )
   W( 98 ) = JVS( 191 )
  a = -W( 29 ) / JVS(          187  )
  W( 29 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 188 )
  JVS( 189) = W( 29 )
  JVS( 190) = W( 30 )
  JVS( 191) = W( 98 )
  IF ( ABS(  JVS( 192 )) < TINY(a) ) THEN
         IER = 31                                      
         RETURN
  END IF
   W( 31 ) = JVS( 192 )
   W( 98 ) = JVS( 193 )
  JVS( 192) = W( 31 )
  JVS( 193) = W( 98 )
  IF ( ABS(  JVS( 194 )) < TINY(a) ) THEN
         IER = 32                                      
         RETURN
  END IF
   W( 32 ) = JVS( 194 )
   W( 98 ) = JVS( 195 )
  JVS( 194) = W( 32 )
  JVS( 195) = W( 98 )
  IF ( ABS(  JVS( 196 )) < TINY(a) ) THEN
         IER = 33                                      
         RETURN
  END IF
   W( 33 ) = JVS( 196 )
   W( 98 ) = JVS( 197 )
  JVS( 196) = W( 33 )
  JVS( 197) = W( 98 )
  IF ( ABS(  JVS( 198 )) < TINY(a) ) THEN
         IER = 34                                      
         RETURN
  END IF
   W( 34 ) = JVS( 198 )
   W( 89 ) = JVS( 199 )
   W( 92 ) = JVS( 200 )
  JVS( 198) = W( 34 )
  JVS( 199) = W( 89 )
  JVS( 200) = W( 92 )
  IF ( ABS(  JVS( 201 )) < TINY(a) ) THEN
         IER = 35                                      
         RETURN
  END IF
   W( 35 ) = JVS( 201 )
   W( 92 ) = JVS( 202 )
   W( 99 ) = JVS( 203 )
  JVS( 201) = W( 35 )
  JVS( 202) = W( 92 )
  JVS( 203) = W( 99 )
  IF ( ABS(  JVS( 204 )) < TINY(a) ) THEN
         IER = 36                                      
         RETURN
  END IF
   W( 36 ) = JVS( 204 )
   W( 92 ) = JVS( 205 )
   W( 95 ) = JVS( 206 )
  JVS( 204) = W( 36 )
  JVS( 205) = W( 92 )
  JVS( 206) = W( 95 )
  IF ( ABS(  JVS( 207 )) < TINY(a) ) THEN
         IER = 37                                      
         RETURN
  END IF
   W( 37 ) = JVS( 207 )
   W( 92 ) = JVS( 208 )
   W( 96 ) = JVS( 209 )
  JVS( 207) = W( 37 )
  JVS( 208) = W( 92 )
  JVS( 209) = W( 96 )
  IF ( ABS(  JVS( 210 )) < TINY(a) ) THEN
         IER = 38                                      
         RETURN
  END IF
   W( 38 ) = JVS( 210 )
   W( 97 ) = JVS( 211 )
   W( 98 ) = JVS( 212 )
  JVS( 210) = W( 38 )
  JVS( 211) = W( 97 )
  JVS( 212) = W( 98 )
  IF ( ABS(  JVS( 213 )) < TINY(a) ) THEN
         IER = 39                                      
         RETURN
  END IF
   W( 39 ) = JVS( 213 )
   W( 98 ) = JVS( 214 )
  JVS( 213) = W( 39 )
  JVS( 214) = W( 98 )
  IF ( ABS(  JVS( 215 )) < TINY(a) ) THEN
         IER = 40                                      
         RETURN
  END IF
   W( 40 ) = JVS( 215 )
   W( 50 ) = JVS( 216 )
   W( 75 ) = JVS( 217 )
   W( 76 ) = JVS( 218 )
   W( 88 ) = JVS( 219 )
   W( 98 ) = JVS( 220 )
  JVS( 215) = W( 40 )
  JVS( 216) = W( 50 )
  JVS( 217) = W( 75 )
  JVS( 218) = W( 76 )
  JVS( 219) = W( 88 )
  JVS( 220) = W( 98 )
  IF ( ABS(  JVS( 221 )) < TINY(a) ) THEN
         IER = 41                                      
         RETURN
  END IF
   W( 41 ) = JVS( 221 )
   W( 98 ) = JVS( 222 )
  JVS( 221) = W( 41 )
  JVS( 222) = W( 98 )
  IF ( ABS(  JVS( 223 )) < TINY(a) ) THEN
         IER = 42                                      
         RETURN
  END IF
   W( 42 ) = JVS( 223 )
   W( 92 ) = JVS( 224 )
   W( 93 ) = JVS( 225 )
  JVS( 223) = W( 42 )
  JVS( 224) = W( 92 )
  JVS( 225) = W( 93 )
  IF ( ABS(  JVS( 226 )) < TINY(a) ) THEN
         IER = 43                                      
         RETURN
  END IF
   W( 43 ) = JVS( 226 )
   W( 98 ) = JVS( 227 )
   W( 100 ) = JVS( 228 )
  JVS( 226) = W( 43 )
  JVS( 227) = W( 98 )
  JVS( 228) = W( 100 )
  IF ( ABS(  JVS( 229 )) < TINY(a) ) THEN
         IER = 44                                      
         RETURN
  END IF
   W( 44 ) = JVS( 229 )
   W( 98 ) = JVS( 230 )
  JVS( 229) = W( 44 )
  JVS( 230) = W( 98 )
  IF ( ABS(  JVS( 232 )) < TINY(a) ) THEN
         IER = 45                                      
         RETURN
  END IF
   W( 44 ) = JVS( 231 )
   W( 45 ) = JVS( 232 )
   W( 92 ) = JVS( 233 )
   W( 98 ) = JVS( 234 )
  a = -W( 44 ) / JVS(          229  )
  W( 44 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 230 )
  JVS( 231) = W( 44 )
  JVS( 232) = W( 45 )
  JVS( 233) = W( 92 )
  JVS( 234) = W( 98 )
  IF ( ABS(  JVS( 235 )) < TINY(a) ) THEN
         IER = 46                                      
         RETURN
  END IF
   W( 46 ) = JVS( 235 )
   W( 98 ) = JVS( 236 )
  JVS( 235) = W( 46 )
  JVS( 236) = W( 98 )
  IF ( ABS(  JVS( 237 )) < TINY(a) ) THEN
         IER = 47                                      
         RETURN
  END IF
   W( 47 ) = JVS( 237 )
   W( 94 ) = JVS( 238 )
   W( 97 ) = JVS( 239 )
   W( 98 ) = JVS( 240 )
  JVS( 237) = W( 47 )
  JVS( 238) = W( 94 )
  JVS( 239) = W( 97 )
  JVS( 240) = W( 98 )
  IF ( ABS(  JVS( 241 )) < TINY(a) ) THEN
         IER = 48                                      
         RETURN
  END IF
   W( 48 ) = JVS( 241 )
   W( 82 ) = JVS( 242 )
   W( 97 ) = JVS( 243 )
   W( 100 ) = JVS( 244 )
  JVS( 241) = W( 48 )
  JVS( 242) = W( 82 )
  JVS( 243) = W( 97 )
  JVS( 244) = W( 100 )
  IF ( ABS(  JVS( 245 )) < TINY(a) ) THEN
         IER = 49                                      
         RETURN
  END IF
   W( 49 ) = JVS( 245 )
   W( 61 ) = JVS( 246 )
   W( 92 ) = JVS( 247 )
   W( 93 ) = JVS( 248 )
   W( 97 ) = JVS( 249 )
  JVS( 245) = W( 49 )
  JVS( 246) = W( 61 )
  JVS( 247) = W( 92 )
  JVS( 248) = W( 93 )
  JVS( 249) = W( 97 )
  IF ( ABS(  JVS( 250 )) < TINY(a) ) THEN
         IER = 50                                      
         RETURN
  END IF
   W( 50 ) = JVS( 250 )
   W( 98 ) = JVS( 251 )
  JVS( 250) = W( 50 )
  JVS( 251) = W( 98 )
  IF ( ABS(  JVS( 252 )) < TINY(a) ) THEN
         IER = 51                                      
         RETURN
  END IF
   W( 51 ) = JVS( 252 )
   W( 92 ) = JVS( 253 )
   W( 97 ) = JVS( 254 )
   W( 98 ) = JVS( 255 )
  JVS( 252) = W( 51 )
  JVS( 253) = W( 92 )
  JVS( 254) = W( 97 )
  JVS( 255) = W( 98 )
  IF ( ABS(  JVS( 256 )) < TINY(a) ) THEN
         IER = 52                                      
         RETURN
  END IF
   W( 52 ) = JVS( 256 )
   W( 90 ) = JVS( 257 )
   W( 91 ) = JVS( 258 )
   W( 94 ) = JVS( 259 )
   W( 98 ) = JVS( 260 )
  JVS( 256) = W( 52 )
  JVS( 257) = W( 90 )
  JVS( 258) = W( 91 )
  JVS( 259) = W( 94 )
  JVS( 260) = W( 98 )
  IF ( ABS(  JVS( 261 )) < TINY(a) ) THEN
         IER = 53                                      
         RETURN
  END IF
   W( 53 ) = JVS( 261 )
   W( 98 ) = JVS( 262 )
  JVS( 261) = W( 53 )
  JVS( 262) = W( 98 )
  IF ( ABS(  JVS( 263 )) < TINY(a) ) THEN
         IER = 54                                      
         RETURN
  END IF
   W( 54 ) = JVS( 263 )
   W( 98 ) = JVS( 264 )
  JVS( 263) = W( 54 )
  JVS( 264) = W( 98 )
  IF ( ABS(  JVS( 267 )) < TINY(a) ) THEN
         IER = 55                                      
         RETURN
  END IF
   W( 50 ) = JVS( 265 )
   W( 54 ) = JVS( 266 )
   W( 55 ) = JVS( 267 )
   W( 98 ) = JVS( 268 )
  a = -W( 50 ) / JVS(          250  )
  W( 50 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 251 )
  a = -W( 54 ) / JVS(          263  )
  W( 54 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 264 )
  JVS( 265) = W( 50 )
  JVS( 266) = W( 54 )
  JVS( 267) = W( 55 )
  JVS( 268) = W( 98 )
  IF ( ABS(  JVS( 271 )) < TINY(a) ) THEN
         IER = 56                                      
         RETURN
  END IF
   W( 50 ) = JVS( 269 )
   W( 54 ) = JVS( 270 )
   W( 56 ) = JVS( 271 )
   W( 98 ) = JVS( 272 )
  a = -W( 50 ) / JVS(          250  )
  W( 50 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 251 )
  a = -W( 54 ) / JVS(          263  )
  W( 54 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 264 )
  JVS( 269) = W( 50 )
  JVS( 270) = W( 54 )
  JVS( 271) = W( 56 )
  JVS( 272) = W( 98 )
  IF ( ABS(  JVS( 275 )) < TINY(a) ) THEN
         IER = 57                                      
         RETURN
  END IF
   W( 50 ) = JVS( 273 )
   W( 54 ) = JVS( 274 )
   W( 57 ) = JVS( 275 )
   W( 93 ) = JVS( 276 )
   W( 98 ) = JVS( 277 )
  a = -W( 50 ) / JVS(          250  )
  W( 50 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 251 )
  a = -W( 54 ) / JVS(          263  )
  W( 54 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 264 )
  JVS( 273) = W( 50 )
  JVS( 274) = W( 54 )
  JVS( 275) = W( 57 )
  JVS( 276) = W( 93 )
  JVS( 277) = W( 98 )
  IF ( ABS(  JVS( 280 )) < TINY(a) ) THEN
         IER = 58                                      
         RETURN
  END IF
   W( 50 ) = JVS( 278 )
   W( 54 ) = JVS( 279 )
   W( 58 ) = JVS( 280 )
   W( 88 ) = JVS( 281 )
   W( 98 ) = JVS( 282 )
  a = -W( 50 ) / JVS(          250  )
  W( 50 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 251 )
  a = -W( 54 ) / JVS(          263  )
  W( 54 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 264 )
  JVS( 278) = W( 50 )
  JVS( 279) = W( 54 )
  JVS( 280) = W( 58 )
  JVS( 281) = W( 88 )
  JVS( 282) = W( 98 )
  IF ( ABS(  JVS( 283 )) < TINY(a) ) THEN
         IER = 59                                      
         RETURN
  END IF
   W( 59 ) = JVS( 283 )
   W( 88 ) = JVS( 284 )
   W( 98 ) = JVS( 285 )
  JVS( 283) = W( 59 )
  JVS( 284) = W( 88 )
  JVS( 285) = W( 98 )
  IF ( ABS(  JVS( 286 )) < TINY(a) ) THEN
         IER = 60                                      
         RETURN
  END IF
   W( 60 ) = JVS( 286 )
   W( 90 ) = JVS( 287 )
   W( 91 ) = JVS( 288 )
   W( 97 ) = JVS( 289 )
   W( 98 ) = JVS( 290 )
  JVS( 286) = W( 60 )
  JVS( 287) = W( 90 )
  JVS( 288) = W( 91 )
  JVS( 289) = W( 97 )
  JVS( 290) = W( 98 )
  IF ( ABS(  JVS( 292 )) < TINY(a) ) THEN
         IER = 61                                      
         RETURN
  END IF
   W( 49 ) = JVS( 291 )
   W( 61 ) = JVS( 292 )
   W( 70 ) = JVS( 293 )
   W( 92 ) = JVS( 294 )
   W( 93 ) = JVS( 295 )
   W( 97 ) = JVS( 296 )
  a = -W( 49 ) / JVS(          245  )
  W( 49 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 246 )
  W( 92 ) = W( 92 ) + a*JVS( 247 )
  W( 93 ) = W( 93 ) + a*JVS( 248 )
  W( 97 ) = W( 97 ) + a*JVS( 249 )
  JVS( 291) = W( 49 )
  JVS( 292) = W( 61 )
  JVS( 293) = W( 70 )
  JVS( 294) = W( 92 )
  JVS( 295) = W( 93 )
  JVS( 296) = W( 97 )
  IF ( ABS(  JVS( 299 )) < TINY(a) ) THEN
         IER = 62                                      
         RETURN
  END IF
   W( 50 ) = JVS( 297 )
   W( 54 ) = JVS( 298 )
   W( 62 ) = JVS( 299 )
   W( 80 ) = JVS( 300 )
   W( 88 ) = JVS( 301 )
   W( 93 ) = JVS( 302 )
   W( 98 ) = JVS( 303 )
  a = -W( 50 ) / JVS(          250  )
  W( 50 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 251 )
  a = -W( 54 ) / JVS(          263  )
  W( 54 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 264 )
  JVS( 297) = W( 50 )
  JVS( 298) = W( 54 )
  JVS( 299) = W( 62 )
  JVS( 300) = W( 80 )
  JVS( 301) = W( 88 )
  JVS( 302) = W( 93 )
  JVS( 303) = W( 98 )
  IF ( ABS(  JVS( 305 )) < TINY(a) ) THEN
         IER = 63                                      
         RETURN
  END IF
   W( 54 ) = JVS( 304 )
   W( 63 ) = JVS( 305 )
   W( 70 ) = JVS( 306 )
   W( 93 ) = JVS( 307 )
   W( 97 ) = JVS( 308 )
   W( 98 ) = JVS( 309 )
  a = -W( 54 ) / JVS(          263  )
  W( 54 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 264 )
  JVS( 304) = W( 54 )
  JVS( 305) = W( 63 )
  JVS( 306) = W( 70 )
  JVS( 307) = W( 93 )
  JVS( 308) = W( 97 )
  JVS( 309) = W( 98 )
  IF ( ABS(  JVS( 315 )) < TINY(a) ) THEN
         IER = 64                                      
         RETURN
  END IF
   W( 50 ) = JVS( 310 )
   W( 54 ) = JVS( 311 )
   W( 55 ) = JVS( 312 )
   W( 56 ) = JVS( 313 )
   W( 57 ) = JVS( 314 )
   W( 64 ) = JVS( 315 )
   W( 74 ) = JVS( 316 )
   W( 78 ) = JVS( 317 )
   W( 81 ) = JVS( 318 )
   W( 88 ) = JVS( 319 )
   W( 93 ) = JVS( 320 )
   W( 98 ) = JVS( 321 )
  a = -W( 50 ) / JVS(          250  )
  W( 50 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 251 )
  a = -W( 54 ) / JVS(          263  )
  W( 54 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 264 )
  a = -W( 55 ) / JVS(          267  )
  W( 55 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 268 )
  a = -W( 56 ) / JVS(          271  )
  W( 56 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 272 )
  a = -W( 57 ) / JVS(          275  )
  W( 57 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 276 )
  W( 98 ) = W( 98 ) + a*JVS( 277 )
  JVS( 310) = W( 50 )
  JVS( 311) = W( 54 )
  JVS( 312) = W( 55 )
  JVS( 313) = W( 56 )
  JVS( 314) = W( 57 )
  JVS( 315) = W( 64 )
  JVS( 316) = W( 74 )
  JVS( 317) = W( 78 )
  JVS( 318) = W( 81 )
  JVS( 319) = W( 88 )
  JVS( 320) = W( 93 )
  JVS( 321) = W( 98 )
  IF ( ABS(  JVS( 328 )) < TINY(a) ) THEN
         IER = 65                                      
         RETURN
  END IF
   W( 53 ) = JVS( 322 )
   W( 55 ) = JVS( 323 )
   W( 56 ) = JVS( 324 )
   W( 58 ) = JVS( 325 )
   W( 59 ) = JVS( 326 )
   W( 64 ) = JVS( 327 )
   W( 65 ) = JVS( 328 )
   W( 67 ) = JVS( 329 )
   W( 69 ) = JVS( 330 )
   W( 71 ) = JVS( 331 )
   W( 72 ) = JVS( 332 )
   W( 74 ) = JVS( 333 )
   W( 75 ) = JVS( 334 )
   W( 76 ) = JVS( 335 )
   W( 77 ) = JVS( 336 )
   W( 78 ) = JVS( 337 )
   W( 79 ) = JVS( 338 )
   W( 80 ) = JVS( 339 )
   W( 81 ) = JVS( 340 )
   W( 82 ) = JVS( 341 )
   W( 84 ) = JVS( 342 )
   W( 85 ) = JVS( 343 )
   W( 88 ) = JVS( 344 )
   W( 93 ) = JVS( 345 )
   W( 98 ) = JVS( 346 )
  a = -W( 53 ) / JVS(          261  )
  W( 53 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 262 )
  a = -W( 55 ) / JVS(          267  )
  W( 55 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 268 )
  a = -W( 56 ) / JVS(          271  )
  W( 56 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 272 )
  a = -W( 58 ) / JVS(          280  )
  W( 58 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 281 )
  W( 98 ) = W( 98 ) + a*JVS( 282 )
  a = -W( 59 ) / JVS(          283  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 284 )
  W( 98 ) = W( 98 ) + a*JVS( 285 )
  a = -W( 64 ) / JVS(          315  )
  W( 64 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 316 )
  W( 78 ) = W( 78 ) + a*JVS( 317 )
  W( 81 ) = W( 81 ) + a*JVS( 318 )
  W( 88 ) = W( 88 ) + a*JVS( 319 )
  W( 93 ) = W( 93 ) + a*JVS( 320 )
  W( 98 ) = W( 98 ) + a*JVS( 321 )
  JVS( 322) = W( 53 )
  JVS( 323) = W( 55 )
  JVS( 324) = W( 56 )
  JVS( 325) = W( 58 )
  JVS( 326) = W( 59 )
  JVS( 327) = W( 64 )
  JVS( 328) = W( 65 )
  JVS( 329) = W( 67 )
  JVS( 330) = W( 69 )
  JVS( 331) = W( 71 )
  JVS( 332) = W( 72 )
  JVS( 333) = W( 74 )
  JVS( 334) = W( 75 )
  JVS( 335) = W( 76 )
  JVS( 336) = W( 77 )
  JVS( 337) = W( 78 )
  JVS( 338) = W( 79 )
  JVS( 339) = W( 80 )
  JVS( 340) = W( 81 )
  JVS( 341) = W( 82 )
  JVS( 342) = W( 84 )
  JVS( 343) = W( 85 )
  JVS( 344) = W( 88 )
  JVS( 345) = W( 93 )
  JVS( 346) = W( 98 )
  IF ( ABS(  JVS( 353 )) < TINY(a) ) THEN
         IER = 66                                      
         RETURN
  END IF
   W( 42 ) = JVS( 347 )
   W( 57 ) = JVS( 348 )
   W( 61 ) = JVS( 349 )
   W( 62 ) = JVS( 350 )
   W( 63 ) = JVS( 351 )
   W( 64 ) = JVS( 352 )
   W( 66 ) = JVS( 353 )
   W( 70 ) = JVS( 354 )
   W( 74 ) = JVS( 355 )
   W( 77 ) = JVS( 356 )
   W( 78 ) = JVS( 357 )
   W( 79 ) = JVS( 358 )
   W( 80 ) = JVS( 359 )
   W( 81 ) = JVS( 360 )
   W( 82 ) = JVS( 361 )
   W( 85 ) = JVS( 362 )
   W( 88 ) = JVS( 363 )
   W( 92 ) = JVS( 364 )
   W( 93 ) = JVS( 365 )
   W( 97 ) = JVS( 366 )
   W( 98 ) = JVS( 367 )
  a = -W( 42 ) / JVS(          223  )
  W( 42 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 224 )
  W( 93 ) = W( 93 ) + a*JVS( 225 )
  a = -W( 57 ) / JVS(          275  )
  W( 57 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 276 )
  W( 98 ) = W( 98 ) + a*JVS( 277 )
  a = -W( 61 ) / JVS(          292  )
  W( 61 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 293 )
  W( 92 ) = W( 92 ) + a*JVS( 294 )
  W( 93 ) = W( 93 ) + a*JVS( 295 )
  W( 97 ) = W( 97 ) + a*JVS( 296 )
  a = -W( 62 ) / JVS(          299  )
  W( 62 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 300 )
  W( 88 ) = W( 88 ) + a*JVS( 301 )
  W( 93 ) = W( 93 ) + a*JVS( 302 )
  W( 98 ) = W( 98 ) + a*JVS( 303 )
  a = -W( 63 ) / JVS(          305  )
  W( 63 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 306 )
  W( 93 ) = W( 93 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  W( 98 ) = W( 98 ) + a*JVS( 309 )
  a = -W( 64 ) / JVS(          315  )
  W( 64 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 316 )
  W( 78 ) = W( 78 ) + a*JVS( 317 )
  W( 81 ) = W( 81 ) + a*JVS( 318 )
  W( 88 ) = W( 88 ) + a*JVS( 319 )
  W( 93 ) = W( 93 ) + a*JVS( 320 )
  W( 98 ) = W( 98 ) + a*JVS( 321 )
  JVS( 347) = W( 42 )
  JVS( 348) = W( 57 )
  JVS( 349) = W( 61 )
  JVS( 350) = W( 62 )
  JVS( 351) = W( 63 )
  JVS( 352) = W( 64 )
  JVS( 353) = W( 66 )
  JVS( 354) = W( 70 )
  JVS( 355) = W( 74 )
  JVS( 356) = W( 77 )
  JVS( 357) = W( 78 )
  JVS( 358) = W( 79 )
  JVS( 359) = W( 80 )
  JVS( 360) = W( 81 )
  JVS( 361) = W( 82 )
  JVS( 362) = W( 85 )
  JVS( 363) = W( 88 )
  JVS( 364) = W( 92 )
  JVS( 365) = W( 93 )
  JVS( 366) = W( 97 )
  JVS( 367) = W( 98 )
  IF ( ABS(  JVS( 368 )) < TINY(a) ) THEN
         IER = 67                                      
         RETURN
  END IF
   W( 67 ) = JVS( 368 )
   W( 84 ) = JVS( 369 )
   W( 88 ) = JVS( 370 )
   W( 93 ) = JVS( 371 )
   W( 98 ) = JVS( 372 )
  JVS( 368) = W( 67 )
  JVS( 369) = W( 84 )
  JVS( 370) = W( 88 )
  JVS( 371) = W( 93 )
  JVS( 372) = W( 98 )
  IF ( ABS(  JVS( 378 )) < TINY(a) ) THEN
         IER = 68                                      
         RETURN
  END IF
   W( 39 ) = JVS( 373 )
   W( 44 ) = JVS( 374 )
   W( 45 ) = JVS( 375 )
   W( 46 ) = JVS( 376 )
   W( 53 ) = JVS( 377 )
   W( 68 ) = JVS( 378 )
   W( 72 ) = JVS( 379 )
   W( 75 ) = JVS( 380 )
   W( 76 ) = JVS( 381 )
   W( 80 ) = JVS( 382 )
   W( 83 ) = JVS( 383 )
   W( 88 ) = JVS( 384 )
   W( 92 ) = JVS( 385 )
   W( 93 ) = JVS( 386 )
   W( 98 ) = JVS( 387 )
  a = -W( 39 ) / JVS(          213  )
  W( 39 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 214 )
  a = -W( 44 ) / JVS(          229  )
  W( 44 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 230 )
  a = -W( 45 ) / JVS(          232  )
  W( 45 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 233 )
  W( 98 ) = W( 98 ) + a*JVS( 234 )
  a = -W( 46 ) / JVS(          235  )
  W( 46 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 236 )
  a = -W( 53 ) / JVS(          261  )
  W( 53 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 262 )
  JVS( 373) = W( 39 )
  JVS( 374) = W( 44 )
  JVS( 375) = W( 45 )
  JVS( 376) = W( 46 )
  JVS( 377) = W( 53 )
  JVS( 378) = W( 68 )
  JVS( 379) = W( 72 )
  JVS( 380) = W( 75 )
  JVS( 381) = W( 76 )
  JVS( 382) = W( 80 )
  JVS( 383) = W( 83 )
  JVS( 384) = W( 88 )
  JVS( 385) = W( 92 )
  JVS( 386) = W( 93 )
  JVS( 387) = W( 98 )
  IF ( ABS(  JVS( 388 )) < TINY(a) ) THEN
         IER = 69                                      
         RETURN
  END IF
   W( 69 ) = JVS( 388 )
   W( 84 ) = JVS( 389 )
   W( 88 ) = JVS( 390 )
   W( 93 ) = JVS( 391 )
   W( 98 ) = JVS( 392 )
  JVS( 388) = W( 69 )
  JVS( 389) = W( 84 )
  JVS( 390) = W( 88 )
  JVS( 391) = W( 93 )
  JVS( 392) = W( 98 )
  IF ( ABS(  JVS( 395 )) < TINY(a) ) THEN
         IER = 70                                      
         RETURN
  END IF
   W( 57 ) = JVS( 393 )
   W( 63 ) = JVS( 394 )
   W( 70 ) = JVS( 395 )
   W( 89 ) = JVS( 396 )
   W( 92 ) = JVS( 397 )
   W( 93 ) = JVS( 398 )
   W( 95 ) = JVS( 399 )
   W( 96 ) = JVS( 400 )
   W( 97 ) = JVS( 401 )
   W( 98 ) = JVS( 402 )
   W( 99 ) = JVS( 403 )
   W( 100 ) = JVS( 404 )
  a = -W( 57 ) / JVS(          275  )
  W( 57 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 276 )
  W( 98 ) = W( 98 ) + a*JVS( 277 )
  a = -W( 63 ) / JVS(          305  )
  W( 63 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 306 )
  W( 93 ) = W( 93 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  W( 98 ) = W( 98 ) + a*JVS( 309 )
  JVS( 393) = W( 57 )
  JVS( 394) = W( 63 )
  JVS( 395) = W( 70 )
  JVS( 396) = W( 89 )
  JVS( 397) = W( 92 )
  JVS( 398) = W( 93 )
  JVS( 399) = W( 95 )
  JVS( 400) = W( 96 )
  JVS( 401) = W( 97 )
  JVS( 402) = W( 98 )
  JVS( 403) = W( 99 )
  JVS( 404) = W( 100 )
  IF ( ABS(  JVS( 405 )) < TINY(a) ) THEN
         IER = 71                                      
         RETURN
  END IF
   W( 71 ) = JVS( 405 )
   W( 84 ) = JVS( 406 )
   W( 88 ) = JVS( 407 )
   W( 93 ) = JVS( 408 )
   W( 98 ) = JVS( 409 )
  JVS( 405) = W( 71 )
  JVS( 406) = W( 84 )
  JVS( 407) = W( 88 )
  JVS( 408) = W( 93 )
  JVS( 409) = W( 98 )
  IF ( ABS(  JVS( 410 )) < TINY(a) ) THEN
         IER = 72                                      
         RETURN
  END IF
   W( 72 ) = JVS( 410 )
   W( 84 ) = JVS( 411 )
   W( 88 ) = JVS( 412 )
   W( 93 ) = JVS( 413 )
   W( 98 ) = JVS( 414 )
  JVS( 410) = W( 72 )
  JVS( 411) = W( 84 )
  JVS( 412) = W( 88 )
  JVS( 413) = W( 93 )
  JVS( 414) = W( 98 )
  IF ( ABS(  JVS( 423 )) < TINY(a) ) THEN
         IER = 73                                      
         RETURN
  END IF
   W( 44 ) = JVS( 415 )
   W( 46 ) = JVS( 416 )
   W( 53 ) = JVS( 417 )
   W( 55 ) = JVS( 418 )
   W( 56 ) = JVS( 419 )
   W( 68 ) = JVS( 420 )
   W( 71 ) = JVS( 421 )
   W( 72 ) = JVS( 422 )
   W( 73 ) = JVS( 423 )
   W( 75 ) = JVS( 424 )
   W( 76 ) = JVS( 425 )
   W( 80 ) = JVS( 426 )
   W( 81 ) = JVS( 427 )
   W( 83 ) = JVS( 428 )
   W( 84 ) = JVS( 429 )
   W( 86 ) = JVS( 430 )
   W( 87 ) = JVS( 431 )
   W( 88 ) = JVS( 432 )
   W( 89 ) = JVS( 433 )
   W( 90 ) = JVS( 434 )
   W( 91 ) = JVS( 435 )
   W( 92 ) = JVS( 436 )
   W( 93 ) = JVS( 437 )
   W( 94 ) = JVS( 438 )
   W( 95 ) = JVS( 439 )
   W( 96 ) = JVS( 440 )
   W( 97 ) = JVS( 441 )
   W( 98 ) = JVS( 442 )
   W( 99 ) = JVS( 443 )
   W( 100 ) = JVS( 444 )
  a = -W( 44 ) / JVS(          229  )
  W( 44 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 230 )
  a = -W( 46 ) / JVS(          235  )
  W( 46 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 236 )
  a = -W( 53 ) / JVS(          261  )
  W( 53 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 262 )
  a = -W( 55 ) / JVS(          267  )
  W( 55 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 268 )
  a = -W( 56 ) / JVS(          271  )
  W( 56 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 272 )
  a = -W( 68 ) / JVS(          378  )
  W( 68 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 80 ) = W( 80 ) + a*JVS( 382 )
  W( 83 ) = W( 83 ) + a*JVS( 383 )
  W( 88 ) = W( 88 ) + a*JVS( 384 )
  W( 92 ) = W( 92 ) + a*JVS( 385 )
  W( 93 ) = W( 93 ) + a*JVS( 386 )
  W( 98 ) = W( 98 ) + a*JVS( 387 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  JVS( 415) = W( 44 )
  JVS( 416) = W( 46 )
  JVS( 417) = W( 53 )
  JVS( 418) = W( 55 )
  JVS( 419) = W( 56 )
  JVS( 420) = W( 68 )
  JVS( 421) = W( 71 )
  JVS( 422) = W( 72 )
  JVS( 423) = W( 73 )
  JVS( 424) = W( 75 )
  JVS( 425) = W( 76 )
  JVS( 426) = W( 80 )
  JVS( 427) = W( 81 )
  JVS( 428) = W( 83 )
  JVS( 429) = W( 84 )
  JVS( 430) = W( 86 )
  JVS( 431) = W( 87 )
  JVS( 432) = W( 88 )
  JVS( 433) = W( 89 )
  JVS( 434) = W( 90 )
  JVS( 435) = W( 91 )
  JVS( 436) = W( 92 )
  JVS( 437) = W( 93 )
  JVS( 438) = W( 94 )
  JVS( 439) = W( 95 )
  JVS( 440) = W( 96 )
  JVS( 441) = W( 97 )
  JVS( 442) = W( 98 )
  JVS( 443) = W( 99 )
  JVS( 444) = W( 100 )
  IF ( ABS(  JVS( 446 )) < TINY(a) ) THEN
         IER = 74                                      
         RETURN
  END IF
   W( 71 ) = JVS( 445 )
   W( 74 ) = JVS( 446 )
   W( 80 ) = JVS( 447 )
   W( 84 ) = JVS( 448 )
   W( 88 ) = JVS( 449 )
   W( 93 ) = JVS( 450 )
   W( 98 ) = JVS( 451 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  JVS( 445) = W( 71 )
  JVS( 446) = W( 74 )
  JVS( 447) = W( 80 )
  JVS( 448) = W( 84 )
  JVS( 449) = W( 88 )
  JVS( 450) = W( 93 )
  JVS( 451) = W( 98 )
  IF ( ABS(  JVS( 452 )) < TINY(a) ) THEN
         IER = 75                                      
         RETURN
  END IF
   W( 75 ) = JVS( 452 )
   W( 84 ) = JVS( 453 )
   W( 88 ) = JVS( 454 )
   W( 93 ) = JVS( 455 )
   W( 98 ) = JVS( 456 )
  JVS( 452) = W( 75 )
  JVS( 453) = W( 84 )
  JVS( 454) = W( 88 )
  JVS( 455) = W( 93 )
  JVS( 456) = W( 98 )
  IF ( ABS(  JVS( 457 )) < TINY(a) ) THEN
         IER = 76                                      
         RETURN
  END IF
   W( 76 ) = JVS( 457 )
   W( 84 ) = JVS( 458 )
   W( 88 ) = JVS( 459 )
   W( 93 ) = JVS( 460 )
   W( 98 ) = JVS( 461 )
  JVS( 457) = W( 76 )
  JVS( 458) = W( 84 )
  JVS( 459) = W( 88 )
  JVS( 460) = W( 93 )
  JVS( 461) = W( 98 )
  IF ( ABS(  JVS( 473 )) < TINY(a) ) THEN
         IER = 77                                      
         RETURN
  END IF
   W( 50 ) = JVS( 462 )
   W( 54 ) = JVS( 463 )
   W( 55 ) = JVS( 464 )
   W( 56 ) = JVS( 465 )
   W( 58 ) = JVS( 466 )
   W( 59 ) = JVS( 467 )
   W( 63 ) = JVS( 468 )
   W( 67 ) = JVS( 469 )
   W( 70 ) = JVS( 470 )
   W( 75 ) = JVS( 471 )
   W( 76 ) = JVS( 472 )
   W( 77 ) = JVS( 473 )
   W( 78 ) = JVS( 474 )
   W( 84 ) = JVS( 475 )
   W( 88 ) = JVS( 476 )
   W( 89 ) = JVS( 477 )
   W( 92 ) = JVS( 478 )
   W( 93 ) = JVS( 479 )
   W( 95 ) = JVS( 480 )
   W( 96 ) = JVS( 481 )
   W( 97 ) = JVS( 482 )
   W( 98 ) = JVS( 483 )
   W( 99 ) = JVS( 484 )
   W( 100 ) = JVS( 485 )
  a = -W( 50 ) / JVS(          250  )
  W( 50 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 251 )
  a = -W( 54 ) / JVS(          263  )
  W( 54 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 264 )
  a = -W( 55 ) / JVS(          267  )
  W( 55 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 268 )
  a = -W( 56 ) / JVS(          271  )
  W( 56 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 272 )
  a = -W( 58 ) / JVS(          280  )
  W( 58 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 281 )
  W( 98 ) = W( 98 ) + a*JVS( 282 )
  a = -W( 59 ) / JVS(          283  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 284 )
  W( 98 ) = W( 98 ) + a*JVS( 285 )
  a = -W( 63 ) / JVS(          305  )
  W( 63 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 306 )
  W( 93 ) = W( 93 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  W( 98 ) = W( 98 ) + a*JVS( 309 )
  a = -W( 67 ) / JVS(          368  )
  W( 67 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 369 )
  W( 88 ) = W( 88 ) + a*JVS( 370 )
  W( 93 ) = W( 93 ) + a*JVS( 371 )
  W( 98 ) = W( 98 ) + a*JVS( 372 )
  a = -W( 70 ) / JVS(          395  )
  W( 70 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 396 )
  W( 92 ) = W( 92 ) + a*JVS( 397 )
  W( 93 ) = W( 93 ) + a*JVS( 398 )
  W( 95 ) = W( 95 ) + a*JVS( 399 )
  W( 96 ) = W( 96 ) + a*JVS( 400 )
  W( 97 ) = W( 97 ) + a*JVS( 401 )
  W( 98 ) = W( 98 ) + a*JVS( 402 )
  W( 99 ) = W( 99 ) + a*JVS( 403 )
  W( 100 ) = W( 100 ) + a*JVS( 404 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  JVS( 462) = W( 50 )
  JVS( 463) = W( 54 )
  JVS( 464) = W( 55 )
  JVS( 465) = W( 56 )
  JVS( 466) = W( 58 )
  JVS( 467) = W( 59 )
  JVS( 468) = W( 63 )
  JVS( 469) = W( 67 )
  JVS( 470) = W( 70 )
  JVS( 471) = W( 75 )
  JVS( 472) = W( 76 )
  JVS( 473) = W( 77 )
  JVS( 474) = W( 78 )
  JVS( 475) = W( 84 )
  JVS( 476) = W( 88 )
  JVS( 477) = W( 89 )
  JVS( 478) = W( 92 )
  JVS( 479) = W( 93 )
  JVS( 480) = W( 95 )
  JVS( 481) = W( 96 )
  JVS( 482) = W( 97 )
  JVS( 483) = W( 98 )
  JVS( 484) = W( 99 )
  JVS( 485) = W( 100 )
  IF ( ABS(  JVS( 487 )) < TINY(a) ) THEN
         IER = 78                                      
         RETURN
  END IF
   W( 71 ) = JVS( 486 )
   W( 78 ) = JVS( 487 )
   W( 80 ) = JVS( 488 )
   W( 84 ) = JVS( 489 )
   W( 88 ) = JVS( 490 )
   W( 93 ) = JVS( 491 )
   W( 98 ) = JVS( 492 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  JVS( 486) = W( 71 )
  JVS( 487) = W( 78 )
  JVS( 488) = W( 80 )
  JVS( 489) = W( 84 )
  JVS( 490) = W( 88 )
  JVS( 491) = W( 93 )
  JVS( 492) = W( 98 )
  IF ( ABS(  JVS( 502 )) < TINY(a) ) THEN
         IER = 79                                      
         RETURN
  END IF
   W( 33 ) = JVS( 493 )
   W( 41 ) = JVS( 494 )
   W( 44 ) = JVS( 495 )
   W( 46 ) = JVS( 496 )
   W( 53 ) = JVS( 497 )
   W( 67 ) = JVS( 498 )
   W( 69 ) = JVS( 499 )
   W( 72 ) = JVS( 500 )
   W( 78 ) = JVS( 501 )
   W( 79 ) = JVS( 502 )
   W( 80 ) = JVS( 503 )
   W( 83 ) = JVS( 504 )
   W( 84 ) = JVS( 505 )
   W( 85 ) = JVS( 506 )
   W( 86 ) = JVS( 507 )
   W( 87 ) = JVS( 508 )
   W( 88 ) = JVS( 509 )
   W( 89 ) = JVS( 510 )
   W( 93 ) = JVS( 511 )
   W( 95 ) = JVS( 512 )
   W( 96 ) = JVS( 513 )
   W( 98 ) = JVS( 514 )
   W( 99 ) = JVS( 515 )
   W( 100 ) = JVS( 516 )
  a = -W( 33 ) / JVS(          196  )
  W( 33 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 197 )
  a = -W( 41 ) / JVS(          221  )
  W( 41 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 222 )
  a = -W( 44 ) / JVS(          229  )
  W( 44 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 230 )
  a = -W( 46 ) / JVS(          235  )
  W( 46 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 236 )
  a = -W( 53 ) / JVS(          261  )
  W( 53 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 262 )
  a = -W( 67 ) / JVS(          368  )
  W( 67 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 369 )
  W( 88 ) = W( 88 ) + a*JVS( 370 )
  W( 93 ) = W( 93 ) + a*JVS( 371 )
  W( 98 ) = W( 98 ) + a*JVS( 372 )
  a = -W( 69 ) / JVS(          388  )
  W( 69 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 389 )
  W( 88 ) = W( 88 ) + a*JVS( 390 )
  W( 93 ) = W( 93 ) + a*JVS( 391 )
  W( 98 ) = W( 98 ) + a*JVS( 392 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  JVS( 493) = W( 33 )
  JVS( 494) = W( 41 )
  JVS( 495) = W( 44 )
  JVS( 496) = W( 46 )
  JVS( 497) = W( 53 )
  JVS( 498) = W( 67 )
  JVS( 499) = W( 69 )
  JVS( 500) = W( 72 )
  JVS( 501) = W( 78 )
  JVS( 502) = W( 79 )
  JVS( 503) = W( 80 )
  JVS( 504) = W( 83 )
  JVS( 505) = W( 84 )
  JVS( 506) = W( 85 )
  JVS( 507) = W( 86 )
  JVS( 508) = W( 87 )
  JVS( 509) = W( 88 )
  JVS( 510) = W( 89 )
  JVS( 511) = W( 93 )
  JVS( 512) = W( 95 )
  JVS( 513) = W( 96 )
  JVS( 514) = W( 98 )
  JVS( 515) = W( 99 )
  JVS( 516) = W( 100 )
  IF ( ABS(  JVS( 517 )) < TINY(a) ) THEN
         IER = 80                                      
         RETURN
  END IF
   W( 80 ) = JVS( 517 )
   W( 84 ) = JVS( 518 )
   W( 88 ) = JVS( 519 )
   W( 93 ) = JVS( 520 )
   W( 98 ) = JVS( 521 )
  JVS( 517) = W( 80 )
  JVS( 518) = W( 84 )
  JVS( 519) = W( 88 )
  JVS( 520) = W( 93 )
  JVS( 521) = W( 98 )
  IF ( ABS(  JVS( 524 )) < TINY(a) ) THEN
         IER = 81                                      
         RETURN
  END IF
   W( 71 ) = JVS( 522 )
   W( 80 ) = JVS( 523 )
   W( 81 ) = JVS( 524 )
   W( 84 ) = JVS( 525 )
   W( 88 ) = JVS( 526 )
   W( 93 ) = JVS( 527 )
   W( 98 ) = JVS( 528 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  JVS( 522) = W( 71 )
  JVS( 523) = W( 80 )
  JVS( 524) = W( 81 )
  JVS( 525) = W( 84 )
  JVS( 526) = W( 88 )
  JVS( 527) = W( 93 )
  JVS( 528) = W( 98 )
  IF ( ABS(  JVS( 549 )) < TINY(a) ) THEN
         IER = 82                                      
         RETURN
  END IF
   W( 41 ) = JVS( 529 )
   W( 44 ) = JVS( 530 )
   W( 46 ) = JVS( 531 )
   W( 47 ) = JVS( 532 )
   W( 48 ) = JVS( 533 )
   W( 52 ) = JVS( 534 )
   W( 53 ) = JVS( 535 )
   W( 59 ) = JVS( 536 )
   W( 67 ) = JVS( 537 )
   W( 68 ) = JVS( 538 )
   W( 69 ) = JVS( 539 )
   W( 71 ) = JVS( 540 )
   W( 72 ) = JVS( 541 )
   W( 74 ) = JVS( 542 )
   W( 75 ) = JVS( 543 )
   W( 76 ) = JVS( 544 )
   W( 77 ) = JVS( 545 )
   W( 78 ) = JVS( 546 )
   W( 80 ) = JVS( 547 )
   W( 81 ) = JVS( 548 )
   W( 82 ) = JVS( 549 )
   W( 83 ) = JVS( 550 )
   W( 84 ) = JVS( 551 )
   W( 86 ) = JVS( 552 )
   W( 87 ) = JVS( 553 )
   W( 88 ) = JVS( 554 )
   W( 89 ) = JVS( 555 )
   W( 90 ) = JVS( 556 )
   W( 91 ) = JVS( 557 )
   W( 92 ) = JVS( 558 )
   W( 93 ) = JVS( 559 )
   W( 94 ) = JVS( 560 )
   W( 95 ) = JVS( 561 )
   W( 96 ) = JVS( 562 )
   W( 97 ) = JVS( 563 )
   W( 98 ) = JVS( 564 )
   W( 99 ) = JVS( 565 )
   W( 100 ) = JVS( 566 )
  a = -W( 41 ) / JVS(          221  )
  W( 41 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 222 )
  a = -W( 44 ) / JVS(          229  )
  W( 44 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 230 )
  a = -W( 46 ) / JVS(          235  )
  W( 46 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 236 )
  a = -W( 47 ) / JVS(          237  )
  W( 47 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 238 )
  W( 97 ) = W( 97 ) + a*JVS( 239 )
  W( 98 ) = W( 98 ) + a*JVS( 240 )
  a = -W( 48 ) / JVS(          241  )
  W( 48 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 242 )
  W( 97 ) = W( 97 ) + a*JVS( 243 )
  W( 100 ) = W( 100 ) + a*JVS( 244 )
  a = -W( 52 ) / JVS(          256  )
  W( 52 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 257 )
  W( 91 ) = W( 91 ) + a*JVS( 258 )
  W( 94 ) = W( 94 ) + a*JVS( 259 )
  W( 98 ) = W( 98 ) + a*JVS( 260 )
  a = -W( 53 ) / JVS(          261  )
  W( 53 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 262 )
  a = -W( 59 ) / JVS(          283  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 284 )
  W( 98 ) = W( 98 ) + a*JVS( 285 )
  a = -W( 67 ) / JVS(          368  )
  W( 67 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 369 )
  W( 88 ) = W( 88 ) + a*JVS( 370 )
  W( 93 ) = W( 93 ) + a*JVS( 371 )
  W( 98 ) = W( 98 ) + a*JVS( 372 )
  a = -W( 68 ) / JVS(          378  )
  W( 68 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 80 ) = W( 80 ) + a*JVS( 382 )
  W( 83 ) = W( 83 ) + a*JVS( 383 )
  W( 88 ) = W( 88 ) + a*JVS( 384 )
  W( 92 ) = W( 92 ) + a*JVS( 385 )
  W( 93 ) = W( 93 ) + a*JVS( 386 )
  W( 98 ) = W( 98 ) + a*JVS( 387 )
  a = -W( 69 ) / JVS(          388  )
  W( 69 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 389 )
  W( 88 ) = W( 88 ) + a*JVS( 390 )
  W( 93 ) = W( 93 ) + a*JVS( 391 )
  W( 98 ) = W( 98 ) + a*JVS( 392 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 74 ) / JVS(          446  )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 447 )
  W( 84 ) = W( 84 ) + a*JVS( 448 )
  W( 88 ) = W( 88 ) + a*JVS( 449 )
  W( 93 ) = W( 93 ) + a*JVS( 450 )
  W( 98 ) = W( 98 ) + a*JVS( 451 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 77 ) / JVS(          473  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 474 )
  W( 84 ) = W( 84 ) + a*JVS( 475 )
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 89 ) = W( 89 ) + a*JVS( 477 )
  W( 92 ) = W( 92 ) + a*JVS( 478 )
  W( 93 ) = W( 93 ) + a*JVS( 479 )
  W( 95 ) = W( 95 ) + a*JVS( 480 )
  W( 96 ) = W( 96 ) + a*JVS( 481 )
  W( 97 ) = W( 97 ) + a*JVS( 482 )
  W( 98 ) = W( 98 ) + a*JVS( 483 )
  W( 99 ) = W( 99 ) + a*JVS( 484 )
  W( 100 ) = W( 100 ) + a*JVS( 485 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  JVS( 529) = W( 41 )
  JVS( 530) = W( 44 )
  JVS( 531) = W( 46 )
  JVS( 532) = W( 47 )
  JVS( 533) = W( 48 )
  JVS( 534) = W( 52 )
  JVS( 535) = W( 53 )
  JVS( 536) = W( 59 )
  JVS( 537) = W( 67 )
  JVS( 538) = W( 68 )
  JVS( 539) = W( 69 )
  JVS( 540) = W( 71 )
  JVS( 541) = W( 72 )
  JVS( 542) = W( 74 )
  JVS( 543) = W( 75 )
  JVS( 544) = W( 76 )
  JVS( 545) = W( 77 )
  JVS( 546) = W( 78 )
  JVS( 547) = W( 80 )
  JVS( 548) = W( 81 )
  JVS( 549) = W( 82 )
  JVS( 550) = W( 83 )
  JVS( 551) = W( 84 )
  JVS( 552) = W( 86 )
  JVS( 553) = W( 87 )
  JVS( 554) = W( 88 )
  JVS( 555) = W( 89 )
  JVS( 556) = W( 90 )
  JVS( 557) = W( 91 )
  JVS( 558) = W( 92 )
  JVS( 559) = W( 93 )
  JVS( 560) = W( 94 )
  JVS( 561) = W( 95 )
  JVS( 562) = W( 96 )
  JVS( 563) = W( 97 )
  JVS( 564) = W( 98 )
  JVS( 565) = W( 99 )
  JVS( 566) = W( 100 )
  IF ( ABS(  JVS( 573 )) < TINY(a) ) THEN
         IER = 83                                      
         RETURN
  END IF
   W( 45 ) = JVS( 567 )
   W( 72 ) = JVS( 568 )
   W( 75 ) = JVS( 569 )
   W( 76 ) = JVS( 570 )
   W( 78 ) = JVS( 571 )
   W( 80 ) = JVS( 572 )
   W( 83 ) = JVS( 573 )
   W( 84 ) = JVS( 574 )
   W( 88 ) = JVS( 575 )
   W( 91 ) = JVS( 576 )
   W( 92 ) = JVS( 577 )
   W( 93 ) = JVS( 578 )
   W( 98 ) = JVS( 579 )
   W( 100 ) = JVS( 580 )
  a = -W( 45 ) / JVS(          232  )
  W( 45 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 233 )
  W( 98 ) = W( 98 ) + a*JVS( 234 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  JVS( 567) = W( 45 )
  JVS( 568) = W( 72 )
  JVS( 569) = W( 75 )
  JVS( 570) = W( 76 )
  JVS( 571) = W( 78 )
  JVS( 572) = W( 80 )
  JVS( 573) = W( 83 )
  JVS( 574) = W( 84 )
  JVS( 575) = W( 88 )
  JVS( 576) = W( 91 )
  JVS( 577) = W( 92 )
  JVS( 578) = W( 93 )
  JVS( 579) = W( 98 )
  JVS( 580) = W( 100 )
  IF ( ABS(  JVS( 591 )) < TINY(a) ) THEN
         IER = 84                                      
         RETURN
  END IF
   W( 26 ) = JVS( 581 )
   W( 67 ) = JVS( 582 )
   W( 69 ) = JVS( 583 )
   W( 71 ) = JVS( 584 )
   W( 72 ) = JVS( 585 )
   W( 74 ) = JVS( 586 )
   W( 75 ) = JVS( 587 )
   W( 76 ) = JVS( 588 )
   W( 80 ) = JVS( 589 )
   W( 81 ) = JVS( 590 )
   W( 84 ) = JVS( 591 )
   W( 88 ) = JVS( 592 )
   W( 92 ) = JVS( 593 )
   W( 93 ) = JVS( 594 )
   W( 98 ) = JVS( 595 )
   W( 100 ) = JVS( 596 )
  a = -W( 26 ) / JVS(          180  )
  W( 26 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 181 )
  a = -W( 67 ) / JVS(          368  )
  W( 67 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 369 )
  W( 88 ) = W( 88 ) + a*JVS( 370 )
  W( 93 ) = W( 93 ) + a*JVS( 371 )
  W( 98 ) = W( 98 ) + a*JVS( 372 )
  a = -W( 69 ) / JVS(          388  )
  W( 69 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 389 )
  W( 88 ) = W( 88 ) + a*JVS( 390 )
  W( 93 ) = W( 93 ) + a*JVS( 391 )
  W( 98 ) = W( 98 ) + a*JVS( 392 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 74 ) / JVS(          446  )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 447 )
  W( 84 ) = W( 84 ) + a*JVS( 448 )
  W( 88 ) = W( 88 ) + a*JVS( 449 )
  W( 93 ) = W( 93 ) + a*JVS( 450 )
  W( 98 ) = W( 98 ) + a*JVS( 451 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  JVS( 581) = W( 26 )
  JVS( 582) = W( 67 )
  JVS( 583) = W( 69 )
  JVS( 584) = W( 71 )
  JVS( 585) = W( 72 )
  JVS( 586) = W( 74 )
  JVS( 587) = W( 75 )
  JVS( 588) = W( 76 )
  JVS( 589) = W( 80 )
  JVS( 590) = W( 81 )
  JVS( 591) = W( 84 )
  JVS( 592) = W( 88 )
  JVS( 593) = W( 92 )
  JVS( 594) = W( 93 )
  JVS( 595) = W( 98 )
  JVS( 596) = W( 100 )
  IF ( ABS(  JVS( 616 )) < TINY(a) ) THEN
         IER = 85                                      
         RETURN
  END IF
   W( 39 ) = JVS( 597 )
   W( 44 ) = JVS( 598 )
   W( 46 ) = JVS( 599 )
   W( 53 ) = JVS( 600 )
   W( 55 ) = JVS( 601 )
   W( 56 ) = JVS( 602 )
   W( 58 ) = JVS( 603 )
   W( 60 ) = JVS( 604 )
   W( 67 ) = JVS( 605 )
   W( 69 ) = JVS( 606 )
   W( 72 ) = JVS( 607 )
   W( 74 ) = JVS( 608 )
   W( 75 ) = JVS( 609 )
   W( 76 ) = JVS( 610 )
   W( 78 ) = JVS( 611 )
   W( 80 ) = JVS( 612 )
   W( 81 ) = JVS( 613 )
   W( 83 ) = JVS( 614 )
   W( 84 ) = JVS( 615 )
   W( 85 ) = JVS( 616 )
   W( 86 ) = JVS( 617 )
   W( 87 ) = JVS( 618 )
   W( 88 ) = JVS( 619 )
   W( 90 ) = JVS( 620 )
   W( 91 ) = JVS( 621 )
   W( 92 ) = JVS( 622 )
   W( 93 ) = JVS( 623 )
   W( 97 ) = JVS( 624 )
   W( 98 ) = JVS( 625 )
   W( 100 ) = JVS( 626 )
  a = -W( 39 ) / JVS(          213  )
  W( 39 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 214 )
  a = -W( 44 ) / JVS(          229  )
  W( 44 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 230 )
  a = -W( 46 ) / JVS(          235  )
  W( 46 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 236 )
  a = -W( 53 ) / JVS(          261  )
  W( 53 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 262 )
  a = -W( 55 ) / JVS(          267  )
  W( 55 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 268 )
  a = -W( 56 ) / JVS(          271  )
  W( 56 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 272 )
  a = -W( 58 ) / JVS(          280  )
  W( 58 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 281 )
  W( 98 ) = W( 98 ) + a*JVS( 282 )
  a = -W( 60 ) / JVS(          286  )
  W( 60 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 287 )
  W( 91 ) = W( 91 ) + a*JVS( 288 )
  W( 97 ) = W( 97 ) + a*JVS( 289 )
  W( 98 ) = W( 98 ) + a*JVS( 290 )
  a = -W( 67 ) / JVS(          368  )
  W( 67 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 369 )
  W( 88 ) = W( 88 ) + a*JVS( 370 )
  W( 93 ) = W( 93 ) + a*JVS( 371 )
  W( 98 ) = W( 98 ) + a*JVS( 372 )
  a = -W( 69 ) / JVS(          388  )
  W( 69 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 389 )
  W( 88 ) = W( 88 ) + a*JVS( 390 )
  W( 93 ) = W( 93 ) + a*JVS( 391 )
  W( 98 ) = W( 98 ) + a*JVS( 392 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 74 ) / JVS(          446  )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 447 )
  W( 84 ) = W( 84 ) + a*JVS( 448 )
  W( 88 ) = W( 88 ) + a*JVS( 449 )
  W( 93 ) = W( 93 ) + a*JVS( 450 )
  W( 98 ) = W( 98 ) + a*JVS( 451 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 83 ) / JVS(          573  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 574 )
  W( 88 ) = W( 88 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 98 ) = W( 98 ) + a*JVS( 579 )
  W( 100 ) = W( 100 ) + a*JVS( 580 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  JVS( 597) = W( 39 )
  JVS( 598) = W( 44 )
  JVS( 599) = W( 46 )
  JVS( 600) = W( 53 )
  JVS( 601) = W( 55 )
  JVS( 602) = W( 56 )
  JVS( 603) = W( 58 )
  JVS( 604) = W( 60 )
  JVS( 605) = W( 67 )
  JVS( 606) = W( 69 )
  JVS( 607) = W( 72 )
  JVS( 608) = W( 74 )
  JVS( 609) = W( 75 )
  JVS( 610) = W( 76 )
  JVS( 611) = W( 78 )
  JVS( 612) = W( 80 )
  JVS( 613) = W( 81 )
  JVS( 614) = W( 83 )
  JVS( 615) = W( 84 )
  JVS( 616) = W( 85 )
  JVS( 617) = W( 86 )
  JVS( 618) = W( 87 )
  JVS( 619) = W( 88 )
  JVS( 620) = W( 90 )
  JVS( 621) = W( 91 )
  JVS( 622) = W( 92 )
  JVS( 623) = W( 93 )
  JVS( 624) = W( 97 )
  JVS( 625) = W( 98 )
  JVS( 626) = W( 100 )
  IF ( ABS(  JVS( 638 )) < TINY(a) ) THEN
         IER = 86                                      
         RETURN
  END IF
   W( 44 ) = JVS( 627 )
   W( 46 ) = JVS( 628 )
   W( 53 ) = JVS( 629 )
   W( 69 ) = JVS( 630 )
   W( 72 ) = JVS( 631 )
   W( 74 ) = JVS( 632 )
   W( 78 ) = JVS( 633 )
   W( 80 ) = JVS( 634 )
   W( 81 ) = JVS( 635 )
   W( 83 ) = JVS( 636 )
   W( 84 ) = JVS( 637 )
   W( 86 ) = JVS( 638 )
   W( 87 ) = JVS( 639 )
   W( 88 ) = JVS( 640 )
   W( 90 ) = JVS( 641 )
   W( 91 ) = JVS( 642 )
   W( 92 ) = JVS( 643 )
   W( 93 ) = JVS( 644 )
   W( 94 ) = JVS( 645 )
   W( 98 ) = JVS( 646 )
   W( 100 ) = JVS( 647 )
  a = -W( 44 ) / JVS(          229  )
  W( 44 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 230 )
  a = -W( 46 ) / JVS(          235  )
  W( 46 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 236 )
  a = -W( 53 ) / JVS(          261  )
  W( 53 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 262 )
  a = -W( 69 ) / JVS(          388  )
  W( 69 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 389 )
  W( 88 ) = W( 88 ) + a*JVS( 390 )
  W( 93 ) = W( 93 ) + a*JVS( 391 )
  W( 98 ) = W( 98 ) + a*JVS( 392 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 74 ) / JVS(          446  )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 447 )
  W( 84 ) = W( 84 ) + a*JVS( 448 )
  W( 88 ) = W( 88 ) + a*JVS( 449 )
  W( 93 ) = W( 93 ) + a*JVS( 450 )
  W( 98 ) = W( 98 ) + a*JVS( 451 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 83 ) / JVS(          573  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 574 )
  W( 88 ) = W( 88 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 98 ) = W( 98 ) + a*JVS( 579 )
  W( 100 ) = W( 100 ) + a*JVS( 580 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  JVS( 627) = W( 44 )
  JVS( 628) = W( 46 )
  JVS( 629) = W( 53 )
  JVS( 630) = W( 69 )
  JVS( 631) = W( 72 )
  JVS( 632) = W( 74 )
  JVS( 633) = W( 78 )
  JVS( 634) = W( 80 )
  JVS( 635) = W( 81 )
  JVS( 636) = W( 83 )
  JVS( 637) = W( 84 )
  JVS( 638) = W( 86 )
  JVS( 639) = W( 87 )
  JVS( 640) = W( 88 )
  JVS( 641) = W( 90 )
  JVS( 642) = W( 91 )
  JVS( 643) = W( 92 )
  JVS( 644) = W( 93 )
  JVS( 645) = W( 94 )
  JVS( 646) = W( 98 )
  JVS( 647) = W( 100 )
  IF ( ABS(  JVS( 660 )) < TINY(a) ) THEN
         IER = 87                                      
         RETURN
  END IF
   W( 46 ) = JVS( 648 )
   W( 53 ) = JVS( 649 )
   W( 54 ) = JVS( 650 )
   W( 71 ) = JVS( 651 )
   W( 72 ) = JVS( 652 )
   W( 75 ) = JVS( 653 )
   W( 76 ) = JVS( 654 )
   W( 78 ) = JVS( 655 )
   W( 80 ) = JVS( 656 )
   W( 81 ) = JVS( 657 )
   W( 83 ) = JVS( 658 )
   W( 84 ) = JVS( 659 )
   W( 87 ) = JVS( 660 )
   W( 88 ) = JVS( 661 )
   W( 89 ) = JVS( 662 )
   W( 90 ) = JVS( 663 )
   W( 91 ) = JVS( 664 )
   W( 92 ) = JVS( 665 )
   W( 93 ) = JVS( 666 )
   W( 94 ) = JVS( 667 )
   W( 95 ) = JVS( 668 )
   W( 98 ) = JVS( 669 )
   W( 99 ) = JVS( 670 )
   W( 100 ) = JVS( 671 )
  a = -W( 46 ) / JVS(          235  )
  W( 46 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 236 )
  a = -W( 53 ) / JVS(          261  )
  W( 53 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 262 )
  a = -W( 54 ) / JVS(          263  )
  W( 54 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 264 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 83 ) / JVS(          573  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 574 )
  W( 88 ) = W( 88 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 98 ) = W( 98 ) + a*JVS( 579 )
  W( 100 ) = W( 100 ) + a*JVS( 580 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  JVS( 648) = W( 46 )
  JVS( 649) = W( 53 )
  JVS( 650) = W( 54 )
  JVS( 651) = W( 71 )
  JVS( 652) = W( 72 )
  JVS( 653) = W( 75 )
  JVS( 654) = W( 76 )
  JVS( 655) = W( 78 )
  JVS( 656) = W( 80 )
  JVS( 657) = W( 81 )
  JVS( 658) = W( 83 )
  JVS( 659) = W( 84 )
  JVS( 660) = W( 87 )
  JVS( 661) = W( 88 )
  JVS( 662) = W( 89 )
  JVS( 663) = W( 90 )
  JVS( 664) = W( 91 )
  JVS( 665) = W( 92 )
  JVS( 666) = W( 93 )
  JVS( 667) = W( 94 )
  JVS( 668) = W( 95 )
  JVS( 669) = W( 98 )
  JVS( 670) = W( 99 )
  JVS( 671) = W( 100 )
  IF ( ABS(  JVS( 685 )) < TINY(a) ) THEN
         IER = 88                                      
         RETURN
  END IF
   W( 58 ) = JVS( 672 )
   W( 59 ) = JVS( 673 )
   W( 67 ) = JVS( 674 )
   W( 69 ) = JVS( 675 )
   W( 71 ) = JVS( 676 )
   W( 72 ) = JVS( 677 )
   W( 74 ) = JVS( 678 )
   W( 75 ) = JVS( 679 )
   W( 76 ) = JVS( 680 )
   W( 78 ) = JVS( 681 )
   W( 80 ) = JVS( 682 )
   W( 81 ) = JVS( 683 )
   W( 84 ) = JVS( 684 )
   W( 88 ) = JVS( 685 )
   W( 89 ) = JVS( 686 )
   W( 92 ) = JVS( 687 )
   W( 93 ) = JVS( 688 )
   W( 95 ) = JVS( 689 )
   W( 96 ) = JVS( 690 )
   W( 97 ) = JVS( 691 )
   W( 98 ) = JVS( 692 )
   W( 99 ) = JVS( 693 )
   W( 100 ) = JVS( 694 )
  a = -W( 58 ) / JVS(          280  )
  W( 58 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 281 )
  W( 98 ) = W( 98 ) + a*JVS( 282 )
  a = -W( 59 ) / JVS(          283  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 284 )
  W( 98 ) = W( 98 ) + a*JVS( 285 )
  a = -W( 67 ) / JVS(          368  )
  W( 67 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 369 )
  W( 88 ) = W( 88 ) + a*JVS( 370 )
  W( 93 ) = W( 93 ) + a*JVS( 371 )
  W( 98 ) = W( 98 ) + a*JVS( 372 )
  a = -W( 69 ) / JVS(          388  )
  W( 69 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 389 )
  W( 88 ) = W( 88 ) + a*JVS( 390 )
  W( 93 ) = W( 93 ) + a*JVS( 391 )
  W( 98 ) = W( 98 ) + a*JVS( 392 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 74 ) / JVS(          446  )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 447 )
  W( 84 ) = W( 84 ) + a*JVS( 448 )
  W( 88 ) = W( 88 ) + a*JVS( 449 )
  W( 93 ) = W( 93 ) + a*JVS( 450 )
  W( 98 ) = W( 98 ) + a*JVS( 451 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  JVS( 672) = W( 58 )
  JVS( 673) = W( 59 )
  JVS( 674) = W( 67 )
  JVS( 675) = W( 69 )
  JVS( 676) = W( 71 )
  JVS( 677) = W( 72 )
  JVS( 678) = W( 74 )
  JVS( 679) = W( 75 )
  JVS( 680) = W( 76 )
  JVS( 681) = W( 78 )
  JVS( 682) = W( 80 )
  JVS( 683) = W( 81 )
  JVS( 684) = W( 84 )
  JVS( 685) = W( 88 )
  JVS( 686) = W( 89 )
  JVS( 687) = W( 92 )
  JVS( 688) = W( 93 )
  JVS( 689) = W( 95 )
  JVS( 690) = W( 96 )
  JVS( 691) = W( 97 )
  JVS( 692) = W( 98 )
  JVS( 693) = W( 99 )
  JVS( 694) = W( 100 )
  IF ( ABS(  JVS( 717 )) < TINY(a) ) THEN
         IER = 89                                      
         RETURN
  END IF
   W( 34 ) = JVS( 695 )
   W( 40 ) = JVS( 696 )
   W( 50 ) = JVS( 697 )
   W( 53 ) = JVS( 698 )
   W( 55 ) = JVS( 699 )
   W( 56 ) = JVS( 700 )
   W( 64 ) = JVS( 701 )
   W( 68 ) = JVS( 702 )
   W( 72 ) = JVS( 703 )
   W( 74 ) = JVS( 704 )
   W( 75 ) = JVS( 705 )
   W( 76 ) = JVS( 706 )
   W( 78 ) = JVS( 707 )
   W( 79 ) = JVS( 708 )
   W( 80 ) = JVS( 709 )
   W( 81 ) = JVS( 710 )
   W( 83 ) = JVS( 711 )
   W( 84 ) = JVS( 712 )
   W( 85 ) = JVS( 713 )
   W( 86 ) = JVS( 714 )
   W( 87 ) = JVS( 715 )
   W( 88 ) = JVS( 716 )
   W( 89 ) = JVS( 717 )
   W( 90 ) = JVS( 718 )
   W( 91 ) = JVS( 719 )
   W( 92 ) = JVS( 720 )
   W( 93 ) = JVS( 721 )
   W( 94 ) = JVS( 722 )
   W( 95 ) = JVS( 723 )
   W( 96 ) = JVS( 724 )
   W( 97 ) = JVS( 725 )
   W( 98 ) = JVS( 726 )
   W( 99 ) = JVS( 727 )
   W( 100 ) = JVS( 728 )
  a = -W( 34 ) / JVS(          198  )
  W( 34 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 199 )
  W( 92 ) = W( 92 ) + a*JVS( 200 )
  a = -W( 40 ) / JVS(          215  )
  W( 40 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 216 )
  W( 75 ) = W( 75 ) + a*JVS( 217 )
  W( 76 ) = W( 76 ) + a*JVS( 218 )
  W( 88 ) = W( 88 ) + a*JVS( 219 )
  W( 98 ) = W( 98 ) + a*JVS( 220 )
  a = -W( 50 ) / JVS(          250  )
  W( 50 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 251 )
  a = -W( 53 ) / JVS(          261  )
  W( 53 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 262 )
  a = -W( 55 ) / JVS(          267  )
  W( 55 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 268 )
  a = -W( 56 ) / JVS(          271  )
  W( 56 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 272 )
  a = -W( 64 ) / JVS(          315  )
  W( 64 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 316 )
  W( 78 ) = W( 78 ) + a*JVS( 317 )
  W( 81 ) = W( 81 ) + a*JVS( 318 )
  W( 88 ) = W( 88 ) + a*JVS( 319 )
  W( 93 ) = W( 93 ) + a*JVS( 320 )
  W( 98 ) = W( 98 ) + a*JVS( 321 )
  a = -W( 68 ) / JVS(          378  )
  W( 68 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 80 ) = W( 80 ) + a*JVS( 382 )
  W( 83 ) = W( 83 ) + a*JVS( 383 )
  W( 88 ) = W( 88 ) + a*JVS( 384 )
  W( 92 ) = W( 92 ) + a*JVS( 385 )
  W( 93 ) = W( 93 ) + a*JVS( 386 )
  W( 98 ) = W( 98 ) + a*JVS( 387 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 74 ) / JVS(          446  )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 447 )
  W( 84 ) = W( 84 ) + a*JVS( 448 )
  W( 88 ) = W( 88 ) + a*JVS( 449 )
  W( 93 ) = W( 93 ) + a*JVS( 450 )
  W( 98 ) = W( 98 ) + a*JVS( 451 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 79 ) / JVS(          502  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 503 )
  W( 83 ) = W( 83 ) + a*JVS( 504 )
  W( 84 ) = W( 84 ) + a*JVS( 505 )
  W( 85 ) = W( 85 ) + a*JVS( 506 )
  W( 86 ) = W( 86 ) + a*JVS( 507 )
  W( 87 ) = W( 87 ) + a*JVS( 508 )
  W( 88 ) = W( 88 ) + a*JVS( 509 )
  W( 89 ) = W( 89 ) + a*JVS( 510 )
  W( 93 ) = W( 93 ) + a*JVS( 511 )
  W( 95 ) = W( 95 ) + a*JVS( 512 )
  W( 96 ) = W( 96 ) + a*JVS( 513 )
  W( 98 ) = W( 98 ) + a*JVS( 514 )
  W( 99 ) = W( 99 ) + a*JVS( 515 )
  W( 100 ) = W( 100 ) + a*JVS( 516 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 83 ) / JVS(          573  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 574 )
  W( 88 ) = W( 88 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 98 ) = W( 98 ) + a*JVS( 579 )
  W( 100 ) = W( 100 ) + a*JVS( 580 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  a = -W( 85 ) / JVS(          616  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 617 )
  W( 87 ) = W( 87 ) + a*JVS( 618 )
  W( 88 ) = W( 88 ) + a*JVS( 619 )
  W( 90 ) = W( 90 ) + a*JVS( 620 )
  W( 91 ) = W( 91 ) + a*JVS( 621 )
  W( 92 ) = W( 92 ) + a*JVS( 622 )
  W( 93 ) = W( 93 ) + a*JVS( 623 )
  W( 97 ) = W( 97 ) + a*JVS( 624 )
  W( 98 ) = W( 98 ) + a*JVS( 625 )
  W( 100 ) = W( 100 ) + a*JVS( 626 )
  a = -W( 86 ) / JVS(          638  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 639 )
  W( 88 ) = W( 88 ) + a*JVS( 640 )
  W( 90 ) = W( 90 ) + a*JVS( 641 )
  W( 91 ) = W( 91 ) + a*JVS( 642 )
  W( 92 ) = W( 92 ) + a*JVS( 643 )
  W( 93 ) = W( 93 ) + a*JVS( 644 )
  W( 94 ) = W( 94 ) + a*JVS( 645 )
  W( 98 ) = W( 98 ) + a*JVS( 646 )
  W( 100 ) = W( 100 ) + a*JVS( 647 )
  a = -W( 87 ) / JVS(          660  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 661 )
  W( 89 ) = W( 89 ) + a*JVS( 662 )
  W( 90 ) = W( 90 ) + a*JVS( 663 )
  W( 91 ) = W( 91 ) + a*JVS( 664 )
  W( 92 ) = W( 92 ) + a*JVS( 665 )
  W( 93 ) = W( 93 ) + a*JVS( 666 )
  W( 94 ) = W( 94 ) + a*JVS( 667 )
  W( 95 ) = W( 95 ) + a*JVS( 668 )
  W( 98 ) = W( 98 ) + a*JVS( 669 )
  W( 99 ) = W( 99 ) + a*JVS( 670 )
  W( 100 ) = W( 100 ) + a*JVS( 671 )
  a = -W( 88 ) / JVS(          685  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 686 )
  W( 92 ) = W( 92 ) + a*JVS( 687 )
  W( 93 ) = W( 93 ) + a*JVS( 688 )
  W( 95 ) = W( 95 ) + a*JVS( 689 )
  W( 96 ) = W( 96 ) + a*JVS( 690 )
  W( 97 ) = W( 97 ) + a*JVS( 691 )
  W( 98 ) = W( 98 ) + a*JVS( 692 )
  W( 99 ) = W( 99 ) + a*JVS( 693 )
  W( 100 ) = W( 100 ) + a*JVS( 694 )
  JVS( 695) = W( 34 )
  JVS( 696) = W( 40 )
  JVS( 697) = W( 50 )
  JVS( 698) = W( 53 )
  JVS( 699) = W( 55 )
  JVS( 700) = W( 56 )
  JVS( 701) = W( 64 )
  JVS( 702) = W( 68 )
  JVS( 703) = W( 72 )
  JVS( 704) = W( 74 )
  JVS( 705) = W( 75 )
  JVS( 706) = W( 76 )
  JVS( 707) = W( 78 )
  JVS( 708) = W( 79 )
  JVS( 709) = W( 80 )
  JVS( 710) = W( 81 )
  JVS( 711) = W( 83 )
  JVS( 712) = W( 84 )
  JVS( 713) = W( 85 )
  JVS( 714) = W( 86 )
  JVS( 715) = W( 87 )
  JVS( 716) = W( 88 )
  JVS( 717) = W( 89 )
  JVS( 718) = W( 90 )
  JVS( 719) = W( 91 )
  JVS( 720) = W( 92 )
  JVS( 721) = W( 93 )
  JVS( 722) = W( 94 )
  JVS( 723) = W( 95 )
  JVS( 724) = W( 96 )
  JVS( 725) = W( 97 )
  JVS( 726) = W( 98 )
  JVS( 727) = W( 99 )
  JVS( 728) = W( 100 )
  IF ( ABS(  JVS( 762 )) < TINY(a) ) THEN
         IER = 90                                      
         RETURN
  END IF
   W( 33 ) = JVS( 729 )
   W( 39 ) = JVS( 730 )
   W( 41 ) = JVS( 731 )
   W( 44 ) = JVS( 732 )
   W( 46 ) = JVS( 733 )
   W( 50 ) = JVS( 734 )
   W( 53 ) = JVS( 735 )
   W( 54 ) = JVS( 736 )
   W( 55 ) = JVS( 737 )
   W( 56 ) = JVS( 738 )
   W( 57 ) = JVS( 739 )
   W( 58 ) = JVS( 740 )
   W( 59 ) = JVS( 741 )
   W( 60 ) = JVS( 742 )
   W( 63 ) = JVS( 743 )
   W( 67 ) = JVS( 744 )
   W( 69 ) = JVS( 745 )
   W( 70 ) = JVS( 746 )
   W( 71 ) = JVS( 747 )
   W( 72 ) = JVS( 748 )
   W( 74 ) = JVS( 749 )
   W( 75 ) = JVS( 750 )
   W( 76 ) = JVS( 751 )
   W( 78 ) = JVS( 752 )
   W( 80 ) = JVS( 753 )
   W( 81 ) = JVS( 754 )
   W( 83 ) = JVS( 755 )
   W( 84 ) = JVS( 756 )
   W( 85 ) = JVS( 757 )
   W( 86 ) = JVS( 758 )
   W( 87 ) = JVS( 759 )
   W( 88 ) = JVS( 760 )
   W( 89 ) = JVS( 761 )
   W( 90 ) = JVS( 762 )
   W( 91 ) = JVS( 763 )
   W( 92 ) = JVS( 764 )
   W( 93 ) = JVS( 765 )
   W( 94 ) = JVS( 766 )
   W( 95 ) = JVS( 767 )
   W( 96 ) = JVS( 768 )
   W( 97 ) = JVS( 769 )
   W( 98 ) = JVS( 770 )
   W( 99 ) = JVS( 771 )
   W( 100 ) = JVS( 772 )
  a = -W( 33 ) / JVS(          196  )
  W( 33 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 197 )
  a = -W( 39 ) / JVS(          213  )
  W( 39 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 214 )
  a = -W( 41 ) / JVS(          221  )
  W( 41 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 222 )
  a = -W( 44 ) / JVS(          229  )
  W( 44 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 230 )
  a = -W( 46 ) / JVS(          235  )
  W( 46 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 236 )
  a = -W( 50 ) / JVS(          250  )
  W( 50 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 251 )
  a = -W( 53 ) / JVS(          261  )
  W( 53 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 262 )
  a = -W( 54 ) / JVS(          263  )
  W( 54 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 264 )
  a = -W( 55 ) / JVS(          267  )
  W( 55 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 268 )
  a = -W( 56 ) / JVS(          271  )
  W( 56 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 272 )
  a = -W( 57 ) / JVS(          275  )
  W( 57 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 276 )
  W( 98 ) = W( 98 ) + a*JVS( 277 )
  a = -W( 58 ) / JVS(          280  )
  W( 58 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 281 )
  W( 98 ) = W( 98 ) + a*JVS( 282 )
  a = -W( 59 ) / JVS(          283  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 284 )
  W( 98 ) = W( 98 ) + a*JVS( 285 )
  a = -W( 60 ) / JVS(          286  )
  W( 60 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 287 )
  W( 91 ) = W( 91 ) + a*JVS( 288 )
  W( 97 ) = W( 97 ) + a*JVS( 289 )
  W( 98 ) = W( 98 ) + a*JVS( 290 )
  a = -W( 63 ) / JVS(          305  )
  W( 63 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 306 )
  W( 93 ) = W( 93 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  W( 98 ) = W( 98 ) + a*JVS( 309 )
  a = -W( 67 ) / JVS(          368  )
  W( 67 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 369 )
  W( 88 ) = W( 88 ) + a*JVS( 370 )
  W( 93 ) = W( 93 ) + a*JVS( 371 )
  W( 98 ) = W( 98 ) + a*JVS( 372 )
  a = -W( 69 ) / JVS(          388  )
  W( 69 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 389 )
  W( 88 ) = W( 88 ) + a*JVS( 390 )
  W( 93 ) = W( 93 ) + a*JVS( 391 )
  W( 98 ) = W( 98 ) + a*JVS( 392 )
  a = -W( 70 ) / JVS(          395  )
  W( 70 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 396 )
  W( 92 ) = W( 92 ) + a*JVS( 397 )
  W( 93 ) = W( 93 ) + a*JVS( 398 )
  W( 95 ) = W( 95 ) + a*JVS( 399 )
  W( 96 ) = W( 96 ) + a*JVS( 400 )
  W( 97 ) = W( 97 ) + a*JVS( 401 )
  W( 98 ) = W( 98 ) + a*JVS( 402 )
  W( 99 ) = W( 99 ) + a*JVS( 403 )
  W( 100 ) = W( 100 ) + a*JVS( 404 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 74 ) / JVS(          446  )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 447 )
  W( 84 ) = W( 84 ) + a*JVS( 448 )
  W( 88 ) = W( 88 ) + a*JVS( 449 )
  W( 93 ) = W( 93 ) + a*JVS( 450 )
  W( 98 ) = W( 98 ) + a*JVS( 451 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 83 ) / JVS(          573  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 574 )
  W( 88 ) = W( 88 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 98 ) = W( 98 ) + a*JVS( 579 )
  W( 100 ) = W( 100 ) + a*JVS( 580 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  a = -W( 85 ) / JVS(          616  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 617 )
  W( 87 ) = W( 87 ) + a*JVS( 618 )
  W( 88 ) = W( 88 ) + a*JVS( 619 )
  W( 90 ) = W( 90 ) + a*JVS( 620 )
  W( 91 ) = W( 91 ) + a*JVS( 621 )
  W( 92 ) = W( 92 ) + a*JVS( 622 )
  W( 93 ) = W( 93 ) + a*JVS( 623 )
  W( 97 ) = W( 97 ) + a*JVS( 624 )
  W( 98 ) = W( 98 ) + a*JVS( 625 )
  W( 100 ) = W( 100 ) + a*JVS( 626 )
  a = -W( 86 ) / JVS(          638  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 639 )
  W( 88 ) = W( 88 ) + a*JVS( 640 )
  W( 90 ) = W( 90 ) + a*JVS( 641 )
  W( 91 ) = W( 91 ) + a*JVS( 642 )
  W( 92 ) = W( 92 ) + a*JVS( 643 )
  W( 93 ) = W( 93 ) + a*JVS( 644 )
  W( 94 ) = W( 94 ) + a*JVS( 645 )
  W( 98 ) = W( 98 ) + a*JVS( 646 )
  W( 100 ) = W( 100 ) + a*JVS( 647 )
  a = -W( 87 ) / JVS(          660  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 661 )
  W( 89 ) = W( 89 ) + a*JVS( 662 )
  W( 90 ) = W( 90 ) + a*JVS( 663 )
  W( 91 ) = W( 91 ) + a*JVS( 664 )
  W( 92 ) = W( 92 ) + a*JVS( 665 )
  W( 93 ) = W( 93 ) + a*JVS( 666 )
  W( 94 ) = W( 94 ) + a*JVS( 667 )
  W( 95 ) = W( 95 ) + a*JVS( 668 )
  W( 98 ) = W( 98 ) + a*JVS( 669 )
  W( 99 ) = W( 99 ) + a*JVS( 670 )
  W( 100 ) = W( 100 ) + a*JVS( 671 )
  a = -W( 88 ) / JVS(          685  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 686 )
  W( 92 ) = W( 92 ) + a*JVS( 687 )
  W( 93 ) = W( 93 ) + a*JVS( 688 )
  W( 95 ) = W( 95 ) + a*JVS( 689 )
  W( 96 ) = W( 96 ) + a*JVS( 690 )
  W( 97 ) = W( 97 ) + a*JVS( 691 )
  W( 98 ) = W( 98 ) + a*JVS( 692 )
  W( 99 ) = W( 99 ) + a*JVS( 693 )
  W( 100 ) = W( 100 ) + a*JVS( 694 )
  a = -W( 89 ) / JVS(          717  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 718 )
  W( 91 ) = W( 91 ) + a*JVS( 719 )
  W( 92 ) = W( 92 ) + a*JVS( 720 )
  W( 93 ) = W( 93 ) + a*JVS( 721 )
  W( 94 ) = W( 94 ) + a*JVS( 722 )
  W( 95 ) = W( 95 ) + a*JVS( 723 )
  W( 96 ) = W( 96 ) + a*JVS( 724 )
  W( 97 ) = W( 97 ) + a*JVS( 725 )
  W( 98 ) = W( 98 ) + a*JVS( 726 )
  W( 99 ) = W( 99 ) + a*JVS( 727 )
  W( 100 ) = W( 100 ) + a*JVS( 728 )
  JVS( 729) = W( 33 )
  JVS( 730) = W( 39 )
  JVS( 731) = W( 41 )
  JVS( 732) = W( 44 )
  JVS( 733) = W( 46 )
  JVS( 734) = W( 50 )
  JVS( 735) = W( 53 )
  JVS( 736) = W( 54 )
  JVS( 737) = W( 55 )
  JVS( 738) = W( 56 )
  JVS( 739) = W( 57 )
  JVS( 740) = W( 58 )
  JVS( 741) = W( 59 )
  JVS( 742) = W( 60 )
  JVS( 743) = W( 63 )
  JVS( 744) = W( 67 )
  JVS( 745) = W( 69 )
  JVS( 746) = W( 70 )
  JVS( 747) = W( 71 )
  JVS( 748) = W( 72 )
  JVS( 749) = W( 74 )
  JVS( 750) = W( 75 )
  JVS( 751) = W( 76 )
  JVS( 752) = W( 78 )
  JVS( 753) = W( 80 )
  JVS( 754) = W( 81 )
  JVS( 755) = W( 83 )
  JVS( 756) = W( 84 )
  JVS( 757) = W( 85 )
  JVS( 758) = W( 86 )
  JVS( 759) = W( 87 )
  JVS( 760) = W( 88 )
  JVS( 761) = W( 89 )
  JVS( 762) = W( 90 )
  JVS( 763) = W( 91 )
  JVS( 764) = W( 92 )
  JVS( 765) = W( 93 )
  JVS( 766) = W( 94 )
  JVS( 767) = W( 95 )
  JVS( 768) = W( 96 )
  JVS( 769) = W( 97 )
  JVS( 770) = W( 98 )
  JVS( 771) = W( 99 )
  JVS( 772) = W( 100 )
  IF ( ABS(  JVS( 795 )) < TINY(a) ) THEN
         IER = 91                                      
         RETURN
  END IF
   W( 39 ) = JVS( 773 )
   W( 44 ) = JVS( 774 )
   W( 46 ) = JVS( 775 )
   W( 50 ) = JVS( 776 )
   W( 53 ) = JVS( 777 )
   W( 54 ) = JVS( 778 )
   W( 69 ) = JVS( 779 )
   W( 71 ) = JVS( 780 )
   W( 72 ) = JVS( 781 )
   W( 75 ) = JVS( 782 )
   W( 76 ) = JVS( 783 )
   W( 78 ) = JVS( 784 )
   W( 80 ) = JVS( 785 )
   W( 81 ) = JVS( 786 )
   W( 83 ) = JVS( 787 )
   W( 84 ) = JVS( 788 )
   W( 85 ) = JVS( 789 )
   W( 86 ) = JVS( 790 )
   W( 87 ) = JVS( 791 )
   W( 88 ) = JVS( 792 )
   W( 89 ) = JVS( 793 )
   W( 90 ) = JVS( 794 )
   W( 91 ) = JVS( 795 )
   W( 92 ) = JVS( 796 )
   W( 93 ) = JVS( 797 )
   W( 94 ) = JVS( 798 )
   W( 95 ) = JVS( 799 )
   W( 96 ) = JVS( 800 )
   W( 97 ) = JVS( 801 )
   W( 98 ) = JVS( 802 )
   W( 99 ) = JVS( 803 )
   W( 100 ) = JVS( 804 )
  a = -W( 39 ) / JVS(          213  )
  W( 39 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 214 )
  a = -W( 44 ) / JVS(          229  )
  W( 44 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 230 )
  a = -W( 46 ) / JVS(          235  )
  W( 46 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 236 )
  a = -W( 50 ) / JVS(          250  )
  W( 50 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 251 )
  a = -W( 53 ) / JVS(          261  )
  W( 53 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 262 )
  a = -W( 54 ) / JVS(          263  )
  W( 54 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 264 )
  a = -W( 69 ) / JVS(          388  )
  W( 69 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 389 )
  W( 88 ) = W( 88 ) + a*JVS( 390 )
  W( 93 ) = W( 93 ) + a*JVS( 391 )
  W( 98 ) = W( 98 ) + a*JVS( 392 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 83 ) / JVS(          573  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 574 )
  W( 88 ) = W( 88 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 98 ) = W( 98 ) + a*JVS( 579 )
  W( 100 ) = W( 100 ) + a*JVS( 580 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  a = -W( 85 ) / JVS(          616  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 617 )
  W( 87 ) = W( 87 ) + a*JVS( 618 )
  W( 88 ) = W( 88 ) + a*JVS( 619 )
  W( 90 ) = W( 90 ) + a*JVS( 620 )
  W( 91 ) = W( 91 ) + a*JVS( 621 )
  W( 92 ) = W( 92 ) + a*JVS( 622 )
  W( 93 ) = W( 93 ) + a*JVS( 623 )
  W( 97 ) = W( 97 ) + a*JVS( 624 )
  W( 98 ) = W( 98 ) + a*JVS( 625 )
  W( 100 ) = W( 100 ) + a*JVS( 626 )
  a = -W( 86 ) / JVS(          638  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 639 )
  W( 88 ) = W( 88 ) + a*JVS( 640 )
  W( 90 ) = W( 90 ) + a*JVS( 641 )
  W( 91 ) = W( 91 ) + a*JVS( 642 )
  W( 92 ) = W( 92 ) + a*JVS( 643 )
  W( 93 ) = W( 93 ) + a*JVS( 644 )
  W( 94 ) = W( 94 ) + a*JVS( 645 )
  W( 98 ) = W( 98 ) + a*JVS( 646 )
  W( 100 ) = W( 100 ) + a*JVS( 647 )
  a = -W( 87 ) / JVS(          660  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 661 )
  W( 89 ) = W( 89 ) + a*JVS( 662 )
  W( 90 ) = W( 90 ) + a*JVS( 663 )
  W( 91 ) = W( 91 ) + a*JVS( 664 )
  W( 92 ) = W( 92 ) + a*JVS( 665 )
  W( 93 ) = W( 93 ) + a*JVS( 666 )
  W( 94 ) = W( 94 ) + a*JVS( 667 )
  W( 95 ) = W( 95 ) + a*JVS( 668 )
  W( 98 ) = W( 98 ) + a*JVS( 669 )
  W( 99 ) = W( 99 ) + a*JVS( 670 )
  W( 100 ) = W( 100 ) + a*JVS( 671 )
  a = -W( 88 ) / JVS(          685  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 686 )
  W( 92 ) = W( 92 ) + a*JVS( 687 )
  W( 93 ) = W( 93 ) + a*JVS( 688 )
  W( 95 ) = W( 95 ) + a*JVS( 689 )
  W( 96 ) = W( 96 ) + a*JVS( 690 )
  W( 97 ) = W( 97 ) + a*JVS( 691 )
  W( 98 ) = W( 98 ) + a*JVS( 692 )
  W( 99 ) = W( 99 ) + a*JVS( 693 )
  W( 100 ) = W( 100 ) + a*JVS( 694 )
  a = -W( 89 ) / JVS(          717  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 718 )
  W( 91 ) = W( 91 ) + a*JVS( 719 )
  W( 92 ) = W( 92 ) + a*JVS( 720 )
  W( 93 ) = W( 93 ) + a*JVS( 721 )
  W( 94 ) = W( 94 ) + a*JVS( 722 )
  W( 95 ) = W( 95 ) + a*JVS( 723 )
  W( 96 ) = W( 96 ) + a*JVS( 724 )
  W( 97 ) = W( 97 ) + a*JVS( 725 )
  W( 98 ) = W( 98 ) + a*JVS( 726 )
  W( 99 ) = W( 99 ) + a*JVS( 727 )
  W( 100 ) = W( 100 ) + a*JVS( 728 )
  a = -W( 90 ) / JVS(          762  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 92 ) = W( 92 ) + a*JVS( 764 )
  W( 93 ) = W( 93 ) + a*JVS( 765 )
  W( 94 ) = W( 94 ) + a*JVS( 766 )
  W( 95 ) = W( 95 ) + a*JVS( 767 )
  W( 96 ) = W( 96 ) + a*JVS( 768 )
  W( 97 ) = W( 97 ) + a*JVS( 769 )
  W( 98 ) = W( 98 ) + a*JVS( 770 )
  W( 99 ) = W( 99 ) + a*JVS( 771 )
  W( 100 ) = W( 100 ) + a*JVS( 772 )
  JVS( 773) = W( 39 )
  JVS( 774) = W( 44 )
  JVS( 775) = W( 46 )
  JVS( 776) = W( 50 )
  JVS( 777) = W( 53 )
  JVS( 778) = W( 54 )
  JVS( 779) = W( 69 )
  JVS( 780) = W( 71 )
  JVS( 781) = W( 72 )
  JVS( 782) = W( 75 )
  JVS( 783) = W( 76 )
  JVS( 784) = W( 78 )
  JVS( 785) = W( 80 )
  JVS( 786) = W( 81 )
  JVS( 787) = W( 83 )
  JVS( 788) = W( 84 )
  JVS( 789) = W( 85 )
  JVS( 790) = W( 86 )
  JVS( 791) = W( 87 )
  JVS( 792) = W( 88 )
  JVS( 793) = W( 89 )
  JVS( 794) = W( 90 )
  JVS( 795) = W( 91 )
  JVS( 796) = W( 92 )
  JVS( 797) = W( 93 )
  JVS( 798) = W( 94 )
  JVS( 799) = W( 95 )
  JVS( 800) = W( 96 )
  JVS( 801) = W( 97 )
  JVS( 802) = W( 98 )
  JVS( 803) = W( 99 )
  JVS( 804) = W( 100 )
  IF ( ABS(  JVS( 838 )) < TINY(a) ) THEN
         IER = 92                                      
         RETURN
  END IF
   W( 34 ) = JVS( 805 )
   W( 35 ) = JVS( 806 )
   W( 36 ) = JVS( 807 )
   W( 37 ) = JVS( 808 )
   W( 42 ) = JVS( 809 )
   W( 43 ) = JVS( 810 )
   W( 45 ) = JVS( 811 )
   W( 48 ) = JVS( 812 )
   W( 49 ) = JVS( 813 )
   W( 51 ) = JVS( 814 )
   W( 61 ) = JVS( 815 )
   W( 66 ) = JVS( 816 )
   W( 70 ) = JVS( 817 )
   W( 71 ) = JVS( 818 )
   W( 73 ) = JVS( 819 )
   W( 74 ) = JVS( 820 )
   W( 75 ) = JVS( 821 )
   W( 76 ) = JVS( 822 )
   W( 77 ) = JVS( 823 )
   W( 78 ) = JVS( 824 )
   W( 79 ) = JVS( 825 )
   W( 80 ) = JVS( 826 )
   W( 81 ) = JVS( 827 )
   W( 82 ) = JVS( 828 )
   W( 83 ) = JVS( 829 )
   W( 84 ) = JVS( 830 )
   W( 85 ) = JVS( 831 )
   W( 86 ) = JVS( 832 )
   W( 87 ) = JVS( 833 )
   W( 88 ) = JVS( 834 )
   W( 89 ) = JVS( 835 )
   W( 90 ) = JVS( 836 )
   W( 91 ) = JVS( 837 )
   W( 92 ) = JVS( 838 )
   W( 93 ) = JVS( 839 )
   W( 94 ) = JVS( 840 )
   W( 95 ) = JVS( 841 )
   W( 96 ) = JVS( 842 )
   W( 97 ) = JVS( 843 )
   W( 98 ) = JVS( 844 )
   W( 99 ) = JVS( 845 )
   W( 100 ) = JVS( 846 )
  a = -W( 34 ) / JVS(          198  )
  W( 34 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 199 )
  W( 92 ) = W( 92 ) + a*JVS( 200 )
  a = -W( 35 ) / JVS(          201  )
  W( 35 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 202 )
  W( 99 ) = W( 99 ) + a*JVS( 203 )
  a = -W( 36 ) / JVS(          204  )
  W( 36 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 205 )
  W( 95 ) = W( 95 ) + a*JVS( 206 )
  a = -W( 37 ) / JVS(          207  )
  W( 37 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 208 )
  W( 96 ) = W( 96 ) + a*JVS( 209 )
  a = -W( 42 ) / JVS(          223  )
  W( 42 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 224 )
  W( 93 ) = W( 93 ) + a*JVS( 225 )
  a = -W( 43 ) / JVS(          226  )
  W( 43 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 227 )
  W( 100 ) = W( 100 ) + a*JVS( 228 )
  a = -W( 45 ) / JVS(          232  )
  W( 45 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 233 )
  W( 98 ) = W( 98 ) + a*JVS( 234 )
  a = -W( 48 ) / JVS(          241  )
  W( 48 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 242 )
  W( 97 ) = W( 97 ) + a*JVS( 243 )
  W( 100 ) = W( 100 ) + a*JVS( 244 )
  a = -W( 49 ) / JVS(          245  )
  W( 49 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 246 )
  W( 92 ) = W( 92 ) + a*JVS( 247 )
  W( 93 ) = W( 93 ) + a*JVS( 248 )
  W( 97 ) = W( 97 ) + a*JVS( 249 )
  a = -W( 51 ) / JVS(          252  )
  W( 51 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 253 )
  W( 97 ) = W( 97 ) + a*JVS( 254 )
  W( 98 ) = W( 98 ) + a*JVS( 255 )
  a = -W( 61 ) / JVS(          292  )
  W( 61 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 293 )
  W( 92 ) = W( 92 ) + a*JVS( 294 )
  W( 93 ) = W( 93 ) + a*JVS( 295 )
  W( 97 ) = W( 97 ) + a*JVS( 296 )
  a = -W( 66 ) / JVS(          353  )
  W( 66 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 354 )
  W( 74 ) = W( 74 ) + a*JVS( 355 )
  W( 77 ) = W( 77 ) + a*JVS( 356 )
  W( 78 ) = W( 78 ) + a*JVS( 357 )
  W( 79 ) = W( 79 ) + a*JVS( 358 )
  W( 80 ) = W( 80 ) + a*JVS( 359 )
  W( 81 ) = W( 81 ) + a*JVS( 360 )
  W( 82 ) = W( 82 ) + a*JVS( 361 )
  W( 85 ) = W( 85 ) + a*JVS( 362 )
  W( 88 ) = W( 88 ) + a*JVS( 363 )
  W( 92 ) = W( 92 ) + a*JVS( 364 )
  W( 93 ) = W( 93 ) + a*JVS( 365 )
  W( 97 ) = W( 97 ) + a*JVS( 366 )
  W( 98 ) = W( 98 ) + a*JVS( 367 )
  a = -W( 70 ) / JVS(          395  )
  W( 70 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 396 )
  W( 92 ) = W( 92 ) + a*JVS( 397 )
  W( 93 ) = W( 93 ) + a*JVS( 398 )
  W( 95 ) = W( 95 ) + a*JVS( 399 )
  W( 96 ) = W( 96 ) + a*JVS( 400 )
  W( 97 ) = W( 97 ) + a*JVS( 401 )
  W( 98 ) = W( 98 ) + a*JVS( 402 )
  W( 99 ) = W( 99 ) + a*JVS( 403 )
  W( 100 ) = W( 100 ) + a*JVS( 404 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  a = -W( 73 ) / JVS(          423  )
  W( 73 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 424 )
  W( 76 ) = W( 76 ) + a*JVS( 425 )
  W( 80 ) = W( 80 ) + a*JVS( 426 )
  W( 81 ) = W( 81 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 86 ) = W( 86 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 88 ) = W( 88 ) + a*JVS( 432 )
  W( 89 ) = W( 89 ) + a*JVS( 433 )
  W( 90 ) = W( 90 ) + a*JVS( 434 )
  W( 91 ) = W( 91 ) + a*JVS( 435 )
  W( 92 ) = W( 92 ) + a*JVS( 436 )
  W( 93 ) = W( 93 ) + a*JVS( 437 )
  W( 94 ) = W( 94 ) + a*JVS( 438 )
  W( 95 ) = W( 95 ) + a*JVS( 439 )
  W( 96 ) = W( 96 ) + a*JVS( 440 )
  W( 97 ) = W( 97 ) + a*JVS( 441 )
  W( 98 ) = W( 98 ) + a*JVS( 442 )
  W( 99 ) = W( 99 ) + a*JVS( 443 )
  W( 100 ) = W( 100 ) + a*JVS( 444 )
  a = -W( 74 ) / JVS(          446  )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 447 )
  W( 84 ) = W( 84 ) + a*JVS( 448 )
  W( 88 ) = W( 88 ) + a*JVS( 449 )
  W( 93 ) = W( 93 ) + a*JVS( 450 )
  W( 98 ) = W( 98 ) + a*JVS( 451 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 77 ) / JVS(          473  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 474 )
  W( 84 ) = W( 84 ) + a*JVS( 475 )
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 89 ) = W( 89 ) + a*JVS( 477 )
  W( 92 ) = W( 92 ) + a*JVS( 478 )
  W( 93 ) = W( 93 ) + a*JVS( 479 )
  W( 95 ) = W( 95 ) + a*JVS( 480 )
  W( 96 ) = W( 96 ) + a*JVS( 481 )
  W( 97 ) = W( 97 ) + a*JVS( 482 )
  W( 98 ) = W( 98 ) + a*JVS( 483 )
  W( 99 ) = W( 99 ) + a*JVS( 484 )
  W( 100 ) = W( 100 ) + a*JVS( 485 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 79 ) / JVS(          502  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 503 )
  W( 83 ) = W( 83 ) + a*JVS( 504 )
  W( 84 ) = W( 84 ) + a*JVS( 505 )
  W( 85 ) = W( 85 ) + a*JVS( 506 )
  W( 86 ) = W( 86 ) + a*JVS( 507 )
  W( 87 ) = W( 87 ) + a*JVS( 508 )
  W( 88 ) = W( 88 ) + a*JVS( 509 )
  W( 89 ) = W( 89 ) + a*JVS( 510 )
  W( 93 ) = W( 93 ) + a*JVS( 511 )
  W( 95 ) = W( 95 ) + a*JVS( 512 )
  W( 96 ) = W( 96 ) + a*JVS( 513 )
  W( 98 ) = W( 98 ) + a*JVS( 514 )
  W( 99 ) = W( 99 ) + a*JVS( 515 )
  W( 100 ) = W( 100 ) + a*JVS( 516 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 82 ) / JVS(          549  )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 550 )
  W( 84 ) = W( 84 ) + a*JVS( 551 )
  W( 86 ) = W( 86 ) + a*JVS( 552 )
  W( 87 ) = W( 87 ) + a*JVS( 553 )
  W( 88 ) = W( 88 ) + a*JVS( 554 )
  W( 89 ) = W( 89 ) + a*JVS( 555 )
  W( 90 ) = W( 90 ) + a*JVS( 556 )
  W( 91 ) = W( 91 ) + a*JVS( 557 )
  W( 92 ) = W( 92 ) + a*JVS( 558 )
  W( 93 ) = W( 93 ) + a*JVS( 559 )
  W( 94 ) = W( 94 ) + a*JVS( 560 )
  W( 95 ) = W( 95 ) + a*JVS( 561 )
  W( 96 ) = W( 96 ) + a*JVS( 562 )
  W( 97 ) = W( 97 ) + a*JVS( 563 )
  W( 98 ) = W( 98 ) + a*JVS( 564 )
  W( 99 ) = W( 99 ) + a*JVS( 565 )
  W( 100 ) = W( 100 ) + a*JVS( 566 )
  a = -W( 83 ) / JVS(          573  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 574 )
  W( 88 ) = W( 88 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 98 ) = W( 98 ) + a*JVS( 579 )
  W( 100 ) = W( 100 ) + a*JVS( 580 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  a = -W( 85 ) / JVS(          616  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 617 )
  W( 87 ) = W( 87 ) + a*JVS( 618 )
  W( 88 ) = W( 88 ) + a*JVS( 619 )
  W( 90 ) = W( 90 ) + a*JVS( 620 )
  W( 91 ) = W( 91 ) + a*JVS( 621 )
  W( 92 ) = W( 92 ) + a*JVS( 622 )
  W( 93 ) = W( 93 ) + a*JVS( 623 )
  W( 97 ) = W( 97 ) + a*JVS( 624 )
  W( 98 ) = W( 98 ) + a*JVS( 625 )
  W( 100 ) = W( 100 ) + a*JVS( 626 )
  a = -W( 86 ) / JVS(          638  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 639 )
  W( 88 ) = W( 88 ) + a*JVS( 640 )
  W( 90 ) = W( 90 ) + a*JVS( 641 )
  W( 91 ) = W( 91 ) + a*JVS( 642 )
  W( 92 ) = W( 92 ) + a*JVS( 643 )
  W( 93 ) = W( 93 ) + a*JVS( 644 )
  W( 94 ) = W( 94 ) + a*JVS( 645 )
  W( 98 ) = W( 98 ) + a*JVS( 646 )
  W( 100 ) = W( 100 ) + a*JVS( 647 )
  a = -W( 87 ) / JVS(          660  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 661 )
  W( 89 ) = W( 89 ) + a*JVS( 662 )
  W( 90 ) = W( 90 ) + a*JVS( 663 )
  W( 91 ) = W( 91 ) + a*JVS( 664 )
  W( 92 ) = W( 92 ) + a*JVS( 665 )
  W( 93 ) = W( 93 ) + a*JVS( 666 )
  W( 94 ) = W( 94 ) + a*JVS( 667 )
  W( 95 ) = W( 95 ) + a*JVS( 668 )
  W( 98 ) = W( 98 ) + a*JVS( 669 )
  W( 99 ) = W( 99 ) + a*JVS( 670 )
  W( 100 ) = W( 100 ) + a*JVS( 671 )
  a = -W( 88 ) / JVS(          685  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 686 )
  W( 92 ) = W( 92 ) + a*JVS( 687 )
  W( 93 ) = W( 93 ) + a*JVS( 688 )
  W( 95 ) = W( 95 ) + a*JVS( 689 )
  W( 96 ) = W( 96 ) + a*JVS( 690 )
  W( 97 ) = W( 97 ) + a*JVS( 691 )
  W( 98 ) = W( 98 ) + a*JVS( 692 )
  W( 99 ) = W( 99 ) + a*JVS( 693 )
  W( 100 ) = W( 100 ) + a*JVS( 694 )
  a = -W( 89 ) / JVS(          717  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 718 )
  W( 91 ) = W( 91 ) + a*JVS( 719 )
  W( 92 ) = W( 92 ) + a*JVS( 720 )
  W( 93 ) = W( 93 ) + a*JVS( 721 )
  W( 94 ) = W( 94 ) + a*JVS( 722 )
  W( 95 ) = W( 95 ) + a*JVS( 723 )
  W( 96 ) = W( 96 ) + a*JVS( 724 )
  W( 97 ) = W( 97 ) + a*JVS( 725 )
  W( 98 ) = W( 98 ) + a*JVS( 726 )
  W( 99 ) = W( 99 ) + a*JVS( 727 )
  W( 100 ) = W( 100 ) + a*JVS( 728 )
  a = -W( 90 ) / JVS(          762  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 92 ) = W( 92 ) + a*JVS( 764 )
  W( 93 ) = W( 93 ) + a*JVS( 765 )
  W( 94 ) = W( 94 ) + a*JVS( 766 )
  W( 95 ) = W( 95 ) + a*JVS( 767 )
  W( 96 ) = W( 96 ) + a*JVS( 768 )
  W( 97 ) = W( 97 ) + a*JVS( 769 )
  W( 98 ) = W( 98 ) + a*JVS( 770 )
  W( 99 ) = W( 99 ) + a*JVS( 771 )
  W( 100 ) = W( 100 ) + a*JVS( 772 )
  a = -W( 91 ) / JVS(          795  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 796 )
  W( 93 ) = W( 93 ) + a*JVS( 797 )
  W( 94 ) = W( 94 ) + a*JVS( 798 )
  W( 95 ) = W( 95 ) + a*JVS( 799 )
  W( 96 ) = W( 96 ) + a*JVS( 800 )
  W( 97 ) = W( 97 ) + a*JVS( 801 )
  W( 98 ) = W( 98 ) + a*JVS( 802 )
  W( 99 ) = W( 99 ) + a*JVS( 803 )
  W( 100 ) = W( 100 ) + a*JVS( 804 )
  JVS( 805) = W( 34 )
  JVS( 806) = W( 35 )
  JVS( 807) = W( 36 )
  JVS( 808) = W( 37 )
  JVS( 809) = W( 42 )
  JVS( 810) = W( 43 )
  JVS( 811) = W( 45 )
  JVS( 812) = W( 48 )
  JVS( 813) = W( 49 )
  JVS( 814) = W( 51 )
  JVS( 815) = W( 61 )
  JVS( 816) = W( 66 )
  JVS( 817) = W( 70 )
  JVS( 818) = W( 71 )
  JVS( 819) = W( 73 )
  JVS( 820) = W( 74 )
  JVS( 821) = W( 75 )
  JVS( 822) = W( 76 )
  JVS( 823) = W( 77 )
  JVS( 824) = W( 78 )
  JVS( 825) = W( 79 )
  JVS( 826) = W( 80 )
  JVS( 827) = W( 81 )
  JVS( 828) = W( 82 )
  JVS( 829) = W( 83 )
  JVS( 830) = W( 84 )
  JVS( 831) = W( 85 )
  JVS( 832) = W( 86 )
  JVS( 833) = W( 87 )
  JVS( 834) = W( 88 )
  JVS( 835) = W( 89 )
  JVS( 836) = W( 90 )
  JVS( 837) = W( 91 )
  JVS( 838) = W( 92 )
  JVS( 839) = W( 93 )
  JVS( 840) = W( 94 )
  JVS( 841) = W( 95 )
  JVS( 842) = W( 96 )
  JVS( 843) = W( 97 )
  JVS( 844) = W( 98 )
  JVS( 845) = W( 99 )
  JVS( 846) = W( 100 )
  IF ( ABS(  JVS( 880 )) < TINY(a) ) THEN
         IER = 93                                      
         RETURN
  END IF
   W( 42 ) = JVS( 847 )
   W( 51 ) = JVS( 848 )
   W( 57 ) = JVS( 849 )
   W( 61 ) = JVS( 850 )
   W( 62 ) = JVS( 851 )
   W( 63 ) = JVS( 852 )
   W( 64 ) = JVS( 853 )
   W( 66 ) = JVS( 854 )
   W( 67 ) = JVS( 855 )
   W( 69 ) = JVS( 856 )
   W( 70 ) = JVS( 857 )
   W( 71 ) = JVS( 858 )
   W( 72 ) = JVS( 859 )
   W( 73 ) = JVS( 860 )
   W( 74 ) = JVS( 861 )
   W( 75 ) = JVS( 862 )
   W( 76 ) = JVS( 863 )
   W( 77 ) = JVS( 864 )
   W( 78 ) = JVS( 865 )
   W( 79 ) = JVS( 866 )
   W( 80 ) = JVS( 867 )
   W( 81 ) = JVS( 868 )
   W( 82 ) = JVS( 869 )
   W( 83 ) = JVS( 870 )
   W( 84 ) = JVS( 871 )
   W( 85 ) = JVS( 872 )
   W( 86 ) = JVS( 873 )
   W( 87 ) = JVS( 874 )
   W( 88 ) = JVS( 875 )
   W( 89 ) = JVS( 876 )
   W( 90 ) = JVS( 877 )
   W( 91 ) = JVS( 878 )
   W( 92 ) = JVS( 879 )
   W( 93 ) = JVS( 880 )
   W( 94 ) = JVS( 881 )
   W( 95 ) = JVS( 882 )
   W( 96 ) = JVS( 883 )
   W( 97 ) = JVS( 884 )
   W( 98 ) = JVS( 885 )
   W( 99 ) = JVS( 886 )
   W( 100 ) = JVS( 887 )
  a = -W( 42 ) / JVS(          223  )
  W( 42 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 224 )
  W( 93 ) = W( 93 ) + a*JVS( 225 )
  a = -W( 51 ) / JVS(          252  )
  W( 51 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 253 )
  W( 97 ) = W( 97 ) + a*JVS( 254 )
  W( 98 ) = W( 98 ) + a*JVS( 255 )
  a = -W( 57 ) / JVS(          275  )
  W( 57 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 276 )
  W( 98 ) = W( 98 ) + a*JVS( 277 )
  a = -W( 61 ) / JVS(          292  )
  W( 61 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 293 )
  W( 92 ) = W( 92 ) + a*JVS( 294 )
  W( 93 ) = W( 93 ) + a*JVS( 295 )
  W( 97 ) = W( 97 ) + a*JVS( 296 )
  a = -W( 62 ) / JVS(          299  )
  W( 62 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 300 )
  W( 88 ) = W( 88 ) + a*JVS( 301 )
  W( 93 ) = W( 93 ) + a*JVS( 302 )
  W( 98 ) = W( 98 ) + a*JVS( 303 )
  a = -W( 63 ) / JVS(          305  )
  W( 63 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 306 )
  W( 93 ) = W( 93 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  W( 98 ) = W( 98 ) + a*JVS( 309 )
  a = -W( 64 ) / JVS(          315  )
  W( 64 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 316 )
  W( 78 ) = W( 78 ) + a*JVS( 317 )
  W( 81 ) = W( 81 ) + a*JVS( 318 )
  W( 88 ) = W( 88 ) + a*JVS( 319 )
  W( 93 ) = W( 93 ) + a*JVS( 320 )
  W( 98 ) = W( 98 ) + a*JVS( 321 )
  a = -W( 66 ) / JVS(          353  )
  W( 66 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 354 )
  W( 74 ) = W( 74 ) + a*JVS( 355 )
  W( 77 ) = W( 77 ) + a*JVS( 356 )
  W( 78 ) = W( 78 ) + a*JVS( 357 )
  W( 79 ) = W( 79 ) + a*JVS( 358 )
  W( 80 ) = W( 80 ) + a*JVS( 359 )
  W( 81 ) = W( 81 ) + a*JVS( 360 )
  W( 82 ) = W( 82 ) + a*JVS( 361 )
  W( 85 ) = W( 85 ) + a*JVS( 362 )
  W( 88 ) = W( 88 ) + a*JVS( 363 )
  W( 92 ) = W( 92 ) + a*JVS( 364 )
  W( 93 ) = W( 93 ) + a*JVS( 365 )
  W( 97 ) = W( 97 ) + a*JVS( 366 )
  W( 98 ) = W( 98 ) + a*JVS( 367 )
  a = -W( 67 ) / JVS(          368  )
  W( 67 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 369 )
  W( 88 ) = W( 88 ) + a*JVS( 370 )
  W( 93 ) = W( 93 ) + a*JVS( 371 )
  W( 98 ) = W( 98 ) + a*JVS( 372 )
  a = -W( 69 ) / JVS(          388  )
  W( 69 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 389 )
  W( 88 ) = W( 88 ) + a*JVS( 390 )
  W( 93 ) = W( 93 ) + a*JVS( 391 )
  W( 98 ) = W( 98 ) + a*JVS( 392 )
  a = -W( 70 ) / JVS(          395  )
  W( 70 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 396 )
  W( 92 ) = W( 92 ) + a*JVS( 397 )
  W( 93 ) = W( 93 ) + a*JVS( 398 )
  W( 95 ) = W( 95 ) + a*JVS( 399 )
  W( 96 ) = W( 96 ) + a*JVS( 400 )
  W( 97 ) = W( 97 ) + a*JVS( 401 )
  W( 98 ) = W( 98 ) + a*JVS( 402 )
  W( 99 ) = W( 99 ) + a*JVS( 403 )
  W( 100 ) = W( 100 ) + a*JVS( 404 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 73 ) / JVS(          423  )
  W( 73 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 424 )
  W( 76 ) = W( 76 ) + a*JVS( 425 )
  W( 80 ) = W( 80 ) + a*JVS( 426 )
  W( 81 ) = W( 81 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 86 ) = W( 86 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 88 ) = W( 88 ) + a*JVS( 432 )
  W( 89 ) = W( 89 ) + a*JVS( 433 )
  W( 90 ) = W( 90 ) + a*JVS( 434 )
  W( 91 ) = W( 91 ) + a*JVS( 435 )
  W( 92 ) = W( 92 ) + a*JVS( 436 )
  W( 93 ) = W( 93 ) + a*JVS( 437 )
  W( 94 ) = W( 94 ) + a*JVS( 438 )
  W( 95 ) = W( 95 ) + a*JVS( 439 )
  W( 96 ) = W( 96 ) + a*JVS( 440 )
  W( 97 ) = W( 97 ) + a*JVS( 441 )
  W( 98 ) = W( 98 ) + a*JVS( 442 )
  W( 99 ) = W( 99 ) + a*JVS( 443 )
  W( 100 ) = W( 100 ) + a*JVS( 444 )
  a = -W( 74 ) / JVS(          446  )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 447 )
  W( 84 ) = W( 84 ) + a*JVS( 448 )
  W( 88 ) = W( 88 ) + a*JVS( 449 )
  W( 93 ) = W( 93 ) + a*JVS( 450 )
  W( 98 ) = W( 98 ) + a*JVS( 451 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 77 ) / JVS(          473  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 474 )
  W( 84 ) = W( 84 ) + a*JVS( 475 )
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 89 ) = W( 89 ) + a*JVS( 477 )
  W( 92 ) = W( 92 ) + a*JVS( 478 )
  W( 93 ) = W( 93 ) + a*JVS( 479 )
  W( 95 ) = W( 95 ) + a*JVS( 480 )
  W( 96 ) = W( 96 ) + a*JVS( 481 )
  W( 97 ) = W( 97 ) + a*JVS( 482 )
  W( 98 ) = W( 98 ) + a*JVS( 483 )
  W( 99 ) = W( 99 ) + a*JVS( 484 )
  W( 100 ) = W( 100 ) + a*JVS( 485 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 79 ) / JVS(          502  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 503 )
  W( 83 ) = W( 83 ) + a*JVS( 504 )
  W( 84 ) = W( 84 ) + a*JVS( 505 )
  W( 85 ) = W( 85 ) + a*JVS( 506 )
  W( 86 ) = W( 86 ) + a*JVS( 507 )
  W( 87 ) = W( 87 ) + a*JVS( 508 )
  W( 88 ) = W( 88 ) + a*JVS( 509 )
  W( 89 ) = W( 89 ) + a*JVS( 510 )
  W( 93 ) = W( 93 ) + a*JVS( 511 )
  W( 95 ) = W( 95 ) + a*JVS( 512 )
  W( 96 ) = W( 96 ) + a*JVS( 513 )
  W( 98 ) = W( 98 ) + a*JVS( 514 )
  W( 99 ) = W( 99 ) + a*JVS( 515 )
  W( 100 ) = W( 100 ) + a*JVS( 516 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 82 ) / JVS(          549  )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 550 )
  W( 84 ) = W( 84 ) + a*JVS( 551 )
  W( 86 ) = W( 86 ) + a*JVS( 552 )
  W( 87 ) = W( 87 ) + a*JVS( 553 )
  W( 88 ) = W( 88 ) + a*JVS( 554 )
  W( 89 ) = W( 89 ) + a*JVS( 555 )
  W( 90 ) = W( 90 ) + a*JVS( 556 )
  W( 91 ) = W( 91 ) + a*JVS( 557 )
  W( 92 ) = W( 92 ) + a*JVS( 558 )
  W( 93 ) = W( 93 ) + a*JVS( 559 )
  W( 94 ) = W( 94 ) + a*JVS( 560 )
  W( 95 ) = W( 95 ) + a*JVS( 561 )
  W( 96 ) = W( 96 ) + a*JVS( 562 )
  W( 97 ) = W( 97 ) + a*JVS( 563 )
  W( 98 ) = W( 98 ) + a*JVS( 564 )
  W( 99 ) = W( 99 ) + a*JVS( 565 )
  W( 100 ) = W( 100 ) + a*JVS( 566 )
  a = -W( 83 ) / JVS(          573  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 574 )
  W( 88 ) = W( 88 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 98 ) = W( 98 ) + a*JVS( 579 )
  W( 100 ) = W( 100 ) + a*JVS( 580 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  a = -W( 85 ) / JVS(          616  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 617 )
  W( 87 ) = W( 87 ) + a*JVS( 618 )
  W( 88 ) = W( 88 ) + a*JVS( 619 )
  W( 90 ) = W( 90 ) + a*JVS( 620 )
  W( 91 ) = W( 91 ) + a*JVS( 621 )
  W( 92 ) = W( 92 ) + a*JVS( 622 )
  W( 93 ) = W( 93 ) + a*JVS( 623 )
  W( 97 ) = W( 97 ) + a*JVS( 624 )
  W( 98 ) = W( 98 ) + a*JVS( 625 )
  W( 100 ) = W( 100 ) + a*JVS( 626 )
  a = -W( 86 ) / JVS(          638  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 639 )
  W( 88 ) = W( 88 ) + a*JVS( 640 )
  W( 90 ) = W( 90 ) + a*JVS( 641 )
  W( 91 ) = W( 91 ) + a*JVS( 642 )
  W( 92 ) = W( 92 ) + a*JVS( 643 )
  W( 93 ) = W( 93 ) + a*JVS( 644 )
  W( 94 ) = W( 94 ) + a*JVS( 645 )
  W( 98 ) = W( 98 ) + a*JVS( 646 )
  W( 100 ) = W( 100 ) + a*JVS( 647 )
  a = -W( 87 ) / JVS(          660  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 661 )
  W( 89 ) = W( 89 ) + a*JVS( 662 )
  W( 90 ) = W( 90 ) + a*JVS( 663 )
  W( 91 ) = W( 91 ) + a*JVS( 664 )
  W( 92 ) = W( 92 ) + a*JVS( 665 )
  W( 93 ) = W( 93 ) + a*JVS( 666 )
  W( 94 ) = W( 94 ) + a*JVS( 667 )
  W( 95 ) = W( 95 ) + a*JVS( 668 )
  W( 98 ) = W( 98 ) + a*JVS( 669 )
  W( 99 ) = W( 99 ) + a*JVS( 670 )
  W( 100 ) = W( 100 ) + a*JVS( 671 )
  a = -W( 88 ) / JVS(          685  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 686 )
  W( 92 ) = W( 92 ) + a*JVS( 687 )
  W( 93 ) = W( 93 ) + a*JVS( 688 )
  W( 95 ) = W( 95 ) + a*JVS( 689 )
  W( 96 ) = W( 96 ) + a*JVS( 690 )
  W( 97 ) = W( 97 ) + a*JVS( 691 )
  W( 98 ) = W( 98 ) + a*JVS( 692 )
  W( 99 ) = W( 99 ) + a*JVS( 693 )
  W( 100 ) = W( 100 ) + a*JVS( 694 )
  a = -W( 89 ) / JVS(          717  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 718 )
  W( 91 ) = W( 91 ) + a*JVS( 719 )
  W( 92 ) = W( 92 ) + a*JVS( 720 )
  W( 93 ) = W( 93 ) + a*JVS( 721 )
  W( 94 ) = W( 94 ) + a*JVS( 722 )
  W( 95 ) = W( 95 ) + a*JVS( 723 )
  W( 96 ) = W( 96 ) + a*JVS( 724 )
  W( 97 ) = W( 97 ) + a*JVS( 725 )
  W( 98 ) = W( 98 ) + a*JVS( 726 )
  W( 99 ) = W( 99 ) + a*JVS( 727 )
  W( 100 ) = W( 100 ) + a*JVS( 728 )
  a = -W( 90 ) / JVS(          762  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 92 ) = W( 92 ) + a*JVS( 764 )
  W( 93 ) = W( 93 ) + a*JVS( 765 )
  W( 94 ) = W( 94 ) + a*JVS( 766 )
  W( 95 ) = W( 95 ) + a*JVS( 767 )
  W( 96 ) = W( 96 ) + a*JVS( 768 )
  W( 97 ) = W( 97 ) + a*JVS( 769 )
  W( 98 ) = W( 98 ) + a*JVS( 770 )
  W( 99 ) = W( 99 ) + a*JVS( 771 )
  W( 100 ) = W( 100 ) + a*JVS( 772 )
  a = -W( 91 ) / JVS(          795  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 796 )
  W( 93 ) = W( 93 ) + a*JVS( 797 )
  W( 94 ) = W( 94 ) + a*JVS( 798 )
  W( 95 ) = W( 95 ) + a*JVS( 799 )
  W( 96 ) = W( 96 ) + a*JVS( 800 )
  W( 97 ) = W( 97 ) + a*JVS( 801 )
  W( 98 ) = W( 98 ) + a*JVS( 802 )
  W( 99 ) = W( 99 ) + a*JVS( 803 )
  W( 100 ) = W( 100 ) + a*JVS( 804 )
  a = -W( 92 ) / JVS(          838  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 839 )
  W( 94 ) = W( 94 ) + a*JVS( 840 )
  W( 95 ) = W( 95 ) + a*JVS( 841 )
  W( 96 ) = W( 96 ) + a*JVS( 842 )
  W( 97 ) = W( 97 ) + a*JVS( 843 )
  W( 98 ) = W( 98 ) + a*JVS( 844 )
  W( 99 ) = W( 99 ) + a*JVS( 845 )
  W( 100 ) = W( 100 ) + a*JVS( 846 )
  JVS( 847) = W( 42 )
  JVS( 848) = W( 51 )
  JVS( 849) = W( 57 )
  JVS( 850) = W( 61 )
  JVS( 851) = W( 62 )
  JVS( 852) = W( 63 )
  JVS( 853) = W( 64 )
  JVS( 854) = W( 66 )
  JVS( 855) = W( 67 )
  JVS( 856) = W( 69 )
  JVS( 857) = W( 70 )
  JVS( 858) = W( 71 )
  JVS( 859) = W( 72 )
  JVS( 860) = W( 73 )
  JVS( 861) = W( 74 )
  JVS( 862) = W( 75 )
  JVS( 863) = W( 76 )
  JVS( 864) = W( 77 )
  JVS( 865) = W( 78 )
  JVS( 866) = W( 79 )
  JVS( 867) = W( 80 )
  JVS( 868) = W( 81 )
  JVS( 869) = W( 82 )
  JVS( 870) = W( 83 )
  JVS( 871) = W( 84 )
  JVS( 872) = W( 85 )
  JVS( 873) = W( 86 )
  JVS( 874) = W( 87 )
  JVS( 875) = W( 88 )
  JVS( 876) = W( 89 )
  JVS( 877) = W( 90 )
  JVS( 878) = W( 91 )
  JVS( 879) = W( 92 )
  JVS( 880) = W( 93 )
  JVS( 881) = W( 94 )
  JVS( 882) = W( 95 )
  JVS( 883) = W( 96 )
  JVS( 884) = W( 97 )
  JVS( 885) = W( 98 )
  JVS( 886) = W( 99 )
  JVS( 887) = W( 100 )
  IF ( ABS(  JVS( 913 )) < TINY(a) ) THEN
         IER = 94                                      
         RETURN
  END IF
   W( 31 ) = JVS( 888 )
   W( 45 ) = JVS( 889 )
   W( 47 ) = JVS( 890 )
   W( 53 ) = JVS( 891 )
   W( 67 ) = JVS( 892 )
   W( 68 ) = JVS( 893 )
   W( 69 ) = JVS( 894 )
   W( 71 ) = JVS( 895 )
   W( 72 ) = JVS( 896 )
   W( 75 ) = JVS( 897 )
   W( 76 ) = JVS( 898 )
   W( 79 ) = JVS( 899 )
   W( 80 ) = JVS( 900 )
   W( 81 ) = JVS( 901 )
   W( 83 ) = JVS( 902 )
   W( 84 ) = JVS( 903 )
   W( 85 ) = JVS( 904 )
   W( 86 ) = JVS( 905 )
   W( 87 ) = JVS( 906 )
   W( 88 ) = JVS( 907 )
   W( 89 ) = JVS( 908 )
   W( 90 ) = JVS( 909 )
   W( 91 ) = JVS( 910 )
   W( 92 ) = JVS( 911 )
   W( 93 ) = JVS( 912 )
   W( 94 ) = JVS( 913 )
   W( 95 ) = JVS( 914 )
   W( 96 ) = JVS( 915 )
   W( 97 ) = JVS( 916 )
   W( 98 ) = JVS( 917 )
   W( 99 ) = JVS( 918 )
   W( 100 ) = JVS( 919 )
  a = -W( 31 ) / JVS(          192  )
  W( 31 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 193 )
  a = -W( 45 ) / JVS(          232  )
  W( 45 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 233 )
  W( 98 ) = W( 98 ) + a*JVS( 234 )
  a = -W( 47 ) / JVS(          237  )
  W( 47 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 238 )
  W( 97 ) = W( 97 ) + a*JVS( 239 )
  W( 98 ) = W( 98 ) + a*JVS( 240 )
  a = -W( 53 ) / JVS(          261  )
  W( 53 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 262 )
  a = -W( 67 ) / JVS(          368  )
  W( 67 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 369 )
  W( 88 ) = W( 88 ) + a*JVS( 370 )
  W( 93 ) = W( 93 ) + a*JVS( 371 )
  W( 98 ) = W( 98 ) + a*JVS( 372 )
  a = -W( 68 ) / JVS(          378  )
  W( 68 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 80 ) = W( 80 ) + a*JVS( 382 )
  W( 83 ) = W( 83 ) + a*JVS( 383 )
  W( 88 ) = W( 88 ) + a*JVS( 384 )
  W( 92 ) = W( 92 ) + a*JVS( 385 )
  W( 93 ) = W( 93 ) + a*JVS( 386 )
  W( 98 ) = W( 98 ) + a*JVS( 387 )
  a = -W( 69 ) / JVS(          388  )
  W( 69 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 389 )
  W( 88 ) = W( 88 ) + a*JVS( 390 )
  W( 93 ) = W( 93 ) + a*JVS( 391 )
  W( 98 ) = W( 98 ) + a*JVS( 392 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 79 ) / JVS(          502  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 503 )
  W( 83 ) = W( 83 ) + a*JVS( 504 )
  W( 84 ) = W( 84 ) + a*JVS( 505 )
  W( 85 ) = W( 85 ) + a*JVS( 506 )
  W( 86 ) = W( 86 ) + a*JVS( 507 )
  W( 87 ) = W( 87 ) + a*JVS( 508 )
  W( 88 ) = W( 88 ) + a*JVS( 509 )
  W( 89 ) = W( 89 ) + a*JVS( 510 )
  W( 93 ) = W( 93 ) + a*JVS( 511 )
  W( 95 ) = W( 95 ) + a*JVS( 512 )
  W( 96 ) = W( 96 ) + a*JVS( 513 )
  W( 98 ) = W( 98 ) + a*JVS( 514 )
  W( 99 ) = W( 99 ) + a*JVS( 515 )
  W( 100 ) = W( 100 ) + a*JVS( 516 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 83 ) / JVS(          573  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 574 )
  W( 88 ) = W( 88 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 98 ) = W( 98 ) + a*JVS( 579 )
  W( 100 ) = W( 100 ) + a*JVS( 580 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  a = -W( 85 ) / JVS(          616  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 617 )
  W( 87 ) = W( 87 ) + a*JVS( 618 )
  W( 88 ) = W( 88 ) + a*JVS( 619 )
  W( 90 ) = W( 90 ) + a*JVS( 620 )
  W( 91 ) = W( 91 ) + a*JVS( 621 )
  W( 92 ) = W( 92 ) + a*JVS( 622 )
  W( 93 ) = W( 93 ) + a*JVS( 623 )
  W( 97 ) = W( 97 ) + a*JVS( 624 )
  W( 98 ) = W( 98 ) + a*JVS( 625 )
  W( 100 ) = W( 100 ) + a*JVS( 626 )
  a = -W( 86 ) / JVS(          638  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 639 )
  W( 88 ) = W( 88 ) + a*JVS( 640 )
  W( 90 ) = W( 90 ) + a*JVS( 641 )
  W( 91 ) = W( 91 ) + a*JVS( 642 )
  W( 92 ) = W( 92 ) + a*JVS( 643 )
  W( 93 ) = W( 93 ) + a*JVS( 644 )
  W( 94 ) = W( 94 ) + a*JVS( 645 )
  W( 98 ) = W( 98 ) + a*JVS( 646 )
  W( 100 ) = W( 100 ) + a*JVS( 647 )
  a = -W( 87 ) / JVS(          660  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 661 )
  W( 89 ) = W( 89 ) + a*JVS( 662 )
  W( 90 ) = W( 90 ) + a*JVS( 663 )
  W( 91 ) = W( 91 ) + a*JVS( 664 )
  W( 92 ) = W( 92 ) + a*JVS( 665 )
  W( 93 ) = W( 93 ) + a*JVS( 666 )
  W( 94 ) = W( 94 ) + a*JVS( 667 )
  W( 95 ) = W( 95 ) + a*JVS( 668 )
  W( 98 ) = W( 98 ) + a*JVS( 669 )
  W( 99 ) = W( 99 ) + a*JVS( 670 )
  W( 100 ) = W( 100 ) + a*JVS( 671 )
  a = -W( 88 ) / JVS(          685  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 686 )
  W( 92 ) = W( 92 ) + a*JVS( 687 )
  W( 93 ) = W( 93 ) + a*JVS( 688 )
  W( 95 ) = W( 95 ) + a*JVS( 689 )
  W( 96 ) = W( 96 ) + a*JVS( 690 )
  W( 97 ) = W( 97 ) + a*JVS( 691 )
  W( 98 ) = W( 98 ) + a*JVS( 692 )
  W( 99 ) = W( 99 ) + a*JVS( 693 )
  W( 100 ) = W( 100 ) + a*JVS( 694 )
  a = -W( 89 ) / JVS(          717  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 718 )
  W( 91 ) = W( 91 ) + a*JVS( 719 )
  W( 92 ) = W( 92 ) + a*JVS( 720 )
  W( 93 ) = W( 93 ) + a*JVS( 721 )
  W( 94 ) = W( 94 ) + a*JVS( 722 )
  W( 95 ) = W( 95 ) + a*JVS( 723 )
  W( 96 ) = W( 96 ) + a*JVS( 724 )
  W( 97 ) = W( 97 ) + a*JVS( 725 )
  W( 98 ) = W( 98 ) + a*JVS( 726 )
  W( 99 ) = W( 99 ) + a*JVS( 727 )
  W( 100 ) = W( 100 ) + a*JVS( 728 )
  a = -W( 90 ) / JVS(          762  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 92 ) = W( 92 ) + a*JVS( 764 )
  W( 93 ) = W( 93 ) + a*JVS( 765 )
  W( 94 ) = W( 94 ) + a*JVS( 766 )
  W( 95 ) = W( 95 ) + a*JVS( 767 )
  W( 96 ) = W( 96 ) + a*JVS( 768 )
  W( 97 ) = W( 97 ) + a*JVS( 769 )
  W( 98 ) = W( 98 ) + a*JVS( 770 )
  W( 99 ) = W( 99 ) + a*JVS( 771 )
  W( 100 ) = W( 100 ) + a*JVS( 772 )
  a = -W( 91 ) / JVS(          795  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 796 )
  W( 93 ) = W( 93 ) + a*JVS( 797 )
  W( 94 ) = W( 94 ) + a*JVS( 798 )
  W( 95 ) = W( 95 ) + a*JVS( 799 )
  W( 96 ) = W( 96 ) + a*JVS( 800 )
  W( 97 ) = W( 97 ) + a*JVS( 801 )
  W( 98 ) = W( 98 ) + a*JVS( 802 )
  W( 99 ) = W( 99 ) + a*JVS( 803 )
  W( 100 ) = W( 100 ) + a*JVS( 804 )
  a = -W( 92 ) / JVS(          838  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 839 )
  W( 94 ) = W( 94 ) + a*JVS( 840 )
  W( 95 ) = W( 95 ) + a*JVS( 841 )
  W( 96 ) = W( 96 ) + a*JVS( 842 )
  W( 97 ) = W( 97 ) + a*JVS( 843 )
  W( 98 ) = W( 98 ) + a*JVS( 844 )
  W( 99 ) = W( 99 ) + a*JVS( 845 )
  W( 100 ) = W( 100 ) + a*JVS( 846 )
  a = -W( 93 ) / JVS(          880  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 881 )
  W( 95 ) = W( 95 ) + a*JVS( 882 )
  W( 96 ) = W( 96 ) + a*JVS( 883 )
  W( 97 ) = W( 97 ) + a*JVS( 884 )
  W( 98 ) = W( 98 ) + a*JVS( 885 )
  W( 99 ) = W( 99 ) + a*JVS( 886 )
  W( 100 ) = W( 100 ) + a*JVS( 887 )
  JVS( 888) = W( 31 )
  JVS( 889) = W( 45 )
  JVS( 890) = W( 47 )
  JVS( 891) = W( 53 )
  JVS( 892) = W( 67 )
  JVS( 893) = W( 68 )
  JVS( 894) = W( 69 )
  JVS( 895) = W( 71 )
  JVS( 896) = W( 72 )
  JVS( 897) = W( 75 )
  JVS( 898) = W( 76 )
  JVS( 899) = W( 79 )
  JVS( 900) = W( 80 )
  JVS( 901) = W( 81 )
  JVS( 902) = W( 83 )
  JVS( 903) = W( 84 )
  JVS( 904) = W( 85 )
  JVS( 905) = W( 86 )
  JVS( 906) = W( 87 )
  JVS( 907) = W( 88 )
  JVS( 908) = W( 89 )
  JVS( 909) = W( 90 )
  JVS( 910) = W( 91 )
  JVS( 911) = W( 92 )
  JVS( 912) = W( 93 )
  JVS( 913) = W( 94 )
  JVS( 914) = W( 95 )
  JVS( 915) = W( 96 )
  JVS( 916) = W( 97 )
  JVS( 917) = W( 98 )
  JVS( 918) = W( 99 )
  JVS( 919) = W( 100 )
  IF ( ABS(  JVS( 931 )) < TINY(a) ) THEN
         IER = 95                                      
         RETURN
  END IF
   W( 36 ) = JVS( 920 )
   W( 62 ) = JVS( 921 )
   W( 80 ) = JVS( 922 )
   W( 84 ) = JVS( 923 )
   W( 88 ) = JVS( 924 )
   W( 89 ) = JVS( 925 )
   W( 90 ) = JVS( 926 )
   W( 91 ) = JVS( 927 )
   W( 92 ) = JVS( 928 )
   W( 93 ) = JVS( 929 )
   W( 94 ) = JVS( 930 )
   W( 95 ) = JVS( 931 )
   W( 96 ) = JVS( 932 )
   W( 97 ) = JVS( 933 )
   W( 98 ) = JVS( 934 )
   W( 99 ) = JVS( 935 )
   W( 100 ) = JVS( 936 )
  a = -W( 36 ) / JVS(          204  )
  W( 36 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 205 )
  W( 95 ) = W( 95 ) + a*JVS( 206 )
  a = -W( 62 ) / JVS(          299  )
  W( 62 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 300 )
  W( 88 ) = W( 88 ) + a*JVS( 301 )
  W( 93 ) = W( 93 ) + a*JVS( 302 )
  W( 98 ) = W( 98 ) + a*JVS( 303 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  a = -W( 88 ) / JVS(          685  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 686 )
  W( 92 ) = W( 92 ) + a*JVS( 687 )
  W( 93 ) = W( 93 ) + a*JVS( 688 )
  W( 95 ) = W( 95 ) + a*JVS( 689 )
  W( 96 ) = W( 96 ) + a*JVS( 690 )
  W( 97 ) = W( 97 ) + a*JVS( 691 )
  W( 98 ) = W( 98 ) + a*JVS( 692 )
  W( 99 ) = W( 99 ) + a*JVS( 693 )
  W( 100 ) = W( 100 ) + a*JVS( 694 )
  a = -W( 89 ) / JVS(          717  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 718 )
  W( 91 ) = W( 91 ) + a*JVS( 719 )
  W( 92 ) = W( 92 ) + a*JVS( 720 )
  W( 93 ) = W( 93 ) + a*JVS( 721 )
  W( 94 ) = W( 94 ) + a*JVS( 722 )
  W( 95 ) = W( 95 ) + a*JVS( 723 )
  W( 96 ) = W( 96 ) + a*JVS( 724 )
  W( 97 ) = W( 97 ) + a*JVS( 725 )
  W( 98 ) = W( 98 ) + a*JVS( 726 )
  W( 99 ) = W( 99 ) + a*JVS( 727 )
  W( 100 ) = W( 100 ) + a*JVS( 728 )
  a = -W( 90 ) / JVS(          762  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 92 ) = W( 92 ) + a*JVS( 764 )
  W( 93 ) = W( 93 ) + a*JVS( 765 )
  W( 94 ) = W( 94 ) + a*JVS( 766 )
  W( 95 ) = W( 95 ) + a*JVS( 767 )
  W( 96 ) = W( 96 ) + a*JVS( 768 )
  W( 97 ) = W( 97 ) + a*JVS( 769 )
  W( 98 ) = W( 98 ) + a*JVS( 770 )
  W( 99 ) = W( 99 ) + a*JVS( 771 )
  W( 100 ) = W( 100 ) + a*JVS( 772 )
  a = -W( 91 ) / JVS(          795  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 796 )
  W( 93 ) = W( 93 ) + a*JVS( 797 )
  W( 94 ) = W( 94 ) + a*JVS( 798 )
  W( 95 ) = W( 95 ) + a*JVS( 799 )
  W( 96 ) = W( 96 ) + a*JVS( 800 )
  W( 97 ) = W( 97 ) + a*JVS( 801 )
  W( 98 ) = W( 98 ) + a*JVS( 802 )
  W( 99 ) = W( 99 ) + a*JVS( 803 )
  W( 100 ) = W( 100 ) + a*JVS( 804 )
  a = -W( 92 ) / JVS(          838  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 839 )
  W( 94 ) = W( 94 ) + a*JVS( 840 )
  W( 95 ) = W( 95 ) + a*JVS( 841 )
  W( 96 ) = W( 96 ) + a*JVS( 842 )
  W( 97 ) = W( 97 ) + a*JVS( 843 )
  W( 98 ) = W( 98 ) + a*JVS( 844 )
  W( 99 ) = W( 99 ) + a*JVS( 845 )
  W( 100 ) = W( 100 ) + a*JVS( 846 )
  a = -W( 93 ) / JVS(          880  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 881 )
  W( 95 ) = W( 95 ) + a*JVS( 882 )
  W( 96 ) = W( 96 ) + a*JVS( 883 )
  W( 97 ) = W( 97 ) + a*JVS( 884 )
  W( 98 ) = W( 98 ) + a*JVS( 885 )
  W( 99 ) = W( 99 ) + a*JVS( 886 )
  W( 100 ) = W( 100 ) + a*JVS( 887 )
  a = -W( 94 ) / JVS(          913  )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 914 )
  W( 96 ) = W( 96 ) + a*JVS( 915 )
  W( 97 ) = W( 97 ) + a*JVS( 916 )
  W( 98 ) = W( 98 ) + a*JVS( 917 )
  W( 99 ) = W( 99 ) + a*JVS( 918 )
  W( 100 ) = W( 100 ) + a*JVS( 919 )
  JVS( 920) = W( 36 )
  JVS( 921) = W( 62 )
  JVS( 922) = W( 80 )
  JVS( 923) = W( 84 )
  JVS( 924) = W( 88 )
  JVS( 925) = W( 89 )
  JVS( 926) = W( 90 )
  JVS( 927) = W( 91 )
  JVS( 928) = W( 92 )
  JVS( 929) = W( 93 )
  JVS( 930) = W( 94 )
  JVS( 931) = W( 95 )
  JVS( 932) = W( 96 )
  JVS( 933) = W( 97 )
  JVS( 934) = W( 98 )
  JVS( 935) = W( 99 )
  JVS( 936) = W( 100 )
  IF ( ABS(  JVS( 952 )) < TINY(a) ) THEN
         IER = 96                                      
         RETURN
  END IF
   W( 37 ) = JVS( 937 )
   W( 71 ) = JVS( 938 )
   W( 74 ) = JVS( 939 )
   W( 78 ) = JVS( 940 )
   W( 80 ) = JVS( 941 )
   W( 81 ) = JVS( 942 )
   W( 84 ) = JVS( 943 )
   W( 88 ) = JVS( 944 )
   W( 89 ) = JVS( 945 )
   W( 90 ) = JVS( 946 )
   W( 91 ) = JVS( 947 )
   W( 92 ) = JVS( 948 )
   W( 93 ) = JVS( 949 )
   W( 94 ) = JVS( 950 )
   W( 95 ) = JVS( 951 )
   W( 96 ) = JVS( 952 )
   W( 97 ) = JVS( 953 )
   W( 98 ) = JVS( 954 )
   W( 99 ) = JVS( 955 )
   W( 100 ) = JVS( 956 )
  a = -W( 37 ) / JVS(          207  )
  W( 37 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 208 )
  W( 96 ) = W( 96 ) + a*JVS( 209 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  a = -W( 74 ) / JVS(          446  )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 447 )
  W( 84 ) = W( 84 ) + a*JVS( 448 )
  W( 88 ) = W( 88 ) + a*JVS( 449 )
  W( 93 ) = W( 93 ) + a*JVS( 450 )
  W( 98 ) = W( 98 ) + a*JVS( 451 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  a = -W( 88 ) / JVS(          685  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 686 )
  W( 92 ) = W( 92 ) + a*JVS( 687 )
  W( 93 ) = W( 93 ) + a*JVS( 688 )
  W( 95 ) = W( 95 ) + a*JVS( 689 )
  W( 96 ) = W( 96 ) + a*JVS( 690 )
  W( 97 ) = W( 97 ) + a*JVS( 691 )
  W( 98 ) = W( 98 ) + a*JVS( 692 )
  W( 99 ) = W( 99 ) + a*JVS( 693 )
  W( 100 ) = W( 100 ) + a*JVS( 694 )
  a = -W( 89 ) / JVS(          717  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 718 )
  W( 91 ) = W( 91 ) + a*JVS( 719 )
  W( 92 ) = W( 92 ) + a*JVS( 720 )
  W( 93 ) = W( 93 ) + a*JVS( 721 )
  W( 94 ) = W( 94 ) + a*JVS( 722 )
  W( 95 ) = W( 95 ) + a*JVS( 723 )
  W( 96 ) = W( 96 ) + a*JVS( 724 )
  W( 97 ) = W( 97 ) + a*JVS( 725 )
  W( 98 ) = W( 98 ) + a*JVS( 726 )
  W( 99 ) = W( 99 ) + a*JVS( 727 )
  W( 100 ) = W( 100 ) + a*JVS( 728 )
  a = -W( 90 ) / JVS(          762  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 92 ) = W( 92 ) + a*JVS( 764 )
  W( 93 ) = W( 93 ) + a*JVS( 765 )
  W( 94 ) = W( 94 ) + a*JVS( 766 )
  W( 95 ) = W( 95 ) + a*JVS( 767 )
  W( 96 ) = W( 96 ) + a*JVS( 768 )
  W( 97 ) = W( 97 ) + a*JVS( 769 )
  W( 98 ) = W( 98 ) + a*JVS( 770 )
  W( 99 ) = W( 99 ) + a*JVS( 771 )
  W( 100 ) = W( 100 ) + a*JVS( 772 )
  a = -W( 91 ) / JVS(          795  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 796 )
  W( 93 ) = W( 93 ) + a*JVS( 797 )
  W( 94 ) = W( 94 ) + a*JVS( 798 )
  W( 95 ) = W( 95 ) + a*JVS( 799 )
  W( 96 ) = W( 96 ) + a*JVS( 800 )
  W( 97 ) = W( 97 ) + a*JVS( 801 )
  W( 98 ) = W( 98 ) + a*JVS( 802 )
  W( 99 ) = W( 99 ) + a*JVS( 803 )
  W( 100 ) = W( 100 ) + a*JVS( 804 )
  a = -W( 92 ) / JVS(          838  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 839 )
  W( 94 ) = W( 94 ) + a*JVS( 840 )
  W( 95 ) = W( 95 ) + a*JVS( 841 )
  W( 96 ) = W( 96 ) + a*JVS( 842 )
  W( 97 ) = W( 97 ) + a*JVS( 843 )
  W( 98 ) = W( 98 ) + a*JVS( 844 )
  W( 99 ) = W( 99 ) + a*JVS( 845 )
  W( 100 ) = W( 100 ) + a*JVS( 846 )
  a = -W( 93 ) / JVS(          880  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 881 )
  W( 95 ) = W( 95 ) + a*JVS( 882 )
  W( 96 ) = W( 96 ) + a*JVS( 883 )
  W( 97 ) = W( 97 ) + a*JVS( 884 )
  W( 98 ) = W( 98 ) + a*JVS( 885 )
  W( 99 ) = W( 99 ) + a*JVS( 886 )
  W( 100 ) = W( 100 ) + a*JVS( 887 )
  a = -W( 94 ) / JVS(          913  )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 914 )
  W( 96 ) = W( 96 ) + a*JVS( 915 )
  W( 97 ) = W( 97 ) + a*JVS( 916 )
  W( 98 ) = W( 98 ) + a*JVS( 917 )
  W( 99 ) = W( 99 ) + a*JVS( 918 )
  W( 100 ) = W( 100 ) + a*JVS( 919 )
  a = -W( 95 ) / JVS(          931  )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 932 )
  W( 97 ) = W( 97 ) + a*JVS( 933 )
  W( 98 ) = W( 98 ) + a*JVS( 934 )
  W( 99 ) = W( 99 ) + a*JVS( 935 )
  W( 100 ) = W( 100 ) + a*JVS( 936 )
  JVS( 937) = W( 37 )
  JVS( 938) = W( 71 )
  JVS( 939) = W( 74 )
  JVS( 940) = W( 78 )
  JVS( 941) = W( 80 )
  JVS( 942) = W( 81 )
  JVS( 943) = W( 84 )
  JVS( 944) = W( 88 )
  JVS( 945) = W( 89 )
  JVS( 946) = W( 90 )
  JVS( 947) = W( 91 )
  JVS( 948) = W( 92 )
  JVS( 949) = W( 93 )
  JVS( 950) = W( 94 )
  JVS( 951) = W( 95 )
  JVS( 952) = W( 96 )
  JVS( 953) = W( 97 )
  JVS( 954) = W( 98 )
  JVS( 955) = W( 99 )
  JVS( 956) = W( 100 )
  IF ( ABS(  JVS( 1004 )) < TINY(a) ) THEN
         IER = 97                                      
         RETURN
  END IF
   W( 32 ) = JVS( 957 )
   W( 38 ) = JVS( 958 )
   W( 41 ) = JVS( 959 )
   W( 43 ) = JVS( 960 )
   W( 47 ) = JVS( 961 )
   W( 48 ) = JVS( 962 )
   W( 49 ) = JVS( 963 )
   W( 50 ) = JVS( 964 )
   W( 51 ) = JVS( 965 )
   W( 52 ) = JVS( 966 )
   W( 54 ) = JVS( 967 )
   W( 55 ) = JVS( 968 )
   W( 56 ) = JVS( 969 )
   W( 58 ) = JVS( 970 )
   W( 59 ) = JVS( 971 )
   W( 60 ) = JVS( 972 )
   W( 61 ) = JVS( 973 )
   W( 64 ) = JVS( 974 )
   W( 65 ) = JVS( 975 )
   W( 67 ) = JVS( 976 )
   W( 69 ) = JVS( 977 )
   W( 70 ) = JVS( 978 )
   W( 71 ) = JVS( 979 )
   W( 72 ) = JVS( 980 )
   W( 74 ) = JVS( 981 )
   W( 75 ) = JVS( 982 )
   W( 76 ) = JVS( 983 )
   W( 77 ) = JVS( 984 )
   W( 78 ) = JVS( 985 )
   W( 79 ) = JVS( 986 )
   W( 80 ) = JVS( 987 )
   W( 81 ) = JVS( 988 )
   W( 82 ) = JVS( 989 )
   W( 83 ) = JVS( 990 )
   W( 84 ) = JVS( 991 )
   W( 85 ) = JVS( 992 )
   W( 86 ) = JVS( 993 )
   W( 87 ) = JVS( 994 )
   W( 88 ) = JVS( 995 )
   W( 89 ) = JVS( 996 )
   W( 90 ) = JVS( 997 )
   W( 91 ) = JVS( 998 )
   W( 92 ) = JVS( 999 )
   W( 93 ) = JVS( 1000 )
   W( 94 ) = JVS( 1001 )
   W( 95 ) = JVS( 1002 )
   W( 96 ) = JVS( 1003 )
   W( 97 ) = JVS( 1004 )
   W( 98 ) = JVS( 1005 )
   W( 99 ) = JVS( 1006 )
   W( 100 ) = JVS( 1007 )
  a = -W( 32 ) / JVS(          194  )
  W( 32 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 195 )
  a = -W( 38 ) / JVS(          210  )
  W( 38 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 211 )
  W( 98 ) = W( 98 ) + a*JVS( 212 )
  a = -W( 41 ) / JVS(          221  )
  W( 41 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 222 )
  a = -W( 43 ) / JVS(          226  )
  W( 43 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 227 )
  W( 100 ) = W( 100 ) + a*JVS( 228 )
  a = -W( 47 ) / JVS(          237  )
  W( 47 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 238 )
  W( 97 ) = W( 97 ) + a*JVS( 239 )
  W( 98 ) = W( 98 ) + a*JVS( 240 )
  a = -W( 48 ) / JVS(          241  )
  W( 48 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 242 )
  W( 97 ) = W( 97 ) + a*JVS( 243 )
  W( 100 ) = W( 100 ) + a*JVS( 244 )
  a = -W( 49 ) / JVS(          245  )
  W( 49 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 246 )
  W( 92 ) = W( 92 ) + a*JVS( 247 )
  W( 93 ) = W( 93 ) + a*JVS( 248 )
  W( 97 ) = W( 97 ) + a*JVS( 249 )
  a = -W( 50 ) / JVS(          250  )
  W( 50 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 251 )
  a = -W( 51 ) / JVS(          252  )
  W( 51 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 253 )
  W( 97 ) = W( 97 ) + a*JVS( 254 )
  W( 98 ) = W( 98 ) + a*JVS( 255 )
  a = -W( 52 ) / JVS(          256  )
  W( 52 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 257 )
  W( 91 ) = W( 91 ) + a*JVS( 258 )
  W( 94 ) = W( 94 ) + a*JVS( 259 )
  W( 98 ) = W( 98 ) + a*JVS( 260 )
  a = -W( 54 ) / JVS(          263  )
  W( 54 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 264 )
  a = -W( 55 ) / JVS(          267  )
  W( 55 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 268 )
  a = -W( 56 ) / JVS(          271  )
  W( 56 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 272 )
  a = -W( 58 ) / JVS(          280  )
  W( 58 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 281 )
  W( 98 ) = W( 98 ) + a*JVS( 282 )
  a = -W( 59 ) / JVS(          283  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 284 )
  W( 98 ) = W( 98 ) + a*JVS( 285 )
  a = -W( 60 ) / JVS(          286  )
  W( 60 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 287 )
  W( 91 ) = W( 91 ) + a*JVS( 288 )
  W( 97 ) = W( 97 ) + a*JVS( 289 )
  W( 98 ) = W( 98 ) + a*JVS( 290 )
  a = -W( 61 ) / JVS(          292  )
  W( 61 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 293 )
  W( 92 ) = W( 92 ) + a*JVS( 294 )
  W( 93 ) = W( 93 ) + a*JVS( 295 )
  W( 97 ) = W( 97 ) + a*JVS( 296 )
  a = -W( 64 ) / JVS(          315  )
  W( 64 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 316 )
  W( 78 ) = W( 78 ) + a*JVS( 317 )
  W( 81 ) = W( 81 ) + a*JVS( 318 )
  W( 88 ) = W( 88 ) + a*JVS( 319 )
  W( 93 ) = W( 93 ) + a*JVS( 320 )
  W( 98 ) = W( 98 ) + a*JVS( 321 )
  a = -W( 65 ) / JVS(          328  )
  W( 65 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 329 )
  W( 69 ) = W( 69 ) + a*JVS( 330 )
  W( 71 ) = W( 71 ) + a*JVS( 331 )
  W( 72 ) = W( 72 ) + a*JVS( 332 )
  W( 74 ) = W( 74 ) + a*JVS( 333 )
  W( 75 ) = W( 75 ) + a*JVS( 334 )
  W( 76 ) = W( 76 ) + a*JVS( 335 )
  W( 77 ) = W( 77 ) + a*JVS( 336 )
  W( 78 ) = W( 78 ) + a*JVS( 337 )
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 81 ) = W( 81 ) + a*JVS( 340 )
  W( 82 ) = W( 82 ) + a*JVS( 341 )
  W( 84 ) = W( 84 ) + a*JVS( 342 )
  W( 85 ) = W( 85 ) + a*JVS( 343 )
  W( 88 ) = W( 88 ) + a*JVS( 344 )
  W( 93 ) = W( 93 ) + a*JVS( 345 )
  W( 98 ) = W( 98 ) + a*JVS( 346 )
  a = -W( 67 ) / JVS(          368  )
  W( 67 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 369 )
  W( 88 ) = W( 88 ) + a*JVS( 370 )
  W( 93 ) = W( 93 ) + a*JVS( 371 )
  W( 98 ) = W( 98 ) + a*JVS( 372 )
  a = -W( 69 ) / JVS(          388  )
  W( 69 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 389 )
  W( 88 ) = W( 88 ) + a*JVS( 390 )
  W( 93 ) = W( 93 ) + a*JVS( 391 )
  W( 98 ) = W( 98 ) + a*JVS( 392 )
  a = -W( 70 ) / JVS(          395  )
  W( 70 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 396 )
  W( 92 ) = W( 92 ) + a*JVS( 397 )
  W( 93 ) = W( 93 ) + a*JVS( 398 )
  W( 95 ) = W( 95 ) + a*JVS( 399 )
  W( 96 ) = W( 96 ) + a*JVS( 400 )
  W( 97 ) = W( 97 ) + a*JVS( 401 )
  W( 98 ) = W( 98 ) + a*JVS( 402 )
  W( 99 ) = W( 99 ) + a*JVS( 403 )
  W( 100 ) = W( 100 ) + a*JVS( 404 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 74 ) / JVS(          446  )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 447 )
  W( 84 ) = W( 84 ) + a*JVS( 448 )
  W( 88 ) = W( 88 ) + a*JVS( 449 )
  W( 93 ) = W( 93 ) + a*JVS( 450 )
  W( 98 ) = W( 98 ) + a*JVS( 451 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 77 ) / JVS(          473  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 474 )
  W( 84 ) = W( 84 ) + a*JVS( 475 )
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 89 ) = W( 89 ) + a*JVS( 477 )
  W( 92 ) = W( 92 ) + a*JVS( 478 )
  W( 93 ) = W( 93 ) + a*JVS( 479 )
  W( 95 ) = W( 95 ) + a*JVS( 480 )
  W( 96 ) = W( 96 ) + a*JVS( 481 )
  W( 97 ) = W( 97 ) + a*JVS( 482 )
  W( 98 ) = W( 98 ) + a*JVS( 483 )
  W( 99 ) = W( 99 ) + a*JVS( 484 )
  W( 100 ) = W( 100 ) + a*JVS( 485 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 79 ) / JVS(          502  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 503 )
  W( 83 ) = W( 83 ) + a*JVS( 504 )
  W( 84 ) = W( 84 ) + a*JVS( 505 )
  W( 85 ) = W( 85 ) + a*JVS( 506 )
  W( 86 ) = W( 86 ) + a*JVS( 507 )
  W( 87 ) = W( 87 ) + a*JVS( 508 )
  W( 88 ) = W( 88 ) + a*JVS( 509 )
  W( 89 ) = W( 89 ) + a*JVS( 510 )
  W( 93 ) = W( 93 ) + a*JVS( 511 )
  W( 95 ) = W( 95 ) + a*JVS( 512 )
  W( 96 ) = W( 96 ) + a*JVS( 513 )
  W( 98 ) = W( 98 ) + a*JVS( 514 )
  W( 99 ) = W( 99 ) + a*JVS( 515 )
  W( 100 ) = W( 100 ) + a*JVS( 516 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 82 ) / JVS(          549  )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 550 )
  W( 84 ) = W( 84 ) + a*JVS( 551 )
  W( 86 ) = W( 86 ) + a*JVS( 552 )
  W( 87 ) = W( 87 ) + a*JVS( 553 )
  W( 88 ) = W( 88 ) + a*JVS( 554 )
  W( 89 ) = W( 89 ) + a*JVS( 555 )
  W( 90 ) = W( 90 ) + a*JVS( 556 )
  W( 91 ) = W( 91 ) + a*JVS( 557 )
  W( 92 ) = W( 92 ) + a*JVS( 558 )
  W( 93 ) = W( 93 ) + a*JVS( 559 )
  W( 94 ) = W( 94 ) + a*JVS( 560 )
  W( 95 ) = W( 95 ) + a*JVS( 561 )
  W( 96 ) = W( 96 ) + a*JVS( 562 )
  W( 97 ) = W( 97 ) + a*JVS( 563 )
  W( 98 ) = W( 98 ) + a*JVS( 564 )
  W( 99 ) = W( 99 ) + a*JVS( 565 )
  W( 100 ) = W( 100 ) + a*JVS( 566 )
  a = -W( 83 ) / JVS(          573  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 574 )
  W( 88 ) = W( 88 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 98 ) = W( 98 ) + a*JVS( 579 )
  W( 100 ) = W( 100 ) + a*JVS( 580 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  a = -W( 85 ) / JVS(          616  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 617 )
  W( 87 ) = W( 87 ) + a*JVS( 618 )
  W( 88 ) = W( 88 ) + a*JVS( 619 )
  W( 90 ) = W( 90 ) + a*JVS( 620 )
  W( 91 ) = W( 91 ) + a*JVS( 621 )
  W( 92 ) = W( 92 ) + a*JVS( 622 )
  W( 93 ) = W( 93 ) + a*JVS( 623 )
  W( 97 ) = W( 97 ) + a*JVS( 624 )
  W( 98 ) = W( 98 ) + a*JVS( 625 )
  W( 100 ) = W( 100 ) + a*JVS( 626 )
  a = -W( 86 ) / JVS(          638  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 639 )
  W( 88 ) = W( 88 ) + a*JVS( 640 )
  W( 90 ) = W( 90 ) + a*JVS( 641 )
  W( 91 ) = W( 91 ) + a*JVS( 642 )
  W( 92 ) = W( 92 ) + a*JVS( 643 )
  W( 93 ) = W( 93 ) + a*JVS( 644 )
  W( 94 ) = W( 94 ) + a*JVS( 645 )
  W( 98 ) = W( 98 ) + a*JVS( 646 )
  W( 100 ) = W( 100 ) + a*JVS( 647 )
  a = -W( 87 ) / JVS(          660  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 661 )
  W( 89 ) = W( 89 ) + a*JVS( 662 )
  W( 90 ) = W( 90 ) + a*JVS( 663 )
  W( 91 ) = W( 91 ) + a*JVS( 664 )
  W( 92 ) = W( 92 ) + a*JVS( 665 )
  W( 93 ) = W( 93 ) + a*JVS( 666 )
  W( 94 ) = W( 94 ) + a*JVS( 667 )
  W( 95 ) = W( 95 ) + a*JVS( 668 )
  W( 98 ) = W( 98 ) + a*JVS( 669 )
  W( 99 ) = W( 99 ) + a*JVS( 670 )
  W( 100 ) = W( 100 ) + a*JVS( 671 )
  a = -W( 88 ) / JVS(          685  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 686 )
  W( 92 ) = W( 92 ) + a*JVS( 687 )
  W( 93 ) = W( 93 ) + a*JVS( 688 )
  W( 95 ) = W( 95 ) + a*JVS( 689 )
  W( 96 ) = W( 96 ) + a*JVS( 690 )
  W( 97 ) = W( 97 ) + a*JVS( 691 )
  W( 98 ) = W( 98 ) + a*JVS( 692 )
  W( 99 ) = W( 99 ) + a*JVS( 693 )
  W( 100 ) = W( 100 ) + a*JVS( 694 )
  a = -W( 89 ) / JVS(          717  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 718 )
  W( 91 ) = W( 91 ) + a*JVS( 719 )
  W( 92 ) = W( 92 ) + a*JVS( 720 )
  W( 93 ) = W( 93 ) + a*JVS( 721 )
  W( 94 ) = W( 94 ) + a*JVS( 722 )
  W( 95 ) = W( 95 ) + a*JVS( 723 )
  W( 96 ) = W( 96 ) + a*JVS( 724 )
  W( 97 ) = W( 97 ) + a*JVS( 725 )
  W( 98 ) = W( 98 ) + a*JVS( 726 )
  W( 99 ) = W( 99 ) + a*JVS( 727 )
  W( 100 ) = W( 100 ) + a*JVS( 728 )
  a = -W( 90 ) / JVS(          762  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 92 ) = W( 92 ) + a*JVS( 764 )
  W( 93 ) = W( 93 ) + a*JVS( 765 )
  W( 94 ) = W( 94 ) + a*JVS( 766 )
  W( 95 ) = W( 95 ) + a*JVS( 767 )
  W( 96 ) = W( 96 ) + a*JVS( 768 )
  W( 97 ) = W( 97 ) + a*JVS( 769 )
  W( 98 ) = W( 98 ) + a*JVS( 770 )
  W( 99 ) = W( 99 ) + a*JVS( 771 )
  W( 100 ) = W( 100 ) + a*JVS( 772 )
  a = -W( 91 ) / JVS(          795  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 796 )
  W( 93 ) = W( 93 ) + a*JVS( 797 )
  W( 94 ) = W( 94 ) + a*JVS( 798 )
  W( 95 ) = W( 95 ) + a*JVS( 799 )
  W( 96 ) = W( 96 ) + a*JVS( 800 )
  W( 97 ) = W( 97 ) + a*JVS( 801 )
  W( 98 ) = W( 98 ) + a*JVS( 802 )
  W( 99 ) = W( 99 ) + a*JVS( 803 )
  W( 100 ) = W( 100 ) + a*JVS( 804 )
  a = -W( 92 ) / JVS(          838  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 839 )
  W( 94 ) = W( 94 ) + a*JVS( 840 )
  W( 95 ) = W( 95 ) + a*JVS( 841 )
  W( 96 ) = W( 96 ) + a*JVS( 842 )
  W( 97 ) = W( 97 ) + a*JVS( 843 )
  W( 98 ) = W( 98 ) + a*JVS( 844 )
  W( 99 ) = W( 99 ) + a*JVS( 845 )
  W( 100 ) = W( 100 ) + a*JVS( 846 )
  a = -W( 93 ) / JVS(          880  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 881 )
  W( 95 ) = W( 95 ) + a*JVS( 882 )
  W( 96 ) = W( 96 ) + a*JVS( 883 )
  W( 97 ) = W( 97 ) + a*JVS( 884 )
  W( 98 ) = W( 98 ) + a*JVS( 885 )
  W( 99 ) = W( 99 ) + a*JVS( 886 )
  W( 100 ) = W( 100 ) + a*JVS( 887 )
  a = -W( 94 ) / JVS(          913  )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 914 )
  W( 96 ) = W( 96 ) + a*JVS( 915 )
  W( 97 ) = W( 97 ) + a*JVS( 916 )
  W( 98 ) = W( 98 ) + a*JVS( 917 )
  W( 99 ) = W( 99 ) + a*JVS( 918 )
  W( 100 ) = W( 100 ) + a*JVS( 919 )
  a = -W( 95 ) / JVS(          931  )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 932 )
  W( 97 ) = W( 97 ) + a*JVS( 933 )
  W( 98 ) = W( 98 ) + a*JVS( 934 )
  W( 99 ) = W( 99 ) + a*JVS( 935 )
  W( 100 ) = W( 100 ) + a*JVS( 936 )
  a = -W( 96 ) / JVS(          952  )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 953 )
  W( 98 ) = W( 98 ) + a*JVS( 954 )
  W( 99 ) = W( 99 ) + a*JVS( 955 )
  W( 100 ) = W( 100 ) + a*JVS( 956 )
  JVS( 957) = W( 32 )
  JVS( 958) = W( 38 )
  JVS( 959) = W( 41 )
  JVS( 960) = W( 43 )
  JVS( 961) = W( 47 )
  JVS( 962) = W( 48 )
  JVS( 963) = W( 49 )
  JVS( 964) = W( 50 )
  JVS( 965) = W( 51 )
  JVS( 966) = W( 52 )
  JVS( 967) = W( 54 )
  JVS( 968) = W( 55 )
  JVS( 969) = W( 56 )
  JVS( 970) = W( 58 )
  JVS( 971) = W( 59 )
  JVS( 972) = W( 60 )
  JVS( 973) = W( 61 )
  JVS( 974) = W( 64 )
  JVS( 975) = W( 65 )
  JVS( 976) = W( 67 )
  JVS( 977) = W( 69 )
  JVS( 978) = W( 70 )
  JVS( 979) = W( 71 )
  JVS( 980) = W( 72 )
  JVS( 981) = W( 74 )
  JVS( 982) = W( 75 )
  JVS( 983) = W( 76 )
  JVS( 984) = W( 77 )
  JVS( 985) = W( 78 )
  JVS( 986) = W( 79 )
  JVS( 987) = W( 80 )
  JVS( 988) = W( 81 )
  JVS( 989) = W( 82 )
  JVS( 990) = W( 83 )
  JVS( 991) = W( 84 )
  JVS( 992) = W( 85 )
  JVS( 993) = W( 86 )
  JVS( 994) = W( 87 )
  JVS( 995) = W( 88 )
  JVS( 996) = W( 89 )
  JVS( 997) = W( 90 )
  JVS( 998) = W( 91 )
  JVS( 999) = W( 92 )
  JVS( 1000) = W( 93 )
  JVS( 1001) = W( 94 )
  JVS( 1002) = W( 95 )
  JVS( 1003) = W( 96 )
  JVS( 1004) = W( 97 )
  JVS( 1005) = W( 98 )
  JVS( 1006) = W( 99 )
  JVS( 1007) = W( 100 )
  IF ( ABS(  JVS( 1071 )) < TINY(a) ) THEN
         IER = 98                                      
         RETURN
  END IF
   W( 3 ) = JVS( 1008 )
   W( 5 ) = JVS( 1009 )
   W( 26 ) = JVS( 1010 )
   W( 27 ) = JVS( 1011 )
   W( 28 ) = JVS( 1012 )
   W( 29 ) = JVS( 1013 )
   W( 30 ) = JVS( 1014 )
   W( 31 ) = JVS( 1015 )
   W( 32 ) = JVS( 1016 )
   W( 33 ) = JVS( 1017 )
   W( 38 ) = JVS( 1018 )
   W( 39 ) = JVS( 1019 )
   W( 41 ) = JVS( 1020 )
   W( 43 ) = JVS( 1021 )
   W( 44 ) = JVS( 1022 )
   W( 46 ) = JVS( 1023 )
   W( 47 ) = JVS( 1024 )
   W( 50 ) = JVS( 1025 )
   W( 51 ) = JVS( 1026 )
   W( 52 ) = JVS( 1027 )
   W( 53 ) = JVS( 1028 )
   W( 54 ) = JVS( 1029 )
   W( 55 ) = JVS( 1030 )
   W( 56 ) = JVS( 1031 )
   W( 57 ) = JVS( 1032 )
   W( 58 ) = JVS( 1033 )
   W( 59 ) = JVS( 1034 )
   W( 60 ) = JVS( 1035 )
   W( 62 ) = JVS( 1036 )
   W( 63 ) = JVS( 1037 )
   W( 64 ) = JVS( 1038 )
   W( 65 ) = JVS( 1039 )
   W( 66 ) = JVS( 1040 )
   W( 67 ) = JVS( 1041 )
   W( 68 ) = JVS( 1042 )
   W( 69 ) = JVS( 1043 )
   W( 70 ) = JVS( 1044 )
   W( 71 ) = JVS( 1045 )
   W( 72 ) = JVS( 1046 )
   W( 74 ) = JVS( 1047 )
   W( 75 ) = JVS( 1048 )
   W( 76 ) = JVS( 1049 )
   W( 77 ) = JVS( 1050 )
   W( 78 ) = JVS( 1051 )
   W( 79 ) = JVS( 1052 )
   W( 80 ) = JVS( 1053 )
   W( 81 ) = JVS( 1054 )
   W( 82 ) = JVS( 1055 )
   W( 83 ) = JVS( 1056 )
   W( 84 ) = JVS( 1057 )
   W( 85 ) = JVS( 1058 )
   W( 86 ) = JVS( 1059 )
   W( 87 ) = JVS( 1060 )
   W( 88 ) = JVS( 1061 )
   W( 89 ) = JVS( 1062 )
   W( 90 ) = JVS( 1063 )
   W( 91 ) = JVS( 1064 )
   W( 92 ) = JVS( 1065 )
   W( 93 ) = JVS( 1066 )
   W( 94 ) = JVS( 1067 )
   W( 95 ) = JVS( 1068 )
   W( 96 ) = JVS( 1069 )
   W( 97 ) = JVS( 1070 )
   W( 98 ) = JVS( 1071 )
   W( 99 ) = JVS( 1072 )
   W( 100 ) = JVS( 1073 )
  a = -W( 3 ) / JVS(            8  )
  W( 3 ) = -a
  a = -W( 5 ) / JVS(           10  )
  W( 5 ) = -a
  a = -W( 26 ) / JVS(          180  )
  W( 26 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 181 )
  a = -W( 27 ) / JVS(          182  )
  W( 27 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 183 )
  a = -W( 28 ) / JVS(          185  )
  W( 28 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 186 )
  a = -W( 29 ) / JVS(          187  )
  W( 29 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 188 )
  a = -W( 30 ) / JVS(          190  )
  W( 30 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 191 )
  a = -W( 31 ) / JVS(          192  )
  W( 31 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 193 )
  a = -W( 32 ) / JVS(          194  )
  W( 32 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 195 )
  a = -W( 33 ) / JVS(          196  )
  W( 33 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 197 )
  a = -W( 38 ) / JVS(          210  )
  W( 38 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 211 )
  W( 98 ) = W( 98 ) + a*JVS( 212 )
  a = -W( 39 ) / JVS(          213  )
  W( 39 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 214 )
  a = -W( 41 ) / JVS(          221  )
  W( 41 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 222 )
  a = -W( 43 ) / JVS(          226  )
  W( 43 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 227 )
  W( 100 ) = W( 100 ) + a*JVS( 228 )
  a = -W( 44 ) / JVS(          229  )
  W( 44 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 230 )
  a = -W( 46 ) / JVS(          235  )
  W( 46 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 236 )
  a = -W( 47 ) / JVS(          237  )
  W( 47 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 238 )
  W( 97 ) = W( 97 ) + a*JVS( 239 )
  W( 98 ) = W( 98 ) + a*JVS( 240 )
  a = -W( 50 ) / JVS(          250  )
  W( 50 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 251 )
  a = -W( 51 ) / JVS(          252  )
  W( 51 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 253 )
  W( 97 ) = W( 97 ) + a*JVS( 254 )
  W( 98 ) = W( 98 ) + a*JVS( 255 )
  a = -W( 52 ) / JVS(          256  )
  W( 52 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 257 )
  W( 91 ) = W( 91 ) + a*JVS( 258 )
  W( 94 ) = W( 94 ) + a*JVS( 259 )
  W( 98 ) = W( 98 ) + a*JVS( 260 )
  a = -W( 53 ) / JVS(          261  )
  W( 53 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 262 )
  a = -W( 54 ) / JVS(          263  )
  W( 54 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 264 )
  a = -W( 55 ) / JVS(          267  )
  W( 55 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 268 )
  a = -W( 56 ) / JVS(          271  )
  W( 56 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 272 )
  a = -W( 57 ) / JVS(          275  )
  W( 57 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 276 )
  W( 98 ) = W( 98 ) + a*JVS( 277 )
  a = -W( 58 ) / JVS(          280  )
  W( 58 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 281 )
  W( 98 ) = W( 98 ) + a*JVS( 282 )
  a = -W( 59 ) / JVS(          283  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 284 )
  W( 98 ) = W( 98 ) + a*JVS( 285 )
  a = -W( 60 ) / JVS(          286  )
  W( 60 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 287 )
  W( 91 ) = W( 91 ) + a*JVS( 288 )
  W( 97 ) = W( 97 ) + a*JVS( 289 )
  W( 98 ) = W( 98 ) + a*JVS( 290 )
  a = -W( 62 ) / JVS(          299  )
  W( 62 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 300 )
  W( 88 ) = W( 88 ) + a*JVS( 301 )
  W( 93 ) = W( 93 ) + a*JVS( 302 )
  W( 98 ) = W( 98 ) + a*JVS( 303 )
  a = -W( 63 ) / JVS(          305  )
  W( 63 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 306 )
  W( 93 ) = W( 93 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  W( 98 ) = W( 98 ) + a*JVS( 309 )
  a = -W( 64 ) / JVS(          315  )
  W( 64 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 316 )
  W( 78 ) = W( 78 ) + a*JVS( 317 )
  W( 81 ) = W( 81 ) + a*JVS( 318 )
  W( 88 ) = W( 88 ) + a*JVS( 319 )
  W( 93 ) = W( 93 ) + a*JVS( 320 )
  W( 98 ) = W( 98 ) + a*JVS( 321 )
  a = -W( 65 ) / JVS(          328  )
  W( 65 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 329 )
  W( 69 ) = W( 69 ) + a*JVS( 330 )
  W( 71 ) = W( 71 ) + a*JVS( 331 )
  W( 72 ) = W( 72 ) + a*JVS( 332 )
  W( 74 ) = W( 74 ) + a*JVS( 333 )
  W( 75 ) = W( 75 ) + a*JVS( 334 )
  W( 76 ) = W( 76 ) + a*JVS( 335 )
  W( 77 ) = W( 77 ) + a*JVS( 336 )
  W( 78 ) = W( 78 ) + a*JVS( 337 )
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 81 ) = W( 81 ) + a*JVS( 340 )
  W( 82 ) = W( 82 ) + a*JVS( 341 )
  W( 84 ) = W( 84 ) + a*JVS( 342 )
  W( 85 ) = W( 85 ) + a*JVS( 343 )
  W( 88 ) = W( 88 ) + a*JVS( 344 )
  W( 93 ) = W( 93 ) + a*JVS( 345 )
  W( 98 ) = W( 98 ) + a*JVS( 346 )
  a = -W( 66 ) / JVS(          353  )
  W( 66 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 354 )
  W( 74 ) = W( 74 ) + a*JVS( 355 )
  W( 77 ) = W( 77 ) + a*JVS( 356 )
  W( 78 ) = W( 78 ) + a*JVS( 357 )
  W( 79 ) = W( 79 ) + a*JVS( 358 )
  W( 80 ) = W( 80 ) + a*JVS( 359 )
  W( 81 ) = W( 81 ) + a*JVS( 360 )
  W( 82 ) = W( 82 ) + a*JVS( 361 )
  W( 85 ) = W( 85 ) + a*JVS( 362 )
  W( 88 ) = W( 88 ) + a*JVS( 363 )
  W( 92 ) = W( 92 ) + a*JVS( 364 )
  W( 93 ) = W( 93 ) + a*JVS( 365 )
  W( 97 ) = W( 97 ) + a*JVS( 366 )
  W( 98 ) = W( 98 ) + a*JVS( 367 )
  a = -W( 67 ) / JVS(          368  )
  W( 67 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 369 )
  W( 88 ) = W( 88 ) + a*JVS( 370 )
  W( 93 ) = W( 93 ) + a*JVS( 371 )
  W( 98 ) = W( 98 ) + a*JVS( 372 )
  a = -W( 68 ) / JVS(          378  )
  W( 68 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 80 ) = W( 80 ) + a*JVS( 382 )
  W( 83 ) = W( 83 ) + a*JVS( 383 )
  W( 88 ) = W( 88 ) + a*JVS( 384 )
  W( 92 ) = W( 92 ) + a*JVS( 385 )
  W( 93 ) = W( 93 ) + a*JVS( 386 )
  W( 98 ) = W( 98 ) + a*JVS( 387 )
  a = -W( 69 ) / JVS(          388  )
  W( 69 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 389 )
  W( 88 ) = W( 88 ) + a*JVS( 390 )
  W( 93 ) = W( 93 ) + a*JVS( 391 )
  W( 98 ) = W( 98 ) + a*JVS( 392 )
  a = -W( 70 ) / JVS(          395  )
  W( 70 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 396 )
  W( 92 ) = W( 92 ) + a*JVS( 397 )
  W( 93 ) = W( 93 ) + a*JVS( 398 )
  W( 95 ) = W( 95 ) + a*JVS( 399 )
  W( 96 ) = W( 96 ) + a*JVS( 400 )
  W( 97 ) = W( 97 ) + a*JVS( 401 )
  W( 98 ) = W( 98 ) + a*JVS( 402 )
  W( 99 ) = W( 99 ) + a*JVS( 403 )
  W( 100 ) = W( 100 ) + a*JVS( 404 )
  a = -W( 71 ) / JVS(          405  )
  W( 71 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 406 )
  W( 88 ) = W( 88 ) + a*JVS( 407 )
  W( 93 ) = W( 93 ) + a*JVS( 408 )
  W( 98 ) = W( 98 ) + a*JVS( 409 )
  a = -W( 72 ) / JVS(          410  )
  W( 72 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 88 ) = W( 88 ) + a*JVS( 412 )
  W( 93 ) = W( 93 ) + a*JVS( 413 )
  W( 98 ) = W( 98 ) + a*JVS( 414 )
  a = -W( 74 ) / JVS(          446  )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 447 )
  W( 84 ) = W( 84 ) + a*JVS( 448 )
  W( 88 ) = W( 88 ) + a*JVS( 449 )
  W( 93 ) = W( 93 ) + a*JVS( 450 )
  W( 98 ) = W( 98 ) + a*JVS( 451 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 77 ) / JVS(          473  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 474 )
  W( 84 ) = W( 84 ) + a*JVS( 475 )
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 89 ) = W( 89 ) + a*JVS( 477 )
  W( 92 ) = W( 92 ) + a*JVS( 478 )
  W( 93 ) = W( 93 ) + a*JVS( 479 )
  W( 95 ) = W( 95 ) + a*JVS( 480 )
  W( 96 ) = W( 96 ) + a*JVS( 481 )
  W( 97 ) = W( 97 ) + a*JVS( 482 )
  W( 98 ) = W( 98 ) + a*JVS( 483 )
  W( 99 ) = W( 99 ) + a*JVS( 484 )
  W( 100 ) = W( 100 ) + a*JVS( 485 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 79 ) / JVS(          502  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 503 )
  W( 83 ) = W( 83 ) + a*JVS( 504 )
  W( 84 ) = W( 84 ) + a*JVS( 505 )
  W( 85 ) = W( 85 ) + a*JVS( 506 )
  W( 86 ) = W( 86 ) + a*JVS( 507 )
  W( 87 ) = W( 87 ) + a*JVS( 508 )
  W( 88 ) = W( 88 ) + a*JVS( 509 )
  W( 89 ) = W( 89 ) + a*JVS( 510 )
  W( 93 ) = W( 93 ) + a*JVS( 511 )
  W( 95 ) = W( 95 ) + a*JVS( 512 )
  W( 96 ) = W( 96 ) + a*JVS( 513 )
  W( 98 ) = W( 98 ) + a*JVS( 514 )
  W( 99 ) = W( 99 ) + a*JVS( 515 )
  W( 100 ) = W( 100 ) + a*JVS( 516 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 82 ) / JVS(          549  )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 550 )
  W( 84 ) = W( 84 ) + a*JVS( 551 )
  W( 86 ) = W( 86 ) + a*JVS( 552 )
  W( 87 ) = W( 87 ) + a*JVS( 553 )
  W( 88 ) = W( 88 ) + a*JVS( 554 )
  W( 89 ) = W( 89 ) + a*JVS( 555 )
  W( 90 ) = W( 90 ) + a*JVS( 556 )
  W( 91 ) = W( 91 ) + a*JVS( 557 )
  W( 92 ) = W( 92 ) + a*JVS( 558 )
  W( 93 ) = W( 93 ) + a*JVS( 559 )
  W( 94 ) = W( 94 ) + a*JVS( 560 )
  W( 95 ) = W( 95 ) + a*JVS( 561 )
  W( 96 ) = W( 96 ) + a*JVS( 562 )
  W( 97 ) = W( 97 ) + a*JVS( 563 )
  W( 98 ) = W( 98 ) + a*JVS( 564 )
  W( 99 ) = W( 99 ) + a*JVS( 565 )
  W( 100 ) = W( 100 ) + a*JVS( 566 )
  a = -W( 83 ) / JVS(          573  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 574 )
  W( 88 ) = W( 88 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 98 ) = W( 98 ) + a*JVS( 579 )
  W( 100 ) = W( 100 ) + a*JVS( 580 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  a = -W( 85 ) / JVS(          616  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 617 )
  W( 87 ) = W( 87 ) + a*JVS( 618 )
  W( 88 ) = W( 88 ) + a*JVS( 619 )
  W( 90 ) = W( 90 ) + a*JVS( 620 )
  W( 91 ) = W( 91 ) + a*JVS( 621 )
  W( 92 ) = W( 92 ) + a*JVS( 622 )
  W( 93 ) = W( 93 ) + a*JVS( 623 )
  W( 97 ) = W( 97 ) + a*JVS( 624 )
  W( 98 ) = W( 98 ) + a*JVS( 625 )
  W( 100 ) = W( 100 ) + a*JVS( 626 )
  a = -W( 86 ) / JVS(          638  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 639 )
  W( 88 ) = W( 88 ) + a*JVS( 640 )
  W( 90 ) = W( 90 ) + a*JVS( 641 )
  W( 91 ) = W( 91 ) + a*JVS( 642 )
  W( 92 ) = W( 92 ) + a*JVS( 643 )
  W( 93 ) = W( 93 ) + a*JVS( 644 )
  W( 94 ) = W( 94 ) + a*JVS( 645 )
  W( 98 ) = W( 98 ) + a*JVS( 646 )
  W( 100 ) = W( 100 ) + a*JVS( 647 )
  a = -W( 87 ) / JVS(          660  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 661 )
  W( 89 ) = W( 89 ) + a*JVS( 662 )
  W( 90 ) = W( 90 ) + a*JVS( 663 )
  W( 91 ) = W( 91 ) + a*JVS( 664 )
  W( 92 ) = W( 92 ) + a*JVS( 665 )
  W( 93 ) = W( 93 ) + a*JVS( 666 )
  W( 94 ) = W( 94 ) + a*JVS( 667 )
  W( 95 ) = W( 95 ) + a*JVS( 668 )
  W( 98 ) = W( 98 ) + a*JVS( 669 )
  W( 99 ) = W( 99 ) + a*JVS( 670 )
  W( 100 ) = W( 100 ) + a*JVS( 671 )
  a = -W( 88 ) / JVS(          685  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 686 )
  W( 92 ) = W( 92 ) + a*JVS( 687 )
  W( 93 ) = W( 93 ) + a*JVS( 688 )
  W( 95 ) = W( 95 ) + a*JVS( 689 )
  W( 96 ) = W( 96 ) + a*JVS( 690 )
  W( 97 ) = W( 97 ) + a*JVS( 691 )
  W( 98 ) = W( 98 ) + a*JVS( 692 )
  W( 99 ) = W( 99 ) + a*JVS( 693 )
  W( 100 ) = W( 100 ) + a*JVS( 694 )
  a = -W( 89 ) / JVS(          717  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 718 )
  W( 91 ) = W( 91 ) + a*JVS( 719 )
  W( 92 ) = W( 92 ) + a*JVS( 720 )
  W( 93 ) = W( 93 ) + a*JVS( 721 )
  W( 94 ) = W( 94 ) + a*JVS( 722 )
  W( 95 ) = W( 95 ) + a*JVS( 723 )
  W( 96 ) = W( 96 ) + a*JVS( 724 )
  W( 97 ) = W( 97 ) + a*JVS( 725 )
  W( 98 ) = W( 98 ) + a*JVS( 726 )
  W( 99 ) = W( 99 ) + a*JVS( 727 )
  W( 100 ) = W( 100 ) + a*JVS( 728 )
  a = -W( 90 ) / JVS(          762  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 92 ) = W( 92 ) + a*JVS( 764 )
  W( 93 ) = W( 93 ) + a*JVS( 765 )
  W( 94 ) = W( 94 ) + a*JVS( 766 )
  W( 95 ) = W( 95 ) + a*JVS( 767 )
  W( 96 ) = W( 96 ) + a*JVS( 768 )
  W( 97 ) = W( 97 ) + a*JVS( 769 )
  W( 98 ) = W( 98 ) + a*JVS( 770 )
  W( 99 ) = W( 99 ) + a*JVS( 771 )
  W( 100 ) = W( 100 ) + a*JVS( 772 )
  a = -W( 91 ) / JVS(          795  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 796 )
  W( 93 ) = W( 93 ) + a*JVS( 797 )
  W( 94 ) = W( 94 ) + a*JVS( 798 )
  W( 95 ) = W( 95 ) + a*JVS( 799 )
  W( 96 ) = W( 96 ) + a*JVS( 800 )
  W( 97 ) = W( 97 ) + a*JVS( 801 )
  W( 98 ) = W( 98 ) + a*JVS( 802 )
  W( 99 ) = W( 99 ) + a*JVS( 803 )
  W( 100 ) = W( 100 ) + a*JVS( 804 )
  a = -W( 92 ) / JVS(          838  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 839 )
  W( 94 ) = W( 94 ) + a*JVS( 840 )
  W( 95 ) = W( 95 ) + a*JVS( 841 )
  W( 96 ) = W( 96 ) + a*JVS( 842 )
  W( 97 ) = W( 97 ) + a*JVS( 843 )
  W( 98 ) = W( 98 ) + a*JVS( 844 )
  W( 99 ) = W( 99 ) + a*JVS( 845 )
  W( 100 ) = W( 100 ) + a*JVS( 846 )
  a = -W( 93 ) / JVS(          880  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 881 )
  W( 95 ) = W( 95 ) + a*JVS( 882 )
  W( 96 ) = W( 96 ) + a*JVS( 883 )
  W( 97 ) = W( 97 ) + a*JVS( 884 )
  W( 98 ) = W( 98 ) + a*JVS( 885 )
  W( 99 ) = W( 99 ) + a*JVS( 886 )
  W( 100 ) = W( 100 ) + a*JVS( 887 )
  a = -W( 94 ) / JVS(          913  )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 914 )
  W( 96 ) = W( 96 ) + a*JVS( 915 )
  W( 97 ) = W( 97 ) + a*JVS( 916 )
  W( 98 ) = W( 98 ) + a*JVS( 917 )
  W( 99 ) = W( 99 ) + a*JVS( 918 )
  W( 100 ) = W( 100 ) + a*JVS( 919 )
  a = -W( 95 ) / JVS(          931  )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 932 )
  W( 97 ) = W( 97 ) + a*JVS( 933 )
  W( 98 ) = W( 98 ) + a*JVS( 934 )
  W( 99 ) = W( 99 ) + a*JVS( 935 )
  W( 100 ) = W( 100 ) + a*JVS( 936 )
  a = -W( 96 ) / JVS(          952  )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 953 )
  W( 98 ) = W( 98 ) + a*JVS( 954 )
  W( 99 ) = W( 99 ) + a*JVS( 955 )
  W( 100 ) = W( 100 ) + a*JVS( 956 )
  a = -W( 97 ) / JVS(         1004  )
  W( 97 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 1005 )
  W( 99 ) = W( 99 ) + a*JVS( 1006 )
  W( 100 ) = W( 100 ) + a*JVS( 1007 )
  JVS( 1008) = W( 3 )
  JVS( 1009) = W( 5 )
  JVS( 1010) = W( 26 )
  JVS( 1011) = W( 27 )
  JVS( 1012) = W( 28 )
  JVS( 1013) = W( 29 )
  JVS( 1014) = W( 30 )
  JVS( 1015) = W( 31 )
  JVS( 1016) = W( 32 )
  JVS( 1017) = W( 33 )
  JVS( 1018) = W( 38 )
  JVS( 1019) = W( 39 )
  JVS( 1020) = W( 41 )
  JVS( 1021) = W( 43 )
  JVS( 1022) = W( 44 )
  JVS( 1023) = W( 46 )
  JVS( 1024) = W( 47 )
  JVS( 1025) = W( 50 )
  JVS( 1026) = W( 51 )
  JVS( 1027) = W( 52 )
  JVS( 1028) = W( 53 )
  JVS( 1029) = W( 54 )
  JVS( 1030) = W( 55 )
  JVS( 1031) = W( 56 )
  JVS( 1032) = W( 57 )
  JVS( 1033) = W( 58 )
  JVS( 1034) = W( 59 )
  JVS( 1035) = W( 60 )
  JVS( 1036) = W( 62 )
  JVS( 1037) = W( 63 )
  JVS( 1038) = W( 64 )
  JVS( 1039) = W( 65 )
  JVS( 1040) = W( 66 )
  JVS( 1041) = W( 67 )
  JVS( 1042) = W( 68 )
  JVS( 1043) = W( 69 )
  JVS( 1044) = W( 70 )
  JVS( 1045) = W( 71 )
  JVS( 1046) = W( 72 )
  JVS( 1047) = W( 74 )
  JVS( 1048) = W( 75 )
  JVS( 1049) = W( 76 )
  JVS( 1050) = W( 77 )
  JVS( 1051) = W( 78 )
  JVS( 1052) = W( 79 )
  JVS( 1053) = W( 80 )
  JVS( 1054) = W( 81 )
  JVS( 1055) = W( 82 )
  JVS( 1056) = W( 83 )
  JVS( 1057) = W( 84 )
  JVS( 1058) = W( 85 )
  JVS( 1059) = W( 86 )
  JVS( 1060) = W( 87 )
  JVS( 1061) = W( 88 )
  JVS( 1062) = W( 89 )
  JVS( 1063) = W( 90 )
  JVS( 1064) = W( 91 )
  JVS( 1065) = W( 92 )
  JVS( 1066) = W( 93 )
  JVS( 1067) = W( 94 )
  JVS( 1068) = W( 95 )
  JVS( 1069) = W( 96 )
  JVS( 1070) = W( 97 )
  JVS( 1071) = W( 98 )
  JVS( 1072) = W( 99 )
  JVS( 1073) = W( 100 )
  IF ( ABS(  JVS( 1097 )) < TINY(a) ) THEN
         IER = 99                                      
         RETURN
  END IF
   W( 35 ) = JVS( 1074 )
   W( 74 ) = JVS( 1075 )
   W( 75 ) = JVS( 1076 )
   W( 76 ) = JVS( 1077 )
   W( 77 ) = JVS( 1078 )
   W( 78 ) = JVS( 1079 )
   W( 80 ) = JVS( 1080 )
   W( 81 ) = JVS( 1081 )
   W( 84 ) = JVS( 1082 )
   W( 85 ) = JVS( 1083 )
   W( 86 ) = JVS( 1084 )
   W( 87 ) = JVS( 1085 )
   W( 88 ) = JVS( 1086 )
   W( 89 ) = JVS( 1087 )
   W( 90 ) = JVS( 1088 )
   W( 91 ) = JVS( 1089 )
   W( 92 ) = JVS( 1090 )
   W( 93 ) = JVS( 1091 )
   W( 94 ) = JVS( 1092 )
   W( 95 ) = JVS( 1093 )
   W( 96 ) = JVS( 1094 )
   W( 97 ) = JVS( 1095 )
   W( 98 ) = JVS( 1096 )
   W( 99 ) = JVS( 1097 )
   W( 100 ) = JVS( 1098 )
  a = -W( 35 ) / JVS(          201  )
  W( 35 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 202 )
  W( 99 ) = W( 99 ) + a*JVS( 203 )
  a = -W( 74 ) / JVS(          446  )
  W( 74 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 447 )
  W( 84 ) = W( 84 ) + a*JVS( 448 )
  W( 88 ) = W( 88 ) + a*JVS( 449 )
  W( 93 ) = W( 93 ) + a*JVS( 450 )
  W( 98 ) = W( 98 ) + a*JVS( 451 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 77 ) / JVS(          473  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 474 )
  W( 84 ) = W( 84 ) + a*JVS( 475 )
  W( 88 ) = W( 88 ) + a*JVS( 476 )
  W( 89 ) = W( 89 ) + a*JVS( 477 )
  W( 92 ) = W( 92 ) + a*JVS( 478 )
  W( 93 ) = W( 93 ) + a*JVS( 479 )
  W( 95 ) = W( 95 ) + a*JVS( 480 )
  W( 96 ) = W( 96 ) + a*JVS( 481 )
  W( 97 ) = W( 97 ) + a*JVS( 482 )
  W( 98 ) = W( 98 ) + a*JVS( 483 )
  W( 99 ) = W( 99 ) + a*JVS( 484 )
  W( 100 ) = W( 100 ) + a*JVS( 485 )
  a = -W( 78 ) / JVS(          487  )
  W( 78 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 488 )
  W( 84 ) = W( 84 ) + a*JVS( 489 )
  W( 88 ) = W( 88 ) + a*JVS( 490 )
  W( 93 ) = W( 93 ) + a*JVS( 491 )
  W( 98 ) = W( 98 ) + a*JVS( 492 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  a = -W( 85 ) / JVS(          616  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 617 )
  W( 87 ) = W( 87 ) + a*JVS( 618 )
  W( 88 ) = W( 88 ) + a*JVS( 619 )
  W( 90 ) = W( 90 ) + a*JVS( 620 )
  W( 91 ) = W( 91 ) + a*JVS( 621 )
  W( 92 ) = W( 92 ) + a*JVS( 622 )
  W( 93 ) = W( 93 ) + a*JVS( 623 )
  W( 97 ) = W( 97 ) + a*JVS( 624 )
  W( 98 ) = W( 98 ) + a*JVS( 625 )
  W( 100 ) = W( 100 ) + a*JVS( 626 )
  a = -W( 86 ) / JVS(          638  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 639 )
  W( 88 ) = W( 88 ) + a*JVS( 640 )
  W( 90 ) = W( 90 ) + a*JVS( 641 )
  W( 91 ) = W( 91 ) + a*JVS( 642 )
  W( 92 ) = W( 92 ) + a*JVS( 643 )
  W( 93 ) = W( 93 ) + a*JVS( 644 )
  W( 94 ) = W( 94 ) + a*JVS( 645 )
  W( 98 ) = W( 98 ) + a*JVS( 646 )
  W( 100 ) = W( 100 ) + a*JVS( 647 )
  a = -W( 87 ) / JVS(          660  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 661 )
  W( 89 ) = W( 89 ) + a*JVS( 662 )
  W( 90 ) = W( 90 ) + a*JVS( 663 )
  W( 91 ) = W( 91 ) + a*JVS( 664 )
  W( 92 ) = W( 92 ) + a*JVS( 665 )
  W( 93 ) = W( 93 ) + a*JVS( 666 )
  W( 94 ) = W( 94 ) + a*JVS( 667 )
  W( 95 ) = W( 95 ) + a*JVS( 668 )
  W( 98 ) = W( 98 ) + a*JVS( 669 )
  W( 99 ) = W( 99 ) + a*JVS( 670 )
  W( 100 ) = W( 100 ) + a*JVS( 671 )
  a = -W( 88 ) / JVS(          685  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 686 )
  W( 92 ) = W( 92 ) + a*JVS( 687 )
  W( 93 ) = W( 93 ) + a*JVS( 688 )
  W( 95 ) = W( 95 ) + a*JVS( 689 )
  W( 96 ) = W( 96 ) + a*JVS( 690 )
  W( 97 ) = W( 97 ) + a*JVS( 691 )
  W( 98 ) = W( 98 ) + a*JVS( 692 )
  W( 99 ) = W( 99 ) + a*JVS( 693 )
  W( 100 ) = W( 100 ) + a*JVS( 694 )
  a = -W( 89 ) / JVS(          717  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 718 )
  W( 91 ) = W( 91 ) + a*JVS( 719 )
  W( 92 ) = W( 92 ) + a*JVS( 720 )
  W( 93 ) = W( 93 ) + a*JVS( 721 )
  W( 94 ) = W( 94 ) + a*JVS( 722 )
  W( 95 ) = W( 95 ) + a*JVS( 723 )
  W( 96 ) = W( 96 ) + a*JVS( 724 )
  W( 97 ) = W( 97 ) + a*JVS( 725 )
  W( 98 ) = W( 98 ) + a*JVS( 726 )
  W( 99 ) = W( 99 ) + a*JVS( 727 )
  W( 100 ) = W( 100 ) + a*JVS( 728 )
  a = -W( 90 ) / JVS(          762  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 92 ) = W( 92 ) + a*JVS( 764 )
  W( 93 ) = W( 93 ) + a*JVS( 765 )
  W( 94 ) = W( 94 ) + a*JVS( 766 )
  W( 95 ) = W( 95 ) + a*JVS( 767 )
  W( 96 ) = W( 96 ) + a*JVS( 768 )
  W( 97 ) = W( 97 ) + a*JVS( 769 )
  W( 98 ) = W( 98 ) + a*JVS( 770 )
  W( 99 ) = W( 99 ) + a*JVS( 771 )
  W( 100 ) = W( 100 ) + a*JVS( 772 )
  a = -W( 91 ) / JVS(          795  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 796 )
  W( 93 ) = W( 93 ) + a*JVS( 797 )
  W( 94 ) = W( 94 ) + a*JVS( 798 )
  W( 95 ) = W( 95 ) + a*JVS( 799 )
  W( 96 ) = W( 96 ) + a*JVS( 800 )
  W( 97 ) = W( 97 ) + a*JVS( 801 )
  W( 98 ) = W( 98 ) + a*JVS( 802 )
  W( 99 ) = W( 99 ) + a*JVS( 803 )
  W( 100 ) = W( 100 ) + a*JVS( 804 )
  a = -W( 92 ) / JVS(          838  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 839 )
  W( 94 ) = W( 94 ) + a*JVS( 840 )
  W( 95 ) = W( 95 ) + a*JVS( 841 )
  W( 96 ) = W( 96 ) + a*JVS( 842 )
  W( 97 ) = W( 97 ) + a*JVS( 843 )
  W( 98 ) = W( 98 ) + a*JVS( 844 )
  W( 99 ) = W( 99 ) + a*JVS( 845 )
  W( 100 ) = W( 100 ) + a*JVS( 846 )
  a = -W( 93 ) / JVS(          880  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 881 )
  W( 95 ) = W( 95 ) + a*JVS( 882 )
  W( 96 ) = W( 96 ) + a*JVS( 883 )
  W( 97 ) = W( 97 ) + a*JVS( 884 )
  W( 98 ) = W( 98 ) + a*JVS( 885 )
  W( 99 ) = W( 99 ) + a*JVS( 886 )
  W( 100 ) = W( 100 ) + a*JVS( 887 )
  a = -W( 94 ) / JVS(          913  )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 914 )
  W( 96 ) = W( 96 ) + a*JVS( 915 )
  W( 97 ) = W( 97 ) + a*JVS( 916 )
  W( 98 ) = W( 98 ) + a*JVS( 917 )
  W( 99 ) = W( 99 ) + a*JVS( 918 )
  W( 100 ) = W( 100 ) + a*JVS( 919 )
  a = -W( 95 ) / JVS(          931  )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 932 )
  W( 97 ) = W( 97 ) + a*JVS( 933 )
  W( 98 ) = W( 98 ) + a*JVS( 934 )
  W( 99 ) = W( 99 ) + a*JVS( 935 )
  W( 100 ) = W( 100 ) + a*JVS( 936 )
  a = -W( 96 ) / JVS(          952  )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 953 )
  W( 98 ) = W( 98 ) + a*JVS( 954 )
  W( 99 ) = W( 99 ) + a*JVS( 955 )
  W( 100 ) = W( 100 ) + a*JVS( 956 )
  a = -W( 97 ) / JVS(         1004  )
  W( 97 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 1005 )
  W( 99 ) = W( 99 ) + a*JVS( 1006 )
  W( 100 ) = W( 100 ) + a*JVS( 1007 )
  a = -W( 98 ) / JVS(         1071  )
  W( 98 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 1072 )
  W( 100 ) = W( 100 ) + a*JVS( 1073 )
  JVS( 1074) = W( 35 )
  JVS( 1075) = W( 74 )
  JVS( 1076) = W( 75 )
  JVS( 1077) = W( 76 )
  JVS( 1078) = W( 77 )
  JVS( 1079) = W( 78 )
  JVS( 1080) = W( 80 )
  JVS( 1081) = W( 81 )
  JVS( 1082) = W( 84 )
  JVS( 1083) = W( 85 )
  JVS( 1084) = W( 86 )
  JVS( 1085) = W( 87 )
  JVS( 1086) = W( 88 )
  JVS( 1087) = W( 89 )
  JVS( 1088) = W( 90 )
  JVS( 1089) = W( 91 )
  JVS( 1090) = W( 92 )
  JVS( 1091) = W( 93 )
  JVS( 1092) = W( 94 )
  JVS( 1093) = W( 95 )
  JVS( 1094) = W( 96 )
  JVS( 1095) = W( 97 )
  JVS( 1096) = W( 98 )
  JVS( 1097) = W( 99 )
  JVS( 1098) = W( 100 )
  IF ( ABS(  JVS( 1123 )) < TINY(a) ) THEN
         IER = 100                                     
         RETURN
  END IF
   W( 43 ) = JVS( 1099 )
   W( 48 ) = JVS( 1100 )
   W( 73 ) = JVS( 1101 )
   W( 75 ) = JVS( 1102 )
   W( 76 ) = JVS( 1103 )
   W( 80 ) = JVS( 1104 )
   W( 81 ) = JVS( 1105 )
   W( 82 ) = JVS( 1106 )
   W( 83 ) = JVS( 1107 )
   W( 84 ) = JVS( 1108 )
   W( 86 ) = JVS( 1109 )
   W( 87 ) = JVS( 1110 )
   W( 88 ) = JVS( 1111 )
   W( 89 ) = JVS( 1112 )
   W( 90 ) = JVS( 1113 )
   W( 91 ) = JVS( 1114 )
   W( 92 ) = JVS( 1115 )
   W( 93 ) = JVS( 1116 )
   W( 94 ) = JVS( 1117 )
   W( 95 ) = JVS( 1118 )
   W( 96 ) = JVS( 1119 )
   W( 97 ) = JVS( 1120 )
   W( 98 ) = JVS( 1121 )
   W( 99 ) = JVS( 1122 )
   W( 100 ) = JVS( 1123 )
  a = -W( 43 ) / JVS(          226  )
  W( 43 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 227 )
  W( 100 ) = W( 100 ) + a*JVS( 228 )
  a = -W( 48 ) / JVS(          241  )
  W( 48 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 242 )
  W( 97 ) = W( 97 ) + a*JVS( 243 )
  W( 100 ) = W( 100 ) + a*JVS( 244 )
  a = -W( 73 ) / JVS(          423  )
  W( 73 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 424 )
  W( 76 ) = W( 76 ) + a*JVS( 425 )
  W( 80 ) = W( 80 ) + a*JVS( 426 )
  W( 81 ) = W( 81 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 86 ) = W( 86 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 88 ) = W( 88 ) + a*JVS( 432 )
  W( 89 ) = W( 89 ) + a*JVS( 433 )
  W( 90 ) = W( 90 ) + a*JVS( 434 )
  W( 91 ) = W( 91 ) + a*JVS( 435 )
  W( 92 ) = W( 92 ) + a*JVS( 436 )
  W( 93 ) = W( 93 ) + a*JVS( 437 )
  W( 94 ) = W( 94 ) + a*JVS( 438 )
  W( 95 ) = W( 95 ) + a*JVS( 439 )
  W( 96 ) = W( 96 ) + a*JVS( 440 )
  W( 97 ) = W( 97 ) + a*JVS( 441 )
  W( 98 ) = W( 98 ) + a*JVS( 442 )
  W( 99 ) = W( 99 ) + a*JVS( 443 )
  W( 100 ) = W( 100 ) + a*JVS( 444 )
  a = -W( 75 ) / JVS(          452  )
  W( 75 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 453 )
  W( 88 ) = W( 88 ) + a*JVS( 454 )
  W( 93 ) = W( 93 ) + a*JVS( 455 )
  W( 98 ) = W( 98 ) + a*JVS( 456 )
  a = -W( 76 ) / JVS(          457  )
  W( 76 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 458 )
  W( 88 ) = W( 88 ) + a*JVS( 459 )
  W( 93 ) = W( 93 ) + a*JVS( 460 )
  W( 98 ) = W( 98 ) + a*JVS( 461 )
  a = -W( 80 ) / JVS(          517  )
  W( 80 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 518 )
  W( 88 ) = W( 88 ) + a*JVS( 519 )
  W( 93 ) = W( 93 ) + a*JVS( 520 )
  W( 98 ) = W( 98 ) + a*JVS( 521 )
  a = -W( 81 ) / JVS(          524  )
  W( 81 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 525 )
  W( 88 ) = W( 88 ) + a*JVS( 526 )
  W( 93 ) = W( 93 ) + a*JVS( 527 )
  W( 98 ) = W( 98 ) + a*JVS( 528 )
  a = -W( 82 ) / JVS(          549  )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 550 )
  W( 84 ) = W( 84 ) + a*JVS( 551 )
  W( 86 ) = W( 86 ) + a*JVS( 552 )
  W( 87 ) = W( 87 ) + a*JVS( 553 )
  W( 88 ) = W( 88 ) + a*JVS( 554 )
  W( 89 ) = W( 89 ) + a*JVS( 555 )
  W( 90 ) = W( 90 ) + a*JVS( 556 )
  W( 91 ) = W( 91 ) + a*JVS( 557 )
  W( 92 ) = W( 92 ) + a*JVS( 558 )
  W( 93 ) = W( 93 ) + a*JVS( 559 )
  W( 94 ) = W( 94 ) + a*JVS( 560 )
  W( 95 ) = W( 95 ) + a*JVS( 561 )
  W( 96 ) = W( 96 ) + a*JVS( 562 )
  W( 97 ) = W( 97 ) + a*JVS( 563 )
  W( 98 ) = W( 98 ) + a*JVS( 564 )
  W( 99 ) = W( 99 ) + a*JVS( 565 )
  W( 100 ) = W( 100 ) + a*JVS( 566 )
  a = -W( 83 ) / JVS(          573  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 574 )
  W( 88 ) = W( 88 ) + a*JVS( 575 )
  W( 91 ) = W( 91 ) + a*JVS( 576 )
  W( 92 ) = W( 92 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 98 ) = W( 98 ) + a*JVS( 579 )
  W( 100 ) = W( 100 ) + a*JVS( 580 )
  a = -W( 84 ) / JVS(          591  )
  W( 84 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 592 )
  W( 92 ) = W( 92 ) + a*JVS( 593 )
  W( 93 ) = W( 93 ) + a*JVS( 594 )
  W( 98 ) = W( 98 ) + a*JVS( 595 )
  W( 100 ) = W( 100 ) + a*JVS( 596 )
  a = -W( 86 ) / JVS(          638  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 639 )
  W( 88 ) = W( 88 ) + a*JVS( 640 )
  W( 90 ) = W( 90 ) + a*JVS( 641 )
  W( 91 ) = W( 91 ) + a*JVS( 642 )
  W( 92 ) = W( 92 ) + a*JVS( 643 )
  W( 93 ) = W( 93 ) + a*JVS( 644 )
  W( 94 ) = W( 94 ) + a*JVS( 645 )
  W( 98 ) = W( 98 ) + a*JVS( 646 )
  W( 100 ) = W( 100 ) + a*JVS( 647 )
  a = -W( 87 ) / JVS(          660  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 661 )
  W( 89 ) = W( 89 ) + a*JVS( 662 )
  W( 90 ) = W( 90 ) + a*JVS( 663 )
  W( 91 ) = W( 91 ) + a*JVS( 664 )
  W( 92 ) = W( 92 ) + a*JVS( 665 )
  W( 93 ) = W( 93 ) + a*JVS( 666 )
  W( 94 ) = W( 94 ) + a*JVS( 667 )
  W( 95 ) = W( 95 ) + a*JVS( 668 )
  W( 98 ) = W( 98 ) + a*JVS( 669 )
  W( 99 ) = W( 99 ) + a*JVS( 670 )
  W( 100 ) = W( 100 ) + a*JVS( 671 )
  a = -W( 88 ) / JVS(          685  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 686 )
  W( 92 ) = W( 92 ) + a*JVS( 687 )
  W( 93 ) = W( 93 ) + a*JVS( 688 )
  W( 95 ) = W( 95 ) + a*JVS( 689 )
  W( 96 ) = W( 96 ) + a*JVS( 690 )
  W( 97 ) = W( 97 ) + a*JVS( 691 )
  W( 98 ) = W( 98 ) + a*JVS( 692 )
  W( 99 ) = W( 99 ) + a*JVS( 693 )
  W( 100 ) = W( 100 ) + a*JVS( 694 )
  a = -W( 89 ) / JVS(          717  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 718 )
  W( 91 ) = W( 91 ) + a*JVS( 719 )
  W( 92 ) = W( 92 ) + a*JVS( 720 )
  W( 93 ) = W( 93 ) + a*JVS( 721 )
  W( 94 ) = W( 94 ) + a*JVS( 722 )
  W( 95 ) = W( 95 ) + a*JVS( 723 )
  W( 96 ) = W( 96 ) + a*JVS( 724 )
  W( 97 ) = W( 97 ) + a*JVS( 725 )
  W( 98 ) = W( 98 ) + a*JVS( 726 )
  W( 99 ) = W( 99 ) + a*JVS( 727 )
  W( 100 ) = W( 100 ) + a*JVS( 728 )
  a = -W( 90 ) / JVS(          762  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 92 ) = W( 92 ) + a*JVS( 764 )
  W( 93 ) = W( 93 ) + a*JVS( 765 )
  W( 94 ) = W( 94 ) + a*JVS( 766 )
  W( 95 ) = W( 95 ) + a*JVS( 767 )
  W( 96 ) = W( 96 ) + a*JVS( 768 )
  W( 97 ) = W( 97 ) + a*JVS( 769 )
  W( 98 ) = W( 98 ) + a*JVS( 770 )
  W( 99 ) = W( 99 ) + a*JVS( 771 )
  W( 100 ) = W( 100 ) + a*JVS( 772 )
  a = -W( 91 ) / JVS(          795  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 796 )
  W( 93 ) = W( 93 ) + a*JVS( 797 )
  W( 94 ) = W( 94 ) + a*JVS( 798 )
  W( 95 ) = W( 95 ) + a*JVS( 799 )
  W( 96 ) = W( 96 ) + a*JVS( 800 )
  W( 97 ) = W( 97 ) + a*JVS( 801 )
  W( 98 ) = W( 98 ) + a*JVS( 802 )
  W( 99 ) = W( 99 ) + a*JVS( 803 )
  W( 100 ) = W( 100 ) + a*JVS( 804 )
  a = -W( 92 ) / JVS(          838  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 839 )
  W( 94 ) = W( 94 ) + a*JVS( 840 )
  W( 95 ) = W( 95 ) + a*JVS( 841 )
  W( 96 ) = W( 96 ) + a*JVS( 842 )
  W( 97 ) = W( 97 ) + a*JVS( 843 )
  W( 98 ) = W( 98 ) + a*JVS( 844 )
  W( 99 ) = W( 99 ) + a*JVS( 845 )
  W( 100 ) = W( 100 ) + a*JVS( 846 )
  a = -W( 93 ) / JVS(          880  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 881 )
  W( 95 ) = W( 95 ) + a*JVS( 882 )
  W( 96 ) = W( 96 ) + a*JVS( 883 )
  W( 97 ) = W( 97 ) + a*JVS( 884 )
  W( 98 ) = W( 98 ) + a*JVS( 885 )
  W( 99 ) = W( 99 ) + a*JVS( 886 )
  W( 100 ) = W( 100 ) + a*JVS( 887 )
  a = -W( 94 ) / JVS(          913  )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 914 )
  W( 96 ) = W( 96 ) + a*JVS( 915 )
  W( 97 ) = W( 97 ) + a*JVS( 916 )
  W( 98 ) = W( 98 ) + a*JVS( 917 )
  W( 99 ) = W( 99 ) + a*JVS( 918 )
  W( 100 ) = W( 100 ) + a*JVS( 919 )
  a = -W( 95 ) / JVS(          931  )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 932 )
  W( 97 ) = W( 97 ) + a*JVS( 933 )
  W( 98 ) = W( 98 ) + a*JVS( 934 )
  W( 99 ) = W( 99 ) + a*JVS( 935 )
  W( 100 ) = W( 100 ) + a*JVS( 936 )
  a = -W( 96 ) / JVS(          952  )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 953 )
  W( 98 ) = W( 98 ) + a*JVS( 954 )
  W( 99 ) = W( 99 ) + a*JVS( 955 )
  W( 100 ) = W( 100 ) + a*JVS( 956 )
  a = -W( 97 ) / JVS(         1004  )
  W( 97 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 1005 )
  W( 99 ) = W( 99 ) + a*JVS( 1006 )
  W( 100 ) = W( 100 ) + a*JVS( 1007 )
  a = -W( 98 ) / JVS(         1071  )
  W( 98 ) = -a
  W( 99 ) = W( 99 ) + a*JVS( 1072 )
  W( 100 ) = W( 100 ) + a*JVS( 1073 )
  a = -W( 99 ) / JVS(         1097  )
  W( 99 ) = -a
  W( 100 ) = W( 100 ) + a*JVS( 1098 )
  JVS( 1099) = W( 43 )
  JVS( 1100) = W( 48 )
  JVS( 1101) = W( 73 )
  JVS( 1102) = W( 75 )
  JVS( 1103) = W( 76 )
  JVS( 1104) = W( 80 )
  JVS( 1105) = W( 81 )
  JVS( 1106) = W( 82 )
  JVS( 1107) = W( 83 )
  JVS( 1108) = W( 84 )
  JVS( 1109) = W( 86 )
  JVS( 1110) = W( 87 )
  JVS( 1111) = W( 88 )
  JVS( 1112) = W( 89 )
  JVS( 1113) = W( 90 )
  JVS( 1114) = W( 91 )
  JVS( 1115) = W( 92 )
  JVS( 1116) = W( 93 )
  JVS( 1117) = W( 94 )
  JVS( 1118) = W( 95 )
  JVS( 1119) = W( 96 )
  JVS( 1120) = W( 97 )
  JVS( 1121) = W( 98 )
  JVS( 1122) = W( 99 )
  JVS( 1123) = W( 100 )
   
   END SUBROUTINE decomp_saprc99_mosaic_4bin_vbs2
 


END MODULE saprc99_mosaic_4bin_vbs2_Integrator
