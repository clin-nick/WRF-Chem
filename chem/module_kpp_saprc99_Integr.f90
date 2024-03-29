
































MODULE saprc99_Integrator

 USE saprc99_Parameters
 USE saprc99_Precision
 USE saprc99_JacobianSP

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

SUBROUTINE  saprc99_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE saprc99_Parameters

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

   CALL saprc99_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  saprc99_INTEGRATE


SUBROUTINE  saprc99_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE saprc99_Parameters

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
      CALL saprc99_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL saprc99_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = saprc99_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL saprc99_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL saprc99_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL saprc99_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL saprc99_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL saprc99_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL saprc99_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL saprc99_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL saprc99_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL saprc99_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL saprc99_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL saprc99_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL saprc99_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL saprc99_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL saprc99_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL saprc99_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  saprc99_ros_ErrorMsg(Code,T,H,IERR)



   USE saprc99_Precision

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

 END SUBROUTINE  saprc99_ros_ErrorMsg


 SUBROUTINE  saprc99_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL saprc99_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL saprc99_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL saprc99_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL saprc99_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL saprc99_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL saprc99_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL saprc99_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL saprc99_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL saprc99_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL saprc99_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL saprc99_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL saprc99_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL saprc99_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL saprc99_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL saprc99_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL saprc99_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL saprc99_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL saprc99_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL saprc99_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL saprc99_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = saprc99_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL saprc99_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  saprc99_ros_Integrator



  REAL(kind=dp) FUNCTION  saprc99_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    saprc99_ros_ErrorNorm = Err

  END FUNCTION  saprc99_ros_ErrorNorm



  SUBROUTINE saprc99_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL saprc99_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL saprc99_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL saprc99_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  saprc99_ros_FunTimeDeriv



  SUBROUTINE  saprc99_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL saprc99_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL saprc99_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL saprc99_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  saprc99_ros_PrepareMatrix



  SUBROUTINE  saprc99_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 



CALL decomp_saprc99 ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  saprc99_ros_Decomp



  SUBROUTINE  saprc99_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL saprc99_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  saprc99_ros_Solve




  SUBROUTINE  saprc99_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  saprc99_Ros2



  SUBROUTINE  saprc99_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_Ros3





  SUBROUTINE  saprc99_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_Ros4


  SUBROUTINE  saprc99_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_Rodas3


  SUBROUTINE  saprc99_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  saprc99_Rodas4




END SUBROUTINE  saprc99_Rosenbrock




SUBROUTINE  saprc99_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE saprc99_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL saprc99_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  saprc99_FunTemplate



SUBROUTINE  saprc99_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE saprc99_Parameters
 
 USE saprc99_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL saprc99_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  saprc99_JacTemplate

















SUBROUTINE saprc99_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(74)
  A(2) = RCT(2)*V(63)*F(2)
  A(3) = RCT(3)*V(63)*V(67)
  A(4) = RCT(4)*V(63)*V(73)*F(2)
  A(5) = RCT(5)*V(63)*V(74)
  A(6) = RCT(6)*V(63)*V(74)
  A(7) = RCT(7)*V(67)*V(73)
  A(8) = RCT(8)*V(67)*V(74)
  A(9) = RCT(9)*V(73)*V(75)
  A(10) = RCT(10)*V(73)*V(73)*F(2)
  A(11) = RCT(11)*V(74)*V(75)
  A(12) = RCT(12)*V(22)
  A(13) = RCT(13)*V(22)*F(1)
  A(14) = RCT(14)*V(74)*V(75)
  A(15) = RCT(15)*V(75)
  A(16) = RCT(16)*V(75)
  A(17) = RCT(17)*V(67)
  A(18) = RCT(18)*V(67)
  A(19) = RCT(19)*V(11)*F(1)
  A(20) = RCT(20)*V(11)*F(2)
  A(21) = RCT(21)*V(71)*V(73)
  A(22) = RCT(22)*V(23)
  A(23) = RCT(23)*V(23)
  A(24) = RCT(24)*V(23)*V(71)
  A(25) = RCT(25)*V(71)*V(74)
  A(26) = RCT(26)*V(71)*V(75)
  A(27) = RCT(27)*V(47)*V(71)
  A(28) = RCT(28)*V(47)
  A(29) = RCT(29)*V(45)*V(71)
  A(30) = RCT(30)*V(67)*V(71)
  A(31) = RCT(31)*V(73)*V(78)
  A(32) = RCT(32)*V(74)*V(78)
  A(33) = RCT(33)*V(32)
  A(34) = RCT(34)*V(32)
  A(35) = RCT(35)*V(32)*V(71)
  A(36) = RCT(36)*V(67)*V(78)
  A(37) = RCT(37)*V(78)*V(78)
  A(38) = RCT(38)*V(78)*V(78)*F(1)
  A(39) = RCT(39)*V(75)*V(78)
  A(40) = RCT(40)*V(75)*V(75)
  A(41) = RCT(41)*V(19)
  A(42) = RCT(42)*V(19)*V(71)
  A(43) = RCT(43)*V(71)*V(78)
  A(44) = RCT(44)*V(10)*V(71)
  A(45) = RCT(45)*V(71)*F(2)
  A(46) = RCT(46)*V(73)*V(76)
  A(47) = RCT(47)*V(76)*V(78)
  A(48) = RCT(48)*V(75)*V(76)
  A(49) = RCT(49)*V(76)*V(76)
  A(50) = RCT(50)*V(76)*V(76)
  A(51) = RCT(51)*V(73)*V(77)
  A(52) = RCT(52)*V(77)*V(78)
  A(53) = RCT(53)*V(75)*V(77)
  A(54) = RCT(54)*V(76)*V(77)
  A(55) = RCT(55)*V(77)*V(77)
  A(56) = RCT(56)*V(53)*V(73)
  A(57) = RCT(57)*V(53)*V(78)
  A(58) = RCT(58)*V(53)*V(75)
  A(59) = RCT(59)*V(53)*V(76)
  A(60) = RCT(60)*V(53)*V(77)
  A(62) = RCT(62)*V(69)*V(73)
  A(63) = RCT(63)*V(69)*V(78)
  A(64) = RCT(64)*V(69)*V(76)
  A(65) = RCT(65)*V(69)*V(75)
  A(66) = RCT(66)*V(69)*V(77)
  A(67) = RCT(67)*V(53)*V(69)
  A(68) = RCT(68)*V(69)*V(69)
  A(69) = RCT(69)*V(68)*V(74)
  A(70) = RCT(70)*V(15)
  A(71) = RCT(71)*V(68)*V(73)
  A(72) = RCT(72)*V(68)*V(78)
  A(73) = RCT(73)*V(68)*V(75)
  A(74) = RCT(74)*V(68)*V(76)
  A(75) = RCT(75)*V(68)*V(77)
  A(76) = RCT(76)*V(53)*V(68)
  A(77) = RCT(77)*V(68)*V(69)
  A(78) = RCT(78)*V(68)*V(68)
  A(79) = RCT(79)*V(72)*V(74)
  A(80) = RCT(80)*V(16)
  A(81) = RCT(81)*V(72)*V(73)
  A(82) = RCT(82)*V(72)*V(78)
  A(83) = RCT(83)*V(72)*V(75)
  A(84) = RCT(84)*V(72)*V(76)
  A(85) = RCT(85)*V(72)*V(77)
  A(86) = RCT(86)*V(53)*V(72)
  A(87) = RCT(87)*V(69)*V(72)
  A(88) = RCT(88)*V(68)*V(72)
  A(89) = RCT(89)*V(72)*V(72)
  A(90) = RCT(90)*V(74)*V(79)
  A(91) = RCT(91)*V(17)
  A(92) = RCT(92)*V(73)*V(79)
  A(93) = RCT(93)*V(78)*V(79)
  A(94) = RCT(94)*V(75)*V(79)
  A(95) = RCT(95)*V(76)*V(79)
  A(96) = RCT(96)*V(77)*V(79)
  A(97) = RCT(97)*V(53)*V(79)
  A(98) = RCT(98)*V(69)*V(79)
  A(99) = RCT(99)*V(68)*V(79)
  A(100) = RCT(100)*V(72)*V(79)
  A(101) = RCT(101)*V(79)*V(79)
  A(102) = RCT(102)*V(70)*V(74)
  A(103) = RCT(103)*V(18)
  A(104) = RCT(104)*V(70)*V(73)
  A(105) = RCT(105)*V(70)*V(78)
  A(106) = RCT(106)*V(70)*V(75)
  A(107) = RCT(107)*V(70)*V(76)
  A(108) = RCT(108)*V(70)*V(77)
  A(109) = RCT(109)*V(53)*V(70)
  A(110) = RCT(110)*V(69)*V(70)
  A(111) = RCT(111)*V(68)*V(70)
  A(112) = RCT(112)*V(70)*V(72)
  A(113) = RCT(113)*V(70)*V(79)
  A(114) = RCT(114)*V(70)*V(70)
  A(115) = RCT(115)*V(25)*V(74)
  A(116) = RCT(116)*V(25)
  A(117) = RCT(117)*V(51)*V(74)
  A(118) = RCT(118)*V(51)*V(78)
  A(119) = RCT(119)*V(51)
  A(120) = RCT(120)*V(31)*V(74)
  A(121) = RCT(121)*V(31)*V(78)
  A(122) = RCT(122)*V(31)
  A(123) = RCT(123)*V(61)
  A(124) = RCT(124)*V(61)
  A(125) = RCT(125)*V(61)*V(71)
  A(126) = RCT(126)*V(61)*V(78)
  A(127) = RCT(127)*V(30)
  A(128) = RCT(128)*V(30)*V(73)
  A(129) = RCT(129)*V(61)*V(75)
  A(130) = RCT(130)*V(60)*V(71)
  A(131) = RCT(131)*V(60)
  A(132) = RCT(132)*V(60)*V(75)
  A(133) = RCT(133)*V(64)*V(71)
  A(134) = RCT(134)*V(64)
  A(135) = RCT(135)*V(64)*V(75)
  A(136) = RCT(136)*V(46)*V(71)
  A(137) = RCT(137)*V(46)
  A(138) = RCT(138)*V(65)*V(71)
  A(139) = RCT(139)*V(65)
  A(140) = RCT(140)*V(28)*V(71)
  A(141) = RCT(141)*V(21)*V(71)
  A(142) = RCT(142)*V(29)*V(71)
  A(143) = RCT(143)*V(29)
  A(144) = RCT(144)*V(42)*V(71)
  A(145) = RCT(145)*V(42)
  A(146) = RCT(146)*V(49)
  A(147) = RCT(147)*V(49)
  A(148) = RCT(148)*V(49)*V(71)
  A(149) = RCT(149)*V(49)*V(75)
  A(150) = RCT(150)*V(44)
  A(151) = RCT(151)*V(44)*V(71)
  A(152) = RCT(152)*V(44)*V(75)
  A(153) = RCT(153)*V(14)
  A(154) = RCT(154)*V(43)*V(71)
  A(155) = RCT(155)*V(43)*V(75)
  A(156) = RCT(156)*V(37)*V(71)
  A(157) = RCT(157)*V(37)*V(75)
  A(158) = RCT(158)*V(40)*V(75)
  A(159) = RCT(159)*V(41)*V(71)
  A(160) = RCT(160)*V(41)
  A(161) = RCT(161)*V(41)*V(75)
  A(162) = RCT(162)*V(55)*V(71)
  A(163) = RCT(163)*V(55)*V(67)
  A(164) = RCT(164)*V(55)*V(75)
  A(165) = RCT(165)*V(55)*V(63)
  A(166) = RCT(166)*V(55)
  A(167) = RCT(167)*V(59)*V(71)
  A(168) = RCT(168)*V(59)*V(67)
  A(169) = RCT(169)*V(59)*V(63)
  A(170) = RCT(170)*V(59)
  A(171) = RCT(171)*V(57)*V(71)
  A(172) = RCT(172)*V(57)*V(67)
  A(173) = RCT(173)*V(57)*V(75)
  A(174) = RCT(174)*V(57)
  A(175) = RCT(175)*V(66)*V(71)
  A(176) = RCT(176)*V(66)
  A(177) = RCT(177)*V(62)*V(71)
  A(178) = RCT(178)*V(62)
  A(179) = RCT(179)*V(38)*V(71)
  A(180) = RCT(180)*V(38)*V(67)
  A(181) = RCT(181)*V(35)*V(71)
  A(182) = RCT(182)*V(35)
  A(183) = RCT(183)*V(36)*V(71)
  A(184) = RCT(184)*V(36)
  A(185) = RCT(185)*V(12)*V(71)
  A(186) = RCT(186)*V(48)*V(71)
  A(187) = RCT(187)*V(48)*V(67)
  A(188) = RCT(188)*V(48)*V(75)
  A(189) = RCT(189)*V(48)*V(63)
  A(190) = RCT(190)*V(52)*V(71)
  A(191) = RCT(191)*V(52)*V(67)
  A(192) = RCT(192)*V(52)*V(75)
  A(193) = RCT(193)*V(52)*V(63)
  A(194) = RCT(194)*V(54)*V(71)
  A(195) = RCT(195)*V(54)*V(67)
  A(196) = RCT(196)*V(54)*V(75)
  A(197) = RCT(197)*V(54)*V(63)
  A(198) = RCT(198)*V(13)*V(71)
  A(199) = RCT(199)*V(20)*V(71)
  A(200) = RCT(200)*V(39)*V(71)
  A(201) = RCT(201)*V(24)*V(71)
  A(202) = RCT(202)*V(33)*V(71)
  A(203) = RCT(203)*V(26)*V(71)
  A(204) = RCT(204)*V(34)*V(71)
  A(205) = RCT(205)*V(27)*V(71)
  A(206) = RCT(206)*V(56)*V(71)
  A(207) = RCT(207)*V(56)*V(67)
  A(208) = RCT(208)*V(56)*V(75)
  A(209) = RCT(209)*V(56)*V(63)
  A(210) = RCT(210)*V(58)*V(71)
  A(211) = RCT(211)*V(58)*V(67)
  A(212) = RCT(212)*V(58)*V(75)
  A(213) = RCT(213)*V(58)*V(63)
  A(214) = RCT(214)*V(39)*V(67)
  A(215) = RCT(215)*V(50)*V(71)
  A(216) = RCT(216)*V(50)*V(67)
  A(217) = RCT(217)*V(50)*V(75)
  A(218) = RCT(218)*V(50)*V(63)
  A(219) = RCT(219)*V(10)
  A(220) = RCT(220)*V(78)
  A(221) = RCT(221)*V(10)
  A(222) = RCT(222)*V(1)
  A(223) = RCT(223)*V(47)
  A(224) = RCT(224)*V(19)
  A(225) = RCT(225)*V(2)


  Vdot(1) = A(44)+A(219)-A(222)
  Vdot(2) = 0.5*A(214)+0.135*A(216)-A(225)
  Vdot(3) = A(128)+0.333*A(163)+0.351*A(168)+0.1*A(172)+0.37*A(187)+0.204*A(191)+0.103*A(195)+0.297*A(200)+0.185*A(207)&
              &+0.073*A(211)+0.185*A(216)
  Vdot(4) = 0.25*A(72)+A(74)+A(75)+A(77)+0.05*A(207)+0.129*A(211)+0.17*A(216)
  Vdot(5) = 0.25*A(82)+A(84)+A(85)+A(87)+0.25*A(93)+A(95)+A(96)+A(98)+0.25*A(105)+A(107)+A(108)+2*A(110)+0.372*A(172)&
              &+0.15*A(191)+0.189*A(195)+0.119*A(207)+0.247*A(211)
  Vdot(6) = 0.75*A(72)
  Vdot(7) = 0.75*A(82)+0.75*A(93)+0.75*A(105)
  Vdot(8) = 2*A(120)+A(217)
  Vdot(9) = 6*A(120)+7*A(160)+0.048*A(215)+0.07*A(216)+2.693*A(217)+0.55*A(218)
  Vdot(10) = -A(44)-A(219)-A(221)
  Vdot(11) = A(18)-A(19)-A(20)
  Vdot(12) = -A(185)
  Vdot(13) = -A(198)
  Vdot(14) = -A(153)+0.031*A(195)+0.087*A(205)
  Vdot(15) = A(69)-A(70)
  Vdot(16) = A(79)-A(80)
  Vdot(17) = A(90)-A(91)
  Vdot(18) = A(102)-A(103)
  Vdot(19) = A(37)+A(38)-A(41)-A(42)-A(224)
  Vdot(20) = -A(199)
  Vdot(21) = -A(141)
  Vdot(22) = A(11)-A(12)-A(13)
  Vdot(23) = A(21)-A(22)-A(23)-A(24)
  Vdot(24) = -A(201)
  Vdot(25) = -A(115)-A(116)+0.236*A(201)
  Vdot(26) = -A(203)
  Vdot(27) = -A(205)
  Vdot(28) = A(49)+0.25*A(54)+0.25*A(64)-A(140)
  Vdot(29) = A(47)-A(142)-A(143)
  Vdot(30) = A(126)-A(127)-A(128)
  Vdot(31) = -A(120)-A(121)-A(122)+A(158)
  Vdot(32) = A(32)-A(33)-A(34)-A(35)
  Vdot(33) = -A(202)
  Vdot(34) = -A(204)
  Vdot(35) = -A(181)-A(182)+0.108*A(204)+0.099*A(205)
  Vdot(36) = -A(183)-A(184)+0.051*A(204)+0.093*A(205)
  Vdot(37) = -A(156)-A(157)+0.207*A(204)+0.187*A(205)
  Vdot(38) = -A(179)-A(180)+0.491*A(204)+0.561*A(205)
  Vdot(39) = -A(200)-A(214)
  Vdot(40) = A(117)+A(121)+A(122)-A(158)
  Vdot(41) = -A(159)-A(160)-A(161)+0.059*A(204)+0.05*A(205)+0.061*A(210)+0.042*A(211)+0.015*A(212)
  Vdot(42) = A(52)+A(63)-A(144)-A(145)
  Vdot(43) = A(118)+A(119)-A(154)-A(155)+0.017*A(204)
  Vdot(44) = -A(150)-A(151)-A(152)+0.23*A(156)+0.084*A(162)+0.9*A(163)+0.3*A(167)+0.95*A(168)+0.174*A(171)+0.742*A(172)&
               &+0.008*A(173)+0.5*A(182)+0.5*A(184)+0.119*A(204)+0.287*A(205)
  Vdot(45) = -A(29)+A(123)+A(124)+A(125)+A(129)+A(131)+0.034*A(133)+A(134)+2*A(146)+A(147)+1.26*A(148)+1.26*A(149)&
               &+A(150)+A(151)+A(152)+0.416*A(162)+0.45*A(163)+0.5*A(164)+0.67*A(166)+0.475*A(168)+0.7*A(170)+0.336*A(171)&
               &+0.498*A(172)+0.572*A(173)+1.233*A(174)+A(179)+1.5*A(180)+A(182)+A(184)+0.5*A(187)+0.491*A(189)+0.275*A(191)&
               &+0.157*A(195)+0.393*A(200)+0.002*A(202)+0.345*A(207)+0.265*A(211)+0.012*A(213)+1.5*A(214)+0.51*A(216)
  Vdot(46) = A(116)-A(136)-A(137)+0.006*A(177)+0.02*A(178)+0.13*A(195)+0.704*A(199)+0.024*A(201)+0.452*A(202)+0.072&
               &*A(203)+0.005*A(206)+0.001*A(207)+0.024*A(208)+0.127*A(210)+0.045*A(211)+0.102*A(212)
  Vdot(47) = 2*A(13)+A(25)-A(27)-A(28)+0.2*A(39)+A(129)+A(132)+A(135)+A(149)+A(152)+A(155)+A(157)+A(158)+A(161)+0.5&
               &*A(164)+0.15*A(173)-A(223)
  Vdot(48) = -A(186)-A(187)-A(188)-A(189)
  Vdot(49) = -A(146)-A(147)-A(148)-A(149)+0.23*A(154)+0.15*A(171)+0.023*A(172)+A(180)+0.5*A(182)+0.5*A(184)+0.009*A(189)&
               &+0.001*A(195)+0.607*A(200)+0.118*A(204)+0.097*A(205)
  Vdot(50) = -A(215)-A(216)-A(217)-A(218)
  Vdot(51) = A(92)+A(94)+A(99)+A(100)+2*A(101)+A(113)-A(117)-A(118)-A(119)+0.24*A(154)+A(155)+0.24*A(156)+A(157)
  Vdot(52) = -A(190)-A(191)-A(192)-A(193)
  Vdot(53) = -A(56)-A(57)-A(58)-A(59)-A(60)-A(67)-A(76)-A(86)+A(92)+A(94)-A(97)+A(99)+A(100)+2*A(101)-A(109)+A(113)&
               &+A(136)+0.616*A(138)+0.675*A(167)+0.515*A(176)+0.596*A(177)+0.152*A(178)+A(181)+A(182)+A(183)+A(184)+0.079&
               &*A(190)+0.126*A(191)+0.187*A(192)+0.24*A(193)+0.5*A(194)+0.729*A(195)+0.75*A(196)+0.559*A(201)+0.936*A(202)&
               &+0.948*A(203)+0.205*A(206)+0.488*A(208)+0.001*A(210)+0.137*A(211)+0.711*A(212)
  Vdot(54) = -A(194)-A(195)-A(196)-A(197)
  Vdot(55) = -A(162)-A(163)-A(164)-A(165)-A(166)+0.23*A(190)+0.39*A(191)+0.025*A(210)+0.026*A(211)+0.012*A(213)
  Vdot(56) = -A(206)-A(207)-A(208)-A(209)
  Vdot(57) = -A(171)-A(172)-A(173)-A(174)+0.357*A(190)+0.936*A(192)+0.025*A(210)
  Vdot(58) = -A(210)-A(211)-A(212)-A(213)
  Vdot(59) = -A(167)-A(168)-A(169)-A(170)+0.32*A(190)+0.16*A(191)+0.019*A(211)+0.048*A(212)
  Vdot(60) = A(81)+A(83)+A(88)+2*A(89)+A(100)+A(112)-A(130)-A(131)-A(132)+0.034*A(133)+A(134)+0.482*A(138)+A(139)+0.96&
               &*A(141)+0.129*A(171)+0.047*A(172)+0.467*A(174)+0.084*A(175)+0.246*A(176)+0.439*A(177)+0.431*A(178)+0.195&
               &*A(186)+0.25*A(189)+A(198)+0.445*A(201)+0.455*A(202)+0.099*A(203)+0.294*A(206)+0.154*A(207)+0.009*A(208)&
               &+0.732*A(210)+0.456*A(211)+0.507*A(212)+0.984*A(215)+0.5*A(216)
  Vdot(61) = A(46)+A(48)+A(49)+2*A(50)+0.75*A(54)+0.75*A(64)+A(74)+A(84)+A(95)+A(104)+A(106)+A(107)+A(111)+A(112)+A(113)&
               &+2*A(114)-A(123)-A(124)-A(125)-A(126)+A(127)-A(129)+A(136)+0.115*A(138)+A(140)+0.081*A(141)+0.35*A(142)&
               &+A(143)+A(147)+0.084*A(162)+0.2*A(163)+0.67*A(166)+0.3*A(167)+0.1*A(168)+0.055*A(171)+0.125*A(172)+0.227&
               &*A(173)+0.3*A(174)+0.213*A(175)+0.506*A(176)+0.01*A(177)+0.134*A(178)+1.61*A(186)+A(187)+0.191*A(189)+0.624&
               &*A(190)+0.592*A(191)+0.24*A(193)+0.276*A(194)+0.235*A(195)+0.096*A(200)+0.026*A(201)+0.024*A(202)+0.026&
               &*A(203)+0.732*A(206)+0.5*A(207)+0.244*A(210)+0.269*A(211)+0.079*A(212)+0.984*A(215)+0.5*A(216)
  Vdot(62) = A(62)+A(115)+0.572*A(173)-0.69*A(177)-A(178)+0.276*A(196)+0.511*A(208)+0.321*A(212)
  Vdot(63) = A(1)-A(2)-A(3)-A(4)-A(5)-A(6)+A(16)+A(17)+A(20)-A(165)-A(169)-A(189)-A(193)-A(197)-A(209)-A(213)-A(218)
  Vdot(64) = -A(133)-A(134)-A(135)+0.37*A(138)+A(144)+A(145)+A(165)+0.675*A(167)+0.45*A(169)+0.013*A(171)+0.218*A(173)&
               &+0.558*A(175)+0.71*A(176)+0.213*A(177)+0.147*A(178)+A(179)+A(181)+A(183)+A(188)+0.474*A(194)+0.205*A(195)&
               &+0.474*A(196)+0.147*A(197)+0.261*A(199)+0.122*A(201)+0.244*A(202)+0.204*A(203)+0.497*A(206)+0.363*A(207)&
               &+0.037*A(208)+0.45*A(209)+0.511*A(210)+0.305*A(211)+0.151*A(212)+0.069*A(213)+0.45*A(218)
  Vdot(65) = 0.5*A(64)+A(65)+0.5*A(66)+A(68)-A(138)-A(139)+0.416*A(162)+0.55*A(169)+0.15*A(171)+0.21*A(172)+0.233*A(174)&
               &+0.115*A(175)+0.177*A(177)+0.243*A(178)+0.332*A(201)+0.11*A(202)+0.089*A(203)+0.437*A(209)+0.072*A(210)&
               &+0.026*A(211)+0.001*A(212)+0.659*A(213)+0.55*A(218)
  Vdot(66) = 0.5*A(64)+0.5*A(66)+A(68)+A(77)+A(87)+A(98)+0.7*A(170)+0.332*A(171)-0.671*A(175)-A(176)+0.048*A(177)+0.435&
               &*A(178)+0.1*A(191)+0.75*A(193)+0.276*A(194)+0.276*A(195)+0.853*A(197)+0.125*A(202)+0.417*A(203)+0.055*A(204)&
               &+0.119*A(206)+0.215*A(207)+0.113*A(209)+0.043*A(211)+0.259*A(213)
  Vdot(67) = A(2)-A(3)-A(7)-A(8)-A(17)-A(18)-A(30)-A(36)+0.25*A(72)+0.25*A(82)+0.25*A(93)+0.25*A(105)-A(163)-A(168)&
               &-A(172)-A(180)-A(187)-A(191)-A(195)-A(207)-A(211)-A(214)-A(216)
  Vdot(68) = -A(69)+A(70)-A(71)-A(72)-A(73)-A(74)-A(75)-A(77)-2*A(78)-A(88)-A(99)+A(104)+A(106)+A(112)+A(113)+2*A(114)&
               &+A(130)+A(132)+A(136)+A(137)+0.492*A(138)+A(139)+A(150)+A(151)+A(152)+2*A(153)+0.67*A(166)+0.675*A(167)&
               &+0.467*A(174)+0.029*A(175)+0.667*A(176)+A(181)+0.5*A(182)+A(183)+0.5*A(184)+0.123*A(195)+0.011*A(202)+0.137&
               &*A(211)
  Vdot(69) = -A(62)-A(63)-A(64)-A(65)-A(66)-2*A(68)-A(77)-A(87)-A(98)-A(110)+0.001*A(133)+0.042*A(138)+0.025*A(167)&
               &+0.041*A(171)+0.051*A(173)+0.07*A(175)+0.04*A(176)+0.173*A(177)+0.095*A(178)+0.093*A(190)+0.008*A(191)+0.064&
               &*A(192)+0.01*A(193)+0.25*A(194)+0.18*A(195)+0.25*A(196)+0.035*A(199)+0.07*A(201)+0.143*A(202)+0.347*A(203)&
               &+0.011*A(204)+0.009*A(205)+0.09*A(206)+0.001*A(207)+0.176*A(208)+0.082*A(210)+0.002*A(211)+0.136*A(212)&
               &+0.001*A(213)+0.016*A(215)+0.051*A(217)
  Vdot(70) = -A(102)+A(103)-A(104)-A(105)-A(106)-A(107)-A(108)-A(110)-A(111)-A(112)-A(113)-2*A(114)+0.5*A(162)+0.5&
               &*A(164)+0.33*A(166)+0.3*A(170)+0.289*A(171)+0.15*A(173)+0.192*A(191)+0.24*A(193)
  Vdot(71) = 2*A(19)-A(21)+A(22)-A(24)-A(25)-A(26)-A(27)+A(28)-A(29)-A(30)+A(31)+0.39*A(34)-A(35)+A(36)+0.8*A(39)+2&
               &*A(41)-A(42)-A(43)-A(44)-A(45)-A(125)-A(130)-A(133)-A(136)-A(138)-A(140)-A(141)-0.65*A(142)+A(143)-0.34&
               &*A(144)+A(145)-A(148)-A(151)-A(154)-A(156)-A(159)-A(162)+0.208*A(163)+0.33*A(166)-A(167)+0.164*A(168)-A(171)&
               &+0.285*A(172)-A(175)-A(177)-A(179)+0.5*A(180)-A(181)-A(183)-A(185)-A(186)+0.12*A(187)-A(190)+0.266*A(191)&
               &-A(194)+0.567*A(195)-A(198)-A(199)-0.397*A(200)-A(201)-A(202)-A(203)-A(204)-A(205)-A(206)+0.155*A(207)&
               &-A(210)+0.378*A(211)+0.5*A(214)-A(215)+0.32*A(216)
  Vdot(72) = -A(79)+A(80)-A(81)-A(82)-A(83)-A(84)-A(85)-A(87)-A(88)-2*A(89)-A(100)-A(112)+0.965*A(133)+A(135)+0.096&
               &*A(138)+0.37*A(148)+0.37*A(149)+0.1*A(163)+0.05*A(168)+0.048*A(172)+0.3*A(174)+0.049*A(175)+0.333*A(176)&
               &+0.201*A(195)+0.006*A(211)
  Vdot(73) = A(1)-A(4)+A(5)-A(7)-A(9)-2*A(10)+A(14)+A(15)-A(21)+A(22)-A(31)-A(46)-A(51)-A(56)-A(62)-A(71)-A(81)-A(92)&
               &-A(104)-A(128)
  Vdot(74) = -A(1)+A(4)-A(5)-A(6)+A(7)-A(8)+2*A(9)+2*A(10)-A(11)+A(12)+A(16)+A(23)+A(24)-A(25)+A(26)+A(28)+A(31)-A(32)&
               &+A(33)+0.61*A(34)+A(35)+0.8*A(39)+2*A(40)+A(46)+A(48)+A(51)+A(53)+A(56)+A(58)+A(65)-A(69)+A(70)+A(71)+A(73)&
               &-A(79)+A(80)+A(81)+A(83)-A(90)+A(91)+A(92)+A(94)-A(102)+A(103)+A(104)+A(106)-A(115)-A(117)-A(120)+A(128)&
               &+0.338*A(177)+A(178)+0.187*A(192)+0.474*A(196)+0.391*A(212)
  Vdot(75) = A(6)+A(8)-A(9)-A(11)+A(12)-A(14)-A(15)-A(16)-A(26)+A(27)+0.39*A(34)-A(39)-2*A(40)-A(48)-A(53)-A(58)-A(65)&
               &-A(73)-A(83)-A(94)-A(106)-A(129)-A(132)-A(135)-A(149)-A(152)-A(155)-A(157)-A(158)-A(161)-A(164)-A(173)&
               &-A(188)-A(192)-A(196)-A(208)-A(212)-A(217)
  Vdot(76) = -A(46)-A(47)-A(48)-2*A(49)-2*A(50)-A(54)-A(64)+A(71)+A(73)-A(74)+2*A(78)-A(84)+A(88)-A(95)+A(99)-A(107)&
               &+A(111)+A(116)+A(131)+A(137)+0.65*A(142)+0.3*A(170)+A(185)+0.3*A(189)+0.25*A(193)+0.011*A(202)+0.076*A(207)&
               &+0.197*A(211)+0.03*A(212)+0.26*A(216)
  Vdot(77) = -A(51)-A(52)-A(53)-A(54)-2*A(55)-A(66)-A(75)+A(81)+A(83)-A(85)+A(88)+2*A(89)-A(96)+A(100)-A(108)+A(112)&
               &+0.034*A(133)+A(134)+0.37*A(138)+A(139)+0.05*A(141)+0.34*A(144)+0.76*A(154)+0.76*A(156)+0.5*A(162)+0.1&
               &*A(163)+0.5*A(164)+0.33*A(166)+0.3*A(167)+0.05*A(168)+0.67*A(171)+0.048*A(172)+0.799*A(173)+0.473*A(175)&
               &+0.96*A(176)+0.376*A(177)+0.564*A(178)+A(179)+A(182)+A(184)+A(186)+A(188)+0.2*A(189)+0.907*A(190)+0.066&
               &*A(191)+0.749*A(192)+0.75*A(194)+0.031*A(195)+0.276*A(196)+A(198)+0.965*A(199)+0.1*A(200)+0.695*A(201)+0.835&
               &*A(202)+0.653*A(203)+0.765*A(204)+0.804*A(205)+0.91*A(206)+0.022*A(207)+0.824*A(208)+0.918*A(210)+0.033&
               &*A(211)+0.442*A(212)+0.012*A(213)+0.984*A(215)+0.949*A(217)
  Vdot(78) = A(23)+A(26)+A(29)+A(30)-A(31)-A(32)+A(33)+0.61*A(34)-A(36)-2*A(37)-2*A(38)-A(39)+A(42)-A(43)+A(44)+A(45)&
               &+A(46)-A(47)+A(48)+2*A(50)+A(51)-A(52)+A(53)+A(54)+A(55)-A(63)+A(64)+A(65)+A(66)+A(68)-A(72)-A(82)-A(93)&
               &-A(105)-A(118)-A(121)+2*A(123)+A(125)-A(126)+A(127)+A(128)+A(129)+A(131)+A(134)+A(140)+0.95*A(141)+A(143)&
               &+A(145)+2*A(146)+0.63*A(148)+0.63*A(149)+A(150)+0.008*A(163)+0.34*A(166)+0.064*A(168)+0.4*A(172)+1.233&
               &*A(174)+0.379*A(175)+0.113*A(177)+0.341*A(178)+1.5*A(180)+0.5*A(182)+0.5*A(184)+0.12*A(187)+0.5*A(189)+0.033&
               &*A(195)+0.297*A(200)+0.224*A(204)+0.187*A(205)+0.056*A(207)+0.003*A(211)+0.013*A(213)+1.5*A(214)+0.06*A(216)&
               &-A(220)
  Vdot(79) = -A(90)+A(91)-A(92)-A(93)-A(94)-A(95)-A(96)-A(98)-A(99)-A(100)-2*A(101)-A(113)+A(159)+A(161)
      
END SUBROUTINE saprc99_Fun
















SUBROUTINE saprc99_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(74)
  IRR(2) = RCT(2)*V(63)*F(2)
  IRR(3) = RCT(3)*V(63)*V(67)
  IRR(4) = RCT(4)*V(63)*V(73)*F(2)
  IRR(5) = RCT(5)*V(63)*V(74)
  IRR(6) = RCT(6)*V(63)*V(74)
  IRR(7) = RCT(7)*V(67)*V(73)
  IRR(8) = RCT(8)*V(67)*V(74)
  IRR(9) = RCT(9)*V(73)*V(75)
  IRR(10) = RCT(10)*V(73)*V(73)*F(2)
  IRR(11) = RCT(11)*V(74)*V(75)
  IRR(12) = RCT(12)*V(22)
  IRR(13) = RCT(13)*V(22)*F(1)
  IRR(14) = RCT(14)*V(74)*V(75)
  IRR(15) = RCT(15)*V(75)
  IRR(16) = RCT(16)*V(75)
  IRR(17) = RCT(17)*V(67)
  IRR(18) = RCT(18)*V(67)
  IRR(19) = RCT(19)*V(11)*F(1)
  IRR(20) = RCT(20)*V(11)*F(2)
  IRR(21) = RCT(21)*V(71)*V(73)
  IRR(22) = RCT(22)*V(23)
  IRR(23) = RCT(23)*V(23)
  IRR(24) = RCT(24)*V(23)*V(71)
  IRR(25) = RCT(25)*V(71)*V(74)
  IRR(26) = RCT(26)*V(71)*V(75)
  IRR(27) = RCT(27)*V(47)*V(71)
  IRR(28) = RCT(28)*V(47)
  IRR(29) = RCT(29)*V(45)*V(71)
  IRR(30) = RCT(30)*V(67)*V(71)
  IRR(31) = RCT(31)*V(73)*V(78)
  IRR(32) = RCT(32)*V(74)*V(78)
  IRR(33) = RCT(33)*V(32)
  IRR(34) = RCT(34)*V(32)
  IRR(35) = RCT(35)*V(32)*V(71)
  IRR(36) = RCT(36)*V(67)*V(78)
  IRR(37) = RCT(37)*V(78)*V(78)
  IRR(38) = RCT(38)*V(78)*V(78)*F(1)
  IRR(39) = RCT(39)*V(75)*V(78)
  IRR(40) = RCT(40)*V(75)*V(75)
  IRR(41) = RCT(41)*V(19)
  IRR(42) = RCT(42)*V(19)*V(71)
  IRR(43) = RCT(43)*V(71)*V(78)
  IRR(44) = RCT(44)*V(10)*V(71)
  IRR(45) = RCT(45)*V(71)*F(2)
  IRR(46) = RCT(46)*V(73)*V(76)
  IRR(47) = RCT(47)*V(76)*V(78)
  IRR(48) = RCT(48)*V(75)*V(76)
  IRR(49) = RCT(49)*V(76)*V(76)
  IRR(50) = RCT(50)*V(76)*V(76)
  IRR(51) = RCT(51)*V(73)*V(77)
  IRR(52) = RCT(52)*V(77)*V(78)
  IRR(53) = RCT(53)*V(75)*V(77)
  IRR(54) = RCT(54)*V(76)*V(77)
  IRR(55) = RCT(55)*V(77)*V(77)
  IRR(56) = RCT(56)*V(53)*V(73)
  IRR(57) = RCT(57)*V(53)*V(78)
  IRR(58) = RCT(58)*V(53)*V(75)
  IRR(59) = RCT(59)*V(53)*V(76)
  IRR(60) = RCT(60)*V(53)*V(77)
  IRR(62) = RCT(62)*V(69)*V(73)
  IRR(63) = RCT(63)*V(69)*V(78)
  IRR(64) = RCT(64)*V(69)*V(76)
  IRR(65) = RCT(65)*V(69)*V(75)
  IRR(66) = RCT(66)*V(69)*V(77)
  IRR(67) = RCT(67)*V(53)*V(69)
  IRR(68) = RCT(68)*V(69)*V(69)
  IRR(69) = RCT(69)*V(68)*V(74)
  IRR(70) = RCT(70)*V(15)
  IRR(71) = RCT(71)*V(68)*V(73)
  IRR(72) = RCT(72)*V(68)*V(78)
  IRR(73) = RCT(73)*V(68)*V(75)
  IRR(74) = RCT(74)*V(68)*V(76)
  IRR(75) = RCT(75)*V(68)*V(77)
  IRR(76) = RCT(76)*V(53)*V(68)
  IRR(77) = RCT(77)*V(68)*V(69)
  IRR(78) = RCT(78)*V(68)*V(68)
  IRR(79) = RCT(79)*V(72)*V(74)
  IRR(80) = RCT(80)*V(16)
  IRR(81) = RCT(81)*V(72)*V(73)
  IRR(82) = RCT(82)*V(72)*V(78)
  IRR(83) = RCT(83)*V(72)*V(75)
  IRR(84) = RCT(84)*V(72)*V(76)
  IRR(85) = RCT(85)*V(72)*V(77)
  IRR(86) = RCT(86)*V(53)*V(72)
  IRR(87) = RCT(87)*V(69)*V(72)
  IRR(88) = RCT(88)*V(68)*V(72)
  IRR(89) = RCT(89)*V(72)*V(72)
  IRR(90) = RCT(90)*V(74)*V(79)
  IRR(91) = RCT(91)*V(17)
  IRR(92) = RCT(92)*V(73)*V(79)
  IRR(93) = RCT(93)*V(78)*V(79)
  IRR(94) = RCT(94)*V(75)*V(79)
  IRR(95) = RCT(95)*V(76)*V(79)
  IRR(96) = RCT(96)*V(77)*V(79)
  IRR(97) = RCT(97)*V(53)*V(79)
  IRR(98) = RCT(98)*V(69)*V(79)
  IRR(99) = RCT(99)*V(68)*V(79)
  IRR(100) = RCT(100)*V(72)*V(79)
  IRR(101) = RCT(101)*V(79)*V(79)
  IRR(102) = RCT(102)*V(70)*V(74)
  IRR(103) = RCT(103)*V(18)
  IRR(104) = RCT(104)*V(70)*V(73)
  IRR(105) = RCT(105)*V(70)*V(78)
  IRR(106) = RCT(106)*V(70)*V(75)
  IRR(107) = RCT(107)*V(70)*V(76)
  IRR(108) = RCT(108)*V(70)*V(77)
  IRR(109) = RCT(109)*V(53)*V(70)
  IRR(110) = RCT(110)*V(69)*V(70)
  IRR(111) = RCT(111)*V(68)*V(70)
  IRR(112) = RCT(112)*V(70)*V(72)
  IRR(113) = RCT(113)*V(70)*V(79)
  IRR(114) = RCT(114)*V(70)*V(70)
  IRR(115) = RCT(115)*V(25)*V(74)
  IRR(116) = RCT(116)*V(25)
  IRR(117) = RCT(117)*V(51)*V(74)
  IRR(118) = RCT(118)*V(51)*V(78)
  IRR(119) = RCT(119)*V(51)
  IRR(120) = RCT(120)*V(31)*V(74)
  IRR(121) = RCT(121)*V(31)*V(78)
  IRR(122) = RCT(122)*V(31)
  IRR(123) = RCT(123)*V(61)
  IRR(124) = RCT(124)*V(61)
  IRR(125) = RCT(125)*V(61)*V(71)
  IRR(126) = RCT(126)*V(61)*V(78)
  IRR(127) = RCT(127)*V(30)
  IRR(128) = RCT(128)*V(30)*V(73)
  IRR(129) = RCT(129)*V(61)*V(75)
  IRR(130) = RCT(130)*V(60)*V(71)
  IRR(131) = RCT(131)*V(60)
  IRR(132) = RCT(132)*V(60)*V(75)
  IRR(133) = RCT(133)*V(64)*V(71)
  IRR(134) = RCT(134)*V(64)
  IRR(135) = RCT(135)*V(64)*V(75)
  IRR(136) = RCT(136)*V(46)*V(71)
  IRR(137) = RCT(137)*V(46)
  IRR(138) = RCT(138)*V(65)*V(71)
  IRR(139) = RCT(139)*V(65)
  IRR(140) = RCT(140)*V(28)*V(71)
  IRR(141) = RCT(141)*V(21)*V(71)
  IRR(142) = RCT(142)*V(29)*V(71)
  IRR(143) = RCT(143)*V(29)
  IRR(144) = RCT(144)*V(42)*V(71)
  IRR(145) = RCT(145)*V(42)
  IRR(146) = RCT(146)*V(49)
  IRR(147) = RCT(147)*V(49)
  IRR(148) = RCT(148)*V(49)*V(71)
  IRR(149) = RCT(149)*V(49)*V(75)
  IRR(150) = RCT(150)*V(44)
  IRR(151) = RCT(151)*V(44)*V(71)
  IRR(152) = RCT(152)*V(44)*V(75)
  IRR(153) = RCT(153)*V(14)
  IRR(154) = RCT(154)*V(43)*V(71)
  IRR(155) = RCT(155)*V(43)*V(75)
  IRR(156) = RCT(156)*V(37)*V(71)
  IRR(157) = RCT(157)*V(37)*V(75)
  IRR(158) = RCT(158)*V(40)*V(75)
  IRR(159) = RCT(159)*V(41)*V(71)
  IRR(160) = RCT(160)*V(41)
  IRR(161) = RCT(161)*V(41)*V(75)
  IRR(162) = RCT(162)*V(55)*V(71)
  IRR(163) = RCT(163)*V(55)*V(67)
  IRR(164) = RCT(164)*V(55)*V(75)
  IRR(165) = RCT(165)*V(55)*V(63)
  IRR(166) = RCT(166)*V(55)
  IRR(167) = RCT(167)*V(59)*V(71)
  IRR(168) = RCT(168)*V(59)*V(67)
  IRR(169) = RCT(169)*V(59)*V(63)
  IRR(170) = RCT(170)*V(59)
  IRR(171) = RCT(171)*V(57)*V(71)
  IRR(172) = RCT(172)*V(57)*V(67)
  IRR(173) = RCT(173)*V(57)*V(75)
  IRR(174) = RCT(174)*V(57)
  IRR(175) = RCT(175)*V(66)*V(71)
  IRR(176) = RCT(176)*V(66)
  IRR(177) = RCT(177)*V(62)*V(71)
  IRR(178) = RCT(178)*V(62)
  IRR(179) = RCT(179)*V(38)*V(71)
  IRR(180) = RCT(180)*V(38)*V(67)
  IRR(181) = RCT(181)*V(35)*V(71)
  IRR(182) = RCT(182)*V(35)
  IRR(183) = RCT(183)*V(36)*V(71)
  IRR(184) = RCT(184)*V(36)
  IRR(185) = RCT(185)*V(12)*V(71)
  IRR(186) = RCT(186)*V(48)*V(71)
  IRR(187) = RCT(187)*V(48)*V(67)
  IRR(188) = RCT(188)*V(48)*V(75)
  IRR(189) = RCT(189)*V(48)*V(63)
  IRR(190) = RCT(190)*V(52)*V(71)
  IRR(191) = RCT(191)*V(52)*V(67)
  IRR(192) = RCT(192)*V(52)*V(75)
  IRR(193) = RCT(193)*V(52)*V(63)
  IRR(194) = RCT(194)*V(54)*V(71)
  IRR(195) = RCT(195)*V(54)*V(67)
  IRR(196) = RCT(196)*V(54)*V(75)
  IRR(197) = RCT(197)*V(54)*V(63)
  IRR(198) = RCT(198)*V(13)*V(71)
  IRR(199) = RCT(199)*V(20)*V(71)
  IRR(200) = RCT(200)*V(39)*V(71)
  IRR(201) = RCT(201)*V(24)*V(71)
  IRR(202) = RCT(202)*V(33)*V(71)
  IRR(203) = RCT(203)*V(26)*V(71)
  IRR(204) = RCT(204)*V(34)*V(71)
  IRR(205) = RCT(205)*V(27)*V(71)
  IRR(206) = RCT(206)*V(56)*V(71)
  IRR(207) = RCT(207)*V(56)*V(67)
  IRR(208) = RCT(208)*V(56)*V(75)
  IRR(209) = RCT(209)*V(56)*V(63)
  IRR(210) = RCT(210)*V(58)*V(71)
  IRR(211) = RCT(211)*V(58)*V(67)
  IRR(212) = RCT(212)*V(58)*V(75)
  IRR(213) = RCT(213)*V(58)*V(63)
  IRR(214) = RCT(214)*V(39)*V(67)
  IRR(215) = RCT(215)*V(50)*V(71)
  IRR(216) = RCT(216)*V(50)*V(67)
  IRR(217) = RCT(217)*V(50)*V(75)
  IRR(218) = RCT(218)*V(50)*V(63)
  IRR(219) = RCT(219)*V(10)
  IRR(220) = RCT(220)*V(78)
  IRR(221) = RCT(221)*V(10)
  IRR(222) = RCT(222)*V(1)
  IRR(223) = RCT(223)*V(47)
  IRR(224) = RCT(224)*V(19)
  IRR(225) = RCT(225)*V(2)
      
END SUBROUTINE saprc99_IRRFun
















SUBROUTINE saprc99_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(393)


  B(1) = RCT(1)

  B(2) = RCT(2)*F(2)

  B(4) = RCT(3)*V(67)

  B(5) = RCT(3)*V(63)

  B(6) = RCT(4)*V(73)*F(2)

  B(7) = RCT(4)*V(63)*F(2)

  B(9) = RCT(5)*V(74)

  B(10) = RCT(5)*V(63)

  B(11) = RCT(6)*V(74)

  B(12) = RCT(6)*V(63)

  B(13) = RCT(7)*V(73)

  B(14) = RCT(7)*V(67)

  B(15) = RCT(8)*V(74)

  B(16) = RCT(8)*V(67)

  B(17) = RCT(9)*V(75)

  B(18) = RCT(9)*V(73)

  B(19) = RCT(10)*2*V(73)*F(2)

  B(21) = RCT(11)*V(75)

  B(22) = RCT(11)*V(74)

  B(23) = RCT(12)

  B(24) = RCT(13)*F(1)

  B(26) = RCT(14)*V(75)

  B(27) = RCT(14)*V(74)

  B(28) = RCT(15)

  B(29) = RCT(16)

  B(30) = RCT(17)

  B(31) = RCT(18)

  B(32) = RCT(19)*F(1)

  B(34) = RCT(20)*F(2)

  B(36) = RCT(21)*V(73)

  B(37) = RCT(21)*V(71)

  B(38) = RCT(22)

  B(39) = RCT(23)

  B(40) = RCT(24)*V(71)

  B(41) = RCT(24)*V(23)

  B(42) = RCT(25)*V(74)

  B(43) = RCT(25)*V(71)

  B(44) = RCT(26)*V(75)

  B(45) = RCT(26)*V(71)

  B(46) = RCT(27)*V(71)

  B(47) = RCT(27)*V(47)

  B(48) = RCT(28)

  B(49) = RCT(29)*V(71)

  B(50) = RCT(29)*V(45)

  B(51) = RCT(30)*V(71)

  B(52) = RCT(30)*V(67)

  B(53) = RCT(31)*V(78)

  B(54) = RCT(31)*V(73)

  B(55) = RCT(32)*V(78)

  B(56) = RCT(32)*V(74)

  B(57) = RCT(33)

  B(58) = RCT(34)

  B(59) = RCT(35)*V(71)

  B(60) = RCT(35)*V(32)

  B(61) = RCT(36)*V(78)

  B(62) = RCT(36)*V(67)

  B(63) = RCT(37)*2*V(78)

  B(64) = RCT(38)*2*V(78)*F(1)

  B(66) = RCT(39)*V(78)

  B(67) = RCT(39)*V(75)

  B(68) = RCT(40)*2*V(75)

  B(69) = RCT(41)

  B(70) = RCT(42)*V(71)

  B(71) = RCT(42)*V(19)

  B(72) = RCT(43)*V(78)

  B(73) = RCT(43)*V(71)

  B(74) = RCT(44)*V(71)

  B(75) = RCT(44)*V(10)

  B(76) = RCT(45)*F(2)

  B(78) = RCT(46)*V(76)

  B(79) = RCT(46)*V(73)

  B(80) = RCT(47)*V(78)

  B(81) = RCT(47)*V(76)

  B(82) = RCT(48)*V(76)

  B(83) = RCT(48)*V(75)

  B(84) = RCT(49)*2*V(76)

  B(85) = RCT(50)*2*V(76)

  B(86) = RCT(51)*V(77)

  B(87) = RCT(51)*V(73)

  B(88) = RCT(52)*V(78)

  B(89) = RCT(52)*V(77)

  B(90) = RCT(53)*V(77)

  B(91) = RCT(53)*V(75)

  B(92) = RCT(54)*V(77)

  B(93) = RCT(54)*V(76)

  B(94) = RCT(55)*2*V(77)

  B(95) = RCT(56)*V(73)

  B(96) = RCT(56)*V(53)

  B(97) = RCT(57)*V(78)

  B(98) = RCT(57)*V(53)

  B(99) = RCT(58)*V(75)

  B(100) = RCT(58)*V(53)

  B(101) = RCT(59)*V(76)

  B(102) = RCT(59)*V(53)

  B(103) = RCT(60)*V(77)

  B(104) = RCT(60)*V(53)

  B(105) = RCT(61)*2*V(53)

  B(106) = RCT(62)*V(73)

  B(107) = RCT(62)*V(69)

  B(108) = RCT(63)*V(78)

  B(109) = RCT(63)*V(69)

  B(110) = RCT(64)*V(76)

  B(111) = RCT(64)*V(69)

  B(112) = RCT(65)*V(75)

  B(113) = RCT(65)*V(69)

  B(114) = RCT(66)*V(77)

  B(115) = RCT(66)*V(69)

  B(116) = RCT(67)*V(69)

  B(117) = RCT(67)*V(53)

  B(118) = RCT(68)*2*V(69)

  B(119) = RCT(69)*V(74)

  B(120) = RCT(69)*V(68)

  B(121) = RCT(70)

  B(122) = RCT(71)*V(73)

  B(123) = RCT(71)*V(68)

  B(124) = RCT(72)*V(78)

  B(125) = RCT(72)*V(68)

  B(126) = RCT(73)*V(75)

  B(127) = RCT(73)*V(68)

  B(128) = RCT(74)*V(76)

  B(129) = RCT(74)*V(68)

  B(130) = RCT(75)*V(77)

  B(131) = RCT(75)*V(68)

  B(132) = RCT(76)*V(68)

  B(133) = RCT(76)*V(53)

  B(134) = RCT(77)*V(69)

  B(135) = RCT(77)*V(68)

  B(136) = RCT(78)*2*V(68)

  B(137) = RCT(79)*V(74)

  B(138) = RCT(79)*V(72)

  B(139) = RCT(80)

  B(140) = RCT(81)*V(73)

  B(141) = RCT(81)*V(72)

  B(142) = RCT(82)*V(78)

  B(143) = RCT(82)*V(72)

  B(144) = RCT(83)*V(75)

  B(145) = RCT(83)*V(72)

  B(146) = RCT(84)*V(76)

  B(147) = RCT(84)*V(72)

  B(148) = RCT(85)*V(77)

  B(149) = RCT(85)*V(72)

  B(150) = RCT(86)*V(72)

  B(151) = RCT(86)*V(53)

  B(152) = RCT(87)*V(72)

  B(153) = RCT(87)*V(69)

  B(154) = RCT(88)*V(72)

  B(155) = RCT(88)*V(68)

  B(156) = RCT(89)*2*V(72)

  B(157) = RCT(90)*V(79)

  B(158) = RCT(90)*V(74)

  B(159) = RCT(91)

  B(160) = RCT(92)*V(79)

  B(161) = RCT(92)*V(73)

  B(162) = RCT(93)*V(79)

  B(163) = RCT(93)*V(78)

  B(164) = RCT(94)*V(79)

  B(165) = RCT(94)*V(75)

  B(166) = RCT(95)*V(79)

  B(167) = RCT(95)*V(76)

  B(168) = RCT(96)*V(79)

  B(169) = RCT(96)*V(77)

  B(170) = RCT(97)*V(79)

  B(171) = RCT(97)*V(53)

  B(172) = RCT(98)*V(79)

  B(173) = RCT(98)*V(69)

  B(174) = RCT(99)*V(79)

  B(175) = RCT(99)*V(68)

  B(176) = RCT(100)*V(79)

  B(177) = RCT(100)*V(72)

  B(178) = RCT(101)*2*V(79)

  B(179) = RCT(102)*V(74)

  B(180) = RCT(102)*V(70)

  B(181) = RCT(103)

  B(182) = RCT(104)*V(73)

  B(183) = RCT(104)*V(70)

  B(184) = RCT(105)*V(78)

  B(185) = RCT(105)*V(70)

  B(186) = RCT(106)*V(75)

  B(187) = RCT(106)*V(70)

  B(188) = RCT(107)*V(76)

  B(189) = RCT(107)*V(70)

  B(190) = RCT(108)*V(77)

  B(191) = RCT(108)*V(70)

  B(192) = RCT(109)*V(70)

  B(193) = RCT(109)*V(53)

  B(194) = RCT(110)*V(70)

  B(195) = RCT(110)*V(69)

  B(196) = RCT(111)*V(70)

  B(197) = RCT(111)*V(68)

  B(198) = RCT(112)*V(72)

  B(199) = RCT(112)*V(70)

  B(200) = RCT(113)*V(79)

  B(201) = RCT(113)*V(70)

  B(202) = RCT(114)*2*V(70)

  B(203) = RCT(115)*V(74)

  B(204) = RCT(115)*V(25)

  B(205) = RCT(116)

  B(206) = RCT(117)*V(74)

  B(207) = RCT(117)*V(51)

  B(208) = RCT(118)*V(78)

  B(209) = RCT(118)*V(51)

  B(210) = RCT(119)

  B(211) = RCT(120)*V(74)

  B(212) = RCT(120)*V(31)

  B(213) = RCT(121)*V(78)

  B(214) = RCT(121)*V(31)

  B(215) = RCT(122)

  B(216) = RCT(123)

  B(217) = RCT(124)

  B(218) = RCT(125)*V(71)

  B(219) = RCT(125)*V(61)

  B(220) = RCT(126)*V(78)

  B(221) = RCT(126)*V(61)

  B(222) = RCT(127)

  B(223) = RCT(128)*V(73)

  B(224) = RCT(128)*V(30)

  B(225) = RCT(129)*V(75)

  B(226) = RCT(129)*V(61)

  B(227) = RCT(130)*V(71)

  B(228) = RCT(130)*V(60)

  B(229) = RCT(131)

  B(230) = RCT(132)*V(75)

  B(231) = RCT(132)*V(60)

  B(232) = RCT(133)*V(71)

  B(233) = RCT(133)*V(64)

  B(234) = RCT(134)

  B(235) = RCT(135)*V(75)

  B(236) = RCT(135)*V(64)

  B(237) = RCT(136)*V(71)

  B(238) = RCT(136)*V(46)

  B(239) = RCT(137)

  B(240) = RCT(138)*V(71)

  B(241) = RCT(138)*V(65)

  B(242) = RCT(139)

  B(243) = RCT(140)*V(71)

  B(244) = RCT(140)*V(28)

  B(245) = RCT(141)*V(71)

  B(246) = RCT(141)*V(21)

  B(247) = RCT(142)*V(71)

  B(248) = RCT(142)*V(29)

  B(249) = RCT(143)

  B(250) = RCT(144)*V(71)

  B(251) = RCT(144)*V(42)

  B(252) = RCT(145)

  B(253) = RCT(146)

  B(254) = RCT(147)

  B(255) = RCT(148)*V(71)

  B(256) = RCT(148)*V(49)

  B(257) = RCT(149)*V(75)

  B(258) = RCT(149)*V(49)

  B(259) = RCT(150)

  B(260) = RCT(151)*V(71)

  B(261) = RCT(151)*V(44)

  B(262) = RCT(152)*V(75)

  B(263) = RCT(152)*V(44)

  B(264) = RCT(153)

  B(265) = RCT(154)*V(71)

  B(266) = RCT(154)*V(43)

  B(267) = RCT(155)*V(75)

  B(268) = RCT(155)*V(43)

  B(269) = RCT(156)*V(71)

  B(270) = RCT(156)*V(37)

  B(271) = RCT(157)*V(75)

  B(272) = RCT(157)*V(37)

  B(273) = RCT(158)*V(75)

  B(274) = RCT(158)*V(40)

  B(275) = RCT(159)*V(71)

  B(276) = RCT(159)*V(41)

  B(277) = RCT(160)

  B(278) = RCT(161)*V(75)

  B(279) = RCT(161)*V(41)

  B(280) = RCT(162)*V(71)

  B(281) = RCT(162)*V(55)

  B(282) = RCT(163)*V(67)

  B(283) = RCT(163)*V(55)

  B(284) = RCT(164)*V(75)

  B(285) = RCT(164)*V(55)

  B(286) = RCT(165)*V(63)

  B(287) = RCT(165)*V(55)

  B(288) = RCT(166)

  B(289) = RCT(167)*V(71)

  B(290) = RCT(167)*V(59)

  B(291) = RCT(168)*V(67)

  B(292) = RCT(168)*V(59)

  B(293) = RCT(169)*V(63)

  B(294) = RCT(169)*V(59)

  B(295) = RCT(170)

  B(296) = RCT(171)*V(71)

  B(297) = RCT(171)*V(57)

  B(298) = RCT(172)*V(67)

  B(299) = RCT(172)*V(57)

  B(300) = RCT(173)*V(75)

  B(301) = RCT(173)*V(57)

  B(302) = RCT(174)

  B(303) = RCT(175)*V(71)

  B(304) = RCT(175)*V(66)

  B(305) = RCT(176)

  B(306) = RCT(177)*V(71)

  B(307) = RCT(177)*V(62)

  B(308) = RCT(178)

  B(309) = RCT(179)*V(71)

  B(310) = RCT(179)*V(38)

  B(311) = RCT(180)*V(67)

  B(312) = RCT(180)*V(38)

  B(313) = RCT(181)*V(71)

  B(314) = RCT(181)*V(35)

  B(315) = RCT(182)

  B(316) = RCT(183)*V(71)

  B(317) = RCT(183)*V(36)

  B(318) = RCT(184)

  B(319) = RCT(185)*V(71)

  B(320) = RCT(185)*V(12)

  B(321) = RCT(186)*V(71)

  B(322) = RCT(186)*V(48)

  B(323) = RCT(187)*V(67)

  B(324) = RCT(187)*V(48)

  B(325) = RCT(188)*V(75)

  B(326) = RCT(188)*V(48)

  B(327) = RCT(189)*V(63)

  B(328) = RCT(189)*V(48)

  B(329) = RCT(190)*V(71)

  B(330) = RCT(190)*V(52)

  B(331) = RCT(191)*V(67)

  B(332) = RCT(191)*V(52)

  B(333) = RCT(192)*V(75)

  B(334) = RCT(192)*V(52)

  B(335) = RCT(193)*V(63)

  B(336) = RCT(193)*V(52)

  B(337) = RCT(194)*V(71)

  B(338) = RCT(194)*V(54)

  B(339) = RCT(195)*V(67)

  B(340) = RCT(195)*V(54)

  B(341) = RCT(196)*V(75)

  B(342) = RCT(196)*V(54)

  B(343) = RCT(197)*V(63)

  B(344) = RCT(197)*V(54)

  B(345) = RCT(198)*V(71)

  B(346) = RCT(198)*V(13)

  B(347) = RCT(199)*V(71)

  B(348) = RCT(199)*V(20)

  B(349) = RCT(200)*V(71)

  B(350) = RCT(200)*V(39)

  B(351) = RCT(201)*V(71)

  B(352) = RCT(201)*V(24)

  B(353) = RCT(202)*V(71)

  B(354) = RCT(202)*V(33)

  B(355) = RCT(203)*V(71)

  B(356) = RCT(203)*V(26)

  B(357) = RCT(204)*V(71)

  B(358) = RCT(204)*V(34)

  B(359) = RCT(205)*V(71)

  B(360) = RCT(205)*V(27)

  B(361) = RCT(206)*V(71)

  B(362) = RCT(206)*V(56)

  B(363) = RCT(207)*V(67)

  B(364) = RCT(207)*V(56)

  B(365) = RCT(208)*V(75)

  B(366) = RCT(208)*V(56)

  B(367) = RCT(209)*V(63)

  B(368) = RCT(209)*V(56)

  B(369) = RCT(210)*V(71)

  B(370) = RCT(210)*V(58)

  B(371) = RCT(211)*V(67)

  B(372) = RCT(211)*V(58)

  B(373) = RCT(212)*V(75)

  B(374) = RCT(212)*V(58)

  B(375) = RCT(213)*V(63)

  B(376) = RCT(213)*V(58)

  B(377) = RCT(214)*V(67)

  B(378) = RCT(214)*V(39)

  B(379) = RCT(215)*V(71)

  B(380) = RCT(215)*V(50)

  B(381) = RCT(216)*V(67)

  B(382) = RCT(216)*V(50)

  B(383) = RCT(217)*V(75)

  B(384) = RCT(217)*V(50)

  B(385) = RCT(218)*V(63)

  B(386) = RCT(218)*V(50)

  B(387) = RCT(219)

  B(388) = RCT(220)

  B(389) = RCT(221)

  B(390) = RCT(222)

  B(391) = RCT(223)

  B(392) = RCT(224)

  B(393) = RCT(225)



  JVS(1) = -B(390)

  JVS(2) = B(74)+B(387)

  JVS(3) = B(75)

  JVS(4) = -B(393)

  JVS(5) = 0.5*B(377)

  JVS(6) = 0.135*B(381)

  JVS(7) = 0.5*B(378)+0.135*B(382)

  JVS(8) = 0

  JVS(9) = B(223)

  JVS(10) = 0.297*B(349)

  JVS(11) = 0.37*B(323)

  JVS(12) = 0.185*B(381)

  JVS(13) = 0.204*B(331)

  JVS(14) = 0.103*B(339)

  JVS(15) = 0.333*B(282)

  JVS(16) = 0.185*B(363)

  JVS(17) = 0.1*B(298)

  JVS(18) = 0.073*B(371)

  JVS(19) = 0.351*B(291)

  JVS(20) = 0.333*B(283)+0.351*B(292)+0.1*B(299)+0.37*B(324)+0.204*B(332)+0.103*B(340)+0.185*B(364)+0.073*B(372)+0.185&
              &*B(382)

  JVS(21) = 0.297*B(350)

  JVS(22) = B(224)

  JVS(23) = 0

  JVS(24) = 0.17*B(381)

  JVS(25) = 0.05*B(363)

  JVS(26) = 0.129*B(371)

  JVS(27) = 0.05*B(364)+0.129*B(372)+0.17*B(382)

  JVS(28) = 0.25*B(124)+B(128)+B(130)+B(134)

  JVS(29) = B(135)

  JVS(30) = B(129)

  JVS(31) = B(131)

  JVS(32) = 0.25*B(125)

  JVS(33) = 0

  JVS(34) = 0.15*B(331)

  JVS(35) = 0.189*B(339)

  JVS(36) = 0.119*B(363)

  JVS(37) = 0.372*B(298)

  JVS(38) = 0.247*B(371)

  JVS(39) = 0.372*B(299)+0.15*B(332)+0.189*B(340)+0.119*B(364)+0.247*B(372)

  JVS(40) = B(152)+B(172)+2*B(194)

  JVS(41) = 0.25*B(184)+B(188)+B(190)+2*B(195)

  JVS(42) = 0.25*B(142)+B(146)+B(148)+B(153)

  JVS(43) = B(147)+B(166)+B(189)

  JVS(44) = B(149)+B(168)+B(191)

  JVS(45) = 0.25*B(143)+0.25*B(162)+0.25*B(185)

  JVS(46) = 0.25*B(163)+B(167)+B(169)+B(173)

  JVS(47) = 0

  JVS(48) = 0.75*B(124)

  JVS(49) = 0.75*B(125)

  JVS(50) = 0

  JVS(51) = 0.75*B(184)

  JVS(52) = 0.75*B(142)

  JVS(53) = 0.75*B(143)+0.75*B(162)+0.75*B(185)

  JVS(54) = 0.75*B(163)

  JVS(55) = 0

  JVS(56) = 2*B(211)

  JVS(57) = B(383)

  JVS(58) = 2*B(212)

  JVS(59) = B(384)

  JVS(60) = 0

  JVS(61) = 6*B(211)

  JVS(62) = 7*B(277)

  JVS(63) = 0.048*B(379)+0.07*B(381)+2.693*B(383)+0.55*B(385)

  JVS(64) = 0.55*B(386)

  JVS(65) = 0.07*B(382)

  JVS(66) = 0.048*B(380)

  JVS(67) = 6*B(212)

  JVS(68) = 2.693*B(384)

  JVS(69) = -B(74)-B(387)-B(389)

  JVS(70) = -B(75)

  JVS(71) = -B(32)-B(34)

  JVS(72) = B(31)

  JVS(73) = -B(319)

  JVS(74) = -B(320)

  JVS(75) = -B(345)

  JVS(76) = -B(346)

  JVS(77) = -B(264)

  JVS(78) = 0.087*B(359)

  JVS(79) = 0.031*B(339)

  JVS(80) = 0.031*B(340)

  JVS(81) = 0.087*B(360)

  JVS(82) = -B(121)

  JVS(83) = B(119)

  JVS(84) = B(120)

  JVS(85) = -B(139)

  JVS(86) = B(137)

  JVS(87) = B(138)

  JVS(88) = -B(159)

  JVS(89) = B(157)

  JVS(90) = B(158)

  JVS(91) = -B(181)

  JVS(92) = B(179)

  JVS(93) = B(180)

  JVS(94) = -B(69)-B(70)-B(392)

  JVS(95) = -B(71)

  JVS(96) = B(63)+B(64)

  JVS(97) = -B(347)

  JVS(98) = -B(348)

  JVS(99) = -B(245)

  JVS(100) = -B(246)

  JVS(101) = -B(23)-B(24)

  JVS(102) = B(21)

  JVS(103) = B(22)

  JVS(104) = -B(38)-B(39)-B(40)

  JVS(105) = B(36)-B(41)

  JVS(106) = B(37)

  JVS(107) = -B(351)

  JVS(108) = -B(352)

  JVS(109) = 0.236*B(351)

  JVS(110) = -B(203)-B(205)

  JVS(111) = 0.236*B(352)

  JVS(112) = -B(204)

  JVS(113) = -B(355)

  JVS(114) = -B(356)

  JVS(115) = -B(359)

  JVS(116) = -B(360)

  JVS(117) = -B(243)

  JVS(118) = 0.25*B(110)

  JVS(119) = -B(244)

  JVS(120) = B(84)+0.25*B(92)+0.25*B(111)

  JVS(121) = 0.25*B(93)

  JVS(122) = -B(247)-B(249)

  JVS(123) = -B(248)

  JVS(124) = B(80)

  JVS(125) = B(81)

  JVS(126) = -B(222)-B(223)

  JVS(127) = B(220)

  JVS(128) = -B(224)

  JVS(129) = B(221)

  JVS(130) = -B(211)-B(213)-B(215)

  JVS(131) = B(273)

  JVS(132) = -B(212)

  JVS(133) = B(274)

  JVS(134) = -B(214)

  JVS(135) = -B(57)-B(58)-B(59)

  JVS(136) = -B(60)

  JVS(137) = B(55)

  JVS(138) = B(56)

  JVS(139) = -B(353)

  JVS(140) = -B(354)

  JVS(141) = -B(357)

  JVS(142) = -B(358)

  JVS(143) = 0.099*B(359)

  JVS(144) = 0.108*B(357)

  JVS(145) = -B(313)-B(315)

  JVS(146) = -B(314)+0.108*B(358)+0.099*B(360)

  JVS(147) = 0.093*B(359)

  JVS(148) = 0.051*B(357)

  JVS(149) = -B(316)-B(318)

  JVS(150) = -B(317)+0.051*B(358)+0.093*B(360)

  JVS(151) = 0.187*B(359)

  JVS(152) = 0.207*B(357)

  JVS(153) = -B(269)-B(271)

  JVS(154) = -B(270)+0.207*B(358)+0.187*B(360)

  JVS(155) = -B(272)

  JVS(156) = 0.561*B(359)

  JVS(157) = 0.491*B(357)

  JVS(158) = -B(309)-B(311)

  JVS(159) = -B(312)

  JVS(160) = -B(310)+0.491*B(358)+0.561*B(360)

  JVS(161) = -B(349)-B(377)

  JVS(162) = -B(378)

  JVS(163) = -B(350)

  JVS(164) = B(213)+B(215)

  JVS(165) = -B(273)

  JVS(166) = B(206)

  JVS(167) = B(207)

  JVS(168) = -B(274)

  JVS(169) = B(214)

  JVS(170) = 0.05*B(359)

  JVS(171) = 0.059*B(357)

  JVS(172) = -B(275)-B(277)-B(278)

  JVS(173) = 0.061*B(369)+0.042*B(371)+0.015*B(373)

  JVS(174) = 0.042*B(372)

  JVS(175) = -B(276)+0.059*B(358)+0.05*B(360)+0.061*B(370)

  JVS(176) = -B(279)+0.015*B(374)

  JVS(177) = -B(250)-B(252)

  JVS(178) = B(108)

  JVS(179) = -B(251)

  JVS(180) = B(88)

  JVS(181) = B(89)+B(109)

  JVS(182) = 0.017*B(357)

  JVS(183) = -B(265)-B(267)

  JVS(184) = B(208)+B(210)

  JVS(185) = -B(266)+0.017*B(358)

  JVS(186) = -B(268)

  JVS(187) = B(209)

  JVS(188) = 0.287*B(359)

  JVS(189) = 0.119*B(357)

  JVS(190) = 0.5*B(315)

  JVS(191) = 0.5*B(318)

  JVS(192) = 0.23*B(269)

  JVS(193) = -B(259)-B(260)-B(262)

  JVS(194) = 0.084*B(280)+0.9*B(282)

  JVS(195) = 0.174*B(296)+0.742*B(298)+0.008*B(300)

  JVS(196) = 0.3*B(289)+0.95*B(291)

  JVS(197) = 0.9*B(283)+0.95*B(292)+0.742*B(299)

  JVS(198) = -B(261)+0.23*B(270)+0.084*B(281)+0.3*B(290)+0.174*B(297)+0.119*B(358)+0.287*B(360)

  JVS(199) = -B(263)+0.008*B(301)

  JVS(200) = 0.002*B(353)

  JVS(201) = B(315)

  JVS(202) = B(318)

  JVS(203) = B(309)+1.5*B(311)

  JVS(204) = 0.393*B(349)+1.5*B(377)

  JVS(205) = B(259)+B(260)+B(262)

  JVS(206) = -B(49)

  JVS(207) = 0.5*B(323)+0.491*B(327)

  JVS(208) = 2*B(253)+B(254)+1.26*B(255)+1.26*B(257)

  JVS(209) = 0.51*B(381)

  JVS(210) = 0.275*B(331)

  JVS(211) = 0.157*B(339)

  JVS(212) = 0.416*B(280)+0.45*B(282)+0.5*B(284)+0.67*B(288)

  JVS(213) = 0.345*B(363)

  JVS(214) = 0.336*B(296)+0.498*B(298)+0.572*B(300)+1.233*B(302)

  JVS(215) = 0.265*B(371)+0.012*B(375)

  JVS(216) = 0.475*B(291)+0.7*B(295)

  JVS(217) = B(229)

  JVS(218) = B(216)+B(217)+B(218)+B(225)

  JVS(219) = 0.491*B(328)+0.012*B(376)

  JVS(220) = 0.034*B(232)+B(234)

  JVS(221) = 0.45*B(283)+0.475*B(292)+0.498*B(299)+1.5*B(312)+0.5*B(324)+0.275*B(332)+0.157*B(340)+0.345*B(364)+0.265&
               &*B(372)+1.5*B(378)+0.51*B(382)

  JVS(222) = -B(50)+B(219)+0.034*B(233)+1.26*B(256)+B(261)+0.416*B(281)+0.336*B(297)+B(310)+0.393*B(350)+0.002*B(354)

  JVS(223) = B(226)+1.26*B(258)+B(263)+0.5*B(285)+0.572*B(301)

  JVS(224) = 0.704*B(347)

  JVS(225) = 0.024*B(351)

  JVS(226) = B(205)

  JVS(227) = 0.072*B(355)

  JVS(228) = 0.452*B(353)

  JVS(229) = -B(237)-B(239)

  JVS(230) = 0.13*B(339)

  JVS(231) = 0.005*B(361)+0.001*B(363)+0.024*B(365)

  JVS(232) = 0.127*B(369)+0.045*B(371)+0.102*B(373)

  JVS(233) = 0.006*B(306)+0.02*B(308)

  JVS(234) = 0.13*B(340)+0.001*B(364)+0.045*B(372)

  JVS(235) = -B(238)+0.006*B(307)+0.704*B(348)+0.024*B(352)+0.452*B(354)+0.072*B(356)+0.005*B(362)+0.127*B(370)

  JVS(236) = 0

  JVS(237) = 0.024*B(366)+0.102*B(374)

  JVS(238) = 2*B(24)

  JVS(239) = B(271)

  JVS(240) = B(273)

  JVS(241) = B(278)

  JVS(242) = B(267)

  JVS(243) = B(262)

  JVS(244) = -B(46)-B(48)-B(391)

  JVS(245) = B(257)

  JVS(246) = 0

  JVS(247) = 0.5*B(284)

  JVS(248) = 0.15*B(300)

  JVS(249) = 0

  JVS(250) = 0

  JVS(251) = B(230)

  JVS(252) = B(225)

  JVS(253) = B(235)

  JVS(254) = 0

  JVS(255) = B(42)-B(47)

  JVS(256) = B(43)

  JVS(257) = 0.2*B(66)+B(226)+B(231)+B(236)+B(258)+B(263)+B(268)+B(272)+B(274)+B(279)+0.5*B(285)+0.15*B(301)

  JVS(258) = 0.2*B(67)

  JVS(259) = -B(321)-B(323)-B(325)-B(327)

  JVS(260) = -B(328)

  JVS(261) = -B(324)

  JVS(262) = -B(322)

  JVS(263) = -B(326)

  JVS(264) = 0.097*B(359)

  JVS(265) = 0.118*B(357)

  JVS(266) = 0.5*B(315)

  JVS(267) = 0.5*B(318)

  JVS(268) = B(311)

  JVS(269) = 0.607*B(349)

  JVS(270) = 0.23*B(265)

  JVS(271) = 0.009*B(327)

  JVS(272) = -B(253)-B(254)-B(255)-B(257)

  JVS(273) = 0

  JVS(274) = 0.001*B(339)

  JVS(275) = 0.15*B(296)+0.023*B(298)

  JVS(276) = 0.009*B(328)

  JVS(277) = 0.023*B(299)+B(312)+0.001*B(340)

  JVS(278) = -B(256)+0.23*B(266)+0.15*B(297)+0.607*B(350)+0.118*B(358)+0.097*B(360)

  JVS(279) = -B(258)

  JVS(280) = 0

  JVS(281) = -B(379)-B(381)-B(383)-B(385)

  JVS(282) = -B(386)

  JVS(283) = -B(382)

  JVS(284) = -B(380)

  JVS(285) = -B(384)

  JVS(286) = 0.24*B(269)+B(271)

  JVS(287) = 0.24*B(265)+B(267)

  JVS(288) = -B(206)-B(208)-B(210)

  JVS(289) = B(174)

  JVS(290) = B(200)

  JVS(291) = 0.24*B(266)+0.24*B(270)

  JVS(292) = B(176)

  JVS(293) = B(160)

  JVS(294) = -B(207)

  JVS(295) = B(164)+B(268)+B(272)

  JVS(296) = -B(209)

  JVS(297) = B(161)+B(165)+B(175)+B(177)+2*B(178)+B(201)

  JVS(298) = -B(329)-B(331)-B(333)-B(335)

  JVS(299) = -B(336)

  JVS(300) = -B(332)

  JVS(301) = -B(330)

  JVS(302) = -B(334)

  JVS(303) = 0.559*B(351)

  JVS(304) = 0.948*B(355)

  JVS(305) = 0.936*B(353)

  JVS(306) = B(313)+B(315)

  JVS(307) = B(316)+B(318)

  JVS(308) = B(237)

  JVS(309) = 0.079*B(329)+0.126*B(331)+0.187*B(333)+0.24*B(335)

  JVS(310) = -B(95)-B(97)-B(99)-B(101)-B(103)-B(116)-B(132)-B(150)-B(170)-B(192)

  JVS(311) = 0.5*B(337)+0.729*B(339)+0.75*B(341)

  JVS(312) = 0.205*B(361)+0.488*B(365)

  JVS(313) = 0.001*B(369)+0.137*B(371)+0.711*B(373)

  JVS(314) = 0.675*B(289)

  JVS(315) = 0.596*B(306)+0.152*B(308)

  JVS(316) = 0.24*B(336)

  JVS(317) = 0.616*B(240)

  JVS(318) = 0.515*B(305)

  JVS(319) = 0.126*B(332)+0.729*B(340)+0.137*B(372)

  JVS(320) = -B(133)+B(174)

  JVS(321) = -B(117)

  JVS(322) = -B(193)+B(200)

  JVS(323) = B(238)+0.616*B(241)+0.675*B(290)+0.596*B(307)+B(314)+B(317)+0.079*B(330)+0.5*B(338)+0.559*B(352)+0.936&
               &*B(354)+0.948*B(356)+0.205*B(362)+0.001*B(370)

  JVS(324) = -B(151)+B(176)

  JVS(325) = -B(96)+B(160)

  JVS(326) = 0

  JVS(327) = -B(100)+B(164)+0.187*B(334)+0.75*B(342)+0.488*B(366)+0.711*B(374)

  JVS(328) = -B(102)

  JVS(329) = -B(104)

  JVS(330) = -B(98)

  JVS(331) = B(161)+B(165)-B(171)+B(175)+B(177)+2*B(178)+B(201)

  JVS(332) = -B(337)-B(339)-B(341)-B(343)

  JVS(333) = -B(344)

  JVS(334) = -B(340)

  JVS(335) = -B(338)

  JVS(336) = -B(342)

  JVS(337) = 0.23*B(329)+0.39*B(331)

  JVS(338) = -B(280)-B(282)-B(284)-B(286)-B(288)

  JVS(339) = 0.025*B(369)+0.026*B(371)+0.012*B(375)

  JVS(340) = -B(287)+0.012*B(376)

  JVS(341) = -B(283)+0.39*B(332)+0.026*B(372)

  JVS(342) = -B(281)+0.23*B(330)+0.025*B(370)

  JVS(343) = -B(285)

  JVS(344) = -B(361)-B(363)-B(365)-B(367)

  JVS(345) = -B(368)

  JVS(346) = -B(364)

  JVS(347) = -B(362)

  JVS(348) = -B(366)

  JVS(349) = 0.357*B(329)+0.936*B(333)

  JVS(350) = -B(296)-B(298)-B(300)-B(302)

  JVS(351) = 0.025*B(369)

  JVS(352) = 0

  JVS(353) = -B(299)

  JVS(354) = -B(297)+0.357*B(330)+0.025*B(370)

  JVS(355) = -B(301)+0.936*B(334)

  JVS(356) = -B(369)-B(371)-B(373)-B(375)

  JVS(357) = -B(376)

  JVS(358) = -B(372)

  JVS(359) = -B(370)

  JVS(360) = -B(374)

  JVS(361) = 0.32*B(329)+0.16*B(331)

  JVS(362) = 0.019*B(371)+0.048*B(373)

  JVS(363) = -B(289)-B(291)-B(293)-B(295)

  JVS(364) = -B(294)

  JVS(365) = -B(292)+0.16*B(332)+0.019*B(372)

  JVS(366) = -B(290)+0.32*B(330)

  JVS(367) = 0.048*B(374)

  JVS(368) = B(345)

  JVS(369) = 0.96*B(245)

  JVS(370) = 0.445*B(351)

  JVS(371) = 0.099*B(355)

  JVS(372) = 0.455*B(353)

  JVS(373) = 0.195*B(321)+0.25*B(327)

  JVS(374) = 0.984*B(379)+0.5*B(381)

  JVS(375) = 0.294*B(361)+0.154*B(363)+0.009*B(365)

  JVS(376) = 0.129*B(296)+0.047*B(298)+0.467*B(302)

  JVS(377) = 0.732*B(369)+0.456*B(371)+0.507*B(373)

  JVS(378) = -B(227)-B(229)-B(230)

  JVS(379) = 0.439*B(306)+0.431*B(308)

  JVS(380) = 0.25*B(328)

  JVS(381) = 0.034*B(232)+B(234)

  JVS(382) = 0.482*B(240)+B(242)

  JVS(383) = 0.084*B(303)+0.246*B(305)

  JVS(384) = 0.047*B(299)+0.154*B(364)+0.456*B(372)+0.5*B(382)

  JVS(385) = B(154)

  JVS(386) = B(198)

  JVS(387) = -B(228)+0.034*B(233)+0.482*B(241)+0.96*B(246)+0.129*B(297)+0.084*B(304)+0.439*B(307)+0.195*B(322)+B(346)&
               &+0.445*B(352)+0.455*B(354)+0.099*B(356)+0.294*B(362)+0.732*B(370)+0.984*B(380)

  JVS(388) = B(140)+B(144)+B(155)+2*B(156)+B(176)+B(199)

  JVS(389) = B(141)

  JVS(390) = B(145)-B(231)+0.009*B(366)+0.507*B(374)

  JVS(391) = B(177)

  JVS(392) = 0.081*B(245)

  JVS(393) = 0.026*B(351)

  JVS(394) = 0.026*B(355)

  JVS(395) = B(243)

  JVS(396) = 0.35*B(247)+B(249)

  JVS(397) = B(222)

  JVS(398) = 0.024*B(353)

  JVS(399) = 0.096*B(349)

  JVS(400) = B(237)

  JVS(401) = 1.61*B(321)+B(323)+0.191*B(327)

  JVS(402) = B(254)

  JVS(403) = 0.984*B(379)+0.5*B(381)

  JVS(404) = 0

  JVS(405) = 0.624*B(329)+0.592*B(331)+0.24*B(335)

  JVS(406) = 0.276*B(337)+0.235*B(339)

  JVS(407) = 0.084*B(280)+0.2*B(282)+0.67*B(288)

  JVS(408) = 0.732*B(361)+0.5*B(363)

  JVS(409) = 0.055*B(296)+0.125*B(298)+0.227*B(300)+0.3*B(302)

  JVS(410) = 0.244*B(369)+0.269*B(371)+0.079*B(373)

  JVS(411) = 0.3*B(289)+0.1*B(291)

  JVS(412) = -B(216)-B(217)-B(218)-B(220)-B(225)

  JVS(413) = 0.01*B(306)+0.134*B(308)

  JVS(414) = 0.191*B(328)+0.24*B(336)

  JVS(415) = 0.115*B(240)

  JVS(416) = 0.213*B(303)+0.506*B(305)

  JVS(417) = 0.2*B(283)+0.1*B(292)+0.125*B(299)+B(324)+0.592*B(332)+0.235*B(340)+0.5*B(364)+0.269*B(372)+0.5*B(382)

  JVS(418) = B(128)+B(196)

  JVS(419) = 0.75*B(110)

  JVS(420) = B(182)+B(186)+B(188)+B(197)+B(198)+B(200)+2*B(202)

  JVS(421) = -B(219)+B(238)+0.115*B(241)+B(244)+0.081*B(246)+0.35*B(248)+0.084*B(281)+0.3*B(290)+0.055*B(297)+0.213&
               &*B(304)+0.01*B(307)+1.61*B(322)+0.624*B(330)+0.276*B(338)+0.096*B(350)+0.026*B(352)+0.024*B(354)+0.026&
               &*B(356)+0.732*B(362)+0.244*B(370)+0.984*B(380)

  JVS(422) = B(146)+B(199)

  JVS(423) = B(78)+B(183)

  JVS(424) = 0

  JVS(425) = B(82)+B(187)-B(226)+0.227*B(301)+0.079*B(374)

  JVS(426) = B(79)+B(83)+B(84)+2*B(85)+0.75*B(92)+0.75*B(111)+B(129)+B(147)+B(166)+B(189)

  JVS(427) = 0.75*B(93)

  JVS(428) = -B(221)

  JVS(429) = B(167)+B(201)

  JVS(430) = B(203)

  JVS(431) = 0.276*B(341)

  JVS(432) = 0.511*B(365)

  JVS(433) = 0.572*B(300)

  JVS(434) = 0.321*B(373)

  JVS(435) = -0.69*B(306)-B(308)

  JVS(436) = 0

  JVS(437) = 0

  JVS(438) = B(106)

  JVS(439) = -0.69*B(307)

  JVS(440) = B(107)

  JVS(441) = B(204)

  JVS(442) = 0.572*B(301)+0.276*B(342)+0.511*B(366)+0.321*B(374)

  JVS(443) = B(34)

  JVS(444) = -B(327)

  JVS(445) = -B(385)

  JVS(446) = -B(335)

  JVS(447) = -B(343)

  JVS(448) = -B(286)

  JVS(449) = -B(367)

  JVS(450) = -B(375)

  JVS(451) = -B(293)

  JVS(452) = -B(2)-B(4)-B(6)-B(9)-B(11)-B(287)-B(294)-B(328)-B(336)-B(344)-B(368)-B(376)-B(386)

  JVS(453) = -B(5)+B(30)

  JVS(454) = 0

  JVS(455) = -B(7)

  JVS(456) = B(1)-B(10)-B(12)

  JVS(457) = B(29)

  JVS(458) = 0.261*B(347)

  JVS(459) = 0.122*B(351)

  JVS(460) = 0.204*B(355)

  JVS(461) = 0.244*B(353)

  JVS(462) = B(313)

  JVS(463) = B(316)

  JVS(464) = B(309)

  JVS(465) = B(250)+B(252)

  JVS(466) = B(325)

  JVS(467) = 0.45*B(385)

  JVS(468) = 0.474*B(337)+0.205*B(339)+0.474*B(341)+0.147*B(343)

  JVS(469) = B(286)

  JVS(470) = 0.497*B(361)+0.363*B(363)+0.037*B(365)+0.45*B(367)

  JVS(471) = 0.013*B(296)+0.218*B(300)

  JVS(472) = 0.511*B(369)+0.305*B(371)+0.151*B(373)+0.069*B(375)

  JVS(473) = 0.675*B(289)+0.45*B(293)

  JVS(474) = 0.213*B(306)+0.147*B(308)

  JVS(475) = B(287)+0.45*B(294)+0.147*B(344)+0.45*B(368)+0.069*B(376)+0.45*B(386)

  JVS(476) = -B(232)-B(234)-B(235)

  JVS(477) = 0.37*B(240)

  JVS(478) = 0.558*B(303)+0.71*B(305)

  JVS(479) = 0.205*B(340)+0.363*B(364)+0.305*B(372)

  JVS(480) = 0

  JVS(481) = -B(233)+0.37*B(241)+B(251)+0.675*B(290)+0.013*B(297)+0.558*B(304)+0.213*B(307)+B(310)+B(314)+B(317)+0.474&
               &*B(338)+0.261*B(348)+0.122*B(352)+0.244*B(354)+0.204*B(356)+0.497*B(362)+0.511*B(370)

  JVS(482) = 0

  JVS(483) = 0

  JVS(484) = -B(236)+0.218*B(301)+B(326)+0.474*B(342)+0.037*B(366)+0.151*B(374)

  JVS(485) = 0

  JVS(486) = 0

  JVS(487) = 0.332*B(351)

  JVS(488) = 0.089*B(355)

  JVS(489) = 0.11*B(353)

  JVS(490) = 0.55*B(385)

  JVS(491) = 0.416*B(280)

  JVS(492) = 0.437*B(367)

  JVS(493) = 0.15*B(296)+0.21*B(298)+0.233*B(302)

  JVS(494) = 0.072*B(369)+0.026*B(371)+0.001*B(373)+0.659*B(375)

  JVS(495) = 0.55*B(293)

  JVS(496) = 0.177*B(306)+0.243*B(308)

  JVS(497) = 0.55*B(294)+0.437*B(368)+0.659*B(376)+0.55*B(386)

  JVS(498) = -B(240)-B(242)

  JVS(499) = 0.115*B(303)

  JVS(500) = 0.21*B(299)+0.026*B(372)

  JVS(501) = 0.5*B(110)+B(112)+0.5*B(114)+B(118)

  JVS(502) = -B(241)+0.416*B(281)+0.15*B(297)+0.115*B(304)+0.177*B(307)+0.332*B(352)+0.11*B(354)+0.089*B(356)+0.072&
               &*B(370)

  JVS(503) = 0

  JVS(504) = 0

  JVS(505) = B(113)+0.001*B(374)

  JVS(506) = 0.5*B(111)

  JVS(507) = 0.5*B(115)

  JVS(508) = 0.417*B(355)

  JVS(509) = 0.125*B(353)

  JVS(510) = 0.055*B(357)

  JVS(511) = 0.1*B(331)+0.75*B(335)

  JVS(512) = 0.276*B(337)+0.276*B(339)+0.853*B(343)

  JVS(513) = 0.119*B(361)+0.215*B(363)+0.113*B(367)

  JVS(514) = 0.332*B(296)

  JVS(515) = 0.043*B(371)+0.259*B(375)

  JVS(516) = 0.7*B(295)

  JVS(517) = 0.048*B(306)+0.435*B(308)

  JVS(518) = 0.75*B(336)+0.853*B(344)+0.113*B(368)+0.259*B(376)

  JVS(519) = -0.671*B(303)-B(305)

  JVS(520) = 0.1*B(332)+0.276*B(340)+0.215*B(364)+0.043*B(372)

  JVS(521) = B(134)

  JVS(522) = 0.5*B(110)+0.5*B(114)+B(118)+B(135)+B(152)+B(172)

  JVS(523) = 0.332*B(297)-0.671*B(304)+0.048*B(307)+0.276*B(338)+0.125*B(354)+0.417*B(356)+0.055*B(358)+0.119*B(362)

  JVS(524) = B(153)

  JVS(525) = 0

  JVS(526) = 0

  JVS(527) = 0

  JVS(528) = 0.5*B(111)

  JVS(529) = 0.5*B(115)

  JVS(530) = B(173)

  JVS(531) = -B(311)

  JVS(532) = -B(377)

  JVS(533) = -B(323)

  JVS(534) = -B(381)

  JVS(535) = -B(331)

  JVS(536) = -B(339)

  JVS(537) = -B(282)

  JVS(538) = -B(363)

  JVS(539) = -B(298)

  JVS(540) = -B(371)

  JVS(541) = -B(291)

  JVS(542) = B(2)-B(4)

  JVS(543) = -B(5)-B(13)-B(15)-B(30)-B(31)-B(51)-B(61)-B(283)-B(292)-B(299)-B(312)-B(324)-B(332)-B(340)-B(364)-B(372)&
               &-B(378)-B(382)

  JVS(544) = 0.25*B(124)

  JVS(545) = 0.25*B(184)

  JVS(546) = -B(52)

  JVS(547) = 0.25*B(142)

  JVS(548) = -B(14)

  JVS(549) = -B(16)

  JVS(550) = 0

  JVS(551) = -B(62)+0.25*B(125)+0.25*B(143)+0.25*B(162)+0.25*B(185)

  JVS(552) = 0.25*B(163)

  JVS(553) = 2*B(264)

  JVS(554) = B(121)

  JVS(555) = 0

  JVS(556) = 0.011*B(353)

  JVS(557) = B(313)+0.5*B(315)

  JVS(558) = B(316)+0.5*B(318)

  JVS(559) = B(259)+B(260)+B(262)

  JVS(560) = B(237)+B(239)

  JVS(561) = 0.123*B(339)

  JVS(562) = 0.67*B(288)

  JVS(563) = 0

  JVS(564) = 0.467*B(302)

  JVS(565) = 0.137*B(371)

  JVS(566) = 0.675*B(289)

  JVS(567) = B(227)+B(230)

  JVS(568) = 0

  JVS(569) = 0

  JVS(570) = 0

  JVS(571) = 0.492*B(240)+B(242)

  JVS(572) = 0.029*B(303)+0.667*B(305)

  JVS(573) = 0.123*B(340)+0.137*B(372)

  JVS(574) = -B(119)-B(122)-B(124)-B(126)-B(128)-B(130)-B(134)-2*B(136)-B(154)-B(174)

  JVS(575) = -B(135)

  JVS(576) = B(182)+B(186)+B(198)+B(200)+2*B(202)

  JVS(577) = B(228)+B(238)+0.492*B(241)+B(261)+0.675*B(290)+0.029*B(304)+B(314)+B(317)+0.011*B(354)

  JVS(578) = -B(155)+B(199)

  JVS(579) = -B(123)+B(183)

  JVS(580) = -B(120)

  JVS(581) = -B(127)+B(187)+B(231)+B(263)

  JVS(582) = -B(129)

  JVS(583) = -B(131)

  JVS(584) = -B(125)

  JVS(585) = -B(175)+B(201)

  JVS(586) = 0.035*B(347)

  JVS(587) = 0.07*B(351)

  JVS(588) = 0.347*B(355)

  JVS(589) = 0.009*B(359)

  JVS(590) = 0.143*B(353)

  JVS(591) = 0.011*B(357)

  JVS(592) = 0.016*B(379)+0.051*B(383)

  JVS(593) = 0.093*B(329)+0.008*B(331)+0.064*B(333)+0.01*B(335)

  JVS(594) = 0.25*B(337)+0.18*B(339)+0.25*B(341)

  JVS(595) = 0.09*B(361)+0.001*B(363)+0.176*B(365)

  JVS(596) = 0.041*B(296)+0.051*B(300)

  JVS(597) = 0.082*B(369)+0.002*B(371)+0.136*B(373)+0.001*B(375)

  JVS(598) = 0.025*B(289)

  JVS(599) = 0.173*B(306)+0.095*B(308)

  JVS(600) = 0.01*B(336)+0.001*B(376)

  JVS(601) = 0.001*B(232)

  JVS(602) = 0.042*B(240)

  JVS(603) = 0.07*B(303)+0.04*B(305)

  JVS(604) = 0.008*B(332)+0.18*B(340)+0.001*B(364)+0.002*B(372)

  JVS(605) = -B(134)

  JVS(606) = -B(106)-B(108)-B(110)-B(112)-B(114)-2*B(118)-B(135)-B(152)-B(172)-B(194)

  JVS(607) = -B(195)

  JVS(608) = 0.001*B(233)+0.042*B(241)+0.025*B(290)+0.041*B(297)+0.07*B(304)+0.173*B(307)+0.093*B(330)+0.25*B(338)+0.035&
               &*B(348)+0.07*B(352)+0.143*B(354)+0.347*B(356)+0.011*B(358)+0.009*B(360)+0.09*B(362)+0.082*B(370)+0.016&
               &*B(380)

  JVS(609) = -B(153)

  JVS(610) = -B(107)

  JVS(611) = 0

  JVS(612) = -B(113)+0.051*B(301)+0.064*B(334)+0.25*B(342)+0.176*B(366)+0.136*B(374)+0.051*B(384)

  JVS(613) = -B(111)

  JVS(614) = -B(115)

  JVS(615) = -B(109)

  JVS(616) = -B(173)

  JVS(617) = B(181)

  JVS(618) = 0.192*B(331)+0.24*B(335)

  JVS(619) = 0.5*B(280)+0.5*B(284)+0.33*B(288)

  JVS(620) = 0.289*B(296)+0.15*B(300)

  JVS(621) = 0

  JVS(622) = 0.3*B(295)

  JVS(623) = 0.24*B(336)

  JVS(624) = 0.192*B(332)

  JVS(625) = -B(196)

  JVS(626) = -B(194)

  JVS(627) = -B(179)-B(182)-B(184)-B(186)-B(188)-B(190)-B(195)-B(197)-B(198)-B(200)-2*B(202)

  JVS(628) = 0.5*B(281)+0.289*B(297)

  JVS(629) = -B(199)

  JVS(630) = -B(183)

  JVS(631) = -B(180)

  JVS(632) = -B(187)+0.5*B(285)+0.15*B(301)

  JVS(633) = -B(189)

  JVS(634) = -B(191)

  JVS(635) = -B(185)

  JVS(636) = -B(201)

  JVS(637) = -B(74)

  JVS(638) = 2*B(32)

  JVS(639) = -B(319)

  JVS(640) = -B(345)

  JVS(641) = 2*B(69)-B(70)

  JVS(642) = -B(347)

  JVS(643) = -B(245)

  JVS(644) = B(38)-B(40)

  JVS(645) = -B(351)

  JVS(646) = -B(355)

  JVS(647) = -B(359)

  JVS(648) = -B(243)

  JVS(649) = -0.65*B(247)+B(249)

  JVS(650) = 0.39*B(58)-B(59)

  JVS(651) = -B(353)

  JVS(652) = -B(357)

  JVS(653) = -B(313)

  JVS(654) = -B(316)

  JVS(655) = -B(269)

  JVS(656) = -B(309)+0.5*B(311)

  JVS(657) = -0.397*B(349)+0.5*B(377)

  JVS(658) = -B(275)

  JVS(659) = -0.34*B(250)+B(252)

  JVS(660) = -B(265)

  JVS(661) = -B(260)

  JVS(662) = -B(49)

  JVS(663) = -B(237)

  JVS(664) = -B(46)+B(48)

  JVS(665) = -B(321)+0.12*B(323)

  JVS(666) = -B(255)

  JVS(667) = -B(379)+0.32*B(381)

  JVS(668) = 0

  JVS(669) = -B(329)+0.266*B(331)

  JVS(670) = -B(337)+0.567*B(339)

  JVS(671) = -B(280)+0.208*B(282)+0.33*B(288)

  JVS(672) = -B(361)+0.155*B(363)

  JVS(673) = -B(296)+0.285*B(298)

  JVS(674) = -B(369)+0.378*B(371)

  JVS(675) = -B(289)+0.164*B(291)

  JVS(676) = -B(227)

  JVS(677) = -B(218)

  JVS(678) = -B(306)

  JVS(679) = 0

  JVS(680) = -B(232)

  JVS(681) = -B(240)

  JVS(682) = -B(303)

  JVS(683) = -B(51)+B(61)+0.208*B(283)+0.164*B(292)+0.285*B(299)+0.5*B(312)+0.12*B(324)+0.266*B(332)+0.567*B(340)+0.155&
               &*B(364)+0.378*B(372)+0.5*B(378)+0.32*B(382)

  JVS(684) = 0

  JVS(685) = 0

  JVS(686) = 0

  JVS(687) = -B(36)-B(41)-B(42)-B(44)-B(47)-B(50)-B(52)-B(60)-B(71)-B(72)-B(75)-B(76)-B(219)-B(228)-B(233)-B(238)-B(241)&
               &-B(244)-B(246)-0.65*B(248)-0.34*B(251)-B(256)-B(261)-B(266)-B(270)-B(276)-B(281)-B(290)-B(297)-B(304)-B(307)&
               &-B(310)-B(314)-B(317)-B(320)-B(322)-B(330)-B(338)-B(346)-B(348)-0.397*B(350)-B(352)-B(354)-B(356)-B(358)&
               &-B(360)-B(362)-B(370)-B(380)

  JVS(688) = 0

  JVS(689) = -B(37)+B(53)

  JVS(690) = -B(43)

  JVS(691) = -B(45)+0.8*B(66)

  JVS(692) = 0

  JVS(693) = 0

  JVS(694) = B(54)+B(62)+0.8*B(67)-B(73)

  JVS(695) = 0

  JVS(696) = B(139)

  JVS(697) = 0.37*B(255)+0.37*B(257)

  JVS(698) = 0

  JVS(699) = 0.201*B(339)

  JVS(700) = 0.1*B(282)

  JVS(701) = 0.048*B(298)+0.3*B(302)

  JVS(702) = 0.006*B(371)

  JVS(703) = 0.05*B(291)

  JVS(704) = 0

  JVS(705) = 0.965*B(232)+B(235)

  JVS(706) = 0.096*B(240)

  JVS(707) = 0.049*B(303)+0.333*B(305)

  JVS(708) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.201*B(340)+0.006*B(372)

  JVS(709) = -B(154)

  JVS(710) = -B(152)

  JVS(711) = -B(198)

  JVS(712) = 0.965*B(233)+0.096*B(241)+0.37*B(256)+0.049*B(304)

  JVS(713) = -B(137)-B(140)-B(142)-B(144)-B(146)-B(148)-B(153)-B(155)-2*B(156)-B(176)-B(199)

  JVS(714) = -B(141)

  JVS(715) = -B(138)

  JVS(716) = -B(145)+B(236)+0.37*B(258)

  JVS(717) = -B(147)

  JVS(718) = -B(149)

  JVS(719) = -B(143)

  JVS(720) = -B(177)

  JVS(721) = B(38)

  JVS(722) = -B(223)

  JVS(723) = -B(95)

  JVS(724) = 0

  JVS(725) = 0

  JVS(726) = 0

  JVS(727) = 0

  JVS(728) = 0

  JVS(729) = 0

  JVS(730) = -B(6)+B(9)

  JVS(731) = 0

  JVS(732) = 0

  JVS(733) = -B(13)

  JVS(734) = -B(122)

  JVS(735) = -B(106)

  JVS(736) = -B(182)

  JVS(737) = -B(36)

  JVS(738) = -B(140)

  JVS(739) = -B(7)-B(14)-B(17)-2*B(19)-B(37)-B(53)-B(78)-B(86)-B(96)-B(107)-B(123)-B(141)-B(160)-B(183)-B(224)

  JVS(740) = B(1)+B(10)+B(26)

  JVS(741) = -B(18)+B(27)+B(28)

  JVS(742) = -B(79)

  JVS(743) = -B(87)

  JVS(744) = -B(54)

  JVS(745) = -B(161)

  JVS(746) = B(121)

  JVS(747) = B(139)

  JVS(748) = B(159)

  JVS(749) = B(181)

  JVS(750) = B(23)

  JVS(751) = B(39)+B(40)

  JVS(752) = -B(203)

  JVS(753) = B(223)

  JVS(754) = -B(211)

  JVS(755) = B(57)+0.61*B(58)+B(59)

  JVS(756) = 0

  JVS(757) = B(48)

  JVS(758) = 0

  JVS(759) = -B(206)

  JVS(760) = 0.187*B(333)

  JVS(761) = B(95)+B(99)

  JVS(762) = 0.474*B(341)

  JVS(763) = 0

  JVS(764) = 0

  JVS(765) = 0

  JVS(766) = 0.391*B(373)

  JVS(767) = 0

  JVS(768) = 0

  JVS(769) = 0

  JVS(770) = 0.338*B(306)+B(308)

  JVS(771) = B(6)-B(9)-B(11)

  JVS(772) = 0

  JVS(773) = 0

  JVS(774) = 0

  JVS(775) = B(13)-B(15)

  JVS(776) = -B(119)+B(122)+B(126)

  JVS(777) = B(112)

  JVS(778) = -B(179)+B(182)+B(186)

  JVS(779) = B(41)-B(42)+B(44)+B(60)+0.338*B(307)

  JVS(780) = -B(137)+B(140)+B(144)

  JVS(781) = B(7)+B(14)+2*B(17)+2*B(19)+B(53)+B(78)+B(86)+B(96)+B(123)+B(141)+B(160)+B(183)+B(224)

  JVS(782) = -B(1)-B(10)-B(12)-B(16)-B(21)-B(43)-B(55)-B(120)-B(138)-B(157)-B(180)-B(204)-B(207)-B(212)

  JVS(783) = 2*B(18)-B(22)+B(29)+B(45)+0.8*B(66)+2*B(68)+B(82)+B(90)+B(100)+B(113)+B(127)+B(145)+B(164)+B(187)+0.187&
               &*B(334)+0.474*B(342)+0.391*B(374)

  JVS(784) = B(79)+B(83)

  JVS(785) = B(87)+B(91)

  JVS(786) = B(54)-B(56)+0.8*B(67)

  JVS(787) = -B(158)+B(161)+B(165)

  JVS(788) = B(23)

  JVS(789) = 0.39*B(58)

  JVS(790) = -B(271)

  JVS(791) = -B(273)

  JVS(792) = -B(278)

  JVS(793) = -B(267)

  JVS(794) = -B(262)

  JVS(795) = B(46)

  JVS(796) = -B(325)

  JVS(797) = -B(257)

  JVS(798) = -B(383)

  JVS(799) = 0

  JVS(800) = -B(333)

  JVS(801) = -B(99)

  JVS(802) = -B(341)

  JVS(803) = -B(284)

  JVS(804) = -B(365)

  JVS(805) = -B(300)

  JVS(806) = -B(373)

  JVS(807) = 0

  JVS(808) = -B(230)

  JVS(809) = -B(225)

  JVS(810) = 0

  JVS(811) = B(11)

  JVS(812) = -B(235)

  JVS(813) = 0

  JVS(814) = 0

  JVS(815) = B(15)

  JVS(816) = -B(126)

  JVS(817) = -B(112)

  JVS(818) = -B(186)

  JVS(819) = -B(44)+B(47)

  JVS(820) = -B(144)

  JVS(821) = -B(17)

  JVS(822) = B(12)+B(16)-B(21)-B(26)

  JVS(823) = -B(18)-B(22)-B(27)-B(28)-B(29)-B(45)-B(66)-2*B(68)-B(82)-B(90)-B(100)-B(113)-B(127)-B(145)-B(164)-B(187)&
               &-B(226)-B(231)-B(236)-B(258)-B(263)-B(268)-B(272)-B(274)-B(279)-B(285)-B(301)-B(326)-B(334)-B(342)-B(366)&
               &-B(374)-B(384)

  JVS(824) = -B(83)

  JVS(825) = -B(91)

  JVS(826) = -B(67)

  JVS(827) = -B(165)

  JVS(828) = B(319)

  JVS(829) = B(205)

  JVS(830) = 0.65*B(247)

  JVS(831) = 0.011*B(353)

  JVS(832) = B(239)

  JVS(833) = 0.3*B(327)

  JVS(834) = 0.26*B(381)

  JVS(835) = 0.25*B(335)

  JVS(836) = 0

  JVS(837) = 0.076*B(363)

  JVS(838) = 0.197*B(371)+0.03*B(373)

  JVS(839) = 0.3*B(295)

  JVS(840) = B(229)

  JVS(841) = 0

  JVS(842) = 0.3*B(328)+0.25*B(336)

  JVS(843) = 0

  JVS(844) = 0

  JVS(845) = 0

  JVS(846) = 0.076*B(364)+0.197*B(372)+0.26*B(382)

  JVS(847) = B(122)+B(126)-B(128)+2*B(136)+B(154)+B(174)+B(196)

  JVS(848) = -B(110)

  JVS(849) = -B(188)+B(197)

  JVS(850) = 0.65*B(248)+B(320)+0.011*B(354)

  JVS(851) = -B(146)+B(155)

  JVS(852) = -B(78)+B(123)

  JVS(853) = 0

  JVS(854) = -B(82)+B(127)+0.03*B(374)

  JVS(855) = -B(79)-B(80)-B(83)-2*B(84)-2*B(85)-B(92)-B(111)-B(129)-B(147)-B(166)-B(189)

  JVS(856) = -B(93)

  JVS(857) = -B(81)

  JVS(858) = -B(167)+B(175)

  JVS(859) = B(345)

  JVS(860) = 0.965*B(347)

  JVS(861) = 0.05*B(245)

  JVS(862) = 0.695*B(351)

  JVS(863) = 0.653*B(355)

  JVS(864) = 0.804*B(359)

  JVS(865) = 0.835*B(353)

  JVS(866) = 0.765*B(357)

  JVS(867) = B(315)

  JVS(868) = B(318)

  JVS(869) = 0.76*B(269)

  JVS(870) = B(309)

  JVS(871) = 0.1*B(349)

  JVS(872) = 0.34*B(250)

  JVS(873) = 0.76*B(265)

  JVS(874) = B(321)+B(325)+0.2*B(327)

  JVS(875) = 0.984*B(379)+0.949*B(383)

  JVS(876) = 0

  JVS(877) = 0.907*B(329)+0.066*B(331)+0.749*B(333)

  JVS(878) = 0.75*B(337)+0.031*B(339)+0.276*B(341)

  JVS(879) = 0.5*B(280)+0.1*B(282)+0.5*B(284)+0.33*B(288)

  JVS(880) = 0.91*B(361)+0.022*B(363)+0.824*B(365)

  JVS(881) = 0.67*B(296)+0.048*B(298)+0.799*B(300)

  JVS(882) = 0.918*B(369)+0.033*B(371)+0.442*B(373)+0.012*B(375)

  JVS(883) = 0.3*B(289)+0.05*B(291)

  JVS(884) = 0.376*B(306)+0.564*B(308)

  JVS(885) = 0.2*B(328)+0.012*B(376)

  JVS(886) = 0.034*B(232)+B(234)

  JVS(887) = 0.37*B(240)+B(242)

  JVS(888) = 0.473*B(303)+0.96*B(305)

  JVS(889) = 0.1*B(283)+0.05*B(292)+0.048*B(299)+0.066*B(332)+0.031*B(340)+0.022*B(364)+0.033*B(372)

  JVS(890) = -B(130)+B(154)

  JVS(891) = -B(114)

  JVS(892) = -B(190)+B(198)

  JVS(893) = 0.034*B(233)+0.37*B(241)+0.05*B(246)+0.34*B(251)+0.76*B(266)+0.76*B(270)+0.5*B(281)+0.3*B(290)+0.67*B(297)&
               &+0.473*B(304)+0.376*B(307)+B(310)+B(322)+0.907*B(330)+0.75*B(338)+B(346)+0.965*B(348)+0.1*B(350)+0.695&
               &*B(352)+0.835*B(354)+0.653*B(356)+0.765*B(358)+0.804*B(360)+0.91*B(362)+0.918*B(370)+0.984*B(380)

  JVS(894) = B(140)+B(144)-B(148)+B(155)+2*B(156)+B(176)+B(199)

  JVS(895) = -B(86)+B(141)

  JVS(896) = 0

  JVS(897) = -B(90)+B(145)+0.5*B(285)+0.799*B(301)+B(326)+0.749*B(334)+0.276*B(342)+0.824*B(366)+0.442*B(374)+0.949&
               &*B(384)

  JVS(898) = -B(92)

  JVS(899) = -B(87)-B(88)-B(91)-B(93)-2*B(94)-B(115)-B(131)-B(149)-B(168)-B(191)

  JVS(900) = -B(89)

  JVS(901) = -B(169)+B(177)

  JVS(902) = B(74)

  JVS(903) = B(70)

  JVS(904) = 0.95*B(245)

  JVS(905) = B(39)

  JVS(906) = 0.187*B(359)

  JVS(907) = B(243)

  JVS(908) = B(249)

  JVS(909) = B(222)+B(223)

  JVS(910) = -B(213)

  JVS(911) = B(57)+0.61*B(58)

  JVS(912) = 0.224*B(357)

  JVS(913) = 0.5*B(315)

  JVS(914) = 0.5*B(318)

  JVS(915) = 1.5*B(311)

  JVS(916) = 0.297*B(349)+1.5*B(377)

  JVS(917) = 0

  JVS(918) = B(252)

  JVS(919) = B(259)

  JVS(920) = B(49)

  JVS(921) = 0.12*B(323)+0.5*B(327)

  JVS(922) = 2*B(253)+0.63*B(255)+0.63*B(257)

  JVS(923) = 0.06*B(381)

  JVS(924) = -B(208)

  JVS(925) = 0

  JVS(926) = 0.033*B(339)

  JVS(927) = 0.008*B(282)+0.34*B(288)

  JVS(928) = 0.056*B(363)

  JVS(929) = 0.4*B(298)+1.233*B(302)

  JVS(930) = 0.003*B(371)+0.013*B(375)

  JVS(931) = 0.064*B(291)

  JVS(932) = B(229)

  JVS(933) = 2*B(216)+B(218)-B(220)+B(225)

  JVS(934) = 0.113*B(306)+0.341*B(308)

  JVS(935) = 0.5*B(328)+0.013*B(376)

  JVS(936) = B(234)

  JVS(937) = 0

  JVS(938) = 0.379*B(303)

  JVS(939) = B(51)-B(61)+0.008*B(283)+0.064*B(292)+0.4*B(299)+1.5*B(312)+0.12*B(324)+0.033*B(340)+0.056*B(364)+0.003&
               &*B(372)+1.5*B(378)+0.06*B(382)

  JVS(940) = -B(124)

  JVS(941) = -B(108)+B(110)+B(112)+B(114)+B(118)

  JVS(942) = -B(184)

  JVS(943) = B(44)+B(50)+B(52)+B(71)-B(72)+B(75)+B(76)+B(219)+B(244)+0.95*B(246)+0.63*B(256)+0.379*B(304)+0.113*B(307)&
               &+0.297*B(350)+0.224*B(358)+0.187*B(360)

  JVS(944) = -B(142)

  JVS(945) = -B(53)+B(78)+B(86)+B(224)

  JVS(946) = -B(55)

  JVS(947) = B(45)-B(66)+B(82)+B(90)+B(113)+B(226)+0.63*B(258)

  JVS(948) = B(79)-B(80)+B(83)+2*B(85)+B(92)+B(111)

  JVS(949) = B(87)-B(88)+B(91)+B(93)+B(94)+B(115)

  JVS(950) = -B(54)-B(56)-B(62)-2*B(63)-2*B(64)-B(67)-B(73)-B(81)-B(89)-B(109)-B(125)-B(143)-B(162)-B(185)-B(209)-B(214)&
               &-B(221)-B(388)

  JVS(951) = -B(163)

  JVS(952) = B(159)

  JVS(953) = B(275)+B(278)

  JVS(954) = 0

  JVS(955) = 0

  JVS(956) = 0

  JVS(957) = -B(174)

  JVS(958) = -B(172)

  JVS(959) = -B(200)

  JVS(960) = B(276)

  JVS(961) = -B(176)

  JVS(962) = -B(160)

  JVS(963) = -B(157)

  JVS(964) = -B(164)+B(279)

  JVS(965) = -B(166)

  JVS(966) = -B(168)

  JVS(967) = -B(162)

  JVS(968) = -B(158)-B(161)-B(163)-B(165)-B(167)-B(169)-B(173)-B(175)-B(177)-2*B(178)-B(201)
      
END SUBROUTINE saprc99_Jac_SP














SUBROUTINE saprc99_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(968), W(79), a
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
      
END SUBROUTINE saprc99_KppDecomp



SUBROUTINE saprc99_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(968), W(79), a
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
      
END SUBROUTINE saprc99_KppDecompCmplx


SUBROUTINE saprc99_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(968), X(79), sum

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
      
END SUBROUTINE saprc99_KppSolveIndirect


SUBROUTINE saprc99_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(968), X(79), sum

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
      
END SUBROUTINE saprc99_KppSolveCmplx













SUBROUTINE saprc99_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(25) = X(25)-JVS(109)*X(24)
  X(35) = X(35)-JVS(143)*X(27)-JVS(144)*X(34)
  X(36) = X(36)-JVS(147)*X(27)-JVS(148)*X(34)
  X(37) = X(37)-JVS(151)*X(27)-JVS(152)*X(34)
  X(38) = X(38)-JVS(156)*X(27)-JVS(157)*X(34)
  X(40) = X(40)-JVS(164)*X(31)
  X(41) = X(41)-JVS(170)*X(27)-JVS(171)*X(34)
  X(43) = X(43)-JVS(182)*X(34)
  X(44) = X(44)-JVS(188)*X(27)-JVS(189)*X(34)-JVS(190)*X(35)-JVS(191)*X(36)-JVS(192)*X(37)
  X(45) = X(45)-JVS(200)*X(33)-JVS(201)*X(35)-JVS(202)*X(36)-JVS(203)*X(38)-JVS(204)*X(39)-JVS(205)*X(44)
  X(46) = X(46)-JVS(224)*X(20)-JVS(225)*X(24)-JVS(226)*X(25)-JVS(227)*X(26)-JVS(228)*X(33)
  X(47) = X(47)-JVS(238)*X(22)-JVS(239)*X(37)-JVS(240)*X(40)-JVS(241)*X(41)-JVS(242)*X(43)-JVS(243)*X(44)
  X(49) = X(49)-JVS(264)*X(27)-JVS(265)*X(34)-JVS(266)*X(35)-JVS(267)*X(36)-JVS(268)*X(38)-JVS(269)*X(39)-JVS(270)*X(43)&
            &-JVS(271)*X(48)
  X(51) = X(51)-JVS(286)*X(37)-JVS(287)*X(43)
  X(53) = X(53)-JVS(303)*X(24)-JVS(304)*X(26)-JVS(305)*X(33)-JVS(306)*X(35)-JVS(307)*X(36)-JVS(308)*X(46)-JVS(309)*X(52)
  X(55) = X(55)-JVS(337)*X(52)
  X(57) = X(57)-JVS(349)*X(52)
  X(59) = X(59)-JVS(361)*X(52)-JVS(362)*X(58)
  X(60) = X(60)-JVS(368)*X(13)-JVS(369)*X(21)-JVS(370)*X(24)-JVS(371)*X(26)-JVS(372)*X(33)-JVS(373)*X(48)-JVS(374)*X(50)&
            &-JVS(375)*X(56)-JVS(376)*X(57)-JVS(377)*X(58)
  X(61) = X(61)-JVS(392)*X(21)-JVS(393)*X(24)-JVS(394)*X(26)-JVS(395)*X(28)-JVS(396)*X(29)-JVS(397)*X(30)-JVS(398)*X(33)&
            &-JVS(399)*X(39)-JVS(400)*X(46)-JVS(401)*X(48)-JVS(402)*X(49)-JVS(403)*X(50)-JVS(404)*X(51)-JVS(405)*X(52)&
            &-JVS(406)*X(54)-JVS(407)*X(55)-JVS(408)*X(56)-JVS(409)*X(57)-JVS(410)*X(58)-JVS(411)*X(59)
  X(62) = X(62)-JVS(430)*X(25)-JVS(431)*X(54)-JVS(432)*X(56)-JVS(433)*X(57)-JVS(434)*X(58)
  X(63) = X(63)-JVS(443)*X(11)-JVS(444)*X(48)-JVS(445)*X(50)-JVS(446)*X(52)-JVS(447)*X(54)-JVS(448)*X(55)-JVS(449)*X(56)&
            &-JVS(450)*X(58)-JVS(451)*X(59)
  X(64) = X(64)-JVS(458)*X(20)-JVS(459)*X(24)-JVS(460)*X(26)-JVS(461)*X(33)-JVS(462)*X(35)-JVS(463)*X(36)-JVS(464)*X(38)&
            &-JVS(465)*X(42)-JVS(466)*X(48)-JVS(467)*X(50)-JVS(468)*X(54)-JVS(469)*X(55)-JVS(470)*X(56)-JVS(471)*X(57)&
            &-JVS(472)*X(58)-JVS(473)*X(59)-JVS(474)*X(62)-JVS(475)*X(63)
  X(65) = X(65)-JVS(487)*X(24)-JVS(488)*X(26)-JVS(489)*X(33)-JVS(490)*X(50)-JVS(491)*X(55)-JVS(492)*X(56)-JVS(493)*X(57)&
            &-JVS(494)*X(58)-JVS(495)*X(59)-JVS(496)*X(62)-JVS(497)*X(63)
  X(66) = X(66)-JVS(508)*X(26)-JVS(509)*X(33)-JVS(510)*X(34)-JVS(511)*X(52)-JVS(512)*X(54)-JVS(513)*X(56)-JVS(514)*X(57)&
            &-JVS(515)*X(58)-JVS(516)*X(59)-JVS(517)*X(62)-JVS(518)*X(63)
  X(67) = X(67)-JVS(531)*X(38)-JVS(532)*X(39)-JVS(533)*X(48)-JVS(534)*X(50)-JVS(535)*X(52)-JVS(536)*X(54)-JVS(537)*X(55)&
            &-JVS(538)*X(56)-JVS(539)*X(57)-JVS(540)*X(58)-JVS(541)*X(59)-JVS(542)*X(63)
  X(68) = X(68)-JVS(553)*X(14)-JVS(554)*X(15)-JVS(555)*X(27)-JVS(556)*X(33)-JVS(557)*X(35)-JVS(558)*X(36)-JVS(559)*X(44)&
            &-JVS(560)*X(46)-JVS(561)*X(54)-JVS(562)*X(55)-JVS(563)*X(56)-JVS(564)*X(57)-JVS(565)*X(58)-JVS(566)*X(59)&
            &-JVS(567)*X(60)-JVS(568)*X(62)-JVS(569)*X(63)-JVS(570)*X(64)-JVS(571)*X(65)-JVS(572)*X(66)-JVS(573)*X(67)
  X(69) = X(69)-JVS(586)*X(20)-JVS(587)*X(24)-JVS(588)*X(26)-JVS(589)*X(27)-JVS(590)*X(33)-JVS(591)*X(34)-JVS(592)*X(50)&
            &-JVS(593)*X(52)-JVS(594)*X(54)-JVS(595)*X(56)-JVS(596)*X(57)-JVS(597)*X(58)-JVS(598)*X(59)-JVS(599)*X(62)&
            &-JVS(600)*X(63)-JVS(601)*X(64)-JVS(602)*X(65)-JVS(603)*X(66)-JVS(604)*X(67)-JVS(605)*X(68)
  X(70) = X(70)-JVS(617)*X(18)-JVS(618)*X(52)-JVS(619)*X(55)-JVS(620)*X(57)-JVS(621)*X(58)-JVS(622)*X(59)-JVS(623)*X(63)&
            &-JVS(624)*X(67)-JVS(625)*X(68)-JVS(626)*X(69)
  X(71) = X(71)-JVS(637)*X(10)-JVS(638)*X(11)-JVS(639)*X(12)-JVS(640)*X(13)-JVS(641)*X(19)-JVS(642)*X(20)-JVS(643)*X(21)&
            &-JVS(644)*X(23)-JVS(645)*X(24)-JVS(646)*X(26)-JVS(647)*X(27)-JVS(648)*X(28)-JVS(649)*X(29)-JVS(650)*X(32)&
            &-JVS(651)*X(33)-JVS(652)*X(34)-JVS(653)*X(35)-JVS(654)*X(36)-JVS(655)*X(37)-JVS(656)*X(38)-JVS(657)*X(39)&
            &-JVS(658)*X(41)-JVS(659)*X(42)-JVS(660)*X(43)-JVS(661)*X(44)-JVS(662)*X(45)-JVS(663)*X(46)-JVS(664)*X(47)&
            &-JVS(665)*X(48)-JVS(666)*X(49)-JVS(667)*X(50)-JVS(668)*X(51)-JVS(669)*X(52)-JVS(670)*X(54)-JVS(671)*X(55)&
            &-JVS(672)*X(56)-JVS(673)*X(57)-JVS(674)*X(58)-JVS(675)*X(59)-JVS(676)*X(60)-JVS(677)*X(61)-JVS(678)*X(62)&
            &-JVS(679)*X(63)-JVS(680)*X(64)-JVS(681)*X(65)-JVS(682)*X(66)-JVS(683)*X(67)-JVS(684)*X(68)-JVS(685)*X(69)&
            &-JVS(686)*X(70)
  X(72) = X(72)-JVS(696)*X(16)-JVS(697)*X(49)-JVS(698)*X(51)-JVS(699)*X(54)-JVS(700)*X(55)-JVS(701)*X(57)-JVS(702)*X(58)&
            &-JVS(703)*X(59)-JVS(704)*X(63)-JVS(705)*X(64)-JVS(706)*X(65)-JVS(707)*X(66)-JVS(708)*X(67)-JVS(709)*X(68)&
            &-JVS(710)*X(69)-JVS(711)*X(70)-JVS(712)*X(71)
  X(73) = X(73)-JVS(721)*X(23)-JVS(722)*X(30)-JVS(723)*X(53)-JVS(724)*X(54)-JVS(725)*X(56)-JVS(726)*X(58)-JVS(727)*X(59)&
            &-JVS(728)*X(61)-JVS(729)*X(62)-JVS(730)*X(63)-JVS(731)*X(65)-JVS(732)*X(66)-JVS(733)*X(67)-JVS(734)*X(68)&
            &-JVS(735)*X(69)-JVS(736)*X(70)-JVS(737)*X(71)-JVS(738)*X(72)
  X(74) = X(74)-JVS(746)*X(15)-JVS(747)*X(16)-JVS(748)*X(17)-JVS(749)*X(18)-JVS(750)*X(22)-JVS(751)*X(23)-JVS(752)*X(25)&
            &-JVS(753)*X(30)-JVS(754)*X(31)-JVS(755)*X(32)-JVS(756)*X(40)-JVS(757)*X(47)-JVS(758)*X(49)-JVS(759)*X(51)&
            &-JVS(760)*X(52)-JVS(761)*X(53)-JVS(762)*X(54)-JVS(763)*X(55)-JVS(764)*X(56)-JVS(765)*X(57)-JVS(766)*X(58)&
            &-JVS(767)*X(59)-JVS(768)*X(60)-JVS(769)*X(61)-JVS(770)*X(62)-JVS(771)*X(63)-JVS(772)*X(64)-JVS(773)*X(65)&
            &-JVS(774)*X(66)-JVS(775)*X(67)-JVS(776)*X(68)-JVS(777)*X(69)-JVS(778)*X(70)-JVS(779)*X(71)-JVS(780)*X(72)&
            &-JVS(781)*X(73)
  X(75) = X(75)-JVS(788)*X(22)-JVS(789)*X(32)-JVS(790)*X(37)-JVS(791)*X(40)-JVS(792)*X(41)-JVS(793)*X(43)-JVS(794)*X(44)&
            &-JVS(795)*X(47)-JVS(796)*X(48)-JVS(797)*X(49)-JVS(798)*X(50)-JVS(799)*X(51)-JVS(800)*X(52)-JVS(801)*X(53)&
            &-JVS(802)*X(54)-JVS(803)*X(55)-JVS(804)*X(56)-JVS(805)*X(57)-JVS(806)*X(58)-JVS(807)*X(59)-JVS(808)*X(60)&
            &-JVS(809)*X(61)-JVS(810)*X(62)-JVS(811)*X(63)-JVS(812)*X(64)-JVS(813)*X(65)-JVS(814)*X(66)-JVS(815)*X(67)&
            &-JVS(816)*X(68)-JVS(817)*X(69)-JVS(818)*X(70)-JVS(819)*X(71)-JVS(820)*X(72)-JVS(821)*X(73)-JVS(822)*X(74)
  X(76) = X(76)-JVS(828)*X(12)-JVS(829)*X(25)-JVS(830)*X(29)-JVS(831)*X(33)-JVS(832)*X(46)-JVS(833)*X(48)-JVS(834)*X(50)&
            &-JVS(835)*X(52)-JVS(836)*X(54)-JVS(837)*X(56)-JVS(838)*X(58)-JVS(839)*X(59)-JVS(840)*X(60)-JVS(841)*X(62)&
            &-JVS(842)*X(63)-JVS(843)*X(64)-JVS(844)*X(65)-JVS(845)*X(66)-JVS(846)*X(67)-JVS(847)*X(68)-JVS(848)*X(69)&
            &-JVS(849)*X(70)-JVS(850)*X(71)-JVS(851)*X(72)-JVS(852)*X(73)-JVS(853)*X(74)-JVS(854)*X(75)
  X(77) = X(77)-JVS(859)*X(13)-JVS(860)*X(20)-JVS(861)*X(21)-JVS(862)*X(24)-JVS(863)*X(26)-JVS(864)*X(27)-JVS(865)*X(33)&
            &-JVS(866)*X(34)-JVS(867)*X(35)-JVS(868)*X(36)-JVS(869)*X(37)-JVS(870)*X(38)-JVS(871)*X(39)-JVS(872)*X(42)&
            &-JVS(873)*X(43)-JVS(874)*X(48)-JVS(875)*X(50)-JVS(876)*X(51)-JVS(877)*X(52)-JVS(878)*X(54)-JVS(879)*X(55)&
            &-JVS(880)*X(56)-JVS(881)*X(57)-JVS(882)*X(58)-JVS(883)*X(59)-JVS(884)*X(62)-JVS(885)*X(63)-JVS(886)*X(64)&
            &-JVS(887)*X(65)-JVS(888)*X(66)-JVS(889)*X(67)-JVS(890)*X(68)-JVS(891)*X(69)-JVS(892)*X(70)-JVS(893)*X(71)&
            &-JVS(894)*X(72)-JVS(895)*X(73)-JVS(896)*X(74)-JVS(897)*X(75)-JVS(898)*X(76)
  X(78) = X(78)-JVS(902)*X(10)-JVS(903)*X(19)-JVS(904)*X(21)-JVS(905)*X(23)-JVS(906)*X(27)-JVS(907)*X(28)-JVS(908)*X(29)&
            &-JVS(909)*X(30)-JVS(910)*X(31)-JVS(911)*X(32)-JVS(912)*X(34)-JVS(913)*X(35)-JVS(914)*X(36)-JVS(915)*X(38)&
            &-JVS(916)*X(39)-JVS(917)*X(40)-JVS(918)*X(42)-JVS(919)*X(44)-JVS(920)*X(45)-JVS(921)*X(48)-JVS(922)*X(49)&
            &-JVS(923)*X(50)-JVS(924)*X(51)-JVS(925)*X(52)-JVS(926)*X(54)-JVS(927)*X(55)-JVS(928)*X(56)-JVS(929)*X(57)&
            &-JVS(930)*X(58)-JVS(931)*X(59)-JVS(932)*X(60)-JVS(933)*X(61)-JVS(934)*X(62)-JVS(935)*X(63)-JVS(936)*X(64)&
            &-JVS(937)*X(65)-JVS(938)*X(66)-JVS(939)*X(67)-JVS(940)*X(68)-JVS(941)*X(69)-JVS(942)*X(70)-JVS(943)*X(71)&
            &-JVS(944)*X(72)-JVS(945)*X(73)-JVS(946)*X(74)-JVS(947)*X(75)-JVS(948)*X(76)-JVS(949)*X(77)
  X(79) = X(79)-JVS(952)*X(17)-JVS(953)*X(41)-JVS(954)*X(58)-JVS(955)*X(63)-JVS(956)*X(67)-JVS(957)*X(68)-JVS(958)*X(69)&
            &-JVS(959)*X(70)-JVS(960)*X(71)-JVS(961)*X(72)-JVS(962)*X(73)-JVS(963)*X(74)-JVS(964)*X(75)-JVS(965)*X(76)&
            &-JVS(966)*X(77)-JVS(967)*X(78)
  X(79) = X(79)/JVS(968)
  X(78) = (X(78)-JVS(951)*X(79))/(JVS(950))
  X(77) = (X(77)-JVS(900)*X(78)-JVS(901)*X(79))/(JVS(899))
  X(76) = (X(76)-JVS(856)*X(77)-JVS(857)*X(78)-JVS(858)*X(79))/(JVS(855))
  X(75) = (X(75)-JVS(824)*X(76)-JVS(825)*X(77)-JVS(826)*X(78)-JVS(827)*X(79))/(JVS(823))
  X(74) = (X(74)-JVS(783)*X(75)-JVS(784)*X(76)-JVS(785)*X(77)-JVS(786)*X(78)-JVS(787)*X(79))/(JVS(782))
  X(73) = (X(73)-JVS(740)*X(74)-JVS(741)*X(75)-JVS(742)*X(76)-JVS(743)*X(77)-JVS(744)*X(78)-JVS(745)*X(79))/(JVS(739))
  X(72) = (X(72)-JVS(714)*X(73)-JVS(715)*X(74)-JVS(716)*X(75)-JVS(717)*X(76)-JVS(718)*X(77)-JVS(719)*X(78)-JVS(720)&
            &*X(79))/(JVS(713))
  X(71) = (X(71)-JVS(688)*X(72)-JVS(689)*X(73)-JVS(690)*X(74)-JVS(691)*X(75)-JVS(692)*X(76)-JVS(693)*X(77)-JVS(694)&
            &*X(78)-JVS(695)*X(79))/(JVS(687))
  X(70) = (X(70)-JVS(628)*X(71)-JVS(629)*X(72)-JVS(630)*X(73)-JVS(631)*X(74)-JVS(632)*X(75)-JVS(633)*X(76)-JVS(634)&
            &*X(77)-JVS(635)*X(78)-JVS(636)*X(79))/(JVS(627))
  X(69) = (X(69)-JVS(607)*X(70)-JVS(608)*X(71)-JVS(609)*X(72)-JVS(610)*X(73)-JVS(611)*X(74)-JVS(612)*X(75)-JVS(613)&
            &*X(76)-JVS(614)*X(77)-JVS(615)*X(78)-JVS(616)*X(79))/(JVS(606))
  X(68) = (X(68)-JVS(575)*X(69)-JVS(576)*X(70)-JVS(577)*X(71)-JVS(578)*X(72)-JVS(579)*X(73)-JVS(580)*X(74)-JVS(581)&
            &*X(75)-JVS(582)*X(76)-JVS(583)*X(77)-JVS(584)*X(78)-JVS(585)*X(79))/(JVS(574))
  X(67) = (X(67)-JVS(544)*X(68)-JVS(545)*X(70)-JVS(546)*X(71)-JVS(547)*X(72)-JVS(548)*X(73)-JVS(549)*X(74)-JVS(550)&
            &*X(75)-JVS(551)*X(78)-JVS(552)*X(79))/(JVS(543))
  X(66) = (X(66)-JVS(520)*X(67)-JVS(521)*X(68)-JVS(522)*X(69)-JVS(523)*X(71)-JVS(524)*X(72)-JVS(525)*X(73)-JVS(526)&
            &*X(74)-JVS(527)*X(75)-JVS(528)*X(76)-JVS(529)*X(77)-JVS(530)*X(79))/(JVS(519))
  X(65) = (X(65)-JVS(499)*X(66)-JVS(500)*X(67)-JVS(501)*X(69)-JVS(502)*X(71)-JVS(503)*X(73)-JVS(504)*X(74)-JVS(505)&
            &*X(75)-JVS(506)*X(76)-JVS(507)*X(77))/(JVS(498))
  X(64) = (X(64)-JVS(477)*X(65)-JVS(478)*X(66)-JVS(479)*X(67)-JVS(480)*X(69)-JVS(481)*X(71)-JVS(482)*X(73)-JVS(483)&
            &*X(74)-JVS(484)*X(75)-JVS(485)*X(77)-JVS(486)*X(78))/(JVS(476))
  X(63) = (X(63)-JVS(453)*X(67)-JVS(454)*X(71)-JVS(455)*X(73)-JVS(456)*X(74)-JVS(457)*X(75))/(JVS(452))
  X(62) = (X(62)-JVS(436)*X(63)-JVS(437)*X(67)-JVS(438)*X(69)-JVS(439)*X(71)-JVS(440)*X(73)-JVS(441)*X(74)-JVS(442)&
            &*X(75))/(JVS(435))
  X(61) = (X(61)-JVS(413)*X(62)-JVS(414)*X(63)-JVS(415)*X(65)-JVS(416)*X(66)-JVS(417)*X(67)-JVS(418)*X(68)-JVS(419)&
            &*X(69)-JVS(420)*X(70)-JVS(421)*X(71)-JVS(422)*X(72)-JVS(423)*X(73)-JVS(424)*X(74)-JVS(425)*X(75)-JVS(426)*X(76)&
            &-JVS(427)*X(77)-JVS(428)*X(78)-JVS(429)*X(79))/(JVS(412))
  X(60) = (X(60)-JVS(379)*X(62)-JVS(380)*X(63)-JVS(381)*X(64)-JVS(382)*X(65)-JVS(383)*X(66)-JVS(384)*X(67)-JVS(385)&
            &*X(68)-JVS(386)*X(70)-JVS(387)*X(71)-JVS(388)*X(72)-JVS(389)*X(73)-JVS(390)*X(75)-JVS(391)*X(79))/(JVS(378))
  X(59) = (X(59)-JVS(364)*X(63)-JVS(365)*X(67)-JVS(366)*X(71)-JVS(367)*X(75))/(JVS(363))
  X(58) = (X(58)-JVS(357)*X(63)-JVS(358)*X(67)-JVS(359)*X(71)-JVS(360)*X(75))/(JVS(356))
  X(57) = (X(57)-JVS(351)*X(58)-JVS(352)*X(63)-JVS(353)*X(67)-JVS(354)*X(71)-JVS(355)*X(75))/(JVS(350))
  X(56) = (X(56)-JVS(345)*X(63)-JVS(346)*X(67)-JVS(347)*X(71)-JVS(348)*X(75))/(JVS(344))
  X(55) = (X(55)-JVS(339)*X(58)-JVS(340)*X(63)-JVS(341)*X(67)-JVS(342)*X(71)-JVS(343)*X(75))/(JVS(338))
  X(54) = (X(54)-JVS(333)*X(63)-JVS(334)*X(67)-JVS(335)*X(71)-JVS(336)*X(75))/(JVS(332))
  X(53) = (X(53)-JVS(311)*X(54)-JVS(312)*X(56)-JVS(313)*X(58)-JVS(314)*X(59)-JVS(315)*X(62)-JVS(316)*X(63)-JVS(317)&
            &*X(65)-JVS(318)*X(66)-JVS(319)*X(67)-JVS(320)*X(68)-JVS(321)*X(69)-JVS(322)*X(70)-JVS(323)*X(71)-JVS(324)*X(72)&
            &-JVS(325)*X(73)-JVS(326)*X(74)-JVS(327)*X(75)-JVS(328)*X(76)-JVS(329)*X(77)-JVS(330)*X(78)-JVS(331)*X(79))&
            &/(JVS(310))
  X(52) = (X(52)-JVS(299)*X(63)-JVS(300)*X(67)-JVS(301)*X(71)-JVS(302)*X(75))/(JVS(298))
  X(51) = (X(51)-JVS(289)*X(68)-JVS(290)*X(70)-JVS(291)*X(71)-JVS(292)*X(72)-JVS(293)*X(73)-JVS(294)*X(74)-JVS(295)&
            &*X(75)-JVS(296)*X(78)-JVS(297)*X(79))/(JVS(288))
  X(50) = (X(50)-JVS(282)*X(63)-JVS(283)*X(67)-JVS(284)*X(71)-JVS(285)*X(75))/(JVS(281))
  X(49) = (X(49)-JVS(273)*X(51)-JVS(274)*X(54)-JVS(275)*X(57)-JVS(276)*X(63)-JVS(277)*X(67)-JVS(278)*X(71)-JVS(279)&
            &*X(75)-JVS(280)*X(78))/(JVS(272))
  X(48) = (X(48)-JVS(260)*X(63)-JVS(261)*X(67)-JVS(262)*X(71)-JVS(263)*X(75))/(JVS(259))
  X(47) = (X(47)-JVS(245)*X(49)-JVS(246)*X(51)-JVS(247)*X(55)-JVS(248)*X(57)-JVS(249)*X(58)-JVS(250)*X(59)-JVS(251)&
            &*X(60)-JVS(252)*X(61)-JVS(253)*X(64)-JVS(254)*X(67)-JVS(255)*X(71)-JVS(256)*X(74)-JVS(257)*X(75)-JVS(258)&
            &*X(78))/(JVS(244))
  X(46) = (X(46)-JVS(230)*X(54)-JVS(231)*X(56)-JVS(232)*X(58)-JVS(233)*X(62)-JVS(234)*X(67)-JVS(235)*X(71)-JVS(236)&
            &*X(74)-JVS(237)*X(75))/(JVS(229))
  X(45) = (X(45)-JVS(207)*X(48)-JVS(208)*X(49)-JVS(209)*X(50)-JVS(210)*X(52)-JVS(211)*X(54)-JVS(212)*X(55)-JVS(213)&
            &*X(56)-JVS(214)*X(57)-JVS(215)*X(58)-JVS(216)*X(59)-JVS(217)*X(60)-JVS(218)*X(61)-JVS(219)*X(63)-JVS(220)*X(64)&
            &-JVS(221)*X(67)-JVS(222)*X(71)-JVS(223)*X(75))/(JVS(206))
  X(44) = (X(44)-JVS(194)*X(55)-JVS(195)*X(57)-JVS(196)*X(59)-JVS(197)*X(67)-JVS(198)*X(71)-JVS(199)*X(75))/(JVS(193))
  X(43) = (X(43)-JVS(184)*X(51)-JVS(185)*X(71)-JVS(186)*X(75)-JVS(187)*X(78))/(JVS(183))
  X(42) = (X(42)-JVS(178)*X(69)-JVS(179)*X(71)-JVS(180)*X(77)-JVS(181)*X(78))/(JVS(177))
  X(41) = (X(41)-JVS(173)*X(58)-JVS(174)*X(67)-JVS(175)*X(71)-JVS(176)*X(75))/(JVS(172))
  X(40) = (X(40)-JVS(166)*X(51)-JVS(167)*X(74)-JVS(168)*X(75)-JVS(169)*X(78))/(JVS(165))
  X(39) = (X(39)-JVS(162)*X(67)-JVS(163)*X(71))/(JVS(161))
  X(38) = (X(38)-JVS(159)*X(67)-JVS(160)*X(71))/(JVS(158))
  X(37) = (X(37)-JVS(154)*X(71)-JVS(155)*X(75))/(JVS(153))
  X(36) = (X(36)-JVS(150)*X(71))/(JVS(149))
  X(35) = (X(35)-JVS(146)*X(71))/(JVS(145))
  X(34) = (X(34)-JVS(142)*X(71))/(JVS(141))
  X(33) = (X(33)-JVS(140)*X(71))/(JVS(139))
  X(32) = (X(32)-JVS(136)*X(71)-JVS(137)*X(74)-JVS(138)*X(78))/(JVS(135))
  X(31) = (X(31)-JVS(131)*X(40)-JVS(132)*X(74)-JVS(133)*X(75)-JVS(134)*X(78))/(JVS(130))
  X(30) = (X(30)-JVS(127)*X(61)-JVS(128)*X(73)-JVS(129)*X(78))/(JVS(126))
  X(29) = (X(29)-JVS(123)*X(71)-JVS(124)*X(76)-JVS(125)*X(78))/(JVS(122))
  X(28) = (X(28)-JVS(118)*X(69)-JVS(119)*X(71)-JVS(120)*X(76)-JVS(121)*X(77))/(JVS(117))
  X(27) = (X(27)-JVS(116)*X(71))/(JVS(115))
  X(26) = (X(26)-JVS(114)*X(71))/(JVS(113))
  X(25) = (X(25)-JVS(111)*X(71)-JVS(112)*X(74))/(JVS(110))
  X(24) = (X(24)-JVS(108)*X(71))/(JVS(107))
  X(23) = (X(23)-JVS(105)*X(71)-JVS(106)*X(73))/(JVS(104))
  X(22) = (X(22)-JVS(102)*X(74)-JVS(103)*X(75))/(JVS(101))
  X(21) = (X(21)-JVS(100)*X(71))/(JVS(99))
  X(20) = (X(20)-JVS(98)*X(71))/(JVS(97))
  X(19) = (X(19)-JVS(95)*X(71)-JVS(96)*X(78))/(JVS(94))
  X(18) = (X(18)-JVS(92)*X(70)-JVS(93)*X(74))/(JVS(91))
  X(17) = (X(17)-JVS(89)*X(74)-JVS(90)*X(79))/(JVS(88))
  X(16) = (X(16)-JVS(86)*X(72)-JVS(87)*X(74))/(JVS(85))
  X(15) = (X(15)-JVS(83)*X(68)-JVS(84)*X(74))/(JVS(82))
  X(14) = (X(14)-JVS(78)*X(27)-JVS(79)*X(54)-JVS(80)*X(67)-JVS(81)*X(71))/(JVS(77))
  X(13) = (X(13)-JVS(76)*X(71))/(JVS(75))
  X(12) = (X(12)-JVS(74)*X(71))/(JVS(73))
  X(11) = (X(11)-JVS(72)*X(67))/(JVS(71))
  X(10) = (X(10)-JVS(70)*X(71))/(JVS(69))
  X(9) = (X(9)-JVS(61)*X(31)-JVS(62)*X(41)-JVS(63)*X(50)-JVS(64)*X(63)-JVS(65)*X(67)-JVS(66)*X(71)-JVS(67)*X(74)-JVS(68)&
           &*X(75))/(JVS(60))
  X(8) = (X(8)-JVS(56)*X(31)-JVS(57)*X(50)-JVS(58)*X(74)-JVS(59)*X(75))/(JVS(55))
  X(7) = (X(7)-JVS(51)*X(70)-JVS(52)*X(72)-JVS(53)*X(78)-JVS(54)*X(79))/(JVS(50))
  X(6) = (X(6)-JVS(48)*X(68)-JVS(49)*X(78))/(JVS(47))
  X(5) = (X(5)-JVS(34)*X(52)-JVS(35)*X(54)-JVS(36)*X(56)-JVS(37)*X(57)-JVS(38)*X(58)-JVS(39)*X(67)-JVS(40)*X(69)-JVS(41)&
           &*X(70)-JVS(42)*X(72)-JVS(43)*X(76)-JVS(44)*X(77)-JVS(45)*X(78)-JVS(46)*X(79))/(JVS(33))
  X(4) = (X(4)-JVS(24)*X(50)-JVS(25)*X(56)-JVS(26)*X(58)-JVS(27)*X(67)-JVS(28)*X(68)-JVS(29)*X(69)-JVS(30)*X(76)-JVS(31)&
           &*X(77)-JVS(32)*X(78))/(JVS(23))
  X(3) = (X(3)-JVS(9)*X(30)-JVS(10)*X(39)-JVS(11)*X(48)-JVS(12)*X(50)-JVS(13)*X(52)-JVS(14)*X(54)-JVS(15)*X(55)-JVS(16)&
           &*X(56)-JVS(17)*X(57)-JVS(18)*X(58)-JVS(19)*X(59)-JVS(20)*X(67)-JVS(21)*X(71)-JVS(22)*X(73))/(JVS(8))
  X(2) = (X(2)-JVS(5)*X(39)-JVS(6)*X(50)-JVS(7)*X(67))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(10)-JVS(3)*X(71))/(JVS(1))
      
END SUBROUTINE saprc99_KppSolve
























      SUBROUTINE saprc99_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE saprc99_WCOPY



      SUBROUTINE saprc99_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE saprc99_WAXPY




      SUBROUTINE saprc99_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE saprc99_WSCAL


      REAL(kind=dp) FUNCTION saprc99_WLAMCH( C )








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
          CALL saprc99_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      saprc99_WLAMCH = Eps

      END FUNCTION saprc99_WLAMCH
     
      SUBROUTINE saprc99_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE saprc99_WLAMCH_ADD




      SUBROUTINE saprc99_SET2ZERO(N,Y)




      
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

      END SUBROUTINE saprc99_SET2ZERO



      REAL(kind=dp) FUNCTION saprc99_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      saprc99_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        saprc99_WDOT = saprc99_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         saprc99_WDOT = saprc99_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          saprc99_WDOT = saprc99_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        saprc99_WDOT = saprc99_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION saprc99_WDOT                                          




   SUBROUTINE decomp_saprc99( JVS, IER )
   
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
   W( 10 ) = JVS( 2 )
   W( 71 ) = JVS( 3 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 10 )
  JVS( 3) = W( 71 )
  IF ( ABS(  JVS( 4 )) < TINY(a) ) THEN
         IER = 2                                       
         RETURN
  END IF
   W( 2 ) = JVS( 4 )
   W( 39 ) = JVS( 5 )
   W( 50 ) = JVS( 6 )
   W( 67 ) = JVS( 7 )
  JVS( 4) = W( 2 )
  JVS( 5) = W( 39 )
  JVS( 6) = W( 50 )
  JVS( 7) = W( 67 )
  IF ( ABS(  JVS( 8 )) < TINY(a) ) THEN
         IER = 3                                       
         RETURN
  END IF
   W( 3 ) = JVS( 8 )
   W( 30 ) = JVS( 9 )
   W( 39 ) = JVS( 10 )
   W( 48 ) = JVS( 11 )
   W( 50 ) = JVS( 12 )
   W( 52 ) = JVS( 13 )
   W( 54 ) = JVS( 14 )
   W( 55 ) = JVS( 15 )
   W( 56 ) = JVS( 16 )
   W( 57 ) = JVS( 17 )
   W( 58 ) = JVS( 18 )
   W( 59 ) = JVS( 19 )
   W( 67 ) = JVS( 20 )
   W( 71 ) = JVS( 21 )
   W( 73 ) = JVS( 22 )
  JVS( 8) = W( 3 )
  JVS( 9) = W( 30 )
  JVS( 10) = W( 39 )
  JVS( 11) = W( 48 )
  JVS( 12) = W( 50 )
  JVS( 13) = W( 52 )
  JVS( 14) = W( 54 )
  JVS( 15) = W( 55 )
  JVS( 16) = W( 56 )
  JVS( 17) = W( 57 )
  JVS( 18) = W( 58 )
  JVS( 19) = W( 59 )
  JVS( 20) = W( 67 )
  JVS( 21) = W( 71 )
  JVS( 22) = W( 73 )
  IF ( ABS(  JVS( 23 )) < TINY(a) ) THEN
         IER = 4                                       
         RETURN
  END IF
   W( 4 ) = JVS( 23 )
   W( 50 ) = JVS( 24 )
   W( 56 ) = JVS( 25 )
   W( 58 ) = JVS( 26 )
   W( 67 ) = JVS( 27 )
   W( 68 ) = JVS( 28 )
   W( 69 ) = JVS( 29 )
   W( 76 ) = JVS( 30 )
   W( 77 ) = JVS( 31 )
   W( 78 ) = JVS( 32 )
  JVS( 23) = W( 4 )
  JVS( 24) = W( 50 )
  JVS( 25) = W( 56 )
  JVS( 26) = W( 58 )
  JVS( 27) = W( 67 )
  JVS( 28) = W( 68 )
  JVS( 29) = W( 69 )
  JVS( 30) = W( 76 )
  JVS( 31) = W( 77 )
  JVS( 32) = W( 78 )
  IF ( ABS(  JVS( 33 )) < TINY(a) ) THEN
         IER = 5                                       
         RETURN
  END IF
   W( 5 ) = JVS( 33 )
   W( 52 ) = JVS( 34 )
   W( 54 ) = JVS( 35 )
   W( 56 ) = JVS( 36 )
   W( 57 ) = JVS( 37 )
   W( 58 ) = JVS( 38 )
   W( 67 ) = JVS( 39 )
   W( 69 ) = JVS( 40 )
   W( 70 ) = JVS( 41 )
   W( 72 ) = JVS( 42 )
   W( 76 ) = JVS( 43 )
   W( 77 ) = JVS( 44 )
   W( 78 ) = JVS( 45 )
   W( 79 ) = JVS( 46 )
  JVS( 33) = W( 5 )
  JVS( 34) = W( 52 )
  JVS( 35) = W( 54 )
  JVS( 36) = W( 56 )
  JVS( 37) = W( 57 )
  JVS( 38) = W( 58 )
  JVS( 39) = W( 67 )
  JVS( 40) = W( 69 )
  JVS( 41) = W( 70 )
  JVS( 42) = W( 72 )
  JVS( 43) = W( 76 )
  JVS( 44) = W( 77 )
  JVS( 45) = W( 78 )
  JVS( 46) = W( 79 )
  IF ( ABS(  JVS( 47 )) < TINY(a) ) THEN
         IER = 6                                       
         RETURN
  END IF
   W( 6 ) = JVS( 47 )
   W( 68 ) = JVS( 48 )
   W( 78 ) = JVS( 49 )
  JVS( 47) = W( 6 )
  JVS( 48) = W( 68 )
  JVS( 49) = W( 78 )
  IF ( ABS(  JVS( 50 )) < TINY(a) ) THEN
         IER = 7                                       
         RETURN
  END IF
   W( 7 ) = JVS( 50 )
   W( 70 ) = JVS( 51 )
   W( 72 ) = JVS( 52 )
   W( 78 ) = JVS( 53 )
   W( 79 ) = JVS( 54 )
  JVS( 50) = W( 7 )
  JVS( 51) = W( 70 )
  JVS( 52) = W( 72 )
  JVS( 53) = W( 78 )
  JVS( 54) = W( 79 )
  IF ( ABS(  JVS( 55 )) < TINY(a) ) THEN
         IER = 8                                       
         RETURN
  END IF
   W( 8 ) = JVS( 55 )
   W( 31 ) = JVS( 56 )
   W( 50 ) = JVS( 57 )
   W( 74 ) = JVS( 58 )
   W( 75 ) = JVS( 59 )
  JVS( 55) = W( 8 )
  JVS( 56) = W( 31 )
  JVS( 57) = W( 50 )
  JVS( 58) = W( 74 )
  JVS( 59) = W( 75 )
  IF ( ABS(  JVS( 60 )) < TINY(a) ) THEN
         IER = 9                                       
         RETURN
  END IF
   W( 9 ) = JVS( 60 )
   W( 31 ) = JVS( 61 )
   W( 41 ) = JVS( 62 )
   W( 50 ) = JVS( 63 )
   W( 63 ) = JVS( 64 )
   W( 67 ) = JVS( 65 )
   W( 71 ) = JVS( 66 )
   W( 74 ) = JVS( 67 )
   W( 75 ) = JVS( 68 )
  JVS( 60) = W( 9 )
  JVS( 61) = W( 31 )
  JVS( 62) = W( 41 )
  JVS( 63) = W( 50 )
  JVS( 64) = W( 63 )
  JVS( 65) = W( 67 )
  JVS( 66) = W( 71 )
  JVS( 67) = W( 74 )
  JVS( 68) = W( 75 )
  IF ( ABS(  JVS( 69 )) < TINY(a) ) THEN
         IER = 10                                      
         RETURN
  END IF
   W( 10 ) = JVS( 69 )
   W( 71 ) = JVS( 70 )
  JVS( 69) = W( 10 )
  JVS( 70) = W( 71 )
  IF ( ABS(  JVS( 71 )) < TINY(a) ) THEN
         IER = 11                                      
         RETURN
  END IF
   W( 11 ) = JVS( 71 )
   W( 67 ) = JVS( 72 )
  JVS( 71) = W( 11 )
  JVS( 72) = W( 67 )
  IF ( ABS(  JVS( 73 )) < TINY(a) ) THEN
         IER = 12                                      
         RETURN
  END IF
   W( 12 ) = JVS( 73 )
   W( 71 ) = JVS( 74 )
  JVS( 73) = W( 12 )
  JVS( 74) = W( 71 )
  IF ( ABS(  JVS( 75 )) < TINY(a) ) THEN
         IER = 13                                      
         RETURN
  END IF
   W( 13 ) = JVS( 75 )
   W( 71 ) = JVS( 76 )
  JVS( 75) = W( 13 )
  JVS( 76) = W( 71 )
  IF ( ABS(  JVS( 77 )) < TINY(a) ) THEN
         IER = 14                                      
         RETURN
  END IF
   W( 14 ) = JVS( 77 )
   W( 27 ) = JVS( 78 )
   W( 54 ) = JVS( 79 )
   W( 67 ) = JVS( 80 )
   W( 71 ) = JVS( 81 )
  JVS( 77) = W( 14 )
  JVS( 78) = W( 27 )
  JVS( 79) = W( 54 )
  JVS( 80) = W( 67 )
  JVS( 81) = W( 71 )
  IF ( ABS(  JVS( 82 )) < TINY(a) ) THEN
         IER = 15                                      
         RETURN
  END IF
   W( 15 ) = JVS( 82 )
   W( 68 ) = JVS( 83 )
   W( 74 ) = JVS( 84 )
  JVS( 82) = W( 15 )
  JVS( 83) = W( 68 )
  JVS( 84) = W( 74 )
  IF ( ABS(  JVS( 85 )) < TINY(a) ) THEN
         IER = 16                                      
         RETURN
  END IF
   W( 16 ) = JVS( 85 )
   W( 72 ) = JVS( 86 )
   W( 74 ) = JVS( 87 )
  JVS( 85) = W( 16 )
  JVS( 86) = W( 72 )
  JVS( 87) = W( 74 )
  IF ( ABS(  JVS( 88 )) < TINY(a) ) THEN
         IER = 17                                      
         RETURN
  END IF
   W( 17 ) = JVS( 88 )
   W( 74 ) = JVS( 89 )
   W( 79 ) = JVS( 90 )
  JVS( 88) = W( 17 )
  JVS( 89) = W( 74 )
  JVS( 90) = W( 79 )
  IF ( ABS(  JVS( 91 )) < TINY(a) ) THEN
         IER = 18                                      
         RETURN
  END IF
   W( 18 ) = JVS( 91 )
   W( 70 ) = JVS( 92 )
   W( 74 ) = JVS( 93 )
  JVS( 91) = W( 18 )
  JVS( 92) = W( 70 )
  JVS( 93) = W( 74 )
  IF ( ABS(  JVS( 94 )) < TINY(a) ) THEN
         IER = 19                                      
         RETURN
  END IF
   W( 19 ) = JVS( 94 )
   W( 71 ) = JVS( 95 )
   W( 78 ) = JVS( 96 )
  JVS( 94) = W( 19 )
  JVS( 95) = W( 71 )
  JVS( 96) = W( 78 )
  IF ( ABS(  JVS( 97 )) < TINY(a) ) THEN
         IER = 20                                      
         RETURN
  END IF
   W( 20 ) = JVS( 97 )
   W( 71 ) = JVS( 98 )
  JVS( 97) = W( 20 )
  JVS( 98) = W( 71 )
  IF ( ABS(  JVS( 99 )) < TINY(a) ) THEN
         IER = 21                                      
         RETURN
  END IF
   W( 21 ) = JVS( 99 )
   W( 71 ) = JVS( 100 )
  JVS( 99) = W( 21 )
  JVS( 100) = W( 71 )
  IF ( ABS(  JVS( 101 )) < TINY(a) ) THEN
         IER = 22                                      
         RETURN
  END IF
   W( 22 ) = JVS( 101 )
   W( 74 ) = JVS( 102 )
   W( 75 ) = JVS( 103 )
  JVS( 101) = W( 22 )
  JVS( 102) = W( 74 )
  JVS( 103) = W( 75 )
  IF ( ABS(  JVS( 104 )) < TINY(a) ) THEN
         IER = 23                                      
         RETURN
  END IF
   W( 23 ) = JVS( 104 )
   W( 71 ) = JVS( 105 )
   W( 73 ) = JVS( 106 )
  JVS( 104) = W( 23 )
  JVS( 105) = W( 71 )
  JVS( 106) = W( 73 )
  IF ( ABS(  JVS( 107 )) < TINY(a) ) THEN
         IER = 24                                      
         RETURN
  END IF
   W( 24 ) = JVS( 107 )
   W( 71 ) = JVS( 108 )
  JVS( 107) = W( 24 )
  JVS( 108) = W( 71 )
  IF ( ABS(  JVS( 110 )) < TINY(a) ) THEN
         IER = 25                                      
         RETURN
  END IF
   W( 24 ) = JVS( 109 )
   W( 25 ) = JVS( 110 )
   W( 71 ) = JVS( 111 )
   W( 74 ) = JVS( 112 )
  a = -W( 24 ) / JVS(          107  )
  W( 24 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 108 )
  JVS( 109) = W( 24 )
  JVS( 110) = W( 25 )
  JVS( 111) = W( 71 )
  JVS( 112) = W( 74 )
  IF ( ABS(  JVS( 113 )) < TINY(a) ) THEN
         IER = 26                                      
         RETURN
  END IF
   W( 26 ) = JVS( 113 )
   W( 71 ) = JVS( 114 )
  JVS( 113) = W( 26 )
  JVS( 114) = W( 71 )
  IF ( ABS(  JVS( 115 )) < TINY(a) ) THEN
         IER = 27                                      
         RETURN
  END IF
   W( 27 ) = JVS( 115 )
   W( 71 ) = JVS( 116 )
  JVS( 115) = W( 27 )
  JVS( 116) = W( 71 )
  IF ( ABS(  JVS( 117 )) < TINY(a) ) THEN
         IER = 28                                      
         RETURN
  END IF
   W( 28 ) = JVS( 117 )
   W( 69 ) = JVS( 118 )
   W( 71 ) = JVS( 119 )
   W( 76 ) = JVS( 120 )
   W( 77 ) = JVS( 121 )
  JVS( 117) = W( 28 )
  JVS( 118) = W( 69 )
  JVS( 119) = W( 71 )
  JVS( 120) = W( 76 )
  JVS( 121) = W( 77 )
  IF ( ABS(  JVS( 122 )) < TINY(a) ) THEN
         IER = 29                                      
         RETURN
  END IF
   W( 29 ) = JVS( 122 )
   W( 71 ) = JVS( 123 )
   W( 76 ) = JVS( 124 )
   W( 78 ) = JVS( 125 )
  JVS( 122) = W( 29 )
  JVS( 123) = W( 71 )
  JVS( 124) = W( 76 )
  JVS( 125) = W( 78 )
  IF ( ABS(  JVS( 126 )) < TINY(a) ) THEN
         IER = 30                                      
         RETURN
  END IF
   W( 30 ) = JVS( 126 )
   W( 61 ) = JVS( 127 )
   W( 73 ) = JVS( 128 )
   W( 78 ) = JVS( 129 )
  JVS( 126) = W( 30 )
  JVS( 127) = W( 61 )
  JVS( 128) = W( 73 )
  JVS( 129) = W( 78 )
  IF ( ABS(  JVS( 130 )) < TINY(a) ) THEN
         IER = 31                                      
         RETURN
  END IF
   W( 31 ) = JVS( 130 )
   W( 40 ) = JVS( 131 )
   W( 74 ) = JVS( 132 )
   W( 75 ) = JVS( 133 )
   W( 78 ) = JVS( 134 )
  JVS( 130) = W( 31 )
  JVS( 131) = W( 40 )
  JVS( 132) = W( 74 )
  JVS( 133) = W( 75 )
  JVS( 134) = W( 78 )
  IF ( ABS(  JVS( 135 )) < TINY(a) ) THEN
         IER = 32                                      
         RETURN
  END IF
   W( 32 ) = JVS( 135 )
   W( 71 ) = JVS( 136 )
   W( 74 ) = JVS( 137 )
   W( 78 ) = JVS( 138 )
  JVS( 135) = W( 32 )
  JVS( 136) = W( 71 )
  JVS( 137) = W( 74 )
  JVS( 138) = W( 78 )
  IF ( ABS(  JVS( 139 )) < TINY(a) ) THEN
         IER = 33                                      
         RETURN
  END IF
   W( 33 ) = JVS( 139 )
   W( 71 ) = JVS( 140 )
  JVS( 139) = W( 33 )
  JVS( 140) = W( 71 )
  IF ( ABS(  JVS( 141 )) < TINY(a) ) THEN
         IER = 34                                      
         RETURN
  END IF
   W( 34 ) = JVS( 141 )
   W( 71 ) = JVS( 142 )
  JVS( 141) = W( 34 )
  JVS( 142) = W( 71 )
  IF ( ABS(  JVS( 145 )) < TINY(a) ) THEN
         IER = 35                                      
         RETURN
  END IF
   W( 27 ) = JVS( 143 )
   W( 34 ) = JVS( 144 )
   W( 35 ) = JVS( 145 )
   W( 71 ) = JVS( 146 )
  a = -W( 27 ) / JVS(          115  )
  W( 27 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 116 )
  a = -W( 34 ) / JVS(          141  )
  W( 34 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 142 )
  JVS( 143) = W( 27 )
  JVS( 144) = W( 34 )
  JVS( 145) = W( 35 )
  JVS( 146) = W( 71 )
  IF ( ABS(  JVS( 149 )) < TINY(a) ) THEN
         IER = 36                                      
         RETURN
  END IF
   W( 27 ) = JVS( 147 )
   W( 34 ) = JVS( 148 )
   W( 36 ) = JVS( 149 )
   W( 71 ) = JVS( 150 )
  a = -W( 27 ) / JVS(          115  )
  W( 27 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 116 )
  a = -W( 34 ) / JVS(          141  )
  W( 34 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 142 )
  JVS( 147) = W( 27 )
  JVS( 148) = W( 34 )
  JVS( 149) = W( 36 )
  JVS( 150) = W( 71 )
  IF ( ABS(  JVS( 153 )) < TINY(a) ) THEN
         IER = 37                                      
         RETURN
  END IF
   W( 27 ) = JVS( 151 )
   W( 34 ) = JVS( 152 )
   W( 37 ) = JVS( 153 )
   W( 71 ) = JVS( 154 )
   W( 75 ) = JVS( 155 )
  a = -W( 27 ) / JVS(          115  )
  W( 27 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 116 )
  a = -W( 34 ) / JVS(          141  )
  W( 34 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 142 )
  JVS( 151) = W( 27 )
  JVS( 152) = W( 34 )
  JVS( 153) = W( 37 )
  JVS( 154) = W( 71 )
  JVS( 155) = W( 75 )
  IF ( ABS(  JVS( 158 )) < TINY(a) ) THEN
         IER = 38                                      
         RETURN
  END IF
   W( 27 ) = JVS( 156 )
   W( 34 ) = JVS( 157 )
   W( 38 ) = JVS( 158 )
   W( 67 ) = JVS( 159 )
   W( 71 ) = JVS( 160 )
  a = -W( 27 ) / JVS(          115  )
  W( 27 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 116 )
  a = -W( 34 ) / JVS(          141  )
  W( 34 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 142 )
  JVS( 156) = W( 27 )
  JVS( 157) = W( 34 )
  JVS( 158) = W( 38 )
  JVS( 159) = W( 67 )
  JVS( 160) = W( 71 )
  IF ( ABS(  JVS( 161 )) < TINY(a) ) THEN
         IER = 39                                      
         RETURN
  END IF
   W( 39 ) = JVS( 161 )
   W( 67 ) = JVS( 162 )
   W( 71 ) = JVS( 163 )
  JVS( 161) = W( 39 )
  JVS( 162) = W( 67 )
  JVS( 163) = W( 71 )
  IF ( ABS(  JVS( 165 )) < TINY(a) ) THEN
         IER = 40                                      
         RETURN
  END IF
   W( 31 ) = JVS( 164 )
   W( 40 ) = JVS( 165 )
   W( 51 ) = JVS( 166 )
   W( 74 ) = JVS( 167 )
   W( 75 ) = JVS( 168 )
   W( 78 ) = JVS( 169 )
  a = -W( 31 ) / JVS(          130  )
  W( 31 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 131 )
  W( 74 ) = W( 74 ) + a*JVS( 132 )
  W( 75 ) = W( 75 ) + a*JVS( 133 )
  W( 78 ) = W( 78 ) + a*JVS( 134 )
  JVS( 164) = W( 31 )
  JVS( 165) = W( 40 )
  JVS( 166) = W( 51 )
  JVS( 167) = W( 74 )
  JVS( 168) = W( 75 )
  JVS( 169) = W( 78 )
  IF ( ABS(  JVS( 172 )) < TINY(a) ) THEN
         IER = 41                                      
         RETURN
  END IF
   W( 27 ) = JVS( 170 )
   W( 34 ) = JVS( 171 )
   W( 41 ) = JVS( 172 )
   W( 58 ) = JVS( 173 )
   W( 67 ) = JVS( 174 )
   W( 71 ) = JVS( 175 )
   W( 75 ) = JVS( 176 )
  a = -W( 27 ) / JVS(          115  )
  W( 27 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 116 )
  a = -W( 34 ) / JVS(          141  )
  W( 34 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 142 )
  JVS( 170) = W( 27 )
  JVS( 171) = W( 34 )
  JVS( 172) = W( 41 )
  JVS( 173) = W( 58 )
  JVS( 174) = W( 67 )
  JVS( 175) = W( 71 )
  JVS( 176) = W( 75 )
  IF ( ABS(  JVS( 177 )) < TINY(a) ) THEN
         IER = 42                                      
         RETURN
  END IF
   W( 42 ) = JVS( 177 )
   W( 69 ) = JVS( 178 )
   W( 71 ) = JVS( 179 )
   W( 77 ) = JVS( 180 )
   W( 78 ) = JVS( 181 )
  JVS( 177) = W( 42 )
  JVS( 178) = W( 69 )
  JVS( 179) = W( 71 )
  JVS( 180) = W( 77 )
  JVS( 181) = W( 78 )
  IF ( ABS(  JVS( 183 )) < TINY(a) ) THEN
         IER = 43                                      
         RETURN
  END IF
   W( 34 ) = JVS( 182 )
   W( 43 ) = JVS( 183 )
   W( 51 ) = JVS( 184 )
   W( 71 ) = JVS( 185 )
   W( 75 ) = JVS( 186 )
   W( 78 ) = JVS( 187 )
  a = -W( 34 ) / JVS(          141  )
  W( 34 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 142 )
  JVS( 182) = W( 34 )
  JVS( 183) = W( 43 )
  JVS( 184) = W( 51 )
  JVS( 185) = W( 71 )
  JVS( 186) = W( 75 )
  JVS( 187) = W( 78 )
  IF ( ABS(  JVS( 193 )) < TINY(a) ) THEN
         IER = 44                                      
         RETURN
  END IF
   W( 27 ) = JVS( 188 )
   W( 34 ) = JVS( 189 )
   W( 35 ) = JVS( 190 )
   W( 36 ) = JVS( 191 )
   W( 37 ) = JVS( 192 )
   W( 44 ) = JVS( 193 )
   W( 55 ) = JVS( 194 )
   W( 57 ) = JVS( 195 )
   W( 59 ) = JVS( 196 )
   W( 67 ) = JVS( 197 )
   W( 71 ) = JVS( 198 )
   W( 75 ) = JVS( 199 )
  a = -W( 27 ) / JVS(          115  )
  W( 27 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 116 )
  a = -W( 34 ) / JVS(          141  )
  W( 34 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 142 )
  a = -W( 35 ) / JVS(          145  )
  W( 35 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 146 )
  a = -W( 36 ) / JVS(          149  )
  W( 36 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 150 )
  a = -W( 37 ) / JVS(          153  )
  W( 37 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 154 )
  W( 75 ) = W( 75 ) + a*JVS( 155 )
  JVS( 188) = W( 27 )
  JVS( 189) = W( 34 )
  JVS( 190) = W( 35 )
  JVS( 191) = W( 36 )
  JVS( 192) = W( 37 )
  JVS( 193) = W( 44 )
  JVS( 194) = W( 55 )
  JVS( 195) = W( 57 )
  JVS( 196) = W( 59 )
  JVS( 197) = W( 67 )
  JVS( 198) = W( 71 )
  JVS( 199) = W( 75 )
  IF ( ABS(  JVS( 206 )) < TINY(a) ) THEN
         IER = 45                                      
         RETURN
  END IF
   W( 33 ) = JVS( 200 )
   W( 35 ) = JVS( 201 )
   W( 36 ) = JVS( 202 )
   W( 38 ) = JVS( 203 )
   W( 39 ) = JVS( 204 )
   W( 44 ) = JVS( 205 )
   W( 45 ) = JVS( 206 )
   W( 48 ) = JVS( 207 )
   W( 49 ) = JVS( 208 )
   W( 50 ) = JVS( 209 )
   W( 52 ) = JVS( 210 )
   W( 54 ) = JVS( 211 )
   W( 55 ) = JVS( 212 )
   W( 56 ) = JVS( 213 )
   W( 57 ) = JVS( 214 )
   W( 58 ) = JVS( 215 )
   W( 59 ) = JVS( 216 )
   W( 60 ) = JVS( 217 )
   W( 61 ) = JVS( 218 )
   W( 63 ) = JVS( 219 )
   W( 64 ) = JVS( 220 )
   W( 67 ) = JVS( 221 )
   W( 71 ) = JVS( 222 )
   W( 75 ) = JVS( 223 )
  a = -W( 33 ) / JVS(          139  )
  W( 33 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 140 )
  a = -W( 35 ) / JVS(          145  )
  W( 35 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 146 )
  a = -W( 36 ) / JVS(          149  )
  W( 36 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 150 )
  a = -W( 38 ) / JVS(          158  )
  W( 38 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 159 )
  W( 71 ) = W( 71 ) + a*JVS( 160 )
  a = -W( 39 ) / JVS(          161  )
  W( 39 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 162 )
  W( 71 ) = W( 71 ) + a*JVS( 163 )
  a = -W( 44 ) / JVS(          193  )
  W( 44 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  W( 57 ) = W( 57 ) + a*JVS( 195 )
  W( 59 ) = W( 59 ) + a*JVS( 196 )
  W( 67 ) = W( 67 ) + a*JVS( 197 )
  W( 71 ) = W( 71 ) + a*JVS( 198 )
  W( 75 ) = W( 75 ) + a*JVS( 199 )
  JVS( 200) = W( 33 )
  JVS( 201) = W( 35 )
  JVS( 202) = W( 36 )
  JVS( 203) = W( 38 )
  JVS( 204) = W( 39 )
  JVS( 205) = W( 44 )
  JVS( 206) = W( 45 )
  JVS( 207) = W( 48 )
  JVS( 208) = W( 49 )
  JVS( 209) = W( 50 )
  JVS( 210) = W( 52 )
  JVS( 211) = W( 54 )
  JVS( 212) = W( 55 )
  JVS( 213) = W( 56 )
  JVS( 214) = W( 57 )
  JVS( 215) = W( 58 )
  JVS( 216) = W( 59 )
  JVS( 217) = W( 60 )
  JVS( 218) = W( 61 )
  JVS( 219) = W( 63 )
  JVS( 220) = W( 64 )
  JVS( 221) = W( 67 )
  JVS( 222) = W( 71 )
  JVS( 223) = W( 75 )
  IF ( ABS(  JVS( 229 )) < TINY(a) ) THEN
         IER = 46                                      
         RETURN
  END IF
   W( 20 ) = JVS( 224 )
   W( 24 ) = JVS( 225 )
   W( 25 ) = JVS( 226 )
   W( 26 ) = JVS( 227 )
   W( 33 ) = JVS( 228 )
   W( 46 ) = JVS( 229 )
   W( 54 ) = JVS( 230 )
   W( 56 ) = JVS( 231 )
   W( 58 ) = JVS( 232 )
   W( 62 ) = JVS( 233 )
   W( 67 ) = JVS( 234 )
   W( 71 ) = JVS( 235 )
   W( 74 ) = JVS( 236 )
   W( 75 ) = JVS( 237 )
  a = -W( 20 ) / JVS(           97  )
  W( 20 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 98 )
  a = -W( 24 ) / JVS(          107  )
  W( 24 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 108 )
  a = -W( 25 ) / JVS(          110  )
  W( 25 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 111 )
  W( 74 ) = W( 74 ) + a*JVS( 112 )
  a = -W( 26 ) / JVS(          113  )
  W( 26 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 114 )
  a = -W( 33 ) / JVS(          139  )
  W( 33 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 140 )
  JVS( 224) = W( 20 )
  JVS( 225) = W( 24 )
  JVS( 226) = W( 25 )
  JVS( 227) = W( 26 )
  JVS( 228) = W( 33 )
  JVS( 229) = W( 46 )
  JVS( 230) = W( 54 )
  JVS( 231) = W( 56 )
  JVS( 232) = W( 58 )
  JVS( 233) = W( 62 )
  JVS( 234) = W( 67 )
  JVS( 235) = W( 71 )
  JVS( 236) = W( 74 )
  JVS( 237) = W( 75 )
  IF ( ABS(  JVS( 244 )) < TINY(a) ) THEN
         IER = 47                                      
         RETURN
  END IF
   W( 22 ) = JVS( 238 )
   W( 37 ) = JVS( 239 )
   W( 40 ) = JVS( 240 )
   W( 41 ) = JVS( 241 )
   W( 43 ) = JVS( 242 )
   W( 44 ) = JVS( 243 )
   W( 47 ) = JVS( 244 )
   W( 49 ) = JVS( 245 )
   W( 51 ) = JVS( 246 )
   W( 55 ) = JVS( 247 )
   W( 57 ) = JVS( 248 )
   W( 58 ) = JVS( 249 )
   W( 59 ) = JVS( 250 )
   W( 60 ) = JVS( 251 )
   W( 61 ) = JVS( 252 )
   W( 64 ) = JVS( 253 )
   W( 67 ) = JVS( 254 )
   W( 71 ) = JVS( 255 )
   W( 74 ) = JVS( 256 )
   W( 75 ) = JVS( 257 )
   W( 78 ) = JVS( 258 )
  a = -W( 22 ) / JVS(          101  )
  W( 22 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 102 )
  W( 75 ) = W( 75 ) + a*JVS( 103 )
  a = -W( 37 ) / JVS(          153  )
  W( 37 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 154 )
  W( 75 ) = W( 75 ) + a*JVS( 155 )
  a = -W( 40 ) / JVS(          165  )
  W( 40 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 166 )
  W( 74 ) = W( 74 ) + a*JVS( 167 )
  W( 75 ) = W( 75 ) + a*JVS( 168 )
  W( 78 ) = W( 78 ) + a*JVS( 169 )
  a = -W( 41 ) / JVS(          172  )
  W( 41 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 173 )
  W( 67 ) = W( 67 ) + a*JVS( 174 )
  W( 71 ) = W( 71 ) + a*JVS( 175 )
  W( 75 ) = W( 75 ) + a*JVS( 176 )
  a = -W( 43 ) / JVS(          183  )
  W( 43 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 184 )
  W( 71 ) = W( 71 ) + a*JVS( 185 )
  W( 75 ) = W( 75 ) + a*JVS( 186 )
  W( 78 ) = W( 78 ) + a*JVS( 187 )
  a = -W( 44 ) / JVS(          193  )
  W( 44 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  W( 57 ) = W( 57 ) + a*JVS( 195 )
  W( 59 ) = W( 59 ) + a*JVS( 196 )
  W( 67 ) = W( 67 ) + a*JVS( 197 )
  W( 71 ) = W( 71 ) + a*JVS( 198 )
  W( 75 ) = W( 75 ) + a*JVS( 199 )
  JVS( 238) = W( 22 )
  JVS( 239) = W( 37 )
  JVS( 240) = W( 40 )
  JVS( 241) = W( 41 )
  JVS( 242) = W( 43 )
  JVS( 243) = W( 44 )
  JVS( 244) = W( 47 )
  JVS( 245) = W( 49 )
  JVS( 246) = W( 51 )
  JVS( 247) = W( 55 )
  JVS( 248) = W( 57 )
  JVS( 249) = W( 58 )
  JVS( 250) = W( 59 )
  JVS( 251) = W( 60 )
  JVS( 252) = W( 61 )
  JVS( 253) = W( 64 )
  JVS( 254) = W( 67 )
  JVS( 255) = W( 71 )
  JVS( 256) = W( 74 )
  JVS( 257) = W( 75 )
  JVS( 258) = W( 78 )
  IF ( ABS(  JVS( 259 )) < TINY(a) ) THEN
         IER = 48                                      
         RETURN
  END IF
   W( 48 ) = JVS( 259 )
   W( 63 ) = JVS( 260 )
   W( 67 ) = JVS( 261 )
   W( 71 ) = JVS( 262 )
   W( 75 ) = JVS( 263 )
  JVS( 259) = W( 48 )
  JVS( 260) = W( 63 )
  JVS( 261) = W( 67 )
  JVS( 262) = W( 71 )
  JVS( 263) = W( 75 )
  IF ( ABS(  JVS( 272 )) < TINY(a) ) THEN
         IER = 49                                      
         RETURN
  END IF
   W( 27 ) = JVS( 264 )
   W( 34 ) = JVS( 265 )
   W( 35 ) = JVS( 266 )
   W( 36 ) = JVS( 267 )
   W( 38 ) = JVS( 268 )
   W( 39 ) = JVS( 269 )
   W( 43 ) = JVS( 270 )
   W( 48 ) = JVS( 271 )
   W( 49 ) = JVS( 272 )
   W( 51 ) = JVS( 273 )
   W( 54 ) = JVS( 274 )
   W( 57 ) = JVS( 275 )
   W( 63 ) = JVS( 276 )
   W( 67 ) = JVS( 277 )
   W( 71 ) = JVS( 278 )
   W( 75 ) = JVS( 279 )
   W( 78 ) = JVS( 280 )
  a = -W( 27 ) / JVS(          115  )
  W( 27 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 116 )
  a = -W( 34 ) / JVS(          141  )
  W( 34 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 142 )
  a = -W( 35 ) / JVS(          145  )
  W( 35 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 146 )
  a = -W( 36 ) / JVS(          149  )
  W( 36 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 150 )
  a = -W( 38 ) / JVS(          158  )
  W( 38 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 159 )
  W( 71 ) = W( 71 ) + a*JVS( 160 )
  a = -W( 39 ) / JVS(          161  )
  W( 39 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 162 )
  W( 71 ) = W( 71 ) + a*JVS( 163 )
  a = -W( 43 ) / JVS(          183  )
  W( 43 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 184 )
  W( 71 ) = W( 71 ) + a*JVS( 185 )
  W( 75 ) = W( 75 ) + a*JVS( 186 )
  W( 78 ) = W( 78 ) + a*JVS( 187 )
  a = -W( 48 ) / JVS(          259  )
  W( 48 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 260 )
  W( 67 ) = W( 67 ) + a*JVS( 261 )
  W( 71 ) = W( 71 ) + a*JVS( 262 )
  W( 75 ) = W( 75 ) + a*JVS( 263 )
  JVS( 264) = W( 27 )
  JVS( 265) = W( 34 )
  JVS( 266) = W( 35 )
  JVS( 267) = W( 36 )
  JVS( 268) = W( 38 )
  JVS( 269) = W( 39 )
  JVS( 270) = W( 43 )
  JVS( 271) = W( 48 )
  JVS( 272) = W( 49 )
  JVS( 273) = W( 51 )
  JVS( 274) = W( 54 )
  JVS( 275) = W( 57 )
  JVS( 276) = W( 63 )
  JVS( 277) = W( 67 )
  JVS( 278) = W( 71 )
  JVS( 279) = W( 75 )
  JVS( 280) = W( 78 )
  IF ( ABS(  JVS( 281 )) < TINY(a) ) THEN
         IER = 50                                      
         RETURN
  END IF
   W( 50 ) = JVS( 281 )
   W( 63 ) = JVS( 282 )
   W( 67 ) = JVS( 283 )
   W( 71 ) = JVS( 284 )
   W( 75 ) = JVS( 285 )
  JVS( 281) = W( 50 )
  JVS( 282) = W( 63 )
  JVS( 283) = W( 67 )
  JVS( 284) = W( 71 )
  JVS( 285) = W( 75 )
  IF ( ABS(  JVS( 288 )) < TINY(a) ) THEN
         IER = 51                                      
         RETURN
  END IF
   W( 37 ) = JVS( 286 )
   W( 43 ) = JVS( 287 )
   W( 51 ) = JVS( 288 )
   W( 68 ) = JVS( 289 )
   W( 70 ) = JVS( 290 )
   W( 71 ) = JVS( 291 )
   W( 72 ) = JVS( 292 )
   W( 73 ) = JVS( 293 )
   W( 74 ) = JVS( 294 )
   W( 75 ) = JVS( 295 )
   W( 78 ) = JVS( 296 )
   W( 79 ) = JVS( 297 )
  a = -W( 37 ) / JVS(          153  )
  W( 37 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 154 )
  W( 75 ) = W( 75 ) + a*JVS( 155 )
  a = -W( 43 ) / JVS(          183  )
  W( 43 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 184 )
  W( 71 ) = W( 71 ) + a*JVS( 185 )
  W( 75 ) = W( 75 ) + a*JVS( 186 )
  W( 78 ) = W( 78 ) + a*JVS( 187 )
  JVS( 286) = W( 37 )
  JVS( 287) = W( 43 )
  JVS( 288) = W( 51 )
  JVS( 289) = W( 68 )
  JVS( 290) = W( 70 )
  JVS( 291) = W( 71 )
  JVS( 292) = W( 72 )
  JVS( 293) = W( 73 )
  JVS( 294) = W( 74 )
  JVS( 295) = W( 75 )
  JVS( 296) = W( 78 )
  JVS( 297) = W( 79 )
  IF ( ABS(  JVS( 298 )) < TINY(a) ) THEN
         IER = 52                                      
         RETURN
  END IF
   W( 52 ) = JVS( 298 )
   W( 63 ) = JVS( 299 )
   W( 67 ) = JVS( 300 )
   W( 71 ) = JVS( 301 )
   W( 75 ) = JVS( 302 )
  JVS( 298) = W( 52 )
  JVS( 299) = W( 63 )
  JVS( 300) = W( 67 )
  JVS( 301) = W( 71 )
  JVS( 302) = W( 75 )
  IF ( ABS(  JVS( 310 )) < TINY(a) ) THEN
         IER = 53                                      
         RETURN
  END IF
   W( 24 ) = JVS( 303 )
   W( 26 ) = JVS( 304 )
   W( 33 ) = JVS( 305 )
   W( 35 ) = JVS( 306 )
   W( 36 ) = JVS( 307 )
   W( 46 ) = JVS( 308 )
   W( 52 ) = JVS( 309 )
   W( 53 ) = JVS( 310 )
   W( 54 ) = JVS( 311 )
   W( 56 ) = JVS( 312 )
   W( 58 ) = JVS( 313 )
   W( 59 ) = JVS( 314 )
   W( 62 ) = JVS( 315 )
   W( 63 ) = JVS( 316 )
   W( 65 ) = JVS( 317 )
   W( 66 ) = JVS( 318 )
   W( 67 ) = JVS( 319 )
   W( 68 ) = JVS( 320 )
   W( 69 ) = JVS( 321 )
   W( 70 ) = JVS( 322 )
   W( 71 ) = JVS( 323 )
   W( 72 ) = JVS( 324 )
   W( 73 ) = JVS( 325 )
   W( 74 ) = JVS( 326 )
   W( 75 ) = JVS( 327 )
   W( 76 ) = JVS( 328 )
   W( 77 ) = JVS( 329 )
   W( 78 ) = JVS( 330 )
   W( 79 ) = JVS( 331 )
  a = -W( 24 ) / JVS(          107  )
  W( 24 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 108 )
  a = -W( 26 ) / JVS(          113  )
  W( 26 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 114 )
  a = -W( 33 ) / JVS(          139  )
  W( 33 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 140 )
  a = -W( 35 ) / JVS(          145  )
  W( 35 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 146 )
  a = -W( 36 ) / JVS(          149  )
  W( 36 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 150 )
  a = -W( 46 ) / JVS(          229  )
  W( 46 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 230 )
  W( 56 ) = W( 56 ) + a*JVS( 231 )
  W( 58 ) = W( 58 ) + a*JVS( 232 )
  W( 62 ) = W( 62 ) + a*JVS( 233 )
  W( 67 ) = W( 67 ) + a*JVS( 234 )
  W( 71 ) = W( 71 ) + a*JVS( 235 )
  W( 74 ) = W( 74 ) + a*JVS( 236 )
  W( 75 ) = W( 75 ) + a*JVS( 237 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  JVS( 303) = W( 24 )
  JVS( 304) = W( 26 )
  JVS( 305) = W( 33 )
  JVS( 306) = W( 35 )
  JVS( 307) = W( 36 )
  JVS( 308) = W( 46 )
  JVS( 309) = W( 52 )
  JVS( 310) = W( 53 )
  JVS( 311) = W( 54 )
  JVS( 312) = W( 56 )
  JVS( 313) = W( 58 )
  JVS( 314) = W( 59 )
  JVS( 315) = W( 62 )
  JVS( 316) = W( 63 )
  JVS( 317) = W( 65 )
  JVS( 318) = W( 66 )
  JVS( 319) = W( 67 )
  JVS( 320) = W( 68 )
  JVS( 321) = W( 69 )
  JVS( 322) = W( 70 )
  JVS( 323) = W( 71 )
  JVS( 324) = W( 72 )
  JVS( 325) = W( 73 )
  JVS( 326) = W( 74 )
  JVS( 327) = W( 75 )
  JVS( 328) = W( 76 )
  JVS( 329) = W( 77 )
  JVS( 330) = W( 78 )
  JVS( 331) = W( 79 )
  IF ( ABS(  JVS( 332 )) < TINY(a) ) THEN
         IER = 54                                      
         RETURN
  END IF
   W( 54 ) = JVS( 332 )
   W( 63 ) = JVS( 333 )
   W( 67 ) = JVS( 334 )
   W( 71 ) = JVS( 335 )
   W( 75 ) = JVS( 336 )
  JVS( 332) = W( 54 )
  JVS( 333) = W( 63 )
  JVS( 334) = W( 67 )
  JVS( 335) = W( 71 )
  JVS( 336) = W( 75 )
  IF ( ABS(  JVS( 338 )) < TINY(a) ) THEN
         IER = 55                                      
         RETURN
  END IF
   W( 52 ) = JVS( 337 )
   W( 55 ) = JVS( 338 )
   W( 58 ) = JVS( 339 )
   W( 63 ) = JVS( 340 )
   W( 67 ) = JVS( 341 )
   W( 71 ) = JVS( 342 )
   W( 75 ) = JVS( 343 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  JVS( 337) = W( 52 )
  JVS( 338) = W( 55 )
  JVS( 339) = W( 58 )
  JVS( 340) = W( 63 )
  JVS( 341) = W( 67 )
  JVS( 342) = W( 71 )
  JVS( 343) = W( 75 )
  IF ( ABS(  JVS( 344 )) < TINY(a) ) THEN
         IER = 56                                      
         RETURN
  END IF
   W( 56 ) = JVS( 344 )
   W( 63 ) = JVS( 345 )
   W( 67 ) = JVS( 346 )
   W( 71 ) = JVS( 347 )
   W( 75 ) = JVS( 348 )
  JVS( 344) = W( 56 )
  JVS( 345) = W( 63 )
  JVS( 346) = W( 67 )
  JVS( 347) = W( 71 )
  JVS( 348) = W( 75 )
  IF ( ABS(  JVS( 350 )) < TINY(a) ) THEN
         IER = 57                                      
         RETURN
  END IF
   W( 52 ) = JVS( 349 )
   W( 57 ) = JVS( 350 )
   W( 58 ) = JVS( 351 )
   W( 63 ) = JVS( 352 )
   W( 67 ) = JVS( 353 )
   W( 71 ) = JVS( 354 )
   W( 75 ) = JVS( 355 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  JVS( 349) = W( 52 )
  JVS( 350) = W( 57 )
  JVS( 351) = W( 58 )
  JVS( 352) = W( 63 )
  JVS( 353) = W( 67 )
  JVS( 354) = W( 71 )
  JVS( 355) = W( 75 )
  IF ( ABS(  JVS( 356 )) < TINY(a) ) THEN
         IER = 58                                      
         RETURN
  END IF
   W( 58 ) = JVS( 356 )
   W( 63 ) = JVS( 357 )
   W( 67 ) = JVS( 358 )
   W( 71 ) = JVS( 359 )
   W( 75 ) = JVS( 360 )
  JVS( 356) = W( 58 )
  JVS( 357) = W( 63 )
  JVS( 358) = W( 67 )
  JVS( 359) = W( 71 )
  JVS( 360) = W( 75 )
  IF ( ABS(  JVS( 363 )) < TINY(a) ) THEN
         IER = 59                                      
         RETURN
  END IF
   W( 52 ) = JVS( 361 )
   W( 58 ) = JVS( 362 )
   W( 59 ) = JVS( 363 )
   W( 63 ) = JVS( 364 )
   W( 67 ) = JVS( 365 )
   W( 71 ) = JVS( 366 )
   W( 75 ) = JVS( 367 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  JVS( 361) = W( 52 )
  JVS( 362) = W( 58 )
  JVS( 363) = W( 59 )
  JVS( 364) = W( 63 )
  JVS( 365) = W( 67 )
  JVS( 366) = W( 71 )
  JVS( 367) = W( 75 )
  IF ( ABS(  JVS( 378 )) < TINY(a) ) THEN
         IER = 60                                      
         RETURN
  END IF
   W( 13 ) = JVS( 368 )
   W( 21 ) = JVS( 369 )
   W( 24 ) = JVS( 370 )
   W( 26 ) = JVS( 371 )
   W( 33 ) = JVS( 372 )
   W( 48 ) = JVS( 373 )
   W( 50 ) = JVS( 374 )
   W( 56 ) = JVS( 375 )
   W( 57 ) = JVS( 376 )
   W( 58 ) = JVS( 377 )
   W( 60 ) = JVS( 378 )
   W( 62 ) = JVS( 379 )
   W( 63 ) = JVS( 380 )
   W( 64 ) = JVS( 381 )
   W( 65 ) = JVS( 382 )
   W( 66 ) = JVS( 383 )
   W( 67 ) = JVS( 384 )
   W( 68 ) = JVS( 385 )
   W( 70 ) = JVS( 386 )
   W( 71 ) = JVS( 387 )
   W( 72 ) = JVS( 388 )
   W( 73 ) = JVS( 389 )
   W( 75 ) = JVS( 390 )
   W( 79 ) = JVS( 391 )
  a = -W( 13 ) / JVS(           75  )
  W( 13 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 76 )
  a = -W( 21 ) / JVS(           99  )
  W( 21 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 100 )
  a = -W( 24 ) / JVS(          107  )
  W( 24 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 108 )
  a = -W( 26 ) / JVS(          113  )
  W( 26 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 114 )
  a = -W( 33 ) / JVS(          139  )
  W( 33 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 140 )
  a = -W( 48 ) / JVS(          259  )
  W( 48 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 260 )
  W( 67 ) = W( 67 ) + a*JVS( 261 )
  W( 71 ) = W( 71 ) + a*JVS( 262 )
  W( 75 ) = W( 75 ) + a*JVS( 263 )
  a = -W( 50 ) / JVS(          281  )
  W( 50 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 282 )
  W( 67 ) = W( 67 ) + a*JVS( 283 )
  W( 71 ) = W( 71 ) + a*JVS( 284 )
  W( 75 ) = W( 75 ) + a*JVS( 285 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  JVS( 368) = W( 13 )
  JVS( 369) = W( 21 )
  JVS( 370) = W( 24 )
  JVS( 371) = W( 26 )
  JVS( 372) = W( 33 )
  JVS( 373) = W( 48 )
  JVS( 374) = W( 50 )
  JVS( 375) = W( 56 )
  JVS( 376) = W( 57 )
  JVS( 377) = W( 58 )
  JVS( 378) = W( 60 )
  JVS( 379) = W( 62 )
  JVS( 380) = W( 63 )
  JVS( 381) = W( 64 )
  JVS( 382) = W( 65 )
  JVS( 383) = W( 66 )
  JVS( 384) = W( 67 )
  JVS( 385) = W( 68 )
  JVS( 386) = W( 70 )
  JVS( 387) = W( 71 )
  JVS( 388) = W( 72 )
  JVS( 389) = W( 73 )
  JVS( 390) = W( 75 )
  JVS( 391) = W( 79 )
  IF ( ABS(  JVS( 412 )) < TINY(a) ) THEN
         IER = 61                                      
         RETURN
  END IF
   W( 21 ) = JVS( 392 )
   W( 24 ) = JVS( 393 )
   W( 26 ) = JVS( 394 )
   W( 28 ) = JVS( 395 )
   W( 29 ) = JVS( 396 )
   W( 30 ) = JVS( 397 )
   W( 33 ) = JVS( 398 )
   W( 39 ) = JVS( 399 )
   W( 46 ) = JVS( 400 )
   W( 48 ) = JVS( 401 )
   W( 49 ) = JVS( 402 )
   W( 50 ) = JVS( 403 )
   W( 51 ) = JVS( 404 )
   W( 52 ) = JVS( 405 )
   W( 54 ) = JVS( 406 )
   W( 55 ) = JVS( 407 )
   W( 56 ) = JVS( 408 )
   W( 57 ) = JVS( 409 )
   W( 58 ) = JVS( 410 )
   W( 59 ) = JVS( 411 )
   W( 61 ) = JVS( 412 )
   W( 62 ) = JVS( 413 )
   W( 63 ) = JVS( 414 )
   W( 65 ) = JVS( 415 )
   W( 66 ) = JVS( 416 )
   W( 67 ) = JVS( 417 )
   W( 68 ) = JVS( 418 )
   W( 69 ) = JVS( 419 )
   W( 70 ) = JVS( 420 )
   W( 71 ) = JVS( 421 )
   W( 72 ) = JVS( 422 )
   W( 73 ) = JVS( 423 )
   W( 74 ) = JVS( 424 )
   W( 75 ) = JVS( 425 )
   W( 76 ) = JVS( 426 )
   W( 77 ) = JVS( 427 )
   W( 78 ) = JVS( 428 )
   W( 79 ) = JVS( 429 )
  a = -W( 21 ) / JVS(           99  )
  W( 21 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 100 )
  a = -W( 24 ) / JVS(          107  )
  W( 24 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 108 )
  a = -W( 26 ) / JVS(          113  )
  W( 26 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 114 )
  a = -W( 28 ) / JVS(          117  )
  W( 28 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 118 )
  W( 71 ) = W( 71 ) + a*JVS( 119 )
  W( 76 ) = W( 76 ) + a*JVS( 120 )
  W( 77 ) = W( 77 ) + a*JVS( 121 )
  a = -W( 29 ) / JVS(          122  )
  W( 29 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 123 )
  W( 76 ) = W( 76 ) + a*JVS( 124 )
  W( 78 ) = W( 78 ) + a*JVS( 125 )
  a = -W( 30 ) / JVS(          126  )
  W( 30 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 127 )
  W( 73 ) = W( 73 ) + a*JVS( 128 )
  W( 78 ) = W( 78 ) + a*JVS( 129 )
  a = -W( 33 ) / JVS(          139  )
  W( 33 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 140 )
  a = -W( 39 ) / JVS(          161  )
  W( 39 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 162 )
  W( 71 ) = W( 71 ) + a*JVS( 163 )
  a = -W( 46 ) / JVS(          229  )
  W( 46 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 230 )
  W( 56 ) = W( 56 ) + a*JVS( 231 )
  W( 58 ) = W( 58 ) + a*JVS( 232 )
  W( 62 ) = W( 62 ) + a*JVS( 233 )
  W( 67 ) = W( 67 ) + a*JVS( 234 )
  W( 71 ) = W( 71 ) + a*JVS( 235 )
  W( 74 ) = W( 74 ) + a*JVS( 236 )
  W( 75 ) = W( 75 ) + a*JVS( 237 )
  a = -W( 48 ) / JVS(          259  )
  W( 48 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 260 )
  W( 67 ) = W( 67 ) + a*JVS( 261 )
  W( 71 ) = W( 71 ) + a*JVS( 262 )
  W( 75 ) = W( 75 ) + a*JVS( 263 )
  a = -W( 49 ) / JVS(          272  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 273 )
  W( 54 ) = W( 54 ) + a*JVS( 274 )
  W( 57 ) = W( 57 ) + a*JVS( 275 )
  W( 63 ) = W( 63 ) + a*JVS( 276 )
  W( 67 ) = W( 67 ) + a*JVS( 277 )
  W( 71 ) = W( 71 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 78 ) = W( 78 ) + a*JVS( 280 )
  a = -W( 50 ) / JVS(          281  )
  W( 50 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 282 )
  W( 67 ) = W( 67 ) + a*JVS( 283 )
  W( 71 ) = W( 71 ) + a*JVS( 284 )
  W( 75 ) = W( 75 ) + a*JVS( 285 )
  a = -W( 51 ) / JVS(          288  )
  W( 51 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 289 )
  W( 70 ) = W( 70 ) + a*JVS( 290 )
  W( 71 ) = W( 71 ) + a*JVS( 291 )
  W( 72 ) = W( 72 ) + a*JVS( 292 )
  W( 73 ) = W( 73 ) + a*JVS( 293 )
  W( 74 ) = W( 74 ) + a*JVS( 294 )
  W( 75 ) = W( 75 ) + a*JVS( 295 )
  W( 78 ) = W( 78 ) + a*JVS( 296 )
  W( 79 ) = W( 79 ) + a*JVS( 297 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 55 ) / JVS(          338  )
  W( 55 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 339 )
  W( 63 ) = W( 63 ) + a*JVS( 340 )
  W( 67 ) = W( 67 ) + a*JVS( 341 )
  W( 71 ) = W( 71 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  JVS( 392) = W( 21 )
  JVS( 393) = W( 24 )
  JVS( 394) = W( 26 )
  JVS( 395) = W( 28 )
  JVS( 396) = W( 29 )
  JVS( 397) = W( 30 )
  JVS( 398) = W( 33 )
  JVS( 399) = W( 39 )
  JVS( 400) = W( 46 )
  JVS( 401) = W( 48 )
  JVS( 402) = W( 49 )
  JVS( 403) = W( 50 )
  JVS( 404) = W( 51 )
  JVS( 405) = W( 52 )
  JVS( 406) = W( 54 )
  JVS( 407) = W( 55 )
  JVS( 408) = W( 56 )
  JVS( 409) = W( 57 )
  JVS( 410) = W( 58 )
  JVS( 411) = W( 59 )
  JVS( 412) = W( 61 )
  JVS( 413) = W( 62 )
  JVS( 414) = W( 63 )
  JVS( 415) = W( 65 )
  JVS( 416) = W( 66 )
  JVS( 417) = W( 67 )
  JVS( 418) = W( 68 )
  JVS( 419) = W( 69 )
  JVS( 420) = W( 70 )
  JVS( 421) = W( 71 )
  JVS( 422) = W( 72 )
  JVS( 423) = W( 73 )
  JVS( 424) = W( 74 )
  JVS( 425) = W( 75 )
  JVS( 426) = W( 76 )
  JVS( 427) = W( 77 )
  JVS( 428) = W( 78 )
  JVS( 429) = W( 79 )
  IF ( ABS(  JVS( 435 )) < TINY(a) ) THEN
         IER = 62                                      
         RETURN
  END IF
   W( 25 ) = JVS( 430 )
   W( 54 ) = JVS( 431 )
   W( 56 ) = JVS( 432 )
   W( 57 ) = JVS( 433 )
   W( 58 ) = JVS( 434 )
   W( 62 ) = JVS( 435 )
   W( 63 ) = JVS( 436 )
   W( 67 ) = JVS( 437 )
   W( 69 ) = JVS( 438 )
   W( 71 ) = JVS( 439 )
   W( 73 ) = JVS( 440 )
   W( 74 ) = JVS( 441 )
   W( 75 ) = JVS( 442 )
  a = -W( 25 ) / JVS(          110  )
  W( 25 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 111 )
  W( 74 ) = W( 74 ) + a*JVS( 112 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  JVS( 430) = W( 25 )
  JVS( 431) = W( 54 )
  JVS( 432) = W( 56 )
  JVS( 433) = W( 57 )
  JVS( 434) = W( 58 )
  JVS( 435) = W( 62 )
  JVS( 436) = W( 63 )
  JVS( 437) = W( 67 )
  JVS( 438) = W( 69 )
  JVS( 439) = W( 71 )
  JVS( 440) = W( 73 )
  JVS( 441) = W( 74 )
  JVS( 442) = W( 75 )
  IF ( ABS(  JVS( 452 )) < TINY(a) ) THEN
         IER = 63                                      
         RETURN
  END IF
   W( 11 ) = JVS( 443 )
   W( 48 ) = JVS( 444 )
   W( 50 ) = JVS( 445 )
   W( 52 ) = JVS( 446 )
   W( 54 ) = JVS( 447 )
   W( 55 ) = JVS( 448 )
   W( 56 ) = JVS( 449 )
   W( 58 ) = JVS( 450 )
   W( 59 ) = JVS( 451 )
   W( 63 ) = JVS( 452 )
   W( 67 ) = JVS( 453 )
   W( 71 ) = JVS( 454 )
   W( 73 ) = JVS( 455 )
   W( 74 ) = JVS( 456 )
   W( 75 ) = JVS( 457 )
  a = -W( 11 ) / JVS(           71  )
  W( 11 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 72 )
  a = -W( 48 ) / JVS(          259  )
  W( 48 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 260 )
  W( 67 ) = W( 67 ) + a*JVS( 261 )
  W( 71 ) = W( 71 ) + a*JVS( 262 )
  W( 75 ) = W( 75 ) + a*JVS( 263 )
  a = -W( 50 ) / JVS(          281  )
  W( 50 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 282 )
  W( 67 ) = W( 67 ) + a*JVS( 283 )
  W( 71 ) = W( 71 ) + a*JVS( 284 )
  W( 75 ) = W( 75 ) + a*JVS( 285 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 55 ) / JVS(          338  )
  W( 55 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 339 )
  W( 63 ) = W( 63 ) + a*JVS( 340 )
  W( 67 ) = W( 67 ) + a*JVS( 341 )
  W( 71 ) = W( 71 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  JVS( 443) = W( 11 )
  JVS( 444) = W( 48 )
  JVS( 445) = W( 50 )
  JVS( 446) = W( 52 )
  JVS( 447) = W( 54 )
  JVS( 448) = W( 55 )
  JVS( 449) = W( 56 )
  JVS( 450) = W( 58 )
  JVS( 451) = W( 59 )
  JVS( 452) = W( 63 )
  JVS( 453) = W( 67 )
  JVS( 454) = W( 71 )
  JVS( 455) = W( 73 )
  JVS( 456) = W( 74 )
  JVS( 457) = W( 75 )
  IF ( ABS(  JVS( 476 )) < TINY(a) ) THEN
         IER = 64                                      
         RETURN
  END IF
   W( 20 ) = JVS( 458 )
   W( 24 ) = JVS( 459 )
   W( 26 ) = JVS( 460 )
   W( 33 ) = JVS( 461 )
   W( 35 ) = JVS( 462 )
   W( 36 ) = JVS( 463 )
   W( 38 ) = JVS( 464 )
   W( 42 ) = JVS( 465 )
   W( 48 ) = JVS( 466 )
   W( 50 ) = JVS( 467 )
   W( 54 ) = JVS( 468 )
   W( 55 ) = JVS( 469 )
   W( 56 ) = JVS( 470 )
   W( 57 ) = JVS( 471 )
   W( 58 ) = JVS( 472 )
   W( 59 ) = JVS( 473 )
   W( 62 ) = JVS( 474 )
   W( 63 ) = JVS( 475 )
   W( 64 ) = JVS( 476 )
   W( 65 ) = JVS( 477 )
   W( 66 ) = JVS( 478 )
   W( 67 ) = JVS( 479 )
   W( 69 ) = JVS( 480 )
   W( 71 ) = JVS( 481 )
   W( 73 ) = JVS( 482 )
   W( 74 ) = JVS( 483 )
   W( 75 ) = JVS( 484 )
   W( 77 ) = JVS( 485 )
   W( 78 ) = JVS( 486 )
  a = -W( 20 ) / JVS(           97  )
  W( 20 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 98 )
  a = -W( 24 ) / JVS(          107  )
  W( 24 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 108 )
  a = -W( 26 ) / JVS(          113  )
  W( 26 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 114 )
  a = -W( 33 ) / JVS(          139  )
  W( 33 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 140 )
  a = -W( 35 ) / JVS(          145  )
  W( 35 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 146 )
  a = -W( 36 ) / JVS(          149  )
  W( 36 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 150 )
  a = -W( 38 ) / JVS(          158  )
  W( 38 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 159 )
  W( 71 ) = W( 71 ) + a*JVS( 160 )
  a = -W( 42 ) / JVS(          177  )
  W( 42 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 178 )
  W( 71 ) = W( 71 ) + a*JVS( 179 )
  W( 77 ) = W( 77 ) + a*JVS( 180 )
  W( 78 ) = W( 78 ) + a*JVS( 181 )
  a = -W( 48 ) / JVS(          259  )
  W( 48 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 260 )
  W( 67 ) = W( 67 ) + a*JVS( 261 )
  W( 71 ) = W( 71 ) + a*JVS( 262 )
  W( 75 ) = W( 75 ) + a*JVS( 263 )
  a = -W( 50 ) / JVS(          281  )
  W( 50 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 282 )
  W( 67 ) = W( 67 ) + a*JVS( 283 )
  W( 71 ) = W( 71 ) + a*JVS( 284 )
  W( 75 ) = W( 75 ) + a*JVS( 285 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 55 ) / JVS(          338  )
  W( 55 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 339 )
  W( 63 ) = W( 63 ) + a*JVS( 340 )
  W( 67 ) = W( 67 ) + a*JVS( 341 )
  W( 71 ) = W( 71 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 62 ) / JVS(          435  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 436 )
  W( 67 ) = W( 67 ) + a*JVS( 437 )
  W( 69 ) = W( 69 ) + a*JVS( 438 )
  W( 71 ) = W( 71 ) + a*JVS( 439 )
  W( 73 ) = W( 73 ) + a*JVS( 440 )
  W( 74 ) = W( 74 ) + a*JVS( 441 )
  W( 75 ) = W( 75 ) + a*JVS( 442 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  JVS( 458) = W( 20 )
  JVS( 459) = W( 24 )
  JVS( 460) = W( 26 )
  JVS( 461) = W( 33 )
  JVS( 462) = W( 35 )
  JVS( 463) = W( 36 )
  JVS( 464) = W( 38 )
  JVS( 465) = W( 42 )
  JVS( 466) = W( 48 )
  JVS( 467) = W( 50 )
  JVS( 468) = W( 54 )
  JVS( 469) = W( 55 )
  JVS( 470) = W( 56 )
  JVS( 471) = W( 57 )
  JVS( 472) = W( 58 )
  JVS( 473) = W( 59 )
  JVS( 474) = W( 62 )
  JVS( 475) = W( 63 )
  JVS( 476) = W( 64 )
  JVS( 477) = W( 65 )
  JVS( 478) = W( 66 )
  JVS( 479) = W( 67 )
  JVS( 480) = W( 69 )
  JVS( 481) = W( 71 )
  JVS( 482) = W( 73 )
  JVS( 483) = W( 74 )
  JVS( 484) = W( 75 )
  JVS( 485) = W( 77 )
  JVS( 486) = W( 78 )
  IF ( ABS(  JVS( 498 )) < TINY(a) ) THEN
         IER = 65                                      
         RETURN
  END IF
   W( 24 ) = JVS( 487 )
   W( 26 ) = JVS( 488 )
   W( 33 ) = JVS( 489 )
   W( 50 ) = JVS( 490 )
   W( 55 ) = JVS( 491 )
   W( 56 ) = JVS( 492 )
   W( 57 ) = JVS( 493 )
   W( 58 ) = JVS( 494 )
   W( 59 ) = JVS( 495 )
   W( 62 ) = JVS( 496 )
   W( 63 ) = JVS( 497 )
   W( 65 ) = JVS( 498 )
   W( 66 ) = JVS( 499 )
   W( 67 ) = JVS( 500 )
   W( 69 ) = JVS( 501 )
   W( 71 ) = JVS( 502 )
   W( 73 ) = JVS( 503 )
   W( 74 ) = JVS( 504 )
   W( 75 ) = JVS( 505 )
   W( 76 ) = JVS( 506 )
   W( 77 ) = JVS( 507 )
  a = -W( 24 ) / JVS(          107  )
  W( 24 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 108 )
  a = -W( 26 ) / JVS(          113  )
  W( 26 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 114 )
  a = -W( 33 ) / JVS(          139  )
  W( 33 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 140 )
  a = -W( 50 ) / JVS(          281  )
  W( 50 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 282 )
  W( 67 ) = W( 67 ) + a*JVS( 283 )
  W( 71 ) = W( 71 ) + a*JVS( 284 )
  W( 75 ) = W( 75 ) + a*JVS( 285 )
  a = -W( 55 ) / JVS(          338  )
  W( 55 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 339 )
  W( 63 ) = W( 63 ) + a*JVS( 340 )
  W( 67 ) = W( 67 ) + a*JVS( 341 )
  W( 71 ) = W( 71 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 62 ) / JVS(          435  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 436 )
  W( 67 ) = W( 67 ) + a*JVS( 437 )
  W( 69 ) = W( 69 ) + a*JVS( 438 )
  W( 71 ) = W( 71 ) + a*JVS( 439 )
  W( 73 ) = W( 73 ) + a*JVS( 440 )
  W( 74 ) = W( 74 ) + a*JVS( 441 )
  W( 75 ) = W( 75 ) + a*JVS( 442 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  JVS( 487) = W( 24 )
  JVS( 488) = W( 26 )
  JVS( 489) = W( 33 )
  JVS( 490) = W( 50 )
  JVS( 491) = W( 55 )
  JVS( 492) = W( 56 )
  JVS( 493) = W( 57 )
  JVS( 494) = W( 58 )
  JVS( 495) = W( 59 )
  JVS( 496) = W( 62 )
  JVS( 497) = W( 63 )
  JVS( 498) = W( 65 )
  JVS( 499) = W( 66 )
  JVS( 500) = W( 67 )
  JVS( 501) = W( 69 )
  JVS( 502) = W( 71 )
  JVS( 503) = W( 73 )
  JVS( 504) = W( 74 )
  JVS( 505) = W( 75 )
  JVS( 506) = W( 76 )
  JVS( 507) = W( 77 )
  IF ( ABS(  JVS( 519 )) < TINY(a) ) THEN
         IER = 66                                      
         RETURN
  END IF
   W( 26 ) = JVS( 508 )
   W( 33 ) = JVS( 509 )
   W( 34 ) = JVS( 510 )
   W( 52 ) = JVS( 511 )
   W( 54 ) = JVS( 512 )
   W( 56 ) = JVS( 513 )
   W( 57 ) = JVS( 514 )
   W( 58 ) = JVS( 515 )
   W( 59 ) = JVS( 516 )
   W( 62 ) = JVS( 517 )
   W( 63 ) = JVS( 518 )
   W( 66 ) = JVS( 519 )
   W( 67 ) = JVS( 520 )
   W( 68 ) = JVS( 521 )
   W( 69 ) = JVS( 522 )
   W( 71 ) = JVS( 523 )
   W( 72 ) = JVS( 524 )
   W( 73 ) = JVS( 525 )
   W( 74 ) = JVS( 526 )
   W( 75 ) = JVS( 527 )
   W( 76 ) = JVS( 528 )
   W( 77 ) = JVS( 529 )
   W( 79 ) = JVS( 530 )
  a = -W( 26 ) / JVS(          113  )
  W( 26 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 114 )
  a = -W( 33 ) / JVS(          139  )
  W( 33 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 140 )
  a = -W( 34 ) / JVS(          141  )
  W( 34 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 142 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 62 ) / JVS(          435  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 436 )
  W( 67 ) = W( 67 ) + a*JVS( 437 )
  W( 69 ) = W( 69 ) + a*JVS( 438 )
  W( 71 ) = W( 71 ) + a*JVS( 439 )
  W( 73 ) = W( 73 ) + a*JVS( 440 )
  W( 74 ) = W( 74 ) + a*JVS( 441 )
  W( 75 ) = W( 75 ) + a*JVS( 442 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  JVS( 508) = W( 26 )
  JVS( 509) = W( 33 )
  JVS( 510) = W( 34 )
  JVS( 511) = W( 52 )
  JVS( 512) = W( 54 )
  JVS( 513) = W( 56 )
  JVS( 514) = W( 57 )
  JVS( 515) = W( 58 )
  JVS( 516) = W( 59 )
  JVS( 517) = W( 62 )
  JVS( 518) = W( 63 )
  JVS( 519) = W( 66 )
  JVS( 520) = W( 67 )
  JVS( 521) = W( 68 )
  JVS( 522) = W( 69 )
  JVS( 523) = W( 71 )
  JVS( 524) = W( 72 )
  JVS( 525) = W( 73 )
  JVS( 526) = W( 74 )
  JVS( 527) = W( 75 )
  JVS( 528) = W( 76 )
  JVS( 529) = W( 77 )
  JVS( 530) = W( 79 )
  IF ( ABS(  JVS( 543 )) < TINY(a) ) THEN
         IER = 67                                      
         RETURN
  END IF
   W( 38 ) = JVS( 531 )
   W( 39 ) = JVS( 532 )
   W( 48 ) = JVS( 533 )
   W( 50 ) = JVS( 534 )
   W( 52 ) = JVS( 535 )
   W( 54 ) = JVS( 536 )
   W( 55 ) = JVS( 537 )
   W( 56 ) = JVS( 538 )
   W( 57 ) = JVS( 539 )
   W( 58 ) = JVS( 540 )
   W( 59 ) = JVS( 541 )
   W( 63 ) = JVS( 542 )
   W( 67 ) = JVS( 543 )
   W( 68 ) = JVS( 544 )
   W( 70 ) = JVS( 545 )
   W( 71 ) = JVS( 546 )
   W( 72 ) = JVS( 547 )
   W( 73 ) = JVS( 548 )
   W( 74 ) = JVS( 549 )
   W( 75 ) = JVS( 550 )
   W( 78 ) = JVS( 551 )
   W( 79 ) = JVS( 552 )
  a = -W( 38 ) / JVS(          158  )
  W( 38 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 159 )
  W( 71 ) = W( 71 ) + a*JVS( 160 )
  a = -W( 39 ) / JVS(          161  )
  W( 39 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 162 )
  W( 71 ) = W( 71 ) + a*JVS( 163 )
  a = -W( 48 ) / JVS(          259  )
  W( 48 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 260 )
  W( 67 ) = W( 67 ) + a*JVS( 261 )
  W( 71 ) = W( 71 ) + a*JVS( 262 )
  W( 75 ) = W( 75 ) + a*JVS( 263 )
  a = -W( 50 ) / JVS(          281  )
  W( 50 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 282 )
  W( 67 ) = W( 67 ) + a*JVS( 283 )
  W( 71 ) = W( 71 ) + a*JVS( 284 )
  W( 75 ) = W( 75 ) + a*JVS( 285 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 55 ) / JVS(          338  )
  W( 55 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 339 )
  W( 63 ) = W( 63 ) + a*JVS( 340 )
  W( 67 ) = W( 67 ) + a*JVS( 341 )
  W( 71 ) = W( 71 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  JVS( 531) = W( 38 )
  JVS( 532) = W( 39 )
  JVS( 533) = W( 48 )
  JVS( 534) = W( 50 )
  JVS( 535) = W( 52 )
  JVS( 536) = W( 54 )
  JVS( 537) = W( 55 )
  JVS( 538) = W( 56 )
  JVS( 539) = W( 57 )
  JVS( 540) = W( 58 )
  JVS( 541) = W( 59 )
  JVS( 542) = W( 63 )
  JVS( 543) = W( 67 )
  JVS( 544) = W( 68 )
  JVS( 545) = W( 70 )
  JVS( 546) = W( 71 )
  JVS( 547) = W( 72 )
  JVS( 548) = W( 73 )
  JVS( 549) = W( 74 )
  JVS( 550) = W( 75 )
  JVS( 551) = W( 78 )
  JVS( 552) = W( 79 )
  IF ( ABS(  JVS( 574 )) < TINY(a) ) THEN
         IER = 68                                      
         RETURN
  END IF
   W( 14 ) = JVS( 553 )
   W( 15 ) = JVS( 554 )
   W( 27 ) = JVS( 555 )
   W( 33 ) = JVS( 556 )
   W( 35 ) = JVS( 557 )
   W( 36 ) = JVS( 558 )
   W( 44 ) = JVS( 559 )
   W( 46 ) = JVS( 560 )
   W( 54 ) = JVS( 561 )
   W( 55 ) = JVS( 562 )
   W( 56 ) = JVS( 563 )
   W( 57 ) = JVS( 564 )
   W( 58 ) = JVS( 565 )
   W( 59 ) = JVS( 566 )
   W( 60 ) = JVS( 567 )
   W( 62 ) = JVS( 568 )
   W( 63 ) = JVS( 569 )
   W( 64 ) = JVS( 570 )
   W( 65 ) = JVS( 571 )
   W( 66 ) = JVS( 572 )
   W( 67 ) = JVS( 573 )
   W( 68 ) = JVS( 574 )
   W( 69 ) = JVS( 575 )
   W( 70 ) = JVS( 576 )
   W( 71 ) = JVS( 577 )
   W( 72 ) = JVS( 578 )
   W( 73 ) = JVS( 579 )
   W( 74 ) = JVS( 580 )
   W( 75 ) = JVS( 581 )
   W( 76 ) = JVS( 582 )
   W( 77 ) = JVS( 583 )
   W( 78 ) = JVS( 584 )
   W( 79 ) = JVS( 585 )
  a = -W( 14 ) / JVS(           77  )
  W( 14 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 78 )
  W( 54 ) = W( 54 ) + a*JVS( 79 )
  W( 67 ) = W( 67 ) + a*JVS( 80 )
  W( 71 ) = W( 71 ) + a*JVS( 81 )
  a = -W( 15 ) / JVS(           82  )
  W( 15 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 83 )
  W( 74 ) = W( 74 ) + a*JVS( 84 )
  a = -W( 27 ) / JVS(          115  )
  W( 27 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 116 )
  a = -W( 33 ) / JVS(          139  )
  W( 33 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 140 )
  a = -W( 35 ) / JVS(          145  )
  W( 35 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 146 )
  a = -W( 36 ) / JVS(          149  )
  W( 36 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 150 )
  a = -W( 44 ) / JVS(          193  )
  W( 44 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  W( 57 ) = W( 57 ) + a*JVS( 195 )
  W( 59 ) = W( 59 ) + a*JVS( 196 )
  W( 67 ) = W( 67 ) + a*JVS( 197 )
  W( 71 ) = W( 71 ) + a*JVS( 198 )
  W( 75 ) = W( 75 ) + a*JVS( 199 )
  a = -W( 46 ) / JVS(          229  )
  W( 46 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 230 )
  W( 56 ) = W( 56 ) + a*JVS( 231 )
  W( 58 ) = W( 58 ) + a*JVS( 232 )
  W( 62 ) = W( 62 ) + a*JVS( 233 )
  W( 67 ) = W( 67 ) + a*JVS( 234 )
  W( 71 ) = W( 71 ) + a*JVS( 235 )
  W( 74 ) = W( 74 ) + a*JVS( 236 )
  W( 75 ) = W( 75 ) + a*JVS( 237 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 55 ) / JVS(          338  )
  W( 55 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 339 )
  W( 63 ) = W( 63 ) + a*JVS( 340 )
  W( 67 ) = W( 67 ) + a*JVS( 341 )
  W( 71 ) = W( 71 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 60 ) / JVS(          378  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 379 )
  W( 63 ) = W( 63 ) + a*JVS( 380 )
  W( 64 ) = W( 64 ) + a*JVS( 381 )
  W( 65 ) = W( 65 ) + a*JVS( 382 )
  W( 66 ) = W( 66 ) + a*JVS( 383 )
  W( 67 ) = W( 67 ) + a*JVS( 384 )
  W( 68 ) = W( 68 ) + a*JVS( 385 )
  W( 70 ) = W( 70 ) + a*JVS( 386 )
  W( 71 ) = W( 71 ) + a*JVS( 387 )
  W( 72 ) = W( 72 ) + a*JVS( 388 )
  W( 73 ) = W( 73 ) + a*JVS( 389 )
  W( 75 ) = W( 75 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 62 ) / JVS(          435  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 436 )
  W( 67 ) = W( 67 ) + a*JVS( 437 )
  W( 69 ) = W( 69 ) + a*JVS( 438 )
  W( 71 ) = W( 71 ) + a*JVS( 439 )
  W( 73 ) = W( 73 ) + a*JVS( 440 )
  W( 74 ) = W( 74 ) + a*JVS( 441 )
  W( 75 ) = W( 75 ) + a*JVS( 442 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  a = -W( 64 ) / JVS(          476  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 477 )
  W( 66 ) = W( 66 ) + a*JVS( 478 )
  W( 67 ) = W( 67 ) + a*JVS( 479 )
  W( 69 ) = W( 69 ) + a*JVS( 480 )
  W( 71 ) = W( 71 ) + a*JVS( 481 )
  W( 73 ) = W( 73 ) + a*JVS( 482 )
  W( 74 ) = W( 74 ) + a*JVS( 483 )
  W( 75 ) = W( 75 ) + a*JVS( 484 )
  W( 77 ) = W( 77 ) + a*JVS( 485 )
  W( 78 ) = W( 78 ) + a*JVS( 486 )
  a = -W( 65 ) / JVS(          498  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 499 )
  W( 67 ) = W( 67 ) + a*JVS( 500 )
  W( 69 ) = W( 69 ) + a*JVS( 501 )
  W( 71 ) = W( 71 ) + a*JVS( 502 )
  W( 73 ) = W( 73 ) + a*JVS( 503 )
  W( 74 ) = W( 74 ) + a*JVS( 504 )
  W( 75 ) = W( 75 ) + a*JVS( 505 )
  W( 76 ) = W( 76 ) + a*JVS( 506 )
  W( 77 ) = W( 77 ) + a*JVS( 507 )
  a = -W( 66 ) / JVS(          519  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 520 )
  W( 68 ) = W( 68 ) + a*JVS( 521 )
  W( 69 ) = W( 69 ) + a*JVS( 522 )
  W( 71 ) = W( 71 ) + a*JVS( 523 )
  W( 72 ) = W( 72 ) + a*JVS( 524 )
  W( 73 ) = W( 73 ) + a*JVS( 525 )
  W( 74 ) = W( 74 ) + a*JVS( 526 )
  W( 75 ) = W( 75 ) + a*JVS( 527 )
  W( 76 ) = W( 76 ) + a*JVS( 528 )
  W( 77 ) = W( 77 ) + a*JVS( 529 )
  W( 79 ) = W( 79 ) + a*JVS( 530 )
  a = -W( 67 ) / JVS(          543  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 544 )
  W( 70 ) = W( 70 ) + a*JVS( 545 )
  W( 71 ) = W( 71 ) + a*JVS( 546 )
  W( 72 ) = W( 72 ) + a*JVS( 547 )
  W( 73 ) = W( 73 ) + a*JVS( 548 )
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 78 ) = W( 78 ) + a*JVS( 551 )
  W( 79 ) = W( 79 ) + a*JVS( 552 )
  JVS( 553) = W( 14 )
  JVS( 554) = W( 15 )
  JVS( 555) = W( 27 )
  JVS( 556) = W( 33 )
  JVS( 557) = W( 35 )
  JVS( 558) = W( 36 )
  JVS( 559) = W( 44 )
  JVS( 560) = W( 46 )
  JVS( 561) = W( 54 )
  JVS( 562) = W( 55 )
  JVS( 563) = W( 56 )
  JVS( 564) = W( 57 )
  JVS( 565) = W( 58 )
  JVS( 566) = W( 59 )
  JVS( 567) = W( 60 )
  JVS( 568) = W( 62 )
  JVS( 569) = W( 63 )
  JVS( 570) = W( 64 )
  JVS( 571) = W( 65 )
  JVS( 572) = W( 66 )
  JVS( 573) = W( 67 )
  JVS( 574) = W( 68 )
  JVS( 575) = W( 69 )
  JVS( 576) = W( 70 )
  JVS( 577) = W( 71 )
  JVS( 578) = W( 72 )
  JVS( 579) = W( 73 )
  JVS( 580) = W( 74 )
  JVS( 581) = W( 75 )
  JVS( 582) = W( 76 )
  JVS( 583) = W( 77 )
  JVS( 584) = W( 78 )
  JVS( 585) = W( 79 )
  IF ( ABS(  JVS( 606 )) < TINY(a) ) THEN
         IER = 69                                      
         RETURN
  END IF
   W( 20 ) = JVS( 586 )
   W( 24 ) = JVS( 587 )
   W( 26 ) = JVS( 588 )
   W( 27 ) = JVS( 589 )
   W( 33 ) = JVS( 590 )
   W( 34 ) = JVS( 591 )
   W( 50 ) = JVS( 592 )
   W( 52 ) = JVS( 593 )
   W( 54 ) = JVS( 594 )
   W( 56 ) = JVS( 595 )
   W( 57 ) = JVS( 596 )
   W( 58 ) = JVS( 597 )
   W( 59 ) = JVS( 598 )
   W( 62 ) = JVS( 599 )
   W( 63 ) = JVS( 600 )
   W( 64 ) = JVS( 601 )
   W( 65 ) = JVS( 602 )
   W( 66 ) = JVS( 603 )
   W( 67 ) = JVS( 604 )
   W( 68 ) = JVS( 605 )
   W( 69 ) = JVS( 606 )
   W( 70 ) = JVS( 607 )
   W( 71 ) = JVS( 608 )
   W( 72 ) = JVS( 609 )
   W( 73 ) = JVS( 610 )
   W( 74 ) = JVS( 611 )
   W( 75 ) = JVS( 612 )
   W( 76 ) = JVS( 613 )
   W( 77 ) = JVS( 614 )
   W( 78 ) = JVS( 615 )
   W( 79 ) = JVS( 616 )
  a = -W( 20 ) / JVS(           97  )
  W( 20 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 98 )
  a = -W( 24 ) / JVS(          107  )
  W( 24 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 108 )
  a = -W( 26 ) / JVS(          113  )
  W( 26 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 114 )
  a = -W( 27 ) / JVS(          115  )
  W( 27 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 116 )
  a = -W( 33 ) / JVS(          139  )
  W( 33 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 140 )
  a = -W( 34 ) / JVS(          141  )
  W( 34 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 142 )
  a = -W( 50 ) / JVS(          281  )
  W( 50 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 282 )
  W( 67 ) = W( 67 ) + a*JVS( 283 )
  W( 71 ) = W( 71 ) + a*JVS( 284 )
  W( 75 ) = W( 75 ) + a*JVS( 285 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 62 ) / JVS(          435  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 436 )
  W( 67 ) = W( 67 ) + a*JVS( 437 )
  W( 69 ) = W( 69 ) + a*JVS( 438 )
  W( 71 ) = W( 71 ) + a*JVS( 439 )
  W( 73 ) = W( 73 ) + a*JVS( 440 )
  W( 74 ) = W( 74 ) + a*JVS( 441 )
  W( 75 ) = W( 75 ) + a*JVS( 442 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  a = -W( 64 ) / JVS(          476  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 477 )
  W( 66 ) = W( 66 ) + a*JVS( 478 )
  W( 67 ) = W( 67 ) + a*JVS( 479 )
  W( 69 ) = W( 69 ) + a*JVS( 480 )
  W( 71 ) = W( 71 ) + a*JVS( 481 )
  W( 73 ) = W( 73 ) + a*JVS( 482 )
  W( 74 ) = W( 74 ) + a*JVS( 483 )
  W( 75 ) = W( 75 ) + a*JVS( 484 )
  W( 77 ) = W( 77 ) + a*JVS( 485 )
  W( 78 ) = W( 78 ) + a*JVS( 486 )
  a = -W( 65 ) / JVS(          498  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 499 )
  W( 67 ) = W( 67 ) + a*JVS( 500 )
  W( 69 ) = W( 69 ) + a*JVS( 501 )
  W( 71 ) = W( 71 ) + a*JVS( 502 )
  W( 73 ) = W( 73 ) + a*JVS( 503 )
  W( 74 ) = W( 74 ) + a*JVS( 504 )
  W( 75 ) = W( 75 ) + a*JVS( 505 )
  W( 76 ) = W( 76 ) + a*JVS( 506 )
  W( 77 ) = W( 77 ) + a*JVS( 507 )
  a = -W( 66 ) / JVS(          519  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 520 )
  W( 68 ) = W( 68 ) + a*JVS( 521 )
  W( 69 ) = W( 69 ) + a*JVS( 522 )
  W( 71 ) = W( 71 ) + a*JVS( 523 )
  W( 72 ) = W( 72 ) + a*JVS( 524 )
  W( 73 ) = W( 73 ) + a*JVS( 525 )
  W( 74 ) = W( 74 ) + a*JVS( 526 )
  W( 75 ) = W( 75 ) + a*JVS( 527 )
  W( 76 ) = W( 76 ) + a*JVS( 528 )
  W( 77 ) = W( 77 ) + a*JVS( 529 )
  W( 79 ) = W( 79 ) + a*JVS( 530 )
  a = -W( 67 ) / JVS(          543  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 544 )
  W( 70 ) = W( 70 ) + a*JVS( 545 )
  W( 71 ) = W( 71 ) + a*JVS( 546 )
  W( 72 ) = W( 72 ) + a*JVS( 547 )
  W( 73 ) = W( 73 ) + a*JVS( 548 )
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 78 ) = W( 78 ) + a*JVS( 551 )
  W( 79 ) = W( 79 ) + a*JVS( 552 )
  a = -W( 68 ) / JVS(          574  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 575 )
  W( 70 ) = W( 70 ) + a*JVS( 576 )
  W( 71 ) = W( 71 ) + a*JVS( 577 )
  W( 72 ) = W( 72 ) + a*JVS( 578 )
  W( 73 ) = W( 73 ) + a*JVS( 579 )
  W( 74 ) = W( 74 ) + a*JVS( 580 )
  W( 75 ) = W( 75 ) + a*JVS( 581 )
  W( 76 ) = W( 76 ) + a*JVS( 582 )
  W( 77 ) = W( 77 ) + a*JVS( 583 )
  W( 78 ) = W( 78 ) + a*JVS( 584 )
  W( 79 ) = W( 79 ) + a*JVS( 585 )
  JVS( 586) = W( 20 )
  JVS( 587) = W( 24 )
  JVS( 588) = W( 26 )
  JVS( 589) = W( 27 )
  JVS( 590) = W( 33 )
  JVS( 591) = W( 34 )
  JVS( 592) = W( 50 )
  JVS( 593) = W( 52 )
  JVS( 594) = W( 54 )
  JVS( 595) = W( 56 )
  JVS( 596) = W( 57 )
  JVS( 597) = W( 58 )
  JVS( 598) = W( 59 )
  JVS( 599) = W( 62 )
  JVS( 600) = W( 63 )
  JVS( 601) = W( 64 )
  JVS( 602) = W( 65 )
  JVS( 603) = W( 66 )
  JVS( 604) = W( 67 )
  JVS( 605) = W( 68 )
  JVS( 606) = W( 69 )
  JVS( 607) = W( 70 )
  JVS( 608) = W( 71 )
  JVS( 609) = W( 72 )
  JVS( 610) = W( 73 )
  JVS( 611) = W( 74 )
  JVS( 612) = W( 75 )
  JVS( 613) = W( 76 )
  JVS( 614) = W( 77 )
  JVS( 615) = W( 78 )
  JVS( 616) = W( 79 )
  IF ( ABS(  JVS( 627 )) < TINY(a) ) THEN
         IER = 70                                      
         RETURN
  END IF
   W( 18 ) = JVS( 617 )
   W( 52 ) = JVS( 618 )
   W( 55 ) = JVS( 619 )
   W( 57 ) = JVS( 620 )
   W( 58 ) = JVS( 621 )
   W( 59 ) = JVS( 622 )
   W( 63 ) = JVS( 623 )
   W( 67 ) = JVS( 624 )
   W( 68 ) = JVS( 625 )
   W( 69 ) = JVS( 626 )
   W( 70 ) = JVS( 627 )
   W( 71 ) = JVS( 628 )
   W( 72 ) = JVS( 629 )
   W( 73 ) = JVS( 630 )
   W( 74 ) = JVS( 631 )
   W( 75 ) = JVS( 632 )
   W( 76 ) = JVS( 633 )
   W( 77 ) = JVS( 634 )
   W( 78 ) = JVS( 635 )
   W( 79 ) = JVS( 636 )
  a = -W( 18 ) / JVS(           91  )
  W( 18 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 92 )
  W( 74 ) = W( 74 ) + a*JVS( 93 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  a = -W( 55 ) / JVS(          338  )
  W( 55 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 339 )
  W( 63 ) = W( 63 ) + a*JVS( 340 )
  W( 67 ) = W( 67 ) + a*JVS( 341 )
  W( 71 ) = W( 71 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  a = -W( 67 ) / JVS(          543  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 544 )
  W( 70 ) = W( 70 ) + a*JVS( 545 )
  W( 71 ) = W( 71 ) + a*JVS( 546 )
  W( 72 ) = W( 72 ) + a*JVS( 547 )
  W( 73 ) = W( 73 ) + a*JVS( 548 )
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 78 ) = W( 78 ) + a*JVS( 551 )
  W( 79 ) = W( 79 ) + a*JVS( 552 )
  a = -W( 68 ) / JVS(          574  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 575 )
  W( 70 ) = W( 70 ) + a*JVS( 576 )
  W( 71 ) = W( 71 ) + a*JVS( 577 )
  W( 72 ) = W( 72 ) + a*JVS( 578 )
  W( 73 ) = W( 73 ) + a*JVS( 579 )
  W( 74 ) = W( 74 ) + a*JVS( 580 )
  W( 75 ) = W( 75 ) + a*JVS( 581 )
  W( 76 ) = W( 76 ) + a*JVS( 582 )
  W( 77 ) = W( 77 ) + a*JVS( 583 )
  W( 78 ) = W( 78 ) + a*JVS( 584 )
  W( 79 ) = W( 79 ) + a*JVS( 585 )
  a = -W( 69 ) / JVS(          606  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 607 )
  W( 71 ) = W( 71 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  W( 74 ) = W( 74 ) + a*JVS( 611 )
  W( 75 ) = W( 75 ) + a*JVS( 612 )
  W( 76 ) = W( 76 ) + a*JVS( 613 )
  W( 77 ) = W( 77 ) + a*JVS( 614 )
  W( 78 ) = W( 78 ) + a*JVS( 615 )
  W( 79 ) = W( 79 ) + a*JVS( 616 )
  JVS( 617) = W( 18 )
  JVS( 618) = W( 52 )
  JVS( 619) = W( 55 )
  JVS( 620) = W( 57 )
  JVS( 621) = W( 58 )
  JVS( 622) = W( 59 )
  JVS( 623) = W( 63 )
  JVS( 624) = W( 67 )
  JVS( 625) = W( 68 )
  JVS( 626) = W( 69 )
  JVS( 627) = W( 70 )
  JVS( 628) = W( 71 )
  JVS( 629) = W( 72 )
  JVS( 630) = W( 73 )
  JVS( 631) = W( 74 )
  JVS( 632) = W( 75 )
  JVS( 633) = W( 76 )
  JVS( 634) = W( 77 )
  JVS( 635) = W( 78 )
  JVS( 636) = W( 79 )
  IF ( ABS(  JVS( 687 )) < TINY(a) ) THEN
         IER = 71                                      
         RETURN
  END IF
   W( 10 ) = JVS( 637 )
   W( 11 ) = JVS( 638 )
   W( 12 ) = JVS( 639 )
   W( 13 ) = JVS( 640 )
   W( 19 ) = JVS( 641 )
   W( 20 ) = JVS( 642 )
   W( 21 ) = JVS( 643 )
   W( 23 ) = JVS( 644 )
   W( 24 ) = JVS( 645 )
   W( 26 ) = JVS( 646 )
   W( 27 ) = JVS( 647 )
   W( 28 ) = JVS( 648 )
   W( 29 ) = JVS( 649 )
   W( 32 ) = JVS( 650 )
   W( 33 ) = JVS( 651 )
   W( 34 ) = JVS( 652 )
   W( 35 ) = JVS( 653 )
   W( 36 ) = JVS( 654 )
   W( 37 ) = JVS( 655 )
   W( 38 ) = JVS( 656 )
   W( 39 ) = JVS( 657 )
   W( 41 ) = JVS( 658 )
   W( 42 ) = JVS( 659 )
   W( 43 ) = JVS( 660 )
   W( 44 ) = JVS( 661 )
   W( 45 ) = JVS( 662 )
   W( 46 ) = JVS( 663 )
   W( 47 ) = JVS( 664 )
   W( 48 ) = JVS( 665 )
   W( 49 ) = JVS( 666 )
   W( 50 ) = JVS( 667 )
   W( 51 ) = JVS( 668 )
   W( 52 ) = JVS( 669 )
   W( 54 ) = JVS( 670 )
   W( 55 ) = JVS( 671 )
   W( 56 ) = JVS( 672 )
   W( 57 ) = JVS( 673 )
   W( 58 ) = JVS( 674 )
   W( 59 ) = JVS( 675 )
   W( 60 ) = JVS( 676 )
   W( 61 ) = JVS( 677 )
   W( 62 ) = JVS( 678 )
   W( 63 ) = JVS( 679 )
   W( 64 ) = JVS( 680 )
   W( 65 ) = JVS( 681 )
   W( 66 ) = JVS( 682 )
   W( 67 ) = JVS( 683 )
   W( 68 ) = JVS( 684 )
   W( 69 ) = JVS( 685 )
   W( 70 ) = JVS( 686 )
   W( 71 ) = JVS( 687 )
   W( 72 ) = JVS( 688 )
   W( 73 ) = JVS( 689 )
   W( 74 ) = JVS( 690 )
   W( 75 ) = JVS( 691 )
   W( 76 ) = JVS( 692 )
   W( 77 ) = JVS( 693 )
   W( 78 ) = JVS( 694 )
   W( 79 ) = JVS( 695 )
  a = -W( 10 ) / JVS(           69  )
  W( 10 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 70 )
  a = -W( 11 ) / JVS(           71  )
  W( 11 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 72 )
  a = -W( 12 ) / JVS(           73  )
  W( 12 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 74 )
  a = -W( 13 ) / JVS(           75  )
  W( 13 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 76 )
  a = -W( 19 ) / JVS(           94  )
  W( 19 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 95 )
  W( 78 ) = W( 78 ) + a*JVS( 96 )
  a = -W( 20 ) / JVS(           97  )
  W( 20 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 98 )
  a = -W( 21 ) / JVS(           99  )
  W( 21 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 100 )
  a = -W( 23 ) / JVS(          104  )
  W( 23 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 105 )
  W( 73 ) = W( 73 ) + a*JVS( 106 )
  a = -W( 24 ) / JVS(          107  )
  W( 24 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 108 )
  a = -W( 26 ) / JVS(          113  )
  W( 26 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 114 )
  a = -W( 27 ) / JVS(          115  )
  W( 27 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 116 )
  a = -W( 28 ) / JVS(          117  )
  W( 28 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 118 )
  W( 71 ) = W( 71 ) + a*JVS( 119 )
  W( 76 ) = W( 76 ) + a*JVS( 120 )
  W( 77 ) = W( 77 ) + a*JVS( 121 )
  a = -W( 29 ) / JVS(          122  )
  W( 29 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 123 )
  W( 76 ) = W( 76 ) + a*JVS( 124 )
  W( 78 ) = W( 78 ) + a*JVS( 125 )
  a = -W( 32 ) / JVS(          135  )
  W( 32 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 136 )
  W( 74 ) = W( 74 ) + a*JVS( 137 )
  W( 78 ) = W( 78 ) + a*JVS( 138 )
  a = -W( 33 ) / JVS(          139  )
  W( 33 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 140 )
  a = -W( 34 ) / JVS(          141  )
  W( 34 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 142 )
  a = -W( 35 ) / JVS(          145  )
  W( 35 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 146 )
  a = -W( 36 ) / JVS(          149  )
  W( 36 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 150 )
  a = -W( 37 ) / JVS(          153  )
  W( 37 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 154 )
  W( 75 ) = W( 75 ) + a*JVS( 155 )
  a = -W( 38 ) / JVS(          158  )
  W( 38 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 159 )
  W( 71 ) = W( 71 ) + a*JVS( 160 )
  a = -W( 39 ) / JVS(          161  )
  W( 39 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 162 )
  W( 71 ) = W( 71 ) + a*JVS( 163 )
  a = -W( 41 ) / JVS(          172  )
  W( 41 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 173 )
  W( 67 ) = W( 67 ) + a*JVS( 174 )
  W( 71 ) = W( 71 ) + a*JVS( 175 )
  W( 75 ) = W( 75 ) + a*JVS( 176 )
  a = -W( 42 ) / JVS(          177  )
  W( 42 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 178 )
  W( 71 ) = W( 71 ) + a*JVS( 179 )
  W( 77 ) = W( 77 ) + a*JVS( 180 )
  W( 78 ) = W( 78 ) + a*JVS( 181 )
  a = -W( 43 ) / JVS(          183  )
  W( 43 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 184 )
  W( 71 ) = W( 71 ) + a*JVS( 185 )
  W( 75 ) = W( 75 ) + a*JVS( 186 )
  W( 78 ) = W( 78 ) + a*JVS( 187 )
  a = -W( 44 ) / JVS(          193  )
  W( 44 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  W( 57 ) = W( 57 ) + a*JVS( 195 )
  W( 59 ) = W( 59 ) + a*JVS( 196 )
  W( 67 ) = W( 67 ) + a*JVS( 197 )
  W( 71 ) = W( 71 ) + a*JVS( 198 )
  W( 75 ) = W( 75 ) + a*JVS( 199 )
  a = -W( 45 ) / JVS(          206  )
  W( 45 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 207 )
  W( 49 ) = W( 49 ) + a*JVS( 208 )
  W( 50 ) = W( 50 ) + a*JVS( 209 )
  W( 52 ) = W( 52 ) + a*JVS( 210 )
  W( 54 ) = W( 54 ) + a*JVS( 211 )
  W( 55 ) = W( 55 ) + a*JVS( 212 )
  W( 56 ) = W( 56 ) + a*JVS( 213 )
  W( 57 ) = W( 57 ) + a*JVS( 214 )
  W( 58 ) = W( 58 ) + a*JVS( 215 )
  W( 59 ) = W( 59 ) + a*JVS( 216 )
  W( 60 ) = W( 60 ) + a*JVS( 217 )
  W( 61 ) = W( 61 ) + a*JVS( 218 )
  W( 63 ) = W( 63 ) + a*JVS( 219 )
  W( 64 ) = W( 64 ) + a*JVS( 220 )
  W( 67 ) = W( 67 ) + a*JVS( 221 )
  W( 71 ) = W( 71 ) + a*JVS( 222 )
  W( 75 ) = W( 75 ) + a*JVS( 223 )
  a = -W( 46 ) / JVS(          229  )
  W( 46 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 230 )
  W( 56 ) = W( 56 ) + a*JVS( 231 )
  W( 58 ) = W( 58 ) + a*JVS( 232 )
  W( 62 ) = W( 62 ) + a*JVS( 233 )
  W( 67 ) = W( 67 ) + a*JVS( 234 )
  W( 71 ) = W( 71 ) + a*JVS( 235 )
  W( 74 ) = W( 74 ) + a*JVS( 236 )
  W( 75 ) = W( 75 ) + a*JVS( 237 )
  a = -W( 47 ) / JVS(          244  )
  W( 47 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 245 )
  W( 51 ) = W( 51 ) + a*JVS( 246 )
  W( 55 ) = W( 55 ) + a*JVS( 247 )
  W( 57 ) = W( 57 ) + a*JVS( 248 )
  W( 58 ) = W( 58 ) + a*JVS( 249 )
  W( 59 ) = W( 59 ) + a*JVS( 250 )
  W( 60 ) = W( 60 ) + a*JVS( 251 )
  W( 61 ) = W( 61 ) + a*JVS( 252 )
  W( 64 ) = W( 64 ) + a*JVS( 253 )
  W( 67 ) = W( 67 ) + a*JVS( 254 )
  W( 71 ) = W( 71 ) + a*JVS( 255 )
  W( 74 ) = W( 74 ) + a*JVS( 256 )
  W( 75 ) = W( 75 ) + a*JVS( 257 )
  W( 78 ) = W( 78 ) + a*JVS( 258 )
  a = -W( 48 ) / JVS(          259  )
  W( 48 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 260 )
  W( 67 ) = W( 67 ) + a*JVS( 261 )
  W( 71 ) = W( 71 ) + a*JVS( 262 )
  W( 75 ) = W( 75 ) + a*JVS( 263 )
  a = -W( 49 ) / JVS(          272  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 273 )
  W( 54 ) = W( 54 ) + a*JVS( 274 )
  W( 57 ) = W( 57 ) + a*JVS( 275 )
  W( 63 ) = W( 63 ) + a*JVS( 276 )
  W( 67 ) = W( 67 ) + a*JVS( 277 )
  W( 71 ) = W( 71 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 78 ) = W( 78 ) + a*JVS( 280 )
  a = -W( 50 ) / JVS(          281  )
  W( 50 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 282 )
  W( 67 ) = W( 67 ) + a*JVS( 283 )
  W( 71 ) = W( 71 ) + a*JVS( 284 )
  W( 75 ) = W( 75 ) + a*JVS( 285 )
  a = -W( 51 ) / JVS(          288  )
  W( 51 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 289 )
  W( 70 ) = W( 70 ) + a*JVS( 290 )
  W( 71 ) = W( 71 ) + a*JVS( 291 )
  W( 72 ) = W( 72 ) + a*JVS( 292 )
  W( 73 ) = W( 73 ) + a*JVS( 293 )
  W( 74 ) = W( 74 ) + a*JVS( 294 )
  W( 75 ) = W( 75 ) + a*JVS( 295 )
  W( 78 ) = W( 78 ) + a*JVS( 296 )
  W( 79 ) = W( 79 ) + a*JVS( 297 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 55 ) / JVS(          338  )
  W( 55 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 339 )
  W( 63 ) = W( 63 ) + a*JVS( 340 )
  W( 67 ) = W( 67 ) + a*JVS( 341 )
  W( 71 ) = W( 71 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 60 ) / JVS(          378  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 379 )
  W( 63 ) = W( 63 ) + a*JVS( 380 )
  W( 64 ) = W( 64 ) + a*JVS( 381 )
  W( 65 ) = W( 65 ) + a*JVS( 382 )
  W( 66 ) = W( 66 ) + a*JVS( 383 )
  W( 67 ) = W( 67 ) + a*JVS( 384 )
  W( 68 ) = W( 68 ) + a*JVS( 385 )
  W( 70 ) = W( 70 ) + a*JVS( 386 )
  W( 71 ) = W( 71 ) + a*JVS( 387 )
  W( 72 ) = W( 72 ) + a*JVS( 388 )
  W( 73 ) = W( 73 ) + a*JVS( 389 )
  W( 75 ) = W( 75 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 61 ) / JVS(          412  )
  W( 61 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 413 )
  W( 63 ) = W( 63 ) + a*JVS( 414 )
  W( 65 ) = W( 65 ) + a*JVS( 415 )
  W( 66 ) = W( 66 ) + a*JVS( 416 )
  W( 67 ) = W( 67 ) + a*JVS( 417 )
  W( 68 ) = W( 68 ) + a*JVS( 418 )
  W( 69 ) = W( 69 ) + a*JVS( 419 )
  W( 70 ) = W( 70 ) + a*JVS( 420 )
  W( 71 ) = W( 71 ) + a*JVS( 421 )
  W( 72 ) = W( 72 ) + a*JVS( 422 )
  W( 73 ) = W( 73 ) + a*JVS( 423 )
  W( 74 ) = W( 74 ) + a*JVS( 424 )
  W( 75 ) = W( 75 ) + a*JVS( 425 )
  W( 76 ) = W( 76 ) + a*JVS( 426 )
  W( 77 ) = W( 77 ) + a*JVS( 427 )
  W( 78 ) = W( 78 ) + a*JVS( 428 )
  W( 79 ) = W( 79 ) + a*JVS( 429 )
  a = -W( 62 ) / JVS(          435  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 436 )
  W( 67 ) = W( 67 ) + a*JVS( 437 )
  W( 69 ) = W( 69 ) + a*JVS( 438 )
  W( 71 ) = W( 71 ) + a*JVS( 439 )
  W( 73 ) = W( 73 ) + a*JVS( 440 )
  W( 74 ) = W( 74 ) + a*JVS( 441 )
  W( 75 ) = W( 75 ) + a*JVS( 442 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  a = -W( 64 ) / JVS(          476  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 477 )
  W( 66 ) = W( 66 ) + a*JVS( 478 )
  W( 67 ) = W( 67 ) + a*JVS( 479 )
  W( 69 ) = W( 69 ) + a*JVS( 480 )
  W( 71 ) = W( 71 ) + a*JVS( 481 )
  W( 73 ) = W( 73 ) + a*JVS( 482 )
  W( 74 ) = W( 74 ) + a*JVS( 483 )
  W( 75 ) = W( 75 ) + a*JVS( 484 )
  W( 77 ) = W( 77 ) + a*JVS( 485 )
  W( 78 ) = W( 78 ) + a*JVS( 486 )
  a = -W( 65 ) / JVS(          498  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 499 )
  W( 67 ) = W( 67 ) + a*JVS( 500 )
  W( 69 ) = W( 69 ) + a*JVS( 501 )
  W( 71 ) = W( 71 ) + a*JVS( 502 )
  W( 73 ) = W( 73 ) + a*JVS( 503 )
  W( 74 ) = W( 74 ) + a*JVS( 504 )
  W( 75 ) = W( 75 ) + a*JVS( 505 )
  W( 76 ) = W( 76 ) + a*JVS( 506 )
  W( 77 ) = W( 77 ) + a*JVS( 507 )
  a = -W( 66 ) / JVS(          519  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 520 )
  W( 68 ) = W( 68 ) + a*JVS( 521 )
  W( 69 ) = W( 69 ) + a*JVS( 522 )
  W( 71 ) = W( 71 ) + a*JVS( 523 )
  W( 72 ) = W( 72 ) + a*JVS( 524 )
  W( 73 ) = W( 73 ) + a*JVS( 525 )
  W( 74 ) = W( 74 ) + a*JVS( 526 )
  W( 75 ) = W( 75 ) + a*JVS( 527 )
  W( 76 ) = W( 76 ) + a*JVS( 528 )
  W( 77 ) = W( 77 ) + a*JVS( 529 )
  W( 79 ) = W( 79 ) + a*JVS( 530 )
  a = -W( 67 ) / JVS(          543  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 544 )
  W( 70 ) = W( 70 ) + a*JVS( 545 )
  W( 71 ) = W( 71 ) + a*JVS( 546 )
  W( 72 ) = W( 72 ) + a*JVS( 547 )
  W( 73 ) = W( 73 ) + a*JVS( 548 )
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 78 ) = W( 78 ) + a*JVS( 551 )
  W( 79 ) = W( 79 ) + a*JVS( 552 )
  a = -W( 68 ) / JVS(          574  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 575 )
  W( 70 ) = W( 70 ) + a*JVS( 576 )
  W( 71 ) = W( 71 ) + a*JVS( 577 )
  W( 72 ) = W( 72 ) + a*JVS( 578 )
  W( 73 ) = W( 73 ) + a*JVS( 579 )
  W( 74 ) = W( 74 ) + a*JVS( 580 )
  W( 75 ) = W( 75 ) + a*JVS( 581 )
  W( 76 ) = W( 76 ) + a*JVS( 582 )
  W( 77 ) = W( 77 ) + a*JVS( 583 )
  W( 78 ) = W( 78 ) + a*JVS( 584 )
  W( 79 ) = W( 79 ) + a*JVS( 585 )
  a = -W( 69 ) / JVS(          606  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 607 )
  W( 71 ) = W( 71 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  W( 74 ) = W( 74 ) + a*JVS( 611 )
  W( 75 ) = W( 75 ) + a*JVS( 612 )
  W( 76 ) = W( 76 ) + a*JVS( 613 )
  W( 77 ) = W( 77 ) + a*JVS( 614 )
  W( 78 ) = W( 78 ) + a*JVS( 615 )
  W( 79 ) = W( 79 ) + a*JVS( 616 )
  a = -W( 70 ) / JVS(          627  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 628 )
  W( 72 ) = W( 72 ) + a*JVS( 629 )
  W( 73 ) = W( 73 ) + a*JVS( 630 )
  W( 74 ) = W( 74 ) + a*JVS( 631 )
  W( 75 ) = W( 75 ) + a*JVS( 632 )
  W( 76 ) = W( 76 ) + a*JVS( 633 )
  W( 77 ) = W( 77 ) + a*JVS( 634 )
  W( 78 ) = W( 78 ) + a*JVS( 635 )
  W( 79 ) = W( 79 ) + a*JVS( 636 )
  JVS( 637) = W( 10 )
  JVS( 638) = W( 11 )
  JVS( 639) = W( 12 )
  JVS( 640) = W( 13 )
  JVS( 641) = W( 19 )
  JVS( 642) = W( 20 )
  JVS( 643) = W( 21 )
  JVS( 644) = W( 23 )
  JVS( 645) = W( 24 )
  JVS( 646) = W( 26 )
  JVS( 647) = W( 27 )
  JVS( 648) = W( 28 )
  JVS( 649) = W( 29 )
  JVS( 650) = W( 32 )
  JVS( 651) = W( 33 )
  JVS( 652) = W( 34 )
  JVS( 653) = W( 35 )
  JVS( 654) = W( 36 )
  JVS( 655) = W( 37 )
  JVS( 656) = W( 38 )
  JVS( 657) = W( 39 )
  JVS( 658) = W( 41 )
  JVS( 659) = W( 42 )
  JVS( 660) = W( 43 )
  JVS( 661) = W( 44 )
  JVS( 662) = W( 45 )
  JVS( 663) = W( 46 )
  JVS( 664) = W( 47 )
  JVS( 665) = W( 48 )
  JVS( 666) = W( 49 )
  JVS( 667) = W( 50 )
  JVS( 668) = W( 51 )
  JVS( 669) = W( 52 )
  JVS( 670) = W( 54 )
  JVS( 671) = W( 55 )
  JVS( 672) = W( 56 )
  JVS( 673) = W( 57 )
  JVS( 674) = W( 58 )
  JVS( 675) = W( 59 )
  JVS( 676) = W( 60 )
  JVS( 677) = W( 61 )
  JVS( 678) = W( 62 )
  JVS( 679) = W( 63 )
  JVS( 680) = W( 64 )
  JVS( 681) = W( 65 )
  JVS( 682) = W( 66 )
  JVS( 683) = W( 67 )
  JVS( 684) = W( 68 )
  JVS( 685) = W( 69 )
  JVS( 686) = W( 70 )
  JVS( 687) = W( 71 )
  JVS( 688) = W( 72 )
  JVS( 689) = W( 73 )
  JVS( 690) = W( 74 )
  JVS( 691) = W( 75 )
  JVS( 692) = W( 76 )
  JVS( 693) = W( 77 )
  JVS( 694) = W( 78 )
  JVS( 695) = W( 79 )
  IF ( ABS(  JVS( 713 )) < TINY(a) ) THEN
         IER = 72                                      
         RETURN
  END IF
   W( 16 ) = JVS( 696 )
   W( 49 ) = JVS( 697 )
   W( 51 ) = JVS( 698 )
   W( 54 ) = JVS( 699 )
   W( 55 ) = JVS( 700 )
   W( 57 ) = JVS( 701 )
   W( 58 ) = JVS( 702 )
   W( 59 ) = JVS( 703 )
   W( 63 ) = JVS( 704 )
   W( 64 ) = JVS( 705 )
   W( 65 ) = JVS( 706 )
   W( 66 ) = JVS( 707 )
   W( 67 ) = JVS( 708 )
   W( 68 ) = JVS( 709 )
   W( 69 ) = JVS( 710 )
   W( 70 ) = JVS( 711 )
   W( 71 ) = JVS( 712 )
   W( 72 ) = JVS( 713 )
   W( 73 ) = JVS( 714 )
   W( 74 ) = JVS( 715 )
   W( 75 ) = JVS( 716 )
   W( 76 ) = JVS( 717 )
   W( 77 ) = JVS( 718 )
   W( 78 ) = JVS( 719 )
   W( 79 ) = JVS( 720 )
  a = -W( 16 ) / JVS(           85  )
  W( 16 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 86 )
  W( 74 ) = W( 74 ) + a*JVS( 87 )
  a = -W( 49 ) / JVS(          272  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 273 )
  W( 54 ) = W( 54 ) + a*JVS( 274 )
  W( 57 ) = W( 57 ) + a*JVS( 275 )
  W( 63 ) = W( 63 ) + a*JVS( 276 )
  W( 67 ) = W( 67 ) + a*JVS( 277 )
  W( 71 ) = W( 71 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 78 ) = W( 78 ) + a*JVS( 280 )
  a = -W( 51 ) / JVS(          288  )
  W( 51 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 289 )
  W( 70 ) = W( 70 ) + a*JVS( 290 )
  W( 71 ) = W( 71 ) + a*JVS( 291 )
  W( 72 ) = W( 72 ) + a*JVS( 292 )
  W( 73 ) = W( 73 ) + a*JVS( 293 )
  W( 74 ) = W( 74 ) + a*JVS( 294 )
  W( 75 ) = W( 75 ) + a*JVS( 295 )
  W( 78 ) = W( 78 ) + a*JVS( 296 )
  W( 79 ) = W( 79 ) + a*JVS( 297 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 55 ) / JVS(          338  )
  W( 55 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 339 )
  W( 63 ) = W( 63 ) + a*JVS( 340 )
  W( 67 ) = W( 67 ) + a*JVS( 341 )
  W( 71 ) = W( 71 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  a = -W( 64 ) / JVS(          476  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 477 )
  W( 66 ) = W( 66 ) + a*JVS( 478 )
  W( 67 ) = W( 67 ) + a*JVS( 479 )
  W( 69 ) = W( 69 ) + a*JVS( 480 )
  W( 71 ) = W( 71 ) + a*JVS( 481 )
  W( 73 ) = W( 73 ) + a*JVS( 482 )
  W( 74 ) = W( 74 ) + a*JVS( 483 )
  W( 75 ) = W( 75 ) + a*JVS( 484 )
  W( 77 ) = W( 77 ) + a*JVS( 485 )
  W( 78 ) = W( 78 ) + a*JVS( 486 )
  a = -W( 65 ) / JVS(          498  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 499 )
  W( 67 ) = W( 67 ) + a*JVS( 500 )
  W( 69 ) = W( 69 ) + a*JVS( 501 )
  W( 71 ) = W( 71 ) + a*JVS( 502 )
  W( 73 ) = W( 73 ) + a*JVS( 503 )
  W( 74 ) = W( 74 ) + a*JVS( 504 )
  W( 75 ) = W( 75 ) + a*JVS( 505 )
  W( 76 ) = W( 76 ) + a*JVS( 506 )
  W( 77 ) = W( 77 ) + a*JVS( 507 )
  a = -W( 66 ) / JVS(          519  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 520 )
  W( 68 ) = W( 68 ) + a*JVS( 521 )
  W( 69 ) = W( 69 ) + a*JVS( 522 )
  W( 71 ) = W( 71 ) + a*JVS( 523 )
  W( 72 ) = W( 72 ) + a*JVS( 524 )
  W( 73 ) = W( 73 ) + a*JVS( 525 )
  W( 74 ) = W( 74 ) + a*JVS( 526 )
  W( 75 ) = W( 75 ) + a*JVS( 527 )
  W( 76 ) = W( 76 ) + a*JVS( 528 )
  W( 77 ) = W( 77 ) + a*JVS( 529 )
  W( 79 ) = W( 79 ) + a*JVS( 530 )
  a = -W( 67 ) / JVS(          543  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 544 )
  W( 70 ) = W( 70 ) + a*JVS( 545 )
  W( 71 ) = W( 71 ) + a*JVS( 546 )
  W( 72 ) = W( 72 ) + a*JVS( 547 )
  W( 73 ) = W( 73 ) + a*JVS( 548 )
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 78 ) = W( 78 ) + a*JVS( 551 )
  W( 79 ) = W( 79 ) + a*JVS( 552 )
  a = -W( 68 ) / JVS(          574  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 575 )
  W( 70 ) = W( 70 ) + a*JVS( 576 )
  W( 71 ) = W( 71 ) + a*JVS( 577 )
  W( 72 ) = W( 72 ) + a*JVS( 578 )
  W( 73 ) = W( 73 ) + a*JVS( 579 )
  W( 74 ) = W( 74 ) + a*JVS( 580 )
  W( 75 ) = W( 75 ) + a*JVS( 581 )
  W( 76 ) = W( 76 ) + a*JVS( 582 )
  W( 77 ) = W( 77 ) + a*JVS( 583 )
  W( 78 ) = W( 78 ) + a*JVS( 584 )
  W( 79 ) = W( 79 ) + a*JVS( 585 )
  a = -W( 69 ) / JVS(          606  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 607 )
  W( 71 ) = W( 71 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  W( 74 ) = W( 74 ) + a*JVS( 611 )
  W( 75 ) = W( 75 ) + a*JVS( 612 )
  W( 76 ) = W( 76 ) + a*JVS( 613 )
  W( 77 ) = W( 77 ) + a*JVS( 614 )
  W( 78 ) = W( 78 ) + a*JVS( 615 )
  W( 79 ) = W( 79 ) + a*JVS( 616 )
  a = -W( 70 ) / JVS(          627  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 628 )
  W( 72 ) = W( 72 ) + a*JVS( 629 )
  W( 73 ) = W( 73 ) + a*JVS( 630 )
  W( 74 ) = W( 74 ) + a*JVS( 631 )
  W( 75 ) = W( 75 ) + a*JVS( 632 )
  W( 76 ) = W( 76 ) + a*JVS( 633 )
  W( 77 ) = W( 77 ) + a*JVS( 634 )
  W( 78 ) = W( 78 ) + a*JVS( 635 )
  W( 79 ) = W( 79 ) + a*JVS( 636 )
  a = -W( 71 ) / JVS(          687  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 688 )
  W( 73 ) = W( 73 ) + a*JVS( 689 )
  W( 74 ) = W( 74 ) + a*JVS( 690 )
  W( 75 ) = W( 75 ) + a*JVS( 691 )
  W( 76 ) = W( 76 ) + a*JVS( 692 )
  W( 77 ) = W( 77 ) + a*JVS( 693 )
  W( 78 ) = W( 78 ) + a*JVS( 694 )
  W( 79 ) = W( 79 ) + a*JVS( 695 )
  JVS( 696) = W( 16 )
  JVS( 697) = W( 49 )
  JVS( 698) = W( 51 )
  JVS( 699) = W( 54 )
  JVS( 700) = W( 55 )
  JVS( 701) = W( 57 )
  JVS( 702) = W( 58 )
  JVS( 703) = W( 59 )
  JVS( 704) = W( 63 )
  JVS( 705) = W( 64 )
  JVS( 706) = W( 65 )
  JVS( 707) = W( 66 )
  JVS( 708) = W( 67 )
  JVS( 709) = W( 68 )
  JVS( 710) = W( 69 )
  JVS( 711) = W( 70 )
  JVS( 712) = W( 71 )
  JVS( 713) = W( 72 )
  JVS( 714) = W( 73 )
  JVS( 715) = W( 74 )
  JVS( 716) = W( 75 )
  JVS( 717) = W( 76 )
  JVS( 718) = W( 77 )
  JVS( 719) = W( 78 )
  JVS( 720) = W( 79 )
  IF ( ABS(  JVS( 739 )) < TINY(a) ) THEN
         IER = 73                                      
         RETURN
  END IF
   W( 23 ) = JVS( 721 )
   W( 30 ) = JVS( 722 )
   W( 53 ) = JVS( 723 )
   W( 54 ) = JVS( 724 )
   W( 56 ) = JVS( 725 )
   W( 58 ) = JVS( 726 )
   W( 59 ) = JVS( 727 )
   W( 61 ) = JVS( 728 )
   W( 62 ) = JVS( 729 )
   W( 63 ) = JVS( 730 )
   W( 65 ) = JVS( 731 )
   W( 66 ) = JVS( 732 )
   W( 67 ) = JVS( 733 )
   W( 68 ) = JVS( 734 )
   W( 69 ) = JVS( 735 )
   W( 70 ) = JVS( 736 )
   W( 71 ) = JVS( 737 )
   W( 72 ) = JVS( 738 )
   W( 73 ) = JVS( 739 )
   W( 74 ) = JVS( 740 )
   W( 75 ) = JVS( 741 )
   W( 76 ) = JVS( 742 )
   W( 77 ) = JVS( 743 )
   W( 78 ) = JVS( 744 )
   W( 79 ) = JVS( 745 )
  a = -W( 23 ) / JVS(          104  )
  W( 23 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 105 )
  W( 73 ) = W( 73 ) + a*JVS( 106 )
  a = -W( 30 ) / JVS(          126  )
  W( 30 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 127 )
  W( 73 ) = W( 73 ) + a*JVS( 128 )
  W( 78 ) = W( 78 ) + a*JVS( 129 )
  a = -W( 53 ) / JVS(          310  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 311 )
  W( 56 ) = W( 56 ) + a*JVS( 312 )
  W( 58 ) = W( 58 ) + a*JVS( 313 )
  W( 59 ) = W( 59 ) + a*JVS( 314 )
  W( 62 ) = W( 62 ) + a*JVS( 315 )
  W( 63 ) = W( 63 ) + a*JVS( 316 )
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 69 ) = W( 69 ) + a*JVS( 321 )
  W( 70 ) = W( 70 ) + a*JVS( 322 )
  W( 71 ) = W( 71 ) + a*JVS( 323 )
  W( 72 ) = W( 72 ) + a*JVS( 324 )
  W( 73 ) = W( 73 ) + a*JVS( 325 )
  W( 74 ) = W( 74 ) + a*JVS( 326 )
  W( 75 ) = W( 75 ) + a*JVS( 327 )
  W( 76 ) = W( 76 ) + a*JVS( 328 )
  W( 77 ) = W( 77 ) + a*JVS( 329 )
  W( 78 ) = W( 78 ) + a*JVS( 330 )
  W( 79 ) = W( 79 ) + a*JVS( 331 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 61 ) / JVS(          412  )
  W( 61 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 413 )
  W( 63 ) = W( 63 ) + a*JVS( 414 )
  W( 65 ) = W( 65 ) + a*JVS( 415 )
  W( 66 ) = W( 66 ) + a*JVS( 416 )
  W( 67 ) = W( 67 ) + a*JVS( 417 )
  W( 68 ) = W( 68 ) + a*JVS( 418 )
  W( 69 ) = W( 69 ) + a*JVS( 419 )
  W( 70 ) = W( 70 ) + a*JVS( 420 )
  W( 71 ) = W( 71 ) + a*JVS( 421 )
  W( 72 ) = W( 72 ) + a*JVS( 422 )
  W( 73 ) = W( 73 ) + a*JVS( 423 )
  W( 74 ) = W( 74 ) + a*JVS( 424 )
  W( 75 ) = W( 75 ) + a*JVS( 425 )
  W( 76 ) = W( 76 ) + a*JVS( 426 )
  W( 77 ) = W( 77 ) + a*JVS( 427 )
  W( 78 ) = W( 78 ) + a*JVS( 428 )
  W( 79 ) = W( 79 ) + a*JVS( 429 )
  a = -W( 62 ) / JVS(          435  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 436 )
  W( 67 ) = W( 67 ) + a*JVS( 437 )
  W( 69 ) = W( 69 ) + a*JVS( 438 )
  W( 71 ) = W( 71 ) + a*JVS( 439 )
  W( 73 ) = W( 73 ) + a*JVS( 440 )
  W( 74 ) = W( 74 ) + a*JVS( 441 )
  W( 75 ) = W( 75 ) + a*JVS( 442 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  a = -W( 65 ) / JVS(          498  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 499 )
  W( 67 ) = W( 67 ) + a*JVS( 500 )
  W( 69 ) = W( 69 ) + a*JVS( 501 )
  W( 71 ) = W( 71 ) + a*JVS( 502 )
  W( 73 ) = W( 73 ) + a*JVS( 503 )
  W( 74 ) = W( 74 ) + a*JVS( 504 )
  W( 75 ) = W( 75 ) + a*JVS( 505 )
  W( 76 ) = W( 76 ) + a*JVS( 506 )
  W( 77 ) = W( 77 ) + a*JVS( 507 )
  a = -W( 66 ) / JVS(          519  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 520 )
  W( 68 ) = W( 68 ) + a*JVS( 521 )
  W( 69 ) = W( 69 ) + a*JVS( 522 )
  W( 71 ) = W( 71 ) + a*JVS( 523 )
  W( 72 ) = W( 72 ) + a*JVS( 524 )
  W( 73 ) = W( 73 ) + a*JVS( 525 )
  W( 74 ) = W( 74 ) + a*JVS( 526 )
  W( 75 ) = W( 75 ) + a*JVS( 527 )
  W( 76 ) = W( 76 ) + a*JVS( 528 )
  W( 77 ) = W( 77 ) + a*JVS( 529 )
  W( 79 ) = W( 79 ) + a*JVS( 530 )
  a = -W( 67 ) / JVS(          543  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 544 )
  W( 70 ) = W( 70 ) + a*JVS( 545 )
  W( 71 ) = W( 71 ) + a*JVS( 546 )
  W( 72 ) = W( 72 ) + a*JVS( 547 )
  W( 73 ) = W( 73 ) + a*JVS( 548 )
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 78 ) = W( 78 ) + a*JVS( 551 )
  W( 79 ) = W( 79 ) + a*JVS( 552 )
  a = -W( 68 ) / JVS(          574  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 575 )
  W( 70 ) = W( 70 ) + a*JVS( 576 )
  W( 71 ) = W( 71 ) + a*JVS( 577 )
  W( 72 ) = W( 72 ) + a*JVS( 578 )
  W( 73 ) = W( 73 ) + a*JVS( 579 )
  W( 74 ) = W( 74 ) + a*JVS( 580 )
  W( 75 ) = W( 75 ) + a*JVS( 581 )
  W( 76 ) = W( 76 ) + a*JVS( 582 )
  W( 77 ) = W( 77 ) + a*JVS( 583 )
  W( 78 ) = W( 78 ) + a*JVS( 584 )
  W( 79 ) = W( 79 ) + a*JVS( 585 )
  a = -W( 69 ) / JVS(          606  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 607 )
  W( 71 ) = W( 71 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  W( 74 ) = W( 74 ) + a*JVS( 611 )
  W( 75 ) = W( 75 ) + a*JVS( 612 )
  W( 76 ) = W( 76 ) + a*JVS( 613 )
  W( 77 ) = W( 77 ) + a*JVS( 614 )
  W( 78 ) = W( 78 ) + a*JVS( 615 )
  W( 79 ) = W( 79 ) + a*JVS( 616 )
  a = -W( 70 ) / JVS(          627  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 628 )
  W( 72 ) = W( 72 ) + a*JVS( 629 )
  W( 73 ) = W( 73 ) + a*JVS( 630 )
  W( 74 ) = W( 74 ) + a*JVS( 631 )
  W( 75 ) = W( 75 ) + a*JVS( 632 )
  W( 76 ) = W( 76 ) + a*JVS( 633 )
  W( 77 ) = W( 77 ) + a*JVS( 634 )
  W( 78 ) = W( 78 ) + a*JVS( 635 )
  W( 79 ) = W( 79 ) + a*JVS( 636 )
  a = -W( 71 ) / JVS(          687  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 688 )
  W( 73 ) = W( 73 ) + a*JVS( 689 )
  W( 74 ) = W( 74 ) + a*JVS( 690 )
  W( 75 ) = W( 75 ) + a*JVS( 691 )
  W( 76 ) = W( 76 ) + a*JVS( 692 )
  W( 77 ) = W( 77 ) + a*JVS( 693 )
  W( 78 ) = W( 78 ) + a*JVS( 694 )
  W( 79 ) = W( 79 ) + a*JVS( 695 )
  a = -W( 72 ) / JVS(          713  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 714 )
  W( 74 ) = W( 74 ) + a*JVS( 715 )
  W( 75 ) = W( 75 ) + a*JVS( 716 )
  W( 76 ) = W( 76 ) + a*JVS( 717 )
  W( 77 ) = W( 77 ) + a*JVS( 718 )
  W( 78 ) = W( 78 ) + a*JVS( 719 )
  W( 79 ) = W( 79 ) + a*JVS( 720 )
  JVS( 721) = W( 23 )
  JVS( 722) = W( 30 )
  JVS( 723) = W( 53 )
  JVS( 724) = W( 54 )
  JVS( 725) = W( 56 )
  JVS( 726) = W( 58 )
  JVS( 727) = W( 59 )
  JVS( 728) = W( 61 )
  JVS( 729) = W( 62 )
  JVS( 730) = W( 63 )
  JVS( 731) = W( 65 )
  JVS( 732) = W( 66 )
  JVS( 733) = W( 67 )
  JVS( 734) = W( 68 )
  JVS( 735) = W( 69 )
  JVS( 736) = W( 70 )
  JVS( 737) = W( 71 )
  JVS( 738) = W( 72 )
  JVS( 739) = W( 73 )
  JVS( 740) = W( 74 )
  JVS( 741) = W( 75 )
  JVS( 742) = W( 76 )
  JVS( 743) = W( 77 )
  JVS( 744) = W( 78 )
  JVS( 745) = W( 79 )
  IF ( ABS(  JVS( 782 )) < TINY(a) ) THEN
         IER = 74                                      
         RETURN
  END IF
   W( 15 ) = JVS( 746 )
   W( 16 ) = JVS( 747 )
   W( 17 ) = JVS( 748 )
   W( 18 ) = JVS( 749 )
   W( 22 ) = JVS( 750 )
   W( 23 ) = JVS( 751 )
   W( 25 ) = JVS( 752 )
   W( 30 ) = JVS( 753 )
   W( 31 ) = JVS( 754 )
   W( 32 ) = JVS( 755 )
   W( 40 ) = JVS( 756 )
   W( 47 ) = JVS( 757 )
   W( 49 ) = JVS( 758 )
   W( 51 ) = JVS( 759 )
   W( 52 ) = JVS( 760 )
   W( 53 ) = JVS( 761 )
   W( 54 ) = JVS( 762 )
   W( 55 ) = JVS( 763 )
   W( 56 ) = JVS( 764 )
   W( 57 ) = JVS( 765 )
   W( 58 ) = JVS( 766 )
   W( 59 ) = JVS( 767 )
   W( 60 ) = JVS( 768 )
   W( 61 ) = JVS( 769 )
   W( 62 ) = JVS( 770 )
   W( 63 ) = JVS( 771 )
   W( 64 ) = JVS( 772 )
   W( 65 ) = JVS( 773 )
   W( 66 ) = JVS( 774 )
   W( 67 ) = JVS( 775 )
   W( 68 ) = JVS( 776 )
   W( 69 ) = JVS( 777 )
   W( 70 ) = JVS( 778 )
   W( 71 ) = JVS( 779 )
   W( 72 ) = JVS( 780 )
   W( 73 ) = JVS( 781 )
   W( 74 ) = JVS( 782 )
   W( 75 ) = JVS( 783 )
   W( 76 ) = JVS( 784 )
   W( 77 ) = JVS( 785 )
   W( 78 ) = JVS( 786 )
   W( 79 ) = JVS( 787 )
  a = -W( 15 ) / JVS(           82  )
  W( 15 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 83 )
  W( 74 ) = W( 74 ) + a*JVS( 84 )
  a = -W( 16 ) / JVS(           85  )
  W( 16 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 86 )
  W( 74 ) = W( 74 ) + a*JVS( 87 )
  a = -W( 17 ) / JVS(           88  )
  W( 17 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 89 )
  W( 79 ) = W( 79 ) + a*JVS( 90 )
  a = -W( 18 ) / JVS(           91  )
  W( 18 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 92 )
  W( 74 ) = W( 74 ) + a*JVS( 93 )
  a = -W( 22 ) / JVS(          101  )
  W( 22 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 102 )
  W( 75 ) = W( 75 ) + a*JVS( 103 )
  a = -W( 23 ) / JVS(          104  )
  W( 23 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 105 )
  W( 73 ) = W( 73 ) + a*JVS( 106 )
  a = -W( 25 ) / JVS(          110  )
  W( 25 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 111 )
  W( 74 ) = W( 74 ) + a*JVS( 112 )
  a = -W( 30 ) / JVS(          126  )
  W( 30 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 127 )
  W( 73 ) = W( 73 ) + a*JVS( 128 )
  W( 78 ) = W( 78 ) + a*JVS( 129 )
  a = -W( 31 ) / JVS(          130  )
  W( 31 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 131 )
  W( 74 ) = W( 74 ) + a*JVS( 132 )
  W( 75 ) = W( 75 ) + a*JVS( 133 )
  W( 78 ) = W( 78 ) + a*JVS( 134 )
  a = -W( 32 ) / JVS(          135  )
  W( 32 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 136 )
  W( 74 ) = W( 74 ) + a*JVS( 137 )
  W( 78 ) = W( 78 ) + a*JVS( 138 )
  a = -W( 40 ) / JVS(          165  )
  W( 40 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 166 )
  W( 74 ) = W( 74 ) + a*JVS( 167 )
  W( 75 ) = W( 75 ) + a*JVS( 168 )
  W( 78 ) = W( 78 ) + a*JVS( 169 )
  a = -W( 47 ) / JVS(          244  )
  W( 47 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 245 )
  W( 51 ) = W( 51 ) + a*JVS( 246 )
  W( 55 ) = W( 55 ) + a*JVS( 247 )
  W( 57 ) = W( 57 ) + a*JVS( 248 )
  W( 58 ) = W( 58 ) + a*JVS( 249 )
  W( 59 ) = W( 59 ) + a*JVS( 250 )
  W( 60 ) = W( 60 ) + a*JVS( 251 )
  W( 61 ) = W( 61 ) + a*JVS( 252 )
  W( 64 ) = W( 64 ) + a*JVS( 253 )
  W( 67 ) = W( 67 ) + a*JVS( 254 )
  W( 71 ) = W( 71 ) + a*JVS( 255 )
  W( 74 ) = W( 74 ) + a*JVS( 256 )
  W( 75 ) = W( 75 ) + a*JVS( 257 )
  W( 78 ) = W( 78 ) + a*JVS( 258 )
  a = -W( 49 ) / JVS(          272  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 273 )
  W( 54 ) = W( 54 ) + a*JVS( 274 )
  W( 57 ) = W( 57 ) + a*JVS( 275 )
  W( 63 ) = W( 63 ) + a*JVS( 276 )
  W( 67 ) = W( 67 ) + a*JVS( 277 )
  W( 71 ) = W( 71 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 78 ) = W( 78 ) + a*JVS( 280 )
  a = -W( 51 ) / JVS(          288  )
  W( 51 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 289 )
  W( 70 ) = W( 70 ) + a*JVS( 290 )
  W( 71 ) = W( 71 ) + a*JVS( 291 )
  W( 72 ) = W( 72 ) + a*JVS( 292 )
  W( 73 ) = W( 73 ) + a*JVS( 293 )
  W( 74 ) = W( 74 ) + a*JVS( 294 )
  W( 75 ) = W( 75 ) + a*JVS( 295 )
  W( 78 ) = W( 78 ) + a*JVS( 296 )
  W( 79 ) = W( 79 ) + a*JVS( 297 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  a = -W( 53 ) / JVS(          310  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 311 )
  W( 56 ) = W( 56 ) + a*JVS( 312 )
  W( 58 ) = W( 58 ) + a*JVS( 313 )
  W( 59 ) = W( 59 ) + a*JVS( 314 )
  W( 62 ) = W( 62 ) + a*JVS( 315 )
  W( 63 ) = W( 63 ) + a*JVS( 316 )
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 69 ) = W( 69 ) + a*JVS( 321 )
  W( 70 ) = W( 70 ) + a*JVS( 322 )
  W( 71 ) = W( 71 ) + a*JVS( 323 )
  W( 72 ) = W( 72 ) + a*JVS( 324 )
  W( 73 ) = W( 73 ) + a*JVS( 325 )
  W( 74 ) = W( 74 ) + a*JVS( 326 )
  W( 75 ) = W( 75 ) + a*JVS( 327 )
  W( 76 ) = W( 76 ) + a*JVS( 328 )
  W( 77 ) = W( 77 ) + a*JVS( 329 )
  W( 78 ) = W( 78 ) + a*JVS( 330 )
  W( 79 ) = W( 79 ) + a*JVS( 331 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 55 ) / JVS(          338  )
  W( 55 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 339 )
  W( 63 ) = W( 63 ) + a*JVS( 340 )
  W( 67 ) = W( 67 ) + a*JVS( 341 )
  W( 71 ) = W( 71 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 60 ) / JVS(          378  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 379 )
  W( 63 ) = W( 63 ) + a*JVS( 380 )
  W( 64 ) = W( 64 ) + a*JVS( 381 )
  W( 65 ) = W( 65 ) + a*JVS( 382 )
  W( 66 ) = W( 66 ) + a*JVS( 383 )
  W( 67 ) = W( 67 ) + a*JVS( 384 )
  W( 68 ) = W( 68 ) + a*JVS( 385 )
  W( 70 ) = W( 70 ) + a*JVS( 386 )
  W( 71 ) = W( 71 ) + a*JVS( 387 )
  W( 72 ) = W( 72 ) + a*JVS( 388 )
  W( 73 ) = W( 73 ) + a*JVS( 389 )
  W( 75 ) = W( 75 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 61 ) / JVS(          412  )
  W( 61 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 413 )
  W( 63 ) = W( 63 ) + a*JVS( 414 )
  W( 65 ) = W( 65 ) + a*JVS( 415 )
  W( 66 ) = W( 66 ) + a*JVS( 416 )
  W( 67 ) = W( 67 ) + a*JVS( 417 )
  W( 68 ) = W( 68 ) + a*JVS( 418 )
  W( 69 ) = W( 69 ) + a*JVS( 419 )
  W( 70 ) = W( 70 ) + a*JVS( 420 )
  W( 71 ) = W( 71 ) + a*JVS( 421 )
  W( 72 ) = W( 72 ) + a*JVS( 422 )
  W( 73 ) = W( 73 ) + a*JVS( 423 )
  W( 74 ) = W( 74 ) + a*JVS( 424 )
  W( 75 ) = W( 75 ) + a*JVS( 425 )
  W( 76 ) = W( 76 ) + a*JVS( 426 )
  W( 77 ) = W( 77 ) + a*JVS( 427 )
  W( 78 ) = W( 78 ) + a*JVS( 428 )
  W( 79 ) = W( 79 ) + a*JVS( 429 )
  a = -W( 62 ) / JVS(          435  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 436 )
  W( 67 ) = W( 67 ) + a*JVS( 437 )
  W( 69 ) = W( 69 ) + a*JVS( 438 )
  W( 71 ) = W( 71 ) + a*JVS( 439 )
  W( 73 ) = W( 73 ) + a*JVS( 440 )
  W( 74 ) = W( 74 ) + a*JVS( 441 )
  W( 75 ) = W( 75 ) + a*JVS( 442 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  a = -W( 64 ) / JVS(          476  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 477 )
  W( 66 ) = W( 66 ) + a*JVS( 478 )
  W( 67 ) = W( 67 ) + a*JVS( 479 )
  W( 69 ) = W( 69 ) + a*JVS( 480 )
  W( 71 ) = W( 71 ) + a*JVS( 481 )
  W( 73 ) = W( 73 ) + a*JVS( 482 )
  W( 74 ) = W( 74 ) + a*JVS( 483 )
  W( 75 ) = W( 75 ) + a*JVS( 484 )
  W( 77 ) = W( 77 ) + a*JVS( 485 )
  W( 78 ) = W( 78 ) + a*JVS( 486 )
  a = -W( 65 ) / JVS(          498  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 499 )
  W( 67 ) = W( 67 ) + a*JVS( 500 )
  W( 69 ) = W( 69 ) + a*JVS( 501 )
  W( 71 ) = W( 71 ) + a*JVS( 502 )
  W( 73 ) = W( 73 ) + a*JVS( 503 )
  W( 74 ) = W( 74 ) + a*JVS( 504 )
  W( 75 ) = W( 75 ) + a*JVS( 505 )
  W( 76 ) = W( 76 ) + a*JVS( 506 )
  W( 77 ) = W( 77 ) + a*JVS( 507 )
  a = -W( 66 ) / JVS(          519  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 520 )
  W( 68 ) = W( 68 ) + a*JVS( 521 )
  W( 69 ) = W( 69 ) + a*JVS( 522 )
  W( 71 ) = W( 71 ) + a*JVS( 523 )
  W( 72 ) = W( 72 ) + a*JVS( 524 )
  W( 73 ) = W( 73 ) + a*JVS( 525 )
  W( 74 ) = W( 74 ) + a*JVS( 526 )
  W( 75 ) = W( 75 ) + a*JVS( 527 )
  W( 76 ) = W( 76 ) + a*JVS( 528 )
  W( 77 ) = W( 77 ) + a*JVS( 529 )
  W( 79 ) = W( 79 ) + a*JVS( 530 )
  a = -W( 67 ) / JVS(          543  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 544 )
  W( 70 ) = W( 70 ) + a*JVS( 545 )
  W( 71 ) = W( 71 ) + a*JVS( 546 )
  W( 72 ) = W( 72 ) + a*JVS( 547 )
  W( 73 ) = W( 73 ) + a*JVS( 548 )
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 78 ) = W( 78 ) + a*JVS( 551 )
  W( 79 ) = W( 79 ) + a*JVS( 552 )
  a = -W( 68 ) / JVS(          574  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 575 )
  W( 70 ) = W( 70 ) + a*JVS( 576 )
  W( 71 ) = W( 71 ) + a*JVS( 577 )
  W( 72 ) = W( 72 ) + a*JVS( 578 )
  W( 73 ) = W( 73 ) + a*JVS( 579 )
  W( 74 ) = W( 74 ) + a*JVS( 580 )
  W( 75 ) = W( 75 ) + a*JVS( 581 )
  W( 76 ) = W( 76 ) + a*JVS( 582 )
  W( 77 ) = W( 77 ) + a*JVS( 583 )
  W( 78 ) = W( 78 ) + a*JVS( 584 )
  W( 79 ) = W( 79 ) + a*JVS( 585 )
  a = -W( 69 ) / JVS(          606  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 607 )
  W( 71 ) = W( 71 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  W( 74 ) = W( 74 ) + a*JVS( 611 )
  W( 75 ) = W( 75 ) + a*JVS( 612 )
  W( 76 ) = W( 76 ) + a*JVS( 613 )
  W( 77 ) = W( 77 ) + a*JVS( 614 )
  W( 78 ) = W( 78 ) + a*JVS( 615 )
  W( 79 ) = W( 79 ) + a*JVS( 616 )
  a = -W( 70 ) / JVS(          627  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 628 )
  W( 72 ) = W( 72 ) + a*JVS( 629 )
  W( 73 ) = W( 73 ) + a*JVS( 630 )
  W( 74 ) = W( 74 ) + a*JVS( 631 )
  W( 75 ) = W( 75 ) + a*JVS( 632 )
  W( 76 ) = W( 76 ) + a*JVS( 633 )
  W( 77 ) = W( 77 ) + a*JVS( 634 )
  W( 78 ) = W( 78 ) + a*JVS( 635 )
  W( 79 ) = W( 79 ) + a*JVS( 636 )
  a = -W( 71 ) / JVS(          687  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 688 )
  W( 73 ) = W( 73 ) + a*JVS( 689 )
  W( 74 ) = W( 74 ) + a*JVS( 690 )
  W( 75 ) = W( 75 ) + a*JVS( 691 )
  W( 76 ) = W( 76 ) + a*JVS( 692 )
  W( 77 ) = W( 77 ) + a*JVS( 693 )
  W( 78 ) = W( 78 ) + a*JVS( 694 )
  W( 79 ) = W( 79 ) + a*JVS( 695 )
  a = -W( 72 ) / JVS(          713  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 714 )
  W( 74 ) = W( 74 ) + a*JVS( 715 )
  W( 75 ) = W( 75 ) + a*JVS( 716 )
  W( 76 ) = W( 76 ) + a*JVS( 717 )
  W( 77 ) = W( 77 ) + a*JVS( 718 )
  W( 78 ) = W( 78 ) + a*JVS( 719 )
  W( 79 ) = W( 79 ) + a*JVS( 720 )
  a = -W( 73 ) / JVS(          739  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 740 )
  W( 75 ) = W( 75 ) + a*JVS( 741 )
  W( 76 ) = W( 76 ) + a*JVS( 742 )
  W( 77 ) = W( 77 ) + a*JVS( 743 )
  W( 78 ) = W( 78 ) + a*JVS( 744 )
  W( 79 ) = W( 79 ) + a*JVS( 745 )
  JVS( 746) = W( 15 )
  JVS( 747) = W( 16 )
  JVS( 748) = W( 17 )
  JVS( 749) = W( 18 )
  JVS( 750) = W( 22 )
  JVS( 751) = W( 23 )
  JVS( 752) = W( 25 )
  JVS( 753) = W( 30 )
  JVS( 754) = W( 31 )
  JVS( 755) = W( 32 )
  JVS( 756) = W( 40 )
  JVS( 757) = W( 47 )
  JVS( 758) = W( 49 )
  JVS( 759) = W( 51 )
  JVS( 760) = W( 52 )
  JVS( 761) = W( 53 )
  JVS( 762) = W( 54 )
  JVS( 763) = W( 55 )
  JVS( 764) = W( 56 )
  JVS( 765) = W( 57 )
  JVS( 766) = W( 58 )
  JVS( 767) = W( 59 )
  JVS( 768) = W( 60 )
  JVS( 769) = W( 61 )
  JVS( 770) = W( 62 )
  JVS( 771) = W( 63 )
  JVS( 772) = W( 64 )
  JVS( 773) = W( 65 )
  JVS( 774) = W( 66 )
  JVS( 775) = W( 67 )
  JVS( 776) = W( 68 )
  JVS( 777) = W( 69 )
  JVS( 778) = W( 70 )
  JVS( 779) = W( 71 )
  JVS( 780) = W( 72 )
  JVS( 781) = W( 73 )
  JVS( 782) = W( 74 )
  JVS( 783) = W( 75 )
  JVS( 784) = W( 76 )
  JVS( 785) = W( 77 )
  JVS( 786) = W( 78 )
  JVS( 787) = W( 79 )
  IF ( ABS(  JVS( 823 )) < TINY(a) ) THEN
         IER = 75                                      
         RETURN
  END IF
   W( 22 ) = JVS( 788 )
   W( 32 ) = JVS( 789 )
   W( 37 ) = JVS( 790 )
   W( 40 ) = JVS( 791 )
   W( 41 ) = JVS( 792 )
   W( 43 ) = JVS( 793 )
   W( 44 ) = JVS( 794 )
   W( 47 ) = JVS( 795 )
   W( 48 ) = JVS( 796 )
   W( 49 ) = JVS( 797 )
   W( 50 ) = JVS( 798 )
   W( 51 ) = JVS( 799 )
   W( 52 ) = JVS( 800 )
   W( 53 ) = JVS( 801 )
   W( 54 ) = JVS( 802 )
   W( 55 ) = JVS( 803 )
   W( 56 ) = JVS( 804 )
   W( 57 ) = JVS( 805 )
   W( 58 ) = JVS( 806 )
   W( 59 ) = JVS( 807 )
   W( 60 ) = JVS( 808 )
   W( 61 ) = JVS( 809 )
   W( 62 ) = JVS( 810 )
   W( 63 ) = JVS( 811 )
   W( 64 ) = JVS( 812 )
   W( 65 ) = JVS( 813 )
   W( 66 ) = JVS( 814 )
   W( 67 ) = JVS( 815 )
   W( 68 ) = JVS( 816 )
   W( 69 ) = JVS( 817 )
   W( 70 ) = JVS( 818 )
   W( 71 ) = JVS( 819 )
   W( 72 ) = JVS( 820 )
   W( 73 ) = JVS( 821 )
   W( 74 ) = JVS( 822 )
   W( 75 ) = JVS( 823 )
   W( 76 ) = JVS( 824 )
   W( 77 ) = JVS( 825 )
   W( 78 ) = JVS( 826 )
   W( 79 ) = JVS( 827 )
  a = -W( 22 ) / JVS(          101  )
  W( 22 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 102 )
  W( 75 ) = W( 75 ) + a*JVS( 103 )
  a = -W( 32 ) / JVS(          135  )
  W( 32 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 136 )
  W( 74 ) = W( 74 ) + a*JVS( 137 )
  W( 78 ) = W( 78 ) + a*JVS( 138 )
  a = -W( 37 ) / JVS(          153  )
  W( 37 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 154 )
  W( 75 ) = W( 75 ) + a*JVS( 155 )
  a = -W( 40 ) / JVS(          165  )
  W( 40 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 166 )
  W( 74 ) = W( 74 ) + a*JVS( 167 )
  W( 75 ) = W( 75 ) + a*JVS( 168 )
  W( 78 ) = W( 78 ) + a*JVS( 169 )
  a = -W( 41 ) / JVS(          172  )
  W( 41 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 173 )
  W( 67 ) = W( 67 ) + a*JVS( 174 )
  W( 71 ) = W( 71 ) + a*JVS( 175 )
  W( 75 ) = W( 75 ) + a*JVS( 176 )
  a = -W( 43 ) / JVS(          183  )
  W( 43 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 184 )
  W( 71 ) = W( 71 ) + a*JVS( 185 )
  W( 75 ) = W( 75 ) + a*JVS( 186 )
  W( 78 ) = W( 78 ) + a*JVS( 187 )
  a = -W( 44 ) / JVS(          193  )
  W( 44 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  W( 57 ) = W( 57 ) + a*JVS( 195 )
  W( 59 ) = W( 59 ) + a*JVS( 196 )
  W( 67 ) = W( 67 ) + a*JVS( 197 )
  W( 71 ) = W( 71 ) + a*JVS( 198 )
  W( 75 ) = W( 75 ) + a*JVS( 199 )
  a = -W( 47 ) / JVS(          244  )
  W( 47 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 245 )
  W( 51 ) = W( 51 ) + a*JVS( 246 )
  W( 55 ) = W( 55 ) + a*JVS( 247 )
  W( 57 ) = W( 57 ) + a*JVS( 248 )
  W( 58 ) = W( 58 ) + a*JVS( 249 )
  W( 59 ) = W( 59 ) + a*JVS( 250 )
  W( 60 ) = W( 60 ) + a*JVS( 251 )
  W( 61 ) = W( 61 ) + a*JVS( 252 )
  W( 64 ) = W( 64 ) + a*JVS( 253 )
  W( 67 ) = W( 67 ) + a*JVS( 254 )
  W( 71 ) = W( 71 ) + a*JVS( 255 )
  W( 74 ) = W( 74 ) + a*JVS( 256 )
  W( 75 ) = W( 75 ) + a*JVS( 257 )
  W( 78 ) = W( 78 ) + a*JVS( 258 )
  a = -W( 48 ) / JVS(          259  )
  W( 48 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 260 )
  W( 67 ) = W( 67 ) + a*JVS( 261 )
  W( 71 ) = W( 71 ) + a*JVS( 262 )
  W( 75 ) = W( 75 ) + a*JVS( 263 )
  a = -W( 49 ) / JVS(          272  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 273 )
  W( 54 ) = W( 54 ) + a*JVS( 274 )
  W( 57 ) = W( 57 ) + a*JVS( 275 )
  W( 63 ) = W( 63 ) + a*JVS( 276 )
  W( 67 ) = W( 67 ) + a*JVS( 277 )
  W( 71 ) = W( 71 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 78 ) = W( 78 ) + a*JVS( 280 )
  a = -W( 50 ) / JVS(          281  )
  W( 50 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 282 )
  W( 67 ) = W( 67 ) + a*JVS( 283 )
  W( 71 ) = W( 71 ) + a*JVS( 284 )
  W( 75 ) = W( 75 ) + a*JVS( 285 )
  a = -W( 51 ) / JVS(          288  )
  W( 51 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 289 )
  W( 70 ) = W( 70 ) + a*JVS( 290 )
  W( 71 ) = W( 71 ) + a*JVS( 291 )
  W( 72 ) = W( 72 ) + a*JVS( 292 )
  W( 73 ) = W( 73 ) + a*JVS( 293 )
  W( 74 ) = W( 74 ) + a*JVS( 294 )
  W( 75 ) = W( 75 ) + a*JVS( 295 )
  W( 78 ) = W( 78 ) + a*JVS( 296 )
  W( 79 ) = W( 79 ) + a*JVS( 297 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  a = -W( 53 ) / JVS(          310  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 311 )
  W( 56 ) = W( 56 ) + a*JVS( 312 )
  W( 58 ) = W( 58 ) + a*JVS( 313 )
  W( 59 ) = W( 59 ) + a*JVS( 314 )
  W( 62 ) = W( 62 ) + a*JVS( 315 )
  W( 63 ) = W( 63 ) + a*JVS( 316 )
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 69 ) = W( 69 ) + a*JVS( 321 )
  W( 70 ) = W( 70 ) + a*JVS( 322 )
  W( 71 ) = W( 71 ) + a*JVS( 323 )
  W( 72 ) = W( 72 ) + a*JVS( 324 )
  W( 73 ) = W( 73 ) + a*JVS( 325 )
  W( 74 ) = W( 74 ) + a*JVS( 326 )
  W( 75 ) = W( 75 ) + a*JVS( 327 )
  W( 76 ) = W( 76 ) + a*JVS( 328 )
  W( 77 ) = W( 77 ) + a*JVS( 329 )
  W( 78 ) = W( 78 ) + a*JVS( 330 )
  W( 79 ) = W( 79 ) + a*JVS( 331 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 55 ) / JVS(          338  )
  W( 55 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 339 )
  W( 63 ) = W( 63 ) + a*JVS( 340 )
  W( 67 ) = W( 67 ) + a*JVS( 341 )
  W( 71 ) = W( 71 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 60 ) / JVS(          378  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 379 )
  W( 63 ) = W( 63 ) + a*JVS( 380 )
  W( 64 ) = W( 64 ) + a*JVS( 381 )
  W( 65 ) = W( 65 ) + a*JVS( 382 )
  W( 66 ) = W( 66 ) + a*JVS( 383 )
  W( 67 ) = W( 67 ) + a*JVS( 384 )
  W( 68 ) = W( 68 ) + a*JVS( 385 )
  W( 70 ) = W( 70 ) + a*JVS( 386 )
  W( 71 ) = W( 71 ) + a*JVS( 387 )
  W( 72 ) = W( 72 ) + a*JVS( 388 )
  W( 73 ) = W( 73 ) + a*JVS( 389 )
  W( 75 ) = W( 75 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 61 ) / JVS(          412  )
  W( 61 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 413 )
  W( 63 ) = W( 63 ) + a*JVS( 414 )
  W( 65 ) = W( 65 ) + a*JVS( 415 )
  W( 66 ) = W( 66 ) + a*JVS( 416 )
  W( 67 ) = W( 67 ) + a*JVS( 417 )
  W( 68 ) = W( 68 ) + a*JVS( 418 )
  W( 69 ) = W( 69 ) + a*JVS( 419 )
  W( 70 ) = W( 70 ) + a*JVS( 420 )
  W( 71 ) = W( 71 ) + a*JVS( 421 )
  W( 72 ) = W( 72 ) + a*JVS( 422 )
  W( 73 ) = W( 73 ) + a*JVS( 423 )
  W( 74 ) = W( 74 ) + a*JVS( 424 )
  W( 75 ) = W( 75 ) + a*JVS( 425 )
  W( 76 ) = W( 76 ) + a*JVS( 426 )
  W( 77 ) = W( 77 ) + a*JVS( 427 )
  W( 78 ) = W( 78 ) + a*JVS( 428 )
  W( 79 ) = W( 79 ) + a*JVS( 429 )
  a = -W( 62 ) / JVS(          435  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 436 )
  W( 67 ) = W( 67 ) + a*JVS( 437 )
  W( 69 ) = W( 69 ) + a*JVS( 438 )
  W( 71 ) = W( 71 ) + a*JVS( 439 )
  W( 73 ) = W( 73 ) + a*JVS( 440 )
  W( 74 ) = W( 74 ) + a*JVS( 441 )
  W( 75 ) = W( 75 ) + a*JVS( 442 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  a = -W( 64 ) / JVS(          476  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 477 )
  W( 66 ) = W( 66 ) + a*JVS( 478 )
  W( 67 ) = W( 67 ) + a*JVS( 479 )
  W( 69 ) = W( 69 ) + a*JVS( 480 )
  W( 71 ) = W( 71 ) + a*JVS( 481 )
  W( 73 ) = W( 73 ) + a*JVS( 482 )
  W( 74 ) = W( 74 ) + a*JVS( 483 )
  W( 75 ) = W( 75 ) + a*JVS( 484 )
  W( 77 ) = W( 77 ) + a*JVS( 485 )
  W( 78 ) = W( 78 ) + a*JVS( 486 )
  a = -W( 65 ) / JVS(          498  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 499 )
  W( 67 ) = W( 67 ) + a*JVS( 500 )
  W( 69 ) = W( 69 ) + a*JVS( 501 )
  W( 71 ) = W( 71 ) + a*JVS( 502 )
  W( 73 ) = W( 73 ) + a*JVS( 503 )
  W( 74 ) = W( 74 ) + a*JVS( 504 )
  W( 75 ) = W( 75 ) + a*JVS( 505 )
  W( 76 ) = W( 76 ) + a*JVS( 506 )
  W( 77 ) = W( 77 ) + a*JVS( 507 )
  a = -W( 66 ) / JVS(          519  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 520 )
  W( 68 ) = W( 68 ) + a*JVS( 521 )
  W( 69 ) = W( 69 ) + a*JVS( 522 )
  W( 71 ) = W( 71 ) + a*JVS( 523 )
  W( 72 ) = W( 72 ) + a*JVS( 524 )
  W( 73 ) = W( 73 ) + a*JVS( 525 )
  W( 74 ) = W( 74 ) + a*JVS( 526 )
  W( 75 ) = W( 75 ) + a*JVS( 527 )
  W( 76 ) = W( 76 ) + a*JVS( 528 )
  W( 77 ) = W( 77 ) + a*JVS( 529 )
  W( 79 ) = W( 79 ) + a*JVS( 530 )
  a = -W( 67 ) / JVS(          543  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 544 )
  W( 70 ) = W( 70 ) + a*JVS( 545 )
  W( 71 ) = W( 71 ) + a*JVS( 546 )
  W( 72 ) = W( 72 ) + a*JVS( 547 )
  W( 73 ) = W( 73 ) + a*JVS( 548 )
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 78 ) = W( 78 ) + a*JVS( 551 )
  W( 79 ) = W( 79 ) + a*JVS( 552 )
  a = -W( 68 ) / JVS(          574  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 575 )
  W( 70 ) = W( 70 ) + a*JVS( 576 )
  W( 71 ) = W( 71 ) + a*JVS( 577 )
  W( 72 ) = W( 72 ) + a*JVS( 578 )
  W( 73 ) = W( 73 ) + a*JVS( 579 )
  W( 74 ) = W( 74 ) + a*JVS( 580 )
  W( 75 ) = W( 75 ) + a*JVS( 581 )
  W( 76 ) = W( 76 ) + a*JVS( 582 )
  W( 77 ) = W( 77 ) + a*JVS( 583 )
  W( 78 ) = W( 78 ) + a*JVS( 584 )
  W( 79 ) = W( 79 ) + a*JVS( 585 )
  a = -W( 69 ) / JVS(          606  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 607 )
  W( 71 ) = W( 71 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  W( 74 ) = W( 74 ) + a*JVS( 611 )
  W( 75 ) = W( 75 ) + a*JVS( 612 )
  W( 76 ) = W( 76 ) + a*JVS( 613 )
  W( 77 ) = W( 77 ) + a*JVS( 614 )
  W( 78 ) = W( 78 ) + a*JVS( 615 )
  W( 79 ) = W( 79 ) + a*JVS( 616 )
  a = -W( 70 ) / JVS(          627  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 628 )
  W( 72 ) = W( 72 ) + a*JVS( 629 )
  W( 73 ) = W( 73 ) + a*JVS( 630 )
  W( 74 ) = W( 74 ) + a*JVS( 631 )
  W( 75 ) = W( 75 ) + a*JVS( 632 )
  W( 76 ) = W( 76 ) + a*JVS( 633 )
  W( 77 ) = W( 77 ) + a*JVS( 634 )
  W( 78 ) = W( 78 ) + a*JVS( 635 )
  W( 79 ) = W( 79 ) + a*JVS( 636 )
  a = -W( 71 ) / JVS(          687  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 688 )
  W( 73 ) = W( 73 ) + a*JVS( 689 )
  W( 74 ) = W( 74 ) + a*JVS( 690 )
  W( 75 ) = W( 75 ) + a*JVS( 691 )
  W( 76 ) = W( 76 ) + a*JVS( 692 )
  W( 77 ) = W( 77 ) + a*JVS( 693 )
  W( 78 ) = W( 78 ) + a*JVS( 694 )
  W( 79 ) = W( 79 ) + a*JVS( 695 )
  a = -W( 72 ) / JVS(          713  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 714 )
  W( 74 ) = W( 74 ) + a*JVS( 715 )
  W( 75 ) = W( 75 ) + a*JVS( 716 )
  W( 76 ) = W( 76 ) + a*JVS( 717 )
  W( 77 ) = W( 77 ) + a*JVS( 718 )
  W( 78 ) = W( 78 ) + a*JVS( 719 )
  W( 79 ) = W( 79 ) + a*JVS( 720 )
  a = -W( 73 ) / JVS(          739  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 740 )
  W( 75 ) = W( 75 ) + a*JVS( 741 )
  W( 76 ) = W( 76 ) + a*JVS( 742 )
  W( 77 ) = W( 77 ) + a*JVS( 743 )
  W( 78 ) = W( 78 ) + a*JVS( 744 )
  W( 79 ) = W( 79 ) + a*JVS( 745 )
  a = -W( 74 ) / JVS(          782  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 783 )
  W( 76 ) = W( 76 ) + a*JVS( 784 )
  W( 77 ) = W( 77 ) + a*JVS( 785 )
  W( 78 ) = W( 78 ) + a*JVS( 786 )
  W( 79 ) = W( 79 ) + a*JVS( 787 )
  JVS( 788) = W( 22 )
  JVS( 789) = W( 32 )
  JVS( 790) = W( 37 )
  JVS( 791) = W( 40 )
  JVS( 792) = W( 41 )
  JVS( 793) = W( 43 )
  JVS( 794) = W( 44 )
  JVS( 795) = W( 47 )
  JVS( 796) = W( 48 )
  JVS( 797) = W( 49 )
  JVS( 798) = W( 50 )
  JVS( 799) = W( 51 )
  JVS( 800) = W( 52 )
  JVS( 801) = W( 53 )
  JVS( 802) = W( 54 )
  JVS( 803) = W( 55 )
  JVS( 804) = W( 56 )
  JVS( 805) = W( 57 )
  JVS( 806) = W( 58 )
  JVS( 807) = W( 59 )
  JVS( 808) = W( 60 )
  JVS( 809) = W( 61 )
  JVS( 810) = W( 62 )
  JVS( 811) = W( 63 )
  JVS( 812) = W( 64 )
  JVS( 813) = W( 65 )
  JVS( 814) = W( 66 )
  JVS( 815) = W( 67 )
  JVS( 816) = W( 68 )
  JVS( 817) = W( 69 )
  JVS( 818) = W( 70 )
  JVS( 819) = W( 71 )
  JVS( 820) = W( 72 )
  JVS( 821) = W( 73 )
  JVS( 822) = W( 74 )
  JVS( 823) = W( 75 )
  JVS( 824) = W( 76 )
  JVS( 825) = W( 77 )
  JVS( 826) = W( 78 )
  JVS( 827) = W( 79 )
  IF ( ABS(  JVS( 855 )) < TINY(a) ) THEN
         IER = 76                                      
         RETURN
  END IF
   W( 12 ) = JVS( 828 )
   W( 25 ) = JVS( 829 )
   W( 29 ) = JVS( 830 )
   W( 33 ) = JVS( 831 )
   W( 46 ) = JVS( 832 )
   W( 48 ) = JVS( 833 )
   W( 50 ) = JVS( 834 )
   W( 52 ) = JVS( 835 )
   W( 54 ) = JVS( 836 )
   W( 56 ) = JVS( 837 )
   W( 58 ) = JVS( 838 )
   W( 59 ) = JVS( 839 )
   W( 60 ) = JVS( 840 )
   W( 62 ) = JVS( 841 )
   W( 63 ) = JVS( 842 )
   W( 64 ) = JVS( 843 )
   W( 65 ) = JVS( 844 )
   W( 66 ) = JVS( 845 )
   W( 67 ) = JVS( 846 )
   W( 68 ) = JVS( 847 )
   W( 69 ) = JVS( 848 )
   W( 70 ) = JVS( 849 )
   W( 71 ) = JVS( 850 )
   W( 72 ) = JVS( 851 )
   W( 73 ) = JVS( 852 )
   W( 74 ) = JVS( 853 )
   W( 75 ) = JVS( 854 )
   W( 76 ) = JVS( 855 )
   W( 77 ) = JVS( 856 )
   W( 78 ) = JVS( 857 )
   W( 79 ) = JVS( 858 )
  a = -W( 12 ) / JVS(           73  )
  W( 12 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 74 )
  a = -W( 25 ) / JVS(          110  )
  W( 25 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 111 )
  W( 74 ) = W( 74 ) + a*JVS( 112 )
  a = -W( 29 ) / JVS(          122  )
  W( 29 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 123 )
  W( 76 ) = W( 76 ) + a*JVS( 124 )
  W( 78 ) = W( 78 ) + a*JVS( 125 )
  a = -W( 33 ) / JVS(          139  )
  W( 33 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 140 )
  a = -W( 46 ) / JVS(          229  )
  W( 46 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 230 )
  W( 56 ) = W( 56 ) + a*JVS( 231 )
  W( 58 ) = W( 58 ) + a*JVS( 232 )
  W( 62 ) = W( 62 ) + a*JVS( 233 )
  W( 67 ) = W( 67 ) + a*JVS( 234 )
  W( 71 ) = W( 71 ) + a*JVS( 235 )
  W( 74 ) = W( 74 ) + a*JVS( 236 )
  W( 75 ) = W( 75 ) + a*JVS( 237 )
  a = -W( 48 ) / JVS(          259  )
  W( 48 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 260 )
  W( 67 ) = W( 67 ) + a*JVS( 261 )
  W( 71 ) = W( 71 ) + a*JVS( 262 )
  W( 75 ) = W( 75 ) + a*JVS( 263 )
  a = -W( 50 ) / JVS(          281  )
  W( 50 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 282 )
  W( 67 ) = W( 67 ) + a*JVS( 283 )
  W( 71 ) = W( 71 ) + a*JVS( 284 )
  W( 75 ) = W( 75 ) + a*JVS( 285 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 60 ) / JVS(          378  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 379 )
  W( 63 ) = W( 63 ) + a*JVS( 380 )
  W( 64 ) = W( 64 ) + a*JVS( 381 )
  W( 65 ) = W( 65 ) + a*JVS( 382 )
  W( 66 ) = W( 66 ) + a*JVS( 383 )
  W( 67 ) = W( 67 ) + a*JVS( 384 )
  W( 68 ) = W( 68 ) + a*JVS( 385 )
  W( 70 ) = W( 70 ) + a*JVS( 386 )
  W( 71 ) = W( 71 ) + a*JVS( 387 )
  W( 72 ) = W( 72 ) + a*JVS( 388 )
  W( 73 ) = W( 73 ) + a*JVS( 389 )
  W( 75 ) = W( 75 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 62 ) / JVS(          435  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 436 )
  W( 67 ) = W( 67 ) + a*JVS( 437 )
  W( 69 ) = W( 69 ) + a*JVS( 438 )
  W( 71 ) = W( 71 ) + a*JVS( 439 )
  W( 73 ) = W( 73 ) + a*JVS( 440 )
  W( 74 ) = W( 74 ) + a*JVS( 441 )
  W( 75 ) = W( 75 ) + a*JVS( 442 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  a = -W( 64 ) / JVS(          476  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 477 )
  W( 66 ) = W( 66 ) + a*JVS( 478 )
  W( 67 ) = W( 67 ) + a*JVS( 479 )
  W( 69 ) = W( 69 ) + a*JVS( 480 )
  W( 71 ) = W( 71 ) + a*JVS( 481 )
  W( 73 ) = W( 73 ) + a*JVS( 482 )
  W( 74 ) = W( 74 ) + a*JVS( 483 )
  W( 75 ) = W( 75 ) + a*JVS( 484 )
  W( 77 ) = W( 77 ) + a*JVS( 485 )
  W( 78 ) = W( 78 ) + a*JVS( 486 )
  a = -W( 65 ) / JVS(          498  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 499 )
  W( 67 ) = W( 67 ) + a*JVS( 500 )
  W( 69 ) = W( 69 ) + a*JVS( 501 )
  W( 71 ) = W( 71 ) + a*JVS( 502 )
  W( 73 ) = W( 73 ) + a*JVS( 503 )
  W( 74 ) = W( 74 ) + a*JVS( 504 )
  W( 75 ) = W( 75 ) + a*JVS( 505 )
  W( 76 ) = W( 76 ) + a*JVS( 506 )
  W( 77 ) = W( 77 ) + a*JVS( 507 )
  a = -W( 66 ) / JVS(          519  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 520 )
  W( 68 ) = W( 68 ) + a*JVS( 521 )
  W( 69 ) = W( 69 ) + a*JVS( 522 )
  W( 71 ) = W( 71 ) + a*JVS( 523 )
  W( 72 ) = W( 72 ) + a*JVS( 524 )
  W( 73 ) = W( 73 ) + a*JVS( 525 )
  W( 74 ) = W( 74 ) + a*JVS( 526 )
  W( 75 ) = W( 75 ) + a*JVS( 527 )
  W( 76 ) = W( 76 ) + a*JVS( 528 )
  W( 77 ) = W( 77 ) + a*JVS( 529 )
  W( 79 ) = W( 79 ) + a*JVS( 530 )
  a = -W( 67 ) / JVS(          543  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 544 )
  W( 70 ) = W( 70 ) + a*JVS( 545 )
  W( 71 ) = W( 71 ) + a*JVS( 546 )
  W( 72 ) = W( 72 ) + a*JVS( 547 )
  W( 73 ) = W( 73 ) + a*JVS( 548 )
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 78 ) = W( 78 ) + a*JVS( 551 )
  W( 79 ) = W( 79 ) + a*JVS( 552 )
  a = -W( 68 ) / JVS(          574  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 575 )
  W( 70 ) = W( 70 ) + a*JVS( 576 )
  W( 71 ) = W( 71 ) + a*JVS( 577 )
  W( 72 ) = W( 72 ) + a*JVS( 578 )
  W( 73 ) = W( 73 ) + a*JVS( 579 )
  W( 74 ) = W( 74 ) + a*JVS( 580 )
  W( 75 ) = W( 75 ) + a*JVS( 581 )
  W( 76 ) = W( 76 ) + a*JVS( 582 )
  W( 77 ) = W( 77 ) + a*JVS( 583 )
  W( 78 ) = W( 78 ) + a*JVS( 584 )
  W( 79 ) = W( 79 ) + a*JVS( 585 )
  a = -W( 69 ) / JVS(          606  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 607 )
  W( 71 ) = W( 71 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  W( 74 ) = W( 74 ) + a*JVS( 611 )
  W( 75 ) = W( 75 ) + a*JVS( 612 )
  W( 76 ) = W( 76 ) + a*JVS( 613 )
  W( 77 ) = W( 77 ) + a*JVS( 614 )
  W( 78 ) = W( 78 ) + a*JVS( 615 )
  W( 79 ) = W( 79 ) + a*JVS( 616 )
  a = -W( 70 ) / JVS(          627  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 628 )
  W( 72 ) = W( 72 ) + a*JVS( 629 )
  W( 73 ) = W( 73 ) + a*JVS( 630 )
  W( 74 ) = W( 74 ) + a*JVS( 631 )
  W( 75 ) = W( 75 ) + a*JVS( 632 )
  W( 76 ) = W( 76 ) + a*JVS( 633 )
  W( 77 ) = W( 77 ) + a*JVS( 634 )
  W( 78 ) = W( 78 ) + a*JVS( 635 )
  W( 79 ) = W( 79 ) + a*JVS( 636 )
  a = -W( 71 ) / JVS(          687  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 688 )
  W( 73 ) = W( 73 ) + a*JVS( 689 )
  W( 74 ) = W( 74 ) + a*JVS( 690 )
  W( 75 ) = W( 75 ) + a*JVS( 691 )
  W( 76 ) = W( 76 ) + a*JVS( 692 )
  W( 77 ) = W( 77 ) + a*JVS( 693 )
  W( 78 ) = W( 78 ) + a*JVS( 694 )
  W( 79 ) = W( 79 ) + a*JVS( 695 )
  a = -W( 72 ) / JVS(          713  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 714 )
  W( 74 ) = W( 74 ) + a*JVS( 715 )
  W( 75 ) = W( 75 ) + a*JVS( 716 )
  W( 76 ) = W( 76 ) + a*JVS( 717 )
  W( 77 ) = W( 77 ) + a*JVS( 718 )
  W( 78 ) = W( 78 ) + a*JVS( 719 )
  W( 79 ) = W( 79 ) + a*JVS( 720 )
  a = -W( 73 ) / JVS(          739  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 740 )
  W( 75 ) = W( 75 ) + a*JVS( 741 )
  W( 76 ) = W( 76 ) + a*JVS( 742 )
  W( 77 ) = W( 77 ) + a*JVS( 743 )
  W( 78 ) = W( 78 ) + a*JVS( 744 )
  W( 79 ) = W( 79 ) + a*JVS( 745 )
  a = -W( 74 ) / JVS(          782  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 783 )
  W( 76 ) = W( 76 ) + a*JVS( 784 )
  W( 77 ) = W( 77 ) + a*JVS( 785 )
  W( 78 ) = W( 78 ) + a*JVS( 786 )
  W( 79 ) = W( 79 ) + a*JVS( 787 )
  a = -W( 75 ) / JVS(          823  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 824 )
  W( 77 ) = W( 77 ) + a*JVS( 825 )
  W( 78 ) = W( 78 ) + a*JVS( 826 )
  W( 79 ) = W( 79 ) + a*JVS( 827 )
  JVS( 828) = W( 12 )
  JVS( 829) = W( 25 )
  JVS( 830) = W( 29 )
  JVS( 831) = W( 33 )
  JVS( 832) = W( 46 )
  JVS( 833) = W( 48 )
  JVS( 834) = W( 50 )
  JVS( 835) = W( 52 )
  JVS( 836) = W( 54 )
  JVS( 837) = W( 56 )
  JVS( 838) = W( 58 )
  JVS( 839) = W( 59 )
  JVS( 840) = W( 60 )
  JVS( 841) = W( 62 )
  JVS( 842) = W( 63 )
  JVS( 843) = W( 64 )
  JVS( 844) = W( 65 )
  JVS( 845) = W( 66 )
  JVS( 846) = W( 67 )
  JVS( 847) = W( 68 )
  JVS( 848) = W( 69 )
  JVS( 849) = W( 70 )
  JVS( 850) = W( 71 )
  JVS( 851) = W( 72 )
  JVS( 852) = W( 73 )
  JVS( 853) = W( 74 )
  JVS( 854) = W( 75 )
  JVS( 855) = W( 76 )
  JVS( 856) = W( 77 )
  JVS( 857) = W( 78 )
  JVS( 858) = W( 79 )
  IF ( ABS(  JVS( 899 )) < TINY(a) ) THEN
         IER = 77                                      
         RETURN
  END IF
   W( 13 ) = JVS( 859 )
   W( 20 ) = JVS( 860 )
   W( 21 ) = JVS( 861 )
   W( 24 ) = JVS( 862 )
   W( 26 ) = JVS( 863 )
   W( 27 ) = JVS( 864 )
   W( 33 ) = JVS( 865 )
   W( 34 ) = JVS( 866 )
   W( 35 ) = JVS( 867 )
   W( 36 ) = JVS( 868 )
   W( 37 ) = JVS( 869 )
   W( 38 ) = JVS( 870 )
   W( 39 ) = JVS( 871 )
   W( 42 ) = JVS( 872 )
   W( 43 ) = JVS( 873 )
   W( 48 ) = JVS( 874 )
   W( 50 ) = JVS( 875 )
   W( 51 ) = JVS( 876 )
   W( 52 ) = JVS( 877 )
   W( 54 ) = JVS( 878 )
   W( 55 ) = JVS( 879 )
   W( 56 ) = JVS( 880 )
   W( 57 ) = JVS( 881 )
   W( 58 ) = JVS( 882 )
   W( 59 ) = JVS( 883 )
   W( 62 ) = JVS( 884 )
   W( 63 ) = JVS( 885 )
   W( 64 ) = JVS( 886 )
   W( 65 ) = JVS( 887 )
   W( 66 ) = JVS( 888 )
   W( 67 ) = JVS( 889 )
   W( 68 ) = JVS( 890 )
   W( 69 ) = JVS( 891 )
   W( 70 ) = JVS( 892 )
   W( 71 ) = JVS( 893 )
   W( 72 ) = JVS( 894 )
   W( 73 ) = JVS( 895 )
   W( 74 ) = JVS( 896 )
   W( 75 ) = JVS( 897 )
   W( 76 ) = JVS( 898 )
   W( 77 ) = JVS( 899 )
   W( 78 ) = JVS( 900 )
   W( 79 ) = JVS( 901 )
  a = -W( 13 ) / JVS(           75  )
  W( 13 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 76 )
  a = -W( 20 ) / JVS(           97  )
  W( 20 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 98 )
  a = -W( 21 ) / JVS(           99  )
  W( 21 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 100 )
  a = -W( 24 ) / JVS(          107  )
  W( 24 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 108 )
  a = -W( 26 ) / JVS(          113  )
  W( 26 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 114 )
  a = -W( 27 ) / JVS(          115  )
  W( 27 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 116 )
  a = -W( 33 ) / JVS(          139  )
  W( 33 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 140 )
  a = -W( 34 ) / JVS(          141  )
  W( 34 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 142 )
  a = -W( 35 ) / JVS(          145  )
  W( 35 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 146 )
  a = -W( 36 ) / JVS(          149  )
  W( 36 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 150 )
  a = -W( 37 ) / JVS(          153  )
  W( 37 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 154 )
  W( 75 ) = W( 75 ) + a*JVS( 155 )
  a = -W( 38 ) / JVS(          158  )
  W( 38 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 159 )
  W( 71 ) = W( 71 ) + a*JVS( 160 )
  a = -W( 39 ) / JVS(          161  )
  W( 39 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 162 )
  W( 71 ) = W( 71 ) + a*JVS( 163 )
  a = -W( 42 ) / JVS(          177  )
  W( 42 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 178 )
  W( 71 ) = W( 71 ) + a*JVS( 179 )
  W( 77 ) = W( 77 ) + a*JVS( 180 )
  W( 78 ) = W( 78 ) + a*JVS( 181 )
  a = -W( 43 ) / JVS(          183  )
  W( 43 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 184 )
  W( 71 ) = W( 71 ) + a*JVS( 185 )
  W( 75 ) = W( 75 ) + a*JVS( 186 )
  W( 78 ) = W( 78 ) + a*JVS( 187 )
  a = -W( 48 ) / JVS(          259  )
  W( 48 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 260 )
  W( 67 ) = W( 67 ) + a*JVS( 261 )
  W( 71 ) = W( 71 ) + a*JVS( 262 )
  W( 75 ) = W( 75 ) + a*JVS( 263 )
  a = -W( 50 ) / JVS(          281  )
  W( 50 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 282 )
  W( 67 ) = W( 67 ) + a*JVS( 283 )
  W( 71 ) = W( 71 ) + a*JVS( 284 )
  W( 75 ) = W( 75 ) + a*JVS( 285 )
  a = -W( 51 ) / JVS(          288  )
  W( 51 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 289 )
  W( 70 ) = W( 70 ) + a*JVS( 290 )
  W( 71 ) = W( 71 ) + a*JVS( 291 )
  W( 72 ) = W( 72 ) + a*JVS( 292 )
  W( 73 ) = W( 73 ) + a*JVS( 293 )
  W( 74 ) = W( 74 ) + a*JVS( 294 )
  W( 75 ) = W( 75 ) + a*JVS( 295 )
  W( 78 ) = W( 78 ) + a*JVS( 296 )
  W( 79 ) = W( 79 ) + a*JVS( 297 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 55 ) / JVS(          338  )
  W( 55 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 339 )
  W( 63 ) = W( 63 ) + a*JVS( 340 )
  W( 67 ) = W( 67 ) + a*JVS( 341 )
  W( 71 ) = W( 71 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 62 ) / JVS(          435  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 436 )
  W( 67 ) = W( 67 ) + a*JVS( 437 )
  W( 69 ) = W( 69 ) + a*JVS( 438 )
  W( 71 ) = W( 71 ) + a*JVS( 439 )
  W( 73 ) = W( 73 ) + a*JVS( 440 )
  W( 74 ) = W( 74 ) + a*JVS( 441 )
  W( 75 ) = W( 75 ) + a*JVS( 442 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  a = -W( 64 ) / JVS(          476  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 477 )
  W( 66 ) = W( 66 ) + a*JVS( 478 )
  W( 67 ) = W( 67 ) + a*JVS( 479 )
  W( 69 ) = W( 69 ) + a*JVS( 480 )
  W( 71 ) = W( 71 ) + a*JVS( 481 )
  W( 73 ) = W( 73 ) + a*JVS( 482 )
  W( 74 ) = W( 74 ) + a*JVS( 483 )
  W( 75 ) = W( 75 ) + a*JVS( 484 )
  W( 77 ) = W( 77 ) + a*JVS( 485 )
  W( 78 ) = W( 78 ) + a*JVS( 486 )
  a = -W( 65 ) / JVS(          498  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 499 )
  W( 67 ) = W( 67 ) + a*JVS( 500 )
  W( 69 ) = W( 69 ) + a*JVS( 501 )
  W( 71 ) = W( 71 ) + a*JVS( 502 )
  W( 73 ) = W( 73 ) + a*JVS( 503 )
  W( 74 ) = W( 74 ) + a*JVS( 504 )
  W( 75 ) = W( 75 ) + a*JVS( 505 )
  W( 76 ) = W( 76 ) + a*JVS( 506 )
  W( 77 ) = W( 77 ) + a*JVS( 507 )
  a = -W( 66 ) / JVS(          519  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 520 )
  W( 68 ) = W( 68 ) + a*JVS( 521 )
  W( 69 ) = W( 69 ) + a*JVS( 522 )
  W( 71 ) = W( 71 ) + a*JVS( 523 )
  W( 72 ) = W( 72 ) + a*JVS( 524 )
  W( 73 ) = W( 73 ) + a*JVS( 525 )
  W( 74 ) = W( 74 ) + a*JVS( 526 )
  W( 75 ) = W( 75 ) + a*JVS( 527 )
  W( 76 ) = W( 76 ) + a*JVS( 528 )
  W( 77 ) = W( 77 ) + a*JVS( 529 )
  W( 79 ) = W( 79 ) + a*JVS( 530 )
  a = -W( 67 ) / JVS(          543  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 544 )
  W( 70 ) = W( 70 ) + a*JVS( 545 )
  W( 71 ) = W( 71 ) + a*JVS( 546 )
  W( 72 ) = W( 72 ) + a*JVS( 547 )
  W( 73 ) = W( 73 ) + a*JVS( 548 )
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 78 ) = W( 78 ) + a*JVS( 551 )
  W( 79 ) = W( 79 ) + a*JVS( 552 )
  a = -W( 68 ) / JVS(          574  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 575 )
  W( 70 ) = W( 70 ) + a*JVS( 576 )
  W( 71 ) = W( 71 ) + a*JVS( 577 )
  W( 72 ) = W( 72 ) + a*JVS( 578 )
  W( 73 ) = W( 73 ) + a*JVS( 579 )
  W( 74 ) = W( 74 ) + a*JVS( 580 )
  W( 75 ) = W( 75 ) + a*JVS( 581 )
  W( 76 ) = W( 76 ) + a*JVS( 582 )
  W( 77 ) = W( 77 ) + a*JVS( 583 )
  W( 78 ) = W( 78 ) + a*JVS( 584 )
  W( 79 ) = W( 79 ) + a*JVS( 585 )
  a = -W( 69 ) / JVS(          606  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 607 )
  W( 71 ) = W( 71 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  W( 74 ) = W( 74 ) + a*JVS( 611 )
  W( 75 ) = W( 75 ) + a*JVS( 612 )
  W( 76 ) = W( 76 ) + a*JVS( 613 )
  W( 77 ) = W( 77 ) + a*JVS( 614 )
  W( 78 ) = W( 78 ) + a*JVS( 615 )
  W( 79 ) = W( 79 ) + a*JVS( 616 )
  a = -W( 70 ) / JVS(          627  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 628 )
  W( 72 ) = W( 72 ) + a*JVS( 629 )
  W( 73 ) = W( 73 ) + a*JVS( 630 )
  W( 74 ) = W( 74 ) + a*JVS( 631 )
  W( 75 ) = W( 75 ) + a*JVS( 632 )
  W( 76 ) = W( 76 ) + a*JVS( 633 )
  W( 77 ) = W( 77 ) + a*JVS( 634 )
  W( 78 ) = W( 78 ) + a*JVS( 635 )
  W( 79 ) = W( 79 ) + a*JVS( 636 )
  a = -W( 71 ) / JVS(          687  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 688 )
  W( 73 ) = W( 73 ) + a*JVS( 689 )
  W( 74 ) = W( 74 ) + a*JVS( 690 )
  W( 75 ) = W( 75 ) + a*JVS( 691 )
  W( 76 ) = W( 76 ) + a*JVS( 692 )
  W( 77 ) = W( 77 ) + a*JVS( 693 )
  W( 78 ) = W( 78 ) + a*JVS( 694 )
  W( 79 ) = W( 79 ) + a*JVS( 695 )
  a = -W( 72 ) / JVS(          713  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 714 )
  W( 74 ) = W( 74 ) + a*JVS( 715 )
  W( 75 ) = W( 75 ) + a*JVS( 716 )
  W( 76 ) = W( 76 ) + a*JVS( 717 )
  W( 77 ) = W( 77 ) + a*JVS( 718 )
  W( 78 ) = W( 78 ) + a*JVS( 719 )
  W( 79 ) = W( 79 ) + a*JVS( 720 )
  a = -W( 73 ) / JVS(          739  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 740 )
  W( 75 ) = W( 75 ) + a*JVS( 741 )
  W( 76 ) = W( 76 ) + a*JVS( 742 )
  W( 77 ) = W( 77 ) + a*JVS( 743 )
  W( 78 ) = W( 78 ) + a*JVS( 744 )
  W( 79 ) = W( 79 ) + a*JVS( 745 )
  a = -W( 74 ) / JVS(          782  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 783 )
  W( 76 ) = W( 76 ) + a*JVS( 784 )
  W( 77 ) = W( 77 ) + a*JVS( 785 )
  W( 78 ) = W( 78 ) + a*JVS( 786 )
  W( 79 ) = W( 79 ) + a*JVS( 787 )
  a = -W( 75 ) / JVS(          823  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 824 )
  W( 77 ) = W( 77 ) + a*JVS( 825 )
  W( 78 ) = W( 78 ) + a*JVS( 826 )
  W( 79 ) = W( 79 ) + a*JVS( 827 )
  a = -W( 76 ) / JVS(          855  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 856 )
  W( 78 ) = W( 78 ) + a*JVS( 857 )
  W( 79 ) = W( 79 ) + a*JVS( 858 )
  JVS( 859) = W( 13 )
  JVS( 860) = W( 20 )
  JVS( 861) = W( 21 )
  JVS( 862) = W( 24 )
  JVS( 863) = W( 26 )
  JVS( 864) = W( 27 )
  JVS( 865) = W( 33 )
  JVS( 866) = W( 34 )
  JVS( 867) = W( 35 )
  JVS( 868) = W( 36 )
  JVS( 869) = W( 37 )
  JVS( 870) = W( 38 )
  JVS( 871) = W( 39 )
  JVS( 872) = W( 42 )
  JVS( 873) = W( 43 )
  JVS( 874) = W( 48 )
  JVS( 875) = W( 50 )
  JVS( 876) = W( 51 )
  JVS( 877) = W( 52 )
  JVS( 878) = W( 54 )
  JVS( 879) = W( 55 )
  JVS( 880) = W( 56 )
  JVS( 881) = W( 57 )
  JVS( 882) = W( 58 )
  JVS( 883) = W( 59 )
  JVS( 884) = W( 62 )
  JVS( 885) = W( 63 )
  JVS( 886) = W( 64 )
  JVS( 887) = W( 65 )
  JVS( 888) = W( 66 )
  JVS( 889) = W( 67 )
  JVS( 890) = W( 68 )
  JVS( 891) = W( 69 )
  JVS( 892) = W( 70 )
  JVS( 893) = W( 71 )
  JVS( 894) = W( 72 )
  JVS( 895) = W( 73 )
  JVS( 896) = W( 74 )
  JVS( 897) = W( 75 )
  JVS( 898) = W( 76 )
  JVS( 899) = W( 77 )
  JVS( 900) = W( 78 )
  JVS( 901) = W( 79 )
  IF ( ABS(  JVS( 950 )) < TINY(a) ) THEN
         IER = 78                                      
         RETURN
  END IF
   W( 10 ) = JVS( 902 )
   W( 19 ) = JVS( 903 )
   W( 21 ) = JVS( 904 )
   W( 23 ) = JVS( 905 )
   W( 27 ) = JVS( 906 )
   W( 28 ) = JVS( 907 )
   W( 29 ) = JVS( 908 )
   W( 30 ) = JVS( 909 )
   W( 31 ) = JVS( 910 )
   W( 32 ) = JVS( 911 )
   W( 34 ) = JVS( 912 )
   W( 35 ) = JVS( 913 )
   W( 36 ) = JVS( 914 )
   W( 38 ) = JVS( 915 )
   W( 39 ) = JVS( 916 )
   W( 40 ) = JVS( 917 )
   W( 42 ) = JVS( 918 )
   W( 44 ) = JVS( 919 )
   W( 45 ) = JVS( 920 )
   W( 48 ) = JVS( 921 )
   W( 49 ) = JVS( 922 )
   W( 50 ) = JVS( 923 )
   W( 51 ) = JVS( 924 )
   W( 52 ) = JVS( 925 )
   W( 54 ) = JVS( 926 )
   W( 55 ) = JVS( 927 )
   W( 56 ) = JVS( 928 )
   W( 57 ) = JVS( 929 )
   W( 58 ) = JVS( 930 )
   W( 59 ) = JVS( 931 )
   W( 60 ) = JVS( 932 )
   W( 61 ) = JVS( 933 )
   W( 62 ) = JVS( 934 )
   W( 63 ) = JVS( 935 )
   W( 64 ) = JVS( 936 )
   W( 65 ) = JVS( 937 )
   W( 66 ) = JVS( 938 )
   W( 67 ) = JVS( 939 )
   W( 68 ) = JVS( 940 )
   W( 69 ) = JVS( 941 )
   W( 70 ) = JVS( 942 )
   W( 71 ) = JVS( 943 )
   W( 72 ) = JVS( 944 )
   W( 73 ) = JVS( 945 )
   W( 74 ) = JVS( 946 )
   W( 75 ) = JVS( 947 )
   W( 76 ) = JVS( 948 )
   W( 77 ) = JVS( 949 )
   W( 78 ) = JVS( 950 )
   W( 79 ) = JVS( 951 )
  a = -W( 10 ) / JVS(           69  )
  W( 10 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 70 )
  a = -W( 19 ) / JVS(           94  )
  W( 19 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 95 )
  W( 78 ) = W( 78 ) + a*JVS( 96 )
  a = -W( 21 ) / JVS(           99  )
  W( 21 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 100 )
  a = -W( 23 ) / JVS(          104  )
  W( 23 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 105 )
  W( 73 ) = W( 73 ) + a*JVS( 106 )
  a = -W( 27 ) / JVS(          115  )
  W( 27 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 116 )
  a = -W( 28 ) / JVS(          117  )
  W( 28 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 118 )
  W( 71 ) = W( 71 ) + a*JVS( 119 )
  W( 76 ) = W( 76 ) + a*JVS( 120 )
  W( 77 ) = W( 77 ) + a*JVS( 121 )
  a = -W( 29 ) / JVS(          122  )
  W( 29 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 123 )
  W( 76 ) = W( 76 ) + a*JVS( 124 )
  W( 78 ) = W( 78 ) + a*JVS( 125 )
  a = -W( 30 ) / JVS(          126  )
  W( 30 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 127 )
  W( 73 ) = W( 73 ) + a*JVS( 128 )
  W( 78 ) = W( 78 ) + a*JVS( 129 )
  a = -W( 31 ) / JVS(          130  )
  W( 31 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 131 )
  W( 74 ) = W( 74 ) + a*JVS( 132 )
  W( 75 ) = W( 75 ) + a*JVS( 133 )
  W( 78 ) = W( 78 ) + a*JVS( 134 )
  a = -W( 32 ) / JVS(          135  )
  W( 32 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 136 )
  W( 74 ) = W( 74 ) + a*JVS( 137 )
  W( 78 ) = W( 78 ) + a*JVS( 138 )
  a = -W( 34 ) / JVS(          141  )
  W( 34 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 142 )
  a = -W( 35 ) / JVS(          145  )
  W( 35 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 146 )
  a = -W( 36 ) / JVS(          149  )
  W( 36 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 150 )
  a = -W( 38 ) / JVS(          158  )
  W( 38 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 159 )
  W( 71 ) = W( 71 ) + a*JVS( 160 )
  a = -W( 39 ) / JVS(          161  )
  W( 39 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 162 )
  W( 71 ) = W( 71 ) + a*JVS( 163 )
  a = -W( 40 ) / JVS(          165  )
  W( 40 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 166 )
  W( 74 ) = W( 74 ) + a*JVS( 167 )
  W( 75 ) = W( 75 ) + a*JVS( 168 )
  W( 78 ) = W( 78 ) + a*JVS( 169 )
  a = -W( 42 ) / JVS(          177  )
  W( 42 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 178 )
  W( 71 ) = W( 71 ) + a*JVS( 179 )
  W( 77 ) = W( 77 ) + a*JVS( 180 )
  W( 78 ) = W( 78 ) + a*JVS( 181 )
  a = -W( 44 ) / JVS(          193  )
  W( 44 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  W( 57 ) = W( 57 ) + a*JVS( 195 )
  W( 59 ) = W( 59 ) + a*JVS( 196 )
  W( 67 ) = W( 67 ) + a*JVS( 197 )
  W( 71 ) = W( 71 ) + a*JVS( 198 )
  W( 75 ) = W( 75 ) + a*JVS( 199 )
  a = -W( 45 ) / JVS(          206  )
  W( 45 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 207 )
  W( 49 ) = W( 49 ) + a*JVS( 208 )
  W( 50 ) = W( 50 ) + a*JVS( 209 )
  W( 52 ) = W( 52 ) + a*JVS( 210 )
  W( 54 ) = W( 54 ) + a*JVS( 211 )
  W( 55 ) = W( 55 ) + a*JVS( 212 )
  W( 56 ) = W( 56 ) + a*JVS( 213 )
  W( 57 ) = W( 57 ) + a*JVS( 214 )
  W( 58 ) = W( 58 ) + a*JVS( 215 )
  W( 59 ) = W( 59 ) + a*JVS( 216 )
  W( 60 ) = W( 60 ) + a*JVS( 217 )
  W( 61 ) = W( 61 ) + a*JVS( 218 )
  W( 63 ) = W( 63 ) + a*JVS( 219 )
  W( 64 ) = W( 64 ) + a*JVS( 220 )
  W( 67 ) = W( 67 ) + a*JVS( 221 )
  W( 71 ) = W( 71 ) + a*JVS( 222 )
  W( 75 ) = W( 75 ) + a*JVS( 223 )
  a = -W( 48 ) / JVS(          259  )
  W( 48 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 260 )
  W( 67 ) = W( 67 ) + a*JVS( 261 )
  W( 71 ) = W( 71 ) + a*JVS( 262 )
  W( 75 ) = W( 75 ) + a*JVS( 263 )
  a = -W( 49 ) / JVS(          272  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 273 )
  W( 54 ) = W( 54 ) + a*JVS( 274 )
  W( 57 ) = W( 57 ) + a*JVS( 275 )
  W( 63 ) = W( 63 ) + a*JVS( 276 )
  W( 67 ) = W( 67 ) + a*JVS( 277 )
  W( 71 ) = W( 71 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 78 ) = W( 78 ) + a*JVS( 280 )
  a = -W( 50 ) / JVS(          281  )
  W( 50 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 282 )
  W( 67 ) = W( 67 ) + a*JVS( 283 )
  W( 71 ) = W( 71 ) + a*JVS( 284 )
  W( 75 ) = W( 75 ) + a*JVS( 285 )
  a = -W( 51 ) / JVS(          288  )
  W( 51 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 289 )
  W( 70 ) = W( 70 ) + a*JVS( 290 )
  W( 71 ) = W( 71 ) + a*JVS( 291 )
  W( 72 ) = W( 72 ) + a*JVS( 292 )
  W( 73 ) = W( 73 ) + a*JVS( 293 )
  W( 74 ) = W( 74 ) + a*JVS( 294 )
  W( 75 ) = W( 75 ) + a*JVS( 295 )
  W( 78 ) = W( 78 ) + a*JVS( 296 )
  W( 79 ) = W( 79 ) + a*JVS( 297 )
  a = -W( 52 ) / JVS(          298  )
  W( 52 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 299 )
  W( 67 ) = W( 67 ) + a*JVS( 300 )
  W( 71 ) = W( 71 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  a = -W( 54 ) / JVS(          332  )
  W( 54 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 333 )
  W( 67 ) = W( 67 ) + a*JVS( 334 )
  W( 71 ) = W( 71 ) + a*JVS( 335 )
  W( 75 ) = W( 75 ) + a*JVS( 336 )
  a = -W( 55 ) / JVS(          338  )
  W( 55 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 339 )
  W( 63 ) = W( 63 ) + a*JVS( 340 )
  W( 67 ) = W( 67 ) + a*JVS( 341 )
  W( 71 ) = W( 71 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  a = -W( 56 ) / JVS(          344  )
  W( 56 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 345 )
  W( 67 ) = W( 67 ) + a*JVS( 346 )
  W( 71 ) = W( 71 ) + a*JVS( 347 )
  W( 75 ) = W( 75 ) + a*JVS( 348 )
  a = -W( 57 ) / JVS(          350  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 63 ) = W( 63 ) + a*JVS( 352 )
  W( 67 ) = W( 67 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 59 ) / JVS(          363  )
  W( 59 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 364 )
  W( 67 ) = W( 67 ) + a*JVS( 365 )
  W( 71 ) = W( 71 ) + a*JVS( 366 )
  W( 75 ) = W( 75 ) + a*JVS( 367 )
  a = -W( 60 ) / JVS(          378  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 379 )
  W( 63 ) = W( 63 ) + a*JVS( 380 )
  W( 64 ) = W( 64 ) + a*JVS( 381 )
  W( 65 ) = W( 65 ) + a*JVS( 382 )
  W( 66 ) = W( 66 ) + a*JVS( 383 )
  W( 67 ) = W( 67 ) + a*JVS( 384 )
  W( 68 ) = W( 68 ) + a*JVS( 385 )
  W( 70 ) = W( 70 ) + a*JVS( 386 )
  W( 71 ) = W( 71 ) + a*JVS( 387 )
  W( 72 ) = W( 72 ) + a*JVS( 388 )
  W( 73 ) = W( 73 ) + a*JVS( 389 )
  W( 75 ) = W( 75 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 61 ) / JVS(          412  )
  W( 61 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 413 )
  W( 63 ) = W( 63 ) + a*JVS( 414 )
  W( 65 ) = W( 65 ) + a*JVS( 415 )
  W( 66 ) = W( 66 ) + a*JVS( 416 )
  W( 67 ) = W( 67 ) + a*JVS( 417 )
  W( 68 ) = W( 68 ) + a*JVS( 418 )
  W( 69 ) = W( 69 ) + a*JVS( 419 )
  W( 70 ) = W( 70 ) + a*JVS( 420 )
  W( 71 ) = W( 71 ) + a*JVS( 421 )
  W( 72 ) = W( 72 ) + a*JVS( 422 )
  W( 73 ) = W( 73 ) + a*JVS( 423 )
  W( 74 ) = W( 74 ) + a*JVS( 424 )
  W( 75 ) = W( 75 ) + a*JVS( 425 )
  W( 76 ) = W( 76 ) + a*JVS( 426 )
  W( 77 ) = W( 77 ) + a*JVS( 427 )
  W( 78 ) = W( 78 ) + a*JVS( 428 )
  W( 79 ) = W( 79 ) + a*JVS( 429 )
  a = -W( 62 ) / JVS(          435  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 436 )
  W( 67 ) = W( 67 ) + a*JVS( 437 )
  W( 69 ) = W( 69 ) + a*JVS( 438 )
  W( 71 ) = W( 71 ) + a*JVS( 439 )
  W( 73 ) = W( 73 ) + a*JVS( 440 )
  W( 74 ) = W( 74 ) + a*JVS( 441 )
  W( 75 ) = W( 75 ) + a*JVS( 442 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  a = -W( 64 ) / JVS(          476  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 477 )
  W( 66 ) = W( 66 ) + a*JVS( 478 )
  W( 67 ) = W( 67 ) + a*JVS( 479 )
  W( 69 ) = W( 69 ) + a*JVS( 480 )
  W( 71 ) = W( 71 ) + a*JVS( 481 )
  W( 73 ) = W( 73 ) + a*JVS( 482 )
  W( 74 ) = W( 74 ) + a*JVS( 483 )
  W( 75 ) = W( 75 ) + a*JVS( 484 )
  W( 77 ) = W( 77 ) + a*JVS( 485 )
  W( 78 ) = W( 78 ) + a*JVS( 486 )
  a = -W( 65 ) / JVS(          498  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 499 )
  W( 67 ) = W( 67 ) + a*JVS( 500 )
  W( 69 ) = W( 69 ) + a*JVS( 501 )
  W( 71 ) = W( 71 ) + a*JVS( 502 )
  W( 73 ) = W( 73 ) + a*JVS( 503 )
  W( 74 ) = W( 74 ) + a*JVS( 504 )
  W( 75 ) = W( 75 ) + a*JVS( 505 )
  W( 76 ) = W( 76 ) + a*JVS( 506 )
  W( 77 ) = W( 77 ) + a*JVS( 507 )
  a = -W( 66 ) / JVS(          519  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 520 )
  W( 68 ) = W( 68 ) + a*JVS( 521 )
  W( 69 ) = W( 69 ) + a*JVS( 522 )
  W( 71 ) = W( 71 ) + a*JVS( 523 )
  W( 72 ) = W( 72 ) + a*JVS( 524 )
  W( 73 ) = W( 73 ) + a*JVS( 525 )
  W( 74 ) = W( 74 ) + a*JVS( 526 )
  W( 75 ) = W( 75 ) + a*JVS( 527 )
  W( 76 ) = W( 76 ) + a*JVS( 528 )
  W( 77 ) = W( 77 ) + a*JVS( 529 )
  W( 79 ) = W( 79 ) + a*JVS( 530 )
  a = -W( 67 ) / JVS(          543  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 544 )
  W( 70 ) = W( 70 ) + a*JVS( 545 )
  W( 71 ) = W( 71 ) + a*JVS( 546 )
  W( 72 ) = W( 72 ) + a*JVS( 547 )
  W( 73 ) = W( 73 ) + a*JVS( 548 )
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 78 ) = W( 78 ) + a*JVS( 551 )
  W( 79 ) = W( 79 ) + a*JVS( 552 )
  a = -W( 68 ) / JVS(          574  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 575 )
  W( 70 ) = W( 70 ) + a*JVS( 576 )
  W( 71 ) = W( 71 ) + a*JVS( 577 )
  W( 72 ) = W( 72 ) + a*JVS( 578 )
  W( 73 ) = W( 73 ) + a*JVS( 579 )
  W( 74 ) = W( 74 ) + a*JVS( 580 )
  W( 75 ) = W( 75 ) + a*JVS( 581 )
  W( 76 ) = W( 76 ) + a*JVS( 582 )
  W( 77 ) = W( 77 ) + a*JVS( 583 )
  W( 78 ) = W( 78 ) + a*JVS( 584 )
  W( 79 ) = W( 79 ) + a*JVS( 585 )
  a = -W( 69 ) / JVS(          606  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 607 )
  W( 71 ) = W( 71 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  W( 74 ) = W( 74 ) + a*JVS( 611 )
  W( 75 ) = W( 75 ) + a*JVS( 612 )
  W( 76 ) = W( 76 ) + a*JVS( 613 )
  W( 77 ) = W( 77 ) + a*JVS( 614 )
  W( 78 ) = W( 78 ) + a*JVS( 615 )
  W( 79 ) = W( 79 ) + a*JVS( 616 )
  a = -W( 70 ) / JVS(          627  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 628 )
  W( 72 ) = W( 72 ) + a*JVS( 629 )
  W( 73 ) = W( 73 ) + a*JVS( 630 )
  W( 74 ) = W( 74 ) + a*JVS( 631 )
  W( 75 ) = W( 75 ) + a*JVS( 632 )
  W( 76 ) = W( 76 ) + a*JVS( 633 )
  W( 77 ) = W( 77 ) + a*JVS( 634 )
  W( 78 ) = W( 78 ) + a*JVS( 635 )
  W( 79 ) = W( 79 ) + a*JVS( 636 )
  a = -W( 71 ) / JVS(          687  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 688 )
  W( 73 ) = W( 73 ) + a*JVS( 689 )
  W( 74 ) = W( 74 ) + a*JVS( 690 )
  W( 75 ) = W( 75 ) + a*JVS( 691 )
  W( 76 ) = W( 76 ) + a*JVS( 692 )
  W( 77 ) = W( 77 ) + a*JVS( 693 )
  W( 78 ) = W( 78 ) + a*JVS( 694 )
  W( 79 ) = W( 79 ) + a*JVS( 695 )
  a = -W( 72 ) / JVS(          713  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 714 )
  W( 74 ) = W( 74 ) + a*JVS( 715 )
  W( 75 ) = W( 75 ) + a*JVS( 716 )
  W( 76 ) = W( 76 ) + a*JVS( 717 )
  W( 77 ) = W( 77 ) + a*JVS( 718 )
  W( 78 ) = W( 78 ) + a*JVS( 719 )
  W( 79 ) = W( 79 ) + a*JVS( 720 )
  a = -W( 73 ) / JVS(          739  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 740 )
  W( 75 ) = W( 75 ) + a*JVS( 741 )
  W( 76 ) = W( 76 ) + a*JVS( 742 )
  W( 77 ) = W( 77 ) + a*JVS( 743 )
  W( 78 ) = W( 78 ) + a*JVS( 744 )
  W( 79 ) = W( 79 ) + a*JVS( 745 )
  a = -W( 74 ) / JVS(          782  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 783 )
  W( 76 ) = W( 76 ) + a*JVS( 784 )
  W( 77 ) = W( 77 ) + a*JVS( 785 )
  W( 78 ) = W( 78 ) + a*JVS( 786 )
  W( 79 ) = W( 79 ) + a*JVS( 787 )
  a = -W( 75 ) / JVS(          823  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 824 )
  W( 77 ) = W( 77 ) + a*JVS( 825 )
  W( 78 ) = W( 78 ) + a*JVS( 826 )
  W( 79 ) = W( 79 ) + a*JVS( 827 )
  a = -W( 76 ) / JVS(          855  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 856 )
  W( 78 ) = W( 78 ) + a*JVS( 857 )
  W( 79 ) = W( 79 ) + a*JVS( 858 )
  a = -W( 77 ) / JVS(          899  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 900 )
  W( 79 ) = W( 79 ) + a*JVS( 901 )
  JVS( 902) = W( 10 )
  JVS( 903) = W( 19 )
  JVS( 904) = W( 21 )
  JVS( 905) = W( 23 )
  JVS( 906) = W( 27 )
  JVS( 907) = W( 28 )
  JVS( 908) = W( 29 )
  JVS( 909) = W( 30 )
  JVS( 910) = W( 31 )
  JVS( 911) = W( 32 )
  JVS( 912) = W( 34 )
  JVS( 913) = W( 35 )
  JVS( 914) = W( 36 )
  JVS( 915) = W( 38 )
  JVS( 916) = W( 39 )
  JVS( 917) = W( 40 )
  JVS( 918) = W( 42 )
  JVS( 919) = W( 44 )
  JVS( 920) = W( 45 )
  JVS( 921) = W( 48 )
  JVS( 922) = W( 49 )
  JVS( 923) = W( 50 )
  JVS( 924) = W( 51 )
  JVS( 925) = W( 52 )
  JVS( 926) = W( 54 )
  JVS( 927) = W( 55 )
  JVS( 928) = W( 56 )
  JVS( 929) = W( 57 )
  JVS( 930) = W( 58 )
  JVS( 931) = W( 59 )
  JVS( 932) = W( 60 )
  JVS( 933) = W( 61 )
  JVS( 934) = W( 62 )
  JVS( 935) = W( 63 )
  JVS( 936) = W( 64 )
  JVS( 937) = W( 65 )
  JVS( 938) = W( 66 )
  JVS( 939) = W( 67 )
  JVS( 940) = W( 68 )
  JVS( 941) = W( 69 )
  JVS( 942) = W( 70 )
  JVS( 943) = W( 71 )
  JVS( 944) = W( 72 )
  JVS( 945) = W( 73 )
  JVS( 946) = W( 74 )
  JVS( 947) = W( 75 )
  JVS( 948) = W( 76 )
  JVS( 949) = W( 77 )
  JVS( 950) = W( 78 )
  JVS( 951) = W( 79 )
  IF ( ABS(  JVS( 968 )) < TINY(a) ) THEN
         IER = 79                                      
         RETURN
  END IF
   W( 17 ) = JVS( 952 )
   W( 41 ) = JVS( 953 )
   W( 58 ) = JVS( 954 )
   W( 63 ) = JVS( 955 )
   W( 67 ) = JVS( 956 )
   W( 68 ) = JVS( 957 )
   W( 69 ) = JVS( 958 )
   W( 70 ) = JVS( 959 )
   W( 71 ) = JVS( 960 )
   W( 72 ) = JVS( 961 )
   W( 73 ) = JVS( 962 )
   W( 74 ) = JVS( 963 )
   W( 75 ) = JVS( 964 )
   W( 76 ) = JVS( 965 )
   W( 77 ) = JVS( 966 )
   W( 78 ) = JVS( 967 )
   W( 79 ) = JVS( 968 )
  a = -W( 17 ) / JVS(           88  )
  W( 17 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 89 )
  W( 79 ) = W( 79 ) + a*JVS( 90 )
  a = -W( 41 ) / JVS(          172  )
  W( 41 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 173 )
  W( 67 ) = W( 67 ) + a*JVS( 174 )
  W( 71 ) = W( 71 ) + a*JVS( 175 )
  W( 75 ) = W( 75 ) + a*JVS( 176 )
  a = -W( 58 ) / JVS(          356  )
  W( 58 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 357 )
  W( 67 ) = W( 67 ) + a*JVS( 358 )
  W( 71 ) = W( 71 ) + a*JVS( 359 )
  W( 75 ) = W( 75 ) + a*JVS( 360 )
  a = -W( 63 ) / JVS(          452  )
  W( 63 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 453 )
  W( 71 ) = W( 71 ) + a*JVS( 454 )
  W( 73 ) = W( 73 ) + a*JVS( 455 )
  W( 74 ) = W( 74 ) + a*JVS( 456 )
  W( 75 ) = W( 75 ) + a*JVS( 457 )
  a = -W( 67 ) / JVS(          543  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 544 )
  W( 70 ) = W( 70 ) + a*JVS( 545 )
  W( 71 ) = W( 71 ) + a*JVS( 546 )
  W( 72 ) = W( 72 ) + a*JVS( 547 )
  W( 73 ) = W( 73 ) + a*JVS( 548 )
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 78 ) = W( 78 ) + a*JVS( 551 )
  W( 79 ) = W( 79 ) + a*JVS( 552 )
  a = -W( 68 ) / JVS(          574  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 575 )
  W( 70 ) = W( 70 ) + a*JVS( 576 )
  W( 71 ) = W( 71 ) + a*JVS( 577 )
  W( 72 ) = W( 72 ) + a*JVS( 578 )
  W( 73 ) = W( 73 ) + a*JVS( 579 )
  W( 74 ) = W( 74 ) + a*JVS( 580 )
  W( 75 ) = W( 75 ) + a*JVS( 581 )
  W( 76 ) = W( 76 ) + a*JVS( 582 )
  W( 77 ) = W( 77 ) + a*JVS( 583 )
  W( 78 ) = W( 78 ) + a*JVS( 584 )
  W( 79 ) = W( 79 ) + a*JVS( 585 )
  a = -W( 69 ) / JVS(          606  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 607 )
  W( 71 ) = W( 71 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  W( 74 ) = W( 74 ) + a*JVS( 611 )
  W( 75 ) = W( 75 ) + a*JVS( 612 )
  W( 76 ) = W( 76 ) + a*JVS( 613 )
  W( 77 ) = W( 77 ) + a*JVS( 614 )
  W( 78 ) = W( 78 ) + a*JVS( 615 )
  W( 79 ) = W( 79 ) + a*JVS( 616 )
  a = -W( 70 ) / JVS(          627  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 628 )
  W( 72 ) = W( 72 ) + a*JVS( 629 )
  W( 73 ) = W( 73 ) + a*JVS( 630 )
  W( 74 ) = W( 74 ) + a*JVS( 631 )
  W( 75 ) = W( 75 ) + a*JVS( 632 )
  W( 76 ) = W( 76 ) + a*JVS( 633 )
  W( 77 ) = W( 77 ) + a*JVS( 634 )
  W( 78 ) = W( 78 ) + a*JVS( 635 )
  W( 79 ) = W( 79 ) + a*JVS( 636 )
  a = -W( 71 ) / JVS(          687  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 688 )
  W( 73 ) = W( 73 ) + a*JVS( 689 )
  W( 74 ) = W( 74 ) + a*JVS( 690 )
  W( 75 ) = W( 75 ) + a*JVS( 691 )
  W( 76 ) = W( 76 ) + a*JVS( 692 )
  W( 77 ) = W( 77 ) + a*JVS( 693 )
  W( 78 ) = W( 78 ) + a*JVS( 694 )
  W( 79 ) = W( 79 ) + a*JVS( 695 )
  a = -W( 72 ) / JVS(          713  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 714 )
  W( 74 ) = W( 74 ) + a*JVS( 715 )
  W( 75 ) = W( 75 ) + a*JVS( 716 )
  W( 76 ) = W( 76 ) + a*JVS( 717 )
  W( 77 ) = W( 77 ) + a*JVS( 718 )
  W( 78 ) = W( 78 ) + a*JVS( 719 )
  W( 79 ) = W( 79 ) + a*JVS( 720 )
  a = -W( 73 ) / JVS(          739  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 740 )
  W( 75 ) = W( 75 ) + a*JVS( 741 )
  W( 76 ) = W( 76 ) + a*JVS( 742 )
  W( 77 ) = W( 77 ) + a*JVS( 743 )
  W( 78 ) = W( 78 ) + a*JVS( 744 )
  W( 79 ) = W( 79 ) + a*JVS( 745 )
  a = -W( 74 ) / JVS(          782  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 783 )
  W( 76 ) = W( 76 ) + a*JVS( 784 )
  W( 77 ) = W( 77 ) + a*JVS( 785 )
  W( 78 ) = W( 78 ) + a*JVS( 786 )
  W( 79 ) = W( 79 ) + a*JVS( 787 )
  a = -W( 75 ) / JVS(          823  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 824 )
  W( 77 ) = W( 77 ) + a*JVS( 825 )
  W( 78 ) = W( 78 ) + a*JVS( 826 )
  W( 79 ) = W( 79 ) + a*JVS( 827 )
  a = -W( 76 ) / JVS(          855  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 856 )
  W( 78 ) = W( 78 ) + a*JVS( 857 )
  W( 79 ) = W( 79 ) + a*JVS( 858 )
  a = -W( 77 ) / JVS(          899  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 900 )
  W( 79 ) = W( 79 ) + a*JVS( 901 )
  a = -W( 78 ) / JVS(          950  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 951 )
  JVS( 952) = W( 17 )
  JVS( 953) = W( 41 )
  JVS( 954) = W( 58 )
  JVS( 955) = W( 63 )
  JVS( 956) = W( 67 )
  JVS( 957) = W( 68 )
  JVS( 958) = W( 69 )
  JVS( 959) = W( 70 )
  JVS( 960) = W( 71 )
  JVS( 961) = W( 72 )
  JVS( 962) = W( 73 )
  JVS( 963) = W( 74 )
  JVS( 964) = W( 75 )
  JVS( 965) = W( 76 )
  JVS( 966) = W( 77 )
  JVS( 967) = W( 78 )
  JVS( 968) = W( 79 )
   
   END SUBROUTINE decomp_saprc99
 


END MODULE saprc99_Integrator
