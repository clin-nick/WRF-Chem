
































MODULE racmpm_Integrator

 USE racmpm_Parameters
 USE racmpm_Precision
 USE racmpm_JacobianSP

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

SUBROUTINE  racmpm_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE racmpm_Parameters

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

   CALL racmpm_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  racmpm_INTEGRATE


SUBROUTINE  racmpm_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE racmpm_Parameters

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
      CALL racmpm_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL racmpm_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = racmpm_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL racmpm_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL racmpm_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL racmpm_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL racmpm_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL racmpm_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL racmpm_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL racmpm_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL racmpm_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL racmpm_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL racmpm_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL racmpm_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL racmpm_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL racmpm_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL racmpm_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL racmpm_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  racmpm_ros_ErrorMsg(Code,T,H,IERR)



   USE racmpm_Precision

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

 END SUBROUTINE  racmpm_ros_ErrorMsg


 SUBROUTINE  racmpm_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL racmpm_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL racmpm_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL racmpm_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL racmpm_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL racmpm_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL racmpm_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL racmpm_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL racmpm_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL racmpm_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL racmpm_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL racmpm_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL racmpm_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL racmpm_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL racmpm_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL racmpm_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL racmpm_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL racmpm_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL racmpm_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL racmpm_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL racmpm_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = racmpm_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL racmpm_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  racmpm_ros_Integrator



  REAL(kind=dp) FUNCTION  racmpm_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    racmpm_ros_ErrorNorm = Err

  END FUNCTION  racmpm_ros_ErrorNorm



  SUBROUTINE racmpm_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL racmpm_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL racmpm_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL racmpm_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  racmpm_ros_FunTimeDeriv



  SUBROUTINE  racmpm_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL racmpm_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL racmpm_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL racmpm_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  racmpm_ros_PrepareMatrix



  SUBROUTINE  racmpm_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 



CALL decomp_racmpm ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  racmpm_ros_Decomp



  SUBROUTINE  racmpm_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL racmpm_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  racmpm_ros_Solve




  SUBROUTINE  racmpm_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  racmpm_Ros2



  SUBROUTINE  racmpm_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racmpm_Ros3





  SUBROUTINE  racmpm_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racmpm_Ros4


  SUBROUTINE  racmpm_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racmpm_Rodas3


  SUBROUTINE  racmpm_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racmpm_Rodas4




END SUBROUTINE  racmpm_Rosenbrock




SUBROUTINE  racmpm_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE racmpm_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL racmpm_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  racmpm_FunTemplate



SUBROUTINE  racmpm_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE racmpm_Parameters
 
 USE racmpm_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL racmpm_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  racmpm_JacTemplate

















SUBROUTINE racmpm_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(73)
  A(2) = RCT(2)*V(63)
  A(3) = RCT(3)*V(63)
  A(4) = RCT(4)*V(18)
  A(5) = RCT(5)*V(27)
  A(6) = RCT(6)*V(16)
  A(7) = RCT(7)*V(70)
  A(8) = RCT(8)*V(70)
  A(9) = RCT(9)*V(19)
  A(10) = RCT(10)*V(48)
  A(11) = RCT(11)*V(48)
  A(12) = RCT(12)*V(53)
  A(13) = RCT(13)*V(17)
  A(14) = RCT(14)*V(71)
  A(15) = RCT(15)*V(26)
  A(16) = RCT(16)*V(42)
  A(17) = RCT(17)*V(34)
  A(18) = RCT(18)*V(34)
  A(19) = RCT(19)*V(40)
  A(20) = RCT(20)*V(54)
  A(21) = RCT(21)*V(69)
  A(22) = RCT(22)*V(51)
  A(23) = RCT(23)*V(23)
  A(24) = RCT(24)*V(62)*F(2)
  A(25) = RCT(25)*V(62)*V(63)
  A(26) = RCT(26)*V(6)*F(2)
  A(27) = RCT(27)*V(6)*F(2)
  A(28) = RCT(28)*V(6)*F(1)
  A(29) = RCT(29)*V(63)*V(68)
  A(30) = RCT(30)*V(63)*V(67)
  A(31) = RCT(31)*V(67)*V(68)
  A(32) = RCT(32)*V(19)*V(68)
  A(33) = RCT(33)*V(67)*V(67)
  A(34) = RCT(34)*V(67)*V(67)*F(1)
  A(35) = RCT(35)*V(62)*V(72)
  A(36) = RCT(36)*V(62)*V(73)
  A(37) = RCT(37)*V(62)*V(73)
  A(38) = RCT(38)*V(68)*V(72)
  A(39) = RCT(39)*V(68)*V(73)
  A(40) = RCT(40)*V(68)*V(70)
  A(41) = RCT(41)*V(67)*V(72)
  A(42) = RCT(42)*V(67)*V(73)
  A(43) = RCT(43)*V(16)
  A(44) = RCT(44)*V(67)*V(70)
  A(45) = RCT(45)*V(18)*V(68)
  A(46) = RCT(46)*V(27)*V(68)
  A(47) = RCT(47)*V(16)*V(68)
  A(48) = RCT(48)*V(63)*V(72)
  A(49) = RCT(49)*V(63)*V(73)
  A(50) = RCT(50)*V(72)*V(72)*F(2)
  A(51) = RCT(51)*V(70)*V(72)
  A(52) = RCT(52)*V(70)*V(73)
  A(53) = RCT(53)*V(70)*V(73)
  A(54) = RCT(54)*V(10)
  A(55) = RCT(55)*V(70)*V(70)
  A(56) = RCT(56)*V(68)*F(2)
  A(57) = RCT(57)*V(5)*V(68)
  A(58) = RCT(58)*V(28)*V(68)
  A(59) = RCT(59)*V(37)*V(62)
  A(60) = RCT(60)*V(51)*V(62)
  A(61) = RCT(61)*V(14)*V(68)
  A(62) = RCT(62)*V(13)*V(68)
  A(63) = RCT(63)*V(12)*V(68)
  A(64) = RCT(64)*V(7)*V(68)
  A(65) = RCT(65)*V(11)*V(68)
  A(66) = RCT(66)*V(24)*V(68)
  A(67) = RCT(67)*V(57)*V(68)
  A(68) = RCT(68)*V(58)*V(68)
  A(69) = RCT(69)*V(33)*V(68)
  A(70) = RCT(70)*V(37)*V(68)
  A(71) = RCT(71)*V(29)*V(68)
  A(72) = RCT(72)*V(30)*V(68)
  A(73) = RCT(73)*V(8)*V(68)
  A(74) = RCT(74)*V(9)*V(68)
  A(75) = RCT(75)*V(32)*V(68)
  A(76) = RCT(76)*V(48)*V(68)
  A(77) = RCT(77)*V(53)*V(68)
  A(78) = RCT(78)*V(42)*V(68)
  A(79) = RCT(79)*V(23)*V(68)
  A(80) = RCT(80)*V(34)*V(68)
  A(81) = RCT(81)*V(40)*V(68)
  A(82) = RCT(82)*V(51)*V(68)
  A(83) = RCT(83)*V(54)*V(68)
  A(84) = RCT(84)*V(15)*V(68)
  A(85) = RCT(85)*V(17)*V(68)
  A(86) = RCT(86)*V(68)*V(71)
  A(87) = RCT(87)*V(26)*V(68)
  A(88) = RCT(88)*V(31)*V(68)
  A(89) = RCT(89)*V(35)*V(68)
  A(90) = RCT(90)*V(68)*V(69)
  A(91) = RCT(91)*V(48)*V(70)
  A(92) = RCT(92)*V(53)*V(70)
  A(93) = RCT(93)*V(34)*V(70)
  A(94) = RCT(94)*V(40)*V(70)
  A(95) = RCT(95)*V(51)*V(70)
  A(96) = RCT(96)*V(54)*V(70)
  A(97) = RCT(97)*V(32)*V(70)
  A(98) = RCT(98)*V(24)*V(70)
  A(99) = RCT(99)*V(57)*V(70)
  A(100) = RCT(100)*V(58)*V(70)
  A(101) = RCT(101)*V(33)*V(70)
  A(102) = RCT(102)*V(37)*V(70)
  A(103) = RCT(103)*V(29)*V(70)
  A(104) = RCT(104)*V(30)*V(70)
  A(105) = RCT(105)*V(35)*V(70)
  A(106) = RCT(106)*V(24)*V(63)
  A(107) = RCT(107)*V(57)*V(63)
  A(108) = RCT(108)*V(58)*V(63)
  A(109) = RCT(109)*V(33)*V(63)
  A(110) = RCT(110)*V(37)*V(63)
  A(111) = RCT(111)*V(29)*V(63)
  A(112) = RCT(112)*V(30)*V(63)
  A(113) = RCT(113)*V(51)*V(63)
  A(114) = RCT(114)*V(54)*V(63)
  A(115) = RCT(115)*V(35)*V(63)
  A(116) = RCT(116)*V(20)*V(73)
  A(117) = RCT(117)*V(20)*V(67)
  A(118) = RCT(118)*V(21)*V(73)
  A(119) = RCT(119)*V(21)*F(2)
  A(120) = RCT(120)*V(21)*V(63)
  A(121) = RCT(121)*V(22)*V(73)
  A(122) = RCT(122)*V(22)*F(2)
  A(123) = RCT(123)*V(22)*V(63)
  A(124) = RCT(124)*V(25)*V(73)
  A(125) = RCT(125)*V(25)*F(2)
  A(126) = RCT(126)*V(25)*V(63)
  A(127) = RCT(127)*V(66)*V(73)
  A(128) = RCT(128)*V(31)
  A(129) = RCT(129)*V(55)*V(73)
  A(130) = RCT(130)*V(35)
  A(131) = RCT(131)*V(65)*V(72)
  A(132) = RCT(132)*V(61)*V(72)
  A(133) = RCT(133)*V(52)*V(72)
  A(134) = RCT(134)*V(44)*V(72)
  A(135) = RCT(135)*V(45)*V(72)
  A(136) = RCT(136)*V(36)*V(72)
  A(137) = RCT(137)*V(38)*V(72)
  A(138) = RCT(138)*V(39)*V(72)
  A(139) = RCT(139)*V(50)*V(72)
  A(140) = RCT(140)*V(49)*V(72)
  A(141) = RCT(141)*V(43)*V(72)
  A(142) = RCT(142)*V(46)*V(72)
  A(143) = RCT(143)*V(47)*V(72)
  A(144) = RCT(144)*V(41)*V(72)
  A(145) = RCT(145)*V(66)*V(72)
  A(146) = RCT(146)*V(55)*V(72)
  A(147) = RCT(147)*V(64)*V(72)
  A(148) = RCT(148)*V(59)*V(72)
  A(149) = RCT(149)*V(60)*V(72)
  A(150) = RCT(150)*V(65)*V(67)
  A(151) = RCT(151)*V(61)*V(67)
  A(152) = RCT(152)*V(52)*V(67)
  A(153) = RCT(153)*V(44)*V(67)
  A(154) = RCT(154)*V(45)*V(67)
  A(155) = RCT(155)*V(36)*V(67)
  A(156) = RCT(156)*V(38)*V(67)
  A(157) = RCT(157)*V(39)*V(67)
  A(158) = RCT(158)*V(50)*V(67)
  A(159) = RCT(159)*V(49)*V(67)
  A(160) = RCT(160)*V(43)*V(67)
  A(161) = RCT(161)*V(46)*V(67)
  A(162) = RCT(162)*V(47)*V(67)
  A(163) = RCT(163)*V(41)*V(67)
  A(164) = RCT(164)*V(66)*V(67)
  A(165) = RCT(165)*V(66)*V(67)
  A(166) = RCT(166)*V(55)*V(67)
  A(167) = RCT(167)*V(55)*V(67)
  A(168) = RCT(168)*V(64)*V(67)
  A(169) = RCT(169)*V(59)*V(67)
  A(170) = RCT(170)*V(60)*V(67)
  A(171) = RCT(171)*V(65)*V(65)
  A(172) = RCT(172)*V(61)*V(65)
  A(173) = RCT(173)*V(52)*V(65)
  A(174) = RCT(174)*V(44)*V(65)
  A(175) = RCT(175)*V(45)*V(65)
  A(176) = RCT(176)*V(36)*V(65)
  A(177) = RCT(177)*V(38)*V(65)
  A(178) = RCT(178)*V(39)*V(65)
  A(179) = RCT(179)*V(50)*V(65)
  A(180) = RCT(180)*V(49)*V(65)
  A(181) = RCT(181)*V(43)*V(65)
  A(182) = RCT(182)*V(46)*V(65)
  A(183) = RCT(183)*V(47)*V(65)
  A(184) = RCT(184)*V(41)*V(65)
  A(185) = RCT(185)*V(65)*V(66)
  A(186) = RCT(186)*V(65)*V(66)
  A(187) = RCT(187)*V(55)*V(65)
  A(188) = RCT(188)*V(55)*V(65)
  A(189) = RCT(189)*V(64)*V(65)
  A(190) = RCT(190)*V(59)*V(65)
  A(191) = RCT(191)*V(60)*V(65)
  A(192) = RCT(192)*V(61)*V(66)
  A(193) = RCT(193)*V(52)*V(66)
  A(194) = RCT(194)*V(44)*V(66)
  A(195) = RCT(195)*V(45)*V(66)
  A(196) = RCT(196)*V(36)*V(66)
  A(197) = RCT(197)*V(38)*V(66)
  A(198) = RCT(198)*V(39)*V(66)
  A(199) = RCT(199)*V(50)*V(66)
  A(200) = RCT(200)*V(49)*V(66)
  A(201) = RCT(201)*V(43)*V(66)
  A(202) = RCT(202)*V(46)*V(66)
  A(203) = RCT(203)*V(47)*V(66)
  A(204) = RCT(204)*V(41)*V(66)
  A(205) = RCT(205)*V(66)*V(66)
  A(206) = RCT(206)*V(55)*V(66)
  A(207) = RCT(207)*V(64)*V(66)
  A(208) = RCT(208)*V(59)*V(66)
  A(209) = RCT(209)*V(60)*V(66)
  A(210) = RCT(210)*V(59)*V(59)
  A(211) = RCT(211)*V(59)*V(60)
  A(212) = RCT(212)*V(60)*V(60)
  A(213) = RCT(213)*V(65)*V(70)
  A(214) = RCT(214)*V(61)*V(70)
  A(215) = RCT(215)*V(52)*V(70)
  A(216) = RCT(216)*V(44)*V(70)
  A(217) = RCT(217)*V(45)*V(70)
  A(218) = RCT(218)*V(36)*V(70)
  A(219) = RCT(219)*V(38)*V(70)
  A(220) = RCT(220)*V(39)*V(70)
  A(221) = RCT(221)*V(50)*V(70)
  A(222) = RCT(222)*V(49)*V(70)
  A(223) = RCT(223)*V(43)*V(70)
  A(224) = RCT(224)*V(46)*V(70)
  A(225) = RCT(225)*V(47)*V(70)
  A(226) = RCT(226)*V(41)*V(70)
  A(227) = RCT(227)*V(66)*V(70)
  A(228) = RCT(228)*V(55)*V(70)
  A(229) = RCT(229)*V(64)*V(70)
  A(230) = RCT(230)*V(59)*V(70)
  A(231) = RCT(231)*V(60)*V(70)
  A(232) = RCT(232)*V(56)*V(67)
  A(233) = RCT(233)*V(56)*V(65)
  A(234) = RCT(234)*V(56)*V(66)
  A(235) = RCT(235)*V(56)*V(56)
  A(236) = RCT(236)*V(56)*V(72)
  A(237) = RCT(237)*V(56)*V(70)


  Vdot(1) = A(57)
  Vdot(2) = A(58)
  Vdot(3) = 0.036*A(63)+0.37*A(106)+0.14*A(107)+0.15*A(109)+0.15*A(110)+0.01*A(112)+0.22*A(113)+0.11*A(114)+0.11*A(115)
  Vdot(4) = 0.1*A(107)+0.14*A(108)+0.07*A(112)+0.13*A(113)+0.21*A(114)+A(165)+A(167)+A(186)+A(188)+0.5*A(192)+0.499&
              &*A(193)+0.495*A(194)+0.495*A(195)+0.5*A(196)+0.499*A(197)+0.49*A(198)+0.494*A(199)+0.5*A(207)+0.5*A(208)&
              &+0.484*A(209)
  Vdot(5) = -A(57)
  Vdot(6) = A(2)-A(26)-A(27)-A(28)
  Vdot(7) = -A(64)
  Vdot(8) = -A(73)
  Vdot(9) = -A(74)
  Vdot(10) = A(53)-A(54)
  Vdot(11) = -A(65)
  Vdot(12) = -A(63)
  Vdot(13) = -A(62)+0.03*A(107)+0.06*A(108)
  Vdot(14) = -A(61)+0.06*A(107)+0.07*A(108)
  Vdot(15) = 0.35*A(83)-A(84)
  Vdot(16) = -A(6)+A(42)-A(43)-A(47)
  Vdot(17) = -A(13)-A(85)+A(150)
  Vdot(18) = -A(4)+A(38)-A(45)+A(118)+A(121)+A(124)
  Vdot(19) = -A(9)-A(32)+A(33)+A(34)+0.006*A(107)+0.011*A(108)+0.001*A(109)+0.001*A(110)+0.02*A(111)+0.02*A(112)
  Vdot(20) = 0.1*A(75)+A(97)-A(116)-A(117)
  Vdot(21) = 0.9*A(73)-A(118)-A(119)-A(120)
  Vdot(22) = 0.9*A(74)-A(121)-A(122)-A(123)
  Vdot(23) = -A(23)+0.024*A(65)-A(79)+0.41*A(82)+0.6*A(89)+0.3*A(189)
  Vdot(24) = -A(66)-A(98)-A(106)
  Vdot(25) = 0.85*A(75)-A(124)-A(125)-A(126)
  Vdot(26) = -A(15)-A(87)+0.11*A(114)+A(164)
  Vdot(27) = -A(5)+A(39)+0.3*A(44)-A(46)+A(91)+A(92)+A(93)+A(94)+0.2*A(95)+0.5*A(96)+A(97)
  Vdot(28) = A(10)+A(11)+A(12)+1.87*A(17)+1.55*A(18)+A(19)+A(22)-A(58)+0.01*A(59)+0.036*A(63)+A(76)+2*A(80)+A(81)+0.41&
               &*A(82)+A(91)+2*A(93)+A(94)+0.8*A(95)+0.43*A(106)+0.37*A(107)+0.3*A(108)+0.36*A(109)+0.36*A(110)+0.14*A(111)&
               &+0.14*A(112)+0.54*A(113)+0.66*A(114)+0.13*A(115)
  Vdot(29) = -A(71)-A(103)-A(111)
  Vdot(30) = -A(72)-A(104)-A(112)
  Vdot(31) = -A(88)+0.4*A(89)+0.4*A(105)+0.3*A(115)+A(127)-A(128)
  Vdot(32) = -A(75)-A(97)+0.1*A(116)+A(117)+A(118)+0.02*A(119)+A(120)+A(121)+0.02*A(122)+A(123)+A(124)+0.02*A(125)&
               &+A(126)
  Vdot(33) = -A(69)-A(101)-A(109)
  Vdot(34) = -A(17)-A(18)+0.036*A(63)-A(80)+0.15*A(83)-A(93)+0.25*A(96)+0.5*A(114)+0.063*A(133)+1.2*A(142)+0.35*A(143)&
               &+A(144)+0.119*A(173)+0.65*A(182)+0.37*A(183)+A(184)+0.1*A(193)+0.65*A(202)+0.37*A(203)+A(204)+0.063*A(215)&
               &+1.3*A(224)+0.74*A(225)+A(226)
  Vdot(35) = -A(89)-A(105)-A(115)+A(129)-A(130)
  Vdot(36) = A(66)-A(136)-A(155)-A(176)-A(196)-A(218)
  Vdot(37) = -A(59)-A(70)-A(102)-A(110)
  Vdot(38) = A(67)-A(137)-A(156)-A(177)-A(197)-A(219)
  Vdot(39) = A(68)-A(138)-A(157)-A(178)-A(198)-A(220)
  Vdot(40) = -A(19)+A(79)-A(81)+0.08*A(82)+0.15*A(83)-A(94)+0.25*A(96)+0.6*A(113)+0.62*A(114)+0.65*A(142)+0.6*A(143)&
               &+A(144)+0.54*A(147)+0.005*A(173)+0.35*A(182)+0.63*A(183)+A(184)+0.4*A(189)+0.004*A(193)+0.35*A(202)+0.63&
               &*A(203)+A(204)+0.54*A(207)+0.7*A(224)+1.26*A(225)+A(226)+0.54*A(229)
  Vdot(41) = 0.98*A(125)-A(144)-A(163)-A(184)-A(204)-A(226)
  Vdot(42) = -A(16)+0.8*A(21)+0.25*A(64)-A(78)+0.12*A(84)+0.41*A(86)+0.03*A(96)+0.03*A(107)+0.16*A(108)+0.53*A(111)&
               &+0.623*A(133)+0.722*A(134)+0.642*A(135)+0.06*A(137)+0.29*A(138)+0.8*A(140)+0.464*A(149)+0.018*A(173)+0.24&
               &*A(174)+0.419*A(175)+0.081*A(177)+0.313*A(178)+A(180)+0.149*A(191)+0.127*A(193)+0.33*A(194)+0.581*A(195)&
               &+0.141*A(197)+0.569*A(198)+A(200)+0.11*A(207)+0.167*A(209)+0.149*A(211)+0.285*A(212)+0.67*A(215)+0.828&
               &*A(216)+0.88*A(217)+0.06*A(219)+0.29*A(220)+A(222)+0.469*A(231)
  Vdot(43) = A(72)-A(141)-A(160)-A(181)-A(201)-A(223)
  Vdot(44) = 0.75*A(64)-A(134)-A(153)-A(174)-A(194)-A(216)
  Vdot(45) = 0.9511*A(65)-A(135)-A(154)-A(175)-A(195)-A(217)
  Vdot(46) = 0.98*A(119)-A(142)-A(161)-A(182)-A(202)-A(224)
  Vdot(47) = 0.98*A(122)-A(143)-A(162)-A(183)-A(203)-A(225)
  Vdot(48) = -A(10)-A(11)+A(13)+0.13*A(17)+0.45*A(18)+A(22)+A(23)+0.05*A(59)+0.01*A(63)-A(76)+0.08*A(82)+0.35*A(85)+0.35&
               &*A(87)+A(88)+0.4*A(89)-A(91)+0.4*A(105)+A(106)+0.64*A(107)+0.02*A(108)+0.9*A(109)+0.9*A(110)+0.04*A(112)+0.4&
               &*A(113)+0.7*A(115)+A(131)+0.047*A(133)+0.021*A(134)+1.6*A(136)+A(137)+0.606*A(139)+0.25*A(141)+A(146)+0.287&
               &*A(149)+1.33*A(171)+0.75*A(172)+0.81*A(173)+0.829*A(174)+0.753*A(175)+1.55*A(176)+1.25*A(177)+0.755*A(178)&
               &+1.09*A(179)+A(180)+1.4*A(181)+A(182)+A(183)+A(184)+A(185)+A(186)+2*A(187)+A(188)+0.75*A(189)+0.75*A(190)&
               &+0.96*A(191)+0.091*A(193)+0.076*A(194)+0.8*A(196)+0.501*A(197)+0.34*A(199)+0.4*A(201)+A(206)+0.207*A(209)&
               &+0.202*A(211)+0.504*A(212)+A(213)+0.048*A(215)+0.021*A(216)+1.6*A(218)+A(219)+0.686*A(221)+0.4*A(223)+A(228)&
               &+0.28*A(231)+A(233)
  Vdot(49) = A(71)-A(140)-A(159)-A(180)-A(200)-A(222)
  Vdot(50) = A(69)+A(70)-A(139)-A(158)-A(179)-A(199)-A(221)
  Vdot(51) = -A(22)-A(60)-A(82)-A(95)+0.9*A(101)+0.9*A(102)+0.39*A(109)+0.39*A(110)+0.79*A(112)-A(113)+0.446*A(139)+0.4&
               &*A(141)+0.55*A(179)+0.6*A(181)+0.771*A(199)+0.6*A(201)+0.6*A(221)+0.6*A(223)
  Vdot(52) = 0.583*A(63)+0.44*A(86)+A(90)-A(133)-A(152)-A(173)-A(193)-A(215)
  Vdot(53) = -A(12)+A(14)+0.2*A(21)+A(60)+0.335*A(63)+0.025*A(65)-A(77)+0.88*A(84)+0.08*A(86)-A(92)+0.25*A(96)+0.44&
               &*A(107)+0.99*A(108)+0.65*A(111)+0.16*A(114)+A(132)+0.233*A(133)+0.211*A(134)+0.15*A(135)+0.2*A(136)+0.94&
               &*A(137)+1.71*A(138)+0.8*A(140)+0.46*A(147)+1.24*A(149)+0.75*A(172)+0.58*A(173)+0.523*A(174)+0.411*A(175)&
               &+0.35*A(176)+0.669*A(177)+0.932*A(178)+A(180)+0.3*A(189)+0.64*A(191)+A(192)+0.724*A(193)+0.677*A(194)+0.497&
               &*A(195)+0.6*A(196)+0.859*A(197)+0.941*A(198)+A(200)+0.35*A(207)+0.65*A(209)+0.64*A(211)+1.21*A(212)+A(214)&
               &+0.243*A(215)+0.239*A(216)+0.187*A(217)+0.2*A(218)+0.94*A(219)+1.71*A(220)+A(222)+0.46*A(229)+1.24*A(231)
  Vdot(54) = -A(20)+0.13*A(59)-A(83)-A(96)-A(114)+0.5*A(142)+0.95*A(143)+A(182)+A(183)+A(202)+A(203)+0.5*A(224)+A(225)
  Vdot(55) = A(20)+0.51*A(82)+0.5*A(83)+0.2*A(95)+0.5*A(96)-A(129)+A(130)-A(146)-A(166)-A(167)-A(187)-A(188)-A(206)&
               &-A(228)
  Vdot(56) = 0.15*A(59)+0.1*A(73)+0.1*A(74)+0.05*A(75)+0.49*A(82)+0.5*A(83)+0.07*A(86)+0.35*A(87)+A(88)+A(89)+0.5*A(96)&
               &+A(105)+0.13*A(109)+0.13*A(110)+0.048*A(133)+0.334*A(134)+0.416*A(135)+0.16*A(147)+0.085*A(173)+0.245*A(174)&
               &+0.322*A(175)+0.08*A(189)+0.071*A(193)+0.237*A(194)+0.318*A(195)+0.08*A(207)+0.051*A(215)+0.391*A(216)+0.587&
               &*A(217)+0.16*A(229)-A(232)-A(233)-A(234)-2*A(235)-A(236)-A(237)
  Vdot(57) = 0.86*A(59)-A(67)-A(99)-A(107)+0.35*A(109)+0.35*A(110)+0.46*A(112)+0.354*A(139)+0.37*A(179)+0.229*A(199)+0.4&
               &*A(221)
  Vdot(58) = -A(68)-A(100)-A(108)+0.25*A(141)+0.08*A(179)+0.4*A(181)+0.4*A(201)+0.4*A(223)
  Vdot(59) = 0.8*A(95)+0.8*A(98)+0.43*A(99)+0.11*A(100)+0.9*A(101)+0.9*A(102)+0.1*A(103)+0.13*A(104)-A(148)-A(169)&
               &-A(190)-A(208)-2*A(210)-A(211)-A(230)
  Vdot(60) = 0.2*A(98)+0.57*A(99)+0.89*A(100)+0.1*A(101)+0.1*A(102)+0.9*A(103)+0.87*A(104)-A(149)-A(170)-A(191)-A(209)&
               &-A(211)-2*A(212)-A(231)
  Vdot(61) = A(16)+A(62)+0.1*A(107)+0.18*A(108)+0.2*A(111)+0.16*A(112)-A(132)+0.048*A(133)+0.245*A(134)+0.133*A(135)&
               &-A(151)-A(172)+0.014*A(174)+0.013*A(175)-A(192)+0.006*A(193)+0.018*A(194)+0.015*A(195)-A(214)+0.053*A(215)&
               &+0.262*A(216)+0.155*A(217)
  Vdot(62) = A(1)+A(3)+A(8)-A(24)-A(25)+A(27)-A(35)-A(36)-A(37)-A(59)-A(60)+0.09*A(109)+0.09*A(110)
  Vdot(63) = -A(2)-A(3)+A(24)-A(25)+A(26)-A(29)-A(30)-A(48)-A(49)-A(106)-A(107)-A(108)-A(109)-A(110)-A(111)-A(112)&
               &-A(113)-A(114)-A(115)-A(120)-A(123)-A(126)+A(165)+A(167)
  Vdot(64) = A(78)+0.03*A(107)+0.12*A(108)+0.02*A(109)+0.02*A(110)+0.42*A(111)+0.42*A(112)-A(147)-A(168)-A(189)-A(207)&
               &-A(229)
  Vdot(65) = A(12)+A(15)+A(61)+0.65*A(85)+0.19*A(107)+0.23*A(108)+0.03*A(109)+0.03*A(110)-A(131)+0.15*A(133)+0.031&
               &*A(134)+A(145)-A(150)-2*A(171)-A(172)-0.993*A(173)-0.951*A(174)-A(175)-A(176)-A(177)-A(178)-A(179)-A(180)&
               &-A(181)-A(182)-A(183)-A(184)-A(186)-A(187)-A(188)-A(189)-A(190)-A(191)+0.5*A(192)+0.508*A(193)+0.554*A(194)&
               &+0.507*A(195)+0.5*A(196)+0.501*A(197)+0.51*A(198)+0.506*A(199)+A(200)+A(201)+A(202)+A(203)+A(204)+2*A(205)&
               &+A(206)+0.5*A(207)+0.5*A(208)+0.516*A(209)-A(213)+0.155*A(215)+0.04*A(216)+A(227)-A(233)+A(234)
  Vdot(66) = A(16)+A(19)+A(22)+A(23)+A(77)+A(81)+0.65*A(87)+A(92)+A(94)+0.15*A(109)+0.15*A(110)+0.13*A(113)+0.28*A(114)&
               &+0.7*A(115)-A(127)+A(128)-A(145)+A(146)+0.23*A(147)-A(164)-A(165)-A(185)-A(186)+A(187)+0.12*A(189)-A(192)&
               &-A(193)-A(194)-A(195)-A(196)-A(197)-A(198)-A(199)-A(200)-A(201)-A(202)-A(203)-A(204)-2*A(205)-0.88*A(207)&
               &-A(208)-A(209)-A(227)+A(228)+0.23*A(229)-A(234)
  Vdot(67) = 0.65*A(6)+2*A(11)+A(12)+A(13)+A(14)+0.8*A(18)+A(19)+A(20)+A(21)+A(22)+A(23)+A(29)-A(30)-A(31)+A(32)-2*A(33)&
               &-2*A(34)+A(40)-A(41)-A(42)+A(43)-A(44)+A(56)+A(57)+A(58)+0.28*A(59)+0.381*A(63)+0.25*A(64)+0.049*A(65)+0.1&
               &*A(73)+0.1*A(74)+0.05*A(75)+A(76)+A(79)+A(80)+0.49*A(82)+0.5*A(83)+A(84)+0.35*A(87)+0.4*A(89)+A(91)+A(93)&
               &+0.5*A(96)+0.26*A(106)+0.25*A(107)+0.22*A(108)+0.3*A(109)+0.3*A(110)+0.1*A(111)+0.1*A(112)+0.29*A(113)+0.29&
               &*A(114)+0.08*A(115)-A(117)+0.02*A(119)+0.02*A(122)+0.02*A(125)+A(131)+A(132)+0.742*A(133)+0.599*A(134)+0.606&
               &*A(135)+A(136)+A(137)+A(138)+0.847*A(139)+0.8*A(140)+0.65*A(141)+0.95*A(142)+0.95*A(143)+A(144)+0.77*A(147)&
               &+A(148)-A(150)-A(151)-A(152)-A(153)-A(154)-A(155)-A(156)-A(157)-A(158)-A(159)-A(160)-A(161)-A(162)-A(163)&
               &-A(164)-A(165)-A(166)-A(167)-A(168)-A(169)-A(170)+0.66*A(171)+A(172)+0.992*A(173)+0.946*A(174)+0.993*A(175)&
               &+A(176)+A(177)+A(178)+A(179)+2*A(180)+2*A(181)+A(182)+A(183)+2*A(184)+A(185)+A(187)+0.88*A(189)+A(190)+0.5&
               &*A(191)+0.5*A(192)+0.488*A(193)+0.438*A(194)+0.489*A(195)+0.5*A(196)+0.501*A(197)+0.51*A(198)+0.506*A(199)&
               &+A(200)+A(201)+A(202)+A(203)+A(204)+0.38*A(207)+0.5*A(208)+A(210)+0.5*A(211)+A(213)+A(214)+0.792*A(215)&
               &+0.699*A(216)+0.845*A(217)+A(218)+A(219)+A(220)+A(221)+A(222)+A(223)+A(224)+A(225)+A(226)+0.77*A(229)+A(230)&
               &-A(232)+A(233)
  Vdot(68) = A(4)+A(5)+0.35*A(6)+2*A(9)+A(13)+A(14)+A(15)+2*A(28)-A(29)+A(30)-A(31)-A(32)-A(38)-A(39)-A(40)+A(41)+0.7&
               &*A(44)-A(45)-A(46)-A(47)-A(56)-A(57)-A(58)+0.02*A(59)-A(61)-A(62)-0.964*A(63)-A(64)-A(65)-A(66)-A(67)-A(68)&
               &-A(69)-A(70)-A(71)-A(72)-A(73)-A(74)-A(75)-A(76)-A(77)-A(78)-A(79)-A(80)-A(81)-A(82)-A(83)-A(84)-0.65*A(85)&
               &-0.51*A(86)-A(87)-A(88)-A(89)-A(90)+0.12*A(106)+0.4*A(107)+0.63*A(108)+0.28*A(109)+0.28*A(110)+0.85*A(111)&
               &+0.85*A(112)+0.07*A(113)+0.21*A(114)+0.036*A(115)+A(120)+A(123)+A(126)
  Vdot(69) = -A(21)-A(90)+0.6*A(105)+A(116)+0.059*A(133)+0.124*A(134)+0.261*A(135)+0.153*A(139)+0.2*A(140)+0.35*A(141)&
               &+0.05*A(142)+0.05*A(143)+A(148)+A(169)+A(170)+A(190)+0.5*A(191)+A(208)+0.484*A(209)+2*A(210)+1.5*A(211)&
               &+A(212)+A(230)
  Vdot(70) = 0.35*A(6)-A(7)-A(8)+A(37)-A(40)-A(44)+A(46)+A(49)-A(51)-A(52)-A(53)+A(54)-2*A(55)+A(88)+0.6*A(89)-A(91)&
               &-A(92)-A(93)-A(94)-A(95)-A(96)-A(97)-A(98)-A(99)-A(100)-A(101)-A(102)-A(103)-A(104)-0.4*A(105)-A(213)-A(214)&
               &-A(215)-A(216)-A(217)-A(218)-A(219)-A(220)-A(221)-A(222)-A(223)-A(224)-A(225)-A(226)-A(227)-A(228)-A(229)&
               &-A(230)-A(231)-A(237)
  Vdot(71) = -A(14)-A(86)+0.13*A(113)+A(151)+A(152)+A(153)+A(154)+A(155)+A(156)+A(157)+A(158)+A(159)+A(160)+A(161)&
               &+A(162)+A(163)+A(166)+A(168)+A(232)
  Vdot(72) = A(1)+A(4)+A(7)-A(35)+A(36)-A(38)-A(41)-A(48)-2*A(50)-A(51)+A(52)-A(131)-A(132)-A(133)-A(134)-A(135)-A(136)&
               &-A(137)-A(138)-A(139)-A(140)-A(141)-A(142)-A(143)-A(144)-A(145)-A(146)-A(147)-A(148)-A(149)-A(236)
  Vdot(73) = -A(1)+A(5)+0.65*A(6)+A(8)+A(21)+A(35)-A(36)-A(37)-A(39)+A(40)+A(41)-A(42)+A(43)+0.7*A(44)+A(45)+A(47)+A(48)&
               &-A(49)+2*A(50)+2*A(51)-A(53)+A(54)+2*A(55)+A(90)+0.5*A(96)+0.4*A(105)+0.7*A(115)-A(116)-A(118)-A(121)-A(124)&
               &-A(127)+A(128)-A(129)+A(130)+A(131)+A(132)+0.941*A(133)+0.876*A(134)+0.739*A(135)+A(136)+A(137)+A(138)+0.847&
               &*A(139)+0.8*A(140)+0.65*A(141)+0.95*A(142)+0.95*A(143)+A(144)+A(145)+A(146)+A(147)+A(148)+2*A(149)+0.5&
               &*A(191)+0.516*A(209)+0.5*A(211)+A(212)+A(213)+A(214)+A(215)+A(216)+A(217)+A(218)+A(219)+A(220)+A(221)+A(222)&
               &+A(223)+A(224)+A(225)+A(226)+A(227)+A(228)+A(229)+A(230)+2*A(231)+A(236)+A(237)
      
END SUBROUTINE racmpm_Fun
















SUBROUTINE racmpm_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(73)
  IRR(2) = RCT(2)*V(63)
  IRR(3) = RCT(3)*V(63)
  IRR(4) = RCT(4)*V(18)
  IRR(5) = RCT(5)*V(27)
  IRR(6) = RCT(6)*V(16)
  IRR(7) = RCT(7)*V(70)
  IRR(8) = RCT(8)*V(70)
  IRR(9) = RCT(9)*V(19)
  IRR(10) = RCT(10)*V(48)
  IRR(11) = RCT(11)*V(48)
  IRR(12) = RCT(12)*V(53)
  IRR(13) = RCT(13)*V(17)
  IRR(14) = RCT(14)*V(71)
  IRR(15) = RCT(15)*V(26)
  IRR(16) = RCT(16)*V(42)
  IRR(17) = RCT(17)*V(34)
  IRR(18) = RCT(18)*V(34)
  IRR(19) = RCT(19)*V(40)
  IRR(20) = RCT(20)*V(54)
  IRR(21) = RCT(21)*V(69)
  IRR(22) = RCT(22)*V(51)
  IRR(23) = RCT(23)*V(23)
  IRR(24) = RCT(24)*V(62)*F(2)
  IRR(25) = RCT(25)*V(62)*V(63)
  IRR(26) = RCT(26)*V(6)*F(2)
  IRR(27) = RCT(27)*V(6)*F(2)
  IRR(28) = RCT(28)*V(6)*F(1)
  IRR(29) = RCT(29)*V(63)*V(68)
  IRR(30) = RCT(30)*V(63)*V(67)
  IRR(31) = RCT(31)*V(67)*V(68)
  IRR(32) = RCT(32)*V(19)*V(68)
  IRR(33) = RCT(33)*V(67)*V(67)
  IRR(34) = RCT(34)*V(67)*V(67)*F(1)
  IRR(35) = RCT(35)*V(62)*V(72)
  IRR(36) = RCT(36)*V(62)*V(73)
  IRR(37) = RCT(37)*V(62)*V(73)
  IRR(38) = RCT(38)*V(68)*V(72)
  IRR(39) = RCT(39)*V(68)*V(73)
  IRR(40) = RCT(40)*V(68)*V(70)
  IRR(41) = RCT(41)*V(67)*V(72)
  IRR(42) = RCT(42)*V(67)*V(73)
  IRR(43) = RCT(43)*V(16)
  IRR(44) = RCT(44)*V(67)*V(70)
  IRR(45) = RCT(45)*V(18)*V(68)
  IRR(46) = RCT(46)*V(27)*V(68)
  IRR(47) = RCT(47)*V(16)*V(68)
  IRR(48) = RCT(48)*V(63)*V(72)
  IRR(49) = RCT(49)*V(63)*V(73)
  IRR(50) = RCT(50)*V(72)*V(72)*F(2)
  IRR(51) = RCT(51)*V(70)*V(72)
  IRR(52) = RCT(52)*V(70)*V(73)
  IRR(53) = RCT(53)*V(70)*V(73)
  IRR(54) = RCT(54)*V(10)
  IRR(55) = RCT(55)*V(70)*V(70)
  IRR(56) = RCT(56)*V(68)*F(2)
  IRR(57) = RCT(57)*V(5)*V(68)
  IRR(58) = RCT(58)*V(28)*V(68)
  IRR(59) = RCT(59)*V(37)*V(62)
  IRR(60) = RCT(60)*V(51)*V(62)
  IRR(61) = RCT(61)*V(14)*V(68)
  IRR(62) = RCT(62)*V(13)*V(68)
  IRR(63) = RCT(63)*V(12)*V(68)
  IRR(64) = RCT(64)*V(7)*V(68)
  IRR(65) = RCT(65)*V(11)*V(68)
  IRR(66) = RCT(66)*V(24)*V(68)
  IRR(67) = RCT(67)*V(57)*V(68)
  IRR(68) = RCT(68)*V(58)*V(68)
  IRR(69) = RCT(69)*V(33)*V(68)
  IRR(70) = RCT(70)*V(37)*V(68)
  IRR(71) = RCT(71)*V(29)*V(68)
  IRR(72) = RCT(72)*V(30)*V(68)
  IRR(73) = RCT(73)*V(8)*V(68)
  IRR(74) = RCT(74)*V(9)*V(68)
  IRR(75) = RCT(75)*V(32)*V(68)
  IRR(76) = RCT(76)*V(48)*V(68)
  IRR(77) = RCT(77)*V(53)*V(68)
  IRR(78) = RCT(78)*V(42)*V(68)
  IRR(79) = RCT(79)*V(23)*V(68)
  IRR(80) = RCT(80)*V(34)*V(68)
  IRR(81) = RCT(81)*V(40)*V(68)
  IRR(82) = RCT(82)*V(51)*V(68)
  IRR(83) = RCT(83)*V(54)*V(68)
  IRR(84) = RCT(84)*V(15)*V(68)
  IRR(85) = RCT(85)*V(17)*V(68)
  IRR(86) = RCT(86)*V(68)*V(71)
  IRR(87) = RCT(87)*V(26)*V(68)
  IRR(88) = RCT(88)*V(31)*V(68)
  IRR(89) = RCT(89)*V(35)*V(68)
  IRR(90) = RCT(90)*V(68)*V(69)
  IRR(91) = RCT(91)*V(48)*V(70)
  IRR(92) = RCT(92)*V(53)*V(70)
  IRR(93) = RCT(93)*V(34)*V(70)
  IRR(94) = RCT(94)*V(40)*V(70)
  IRR(95) = RCT(95)*V(51)*V(70)
  IRR(96) = RCT(96)*V(54)*V(70)
  IRR(97) = RCT(97)*V(32)*V(70)
  IRR(98) = RCT(98)*V(24)*V(70)
  IRR(99) = RCT(99)*V(57)*V(70)
  IRR(100) = RCT(100)*V(58)*V(70)
  IRR(101) = RCT(101)*V(33)*V(70)
  IRR(102) = RCT(102)*V(37)*V(70)
  IRR(103) = RCT(103)*V(29)*V(70)
  IRR(104) = RCT(104)*V(30)*V(70)
  IRR(105) = RCT(105)*V(35)*V(70)
  IRR(106) = RCT(106)*V(24)*V(63)
  IRR(107) = RCT(107)*V(57)*V(63)
  IRR(108) = RCT(108)*V(58)*V(63)
  IRR(109) = RCT(109)*V(33)*V(63)
  IRR(110) = RCT(110)*V(37)*V(63)
  IRR(111) = RCT(111)*V(29)*V(63)
  IRR(112) = RCT(112)*V(30)*V(63)
  IRR(113) = RCT(113)*V(51)*V(63)
  IRR(114) = RCT(114)*V(54)*V(63)
  IRR(115) = RCT(115)*V(35)*V(63)
  IRR(116) = RCT(116)*V(20)*V(73)
  IRR(117) = RCT(117)*V(20)*V(67)
  IRR(118) = RCT(118)*V(21)*V(73)
  IRR(119) = RCT(119)*V(21)*F(2)
  IRR(120) = RCT(120)*V(21)*V(63)
  IRR(121) = RCT(121)*V(22)*V(73)
  IRR(122) = RCT(122)*V(22)*F(2)
  IRR(123) = RCT(123)*V(22)*V(63)
  IRR(124) = RCT(124)*V(25)*V(73)
  IRR(125) = RCT(125)*V(25)*F(2)
  IRR(126) = RCT(126)*V(25)*V(63)
  IRR(127) = RCT(127)*V(66)*V(73)
  IRR(128) = RCT(128)*V(31)
  IRR(129) = RCT(129)*V(55)*V(73)
  IRR(130) = RCT(130)*V(35)
  IRR(131) = RCT(131)*V(65)*V(72)
  IRR(132) = RCT(132)*V(61)*V(72)
  IRR(133) = RCT(133)*V(52)*V(72)
  IRR(134) = RCT(134)*V(44)*V(72)
  IRR(135) = RCT(135)*V(45)*V(72)
  IRR(136) = RCT(136)*V(36)*V(72)
  IRR(137) = RCT(137)*V(38)*V(72)
  IRR(138) = RCT(138)*V(39)*V(72)
  IRR(139) = RCT(139)*V(50)*V(72)
  IRR(140) = RCT(140)*V(49)*V(72)
  IRR(141) = RCT(141)*V(43)*V(72)
  IRR(142) = RCT(142)*V(46)*V(72)
  IRR(143) = RCT(143)*V(47)*V(72)
  IRR(144) = RCT(144)*V(41)*V(72)
  IRR(145) = RCT(145)*V(66)*V(72)
  IRR(146) = RCT(146)*V(55)*V(72)
  IRR(147) = RCT(147)*V(64)*V(72)
  IRR(148) = RCT(148)*V(59)*V(72)
  IRR(149) = RCT(149)*V(60)*V(72)
  IRR(150) = RCT(150)*V(65)*V(67)
  IRR(151) = RCT(151)*V(61)*V(67)
  IRR(152) = RCT(152)*V(52)*V(67)
  IRR(153) = RCT(153)*V(44)*V(67)
  IRR(154) = RCT(154)*V(45)*V(67)
  IRR(155) = RCT(155)*V(36)*V(67)
  IRR(156) = RCT(156)*V(38)*V(67)
  IRR(157) = RCT(157)*V(39)*V(67)
  IRR(158) = RCT(158)*V(50)*V(67)
  IRR(159) = RCT(159)*V(49)*V(67)
  IRR(160) = RCT(160)*V(43)*V(67)
  IRR(161) = RCT(161)*V(46)*V(67)
  IRR(162) = RCT(162)*V(47)*V(67)
  IRR(163) = RCT(163)*V(41)*V(67)
  IRR(164) = RCT(164)*V(66)*V(67)
  IRR(165) = RCT(165)*V(66)*V(67)
  IRR(166) = RCT(166)*V(55)*V(67)
  IRR(167) = RCT(167)*V(55)*V(67)
  IRR(168) = RCT(168)*V(64)*V(67)
  IRR(169) = RCT(169)*V(59)*V(67)
  IRR(170) = RCT(170)*V(60)*V(67)
  IRR(171) = RCT(171)*V(65)*V(65)
  IRR(172) = RCT(172)*V(61)*V(65)
  IRR(173) = RCT(173)*V(52)*V(65)
  IRR(174) = RCT(174)*V(44)*V(65)
  IRR(175) = RCT(175)*V(45)*V(65)
  IRR(176) = RCT(176)*V(36)*V(65)
  IRR(177) = RCT(177)*V(38)*V(65)
  IRR(178) = RCT(178)*V(39)*V(65)
  IRR(179) = RCT(179)*V(50)*V(65)
  IRR(180) = RCT(180)*V(49)*V(65)
  IRR(181) = RCT(181)*V(43)*V(65)
  IRR(182) = RCT(182)*V(46)*V(65)
  IRR(183) = RCT(183)*V(47)*V(65)
  IRR(184) = RCT(184)*V(41)*V(65)
  IRR(185) = RCT(185)*V(65)*V(66)
  IRR(186) = RCT(186)*V(65)*V(66)
  IRR(187) = RCT(187)*V(55)*V(65)
  IRR(188) = RCT(188)*V(55)*V(65)
  IRR(189) = RCT(189)*V(64)*V(65)
  IRR(190) = RCT(190)*V(59)*V(65)
  IRR(191) = RCT(191)*V(60)*V(65)
  IRR(192) = RCT(192)*V(61)*V(66)
  IRR(193) = RCT(193)*V(52)*V(66)
  IRR(194) = RCT(194)*V(44)*V(66)
  IRR(195) = RCT(195)*V(45)*V(66)
  IRR(196) = RCT(196)*V(36)*V(66)
  IRR(197) = RCT(197)*V(38)*V(66)
  IRR(198) = RCT(198)*V(39)*V(66)
  IRR(199) = RCT(199)*V(50)*V(66)
  IRR(200) = RCT(200)*V(49)*V(66)
  IRR(201) = RCT(201)*V(43)*V(66)
  IRR(202) = RCT(202)*V(46)*V(66)
  IRR(203) = RCT(203)*V(47)*V(66)
  IRR(204) = RCT(204)*V(41)*V(66)
  IRR(205) = RCT(205)*V(66)*V(66)
  IRR(206) = RCT(206)*V(55)*V(66)
  IRR(207) = RCT(207)*V(64)*V(66)
  IRR(208) = RCT(208)*V(59)*V(66)
  IRR(209) = RCT(209)*V(60)*V(66)
  IRR(210) = RCT(210)*V(59)*V(59)
  IRR(211) = RCT(211)*V(59)*V(60)
  IRR(212) = RCT(212)*V(60)*V(60)
  IRR(213) = RCT(213)*V(65)*V(70)
  IRR(214) = RCT(214)*V(61)*V(70)
  IRR(215) = RCT(215)*V(52)*V(70)
  IRR(216) = RCT(216)*V(44)*V(70)
  IRR(217) = RCT(217)*V(45)*V(70)
  IRR(218) = RCT(218)*V(36)*V(70)
  IRR(219) = RCT(219)*V(38)*V(70)
  IRR(220) = RCT(220)*V(39)*V(70)
  IRR(221) = RCT(221)*V(50)*V(70)
  IRR(222) = RCT(222)*V(49)*V(70)
  IRR(223) = RCT(223)*V(43)*V(70)
  IRR(224) = RCT(224)*V(46)*V(70)
  IRR(225) = RCT(225)*V(47)*V(70)
  IRR(226) = RCT(226)*V(41)*V(70)
  IRR(227) = RCT(227)*V(66)*V(70)
  IRR(228) = RCT(228)*V(55)*V(70)
  IRR(229) = RCT(229)*V(64)*V(70)
  IRR(230) = RCT(230)*V(59)*V(70)
  IRR(231) = RCT(231)*V(60)*V(70)
  IRR(232) = RCT(232)*V(56)*V(67)
  IRR(233) = RCT(233)*V(56)*V(65)
  IRR(234) = RCT(234)*V(56)*V(66)
  IRR(235) = RCT(235)*V(56)*V(56)
  IRR(236) = RCT(236)*V(56)*V(72)
  IRR(237) = RCT(237)*V(56)*V(70)
      
END SUBROUTINE racmpm_IRRFun
















SUBROUTINE racmpm_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(440)


  B(1) = RCT(1)

  B(2) = RCT(2)

  B(3) = RCT(3)

  B(4) = RCT(4)

  B(5) = RCT(5)

  B(6) = RCT(6)

  B(7) = RCT(7)

  B(8) = RCT(8)

  B(9) = RCT(9)

  B(10) = RCT(10)

  B(11) = RCT(11)

  B(12) = RCT(12)

  B(13) = RCT(13)

  B(14) = RCT(14)

  B(15) = RCT(15)

  B(16) = RCT(16)

  B(17) = RCT(17)

  B(18) = RCT(18)

  B(19) = RCT(19)

  B(20) = RCT(20)

  B(21) = RCT(21)

  B(22) = RCT(22)

  B(23) = RCT(23)

  B(24) = RCT(24)*F(2)

  B(26) = RCT(25)*V(63)

  B(27) = RCT(25)*V(62)

  B(28) = RCT(26)*F(2)

  B(30) = RCT(27)*F(2)

  B(32) = RCT(28)*F(1)

  B(34) = RCT(29)*V(68)

  B(35) = RCT(29)*V(63)

  B(36) = RCT(30)*V(67)

  B(37) = RCT(30)*V(63)

  B(38) = RCT(31)*V(68)

  B(39) = RCT(31)*V(67)

  B(40) = RCT(32)*V(68)

  B(41) = RCT(32)*V(19)

  B(42) = RCT(33)*2*V(67)

  B(43) = RCT(34)*2*V(67)*F(1)

  B(45) = RCT(35)*V(72)

  B(46) = RCT(35)*V(62)

  B(47) = RCT(36)*V(73)

  B(48) = RCT(36)*V(62)

  B(49) = RCT(37)*V(73)

  B(50) = RCT(37)*V(62)

  B(51) = RCT(38)*V(72)

  B(52) = RCT(38)*V(68)

  B(53) = RCT(39)*V(73)

  B(54) = RCT(39)*V(68)

  B(55) = RCT(40)*V(70)

  B(56) = RCT(40)*V(68)

  B(57) = RCT(41)*V(72)

  B(58) = RCT(41)*V(67)

  B(59) = RCT(42)*V(73)

  B(60) = RCT(42)*V(67)

  B(61) = RCT(43)

  B(62) = RCT(44)*V(70)

  B(63) = RCT(44)*V(67)

  B(64) = RCT(45)*V(68)

  B(65) = RCT(45)*V(18)

  B(66) = RCT(46)*V(68)

  B(67) = RCT(46)*V(27)

  B(68) = RCT(47)*V(68)

  B(69) = RCT(47)*V(16)

  B(70) = RCT(48)*V(72)

  B(71) = RCT(48)*V(63)

  B(72) = RCT(49)*V(73)

  B(73) = RCT(49)*V(63)

  B(74) = RCT(50)*2*V(72)*F(2)

  B(76) = RCT(51)*V(72)

  B(77) = RCT(51)*V(70)

  B(78) = RCT(52)*V(73)

  B(79) = RCT(52)*V(70)

  B(80) = RCT(53)*V(73)

  B(81) = RCT(53)*V(70)

  B(82) = RCT(54)

  B(83) = RCT(55)*2*V(70)

  B(84) = RCT(56)*F(2)

  B(86) = RCT(57)*V(68)

  B(87) = RCT(57)*V(5)

  B(88) = RCT(58)*V(68)

  B(89) = RCT(58)*V(28)

  B(90) = RCT(59)*V(62)

  B(91) = RCT(59)*V(37)

  B(92) = RCT(60)*V(62)

  B(93) = RCT(60)*V(51)

  B(94) = RCT(61)*V(68)

  B(95) = RCT(61)*V(14)

  B(96) = RCT(62)*V(68)

  B(97) = RCT(62)*V(13)

  B(98) = RCT(63)*V(68)

  B(99) = RCT(63)*V(12)

  B(100) = RCT(64)*V(68)

  B(101) = RCT(64)*V(7)

  B(102) = RCT(65)*V(68)

  B(103) = RCT(65)*V(11)

  B(104) = RCT(66)*V(68)

  B(105) = RCT(66)*V(24)

  B(106) = RCT(67)*V(68)

  B(107) = RCT(67)*V(57)

  B(108) = RCT(68)*V(68)

  B(109) = RCT(68)*V(58)

  B(110) = RCT(69)*V(68)

  B(111) = RCT(69)*V(33)

  B(112) = RCT(70)*V(68)

  B(113) = RCT(70)*V(37)

  B(114) = RCT(71)*V(68)

  B(115) = RCT(71)*V(29)

  B(116) = RCT(72)*V(68)

  B(117) = RCT(72)*V(30)

  B(118) = RCT(73)*V(68)

  B(119) = RCT(73)*V(8)

  B(120) = RCT(74)*V(68)

  B(121) = RCT(74)*V(9)

  B(122) = RCT(75)*V(68)

  B(123) = RCT(75)*V(32)

  B(124) = RCT(76)*V(68)

  B(125) = RCT(76)*V(48)

  B(126) = RCT(77)*V(68)

  B(127) = RCT(77)*V(53)

  B(128) = RCT(78)*V(68)

  B(129) = RCT(78)*V(42)

  B(130) = RCT(79)*V(68)

  B(131) = RCT(79)*V(23)

  B(132) = RCT(80)*V(68)

  B(133) = RCT(80)*V(34)

  B(134) = RCT(81)*V(68)

  B(135) = RCT(81)*V(40)

  B(136) = RCT(82)*V(68)

  B(137) = RCT(82)*V(51)

  B(138) = RCT(83)*V(68)

  B(139) = RCT(83)*V(54)

  B(140) = RCT(84)*V(68)

  B(141) = RCT(84)*V(15)

  B(142) = RCT(85)*V(68)

  B(143) = RCT(85)*V(17)

  B(144) = RCT(86)*V(71)

  B(145) = RCT(86)*V(68)

  B(146) = RCT(87)*V(68)

  B(147) = RCT(87)*V(26)

  B(148) = RCT(88)*V(68)

  B(149) = RCT(88)*V(31)

  B(150) = RCT(89)*V(68)

  B(151) = RCT(89)*V(35)

  B(152) = RCT(90)*V(69)

  B(153) = RCT(90)*V(68)

  B(154) = RCT(91)*V(70)

  B(155) = RCT(91)*V(48)

  B(156) = RCT(92)*V(70)

  B(157) = RCT(92)*V(53)

  B(158) = RCT(93)*V(70)

  B(159) = RCT(93)*V(34)

  B(160) = RCT(94)*V(70)

  B(161) = RCT(94)*V(40)

  B(162) = RCT(95)*V(70)

  B(163) = RCT(95)*V(51)

  B(164) = RCT(96)*V(70)

  B(165) = RCT(96)*V(54)

  B(166) = RCT(97)*V(70)

  B(167) = RCT(97)*V(32)

  B(168) = RCT(98)*V(70)

  B(169) = RCT(98)*V(24)

  B(170) = RCT(99)*V(70)

  B(171) = RCT(99)*V(57)

  B(172) = RCT(100)*V(70)

  B(173) = RCT(100)*V(58)

  B(174) = RCT(101)*V(70)

  B(175) = RCT(101)*V(33)

  B(176) = RCT(102)*V(70)

  B(177) = RCT(102)*V(37)

  B(178) = RCT(103)*V(70)

  B(179) = RCT(103)*V(29)

  B(180) = RCT(104)*V(70)

  B(181) = RCT(104)*V(30)

  B(182) = RCT(105)*V(70)

  B(183) = RCT(105)*V(35)

  B(184) = RCT(106)*V(63)

  B(185) = RCT(106)*V(24)

  B(186) = RCT(107)*V(63)

  B(187) = RCT(107)*V(57)

  B(188) = RCT(108)*V(63)

  B(189) = RCT(108)*V(58)

  B(190) = RCT(109)*V(63)

  B(191) = RCT(109)*V(33)

  B(192) = RCT(110)*V(63)

  B(193) = RCT(110)*V(37)

  B(194) = RCT(111)*V(63)

  B(195) = RCT(111)*V(29)

  B(196) = RCT(112)*V(63)

  B(197) = RCT(112)*V(30)

  B(198) = RCT(113)*V(63)

  B(199) = RCT(113)*V(51)

  B(200) = RCT(114)*V(63)

  B(201) = RCT(114)*V(54)

  B(202) = RCT(115)*V(63)

  B(203) = RCT(115)*V(35)

  B(204) = RCT(116)*V(73)

  B(205) = RCT(116)*V(20)

  B(206) = RCT(117)*V(67)

  B(207) = RCT(117)*V(20)

  B(208) = RCT(118)*V(73)

  B(209) = RCT(118)*V(21)

  B(210) = RCT(119)*F(2)

  B(212) = RCT(120)*V(63)

  B(213) = RCT(120)*V(21)

  B(214) = RCT(121)*V(73)

  B(215) = RCT(121)*V(22)

  B(216) = RCT(122)*F(2)

  B(218) = RCT(123)*V(63)

  B(219) = RCT(123)*V(22)

  B(220) = RCT(124)*V(73)

  B(221) = RCT(124)*V(25)

  B(222) = RCT(125)*F(2)

  B(224) = RCT(126)*V(63)

  B(225) = RCT(126)*V(25)

  B(226) = RCT(127)*V(73)

  B(227) = RCT(127)*V(66)

  B(228) = RCT(128)

  B(229) = RCT(129)*V(73)

  B(230) = RCT(129)*V(55)

  B(231) = RCT(130)

  B(232) = RCT(131)*V(72)

  B(233) = RCT(131)*V(65)

  B(234) = RCT(132)*V(72)

  B(235) = RCT(132)*V(61)

  B(236) = RCT(133)*V(72)

  B(237) = RCT(133)*V(52)

  B(238) = RCT(134)*V(72)

  B(239) = RCT(134)*V(44)

  B(240) = RCT(135)*V(72)

  B(241) = RCT(135)*V(45)

  B(242) = RCT(136)*V(72)

  B(243) = RCT(136)*V(36)

  B(244) = RCT(137)*V(72)

  B(245) = RCT(137)*V(38)

  B(246) = RCT(138)*V(72)

  B(247) = RCT(138)*V(39)

  B(248) = RCT(139)*V(72)

  B(249) = RCT(139)*V(50)

  B(250) = RCT(140)*V(72)

  B(251) = RCT(140)*V(49)

  B(252) = RCT(141)*V(72)

  B(253) = RCT(141)*V(43)

  B(254) = RCT(142)*V(72)

  B(255) = RCT(142)*V(46)

  B(256) = RCT(143)*V(72)

  B(257) = RCT(143)*V(47)

  B(258) = RCT(144)*V(72)

  B(259) = RCT(144)*V(41)

  B(260) = RCT(145)*V(72)

  B(261) = RCT(145)*V(66)

  B(262) = RCT(146)*V(72)

  B(263) = RCT(146)*V(55)

  B(264) = RCT(147)*V(72)

  B(265) = RCT(147)*V(64)

  B(266) = RCT(148)*V(72)

  B(267) = RCT(148)*V(59)

  B(268) = RCT(149)*V(72)

  B(269) = RCT(149)*V(60)

  B(270) = RCT(150)*V(67)

  B(271) = RCT(150)*V(65)

  B(272) = RCT(151)*V(67)

  B(273) = RCT(151)*V(61)

  B(274) = RCT(152)*V(67)

  B(275) = RCT(152)*V(52)

  B(276) = RCT(153)*V(67)

  B(277) = RCT(153)*V(44)

  B(278) = RCT(154)*V(67)

  B(279) = RCT(154)*V(45)

  B(280) = RCT(155)*V(67)

  B(281) = RCT(155)*V(36)

  B(282) = RCT(156)*V(67)

  B(283) = RCT(156)*V(38)

  B(284) = RCT(157)*V(67)

  B(285) = RCT(157)*V(39)

  B(286) = RCT(158)*V(67)

  B(287) = RCT(158)*V(50)

  B(288) = RCT(159)*V(67)

  B(289) = RCT(159)*V(49)

  B(290) = RCT(160)*V(67)

  B(291) = RCT(160)*V(43)

  B(292) = RCT(161)*V(67)

  B(293) = RCT(161)*V(46)

  B(294) = RCT(162)*V(67)

  B(295) = RCT(162)*V(47)

  B(296) = RCT(163)*V(67)

  B(297) = RCT(163)*V(41)

  B(298) = RCT(164)*V(67)

  B(299) = RCT(164)*V(66)

  B(300) = RCT(165)*V(67)

  B(301) = RCT(165)*V(66)

  B(302) = RCT(166)*V(67)

  B(303) = RCT(166)*V(55)

  B(304) = RCT(167)*V(67)

  B(305) = RCT(167)*V(55)

  B(306) = RCT(168)*V(67)

  B(307) = RCT(168)*V(64)

  B(308) = RCT(169)*V(67)

  B(309) = RCT(169)*V(59)

  B(310) = RCT(170)*V(67)

  B(311) = RCT(170)*V(60)

  B(312) = RCT(171)*2*V(65)

  B(313) = RCT(172)*V(65)

  B(314) = RCT(172)*V(61)

  B(315) = RCT(173)*V(65)

  B(316) = RCT(173)*V(52)

  B(317) = RCT(174)*V(65)

  B(318) = RCT(174)*V(44)

  B(319) = RCT(175)*V(65)

  B(320) = RCT(175)*V(45)

  B(321) = RCT(176)*V(65)

  B(322) = RCT(176)*V(36)

  B(323) = RCT(177)*V(65)

  B(324) = RCT(177)*V(38)

  B(325) = RCT(178)*V(65)

  B(326) = RCT(178)*V(39)

  B(327) = RCT(179)*V(65)

  B(328) = RCT(179)*V(50)

  B(329) = RCT(180)*V(65)

  B(330) = RCT(180)*V(49)

  B(331) = RCT(181)*V(65)

  B(332) = RCT(181)*V(43)

  B(333) = RCT(182)*V(65)

  B(334) = RCT(182)*V(46)

  B(335) = RCT(183)*V(65)

  B(336) = RCT(183)*V(47)

  B(337) = RCT(184)*V(65)

  B(338) = RCT(184)*V(41)

  B(339) = RCT(185)*V(66)

  B(340) = RCT(185)*V(65)

  B(341) = RCT(186)*V(66)

  B(342) = RCT(186)*V(65)

  B(343) = RCT(187)*V(65)

  B(344) = RCT(187)*V(55)

  B(345) = RCT(188)*V(65)

  B(346) = RCT(188)*V(55)

  B(347) = RCT(189)*V(65)

  B(348) = RCT(189)*V(64)

  B(349) = RCT(190)*V(65)

  B(350) = RCT(190)*V(59)

  B(351) = RCT(191)*V(65)

  B(352) = RCT(191)*V(60)

  B(353) = RCT(192)*V(66)

  B(354) = RCT(192)*V(61)

  B(355) = RCT(193)*V(66)

  B(356) = RCT(193)*V(52)

  B(357) = RCT(194)*V(66)

  B(358) = RCT(194)*V(44)

  B(359) = RCT(195)*V(66)

  B(360) = RCT(195)*V(45)

  B(361) = RCT(196)*V(66)

  B(362) = RCT(196)*V(36)

  B(363) = RCT(197)*V(66)

  B(364) = RCT(197)*V(38)

  B(365) = RCT(198)*V(66)

  B(366) = RCT(198)*V(39)

  B(367) = RCT(199)*V(66)

  B(368) = RCT(199)*V(50)

  B(369) = RCT(200)*V(66)

  B(370) = RCT(200)*V(49)

  B(371) = RCT(201)*V(66)

  B(372) = RCT(201)*V(43)

  B(373) = RCT(202)*V(66)

  B(374) = RCT(202)*V(46)

  B(375) = RCT(203)*V(66)

  B(376) = RCT(203)*V(47)

  B(377) = RCT(204)*V(66)

  B(378) = RCT(204)*V(41)

  B(379) = RCT(205)*2*V(66)

  B(380) = RCT(206)*V(66)

  B(381) = RCT(206)*V(55)

  B(382) = RCT(207)*V(66)

  B(383) = RCT(207)*V(64)

  B(384) = RCT(208)*V(66)

  B(385) = RCT(208)*V(59)

  B(386) = RCT(209)*V(66)

  B(387) = RCT(209)*V(60)

  B(388) = RCT(210)*2*V(59)

  B(389) = RCT(211)*V(60)

  B(390) = RCT(211)*V(59)

  B(391) = RCT(212)*2*V(60)

  B(392) = RCT(213)*V(70)

  B(393) = RCT(213)*V(65)

  B(394) = RCT(214)*V(70)

  B(395) = RCT(214)*V(61)

  B(396) = RCT(215)*V(70)

  B(397) = RCT(215)*V(52)

  B(398) = RCT(216)*V(70)

  B(399) = RCT(216)*V(44)

  B(400) = RCT(217)*V(70)

  B(401) = RCT(217)*V(45)

  B(402) = RCT(218)*V(70)

  B(403) = RCT(218)*V(36)

  B(404) = RCT(219)*V(70)

  B(405) = RCT(219)*V(38)

  B(406) = RCT(220)*V(70)

  B(407) = RCT(220)*V(39)

  B(408) = RCT(221)*V(70)

  B(409) = RCT(221)*V(50)

  B(410) = RCT(222)*V(70)

  B(411) = RCT(222)*V(49)

  B(412) = RCT(223)*V(70)

  B(413) = RCT(223)*V(43)

  B(414) = RCT(224)*V(70)

  B(415) = RCT(224)*V(46)

  B(416) = RCT(225)*V(70)

  B(417) = RCT(225)*V(47)

  B(418) = RCT(226)*V(70)

  B(419) = RCT(226)*V(41)

  B(420) = RCT(227)*V(70)

  B(421) = RCT(227)*V(66)

  B(422) = RCT(228)*V(70)

  B(423) = RCT(228)*V(55)

  B(424) = RCT(229)*V(70)

  B(425) = RCT(229)*V(64)

  B(426) = RCT(230)*V(70)

  B(427) = RCT(230)*V(59)

  B(428) = RCT(231)*V(70)

  B(429) = RCT(231)*V(60)

  B(430) = RCT(232)*V(67)

  B(431) = RCT(232)*V(56)

  B(432) = RCT(233)*V(65)

  B(433) = RCT(233)*V(56)

  B(434) = RCT(234)*V(66)

  B(435) = RCT(234)*V(56)

  B(436) = RCT(235)*2*V(56)

  B(437) = RCT(236)*V(72)

  B(438) = RCT(236)*V(56)

  B(439) = RCT(237)*V(70)

  B(440) = RCT(237)*V(56)



  JVS(1) = 0

  JVS(2) = B(86)

  JVS(3) = B(87)

  JVS(4) = 0

  JVS(5) = B(88)

  JVS(6) = B(89)

  JVS(7) = 0

  JVS(8) = 0.036*B(98)

  JVS(9) = 0.37*B(184)

  JVS(10) = 0.01*B(196)

  JVS(11) = 0.15*B(190)

  JVS(12) = 0.11*B(202)

  JVS(13) = 0.15*B(192)

  JVS(14) = 0.22*B(198)

  JVS(15) = 0.11*B(200)

  JVS(16) = 0.14*B(186)

  JVS(17) = 0.37*B(185)+0.14*B(187)+0.15*B(191)+0.15*B(193)+0.01*B(197)+0.22*B(199)+0.11*B(201)+0.11*B(203)

  JVS(18) = 0.036*B(99)

  JVS(19) = 0

  JVS(20) = 0.07*B(196)

  JVS(21) = 0.5*B(361)

  JVS(22) = 0.499*B(363)

  JVS(23) = 0.49*B(365)

  JVS(24) = 0.495*B(357)

  JVS(25) = 0.495*B(359)

  JVS(26) = 0.494*B(367)

  JVS(27) = 0.13*B(198)

  JVS(28) = 0.499*B(355)

  JVS(29) = 0.21*B(200)

  JVS(30) = B(304)+B(345)

  JVS(31) = 0.1*B(186)

  JVS(32) = 0.14*B(188)

  JVS(33) = 0.5*B(384)

  JVS(34) = 0.484*B(386)

  JVS(35) = 0.5*B(353)

  JVS(36) = 0.1*B(187)+0.14*B(189)+0.07*B(197)+0.13*B(199)+0.21*B(201)

  JVS(37) = 0.5*B(382)

  JVS(38) = B(341)+B(346)

  JVS(39) = B(300)+B(342)+0.5*B(354)+0.499*B(356)+0.495*B(358)+0.495*B(360)+0.5*B(362)+0.499*B(364)+0.49*B(366)+0.494&
              &*B(368)+0.5*B(383)+0.5*B(385)+0.484*B(387)

  JVS(40) = B(301)+B(305)

  JVS(41) = -B(86)

  JVS(42) = -B(87)

  JVS(43) = -B(28)-B(30)-B(32)

  JVS(44) = B(2)

  JVS(45) = -B(100)

  JVS(46) = -B(101)

  JVS(47) = -B(118)

  JVS(48) = -B(119)

  JVS(49) = -B(120)

  JVS(50) = -B(121)

  JVS(51) = -B(82)

  JVS(52) = B(80)

  JVS(53) = B(81)

  JVS(54) = -B(102)

  JVS(55) = -B(103)

  JVS(56) = -B(98)

  JVS(57) = -B(99)

  JVS(58) = -B(96)

  JVS(59) = 0.03*B(186)

  JVS(60) = 0.06*B(188)

  JVS(61) = 0.03*B(187)+0.06*B(189)

  JVS(62) = -B(97)

  JVS(63) = -B(94)

  JVS(64) = 0.06*B(186)

  JVS(65) = 0.07*B(188)

  JVS(66) = 0.06*B(187)+0.07*B(189)

  JVS(67) = -B(95)

  JVS(68) = -B(140)

  JVS(69) = 0.35*B(138)

  JVS(70) = 0.35*B(139)-B(141)

  JVS(71) = -B(6)-B(61)-B(68)

  JVS(72) = B(59)

  JVS(73) = -B(69)

  JVS(74) = B(60)

  JVS(75) = -B(13)-B(142)

  JVS(76) = B(270)

  JVS(77) = B(271)

  JVS(78) = -B(143)

  JVS(79) = -B(4)-B(64)

  JVS(80) = B(208)

  JVS(81) = B(214)

  JVS(82) = B(220)

  JVS(83) = B(51)-B(65)

  JVS(84) = B(52)

  JVS(85) = B(209)+B(215)+B(221)

  JVS(86) = -B(9)-B(40)

  JVS(87) = 0.02*B(194)

  JVS(88) = 0.02*B(196)

  JVS(89) = 0.001*B(190)

  JVS(90) = 0.001*B(192)

  JVS(91) = 0.006*B(186)

  JVS(92) = 0.011*B(188)

  JVS(93) = 0.006*B(187)+0.011*B(189)+0.001*B(191)+0.001*B(193)+0.02*B(195)+0.02*B(197)

  JVS(94) = B(42)+B(43)

  JVS(95) = -B(41)

  JVS(96) = -B(204)-B(206)

  JVS(97) = 0.1*B(122)+B(166)

  JVS(98) = -B(207)

  JVS(99) = 0.1*B(123)

  JVS(100) = B(167)

  JVS(101) = -B(205)

  JVS(102) = 0.9*B(118)

  JVS(103) = -B(208)-B(210)-B(212)

  JVS(104) = -B(213)

  JVS(105) = 0.9*B(119)

  JVS(106) = -B(209)

  JVS(107) = 0.9*B(120)

  JVS(108) = -B(214)-B(216)-B(218)

  JVS(109) = -B(219)

  JVS(110) = 0.9*B(121)

  JVS(111) = -B(215)

  JVS(112) = 0.024*B(102)

  JVS(113) = -B(23)-B(130)

  JVS(114) = 0.6*B(150)

  JVS(115) = 0.41*B(136)

  JVS(116) = 0.3*B(347)

  JVS(117) = 0.3*B(348)

  JVS(118) = 0.024*B(103)-B(131)+0.41*B(137)+0.6*B(151)

  JVS(119) = -B(104)-B(168)-B(184)

  JVS(120) = -B(185)

  JVS(121) = -B(105)

  JVS(122) = -B(169)

  JVS(123) = -B(220)-B(222)-B(224)

  JVS(124) = 0.85*B(122)

  JVS(125) = -B(225)

  JVS(126) = 0.85*B(123)

  JVS(127) = -B(221)

  JVS(128) = -B(15)-B(146)

  JVS(129) = 0.11*B(200)

  JVS(130) = 0.11*B(201)

  JVS(131) = B(298)

  JVS(132) = B(299)

  JVS(133) = -B(147)

  JVS(134) = -B(5)-B(66)

  JVS(135) = B(166)

  JVS(136) = B(158)

  JVS(137) = B(160)

  JVS(138) = B(154)

  JVS(139) = 0.2*B(162)

  JVS(140) = B(156)

  JVS(141) = 0.5*B(164)

  JVS(142) = 0.3*B(62)

  JVS(143) = B(53)-B(67)

  JVS(144) = 0.3*B(63)+B(155)+B(157)+B(159)+B(161)+0.2*B(163)+0.5*B(165)+B(167)

  JVS(145) = B(54)

  JVS(146) = 0.036*B(98)

  JVS(147) = 0.43*B(184)

  JVS(148) = -B(88)

  JVS(149) = 0.14*B(194)

  JVS(150) = 0.14*B(196)

  JVS(151) = 0.36*B(190)

  JVS(152) = 1.87*B(17)+1.55*B(18)+2*B(132)+2*B(158)

  JVS(153) = 0.13*B(202)

  JVS(154) = 0.01*B(90)+0.36*B(192)

  JVS(155) = B(19)+B(134)+B(160)

  JVS(156) = B(10)+B(11)+B(124)+B(154)

  JVS(157) = B(22)+0.41*B(136)+0.8*B(162)+0.54*B(198)

  JVS(158) = B(12)

  JVS(159) = 0.66*B(200)

  JVS(160) = 0.37*B(186)

  JVS(161) = 0.3*B(188)

  JVS(162) = 0.01*B(91)

  JVS(163) = 0.43*B(185)+0.37*B(187)+0.3*B(189)+0.36*B(191)+0.36*B(193)+0.14*B(195)+0.14*B(197)+0.54*B(199)+0.66*B(201)&
               &+0.13*B(203)

  JVS(164) = -B(89)+0.036*B(99)+B(125)+2*B(133)+B(135)+0.41*B(137)

  JVS(165) = B(155)+2*B(159)+B(161)+0.8*B(163)

  JVS(166) = -B(114)-B(178)-B(194)

  JVS(167) = -B(195)

  JVS(168) = -B(115)

  JVS(169) = -B(179)

  JVS(170) = -B(116)-B(180)-B(196)

  JVS(171) = -B(197)

  JVS(172) = -B(117)

  JVS(173) = -B(181)

  JVS(174) = -B(148)-B(228)

  JVS(175) = 0.4*B(150)+0.4*B(182)+0.3*B(202)

  JVS(176) = 0.3*B(203)

  JVS(177) = B(226)

  JVS(178) = -B(149)+0.4*B(151)

  JVS(179) = 0.4*B(183)

  JVS(180) = B(227)

  JVS(181) = 0.1*B(204)+B(206)

  JVS(182) = B(208)+0.02*B(210)+B(212)

  JVS(183) = B(214)+0.02*B(216)+B(218)

  JVS(184) = B(220)+0.02*B(222)+B(224)

  JVS(185) = -B(122)-B(166)

  JVS(186) = B(213)+B(219)+B(225)

  JVS(187) = B(207)

  JVS(188) = -B(123)

  JVS(189) = -B(167)

  JVS(190) = 0.1*B(205)+B(209)+B(215)+B(221)

  JVS(191) = -B(110)-B(174)-B(190)

  JVS(192) = -B(191)

  JVS(193) = -B(111)

  JVS(194) = -B(175)

  JVS(195) = 0.036*B(98)

  JVS(196) = -B(17)-B(18)-B(132)-B(158)

  JVS(197) = B(258)+B(337)+B(377)+B(418)

  JVS(198) = 1.2*B(254)+0.65*B(333)+0.65*B(373)+1.3*B(414)

  JVS(199) = 0.35*B(256)+0.37*B(335)+0.37*B(375)+0.74*B(416)

  JVS(200) = 0.063*B(236)+0.119*B(315)+0.1*B(355)+0.063*B(396)

  JVS(201) = 0.15*B(138)+0.25*B(164)+0.5*B(200)

  JVS(202) = 0.5*B(201)

  JVS(203) = 0.119*B(316)+0.65*B(334)+0.37*B(336)+B(338)

  JVS(204) = 0.1*B(356)+0.65*B(374)+0.37*B(376)+B(378)

  JVS(205) = 0.036*B(99)-B(133)+0.15*B(139)

  JVS(206) = -B(159)+0.25*B(165)+0.063*B(397)+1.3*B(415)+0.74*B(417)+B(419)

  JVS(207) = 0.063*B(237)+1.2*B(255)+0.35*B(257)+B(259)

  JVS(208) = -B(150)-B(182)-B(202)-B(231)

  JVS(209) = B(229)

  JVS(210) = -B(203)

  JVS(211) = -B(151)

  JVS(212) = -B(183)

  JVS(213) = B(230)

  JVS(214) = B(104)

  JVS(215) = -B(242)-B(280)-B(321)-B(361)-B(402)

  JVS(216) = 0

  JVS(217) = -B(322)

  JVS(218) = -B(362)

  JVS(219) = -B(281)

  JVS(220) = B(105)

  JVS(221) = -B(403)

  JVS(222) = -B(243)

  JVS(223) = -B(90)-B(112)-B(176)-B(192)

  JVS(224) = -B(91)

  JVS(225) = -B(193)

  JVS(226) = -B(113)

  JVS(227) = -B(177)

  JVS(228) = -B(244)-B(282)-B(323)-B(363)-B(404)

  JVS(229) = B(106)

  JVS(230) = -B(324)

  JVS(231) = -B(364)

  JVS(232) = -B(283)

  JVS(233) = B(107)

  JVS(234) = -B(405)

  JVS(235) = -B(245)

  JVS(236) = -B(246)-B(284)-B(325)-B(365)-B(406)

  JVS(237) = B(108)

  JVS(238) = -B(326)

  JVS(239) = -B(366)

  JVS(240) = -B(285)

  JVS(241) = B(109)

  JVS(242) = -B(407)

  JVS(243) = -B(247)

  JVS(244) = B(130)

  JVS(245) = 0

  JVS(246) = -B(19)-B(134)-B(160)

  JVS(247) = B(258)+B(337)+B(377)+B(418)

  JVS(248) = 0.65*B(254)+0.35*B(333)+0.35*B(373)+0.7*B(414)

  JVS(249) = 0.6*B(256)+0.63*B(335)+0.63*B(375)+1.26*B(416)

  JVS(250) = 0.08*B(136)+0.6*B(198)

  JVS(251) = 0.005*B(315)+0.004*B(355)

  JVS(252) = 0.15*B(138)+0.25*B(164)+0.62*B(200)

  JVS(253) = 0

  JVS(254) = 0.6*B(199)+0.62*B(201)

  JVS(255) = 0.54*B(264)+0.4*B(347)+0.54*B(382)+0.54*B(424)

  JVS(256) = 0.005*B(316)+0.35*B(334)+0.63*B(336)+B(338)+0.4*B(348)

  JVS(257) = 0.004*B(356)+0.35*B(374)+0.63*B(376)+B(378)+0.54*B(383)

  JVS(258) = B(131)-B(135)+0.08*B(137)+0.15*B(139)

  JVS(259) = -B(161)+0.25*B(165)+0.7*B(415)+1.26*B(417)+B(419)+0.54*B(425)

  JVS(260) = 0.65*B(255)+0.6*B(257)+B(259)+0.54*B(265)

  JVS(261) = 0

  JVS(262) = 0.98*B(222)

  JVS(263) = 0

  JVS(264) = -B(258)-B(296)-B(337)-B(377)-B(418)

  JVS(265) = 0

  JVS(266) = -B(338)

  JVS(267) = -B(378)

  JVS(268) = -B(297)

  JVS(269) = 0

  JVS(270) = -B(419)

  JVS(271) = -B(259)

  JVS(272) = 0

  JVS(273) = 0.25*B(100)

  JVS(274) = 0.12*B(140)

  JVS(275) = 0.53*B(194)

  JVS(276) = 0.06*B(244)+0.081*B(323)+0.141*B(363)+0.06*B(404)

  JVS(277) = 0.29*B(246)+0.313*B(325)+0.569*B(365)+0.29*B(406)

  JVS(278) = -B(16)-B(128)

  JVS(279) = 0.722*B(238)+0.24*B(317)+0.33*B(357)+0.828*B(398)

  JVS(280) = 0.642*B(240)+0.419*B(319)+0.581*B(359)+0.88*B(400)

  JVS(281) = 0.8*B(250)+B(329)+B(369)+B(410)

  JVS(282) = 0.623*B(236)+0.018*B(315)+0.127*B(355)+0.67*B(396)

  JVS(283) = 0.03*B(164)

  JVS(284) = 0.03*B(186)

  JVS(285) = 0.16*B(188)

  JVS(286) = 0.149*B(389)

  JVS(287) = 0.464*B(268)+0.149*B(351)+0.167*B(386)+0.149*B(390)+0.285*B(391)+0.469*B(428)

  JVS(288) = 0.03*B(187)+0.16*B(189)+0.53*B(195)

  JVS(289) = 0.11*B(382)

  JVS(290) = 0.018*B(316)+0.24*B(318)+0.419*B(320)+0.081*B(324)+0.313*B(326)+B(330)+0.149*B(352)

  JVS(291) = 0.127*B(356)+0.33*B(358)+0.581*B(360)+0.141*B(364)+0.569*B(366)+B(370)+0.11*B(383)+0.167*B(387)

  JVS(292) = 0

  JVS(293) = 0.25*B(101)-B(129)+0.12*B(141)+0.41*B(144)

  JVS(294) = 0.8*B(21)

  JVS(295) = 0.03*B(165)+0.67*B(397)+0.828*B(399)+0.88*B(401)+0.06*B(405)+0.29*B(407)+B(411)+0.469*B(429)

  JVS(296) = 0.41*B(145)

  JVS(297) = 0.623*B(237)+0.722*B(239)+0.642*B(241)+0.06*B(245)+0.29*B(247)+0.8*B(251)+0.464*B(269)

  JVS(298) = B(116)

  JVS(299) = -B(252)-B(290)-B(331)-B(371)-B(412)

  JVS(300) = 0

  JVS(301) = -B(332)

  JVS(302) = -B(372)

  JVS(303) = -B(291)

  JVS(304) = B(117)

  JVS(305) = -B(413)

  JVS(306) = -B(253)

  JVS(307) = 0.75*B(100)

  JVS(308) = -B(238)-B(276)-B(317)-B(357)-B(398)

  JVS(309) = -B(318)

  JVS(310) = -B(358)

  JVS(311) = -B(277)

  JVS(312) = 0.75*B(101)

  JVS(313) = -B(399)

  JVS(314) = -B(239)

  JVS(315) = 0.9511*B(102)

  JVS(316) = -B(240)-B(278)-B(319)-B(359)-B(400)

  JVS(317) = -B(320)

  JVS(318) = -B(360)

  JVS(319) = -B(279)

  JVS(320) = 0.9511*B(103)

  JVS(321) = -B(401)

  JVS(322) = -B(241)

  JVS(323) = 0.98*B(210)

  JVS(324) = -B(254)-B(292)-B(333)-B(373)-B(414)

  JVS(325) = 0

  JVS(326) = -B(334)

  JVS(327) = -B(374)

  JVS(328) = -B(293)

  JVS(329) = 0

  JVS(330) = -B(415)

  JVS(331) = -B(255)

  JVS(332) = 0

  JVS(333) = 0.98*B(216)

  JVS(334) = -B(256)-B(294)-B(335)-B(375)-B(416)

  JVS(335) = 0

  JVS(336) = -B(336)

  JVS(337) = -B(376)

  JVS(338) = -B(295)

  JVS(339) = 0

  JVS(340) = -B(417)

  JVS(341) = -B(257)

  JVS(342) = 0

  JVS(343) = 0.01*B(98)

  JVS(344) = B(13)+0.35*B(142)

  JVS(345) = B(23)

  JVS(346) = B(184)

  JVS(347) = 0.35*B(146)

  JVS(348) = 0.04*B(196)

  JVS(349) = B(148)

  JVS(350) = 0.9*B(190)

  JVS(351) = 0.13*B(17)+0.45*B(18)

  JVS(352) = 0.4*B(150)+0.4*B(182)+0.7*B(202)

  JVS(353) = 1.6*B(242)+1.55*B(321)+0.8*B(361)+1.6*B(402)

  JVS(354) = 0.05*B(90)+0.9*B(192)

  JVS(355) = B(244)+1.25*B(323)+0.501*B(363)+B(404)

  JVS(356) = 0.755*B(325)

  JVS(357) = B(337)

  JVS(358) = 0.25*B(252)+1.4*B(331)+0.4*B(371)+0.4*B(412)

  JVS(359) = 0.021*B(238)+0.829*B(317)+0.076*B(357)+0.021*B(398)

  JVS(360) = 0.753*B(319)

  JVS(361) = B(333)

  JVS(362) = B(335)

  JVS(363) = -B(10)-B(11)-B(124)-B(154)

  JVS(364) = B(329)

  JVS(365) = 0.606*B(248)+1.09*B(327)+0.34*B(367)+0.686*B(408)

  JVS(366) = B(22)+0.08*B(136)+0.4*B(198)

  JVS(367) = 0.047*B(236)+0.81*B(315)+0.091*B(355)+0.048*B(396)

  JVS(368) = 0

  JVS(369) = B(262)+2*B(343)+B(345)+B(380)+B(422)

  JVS(370) = B(432)

  JVS(371) = 0.64*B(186)

  JVS(372) = 0.02*B(188)

  JVS(373) = 0.75*B(349)+0.202*B(389)

  JVS(374) = 0.287*B(268)+0.96*B(351)+0.207*B(386)+0.202*B(390)+0.504*B(391)+0.28*B(428)

  JVS(375) = 0.75*B(313)

  JVS(376) = 0.05*B(91)

  JVS(377) = B(185)+0.64*B(187)+0.02*B(189)+0.9*B(191)+0.9*B(193)+0.04*B(197)+0.4*B(199)+0.7*B(203)

  JVS(378) = 0.75*B(347)

  JVS(379) = B(232)+1.33*B(312)+0.75*B(314)+0.81*B(316)+0.829*B(318)+0.753*B(320)+1.55*B(322)+1.25*B(324)+0.755*B(326)&
               &+1.09*B(328)+B(330)+1.4*B(332)+B(334)+B(336)+B(338)+B(339)+B(341)+2*B(344)+B(346)+0.75*B(348)+0.75*B(350)&
               &+0.96*B(352)+B(392)+B(433)

  JVS(380) = B(340)+B(342)+0.091*B(356)+0.076*B(358)+0.8*B(362)+0.501*B(364)+0.34*B(368)+0.4*B(372)+B(381)+0.207*B(387)

  JVS(381) = 0

  JVS(382) = 0.01*B(99)-B(125)+0.08*B(137)+0.35*B(143)+0.35*B(147)+B(149)+0.4*B(151)

  JVS(383) = -B(155)+0.4*B(183)+B(393)+0.048*B(397)+0.021*B(399)+1.6*B(403)+B(405)+0.686*B(409)+0.4*B(413)+B(423)+0.28&
               &*B(429)

  JVS(384) = B(233)+0.047*B(237)+0.021*B(239)+1.6*B(243)+B(245)+0.606*B(249)+0.25*B(253)+B(263)+0.287*B(269)

  JVS(385) = 0

  JVS(386) = B(114)

  JVS(387) = -B(250)-B(288)-B(329)-B(369)-B(410)

  JVS(388) = 0

  JVS(389) = -B(330)

  JVS(390) = -B(370)

  JVS(391) = -B(289)

  JVS(392) = B(115)

  JVS(393) = -B(411)

  JVS(394) = -B(251)

  JVS(395) = B(110)

  JVS(396) = B(112)

  JVS(397) = -B(248)-B(286)-B(327)-B(367)-B(408)

  JVS(398) = 0

  JVS(399) = 0

  JVS(400) = -B(328)

  JVS(401) = -B(368)

  JVS(402) = -B(287)

  JVS(403) = B(111)+B(113)

  JVS(404) = -B(409)

  JVS(405) = -B(249)

  JVS(406) = 0.79*B(196)

  JVS(407) = 0.9*B(174)+0.39*B(190)

  JVS(408) = 0.9*B(176)+0.39*B(192)

  JVS(409) = 0.4*B(252)+0.6*B(331)+0.6*B(371)+0.6*B(412)

  JVS(410) = 0.446*B(248)+0.55*B(327)+0.771*B(367)+0.6*B(408)

  JVS(411) = -B(22)-B(92)-B(136)-B(162)-B(198)

  JVS(412) = -B(93)

  JVS(413) = 0.39*B(191)+0.39*B(193)+0.79*B(197)-B(199)

  JVS(414) = 0.55*B(328)+0.6*B(332)

  JVS(415) = 0.771*B(368)+0.6*B(372)

  JVS(416) = 0

  JVS(417) = -B(137)

  JVS(418) = -B(163)+0.9*B(175)+0.9*B(177)+0.6*B(409)+0.6*B(413)

  JVS(419) = 0.446*B(249)+0.4*B(253)

  JVS(420) = 0.583*B(98)

  JVS(421) = -B(236)-B(274)-B(315)-B(355)-B(396)

  JVS(422) = -B(316)

  JVS(423) = -B(356)

  JVS(424) = -B(275)

  JVS(425) = 0.583*B(99)+0.44*B(144)+B(152)

  JVS(426) = B(153)

  JVS(427) = -B(397)

  JVS(428) = 0.44*B(145)

  JVS(429) = -B(237)

  JVS(430) = 0.025*B(102)

  JVS(431) = 0.335*B(98)

  JVS(432) = 0.88*B(140)

  JVS(433) = 0.65*B(194)

  JVS(434) = 0.2*B(242)+0.35*B(321)+0.6*B(361)+0.2*B(402)

  JVS(435) = 0.94*B(244)+0.669*B(323)+0.859*B(363)+0.94*B(404)

  JVS(436) = 1.71*B(246)+0.932*B(325)+0.941*B(365)+1.71*B(406)

  JVS(437) = 0.211*B(238)+0.523*B(317)+0.677*B(357)+0.239*B(398)

  JVS(438) = 0.15*B(240)+0.411*B(319)+0.497*B(359)+0.187*B(400)

  JVS(439) = 0.8*B(250)+B(329)+B(369)+B(410)

  JVS(440) = B(92)

  JVS(441) = 0.233*B(236)+0.58*B(315)+0.724*B(355)+0.243*B(396)

  JVS(442) = -B(12)-B(126)-B(156)

  JVS(443) = 0.25*B(164)+0.16*B(200)

  JVS(444) = 0.44*B(186)

  JVS(445) = 0.99*B(188)

  JVS(446) = 0.64*B(389)

  JVS(447) = 1.24*B(268)+0.64*B(351)+0.65*B(386)+0.64*B(390)+1.21*B(391)+1.24*B(428)

  JVS(448) = B(234)+0.75*B(313)+B(353)+B(394)

  JVS(449) = B(93)

  JVS(450) = 0.44*B(187)+0.99*B(189)+0.65*B(195)+0.16*B(201)

  JVS(451) = 0.46*B(264)+0.3*B(347)+0.35*B(382)+0.46*B(424)

  JVS(452) = 0.75*B(314)+0.58*B(316)+0.523*B(318)+0.411*B(320)+0.35*B(322)+0.669*B(324)+0.932*B(326)+B(330)+0.3*B(348)&
               &+0.64*B(352)

  JVS(453) = B(354)+0.724*B(356)+0.677*B(358)+0.497*B(360)+0.6*B(362)+0.859*B(364)+0.941*B(366)+B(370)+0.35*B(383)+0.65&
               &*B(387)

  JVS(454) = 0

  JVS(455) = 0.335*B(99)+0.025*B(103)-B(127)+0.88*B(141)+0.08*B(144)

  JVS(456) = 0.2*B(21)

  JVS(457) = -B(157)+0.25*B(165)+B(395)+0.243*B(397)+0.239*B(399)+0.187*B(401)+0.2*B(403)+0.94*B(405)+1.71*B(407)+B(411)&
               &+0.46*B(425)+1.24*B(429)

  JVS(458) = B(14)+0.08*B(145)

  JVS(459) = B(235)+0.233*B(237)+0.211*B(239)+0.15*B(241)+0.2*B(243)+0.94*B(245)+1.71*B(247)+0.8*B(251)+0.46*B(265)+1.24&
               &*B(269)

  JVS(460) = 0.13*B(90)

  JVS(461) = 0.5*B(254)+B(333)+B(373)+0.5*B(414)

  JVS(462) = 0.95*B(256)+B(335)+B(375)+B(416)

  JVS(463) = -B(20)-B(138)-B(164)-B(200)

  JVS(464) = 0.13*B(91)

  JVS(465) = -B(201)

  JVS(466) = B(334)+B(336)

  JVS(467) = B(374)+B(376)

  JVS(468) = 0

  JVS(469) = -B(139)

  JVS(470) = -B(165)+0.5*B(415)+B(417)

  JVS(471) = 0.5*B(255)+0.95*B(257)

  JVS(472) = 0

  JVS(473) = B(231)

  JVS(474) = 0.51*B(136)+0.2*B(162)

  JVS(475) = B(20)+0.5*B(138)+0.5*B(164)

  JVS(476) = -B(229)-B(262)-B(302)-B(304)-B(343)-B(345)-B(380)-B(422)

  JVS(477) = 0

  JVS(478) = 0

  JVS(479) = -B(344)-B(346)

  JVS(480) = -B(381)

  JVS(481) = -B(303)-B(305)

  JVS(482) = 0.51*B(137)+0.5*B(139)

  JVS(483) = 0.2*B(163)+0.5*B(165)-B(423)

  JVS(484) = -B(263)

  JVS(485) = -B(230)

  JVS(486) = 0.1*B(118)

  JVS(487) = 0.1*B(120)

  JVS(488) = 0.35*B(146)

  JVS(489) = B(148)

  JVS(490) = 0.05*B(122)

  JVS(491) = 0.13*B(190)

  JVS(492) = B(150)+B(182)

  JVS(493) = 0.15*B(90)+0.13*B(192)

  JVS(494) = 0.334*B(238)+0.245*B(317)+0.237*B(357)+0.391*B(398)

  JVS(495) = 0.416*B(240)+0.322*B(319)+0.318*B(359)+0.587*B(400)

  JVS(496) = 0.49*B(136)

  JVS(497) = 0.048*B(236)+0.085*B(315)+0.071*B(355)+0.051*B(396)

  JVS(498) = 0.5*B(138)+0.5*B(164)

  JVS(499) = 0

  JVS(500) = -B(430)-B(432)-B(434)-2*B(436)-B(437)-B(439)

  JVS(501) = 0.15*B(91)

  JVS(502) = 0.13*B(191)+0.13*B(193)

  JVS(503) = 0.16*B(264)+0.08*B(347)+0.08*B(382)+0.16*B(424)

  JVS(504) = 0.085*B(316)+0.245*B(318)+0.322*B(320)+0.08*B(348)-B(433)

  JVS(505) = 0.071*B(356)+0.237*B(358)+0.318*B(360)+0.08*B(383)-B(435)

  JVS(506) = -B(431)

  JVS(507) = 0.1*B(119)+0.1*B(121)+0.05*B(123)+0.49*B(137)+0.5*B(139)+0.07*B(144)+0.35*B(147)+B(149)+B(151)

  JVS(508) = 0

  JVS(509) = 0.5*B(165)+B(183)+0.051*B(397)+0.391*B(399)+0.587*B(401)+0.16*B(425)-B(440)

  JVS(510) = 0.07*B(145)

  JVS(511) = 0.048*B(237)+0.334*B(239)+0.416*B(241)+0.16*B(265)-B(438)

  JVS(512) = 0

  JVS(513) = 0.46*B(196)

  JVS(514) = 0.35*B(190)

  JVS(515) = 0.86*B(90)+0.35*B(192)

  JVS(516) = 0.354*B(248)+0.37*B(327)+0.229*B(367)+0.4*B(408)

  JVS(517) = -B(106)-B(170)-B(186)

  JVS(518) = 0.86*B(91)

  JVS(519) = -B(187)+0.35*B(191)+0.35*B(193)+0.46*B(197)

  JVS(520) = 0.37*B(328)

  JVS(521) = 0.229*B(368)

  JVS(522) = 0

  JVS(523) = -B(107)

  JVS(524) = -B(171)+0.4*B(409)

  JVS(525) = 0.354*B(249)

  JVS(526) = 0.25*B(252)+0.4*B(331)+0.4*B(371)+0.4*B(412)

  JVS(527) = 0.08*B(327)

  JVS(528) = -B(108)-B(172)-B(188)

  JVS(529) = 0

  JVS(530) = -B(189)

  JVS(531) = 0.08*B(328)+0.4*B(332)

  JVS(532) = 0.4*B(372)

  JVS(533) = 0

  JVS(534) = -B(109)

  JVS(535) = -B(173)+0.4*B(413)

  JVS(536) = 0.25*B(253)

  JVS(537) = 0.8*B(168)

  JVS(538) = 0.1*B(178)

  JVS(539) = 0.13*B(180)

  JVS(540) = 0.9*B(174)

  JVS(541) = 0.9*B(176)

  JVS(542) = 0.8*B(162)

  JVS(543) = 0.43*B(170)

  JVS(544) = 0.11*B(172)

  JVS(545) = -B(266)-B(308)-B(349)-B(384)-2*B(388)-B(389)-B(426)

  JVS(546) = -B(390)

  JVS(547) = 0

  JVS(548) = 0

  JVS(549) = -B(350)

  JVS(550) = -B(385)

  JVS(551) = -B(309)

  JVS(552) = 0

  JVS(553) = 0.8*B(163)+0.8*B(169)+0.43*B(171)+0.11*B(173)+0.9*B(175)+0.9*B(177)+0.1*B(179)+0.13*B(181)-B(427)

  JVS(554) = -B(267)

  JVS(555) = 0.2*B(168)

  JVS(556) = 0.9*B(178)

  JVS(557) = 0.87*B(180)

  JVS(558) = 0.1*B(174)

  JVS(559) = 0.1*B(176)

  JVS(560) = 0.57*B(170)

  JVS(561) = 0.89*B(172)

  JVS(562) = -B(389)

  JVS(563) = -B(268)-B(310)-B(351)-B(386)-B(390)-2*B(391)-B(428)

  JVS(564) = 0

  JVS(565) = 0

  JVS(566) = -B(352)

  JVS(567) = -B(387)

  JVS(568) = -B(311)

  JVS(569) = 0

  JVS(570) = 0.2*B(169)+0.57*B(171)+0.89*B(173)+0.1*B(175)+0.1*B(177)+0.9*B(179)+0.87*B(181)-B(429)

  JVS(571) = -B(269)

  JVS(572) = B(96)

  JVS(573) = 0.2*B(194)

  JVS(574) = 0.16*B(196)

  JVS(575) = B(16)

  JVS(576) = 0.245*B(238)+0.014*B(317)+0.018*B(357)+0.262*B(398)

  JVS(577) = 0.133*B(240)+0.013*B(319)+0.015*B(359)+0.155*B(400)

  JVS(578) = 0

  JVS(579) = 0.048*B(236)+0.006*B(355)+0.053*B(396)

  JVS(580) = 0

  JVS(581) = 0.1*B(186)

  JVS(582) = 0.18*B(188)

  JVS(583) = 0

  JVS(584) = 0

  JVS(585) = -B(234)-B(272)-B(313)-B(353)-B(394)

  JVS(586) = 0

  JVS(587) = 0.1*B(187)+0.18*B(189)+0.2*B(195)+0.16*B(197)

  JVS(588) = 0

  JVS(589) = -B(314)+0.014*B(318)+0.013*B(320)

  JVS(590) = -B(354)+0.006*B(356)+0.018*B(358)+0.015*B(360)

  JVS(591) = -B(273)

  JVS(592) = B(97)

  JVS(593) = 0

  JVS(594) = -B(395)+0.053*B(397)+0.262*B(399)+0.155*B(401)

  JVS(595) = 0

  JVS(596) = -B(235)+0.048*B(237)+0.245*B(239)+0.133*B(241)

  JVS(597) = 0

  JVS(598) = B(30)

  JVS(599) = 0.09*B(190)

  JVS(600) = -B(90)+0.09*B(192)

  JVS(601) = -B(92)

  JVS(602) = -B(24)-B(26)-B(45)-B(47)-B(49)-B(91)-B(93)

  JVS(603) = B(3)-B(27)+0.09*B(191)+0.09*B(193)

  JVS(604) = 0

  JVS(605) = 0

  JVS(606) = 0

  JVS(607) = 0

  JVS(608) = B(8)

  JVS(609) = -B(46)

  JVS(610) = B(1)-B(48)-B(50)

  JVS(611) = B(28)

  JVS(612) = -B(212)

  JVS(613) = -B(218)

  JVS(614) = -B(184)

  JVS(615) = -B(224)

  JVS(616) = -B(194)

  JVS(617) = -B(196)

  JVS(618) = 0

  JVS(619) = -B(190)

  JVS(620) = -B(202)

  JVS(621) = -B(192)

  JVS(622) = -B(198)

  JVS(623) = -B(200)

  JVS(624) = B(304)

  JVS(625) = -B(186)

  JVS(626) = -B(188)

  JVS(627) = B(24)-B(26)

  JVS(628) = -B(2)-B(3)-B(27)-B(34)-B(36)-B(70)-B(72)-B(185)-B(187)-B(189)-B(191)-B(193)-B(195)-B(197)-B(199)-B(201)&
               &-B(203)-B(213)-B(219)-B(225)

  JVS(629) = 0

  JVS(630) = B(300)

  JVS(631) = -B(37)+B(301)+B(305)

  JVS(632) = -B(35)

  JVS(633) = 0

  JVS(634) = -B(71)

  JVS(635) = -B(73)

  JVS(636) = 0.42*B(194)

  JVS(637) = 0.42*B(196)

  JVS(638) = 0.02*B(190)

  JVS(639) = 0.02*B(192)

  JVS(640) = B(128)

  JVS(641) = 0

  JVS(642) = 0

  JVS(643) = 0

  JVS(644) = 0

  JVS(645) = 0

  JVS(646) = 0.03*B(186)

  JVS(647) = 0.12*B(188)

  JVS(648) = 0

  JVS(649) = 0

  JVS(650) = 0

  JVS(651) = 0.03*B(187)+0.12*B(189)+0.02*B(191)+0.02*B(193)+0.42*B(195)+0.42*B(197)

  JVS(652) = -B(264)-B(306)-B(347)-B(382)-B(424)

  JVS(653) = -B(348)

  JVS(654) = -B(383)

  JVS(655) = -B(307)

  JVS(656) = B(129)

  JVS(657) = 0

  JVS(658) = -B(425)

  JVS(659) = 0

  JVS(660) = -B(265)

  JVS(661) = 0

  JVS(662) = B(94)

  JVS(663) = 0.65*B(142)

  JVS(664) = B(15)

  JVS(665) = 0.03*B(190)

  JVS(666) = -B(321)+0.5*B(361)

  JVS(667) = 0.03*B(192)

  JVS(668) = -B(323)+0.501*B(363)

  JVS(669) = -B(325)+0.51*B(365)

  JVS(670) = -B(337)+B(377)

  JVS(671) = -B(331)+B(371)

  JVS(672) = 0.031*B(238)-0.951*B(317)+0.554*B(357)+0.04*B(398)

  JVS(673) = -B(319)+0.507*B(359)

  JVS(674) = -B(333)+B(373)

  JVS(675) = -B(335)+B(375)

  JVS(676) = -B(329)+B(369)

  JVS(677) = -B(327)+0.506*B(367)

  JVS(678) = 0.15*B(236)-0.993*B(315)+0.508*B(355)+0.155*B(396)

  JVS(679) = B(12)

  JVS(680) = 0

  JVS(681) = -B(343)-B(345)+B(380)

  JVS(682) = -B(432)+B(434)

  JVS(683) = 0.19*B(186)

  JVS(684) = 0.23*B(188)

  JVS(685) = -B(349)+0.5*B(384)

  JVS(686) = -B(351)+0.516*B(386)

  JVS(687) = -B(313)+0.5*B(353)

  JVS(688) = 0

  JVS(689) = 0.19*B(187)+0.23*B(189)+0.03*B(191)+0.03*B(193)

  JVS(690) = -B(347)+0.5*B(382)

  JVS(691) = -B(232)-B(270)-2*B(312)-B(314)-0.993*B(316)-0.951*B(318)-B(320)-B(322)-B(324)-B(326)-B(328)-B(330)-B(332)&
               &-B(334)-B(336)-B(338)-B(341)-B(344)-B(346)-B(348)-B(350)-B(352)-B(392)-B(433)

  JVS(692) = B(260)-B(342)+0.5*B(354)+0.508*B(356)+0.554*B(358)+0.507*B(360)+0.5*B(362)+0.501*B(364)+0.51*B(366)+0.506&
               &*B(368)+B(370)+B(372)+B(374)+B(376)+B(378)+2*B(379)+B(381)+0.5*B(383)+0.5*B(385)+0.516*B(387)+B(420)+B(435)

  JVS(693) = -B(271)

  JVS(694) = B(95)+0.65*B(143)

  JVS(695) = 0

  JVS(696) = -B(393)+0.155*B(397)+0.04*B(399)+B(421)

  JVS(697) = 0

  JVS(698) = -B(233)+0.15*B(237)+0.031*B(239)+B(261)

  JVS(699) = 0

  JVS(700) = B(23)

  JVS(701) = 0.65*B(146)

  JVS(702) = B(228)

  JVS(703) = 0.15*B(190)

  JVS(704) = 0.7*B(202)

  JVS(705) = -B(361)

  JVS(706) = 0.15*B(192)

  JVS(707) = -B(363)

  JVS(708) = -B(365)

  JVS(709) = B(19)+B(134)+B(160)

  JVS(710) = -B(377)

  JVS(711) = B(16)

  JVS(712) = -B(371)

  JVS(713) = -B(357)

  JVS(714) = -B(359)

  JVS(715) = -B(373)

  JVS(716) = -B(375)

  JVS(717) = -B(369)

  JVS(718) = -B(367)

  JVS(719) = B(22)+0.13*B(198)

  JVS(720) = -B(355)

  JVS(721) = B(126)+B(156)

  JVS(722) = 0.28*B(200)

  JVS(723) = B(262)+B(343)+B(422)

  JVS(724) = -B(434)

  JVS(725) = 0

  JVS(726) = 0

  JVS(727) = -B(384)

  JVS(728) = -B(386)

  JVS(729) = -B(353)

  JVS(730) = 0

  JVS(731) = 0.15*B(191)+0.15*B(193)+0.13*B(199)+0.28*B(201)+0.7*B(203)

  JVS(732) = 0.23*B(264)+0.12*B(347)-0.88*B(382)+0.23*B(424)

  JVS(733) = -B(339)-B(341)+B(344)+0.12*B(348)

  JVS(734) = -B(226)-B(260)-B(298)-B(300)-B(340)-B(342)-B(354)-B(356)-B(358)-B(360)-B(362)-B(364)-B(366)-B(368)-B(370)&
               &-B(372)-B(374)-B(376)-B(378)-2*B(379)-0.88*B(383)-B(385)-B(387)-B(420)-B(435)

  JVS(735) = -B(299)-B(301)

  JVS(736) = B(127)+B(135)+0.65*B(147)

  JVS(737) = 0

  JVS(738) = B(157)+B(161)-B(421)+B(423)+0.23*B(425)

  JVS(739) = 0

  JVS(740) = -B(261)+B(263)+0.23*B(265)

  JVS(741) = -B(227)

  JVS(742) = B(86)

  JVS(743) = 0.25*B(100)

  JVS(744) = 0.1*B(118)

  JVS(745) = 0.1*B(120)

  JVS(746) = 0.049*B(102)

  JVS(747) = 0.381*B(98)

  JVS(748) = B(140)

  JVS(749) = 0.65*B(6)+B(61)

  JVS(750) = B(13)

  JVS(751) = B(40)

  JVS(752) = -B(206)

  JVS(753) = 0.02*B(210)

  JVS(754) = 0.02*B(216)

  JVS(755) = B(23)+B(130)

  JVS(756) = 0.26*B(184)

  JVS(757) = 0.02*B(222)

  JVS(758) = 0.35*B(146)

  JVS(759) = B(88)

  JVS(760) = 0.1*B(194)

  JVS(761) = 0.1*B(196)

  JVS(762) = 0.05*B(122)

  JVS(763) = 0.3*B(190)

  JVS(764) = 0.8*B(18)+B(132)+B(158)

  JVS(765) = 0.4*B(150)+0.08*B(202)

  JVS(766) = B(242)-B(280)+B(321)+0.5*B(361)+B(402)

  JVS(767) = 0.28*B(90)+0.3*B(192)

  JVS(768) = B(244)-B(282)+B(323)+0.501*B(363)+B(404)

  JVS(769) = B(246)-B(284)+B(325)+0.51*B(365)+B(406)

  JVS(770) = B(19)

  JVS(771) = B(258)-B(296)+2*B(337)+B(377)+B(418)

  JVS(772) = 0.65*B(252)-B(290)+2*B(331)+B(371)+B(412)

  JVS(773) = 0.599*B(238)-B(276)+0.946*B(317)+0.438*B(357)+0.699*B(398)

  JVS(774) = 0.606*B(240)-B(278)+0.993*B(319)+0.489*B(359)+0.845*B(400)

  JVS(775) = 0.95*B(254)-B(292)+B(333)+B(373)+B(414)

  JVS(776) = 0.95*B(256)-B(294)+B(335)+B(375)+B(416)

  JVS(777) = 2*B(11)+B(124)+B(154)

  JVS(778) = 0.8*B(250)-B(288)+2*B(329)+B(369)+B(410)

  JVS(779) = 0.847*B(248)-B(286)+B(327)+0.506*B(367)+B(408)

  JVS(780) = B(22)+0.49*B(136)+0.29*B(198)

  JVS(781) = 0.742*B(236)-B(274)+0.992*B(315)+0.488*B(355)+0.792*B(396)

  JVS(782) = B(12)

  JVS(783) = B(20)+0.5*B(138)+0.5*B(164)+0.29*B(200)

  JVS(784) = -B(302)-B(304)+B(343)

  JVS(785) = -B(430)+B(432)

  JVS(786) = 0.25*B(186)

  JVS(787) = 0.22*B(188)

  JVS(788) = B(266)-B(308)+B(349)+0.5*B(384)+B(388)+0.5*B(389)+B(426)

  JVS(789) = -B(310)+0.5*B(351)+0.5*B(390)

  JVS(790) = B(234)-B(272)+B(313)+0.5*B(353)+B(394)

  JVS(791) = 0.28*B(91)

  JVS(792) = B(34)-B(36)+0.26*B(185)+0.25*B(187)+0.22*B(189)+0.3*B(191)+0.3*B(193)+0.1*B(195)+0.1*B(197)+0.29*B(199)&
               &+0.29*B(201)+0.08*B(203)

  JVS(793) = 0.77*B(264)-B(306)+0.88*B(347)+0.38*B(382)+0.77*B(424)

  JVS(794) = B(232)-B(270)+0.66*B(312)+B(314)+0.992*B(316)+0.946*B(318)+0.993*B(320)+B(322)+B(324)+B(326)+B(328)+2&
               &*B(330)+2*B(332)+B(334)+B(336)+2*B(338)+B(339)+B(344)+0.88*B(348)+B(350)+0.5*B(352)+B(392)+B(433)

  JVS(795) = -B(298)-B(300)+B(340)+0.5*B(354)+0.488*B(356)+0.438*B(358)+0.489*B(360)+0.5*B(362)+0.501*B(364)+0.51*B(366)&
               &+0.506*B(368)+B(370)+B(372)+B(374)+B(376)+B(378)+0.38*B(383)+0.5*B(385)

  JVS(796) = -B(37)-B(38)-2*B(42)-2*B(43)-B(57)-B(59)-B(62)-B(207)-B(271)-B(273)-B(275)-B(277)-B(279)-B(281)-B(283)&
               &-B(285)-B(287)-B(289)-B(291)-B(293)-B(295)-B(297)-B(299)-B(301)-B(303)-B(305)-B(307)-B(309)-B(311)-B(431)

  JVS(797) = B(35)-B(39)+B(41)+B(55)+B(84)+B(87)+B(89)+0.381*B(99)+0.25*B(101)+0.049*B(103)+0.1*B(119)+0.1*B(121)+0.05&
               &*B(123)+B(125)+B(131)+B(133)+0.49*B(137)+0.5*B(139)+B(141)+0.35*B(147)+0.4*B(151)

  JVS(798) = B(21)

  JVS(799) = B(56)-B(63)+B(155)+B(159)+0.5*B(165)+B(393)+B(395)+0.792*B(397)+0.699*B(399)+0.845*B(401)+B(403)+B(405)&
               &+B(407)+B(409)+B(411)+B(413)+B(415)+B(417)+B(419)+0.77*B(425)+B(427)

  JVS(800) = B(14)

  JVS(801) = -B(58)+B(233)+B(235)+0.742*B(237)+0.599*B(239)+0.606*B(241)+B(243)+B(245)+B(247)+0.847*B(249)+0.8*B(251)&
               &+0.65*B(253)+0.95*B(255)+0.95*B(257)+B(259)+0.77*B(265)+B(267)

  JVS(802) = -B(60)

  JVS(803) = -B(86)

  JVS(804) = 2*B(32)

  JVS(805) = -B(100)

  JVS(806) = -B(118)

  JVS(807) = -B(120)

  JVS(808) = -B(102)

  JVS(809) = -0.964*B(98)

  JVS(810) = -B(96)

  JVS(811) = -B(94)

  JVS(812) = -B(140)

  JVS(813) = 0.35*B(6)-B(68)

  JVS(814) = B(13)-0.65*B(142)

  JVS(815) = B(4)-B(64)

  JVS(816) = 2*B(9)-B(40)

  JVS(817) = B(212)

  JVS(818) = B(218)

  JVS(819) = -B(130)

  JVS(820) = -B(104)+0.12*B(184)

  JVS(821) = B(224)

  JVS(822) = B(15)-B(146)

  JVS(823) = B(5)-B(66)

  JVS(824) = -B(88)

  JVS(825) = -B(114)+0.85*B(194)

  JVS(826) = -B(116)+0.85*B(196)

  JVS(827) = -B(148)

  JVS(828) = -B(122)

  JVS(829) = -B(110)+0.28*B(190)

  JVS(830) = -B(132)

  JVS(831) = -B(150)+0.036*B(202)

  JVS(832) = 0.02*B(90)-B(112)+0.28*B(192)

  JVS(833) = -B(134)

  JVS(834) = 0

  JVS(835) = -B(128)

  JVS(836) = 0

  JVS(837) = 0

  JVS(838) = 0

  JVS(839) = 0

  JVS(840) = -B(124)

  JVS(841) = 0

  JVS(842) = 0

  JVS(843) = -B(136)+0.07*B(198)

  JVS(844) = 0

  JVS(845) = -B(126)

  JVS(846) = -B(138)+0.21*B(200)

  JVS(847) = 0

  JVS(848) = 0

  JVS(849) = -B(106)+0.4*B(186)

  JVS(850) = -B(108)+0.63*B(188)

  JVS(851) = 0

  JVS(852) = 0

  JVS(853) = 0

  JVS(854) = 0.02*B(91)

  JVS(855) = -B(34)+B(36)+0.12*B(185)+0.4*B(187)+0.63*B(189)+0.28*B(191)+0.28*B(193)+0.85*B(195)+0.85*B(197)+0.07*B(199)&
               &+0.21*B(201)+0.036*B(203)+B(213)+B(219)+B(225)

  JVS(856) = 0

  JVS(857) = 0

  JVS(858) = 0

  JVS(859) = B(37)-B(38)+B(57)+0.7*B(62)

  JVS(860) = -B(35)-B(39)-B(41)-B(51)-B(53)-B(55)-B(65)-B(67)-B(69)-B(84)-B(87)-B(89)-B(95)-B(97)-0.964*B(99)-B(101)&
               &-B(103)-B(105)-B(107)-B(109)-B(111)-B(113)-B(115)-B(117)-B(119)-B(121)-B(123)-B(125)-B(127)-B(129)-B(131)&
               &-B(133)-B(135)-B(137)-B(139)-B(141)-0.65*B(143)-0.51*B(144)-B(147)-B(149)-B(151)-B(152)

  JVS(861) = -B(153)

  JVS(862) = -B(56)+0.7*B(63)

  JVS(863) = B(14)-0.51*B(145)

  JVS(864) = -B(52)+B(58)

  JVS(865) = -B(54)

  JVS(866) = B(204)

  JVS(867) = 0

  JVS(868) = 0.6*B(182)

  JVS(869) = 0.35*B(252)

  JVS(870) = 0.124*B(238)

  JVS(871) = 0.261*B(240)

  JVS(872) = 0.05*B(254)

  JVS(873) = 0.05*B(256)

  JVS(874) = 0.2*B(250)

  JVS(875) = 0.153*B(248)

  JVS(876) = 0.059*B(236)

  JVS(877) = 0

  JVS(878) = B(266)+B(308)+B(349)+B(384)+2*B(388)+1.5*B(389)+B(426)

  JVS(879) = B(310)+0.5*B(351)+0.484*B(386)+1.5*B(390)+B(391)

  JVS(880) = 0

  JVS(881) = 0

  JVS(882) = B(350)+0.5*B(352)

  JVS(883) = B(385)+0.484*B(387)

  JVS(884) = B(309)+B(311)

  JVS(885) = -B(152)

  JVS(886) = -B(21)-B(153)

  JVS(887) = 0.6*B(183)+B(427)

  JVS(888) = 0

  JVS(889) = 0.059*B(237)+0.124*B(239)+0.261*B(241)+0.153*B(249)+0.2*B(251)+0.35*B(253)+0.05*B(255)+0.05*B(257)+B(267)

  JVS(890) = B(205)

  JVS(891) = B(82)

  JVS(892) = 0.35*B(6)

  JVS(893) = -B(168)

  JVS(894) = B(66)

  JVS(895) = -B(178)

  JVS(896) = -B(180)

  JVS(897) = B(148)

  JVS(898) = -B(166)

  JVS(899) = -B(174)

  JVS(900) = -B(158)

  JVS(901) = 0.6*B(150)-0.4*B(182)

  JVS(902) = -B(402)

  JVS(903) = -B(176)

  JVS(904) = -B(404)

  JVS(905) = -B(406)

  JVS(906) = -B(160)

  JVS(907) = -B(418)

  JVS(908) = -B(412)

  JVS(909) = -B(398)

  JVS(910) = -B(400)

  JVS(911) = -B(414)

  JVS(912) = -B(416)

  JVS(913) = -B(154)

  JVS(914) = -B(410)

  JVS(915) = -B(408)

  JVS(916) = -B(162)

  JVS(917) = -B(396)

  JVS(918) = -B(156)

  JVS(919) = -B(164)

  JVS(920) = -B(422)

  JVS(921) = -B(439)

  JVS(922) = -B(170)

  JVS(923) = -B(172)

  JVS(924) = -B(426)

  JVS(925) = -B(428)

  JVS(926) = -B(394)

  JVS(927) = B(49)

  JVS(928) = B(72)

  JVS(929) = -B(424)

  JVS(930) = -B(392)

  JVS(931) = -B(420)

  JVS(932) = -B(62)

  JVS(933) = -B(55)+B(67)+B(149)+0.6*B(151)

  JVS(934) = 0

  JVS(935) = -B(7)-B(8)-B(56)-B(63)-B(76)-B(78)-B(80)-2*B(83)-B(155)-B(157)-B(159)-B(161)-B(163)-B(165)-B(167)-B(169)&
               &-B(171)-B(173)-B(175)-B(177)-B(179)-B(181)-0.4*B(183)-B(393)-B(395)-B(397)-B(399)-B(401)-B(403)-B(405)&
               &-B(407)-B(409)-B(411)-B(413)-B(415)-B(417)-B(419)-B(421)-B(423)-B(425)-B(427)-B(429)-B(440)

  JVS(936) = 0

  JVS(937) = -B(77)

  JVS(938) = B(50)+B(73)-B(79)-B(81)

  JVS(939) = B(280)

  JVS(940) = B(282)

  JVS(941) = B(284)

  JVS(942) = B(296)

  JVS(943) = B(290)

  JVS(944) = B(276)

  JVS(945) = B(278)

  JVS(946) = B(292)

  JVS(947) = B(294)

  JVS(948) = B(288)

  JVS(949) = B(286)

  JVS(950) = 0.13*B(198)

  JVS(951) = B(274)

  JVS(952) = B(302)

  JVS(953) = B(430)

  JVS(954) = 0

  JVS(955) = 0

  JVS(956) = B(272)

  JVS(957) = 0

  JVS(958) = 0.13*B(199)

  JVS(959) = B(306)

  JVS(960) = 0

  JVS(961) = 0

  JVS(962) = B(273)+B(275)+B(277)+B(279)+B(281)+B(283)+B(285)+B(287)+B(289)+B(291)+B(293)+B(295)+B(297)+B(303)+B(307)&
               &+B(431)

  JVS(963) = -B(144)

  JVS(964) = 0

  JVS(965) = 0

  JVS(966) = -B(14)-B(145)

  JVS(967) = 0

  JVS(968) = 0

  JVS(969) = B(4)

  JVS(970) = 0

  JVS(971) = 0

  JVS(972) = 0

  JVS(973) = 0

  JVS(974) = -B(242)

  JVS(975) = -B(244)

  JVS(976) = -B(246)

  JVS(977) = -B(258)

  JVS(978) = -B(252)

  JVS(979) = -B(238)

  JVS(980) = -B(240)

  JVS(981) = -B(254)

  JVS(982) = -B(256)

  JVS(983) = -B(250)

  JVS(984) = -B(248)

  JVS(985) = -B(236)

  JVS(986) = -B(262)

  JVS(987) = -B(437)

  JVS(988) = 0

  JVS(989) = 0

  JVS(990) = -B(266)

  JVS(991) = -B(268)

  JVS(992) = -B(234)

  JVS(993) = -B(45)+B(47)

  JVS(994) = -B(70)

  JVS(995) = -B(264)

  JVS(996) = -B(232)

  JVS(997) = -B(260)

  JVS(998) = -B(57)

  JVS(999) = -B(51)

  JVS(1000) = 0

  JVS(1001) = B(7)-B(76)+B(78)

  JVS(1002) = 0

  JVS(1003) = -B(46)-B(52)-B(58)-B(71)-2*B(74)-B(77)-B(233)-B(235)-B(237)-B(239)-B(241)-B(243)-B(245)-B(247)-B(249)&
                &-B(251)-B(253)-B(255)-B(257)-B(259)-B(261)-B(263)-B(265)-B(267)-B(269)-B(438)

  JVS(1004) = B(1)+B(48)+B(79)

  JVS(1005) = B(82)

  JVS(1006) = 0.65*B(6)+B(61)+B(68)

  JVS(1007) = B(64)

  JVS(1008) = -B(204)

  JVS(1009) = -B(208)

  JVS(1010) = -B(214)

  JVS(1011) = -B(220)

  JVS(1012) = B(5)

  JVS(1013) = B(228)

  JVS(1014) = 0

  JVS(1015) = 0

  JVS(1016) = 0.4*B(182)+0.7*B(202)+B(231)

  JVS(1017) = B(242)+B(402)

  JVS(1018) = B(244)+B(404)

  JVS(1019) = B(246)+B(406)

  JVS(1020) = 0

  JVS(1021) = B(258)+B(418)

  JVS(1022) = 0.65*B(252)+B(412)

  JVS(1023) = 0.876*B(238)+B(398)

  JVS(1024) = 0.739*B(240)+B(400)

  JVS(1025) = 0.95*B(254)+B(414)

  JVS(1026) = 0.95*B(256)+B(416)

  JVS(1027) = 0

  JVS(1028) = 0.8*B(250)+B(410)

  JVS(1029) = 0.847*B(248)+B(408)

  JVS(1030) = 0

  JVS(1031) = 0.941*B(236)+B(396)

  JVS(1032) = 0

  JVS(1033) = 0.5*B(164)

  JVS(1034) = -B(229)+B(262)+B(422)

  JVS(1035) = B(437)+B(439)

  JVS(1036) = 0

  JVS(1037) = 0

  JVS(1038) = B(266)+0.5*B(389)+B(426)

  JVS(1039) = 2*B(268)+0.5*B(351)+0.516*B(386)+0.5*B(390)+B(391)+2*B(428)

  JVS(1040) = B(234)+B(394)

  JVS(1041) = B(45)-B(47)-B(49)

  JVS(1042) = B(70)-B(72)+0.7*B(203)

  JVS(1043) = B(264)+B(424)

  JVS(1044) = B(232)+0.5*B(352)+B(392)

  JVS(1045) = -B(226)+B(260)+0.516*B(387)+B(420)

  JVS(1046) = B(57)-B(59)+0.7*B(62)

  JVS(1047) = -B(53)+B(55)+B(65)+B(69)+B(152)

  JVS(1048) = B(21)+B(153)

  JVS(1049) = B(8)+B(56)+0.7*B(63)+2*B(76)-B(80)+2*B(83)+0.5*B(165)+0.4*B(183)+B(393)+B(395)+B(397)+B(399)+B(401)+B(403)&
                &+B(405)+B(407)+B(409)+B(411)+B(413)+B(415)+B(417)+B(419)+B(421)+B(423)+B(425)+B(427)+2*B(429)+B(440)

  JVS(1050) = 0

  JVS(1051) = B(46)+B(58)+B(71)+2*B(74)+2*B(77)+B(233)+B(235)+0.941*B(237)+0.876*B(239)+0.739*B(241)+B(243)+B(245)&
                &+B(247)+0.847*B(249)+0.8*B(251)+0.65*B(253)+0.95*B(255)+0.95*B(257)+B(259)+B(261)+B(263)+B(265)+B(267)+2&
                &*B(269)+B(438)

  JVS(1052) = -B(1)-B(48)-B(50)-B(54)-B(60)-B(73)-B(81)-B(205)-B(209)-B(215)-B(221)-B(227)-B(230)
      
END SUBROUTINE racmpm_Jac_SP














SUBROUTINE racmpm_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(1052), W(73), a
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
      
END SUBROUTINE racmpm_KppDecomp



SUBROUTINE racmpm_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(1052), W(73), a
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
      
END SUBROUTINE racmpm_KppDecompCmplx


SUBROUTINE racmpm_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(1052), X(73), sum

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
      
END SUBROUTINE racmpm_KppSolveIndirect


SUBROUTINE racmpm_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(1052), X(73), sum

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
      
END SUBROUTINE racmpm_KppSolveCmplx













SUBROUTINE racmpm_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(21) = X(21)-JVS(102)*X(8)
  X(22) = X(22)-JVS(107)*X(9)
  X(23) = X(23)-JVS(112)*X(11)
  X(28) = X(28)-JVS(146)*X(12)-JVS(147)*X(24)
  X(32) = X(32)-JVS(181)*X(20)-JVS(182)*X(21)-JVS(183)*X(22)-JVS(184)*X(25)
  X(34) = X(34)-JVS(195)*X(12)
  X(36) = X(36)-JVS(214)*X(24)
  X(40) = X(40)-JVS(244)*X(23)-JVS(245)*X(35)
  X(41) = X(41)-JVS(262)*X(25)-JVS(263)*X(32)
  X(42) = X(42)-JVS(273)*X(7)-JVS(274)*X(15)-JVS(275)*X(29)-JVS(276)*X(38)-JVS(277)*X(39)
  X(43) = X(43)-JVS(298)*X(30)
  X(44) = X(44)-JVS(307)*X(7)
  X(45) = X(45)-JVS(315)*X(11)
  X(46) = X(46)-JVS(323)*X(21)
  X(47) = X(47)-JVS(333)*X(22)
  X(48) = X(48)-JVS(343)*X(12)-JVS(344)*X(17)-JVS(345)*X(23)-JVS(346)*X(24)-JVS(347)*X(26)-JVS(348)*X(30)-JVS(349)*X(31)&
            &-JVS(350)*X(33)-JVS(351)*X(34)-JVS(352)*X(35)-JVS(353)*X(36)-JVS(354)*X(37)-JVS(355)*X(38)-JVS(356)*X(39)&
            &-JVS(357)*X(41)-JVS(358)*X(43)-JVS(359)*X(44)-JVS(360)*X(45)-JVS(361)*X(46)-JVS(362)*X(47)
  X(49) = X(49)-JVS(386)*X(29)
  X(50) = X(50)-JVS(395)*X(33)-JVS(396)*X(37)
  X(51) = X(51)-JVS(406)*X(30)-JVS(407)*X(33)-JVS(408)*X(37)-JVS(409)*X(43)-JVS(410)*X(50)
  X(52) = X(52)-JVS(420)*X(12)
  X(53) = X(53)-JVS(430)*X(11)-JVS(431)*X(12)-JVS(432)*X(15)-JVS(433)*X(29)-JVS(434)*X(36)-JVS(435)*X(38)-JVS(436)*X(39)&
            &-JVS(437)*X(44)-JVS(438)*X(45)-JVS(439)*X(49)-JVS(440)*X(51)-JVS(441)*X(52)
  X(54) = X(54)-JVS(460)*X(37)-JVS(461)*X(46)-JVS(462)*X(47)
  X(55) = X(55)-JVS(473)*X(35)-JVS(474)*X(51)-JVS(475)*X(54)
  X(56) = X(56)-JVS(486)*X(8)-JVS(487)*X(9)-JVS(488)*X(26)-JVS(489)*X(31)-JVS(490)*X(32)-JVS(491)*X(33)-JVS(492)*X(35)&
            &-JVS(493)*X(37)-JVS(494)*X(44)-JVS(495)*X(45)-JVS(496)*X(51)-JVS(497)*X(52)-JVS(498)*X(54)-JVS(499)*X(55)
  X(57) = X(57)-JVS(513)*X(30)-JVS(514)*X(33)-JVS(515)*X(37)-JVS(516)*X(50)
  X(58) = X(58)-JVS(526)*X(43)-JVS(527)*X(50)
  X(59) = X(59)-JVS(537)*X(24)-JVS(538)*X(29)-JVS(539)*X(30)-JVS(540)*X(33)-JVS(541)*X(37)-JVS(542)*X(51)-JVS(543)*X(57)&
            &-JVS(544)*X(58)
  X(60) = X(60)-JVS(555)*X(24)-JVS(556)*X(29)-JVS(557)*X(30)-JVS(558)*X(33)-JVS(559)*X(37)-JVS(560)*X(57)-JVS(561)*X(58)&
            &-JVS(562)*X(59)
  X(61) = X(61)-JVS(572)*X(13)-JVS(573)*X(29)-JVS(574)*X(30)-JVS(575)*X(42)-JVS(576)*X(44)-JVS(577)*X(45)-JVS(578)*X(49)&
            &-JVS(579)*X(52)-JVS(580)*X(54)-JVS(581)*X(57)-JVS(582)*X(58)-JVS(583)*X(59)-JVS(584)*X(60)
  X(62) = X(62)-JVS(598)*X(6)-JVS(599)*X(33)-JVS(600)*X(37)-JVS(601)*X(51)
  X(63) = X(63)-JVS(611)*X(6)-JVS(612)*X(21)-JVS(613)*X(22)-JVS(614)*X(24)-JVS(615)*X(25)-JVS(616)*X(29)-JVS(617)*X(30)&
            &-JVS(618)*X(32)-JVS(619)*X(33)-JVS(620)*X(35)-JVS(621)*X(37)-JVS(622)*X(51)-JVS(623)*X(54)-JVS(624)*X(55)&
            &-JVS(625)*X(57)-JVS(626)*X(58)-JVS(627)*X(62)
  X(64) = X(64)-JVS(636)*X(29)-JVS(637)*X(30)-JVS(638)*X(33)-JVS(639)*X(37)-JVS(640)*X(42)-JVS(641)*X(44)-JVS(642)*X(45)&
            &-JVS(643)*X(49)-JVS(644)*X(52)-JVS(645)*X(54)-JVS(646)*X(57)-JVS(647)*X(58)-JVS(648)*X(59)-JVS(649)*X(60)&
            &-JVS(650)*X(62)-JVS(651)*X(63)
  X(65) = X(65)-JVS(662)*X(14)-JVS(663)*X(17)-JVS(664)*X(26)-JVS(665)*X(33)-JVS(666)*X(36)-JVS(667)*X(37)-JVS(668)*X(38)&
            &-JVS(669)*X(39)-JVS(670)*X(41)-JVS(671)*X(43)-JVS(672)*X(44)-JVS(673)*X(45)-JVS(674)*X(46)-JVS(675)*X(47)&
            &-JVS(676)*X(49)-JVS(677)*X(50)-JVS(678)*X(52)-JVS(679)*X(53)-JVS(680)*X(54)-JVS(681)*X(55)-JVS(682)*X(56)&
            &-JVS(683)*X(57)-JVS(684)*X(58)-JVS(685)*X(59)-JVS(686)*X(60)-JVS(687)*X(61)-JVS(688)*X(62)-JVS(689)*X(63)&
            &-JVS(690)*X(64)
  X(66) = X(66)-JVS(700)*X(23)-JVS(701)*X(26)-JVS(702)*X(31)-JVS(703)*X(33)-JVS(704)*X(35)-JVS(705)*X(36)-JVS(706)*X(37)&
            &-JVS(707)*X(38)-JVS(708)*X(39)-JVS(709)*X(40)-JVS(710)*X(41)-JVS(711)*X(42)-JVS(712)*X(43)-JVS(713)*X(44)&
            &-JVS(714)*X(45)-JVS(715)*X(46)-JVS(716)*X(47)-JVS(717)*X(49)-JVS(718)*X(50)-JVS(719)*X(51)-JVS(720)*X(52)&
            &-JVS(721)*X(53)-JVS(722)*X(54)-JVS(723)*X(55)-JVS(724)*X(56)-JVS(725)*X(57)-JVS(726)*X(58)-JVS(727)*X(59)&
            &-JVS(728)*X(60)-JVS(729)*X(61)-JVS(730)*X(62)-JVS(731)*X(63)-JVS(732)*X(64)-JVS(733)*X(65)
  X(67) = X(67)-JVS(742)*X(5)-JVS(743)*X(7)-JVS(744)*X(8)-JVS(745)*X(9)-JVS(746)*X(11)-JVS(747)*X(12)-JVS(748)*X(15)&
            &-JVS(749)*X(16)-JVS(750)*X(17)-JVS(751)*X(19)-JVS(752)*X(20)-JVS(753)*X(21)-JVS(754)*X(22)-JVS(755)*X(23)&
            &-JVS(756)*X(24)-JVS(757)*X(25)-JVS(758)*X(26)-JVS(759)*X(28)-JVS(760)*X(29)-JVS(761)*X(30)-JVS(762)*X(32)&
            &-JVS(763)*X(33)-JVS(764)*X(34)-JVS(765)*X(35)-JVS(766)*X(36)-JVS(767)*X(37)-JVS(768)*X(38)-JVS(769)*X(39)&
            &-JVS(770)*X(40)-JVS(771)*X(41)-JVS(772)*X(43)-JVS(773)*X(44)-JVS(774)*X(45)-JVS(775)*X(46)-JVS(776)*X(47)&
            &-JVS(777)*X(48)-JVS(778)*X(49)-JVS(779)*X(50)-JVS(780)*X(51)-JVS(781)*X(52)-JVS(782)*X(53)-JVS(783)*X(54)&
            &-JVS(784)*X(55)-JVS(785)*X(56)-JVS(786)*X(57)-JVS(787)*X(58)-JVS(788)*X(59)-JVS(789)*X(60)-JVS(790)*X(61)&
            &-JVS(791)*X(62)-JVS(792)*X(63)-JVS(793)*X(64)-JVS(794)*X(65)-JVS(795)*X(66)
  X(68) = X(68)-JVS(803)*X(5)-JVS(804)*X(6)-JVS(805)*X(7)-JVS(806)*X(8)-JVS(807)*X(9)-JVS(808)*X(11)-JVS(809)*X(12)&
            &-JVS(810)*X(13)-JVS(811)*X(14)-JVS(812)*X(15)-JVS(813)*X(16)-JVS(814)*X(17)-JVS(815)*X(18)-JVS(816)*X(19)&
            &-JVS(817)*X(21)-JVS(818)*X(22)-JVS(819)*X(23)-JVS(820)*X(24)-JVS(821)*X(25)-JVS(822)*X(26)-JVS(823)*X(27)&
            &-JVS(824)*X(28)-JVS(825)*X(29)-JVS(826)*X(30)-JVS(827)*X(31)-JVS(828)*X(32)-JVS(829)*X(33)-JVS(830)*X(34)&
            &-JVS(831)*X(35)-JVS(832)*X(37)-JVS(833)*X(40)-JVS(834)*X(41)-JVS(835)*X(42)-JVS(836)*X(44)-JVS(837)*X(45)&
            &-JVS(838)*X(46)-JVS(839)*X(47)-JVS(840)*X(48)-JVS(841)*X(49)-JVS(842)*X(50)-JVS(843)*X(51)-JVS(844)*X(52)&
            &-JVS(845)*X(53)-JVS(846)*X(54)-JVS(847)*X(55)-JVS(848)*X(56)-JVS(849)*X(57)-JVS(850)*X(58)-JVS(851)*X(59)&
            &-JVS(852)*X(60)-JVS(853)*X(61)-JVS(854)*X(62)-JVS(855)*X(63)-JVS(856)*X(64)-JVS(857)*X(65)-JVS(858)*X(66)&
            &-JVS(859)*X(67)
  X(69) = X(69)-JVS(866)*X(20)-JVS(867)*X(32)-JVS(868)*X(35)-JVS(869)*X(43)-JVS(870)*X(44)-JVS(871)*X(45)-JVS(872)*X(46)&
            &-JVS(873)*X(47)-JVS(874)*X(49)-JVS(875)*X(50)-JVS(876)*X(52)-JVS(877)*X(55)-JVS(878)*X(59)-JVS(879)*X(60)&
            &-JVS(880)*X(62)-JVS(881)*X(63)-JVS(882)*X(65)-JVS(883)*X(66)-JVS(884)*X(67)-JVS(885)*X(68)
  X(70) = X(70)-JVS(891)*X(10)-JVS(892)*X(16)-JVS(893)*X(24)-JVS(894)*X(27)-JVS(895)*X(29)-JVS(896)*X(30)-JVS(897)*X(31)&
            &-JVS(898)*X(32)-JVS(899)*X(33)-JVS(900)*X(34)-JVS(901)*X(35)-JVS(902)*X(36)-JVS(903)*X(37)-JVS(904)*X(38)&
            &-JVS(905)*X(39)-JVS(906)*X(40)-JVS(907)*X(41)-JVS(908)*X(43)-JVS(909)*X(44)-JVS(910)*X(45)-JVS(911)*X(46)&
            &-JVS(912)*X(47)-JVS(913)*X(48)-JVS(914)*X(49)-JVS(915)*X(50)-JVS(916)*X(51)-JVS(917)*X(52)-JVS(918)*X(53)&
            &-JVS(919)*X(54)-JVS(920)*X(55)-JVS(921)*X(56)-JVS(922)*X(57)-JVS(923)*X(58)-JVS(924)*X(59)-JVS(925)*X(60)&
            &-JVS(926)*X(61)-JVS(927)*X(62)-JVS(928)*X(63)-JVS(929)*X(64)-JVS(930)*X(65)-JVS(931)*X(66)-JVS(932)*X(67)&
            &-JVS(933)*X(68)-JVS(934)*X(69)
  X(71) = X(71)-JVS(939)*X(36)-JVS(940)*X(38)-JVS(941)*X(39)-JVS(942)*X(41)-JVS(943)*X(43)-JVS(944)*X(44)-JVS(945)*X(45)&
            &-JVS(946)*X(46)-JVS(947)*X(47)-JVS(948)*X(49)-JVS(949)*X(50)-JVS(950)*X(51)-JVS(951)*X(52)-JVS(952)*X(55)&
            &-JVS(953)*X(56)-JVS(954)*X(57)-JVS(955)*X(58)-JVS(956)*X(61)-JVS(957)*X(62)-JVS(958)*X(63)-JVS(959)*X(64)&
            &-JVS(960)*X(65)-JVS(961)*X(66)-JVS(962)*X(67)-JVS(963)*X(68)-JVS(964)*X(69)-JVS(965)*X(70)
  X(72) = X(72)-JVS(969)*X(18)-JVS(970)*X(21)-JVS(971)*X(22)-JVS(972)*X(25)-JVS(973)*X(32)-JVS(974)*X(36)-JVS(975)*X(38)&
            &-JVS(976)*X(39)-JVS(977)*X(41)-JVS(978)*X(43)-JVS(979)*X(44)-JVS(980)*X(45)-JVS(981)*X(46)-JVS(982)*X(47)&
            &-JVS(983)*X(49)-JVS(984)*X(50)-JVS(985)*X(52)-JVS(986)*X(55)-JVS(987)*X(56)-JVS(988)*X(57)-JVS(989)*X(58)&
            &-JVS(990)*X(59)-JVS(991)*X(60)-JVS(992)*X(61)-JVS(993)*X(62)-JVS(994)*X(63)-JVS(995)*X(64)-JVS(996)*X(65)&
            &-JVS(997)*X(66)-JVS(998)*X(67)-JVS(999)*X(68)-JVS(1000)*X(69)-JVS(1001)*X(70)-JVS(1002)*X(71)
  X(73) = X(73)-JVS(1005)*X(10)-JVS(1006)*X(16)-JVS(1007)*X(18)-JVS(1008)*X(20)-JVS(1009)*X(21)-JVS(1010)*X(22)&
            &-JVS(1011)*X(25)-JVS(1012)*X(27)-JVS(1013)*X(31)-JVS(1014)*X(32)-JVS(1015)*X(34)-JVS(1016)*X(35)-JVS(1017)&
            &*X(36)-JVS(1018)*X(38)-JVS(1019)*X(39)-JVS(1020)*X(40)-JVS(1021)*X(41)-JVS(1022)*X(43)-JVS(1023)*X(44)&
            &-JVS(1024)*X(45)-JVS(1025)*X(46)-JVS(1026)*X(47)-JVS(1027)*X(48)-JVS(1028)*X(49)-JVS(1029)*X(50)-JVS(1030)&
            &*X(51)-JVS(1031)*X(52)-JVS(1032)*X(53)-JVS(1033)*X(54)-JVS(1034)*X(55)-JVS(1035)*X(56)-JVS(1036)*X(57)&
            &-JVS(1037)*X(58)-JVS(1038)*X(59)-JVS(1039)*X(60)-JVS(1040)*X(61)-JVS(1041)*X(62)-JVS(1042)*X(63)-JVS(1043)&
            &*X(64)-JVS(1044)*X(65)-JVS(1045)*X(66)-JVS(1046)*X(67)-JVS(1047)*X(68)-JVS(1048)*X(69)-JVS(1049)*X(70)&
            &-JVS(1050)*X(71)-JVS(1051)*X(72)
  X(73) = X(73)/JVS(1052)
  X(72) = (X(72)-JVS(1004)*X(73))/(JVS(1003))
  X(71) = (X(71)-JVS(967)*X(72)-JVS(968)*X(73))/(JVS(966))
  X(70) = (X(70)-JVS(936)*X(71)-JVS(937)*X(72)-JVS(938)*X(73))/(JVS(935))
  X(69) = (X(69)-JVS(887)*X(70)-JVS(888)*X(71)-JVS(889)*X(72)-JVS(890)*X(73))/(JVS(886))
  X(68) = (X(68)-JVS(861)*X(69)-JVS(862)*X(70)-JVS(863)*X(71)-JVS(864)*X(72)-JVS(865)*X(73))/(JVS(860))
  X(67) = (X(67)-JVS(797)*X(68)-JVS(798)*X(69)-JVS(799)*X(70)-JVS(800)*X(71)-JVS(801)*X(72)-JVS(802)*X(73))/(JVS(796))
  X(66) = (X(66)-JVS(735)*X(67)-JVS(736)*X(68)-JVS(737)*X(69)-JVS(738)*X(70)-JVS(739)*X(71)-JVS(740)*X(72)-JVS(741)&
            &*X(73))/(JVS(734))
  X(65) = (X(65)-JVS(692)*X(66)-JVS(693)*X(67)-JVS(694)*X(68)-JVS(695)*X(69)-JVS(696)*X(70)-JVS(697)*X(71)-JVS(698)&
            &*X(72)-JVS(699)*X(73))/(JVS(691))
  X(64) = (X(64)-JVS(653)*X(65)-JVS(654)*X(66)-JVS(655)*X(67)-JVS(656)*X(68)-JVS(657)*X(69)-JVS(658)*X(70)-JVS(659)&
            &*X(71)-JVS(660)*X(72)-JVS(661)*X(73))/(JVS(652))
  X(63) = (X(63)-JVS(629)*X(65)-JVS(630)*X(66)-JVS(631)*X(67)-JVS(632)*X(68)-JVS(633)*X(70)-JVS(634)*X(72)-JVS(635)&
            &*X(73))/(JVS(628))
  X(62) = (X(62)-JVS(603)*X(63)-JVS(604)*X(65)-JVS(605)*X(66)-JVS(606)*X(67)-JVS(607)*X(68)-JVS(608)*X(70)-JVS(609)&
            &*X(72)-JVS(610)*X(73))/(JVS(602))
  X(61) = (X(61)-JVS(586)*X(62)-JVS(587)*X(63)-JVS(588)*X(64)-JVS(589)*X(65)-JVS(590)*X(66)-JVS(591)*X(67)-JVS(592)&
            &*X(68)-JVS(593)*X(69)-JVS(594)*X(70)-JVS(595)*X(71)-JVS(596)*X(72)-JVS(597)*X(73))/(JVS(585))
  X(60) = (X(60)-JVS(564)*X(62)-JVS(565)*X(63)-JVS(566)*X(65)-JVS(567)*X(66)-JVS(568)*X(67)-JVS(569)*X(68)-JVS(570)&
            &*X(70)-JVS(571)*X(72))/(JVS(563))
  X(59) = (X(59)-JVS(546)*X(60)-JVS(547)*X(62)-JVS(548)*X(63)-JVS(549)*X(65)-JVS(550)*X(66)-JVS(551)*X(67)-JVS(552)&
            &*X(68)-JVS(553)*X(70)-JVS(554)*X(72))/(JVS(545))
  X(58) = (X(58)-JVS(529)*X(62)-JVS(530)*X(63)-JVS(531)*X(65)-JVS(532)*X(66)-JVS(533)*X(67)-JVS(534)*X(68)-JVS(535)&
            &*X(70)-JVS(536)*X(72))/(JVS(528))
  X(57) = (X(57)-JVS(518)*X(62)-JVS(519)*X(63)-JVS(520)*X(65)-JVS(521)*X(66)-JVS(522)*X(67)-JVS(523)*X(68)-JVS(524)&
            &*X(70)-JVS(525)*X(72))/(JVS(517))
  X(56) = (X(56)-JVS(501)*X(62)-JVS(502)*X(63)-JVS(503)*X(64)-JVS(504)*X(65)-JVS(505)*X(66)-JVS(506)*X(67)-JVS(507)&
            &*X(68)-JVS(508)*X(69)-JVS(509)*X(70)-JVS(510)*X(71)-JVS(511)*X(72)-JVS(512)*X(73))/(JVS(500))
  X(55) = (X(55)-JVS(477)*X(62)-JVS(478)*X(63)-JVS(479)*X(65)-JVS(480)*X(66)-JVS(481)*X(67)-JVS(482)*X(68)-JVS(483)&
            &*X(70)-JVS(484)*X(72)-JVS(485)*X(73))/(JVS(476))
  X(54) = (X(54)-JVS(464)*X(62)-JVS(465)*X(63)-JVS(466)*X(65)-JVS(467)*X(66)-JVS(468)*X(67)-JVS(469)*X(68)-JVS(470)&
            &*X(70)-JVS(471)*X(72)-JVS(472)*X(73))/(JVS(463))
  X(53) = (X(53)-JVS(443)*X(54)-JVS(444)*X(57)-JVS(445)*X(58)-JVS(446)*X(59)-JVS(447)*X(60)-JVS(448)*X(61)-JVS(449)&
            &*X(62)-JVS(450)*X(63)-JVS(451)*X(64)-JVS(452)*X(65)-JVS(453)*X(66)-JVS(454)*X(67)-JVS(455)*X(68)-JVS(456)*X(69)&
            &-JVS(457)*X(70)-JVS(458)*X(71)-JVS(459)*X(72))/(JVS(442))
  X(52) = (X(52)-JVS(422)*X(65)-JVS(423)*X(66)-JVS(424)*X(67)-JVS(425)*X(68)-JVS(426)*X(69)-JVS(427)*X(70)-JVS(428)&
            &*X(71)-JVS(429)*X(72))/(JVS(421))
  X(51) = (X(51)-JVS(412)*X(62)-JVS(413)*X(63)-JVS(414)*X(65)-JVS(415)*X(66)-JVS(416)*X(67)-JVS(417)*X(68)-JVS(418)&
            &*X(70)-JVS(419)*X(72))/(JVS(411))
  X(50) = (X(50)-JVS(398)*X(62)-JVS(399)*X(63)-JVS(400)*X(65)-JVS(401)*X(66)-JVS(402)*X(67)-JVS(403)*X(68)-JVS(404)&
            &*X(70)-JVS(405)*X(72))/(JVS(397))
  X(49) = (X(49)-JVS(388)*X(63)-JVS(389)*X(65)-JVS(390)*X(66)-JVS(391)*X(67)-JVS(392)*X(68)-JVS(393)*X(70)-JVS(394)&
            &*X(72))/(JVS(387))
  X(48) = (X(48)-JVS(364)*X(49)-JVS(365)*X(50)-JVS(366)*X(51)-JVS(367)*X(52)-JVS(368)*X(54)-JVS(369)*X(55)-JVS(370)&
            &*X(56)-JVS(371)*X(57)-JVS(372)*X(58)-JVS(373)*X(59)-JVS(374)*X(60)-JVS(375)*X(61)-JVS(376)*X(62)-JVS(377)*X(63)&
            &-JVS(378)*X(64)-JVS(379)*X(65)-JVS(380)*X(66)-JVS(381)*X(67)-JVS(382)*X(68)-JVS(383)*X(70)-JVS(384)*X(72)&
            &-JVS(385)*X(73))/(JVS(363))
  X(47) = (X(47)-JVS(335)*X(63)-JVS(336)*X(65)-JVS(337)*X(66)-JVS(338)*X(67)-JVS(339)*X(68)-JVS(340)*X(70)-JVS(341)&
            &*X(72)-JVS(342)*X(73))/(JVS(334))
  X(46) = (X(46)-JVS(325)*X(63)-JVS(326)*X(65)-JVS(327)*X(66)-JVS(328)*X(67)-JVS(329)*X(68)-JVS(330)*X(70)-JVS(331)&
            &*X(72)-JVS(332)*X(73))/(JVS(324))
  X(45) = (X(45)-JVS(317)*X(65)-JVS(318)*X(66)-JVS(319)*X(67)-JVS(320)*X(68)-JVS(321)*X(70)-JVS(322)*X(72))/(JVS(316))
  X(44) = (X(44)-JVS(309)*X(65)-JVS(310)*X(66)-JVS(311)*X(67)-JVS(312)*X(68)-JVS(313)*X(70)-JVS(314)*X(72))/(JVS(308))
  X(43) = (X(43)-JVS(300)*X(63)-JVS(301)*X(65)-JVS(302)*X(66)-JVS(303)*X(67)-JVS(304)*X(68)-JVS(305)*X(70)-JVS(306)&
            &*X(72))/(JVS(299))
  X(42) = (X(42)-JVS(279)*X(44)-JVS(280)*X(45)-JVS(281)*X(49)-JVS(282)*X(52)-JVS(283)*X(54)-JVS(284)*X(57)-JVS(285)&
            &*X(58)-JVS(286)*X(59)-JVS(287)*X(60)-JVS(288)*X(63)-JVS(289)*X(64)-JVS(290)*X(65)-JVS(291)*X(66)-JVS(292)*X(67)&
            &-JVS(293)*X(68)-JVS(294)*X(69)-JVS(295)*X(70)-JVS(296)*X(71)-JVS(297)*X(72))/(JVS(278))
  X(41) = (X(41)-JVS(265)*X(63)-JVS(266)*X(65)-JVS(267)*X(66)-JVS(268)*X(67)-JVS(269)*X(68)-JVS(270)*X(70)-JVS(271)&
            &*X(72)-JVS(272)*X(73))/(JVS(264))
  X(40) = (X(40)-JVS(247)*X(41)-JVS(248)*X(46)-JVS(249)*X(47)-JVS(250)*X(51)-JVS(251)*X(52)-JVS(252)*X(54)-JVS(253)&
            &*X(55)-JVS(254)*X(63)-JVS(255)*X(64)-JVS(256)*X(65)-JVS(257)*X(66)-JVS(258)*X(68)-JVS(259)*X(70)-JVS(260)*X(72)&
            &-JVS(261)*X(73))/(JVS(246))
  X(39) = (X(39)-JVS(237)*X(58)-JVS(238)*X(65)-JVS(239)*X(66)-JVS(240)*X(67)-JVS(241)*X(68)-JVS(242)*X(70)-JVS(243)&
            &*X(72))/(JVS(236))
  X(38) = (X(38)-JVS(229)*X(57)-JVS(230)*X(65)-JVS(231)*X(66)-JVS(232)*X(67)-JVS(233)*X(68)-JVS(234)*X(70)-JVS(235)&
            &*X(72))/(JVS(228))
  X(37) = (X(37)-JVS(224)*X(62)-JVS(225)*X(63)-JVS(226)*X(68)-JVS(227)*X(70))/(JVS(223))
  X(36) = (X(36)-JVS(216)*X(63)-JVS(217)*X(65)-JVS(218)*X(66)-JVS(219)*X(67)-JVS(220)*X(68)-JVS(221)*X(70)-JVS(222)&
            &*X(72))/(JVS(215))
  X(35) = (X(35)-JVS(209)*X(55)-JVS(210)*X(63)-JVS(211)*X(68)-JVS(212)*X(70)-JVS(213)*X(73))/(JVS(208))
  X(34) = (X(34)-JVS(197)*X(41)-JVS(198)*X(46)-JVS(199)*X(47)-JVS(200)*X(52)-JVS(201)*X(54)-JVS(202)*X(63)-JVS(203)&
            &*X(65)-JVS(204)*X(66)-JVS(205)*X(68)-JVS(206)*X(70)-JVS(207)*X(72))/(JVS(196))
  X(33) = (X(33)-JVS(192)*X(63)-JVS(193)*X(68)-JVS(194)*X(70))/(JVS(191))
  X(32) = (X(32)-JVS(186)*X(63)-JVS(187)*X(67)-JVS(188)*X(68)-JVS(189)*X(70)-JVS(190)*X(73))/(JVS(185))
  X(31) = (X(31)-JVS(175)*X(35)-JVS(176)*X(63)-JVS(177)*X(66)-JVS(178)*X(68)-JVS(179)*X(70)-JVS(180)*X(73))/(JVS(174))
  X(30) = (X(30)-JVS(171)*X(63)-JVS(172)*X(68)-JVS(173)*X(70))/(JVS(170))
  X(29) = (X(29)-JVS(167)*X(63)-JVS(168)*X(68)-JVS(169)*X(70))/(JVS(166))
  X(28) = (X(28)-JVS(149)*X(29)-JVS(150)*X(30)-JVS(151)*X(33)-JVS(152)*X(34)-JVS(153)*X(35)-JVS(154)*X(37)-JVS(155)&
            &*X(40)-JVS(156)*X(48)-JVS(157)*X(51)-JVS(158)*X(53)-JVS(159)*X(54)-JVS(160)*X(57)-JVS(161)*X(58)-JVS(162)*X(62)&
            &-JVS(163)*X(63)-JVS(164)*X(68)-JVS(165)*X(70))/(JVS(148))
  X(27) = (X(27)-JVS(135)*X(32)-JVS(136)*X(34)-JVS(137)*X(40)-JVS(138)*X(48)-JVS(139)*X(51)-JVS(140)*X(53)-JVS(141)&
            &*X(54)-JVS(142)*X(67)-JVS(143)*X(68)-JVS(144)*X(70)-JVS(145)*X(73))/(JVS(134))
  X(26) = (X(26)-JVS(129)*X(54)-JVS(130)*X(63)-JVS(131)*X(66)-JVS(132)*X(67)-JVS(133)*X(68))/(JVS(128))
  X(25) = (X(25)-JVS(124)*X(32)-JVS(125)*X(63)-JVS(126)*X(68)-JVS(127)*X(73))/(JVS(123))
  X(24) = (X(24)-JVS(120)*X(63)-JVS(121)*X(68)-JVS(122)*X(70))/(JVS(119))
  X(23) = (X(23)-JVS(114)*X(35)-JVS(115)*X(51)-JVS(116)*X(64)-JVS(117)*X(65)-JVS(118)*X(68))/(JVS(113))
  X(22) = (X(22)-JVS(109)*X(63)-JVS(110)*X(68)-JVS(111)*X(73))/(JVS(108))
  X(21) = (X(21)-JVS(104)*X(63)-JVS(105)*X(68)-JVS(106)*X(73))/(JVS(103))
  X(20) = (X(20)-JVS(97)*X(32)-JVS(98)*X(67)-JVS(99)*X(68)-JVS(100)*X(70)-JVS(101)*X(73))/(JVS(96))
  X(19) = (X(19)-JVS(87)*X(29)-JVS(88)*X(30)-JVS(89)*X(33)-JVS(90)*X(37)-JVS(91)*X(57)-JVS(92)*X(58)-JVS(93)*X(63)&
            &-JVS(94)*X(67)-JVS(95)*X(68))/(JVS(86))
  X(18) = (X(18)-JVS(80)*X(21)-JVS(81)*X(22)-JVS(82)*X(25)-JVS(83)*X(68)-JVS(84)*X(72)-JVS(85)*X(73))/(JVS(79))
  X(17) = (X(17)-JVS(76)*X(65)-JVS(77)*X(67)-JVS(78)*X(68))/(JVS(75))
  X(16) = (X(16)-JVS(72)*X(67)-JVS(73)*X(68)-JVS(74)*X(73))/(JVS(71))
  X(15) = (X(15)-JVS(69)*X(54)-JVS(70)*X(68))/(JVS(68))
  X(14) = (X(14)-JVS(64)*X(57)-JVS(65)*X(58)-JVS(66)*X(63)-JVS(67)*X(68))/(JVS(63))
  X(13) = (X(13)-JVS(59)*X(57)-JVS(60)*X(58)-JVS(61)*X(63)-JVS(62)*X(68))/(JVS(58))
  X(12) = (X(12)-JVS(57)*X(68))/(JVS(56))
  X(11) = (X(11)-JVS(55)*X(68))/(JVS(54))
  X(10) = (X(10)-JVS(52)*X(70)-JVS(53)*X(73))/(JVS(51))
  X(9) = (X(9)-JVS(50)*X(68))/(JVS(49))
  X(8) = (X(8)-JVS(48)*X(68))/(JVS(47))
  X(7) = (X(7)-JVS(46)*X(68))/(JVS(45))
  X(6) = (X(6)-JVS(44)*X(63))/(JVS(43))
  X(5) = (X(5)-JVS(42)*X(68))/(JVS(41))
  X(4) = (X(4)-JVS(20)*X(30)-JVS(21)*X(36)-JVS(22)*X(38)-JVS(23)*X(39)-JVS(24)*X(44)-JVS(25)*X(45)-JVS(26)*X(50)-JVS(27)&
           &*X(51)-JVS(28)*X(52)-JVS(29)*X(54)-JVS(30)*X(55)-JVS(31)*X(57)-JVS(32)*X(58)-JVS(33)*X(59)-JVS(34)*X(60)-JVS(35)&
           &*X(61)-JVS(36)*X(63)-JVS(37)*X(64)-JVS(38)*X(65)-JVS(39)*X(66)-JVS(40)*X(67))/(JVS(19))
  X(3) = (X(3)-JVS(8)*X(12)-JVS(9)*X(24)-JVS(10)*X(30)-JVS(11)*X(33)-JVS(12)*X(35)-JVS(13)*X(37)-JVS(14)*X(51)-JVS(15)&
           &*X(54)-JVS(16)*X(57)-JVS(17)*X(63)-JVS(18)*X(68))/(JVS(7))
  X(2) = (X(2)-JVS(5)*X(28)-JVS(6)*X(68))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(5)-JVS(3)*X(68))/(JVS(1))
      
END SUBROUTINE racmpm_KppSolve
























      SUBROUTINE racmpm_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE racmpm_WCOPY



      SUBROUTINE racmpm_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE racmpm_WAXPY




      SUBROUTINE racmpm_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE racmpm_WSCAL


      REAL(kind=dp) FUNCTION racmpm_WLAMCH( C )








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
          CALL racmpm_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      racmpm_WLAMCH = Eps

      END FUNCTION racmpm_WLAMCH
     
      SUBROUTINE racmpm_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE racmpm_WLAMCH_ADD




      SUBROUTINE racmpm_SET2ZERO(N,Y)




      
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

      END SUBROUTINE racmpm_SET2ZERO



      REAL(kind=dp) FUNCTION racmpm_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      racmpm_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        racmpm_WDOT = racmpm_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         racmpm_WDOT = racmpm_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          racmpm_WDOT = racmpm_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        racmpm_WDOT = racmpm_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION racmpm_WDOT                                          




   SUBROUTINE decomp_racmpm( JVS, IER )
   
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
   W( 5 ) = JVS( 2 )
   W( 68 ) = JVS( 3 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 5 )
  JVS( 3) = W( 68 )
  IF ( ABS(  JVS( 4 )) < TINY(a) ) THEN
         IER = 2                                       
         RETURN
  END IF
   W( 2 ) = JVS( 4 )
   W( 28 ) = JVS( 5 )
   W( 68 ) = JVS( 6 )
  JVS( 4) = W( 2 )
  JVS( 5) = W( 28 )
  JVS( 6) = W( 68 )
  IF ( ABS(  JVS( 7 )) < TINY(a) ) THEN
         IER = 3                                       
         RETURN
  END IF
   W( 3 ) = JVS( 7 )
   W( 12 ) = JVS( 8 )
   W( 24 ) = JVS( 9 )
   W( 30 ) = JVS( 10 )
   W( 33 ) = JVS( 11 )
   W( 35 ) = JVS( 12 )
   W( 37 ) = JVS( 13 )
   W( 51 ) = JVS( 14 )
   W( 54 ) = JVS( 15 )
   W( 57 ) = JVS( 16 )
   W( 63 ) = JVS( 17 )
   W( 68 ) = JVS( 18 )
  JVS( 7) = W( 3 )
  JVS( 8) = W( 12 )
  JVS( 9) = W( 24 )
  JVS( 10) = W( 30 )
  JVS( 11) = W( 33 )
  JVS( 12) = W( 35 )
  JVS( 13) = W( 37 )
  JVS( 14) = W( 51 )
  JVS( 15) = W( 54 )
  JVS( 16) = W( 57 )
  JVS( 17) = W( 63 )
  JVS( 18) = W( 68 )
  IF ( ABS(  JVS( 19 )) < TINY(a) ) THEN
         IER = 4                                       
         RETURN
  END IF
   W( 4 ) = JVS( 19 )
   W( 30 ) = JVS( 20 )
   W( 36 ) = JVS( 21 )
   W( 38 ) = JVS( 22 )
   W( 39 ) = JVS( 23 )
   W( 44 ) = JVS( 24 )
   W( 45 ) = JVS( 25 )
   W( 50 ) = JVS( 26 )
   W( 51 ) = JVS( 27 )
   W( 52 ) = JVS( 28 )
   W( 54 ) = JVS( 29 )
   W( 55 ) = JVS( 30 )
   W( 57 ) = JVS( 31 )
   W( 58 ) = JVS( 32 )
   W( 59 ) = JVS( 33 )
   W( 60 ) = JVS( 34 )
   W( 61 ) = JVS( 35 )
   W( 63 ) = JVS( 36 )
   W( 64 ) = JVS( 37 )
   W( 65 ) = JVS( 38 )
   W( 66 ) = JVS( 39 )
   W( 67 ) = JVS( 40 )
  JVS( 19) = W( 4 )
  JVS( 20) = W( 30 )
  JVS( 21) = W( 36 )
  JVS( 22) = W( 38 )
  JVS( 23) = W( 39 )
  JVS( 24) = W( 44 )
  JVS( 25) = W( 45 )
  JVS( 26) = W( 50 )
  JVS( 27) = W( 51 )
  JVS( 28) = W( 52 )
  JVS( 29) = W( 54 )
  JVS( 30) = W( 55 )
  JVS( 31) = W( 57 )
  JVS( 32) = W( 58 )
  JVS( 33) = W( 59 )
  JVS( 34) = W( 60 )
  JVS( 35) = W( 61 )
  JVS( 36) = W( 63 )
  JVS( 37) = W( 64 )
  JVS( 38) = W( 65 )
  JVS( 39) = W( 66 )
  JVS( 40) = W( 67 )
  IF ( ABS(  JVS( 41 )) < TINY(a) ) THEN
         IER = 5                                       
         RETURN
  END IF
   W( 5 ) = JVS( 41 )
   W( 68 ) = JVS( 42 )
  JVS( 41) = W( 5 )
  JVS( 42) = W( 68 )
  IF ( ABS(  JVS( 43 )) < TINY(a) ) THEN
         IER = 6                                       
         RETURN
  END IF
   W( 6 ) = JVS( 43 )
   W( 63 ) = JVS( 44 )
  JVS( 43) = W( 6 )
  JVS( 44) = W( 63 )
  IF ( ABS(  JVS( 45 )) < TINY(a) ) THEN
         IER = 7                                       
         RETURN
  END IF
   W( 7 ) = JVS( 45 )
   W( 68 ) = JVS( 46 )
  JVS( 45) = W( 7 )
  JVS( 46) = W( 68 )
  IF ( ABS(  JVS( 47 )) < TINY(a) ) THEN
         IER = 8                                       
         RETURN
  END IF
   W( 8 ) = JVS( 47 )
   W( 68 ) = JVS( 48 )
  JVS( 47) = W( 8 )
  JVS( 48) = W( 68 )
  IF ( ABS(  JVS( 49 )) < TINY(a) ) THEN
         IER = 9                                       
         RETURN
  END IF
   W( 9 ) = JVS( 49 )
   W( 68 ) = JVS( 50 )
  JVS( 49) = W( 9 )
  JVS( 50) = W( 68 )
  IF ( ABS(  JVS( 51 )) < TINY(a) ) THEN
         IER = 10                                      
         RETURN
  END IF
   W( 10 ) = JVS( 51 )
   W( 70 ) = JVS( 52 )
   W( 73 ) = JVS( 53 )
  JVS( 51) = W( 10 )
  JVS( 52) = W( 70 )
  JVS( 53) = W( 73 )
  IF ( ABS(  JVS( 54 )) < TINY(a) ) THEN
         IER = 11                                      
         RETURN
  END IF
   W( 11 ) = JVS( 54 )
   W( 68 ) = JVS( 55 )
  JVS( 54) = W( 11 )
  JVS( 55) = W( 68 )
  IF ( ABS(  JVS( 56 )) < TINY(a) ) THEN
         IER = 12                                      
         RETURN
  END IF
   W( 12 ) = JVS( 56 )
   W( 68 ) = JVS( 57 )
  JVS( 56) = W( 12 )
  JVS( 57) = W( 68 )
  IF ( ABS(  JVS( 58 )) < TINY(a) ) THEN
         IER = 13                                      
         RETURN
  END IF
   W( 13 ) = JVS( 58 )
   W( 57 ) = JVS( 59 )
   W( 58 ) = JVS( 60 )
   W( 63 ) = JVS( 61 )
   W( 68 ) = JVS( 62 )
  JVS( 58) = W( 13 )
  JVS( 59) = W( 57 )
  JVS( 60) = W( 58 )
  JVS( 61) = W( 63 )
  JVS( 62) = W( 68 )
  IF ( ABS(  JVS( 63 )) < TINY(a) ) THEN
         IER = 14                                      
         RETURN
  END IF
   W( 14 ) = JVS( 63 )
   W( 57 ) = JVS( 64 )
   W( 58 ) = JVS( 65 )
   W( 63 ) = JVS( 66 )
   W( 68 ) = JVS( 67 )
  JVS( 63) = W( 14 )
  JVS( 64) = W( 57 )
  JVS( 65) = W( 58 )
  JVS( 66) = W( 63 )
  JVS( 67) = W( 68 )
  IF ( ABS(  JVS( 68 )) < TINY(a) ) THEN
         IER = 15                                      
         RETURN
  END IF
   W( 15 ) = JVS( 68 )
   W( 54 ) = JVS( 69 )
   W( 68 ) = JVS( 70 )
  JVS( 68) = W( 15 )
  JVS( 69) = W( 54 )
  JVS( 70) = W( 68 )
  IF ( ABS(  JVS( 71 )) < TINY(a) ) THEN
         IER = 16                                      
         RETURN
  END IF
   W( 16 ) = JVS( 71 )
   W( 67 ) = JVS( 72 )
   W( 68 ) = JVS( 73 )
   W( 73 ) = JVS( 74 )
  JVS( 71) = W( 16 )
  JVS( 72) = W( 67 )
  JVS( 73) = W( 68 )
  JVS( 74) = W( 73 )
  IF ( ABS(  JVS( 75 )) < TINY(a) ) THEN
         IER = 17                                      
         RETURN
  END IF
   W( 17 ) = JVS( 75 )
   W( 65 ) = JVS( 76 )
   W( 67 ) = JVS( 77 )
   W( 68 ) = JVS( 78 )
  JVS( 75) = W( 17 )
  JVS( 76) = W( 65 )
  JVS( 77) = W( 67 )
  JVS( 78) = W( 68 )
  IF ( ABS(  JVS( 79 )) < TINY(a) ) THEN
         IER = 18                                      
         RETURN
  END IF
   W( 18 ) = JVS( 79 )
   W( 21 ) = JVS( 80 )
   W( 22 ) = JVS( 81 )
   W( 25 ) = JVS( 82 )
   W( 68 ) = JVS( 83 )
   W( 72 ) = JVS( 84 )
   W( 73 ) = JVS( 85 )
  JVS( 79) = W( 18 )
  JVS( 80) = W( 21 )
  JVS( 81) = W( 22 )
  JVS( 82) = W( 25 )
  JVS( 83) = W( 68 )
  JVS( 84) = W( 72 )
  JVS( 85) = W( 73 )
  IF ( ABS(  JVS( 86 )) < TINY(a) ) THEN
         IER = 19                                      
         RETURN
  END IF
   W( 19 ) = JVS( 86 )
   W( 29 ) = JVS( 87 )
   W( 30 ) = JVS( 88 )
   W( 33 ) = JVS( 89 )
   W( 37 ) = JVS( 90 )
   W( 57 ) = JVS( 91 )
   W( 58 ) = JVS( 92 )
   W( 63 ) = JVS( 93 )
   W( 67 ) = JVS( 94 )
   W( 68 ) = JVS( 95 )
  JVS( 86) = W( 19 )
  JVS( 87) = W( 29 )
  JVS( 88) = W( 30 )
  JVS( 89) = W( 33 )
  JVS( 90) = W( 37 )
  JVS( 91) = W( 57 )
  JVS( 92) = W( 58 )
  JVS( 93) = W( 63 )
  JVS( 94) = W( 67 )
  JVS( 95) = W( 68 )
  IF ( ABS(  JVS( 96 )) < TINY(a) ) THEN
         IER = 20                                      
         RETURN
  END IF
   W( 20 ) = JVS( 96 )
   W( 32 ) = JVS( 97 )
   W( 67 ) = JVS( 98 )
   W( 68 ) = JVS( 99 )
   W( 70 ) = JVS( 100 )
   W( 73 ) = JVS( 101 )
  JVS( 96) = W( 20 )
  JVS( 97) = W( 32 )
  JVS( 98) = W( 67 )
  JVS( 99) = W( 68 )
  JVS( 100) = W( 70 )
  JVS( 101) = W( 73 )
  IF ( ABS(  JVS( 103 )) < TINY(a) ) THEN
         IER = 21                                      
         RETURN
  END IF
   W( 8 ) = JVS( 102 )
   W( 21 ) = JVS( 103 )
   W( 63 ) = JVS( 104 )
   W( 68 ) = JVS( 105 )
   W( 73 ) = JVS( 106 )
  a = -W( 8 ) / JVS(           47  )
  W( 8 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 48 )
  JVS( 102) = W( 8 )
  JVS( 103) = W( 21 )
  JVS( 104) = W( 63 )
  JVS( 105) = W( 68 )
  JVS( 106) = W( 73 )
  IF ( ABS(  JVS( 108 )) < TINY(a) ) THEN
         IER = 22                                      
         RETURN
  END IF
   W( 9 ) = JVS( 107 )
   W( 22 ) = JVS( 108 )
   W( 63 ) = JVS( 109 )
   W( 68 ) = JVS( 110 )
   W( 73 ) = JVS( 111 )
  a = -W( 9 ) / JVS(           49  )
  W( 9 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 50 )
  JVS( 107) = W( 9 )
  JVS( 108) = W( 22 )
  JVS( 109) = W( 63 )
  JVS( 110) = W( 68 )
  JVS( 111) = W( 73 )
  IF ( ABS(  JVS( 113 )) < TINY(a) ) THEN
         IER = 23                                      
         RETURN
  END IF
   W( 11 ) = JVS( 112 )
   W( 23 ) = JVS( 113 )
   W( 35 ) = JVS( 114 )
   W( 51 ) = JVS( 115 )
   W( 64 ) = JVS( 116 )
   W( 65 ) = JVS( 117 )
   W( 68 ) = JVS( 118 )
  a = -W( 11 ) / JVS(           54  )
  W( 11 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 55 )
  JVS( 112) = W( 11 )
  JVS( 113) = W( 23 )
  JVS( 114) = W( 35 )
  JVS( 115) = W( 51 )
  JVS( 116) = W( 64 )
  JVS( 117) = W( 65 )
  JVS( 118) = W( 68 )
  IF ( ABS(  JVS( 119 )) < TINY(a) ) THEN
         IER = 24                                      
         RETURN
  END IF
   W( 24 ) = JVS( 119 )
   W( 63 ) = JVS( 120 )
   W( 68 ) = JVS( 121 )
   W( 70 ) = JVS( 122 )
  JVS( 119) = W( 24 )
  JVS( 120) = W( 63 )
  JVS( 121) = W( 68 )
  JVS( 122) = W( 70 )
  IF ( ABS(  JVS( 123 )) < TINY(a) ) THEN
         IER = 25                                      
         RETURN
  END IF
   W( 25 ) = JVS( 123 )
   W( 32 ) = JVS( 124 )
   W( 63 ) = JVS( 125 )
   W( 68 ) = JVS( 126 )
   W( 73 ) = JVS( 127 )
  JVS( 123) = W( 25 )
  JVS( 124) = W( 32 )
  JVS( 125) = W( 63 )
  JVS( 126) = W( 68 )
  JVS( 127) = W( 73 )
  IF ( ABS(  JVS( 128 )) < TINY(a) ) THEN
         IER = 26                                      
         RETURN
  END IF
   W( 26 ) = JVS( 128 )
   W( 54 ) = JVS( 129 )
   W( 63 ) = JVS( 130 )
   W( 66 ) = JVS( 131 )
   W( 67 ) = JVS( 132 )
   W( 68 ) = JVS( 133 )
  JVS( 128) = W( 26 )
  JVS( 129) = W( 54 )
  JVS( 130) = W( 63 )
  JVS( 131) = W( 66 )
  JVS( 132) = W( 67 )
  JVS( 133) = W( 68 )
  IF ( ABS(  JVS( 134 )) < TINY(a) ) THEN
         IER = 27                                      
         RETURN
  END IF
   W( 27 ) = JVS( 134 )
   W( 32 ) = JVS( 135 )
   W( 34 ) = JVS( 136 )
   W( 40 ) = JVS( 137 )
   W( 48 ) = JVS( 138 )
   W( 51 ) = JVS( 139 )
   W( 53 ) = JVS( 140 )
   W( 54 ) = JVS( 141 )
   W( 67 ) = JVS( 142 )
   W( 68 ) = JVS( 143 )
   W( 70 ) = JVS( 144 )
   W( 73 ) = JVS( 145 )
  JVS( 134) = W( 27 )
  JVS( 135) = W( 32 )
  JVS( 136) = W( 34 )
  JVS( 137) = W( 40 )
  JVS( 138) = W( 48 )
  JVS( 139) = W( 51 )
  JVS( 140) = W( 53 )
  JVS( 141) = W( 54 )
  JVS( 142) = W( 67 )
  JVS( 143) = W( 68 )
  JVS( 144) = W( 70 )
  JVS( 145) = W( 73 )
  IF ( ABS(  JVS( 148 )) < TINY(a) ) THEN
         IER = 28                                      
         RETURN
  END IF
   W( 12 ) = JVS( 146 )
   W( 24 ) = JVS( 147 )
   W( 28 ) = JVS( 148 )
   W( 29 ) = JVS( 149 )
   W( 30 ) = JVS( 150 )
   W( 33 ) = JVS( 151 )
   W( 34 ) = JVS( 152 )
   W( 35 ) = JVS( 153 )
   W( 37 ) = JVS( 154 )
   W( 40 ) = JVS( 155 )
   W( 48 ) = JVS( 156 )
   W( 51 ) = JVS( 157 )
   W( 53 ) = JVS( 158 )
   W( 54 ) = JVS( 159 )
   W( 57 ) = JVS( 160 )
   W( 58 ) = JVS( 161 )
   W( 62 ) = JVS( 162 )
   W( 63 ) = JVS( 163 )
   W( 68 ) = JVS( 164 )
   W( 70 ) = JVS( 165 )
  a = -W( 12 ) / JVS(           56  )
  W( 12 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 57 )
  a = -W( 24 ) / JVS(          119  )
  W( 24 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 120 )
  W( 68 ) = W( 68 ) + a*JVS( 121 )
  W( 70 ) = W( 70 ) + a*JVS( 122 )
  JVS( 146) = W( 12 )
  JVS( 147) = W( 24 )
  JVS( 148) = W( 28 )
  JVS( 149) = W( 29 )
  JVS( 150) = W( 30 )
  JVS( 151) = W( 33 )
  JVS( 152) = W( 34 )
  JVS( 153) = W( 35 )
  JVS( 154) = W( 37 )
  JVS( 155) = W( 40 )
  JVS( 156) = W( 48 )
  JVS( 157) = W( 51 )
  JVS( 158) = W( 53 )
  JVS( 159) = W( 54 )
  JVS( 160) = W( 57 )
  JVS( 161) = W( 58 )
  JVS( 162) = W( 62 )
  JVS( 163) = W( 63 )
  JVS( 164) = W( 68 )
  JVS( 165) = W( 70 )
  IF ( ABS(  JVS( 166 )) < TINY(a) ) THEN
         IER = 29                                      
         RETURN
  END IF
   W( 29 ) = JVS( 166 )
   W( 63 ) = JVS( 167 )
   W( 68 ) = JVS( 168 )
   W( 70 ) = JVS( 169 )
  JVS( 166) = W( 29 )
  JVS( 167) = W( 63 )
  JVS( 168) = W( 68 )
  JVS( 169) = W( 70 )
  IF ( ABS(  JVS( 170 )) < TINY(a) ) THEN
         IER = 30                                      
         RETURN
  END IF
   W( 30 ) = JVS( 170 )
   W( 63 ) = JVS( 171 )
   W( 68 ) = JVS( 172 )
   W( 70 ) = JVS( 173 )
  JVS( 170) = W( 30 )
  JVS( 171) = W( 63 )
  JVS( 172) = W( 68 )
  JVS( 173) = W( 70 )
  IF ( ABS(  JVS( 174 )) < TINY(a) ) THEN
         IER = 31                                      
         RETURN
  END IF
   W( 31 ) = JVS( 174 )
   W( 35 ) = JVS( 175 )
   W( 63 ) = JVS( 176 )
   W( 66 ) = JVS( 177 )
   W( 68 ) = JVS( 178 )
   W( 70 ) = JVS( 179 )
   W( 73 ) = JVS( 180 )
  JVS( 174) = W( 31 )
  JVS( 175) = W( 35 )
  JVS( 176) = W( 63 )
  JVS( 177) = W( 66 )
  JVS( 178) = W( 68 )
  JVS( 179) = W( 70 )
  JVS( 180) = W( 73 )
  IF ( ABS(  JVS( 185 )) < TINY(a) ) THEN
         IER = 32                                      
         RETURN
  END IF
   W( 20 ) = JVS( 181 )
   W( 21 ) = JVS( 182 )
   W( 22 ) = JVS( 183 )
   W( 25 ) = JVS( 184 )
   W( 32 ) = JVS( 185 )
   W( 63 ) = JVS( 186 )
   W( 67 ) = JVS( 187 )
   W( 68 ) = JVS( 188 )
   W( 70 ) = JVS( 189 )
   W( 73 ) = JVS( 190 )
  a = -W( 20 ) / JVS(           96  )
  W( 20 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 97 )
  W( 67 ) = W( 67 ) + a*JVS( 98 )
  W( 68 ) = W( 68 ) + a*JVS( 99 )
  W( 70 ) = W( 70 ) + a*JVS( 100 )
  W( 73 ) = W( 73 ) + a*JVS( 101 )
  a = -W( 21 ) / JVS(          103  )
  W( 21 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 104 )
  W( 68 ) = W( 68 ) + a*JVS( 105 )
  W( 73 ) = W( 73 ) + a*JVS( 106 )
  a = -W( 22 ) / JVS(          108  )
  W( 22 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 109 )
  W( 68 ) = W( 68 ) + a*JVS( 110 )
  W( 73 ) = W( 73 ) + a*JVS( 111 )
  a = -W( 25 ) / JVS(          123  )
  W( 25 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 124 )
  W( 63 ) = W( 63 ) + a*JVS( 125 )
  W( 68 ) = W( 68 ) + a*JVS( 126 )
  W( 73 ) = W( 73 ) + a*JVS( 127 )
  JVS( 181) = W( 20 )
  JVS( 182) = W( 21 )
  JVS( 183) = W( 22 )
  JVS( 184) = W( 25 )
  JVS( 185) = W( 32 )
  JVS( 186) = W( 63 )
  JVS( 187) = W( 67 )
  JVS( 188) = W( 68 )
  JVS( 189) = W( 70 )
  JVS( 190) = W( 73 )
  IF ( ABS(  JVS( 191 )) < TINY(a) ) THEN
         IER = 33                                      
         RETURN
  END IF
   W( 33 ) = JVS( 191 )
   W( 63 ) = JVS( 192 )
   W( 68 ) = JVS( 193 )
   W( 70 ) = JVS( 194 )
  JVS( 191) = W( 33 )
  JVS( 192) = W( 63 )
  JVS( 193) = W( 68 )
  JVS( 194) = W( 70 )
  IF ( ABS(  JVS( 196 )) < TINY(a) ) THEN
         IER = 34                                      
         RETURN
  END IF
   W( 12 ) = JVS( 195 )
   W( 34 ) = JVS( 196 )
   W( 41 ) = JVS( 197 )
   W( 46 ) = JVS( 198 )
   W( 47 ) = JVS( 199 )
   W( 52 ) = JVS( 200 )
   W( 54 ) = JVS( 201 )
   W( 63 ) = JVS( 202 )
   W( 65 ) = JVS( 203 )
   W( 66 ) = JVS( 204 )
   W( 68 ) = JVS( 205 )
   W( 70 ) = JVS( 206 )
   W( 72 ) = JVS( 207 )
  a = -W( 12 ) / JVS(           56  )
  W( 12 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 57 )
  JVS( 195) = W( 12 )
  JVS( 196) = W( 34 )
  JVS( 197) = W( 41 )
  JVS( 198) = W( 46 )
  JVS( 199) = W( 47 )
  JVS( 200) = W( 52 )
  JVS( 201) = W( 54 )
  JVS( 202) = W( 63 )
  JVS( 203) = W( 65 )
  JVS( 204) = W( 66 )
  JVS( 205) = W( 68 )
  JVS( 206) = W( 70 )
  JVS( 207) = W( 72 )
  IF ( ABS(  JVS( 208 )) < TINY(a) ) THEN
         IER = 35                                      
         RETURN
  END IF
   W( 35 ) = JVS( 208 )
   W( 55 ) = JVS( 209 )
   W( 63 ) = JVS( 210 )
   W( 68 ) = JVS( 211 )
   W( 70 ) = JVS( 212 )
   W( 73 ) = JVS( 213 )
  JVS( 208) = W( 35 )
  JVS( 209) = W( 55 )
  JVS( 210) = W( 63 )
  JVS( 211) = W( 68 )
  JVS( 212) = W( 70 )
  JVS( 213) = W( 73 )
  IF ( ABS(  JVS( 215 )) < TINY(a) ) THEN
         IER = 36                                      
         RETURN
  END IF
   W( 24 ) = JVS( 214 )
   W( 36 ) = JVS( 215 )
   W( 63 ) = JVS( 216 )
   W( 65 ) = JVS( 217 )
   W( 66 ) = JVS( 218 )
   W( 67 ) = JVS( 219 )
   W( 68 ) = JVS( 220 )
   W( 70 ) = JVS( 221 )
   W( 72 ) = JVS( 222 )
  a = -W( 24 ) / JVS(          119  )
  W( 24 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 120 )
  W( 68 ) = W( 68 ) + a*JVS( 121 )
  W( 70 ) = W( 70 ) + a*JVS( 122 )
  JVS( 214) = W( 24 )
  JVS( 215) = W( 36 )
  JVS( 216) = W( 63 )
  JVS( 217) = W( 65 )
  JVS( 218) = W( 66 )
  JVS( 219) = W( 67 )
  JVS( 220) = W( 68 )
  JVS( 221) = W( 70 )
  JVS( 222) = W( 72 )
  IF ( ABS(  JVS( 223 )) < TINY(a) ) THEN
         IER = 37                                      
         RETURN
  END IF
   W( 37 ) = JVS( 223 )
   W( 62 ) = JVS( 224 )
   W( 63 ) = JVS( 225 )
   W( 68 ) = JVS( 226 )
   W( 70 ) = JVS( 227 )
  JVS( 223) = W( 37 )
  JVS( 224) = W( 62 )
  JVS( 225) = W( 63 )
  JVS( 226) = W( 68 )
  JVS( 227) = W( 70 )
  IF ( ABS(  JVS( 228 )) < TINY(a) ) THEN
         IER = 38                                      
         RETURN
  END IF
   W( 38 ) = JVS( 228 )
   W( 57 ) = JVS( 229 )
   W( 65 ) = JVS( 230 )
   W( 66 ) = JVS( 231 )
   W( 67 ) = JVS( 232 )
   W( 68 ) = JVS( 233 )
   W( 70 ) = JVS( 234 )
   W( 72 ) = JVS( 235 )
  JVS( 228) = W( 38 )
  JVS( 229) = W( 57 )
  JVS( 230) = W( 65 )
  JVS( 231) = W( 66 )
  JVS( 232) = W( 67 )
  JVS( 233) = W( 68 )
  JVS( 234) = W( 70 )
  JVS( 235) = W( 72 )
  IF ( ABS(  JVS( 236 )) < TINY(a) ) THEN
         IER = 39                                      
         RETURN
  END IF
   W( 39 ) = JVS( 236 )
   W( 58 ) = JVS( 237 )
   W( 65 ) = JVS( 238 )
   W( 66 ) = JVS( 239 )
   W( 67 ) = JVS( 240 )
   W( 68 ) = JVS( 241 )
   W( 70 ) = JVS( 242 )
   W( 72 ) = JVS( 243 )
  JVS( 236) = W( 39 )
  JVS( 237) = W( 58 )
  JVS( 238) = W( 65 )
  JVS( 239) = W( 66 )
  JVS( 240) = W( 67 )
  JVS( 241) = W( 68 )
  JVS( 242) = W( 70 )
  JVS( 243) = W( 72 )
  IF ( ABS(  JVS( 246 )) < TINY(a) ) THEN
         IER = 40                                      
         RETURN
  END IF
   W( 23 ) = JVS( 244 )
   W( 35 ) = JVS( 245 )
   W( 40 ) = JVS( 246 )
   W( 41 ) = JVS( 247 )
   W( 46 ) = JVS( 248 )
   W( 47 ) = JVS( 249 )
   W( 51 ) = JVS( 250 )
   W( 52 ) = JVS( 251 )
   W( 54 ) = JVS( 252 )
   W( 55 ) = JVS( 253 )
   W( 63 ) = JVS( 254 )
   W( 64 ) = JVS( 255 )
   W( 65 ) = JVS( 256 )
   W( 66 ) = JVS( 257 )
   W( 68 ) = JVS( 258 )
   W( 70 ) = JVS( 259 )
   W( 72 ) = JVS( 260 )
   W( 73 ) = JVS( 261 )
  a = -W( 23 ) / JVS(          113  )
  W( 23 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 114 )
  W( 51 ) = W( 51 ) + a*JVS( 115 )
  W( 64 ) = W( 64 ) + a*JVS( 116 )
  W( 65 ) = W( 65 ) + a*JVS( 117 )
  W( 68 ) = W( 68 ) + a*JVS( 118 )
  a = -W( 35 ) / JVS(          208  )
  W( 35 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 209 )
  W( 63 ) = W( 63 ) + a*JVS( 210 )
  W( 68 ) = W( 68 ) + a*JVS( 211 )
  W( 70 ) = W( 70 ) + a*JVS( 212 )
  W( 73 ) = W( 73 ) + a*JVS( 213 )
  JVS( 244) = W( 23 )
  JVS( 245) = W( 35 )
  JVS( 246) = W( 40 )
  JVS( 247) = W( 41 )
  JVS( 248) = W( 46 )
  JVS( 249) = W( 47 )
  JVS( 250) = W( 51 )
  JVS( 251) = W( 52 )
  JVS( 252) = W( 54 )
  JVS( 253) = W( 55 )
  JVS( 254) = W( 63 )
  JVS( 255) = W( 64 )
  JVS( 256) = W( 65 )
  JVS( 257) = W( 66 )
  JVS( 258) = W( 68 )
  JVS( 259) = W( 70 )
  JVS( 260) = W( 72 )
  JVS( 261) = W( 73 )
  IF ( ABS(  JVS( 264 )) < TINY(a) ) THEN
         IER = 41                                      
         RETURN
  END IF
   W( 25 ) = JVS( 262 )
   W( 32 ) = JVS( 263 )
   W( 41 ) = JVS( 264 )
   W( 63 ) = JVS( 265 )
   W( 65 ) = JVS( 266 )
   W( 66 ) = JVS( 267 )
   W( 67 ) = JVS( 268 )
   W( 68 ) = JVS( 269 )
   W( 70 ) = JVS( 270 )
   W( 72 ) = JVS( 271 )
   W( 73 ) = JVS( 272 )
  a = -W( 25 ) / JVS(          123  )
  W( 25 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 124 )
  W( 63 ) = W( 63 ) + a*JVS( 125 )
  W( 68 ) = W( 68 ) + a*JVS( 126 )
  W( 73 ) = W( 73 ) + a*JVS( 127 )
  a = -W( 32 ) / JVS(          185  )
  W( 32 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 186 )
  W( 67 ) = W( 67 ) + a*JVS( 187 )
  W( 68 ) = W( 68 ) + a*JVS( 188 )
  W( 70 ) = W( 70 ) + a*JVS( 189 )
  W( 73 ) = W( 73 ) + a*JVS( 190 )
  JVS( 262) = W( 25 )
  JVS( 263) = W( 32 )
  JVS( 264) = W( 41 )
  JVS( 265) = W( 63 )
  JVS( 266) = W( 65 )
  JVS( 267) = W( 66 )
  JVS( 268) = W( 67 )
  JVS( 269) = W( 68 )
  JVS( 270) = W( 70 )
  JVS( 271) = W( 72 )
  JVS( 272) = W( 73 )
  IF ( ABS(  JVS( 278 )) < TINY(a) ) THEN
         IER = 42                                      
         RETURN
  END IF
   W( 7 ) = JVS( 273 )
   W( 15 ) = JVS( 274 )
   W( 29 ) = JVS( 275 )
   W( 38 ) = JVS( 276 )
   W( 39 ) = JVS( 277 )
   W( 42 ) = JVS( 278 )
   W( 44 ) = JVS( 279 )
   W( 45 ) = JVS( 280 )
   W( 49 ) = JVS( 281 )
   W( 52 ) = JVS( 282 )
   W( 54 ) = JVS( 283 )
   W( 57 ) = JVS( 284 )
   W( 58 ) = JVS( 285 )
   W( 59 ) = JVS( 286 )
   W( 60 ) = JVS( 287 )
   W( 63 ) = JVS( 288 )
   W( 64 ) = JVS( 289 )
   W( 65 ) = JVS( 290 )
   W( 66 ) = JVS( 291 )
   W( 67 ) = JVS( 292 )
   W( 68 ) = JVS( 293 )
   W( 69 ) = JVS( 294 )
   W( 70 ) = JVS( 295 )
   W( 71 ) = JVS( 296 )
   W( 72 ) = JVS( 297 )
  a = -W( 7 ) / JVS(           45  )
  W( 7 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 46 )
  a = -W( 15 ) / JVS(           68  )
  W( 15 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 69 )
  W( 68 ) = W( 68 ) + a*JVS( 70 )
  a = -W( 29 ) / JVS(          166  )
  W( 29 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 167 )
  W( 68 ) = W( 68 ) + a*JVS( 168 )
  W( 70 ) = W( 70 ) + a*JVS( 169 )
  a = -W( 38 ) / JVS(          228  )
  W( 38 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 229 )
  W( 65 ) = W( 65 ) + a*JVS( 230 )
  W( 66 ) = W( 66 ) + a*JVS( 231 )
  W( 67 ) = W( 67 ) + a*JVS( 232 )
  W( 68 ) = W( 68 ) + a*JVS( 233 )
  W( 70 ) = W( 70 ) + a*JVS( 234 )
  W( 72 ) = W( 72 ) + a*JVS( 235 )
  a = -W( 39 ) / JVS(          236  )
  W( 39 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 237 )
  W( 65 ) = W( 65 ) + a*JVS( 238 )
  W( 66 ) = W( 66 ) + a*JVS( 239 )
  W( 67 ) = W( 67 ) + a*JVS( 240 )
  W( 68 ) = W( 68 ) + a*JVS( 241 )
  W( 70 ) = W( 70 ) + a*JVS( 242 )
  W( 72 ) = W( 72 ) + a*JVS( 243 )
  JVS( 273) = W( 7 )
  JVS( 274) = W( 15 )
  JVS( 275) = W( 29 )
  JVS( 276) = W( 38 )
  JVS( 277) = W( 39 )
  JVS( 278) = W( 42 )
  JVS( 279) = W( 44 )
  JVS( 280) = W( 45 )
  JVS( 281) = W( 49 )
  JVS( 282) = W( 52 )
  JVS( 283) = W( 54 )
  JVS( 284) = W( 57 )
  JVS( 285) = W( 58 )
  JVS( 286) = W( 59 )
  JVS( 287) = W( 60 )
  JVS( 288) = W( 63 )
  JVS( 289) = W( 64 )
  JVS( 290) = W( 65 )
  JVS( 291) = W( 66 )
  JVS( 292) = W( 67 )
  JVS( 293) = W( 68 )
  JVS( 294) = W( 69 )
  JVS( 295) = W( 70 )
  JVS( 296) = W( 71 )
  JVS( 297) = W( 72 )
  IF ( ABS(  JVS( 299 )) < TINY(a) ) THEN
         IER = 43                                      
         RETURN
  END IF
   W( 30 ) = JVS( 298 )
   W( 43 ) = JVS( 299 )
   W( 63 ) = JVS( 300 )
   W( 65 ) = JVS( 301 )
   W( 66 ) = JVS( 302 )
   W( 67 ) = JVS( 303 )
   W( 68 ) = JVS( 304 )
   W( 70 ) = JVS( 305 )
   W( 72 ) = JVS( 306 )
  a = -W( 30 ) / JVS(          170  )
  W( 30 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 171 )
  W( 68 ) = W( 68 ) + a*JVS( 172 )
  W( 70 ) = W( 70 ) + a*JVS( 173 )
  JVS( 298) = W( 30 )
  JVS( 299) = W( 43 )
  JVS( 300) = W( 63 )
  JVS( 301) = W( 65 )
  JVS( 302) = W( 66 )
  JVS( 303) = W( 67 )
  JVS( 304) = W( 68 )
  JVS( 305) = W( 70 )
  JVS( 306) = W( 72 )
  IF ( ABS(  JVS( 308 )) < TINY(a) ) THEN
         IER = 44                                      
         RETURN
  END IF
   W( 7 ) = JVS( 307 )
   W( 44 ) = JVS( 308 )
   W( 65 ) = JVS( 309 )
   W( 66 ) = JVS( 310 )
   W( 67 ) = JVS( 311 )
   W( 68 ) = JVS( 312 )
   W( 70 ) = JVS( 313 )
   W( 72 ) = JVS( 314 )
  a = -W( 7 ) / JVS(           45  )
  W( 7 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 46 )
  JVS( 307) = W( 7 )
  JVS( 308) = W( 44 )
  JVS( 309) = W( 65 )
  JVS( 310) = W( 66 )
  JVS( 311) = W( 67 )
  JVS( 312) = W( 68 )
  JVS( 313) = W( 70 )
  JVS( 314) = W( 72 )
  IF ( ABS(  JVS( 316 )) < TINY(a) ) THEN
         IER = 45                                      
         RETURN
  END IF
   W( 11 ) = JVS( 315 )
   W( 45 ) = JVS( 316 )
   W( 65 ) = JVS( 317 )
   W( 66 ) = JVS( 318 )
   W( 67 ) = JVS( 319 )
   W( 68 ) = JVS( 320 )
   W( 70 ) = JVS( 321 )
   W( 72 ) = JVS( 322 )
  a = -W( 11 ) / JVS(           54  )
  W( 11 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 55 )
  JVS( 315) = W( 11 )
  JVS( 316) = W( 45 )
  JVS( 317) = W( 65 )
  JVS( 318) = W( 66 )
  JVS( 319) = W( 67 )
  JVS( 320) = W( 68 )
  JVS( 321) = W( 70 )
  JVS( 322) = W( 72 )
  IF ( ABS(  JVS( 324 )) < TINY(a) ) THEN
         IER = 46                                      
         RETURN
  END IF
   W( 21 ) = JVS( 323 )
   W( 46 ) = JVS( 324 )
   W( 63 ) = JVS( 325 )
   W( 65 ) = JVS( 326 )
   W( 66 ) = JVS( 327 )
   W( 67 ) = JVS( 328 )
   W( 68 ) = JVS( 329 )
   W( 70 ) = JVS( 330 )
   W( 72 ) = JVS( 331 )
   W( 73 ) = JVS( 332 )
  a = -W( 21 ) / JVS(          103  )
  W( 21 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 104 )
  W( 68 ) = W( 68 ) + a*JVS( 105 )
  W( 73 ) = W( 73 ) + a*JVS( 106 )
  JVS( 323) = W( 21 )
  JVS( 324) = W( 46 )
  JVS( 325) = W( 63 )
  JVS( 326) = W( 65 )
  JVS( 327) = W( 66 )
  JVS( 328) = W( 67 )
  JVS( 329) = W( 68 )
  JVS( 330) = W( 70 )
  JVS( 331) = W( 72 )
  JVS( 332) = W( 73 )
  IF ( ABS(  JVS( 334 )) < TINY(a) ) THEN
         IER = 47                                      
         RETURN
  END IF
   W( 22 ) = JVS( 333 )
   W( 47 ) = JVS( 334 )
   W( 63 ) = JVS( 335 )
   W( 65 ) = JVS( 336 )
   W( 66 ) = JVS( 337 )
   W( 67 ) = JVS( 338 )
   W( 68 ) = JVS( 339 )
   W( 70 ) = JVS( 340 )
   W( 72 ) = JVS( 341 )
   W( 73 ) = JVS( 342 )
  a = -W( 22 ) / JVS(          108  )
  W( 22 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 109 )
  W( 68 ) = W( 68 ) + a*JVS( 110 )
  W( 73 ) = W( 73 ) + a*JVS( 111 )
  JVS( 333) = W( 22 )
  JVS( 334) = W( 47 )
  JVS( 335) = W( 63 )
  JVS( 336) = W( 65 )
  JVS( 337) = W( 66 )
  JVS( 338) = W( 67 )
  JVS( 339) = W( 68 )
  JVS( 340) = W( 70 )
  JVS( 341) = W( 72 )
  JVS( 342) = W( 73 )
  IF ( ABS(  JVS( 363 )) < TINY(a) ) THEN
         IER = 48                                      
         RETURN
  END IF
   W( 12 ) = JVS( 343 )
   W( 17 ) = JVS( 344 )
   W( 23 ) = JVS( 345 )
   W( 24 ) = JVS( 346 )
   W( 26 ) = JVS( 347 )
   W( 30 ) = JVS( 348 )
   W( 31 ) = JVS( 349 )
   W( 33 ) = JVS( 350 )
   W( 34 ) = JVS( 351 )
   W( 35 ) = JVS( 352 )
   W( 36 ) = JVS( 353 )
   W( 37 ) = JVS( 354 )
   W( 38 ) = JVS( 355 )
   W( 39 ) = JVS( 356 )
   W( 41 ) = JVS( 357 )
   W( 43 ) = JVS( 358 )
   W( 44 ) = JVS( 359 )
   W( 45 ) = JVS( 360 )
   W( 46 ) = JVS( 361 )
   W( 47 ) = JVS( 362 )
   W( 48 ) = JVS( 363 )
   W( 49 ) = JVS( 364 )
   W( 50 ) = JVS( 365 )
   W( 51 ) = JVS( 366 )
   W( 52 ) = JVS( 367 )
   W( 54 ) = JVS( 368 )
   W( 55 ) = JVS( 369 )
   W( 56 ) = JVS( 370 )
   W( 57 ) = JVS( 371 )
   W( 58 ) = JVS( 372 )
   W( 59 ) = JVS( 373 )
   W( 60 ) = JVS( 374 )
   W( 61 ) = JVS( 375 )
   W( 62 ) = JVS( 376 )
   W( 63 ) = JVS( 377 )
   W( 64 ) = JVS( 378 )
   W( 65 ) = JVS( 379 )
   W( 66 ) = JVS( 380 )
   W( 67 ) = JVS( 381 )
   W( 68 ) = JVS( 382 )
   W( 70 ) = JVS( 383 )
   W( 72 ) = JVS( 384 )
   W( 73 ) = JVS( 385 )
  a = -W( 12 ) / JVS(           56  )
  W( 12 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 57 )
  a = -W( 17 ) / JVS(           75  )
  W( 17 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 76 )
  W( 67 ) = W( 67 ) + a*JVS( 77 )
  W( 68 ) = W( 68 ) + a*JVS( 78 )
  a = -W( 23 ) / JVS(          113  )
  W( 23 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 114 )
  W( 51 ) = W( 51 ) + a*JVS( 115 )
  W( 64 ) = W( 64 ) + a*JVS( 116 )
  W( 65 ) = W( 65 ) + a*JVS( 117 )
  W( 68 ) = W( 68 ) + a*JVS( 118 )
  a = -W( 24 ) / JVS(          119  )
  W( 24 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 120 )
  W( 68 ) = W( 68 ) + a*JVS( 121 )
  W( 70 ) = W( 70 ) + a*JVS( 122 )
  a = -W( 26 ) / JVS(          128  )
  W( 26 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 129 )
  W( 63 ) = W( 63 ) + a*JVS( 130 )
  W( 66 ) = W( 66 ) + a*JVS( 131 )
  W( 67 ) = W( 67 ) + a*JVS( 132 )
  W( 68 ) = W( 68 ) + a*JVS( 133 )
  a = -W( 30 ) / JVS(          170  )
  W( 30 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 171 )
  W( 68 ) = W( 68 ) + a*JVS( 172 )
  W( 70 ) = W( 70 ) + a*JVS( 173 )
  a = -W( 31 ) / JVS(          174  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 175 )
  W( 63 ) = W( 63 ) + a*JVS( 176 )
  W( 66 ) = W( 66 ) + a*JVS( 177 )
  W( 68 ) = W( 68 ) + a*JVS( 178 )
  W( 70 ) = W( 70 ) + a*JVS( 179 )
  W( 73 ) = W( 73 ) + a*JVS( 180 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 34 ) / JVS(          196  )
  W( 34 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 197 )
  W( 46 ) = W( 46 ) + a*JVS( 198 )
  W( 47 ) = W( 47 ) + a*JVS( 199 )
  W( 52 ) = W( 52 ) + a*JVS( 200 )
  W( 54 ) = W( 54 ) + a*JVS( 201 )
  W( 63 ) = W( 63 ) + a*JVS( 202 )
  W( 65 ) = W( 65 ) + a*JVS( 203 )
  W( 66 ) = W( 66 ) + a*JVS( 204 )
  W( 68 ) = W( 68 ) + a*JVS( 205 )
  W( 70 ) = W( 70 ) + a*JVS( 206 )
  W( 72 ) = W( 72 ) + a*JVS( 207 )
  a = -W( 35 ) / JVS(          208  )
  W( 35 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 209 )
  W( 63 ) = W( 63 ) + a*JVS( 210 )
  W( 68 ) = W( 68 ) + a*JVS( 211 )
  W( 70 ) = W( 70 ) + a*JVS( 212 )
  W( 73 ) = W( 73 ) + a*JVS( 213 )
  a = -W( 36 ) / JVS(          215  )
  W( 36 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 216 )
  W( 65 ) = W( 65 ) + a*JVS( 217 )
  W( 66 ) = W( 66 ) + a*JVS( 218 )
  W( 67 ) = W( 67 ) + a*JVS( 219 )
  W( 68 ) = W( 68 ) + a*JVS( 220 )
  W( 70 ) = W( 70 ) + a*JVS( 221 )
  W( 72 ) = W( 72 ) + a*JVS( 222 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 38 ) / JVS(          228  )
  W( 38 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 229 )
  W( 65 ) = W( 65 ) + a*JVS( 230 )
  W( 66 ) = W( 66 ) + a*JVS( 231 )
  W( 67 ) = W( 67 ) + a*JVS( 232 )
  W( 68 ) = W( 68 ) + a*JVS( 233 )
  W( 70 ) = W( 70 ) + a*JVS( 234 )
  W( 72 ) = W( 72 ) + a*JVS( 235 )
  a = -W( 39 ) / JVS(          236  )
  W( 39 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 237 )
  W( 65 ) = W( 65 ) + a*JVS( 238 )
  W( 66 ) = W( 66 ) + a*JVS( 239 )
  W( 67 ) = W( 67 ) + a*JVS( 240 )
  W( 68 ) = W( 68 ) + a*JVS( 241 )
  W( 70 ) = W( 70 ) + a*JVS( 242 )
  W( 72 ) = W( 72 ) + a*JVS( 243 )
  a = -W( 41 ) / JVS(          264  )
  W( 41 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 265 )
  W( 65 ) = W( 65 ) + a*JVS( 266 )
  W( 66 ) = W( 66 ) + a*JVS( 267 )
  W( 67 ) = W( 67 ) + a*JVS( 268 )
  W( 68 ) = W( 68 ) + a*JVS( 269 )
  W( 70 ) = W( 70 ) + a*JVS( 270 )
  W( 72 ) = W( 72 ) + a*JVS( 271 )
  W( 73 ) = W( 73 ) + a*JVS( 272 )
  a = -W( 43 ) / JVS(          299  )
  W( 43 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 300 )
  W( 65 ) = W( 65 ) + a*JVS( 301 )
  W( 66 ) = W( 66 ) + a*JVS( 302 )
  W( 67 ) = W( 67 ) + a*JVS( 303 )
  W( 68 ) = W( 68 ) + a*JVS( 304 )
  W( 70 ) = W( 70 ) + a*JVS( 305 )
  W( 72 ) = W( 72 ) + a*JVS( 306 )
  a = -W( 44 ) / JVS(          308  )
  W( 44 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 309 )
  W( 66 ) = W( 66 ) + a*JVS( 310 )
  W( 67 ) = W( 67 ) + a*JVS( 311 )
  W( 68 ) = W( 68 ) + a*JVS( 312 )
  W( 70 ) = W( 70 ) + a*JVS( 313 )
  W( 72 ) = W( 72 ) + a*JVS( 314 )
  a = -W( 45 ) / JVS(          316  )
  W( 45 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 70 ) = W( 70 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 325 )
  W( 65 ) = W( 65 ) + a*JVS( 326 )
  W( 66 ) = W( 66 ) + a*JVS( 327 )
  W( 67 ) = W( 67 ) + a*JVS( 328 )
  W( 68 ) = W( 68 ) + a*JVS( 329 )
  W( 70 ) = W( 70 ) + a*JVS( 330 )
  W( 72 ) = W( 72 ) + a*JVS( 331 )
  W( 73 ) = W( 73 ) + a*JVS( 332 )
  a = -W( 47 ) / JVS(          334  )
  W( 47 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 335 )
  W( 65 ) = W( 65 ) + a*JVS( 336 )
  W( 66 ) = W( 66 ) + a*JVS( 337 )
  W( 67 ) = W( 67 ) + a*JVS( 338 )
  W( 68 ) = W( 68 ) + a*JVS( 339 )
  W( 70 ) = W( 70 ) + a*JVS( 340 )
  W( 72 ) = W( 72 ) + a*JVS( 341 )
  W( 73 ) = W( 73 ) + a*JVS( 342 )
  JVS( 343) = W( 12 )
  JVS( 344) = W( 17 )
  JVS( 345) = W( 23 )
  JVS( 346) = W( 24 )
  JVS( 347) = W( 26 )
  JVS( 348) = W( 30 )
  JVS( 349) = W( 31 )
  JVS( 350) = W( 33 )
  JVS( 351) = W( 34 )
  JVS( 352) = W( 35 )
  JVS( 353) = W( 36 )
  JVS( 354) = W( 37 )
  JVS( 355) = W( 38 )
  JVS( 356) = W( 39 )
  JVS( 357) = W( 41 )
  JVS( 358) = W( 43 )
  JVS( 359) = W( 44 )
  JVS( 360) = W( 45 )
  JVS( 361) = W( 46 )
  JVS( 362) = W( 47 )
  JVS( 363) = W( 48 )
  JVS( 364) = W( 49 )
  JVS( 365) = W( 50 )
  JVS( 366) = W( 51 )
  JVS( 367) = W( 52 )
  JVS( 368) = W( 54 )
  JVS( 369) = W( 55 )
  JVS( 370) = W( 56 )
  JVS( 371) = W( 57 )
  JVS( 372) = W( 58 )
  JVS( 373) = W( 59 )
  JVS( 374) = W( 60 )
  JVS( 375) = W( 61 )
  JVS( 376) = W( 62 )
  JVS( 377) = W( 63 )
  JVS( 378) = W( 64 )
  JVS( 379) = W( 65 )
  JVS( 380) = W( 66 )
  JVS( 381) = W( 67 )
  JVS( 382) = W( 68 )
  JVS( 383) = W( 70 )
  JVS( 384) = W( 72 )
  JVS( 385) = W( 73 )
  IF ( ABS(  JVS( 387 )) < TINY(a) ) THEN
         IER = 49                                      
         RETURN
  END IF
   W( 29 ) = JVS( 386 )
   W( 49 ) = JVS( 387 )
   W( 63 ) = JVS( 388 )
   W( 65 ) = JVS( 389 )
   W( 66 ) = JVS( 390 )
   W( 67 ) = JVS( 391 )
   W( 68 ) = JVS( 392 )
   W( 70 ) = JVS( 393 )
   W( 72 ) = JVS( 394 )
  a = -W( 29 ) / JVS(          166  )
  W( 29 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 167 )
  W( 68 ) = W( 68 ) + a*JVS( 168 )
  W( 70 ) = W( 70 ) + a*JVS( 169 )
  JVS( 386) = W( 29 )
  JVS( 387) = W( 49 )
  JVS( 388) = W( 63 )
  JVS( 389) = W( 65 )
  JVS( 390) = W( 66 )
  JVS( 391) = W( 67 )
  JVS( 392) = W( 68 )
  JVS( 393) = W( 70 )
  JVS( 394) = W( 72 )
  IF ( ABS(  JVS( 397 )) < TINY(a) ) THEN
         IER = 50                                      
         RETURN
  END IF
   W( 33 ) = JVS( 395 )
   W( 37 ) = JVS( 396 )
   W( 50 ) = JVS( 397 )
   W( 62 ) = JVS( 398 )
   W( 63 ) = JVS( 399 )
   W( 65 ) = JVS( 400 )
   W( 66 ) = JVS( 401 )
   W( 67 ) = JVS( 402 )
   W( 68 ) = JVS( 403 )
   W( 70 ) = JVS( 404 )
   W( 72 ) = JVS( 405 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  JVS( 395) = W( 33 )
  JVS( 396) = W( 37 )
  JVS( 397) = W( 50 )
  JVS( 398) = W( 62 )
  JVS( 399) = W( 63 )
  JVS( 400) = W( 65 )
  JVS( 401) = W( 66 )
  JVS( 402) = W( 67 )
  JVS( 403) = W( 68 )
  JVS( 404) = W( 70 )
  JVS( 405) = W( 72 )
  IF ( ABS(  JVS( 411 )) < TINY(a) ) THEN
         IER = 51                                      
         RETURN
  END IF
   W( 30 ) = JVS( 406 )
   W( 33 ) = JVS( 407 )
   W( 37 ) = JVS( 408 )
   W( 43 ) = JVS( 409 )
   W( 50 ) = JVS( 410 )
   W( 51 ) = JVS( 411 )
   W( 62 ) = JVS( 412 )
   W( 63 ) = JVS( 413 )
   W( 65 ) = JVS( 414 )
   W( 66 ) = JVS( 415 )
   W( 67 ) = JVS( 416 )
   W( 68 ) = JVS( 417 )
   W( 70 ) = JVS( 418 )
   W( 72 ) = JVS( 419 )
  a = -W( 30 ) / JVS(          170  )
  W( 30 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 171 )
  W( 68 ) = W( 68 ) + a*JVS( 172 )
  W( 70 ) = W( 70 ) + a*JVS( 173 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 43 ) / JVS(          299  )
  W( 43 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 300 )
  W( 65 ) = W( 65 ) + a*JVS( 301 )
  W( 66 ) = W( 66 ) + a*JVS( 302 )
  W( 67 ) = W( 67 ) + a*JVS( 303 )
  W( 68 ) = W( 68 ) + a*JVS( 304 )
  W( 70 ) = W( 70 ) + a*JVS( 305 )
  W( 72 ) = W( 72 ) + a*JVS( 306 )
  a = -W( 50 ) / JVS(          397  )
  W( 50 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 398 )
  W( 63 ) = W( 63 ) + a*JVS( 399 )
  W( 65 ) = W( 65 ) + a*JVS( 400 )
  W( 66 ) = W( 66 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 70 ) = W( 70 ) + a*JVS( 404 )
  W( 72 ) = W( 72 ) + a*JVS( 405 )
  JVS( 406) = W( 30 )
  JVS( 407) = W( 33 )
  JVS( 408) = W( 37 )
  JVS( 409) = W( 43 )
  JVS( 410) = W( 50 )
  JVS( 411) = W( 51 )
  JVS( 412) = W( 62 )
  JVS( 413) = W( 63 )
  JVS( 414) = W( 65 )
  JVS( 415) = W( 66 )
  JVS( 416) = W( 67 )
  JVS( 417) = W( 68 )
  JVS( 418) = W( 70 )
  JVS( 419) = W( 72 )
  IF ( ABS(  JVS( 421 )) < TINY(a) ) THEN
         IER = 52                                      
         RETURN
  END IF
   W( 12 ) = JVS( 420 )
   W( 52 ) = JVS( 421 )
   W( 65 ) = JVS( 422 )
   W( 66 ) = JVS( 423 )
   W( 67 ) = JVS( 424 )
   W( 68 ) = JVS( 425 )
   W( 69 ) = JVS( 426 )
   W( 70 ) = JVS( 427 )
   W( 71 ) = JVS( 428 )
   W( 72 ) = JVS( 429 )
  a = -W( 12 ) / JVS(           56  )
  W( 12 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 57 )
  JVS( 420) = W( 12 )
  JVS( 421) = W( 52 )
  JVS( 422) = W( 65 )
  JVS( 423) = W( 66 )
  JVS( 424) = W( 67 )
  JVS( 425) = W( 68 )
  JVS( 426) = W( 69 )
  JVS( 427) = W( 70 )
  JVS( 428) = W( 71 )
  JVS( 429) = W( 72 )
  IF ( ABS(  JVS( 442 )) < TINY(a) ) THEN
         IER = 53                                      
         RETURN
  END IF
   W( 11 ) = JVS( 430 )
   W( 12 ) = JVS( 431 )
   W( 15 ) = JVS( 432 )
   W( 29 ) = JVS( 433 )
   W( 36 ) = JVS( 434 )
   W( 38 ) = JVS( 435 )
   W( 39 ) = JVS( 436 )
   W( 44 ) = JVS( 437 )
   W( 45 ) = JVS( 438 )
   W( 49 ) = JVS( 439 )
   W( 51 ) = JVS( 440 )
   W( 52 ) = JVS( 441 )
   W( 53 ) = JVS( 442 )
   W( 54 ) = JVS( 443 )
   W( 57 ) = JVS( 444 )
   W( 58 ) = JVS( 445 )
   W( 59 ) = JVS( 446 )
   W( 60 ) = JVS( 447 )
   W( 61 ) = JVS( 448 )
   W( 62 ) = JVS( 449 )
   W( 63 ) = JVS( 450 )
   W( 64 ) = JVS( 451 )
   W( 65 ) = JVS( 452 )
   W( 66 ) = JVS( 453 )
   W( 67 ) = JVS( 454 )
   W( 68 ) = JVS( 455 )
   W( 69 ) = JVS( 456 )
   W( 70 ) = JVS( 457 )
   W( 71 ) = JVS( 458 )
   W( 72 ) = JVS( 459 )
  a = -W( 11 ) / JVS(           54  )
  W( 11 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 55 )
  a = -W( 12 ) / JVS(           56  )
  W( 12 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 57 )
  a = -W( 15 ) / JVS(           68  )
  W( 15 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 69 )
  W( 68 ) = W( 68 ) + a*JVS( 70 )
  a = -W( 29 ) / JVS(          166  )
  W( 29 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 167 )
  W( 68 ) = W( 68 ) + a*JVS( 168 )
  W( 70 ) = W( 70 ) + a*JVS( 169 )
  a = -W( 36 ) / JVS(          215  )
  W( 36 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 216 )
  W( 65 ) = W( 65 ) + a*JVS( 217 )
  W( 66 ) = W( 66 ) + a*JVS( 218 )
  W( 67 ) = W( 67 ) + a*JVS( 219 )
  W( 68 ) = W( 68 ) + a*JVS( 220 )
  W( 70 ) = W( 70 ) + a*JVS( 221 )
  W( 72 ) = W( 72 ) + a*JVS( 222 )
  a = -W( 38 ) / JVS(          228  )
  W( 38 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 229 )
  W( 65 ) = W( 65 ) + a*JVS( 230 )
  W( 66 ) = W( 66 ) + a*JVS( 231 )
  W( 67 ) = W( 67 ) + a*JVS( 232 )
  W( 68 ) = W( 68 ) + a*JVS( 233 )
  W( 70 ) = W( 70 ) + a*JVS( 234 )
  W( 72 ) = W( 72 ) + a*JVS( 235 )
  a = -W( 39 ) / JVS(          236  )
  W( 39 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 237 )
  W( 65 ) = W( 65 ) + a*JVS( 238 )
  W( 66 ) = W( 66 ) + a*JVS( 239 )
  W( 67 ) = W( 67 ) + a*JVS( 240 )
  W( 68 ) = W( 68 ) + a*JVS( 241 )
  W( 70 ) = W( 70 ) + a*JVS( 242 )
  W( 72 ) = W( 72 ) + a*JVS( 243 )
  a = -W( 44 ) / JVS(          308  )
  W( 44 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 309 )
  W( 66 ) = W( 66 ) + a*JVS( 310 )
  W( 67 ) = W( 67 ) + a*JVS( 311 )
  W( 68 ) = W( 68 ) + a*JVS( 312 )
  W( 70 ) = W( 70 ) + a*JVS( 313 )
  W( 72 ) = W( 72 ) + a*JVS( 314 )
  a = -W( 45 ) / JVS(          316  )
  W( 45 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 70 ) = W( 70 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  a = -W( 49 ) / JVS(          387  )
  W( 49 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 388 )
  W( 65 ) = W( 65 ) + a*JVS( 389 )
  W( 66 ) = W( 66 ) + a*JVS( 390 )
  W( 67 ) = W( 67 ) + a*JVS( 391 )
  W( 68 ) = W( 68 ) + a*JVS( 392 )
  W( 70 ) = W( 70 ) + a*JVS( 393 )
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  a = -W( 51 ) / JVS(          411  )
  W( 51 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 412 )
  W( 63 ) = W( 63 ) + a*JVS( 413 )
  W( 65 ) = W( 65 ) + a*JVS( 414 )
  W( 66 ) = W( 66 ) + a*JVS( 415 )
  W( 67 ) = W( 67 ) + a*JVS( 416 )
  W( 68 ) = W( 68 ) + a*JVS( 417 )
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 72 ) = W( 72 ) + a*JVS( 419 )
  a = -W( 52 ) / JVS(          421  )
  W( 52 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 422 )
  W( 66 ) = W( 66 ) + a*JVS( 423 )
  W( 67 ) = W( 67 ) + a*JVS( 424 )
  W( 68 ) = W( 68 ) + a*JVS( 425 )
  W( 69 ) = W( 69 ) + a*JVS( 426 )
  W( 70 ) = W( 70 ) + a*JVS( 427 )
  W( 71 ) = W( 71 ) + a*JVS( 428 )
  W( 72 ) = W( 72 ) + a*JVS( 429 )
  JVS( 430) = W( 11 )
  JVS( 431) = W( 12 )
  JVS( 432) = W( 15 )
  JVS( 433) = W( 29 )
  JVS( 434) = W( 36 )
  JVS( 435) = W( 38 )
  JVS( 436) = W( 39 )
  JVS( 437) = W( 44 )
  JVS( 438) = W( 45 )
  JVS( 439) = W( 49 )
  JVS( 440) = W( 51 )
  JVS( 441) = W( 52 )
  JVS( 442) = W( 53 )
  JVS( 443) = W( 54 )
  JVS( 444) = W( 57 )
  JVS( 445) = W( 58 )
  JVS( 446) = W( 59 )
  JVS( 447) = W( 60 )
  JVS( 448) = W( 61 )
  JVS( 449) = W( 62 )
  JVS( 450) = W( 63 )
  JVS( 451) = W( 64 )
  JVS( 452) = W( 65 )
  JVS( 453) = W( 66 )
  JVS( 454) = W( 67 )
  JVS( 455) = W( 68 )
  JVS( 456) = W( 69 )
  JVS( 457) = W( 70 )
  JVS( 458) = W( 71 )
  JVS( 459) = W( 72 )
  IF ( ABS(  JVS( 463 )) < TINY(a) ) THEN
         IER = 54                                      
         RETURN
  END IF
   W( 37 ) = JVS( 460 )
   W( 46 ) = JVS( 461 )
   W( 47 ) = JVS( 462 )
   W( 54 ) = JVS( 463 )
   W( 62 ) = JVS( 464 )
   W( 63 ) = JVS( 465 )
   W( 65 ) = JVS( 466 )
   W( 66 ) = JVS( 467 )
   W( 67 ) = JVS( 468 )
   W( 68 ) = JVS( 469 )
   W( 70 ) = JVS( 470 )
   W( 72 ) = JVS( 471 )
   W( 73 ) = JVS( 472 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 325 )
  W( 65 ) = W( 65 ) + a*JVS( 326 )
  W( 66 ) = W( 66 ) + a*JVS( 327 )
  W( 67 ) = W( 67 ) + a*JVS( 328 )
  W( 68 ) = W( 68 ) + a*JVS( 329 )
  W( 70 ) = W( 70 ) + a*JVS( 330 )
  W( 72 ) = W( 72 ) + a*JVS( 331 )
  W( 73 ) = W( 73 ) + a*JVS( 332 )
  a = -W( 47 ) / JVS(          334  )
  W( 47 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 335 )
  W( 65 ) = W( 65 ) + a*JVS( 336 )
  W( 66 ) = W( 66 ) + a*JVS( 337 )
  W( 67 ) = W( 67 ) + a*JVS( 338 )
  W( 68 ) = W( 68 ) + a*JVS( 339 )
  W( 70 ) = W( 70 ) + a*JVS( 340 )
  W( 72 ) = W( 72 ) + a*JVS( 341 )
  W( 73 ) = W( 73 ) + a*JVS( 342 )
  JVS( 460) = W( 37 )
  JVS( 461) = W( 46 )
  JVS( 462) = W( 47 )
  JVS( 463) = W( 54 )
  JVS( 464) = W( 62 )
  JVS( 465) = W( 63 )
  JVS( 466) = W( 65 )
  JVS( 467) = W( 66 )
  JVS( 468) = W( 67 )
  JVS( 469) = W( 68 )
  JVS( 470) = W( 70 )
  JVS( 471) = W( 72 )
  JVS( 472) = W( 73 )
  IF ( ABS(  JVS( 476 )) < TINY(a) ) THEN
         IER = 55                                      
         RETURN
  END IF
   W( 35 ) = JVS( 473 )
   W( 51 ) = JVS( 474 )
   W( 54 ) = JVS( 475 )
   W( 55 ) = JVS( 476 )
   W( 62 ) = JVS( 477 )
   W( 63 ) = JVS( 478 )
   W( 65 ) = JVS( 479 )
   W( 66 ) = JVS( 480 )
   W( 67 ) = JVS( 481 )
   W( 68 ) = JVS( 482 )
   W( 70 ) = JVS( 483 )
   W( 72 ) = JVS( 484 )
   W( 73 ) = JVS( 485 )
  a = -W( 35 ) / JVS(          208  )
  W( 35 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 209 )
  W( 63 ) = W( 63 ) + a*JVS( 210 )
  W( 68 ) = W( 68 ) + a*JVS( 211 )
  W( 70 ) = W( 70 ) + a*JVS( 212 )
  W( 73 ) = W( 73 ) + a*JVS( 213 )
  a = -W( 51 ) / JVS(          411  )
  W( 51 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 412 )
  W( 63 ) = W( 63 ) + a*JVS( 413 )
  W( 65 ) = W( 65 ) + a*JVS( 414 )
  W( 66 ) = W( 66 ) + a*JVS( 415 )
  W( 67 ) = W( 67 ) + a*JVS( 416 )
  W( 68 ) = W( 68 ) + a*JVS( 417 )
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 72 ) = W( 72 ) + a*JVS( 419 )
  a = -W( 54 ) / JVS(          463  )
  W( 54 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 464 )
  W( 63 ) = W( 63 ) + a*JVS( 465 )
  W( 65 ) = W( 65 ) + a*JVS( 466 )
  W( 66 ) = W( 66 ) + a*JVS( 467 )
  W( 67 ) = W( 67 ) + a*JVS( 468 )
  W( 68 ) = W( 68 ) + a*JVS( 469 )
  W( 70 ) = W( 70 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  JVS( 473) = W( 35 )
  JVS( 474) = W( 51 )
  JVS( 475) = W( 54 )
  JVS( 476) = W( 55 )
  JVS( 477) = W( 62 )
  JVS( 478) = W( 63 )
  JVS( 479) = W( 65 )
  JVS( 480) = W( 66 )
  JVS( 481) = W( 67 )
  JVS( 482) = W( 68 )
  JVS( 483) = W( 70 )
  JVS( 484) = W( 72 )
  JVS( 485) = W( 73 )
  IF ( ABS(  JVS( 500 )) < TINY(a) ) THEN
         IER = 56                                      
         RETURN
  END IF
   W( 8 ) = JVS( 486 )
   W( 9 ) = JVS( 487 )
   W( 26 ) = JVS( 488 )
   W( 31 ) = JVS( 489 )
   W( 32 ) = JVS( 490 )
   W( 33 ) = JVS( 491 )
   W( 35 ) = JVS( 492 )
   W( 37 ) = JVS( 493 )
   W( 44 ) = JVS( 494 )
   W( 45 ) = JVS( 495 )
   W( 51 ) = JVS( 496 )
   W( 52 ) = JVS( 497 )
   W( 54 ) = JVS( 498 )
   W( 55 ) = JVS( 499 )
   W( 56 ) = JVS( 500 )
   W( 62 ) = JVS( 501 )
   W( 63 ) = JVS( 502 )
   W( 64 ) = JVS( 503 )
   W( 65 ) = JVS( 504 )
   W( 66 ) = JVS( 505 )
   W( 67 ) = JVS( 506 )
   W( 68 ) = JVS( 507 )
   W( 69 ) = JVS( 508 )
   W( 70 ) = JVS( 509 )
   W( 71 ) = JVS( 510 )
   W( 72 ) = JVS( 511 )
   W( 73 ) = JVS( 512 )
  a = -W( 8 ) / JVS(           47  )
  W( 8 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 48 )
  a = -W( 9 ) / JVS(           49  )
  W( 9 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 50 )
  a = -W( 26 ) / JVS(          128  )
  W( 26 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 129 )
  W( 63 ) = W( 63 ) + a*JVS( 130 )
  W( 66 ) = W( 66 ) + a*JVS( 131 )
  W( 67 ) = W( 67 ) + a*JVS( 132 )
  W( 68 ) = W( 68 ) + a*JVS( 133 )
  a = -W( 31 ) / JVS(          174  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 175 )
  W( 63 ) = W( 63 ) + a*JVS( 176 )
  W( 66 ) = W( 66 ) + a*JVS( 177 )
  W( 68 ) = W( 68 ) + a*JVS( 178 )
  W( 70 ) = W( 70 ) + a*JVS( 179 )
  W( 73 ) = W( 73 ) + a*JVS( 180 )
  a = -W( 32 ) / JVS(          185  )
  W( 32 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 186 )
  W( 67 ) = W( 67 ) + a*JVS( 187 )
  W( 68 ) = W( 68 ) + a*JVS( 188 )
  W( 70 ) = W( 70 ) + a*JVS( 189 )
  W( 73 ) = W( 73 ) + a*JVS( 190 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 35 ) / JVS(          208  )
  W( 35 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 209 )
  W( 63 ) = W( 63 ) + a*JVS( 210 )
  W( 68 ) = W( 68 ) + a*JVS( 211 )
  W( 70 ) = W( 70 ) + a*JVS( 212 )
  W( 73 ) = W( 73 ) + a*JVS( 213 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 44 ) / JVS(          308  )
  W( 44 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 309 )
  W( 66 ) = W( 66 ) + a*JVS( 310 )
  W( 67 ) = W( 67 ) + a*JVS( 311 )
  W( 68 ) = W( 68 ) + a*JVS( 312 )
  W( 70 ) = W( 70 ) + a*JVS( 313 )
  W( 72 ) = W( 72 ) + a*JVS( 314 )
  a = -W( 45 ) / JVS(          316  )
  W( 45 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 70 ) = W( 70 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  a = -W( 51 ) / JVS(          411  )
  W( 51 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 412 )
  W( 63 ) = W( 63 ) + a*JVS( 413 )
  W( 65 ) = W( 65 ) + a*JVS( 414 )
  W( 66 ) = W( 66 ) + a*JVS( 415 )
  W( 67 ) = W( 67 ) + a*JVS( 416 )
  W( 68 ) = W( 68 ) + a*JVS( 417 )
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 72 ) = W( 72 ) + a*JVS( 419 )
  a = -W( 52 ) / JVS(          421  )
  W( 52 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 422 )
  W( 66 ) = W( 66 ) + a*JVS( 423 )
  W( 67 ) = W( 67 ) + a*JVS( 424 )
  W( 68 ) = W( 68 ) + a*JVS( 425 )
  W( 69 ) = W( 69 ) + a*JVS( 426 )
  W( 70 ) = W( 70 ) + a*JVS( 427 )
  W( 71 ) = W( 71 ) + a*JVS( 428 )
  W( 72 ) = W( 72 ) + a*JVS( 429 )
  a = -W( 54 ) / JVS(          463  )
  W( 54 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 464 )
  W( 63 ) = W( 63 ) + a*JVS( 465 )
  W( 65 ) = W( 65 ) + a*JVS( 466 )
  W( 66 ) = W( 66 ) + a*JVS( 467 )
  W( 67 ) = W( 67 ) + a*JVS( 468 )
  W( 68 ) = W( 68 ) + a*JVS( 469 )
  W( 70 ) = W( 70 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  a = -W( 55 ) / JVS(          476  )
  W( 55 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 477 )
  W( 63 ) = W( 63 ) + a*JVS( 478 )
  W( 65 ) = W( 65 ) + a*JVS( 479 )
  W( 66 ) = W( 66 ) + a*JVS( 480 )
  W( 67 ) = W( 67 ) + a*JVS( 481 )
  W( 68 ) = W( 68 ) + a*JVS( 482 )
  W( 70 ) = W( 70 ) + a*JVS( 483 )
  W( 72 ) = W( 72 ) + a*JVS( 484 )
  W( 73 ) = W( 73 ) + a*JVS( 485 )
  JVS( 486) = W( 8 )
  JVS( 487) = W( 9 )
  JVS( 488) = W( 26 )
  JVS( 489) = W( 31 )
  JVS( 490) = W( 32 )
  JVS( 491) = W( 33 )
  JVS( 492) = W( 35 )
  JVS( 493) = W( 37 )
  JVS( 494) = W( 44 )
  JVS( 495) = W( 45 )
  JVS( 496) = W( 51 )
  JVS( 497) = W( 52 )
  JVS( 498) = W( 54 )
  JVS( 499) = W( 55 )
  JVS( 500) = W( 56 )
  JVS( 501) = W( 62 )
  JVS( 502) = W( 63 )
  JVS( 503) = W( 64 )
  JVS( 504) = W( 65 )
  JVS( 505) = W( 66 )
  JVS( 506) = W( 67 )
  JVS( 507) = W( 68 )
  JVS( 508) = W( 69 )
  JVS( 509) = W( 70 )
  JVS( 510) = W( 71 )
  JVS( 511) = W( 72 )
  JVS( 512) = W( 73 )
  IF ( ABS(  JVS( 517 )) < TINY(a) ) THEN
         IER = 57                                      
         RETURN
  END IF
   W( 30 ) = JVS( 513 )
   W( 33 ) = JVS( 514 )
   W( 37 ) = JVS( 515 )
   W( 50 ) = JVS( 516 )
   W( 57 ) = JVS( 517 )
   W( 62 ) = JVS( 518 )
   W( 63 ) = JVS( 519 )
   W( 65 ) = JVS( 520 )
   W( 66 ) = JVS( 521 )
   W( 67 ) = JVS( 522 )
   W( 68 ) = JVS( 523 )
   W( 70 ) = JVS( 524 )
   W( 72 ) = JVS( 525 )
  a = -W( 30 ) / JVS(          170  )
  W( 30 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 171 )
  W( 68 ) = W( 68 ) + a*JVS( 172 )
  W( 70 ) = W( 70 ) + a*JVS( 173 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 50 ) / JVS(          397  )
  W( 50 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 398 )
  W( 63 ) = W( 63 ) + a*JVS( 399 )
  W( 65 ) = W( 65 ) + a*JVS( 400 )
  W( 66 ) = W( 66 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 70 ) = W( 70 ) + a*JVS( 404 )
  W( 72 ) = W( 72 ) + a*JVS( 405 )
  JVS( 513) = W( 30 )
  JVS( 514) = W( 33 )
  JVS( 515) = W( 37 )
  JVS( 516) = W( 50 )
  JVS( 517) = W( 57 )
  JVS( 518) = W( 62 )
  JVS( 519) = W( 63 )
  JVS( 520) = W( 65 )
  JVS( 521) = W( 66 )
  JVS( 522) = W( 67 )
  JVS( 523) = W( 68 )
  JVS( 524) = W( 70 )
  JVS( 525) = W( 72 )
  IF ( ABS(  JVS( 528 )) < TINY(a) ) THEN
         IER = 58                                      
         RETURN
  END IF
   W( 43 ) = JVS( 526 )
   W( 50 ) = JVS( 527 )
   W( 58 ) = JVS( 528 )
   W( 62 ) = JVS( 529 )
   W( 63 ) = JVS( 530 )
   W( 65 ) = JVS( 531 )
   W( 66 ) = JVS( 532 )
   W( 67 ) = JVS( 533 )
   W( 68 ) = JVS( 534 )
   W( 70 ) = JVS( 535 )
   W( 72 ) = JVS( 536 )
  a = -W( 43 ) / JVS(          299  )
  W( 43 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 300 )
  W( 65 ) = W( 65 ) + a*JVS( 301 )
  W( 66 ) = W( 66 ) + a*JVS( 302 )
  W( 67 ) = W( 67 ) + a*JVS( 303 )
  W( 68 ) = W( 68 ) + a*JVS( 304 )
  W( 70 ) = W( 70 ) + a*JVS( 305 )
  W( 72 ) = W( 72 ) + a*JVS( 306 )
  a = -W( 50 ) / JVS(          397  )
  W( 50 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 398 )
  W( 63 ) = W( 63 ) + a*JVS( 399 )
  W( 65 ) = W( 65 ) + a*JVS( 400 )
  W( 66 ) = W( 66 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 70 ) = W( 70 ) + a*JVS( 404 )
  W( 72 ) = W( 72 ) + a*JVS( 405 )
  JVS( 526) = W( 43 )
  JVS( 527) = W( 50 )
  JVS( 528) = W( 58 )
  JVS( 529) = W( 62 )
  JVS( 530) = W( 63 )
  JVS( 531) = W( 65 )
  JVS( 532) = W( 66 )
  JVS( 533) = W( 67 )
  JVS( 534) = W( 68 )
  JVS( 535) = W( 70 )
  JVS( 536) = W( 72 )
  IF ( ABS(  JVS( 545 )) < TINY(a) ) THEN
         IER = 59                                      
         RETURN
  END IF
   W( 24 ) = JVS( 537 )
   W( 29 ) = JVS( 538 )
   W( 30 ) = JVS( 539 )
   W( 33 ) = JVS( 540 )
   W( 37 ) = JVS( 541 )
   W( 51 ) = JVS( 542 )
   W( 57 ) = JVS( 543 )
   W( 58 ) = JVS( 544 )
   W( 59 ) = JVS( 545 )
   W( 60 ) = JVS( 546 )
   W( 62 ) = JVS( 547 )
   W( 63 ) = JVS( 548 )
   W( 65 ) = JVS( 549 )
   W( 66 ) = JVS( 550 )
   W( 67 ) = JVS( 551 )
   W( 68 ) = JVS( 552 )
   W( 70 ) = JVS( 553 )
   W( 72 ) = JVS( 554 )
  a = -W( 24 ) / JVS(          119  )
  W( 24 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 120 )
  W( 68 ) = W( 68 ) + a*JVS( 121 )
  W( 70 ) = W( 70 ) + a*JVS( 122 )
  a = -W( 29 ) / JVS(          166  )
  W( 29 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 167 )
  W( 68 ) = W( 68 ) + a*JVS( 168 )
  W( 70 ) = W( 70 ) + a*JVS( 169 )
  a = -W( 30 ) / JVS(          170  )
  W( 30 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 171 )
  W( 68 ) = W( 68 ) + a*JVS( 172 )
  W( 70 ) = W( 70 ) + a*JVS( 173 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 51 ) / JVS(          411  )
  W( 51 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 412 )
  W( 63 ) = W( 63 ) + a*JVS( 413 )
  W( 65 ) = W( 65 ) + a*JVS( 414 )
  W( 66 ) = W( 66 ) + a*JVS( 415 )
  W( 67 ) = W( 67 ) + a*JVS( 416 )
  W( 68 ) = W( 68 ) + a*JVS( 417 )
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 72 ) = W( 72 ) + a*JVS( 419 )
  a = -W( 57 ) / JVS(          517  )
  W( 57 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 518 )
  W( 63 ) = W( 63 ) + a*JVS( 519 )
  W( 65 ) = W( 65 ) + a*JVS( 520 )
  W( 66 ) = W( 66 ) + a*JVS( 521 )
  W( 67 ) = W( 67 ) + a*JVS( 522 )
  W( 68 ) = W( 68 ) + a*JVS( 523 )
  W( 70 ) = W( 70 ) + a*JVS( 524 )
  W( 72 ) = W( 72 ) + a*JVS( 525 )
  a = -W( 58 ) / JVS(          528  )
  W( 58 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 529 )
  W( 63 ) = W( 63 ) + a*JVS( 530 )
  W( 65 ) = W( 65 ) + a*JVS( 531 )
  W( 66 ) = W( 66 ) + a*JVS( 532 )
  W( 67 ) = W( 67 ) + a*JVS( 533 )
  W( 68 ) = W( 68 ) + a*JVS( 534 )
  W( 70 ) = W( 70 ) + a*JVS( 535 )
  W( 72 ) = W( 72 ) + a*JVS( 536 )
  JVS( 537) = W( 24 )
  JVS( 538) = W( 29 )
  JVS( 539) = W( 30 )
  JVS( 540) = W( 33 )
  JVS( 541) = W( 37 )
  JVS( 542) = W( 51 )
  JVS( 543) = W( 57 )
  JVS( 544) = W( 58 )
  JVS( 545) = W( 59 )
  JVS( 546) = W( 60 )
  JVS( 547) = W( 62 )
  JVS( 548) = W( 63 )
  JVS( 549) = W( 65 )
  JVS( 550) = W( 66 )
  JVS( 551) = W( 67 )
  JVS( 552) = W( 68 )
  JVS( 553) = W( 70 )
  JVS( 554) = W( 72 )
  IF ( ABS(  JVS( 563 )) < TINY(a) ) THEN
         IER = 60                                      
         RETURN
  END IF
   W( 24 ) = JVS( 555 )
   W( 29 ) = JVS( 556 )
   W( 30 ) = JVS( 557 )
   W( 33 ) = JVS( 558 )
   W( 37 ) = JVS( 559 )
   W( 57 ) = JVS( 560 )
   W( 58 ) = JVS( 561 )
   W( 59 ) = JVS( 562 )
   W( 60 ) = JVS( 563 )
   W( 62 ) = JVS( 564 )
   W( 63 ) = JVS( 565 )
   W( 65 ) = JVS( 566 )
   W( 66 ) = JVS( 567 )
   W( 67 ) = JVS( 568 )
   W( 68 ) = JVS( 569 )
   W( 70 ) = JVS( 570 )
   W( 72 ) = JVS( 571 )
  a = -W( 24 ) / JVS(          119  )
  W( 24 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 120 )
  W( 68 ) = W( 68 ) + a*JVS( 121 )
  W( 70 ) = W( 70 ) + a*JVS( 122 )
  a = -W( 29 ) / JVS(          166  )
  W( 29 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 167 )
  W( 68 ) = W( 68 ) + a*JVS( 168 )
  W( 70 ) = W( 70 ) + a*JVS( 169 )
  a = -W( 30 ) / JVS(          170  )
  W( 30 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 171 )
  W( 68 ) = W( 68 ) + a*JVS( 172 )
  W( 70 ) = W( 70 ) + a*JVS( 173 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 57 ) / JVS(          517  )
  W( 57 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 518 )
  W( 63 ) = W( 63 ) + a*JVS( 519 )
  W( 65 ) = W( 65 ) + a*JVS( 520 )
  W( 66 ) = W( 66 ) + a*JVS( 521 )
  W( 67 ) = W( 67 ) + a*JVS( 522 )
  W( 68 ) = W( 68 ) + a*JVS( 523 )
  W( 70 ) = W( 70 ) + a*JVS( 524 )
  W( 72 ) = W( 72 ) + a*JVS( 525 )
  a = -W( 58 ) / JVS(          528  )
  W( 58 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 529 )
  W( 63 ) = W( 63 ) + a*JVS( 530 )
  W( 65 ) = W( 65 ) + a*JVS( 531 )
  W( 66 ) = W( 66 ) + a*JVS( 532 )
  W( 67 ) = W( 67 ) + a*JVS( 533 )
  W( 68 ) = W( 68 ) + a*JVS( 534 )
  W( 70 ) = W( 70 ) + a*JVS( 535 )
  W( 72 ) = W( 72 ) + a*JVS( 536 )
  a = -W( 59 ) / JVS(          545  )
  W( 59 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 546 )
  W( 62 ) = W( 62 ) + a*JVS( 547 )
  W( 63 ) = W( 63 ) + a*JVS( 548 )
  W( 65 ) = W( 65 ) + a*JVS( 549 )
  W( 66 ) = W( 66 ) + a*JVS( 550 )
  W( 67 ) = W( 67 ) + a*JVS( 551 )
  W( 68 ) = W( 68 ) + a*JVS( 552 )
  W( 70 ) = W( 70 ) + a*JVS( 553 )
  W( 72 ) = W( 72 ) + a*JVS( 554 )
  JVS( 555) = W( 24 )
  JVS( 556) = W( 29 )
  JVS( 557) = W( 30 )
  JVS( 558) = W( 33 )
  JVS( 559) = W( 37 )
  JVS( 560) = W( 57 )
  JVS( 561) = W( 58 )
  JVS( 562) = W( 59 )
  JVS( 563) = W( 60 )
  JVS( 564) = W( 62 )
  JVS( 565) = W( 63 )
  JVS( 566) = W( 65 )
  JVS( 567) = W( 66 )
  JVS( 568) = W( 67 )
  JVS( 569) = W( 68 )
  JVS( 570) = W( 70 )
  JVS( 571) = W( 72 )
  IF ( ABS(  JVS( 585 )) < TINY(a) ) THEN
         IER = 61                                      
         RETURN
  END IF
   W( 13 ) = JVS( 572 )
   W( 29 ) = JVS( 573 )
   W( 30 ) = JVS( 574 )
   W( 42 ) = JVS( 575 )
   W( 44 ) = JVS( 576 )
   W( 45 ) = JVS( 577 )
   W( 49 ) = JVS( 578 )
   W( 52 ) = JVS( 579 )
   W( 54 ) = JVS( 580 )
   W( 57 ) = JVS( 581 )
   W( 58 ) = JVS( 582 )
   W( 59 ) = JVS( 583 )
   W( 60 ) = JVS( 584 )
   W( 61 ) = JVS( 585 )
   W( 62 ) = JVS( 586 )
   W( 63 ) = JVS( 587 )
   W( 64 ) = JVS( 588 )
   W( 65 ) = JVS( 589 )
   W( 66 ) = JVS( 590 )
   W( 67 ) = JVS( 591 )
   W( 68 ) = JVS( 592 )
   W( 69 ) = JVS( 593 )
   W( 70 ) = JVS( 594 )
   W( 71 ) = JVS( 595 )
   W( 72 ) = JVS( 596 )
   W( 73 ) = JVS( 597 )
  a = -W( 13 ) / JVS(           58  )
  W( 13 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 59 )
  W( 58 ) = W( 58 ) + a*JVS( 60 )
  W( 63 ) = W( 63 ) + a*JVS( 61 )
  W( 68 ) = W( 68 ) + a*JVS( 62 )
  a = -W( 29 ) / JVS(          166  )
  W( 29 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 167 )
  W( 68 ) = W( 68 ) + a*JVS( 168 )
  W( 70 ) = W( 70 ) + a*JVS( 169 )
  a = -W( 30 ) / JVS(          170  )
  W( 30 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 171 )
  W( 68 ) = W( 68 ) + a*JVS( 172 )
  W( 70 ) = W( 70 ) + a*JVS( 173 )
  a = -W( 42 ) / JVS(          278  )
  W( 42 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 279 )
  W( 45 ) = W( 45 ) + a*JVS( 280 )
  W( 49 ) = W( 49 ) + a*JVS( 281 )
  W( 52 ) = W( 52 ) + a*JVS( 282 )
  W( 54 ) = W( 54 ) + a*JVS( 283 )
  W( 57 ) = W( 57 ) + a*JVS( 284 )
  W( 58 ) = W( 58 ) + a*JVS( 285 )
  W( 59 ) = W( 59 ) + a*JVS( 286 )
  W( 60 ) = W( 60 ) + a*JVS( 287 )
  W( 63 ) = W( 63 ) + a*JVS( 288 )
  W( 64 ) = W( 64 ) + a*JVS( 289 )
  W( 65 ) = W( 65 ) + a*JVS( 290 )
  W( 66 ) = W( 66 ) + a*JVS( 291 )
  W( 67 ) = W( 67 ) + a*JVS( 292 )
  W( 68 ) = W( 68 ) + a*JVS( 293 )
  W( 69 ) = W( 69 ) + a*JVS( 294 )
  W( 70 ) = W( 70 ) + a*JVS( 295 )
  W( 71 ) = W( 71 ) + a*JVS( 296 )
  W( 72 ) = W( 72 ) + a*JVS( 297 )
  a = -W( 44 ) / JVS(          308  )
  W( 44 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 309 )
  W( 66 ) = W( 66 ) + a*JVS( 310 )
  W( 67 ) = W( 67 ) + a*JVS( 311 )
  W( 68 ) = W( 68 ) + a*JVS( 312 )
  W( 70 ) = W( 70 ) + a*JVS( 313 )
  W( 72 ) = W( 72 ) + a*JVS( 314 )
  a = -W( 45 ) / JVS(          316  )
  W( 45 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 70 ) = W( 70 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  a = -W( 49 ) / JVS(          387  )
  W( 49 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 388 )
  W( 65 ) = W( 65 ) + a*JVS( 389 )
  W( 66 ) = W( 66 ) + a*JVS( 390 )
  W( 67 ) = W( 67 ) + a*JVS( 391 )
  W( 68 ) = W( 68 ) + a*JVS( 392 )
  W( 70 ) = W( 70 ) + a*JVS( 393 )
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  a = -W( 52 ) / JVS(          421  )
  W( 52 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 422 )
  W( 66 ) = W( 66 ) + a*JVS( 423 )
  W( 67 ) = W( 67 ) + a*JVS( 424 )
  W( 68 ) = W( 68 ) + a*JVS( 425 )
  W( 69 ) = W( 69 ) + a*JVS( 426 )
  W( 70 ) = W( 70 ) + a*JVS( 427 )
  W( 71 ) = W( 71 ) + a*JVS( 428 )
  W( 72 ) = W( 72 ) + a*JVS( 429 )
  a = -W( 54 ) / JVS(          463  )
  W( 54 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 464 )
  W( 63 ) = W( 63 ) + a*JVS( 465 )
  W( 65 ) = W( 65 ) + a*JVS( 466 )
  W( 66 ) = W( 66 ) + a*JVS( 467 )
  W( 67 ) = W( 67 ) + a*JVS( 468 )
  W( 68 ) = W( 68 ) + a*JVS( 469 )
  W( 70 ) = W( 70 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  a = -W( 57 ) / JVS(          517  )
  W( 57 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 518 )
  W( 63 ) = W( 63 ) + a*JVS( 519 )
  W( 65 ) = W( 65 ) + a*JVS( 520 )
  W( 66 ) = W( 66 ) + a*JVS( 521 )
  W( 67 ) = W( 67 ) + a*JVS( 522 )
  W( 68 ) = W( 68 ) + a*JVS( 523 )
  W( 70 ) = W( 70 ) + a*JVS( 524 )
  W( 72 ) = W( 72 ) + a*JVS( 525 )
  a = -W( 58 ) / JVS(          528  )
  W( 58 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 529 )
  W( 63 ) = W( 63 ) + a*JVS( 530 )
  W( 65 ) = W( 65 ) + a*JVS( 531 )
  W( 66 ) = W( 66 ) + a*JVS( 532 )
  W( 67 ) = W( 67 ) + a*JVS( 533 )
  W( 68 ) = W( 68 ) + a*JVS( 534 )
  W( 70 ) = W( 70 ) + a*JVS( 535 )
  W( 72 ) = W( 72 ) + a*JVS( 536 )
  a = -W( 59 ) / JVS(          545  )
  W( 59 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 546 )
  W( 62 ) = W( 62 ) + a*JVS( 547 )
  W( 63 ) = W( 63 ) + a*JVS( 548 )
  W( 65 ) = W( 65 ) + a*JVS( 549 )
  W( 66 ) = W( 66 ) + a*JVS( 550 )
  W( 67 ) = W( 67 ) + a*JVS( 551 )
  W( 68 ) = W( 68 ) + a*JVS( 552 )
  W( 70 ) = W( 70 ) + a*JVS( 553 )
  W( 72 ) = W( 72 ) + a*JVS( 554 )
  a = -W( 60 ) / JVS(          563  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 564 )
  W( 63 ) = W( 63 ) + a*JVS( 565 )
  W( 65 ) = W( 65 ) + a*JVS( 566 )
  W( 66 ) = W( 66 ) + a*JVS( 567 )
  W( 67 ) = W( 67 ) + a*JVS( 568 )
  W( 68 ) = W( 68 ) + a*JVS( 569 )
  W( 70 ) = W( 70 ) + a*JVS( 570 )
  W( 72 ) = W( 72 ) + a*JVS( 571 )
  JVS( 572) = W( 13 )
  JVS( 573) = W( 29 )
  JVS( 574) = W( 30 )
  JVS( 575) = W( 42 )
  JVS( 576) = W( 44 )
  JVS( 577) = W( 45 )
  JVS( 578) = W( 49 )
  JVS( 579) = W( 52 )
  JVS( 580) = W( 54 )
  JVS( 581) = W( 57 )
  JVS( 582) = W( 58 )
  JVS( 583) = W( 59 )
  JVS( 584) = W( 60 )
  JVS( 585) = W( 61 )
  JVS( 586) = W( 62 )
  JVS( 587) = W( 63 )
  JVS( 588) = W( 64 )
  JVS( 589) = W( 65 )
  JVS( 590) = W( 66 )
  JVS( 591) = W( 67 )
  JVS( 592) = W( 68 )
  JVS( 593) = W( 69 )
  JVS( 594) = W( 70 )
  JVS( 595) = W( 71 )
  JVS( 596) = W( 72 )
  JVS( 597) = W( 73 )
  IF ( ABS(  JVS( 602 )) < TINY(a) ) THEN
         IER = 62                                      
         RETURN
  END IF
   W( 6 ) = JVS( 598 )
   W( 33 ) = JVS( 599 )
   W( 37 ) = JVS( 600 )
   W( 51 ) = JVS( 601 )
   W( 62 ) = JVS( 602 )
   W( 63 ) = JVS( 603 )
   W( 65 ) = JVS( 604 )
   W( 66 ) = JVS( 605 )
   W( 67 ) = JVS( 606 )
   W( 68 ) = JVS( 607 )
   W( 70 ) = JVS( 608 )
   W( 72 ) = JVS( 609 )
   W( 73 ) = JVS( 610 )
  a = -W( 6 ) / JVS(           43  )
  W( 6 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 44 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 51 ) / JVS(          411  )
  W( 51 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 412 )
  W( 63 ) = W( 63 ) + a*JVS( 413 )
  W( 65 ) = W( 65 ) + a*JVS( 414 )
  W( 66 ) = W( 66 ) + a*JVS( 415 )
  W( 67 ) = W( 67 ) + a*JVS( 416 )
  W( 68 ) = W( 68 ) + a*JVS( 417 )
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 72 ) = W( 72 ) + a*JVS( 419 )
  JVS( 598) = W( 6 )
  JVS( 599) = W( 33 )
  JVS( 600) = W( 37 )
  JVS( 601) = W( 51 )
  JVS( 602) = W( 62 )
  JVS( 603) = W( 63 )
  JVS( 604) = W( 65 )
  JVS( 605) = W( 66 )
  JVS( 606) = W( 67 )
  JVS( 607) = W( 68 )
  JVS( 608) = W( 70 )
  JVS( 609) = W( 72 )
  JVS( 610) = W( 73 )
  IF ( ABS(  JVS( 628 )) < TINY(a) ) THEN
         IER = 63                                      
         RETURN
  END IF
   W( 6 ) = JVS( 611 )
   W( 21 ) = JVS( 612 )
   W( 22 ) = JVS( 613 )
   W( 24 ) = JVS( 614 )
   W( 25 ) = JVS( 615 )
   W( 29 ) = JVS( 616 )
   W( 30 ) = JVS( 617 )
   W( 32 ) = JVS( 618 )
   W( 33 ) = JVS( 619 )
   W( 35 ) = JVS( 620 )
   W( 37 ) = JVS( 621 )
   W( 51 ) = JVS( 622 )
   W( 54 ) = JVS( 623 )
   W( 55 ) = JVS( 624 )
   W( 57 ) = JVS( 625 )
   W( 58 ) = JVS( 626 )
   W( 62 ) = JVS( 627 )
   W( 63 ) = JVS( 628 )
   W( 65 ) = JVS( 629 )
   W( 66 ) = JVS( 630 )
   W( 67 ) = JVS( 631 )
   W( 68 ) = JVS( 632 )
   W( 70 ) = JVS( 633 )
   W( 72 ) = JVS( 634 )
   W( 73 ) = JVS( 635 )
  a = -W( 6 ) / JVS(           43  )
  W( 6 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 44 )
  a = -W( 21 ) / JVS(          103  )
  W( 21 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 104 )
  W( 68 ) = W( 68 ) + a*JVS( 105 )
  W( 73 ) = W( 73 ) + a*JVS( 106 )
  a = -W( 22 ) / JVS(          108  )
  W( 22 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 109 )
  W( 68 ) = W( 68 ) + a*JVS( 110 )
  W( 73 ) = W( 73 ) + a*JVS( 111 )
  a = -W( 24 ) / JVS(          119  )
  W( 24 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 120 )
  W( 68 ) = W( 68 ) + a*JVS( 121 )
  W( 70 ) = W( 70 ) + a*JVS( 122 )
  a = -W( 25 ) / JVS(          123  )
  W( 25 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 124 )
  W( 63 ) = W( 63 ) + a*JVS( 125 )
  W( 68 ) = W( 68 ) + a*JVS( 126 )
  W( 73 ) = W( 73 ) + a*JVS( 127 )
  a = -W( 29 ) / JVS(          166  )
  W( 29 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 167 )
  W( 68 ) = W( 68 ) + a*JVS( 168 )
  W( 70 ) = W( 70 ) + a*JVS( 169 )
  a = -W( 30 ) / JVS(          170  )
  W( 30 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 171 )
  W( 68 ) = W( 68 ) + a*JVS( 172 )
  W( 70 ) = W( 70 ) + a*JVS( 173 )
  a = -W( 32 ) / JVS(          185  )
  W( 32 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 186 )
  W( 67 ) = W( 67 ) + a*JVS( 187 )
  W( 68 ) = W( 68 ) + a*JVS( 188 )
  W( 70 ) = W( 70 ) + a*JVS( 189 )
  W( 73 ) = W( 73 ) + a*JVS( 190 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 35 ) / JVS(          208  )
  W( 35 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 209 )
  W( 63 ) = W( 63 ) + a*JVS( 210 )
  W( 68 ) = W( 68 ) + a*JVS( 211 )
  W( 70 ) = W( 70 ) + a*JVS( 212 )
  W( 73 ) = W( 73 ) + a*JVS( 213 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 51 ) / JVS(          411  )
  W( 51 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 412 )
  W( 63 ) = W( 63 ) + a*JVS( 413 )
  W( 65 ) = W( 65 ) + a*JVS( 414 )
  W( 66 ) = W( 66 ) + a*JVS( 415 )
  W( 67 ) = W( 67 ) + a*JVS( 416 )
  W( 68 ) = W( 68 ) + a*JVS( 417 )
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 72 ) = W( 72 ) + a*JVS( 419 )
  a = -W( 54 ) / JVS(          463  )
  W( 54 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 464 )
  W( 63 ) = W( 63 ) + a*JVS( 465 )
  W( 65 ) = W( 65 ) + a*JVS( 466 )
  W( 66 ) = W( 66 ) + a*JVS( 467 )
  W( 67 ) = W( 67 ) + a*JVS( 468 )
  W( 68 ) = W( 68 ) + a*JVS( 469 )
  W( 70 ) = W( 70 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  a = -W( 55 ) / JVS(          476  )
  W( 55 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 477 )
  W( 63 ) = W( 63 ) + a*JVS( 478 )
  W( 65 ) = W( 65 ) + a*JVS( 479 )
  W( 66 ) = W( 66 ) + a*JVS( 480 )
  W( 67 ) = W( 67 ) + a*JVS( 481 )
  W( 68 ) = W( 68 ) + a*JVS( 482 )
  W( 70 ) = W( 70 ) + a*JVS( 483 )
  W( 72 ) = W( 72 ) + a*JVS( 484 )
  W( 73 ) = W( 73 ) + a*JVS( 485 )
  a = -W( 57 ) / JVS(          517  )
  W( 57 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 518 )
  W( 63 ) = W( 63 ) + a*JVS( 519 )
  W( 65 ) = W( 65 ) + a*JVS( 520 )
  W( 66 ) = W( 66 ) + a*JVS( 521 )
  W( 67 ) = W( 67 ) + a*JVS( 522 )
  W( 68 ) = W( 68 ) + a*JVS( 523 )
  W( 70 ) = W( 70 ) + a*JVS( 524 )
  W( 72 ) = W( 72 ) + a*JVS( 525 )
  a = -W( 58 ) / JVS(          528  )
  W( 58 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 529 )
  W( 63 ) = W( 63 ) + a*JVS( 530 )
  W( 65 ) = W( 65 ) + a*JVS( 531 )
  W( 66 ) = W( 66 ) + a*JVS( 532 )
  W( 67 ) = W( 67 ) + a*JVS( 533 )
  W( 68 ) = W( 68 ) + a*JVS( 534 )
  W( 70 ) = W( 70 ) + a*JVS( 535 )
  W( 72 ) = W( 72 ) + a*JVS( 536 )
  a = -W( 62 ) / JVS(          602  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 603 )
  W( 65 ) = W( 65 ) + a*JVS( 604 )
  W( 66 ) = W( 66 ) + a*JVS( 605 )
  W( 67 ) = W( 67 ) + a*JVS( 606 )
  W( 68 ) = W( 68 ) + a*JVS( 607 )
  W( 70 ) = W( 70 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  JVS( 611) = W( 6 )
  JVS( 612) = W( 21 )
  JVS( 613) = W( 22 )
  JVS( 614) = W( 24 )
  JVS( 615) = W( 25 )
  JVS( 616) = W( 29 )
  JVS( 617) = W( 30 )
  JVS( 618) = W( 32 )
  JVS( 619) = W( 33 )
  JVS( 620) = W( 35 )
  JVS( 621) = W( 37 )
  JVS( 622) = W( 51 )
  JVS( 623) = W( 54 )
  JVS( 624) = W( 55 )
  JVS( 625) = W( 57 )
  JVS( 626) = W( 58 )
  JVS( 627) = W( 62 )
  JVS( 628) = W( 63 )
  JVS( 629) = W( 65 )
  JVS( 630) = W( 66 )
  JVS( 631) = W( 67 )
  JVS( 632) = W( 68 )
  JVS( 633) = W( 70 )
  JVS( 634) = W( 72 )
  JVS( 635) = W( 73 )
  IF ( ABS(  JVS( 652 )) < TINY(a) ) THEN
         IER = 64                                      
         RETURN
  END IF
   W( 29 ) = JVS( 636 )
   W( 30 ) = JVS( 637 )
   W( 33 ) = JVS( 638 )
   W( 37 ) = JVS( 639 )
   W( 42 ) = JVS( 640 )
   W( 44 ) = JVS( 641 )
   W( 45 ) = JVS( 642 )
   W( 49 ) = JVS( 643 )
   W( 52 ) = JVS( 644 )
   W( 54 ) = JVS( 645 )
   W( 57 ) = JVS( 646 )
   W( 58 ) = JVS( 647 )
   W( 59 ) = JVS( 648 )
   W( 60 ) = JVS( 649 )
   W( 62 ) = JVS( 650 )
   W( 63 ) = JVS( 651 )
   W( 64 ) = JVS( 652 )
   W( 65 ) = JVS( 653 )
   W( 66 ) = JVS( 654 )
   W( 67 ) = JVS( 655 )
   W( 68 ) = JVS( 656 )
   W( 69 ) = JVS( 657 )
   W( 70 ) = JVS( 658 )
   W( 71 ) = JVS( 659 )
   W( 72 ) = JVS( 660 )
   W( 73 ) = JVS( 661 )
  a = -W( 29 ) / JVS(          166  )
  W( 29 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 167 )
  W( 68 ) = W( 68 ) + a*JVS( 168 )
  W( 70 ) = W( 70 ) + a*JVS( 169 )
  a = -W( 30 ) / JVS(          170  )
  W( 30 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 171 )
  W( 68 ) = W( 68 ) + a*JVS( 172 )
  W( 70 ) = W( 70 ) + a*JVS( 173 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 42 ) / JVS(          278  )
  W( 42 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 279 )
  W( 45 ) = W( 45 ) + a*JVS( 280 )
  W( 49 ) = W( 49 ) + a*JVS( 281 )
  W( 52 ) = W( 52 ) + a*JVS( 282 )
  W( 54 ) = W( 54 ) + a*JVS( 283 )
  W( 57 ) = W( 57 ) + a*JVS( 284 )
  W( 58 ) = W( 58 ) + a*JVS( 285 )
  W( 59 ) = W( 59 ) + a*JVS( 286 )
  W( 60 ) = W( 60 ) + a*JVS( 287 )
  W( 63 ) = W( 63 ) + a*JVS( 288 )
  W( 64 ) = W( 64 ) + a*JVS( 289 )
  W( 65 ) = W( 65 ) + a*JVS( 290 )
  W( 66 ) = W( 66 ) + a*JVS( 291 )
  W( 67 ) = W( 67 ) + a*JVS( 292 )
  W( 68 ) = W( 68 ) + a*JVS( 293 )
  W( 69 ) = W( 69 ) + a*JVS( 294 )
  W( 70 ) = W( 70 ) + a*JVS( 295 )
  W( 71 ) = W( 71 ) + a*JVS( 296 )
  W( 72 ) = W( 72 ) + a*JVS( 297 )
  a = -W( 44 ) / JVS(          308  )
  W( 44 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 309 )
  W( 66 ) = W( 66 ) + a*JVS( 310 )
  W( 67 ) = W( 67 ) + a*JVS( 311 )
  W( 68 ) = W( 68 ) + a*JVS( 312 )
  W( 70 ) = W( 70 ) + a*JVS( 313 )
  W( 72 ) = W( 72 ) + a*JVS( 314 )
  a = -W( 45 ) / JVS(          316  )
  W( 45 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 70 ) = W( 70 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  a = -W( 49 ) / JVS(          387  )
  W( 49 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 388 )
  W( 65 ) = W( 65 ) + a*JVS( 389 )
  W( 66 ) = W( 66 ) + a*JVS( 390 )
  W( 67 ) = W( 67 ) + a*JVS( 391 )
  W( 68 ) = W( 68 ) + a*JVS( 392 )
  W( 70 ) = W( 70 ) + a*JVS( 393 )
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  a = -W( 52 ) / JVS(          421  )
  W( 52 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 422 )
  W( 66 ) = W( 66 ) + a*JVS( 423 )
  W( 67 ) = W( 67 ) + a*JVS( 424 )
  W( 68 ) = W( 68 ) + a*JVS( 425 )
  W( 69 ) = W( 69 ) + a*JVS( 426 )
  W( 70 ) = W( 70 ) + a*JVS( 427 )
  W( 71 ) = W( 71 ) + a*JVS( 428 )
  W( 72 ) = W( 72 ) + a*JVS( 429 )
  a = -W( 54 ) / JVS(          463  )
  W( 54 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 464 )
  W( 63 ) = W( 63 ) + a*JVS( 465 )
  W( 65 ) = W( 65 ) + a*JVS( 466 )
  W( 66 ) = W( 66 ) + a*JVS( 467 )
  W( 67 ) = W( 67 ) + a*JVS( 468 )
  W( 68 ) = W( 68 ) + a*JVS( 469 )
  W( 70 ) = W( 70 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  a = -W( 57 ) / JVS(          517  )
  W( 57 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 518 )
  W( 63 ) = W( 63 ) + a*JVS( 519 )
  W( 65 ) = W( 65 ) + a*JVS( 520 )
  W( 66 ) = W( 66 ) + a*JVS( 521 )
  W( 67 ) = W( 67 ) + a*JVS( 522 )
  W( 68 ) = W( 68 ) + a*JVS( 523 )
  W( 70 ) = W( 70 ) + a*JVS( 524 )
  W( 72 ) = W( 72 ) + a*JVS( 525 )
  a = -W( 58 ) / JVS(          528  )
  W( 58 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 529 )
  W( 63 ) = W( 63 ) + a*JVS( 530 )
  W( 65 ) = W( 65 ) + a*JVS( 531 )
  W( 66 ) = W( 66 ) + a*JVS( 532 )
  W( 67 ) = W( 67 ) + a*JVS( 533 )
  W( 68 ) = W( 68 ) + a*JVS( 534 )
  W( 70 ) = W( 70 ) + a*JVS( 535 )
  W( 72 ) = W( 72 ) + a*JVS( 536 )
  a = -W( 59 ) / JVS(          545  )
  W( 59 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 546 )
  W( 62 ) = W( 62 ) + a*JVS( 547 )
  W( 63 ) = W( 63 ) + a*JVS( 548 )
  W( 65 ) = W( 65 ) + a*JVS( 549 )
  W( 66 ) = W( 66 ) + a*JVS( 550 )
  W( 67 ) = W( 67 ) + a*JVS( 551 )
  W( 68 ) = W( 68 ) + a*JVS( 552 )
  W( 70 ) = W( 70 ) + a*JVS( 553 )
  W( 72 ) = W( 72 ) + a*JVS( 554 )
  a = -W( 60 ) / JVS(          563  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 564 )
  W( 63 ) = W( 63 ) + a*JVS( 565 )
  W( 65 ) = W( 65 ) + a*JVS( 566 )
  W( 66 ) = W( 66 ) + a*JVS( 567 )
  W( 67 ) = W( 67 ) + a*JVS( 568 )
  W( 68 ) = W( 68 ) + a*JVS( 569 )
  W( 70 ) = W( 70 ) + a*JVS( 570 )
  W( 72 ) = W( 72 ) + a*JVS( 571 )
  a = -W( 62 ) / JVS(          602  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 603 )
  W( 65 ) = W( 65 ) + a*JVS( 604 )
  W( 66 ) = W( 66 ) + a*JVS( 605 )
  W( 67 ) = W( 67 ) + a*JVS( 606 )
  W( 68 ) = W( 68 ) + a*JVS( 607 )
  W( 70 ) = W( 70 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  a = -W( 63 ) / JVS(          628  )
  W( 63 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 629 )
  W( 66 ) = W( 66 ) + a*JVS( 630 )
  W( 67 ) = W( 67 ) + a*JVS( 631 )
  W( 68 ) = W( 68 ) + a*JVS( 632 )
  W( 70 ) = W( 70 ) + a*JVS( 633 )
  W( 72 ) = W( 72 ) + a*JVS( 634 )
  W( 73 ) = W( 73 ) + a*JVS( 635 )
  JVS( 636) = W( 29 )
  JVS( 637) = W( 30 )
  JVS( 638) = W( 33 )
  JVS( 639) = W( 37 )
  JVS( 640) = W( 42 )
  JVS( 641) = W( 44 )
  JVS( 642) = W( 45 )
  JVS( 643) = W( 49 )
  JVS( 644) = W( 52 )
  JVS( 645) = W( 54 )
  JVS( 646) = W( 57 )
  JVS( 647) = W( 58 )
  JVS( 648) = W( 59 )
  JVS( 649) = W( 60 )
  JVS( 650) = W( 62 )
  JVS( 651) = W( 63 )
  JVS( 652) = W( 64 )
  JVS( 653) = W( 65 )
  JVS( 654) = W( 66 )
  JVS( 655) = W( 67 )
  JVS( 656) = W( 68 )
  JVS( 657) = W( 69 )
  JVS( 658) = W( 70 )
  JVS( 659) = W( 71 )
  JVS( 660) = W( 72 )
  JVS( 661) = W( 73 )
  IF ( ABS(  JVS( 691 )) < TINY(a) ) THEN
         IER = 65                                      
         RETURN
  END IF
   W( 14 ) = JVS( 662 )
   W( 17 ) = JVS( 663 )
   W( 26 ) = JVS( 664 )
   W( 33 ) = JVS( 665 )
   W( 36 ) = JVS( 666 )
   W( 37 ) = JVS( 667 )
   W( 38 ) = JVS( 668 )
   W( 39 ) = JVS( 669 )
   W( 41 ) = JVS( 670 )
   W( 43 ) = JVS( 671 )
   W( 44 ) = JVS( 672 )
   W( 45 ) = JVS( 673 )
   W( 46 ) = JVS( 674 )
   W( 47 ) = JVS( 675 )
   W( 49 ) = JVS( 676 )
   W( 50 ) = JVS( 677 )
   W( 52 ) = JVS( 678 )
   W( 53 ) = JVS( 679 )
   W( 54 ) = JVS( 680 )
   W( 55 ) = JVS( 681 )
   W( 56 ) = JVS( 682 )
   W( 57 ) = JVS( 683 )
   W( 58 ) = JVS( 684 )
   W( 59 ) = JVS( 685 )
   W( 60 ) = JVS( 686 )
   W( 61 ) = JVS( 687 )
   W( 62 ) = JVS( 688 )
   W( 63 ) = JVS( 689 )
   W( 64 ) = JVS( 690 )
   W( 65 ) = JVS( 691 )
   W( 66 ) = JVS( 692 )
   W( 67 ) = JVS( 693 )
   W( 68 ) = JVS( 694 )
   W( 69 ) = JVS( 695 )
   W( 70 ) = JVS( 696 )
   W( 71 ) = JVS( 697 )
   W( 72 ) = JVS( 698 )
   W( 73 ) = JVS( 699 )
  a = -W( 14 ) / JVS(           63  )
  W( 14 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 64 )
  W( 58 ) = W( 58 ) + a*JVS( 65 )
  W( 63 ) = W( 63 ) + a*JVS( 66 )
  W( 68 ) = W( 68 ) + a*JVS( 67 )
  a = -W( 17 ) / JVS(           75  )
  W( 17 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 76 )
  W( 67 ) = W( 67 ) + a*JVS( 77 )
  W( 68 ) = W( 68 ) + a*JVS( 78 )
  a = -W( 26 ) / JVS(          128  )
  W( 26 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 129 )
  W( 63 ) = W( 63 ) + a*JVS( 130 )
  W( 66 ) = W( 66 ) + a*JVS( 131 )
  W( 67 ) = W( 67 ) + a*JVS( 132 )
  W( 68 ) = W( 68 ) + a*JVS( 133 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 36 ) / JVS(          215  )
  W( 36 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 216 )
  W( 65 ) = W( 65 ) + a*JVS( 217 )
  W( 66 ) = W( 66 ) + a*JVS( 218 )
  W( 67 ) = W( 67 ) + a*JVS( 219 )
  W( 68 ) = W( 68 ) + a*JVS( 220 )
  W( 70 ) = W( 70 ) + a*JVS( 221 )
  W( 72 ) = W( 72 ) + a*JVS( 222 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 38 ) / JVS(          228  )
  W( 38 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 229 )
  W( 65 ) = W( 65 ) + a*JVS( 230 )
  W( 66 ) = W( 66 ) + a*JVS( 231 )
  W( 67 ) = W( 67 ) + a*JVS( 232 )
  W( 68 ) = W( 68 ) + a*JVS( 233 )
  W( 70 ) = W( 70 ) + a*JVS( 234 )
  W( 72 ) = W( 72 ) + a*JVS( 235 )
  a = -W( 39 ) / JVS(          236  )
  W( 39 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 237 )
  W( 65 ) = W( 65 ) + a*JVS( 238 )
  W( 66 ) = W( 66 ) + a*JVS( 239 )
  W( 67 ) = W( 67 ) + a*JVS( 240 )
  W( 68 ) = W( 68 ) + a*JVS( 241 )
  W( 70 ) = W( 70 ) + a*JVS( 242 )
  W( 72 ) = W( 72 ) + a*JVS( 243 )
  a = -W( 41 ) / JVS(          264  )
  W( 41 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 265 )
  W( 65 ) = W( 65 ) + a*JVS( 266 )
  W( 66 ) = W( 66 ) + a*JVS( 267 )
  W( 67 ) = W( 67 ) + a*JVS( 268 )
  W( 68 ) = W( 68 ) + a*JVS( 269 )
  W( 70 ) = W( 70 ) + a*JVS( 270 )
  W( 72 ) = W( 72 ) + a*JVS( 271 )
  W( 73 ) = W( 73 ) + a*JVS( 272 )
  a = -W( 43 ) / JVS(          299  )
  W( 43 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 300 )
  W( 65 ) = W( 65 ) + a*JVS( 301 )
  W( 66 ) = W( 66 ) + a*JVS( 302 )
  W( 67 ) = W( 67 ) + a*JVS( 303 )
  W( 68 ) = W( 68 ) + a*JVS( 304 )
  W( 70 ) = W( 70 ) + a*JVS( 305 )
  W( 72 ) = W( 72 ) + a*JVS( 306 )
  a = -W( 44 ) / JVS(          308  )
  W( 44 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 309 )
  W( 66 ) = W( 66 ) + a*JVS( 310 )
  W( 67 ) = W( 67 ) + a*JVS( 311 )
  W( 68 ) = W( 68 ) + a*JVS( 312 )
  W( 70 ) = W( 70 ) + a*JVS( 313 )
  W( 72 ) = W( 72 ) + a*JVS( 314 )
  a = -W( 45 ) / JVS(          316  )
  W( 45 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 70 ) = W( 70 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 325 )
  W( 65 ) = W( 65 ) + a*JVS( 326 )
  W( 66 ) = W( 66 ) + a*JVS( 327 )
  W( 67 ) = W( 67 ) + a*JVS( 328 )
  W( 68 ) = W( 68 ) + a*JVS( 329 )
  W( 70 ) = W( 70 ) + a*JVS( 330 )
  W( 72 ) = W( 72 ) + a*JVS( 331 )
  W( 73 ) = W( 73 ) + a*JVS( 332 )
  a = -W( 47 ) / JVS(          334  )
  W( 47 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 335 )
  W( 65 ) = W( 65 ) + a*JVS( 336 )
  W( 66 ) = W( 66 ) + a*JVS( 337 )
  W( 67 ) = W( 67 ) + a*JVS( 338 )
  W( 68 ) = W( 68 ) + a*JVS( 339 )
  W( 70 ) = W( 70 ) + a*JVS( 340 )
  W( 72 ) = W( 72 ) + a*JVS( 341 )
  W( 73 ) = W( 73 ) + a*JVS( 342 )
  a = -W( 49 ) / JVS(          387  )
  W( 49 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 388 )
  W( 65 ) = W( 65 ) + a*JVS( 389 )
  W( 66 ) = W( 66 ) + a*JVS( 390 )
  W( 67 ) = W( 67 ) + a*JVS( 391 )
  W( 68 ) = W( 68 ) + a*JVS( 392 )
  W( 70 ) = W( 70 ) + a*JVS( 393 )
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  a = -W( 50 ) / JVS(          397  )
  W( 50 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 398 )
  W( 63 ) = W( 63 ) + a*JVS( 399 )
  W( 65 ) = W( 65 ) + a*JVS( 400 )
  W( 66 ) = W( 66 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 70 ) = W( 70 ) + a*JVS( 404 )
  W( 72 ) = W( 72 ) + a*JVS( 405 )
  a = -W( 52 ) / JVS(          421  )
  W( 52 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 422 )
  W( 66 ) = W( 66 ) + a*JVS( 423 )
  W( 67 ) = W( 67 ) + a*JVS( 424 )
  W( 68 ) = W( 68 ) + a*JVS( 425 )
  W( 69 ) = W( 69 ) + a*JVS( 426 )
  W( 70 ) = W( 70 ) + a*JVS( 427 )
  W( 71 ) = W( 71 ) + a*JVS( 428 )
  W( 72 ) = W( 72 ) + a*JVS( 429 )
  a = -W( 53 ) / JVS(          442  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 443 )
  W( 57 ) = W( 57 ) + a*JVS( 444 )
  W( 58 ) = W( 58 ) + a*JVS( 445 )
  W( 59 ) = W( 59 ) + a*JVS( 446 )
  W( 60 ) = W( 60 ) + a*JVS( 447 )
  W( 61 ) = W( 61 ) + a*JVS( 448 )
  W( 62 ) = W( 62 ) + a*JVS( 449 )
  W( 63 ) = W( 63 ) + a*JVS( 450 )
  W( 64 ) = W( 64 ) + a*JVS( 451 )
  W( 65 ) = W( 65 ) + a*JVS( 452 )
  W( 66 ) = W( 66 ) + a*JVS( 453 )
  W( 67 ) = W( 67 ) + a*JVS( 454 )
  W( 68 ) = W( 68 ) + a*JVS( 455 )
  W( 69 ) = W( 69 ) + a*JVS( 456 )
  W( 70 ) = W( 70 ) + a*JVS( 457 )
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 72 ) = W( 72 ) + a*JVS( 459 )
  a = -W( 54 ) / JVS(          463  )
  W( 54 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 464 )
  W( 63 ) = W( 63 ) + a*JVS( 465 )
  W( 65 ) = W( 65 ) + a*JVS( 466 )
  W( 66 ) = W( 66 ) + a*JVS( 467 )
  W( 67 ) = W( 67 ) + a*JVS( 468 )
  W( 68 ) = W( 68 ) + a*JVS( 469 )
  W( 70 ) = W( 70 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  a = -W( 55 ) / JVS(          476  )
  W( 55 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 477 )
  W( 63 ) = W( 63 ) + a*JVS( 478 )
  W( 65 ) = W( 65 ) + a*JVS( 479 )
  W( 66 ) = W( 66 ) + a*JVS( 480 )
  W( 67 ) = W( 67 ) + a*JVS( 481 )
  W( 68 ) = W( 68 ) + a*JVS( 482 )
  W( 70 ) = W( 70 ) + a*JVS( 483 )
  W( 72 ) = W( 72 ) + a*JVS( 484 )
  W( 73 ) = W( 73 ) + a*JVS( 485 )
  a = -W( 56 ) / JVS(          500  )
  W( 56 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 501 )
  W( 63 ) = W( 63 ) + a*JVS( 502 )
  W( 64 ) = W( 64 ) + a*JVS( 503 )
  W( 65 ) = W( 65 ) + a*JVS( 504 )
  W( 66 ) = W( 66 ) + a*JVS( 505 )
  W( 67 ) = W( 67 ) + a*JVS( 506 )
  W( 68 ) = W( 68 ) + a*JVS( 507 )
  W( 69 ) = W( 69 ) + a*JVS( 508 )
  W( 70 ) = W( 70 ) + a*JVS( 509 )
  W( 71 ) = W( 71 ) + a*JVS( 510 )
  W( 72 ) = W( 72 ) + a*JVS( 511 )
  W( 73 ) = W( 73 ) + a*JVS( 512 )
  a = -W( 57 ) / JVS(          517  )
  W( 57 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 518 )
  W( 63 ) = W( 63 ) + a*JVS( 519 )
  W( 65 ) = W( 65 ) + a*JVS( 520 )
  W( 66 ) = W( 66 ) + a*JVS( 521 )
  W( 67 ) = W( 67 ) + a*JVS( 522 )
  W( 68 ) = W( 68 ) + a*JVS( 523 )
  W( 70 ) = W( 70 ) + a*JVS( 524 )
  W( 72 ) = W( 72 ) + a*JVS( 525 )
  a = -W( 58 ) / JVS(          528  )
  W( 58 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 529 )
  W( 63 ) = W( 63 ) + a*JVS( 530 )
  W( 65 ) = W( 65 ) + a*JVS( 531 )
  W( 66 ) = W( 66 ) + a*JVS( 532 )
  W( 67 ) = W( 67 ) + a*JVS( 533 )
  W( 68 ) = W( 68 ) + a*JVS( 534 )
  W( 70 ) = W( 70 ) + a*JVS( 535 )
  W( 72 ) = W( 72 ) + a*JVS( 536 )
  a = -W( 59 ) / JVS(          545  )
  W( 59 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 546 )
  W( 62 ) = W( 62 ) + a*JVS( 547 )
  W( 63 ) = W( 63 ) + a*JVS( 548 )
  W( 65 ) = W( 65 ) + a*JVS( 549 )
  W( 66 ) = W( 66 ) + a*JVS( 550 )
  W( 67 ) = W( 67 ) + a*JVS( 551 )
  W( 68 ) = W( 68 ) + a*JVS( 552 )
  W( 70 ) = W( 70 ) + a*JVS( 553 )
  W( 72 ) = W( 72 ) + a*JVS( 554 )
  a = -W( 60 ) / JVS(          563  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 564 )
  W( 63 ) = W( 63 ) + a*JVS( 565 )
  W( 65 ) = W( 65 ) + a*JVS( 566 )
  W( 66 ) = W( 66 ) + a*JVS( 567 )
  W( 67 ) = W( 67 ) + a*JVS( 568 )
  W( 68 ) = W( 68 ) + a*JVS( 569 )
  W( 70 ) = W( 70 ) + a*JVS( 570 )
  W( 72 ) = W( 72 ) + a*JVS( 571 )
  a = -W( 61 ) / JVS(          585  )
  W( 61 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 586 )
  W( 63 ) = W( 63 ) + a*JVS( 587 )
  W( 64 ) = W( 64 ) + a*JVS( 588 )
  W( 65 ) = W( 65 ) + a*JVS( 589 )
  W( 66 ) = W( 66 ) + a*JVS( 590 )
  W( 67 ) = W( 67 ) + a*JVS( 591 )
  W( 68 ) = W( 68 ) + a*JVS( 592 )
  W( 69 ) = W( 69 ) + a*JVS( 593 )
  W( 70 ) = W( 70 ) + a*JVS( 594 )
  W( 71 ) = W( 71 ) + a*JVS( 595 )
  W( 72 ) = W( 72 ) + a*JVS( 596 )
  W( 73 ) = W( 73 ) + a*JVS( 597 )
  a = -W( 62 ) / JVS(          602  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 603 )
  W( 65 ) = W( 65 ) + a*JVS( 604 )
  W( 66 ) = W( 66 ) + a*JVS( 605 )
  W( 67 ) = W( 67 ) + a*JVS( 606 )
  W( 68 ) = W( 68 ) + a*JVS( 607 )
  W( 70 ) = W( 70 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  a = -W( 63 ) / JVS(          628  )
  W( 63 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 629 )
  W( 66 ) = W( 66 ) + a*JVS( 630 )
  W( 67 ) = W( 67 ) + a*JVS( 631 )
  W( 68 ) = W( 68 ) + a*JVS( 632 )
  W( 70 ) = W( 70 ) + a*JVS( 633 )
  W( 72 ) = W( 72 ) + a*JVS( 634 )
  W( 73 ) = W( 73 ) + a*JVS( 635 )
  a = -W( 64 ) / JVS(          652  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 653 )
  W( 66 ) = W( 66 ) + a*JVS( 654 )
  W( 67 ) = W( 67 ) + a*JVS( 655 )
  W( 68 ) = W( 68 ) + a*JVS( 656 )
  W( 69 ) = W( 69 ) + a*JVS( 657 )
  W( 70 ) = W( 70 ) + a*JVS( 658 )
  W( 71 ) = W( 71 ) + a*JVS( 659 )
  W( 72 ) = W( 72 ) + a*JVS( 660 )
  W( 73 ) = W( 73 ) + a*JVS( 661 )
  JVS( 662) = W( 14 )
  JVS( 663) = W( 17 )
  JVS( 664) = W( 26 )
  JVS( 665) = W( 33 )
  JVS( 666) = W( 36 )
  JVS( 667) = W( 37 )
  JVS( 668) = W( 38 )
  JVS( 669) = W( 39 )
  JVS( 670) = W( 41 )
  JVS( 671) = W( 43 )
  JVS( 672) = W( 44 )
  JVS( 673) = W( 45 )
  JVS( 674) = W( 46 )
  JVS( 675) = W( 47 )
  JVS( 676) = W( 49 )
  JVS( 677) = W( 50 )
  JVS( 678) = W( 52 )
  JVS( 679) = W( 53 )
  JVS( 680) = W( 54 )
  JVS( 681) = W( 55 )
  JVS( 682) = W( 56 )
  JVS( 683) = W( 57 )
  JVS( 684) = W( 58 )
  JVS( 685) = W( 59 )
  JVS( 686) = W( 60 )
  JVS( 687) = W( 61 )
  JVS( 688) = W( 62 )
  JVS( 689) = W( 63 )
  JVS( 690) = W( 64 )
  JVS( 691) = W( 65 )
  JVS( 692) = W( 66 )
  JVS( 693) = W( 67 )
  JVS( 694) = W( 68 )
  JVS( 695) = W( 69 )
  JVS( 696) = W( 70 )
  JVS( 697) = W( 71 )
  JVS( 698) = W( 72 )
  JVS( 699) = W( 73 )
  IF ( ABS(  JVS( 734 )) < TINY(a) ) THEN
         IER = 66                                      
         RETURN
  END IF
   W( 23 ) = JVS( 700 )
   W( 26 ) = JVS( 701 )
   W( 31 ) = JVS( 702 )
   W( 33 ) = JVS( 703 )
   W( 35 ) = JVS( 704 )
   W( 36 ) = JVS( 705 )
   W( 37 ) = JVS( 706 )
   W( 38 ) = JVS( 707 )
   W( 39 ) = JVS( 708 )
   W( 40 ) = JVS( 709 )
   W( 41 ) = JVS( 710 )
   W( 42 ) = JVS( 711 )
   W( 43 ) = JVS( 712 )
   W( 44 ) = JVS( 713 )
   W( 45 ) = JVS( 714 )
   W( 46 ) = JVS( 715 )
   W( 47 ) = JVS( 716 )
   W( 49 ) = JVS( 717 )
   W( 50 ) = JVS( 718 )
   W( 51 ) = JVS( 719 )
   W( 52 ) = JVS( 720 )
   W( 53 ) = JVS( 721 )
   W( 54 ) = JVS( 722 )
   W( 55 ) = JVS( 723 )
   W( 56 ) = JVS( 724 )
   W( 57 ) = JVS( 725 )
   W( 58 ) = JVS( 726 )
   W( 59 ) = JVS( 727 )
   W( 60 ) = JVS( 728 )
   W( 61 ) = JVS( 729 )
   W( 62 ) = JVS( 730 )
   W( 63 ) = JVS( 731 )
   W( 64 ) = JVS( 732 )
   W( 65 ) = JVS( 733 )
   W( 66 ) = JVS( 734 )
   W( 67 ) = JVS( 735 )
   W( 68 ) = JVS( 736 )
   W( 69 ) = JVS( 737 )
   W( 70 ) = JVS( 738 )
   W( 71 ) = JVS( 739 )
   W( 72 ) = JVS( 740 )
   W( 73 ) = JVS( 741 )
  a = -W( 23 ) / JVS(          113  )
  W( 23 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 114 )
  W( 51 ) = W( 51 ) + a*JVS( 115 )
  W( 64 ) = W( 64 ) + a*JVS( 116 )
  W( 65 ) = W( 65 ) + a*JVS( 117 )
  W( 68 ) = W( 68 ) + a*JVS( 118 )
  a = -W( 26 ) / JVS(          128  )
  W( 26 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 129 )
  W( 63 ) = W( 63 ) + a*JVS( 130 )
  W( 66 ) = W( 66 ) + a*JVS( 131 )
  W( 67 ) = W( 67 ) + a*JVS( 132 )
  W( 68 ) = W( 68 ) + a*JVS( 133 )
  a = -W( 31 ) / JVS(          174  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 175 )
  W( 63 ) = W( 63 ) + a*JVS( 176 )
  W( 66 ) = W( 66 ) + a*JVS( 177 )
  W( 68 ) = W( 68 ) + a*JVS( 178 )
  W( 70 ) = W( 70 ) + a*JVS( 179 )
  W( 73 ) = W( 73 ) + a*JVS( 180 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 35 ) / JVS(          208  )
  W( 35 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 209 )
  W( 63 ) = W( 63 ) + a*JVS( 210 )
  W( 68 ) = W( 68 ) + a*JVS( 211 )
  W( 70 ) = W( 70 ) + a*JVS( 212 )
  W( 73 ) = W( 73 ) + a*JVS( 213 )
  a = -W( 36 ) / JVS(          215  )
  W( 36 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 216 )
  W( 65 ) = W( 65 ) + a*JVS( 217 )
  W( 66 ) = W( 66 ) + a*JVS( 218 )
  W( 67 ) = W( 67 ) + a*JVS( 219 )
  W( 68 ) = W( 68 ) + a*JVS( 220 )
  W( 70 ) = W( 70 ) + a*JVS( 221 )
  W( 72 ) = W( 72 ) + a*JVS( 222 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 38 ) / JVS(          228  )
  W( 38 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 229 )
  W( 65 ) = W( 65 ) + a*JVS( 230 )
  W( 66 ) = W( 66 ) + a*JVS( 231 )
  W( 67 ) = W( 67 ) + a*JVS( 232 )
  W( 68 ) = W( 68 ) + a*JVS( 233 )
  W( 70 ) = W( 70 ) + a*JVS( 234 )
  W( 72 ) = W( 72 ) + a*JVS( 235 )
  a = -W( 39 ) / JVS(          236  )
  W( 39 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 237 )
  W( 65 ) = W( 65 ) + a*JVS( 238 )
  W( 66 ) = W( 66 ) + a*JVS( 239 )
  W( 67 ) = W( 67 ) + a*JVS( 240 )
  W( 68 ) = W( 68 ) + a*JVS( 241 )
  W( 70 ) = W( 70 ) + a*JVS( 242 )
  W( 72 ) = W( 72 ) + a*JVS( 243 )
  a = -W( 40 ) / JVS(          246  )
  W( 40 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 247 )
  W( 46 ) = W( 46 ) + a*JVS( 248 )
  W( 47 ) = W( 47 ) + a*JVS( 249 )
  W( 51 ) = W( 51 ) + a*JVS( 250 )
  W( 52 ) = W( 52 ) + a*JVS( 251 )
  W( 54 ) = W( 54 ) + a*JVS( 252 )
  W( 55 ) = W( 55 ) + a*JVS( 253 )
  W( 63 ) = W( 63 ) + a*JVS( 254 )
  W( 64 ) = W( 64 ) + a*JVS( 255 )
  W( 65 ) = W( 65 ) + a*JVS( 256 )
  W( 66 ) = W( 66 ) + a*JVS( 257 )
  W( 68 ) = W( 68 ) + a*JVS( 258 )
  W( 70 ) = W( 70 ) + a*JVS( 259 )
  W( 72 ) = W( 72 ) + a*JVS( 260 )
  W( 73 ) = W( 73 ) + a*JVS( 261 )
  a = -W( 41 ) / JVS(          264  )
  W( 41 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 265 )
  W( 65 ) = W( 65 ) + a*JVS( 266 )
  W( 66 ) = W( 66 ) + a*JVS( 267 )
  W( 67 ) = W( 67 ) + a*JVS( 268 )
  W( 68 ) = W( 68 ) + a*JVS( 269 )
  W( 70 ) = W( 70 ) + a*JVS( 270 )
  W( 72 ) = W( 72 ) + a*JVS( 271 )
  W( 73 ) = W( 73 ) + a*JVS( 272 )
  a = -W( 42 ) / JVS(          278  )
  W( 42 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 279 )
  W( 45 ) = W( 45 ) + a*JVS( 280 )
  W( 49 ) = W( 49 ) + a*JVS( 281 )
  W( 52 ) = W( 52 ) + a*JVS( 282 )
  W( 54 ) = W( 54 ) + a*JVS( 283 )
  W( 57 ) = W( 57 ) + a*JVS( 284 )
  W( 58 ) = W( 58 ) + a*JVS( 285 )
  W( 59 ) = W( 59 ) + a*JVS( 286 )
  W( 60 ) = W( 60 ) + a*JVS( 287 )
  W( 63 ) = W( 63 ) + a*JVS( 288 )
  W( 64 ) = W( 64 ) + a*JVS( 289 )
  W( 65 ) = W( 65 ) + a*JVS( 290 )
  W( 66 ) = W( 66 ) + a*JVS( 291 )
  W( 67 ) = W( 67 ) + a*JVS( 292 )
  W( 68 ) = W( 68 ) + a*JVS( 293 )
  W( 69 ) = W( 69 ) + a*JVS( 294 )
  W( 70 ) = W( 70 ) + a*JVS( 295 )
  W( 71 ) = W( 71 ) + a*JVS( 296 )
  W( 72 ) = W( 72 ) + a*JVS( 297 )
  a = -W( 43 ) / JVS(          299  )
  W( 43 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 300 )
  W( 65 ) = W( 65 ) + a*JVS( 301 )
  W( 66 ) = W( 66 ) + a*JVS( 302 )
  W( 67 ) = W( 67 ) + a*JVS( 303 )
  W( 68 ) = W( 68 ) + a*JVS( 304 )
  W( 70 ) = W( 70 ) + a*JVS( 305 )
  W( 72 ) = W( 72 ) + a*JVS( 306 )
  a = -W( 44 ) / JVS(          308  )
  W( 44 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 309 )
  W( 66 ) = W( 66 ) + a*JVS( 310 )
  W( 67 ) = W( 67 ) + a*JVS( 311 )
  W( 68 ) = W( 68 ) + a*JVS( 312 )
  W( 70 ) = W( 70 ) + a*JVS( 313 )
  W( 72 ) = W( 72 ) + a*JVS( 314 )
  a = -W( 45 ) / JVS(          316  )
  W( 45 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 70 ) = W( 70 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 325 )
  W( 65 ) = W( 65 ) + a*JVS( 326 )
  W( 66 ) = W( 66 ) + a*JVS( 327 )
  W( 67 ) = W( 67 ) + a*JVS( 328 )
  W( 68 ) = W( 68 ) + a*JVS( 329 )
  W( 70 ) = W( 70 ) + a*JVS( 330 )
  W( 72 ) = W( 72 ) + a*JVS( 331 )
  W( 73 ) = W( 73 ) + a*JVS( 332 )
  a = -W( 47 ) / JVS(          334  )
  W( 47 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 335 )
  W( 65 ) = W( 65 ) + a*JVS( 336 )
  W( 66 ) = W( 66 ) + a*JVS( 337 )
  W( 67 ) = W( 67 ) + a*JVS( 338 )
  W( 68 ) = W( 68 ) + a*JVS( 339 )
  W( 70 ) = W( 70 ) + a*JVS( 340 )
  W( 72 ) = W( 72 ) + a*JVS( 341 )
  W( 73 ) = W( 73 ) + a*JVS( 342 )
  a = -W( 49 ) / JVS(          387  )
  W( 49 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 388 )
  W( 65 ) = W( 65 ) + a*JVS( 389 )
  W( 66 ) = W( 66 ) + a*JVS( 390 )
  W( 67 ) = W( 67 ) + a*JVS( 391 )
  W( 68 ) = W( 68 ) + a*JVS( 392 )
  W( 70 ) = W( 70 ) + a*JVS( 393 )
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  a = -W( 50 ) / JVS(          397  )
  W( 50 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 398 )
  W( 63 ) = W( 63 ) + a*JVS( 399 )
  W( 65 ) = W( 65 ) + a*JVS( 400 )
  W( 66 ) = W( 66 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 70 ) = W( 70 ) + a*JVS( 404 )
  W( 72 ) = W( 72 ) + a*JVS( 405 )
  a = -W( 51 ) / JVS(          411  )
  W( 51 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 412 )
  W( 63 ) = W( 63 ) + a*JVS( 413 )
  W( 65 ) = W( 65 ) + a*JVS( 414 )
  W( 66 ) = W( 66 ) + a*JVS( 415 )
  W( 67 ) = W( 67 ) + a*JVS( 416 )
  W( 68 ) = W( 68 ) + a*JVS( 417 )
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 72 ) = W( 72 ) + a*JVS( 419 )
  a = -W( 52 ) / JVS(          421  )
  W( 52 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 422 )
  W( 66 ) = W( 66 ) + a*JVS( 423 )
  W( 67 ) = W( 67 ) + a*JVS( 424 )
  W( 68 ) = W( 68 ) + a*JVS( 425 )
  W( 69 ) = W( 69 ) + a*JVS( 426 )
  W( 70 ) = W( 70 ) + a*JVS( 427 )
  W( 71 ) = W( 71 ) + a*JVS( 428 )
  W( 72 ) = W( 72 ) + a*JVS( 429 )
  a = -W( 53 ) / JVS(          442  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 443 )
  W( 57 ) = W( 57 ) + a*JVS( 444 )
  W( 58 ) = W( 58 ) + a*JVS( 445 )
  W( 59 ) = W( 59 ) + a*JVS( 446 )
  W( 60 ) = W( 60 ) + a*JVS( 447 )
  W( 61 ) = W( 61 ) + a*JVS( 448 )
  W( 62 ) = W( 62 ) + a*JVS( 449 )
  W( 63 ) = W( 63 ) + a*JVS( 450 )
  W( 64 ) = W( 64 ) + a*JVS( 451 )
  W( 65 ) = W( 65 ) + a*JVS( 452 )
  W( 66 ) = W( 66 ) + a*JVS( 453 )
  W( 67 ) = W( 67 ) + a*JVS( 454 )
  W( 68 ) = W( 68 ) + a*JVS( 455 )
  W( 69 ) = W( 69 ) + a*JVS( 456 )
  W( 70 ) = W( 70 ) + a*JVS( 457 )
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 72 ) = W( 72 ) + a*JVS( 459 )
  a = -W( 54 ) / JVS(          463  )
  W( 54 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 464 )
  W( 63 ) = W( 63 ) + a*JVS( 465 )
  W( 65 ) = W( 65 ) + a*JVS( 466 )
  W( 66 ) = W( 66 ) + a*JVS( 467 )
  W( 67 ) = W( 67 ) + a*JVS( 468 )
  W( 68 ) = W( 68 ) + a*JVS( 469 )
  W( 70 ) = W( 70 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  a = -W( 55 ) / JVS(          476  )
  W( 55 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 477 )
  W( 63 ) = W( 63 ) + a*JVS( 478 )
  W( 65 ) = W( 65 ) + a*JVS( 479 )
  W( 66 ) = W( 66 ) + a*JVS( 480 )
  W( 67 ) = W( 67 ) + a*JVS( 481 )
  W( 68 ) = W( 68 ) + a*JVS( 482 )
  W( 70 ) = W( 70 ) + a*JVS( 483 )
  W( 72 ) = W( 72 ) + a*JVS( 484 )
  W( 73 ) = W( 73 ) + a*JVS( 485 )
  a = -W( 56 ) / JVS(          500  )
  W( 56 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 501 )
  W( 63 ) = W( 63 ) + a*JVS( 502 )
  W( 64 ) = W( 64 ) + a*JVS( 503 )
  W( 65 ) = W( 65 ) + a*JVS( 504 )
  W( 66 ) = W( 66 ) + a*JVS( 505 )
  W( 67 ) = W( 67 ) + a*JVS( 506 )
  W( 68 ) = W( 68 ) + a*JVS( 507 )
  W( 69 ) = W( 69 ) + a*JVS( 508 )
  W( 70 ) = W( 70 ) + a*JVS( 509 )
  W( 71 ) = W( 71 ) + a*JVS( 510 )
  W( 72 ) = W( 72 ) + a*JVS( 511 )
  W( 73 ) = W( 73 ) + a*JVS( 512 )
  a = -W( 57 ) / JVS(          517  )
  W( 57 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 518 )
  W( 63 ) = W( 63 ) + a*JVS( 519 )
  W( 65 ) = W( 65 ) + a*JVS( 520 )
  W( 66 ) = W( 66 ) + a*JVS( 521 )
  W( 67 ) = W( 67 ) + a*JVS( 522 )
  W( 68 ) = W( 68 ) + a*JVS( 523 )
  W( 70 ) = W( 70 ) + a*JVS( 524 )
  W( 72 ) = W( 72 ) + a*JVS( 525 )
  a = -W( 58 ) / JVS(          528  )
  W( 58 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 529 )
  W( 63 ) = W( 63 ) + a*JVS( 530 )
  W( 65 ) = W( 65 ) + a*JVS( 531 )
  W( 66 ) = W( 66 ) + a*JVS( 532 )
  W( 67 ) = W( 67 ) + a*JVS( 533 )
  W( 68 ) = W( 68 ) + a*JVS( 534 )
  W( 70 ) = W( 70 ) + a*JVS( 535 )
  W( 72 ) = W( 72 ) + a*JVS( 536 )
  a = -W( 59 ) / JVS(          545  )
  W( 59 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 546 )
  W( 62 ) = W( 62 ) + a*JVS( 547 )
  W( 63 ) = W( 63 ) + a*JVS( 548 )
  W( 65 ) = W( 65 ) + a*JVS( 549 )
  W( 66 ) = W( 66 ) + a*JVS( 550 )
  W( 67 ) = W( 67 ) + a*JVS( 551 )
  W( 68 ) = W( 68 ) + a*JVS( 552 )
  W( 70 ) = W( 70 ) + a*JVS( 553 )
  W( 72 ) = W( 72 ) + a*JVS( 554 )
  a = -W( 60 ) / JVS(          563  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 564 )
  W( 63 ) = W( 63 ) + a*JVS( 565 )
  W( 65 ) = W( 65 ) + a*JVS( 566 )
  W( 66 ) = W( 66 ) + a*JVS( 567 )
  W( 67 ) = W( 67 ) + a*JVS( 568 )
  W( 68 ) = W( 68 ) + a*JVS( 569 )
  W( 70 ) = W( 70 ) + a*JVS( 570 )
  W( 72 ) = W( 72 ) + a*JVS( 571 )
  a = -W( 61 ) / JVS(          585  )
  W( 61 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 586 )
  W( 63 ) = W( 63 ) + a*JVS( 587 )
  W( 64 ) = W( 64 ) + a*JVS( 588 )
  W( 65 ) = W( 65 ) + a*JVS( 589 )
  W( 66 ) = W( 66 ) + a*JVS( 590 )
  W( 67 ) = W( 67 ) + a*JVS( 591 )
  W( 68 ) = W( 68 ) + a*JVS( 592 )
  W( 69 ) = W( 69 ) + a*JVS( 593 )
  W( 70 ) = W( 70 ) + a*JVS( 594 )
  W( 71 ) = W( 71 ) + a*JVS( 595 )
  W( 72 ) = W( 72 ) + a*JVS( 596 )
  W( 73 ) = W( 73 ) + a*JVS( 597 )
  a = -W( 62 ) / JVS(          602  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 603 )
  W( 65 ) = W( 65 ) + a*JVS( 604 )
  W( 66 ) = W( 66 ) + a*JVS( 605 )
  W( 67 ) = W( 67 ) + a*JVS( 606 )
  W( 68 ) = W( 68 ) + a*JVS( 607 )
  W( 70 ) = W( 70 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  a = -W( 63 ) / JVS(          628  )
  W( 63 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 629 )
  W( 66 ) = W( 66 ) + a*JVS( 630 )
  W( 67 ) = W( 67 ) + a*JVS( 631 )
  W( 68 ) = W( 68 ) + a*JVS( 632 )
  W( 70 ) = W( 70 ) + a*JVS( 633 )
  W( 72 ) = W( 72 ) + a*JVS( 634 )
  W( 73 ) = W( 73 ) + a*JVS( 635 )
  a = -W( 64 ) / JVS(          652  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 653 )
  W( 66 ) = W( 66 ) + a*JVS( 654 )
  W( 67 ) = W( 67 ) + a*JVS( 655 )
  W( 68 ) = W( 68 ) + a*JVS( 656 )
  W( 69 ) = W( 69 ) + a*JVS( 657 )
  W( 70 ) = W( 70 ) + a*JVS( 658 )
  W( 71 ) = W( 71 ) + a*JVS( 659 )
  W( 72 ) = W( 72 ) + a*JVS( 660 )
  W( 73 ) = W( 73 ) + a*JVS( 661 )
  a = -W( 65 ) / JVS(          691  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 692 )
  W( 67 ) = W( 67 ) + a*JVS( 693 )
  W( 68 ) = W( 68 ) + a*JVS( 694 )
  W( 69 ) = W( 69 ) + a*JVS( 695 )
  W( 70 ) = W( 70 ) + a*JVS( 696 )
  W( 71 ) = W( 71 ) + a*JVS( 697 )
  W( 72 ) = W( 72 ) + a*JVS( 698 )
  W( 73 ) = W( 73 ) + a*JVS( 699 )
  JVS( 700) = W( 23 )
  JVS( 701) = W( 26 )
  JVS( 702) = W( 31 )
  JVS( 703) = W( 33 )
  JVS( 704) = W( 35 )
  JVS( 705) = W( 36 )
  JVS( 706) = W( 37 )
  JVS( 707) = W( 38 )
  JVS( 708) = W( 39 )
  JVS( 709) = W( 40 )
  JVS( 710) = W( 41 )
  JVS( 711) = W( 42 )
  JVS( 712) = W( 43 )
  JVS( 713) = W( 44 )
  JVS( 714) = W( 45 )
  JVS( 715) = W( 46 )
  JVS( 716) = W( 47 )
  JVS( 717) = W( 49 )
  JVS( 718) = W( 50 )
  JVS( 719) = W( 51 )
  JVS( 720) = W( 52 )
  JVS( 721) = W( 53 )
  JVS( 722) = W( 54 )
  JVS( 723) = W( 55 )
  JVS( 724) = W( 56 )
  JVS( 725) = W( 57 )
  JVS( 726) = W( 58 )
  JVS( 727) = W( 59 )
  JVS( 728) = W( 60 )
  JVS( 729) = W( 61 )
  JVS( 730) = W( 62 )
  JVS( 731) = W( 63 )
  JVS( 732) = W( 64 )
  JVS( 733) = W( 65 )
  JVS( 734) = W( 66 )
  JVS( 735) = W( 67 )
  JVS( 736) = W( 68 )
  JVS( 737) = W( 69 )
  JVS( 738) = W( 70 )
  JVS( 739) = W( 71 )
  JVS( 740) = W( 72 )
  JVS( 741) = W( 73 )
  IF ( ABS(  JVS( 796 )) < TINY(a) ) THEN
         IER = 67                                      
         RETURN
  END IF
   W( 5 ) = JVS( 742 )
   W( 7 ) = JVS( 743 )
   W( 8 ) = JVS( 744 )
   W( 9 ) = JVS( 745 )
   W( 11 ) = JVS( 746 )
   W( 12 ) = JVS( 747 )
   W( 15 ) = JVS( 748 )
   W( 16 ) = JVS( 749 )
   W( 17 ) = JVS( 750 )
   W( 19 ) = JVS( 751 )
   W( 20 ) = JVS( 752 )
   W( 21 ) = JVS( 753 )
   W( 22 ) = JVS( 754 )
   W( 23 ) = JVS( 755 )
   W( 24 ) = JVS( 756 )
   W( 25 ) = JVS( 757 )
   W( 26 ) = JVS( 758 )
   W( 28 ) = JVS( 759 )
   W( 29 ) = JVS( 760 )
   W( 30 ) = JVS( 761 )
   W( 32 ) = JVS( 762 )
   W( 33 ) = JVS( 763 )
   W( 34 ) = JVS( 764 )
   W( 35 ) = JVS( 765 )
   W( 36 ) = JVS( 766 )
   W( 37 ) = JVS( 767 )
   W( 38 ) = JVS( 768 )
   W( 39 ) = JVS( 769 )
   W( 40 ) = JVS( 770 )
   W( 41 ) = JVS( 771 )
   W( 43 ) = JVS( 772 )
   W( 44 ) = JVS( 773 )
   W( 45 ) = JVS( 774 )
   W( 46 ) = JVS( 775 )
   W( 47 ) = JVS( 776 )
   W( 48 ) = JVS( 777 )
   W( 49 ) = JVS( 778 )
   W( 50 ) = JVS( 779 )
   W( 51 ) = JVS( 780 )
   W( 52 ) = JVS( 781 )
   W( 53 ) = JVS( 782 )
   W( 54 ) = JVS( 783 )
   W( 55 ) = JVS( 784 )
   W( 56 ) = JVS( 785 )
   W( 57 ) = JVS( 786 )
   W( 58 ) = JVS( 787 )
   W( 59 ) = JVS( 788 )
   W( 60 ) = JVS( 789 )
   W( 61 ) = JVS( 790 )
   W( 62 ) = JVS( 791 )
   W( 63 ) = JVS( 792 )
   W( 64 ) = JVS( 793 )
   W( 65 ) = JVS( 794 )
   W( 66 ) = JVS( 795 )
   W( 67 ) = JVS( 796 )
   W( 68 ) = JVS( 797 )
   W( 69 ) = JVS( 798 )
   W( 70 ) = JVS( 799 )
   W( 71 ) = JVS( 800 )
   W( 72 ) = JVS( 801 )
   W( 73 ) = JVS( 802 )
  a = -W( 5 ) / JVS(           41  )
  W( 5 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 42 )
  a = -W( 7 ) / JVS(           45  )
  W( 7 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 46 )
  a = -W( 8 ) / JVS(           47  )
  W( 8 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 48 )
  a = -W( 9 ) / JVS(           49  )
  W( 9 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 50 )
  a = -W( 11 ) / JVS(           54  )
  W( 11 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 55 )
  a = -W( 12 ) / JVS(           56  )
  W( 12 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 57 )
  a = -W( 15 ) / JVS(           68  )
  W( 15 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 69 )
  W( 68 ) = W( 68 ) + a*JVS( 70 )
  a = -W( 16 ) / JVS(           71  )
  W( 16 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 72 )
  W( 68 ) = W( 68 ) + a*JVS( 73 )
  W( 73 ) = W( 73 ) + a*JVS( 74 )
  a = -W( 17 ) / JVS(           75  )
  W( 17 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 76 )
  W( 67 ) = W( 67 ) + a*JVS( 77 )
  W( 68 ) = W( 68 ) + a*JVS( 78 )
  a = -W( 19 ) / JVS(           86  )
  W( 19 ) = -a
  W( 29 ) = W( 29 ) + a*JVS( 87 )
  W( 30 ) = W( 30 ) + a*JVS( 88 )
  W( 33 ) = W( 33 ) + a*JVS( 89 )
  W( 37 ) = W( 37 ) + a*JVS( 90 )
  W( 57 ) = W( 57 ) + a*JVS( 91 )
  W( 58 ) = W( 58 ) + a*JVS( 92 )
  W( 63 ) = W( 63 ) + a*JVS( 93 )
  W( 67 ) = W( 67 ) + a*JVS( 94 )
  W( 68 ) = W( 68 ) + a*JVS( 95 )
  a = -W( 20 ) / JVS(           96  )
  W( 20 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 97 )
  W( 67 ) = W( 67 ) + a*JVS( 98 )
  W( 68 ) = W( 68 ) + a*JVS( 99 )
  W( 70 ) = W( 70 ) + a*JVS( 100 )
  W( 73 ) = W( 73 ) + a*JVS( 101 )
  a = -W( 21 ) / JVS(          103  )
  W( 21 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 104 )
  W( 68 ) = W( 68 ) + a*JVS( 105 )
  W( 73 ) = W( 73 ) + a*JVS( 106 )
  a = -W( 22 ) / JVS(          108  )
  W( 22 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 109 )
  W( 68 ) = W( 68 ) + a*JVS( 110 )
  W( 73 ) = W( 73 ) + a*JVS( 111 )
  a = -W( 23 ) / JVS(          113  )
  W( 23 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 114 )
  W( 51 ) = W( 51 ) + a*JVS( 115 )
  W( 64 ) = W( 64 ) + a*JVS( 116 )
  W( 65 ) = W( 65 ) + a*JVS( 117 )
  W( 68 ) = W( 68 ) + a*JVS( 118 )
  a = -W( 24 ) / JVS(          119  )
  W( 24 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 120 )
  W( 68 ) = W( 68 ) + a*JVS( 121 )
  W( 70 ) = W( 70 ) + a*JVS( 122 )
  a = -W( 25 ) / JVS(          123  )
  W( 25 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 124 )
  W( 63 ) = W( 63 ) + a*JVS( 125 )
  W( 68 ) = W( 68 ) + a*JVS( 126 )
  W( 73 ) = W( 73 ) + a*JVS( 127 )
  a = -W( 26 ) / JVS(          128  )
  W( 26 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 129 )
  W( 63 ) = W( 63 ) + a*JVS( 130 )
  W( 66 ) = W( 66 ) + a*JVS( 131 )
  W( 67 ) = W( 67 ) + a*JVS( 132 )
  W( 68 ) = W( 68 ) + a*JVS( 133 )
  a = -W( 28 ) / JVS(          148  )
  W( 28 ) = -a
  W( 29 ) = W( 29 ) + a*JVS( 149 )
  W( 30 ) = W( 30 ) + a*JVS( 150 )
  W( 33 ) = W( 33 ) + a*JVS( 151 )
  W( 34 ) = W( 34 ) + a*JVS( 152 )
  W( 35 ) = W( 35 ) + a*JVS( 153 )
  W( 37 ) = W( 37 ) + a*JVS( 154 )
  W( 40 ) = W( 40 ) + a*JVS( 155 )
  W( 48 ) = W( 48 ) + a*JVS( 156 )
  W( 51 ) = W( 51 ) + a*JVS( 157 )
  W( 53 ) = W( 53 ) + a*JVS( 158 )
  W( 54 ) = W( 54 ) + a*JVS( 159 )
  W( 57 ) = W( 57 ) + a*JVS( 160 )
  W( 58 ) = W( 58 ) + a*JVS( 161 )
  W( 62 ) = W( 62 ) + a*JVS( 162 )
  W( 63 ) = W( 63 ) + a*JVS( 163 )
  W( 68 ) = W( 68 ) + a*JVS( 164 )
  W( 70 ) = W( 70 ) + a*JVS( 165 )
  a = -W( 29 ) / JVS(          166  )
  W( 29 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 167 )
  W( 68 ) = W( 68 ) + a*JVS( 168 )
  W( 70 ) = W( 70 ) + a*JVS( 169 )
  a = -W( 30 ) / JVS(          170  )
  W( 30 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 171 )
  W( 68 ) = W( 68 ) + a*JVS( 172 )
  W( 70 ) = W( 70 ) + a*JVS( 173 )
  a = -W( 32 ) / JVS(          185  )
  W( 32 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 186 )
  W( 67 ) = W( 67 ) + a*JVS( 187 )
  W( 68 ) = W( 68 ) + a*JVS( 188 )
  W( 70 ) = W( 70 ) + a*JVS( 189 )
  W( 73 ) = W( 73 ) + a*JVS( 190 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 34 ) / JVS(          196  )
  W( 34 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 197 )
  W( 46 ) = W( 46 ) + a*JVS( 198 )
  W( 47 ) = W( 47 ) + a*JVS( 199 )
  W( 52 ) = W( 52 ) + a*JVS( 200 )
  W( 54 ) = W( 54 ) + a*JVS( 201 )
  W( 63 ) = W( 63 ) + a*JVS( 202 )
  W( 65 ) = W( 65 ) + a*JVS( 203 )
  W( 66 ) = W( 66 ) + a*JVS( 204 )
  W( 68 ) = W( 68 ) + a*JVS( 205 )
  W( 70 ) = W( 70 ) + a*JVS( 206 )
  W( 72 ) = W( 72 ) + a*JVS( 207 )
  a = -W( 35 ) / JVS(          208  )
  W( 35 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 209 )
  W( 63 ) = W( 63 ) + a*JVS( 210 )
  W( 68 ) = W( 68 ) + a*JVS( 211 )
  W( 70 ) = W( 70 ) + a*JVS( 212 )
  W( 73 ) = W( 73 ) + a*JVS( 213 )
  a = -W( 36 ) / JVS(          215  )
  W( 36 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 216 )
  W( 65 ) = W( 65 ) + a*JVS( 217 )
  W( 66 ) = W( 66 ) + a*JVS( 218 )
  W( 67 ) = W( 67 ) + a*JVS( 219 )
  W( 68 ) = W( 68 ) + a*JVS( 220 )
  W( 70 ) = W( 70 ) + a*JVS( 221 )
  W( 72 ) = W( 72 ) + a*JVS( 222 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 38 ) / JVS(          228  )
  W( 38 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 229 )
  W( 65 ) = W( 65 ) + a*JVS( 230 )
  W( 66 ) = W( 66 ) + a*JVS( 231 )
  W( 67 ) = W( 67 ) + a*JVS( 232 )
  W( 68 ) = W( 68 ) + a*JVS( 233 )
  W( 70 ) = W( 70 ) + a*JVS( 234 )
  W( 72 ) = W( 72 ) + a*JVS( 235 )
  a = -W( 39 ) / JVS(          236  )
  W( 39 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 237 )
  W( 65 ) = W( 65 ) + a*JVS( 238 )
  W( 66 ) = W( 66 ) + a*JVS( 239 )
  W( 67 ) = W( 67 ) + a*JVS( 240 )
  W( 68 ) = W( 68 ) + a*JVS( 241 )
  W( 70 ) = W( 70 ) + a*JVS( 242 )
  W( 72 ) = W( 72 ) + a*JVS( 243 )
  a = -W( 40 ) / JVS(          246  )
  W( 40 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 247 )
  W( 46 ) = W( 46 ) + a*JVS( 248 )
  W( 47 ) = W( 47 ) + a*JVS( 249 )
  W( 51 ) = W( 51 ) + a*JVS( 250 )
  W( 52 ) = W( 52 ) + a*JVS( 251 )
  W( 54 ) = W( 54 ) + a*JVS( 252 )
  W( 55 ) = W( 55 ) + a*JVS( 253 )
  W( 63 ) = W( 63 ) + a*JVS( 254 )
  W( 64 ) = W( 64 ) + a*JVS( 255 )
  W( 65 ) = W( 65 ) + a*JVS( 256 )
  W( 66 ) = W( 66 ) + a*JVS( 257 )
  W( 68 ) = W( 68 ) + a*JVS( 258 )
  W( 70 ) = W( 70 ) + a*JVS( 259 )
  W( 72 ) = W( 72 ) + a*JVS( 260 )
  W( 73 ) = W( 73 ) + a*JVS( 261 )
  a = -W( 41 ) / JVS(          264  )
  W( 41 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 265 )
  W( 65 ) = W( 65 ) + a*JVS( 266 )
  W( 66 ) = W( 66 ) + a*JVS( 267 )
  W( 67 ) = W( 67 ) + a*JVS( 268 )
  W( 68 ) = W( 68 ) + a*JVS( 269 )
  W( 70 ) = W( 70 ) + a*JVS( 270 )
  W( 72 ) = W( 72 ) + a*JVS( 271 )
  W( 73 ) = W( 73 ) + a*JVS( 272 )
  a = -W( 43 ) / JVS(          299  )
  W( 43 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 300 )
  W( 65 ) = W( 65 ) + a*JVS( 301 )
  W( 66 ) = W( 66 ) + a*JVS( 302 )
  W( 67 ) = W( 67 ) + a*JVS( 303 )
  W( 68 ) = W( 68 ) + a*JVS( 304 )
  W( 70 ) = W( 70 ) + a*JVS( 305 )
  W( 72 ) = W( 72 ) + a*JVS( 306 )
  a = -W( 44 ) / JVS(          308  )
  W( 44 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 309 )
  W( 66 ) = W( 66 ) + a*JVS( 310 )
  W( 67 ) = W( 67 ) + a*JVS( 311 )
  W( 68 ) = W( 68 ) + a*JVS( 312 )
  W( 70 ) = W( 70 ) + a*JVS( 313 )
  W( 72 ) = W( 72 ) + a*JVS( 314 )
  a = -W( 45 ) / JVS(          316  )
  W( 45 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 70 ) = W( 70 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 325 )
  W( 65 ) = W( 65 ) + a*JVS( 326 )
  W( 66 ) = W( 66 ) + a*JVS( 327 )
  W( 67 ) = W( 67 ) + a*JVS( 328 )
  W( 68 ) = W( 68 ) + a*JVS( 329 )
  W( 70 ) = W( 70 ) + a*JVS( 330 )
  W( 72 ) = W( 72 ) + a*JVS( 331 )
  W( 73 ) = W( 73 ) + a*JVS( 332 )
  a = -W( 47 ) / JVS(          334  )
  W( 47 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 335 )
  W( 65 ) = W( 65 ) + a*JVS( 336 )
  W( 66 ) = W( 66 ) + a*JVS( 337 )
  W( 67 ) = W( 67 ) + a*JVS( 338 )
  W( 68 ) = W( 68 ) + a*JVS( 339 )
  W( 70 ) = W( 70 ) + a*JVS( 340 )
  W( 72 ) = W( 72 ) + a*JVS( 341 )
  W( 73 ) = W( 73 ) + a*JVS( 342 )
  a = -W( 48 ) / JVS(          363  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 364 )
  W( 50 ) = W( 50 ) + a*JVS( 365 )
  W( 51 ) = W( 51 ) + a*JVS( 366 )
  W( 52 ) = W( 52 ) + a*JVS( 367 )
  W( 54 ) = W( 54 ) + a*JVS( 368 )
  W( 55 ) = W( 55 ) + a*JVS( 369 )
  W( 56 ) = W( 56 ) + a*JVS( 370 )
  W( 57 ) = W( 57 ) + a*JVS( 371 )
  W( 58 ) = W( 58 ) + a*JVS( 372 )
  W( 59 ) = W( 59 ) + a*JVS( 373 )
  W( 60 ) = W( 60 ) + a*JVS( 374 )
  W( 61 ) = W( 61 ) + a*JVS( 375 )
  W( 62 ) = W( 62 ) + a*JVS( 376 )
  W( 63 ) = W( 63 ) + a*JVS( 377 )
  W( 64 ) = W( 64 ) + a*JVS( 378 )
  W( 65 ) = W( 65 ) + a*JVS( 379 )
  W( 66 ) = W( 66 ) + a*JVS( 380 )
  W( 67 ) = W( 67 ) + a*JVS( 381 )
  W( 68 ) = W( 68 ) + a*JVS( 382 )
  W( 70 ) = W( 70 ) + a*JVS( 383 )
  W( 72 ) = W( 72 ) + a*JVS( 384 )
  W( 73 ) = W( 73 ) + a*JVS( 385 )
  a = -W( 49 ) / JVS(          387  )
  W( 49 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 388 )
  W( 65 ) = W( 65 ) + a*JVS( 389 )
  W( 66 ) = W( 66 ) + a*JVS( 390 )
  W( 67 ) = W( 67 ) + a*JVS( 391 )
  W( 68 ) = W( 68 ) + a*JVS( 392 )
  W( 70 ) = W( 70 ) + a*JVS( 393 )
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  a = -W( 50 ) / JVS(          397  )
  W( 50 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 398 )
  W( 63 ) = W( 63 ) + a*JVS( 399 )
  W( 65 ) = W( 65 ) + a*JVS( 400 )
  W( 66 ) = W( 66 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 70 ) = W( 70 ) + a*JVS( 404 )
  W( 72 ) = W( 72 ) + a*JVS( 405 )
  a = -W( 51 ) / JVS(          411  )
  W( 51 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 412 )
  W( 63 ) = W( 63 ) + a*JVS( 413 )
  W( 65 ) = W( 65 ) + a*JVS( 414 )
  W( 66 ) = W( 66 ) + a*JVS( 415 )
  W( 67 ) = W( 67 ) + a*JVS( 416 )
  W( 68 ) = W( 68 ) + a*JVS( 417 )
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 72 ) = W( 72 ) + a*JVS( 419 )
  a = -W( 52 ) / JVS(          421  )
  W( 52 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 422 )
  W( 66 ) = W( 66 ) + a*JVS( 423 )
  W( 67 ) = W( 67 ) + a*JVS( 424 )
  W( 68 ) = W( 68 ) + a*JVS( 425 )
  W( 69 ) = W( 69 ) + a*JVS( 426 )
  W( 70 ) = W( 70 ) + a*JVS( 427 )
  W( 71 ) = W( 71 ) + a*JVS( 428 )
  W( 72 ) = W( 72 ) + a*JVS( 429 )
  a = -W( 53 ) / JVS(          442  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 443 )
  W( 57 ) = W( 57 ) + a*JVS( 444 )
  W( 58 ) = W( 58 ) + a*JVS( 445 )
  W( 59 ) = W( 59 ) + a*JVS( 446 )
  W( 60 ) = W( 60 ) + a*JVS( 447 )
  W( 61 ) = W( 61 ) + a*JVS( 448 )
  W( 62 ) = W( 62 ) + a*JVS( 449 )
  W( 63 ) = W( 63 ) + a*JVS( 450 )
  W( 64 ) = W( 64 ) + a*JVS( 451 )
  W( 65 ) = W( 65 ) + a*JVS( 452 )
  W( 66 ) = W( 66 ) + a*JVS( 453 )
  W( 67 ) = W( 67 ) + a*JVS( 454 )
  W( 68 ) = W( 68 ) + a*JVS( 455 )
  W( 69 ) = W( 69 ) + a*JVS( 456 )
  W( 70 ) = W( 70 ) + a*JVS( 457 )
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 72 ) = W( 72 ) + a*JVS( 459 )
  a = -W( 54 ) / JVS(          463  )
  W( 54 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 464 )
  W( 63 ) = W( 63 ) + a*JVS( 465 )
  W( 65 ) = W( 65 ) + a*JVS( 466 )
  W( 66 ) = W( 66 ) + a*JVS( 467 )
  W( 67 ) = W( 67 ) + a*JVS( 468 )
  W( 68 ) = W( 68 ) + a*JVS( 469 )
  W( 70 ) = W( 70 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  a = -W( 55 ) / JVS(          476  )
  W( 55 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 477 )
  W( 63 ) = W( 63 ) + a*JVS( 478 )
  W( 65 ) = W( 65 ) + a*JVS( 479 )
  W( 66 ) = W( 66 ) + a*JVS( 480 )
  W( 67 ) = W( 67 ) + a*JVS( 481 )
  W( 68 ) = W( 68 ) + a*JVS( 482 )
  W( 70 ) = W( 70 ) + a*JVS( 483 )
  W( 72 ) = W( 72 ) + a*JVS( 484 )
  W( 73 ) = W( 73 ) + a*JVS( 485 )
  a = -W( 56 ) / JVS(          500  )
  W( 56 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 501 )
  W( 63 ) = W( 63 ) + a*JVS( 502 )
  W( 64 ) = W( 64 ) + a*JVS( 503 )
  W( 65 ) = W( 65 ) + a*JVS( 504 )
  W( 66 ) = W( 66 ) + a*JVS( 505 )
  W( 67 ) = W( 67 ) + a*JVS( 506 )
  W( 68 ) = W( 68 ) + a*JVS( 507 )
  W( 69 ) = W( 69 ) + a*JVS( 508 )
  W( 70 ) = W( 70 ) + a*JVS( 509 )
  W( 71 ) = W( 71 ) + a*JVS( 510 )
  W( 72 ) = W( 72 ) + a*JVS( 511 )
  W( 73 ) = W( 73 ) + a*JVS( 512 )
  a = -W( 57 ) / JVS(          517  )
  W( 57 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 518 )
  W( 63 ) = W( 63 ) + a*JVS( 519 )
  W( 65 ) = W( 65 ) + a*JVS( 520 )
  W( 66 ) = W( 66 ) + a*JVS( 521 )
  W( 67 ) = W( 67 ) + a*JVS( 522 )
  W( 68 ) = W( 68 ) + a*JVS( 523 )
  W( 70 ) = W( 70 ) + a*JVS( 524 )
  W( 72 ) = W( 72 ) + a*JVS( 525 )
  a = -W( 58 ) / JVS(          528  )
  W( 58 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 529 )
  W( 63 ) = W( 63 ) + a*JVS( 530 )
  W( 65 ) = W( 65 ) + a*JVS( 531 )
  W( 66 ) = W( 66 ) + a*JVS( 532 )
  W( 67 ) = W( 67 ) + a*JVS( 533 )
  W( 68 ) = W( 68 ) + a*JVS( 534 )
  W( 70 ) = W( 70 ) + a*JVS( 535 )
  W( 72 ) = W( 72 ) + a*JVS( 536 )
  a = -W( 59 ) / JVS(          545  )
  W( 59 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 546 )
  W( 62 ) = W( 62 ) + a*JVS( 547 )
  W( 63 ) = W( 63 ) + a*JVS( 548 )
  W( 65 ) = W( 65 ) + a*JVS( 549 )
  W( 66 ) = W( 66 ) + a*JVS( 550 )
  W( 67 ) = W( 67 ) + a*JVS( 551 )
  W( 68 ) = W( 68 ) + a*JVS( 552 )
  W( 70 ) = W( 70 ) + a*JVS( 553 )
  W( 72 ) = W( 72 ) + a*JVS( 554 )
  a = -W( 60 ) / JVS(          563  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 564 )
  W( 63 ) = W( 63 ) + a*JVS( 565 )
  W( 65 ) = W( 65 ) + a*JVS( 566 )
  W( 66 ) = W( 66 ) + a*JVS( 567 )
  W( 67 ) = W( 67 ) + a*JVS( 568 )
  W( 68 ) = W( 68 ) + a*JVS( 569 )
  W( 70 ) = W( 70 ) + a*JVS( 570 )
  W( 72 ) = W( 72 ) + a*JVS( 571 )
  a = -W( 61 ) / JVS(          585  )
  W( 61 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 586 )
  W( 63 ) = W( 63 ) + a*JVS( 587 )
  W( 64 ) = W( 64 ) + a*JVS( 588 )
  W( 65 ) = W( 65 ) + a*JVS( 589 )
  W( 66 ) = W( 66 ) + a*JVS( 590 )
  W( 67 ) = W( 67 ) + a*JVS( 591 )
  W( 68 ) = W( 68 ) + a*JVS( 592 )
  W( 69 ) = W( 69 ) + a*JVS( 593 )
  W( 70 ) = W( 70 ) + a*JVS( 594 )
  W( 71 ) = W( 71 ) + a*JVS( 595 )
  W( 72 ) = W( 72 ) + a*JVS( 596 )
  W( 73 ) = W( 73 ) + a*JVS( 597 )
  a = -W( 62 ) / JVS(          602  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 603 )
  W( 65 ) = W( 65 ) + a*JVS( 604 )
  W( 66 ) = W( 66 ) + a*JVS( 605 )
  W( 67 ) = W( 67 ) + a*JVS( 606 )
  W( 68 ) = W( 68 ) + a*JVS( 607 )
  W( 70 ) = W( 70 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  a = -W( 63 ) / JVS(          628  )
  W( 63 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 629 )
  W( 66 ) = W( 66 ) + a*JVS( 630 )
  W( 67 ) = W( 67 ) + a*JVS( 631 )
  W( 68 ) = W( 68 ) + a*JVS( 632 )
  W( 70 ) = W( 70 ) + a*JVS( 633 )
  W( 72 ) = W( 72 ) + a*JVS( 634 )
  W( 73 ) = W( 73 ) + a*JVS( 635 )
  a = -W( 64 ) / JVS(          652  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 653 )
  W( 66 ) = W( 66 ) + a*JVS( 654 )
  W( 67 ) = W( 67 ) + a*JVS( 655 )
  W( 68 ) = W( 68 ) + a*JVS( 656 )
  W( 69 ) = W( 69 ) + a*JVS( 657 )
  W( 70 ) = W( 70 ) + a*JVS( 658 )
  W( 71 ) = W( 71 ) + a*JVS( 659 )
  W( 72 ) = W( 72 ) + a*JVS( 660 )
  W( 73 ) = W( 73 ) + a*JVS( 661 )
  a = -W( 65 ) / JVS(          691  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 692 )
  W( 67 ) = W( 67 ) + a*JVS( 693 )
  W( 68 ) = W( 68 ) + a*JVS( 694 )
  W( 69 ) = W( 69 ) + a*JVS( 695 )
  W( 70 ) = W( 70 ) + a*JVS( 696 )
  W( 71 ) = W( 71 ) + a*JVS( 697 )
  W( 72 ) = W( 72 ) + a*JVS( 698 )
  W( 73 ) = W( 73 ) + a*JVS( 699 )
  a = -W( 66 ) / JVS(          734  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 735 )
  W( 68 ) = W( 68 ) + a*JVS( 736 )
  W( 69 ) = W( 69 ) + a*JVS( 737 )
  W( 70 ) = W( 70 ) + a*JVS( 738 )
  W( 71 ) = W( 71 ) + a*JVS( 739 )
  W( 72 ) = W( 72 ) + a*JVS( 740 )
  W( 73 ) = W( 73 ) + a*JVS( 741 )
  JVS( 742) = W( 5 )
  JVS( 743) = W( 7 )
  JVS( 744) = W( 8 )
  JVS( 745) = W( 9 )
  JVS( 746) = W( 11 )
  JVS( 747) = W( 12 )
  JVS( 748) = W( 15 )
  JVS( 749) = W( 16 )
  JVS( 750) = W( 17 )
  JVS( 751) = W( 19 )
  JVS( 752) = W( 20 )
  JVS( 753) = W( 21 )
  JVS( 754) = W( 22 )
  JVS( 755) = W( 23 )
  JVS( 756) = W( 24 )
  JVS( 757) = W( 25 )
  JVS( 758) = W( 26 )
  JVS( 759) = W( 28 )
  JVS( 760) = W( 29 )
  JVS( 761) = W( 30 )
  JVS( 762) = W( 32 )
  JVS( 763) = W( 33 )
  JVS( 764) = W( 34 )
  JVS( 765) = W( 35 )
  JVS( 766) = W( 36 )
  JVS( 767) = W( 37 )
  JVS( 768) = W( 38 )
  JVS( 769) = W( 39 )
  JVS( 770) = W( 40 )
  JVS( 771) = W( 41 )
  JVS( 772) = W( 43 )
  JVS( 773) = W( 44 )
  JVS( 774) = W( 45 )
  JVS( 775) = W( 46 )
  JVS( 776) = W( 47 )
  JVS( 777) = W( 48 )
  JVS( 778) = W( 49 )
  JVS( 779) = W( 50 )
  JVS( 780) = W( 51 )
  JVS( 781) = W( 52 )
  JVS( 782) = W( 53 )
  JVS( 783) = W( 54 )
  JVS( 784) = W( 55 )
  JVS( 785) = W( 56 )
  JVS( 786) = W( 57 )
  JVS( 787) = W( 58 )
  JVS( 788) = W( 59 )
  JVS( 789) = W( 60 )
  JVS( 790) = W( 61 )
  JVS( 791) = W( 62 )
  JVS( 792) = W( 63 )
  JVS( 793) = W( 64 )
  JVS( 794) = W( 65 )
  JVS( 795) = W( 66 )
  JVS( 796) = W( 67 )
  JVS( 797) = W( 68 )
  JVS( 798) = W( 69 )
  JVS( 799) = W( 70 )
  JVS( 800) = W( 71 )
  JVS( 801) = W( 72 )
  JVS( 802) = W( 73 )
  IF ( ABS(  JVS( 860 )) < TINY(a) ) THEN
         IER = 68                                      
         RETURN
  END IF
   W( 5 ) = JVS( 803 )
   W( 6 ) = JVS( 804 )
   W( 7 ) = JVS( 805 )
   W( 8 ) = JVS( 806 )
   W( 9 ) = JVS( 807 )
   W( 11 ) = JVS( 808 )
   W( 12 ) = JVS( 809 )
   W( 13 ) = JVS( 810 )
   W( 14 ) = JVS( 811 )
   W( 15 ) = JVS( 812 )
   W( 16 ) = JVS( 813 )
   W( 17 ) = JVS( 814 )
   W( 18 ) = JVS( 815 )
   W( 19 ) = JVS( 816 )
   W( 21 ) = JVS( 817 )
   W( 22 ) = JVS( 818 )
   W( 23 ) = JVS( 819 )
   W( 24 ) = JVS( 820 )
   W( 25 ) = JVS( 821 )
   W( 26 ) = JVS( 822 )
   W( 27 ) = JVS( 823 )
   W( 28 ) = JVS( 824 )
   W( 29 ) = JVS( 825 )
   W( 30 ) = JVS( 826 )
   W( 31 ) = JVS( 827 )
   W( 32 ) = JVS( 828 )
   W( 33 ) = JVS( 829 )
   W( 34 ) = JVS( 830 )
   W( 35 ) = JVS( 831 )
   W( 37 ) = JVS( 832 )
   W( 40 ) = JVS( 833 )
   W( 41 ) = JVS( 834 )
   W( 42 ) = JVS( 835 )
   W( 44 ) = JVS( 836 )
   W( 45 ) = JVS( 837 )
   W( 46 ) = JVS( 838 )
   W( 47 ) = JVS( 839 )
   W( 48 ) = JVS( 840 )
   W( 49 ) = JVS( 841 )
   W( 50 ) = JVS( 842 )
   W( 51 ) = JVS( 843 )
   W( 52 ) = JVS( 844 )
   W( 53 ) = JVS( 845 )
   W( 54 ) = JVS( 846 )
   W( 55 ) = JVS( 847 )
   W( 56 ) = JVS( 848 )
   W( 57 ) = JVS( 849 )
   W( 58 ) = JVS( 850 )
   W( 59 ) = JVS( 851 )
   W( 60 ) = JVS( 852 )
   W( 61 ) = JVS( 853 )
   W( 62 ) = JVS( 854 )
   W( 63 ) = JVS( 855 )
   W( 64 ) = JVS( 856 )
   W( 65 ) = JVS( 857 )
   W( 66 ) = JVS( 858 )
   W( 67 ) = JVS( 859 )
   W( 68 ) = JVS( 860 )
   W( 69 ) = JVS( 861 )
   W( 70 ) = JVS( 862 )
   W( 71 ) = JVS( 863 )
   W( 72 ) = JVS( 864 )
   W( 73 ) = JVS( 865 )
  a = -W( 5 ) / JVS(           41  )
  W( 5 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 42 )
  a = -W( 6 ) / JVS(           43  )
  W( 6 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 44 )
  a = -W( 7 ) / JVS(           45  )
  W( 7 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 46 )
  a = -W( 8 ) / JVS(           47  )
  W( 8 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 48 )
  a = -W( 9 ) / JVS(           49  )
  W( 9 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 50 )
  a = -W( 11 ) / JVS(           54  )
  W( 11 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 55 )
  a = -W( 12 ) / JVS(           56  )
  W( 12 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 57 )
  a = -W( 13 ) / JVS(           58  )
  W( 13 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 59 )
  W( 58 ) = W( 58 ) + a*JVS( 60 )
  W( 63 ) = W( 63 ) + a*JVS( 61 )
  W( 68 ) = W( 68 ) + a*JVS( 62 )
  a = -W( 14 ) / JVS(           63  )
  W( 14 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 64 )
  W( 58 ) = W( 58 ) + a*JVS( 65 )
  W( 63 ) = W( 63 ) + a*JVS( 66 )
  W( 68 ) = W( 68 ) + a*JVS( 67 )
  a = -W( 15 ) / JVS(           68  )
  W( 15 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 69 )
  W( 68 ) = W( 68 ) + a*JVS( 70 )
  a = -W( 16 ) / JVS(           71  )
  W( 16 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 72 )
  W( 68 ) = W( 68 ) + a*JVS( 73 )
  W( 73 ) = W( 73 ) + a*JVS( 74 )
  a = -W( 17 ) / JVS(           75  )
  W( 17 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 76 )
  W( 67 ) = W( 67 ) + a*JVS( 77 )
  W( 68 ) = W( 68 ) + a*JVS( 78 )
  a = -W( 18 ) / JVS(           79  )
  W( 18 ) = -a
  W( 21 ) = W( 21 ) + a*JVS( 80 )
  W( 22 ) = W( 22 ) + a*JVS( 81 )
  W( 25 ) = W( 25 ) + a*JVS( 82 )
  W( 68 ) = W( 68 ) + a*JVS( 83 )
  W( 72 ) = W( 72 ) + a*JVS( 84 )
  W( 73 ) = W( 73 ) + a*JVS( 85 )
  a = -W( 19 ) / JVS(           86  )
  W( 19 ) = -a
  W( 29 ) = W( 29 ) + a*JVS( 87 )
  W( 30 ) = W( 30 ) + a*JVS( 88 )
  W( 33 ) = W( 33 ) + a*JVS( 89 )
  W( 37 ) = W( 37 ) + a*JVS( 90 )
  W( 57 ) = W( 57 ) + a*JVS( 91 )
  W( 58 ) = W( 58 ) + a*JVS( 92 )
  W( 63 ) = W( 63 ) + a*JVS( 93 )
  W( 67 ) = W( 67 ) + a*JVS( 94 )
  W( 68 ) = W( 68 ) + a*JVS( 95 )
  a = -W( 21 ) / JVS(          103  )
  W( 21 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 104 )
  W( 68 ) = W( 68 ) + a*JVS( 105 )
  W( 73 ) = W( 73 ) + a*JVS( 106 )
  a = -W( 22 ) / JVS(          108  )
  W( 22 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 109 )
  W( 68 ) = W( 68 ) + a*JVS( 110 )
  W( 73 ) = W( 73 ) + a*JVS( 111 )
  a = -W( 23 ) / JVS(          113  )
  W( 23 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 114 )
  W( 51 ) = W( 51 ) + a*JVS( 115 )
  W( 64 ) = W( 64 ) + a*JVS( 116 )
  W( 65 ) = W( 65 ) + a*JVS( 117 )
  W( 68 ) = W( 68 ) + a*JVS( 118 )
  a = -W( 24 ) / JVS(          119  )
  W( 24 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 120 )
  W( 68 ) = W( 68 ) + a*JVS( 121 )
  W( 70 ) = W( 70 ) + a*JVS( 122 )
  a = -W( 25 ) / JVS(          123  )
  W( 25 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 124 )
  W( 63 ) = W( 63 ) + a*JVS( 125 )
  W( 68 ) = W( 68 ) + a*JVS( 126 )
  W( 73 ) = W( 73 ) + a*JVS( 127 )
  a = -W( 26 ) / JVS(          128  )
  W( 26 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 129 )
  W( 63 ) = W( 63 ) + a*JVS( 130 )
  W( 66 ) = W( 66 ) + a*JVS( 131 )
  W( 67 ) = W( 67 ) + a*JVS( 132 )
  W( 68 ) = W( 68 ) + a*JVS( 133 )
  a = -W( 27 ) / JVS(          134  )
  W( 27 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 135 )
  W( 34 ) = W( 34 ) + a*JVS( 136 )
  W( 40 ) = W( 40 ) + a*JVS( 137 )
  W( 48 ) = W( 48 ) + a*JVS( 138 )
  W( 51 ) = W( 51 ) + a*JVS( 139 )
  W( 53 ) = W( 53 ) + a*JVS( 140 )
  W( 54 ) = W( 54 ) + a*JVS( 141 )
  W( 67 ) = W( 67 ) + a*JVS( 142 )
  W( 68 ) = W( 68 ) + a*JVS( 143 )
  W( 70 ) = W( 70 ) + a*JVS( 144 )
  W( 73 ) = W( 73 ) + a*JVS( 145 )
  a = -W( 28 ) / JVS(          148  )
  W( 28 ) = -a
  W( 29 ) = W( 29 ) + a*JVS( 149 )
  W( 30 ) = W( 30 ) + a*JVS( 150 )
  W( 33 ) = W( 33 ) + a*JVS( 151 )
  W( 34 ) = W( 34 ) + a*JVS( 152 )
  W( 35 ) = W( 35 ) + a*JVS( 153 )
  W( 37 ) = W( 37 ) + a*JVS( 154 )
  W( 40 ) = W( 40 ) + a*JVS( 155 )
  W( 48 ) = W( 48 ) + a*JVS( 156 )
  W( 51 ) = W( 51 ) + a*JVS( 157 )
  W( 53 ) = W( 53 ) + a*JVS( 158 )
  W( 54 ) = W( 54 ) + a*JVS( 159 )
  W( 57 ) = W( 57 ) + a*JVS( 160 )
  W( 58 ) = W( 58 ) + a*JVS( 161 )
  W( 62 ) = W( 62 ) + a*JVS( 162 )
  W( 63 ) = W( 63 ) + a*JVS( 163 )
  W( 68 ) = W( 68 ) + a*JVS( 164 )
  W( 70 ) = W( 70 ) + a*JVS( 165 )
  a = -W( 29 ) / JVS(          166  )
  W( 29 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 167 )
  W( 68 ) = W( 68 ) + a*JVS( 168 )
  W( 70 ) = W( 70 ) + a*JVS( 169 )
  a = -W( 30 ) / JVS(          170  )
  W( 30 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 171 )
  W( 68 ) = W( 68 ) + a*JVS( 172 )
  W( 70 ) = W( 70 ) + a*JVS( 173 )
  a = -W( 31 ) / JVS(          174  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 175 )
  W( 63 ) = W( 63 ) + a*JVS( 176 )
  W( 66 ) = W( 66 ) + a*JVS( 177 )
  W( 68 ) = W( 68 ) + a*JVS( 178 )
  W( 70 ) = W( 70 ) + a*JVS( 179 )
  W( 73 ) = W( 73 ) + a*JVS( 180 )
  a = -W( 32 ) / JVS(          185  )
  W( 32 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 186 )
  W( 67 ) = W( 67 ) + a*JVS( 187 )
  W( 68 ) = W( 68 ) + a*JVS( 188 )
  W( 70 ) = W( 70 ) + a*JVS( 189 )
  W( 73 ) = W( 73 ) + a*JVS( 190 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 34 ) / JVS(          196  )
  W( 34 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 197 )
  W( 46 ) = W( 46 ) + a*JVS( 198 )
  W( 47 ) = W( 47 ) + a*JVS( 199 )
  W( 52 ) = W( 52 ) + a*JVS( 200 )
  W( 54 ) = W( 54 ) + a*JVS( 201 )
  W( 63 ) = W( 63 ) + a*JVS( 202 )
  W( 65 ) = W( 65 ) + a*JVS( 203 )
  W( 66 ) = W( 66 ) + a*JVS( 204 )
  W( 68 ) = W( 68 ) + a*JVS( 205 )
  W( 70 ) = W( 70 ) + a*JVS( 206 )
  W( 72 ) = W( 72 ) + a*JVS( 207 )
  a = -W( 35 ) / JVS(          208  )
  W( 35 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 209 )
  W( 63 ) = W( 63 ) + a*JVS( 210 )
  W( 68 ) = W( 68 ) + a*JVS( 211 )
  W( 70 ) = W( 70 ) + a*JVS( 212 )
  W( 73 ) = W( 73 ) + a*JVS( 213 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 40 ) / JVS(          246  )
  W( 40 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 247 )
  W( 46 ) = W( 46 ) + a*JVS( 248 )
  W( 47 ) = W( 47 ) + a*JVS( 249 )
  W( 51 ) = W( 51 ) + a*JVS( 250 )
  W( 52 ) = W( 52 ) + a*JVS( 251 )
  W( 54 ) = W( 54 ) + a*JVS( 252 )
  W( 55 ) = W( 55 ) + a*JVS( 253 )
  W( 63 ) = W( 63 ) + a*JVS( 254 )
  W( 64 ) = W( 64 ) + a*JVS( 255 )
  W( 65 ) = W( 65 ) + a*JVS( 256 )
  W( 66 ) = W( 66 ) + a*JVS( 257 )
  W( 68 ) = W( 68 ) + a*JVS( 258 )
  W( 70 ) = W( 70 ) + a*JVS( 259 )
  W( 72 ) = W( 72 ) + a*JVS( 260 )
  W( 73 ) = W( 73 ) + a*JVS( 261 )
  a = -W( 41 ) / JVS(          264  )
  W( 41 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 265 )
  W( 65 ) = W( 65 ) + a*JVS( 266 )
  W( 66 ) = W( 66 ) + a*JVS( 267 )
  W( 67 ) = W( 67 ) + a*JVS( 268 )
  W( 68 ) = W( 68 ) + a*JVS( 269 )
  W( 70 ) = W( 70 ) + a*JVS( 270 )
  W( 72 ) = W( 72 ) + a*JVS( 271 )
  W( 73 ) = W( 73 ) + a*JVS( 272 )
  a = -W( 42 ) / JVS(          278  )
  W( 42 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 279 )
  W( 45 ) = W( 45 ) + a*JVS( 280 )
  W( 49 ) = W( 49 ) + a*JVS( 281 )
  W( 52 ) = W( 52 ) + a*JVS( 282 )
  W( 54 ) = W( 54 ) + a*JVS( 283 )
  W( 57 ) = W( 57 ) + a*JVS( 284 )
  W( 58 ) = W( 58 ) + a*JVS( 285 )
  W( 59 ) = W( 59 ) + a*JVS( 286 )
  W( 60 ) = W( 60 ) + a*JVS( 287 )
  W( 63 ) = W( 63 ) + a*JVS( 288 )
  W( 64 ) = W( 64 ) + a*JVS( 289 )
  W( 65 ) = W( 65 ) + a*JVS( 290 )
  W( 66 ) = W( 66 ) + a*JVS( 291 )
  W( 67 ) = W( 67 ) + a*JVS( 292 )
  W( 68 ) = W( 68 ) + a*JVS( 293 )
  W( 69 ) = W( 69 ) + a*JVS( 294 )
  W( 70 ) = W( 70 ) + a*JVS( 295 )
  W( 71 ) = W( 71 ) + a*JVS( 296 )
  W( 72 ) = W( 72 ) + a*JVS( 297 )
  a = -W( 44 ) / JVS(          308  )
  W( 44 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 309 )
  W( 66 ) = W( 66 ) + a*JVS( 310 )
  W( 67 ) = W( 67 ) + a*JVS( 311 )
  W( 68 ) = W( 68 ) + a*JVS( 312 )
  W( 70 ) = W( 70 ) + a*JVS( 313 )
  W( 72 ) = W( 72 ) + a*JVS( 314 )
  a = -W( 45 ) / JVS(          316  )
  W( 45 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 70 ) = W( 70 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 325 )
  W( 65 ) = W( 65 ) + a*JVS( 326 )
  W( 66 ) = W( 66 ) + a*JVS( 327 )
  W( 67 ) = W( 67 ) + a*JVS( 328 )
  W( 68 ) = W( 68 ) + a*JVS( 329 )
  W( 70 ) = W( 70 ) + a*JVS( 330 )
  W( 72 ) = W( 72 ) + a*JVS( 331 )
  W( 73 ) = W( 73 ) + a*JVS( 332 )
  a = -W( 47 ) / JVS(          334  )
  W( 47 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 335 )
  W( 65 ) = W( 65 ) + a*JVS( 336 )
  W( 66 ) = W( 66 ) + a*JVS( 337 )
  W( 67 ) = W( 67 ) + a*JVS( 338 )
  W( 68 ) = W( 68 ) + a*JVS( 339 )
  W( 70 ) = W( 70 ) + a*JVS( 340 )
  W( 72 ) = W( 72 ) + a*JVS( 341 )
  W( 73 ) = W( 73 ) + a*JVS( 342 )
  a = -W( 48 ) / JVS(          363  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 364 )
  W( 50 ) = W( 50 ) + a*JVS( 365 )
  W( 51 ) = W( 51 ) + a*JVS( 366 )
  W( 52 ) = W( 52 ) + a*JVS( 367 )
  W( 54 ) = W( 54 ) + a*JVS( 368 )
  W( 55 ) = W( 55 ) + a*JVS( 369 )
  W( 56 ) = W( 56 ) + a*JVS( 370 )
  W( 57 ) = W( 57 ) + a*JVS( 371 )
  W( 58 ) = W( 58 ) + a*JVS( 372 )
  W( 59 ) = W( 59 ) + a*JVS( 373 )
  W( 60 ) = W( 60 ) + a*JVS( 374 )
  W( 61 ) = W( 61 ) + a*JVS( 375 )
  W( 62 ) = W( 62 ) + a*JVS( 376 )
  W( 63 ) = W( 63 ) + a*JVS( 377 )
  W( 64 ) = W( 64 ) + a*JVS( 378 )
  W( 65 ) = W( 65 ) + a*JVS( 379 )
  W( 66 ) = W( 66 ) + a*JVS( 380 )
  W( 67 ) = W( 67 ) + a*JVS( 381 )
  W( 68 ) = W( 68 ) + a*JVS( 382 )
  W( 70 ) = W( 70 ) + a*JVS( 383 )
  W( 72 ) = W( 72 ) + a*JVS( 384 )
  W( 73 ) = W( 73 ) + a*JVS( 385 )
  a = -W( 49 ) / JVS(          387  )
  W( 49 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 388 )
  W( 65 ) = W( 65 ) + a*JVS( 389 )
  W( 66 ) = W( 66 ) + a*JVS( 390 )
  W( 67 ) = W( 67 ) + a*JVS( 391 )
  W( 68 ) = W( 68 ) + a*JVS( 392 )
  W( 70 ) = W( 70 ) + a*JVS( 393 )
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  a = -W( 50 ) / JVS(          397  )
  W( 50 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 398 )
  W( 63 ) = W( 63 ) + a*JVS( 399 )
  W( 65 ) = W( 65 ) + a*JVS( 400 )
  W( 66 ) = W( 66 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 70 ) = W( 70 ) + a*JVS( 404 )
  W( 72 ) = W( 72 ) + a*JVS( 405 )
  a = -W( 51 ) / JVS(          411  )
  W( 51 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 412 )
  W( 63 ) = W( 63 ) + a*JVS( 413 )
  W( 65 ) = W( 65 ) + a*JVS( 414 )
  W( 66 ) = W( 66 ) + a*JVS( 415 )
  W( 67 ) = W( 67 ) + a*JVS( 416 )
  W( 68 ) = W( 68 ) + a*JVS( 417 )
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 72 ) = W( 72 ) + a*JVS( 419 )
  a = -W( 52 ) / JVS(          421  )
  W( 52 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 422 )
  W( 66 ) = W( 66 ) + a*JVS( 423 )
  W( 67 ) = W( 67 ) + a*JVS( 424 )
  W( 68 ) = W( 68 ) + a*JVS( 425 )
  W( 69 ) = W( 69 ) + a*JVS( 426 )
  W( 70 ) = W( 70 ) + a*JVS( 427 )
  W( 71 ) = W( 71 ) + a*JVS( 428 )
  W( 72 ) = W( 72 ) + a*JVS( 429 )
  a = -W( 53 ) / JVS(          442  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 443 )
  W( 57 ) = W( 57 ) + a*JVS( 444 )
  W( 58 ) = W( 58 ) + a*JVS( 445 )
  W( 59 ) = W( 59 ) + a*JVS( 446 )
  W( 60 ) = W( 60 ) + a*JVS( 447 )
  W( 61 ) = W( 61 ) + a*JVS( 448 )
  W( 62 ) = W( 62 ) + a*JVS( 449 )
  W( 63 ) = W( 63 ) + a*JVS( 450 )
  W( 64 ) = W( 64 ) + a*JVS( 451 )
  W( 65 ) = W( 65 ) + a*JVS( 452 )
  W( 66 ) = W( 66 ) + a*JVS( 453 )
  W( 67 ) = W( 67 ) + a*JVS( 454 )
  W( 68 ) = W( 68 ) + a*JVS( 455 )
  W( 69 ) = W( 69 ) + a*JVS( 456 )
  W( 70 ) = W( 70 ) + a*JVS( 457 )
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 72 ) = W( 72 ) + a*JVS( 459 )
  a = -W( 54 ) / JVS(          463  )
  W( 54 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 464 )
  W( 63 ) = W( 63 ) + a*JVS( 465 )
  W( 65 ) = W( 65 ) + a*JVS( 466 )
  W( 66 ) = W( 66 ) + a*JVS( 467 )
  W( 67 ) = W( 67 ) + a*JVS( 468 )
  W( 68 ) = W( 68 ) + a*JVS( 469 )
  W( 70 ) = W( 70 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  a = -W( 55 ) / JVS(          476  )
  W( 55 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 477 )
  W( 63 ) = W( 63 ) + a*JVS( 478 )
  W( 65 ) = W( 65 ) + a*JVS( 479 )
  W( 66 ) = W( 66 ) + a*JVS( 480 )
  W( 67 ) = W( 67 ) + a*JVS( 481 )
  W( 68 ) = W( 68 ) + a*JVS( 482 )
  W( 70 ) = W( 70 ) + a*JVS( 483 )
  W( 72 ) = W( 72 ) + a*JVS( 484 )
  W( 73 ) = W( 73 ) + a*JVS( 485 )
  a = -W( 56 ) / JVS(          500  )
  W( 56 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 501 )
  W( 63 ) = W( 63 ) + a*JVS( 502 )
  W( 64 ) = W( 64 ) + a*JVS( 503 )
  W( 65 ) = W( 65 ) + a*JVS( 504 )
  W( 66 ) = W( 66 ) + a*JVS( 505 )
  W( 67 ) = W( 67 ) + a*JVS( 506 )
  W( 68 ) = W( 68 ) + a*JVS( 507 )
  W( 69 ) = W( 69 ) + a*JVS( 508 )
  W( 70 ) = W( 70 ) + a*JVS( 509 )
  W( 71 ) = W( 71 ) + a*JVS( 510 )
  W( 72 ) = W( 72 ) + a*JVS( 511 )
  W( 73 ) = W( 73 ) + a*JVS( 512 )
  a = -W( 57 ) / JVS(          517  )
  W( 57 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 518 )
  W( 63 ) = W( 63 ) + a*JVS( 519 )
  W( 65 ) = W( 65 ) + a*JVS( 520 )
  W( 66 ) = W( 66 ) + a*JVS( 521 )
  W( 67 ) = W( 67 ) + a*JVS( 522 )
  W( 68 ) = W( 68 ) + a*JVS( 523 )
  W( 70 ) = W( 70 ) + a*JVS( 524 )
  W( 72 ) = W( 72 ) + a*JVS( 525 )
  a = -W( 58 ) / JVS(          528  )
  W( 58 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 529 )
  W( 63 ) = W( 63 ) + a*JVS( 530 )
  W( 65 ) = W( 65 ) + a*JVS( 531 )
  W( 66 ) = W( 66 ) + a*JVS( 532 )
  W( 67 ) = W( 67 ) + a*JVS( 533 )
  W( 68 ) = W( 68 ) + a*JVS( 534 )
  W( 70 ) = W( 70 ) + a*JVS( 535 )
  W( 72 ) = W( 72 ) + a*JVS( 536 )
  a = -W( 59 ) / JVS(          545  )
  W( 59 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 546 )
  W( 62 ) = W( 62 ) + a*JVS( 547 )
  W( 63 ) = W( 63 ) + a*JVS( 548 )
  W( 65 ) = W( 65 ) + a*JVS( 549 )
  W( 66 ) = W( 66 ) + a*JVS( 550 )
  W( 67 ) = W( 67 ) + a*JVS( 551 )
  W( 68 ) = W( 68 ) + a*JVS( 552 )
  W( 70 ) = W( 70 ) + a*JVS( 553 )
  W( 72 ) = W( 72 ) + a*JVS( 554 )
  a = -W( 60 ) / JVS(          563  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 564 )
  W( 63 ) = W( 63 ) + a*JVS( 565 )
  W( 65 ) = W( 65 ) + a*JVS( 566 )
  W( 66 ) = W( 66 ) + a*JVS( 567 )
  W( 67 ) = W( 67 ) + a*JVS( 568 )
  W( 68 ) = W( 68 ) + a*JVS( 569 )
  W( 70 ) = W( 70 ) + a*JVS( 570 )
  W( 72 ) = W( 72 ) + a*JVS( 571 )
  a = -W( 61 ) / JVS(          585  )
  W( 61 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 586 )
  W( 63 ) = W( 63 ) + a*JVS( 587 )
  W( 64 ) = W( 64 ) + a*JVS( 588 )
  W( 65 ) = W( 65 ) + a*JVS( 589 )
  W( 66 ) = W( 66 ) + a*JVS( 590 )
  W( 67 ) = W( 67 ) + a*JVS( 591 )
  W( 68 ) = W( 68 ) + a*JVS( 592 )
  W( 69 ) = W( 69 ) + a*JVS( 593 )
  W( 70 ) = W( 70 ) + a*JVS( 594 )
  W( 71 ) = W( 71 ) + a*JVS( 595 )
  W( 72 ) = W( 72 ) + a*JVS( 596 )
  W( 73 ) = W( 73 ) + a*JVS( 597 )
  a = -W( 62 ) / JVS(          602  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 603 )
  W( 65 ) = W( 65 ) + a*JVS( 604 )
  W( 66 ) = W( 66 ) + a*JVS( 605 )
  W( 67 ) = W( 67 ) + a*JVS( 606 )
  W( 68 ) = W( 68 ) + a*JVS( 607 )
  W( 70 ) = W( 70 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  a = -W( 63 ) / JVS(          628  )
  W( 63 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 629 )
  W( 66 ) = W( 66 ) + a*JVS( 630 )
  W( 67 ) = W( 67 ) + a*JVS( 631 )
  W( 68 ) = W( 68 ) + a*JVS( 632 )
  W( 70 ) = W( 70 ) + a*JVS( 633 )
  W( 72 ) = W( 72 ) + a*JVS( 634 )
  W( 73 ) = W( 73 ) + a*JVS( 635 )
  a = -W( 64 ) / JVS(          652  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 653 )
  W( 66 ) = W( 66 ) + a*JVS( 654 )
  W( 67 ) = W( 67 ) + a*JVS( 655 )
  W( 68 ) = W( 68 ) + a*JVS( 656 )
  W( 69 ) = W( 69 ) + a*JVS( 657 )
  W( 70 ) = W( 70 ) + a*JVS( 658 )
  W( 71 ) = W( 71 ) + a*JVS( 659 )
  W( 72 ) = W( 72 ) + a*JVS( 660 )
  W( 73 ) = W( 73 ) + a*JVS( 661 )
  a = -W( 65 ) / JVS(          691  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 692 )
  W( 67 ) = W( 67 ) + a*JVS( 693 )
  W( 68 ) = W( 68 ) + a*JVS( 694 )
  W( 69 ) = W( 69 ) + a*JVS( 695 )
  W( 70 ) = W( 70 ) + a*JVS( 696 )
  W( 71 ) = W( 71 ) + a*JVS( 697 )
  W( 72 ) = W( 72 ) + a*JVS( 698 )
  W( 73 ) = W( 73 ) + a*JVS( 699 )
  a = -W( 66 ) / JVS(          734  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 735 )
  W( 68 ) = W( 68 ) + a*JVS( 736 )
  W( 69 ) = W( 69 ) + a*JVS( 737 )
  W( 70 ) = W( 70 ) + a*JVS( 738 )
  W( 71 ) = W( 71 ) + a*JVS( 739 )
  W( 72 ) = W( 72 ) + a*JVS( 740 )
  W( 73 ) = W( 73 ) + a*JVS( 741 )
  a = -W( 67 ) / JVS(          796  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 797 )
  W( 69 ) = W( 69 ) + a*JVS( 798 )
  W( 70 ) = W( 70 ) + a*JVS( 799 )
  W( 71 ) = W( 71 ) + a*JVS( 800 )
  W( 72 ) = W( 72 ) + a*JVS( 801 )
  W( 73 ) = W( 73 ) + a*JVS( 802 )
  JVS( 803) = W( 5 )
  JVS( 804) = W( 6 )
  JVS( 805) = W( 7 )
  JVS( 806) = W( 8 )
  JVS( 807) = W( 9 )
  JVS( 808) = W( 11 )
  JVS( 809) = W( 12 )
  JVS( 810) = W( 13 )
  JVS( 811) = W( 14 )
  JVS( 812) = W( 15 )
  JVS( 813) = W( 16 )
  JVS( 814) = W( 17 )
  JVS( 815) = W( 18 )
  JVS( 816) = W( 19 )
  JVS( 817) = W( 21 )
  JVS( 818) = W( 22 )
  JVS( 819) = W( 23 )
  JVS( 820) = W( 24 )
  JVS( 821) = W( 25 )
  JVS( 822) = W( 26 )
  JVS( 823) = W( 27 )
  JVS( 824) = W( 28 )
  JVS( 825) = W( 29 )
  JVS( 826) = W( 30 )
  JVS( 827) = W( 31 )
  JVS( 828) = W( 32 )
  JVS( 829) = W( 33 )
  JVS( 830) = W( 34 )
  JVS( 831) = W( 35 )
  JVS( 832) = W( 37 )
  JVS( 833) = W( 40 )
  JVS( 834) = W( 41 )
  JVS( 835) = W( 42 )
  JVS( 836) = W( 44 )
  JVS( 837) = W( 45 )
  JVS( 838) = W( 46 )
  JVS( 839) = W( 47 )
  JVS( 840) = W( 48 )
  JVS( 841) = W( 49 )
  JVS( 842) = W( 50 )
  JVS( 843) = W( 51 )
  JVS( 844) = W( 52 )
  JVS( 845) = W( 53 )
  JVS( 846) = W( 54 )
  JVS( 847) = W( 55 )
  JVS( 848) = W( 56 )
  JVS( 849) = W( 57 )
  JVS( 850) = W( 58 )
  JVS( 851) = W( 59 )
  JVS( 852) = W( 60 )
  JVS( 853) = W( 61 )
  JVS( 854) = W( 62 )
  JVS( 855) = W( 63 )
  JVS( 856) = W( 64 )
  JVS( 857) = W( 65 )
  JVS( 858) = W( 66 )
  JVS( 859) = W( 67 )
  JVS( 860) = W( 68 )
  JVS( 861) = W( 69 )
  JVS( 862) = W( 70 )
  JVS( 863) = W( 71 )
  JVS( 864) = W( 72 )
  JVS( 865) = W( 73 )
  IF ( ABS(  JVS( 886 )) < TINY(a) ) THEN
         IER = 69                                      
         RETURN
  END IF
   W( 20 ) = JVS( 866 )
   W( 32 ) = JVS( 867 )
   W( 35 ) = JVS( 868 )
   W( 43 ) = JVS( 869 )
   W( 44 ) = JVS( 870 )
   W( 45 ) = JVS( 871 )
   W( 46 ) = JVS( 872 )
   W( 47 ) = JVS( 873 )
   W( 49 ) = JVS( 874 )
   W( 50 ) = JVS( 875 )
   W( 52 ) = JVS( 876 )
   W( 55 ) = JVS( 877 )
   W( 59 ) = JVS( 878 )
   W( 60 ) = JVS( 879 )
   W( 62 ) = JVS( 880 )
   W( 63 ) = JVS( 881 )
   W( 65 ) = JVS( 882 )
   W( 66 ) = JVS( 883 )
   W( 67 ) = JVS( 884 )
   W( 68 ) = JVS( 885 )
   W( 69 ) = JVS( 886 )
   W( 70 ) = JVS( 887 )
   W( 71 ) = JVS( 888 )
   W( 72 ) = JVS( 889 )
   W( 73 ) = JVS( 890 )
  a = -W( 20 ) / JVS(           96  )
  W( 20 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 97 )
  W( 67 ) = W( 67 ) + a*JVS( 98 )
  W( 68 ) = W( 68 ) + a*JVS( 99 )
  W( 70 ) = W( 70 ) + a*JVS( 100 )
  W( 73 ) = W( 73 ) + a*JVS( 101 )
  a = -W( 32 ) / JVS(          185  )
  W( 32 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 186 )
  W( 67 ) = W( 67 ) + a*JVS( 187 )
  W( 68 ) = W( 68 ) + a*JVS( 188 )
  W( 70 ) = W( 70 ) + a*JVS( 189 )
  W( 73 ) = W( 73 ) + a*JVS( 190 )
  a = -W( 35 ) / JVS(          208  )
  W( 35 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 209 )
  W( 63 ) = W( 63 ) + a*JVS( 210 )
  W( 68 ) = W( 68 ) + a*JVS( 211 )
  W( 70 ) = W( 70 ) + a*JVS( 212 )
  W( 73 ) = W( 73 ) + a*JVS( 213 )
  a = -W( 43 ) / JVS(          299  )
  W( 43 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 300 )
  W( 65 ) = W( 65 ) + a*JVS( 301 )
  W( 66 ) = W( 66 ) + a*JVS( 302 )
  W( 67 ) = W( 67 ) + a*JVS( 303 )
  W( 68 ) = W( 68 ) + a*JVS( 304 )
  W( 70 ) = W( 70 ) + a*JVS( 305 )
  W( 72 ) = W( 72 ) + a*JVS( 306 )
  a = -W( 44 ) / JVS(          308  )
  W( 44 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 309 )
  W( 66 ) = W( 66 ) + a*JVS( 310 )
  W( 67 ) = W( 67 ) + a*JVS( 311 )
  W( 68 ) = W( 68 ) + a*JVS( 312 )
  W( 70 ) = W( 70 ) + a*JVS( 313 )
  W( 72 ) = W( 72 ) + a*JVS( 314 )
  a = -W( 45 ) / JVS(          316  )
  W( 45 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 70 ) = W( 70 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 325 )
  W( 65 ) = W( 65 ) + a*JVS( 326 )
  W( 66 ) = W( 66 ) + a*JVS( 327 )
  W( 67 ) = W( 67 ) + a*JVS( 328 )
  W( 68 ) = W( 68 ) + a*JVS( 329 )
  W( 70 ) = W( 70 ) + a*JVS( 330 )
  W( 72 ) = W( 72 ) + a*JVS( 331 )
  W( 73 ) = W( 73 ) + a*JVS( 332 )
  a = -W( 47 ) / JVS(          334  )
  W( 47 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 335 )
  W( 65 ) = W( 65 ) + a*JVS( 336 )
  W( 66 ) = W( 66 ) + a*JVS( 337 )
  W( 67 ) = W( 67 ) + a*JVS( 338 )
  W( 68 ) = W( 68 ) + a*JVS( 339 )
  W( 70 ) = W( 70 ) + a*JVS( 340 )
  W( 72 ) = W( 72 ) + a*JVS( 341 )
  W( 73 ) = W( 73 ) + a*JVS( 342 )
  a = -W( 49 ) / JVS(          387  )
  W( 49 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 388 )
  W( 65 ) = W( 65 ) + a*JVS( 389 )
  W( 66 ) = W( 66 ) + a*JVS( 390 )
  W( 67 ) = W( 67 ) + a*JVS( 391 )
  W( 68 ) = W( 68 ) + a*JVS( 392 )
  W( 70 ) = W( 70 ) + a*JVS( 393 )
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  a = -W( 50 ) / JVS(          397  )
  W( 50 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 398 )
  W( 63 ) = W( 63 ) + a*JVS( 399 )
  W( 65 ) = W( 65 ) + a*JVS( 400 )
  W( 66 ) = W( 66 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 70 ) = W( 70 ) + a*JVS( 404 )
  W( 72 ) = W( 72 ) + a*JVS( 405 )
  a = -W( 52 ) / JVS(          421  )
  W( 52 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 422 )
  W( 66 ) = W( 66 ) + a*JVS( 423 )
  W( 67 ) = W( 67 ) + a*JVS( 424 )
  W( 68 ) = W( 68 ) + a*JVS( 425 )
  W( 69 ) = W( 69 ) + a*JVS( 426 )
  W( 70 ) = W( 70 ) + a*JVS( 427 )
  W( 71 ) = W( 71 ) + a*JVS( 428 )
  W( 72 ) = W( 72 ) + a*JVS( 429 )
  a = -W( 55 ) / JVS(          476  )
  W( 55 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 477 )
  W( 63 ) = W( 63 ) + a*JVS( 478 )
  W( 65 ) = W( 65 ) + a*JVS( 479 )
  W( 66 ) = W( 66 ) + a*JVS( 480 )
  W( 67 ) = W( 67 ) + a*JVS( 481 )
  W( 68 ) = W( 68 ) + a*JVS( 482 )
  W( 70 ) = W( 70 ) + a*JVS( 483 )
  W( 72 ) = W( 72 ) + a*JVS( 484 )
  W( 73 ) = W( 73 ) + a*JVS( 485 )
  a = -W( 59 ) / JVS(          545  )
  W( 59 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 546 )
  W( 62 ) = W( 62 ) + a*JVS( 547 )
  W( 63 ) = W( 63 ) + a*JVS( 548 )
  W( 65 ) = W( 65 ) + a*JVS( 549 )
  W( 66 ) = W( 66 ) + a*JVS( 550 )
  W( 67 ) = W( 67 ) + a*JVS( 551 )
  W( 68 ) = W( 68 ) + a*JVS( 552 )
  W( 70 ) = W( 70 ) + a*JVS( 553 )
  W( 72 ) = W( 72 ) + a*JVS( 554 )
  a = -W( 60 ) / JVS(          563  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 564 )
  W( 63 ) = W( 63 ) + a*JVS( 565 )
  W( 65 ) = W( 65 ) + a*JVS( 566 )
  W( 66 ) = W( 66 ) + a*JVS( 567 )
  W( 67 ) = W( 67 ) + a*JVS( 568 )
  W( 68 ) = W( 68 ) + a*JVS( 569 )
  W( 70 ) = W( 70 ) + a*JVS( 570 )
  W( 72 ) = W( 72 ) + a*JVS( 571 )
  a = -W( 62 ) / JVS(          602  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 603 )
  W( 65 ) = W( 65 ) + a*JVS( 604 )
  W( 66 ) = W( 66 ) + a*JVS( 605 )
  W( 67 ) = W( 67 ) + a*JVS( 606 )
  W( 68 ) = W( 68 ) + a*JVS( 607 )
  W( 70 ) = W( 70 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  a = -W( 63 ) / JVS(          628  )
  W( 63 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 629 )
  W( 66 ) = W( 66 ) + a*JVS( 630 )
  W( 67 ) = W( 67 ) + a*JVS( 631 )
  W( 68 ) = W( 68 ) + a*JVS( 632 )
  W( 70 ) = W( 70 ) + a*JVS( 633 )
  W( 72 ) = W( 72 ) + a*JVS( 634 )
  W( 73 ) = W( 73 ) + a*JVS( 635 )
  a = -W( 65 ) / JVS(          691  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 692 )
  W( 67 ) = W( 67 ) + a*JVS( 693 )
  W( 68 ) = W( 68 ) + a*JVS( 694 )
  W( 69 ) = W( 69 ) + a*JVS( 695 )
  W( 70 ) = W( 70 ) + a*JVS( 696 )
  W( 71 ) = W( 71 ) + a*JVS( 697 )
  W( 72 ) = W( 72 ) + a*JVS( 698 )
  W( 73 ) = W( 73 ) + a*JVS( 699 )
  a = -W( 66 ) / JVS(          734  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 735 )
  W( 68 ) = W( 68 ) + a*JVS( 736 )
  W( 69 ) = W( 69 ) + a*JVS( 737 )
  W( 70 ) = W( 70 ) + a*JVS( 738 )
  W( 71 ) = W( 71 ) + a*JVS( 739 )
  W( 72 ) = W( 72 ) + a*JVS( 740 )
  W( 73 ) = W( 73 ) + a*JVS( 741 )
  a = -W( 67 ) / JVS(          796  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 797 )
  W( 69 ) = W( 69 ) + a*JVS( 798 )
  W( 70 ) = W( 70 ) + a*JVS( 799 )
  W( 71 ) = W( 71 ) + a*JVS( 800 )
  W( 72 ) = W( 72 ) + a*JVS( 801 )
  W( 73 ) = W( 73 ) + a*JVS( 802 )
  a = -W( 68 ) / JVS(          860  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 861 )
  W( 70 ) = W( 70 ) + a*JVS( 862 )
  W( 71 ) = W( 71 ) + a*JVS( 863 )
  W( 72 ) = W( 72 ) + a*JVS( 864 )
  W( 73 ) = W( 73 ) + a*JVS( 865 )
  JVS( 866) = W( 20 )
  JVS( 867) = W( 32 )
  JVS( 868) = W( 35 )
  JVS( 869) = W( 43 )
  JVS( 870) = W( 44 )
  JVS( 871) = W( 45 )
  JVS( 872) = W( 46 )
  JVS( 873) = W( 47 )
  JVS( 874) = W( 49 )
  JVS( 875) = W( 50 )
  JVS( 876) = W( 52 )
  JVS( 877) = W( 55 )
  JVS( 878) = W( 59 )
  JVS( 879) = W( 60 )
  JVS( 880) = W( 62 )
  JVS( 881) = W( 63 )
  JVS( 882) = W( 65 )
  JVS( 883) = W( 66 )
  JVS( 884) = W( 67 )
  JVS( 885) = W( 68 )
  JVS( 886) = W( 69 )
  JVS( 887) = W( 70 )
  JVS( 888) = W( 71 )
  JVS( 889) = W( 72 )
  JVS( 890) = W( 73 )
  IF ( ABS(  JVS( 935 )) < TINY(a) ) THEN
         IER = 70                                      
         RETURN
  END IF
   W( 10 ) = JVS( 891 )
   W( 16 ) = JVS( 892 )
   W( 24 ) = JVS( 893 )
   W( 27 ) = JVS( 894 )
   W( 29 ) = JVS( 895 )
   W( 30 ) = JVS( 896 )
   W( 31 ) = JVS( 897 )
   W( 32 ) = JVS( 898 )
   W( 33 ) = JVS( 899 )
   W( 34 ) = JVS( 900 )
   W( 35 ) = JVS( 901 )
   W( 36 ) = JVS( 902 )
   W( 37 ) = JVS( 903 )
   W( 38 ) = JVS( 904 )
   W( 39 ) = JVS( 905 )
   W( 40 ) = JVS( 906 )
   W( 41 ) = JVS( 907 )
   W( 43 ) = JVS( 908 )
   W( 44 ) = JVS( 909 )
   W( 45 ) = JVS( 910 )
   W( 46 ) = JVS( 911 )
   W( 47 ) = JVS( 912 )
   W( 48 ) = JVS( 913 )
   W( 49 ) = JVS( 914 )
   W( 50 ) = JVS( 915 )
   W( 51 ) = JVS( 916 )
   W( 52 ) = JVS( 917 )
   W( 53 ) = JVS( 918 )
   W( 54 ) = JVS( 919 )
   W( 55 ) = JVS( 920 )
   W( 56 ) = JVS( 921 )
   W( 57 ) = JVS( 922 )
   W( 58 ) = JVS( 923 )
   W( 59 ) = JVS( 924 )
   W( 60 ) = JVS( 925 )
   W( 61 ) = JVS( 926 )
   W( 62 ) = JVS( 927 )
   W( 63 ) = JVS( 928 )
   W( 64 ) = JVS( 929 )
   W( 65 ) = JVS( 930 )
   W( 66 ) = JVS( 931 )
   W( 67 ) = JVS( 932 )
   W( 68 ) = JVS( 933 )
   W( 69 ) = JVS( 934 )
   W( 70 ) = JVS( 935 )
   W( 71 ) = JVS( 936 )
   W( 72 ) = JVS( 937 )
   W( 73 ) = JVS( 938 )
  a = -W( 10 ) / JVS(           51  )
  W( 10 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 52 )
  W( 73 ) = W( 73 ) + a*JVS( 53 )
  a = -W( 16 ) / JVS(           71  )
  W( 16 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 72 )
  W( 68 ) = W( 68 ) + a*JVS( 73 )
  W( 73 ) = W( 73 ) + a*JVS( 74 )
  a = -W( 24 ) / JVS(          119  )
  W( 24 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 120 )
  W( 68 ) = W( 68 ) + a*JVS( 121 )
  W( 70 ) = W( 70 ) + a*JVS( 122 )
  a = -W( 27 ) / JVS(          134  )
  W( 27 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 135 )
  W( 34 ) = W( 34 ) + a*JVS( 136 )
  W( 40 ) = W( 40 ) + a*JVS( 137 )
  W( 48 ) = W( 48 ) + a*JVS( 138 )
  W( 51 ) = W( 51 ) + a*JVS( 139 )
  W( 53 ) = W( 53 ) + a*JVS( 140 )
  W( 54 ) = W( 54 ) + a*JVS( 141 )
  W( 67 ) = W( 67 ) + a*JVS( 142 )
  W( 68 ) = W( 68 ) + a*JVS( 143 )
  W( 70 ) = W( 70 ) + a*JVS( 144 )
  W( 73 ) = W( 73 ) + a*JVS( 145 )
  a = -W( 29 ) / JVS(          166  )
  W( 29 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 167 )
  W( 68 ) = W( 68 ) + a*JVS( 168 )
  W( 70 ) = W( 70 ) + a*JVS( 169 )
  a = -W( 30 ) / JVS(          170  )
  W( 30 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 171 )
  W( 68 ) = W( 68 ) + a*JVS( 172 )
  W( 70 ) = W( 70 ) + a*JVS( 173 )
  a = -W( 31 ) / JVS(          174  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 175 )
  W( 63 ) = W( 63 ) + a*JVS( 176 )
  W( 66 ) = W( 66 ) + a*JVS( 177 )
  W( 68 ) = W( 68 ) + a*JVS( 178 )
  W( 70 ) = W( 70 ) + a*JVS( 179 )
  W( 73 ) = W( 73 ) + a*JVS( 180 )
  a = -W( 32 ) / JVS(          185  )
  W( 32 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 186 )
  W( 67 ) = W( 67 ) + a*JVS( 187 )
  W( 68 ) = W( 68 ) + a*JVS( 188 )
  W( 70 ) = W( 70 ) + a*JVS( 189 )
  W( 73 ) = W( 73 ) + a*JVS( 190 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 192 )
  W( 68 ) = W( 68 ) + a*JVS( 193 )
  W( 70 ) = W( 70 ) + a*JVS( 194 )
  a = -W( 34 ) / JVS(          196  )
  W( 34 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 197 )
  W( 46 ) = W( 46 ) + a*JVS( 198 )
  W( 47 ) = W( 47 ) + a*JVS( 199 )
  W( 52 ) = W( 52 ) + a*JVS( 200 )
  W( 54 ) = W( 54 ) + a*JVS( 201 )
  W( 63 ) = W( 63 ) + a*JVS( 202 )
  W( 65 ) = W( 65 ) + a*JVS( 203 )
  W( 66 ) = W( 66 ) + a*JVS( 204 )
  W( 68 ) = W( 68 ) + a*JVS( 205 )
  W( 70 ) = W( 70 ) + a*JVS( 206 )
  W( 72 ) = W( 72 ) + a*JVS( 207 )
  a = -W( 35 ) / JVS(          208  )
  W( 35 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 209 )
  W( 63 ) = W( 63 ) + a*JVS( 210 )
  W( 68 ) = W( 68 ) + a*JVS( 211 )
  W( 70 ) = W( 70 ) + a*JVS( 212 )
  W( 73 ) = W( 73 ) + a*JVS( 213 )
  a = -W( 36 ) / JVS(          215  )
  W( 36 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 216 )
  W( 65 ) = W( 65 ) + a*JVS( 217 )
  W( 66 ) = W( 66 ) + a*JVS( 218 )
  W( 67 ) = W( 67 ) + a*JVS( 219 )
  W( 68 ) = W( 68 ) + a*JVS( 220 )
  W( 70 ) = W( 70 ) + a*JVS( 221 )
  W( 72 ) = W( 72 ) + a*JVS( 222 )
  a = -W( 37 ) / JVS(          223  )
  W( 37 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 224 )
  W( 63 ) = W( 63 ) + a*JVS( 225 )
  W( 68 ) = W( 68 ) + a*JVS( 226 )
  W( 70 ) = W( 70 ) + a*JVS( 227 )
  a = -W( 38 ) / JVS(          228  )
  W( 38 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 229 )
  W( 65 ) = W( 65 ) + a*JVS( 230 )
  W( 66 ) = W( 66 ) + a*JVS( 231 )
  W( 67 ) = W( 67 ) + a*JVS( 232 )
  W( 68 ) = W( 68 ) + a*JVS( 233 )
  W( 70 ) = W( 70 ) + a*JVS( 234 )
  W( 72 ) = W( 72 ) + a*JVS( 235 )
  a = -W( 39 ) / JVS(          236  )
  W( 39 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 237 )
  W( 65 ) = W( 65 ) + a*JVS( 238 )
  W( 66 ) = W( 66 ) + a*JVS( 239 )
  W( 67 ) = W( 67 ) + a*JVS( 240 )
  W( 68 ) = W( 68 ) + a*JVS( 241 )
  W( 70 ) = W( 70 ) + a*JVS( 242 )
  W( 72 ) = W( 72 ) + a*JVS( 243 )
  a = -W( 40 ) / JVS(          246  )
  W( 40 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 247 )
  W( 46 ) = W( 46 ) + a*JVS( 248 )
  W( 47 ) = W( 47 ) + a*JVS( 249 )
  W( 51 ) = W( 51 ) + a*JVS( 250 )
  W( 52 ) = W( 52 ) + a*JVS( 251 )
  W( 54 ) = W( 54 ) + a*JVS( 252 )
  W( 55 ) = W( 55 ) + a*JVS( 253 )
  W( 63 ) = W( 63 ) + a*JVS( 254 )
  W( 64 ) = W( 64 ) + a*JVS( 255 )
  W( 65 ) = W( 65 ) + a*JVS( 256 )
  W( 66 ) = W( 66 ) + a*JVS( 257 )
  W( 68 ) = W( 68 ) + a*JVS( 258 )
  W( 70 ) = W( 70 ) + a*JVS( 259 )
  W( 72 ) = W( 72 ) + a*JVS( 260 )
  W( 73 ) = W( 73 ) + a*JVS( 261 )
  a = -W( 41 ) / JVS(          264  )
  W( 41 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 265 )
  W( 65 ) = W( 65 ) + a*JVS( 266 )
  W( 66 ) = W( 66 ) + a*JVS( 267 )
  W( 67 ) = W( 67 ) + a*JVS( 268 )
  W( 68 ) = W( 68 ) + a*JVS( 269 )
  W( 70 ) = W( 70 ) + a*JVS( 270 )
  W( 72 ) = W( 72 ) + a*JVS( 271 )
  W( 73 ) = W( 73 ) + a*JVS( 272 )
  a = -W( 43 ) / JVS(          299  )
  W( 43 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 300 )
  W( 65 ) = W( 65 ) + a*JVS( 301 )
  W( 66 ) = W( 66 ) + a*JVS( 302 )
  W( 67 ) = W( 67 ) + a*JVS( 303 )
  W( 68 ) = W( 68 ) + a*JVS( 304 )
  W( 70 ) = W( 70 ) + a*JVS( 305 )
  W( 72 ) = W( 72 ) + a*JVS( 306 )
  a = -W( 44 ) / JVS(          308  )
  W( 44 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 309 )
  W( 66 ) = W( 66 ) + a*JVS( 310 )
  W( 67 ) = W( 67 ) + a*JVS( 311 )
  W( 68 ) = W( 68 ) + a*JVS( 312 )
  W( 70 ) = W( 70 ) + a*JVS( 313 )
  W( 72 ) = W( 72 ) + a*JVS( 314 )
  a = -W( 45 ) / JVS(          316  )
  W( 45 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 70 ) = W( 70 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 325 )
  W( 65 ) = W( 65 ) + a*JVS( 326 )
  W( 66 ) = W( 66 ) + a*JVS( 327 )
  W( 67 ) = W( 67 ) + a*JVS( 328 )
  W( 68 ) = W( 68 ) + a*JVS( 329 )
  W( 70 ) = W( 70 ) + a*JVS( 330 )
  W( 72 ) = W( 72 ) + a*JVS( 331 )
  W( 73 ) = W( 73 ) + a*JVS( 332 )
  a = -W( 47 ) / JVS(          334  )
  W( 47 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 335 )
  W( 65 ) = W( 65 ) + a*JVS( 336 )
  W( 66 ) = W( 66 ) + a*JVS( 337 )
  W( 67 ) = W( 67 ) + a*JVS( 338 )
  W( 68 ) = W( 68 ) + a*JVS( 339 )
  W( 70 ) = W( 70 ) + a*JVS( 340 )
  W( 72 ) = W( 72 ) + a*JVS( 341 )
  W( 73 ) = W( 73 ) + a*JVS( 342 )
  a = -W( 48 ) / JVS(          363  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 364 )
  W( 50 ) = W( 50 ) + a*JVS( 365 )
  W( 51 ) = W( 51 ) + a*JVS( 366 )
  W( 52 ) = W( 52 ) + a*JVS( 367 )
  W( 54 ) = W( 54 ) + a*JVS( 368 )
  W( 55 ) = W( 55 ) + a*JVS( 369 )
  W( 56 ) = W( 56 ) + a*JVS( 370 )
  W( 57 ) = W( 57 ) + a*JVS( 371 )
  W( 58 ) = W( 58 ) + a*JVS( 372 )
  W( 59 ) = W( 59 ) + a*JVS( 373 )
  W( 60 ) = W( 60 ) + a*JVS( 374 )
  W( 61 ) = W( 61 ) + a*JVS( 375 )
  W( 62 ) = W( 62 ) + a*JVS( 376 )
  W( 63 ) = W( 63 ) + a*JVS( 377 )
  W( 64 ) = W( 64 ) + a*JVS( 378 )
  W( 65 ) = W( 65 ) + a*JVS( 379 )
  W( 66 ) = W( 66 ) + a*JVS( 380 )
  W( 67 ) = W( 67 ) + a*JVS( 381 )
  W( 68 ) = W( 68 ) + a*JVS( 382 )
  W( 70 ) = W( 70 ) + a*JVS( 383 )
  W( 72 ) = W( 72 ) + a*JVS( 384 )
  W( 73 ) = W( 73 ) + a*JVS( 385 )
  a = -W( 49 ) / JVS(          387  )
  W( 49 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 388 )
  W( 65 ) = W( 65 ) + a*JVS( 389 )
  W( 66 ) = W( 66 ) + a*JVS( 390 )
  W( 67 ) = W( 67 ) + a*JVS( 391 )
  W( 68 ) = W( 68 ) + a*JVS( 392 )
  W( 70 ) = W( 70 ) + a*JVS( 393 )
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  a = -W( 50 ) / JVS(          397  )
  W( 50 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 398 )
  W( 63 ) = W( 63 ) + a*JVS( 399 )
  W( 65 ) = W( 65 ) + a*JVS( 400 )
  W( 66 ) = W( 66 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 70 ) = W( 70 ) + a*JVS( 404 )
  W( 72 ) = W( 72 ) + a*JVS( 405 )
  a = -W( 51 ) / JVS(          411  )
  W( 51 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 412 )
  W( 63 ) = W( 63 ) + a*JVS( 413 )
  W( 65 ) = W( 65 ) + a*JVS( 414 )
  W( 66 ) = W( 66 ) + a*JVS( 415 )
  W( 67 ) = W( 67 ) + a*JVS( 416 )
  W( 68 ) = W( 68 ) + a*JVS( 417 )
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 72 ) = W( 72 ) + a*JVS( 419 )
  a = -W( 52 ) / JVS(          421  )
  W( 52 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 422 )
  W( 66 ) = W( 66 ) + a*JVS( 423 )
  W( 67 ) = W( 67 ) + a*JVS( 424 )
  W( 68 ) = W( 68 ) + a*JVS( 425 )
  W( 69 ) = W( 69 ) + a*JVS( 426 )
  W( 70 ) = W( 70 ) + a*JVS( 427 )
  W( 71 ) = W( 71 ) + a*JVS( 428 )
  W( 72 ) = W( 72 ) + a*JVS( 429 )
  a = -W( 53 ) / JVS(          442  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 443 )
  W( 57 ) = W( 57 ) + a*JVS( 444 )
  W( 58 ) = W( 58 ) + a*JVS( 445 )
  W( 59 ) = W( 59 ) + a*JVS( 446 )
  W( 60 ) = W( 60 ) + a*JVS( 447 )
  W( 61 ) = W( 61 ) + a*JVS( 448 )
  W( 62 ) = W( 62 ) + a*JVS( 449 )
  W( 63 ) = W( 63 ) + a*JVS( 450 )
  W( 64 ) = W( 64 ) + a*JVS( 451 )
  W( 65 ) = W( 65 ) + a*JVS( 452 )
  W( 66 ) = W( 66 ) + a*JVS( 453 )
  W( 67 ) = W( 67 ) + a*JVS( 454 )
  W( 68 ) = W( 68 ) + a*JVS( 455 )
  W( 69 ) = W( 69 ) + a*JVS( 456 )
  W( 70 ) = W( 70 ) + a*JVS( 457 )
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 72 ) = W( 72 ) + a*JVS( 459 )
  a = -W( 54 ) / JVS(          463  )
  W( 54 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 464 )
  W( 63 ) = W( 63 ) + a*JVS( 465 )
  W( 65 ) = W( 65 ) + a*JVS( 466 )
  W( 66 ) = W( 66 ) + a*JVS( 467 )
  W( 67 ) = W( 67 ) + a*JVS( 468 )
  W( 68 ) = W( 68 ) + a*JVS( 469 )
  W( 70 ) = W( 70 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  a = -W( 55 ) / JVS(          476  )
  W( 55 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 477 )
  W( 63 ) = W( 63 ) + a*JVS( 478 )
  W( 65 ) = W( 65 ) + a*JVS( 479 )
  W( 66 ) = W( 66 ) + a*JVS( 480 )
  W( 67 ) = W( 67 ) + a*JVS( 481 )
  W( 68 ) = W( 68 ) + a*JVS( 482 )
  W( 70 ) = W( 70 ) + a*JVS( 483 )
  W( 72 ) = W( 72 ) + a*JVS( 484 )
  W( 73 ) = W( 73 ) + a*JVS( 485 )
  a = -W( 56 ) / JVS(          500  )
  W( 56 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 501 )
  W( 63 ) = W( 63 ) + a*JVS( 502 )
  W( 64 ) = W( 64 ) + a*JVS( 503 )
  W( 65 ) = W( 65 ) + a*JVS( 504 )
  W( 66 ) = W( 66 ) + a*JVS( 505 )
  W( 67 ) = W( 67 ) + a*JVS( 506 )
  W( 68 ) = W( 68 ) + a*JVS( 507 )
  W( 69 ) = W( 69 ) + a*JVS( 508 )
  W( 70 ) = W( 70 ) + a*JVS( 509 )
  W( 71 ) = W( 71 ) + a*JVS( 510 )
  W( 72 ) = W( 72 ) + a*JVS( 511 )
  W( 73 ) = W( 73 ) + a*JVS( 512 )
  a = -W( 57 ) / JVS(          517  )
  W( 57 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 518 )
  W( 63 ) = W( 63 ) + a*JVS( 519 )
  W( 65 ) = W( 65 ) + a*JVS( 520 )
  W( 66 ) = W( 66 ) + a*JVS( 521 )
  W( 67 ) = W( 67 ) + a*JVS( 522 )
  W( 68 ) = W( 68 ) + a*JVS( 523 )
  W( 70 ) = W( 70 ) + a*JVS( 524 )
  W( 72 ) = W( 72 ) + a*JVS( 525 )
  a = -W( 58 ) / JVS(          528  )
  W( 58 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 529 )
  W( 63 ) = W( 63 ) + a*JVS( 530 )
  W( 65 ) = W( 65 ) + a*JVS( 531 )
  W( 66 ) = W( 66 ) + a*JVS( 532 )
  W( 67 ) = W( 67 ) + a*JVS( 533 )
  W( 68 ) = W( 68 ) + a*JVS( 534 )
  W( 70 ) = W( 70 ) + a*JVS( 535 )
  W( 72 ) = W( 72 ) + a*JVS( 536 )
  a = -W( 59 ) / JVS(          545  )
  W( 59 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 546 )
  W( 62 ) = W( 62 ) + a*JVS( 547 )
  W( 63 ) = W( 63 ) + a*JVS( 548 )
  W( 65 ) = W( 65 ) + a*JVS( 549 )
  W( 66 ) = W( 66 ) + a*JVS( 550 )
  W( 67 ) = W( 67 ) + a*JVS( 551 )
  W( 68 ) = W( 68 ) + a*JVS( 552 )
  W( 70 ) = W( 70 ) + a*JVS( 553 )
  W( 72 ) = W( 72 ) + a*JVS( 554 )
  a = -W( 60 ) / JVS(          563  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 564 )
  W( 63 ) = W( 63 ) + a*JVS( 565 )
  W( 65 ) = W( 65 ) + a*JVS( 566 )
  W( 66 ) = W( 66 ) + a*JVS( 567 )
  W( 67 ) = W( 67 ) + a*JVS( 568 )
  W( 68 ) = W( 68 ) + a*JVS( 569 )
  W( 70 ) = W( 70 ) + a*JVS( 570 )
  W( 72 ) = W( 72 ) + a*JVS( 571 )
  a = -W( 61 ) / JVS(          585  )
  W( 61 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 586 )
  W( 63 ) = W( 63 ) + a*JVS( 587 )
  W( 64 ) = W( 64 ) + a*JVS( 588 )
  W( 65 ) = W( 65 ) + a*JVS( 589 )
  W( 66 ) = W( 66 ) + a*JVS( 590 )
  W( 67 ) = W( 67 ) + a*JVS( 591 )
  W( 68 ) = W( 68 ) + a*JVS( 592 )
  W( 69 ) = W( 69 ) + a*JVS( 593 )
  W( 70 ) = W( 70 ) + a*JVS( 594 )
  W( 71 ) = W( 71 ) + a*JVS( 595 )
  W( 72 ) = W( 72 ) + a*JVS( 596 )
  W( 73 ) = W( 73 ) + a*JVS( 597 )
  a = -W( 62 ) / JVS(          602  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 603 )
  W( 65 ) = W( 65 ) + a*JVS( 604 )
  W( 66 ) = W( 66 ) + a*JVS( 605 )
  W( 67 ) = W( 67 ) + a*JVS( 606 )
  W( 68 ) = W( 68 ) + a*JVS( 607 )
  W( 70 ) = W( 70 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  a = -W( 63 ) / JVS(          628  )
  W( 63 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 629 )
  W( 66 ) = W( 66 ) + a*JVS( 630 )
  W( 67 ) = W( 67 ) + a*JVS( 631 )
  W( 68 ) = W( 68 ) + a*JVS( 632 )
  W( 70 ) = W( 70 ) + a*JVS( 633 )
  W( 72 ) = W( 72 ) + a*JVS( 634 )
  W( 73 ) = W( 73 ) + a*JVS( 635 )
  a = -W( 64 ) / JVS(          652  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 653 )
  W( 66 ) = W( 66 ) + a*JVS( 654 )
  W( 67 ) = W( 67 ) + a*JVS( 655 )
  W( 68 ) = W( 68 ) + a*JVS( 656 )
  W( 69 ) = W( 69 ) + a*JVS( 657 )
  W( 70 ) = W( 70 ) + a*JVS( 658 )
  W( 71 ) = W( 71 ) + a*JVS( 659 )
  W( 72 ) = W( 72 ) + a*JVS( 660 )
  W( 73 ) = W( 73 ) + a*JVS( 661 )
  a = -W( 65 ) / JVS(          691  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 692 )
  W( 67 ) = W( 67 ) + a*JVS( 693 )
  W( 68 ) = W( 68 ) + a*JVS( 694 )
  W( 69 ) = W( 69 ) + a*JVS( 695 )
  W( 70 ) = W( 70 ) + a*JVS( 696 )
  W( 71 ) = W( 71 ) + a*JVS( 697 )
  W( 72 ) = W( 72 ) + a*JVS( 698 )
  W( 73 ) = W( 73 ) + a*JVS( 699 )
  a = -W( 66 ) / JVS(          734  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 735 )
  W( 68 ) = W( 68 ) + a*JVS( 736 )
  W( 69 ) = W( 69 ) + a*JVS( 737 )
  W( 70 ) = W( 70 ) + a*JVS( 738 )
  W( 71 ) = W( 71 ) + a*JVS( 739 )
  W( 72 ) = W( 72 ) + a*JVS( 740 )
  W( 73 ) = W( 73 ) + a*JVS( 741 )
  a = -W( 67 ) / JVS(          796  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 797 )
  W( 69 ) = W( 69 ) + a*JVS( 798 )
  W( 70 ) = W( 70 ) + a*JVS( 799 )
  W( 71 ) = W( 71 ) + a*JVS( 800 )
  W( 72 ) = W( 72 ) + a*JVS( 801 )
  W( 73 ) = W( 73 ) + a*JVS( 802 )
  a = -W( 68 ) / JVS(          860  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 861 )
  W( 70 ) = W( 70 ) + a*JVS( 862 )
  W( 71 ) = W( 71 ) + a*JVS( 863 )
  W( 72 ) = W( 72 ) + a*JVS( 864 )
  W( 73 ) = W( 73 ) + a*JVS( 865 )
  a = -W( 69 ) / JVS(          886  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 887 )
  W( 71 ) = W( 71 ) + a*JVS( 888 )
  W( 72 ) = W( 72 ) + a*JVS( 889 )
  W( 73 ) = W( 73 ) + a*JVS( 890 )
  JVS( 891) = W( 10 )
  JVS( 892) = W( 16 )
  JVS( 893) = W( 24 )
  JVS( 894) = W( 27 )
  JVS( 895) = W( 29 )
  JVS( 896) = W( 30 )
  JVS( 897) = W( 31 )
  JVS( 898) = W( 32 )
  JVS( 899) = W( 33 )
  JVS( 900) = W( 34 )
  JVS( 901) = W( 35 )
  JVS( 902) = W( 36 )
  JVS( 903) = W( 37 )
  JVS( 904) = W( 38 )
  JVS( 905) = W( 39 )
  JVS( 906) = W( 40 )
  JVS( 907) = W( 41 )
  JVS( 908) = W( 43 )
  JVS( 909) = W( 44 )
  JVS( 910) = W( 45 )
  JVS( 911) = W( 46 )
  JVS( 912) = W( 47 )
  JVS( 913) = W( 48 )
  JVS( 914) = W( 49 )
  JVS( 915) = W( 50 )
  JVS( 916) = W( 51 )
  JVS( 917) = W( 52 )
  JVS( 918) = W( 53 )
  JVS( 919) = W( 54 )
  JVS( 920) = W( 55 )
  JVS( 921) = W( 56 )
  JVS( 922) = W( 57 )
  JVS( 923) = W( 58 )
  JVS( 924) = W( 59 )
  JVS( 925) = W( 60 )
  JVS( 926) = W( 61 )
  JVS( 927) = W( 62 )
  JVS( 928) = W( 63 )
  JVS( 929) = W( 64 )
  JVS( 930) = W( 65 )
  JVS( 931) = W( 66 )
  JVS( 932) = W( 67 )
  JVS( 933) = W( 68 )
  JVS( 934) = W( 69 )
  JVS( 935) = W( 70 )
  JVS( 936) = W( 71 )
  JVS( 937) = W( 72 )
  JVS( 938) = W( 73 )
  IF ( ABS(  JVS( 966 )) < TINY(a) ) THEN
         IER = 71                                      
         RETURN
  END IF
   W( 36 ) = JVS( 939 )
   W( 38 ) = JVS( 940 )
   W( 39 ) = JVS( 941 )
   W( 41 ) = JVS( 942 )
   W( 43 ) = JVS( 943 )
   W( 44 ) = JVS( 944 )
   W( 45 ) = JVS( 945 )
   W( 46 ) = JVS( 946 )
   W( 47 ) = JVS( 947 )
   W( 49 ) = JVS( 948 )
   W( 50 ) = JVS( 949 )
   W( 51 ) = JVS( 950 )
   W( 52 ) = JVS( 951 )
   W( 55 ) = JVS( 952 )
   W( 56 ) = JVS( 953 )
   W( 57 ) = JVS( 954 )
   W( 58 ) = JVS( 955 )
   W( 61 ) = JVS( 956 )
   W( 62 ) = JVS( 957 )
   W( 63 ) = JVS( 958 )
   W( 64 ) = JVS( 959 )
   W( 65 ) = JVS( 960 )
   W( 66 ) = JVS( 961 )
   W( 67 ) = JVS( 962 )
   W( 68 ) = JVS( 963 )
   W( 69 ) = JVS( 964 )
   W( 70 ) = JVS( 965 )
   W( 71 ) = JVS( 966 )
   W( 72 ) = JVS( 967 )
   W( 73 ) = JVS( 968 )
  a = -W( 36 ) / JVS(          215  )
  W( 36 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 216 )
  W( 65 ) = W( 65 ) + a*JVS( 217 )
  W( 66 ) = W( 66 ) + a*JVS( 218 )
  W( 67 ) = W( 67 ) + a*JVS( 219 )
  W( 68 ) = W( 68 ) + a*JVS( 220 )
  W( 70 ) = W( 70 ) + a*JVS( 221 )
  W( 72 ) = W( 72 ) + a*JVS( 222 )
  a = -W( 38 ) / JVS(          228  )
  W( 38 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 229 )
  W( 65 ) = W( 65 ) + a*JVS( 230 )
  W( 66 ) = W( 66 ) + a*JVS( 231 )
  W( 67 ) = W( 67 ) + a*JVS( 232 )
  W( 68 ) = W( 68 ) + a*JVS( 233 )
  W( 70 ) = W( 70 ) + a*JVS( 234 )
  W( 72 ) = W( 72 ) + a*JVS( 235 )
  a = -W( 39 ) / JVS(          236  )
  W( 39 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 237 )
  W( 65 ) = W( 65 ) + a*JVS( 238 )
  W( 66 ) = W( 66 ) + a*JVS( 239 )
  W( 67 ) = W( 67 ) + a*JVS( 240 )
  W( 68 ) = W( 68 ) + a*JVS( 241 )
  W( 70 ) = W( 70 ) + a*JVS( 242 )
  W( 72 ) = W( 72 ) + a*JVS( 243 )
  a = -W( 41 ) / JVS(          264  )
  W( 41 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 265 )
  W( 65 ) = W( 65 ) + a*JVS( 266 )
  W( 66 ) = W( 66 ) + a*JVS( 267 )
  W( 67 ) = W( 67 ) + a*JVS( 268 )
  W( 68 ) = W( 68 ) + a*JVS( 269 )
  W( 70 ) = W( 70 ) + a*JVS( 270 )
  W( 72 ) = W( 72 ) + a*JVS( 271 )
  W( 73 ) = W( 73 ) + a*JVS( 272 )
  a = -W( 43 ) / JVS(          299  )
  W( 43 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 300 )
  W( 65 ) = W( 65 ) + a*JVS( 301 )
  W( 66 ) = W( 66 ) + a*JVS( 302 )
  W( 67 ) = W( 67 ) + a*JVS( 303 )
  W( 68 ) = W( 68 ) + a*JVS( 304 )
  W( 70 ) = W( 70 ) + a*JVS( 305 )
  W( 72 ) = W( 72 ) + a*JVS( 306 )
  a = -W( 44 ) / JVS(          308  )
  W( 44 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 309 )
  W( 66 ) = W( 66 ) + a*JVS( 310 )
  W( 67 ) = W( 67 ) + a*JVS( 311 )
  W( 68 ) = W( 68 ) + a*JVS( 312 )
  W( 70 ) = W( 70 ) + a*JVS( 313 )
  W( 72 ) = W( 72 ) + a*JVS( 314 )
  a = -W( 45 ) / JVS(          316  )
  W( 45 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 70 ) = W( 70 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 325 )
  W( 65 ) = W( 65 ) + a*JVS( 326 )
  W( 66 ) = W( 66 ) + a*JVS( 327 )
  W( 67 ) = W( 67 ) + a*JVS( 328 )
  W( 68 ) = W( 68 ) + a*JVS( 329 )
  W( 70 ) = W( 70 ) + a*JVS( 330 )
  W( 72 ) = W( 72 ) + a*JVS( 331 )
  W( 73 ) = W( 73 ) + a*JVS( 332 )
  a = -W( 47 ) / JVS(          334  )
  W( 47 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 335 )
  W( 65 ) = W( 65 ) + a*JVS( 336 )
  W( 66 ) = W( 66 ) + a*JVS( 337 )
  W( 67 ) = W( 67 ) + a*JVS( 338 )
  W( 68 ) = W( 68 ) + a*JVS( 339 )
  W( 70 ) = W( 70 ) + a*JVS( 340 )
  W( 72 ) = W( 72 ) + a*JVS( 341 )
  W( 73 ) = W( 73 ) + a*JVS( 342 )
  a = -W( 49 ) / JVS(          387  )
  W( 49 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 388 )
  W( 65 ) = W( 65 ) + a*JVS( 389 )
  W( 66 ) = W( 66 ) + a*JVS( 390 )
  W( 67 ) = W( 67 ) + a*JVS( 391 )
  W( 68 ) = W( 68 ) + a*JVS( 392 )
  W( 70 ) = W( 70 ) + a*JVS( 393 )
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  a = -W( 50 ) / JVS(          397  )
  W( 50 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 398 )
  W( 63 ) = W( 63 ) + a*JVS( 399 )
  W( 65 ) = W( 65 ) + a*JVS( 400 )
  W( 66 ) = W( 66 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 70 ) = W( 70 ) + a*JVS( 404 )
  W( 72 ) = W( 72 ) + a*JVS( 405 )
  a = -W( 51 ) / JVS(          411  )
  W( 51 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 412 )
  W( 63 ) = W( 63 ) + a*JVS( 413 )
  W( 65 ) = W( 65 ) + a*JVS( 414 )
  W( 66 ) = W( 66 ) + a*JVS( 415 )
  W( 67 ) = W( 67 ) + a*JVS( 416 )
  W( 68 ) = W( 68 ) + a*JVS( 417 )
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 72 ) = W( 72 ) + a*JVS( 419 )
  a = -W( 52 ) / JVS(          421  )
  W( 52 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 422 )
  W( 66 ) = W( 66 ) + a*JVS( 423 )
  W( 67 ) = W( 67 ) + a*JVS( 424 )
  W( 68 ) = W( 68 ) + a*JVS( 425 )
  W( 69 ) = W( 69 ) + a*JVS( 426 )
  W( 70 ) = W( 70 ) + a*JVS( 427 )
  W( 71 ) = W( 71 ) + a*JVS( 428 )
  W( 72 ) = W( 72 ) + a*JVS( 429 )
  a = -W( 55 ) / JVS(          476  )
  W( 55 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 477 )
  W( 63 ) = W( 63 ) + a*JVS( 478 )
  W( 65 ) = W( 65 ) + a*JVS( 479 )
  W( 66 ) = W( 66 ) + a*JVS( 480 )
  W( 67 ) = W( 67 ) + a*JVS( 481 )
  W( 68 ) = W( 68 ) + a*JVS( 482 )
  W( 70 ) = W( 70 ) + a*JVS( 483 )
  W( 72 ) = W( 72 ) + a*JVS( 484 )
  W( 73 ) = W( 73 ) + a*JVS( 485 )
  a = -W( 56 ) / JVS(          500  )
  W( 56 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 501 )
  W( 63 ) = W( 63 ) + a*JVS( 502 )
  W( 64 ) = W( 64 ) + a*JVS( 503 )
  W( 65 ) = W( 65 ) + a*JVS( 504 )
  W( 66 ) = W( 66 ) + a*JVS( 505 )
  W( 67 ) = W( 67 ) + a*JVS( 506 )
  W( 68 ) = W( 68 ) + a*JVS( 507 )
  W( 69 ) = W( 69 ) + a*JVS( 508 )
  W( 70 ) = W( 70 ) + a*JVS( 509 )
  W( 71 ) = W( 71 ) + a*JVS( 510 )
  W( 72 ) = W( 72 ) + a*JVS( 511 )
  W( 73 ) = W( 73 ) + a*JVS( 512 )
  a = -W( 57 ) / JVS(          517  )
  W( 57 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 518 )
  W( 63 ) = W( 63 ) + a*JVS( 519 )
  W( 65 ) = W( 65 ) + a*JVS( 520 )
  W( 66 ) = W( 66 ) + a*JVS( 521 )
  W( 67 ) = W( 67 ) + a*JVS( 522 )
  W( 68 ) = W( 68 ) + a*JVS( 523 )
  W( 70 ) = W( 70 ) + a*JVS( 524 )
  W( 72 ) = W( 72 ) + a*JVS( 525 )
  a = -W( 58 ) / JVS(          528  )
  W( 58 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 529 )
  W( 63 ) = W( 63 ) + a*JVS( 530 )
  W( 65 ) = W( 65 ) + a*JVS( 531 )
  W( 66 ) = W( 66 ) + a*JVS( 532 )
  W( 67 ) = W( 67 ) + a*JVS( 533 )
  W( 68 ) = W( 68 ) + a*JVS( 534 )
  W( 70 ) = W( 70 ) + a*JVS( 535 )
  W( 72 ) = W( 72 ) + a*JVS( 536 )
  a = -W( 61 ) / JVS(          585  )
  W( 61 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 586 )
  W( 63 ) = W( 63 ) + a*JVS( 587 )
  W( 64 ) = W( 64 ) + a*JVS( 588 )
  W( 65 ) = W( 65 ) + a*JVS( 589 )
  W( 66 ) = W( 66 ) + a*JVS( 590 )
  W( 67 ) = W( 67 ) + a*JVS( 591 )
  W( 68 ) = W( 68 ) + a*JVS( 592 )
  W( 69 ) = W( 69 ) + a*JVS( 593 )
  W( 70 ) = W( 70 ) + a*JVS( 594 )
  W( 71 ) = W( 71 ) + a*JVS( 595 )
  W( 72 ) = W( 72 ) + a*JVS( 596 )
  W( 73 ) = W( 73 ) + a*JVS( 597 )
  a = -W( 62 ) / JVS(          602  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 603 )
  W( 65 ) = W( 65 ) + a*JVS( 604 )
  W( 66 ) = W( 66 ) + a*JVS( 605 )
  W( 67 ) = W( 67 ) + a*JVS( 606 )
  W( 68 ) = W( 68 ) + a*JVS( 607 )
  W( 70 ) = W( 70 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  a = -W( 63 ) / JVS(          628  )
  W( 63 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 629 )
  W( 66 ) = W( 66 ) + a*JVS( 630 )
  W( 67 ) = W( 67 ) + a*JVS( 631 )
  W( 68 ) = W( 68 ) + a*JVS( 632 )
  W( 70 ) = W( 70 ) + a*JVS( 633 )
  W( 72 ) = W( 72 ) + a*JVS( 634 )
  W( 73 ) = W( 73 ) + a*JVS( 635 )
  a = -W( 64 ) / JVS(          652  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 653 )
  W( 66 ) = W( 66 ) + a*JVS( 654 )
  W( 67 ) = W( 67 ) + a*JVS( 655 )
  W( 68 ) = W( 68 ) + a*JVS( 656 )
  W( 69 ) = W( 69 ) + a*JVS( 657 )
  W( 70 ) = W( 70 ) + a*JVS( 658 )
  W( 71 ) = W( 71 ) + a*JVS( 659 )
  W( 72 ) = W( 72 ) + a*JVS( 660 )
  W( 73 ) = W( 73 ) + a*JVS( 661 )
  a = -W( 65 ) / JVS(          691  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 692 )
  W( 67 ) = W( 67 ) + a*JVS( 693 )
  W( 68 ) = W( 68 ) + a*JVS( 694 )
  W( 69 ) = W( 69 ) + a*JVS( 695 )
  W( 70 ) = W( 70 ) + a*JVS( 696 )
  W( 71 ) = W( 71 ) + a*JVS( 697 )
  W( 72 ) = W( 72 ) + a*JVS( 698 )
  W( 73 ) = W( 73 ) + a*JVS( 699 )
  a = -W( 66 ) / JVS(          734  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 735 )
  W( 68 ) = W( 68 ) + a*JVS( 736 )
  W( 69 ) = W( 69 ) + a*JVS( 737 )
  W( 70 ) = W( 70 ) + a*JVS( 738 )
  W( 71 ) = W( 71 ) + a*JVS( 739 )
  W( 72 ) = W( 72 ) + a*JVS( 740 )
  W( 73 ) = W( 73 ) + a*JVS( 741 )
  a = -W( 67 ) / JVS(          796  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 797 )
  W( 69 ) = W( 69 ) + a*JVS( 798 )
  W( 70 ) = W( 70 ) + a*JVS( 799 )
  W( 71 ) = W( 71 ) + a*JVS( 800 )
  W( 72 ) = W( 72 ) + a*JVS( 801 )
  W( 73 ) = W( 73 ) + a*JVS( 802 )
  a = -W( 68 ) / JVS(          860  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 861 )
  W( 70 ) = W( 70 ) + a*JVS( 862 )
  W( 71 ) = W( 71 ) + a*JVS( 863 )
  W( 72 ) = W( 72 ) + a*JVS( 864 )
  W( 73 ) = W( 73 ) + a*JVS( 865 )
  a = -W( 69 ) / JVS(          886  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 887 )
  W( 71 ) = W( 71 ) + a*JVS( 888 )
  W( 72 ) = W( 72 ) + a*JVS( 889 )
  W( 73 ) = W( 73 ) + a*JVS( 890 )
  a = -W( 70 ) / JVS(          935  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 936 )
  W( 72 ) = W( 72 ) + a*JVS( 937 )
  W( 73 ) = W( 73 ) + a*JVS( 938 )
  JVS( 939) = W( 36 )
  JVS( 940) = W( 38 )
  JVS( 941) = W( 39 )
  JVS( 942) = W( 41 )
  JVS( 943) = W( 43 )
  JVS( 944) = W( 44 )
  JVS( 945) = W( 45 )
  JVS( 946) = W( 46 )
  JVS( 947) = W( 47 )
  JVS( 948) = W( 49 )
  JVS( 949) = W( 50 )
  JVS( 950) = W( 51 )
  JVS( 951) = W( 52 )
  JVS( 952) = W( 55 )
  JVS( 953) = W( 56 )
  JVS( 954) = W( 57 )
  JVS( 955) = W( 58 )
  JVS( 956) = W( 61 )
  JVS( 957) = W( 62 )
  JVS( 958) = W( 63 )
  JVS( 959) = W( 64 )
  JVS( 960) = W( 65 )
  JVS( 961) = W( 66 )
  JVS( 962) = W( 67 )
  JVS( 963) = W( 68 )
  JVS( 964) = W( 69 )
  JVS( 965) = W( 70 )
  JVS( 966) = W( 71 )
  JVS( 967) = W( 72 )
  JVS( 968) = W( 73 )
  IF ( ABS(  JVS( 1003 )) < TINY(a) ) THEN
         IER = 72                                      
         RETURN
  END IF
   W( 18 ) = JVS( 969 )
   W( 21 ) = JVS( 970 )
   W( 22 ) = JVS( 971 )
   W( 25 ) = JVS( 972 )
   W( 32 ) = JVS( 973 )
   W( 36 ) = JVS( 974 )
   W( 38 ) = JVS( 975 )
   W( 39 ) = JVS( 976 )
   W( 41 ) = JVS( 977 )
   W( 43 ) = JVS( 978 )
   W( 44 ) = JVS( 979 )
   W( 45 ) = JVS( 980 )
   W( 46 ) = JVS( 981 )
   W( 47 ) = JVS( 982 )
   W( 49 ) = JVS( 983 )
   W( 50 ) = JVS( 984 )
   W( 52 ) = JVS( 985 )
   W( 55 ) = JVS( 986 )
   W( 56 ) = JVS( 987 )
   W( 57 ) = JVS( 988 )
   W( 58 ) = JVS( 989 )
   W( 59 ) = JVS( 990 )
   W( 60 ) = JVS( 991 )
   W( 61 ) = JVS( 992 )
   W( 62 ) = JVS( 993 )
   W( 63 ) = JVS( 994 )
   W( 64 ) = JVS( 995 )
   W( 65 ) = JVS( 996 )
   W( 66 ) = JVS( 997 )
   W( 67 ) = JVS( 998 )
   W( 68 ) = JVS( 999 )
   W( 69 ) = JVS( 1000 )
   W( 70 ) = JVS( 1001 )
   W( 71 ) = JVS( 1002 )
   W( 72 ) = JVS( 1003 )
   W( 73 ) = JVS( 1004 )
  a = -W( 18 ) / JVS(           79  )
  W( 18 ) = -a
  W( 21 ) = W( 21 ) + a*JVS( 80 )
  W( 22 ) = W( 22 ) + a*JVS( 81 )
  W( 25 ) = W( 25 ) + a*JVS( 82 )
  W( 68 ) = W( 68 ) + a*JVS( 83 )
  W( 72 ) = W( 72 ) + a*JVS( 84 )
  W( 73 ) = W( 73 ) + a*JVS( 85 )
  a = -W( 21 ) / JVS(          103  )
  W( 21 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 104 )
  W( 68 ) = W( 68 ) + a*JVS( 105 )
  W( 73 ) = W( 73 ) + a*JVS( 106 )
  a = -W( 22 ) / JVS(          108  )
  W( 22 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 109 )
  W( 68 ) = W( 68 ) + a*JVS( 110 )
  W( 73 ) = W( 73 ) + a*JVS( 111 )
  a = -W( 25 ) / JVS(          123  )
  W( 25 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 124 )
  W( 63 ) = W( 63 ) + a*JVS( 125 )
  W( 68 ) = W( 68 ) + a*JVS( 126 )
  W( 73 ) = W( 73 ) + a*JVS( 127 )
  a = -W( 32 ) / JVS(          185  )
  W( 32 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 186 )
  W( 67 ) = W( 67 ) + a*JVS( 187 )
  W( 68 ) = W( 68 ) + a*JVS( 188 )
  W( 70 ) = W( 70 ) + a*JVS( 189 )
  W( 73 ) = W( 73 ) + a*JVS( 190 )
  a = -W( 36 ) / JVS(          215  )
  W( 36 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 216 )
  W( 65 ) = W( 65 ) + a*JVS( 217 )
  W( 66 ) = W( 66 ) + a*JVS( 218 )
  W( 67 ) = W( 67 ) + a*JVS( 219 )
  W( 68 ) = W( 68 ) + a*JVS( 220 )
  W( 70 ) = W( 70 ) + a*JVS( 221 )
  W( 72 ) = W( 72 ) + a*JVS( 222 )
  a = -W( 38 ) / JVS(          228  )
  W( 38 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 229 )
  W( 65 ) = W( 65 ) + a*JVS( 230 )
  W( 66 ) = W( 66 ) + a*JVS( 231 )
  W( 67 ) = W( 67 ) + a*JVS( 232 )
  W( 68 ) = W( 68 ) + a*JVS( 233 )
  W( 70 ) = W( 70 ) + a*JVS( 234 )
  W( 72 ) = W( 72 ) + a*JVS( 235 )
  a = -W( 39 ) / JVS(          236  )
  W( 39 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 237 )
  W( 65 ) = W( 65 ) + a*JVS( 238 )
  W( 66 ) = W( 66 ) + a*JVS( 239 )
  W( 67 ) = W( 67 ) + a*JVS( 240 )
  W( 68 ) = W( 68 ) + a*JVS( 241 )
  W( 70 ) = W( 70 ) + a*JVS( 242 )
  W( 72 ) = W( 72 ) + a*JVS( 243 )
  a = -W( 41 ) / JVS(          264  )
  W( 41 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 265 )
  W( 65 ) = W( 65 ) + a*JVS( 266 )
  W( 66 ) = W( 66 ) + a*JVS( 267 )
  W( 67 ) = W( 67 ) + a*JVS( 268 )
  W( 68 ) = W( 68 ) + a*JVS( 269 )
  W( 70 ) = W( 70 ) + a*JVS( 270 )
  W( 72 ) = W( 72 ) + a*JVS( 271 )
  W( 73 ) = W( 73 ) + a*JVS( 272 )
  a = -W( 43 ) / JVS(          299  )
  W( 43 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 300 )
  W( 65 ) = W( 65 ) + a*JVS( 301 )
  W( 66 ) = W( 66 ) + a*JVS( 302 )
  W( 67 ) = W( 67 ) + a*JVS( 303 )
  W( 68 ) = W( 68 ) + a*JVS( 304 )
  W( 70 ) = W( 70 ) + a*JVS( 305 )
  W( 72 ) = W( 72 ) + a*JVS( 306 )
  a = -W( 44 ) / JVS(          308  )
  W( 44 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 309 )
  W( 66 ) = W( 66 ) + a*JVS( 310 )
  W( 67 ) = W( 67 ) + a*JVS( 311 )
  W( 68 ) = W( 68 ) + a*JVS( 312 )
  W( 70 ) = W( 70 ) + a*JVS( 313 )
  W( 72 ) = W( 72 ) + a*JVS( 314 )
  a = -W( 45 ) / JVS(          316  )
  W( 45 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 70 ) = W( 70 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 325 )
  W( 65 ) = W( 65 ) + a*JVS( 326 )
  W( 66 ) = W( 66 ) + a*JVS( 327 )
  W( 67 ) = W( 67 ) + a*JVS( 328 )
  W( 68 ) = W( 68 ) + a*JVS( 329 )
  W( 70 ) = W( 70 ) + a*JVS( 330 )
  W( 72 ) = W( 72 ) + a*JVS( 331 )
  W( 73 ) = W( 73 ) + a*JVS( 332 )
  a = -W( 47 ) / JVS(          334  )
  W( 47 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 335 )
  W( 65 ) = W( 65 ) + a*JVS( 336 )
  W( 66 ) = W( 66 ) + a*JVS( 337 )
  W( 67 ) = W( 67 ) + a*JVS( 338 )
  W( 68 ) = W( 68 ) + a*JVS( 339 )
  W( 70 ) = W( 70 ) + a*JVS( 340 )
  W( 72 ) = W( 72 ) + a*JVS( 341 )
  W( 73 ) = W( 73 ) + a*JVS( 342 )
  a = -W( 49 ) / JVS(          387  )
  W( 49 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 388 )
  W( 65 ) = W( 65 ) + a*JVS( 389 )
  W( 66 ) = W( 66 ) + a*JVS( 390 )
  W( 67 ) = W( 67 ) + a*JVS( 391 )
  W( 68 ) = W( 68 ) + a*JVS( 392 )
  W( 70 ) = W( 70 ) + a*JVS( 393 )
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  a = -W( 50 ) / JVS(          397  )
  W( 50 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 398 )
  W( 63 ) = W( 63 ) + a*JVS( 399 )
  W( 65 ) = W( 65 ) + a*JVS( 400 )
  W( 66 ) = W( 66 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 70 ) = W( 70 ) + a*JVS( 404 )
  W( 72 ) = W( 72 ) + a*JVS( 405 )
  a = -W( 52 ) / JVS(          421  )
  W( 52 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 422 )
  W( 66 ) = W( 66 ) + a*JVS( 423 )
  W( 67 ) = W( 67 ) + a*JVS( 424 )
  W( 68 ) = W( 68 ) + a*JVS( 425 )
  W( 69 ) = W( 69 ) + a*JVS( 426 )
  W( 70 ) = W( 70 ) + a*JVS( 427 )
  W( 71 ) = W( 71 ) + a*JVS( 428 )
  W( 72 ) = W( 72 ) + a*JVS( 429 )
  a = -W( 55 ) / JVS(          476  )
  W( 55 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 477 )
  W( 63 ) = W( 63 ) + a*JVS( 478 )
  W( 65 ) = W( 65 ) + a*JVS( 479 )
  W( 66 ) = W( 66 ) + a*JVS( 480 )
  W( 67 ) = W( 67 ) + a*JVS( 481 )
  W( 68 ) = W( 68 ) + a*JVS( 482 )
  W( 70 ) = W( 70 ) + a*JVS( 483 )
  W( 72 ) = W( 72 ) + a*JVS( 484 )
  W( 73 ) = W( 73 ) + a*JVS( 485 )
  a = -W( 56 ) / JVS(          500  )
  W( 56 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 501 )
  W( 63 ) = W( 63 ) + a*JVS( 502 )
  W( 64 ) = W( 64 ) + a*JVS( 503 )
  W( 65 ) = W( 65 ) + a*JVS( 504 )
  W( 66 ) = W( 66 ) + a*JVS( 505 )
  W( 67 ) = W( 67 ) + a*JVS( 506 )
  W( 68 ) = W( 68 ) + a*JVS( 507 )
  W( 69 ) = W( 69 ) + a*JVS( 508 )
  W( 70 ) = W( 70 ) + a*JVS( 509 )
  W( 71 ) = W( 71 ) + a*JVS( 510 )
  W( 72 ) = W( 72 ) + a*JVS( 511 )
  W( 73 ) = W( 73 ) + a*JVS( 512 )
  a = -W( 57 ) / JVS(          517  )
  W( 57 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 518 )
  W( 63 ) = W( 63 ) + a*JVS( 519 )
  W( 65 ) = W( 65 ) + a*JVS( 520 )
  W( 66 ) = W( 66 ) + a*JVS( 521 )
  W( 67 ) = W( 67 ) + a*JVS( 522 )
  W( 68 ) = W( 68 ) + a*JVS( 523 )
  W( 70 ) = W( 70 ) + a*JVS( 524 )
  W( 72 ) = W( 72 ) + a*JVS( 525 )
  a = -W( 58 ) / JVS(          528  )
  W( 58 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 529 )
  W( 63 ) = W( 63 ) + a*JVS( 530 )
  W( 65 ) = W( 65 ) + a*JVS( 531 )
  W( 66 ) = W( 66 ) + a*JVS( 532 )
  W( 67 ) = W( 67 ) + a*JVS( 533 )
  W( 68 ) = W( 68 ) + a*JVS( 534 )
  W( 70 ) = W( 70 ) + a*JVS( 535 )
  W( 72 ) = W( 72 ) + a*JVS( 536 )
  a = -W( 59 ) / JVS(          545  )
  W( 59 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 546 )
  W( 62 ) = W( 62 ) + a*JVS( 547 )
  W( 63 ) = W( 63 ) + a*JVS( 548 )
  W( 65 ) = W( 65 ) + a*JVS( 549 )
  W( 66 ) = W( 66 ) + a*JVS( 550 )
  W( 67 ) = W( 67 ) + a*JVS( 551 )
  W( 68 ) = W( 68 ) + a*JVS( 552 )
  W( 70 ) = W( 70 ) + a*JVS( 553 )
  W( 72 ) = W( 72 ) + a*JVS( 554 )
  a = -W( 60 ) / JVS(          563  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 564 )
  W( 63 ) = W( 63 ) + a*JVS( 565 )
  W( 65 ) = W( 65 ) + a*JVS( 566 )
  W( 66 ) = W( 66 ) + a*JVS( 567 )
  W( 67 ) = W( 67 ) + a*JVS( 568 )
  W( 68 ) = W( 68 ) + a*JVS( 569 )
  W( 70 ) = W( 70 ) + a*JVS( 570 )
  W( 72 ) = W( 72 ) + a*JVS( 571 )
  a = -W( 61 ) / JVS(          585  )
  W( 61 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 586 )
  W( 63 ) = W( 63 ) + a*JVS( 587 )
  W( 64 ) = W( 64 ) + a*JVS( 588 )
  W( 65 ) = W( 65 ) + a*JVS( 589 )
  W( 66 ) = W( 66 ) + a*JVS( 590 )
  W( 67 ) = W( 67 ) + a*JVS( 591 )
  W( 68 ) = W( 68 ) + a*JVS( 592 )
  W( 69 ) = W( 69 ) + a*JVS( 593 )
  W( 70 ) = W( 70 ) + a*JVS( 594 )
  W( 71 ) = W( 71 ) + a*JVS( 595 )
  W( 72 ) = W( 72 ) + a*JVS( 596 )
  W( 73 ) = W( 73 ) + a*JVS( 597 )
  a = -W( 62 ) / JVS(          602  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 603 )
  W( 65 ) = W( 65 ) + a*JVS( 604 )
  W( 66 ) = W( 66 ) + a*JVS( 605 )
  W( 67 ) = W( 67 ) + a*JVS( 606 )
  W( 68 ) = W( 68 ) + a*JVS( 607 )
  W( 70 ) = W( 70 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  a = -W( 63 ) / JVS(          628  )
  W( 63 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 629 )
  W( 66 ) = W( 66 ) + a*JVS( 630 )
  W( 67 ) = W( 67 ) + a*JVS( 631 )
  W( 68 ) = W( 68 ) + a*JVS( 632 )
  W( 70 ) = W( 70 ) + a*JVS( 633 )
  W( 72 ) = W( 72 ) + a*JVS( 634 )
  W( 73 ) = W( 73 ) + a*JVS( 635 )
  a = -W( 64 ) / JVS(          652  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 653 )
  W( 66 ) = W( 66 ) + a*JVS( 654 )
  W( 67 ) = W( 67 ) + a*JVS( 655 )
  W( 68 ) = W( 68 ) + a*JVS( 656 )
  W( 69 ) = W( 69 ) + a*JVS( 657 )
  W( 70 ) = W( 70 ) + a*JVS( 658 )
  W( 71 ) = W( 71 ) + a*JVS( 659 )
  W( 72 ) = W( 72 ) + a*JVS( 660 )
  W( 73 ) = W( 73 ) + a*JVS( 661 )
  a = -W( 65 ) / JVS(          691  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 692 )
  W( 67 ) = W( 67 ) + a*JVS( 693 )
  W( 68 ) = W( 68 ) + a*JVS( 694 )
  W( 69 ) = W( 69 ) + a*JVS( 695 )
  W( 70 ) = W( 70 ) + a*JVS( 696 )
  W( 71 ) = W( 71 ) + a*JVS( 697 )
  W( 72 ) = W( 72 ) + a*JVS( 698 )
  W( 73 ) = W( 73 ) + a*JVS( 699 )
  a = -W( 66 ) / JVS(          734  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 735 )
  W( 68 ) = W( 68 ) + a*JVS( 736 )
  W( 69 ) = W( 69 ) + a*JVS( 737 )
  W( 70 ) = W( 70 ) + a*JVS( 738 )
  W( 71 ) = W( 71 ) + a*JVS( 739 )
  W( 72 ) = W( 72 ) + a*JVS( 740 )
  W( 73 ) = W( 73 ) + a*JVS( 741 )
  a = -W( 67 ) / JVS(          796  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 797 )
  W( 69 ) = W( 69 ) + a*JVS( 798 )
  W( 70 ) = W( 70 ) + a*JVS( 799 )
  W( 71 ) = W( 71 ) + a*JVS( 800 )
  W( 72 ) = W( 72 ) + a*JVS( 801 )
  W( 73 ) = W( 73 ) + a*JVS( 802 )
  a = -W( 68 ) / JVS(          860  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 861 )
  W( 70 ) = W( 70 ) + a*JVS( 862 )
  W( 71 ) = W( 71 ) + a*JVS( 863 )
  W( 72 ) = W( 72 ) + a*JVS( 864 )
  W( 73 ) = W( 73 ) + a*JVS( 865 )
  a = -W( 69 ) / JVS(          886  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 887 )
  W( 71 ) = W( 71 ) + a*JVS( 888 )
  W( 72 ) = W( 72 ) + a*JVS( 889 )
  W( 73 ) = W( 73 ) + a*JVS( 890 )
  a = -W( 70 ) / JVS(          935  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 936 )
  W( 72 ) = W( 72 ) + a*JVS( 937 )
  W( 73 ) = W( 73 ) + a*JVS( 938 )
  a = -W( 71 ) / JVS(          966  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 967 )
  W( 73 ) = W( 73 ) + a*JVS( 968 )
  JVS( 969) = W( 18 )
  JVS( 970) = W( 21 )
  JVS( 971) = W( 22 )
  JVS( 972) = W( 25 )
  JVS( 973) = W( 32 )
  JVS( 974) = W( 36 )
  JVS( 975) = W( 38 )
  JVS( 976) = W( 39 )
  JVS( 977) = W( 41 )
  JVS( 978) = W( 43 )
  JVS( 979) = W( 44 )
  JVS( 980) = W( 45 )
  JVS( 981) = W( 46 )
  JVS( 982) = W( 47 )
  JVS( 983) = W( 49 )
  JVS( 984) = W( 50 )
  JVS( 985) = W( 52 )
  JVS( 986) = W( 55 )
  JVS( 987) = W( 56 )
  JVS( 988) = W( 57 )
  JVS( 989) = W( 58 )
  JVS( 990) = W( 59 )
  JVS( 991) = W( 60 )
  JVS( 992) = W( 61 )
  JVS( 993) = W( 62 )
  JVS( 994) = W( 63 )
  JVS( 995) = W( 64 )
  JVS( 996) = W( 65 )
  JVS( 997) = W( 66 )
  JVS( 998) = W( 67 )
  JVS( 999) = W( 68 )
  JVS( 1000) = W( 69 )
  JVS( 1001) = W( 70 )
  JVS( 1002) = W( 71 )
  JVS( 1003) = W( 72 )
  JVS( 1004) = W( 73 )
  IF ( ABS(  JVS( 1052 )) < TINY(a) ) THEN
         IER = 73                                      
         RETURN
  END IF
   W( 10 ) = JVS( 1005 )
   W( 16 ) = JVS( 1006 )
   W( 18 ) = JVS( 1007 )
   W( 20 ) = JVS( 1008 )
   W( 21 ) = JVS( 1009 )
   W( 22 ) = JVS( 1010 )
   W( 25 ) = JVS( 1011 )
   W( 27 ) = JVS( 1012 )
   W( 31 ) = JVS( 1013 )
   W( 32 ) = JVS( 1014 )
   W( 34 ) = JVS( 1015 )
   W( 35 ) = JVS( 1016 )
   W( 36 ) = JVS( 1017 )
   W( 38 ) = JVS( 1018 )
   W( 39 ) = JVS( 1019 )
   W( 40 ) = JVS( 1020 )
   W( 41 ) = JVS( 1021 )
   W( 43 ) = JVS( 1022 )
   W( 44 ) = JVS( 1023 )
   W( 45 ) = JVS( 1024 )
   W( 46 ) = JVS( 1025 )
   W( 47 ) = JVS( 1026 )
   W( 48 ) = JVS( 1027 )
   W( 49 ) = JVS( 1028 )
   W( 50 ) = JVS( 1029 )
   W( 51 ) = JVS( 1030 )
   W( 52 ) = JVS( 1031 )
   W( 53 ) = JVS( 1032 )
   W( 54 ) = JVS( 1033 )
   W( 55 ) = JVS( 1034 )
   W( 56 ) = JVS( 1035 )
   W( 57 ) = JVS( 1036 )
   W( 58 ) = JVS( 1037 )
   W( 59 ) = JVS( 1038 )
   W( 60 ) = JVS( 1039 )
   W( 61 ) = JVS( 1040 )
   W( 62 ) = JVS( 1041 )
   W( 63 ) = JVS( 1042 )
   W( 64 ) = JVS( 1043 )
   W( 65 ) = JVS( 1044 )
   W( 66 ) = JVS( 1045 )
   W( 67 ) = JVS( 1046 )
   W( 68 ) = JVS( 1047 )
   W( 69 ) = JVS( 1048 )
   W( 70 ) = JVS( 1049 )
   W( 71 ) = JVS( 1050 )
   W( 72 ) = JVS( 1051 )
   W( 73 ) = JVS( 1052 )
  a = -W( 10 ) / JVS(           51  )
  W( 10 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 52 )
  W( 73 ) = W( 73 ) + a*JVS( 53 )
  a = -W( 16 ) / JVS(           71  )
  W( 16 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 72 )
  W( 68 ) = W( 68 ) + a*JVS( 73 )
  W( 73 ) = W( 73 ) + a*JVS( 74 )
  a = -W( 18 ) / JVS(           79  )
  W( 18 ) = -a
  W( 21 ) = W( 21 ) + a*JVS( 80 )
  W( 22 ) = W( 22 ) + a*JVS( 81 )
  W( 25 ) = W( 25 ) + a*JVS( 82 )
  W( 68 ) = W( 68 ) + a*JVS( 83 )
  W( 72 ) = W( 72 ) + a*JVS( 84 )
  W( 73 ) = W( 73 ) + a*JVS( 85 )
  a = -W( 20 ) / JVS(           96  )
  W( 20 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 97 )
  W( 67 ) = W( 67 ) + a*JVS( 98 )
  W( 68 ) = W( 68 ) + a*JVS( 99 )
  W( 70 ) = W( 70 ) + a*JVS( 100 )
  W( 73 ) = W( 73 ) + a*JVS( 101 )
  a = -W( 21 ) / JVS(          103  )
  W( 21 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 104 )
  W( 68 ) = W( 68 ) + a*JVS( 105 )
  W( 73 ) = W( 73 ) + a*JVS( 106 )
  a = -W( 22 ) / JVS(          108  )
  W( 22 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 109 )
  W( 68 ) = W( 68 ) + a*JVS( 110 )
  W( 73 ) = W( 73 ) + a*JVS( 111 )
  a = -W( 25 ) / JVS(          123  )
  W( 25 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 124 )
  W( 63 ) = W( 63 ) + a*JVS( 125 )
  W( 68 ) = W( 68 ) + a*JVS( 126 )
  W( 73 ) = W( 73 ) + a*JVS( 127 )
  a = -W( 27 ) / JVS(          134  )
  W( 27 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 135 )
  W( 34 ) = W( 34 ) + a*JVS( 136 )
  W( 40 ) = W( 40 ) + a*JVS( 137 )
  W( 48 ) = W( 48 ) + a*JVS( 138 )
  W( 51 ) = W( 51 ) + a*JVS( 139 )
  W( 53 ) = W( 53 ) + a*JVS( 140 )
  W( 54 ) = W( 54 ) + a*JVS( 141 )
  W( 67 ) = W( 67 ) + a*JVS( 142 )
  W( 68 ) = W( 68 ) + a*JVS( 143 )
  W( 70 ) = W( 70 ) + a*JVS( 144 )
  W( 73 ) = W( 73 ) + a*JVS( 145 )
  a = -W( 31 ) / JVS(          174  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 175 )
  W( 63 ) = W( 63 ) + a*JVS( 176 )
  W( 66 ) = W( 66 ) + a*JVS( 177 )
  W( 68 ) = W( 68 ) + a*JVS( 178 )
  W( 70 ) = W( 70 ) + a*JVS( 179 )
  W( 73 ) = W( 73 ) + a*JVS( 180 )
  a = -W( 32 ) / JVS(          185  )
  W( 32 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 186 )
  W( 67 ) = W( 67 ) + a*JVS( 187 )
  W( 68 ) = W( 68 ) + a*JVS( 188 )
  W( 70 ) = W( 70 ) + a*JVS( 189 )
  W( 73 ) = W( 73 ) + a*JVS( 190 )
  a = -W( 34 ) / JVS(          196  )
  W( 34 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 197 )
  W( 46 ) = W( 46 ) + a*JVS( 198 )
  W( 47 ) = W( 47 ) + a*JVS( 199 )
  W( 52 ) = W( 52 ) + a*JVS( 200 )
  W( 54 ) = W( 54 ) + a*JVS( 201 )
  W( 63 ) = W( 63 ) + a*JVS( 202 )
  W( 65 ) = W( 65 ) + a*JVS( 203 )
  W( 66 ) = W( 66 ) + a*JVS( 204 )
  W( 68 ) = W( 68 ) + a*JVS( 205 )
  W( 70 ) = W( 70 ) + a*JVS( 206 )
  W( 72 ) = W( 72 ) + a*JVS( 207 )
  a = -W( 35 ) / JVS(          208  )
  W( 35 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 209 )
  W( 63 ) = W( 63 ) + a*JVS( 210 )
  W( 68 ) = W( 68 ) + a*JVS( 211 )
  W( 70 ) = W( 70 ) + a*JVS( 212 )
  W( 73 ) = W( 73 ) + a*JVS( 213 )
  a = -W( 36 ) / JVS(          215  )
  W( 36 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 216 )
  W( 65 ) = W( 65 ) + a*JVS( 217 )
  W( 66 ) = W( 66 ) + a*JVS( 218 )
  W( 67 ) = W( 67 ) + a*JVS( 219 )
  W( 68 ) = W( 68 ) + a*JVS( 220 )
  W( 70 ) = W( 70 ) + a*JVS( 221 )
  W( 72 ) = W( 72 ) + a*JVS( 222 )
  a = -W( 38 ) / JVS(          228  )
  W( 38 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 229 )
  W( 65 ) = W( 65 ) + a*JVS( 230 )
  W( 66 ) = W( 66 ) + a*JVS( 231 )
  W( 67 ) = W( 67 ) + a*JVS( 232 )
  W( 68 ) = W( 68 ) + a*JVS( 233 )
  W( 70 ) = W( 70 ) + a*JVS( 234 )
  W( 72 ) = W( 72 ) + a*JVS( 235 )
  a = -W( 39 ) / JVS(          236  )
  W( 39 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 237 )
  W( 65 ) = W( 65 ) + a*JVS( 238 )
  W( 66 ) = W( 66 ) + a*JVS( 239 )
  W( 67 ) = W( 67 ) + a*JVS( 240 )
  W( 68 ) = W( 68 ) + a*JVS( 241 )
  W( 70 ) = W( 70 ) + a*JVS( 242 )
  W( 72 ) = W( 72 ) + a*JVS( 243 )
  a = -W( 40 ) / JVS(          246  )
  W( 40 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 247 )
  W( 46 ) = W( 46 ) + a*JVS( 248 )
  W( 47 ) = W( 47 ) + a*JVS( 249 )
  W( 51 ) = W( 51 ) + a*JVS( 250 )
  W( 52 ) = W( 52 ) + a*JVS( 251 )
  W( 54 ) = W( 54 ) + a*JVS( 252 )
  W( 55 ) = W( 55 ) + a*JVS( 253 )
  W( 63 ) = W( 63 ) + a*JVS( 254 )
  W( 64 ) = W( 64 ) + a*JVS( 255 )
  W( 65 ) = W( 65 ) + a*JVS( 256 )
  W( 66 ) = W( 66 ) + a*JVS( 257 )
  W( 68 ) = W( 68 ) + a*JVS( 258 )
  W( 70 ) = W( 70 ) + a*JVS( 259 )
  W( 72 ) = W( 72 ) + a*JVS( 260 )
  W( 73 ) = W( 73 ) + a*JVS( 261 )
  a = -W( 41 ) / JVS(          264  )
  W( 41 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 265 )
  W( 65 ) = W( 65 ) + a*JVS( 266 )
  W( 66 ) = W( 66 ) + a*JVS( 267 )
  W( 67 ) = W( 67 ) + a*JVS( 268 )
  W( 68 ) = W( 68 ) + a*JVS( 269 )
  W( 70 ) = W( 70 ) + a*JVS( 270 )
  W( 72 ) = W( 72 ) + a*JVS( 271 )
  W( 73 ) = W( 73 ) + a*JVS( 272 )
  a = -W( 43 ) / JVS(          299  )
  W( 43 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 300 )
  W( 65 ) = W( 65 ) + a*JVS( 301 )
  W( 66 ) = W( 66 ) + a*JVS( 302 )
  W( 67 ) = W( 67 ) + a*JVS( 303 )
  W( 68 ) = W( 68 ) + a*JVS( 304 )
  W( 70 ) = W( 70 ) + a*JVS( 305 )
  W( 72 ) = W( 72 ) + a*JVS( 306 )
  a = -W( 44 ) / JVS(          308  )
  W( 44 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 309 )
  W( 66 ) = W( 66 ) + a*JVS( 310 )
  W( 67 ) = W( 67 ) + a*JVS( 311 )
  W( 68 ) = W( 68 ) + a*JVS( 312 )
  W( 70 ) = W( 70 ) + a*JVS( 313 )
  W( 72 ) = W( 72 ) + a*JVS( 314 )
  a = -W( 45 ) / JVS(          316  )
  W( 45 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 317 )
  W( 66 ) = W( 66 ) + a*JVS( 318 )
  W( 67 ) = W( 67 ) + a*JVS( 319 )
  W( 68 ) = W( 68 ) + a*JVS( 320 )
  W( 70 ) = W( 70 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 325 )
  W( 65 ) = W( 65 ) + a*JVS( 326 )
  W( 66 ) = W( 66 ) + a*JVS( 327 )
  W( 67 ) = W( 67 ) + a*JVS( 328 )
  W( 68 ) = W( 68 ) + a*JVS( 329 )
  W( 70 ) = W( 70 ) + a*JVS( 330 )
  W( 72 ) = W( 72 ) + a*JVS( 331 )
  W( 73 ) = W( 73 ) + a*JVS( 332 )
  a = -W( 47 ) / JVS(          334  )
  W( 47 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 335 )
  W( 65 ) = W( 65 ) + a*JVS( 336 )
  W( 66 ) = W( 66 ) + a*JVS( 337 )
  W( 67 ) = W( 67 ) + a*JVS( 338 )
  W( 68 ) = W( 68 ) + a*JVS( 339 )
  W( 70 ) = W( 70 ) + a*JVS( 340 )
  W( 72 ) = W( 72 ) + a*JVS( 341 )
  W( 73 ) = W( 73 ) + a*JVS( 342 )
  a = -W( 48 ) / JVS(          363  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 364 )
  W( 50 ) = W( 50 ) + a*JVS( 365 )
  W( 51 ) = W( 51 ) + a*JVS( 366 )
  W( 52 ) = W( 52 ) + a*JVS( 367 )
  W( 54 ) = W( 54 ) + a*JVS( 368 )
  W( 55 ) = W( 55 ) + a*JVS( 369 )
  W( 56 ) = W( 56 ) + a*JVS( 370 )
  W( 57 ) = W( 57 ) + a*JVS( 371 )
  W( 58 ) = W( 58 ) + a*JVS( 372 )
  W( 59 ) = W( 59 ) + a*JVS( 373 )
  W( 60 ) = W( 60 ) + a*JVS( 374 )
  W( 61 ) = W( 61 ) + a*JVS( 375 )
  W( 62 ) = W( 62 ) + a*JVS( 376 )
  W( 63 ) = W( 63 ) + a*JVS( 377 )
  W( 64 ) = W( 64 ) + a*JVS( 378 )
  W( 65 ) = W( 65 ) + a*JVS( 379 )
  W( 66 ) = W( 66 ) + a*JVS( 380 )
  W( 67 ) = W( 67 ) + a*JVS( 381 )
  W( 68 ) = W( 68 ) + a*JVS( 382 )
  W( 70 ) = W( 70 ) + a*JVS( 383 )
  W( 72 ) = W( 72 ) + a*JVS( 384 )
  W( 73 ) = W( 73 ) + a*JVS( 385 )
  a = -W( 49 ) / JVS(          387  )
  W( 49 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 388 )
  W( 65 ) = W( 65 ) + a*JVS( 389 )
  W( 66 ) = W( 66 ) + a*JVS( 390 )
  W( 67 ) = W( 67 ) + a*JVS( 391 )
  W( 68 ) = W( 68 ) + a*JVS( 392 )
  W( 70 ) = W( 70 ) + a*JVS( 393 )
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  a = -W( 50 ) / JVS(          397  )
  W( 50 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 398 )
  W( 63 ) = W( 63 ) + a*JVS( 399 )
  W( 65 ) = W( 65 ) + a*JVS( 400 )
  W( 66 ) = W( 66 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 70 ) = W( 70 ) + a*JVS( 404 )
  W( 72 ) = W( 72 ) + a*JVS( 405 )
  a = -W( 51 ) / JVS(          411  )
  W( 51 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 412 )
  W( 63 ) = W( 63 ) + a*JVS( 413 )
  W( 65 ) = W( 65 ) + a*JVS( 414 )
  W( 66 ) = W( 66 ) + a*JVS( 415 )
  W( 67 ) = W( 67 ) + a*JVS( 416 )
  W( 68 ) = W( 68 ) + a*JVS( 417 )
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 72 ) = W( 72 ) + a*JVS( 419 )
  a = -W( 52 ) / JVS(          421  )
  W( 52 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 422 )
  W( 66 ) = W( 66 ) + a*JVS( 423 )
  W( 67 ) = W( 67 ) + a*JVS( 424 )
  W( 68 ) = W( 68 ) + a*JVS( 425 )
  W( 69 ) = W( 69 ) + a*JVS( 426 )
  W( 70 ) = W( 70 ) + a*JVS( 427 )
  W( 71 ) = W( 71 ) + a*JVS( 428 )
  W( 72 ) = W( 72 ) + a*JVS( 429 )
  a = -W( 53 ) / JVS(          442  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 443 )
  W( 57 ) = W( 57 ) + a*JVS( 444 )
  W( 58 ) = W( 58 ) + a*JVS( 445 )
  W( 59 ) = W( 59 ) + a*JVS( 446 )
  W( 60 ) = W( 60 ) + a*JVS( 447 )
  W( 61 ) = W( 61 ) + a*JVS( 448 )
  W( 62 ) = W( 62 ) + a*JVS( 449 )
  W( 63 ) = W( 63 ) + a*JVS( 450 )
  W( 64 ) = W( 64 ) + a*JVS( 451 )
  W( 65 ) = W( 65 ) + a*JVS( 452 )
  W( 66 ) = W( 66 ) + a*JVS( 453 )
  W( 67 ) = W( 67 ) + a*JVS( 454 )
  W( 68 ) = W( 68 ) + a*JVS( 455 )
  W( 69 ) = W( 69 ) + a*JVS( 456 )
  W( 70 ) = W( 70 ) + a*JVS( 457 )
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 72 ) = W( 72 ) + a*JVS( 459 )
  a = -W( 54 ) / JVS(          463  )
  W( 54 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 464 )
  W( 63 ) = W( 63 ) + a*JVS( 465 )
  W( 65 ) = W( 65 ) + a*JVS( 466 )
  W( 66 ) = W( 66 ) + a*JVS( 467 )
  W( 67 ) = W( 67 ) + a*JVS( 468 )
  W( 68 ) = W( 68 ) + a*JVS( 469 )
  W( 70 ) = W( 70 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  a = -W( 55 ) / JVS(          476  )
  W( 55 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 477 )
  W( 63 ) = W( 63 ) + a*JVS( 478 )
  W( 65 ) = W( 65 ) + a*JVS( 479 )
  W( 66 ) = W( 66 ) + a*JVS( 480 )
  W( 67 ) = W( 67 ) + a*JVS( 481 )
  W( 68 ) = W( 68 ) + a*JVS( 482 )
  W( 70 ) = W( 70 ) + a*JVS( 483 )
  W( 72 ) = W( 72 ) + a*JVS( 484 )
  W( 73 ) = W( 73 ) + a*JVS( 485 )
  a = -W( 56 ) / JVS(          500  )
  W( 56 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 501 )
  W( 63 ) = W( 63 ) + a*JVS( 502 )
  W( 64 ) = W( 64 ) + a*JVS( 503 )
  W( 65 ) = W( 65 ) + a*JVS( 504 )
  W( 66 ) = W( 66 ) + a*JVS( 505 )
  W( 67 ) = W( 67 ) + a*JVS( 506 )
  W( 68 ) = W( 68 ) + a*JVS( 507 )
  W( 69 ) = W( 69 ) + a*JVS( 508 )
  W( 70 ) = W( 70 ) + a*JVS( 509 )
  W( 71 ) = W( 71 ) + a*JVS( 510 )
  W( 72 ) = W( 72 ) + a*JVS( 511 )
  W( 73 ) = W( 73 ) + a*JVS( 512 )
  a = -W( 57 ) / JVS(          517  )
  W( 57 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 518 )
  W( 63 ) = W( 63 ) + a*JVS( 519 )
  W( 65 ) = W( 65 ) + a*JVS( 520 )
  W( 66 ) = W( 66 ) + a*JVS( 521 )
  W( 67 ) = W( 67 ) + a*JVS( 522 )
  W( 68 ) = W( 68 ) + a*JVS( 523 )
  W( 70 ) = W( 70 ) + a*JVS( 524 )
  W( 72 ) = W( 72 ) + a*JVS( 525 )
  a = -W( 58 ) / JVS(          528  )
  W( 58 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 529 )
  W( 63 ) = W( 63 ) + a*JVS( 530 )
  W( 65 ) = W( 65 ) + a*JVS( 531 )
  W( 66 ) = W( 66 ) + a*JVS( 532 )
  W( 67 ) = W( 67 ) + a*JVS( 533 )
  W( 68 ) = W( 68 ) + a*JVS( 534 )
  W( 70 ) = W( 70 ) + a*JVS( 535 )
  W( 72 ) = W( 72 ) + a*JVS( 536 )
  a = -W( 59 ) / JVS(          545  )
  W( 59 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 546 )
  W( 62 ) = W( 62 ) + a*JVS( 547 )
  W( 63 ) = W( 63 ) + a*JVS( 548 )
  W( 65 ) = W( 65 ) + a*JVS( 549 )
  W( 66 ) = W( 66 ) + a*JVS( 550 )
  W( 67 ) = W( 67 ) + a*JVS( 551 )
  W( 68 ) = W( 68 ) + a*JVS( 552 )
  W( 70 ) = W( 70 ) + a*JVS( 553 )
  W( 72 ) = W( 72 ) + a*JVS( 554 )
  a = -W( 60 ) / JVS(          563  )
  W( 60 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 564 )
  W( 63 ) = W( 63 ) + a*JVS( 565 )
  W( 65 ) = W( 65 ) + a*JVS( 566 )
  W( 66 ) = W( 66 ) + a*JVS( 567 )
  W( 67 ) = W( 67 ) + a*JVS( 568 )
  W( 68 ) = W( 68 ) + a*JVS( 569 )
  W( 70 ) = W( 70 ) + a*JVS( 570 )
  W( 72 ) = W( 72 ) + a*JVS( 571 )
  a = -W( 61 ) / JVS(          585  )
  W( 61 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 586 )
  W( 63 ) = W( 63 ) + a*JVS( 587 )
  W( 64 ) = W( 64 ) + a*JVS( 588 )
  W( 65 ) = W( 65 ) + a*JVS( 589 )
  W( 66 ) = W( 66 ) + a*JVS( 590 )
  W( 67 ) = W( 67 ) + a*JVS( 591 )
  W( 68 ) = W( 68 ) + a*JVS( 592 )
  W( 69 ) = W( 69 ) + a*JVS( 593 )
  W( 70 ) = W( 70 ) + a*JVS( 594 )
  W( 71 ) = W( 71 ) + a*JVS( 595 )
  W( 72 ) = W( 72 ) + a*JVS( 596 )
  W( 73 ) = W( 73 ) + a*JVS( 597 )
  a = -W( 62 ) / JVS(          602  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 603 )
  W( 65 ) = W( 65 ) + a*JVS( 604 )
  W( 66 ) = W( 66 ) + a*JVS( 605 )
  W( 67 ) = W( 67 ) + a*JVS( 606 )
  W( 68 ) = W( 68 ) + a*JVS( 607 )
  W( 70 ) = W( 70 ) + a*JVS( 608 )
  W( 72 ) = W( 72 ) + a*JVS( 609 )
  W( 73 ) = W( 73 ) + a*JVS( 610 )
  a = -W( 63 ) / JVS(          628  )
  W( 63 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 629 )
  W( 66 ) = W( 66 ) + a*JVS( 630 )
  W( 67 ) = W( 67 ) + a*JVS( 631 )
  W( 68 ) = W( 68 ) + a*JVS( 632 )
  W( 70 ) = W( 70 ) + a*JVS( 633 )
  W( 72 ) = W( 72 ) + a*JVS( 634 )
  W( 73 ) = W( 73 ) + a*JVS( 635 )
  a = -W( 64 ) / JVS(          652  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 653 )
  W( 66 ) = W( 66 ) + a*JVS( 654 )
  W( 67 ) = W( 67 ) + a*JVS( 655 )
  W( 68 ) = W( 68 ) + a*JVS( 656 )
  W( 69 ) = W( 69 ) + a*JVS( 657 )
  W( 70 ) = W( 70 ) + a*JVS( 658 )
  W( 71 ) = W( 71 ) + a*JVS( 659 )
  W( 72 ) = W( 72 ) + a*JVS( 660 )
  W( 73 ) = W( 73 ) + a*JVS( 661 )
  a = -W( 65 ) / JVS(          691  )
  W( 65 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 692 )
  W( 67 ) = W( 67 ) + a*JVS( 693 )
  W( 68 ) = W( 68 ) + a*JVS( 694 )
  W( 69 ) = W( 69 ) + a*JVS( 695 )
  W( 70 ) = W( 70 ) + a*JVS( 696 )
  W( 71 ) = W( 71 ) + a*JVS( 697 )
  W( 72 ) = W( 72 ) + a*JVS( 698 )
  W( 73 ) = W( 73 ) + a*JVS( 699 )
  a = -W( 66 ) / JVS(          734  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 735 )
  W( 68 ) = W( 68 ) + a*JVS( 736 )
  W( 69 ) = W( 69 ) + a*JVS( 737 )
  W( 70 ) = W( 70 ) + a*JVS( 738 )
  W( 71 ) = W( 71 ) + a*JVS( 739 )
  W( 72 ) = W( 72 ) + a*JVS( 740 )
  W( 73 ) = W( 73 ) + a*JVS( 741 )
  a = -W( 67 ) / JVS(          796  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 797 )
  W( 69 ) = W( 69 ) + a*JVS( 798 )
  W( 70 ) = W( 70 ) + a*JVS( 799 )
  W( 71 ) = W( 71 ) + a*JVS( 800 )
  W( 72 ) = W( 72 ) + a*JVS( 801 )
  W( 73 ) = W( 73 ) + a*JVS( 802 )
  a = -W( 68 ) / JVS(          860  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 861 )
  W( 70 ) = W( 70 ) + a*JVS( 862 )
  W( 71 ) = W( 71 ) + a*JVS( 863 )
  W( 72 ) = W( 72 ) + a*JVS( 864 )
  W( 73 ) = W( 73 ) + a*JVS( 865 )
  a = -W( 69 ) / JVS(          886  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 887 )
  W( 71 ) = W( 71 ) + a*JVS( 888 )
  W( 72 ) = W( 72 ) + a*JVS( 889 )
  W( 73 ) = W( 73 ) + a*JVS( 890 )
  a = -W( 70 ) / JVS(          935  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 936 )
  W( 72 ) = W( 72 ) + a*JVS( 937 )
  W( 73 ) = W( 73 ) + a*JVS( 938 )
  a = -W( 71 ) / JVS(          966  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 967 )
  W( 73 ) = W( 73 ) + a*JVS( 968 )
  a = -W( 72 ) / JVS(         1003  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 1004 )
  JVS( 1005) = W( 10 )
  JVS( 1006) = W( 16 )
  JVS( 1007) = W( 18 )
  JVS( 1008) = W( 20 )
  JVS( 1009) = W( 21 )
  JVS( 1010) = W( 22 )
  JVS( 1011) = W( 25 )
  JVS( 1012) = W( 27 )
  JVS( 1013) = W( 31 )
  JVS( 1014) = W( 32 )
  JVS( 1015) = W( 34 )
  JVS( 1016) = W( 35 )
  JVS( 1017) = W( 36 )
  JVS( 1018) = W( 38 )
  JVS( 1019) = W( 39 )
  JVS( 1020) = W( 40 )
  JVS( 1021) = W( 41 )
  JVS( 1022) = W( 43 )
  JVS( 1023) = W( 44 )
  JVS( 1024) = W( 45 )
  JVS( 1025) = W( 46 )
  JVS( 1026) = W( 47 )
  JVS( 1027) = W( 48 )
  JVS( 1028) = W( 49 )
  JVS( 1029) = W( 50 )
  JVS( 1030) = W( 51 )
  JVS( 1031) = W( 52 )
  JVS( 1032) = W( 53 )
  JVS( 1033) = W( 54 )
  JVS( 1034) = W( 55 )
  JVS( 1035) = W( 56 )
  JVS( 1036) = W( 57 )
  JVS( 1037) = W( 58 )
  JVS( 1038) = W( 59 )
  JVS( 1039) = W( 60 )
  JVS( 1040) = W( 61 )
  JVS( 1041) = W( 62 )
  JVS( 1042) = W( 63 )
  JVS( 1043) = W( 64 )
  JVS( 1044) = W( 65 )
  JVS( 1045) = W( 66 )
  JVS( 1046) = W( 67 )
  JVS( 1047) = W( 68 )
  JVS( 1048) = W( 69 )
  JVS( 1049) = W( 70 )
  JVS( 1050) = W( 71 )
  JVS( 1051) = W( 72 )
  JVS( 1052) = W( 73 )
   
   END SUBROUTINE decomp_racmpm
 


END MODULE racmpm_Integrator
