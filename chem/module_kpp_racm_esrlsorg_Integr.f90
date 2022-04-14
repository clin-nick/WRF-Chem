
































MODULE racm_esrlsorg_Integrator

 USE racm_esrlsorg_Parameters
 USE racm_esrlsorg_Precision
 USE racm_esrlsorg_JacobianSP

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

SUBROUTINE  racm_esrlsorg_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE racm_esrlsorg_Parameters

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

   CALL racm_esrlsorg_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  racm_esrlsorg_INTEGRATE


SUBROUTINE  racm_esrlsorg_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE racm_esrlsorg_Parameters

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
      CALL racm_esrlsorg_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL racm_esrlsorg_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = racm_esrlsorg_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL racm_esrlsorg_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL racm_esrlsorg_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL racm_esrlsorg_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL racm_esrlsorg_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL racm_esrlsorg_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL racm_esrlsorg_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL racm_esrlsorg_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL racm_esrlsorg_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL racm_esrlsorg_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL racm_esrlsorg_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL racm_esrlsorg_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL racm_esrlsorg_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL racm_esrlsorg_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL racm_esrlsorg_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL racm_esrlsorg_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  racm_esrlsorg_ros_ErrorMsg(Code,T,H,IERR)



   USE racm_esrlsorg_Precision

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

 END SUBROUTINE  racm_esrlsorg_ros_ErrorMsg


 SUBROUTINE  racm_esrlsorg_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL racm_esrlsorg_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL racm_esrlsorg_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL racm_esrlsorg_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL racm_esrlsorg_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL racm_esrlsorg_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL racm_esrlsorg_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL racm_esrlsorg_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL racm_esrlsorg_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL racm_esrlsorg_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL racm_esrlsorg_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL racm_esrlsorg_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL racm_esrlsorg_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL racm_esrlsorg_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL racm_esrlsorg_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL racm_esrlsorg_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL racm_esrlsorg_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL racm_esrlsorg_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL racm_esrlsorg_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL racm_esrlsorg_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL racm_esrlsorg_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = racm_esrlsorg_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL racm_esrlsorg_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  racm_esrlsorg_ros_Integrator



  REAL(kind=dp) FUNCTION  racm_esrlsorg_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    racm_esrlsorg_ros_ErrorNorm = Err

  END FUNCTION  racm_esrlsorg_ros_ErrorNorm



  SUBROUTINE racm_esrlsorg_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL racm_esrlsorg_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL racm_esrlsorg_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL racm_esrlsorg_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  racm_esrlsorg_ros_FunTimeDeriv



  SUBROUTINE  racm_esrlsorg_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL racm_esrlsorg_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL racm_esrlsorg_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL racm_esrlsorg_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  racm_esrlsorg_ros_PrepareMatrix



  SUBROUTINE  racm_esrlsorg_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 



CALL decomp_racm_esrlsorg ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  racm_esrlsorg_ros_Decomp



  SUBROUTINE  racm_esrlsorg_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL racm_esrlsorg_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  racm_esrlsorg_ros_Solve




  SUBROUTINE  racm_esrlsorg_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  racm_esrlsorg_Ros2



  SUBROUTINE  racm_esrlsorg_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racm_esrlsorg_Ros3





  SUBROUTINE  racm_esrlsorg_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racm_esrlsorg_Ros4


  SUBROUTINE  racm_esrlsorg_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racm_esrlsorg_Rodas3


  SUBROUTINE  racm_esrlsorg_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racm_esrlsorg_Rodas4




END SUBROUTINE  racm_esrlsorg_Rosenbrock




SUBROUTINE  racm_esrlsorg_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE racm_esrlsorg_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL racm_esrlsorg_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  racm_esrlsorg_FunTemplate



SUBROUTINE  racm_esrlsorg_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE racm_esrlsorg_Parameters
 
 USE racm_esrlsorg_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL racm_esrlsorg_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  racm_esrlsorg_JacTemplate

















SUBROUTINE racm_esrlsorg_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(80)
  A(2) = RCT(2)*V(77)
  A(3) = RCT(3)*V(77)
  A(4) = RCT(4)*V(23)
  A(5) = RCT(5)*V(32)
  A(6) = RCT(6)*V(19)
  A(7) = RCT(7)*V(74)
  A(8) = RCT(8)*V(74)
  A(9) = RCT(9)*V(24)
  A(10) = RCT(10)*V(56)
  A(11) = RCT(11)*V(56)
  A(12) = RCT(12)*V(64)
  A(13) = RCT(13)*V(20)
  A(14) = RCT(14)*V(70)
  A(15) = RCT(15)*V(33)
  A(16) = RCT(16)*V(49)
  A(17) = RCT(17)*V(44)
  A(18) = RCT(18)*V(44)
  A(19) = RCT(19)*V(51)
  A(20) = RCT(20)*V(61)
  A(21) = RCT(21)*V(60)
  A(22) = RCT(22)*V(48)
  A(23) = RCT(23)*V(25)
  A(24) = RCT(24)*V(26)*F(2)
  A(25) = RCT(25)*V(26)*V(77)
  A(26) = RCT(26)*V(6)*F(2)
  A(27) = RCT(27)*V(6)*F(1)
  A(28) = RCT(28)*V(73)*V(77)
  A(29) = RCT(29)*V(75)*V(77)
  A(30) = RCT(30)*V(73)*V(75)
  A(31) = RCT(31)*V(24)*V(73)
  A(32) = RCT(32)*V(75)*V(75)
  A(33) = RCT(33)*V(75)*V(75)*F(1)
  A(34) = RCT(34)*V(26)*V(79)
  A(35) = RCT(35)*V(26)*V(80)
  A(36) = RCT(36)*V(26)*V(80)
  A(37) = RCT(37)*V(73)*V(79)
  A(38) = RCT(38)*V(73)*V(80)
  A(39) = RCT(39)*V(73)*V(74)
  A(40) = RCT(40)*V(75)*V(79)
  A(41) = RCT(41)*V(75)*V(80)
  A(42) = RCT(42)*V(19)
  A(43) = RCT(43)*V(74)*V(75)
  A(44) = RCT(44)*V(23)*V(73)
  A(45) = RCT(45)*V(32)*V(73)
  A(46) = RCT(46)*V(19)*V(73)
  A(47) = RCT(47)*V(77)*V(79)
  A(48) = RCT(48)*V(77)*V(80)
  A(49) = RCT(49)*V(79)*V(79)*F(2)
  A(50) = RCT(50)*V(74)*V(79)
  A(51) = RCT(51)*V(74)*V(80)
  A(52) = RCT(52)*V(74)*V(80)
  A(53) = RCT(53)*V(11)
  A(54) = RCT(54)*V(74)*V(74)
  A(55) = RCT(55)*V(73)*F(2)
  A(56) = RCT(56)*V(5)*V(73)
  A(57) = RCT(57)*V(37)*V(73)
  A(58) = RCT(58)*V(18)*V(73)
  A(59) = RCT(59)*V(22)*V(73)
  A(60) = RCT(60)*V(15)*V(73)
  A(61) = RCT(61)*V(17)*V(73)
  A(62) = RCT(62)*V(14)*V(73)
  A(63) = RCT(63)*V(8)*V(73)
  A(64) = RCT(64)*V(12)*V(73)
  A(65) = RCT(65)*V(30)*V(73)
  A(66) = RCT(66)*V(63)*V(73)
  A(67) = RCT(67)*V(65)*V(73)
  A(68) = RCT(68)*V(41)*V(73)
  A(69) = RCT(69)*V(39)*V(73)
  A(70) = RCT(70)*V(36)*V(73)
  A(71) = RCT(71)*V(38)*V(73)
  A(72) = RCT(72)*V(9)*V(73)
  A(73) = RCT(73)*V(10)*V(73)
  A(74) = RCT(74)*V(40)*V(73)
  A(75) = RCT(75)*V(56)*V(73)
  A(76) = RCT(76)*V(64)*V(73)
  A(77) = RCT(77)*V(49)*V(73)
  A(78) = RCT(78)*V(25)*V(73)
  A(79) = RCT(79)*V(44)*V(73)
  A(80) = RCT(80)*V(51)*V(73)
  A(81) = RCT(81)*V(48)*V(73)
  A(82) = RCT(82)*V(61)*V(73)
  A(83) = RCT(83)*V(16)*V(73)
  A(84) = RCT(84)*V(20)*V(73)
  A(85) = RCT(85)*V(70)*V(73)
  A(86) = RCT(86)*V(33)*V(73)
  A(87) = RCT(87)*V(35)*V(73)
  A(88) = RCT(88)*V(43)*V(73)
  A(89) = RCT(89)*V(60)*V(73)
  A(90) = RCT(90)*V(56)*V(74)
  A(91) = RCT(91)*V(64)*V(74)
  A(92) = RCT(92)*V(44)*V(74)
  A(93) = RCT(93)*V(51)*V(74)
  A(94) = RCT(94)*V(13)*V(73)
  A(95) = RCT(95)*V(61)*V(74)
  A(96) = RCT(96)*V(40)*V(74)
  A(97) = RCT(97)*V(30)*V(74)
  A(98) = RCT(98)*V(63)*V(74)
  A(99) = RCT(99)*V(65)*V(74)
  A(100) = RCT(100)*V(41)*V(74)
  A(101) = RCT(101)*V(39)*V(74)
  A(102) = RCT(102)*V(36)*V(74)
  A(103) = RCT(103)*V(38)*V(74)
  A(104) = RCT(104)*V(43)*V(74)
  A(105) = RCT(105)*V(30)*V(77)
  A(106) = RCT(106)*V(63)*V(77)
  A(107) = RCT(107)*V(65)*V(77)
  A(108) = RCT(108)*V(41)*V(77)
  A(109) = RCT(109)*V(39)*V(77)
  A(110) = RCT(110)*V(36)*V(77)
  A(111) = RCT(111)*V(38)*V(77)
  A(112) = RCT(112)*V(48)*V(77)
  A(113) = RCT(113)*V(61)*V(77)
  A(114) = RCT(114)*V(43)*V(77)
  A(115) = RCT(115)*V(27)*V(80)
  A(116) = RCT(116)*V(27)*V(75)
  A(117) = RCT(117)*V(28)*V(80)
  A(118) = RCT(118)*V(28)*F(2)
  A(119) = RCT(119)*V(28)*V(77)
  A(120) = RCT(120)*V(29)*V(80)
  A(121) = RCT(121)*V(29)*F(2)
  A(122) = RCT(122)*V(29)*V(77)
  A(123) = RCT(123)*V(31)*V(80)
  A(124) = RCT(124)*V(31)*F(2)
  A(125) = RCT(125)*V(31)*V(77)
  A(126) = RCT(126)*V(72)*V(80)
  A(127) = RCT(127)*V(35)
  A(128) = RCT(128)*V(78)*V(80)
  A(129) = RCT(129)*V(43)
  A(130) = RCT(130)*V(76)*V(79)
  A(131) = RCT(131)*V(68)*V(79)
  A(132) = RCT(132)*V(71)*V(79)
  A(133) = RCT(133)*V(55)*V(79)
  A(134) = RCT(134)*V(53)*V(79)
  A(135) = RCT(135)*V(45)*V(79)
  A(136) = RCT(136)*V(46)*V(79)
  A(137) = RCT(137)*V(47)*V(79)
  A(138) = RCT(138)*V(50)*V(79)
  A(139) = RCT(139)*V(59)*V(79)
  A(140) = RCT(140)*V(54)*V(79)
  A(141) = RCT(141)*V(57)*V(79)
  A(142) = RCT(142)*V(58)*V(79)
  A(143) = RCT(143)*V(52)*V(79)
  A(144) = RCT(144)*V(72)*V(79)
  A(145) = RCT(145)*V(78)*V(79)
  A(146) = RCT(146)*V(69)*V(79)
  A(147) = RCT(147)*V(66)*V(79)
  A(148) = RCT(148)*V(67)*V(79)
  A(149) = RCT(149)*V(75)*V(76)
  A(150) = RCT(150)*V(68)*V(75)
  A(151) = RCT(151)*V(71)*V(75)
  A(152) = RCT(152)*V(55)*V(75)
  A(153) = RCT(153)*V(53)*V(75)
  A(154) = RCT(154)*V(45)*V(75)
  A(155) = RCT(155)*V(46)*V(75)
  A(156) = RCT(156)*V(47)*V(75)
  A(157) = RCT(157)*V(50)*V(75)
  A(158) = RCT(158)*V(59)*V(75)
  A(159) = RCT(159)*V(54)*V(75)
  A(160) = RCT(160)*V(57)*V(75)
  A(161) = RCT(161)*V(58)*V(75)
  A(162) = RCT(162)*V(52)*V(75)
  A(163) = RCT(163)*V(72)*V(75)
  A(164) = RCT(164)*V(72)*V(75)
  A(165) = RCT(165)*V(75)*V(78)
  A(166) = RCT(166)*V(75)*V(78)
  A(167) = RCT(167)*V(69)*V(75)
  A(168) = RCT(168)*V(66)*V(75)
  A(169) = RCT(169)*V(67)*V(75)
  A(170) = RCT(170)*V(76)*V(76)
  A(171) = RCT(171)*V(68)*V(76)
  A(172) = RCT(172)*V(71)*V(76)
  A(173) = RCT(173)*V(55)*V(76)
  A(174) = RCT(174)*V(53)*V(76)
  A(175) = RCT(175)*V(45)*V(76)
  A(176) = RCT(176)*V(46)*V(76)
  A(177) = RCT(177)*V(47)*V(76)
  A(178) = RCT(178)*V(50)*V(76)
  A(179) = RCT(179)*V(59)*V(76)
  A(180) = RCT(180)*V(54)*V(76)
  A(181) = RCT(181)*V(57)*V(76)
  A(182) = RCT(182)*V(58)*V(76)
  A(183) = RCT(183)*V(52)*V(76)
  A(184) = RCT(184)*V(72)*V(76)
  A(185) = RCT(185)*V(72)*V(76)
  A(186) = RCT(186)*V(76)*V(78)
  A(187) = RCT(187)*V(76)*V(78)
  A(188) = RCT(188)*V(69)*V(76)
  A(189) = RCT(189)*V(66)*V(76)
  A(190) = RCT(190)*V(67)*V(76)
  A(191) = RCT(191)*V(68)*V(72)
  A(192) = RCT(192)*V(71)*V(72)
  A(193) = RCT(193)*V(55)*V(72)
  A(194) = RCT(194)*V(53)*V(72)
  A(195) = RCT(195)*V(45)*V(72)
  A(196) = RCT(196)*V(46)*V(72)
  A(197) = RCT(197)*V(47)*V(72)
  A(198) = RCT(198)*V(50)*V(72)
  A(199) = RCT(199)*V(59)*V(72)
  A(200) = RCT(200)*V(54)*V(72)
  A(201) = RCT(201)*V(57)*V(72)
  A(202) = RCT(202)*V(58)*V(72)
  A(203) = RCT(203)*V(52)*V(72)
  A(204) = RCT(204)*V(72)*V(72)
  A(205) = RCT(205)*V(72)*V(78)
  A(206) = RCT(206)*V(69)*V(72)
  A(207) = RCT(207)*V(66)*V(72)
  A(208) = RCT(208)*V(67)*V(72)
  A(209) = RCT(209)*V(66)*V(66)
  A(210) = RCT(210)*V(66)*V(67)
  A(211) = RCT(211)*V(67)*V(67)
  A(212) = RCT(212)*V(74)*V(76)
  A(213) = RCT(213)*V(68)*V(74)
  A(214) = RCT(214)*V(71)*V(74)
  A(215) = RCT(215)*V(55)*V(74)
  A(216) = RCT(216)*V(53)*V(74)
  A(217) = RCT(217)*V(45)*V(74)
  A(218) = RCT(218)*V(46)*V(74)
  A(219) = RCT(219)*V(47)*V(74)
  A(220) = RCT(220)*V(21)*V(73)
  A(221) = RCT(221)*V(59)*V(74)
  A(222) = RCT(222)*V(54)*V(74)
  A(223) = RCT(223)*V(57)*V(74)
  A(224) = RCT(224)*V(58)*V(74)
  A(225) = RCT(225)*V(52)*V(74)
  A(226) = RCT(226)*V(72)*V(74)
  A(227) = RCT(227)*V(74)*V(78)
  A(228) = RCT(228)*V(69)*V(74)
  A(229) = RCT(229)*V(66)*V(74)
  A(230) = RCT(230)*V(67)*V(74)
  A(231) = RCT(231)*V(62)*V(75)
  A(232) = RCT(232)*V(62)*V(76)
  A(233) = RCT(233)*V(62)*V(72)
  A(234) = RCT(234)*V(62)*V(62)
  A(235) = RCT(235)*V(62)*V(79)
  A(236) = RCT(236)*V(62)*V(74)
  A(237) = RCT(237)*V(50)*V(50)
  A(238) = RCT(238)*V(7)*V(73)
  A(239) = RCT(239)*V(34)*V(73)
  A(240) = RCT(240)*V(42)*V(79)
  A(241) = RCT(241)*V(42)*V(75)
  A(242) = RCT(242)*V(42)*V(42)
  A(243) = RCT(243)*V(42)*V(80)
  A(244) = RCT(244)*V(21)


  Vdot(1) = A(56)
  Vdot(2) = A(57)
  Vdot(3) = 0.036*A(62)+0.37*A(105)+0.14*A(106)+0.15*A(108)+0.28*A(109)+0.01*A(111)+0.45*A(112)+0.11*A(113)+0.11*A(114)
  Vdot(4) = 0.1*A(106)+0.14*A(107)+0.07*A(111)+0.21*A(113)+A(164)+A(166)+A(185)+A(187)+0.5*A(191)+0.499*A(192)+0.495&
              &*A(193)+0.495*A(194)+0.5*A(195)+0.499*A(196)+0.49*A(197)+0.494*A(198)+0.5*A(206)+0.5*A(207)+0.484*A(208)
  Vdot(5) = -A(56)
  Vdot(6) = A(2)-A(26)-A(27)
  Vdot(7) = A(157)-A(238)
  Vdot(8) = -A(63)
  Vdot(9) = -A(72)
  Vdot(10) = -A(73)
  Vdot(11) = A(52)-A(53)
  Vdot(12) = -A(64)
  Vdot(13) = -A(94)+A(241)
  Vdot(14) = -A(62)
  Vdot(15) = -A(60)+0.06*A(106)+0.07*A(107)
  Vdot(16) = 0.35*A(82)-A(83)
  Vdot(17) = -A(61)+0.03*A(106)+0.06*A(107)
  Vdot(18) = -A(58)+A(239)
  Vdot(19) = -A(6)+A(41)-A(42)-A(46)
  Vdot(20) = -A(13)-A(84)+A(149)
  Vdot(21) = -A(220)+A(243)-A(244)
  Vdot(22) = -A(59)+A(220)+A(239)+0.25*A(240)+A(242)
  Vdot(23) = -A(4)+A(37)-A(44)+A(117)+A(120)+A(123)
  Vdot(24) = -A(9)-A(31)+A(32)+A(33)+0.006*A(106)+0.011*A(107)+0.001*A(108)+0.09*A(109)+0.02*A(110)+0.02*A(111)
  Vdot(25) = -A(23)+0.024*A(64)-A(78)+0.6*A(88)+0.3*A(188)
  Vdot(26) = A(1)+A(3)+A(8)-A(24)-A(25)+A(26)-A(34)-A(35)-A(36)+0.09*A(108)
  Vdot(27) = 0.1*A(74)+A(96)-A(115)-A(116)
  Vdot(28) = 0.9*A(72)-A(117)-A(118)-A(119)
  Vdot(29) = 0.9*A(73)-A(120)-A(121)-A(122)
  Vdot(30) = -A(65)-A(97)-A(105)
  Vdot(31) = 0.85*A(74)-A(123)-A(124)-A(125)
  Vdot(32) = -A(5)+A(38)+0.3*A(43)-A(45)+A(90)+A(91)+A(92)+A(93)+0.5*A(95)+A(96)
  Vdot(33) = -A(15)-A(86)+0.11*A(113)+A(163)
  Vdot(34) = A(101)+0.046*A(138)-A(239)
  Vdot(35) = -A(87)+0.4*A(88)+0.4*A(104)+0.3*A(114)+A(126)-A(127)
  Vdot(36) = -A(70)-A(102)-A(110)
  Vdot(37) = A(10)+A(11)+A(12)+1.87*A(17)+1.55*A(18)+A(19)+A(22)-A(57)+A(58)+0.036*A(62)+A(75)+2*A(79)+A(80)+A(90)+2&
               &*A(92)+A(93)+0.43*A(105)+0.37*A(106)+0.3*A(107)+0.36*A(108)+0.14*A(109)+0.14*A(110)+0.14*A(111)+0.22*A(112)&
               &+0.66*A(113)+0.13*A(114)+0.25*A(240)+0.5*A(242)
  Vdot(38) = -A(71)-A(103)-A(111)
  Vdot(39) = -A(69)-A(101)-A(109)
  Vdot(40) = -A(74)-A(96)+0.1*A(115)+A(116)+A(117)+0.02*A(118)+A(119)+A(120)+0.02*A(121)+A(122)+A(123)+0.02*A(124)&
               &+A(125)
  Vdot(41) = -A(68)-A(100)-A(108)
  Vdot(42) = A(81)+A(94)+0.1*A(109)-A(240)-A(241)-2*A(242)-A(243)+A(244)
  Vdot(43) = -A(88)-A(104)-A(114)+A(128)-A(129)
  Vdot(44) = -A(17)-A(18)+0.036*A(62)-A(79)+0.15*A(82)-A(92)+0.25*A(95)+0.5*A(113)+0.063*A(132)+1.2*A(141)+0.35*A(142)&
               &+A(143)+0.119*A(172)+0.65*A(181)+0.37*A(182)+A(183)+0.1*A(192)+0.65*A(201)+0.37*A(202)+A(203)+0.063*A(214)&
               &+1.3*A(223)+0.74*A(224)+A(225)
  Vdot(45) = A(65)-A(135)-A(154)-A(175)-A(195)-A(217)
  Vdot(46) = A(66)-A(136)-A(155)-A(176)-A(196)-A(218)
  Vdot(47) = A(67)-A(137)-A(156)-A(177)-A(197)-A(219)
  Vdot(48) = -A(22)-A(81)+0.9*A(100)+0.39*A(108)+0.65*A(109)+0.79*A(111)-A(112)+A(138)+0.4*A(140)+0.55*A(178)+0.6*A(180)&
               &+0.771*A(198)+0.6*A(200)+0.6*A(222)+2*A(237)+A(238)
  Vdot(49) = -A(16)+0.8*A(21)+0.25*A(63)-A(77)+0.12*A(83)+0.41*A(85)+0.03*A(95)+0.03*A(106)+0.16*A(107)+0.53*A(110)&
               &+0.623*A(132)+0.722*A(133)+0.642*A(134)+0.06*A(136)+0.29*A(137)+0.8*A(139)+0.464*A(148)+0.018*A(172)+0.24&
               &*A(173)+0.419*A(174)+0.081*A(176)+0.313*A(177)+A(179)+0.149*A(190)+0.127*A(192)+0.33*A(193)+0.581*A(194)&
               &+0.141*A(196)+0.569*A(197)+A(199)+0.11*A(206)+0.167*A(208)+0.149*A(210)+0.285*A(211)+0.67*A(214)+0.828&
               &*A(215)+0.88*A(216)+0.06*A(218)+0.29*A(219)+A(221)+0.469*A(230)
  Vdot(50) = A(68)+A(69)-A(138)-A(157)-A(178)-A(198)-2*A(237)
  Vdot(51) = -A(19)+A(59)+A(78)-A(80)+0.15*A(82)-A(93)+0.25*A(95)+0.9*A(112)+0.62*A(113)+0.65*A(141)+0.6*A(142)+A(143)&
               &+0.54*A(146)+0.005*A(172)+0.35*A(181)+0.63*A(182)+A(183)+0.4*A(188)+0.004*A(192)+0.35*A(201)+0.63*A(202)&
               &+A(203)+0.54*A(206)+0.7*A(223)+1.26*A(224)+A(225)+0.54*A(228)+0.5*A(240)+A(242)
  Vdot(52) = 0.98*A(124)-A(143)-A(162)-A(183)-A(203)-A(225)
  Vdot(53) = 0.9511*A(64)-A(134)-A(153)-A(174)-A(194)-A(216)
  Vdot(54) = A(71)-A(140)-A(159)-A(180)-A(200)-A(222)
  Vdot(55) = 0.75*A(63)-A(133)-A(152)-A(173)-A(193)-A(215)
  Vdot(56) = -A(10)-A(11)+A(13)+0.13*A(17)+0.45*A(18)+A(22)+A(23)+A(58)+0.01*A(62)-A(75)+0.35*A(84)+0.35*A(86)+A(87)+0.4&
               &*A(88)-A(90)+0.4*A(104)+A(105)+0.64*A(106)+0.02*A(107)+0.9*A(108)+0.58*A(109)+0.04*A(111)+0.7*A(114)+A(130)&
               &+0.047*A(132)+0.021*A(133)+1.6*A(135)+A(136)+A(138)+0.25*A(140)+A(145)+0.287*A(148)+1.33*A(170)+0.75*A(171)&
               &+0.81*A(172)+0.829*A(173)+0.753*A(174)+1.55*A(175)+1.25*A(176)+0.755*A(177)+1.09*A(178)+A(179)+1.4*A(180)&
               &+A(181)+A(182)+A(183)+A(184)+A(185)+2*A(186)+A(187)+0.75*A(188)+0.75*A(189)+0.96*A(190)+0.091*A(192)+0.076&
               &*A(193)+0.8*A(195)+0.501*A(196)+0.34*A(198)+0.4*A(200)+A(205)+0.207*A(208)+0.202*A(210)+0.504*A(211)+A(212)&
               &+0.048*A(214)+0.021*A(215)+1.6*A(217)+A(218)+0.4*A(222)+A(227)+0.28*A(230)+A(232)+A(237)+0.75*A(240)+0.5&
               &*A(242)
  Vdot(57) = 0.98*A(118)-A(141)-A(160)-A(181)-A(201)-A(223)
  Vdot(58) = 0.98*A(121)-A(142)-A(161)-A(182)-A(202)-A(224)
  Vdot(59) = A(70)-A(139)-A(158)-A(179)-A(199)-A(221)
  Vdot(60) = -A(21)-A(89)+0.6*A(104)+A(115)+0.059*A(132)+0.124*A(133)+0.261*A(134)+0.2*A(139)+0.35*A(140)+0.05*A(141)&
               &+0.05*A(142)+A(147)+A(168)+A(169)+A(189)+0.5*A(190)+A(207)+0.484*A(208)+2*A(209)+1.5*A(210)+A(211)+A(229)
  Vdot(61) = -A(20)-A(82)-A(95)-A(113)+0.5*A(141)+0.95*A(142)+A(181)+A(182)+A(201)+A(202)+0.5*A(223)+A(224)
  Vdot(62) = 0.1*A(72)+0.1*A(73)+0.05*A(74)+0.5*A(82)+0.07*A(85)+0.35*A(86)+A(87)+A(88)+0.5*A(95)+A(104)+0.13*A(108)&
               &+0.048*A(132)+0.334*A(133)+0.416*A(134)+0.16*A(146)+0.085*A(172)+0.245*A(173)+0.322*A(174)+0.08*A(188)+0.071&
               &*A(192)+0.237*A(193)+0.318*A(194)+0.08*A(206)+0.051*A(214)+0.391*A(215)+0.587*A(216)+0.16*A(228)-A(231)&
               &-A(232)-A(233)-2*A(234)-A(235)-A(236)
  Vdot(63) = -A(66)-A(98)-A(106)+0.35*A(108)+0.46*A(111)+0.37*A(178)+0.229*A(198)
  Vdot(64) = -A(12)+A(14)+0.2*A(21)+0.335*A(62)+0.025*A(64)-A(76)+0.88*A(83)+0.08*A(85)-A(91)+0.25*A(95)+0.44*A(106)&
               &+0.99*A(107)+0.65*A(110)+0.16*A(113)+A(131)+0.233*A(132)+0.211*A(133)+0.15*A(134)+0.2*A(135)+0.94*A(136)&
               &+1.71*A(137)+0.8*A(139)+0.46*A(146)+1.24*A(148)+0.75*A(171)+0.58*A(172)+0.523*A(173)+0.411*A(174)+0.35&
               &*A(175)+0.669*A(176)+0.932*A(177)+A(179)+0.3*A(188)+0.64*A(190)+A(191)+0.724*A(192)+0.677*A(193)+0.497&
               &*A(194)+0.6*A(195)+0.859*A(196)+0.941*A(197)+A(199)+0.35*A(206)+0.65*A(208)+0.64*A(210)+1.21*A(211)+A(213)&
               &+0.243*A(214)+0.239*A(215)+0.187*A(216)+0.2*A(217)+0.94*A(218)+1.71*A(219)+A(221)+0.46*A(228)+1.24*A(230)
  Vdot(65) = -A(67)-A(99)-A(107)+0.25*A(140)+0.08*A(178)+0.4*A(180)+0.4*A(200)+0.4*A(222)
  Vdot(66) = 0.8*A(97)+0.43*A(98)+0.11*A(99)+0.9*A(100)+0.1*A(102)+0.13*A(103)-A(147)-A(168)-A(189)-A(207)-2*A(209)&
               &-A(210)-A(229)
  Vdot(67) = 0.2*A(97)+0.57*A(98)+0.89*A(99)+0.1*A(100)+0.9*A(102)+0.87*A(103)-A(148)-A(169)-A(190)-A(208)-A(210)-2&
               &*A(211)-A(230)
  Vdot(68) = A(16)+A(61)+0.1*A(106)+0.18*A(107)+0.2*A(110)+0.16*A(111)-A(131)+0.048*A(132)+0.245*A(133)+0.133*A(134)&
               &-A(150)-A(171)+0.014*A(173)+0.013*A(174)-A(191)+0.006*A(192)+0.018*A(193)+0.015*A(194)-A(213)+0.053*A(214)&
               &+0.262*A(215)+0.155*A(216)
  Vdot(69) = A(77)+0.03*A(106)+0.12*A(107)+0.02*A(108)+0.42*A(110)+0.42*A(111)-A(146)-A(167)-A(188)-A(206)-A(228)
  Vdot(70) = -A(14)-A(85)+A(150)+A(151)+A(152)+A(153)+A(154)+A(155)+A(156)+A(158)+A(159)+A(160)+A(161)+A(162)+A(165)&
               &+A(167)+A(231)
  Vdot(71) = 0.583*A(62)+0.44*A(85)+A(89)-A(132)-A(151)-A(172)-A(192)-A(214)
  Vdot(72) = A(16)+A(19)+A(22)+A(23)+A(76)+A(80)+0.65*A(86)+A(91)+A(93)+0.15*A(108)+0.1*A(109)+0.1*A(112)+0.28*A(113)&
               &+0.7*A(114)-A(126)+A(127)-A(144)+A(145)+0.23*A(146)-A(163)-A(164)-A(184)-A(185)+A(186)+0.12*A(188)-A(191)&
               &-A(192)-A(193)-A(194)-A(195)-A(196)-A(197)-A(198)-A(199)-A(200)-A(201)-A(202)-A(203)-2*A(204)-0.88*A(206)&
               &-A(207)-A(208)-A(226)+A(227)+0.23*A(228)-A(233)+0.25*A(240)
  Vdot(73) = A(4)+A(5)+0.35*A(6)+2*A(9)+A(13)+A(14)+A(15)+2*A(27)-A(28)+A(29)-A(30)-A(31)-A(37)-A(38)-A(39)+A(40)+0.7&
               &*A(43)-A(44)-A(45)-A(46)-A(55)-A(56)-A(57)-A(58)-A(59)-A(60)-A(61)-0.964*A(62)-A(63)-A(64)-A(65)-A(66)-A(67)&
               &-A(68)-A(69)-A(70)-A(71)-A(72)-A(73)-A(74)-A(75)-A(76)-A(77)-A(78)-A(79)-A(80)-A(81)-A(82)-A(83)-0.65*A(84)&
               &-0.51*A(85)-A(86)-A(87)-A(88)-A(89)-A(94)+0.12*A(105)+0.4*A(106)+0.63*A(107)+0.28*A(108)+0.25*A(109)+0.85&
               &*A(110)+0.85*A(111)+0.19*A(112)+0.21*A(113)+0.036*A(114)+A(119)+A(122)+A(125)-A(220)-A(239)
  Vdot(74) = 0.35*A(6)-A(7)-A(8)+A(36)-A(39)-A(43)+A(45)+A(48)-A(50)-A(51)-A(52)+A(53)-2*A(54)+A(87)+0.6*A(88)-A(90)&
               &-A(91)-A(92)-A(93)-A(95)-A(96)-A(97)-A(98)-A(99)-A(100)-A(101)-A(102)-A(103)-0.4*A(104)-A(212)-A(213)-A(214)&
               &-A(215)-A(216)-A(217)-A(218)-A(219)-A(221)-A(222)-A(223)-A(224)-A(225)-A(226)-A(227)-A(228)-A(229)-A(230)&
               &-A(236)
  Vdot(75) = 0.65*A(6)+2*A(11)+A(12)+A(13)+A(14)+0.8*A(18)+A(19)+A(20)+A(21)+A(22)+A(23)+A(28)-A(29)-A(30)+A(31)-2*A(32)&
               &-2*A(33)+A(39)-A(40)-A(41)+A(42)-A(43)+A(55)+A(56)+A(57)+A(59)+0.381*A(62)+0.25*A(63)+0.049*A(64)+0.1*A(72)&
               &+0.1*A(73)+0.05*A(74)+A(75)+A(78)+A(79)+0.5*A(82)+A(83)+0.35*A(86)+0.4*A(88)+A(90)+A(92)+0.5*A(95)+0.26&
               &*A(105)+0.25*A(106)+0.22*A(107)+0.3*A(108)+0.25*A(109)+0.1*A(110)+0.1*A(111)+0.32*A(112)+0.29*A(113)+0.08&
               &*A(114)-A(116)+0.02*A(118)+0.02*A(121)+0.02*A(124)+A(130)+A(131)+0.742*A(132)+0.599*A(133)+0.606*A(134)&
               &+A(135)+A(136)+A(137)+A(138)+0.8*A(139)+0.65*A(140)+0.95*A(141)+0.95*A(142)+A(143)+0.77*A(146)+A(147)-A(149)&
               &-A(150)-A(151)-A(152)-A(153)-A(154)-A(155)-A(156)-A(157)-A(158)-A(159)-A(160)-A(161)-A(162)-A(163)-A(164)&
               &-A(165)-A(166)-A(167)-A(168)-A(169)+0.66*A(170)+A(171)+0.992*A(172)+0.946*A(173)+0.993*A(174)+A(175)+A(176)&
               &+A(177)+A(178)+2*A(179)+2*A(180)+A(181)+A(182)+2*A(183)+A(184)+A(186)+0.88*A(188)+A(189)+0.5*A(190)+0.5&
               &*A(191)+0.488*A(192)+0.438*A(193)+0.489*A(194)+0.5*A(195)+0.501*A(196)+0.51*A(197)+0.506*A(198)+A(199)&
               &+A(200)+A(201)+A(202)+A(203)+0.38*A(206)+0.5*A(207)+A(209)+0.5*A(210)+A(212)+A(213)+0.792*A(214)+0.699&
               &*A(215)+0.845*A(216)+A(217)+A(218)+A(219)+A(221)+A(222)+A(223)+A(224)+A(225)+0.77*A(228)+A(229)-A(231)&
               &+A(232)+A(237)+0.75*A(240)-A(241)+A(242)
  Vdot(76) = A(12)+A(15)+A(60)+0.65*A(84)+0.19*A(106)+0.23*A(107)+0.03*A(108)+0.08*A(109)-A(130)+0.15*A(132)+0.031&
               &*A(133)+A(144)-A(149)-2*A(170)-A(171)-0.993*A(172)-0.951*A(173)-A(174)-A(175)-A(176)-A(177)-A(178)-A(179)&
               &-A(180)-A(181)-A(182)-A(183)-A(185)-A(186)-A(187)-A(188)-A(189)-A(190)+0.5*A(191)+0.508*A(192)+0.554*A(193)&
               &+0.507*A(194)+0.5*A(195)+0.501*A(196)+0.51*A(197)+0.506*A(198)+A(199)+A(200)+A(201)+A(202)+A(203)+2*A(204)&
               &+A(205)+0.5*A(206)+0.5*A(207)+0.516*A(208)-A(212)+0.155*A(214)+0.04*A(215)+A(226)-A(232)+A(233)
  Vdot(77) = -A(2)-A(3)+A(24)-A(25)-A(28)-A(29)-A(47)-A(48)-A(105)-A(106)-A(107)-A(108)-A(109)-A(110)-A(111)-A(112)&
               &-A(113)-A(114)-A(119)-A(122)-A(125)+A(164)+A(166)
  Vdot(78) = A(20)+0.5*A(82)+0.5*A(95)-A(128)+A(129)-A(145)-A(165)-A(166)-A(186)-A(187)-A(205)-A(227)
  Vdot(79) = A(1)+A(4)+A(7)-A(34)+A(35)-A(37)-A(40)-A(47)-2*A(49)-A(50)+A(51)-A(130)-A(131)-A(132)-A(133)-A(134)-A(135)&
               &-A(136)-A(137)-A(138)-A(139)-A(140)-A(141)-A(142)-A(143)-A(144)-A(145)-A(146)-A(147)-A(148)-A(235)-A(240)
  Vdot(80) = -A(1)+A(5)+0.65*A(6)+A(8)+A(21)+A(34)-A(35)-A(36)-A(38)+A(39)+A(40)-A(41)+A(42)+0.7*A(43)+A(44)+A(46)+A(47)&
               &-A(48)+2*A(49)+2*A(50)-A(52)+A(53)+2*A(54)+A(58)+A(89)+0.5*A(95)+0.4*A(104)+0.7*A(114)-A(115)-A(117)-A(120)&
               &-A(123)-A(126)+A(127)-A(128)+A(129)+A(130)+A(131)+0.941*A(132)+0.876*A(133)+0.739*A(134)+A(135)+A(136)&
               &+A(137)+A(138)+0.8*A(139)+0.65*A(140)+0.95*A(141)+0.95*A(142)+A(143)+A(144)+A(145)+A(146)+A(147)+2*A(148)&
               &+0.5*A(190)+0.516*A(208)+0.5*A(210)+A(211)+A(212)+A(213)+A(214)+A(215)+A(216)+A(217)+A(218)+A(219)+A(220)&
               &+A(221)+A(222)+A(223)+A(224)+A(225)+A(226)+A(227)+A(228)+A(229)+2*A(230)+A(235)+A(236)+A(240)-A(243)+A(244)
      
END SUBROUTINE racm_esrlsorg_Fun
















SUBROUTINE racm_esrlsorg_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(80)
  IRR(2) = RCT(2)*V(77)
  IRR(3) = RCT(3)*V(77)
  IRR(4) = RCT(4)*V(23)
  IRR(5) = RCT(5)*V(32)
  IRR(6) = RCT(6)*V(19)
  IRR(7) = RCT(7)*V(74)
  IRR(8) = RCT(8)*V(74)
  IRR(9) = RCT(9)*V(24)
  IRR(10) = RCT(10)*V(56)
  IRR(11) = RCT(11)*V(56)
  IRR(12) = RCT(12)*V(64)
  IRR(13) = RCT(13)*V(20)
  IRR(14) = RCT(14)*V(70)
  IRR(15) = RCT(15)*V(33)
  IRR(16) = RCT(16)*V(49)
  IRR(17) = RCT(17)*V(44)
  IRR(18) = RCT(18)*V(44)
  IRR(19) = RCT(19)*V(51)
  IRR(20) = RCT(20)*V(61)
  IRR(21) = RCT(21)*V(60)
  IRR(22) = RCT(22)*V(48)
  IRR(23) = RCT(23)*V(25)
  IRR(24) = RCT(24)*V(26)*F(2)
  IRR(25) = RCT(25)*V(26)*V(77)
  IRR(26) = RCT(26)*V(6)*F(2)
  IRR(27) = RCT(27)*V(6)*F(1)
  IRR(28) = RCT(28)*V(73)*V(77)
  IRR(29) = RCT(29)*V(75)*V(77)
  IRR(30) = RCT(30)*V(73)*V(75)
  IRR(31) = RCT(31)*V(24)*V(73)
  IRR(32) = RCT(32)*V(75)*V(75)
  IRR(33) = RCT(33)*V(75)*V(75)*F(1)
  IRR(34) = RCT(34)*V(26)*V(79)
  IRR(35) = RCT(35)*V(26)*V(80)
  IRR(36) = RCT(36)*V(26)*V(80)
  IRR(37) = RCT(37)*V(73)*V(79)
  IRR(38) = RCT(38)*V(73)*V(80)
  IRR(39) = RCT(39)*V(73)*V(74)
  IRR(40) = RCT(40)*V(75)*V(79)
  IRR(41) = RCT(41)*V(75)*V(80)
  IRR(42) = RCT(42)*V(19)
  IRR(43) = RCT(43)*V(74)*V(75)
  IRR(44) = RCT(44)*V(23)*V(73)
  IRR(45) = RCT(45)*V(32)*V(73)
  IRR(46) = RCT(46)*V(19)*V(73)
  IRR(47) = RCT(47)*V(77)*V(79)
  IRR(48) = RCT(48)*V(77)*V(80)
  IRR(49) = RCT(49)*V(79)*V(79)*F(2)
  IRR(50) = RCT(50)*V(74)*V(79)
  IRR(51) = RCT(51)*V(74)*V(80)
  IRR(52) = RCT(52)*V(74)*V(80)
  IRR(53) = RCT(53)*V(11)
  IRR(54) = RCT(54)*V(74)*V(74)
  IRR(55) = RCT(55)*V(73)*F(2)
  IRR(56) = RCT(56)*V(5)*V(73)
  IRR(57) = RCT(57)*V(37)*V(73)
  IRR(58) = RCT(58)*V(18)*V(73)
  IRR(59) = RCT(59)*V(22)*V(73)
  IRR(60) = RCT(60)*V(15)*V(73)
  IRR(61) = RCT(61)*V(17)*V(73)
  IRR(62) = RCT(62)*V(14)*V(73)
  IRR(63) = RCT(63)*V(8)*V(73)
  IRR(64) = RCT(64)*V(12)*V(73)
  IRR(65) = RCT(65)*V(30)*V(73)
  IRR(66) = RCT(66)*V(63)*V(73)
  IRR(67) = RCT(67)*V(65)*V(73)
  IRR(68) = RCT(68)*V(41)*V(73)
  IRR(69) = RCT(69)*V(39)*V(73)
  IRR(70) = RCT(70)*V(36)*V(73)
  IRR(71) = RCT(71)*V(38)*V(73)
  IRR(72) = RCT(72)*V(9)*V(73)
  IRR(73) = RCT(73)*V(10)*V(73)
  IRR(74) = RCT(74)*V(40)*V(73)
  IRR(75) = RCT(75)*V(56)*V(73)
  IRR(76) = RCT(76)*V(64)*V(73)
  IRR(77) = RCT(77)*V(49)*V(73)
  IRR(78) = RCT(78)*V(25)*V(73)
  IRR(79) = RCT(79)*V(44)*V(73)
  IRR(80) = RCT(80)*V(51)*V(73)
  IRR(81) = RCT(81)*V(48)*V(73)
  IRR(82) = RCT(82)*V(61)*V(73)
  IRR(83) = RCT(83)*V(16)*V(73)
  IRR(84) = RCT(84)*V(20)*V(73)
  IRR(85) = RCT(85)*V(70)*V(73)
  IRR(86) = RCT(86)*V(33)*V(73)
  IRR(87) = RCT(87)*V(35)*V(73)
  IRR(88) = RCT(88)*V(43)*V(73)
  IRR(89) = RCT(89)*V(60)*V(73)
  IRR(90) = RCT(90)*V(56)*V(74)
  IRR(91) = RCT(91)*V(64)*V(74)
  IRR(92) = RCT(92)*V(44)*V(74)
  IRR(93) = RCT(93)*V(51)*V(74)
  IRR(94) = RCT(94)*V(13)*V(73)
  IRR(95) = RCT(95)*V(61)*V(74)
  IRR(96) = RCT(96)*V(40)*V(74)
  IRR(97) = RCT(97)*V(30)*V(74)
  IRR(98) = RCT(98)*V(63)*V(74)
  IRR(99) = RCT(99)*V(65)*V(74)
  IRR(100) = RCT(100)*V(41)*V(74)
  IRR(101) = RCT(101)*V(39)*V(74)
  IRR(102) = RCT(102)*V(36)*V(74)
  IRR(103) = RCT(103)*V(38)*V(74)
  IRR(104) = RCT(104)*V(43)*V(74)
  IRR(105) = RCT(105)*V(30)*V(77)
  IRR(106) = RCT(106)*V(63)*V(77)
  IRR(107) = RCT(107)*V(65)*V(77)
  IRR(108) = RCT(108)*V(41)*V(77)
  IRR(109) = RCT(109)*V(39)*V(77)
  IRR(110) = RCT(110)*V(36)*V(77)
  IRR(111) = RCT(111)*V(38)*V(77)
  IRR(112) = RCT(112)*V(48)*V(77)
  IRR(113) = RCT(113)*V(61)*V(77)
  IRR(114) = RCT(114)*V(43)*V(77)
  IRR(115) = RCT(115)*V(27)*V(80)
  IRR(116) = RCT(116)*V(27)*V(75)
  IRR(117) = RCT(117)*V(28)*V(80)
  IRR(118) = RCT(118)*V(28)*F(2)
  IRR(119) = RCT(119)*V(28)*V(77)
  IRR(120) = RCT(120)*V(29)*V(80)
  IRR(121) = RCT(121)*V(29)*F(2)
  IRR(122) = RCT(122)*V(29)*V(77)
  IRR(123) = RCT(123)*V(31)*V(80)
  IRR(124) = RCT(124)*V(31)*F(2)
  IRR(125) = RCT(125)*V(31)*V(77)
  IRR(126) = RCT(126)*V(72)*V(80)
  IRR(127) = RCT(127)*V(35)
  IRR(128) = RCT(128)*V(78)*V(80)
  IRR(129) = RCT(129)*V(43)
  IRR(130) = RCT(130)*V(76)*V(79)
  IRR(131) = RCT(131)*V(68)*V(79)
  IRR(132) = RCT(132)*V(71)*V(79)
  IRR(133) = RCT(133)*V(55)*V(79)
  IRR(134) = RCT(134)*V(53)*V(79)
  IRR(135) = RCT(135)*V(45)*V(79)
  IRR(136) = RCT(136)*V(46)*V(79)
  IRR(137) = RCT(137)*V(47)*V(79)
  IRR(138) = RCT(138)*V(50)*V(79)
  IRR(139) = RCT(139)*V(59)*V(79)
  IRR(140) = RCT(140)*V(54)*V(79)
  IRR(141) = RCT(141)*V(57)*V(79)
  IRR(142) = RCT(142)*V(58)*V(79)
  IRR(143) = RCT(143)*V(52)*V(79)
  IRR(144) = RCT(144)*V(72)*V(79)
  IRR(145) = RCT(145)*V(78)*V(79)
  IRR(146) = RCT(146)*V(69)*V(79)
  IRR(147) = RCT(147)*V(66)*V(79)
  IRR(148) = RCT(148)*V(67)*V(79)
  IRR(149) = RCT(149)*V(75)*V(76)
  IRR(150) = RCT(150)*V(68)*V(75)
  IRR(151) = RCT(151)*V(71)*V(75)
  IRR(152) = RCT(152)*V(55)*V(75)
  IRR(153) = RCT(153)*V(53)*V(75)
  IRR(154) = RCT(154)*V(45)*V(75)
  IRR(155) = RCT(155)*V(46)*V(75)
  IRR(156) = RCT(156)*V(47)*V(75)
  IRR(157) = RCT(157)*V(50)*V(75)
  IRR(158) = RCT(158)*V(59)*V(75)
  IRR(159) = RCT(159)*V(54)*V(75)
  IRR(160) = RCT(160)*V(57)*V(75)
  IRR(161) = RCT(161)*V(58)*V(75)
  IRR(162) = RCT(162)*V(52)*V(75)
  IRR(163) = RCT(163)*V(72)*V(75)
  IRR(164) = RCT(164)*V(72)*V(75)
  IRR(165) = RCT(165)*V(75)*V(78)
  IRR(166) = RCT(166)*V(75)*V(78)
  IRR(167) = RCT(167)*V(69)*V(75)
  IRR(168) = RCT(168)*V(66)*V(75)
  IRR(169) = RCT(169)*V(67)*V(75)
  IRR(170) = RCT(170)*V(76)*V(76)
  IRR(171) = RCT(171)*V(68)*V(76)
  IRR(172) = RCT(172)*V(71)*V(76)
  IRR(173) = RCT(173)*V(55)*V(76)
  IRR(174) = RCT(174)*V(53)*V(76)
  IRR(175) = RCT(175)*V(45)*V(76)
  IRR(176) = RCT(176)*V(46)*V(76)
  IRR(177) = RCT(177)*V(47)*V(76)
  IRR(178) = RCT(178)*V(50)*V(76)
  IRR(179) = RCT(179)*V(59)*V(76)
  IRR(180) = RCT(180)*V(54)*V(76)
  IRR(181) = RCT(181)*V(57)*V(76)
  IRR(182) = RCT(182)*V(58)*V(76)
  IRR(183) = RCT(183)*V(52)*V(76)
  IRR(184) = RCT(184)*V(72)*V(76)
  IRR(185) = RCT(185)*V(72)*V(76)
  IRR(186) = RCT(186)*V(76)*V(78)
  IRR(187) = RCT(187)*V(76)*V(78)
  IRR(188) = RCT(188)*V(69)*V(76)
  IRR(189) = RCT(189)*V(66)*V(76)
  IRR(190) = RCT(190)*V(67)*V(76)
  IRR(191) = RCT(191)*V(68)*V(72)
  IRR(192) = RCT(192)*V(71)*V(72)
  IRR(193) = RCT(193)*V(55)*V(72)
  IRR(194) = RCT(194)*V(53)*V(72)
  IRR(195) = RCT(195)*V(45)*V(72)
  IRR(196) = RCT(196)*V(46)*V(72)
  IRR(197) = RCT(197)*V(47)*V(72)
  IRR(198) = RCT(198)*V(50)*V(72)
  IRR(199) = RCT(199)*V(59)*V(72)
  IRR(200) = RCT(200)*V(54)*V(72)
  IRR(201) = RCT(201)*V(57)*V(72)
  IRR(202) = RCT(202)*V(58)*V(72)
  IRR(203) = RCT(203)*V(52)*V(72)
  IRR(204) = RCT(204)*V(72)*V(72)
  IRR(205) = RCT(205)*V(72)*V(78)
  IRR(206) = RCT(206)*V(69)*V(72)
  IRR(207) = RCT(207)*V(66)*V(72)
  IRR(208) = RCT(208)*V(67)*V(72)
  IRR(209) = RCT(209)*V(66)*V(66)
  IRR(210) = RCT(210)*V(66)*V(67)
  IRR(211) = RCT(211)*V(67)*V(67)
  IRR(212) = RCT(212)*V(74)*V(76)
  IRR(213) = RCT(213)*V(68)*V(74)
  IRR(214) = RCT(214)*V(71)*V(74)
  IRR(215) = RCT(215)*V(55)*V(74)
  IRR(216) = RCT(216)*V(53)*V(74)
  IRR(217) = RCT(217)*V(45)*V(74)
  IRR(218) = RCT(218)*V(46)*V(74)
  IRR(219) = RCT(219)*V(47)*V(74)
  IRR(220) = RCT(220)*V(21)*V(73)
  IRR(221) = RCT(221)*V(59)*V(74)
  IRR(222) = RCT(222)*V(54)*V(74)
  IRR(223) = RCT(223)*V(57)*V(74)
  IRR(224) = RCT(224)*V(58)*V(74)
  IRR(225) = RCT(225)*V(52)*V(74)
  IRR(226) = RCT(226)*V(72)*V(74)
  IRR(227) = RCT(227)*V(74)*V(78)
  IRR(228) = RCT(228)*V(69)*V(74)
  IRR(229) = RCT(229)*V(66)*V(74)
  IRR(230) = RCT(230)*V(67)*V(74)
  IRR(231) = RCT(231)*V(62)*V(75)
  IRR(232) = RCT(232)*V(62)*V(76)
  IRR(233) = RCT(233)*V(62)*V(72)
  IRR(234) = RCT(234)*V(62)*V(62)
  IRR(235) = RCT(235)*V(62)*V(79)
  IRR(236) = RCT(236)*V(62)*V(74)
  IRR(237) = RCT(237)*V(50)*V(50)
  IRR(238) = RCT(238)*V(7)*V(73)
  IRR(239) = RCT(239)*V(34)*V(73)
  IRR(240) = RCT(240)*V(42)*V(79)
  IRR(241) = RCT(241)*V(42)*V(75)
  IRR(242) = RCT(242)*V(42)*V(42)
  IRR(243) = RCT(243)*V(42)*V(80)
  IRR(244) = RCT(244)*V(21)
      
END SUBROUTINE racm_esrlsorg_IRRFun
















SUBROUTINE racm_esrlsorg_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(451)


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

  B(26) = RCT(25)*V(77)

  B(27) = RCT(25)*V(26)

  B(28) = RCT(26)*F(2)

  B(30) = RCT(27)*F(1)

  B(32) = RCT(28)*V(77)

  B(33) = RCT(28)*V(73)

  B(34) = RCT(29)*V(77)

  B(35) = RCT(29)*V(75)

  B(36) = RCT(30)*V(75)

  B(37) = RCT(30)*V(73)

  B(38) = RCT(31)*V(73)

  B(39) = RCT(31)*V(24)

  B(40) = RCT(32)*2*V(75)

  B(41) = RCT(33)*2*V(75)*F(1)

  B(43) = RCT(34)*V(79)

  B(44) = RCT(34)*V(26)

  B(45) = RCT(35)*V(80)

  B(46) = RCT(35)*V(26)

  B(47) = RCT(36)*V(80)

  B(48) = RCT(36)*V(26)

  B(49) = RCT(37)*V(79)

  B(50) = RCT(37)*V(73)

  B(51) = RCT(38)*V(80)

  B(52) = RCT(38)*V(73)

  B(53) = RCT(39)*V(74)

  B(54) = RCT(39)*V(73)

  B(55) = RCT(40)*V(79)

  B(56) = RCT(40)*V(75)

  B(57) = RCT(41)*V(80)

  B(58) = RCT(41)*V(75)

  B(59) = RCT(42)

  B(60) = RCT(43)*V(75)

  B(61) = RCT(43)*V(74)

  B(62) = RCT(44)*V(73)

  B(63) = RCT(44)*V(23)

  B(64) = RCT(45)*V(73)

  B(65) = RCT(45)*V(32)

  B(66) = RCT(46)*V(73)

  B(67) = RCT(46)*V(19)

  B(68) = RCT(47)*V(79)

  B(69) = RCT(47)*V(77)

  B(70) = RCT(48)*V(80)

  B(71) = RCT(48)*V(77)

  B(72) = RCT(49)*2*V(79)*F(2)

  B(74) = RCT(50)*V(79)

  B(75) = RCT(50)*V(74)

  B(76) = RCT(51)*V(80)

  B(77) = RCT(51)*V(74)

  B(78) = RCT(52)*V(80)

  B(79) = RCT(52)*V(74)

  B(80) = RCT(53)

  B(81) = RCT(54)*2*V(74)

  B(82) = RCT(55)*F(2)

  B(84) = RCT(56)*V(73)

  B(85) = RCT(56)*V(5)

  B(86) = RCT(57)*V(73)

  B(87) = RCT(57)*V(37)

  B(88) = RCT(58)*V(73)

  B(89) = RCT(58)*V(18)

  B(90) = RCT(59)*V(73)

  B(91) = RCT(59)*V(22)

  B(92) = RCT(60)*V(73)

  B(93) = RCT(60)*V(15)

  B(94) = RCT(61)*V(73)

  B(95) = RCT(61)*V(17)

  B(96) = RCT(62)*V(73)

  B(97) = RCT(62)*V(14)

  B(98) = RCT(63)*V(73)

  B(99) = RCT(63)*V(8)

  B(100) = RCT(64)*V(73)

  B(101) = RCT(64)*V(12)

  B(102) = RCT(65)*V(73)

  B(103) = RCT(65)*V(30)

  B(104) = RCT(66)*V(73)

  B(105) = RCT(66)*V(63)

  B(106) = RCT(67)*V(73)

  B(107) = RCT(67)*V(65)

  B(108) = RCT(68)*V(73)

  B(109) = RCT(68)*V(41)

  B(110) = RCT(69)*V(73)

  B(111) = RCT(69)*V(39)

  B(112) = RCT(70)*V(73)

  B(113) = RCT(70)*V(36)

  B(114) = RCT(71)*V(73)

  B(115) = RCT(71)*V(38)

  B(116) = RCT(72)*V(73)

  B(117) = RCT(72)*V(9)

  B(118) = RCT(73)*V(73)

  B(119) = RCT(73)*V(10)

  B(120) = RCT(74)*V(73)

  B(121) = RCT(74)*V(40)

  B(122) = RCT(75)*V(73)

  B(123) = RCT(75)*V(56)

  B(124) = RCT(76)*V(73)

  B(125) = RCT(76)*V(64)

  B(126) = RCT(77)*V(73)

  B(127) = RCT(77)*V(49)

  B(128) = RCT(78)*V(73)

  B(129) = RCT(78)*V(25)

  B(130) = RCT(79)*V(73)

  B(131) = RCT(79)*V(44)

  B(132) = RCT(80)*V(73)

  B(133) = RCT(80)*V(51)

  B(134) = RCT(81)*V(73)

  B(135) = RCT(81)*V(48)

  B(136) = RCT(82)*V(73)

  B(137) = RCT(82)*V(61)

  B(138) = RCT(83)*V(73)

  B(139) = RCT(83)*V(16)

  B(140) = RCT(84)*V(73)

  B(141) = RCT(84)*V(20)

  B(142) = RCT(85)*V(73)

  B(143) = RCT(85)*V(70)

  B(144) = RCT(86)*V(73)

  B(145) = RCT(86)*V(33)

  B(146) = RCT(87)*V(73)

  B(147) = RCT(87)*V(35)

  B(148) = RCT(88)*V(73)

  B(149) = RCT(88)*V(43)

  B(150) = RCT(89)*V(73)

  B(151) = RCT(89)*V(60)

  B(152) = RCT(90)*V(74)

  B(153) = RCT(90)*V(56)

  B(154) = RCT(91)*V(74)

  B(155) = RCT(91)*V(64)

  B(156) = RCT(92)*V(74)

  B(157) = RCT(92)*V(44)

  B(158) = RCT(93)*V(74)

  B(159) = RCT(93)*V(51)

  B(160) = RCT(94)*V(73)

  B(161) = RCT(94)*V(13)

  B(162) = RCT(95)*V(74)

  B(163) = RCT(95)*V(61)

  B(164) = RCT(96)*V(74)

  B(165) = RCT(96)*V(40)

  B(166) = RCT(97)*V(74)

  B(167) = RCT(97)*V(30)

  B(168) = RCT(98)*V(74)

  B(169) = RCT(98)*V(63)

  B(170) = RCT(99)*V(74)

  B(171) = RCT(99)*V(65)

  B(172) = RCT(100)*V(74)

  B(173) = RCT(100)*V(41)

  B(174) = RCT(101)*V(74)

  B(175) = RCT(101)*V(39)

  B(176) = RCT(102)*V(74)

  B(177) = RCT(102)*V(36)

  B(178) = RCT(103)*V(74)

  B(179) = RCT(103)*V(38)

  B(180) = RCT(104)*V(74)

  B(181) = RCT(104)*V(43)

  B(182) = RCT(105)*V(77)

  B(183) = RCT(105)*V(30)

  B(184) = RCT(106)*V(77)

  B(185) = RCT(106)*V(63)

  B(186) = RCT(107)*V(77)

  B(187) = RCT(107)*V(65)

  B(188) = RCT(108)*V(77)

  B(189) = RCT(108)*V(41)

  B(190) = RCT(109)*V(77)

  B(191) = RCT(109)*V(39)

  B(192) = RCT(110)*V(77)

  B(193) = RCT(110)*V(36)

  B(194) = RCT(111)*V(77)

  B(195) = RCT(111)*V(38)

  B(196) = RCT(112)*V(77)

  B(197) = RCT(112)*V(48)

  B(198) = RCT(113)*V(77)

  B(199) = RCT(113)*V(61)

  B(200) = RCT(114)*V(77)

  B(201) = RCT(114)*V(43)

  B(202) = RCT(115)*V(80)

  B(203) = RCT(115)*V(27)

  B(204) = RCT(116)*V(75)

  B(205) = RCT(116)*V(27)

  B(206) = RCT(117)*V(80)

  B(207) = RCT(117)*V(28)

  B(208) = RCT(118)*F(2)

  B(210) = RCT(119)*V(77)

  B(211) = RCT(119)*V(28)

  B(212) = RCT(120)*V(80)

  B(213) = RCT(120)*V(29)

  B(214) = RCT(121)*F(2)

  B(216) = RCT(122)*V(77)

  B(217) = RCT(122)*V(29)

  B(218) = RCT(123)*V(80)

  B(219) = RCT(123)*V(31)

  B(220) = RCT(124)*F(2)

  B(222) = RCT(125)*V(77)

  B(223) = RCT(125)*V(31)

  B(224) = RCT(126)*V(80)

  B(225) = RCT(126)*V(72)

  B(226) = RCT(127)

  B(227) = RCT(128)*V(80)

  B(228) = RCT(128)*V(78)

  B(229) = RCT(129)

  B(230) = RCT(130)*V(79)

  B(231) = RCT(130)*V(76)

  B(232) = RCT(131)*V(79)

  B(233) = RCT(131)*V(68)

  B(234) = RCT(132)*V(79)

  B(235) = RCT(132)*V(71)

  B(236) = RCT(133)*V(79)

  B(237) = RCT(133)*V(55)

  B(238) = RCT(134)*V(79)

  B(239) = RCT(134)*V(53)

  B(240) = RCT(135)*V(79)

  B(241) = RCT(135)*V(45)

  B(242) = RCT(136)*V(79)

  B(243) = RCT(136)*V(46)

  B(244) = RCT(137)*V(79)

  B(245) = RCT(137)*V(47)

  B(246) = RCT(138)*V(79)

  B(247) = RCT(138)*V(50)

  B(248) = RCT(139)*V(79)

  B(249) = RCT(139)*V(59)

  B(250) = RCT(140)*V(79)

  B(251) = RCT(140)*V(54)

  B(252) = RCT(141)*V(79)

  B(253) = RCT(141)*V(57)

  B(254) = RCT(142)*V(79)

  B(255) = RCT(142)*V(58)

  B(256) = RCT(143)*V(79)

  B(257) = RCT(143)*V(52)

  B(258) = RCT(144)*V(79)

  B(259) = RCT(144)*V(72)

  B(260) = RCT(145)*V(79)

  B(261) = RCT(145)*V(78)

  B(262) = RCT(146)*V(79)

  B(263) = RCT(146)*V(69)

  B(264) = RCT(147)*V(79)

  B(265) = RCT(147)*V(66)

  B(266) = RCT(148)*V(79)

  B(267) = RCT(148)*V(67)

  B(268) = RCT(149)*V(76)

  B(269) = RCT(149)*V(75)

  B(270) = RCT(150)*V(75)

  B(271) = RCT(150)*V(68)

  B(272) = RCT(151)*V(75)

  B(273) = RCT(151)*V(71)

  B(274) = RCT(152)*V(75)

  B(275) = RCT(152)*V(55)

  B(276) = RCT(153)*V(75)

  B(277) = RCT(153)*V(53)

  B(278) = RCT(154)*V(75)

  B(279) = RCT(154)*V(45)

  B(280) = RCT(155)*V(75)

  B(281) = RCT(155)*V(46)

  B(282) = RCT(156)*V(75)

  B(283) = RCT(156)*V(47)

  B(284) = RCT(157)*V(75)

  B(285) = RCT(157)*V(50)

  B(286) = RCT(158)*V(75)

  B(287) = RCT(158)*V(59)

  B(288) = RCT(159)*V(75)

  B(289) = RCT(159)*V(54)

  B(290) = RCT(160)*V(75)

  B(291) = RCT(160)*V(57)

  B(292) = RCT(161)*V(75)

  B(293) = RCT(161)*V(58)

  B(294) = RCT(162)*V(75)

  B(295) = RCT(162)*V(52)

  B(296) = RCT(163)*V(75)

  B(297) = RCT(163)*V(72)

  B(298) = RCT(164)*V(75)

  B(299) = RCT(164)*V(72)

  B(300) = RCT(165)*V(78)

  B(301) = RCT(165)*V(75)

  B(302) = RCT(166)*V(78)

  B(303) = RCT(166)*V(75)

  B(304) = RCT(167)*V(75)

  B(305) = RCT(167)*V(69)

  B(306) = RCT(168)*V(75)

  B(307) = RCT(168)*V(66)

  B(308) = RCT(169)*V(75)

  B(309) = RCT(169)*V(67)

  B(310) = RCT(170)*2*V(76)

  B(311) = RCT(171)*V(76)

  B(312) = RCT(171)*V(68)

  B(313) = RCT(172)*V(76)

  B(314) = RCT(172)*V(71)

  B(315) = RCT(173)*V(76)

  B(316) = RCT(173)*V(55)

  B(317) = RCT(174)*V(76)

  B(318) = RCT(174)*V(53)

  B(319) = RCT(175)*V(76)

  B(320) = RCT(175)*V(45)

  B(321) = RCT(176)*V(76)

  B(322) = RCT(176)*V(46)

  B(323) = RCT(177)*V(76)

  B(324) = RCT(177)*V(47)

  B(325) = RCT(178)*V(76)

  B(326) = RCT(178)*V(50)

  B(327) = RCT(179)*V(76)

  B(328) = RCT(179)*V(59)

  B(329) = RCT(180)*V(76)

  B(330) = RCT(180)*V(54)

  B(331) = RCT(181)*V(76)

  B(332) = RCT(181)*V(57)

  B(333) = RCT(182)*V(76)

  B(334) = RCT(182)*V(58)

  B(335) = RCT(183)*V(76)

  B(336) = RCT(183)*V(52)

  B(337) = RCT(184)*V(76)

  B(338) = RCT(184)*V(72)

  B(339) = RCT(185)*V(76)

  B(340) = RCT(185)*V(72)

  B(341) = RCT(186)*V(78)

  B(342) = RCT(186)*V(76)

  B(343) = RCT(187)*V(78)

  B(344) = RCT(187)*V(76)

  B(345) = RCT(188)*V(76)

  B(346) = RCT(188)*V(69)

  B(347) = RCT(189)*V(76)

  B(348) = RCT(189)*V(66)

  B(349) = RCT(190)*V(76)

  B(350) = RCT(190)*V(67)

  B(351) = RCT(191)*V(72)

  B(352) = RCT(191)*V(68)

  B(353) = RCT(192)*V(72)

  B(354) = RCT(192)*V(71)

  B(355) = RCT(193)*V(72)

  B(356) = RCT(193)*V(55)

  B(357) = RCT(194)*V(72)

  B(358) = RCT(194)*V(53)

  B(359) = RCT(195)*V(72)

  B(360) = RCT(195)*V(45)

  B(361) = RCT(196)*V(72)

  B(362) = RCT(196)*V(46)

  B(363) = RCT(197)*V(72)

  B(364) = RCT(197)*V(47)

  B(365) = RCT(198)*V(72)

  B(366) = RCT(198)*V(50)

  B(367) = RCT(199)*V(72)

  B(368) = RCT(199)*V(59)

  B(369) = RCT(200)*V(72)

  B(370) = RCT(200)*V(54)

  B(371) = RCT(201)*V(72)

  B(372) = RCT(201)*V(57)

  B(373) = RCT(202)*V(72)

  B(374) = RCT(202)*V(58)

  B(375) = RCT(203)*V(72)

  B(376) = RCT(203)*V(52)

  B(377) = RCT(204)*2*V(72)

  B(378) = RCT(205)*V(78)

  B(379) = RCT(205)*V(72)

  B(380) = RCT(206)*V(72)

  B(381) = RCT(206)*V(69)

  B(382) = RCT(207)*V(72)

  B(383) = RCT(207)*V(66)

  B(384) = RCT(208)*V(72)

  B(385) = RCT(208)*V(67)

  B(386) = RCT(209)*2*V(66)

  B(387) = RCT(210)*V(67)

  B(388) = RCT(210)*V(66)

  B(389) = RCT(211)*2*V(67)

  B(390) = RCT(212)*V(76)

  B(391) = RCT(212)*V(74)

  B(392) = RCT(213)*V(74)

  B(393) = RCT(213)*V(68)

  B(394) = RCT(214)*V(74)

  B(395) = RCT(214)*V(71)

  B(396) = RCT(215)*V(74)

  B(397) = RCT(215)*V(55)

  B(398) = RCT(216)*V(74)

  B(399) = RCT(216)*V(53)

  B(400) = RCT(217)*V(74)

  B(401) = RCT(217)*V(45)

  B(402) = RCT(218)*V(74)

  B(403) = RCT(218)*V(46)

  B(404) = RCT(219)*V(74)

  B(405) = RCT(219)*V(47)

  B(406) = RCT(220)*V(73)

  B(407) = RCT(220)*V(21)

  B(408) = RCT(221)*V(74)

  B(409) = RCT(221)*V(59)

  B(410) = RCT(222)*V(74)

  B(411) = RCT(222)*V(54)

  B(412) = RCT(223)*V(74)

  B(413) = RCT(223)*V(57)

  B(414) = RCT(224)*V(74)

  B(415) = RCT(224)*V(58)

  B(416) = RCT(225)*V(74)

  B(417) = RCT(225)*V(52)

  B(418) = RCT(226)*V(74)

  B(419) = RCT(226)*V(72)

  B(420) = RCT(227)*V(78)

  B(421) = RCT(227)*V(74)

  B(422) = RCT(228)*V(74)

  B(423) = RCT(228)*V(69)

  B(424) = RCT(229)*V(74)

  B(425) = RCT(229)*V(66)

  B(426) = RCT(230)*V(74)

  B(427) = RCT(230)*V(67)

  B(428) = RCT(231)*V(75)

  B(429) = RCT(231)*V(62)

  B(430) = RCT(232)*V(76)

  B(431) = RCT(232)*V(62)

  B(432) = RCT(233)*V(72)

  B(433) = RCT(233)*V(62)

  B(434) = RCT(234)*2*V(62)

  B(435) = RCT(235)*V(79)

  B(436) = RCT(235)*V(62)

  B(437) = RCT(236)*V(74)

  B(438) = RCT(236)*V(62)

  B(439) = RCT(237)*2*V(50)

  B(440) = RCT(238)*V(73)

  B(441) = RCT(238)*V(7)

  B(442) = RCT(239)*V(73)

  B(443) = RCT(239)*V(34)

  B(444) = RCT(240)*V(79)

  B(445) = RCT(240)*V(42)

  B(446) = RCT(241)*V(75)

  B(447) = RCT(241)*V(42)

  B(448) = RCT(242)*2*V(42)

  B(449) = RCT(243)*V(80)

  B(450) = RCT(243)*V(42)

  B(451) = RCT(244)



  JVS(1) = 0

  JVS(2) = B(84)

  JVS(3) = B(85)

  JVS(4) = 0

  JVS(5) = B(86)

  JVS(6) = B(87)

  JVS(7) = 0

  JVS(8) = 0.036*B(96)

  JVS(9) = 0.37*B(182)

  JVS(10) = 0.01*B(194)

  JVS(11) = 0.28*B(190)

  JVS(12) = 0.15*B(188)

  JVS(13) = 0.11*B(200)

  JVS(14) = 0.45*B(196)

  JVS(15) = 0.11*B(198)

  JVS(16) = 0.14*B(184)

  JVS(17) = 0.036*B(97)

  JVS(18) = 0.37*B(183)+0.14*B(185)+0.15*B(189)+0.28*B(191)+0.01*B(195)+0.45*B(197)+0.11*B(199)+0.11*B(201)

  JVS(19) = 0

  JVS(20) = 0.07*B(194)

  JVS(21) = 0.5*B(359)

  JVS(22) = 0.499*B(361)

  JVS(23) = 0.49*B(363)

  JVS(24) = 0.494*B(365)

  JVS(25) = 0.495*B(357)

  JVS(26) = 0.495*B(355)

  JVS(27) = 0.21*B(198)

  JVS(28) = 0.1*B(184)

  JVS(29) = 0.14*B(186)

  JVS(30) = 0.5*B(382)

  JVS(31) = 0.484*B(384)

  JVS(32) = 0.5*B(351)

  JVS(33) = 0.5*B(380)

  JVS(34) = 0.499*B(353)

  JVS(35) = B(298)+B(339)+0.5*B(352)+0.499*B(354)+0.495*B(356)+0.495*B(358)+0.5*B(360)+0.499*B(362)+0.49*B(364)+0.494&
              &*B(366)+0.5*B(381)+0.5*B(383)+0.484*B(385)

  JVS(36) = B(299)+B(302)

  JVS(37) = B(340)+B(343)

  JVS(38) = 0.1*B(185)+0.14*B(187)+0.07*B(195)+0.21*B(199)

  JVS(39) = B(303)+B(344)

  JVS(40) = -B(84)

  JVS(41) = -B(85)

  JVS(42) = -B(28)-B(30)

  JVS(43) = B(2)

  JVS(44) = -B(440)

  JVS(45) = B(284)

  JVS(46) = -B(441)

  JVS(47) = B(285)

  JVS(48) = -B(98)

  JVS(49) = -B(99)

  JVS(50) = -B(116)

  JVS(51) = -B(117)

  JVS(52) = -B(118)

  JVS(53) = -B(119)

  JVS(54) = -B(80)

  JVS(55) = B(78)

  JVS(56) = B(79)

  JVS(57) = -B(100)

  JVS(58) = -B(101)

  JVS(59) = -B(160)

  JVS(60) = B(446)

  JVS(61) = -B(161)

  JVS(62) = B(447)

  JVS(63) = -B(96)

  JVS(64) = -B(97)

  JVS(65) = -B(92)

  JVS(66) = 0.06*B(184)

  JVS(67) = 0.07*B(186)

  JVS(68) = -B(93)

  JVS(69) = 0.06*B(185)+0.07*B(187)

  JVS(70) = -B(138)

  JVS(71) = 0.35*B(136)

  JVS(72) = 0.35*B(137)-B(139)

  JVS(73) = -B(94)

  JVS(74) = 0.03*B(184)

  JVS(75) = 0.06*B(186)

  JVS(76) = -B(95)

  JVS(77) = 0.03*B(185)+0.06*B(187)

  JVS(78) = -B(88)

  JVS(79) = B(442)

  JVS(80) = -B(89)+B(443)

  JVS(81) = -B(6)-B(59)-B(66)

  JVS(82) = -B(67)

  JVS(83) = B(57)

  JVS(84) = B(58)

  JVS(85) = -B(13)-B(140)

  JVS(86) = -B(141)

  JVS(87) = B(268)

  JVS(88) = B(269)

  JVS(89) = -B(406)-B(451)

  JVS(90) = B(449)

  JVS(91) = -B(407)

  JVS(92) = B(450)

  JVS(93) = B(406)

  JVS(94) = -B(90)

  JVS(95) = B(442)

  JVS(96) = 0.25*B(444)+B(448)

  JVS(97) = -B(91)+B(407)+B(443)

  JVS(98) = 0.25*B(445)

  JVS(99) = 0

  JVS(100) = -B(4)-B(62)

  JVS(101) = B(206)

  JVS(102) = B(212)

  JVS(103) = B(218)

  JVS(104) = B(49)-B(63)

  JVS(105) = B(50)

  JVS(106) = B(207)+B(213)+B(219)

  JVS(107) = -B(9)-B(38)

  JVS(108) = 0.02*B(192)

  JVS(109) = 0.02*B(194)

  JVS(110) = 0.09*B(190)

  JVS(111) = 0.001*B(188)

  JVS(112) = 0.006*B(184)

  JVS(113) = 0.011*B(186)

  JVS(114) = -B(39)

  JVS(115) = B(40)+B(41)

  JVS(116) = 0.006*B(185)+0.011*B(187)+0.001*B(189)+0.09*B(191)+0.02*B(193)+0.02*B(195)

  JVS(117) = 0.024*B(100)

  JVS(118) = -B(23)-B(128)

  JVS(119) = 0.6*B(148)

  JVS(120) = 0.3*B(345)

  JVS(121) = 0.024*B(101)-B(129)+0.6*B(149)

  JVS(122) = 0.3*B(346)

  JVS(123) = B(28)

  JVS(124) = -B(24)-B(26)-B(43)-B(45)-B(47)

  JVS(125) = 0.09*B(188)

  JVS(126) = B(8)

  JVS(127) = B(3)-B(27)+0.09*B(189)

  JVS(128) = -B(44)

  JVS(129) = B(1)-B(46)-B(48)

  JVS(130) = -B(202)-B(204)

  JVS(131) = 0.1*B(120)+B(164)

  JVS(132) = 0.1*B(121)

  JVS(133) = B(165)

  JVS(134) = -B(205)

  JVS(135) = -B(203)

  JVS(136) = 0.9*B(116)

  JVS(137) = -B(206)-B(208)-B(210)

  JVS(138) = 0.9*B(117)

  JVS(139) = -B(211)

  JVS(140) = -B(207)

  JVS(141) = 0.9*B(118)

  JVS(142) = -B(212)-B(214)-B(216)

  JVS(143) = 0.9*B(119)

  JVS(144) = -B(217)

  JVS(145) = -B(213)

  JVS(146) = -B(102)-B(166)-B(182)

  JVS(147) = -B(103)

  JVS(148) = -B(167)

  JVS(149) = -B(183)

  JVS(150) = -B(218)-B(220)-B(222)

  JVS(151) = 0.85*B(120)

  JVS(152) = 0.85*B(121)

  JVS(153) = -B(223)

  JVS(154) = -B(219)

  JVS(155) = -B(5)-B(64)

  JVS(156) = B(164)

  JVS(157) = B(156)

  JVS(158) = B(158)

  JVS(159) = B(152)

  JVS(160) = 0.5*B(162)

  JVS(161) = B(154)

  JVS(162) = B(51)-B(65)

  JVS(163) = 0.3*B(60)+B(153)+B(155)+B(157)+B(159)+0.5*B(163)+B(165)

  JVS(164) = 0.3*B(61)

  JVS(165) = B(52)

  JVS(166) = -B(15)-B(144)

  JVS(167) = 0.11*B(198)

  JVS(168) = B(296)

  JVS(169) = -B(145)

  JVS(170) = B(297)

  JVS(171) = 0.11*B(199)

  JVS(172) = -B(442)

  JVS(173) = B(174)

  JVS(174) = 0.046*B(246)

  JVS(175) = -B(443)

  JVS(176) = B(175)

  JVS(177) = 0.046*B(247)

  JVS(178) = -B(146)-B(226)

  JVS(179) = 0.4*B(148)+0.4*B(180)+0.3*B(200)

  JVS(180) = B(224)

  JVS(181) = -B(147)+0.4*B(149)

  JVS(182) = 0.4*B(181)

  JVS(183) = 0.3*B(201)

  JVS(184) = B(225)

  JVS(185) = -B(112)-B(176)-B(192)

  JVS(186) = -B(113)

  JVS(187) = -B(177)

  JVS(188) = -B(193)

  JVS(189) = 0.036*B(96)

  JVS(190) = B(88)

  JVS(191) = 0.43*B(182)

  JVS(192) = 0

  JVS(193) = 0.14*B(192)

  JVS(194) = -B(86)

  JVS(195) = 0.14*B(194)

  JVS(196) = 0.14*B(190)

  JVS(197) = 0.36*B(188)

  JVS(198) = 0.25*B(444)+0.5*B(448)

  JVS(199) = 0.13*B(200)

  JVS(200) = 1.87*B(17)+1.55*B(18)+2*B(130)+2*B(156)

  JVS(201) = B(22)+0.22*B(196)

  JVS(202) = 0

  JVS(203) = B(19)+B(132)+B(158)

  JVS(204) = B(10)+B(11)+B(122)+B(152)

  JVS(205) = 0.66*B(198)

  JVS(206) = 0.37*B(184)

  JVS(207) = B(12)

  JVS(208) = 0.3*B(186)

  JVS(209) = -B(87)+B(89)+0.036*B(97)+B(123)+2*B(131)+B(133)

  JVS(210) = B(153)+2*B(157)+B(159)

  JVS(211) = 0.43*B(183)+0.37*B(185)+0.3*B(187)+0.36*B(189)+0.14*B(191)+0.14*B(193)+0.14*B(195)+0.22*B(197)+0.66*B(199)&
               &+0.13*B(201)

  JVS(212) = 0.25*B(445)

  JVS(213) = -B(114)-B(178)-B(194)

  JVS(214) = -B(115)

  JVS(215) = -B(179)

  JVS(216) = -B(195)

  JVS(217) = -B(110)-B(174)-B(190)

  JVS(218) = -B(111)

  JVS(219) = -B(175)

  JVS(220) = -B(191)

  JVS(221) = 0.1*B(202)+B(204)

  JVS(222) = B(206)+0.02*B(208)+B(210)

  JVS(223) = B(212)+0.02*B(214)+B(216)

  JVS(224) = B(218)+0.02*B(220)+B(222)

  JVS(225) = -B(120)-B(164)

  JVS(226) = -B(121)

  JVS(227) = -B(165)

  JVS(228) = B(205)

  JVS(229) = B(211)+B(217)+B(223)

  JVS(230) = 0.1*B(203)+B(207)+B(213)+B(219)

  JVS(231) = -B(108)-B(172)-B(188)

  JVS(232) = -B(109)

  JVS(233) = -B(173)

  JVS(234) = -B(189)

  JVS(235) = B(160)

  JVS(236) = B(451)

  JVS(237) = 0.1*B(190)

  JVS(238) = -B(444)-B(446)-2*B(448)-B(449)

  JVS(239) = B(134)

  JVS(240) = B(135)+B(161)

  JVS(241) = 0

  JVS(242) = -B(447)

  JVS(243) = 0.1*B(191)

  JVS(244) = -B(445)

  JVS(245) = -B(450)

  JVS(246) = -B(148)-B(180)-B(200)-B(229)

  JVS(247) = -B(149)

  JVS(248) = -B(181)

  JVS(249) = -B(201)

  JVS(250) = B(227)

  JVS(251) = B(228)

  JVS(252) = 0.036*B(96)

  JVS(253) = -B(17)-B(18)-B(130)-B(156)

  JVS(254) = B(256)+B(335)+B(375)+B(416)

  JVS(255) = 1.2*B(252)+0.65*B(331)+0.65*B(371)+1.3*B(412)

  JVS(256) = 0.35*B(254)+0.37*B(333)+0.37*B(373)+0.74*B(414)

  JVS(257) = 0.15*B(136)+0.25*B(162)+0.5*B(198)

  JVS(258) = 0.063*B(234)+0.119*B(313)+0.1*B(353)+0.063*B(394)

  JVS(259) = 0.1*B(354)+0.65*B(372)+0.37*B(374)+B(376)

  JVS(260) = 0.036*B(97)-B(131)+0.15*B(137)

  JVS(261) = -B(157)+0.25*B(163)+0.063*B(395)+1.3*B(413)+0.74*B(415)+B(417)

  JVS(262) = 0.119*B(314)+0.65*B(332)+0.37*B(334)+B(336)

  JVS(263) = 0.5*B(199)

  JVS(264) = 0.063*B(235)+1.2*B(253)+0.35*B(255)+B(257)

  JVS(265) = B(102)

  JVS(266) = -B(240)-B(278)-B(319)-B(359)-B(400)

  JVS(267) = -B(360)

  JVS(268) = B(103)

  JVS(269) = -B(401)

  JVS(270) = -B(279)

  JVS(271) = -B(320)

  JVS(272) = 0

  JVS(273) = -B(241)

  JVS(274) = -B(242)-B(280)-B(321)-B(361)-B(402)

  JVS(275) = B(104)

  JVS(276) = -B(362)

  JVS(277) = B(105)

  JVS(278) = -B(403)

  JVS(279) = -B(281)

  JVS(280) = -B(322)

  JVS(281) = -B(243)

  JVS(282) = -B(244)-B(282)-B(323)-B(363)-B(404)

  JVS(283) = B(106)

  JVS(284) = -B(364)

  JVS(285) = B(107)

  JVS(286) = -B(405)

  JVS(287) = -B(283)

  JVS(288) = -B(324)

  JVS(289) = -B(245)

  JVS(290) = B(440)

  JVS(291) = 0.79*B(194)

  JVS(292) = 0.65*B(190)

  JVS(293) = 0.9*B(172)+0.39*B(188)

  JVS(294) = -B(22)-B(134)-B(196)

  JVS(295) = B(246)+0.55*B(325)+0.771*B(365)+2*B(439)

  JVS(296) = 0.4*B(250)+0.6*B(329)+0.6*B(369)+0.6*B(410)

  JVS(297) = 0.771*B(366)+0.6*B(370)

  JVS(298) = -B(135)+B(441)

  JVS(299) = 0.9*B(173)+0.6*B(411)

  JVS(300) = 0

  JVS(301) = 0.55*B(326)+0.6*B(330)

  JVS(302) = 0.39*B(189)+0.65*B(191)+0.79*B(195)-B(197)

  JVS(303) = B(247)+0.4*B(251)

  JVS(304) = 0.25*B(98)

  JVS(305) = 0.12*B(138)

  JVS(306) = 0.53*B(192)

  JVS(307) = 0.06*B(242)+0.081*B(321)+0.141*B(361)+0.06*B(402)

  JVS(308) = 0.29*B(244)+0.313*B(323)+0.569*B(363)+0.29*B(404)

  JVS(309) = -B(16)-B(126)

  JVS(310) = 0.642*B(238)+0.419*B(317)+0.581*B(357)+0.88*B(398)

  JVS(311) = 0.722*B(236)+0.24*B(315)+0.33*B(355)+0.828*B(396)

  JVS(312) = 0.8*B(248)+B(327)+B(367)+B(408)

  JVS(313) = 0.8*B(21)

  JVS(314) = 0.03*B(162)

  JVS(315) = 0.03*B(184)

  JVS(316) = 0.16*B(186)

  JVS(317) = 0.149*B(387)

  JVS(318) = 0.464*B(266)+0.149*B(349)+0.167*B(384)+0.149*B(388)+0.285*B(389)+0.469*B(426)

  JVS(319) = 0.11*B(380)

  JVS(320) = 0.41*B(142)

  JVS(321) = 0.623*B(234)+0.018*B(313)+0.127*B(353)+0.67*B(394)

  JVS(322) = 0.127*B(354)+0.33*B(356)+0.581*B(358)+0.141*B(362)+0.569*B(364)+B(368)+0.11*B(381)+0.167*B(385)

  JVS(323) = 0.25*B(99)-B(127)+0.12*B(139)+0.41*B(143)

  JVS(324) = 0.03*B(163)+0.67*B(395)+0.828*B(397)+0.88*B(399)+0.06*B(403)+0.29*B(405)+B(409)+0.469*B(427)

  JVS(325) = 0

  JVS(326) = 0.018*B(314)+0.24*B(316)+0.419*B(318)+0.081*B(322)+0.313*B(324)+B(328)+0.149*B(350)

  JVS(327) = 0.03*B(185)+0.16*B(187)+0.53*B(193)

  JVS(328) = 0.623*B(235)+0.722*B(237)+0.642*B(239)+0.06*B(243)+0.29*B(245)+0.8*B(249)+0.464*B(267)

  JVS(329) = B(110)

  JVS(330) = B(108)

  JVS(331) = -B(246)-B(284)-B(325)-B(365)-2*B(439)

  JVS(332) = -B(366)

  JVS(333) = B(109)+B(111)

  JVS(334) = 0

  JVS(335) = -B(285)

  JVS(336) = -B(326)

  JVS(337) = 0

  JVS(338) = -B(247)

  JVS(339) = B(90)

  JVS(340) = B(128)

  JVS(341) = 0

  JVS(342) = 0

  JVS(343) = 0.5*B(444)+B(448)

  JVS(344) = 0

  JVS(345) = 0.9*B(196)

  JVS(346) = 0

  JVS(347) = -B(19)-B(132)-B(158)

  JVS(348) = B(256)+B(335)+B(375)+B(416)

  JVS(349) = 0

  JVS(350) = 0.65*B(252)+0.35*B(331)+0.35*B(371)+0.7*B(412)

  JVS(351) = 0.6*B(254)+0.63*B(333)+0.63*B(373)+1.26*B(414)

  JVS(352) = 0.15*B(136)+0.25*B(162)+0.62*B(198)

  JVS(353) = 0.54*B(262)+0.4*B(345)+0.54*B(380)+0.54*B(422)

  JVS(354) = 0.005*B(313)+0.004*B(353)

  JVS(355) = 0.004*B(354)+0.35*B(372)+0.63*B(374)+B(376)+0.54*B(381)

  JVS(356) = B(91)+B(129)-B(133)+0.15*B(137)

  JVS(357) = -B(159)+0.25*B(163)+0.7*B(413)+1.26*B(415)+B(417)+0.54*B(423)

  JVS(358) = 0

  JVS(359) = 0.005*B(314)+0.35*B(332)+0.63*B(334)+B(336)+0.4*B(346)

  JVS(360) = 0.9*B(197)+0.62*B(199)

  JVS(361) = 0

  JVS(362) = 0.65*B(253)+0.6*B(255)+B(257)+0.54*B(263)+0.5*B(445)

  JVS(363) = 0

  JVS(364) = 0.98*B(220)

  JVS(365) = 0

  JVS(366) = -B(256)-B(294)-B(335)-B(375)-B(416)

  JVS(367) = -B(376)

  JVS(368) = 0

  JVS(369) = -B(417)

  JVS(370) = -B(295)

  JVS(371) = -B(336)

  JVS(372) = 0

  JVS(373) = -B(257)

  JVS(374) = 0

  JVS(375) = 0.9511*B(100)

  JVS(376) = -B(238)-B(276)-B(317)-B(357)-B(398)

  JVS(377) = -B(358)

  JVS(378) = 0.9511*B(101)

  JVS(379) = -B(399)

  JVS(380) = -B(277)

  JVS(381) = -B(318)

  JVS(382) = -B(239)

  JVS(383) = B(114)

  JVS(384) = -B(250)-B(288)-B(329)-B(369)-B(410)

  JVS(385) = -B(370)

  JVS(386) = B(115)

  JVS(387) = -B(411)

  JVS(388) = -B(289)

  JVS(389) = -B(330)

  JVS(390) = 0

  JVS(391) = -B(251)

  JVS(392) = 0.75*B(98)

  JVS(393) = -B(236)-B(274)-B(315)-B(355)-B(396)

  JVS(394) = -B(356)

  JVS(395) = 0.75*B(99)

  JVS(396) = -B(397)

  JVS(397) = -B(275)

  JVS(398) = -B(316)

  JVS(399) = -B(237)

  JVS(400) = 0.01*B(96)

  JVS(401) = B(88)

  JVS(402) = B(13)+0.35*B(140)

  JVS(403) = B(23)

  JVS(404) = B(182)

  JVS(405) = 0.35*B(144)

  JVS(406) = 0

  JVS(407) = B(146)

  JVS(408) = 0.04*B(194)

  JVS(409) = 0.58*B(190)

  JVS(410) = 0.9*B(188)

  JVS(411) = 0.75*B(444)+0.5*B(448)

  JVS(412) = 0.4*B(148)+0.4*B(180)+0.7*B(200)

  JVS(413) = 0.13*B(17)+0.45*B(18)

  JVS(414) = 1.6*B(240)+1.55*B(319)+0.8*B(359)+1.6*B(400)

  JVS(415) = B(242)+1.25*B(321)+0.501*B(361)+B(402)

  JVS(416) = 0.755*B(323)

  JVS(417) = B(22)

  JVS(418) = B(246)+1.09*B(325)+0.34*B(365)+B(439)

  JVS(419) = B(335)

  JVS(420) = 0.753*B(317)

  JVS(421) = 0.25*B(250)+1.4*B(329)+0.4*B(369)+0.4*B(410)

  JVS(422) = 0.021*B(236)+0.829*B(315)+0.076*B(355)+0.021*B(396)

  JVS(423) = -B(10)-B(11)-B(122)-B(152)

  JVS(424) = B(331)

  JVS(425) = B(333)

  JVS(426) = B(327)

  JVS(427) = 0

  JVS(428) = B(430)

  JVS(429) = 0.64*B(184)

  JVS(430) = 0.02*B(186)

  JVS(431) = 0.75*B(347)+0.202*B(387)

  JVS(432) = 0.287*B(266)+0.96*B(349)+0.207*B(384)+0.202*B(388)+0.504*B(389)+0.28*B(426)

  JVS(433) = 0.75*B(311)

  JVS(434) = 0.75*B(345)

  JVS(435) = 0.047*B(234)+0.81*B(313)+0.091*B(353)+0.048*B(394)

  JVS(436) = B(337)+B(339)+0.091*B(354)+0.076*B(356)+0.8*B(360)+0.501*B(362)+0.34*B(366)+0.4*B(370)+B(378)+0.207*B(385)

  JVS(437) = B(89)+0.01*B(97)-B(123)+0.35*B(141)+0.35*B(145)+B(147)+0.4*B(149)

  JVS(438) = -B(153)+0.4*B(181)+B(390)+0.048*B(395)+0.021*B(397)+1.6*B(401)+B(403)+0.4*B(411)+B(420)+0.28*B(427)

  JVS(439) = 0

  JVS(440) = B(230)+1.33*B(310)+0.75*B(312)+0.81*B(314)+0.829*B(316)+0.753*B(318)+1.55*B(320)+1.25*B(322)+0.755*B(324)&
               &+1.09*B(326)+B(328)+1.4*B(330)+B(332)+B(334)+B(336)+B(338)+B(340)+2*B(341)+B(343)+0.75*B(346)+0.75*B(348)&
               &+0.96*B(350)+B(391)+B(431)

  JVS(441) = B(183)+0.64*B(185)+0.02*B(187)+0.9*B(189)+0.58*B(191)+0.04*B(195)+0.7*B(201)

  JVS(442) = B(260)+2*B(342)+B(344)+B(379)+B(421)

  JVS(443) = B(231)+0.047*B(235)+0.021*B(237)+1.6*B(241)+B(243)+B(247)+0.25*B(251)+B(261)+0.287*B(267)+0.75*B(445)

  JVS(444) = 0

  JVS(445) = 0.98*B(208)

  JVS(446) = -B(252)-B(290)-B(331)-B(371)-B(412)

  JVS(447) = -B(372)

  JVS(448) = 0

  JVS(449) = -B(413)

  JVS(450) = -B(291)

  JVS(451) = -B(332)

  JVS(452) = 0

  JVS(453) = -B(253)

  JVS(454) = 0

  JVS(455) = 0.98*B(214)

  JVS(456) = -B(254)-B(292)-B(333)-B(373)-B(414)

  JVS(457) = -B(374)

  JVS(458) = 0

  JVS(459) = -B(415)

  JVS(460) = -B(293)

  JVS(461) = -B(334)

  JVS(462) = 0

  JVS(463) = -B(255)

  JVS(464) = 0

  JVS(465) = B(112)

  JVS(466) = -B(248)-B(286)-B(327)-B(367)-B(408)

  JVS(467) = -B(368)

  JVS(468) = B(113)

  JVS(469) = -B(409)

  JVS(470) = -B(287)

  JVS(471) = -B(328)

  JVS(472) = 0

  JVS(473) = -B(249)

  JVS(474) = B(202)

  JVS(475) = 0

  JVS(476) = 0.6*B(180)

  JVS(477) = 0.261*B(238)

  JVS(478) = 0.35*B(250)

  JVS(479) = 0.124*B(236)

  JVS(480) = 0.05*B(252)

  JVS(481) = 0.05*B(254)

  JVS(482) = 0.2*B(248)

  JVS(483) = -B(21)-B(150)

  JVS(484) = B(264)+B(306)+B(347)+B(382)+2*B(386)+1.5*B(387)+B(424)

  JVS(485) = B(308)+0.5*B(349)+0.484*B(384)+1.5*B(388)+B(389)

  JVS(486) = 0.059*B(234)

  JVS(487) = B(383)+0.484*B(385)

  JVS(488) = -B(151)

  JVS(489) = 0.6*B(181)+B(425)

  JVS(490) = B(307)+B(309)

  JVS(491) = B(348)+0.5*B(350)

  JVS(492) = 0

  JVS(493) = 0

  JVS(494) = 0.059*B(235)+0.124*B(237)+0.261*B(239)+0.2*B(249)+0.35*B(251)+0.05*B(253)+0.05*B(255)+B(265)

  JVS(495) = B(203)

  JVS(496) = 0.5*B(252)+B(331)+B(371)+0.5*B(412)

  JVS(497) = 0.95*B(254)+B(333)+B(373)+B(414)

  JVS(498) = -B(20)-B(136)-B(162)-B(198)

  JVS(499) = B(372)+B(374)

  JVS(500) = -B(137)

  JVS(501) = -B(163)+0.5*B(413)+B(415)

  JVS(502) = 0

  JVS(503) = B(332)+B(334)

  JVS(504) = -B(199)

  JVS(505) = 0.5*B(253)+0.95*B(255)

  JVS(506) = 0

  JVS(507) = 0.1*B(116)

  JVS(508) = 0.1*B(118)

  JVS(509) = 0.35*B(144)

  JVS(510) = B(146)

  JVS(511) = 0.05*B(120)

  JVS(512) = 0.13*B(188)

  JVS(513) = B(148)+B(180)

  JVS(514) = 0.416*B(238)+0.322*B(317)+0.318*B(357)+0.587*B(398)

  JVS(515) = 0.334*B(236)+0.245*B(315)+0.237*B(355)+0.391*B(396)

  JVS(516) = 0.5*B(136)+0.5*B(162)

  JVS(517) = -B(428)-B(430)-B(432)-2*B(434)-B(435)-B(437)

  JVS(518) = 0.16*B(262)+0.08*B(345)+0.08*B(380)+0.16*B(422)

  JVS(519) = 0.07*B(142)

  JVS(520) = 0.048*B(234)+0.085*B(313)+0.071*B(353)+0.051*B(394)

  JVS(521) = 0.071*B(354)+0.237*B(356)+0.318*B(358)+0.08*B(381)-B(433)

  JVS(522) = 0.1*B(117)+0.1*B(119)+0.05*B(121)+0.5*B(137)+0.07*B(143)+0.35*B(145)+B(147)+B(149)

  JVS(523) = 0.5*B(163)+B(181)+0.051*B(395)+0.391*B(397)+0.587*B(399)+0.16*B(423)-B(438)

  JVS(524) = -B(429)

  JVS(525) = 0.085*B(314)+0.245*B(316)+0.322*B(318)+0.08*B(346)-B(431)

  JVS(526) = 0.13*B(189)

  JVS(527) = 0

  JVS(528) = 0.048*B(235)+0.334*B(237)+0.416*B(239)+0.16*B(263)-B(436)

  JVS(529) = 0

  JVS(530) = 0.46*B(194)

  JVS(531) = 0.35*B(188)

  JVS(532) = 0.37*B(325)+0.229*B(365)

  JVS(533) = -B(104)-B(168)-B(184)

  JVS(534) = 0.229*B(366)

  JVS(535) = -B(105)

  JVS(536) = -B(169)

  JVS(537) = 0

  JVS(538) = 0.37*B(326)

  JVS(539) = -B(185)+0.35*B(189)+0.46*B(195)

  JVS(540) = 0

  JVS(541) = 0.025*B(100)

  JVS(542) = 0.335*B(96)

  JVS(543) = 0.88*B(138)

  JVS(544) = 0.65*B(192)

  JVS(545) = 0.2*B(240)+0.35*B(319)+0.6*B(359)+0.2*B(400)

  JVS(546) = 0.94*B(242)+0.669*B(321)+0.859*B(361)+0.94*B(402)

  JVS(547) = 1.71*B(244)+0.932*B(323)+0.941*B(363)+1.71*B(404)

  JVS(548) = 0.15*B(238)+0.411*B(317)+0.497*B(357)+0.187*B(398)

  JVS(549) = 0.211*B(236)+0.523*B(315)+0.677*B(355)+0.239*B(396)

  JVS(550) = 0.8*B(248)+B(327)+B(367)+B(408)

  JVS(551) = 0.2*B(21)

  JVS(552) = 0.25*B(162)+0.16*B(198)

  JVS(553) = 0.44*B(184)

  JVS(554) = -B(12)-B(124)-B(154)

  JVS(555) = 0.99*B(186)

  JVS(556) = 0.64*B(387)

  JVS(557) = 1.24*B(266)+0.64*B(349)+0.65*B(384)+0.64*B(388)+1.21*B(389)+1.24*B(426)

  JVS(558) = B(232)+0.75*B(311)+B(351)+B(392)

  JVS(559) = 0.46*B(262)+0.3*B(345)+0.35*B(380)+0.46*B(422)

  JVS(560) = B(14)+0.08*B(142)

  JVS(561) = 0.233*B(234)+0.58*B(313)+0.724*B(353)+0.243*B(394)

  JVS(562) = B(352)+0.724*B(354)+0.677*B(356)+0.497*B(358)+0.6*B(360)+0.859*B(362)+0.941*B(364)+B(368)+0.35*B(381)+0.65&
               &*B(385)

  JVS(563) = 0.335*B(97)+0.025*B(101)-B(125)+0.88*B(139)+0.08*B(143)

  JVS(564) = -B(155)+0.25*B(163)+B(393)+0.243*B(395)+0.239*B(397)+0.187*B(399)+0.2*B(401)+0.94*B(403)+1.71*B(405)+B(409)&
               &+0.46*B(423)+1.24*B(427)

  JVS(565) = 0

  JVS(566) = 0.75*B(312)+0.58*B(314)+0.523*B(316)+0.411*B(318)+0.35*B(320)+0.669*B(322)+0.932*B(324)+B(328)+0.3*B(346)&
               &+0.64*B(350)

  JVS(567) = 0.44*B(185)+0.99*B(187)+0.65*B(193)+0.16*B(199)

  JVS(568) = 0

  JVS(569) = B(233)+0.233*B(235)+0.211*B(237)+0.15*B(239)+0.2*B(241)+0.94*B(243)+1.71*B(245)+0.8*B(249)+0.46*B(263)+1.24&
               &*B(267)

  JVS(570) = 0

  JVS(571) = 0.08*B(325)

  JVS(572) = 0.25*B(250)+0.4*B(329)+0.4*B(369)+0.4*B(410)

  JVS(573) = -B(106)-B(170)-B(186)

  JVS(574) = 0.4*B(370)

  JVS(575) = -B(107)

  JVS(576) = -B(171)+0.4*B(411)

  JVS(577) = 0

  JVS(578) = 0.08*B(326)+0.4*B(330)

  JVS(579) = -B(187)

  JVS(580) = 0.25*B(251)

  JVS(581) = 0.8*B(166)

  JVS(582) = 0.1*B(176)

  JVS(583) = 0.13*B(178)

  JVS(584) = 0.9*B(172)

  JVS(585) = 0.43*B(168)

  JVS(586) = 0.11*B(170)

  JVS(587) = -B(264)-B(306)-B(347)-B(382)-2*B(386)-B(387)-B(424)

  JVS(588) = -B(388)

  JVS(589) = -B(383)

  JVS(590) = 0

  JVS(591) = 0.8*B(167)+0.43*B(169)+0.11*B(171)+0.9*B(173)+0.1*B(177)+0.13*B(179)-B(425)

  JVS(592) = -B(307)

  JVS(593) = -B(348)

  JVS(594) = 0

  JVS(595) = -B(265)

  JVS(596) = 0.2*B(166)

  JVS(597) = 0.9*B(176)

  JVS(598) = 0.87*B(178)

  JVS(599) = 0.1*B(172)

  JVS(600) = 0.57*B(168)

  JVS(601) = 0.89*B(170)

  JVS(602) = -B(387)

  JVS(603) = -B(266)-B(308)-B(349)-B(384)-B(388)-2*B(389)-B(426)

  JVS(604) = -B(385)

  JVS(605) = 0

  JVS(606) = 0.2*B(167)+0.57*B(169)+0.89*B(171)+0.1*B(173)+0.9*B(177)+0.87*B(179)-B(427)

  JVS(607) = -B(309)

  JVS(608) = -B(350)

  JVS(609) = 0

  JVS(610) = -B(267)

  JVS(611) = B(94)

  JVS(612) = 0.2*B(192)

  JVS(613) = 0.16*B(194)

  JVS(614) = B(16)

  JVS(615) = 0.133*B(238)+0.013*B(317)+0.015*B(357)+0.155*B(398)

  JVS(616) = 0.245*B(236)+0.014*B(315)+0.018*B(355)+0.262*B(396)

  JVS(617) = 0

  JVS(618) = 0

  JVS(619) = 0

  JVS(620) = 0.1*B(184)

  JVS(621) = 0.18*B(186)

  JVS(622) = 0

  JVS(623) = 0

  JVS(624) = -B(232)-B(270)-B(311)-B(351)-B(392)

  JVS(625) = 0

  JVS(626) = 0

  JVS(627) = 0.048*B(234)+0.006*B(353)+0.053*B(394)

  JVS(628) = -B(352)+0.006*B(354)+0.018*B(356)+0.015*B(358)

  JVS(629) = B(95)

  JVS(630) = -B(393)+0.053*B(395)+0.262*B(397)+0.155*B(399)

  JVS(631) = -B(271)

  JVS(632) = -B(312)+0.014*B(316)+0.013*B(318)

  JVS(633) = 0.1*B(185)+0.18*B(187)+0.2*B(193)+0.16*B(195)

  JVS(634) = 0

  JVS(635) = -B(233)+0.048*B(235)+0.245*B(237)+0.133*B(239)

  JVS(636) = 0

  JVS(637) = 0.42*B(192)

  JVS(638) = 0.42*B(194)

  JVS(639) = 0.02*B(188)

  JVS(640) = B(126)

  JVS(641) = 0

  JVS(642) = 0

  JVS(643) = 0

  JVS(644) = 0

  JVS(645) = 0

  JVS(646) = 0.03*B(184)

  JVS(647) = 0.12*B(186)

  JVS(648) = 0

  JVS(649) = 0

  JVS(650) = -B(262)-B(304)-B(345)-B(380)-B(422)

  JVS(651) = 0

  JVS(652) = 0

  JVS(653) = -B(381)

  JVS(654) = B(127)

  JVS(655) = -B(423)

  JVS(656) = -B(305)

  JVS(657) = -B(346)

  JVS(658) = 0.03*B(185)+0.12*B(187)+0.02*B(189)+0.42*B(193)+0.42*B(195)

  JVS(659) = 0

  JVS(660) = -B(263)

  JVS(661) = 0

  JVS(662) = B(278)

  JVS(663) = B(280)

  JVS(664) = B(282)

  JVS(665) = B(294)

  JVS(666) = B(276)

  JVS(667) = B(288)

  JVS(668) = B(274)

  JVS(669) = B(290)

  JVS(670) = B(292)

  JVS(671) = B(286)

  JVS(672) = B(428)

  JVS(673) = 0

  JVS(674) = 0

  JVS(675) = B(270)

  JVS(676) = B(304)

  JVS(677) = -B(14)-B(142)

  JVS(678) = B(272)

  JVS(679) = 0

  JVS(680) = -B(143)

  JVS(681) = 0

  JVS(682) = B(271)+B(273)+B(275)+B(277)+B(279)+B(281)+B(283)+B(287)+B(289)+B(291)+B(293)+B(295)+B(300)+B(305)+B(429)

  JVS(683) = 0

  JVS(684) = 0

  JVS(685) = B(301)

  JVS(686) = 0

  JVS(687) = 0

  JVS(688) = 0.583*B(96)

  JVS(689) = B(150)

  JVS(690) = 0

  JVS(691) = 0

  JVS(692) = 0.44*B(142)

  JVS(693) = -B(234)-B(272)-B(313)-B(353)-B(394)

  JVS(694) = -B(354)

  JVS(695) = 0.583*B(97)+0.44*B(143)+B(151)

  JVS(696) = -B(395)

  JVS(697) = -B(273)

  JVS(698) = -B(314)

  JVS(699) = 0

  JVS(700) = 0

  JVS(701) = -B(235)

  JVS(702) = 0

  JVS(703) = B(23)

  JVS(704) = 0.65*B(144)

  JVS(705) = B(226)

  JVS(706) = 0.1*B(190)

  JVS(707) = 0.15*B(188)

  JVS(708) = 0.25*B(444)

  JVS(709) = 0.7*B(200)

  JVS(710) = -B(359)

  JVS(711) = -B(361)

  JVS(712) = -B(363)

  JVS(713) = B(22)+0.1*B(196)

  JVS(714) = B(16)

  JVS(715) = -B(365)

  JVS(716) = B(19)+B(132)+B(158)

  JVS(717) = -B(375)

  JVS(718) = -B(357)

  JVS(719) = -B(369)

  JVS(720) = -B(355)

  JVS(721) = -B(371)

  JVS(722) = -B(373)

  JVS(723) = -B(367)

  JVS(724) = 0

  JVS(725) = 0.28*B(198)

  JVS(726) = -B(432)

  JVS(727) = 0

  JVS(728) = B(124)+B(154)

  JVS(729) = 0

  JVS(730) = -B(382)

  JVS(731) = -B(384)

  JVS(732) = -B(351)

  JVS(733) = 0.23*B(262)+0.12*B(345)-0.88*B(380)+0.23*B(422)

  JVS(734) = 0

  JVS(735) = -B(353)

  JVS(736) = -B(224)-B(258)-B(296)-B(298)-B(337)-B(339)-B(352)-B(354)-B(356)-B(358)-B(360)-B(362)-B(364)-B(366)-B(368)&
               &-B(370)-B(372)-B(374)-B(376)-2*B(377)-0.88*B(381)-B(383)-B(385)-B(418)-B(433)

  JVS(737) = B(125)+B(133)+0.65*B(145)

  JVS(738) = B(155)+B(159)-B(419)+B(420)+0.23*B(423)

  JVS(739) = -B(297)-B(299)

  JVS(740) = -B(338)-B(340)+B(341)+0.12*B(346)

  JVS(741) = 0.15*B(189)+0.1*B(191)+0.1*B(197)+0.28*B(199)+0.7*B(201)

  JVS(742) = B(260)+B(342)+B(421)

  JVS(743) = -B(259)+B(261)+0.23*B(263)+0.25*B(445)

  JVS(744) = -B(225)

  JVS(745) = -B(84)

  JVS(746) = 2*B(30)

  JVS(747) = -B(98)

  JVS(748) = -B(116)

  JVS(749) = -B(118)

  JVS(750) = -B(100)

  JVS(751) = -B(160)

  JVS(752) = -0.964*B(96)

  JVS(753) = -B(92)

  JVS(754) = -B(138)

  JVS(755) = -B(94)

  JVS(756) = -B(88)

  JVS(757) = 0.35*B(6)-B(66)

  JVS(758) = B(13)-0.65*B(140)

  JVS(759) = -B(406)

  JVS(760) = -B(90)

  JVS(761) = B(4)-B(62)

  JVS(762) = 2*B(9)-B(38)

  JVS(763) = -B(128)

  JVS(764) = B(210)

  JVS(765) = B(216)

  JVS(766) = -B(102)+0.12*B(182)

  JVS(767) = B(222)

  JVS(768) = B(5)-B(64)

  JVS(769) = B(15)-B(144)

  JVS(770) = -B(442)

  JVS(771) = -B(146)

  JVS(772) = -B(112)+0.85*B(192)

  JVS(773) = -B(86)

  JVS(774) = -B(114)+0.85*B(194)

  JVS(775) = -B(110)+0.25*B(190)

  JVS(776) = -B(120)

  JVS(777) = -B(108)+0.28*B(188)

  JVS(778) = 0

  JVS(779) = -B(148)+0.036*B(200)

  JVS(780) = -B(130)

  JVS(781) = -B(134)+0.19*B(196)

  JVS(782) = -B(126)

  JVS(783) = 0

  JVS(784) = -B(132)

  JVS(785) = 0

  JVS(786) = 0

  JVS(787) = 0

  JVS(788) = 0

  JVS(789) = -B(122)

  JVS(790) = 0

  JVS(791) = 0

  JVS(792) = 0

  JVS(793) = -B(150)

  JVS(794) = -B(136)+0.21*B(198)

  JVS(795) = 0

  JVS(796) = -B(104)+0.4*B(184)

  JVS(797) = -B(124)

  JVS(798) = -B(106)+0.63*B(186)

  JVS(799) = 0

  JVS(800) = 0

  JVS(801) = 0

  JVS(802) = 0

  JVS(803) = B(14)-0.51*B(142)

  JVS(804) = 0

  JVS(805) = 0

  JVS(806) = -B(32)-B(36)-B(39)-B(49)-B(51)-B(53)-B(63)-B(65)-B(67)-B(82)-B(85)-B(87)-B(89)-B(91)-B(93)-B(95)-0.964&
               &*B(97)-B(99)-B(101)-B(103)-B(105)-B(107)-B(109)-B(111)-B(113)-B(115)-B(117)-B(119)-B(121)-B(123)-B(125)&
               &-B(127)-B(129)-B(131)-B(133)-B(135)-B(137)-B(139)-0.65*B(141)-0.51*B(143)-B(145)-B(147)-B(149)-B(151)-B(161)&
               &-B(407)-B(443)

  JVS(807) = -B(54)+0.7*B(60)

  JVS(808) = B(34)-B(37)+B(55)+0.7*B(61)

  JVS(809) = 0

  JVS(810) = -B(33)+B(35)+0.12*B(183)+0.4*B(185)+0.63*B(187)+0.28*B(189)+0.25*B(191)+0.85*B(193)+0.85*B(195)+0.19*B(197)&
               &+0.21*B(199)+0.036*B(201)+B(211)+B(217)+B(223)

  JVS(811) = 0

  JVS(812) = -B(50)+B(56)

  JVS(813) = -B(52)

  JVS(814) = B(80)

  JVS(815) = 0.35*B(6)

  JVS(816) = B(47)

  JVS(817) = -B(166)

  JVS(818) = B(64)

  JVS(819) = B(146)

  JVS(820) = -B(176)

  JVS(821) = -B(178)

  JVS(822) = -B(174)

  JVS(823) = -B(164)

  JVS(824) = -B(172)

  JVS(825) = 0.6*B(148)-0.4*B(180)

  JVS(826) = -B(156)

  JVS(827) = -B(400)

  JVS(828) = -B(402)

  JVS(829) = -B(404)

  JVS(830) = -B(158)

  JVS(831) = -B(416)

  JVS(832) = -B(398)

  JVS(833) = -B(410)

  JVS(834) = -B(396)

  JVS(835) = -B(152)

  JVS(836) = -B(412)

  JVS(837) = -B(414)

  JVS(838) = -B(408)

  JVS(839) = -B(162)

  JVS(840) = -B(437)

  JVS(841) = -B(168)

  JVS(842) = -B(154)

  JVS(843) = -B(170)

  JVS(844) = -B(424)

  JVS(845) = -B(426)

  JVS(846) = -B(392)

  JVS(847) = -B(422)

  JVS(848) = 0

  JVS(849) = -B(394)

  JVS(850) = -B(418)

  JVS(851) = -B(53)+B(65)+B(147)+0.6*B(149)

  JVS(852) = -B(7)-B(8)-B(54)-B(60)-B(74)-B(76)-B(78)-2*B(81)-B(153)-B(155)-B(157)-B(159)-B(163)-B(165)-B(167)-B(169)&
               &-B(171)-B(173)-B(175)-B(177)-B(179)-0.4*B(181)-B(390)-B(393)-B(395)-B(397)-B(399)-B(401)-B(403)-B(405)&
               &-B(409)-B(411)-B(413)-B(415)-B(417)-B(419)-B(420)-B(423)-B(425)-B(427)-B(438)

  JVS(853) = -B(61)

  JVS(854) = -B(391)

  JVS(855) = B(70)

  JVS(856) = -B(421)

  JVS(857) = -B(75)

  JVS(858) = B(48)+B(71)-B(77)-B(79)

  JVS(859) = B(84)

  JVS(860) = 0.25*B(98)

  JVS(861) = 0.1*B(116)

  JVS(862) = 0.1*B(118)

  JVS(863) = 0.049*B(100)

  JVS(864) = 0.381*B(96)

  JVS(865) = B(138)

  JVS(866) = 0.65*B(6)+B(59)

  JVS(867) = B(13)

  JVS(868) = B(90)

  JVS(869) = B(38)

  JVS(870) = B(23)+B(128)

  JVS(871) = -B(204)

  JVS(872) = 0.02*B(208)

  JVS(873) = 0.02*B(214)

  JVS(874) = 0.26*B(182)

  JVS(875) = 0.02*B(220)

  JVS(876) = 0.35*B(144)

  JVS(877) = 0

  JVS(878) = 0.1*B(192)

  JVS(879) = B(86)

  JVS(880) = 0.1*B(194)

  JVS(881) = 0.25*B(190)

  JVS(882) = 0.05*B(120)

  JVS(883) = 0.3*B(188)

  JVS(884) = 0.75*B(444)-B(446)+B(448)

  JVS(885) = 0.4*B(148)+0.08*B(200)

  JVS(886) = 0.8*B(18)+B(130)+B(156)

  JVS(887) = B(240)-B(278)+B(319)+0.5*B(359)+B(400)

  JVS(888) = B(242)-B(280)+B(321)+0.501*B(361)+B(402)

  JVS(889) = B(244)-B(282)+B(323)+0.51*B(363)+B(404)

  JVS(890) = B(22)+0.32*B(196)

  JVS(891) = B(246)-B(284)+B(325)+0.506*B(365)+B(439)

  JVS(892) = B(19)

  JVS(893) = B(256)-B(294)+2*B(335)+B(375)+B(416)

  JVS(894) = 0.606*B(238)-B(276)+0.993*B(317)+0.489*B(357)+0.845*B(398)

  JVS(895) = 0.65*B(250)-B(288)+2*B(329)+B(369)+B(410)

  JVS(896) = 0.599*B(236)-B(274)+0.946*B(315)+0.438*B(355)+0.699*B(396)

  JVS(897) = 2*B(11)+B(122)+B(152)

  JVS(898) = 0.95*B(252)-B(290)+B(331)+B(371)+B(412)

  JVS(899) = 0.95*B(254)-B(292)+B(333)+B(373)+B(414)

  JVS(900) = 0.8*B(248)-B(286)+2*B(327)+B(367)+B(408)

  JVS(901) = B(21)

  JVS(902) = B(20)+0.5*B(136)+0.5*B(162)+0.29*B(198)

  JVS(903) = -B(428)+B(430)

  JVS(904) = 0.25*B(184)

  JVS(905) = B(12)

  JVS(906) = 0.22*B(186)

  JVS(907) = B(264)-B(306)+B(347)+0.5*B(382)+B(386)+0.5*B(387)+B(424)

  JVS(908) = -B(308)+0.5*B(349)+0.5*B(388)

  JVS(909) = B(232)-B(270)+B(311)+0.5*B(351)+B(392)

  JVS(910) = 0.77*B(262)-B(304)+0.88*B(345)+0.38*B(380)+0.77*B(422)

  JVS(911) = B(14)

  JVS(912) = 0.742*B(234)-B(272)+0.992*B(313)+0.488*B(353)+0.792*B(394)

  JVS(913) = -B(296)-B(298)+B(337)+0.5*B(352)+0.488*B(354)+0.438*B(356)+0.489*B(358)+0.5*B(360)+0.501*B(362)+0.51*B(364)&
               &+0.506*B(366)+B(368)+B(370)+B(372)+B(374)+B(376)+0.38*B(381)+0.5*B(383)

  JVS(914) = B(32)-B(36)+B(39)+B(53)+B(82)+B(85)+B(87)+B(91)+0.381*B(97)+0.25*B(99)+0.049*B(101)+0.1*B(117)+0.1*B(119)&
               &+0.05*B(121)+B(123)+B(129)+B(131)+0.5*B(137)+B(139)+0.35*B(145)+0.4*B(149)

  JVS(915) = B(54)-B(60)+B(153)+B(157)+0.5*B(163)+B(390)+B(393)+0.792*B(395)+0.699*B(397)+0.845*B(399)+B(401)+B(403)&
               &+B(405)+B(409)+B(411)+B(413)+B(415)+B(417)+0.77*B(423)+B(425)

  JVS(916) = -B(34)-B(37)-2*B(40)-2*B(41)-B(55)-B(57)-B(61)-B(205)-B(268)-B(271)-B(273)-B(275)-B(277)-B(279)-B(281)&
               &-B(283)-B(285)-B(287)-B(289)-B(291)-B(293)-B(295)-B(297)-B(299)-B(300)-B(302)-B(305)-B(307)-B(309)-B(429)&
               &-B(447)

  JVS(917) = B(230)-B(269)+0.66*B(310)+B(312)+0.992*B(314)+0.946*B(316)+0.993*B(318)+B(320)+B(322)+B(324)+B(326)+2&
               &*B(328)+2*B(330)+B(332)+B(334)+2*B(336)+B(338)+B(341)+0.88*B(346)+B(348)+0.5*B(350)+B(391)+B(431)

  JVS(918) = B(33)-B(35)+0.26*B(183)+0.25*B(185)+0.22*B(187)+0.3*B(189)+0.25*B(191)+0.1*B(193)+0.1*B(195)+0.32*B(197)&
               &+0.29*B(199)+0.08*B(201)

  JVS(919) = -B(301)-B(303)+B(342)

  JVS(920) = -B(56)+B(231)+B(233)+0.742*B(235)+0.599*B(237)+0.606*B(239)+B(241)+B(243)+B(245)+B(247)+0.8*B(249)+0.65&
               &*B(251)+0.95*B(253)+0.95*B(255)+B(257)+0.77*B(263)+B(265)+0.75*B(445)

  JVS(921) = -B(58)

  JVS(922) = B(92)

  JVS(923) = 0.65*B(140)

  JVS(924) = B(15)

  JVS(925) = 0.08*B(190)

  JVS(926) = 0.03*B(188)

  JVS(927) = -B(319)+0.5*B(359)

  JVS(928) = -B(321)+0.501*B(361)

  JVS(929) = -B(323)+0.51*B(363)

  JVS(930) = -B(325)+0.506*B(365)

  JVS(931) = -B(335)+B(375)

  JVS(932) = -B(317)+0.507*B(357)

  JVS(933) = -B(329)+B(369)

  JVS(934) = 0.031*B(236)-0.951*B(315)+0.554*B(355)+0.04*B(396)

  JVS(935) = -B(331)+B(371)

  JVS(936) = -B(333)+B(373)

  JVS(937) = -B(327)+B(367)

  JVS(938) = 0

  JVS(939) = -B(430)+B(432)

  JVS(940) = 0.19*B(184)

  JVS(941) = B(12)

  JVS(942) = 0.23*B(186)

  JVS(943) = -B(347)+0.5*B(382)

  JVS(944) = -B(349)+0.516*B(384)

  JVS(945) = -B(311)+0.5*B(351)

  JVS(946) = -B(345)+0.5*B(380)

  JVS(947) = 0

  JVS(948) = 0.15*B(234)-0.993*B(313)+0.508*B(353)+0.155*B(394)

  JVS(949) = B(258)-B(339)+0.5*B(352)+0.508*B(354)+0.554*B(356)+0.507*B(358)+0.5*B(360)+0.501*B(362)+0.51*B(364)+0.506&
               &*B(366)+B(368)+B(370)+B(372)+B(374)+B(376)+2*B(377)+B(378)+0.5*B(381)+0.5*B(383)+0.516*B(385)+B(418)+B(433)

  JVS(950) = B(93)+0.65*B(141)

  JVS(951) = -B(390)+0.155*B(395)+0.04*B(397)+B(419)

  JVS(952) = -B(268)

  JVS(953) = -B(230)-B(269)-2*B(310)-B(312)-0.993*B(314)-0.951*B(316)-B(318)-B(320)-B(322)-B(324)-B(326)-B(328)-B(330)&
               &-B(332)-B(334)-B(336)-B(340)-B(341)-B(343)-B(346)-B(348)-B(350)-B(391)-B(431)

  JVS(954) = 0.19*B(185)+0.23*B(187)+0.03*B(189)+0.08*B(191)

  JVS(955) = -B(342)-B(344)+B(379)

  JVS(956) = -B(231)+0.15*B(235)+0.031*B(237)+B(259)

  JVS(957) = 0

  JVS(958) = B(24)-B(26)

  JVS(959) = -B(210)

  JVS(960) = -B(216)

  JVS(961) = -B(182)

  JVS(962) = -B(222)

  JVS(963) = -B(192)

  JVS(964) = -B(194)

  JVS(965) = -B(190)

  JVS(966) = 0

  JVS(967) = -B(188)

  JVS(968) = -B(200)

  JVS(969) = -B(196)

  JVS(970) = 0

  JVS(971) = 0

  JVS(972) = -B(198)

  JVS(973) = -B(184)

  JVS(974) = -B(186)

  JVS(975) = B(298)

  JVS(976) = -B(32)

  JVS(977) = 0

  JVS(978) = -B(34)+B(299)+B(302)

  JVS(979) = 0

  JVS(980) = -B(2)-B(3)-B(27)-B(33)-B(35)-B(68)-B(70)-B(183)-B(185)-B(187)-B(189)-B(191)-B(193)-B(195)-B(197)-B(199)&
               &-B(201)-B(211)-B(217)-B(223)

  JVS(981) = B(303)

  JVS(982) = -B(69)

  JVS(983) = -B(71)

  JVS(984) = B(229)

  JVS(985) = B(20)+0.5*B(136)+0.5*B(162)

  JVS(986) = -B(378)

  JVS(987) = 0.5*B(137)

  JVS(988) = 0.5*B(163)-B(420)

  JVS(989) = -B(300)-B(302)

  JVS(990) = -B(341)-B(343)

  JVS(991) = 0

  JVS(992) = -B(227)-B(260)-B(301)-B(303)-B(342)-B(344)-B(379)-B(421)

  JVS(993) = -B(261)

  JVS(994) = -B(228)

  JVS(995) = B(4)

  JVS(996) = -B(43)+B(45)

  JVS(997) = 0

  JVS(998) = 0

  JVS(999) = 0

  JVS(1000) = 0

  JVS(1001) = 0

  JVS(1002) = -B(444)

  JVS(1003) = -B(240)

  JVS(1004) = -B(242)

  JVS(1005) = -B(244)

  JVS(1006) = 0

  JVS(1007) = -B(246)

  JVS(1008) = -B(256)

  JVS(1009) = -B(238)

  JVS(1010) = -B(250)

  JVS(1011) = -B(236)

  JVS(1012) = -B(252)

  JVS(1013) = -B(254)

  JVS(1014) = -B(248)

  JVS(1015) = -B(435)

  JVS(1016) = 0

  JVS(1017) = 0

  JVS(1018) = -B(264)

  JVS(1019) = -B(266)

  JVS(1020) = -B(232)

  JVS(1021) = -B(262)

  JVS(1022) = 0

  JVS(1023) = -B(234)

  JVS(1024) = -B(258)

  JVS(1025) = -B(49)

  JVS(1026) = B(7)-B(74)+B(76)

  JVS(1027) = -B(55)

  JVS(1028) = -B(230)

  JVS(1029) = -B(68)

  JVS(1030) = -B(260)

  JVS(1031) = -B(44)-B(50)-B(56)-B(69)-2*B(72)-B(75)-B(231)-B(233)-B(235)-B(237)-B(239)-B(241)-B(243)-B(245)-B(247)&
                &-B(249)-B(251)-B(253)-B(255)-B(257)-B(259)-B(261)-B(263)-B(265)-B(267)-B(436)-B(445)

  JVS(1032) = B(1)+B(46)+B(77)

  JVS(1033) = B(80)

  JVS(1034) = B(88)

  JVS(1035) = 0.65*B(6)+B(59)+B(66)

  JVS(1036) = B(406)+B(451)

  JVS(1037) = B(62)

  JVS(1038) = B(43)-B(45)-B(47)

  JVS(1039) = -B(202)

  JVS(1040) = -B(206)

  JVS(1041) = -B(212)

  JVS(1042) = -B(218)

  JVS(1043) = B(5)

  JVS(1044) = 0

  JVS(1045) = B(226)

  JVS(1046) = 0

  JVS(1047) = 0

  JVS(1048) = 0

  JVS(1049) = B(444)-B(449)

  JVS(1050) = 0.4*B(180)+0.7*B(200)+B(229)

  JVS(1051) = 0

  JVS(1052) = B(240)+B(400)

  JVS(1053) = B(242)+B(402)

  JVS(1054) = B(244)+B(404)

  JVS(1055) = 0

  JVS(1056) = B(246)

  JVS(1057) = 0

  JVS(1058) = B(256)+B(416)

  JVS(1059) = 0.739*B(238)+B(398)

  JVS(1060) = 0.65*B(250)+B(410)

  JVS(1061) = 0.876*B(236)+B(396)

  JVS(1062) = 0

  JVS(1063) = 0.95*B(252)+B(412)

  JVS(1064) = 0.95*B(254)+B(414)

  JVS(1065) = 0.8*B(248)+B(408)

  JVS(1066) = B(21)+B(150)

  JVS(1067) = 0.5*B(162)

  JVS(1068) = B(435)+B(437)

  JVS(1069) = 0

  JVS(1070) = 0

  JVS(1071) = 0

  JVS(1072) = B(264)+0.5*B(387)+B(424)

  JVS(1073) = 2*B(266)+0.5*B(349)+0.516*B(384)+0.5*B(388)+B(389)+2*B(426)

  JVS(1074) = B(232)+B(392)

  JVS(1075) = B(262)+B(422)

  JVS(1076) = 0

  JVS(1077) = 0.941*B(234)+B(394)

  JVS(1078) = -B(224)+B(258)+0.516*B(385)+B(418)

  JVS(1079) = -B(51)+B(53)+B(63)+B(67)+B(89)+B(151)+B(407)

  JVS(1080) = B(8)+B(54)+0.7*B(60)+2*B(74)-B(78)+2*B(81)+0.5*B(163)+0.4*B(181)+B(390)+B(393)+B(395)+B(397)+B(399)+B(401)&
                &+B(403)+B(405)+B(409)+B(411)+B(413)+B(415)+B(417)+B(419)+B(420)+B(423)+B(425)+2*B(427)+B(438)

  JVS(1081) = B(55)-B(57)+0.7*B(61)

  JVS(1082) = B(230)+0.5*B(350)+B(391)

  JVS(1083) = B(68)-B(70)+0.7*B(201)

  JVS(1084) = -B(227)+B(260)+B(421)

  JVS(1085) = B(44)+B(56)+B(69)+2*B(72)+2*B(75)+B(231)+B(233)+0.941*B(235)+0.876*B(237)+0.739*B(239)+B(241)+B(243)&
                &+B(245)+B(247)+0.8*B(249)+0.65*B(251)+0.95*B(253)+0.95*B(255)+B(257)+B(259)+B(261)+B(263)+B(265)+2*B(267)&
                &+B(436)+B(445)

  JVS(1086) = -B(1)-B(46)-B(48)-B(52)-B(58)-B(71)-B(79)-B(203)-B(207)-B(213)-B(219)-B(225)-B(228)-B(450)
      
END SUBROUTINE racm_esrlsorg_Jac_SP














SUBROUTINE racm_esrlsorg_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(1086), W(80), a
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
      
END SUBROUTINE racm_esrlsorg_KppDecomp



SUBROUTINE racm_esrlsorg_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(1086), W(80), a
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
      
END SUBROUTINE racm_esrlsorg_KppDecompCmplx


SUBROUTINE racm_esrlsorg_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(1086), X(80), sum

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
      
END SUBROUTINE racm_esrlsorg_KppSolveIndirect


SUBROUTINE racm_esrlsorg_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(1086), X(80), sum

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
      
END SUBROUTINE racm_esrlsorg_KppSolveCmplx













SUBROUTINE racm_esrlsorg_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(22) = X(22)-JVS(93)*X(21)
  X(25) = X(25)-JVS(117)*X(12)
  X(26) = X(26)-JVS(123)*X(6)
  X(28) = X(28)-JVS(136)*X(9)
  X(29) = X(29)-JVS(141)*X(10)
  X(37) = X(37)-JVS(189)*X(14)-JVS(190)*X(18)-JVS(191)*X(30)-JVS(192)*X(34)-JVS(193)*X(36)
  X(40) = X(40)-JVS(221)*X(27)-JVS(222)*X(28)-JVS(223)*X(29)-JVS(224)*X(31)
  X(42) = X(42)-JVS(235)*X(13)-JVS(236)*X(21)-JVS(237)*X(39)
  X(44) = X(44)-JVS(252)*X(14)
  X(45) = X(45)-JVS(265)*X(30)
  X(48) = X(48)-JVS(290)*X(7)-JVS(291)*X(38)-JVS(292)*X(39)-JVS(293)*X(41)
  X(49) = X(49)-JVS(304)*X(8)-JVS(305)*X(16)-JVS(306)*X(36)-JVS(307)*X(46)-JVS(308)*X(47)
  X(50) = X(50)-JVS(329)*X(39)-JVS(330)*X(41)
  X(51) = X(51)-JVS(339)*X(22)-JVS(340)*X(25)-JVS(341)*X(34)-JVS(342)*X(39)-JVS(343)*X(42)-JVS(344)*X(43)-JVS(345)*X(48)&
            &-JVS(346)*X(50)
  X(52) = X(52)-JVS(364)*X(31)-JVS(365)*X(40)
  X(53) = X(53)-JVS(375)*X(12)
  X(54) = X(54)-JVS(383)*X(38)
  X(55) = X(55)-JVS(392)*X(8)
  X(56) = X(56)-JVS(400)*X(14)-JVS(401)*X(18)-JVS(402)*X(20)-JVS(403)*X(25)-JVS(404)*X(30)-JVS(405)*X(33)-JVS(406)*X(34)&
            &-JVS(407)*X(35)-JVS(408)*X(38)-JVS(409)*X(39)-JVS(410)*X(41)-JVS(411)*X(42)-JVS(412)*X(43)-JVS(413)*X(44)&
            &-JVS(414)*X(45)-JVS(415)*X(46)-JVS(416)*X(47)-JVS(417)*X(48)-JVS(418)*X(50)-JVS(419)*X(52)-JVS(420)*X(53)&
            &-JVS(421)*X(54)-JVS(422)*X(55)
  X(57) = X(57)-JVS(445)*X(28)
  X(58) = X(58)-JVS(455)*X(29)
  X(59) = X(59)-JVS(465)*X(36)
  X(60) = X(60)-JVS(474)*X(27)-JVS(475)*X(40)-JVS(476)*X(43)-JVS(477)*X(53)-JVS(478)*X(54)-JVS(479)*X(55)-JVS(480)*X(57)&
            &-JVS(481)*X(58)-JVS(482)*X(59)
  X(61) = X(61)-JVS(496)*X(57)-JVS(497)*X(58)
  X(62) = X(62)-JVS(507)*X(9)-JVS(508)*X(10)-JVS(509)*X(33)-JVS(510)*X(35)-JVS(511)*X(40)-JVS(512)*X(41)-JVS(513)*X(43)&
            &-JVS(514)*X(53)-JVS(515)*X(55)-JVS(516)*X(61)
  X(63) = X(63)-JVS(530)*X(38)-JVS(531)*X(41)-JVS(532)*X(50)
  X(64) = X(64)-JVS(541)*X(12)-JVS(542)*X(14)-JVS(543)*X(16)-JVS(544)*X(36)-JVS(545)*X(45)-JVS(546)*X(46)-JVS(547)*X(47)&
            &-JVS(548)*X(53)-JVS(549)*X(55)-JVS(550)*X(59)-JVS(551)*X(60)-JVS(552)*X(61)-JVS(553)*X(63)
  X(65) = X(65)-JVS(571)*X(50)-JVS(572)*X(54)
  X(66) = X(66)-JVS(581)*X(30)-JVS(582)*X(36)-JVS(583)*X(38)-JVS(584)*X(41)-JVS(585)*X(63)-JVS(586)*X(65)
  X(67) = X(67)-JVS(596)*X(30)-JVS(597)*X(36)-JVS(598)*X(38)-JVS(599)*X(41)-JVS(600)*X(63)-JVS(601)*X(65)-JVS(602)*X(66)
  X(68) = X(68)-JVS(611)*X(17)-JVS(612)*X(36)-JVS(613)*X(38)-JVS(614)*X(49)-JVS(615)*X(53)-JVS(616)*X(55)-JVS(617)*X(59)&
            &-JVS(618)*X(60)-JVS(619)*X(61)-JVS(620)*X(63)-JVS(621)*X(65)-JVS(622)*X(66)-JVS(623)*X(67)
  X(69) = X(69)-JVS(637)*X(36)-JVS(638)*X(38)-JVS(639)*X(41)-JVS(640)*X(49)-JVS(641)*X(53)-JVS(642)*X(55)-JVS(643)*X(59)&
            &-JVS(644)*X(60)-JVS(645)*X(61)-JVS(646)*X(63)-JVS(647)*X(65)-JVS(648)*X(66)-JVS(649)*X(67)
  X(70) = X(70)-JVS(662)*X(45)-JVS(663)*X(46)-JVS(664)*X(47)-JVS(665)*X(52)-JVS(666)*X(53)-JVS(667)*X(54)-JVS(668)*X(55)&
            &-JVS(669)*X(57)-JVS(670)*X(58)-JVS(671)*X(59)-JVS(672)*X(62)-JVS(673)*X(63)-JVS(674)*X(65)-JVS(675)*X(68)&
            &-JVS(676)*X(69)
  X(71) = X(71)-JVS(688)*X(14)-JVS(689)*X(60)-JVS(690)*X(66)-JVS(691)*X(67)-JVS(692)*X(70)
  X(72) = X(72)-JVS(703)*X(25)-JVS(704)*X(33)-JVS(705)*X(35)-JVS(706)*X(39)-JVS(707)*X(41)-JVS(708)*X(42)-JVS(709)*X(43)&
            &-JVS(710)*X(45)-JVS(711)*X(46)-JVS(712)*X(47)-JVS(713)*X(48)-JVS(714)*X(49)-JVS(715)*X(50)-JVS(716)*X(51)&
            &-JVS(717)*X(52)-JVS(718)*X(53)-JVS(719)*X(54)-JVS(720)*X(55)-JVS(721)*X(57)-JVS(722)*X(58)-JVS(723)*X(59)&
            &-JVS(724)*X(60)-JVS(725)*X(61)-JVS(726)*X(62)-JVS(727)*X(63)-JVS(728)*X(64)-JVS(729)*X(65)-JVS(730)*X(66)&
            &-JVS(731)*X(67)-JVS(732)*X(68)-JVS(733)*X(69)-JVS(734)*X(70)-JVS(735)*X(71)
  X(73) = X(73)-JVS(745)*X(5)-JVS(746)*X(6)-JVS(747)*X(8)-JVS(748)*X(9)-JVS(749)*X(10)-JVS(750)*X(12)-JVS(751)*X(13)&
            &-JVS(752)*X(14)-JVS(753)*X(15)-JVS(754)*X(16)-JVS(755)*X(17)-JVS(756)*X(18)-JVS(757)*X(19)-JVS(758)*X(20)&
            &-JVS(759)*X(21)-JVS(760)*X(22)-JVS(761)*X(23)-JVS(762)*X(24)-JVS(763)*X(25)-JVS(764)*X(28)-JVS(765)*X(29)&
            &-JVS(766)*X(30)-JVS(767)*X(31)-JVS(768)*X(32)-JVS(769)*X(33)-JVS(770)*X(34)-JVS(771)*X(35)-JVS(772)*X(36)&
            &-JVS(773)*X(37)-JVS(774)*X(38)-JVS(775)*X(39)-JVS(776)*X(40)-JVS(777)*X(41)-JVS(778)*X(42)-JVS(779)*X(43)&
            &-JVS(780)*X(44)-JVS(781)*X(48)-JVS(782)*X(49)-JVS(783)*X(50)-JVS(784)*X(51)-JVS(785)*X(52)-JVS(786)*X(53)&
            &-JVS(787)*X(54)-JVS(788)*X(55)-JVS(789)*X(56)-JVS(790)*X(57)-JVS(791)*X(58)-JVS(792)*X(59)-JVS(793)*X(60)&
            &-JVS(794)*X(61)-JVS(795)*X(62)-JVS(796)*X(63)-JVS(797)*X(64)-JVS(798)*X(65)-JVS(799)*X(66)-JVS(800)*X(67)&
            &-JVS(801)*X(68)-JVS(802)*X(69)-JVS(803)*X(70)-JVS(804)*X(71)-JVS(805)*X(72)
  X(74) = X(74)-JVS(814)*X(11)-JVS(815)*X(19)-JVS(816)*X(26)-JVS(817)*X(30)-JVS(818)*X(32)-JVS(819)*X(35)-JVS(820)*X(36)&
            &-JVS(821)*X(38)-JVS(822)*X(39)-JVS(823)*X(40)-JVS(824)*X(41)-JVS(825)*X(43)-JVS(826)*X(44)-JVS(827)*X(45)&
            &-JVS(828)*X(46)-JVS(829)*X(47)-JVS(830)*X(51)-JVS(831)*X(52)-JVS(832)*X(53)-JVS(833)*X(54)-JVS(834)*X(55)&
            &-JVS(835)*X(56)-JVS(836)*X(57)-JVS(837)*X(58)-JVS(838)*X(59)-JVS(839)*X(61)-JVS(840)*X(62)-JVS(841)*X(63)&
            &-JVS(842)*X(64)-JVS(843)*X(65)-JVS(844)*X(66)-JVS(845)*X(67)-JVS(846)*X(68)-JVS(847)*X(69)-JVS(848)*X(70)&
            &-JVS(849)*X(71)-JVS(850)*X(72)-JVS(851)*X(73)
  X(75) = X(75)-JVS(859)*X(5)-JVS(860)*X(8)-JVS(861)*X(9)-JVS(862)*X(10)-JVS(863)*X(12)-JVS(864)*X(14)-JVS(865)*X(16)&
            &-JVS(866)*X(19)-JVS(867)*X(20)-JVS(868)*X(22)-JVS(869)*X(24)-JVS(870)*X(25)-JVS(871)*X(27)-JVS(872)*X(28)&
            &-JVS(873)*X(29)-JVS(874)*X(30)-JVS(875)*X(31)-JVS(876)*X(33)-JVS(877)*X(34)-JVS(878)*X(36)-JVS(879)*X(37)&
            &-JVS(880)*X(38)-JVS(881)*X(39)-JVS(882)*X(40)-JVS(883)*X(41)-JVS(884)*X(42)-JVS(885)*X(43)-JVS(886)*X(44)&
            &-JVS(887)*X(45)-JVS(888)*X(46)-JVS(889)*X(47)-JVS(890)*X(48)-JVS(891)*X(50)-JVS(892)*X(51)-JVS(893)*X(52)&
            &-JVS(894)*X(53)-JVS(895)*X(54)-JVS(896)*X(55)-JVS(897)*X(56)-JVS(898)*X(57)-JVS(899)*X(58)-JVS(900)*X(59)&
            &-JVS(901)*X(60)-JVS(902)*X(61)-JVS(903)*X(62)-JVS(904)*X(63)-JVS(905)*X(64)-JVS(906)*X(65)-JVS(907)*X(66)&
            &-JVS(908)*X(67)-JVS(909)*X(68)-JVS(910)*X(69)-JVS(911)*X(70)-JVS(912)*X(71)-JVS(913)*X(72)-JVS(914)*X(73)&
            &-JVS(915)*X(74)
  X(76) = X(76)-JVS(922)*X(15)-JVS(923)*X(20)-JVS(924)*X(33)-JVS(925)*X(39)-JVS(926)*X(41)-JVS(927)*X(45)-JVS(928)*X(46)&
            &-JVS(929)*X(47)-JVS(930)*X(50)-JVS(931)*X(52)-JVS(932)*X(53)-JVS(933)*X(54)-JVS(934)*X(55)-JVS(935)*X(57)&
            &-JVS(936)*X(58)-JVS(937)*X(59)-JVS(938)*X(61)-JVS(939)*X(62)-JVS(940)*X(63)-JVS(941)*X(64)-JVS(942)*X(65)&
            &-JVS(943)*X(66)-JVS(944)*X(67)-JVS(945)*X(68)-JVS(946)*X(69)-JVS(947)*X(70)-JVS(948)*X(71)-JVS(949)*X(72)&
            &-JVS(950)*X(73)-JVS(951)*X(74)-JVS(952)*X(75)
  X(77) = X(77)-JVS(958)*X(26)-JVS(959)*X(28)-JVS(960)*X(29)-JVS(961)*X(30)-JVS(962)*X(31)-JVS(963)*X(36)-JVS(964)*X(38)&
            &-JVS(965)*X(39)-JVS(966)*X(40)-JVS(967)*X(41)-JVS(968)*X(43)-JVS(969)*X(48)-JVS(970)*X(50)-JVS(971)*X(54)&
            &-JVS(972)*X(61)-JVS(973)*X(63)-JVS(974)*X(65)-JVS(975)*X(72)-JVS(976)*X(73)-JVS(977)*X(74)-JVS(978)*X(75)&
            &-JVS(979)*X(76)
  X(78) = X(78)-JVS(984)*X(43)-JVS(985)*X(61)-JVS(986)*X(72)-JVS(987)*X(73)-JVS(988)*X(74)-JVS(989)*X(75)-JVS(990)*X(76)&
            &-JVS(991)*X(77)
  X(79) = X(79)-JVS(995)*X(23)-JVS(996)*X(26)-JVS(997)*X(28)-JVS(998)*X(29)-JVS(999)*X(31)-JVS(1000)*X(40)-JVS(1001)&
            &*X(41)-JVS(1002)*X(42)-JVS(1003)*X(45)-JVS(1004)*X(46)-JVS(1005)*X(47)-JVS(1006)*X(48)-JVS(1007)*X(50)&
            &-JVS(1008)*X(52)-JVS(1009)*X(53)-JVS(1010)*X(54)-JVS(1011)*X(55)-JVS(1012)*X(57)-JVS(1013)*X(58)-JVS(1014)&
            &*X(59)-JVS(1015)*X(62)-JVS(1016)*X(63)-JVS(1017)*X(65)-JVS(1018)*X(66)-JVS(1019)*X(67)-JVS(1020)*X(68)&
            &-JVS(1021)*X(69)-JVS(1022)*X(70)-JVS(1023)*X(71)-JVS(1024)*X(72)-JVS(1025)*X(73)-JVS(1026)*X(74)-JVS(1027)&
            &*X(75)-JVS(1028)*X(76)-JVS(1029)*X(77)-JVS(1030)*X(78)
  X(80) = X(80)-JVS(1033)*X(11)-JVS(1034)*X(18)-JVS(1035)*X(19)-JVS(1036)*X(21)-JVS(1037)*X(23)-JVS(1038)*X(26)&
            &-JVS(1039)*X(27)-JVS(1040)*X(28)-JVS(1041)*X(29)-JVS(1042)*X(31)-JVS(1043)*X(32)-JVS(1044)*X(34)-JVS(1045)&
            &*X(35)-JVS(1046)*X(39)-JVS(1047)*X(40)-JVS(1048)*X(41)-JVS(1049)*X(42)-JVS(1050)*X(43)-JVS(1051)*X(44)&
            &-JVS(1052)*X(45)-JVS(1053)*X(46)-JVS(1054)*X(47)-JVS(1055)*X(48)-JVS(1056)*X(50)-JVS(1057)*X(51)-JVS(1058)&
            &*X(52)-JVS(1059)*X(53)-JVS(1060)*X(54)-JVS(1061)*X(55)-JVS(1062)*X(56)-JVS(1063)*X(57)-JVS(1064)*X(58)&
            &-JVS(1065)*X(59)-JVS(1066)*X(60)-JVS(1067)*X(61)-JVS(1068)*X(62)-JVS(1069)*X(63)-JVS(1070)*X(64)-JVS(1071)&
            &*X(65)-JVS(1072)*X(66)-JVS(1073)*X(67)-JVS(1074)*X(68)-JVS(1075)*X(69)-JVS(1076)*X(70)-JVS(1077)*X(71)&
            &-JVS(1078)*X(72)-JVS(1079)*X(73)-JVS(1080)*X(74)-JVS(1081)*X(75)-JVS(1082)*X(76)-JVS(1083)*X(77)-JVS(1084)&
            &*X(78)-JVS(1085)*X(79)
  X(80) = X(80)/JVS(1086)
  X(79) = (X(79)-JVS(1032)*X(80))/(JVS(1031))
  X(78) = (X(78)-JVS(993)*X(79)-JVS(994)*X(80))/(JVS(992))
  X(77) = (X(77)-JVS(981)*X(78)-JVS(982)*X(79)-JVS(983)*X(80))/(JVS(980))
  X(76) = (X(76)-JVS(954)*X(77)-JVS(955)*X(78)-JVS(956)*X(79)-JVS(957)*X(80))/(JVS(953))
  X(75) = (X(75)-JVS(917)*X(76)-JVS(918)*X(77)-JVS(919)*X(78)-JVS(920)*X(79)-JVS(921)*X(80))/(JVS(916))
  X(74) = (X(74)-JVS(853)*X(75)-JVS(854)*X(76)-JVS(855)*X(77)-JVS(856)*X(78)-JVS(857)*X(79)-JVS(858)*X(80))/(JVS(852))
  X(73) = (X(73)-JVS(807)*X(74)-JVS(808)*X(75)-JVS(809)*X(76)-JVS(810)*X(77)-JVS(811)*X(78)-JVS(812)*X(79)-JVS(813)&
            &*X(80))/(JVS(806))
  X(72) = (X(72)-JVS(737)*X(73)-JVS(738)*X(74)-JVS(739)*X(75)-JVS(740)*X(76)-JVS(741)*X(77)-JVS(742)*X(78)-JVS(743)&
            &*X(79)-JVS(744)*X(80))/(JVS(736))
  X(71) = (X(71)-JVS(694)*X(72)-JVS(695)*X(73)-JVS(696)*X(74)-JVS(697)*X(75)-JVS(698)*X(76)-JVS(699)*X(77)-JVS(700)&
            &*X(78)-JVS(701)*X(79)-JVS(702)*X(80))/(JVS(693))
  X(70) = (X(70)-JVS(678)*X(71)-JVS(679)*X(72)-JVS(680)*X(73)-JVS(681)*X(74)-JVS(682)*X(75)-JVS(683)*X(76)-JVS(684)&
            &*X(77)-JVS(685)*X(78)-JVS(686)*X(79)-JVS(687)*X(80))/(JVS(677))
  X(69) = (X(69)-JVS(651)*X(70)-JVS(652)*X(71)-JVS(653)*X(72)-JVS(654)*X(73)-JVS(655)*X(74)-JVS(656)*X(75)-JVS(657)&
            &*X(76)-JVS(658)*X(77)-JVS(659)*X(78)-JVS(660)*X(79)-JVS(661)*X(80))/(JVS(650))
  X(68) = (X(68)-JVS(625)*X(69)-JVS(626)*X(70)-JVS(627)*X(71)-JVS(628)*X(72)-JVS(629)*X(73)-JVS(630)*X(74)-JVS(631)&
            &*X(75)-JVS(632)*X(76)-JVS(633)*X(77)-JVS(634)*X(78)-JVS(635)*X(79)-JVS(636)*X(80))/(JVS(624))
  X(67) = (X(67)-JVS(604)*X(72)-JVS(605)*X(73)-JVS(606)*X(74)-JVS(607)*X(75)-JVS(608)*X(76)-JVS(609)*X(77)-JVS(610)&
            &*X(79))/(JVS(603))
  X(66) = (X(66)-JVS(588)*X(67)-JVS(589)*X(72)-JVS(590)*X(73)-JVS(591)*X(74)-JVS(592)*X(75)-JVS(593)*X(76)-JVS(594)&
            &*X(77)-JVS(595)*X(79))/(JVS(587))
  X(65) = (X(65)-JVS(574)*X(72)-JVS(575)*X(73)-JVS(576)*X(74)-JVS(577)*X(75)-JVS(578)*X(76)-JVS(579)*X(77)-JVS(580)&
            &*X(79))/(JVS(573))
  X(64) = (X(64)-JVS(555)*X(65)-JVS(556)*X(66)-JVS(557)*X(67)-JVS(558)*X(68)-JVS(559)*X(69)-JVS(560)*X(70)-JVS(561)&
            &*X(71)-JVS(562)*X(72)-JVS(563)*X(73)-JVS(564)*X(74)-JVS(565)*X(75)-JVS(566)*X(76)-JVS(567)*X(77)-JVS(568)*X(78)&
            &-JVS(569)*X(79)-JVS(570)*X(80))/(JVS(554))
  X(63) = (X(63)-JVS(534)*X(72)-JVS(535)*X(73)-JVS(536)*X(74)-JVS(537)*X(75)-JVS(538)*X(76)-JVS(539)*X(77)-JVS(540)&
            &*X(79))/(JVS(533))
  X(62) = (X(62)-JVS(518)*X(69)-JVS(519)*X(70)-JVS(520)*X(71)-JVS(521)*X(72)-JVS(522)*X(73)-JVS(523)*X(74)-JVS(524)&
            &*X(75)-JVS(525)*X(76)-JVS(526)*X(77)-JVS(527)*X(78)-JVS(528)*X(79)-JVS(529)*X(80))/(JVS(517))
  X(61) = (X(61)-JVS(499)*X(72)-JVS(500)*X(73)-JVS(501)*X(74)-JVS(502)*X(75)-JVS(503)*X(76)-JVS(504)*X(77)-JVS(505)&
            &*X(79)-JVS(506)*X(80))/(JVS(498))
  X(60) = (X(60)-JVS(484)*X(66)-JVS(485)*X(67)-JVS(486)*X(71)-JVS(487)*X(72)-JVS(488)*X(73)-JVS(489)*X(74)-JVS(490)&
            &*X(75)-JVS(491)*X(76)-JVS(492)*X(77)-JVS(493)*X(78)-JVS(494)*X(79)-JVS(495)*X(80))/(JVS(483))
  X(59) = (X(59)-JVS(467)*X(72)-JVS(468)*X(73)-JVS(469)*X(74)-JVS(470)*X(75)-JVS(471)*X(76)-JVS(472)*X(77)-JVS(473)&
            &*X(79))/(JVS(466))
  X(58) = (X(58)-JVS(457)*X(72)-JVS(458)*X(73)-JVS(459)*X(74)-JVS(460)*X(75)-JVS(461)*X(76)-JVS(462)*X(77)-JVS(463)&
            &*X(79)-JVS(464)*X(80))/(JVS(456))
  X(57) = (X(57)-JVS(447)*X(72)-JVS(448)*X(73)-JVS(449)*X(74)-JVS(450)*X(75)-JVS(451)*X(76)-JVS(452)*X(77)-JVS(453)&
            &*X(79)-JVS(454)*X(80))/(JVS(446))
  X(56) = (X(56)-JVS(424)*X(57)-JVS(425)*X(58)-JVS(426)*X(59)-JVS(427)*X(61)-JVS(428)*X(62)-JVS(429)*X(63)-JVS(430)&
            &*X(65)-JVS(431)*X(66)-JVS(432)*X(67)-JVS(433)*X(68)-JVS(434)*X(69)-JVS(435)*X(71)-JVS(436)*X(72)-JVS(437)*X(73)&
            &-JVS(438)*X(74)-JVS(439)*X(75)-JVS(440)*X(76)-JVS(441)*X(77)-JVS(442)*X(78)-JVS(443)*X(79)-JVS(444)*X(80))&
            &/(JVS(423))
  X(55) = (X(55)-JVS(394)*X(72)-JVS(395)*X(73)-JVS(396)*X(74)-JVS(397)*X(75)-JVS(398)*X(76)-JVS(399)*X(79))/(JVS(393))
  X(54) = (X(54)-JVS(385)*X(72)-JVS(386)*X(73)-JVS(387)*X(74)-JVS(388)*X(75)-JVS(389)*X(76)-JVS(390)*X(77)-JVS(391)&
            &*X(79))/(JVS(384))
  X(53) = (X(53)-JVS(377)*X(72)-JVS(378)*X(73)-JVS(379)*X(74)-JVS(380)*X(75)-JVS(381)*X(76)-JVS(382)*X(79))/(JVS(376))
  X(52) = (X(52)-JVS(367)*X(72)-JVS(368)*X(73)-JVS(369)*X(74)-JVS(370)*X(75)-JVS(371)*X(76)-JVS(372)*X(77)-JVS(373)&
            &*X(79)-JVS(374)*X(80))/(JVS(366))
  X(51) = (X(51)-JVS(348)*X(52)-JVS(349)*X(54)-JVS(350)*X(57)-JVS(351)*X(58)-JVS(352)*X(61)-JVS(353)*X(69)-JVS(354)&
            &*X(71)-JVS(355)*X(72)-JVS(356)*X(73)-JVS(357)*X(74)-JVS(358)*X(75)-JVS(359)*X(76)-JVS(360)*X(77)-JVS(361)*X(78)&
            &-JVS(362)*X(79)-JVS(363)*X(80))/(JVS(347))
  X(50) = (X(50)-JVS(332)*X(72)-JVS(333)*X(73)-JVS(334)*X(74)-JVS(335)*X(75)-JVS(336)*X(76)-JVS(337)*X(77)-JVS(338)&
            &*X(79))/(JVS(331))
  X(49) = (X(49)-JVS(310)*X(53)-JVS(311)*X(55)-JVS(312)*X(59)-JVS(313)*X(60)-JVS(314)*X(61)-JVS(315)*X(63)-JVS(316)&
            &*X(65)-JVS(317)*X(66)-JVS(318)*X(67)-JVS(319)*X(69)-JVS(320)*X(70)-JVS(321)*X(71)-JVS(322)*X(72)-JVS(323)*X(73)&
            &-JVS(324)*X(74)-JVS(325)*X(75)-JVS(326)*X(76)-JVS(327)*X(77)-JVS(328)*X(79))/(JVS(309))
  X(48) = (X(48)-JVS(295)*X(50)-JVS(296)*X(54)-JVS(297)*X(72)-JVS(298)*X(73)-JVS(299)*X(74)-JVS(300)*X(75)-JVS(301)&
            &*X(76)-JVS(302)*X(77)-JVS(303)*X(79))/(JVS(294))
  X(47) = (X(47)-JVS(283)*X(65)-JVS(284)*X(72)-JVS(285)*X(73)-JVS(286)*X(74)-JVS(287)*X(75)-JVS(288)*X(76)-JVS(289)&
            &*X(79))/(JVS(282))
  X(46) = (X(46)-JVS(275)*X(63)-JVS(276)*X(72)-JVS(277)*X(73)-JVS(278)*X(74)-JVS(279)*X(75)-JVS(280)*X(76)-JVS(281)&
            &*X(79))/(JVS(274))
  X(45) = (X(45)-JVS(267)*X(72)-JVS(268)*X(73)-JVS(269)*X(74)-JVS(270)*X(75)-JVS(271)*X(76)-JVS(272)*X(77)-JVS(273)&
            &*X(79))/(JVS(266))
  X(44) = (X(44)-JVS(254)*X(52)-JVS(255)*X(57)-JVS(256)*X(58)-JVS(257)*X(61)-JVS(258)*X(71)-JVS(259)*X(72)-JVS(260)&
            &*X(73)-JVS(261)*X(74)-JVS(262)*X(76)-JVS(263)*X(77)-JVS(264)*X(79))/(JVS(253))
  X(43) = (X(43)-JVS(247)*X(73)-JVS(248)*X(74)-JVS(249)*X(77)-JVS(250)*X(78)-JVS(251)*X(80))/(JVS(246))
  X(42) = (X(42)-JVS(239)*X(48)-JVS(240)*X(73)-JVS(241)*X(74)-JVS(242)*X(75)-JVS(243)*X(77)-JVS(244)*X(79)-JVS(245)&
            &*X(80))/(JVS(238))
  X(41) = (X(41)-JVS(232)*X(73)-JVS(233)*X(74)-JVS(234)*X(77))/(JVS(231))
  X(40) = (X(40)-JVS(226)*X(73)-JVS(227)*X(74)-JVS(228)*X(75)-JVS(229)*X(77)-JVS(230)*X(80))/(JVS(225))
  X(39) = (X(39)-JVS(218)*X(73)-JVS(219)*X(74)-JVS(220)*X(77))/(JVS(217))
  X(38) = (X(38)-JVS(214)*X(73)-JVS(215)*X(74)-JVS(216)*X(77))/(JVS(213))
  X(37) = (X(37)-JVS(195)*X(38)-JVS(196)*X(39)-JVS(197)*X(41)-JVS(198)*X(42)-JVS(199)*X(43)-JVS(200)*X(44)-JVS(201)&
            &*X(48)-JVS(202)*X(50)-JVS(203)*X(51)-JVS(204)*X(56)-JVS(205)*X(61)-JVS(206)*X(63)-JVS(207)*X(64)-JVS(208)*X(65)&
            &-JVS(209)*X(73)-JVS(210)*X(74)-JVS(211)*X(77)-JVS(212)*X(79))/(JVS(194))
  X(36) = (X(36)-JVS(186)*X(73)-JVS(187)*X(74)-JVS(188)*X(77))/(JVS(185))
  X(35) = (X(35)-JVS(179)*X(43)-JVS(180)*X(72)-JVS(181)*X(73)-JVS(182)*X(74)-JVS(183)*X(77)-JVS(184)*X(80))/(JVS(178))
  X(34) = (X(34)-JVS(173)*X(39)-JVS(174)*X(50)-JVS(175)*X(73)-JVS(176)*X(74)-JVS(177)*X(79))/(JVS(172))
  X(33) = (X(33)-JVS(167)*X(61)-JVS(168)*X(72)-JVS(169)*X(73)-JVS(170)*X(75)-JVS(171)*X(77))/(JVS(166))
  X(32) = (X(32)-JVS(156)*X(40)-JVS(157)*X(44)-JVS(158)*X(51)-JVS(159)*X(56)-JVS(160)*X(61)-JVS(161)*X(64)-JVS(162)&
            &*X(73)-JVS(163)*X(74)-JVS(164)*X(75)-JVS(165)*X(80))/(JVS(155))
  X(31) = (X(31)-JVS(151)*X(40)-JVS(152)*X(73)-JVS(153)*X(77)-JVS(154)*X(80))/(JVS(150))
  X(30) = (X(30)-JVS(147)*X(73)-JVS(148)*X(74)-JVS(149)*X(77))/(JVS(146))
  X(29) = (X(29)-JVS(143)*X(73)-JVS(144)*X(77)-JVS(145)*X(80))/(JVS(142))
  X(28) = (X(28)-JVS(138)*X(73)-JVS(139)*X(77)-JVS(140)*X(80))/(JVS(137))
  X(27) = (X(27)-JVS(131)*X(40)-JVS(132)*X(73)-JVS(133)*X(74)-JVS(134)*X(75)-JVS(135)*X(80))/(JVS(130))
  X(26) = (X(26)-JVS(125)*X(41)-JVS(126)*X(74)-JVS(127)*X(77)-JVS(128)*X(79)-JVS(129)*X(80))/(JVS(124))
  X(25) = (X(25)-JVS(119)*X(43)-JVS(120)*X(69)-JVS(121)*X(73)-JVS(122)*X(76))/(JVS(118))
  X(24) = (X(24)-JVS(108)*X(36)-JVS(109)*X(38)-JVS(110)*X(39)-JVS(111)*X(41)-JVS(112)*X(63)-JVS(113)*X(65)-JVS(114)&
            &*X(73)-JVS(115)*X(75)-JVS(116)*X(77))/(JVS(107))
  X(23) = (X(23)-JVS(101)*X(28)-JVS(102)*X(29)-JVS(103)*X(31)-JVS(104)*X(73)-JVS(105)*X(79)-JVS(106)*X(80))/(JVS(100))
  X(22) = (X(22)-JVS(95)*X(34)-JVS(96)*X(42)-JVS(97)*X(73)-JVS(98)*X(79)-JVS(99)*X(80))/(JVS(94))
  X(21) = (X(21)-JVS(90)*X(42)-JVS(91)*X(73)-JVS(92)*X(80))/(JVS(89))
  X(20) = (X(20)-JVS(86)*X(73)-JVS(87)*X(75)-JVS(88)*X(76))/(JVS(85))
  X(19) = (X(19)-JVS(82)*X(73)-JVS(83)*X(75)-JVS(84)*X(80))/(JVS(81))
  X(18) = (X(18)-JVS(79)*X(34)-JVS(80)*X(73))/(JVS(78))
  X(17) = (X(17)-JVS(74)*X(63)-JVS(75)*X(65)-JVS(76)*X(73)-JVS(77)*X(77))/(JVS(73))
  X(16) = (X(16)-JVS(71)*X(61)-JVS(72)*X(73))/(JVS(70))
  X(15) = (X(15)-JVS(66)*X(63)-JVS(67)*X(65)-JVS(68)*X(73)-JVS(69)*X(77))/(JVS(65))
  X(14) = (X(14)-JVS(64)*X(73))/(JVS(63))
  X(13) = (X(13)-JVS(60)*X(42)-JVS(61)*X(73)-JVS(62)*X(75))/(JVS(59))
  X(12) = (X(12)-JVS(58)*X(73))/(JVS(57))
  X(11) = (X(11)-JVS(55)*X(74)-JVS(56)*X(80))/(JVS(54))
  X(10) = (X(10)-JVS(53)*X(73))/(JVS(52))
  X(9) = (X(9)-JVS(51)*X(73))/(JVS(50))
  X(8) = (X(8)-JVS(49)*X(73))/(JVS(48))
  X(7) = (X(7)-JVS(45)*X(50)-JVS(46)*X(73)-JVS(47)*X(75))/(JVS(44))
  X(6) = (X(6)-JVS(43)*X(77))/(JVS(42))
  X(5) = (X(5)-JVS(41)*X(73))/(JVS(40))
  X(4) = (X(4)-JVS(20)*X(38)-JVS(21)*X(45)-JVS(22)*X(46)-JVS(23)*X(47)-JVS(24)*X(50)-JVS(25)*X(53)-JVS(26)*X(55)-JVS(27)&
           &*X(61)-JVS(28)*X(63)-JVS(29)*X(65)-JVS(30)*X(66)-JVS(31)*X(67)-JVS(32)*X(68)-JVS(33)*X(69)-JVS(34)*X(71)-JVS(35)&
           &*X(72)-JVS(36)*X(75)-JVS(37)*X(76)-JVS(38)*X(77)-JVS(39)*X(78))/(JVS(19))
  X(3) = (X(3)-JVS(8)*X(14)-JVS(9)*X(30)-JVS(10)*X(38)-JVS(11)*X(39)-JVS(12)*X(41)-JVS(13)*X(43)-JVS(14)*X(48)-JVS(15)&
           &*X(61)-JVS(16)*X(63)-JVS(17)*X(73)-JVS(18)*X(77))/(JVS(7))
  X(2) = (X(2)-JVS(5)*X(37)-JVS(6)*X(73))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(5)-JVS(3)*X(73))/(JVS(1))
      
END SUBROUTINE racm_esrlsorg_KppSolve
























      SUBROUTINE racm_esrlsorg_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE racm_esrlsorg_WCOPY



      SUBROUTINE racm_esrlsorg_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE racm_esrlsorg_WAXPY




      SUBROUTINE racm_esrlsorg_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE racm_esrlsorg_WSCAL


      REAL(kind=dp) FUNCTION racm_esrlsorg_WLAMCH( C )








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
          CALL racm_esrlsorg_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      racm_esrlsorg_WLAMCH = Eps

      END FUNCTION racm_esrlsorg_WLAMCH
     
      SUBROUTINE racm_esrlsorg_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE racm_esrlsorg_WLAMCH_ADD




      SUBROUTINE racm_esrlsorg_SET2ZERO(N,Y)




      
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

      END SUBROUTINE racm_esrlsorg_SET2ZERO



      REAL(kind=dp) FUNCTION racm_esrlsorg_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      racm_esrlsorg_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        racm_esrlsorg_WDOT = racm_esrlsorg_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         racm_esrlsorg_WDOT = racm_esrlsorg_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          racm_esrlsorg_WDOT = racm_esrlsorg_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        racm_esrlsorg_WDOT = racm_esrlsorg_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION racm_esrlsorg_WDOT                                          




   SUBROUTINE decomp_racm_esrlsorg( JVS, IER )
   
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
   W( 73 ) = JVS( 3 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 5 )
  JVS( 3) = W( 73 )
  IF ( ABS(  JVS( 4 )) < TINY(a) ) THEN
         IER = 2                                       
         RETURN
  END IF
   W( 2 ) = JVS( 4 )
   W( 37 ) = JVS( 5 )
   W( 73 ) = JVS( 6 )
  JVS( 4) = W( 2 )
  JVS( 5) = W( 37 )
  JVS( 6) = W( 73 )
  IF ( ABS(  JVS( 7 )) < TINY(a) ) THEN
         IER = 3                                       
         RETURN
  END IF
   W( 3 ) = JVS( 7 )
   W( 14 ) = JVS( 8 )
   W( 30 ) = JVS( 9 )
   W( 38 ) = JVS( 10 )
   W( 39 ) = JVS( 11 )
   W( 41 ) = JVS( 12 )
   W( 43 ) = JVS( 13 )
   W( 48 ) = JVS( 14 )
   W( 61 ) = JVS( 15 )
   W( 63 ) = JVS( 16 )
   W( 73 ) = JVS( 17 )
   W( 77 ) = JVS( 18 )
  JVS( 7) = W( 3 )
  JVS( 8) = W( 14 )
  JVS( 9) = W( 30 )
  JVS( 10) = W( 38 )
  JVS( 11) = W( 39 )
  JVS( 12) = W( 41 )
  JVS( 13) = W( 43 )
  JVS( 14) = W( 48 )
  JVS( 15) = W( 61 )
  JVS( 16) = W( 63 )
  JVS( 17) = W( 73 )
  JVS( 18) = W( 77 )
  IF ( ABS(  JVS( 19 )) < TINY(a) ) THEN
         IER = 4                                       
         RETURN
  END IF
   W( 4 ) = JVS( 19 )
   W( 38 ) = JVS( 20 )
   W( 45 ) = JVS( 21 )
   W( 46 ) = JVS( 22 )
   W( 47 ) = JVS( 23 )
   W( 50 ) = JVS( 24 )
   W( 53 ) = JVS( 25 )
   W( 55 ) = JVS( 26 )
   W( 61 ) = JVS( 27 )
   W( 63 ) = JVS( 28 )
   W( 65 ) = JVS( 29 )
   W( 66 ) = JVS( 30 )
   W( 67 ) = JVS( 31 )
   W( 68 ) = JVS( 32 )
   W( 69 ) = JVS( 33 )
   W( 71 ) = JVS( 34 )
   W( 72 ) = JVS( 35 )
   W( 75 ) = JVS( 36 )
   W( 76 ) = JVS( 37 )
   W( 77 ) = JVS( 38 )
   W( 78 ) = JVS( 39 )
  JVS( 19) = W( 4 )
  JVS( 20) = W( 38 )
  JVS( 21) = W( 45 )
  JVS( 22) = W( 46 )
  JVS( 23) = W( 47 )
  JVS( 24) = W( 50 )
  JVS( 25) = W( 53 )
  JVS( 26) = W( 55 )
  JVS( 27) = W( 61 )
  JVS( 28) = W( 63 )
  JVS( 29) = W( 65 )
  JVS( 30) = W( 66 )
  JVS( 31) = W( 67 )
  JVS( 32) = W( 68 )
  JVS( 33) = W( 69 )
  JVS( 34) = W( 71 )
  JVS( 35) = W( 72 )
  JVS( 36) = W( 75 )
  JVS( 37) = W( 76 )
  JVS( 38) = W( 77 )
  JVS( 39) = W( 78 )
  IF ( ABS(  JVS( 40 )) < TINY(a) ) THEN
         IER = 5                                       
         RETURN
  END IF
   W( 5 ) = JVS( 40 )
   W( 73 ) = JVS( 41 )
  JVS( 40) = W( 5 )
  JVS( 41) = W( 73 )
  IF ( ABS(  JVS( 42 )) < TINY(a) ) THEN
         IER = 6                                       
         RETURN
  END IF
   W( 6 ) = JVS( 42 )
   W( 77 ) = JVS( 43 )
  JVS( 42) = W( 6 )
  JVS( 43) = W( 77 )
  IF ( ABS(  JVS( 44 )) < TINY(a) ) THEN
         IER = 7                                       
         RETURN
  END IF
   W( 7 ) = JVS( 44 )
   W( 50 ) = JVS( 45 )
   W( 73 ) = JVS( 46 )
   W( 75 ) = JVS( 47 )
  JVS( 44) = W( 7 )
  JVS( 45) = W( 50 )
  JVS( 46) = W( 73 )
  JVS( 47) = W( 75 )
  IF ( ABS(  JVS( 48 )) < TINY(a) ) THEN
         IER = 8                                       
         RETURN
  END IF
   W( 8 ) = JVS( 48 )
   W( 73 ) = JVS( 49 )
  JVS( 48) = W( 8 )
  JVS( 49) = W( 73 )
  IF ( ABS(  JVS( 50 )) < TINY(a) ) THEN
         IER = 9                                       
         RETURN
  END IF
   W( 9 ) = JVS( 50 )
   W( 73 ) = JVS( 51 )
  JVS( 50) = W( 9 )
  JVS( 51) = W( 73 )
  IF ( ABS(  JVS( 52 )) < TINY(a) ) THEN
         IER = 10                                      
         RETURN
  END IF
   W( 10 ) = JVS( 52 )
   W( 73 ) = JVS( 53 )
  JVS( 52) = W( 10 )
  JVS( 53) = W( 73 )
  IF ( ABS(  JVS( 54 )) < TINY(a) ) THEN
         IER = 11                                      
         RETURN
  END IF
   W( 11 ) = JVS( 54 )
   W( 74 ) = JVS( 55 )
   W( 80 ) = JVS( 56 )
  JVS( 54) = W( 11 )
  JVS( 55) = W( 74 )
  JVS( 56) = W( 80 )
  IF ( ABS(  JVS( 57 )) < TINY(a) ) THEN
         IER = 12                                      
         RETURN
  END IF
   W( 12 ) = JVS( 57 )
   W( 73 ) = JVS( 58 )
  JVS( 57) = W( 12 )
  JVS( 58) = W( 73 )
  IF ( ABS(  JVS( 59 )) < TINY(a) ) THEN
         IER = 13                                      
         RETURN
  END IF
   W( 13 ) = JVS( 59 )
   W( 42 ) = JVS( 60 )
   W( 73 ) = JVS( 61 )
   W( 75 ) = JVS( 62 )
  JVS( 59) = W( 13 )
  JVS( 60) = W( 42 )
  JVS( 61) = W( 73 )
  JVS( 62) = W( 75 )
  IF ( ABS(  JVS( 63 )) < TINY(a) ) THEN
         IER = 14                                      
         RETURN
  END IF
   W( 14 ) = JVS( 63 )
   W( 73 ) = JVS( 64 )
  JVS( 63) = W( 14 )
  JVS( 64) = W( 73 )
  IF ( ABS(  JVS( 65 )) < TINY(a) ) THEN
         IER = 15                                      
         RETURN
  END IF
   W( 15 ) = JVS( 65 )
   W( 63 ) = JVS( 66 )
   W( 65 ) = JVS( 67 )
   W( 73 ) = JVS( 68 )
   W( 77 ) = JVS( 69 )
  JVS( 65) = W( 15 )
  JVS( 66) = W( 63 )
  JVS( 67) = W( 65 )
  JVS( 68) = W( 73 )
  JVS( 69) = W( 77 )
  IF ( ABS(  JVS( 70 )) < TINY(a) ) THEN
         IER = 16                                      
         RETURN
  END IF
   W( 16 ) = JVS( 70 )
   W( 61 ) = JVS( 71 )
   W( 73 ) = JVS( 72 )
  JVS( 70) = W( 16 )
  JVS( 71) = W( 61 )
  JVS( 72) = W( 73 )
  IF ( ABS(  JVS( 73 )) < TINY(a) ) THEN
         IER = 17                                      
         RETURN
  END IF
   W( 17 ) = JVS( 73 )
   W( 63 ) = JVS( 74 )
   W( 65 ) = JVS( 75 )
   W( 73 ) = JVS( 76 )
   W( 77 ) = JVS( 77 )
  JVS( 73) = W( 17 )
  JVS( 74) = W( 63 )
  JVS( 75) = W( 65 )
  JVS( 76) = W( 73 )
  JVS( 77) = W( 77 )
  IF ( ABS(  JVS( 78 )) < TINY(a) ) THEN
         IER = 18                                      
         RETURN
  END IF
   W( 18 ) = JVS( 78 )
   W( 34 ) = JVS( 79 )
   W( 73 ) = JVS( 80 )
  JVS( 78) = W( 18 )
  JVS( 79) = W( 34 )
  JVS( 80) = W( 73 )
  IF ( ABS(  JVS( 81 )) < TINY(a) ) THEN
         IER = 19                                      
         RETURN
  END IF
   W( 19 ) = JVS( 81 )
   W( 73 ) = JVS( 82 )
   W( 75 ) = JVS( 83 )
   W( 80 ) = JVS( 84 )
  JVS( 81) = W( 19 )
  JVS( 82) = W( 73 )
  JVS( 83) = W( 75 )
  JVS( 84) = W( 80 )
  IF ( ABS(  JVS( 85 )) < TINY(a) ) THEN
         IER = 20                                      
         RETURN
  END IF
   W( 20 ) = JVS( 85 )
   W( 73 ) = JVS( 86 )
   W( 75 ) = JVS( 87 )
   W( 76 ) = JVS( 88 )
  JVS( 85) = W( 20 )
  JVS( 86) = W( 73 )
  JVS( 87) = W( 75 )
  JVS( 88) = W( 76 )
  IF ( ABS(  JVS( 89 )) < TINY(a) ) THEN
         IER = 21                                      
         RETURN
  END IF
   W( 21 ) = JVS( 89 )
   W( 42 ) = JVS( 90 )
   W( 73 ) = JVS( 91 )
   W( 80 ) = JVS( 92 )
  JVS( 89) = W( 21 )
  JVS( 90) = W( 42 )
  JVS( 91) = W( 73 )
  JVS( 92) = W( 80 )
  IF ( ABS(  JVS( 94 )) < TINY(a) ) THEN
         IER = 22                                      
         RETURN
  END IF
   W( 21 ) = JVS( 93 )
   W( 22 ) = JVS( 94 )
   W( 34 ) = JVS( 95 )
   W( 42 ) = JVS( 96 )
   W( 73 ) = JVS( 97 )
   W( 79 ) = JVS( 98 )
   W( 80 ) = JVS( 99 )
  a = -W( 21 ) / JVS(           89  )
  W( 21 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 90 )
  W( 73 ) = W( 73 ) + a*JVS( 91 )
  W( 80 ) = W( 80 ) + a*JVS( 92 )
  JVS( 93) = W( 21 )
  JVS( 94) = W( 22 )
  JVS( 95) = W( 34 )
  JVS( 96) = W( 42 )
  JVS( 97) = W( 73 )
  JVS( 98) = W( 79 )
  JVS( 99) = W( 80 )
  IF ( ABS(  JVS( 100 )) < TINY(a) ) THEN
         IER = 23                                      
         RETURN
  END IF
   W( 23 ) = JVS( 100 )
   W( 28 ) = JVS( 101 )
   W( 29 ) = JVS( 102 )
   W( 31 ) = JVS( 103 )
   W( 73 ) = JVS( 104 )
   W( 79 ) = JVS( 105 )
   W( 80 ) = JVS( 106 )
  JVS( 100) = W( 23 )
  JVS( 101) = W( 28 )
  JVS( 102) = W( 29 )
  JVS( 103) = W( 31 )
  JVS( 104) = W( 73 )
  JVS( 105) = W( 79 )
  JVS( 106) = W( 80 )
  IF ( ABS(  JVS( 107 )) < TINY(a) ) THEN
         IER = 24                                      
         RETURN
  END IF
   W( 24 ) = JVS( 107 )
   W( 36 ) = JVS( 108 )
   W( 38 ) = JVS( 109 )
   W( 39 ) = JVS( 110 )
   W( 41 ) = JVS( 111 )
   W( 63 ) = JVS( 112 )
   W( 65 ) = JVS( 113 )
   W( 73 ) = JVS( 114 )
   W( 75 ) = JVS( 115 )
   W( 77 ) = JVS( 116 )
  JVS( 107) = W( 24 )
  JVS( 108) = W( 36 )
  JVS( 109) = W( 38 )
  JVS( 110) = W( 39 )
  JVS( 111) = W( 41 )
  JVS( 112) = W( 63 )
  JVS( 113) = W( 65 )
  JVS( 114) = W( 73 )
  JVS( 115) = W( 75 )
  JVS( 116) = W( 77 )
  IF ( ABS(  JVS( 118 )) < TINY(a) ) THEN
         IER = 25                                      
         RETURN
  END IF
   W( 12 ) = JVS( 117 )
   W( 25 ) = JVS( 118 )
   W( 43 ) = JVS( 119 )
   W( 69 ) = JVS( 120 )
   W( 73 ) = JVS( 121 )
   W( 76 ) = JVS( 122 )
  a = -W( 12 ) / JVS(           57  )
  W( 12 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 58 )
  JVS( 117) = W( 12 )
  JVS( 118) = W( 25 )
  JVS( 119) = W( 43 )
  JVS( 120) = W( 69 )
  JVS( 121) = W( 73 )
  JVS( 122) = W( 76 )
  IF ( ABS(  JVS( 124 )) < TINY(a) ) THEN
         IER = 26                                      
         RETURN
  END IF
   W( 6 ) = JVS( 123 )
   W( 26 ) = JVS( 124 )
   W( 41 ) = JVS( 125 )
   W( 74 ) = JVS( 126 )
   W( 77 ) = JVS( 127 )
   W( 79 ) = JVS( 128 )
   W( 80 ) = JVS( 129 )
  a = -W( 6 ) / JVS(           42  )
  W( 6 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 43 )
  JVS( 123) = W( 6 )
  JVS( 124) = W( 26 )
  JVS( 125) = W( 41 )
  JVS( 126) = W( 74 )
  JVS( 127) = W( 77 )
  JVS( 128) = W( 79 )
  JVS( 129) = W( 80 )
  IF ( ABS(  JVS( 130 )) < TINY(a) ) THEN
         IER = 27                                      
         RETURN
  END IF
   W( 27 ) = JVS( 130 )
   W( 40 ) = JVS( 131 )
   W( 73 ) = JVS( 132 )
   W( 74 ) = JVS( 133 )
   W( 75 ) = JVS( 134 )
   W( 80 ) = JVS( 135 )
  JVS( 130) = W( 27 )
  JVS( 131) = W( 40 )
  JVS( 132) = W( 73 )
  JVS( 133) = W( 74 )
  JVS( 134) = W( 75 )
  JVS( 135) = W( 80 )
  IF ( ABS(  JVS( 137 )) < TINY(a) ) THEN
         IER = 28                                      
         RETURN
  END IF
   W( 9 ) = JVS( 136 )
   W( 28 ) = JVS( 137 )
   W( 73 ) = JVS( 138 )
   W( 77 ) = JVS( 139 )
   W( 80 ) = JVS( 140 )
  a = -W( 9 ) / JVS(           50  )
  W( 9 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 51 )
  JVS( 136) = W( 9 )
  JVS( 137) = W( 28 )
  JVS( 138) = W( 73 )
  JVS( 139) = W( 77 )
  JVS( 140) = W( 80 )
  IF ( ABS(  JVS( 142 )) < TINY(a) ) THEN
         IER = 29                                      
         RETURN
  END IF
   W( 10 ) = JVS( 141 )
   W( 29 ) = JVS( 142 )
   W( 73 ) = JVS( 143 )
   W( 77 ) = JVS( 144 )
   W( 80 ) = JVS( 145 )
  a = -W( 10 ) / JVS(           52  )
  W( 10 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 53 )
  JVS( 141) = W( 10 )
  JVS( 142) = W( 29 )
  JVS( 143) = W( 73 )
  JVS( 144) = W( 77 )
  JVS( 145) = W( 80 )
  IF ( ABS(  JVS( 146 )) < TINY(a) ) THEN
         IER = 30                                      
         RETURN
  END IF
   W( 30 ) = JVS( 146 )
   W( 73 ) = JVS( 147 )
   W( 74 ) = JVS( 148 )
   W( 77 ) = JVS( 149 )
  JVS( 146) = W( 30 )
  JVS( 147) = W( 73 )
  JVS( 148) = W( 74 )
  JVS( 149) = W( 77 )
  IF ( ABS(  JVS( 150 )) < TINY(a) ) THEN
         IER = 31                                      
         RETURN
  END IF
   W( 31 ) = JVS( 150 )
   W( 40 ) = JVS( 151 )
   W( 73 ) = JVS( 152 )
   W( 77 ) = JVS( 153 )
   W( 80 ) = JVS( 154 )
  JVS( 150) = W( 31 )
  JVS( 151) = W( 40 )
  JVS( 152) = W( 73 )
  JVS( 153) = W( 77 )
  JVS( 154) = W( 80 )
  IF ( ABS(  JVS( 155 )) < TINY(a) ) THEN
         IER = 32                                      
         RETURN
  END IF
   W( 32 ) = JVS( 155 )
   W( 40 ) = JVS( 156 )
   W( 44 ) = JVS( 157 )
   W( 51 ) = JVS( 158 )
   W( 56 ) = JVS( 159 )
   W( 61 ) = JVS( 160 )
   W( 64 ) = JVS( 161 )
   W( 73 ) = JVS( 162 )
   W( 74 ) = JVS( 163 )
   W( 75 ) = JVS( 164 )
   W( 80 ) = JVS( 165 )
  JVS( 155) = W( 32 )
  JVS( 156) = W( 40 )
  JVS( 157) = W( 44 )
  JVS( 158) = W( 51 )
  JVS( 159) = W( 56 )
  JVS( 160) = W( 61 )
  JVS( 161) = W( 64 )
  JVS( 162) = W( 73 )
  JVS( 163) = W( 74 )
  JVS( 164) = W( 75 )
  JVS( 165) = W( 80 )
  IF ( ABS(  JVS( 166 )) < TINY(a) ) THEN
         IER = 33                                      
         RETURN
  END IF
   W( 33 ) = JVS( 166 )
   W( 61 ) = JVS( 167 )
   W( 72 ) = JVS( 168 )
   W( 73 ) = JVS( 169 )
   W( 75 ) = JVS( 170 )
   W( 77 ) = JVS( 171 )
  JVS( 166) = W( 33 )
  JVS( 167) = W( 61 )
  JVS( 168) = W( 72 )
  JVS( 169) = W( 73 )
  JVS( 170) = W( 75 )
  JVS( 171) = W( 77 )
  IF ( ABS(  JVS( 172 )) < TINY(a) ) THEN
         IER = 34                                      
         RETURN
  END IF
   W( 34 ) = JVS( 172 )
   W( 39 ) = JVS( 173 )
   W( 50 ) = JVS( 174 )
   W( 73 ) = JVS( 175 )
   W( 74 ) = JVS( 176 )
   W( 79 ) = JVS( 177 )
  JVS( 172) = W( 34 )
  JVS( 173) = W( 39 )
  JVS( 174) = W( 50 )
  JVS( 175) = W( 73 )
  JVS( 176) = W( 74 )
  JVS( 177) = W( 79 )
  IF ( ABS(  JVS( 178 )) < TINY(a) ) THEN
         IER = 35                                      
         RETURN
  END IF
   W( 35 ) = JVS( 178 )
   W( 43 ) = JVS( 179 )
   W( 72 ) = JVS( 180 )
   W( 73 ) = JVS( 181 )
   W( 74 ) = JVS( 182 )
   W( 77 ) = JVS( 183 )
   W( 80 ) = JVS( 184 )
  JVS( 178) = W( 35 )
  JVS( 179) = W( 43 )
  JVS( 180) = W( 72 )
  JVS( 181) = W( 73 )
  JVS( 182) = W( 74 )
  JVS( 183) = W( 77 )
  JVS( 184) = W( 80 )
  IF ( ABS(  JVS( 185 )) < TINY(a) ) THEN
         IER = 36                                      
         RETURN
  END IF
   W( 36 ) = JVS( 185 )
   W( 73 ) = JVS( 186 )
   W( 74 ) = JVS( 187 )
   W( 77 ) = JVS( 188 )
  JVS( 185) = W( 36 )
  JVS( 186) = W( 73 )
  JVS( 187) = W( 74 )
  JVS( 188) = W( 77 )
  IF ( ABS(  JVS( 194 )) < TINY(a) ) THEN
         IER = 37                                      
         RETURN
  END IF
   W( 14 ) = JVS( 189 )
   W( 18 ) = JVS( 190 )
   W( 30 ) = JVS( 191 )
   W( 34 ) = JVS( 192 )
   W( 36 ) = JVS( 193 )
   W( 37 ) = JVS( 194 )
   W( 38 ) = JVS( 195 )
   W( 39 ) = JVS( 196 )
   W( 41 ) = JVS( 197 )
   W( 42 ) = JVS( 198 )
   W( 43 ) = JVS( 199 )
   W( 44 ) = JVS( 200 )
   W( 48 ) = JVS( 201 )
   W( 50 ) = JVS( 202 )
   W( 51 ) = JVS( 203 )
   W( 56 ) = JVS( 204 )
   W( 61 ) = JVS( 205 )
   W( 63 ) = JVS( 206 )
   W( 64 ) = JVS( 207 )
   W( 65 ) = JVS( 208 )
   W( 73 ) = JVS( 209 )
   W( 74 ) = JVS( 210 )
   W( 77 ) = JVS( 211 )
   W( 79 ) = JVS( 212 )
  a = -W( 14 ) / JVS(           63  )
  W( 14 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 64 )
  a = -W( 18 ) / JVS(           78  )
  W( 18 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 79 )
  W( 73 ) = W( 73 ) + a*JVS( 80 )
  a = -W( 30 ) / JVS(          146  )
  W( 30 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 147 )
  W( 74 ) = W( 74 ) + a*JVS( 148 )
  W( 77 ) = W( 77 ) + a*JVS( 149 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 173 )
  W( 50 ) = W( 50 ) + a*JVS( 174 )
  W( 73 ) = W( 73 ) + a*JVS( 175 )
  W( 74 ) = W( 74 ) + a*JVS( 176 )
  W( 79 ) = W( 79 ) + a*JVS( 177 )
  a = -W( 36 ) / JVS(          185  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 186 )
  W( 74 ) = W( 74 ) + a*JVS( 187 )
  W( 77 ) = W( 77 ) + a*JVS( 188 )
  JVS( 189) = W( 14 )
  JVS( 190) = W( 18 )
  JVS( 191) = W( 30 )
  JVS( 192) = W( 34 )
  JVS( 193) = W( 36 )
  JVS( 194) = W( 37 )
  JVS( 195) = W( 38 )
  JVS( 196) = W( 39 )
  JVS( 197) = W( 41 )
  JVS( 198) = W( 42 )
  JVS( 199) = W( 43 )
  JVS( 200) = W( 44 )
  JVS( 201) = W( 48 )
  JVS( 202) = W( 50 )
  JVS( 203) = W( 51 )
  JVS( 204) = W( 56 )
  JVS( 205) = W( 61 )
  JVS( 206) = W( 63 )
  JVS( 207) = W( 64 )
  JVS( 208) = W( 65 )
  JVS( 209) = W( 73 )
  JVS( 210) = W( 74 )
  JVS( 211) = W( 77 )
  JVS( 212) = W( 79 )
  IF ( ABS(  JVS( 213 )) < TINY(a) ) THEN
         IER = 38                                      
         RETURN
  END IF
   W( 38 ) = JVS( 213 )
   W( 73 ) = JVS( 214 )
   W( 74 ) = JVS( 215 )
   W( 77 ) = JVS( 216 )
  JVS( 213) = W( 38 )
  JVS( 214) = W( 73 )
  JVS( 215) = W( 74 )
  JVS( 216) = W( 77 )
  IF ( ABS(  JVS( 217 )) < TINY(a) ) THEN
         IER = 39                                      
         RETURN
  END IF
   W( 39 ) = JVS( 217 )
   W( 73 ) = JVS( 218 )
   W( 74 ) = JVS( 219 )
   W( 77 ) = JVS( 220 )
  JVS( 217) = W( 39 )
  JVS( 218) = W( 73 )
  JVS( 219) = W( 74 )
  JVS( 220) = W( 77 )
  IF ( ABS(  JVS( 225 )) < TINY(a) ) THEN
         IER = 40                                      
         RETURN
  END IF
   W( 27 ) = JVS( 221 )
   W( 28 ) = JVS( 222 )
   W( 29 ) = JVS( 223 )
   W( 31 ) = JVS( 224 )
   W( 40 ) = JVS( 225 )
   W( 73 ) = JVS( 226 )
   W( 74 ) = JVS( 227 )
   W( 75 ) = JVS( 228 )
   W( 77 ) = JVS( 229 )
   W( 80 ) = JVS( 230 )
  a = -W( 27 ) / JVS(          130  )
  W( 27 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 131 )
  W( 73 ) = W( 73 ) + a*JVS( 132 )
  W( 74 ) = W( 74 ) + a*JVS( 133 )
  W( 75 ) = W( 75 ) + a*JVS( 134 )
  W( 80 ) = W( 80 ) + a*JVS( 135 )
  a = -W( 28 ) / JVS(          137  )
  W( 28 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 138 )
  W( 77 ) = W( 77 ) + a*JVS( 139 )
  W( 80 ) = W( 80 ) + a*JVS( 140 )
  a = -W( 29 ) / JVS(          142  )
  W( 29 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 143 )
  W( 77 ) = W( 77 ) + a*JVS( 144 )
  W( 80 ) = W( 80 ) + a*JVS( 145 )
  a = -W( 31 ) / JVS(          150  )
  W( 31 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 151 )
  W( 73 ) = W( 73 ) + a*JVS( 152 )
  W( 77 ) = W( 77 ) + a*JVS( 153 )
  W( 80 ) = W( 80 ) + a*JVS( 154 )
  JVS( 221) = W( 27 )
  JVS( 222) = W( 28 )
  JVS( 223) = W( 29 )
  JVS( 224) = W( 31 )
  JVS( 225) = W( 40 )
  JVS( 226) = W( 73 )
  JVS( 227) = W( 74 )
  JVS( 228) = W( 75 )
  JVS( 229) = W( 77 )
  JVS( 230) = W( 80 )
  IF ( ABS(  JVS( 231 )) < TINY(a) ) THEN
         IER = 41                                      
         RETURN
  END IF
   W( 41 ) = JVS( 231 )
   W( 73 ) = JVS( 232 )
   W( 74 ) = JVS( 233 )
   W( 77 ) = JVS( 234 )
  JVS( 231) = W( 41 )
  JVS( 232) = W( 73 )
  JVS( 233) = W( 74 )
  JVS( 234) = W( 77 )
  IF ( ABS(  JVS( 238 )) < TINY(a) ) THEN
         IER = 42                                      
         RETURN
  END IF
   W( 13 ) = JVS( 235 )
   W( 21 ) = JVS( 236 )
   W( 39 ) = JVS( 237 )
   W( 42 ) = JVS( 238 )
   W( 48 ) = JVS( 239 )
   W( 73 ) = JVS( 240 )
   W( 74 ) = JVS( 241 )
   W( 75 ) = JVS( 242 )
   W( 77 ) = JVS( 243 )
   W( 79 ) = JVS( 244 )
   W( 80 ) = JVS( 245 )
  a = -W( 13 ) / JVS(           59  )
  W( 13 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 60 )
  W( 73 ) = W( 73 ) + a*JVS( 61 )
  W( 75 ) = W( 75 ) + a*JVS( 62 )
  a = -W( 21 ) / JVS(           89  )
  W( 21 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 90 )
  W( 73 ) = W( 73 ) + a*JVS( 91 )
  W( 80 ) = W( 80 ) + a*JVS( 92 )
  a = -W( 39 ) / JVS(          217  )
  W( 39 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 218 )
  W( 74 ) = W( 74 ) + a*JVS( 219 )
  W( 77 ) = W( 77 ) + a*JVS( 220 )
  JVS( 235) = W( 13 )
  JVS( 236) = W( 21 )
  JVS( 237) = W( 39 )
  JVS( 238) = W( 42 )
  JVS( 239) = W( 48 )
  JVS( 240) = W( 73 )
  JVS( 241) = W( 74 )
  JVS( 242) = W( 75 )
  JVS( 243) = W( 77 )
  JVS( 244) = W( 79 )
  JVS( 245) = W( 80 )
  IF ( ABS(  JVS( 246 )) < TINY(a) ) THEN
         IER = 43                                      
         RETURN
  END IF
   W( 43 ) = JVS( 246 )
   W( 73 ) = JVS( 247 )
   W( 74 ) = JVS( 248 )
   W( 77 ) = JVS( 249 )
   W( 78 ) = JVS( 250 )
   W( 80 ) = JVS( 251 )
  JVS( 246) = W( 43 )
  JVS( 247) = W( 73 )
  JVS( 248) = W( 74 )
  JVS( 249) = W( 77 )
  JVS( 250) = W( 78 )
  JVS( 251) = W( 80 )
  IF ( ABS(  JVS( 253 )) < TINY(a) ) THEN
         IER = 44                                      
         RETURN
  END IF
   W( 14 ) = JVS( 252 )
   W( 44 ) = JVS( 253 )
   W( 52 ) = JVS( 254 )
   W( 57 ) = JVS( 255 )
   W( 58 ) = JVS( 256 )
   W( 61 ) = JVS( 257 )
   W( 71 ) = JVS( 258 )
   W( 72 ) = JVS( 259 )
   W( 73 ) = JVS( 260 )
   W( 74 ) = JVS( 261 )
   W( 76 ) = JVS( 262 )
   W( 77 ) = JVS( 263 )
   W( 79 ) = JVS( 264 )
  a = -W( 14 ) / JVS(           63  )
  W( 14 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 64 )
  JVS( 252) = W( 14 )
  JVS( 253) = W( 44 )
  JVS( 254) = W( 52 )
  JVS( 255) = W( 57 )
  JVS( 256) = W( 58 )
  JVS( 257) = W( 61 )
  JVS( 258) = W( 71 )
  JVS( 259) = W( 72 )
  JVS( 260) = W( 73 )
  JVS( 261) = W( 74 )
  JVS( 262) = W( 76 )
  JVS( 263) = W( 77 )
  JVS( 264) = W( 79 )
  IF ( ABS(  JVS( 266 )) < TINY(a) ) THEN
         IER = 45                                      
         RETURN
  END IF
   W( 30 ) = JVS( 265 )
   W( 45 ) = JVS( 266 )
   W( 72 ) = JVS( 267 )
   W( 73 ) = JVS( 268 )
   W( 74 ) = JVS( 269 )
   W( 75 ) = JVS( 270 )
   W( 76 ) = JVS( 271 )
   W( 77 ) = JVS( 272 )
   W( 79 ) = JVS( 273 )
  a = -W( 30 ) / JVS(          146  )
  W( 30 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 147 )
  W( 74 ) = W( 74 ) + a*JVS( 148 )
  W( 77 ) = W( 77 ) + a*JVS( 149 )
  JVS( 265) = W( 30 )
  JVS( 266) = W( 45 )
  JVS( 267) = W( 72 )
  JVS( 268) = W( 73 )
  JVS( 269) = W( 74 )
  JVS( 270) = W( 75 )
  JVS( 271) = W( 76 )
  JVS( 272) = W( 77 )
  JVS( 273) = W( 79 )
  IF ( ABS(  JVS( 274 )) < TINY(a) ) THEN
         IER = 46                                      
         RETURN
  END IF
   W( 46 ) = JVS( 274 )
   W( 63 ) = JVS( 275 )
   W( 72 ) = JVS( 276 )
   W( 73 ) = JVS( 277 )
   W( 74 ) = JVS( 278 )
   W( 75 ) = JVS( 279 )
   W( 76 ) = JVS( 280 )
   W( 79 ) = JVS( 281 )
  JVS( 274) = W( 46 )
  JVS( 275) = W( 63 )
  JVS( 276) = W( 72 )
  JVS( 277) = W( 73 )
  JVS( 278) = W( 74 )
  JVS( 279) = W( 75 )
  JVS( 280) = W( 76 )
  JVS( 281) = W( 79 )
  IF ( ABS(  JVS( 282 )) < TINY(a) ) THEN
         IER = 47                                      
         RETURN
  END IF
   W( 47 ) = JVS( 282 )
   W( 65 ) = JVS( 283 )
   W( 72 ) = JVS( 284 )
   W( 73 ) = JVS( 285 )
   W( 74 ) = JVS( 286 )
   W( 75 ) = JVS( 287 )
   W( 76 ) = JVS( 288 )
   W( 79 ) = JVS( 289 )
  JVS( 282) = W( 47 )
  JVS( 283) = W( 65 )
  JVS( 284) = W( 72 )
  JVS( 285) = W( 73 )
  JVS( 286) = W( 74 )
  JVS( 287) = W( 75 )
  JVS( 288) = W( 76 )
  JVS( 289) = W( 79 )
  IF ( ABS(  JVS( 294 )) < TINY(a) ) THEN
         IER = 48                                      
         RETURN
  END IF
   W( 7 ) = JVS( 290 )
   W( 38 ) = JVS( 291 )
   W( 39 ) = JVS( 292 )
   W( 41 ) = JVS( 293 )
   W( 48 ) = JVS( 294 )
   W( 50 ) = JVS( 295 )
   W( 54 ) = JVS( 296 )
   W( 72 ) = JVS( 297 )
   W( 73 ) = JVS( 298 )
   W( 74 ) = JVS( 299 )
   W( 75 ) = JVS( 300 )
   W( 76 ) = JVS( 301 )
   W( 77 ) = JVS( 302 )
   W( 79 ) = JVS( 303 )
  a = -W( 7 ) / JVS(           44  )
  W( 7 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 45 )
  W( 73 ) = W( 73 ) + a*JVS( 46 )
  W( 75 ) = W( 75 ) + a*JVS( 47 )
  a = -W( 38 ) / JVS(          213  )
  W( 38 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 214 )
  W( 74 ) = W( 74 ) + a*JVS( 215 )
  W( 77 ) = W( 77 ) + a*JVS( 216 )
  a = -W( 39 ) / JVS(          217  )
  W( 39 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 218 )
  W( 74 ) = W( 74 ) + a*JVS( 219 )
  W( 77 ) = W( 77 ) + a*JVS( 220 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  JVS( 290) = W( 7 )
  JVS( 291) = W( 38 )
  JVS( 292) = W( 39 )
  JVS( 293) = W( 41 )
  JVS( 294) = W( 48 )
  JVS( 295) = W( 50 )
  JVS( 296) = W( 54 )
  JVS( 297) = W( 72 )
  JVS( 298) = W( 73 )
  JVS( 299) = W( 74 )
  JVS( 300) = W( 75 )
  JVS( 301) = W( 76 )
  JVS( 302) = W( 77 )
  JVS( 303) = W( 79 )
  IF ( ABS(  JVS( 309 )) < TINY(a) ) THEN
         IER = 49                                      
         RETURN
  END IF
   W( 8 ) = JVS( 304 )
   W( 16 ) = JVS( 305 )
   W( 36 ) = JVS( 306 )
   W( 46 ) = JVS( 307 )
   W( 47 ) = JVS( 308 )
   W( 49 ) = JVS( 309 )
   W( 53 ) = JVS( 310 )
   W( 55 ) = JVS( 311 )
   W( 59 ) = JVS( 312 )
   W( 60 ) = JVS( 313 )
   W( 61 ) = JVS( 314 )
   W( 63 ) = JVS( 315 )
   W( 65 ) = JVS( 316 )
   W( 66 ) = JVS( 317 )
   W( 67 ) = JVS( 318 )
   W( 69 ) = JVS( 319 )
   W( 70 ) = JVS( 320 )
   W( 71 ) = JVS( 321 )
   W( 72 ) = JVS( 322 )
   W( 73 ) = JVS( 323 )
   W( 74 ) = JVS( 324 )
   W( 75 ) = JVS( 325 )
   W( 76 ) = JVS( 326 )
   W( 77 ) = JVS( 327 )
   W( 79 ) = JVS( 328 )
  a = -W( 8 ) / JVS(           48  )
  W( 8 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 49 )
  a = -W( 16 ) / JVS(           70  )
  W( 16 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 71 )
  W( 73 ) = W( 73 ) + a*JVS( 72 )
  a = -W( 36 ) / JVS(          185  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 186 )
  W( 74 ) = W( 74 ) + a*JVS( 187 )
  W( 77 ) = W( 77 ) + a*JVS( 188 )
  a = -W( 46 ) / JVS(          274  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 275 )
  W( 72 ) = W( 72 ) + a*JVS( 276 )
  W( 73 ) = W( 73 ) + a*JVS( 277 )
  W( 74 ) = W( 74 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 76 ) = W( 76 ) + a*JVS( 280 )
  W( 79 ) = W( 79 ) + a*JVS( 281 )
  a = -W( 47 ) / JVS(          282  )
  W( 47 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 283 )
  W( 72 ) = W( 72 ) + a*JVS( 284 )
  W( 73 ) = W( 73 ) + a*JVS( 285 )
  W( 74 ) = W( 74 ) + a*JVS( 286 )
  W( 75 ) = W( 75 ) + a*JVS( 287 )
  W( 76 ) = W( 76 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  JVS( 304) = W( 8 )
  JVS( 305) = W( 16 )
  JVS( 306) = W( 36 )
  JVS( 307) = W( 46 )
  JVS( 308) = W( 47 )
  JVS( 309) = W( 49 )
  JVS( 310) = W( 53 )
  JVS( 311) = W( 55 )
  JVS( 312) = W( 59 )
  JVS( 313) = W( 60 )
  JVS( 314) = W( 61 )
  JVS( 315) = W( 63 )
  JVS( 316) = W( 65 )
  JVS( 317) = W( 66 )
  JVS( 318) = W( 67 )
  JVS( 319) = W( 69 )
  JVS( 320) = W( 70 )
  JVS( 321) = W( 71 )
  JVS( 322) = W( 72 )
  JVS( 323) = W( 73 )
  JVS( 324) = W( 74 )
  JVS( 325) = W( 75 )
  JVS( 326) = W( 76 )
  JVS( 327) = W( 77 )
  JVS( 328) = W( 79 )
  IF ( ABS(  JVS( 331 )) < TINY(a) ) THEN
         IER = 50                                      
         RETURN
  END IF
   W( 39 ) = JVS( 329 )
   W( 41 ) = JVS( 330 )
   W( 50 ) = JVS( 331 )
   W( 72 ) = JVS( 332 )
   W( 73 ) = JVS( 333 )
   W( 74 ) = JVS( 334 )
   W( 75 ) = JVS( 335 )
   W( 76 ) = JVS( 336 )
   W( 77 ) = JVS( 337 )
   W( 79 ) = JVS( 338 )
  a = -W( 39 ) / JVS(          217  )
  W( 39 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 218 )
  W( 74 ) = W( 74 ) + a*JVS( 219 )
  W( 77 ) = W( 77 ) + a*JVS( 220 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  JVS( 329) = W( 39 )
  JVS( 330) = W( 41 )
  JVS( 331) = W( 50 )
  JVS( 332) = W( 72 )
  JVS( 333) = W( 73 )
  JVS( 334) = W( 74 )
  JVS( 335) = W( 75 )
  JVS( 336) = W( 76 )
  JVS( 337) = W( 77 )
  JVS( 338) = W( 79 )
  IF ( ABS(  JVS( 347 )) < TINY(a) ) THEN
         IER = 51                                      
         RETURN
  END IF
   W( 22 ) = JVS( 339 )
   W( 25 ) = JVS( 340 )
   W( 34 ) = JVS( 341 )
   W( 39 ) = JVS( 342 )
   W( 42 ) = JVS( 343 )
   W( 43 ) = JVS( 344 )
   W( 48 ) = JVS( 345 )
   W( 50 ) = JVS( 346 )
   W( 51 ) = JVS( 347 )
   W( 52 ) = JVS( 348 )
   W( 54 ) = JVS( 349 )
   W( 57 ) = JVS( 350 )
   W( 58 ) = JVS( 351 )
   W( 61 ) = JVS( 352 )
   W( 69 ) = JVS( 353 )
   W( 71 ) = JVS( 354 )
   W( 72 ) = JVS( 355 )
   W( 73 ) = JVS( 356 )
   W( 74 ) = JVS( 357 )
   W( 75 ) = JVS( 358 )
   W( 76 ) = JVS( 359 )
   W( 77 ) = JVS( 360 )
   W( 78 ) = JVS( 361 )
   W( 79 ) = JVS( 362 )
   W( 80 ) = JVS( 363 )
  a = -W( 22 ) / JVS(           94  )
  W( 22 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 95 )
  W( 42 ) = W( 42 ) + a*JVS( 96 )
  W( 73 ) = W( 73 ) + a*JVS( 97 )
  W( 79 ) = W( 79 ) + a*JVS( 98 )
  W( 80 ) = W( 80 ) + a*JVS( 99 )
  a = -W( 25 ) / JVS(          118  )
  W( 25 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 119 )
  W( 69 ) = W( 69 ) + a*JVS( 120 )
  W( 73 ) = W( 73 ) + a*JVS( 121 )
  W( 76 ) = W( 76 ) + a*JVS( 122 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 173 )
  W( 50 ) = W( 50 ) + a*JVS( 174 )
  W( 73 ) = W( 73 ) + a*JVS( 175 )
  W( 74 ) = W( 74 ) + a*JVS( 176 )
  W( 79 ) = W( 79 ) + a*JVS( 177 )
  a = -W( 39 ) / JVS(          217  )
  W( 39 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 218 )
  W( 74 ) = W( 74 ) + a*JVS( 219 )
  W( 77 ) = W( 77 ) + a*JVS( 220 )
  a = -W( 42 ) / JVS(          238  )
  W( 42 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 239 )
  W( 73 ) = W( 73 ) + a*JVS( 240 )
  W( 74 ) = W( 74 ) + a*JVS( 241 )
  W( 75 ) = W( 75 ) + a*JVS( 242 )
  W( 77 ) = W( 77 ) + a*JVS( 243 )
  W( 79 ) = W( 79 ) + a*JVS( 244 )
  W( 80 ) = W( 80 ) + a*JVS( 245 )
  a = -W( 43 ) / JVS(          246  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 247 )
  W( 74 ) = W( 74 ) + a*JVS( 248 )
  W( 77 ) = W( 77 ) + a*JVS( 249 )
  W( 78 ) = W( 78 ) + a*JVS( 250 )
  W( 80 ) = W( 80 ) + a*JVS( 251 )
  a = -W( 48 ) / JVS(          294  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 295 )
  W( 54 ) = W( 54 ) + a*JVS( 296 )
  W( 72 ) = W( 72 ) + a*JVS( 297 )
  W( 73 ) = W( 73 ) + a*JVS( 298 )
  W( 74 ) = W( 74 ) + a*JVS( 299 )
  W( 75 ) = W( 75 ) + a*JVS( 300 )
  W( 76 ) = W( 76 ) + a*JVS( 301 )
  W( 77 ) = W( 77 ) + a*JVS( 302 )
  W( 79 ) = W( 79 ) + a*JVS( 303 )
  a = -W( 50 ) / JVS(          331  )
  W( 50 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 332 )
  W( 73 ) = W( 73 ) + a*JVS( 333 )
  W( 74 ) = W( 74 ) + a*JVS( 334 )
  W( 75 ) = W( 75 ) + a*JVS( 335 )
  W( 76 ) = W( 76 ) + a*JVS( 336 )
  W( 77 ) = W( 77 ) + a*JVS( 337 )
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  JVS( 339) = W( 22 )
  JVS( 340) = W( 25 )
  JVS( 341) = W( 34 )
  JVS( 342) = W( 39 )
  JVS( 343) = W( 42 )
  JVS( 344) = W( 43 )
  JVS( 345) = W( 48 )
  JVS( 346) = W( 50 )
  JVS( 347) = W( 51 )
  JVS( 348) = W( 52 )
  JVS( 349) = W( 54 )
  JVS( 350) = W( 57 )
  JVS( 351) = W( 58 )
  JVS( 352) = W( 61 )
  JVS( 353) = W( 69 )
  JVS( 354) = W( 71 )
  JVS( 355) = W( 72 )
  JVS( 356) = W( 73 )
  JVS( 357) = W( 74 )
  JVS( 358) = W( 75 )
  JVS( 359) = W( 76 )
  JVS( 360) = W( 77 )
  JVS( 361) = W( 78 )
  JVS( 362) = W( 79 )
  JVS( 363) = W( 80 )
  IF ( ABS(  JVS( 366 )) < TINY(a) ) THEN
         IER = 52                                      
         RETURN
  END IF
   W( 31 ) = JVS( 364 )
   W( 40 ) = JVS( 365 )
   W( 52 ) = JVS( 366 )
   W( 72 ) = JVS( 367 )
   W( 73 ) = JVS( 368 )
   W( 74 ) = JVS( 369 )
   W( 75 ) = JVS( 370 )
   W( 76 ) = JVS( 371 )
   W( 77 ) = JVS( 372 )
   W( 79 ) = JVS( 373 )
   W( 80 ) = JVS( 374 )
  a = -W( 31 ) / JVS(          150  )
  W( 31 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 151 )
  W( 73 ) = W( 73 ) + a*JVS( 152 )
  W( 77 ) = W( 77 ) + a*JVS( 153 )
  W( 80 ) = W( 80 ) + a*JVS( 154 )
  a = -W( 40 ) / JVS(          225  )
  W( 40 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 226 )
  W( 74 ) = W( 74 ) + a*JVS( 227 )
  W( 75 ) = W( 75 ) + a*JVS( 228 )
  W( 77 ) = W( 77 ) + a*JVS( 229 )
  W( 80 ) = W( 80 ) + a*JVS( 230 )
  JVS( 364) = W( 31 )
  JVS( 365) = W( 40 )
  JVS( 366) = W( 52 )
  JVS( 367) = W( 72 )
  JVS( 368) = W( 73 )
  JVS( 369) = W( 74 )
  JVS( 370) = W( 75 )
  JVS( 371) = W( 76 )
  JVS( 372) = W( 77 )
  JVS( 373) = W( 79 )
  JVS( 374) = W( 80 )
  IF ( ABS(  JVS( 376 )) < TINY(a) ) THEN
         IER = 53                                      
         RETURN
  END IF
   W( 12 ) = JVS( 375 )
   W( 53 ) = JVS( 376 )
   W( 72 ) = JVS( 377 )
   W( 73 ) = JVS( 378 )
   W( 74 ) = JVS( 379 )
   W( 75 ) = JVS( 380 )
   W( 76 ) = JVS( 381 )
   W( 79 ) = JVS( 382 )
  a = -W( 12 ) / JVS(           57  )
  W( 12 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 58 )
  JVS( 375) = W( 12 )
  JVS( 376) = W( 53 )
  JVS( 377) = W( 72 )
  JVS( 378) = W( 73 )
  JVS( 379) = W( 74 )
  JVS( 380) = W( 75 )
  JVS( 381) = W( 76 )
  JVS( 382) = W( 79 )
  IF ( ABS(  JVS( 384 )) < TINY(a) ) THEN
         IER = 54                                      
         RETURN
  END IF
   W( 38 ) = JVS( 383 )
   W( 54 ) = JVS( 384 )
   W( 72 ) = JVS( 385 )
   W( 73 ) = JVS( 386 )
   W( 74 ) = JVS( 387 )
   W( 75 ) = JVS( 388 )
   W( 76 ) = JVS( 389 )
   W( 77 ) = JVS( 390 )
   W( 79 ) = JVS( 391 )
  a = -W( 38 ) / JVS(          213  )
  W( 38 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 214 )
  W( 74 ) = W( 74 ) + a*JVS( 215 )
  W( 77 ) = W( 77 ) + a*JVS( 216 )
  JVS( 383) = W( 38 )
  JVS( 384) = W( 54 )
  JVS( 385) = W( 72 )
  JVS( 386) = W( 73 )
  JVS( 387) = W( 74 )
  JVS( 388) = W( 75 )
  JVS( 389) = W( 76 )
  JVS( 390) = W( 77 )
  JVS( 391) = W( 79 )
  IF ( ABS(  JVS( 393 )) < TINY(a) ) THEN
         IER = 55                                      
         RETURN
  END IF
   W( 8 ) = JVS( 392 )
   W( 55 ) = JVS( 393 )
   W( 72 ) = JVS( 394 )
   W( 73 ) = JVS( 395 )
   W( 74 ) = JVS( 396 )
   W( 75 ) = JVS( 397 )
   W( 76 ) = JVS( 398 )
   W( 79 ) = JVS( 399 )
  a = -W( 8 ) / JVS(           48  )
  W( 8 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 49 )
  JVS( 392) = W( 8 )
  JVS( 393) = W( 55 )
  JVS( 394) = W( 72 )
  JVS( 395) = W( 73 )
  JVS( 396) = W( 74 )
  JVS( 397) = W( 75 )
  JVS( 398) = W( 76 )
  JVS( 399) = W( 79 )
  IF ( ABS(  JVS( 423 )) < TINY(a) ) THEN
         IER = 56                                      
         RETURN
  END IF
   W( 14 ) = JVS( 400 )
   W( 18 ) = JVS( 401 )
   W( 20 ) = JVS( 402 )
   W( 25 ) = JVS( 403 )
   W( 30 ) = JVS( 404 )
   W( 33 ) = JVS( 405 )
   W( 34 ) = JVS( 406 )
   W( 35 ) = JVS( 407 )
   W( 38 ) = JVS( 408 )
   W( 39 ) = JVS( 409 )
   W( 41 ) = JVS( 410 )
   W( 42 ) = JVS( 411 )
   W( 43 ) = JVS( 412 )
   W( 44 ) = JVS( 413 )
   W( 45 ) = JVS( 414 )
   W( 46 ) = JVS( 415 )
   W( 47 ) = JVS( 416 )
   W( 48 ) = JVS( 417 )
   W( 50 ) = JVS( 418 )
   W( 52 ) = JVS( 419 )
   W( 53 ) = JVS( 420 )
   W( 54 ) = JVS( 421 )
   W( 55 ) = JVS( 422 )
   W( 56 ) = JVS( 423 )
   W( 57 ) = JVS( 424 )
   W( 58 ) = JVS( 425 )
   W( 59 ) = JVS( 426 )
   W( 61 ) = JVS( 427 )
   W( 62 ) = JVS( 428 )
   W( 63 ) = JVS( 429 )
   W( 65 ) = JVS( 430 )
   W( 66 ) = JVS( 431 )
   W( 67 ) = JVS( 432 )
   W( 68 ) = JVS( 433 )
   W( 69 ) = JVS( 434 )
   W( 71 ) = JVS( 435 )
   W( 72 ) = JVS( 436 )
   W( 73 ) = JVS( 437 )
   W( 74 ) = JVS( 438 )
   W( 75 ) = JVS( 439 )
   W( 76 ) = JVS( 440 )
   W( 77 ) = JVS( 441 )
   W( 78 ) = JVS( 442 )
   W( 79 ) = JVS( 443 )
   W( 80 ) = JVS( 444 )
  a = -W( 14 ) / JVS(           63  )
  W( 14 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 64 )
  a = -W( 18 ) / JVS(           78  )
  W( 18 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 79 )
  W( 73 ) = W( 73 ) + a*JVS( 80 )
  a = -W( 20 ) / JVS(           85  )
  W( 20 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 86 )
  W( 75 ) = W( 75 ) + a*JVS( 87 )
  W( 76 ) = W( 76 ) + a*JVS( 88 )
  a = -W( 25 ) / JVS(          118  )
  W( 25 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 119 )
  W( 69 ) = W( 69 ) + a*JVS( 120 )
  W( 73 ) = W( 73 ) + a*JVS( 121 )
  W( 76 ) = W( 76 ) + a*JVS( 122 )
  a = -W( 30 ) / JVS(          146  )
  W( 30 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 147 )
  W( 74 ) = W( 74 ) + a*JVS( 148 )
  W( 77 ) = W( 77 ) + a*JVS( 149 )
  a = -W( 33 ) / JVS(          166  )
  W( 33 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 167 )
  W( 72 ) = W( 72 ) + a*JVS( 168 )
  W( 73 ) = W( 73 ) + a*JVS( 169 )
  W( 75 ) = W( 75 ) + a*JVS( 170 )
  W( 77 ) = W( 77 ) + a*JVS( 171 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 173 )
  W( 50 ) = W( 50 ) + a*JVS( 174 )
  W( 73 ) = W( 73 ) + a*JVS( 175 )
  W( 74 ) = W( 74 ) + a*JVS( 176 )
  W( 79 ) = W( 79 ) + a*JVS( 177 )
  a = -W( 35 ) / JVS(          178  )
  W( 35 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 179 )
  W( 72 ) = W( 72 ) + a*JVS( 180 )
  W( 73 ) = W( 73 ) + a*JVS( 181 )
  W( 74 ) = W( 74 ) + a*JVS( 182 )
  W( 77 ) = W( 77 ) + a*JVS( 183 )
  W( 80 ) = W( 80 ) + a*JVS( 184 )
  a = -W( 38 ) / JVS(          213  )
  W( 38 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 214 )
  W( 74 ) = W( 74 ) + a*JVS( 215 )
  W( 77 ) = W( 77 ) + a*JVS( 216 )
  a = -W( 39 ) / JVS(          217  )
  W( 39 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 218 )
  W( 74 ) = W( 74 ) + a*JVS( 219 )
  W( 77 ) = W( 77 ) + a*JVS( 220 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  a = -W( 42 ) / JVS(          238  )
  W( 42 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 239 )
  W( 73 ) = W( 73 ) + a*JVS( 240 )
  W( 74 ) = W( 74 ) + a*JVS( 241 )
  W( 75 ) = W( 75 ) + a*JVS( 242 )
  W( 77 ) = W( 77 ) + a*JVS( 243 )
  W( 79 ) = W( 79 ) + a*JVS( 244 )
  W( 80 ) = W( 80 ) + a*JVS( 245 )
  a = -W( 43 ) / JVS(          246  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 247 )
  W( 74 ) = W( 74 ) + a*JVS( 248 )
  W( 77 ) = W( 77 ) + a*JVS( 249 )
  W( 78 ) = W( 78 ) + a*JVS( 250 )
  W( 80 ) = W( 80 ) + a*JVS( 251 )
  a = -W( 44 ) / JVS(          253  )
  W( 44 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 254 )
  W( 57 ) = W( 57 ) + a*JVS( 255 )
  W( 58 ) = W( 58 ) + a*JVS( 256 )
  W( 61 ) = W( 61 ) + a*JVS( 257 )
  W( 71 ) = W( 71 ) + a*JVS( 258 )
  W( 72 ) = W( 72 ) + a*JVS( 259 )
  W( 73 ) = W( 73 ) + a*JVS( 260 )
  W( 74 ) = W( 74 ) + a*JVS( 261 )
  W( 76 ) = W( 76 ) + a*JVS( 262 )
  W( 77 ) = W( 77 ) + a*JVS( 263 )
  W( 79 ) = W( 79 ) + a*JVS( 264 )
  a = -W( 45 ) / JVS(          266  )
  W( 45 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 267 )
  W( 73 ) = W( 73 ) + a*JVS( 268 )
  W( 74 ) = W( 74 ) + a*JVS( 269 )
  W( 75 ) = W( 75 ) + a*JVS( 270 )
  W( 76 ) = W( 76 ) + a*JVS( 271 )
  W( 77 ) = W( 77 ) + a*JVS( 272 )
  W( 79 ) = W( 79 ) + a*JVS( 273 )
  a = -W( 46 ) / JVS(          274  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 275 )
  W( 72 ) = W( 72 ) + a*JVS( 276 )
  W( 73 ) = W( 73 ) + a*JVS( 277 )
  W( 74 ) = W( 74 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 76 ) = W( 76 ) + a*JVS( 280 )
  W( 79 ) = W( 79 ) + a*JVS( 281 )
  a = -W( 47 ) / JVS(          282  )
  W( 47 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 283 )
  W( 72 ) = W( 72 ) + a*JVS( 284 )
  W( 73 ) = W( 73 ) + a*JVS( 285 )
  W( 74 ) = W( 74 ) + a*JVS( 286 )
  W( 75 ) = W( 75 ) + a*JVS( 287 )
  W( 76 ) = W( 76 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  a = -W( 48 ) / JVS(          294  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 295 )
  W( 54 ) = W( 54 ) + a*JVS( 296 )
  W( 72 ) = W( 72 ) + a*JVS( 297 )
  W( 73 ) = W( 73 ) + a*JVS( 298 )
  W( 74 ) = W( 74 ) + a*JVS( 299 )
  W( 75 ) = W( 75 ) + a*JVS( 300 )
  W( 76 ) = W( 76 ) + a*JVS( 301 )
  W( 77 ) = W( 77 ) + a*JVS( 302 )
  W( 79 ) = W( 79 ) + a*JVS( 303 )
  a = -W( 50 ) / JVS(          331  )
  W( 50 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 332 )
  W( 73 ) = W( 73 ) + a*JVS( 333 )
  W( 74 ) = W( 74 ) + a*JVS( 334 )
  W( 75 ) = W( 75 ) + a*JVS( 335 )
  W( 76 ) = W( 76 ) + a*JVS( 336 )
  W( 77 ) = W( 77 ) + a*JVS( 337 )
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  a = -W( 52 ) / JVS(          366  )
  W( 52 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 367 )
  W( 73 ) = W( 73 ) + a*JVS( 368 )
  W( 74 ) = W( 74 ) + a*JVS( 369 )
  W( 75 ) = W( 75 ) + a*JVS( 370 )
  W( 76 ) = W( 76 ) + a*JVS( 371 )
  W( 77 ) = W( 77 ) + a*JVS( 372 )
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 80 ) = W( 80 ) + a*JVS( 374 )
  a = -W( 53 ) / JVS(          376  )
  W( 53 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 377 )
  W( 73 ) = W( 73 ) + a*JVS( 378 )
  W( 74 ) = W( 74 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 79 ) = W( 79 ) + a*JVS( 382 )
  a = -W( 54 ) / JVS(          384  )
  W( 54 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 385 )
  W( 73 ) = W( 73 ) + a*JVS( 386 )
  W( 74 ) = W( 74 ) + a*JVS( 387 )
  W( 75 ) = W( 75 ) + a*JVS( 388 )
  W( 76 ) = W( 76 ) + a*JVS( 389 )
  W( 77 ) = W( 77 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 55 ) / JVS(          393  )
  W( 55 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  W( 73 ) = W( 73 ) + a*JVS( 395 )
  W( 74 ) = W( 74 ) + a*JVS( 396 )
  W( 75 ) = W( 75 ) + a*JVS( 397 )
  W( 76 ) = W( 76 ) + a*JVS( 398 )
  W( 79 ) = W( 79 ) + a*JVS( 399 )
  JVS( 400) = W( 14 )
  JVS( 401) = W( 18 )
  JVS( 402) = W( 20 )
  JVS( 403) = W( 25 )
  JVS( 404) = W( 30 )
  JVS( 405) = W( 33 )
  JVS( 406) = W( 34 )
  JVS( 407) = W( 35 )
  JVS( 408) = W( 38 )
  JVS( 409) = W( 39 )
  JVS( 410) = W( 41 )
  JVS( 411) = W( 42 )
  JVS( 412) = W( 43 )
  JVS( 413) = W( 44 )
  JVS( 414) = W( 45 )
  JVS( 415) = W( 46 )
  JVS( 416) = W( 47 )
  JVS( 417) = W( 48 )
  JVS( 418) = W( 50 )
  JVS( 419) = W( 52 )
  JVS( 420) = W( 53 )
  JVS( 421) = W( 54 )
  JVS( 422) = W( 55 )
  JVS( 423) = W( 56 )
  JVS( 424) = W( 57 )
  JVS( 425) = W( 58 )
  JVS( 426) = W( 59 )
  JVS( 427) = W( 61 )
  JVS( 428) = W( 62 )
  JVS( 429) = W( 63 )
  JVS( 430) = W( 65 )
  JVS( 431) = W( 66 )
  JVS( 432) = W( 67 )
  JVS( 433) = W( 68 )
  JVS( 434) = W( 69 )
  JVS( 435) = W( 71 )
  JVS( 436) = W( 72 )
  JVS( 437) = W( 73 )
  JVS( 438) = W( 74 )
  JVS( 439) = W( 75 )
  JVS( 440) = W( 76 )
  JVS( 441) = W( 77 )
  JVS( 442) = W( 78 )
  JVS( 443) = W( 79 )
  JVS( 444) = W( 80 )
  IF ( ABS(  JVS( 446 )) < TINY(a) ) THEN
         IER = 57                                      
         RETURN
  END IF
   W( 28 ) = JVS( 445 )
   W( 57 ) = JVS( 446 )
   W( 72 ) = JVS( 447 )
   W( 73 ) = JVS( 448 )
   W( 74 ) = JVS( 449 )
   W( 75 ) = JVS( 450 )
   W( 76 ) = JVS( 451 )
   W( 77 ) = JVS( 452 )
   W( 79 ) = JVS( 453 )
   W( 80 ) = JVS( 454 )
  a = -W( 28 ) / JVS(          137  )
  W( 28 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 138 )
  W( 77 ) = W( 77 ) + a*JVS( 139 )
  W( 80 ) = W( 80 ) + a*JVS( 140 )
  JVS( 445) = W( 28 )
  JVS( 446) = W( 57 )
  JVS( 447) = W( 72 )
  JVS( 448) = W( 73 )
  JVS( 449) = W( 74 )
  JVS( 450) = W( 75 )
  JVS( 451) = W( 76 )
  JVS( 452) = W( 77 )
  JVS( 453) = W( 79 )
  JVS( 454) = W( 80 )
  IF ( ABS(  JVS( 456 )) < TINY(a) ) THEN
         IER = 58                                      
         RETURN
  END IF
   W( 29 ) = JVS( 455 )
   W( 58 ) = JVS( 456 )
   W( 72 ) = JVS( 457 )
   W( 73 ) = JVS( 458 )
   W( 74 ) = JVS( 459 )
   W( 75 ) = JVS( 460 )
   W( 76 ) = JVS( 461 )
   W( 77 ) = JVS( 462 )
   W( 79 ) = JVS( 463 )
   W( 80 ) = JVS( 464 )
  a = -W( 29 ) / JVS(          142  )
  W( 29 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 143 )
  W( 77 ) = W( 77 ) + a*JVS( 144 )
  W( 80 ) = W( 80 ) + a*JVS( 145 )
  JVS( 455) = W( 29 )
  JVS( 456) = W( 58 )
  JVS( 457) = W( 72 )
  JVS( 458) = W( 73 )
  JVS( 459) = W( 74 )
  JVS( 460) = W( 75 )
  JVS( 461) = W( 76 )
  JVS( 462) = W( 77 )
  JVS( 463) = W( 79 )
  JVS( 464) = W( 80 )
  IF ( ABS(  JVS( 466 )) < TINY(a) ) THEN
         IER = 59                                      
         RETURN
  END IF
   W( 36 ) = JVS( 465 )
   W( 59 ) = JVS( 466 )
   W( 72 ) = JVS( 467 )
   W( 73 ) = JVS( 468 )
   W( 74 ) = JVS( 469 )
   W( 75 ) = JVS( 470 )
   W( 76 ) = JVS( 471 )
   W( 77 ) = JVS( 472 )
   W( 79 ) = JVS( 473 )
  a = -W( 36 ) / JVS(          185  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 186 )
  W( 74 ) = W( 74 ) + a*JVS( 187 )
  W( 77 ) = W( 77 ) + a*JVS( 188 )
  JVS( 465) = W( 36 )
  JVS( 466) = W( 59 )
  JVS( 467) = W( 72 )
  JVS( 468) = W( 73 )
  JVS( 469) = W( 74 )
  JVS( 470) = W( 75 )
  JVS( 471) = W( 76 )
  JVS( 472) = W( 77 )
  JVS( 473) = W( 79 )
  IF ( ABS(  JVS( 483 )) < TINY(a) ) THEN
         IER = 60                                      
         RETURN
  END IF
   W( 27 ) = JVS( 474 )
   W( 40 ) = JVS( 475 )
   W( 43 ) = JVS( 476 )
   W( 53 ) = JVS( 477 )
   W( 54 ) = JVS( 478 )
   W( 55 ) = JVS( 479 )
   W( 57 ) = JVS( 480 )
   W( 58 ) = JVS( 481 )
   W( 59 ) = JVS( 482 )
   W( 60 ) = JVS( 483 )
   W( 66 ) = JVS( 484 )
   W( 67 ) = JVS( 485 )
   W( 71 ) = JVS( 486 )
   W( 72 ) = JVS( 487 )
   W( 73 ) = JVS( 488 )
   W( 74 ) = JVS( 489 )
   W( 75 ) = JVS( 490 )
   W( 76 ) = JVS( 491 )
   W( 77 ) = JVS( 492 )
   W( 78 ) = JVS( 493 )
   W( 79 ) = JVS( 494 )
   W( 80 ) = JVS( 495 )
  a = -W( 27 ) / JVS(          130  )
  W( 27 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 131 )
  W( 73 ) = W( 73 ) + a*JVS( 132 )
  W( 74 ) = W( 74 ) + a*JVS( 133 )
  W( 75 ) = W( 75 ) + a*JVS( 134 )
  W( 80 ) = W( 80 ) + a*JVS( 135 )
  a = -W( 40 ) / JVS(          225  )
  W( 40 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 226 )
  W( 74 ) = W( 74 ) + a*JVS( 227 )
  W( 75 ) = W( 75 ) + a*JVS( 228 )
  W( 77 ) = W( 77 ) + a*JVS( 229 )
  W( 80 ) = W( 80 ) + a*JVS( 230 )
  a = -W( 43 ) / JVS(          246  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 247 )
  W( 74 ) = W( 74 ) + a*JVS( 248 )
  W( 77 ) = W( 77 ) + a*JVS( 249 )
  W( 78 ) = W( 78 ) + a*JVS( 250 )
  W( 80 ) = W( 80 ) + a*JVS( 251 )
  a = -W( 53 ) / JVS(          376  )
  W( 53 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 377 )
  W( 73 ) = W( 73 ) + a*JVS( 378 )
  W( 74 ) = W( 74 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 79 ) = W( 79 ) + a*JVS( 382 )
  a = -W( 54 ) / JVS(          384  )
  W( 54 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 385 )
  W( 73 ) = W( 73 ) + a*JVS( 386 )
  W( 74 ) = W( 74 ) + a*JVS( 387 )
  W( 75 ) = W( 75 ) + a*JVS( 388 )
  W( 76 ) = W( 76 ) + a*JVS( 389 )
  W( 77 ) = W( 77 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 55 ) / JVS(          393  )
  W( 55 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  W( 73 ) = W( 73 ) + a*JVS( 395 )
  W( 74 ) = W( 74 ) + a*JVS( 396 )
  W( 75 ) = W( 75 ) + a*JVS( 397 )
  W( 76 ) = W( 76 ) + a*JVS( 398 )
  W( 79 ) = W( 79 ) + a*JVS( 399 )
  a = -W( 57 ) / JVS(          446  )
  W( 57 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 447 )
  W( 73 ) = W( 73 ) + a*JVS( 448 )
  W( 74 ) = W( 74 ) + a*JVS( 449 )
  W( 75 ) = W( 75 ) + a*JVS( 450 )
  W( 76 ) = W( 76 ) + a*JVS( 451 )
  W( 77 ) = W( 77 ) + a*JVS( 452 )
  W( 79 ) = W( 79 ) + a*JVS( 453 )
  W( 80 ) = W( 80 ) + a*JVS( 454 )
  a = -W( 58 ) / JVS(          456  )
  W( 58 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 457 )
  W( 73 ) = W( 73 ) + a*JVS( 458 )
  W( 74 ) = W( 74 ) + a*JVS( 459 )
  W( 75 ) = W( 75 ) + a*JVS( 460 )
  W( 76 ) = W( 76 ) + a*JVS( 461 )
  W( 77 ) = W( 77 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  a = -W( 59 ) / JVS(          466  )
  W( 59 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 467 )
  W( 73 ) = W( 73 ) + a*JVS( 468 )
  W( 74 ) = W( 74 ) + a*JVS( 469 )
  W( 75 ) = W( 75 ) + a*JVS( 470 )
  W( 76 ) = W( 76 ) + a*JVS( 471 )
  W( 77 ) = W( 77 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  JVS( 474) = W( 27 )
  JVS( 475) = W( 40 )
  JVS( 476) = W( 43 )
  JVS( 477) = W( 53 )
  JVS( 478) = W( 54 )
  JVS( 479) = W( 55 )
  JVS( 480) = W( 57 )
  JVS( 481) = W( 58 )
  JVS( 482) = W( 59 )
  JVS( 483) = W( 60 )
  JVS( 484) = W( 66 )
  JVS( 485) = W( 67 )
  JVS( 486) = W( 71 )
  JVS( 487) = W( 72 )
  JVS( 488) = W( 73 )
  JVS( 489) = W( 74 )
  JVS( 490) = W( 75 )
  JVS( 491) = W( 76 )
  JVS( 492) = W( 77 )
  JVS( 493) = W( 78 )
  JVS( 494) = W( 79 )
  JVS( 495) = W( 80 )
  IF ( ABS(  JVS( 498 )) < TINY(a) ) THEN
         IER = 61                                      
         RETURN
  END IF
   W( 57 ) = JVS( 496 )
   W( 58 ) = JVS( 497 )
   W( 61 ) = JVS( 498 )
   W( 72 ) = JVS( 499 )
   W( 73 ) = JVS( 500 )
   W( 74 ) = JVS( 501 )
   W( 75 ) = JVS( 502 )
   W( 76 ) = JVS( 503 )
   W( 77 ) = JVS( 504 )
   W( 79 ) = JVS( 505 )
   W( 80 ) = JVS( 506 )
  a = -W( 57 ) / JVS(          446  )
  W( 57 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 447 )
  W( 73 ) = W( 73 ) + a*JVS( 448 )
  W( 74 ) = W( 74 ) + a*JVS( 449 )
  W( 75 ) = W( 75 ) + a*JVS( 450 )
  W( 76 ) = W( 76 ) + a*JVS( 451 )
  W( 77 ) = W( 77 ) + a*JVS( 452 )
  W( 79 ) = W( 79 ) + a*JVS( 453 )
  W( 80 ) = W( 80 ) + a*JVS( 454 )
  a = -W( 58 ) / JVS(          456  )
  W( 58 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 457 )
  W( 73 ) = W( 73 ) + a*JVS( 458 )
  W( 74 ) = W( 74 ) + a*JVS( 459 )
  W( 75 ) = W( 75 ) + a*JVS( 460 )
  W( 76 ) = W( 76 ) + a*JVS( 461 )
  W( 77 ) = W( 77 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  JVS( 496) = W( 57 )
  JVS( 497) = W( 58 )
  JVS( 498) = W( 61 )
  JVS( 499) = W( 72 )
  JVS( 500) = W( 73 )
  JVS( 501) = W( 74 )
  JVS( 502) = W( 75 )
  JVS( 503) = W( 76 )
  JVS( 504) = W( 77 )
  JVS( 505) = W( 79 )
  JVS( 506) = W( 80 )
  IF ( ABS(  JVS( 517 )) < TINY(a) ) THEN
         IER = 62                                      
         RETURN
  END IF
   W( 9 ) = JVS( 507 )
   W( 10 ) = JVS( 508 )
   W( 33 ) = JVS( 509 )
   W( 35 ) = JVS( 510 )
   W( 40 ) = JVS( 511 )
   W( 41 ) = JVS( 512 )
   W( 43 ) = JVS( 513 )
   W( 53 ) = JVS( 514 )
   W( 55 ) = JVS( 515 )
   W( 61 ) = JVS( 516 )
   W( 62 ) = JVS( 517 )
   W( 69 ) = JVS( 518 )
   W( 70 ) = JVS( 519 )
   W( 71 ) = JVS( 520 )
   W( 72 ) = JVS( 521 )
   W( 73 ) = JVS( 522 )
   W( 74 ) = JVS( 523 )
   W( 75 ) = JVS( 524 )
   W( 76 ) = JVS( 525 )
   W( 77 ) = JVS( 526 )
   W( 78 ) = JVS( 527 )
   W( 79 ) = JVS( 528 )
   W( 80 ) = JVS( 529 )
  a = -W( 9 ) / JVS(           50  )
  W( 9 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 51 )
  a = -W( 10 ) / JVS(           52  )
  W( 10 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 53 )
  a = -W( 33 ) / JVS(          166  )
  W( 33 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 167 )
  W( 72 ) = W( 72 ) + a*JVS( 168 )
  W( 73 ) = W( 73 ) + a*JVS( 169 )
  W( 75 ) = W( 75 ) + a*JVS( 170 )
  W( 77 ) = W( 77 ) + a*JVS( 171 )
  a = -W( 35 ) / JVS(          178  )
  W( 35 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 179 )
  W( 72 ) = W( 72 ) + a*JVS( 180 )
  W( 73 ) = W( 73 ) + a*JVS( 181 )
  W( 74 ) = W( 74 ) + a*JVS( 182 )
  W( 77 ) = W( 77 ) + a*JVS( 183 )
  W( 80 ) = W( 80 ) + a*JVS( 184 )
  a = -W( 40 ) / JVS(          225  )
  W( 40 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 226 )
  W( 74 ) = W( 74 ) + a*JVS( 227 )
  W( 75 ) = W( 75 ) + a*JVS( 228 )
  W( 77 ) = W( 77 ) + a*JVS( 229 )
  W( 80 ) = W( 80 ) + a*JVS( 230 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  a = -W( 43 ) / JVS(          246  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 247 )
  W( 74 ) = W( 74 ) + a*JVS( 248 )
  W( 77 ) = W( 77 ) + a*JVS( 249 )
  W( 78 ) = W( 78 ) + a*JVS( 250 )
  W( 80 ) = W( 80 ) + a*JVS( 251 )
  a = -W( 53 ) / JVS(          376  )
  W( 53 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 377 )
  W( 73 ) = W( 73 ) + a*JVS( 378 )
  W( 74 ) = W( 74 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 79 ) = W( 79 ) + a*JVS( 382 )
  a = -W( 55 ) / JVS(          393  )
  W( 55 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  W( 73 ) = W( 73 ) + a*JVS( 395 )
  W( 74 ) = W( 74 ) + a*JVS( 396 )
  W( 75 ) = W( 75 ) + a*JVS( 397 )
  W( 76 ) = W( 76 ) + a*JVS( 398 )
  W( 79 ) = W( 79 ) + a*JVS( 399 )
  a = -W( 61 ) / JVS(          498  )
  W( 61 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 499 )
  W( 73 ) = W( 73 ) + a*JVS( 500 )
  W( 74 ) = W( 74 ) + a*JVS( 501 )
  W( 75 ) = W( 75 ) + a*JVS( 502 )
  W( 76 ) = W( 76 ) + a*JVS( 503 )
  W( 77 ) = W( 77 ) + a*JVS( 504 )
  W( 79 ) = W( 79 ) + a*JVS( 505 )
  W( 80 ) = W( 80 ) + a*JVS( 506 )
  JVS( 507) = W( 9 )
  JVS( 508) = W( 10 )
  JVS( 509) = W( 33 )
  JVS( 510) = W( 35 )
  JVS( 511) = W( 40 )
  JVS( 512) = W( 41 )
  JVS( 513) = W( 43 )
  JVS( 514) = W( 53 )
  JVS( 515) = W( 55 )
  JVS( 516) = W( 61 )
  JVS( 517) = W( 62 )
  JVS( 518) = W( 69 )
  JVS( 519) = W( 70 )
  JVS( 520) = W( 71 )
  JVS( 521) = W( 72 )
  JVS( 522) = W( 73 )
  JVS( 523) = W( 74 )
  JVS( 524) = W( 75 )
  JVS( 525) = W( 76 )
  JVS( 526) = W( 77 )
  JVS( 527) = W( 78 )
  JVS( 528) = W( 79 )
  JVS( 529) = W( 80 )
  IF ( ABS(  JVS( 533 )) < TINY(a) ) THEN
         IER = 63                                      
         RETURN
  END IF
   W( 38 ) = JVS( 530 )
   W( 41 ) = JVS( 531 )
   W( 50 ) = JVS( 532 )
   W( 63 ) = JVS( 533 )
   W( 72 ) = JVS( 534 )
   W( 73 ) = JVS( 535 )
   W( 74 ) = JVS( 536 )
   W( 75 ) = JVS( 537 )
   W( 76 ) = JVS( 538 )
   W( 77 ) = JVS( 539 )
   W( 79 ) = JVS( 540 )
  a = -W( 38 ) / JVS(          213  )
  W( 38 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 214 )
  W( 74 ) = W( 74 ) + a*JVS( 215 )
  W( 77 ) = W( 77 ) + a*JVS( 216 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  a = -W( 50 ) / JVS(          331  )
  W( 50 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 332 )
  W( 73 ) = W( 73 ) + a*JVS( 333 )
  W( 74 ) = W( 74 ) + a*JVS( 334 )
  W( 75 ) = W( 75 ) + a*JVS( 335 )
  W( 76 ) = W( 76 ) + a*JVS( 336 )
  W( 77 ) = W( 77 ) + a*JVS( 337 )
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  JVS( 530) = W( 38 )
  JVS( 531) = W( 41 )
  JVS( 532) = W( 50 )
  JVS( 533) = W( 63 )
  JVS( 534) = W( 72 )
  JVS( 535) = W( 73 )
  JVS( 536) = W( 74 )
  JVS( 537) = W( 75 )
  JVS( 538) = W( 76 )
  JVS( 539) = W( 77 )
  JVS( 540) = W( 79 )
  IF ( ABS(  JVS( 554 )) < TINY(a) ) THEN
         IER = 64                                      
         RETURN
  END IF
   W( 12 ) = JVS( 541 )
   W( 14 ) = JVS( 542 )
   W( 16 ) = JVS( 543 )
   W( 36 ) = JVS( 544 )
   W( 45 ) = JVS( 545 )
   W( 46 ) = JVS( 546 )
   W( 47 ) = JVS( 547 )
   W( 53 ) = JVS( 548 )
   W( 55 ) = JVS( 549 )
   W( 59 ) = JVS( 550 )
   W( 60 ) = JVS( 551 )
   W( 61 ) = JVS( 552 )
   W( 63 ) = JVS( 553 )
   W( 64 ) = JVS( 554 )
   W( 65 ) = JVS( 555 )
   W( 66 ) = JVS( 556 )
   W( 67 ) = JVS( 557 )
   W( 68 ) = JVS( 558 )
   W( 69 ) = JVS( 559 )
   W( 70 ) = JVS( 560 )
   W( 71 ) = JVS( 561 )
   W( 72 ) = JVS( 562 )
   W( 73 ) = JVS( 563 )
   W( 74 ) = JVS( 564 )
   W( 75 ) = JVS( 565 )
   W( 76 ) = JVS( 566 )
   W( 77 ) = JVS( 567 )
   W( 78 ) = JVS( 568 )
   W( 79 ) = JVS( 569 )
   W( 80 ) = JVS( 570 )
  a = -W( 12 ) / JVS(           57  )
  W( 12 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 58 )
  a = -W( 14 ) / JVS(           63  )
  W( 14 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 64 )
  a = -W( 16 ) / JVS(           70  )
  W( 16 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 71 )
  W( 73 ) = W( 73 ) + a*JVS( 72 )
  a = -W( 36 ) / JVS(          185  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 186 )
  W( 74 ) = W( 74 ) + a*JVS( 187 )
  W( 77 ) = W( 77 ) + a*JVS( 188 )
  a = -W( 45 ) / JVS(          266  )
  W( 45 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 267 )
  W( 73 ) = W( 73 ) + a*JVS( 268 )
  W( 74 ) = W( 74 ) + a*JVS( 269 )
  W( 75 ) = W( 75 ) + a*JVS( 270 )
  W( 76 ) = W( 76 ) + a*JVS( 271 )
  W( 77 ) = W( 77 ) + a*JVS( 272 )
  W( 79 ) = W( 79 ) + a*JVS( 273 )
  a = -W( 46 ) / JVS(          274  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 275 )
  W( 72 ) = W( 72 ) + a*JVS( 276 )
  W( 73 ) = W( 73 ) + a*JVS( 277 )
  W( 74 ) = W( 74 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 76 ) = W( 76 ) + a*JVS( 280 )
  W( 79 ) = W( 79 ) + a*JVS( 281 )
  a = -W( 47 ) / JVS(          282  )
  W( 47 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 283 )
  W( 72 ) = W( 72 ) + a*JVS( 284 )
  W( 73 ) = W( 73 ) + a*JVS( 285 )
  W( 74 ) = W( 74 ) + a*JVS( 286 )
  W( 75 ) = W( 75 ) + a*JVS( 287 )
  W( 76 ) = W( 76 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  a = -W( 53 ) / JVS(          376  )
  W( 53 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 377 )
  W( 73 ) = W( 73 ) + a*JVS( 378 )
  W( 74 ) = W( 74 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 79 ) = W( 79 ) + a*JVS( 382 )
  a = -W( 55 ) / JVS(          393  )
  W( 55 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  W( 73 ) = W( 73 ) + a*JVS( 395 )
  W( 74 ) = W( 74 ) + a*JVS( 396 )
  W( 75 ) = W( 75 ) + a*JVS( 397 )
  W( 76 ) = W( 76 ) + a*JVS( 398 )
  W( 79 ) = W( 79 ) + a*JVS( 399 )
  a = -W( 59 ) / JVS(          466  )
  W( 59 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 467 )
  W( 73 ) = W( 73 ) + a*JVS( 468 )
  W( 74 ) = W( 74 ) + a*JVS( 469 )
  W( 75 ) = W( 75 ) + a*JVS( 470 )
  W( 76 ) = W( 76 ) + a*JVS( 471 )
  W( 77 ) = W( 77 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  a = -W( 60 ) / JVS(          483  )
  W( 60 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 484 )
  W( 67 ) = W( 67 ) + a*JVS( 485 )
  W( 71 ) = W( 71 ) + a*JVS( 486 )
  W( 72 ) = W( 72 ) + a*JVS( 487 )
  W( 73 ) = W( 73 ) + a*JVS( 488 )
  W( 74 ) = W( 74 ) + a*JVS( 489 )
  W( 75 ) = W( 75 ) + a*JVS( 490 )
  W( 76 ) = W( 76 ) + a*JVS( 491 )
  W( 77 ) = W( 77 ) + a*JVS( 492 )
  W( 78 ) = W( 78 ) + a*JVS( 493 )
  W( 79 ) = W( 79 ) + a*JVS( 494 )
  W( 80 ) = W( 80 ) + a*JVS( 495 )
  a = -W( 61 ) / JVS(          498  )
  W( 61 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 499 )
  W( 73 ) = W( 73 ) + a*JVS( 500 )
  W( 74 ) = W( 74 ) + a*JVS( 501 )
  W( 75 ) = W( 75 ) + a*JVS( 502 )
  W( 76 ) = W( 76 ) + a*JVS( 503 )
  W( 77 ) = W( 77 ) + a*JVS( 504 )
  W( 79 ) = W( 79 ) + a*JVS( 505 )
  W( 80 ) = W( 80 ) + a*JVS( 506 )
  a = -W( 63 ) / JVS(          533  )
  W( 63 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 534 )
  W( 73 ) = W( 73 ) + a*JVS( 535 )
  W( 74 ) = W( 74 ) + a*JVS( 536 )
  W( 75 ) = W( 75 ) + a*JVS( 537 )
  W( 76 ) = W( 76 ) + a*JVS( 538 )
  W( 77 ) = W( 77 ) + a*JVS( 539 )
  W( 79 ) = W( 79 ) + a*JVS( 540 )
  JVS( 541) = W( 12 )
  JVS( 542) = W( 14 )
  JVS( 543) = W( 16 )
  JVS( 544) = W( 36 )
  JVS( 545) = W( 45 )
  JVS( 546) = W( 46 )
  JVS( 547) = W( 47 )
  JVS( 548) = W( 53 )
  JVS( 549) = W( 55 )
  JVS( 550) = W( 59 )
  JVS( 551) = W( 60 )
  JVS( 552) = W( 61 )
  JVS( 553) = W( 63 )
  JVS( 554) = W( 64 )
  JVS( 555) = W( 65 )
  JVS( 556) = W( 66 )
  JVS( 557) = W( 67 )
  JVS( 558) = W( 68 )
  JVS( 559) = W( 69 )
  JVS( 560) = W( 70 )
  JVS( 561) = W( 71 )
  JVS( 562) = W( 72 )
  JVS( 563) = W( 73 )
  JVS( 564) = W( 74 )
  JVS( 565) = W( 75 )
  JVS( 566) = W( 76 )
  JVS( 567) = W( 77 )
  JVS( 568) = W( 78 )
  JVS( 569) = W( 79 )
  JVS( 570) = W( 80 )
  IF ( ABS(  JVS( 573 )) < TINY(a) ) THEN
         IER = 65                                      
         RETURN
  END IF
   W( 50 ) = JVS( 571 )
   W( 54 ) = JVS( 572 )
   W( 65 ) = JVS( 573 )
   W( 72 ) = JVS( 574 )
   W( 73 ) = JVS( 575 )
   W( 74 ) = JVS( 576 )
   W( 75 ) = JVS( 577 )
   W( 76 ) = JVS( 578 )
   W( 77 ) = JVS( 579 )
   W( 79 ) = JVS( 580 )
  a = -W( 50 ) / JVS(          331  )
  W( 50 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 332 )
  W( 73 ) = W( 73 ) + a*JVS( 333 )
  W( 74 ) = W( 74 ) + a*JVS( 334 )
  W( 75 ) = W( 75 ) + a*JVS( 335 )
  W( 76 ) = W( 76 ) + a*JVS( 336 )
  W( 77 ) = W( 77 ) + a*JVS( 337 )
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  a = -W( 54 ) / JVS(          384  )
  W( 54 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 385 )
  W( 73 ) = W( 73 ) + a*JVS( 386 )
  W( 74 ) = W( 74 ) + a*JVS( 387 )
  W( 75 ) = W( 75 ) + a*JVS( 388 )
  W( 76 ) = W( 76 ) + a*JVS( 389 )
  W( 77 ) = W( 77 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  JVS( 571) = W( 50 )
  JVS( 572) = W( 54 )
  JVS( 573) = W( 65 )
  JVS( 574) = W( 72 )
  JVS( 575) = W( 73 )
  JVS( 576) = W( 74 )
  JVS( 577) = W( 75 )
  JVS( 578) = W( 76 )
  JVS( 579) = W( 77 )
  JVS( 580) = W( 79 )
  IF ( ABS(  JVS( 587 )) < TINY(a) ) THEN
         IER = 66                                      
         RETURN
  END IF
   W( 30 ) = JVS( 581 )
   W( 36 ) = JVS( 582 )
   W( 38 ) = JVS( 583 )
   W( 41 ) = JVS( 584 )
   W( 63 ) = JVS( 585 )
   W( 65 ) = JVS( 586 )
   W( 66 ) = JVS( 587 )
   W( 67 ) = JVS( 588 )
   W( 72 ) = JVS( 589 )
   W( 73 ) = JVS( 590 )
   W( 74 ) = JVS( 591 )
   W( 75 ) = JVS( 592 )
   W( 76 ) = JVS( 593 )
   W( 77 ) = JVS( 594 )
   W( 79 ) = JVS( 595 )
  a = -W( 30 ) / JVS(          146  )
  W( 30 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 147 )
  W( 74 ) = W( 74 ) + a*JVS( 148 )
  W( 77 ) = W( 77 ) + a*JVS( 149 )
  a = -W( 36 ) / JVS(          185  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 186 )
  W( 74 ) = W( 74 ) + a*JVS( 187 )
  W( 77 ) = W( 77 ) + a*JVS( 188 )
  a = -W( 38 ) / JVS(          213  )
  W( 38 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 214 )
  W( 74 ) = W( 74 ) + a*JVS( 215 )
  W( 77 ) = W( 77 ) + a*JVS( 216 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  a = -W( 63 ) / JVS(          533  )
  W( 63 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 534 )
  W( 73 ) = W( 73 ) + a*JVS( 535 )
  W( 74 ) = W( 74 ) + a*JVS( 536 )
  W( 75 ) = W( 75 ) + a*JVS( 537 )
  W( 76 ) = W( 76 ) + a*JVS( 538 )
  W( 77 ) = W( 77 ) + a*JVS( 539 )
  W( 79 ) = W( 79 ) + a*JVS( 540 )
  a = -W( 65 ) / JVS(          573  )
  W( 65 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 574 )
  W( 73 ) = W( 73 ) + a*JVS( 575 )
  W( 74 ) = W( 74 ) + a*JVS( 576 )
  W( 75 ) = W( 75 ) + a*JVS( 577 )
  W( 76 ) = W( 76 ) + a*JVS( 578 )
  W( 77 ) = W( 77 ) + a*JVS( 579 )
  W( 79 ) = W( 79 ) + a*JVS( 580 )
  JVS( 581) = W( 30 )
  JVS( 582) = W( 36 )
  JVS( 583) = W( 38 )
  JVS( 584) = W( 41 )
  JVS( 585) = W( 63 )
  JVS( 586) = W( 65 )
  JVS( 587) = W( 66 )
  JVS( 588) = W( 67 )
  JVS( 589) = W( 72 )
  JVS( 590) = W( 73 )
  JVS( 591) = W( 74 )
  JVS( 592) = W( 75 )
  JVS( 593) = W( 76 )
  JVS( 594) = W( 77 )
  JVS( 595) = W( 79 )
  IF ( ABS(  JVS( 603 )) < TINY(a) ) THEN
         IER = 67                                      
         RETURN
  END IF
   W( 30 ) = JVS( 596 )
   W( 36 ) = JVS( 597 )
   W( 38 ) = JVS( 598 )
   W( 41 ) = JVS( 599 )
   W( 63 ) = JVS( 600 )
   W( 65 ) = JVS( 601 )
   W( 66 ) = JVS( 602 )
   W( 67 ) = JVS( 603 )
   W( 72 ) = JVS( 604 )
   W( 73 ) = JVS( 605 )
   W( 74 ) = JVS( 606 )
   W( 75 ) = JVS( 607 )
   W( 76 ) = JVS( 608 )
   W( 77 ) = JVS( 609 )
   W( 79 ) = JVS( 610 )
  a = -W( 30 ) / JVS(          146  )
  W( 30 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 147 )
  W( 74 ) = W( 74 ) + a*JVS( 148 )
  W( 77 ) = W( 77 ) + a*JVS( 149 )
  a = -W( 36 ) / JVS(          185  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 186 )
  W( 74 ) = W( 74 ) + a*JVS( 187 )
  W( 77 ) = W( 77 ) + a*JVS( 188 )
  a = -W( 38 ) / JVS(          213  )
  W( 38 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 214 )
  W( 74 ) = W( 74 ) + a*JVS( 215 )
  W( 77 ) = W( 77 ) + a*JVS( 216 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  a = -W( 63 ) / JVS(          533  )
  W( 63 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 534 )
  W( 73 ) = W( 73 ) + a*JVS( 535 )
  W( 74 ) = W( 74 ) + a*JVS( 536 )
  W( 75 ) = W( 75 ) + a*JVS( 537 )
  W( 76 ) = W( 76 ) + a*JVS( 538 )
  W( 77 ) = W( 77 ) + a*JVS( 539 )
  W( 79 ) = W( 79 ) + a*JVS( 540 )
  a = -W( 65 ) / JVS(          573  )
  W( 65 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 574 )
  W( 73 ) = W( 73 ) + a*JVS( 575 )
  W( 74 ) = W( 74 ) + a*JVS( 576 )
  W( 75 ) = W( 75 ) + a*JVS( 577 )
  W( 76 ) = W( 76 ) + a*JVS( 578 )
  W( 77 ) = W( 77 ) + a*JVS( 579 )
  W( 79 ) = W( 79 ) + a*JVS( 580 )
  a = -W( 66 ) / JVS(          587  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 588 )
  W( 72 ) = W( 72 ) + a*JVS( 589 )
  W( 73 ) = W( 73 ) + a*JVS( 590 )
  W( 74 ) = W( 74 ) + a*JVS( 591 )
  W( 75 ) = W( 75 ) + a*JVS( 592 )
  W( 76 ) = W( 76 ) + a*JVS( 593 )
  W( 77 ) = W( 77 ) + a*JVS( 594 )
  W( 79 ) = W( 79 ) + a*JVS( 595 )
  JVS( 596) = W( 30 )
  JVS( 597) = W( 36 )
  JVS( 598) = W( 38 )
  JVS( 599) = W( 41 )
  JVS( 600) = W( 63 )
  JVS( 601) = W( 65 )
  JVS( 602) = W( 66 )
  JVS( 603) = W( 67 )
  JVS( 604) = W( 72 )
  JVS( 605) = W( 73 )
  JVS( 606) = W( 74 )
  JVS( 607) = W( 75 )
  JVS( 608) = W( 76 )
  JVS( 609) = W( 77 )
  JVS( 610) = W( 79 )
  IF ( ABS(  JVS( 624 )) < TINY(a) ) THEN
         IER = 68                                      
         RETURN
  END IF
   W( 17 ) = JVS( 611 )
   W( 36 ) = JVS( 612 )
   W( 38 ) = JVS( 613 )
   W( 49 ) = JVS( 614 )
   W( 53 ) = JVS( 615 )
   W( 55 ) = JVS( 616 )
   W( 59 ) = JVS( 617 )
   W( 60 ) = JVS( 618 )
   W( 61 ) = JVS( 619 )
   W( 63 ) = JVS( 620 )
   W( 65 ) = JVS( 621 )
   W( 66 ) = JVS( 622 )
   W( 67 ) = JVS( 623 )
   W( 68 ) = JVS( 624 )
   W( 69 ) = JVS( 625 )
   W( 70 ) = JVS( 626 )
   W( 71 ) = JVS( 627 )
   W( 72 ) = JVS( 628 )
   W( 73 ) = JVS( 629 )
   W( 74 ) = JVS( 630 )
   W( 75 ) = JVS( 631 )
   W( 76 ) = JVS( 632 )
   W( 77 ) = JVS( 633 )
   W( 78 ) = JVS( 634 )
   W( 79 ) = JVS( 635 )
   W( 80 ) = JVS( 636 )
  a = -W( 17 ) / JVS(           73  )
  W( 17 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 74 )
  W( 65 ) = W( 65 ) + a*JVS( 75 )
  W( 73 ) = W( 73 ) + a*JVS( 76 )
  W( 77 ) = W( 77 ) + a*JVS( 77 )
  a = -W( 36 ) / JVS(          185  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 186 )
  W( 74 ) = W( 74 ) + a*JVS( 187 )
  W( 77 ) = W( 77 ) + a*JVS( 188 )
  a = -W( 38 ) / JVS(          213  )
  W( 38 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 214 )
  W( 74 ) = W( 74 ) + a*JVS( 215 )
  W( 77 ) = W( 77 ) + a*JVS( 216 )
  a = -W( 49 ) / JVS(          309  )
  W( 49 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 310 )
  W( 55 ) = W( 55 ) + a*JVS( 311 )
  W( 59 ) = W( 59 ) + a*JVS( 312 )
  W( 60 ) = W( 60 ) + a*JVS( 313 )
  W( 61 ) = W( 61 ) + a*JVS( 314 )
  W( 63 ) = W( 63 ) + a*JVS( 315 )
  W( 65 ) = W( 65 ) + a*JVS( 316 )
  W( 66 ) = W( 66 ) + a*JVS( 317 )
  W( 67 ) = W( 67 ) + a*JVS( 318 )
  W( 69 ) = W( 69 ) + a*JVS( 319 )
  W( 70 ) = W( 70 ) + a*JVS( 320 )
  W( 71 ) = W( 71 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  W( 73 ) = W( 73 ) + a*JVS( 323 )
  W( 74 ) = W( 74 ) + a*JVS( 324 )
  W( 75 ) = W( 75 ) + a*JVS( 325 )
  W( 76 ) = W( 76 ) + a*JVS( 326 )
  W( 77 ) = W( 77 ) + a*JVS( 327 )
  W( 79 ) = W( 79 ) + a*JVS( 328 )
  a = -W( 53 ) / JVS(          376  )
  W( 53 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 377 )
  W( 73 ) = W( 73 ) + a*JVS( 378 )
  W( 74 ) = W( 74 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 79 ) = W( 79 ) + a*JVS( 382 )
  a = -W( 55 ) / JVS(          393  )
  W( 55 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  W( 73 ) = W( 73 ) + a*JVS( 395 )
  W( 74 ) = W( 74 ) + a*JVS( 396 )
  W( 75 ) = W( 75 ) + a*JVS( 397 )
  W( 76 ) = W( 76 ) + a*JVS( 398 )
  W( 79 ) = W( 79 ) + a*JVS( 399 )
  a = -W( 59 ) / JVS(          466  )
  W( 59 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 467 )
  W( 73 ) = W( 73 ) + a*JVS( 468 )
  W( 74 ) = W( 74 ) + a*JVS( 469 )
  W( 75 ) = W( 75 ) + a*JVS( 470 )
  W( 76 ) = W( 76 ) + a*JVS( 471 )
  W( 77 ) = W( 77 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  a = -W( 60 ) / JVS(          483  )
  W( 60 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 484 )
  W( 67 ) = W( 67 ) + a*JVS( 485 )
  W( 71 ) = W( 71 ) + a*JVS( 486 )
  W( 72 ) = W( 72 ) + a*JVS( 487 )
  W( 73 ) = W( 73 ) + a*JVS( 488 )
  W( 74 ) = W( 74 ) + a*JVS( 489 )
  W( 75 ) = W( 75 ) + a*JVS( 490 )
  W( 76 ) = W( 76 ) + a*JVS( 491 )
  W( 77 ) = W( 77 ) + a*JVS( 492 )
  W( 78 ) = W( 78 ) + a*JVS( 493 )
  W( 79 ) = W( 79 ) + a*JVS( 494 )
  W( 80 ) = W( 80 ) + a*JVS( 495 )
  a = -W( 61 ) / JVS(          498  )
  W( 61 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 499 )
  W( 73 ) = W( 73 ) + a*JVS( 500 )
  W( 74 ) = W( 74 ) + a*JVS( 501 )
  W( 75 ) = W( 75 ) + a*JVS( 502 )
  W( 76 ) = W( 76 ) + a*JVS( 503 )
  W( 77 ) = W( 77 ) + a*JVS( 504 )
  W( 79 ) = W( 79 ) + a*JVS( 505 )
  W( 80 ) = W( 80 ) + a*JVS( 506 )
  a = -W( 63 ) / JVS(          533  )
  W( 63 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 534 )
  W( 73 ) = W( 73 ) + a*JVS( 535 )
  W( 74 ) = W( 74 ) + a*JVS( 536 )
  W( 75 ) = W( 75 ) + a*JVS( 537 )
  W( 76 ) = W( 76 ) + a*JVS( 538 )
  W( 77 ) = W( 77 ) + a*JVS( 539 )
  W( 79 ) = W( 79 ) + a*JVS( 540 )
  a = -W( 65 ) / JVS(          573  )
  W( 65 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 574 )
  W( 73 ) = W( 73 ) + a*JVS( 575 )
  W( 74 ) = W( 74 ) + a*JVS( 576 )
  W( 75 ) = W( 75 ) + a*JVS( 577 )
  W( 76 ) = W( 76 ) + a*JVS( 578 )
  W( 77 ) = W( 77 ) + a*JVS( 579 )
  W( 79 ) = W( 79 ) + a*JVS( 580 )
  a = -W( 66 ) / JVS(          587  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 588 )
  W( 72 ) = W( 72 ) + a*JVS( 589 )
  W( 73 ) = W( 73 ) + a*JVS( 590 )
  W( 74 ) = W( 74 ) + a*JVS( 591 )
  W( 75 ) = W( 75 ) + a*JVS( 592 )
  W( 76 ) = W( 76 ) + a*JVS( 593 )
  W( 77 ) = W( 77 ) + a*JVS( 594 )
  W( 79 ) = W( 79 ) + a*JVS( 595 )
  a = -W( 67 ) / JVS(          603  )
  W( 67 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 604 )
  W( 73 ) = W( 73 ) + a*JVS( 605 )
  W( 74 ) = W( 74 ) + a*JVS( 606 )
  W( 75 ) = W( 75 ) + a*JVS( 607 )
  W( 76 ) = W( 76 ) + a*JVS( 608 )
  W( 77 ) = W( 77 ) + a*JVS( 609 )
  W( 79 ) = W( 79 ) + a*JVS( 610 )
  JVS( 611) = W( 17 )
  JVS( 612) = W( 36 )
  JVS( 613) = W( 38 )
  JVS( 614) = W( 49 )
  JVS( 615) = W( 53 )
  JVS( 616) = W( 55 )
  JVS( 617) = W( 59 )
  JVS( 618) = W( 60 )
  JVS( 619) = W( 61 )
  JVS( 620) = W( 63 )
  JVS( 621) = W( 65 )
  JVS( 622) = W( 66 )
  JVS( 623) = W( 67 )
  JVS( 624) = W( 68 )
  JVS( 625) = W( 69 )
  JVS( 626) = W( 70 )
  JVS( 627) = W( 71 )
  JVS( 628) = W( 72 )
  JVS( 629) = W( 73 )
  JVS( 630) = W( 74 )
  JVS( 631) = W( 75 )
  JVS( 632) = W( 76 )
  JVS( 633) = W( 77 )
  JVS( 634) = W( 78 )
  JVS( 635) = W( 79 )
  JVS( 636) = W( 80 )
  IF ( ABS(  JVS( 650 )) < TINY(a) ) THEN
         IER = 69                                      
         RETURN
  END IF
   W( 36 ) = JVS( 637 )
   W( 38 ) = JVS( 638 )
   W( 41 ) = JVS( 639 )
   W( 49 ) = JVS( 640 )
   W( 53 ) = JVS( 641 )
   W( 55 ) = JVS( 642 )
   W( 59 ) = JVS( 643 )
   W( 60 ) = JVS( 644 )
   W( 61 ) = JVS( 645 )
   W( 63 ) = JVS( 646 )
   W( 65 ) = JVS( 647 )
   W( 66 ) = JVS( 648 )
   W( 67 ) = JVS( 649 )
   W( 69 ) = JVS( 650 )
   W( 70 ) = JVS( 651 )
   W( 71 ) = JVS( 652 )
   W( 72 ) = JVS( 653 )
   W( 73 ) = JVS( 654 )
   W( 74 ) = JVS( 655 )
   W( 75 ) = JVS( 656 )
   W( 76 ) = JVS( 657 )
   W( 77 ) = JVS( 658 )
   W( 78 ) = JVS( 659 )
   W( 79 ) = JVS( 660 )
   W( 80 ) = JVS( 661 )
  a = -W( 36 ) / JVS(          185  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 186 )
  W( 74 ) = W( 74 ) + a*JVS( 187 )
  W( 77 ) = W( 77 ) + a*JVS( 188 )
  a = -W( 38 ) / JVS(          213  )
  W( 38 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 214 )
  W( 74 ) = W( 74 ) + a*JVS( 215 )
  W( 77 ) = W( 77 ) + a*JVS( 216 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  a = -W( 49 ) / JVS(          309  )
  W( 49 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 310 )
  W( 55 ) = W( 55 ) + a*JVS( 311 )
  W( 59 ) = W( 59 ) + a*JVS( 312 )
  W( 60 ) = W( 60 ) + a*JVS( 313 )
  W( 61 ) = W( 61 ) + a*JVS( 314 )
  W( 63 ) = W( 63 ) + a*JVS( 315 )
  W( 65 ) = W( 65 ) + a*JVS( 316 )
  W( 66 ) = W( 66 ) + a*JVS( 317 )
  W( 67 ) = W( 67 ) + a*JVS( 318 )
  W( 69 ) = W( 69 ) + a*JVS( 319 )
  W( 70 ) = W( 70 ) + a*JVS( 320 )
  W( 71 ) = W( 71 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  W( 73 ) = W( 73 ) + a*JVS( 323 )
  W( 74 ) = W( 74 ) + a*JVS( 324 )
  W( 75 ) = W( 75 ) + a*JVS( 325 )
  W( 76 ) = W( 76 ) + a*JVS( 326 )
  W( 77 ) = W( 77 ) + a*JVS( 327 )
  W( 79 ) = W( 79 ) + a*JVS( 328 )
  a = -W( 53 ) / JVS(          376  )
  W( 53 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 377 )
  W( 73 ) = W( 73 ) + a*JVS( 378 )
  W( 74 ) = W( 74 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 79 ) = W( 79 ) + a*JVS( 382 )
  a = -W( 55 ) / JVS(          393  )
  W( 55 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  W( 73 ) = W( 73 ) + a*JVS( 395 )
  W( 74 ) = W( 74 ) + a*JVS( 396 )
  W( 75 ) = W( 75 ) + a*JVS( 397 )
  W( 76 ) = W( 76 ) + a*JVS( 398 )
  W( 79 ) = W( 79 ) + a*JVS( 399 )
  a = -W( 59 ) / JVS(          466  )
  W( 59 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 467 )
  W( 73 ) = W( 73 ) + a*JVS( 468 )
  W( 74 ) = W( 74 ) + a*JVS( 469 )
  W( 75 ) = W( 75 ) + a*JVS( 470 )
  W( 76 ) = W( 76 ) + a*JVS( 471 )
  W( 77 ) = W( 77 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  a = -W( 60 ) / JVS(          483  )
  W( 60 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 484 )
  W( 67 ) = W( 67 ) + a*JVS( 485 )
  W( 71 ) = W( 71 ) + a*JVS( 486 )
  W( 72 ) = W( 72 ) + a*JVS( 487 )
  W( 73 ) = W( 73 ) + a*JVS( 488 )
  W( 74 ) = W( 74 ) + a*JVS( 489 )
  W( 75 ) = W( 75 ) + a*JVS( 490 )
  W( 76 ) = W( 76 ) + a*JVS( 491 )
  W( 77 ) = W( 77 ) + a*JVS( 492 )
  W( 78 ) = W( 78 ) + a*JVS( 493 )
  W( 79 ) = W( 79 ) + a*JVS( 494 )
  W( 80 ) = W( 80 ) + a*JVS( 495 )
  a = -W( 61 ) / JVS(          498  )
  W( 61 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 499 )
  W( 73 ) = W( 73 ) + a*JVS( 500 )
  W( 74 ) = W( 74 ) + a*JVS( 501 )
  W( 75 ) = W( 75 ) + a*JVS( 502 )
  W( 76 ) = W( 76 ) + a*JVS( 503 )
  W( 77 ) = W( 77 ) + a*JVS( 504 )
  W( 79 ) = W( 79 ) + a*JVS( 505 )
  W( 80 ) = W( 80 ) + a*JVS( 506 )
  a = -W( 63 ) / JVS(          533  )
  W( 63 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 534 )
  W( 73 ) = W( 73 ) + a*JVS( 535 )
  W( 74 ) = W( 74 ) + a*JVS( 536 )
  W( 75 ) = W( 75 ) + a*JVS( 537 )
  W( 76 ) = W( 76 ) + a*JVS( 538 )
  W( 77 ) = W( 77 ) + a*JVS( 539 )
  W( 79 ) = W( 79 ) + a*JVS( 540 )
  a = -W( 65 ) / JVS(          573  )
  W( 65 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 574 )
  W( 73 ) = W( 73 ) + a*JVS( 575 )
  W( 74 ) = W( 74 ) + a*JVS( 576 )
  W( 75 ) = W( 75 ) + a*JVS( 577 )
  W( 76 ) = W( 76 ) + a*JVS( 578 )
  W( 77 ) = W( 77 ) + a*JVS( 579 )
  W( 79 ) = W( 79 ) + a*JVS( 580 )
  a = -W( 66 ) / JVS(          587  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 588 )
  W( 72 ) = W( 72 ) + a*JVS( 589 )
  W( 73 ) = W( 73 ) + a*JVS( 590 )
  W( 74 ) = W( 74 ) + a*JVS( 591 )
  W( 75 ) = W( 75 ) + a*JVS( 592 )
  W( 76 ) = W( 76 ) + a*JVS( 593 )
  W( 77 ) = W( 77 ) + a*JVS( 594 )
  W( 79 ) = W( 79 ) + a*JVS( 595 )
  a = -W( 67 ) / JVS(          603  )
  W( 67 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 604 )
  W( 73 ) = W( 73 ) + a*JVS( 605 )
  W( 74 ) = W( 74 ) + a*JVS( 606 )
  W( 75 ) = W( 75 ) + a*JVS( 607 )
  W( 76 ) = W( 76 ) + a*JVS( 608 )
  W( 77 ) = W( 77 ) + a*JVS( 609 )
  W( 79 ) = W( 79 ) + a*JVS( 610 )
  JVS( 637) = W( 36 )
  JVS( 638) = W( 38 )
  JVS( 639) = W( 41 )
  JVS( 640) = W( 49 )
  JVS( 641) = W( 53 )
  JVS( 642) = W( 55 )
  JVS( 643) = W( 59 )
  JVS( 644) = W( 60 )
  JVS( 645) = W( 61 )
  JVS( 646) = W( 63 )
  JVS( 647) = W( 65 )
  JVS( 648) = W( 66 )
  JVS( 649) = W( 67 )
  JVS( 650) = W( 69 )
  JVS( 651) = W( 70 )
  JVS( 652) = W( 71 )
  JVS( 653) = W( 72 )
  JVS( 654) = W( 73 )
  JVS( 655) = W( 74 )
  JVS( 656) = W( 75 )
  JVS( 657) = W( 76 )
  JVS( 658) = W( 77 )
  JVS( 659) = W( 78 )
  JVS( 660) = W( 79 )
  JVS( 661) = W( 80 )
  IF ( ABS(  JVS( 677 )) < TINY(a) ) THEN
         IER = 70                                      
         RETURN
  END IF
   W( 45 ) = JVS( 662 )
   W( 46 ) = JVS( 663 )
   W( 47 ) = JVS( 664 )
   W( 52 ) = JVS( 665 )
   W( 53 ) = JVS( 666 )
   W( 54 ) = JVS( 667 )
   W( 55 ) = JVS( 668 )
   W( 57 ) = JVS( 669 )
   W( 58 ) = JVS( 670 )
   W( 59 ) = JVS( 671 )
   W( 62 ) = JVS( 672 )
   W( 63 ) = JVS( 673 )
   W( 65 ) = JVS( 674 )
   W( 68 ) = JVS( 675 )
   W( 69 ) = JVS( 676 )
   W( 70 ) = JVS( 677 )
   W( 71 ) = JVS( 678 )
   W( 72 ) = JVS( 679 )
   W( 73 ) = JVS( 680 )
   W( 74 ) = JVS( 681 )
   W( 75 ) = JVS( 682 )
   W( 76 ) = JVS( 683 )
   W( 77 ) = JVS( 684 )
   W( 78 ) = JVS( 685 )
   W( 79 ) = JVS( 686 )
   W( 80 ) = JVS( 687 )
  a = -W( 45 ) / JVS(          266  )
  W( 45 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 267 )
  W( 73 ) = W( 73 ) + a*JVS( 268 )
  W( 74 ) = W( 74 ) + a*JVS( 269 )
  W( 75 ) = W( 75 ) + a*JVS( 270 )
  W( 76 ) = W( 76 ) + a*JVS( 271 )
  W( 77 ) = W( 77 ) + a*JVS( 272 )
  W( 79 ) = W( 79 ) + a*JVS( 273 )
  a = -W( 46 ) / JVS(          274  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 275 )
  W( 72 ) = W( 72 ) + a*JVS( 276 )
  W( 73 ) = W( 73 ) + a*JVS( 277 )
  W( 74 ) = W( 74 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 76 ) = W( 76 ) + a*JVS( 280 )
  W( 79 ) = W( 79 ) + a*JVS( 281 )
  a = -W( 47 ) / JVS(          282  )
  W( 47 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 283 )
  W( 72 ) = W( 72 ) + a*JVS( 284 )
  W( 73 ) = W( 73 ) + a*JVS( 285 )
  W( 74 ) = W( 74 ) + a*JVS( 286 )
  W( 75 ) = W( 75 ) + a*JVS( 287 )
  W( 76 ) = W( 76 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  a = -W( 52 ) / JVS(          366  )
  W( 52 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 367 )
  W( 73 ) = W( 73 ) + a*JVS( 368 )
  W( 74 ) = W( 74 ) + a*JVS( 369 )
  W( 75 ) = W( 75 ) + a*JVS( 370 )
  W( 76 ) = W( 76 ) + a*JVS( 371 )
  W( 77 ) = W( 77 ) + a*JVS( 372 )
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 80 ) = W( 80 ) + a*JVS( 374 )
  a = -W( 53 ) / JVS(          376  )
  W( 53 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 377 )
  W( 73 ) = W( 73 ) + a*JVS( 378 )
  W( 74 ) = W( 74 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 79 ) = W( 79 ) + a*JVS( 382 )
  a = -W( 54 ) / JVS(          384  )
  W( 54 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 385 )
  W( 73 ) = W( 73 ) + a*JVS( 386 )
  W( 74 ) = W( 74 ) + a*JVS( 387 )
  W( 75 ) = W( 75 ) + a*JVS( 388 )
  W( 76 ) = W( 76 ) + a*JVS( 389 )
  W( 77 ) = W( 77 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 55 ) / JVS(          393  )
  W( 55 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  W( 73 ) = W( 73 ) + a*JVS( 395 )
  W( 74 ) = W( 74 ) + a*JVS( 396 )
  W( 75 ) = W( 75 ) + a*JVS( 397 )
  W( 76 ) = W( 76 ) + a*JVS( 398 )
  W( 79 ) = W( 79 ) + a*JVS( 399 )
  a = -W( 57 ) / JVS(          446  )
  W( 57 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 447 )
  W( 73 ) = W( 73 ) + a*JVS( 448 )
  W( 74 ) = W( 74 ) + a*JVS( 449 )
  W( 75 ) = W( 75 ) + a*JVS( 450 )
  W( 76 ) = W( 76 ) + a*JVS( 451 )
  W( 77 ) = W( 77 ) + a*JVS( 452 )
  W( 79 ) = W( 79 ) + a*JVS( 453 )
  W( 80 ) = W( 80 ) + a*JVS( 454 )
  a = -W( 58 ) / JVS(          456  )
  W( 58 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 457 )
  W( 73 ) = W( 73 ) + a*JVS( 458 )
  W( 74 ) = W( 74 ) + a*JVS( 459 )
  W( 75 ) = W( 75 ) + a*JVS( 460 )
  W( 76 ) = W( 76 ) + a*JVS( 461 )
  W( 77 ) = W( 77 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  a = -W( 59 ) / JVS(          466  )
  W( 59 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 467 )
  W( 73 ) = W( 73 ) + a*JVS( 468 )
  W( 74 ) = W( 74 ) + a*JVS( 469 )
  W( 75 ) = W( 75 ) + a*JVS( 470 )
  W( 76 ) = W( 76 ) + a*JVS( 471 )
  W( 77 ) = W( 77 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  a = -W( 62 ) / JVS(          517  )
  W( 62 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 518 )
  W( 70 ) = W( 70 ) + a*JVS( 519 )
  W( 71 ) = W( 71 ) + a*JVS( 520 )
  W( 72 ) = W( 72 ) + a*JVS( 521 )
  W( 73 ) = W( 73 ) + a*JVS( 522 )
  W( 74 ) = W( 74 ) + a*JVS( 523 )
  W( 75 ) = W( 75 ) + a*JVS( 524 )
  W( 76 ) = W( 76 ) + a*JVS( 525 )
  W( 77 ) = W( 77 ) + a*JVS( 526 )
  W( 78 ) = W( 78 ) + a*JVS( 527 )
  W( 79 ) = W( 79 ) + a*JVS( 528 )
  W( 80 ) = W( 80 ) + a*JVS( 529 )
  a = -W( 63 ) / JVS(          533  )
  W( 63 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 534 )
  W( 73 ) = W( 73 ) + a*JVS( 535 )
  W( 74 ) = W( 74 ) + a*JVS( 536 )
  W( 75 ) = W( 75 ) + a*JVS( 537 )
  W( 76 ) = W( 76 ) + a*JVS( 538 )
  W( 77 ) = W( 77 ) + a*JVS( 539 )
  W( 79 ) = W( 79 ) + a*JVS( 540 )
  a = -W( 65 ) / JVS(          573  )
  W( 65 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 574 )
  W( 73 ) = W( 73 ) + a*JVS( 575 )
  W( 74 ) = W( 74 ) + a*JVS( 576 )
  W( 75 ) = W( 75 ) + a*JVS( 577 )
  W( 76 ) = W( 76 ) + a*JVS( 578 )
  W( 77 ) = W( 77 ) + a*JVS( 579 )
  W( 79 ) = W( 79 ) + a*JVS( 580 )
  a = -W( 68 ) / JVS(          624  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 625 )
  W( 70 ) = W( 70 ) + a*JVS( 626 )
  W( 71 ) = W( 71 ) + a*JVS( 627 )
  W( 72 ) = W( 72 ) + a*JVS( 628 )
  W( 73 ) = W( 73 ) + a*JVS( 629 )
  W( 74 ) = W( 74 ) + a*JVS( 630 )
  W( 75 ) = W( 75 ) + a*JVS( 631 )
  W( 76 ) = W( 76 ) + a*JVS( 632 )
  W( 77 ) = W( 77 ) + a*JVS( 633 )
  W( 78 ) = W( 78 ) + a*JVS( 634 )
  W( 79 ) = W( 79 ) + a*JVS( 635 )
  W( 80 ) = W( 80 ) + a*JVS( 636 )
  a = -W( 69 ) / JVS(          650  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 651 )
  W( 71 ) = W( 71 ) + a*JVS( 652 )
  W( 72 ) = W( 72 ) + a*JVS( 653 )
  W( 73 ) = W( 73 ) + a*JVS( 654 )
  W( 74 ) = W( 74 ) + a*JVS( 655 )
  W( 75 ) = W( 75 ) + a*JVS( 656 )
  W( 76 ) = W( 76 ) + a*JVS( 657 )
  W( 77 ) = W( 77 ) + a*JVS( 658 )
  W( 78 ) = W( 78 ) + a*JVS( 659 )
  W( 79 ) = W( 79 ) + a*JVS( 660 )
  W( 80 ) = W( 80 ) + a*JVS( 661 )
  JVS( 662) = W( 45 )
  JVS( 663) = W( 46 )
  JVS( 664) = W( 47 )
  JVS( 665) = W( 52 )
  JVS( 666) = W( 53 )
  JVS( 667) = W( 54 )
  JVS( 668) = W( 55 )
  JVS( 669) = W( 57 )
  JVS( 670) = W( 58 )
  JVS( 671) = W( 59 )
  JVS( 672) = W( 62 )
  JVS( 673) = W( 63 )
  JVS( 674) = W( 65 )
  JVS( 675) = W( 68 )
  JVS( 676) = W( 69 )
  JVS( 677) = W( 70 )
  JVS( 678) = W( 71 )
  JVS( 679) = W( 72 )
  JVS( 680) = W( 73 )
  JVS( 681) = W( 74 )
  JVS( 682) = W( 75 )
  JVS( 683) = W( 76 )
  JVS( 684) = W( 77 )
  JVS( 685) = W( 78 )
  JVS( 686) = W( 79 )
  JVS( 687) = W( 80 )
  IF ( ABS(  JVS( 693 )) < TINY(a) ) THEN
         IER = 71                                      
         RETURN
  END IF
   W( 14 ) = JVS( 688 )
   W( 60 ) = JVS( 689 )
   W( 66 ) = JVS( 690 )
   W( 67 ) = JVS( 691 )
   W( 70 ) = JVS( 692 )
   W( 71 ) = JVS( 693 )
   W( 72 ) = JVS( 694 )
   W( 73 ) = JVS( 695 )
   W( 74 ) = JVS( 696 )
   W( 75 ) = JVS( 697 )
   W( 76 ) = JVS( 698 )
   W( 77 ) = JVS( 699 )
   W( 78 ) = JVS( 700 )
   W( 79 ) = JVS( 701 )
   W( 80 ) = JVS( 702 )
  a = -W( 14 ) / JVS(           63  )
  W( 14 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 64 )
  a = -W( 60 ) / JVS(          483  )
  W( 60 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 484 )
  W( 67 ) = W( 67 ) + a*JVS( 485 )
  W( 71 ) = W( 71 ) + a*JVS( 486 )
  W( 72 ) = W( 72 ) + a*JVS( 487 )
  W( 73 ) = W( 73 ) + a*JVS( 488 )
  W( 74 ) = W( 74 ) + a*JVS( 489 )
  W( 75 ) = W( 75 ) + a*JVS( 490 )
  W( 76 ) = W( 76 ) + a*JVS( 491 )
  W( 77 ) = W( 77 ) + a*JVS( 492 )
  W( 78 ) = W( 78 ) + a*JVS( 493 )
  W( 79 ) = W( 79 ) + a*JVS( 494 )
  W( 80 ) = W( 80 ) + a*JVS( 495 )
  a = -W( 66 ) / JVS(          587  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 588 )
  W( 72 ) = W( 72 ) + a*JVS( 589 )
  W( 73 ) = W( 73 ) + a*JVS( 590 )
  W( 74 ) = W( 74 ) + a*JVS( 591 )
  W( 75 ) = W( 75 ) + a*JVS( 592 )
  W( 76 ) = W( 76 ) + a*JVS( 593 )
  W( 77 ) = W( 77 ) + a*JVS( 594 )
  W( 79 ) = W( 79 ) + a*JVS( 595 )
  a = -W( 67 ) / JVS(          603  )
  W( 67 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 604 )
  W( 73 ) = W( 73 ) + a*JVS( 605 )
  W( 74 ) = W( 74 ) + a*JVS( 606 )
  W( 75 ) = W( 75 ) + a*JVS( 607 )
  W( 76 ) = W( 76 ) + a*JVS( 608 )
  W( 77 ) = W( 77 ) + a*JVS( 609 )
  W( 79 ) = W( 79 ) + a*JVS( 610 )
  a = -W( 70 ) / JVS(          677  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 678 )
  W( 72 ) = W( 72 ) + a*JVS( 679 )
  W( 73 ) = W( 73 ) + a*JVS( 680 )
  W( 74 ) = W( 74 ) + a*JVS( 681 )
  W( 75 ) = W( 75 ) + a*JVS( 682 )
  W( 76 ) = W( 76 ) + a*JVS( 683 )
  W( 77 ) = W( 77 ) + a*JVS( 684 )
  W( 78 ) = W( 78 ) + a*JVS( 685 )
  W( 79 ) = W( 79 ) + a*JVS( 686 )
  W( 80 ) = W( 80 ) + a*JVS( 687 )
  JVS( 688) = W( 14 )
  JVS( 689) = W( 60 )
  JVS( 690) = W( 66 )
  JVS( 691) = W( 67 )
  JVS( 692) = W( 70 )
  JVS( 693) = W( 71 )
  JVS( 694) = W( 72 )
  JVS( 695) = W( 73 )
  JVS( 696) = W( 74 )
  JVS( 697) = W( 75 )
  JVS( 698) = W( 76 )
  JVS( 699) = W( 77 )
  JVS( 700) = W( 78 )
  JVS( 701) = W( 79 )
  JVS( 702) = W( 80 )
  IF ( ABS(  JVS( 736 )) < TINY(a) ) THEN
         IER = 72                                      
         RETURN
  END IF
   W( 25 ) = JVS( 703 )
   W( 33 ) = JVS( 704 )
   W( 35 ) = JVS( 705 )
   W( 39 ) = JVS( 706 )
   W( 41 ) = JVS( 707 )
   W( 42 ) = JVS( 708 )
   W( 43 ) = JVS( 709 )
   W( 45 ) = JVS( 710 )
   W( 46 ) = JVS( 711 )
   W( 47 ) = JVS( 712 )
   W( 48 ) = JVS( 713 )
   W( 49 ) = JVS( 714 )
   W( 50 ) = JVS( 715 )
   W( 51 ) = JVS( 716 )
   W( 52 ) = JVS( 717 )
   W( 53 ) = JVS( 718 )
   W( 54 ) = JVS( 719 )
   W( 55 ) = JVS( 720 )
   W( 57 ) = JVS( 721 )
   W( 58 ) = JVS( 722 )
   W( 59 ) = JVS( 723 )
   W( 60 ) = JVS( 724 )
   W( 61 ) = JVS( 725 )
   W( 62 ) = JVS( 726 )
   W( 63 ) = JVS( 727 )
   W( 64 ) = JVS( 728 )
   W( 65 ) = JVS( 729 )
   W( 66 ) = JVS( 730 )
   W( 67 ) = JVS( 731 )
   W( 68 ) = JVS( 732 )
   W( 69 ) = JVS( 733 )
   W( 70 ) = JVS( 734 )
   W( 71 ) = JVS( 735 )
   W( 72 ) = JVS( 736 )
   W( 73 ) = JVS( 737 )
   W( 74 ) = JVS( 738 )
   W( 75 ) = JVS( 739 )
   W( 76 ) = JVS( 740 )
   W( 77 ) = JVS( 741 )
   W( 78 ) = JVS( 742 )
   W( 79 ) = JVS( 743 )
   W( 80 ) = JVS( 744 )
  a = -W( 25 ) / JVS(          118  )
  W( 25 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 119 )
  W( 69 ) = W( 69 ) + a*JVS( 120 )
  W( 73 ) = W( 73 ) + a*JVS( 121 )
  W( 76 ) = W( 76 ) + a*JVS( 122 )
  a = -W( 33 ) / JVS(          166  )
  W( 33 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 167 )
  W( 72 ) = W( 72 ) + a*JVS( 168 )
  W( 73 ) = W( 73 ) + a*JVS( 169 )
  W( 75 ) = W( 75 ) + a*JVS( 170 )
  W( 77 ) = W( 77 ) + a*JVS( 171 )
  a = -W( 35 ) / JVS(          178  )
  W( 35 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 179 )
  W( 72 ) = W( 72 ) + a*JVS( 180 )
  W( 73 ) = W( 73 ) + a*JVS( 181 )
  W( 74 ) = W( 74 ) + a*JVS( 182 )
  W( 77 ) = W( 77 ) + a*JVS( 183 )
  W( 80 ) = W( 80 ) + a*JVS( 184 )
  a = -W( 39 ) / JVS(          217  )
  W( 39 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 218 )
  W( 74 ) = W( 74 ) + a*JVS( 219 )
  W( 77 ) = W( 77 ) + a*JVS( 220 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  a = -W( 42 ) / JVS(          238  )
  W( 42 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 239 )
  W( 73 ) = W( 73 ) + a*JVS( 240 )
  W( 74 ) = W( 74 ) + a*JVS( 241 )
  W( 75 ) = W( 75 ) + a*JVS( 242 )
  W( 77 ) = W( 77 ) + a*JVS( 243 )
  W( 79 ) = W( 79 ) + a*JVS( 244 )
  W( 80 ) = W( 80 ) + a*JVS( 245 )
  a = -W( 43 ) / JVS(          246  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 247 )
  W( 74 ) = W( 74 ) + a*JVS( 248 )
  W( 77 ) = W( 77 ) + a*JVS( 249 )
  W( 78 ) = W( 78 ) + a*JVS( 250 )
  W( 80 ) = W( 80 ) + a*JVS( 251 )
  a = -W( 45 ) / JVS(          266  )
  W( 45 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 267 )
  W( 73 ) = W( 73 ) + a*JVS( 268 )
  W( 74 ) = W( 74 ) + a*JVS( 269 )
  W( 75 ) = W( 75 ) + a*JVS( 270 )
  W( 76 ) = W( 76 ) + a*JVS( 271 )
  W( 77 ) = W( 77 ) + a*JVS( 272 )
  W( 79 ) = W( 79 ) + a*JVS( 273 )
  a = -W( 46 ) / JVS(          274  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 275 )
  W( 72 ) = W( 72 ) + a*JVS( 276 )
  W( 73 ) = W( 73 ) + a*JVS( 277 )
  W( 74 ) = W( 74 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 76 ) = W( 76 ) + a*JVS( 280 )
  W( 79 ) = W( 79 ) + a*JVS( 281 )
  a = -W( 47 ) / JVS(          282  )
  W( 47 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 283 )
  W( 72 ) = W( 72 ) + a*JVS( 284 )
  W( 73 ) = W( 73 ) + a*JVS( 285 )
  W( 74 ) = W( 74 ) + a*JVS( 286 )
  W( 75 ) = W( 75 ) + a*JVS( 287 )
  W( 76 ) = W( 76 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  a = -W( 48 ) / JVS(          294  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 295 )
  W( 54 ) = W( 54 ) + a*JVS( 296 )
  W( 72 ) = W( 72 ) + a*JVS( 297 )
  W( 73 ) = W( 73 ) + a*JVS( 298 )
  W( 74 ) = W( 74 ) + a*JVS( 299 )
  W( 75 ) = W( 75 ) + a*JVS( 300 )
  W( 76 ) = W( 76 ) + a*JVS( 301 )
  W( 77 ) = W( 77 ) + a*JVS( 302 )
  W( 79 ) = W( 79 ) + a*JVS( 303 )
  a = -W( 49 ) / JVS(          309  )
  W( 49 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 310 )
  W( 55 ) = W( 55 ) + a*JVS( 311 )
  W( 59 ) = W( 59 ) + a*JVS( 312 )
  W( 60 ) = W( 60 ) + a*JVS( 313 )
  W( 61 ) = W( 61 ) + a*JVS( 314 )
  W( 63 ) = W( 63 ) + a*JVS( 315 )
  W( 65 ) = W( 65 ) + a*JVS( 316 )
  W( 66 ) = W( 66 ) + a*JVS( 317 )
  W( 67 ) = W( 67 ) + a*JVS( 318 )
  W( 69 ) = W( 69 ) + a*JVS( 319 )
  W( 70 ) = W( 70 ) + a*JVS( 320 )
  W( 71 ) = W( 71 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  W( 73 ) = W( 73 ) + a*JVS( 323 )
  W( 74 ) = W( 74 ) + a*JVS( 324 )
  W( 75 ) = W( 75 ) + a*JVS( 325 )
  W( 76 ) = W( 76 ) + a*JVS( 326 )
  W( 77 ) = W( 77 ) + a*JVS( 327 )
  W( 79 ) = W( 79 ) + a*JVS( 328 )
  a = -W( 50 ) / JVS(          331  )
  W( 50 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 332 )
  W( 73 ) = W( 73 ) + a*JVS( 333 )
  W( 74 ) = W( 74 ) + a*JVS( 334 )
  W( 75 ) = W( 75 ) + a*JVS( 335 )
  W( 76 ) = W( 76 ) + a*JVS( 336 )
  W( 77 ) = W( 77 ) + a*JVS( 337 )
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  a = -W( 51 ) / JVS(          347  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 348 )
  W( 54 ) = W( 54 ) + a*JVS( 349 )
  W( 57 ) = W( 57 ) + a*JVS( 350 )
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 61 ) = W( 61 ) + a*JVS( 352 )
  W( 69 ) = W( 69 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 72 ) = W( 72 ) + a*JVS( 355 )
  W( 73 ) = W( 73 ) + a*JVS( 356 )
  W( 74 ) = W( 74 ) + a*JVS( 357 )
  W( 75 ) = W( 75 ) + a*JVS( 358 )
  W( 76 ) = W( 76 ) + a*JVS( 359 )
  W( 77 ) = W( 77 ) + a*JVS( 360 )
  W( 78 ) = W( 78 ) + a*JVS( 361 )
  W( 79 ) = W( 79 ) + a*JVS( 362 )
  W( 80 ) = W( 80 ) + a*JVS( 363 )
  a = -W( 52 ) / JVS(          366  )
  W( 52 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 367 )
  W( 73 ) = W( 73 ) + a*JVS( 368 )
  W( 74 ) = W( 74 ) + a*JVS( 369 )
  W( 75 ) = W( 75 ) + a*JVS( 370 )
  W( 76 ) = W( 76 ) + a*JVS( 371 )
  W( 77 ) = W( 77 ) + a*JVS( 372 )
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 80 ) = W( 80 ) + a*JVS( 374 )
  a = -W( 53 ) / JVS(          376  )
  W( 53 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 377 )
  W( 73 ) = W( 73 ) + a*JVS( 378 )
  W( 74 ) = W( 74 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 79 ) = W( 79 ) + a*JVS( 382 )
  a = -W( 54 ) / JVS(          384  )
  W( 54 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 385 )
  W( 73 ) = W( 73 ) + a*JVS( 386 )
  W( 74 ) = W( 74 ) + a*JVS( 387 )
  W( 75 ) = W( 75 ) + a*JVS( 388 )
  W( 76 ) = W( 76 ) + a*JVS( 389 )
  W( 77 ) = W( 77 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 55 ) / JVS(          393  )
  W( 55 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  W( 73 ) = W( 73 ) + a*JVS( 395 )
  W( 74 ) = W( 74 ) + a*JVS( 396 )
  W( 75 ) = W( 75 ) + a*JVS( 397 )
  W( 76 ) = W( 76 ) + a*JVS( 398 )
  W( 79 ) = W( 79 ) + a*JVS( 399 )
  a = -W( 57 ) / JVS(          446  )
  W( 57 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 447 )
  W( 73 ) = W( 73 ) + a*JVS( 448 )
  W( 74 ) = W( 74 ) + a*JVS( 449 )
  W( 75 ) = W( 75 ) + a*JVS( 450 )
  W( 76 ) = W( 76 ) + a*JVS( 451 )
  W( 77 ) = W( 77 ) + a*JVS( 452 )
  W( 79 ) = W( 79 ) + a*JVS( 453 )
  W( 80 ) = W( 80 ) + a*JVS( 454 )
  a = -W( 58 ) / JVS(          456  )
  W( 58 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 457 )
  W( 73 ) = W( 73 ) + a*JVS( 458 )
  W( 74 ) = W( 74 ) + a*JVS( 459 )
  W( 75 ) = W( 75 ) + a*JVS( 460 )
  W( 76 ) = W( 76 ) + a*JVS( 461 )
  W( 77 ) = W( 77 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  a = -W( 59 ) / JVS(          466  )
  W( 59 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 467 )
  W( 73 ) = W( 73 ) + a*JVS( 468 )
  W( 74 ) = W( 74 ) + a*JVS( 469 )
  W( 75 ) = W( 75 ) + a*JVS( 470 )
  W( 76 ) = W( 76 ) + a*JVS( 471 )
  W( 77 ) = W( 77 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  a = -W( 60 ) / JVS(          483  )
  W( 60 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 484 )
  W( 67 ) = W( 67 ) + a*JVS( 485 )
  W( 71 ) = W( 71 ) + a*JVS( 486 )
  W( 72 ) = W( 72 ) + a*JVS( 487 )
  W( 73 ) = W( 73 ) + a*JVS( 488 )
  W( 74 ) = W( 74 ) + a*JVS( 489 )
  W( 75 ) = W( 75 ) + a*JVS( 490 )
  W( 76 ) = W( 76 ) + a*JVS( 491 )
  W( 77 ) = W( 77 ) + a*JVS( 492 )
  W( 78 ) = W( 78 ) + a*JVS( 493 )
  W( 79 ) = W( 79 ) + a*JVS( 494 )
  W( 80 ) = W( 80 ) + a*JVS( 495 )
  a = -W( 61 ) / JVS(          498  )
  W( 61 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 499 )
  W( 73 ) = W( 73 ) + a*JVS( 500 )
  W( 74 ) = W( 74 ) + a*JVS( 501 )
  W( 75 ) = W( 75 ) + a*JVS( 502 )
  W( 76 ) = W( 76 ) + a*JVS( 503 )
  W( 77 ) = W( 77 ) + a*JVS( 504 )
  W( 79 ) = W( 79 ) + a*JVS( 505 )
  W( 80 ) = W( 80 ) + a*JVS( 506 )
  a = -W( 62 ) / JVS(          517  )
  W( 62 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 518 )
  W( 70 ) = W( 70 ) + a*JVS( 519 )
  W( 71 ) = W( 71 ) + a*JVS( 520 )
  W( 72 ) = W( 72 ) + a*JVS( 521 )
  W( 73 ) = W( 73 ) + a*JVS( 522 )
  W( 74 ) = W( 74 ) + a*JVS( 523 )
  W( 75 ) = W( 75 ) + a*JVS( 524 )
  W( 76 ) = W( 76 ) + a*JVS( 525 )
  W( 77 ) = W( 77 ) + a*JVS( 526 )
  W( 78 ) = W( 78 ) + a*JVS( 527 )
  W( 79 ) = W( 79 ) + a*JVS( 528 )
  W( 80 ) = W( 80 ) + a*JVS( 529 )
  a = -W( 63 ) / JVS(          533  )
  W( 63 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 534 )
  W( 73 ) = W( 73 ) + a*JVS( 535 )
  W( 74 ) = W( 74 ) + a*JVS( 536 )
  W( 75 ) = W( 75 ) + a*JVS( 537 )
  W( 76 ) = W( 76 ) + a*JVS( 538 )
  W( 77 ) = W( 77 ) + a*JVS( 539 )
  W( 79 ) = W( 79 ) + a*JVS( 540 )
  a = -W( 64 ) / JVS(          554  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 555 )
  W( 66 ) = W( 66 ) + a*JVS( 556 )
  W( 67 ) = W( 67 ) + a*JVS( 557 )
  W( 68 ) = W( 68 ) + a*JVS( 558 )
  W( 69 ) = W( 69 ) + a*JVS( 559 )
  W( 70 ) = W( 70 ) + a*JVS( 560 )
  W( 71 ) = W( 71 ) + a*JVS( 561 )
  W( 72 ) = W( 72 ) + a*JVS( 562 )
  W( 73 ) = W( 73 ) + a*JVS( 563 )
  W( 74 ) = W( 74 ) + a*JVS( 564 )
  W( 75 ) = W( 75 ) + a*JVS( 565 )
  W( 76 ) = W( 76 ) + a*JVS( 566 )
  W( 77 ) = W( 77 ) + a*JVS( 567 )
  W( 78 ) = W( 78 ) + a*JVS( 568 )
  W( 79 ) = W( 79 ) + a*JVS( 569 )
  W( 80 ) = W( 80 ) + a*JVS( 570 )
  a = -W( 65 ) / JVS(          573  )
  W( 65 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 574 )
  W( 73 ) = W( 73 ) + a*JVS( 575 )
  W( 74 ) = W( 74 ) + a*JVS( 576 )
  W( 75 ) = W( 75 ) + a*JVS( 577 )
  W( 76 ) = W( 76 ) + a*JVS( 578 )
  W( 77 ) = W( 77 ) + a*JVS( 579 )
  W( 79 ) = W( 79 ) + a*JVS( 580 )
  a = -W( 66 ) / JVS(          587  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 588 )
  W( 72 ) = W( 72 ) + a*JVS( 589 )
  W( 73 ) = W( 73 ) + a*JVS( 590 )
  W( 74 ) = W( 74 ) + a*JVS( 591 )
  W( 75 ) = W( 75 ) + a*JVS( 592 )
  W( 76 ) = W( 76 ) + a*JVS( 593 )
  W( 77 ) = W( 77 ) + a*JVS( 594 )
  W( 79 ) = W( 79 ) + a*JVS( 595 )
  a = -W( 67 ) / JVS(          603  )
  W( 67 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 604 )
  W( 73 ) = W( 73 ) + a*JVS( 605 )
  W( 74 ) = W( 74 ) + a*JVS( 606 )
  W( 75 ) = W( 75 ) + a*JVS( 607 )
  W( 76 ) = W( 76 ) + a*JVS( 608 )
  W( 77 ) = W( 77 ) + a*JVS( 609 )
  W( 79 ) = W( 79 ) + a*JVS( 610 )
  a = -W( 68 ) / JVS(          624  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 625 )
  W( 70 ) = W( 70 ) + a*JVS( 626 )
  W( 71 ) = W( 71 ) + a*JVS( 627 )
  W( 72 ) = W( 72 ) + a*JVS( 628 )
  W( 73 ) = W( 73 ) + a*JVS( 629 )
  W( 74 ) = W( 74 ) + a*JVS( 630 )
  W( 75 ) = W( 75 ) + a*JVS( 631 )
  W( 76 ) = W( 76 ) + a*JVS( 632 )
  W( 77 ) = W( 77 ) + a*JVS( 633 )
  W( 78 ) = W( 78 ) + a*JVS( 634 )
  W( 79 ) = W( 79 ) + a*JVS( 635 )
  W( 80 ) = W( 80 ) + a*JVS( 636 )
  a = -W( 69 ) / JVS(          650  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 651 )
  W( 71 ) = W( 71 ) + a*JVS( 652 )
  W( 72 ) = W( 72 ) + a*JVS( 653 )
  W( 73 ) = W( 73 ) + a*JVS( 654 )
  W( 74 ) = W( 74 ) + a*JVS( 655 )
  W( 75 ) = W( 75 ) + a*JVS( 656 )
  W( 76 ) = W( 76 ) + a*JVS( 657 )
  W( 77 ) = W( 77 ) + a*JVS( 658 )
  W( 78 ) = W( 78 ) + a*JVS( 659 )
  W( 79 ) = W( 79 ) + a*JVS( 660 )
  W( 80 ) = W( 80 ) + a*JVS( 661 )
  a = -W( 70 ) / JVS(          677  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 678 )
  W( 72 ) = W( 72 ) + a*JVS( 679 )
  W( 73 ) = W( 73 ) + a*JVS( 680 )
  W( 74 ) = W( 74 ) + a*JVS( 681 )
  W( 75 ) = W( 75 ) + a*JVS( 682 )
  W( 76 ) = W( 76 ) + a*JVS( 683 )
  W( 77 ) = W( 77 ) + a*JVS( 684 )
  W( 78 ) = W( 78 ) + a*JVS( 685 )
  W( 79 ) = W( 79 ) + a*JVS( 686 )
  W( 80 ) = W( 80 ) + a*JVS( 687 )
  a = -W( 71 ) / JVS(          693  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 694 )
  W( 73 ) = W( 73 ) + a*JVS( 695 )
  W( 74 ) = W( 74 ) + a*JVS( 696 )
  W( 75 ) = W( 75 ) + a*JVS( 697 )
  W( 76 ) = W( 76 ) + a*JVS( 698 )
  W( 77 ) = W( 77 ) + a*JVS( 699 )
  W( 78 ) = W( 78 ) + a*JVS( 700 )
  W( 79 ) = W( 79 ) + a*JVS( 701 )
  W( 80 ) = W( 80 ) + a*JVS( 702 )
  JVS( 703) = W( 25 )
  JVS( 704) = W( 33 )
  JVS( 705) = W( 35 )
  JVS( 706) = W( 39 )
  JVS( 707) = W( 41 )
  JVS( 708) = W( 42 )
  JVS( 709) = W( 43 )
  JVS( 710) = W( 45 )
  JVS( 711) = W( 46 )
  JVS( 712) = W( 47 )
  JVS( 713) = W( 48 )
  JVS( 714) = W( 49 )
  JVS( 715) = W( 50 )
  JVS( 716) = W( 51 )
  JVS( 717) = W( 52 )
  JVS( 718) = W( 53 )
  JVS( 719) = W( 54 )
  JVS( 720) = W( 55 )
  JVS( 721) = W( 57 )
  JVS( 722) = W( 58 )
  JVS( 723) = W( 59 )
  JVS( 724) = W( 60 )
  JVS( 725) = W( 61 )
  JVS( 726) = W( 62 )
  JVS( 727) = W( 63 )
  JVS( 728) = W( 64 )
  JVS( 729) = W( 65 )
  JVS( 730) = W( 66 )
  JVS( 731) = W( 67 )
  JVS( 732) = W( 68 )
  JVS( 733) = W( 69 )
  JVS( 734) = W( 70 )
  JVS( 735) = W( 71 )
  JVS( 736) = W( 72 )
  JVS( 737) = W( 73 )
  JVS( 738) = W( 74 )
  JVS( 739) = W( 75 )
  JVS( 740) = W( 76 )
  JVS( 741) = W( 77 )
  JVS( 742) = W( 78 )
  JVS( 743) = W( 79 )
  JVS( 744) = W( 80 )
  IF ( ABS(  JVS( 806 )) < TINY(a) ) THEN
         IER = 73                                      
         RETURN
  END IF
   W( 5 ) = JVS( 745 )
   W( 6 ) = JVS( 746 )
   W( 8 ) = JVS( 747 )
   W( 9 ) = JVS( 748 )
   W( 10 ) = JVS( 749 )
   W( 12 ) = JVS( 750 )
   W( 13 ) = JVS( 751 )
   W( 14 ) = JVS( 752 )
   W( 15 ) = JVS( 753 )
   W( 16 ) = JVS( 754 )
   W( 17 ) = JVS( 755 )
   W( 18 ) = JVS( 756 )
   W( 19 ) = JVS( 757 )
   W( 20 ) = JVS( 758 )
   W( 21 ) = JVS( 759 )
   W( 22 ) = JVS( 760 )
   W( 23 ) = JVS( 761 )
   W( 24 ) = JVS( 762 )
   W( 25 ) = JVS( 763 )
   W( 28 ) = JVS( 764 )
   W( 29 ) = JVS( 765 )
   W( 30 ) = JVS( 766 )
   W( 31 ) = JVS( 767 )
   W( 32 ) = JVS( 768 )
   W( 33 ) = JVS( 769 )
   W( 34 ) = JVS( 770 )
   W( 35 ) = JVS( 771 )
   W( 36 ) = JVS( 772 )
   W( 37 ) = JVS( 773 )
   W( 38 ) = JVS( 774 )
   W( 39 ) = JVS( 775 )
   W( 40 ) = JVS( 776 )
   W( 41 ) = JVS( 777 )
   W( 42 ) = JVS( 778 )
   W( 43 ) = JVS( 779 )
   W( 44 ) = JVS( 780 )
   W( 48 ) = JVS( 781 )
   W( 49 ) = JVS( 782 )
   W( 50 ) = JVS( 783 )
   W( 51 ) = JVS( 784 )
   W( 52 ) = JVS( 785 )
   W( 53 ) = JVS( 786 )
   W( 54 ) = JVS( 787 )
   W( 55 ) = JVS( 788 )
   W( 56 ) = JVS( 789 )
   W( 57 ) = JVS( 790 )
   W( 58 ) = JVS( 791 )
   W( 59 ) = JVS( 792 )
   W( 60 ) = JVS( 793 )
   W( 61 ) = JVS( 794 )
   W( 62 ) = JVS( 795 )
   W( 63 ) = JVS( 796 )
   W( 64 ) = JVS( 797 )
   W( 65 ) = JVS( 798 )
   W( 66 ) = JVS( 799 )
   W( 67 ) = JVS( 800 )
   W( 68 ) = JVS( 801 )
   W( 69 ) = JVS( 802 )
   W( 70 ) = JVS( 803 )
   W( 71 ) = JVS( 804 )
   W( 72 ) = JVS( 805 )
   W( 73 ) = JVS( 806 )
   W( 74 ) = JVS( 807 )
   W( 75 ) = JVS( 808 )
   W( 76 ) = JVS( 809 )
   W( 77 ) = JVS( 810 )
   W( 78 ) = JVS( 811 )
   W( 79 ) = JVS( 812 )
   W( 80 ) = JVS( 813 )
  a = -W( 5 ) / JVS(           40  )
  W( 5 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 41 )
  a = -W( 6 ) / JVS(           42  )
  W( 6 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 43 )
  a = -W( 8 ) / JVS(           48  )
  W( 8 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 49 )
  a = -W( 9 ) / JVS(           50  )
  W( 9 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 51 )
  a = -W( 10 ) / JVS(           52  )
  W( 10 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 53 )
  a = -W( 12 ) / JVS(           57  )
  W( 12 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 58 )
  a = -W( 13 ) / JVS(           59  )
  W( 13 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 60 )
  W( 73 ) = W( 73 ) + a*JVS( 61 )
  W( 75 ) = W( 75 ) + a*JVS( 62 )
  a = -W( 14 ) / JVS(           63  )
  W( 14 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 64 )
  a = -W( 15 ) / JVS(           65  )
  W( 15 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 66 )
  W( 65 ) = W( 65 ) + a*JVS( 67 )
  W( 73 ) = W( 73 ) + a*JVS( 68 )
  W( 77 ) = W( 77 ) + a*JVS( 69 )
  a = -W( 16 ) / JVS(           70  )
  W( 16 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 71 )
  W( 73 ) = W( 73 ) + a*JVS( 72 )
  a = -W( 17 ) / JVS(           73  )
  W( 17 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 74 )
  W( 65 ) = W( 65 ) + a*JVS( 75 )
  W( 73 ) = W( 73 ) + a*JVS( 76 )
  W( 77 ) = W( 77 ) + a*JVS( 77 )
  a = -W( 18 ) / JVS(           78  )
  W( 18 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 79 )
  W( 73 ) = W( 73 ) + a*JVS( 80 )
  a = -W( 19 ) / JVS(           81  )
  W( 19 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 82 )
  W( 75 ) = W( 75 ) + a*JVS( 83 )
  W( 80 ) = W( 80 ) + a*JVS( 84 )
  a = -W( 20 ) / JVS(           85  )
  W( 20 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 86 )
  W( 75 ) = W( 75 ) + a*JVS( 87 )
  W( 76 ) = W( 76 ) + a*JVS( 88 )
  a = -W( 21 ) / JVS(           89  )
  W( 21 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 90 )
  W( 73 ) = W( 73 ) + a*JVS( 91 )
  W( 80 ) = W( 80 ) + a*JVS( 92 )
  a = -W( 22 ) / JVS(           94  )
  W( 22 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 95 )
  W( 42 ) = W( 42 ) + a*JVS( 96 )
  W( 73 ) = W( 73 ) + a*JVS( 97 )
  W( 79 ) = W( 79 ) + a*JVS( 98 )
  W( 80 ) = W( 80 ) + a*JVS( 99 )
  a = -W( 23 ) / JVS(          100  )
  W( 23 ) = -a
  W( 28 ) = W( 28 ) + a*JVS( 101 )
  W( 29 ) = W( 29 ) + a*JVS( 102 )
  W( 31 ) = W( 31 ) + a*JVS( 103 )
  W( 73 ) = W( 73 ) + a*JVS( 104 )
  W( 79 ) = W( 79 ) + a*JVS( 105 )
  W( 80 ) = W( 80 ) + a*JVS( 106 )
  a = -W( 24 ) / JVS(          107  )
  W( 24 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 108 )
  W( 38 ) = W( 38 ) + a*JVS( 109 )
  W( 39 ) = W( 39 ) + a*JVS( 110 )
  W( 41 ) = W( 41 ) + a*JVS( 111 )
  W( 63 ) = W( 63 ) + a*JVS( 112 )
  W( 65 ) = W( 65 ) + a*JVS( 113 )
  W( 73 ) = W( 73 ) + a*JVS( 114 )
  W( 75 ) = W( 75 ) + a*JVS( 115 )
  W( 77 ) = W( 77 ) + a*JVS( 116 )
  a = -W( 25 ) / JVS(          118  )
  W( 25 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 119 )
  W( 69 ) = W( 69 ) + a*JVS( 120 )
  W( 73 ) = W( 73 ) + a*JVS( 121 )
  W( 76 ) = W( 76 ) + a*JVS( 122 )
  a = -W( 28 ) / JVS(          137  )
  W( 28 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 138 )
  W( 77 ) = W( 77 ) + a*JVS( 139 )
  W( 80 ) = W( 80 ) + a*JVS( 140 )
  a = -W( 29 ) / JVS(          142  )
  W( 29 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 143 )
  W( 77 ) = W( 77 ) + a*JVS( 144 )
  W( 80 ) = W( 80 ) + a*JVS( 145 )
  a = -W( 30 ) / JVS(          146  )
  W( 30 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 147 )
  W( 74 ) = W( 74 ) + a*JVS( 148 )
  W( 77 ) = W( 77 ) + a*JVS( 149 )
  a = -W( 31 ) / JVS(          150  )
  W( 31 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 151 )
  W( 73 ) = W( 73 ) + a*JVS( 152 )
  W( 77 ) = W( 77 ) + a*JVS( 153 )
  W( 80 ) = W( 80 ) + a*JVS( 154 )
  a = -W( 32 ) / JVS(          155  )
  W( 32 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 156 )
  W( 44 ) = W( 44 ) + a*JVS( 157 )
  W( 51 ) = W( 51 ) + a*JVS( 158 )
  W( 56 ) = W( 56 ) + a*JVS( 159 )
  W( 61 ) = W( 61 ) + a*JVS( 160 )
  W( 64 ) = W( 64 ) + a*JVS( 161 )
  W( 73 ) = W( 73 ) + a*JVS( 162 )
  W( 74 ) = W( 74 ) + a*JVS( 163 )
  W( 75 ) = W( 75 ) + a*JVS( 164 )
  W( 80 ) = W( 80 ) + a*JVS( 165 )
  a = -W( 33 ) / JVS(          166  )
  W( 33 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 167 )
  W( 72 ) = W( 72 ) + a*JVS( 168 )
  W( 73 ) = W( 73 ) + a*JVS( 169 )
  W( 75 ) = W( 75 ) + a*JVS( 170 )
  W( 77 ) = W( 77 ) + a*JVS( 171 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 173 )
  W( 50 ) = W( 50 ) + a*JVS( 174 )
  W( 73 ) = W( 73 ) + a*JVS( 175 )
  W( 74 ) = W( 74 ) + a*JVS( 176 )
  W( 79 ) = W( 79 ) + a*JVS( 177 )
  a = -W( 35 ) / JVS(          178  )
  W( 35 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 179 )
  W( 72 ) = W( 72 ) + a*JVS( 180 )
  W( 73 ) = W( 73 ) + a*JVS( 181 )
  W( 74 ) = W( 74 ) + a*JVS( 182 )
  W( 77 ) = W( 77 ) + a*JVS( 183 )
  W( 80 ) = W( 80 ) + a*JVS( 184 )
  a = -W( 36 ) / JVS(          185  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 186 )
  W( 74 ) = W( 74 ) + a*JVS( 187 )
  W( 77 ) = W( 77 ) + a*JVS( 188 )
  a = -W( 37 ) / JVS(          194  )
  W( 37 ) = -a
  W( 38 ) = W( 38 ) + a*JVS( 195 )
  W( 39 ) = W( 39 ) + a*JVS( 196 )
  W( 41 ) = W( 41 ) + a*JVS( 197 )
  W( 42 ) = W( 42 ) + a*JVS( 198 )
  W( 43 ) = W( 43 ) + a*JVS( 199 )
  W( 44 ) = W( 44 ) + a*JVS( 200 )
  W( 48 ) = W( 48 ) + a*JVS( 201 )
  W( 50 ) = W( 50 ) + a*JVS( 202 )
  W( 51 ) = W( 51 ) + a*JVS( 203 )
  W( 56 ) = W( 56 ) + a*JVS( 204 )
  W( 61 ) = W( 61 ) + a*JVS( 205 )
  W( 63 ) = W( 63 ) + a*JVS( 206 )
  W( 64 ) = W( 64 ) + a*JVS( 207 )
  W( 65 ) = W( 65 ) + a*JVS( 208 )
  W( 73 ) = W( 73 ) + a*JVS( 209 )
  W( 74 ) = W( 74 ) + a*JVS( 210 )
  W( 77 ) = W( 77 ) + a*JVS( 211 )
  W( 79 ) = W( 79 ) + a*JVS( 212 )
  a = -W( 38 ) / JVS(          213  )
  W( 38 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 214 )
  W( 74 ) = W( 74 ) + a*JVS( 215 )
  W( 77 ) = W( 77 ) + a*JVS( 216 )
  a = -W( 39 ) / JVS(          217  )
  W( 39 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 218 )
  W( 74 ) = W( 74 ) + a*JVS( 219 )
  W( 77 ) = W( 77 ) + a*JVS( 220 )
  a = -W( 40 ) / JVS(          225  )
  W( 40 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 226 )
  W( 74 ) = W( 74 ) + a*JVS( 227 )
  W( 75 ) = W( 75 ) + a*JVS( 228 )
  W( 77 ) = W( 77 ) + a*JVS( 229 )
  W( 80 ) = W( 80 ) + a*JVS( 230 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  a = -W( 42 ) / JVS(          238  )
  W( 42 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 239 )
  W( 73 ) = W( 73 ) + a*JVS( 240 )
  W( 74 ) = W( 74 ) + a*JVS( 241 )
  W( 75 ) = W( 75 ) + a*JVS( 242 )
  W( 77 ) = W( 77 ) + a*JVS( 243 )
  W( 79 ) = W( 79 ) + a*JVS( 244 )
  W( 80 ) = W( 80 ) + a*JVS( 245 )
  a = -W( 43 ) / JVS(          246  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 247 )
  W( 74 ) = W( 74 ) + a*JVS( 248 )
  W( 77 ) = W( 77 ) + a*JVS( 249 )
  W( 78 ) = W( 78 ) + a*JVS( 250 )
  W( 80 ) = W( 80 ) + a*JVS( 251 )
  a = -W( 44 ) / JVS(          253  )
  W( 44 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 254 )
  W( 57 ) = W( 57 ) + a*JVS( 255 )
  W( 58 ) = W( 58 ) + a*JVS( 256 )
  W( 61 ) = W( 61 ) + a*JVS( 257 )
  W( 71 ) = W( 71 ) + a*JVS( 258 )
  W( 72 ) = W( 72 ) + a*JVS( 259 )
  W( 73 ) = W( 73 ) + a*JVS( 260 )
  W( 74 ) = W( 74 ) + a*JVS( 261 )
  W( 76 ) = W( 76 ) + a*JVS( 262 )
  W( 77 ) = W( 77 ) + a*JVS( 263 )
  W( 79 ) = W( 79 ) + a*JVS( 264 )
  a = -W( 48 ) / JVS(          294  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 295 )
  W( 54 ) = W( 54 ) + a*JVS( 296 )
  W( 72 ) = W( 72 ) + a*JVS( 297 )
  W( 73 ) = W( 73 ) + a*JVS( 298 )
  W( 74 ) = W( 74 ) + a*JVS( 299 )
  W( 75 ) = W( 75 ) + a*JVS( 300 )
  W( 76 ) = W( 76 ) + a*JVS( 301 )
  W( 77 ) = W( 77 ) + a*JVS( 302 )
  W( 79 ) = W( 79 ) + a*JVS( 303 )
  a = -W( 49 ) / JVS(          309  )
  W( 49 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 310 )
  W( 55 ) = W( 55 ) + a*JVS( 311 )
  W( 59 ) = W( 59 ) + a*JVS( 312 )
  W( 60 ) = W( 60 ) + a*JVS( 313 )
  W( 61 ) = W( 61 ) + a*JVS( 314 )
  W( 63 ) = W( 63 ) + a*JVS( 315 )
  W( 65 ) = W( 65 ) + a*JVS( 316 )
  W( 66 ) = W( 66 ) + a*JVS( 317 )
  W( 67 ) = W( 67 ) + a*JVS( 318 )
  W( 69 ) = W( 69 ) + a*JVS( 319 )
  W( 70 ) = W( 70 ) + a*JVS( 320 )
  W( 71 ) = W( 71 ) + a*JVS( 321 )
  W( 72 ) = W( 72 ) + a*JVS( 322 )
  W( 73 ) = W( 73 ) + a*JVS( 323 )
  W( 74 ) = W( 74 ) + a*JVS( 324 )
  W( 75 ) = W( 75 ) + a*JVS( 325 )
  W( 76 ) = W( 76 ) + a*JVS( 326 )
  W( 77 ) = W( 77 ) + a*JVS( 327 )
  W( 79 ) = W( 79 ) + a*JVS( 328 )
  a = -W( 50 ) / JVS(          331  )
  W( 50 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 332 )
  W( 73 ) = W( 73 ) + a*JVS( 333 )
  W( 74 ) = W( 74 ) + a*JVS( 334 )
  W( 75 ) = W( 75 ) + a*JVS( 335 )
  W( 76 ) = W( 76 ) + a*JVS( 336 )
  W( 77 ) = W( 77 ) + a*JVS( 337 )
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  a = -W( 51 ) / JVS(          347  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 348 )
  W( 54 ) = W( 54 ) + a*JVS( 349 )
  W( 57 ) = W( 57 ) + a*JVS( 350 )
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 61 ) = W( 61 ) + a*JVS( 352 )
  W( 69 ) = W( 69 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 72 ) = W( 72 ) + a*JVS( 355 )
  W( 73 ) = W( 73 ) + a*JVS( 356 )
  W( 74 ) = W( 74 ) + a*JVS( 357 )
  W( 75 ) = W( 75 ) + a*JVS( 358 )
  W( 76 ) = W( 76 ) + a*JVS( 359 )
  W( 77 ) = W( 77 ) + a*JVS( 360 )
  W( 78 ) = W( 78 ) + a*JVS( 361 )
  W( 79 ) = W( 79 ) + a*JVS( 362 )
  W( 80 ) = W( 80 ) + a*JVS( 363 )
  a = -W( 52 ) / JVS(          366  )
  W( 52 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 367 )
  W( 73 ) = W( 73 ) + a*JVS( 368 )
  W( 74 ) = W( 74 ) + a*JVS( 369 )
  W( 75 ) = W( 75 ) + a*JVS( 370 )
  W( 76 ) = W( 76 ) + a*JVS( 371 )
  W( 77 ) = W( 77 ) + a*JVS( 372 )
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 80 ) = W( 80 ) + a*JVS( 374 )
  a = -W( 53 ) / JVS(          376  )
  W( 53 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 377 )
  W( 73 ) = W( 73 ) + a*JVS( 378 )
  W( 74 ) = W( 74 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 79 ) = W( 79 ) + a*JVS( 382 )
  a = -W( 54 ) / JVS(          384  )
  W( 54 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 385 )
  W( 73 ) = W( 73 ) + a*JVS( 386 )
  W( 74 ) = W( 74 ) + a*JVS( 387 )
  W( 75 ) = W( 75 ) + a*JVS( 388 )
  W( 76 ) = W( 76 ) + a*JVS( 389 )
  W( 77 ) = W( 77 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 55 ) / JVS(          393  )
  W( 55 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  W( 73 ) = W( 73 ) + a*JVS( 395 )
  W( 74 ) = W( 74 ) + a*JVS( 396 )
  W( 75 ) = W( 75 ) + a*JVS( 397 )
  W( 76 ) = W( 76 ) + a*JVS( 398 )
  W( 79 ) = W( 79 ) + a*JVS( 399 )
  a = -W( 56 ) / JVS(          423  )
  W( 56 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 424 )
  W( 58 ) = W( 58 ) + a*JVS( 425 )
  W( 59 ) = W( 59 ) + a*JVS( 426 )
  W( 61 ) = W( 61 ) + a*JVS( 427 )
  W( 62 ) = W( 62 ) + a*JVS( 428 )
  W( 63 ) = W( 63 ) + a*JVS( 429 )
  W( 65 ) = W( 65 ) + a*JVS( 430 )
  W( 66 ) = W( 66 ) + a*JVS( 431 )
  W( 67 ) = W( 67 ) + a*JVS( 432 )
  W( 68 ) = W( 68 ) + a*JVS( 433 )
  W( 69 ) = W( 69 ) + a*JVS( 434 )
  W( 71 ) = W( 71 ) + a*JVS( 435 )
  W( 72 ) = W( 72 ) + a*JVS( 436 )
  W( 73 ) = W( 73 ) + a*JVS( 437 )
  W( 74 ) = W( 74 ) + a*JVS( 438 )
  W( 75 ) = W( 75 ) + a*JVS( 439 )
  W( 76 ) = W( 76 ) + a*JVS( 440 )
  W( 77 ) = W( 77 ) + a*JVS( 441 )
  W( 78 ) = W( 78 ) + a*JVS( 442 )
  W( 79 ) = W( 79 ) + a*JVS( 443 )
  W( 80 ) = W( 80 ) + a*JVS( 444 )
  a = -W( 57 ) / JVS(          446  )
  W( 57 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 447 )
  W( 73 ) = W( 73 ) + a*JVS( 448 )
  W( 74 ) = W( 74 ) + a*JVS( 449 )
  W( 75 ) = W( 75 ) + a*JVS( 450 )
  W( 76 ) = W( 76 ) + a*JVS( 451 )
  W( 77 ) = W( 77 ) + a*JVS( 452 )
  W( 79 ) = W( 79 ) + a*JVS( 453 )
  W( 80 ) = W( 80 ) + a*JVS( 454 )
  a = -W( 58 ) / JVS(          456  )
  W( 58 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 457 )
  W( 73 ) = W( 73 ) + a*JVS( 458 )
  W( 74 ) = W( 74 ) + a*JVS( 459 )
  W( 75 ) = W( 75 ) + a*JVS( 460 )
  W( 76 ) = W( 76 ) + a*JVS( 461 )
  W( 77 ) = W( 77 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  a = -W( 59 ) / JVS(          466  )
  W( 59 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 467 )
  W( 73 ) = W( 73 ) + a*JVS( 468 )
  W( 74 ) = W( 74 ) + a*JVS( 469 )
  W( 75 ) = W( 75 ) + a*JVS( 470 )
  W( 76 ) = W( 76 ) + a*JVS( 471 )
  W( 77 ) = W( 77 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  a = -W( 60 ) / JVS(          483  )
  W( 60 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 484 )
  W( 67 ) = W( 67 ) + a*JVS( 485 )
  W( 71 ) = W( 71 ) + a*JVS( 486 )
  W( 72 ) = W( 72 ) + a*JVS( 487 )
  W( 73 ) = W( 73 ) + a*JVS( 488 )
  W( 74 ) = W( 74 ) + a*JVS( 489 )
  W( 75 ) = W( 75 ) + a*JVS( 490 )
  W( 76 ) = W( 76 ) + a*JVS( 491 )
  W( 77 ) = W( 77 ) + a*JVS( 492 )
  W( 78 ) = W( 78 ) + a*JVS( 493 )
  W( 79 ) = W( 79 ) + a*JVS( 494 )
  W( 80 ) = W( 80 ) + a*JVS( 495 )
  a = -W( 61 ) / JVS(          498  )
  W( 61 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 499 )
  W( 73 ) = W( 73 ) + a*JVS( 500 )
  W( 74 ) = W( 74 ) + a*JVS( 501 )
  W( 75 ) = W( 75 ) + a*JVS( 502 )
  W( 76 ) = W( 76 ) + a*JVS( 503 )
  W( 77 ) = W( 77 ) + a*JVS( 504 )
  W( 79 ) = W( 79 ) + a*JVS( 505 )
  W( 80 ) = W( 80 ) + a*JVS( 506 )
  a = -W( 62 ) / JVS(          517  )
  W( 62 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 518 )
  W( 70 ) = W( 70 ) + a*JVS( 519 )
  W( 71 ) = W( 71 ) + a*JVS( 520 )
  W( 72 ) = W( 72 ) + a*JVS( 521 )
  W( 73 ) = W( 73 ) + a*JVS( 522 )
  W( 74 ) = W( 74 ) + a*JVS( 523 )
  W( 75 ) = W( 75 ) + a*JVS( 524 )
  W( 76 ) = W( 76 ) + a*JVS( 525 )
  W( 77 ) = W( 77 ) + a*JVS( 526 )
  W( 78 ) = W( 78 ) + a*JVS( 527 )
  W( 79 ) = W( 79 ) + a*JVS( 528 )
  W( 80 ) = W( 80 ) + a*JVS( 529 )
  a = -W( 63 ) / JVS(          533  )
  W( 63 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 534 )
  W( 73 ) = W( 73 ) + a*JVS( 535 )
  W( 74 ) = W( 74 ) + a*JVS( 536 )
  W( 75 ) = W( 75 ) + a*JVS( 537 )
  W( 76 ) = W( 76 ) + a*JVS( 538 )
  W( 77 ) = W( 77 ) + a*JVS( 539 )
  W( 79 ) = W( 79 ) + a*JVS( 540 )
  a = -W( 64 ) / JVS(          554  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 555 )
  W( 66 ) = W( 66 ) + a*JVS( 556 )
  W( 67 ) = W( 67 ) + a*JVS( 557 )
  W( 68 ) = W( 68 ) + a*JVS( 558 )
  W( 69 ) = W( 69 ) + a*JVS( 559 )
  W( 70 ) = W( 70 ) + a*JVS( 560 )
  W( 71 ) = W( 71 ) + a*JVS( 561 )
  W( 72 ) = W( 72 ) + a*JVS( 562 )
  W( 73 ) = W( 73 ) + a*JVS( 563 )
  W( 74 ) = W( 74 ) + a*JVS( 564 )
  W( 75 ) = W( 75 ) + a*JVS( 565 )
  W( 76 ) = W( 76 ) + a*JVS( 566 )
  W( 77 ) = W( 77 ) + a*JVS( 567 )
  W( 78 ) = W( 78 ) + a*JVS( 568 )
  W( 79 ) = W( 79 ) + a*JVS( 569 )
  W( 80 ) = W( 80 ) + a*JVS( 570 )
  a = -W( 65 ) / JVS(          573  )
  W( 65 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 574 )
  W( 73 ) = W( 73 ) + a*JVS( 575 )
  W( 74 ) = W( 74 ) + a*JVS( 576 )
  W( 75 ) = W( 75 ) + a*JVS( 577 )
  W( 76 ) = W( 76 ) + a*JVS( 578 )
  W( 77 ) = W( 77 ) + a*JVS( 579 )
  W( 79 ) = W( 79 ) + a*JVS( 580 )
  a = -W( 66 ) / JVS(          587  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 588 )
  W( 72 ) = W( 72 ) + a*JVS( 589 )
  W( 73 ) = W( 73 ) + a*JVS( 590 )
  W( 74 ) = W( 74 ) + a*JVS( 591 )
  W( 75 ) = W( 75 ) + a*JVS( 592 )
  W( 76 ) = W( 76 ) + a*JVS( 593 )
  W( 77 ) = W( 77 ) + a*JVS( 594 )
  W( 79 ) = W( 79 ) + a*JVS( 595 )
  a = -W( 67 ) / JVS(          603  )
  W( 67 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 604 )
  W( 73 ) = W( 73 ) + a*JVS( 605 )
  W( 74 ) = W( 74 ) + a*JVS( 606 )
  W( 75 ) = W( 75 ) + a*JVS( 607 )
  W( 76 ) = W( 76 ) + a*JVS( 608 )
  W( 77 ) = W( 77 ) + a*JVS( 609 )
  W( 79 ) = W( 79 ) + a*JVS( 610 )
  a = -W( 68 ) / JVS(          624  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 625 )
  W( 70 ) = W( 70 ) + a*JVS( 626 )
  W( 71 ) = W( 71 ) + a*JVS( 627 )
  W( 72 ) = W( 72 ) + a*JVS( 628 )
  W( 73 ) = W( 73 ) + a*JVS( 629 )
  W( 74 ) = W( 74 ) + a*JVS( 630 )
  W( 75 ) = W( 75 ) + a*JVS( 631 )
  W( 76 ) = W( 76 ) + a*JVS( 632 )
  W( 77 ) = W( 77 ) + a*JVS( 633 )
  W( 78 ) = W( 78 ) + a*JVS( 634 )
  W( 79 ) = W( 79 ) + a*JVS( 635 )
  W( 80 ) = W( 80 ) + a*JVS( 636 )
  a = -W( 69 ) / JVS(          650  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 651 )
  W( 71 ) = W( 71 ) + a*JVS( 652 )
  W( 72 ) = W( 72 ) + a*JVS( 653 )
  W( 73 ) = W( 73 ) + a*JVS( 654 )
  W( 74 ) = W( 74 ) + a*JVS( 655 )
  W( 75 ) = W( 75 ) + a*JVS( 656 )
  W( 76 ) = W( 76 ) + a*JVS( 657 )
  W( 77 ) = W( 77 ) + a*JVS( 658 )
  W( 78 ) = W( 78 ) + a*JVS( 659 )
  W( 79 ) = W( 79 ) + a*JVS( 660 )
  W( 80 ) = W( 80 ) + a*JVS( 661 )
  a = -W( 70 ) / JVS(          677  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 678 )
  W( 72 ) = W( 72 ) + a*JVS( 679 )
  W( 73 ) = W( 73 ) + a*JVS( 680 )
  W( 74 ) = W( 74 ) + a*JVS( 681 )
  W( 75 ) = W( 75 ) + a*JVS( 682 )
  W( 76 ) = W( 76 ) + a*JVS( 683 )
  W( 77 ) = W( 77 ) + a*JVS( 684 )
  W( 78 ) = W( 78 ) + a*JVS( 685 )
  W( 79 ) = W( 79 ) + a*JVS( 686 )
  W( 80 ) = W( 80 ) + a*JVS( 687 )
  a = -W( 71 ) / JVS(          693  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 694 )
  W( 73 ) = W( 73 ) + a*JVS( 695 )
  W( 74 ) = W( 74 ) + a*JVS( 696 )
  W( 75 ) = W( 75 ) + a*JVS( 697 )
  W( 76 ) = W( 76 ) + a*JVS( 698 )
  W( 77 ) = W( 77 ) + a*JVS( 699 )
  W( 78 ) = W( 78 ) + a*JVS( 700 )
  W( 79 ) = W( 79 ) + a*JVS( 701 )
  W( 80 ) = W( 80 ) + a*JVS( 702 )
  a = -W( 72 ) / JVS(          736  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 737 )
  W( 74 ) = W( 74 ) + a*JVS( 738 )
  W( 75 ) = W( 75 ) + a*JVS( 739 )
  W( 76 ) = W( 76 ) + a*JVS( 740 )
  W( 77 ) = W( 77 ) + a*JVS( 741 )
  W( 78 ) = W( 78 ) + a*JVS( 742 )
  W( 79 ) = W( 79 ) + a*JVS( 743 )
  W( 80 ) = W( 80 ) + a*JVS( 744 )
  JVS( 745) = W( 5 )
  JVS( 746) = W( 6 )
  JVS( 747) = W( 8 )
  JVS( 748) = W( 9 )
  JVS( 749) = W( 10 )
  JVS( 750) = W( 12 )
  JVS( 751) = W( 13 )
  JVS( 752) = W( 14 )
  JVS( 753) = W( 15 )
  JVS( 754) = W( 16 )
  JVS( 755) = W( 17 )
  JVS( 756) = W( 18 )
  JVS( 757) = W( 19 )
  JVS( 758) = W( 20 )
  JVS( 759) = W( 21 )
  JVS( 760) = W( 22 )
  JVS( 761) = W( 23 )
  JVS( 762) = W( 24 )
  JVS( 763) = W( 25 )
  JVS( 764) = W( 28 )
  JVS( 765) = W( 29 )
  JVS( 766) = W( 30 )
  JVS( 767) = W( 31 )
  JVS( 768) = W( 32 )
  JVS( 769) = W( 33 )
  JVS( 770) = W( 34 )
  JVS( 771) = W( 35 )
  JVS( 772) = W( 36 )
  JVS( 773) = W( 37 )
  JVS( 774) = W( 38 )
  JVS( 775) = W( 39 )
  JVS( 776) = W( 40 )
  JVS( 777) = W( 41 )
  JVS( 778) = W( 42 )
  JVS( 779) = W( 43 )
  JVS( 780) = W( 44 )
  JVS( 781) = W( 48 )
  JVS( 782) = W( 49 )
  JVS( 783) = W( 50 )
  JVS( 784) = W( 51 )
  JVS( 785) = W( 52 )
  JVS( 786) = W( 53 )
  JVS( 787) = W( 54 )
  JVS( 788) = W( 55 )
  JVS( 789) = W( 56 )
  JVS( 790) = W( 57 )
  JVS( 791) = W( 58 )
  JVS( 792) = W( 59 )
  JVS( 793) = W( 60 )
  JVS( 794) = W( 61 )
  JVS( 795) = W( 62 )
  JVS( 796) = W( 63 )
  JVS( 797) = W( 64 )
  JVS( 798) = W( 65 )
  JVS( 799) = W( 66 )
  JVS( 800) = W( 67 )
  JVS( 801) = W( 68 )
  JVS( 802) = W( 69 )
  JVS( 803) = W( 70 )
  JVS( 804) = W( 71 )
  JVS( 805) = W( 72 )
  JVS( 806) = W( 73 )
  JVS( 807) = W( 74 )
  JVS( 808) = W( 75 )
  JVS( 809) = W( 76 )
  JVS( 810) = W( 77 )
  JVS( 811) = W( 78 )
  JVS( 812) = W( 79 )
  JVS( 813) = W( 80 )
  IF ( ABS(  JVS( 852 )) < TINY(a) ) THEN
         IER = 74                                      
         RETURN
  END IF
   W( 11 ) = JVS( 814 )
   W( 19 ) = JVS( 815 )
   W( 26 ) = JVS( 816 )
   W( 30 ) = JVS( 817 )
   W( 32 ) = JVS( 818 )
   W( 35 ) = JVS( 819 )
   W( 36 ) = JVS( 820 )
   W( 38 ) = JVS( 821 )
   W( 39 ) = JVS( 822 )
   W( 40 ) = JVS( 823 )
   W( 41 ) = JVS( 824 )
   W( 43 ) = JVS( 825 )
   W( 44 ) = JVS( 826 )
   W( 45 ) = JVS( 827 )
   W( 46 ) = JVS( 828 )
   W( 47 ) = JVS( 829 )
   W( 51 ) = JVS( 830 )
   W( 52 ) = JVS( 831 )
   W( 53 ) = JVS( 832 )
   W( 54 ) = JVS( 833 )
   W( 55 ) = JVS( 834 )
   W( 56 ) = JVS( 835 )
   W( 57 ) = JVS( 836 )
   W( 58 ) = JVS( 837 )
   W( 59 ) = JVS( 838 )
   W( 61 ) = JVS( 839 )
   W( 62 ) = JVS( 840 )
   W( 63 ) = JVS( 841 )
   W( 64 ) = JVS( 842 )
   W( 65 ) = JVS( 843 )
   W( 66 ) = JVS( 844 )
   W( 67 ) = JVS( 845 )
   W( 68 ) = JVS( 846 )
   W( 69 ) = JVS( 847 )
   W( 70 ) = JVS( 848 )
   W( 71 ) = JVS( 849 )
   W( 72 ) = JVS( 850 )
   W( 73 ) = JVS( 851 )
   W( 74 ) = JVS( 852 )
   W( 75 ) = JVS( 853 )
   W( 76 ) = JVS( 854 )
   W( 77 ) = JVS( 855 )
   W( 78 ) = JVS( 856 )
   W( 79 ) = JVS( 857 )
   W( 80 ) = JVS( 858 )
  a = -W( 11 ) / JVS(           54  )
  W( 11 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 55 )
  W( 80 ) = W( 80 ) + a*JVS( 56 )
  a = -W( 19 ) / JVS(           81  )
  W( 19 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 82 )
  W( 75 ) = W( 75 ) + a*JVS( 83 )
  W( 80 ) = W( 80 ) + a*JVS( 84 )
  a = -W( 26 ) / JVS(          124  )
  W( 26 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 125 )
  W( 74 ) = W( 74 ) + a*JVS( 126 )
  W( 77 ) = W( 77 ) + a*JVS( 127 )
  W( 79 ) = W( 79 ) + a*JVS( 128 )
  W( 80 ) = W( 80 ) + a*JVS( 129 )
  a = -W( 30 ) / JVS(          146  )
  W( 30 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 147 )
  W( 74 ) = W( 74 ) + a*JVS( 148 )
  W( 77 ) = W( 77 ) + a*JVS( 149 )
  a = -W( 32 ) / JVS(          155  )
  W( 32 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 156 )
  W( 44 ) = W( 44 ) + a*JVS( 157 )
  W( 51 ) = W( 51 ) + a*JVS( 158 )
  W( 56 ) = W( 56 ) + a*JVS( 159 )
  W( 61 ) = W( 61 ) + a*JVS( 160 )
  W( 64 ) = W( 64 ) + a*JVS( 161 )
  W( 73 ) = W( 73 ) + a*JVS( 162 )
  W( 74 ) = W( 74 ) + a*JVS( 163 )
  W( 75 ) = W( 75 ) + a*JVS( 164 )
  W( 80 ) = W( 80 ) + a*JVS( 165 )
  a = -W( 35 ) / JVS(          178  )
  W( 35 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 179 )
  W( 72 ) = W( 72 ) + a*JVS( 180 )
  W( 73 ) = W( 73 ) + a*JVS( 181 )
  W( 74 ) = W( 74 ) + a*JVS( 182 )
  W( 77 ) = W( 77 ) + a*JVS( 183 )
  W( 80 ) = W( 80 ) + a*JVS( 184 )
  a = -W( 36 ) / JVS(          185  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 186 )
  W( 74 ) = W( 74 ) + a*JVS( 187 )
  W( 77 ) = W( 77 ) + a*JVS( 188 )
  a = -W( 38 ) / JVS(          213  )
  W( 38 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 214 )
  W( 74 ) = W( 74 ) + a*JVS( 215 )
  W( 77 ) = W( 77 ) + a*JVS( 216 )
  a = -W( 39 ) / JVS(          217  )
  W( 39 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 218 )
  W( 74 ) = W( 74 ) + a*JVS( 219 )
  W( 77 ) = W( 77 ) + a*JVS( 220 )
  a = -W( 40 ) / JVS(          225  )
  W( 40 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 226 )
  W( 74 ) = W( 74 ) + a*JVS( 227 )
  W( 75 ) = W( 75 ) + a*JVS( 228 )
  W( 77 ) = W( 77 ) + a*JVS( 229 )
  W( 80 ) = W( 80 ) + a*JVS( 230 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  a = -W( 43 ) / JVS(          246  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 247 )
  W( 74 ) = W( 74 ) + a*JVS( 248 )
  W( 77 ) = W( 77 ) + a*JVS( 249 )
  W( 78 ) = W( 78 ) + a*JVS( 250 )
  W( 80 ) = W( 80 ) + a*JVS( 251 )
  a = -W( 44 ) / JVS(          253  )
  W( 44 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 254 )
  W( 57 ) = W( 57 ) + a*JVS( 255 )
  W( 58 ) = W( 58 ) + a*JVS( 256 )
  W( 61 ) = W( 61 ) + a*JVS( 257 )
  W( 71 ) = W( 71 ) + a*JVS( 258 )
  W( 72 ) = W( 72 ) + a*JVS( 259 )
  W( 73 ) = W( 73 ) + a*JVS( 260 )
  W( 74 ) = W( 74 ) + a*JVS( 261 )
  W( 76 ) = W( 76 ) + a*JVS( 262 )
  W( 77 ) = W( 77 ) + a*JVS( 263 )
  W( 79 ) = W( 79 ) + a*JVS( 264 )
  a = -W( 45 ) / JVS(          266  )
  W( 45 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 267 )
  W( 73 ) = W( 73 ) + a*JVS( 268 )
  W( 74 ) = W( 74 ) + a*JVS( 269 )
  W( 75 ) = W( 75 ) + a*JVS( 270 )
  W( 76 ) = W( 76 ) + a*JVS( 271 )
  W( 77 ) = W( 77 ) + a*JVS( 272 )
  W( 79 ) = W( 79 ) + a*JVS( 273 )
  a = -W( 46 ) / JVS(          274  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 275 )
  W( 72 ) = W( 72 ) + a*JVS( 276 )
  W( 73 ) = W( 73 ) + a*JVS( 277 )
  W( 74 ) = W( 74 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 76 ) = W( 76 ) + a*JVS( 280 )
  W( 79 ) = W( 79 ) + a*JVS( 281 )
  a = -W( 47 ) / JVS(          282  )
  W( 47 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 283 )
  W( 72 ) = W( 72 ) + a*JVS( 284 )
  W( 73 ) = W( 73 ) + a*JVS( 285 )
  W( 74 ) = W( 74 ) + a*JVS( 286 )
  W( 75 ) = W( 75 ) + a*JVS( 287 )
  W( 76 ) = W( 76 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  a = -W( 51 ) / JVS(          347  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 348 )
  W( 54 ) = W( 54 ) + a*JVS( 349 )
  W( 57 ) = W( 57 ) + a*JVS( 350 )
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 61 ) = W( 61 ) + a*JVS( 352 )
  W( 69 ) = W( 69 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 72 ) = W( 72 ) + a*JVS( 355 )
  W( 73 ) = W( 73 ) + a*JVS( 356 )
  W( 74 ) = W( 74 ) + a*JVS( 357 )
  W( 75 ) = W( 75 ) + a*JVS( 358 )
  W( 76 ) = W( 76 ) + a*JVS( 359 )
  W( 77 ) = W( 77 ) + a*JVS( 360 )
  W( 78 ) = W( 78 ) + a*JVS( 361 )
  W( 79 ) = W( 79 ) + a*JVS( 362 )
  W( 80 ) = W( 80 ) + a*JVS( 363 )
  a = -W( 52 ) / JVS(          366  )
  W( 52 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 367 )
  W( 73 ) = W( 73 ) + a*JVS( 368 )
  W( 74 ) = W( 74 ) + a*JVS( 369 )
  W( 75 ) = W( 75 ) + a*JVS( 370 )
  W( 76 ) = W( 76 ) + a*JVS( 371 )
  W( 77 ) = W( 77 ) + a*JVS( 372 )
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 80 ) = W( 80 ) + a*JVS( 374 )
  a = -W( 53 ) / JVS(          376  )
  W( 53 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 377 )
  W( 73 ) = W( 73 ) + a*JVS( 378 )
  W( 74 ) = W( 74 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 79 ) = W( 79 ) + a*JVS( 382 )
  a = -W( 54 ) / JVS(          384  )
  W( 54 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 385 )
  W( 73 ) = W( 73 ) + a*JVS( 386 )
  W( 74 ) = W( 74 ) + a*JVS( 387 )
  W( 75 ) = W( 75 ) + a*JVS( 388 )
  W( 76 ) = W( 76 ) + a*JVS( 389 )
  W( 77 ) = W( 77 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 55 ) / JVS(          393  )
  W( 55 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  W( 73 ) = W( 73 ) + a*JVS( 395 )
  W( 74 ) = W( 74 ) + a*JVS( 396 )
  W( 75 ) = W( 75 ) + a*JVS( 397 )
  W( 76 ) = W( 76 ) + a*JVS( 398 )
  W( 79 ) = W( 79 ) + a*JVS( 399 )
  a = -W( 56 ) / JVS(          423  )
  W( 56 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 424 )
  W( 58 ) = W( 58 ) + a*JVS( 425 )
  W( 59 ) = W( 59 ) + a*JVS( 426 )
  W( 61 ) = W( 61 ) + a*JVS( 427 )
  W( 62 ) = W( 62 ) + a*JVS( 428 )
  W( 63 ) = W( 63 ) + a*JVS( 429 )
  W( 65 ) = W( 65 ) + a*JVS( 430 )
  W( 66 ) = W( 66 ) + a*JVS( 431 )
  W( 67 ) = W( 67 ) + a*JVS( 432 )
  W( 68 ) = W( 68 ) + a*JVS( 433 )
  W( 69 ) = W( 69 ) + a*JVS( 434 )
  W( 71 ) = W( 71 ) + a*JVS( 435 )
  W( 72 ) = W( 72 ) + a*JVS( 436 )
  W( 73 ) = W( 73 ) + a*JVS( 437 )
  W( 74 ) = W( 74 ) + a*JVS( 438 )
  W( 75 ) = W( 75 ) + a*JVS( 439 )
  W( 76 ) = W( 76 ) + a*JVS( 440 )
  W( 77 ) = W( 77 ) + a*JVS( 441 )
  W( 78 ) = W( 78 ) + a*JVS( 442 )
  W( 79 ) = W( 79 ) + a*JVS( 443 )
  W( 80 ) = W( 80 ) + a*JVS( 444 )
  a = -W( 57 ) / JVS(          446  )
  W( 57 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 447 )
  W( 73 ) = W( 73 ) + a*JVS( 448 )
  W( 74 ) = W( 74 ) + a*JVS( 449 )
  W( 75 ) = W( 75 ) + a*JVS( 450 )
  W( 76 ) = W( 76 ) + a*JVS( 451 )
  W( 77 ) = W( 77 ) + a*JVS( 452 )
  W( 79 ) = W( 79 ) + a*JVS( 453 )
  W( 80 ) = W( 80 ) + a*JVS( 454 )
  a = -W( 58 ) / JVS(          456  )
  W( 58 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 457 )
  W( 73 ) = W( 73 ) + a*JVS( 458 )
  W( 74 ) = W( 74 ) + a*JVS( 459 )
  W( 75 ) = W( 75 ) + a*JVS( 460 )
  W( 76 ) = W( 76 ) + a*JVS( 461 )
  W( 77 ) = W( 77 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  a = -W( 59 ) / JVS(          466  )
  W( 59 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 467 )
  W( 73 ) = W( 73 ) + a*JVS( 468 )
  W( 74 ) = W( 74 ) + a*JVS( 469 )
  W( 75 ) = W( 75 ) + a*JVS( 470 )
  W( 76 ) = W( 76 ) + a*JVS( 471 )
  W( 77 ) = W( 77 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  a = -W( 61 ) / JVS(          498  )
  W( 61 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 499 )
  W( 73 ) = W( 73 ) + a*JVS( 500 )
  W( 74 ) = W( 74 ) + a*JVS( 501 )
  W( 75 ) = W( 75 ) + a*JVS( 502 )
  W( 76 ) = W( 76 ) + a*JVS( 503 )
  W( 77 ) = W( 77 ) + a*JVS( 504 )
  W( 79 ) = W( 79 ) + a*JVS( 505 )
  W( 80 ) = W( 80 ) + a*JVS( 506 )
  a = -W( 62 ) / JVS(          517  )
  W( 62 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 518 )
  W( 70 ) = W( 70 ) + a*JVS( 519 )
  W( 71 ) = W( 71 ) + a*JVS( 520 )
  W( 72 ) = W( 72 ) + a*JVS( 521 )
  W( 73 ) = W( 73 ) + a*JVS( 522 )
  W( 74 ) = W( 74 ) + a*JVS( 523 )
  W( 75 ) = W( 75 ) + a*JVS( 524 )
  W( 76 ) = W( 76 ) + a*JVS( 525 )
  W( 77 ) = W( 77 ) + a*JVS( 526 )
  W( 78 ) = W( 78 ) + a*JVS( 527 )
  W( 79 ) = W( 79 ) + a*JVS( 528 )
  W( 80 ) = W( 80 ) + a*JVS( 529 )
  a = -W( 63 ) / JVS(          533  )
  W( 63 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 534 )
  W( 73 ) = W( 73 ) + a*JVS( 535 )
  W( 74 ) = W( 74 ) + a*JVS( 536 )
  W( 75 ) = W( 75 ) + a*JVS( 537 )
  W( 76 ) = W( 76 ) + a*JVS( 538 )
  W( 77 ) = W( 77 ) + a*JVS( 539 )
  W( 79 ) = W( 79 ) + a*JVS( 540 )
  a = -W( 64 ) / JVS(          554  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 555 )
  W( 66 ) = W( 66 ) + a*JVS( 556 )
  W( 67 ) = W( 67 ) + a*JVS( 557 )
  W( 68 ) = W( 68 ) + a*JVS( 558 )
  W( 69 ) = W( 69 ) + a*JVS( 559 )
  W( 70 ) = W( 70 ) + a*JVS( 560 )
  W( 71 ) = W( 71 ) + a*JVS( 561 )
  W( 72 ) = W( 72 ) + a*JVS( 562 )
  W( 73 ) = W( 73 ) + a*JVS( 563 )
  W( 74 ) = W( 74 ) + a*JVS( 564 )
  W( 75 ) = W( 75 ) + a*JVS( 565 )
  W( 76 ) = W( 76 ) + a*JVS( 566 )
  W( 77 ) = W( 77 ) + a*JVS( 567 )
  W( 78 ) = W( 78 ) + a*JVS( 568 )
  W( 79 ) = W( 79 ) + a*JVS( 569 )
  W( 80 ) = W( 80 ) + a*JVS( 570 )
  a = -W( 65 ) / JVS(          573  )
  W( 65 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 574 )
  W( 73 ) = W( 73 ) + a*JVS( 575 )
  W( 74 ) = W( 74 ) + a*JVS( 576 )
  W( 75 ) = W( 75 ) + a*JVS( 577 )
  W( 76 ) = W( 76 ) + a*JVS( 578 )
  W( 77 ) = W( 77 ) + a*JVS( 579 )
  W( 79 ) = W( 79 ) + a*JVS( 580 )
  a = -W( 66 ) / JVS(          587  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 588 )
  W( 72 ) = W( 72 ) + a*JVS( 589 )
  W( 73 ) = W( 73 ) + a*JVS( 590 )
  W( 74 ) = W( 74 ) + a*JVS( 591 )
  W( 75 ) = W( 75 ) + a*JVS( 592 )
  W( 76 ) = W( 76 ) + a*JVS( 593 )
  W( 77 ) = W( 77 ) + a*JVS( 594 )
  W( 79 ) = W( 79 ) + a*JVS( 595 )
  a = -W( 67 ) / JVS(          603  )
  W( 67 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 604 )
  W( 73 ) = W( 73 ) + a*JVS( 605 )
  W( 74 ) = W( 74 ) + a*JVS( 606 )
  W( 75 ) = W( 75 ) + a*JVS( 607 )
  W( 76 ) = W( 76 ) + a*JVS( 608 )
  W( 77 ) = W( 77 ) + a*JVS( 609 )
  W( 79 ) = W( 79 ) + a*JVS( 610 )
  a = -W( 68 ) / JVS(          624  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 625 )
  W( 70 ) = W( 70 ) + a*JVS( 626 )
  W( 71 ) = W( 71 ) + a*JVS( 627 )
  W( 72 ) = W( 72 ) + a*JVS( 628 )
  W( 73 ) = W( 73 ) + a*JVS( 629 )
  W( 74 ) = W( 74 ) + a*JVS( 630 )
  W( 75 ) = W( 75 ) + a*JVS( 631 )
  W( 76 ) = W( 76 ) + a*JVS( 632 )
  W( 77 ) = W( 77 ) + a*JVS( 633 )
  W( 78 ) = W( 78 ) + a*JVS( 634 )
  W( 79 ) = W( 79 ) + a*JVS( 635 )
  W( 80 ) = W( 80 ) + a*JVS( 636 )
  a = -W( 69 ) / JVS(          650  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 651 )
  W( 71 ) = W( 71 ) + a*JVS( 652 )
  W( 72 ) = W( 72 ) + a*JVS( 653 )
  W( 73 ) = W( 73 ) + a*JVS( 654 )
  W( 74 ) = W( 74 ) + a*JVS( 655 )
  W( 75 ) = W( 75 ) + a*JVS( 656 )
  W( 76 ) = W( 76 ) + a*JVS( 657 )
  W( 77 ) = W( 77 ) + a*JVS( 658 )
  W( 78 ) = W( 78 ) + a*JVS( 659 )
  W( 79 ) = W( 79 ) + a*JVS( 660 )
  W( 80 ) = W( 80 ) + a*JVS( 661 )
  a = -W( 70 ) / JVS(          677  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 678 )
  W( 72 ) = W( 72 ) + a*JVS( 679 )
  W( 73 ) = W( 73 ) + a*JVS( 680 )
  W( 74 ) = W( 74 ) + a*JVS( 681 )
  W( 75 ) = W( 75 ) + a*JVS( 682 )
  W( 76 ) = W( 76 ) + a*JVS( 683 )
  W( 77 ) = W( 77 ) + a*JVS( 684 )
  W( 78 ) = W( 78 ) + a*JVS( 685 )
  W( 79 ) = W( 79 ) + a*JVS( 686 )
  W( 80 ) = W( 80 ) + a*JVS( 687 )
  a = -W( 71 ) / JVS(          693  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 694 )
  W( 73 ) = W( 73 ) + a*JVS( 695 )
  W( 74 ) = W( 74 ) + a*JVS( 696 )
  W( 75 ) = W( 75 ) + a*JVS( 697 )
  W( 76 ) = W( 76 ) + a*JVS( 698 )
  W( 77 ) = W( 77 ) + a*JVS( 699 )
  W( 78 ) = W( 78 ) + a*JVS( 700 )
  W( 79 ) = W( 79 ) + a*JVS( 701 )
  W( 80 ) = W( 80 ) + a*JVS( 702 )
  a = -W( 72 ) / JVS(          736  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 737 )
  W( 74 ) = W( 74 ) + a*JVS( 738 )
  W( 75 ) = W( 75 ) + a*JVS( 739 )
  W( 76 ) = W( 76 ) + a*JVS( 740 )
  W( 77 ) = W( 77 ) + a*JVS( 741 )
  W( 78 ) = W( 78 ) + a*JVS( 742 )
  W( 79 ) = W( 79 ) + a*JVS( 743 )
  W( 80 ) = W( 80 ) + a*JVS( 744 )
  a = -W( 73 ) / JVS(          806  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 807 )
  W( 75 ) = W( 75 ) + a*JVS( 808 )
  W( 76 ) = W( 76 ) + a*JVS( 809 )
  W( 77 ) = W( 77 ) + a*JVS( 810 )
  W( 78 ) = W( 78 ) + a*JVS( 811 )
  W( 79 ) = W( 79 ) + a*JVS( 812 )
  W( 80 ) = W( 80 ) + a*JVS( 813 )
  JVS( 814) = W( 11 )
  JVS( 815) = W( 19 )
  JVS( 816) = W( 26 )
  JVS( 817) = W( 30 )
  JVS( 818) = W( 32 )
  JVS( 819) = W( 35 )
  JVS( 820) = W( 36 )
  JVS( 821) = W( 38 )
  JVS( 822) = W( 39 )
  JVS( 823) = W( 40 )
  JVS( 824) = W( 41 )
  JVS( 825) = W( 43 )
  JVS( 826) = W( 44 )
  JVS( 827) = W( 45 )
  JVS( 828) = W( 46 )
  JVS( 829) = W( 47 )
  JVS( 830) = W( 51 )
  JVS( 831) = W( 52 )
  JVS( 832) = W( 53 )
  JVS( 833) = W( 54 )
  JVS( 834) = W( 55 )
  JVS( 835) = W( 56 )
  JVS( 836) = W( 57 )
  JVS( 837) = W( 58 )
  JVS( 838) = W( 59 )
  JVS( 839) = W( 61 )
  JVS( 840) = W( 62 )
  JVS( 841) = W( 63 )
  JVS( 842) = W( 64 )
  JVS( 843) = W( 65 )
  JVS( 844) = W( 66 )
  JVS( 845) = W( 67 )
  JVS( 846) = W( 68 )
  JVS( 847) = W( 69 )
  JVS( 848) = W( 70 )
  JVS( 849) = W( 71 )
  JVS( 850) = W( 72 )
  JVS( 851) = W( 73 )
  JVS( 852) = W( 74 )
  JVS( 853) = W( 75 )
  JVS( 854) = W( 76 )
  JVS( 855) = W( 77 )
  JVS( 856) = W( 78 )
  JVS( 857) = W( 79 )
  JVS( 858) = W( 80 )
  IF ( ABS(  JVS( 916 )) < TINY(a) ) THEN
         IER = 75                                      
         RETURN
  END IF
   W( 5 ) = JVS( 859 )
   W( 8 ) = JVS( 860 )
   W( 9 ) = JVS( 861 )
   W( 10 ) = JVS( 862 )
   W( 12 ) = JVS( 863 )
   W( 14 ) = JVS( 864 )
   W( 16 ) = JVS( 865 )
   W( 19 ) = JVS( 866 )
   W( 20 ) = JVS( 867 )
   W( 22 ) = JVS( 868 )
   W( 24 ) = JVS( 869 )
   W( 25 ) = JVS( 870 )
   W( 27 ) = JVS( 871 )
   W( 28 ) = JVS( 872 )
   W( 29 ) = JVS( 873 )
   W( 30 ) = JVS( 874 )
   W( 31 ) = JVS( 875 )
   W( 33 ) = JVS( 876 )
   W( 34 ) = JVS( 877 )
   W( 36 ) = JVS( 878 )
   W( 37 ) = JVS( 879 )
   W( 38 ) = JVS( 880 )
   W( 39 ) = JVS( 881 )
   W( 40 ) = JVS( 882 )
   W( 41 ) = JVS( 883 )
   W( 42 ) = JVS( 884 )
   W( 43 ) = JVS( 885 )
   W( 44 ) = JVS( 886 )
   W( 45 ) = JVS( 887 )
   W( 46 ) = JVS( 888 )
   W( 47 ) = JVS( 889 )
   W( 48 ) = JVS( 890 )
   W( 50 ) = JVS( 891 )
   W( 51 ) = JVS( 892 )
   W( 52 ) = JVS( 893 )
   W( 53 ) = JVS( 894 )
   W( 54 ) = JVS( 895 )
   W( 55 ) = JVS( 896 )
   W( 56 ) = JVS( 897 )
   W( 57 ) = JVS( 898 )
   W( 58 ) = JVS( 899 )
   W( 59 ) = JVS( 900 )
   W( 60 ) = JVS( 901 )
   W( 61 ) = JVS( 902 )
   W( 62 ) = JVS( 903 )
   W( 63 ) = JVS( 904 )
   W( 64 ) = JVS( 905 )
   W( 65 ) = JVS( 906 )
   W( 66 ) = JVS( 907 )
   W( 67 ) = JVS( 908 )
   W( 68 ) = JVS( 909 )
   W( 69 ) = JVS( 910 )
   W( 70 ) = JVS( 911 )
   W( 71 ) = JVS( 912 )
   W( 72 ) = JVS( 913 )
   W( 73 ) = JVS( 914 )
   W( 74 ) = JVS( 915 )
   W( 75 ) = JVS( 916 )
   W( 76 ) = JVS( 917 )
   W( 77 ) = JVS( 918 )
   W( 78 ) = JVS( 919 )
   W( 79 ) = JVS( 920 )
   W( 80 ) = JVS( 921 )
  a = -W( 5 ) / JVS(           40  )
  W( 5 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 41 )
  a = -W( 8 ) / JVS(           48  )
  W( 8 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 49 )
  a = -W( 9 ) / JVS(           50  )
  W( 9 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 51 )
  a = -W( 10 ) / JVS(           52  )
  W( 10 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 53 )
  a = -W( 12 ) / JVS(           57  )
  W( 12 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 58 )
  a = -W( 14 ) / JVS(           63  )
  W( 14 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 64 )
  a = -W( 16 ) / JVS(           70  )
  W( 16 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 71 )
  W( 73 ) = W( 73 ) + a*JVS( 72 )
  a = -W( 19 ) / JVS(           81  )
  W( 19 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 82 )
  W( 75 ) = W( 75 ) + a*JVS( 83 )
  W( 80 ) = W( 80 ) + a*JVS( 84 )
  a = -W( 20 ) / JVS(           85  )
  W( 20 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 86 )
  W( 75 ) = W( 75 ) + a*JVS( 87 )
  W( 76 ) = W( 76 ) + a*JVS( 88 )
  a = -W( 22 ) / JVS(           94  )
  W( 22 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 95 )
  W( 42 ) = W( 42 ) + a*JVS( 96 )
  W( 73 ) = W( 73 ) + a*JVS( 97 )
  W( 79 ) = W( 79 ) + a*JVS( 98 )
  W( 80 ) = W( 80 ) + a*JVS( 99 )
  a = -W( 24 ) / JVS(          107  )
  W( 24 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 108 )
  W( 38 ) = W( 38 ) + a*JVS( 109 )
  W( 39 ) = W( 39 ) + a*JVS( 110 )
  W( 41 ) = W( 41 ) + a*JVS( 111 )
  W( 63 ) = W( 63 ) + a*JVS( 112 )
  W( 65 ) = W( 65 ) + a*JVS( 113 )
  W( 73 ) = W( 73 ) + a*JVS( 114 )
  W( 75 ) = W( 75 ) + a*JVS( 115 )
  W( 77 ) = W( 77 ) + a*JVS( 116 )
  a = -W( 25 ) / JVS(          118  )
  W( 25 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 119 )
  W( 69 ) = W( 69 ) + a*JVS( 120 )
  W( 73 ) = W( 73 ) + a*JVS( 121 )
  W( 76 ) = W( 76 ) + a*JVS( 122 )
  a = -W( 27 ) / JVS(          130  )
  W( 27 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 131 )
  W( 73 ) = W( 73 ) + a*JVS( 132 )
  W( 74 ) = W( 74 ) + a*JVS( 133 )
  W( 75 ) = W( 75 ) + a*JVS( 134 )
  W( 80 ) = W( 80 ) + a*JVS( 135 )
  a = -W( 28 ) / JVS(          137  )
  W( 28 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 138 )
  W( 77 ) = W( 77 ) + a*JVS( 139 )
  W( 80 ) = W( 80 ) + a*JVS( 140 )
  a = -W( 29 ) / JVS(          142  )
  W( 29 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 143 )
  W( 77 ) = W( 77 ) + a*JVS( 144 )
  W( 80 ) = W( 80 ) + a*JVS( 145 )
  a = -W( 30 ) / JVS(          146  )
  W( 30 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 147 )
  W( 74 ) = W( 74 ) + a*JVS( 148 )
  W( 77 ) = W( 77 ) + a*JVS( 149 )
  a = -W( 31 ) / JVS(          150  )
  W( 31 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 151 )
  W( 73 ) = W( 73 ) + a*JVS( 152 )
  W( 77 ) = W( 77 ) + a*JVS( 153 )
  W( 80 ) = W( 80 ) + a*JVS( 154 )
  a = -W( 33 ) / JVS(          166  )
  W( 33 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 167 )
  W( 72 ) = W( 72 ) + a*JVS( 168 )
  W( 73 ) = W( 73 ) + a*JVS( 169 )
  W( 75 ) = W( 75 ) + a*JVS( 170 )
  W( 77 ) = W( 77 ) + a*JVS( 171 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 173 )
  W( 50 ) = W( 50 ) + a*JVS( 174 )
  W( 73 ) = W( 73 ) + a*JVS( 175 )
  W( 74 ) = W( 74 ) + a*JVS( 176 )
  W( 79 ) = W( 79 ) + a*JVS( 177 )
  a = -W( 36 ) / JVS(          185  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 186 )
  W( 74 ) = W( 74 ) + a*JVS( 187 )
  W( 77 ) = W( 77 ) + a*JVS( 188 )
  a = -W( 37 ) / JVS(          194  )
  W( 37 ) = -a
  W( 38 ) = W( 38 ) + a*JVS( 195 )
  W( 39 ) = W( 39 ) + a*JVS( 196 )
  W( 41 ) = W( 41 ) + a*JVS( 197 )
  W( 42 ) = W( 42 ) + a*JVS( 198 )
  W( 43 ) = W( 43 ) + a*JVS( 199 )
  W( 44 ) = W( 44 ) + a*JVS( 200 )
  W( 48 ) = W( 48 ) + a*JVS( 201 )
  W( 50 ) = W( 50 ) + a*JVS( 202 )
  W( 51 ) = W( 51 ) + a*JVS( 203 )
  W( 56 ) = W( 56 ) + a*JVS( 204 )
  W( 61 ) = W( 61 ) + a*JVS( 205 )
  W( 63 ) = W( 63 ) + a*JVS( 206 )
  W( 64 ) = W( 64 ) + a*JVS( 207 )
  W( 65 ) = W( 65 ) + a*JVS( 208 )
  W( 73 ) = W( 73 ) + a*JVS( 209 )
  W( 74 ) = W( 74 ) + a*JVS( 210 )
  W( 77 ) = W( 77 ) + a*JVS( 211 )
  W( 79 ) = W( 79 ) + a*JVS( 212 )
  a = -W( 38 ) / JVS(          213  )
  W( 38 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 214 )
  W( 74 ) = W( 74 ) + a*JVS( 215 )
  W( 77 ) = W( 77 ) + a*JVS( 216 )
  a = -W( 39 ) / JVS(          217  )
  W( 39 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 218 )
  W( 74 ) = W( 74 ) + a*JVS( 219 )
  W( 77 ) = W( 77 ) + a*JVS( 220 )
  a = -W( 40 ) / JVS(          225  )
  W( 40 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 226 )
  W( 74 ) = W( 74 ) + a*JVS( 227 )
  W( 75 ) = W( 75 ) + a*JVS( 228 )
  W( 77 ) = W( 77 ) + a*JVS( 229 )
  W( 80 ) = W( 80 ) + a*JVS( 230 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  a = -W( 42 ) / JVS(          238  )
  W( 42 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 239 )
  W( 73 ) = W( 73 ) + a*JVS( 240 )
  W( 74 ) = W( 74 ) + a*JVS( 241 )
  W( 75 ) = W( 75 ) + a*JVS( 242 )
  W( 77 ) = W( 77 ) + a*JVS( 243 )
  W( 79 ) = W( 79 ) + a*JVS( 244 )
  W( 80 ) = W( 80 ) + a*JVS( 245 )
  a = -W( 43 ) / JVS(          246  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 247 )
  W( 74 ) = W( 74 ) + a*JVS( 248 )
  W( 77 ) = W( 77 ) + a*JVS( 249 )
  W( 78 ) = W( 78 ) + a*JVS( 250 )
  W( 80 ) = W( 80 ) + a*JVS( 251 )
  a = -W( 44 ) / JVS(          253  )
  W( 44 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 254 )
  W( 57 ) = W( 57 ) + a*JVS( 255 )
  W( 58 ) = W( 58 ) + a*JVS( 256 )
  W( 61 ) = W( 61 ) + a*JVS( 257 )
  W( 71 ) = W( 71 ) + a*JVS( 258 )
  W( 72 ) = W( 72 ) + a*JVS( 259 )
  W( 73 ) = W( 73 ) + a*JVS( 260 )
  W( 74 ) = W( 74 ) + a*JVS( 261 )
  W( 76 ) = W( 76 ) + a*JVS( 262 )
  W( 77 ) = W( 77 ) + a*JVS( 263 )
  W( 79 ) = W( 79 ) + a*JVS( 264 )
  a = -W( 45 ) / JVS(          266  )
  W( 45 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 267 )
  W( 73 ) = W( 73 ) + a*JVS( 268 )
  W( 74 ) = W( 74 ) + a*JVS( 269 )
  W( 75 ) = W( 75 ) + a*JVS( 270 )
  W( 76 ) = W( 76 ) + a*JVS( 271 )
  W( 77 ) = W( 77 ) + a*JVS( 272 )
  W( 79 ) = W( 79 ) + a*JVS( 273 )
  a = -W( 46 ) / JVS(          274  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 275 )
  W( 72 ) = W( 72 ) + a*JVS( 276 )
  W( 73 ) = W( 73 ) + a*JVS( 277 )
  W( 74 ) = W( 74 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 76 ) = W( 76 ) + a*JVS( 280 )
  W( 79 ) = W( 79 ) + a*JVS( 281 )
  a = -W( 47 ) / JVS(          282  )
  W( 47 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 283 )
  W( 72 ) = W( 72 ) + a*JVS( 284 )
  W( 73 ) = W( 73 ) + a*JVS( 285 )
  W( 74 ) = W( 74 ) + a*JVS( 286 )
  W( 75 ) = W( 75 ) + a*JVS( 287 )
  W( 76 ) = W( 76 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  a = -W( 48 ) / JVS(          294  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 295 )
  W( 54 ) = W( 54 ) + a*JVS( 296 )
  W( 72 ) = W( 72 ) + a*JVS( 297 )
  W( 73 ) = W( 73 ) + a*JVS( 298 )
  W( 74 ) = W( 74 ) + a*JVS( 299 )
  W( 75 ) = W( 75 ) + a*JVS( 300 )
  W( 76 ) = W( 76 ) + a*JVS( 301 )
  W( 77 ) = W( 77 ) + a*JVS( 302 )
  W( 79 ) = W( 79 ) + a*JVS( 303 )
  a = -W( 50 ) / JVS(          331  )
  W( 50 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 332 )
  W( 73 ) = W( 73 ) + a*JVS( 333 )
  W( 74 ) = W( 74 ) + a*JVS( 334 )
  W( 75 ) = W( 75 ) + a*JVS( 335 )
  W( 76 ) = W( 76 ) + a*JVS( 336 )
  W( 77 ) = W( 77 ) + a*JVS( 337 )
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  a = -W( 51 ) / JVS(          347  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 348 )
  W( 54 ) = W( 54 ) + a*JVS( 349 )
  W( 57 ) = W( 57 ) + a*JVS( 350 )
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 61 ) = W( 61 ) + a*JVS( 352 )
  W( 69 ) = W( 69 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 72 ) = W( 72 ) + a*JVS( 355 )
  W( 73 ) = W( 73 ) + a*JVS( 356 )
  W( 74 ) = W( 74 ) + a*JVS( 357 )
  W( 75 ) = W( 75 ) + a*JVS( 358 )
  W( 76 ) = W( 76 ) + a*JVS( 359 )
  W( 77 ) = W( 77 ) + a*JVS( 360 )
  W( 78 ) = W( 78 ) + a*JVS( 361 )
  W( 79 ) = W( 79 ) + a*JVS( 362 )
  W( 80 ) = W( 80 ) + a*JVS( 363 )
  a = -W( 52 ) / JVS(          366  )
  W( 52 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 367 )
  W( 73 ) = W( 73 ) + a*JVS( 368 )
  W( 74 ) = W( 74 ) + a*JVS( 369 )
  W( 75 ) = W( 75 ) + a*JVS( 370 )
  W( 76 ) = W( 76 ) + a*JVS( 371 )
  W( 77 ) = W( 77 ) + a*JVS( 372 )
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 80 ) = W( 80 ) + a*JVS( 374 )
  a = -W( 53 ) / JVS(          376  )
  W( 53 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 377 )
  W( 73 ) = W( 73 ) + a*JVS( 378 )
  W( 74 ) = W( 74 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 79 ) = W( 79 ) + a*JVS( 382 )
  a = -W( 54 ) / JVS(          384  )
  W( 54 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 385 )
  W( 73 ) = W( 73 ) + a*JVS( 386 )
  W( 74 ) = W( 74 ) + a*JVS( 387 )
  W( 75 ) = W( 75 ) + a*JVS( 388 )
  W( 76 ) = W( 76 ) + a*JVS( 389 )
  W( 77 ) = W( 77 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 55 ) / JVS(          393  )
  W( 55 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  W( 73 ) = W( 73 ) + a*JVS( 395 )
  W( 74 ) = W( 74 ) + a*JVS( 396 )
  W( 75 ) = W( 75 ) + a*JVS( 397 )
  W( 76 ) = W( 76 ) + a*JVS( 398 )
  W( 79 ) = W( 79 ) + a*JVS( 399 )
  a = -W( 56 ) / JVS(          423  )
  W( 56 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 424 )
  W( 58 ) = W( 58 ) + a*JVS( 425 )
  W( 59 ) = W( 59 ) + a*JVS( 426 )
  W( 61 ) = W( 61 ) + a*JVS( 427 )
  W( 62 ) = W( 62 ) + a*JVS( 428 )
  W( 63 ) = W( 63 ) + a*JVS( 429 )
  W( 65 ) = W( 65 ) + a*JVS( 430 )
  W( 66 ) = W( 66 ) + a*JVS( 431 )
  W( 67 ) = W( 67 ) + a*JVS( 432 )
  W( 68 ) = W( 68 ) + a*JVS( 433 )
  W( 69 ) = W( 69 ) + a*JVS( 434 )
  W( 71 ) = W( 71 ) + a*JVS( 435 )
  W( 72 ) = W( 72 ) + a*JVS( 436 )
  W( 73 ) = W( 73 ) + a*JVS( 437 )
  W( 74 ) = W( 74 ) + a*JVS( 438 )
  W( 75 ) = W( 75 ) + a*JVS( 439 )
  W( 76 ) = W( 76 ) + a*JVS( 440 )
  W( 77 ) = W( 77 ) + a*JVS( 441 )
  W( 78 ) = W( 78 ) + a*JVS( 442 )
  W( 79 ) = W( 79 ) + a*JVS( 443 )
  W( 80 ) = W( 80 ) + a*JVS( 444 )
  a = -W( 57 ) / JVS(          446  )
  W( 57 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 447 )
  W( 73 ) = W( 73 ) + a*JVS( 448 )
  W( 74 ) = W( 74 ) + a*JVS( 449 )
  W( 75 ) = W( 75 ) + a*JVS( 450 )
  W( 76 ) = W( 76 ) + a*JVS( 451 )
  W( 77 ) = W( 77 ) + a*JVS( 452 )
  W( 79 ) = W( 79 ) + a*JVS( 453 )
  W( 80 ) = W( 80 ) + a*JVS( 454 )
  a = -W( 58 ) / JVS(          456  )
  W( 58 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 457 )
  W( 73 ) = W( 73 ) + a*JVS( 458 )
  W( 74 ) = W( 74 ) + a*JVS( 459 )
  W( 75 ) = W( 75 ) + a*JVS( 460 )
  W( 76 ) = W( 76 ) + a*JVS( 461 )
  W( 77 ) = W( 77 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  a = -W( 59 ) / JVS(          466  )
  W( 59 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 467 )
  W( 73 ) = W( 73 ) + a*JVS( 468 )
  W( 74 ) = W( 74 ) + a*JVS( 469 )
  W( 75 ) = W( 75 ) + a*JVS( 470 )
  W( 76 ) = W( 76 ) + a*JVS( 471 )
  W( 77 ) = W( 77 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  a = -W( 60 ) / JVS(          483  )
  W( 60 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 484 )
  W( 67 ) = W( 67 ) + a*JVS( 485 )
  W( 71 ) = W( 71 ) + a*JVS( 486 )
  W( 72 ) = W( 72 ) + a*JVS( 487 )
  W( 73 ) = W( 73 ) + a*JVS( 488 )
  W( 74 ) = W( 74 ) + a*JVS( 489 )
  W( 75 ) = W( 75 ) + a*JVS( 490 )
  W( 76 ) = W( 76 ) + a*JVS( 491 )
  W( 77 ) = W( 77 ) + a*JVS( 492 )
  W( 78 ) = W( 78 ) + a*JVS( 493 )
  W( 79 ) = W( 79 ) + a*JVS( 494 )
  W( 80 ) = W( 80 ) + a*JVS( 495 )
  a = -W( 61 ) / JVS(          498  )
  W( 61 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 499 )
  W( 73 ) = W( 73 ) + a*JVS( 500 )
  W( 74 ) = W( 74 ) + a*JVS( 501 )
  W( 75 ) = W( 75 ) + a*JVS( 502 )
  W( 76 ) = W( 76 ) + a*JVS( 503 )
  W( 77 ) = W( 77 ) + a*JVS( 504 )
  W( 79 ) = W( 79 ) + a*JVS( 505 )
  W( 80 ) = W( 80 ) + a*JVS( 506 )
  a = -W( 62 ) / JVS(          517  )
  W( 62 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 518 )
  W( 70 ) = W( 70 ) + a*JVS( 519 )
  W( 71 ) = W( 71 ) + a*JVS( 520 )
  W( 72 ) = W( 72 ) + a*JVS( 521 )
  W( 73 ) = W( 73 ) + a*JVS( 522 )
  W( 74 ) = W( 74 ) + a*JVS( 523 )
  W( 75 ) = W( 75 ) + a*JVS( 524 )
  W( 76 ) = W( 76 ) + a*JVS( 525 )
  W( 77 ) = W( 77 ) + a*JVS( 526 )
  W( 78 ) = W( 78 ) + a*JVS( 527 )
  W( 79 ) = W( 79 ) + a*JVS( 528 )
  W( 80 ) = W( 80 ) + a*JVS( 529 )
  a = -W( 63 ) / JVS(          533  )
  W( 63 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 534 )
  W( 73 ) = W( 73 ) + a*JVS( 535 )
  W( 74 ) = W( 74 ) + a*JVS( 536 )
  W( 75 ) = W( 75 ) + a*JVS( 537 )
  W( 76 ) = W( 76 ) + a*JVS( 538 )
  W( 77 ) = W( 77 ) + a*JVS( 539 )
  W( 79 ) = W( 79 ) + a*JVS( 540 )
  a = -W( 64 ) / JVS(          554  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 555 )
  W( 66 ) = W( 66 ) + a*JVS( 556 )
  W( 67 ) = W( 67 ) + a*JVS( 557 )
  W( 68 ) = W( 68 ) + a*JVS( 558 )
  W( 69 ) = W( 69 ) + a*JVS( 559 )
  W( 70 ) = W( 70 ) + a*JVS( 560 )
  W( 71 ) = W( 71 ) + a*JVS( 561 )
  W( 72 ) = W( 72 ) + a*JVS( 562 )
  W( 73 ) = W( 73 ) + a*JVS( 563 )
  W( 74 ) = W( 74 ) + a*JVS( 564 )
  W( 75 ) = W( 75 ) + a*JVS( 565 )
  W( 76 ) = W( 76 ) + a*JVS( 566 )
  W( 77 ) = W( 77 ) + a*JVS( 567 )
  W( 78 ) = W( 78 ) + a*JVS( 568 )
  W( 79 ) = W( 79 ) + a*JVS( 569 )
  W( 80 ) = W( 80 ) + a*JVS( 570 )
  a = -W( 65 ) / JVS(          573  )
  W( 65 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 574 )
  W( 73 ) = W( 73 ) + a*JVS( 575 )
  W( 74 ) = W( 74 ) + a*JVS( 576 )
  W( 75 ) = W( 75 ) + a*JVS( 577 )
  W( 76 ) = W( 76 ) + a*JVS( 578 )
  W( 77 ) = W( 77 ) + a*JVS( 579 )
  W( 79 ) = W( 79 ) + a*JVS( 580 )
  a = -W( 66 ) / JVS(          587  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 588 )
  W( 72 ) = W( 72 ) + a*JVS( 589 )
  W( 73 ) = W( 73 ) + a*JVS( 590 )
  W( 74 ) = W( 74 ) + a*JVS( 591 )
  W( 75 ) = W( 75 ) + a*JVS( 592 )
  W( 76 ) = W( 76 ) + a*JVS( 593 )
  W( 77 ) = W( 77 ) + a*JVS( 594 )
  W( 79 ) = W( 79 ) + a*JVS( 595 )
  a = -W( 67 ) / JVS(          603  )
  W( 67 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 604 )
  W( 73 ) = W( 73 ) + a*JVS( 605 )
  W( 74 ) = W( 74 ) + a*JVS( 606 )
  W( 75 ) = W( 75 ) + a*JVS( 607 )
  W( 76 ) = W( 76 ) + a*JVS( 608 )
  W( 77 ) = W( 77 ) + a*JVS( 609 )
  W( 79 ) = W( 79 ) + a*JVS( 610 )
  a = -W( 68 ) / JVS(          624  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 625 )
  W( 70 ) = W( 70 ) + a*JVS( 626 )
  W( 71 ) = W( 71 ) + a*JVS( 627 )
  W( 72 ) = W( 72 ) + a*JVS( 628 )
  W( 73 ) = W( 73 ) + a*JVS( 629 )
  W( 74 ) = W( 74 ) + a*JVS( 630 )
  W( 75 ) = W( 75 ) + a*JVS( 631 )
  W( 76 ) = W( 76 ) + a*JVS( 632 )
  W( 77 ) = W( 77 ) + a*JVS( 633 )
  W( 78 ) = W( 78 ) + a*JVS( 634 )
  W( 79 ) = W( 79 ) + a*JVS( 635 )
  W( 80 ) = W( 80 ) + a*JVS( 636 )
  a = -W( 69 ) / JVS(          650  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 651 )
  W( 71 ) = W( 71 ) + a*JVS( 652 )
  W( 72 ) = W( 72 ) + a*JVS( 653 )
  W( 73 ) = W( 73 ) + a*JVS( 654 )
  W( 74 ) = W( 74 ) + a*JVS( 655 )
  W( 75 ) = W( 75 ) + a*JVS( 656 )
  W( 76 ) = W( 76 ) + a*JVS( 657 )
  W( 77 ) = W( 77 ) + a*JVS( 658 )
  W( 78 ) = W( 78 ) + a*JVS( 659 )
  W( 79 ) = W( 79 ) + a*JVS( 660 )
  W( 80 ) = W( 80 ) + a*JVS( 661 )
  a = -W( 70 ) / JVS(          677  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 678 )
  W( 72 ) = W( 72 ) + a*JVS( 679 )
  W( 73 ) = W( 73 ) + a*JVS( 680 )
  W( 74 ) = W( 74 ) + a*JVS( 681 )
  W( 75 ) = W( 75 ) + a*JVS( 682 )
  W( 76 ) = W( 76 ) + a*JVS( 683 )
  W( 77 ) = W( 77 ) + a*JVS( 684 )
  W( 78 ) = W( 78 ) + a*JVS( 685 )
  W( 79 ) = W( 79 ) + a*JVS( 686 )
  W( 80 ) = W( 80 ) + a*JVS( 687 )
  a = -W( 71 ) / JVS(          693  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 694 )
  W( 73 ) = W( 73 ) + a*JVS( 695 )
  W( 74 ) = W( 74 ) + a*JVS( 696 )
  W( 75 ) = W( 75 ) + a*JVS( 697 )
  W( 76 ) = W( 76 ) + a*JVS( 698 )
  W( 77 ) = W( 77 ) + a*JVS( 699 )
  W( 78 ) = W( 78 ) + a*JVS( 700 )
  W( 79 ) = W( 79 ) + a*JVS( 701 )
  W( 80 ) = W( 80 ) + a*JVS( 702 )
  a = -W( 72 ) / JVS(          736  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 737 )
  W( 74 ) = W( 74 ) + a*JVS( 738 )
  W( 75 ) = W( 75 ) + a*JVS( 739 )
  W( 76 ) = W( 76 ) + a*JVS( 740 )
  W( 77 ) = W( 77 ) + a*JVS( 741 )
  W( 78 ) = W( 78 ) + a*JVS( 742 )
  W( 79 ) = W( 79 ) + a*JVS( 743 )
  W( 80 ) = W( 80 ) + a*JVS( 744 )
  a = -W( 73 ) / JVS(          806  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 807 )
  W( 75 ) = W( 75 ) + a*JVS( 808 )
  W( 76 ) = W( 76 ) + a*JVS( 809 )
  W( 77 ) = W( 77 ) + a*JVS( 810 )
  W( 78 ) = W( 78 ) + a*JVS( 811 )
  W( 79 ) = W( 79 ) + a*JVS( 812 )
  W( 80 ) = W( 80 ) + a*JVS( 813 )
  a = -W( 74 ) / JVS(          852  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 853 )
  W( 76 ) = W( 76 ) + a*JVS( 854 )
  W( 77 ) = W( 77 ) + a*JVS( 855 )
  W( 78 ) = W( 78 ) + a*JVS( 856 )
  W( 79 ) = W( 79 ) + a*JVS( 857 )
  W( 80 ) = W( 80 ) + a*JVS( 858 )
  JVS( 859) = W( 5 )
  JVS( 860) = W( 8 )
  JVS( 861) = W( 9 )
  JVS( 862) = W( 10 )
  JVS( 863) = W( 12 )
  JVS( 864) = W( 14 )
  JVS( 865) = W( 16 )
  JVS( 866) = W( 19 )
  JVS( 867) = W( 20 )
  JVS( 868) = W( 22 )
  JVS( 869) = W( 24 )
  JVS( 870) = W( 25 )
  JVS( 871) = W( 27 )
  JVS( 872) = W( 28 )
  JVS( 873) = W( 29 )
  JVS( 874) = W( 30 )
  JVS( 875) = W( 31 )
  JVS( 876) = W( 33 )
  JVS( 877) = W( 34 )
  JVS( 878) = W( 36 )
  JVS( 879) = W( 37 )
  JVS( 880) = W( 38 )
  JVS( 881) = W( 39 )
  JVS( 882) = W( 40 )
  JVS( 883) = W( 41 )
  JVS( 884) = W( 42 )
  JVS( 885) = W( 43 )
  JVS( 886) = W( 44 )
  JVS( 887) = W( 45 )
  JVS( 888) = W( 46 )
  JVS( 889) = W( 47 )
  JVS( 890) = W( 48 )
  JVS( 891) = W( 50 )
  JVS( 892) = W( 51 )
  JVS( 893) = W( 52 )
  JVS( 894) = W( 53 )
  JVS( 895) = W( 54 )
  JVS( 896) = W( 55 )
  JVS( 897) = W( 56 )
  JVS( 898) = W( 57 )
  JVS( 899) = W( 58 )
  JVS( 900) = W( 59 )
  JVS( 901) = W( 60 )
  JVS( 902) = W( 61 )
  JVS( 903) = W( 62 )
  JVS( 904) = W( 63 )
  JVS( 905) = W( 64 )
  JVS( 906) = W( 65 )
  JVS( 907) = W( 66 )
  JVS( 908) = W( 67 )
  JVS( 909) = W( 68 )
  JVS( 910) = W( 69 )
  JVS( 911) = W( 70 )
  JVS( 912) = W( 71 )
  JVS( 913) = W( 72 )
  JVS( 914) = W( 73 )
  JVS( 915) = W( 74 )
  JVS( 916) = W( 75 )
  JVS( 917) = W( 76 )
  JVS( 918) = W( 77 )
  JVS( 919) = W( 78 )
  JVS( 920) = W( 79 )
  JVS( 921) = W( 80 )
  IF ( ABS(  JVS( 953 )) < TINY(a) ) THEN
         IER = 76                                      
         RETURN
  END IF
   W( 15 ) = JVS( 922 )
   W( 20 ) = JVS( 923 )
   W( 33 ) = JVS( 924 )
   W( 39 ) = JVS( 925 )
   W( 41 ) = JVS( 926 )
   W( 45 ) = JVS( 927 )
   W( 46 ) = JVS( 928 )
   W( 47 ) = JVS( 929 )
   W( 50 ) = JVS( 930 )
   W( 52 ) = JVS( 931 )
   W( 53 ) = JVS( 932 )
   W( 54 ) = JVS( 933 )
   W( 55 ) = JVS( 934 )
   W( 57 ) = JVS( 935 )
   W( 58 ) = JVS( 936 )
   W( 59 ) = JVS( 937 )
   W( 61 ) = JVS( 938 )
   W( 62 ) = JVS( 939 )
   W( 63 ) = JVS( 940 )
   W( 64 ) = JVS( 941 )
   W( 65 ) = JVS( 942 )
   W( 66 ) = JVS( 943 )
   W( 67 ) = JVS( 944 )
   W( 68 ) = JVS( 945 )
   W( 69 ) = JVS( 946 )
   W( 70 ) = JVS( 947 )
   W( 71 ) = JVS( 948 )
   W( 72 ) = JVS( 949 )
   W( 73 ) = JVS( 950 )
   W( 74 ) = JVS( 951 )
   W( 75 ) = JVS( 952 )
   W( 76 ) = JVS( 953 )
   W( 77 ) = JVS( 954 )
   W( 78 ) = JVS( 955 )
   W( 79 ) = JVS( 956 )
   W( 80 ) = JVS( 957 )
  a = -W( 15 ) / JVS(           65  )
  W( 15 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 66 )
  W( 65 ) = W( 65 ) + a*JVS( 67 )
  W( 73 ) = W( 73 ) + a*JVS( 68 )
  W( 77 ) = W( 77 ) + a*JVS( 69 )
  a = -W( 20 ) / JVS(           85  )
  W( 20 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 86 )
  W( 75 ) = W( 75 ) + a*JVS( 87 )
  W( 76 ) = W( 76 ) + a*JVS( 88 )
  a = -W( 33 ) / JVS(          166  )
  W( 33 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 167 )
  W( 72 ) = W( 72 ) + a*JVS( 168 )
  W( 73 ) = W( 73 ) + a*JVS( 169 )
  W( 75 ) = W( 75 ) + a*JVS( 170 )
  W( 77 ) = W( 77 ) + a*JVS( 171 )
  a = -W( 39 ) / JVS(          217  )
  W( 39 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 218 )
  W( 74 ) = W( 74 ) + a*JVS( 219 )
  W( 77 ) = W( 77 ) + a*JVS( 220 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  a = -W( 45 ) / JVS(          266  )
  W( 45 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 267 )
  W( 73 ) = W( 73 ) + a*JVS( 268 )
  W( 74 ) = W( 74 ) + a*JVS( 269 )
  W( 75 ) = W( 75 ) + a*JVS( 270 )
  W( 76 ) = W( 76 ) + a*JVS( 271 )
  W( 77 ) = W( 77 ) + a*JVS( 272 )
  W( 79 ) = W( 79 ) + a*JVS( 273 )
  a = -W( 46 ) / JVS(          274  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 275 )
  W( 72 ) = W( 72 ) + a*JVS( 276 )
  W( 73 ) = W( 73 ) + a*JVS( 277 )
  W( 74 ) = W( 74 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 76 ) = W( 76 ) + a*JVS( 280 )
  W( 79 ) = W( 79 ) + a*JVS( 281 )
  a = -W( 47 ) / JVS(          282  )
  W( 47 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 283 )
  W( 72 ) = W( 72 ) + a*JVS( 284 )
  W( 73 ) = W( 73 ) + a*JVS( 285 )
  W( 74 ) = W( 74 ) + a*JVS( 286 )
  W( 75 ) = W( 75 ) + a*JVS( 287 )
  W( 76 ) = W( 76 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  a = -W( 50 ) / JVS(          331  )
  W( 50 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 332 )
  W( 73 ) = W( 73 ) + a*JVS( 333 )
  W( 74 ) = W( 74 ) + a*JVS( 334 )
  W( 75 ) = W( 75 ) + a*JVS( 335 )
  W( 76 ) = W( 76 ) + a*JVS( 336 )
  W( 77 ) = W( 77 ) + a*JVS( 337 )
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  a = -W( 52 ) / JVS(          366  )
  W( 52 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 367 )
  W( 73 ) = W( 73 ) + a*JVS( 368 )
  W( 74 ) = W( 74 ) + a*JVS( 369 )
  W( 75 ) = W( 75 ) + a*JVS( 370 )
  W( 76 ) = W( 76 ) + a*JVS( 371 )
  W( 77 ) = W( 77 ) + a*JVS( 372 )
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 80 ) = W( 80 ) + a*JVS( 374 )
  a = -W( 53 ) / JVS(          376  )
  W( 53 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 377 )
  W( 73 ) = W( 73 ) + a*JVS( 378 )
  W( 74 ) = W( 74 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 79 ) = W( 79 ) + a*JVS( 382 )
  a = -W( 54 ) / JVS(          384  )
  W( 54 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 385 )
  W( 73 ) = W( 73 ) + a*JVS( 386 )
  W( 74 ) = W( 74 ) + a*JVS( 387 )
  W( 75 ) = W( 75 ) + a*JVS( 388 )
  W( 76 ) = W( 76 ) + a*JVS( 389 )
  W( 77 ) = W( 77 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 55 ) / JVS(          393  )
  W( 55 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  W( 73 ) = W( 73 ) + a*JVS( 395 )
  W( 74 ) = W( 74 ) + a*JVS( 396 )
  W( 75 ) = W( 75 ) + a*JVS( 397 )
  W( 76 ) = W( 76 ) + a*JVS( 398 )
  W( 79 ) = W( 79 ) + a*JVS( 399 )
  a = -W( 57 ) / JVS(          446  )
  W( 57 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 447 )
  W( 73 ) = W( 73 ) + a*JVS( 448 )
  W( 74 ) = W( 74 ) + a*JVS( 449 )
  W( 75 ) = W( 75 ) + a*JVS( 450 )
  W( 76 ) = W( 76 ) + a*JVS( 451 )
  W( 77 ) = W( 77 ) + a*JVS( 452 )
  W( 79 ) = W( 79 ) + a*JVS( 453 )
  W( 80 ) = W( 80 ) + a*JVS( 454 )
  a = -W( 58 ) / JVS(          456  )
  W( 58 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 457 )
  W( 73 ) = W( 73 ) + a*JVS( 458 )
  W( 74 ) = W( 74 ) + a*JVS( 459 )
  W( 75 ) = W( 75 ) + a*JVS( 460 )
  W( 76 ) = W( 76 ) + a*JVS( 461 )
  W( 77 ) = W( 77 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  a = -W( 59 ) / JVS(          466  )
  W( 59 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 467 )
  W( 73 ) = W( 73 ) + a*JVS( 468 )
  W( 74 ) = W( 74 ) + a*JVS( 469 )
  W( 75 ) = W( 75 ) + a*JVS( 470 )
  W( 76 ) = W( 76 ) + a*JVS( 471 )
  W( 77 ) = W( 77 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  a = -W( 61 ) / JVS(          498  )
  W( 61 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 499 )
  W( 73 ) = W( 73 ) + a*JVS( 500 )
  W( 74 ) = W( 74 ) + a*JVS( 501 )
  W( 75 ) = W( 75 ) + a*JVS( 502 )
  W( 76 ) = W( 76 ) + a*JVS( 503 )
  W( 77 ) = W( 77 ) + a*JVS( 504 )
  W( 79 ) = W( 79 ) + a*JVS( 505 )
  W( 80 ) = W( 80 ) + a*JVS( 506 )
  a = -W( 62 ) / JVS(          517  )
  W( 62 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 518 )
  W( 70 ) = W( 70 ) + a*JVS( 519 )
  W( 71 ) = W( 71 ) + a*JVS( 520 )
  W( 72 ) = W( 72 ) + a*JVS( 521 )
  W( 73 ) = W( 73 ) + a*JVS( 522 )
  W( 74 ) = W( 74 ) + a*JVS( 523 )
  W( 75 ) = W( 75 ) + a*JVS( 524 )
  W( 76 ) = W( 76 ) + a*JVS( 525 )
  W( 77 ) = W( 77 ) + a*JVS( 526 )
  W( 78 ) = W( 78 ) + a*JVS( 527 )
  W( 79 ) = W( 79 ) + a*JVS( 528 )
  W( 80 ) = W( 80 ) + a*JVS( 529 )
  a = -W( 63 ) / JVS(          533  )
  W( 63 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 534 )
  W( 73 ) = W( 73 ) + a*JVS( 535 )
  W( 74 ) = W( 74 ) + a*JVS( 536 )
  W( 75 ) = W( 75 ) + a*JVS( 537 )
  W( 76 ) = W( 76 ) + a*JVS( 538 )
  W( 77 ) = W( 77 ) + a*JVS( 539 )
  W( 79 ) = W( 79 ) + a*JVS( 540 )
  a = -W( 64 ) / JVS(          554  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 555 )
  W( 66 ) = W( 66 ) + a*JVS( 556 )
  W( 67 ) = W( 67 ) + a*JVS( 557 )
  W( 68 ) = W( 68 ) + a*JVS( 558 )
  W( 69 ) = W( 69 ) + a*JVS( 559 )
  W( 70 ) = W( 70 ) + a*JVS( 560 )
  W( 71 ) = W( 71 ) + a*JVS( 561 )
  W( 72 ) = W( 72 ) + a*JVS( 562 )
  W( 73 ) = W( 73 ) + a*JVS( 563 )
  W( 74 ) = W( 74 ) + a*JVS( 564 )
  W( 75 ) = W( 75 ) + a*JVS( 565 )
  W( 76 ) = W( 76 ) + a*JVS( 566 )
  W( 77 ) = W( 77 ) + a*JVS( 567 )
  W( 78 ) = W( 78 ) + a*JVS( 568 )
  W( 79 ) = W( 79 ) + a*JVS( 569 )
  W( 80 ) = W( 80 ) + a*JVS( 570 )
  a = -W( 65 ) / JVS(          573  )
  W( 65 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 574 )
  W( 73 ) = W( 73 ) + a*JVS( 575 )
  W( 74 ) = W( 74 ) + a*JVS( 576 )
  W( 75 ) = W( 75 ) + a*JVS( 577 )
  W( 76 ) = W( 76 ) + a*JVS( 578 )
  W( 77 ) = W( 77 ) + a*JVS( 579 )
  W( 79 ) = W( 79 ) + a*JVS( 580 )
  a = -W( 66 ) / JVS(          587  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 588 )
  W( 72 ) = W( 72 ) + a*JVS( 589 )
  W( 73 ) = W( 73 ) + a*JVS( 590 )
  W( 74 ) = W( 74 ) + a*JVS( 591 )
  W( 75 ) = W( 75 ) + a*JVS( 592 )
  W( 76 ) = W( 76 ) + a*JVS( 593 )
  W( 77 ) = W( 77 ) + a*JVS( 594 )
  W( 79 ) = W( 79 ) + a*JVS( 595 )
  a = -W( 67 ) / JVS(          603  )
  W( 67 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 604 )
  W( 73 ) = W( 73 ) + a*JVS( 605 )
  W( 74 ) = W( 74 ) + a*JVS( 606 )
  W( 75 ) = W( 75 ) + a*JVS( 607 )
  W( 76 ) = W( 76 ) + a*JVS( 608 )
  W( 77 ) = W( 77 ) + a*JVS( 609 )
  W( 79 ) = W( 79 ) + a*JVS( 610 )
  a = -W( 68 ) / JVS(          624  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 625 )
  W( 70 ) = W( 70 ) + a*JVS( 626 )
  W( 71 ) = W( 71 ) + a*JVS( 627 )
  W( 72 ) = W( 72 ) + a*JVS( 628 )
  W( 73 ) = W( 73 ) + a*JVS( 629 )
  W( 74 ) = W( 74 ) + a*JVS( 630 )
  W( 75 ) = W( 75 ) + a*JVS( 631 )
  W( 76 ) = W( 76 ) + a*JVS( 632 )
  W( 77 ) = W( 77 ) + a*JVS( 633 )
  W( 78 ) = W( 78 ) + a*JVS( 634 )
  W( 79 ) = W( 79 ) + a*JVS( 635 )
  W( 80 ) = W( 80 ) + a*JVS( 636 )
  a = -W( 69 ) / JVS(          650  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 651 )
  W( 71 ) = W( 71 ) + a*JVS( 652 )
  W( 72 ) = W( 72 ) + a*JVS( 653 )
  W( 73 ) = W( 73 ) + a*JVS( 654 )
  W( 74 ) = W( 74 ) + a*JVS( 655 )
  W( 75 ) = W( 75 ) + a*JVS( 656 )
  W( 76 ) = W( 76 ) + a*JVS( 657 )
  W( 77 ) = W( 77 ) + a*JVS( 658 )
  W( 78 ) = W( 78 ) + a*JVS( 659 )
  W( 79 ) = W( 79 ) + a*JVS( 660 )
  W( 80 ) = W( 80 ) + a*JVS( 661 )
  a = -W( 70 ) / JVS(          677  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 678 )
  W( 72 ) = W( 72 ) + a*JVS( 679 )
  W( 73 ) = W( 73 ) + a*JVS( 680 )
  W( 74 ) = W( 74 ) + a*JVS( 681 )
  W( 75 ) = W( 75 ) + a*JVS( 682 )
  W( 76 ) = W( 76 ) + a*JVS( 683 )
  W( 77 ) = W( 77 ) + a*JVS( 684 )
  W( 78 ) = W( 78 ) + a*JVS( 685 )
  W( 79 ) = W( 79 ) + a*JVS( 686 )
  W( 80 ) = W( 80 ) + a*JVS( 687 )
  a = -W( 71 ) / JVS(          693  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 694 )
  W( 73 ) = W( 73 ) + a*JVS( 695 )
  W( 74 ) = W( 74 ) + a*JVS( 696 )
  W( 75 ) = W( 75 ) + a*JVS( 697 )
  W( 76 ) = W( 76 ) + a*JVS( 698 )
  W( 77 ) = W( 77 ) + a*JVS( 699 )
  W( 78 ) = W( 78 ) + a*JVS( 700 )
  W( 79 ) = W( 79 ) + a*JVS( 701 )
  W( 80 ) = W( 80 ) + a*JVS( 702 )
  a = -W( 72 ) / JVS(          736  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 737 )
  W( 74 ) = W( 74 ) + a*JVS( 738 )
  W( 75 ) = W( 75 ) + a*JVS( 739 )
  W( 76 ) = W( 76 ) + a*JVS( 740 )
  W( 77 ) = W( 77 ) + a*JVS( 741 )
  W( 78 ) = W( 78 ) + a*JVS( 742 )
  W( 79 ) = W( 79 ) + a*JVS( 743 )
  W( 80 ) = W( 80 ) + a*JVS( 744 )
  a = -W( 73 ) / JVS(          806  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 807 )
  W( 75 ) = W( 75 ) + a*JVS( 808 )
  W( 76 ) = W( 76 ) + a*JVS( 809 )
  W( 77 ) = W( 77 ) + a*JVS( 810 )
  W( 78 ) = W( 78 ) + a*JVS( 811 )
  W( 79 ) = W( 79 ) + a*JVS( 812 )
  W( 80 ) = W( 80 ) + a*JVS( 813 )
  a = -W( 74 ) / JVS(          852  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 853 )
  W( 76 ) = W( 76 ) + a*JVS( 854 )
  W( 77 ) = W( 77 ) + a*JVS( 855 )
  W( 78 ) = W( 78 ) + a*JVS( 856 )
  W( 79 ) = W( 79 ) + a*JVS( 857 )
  W( 80 ) = W( 80 ) + a*JVS( 858 )
  a = -W( 75 ) / JVS(          916  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 917 )
  W( 77 ) = W( 77 ) + a*JVS( 918 )
  W( 78 ) = W( 78 ) + a*JVS( 919 )
  W( 79 ) = W( 79 ) + a*JVS( 920 )
  W( 80 ) = W( 80 ) + a*JVS( 921 )
  JVS( 922) = W( 15 )
  JVS( 923) = W( 20 )
  JVS( 924) = W( 33 )
  JVS( 925) = W( 39 )
  JVS( 926) = W( 41 )
  JVS( 927) = W( 45 )
  JVS( 928) = W( 46 )
  JVS( 929) = W( 47 )
  JVS( 930) = W( 50 )
  JVS( 931) = W( 52 )
  JVS( 932) = W( 53 )
  JVS( 933) = W( 54 )
  JVS( 934) = W( 55 )
  JVS( 935) = W( 57 )
  JVS( 936) = W( 58 )
  JVS( 937) = W( 59 )
  JVS( 938) = W( 61 )
  JVS( 939) = W( 62 )
  JVS( 940) = W( 63 )
  JVS( 941) = W( 64 )
  JVS( 942) = W( 65 )
  JVS( 943) = W( 66 )
  JVS( 944) = W( 67 )
  JVS( 945) = W( 68 )
  JVS( 946) = W( 69 )
  JVS( 947) = W( 70 )
  JVS( 948) = W( 71 )
  JVS( 949) = W( 72 )
  JVS( 950) = W( 73 )
  JVS( 951) = W( 74 )
  JVS( 952) = W( 75 )
  JVS( 953) = W( 76 )
  JVS( 954) = W( 77 )
  JVS( 955) = W( 78 )
  JVS( 956) = W( 79 )
  JVS( 957) = W( 80 )
  IF ( ABS(  JVS( 980 )) < TINY(a) ) THEN
         IER = 77                                      
         RETURN
  END IF
   W( 26 ) = JVS( 958 )
   W( 28 ) = JVS( 959 )
   W( 29 ) = JVS( 960 )
   W( 30 ) = JVS( 961 )
   W( 31 ) = JVS( 962 )
   W( 36 ) = JVS( 963 )
   W( 38 ) = JVS( 964 )
   W( 39 ) = JVS( 965 )
   W( 40 ) = JVS( 966 )
   W( 41 ) = JVS( 967 )
   W( 43 ) = JVS( 968 )
   W( 48 ) = JVS( 969 )
   W( 50 ) = JVS( 970 )
   W( 54 ) = JVS( 971 )
   W( 61 ) = JVS( 972 )
   W( 63 ) = JVS( 973 )
   W( 65 ) = JVS( 974 )
   W( 72 ) = JVS( 975 )
   W( 73 ) = JVS( 976 )
   W( 74 ) = JVS( 977 )
   W( 75 ) = JVS( 978 )
   W( 76 ) = JVS( 979 )
   W( 77 ) = JVS( 980 )
   W( 78 ) = JVS( 981 )
   W( 79 ) = JVS( 982 )
   W( 80 ) = JVS( 983 )
  a = -W( 26 ) / JVS(          124  )
  W( 26 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 125 )
  W( 74 ) = W( 74 ) + a*JVS( 126 )
  W( 77 ) = W( 77 ) + a*JVS( 127 )
  W( 79 ) = W( 79 ) + a*JVS( 128 )
  W( 80 ) = W( 80 ) + a*JVS( 129 )
  a = -W( 28 ) / JVS(          137  )
  W( 28 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 138 )
  W( 77 ) = W( 77 ) + a*JVS( 139 )
  W( 80 ) = W( 80 ) + a*JVS( 140 )
  a = -W( 29 ) / JVS(          142  )
  W( 29 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 143 )
  W( 77 ) = W( 77 ) + a*JVS( 144 )
  W( 80 ) = W( 80 ) + a*JVS( 145 )
  a = -W( 30 ) / JVS(          146  )
  W( 30 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 147 )
  W( 74 ) = W( 74 ) + a*JVS( 148 )
  W( 77 ) = W( 77 ) + a*JVS( 149 )
  a = -W( 31 ) / JVS(          150  )
  W( 31 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 151 )
  W( 73 ) = W( 73 ) + a*JVS( 152 )
  W( 77 ) = W( 77 ) + a*JVS( 153 )
  W( 80 ) = W( 80 ) + a*JVS( 154 )
  a = -W( 36 ) / JVS(          185  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 186 )
  W( 74 ) = W( 74 ) + a*JVS( 187 )
  W( 77 ) = W( 77 ) + a*JVS( 188 )
  a = -W( 38 ) / JVS(          213  )
  W( 38 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 214 )
  W( 74 ) = W( 74 ) + a*JVS( 215 )
  W( 77 ) = W( 77 ) + a*JVS( 216 )
  a = -W( 39 ) / JVS(          217  )
  W( 39 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 218 )
  W( 74 ) = W( 74 ) + a*JVS( 219 )
  W( 77 ) = W( 77 ) + a*JVS( 220 )
  a = -W( 40 ) / JVS(          225  )
  W( 40 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 226 )
  W( 74 ) = W( 74 ) + a*JVS( 227 )
  W( 75 ) = W( 75 ) + a*JVS( 228 )
  W( 77 ) = W( 77 ) + a*JVS( 229 )
  W( 80 ) = W( 80 ) + a*JVS( 230 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  a = -W( 43 ) / JVS(          246  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 247 )
  W( 74 ) = W( 74 ) + a*JVS( 248 )
  W( 77 ) = W( 77 ) + a*JVS( 249 )
  W( 78 ) = W( 78 ) + a*JVS( 250 )
  W( 80 ) = W( 80 ) + a*JVS( 251 )
  a = -W( 48 ) / JVS(          294  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 295 )
  W( 54 ) = W( 54 ) + a*JVS( 296 )
  W( 72 ) = W( 72 ) + a*JVS( 297 )
  W( 73 ) = W( 73 ) + a*JVS( 298 )
  W( 74 ) = W( 74 ) + a*JVS( 299 )
  W( 75 ) = W( 75 ) + a*JVS( 300 )
  W( 76 ) = W( 76 ) + a*JVS( 301 )
  W( 77 ) = W( 77 ) + a*JVS( 302 )
  W( 79 ) = W( 79 ) + a*JVS( 303 )
  a = -W( 50 ) / JVS(          331  )
  W( 50 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 332 )
  W( 73 ) = W( 73 ) + a*JVS( 333 )
  W( 74 ) = W( 74 ) + a*JVS( 334 )
  W( 75 ) = W( 75 ) + a*JVS( 335 )
  W( 76 ) = W( 76 ) + a*JVS( 336 )
  W( 77 ) = W( 77 ) + a*JVS( 337 )
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  a = -W( 54 ) / JVS(          384  )
  W( 54 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 385 )
  W( 73 ) = W( 73 ) + a*JVS( 386 )
  W( 74 ) = W( 74 ) + a*JVS( 387 )
  W( 75 ) = W( 75 ) + a*JVS( 388 )
  W( 76 ) = W( 76 ) + a*JVS( 389 )
  W( 77 ) = W( 77 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 61 ) / JVS(          498  )
  W( 61 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 499 )
  W( 73 ) = W( 73 ) + a*JVS( 500 )
  W( 74 ) = W( 74 ) + a*JVS( 501 )
  W( 75 ) = W( 75 ) + a*JVS( 502 )
  W( 76 ) = W( 76 ) + a*JVS( 503 )
  W( 77 ) = W( 77 ) + a*JVS( 504 )
  W( 79 ) = W( 79 ) + a*JVS( 505 )
  W( 80 ) = W( 80 ) + a*JVS( 506 )
  a = -W( 63 ) / JVS(          533  )
  W( 63 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 534 )
  W( 73 ) = W( 73 ) + a*JVS( 535 )
  W( 74 ) = W( 74 ) + a*JVS( 536 )
  W( 75 ) = W( 75 ) + a*JVS( 537 )
  W( 76 ) = W( 76 ) + a*JVS( 538 )
  W( 77 ) = W( 77 ) + a*JVS( 539 )
  W( 79 ) = W( 79 ) + a*JVS( 540 )
  a = -W( 65 ) / JVS(          573  )
  W( 65 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 574 )
  W( 73 ) = W( 73 ) + a*JVS( 575 )
  W( 74 ) = W( 74 ) + a*JVS( 576 )
  W( 75 ) = W( 75 ) + a*JVS( 577 )
  W( 76 ) = W( 76 ) + a*JVS( 578 )
  W( 77 ) = W( 77 ) + a*JVS( 579 )
  W( 79 ) = W( 79 ) + a*JVS( 580 )
  a = -W( 72 ) / JVS(          736  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 737 )
  W( 74 ) = W( 74 ) + a*JVS( 738 )
  W( 75 ) = W( 75 ) + a*JVS( 739 )
  W( 76 ) = W( 76 ) + a*JVS( 740 )
  W( 77 ) = W( 77 ) + a*JVS( 741 )
  W( 78 ) = W( 78 ) + a*JVS( 742 )
  W( 79 ) = W( 79 ) + a*JVS( 743 )
  W( 80 ) = W( 80 ) + a*JVS( 744 )
  a = -W( 73 ) / JVS(          806  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 807 )
  W( 75 ) = W( 75 ) + a*JVS( 808 )
  W( 76 ) = W( 76 ) + a*JVS( 809 )
  W( 77 ) = W( 77 ) + a*JVS( 810 )
  W( 78 ) = W( 78 ) + a*JVS( 811 )
  W( 79 ) = W( 79 ) + a*JVS( 812 )
  W( 80 ) = W( 80 ) + a*JVS( 813 )
  a = -W( 74 ) / JVS(          852  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 853 )
  W( 76 ) = W( 76 ) + a*JVS( 854 )
  W( 77 ) = W( 77 ) + a*JVS( 855 )
  W( 78 ) = W( 78 ) + a*JVS( 856 )
  W( 79 ) = W( 79 ) + a*JVS( 857 )
  W( 80 ) = W( 80 ) + a*JVS( 858 )
  a = -W( 75 ) / JVS(          916  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 917 )
  W( 77 ) = W( 77 ) + a*JVS( 918 )
  W( 78 ) = W( 78 ) + a*JVS( 919 )
  W( 79 ) = W( 79 ) + a*JVS( 920 )
  W( 80 ) = W( 80 ) + a*JVS( 921 )
  a = -W( 76 ) / JVS(          953  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 954 )
  W( 78 ) = W( 78 ) + a*JVS( 955 )
  W( 79 ) = W( 79 ) + a*JVS( 956 )
  W( 80 ) = W( 80 ) + a*JVS( 957 )
  JVS( 958) = W( 26 )
  JVS( 959) = W( 28 )
  JVS( 960) = W( 29 )
  JVS( 961) = W( 30 )
  JVS( 962) = W( 31 )
  JVS( 963) = W( 36 )
  JVS( 964) = W( 38 )
  JVS( 965) = W( 39 )
  JVS( 966) = W( 40 )
  JVS( 967) = W( 41 )
  JVS( 968) = W( 43 )
  JVS( 969) = W( 48 )
  JVS( 970) = W( 50 )
  JVS( 971) = W( 54 )
  JVS( 972) = W( 61 )
  JVS( 973) = W( 63 )
  JVS( 974) = W( 65 )
  JVS( 975) = W( 72 )
  JVS( 976) = W( 73 )
  JVS( 977) = W( 74 )
  JVS( 978) = W( 75 )
  JVS( 979) = W( 76 )
  JVS( 980) = W( 77 )
  JVS( 981) = W( 78 )
  JVS( 982) = W( 79 )
  JVS( 983) = W( 80 )
  IF ( ABS(  JVS( 992 )) < TINY(a) ) THEN
         IER = 78                                      
         RETURN
  END IF
   W( 43 ) = JVS( 984 )
   W( 61 ) = JVS( 985 )
   W( 72 ) = JVS( 986 )
   W( 73 ) = JVS( 987 )
   W( 74 ) = JVS( 988 )
   W( 75 ) = JVS( 989 )
   W( 76 ) = JVS( 990 )
   W( 77 ) = JVS( 991 )
   W( 78 ) = JVS( 992 )
   W( 79 ) = JVS( 993 )
   W( 80 ) = JVS( 994 )
  a = -W( 43 ) / JVS(          246  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 247 )
  W( 74 ) = W( 74 ) + a*JVS( 248 )
  W( 77 ) = W( 77 ) + a*JVS( 249 )
  W( 78 ) = W( 78 ) + a*JVS( 250 )
  W( 80 ) = W( 80 ) + a*JVS( 251 )
  a = -W( 61 ) / JVS(          498  )
  W( 61 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 499 )
  W( 73 ) = W( 73 ) + a*JVS( 500 )
  W( 74 ) = W( 74 ) + a*JVS( 501 )
  W( 75 ) = W( 75 ) + a*JVS( 502 )
  W( 76 ) = W( 76 ) + a*JVS( 503 )
  W( 77 ) = W( 77 ) + a*JVS( 504 )
  W( 79 ) = W( 79 ) + a*JVS( 505 )
  W( 80 ) = W( 80 ) + a*JVS( 506 )
  a = -W( 72 ) / JVS(          736  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 737 )
  W( 74 ) = W( 74 ) + a*JVS( 738 )
  W( 75 ) = W( 75 ) + a*JVS( 739 )
  W( 76 ) = W( 76 ) + a*JVS( 740 )
  W( 77 ) = W( 77 ) + a*JVS( 741 )
  W( 78 ) = W( 78 ) + a*JVS( 742 )
  W( 79 ) = W( 79 ) + a*JVS( 743 )
  W( 80 ) = W( 80 ) + a*JVS( 744 )
  a = -W( 73 ) / JVS(          806  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 807 )
  W( 75 ) = W( 75 ) + a*JVS( 808 )
  W( 76 ) = W( 76 ) + a*JVS( 809 )
  W( 77 ) = W( 77 ) + a*JVS( 810 )
  W( 78 ) = W( 78 ) + a*JVS( 811 )
  W( 79 ) = W( 79 ) + a*JVS( 812 )
  W( 80 ) = W( 80 ) + a*JVS( 813 )
  a = -W( 74 ) / JVS(          852  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 853 )
  W( 76 ) = W( 76 ) + a*JVS( 854 )
  W( 77 ) = W( 77 ) + a*JVS( 855 )
  W( 78 ) = W( 78 ) + a*JVS( 856 )
  W( 79 ) = W( 79 ) + a*JVS( 857 )
  W( 80 ) = W( 80 ) + a*JVS( 858 )
  a = -W( 75 ) / JVS(          916  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 917 )
  W( 77 ) = W( 77 ) + a*JVS( 918 )
  W( 78 ) = W( 78 ) + a*JVS( 919 )
  W( 79 ) = W( 79 ) + a*JVS( 920 )
  W( 80 ) = W( 80 ) + a*JVS( 921 )
  a = -W( 76 ) / JVS(          953  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 954 )
  W( 78 ) = W( 78 ) + a*JVS( 955 )
  W( 79 ) = W( 79 ) + a*JVS( 956 )
  W( 80 ) = W( 80 ) + a*JVS( 957 )
  a = -W( 77 ) / JVS(          980  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 981 )
  W( 79 ) = W( 79 ) + a*JVS( 982 )
  W( 80 ) = W( 80 ) + a*JVS( 983 )
  JVS( 984) = W( 43 )
  JVS( 985) = W( 61 )
  JVS( 986) = W( 72 )
  JVS( 987) = W( 73 )
  JVS( 988) = W( 74 )
  JVS( 989) = W( 75 )
  JVS( 990) = W( 76 )
  JVS( 991) = W( 77 )
  JVS( 992) = W( 78 )
  JVS( 993) = W( 79 )
  JVS( 994) = W( 80 )
  IF ( ABS(  JVS( 1031 )) < TINY(a) ) THEN
         IER = 79                                      
         RETURN
  END IF
   W( 23 ) = JVS( 995 )
   W( 26 ) = JVS( 996 )
   W( 28 ) = JVS( 997 )
   W( 29 ) = JVS( 998 )
   W( 31 ) = JVS( 999 )
   W( 40 ) = JVS( 1000 )
   W( 41 ) = JVS( 1001 )
   W( 42 ) = JVS( 1002 )
   W( 45 ) = JVS( 1003 )
   W( 46 ) = JVS( 1004 )
   W( 47 ) = JVS( 1005 )
   W( 48 ) = JVS( 1006 )
   W( 50 ) = JVS( 1007 )
   W( 52 ) = JVS( 1008 )
   W( 53 ) = JVS( 1009 )
   W( 54 ) = JVS( 1010 )
   W( 55 ) = JVS( 1011 )
   W( 57 ) = JVS( 1012 )
   W( 58 ) = JVS( 1013 )
   W( 59 ) = JVS( 1014 )
   W( 62 ) = JVS( 1015 )
   W( 63 ) = JVS( 1016 )
   W( 65 ) = JVS( 1017 )
   W( 66 ) = JVS( 1018 )
   W( 67 ) = JVS( 1019 )
   W( 68 ) = JVS( 1020 )
   W( 69 ) = JVS( 1021 )
   W( 70 ) = JVS( 1022 )
   W( 71 ) = JVS( 1023 )
   W( 72 ) = JVS( 1024 )
   W( 73 ) = JVS( 1025 )
   W( 74 ) = JVS( 1026 )
   W( 75 ) = JVS( 1027 )
   W( 76 ) = JVS( 1028 )
   W( 77 ) = JVS( 1029 )
   W( 78 ) = JVS( 1030 )
   W( 79 ) = JVS( 1031 )
   W( 80 ) = JVS( 1032 )
  a = -W( 23 ) / JVS(          100  )
  W( 23 ) = -a
  W( 28 ) = W( 28 ) + a*JVS( 101 )
  W( 29 ) = W( 29 ) + a*JVS( 102 )
  W( 31 ) = W( 31 ) + a*JVS( 103 )
  W( 73 ) = W( 73 ) + a*JVS( 104 )
  W( 79 ) = W( 79 ) + a*JVS( 105 )
  W( 80 ) = W( 80 ) + a*JVS( 106 )
  a = -W( 26 ) / JVS(          124  )
  W( 26 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 125 )
  W( 74 ) = W( 74 ) + a*JVS( 126 )
  W( 77 ) = W( 77 ) + a*JVS( 127 )
  W( 79 ) = W( 79 ) + a*JVS( 128 )
  W( 80 ) = W( 80 ) + a*JVS( 129 )
  a = -W( 28 ) / JVS(          137  )
  W( 28 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 138 )
  W( 77 ) = W( 77 ) + a*JVS( 139 )
  W( 80 ) = W( 80 ) + a*JVS( 140 )
  a = -W( 29 ) / JVS(          142  )
  W( 29 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 143 )
  W( 77 ) = W( 77 ) + a*JVS( 144 )
  W( 80 ) = W( 80 ) + a*JVS( 145 )
  a = -W( 31 ) / JVS(          150  )
  W( 31 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 151 )
  W( 73 ) = W( 73 ) + a*JVS( 152 )
  W( 77 ) = W( 77 ) + a*JVS( 153 )
  W( 80 ) = W( 80 ) + a*JVS( 154 )
  a = -W( 40 ) / JVS(          225  )
  W( 40 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 226 )
  W( 74 ) = W( 74 ) + a*JVS( 227 )
  W( 75 ) = W( 75 ) + a*JVS( 228 )
  W( 77 ) = W( 77 ) + a*JVS( 229 )
  W( 80 ) = W( 80 ) + a*JVS( 230 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  a = -W( 42 ) / JVS(          238  )
  W( 42 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 239 )
  W( 73 ) = W( 73 ) + a*JVS( 240 )
  W( 74 ) = W( 74 ) + a*JVS( 241 )
  W( 75 ) = W( 75 ) + a*JVS( 242 )
  W( 77 ) = W( 77 ) + a*JVS( 243 )
  W( 79 ) = W( 79 ) + a*JVS( 244 )
  W( 80 ) = W( 80 ) + a*JVS( 245 )
  a = -W( 45 ) / JVS(          266  )
  W( 45 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 267 )
  W( 73 ) = W( 73 ) + a*JVS( 268 )
  W( 74 ) = W( 74 ) + a*JVS( 269 )
  W( 75 ) = W( 75 ) + a*JVS( 270 )
  W( 76 ) = W( 76 ) + a*JVS( 271 )
  W( 77 ) = W( 77 ) + a*JVS( 272 )
  W( 79 ) = W( 79 ) + a*JVS( 273 )
  a = -W( 46 ) / JVS(          274  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 275 )
  W( 72 ) = W( 72 ) + a*JVS( 276 )
  W( 73 ) = W( 73 ) + a*JVS( 277 )
  W( 74 ) = W( 74 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 76 ) = W( 76 ) + a*JVS( 280 )
  W( 79 ) = W( 79 ) + a*JVS( 281 )
  a = -W( 47 ) / JVS(          282  )
  W( 47 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 283 )
  W( 72 ) = W( 72 ) + a*JVS( 284 )
  W( 73 ) = W( 73 ) + a*JVS( 285 )
  W( 74 ) = W( 74 ) + a*JVS( 286 )
  W( 75 ) = W( 75 ) + a*JVS( 287 )
  W( 76 ) = W( 76 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  a = -W( 48 ) / JVS(          294  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 295 )
  W( 54 ) = W( 54 ) + a*JVS( 296 )
  W( 72 ) = W( 72 ) + a*JVS( 297 )
  W( 73 ) = W( 73 ) + a*JVS( 298 )
  W( 74 ) = W( 74 ) + a*JVS( 299 )
  W( 75 ) = W( 75 ) + a*JVS( 300 )
  W( 76 ) = W( 76 ) + a*JVS( 301 )
  W( 77 ) = W( 77 ) + a*JVS( 302 )
  W( 79 ) = W( 79 ) + a*JVS( 303 )
  a = -W( 50 ) / JVS(          331  )
  W( 50 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 332 )
  W( 73 ) = W( 73 ) + a*JVS( 333 )
  W( 74 ) = W( 74 ) + a*JVS( 334 )
  W( 75 ) = W( 75 ) + a*JVS( 335 )
  W( 76 ) = W( 76 ) + a*JVS( 336 )
  W( 77 ) = W( 77 ) + a*JVS( 337 )
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  a = -W( 52 ) / JVS(          366  )
  W( 52 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 367 )
  W( 73 ) = W( 73 ) + a*JVS( 368 )
  W( 74 ) = W( 74 ) + a*JVS( 369 )
  W( 75 ) = W( 75 ) + a*JVS( 370 )
  W( 76 ) = W( 76 ) + a*JVS( 371 )
  W( 77 ) = W( 77 ) + a*JVS( 372 )
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 80 ) = W( 80 ) + a*JVS( 374 )
  a = -W( 53 ) / JVS(          376  )
  W( 53 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 377 )
  W( 73 ) = W( 73 ) + a*JVS( 378 )
  W( 74 ) = W( 74 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 79 ) = W( 79 ) + a*JVS( 382 )
  a = -W( 54 ) / JVS(          384  )
  W( 54 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 385 )
  W( 73 ) = W( 73 ) + a*JVS( 386 )
  W( 74 ) = W( 74 ) + a*JVS( 387 )
  W( 75 ) = W( 75 ) + a*JVS( 388 )
  W( 76 ) = W( 76 ) + a*JVS( 389 )
  W( 77 ) = W( 77 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 55 ) / JVS(          393  )
  W( 55 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  W( 73 ) = W( 73 ) + a*JVS( 395 )
  W( 74 ) = W( 74 ) + a*JVS( 396 )
  W( 75 ) = W( 75 ) + a*JVS( 397 )
  W( 76 ) = W( 76 ) + a*JVS( 398 )
  W( 79 ) = W( 79 ) + a*JVS( 399 )
  a = -W( 57 ) / JVS(          446  )
  W( 57 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 447 )
  W( 73 ) = W( 73 ) + a*JVS( 448 )
  W( 74 ) = W( 74 ) + a*JVS( 449 )
  W( 75 ) = W( 75 ) + a*JVS( 450 )
  W( 76 ) = W( 76 ) + a*JVS( 451 )
  W( 77 ) = W( 77 ) + a*JVS( 452 )
  W( 79 ) = W( 79 ) + a*JVS( 453 )
  W( 80 ) = W( 80 ) + a*JVS( 454 )
  a = -W( 58 ) / JVS(          456  )
  W( 58 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 457 )
  W( 73 ) = W( 73 ) + a*JVS( 458 )
  W( 74 ) = W( 74 ) + a*JVS( 459 )
  W( 75 ) = W( 75 ) + a*JVS( 460 )
  W( 76 ) = W( 76 ) + a*JVS( 461 )
  W( 77 ) = W( 77 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  a = -W( 59 ) / JVS(          466  )
  W( 59 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 467 )
  W( 73 ) = W( 73 ) + a*JVS( 468 )
  W( 74 ) = W( 74 ) + a*JVS( 469 )
  W( 75 ) = W( 75 ) + a*JVS( 470 )
  W( 76 ) = W( 76 ) + a*JVS( 471 )
  W( 77 ) = W( 77 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  a = -W( 62 ) / JVS(          517  )
  W( 62 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 518 )
  W( 70 ) = W( 70 ) + a*JVS( 519 )
  W( 71 ) = W( 71 ) + a*JVS( 520 )
  W( 72 ) = W( 72 ) + a*JVS( 521 )
  W( 73 ) = W( 73 ) + a*JVS( 522 )
  W( 74 ) = W( 74 ) + a*JVS( 523 )
  W( 75 ) = W( 75 ) + a*JVS( 524 )
  W( 76 ) = W( 76 ) + a*JVS( 525 )
  W( 77 ) = W( 77 ) + a*JVS( 526 )
  W( 78 ) = W( 78 ) + a*JVS( 527 )
  W( 79 ) = W( 79 ) + a*JVS( 528 )
  W( 80 ) = W( 80 ) + a*JVS( 529 )
  a = -W( 63 ) / JVS(          533  )
  W( 63 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 534 )
  W( 73 ) = W( 73 ) + a*JVS( 535 )
  W( 74 ) = W( 74 ) + a*JVS( 536 )
  W( 75 ) = W( 75 ) + a*JVS( 537 )
  W( 76 ) = W( 76 ) + a*JVS( 538 )
  W( 77 ) = W( 77 ) + a*JVS( 539 )
  W( 79 ) = W( 79 ) + a*JVS( 540 )
  a = -W( 65 ) / JVS(          573  )
  W( 65 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 574 )
  W( 73 ) = W( 73 ) + a*JVS( 575 )
  W( 74 ) = W( 74 ) + a*JVS( 576 )
  W( 75 ) = W( 75 ) + a*JVS( 577 )
  W( 76 ) = W( 76 ) + a*JVS( 578 )
  W( 77 ) = W( 77 ) + a*JVS( 579 )
  W( 79 ) = W( 79 ) + a*JVS( 580 )
  a = -W( 66 ) / JVS(          587  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 588 )
  W( 72 ) = W( 72 ) + a*JVS( 589 )
  W( 73 ) = W( 73 ) + a*JVS( 590 )
  W( 74 ) = W( 74 ) + a*JVS( 591 )
  W( 75 ) = W( 75 ) + a*JVS( 592 )
  W( 76 ) = W( 76 ) + a*JVS( 593 )
  W( 77 ) = W( 77 ) + a*JVS( 594 )
  W( 79 ) = W( 79 ) + a*JVS( 595 )
  a = -W( 67 ) / JVS(          603  )
  W( 67 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 604 )
  W( 73 ) = W( 73 ) + a*JVS( 605 )
  W( 74 ) = W( 74 ) + a*JVS( 606 )
  W( 75 ) = W( 75 ) + a*JVS( 607 )
  W( 76 ) = W( 76 ) + a*JVS( 608 )
  W( 77 ) = W( 77 ) + a*JVS( 609 )
  W( 79 ) = W( 79 ) + a*JVS( 610 )
  a = -W( 68 ) / JVS(          624  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 625 )
  W( 70 ) = W( 70 ) + a*JVS( 626 )
  W( 71 ) = W( 71 ) + a*JVS( 627 )
  W( 72 ) = W( 72 ) + a*JVS( 628 )
  W( 73 ) = W( 73 ) + a*JVS( 629 )
  W( 74 ) = W( 74 ) + a*JVS( 630 )
  W( 75 ) = W( 75 ) + a*JVS( 631 )
  W( 76 ) = W( 76 ) + a*JVS( 632 )
  W( 77 ) = W( 77 ) + a*JVS( 633 )
  W( 78 ) = W( 78 ) + a*JVS( 634 )
  W( 79 ) = W( 79 ) + a*JVS( 635 )
  W( 80 ) = W( 80 ) + a*JVS( 636 )
  a = -W( 69 ) / JVS(          650  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 651 )
  W( 71 ) = W( 71 ) + a*JVS( 652 )
  W( 72 ) = W( 72 ) + a*JVS( 653 )
  W( 73 ) = W( 73 ) + a*JVS( 654 )
  W( 74 ) = W( 74 ) + a*JVS( 655 )
  W( 75 ) = W( 75 ) + a*JVS( 656 )
  W( 76 ) = W( 76 ) + a*JVS( 657 )
  W( 77 ) = W( 77 ) + a*JVS( 658 )
  W( 78 ) = W( 78 ) + a*JVS( 659 )
  W( 79 ) = W( 79 ) + a*JVS( 660 )
  W( 80 ) = W( 80 ) + a*JVS( 661 )
  a = -W( 70 ) / JVS(          677  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 678 )
  W( 72 ) = W( 72 ) + a*JVS( 679 )
  W( 73 ) = W( 73 ) + a*JVS( 680 )
  W( 74 ) = W( 74 ) + a*JVS( 681 )
  W( 75 ) = W( 75 ) + a*JVS( 682 )
  W( 76 ) = W( 76 ) + a*JVS( 683 )
  W( 77 ) = W( 77 ) + a*JVS( 684 )
  W( 78 ) = W( 78 ) + a*JVS( 685 )
  W( 79 ) = W( 79 ) + a*JVS( 686 )
  W( 80 ) = W( 80 ) + a*JVS( 687 )
  a = -W( 71 ) / JVS(          693  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 694 )
  W( 73 ) = W( 73 ) + a*JVS( 695 )
  W( 74 ) = W( 74 ) + a*JVS( 696 )
  W( 75 ) = W( 75 ) + a*JVS( 697 )
  W( 76 ) = W( 76 ) + a*JVS( 698 )
  W( 77 ) = W( 77 ) + a*JVS( 699 )
  W( 78 ) = W( 78 ) + a*JVS( 700 )
  W( 79 ) = W( 79 ) + a*JVS( 701 )
  W( 80 ) = W( 80 ) + a*JVS( 702 )
  a = -W( 72 ) / JVS(          736  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 737 )
  W( 74 ) = W( 74 ) + a*JVS( 738 )
  W( 75 ) = W( 75 ) + a*JVS( 739 )
  W( 76 ) = W( 76 ) + a*JVS( 740 )
  W( 77 ) = W( 77 ) + a*JVS( 741 )
  W( 78 ) = W( 78 ) + a*JVS( 742 )
  W( 79 ) = W( 79 ) + a*JVS( 743 )
  W( 80 ) = W( 80 ) + a*JVS( 744 )
  a = -W( 73 ) / JVS(          806  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 807 )
  W( 75 ) = W( 75 ) + a*JVS( 808 )
  W( 76 ) = W( 76 ) + a*JVS( 809 )
  W( 77 ) = W( 77 ) + a*JVS( 810 )
  W( 78 ) = W( 78 ) + a*JVS( 811 )
  W( 79 ) = W( 79 ) + a*JVS( 812 )
  W( 80 ) = W( 80 ) + a*JVS( 813 )
  a = -W( 74 ) / JVS(          852  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 853 )
  W( 76 ) = W( 76 ) + a*JVS( 854 )
  W( 77 ) = W( 77 ) + a*JVS( 855 )
  W( 78 ) = W( 78 ) + a*JVS( 856 )
  W( 79 ) = W( 79 ) + a*JVS( 857 )
  W( 80 ) = W( 80 ) + a*JVS( 858 )
  a = -W( 75 ) / JVS(          916  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 917 )
  W( 77 ) = W( 77 ) + a*JVS( 918 )
  W( 78 ) = W( 78 ) + a*JVS( 919 )
  W( 79 ) = W( 79 ) + a*JVS( 920 )
  W( 80 ) = W( 80 ) + a*JVS( 921 )
  a = -W( 76 ) / JVS(          953  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 954 )
  W( 78 ) = W( 78 ) + a*JVS( 955 )
  W( 79 ) = W( 79 ) + a*JVS( 956 )
  W( 80 ) = W( 80 ) + a*JVS( 957 )
  a = -W( 77 ) / JVS(          980  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 981 )
  W( 79 ) = W( 79 ) + a*JVS( 982 )
  W( 80 ) = W( 80 ) + a*JVS( 983 )
  a = -W( 78 ) / JVS(          992  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 993 )
  W( 80 ) = W( 80 ) + a*JVS( 994 )
  JVS( 995) = W( 23 )
  JVS( 996) = W( 26 )
  JVS( 997) = W( 28 )
  JVS( 998) = W( 29 )
  JVS( 999) = W( 31 )
  JVS( 1000) = W( 40 )
  JVS( 1001) = W( 41 )
  JVS( 1002) = W( 42 )
  JVS( 1003) = W( 45 )
  JVS( 1004) = W( 46 )
  JVS( 1005) = W( 47 )
  JVS( 1006) = W( 48 )
  JVS( 1007) = W( 50 )
  JVS( 1008) = W( 52 )
  JVS( 1009) = W( 53 )
  JVS( 1010) = W( 54 )
  JVS( 1011) = W( 55 )
  JVS( 1012) = W( 57 )
  JVS( 1013) = W( 58 )
  JVS( 1014) = W( 59 )
  JVS( 1015) = W( 62 )
  JVS( 1016) = W( 63 )
  JVS( 1017) = W( 65 )
  JVS( 1018) = W( 66 )
  JVS( 1019) = W( 67 )
  JVS( 1020) = W( 68 )
  JVS( 1021) = W( 69 )
  JVS( 1022) = W( 70 )
  JVS( 1023) = W( 71 )
  JVS( 1024) = W( 72 )
  JVS( 1025) = W( 73 )
  JVS( 1026) = W( 74 )
  JVS( 1027) = W( 75 )
  JVS( 1028) = W( 76 )
  JVS( 1029) = W( 77 )
  JVS( 1030) = W( 78 )
  JVS( 1031) = W( 79 )
  JVS( 1032) = W( 80 )
  IF ( ABS(  JVS( 1086 )) < TINY(a) ) THEN
         IER = 80                                      
         RETURN
  END IF
   W( 11 ) = JVS( 1033 )
   W( 18 ) = JVS( 1034 )
   W( 19 ) = JVS( 1035 )
   W( 21 ) = JVS( 1036 )
   W( 23 ) = JVS( 1037 )
   W( 26 ) = JVS( 1038 )
   W( 27 ) = JVS( 1039 )
   W( 28 ) = JVS( 1040 )
   W( 29 ) = JVS( 1041 )
   W( 31 ) = JVS( 1042 )
   W( 32 ) = JVS( 1043 )
   W( 34 ) = JVS( 1044 )
   W( 35 ) = JVS( 1045 )
   W( 39 ) = JVS( 1046 )
   W( 40 ) = JVS( 1047 )
   W( 41 ) = JVS( 1048 )
   W( 42 ) = JVS( 1049 )
   W( 43 ) = JVS( 1050 )
   W( 44 ) = JVS( 1051 )
   W( 45 ) = JVS( 1052 )
   W( 46 ) = JVS( 1053 )
   W( 47 ) = JVS( 1054 )
   W( 48 ) = JVS( 1055 )
   W( 50 ) = JVS( 1056 )
   W( 51 ) = JVS( 1057 )
   W( 52 ) = JVS( 1058 )
   W( 53 ) = JVS( 1059 )
   W( 54 ) = JVS( 1060 )
   W( 55 ) = JVS( 1061 )
   W( 56 ) = JVS( 1062 )
   W( 57 ) = JVS( 1063 )
   W( 58 ) = JVS( 1064 )
   W( 59 ) = JVS( 1065 )
   W( 60 ) = JVS( 1066 )
   W( 61 ) = JVS( 1067 )
   W( 62 ) = JVS( 1068 )
   W( 63 ) = JVS( 1069 )
   W( 64 ) = JVS( 1070 )
   W( 65 ) = JVS( 1071 )
   W( 66 ) = JVS( 1072 )
   W( 67 ) = JVS( 1073 )
   W( 68 ) = JVS( 1074 )
   W( 69 ) = JVS( 1075 )
   W( 70 ) = JVS( 1076 )
   W( 71 ) = JVS( 1077 )
   W( 72 ) = JVS( 1078 )
   W( 73 ) = JVS( 1079 )
   W( 74 ) = JVS( 1080 )
   W( 75 ) = JVS( 1081 )
   W( 76 ) = JVS( 1082 )
   W( 77 ) = JVS( 1083 )
   W( 78 ) = JVS( 1084 )
   W( 79 ) = JVS( 1085 )
   W( 80 ) = JVS( 1086 )
  a = -W( 11 ) / JVS(           54  )
  W( 11 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 55 )
  W( 80 ) = W( 80 ) + a*JVS( 56 )
  a = -W( 18 ) / JVS(           78  )
  W( 18 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 79 )
  W( 73 ) = W( 73 ) + a*JVS( 80 )
  a = -W( 19 ) / JVS(           81  )
  W( 19 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 82 )
  W( 75 ) = W( 75 ) + a*JVS( 83 )
  W( 80 ) = W( 80 ) + a*JVS( 84 )
  a = -W( 21 ) / JVS(           89  )
  W( 21 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 90 )
  W( 73 ) = W( 73 ) + a*JVS( 91 )
  W( 80 ) = W( 80 ) + a*JVS( 92 )
  a = -W( 23 ) / JVS(          100  )
  W( 23 ) = -a
  W( 28 ) = W( 28 ) + a*JVS( 101 )
  W( 29 ) = W( 29 ) + a*JVS( 102 )
  W( 31 ) = W( 31 ) + a*JVS( 103 )
  W( 73 ) = W( 73 ) + a*JVS( 104 )
  W( 79 ) = W( 79 ) + a*JVS( 105 )
  W( 80 ) = W( 80 ) + a*JVS( 106 )
  a = -W( 26 ) / JVS(          124  )
  W( 26 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 125 )
  W( 74 ) = W( 74 ) + a*JVS( 126 )
  W( 77 ) = W( 77 ) + a*JVS( 127 )
  W( 79 ) = W( 79 ) + a*JVS( 128 )
  W( 80 ) = W( 80 ) + a*JVS( 129 )
  a = -W( 27 ) / JVS(          130  )
  W( 27 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 131 )
  W( 73 ) = W( 73 ) + a*JVS( 132 )
  W( 74 ) = W( 74 ) + a*JVS( 133 )
  W( 75 ) = W( 75 ) + a*JVS( 134 )
  W( 80 ) = W( 80 ) + a*JVS( 135 )
  a = -W( 28 ) / JVS(          137  )
  W( 28 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 138 )
  W( 77 ) = W( 77 ) + a*JVS( 139 )
  W( 80 ) = W( 80 ) + a*JVS( 140 )
  a = -W( 29 ) / JVS(          142  )
  W( 29 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 143 )
  W( 77 ) = W( 77 ) + a*JVS( 144 )
  W( 80 ) = W( 80 ) + a*JVS( 145 )
  a = -W( 31 ) / JVS(          150  )
  W( 31 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 151 )
  W( 73 ) = W( 73 ) + a*JVS( 152 )
  W( 77 ) = W( 77 ) + a*JVS( 153 )
  W( 80 ) = W( 80 ) + a*JVS( 154 )
  a = -W( 32 ) / JVS(          155  )
  W( 32 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 156 )
  W( 44 ) = W( 44 ) + a*JVS( 157 )
  W( 51 ) = W( 51 ) + a*JVS( 158 )
  W( 56 ) = W( 56 ) + a*JVS( 159 )
  W( 61 ) = W( 61 ) + a*JVS( 160 )
  W( 64 ) = W( 64 ) + a*JVS( 161 )
  W( 73 ) = W( 73 ) + a*JVS( 162 )
  W( 74 ) = W( 74 ) + a*JVS( 163 )
  W( 75 ) = W( 75 ) + a*JVS( 164 )
  W( 80 ) = W( 80 ) + a*JVS( 165 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 173 )
  W( 50 ) = W( 50 ) + a*JVS( 174 )
  W( 73 ) = W( 73 ) + a*JVS( 175 )
  W( 74 ) = W( 74 ) + a*JVS( 176 )
  W( 79 ) = W( 79 ) + a*JVS( 177 )
  a = -W( 35 ) / JVS(          178  )
  W( 35 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 179 )
  W( 72 ) = W( 72 ) + a*JVS( 180 )
  W( 73 ) = W( 73 ) + a*JVS( 181 )
  W( 74 ) = W( 74 ) + a*JVS( 182 )
  W( 77 ) = W( 77 ) + a*JVS( 183 )
  W( 80 ) = W( 80 ) + a*JVS( 184 )
  a = -W( 39 ) / JVS(          217  )
  W( 39 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 218 )
  W( 74 ) = W( 74 ) + a*JVS( 219 )
  W( 77 ) = W( 77 ) + a*JVS( 220 )
  a = -W( 40 ) / JVS(          225  )
  W( 40 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 226 )
  W( 74 ) = W( 74 ) + a*JVS( 227 )
  W( 75 ) = W( 75 ) + a*JVS( 228 )
  W( 77 ) = W( 77 ) + a*JVS( 229 )
  W( 80 ) = W( 80 ) + a*JVS( 230 )
  a = -W( 41 ) / JVS(          231  )
  W( 41 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 232 )
  W( 74 ) = W( 74 ) + a*JVS( 233 )
  W( 77 ) = W( 77 ) + a*JVS( 234 )
  a = -W( 42 ) / JVS(          238  )
  W( 42 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 239 )
  W( 73 ) = W( 73 ) + a*JVS( 240 )
  W( 74 ) = W( 74 ) + a*JVS( 241 )
  W( 75 ) = W( 75 ) + a*JVS( 242 )
  W( 77 ) = W( 77 ) + a*JVS( 243 )
  W( 79 ) = W( 79 ) + a*JVS( 244 )
  W( 80 ) = W( 80 ) + a*JVS( 245 )
  a = -W( 43 ) / JVS(          246  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 247 )
  W( 74 ) = W( 74 ) + a*JVS( 248 )
  W( 77 ) = W( 77 ) + a*JVS( 249 )
  W( 78 ) = W( 78 ) + a*JVS( 250 )
  W( 80 ) = W( 80 ) + a*JVS( 251 )
  a = -W( 44 ) / JVS(          253  )
  W( 44 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 254 )
  W( 57 ) = W( 57 ) + a*JVS( 255 )
  W( 58 ) = W( 58 ) + a*JVS( 256 )
  W( 61 ) = W( 61 ) + a*JVS( 257 )
  W( 71 ) = W( 71 ) + a*JVS( 258 )
  W( 72 ) = W( 72 ) + a*JVS( 259 )
  W( 73 ) = W( 73 ) + a*JVS( 260 )
  W( 74 ) = W( 74 ) + a*JVS( 261 )
  W( 76 ) = W( 76 ) + a*JVS( 262 )
  W( 77 ) = W( 77 ) + a*JVS( 263 )
  W( 79 ) = W( 79 ) + a*JVS( 264 )
  a = -W( 45 ) / JVS(          266  )
  W( 45 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 267 )
  W( 73 ) = W( 73 ) + a*JVS( 268 )
  W( 74 ) = W( 74 ) + a*JVS( 269 )
  W( 75 ) = W( 75 ) + a*JVS( 270 )
  W( 76 ) = W( 76 ) + a*JVS( 271 )
  W( 77 ) = W( 77 ) + a*JVS( 272 )
  W( 79 ) = W( 79 ) + a*JVS( 273 )
  a = -W( 46 ) / JVS(          274  )
  W( 46 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 275 )
  W( 72 ) = W( 72 ) + a*JVS( 276 )
  W( 73 ) = W( 73 ) + a*JVS( 277 )
  W( 74 ) = W( 74 ) + a*JVS( 278 )
  W( 75 ) = W( 75 ) + a*JVS( 279 )
  W( 76 ) = W( 76 ) + a*JVS( 280 )
  W( 79 ) = W( 79 ) + a*JVS( 281 )
  a = -W( 47 ) / JVS(          282  )
  W( 47 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 283 )
  W( 72 ) = W( 72 ) + a*JVS( 284 )
  W( 73 ) = W( 73 ) + a*JVS( 285 )
  W( 74 ) = W( 74 ) + a*JVS( 286 )
  W( 75 ) = W( 75 ) + a*JVS( 287 )
  W( 76 ) = W( 76 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  a = -W( 48 ) / JVS(          294  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 295 )
  W( 54 ) = W( 54 ) + a*JVS( 296 )
  W( 72 ) = W( 72 ) + a*JVS( 297 )
  W( 73 ) = W( 73 ) + a*JVS( 298 )
  W( 74 ) = W( 74 ) + a*JVS( 299 )
  W( 75 ) = W( 75 ) + a*JVS( 300 )
  W( 76 ) = W( 76 ) + a*JVS( 301 )
  W( 77 ) = W( 77 ) + a*JVS( 302 )
  W( 79 ) = W( 79 ) + a*JVS( 303 )
  a = -W( 50 ) / JVS(          331  )
  W( 50 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 332 )
  W( 73 ) = W( 73 ) + a*JVS( 333 )
  W( 74 ) = W( 74 ) + a*JVS( 334 )
  W( 75 ) = W( 75 ) + a*JVS( 335 )
  W( 76 ) = W( 76 ) + a*JVS( 336 )
  W( 77 ) = W( 77 ) + a*JVS( 337 )
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  a = -W( 51 ) / JVS(          347  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 348 )
  W( 54 ) = W( 54 ) + a*JVS( 349 )
  W( 57 ) = W( 57 ) + a*JVS( 350 )
  W( 58 ) = W( 58 ) + a*JVS( 351 )
  W( 61 ) = W( 61 ) + a*JVS( 352 )
  W( 69 ) = W( 69 ) + a*JVS( 353 )
  W( 71 ) = W( 71 ) + a*JVS( 354 )
  W( 72 ) = W( 72 ) + a*JVS( 355 )
  W( 73 ) = W( 73 ) + a*JVS( 356 )
  W( 74 ) = W( 74 ) + a*JVS( 357 )
  W( 75 ) = W( 75 ) + a*JVS( 358 )
  W( 76 ) = W( 76 ) + a*JVS( 359 )
  W( 77 ) = W( 77 ) + a*JVS( 360 )
  W( 78 ) = W( 78 ) + a*JVS( 361 )
  W( 79 ) = W( 79 ) + a*JVS( 362 )
  W( 80 ) = W( 80 ) + a*JVS( 363 )
  a = -W( 52 ) / JVS(          366  )
  W( 52 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 367 )
  W( 73 ) = W( 73 ) + a*JVS( 368 )
  W( 74 ) = W( 74 ) + a*JVS( 369 )
  W( 75 ) = W( 75 ) + a*JVS( 370 )
  W( 76 ) = W( 76 ) + a*JVS( 371 )
  W( 77 ) = W( 77 ) + a*JVS( 372 )
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 80 ) = W( 80 ) + a*JVS( 374 )
  a = -W( 53 ) / JVS(          376  )
  W( 53 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 377 )
  W( 73 ) = W( 73 ) + a*JVS( 378 )
  W( 74 ) = W( 74 ) + a*JVS( 379 )
  W( 75 ) = W( 75 ) + a*JVS( 380 )
  W( 76 ) = W( 76 ) + a*JVS( 381 )
  W( 79 ) = W( 79 ) + a*JVS( 382 )
  a = -W( 54 ) / JVS(          384  )
  W( 54 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 385 )
  W( 73 ) = W( 73 ) + a*JVS( 386 )
  W( 74 ) = W( 74 ) + a*JVS( 387 )
  W( 75 ) = W( 75 ) + a*JVS( 388 )
  W( 76 ) = W( 76 ) + a*JVS( 389 )
  W( 77 ) = W( 77 ) + a*JVS( 390 )
  W( 79 ) = W( 79 ) + a*JVS( 391 )
  a = -W( 55 ) / JVS(          393  )
  W( 55 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 394 )
  W( 73 ) = W( 73 ) + a*JVS( 395 )
  W( 74 ) = W( 74 ) + a*JVS( 396 )
  W( 75 ) = W( 75 ) + a*JVS( 397 )
  W( 76 ) = W( 76 ) + a*JVS( 398 )
  W( 79 ) = W( 79 ) + a*JVS( 399 )
  a = -W( 56 ) / JVS(          423  )
  W( 56 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 424 )
  W( 58 ) = W( 58 ) + a*JVS( 425 )
  W( 59 ) = W( 59 ) + a*JVS( 426 )
  W( 61 ) = W( 61 ) + a*JVS( 427 )
  W( 62 ) = W( 62 ) + a*JVS( 428 )
  W( 63 ) = W( 63 ) + a*JVS( 429 )
  W( 65 ) = W( 65 ) + a*JVS( 430 )
  W( 66 ) = W( 66 ) + a*JVS( 431 )
  W( 67 ) = W( 67 ) + a*JVS( 432 )
  W( 68 ) = W( 68 ) + a*JVS( 433 )
  W( 69 ) = W( 69 ) + a*JVS( 434 )
  W( 71 ) = W( 71 ) + a*JVS( 435 )
  W( 72 ) = W( 72 ) + a*JVS( 436 )
  W( 73 ) = W( 73 ) + a*JVS( 437 )
  W( 74 ) = W( 74 ) + a*JVS( 438 )
  W( 75 ) = W( 75 ) + a*JVS( 439 )
  W( 76 ) = W( 76 ) + a*JVS( 440 )
  W( 77 ) = W( 77 ) + a*JVS( 441 )
  W( 78 ) = W( 78 ) + a*JVS( 442 )
  W( 79 ) = W( 79 ) + a*JVS( 443 )
  W( 80 ) = W( 80 ) + a*JVS( 444 )
  a = -W( 57 ) / JVS(          446  )
  W( 57 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 447 )
  W( 73 ) = W( 73 ) + a*JVS( 448 )
  W( 74 ) = W( 74 ) + a*JVS( 449 )
  W( 75 ) = W( 75 ) + a*JVS( 450 )
  W( 76 ) = W( 76 ) + a*JVS( 451 )
  W( 77 ) = W( 77 ) + a*JVS( 452 )
  W( 79 ) = W( 79 ) + a*JVS( 453 )
  W( 80 ) = W( 80 ) + a*JVS( 454 )
  a = -W( 58 ) / JVS(          456  )
  W( 58 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 457 )
  W( 73 ) = W( 73 ) + a*JVS( 458 )
  W( 74 ) = W( 74 ) + a*JVS( 459 )
  W( 75 ) = W( 75 ) + a*JVS( 460 )
  W( 76 ) = W( 76 ) + a*JVS( 461 )
  W( 77 ) = W( 77 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  a = -W( 59 ) / JVS(          466  )
  W( 59 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 467 )
  W( 73 ) = W( 73 ) + a*JVS( 468 )
  W( 74 ) = W( 74 ) + a*JVS( 469 )
  W( 75 ) = W( 75 ) + a*JVS( 470 )
  W( 76 ) = W( 76 ) + a*JVS( 471 )
  W( 77 ) = W( 77 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  a = -W( 60 ) / JVS(          483  )
  W( 60 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 484 )
  W( 67 ) = W( 67 ) + a*JVS( 485 )
  W( 71 ) = W( 71 ) + a*JVS( 486 )
  W( 72 ) = W( 72 ) + a*JVS( 487 )
  W( 73 ) = W( 73 ) + a*JVS( 488 )
  W( 74 ) = W( 74 ) + a*JVS( 489 )
  W( 75 ) = W( 75 ) + a*JVS( 490 )
  W( 76 ) = W( 76 ) + a*JVS( 491 )
  W( 77 ) = W( 77 ) + a*JVS( 492 )
  W( 78 ) = W( 78 ) + a*JVS( 493 )
  W( 79 ) = W( 79 ) + a*JVS( 494 )
  W( 80 ) = W( 80 ) + a*JVS( 495 )
  a = -W( 61 ) / JVS(          498  )
  W( 61 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 499 )
  W( 73 ) = W( 73 ) + a*JVS( 500 )
  W( 74 ) = W( 74 ) + a*JVS( 501 )
  W( 75 ) = W( 75 ) + a*JVS( 502 )
  W( 76 ) = W( 76 ) + a*JVS( 503 )
  W( 77 ) = W( 77 ) + a*JVS( 504 )
  W( 79 ) = W( 79 ) + a*JVS( 505 )
  W( 80 ) = W( 80 ) + a*JVS( 506 )
  a = -W( 62 ) / JVS(          517  )
  W( 62 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 518 )
  W( 70 ) = W( 70 ) + a*JVS( 519 )
  W( 71 ) = W( 71 ) + a*JVS( 520 )
  W( 72 ) = W( 72 ) + a*JVS( 521 )
  W( 73 ) = W( 73 ) + a*JVS( 522 )
  W( 74 ) = W( 74 ) + a*JVS( 523 )
  W( 75 ) = W( 75 ) + a*JVS( 524 )
  W( 76 ) = W( 76 ) + a*JVS( 525 )
  W( 77 ) = W( 77 ) + a*JVS( 526 )
  W( 78 ) = W( 78 ) + a*JVS( 527 )
  W( 79 ) = W( 79 ) + a*JVS( 528 )
  W( 80 ) = W( 80 ) + a*JVS( 529 )
  a = -W( 63 ) / JVS(          533  )
  W( 63 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 534 )
  W( 73 ) = W( 73 ) + a*JVS( 535 )
  W( 74 ) = W( 74 ) + a*JVS( 536 )
  W( 75 ) = W( 75 ) + a*JVS( 537 )
  W( 76 ) = W( 76 ) + a*JVS( 538 )
  W( 77 ) = W( 77 ) + a*JVS( 539 )
  W( 79 ) = W( 79 ) + a*JVS( 540 )
  a = -W( 64 ) / JVS(          554  )
  W( 64 ) = -a
  W( 65 ) = W( 65 ) + a*JVS( 555 )
  W( 66 ) = W( 66 ) + a*JVS( 556 )
  W( 67 ) = W( 67 ) + a*JVS( 557 )
  W( 68 ) = W( 68 ) + a*JVS( 558 )
  W( 69 ) = W( 69 ) + a*JVS( 559 )
  W( 70 ) = W( 70 ) + a*JVS( 560 )
  W( 71 ) = W( 71 ) + a*JVS( 561 )
  W( 72 ) = W( 72 ) + a*JVS( 562 )
  W( 73 ) = W( 73 ) + a*JVS( 563 )
  W( 74 ) = W( 74 ) + a*JVS( 564 )
  W( 75 ) = W( 75 ) + a*JVS( 565 )
  W( 76 ) = W( 76 ) + a*JVS( 566 )
  W( 77 ) = W( 77 ) + a*JVS( 567 )
  W( 78 ) = W( 78 ) + a*JVS( 568 )
  W( 79 ) = W( 79 ) + a*JVS( 569 )
  W( 80 ) = W( 80 ) + a*JVS( 570 )
  a = -W( 65 ) / JVS(          573  )
  W( 65 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 574 )
  W( 73 ) = W( 73 ) + a*JVS( 575 )
  W( 74 ) = W( 74 ) + a*JVS( 576 )
  W( 75 ) = W( 75 ) + a*JVS( 577 )
  W( 76 ) = W( 76 ) + a*JVS( 578 )
  W( 77 ) = W( 77 ) + a*JVS( 579 )
  W( 79 ) = W( 79 ) + a*JVS( 580 )
  a = -W( 66 ) / JVS(          587  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 588 )
  W( 72 ) = W( 72 ) + a*JVS( 589 )
  W( 73 ) = W( 73 ) + a*JVS( 590 )
  W( 74 ) = W( 74 ) + a*JVS( 591 )
  W( 75 ) = W( 75 ) + a*JVS( 592 )
  W( 76 ) = W( 76 ) + a*JVS( 593 )
  W( 77 ) = W( 77 ) + a*JVS( 594 )
  W( 79 ) = W( 79 ) + a*JVS( 595 )
  a = -W( 67 ) / JVS(          603  )
  W( 67 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 604 )
  W( 73 ) = W( 73 ) + a*JVS( 605 )
  W( 74 ) = W( 74 ) + a*JVS( 606 )
  W( 75 ) = W( 75 ) + a*JVS( 607 )
  W( 76 ) = W( 76 ) + a*JVS( 608 )
  W( 77 ) = W( 77 ) + a*JVS( 609 )
  W( 79 ) = W( 79 ) + a*JVS( 610 )
  a = -W( 68 ) / JVS(          624  )
  W( 68 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 625 )
  W( 70 ) = W( 70 ) + a*JVS( 626 )
  W( 71 ) = W( 71 ) + a*JVS( 627 )
  W( 72 ) = W( 72 ) + a*JVS( 628 )
  W( 73 ) = W( 73 ) + a*JVS( 629 )
  W( 74 ) = W( 74 ) + a*JVS( 630 )
  W( 75 ) = W( 75 ) + a*JVS( 631 )
  W( 76 ) = W( 76 ) + a*JVS( 632 )
  W( 77 ) = W( 77 ) + a*JVS( 633 )
  W( 78 ) = W( 78 ) + a*JVS( 634 )
  W( 79 ) = W( 79 ) + a*JVS( 635 )
  W( 80 ) = W( 80 ) + a*JVS( 636 )
  a = -W( 69 ) / JVS(          650  )
  W( 69 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 651 )
  W( 71 ) = W( 71 ) + a*JVS( 652 )
  W( 72 ) = W( 72 ) + a*JVS( 653 )
  W( 73 ) = W( 73 ) + a*JVS( 654 )
  W( 74 ) = W( 74 ) + a*JVS( 655 )
  W( 75 ) = W( 75 ) + a*JVS( 656 )
  W( 76 ) = W( 76 ) + a*JVS( 657 )
  W( 77 ) = W( 77 ) + a*JVS( 658 )
  W( 78 ) = W( 78 ) + a*JVS( 659 )
  W( 79 ) = W( 79 ) + a*JVS( 660 )
  W( 80 ) = W( 80 ) + a*JVS( 661 )
  a = -W( 70 ) / JVS(          677  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 678 )
  W( 72 ) = W( 72 ) + a*JVS( 679 )
  W( 73 ) = W( 73 ) + a*JVS( 680 )
  W( 74 ) = W( 74 ) + a*JVS( 681 )
  W( 75 ) = W( 75 ) + a*JVS( 682 )
  W( 76 ) = W( 76 ) + a*JVS( 683 )
  W( 77 ) = W( 77 ) + a*JVS( 684 )
  W( 78 ) = W( 78 ) + a*JVS( 685 )
  W( 79 ) = W( 79 ) + a*JVS( 686 )
  W( 80 ) = W( 80 ) + a*JVS( 687 )
  a = -W( 71 ) / JVS(          693  )
  W( 71 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 694 )
  W( 73 ) = W( 73 ) + a*JVS( 695 )
  W( 74 ) = W( 74 ) + a*JVS( 696 )
  W( 75 ) = W( 75 ) + a*JVS( 697 )
  W( 76 ) = W( 76 ) + a*JVS( 698 )
  W( 77 ) = W( 77 ) + a*JVS( 699 )
  W( 78 ) = W( 78 ) + a*JVS( 700 )
  W( 79 ) = W( 79 ) + a*JVS( 701 )
  W( 80 ) = W( 80 ) + a*JVS( 702 )
  a = -W( 72 ) / JVS(          736  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 737 )
  W( 74 ) = W( 74 ) + a*JVS( 738 )
  W( 75 ) = W( 75 ) + a*JVS( 739 )
  W( 76 ) = W( 76 ) + a*JVS( 740 )
  W( 77 ) = W( 77 ) + a*JVS( 741 )
  W( 78 ) = W( 78 ) + a*JVS( 742 )
  W( 79 ) = W( 79 ) + a*JVS( 743 )
  W( 80 ) = W( 80 ) + a*JVS( 744 )
  a = -W( 73 ) / JVS(          806  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 807 )
  W( 75 ) = W( 75 ) + a*JVS( 808 )
  W( 76 ) = W( 76 ) + a*JVS( 809 )
  W( 77 ) = W( 77 ) + a*JVS( 810 )
  W( 78 ) = W( 78 ) + a*JVS( 811 )
  W( 79 ) = W( 79 ) + a*JVS( 812 )
  W( 80 ) = W( 80 ) + a*JVS( 813 )
  a = -W( 74 ) / JVS(          852  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 853 )
  W( 76 ) = W( 76 ) + a*JVS( 854 )
  W( 77 ) = W( 77 ) + a*JVS( 855 )
  W( 78 ) = W( 78 ) + a*JVS( 856 )
  W( 79 ) = W( 79 ) + a*JVS( 857 )
  W( 80 ) = W( 80 ) + a*JVS( 858 )
  a = -W( 75 ) / JVS(          916  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 917 )
  W( 77 ) = W( 77 ) + a*JVS( 918 )
  W( 78 ) = W( 78 ) + a*JVS( 919 )
  W( 79 ) = W( 79 ) + a*JVS( 920 )
  W( 80 ) = W( 80 ) + a*JVS( 921 )
  a = -W( 76 ) / JVS(          953  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 954 )
  W( 78 ) = W( 78 ) + a*JVS( 955 )
  W( 79 ) = W( 79 ) + a*JVS( 956 )
  W( 80 ) = W( 80 ) + a*JVS( 957 )
  a = -W( 77 ) / JVS(          980  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 981 )
  W( 79 ) = W( 79 ) + a*JVS( 982 )
  W( 80 ) = W( 80 ) + a*JVS( 983 )
  a = -W( 78 ) / JVS(          992  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 993 )
  W( 80 ) = W( 80 ) + a*JVS( 994 )
  a = -W( 79 ) / JVS(         1031  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 1032 )
  JVS( 1033) = W( 11 )
  JVS( 1034) = W( 18 )
  JVS( 1035) = W( 19 )
  JVS( 1036) = W( 21 )
  JVS( 1037) = W( 23 )
  JVS( 1038) = W( 26 )
  JVS( 1039) = W( 27 )
  JVS( 1040) = W( 28 )
  JVS( 1041) = W( 29 )
  JVS( 1042) = W( 31 )
  JVS( 1043) = W( 32 )
  JVS( 1044) = W( 34 )
  JVS( 1045) = W( 35 )
  JVS( 1046) = W( 39 )
  JVS( 1047) = W( 40 )
  JVS( 1048) = W( 41 )
  JVS( 1049) = W( 42 )
  JVS( 1050) = W( 43 )
  JVS( 1051) = W( 44 )
  JVS( 1052) = W( 45 )
  JVS( 1053) = W( 46 )
  JVS( 1054) = W( 47 )
  JVS( 1055) = W( 48 )
  JVS( 1056) = W( 50 )
  JVS( 1057) = W( 51 )
  JVS( 1058) = W( 52 )
  JVS( 1059) = W( 53 )
  JVS( 1060) = W( 54 )
  JVS( 1061) = W( 55 )
  JVS( 1062) = W( 56 )
  JVS( 1063) = W( 57 )
  JVS( 1064) = W( 58 )
  JVS( 1065) = W( 59 )
  JVS( 1066) = W( 60 )
  JVS( 1067) = W( 61 )
  JVS( 1068) = W( 62 )
  JVS( 1069) = W( 63 )
  JVS( 1070) = W( 64 )
  JVS( 1071) = W( 65 )
  JVS( 1072) = W( 66 )
  JVS( 1073) = W( 67 )
  JVS( 1074) = W( 68 )
  JVS( 1075) = W( 69 )
  JVS( 1076) = W( 70 )
  JVS( 1077) = W( 71 )
  JVS( 1078) = W( 72 )
  JVS( 1079) = W( 73 )
  JVS( 1080) = W( 74 )
  JVS( 1081) = W( 75 )
  JVS( 1082) = W( 76 )
  JVS( 1083) = W( 77 )
  JVS( 1084) = W( 78 )
  JVS( 1085) = W( 79 )
  JVS( 1086) = W( 80 )
   
   END SUBROUTINE decomp_racm_esrlsorg
 


END MODULE racm_esrlsorg_Integrator
