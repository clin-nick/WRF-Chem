
































MODULE radm2sorg_Integrator

 USE radm2sorg_Parameters
 USE radm2sorg_Precision
 USE radm2sorg_JacobianSP

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

SUBROUTINE  radm2sorg_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE radm2sorg_Parameters

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

   CALL radm2sorg_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  radm2sorg_INTEGRATE


SUBROUTINE  radm2sorg_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE radm2sorg_Parameters

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
      CALL radm2sorg_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL radm2sorg_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = radm2sorg_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL radm2sorg_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL radm2sorg_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL radm2sorg_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL radm2sorg_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL radm2sorg_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL radm2sorg_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL radm2sorg_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL radm2sorg_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL radm2sorg_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL radm2sorg_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL radm2sorg_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL radm2sorg_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL radm2sorg_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL radm2sorg_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL radm2sorg_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  radm2sorg_ros_ErrorMsg(Code,T,H,IERR)



   USE radm2sorg_Precision

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

 END SUBROUTINE  radm2sorg_ros_ErrorMsg


 SUBROUTINE  radm2sorg_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL radm2sorg_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL radm2sorg_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL radm2sorg_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL radm2sorg_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL radm2sorg_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL radm2sorg_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL radm2sorg_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL radm2sorg_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL radm2sorg_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL radm2sorg_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL radm2sorg_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL radm2sorg_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL radm2sorg_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL radm2sorg_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL radm2sorg_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL radm2sorg_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL radm2sorg_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL radm2sorg_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL radm2sorg_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL radm2sorg_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = radm2sorg_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL radm2sorg_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  radm2sorg_ros_Integrator



  REAL(kind=dp) FUNCTION  radm2sorg_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    radm2sorg_ros_ErrorNorm = Err

  END FUNCTION  radm2sorg_ros_ErrorNorm



  SUBROUTINE radm2sorg_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL radm2sorg_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL radm2sorg_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL radm2sorg_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  radm2sorg_ros_FunTimeDeriv



  SUBROUTINE  radm2sorg_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL radm2sorg_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL radm2sorg_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL radm2sorg_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  radm2sorg_ros_PrepareMatrix



  SUBROUTINE  radm2sorg_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 



CALL decomp_radm2sorg ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  radm2sorg_ros_Decomp



  SUBROUTINE  radm2sorg_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL radm2sorg_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  radm2sorg_ros_Solve




  SUBROUTINE  radm2sorg_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  radm2sorg_Ros2



  SUBROUTINE  radm2sorg_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  radm2sorg_Ros3





  SUBROUTINE  radm2sorg_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  radm2sorg_Ros4


  SUBROUTINE  radm2sorg_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  radm2sorg_Rodas3


  SUBROUTINE  radm2sorg_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  radm2sorg_Rodas4




END SUBROUTINE  radm2sorg_Rosenbrock




SUBROUTINE  radm2sorg_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE radm2sorg_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL radm2sorg_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  radm2sorg_FunTemplate



SUBROUTINE  radm2sorg_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE radm2sorg_Parameters
 
 USE radm2sorg_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL radm2sorg_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  radm2sorg_JacTemplate

















SUBROUTINE radm2sorg_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(59)
  A(2) = RCT(2)*V(46)
  A(3) = RCT(3)*V(46)
  A(4) = RCT(4)*V(13)
  A(5) = RCT(5)*V(25)
  A(6) = RCT(6)*V(20)
  A(7) = RCT(7)*V(54)
  A(8) = RCT(8)*V(54)
  A(9) = RCT(9)*V(14)
  A(10) = RCT(10)*V(43)
  A(11) = RCT(11)*V(43)
  A(12) = RCT(12)*V(48)
  A(13) = RCT(13)*V(21)
  A(14) = RCT(14)*V(42)
  A(15) = RCT(15)*V(18)
  A(16) = RCT(16)*V(33)
  A(17) = RCT(17)*V(31)
  A(18) = RCT(18)*V(31)
  A(19) = RCT(19)*V(34)
  A(20) = RCT(20)*V(30)
  A(21) = RCT(21)*V(47)
  A(22) = RCT(22)*V(19)*F(2)
  A(23) = RCT(23)*V(19)*V(59)
  A(24) = RCT(24)*V(7)*F(2)
  A(25) = RCT(25)*V(7)*F(1)
  A(26) = RCT(26)*V(46)*V(58)
  A(27) = RCT(27)*V(46)*V(53)
  A(28) = RCT(28)*V(46)*V(51)
  A(29) = RCT(29)*V(51)*V(58)
  A(30) = RCT(30)*V(51)*V(59)
  A(31) = RCT(31)*V(20)
  A(32) = RCT(32)*V(51)*V(51)
  A(33) = RCT(33)*V(51)*V(51)*F(1)
  A(34) = RCT(34)*V(14)*V(53)
  A(35) = RCT(35)*V(53)*V(58)
  A(36) = RCT(36)*V(58)*V(58)*F(2)
  A(37) = RCT(37)*V(46)*V(59)
  A(38) = RCT(38)*V(54)*V(58)
  A(39) = RCT(39)*V(54)*V(59)
  A(40) = RCT(40)*V(51)*V(54)
  A(41) = RCT(41)*V(54)*V(59)
  A(42) = RCT(42)*V(15)
  A(43) = RCT(43)*V(15)
  A(44) = RCT(44)*V(53)*V(59)
  A(45) = RCT(45)*V(25)*V(53)
  A(46) = RCT(46)*V(20)*V(53)
  A(47) = RCT(47)*V(51)*V(53)
  A(48) = RCT(48)*V(5)*V(53)
  A(49) = RCT(49)*V(26)*V(53)
  A(50) = RCT(50)*V(17)*V(53)
  A(51) = RCT(51)*V(6)*V(53)
  A(52) = RCT(52)*V(16)*V(53)
  A(53) = RCT(53)*V(8)*V(53)
  A(54) = RCT(54)*V(9)*V(53)
  A(55) = RCT(55)*V(24)*V(53)
  A(56) = RCT(56)*V(28)*V(53)
  A(57) = RCT(57)*V(29)*V(53)
  A(58) = RCT(58)*V(10)*V(53)
  A(59) = RCT(59)*V(11)*V(53)
  A(60) = RCT(60)*V(22)*V(53)
  A(61) = RCT(61)*V(43)*V(53)
  A(62) = RCT(62)*V(48)*V(53)
  A(63) = RCT(63)*V(33)*V(53)
  A(64) = RCT(64)*V(31)*V(53)
  A(65) = RCT(65)*V(34)*V(53)
  A(66) = RCT(66)*V(30)*V(53)
  A(67) = RCT(67)*V(21)*V(53)
  A(68) = RCT(68)*V(42)*V(53)
  A(69) = RCT(69)*V(18)*V(53)
  A(70) = RCT(70)*V(23)*V(53)
  A(71) = RCT(71)*V(47)*V(53)
  A(72) = RCT(72)*V(27)*V(53)
  A(73) = RCT(73)*V(55)*V(59)
  A(74) = RCT(74)*V(23)
  A(75) = RCT(75)*V(45)*V(59)
  A(76) = RCT(76)*V(12)
  A(77) = RCT(77)*V(52)*V(58)
  A(78) = RCT(78)*V(56)*V(58)
  A(79) = RCT(79)*V(41)*V(58)
  A(80) = RCT(80)*V(44)*V(58)
  A(81) = RCT(81)*V(40)*V(58)
  A(82) = RCT(82)*V(37)*V(58)
  A(83) = RCT(83)*V(49)*V(58)
  A(84) = RCT(84)*V(55)*V(58)
  A(85) = RCT(85)*V(45)*V(58)
  A(86) = RCT(86)*V(35)*V(58)
  A(87) = RCT(87)*V(36)*V(58)
  A(88) = RCT(88)*V(57)*V(58)
  A(89) = RCT(89)*V(50)*V(58)
  A(90) = RCT(90)*V(38)*V(58)
  A(91) = RCT(91)*V(43)*V(54)
  A(92) = RCT(92)*V(48)*V(54)
  A(93) = RCT(93)*V(31)*V(54)
  A(94) = RCT(94)*V(34)*V(54)
  A(95) = RCT(95)*V(30)*V(54)
  A(96) = RCT(96)*V(22)*V(54)
  A(97) = RCT(97)*V(24)*V(54)
  A(98) = RCT(98)*V(28)*V(54)
  A(99) = RCT(99)*V(29)*V(54)
  A(100) = RCT(100)*V(27)*V(54)
  A(101) = RCT(101)*V(24)*V(46)
  A(102) = RCT(102)*V(28)*V(46)
  A(103) = RCT(103)*V(29)*V(46)
  A(104) = RCT(104)*V(27)*V(46)
  A(105) = RCT(105)*V(51)*V(52)
  A(106) = RCT(106)*V(51)*V(57)
  A(107) = RCT(107)*V(51)*V(56)
  A(108) = RCT(108)*V(41)*V(51)
  A(109) = RCT(109)*V(44)*V(51)
  A(110) = RCT(110)*V(40)*V(51)
  A(111) = RCT(111)*V(37)*V(51)
  A(112) = RCT(112)*V(49)*V(51)
  A(113) = RCT(113)*V(50)*V(51)
  A(114) = RCT(114)*V(51)*V(55)
  A(115) = RCT(115)*V(35)*V(51)
  A(116) = RCT(116)*V(36)*V(51)
  A(117) = RCT(117)*V(45)*V(51)
  A(118) = RCT(118)*V(38)*V(51)
  A(119) = RCT(119)*V(52)*V(52)
  A(120) = RCT(120)*V(52)*V(57)
  A(121) = RCT(121)*V(52)*V(56)
  A(122) = RCT(122)*V(41)*V(52)
  A(123) = RCT(123)*V(44)*V(52)
  A(124) = RCT(124)*V(40)*V(52)
  A(125) = RCT(125)*V(37)*V(52)
  A(126) = RCT(126)*V(49)*V(52)
  A(127) = RCT(127)*V(50)*V(52)
  A(128) = RCT(128)*V(52)*V(55)
  A(129) = RCT(129)*V(35)*V(52)
  A(130) = RCT(130)*V(36)*V(52)
  A(131) = RCT(131)*V(45)*V(52)
  A(132) = RCT(132)*V(55)*V(57)
  A(133) = RCT(133)*V(55)*V(56)
  A(134) = RCT(134)*V(41)*V(55)
  A(135) = RCT(135)*V(44)*V(55)
  A(136) = RCT(136)*V(40)*V(55)
  A(137) = RCT(137)*V(37)*V(55)
  A(138) = RCT(138)*V(49)*V(55)
  A(139) = RCT(139)*V(50)*V(55)
  A(140) = RCT(140)*V(55)*V(55)
  A(141) = RCT(141)*V(35)*V(55)
  A(142) = RCT(142)*V(36)*V(55)
  A(143) = RCT(143)*V(45)*V(55)
  A(144) = RCT(144)*V(39)*V(51)
  A(145) = RCT(145)*V(39)*V(52)
  A(146) = RCT(146)*V(39)*V(55)
  A(147) = RCT(147)*V(39)*V(39)
  A(148) = RCT(148)*V(39)*V(58)
  A(149) = RCT(149)*V(32)*V(59)
  A(150) = RCT(150)*V(32)*V(51)
  A(151) = RCT(151)*V(32)*V(52)
  A(152) = RCT(152)*V(32)*V(55)
  A(153) = RCT(153)*V(32)*V(32)
  A(154) = RCT(154)*V(38)*V(52)
  A(155) = RCT(155)*V(38)*V(55)
  A(156) = RCT(156)*V(38)*V(38)


  Vdot(1) = A(48)
  Vdot(2) = 0.4*A(101)+0.2*A(102)+0.06*A(103)+0.2*A(104)
  Vdot(3) = 0.2*A(102)+0.29*A(103)+0.2*A(104)+0.5*A(128)+0.5*A(131)+0.5*A(132)+0.5*A(133)+0.5*A(134)+0.5*A(135)+0.5&
              &*A(136)+0.5*A(137)+0.5*A(138)+0.5*A(139)+0.5*A(155)
  Vdot(4) = A(49)
  Vdot(5) = -A(48)
  Vdot(6) = -A(51)
  Vdot(7) = A(2)-A(24)-A(25)
  Vdot(8) = -A(53)
  Vdot(9) = -A(54)
  Vdot(10) = -A(58)
  Vdot(11) = -A(59)
  Vdot(12) = A(75)-A(76)
  Vdot(13) = -A(4)+A(35)
  Vdot(14) = -A(9)+A(32)+A(33)-A(34)
  Vdot(15) = A(41)-A(42)-A(43)
  Vdot(16) = -A(52)
  Vdot(17) = -A(50)+0.06*A(102)+0.09*A(103)
  Vdot(18) = -A(15)-A(69)+A(114)
  Vdot(19) = A(1)+A(3)+A(8)-A(22)-A(23)+A(24)
  Vdot(20) = -A(6)+A(30)-A(31)-A(46)
  Vdot(21) = -A(13)-A(67)+A(105)
  Vdot(22) = 0.25*A(58)+0.17*A(59)-A(60)-0.5*A(96)
  Vdot(23) = -A(70)+A(73)-A(74)
  Vdot(24) = -A(55)-A(97)-A(101)
  Vdot(25) = -A(5)+A(40)+2*A(43)+A(44)-A(45)+A(91)+A(92)+A(93)+A(94)+A(95)+A(96)
  Vdot(26) = A(10)+A(11)+A(12)+1.87*A(17)+1.55*A(18)+A(19)-A(49)+A(61)+2*A(64)+A(65)+0.95*A(85)+A(91)+2*A(93)+A(94)+0.42&
               &*A(101)+0.33*A(102)+0.23*A(103)+0.33*A(104)+0.475*A(131)+0.95*A(143)
  Vdot(27) = -A(72)-A(100)-A(104)
  Vdot(28) = -A(56)-A(98)-A(102)
  Vdot(29) = -A(57)-A(99)-A(103)
  Vdot(30) = -A(20)-A(66)+0.7*A(86)+0.806*A(87)-A(95)+0.7*A(129)+0.806*A(130)+A(141)+A(142)
  Vdot(31) = -A(17)-A(18)-A(64)+0.89*A(85)+0.16*A(86)-A(93)+0.16*A(129)+0.445*A(131)+0.2*A(141)+0.89*A(143)
  Vdot(32) = A(96)-A(149)-A(150)-A(151)-A(152)-2*A(153)
  Vdot(33) = -A(16)+0.8*A(21)+0.025*A(52)-A(63)+0.25*A(78)+0.69*A(79)+1.06*A(80)+0.1*A(83)+0.1*A(103)+0.6*A(121)+0.75&
               &*A(122)+1.39*A(123)+0.55*A(126)+0.8*A(133)+0.86*A(134)+0.9*A(135)+0.55*A(138)
  Vdot(34) = -A(19)-A(65)+0.11*A(85)+0.17*A(86)+0.45*A(87)+A(89)-A(94)+0.75*A(127)+0.17*A(129)+0.45*A(130)+0.055*A(131)&
               &+A(139)+0.8*A(141)+A(142)+0.11*A(143)
  Vdot(35) = 0.75*A(58)-A(86)-A(115)-A(129)-A(141)
  Vdot(36) = 0.83*A(59)-A(87)-A(116)-A(130)-A(142)
  Vdot(37) = A(56)+A(72)-A(82)-A(111)-A(125)-A(137)
  Vdot(38) = -A(90)+A(97)+A(98)+A(99)+A(100)-A(118)-A(154)-A(155)-2*A(156)
  Vdot(39) = 0.25*A(53)+0.75*A(54)+0.9*A(60)+A(70)+2*A(85)+A(131)+2*A(143)-A(144)-A(145)-A(146)-2*A(147)-A(148)
  Vdot(40) = A(55)-A(81)-A(110)-A(124)-A(136)
  Vdot(41) = A(53)-A(79)-A(108)-A(122)-A(134)
  Vdot(42) = -A(14)-A(68)+A(106)+A(107)+A(108)+A(109)+A(110)+A(111)+A(112)+A(113)+A(115)+A(116)+A(117)+A(144)+A(150)
  Vdot(43) = -A(10)-A(11)+A(13)+0.13*A(17)+0.45*A(18)+0.009*A(52)-A(61)+0.5*A(67)+A(70)+A(77)+0.09*A(78)+0.04*A(80)+1.6&
               &*A(81)+A(82)+0.28*A(83)+A(90)-A(91)+A(101)+0.53*A(102)+0.18*A(103)+0.53*A(104)+1.5*A(119)+0.75*A(120)+0.75&
               &*A(121)+0.77*A(122)+0.8*A(123)+1.55*A(124)+1.25*A(125)+0.89*A(126)+0.75*A(127)+A(128)+A(129)+A(130)+0.5&
               &*A(131)+0.8*A(136)+0.5*A(137)+0.14*A(138)+A(145)+A(151)+1.75*A(154)+A(155)+2*A(156)
  Vdot(44) = A(54)-A(80)-A(109)-A(123)-A(135)
  Vdot(45) = A(20)+0.9*A(60)+A(66)-A(75)+A(76)-A(85)+A(95)-A(117)-A(131)-A(143)
  Vdot(46) = -A(2)-A(3)+A(22)-A(26)-A(27)-A(28)-A(37)-A(101)-A(102)-A(103)-A(104)
  Vdot(47) = -A(21)-A(71)+0.036*A(78)+0.08*A(79)+0.24*A(80)+A(118)+A(149)
  Vdot(48) = -A(12)+A(14)+0.2*A(21)+0.075*A(52)-A(62)+0.5*A(68)+0.75*A(78)+0.38*A(79)+0.35*A(80)+0.2*A(81)+A(82)+1.45&
               &*A(83)+A(88)+A(90)-A(92)+0.5*A(102)+0.72*A(103)+0.5*A(104)+0.75*A(120)+0.15*A(121)+0.41*A(122)+0.46*A(123)&
               &+0.35*A(124)+0.75*A(125)+0.725*A(126)+A(132)+0.2*A(133)+0.14*A(134)+0.1*A(135)+0.6*A(136)+A(137)+0.725&
               &*A(138)+A(154)+A(155)+2*A(156)
  Vdot(49) = A(57)-A(83)-A(112)-A(126)-A(138)
  Vdot(50) = A(63)-A(89)-A(113)-A(127)-A(139)
  Vdot(51) = 0.65*A(6)+2*A(11)+A(12)+A(13)+A(14)+0.8*A(18)+A(19)+A(20)+A(21)+A(27)-A(28)-A(29)-A(30)+A(31)-2*A(32)-2&
               &*A(33)+A(34)-A(40)-A(47)+A(48)+A(49)+0.17*A(52)+0.25*A(58)+0.17*A(59)+0.1*A(60)+A(61)+A(64)+A(77)+0.964&
               &*A(78)+0.92*A(79)+0.76*A(80)+A(81)+A(82)+A(83)+0.92*A(85)+A(86)+A(87)+A(88)+A(89)+A(91)+A(93)+0.12*A(101)&
               &+0.23*A(102)+0.26*A(103)+0.23*A(104)-A(105)-A(106)-A(107)-A(108)-A(109)-A(110)-A(111)-A(112)-A(113)-A(114)&
               &-A(115)-A(116)-A(117)-A(118)+A(119)+A(120)+A(121)+A(122)+A(123)+A(124)+A(125)+A(126)+A(127)+0.5*A(128)+2&
               &*A(129)+2*A(130)+0.46*A(131)+0.5*A(132)+0.5*A(133)+0.5*A(134)+0.5*A(135)+0.5*A(136)+0.5*A(137)+0.5*A(138)&
               &+0.5*A(139)+A(141)+A(142)+0.92*A(143)-A(144)+A(145)-A(150)+A(151)+0.5*A(154)
  Vdot(52) = A(12)+A(15)+A(50)+0.5*A(67)-A(77)+A(84)+0.22*A(102)+0.31*A(103)+0.22*A(104)-A(105)-2*A(119)-A(120)-A(121)&
               &-A(122)-A(123)-A(124)-A(125)-A(126)-A(127)-0.5*A(128)-A(129)-A(130)-A(131)+0.5*A(132)+0.5*A(133)+0.5*A(134)&
               &+0.5*A(135)+0.5*A(136)+0.5*A(137)+0.5*A(138)+0.5*A(139)+2*A(140)+A(141)+A(142)+A(143)-A(145)+A(146)-A(151)&
               &+A(152)-A(154)+0.5*A(155)
  Vdot(53) = A(4)+A(5)+0.35*A(6)+2*A(9)+A(13)+A(14)+A(15)+2*A(25)-A(27)+A(28)+A(29)-A(34)-A(35)-A(44)-A(45)-A(46)-A(47)&
               &-A(48)-A(49)-A(50)-A(51)-A(52)-A(53)-A(54)-A(55)-A(56)-A(57)-A(58)-A(59)-1.9*A(60)-A(61)-A(62)-A(63)-A(64)&
               &-A(65)-A(66)-0.5*A(67)-0.5*A(68)-A(69)-A(70)-A(71)-A(72)+0.1*A(102)+0.14*A(103)+0.1*A(104)
  Vdot(54) = 0.35*A(6)-A(7)-A(8)+A(37)-A(38)-A(39)-A(40)-A(41)+A(42)+A(45)+A(70)-A(91)-A(92)-A(93)-A(94)-A(95)-A(96)&
               &-A(97)-A(98)-A(99)-A(100)
  Vdot(55) = A(16)+A(19)+A(62)+A(65)+A(69)-A(73)+A(74)-A(84)+0.05*A(85)+A(92)+A(94)-A(114)-A(128)+0.025*A(131)-A(132)&
               &-A(133)-A(134)-A(135)-A(136)-A(137)-A(138)-A(139)-2*A(140)-A(141)-A(142)-0.95*A(143)-A(146)-A(152)-A(155)
  Vdot(56) = 0.83*A(52)+0.5*A(68)+A(71)-A(78)-A(107)-A(121)-A(133)
  Vdot(57) = A(16)+A(51)-A(88)-A(106)-A(120)-A(132)
  Vdot(58) = A(1)+A(4)+A(7)+A(23)-A(26)-A(29)-A(35)-2*A(36)-A(38)+A(39)-A(77)-A(78)-A(79)-A(80)-A(81)-A(82)-A(83)-A(84)&
               &-A(85)-A(86)-A(87)-A(88)-A(89)-A(90)-A(148)
  Vdot(59) = -A(1)+A(5)+0.65*A(6)+A(8)+A(21)-A(23)+A(26)+A(29)-A(30)+A(31)+2*A(36)-A(37)+2*A(38)-A(41)+A(42)-A(44)+A(46)&
               &+A(71)-A(73)+A(74)-A(75)+A(76)+A(77)+0.964*A(78)+0.92*A(79)+0.76*A(80)+A(81)+A(82)+A(83)+A(84)+A(85)+A(86)&
               &+A(87)+A(88)+A(89)+2*A(90)+A(148)-A(149)+A(154)+A(155)+2*A(156)
      
END SUBROUTINE radm2sorg_Fun
















SUBROUTINE radm2sorg_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(59)
  IRR(2) = RCT(2)*V(46)
  IRR(3) = RCT(3)*V(46)
  IRR(4) = RCT(4)*V(13)
  IRR(5) = RCT(5)*V(25)
  IRR(6) = RCT(6)*V(20)
  IRR(7) = RCT(7)*V(54)
  IRR(8) = RCT(8)*V(54)
  IRR(9) = RCT(9)*V(14)
  IRR(10) = RCT(10)*V(43)
  IRR(11) = RCT(11)*V(43)
  IRR(12) = RCT(12)*V(48)
  IRR(13) = RCT(13)*V(21)
  IRR(14) = RCT(14)*V(42)
  IRR(15) = RCT(15)*V(18)
  IRR(16) = RCT(16)*V(33)
  IRR(17) = RCT(17)*V(31)
  IRR(18) = RCT(18)*V(31)
  IRR(19) = RCT(19)*V(34)
  IRR(20) = RCT(20)*V(30)
  IRR(21) = RCT(21)*V(47)
  IRR(22) = RCT(22)*V(19)*F(2)
  IRR(23) = RCT(23)*V(19)*V(59)
  IRR(24) = RCT(24)*V(7)*F(2)
  IRR(25) = RCT(25)*V(7)*F(1)
  IRR(26) = RCT(26)*V(46)*V(58)
  IRR(27) = RCT(27)*V(46)*V(53)
  IRR(28) = RCT(28)*V(46)*V(51)
  IRR(29) = RCT(29)*V(51)*V(58)
  IRR(30) = RCT(30)*V(51)*V(59)
  IRR(31) = RCT(31)*V(20)
  IRR(32) = RCT(32)*V(51)*V(51)
  IRR(33) = RCT(33)*V(51)*V(51)*F(1)
  IRR(34) = RCT(34)*V(14)*V(53)
  IRR(35) = RCT(35)*V(53)*V(58)
  IRR(36) = RCT(36)*V(58)*V(58)*F(2)
  IRR(37) = RCT(37)*V(46)*V(59)
  IRR(38) = RCT(38)*V(54)*V(58)
  IRR(39) = RCT(39)*V(54)*V(59)
  IRR(40) = RCT(40)*V(51)*V(54)
  IRR(41) = RCT(41)*V(54)*V(59)
  IRR(42) = RCT(42)*V(15)
  IRR(43) = RCT(43)*V(15)
  IRR(44) = RCT(44)*V(53)*V(59)
  IRR(45) = RCT(45)*V(25)*V(53)
  IRR(46) = RCT(46)*V(20)*V(53)
  IRR(47) = RCT(47)*V(51)*V(53)
  IRR(48) = RCT(48)*V(5)*V(53)
  IRR(49) = RCT(49)*V(26)*V(53)
  IRR(50) = RCT(50)*V(17)*V(53)
  IRR(51) = RCT(51)*V(6)*V(53)
  IRR(52) = RCT(52)*V(16)*V(53)
  IRR(53) = RCT(53)*V(8)*V(53)
  IRR(54) = RCT(54)*V(9)*V(53)
  IRR(55) = RCT(55)*V(24)*V(53)
  IRR(56) = RCT(56)*V(28)*V(53)
  IRR(57) = RCT(57)*V(29)*V(53)
  IRR(58) = RCT(58)*V(10)*V(53)
  IRR(59) = RCT(59)*V(11)*V(53)
  IRR(60) = RCT(60)*V(22)*V(53)
  IRR(61) = RCT(61)*V(43)*V(53)
  IRR(62) = RCT(62)*V(48)*V(53)
  IRR(63) = RCT(63)*V(33)*V(53)
  IRR(64) = RCT(64)*V(31)*V(53)
  IRR(65) = RCT(65)*V(34)*V(53)
  IRR(66) = RCT(66)*V(30)*V(53)
  IRR(67) = RCT(67)*V(21)*V(53)
  IRR(68) = RCT(68)*V(42)*V(53)
  IRR(69) = RCT(69)*V(18)*V(53)
  IRR(70) = RCT(70)*V(23)*V(53)
  IRR(71) = RCT(71)*V(47)*V(53)
  IRR(72) = RCT(72)*V(27)*V(53)
  IRR(73) = RCT(73)*V(55)*V(59)
  IRR(74) = RCT(74)*V(23)
  IRR(75) = RCT(75)*V(45)*V(59)
  IRR(76) = RCT(76)*V(12)
  IRR(77) = RCT(77)*V(52)*V(58)
  IRR(78) = RCT(78)*V(56)*V(58)
  IRR(79) = RCT(79)*V(41)*V(58)
  IRR(80) = RCT(80)*V(44)*V(58)
  IRR(81) = RCT(81)*V(40)*V(58)
  IRR(82) = RCT(82)*V(37)*V(58)
  IRR(83) = RCT(83)*V(49)*V(58)
  IRR(84) = RCT(84)*V(55)*V(58)
  IRR(85) = RCT(85)*V(45)*V(58)
  IRR(86) = RCT(86)*V(35)*V(58)
  IRR(87) = RCT(87)*V(36)*V(58)
  IRR(88) = RCT(88)*V(57)*V(58)
  IRR(89) = RCT(89)*V(50)*V(58)
  IRR(90) = RCT(90)*V(38)*V(58)
  IRR(91) = RCT(91)*V(43)*V(54)
  IRR(92) = RCT(92)*V(48)*V(54)
  IRR(93) = RCT(93)*V(31)*V(54)
  IRR(94) = RCT(94)*V(34)*V(54)
  IRR(95) = RCT(95)*V(30)*V(54)
  IRR(96) = RCT(96)*V(22)*V(54)
  IRR(97) = RCT(97)*V(24)*V(54)
  IRR(98) = RCT(98)*V(28)*V(54)
  IRR(99) = RCT(99)*V(29)*V(54)
  IRR(100) = RCT(100)*V(27)*V(54)
  IRR(101) = RCT(101)*V(24)*V(46)
  IRR(102) = RCT(102)*V(28)*V(46)
  IRR(103) = RCT(103)*V(29)*V(46)
  IRR(104) = RCT(104)*V(27)*V(46)
  IRR(105) = RCT(105)*V(51)*V(52)
  IRR(106) = RCT(106)*V(51)*V(57)
  IRR(107) = RCT(107)*V(51)*V(56)
  IRR(108) = RCT(108)*V(41)*V(51)
  IRR(109) = RCT(109)*V(44)*V(51)
  IRR(110) = RCT(110)*V(40)*V(51)
  IRR(111) = RCT(111)*V(37)*V(51)
  IRR(112) = RCT(112)*V(49)*V(51)
  IRR(113) = RCT(113)*V(50)*V(51)
  IRR(114) = RCT(114)*V(51)*V(55)
  IRR(115) = RCT(115)*V(35)*V(51)
  IRR(116) = RCT(116)*V(36)*V(51)
  IRR(117) = RCT(117)*V(45)*V(51)
  IRR(118) = RCT(118)*V(38)*V(51)
  IRR(119) = RCT(119)*V(52)*V(52)
  IRR(120) = RCT(120)*V(52)*V(57)
  IRR(121) = RCT(121)*V(52)*V(56)
  IRR(122) = RCT(122)*V(41)*V(52)
  IRR(123) = RCT(123)*V(44)*V(52)
  IRR(124) = RCT(124)*V(40)*V(52)
  IRR(125) = RCT(125)*V(37)*V(52)
  IRR(126) = RCT(126)*V(49)*V(52)
  IRR(127) = RCT(127)*V(50)*V(52)
  IRR(128) = RCT(128)*V(52)*V(55)
  IRR(129) = RCT(129)*V(35)*V(52)
  IRR(130) = RCT(130)*V(36)*V(52)
  IRR(131) = RCT(131)*V(45)*V(52)
  IRR(132) = RCT(132)*V(55)*V(57)
  IRR(133) = RCT(133)*V(55)*V(56)
  IRR(134) = RCT(134)*V(41)*V(55)
  IRR(135) = RCT(135)*V(44)*V(55)
  IRR(136) = RCT(136)*V(40)*V(55)
  IRR(137) = RCT(137)*V(37)*V(55)
  IRR(138) = RCT(138)*V(49)*V(55)
  IRR(139) = RCT(139)*V(50)*V(55)
  IRR(140) = RCT(140)*V(55)*V(55)
  IRR(141) = RCT(141)*V(35)*V(55)
  IRR(142) = RCT(142)*V(36)*V(55)
  IRR(143) = RCT(143)*V(45)*V(55)
  IRR(144) = RCT(144)*V(39)*V(51)
  IRR(145) = RCT(145)*V(39)*V(52)
  IRR(146) = RCT(146)*V(39)*V(55)
  IRR(147) = RCT(147)*V(39)*V(39)
  IRR(148) = RCT(148)*V(39)*V(58)
  IRR(149) = RCT(149)*V(32)*V(59)
  IRR(150) = RCT(150)*V(32)*V(51)
  IRR(151) = RCT(151)*V(32)*V(52)
  IRR(152) = RCT(152)*V(32)*V(55)
  IRR(153) = RCT(153)*V(32)*V(32)
  IRR(154) = RCT(154)*V(38)*V(52)
  IRR(155) = RCT(155)*V(38)*V(55)
  IRR(156) = RCT(156)*V(38)*V(38)
      
END SUBROUTINE radm2sorg_IRRFun
















SUBROUTINE radm2sorg_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(280)


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

  B(22) = RCT(22)*F(2)

  B(24) = RCT(23)*V(59)

  B(25) = RCT(23)*V(19)

  B(26) = RCT(24)*F(2)

  B(28) = RCT(25)*F(1)

  B(30) = RCT(26)*V(58)

  B(31) = RCT(26)*V(46)

  B(32) = RCT(27)*V(53)

  B(33) = RCT(27)*V(46)

  B(34) = RCT(28)*V(51)

  B(35) = RCT(28)*V(46)

  B(36) = RCT(29)*V(58)

  B(37) = RCT(29)*V(51)

  B(38) = RCT(30)*V(59)

  B(39) = RCT(30)*V(51)

  B(40) = RCT(31)

  B(41) = RCT(32)*2*V(51)

  B(42) = RCT(33)*2*V(51)*F(1)

  B(44) = RCT(34)*V(53)

  B(45) = RCT(34)*V(14)

  B(46) = RCT(35)*V(58)

  B(47) = RCT(35)*V(53)

  B(48) = RCT(36)*2*V(58)*F(2)

  B(50) = RCT(37)*V(59)

  B(51) = RCT(37)*V(46)

  B(52) = RCT(38)*V(58)

  B(53) = RCT(38)*V(54)

  B(54) = RCT(39)*V(59)

  B(55) = RCT(39)*V(54)

  B(56) = RCT(40)*V(54)

  B(57) = RCT(40)*V(51)

  B(58) = RCT(41)*V(59)

  B(59) = RCT(41)*V(54)

  B(60) = RCT(42)

  B(61) = RCT(43)

  B(62) = RCT(44)*V(59)

  B(63) = RCT(44)*V(53)

  B(64) = RCT(45)*V(53)

  B(65) = RCT(45)*V(25)

  B(66) = RCT(46)*V(53)

  B(67) = RCT(46)*V(20)

  B(68) = RCT(47)*V(53)

  B(69) = RCT(47)*V(51)

  B(70) = RCT(48)*V(53)

  B(71) = RCT(48)*V(5)

  B(72) = RCT(49)*V(53)

  B(73) = RCT(49)*V(26)

  B(74) = RCT(50)*V(53)

  B(75) = RCT(50)*V(17)

  B(76) = RCT(51)*V(53)

  B(77) = RCT(51)*V(6)

  B(78) = RCT(52)*V(53)

  B(79) = RCT(52)*V(16)

  B(80) = RCT(53)*V(53)

  B(81) = RCT(53)*V(8)

  B(82) = RCT(54)*V(53)

  B(83) = RCT(54)*V(9)

  B(84) = RCT(55)*V(53)

  B(85) = RCT(55)*V(24)

  B(86) = RCT(56)*V(53)

  B(87) = RCT(56)*V(28)

  B(88) = RCT(57)*V(53)

  B(89) = RCT(57)*V(29)

  B(90) = RCT(58)*V(53)

  B(91) = RCT(58)*V(10)

  B(92) = RCT(59)*V(53)

  B(93) = RCT(59)*V(11)

  B(94) = RCT(60)*V(53)

  B(95) = RCT(60)*V(22)

  B(96) = RCT(61)*V(53)

  B(97) = RCT(61)*V(43)

  B(98) = RCT(62)*V(53)

  B(99) = RCT(62)*V(48)

  B(100) = RCT(63)*V(53)

  B(101) = RCT(63)*V(33)

  B(102) = RCT(64)*V(53)

  B(103) = RCT(64)*V(31)

  B(104) = RCT(65)*V(53)

  B(105) = RCT(65)*V(34)

  B(106) = RCT(66)*V(53)

  B(107) = RCT(66)*V(30)

  B(108) = RCT(67)*V(53)

  B(109) = RCT(67)*V(21)

  B(110) = RCT(68)*V(53)

  B(111) = RCT(68)*V(42)

  B(112) = RCT(69)*V(53)

  B(113) = RCT(69)*V(18)

  B(114) = RCT(70)*V(53)

  B(115) = RCT(70)*V(23)

  B(116) = RCT(71)*V(53)

  B(117) = RCT(71)*V(47)

  B(118) = RCT(72)*V(53)

  B(119) = RCT(72)*V(27)

  B(120) = RCT(73)*V(59)

  B(121) = RCT(73)*V(55)

  B(122) = RCT(74)

  B(123) = RCT(75)*V(59)

  B(124) = RCT(75)*V(45)

  B(125) = RCT(76)

  B(126) = RCT(77)*V(58)

  B(127) = RCT(77)*V(52)

  B(128) = RCT(78)*V(58)

  B(129) = RCT(78)*V(56)

  B(130) = RCT(79)*V(58)

  B(131) = RCT(79)*V(41)

  B(132) = RCT(80)*V(58)

  B(133) = RCT(80)*V(44)

  B(134) = RCT(81)*V(58)

  B(135) = RCT(81)*V(40)

  B(136) = RCT(82)*V(58)

  B(137) = RCT(82)*V(37)

  B(138) = RCT(83)*V(58)

  B(139) = RCT(83)*V(49)

  B(140) = RCT(84)*V(58)

  B(141) = RCT(84)*V(55)

  B(142) = RCT(85)*V(58)

  B(143) = RCT(85)*V(45)

  B(144) = RCT(86)*V(58)

  B(145) = RCT(86)*V(35)

  B(146) = RCT(87)*V(58)

  B(147) = RCT(87)*V(36)

  B(148) = RCT(88)*V(58)

  B(149) = RCT(88)*V(57)

  B(150) = RCT(89)*V(58)

  B(151) = RCT(89)*V(50)

  B(152) = RCT(90)*V(58)

  B(153) = RCT(90)*V(38)

  B(154) = RCT(91)*V(54)

  B(155) = RCT(91)*V(43)

  B(156) = RCT(92)*V(54)

  B(157) = RCT(92)*V(48)

  B(158) = RCT(93)*V(54)

  B(159) = RCT(93)*V(31)

  B(160) = RCT(94)*V(54)

  B(161) = RCT(94)*V(34)

  B(162) = RCT(95)*V(54)

  B(163) = RCT(95)*V(30)

  B(164) = RCT(96)*V(54)

  B(165) = RCT(96)*V(22)

  B(166) = RCT(97)*V(54)

  B(167) = RCT(97)*V(24)

  B(168) = RCT(98)*V(54)

  B(169) = RCT(98)*V(28)

  B(170) = RCT(99)*V(54)

  B(171) = RCT(99)*V(29)

  B(172) = RCT(100)*V(54)

  B(173) = RCT(100)*V(27)

  B(174) = RCT(101)*V(46)

  B(175) = RCT(101)*V(24)

  B(176) = RCT(102)*V(46)

  B(177) = RCT(102)*V(28)

  B(178) = RCT(103)*V(46)

  B(179) = RCT(103)*V(29)

  B(180) = RCT(104)*V(46)

  B(181) = RCT(104)*V(27)

  B(182) = RCT(105)*V(52)

  B(183) = RCT(105)*V(51)

  B(184) = RCT(106)*V(57)

  B(185) = RCT(106)*V(51)

  B(186) = RCT(107)*V(56)

  B(187) = RCT(107)*V(51)

  B(188) = RCT(108)*V(51)

  B(189) = RCT(108)*V(41)

  B(190) = RCT(109)*V(51)

  B(191) = RCT(109)*V(44)

  B(192) = RCT(110)*V(51)

  B(193) = RCT(110)*V(40)

  B(194) = RCT(111)*V(51)

  B(195) = RCT(111)*V(37)

  B(196) = RCT(112)*V(51)

  B(197) = RCT(112)*V(49)

  B(198) = RCT(113)*V(51)

  B(199) = RCT(113)*V(50)

  B(200) = RCT(114)*V(55)

  B(201) = RCT(114)*V(51)

  B(202) = RCT(115)*V(51)

  B(203) = RCT(115)*V(35)

  B(204) = RCT(116)*V(51)

  B(205) = RCT(116)*V(36)

  B(206) = RCT(117)*V(51)

  B(207) = RCT(117)*V(45)

  B(208) = RCT(118)*V(51)

  B(209) = RCT(118)*V(38)

  B(210) = RCT(119)*2*V(52)

  B(211) = RCT(120)*V(57)

  B(212) = RCT(120)*V(52)

  B(213) = RCT(121)*V(56)

  B(214) = RCT(121)*V(52)

  B(215) = RCT(122)*V(52)

  B(216) = RCT(122)*V(41)

  B(217) = RCT(123)*V(52)

  B(218) = RCT(123)*V(44)

  B(219) = RCT(124)*V(52)

  B(220) = RCT(124)*V(40)

  B(221) = RCT(125)*V(52)

  B(222) = RCT(125)*V(37)

  B(223) = RCT(126)*V(52)

  B(224) = RCT(126)*V(49)

  B(225) = RCT(127)*V(52)

  B(226) = RCT(127)*V(50)

  B(227) = RCT(128)*V(55)

  B(228) = RCT(128)*V(52)

  B(229) = RCT(129)*V(52)

  B(230) = RCT(129)*V(35)

  B(231) = RCT(130)*V(52)

  B(232) = RCT(130)*V(36)

  B(233) = RCT(131)*V(52)

  B(234) = RCT(131)*V(45)

  B(235) = RCT(132)*V(57)

  B(236) = RCT(132)*V(55)

  B(237) = RCT(133)*V(56)

  B(238) = RCT(133)*V(55)

  B(239) = RCT(134)*V(55)

  B(240) = RCT(134)*V(41)

  B(241) = RCT(135)*V(55)

  B(242) = RCT(135)*V(44)

  B(243) = RCT(136)*V(55)

  B(244) = RCT(136)*V(40)

  B(245) = RCT(137)*V(55)

  B(246) = RCT(137)*V(37)

  B(247) = RCT(138)*V(55)

  B(248) = RCT(138)*V(49)

  B(249) = RCT(139)*V(55)

  B(250) = RCT(139)*V(50)

  B(251) = RCT(140)*2*V(55)

  B(252) = RCT(141)*V(55)

  B(253) = RCT(141)*V(35)

  B(254) = RCT(142)*V(55)

  B(255) = RCT(142)*V(36)

  B(256) = RCT(143)*V(55)

  B(257) = RCT(143)*V(45)

  B(258) = RCT(144)*V(51)

  B(259) = RCT(144)*V(39)

  B(260) = RCT(145)*V(52)

  B(261) = RCT(145)*V(39)

  B(262) = RCT(146)*V(55)

  B(263) = RCT(146)*V(39)

  B(264) = RCT(147)*2*V(39)

  B(265) = RCT(148)*V(58)

  B(266) = RCT(148)*V(39)

  B(267) = RCT(149)*V(59)

  B(268) = RCT(149)*V(32)

  B(269) = RCT(150)*V(51)

  B(270) = RCT(150)*V(32)

  B(271) = RCT(151)*V(52)

  B(272) = RCT(151)*V(32)

  B(273) = RCT(152)*V(55)

  B(274) = RCT(152)*V(32)

  B(275) = RCT(153)*2*V(32)

  B(276) = RCT(154)*V(52)

  B(277) = RCT(154)*V(38)

  B(278) = RCT(155)*V(55)

  B(279) = RCT(155)*V(38)

  B(280) = RCT(156)*2*V(38)



  JVS(1) = 0

  JVS(2) = B(70)

  JVS(3) = B(71)

  JVS(4) = 0

  JVS(5) = 0.4*B(174)

  JVS(6) = 0.2*B(180)

  JVS(7) = 0.2*B(176)

  JVS(8) = 0.06*B(178)

  JVS(9) = 0.4*B(175)+0.2*B(177)+0.06*B(179)+0.2*B(181)

  JVS(10) = 0

  JVS(11) = 0.2*B(180)

  JVS(12) = 0.2*B(176)

  JVS(13) = 0.29*B(178)

  JVS(14) = 0.5*B(245)

  JVS(15) = 0.5*B(278)

  JVS(16) = 0.5*B(243)

  JVS(17) = 0.5*B(239)

  JVS(18) = 0.5*B(241)

  JVS(19) = 0.5*B(233)

  JVS(20) = 0.2*B(177)+0.29*B(179)+0.2*B(181)

  JVS(21) = 0.5*B(247)

  JVS(22) = 0.5*B(249)

  JVS(23) = 0.5*B(227)+0.5*B(234)

  JVS(24) = 0.5*B(228)+0.5*B(235)+0.5*B(237)+0.5*B(240)+0.5*B(242)+0.5*B(244)+0.5*B(246)+0.5*B(248)+0.5*B(250)+0.5&
              &*B(279)

  JVS(25) = 0.5*B(238)

  JVS(26) = 0.5*B(236)

  JVS(27) = 0

  JVS(28) = B(72)

  JVS(29) = B(73)

  JVS(30) = -B(70)

  JVS(31) = -B(71)

  JVS(32) = -B(76)

  JVS(33) = -B(77)

  JVS(34) = -B(26)-B(28)

  JVS(35) = B(2)

  JVS(36) = -B(80)

  JVS(37) = -B(81)

  JVS(38) = -B(82)

  JVS(39) = -B(83)

  JVS(40) = -B(90)

  JVS(41) = -B(91)

  JVS(42) = -B(92)

  JVS(43) = -B(93)

  JVS(44) = -B(125)

  JVS(45) = B(123)

  JVS(46) = B(124)

  JVS(47) = -B(4)

  JVS(48) = B(46)

  JVS(49) = B(47)

  JVS(50) = -B(9)-B(44)

  JVS(51) = B(41)+B(42)

  JVS(52) = -B(45)

  JVS(53) = -B(60)-B(61)

  JVS(54) = B(58)

  JVS(55) = B(59)

  JVS(56) = -B(78)

  JVS(57) = -B(79)

  JVS(58) = -B(74)

  JVS(59) = 0.06*B(176)

  JVS(60) = 0.09*B(178)

  JVS(61) = 0.06*B(177)+0.09*B(179)

  JVS(62) = -B(75)

  JVS(63) = -B(15)-B(112)

  JVS(64) = B(200)

  JVS(65) = -B(113)

  JVS(66) = B(201)

  JVS(67) = B(26)

  JVS(68) = -B(22)-B(24)

  JVS(69) = B(3)

  JVS(70) = B(8)

  JVS(71) = B(1)-B(25)

  JVS(72) = -B(6)-B(40)-B(66)

  JVS(73) = B(38)

  JVS(74) = -B(67)

  JVS(75) = B(39)

  JVS(76) = -B(13)-B(108)

  JVS(77) = B(182)

  JVS(78) = B(183)

  JVS(79) = -B(109)

  JVS(80) = 0.25*B(90)

  JVS(81) = 0.17*B(92)

  JVS(82) = -B(94)-0.5*B(164)

  JVS(83) = 0.25*B(91)+0.17*B(93)-B(95)

  JVS(84) = -0.5*B(165)

  JVS(85) = -B(114)-B(122)

  JVS(86) = -B(115)

  JVS(87) = B(120)

  JVS(88) = B(121)

  JVS(89) = -B(84)-B(166)-B(174)

  JVS(90) = -B(175)

  JVS(91) = -B(85)

  JVS(92) = -B(167)

  JVS(93) = 2*B(61)

  JVS(94) = B(164)

  JVS(95) = -B(5)-B(64)

  JVS(96) = B(162)

  JVS(97) = B(158)

  JVS(98) = B(160)

  JVS(99) = B(154)

  JVS(100) = B(156)

  JVS(101) = B(56)

  JVS(102) = B(62)-B(65)

  JVS(103) = B(57)+B(155)+B(157)+B(159)+B(161)+B(163)+B(165)

  JVS(104) = B(63)

  JVS(105) = 0.42*B(174)

  JVS(106) = -B(72)

  JVS(107) = 0.33*B(180)

  JVS(108) = 0.33*B(176)

  JVS(109) = 0.23*B(178)

  JVS(110) = 1.87*B(17)+1.55*B(18)+2*B(102)+2*B(158)

  JVS(111) = B(19)+B(104)+B(160)

  JVS(112) = B(10)+B(11)+B(96)+B(154)

  JVS(113) = 0.95*B(142)+0.475*B(233)+0.95*B(256)

  JVS(114) = 0.42*B(175)+0.33*B(177)+0.23*B(179)+0.33*B(181)

  JVS(115) = B(12)

  JVS(116) = 0.475*B(234)

  JVS(117) = -B(73)+B(97)+2*B(103)+B(105)

  JVS(118) = B(155)+2*B(159)+B(161)

  JVS(119) = 0.95*B(257)

  JVS(120) = 0.95*B(143)

  JVS(121) = -B(118)-B(172)-B(180)

  JVS(122) = -B(181)

  JVS(123) = -B(119)

  JVS(124) = -B(173)

  JVS(125) = -B(86)-B(168)-B(176)

  JVS(126) = -B(177)

  JVS(127) = -B(87)

  JVS(128) = -B(169)

  JVS(129) = -B(88)-B(170)-B(178)

  JVS(130) = -B(179)

  JVS(131) = -B(89)

  JVS(132) = -B(171)

  JVS(133) = -B(20)-B(106)-B(162)

  JVS(134) = 0.7*B(144)+0.7*B(229)+B(252)

  JVS(135) = 0.806*B(146)+0.806*B(231)+B(254)

  JVS(136) = 0.7*B(230)+0.806*B(232)

  JVS(137) = -B(107)

  JVS(138) = -B(163)

  JVS(139) = B(253)+B(255)

  JVS(140) = 0.7*B(145)+0.806*B(147)

  JVS(141) = -B(17)-B(18)-B(102)-B(158)

  JVS(142) = 0.16*B(144)+0.16*B(229)+0.2*B(252)

  JVS(143) = 0.89*B(142)+0.445*B(233)+0.89*B(256)

  JVS(144) = 0.16*B(230)+0.445*B(234)

  JVS(145) = -B(103)

  JVS(146) = -B(159)

  JVS(147) = 0.2*B(253)+0.89*B(257)

  JVS(148) = 0.89*B(143)+0.16*B(145)

  JVS(149) = B(164)

  JVS(150) = -B(267)-B(269)-B(271)-B(273)-2*B(275)

  JVS(151) = -B(270)

  JVS(152) = -B(272)

  JVS(153) = 0

  JVS(154) = B(165)

  JVS(155) = -B(274)

  JVS(156) = -B(268)

  JVS(157) = 0.025*B(78)

  JVS(158) = 0.1*B(178)

  JVS(159) = -B(16)-B(100)

  JVS(160) = 0.69*B(130)+0.75*B(215)+0.86*B(239)

  JVS(161) = 1.06*B(132)+1.39*B(217)+0.9*B(241)

  JVS(162) = 0.1*B(179)

  JVS(163) = 0.8*B(21)

  JVS(164) = 0.1*B(138)+0.55*B(223)+0.55*B(247)

  JVS(165) = 0.6*B(213)+0.75*B(216)+1.39*B(218)+0.55*B(224)

  JVS(166) = 0.025*B(79)-B(101)

  JVS(167) = 0

  JVS(168) = 0.8*B(237)+0.86*B(240)+0.9*B(242)+0.55*B(248)

  JVS(169) = 0.25*B(128)+0.6*B(214)+0.8*B(238)

  JVS(170) = 0.25*B(129)+0.69*B(131)+1.06*B(133)+0.1*B(139)

  JVS(171) = -B(19)-B(104)-B(160)

  JVS(172) = 0.17*B(144)+0.17*B(229)+0.8*B(252)

  JVS(173) = 0.45*B(146)+0.45*B(231)+B(254)

  JVS(174) = 0.11*B(142)+0.055*B(233)+0.11*B(256)

  JVS(175) = B(150)+0.75*B(225)+B(249)

  JVS(176) = 0.75*B(226)+0.17*B(230)+0.45*B(232)+0.055*B(234)

  JVS(177) = -B(105)

  JVS(178) = -B(161)

  JVS(179) = B(250)+0.8*B(253)+B(255)+0.11*B(257)

  JVS(180) = 0.11*B(143)+0.17*B(145)+0.45*B(147)+B(151)

  JVS(181) = 0.75*B(90)

  JVS(182) = -B(144)-B(202)-B(229)-B(252)

  JVS(183) = -B(203)

  JVS(184) = -B(230)

  JVS(185) = 0.75*B(91)

  JVS(186) = -B(253)

  JVS(187) = -B(145)

  JVS(188) = 0.83*B(92)

  JVS(189) = -B(146)-B(204)-B(231)-B(254)

  JVS(190) = -B(205)

  JVS(191) = -B(232)

  JVS(192) = 0.83*B(93)

  JVS(193) = -B(255)

  JVS(194) = -B(147)

  JVS(195) = B(118)

  JVS(196) = B(86)

  JVS(197) = -B(136)-B(194)-B(221)-B(245)

  JVS(198) = 0

  JVS(199) = -B(195)

  JVS(200) = -B(222)

  JVS(201) = B(87)+B(119)

  JVS(202) = 0

  JVS(203) = -B(246)

  JVS(204) = -B(137)

  JVS(205) = B(166)

  JVS(206) = B(172)

  JVS(207) = B(168)

  JVS(208) = B(170)

  JVS(209) = -B(152)-B(208)-B(276)-B(278)-2*B(280)

  JVS(210) = 0

  JVS(211) = -B(209)

  JVS(212) = -B(277)

  JVS(213) = 0

  JVS(214) = B(167)+B(169)+B(171)+B(173)

  JVS(215) = -B(279)

  JVS(216) = -B(153)

  JVS(217) = 0.25*B(80)

  JVS(218) = 0.75*B(82)

  JVS(219) = 0.9*B(94)

  JVS(220) = B(114)

  JVS(221) = -B(258)-B(260)-B(262)-2*B(264)-B(265)

  JVS(222) = 2*B(142)+B(233)+2*B(256)

  JVS(223) = -B(259)

  JVS(224) = B(234)-B(261)

  JVS(225) = 0.25*B(81)+0.75*B(83)+0.9*B(95)+B(115)

  JVS(226) = 0

  JVS(227) = 2*B(257)-B(263)

  JVS(228) = 2*B(143)-B(266)

  JVS(229) = 0

  JVS(230) = B(84)

  JVS(231) = -B(134)-B(192)-B(219)-B(243)

  JVS(232) = 0

  JVS(233) = -B(193)

  JVS(234) = -B(220)

  JVS(235) = B(85)

  JVS(236) = 0

  JVS(237) = -B(244)

  JVS(238) = -B(135)

  JVS(239) = B(80)

  JVS(240) = -B(130)-B(188)-B(215)-B(239)

  JVS(241) = -B(189)

  JVS(242) = -B(216)

  JVS(243) = B(81)

  JVS(244) = -B(240)

  JVS(245) = -B(131)

  JVS(246) = B(269)

  JVS(247) = B(202)

  JVS(248) = B(204)

  JVS(249) = B(194)

  JVS(250) = B(258)

  JVS(251) = B(192)

  JVS(252) = B(188)

  JVS(253) = -B(14)-B(110)

  JVS(254) = B(190)

  JVS(255) = B(206)

  JVS(256) = 0

  JVS(257) = B(196)

  JVS(258) = B(198)

  JVS(259) = B(184)+B(186)+B(189)+B(191)+B(193)+B(195)+B(197)+B(199)+B(203)+B(205)+B(207)+B(259)+B(270)

  JVS(260) = 0

  JVS(261) = -B(111)

  JVS(262) = 0

  JVS(263) = 0

  JVS(264) = B(187)

  JVS(265) = B(185)

  JVS(266) = 0

  JVS(267) = 0

  JVS(268) = 0.009*B(78)

  JVS(269) = B(13)+0.5*B(108)

  JVS(270) = B(114)

  JVS(271) = B(174)

  JVS(272) = 0.53*B(180)

  JVS(273) = 0.53*B(176)

  JVS(274) = 0.18*B(178)

  JVS(275) = 0.13*B(17)+0.45*B(18)

  JVS(276) = B(271)

  JVS(277) = B(229)

  JVS(278) = B(231)

  JVS(279) = B(136)+1.25*B(221)+0.5*B(245)

  JVS(280) = B(152)+1.75*B(276)+B(278)+2*B(280)

  JVS(281) = B(260)

  JVS(282) = 1.6*B(134)+1.55*B(219)+0.8*B(243)

  JVS(283) = 0.77*B(215)

  JVS(284) = -B(10)-B(11)-B(96)-B(154)

  JVS(285) = 0.04*B(132)+0.8*B(217)

  JVS(286) = 0.5*B(233)

  JVS(287) = B(175)+0.53*B(177)+0.18*B(179)+0.53*B(181)

  JVS(288) = 0.28*B(138)+0.89*B(223)+0.14*B(247)

  JVS(289) = 0.75*B(225)

  JVS(290) = 0

  JVS(291) = B(126)+1.5*B(210)+0.75*B(211)+0.75*B(213)+0.77*B(216)+0.8*B(218)+1.55*B(220)+1.25*B(222)+0.89*B(224)+0.75&
               &*B(226)+B(227)+B(230)+B(232)+0.5*B(234)+B(261)+B(272)+1.75*B(277)

  JVS(292) = 0.009*B(79)-B(97)+0.5*B(109)+B(115)

  JVS(293) = -B(155)

  JVS(294) = B(228)+0.8*B(244)+0.5*B(246)+0.14*B(248)+B(279)

  JVS(295) = 0.09*B(128)+0.75*B(214)

  JVS(296) = 0.75*B(212)

  JVS(297) = B(127)+0.09*B(129)+0.04*B(133)+1.6*B(135)+B(137)+0.28*B(139)+B(153)

  JVS(298) = 0

  JVS(299) = B(82)

  JVS(300) = -B(132)-B(190)-B(217)-B(241)

  JVS(301) = -B(191)

  JVS(302) = -B(218)

  JVS(303) = B(83)

  JVS(304) = -B(242)

  JVS(305) = -B(133)

  JVS(306) = B(125)

  JVS(307) = 0.9*B(94)

  JVS(308) = B(20)+B(106)+B(162)

  JVS(309) = 0

  JVS(310) = 0

  JVS(311) = -B(123)-B(142)-B(206)-B(233)-B(256)

  JVS(312) = -B(207)

  JVS(313) = -B(234)

  JVS(314) = 0.9*B(95)+B(107)

  JVS(315) = B(163)

  JVS(316) = -B(257)

  JVS(317) = -B(143)

  JVS(318) = -B(124)

  JVS(319) = B(22)

  JVS(320) = -B(174)

  JVS(321) = -B(180)

  JVS(322) = -B(176)

  JVS(323) = -B(178)

  JVS(324) = -B(2)-B(3)-B(30)-B(32)-B(34)-B(50)-B(175)-B(177)-B(179)-B(181)

  JVS(325) = -B(35)

  JVS(326) = -B(33)

  JVS(327) = 0

  JVS(328) = -B(31)

  JVS(329) = -B(51)

  JVS(330) = B(267)

  JVS(331) = B(208)

  JVS(332) = 0.08*B(130)

  JVS(333) = 0.24*B(132)

  JVS(334) = 0

  JVS(335) = -B(21)-B(116)

  JVS(336) = B(209)

  JVS(337) = 0

  JVS(338) = -B(117)

  JVS(339) = 0

  JVS(340) = 0

  JVS(341) = 0.036*B(128)

  JVS(342) = 0.036*B(129)+0.08*B(131)+0.24*B(133)

  JVS(343) = B(268)

  JVS(344) = 0.075*B(78)

  JVS(345) = 0.5*B(180)

  JVS(346) = 0.5*B(176)

  JVS(347) = 0.72*B(178)

  JVS(348) = B(136)+0.75*B(221)+B(245)

  JVS(349) = B(152)+B(276)+B(278)+2*B(280)

  JVS(350) = 0.2*B(134)+0.35*B(219)+0.6*B(243)

  JVS(351) = 0.38*B(130)+0.41*B(215)+0.14*B(239)

  JVS(352) = B(14)+0.5*B(110)

  JVS(353) = 0.35*B(132)+0.46*B(217)+0.1*B(241)

  JVS(354) = 0

  JVS(355) = 0.5*B(177)+0.72*B(179)+0.5*B(181)

  JVS(356) = 0.2*B(21)

  JVS(357) = -B(12)-B(98)-B(156)

  JVS(358) = 1.45*B(138)+0.725*B(223)+0.725*B(247)

  JVS(359) = 0

  JVS(360) = 0

  JVS(361) = 0.75*B(211)+0.15*B(213)+0.41*B(216)+0.46*B(218)+0.35*B(220)+0.75*B(222)+0.725*B(224)+B(277)

  JVS(362) = 0.075*B(79)-B(99)+0.5*B(111)

  JVS(363) = -B(157)

  JVS(364) = B(235)+0.2*B(237)+0.14*B(240)+0.1*B(242)+0.6*B(244)+B(246)+0.725*B(248)+B(279)

  JVS(365) = 0.75*B(128)+0.15*B(214)+0.2*B(238)

  JVS(366) = B(148)+0.75*B(212)+B(236)

  JVS(367) = 0.75*B(129)+0.38*B(131)+0.35*B(133)+0.2*B(135)+B(137)+1.45*B(139)+B(149)+B(153)

  JVS(368) = 0

  JVS(369) = B(88)

  JVS(370) = 0

  JVS(371) = -B(138)-B(196)-B(223)-B(247)

  JVS(372) = -B(197)

  JVS(373) = -B(224)

  JVS(374) = B(89)

  JVS(375) = 0

  JVS(376) = -B(248)

  JVS(377) = -B(139)

  JVS(378) = 0

  JVS(379) = B(100)

  JVS(380) = 0

  JVS(381) = 0

  JVS(382) = 0

  JVS(383) = 0

  JVS(384) = 0

  JVS(385) = -B(150)-B(198)-B(225)-B(249)

  JVS(386) = -B(199)

  JVS(387) = -B(226)

  JVS(388) = B(101)

  JVS(389) = 0

  JVS(390) = -B(250)

  JVS(391) = 0

  JVS(392) = -B(151)

  JVS(393) = 0

  JVS(394) = B(70)

  JVS(395) = 0.25*B(90)

  JVS(396) = 0.17*B(92)

  JVS(397) = B(44)

  JVS(398) = 0.17*B(78)

  JVS(399) = 0.65*B(6)+B(40)

  JVS(400) = B(13)

  JVS(401) = 0.1*B(94)

  JVS(402) = 0.12*B(174)

  JVS(403) = B(72)

  JVS(404) = 0.23*B(180)

  JVS(405) = 0.23*B(176)

  JVS(406) = 0.26*B(178)

  JVS(407) = B(20)

  JVS(408) = 0.8*B(18)+B(102)+B(158)

  JVS(409) = -B(269)+B(271)

  JVS(410) = B(19)

  JVS(411) = B(144)-B(202)+2*B(229)+B(252)

  JVS(412) = B(146)-B(204)+2*B(231)+B(254)

  JVS(413) = B(136)-B(194)+B(221)+0.5*B(245)

  JVS(414) = -B(208)+0.5*B(276)

  JVS(415) = -B(258)+B(260)

  JVS(416) = B(134)-B(192)+B(219)+0.5*B(243)

  JVS(417) = 0.92*B(130)-B(188)+B(215)+0.5*B(239)

  JVS(418) = B(14)

  JVS(419) = 2*B(11)+B(96)+B(154)

  JVS(420) = 0.76*B(132)-B(190)+B(217)+0.5*B(241)

  JVS(421) = 0.92*B(142)-B(206)+0.46*B(233)+0.92*B(256)

  JVS(422) = B(32)-B(34)+0.12*B(175)+0.23*B(177)+0.26*B(179)+0.23*B(181)

  JVS(423) = B(21)

  JVS(424) = B(12)

  JVS(425) = B(138)-B(196)+B(223)+0.5*B(247)

  JVS(426) = B(150)-B(198)+B(225)+0.5*B(249)

  JVS(427) = -B(35)-B(36)-B(38)-2*B(41)-2*B(42)-B(56)-B(68)-B(182)-B(184)-B(186)-B(189)-B(191)-B(193)-B(195)-B(197)&
               &-B(199)-B(200)-B(203)-B(205)-B(207)-B(209)-B(259)-B(270)

  JVS(428) = B(126)-B(183)+B(210)+B(211)+B(213)+B(216)+B(218)+B(220)+B(222)+B(224)+B(226)+0.5*B(227)+2*B(230)+2*B(232)&
               &+0.46*B(234)+B(261)+B(272)+0.5*B(277)

  JVS(429) = B(33)+B(45)-B(69)+B(71)+B(73)+0.17*B(79)+0.25*B(91)+0.17*B(93)+0.1*B(95)+B(97)+B(103)

  JVS(430) = -B(57)+B(155)+B(159)

  JVS(431) = -B(201)+0.5*B(228)+0.5*B(235)+0.5*B(237)+0.5*B(240)+0.5*B(242)+0.5*B(244)+0.5*B(246)+0.5*B(248)+0.5*B(250)&
               &+B(253)+B(255)+0.92*B(257)

  JVS(432) = 0.964*B(128)-B(187)+B(214)+0.5*B(238)

  JVS(433) = B(148)-B(185)+B(212)+0.5*B(236)

  JVS(434) = -B(37)+B(127)+0.964*B(129)+0.92*B(131)+0.76*B(133)+B(135)+B(137)+B(139)+0.92*B(143)+B(145)+B(147)+B(149)&
               &+B(151)

  JVS(435) = -B(39)

  JVS(436) = B(74)

  JVS(437) = B(15)

  JVS(438) = 0.5*B(108)

  JVS(439) = 0.22*B(180)

  JVS(440) = 0.22*B(176)

  JVS(441) = 0.31*B(178)

  JVS(442) = -B(271)+B(273)

  JVS(443) = -B(229)+B(252)

  JVS(444) = -B(231)+B(254)

  JVS(445) = -B(221)+0.5*B(245)

  JVS(446) = -B(276)+0.5*B(278)

  JVS(447) = -B(260)+B(262)

  JVS(448) = -B(219)+0.5*B(243)

  JVS(449) = -B(215)+0.5*B(239)

  JVS(450) = -B(217)+0.5*B(241)

  JVS(451) = -B(233)+B(256)

  JVS(452) = 0.22*B(177)+0.31*B(179)+0.22*B(181)

  JVS(453) = B(12)

  JVS(454) = -B(223)+0.5*B(247)

  JVS(455) = -B(225)+0.5*B(249)

  JVS(456) = -B(182)

  JVS(457) = -B(126)-B(183)-2*B(210)-B(211)-B(213)-B(216)-B(218)-B(220)-B(222)-B(224)-B(226)-0.5*B(227)-B(230)-B(232)&
               &-B(234)-B(261)-B(272)-B(277)

  JVS(458) = B(75)+0.5*B(109)

  JVS(459) = 0

  JVS(460) = B(140)-0.5*B(228)+0.5*B(235)+0.5*B(237)+0.5*B(240)+0.5*B(242)+0.5*B(244)+0.5*B(246)+0.5*B(248)+0.5*B(250)+2&
               &*B(251)+B(253)+B(255)+B(257)+B(263)+B(274)+0.5*B(279)

  JVS(461) = -B(214)+0.5*B(238)

  JVS(462) = -B(212)+0.5*B(236)

  JVS(463) = -B(127)+B(141)

  JVS(464) = 0

  JVS(465) = -B(70)

  JVS(466) = -B(76)

  JVS(467) = 2*B(28)

  JVS(468) = -B(80)

  JVS(469) = -B(82)

  JVS(470) = -B(90)

  JVS(471) = -B(92)

  JVS(472) = B(4)

  JVS(473) = 2*B(9)-B(44)

  JVS(474) = -B(78)

  JVS(475) = -B(74)

  JVS(476) = B(15)-B(112)

  JVS(477) = 0.35*B(6)-B(66)

  JVS(478) = B(13)-0.5*B(108)

  JVS(479) = -1.9*B(94)

  JVS(480) = -B(114)

  JVS(481) = -B(84)

  JVS(482) = B(5)-B(64)

  JVS(483) = -B(72)

  JVS(484) = -B(118)+0.1*B(180)

  JVS(485) = -B(86)+0.1*B(176)

  JVS(486) = -B(88)+0.14*B(178)

  JVS(487) = -B(106)

  JVS(488) = -B(102)

  JVS(489) = -B(100)

  JVS(490) = -B(104)

  JVS(491) = 0

  JVS(492) = 0

  JVS(493) = 0

  JVS(494) = B(14)-0.5*B(110)

  JVS(495) = -B(96)

  JVS(496) = 0

  JVS(497) = 0

  JVS(498) = -B(32)+B(34)+0.1*B(177)+0.14*B(179)+0.1*B(181)

  JVS(499) = -B(116)

  JVS(500) = -B(98)

  JVS(501) = 0

  JVS(502) = 0

  JVS(503) = B(35)+B(36)-B(68)

  JVS(504) = 0

  JVS(505) = -B(33)-B(45)-B(46)-B(62)-B(65)-B(67)-B(69)-B(71)-B(73)-B(75)-B(77)-B(79)-B(81)-B(83)-B(85)-B(87)-B(89)&
               &-B(91)-B(93)-1.9*B(95)-B(97)-B(99)-B(101)-B(103)-B(105)-B(107)-0.5*B(109)-0.5*B(111)-B(113)-B(115)-B(117)&
               &-B(119)

  JVS(506) = 0

  JVS(507) = 0

  JVS(508) = 0

  JVS(509) = 0

  JVS(510) = B(37)-B(47)

  JVS(511) = -B(63)

  JVS(512) = B(60)

  JVS(513) = 0.35*B(6)

  JVS(514) = -B(164)

  JVS(515) = B(114)

  JVS(516) = -B(166)

  JVS(517) = B(64)

  JVS(518) = -B(172)

  JVS(519) = -B(168)

  JVS(520) = -B(170)

  JVS(521) = -B(162)

  JVS(522) = -B(158)

  JVS(523) = -B(160)

  JVS(524) = 0

  JVS(525) = 0

  JVS(526) = -B(154)

  JVS(527) = 0

  JVS(528) = 0

  JVS(529) = B(50)

  JVS(530) = -B(156)

  JVS(531) = 0

  JVS(532) = 0

  JVS(533) = -B(56)

  JVS(534) = 0

  JVS(535) = B(65)+B(115)

  JVS(536) = -B(7)-B(8)-B(52)-B(54)-B(57)-B(58)-B(155)-B(157)-B(159)-B(161)-B(163)-B(165)-B(167)-B(169)-B(171)-B(173)

  JVS(537) = 0

  JVS(538) = 0

  JVS(539) = 0

  JVS(540) = -B(53)

  JVS(541) = B(51)-B(55)-B(59)

  JVS(542) = B(112)

  JVS(543) = B(122)

  JVS(544) = -B(273)

  JVS(545) = B(16)

  JVS(546) = B(19)+B(104)+B(160)

  JVS(547) = -B(252)

  JVS(548) = -B(254)

  JVS(549) = -B(245)

  JVS(550) = -B(278)

  JVS(551) = -B(262)

  JVS(552) = -B(243)

  JVS(553) = -B(239)

  JVS(554) = -B(241)

  JVS(555) = 0.05*B(142)+0.025*B(233)-0.95*B(256)

  JVS(556) = 0

  JVS(557) = 0

  JVS(558) = B(98)+B(156)

  JVS(559) = -B(247)

  JVS(560) = -B(249)

  JVS(561) = -B(200)

  JVS(562) = -B(227)+0.025*B(234)

  JVS(563) = B(99)+B(105)+B(113)

  JVS(564) = B(157)+B(161)

  JVS(565) = -B(120)-B(140)-B(201)-B(228)-B(235)-B(237)-B(240)-B(242)-B(244)-B(246)-B(248)-B(250)-2*B(251)-B(253)-B(255)&
               &-0.95*B(257)-B(263)-B(274)-B(279)

  JVS(566) = -B(238)

  JVS(567) = -B(236)

  JVS(568) = -B(141)+0.05*B(143)

  JVS(569) = -B(121)

  JVS(570) = 0.83*B(78)

  JVS(571) = 0.5*B(110)

  JVS(572) = 0

  JVS(573) = 0

  JVS(574) = 0

  JVS(575) = B(116)

  JVS(576) = 0

  JVS(577) = 0

  JVS(578) = -B(186)

  JVS(579) = -B(213)

  JVS(580) = 0.83*B(79)+0.5*B(111)+B(117)

  JVS(581) = 0

  JVS(582) = -B(237)

  JVS(583) = -B(128)-B(187)-B(214)-B(238)

  JVS(584) = 0

  JVS(585) = -B(129)

  JVS(586) = 0

  JVS(587) = B(76)

  JVS(588) = B(16)

  JVS(589) = 0

  JVS(590) = 0

  JVS(591) = 0

  JVS(592) = 0

  JVS(593) = 0

  JVS(594) = -B(184)

  JVS(595) = -B(211)

  JVS(596) = B(77)

  JVS(597) = 0

  JVS(598) = -B(235)

  JVS(599) = 0

  JVS(600) = -B(148)-B(185)-B(212)-B(236)

  JVS(601) = -B(149)

  JVS(602) = 0

  JVS(603) = B(4)

  JVS(604) = B(24)

  JVS(605) = -B(144)

  JVS(606) = -B(146)

  JVS(607) = -B(136)

  JVS(608) = -B(152)

  JVS(609) = -B(265)

  JVS(610) = -B(134)

  JVS(611) = -B(130)

  JVS(612) = -B(132)

  JVS(613) = -B(142)

  JVS(614) = -B(30)

  JVS(615) = -B(138)

  JVS(616) = -B(150)

  JVS(617) = -B(36)

  JVS(618) = -B(126)

  JVS(619) = -B(46)

  JVS(620) = B(7)-B(52)+B(54)

  JVS(621) = -B(140)

  JVS(622) = -B(128)

  JVS(623) = -B(148)

  JVS(624) = -B(31)-B(37)-B(47)-2*B(48)-B(53)-B(127)-B(129)-B(131)-B(133)-B(135)-B(137)-B(139)-B(141)-B(143)-B(145)&
               &-B(147)-B(149)-B(151)-B(153)-B(266)

  JVS(625) = B(1)+B(25)+B(55)

  JVS(626) = B(125)

  JVS(627) = B(60)

  JVS(628) = -B(24)

  JVS(629) = 0.65*B(6)+B(40)+B(66)

  JVS(630) = B(122)

  JVS(631) = B(5)

  JVS(632) = 0

  JVS(633) = 0

  JVS(634) = -B(267)

  JVS(635) = 0

  JVS(636) = B(144)

  JVS(637) = B(146)

  JVS(638) = B(136)

  JVS(639) = 2*B(152)+B(276)+B(278)+2*B(280)

  JVS(640) = B(265)

  JVS(641) = B(134)

  JVS(642) = 0.92*B(130)

  JVS(643) = 0

  JVS(644) = 0.76*B(132)

  JVS(645) = -B(123)+B(142)

  JVS(646) = B(30)-B(50)

  JVS(647) = B(21)+B(116)

  JVS(648) = 0

  JVS(649) = B(138)

  JVS(650) = B(150)

  JVS(651) = B(36)-B(38)

  JVS(652) = B(126)+B(277)

  JVS(653) = -B(62)+B(67)+B(117)

  JVS(654) = B(8)+2*B(52)-B(58)

  JVS(655) = -B(120)+B(140)+B(279)

  JVS(656) = 0.964*B(128)

  JVS(657) = B(148)

  JVS(658) = B(31)+B(37)+2*B(48)+2*B(53)+B(127)+0.964*B(129)+0.92*B(131)+0.76*B(133)+B(135)+B(137)+B(139)+B(141)+B(143)&
               &+B(145)+B(147)+B(149)+B(151)+2*B(153)+B(266)

  JVS(659) = -B(1)-B(25)-B(39)-B(51)-B(59)-B(63)-B(121)-B(124)-B(268)
      
END SUBROUTINE radm2sorg_Jac_SP














SUBROUTINE radm2sorg_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(659), W(59), a
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
      
END SUBROUTINE radm2sorg_KppDecomp



SUBROUTINE radm2sorg_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(659), W(59), a
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
      
END SUBROUTINE radm2sorg_KppDecompCmplx


SUBROUTINE radm2sorg_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(659), X(59), sum

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
      
END SUBROUTINE radm2sorg_KppSolveIndirect


SUBROUTINE radm2sorg_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(659), X(59), sum

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
      
END SUBROUTINE radm2sorg_KppSolveCmplx













SUBROUTINE radm2sorg_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(19) = X(19)-JVS(67)*X(7)
  X(22) = X(22)-JVS(80)*X(10)-JVS(81)*X(11)
  X(25) = X(25)-JVS(93)*X(15)-JVS(94)*X(22)
  X(26) = X(26)-JVS(105)*X(24)
  X(32) = X(32)-JVS(149)*X(22)
  X(33) = X(33)-JVS(157)*X(16)-JVS(158)*X(29)
  X(35) = X(35)-JVS(181)*X(10)
  X(36) = X(36)-JVS(188)*X(11)
  X(37) = X(37)-JVS(195)*X(27)-JVS(196)*X(28)
  X(38) = X(38)-JVS(205)*X(24)-JVS(206)*X(27)-JVS(207)*X(28)-JVS(208)*X(29)
  X(39) = X(39)-JVS(217)*X(8)-JVS(218)*X(9)-JVS(219)*X(22)-JVS(220)*X(23)
  X(40) = X(40)-JVS(230)*X(24)
  X(41) = X(41)-JVS(239)*X(8)
  X(42) = X(42)-JVS(246)*X(32)-JVS(247)*X(35)-JVS(248)*X(36)-JVS(249)*X(37)-JVS(250)*X(39)-JVS(251)*X(40)-JVS(252)*X(41)
  X(43) = X(43)-JVS(268)*X(16)-JVS(269)*X(21)-JVS(270)*X(23)-JVS(271)*X(24)-JVS(272)*X(27)-JVS(273)*X(28)-JVS(274)*X(29)&
            &-JVS(275)*X(31)-JVS(276)*X(32)-JVS(277)*X(35)-JVS(278)*X(36)-JVS(279)*X(37)-JVS(280)*X(38)-JVS(281)*X(39)&
            &-JVS(282)*X(40)-JVS(283)*X(41)
  X(44) = X(44)-JVS(299)*X(9)
  X(45) = X(45)-JVS(306)*X(12)-JVS(307)*X(22)-JVS(308)*X(30)-JVS(309)*X(35)-JVS(310)*X(36)
  X(46) = X(46)-JVS(319)*X(19)-JVS(320)*X(24)-JVS(321)*X(27)-JVS(322)*X(28)-JVS(323)*X(29)
  X(47) = X(47)-JVS(330)*X(32)-JVS(331)*X(38)-JVS(332)*X(41)-JVS(333)*X(44)-JVS(334)*X(46)
  X(48) = X(48)-JVS(344)*X(16)-JVS(345)*X(27)-JVS(346)*X(28)-JVS(347)*X(29)-JVS(348)*X(37)-JVS(349)*X(38)-JVS(350)*X(40)&
            &-JVS(351)*X(41)-JVS(352)*X(42)-JVS(353)*X(44)-JVS(354)*X(45)-JVS(355)*X(46)-JVS(356)*X(47)
  X(49) = X(49)-JVS(369)*X(29)-JVS(370)*X(46)
  X(50) = X(50)-JVS(379)*X(33)-JVS(380)*X(41)-JVS(381)*X(44)-JVS(382)*X(46)-JVS(383)*X(47)-JVS(384)*X(49)
  X(51) = X(51)-JVS(394)*X(5)-JVS(395)*X(10)-JVS(396)*X(11)-JVS(397)*X(14)-JVS(398)*X(16)-JVS(399)*X(20)-JVS(400)*X(21)&
            &-JVS(401)*X(22)-JVS(402)*X(24)-JVS(403)*X(26)-JVS(404)*X(27)-JVS(405)*X(28)-JVS(406)*X(29)-JVS(407)*X(30)&
            &-JVS(408)*X(31)-JVS(409)*X(32)-JVS(410)*X(34)-JVS(411)*X(35)-JVS(412)*X(36)-JVS(413)*X(37)-JVS(414)*X(38)&
            &-JVS(415)*X(39)-JVS(416)*X(40)-JVS(417)*X(41)-JVS(418)*X(42)-JVS(419)*X(43)-JVS(420)*X(44)-JVS(421)*X(45)&
            &-JVS(422)*X(46)-JVS(423)*X(47)-JVS(424)*X(48)-JVS(425)*X(49)-JVS(426)*X(50)
  X(52) = X(52)-JVS(436)*X(17)-JVS(437)*X(18)-JVS(438)*X(21)-JVS(439)*X(27)-JVS(440)*X(28)-JVS(441)*X(29)-JVS(442)*X(32)&
            &-JVS(443)*X(35)-JVS(444)*X(36)-JVS(445)*X(37)-JVS(446)*X(38)-JVS(447)*X(39)-JVS(448)*X(40)-JVS(449)*X(41)&
            &-JVS(450)*X(44)-JVS(451)*X(45)-JVS(452)*X(46)-JVS(453)*X(48)-JVS(454)*X(49)-JVS(455)*X(50)-JVS(456)*X(51)
  X(53) = X(53)-JVS(465)*X(5)-JVS(466)*X(6)-JVS(467)*X(7)-JVS(468)*X(8)-JVS(469)*X(9)-JVS(470)*X(10)-JVS(471)*X(11)&
            &-JVS(472)*X(13)-JVS(473)*X(14)-JVS(474)*X(16)-JVS(475)*X(17)-JVS(476)*X(18)-JVS(477)*X(20)-JVS(478)*X(21)&
            &-JVS(479)*X(22)-JVS(480)*X(23)-JVS(481)*X(24)-JVS(482)*X(25)-JVS(483)*X(26)-JVS(484)*X(27)-JVS(485)*X(28)&
            &-JVS(486)*X(29)-JVS(487)*X(30)-JVS(488)*X(31)-JVS(489)*X(33)-JVS(490)*X(34)-JVS(491)*X(35)-JVS(492)*X(36)&
            &-JVS(493)*X(41)-JVS(494)*X(42)-JVS(495)*X(43)-JVS(496)*X(44)-JVS(497)*X(45)-JVS(498)*X(46)-JVS(499)*X(47)&
            &-JVS(500)*X(48)-JVS(501)*X(49)-JVS(502)*X(50)-JVS(503)*X(51)-JVS(504)*X(52)
  X(54) = X(54)-JVS(512)*X(15)-JVS(513)*X(20)-JVS(514)*X(22)-JVS(515)*X(23)-JVS(516)*X(24)-JVS(517)*X(25)-JVS(518)*X(27)&
            &-JVS(519)*X(28)-JVS(520)*X(29)-JVS(521)*X(30)-JVS(522)*X(31)-JVS(523)*X(34)-JVS(524)*X(35)-JVS(525)*X(36)&
            &-JVS(526)*X(43)-JVS(527)*X(44)-JVS(528)*X(45)-JVS(529)*X(46)-JVS(530)*X(48)-JVS(531)*X(49)-JVS(532)*X(50)&
            &-JVS(533)*X(51)-JVS(534)*X(52)-JVS(535)*X(53)
  X(55) = X(55)-JVS(542)*X(18)-JVS(543)*X(23)-JVS(544)*X(32)-JVS(545)*X(33)-JVS(546)*X(34)-JVS(547)*X(35)-JVS(548)*X(36)&
            &-JVS(549)*X(37)-JVS(550)*X(38)-JVS(551)*X(39)-JVS(552)*X(40)-JVS(553)*X(41)-JVS(554)*X(44)-JVS(555)*X(45)&
            &-JVS(556)*X(46)-JVS(557)*X(47)-JVS(558)*X(48)-JVS(559)*X(49)-JVS(560)*X(50)-JVS(561)*X(51)-JVS(562)*X(52)&
            &-JVS(563)*X(53)-JVS(564)*X(54)
  X(56) = X(56)-JVS(570)*X(16)-JVS(571)*X(42)-JVS(572)*X(44)-JVS(573)*X(45)-JVS(574)*X(46)-JVS(575)*X(47)-JVS(576)*X(49)&
            &-JVS(577)*X(50)-JVS(578)*X(51)-JVS(579)*X(52)-JVS(580)*X(53)-JVS(581)*X(54)-JVS(582)*X(55)
  X(57) = X(57)-JVS(587)*X(6)-JVS(588)*X(33)-JVS(589)*X(41)-JVS(590)*X(44)-JVS(591)*X(46)-JVS(592)*X(47)-JVS(593)*X(49)&
            &-JVS(594)*X(51)-JVS(595)*X(52)-JVS(596)*X(53)-JVS(597)*X(54)-JVS(598)*X(55)-JVS(599)*X(56)
  X(58) = X(58)-JVS(603)*X(13)-JVS(604)*X(19)-JVS(605)*X(35)-JVS(606)*X(36)-JVS(607)*X(37)-JVS(608)*X(38)-JVS(609)*X(39)&
            &-JVS(610)*X(40)-JVS(611)*X(41)-JVS(612)*X(44)-JVS(613)*X(45)-JVS(614)*X(46)-JVS(615)*X(49)-JVS(616)*X(50)&
            &-JVS(617)*X(51)-JVS(618)*X(52)-JVS(619)*X(53)-JVS(620)*X(54)-JVS(621)*X(55)-JVS(622)*X(56)-JVS(623)*X(57)
  X(59) = X(59)-JVS(626)*X(12)-JVS(627)*X(15)-JVS(628)*X(19)-JVS(629)*X(20)-JVS(630)*X(23)-JVS(631)*X(25)-JVS(632)*X(30)&
            &-JVS(633)*X(31)-JVS(634)*X(32)-JVS(635)*X(34)-JVS(636)*X(35)-JVS(637)*X(36)-JVS(638)*X(37)-JVS(639)*X(38)&
            &-JVS(640)*X(39)-JVS(641)*X(40)-JVS(642)*X(41)-JVS(643)*X(43)-JVS(644)*X(44)-JVS(645)*X(45)-JVS(646)*X(46)&
            &-JVS(647)*X(47)-JVS(648)*X(48)-JVS(649)*X(49)-JVS(650)*X(50)-JVS(651)*X(51)-JVS(652)*X(52)-JVS(653)*X(53)&
            &-JVS(654)*X(54)-JVS(655)*X(55)-JVS(656)*X(56)-JVS(657)*X(57)-JVS(658)*X(58)
  X(59) = X(59)/JVS(659)
  X(58) = (X(58)-JVS(625)*X(59))/(JVS(624))
  X(57) = (X(57)-JVS(601)*X(58)-JVS(602)*X(59))/(JVS(600))
  X(56) = (X(56)-JVS(584)*X(57)-JVS(585)*X(58)-JVS(586)*X(59))/(JVS(583))
  X(55) = (X(55)-JVS(566)*X(56)-JVS(567)*X(57)-JVS(568)*X(58)-JVS(569)*X(59))/(JVS(565))
  X(54) = (X(54)-JVS(537)*X(55)-JVS(538)*X(56)-JVS(539)*X(57)-JVS(540)*X(58)-JVS(541)*X(59))/(JVS(536))
  X(53) = (X(53)-JVS(506)*X(54)-JVS(507)*X(55)-JVS(508)*X(56)-JVS(509)*X(57)-JVS(510)*X(58)-JVS(511)*X(59))/(JVS(505))
  X(52) = (X(52)-JVS(458)*X(53)-JVS(459)*X(54)-JVS(460)*X(55)-JVS(461)*X(56)-JVS(462)*X(57)-JVS(463)*X(58)-JVS(464)&
            &*X(59))/(JVS(457))
  X(51) = (X(51)-JVS(428)*X(52)-JVS(429)*X(53)-JVS(430)*X(54)-JVS(431)*X(55)-JVS(432)*X(56)-JVS(433)*X(57)-JVS(434)&
            &*X(58)-JVS(435)*X(59))/(JVS(427))
  X(50) = (X(50)-JVS(386)*X(51)-JVS(387)*X(52)-JVS(388)*X(53)-JVS(389)*X(54)-JVS(390)*X(55)-JVS(391)*X(56)-JVS(392)&
            &*X(58)-JVS(393)*X(59))/(JVS(385))
  X(49) = (X(49)-JVS(372)*X(51)-JVS(373)*X(52)-JVS(374)*X(53)-JVS(375)*X(54)-JVS(376)*X(55)-JVS(377)*X(58)-JVS(378)&
            &*X(59))/(JVS(371))
  X(48) = (X(48)-JVS(358)*X(49)-JVS(359)*X(50)-JVS(360)*X(51)-JVS(361)*X(52)-JVS(362)*X(53)-JVS(363)*X(54)-JVS(364)&
            &*X(55)-JVS(365)*X(56)-JVS(366)*X(57)-JVS(367)*X(58)-JVS(368)*X(59))/(JVS(357))
  X(47) = (X(47)-JVS(336)*X(51)-JVS(337)*X(52)-JVS(338)*X(53)-JVS(339)*X(54)-JVS(340)*X(55)-JVS(341)*X(56)-JVS(342)&
            &*X(58)-JVS(343)*X(59))/(JVS(335))
  X(46) = (X(46)-JVS(325)*X(51)-JVS(326)*X(53)-JVS(327)*X(54)-JVS(328)*X(58)-JVS(329)*X(59))/(JVS(324))
  X(45) = (X(45)-JVS(312)*X(51)-JVS(313)*X(52)-JVS(314)*X(53)-JVS(315)*X(54)-JVS(316)*X(55)-JVS(317)*X(58)-JVS(318)&
            &*X(59))/(JVS(311))
  X(44) = (X(44)-JVS(301)*X(51)-JVS(302)*X(52)-JVS(303)*X(53)-JVS(304)*X(55)-JVS(305)*X(58))/(JVS(300))
  X(43) = (X(43)-JVS(285)*X(44)-JVS(286)*X(45)-JVS(287)*X(46)-JVS(288)*X(49)-JVS(289)*X(50)-JVS(290)*X(51)-JVS(291)&
            &*X(52)-JVS(292)*X(53)-JVS(293)*X(54)-JVS(294)*X(55)-JVS(295)*X(56)-JVS(296)*X(57)-JVS(297)*X(58)-JVS(298)&
            &*X(59))/(JVS(284))
  X(42) = (X(42)-JVS(254)*X(44)-JVS(255)*X(45)-JVS(256)*X(46)-JVS(257)*X(49)-JVS(258)*X(50)-JVS(259)*X(51)-JVS(260)&
            &*X(52)-JVS(261)*X(53)-JVS(262)*X(54)-JVS(263)*X(55)-JVS(264)*X(56)-JVS(265)*X(57)-JVS(266)*X(58)-JVS(267)&
            &*X(59))/(JVS(253))
  X(41) = (X(41)-JVS(241)*X(51)-JVS(242)*X(52)-JVS(243)*X(53)-JVS(244)*X(55)-JVS(245)*X(58))/(JVS(240))
  X(40) = (X(40)-JVS(232)*X(46)-JVS(233)*X(51)-JVS(234)*X(52)-JVS(235)*X(53)-JVS(236)*X(54)-JVS(237)*X(55)-JVS(238)&
            &*X(58))/(JVS(231))
  X(39) = (X(39)-JVS(222)*X(45)-JVS(223)*X(51)-JVS(224)*X(52)-JVS(225)*X(53)-JVS(226)*X(54)-JVS(227)*X(55)-JVS(228)&
            &*X(58)-JVS(229)*X(59))/(JVS(221))
  X(38) = (X(38)-JVS(210)*X(46)-JVS(211)*X(51)-JVS(212)*X(52)-JVS(213)*X(53)-JVS(214)*X(54)-JVS(215)*X(55)-JVS(216)&
            &*X(58))/(JVS(209))
  X(37) = (X(37)-JVS(198)*X(46)-JVS(199)*X(51)-JVS(200)*X(52)-JVS(201)*X(53)-JVS(202)*X(54)-JVS(203)*X(55)-JVS(204)&
            &*X(58))/(JVS(197))
  X(36) = (X(36)-JVS(190)*X(51)-JVS(191)*X(52)-JVS(192)*X(53)-JVS(193)*X(55)-JVS(194)*X(58))/(JVS(189))
  X(35) = (X(35)-JVS(183)*X(51)-JVS(184)*X(52)-JVS(185)*X(53)-JVS(186)*X(55)-JVS(187)*X(58))/(JVS(182))
  X(34) = (X(34)-JVS(172)*X(35)-JVS(173)*X(36)-JVS(174)*X(45)-JVS(175)*X(50)-JVS(176)*X(52)-JVS(177)*X(53)-JVS(178)&
            &*X(54)-JVS(179)*X(55)-JVS(180)*X(58))/(JVS(171))
  X(33) = (X(33)-JVS(160)*X(41)-JVS(161)*X(44)-JVS(162)*X(46)-JVS(163)*X(47)-JVS(164)*X(49)-JVS(165)*X(52)-JVS(166)&
            &*X(53)-JVS(167)*X(54)-JVS(168)*X(55)-JVS(169)*X(56)-JVS(170)*X(58))/(JVS(159))
  X(32) = (X(32)-JVS(151)*X(51)-JVS(152)*X(52)-JVS(153)*X(53)-JVS(154)*X(54)-JVS(155)*X(55)-JVS(156)*X(59))/(JVS(150))
  X(31) = (X(31)-JVS(142)*X(35)-JVS(143)*X(45)-JVS(144)*X(52)-JVS(145)*X(53)-JVS(146)*X(54)-JVS(147)*X(55)-JVS(148)&
            &*X(58))/(JVS(141))
  X(30) = (X(30)-JVS(134)*X(35)-JVS(135)*X(36)-JVS(136)*X(52)-JVS(137)*X(53)-JVS(138)*X(54)-JVS(139)*X(55)-JVS(140)&
            &*X(58))/(JVS(133))
  X(29) = (X(29)-JVS(130)*X(46)-JVS(131)*X(53)-JVS(132)*X(54))/(JVS(129))
  X(28) = (X(28)-JVS(126)*X(46)-JVS(127)*X(53)-JVS(128)*X(54))/(JVS(125))
  X(27) = (X(27)-JVS(122)*X(46)-JVS(123)*X(53)-JVS(124)*X(54))/(JVS(121))
  X(26) = (X(26)-JVS(107)*X(27)-JVS(108)*X(28)-JVS(109)*X(29)-JVS(110)*X(31)-JVS(111)*X(34)-JVS(112)*X(43)-JVS(113)&
            &*X(45)-JVS(114)*X(46)-JVS(115)*X(48)-JVS(116)*X(52)-JVS(117)*X(53)-JVS(118)*X(54)-JVS(119)*X(55)-JVS(120)&
            &*X(58))/(JVS(106))
  X(25) = (X(25)-JVS(96)*X(30)-JVS(97)*X(31)-JVS(98)*X(34)-JVS(99)*X(43)-JVS(100)*X(48)-JVS(101)*X(51)-JVS(102)*X(53)&
            &-JVS(103)*X(54)-JVS(104)*X(59))/(JVS(95))
  X(24) = (X(24)-JVS(90)*X(46)-JVS(91)*X(53)-JVS(92)*X(54))/(JVS(89))
  X(23) = (X(23)-JVS(86)*X(53)-JVS(87)*X(55)-JVS(88)*X(59))/(JVS(85))
  X(22) = (X(22)-JVS(83)*X(53)-JVS(84)*X(54))/(JVS(82))
  X(21) = (X(21)-JVS(77)*X(51)-JVS(78)*X(52)-JVS(79)*X(53))/(JVS(76))
  X(20) = (X(20)-JVS(73)*X(51)-JVS(74)*X(53)-JVS(75)*X(59))/(JVS(72))
  X(19) = (X(19)-JVS(69)*X(46)-JVS(70)*X(54)-JVS(71)*X(59))/(JVS(68))
  X(18) = (X(18)-JVS(64)*X(51)-JVS(65)*X(53)-JVS(66)*X(55))/(JVS(63))
  X(17) = (X(17)-JVS(59)*X(28)-JVS(60)*X(29)-JVS(61)*X(46)-JVS(62)*X(53))/(JVS(58))
  X(16) = (X(16)-JVS(57)*X(53))/(JVS(56))
  X(15) = (X(15)-JVS(54)*X(54)-JVS(55)*X(59))/(JVS(53))
  X(14) = (X(14)-JVS(51)*X(51)-JVS(52)*X(53))/(JVS(50))
  X(13) = (X(13)-JVS(48)*X(53)-JVS(49)*X(58))/(JVS(47))
  X(12) = (X(12)-JVS(45)*X(45)-JVS(46)*X(59))/(JVS(44))
  X(11) = (X(11)-JVS(43)*X(53))/(JVS(42))
  X(10) = (X(10)-JVS(41)*X(53))/(JVS(40))
  X(9) = (X(9)-JVS(39)*X(53))/(JVS(38))
  X(8) = (X(8)-JVS(37)*X(53))/(JVS(36))
  X(7) = (X(7)-JVS(35)*X(46))/(JVS(34))
  X(6) = (X(6)-JVS(33)*X(53))/(JVS(32))
  X(5) = (X(5)-JVS(31)*X(53))/(JVS(30))
  X(4) = (X(4)-JVS(28)*X(26)-JVS(29)*X(53))/(JVS(27))
  X(3) = (X(3)-JVS(11)*X(27)-JVS(12)*X(28)-JVS(13)*X(29)-JVS(14)*X(37)-JVS(15)*X(38)-JVS(16)*X(40)-JVS(17)*X(41)-JVS(18)&
           &*X(44)-JVS(19)*X(45)-JVS(20)*X(46)-JVS(21)*X(49)-JVS(22)*X(50)-JVS(23)*X(52)-JVS(24)*X(55)-JVS(25)*X(56)-JVS(26)&
           &*X(57))/(JVS(10))
  X(2) = (X(2)-JVS(5)*X(24)-JVS(6)*X(27)-JVS(7)*X(28)-JVS(8)*X(29)-JVS(9)*X(46))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(5)-JVS(3)*X(53))/(JVS(1))
      
END SUBROUTINE radm2sorg_KppSolve
























      SUBROUTINE radm2sorg_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE radm2sorg_WCOPY



      SUBROUTINE radm2sorg_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE radm2sorg_WAXPY




      SUBROUTINE radm2sorg_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE radm2sorg_WSCAL


      REAL(kind=dp) FUNCTION radm2sorg_WLAMCH( C )








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
          CALL radm2sorg_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      radm2sorg_WLAMCH = Eps

      END FUNCTION radm2sorg_WLAMCH
     
      SUBROUTINE radm2sorg_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE radm2sorg_WLAMCH_ADD




      SUBROUTINE radm2sorg_SET2ZERO(N,Y)




      
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

      END SUBROUTINE radm2sorg_SET2ZERO



      REAL(kind=dp) FUNCTION radm2sorg_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      radm2sorg_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        radm2sorg_WDOT = radm2sorg_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         radm2sorg_WDOT = radm2sorg_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          radm2sorg_WDOT = radm2sorg_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        radm2sorg_WDOT = radm2sorg_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION radm2sorg_WDOT                                          




   SUBROUTINE decomp_radm2sorg( JVS, IER )
   
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
   W( 53 ) = JVS( 3 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 5 )
  JVS( 3) = W( 53 )
  IF ( ABS(  JVS( 4 )) < TINY(a) ) THEN
         IER = 2                                       
         RETURN
  END IF
   W( 2 ) = JVS( 4 )
   W( 24 ) = JVS( 5 )
   W( 27 ) = JVS( 6 )
   W( 28 ) = JVS( 7 )
   W( 29 ) = JVS( 8 )
   W( 46 ) = JVS( 9 )
  JVS( 4) = W( 2 )
  JVS( 5) = W( 24 )
  JVS( 6) = W( 27 )
  JVS( 7) = W( 28 )
  JVS( 8) = W( 29 )
  JVS( 9) = W( 46 )
  IF ( ABS(  JVS( 10 )) < TINY(a) ) THEN
         IER = 3                                       
         RETURN
  END IF
   W( 3 ) = JVS( 10 )
   W( 27 ) = JVS( 11 )
   W( 28 ) = JVS( 12 )
   W( 29 ) = JVS( 13 )
   W( 37 ) = JVS( 14 )
   W( 38 ) = JVS( 15 )
   W( 40 ) = JVS( 16 )
   W( 41 ) = JVS( 17 )
   W( 44 ) = JVS( 18 )
   W( 45 ) = JVS( 19 )
   W( 46 ) = JVS( 20 )
   W( 49 ) = JVS( 21 )
   W( 50 ) = JVS( 22 )
   W( 52 ) = JVS( 23 )
   W( 55 ) = JVS( 24 )
   W( 56 ) = JVS( 25 )
   W( 57 ) = JVS( 26 )
  JVS( 10) = W( 3 )
  JVS( 11) = W( 27 )
  JVS( 12) = W( 28 )
  JVS( 13) = W( 29 )
  JVS( 14) = W( 37 )
  JVS( 15) = W( 38 )
  JVS( 16) = W( 40 )
  JVS( 17) = W( 41 )
  JVS( 18) = W( 44 )
  JVS( 19) = W( 45 )
  JVS( 20) = W( 46 )
  JVS( 21) = W( 49 )
  JVS( 22) = W( 50 )
  JVS( 23) = W( 52 )
  JVS( 24) = W( 55 )
  JVS( 25) = W( 56 )
  JVS( 26) = W( 57 )
  IF ( ABS(  JVS( 27 )) < TINY(a) ) THEN
         IER = 4                                       
         RETURN
  END IF
   W( 4 ) = JVS( 27 )
   W( 26 ) = JVS( 28 )
   W( 53 ) = JVS( 29 )
  JVS( 27) = W( 4 )
  JVS( 28) = W( 26 )
  JVS( 29) = W( 53 )
  IF ( ABS(  JVS( 30 )) < TINY(a) ) THEN
         IER = 5                                       
         RETURN
  END IF
   W( 5 ) = JVS( 30 )
   W( 53 ) = JVS( 31 )
  JVS( 30) = W( 5 )
  JVS( 31) = W( 53 )
  IF ( ABS(  JVS( 32 )) < TINY(a) ) THEN
         IER = 6                                       
         RETURN
  END IF
   W( 6 ) = JVS( 32 )
   W( 53 ) = JVS( 33 )
  JVS( 32) = W( 6 )
  JVS( 33) = W( 53 )
  IF ( ABS(  JVS( 34 )) < TINY(a) ) THEN
         IER = 7                                       
         RETURN
  END IF
   W( 7 ) = JVS( 34 )
   W( 46 ) = JVS( 35 )
  JVS( 34) = W( 7 )
  JVS( 35) = W( 46 )
  IF ( ABS(  JVS( 36 )) < TINY(a) ) THEN
         IER = 8                                       
         RETURN
  END IF
   W( 8 ) = JVS( 36 )
   W( 53 ) = JVS( 37 )
  JVS( 36) = W( 8 )
  JVS( 37) = W( 53 )
  IF ( ABS(  JVS( 38 )) < TINY(a) ) THEN
         IER = 9                                       
         RETURN
  END IF
   W( 9 ) = JVS( 38 )
   W( 53 ) = JVS( 39 )
  JVS( 38) = W( 9 )
  JVS( 39) = W( 53 )
  IF ( ABS(  JVS( 40 )) < TINY(a) ) THEN
         IER = 10                                      
         RETURN
  END IF
   W( 10 ) = JVS( 40 )
   W( 53 ) = JVS( 41 )
  JVS( 40) = W( 10 )
  JVS( 41) = W( 53 )
  IF ( ABS(  JVS( 42 )) < TINY(a) ) THEN
         IER = 11                                      
         RETURN
  END IF
   W( 11 ) = JVS( 42 )
   W( 53 ) = JVS( 43 )
  JVS( 42) = W( 11 )
  JVS( 43) = W( 53 )
  IF ( ABS(  JVS( 44 )) < TINY(a) ) THEN
         IER = 12                                      
         RETURN
  END IF
   W( 12 ) = JVS( 44 )
   W( 45 ) = JVS( 45 )
   W( 59 ) = JVS( 46 )
  JVS( 44) = W( 12 )
  JVS( 45) = W( 45 )
  JVS( 46) = W( 59 )
  IF ( ABS(  JVS( 47 )) < TINY(a) ) THEN
         IER = 13                                      
         RETURN
  END IF
   W( 13 ) = JVS( 47 )
   W( 53 ) = JVS( 48 )
   W( 58 ) = JVS( 49 )
  JVS( 47) = W( 13 )
  JVS( 48) = W( 53 )
  JVS( 49) = W( 58 )
  IF ( ABS(  JVS( 50 )) < TINY(a) ) THEN
         IER = 14                                      
         RETURN
  END IF
   W( 14 ) = JVS( 50 )
   W( 51 ) = JVS( 51 )
   W( 53 ) = JVS( 52 )
  JVS( 50) = W( 14 )
  JVS( 51) = W( 51 )
  JVS( 52) = W( 53 )
  IF ( ABS(  JVS( 53 )) < TINY(a) ) THEN
         IER = 15                                      
         RETURN
  END IF
   W( 15 ) = JVS( 53 )
   W( 54 ) = JVS( 54 )
   W( 59 ) = JVS( 55 )
  JVS( 53) = W( 15 )
  JVS( 54) = W( 54 )
  JVS( 55) = W( 59 )
  IF ( ABS(  JVS( 56 )) < TINY(a) ) THEN
         IER = 16                                      
         RETURN
  END IF
   W( 16 ) = JVS( 56 )
   W( 53 ) = JVS( 57 )
  JVS( 56) = W( 16 )
  JVS( 57) = W( 53 )
  IF ( ABS(  JVS( 58 )) < TINY(a) ) THEN
         IER = 17                                      
         RETURN
  END IF
   W( 17 ) = JVS( 58 )
   W( 28 ) = JVS( 59 )
   W( 29 ) = JVS( 60 )
   W( 46 ) = JVS( 61 )
   W( 53 ) = JVS( 62 )
  JVS( 58) = W( 17 )
  JVS( 59) = W( 28 )
  JVS( 60) = W( 29 )
  JVS( 61) = W( 46 )
  JVS( 62) = W( 53 )
  IF ( ABS(  JVS( 63 )) < TINY(a) ) THEN
         IER = 18                                      
         RETURN
  END IF
   W( 18 ) = JVS( 63 )
   W( 51 ) = JVS( 64 )
   W( 53 ) = JVS( 65 )
   W( 55 ) = JVS( 66 )
  JVS( 63) = W( 18 )
  JVS( 64) = W( 51 )
  JVS( 65) = W( 53 )
  JVS( 66) = W( 55 )
  IF ( ABS(  JVS( 68 )) < TINY(a) ) THEN
         IER = 19                                      
         RETURN
  END IF
   W( 7 ) = JVS( 67 )
   W( 19 ) = JVS( 68 )
   W( 46 ) = JVS( 69 )
   W( 54 ) = JVS( 70 )
   W( 59 ) = JVS( 71 )
  a = -W( 7 ) / JVS(           34  )
  W( 7 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 35 )
  JVS( 67) = W( 7 )
  JVS( 68) = W( 19 )
  JVS( 69) = W( 46 )
  JVS( 70) = W( 54 )
  JVS( 71) = W( 59 )
  IF ( ABS(  JVS( 72 )) < TINY(a) ) THEN
         IER = 20                                      
         RETURN
  END IF
   W( 20 ) = JVS( 72 )
   W( 51 ) = JVS( 73 )
   W( 53 ) = JVS( 74 )
   W( 59 ) = JVS( 75 )
  JVS( 72) = W( 20 )
  JVS( 73) = W( 51 )
  JVS( 74) = W( 53 )
  JVS( 75) = W( 59 )
  IF ( ABS(  JVS( 76 )) < TINY(a) ) THEN
         IER = 21                                      
         RETURN
  END IF
   W( 21 ) = JVS( 76 )
   W( 51 ) = JVS( 77 )
   W( 52 ) = JVS( 78 )
   W( 53 ) = JVS( 79 )
  JVS( 76) = W( 21 )
  JVS( 77) = W( 51 )
  JVS( 78) = W( 52 )
  JVS( 79) = W( 53 )
  IF ( ABS(  JVS( 82 )) < TINY(a) ) THEN
         IER = 22                                      
         RETURN
  END IF
   W( 10 ) = JVS( 80 )
   W( 11 ) = JVS( 81 )
   W( 22 ) = JVS( 82 )
   W( 53 ) = JVS( 83 )
   W( 54 ) = JVS( 84 )
  a = -W( 10 ) / JVS(           40  )
  W( 10 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 41 )
  a = -W( 11 ) / JVS(           42  )
  W( 11 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 43 )
  JVS( 80) = W( 10 )
  JVS( 81) = W( 11 )
  JVS( 82) = W( 22 )
  JVS( 83) = W( 53 )
  JVS( 84) = W( 54 )
  IF ( ABS(  JVS( 85 )) < TINY(a) ) THEN
         IER = 23                                      
         RETURN
  END IF
   W( 23 ) = JVS( 85 )
   W( 53 ) = JVS( 86 )
   W( 55 ) = JVS( 87 )
   W( 59 ) = JVS( 88 )
  JVS( 85) = W( 23 )
  JVS( 86) = W( 53 )
  JVS( 87) = W( 55 )
  JVS( 88) = W( 59 )
  IF ( ABS(  JVS( 89 )) < TINY(a) ) THEN
         IER = 24                                      
         RETURN
  END IF
   W( 24 ) = JVS( 89 )
   W( 46 ) = JVS( 90 )
   W( 53 ) = JVS( 91 )
   W( 54 ) = JVS( 92 )
  JVS( 89) = W( 24 )
  JVS( 90) = W( 46 )
  JVS( 91) = W( 53 )
  JVS( 92) = W( 54 )
  IF ( ABS(  JVS( 95 )) < TINY(a) ) THEN
         IER = 25                                      
         RETURN
  END IF
   W( 15 ) = JVS( 93 )
   W( 22 ) = JVS( 94 )
   W( 25 ) = JVS( 95 )
   W( 30 ) = JVS( 96 )
   W( 31 ) = JVS( 97 )
   W( 34 ) = JVS( 98 )
   W( 43 ) = JVS( 99 )
   W( 48 ) = JVS( 100 )
   W( 51 ) = JVS( 101 )
   W( 53 ) = JVS( 102 )
   W( 54 ) = JVS( 103 )
   W( 59 ) = JVS( 104 )
  a = -W( 15 ) / JVS(           53  )
  W( 15 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 54 )
  W( 59 ) = W( 59 ) + a*JVS( 55 )
  a = -W( 22 ) / JVS(           82  )
  W( 22 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 83 )
  W( 54 ) = W( 54 ) + a*JVS( 84 )
  JVS( 93) = W( 15 )
  JVS( 94) = W( 22 )
  JVS( 95) = W( 25 )
  JVS( 96) = W( 30 )
  JVS( 97) = W( 31 )
  JVS( 98) = W( 34 )
  JVS( 99) = W( 43 )
  JVS( 100) = W( 48 )
  JVS( 101) = W( 51 )
  JVS( 102) = W( 53 )
  JVS( 103) = W( 54 )
  JVS( 104) = W( 59 )
  IF ( ABS(  JVS( 106 )) < TINY(a) ) THEN
         IER = 26                                      
         RETURN
  END IF
   W( 24 ) = JVS( 105 )
   W( 26 ) = JVS( 106 )
   W( 27 ) = JVS( 107 )
   W( 28 ) = JVS( 108 )
   W( 29 ) = JVS( 109 )
   W( 31 ) = JVS( 110 )
   W( 34 ) = JVS( 111 )
   W( 43 ) = JVS( 112 )
   W( 45 ) = JVS( 113 )
   W( 46 ) = JVS( 114 )
   W( 48 ) = JVS( 115 )
   W( 52 ) = JVS( 116 )
   W( 53 ) = JVS( 117 )
   W( 54 ) = JVS( 118 )
   W( 55 ) = JVS( 119 )
   W( 58 ) = JVS( 120 )
  a = -W( 24 ) / JVS(           89  )
  W( 24 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 90 )
  W( 53 ) = W( 53 ) + a*JVS( 91 )
  W( 54 ) = W( 54 ) + a*JVS( 92 )
  JVS( 105) = W( 24 )
  JVS( 106) = W( 26 )
  JVS( 107) = W( 27 )
  JVS( 108) = W( 28 )
  JVS( 109) = W( 29 )
  JVS( 110) = W( 31 )
  JVS( 111) = W( 34 )
  JVS( 112) = W( 43 )
  JVS( 113) = W( 45 )
  JVS( 114) = W( 46 )
  JVS( 115) = W( 48 )
  JVS( 116) = W( 52 )
  JVS( 117) = W( 53 )
  JVS( 118) = W( 54 )
  JVS( 119) = W( 55 )
  JVS( 120) = W( 58 )
  IF ( ABS(  JVS( 121 )) < TINY(a) ) THEN
         IER = 27                                      
         RETURN
  END IF
   W( 27 ) = JVS( 121 )
   W( 46 ) = JVS( 122 )
   W( 53 ) = JVS( 123 )
   W( 54 ) = JVS( 124 )
  JVS( 121) = W( 27 )
  JVS( 122) = W( 46 )
  JVS( 123) = W( 53 )
  JVS( 124) = W( 54 )
  IF ( ABS(  JVS( 125 )) < TINY(a) ) THEN
         IER = 28                                      
         RETURN
  END IF
   W( 28 ) = JVS( 125 )
   W( 46 ) = JVS( 126 )
   W( 53 ) = JVS( 127 )
   W( 54 ) = JVS( 128 )
  JVS( 125) = W( 28 )
  JVS( 126) = W( 46 )
  JVS( 127) = W( 53 )
  JVS( 128) = W( 54 )
  IF ( ABS(  JVS( 129 )) < TINY(a) ) THEN
         IER = 29                                      
         RETURN
  END IF
   W( 29 ) = JVS( 129 )
   W( 46 ) = JVS( 130 )
   W( 53 ) = JVS( 131 )
   W( 54 ) = JVS( 132 )
  JVS( 129) = W( 29 )
  JVS( 130) = W( 46 )
  JVS( 131) = W( 53 )
  JVS( 132) = W( 54 )
  IF ( ABS(  JVS( 133 )) < TINY(a) ) THEN
         IER = 30                                      
         RETURN
  END IF
   W( 30 ) = JVS( 133 )
   W( 35 ) = JVS( 134 )
   W( 36 ) = JVS( 135 )
   W( 52 ) = JVS( 136 )
   W( 53 ) = JVS( 137 )
   W( 54 ) = JVS( 138 )
   W( 55 ) = JVS( 139 )
   W( 58 ) = JVS( 140 )
  JVS( 133) = W( 30 )
  JVS( 134) = W( 35 )
  JVS( 135) = W( 36 )
  JVS( 136) = W( 52 )
  JVS( 137) = W( 53 )
  JVS( 138) = W( 54 )
  JVS( 139) = W( 55 )
  JVS( 140) = W( 58 )
  IF ( ABS(  JVS( 141 )) < TINY(a) ) THEN
         IER = 31                                      
         RETURN
  END IF
   W( 31 ) = JVS( 141 )
   W( 35 ) = JVS( 142 )
   W( 45 ) = JVS( 143 )
   W( 52 ) = JVS( 144 )
   W( 53 ) = JVS( 145 )
   W( 54 ) = JVS( 146 )
   W( 55 ) = JVS( 147 )
   W( 58 ) = JVS( 148 )
  JVS( 141) = W( 31 )
  JVS( 142) = W( 35 )
  JVS( 143) = W( 45 )
  JVS( 144) = W( 52 )
  JVS( 145) = W( 53 )
  JVS( 146) = W( 54 )
  JVS( 147) = W( 55 )
  JVS( 148) = W( 58 )
  IF ( ABS(  JVS( 150 )) < TINY(a) ) THEN
         IER = 32                                      
         RETURN
  END IF
   W( 22 ) = JVS( 149 )
   W( 32 ) = JVS( 150 )
   W( 51 ) = JVS( 151 )
   W( 52 ) = JVS( 152 )
   W( 53 ) = JVS( 153 )
   W( 54 ) = JVS( 154 )
   W( 55 ) = JVS( 155 )
   W( 59 ) = JVS( 156 )
  a = -W( 22 ) / JVS(           82  )
  W( 22 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 83 )
  W( 54 ) = W( 54 ) + a*JVS( 84 )
  JVS( 149) = W( 22 )
  JVS( 150) = W( 32 )
  JVS( 151) = W( 51 )
  JVS( 152) = W( 52 )
  JVS( 153) = W( 53 )
  JVS( 154) = W( 54 )
  JVS( 155) = W( 55 )
  JVS( 156) = W( 59 )
  IF ( ABS(  JVS( 159 )) < TINY(a) ) THEN
         IER = 33                                      
         RETURN
  END IF
   W( 16 ) = JVS( 157 )
   W( 29 ) = JVS( 158 )
   W( 33 ) = JVS( 159 )
   W( 41 ) = JVS( 160 )
   W( 44 ) = JVS( 161 )
   W( 46 ) = JVS( 162 )
   W( 47 ) = JVS( 163 )
   W( 49 ) = JVS( 164 )
   W( 52 ) = JVS( 165 )
   W( 53 ) = JVS( 166 )
   W( 54 ) = JVS( 167 )
   W( 55 ) = JVS( 168 )
   W( 56 ) = JVS( 169 )
   W( 58 ) = JVS( 170 )
  a = -W( 16 ) / JVS(           56  )
  W( 16 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 57 )
  a = -W( 29 ) / JVS(          129  )
  W( 29 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 130 )
  W( 53 ) = W( 53 ) + a*JVS( 131 )
  W( 54 ) = W( 54 ) + a*JVS( 132 )
  JVS( 157) = W( 16 )
  JVS( 158) = W( 29 )
  JVS( 159) = W( 33 )
  JVS( 160) = W( 41 )
  JVS( 161) = W( 44 )
  JVS( 162) = W( 46 )
  JVS( 163) = W( 47 )
  JVS( 164) = W( 49 )
  JVS( 165) = W( 52 )
  JVS( 166) = W( 53 )
  JVS( 167) = W( 54 )
  JVS( 168) = W( 55 )
  JVS( 169) = W( 56 )
  JVS( 170) = W( 58 )
  IF ( ABS(  JVS( 171 )) < TINY(a) ) THEN
         IER = 34                                      
         RETURN
  END IF
   W( 34 ) = JVS( 171 )
   W( 35 ) = JVS( 172 )
   W( 36 ) = JVS( 173 )
   W( 45 ) = JVS( 174 )
   W( 50 ) = JVS( 175 )
   W( 52 ) = JVS( 176 )
   W( 53 ) = JVS( 177 )
   W( 54 ) = JVS( 178 )
   W( 55 ) = JVS( 179 )
   W( 58 ) = JVS( 180 )
  JVS( 171) = W( 34 )
  JVS( 172) = W( 35 )
  JVS( 173) = W( 36 )
  JVS( 174) = W( 45 )
  JVS( 175) = W( 50 )
  JVS( 176) = W( 52 )
  JVS( 177) = W( 53 )
  JVS( 178) = W( 54 )
  JVS( 179) = W( 55 )
  JVS( 180) = W( 58 )
  IF ( ABS(  JVS( 182 )) < TINY(a) ) THEN
         IER = 35                                      
         RETURN
  END IF
   W( 10 ) = JVS( 181 )
   W( 35 ) = JVS( 182 )
   W( 51 ) = JVS( 183 )
   W( 52 ) = JVS( 184 )
   W( 53 ) = JVS( 185 )
   W( 55 ) = JVS( 186 )
   W( 58 ) = JVS( 187 )
  a = -W( 10 ) / JVS(           40  )
  W( 10 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 41 )
  JVS( 181) = W( 10 )
  JVS( 182) = W( 35 )
  JVS( 183) = W( 51 )
  JVS( 184) = W( 52 )
  JVS( 185) = W( 53 )
  JVS( 186) = W( 55 )
  JVS( 187) = W( 58 )
  IF ( ABS(  JVS( 189 )) < TINY(a) ) THEN
         IER = 36                                      
         RETURN
  END IF
   W( 11 ) = JVS( 188 )
   W( 36 ) = JVS( 189 )
   W( 51 ) = JVS( 190 )
   W( 52 ) = JVS( 191 )
   W( 53 ) = JVS( 192 )
   W( 55 ) = JVS( 193 )
   W( 58 ) = JVS( 194 )
  a = -W( 11 ) / JVS(           42  )
  W( 11 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 43 )
  JVS( 188) = W( 11 )
  JVS( 189) = W( 36 )
  JVS( 190) = W( 51 )
  JVS( 191) = W( 52 )
  JVS( 192) = W( 53 )
  JVS( 193) = W( 55 )
  JVS( 194) = W( 58 )
  IF ( ABS(  JVS( 197 )) < TINY(a) ) THEN
         IER = 37                                      
         RETURN
  END IF
   W( 27 ) = JVS( 195 )
   W( 28 ) = JVS( 196 )
   W( 37 ) = JVS( 197 )
   W( 46 ) = JVS( 198 )
   W( 51 ) = JVS( 199 )
   W( 52 ) = JVS( 200 )
   W( 53 ) = JVS( 201 )
   W( 54 ) = JVS( 202 )
   W( 55 ) = JVS( 203 )
   W( 58 ) = JVS( 204 )
  a = -W( 27 ) / JVS(          121  )
  W( 27 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 122 )
  W( 53 ) = W( 53 ) + a*JVS( 123 )
  W( 54 ) = W( 54 ) + a*JVS( 124 )
  a = -W( 28 ) / JVS(          125  )
  W( 28 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 126 )
  W( 53 ) = W( 53 ) + a*JVS( 127 )
  W( 54 ) = W( 54 ) + a*JVS( 128 )
  JVS( 195) = W( 27 )
  JVS( 196) = W( 28 )
  JVS( 197) = W( 37 )
  JVS( 198) = W( 46 )
  JVS( 199) = W( 51 )
  JVS( 200) = W( 52 )
  JVS( 201) = W( 53 )
  JVS( 202) = W( 54 )
  JVS( 203) = W( 55 )
  JVS( 204) = W( 58 )
  IF ( ABS(  JVS( 209 )) < TINY(a) ) THEN
         IER = 38                                      
         RETURN
  END IF
   W( 24 ) = JVS( 205 )
   W( 27 ) = JVS( 206 )
   W( 28 ) = JVS( 207 )
   W( 29 ) = JVS( 208 )
   W( 38 ) = JVS( 209 )
   W( 46 ) = JVS( 210 )
   W( 51 ) = JVS( 211 )
   W( 52 ) = JVS( 212 )
   W( 53 ) = JVS( 213 )
   W( 54 ) = JVS( 214 )
   W( 55 ) = JVS( 215 )
   W( 58 ) = JVS( 216 )
  a = -W( 24 ) / JVS(           89  )
  W( 24 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 90 )
  W( 53 ) = W( 53 ) + a*JVS( 91 )
  W( 54 ) = W( 54 ) + a*JVS( 92 )
  a = -W( 27 ) / JVS(          121  )
  W( 27 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 122 )
  W( 53 ) = W( 53 ) + a*JVS( 123 )
  W( 54 ) = W( 54 ) + a*JVS( 124 )
  a = -W( 28 ) / JVS(          125  )
  W( 28 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 126 )
  W( 53 ) = W( 53 ) + a*JVS( 127 )
  W( 54 ) = W( 54 ) + a*JVS( 128 )
  a = -W( 29 ) / JVS(          129  )
  W( 29 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 130 )
  W( 53 ) = W( 53 ) + a*JVS( 131 )
  W( 54 ) = W( 54 ) + a*JVS( 132 )
  JVS( 205) = W( 24 )
  JVS( 206) = W( 27 )
  JVS( 207) = W( 28 )
  JVS( 208) = W( 29 )
  JVS( 209) = W( 38 )
  JVS( 210) = W( 46 )
  JVS( 211) = W( 51 )
  JVS( 212) = W( 52 )
  JVS( 213) = W( 53 )
  JVS( 214) = W( 54 )
  JVS( 215) = W( 55 )
  JVS( 216) = W( 58 )
  IF ( ABS(  JVS( 221 )) < TINY(a) ) THEN
         IER = 39                                      
         RETURN
  END IF
   W( 8 ) = JVS( 217 )
   W( 9 ) = JVS( 218 )
   W( 22 ) = JVS( 219 )
   W( 23 ) = JVS( 220 )
   W( 39 ) = JVS( 221 )
   W( 45 ) = JVS( 222 )
   W( 51 ) = JVS( 223 )
   W( 52 ) = JVS( 224 )
   W( 53 ) = JVS( 225 )
   W( 54 ) = JVS( 226 )
   W( 55 ) = JVS( 227 )
   W( 58 ) = JVS( 228 )
   W( 59 ) = JVS( 229 )
  a = -W( 8 ) / JVS(           36  )
  W( 8 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 37 )
  a = -W( 9 ) / JVS(           38  )
  W( 9 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 39 )
  a = -W( 22 ) / JVS(           82  )
  W( 22 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 83 )
  W( 54 ) = W( 54 ) + a*JVS( 84 )
  a = -W( 23 ) / JVS(           85  )
  W( 23 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 86 )
  W( 55 ) = W( 55 ) + a*JVS( 87 )
  W( 59 ) = W( 59 ) + a*JVS( 88 )
  JVS( 217) = W( 8 )
  JVS( 218) = W( 9 )
  JVS( 219) = W( 22 )
  JVS( 220) = W( 23 )
  JVS( 221) = W( 39 )
  JVS( 222) = W( 45 )
  JVS( 223) = W( 51 )
  JVS( 224) = W( 52 )
  JVS( 225) = W( 53 )
  JVS( 226) = W( 54 )
  JVS( 227) = W( 55 )
  JVS( 228) = W( 58 )
  JVS( 229) = W( 59 )
  IF ( ABS(  JVS( 231 )) < TINY(a) ) THEN
         IER = 40                                      
         RETURN
  END IF
   W( 24 ) = JVS( 230 )
   W( 40 ) = JVS( 231 )
   W( 46 ) = JVS( 232 )
   W( 51 ) = JVS( 233 )
   W( 52 ) = JVS( 234 )
   W( 53 ) = JVS( 235 )
   W( 54 ) = JVS( 236 )
   W( 55 ) = JVS( 237 )
   W( 58 ) = JVS( 238 )
  a = -W( 24 ) / JVS(           89  )
  W( 24 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 90 )
  W( 53 ) = W( 53 ) + a*JVS( 91 )
  W( 54 ) = W( 54 ) + a*JVS( 92 )
  JVS( 230) = W( 24 )
  JVS( 231) = W( 40 )
  JVS( 232) = W( 46 )
  JVS( 233) = W( 51 )
  JVS( 234) = W( 52 )
  JVS( 235) = W( 53 )
  JVS( 236) = W( 54 )
  JVS( 237) = W( 55 )
  JVS( 238) = W( 58 )
  IF ( ABS(  JVS( 240 )) < TINY(a) ) THEN
         IER = 41                                      
         RETURN
  END IF
   W( 8 ) = JVS( 239 )
   W( 41 ) = JVS( 240 )
   W( 51 ) = JVS( 241 )
   W( 52 ) = JVS( 242 )
   W( 53 ) = JVS( 243 )
   W( 55 ) = JVS( 244 )
   W( 58 ) = JVS( 245 )
  a = -W( 8 ) / JVS(           36  )
  W( 8 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 37 )
  JVS( 239) = W( 8 )
  JVS( 240) = W( 41 )
  JVS( 241) = W( 51 )
  JVS( 242) = W( 52 )
  JVS( 243) = W( 53 )
  JVS( 244) = W( 55 )
  JVS( 245) = W( 58 )
  IF ( ABS(  JVS( 253 )) < TINY(a) ) THEN
         IER = 42                                      
         RETURN
  END IF
   W( 32 ) = JVS( 246 )
   W( 35 ) = JVS( 247 )
   W( 36 ) = JVS( 248 )
   W( 37 ) = JVS( 249 )
   W( 39 ) = JVS( 250 )
   W( 40 ) = JVS( 251 )
   W( 41 ) = JVS( 252 )
   W( 42 ) = JVS( 253 )
   W( 44 ) = JVS( 254 )
   W( 45 ) = JVS( 255 )
   W( 46 ) = JVS( 256 )
   W( 49 ) = JVS( 257 )
   W( 50 ) = JVS( 258 )
   W( 51 ) = JVS( 259 )
   W( 52 ) = JVS( 260 )
   W( 53 ) = JVS( 261 )
   W( 54 ) = JVS( 262 )
   W( 55 ) = JVS( 263 )
   W( 56 ) = JVS( 264 )
   W( 57 ) = JVS( 265 )
   W( 58 ) = JVS( 266 )
   W( 59 ) = JVS( 267 )
  a = -W( 32 ) / JVS(          150  )
  W( 32 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 151 )
  W( 52 ) = W( 52 ) + a*JVS( 152 )
  W( 53 ) = W( 53 ) + a*JVS( 153 )
  W( 54 ) = W( 54 ) + a*JVS( 154 )
  W( 55 ) = W( 55 ) + a*JVS( 155 )
  W( 59 ) = W( 59 ) + a*JVS( 156 )
  a = -W( 35 ) / JVS(          182  )
  W( 35 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 183 )
  W( 52 ) = W( 52 ) + a*JVS( 184 )
  W( 53 ) = W( 53 ) + a*JVS( 185 )
  W( 55 ) = W( 55 ) + a*JVS( 186 )
  W( 58 ) = W( 58 ) + a*JVS( 187 )
  a = -W( 36 ) / JVS(          189  )
  W( 36 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 190 )
  W( 52 ) = W( 52 ) + a*JVS( 191 )
  W( 53 ) = W( 53 ) + a*JVS( 192 )
  W( 55 ) = W( 55 ) + a*JVS( 193 )
  W( 58 ) = W( 58 ) + a*JVS( 194 )
  a = -W( 37 ) / JVS(          197  )
  W( 37 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 198 )
  W( 51 ) = W( 51 ) + a*JVS( 199 )
  W( 52 ) = W( 52 ) + a*JVS( 200 )
  W( 53 ) = W( 53 ) + a*JVS( 201 )
  W( 54 ) = W( 54 ) + a*JVS( 202 )
  W( 55 ) = W( 55 ) + a*JVS( 203 )
  W( 58 ) = W( 58 ) + a*JVS( 204 )
  a = -W( 39 ) / JVS(          221  )
  W( 39 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 222 )
  W( 51 ) = W( 51 ) + a*JVS( 223 )
  W( 52 ) = W( 52 ) + a*JVS( 224 )
  W( 53 ) = W( 53 ) + a*JVS( 225 )
  W( 54 ) = W( 54 ) + a*JVS( 226 )
  W( 55 ) = W( 55 ) + a*JVS( 227 )
  W( 58 ) = W( 58 ) + a*JVS( 228 )
  W( 59 ) = W( 59 ) + a*JVS( 229 )
  a = -W( 40 ) / JVS(          231  )
  W( 40 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 52 ) = W( 52 ) + a*JVS( 234 )
  W( 53 ) = W( 53 ) + a*JVS( 235 )
  W( 54 ) = W( 54 ) + a*JVS( 236 )
  W( 55 ) = W( 55 ) + a*JVS( 237 )
  W( 58 ) = W( 58 ) + a*JVS( 238 )
  a = -W( 41 ) / JVS(          240  )
  W( 41 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 241 )
  W( 52 ) = W( 52 ) + a*JVS( 242 )
  W( 53 ) = W( 53 ) + a*JVS( 243 )
  W( 55 ) = W( 55 ) + a*JVS( 244 )
  W( 58 ) = W( 58 ) + a*JVS( 245 )
  JVS( 246) = W( 32 )
  JVS( 247) = W( 35 )
  JVS( 248) = W( 36 )
  JVS( 249) = W( 37 )
  JVS( 250) = W( 39 )
  JVS( 251) = W( 40 )
  JVS( 252) = W( 41 )
  JVS( 253) = W( 42 )
  JVS( 254) = W( 44 )
  JVS( 255) = W( 45 )
  JVS( 256) = W( 46 )
  JVS( 257) = W( 49 )
  JVS( 258) = W( 50 )
  JVS( 259) = W( 51 )
  JVS( 260) = W( 52 )
  JVS( 261) = W( 53 )
  JVS( 262) = W( 54 )
  JVS( 263) = W( 55 )
  JVS( 264) = W( 56 )
  JVS( 265) = W( 57 )
  JVS( 266) = W( 58 )
  JVS( 267) = W( 59 )
  IF ( ABS(  JVS( 284 )) < TINY(a) ) THEN
         IER = 43                                      
         RETURN
  END IF
   W( 16 ) = JVS( 268 )
   W( 21 ) = JVS( 269 )
   W( 23 ) = JVS( 270 )
   W( 24 ) = JVS( 271 )
   W( 27 ) = JVS( 272 )
   W( 28 ) = JVS( 273 )
   W( 29 ) = JVS( 274 )
   W( 31 ) = JVS( 275 )
   W( 32 ) = JVS( 276 )
   W( 35 ) = JVS( 277 )
   W( 36 ) = JVS( 278 )
   W( 37 ) = JVS( 279 )
   W( 38 ) = JVS( 280 )
   W( 39 ) = JVS( 281 )
   W( 40 ) = JVS( 282 )
   W( 41 ) = JVS( 283 )
   W( 43 ) = JVS( 284 )
   W( 44 ) = JVS( 285 )
   W( 45 ) = JVS( 286 )
   W( 46 ) = JVS( 287 )
   W( 49 ) = JVS( 288 )
   W( 50 ) = JVS( 289 )
   W( 51 ) = JVS( 290 )
   W( 52 ) = JVS( 291 )
   W( 53 ) = JVS( 292 )
   W( 54 ) = JVS( 293 )
   W( 55 ) = JVS( 294 )
   W( 56 ) = JVS( 295 )
   W( 57 ) = JVS( 296 )
   W( 58 ) = JVS( 297 )
   W( 59 ) = JVS( 298 )
  a = -W( 16 ) / JVS(           56  )
  W( 16 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 57 )
  a = -W( 21 ) / JVS(           76  )
  W( 21 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 77 )
  W( 52 ) = W( 52 ) + a*JVS( 78 )
  W( 53 ) = W( 53 ) + a*JVS( 79 )
  a = -W( 23 ) / JVS(           85  )
  W( 23 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 86 )
  W( 55 ) = W( 55 ) + a*JVS( 87 )
  W( 59 ) = W( 59 ) + a*JVS( 88 )
  a = -W( 24 ) / JVS(           89  )
  W( 24 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 90 )
  W( 53 ) = W( 53 ) + a*JVS( 91 )
  W( 54 ) = W( 54 ) + a*JVS( 92 )
  a = -W( 27 ) / JVS(          121  )
  W( 27 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 122 )
  W( 53 ) = W( 53 ) + a*JVS( 123 )
  W( 54 ) = W( 54 ) + a*JVS( 124 )
  a = -W( 28 ) / JVS(          125  )
  W( 28 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 126 )
  W( 53 ) = W( 53 ) + a*JVS( 127 )
  W( 54 ) = W( 54 ) + a*JVS( 128 )
  a = -W( 29 ) / JVS(          129  )
  W( 29 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 130 )
  W( 53 ) = W( 53 ) + a*JVS( 131 )
  W( 54 ) = W( 54 ) + a*JVS( 132 )
  a = -W( 31 ) / JVS(          141  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 142 )
  W( 45 ) = W( 45 ) + a*JVS( 143 )
  W( 52 ) = W( 52 ) + a*JVS( 144 )
  W( 53 ) = W( 53 ) + a*JVS( 145 )
  W( 54 ) = W( 54 ) + a*JVS( 146 )
  W( 55 ) = W( 55 ) + a*JVS( 147 )
  W( 58 ) = W( 58 ) + a*JVS( 148 )
  a = -W( 32 ) / JVS(          150  )
  W( 32 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 151 )
  W( 52 ) = W( 52 ) + a*JVS( 152 )
  W( 53 ) = W( 53 ) + a*JVS( 153 )
  W( 54 ) = W( 54 ) + a*JVS( 154 )
  W( 55 ) = W( 55 ) + a*JVS( 155 )
  W( 59 ) = W( 59 ) + a*JVS( 156 )
  a = -W( 35 ) / JVS(          182  )
  W( 35 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 183 )
  W( 52 ) = W( 52 ) + a*JVS( 184 )
  W( 53 ) = W( 53 ) + a*JVS( 185 )
  W( 55 ) = W( 55 ) + a*JVS( 186 )
  W( 58 ) = W( 58 ) + a*JVS( 187 )
  a = -W( 36 ) / JVS(          189  )
  W( 36 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 190 )
  W( 52 ) = W( 52 ) + a*JVS( 191 )
  W( 53 ) = W( 53 ) + a*JVS( 192 )
  W( 55 ) = W( 55 ) + a*JVS( 193 )
  W( 58 ) = W( 58 ) + a*JVS( 194 )
  a = -W( 37 ) / JVS(          197  )
  W( 37 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 198 )
  W( 51 ) = W( 51 ) + a*JVS( 199 )
  W( 52 ) = W( 52 ) + a*JVS( 200 )
  W( 53 ) = W( 53 ) + a*JVS( 201 )
  W( 54 ) = W( 54 ) + a*JVS( 202 )
  W( 55 ) = W( 55 ) + a*JVS( 203 )
  W( 58 ) = W( 58 ) + a*JVS( 204 )
  a = -W( 38 ) / JVS(          209  )
  W( 38 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 210 )
  W( 51 ) = W( 51 ) + a*JVS( 211 )
  W( 52 ) = W( 52 ) + a*JVS( 212 )
  W( 53 ) = W( 53 ) + a*JVS( 213 )
  W( 54 ) = W( 54 ) + a*JVS( 214 )
  W( 55 ) = W( 55 ) + a*JVS( 215 )
  W( 58 ) = W( 58 ) + a*JVS( 216 )
  a = -W( 39 ) / JVS(          221  )
  W( 39 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 222 )
  W( 51 ) = W( 51 ) + a*JVS( 223 )
  W( 52 ) = W( 52 ) + a*JVS( 224 )
  W( 53 ) = W( 53 ) + a*JVS( 225 )
  W( 54 ) = W( 54 ) + a*JVS( 226 )
  W( 55 ) = W( 55 ) + a*JVS( 227 )
  W( 58 ) = W( 58 ) + a*JVS( 228 )
  W( 59 ) = W( 59 ) + a*JVS( 229 )
  a = -W( 40 ) / JVS(          231  )
  W( 40 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 52 ) = W( 52 ) + a*JVS( 234 )
  W( 53 ) = W( 53 ) + a*JVS( 235 )
  W( 54 ) = W( 54 ) + a*JVS( 236 )
  W( 55 ) = W( 55 ) + a*JVS( 237 )
  W( 58 ) = W( 58 ) + a*JVS( 238 )
  a = -W( 41 ) / JVS(          240  )
  W( 41 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 241 )
  W( 52 ) = W( 52 ) + a*JVS( 242 )
  W( 53 ) = W( 53 ) + a*JVS( 243 )
  W( 55 ) = W( 55 ) + a*JVS( 244 )
  W( 58 ) = W( 58 ) + a*JVS( 245 )
  JVS( 268) = W( 16 )
  JVS( 269) = W( 21 )
  JVS( 270) = W( 23 )
  JVS( 271) = W( 24 )
  JVS( 272) = W( 27 )
  JVS( 273) = W( 28 )
  JVS( 274) = W( 29 )
  JVS( 275) = W( 31 )
  JVS( 276) = W( 32 )
  JVS( 277) = W( 35 )
  JVS( 278) = W( 36 )
  JVS( 279) = W( 37 )
  JVS( 280) = W( 38 )
  JVS( 281) = W( 39 )
  JVS( 282) = W( 40 )
  JVS( 283) = W( 41 )
  JVS( 284) = W( 43 )
  JVS( 285) = W( 44 )
  JVS( 286) = W( 45 )
  JVS( 287) = W( 46 )
  JVS( 288) = W( 49 )
  JVS( 289) = W( 50 )
  JVS( 290) = W( 51 )
  JVS( 291) = W( 52 )
  JVS( 292) = W( 53 )
  JVS( 293) = W( 54 )
  JVS( 294) = W( 55 )
  JVS( 295) = W( 56 )
  JVS( 296) = W( 57 )
  JVS( 297) = W( 58 )
  JVS( 298) = W( 59 )
  IF ( ABS(  JVS( 300 )) < TINY(a) ) THEN
         IER = 44                                      
         RETURN
  END IF
   W( 9 ) = JVS( 299 )
   W( 44 ) = JVS( 300 )
   W( 51 ) = JVS( 301 )
   W( 52 ) = JVS( 302 )
   W( 53 ) = JVS( 303 )
   W( 55 ) = JVS( 304 )
   W( 58 ) = JVS( 305 )
  a = -W( 9 ) / JVS(           38  )
  W( 9 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 39 )
  JVS( 299) = W( 9 )
  JVS( 300) = W( 44 )
  JVS( 301) = W( 51 )
  JVS( 302) = W( 52 )
  JVS( 303) = W( 53 )
  JVS( 304) = W( 55 )
  JVS( 305) = W( 58 )
  IF ( ABS(  JVS( 311 )) < TINY(a) ) THEN
         IER = 45                                      
         RETURN
  END IF
   W( 12 ) = JVS( 306 )
   W( 22 ) = JVS( 307 )
   W( 30 ) = JVS( 308 )
   W( 35 ) = JVS( 309 )
   W( 36 ) = JVS( 310 )
   W( 45 ) = JVS( 311 )
   W( 51 ) = JVS( 312 )
   W( 52 ) = JVS( 313 )
   W( 53 ) = JVS( 314 )
   W( 54 ) = JVS( 315 )
   W( 55 ) = JVS( 316 )
   W( 58 ) = JVS( 317 )
   W( 59 ) = JVS( 318 )
  a = -W( 12 ) / JVS(           44  )
  W( 12 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 45 )
  W( 59 ) = W( 59 ) + a*JVS( 46 )
  a = -W( 22 ) / JVS(           82  )
  W( 22 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 83 )
  W( 54 ) = W( 54 ) + a*JVS( 84 )
  a = -W( 30 ) / JVS(          133  )
  W( 30 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 134 )
  W( 36 ) = W( 36 ) + a*JVS( 135 )
  W( 52 ) = W( 52 ) + a*JVS( 136 )
  W( 53 ) = W( 53 ) + a*JVS( 137 )
  W( 54 ) = W( 54 ) + a*JVS( 138 )
  W( 55 ) = W( 55 ) + a*JVS( 139 )
  W( 58 ) = W( 58 ) + a*JVS( 140 )
  a = -W( 35 ) / JVS(          182  )
  W( 35 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 183 )
  W( 52 ) = W( 52 ) + a*JVS( 184 )
  W( 53 ) = W( 53 ) + a*JVS( 185 )
  W( 55 ) = W( 55 ) + a*JVS( 186 )
  W( 58 ) = W( 58 ) + a*JVS( 187 )
  a = -W( 36 ) / JVS(          189  )
  W( 36 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 190 )
  W( 52 ) = W( 52 ) + a*JVS( 191 )
  W( 53 ) = W( 53 ) + a*JVS( 192 )
  W( 55 ) = W( 55 ) + a*JVS( 193 )
  W( 58 ) = W( 58 ) + a*JVS( 194 )
  JVS( 306) = W( 12 )
  JVS( 307) = W( 22 )
  JVS( 308) = W( 30 )
  JVS( 309) = W( 35 )
  JVS( 310) = W( 36 )
  JVS( 311) = W( 45 )
  JVS( 312) = W( 51 )
  JVS( 313) = W( 52 )
  JVS( 314) = W( 53 )
  JVS( 315) = W( 54 )
  JVS( 316) = W( 55 )
  JVS( 317) = W( 58 )
  JVS( 318) = W( 59 )
  IF ( ABS(  JVS( 324 )) < TINY(a) ) THEN
         IER = 46                                      
         RETURN
  END IF
   W( 19 ) = JVS( 319 )
   W( 24 ) = JVS( 320 )
   W( 27 ) = JVS( 321 )
   W( 28 ) = JVS( 322 )
   W( 29 ) = JVS( 323 )
   W( 46 ) = JVS( 324 )
   W( 51 ) = JVS( 325 )
   W( 53 ) = JVS( 326 )
   W( 54 ) = JVS( 327 )
   W( 58 ) = JVS( 328 )
   W( 59 ) = JVS( 329 )
  a = -W( 19 ) / JVS(           68  )
  W( 19 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 69 )
  W( 54 ) = W( 54 ) + a*JVS( 70 )
  W( 59 ) = W( 59 ) + a*JVS( 71 )
  a = -W( 24 ) / JVS(           89  )
  W( 24 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 90 )
  W( 53 ) = W( 53 ) + a*JVS( 91 )
  W( 54 ) = W( 54 ) + a*JVS( 92 )
  a = -W( 27 ) / JVS(          121  )
  W( 27 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 122 )
  W( 53 ) = W( 53 ) + a*JVS( 123 )
  W( 54 ) = W( 54 ) + a*JVS( 124 )
  a = -W( 28 ) / JVS(          125  )
  W( 28 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 126 )
  W( 53 ) = W( 53 ) + a*JVS( 127 )
  W( 54 ) = W( 54 ) + a*JVS( 128 )
  a = -W( 29 ) / JVS(          129  )
  W( 29 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 130 )
  W( 53 ) = W( 53 ) + a*JVS( 131 )
  W( 54 ) = W( 54 ) + a*JVS( 132 )
  JVS( 319) = W( 19 )
  JVS( 320) = W( 24 )
  JVS( 321) = W( 27 )
  JVS( 322) = W( 28 )
  JVS( 323) = W( 29 )
  JVS( 324) = W( 46 )
  JVS( 325) = W( 51 )
  JVS( 326) = W( 53 )
  JVS( 327) = W( 54 )
  JVS( 328) = W( 58 )
  JVS( 329) = W( 59 )
  IF ( ABS(  JVS( 335 )) < TINY(a) ) THEN
         IER = 47                                      
         RETURN
  END IF
   W( 32 ) = JVS( 330 )
   W( 38 ) = JVS( 331 )
   W( 41 ) = JVS( 332 )
   W( 44 ) = JVS( 333 )
   W( 46 ) = JVS( 334 )
   W( 47 ) = JVS( 335 )
   W( 51 ) = JVS( 336 )
   W( 52 ) = JVS( 337 )
   W( 53 ) = JVS( 338 )
   W( 54 ) = JVS( 339 )
   W( 55 ) = JVS( 340 )
   W( 56 ) = JVS( 341 )
   W( 58 ) = JVS( 342 )
   W( 59 ) = JVS( 343 )
  a = -W( 32 ) / JVS(          150  )
  W( 32 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 151 )
  W( 52 ) = W( 52 ) + a*JVS( 152 )
  W( 53 ) = W( 53 ) + a*JVS( 153 )
  W( 54 ) = W( 54 ) + a*JVS( 154 )
  W( 55 ) = W( 55 ) + a*JVS( 155 )
  W( 59 ) = W( 59 ) + a*JVS( 156 )
  a = -W( 38 ) / JVS(          209  )
  W( 38 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 210 )
  W( 51 ) = W( 51 ) + a*JVS( 211 )
  W( 52 ) = W( 52 ) + a*JVS( 212 )
  W( 53 ) = W( 53 ) + a*JVS( 213 )
  W( 54 ) = W( 54 ) + a*JVS( 214 )
  W( 55 ) = W( 55 ) + a*JVS( 215 )
  W( 58 ) = W( 58 ) + a*JVS( 216 )
  a = -W( 41 ) / JVS(          240  )
  W( 41 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 241 )
  W( 52 ) = W( 52 ) + a*JVS( 242 )
  W( 53 ) = W( 53 ) + a*JVS( 243 )
  W( 55 ) = W( 55 ) + a*JVS( 244 )
  W( 58 ) = W( 58 ) + a*JVS( 245 )
  a = -W( 44 ) / JVS(          300  )
  W( 44 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 301 )
  W( 52 ) = W( 52 ) + a*JVS( 302 )
  W( 53 ) = W( 53 ) + a*JVS( 303 )
  W( 55 ) = W( 55 ) + a*JVS( 304 )
  W( 58 ) = W( 58 ) + a*JVS( 305 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 325 )
  W( 53 ) = W( 53 ) + a*JVS( 326 )
  W( 54 ) = W( 54 ) + a*JVS( 327 )
  W( 58 ) = W( 58 ) + a*JVS( 328 )
  W( 59 ) = W( 59 ) + a*JVS( 329 )
  JVS( 330) = W( 32 )
  JVS( 331) = W( 38 )
  JVS( 332) = W( 41 )
  JVS( 333) = W( 44 )
  JVS( 334) = W( 46 )
  JVS( 335) = W( 47 )
  JVS( 336) = W( 51 )
  JVS( 337) = W( 52 )
  JVS( 338) = W( 53 )
  JVS( 339) = W( 54 )
  JVS( 340) = W( 55 )
  JVS( 341) = W( 56 )
  JVS( 342) = W( 58 )
  JVS( 343) = W( 59 )
  IF ( ABS(  JVS( 357 )) < TINY(a) ) THEN
         IER = 48                                      
         RETURN
  END IF
   W( 16 ) = JVS( 344 )
   W( 27 ) = JVS( 345 )
   W( 28 ) = JVS( 346 )
   W( 29 ) = JVS( 347 )
   W( 37 ) = JVS( 348 )
   W( 38 ) = JVS( 349 )
   W( 40 ) = JVS( 350 )
   W( 41 ) = JVS( 351 )
   W( 42 ) = JVS( 352 )
   W( 44 ) = JVS( 353 )
   W( 45 ) = JVS( 354 )
   W( 46 ) = JVS( 355 )
   W( 47 ) = JVS( 356 )
   W( 48 ) = JVS( 357 )
   W( 49 ) = JVS( 358 )
   W( 50 ) = JVS( 359 )
   W( 51 ) = JVS( 360 )
   W( 52 ) = JVS( 361 )
   W( 53 ) = JVS( 362 )
   W( 54 ) = JVS( 363 )
   W( 55 ) = JVS( 364 )
   W( 56 ) = JVS( 365 )
   W( 57 ) = JVS( 366 )
   W( 58 ) = JVS( 367 )
   W( 59 ) = JVS( 368 )
  a = -W( 16 ) / JVS(           56  )
  W( 16 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 57 )
  a = -W( 27 ) / JVS(          121  )
  W( 27 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 122 )
  W( 53 ) = W( 53 ) + a*JVS( 123 )
  W( 54 ) = W( 54 ) + a*JVS( 124 )
  a = -W( 28 ) / JVS(          125  )
  W( 28 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 126 )
  W( 53 ) = W( 53 ) + a*JVS( 127 )
  W( 54 ) = W( 54 ) + a*JVS( 128 )
  a = -W( 29 ) / JVS(          129  )
  W( 29 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 130 )
  W( 53 ) = W( 53 ) + a*JVS( 131 )
  W( 54 ) = W( 54 ) + a*JVS( 132 )
  a = -W( 37 ) / JVS(          197  )
  W( 37 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 198 )
  W( 51 ) = W( 51 ) + a*JVS( 199 )
  W( 52 ) = W( 52 ) + a*JVS( 200 )
  W( 53 ) = W( 53 ) + a*JVS( 201 )
  W( 54 ) = W( 54 ) + a*JVS( 202 )
  W( 55 ) = W( 55 ) + a*JVS( 203 )
  W( 58 ) = W( 58 ) + a*JVS( 204 )
  a = -W( 38 ) / JVS(          209  )
  W( 38 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 210 )
  W( 51 ) = W( 51 ) + a*JVS( 211 )
  W( 52 ) = W( 52 ) + a*JVS( 212 )
  W( 53 ) = W( 53 ) + a*JVS( 213 )
  W( 54 ) = W( 54 ) + a*JVS( 214 )
  W( 55 ) = W( 55 ) + a*JVS( 215 )
  W( 58 ) = W( 58 ) + a*JVS( 216 )
  a = -W( 40 ) / JVS(          231  )
  W( 40 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 52 ) = W( 52 ) + a*JVS( 234 )
  W( 53 ) = W( 53 ) + a*JVS( 235 )
  W( 54 ) = W( 54 ) + a*JVS( 236 )
  W( 55 ) = W( 55 ) + a*JVS( 237 )
  W( 58 ) = W( 58 ) + a*JVS( 238 )
  a = -W( 41 ) / JVS(          240  )
  W( 41 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 241 )
  W( 52 ) = W( 52 ) + a*JVS( 242 )
  W( 53 ) = W( 53 ) + a*JVS( 243 )
  W( 55 ) = W( 55 ) + a*JVS( 244 )
  W( 58 ) = W( 58 ) + a*JVS( 245 )
  a = -W( 42 ) / JVS(          253  )
  W( 42 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 254 )
  W( 45 ) = W( 45 ) + a*JVS( 255 )
  W( 46 ) = W( 46 ) + a*JVS( 256 )
  W( 49 ) = W( 49 ) + a*JVS( 257 )
  W( 50 ) = W( 50 ) + a*JVS( 258 )
  W( 51 ) = W( 51 ) + a*JVS( 259 )
  W( 52 ) = W( 52 ) + a*JVS( 260 )
  W( 53 ) = W( 53 ) + a*JVS( 261 )
  W( 54 ) = W( 54 ) + a*JVS( 262 )
  W( 55 ) = W( 55 ) + a*JVS( 263 )
  W( 56 ) = W( 56 ) + a*JVS( 264 )
  W( 57 ) = W( 57 ) + a*JVS( 265 )
  W( 58 ) = W( 58 ) + a*JVS( 266 )
  W( 59 ) = W( 59 ) + a*JVS( 267 )
  a = -W( 44 ) / JVS(          300  )
  W( 44 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 301 )
  W( 52 ) = W( 52 ) + a*JVS( 302 )
  W( 53 ) = W( 53 ) + a*JVS( 303 )
  W( 55 ) = W( 55 ) + a*JVS( 304 )
  W( 58 ) = W( 58 ) + a*JVS( 305 )
  a = -W( 45 ) / JVS(          311  )
  W( 45 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 312 )
  W( 52 ) = W( 52 ) + a*JVS( 313 )
  W( 53 ) = W( 53 ) + a*JVS( 314 )
  W( 54 ) = W( 54 ) + a*JVS( 315 )
  W( 55 ) = W( 55 ) + a*JVS( 316 )
  W( 58 ) = W( 58 ) + a*JVS( 317 )
  W( 59 ) = W( 59 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 325 )
  W( 53 ) = W( 53 ) + a*JVS( 326 )
  W( 54 ) = W( 54 ) + a*JVS( 327 )
  W( 58 ) = W( 58 ) + a*JVS( 328 )
  W( 59 ) = W( 59 ) + a*JVS( 329 )
  a = -W( 47 ) / JVS(          335  )
  W( 47 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 336 )
  W( 52 ) = W( 52 ) + a*JVS( 337 )
  W( 53 ) = W( 53 ) + a*JVS( 338 )
  W( 54 ) = W( 54 ) + a*JVS( 339 )
  W( 55 ) = W( 55 ) + a*JVS( 340 )
  W( 56 ) = W( 56 ) + a*JVS( 341 )
  W( 58 ) = W( 58 ) + a*JVS( 342 )
  W( 59 ) = W( 59 ) + a*JVS( 343 )
  JVS( 344) = W( 16 )
  JVS( 345) = W( 27 )
  JVS( 346) = W( 28 )
  JVS( 347) = W( 29 )
  JVS( 348) = W( 37 )
  JVS( 349) = W( 38 )
  JVS( 350) = W( 40 )
  JVS( 351) = W( 41 )
  JVS( 352) = W( 42 )
  JVS( 353) = W( 44 )
  JVS( 354) = W( 45 )
  JVS( 355) = W( 46 )
  JVS( 356) = W( 47 )
  JVS( 357) = W( 48 )
  JVS( 358) = W( 49 )
  JVS( 359) = W( 50 )
  JVS( 360) = W( 51 )
  JVS( 361) = W( 52 )
  JVS( 362) = W( 53 )
  JVS( 363) = W( 54 )
  JVS( 364) = W( 55 )
  JVS( 365) = W( 56 )
  JVS( 366) = W( 57 )
  JVS( 367) = W( 58 )
  JVS( 368) = W( 59 )
  IF ( ABS(  JVS( 371 )) < TINY(a) ) THEN
         IER = 49                                      
         RETURN
  END IF
   W( 29 ) = JVS( 369 )
   W( 46 ) = JVS( 370 )
   W( 49 ) = JVS( 371 )
   W( 51 ) = JVS( 372 )
   W( 52 ) = JVS( 373 )
   W( 53 ) = JVS( 374 )
   W( 54 ) = JVS( 375 )
   W( 55 ) = JVS( 376 )
   W( 58 ) = JVS( 377 )
   W( 59 ) = JVS( 378 )
  a = -W( 29 ) / JVS(          129  )
  W( 29 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 130 )
  W( 53 ) = W( 53 ) + a*JVS( 131 )
  W( 54 ) = W( 54 ) + a*JVS( 132 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 325 )
  W( 53 ) = W( 53 ) + a*JVS( 326 )
  W( 54 ) = W( 54 ) + a*JVS( 327 )
  W( 58 ) = W( 58 ) + a*JVS( 328 )
  W( 59 ) = W( 59 ) + a*JVS( 329 )
  JVS( 369) = W( 29 )
  JVS( 370) = W( 46 )
  JVS( 371) = W( 49 )
  JVS( 372) = W( 51 )
  JVS( 373) = W( 52 )
  JVS( 374) = W( 53 )
  JVS( 375) = W( 54 )
  JVS( 376) = W( 55 )
  JVS( 377) = W( 58 )
  JVS( 378) = W( 59 )
  IF ( ABS(  JVS( 385 )) < TINY(a) ) THEN
         IER = 50                                      
         RETURN
  END IF
   W( 33 ) = JVS( 379 )
   W( 41 ) = JVS( 380 )
   W( 44 ) = JVS( 381 )
   W( 46 ) = JVS( 382 )
   W( 47 ) = JVS( 383 )
   W( 49 ) = JVS( 384 )
   W( 50 ) = JVS( 385 )
   W( 51 ) = JVS( 386 )
   W( 52 ) = JVS( 387 )
   W( 53 ) = JVS( 388 )
   W( 54 ) = JVS( 389 )
   W( 55 ) = JVS( 390 )
   W( 56 ) = JVS( 391 )
   W( 58 ) = JVS( 392 )
   W( 59 ) = JVS( 393 )
  a = -W( 33 ) / JVS(          159  )
  W( 33 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 160 )
  W( 44 ) = W( 44 ) + a*JVS( 161 )
  W( 46 ) = W( 46 ) + a*JVS( 162 )
  W( 47 ) = W( 47 ) + a*JVS( 163 )
  W( 49 ) = W( 49 ) + a*JVS( 164 )
  W( 52 ) = W( 52 ) + a*JVS( 165 )
  W( 53 ) = W( 53 ) + a*JVS( 166 )
  W( 54 ) = W( 54 ) + a*JVS( 167 )
  W( 55 ) = W( 55 ) + a*JVS( 168 )
  W( 56 ) = W( 56 ) + a*JVS( 169 )
  W( 58 ) = W( 58 ) + a*JVS( 170 )
  a = -W( 41 ) / JVS(          240  )
  W( 41 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 241 )
  W( 52 ) = W( 52 ) + a*JVS( 242 )
  W( 53 ) = W( 53 ) + a*JVS( 243 )
  W( 55 ) = W( 55 ) + a*JVS( 244 )
  W( 58 ) = W( 58 ) + a*JVS( 245 )
  a = -W( 44 ) / JVS(          300  )
  W( 44 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 301 )
  W( 52 ) = W( 52 ) + a*JVS( 302 )
  W( 53 ) = W( 53 ) + a*JVS( 303 )
  W( 55 ) = W( 55 ) + a*JVS( 304 )
  W( 58 ) = W( 58 ) + a*JVS( 305 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 325 )
  W( 53 ) = W( 53 ) + a*JVS( 326 )
  W( 54 ) = W( 54 ) + a*JVS( 327 )
  W( 58 ) = W( 58 ) + a*JVS( 328 )
  W( 59 ) = W( 59 ) + a*JVS( 329 )
  a = -W( 47 ) / JVS(          335  )
  W( 47 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 336 )
  W( 52 ) = W( 52 ) + a*JVS( 337 )
  W( 53 ) = W( 53 ) + a*JVS( 338 )
  W( 54 ) = W( 54 ) + a*JVS( 339 )
  W( 55 ) = W( 55 ) + a*JVS( 340 )
  W( 56 ) = W( 56 ) + a*JVS( 341 )
  W( 58 ) = W( 58 ) + a*JVS( 342 )
  W( 59 ) = W( 59 ) + a*JVS( 343 )
  a = -W( 49 ) / JVS(          371  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 372 )
  W( 52 ) = W( 52 ) + a*JVS( 373 )
  W( 53 ) = W( 53 ) + a*JVS( 374 )
  W( 54 ) = W( 54 ) + a*JVS( 375 )
  W( 55 ) = W( 55 ) + a*JVS( 376 )
  W( 58 ) = W( 58 ) + a*JVS( 377 )
  W( 59 ) = W( 59 ) + a*JVS( 378 )
  JVS( 379) = W( 33 )
  JVS( 380) = W( 41 )
  JVS( 381) = W( 44 )
  JVS( 382) = W( 46 )
  JVS( 383) = W( 47 )
  JVS( 384) = W( 49 )
  JVS( 385) = W( 50 )
  JVS( 386) = W( 51 )
  JVS( 387) = W( 52 )
  JVS( 388) = W( 53 )
  JVS( 389) = W( 54 )
  JVS( 390) = W( 55 )
  JVS( 391) = W( 56 )
  JVS( 392) = W( 58 )
  JVS( 393) = W( 59 )
  IF ( ABS(  JVS( 427 )) < TINY(a) ) THEN
         IER = 51                                      
         RETURN
  END IF
   W( 5 ) = JVS( 394 )
   W( 10 ) = JVS( 395 )
   W( 11 ) = JVS( 396 )
   W( 14 ) = JVS( 397 )
   W( 16 ) = JVS( 398 )
   W( 20 ) = JVS( 399 )
   W( 21 ) = JVS( 400 )
   W( 22 ) = JVS( 401 )
   W( 24 ) = JVS( 402 )
   W( 26 ) = JVS( 403 )
   W( 27 ) = JVS( 404 )
   W( 28 ) = JVS( 405 )
   W( 29 ) = JVS( 406 )
   W( 30 ) = JVS( 407 )
   W( 31 ) = JVS( 408 )
   W( 32 ) = JVS( 409 )
   W( 34 ) = JVS( 410 )
   W( 35 ) = JVS( 411 )
   W( 36 ) = JVS( 412 )
   W( 37 ) = JVS( 413 )
   W( 38 ) = JVS( 414 )
   W( 39 ) = JVS( 415 )
   W( 40 ) = JVS( 416 )
   W( 41 ) = JVS( 417 )
   W( 42 ) = JVS( 418 )
   W( 43 ) = JVS( 419 )
   W( 44 ) = JVS( 420 )
   W( 45 ) = JVS( 421 )
   W( 46 ) = JVS( 422 )
   W( 47 ) = JVS( 423 )
   W( 48 ) = JVS( 424 )
   W( 49 ) = JVS( 425 )
   W( 50 ) = JVS( 426 )
   W( 51 ) = JVS( 427 )
   W( 52 ) = JVS( 428 )
   W( 53 ) = JVS( 429 )
   W( 54 ) = JVS( 430 )
   W( 55 ) = JVS( 431 )
   W( 56 ) = JVS( 432 )
   W( 57 ) = JVS( 433 )
   W( 58 ) = JVS( 434 )
   W( 59 ) = JVS( 435 )
  a = -W( 5 ) / JVS(           30  )
  W( 5 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 31 )
  a = -W( 10 ) / JVS(           40  )
  W( 10 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 41 )
  a = -W( 11 ) / JVS(           42  )
  W( 11 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 43 )
  a = -W( 14 ) / JVS(           50  )
  W( 14 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 51 )
  W( 53 ) = W( 53 ) + a*JVS( 52 )
  a = -W( 16 ) / JVS(           56  )
  W( 16 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 57 )
  a = -W( 20 ) / JVS(           72  )
  W( 20 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 73 )
  W( 53 ) = W( 53 ) + a*JVS( 74 )
  W( 59 ) = W( 59 ) + a*JVS( 75 )
  a = -W( 21 ) / JVS(           76  )
  W( 21 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 77 )
  W( 52 ) = W( 52 ) + a*JVS( 78 )
  W( 53 ) = W( 53 ) + a*JVS( 79 )
  a = -W( 22 ) / JVS(           82  )
  W( 22 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 83 )
  W( 54 ) = W( 54 ) + a*JVS( 84 )
  a = -W( 24 ) / JVS(           89  )
  W( 24 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 90 )
  W( 53 ) = W( 53 ) + a*JVS( 91 )
  W( 54 ) = W( 54 ) + a*JVS( 92 )
  a = -W( 26 ) / JVS(          106  )
  W( 26 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 107 )
  W( 28 ) = W( 28 ) + a*JVS( 108 )
  W( 29 ) = W( 29 ) + a*JVS( 109 )
  W( 31 ) = W( 31 ) + a*JVS( 110 )
  W( 34 ) = W( 34 ) + a*JVS( 111 )
  W( 43 ) = W( 43 ) + a*JVS( 112 )
  W( 45 ) = W( 45 ) + a*JVS( 113 )
  W( 46 ) = W( 46 ) + a*JVS( 114 )
  W( 48 ) = W( 48 ) + a*JVS( 115 )
  W( 52 ) = W( 52 ) + a*JVS( 116 )
  W( 53 ) = W( 53 ) + a*JVS( 117 )
  W( 54 ) = W( 54 ) + a*JVS( 118 )
  W( 55 ) = W( 55 ) + a*JVS( 119 )
  W( 58 ) = W( 58 ) + a*JVS( 120 )
  a = -W( 27 ) / JVS(          121  )
  W( 27 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 122 )
  W( 53 ) = W( 53 ) + a*JVS( 123 )
  W( 54 ) = W( 54 ) + a*JVS( 124 )
  a = -W( 28 ) / JVS(          125  )
  W( 28 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 126 )
  W( 53 ) = W( 53 ) + a*JVS( 127 )
  W( 54 ) = W( 54 ) + a*JVS( 128 )
  a = -W( 29 ) / JVS(          129  )
  W( 29 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 130 )
  W( 53 ) = W( 53 ) + a*JVS( 131 )
  W( 54 ) = W( 54 ) + a*JVS( 132 )
  a = -W( 30 ) / JVS(          133  )
  W( 30 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 134 )
  W( 36 ) = W( 36 ) + a*JVS( 135 )
  W( 52 ) = W( 52 ) + a*JVS( 136 )
  W( 53 ) = W( 53 ) + a*JVS( 137 )
  W( 54 ) = W( 54 ) + a*JVS( 138 )
  W( 55 ) = W( 55 ) + a*JVS( 139 )
  W( 58 ) = W( 58 ) + a*JVS( 140 )
  a = -W( 31 ) / JVS(          141  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 142 )
  W( 45 ) = W( 45 ) + a*JVS( 143 )
  W( 52 ) = W( 52 ) + a*JVS( 144 )
  W( 53 ) = W( 53 ) + a*JVS( 145 )
  W( 54 ) = W( 54 ) + a*JVS( 146 )
  W( 55 ) = W( 55 ) + a*JVS( 147 )
  W( 58 ) = W( 58 ) + a*JVS( 148 )
  a = -W( 32 ) / JVS(          150  )
  W( 32 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 151 )
  W( 52 ) = W( 52 ) + a*JVS( 152 )
  W( 53 ) = W( 53 ) + a*JVS( 153 )
  W( 54 ) = W( 54 ) + a*JVS( 154 )
  W( 55 ) = W( 55 ) + a*JVS( 155 )
  W( 59 ) = W( 59 ) + a*JVS( 156 )
  a = -W( 34 ) / JVS(          171  )
  W( 34 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 172 )
  W( 36 ) = W( 36 ) + a*JVS( 173 )
  W( 45 ) = W( 45 ) + a*JVS( 174 )
  W( 50 ) = W( 50 ) + a*JVS( 175 )
  W( 52 ) = W( 52 ) + a*JVS( 176 )
  W( 53 ) = W( 53 ) + a*JVS( 177 )
  W( 54 ) = W( 54 ) + a*JVS( 178 )
  W( 55 ) = W( 55 ) + a*JVS( 179 )
  W( 58 ) = W( 58 ) + a*JVS( 180 )
  a = -W( 35 ) / JVS(          182  )
  W( 35 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 183 )
  W( 52 ) = W( 52 ) + a*JVS( 184 )
  W( 53 ) = W( 53 ) + a*JVS( 185 )
  W( 55 ) = W( 55 ) + a*JVS( 186 )
  W( 58 ) = W( 58 ) + a*JVS( 187 )
  a = -W( 36 ) / JVS(          189  )
  W( 36 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 190 )
  W( 52 ) = W( 52 ) + a*JVS( 191 )
  W( 53 ) = W( 53 ) + a*JVS( 192 )
  W( 55 ) = W( 55 ) + a*JVS( 193 )
  W( 58 ) = W( 58 ) + a*JVS( 194 )
  a = -W( 37 ) / JVS(          197  )
  W( 37 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 198 )
  W( 51 ) = W( 51 ) + a*JVS( 199 )
  W( 52 ) = W( 52 ) + a*JVS( 200 )
  W( 53 ) = W( 53 ) + a*JVS( 201 )
  W( 54 ) = W( 54 ) + a*JVS( 202 )
  W( 55 ) = W( 55 ) + a*JVS( 203 )
  W( 58 ) = W( 58 ) + a*JVS( 204 )
  a = -W( 38 ) / JVS(          209  )
  W( 38 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 210 )
  W( 51 ) = W( 51 ) + a*JVS( 211 )
  W( 52 ) = W( 52 ) + a*JVS( 212 )
  W( 53 ) = W( 53 ) + a*JVS( 213 )
  W( 54 ) = W( 54 ) + a*JVS( 214 )
  W( 55 ) = W( 55 ) + a*JVS( 215 )
  W( 58 ) = W( 58 ) + a*JVS( 216 )
  a = -W( 39 ) / JVS(          221  )
  W( 39 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 222 )
  W( 51 ) = W( 51 ) + a*JVS( 223 )
  W( 52 ) = W( 52 ) + a*JVS( 224 )
  W( 53 ) = W( 53 ) + a*JVS( 225 )
  W( 54 ) = W( 54 ) + a*JVS( 226 )
  W( 55 ) = W( 55 ) + a*JVS( 227 )
  W( 58 ) = W( 58 ) + a*JVS( 228 )
  W( 59 ) = W( 59 ) + a*JVS( 229 )
  a = -W( 40 ) / JVS(          231  )
  W( 40 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 52 ) = W( 52 ) + a*JVS( 234 )
  W( 53 ) = W( 53 ) + a*JVS( 235 )
  W( 54 ) = W( 54 ) + a*JVS( 236 )
  W( 55 ) = W( 55 ) + a*JVS( 237 )
  W( 58 ) = W( 58 ) + a*JVS( 238 )
  a = -W( 41 ) / JVS(          240  )
  W( 41 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 241 )
  W( 52 ) = W( 52 ) + a*JVS( 242 )
  W( 53 ) = W( 53 ) + a*JVS( 243 )
  W( 55 ) = W( 55 ) + a*JVS( 244 )
  W( 58 ) = W( 58 ) + a*JVS( 245 )
  a = -W( 42 ) / JVS(          253  )
  W( 42 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 254 )
  W( 45 ) = W( 45 ) + a*JVS( 255 )
  W( 46 ) = W( 46 ) + a*JVS( 256 )
  W( 49 ) = W( 49 ) + a*JVS( 257 )
  W( 50 ) = W( 50 ) + a*JVS( 258 )
  W( 51 ) = W( 51 ) + a*JVS( 259 )
  W( 52 ) = W( 52 ) + a*JVS( 260 )
  W( 53 ) = W( 53 ) + a*JVS( 261 )
  W( 54 ) = W( 54 ) + a*JVS( 262 )
  W( 55 ) = W( 55 ) + a*JVS( 263 )
  W( 56 ) = W( 56 ) + a*JVS( 264 )
  W( 57 ) = W( 57 ) + a*JVS( 265 )
  W( 58 ) = W( 58 ) + a*JVS( 266 )
  W( 59 ) = W( 59 ) + a*JVS( 267 )
  a = -W( 43 ) / JVS(          284  )
  W( 43 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 285 )
  W( 45 ) = W( 45 ) + a*JVS( 286 )
  W( 46 ) = W( 46 ) + a*JVS( 287 )
  W( 49 ) = W( 49 ) + a*JVS( 288 )
  W( 50 ) = W( 50 ) + a*JVS( 289 )
  W( 51 ) = W( 51 ) + a*JVS( 290 )
  W( 52 ) = W( 52 ) + a*JVS( 291 )
  W( 53 ) = W( 53 ) + a*JVS( 292 )
  W( 54 ) = W( 54 ) + a*JVS( 293 )
  W( 55 ) = W( 55 ) + a*JVS( 294 )
  W( 56 ) = W( 56 ) + a*JVS( 295 )
  W( 57 ) = W( 57 ) + a*JVS( 296 )
  W( 58 ) = W( 58 ) + a*JVS( 297 )
  W( 59 ) = W( 59 ) + a*JVS( 298 )
  a = -W( 44 ) / JVS(          300  )
  W( 44 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 301 )
  W( 52 ) = W( 52 ) + a*JVS( 302 )
  W( 53 ) = W( 53 ) + a*JVS( 303 )
  W( 55 ) = W( 55 ) + a*JVS( 304 )
  W( 58 ) = W( 58 ) + a*JVS( 305 )
  a = -W( 45 ) / JVS(          311  )
  W( 45 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 312 )
  W( 52 ) = W( 52 ) + a*JVS( 313 )
  W( 53 ) = W( 53 ) + a*JVS( 314 )
  W( 54 ) = W( 54 ) + a*JVS( 315 )
  W( 55 ) = W( 55 ) + a*JVS( 316 )
  W( 58 ) = W( 58 ) + a*JVS( 317 )
  W( 59 ) = W( 59 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 325 )
  W( 53 ) = W( 53 ) + a*JVS( 326 )
  W( 54 ) = W( 54 ) + a*JVS( 327 )
  W( 58 ) = W( 58 ) + a*JVS( 328 )
  W( 59 ) = W( 59 ) + a*JVS( 329 )
  a = -W( 47 ) / JVS(          335  )
  W( 47 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 336 )
  W( 52 ) = W( 52 ) + a*JVS( 337 )
  W( 53 ) = W( 53 ) + a*JVS( 338 )
  W( 54 ) = W( 54 ) + a*JVS( 339 )
  W( 55 ) = W( 55 ) + a*JVS( 340 )
  W( 56 ) = W( 56 ) + a*JVS( 341 )
  W( 58 ) = W( 58 ) + a*JVS( 342 )
  W( 59 ) = W( 59 ) + a*JVS( 343 )
  a = -W( 48 ) / JVS(          357  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 358 )
  W( 50 ) = W( 50 ) + a*JVS( 359 )
  W( 51 ) = W( 51 ) + a*JVS( 360 )
  W( 52 ) = W( 52 ) + a*JVS( 361 )
  W( 53 ) = W( 53 ) + a*JVS( 362 )
  W( 54 ) = W( 54 ) + a*JVS( 363 )
  W( 55 ) = W( 55 ) + a*JVS( 364 )
  W( 56 ) = W( 56 ) + a*JVS( 365 )
  W( 57 ) = W( 57 ) + a*JVS( 366 )
  W( 58 ) = W( 58 ) + a*JVS( 367 )
  W( 59 ) = W( 59 ) + a*JVS( 368 )
  a = -W( 49 ) / JVS(          371  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 372 )
  W( 52 ) = W( 52 ) + a*JVS( 373 )
  W( 53 ) = W( 53 ) + a*JVS( 374 )
  W( 54 ) = W( 54 ) + a*JVS( 375 )
  W( 55 ) = W( 55 ) + a*JVS( 376 )
  W( 58 ) = W( 58 ) + a*JVS( 377 )
  W( 59 ) = W( 59 ) + a*JVS( 378 )
  a = -W( 50 ) / JVS(          385  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 386 )
  W( 52 ) = W( 52 ) + a*JVS( 387 )
  W( 53 ) = W( 53 ) + a*JVS( 388 )
  W( 54 ) = W( 54 ) + a*JVS( 389 )
  W( 55 ) = W( 55 ) + a*JVS( 390 )
  W( 56 ) = W( 56 ) + a*JVS( 391 )
  W( 58 ) = W( 58 ) + a*JVS( 392 )
  W( 59 ) = W( 59 ) + a*JVS( 393 )
  JVS( 394) = W( 5 )
  JVS( 395) = W( 10 )
  JVS( 396) = W( 11 )
  JVS( 397) = W( 14 )
  JVS( 398) = W( 16 )
  JVS( 399) = W( 20 )
  JVS( 400) = W( 21 )
  JVS( 401) = W( 22 )
  JVS( 402) = W( 24 )
  JVS( 403) = W( 26 )
  JVS( 404) = W( 27 )
  JVS( 405) = W( 28 )
  JVS( 406) = W( 29 )
  JVS( 407) = W( 30 )
  JVS( 408) = W( 31 )
  JVS( 409) = W( 32 )
  JVS( 410) = W( 34 )
  JVS( 411) = W( 35 )
  JVS( 412) = W( 36 )
  JVS( 413) = W( 37 )
  JVS( 414) = W( 38 )
  JVS( 415) = W( 39 )
  JVS( 416) = W( 40 )
  JVS( 417) = W( 41 )
  JVS( 418) = W( 42 )
  JVS( 419) = W( 43 )
  JVS( 420) = W( 44 )
  JVS( 421) = W( 45 )
  JVS( 422) = W( 46 )
  JVS( 423) = W( 47 )
  JVS( 424) = W( 48 )
  JVS( 425) = W( 49 )
  JVS( 426) = W( 50 )
  JVS( 427) = W( 51 )
  JVS( 428) = W( 52 )
  JVS( 429) = W( 53 )
  JVS( 430) = W( 54 )
  JVS( 431) = W( 55 )
  JVS( 432) = W( 56 )
  JVS( 433) = W( 57 )
  JVS( 434) = W( 58 )
  JVS( 435) = W( 59 )
  IF ( ABS(  JVS( 457 )) < TINY(a) ) THEN
         IER = 52                                      
         RETURN
  END IF
   W( 17 ) = JVS( 436 )
   W( 18 ) = JVS( 437 )
   W( 21 ) = JVS( 438 )
   W( 27 ) = JVS( 439 )
   W( 28 ) = JVS( 440 )
   W( 29 ) = JVS( 441 )
   W( 32 ) = JVS( 442 )
   W( 35 ) = JVS( 443 )
   W( 36 ) = JVS( 444 )
   W( 37 ) = JVS( 445 )
   W( 38 ) = JVS( 446 )
   W( 39 ) = JVS( 447 )
   W( 40 ) = JVS( 448 )
   W( 41 ) = JVS( 449 )
   W( 44 ) = JVS( 450 )
   W( 45 ) = JVS( 451 )
   W( 46 ) = JVS( 452 )
   W( 48 ) = JVS( 453 )
   W( 49 ) = JVS( 454 )
   W( 50 ) = JVS( 455 )
   W( 51 ) = JVS( 456 )
   W( 52 ) = JVS( 457 )
   W( 53 ) = JVS( 458 )
   W( 54 ) = JVS( 459 )
   W( 55 ) = JVS( 460 )
   W( 56 ) = JVS( 461 )
   W( 57 ) = JVS( 462 )
   W( 58 ) = JVS( 463 )
   W( 59 ) = JVS( 464 )
  a = -W( 17 ) / JVS(           58  )
  W( 17 ) = -a
  W( 28 ) = W( 28 ) + a*JVS( 59 )
  W( 29 ) = W( 29 ) + a*JVS( 60 )
  W( 46 ) = W( 46 ) + a*JVS( 61 )
  W( 53 ) = W( 53 ) + a*JVS( 62 )
  a = -W( 18 ) / JVS(           63  )
  W( 18 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 64 )
  W( 53 ) = W( 53 ) + a*JVS( 65 )
  W( 55 ) = W( 55 ) + a*JVS( 66 )
  a = -W( 21 ) / JVS(           76  )
  W( 21 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 77 )
  W( 52 ) = W( 52 ) + a*JVS( 78 )
  W( 53 ) = W( 53 ) + a*JVS( 79 )
  a = -W( 27 ) / JVS(          121  )
  W( 27 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 122 )
  W( 53 ) = W( 53 ) + a*JVS( 123 )
  W( 54 ) = W( 54 ) + a*JVS( 124 )
  a = -W( 28 ) / JVS(          125  )
  W( 28 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 126 )
  W( 53 ) = W( 53 ) + a*JVS( 127 )
  W( 54 ) = W( 54 ) + a*JVS( 128 )
  a = -W( 29 ) / JVS(          129  )
  W( 29 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 130 )
  W( 53 ) = W( 53 ) + a*JVS( 131 )
  W( 54 ) = W( 54 ) + a*JVS( 132 )
  a = -W( 32 ) / JVS(          150  )
  W( 32 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 151 )
  W( 52 ) = W( 52 ) + a*JVS( 152 )
  W( 53 ) = W( 53 ) + a*JVS( 153 )
  W( 54 ) = W( 54 ) + a*JVS( 154 )
  W( 55 ) = W( 55 ) + a*JVS( 155 )
  W( 59 ) = W( 59 ) + a*JVS( 156 )
  a = -W( 35 ) / JVS(          182  )
  W( 35 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 183 )
  W( 52 ) = W( 52 ) + a*JVS( 184 )
  W( 53 ) = W( 53 ) + a*JVS( 185 )
  W( 55 ) = W( 55 ) + a*JVS( 186 )
  W( 58 ) = W( 58 ) + a*JVS( 187 )
  a = -W( 36 ) / JVS(          189  )
  W( 36 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 190 )
  W( 52 ) = W( 52 ) + a*JVS( 191 )
  W( 53 ) = W( 53 ) + a*JVS( 192 )
  W( 55 ) = W( 55 ) + a*JVS( 193 )
  W( 58 ) = W( 58 ) + a*JVS( 194 )
  a = -W( 37 ) / JVS(          197  )
  W( 37 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 198 )
  W( 51 ) = W( 51 ) + a*JVS( 199 )
  W( 52 ) = W( 52 ) + a*JVS( 200 )
  W( 53 ) = W( 53 ) + a*JVS( 201 )
  W( 54 ) = W( 54 ) + a*JVS( 202 )
  W( 55 ) = W( 55 ) + a*JVS( 203 )
  W( 58 ) = W( 58 ) + a*JVS( 204 )
  a = -W( 38 ) / JVS(          209  )
  W( 38 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 210 )
  W( 51 ) = W( 51 ) + a*JVS( 211 )
  W( 52 ) = W( 52 ) + a*JVS( 212 )
  W( 53 ) = W( 53 ) + a*JVS( 213 )
  W( 54 ) = W( 54 ) + a*JVS( 214 )
  W( 55 ) = W( 55 ) + a*JVS( 215 )
  W( 58 ) = W( 58 ) + a*JVS( 216 )
  a = -W( 39 ) / JVS(          221  )
  W( 39 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 222 )
  W( 51 ) = W( 51 ) + a*JVS( 223 )
  W( 52 ) = W( 52 ) + a*JVS( 224 )
  W( 53 ) = W( 53 ) + a*JVS( 225 )
  W( 54 ) = W( 54 ) + a*JVS( 226 )
  W( 55 ) = W( 55 ) + a*JVS( 227 )
  W( 58 ) = W( 58 ) + a*JVS( 228 )
  W( 59 ) = W( 59 ) + a*JVS( 229 )
  a = -W( 40 ) / JVS(          231  )
  W( 40 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 52 ) = W( 52 ) + a*JVS( 234 )
  W( 53 ) = W( 53 ) + a*JVS( 235 )
  W( 54 ) = W( 54 ) + a*JVS( 236 )
  W( 55 ) = W( 55 ) + a*JVS( 237 )
  W( 58 ) = W( 58 ) + a*JVS( 238 )
  a = -W( 41 ) / JVS(          240  )
  W( 41 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 241 )
  W( 52 ) = W( 52 ) + a*JVS( 242 )
  W( 53 ) = W( 53 ) + a*JVS( 243 )
  W( 55 ) = W( 55 ) + a*JVS( 244 )
  W( 58 ) = W( 58 ) + a*JVS( 245 )
  a = -W( 44 ) / JVS(          300  )
  W( 44 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 301 )
  W( 52 ) = W( 52 ) + a*JVS( 302 )
  W( 53 ) = W( 53 ) + a*JVS( 303 )
  W( 55 ) = W( 55 ) + a*JVS( 304 )
  W( 58 ) = W( 58 ) + a*JVS( 305 )
  a = -W( 45 ) / JVS(          311  )
  W( 45 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 312 )
  W( 52 ) = W( 52 ) + a*JVS( 313 )
  W( 53 ) = W( 53 ) + a*JVS( 314 )
  W( 54 ) = W( 54 ) + a*JVS( 315 )
  W( 55 ) = W( 55 ) + a*JVS( 316 )
  W( 58 ) = W( 58 ) + a*JVS( 317 )
  W( 59 ) = W( 59 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 325 )
  W( 53 ) = W( 53 ) + a*JVS( 326 )
  W( 54 ) = W( 54 ) + a*JVS( 327 )
  W( 58 ) = W( 58 ) + a*JVS( 328 )
  W( 59 ) = W( 59 ) + a*JVS( 329 )
  a = -W( 48 ) / JVS(          357  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 358 )
  W( 50 ) = W( 50 ) + a*JVS( 359 )
  W( 51 ) = W( 51 ) + a*JVS( 360 )
  W( 52 ) = W( 52 ) + a*JVS( 361 )
  W( 53 ) = W( 53 ) + a*JVS( 362 )
  W( 54 ) = W( 54 ) + a*JVS( 363 )
  W( 55 ) = W( 55 ) + a*JVS( 364 )
  W( 56 ) = W( 56 ) + a*JVS( 365 )
  W( 57 ) = W( 57 ) + a*JVS( 366 )
  W( 58 ) = W( 58 ) + a*JVS( 367 )
  W( 59 ) = W( 59 ) + a*JVS( 368 )
  a = -W( 49 ) / JVS(          371  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 372 )
  W( 52 ) = W( 52 ) + a*JVS( 373 )
  W( 53 ) = W( 53 ) + a*JVS( 374 )
  W( 54 ) = W( 54 ) + a*JVS( 375 )
  W( 55 ) = W( 55 ) + a*JVS( 376 )
  W( 58 ) = W( 58 ) + a*JVS( 377 )
  W( 59 ) = W( 59 ) + a*JVS( 378 )
  a = -W( 50 ) / JVS(          385  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 386 )
  W( 52 ) = W( 52 ) + a*JVS( 387 )
  W( 53 ) = W( 53 ) + a*JVS( 388 )
  W( 54 ) = W( 54 ) + a*JVS( 389 )
  W( 55 ) = W( 55 ) + a*JVS( 390 )
  W( 56 ) = W( 56 ) + a*JVS( 391 )
  W( 58 ) = W( 58 ) + a*JVS( 392 )
  W( 59 ) = W( 59 ) + a*JVS( 393 )
  a = -W( 51 ) / JVS(          427  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 428 )
  W( 53 ) = W( 53 ) + a*JVS( 429 )
  W( 54 ) = W( 54 ) + a*JVS( 430 )
  W( 55 ) = W( 55 ) + a*JVS( 431 )
  W( 56 ) = W( 56 ) + a*JVS( 432 )
  W( 57 ) = W( 57 ) + a*JVS( 433 )
  W( 58 ) = W( 58 ) + a*JVS( 434 )
  W( 59 ) = W( 59 ) + a*JVS( 435 )
  JVS( 436) = W( 17 )
  JVS( 437) = W( 18 )
  JVS( 438) = W( 21 )
  JVS( 439) = W( 27 )
  JVS( 440) = W( 28 )
  JVS( 441) = W( 29 )
  JVS( 442) = W( 32 )
  JVS( 443) = W( 35 )
  JVS( 444) = W( 36 )
  JVS( 445) = W( 37 )
  JVS( 446) = W( 38 )
  JVS( 447) = W( 39 )
  JVS( 448) = W( 40 )
  JVS( 449) = W( 41 )
  JVS( 450) = W( 44 )
  JVS( 451) = W( 45 )
  JVS( 452) = W( 46 )
  JVS( 453) = W( 48 )
  JVS( 454) = W( 49 )
  JVS( 455) = W( 50 )
  JVS( 456) = W( 51 )
  JVS( 457) = W( 52 )
  JVS( 458) = W( 53 )
  JVS( 459) = W( 54 )
  JVS( 460) = W( 55 )
  JVS( 461) = W( 56 )
  JVS( 462) = W( 57 )
  JVS( 463) = W( 58 )
  JVS( 464) = W( 59 )
  IF ( ABS(  JVS( 505 )) < TINY(a) ) THEN
         IER = 53                                      
         RETURN
  END IF
   W( 5 ) = JVS( 465 )
   W( 6 ) = JVS( 466 )
   W( 7 ) = JVS( 467 )
   W( 8 ) = JVS( 468 )
   W( 9 ) = JVS( 469 )
   W( 10 ) = JVS( 470 )
   W( 11 ) = JVS( 471 )
   W( 13 ) = JVS( 472 )
   W( 14 ) = JVS( 473 )
   W( 16 ) = JVS( 474 )
   W( 17 ) = JVS( 475 )
   W( 18 ) = JVS( 476 )
   W( 20 ) = JVS( 477 )
   W( 21 ) = JVS( 478 )
   W( 22 ) = JVS( 479 )
   W( 23 ) = JVS( 480 )
   W( 24 ) = JVS( 481 )
   W( 25 ) = JVS( 482 )
   W( 26 ) = JVS( 483 )
   W( 27 ) = JVS( 484 )
   W( 28 ) = JVS( 485 )
   W( 29 ) = JVS( 486 )
   W( 30 ) = JVS( 487 )
   W( 31 ) = JVS( 488 )
   W( 33 ) = JVS( 489 )
   W( 34 ) = JVS( 490 )
   W( 35 ) = JVS( 491 )
   W( 36 ) = JVS( 492 )
   W( 41 ) = JVS( 493 )
   W( 42 ) = JVS( 494 )
   W( 43 ) = JVS( 495 )
   W( 44 ) = JVS( 496 )
   W( 45 ) = JVS( 497 )
   W( 46 ) = JVS( 498 )
   W( 47 ) = JVS( 499 )
   W( 48 ) = JVS( 500 )
   W( 49 ) = JVS( 501 )
   W( 50 ) = JVS( 502 )
   W( 51 ) = JVS( 503 )
   W( 52 ) = JVS( 504 )
   W( 53 ) = JVS( 505 )
   W( 54 ) = JVS( 506 )
   W( 55 ) = JVS( 507 )
   W( 56 ) = JVS( 508 )
   W( 57 ) = JVS( 509 )
   W( 58 ) = JVS( 510 )
   W( 59 ) = JVS( 511 )
  a = -W( 5 ) / JVS(           30  )
  W( 5 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 31 )
  a = -W( 6 ) / JVS(           32  )
  W( 6 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 33 )
  a = -W( 7 ) / JVS(           34  )
  W( 7 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 35 )
  a = -W( 8 ) / JVS(           36  )
  W( 8 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 37 )
  a = -W( 9 ) / JVS(           38  )
  W( 9 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 39 )
  a = -W( 10 ) / JVS(           40  )
  W( 10 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 41 )
  a = -W( 11 ) / JVS(           42  )
  W( 11 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 43 )
  a = -W( 13 ) / JVS(           47  )
  W( 13 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 48 )
  W( 58 ) = W( 58 ) + a*JVS( 49 )
  a = -W( 14 ) / JVS(           50  )
  W( 14 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 51 )
  W( 53 ) = W( 53 ) + a*JVS( 52 )
  a = -W( 16 ) / JVS(           56  )
  W( 16 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 57 )
  a = -W( 17 ) / JVS(           58  )
  W( 17 ) = -a
  W( 28 ) = W( 28 ) + a*JVS( 59 )
  W( 29 ) = W( 29 ) + a*JVS( 60 )
  W( 46 ) = W( 46 ) + a*JVS( 61 )
  W( 53 ) = W( 53 ) + a*JVS( 62 )
  a = -W( 18 ) / JVS(           63  )
  W( 18 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 64 )
  W( 53 ) = W( 53 ) + a*JVS( 65 )
  W( 55 ) = W( 55 ) + a*JVS( 66 )
  a = -W( 20 ) / JVS(           72  )
  W( 20 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 73 )
  W( 53 ) = W( 53 ) + a*JVS( 74 )
  W( 59 ) = W( 59 ) + a*JVS( 75 )
  a = -W( 21 ) / JVS(           76  )
  W( 21 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 77 )
  W( 52 ) = W( 52 ) + a*JVS( 78 )
  W( 53 ) = W( 53 ) + a*JVS( 79 )
  a = -W( 22 ) / JVS(           82  )
  W( 22 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 83 )
  W( 54 ) = W( 54 ) + a*JVS( 84 )
  a = -W( 23 ) / JVS(           85  )
  W( 23 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 86 )
  W( 55 ) = W( 55 ) + a*JVS( 87 )
  W( 59 ) = W( 59 ) + a*JVS( 88 )
  a = -W( 24 ) / JVS(           89  )
  W( 24 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 90 )
  W( 53 ) = W( 53 ) + a*JVS( 91 )
  W( 54 ) = W( 54 ) + a*JVS( 92 )
  a = -W( 25 ) / JVS(           95  )
  W( 25 ) = -a
  W( 30 ) = W( 30 ) + a*JVS( 96 )
  W( 31 ) = W( 31 ) + a*JVS( 97 )
  W( 34 ) = W( 34 ) + a*JVS( 98 )
  W( 43 ) = W( 43 ) + a*JVS( 99 )
  W( 48 ) = W( 48 ) + a*JVS( 100 )
  W( 51 ) = W( 51 ) + a*JVS( 101 )
  W( 53 ) = W( 53 ) + a*JVS( 102 )
  W( 54 ) = W( 54 ) + a*JVS( 103 )
  W( 59 ) = W( 59 ) + a*JVS( 104 )
  a = -W( 26 ) / JVS(          106  )
  W( 26 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 107 )
  W( 28 ) = W( 28 ) + a*JVS( 108 )
  W( 29 ) = W( 29 ) + a*JVS( 109 )
  W( 31 ) = W( 31 ) + a*JVS( 110 )
  W( 34 ) = W( 34 ) + a*JVS( 111 )
  W( 43 ) = W( 43 ) + a*JVS( 112 )
  W( 45 ) = W( 45 ) + a*JVS( 113 )
  W( 46 ) = W( 46 ) + a*JVS( 114 )
  W( 48 ) = W( 48 ) + a*JVS( 115 )
  W( 52 ) = W( 52 ) + a*JVS( 116 )
  W( 53 ) = W( 53 ) + a*JVS( 117 )
  W( 54 ) = W( 54 ) + a*JVS( 118 )
  W( 55 ) = W( 55 ) + a*JVS( 119 )
  W( 58 ) = W( 58 ) + a*JVS( 120 )
  a = -W( 27 ) / JVS(          121  )
  W( 27 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 122 )
  W( 53 ) = W( 53 ) + a*JVS( 123 )
  W( 54 ) = W( 54 ) + a*JVS( 124 )
  a = -W( 28 ) / JVS(          125  )
  W( 28 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 126 )
  W( 53 ) = W( 53 ) + a*JVS( 127 )
  W( 54 ) = W( 54 ) + a*JVS( 128 )
  a = -W( 29 ) / JVS(          129  )
  W( 29 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 130 )
  W( 53 ) = W( 53 ) + a*JVS( 131 )
  W( 54 ) = W( 54 ) + a*JVS( 132 )
  a = -W( 30 ) / JVS(          133  )
  W( 30 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 134 )
  W( 36 ) = W( 36 ) + a*JVS( 135 )
  W( 52 ) = W( 52 ) + a*JVS( 136 )
  W( 53 ) = W( 53 ) + a*JVS( 137 )
  W( 54 ) = W( 54 ) + a*JVS( 138 )
  W( 55 ) = W( 55 ) + a*JVS( 139 )
  W( 58 ) = W( 58 ) + a*JVS( 140 )
  a = -W( 31 ) / JVS(          141  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 142 )
  W( 45 ) = W( 45 ) + a*JVS( 143 )
  W( 52 ) = W( 52 ) + a*JVS( 144 )
  W( 53 ) = W( 53 ) + a*JVS( 145 )
  W( 54 ) = W( 54 ) + a*JVS( 146 )
  W( 55 ) = W( 55 ) + a*JVS( 147 )
  W( 58 ) = W( 58 ) + a*JVS( 148 )
  a = -W( 33 ) / JVS(          159  )
  W( 33 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 160 )
  W( 44 ) = W( 44 ) + a*JVS( 161 )
  W( 46 ) = W( 46 ) + a*JVS( 162 )
  W( 47 ) = W( 47 ) + a*JVS( 163 )
  W( 49 ) = W( 49 ) + a*JVS( 164 )
  W( 52 ) = W( 52 ) + a*JVS( 165 )
  W( 53 ) = W( 53 ) + a*JVS( 166 )
  W( 54 ) = W( 54 ) + a*JVS( 167 )
  W( 55 ) = W( 55 ) + a*JVS( 168 )
  W( 56 ) = W( 56 ) + a*JVS( 169 )
  W( 58 ) = W( 58 ) + a*JVS( 170 )
  a = -W( 34 ) / JVS(          171  )
  W( 34 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 172 )
  W( 36 ) = W( 36 ) + a*JVS( 173 )
  W( 45 ) = W( 45 ) + a*JVS( 174 )
  W( 50 ) = W( 50 ) + a*JVS( 175 )
  W( 52 ) = W( 52 ) + a*JVS( 176 )
  W( 53 ) = W( 53 ) + a*JVS( 177 )
  W( 54 ) = W( 54 ) + a*JVS( 178 )
  W( 55 ) = W( 55 ) + a*JVS( 179 )
  W( 58 ) = W( 58 ) + a*JVS( 180 )
  a = -W( 35 ) / JVS(          182  )
  W( 35 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 183 )
  W( 52 ) = W( 52 ) + a*JVS( 184 )
  W( 53 ) = W( 53 ) + a*JVS( 185 )
  W( 55 ) = W( 55 ) + a*JVS( 186 )
  W( 58 ) = W( 58 ) + a*JVS( 187 )
  a = -W( 36 ) / JVS(          189  )
  W( 36 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 190 )
  W( 52 ) = W( 52 ) + a*JVS( 191 )
  W( 53 ) = W( 53 ) + a*JVS( 192 )
  W( 55 ) = W( 55 ) + a*JVS( 193 )
  W( 58 ) = W( 58 ) + a*JVS( 194 )
  a = -W( 41 ) / JVS(          240  )
  W( 41 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 241 )
  W( 52 ) = W( 52 ) + a*JVS( 242 )
  W( 53 ) = W( 53 ) + a*JVS( 243 )
  W( 55 ) = W( 55 ) + a*JVS( 244 )
  W( 58 ) = W( 58 ) + a*JVS( 245 )
  a = -W( 42 ) / JVS(          253  )
  W( 42 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 254 )
  W( 45 ) = W( 45 ) + a*JVS( 255 )
  W( 46 ) = W( 46 ) + a*JVS( 256 )
  W( 49 ) = W( 49 ) + a*JVS( 257 )
  W( 50 ) = W( 50 ) + a*JVS( 258 )
  W( 51 ) = W( 51 ) + a*JVS( 259 )
  W( 52 ) = W( 52 ) + a*JVS( 260 )
  W( 53 ) = W( 53 ) + a*JVS( 261 )
  W( 54 ) = W( 54 ) + a*JVS( 262 )
  W( 55 ) = W( 55 ) + a*JVS( 263 )
  W( 56 ) = W( 56 ) + a*JVS( 264 )
  W( 57 ) = W( 57 ) + a*JVS( 265 )
  W( 58 ) = W( 58 ) + a*JVS( 266 )
  W( 59 ) = W( 59 ) + a*JVS( 267 )
  a = -W( 43 ) / JVS(          284  )
  W( 43 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 285 )
  W( 45 ) = W( 45 ) + a*JVS( 286 )
  W( 46 ) = W( 46 ) + a*JVS( 287 )
  W( 49 ) = W( 49 ) + a*JVS( 288 )
  W( 50 ) = W( 50 ) + a*JVS( 289 )
  W( 51 ) = W( 51 ) + a*JVS( 290 )
  W( 52 ) = W( 52 ) + a*JVS( 291 )
  W( 53 ) = W( 53 ) + a*JVS( 292 )
  W( 54 ) = W( 54 ) + a*JVS( 293 )
  W( 55 ) = W( 55 ) + a*JVS( 294 )
  W( 56 ) = W( 56 ) + a*JVS( 295 )
  W( 57 ) = W( 57 ) + a*JVS( 296 )
  W( 58 ) = W( 58 ) + a*JVS( 297 )
  W( 59 ) = W( 59 ) + a*JVS( 298 )
  a = -W( 44 ) / JVS(          300  )
  W( 44 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 301 )
  W( 52 ) = W( 52 ) + a*JVS( 302 )
  W( 53 ) = W( 53 ) + a*JVS( 303 )
  W( 55 ) = W( 55 ) + a*JVS( 304 )
  W( 58 ) = W( 58 ) + a*JVS( 305 )
  a = -W( 45 ) / JVS(          311  )
  W( 45 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 312 )
  W( 52 ) = W( 52 ) + a*JVS( 313 )
  W( 53 ) = W( 53 ) + a*JVS( 314 )
  W( 54 ) = W( 54 ) + a*JVS( 315 )
  W( 55 ) = W( 55 ) + a*JVS( 316 )
  W( 58 ) = W( 58 ) + a*JVS( 317 )
  W( 59 ) = W( 59 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 325 )
  W( 53 ) = W( 53 ) + a*JVS( 326 )
  W( 54 ) = W( 54 ) + a*JVS( 327 )
  W( 58 ) = W( 58 ) + a*JVS( 328 )
  W( 59 ) = W( 59 ) + a*JVS( 329 )
  a = -W( 47 ) / JVS(          335  )
  W( 47 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 336 )
  W( 52 ) = W( 52 ) + a*JVS( 337 )
  W( 53 ) = W( 53 ) + a*JVS( 338 )
  W( 54 ) = W( 54 ) + a*JVS( 339 )
  W( 55 ) = W( 55 ) + a*JVS( 340 )
  W( 56 ) = W( 56 ) + a*JVS( 341 )
  W( 58 ) = W( 58 ) + a*JVS( 342 )
  W( 59 ) = W( 59 ) + a*JVS( 343 )
  a = -W( 48 ) / JVS(          357  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 358 )
  W( 50 ) = W( 50 ) + a*JVS( 359 )
  W( 51 ) = W( 51 ) + a*JVS( 360 )
  W( 52 ) = W( 52 ) + a*JVS( 361 )
  W( 53 ) = W( 53 ) + a*JVS( 362 )
  W( 54 ) = W( 54 ) + a*JVS( 363 )
  W( 55 ) = W( 55 ) + a*JVS( 364 )
  W( 56 ) = W( 56 ) + a*JVS( 365 )
  W( 57 ) = W( 57 ) + a*JVS( 366 )
  W( 58 ) = W( 58 ) + a*JVS( 367 )
  W( 59 ) = W( 59 ) + a*JVS( 368 )
  a = -W( 49 ) / JVS(          371  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 372 )
  W( 52 ) = W( 52 ) + a*JVS( 373 )
  W( 53 ) = W( 53 ) + a*JVS( 374 )
  W( 54 ) = W( 54 ) + a*JVS( 375 )
  W( 55 ) = W( 55 ) + a*JVS( 376 )
  W( 58 ) = W( 58 ) + a*JVS( 377 )
  W( 59 ) = W( 59 ) + a*JVS( 378 )
  a = -W( 50 ) / JVS(          385  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 386 )
  W( 52 ) = W( 52 ) + a*JVS( 387 )
  W( 53 ) = W( 53 ) + a*JVS( 388 )
  W( 54 ) = W( 54 ) + a*JVS( 389 )
  W( 55 ) = W( 55 ) + a*JVS( 390 )
  W( 56 ) = W( 56 ) + a*JVS( 391 )
  W( 58 ) = W( 58 ) + a*JVS( 392 )
  W( 59 ) = W( 59 ) + a*JVS( 393 )
  a = -W( 51 ) / JVS(          427  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 428 )
  W( 53 ) = W( 53 ) + a*JVS( 429 )
  W( 54 ) = W( 54 ) + a*JVS( 430 )
  W( 55 ) = W( 55 ) + a*JVS( 431 )
  W( 56 ) = W( 56 ) + a*JVS( 432 )
  W( 57 ) = W( 57 ) + a*JVS( 433 )
  W( 58 ) = W( 58 ) + a*JVS( 434 )
  W( 59 ) = W( 59 ) + a*JVS( 435 )
  a = -W( 52 ) / JVS(          457  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 458 )
  W( 54 ) = W( 54 ) + a*JVS( 459 )
  W( 55 ) = W( 55 ) + a*JVS( 460 )
  W( 56 ) = W( 56 ) + a*JVS( 461 )
  W( 57 ) = W( 57 ) + a*JVS( 462 )
  W( 58 ) = W( 58 ) + a*JVS( 463 )
  W( 59 ) = W( 59 ) + a*JVS( 464 )
  JVS( 465) = W( 5 )
  JVS( 466) = W( 6 )
  JVS( 467) = W( 7 )
  JVS( 468) = W( 8 )
  JVS( 469) = W( 9 )
  JVS( 470) = W( 10 )
  JVS( 471) = W( 11 )
  JVS( 472) = W( 13 )
  JVS( 473) = W( 14 )
  JVS( 474) = W( 16 )
  JVS( 475) = W( 17 )
  JVS( 476) = W( 18 )
  JVS( 477) = W( 20 )
  JVS( 478) = W( 21 )
  JVS( 479) = W( 22 )
  JVS( 480) = W( 23 )
  JVS( 481) = W( 24 )
  JVS( 482) = W( 25 )
  JVS( 483) = W( 26 )
  JVS( 484) = W( 27 )
  JVS( 485) = W( 28 )
  JVS( 486) = W( 29 )
  JVS( 487) = W( 30 )
  JVS( 488) = W( 31 )
  JVS( 489) = W( 33 )
  JVS( 490) = W( 34 )
  JVS( 491) = W( 35 )
  JVS( 492) = W( 36 )
  JVS( 493) = W( 41 )
  JVS( 494) = W( 42 )
  JVS( 495) = W( 43 )
  JVS( 496) = W( 44 )
  JVS( 497) = W( 45 )
  JVS( 498) = W( 46 )
  JVS( 499) = W( 47 )
  JVS( 500) = W( 48 )
  JVS( 501) = W( 49 )
  JVS( 502) = W( 50 )
  JVS( 503) = W( 51 )
  JVS( 504) = W( 52 )
  JVS( 505) = W( 53 )
  JVS( 506) = W( 54 )
  JVS( 507) = W( 55 )
  JVS( 508) = W( 56 )
  JVS( 509) = W( 57 )
  JVS( 510) = W( 58 )
  JVS( 511) = W( 59 )
  IF ( ABS(  JVS( 536 )) < TINY(a) ) THEN
         IER = 54                                      
         RETURN
  END IF
   W( 15 ) = JVS( 512 )
   W( 20 ) = JVS( 513 )
   W( 22 ) = JVS( 514 )
   W( 23 ) = JVS( 515 )
   W( 24 ) = JVS( 516 )
   W( 25 ) = JVS( 517 )
   W( 27 ) = JVS( 518 )
   W( 28 ) = JVS( 519 )
   W( 29 ) = JVS( 520 )
   W( 30 ) = JVS( 521 )
   W( 31 ) = JVS( 522 )
   W( 34 ) = JVS( 523 )
   W( 35 ) = JVS( 524 )
   W( 36 ) = JVS( 525 )
   W( 43 ) = JVS( 526 )
   W( 44 ) = JVS( 527 )
   W( 45 ) = JVS( 528 )
   W( 46 ) = JVS( 529 )
   W( 48 ) = JVS( 530 )
   W( 49 ) = JVS( 531 )
   W( 50 ) = JVS( 532 )
   W( 51 ) = JVS( 533 )
   W( 52 ) = JVS( 534 )
   W( 53 ) = JVS( 535 )
   W( 54 ) = JVS( 536 )
   W( 55 ) = JVS( 537 )
   W( 56 ) = JVS( 538 )
   W( 57 ) = JVS( 539 )
   W( 58 ) = JVS( 540 )
   W( 59 ) = JVS( 541 )
  a = -W( 15 ) / JVS(           53  )
  W( 15 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 54 )
  W( 59 ) = W( 59 ) + a*JVS( 55 )
  a = -W( 20 ) / JVS(           72  )
  W( 20 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 73 )
  W( 53 ) = W( 53 ) + a*JVS( 74 )
  W( 59 ) = W( 59 ) + a*JVS( 75 )
  a = -W( 22 ) / JVS(           82  )
  W( 22 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 83 )
  W( 54 ) = W( 54 ) + a*JVS( 84 )
  a = -W( 23 ) / JVS(           85  )
  W( 23 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 86 )
  W( 55 ) = W( 55 ) + a*JVS( 87 )
  W( 59 ) = W( 59 ) + a*JVS( 88 )
  a = -W( 24 ) / JVS(           89  )
  W( 24 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 90 )
  W( 53 ) = W( 53 ) + a*JVS( 91 )
  W( 54 ) = W( 54 ) + a*JVS( 92 )
  a = -W( 25 ) / JVS(           95  )
  W( 25 ) = -a
  W( 30 ) = W( 30 ) + a*JVS( 96 )
  W( 31 ) = W( 31 ) + a*JVS( 97 )
  W( 34 ) = W( 34 ) + a*JVS( 98 )
  W( 43 ) = W( 43 ) + a*JVS( 99 )
  W( 48 ) = W( 48 ) + a*JVS( 100 )
  W( 51 ) = W( 51 ) + a*JVS( 101 )
  W( 53 ) = W( 53 ) + a*JVS( 102 )
  W( 54 ) = W( 54 ) + a*JVS( 103 )
  W( 59 ) = W( 59 ) + a*JVS( 104 )
  a = -W( 27 ) / JVS(          121  )
  W( 27 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 122 )
  W( 53 ) = W( 53 ) + a*JVS( 123 )
  W( 54 ) = W( 54 ) + a*JVS( 124 )
  a = -W( 28 ) / JVS(          125  )
  W( 28 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 126 )
  W( 53 ) = W( 53 ) + a*JVS( 127 )
  W( 54 ) = W( 54 ) + a*JVS( 128 )
  a = -W( 29 ) / JVS(          129  )
  W( 29 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 130 )
  W( 53 ) = W( 53 ) + a*JVS( 131 )
  W( 54 ) = W( 54 ) + a*JVS( 132 )
  a = -W( 30 ) / JVS(          133  )
  W( 30 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 134 )
  W( 36 ) = W( 36 ) + a*JVS( 135 )
  W( 52 ) = W( 52 ) + a*JVS( 136 )
  W( 53 ) = W( 53 ) + a*JVS( 137 )
  W( 54 ) = W( 54 ) + a*JVS( 138 )
  W( 55 ) = W( 55 ) + a*JVS( 139 )
  W( 58 ) = W( 58 ) + a*JVS( 140 )
  a = -W( 31 ) / JVS(          141  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 142 )
  W( 45 ) = W( 45 ) + a*JVS( 143 )
  W( 52 ) = W( 52 ) + a*JVS( 144 )
  W( 53 ) = W( 53 ) + a*JVS( 145 )
  W( 54 ) = W( 54 ) + a*JVS( 146 )
  W( 55 ) = W( 55 ) + a*JVS( 147 )
  W( 58 ) = W( 58 ) + a*JVS( 148 )
  a = -W( 34 ) / JVS(          171  )
  W( 34 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 172 )
  W( 36 ) = W( 36 ) + a*JVS( 173 )
  W( 45 ) = W( 45 ) + a*JVS( 174 )
  W( 50 ) = W( 50 ) + a*JVS( 175 )
  W( 52 ) = W( 52 ) + a*JVS( 176 )
  W( 53 ) = W( 53 ) + a*JVS( 177 )
  W( 54 ) = W( 54 ) + a*JVS( 178 )
  W( 55 ) = W( 55 ) + a*JVS( 179 )
  W( 58 ) = W( 58 ) + a*JVS( 180 )
  a = -W( 35 ) / JVS(          182  )
  W( 35 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 183 )
  W( 52 ) = W( 52 ) + a*JVS( 184 )
  W( 53 ) = W( 53 ) + a*JVS( 185 )
  W( 55 ) = W( 55 ) + a*JVS( 186 )
  W( 58 ) = W( 58 ) + a*JVS( 187 )
  a = -W( 36 ) / JVS(          189  )
  W( 36 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 190 )
  W( 52 ) = W( 52 ) + a*JVS( 191 )
  W( 53 ) = W( 53 ) + a*JVS( 192 )
  W( 55 ) = W( 55 ) + a*JVS( 193 )
  W( 58 ) = W( 58 ) + a*JVS( 194 )
  a = -W( 43 ) / JVS(          284  )
  W( 43 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 285 )
  W( 45 ) = W( 45 ) + a*JVS( 286 )
  W( 46 ) = W( 46 ) + a*JVS( 287 )
  W( 49 ) = W( 49 ) + a*JVS( 288 )
  W( 50 ) = W( 50 ) + a*JVS( 289 )
  W( 51 ) = W( 51 ) + a*JVS( 290 )
  W( 52 ) = W( 52 ) + a*JVS( 291 )
  W( 53 ) = W( 53 ) + a*JVS( 292 )
  W( 54 ) = W( 54 ) + a*JVS( 293 )
  W( 55 ) = W( 55 ) + a*JVS( 294 )
  W( 56 ) = W( 56 ) + a*JVS( 295 )
  W( 57 ) = W( 57 ) + a*JVS( 296 )
  W( 58 ) = W( 58 ) + a*JVS( 297 )
  W( 59 ) = W( 59 ) + a*JVS( 298 )
  a = -W( 44 ) / JVS(          300  )
  W( 44 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 301 )
  W( 52 ) = W( 52 ) + a*JVS( 302 )
  W( 53 ) = W( 53 ) + a*JVS( 303 )
  W( 55 ) = W( 55 ) + a*JVS( 304 )
  W( 58 ) = W( 58 ) + a*JVS( 305 )
  a = -W( 45 ) / JVS(          311  )
  W( 45 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 312 )
  W( 52 ) = W( 52 ) + a*JVS( 313 )
  W( 53 ) = W( 53 ) + a*JVS( 314 )
  W( 54 ) = W( 54 ) + a*JVS( 315 )
  W( 55 ) = W( 55 ) + a*JVS( 316 )
  W( 58 ) = W( 58 ) + a*JVS( 317 )
  W( 59 ) = W( 59 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 325 )
  W( 53 ) = W( 53 ) + a*JVS( 326 )
  W( 54 ) = W( 54 ) + a*JVS( 327 )
  W( 58 ) = W( 58 ) + a*JVS( 328 )
  W( 59 ) = W( 59 ) + a*JVS( 329 )
  a = -W( 48 ) / JVS(          357  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 358 )
  W( 50 ) = W( 50 ) + a*JVS( 359 )
  W( 51 ) = W( 51 ) + a*JVS( 360 )
  W( 52 ) = W( 52 ) + a*JVS( 361 )
  W( 53 ) = W( 53 ) + a*JVS( 362 )
  W( 54 ) = W( 54 ) + a*JVS( 363 )
  W( 55 ) = W( 55 ) + a*JVS( 364 )
  W( 56 ) = W( 56 ) + a*JVS( 365 )
  W( 57 ) = W( 57 ) + a*JVS( 366 )
  W( 58 ) = W( 58 ) + a*JVS( 367 )
  W( 59 ) = W( 59 ) + a*JVS( 368 )
  a = -W( 49 ) / JVS(          371  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 372 )
  W( 52 ) = W( 52 ) + a*JVS( 373 )
  W( 53 ) = W( 53 ) + a*JVS( 374 )
  W( 54 ) = W( 54 ) + a*JVS( 375 )
  W( 55 ) = W( 55 ) + a*JVS( 376 )
  W( 58 ) = W( 58 ) + a*JVS( 377 )
  W( 59 ) = W( 59 ) + a*JVS( 378 )
  a = -W( 50 ) / JVS(          385  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 386 )
  W( 52 ) = W( 52 ) + a*JVS( 387 )
  W( 53 ) = W( 53 ) + a*JVS( 388 )
  W( 54 ) = W( 54 ) + a*JVS( 389 )
  W( 55 ) = W( 55 ) + a*JVS( 390 )
  W( 56 ) = W( 56 ) + a*JVS( 391 )
  W( 58 ) = W( 58 ) + a*JVS( 392 )
  W( 59 ) = W( 59 ) + a*JVS( 393 )
  a = -W( 51 ) / JVS(          427  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 428 )
  W( 53 ) = W( 53 ) + a*JVS( 429 )
  W( 54 ) = W( 54 ) + a*JVS( 430 )
  W( 55 ) = W( 55 ) + a*JVS( 431 )
  W( 56 ) = W( 56 ) + a*JVS( 432 )
  W( 57 ) = W( 57 ) + a*JVS( 433 )
  W( 58 ) = W( 58 ) + a*JVS( 434 )
  W( 59 ) = W( 59 ) + a*JVS( 435 )
  a = -W( 52 ) / JVS(          457  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 458 )
  W( 54 ) = W( 54 ) + a*JVS( 459 )
  W( 55 ) = W( 55 ) + a*JVS( 460 )
  W( 56 ) = W( 56 ) + a*JVS( 461 )
  W( 57 ) = W( 57 ) + a*JVS( 462 )
  W( 58 ) = W( 58 ) + a*JVS( 463 )
  W( 59 ) = W( 59 ) + a*JVS( 464 )
  a = -W( 53 ) / JVS(          505  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 506 )
  W( 55 ) = W( 55 ) + a*JVS( 507 )
  W( 56 ) = W( 56 ) + a*JVS( 508 )
  W( 57 ) = W( 57 ) + a*JVS( 509 )
  W( 58 ) = W( 58 ) + a*JVS( 510 )
  W( 59 ) = W( 59 ) + a*JVS( 511 )
  JVS( 512) = W( 15 )
  JVS( 513) = W( 20 )
  JVS( 514) = W( 22 )
  JVS( 515) = W( 23 )
  JVS( 516) = W( 24 )
  JVS( 517) = W( 25 )
  JVS( 518) = W( 27 )
  JVS( 519) = W( 28 )
  JVS( 520) = W( 29 )
  JVS( 521) = W( 30 )
  JVS( 522) = W( 31 )
  JVS( 523) = W( 34 )
  JVS( 524) = W( 35 )
  JVS( 525) = W( 36 )
  JVS( 526) = W( 43 )
  JVS( 527) = W( 44 )
  JVS( 528) = W( 45 )
  JVS( 529) = W( 46 )
  JVS( 530) = W( 48 )
  JVS( 531) = W( 49 )
  JVS( 532) = W( 50 )
  JVS( 533) = W( 51 )
  JVS( 534) = W( 52 )
  JVS( 535) = W( 53 )
  JVS( 536) = W( 54 )
  JVS( 537) = W( 55 )
  JVS( 538) = W( 56 )
  JVS( 539) = W( 57 )
  JVS( 540) = W( 58 )
  JVS( 541) = W( 59 )
  IF ( ABS(  JVS( 565 )) < TINY(a) ) THEN
         IER = 55                                      
         RETURN
  END IF
   W( 18 ) = JVS( 542 )
   W( 23 ) = JVS( 543 )
   W( 32 ) = JVS( 544 )
   W( 33 ) = JVS( 545 )
   W( 34 ) = JVS( 546 )
   W( 35 ) = JVS( 547 )
   W( 36 ) = JVS( 548 )
   W( 37 ) = JVS( 549 )
   W( 38 ) = JVS( 550 )
   W( 39 ) = JVS( 551 )
   W( 40 ) = JVS( 552 )
   W( 41 ) = JVS( 553 )
   W( 44 ) = JVS( 554 )
   W( 45 ) = JVS( 555 )
   W( 46 ) = JVS( 556 )
   W( 47 ) = JVS( 557 )
   W( 48 ) = JVS( 558 )
   W( 49 ) = JVS( 559 )
   W( 50 ) = JVS( 560 )
   W( 51 ) = JVS( 561 )
   W( 52 ) = JVS( 562 )
   W( 53 ) = JVS( 563 )
   W( 54 ) = JVS( 564 )
   W( 55 ) = JVS( 565 )
   W( 56 ) = JVS( 566 )
   W( 57 ) = JVS( 567 )
   W( 58 ) = JVS( 568 )
   W( 59 ) = JVS( 569 )
  a = -W( 18 ) / JVS(           63  )
  W( 18 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 64 )
  W( 53 ) = W( 53 ) + a*JVS( 65 )
  W( 55 ) = W( 55 ) + a*JVS( 66 )
  a = -W( 23 ) / JVS(           85  )
  W( 23 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 86 )
  W( 55 ) = W( 55 ) + a*JVS( 87 )
  W( 59 ) = W( 59 ) + a*JVS( 88 )
  a = -W( 32 ) / JVS(          150  )
  W( 32 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 151 )
  W( 52 ) = W( 52 ) + a*JVS( 152 )
  W( 53 ) = W( 53 ) + a*JVS( 153 )
  W( 54 ) = W( 54 ) + a*JVS( 154 )
  W( 55 ) = W( 55 ) + a*JVS( 155 )
  W( 59 ) = W( 59 ) + a*JVS( 156 )
  a = -W( 33 ) / JVS(          159  )
  W( 33 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 160 )
  W( 44 ) = W( 44 ) + a*JVS( 161 )
  W( 46 ) = W( 46 ) + a*JVS( 162 )
  W( 47 ) = W( 47 ) + a*JVS( 163 )
  W( 49 ) = W( 49 ) + a*JVS( 164 )
  W( 52 ) = W( 52 ) + a*JVS( 165 )
  W( 53 ) = W( 53 ) + a*JVS( 166 )
  W( 54 ) = W( 54 ) + a*JVS( 167 )
  W( 55 ) = W( 55 ) + a*JVS( 168 )
  W( 56 ) = W( 56 ) + a*JVS( 169 )
  W( 58 ) = W( 58 ) + a*JVS( 170 )
  a = -W( 34 ) / JVS(          171  )
  W( 34 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 172 )
  W( 36 ) = W( 36 ) + a*JVS( 173 )
  W( 45 ) = W( 45 ) + a*JVS( 174 )
  W( 50 ) = W( 50 ) + a*JVS( 175 )
  W( 52 ) = W( 52 ) + a*JVS( 176 )
  W( 53 ) = W( 53 ) + a*JVS( 177 )
  W( 54 ) = W( 54 ) + a*JVS( 178 )
  W( 55 ) = W( 55 ) + a*JVS( 179 )
  W( 58 ) = W( 58 ) + a*JVS( 180 )
  a = -W( 35 ) / JVS(          182  )
  W( 35 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 183 )
  W( 52 ) = W( 52 ) + a*JVS( 184 )
  W( 53 ) = W( 53 ) + a*JVS( 185 )
  W( 55 ) = W( 55 ) + a*JVS( 186 )
  W( 58 ) = W( 58 ) + a*JVS( 187 )
  a = -W( 36 ) / JVS(          189  )
  W( 36 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 190 )
  W( 52 ) = W( 52 ) + a*JVS( 191 )
  W( 53 ) = W( 53 ) + a*JVS( 192 )
  W( 55 ) = W( 55 ) + a*JVS( 193 )
  W( 58 ) = W( 58 ) + a*JVS( 194 )
  a = -W( 37 ) / JVS(          197  )
  W( 37 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 198 )
  W( 51 ) = W( 51 ) + a*JVS( 199 )
  W( 52 ) = W( 52 ) + a*JVS( 200 )
  W( 53 ) = W( 53 ) + a*JVS( 201 )
  W( 54 ) = W( 54 ) + a*JVS( 202 )
  W( 55 ) = W( 55 ) + a*JVS( 203 )
  W( 58 ) = W( 58 ) + a*JVS( 204 )
  a = -W( 38 ) / JVS(          209  )
  W( 38 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 210 )
  W( 51 ) = W( 51 ) + a*JVS( 211 )
  W( 52 ) = W( 52 ) + a*JVS( 212 )
  W( 53 ) = W( 53 ) + a*JVS( 213 )
  W( 54 ) = W( 54 ) + a*JVS( 214 )
  W( 55 ) = W( 55 ) + a*JVS( 215 )
  W( 58 ) = W( 58 ) + a*JVS( 216 )
  a = -W( 39 ) / JVS(          221  )
  W( 39 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 222 )
  W( 51 ) = W( 51 ) + a*JVS( 223 )
  W( 52 ) = W( 52 ) + a*JVS( 224 )
  W( 53 ) = W( 53 ) + a*JVS( 225 )
  W( 54 ) = W( 54 ) + a*JVS( 226 )
  W( 55 ) = W( 55 ) + a*JVS( 227 )
  W( 58 ) = W( 58 ) + a*JVS( 228 )
  W( 59 ) = W( 59 ) + a*JVS( 229 )
  a = -W( 40 ) / JVS(          231  )
  W( 40 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 52 ) = W( 52 ) + a*JVS( 234 )
  W( 53 ) = W( 53 ) + a*JVS( 235 )
  W( 54 ) = W( 54 ) + a*JVS( 236 )
  W( 55 ) = W( 55 ) + a*JVS( 237 )
  W( 58 ) = W( 58 ) + a*JVS( 238 )
  a = -W( 41 ) / JVS(          240  )
  W( 41 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 241 )
  W( 52 ) = W( 52 ) + a*JVS( 242 )
  W( 53 ) = W( 53 ) + a*JVS( 243 )
  W( 55 ) = W( 55 ) + a*JVS( 244 )
  W( 58 ) = W( 58 ) + a*JVS( 245 )
  a = -W( 44 ) / JVS(          300  )
  W( 44 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 301 )
  W( 52 ) = W( 52 ) + a*JVS( 302 )
  W( 53 ) = W( 53 ) + a*JVS( 303 )
  W( 55 ) = W( 55 ) + a*JVS( 304 )
  W( 58 ) = W( 58 ) + a*JVS( 305 )
  a = -W( 45 ) / JVS(          311  )
  W( 45 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 312 )
  W( 52 ) = W( 52 ) + a*JVS( 313 )
  W( 53 ) = W( 53 ) + a*JVS( 314 )
  W( 54 ) = W( 54 ) + a*JVS( 315 )
  W( 55 ) = W( 55 ) + a*JVS( 316 )
  W( 58 ) = W( 58 ) + a*JVS( 317 )
  W( 59 ) = W( 59 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 325 )
  W( 53 ) = W( 53 ) + a*JVS( 326 )
  W( 54 ) = W( 54 ) + a*JVS( 327 )
  W( 58 ) = W( 58 ) + a*JVS( 328 )
  W( 59 ) = W( 59 ) + a*JVS( 329 )
  a = -W( 47 ) / JVS(          335  )
  W( 47 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 336 )
  W( 52 ) = W( 52 ) + a*JVS( 337 )
  W( 53 ) = W( 53 ) + a*JVS( 338 )
  W( 54 ) = W( 54 ) + a*JVS( 339 )
  W( 55 ) = W( 55 ) + a*JVS( 340 )
  W( 56 ) = W( 56 ) + a*JVS( 341 )
  W( 58 ) = W( 58 ) + a*JVS( 342 )
  W( 59 ) = W( 59 ) + a*JVS( 343 )
  a = -W( 48 ) / JVS(          357  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 358 )
  W( 50 ) = W( 50 ) + a*JVS( 359 )
  W( 51 ) = W( 51 ) + a*JVS( 360 )
  W( 52 ) = W( 52 ) + a*JVS( 361 )
  W( 53 ) = W( 53 ) + a*JVS( 362 )
  W( 54 ) = W( 54 ) + a*JVS( 363 )
  W( 55 ) = W( 55 ) + a*JVS( 364 )
  W( 56 ) = W( 56 ) + a*JVS( 365 )
  W( 57 ) = W( 57 ) + a*JVS( 366 )
  W( 58 ) = W( 58 ) + a*JVS( 367 )
  W( 59 ) = W( 59 ) + a*JVS( 368 )
  a = -W( 49 ) / JVS(          371  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 372 )
  W( 52 ) = W( 52 ) + a*JVS( 373 )
  W( 53 ) = W( 53 ) + a*JVS( 374 )
  W( 54 ) = W( 54 ) + a*JVS( 375 )
  W( 55 ) = W( 55 ) + a*JVS( 376 )
  W( 58 ) = W( 58 ) + a*JVS( 377 )
  W( 59 ) = W( 59 ) + a*JVS( 378 )
  a = -W( 50 ) / JVS(          385  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 386 )
  W( 52 ) = W( 52 ) + a*JVS( 387 )
  W( 53 ) = W( 53 ) + a*JVS( 388 )
  W( 54 ) = W( 54 ) + a*JVS( 389 )
  W( 55 ) = W( 55 ) + a*JVS( 390 )
  W( 56 ) = W( 56 ) + a*JVS( 391 )
  W( 58 ) = W( 58 ) + a*JVS( 392 )
  W( 59 ) = W( 59 ) + a*JVS( 393 )
  a = -W( 51 ) / JVS(          427  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 428 )
  W( 53 ) = W( 53 ) + a*JVS( 429 )
  W( 54 ) = W( 54 ) + a*JVS( 430 )
  W( 55 ) = W( 55 ) + a*JVS( 431 )
  W( 56 ) = W( 56 ) + a*JVS( 432 )
  W( 57 ) = W( 57 ) + a*JVS( 433 )
  W( 58 ) = W( 58 ) + a*JVS( 434 )
  W( 59 ) = W( 59 ) + a*JVS( 435 )
  a = -W( 52 ) / JVS(          457  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 458 )
  W( 54 ) = W( 54 ) + a*JVS( 459 )
  W( 55 ) = W( 55 ) + a*JVS( 460 )
  W( 56 ) = W( 56 ) + a*JVS( 461 )
  W( 57 ) = W( 57 ) + a*JVS( 462 )
  W( 58 ) = W( 58 ) + a*JVS( 463 )
  W( 59 ) = W( 59 ) + a*JVS( 464 )
  a = -W( 53 ) / JVS(          505  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 506 )
  W( 55 ) = W( 55 ) + a*JVS( 507 )
  W( 56 ) = W( 56 ) + a*JVS( 508 )
  W( 57 ) = W( 57 ) + a*JVS( 509 )
  W( 58 ) = W( 58 ) + a*JVS( 510 )
  W( 59 ) = W( 59 ) + a*JVS( 511 )
  a = -W( 54 ) / JVS(          536  )
  W( 54 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 537 )
  W( 56 ) = W( 56 ) + a*JVS( 538 )
  W( 57 ) = W( 57 ) + a*JVS( 539 )
  W( 58 ) = W( 58 ) + a*JVS( 540 )
  W( 59 ) = W( 59 ) + a*JVS( 541 )
  JVS( 542) = W( 18 )
  JVS( 543) = W( 23 )
  JVS( 544) = W( 32 )
  JVS( 545) = W( 33 )
  JVS( 546) = W( 34 )
  JVS( 547) = W( 35 )
  JVS( 548) = W( 36 )
  JVS( 549) = W( 37 )
  JVS( 550) = W( 38 )
  JVS( 551) = W( 39 )
  JVS( 552) = W( 40 )
  JVS( 553) = W( 41 )
  JVS( 554) = W( 44 )
  JVS( 555) = W( 45 )
  JVS( 556) = W( 46 )
  JVS( 557) = W( 47 )
  JVS( 558) = W( 48 )
  JVS( 559) = W( 49 )
  JVS( 560) = W( 50 )
  JVS( 561) = W( 51 )
  JVS( 562) = W( 52 )
  JVS( 563) = W( 53 )
  JVS( 564) = W( 54 )
  JVS( 565) = W( 55 )
  JVS( 566) = W( 56 )
  JVS( 567) = W( 57 )
  JVS( 568) = W( 58 )
  JVS( 569) = W( 59 )
  IF ( ABS(  JVS( 583 )) < TINY(a) ) THEN
         IER = 56                                      
         RETURN
  END IF
   W( 16 ) = JVS( 570 )
   W( 42 ) = JVS( 571 )
   W( 44 ) = JVS( 572 )
   W( 45 ) = JVS( 573 )
   W( 46 ) = JVS( 574 )
   W( 47 ) = JVS( 575 )
   W( 49 ) = JVS( 576 )
   W( 50 ) = JVS( 577 )
   W( 51 ) = JVS( 578 )
   W( 52 ) = JVS( 579 )
   W( 53 ) = JVS( 580 )
   W( 54 ) = JVS( 581 )
   W( 55 ) = JVS( 582 )
   W( 56 ) = JVS( 583 )
   W( 57 ) = JVS( 584 )
   W( 58 ) = JVS( 585 )
   W( 59 ) = JVS( 586 )
  a = -W( 16 ) / JVS(           56  )
  W( 16 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 57 )
  a = -W( 42 ) / JVS(          253  )
  W( 42 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 254 )
  W( 45 ) = W( 45 ) + a*JVS( 255 )
  W( 46 ) = W( 46 ) + a*JVS( 256 )
  W( 49 ) = W( 49 ) + a*JVS( 257 )
  W( 50 ) = W( 50 ) + a*JVS( 258 )
  W( 51 ) = W( 51 ) + a*JVS( 259 )
  W( 52 ) = W( 52 ) + a*JVS( 260 )
  W( 53 ) = W( 53 ) + a*JVS( 261 )
  W( 54 ) = W( 54 ) + a*JVS( 262 )
  W( 55 ) = W( 55 ) + a*JVS( 263 )
  W( 56 ) = W( 56 ) + a*JVS( 264 )
  W( 57 ) = W( 57 ) + a*JVS( 265 )
  W( 58 ) = W( 58 ) + a*JVS( 266 )
  W( 59 ) = W( 59 ) + a*JVS( 267 )
  a = -W( 44 ) / JVS(          300  )
  W( 44 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 301 )
  W( 52 ) = W( 52 ) + a*JVS( 302 )
  W( 53 ) = W( 53 ) + a*JVS( 303 )
  W( 55 ) = W( 55 ) + a*JVS( 304 )
  W( 58 ) = W( 58 ) + a*JVS( 305 )
  a = -W( 45 ) / JVS(          311  )
  W( 45 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 312 )
  W( 52 ) = W( 52 ) + a*JVS( 313 )
  W( 53 ) = W( 53 ) + a*JVS( 314 )
  W( 54 ) = W( 54 ) + a*JVS( 315 )
  W( 55 ) = W( 55 ) + a*JVS( 316 )
  W( 58 ) = W( 58 ) + a*JVS( 317 )
  W( 59 ) = W( 59 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 325 )
  W( 53 ) = W( 53 ) + a*JVS( 326 )
  W( 54 ) = W( 54 ) + a*JVS( 327 )
  W( 58 ) = W( 58 ) + a*JVS( 328 )
  W( 59 ) = W( 59 ) + a*JVS( 329 )
  a = -W( 47 ) / JVS(          335  )
  W( 47 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 336 )
  W( 52 ) = W( 52 ) + a*JVS( 337 )
  W( 53 ) = W( 53 ) + a*JVS( 338 )
  W( 54 ) = W( 54 ) + a*JVS( 339 )
  W( 55 ) = W( 55 ) + a*JVS( 340 )
  W( 56 ) = W( 56 ) + a*JVS( 341 )
  W( 58 ) = W( 58 ) + a*JVS( 342 )
  W( 59 ) = W( 59 ) + a*JVS( 343 )
  a = -W( 49 ) / JVS(          371  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 372 )
  W( 52 ) = W( 52 ) + a*JVS( 373 )
  W( 53 ) = W( 53 ) + a*JVS( 374 )
  W( 54 ) = W( 54 ) + a*JVS( 375 )
  W( 55 ) = W( 55 ) + a*JVS( 376 )
  W( 58 ) = W( 58 ) + a*JVS( 377 )
  W( 59 ) = W( 59 ) + a*JVS( 378 )
  a = -W( 50 ) / JVS(          385  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 386 )
  W( 52 ) = W( 52 ) + a*JVS( 387 )
  W( 53 ) = W( 53 ) + a*JVS( 388 )
  W( 54 ) = W( 54 ) + a*JVS( 389 )
  W( 55 ) = W( 55 ) + a*JVS( 390 )
  W( 56 ) = W( 56 ) + a*JVS( 391 )
  W( 58 ) = W( 58 ) + a*JVS( 392 )
  W( 59 ) = W( 59 ) + a*JVS( 393 )
  a = -W( 51 ) / JVS(          427  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 428 )
  W( 53 ) = W( 53 ) + a*JVS( 429 )
  W( 54 ) = W( 54 ) + a*JVS( 430 )
  W( 55 ) = W( 55 ) + a*JVS( 431 )
  W( 56 ) = W( 56 ) + a*JVS( 432 )
  W( 57 ) = W( 57 ) + a*JVS( 433 )
  W( 58 ) = W( 58 ) + a*JVS( 434 )
  W( 59 ) = W( 59 ) + a*JVS( 435 )
  a = -W( 52 ) / JVS(          457  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 458 )
  W( 54 ) = W( 54 ) + a*JVS( 459 )
  W( 55 ) = W( 55 ) + a*JVS( 460 )
  W( 56 ) = W( 56 ) + a*JVS( 461 )
  W( 57 ) = W( 57 ) + a*JVS( 462 )
  W( 58 ) = W( 58 ) + a*JVS( 463 )
  W( 59 ) = W( 59 ) + a*JVS( 464 )
  a = -W( 53 ) / JVS(          505  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 506 )
  W( 55 ) = W( 55 ) + a*JVS( 507 )
  W( 56 ) = W( 56 ) + a*JVS( 508 )
  W( 57 ) = W( 57 ) + a*JVS( 509 )
  W( 58 ) = W( 58 ) + a*JVS( 510 )
  W( 59 ) = W( 59 ) + a*JVS( 511 )
  a = -W( 54 ) / JVS(          536  )
  W( 54 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 537 )
  W( 56 ) = W( 56 ) + a*JVS( 538 )
  W( 57 ) = W( 57 ) + a*JVS( 539 )
  W( 58 ) = W( 58 ) + a*JVS( 540 )
  W( 59 ) = W( 59 ) + a*JVS( 541 )
  a = -W( 55 ) / JVS(          565  )
  W( 55 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 566 )
  W( 57 ) = W( 57 ) + a*JVS( 567 )
  W( 58 ) = W( 58 ) + a*JVS( 568 )
  W( 59 ) = W( 59 ) + a*JVS( 569 )
  JVS( 570) = W( 16 )
  JVS( 571) = W( 42 )
  JVS( 572) = W( 44 )
  JVS( 573) = W( 45 )
  JVS( 574) = W( 46 )
  JVS( 575) = W( 47 )
  JVS( 576) = W( 49 )
  JVS( 577) = W( 50 )
  JVS( 578) = W( 51 )
  JVS( 579) = W( 52 )
  JVS( 580) = W( 53 )
  JVS( 581) = W( 54 )
  JVS( 582) = W( 55 )
  JVS( 583) = W( 56 )
  JVS( 584) = W( 57 )
  JVS( 585) = W( 58 )
  JVS( 586) = W( 59 )
  IF ( ABS(  JVS( 600 )) < TINY(a) ) THEN
         IER = 57                                      
         RETURN
  END IF
   W( 6 ) = JVS( 587 )
   W( 33 ) = JVS( 588 )
   W( 41 ) = JVS( 589 )
   W( 44 ) = JVS( 590 )
   W( 46 ) = JVS( 591 )
   W( 47 ) = JVS( 592 )
   W( 49 ) = JVS( 593 )
   W( 51 ) = JVS( 594 )
   W( 52 ) = JVS( 595 )
   W( 53 ) = JVS( 596 )
   W( 54 ) = JVS( 597 )
   W( 55 ) = JVS( 598 )
   W( 56 ) = JVS( 599 )
   W( 57 ) = JVS( 600 )
   W( 58 ) = JVS( 601 )
   W( 59 ) = JVS( 602 )
  a = -W( 6 ) / JVS(           32  )
  W( 6 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 33 )
  a = -W( 33 ) / JVS(          159  )
  W( 33 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 160 )
  W( 44 ) = W( 44 ) + a*JVS( 161 )
  W( 46 ) = W( 46 ) + a*JVS( 162 )
  W( 47 ) = W( 47 ) + a*JVS( 163 )
  W( 49 ) = W( 49 ) + a*JVS( 164 )
  W( 52 ) = W( 52 ) + a*JVS( 165 )
  W( 53 ) = W( 53 ) + a*JVS( 166 )
  W( 54 ) = W( 54 ) + a*JVS( 167 )
  W( 55 ) = W( 55 ) + a*JVS( 168 )
  W( 56 ) = W( 56 ) + a*JVS( 169 )
  W( 58 ) = W( 58 ) + a*JVS( 170 )
  a = -W( 41 ) / JVS(          240  )
  W( 41 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 241 )
  W( 52 ) = W( 52 ) + a*JVS( 242 )
  W( 53 ) = W( 53 ) + a*JVS( 243 )
  W( 55 ) = W( 55 ) + a*JVS( 244 )
  W( 58 ) = W( 58 ) + a*JVS( 245 )
  a = -W( 44 ) / JVS(          300  )
  W( 44 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 301 )
  W( 52 ) = W( 52 ) + a*JVS( 302 )
  W( 53 ) = W( 53 ) + a*JVS( 303 )
  W( 55 ) = W( 55 ) + a*JVS( 304 )
  W( 58 ) = W( 58 ) + a*JVS( 305 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 325 )
  W( 53 ) = W( 53 ) + a*JVS( 326 )
  W( 54 ) = W( 54 ) + a*JVS( 327 )
  W( 58 ) = W( 58 ) + a*JVS( 328 )
  W( 59 ) = W( 59 ) + a*JVS( 329 )
  a = -W( 47 ) / JVS(          335  )
  W( 47 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 336 )
  W( 52 ) = W( 52 ) + a*JVS( 337 )
  W( 53 ) = W( 53 ) + a*JVS( 338 )
  W( 54 ) = W( 54 ) + a*JVS( 339 )
  W( 55 ) = W( 55 ) + a*JVS( 340 )
  W( 56 ) = W( 56 ) + a*JVS( 341 )
  W( 58 ) = W( 58 ) + a*JVS( 342 )
  W( 59 ) = W( 59 ) + a*JVS( 343 )
  a = -W( 49 ) / JVS(          371  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 372 )
  W( 52 ) = W( 52 ) + a*JVS( 373 )
  W( 53 ) = W( 53 ) + a*JVS( 374 )
  W( 54 ) = W( 54 ) + a*JVS( 375 )
  W( 55 ) = W( 55 ) + a*JVS( 376 )
  W( 58 ) = W( 58 ) + a*JVS( 377 )
  W( 59 ) = W( 59 ) + a*JVS( 378 )
  a = -W( 51 ) / JVS(          427  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 428 )
  W( 53 ) = W( 53 ) + a*JVS( 429 )
  W( 54 ) = W( 54 ) + a*JVS( 430 )
  W( 55 ) = W( 55 ) + a*JVS( 431 )
  W( 56 ) = W( 56 ) + a*JVS( 432 )
  W( 57 ) = W( 57 ) + a*JVS( 433 )
  W( 58 ) = W( 58 ) + a*JVS( 434 )
  W( 59 ) = W( 59 ) + a*JVS( 435 )
  a = -W( 52 ) / JVS(          457  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 458 )
  W( 54 ) = W( 54 ) + a*JVS( 459 )
  W( 55 ) = W( 55 ) + a*JVS( 460 )
  W( 56 ) = W( 56 ) + a*JVS( 461 )
  W( 57 ) = W( 57 ) + a*JVS( 462 )
  W( 58 ) = W( 58 ) + a*JVS( 463 )
  W( 59 ) = W( 59 ) + a*JVS( 464 )
  a = -W( 53 ) / JVS(          505  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 506 )
  W( 55 ) = W( 55 ) + a*JVS( 507 )
  W( 56 ) = W( 56 ) + a*JVS( 508 )
  W( 57 ) = W( 57 ) + a*JVS( 509 )
  W( 58 ) = W( 58 ) + a*JVS( 510 )
  W( 59 ) = W( 59 ) + a*JVS( 511 )
  a = -W( 54 ) / JVS(          536  )
  W( 54 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 537 )
  W( 56 ) = W( 56 ) + a*JVS( 538 )
  W( 57 ) = W( 57 ) + a*JVS( 539 )
  W( 58 ) = W( 58 ) + a*JVS( 540 )
  W( 59 ) = W( 59 ) + a*JVS( 541 )
  a = -W( 55 ) / JVS(          565  )
  W( 55 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 566 )
  W( 57 ) = W( 57 ) + a*JVS( 567 )
  W( 58 ) = W( 58 ) + a*JVS( 568 )
  W( 59 ) = W( 59 ) + a*JVS( 569 )
  a = -W( 56 ) / JVS(          583  )
  W( 56 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 584 )
  W( 58 ) = W( 58 ) + a*JVS( 585 )
  W( 59 ) = W( 59 ) + a*JVS( 586 )
  JVS( 587) = W( 6 )
  JVS( 588) = W( 33 )
  JVS( 589) = W( 41 )
  JVS( 590) = W( 44 )
  JVS( 591) = W( 46 )
  JVS( 592) = W( 47 )
  JVS( 593) = W( 49 )
  JVS( 594) = W( 51 )
  JVS( 595) = W( 52 )
  JVS( 596) = W( 53 )
  JVS( 597) = W( 54 )
  JVS( 598) = W( 55 )
  JVS( 599) = W( 56 )
  JVS( 600) = W( 57 )
  JVS( 601) = W( 58 )
  JVS( 602) = W( 59 )
  IF ( ABS(  JVS( 624 )) < TINY(a) ) THEN
         IER = 58                                      
         RETURN
  END IF
   W( 13 ) = JVS( 603 )
   W( 19 ) = JVS( 604 )
   W( 35 ) = JVS( 605 )
   W( 36 ) = JVS( 606 )
   W( 37 ) = JVS( 607 )
   W( 38 ) = JVS( 608 )
   W( 39 ) = JVS( 609 )
   W( 40 ) = JVS( 610 )
   W( 41 ) = JVS( 611 )
   W( 44 ) = JVS( 612 )
   W( 45 ) = JVS( 613 )
   W( 46 ) = JVS( 614 )
   W( 49 ) = JVS( 615 )
   W( 50 ) = JVS( 616 )
   W( 51 ) = JVS( 617 )
   W( 52 ) = JVS( 618 )
   W( 53 ) = JVS( 619 )
   W( 54 ) = JVS( 620 )
   W( 55 ) = JVS( 621 )
   W( 56 ) = JVS( 622 )
   W( 57 ) = JVS( 623 )
   W( 58 ) = JVS( 624 )
   W( 59 ) = JVS( 625 )
  a = -W( 13 ) / JVS(           47  )
  W( 13 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 48 )
  W( 58 ) = W( 58 ) + a*JVS( 49 )
  a = -W( 19 ) / JVS(           68  )
  W( 19 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 69 )
  W( 54 ) = W( 54 ) + a*JVS( 70 )
  W( 59 ) = W( 59 ) + a*JVS( 71 )
  a = -W( 35 ) / JVS(          182  )
  W( 35 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 183 )
  W( 52 ) = W( 52 ) + a*JVS( 184 )
  W( 53 ) = W( 53 ) + a*JVS( 185 )
  W( 55 ) = W( 55 ) + a*JVS( 186 )
  W( 58 ) = W( 58 ) + a*JVS( 187 )
  a = -W( 36 ) / JVS(          189  )
  W( 36 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 190 )
  W( 52 ) = W( 52 ) + a*JVS( 191 )
  W( 53 ) = W( 53 ) + a*JVS( 192 )
  W( 55 ) = W( 55 ) + a*JVS( 193 )
  W( 58 ) = W( 58 ) + a*JVS( 194 )
  a = -W( 37 ) / JVS(          197  )
  W( 37 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 198 )
  W( 51 ) = W( 51 ) + a*JVS( 199 )
  W( 52 ) = W( 52 ) + a*JVS( 200 )
  W( 53 ) = W( 53 ) + a*JVS( 201 )
  W( 54 ) = W( 54 ) + a*JVS( 202 )
  W( 55 ) = W( 55 ) + a*JVS( 203 )
  W( 58 ) = W( 58 ) + a*JVS( 204 )
  a = -W( 38 ) / JVS(          209  )
  W( 38 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 210 )
  W( 51 ) = W( 51 ) + a*JVS( 211 )
  W( 52 ) = W( 52 ) + a*JVS( 212 )
  W( 53 ) = W( 53 ) + a*JVS( 213 )
  W( 54 ) = W( 54 ) + a*JVS( 214 )
  W( 55 ) = W( 55 ) + a*JVS( 215 )
  W( 58 ) = W( 58 ) + a*JVS( 216 )
  a = -W( 39 ) / JVS(          221  )
  W( 39 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 222 )
  W( 51 ) = W( 51 ) + a*JVS( 223 )
  W( 52 ) = W( 52 ) + a*JVS( 224 )
  W( 53 ) = W( 53 ) + a*JVS( 225 )
  W( 54 ) = W( 54 ) + a*JVS( 226 )
  W( 55 ) = W( 55 ) + a*JVS( 227 )
  W( 58 ) = W( 58 ) + a*JVS( 228 )
  W( 59 ) = W( 59 ) + a*JVS( 229 )
  a = -W( 40 ) / JVS(          231  )
  W( 40 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 52 ) = W( 52 ) + a*JVS( 234 )
  W( 53 ) = W( 53 ) + a*JVS( 235 )
  W( 54 ) = W( 54 ) + a*JVS( 236 )
  W( 55 ) = W( 55 ) + a*JVS( 237 )
  W( 58 ) = W( 58 ) + a*JVS( 238 )
  a = -W( 41 ) / JVS(          240  )
  W( 41 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 241 )
  W( 52 ) = W( 52 ) + a*JVS( 242 )
  W( 53 ) = W( 53 ) + a*JVS( 243 )
  W( 55 ) = W( 55 ) + a*JVS( 244 )
  W( 58 ) = W( 58 ) + a*JVS( 245 )
  a = -W( 44 ) / JVS(          300  )
  W( 44 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 301 )
  W( 52 ) = W( 52 ) + a*JVS( 302 )
  W( 53 ) = W( 53 ) + a*JVS( 303 )
  W( 55 ) = W( 55 ) + a*JVS( 304 )
  W( 58 ) = W( 58 ) + a*JVS( 305 )
  a = -W( 45 ) / JVS(          311  )
  W( 45 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 312 )
  W( 52 ) = W( 52 ) + a*JVS( 313 )
  W( 53 ) = W( 53 ) + a*JVS( 314 )
  W( 54 ) = W( 54 ) + a*JVS( 315 )
  W( 55 ) = W( 55 ) + a*JVS( 316 )
  W( 58 ) = W( 58 ) + a*JVS( 317 )
  W( 59 ) = W( 59 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 325 )
  W( 53 ) = W( 53 ) + a*JVS( 326 )
  W( 54 ) = W( 54 ) + a*JVS( 327 )
  W( 58 ) = W( 58 ) + a*JVS( 328 )
  W( 59 ) = W( 59 ) + a*JVS( 329 )
  a = -W( 49 ) / JVS(          371  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 372 )
  W( 52 ) = W( 52 ) + a*JVS( 373 )
  W( 53 ) = W( 53 ) + a*JVS( 374 )
  W( 54 ) = W( 54 ) + a*JVS( 375 )
  W( 55 ) = W( 55 ) + a*JVS( 376 )
  W( 58 ) = W( 58 ) + a*JVS( 377 )
  W( 59 ) = W( 59 ) + a*JVS( 378 )
  a = -W( 50 ) / JVS(          385  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 386 )
  W( 52 ) = W( 52 ) + a*JVS( 387 )
  W( 53 ) = W( 53 ) + a*JVS( 388 )
  W( 54 ) = W( 54 ) + a*JVS( 389 )
  W( 55 ) = W( 55 ) + a*JVS( 390 )
  W( 56 ) = W( 56 ) + a*JVS( 391 )
  W( 58 ) = W( 58 ) + a*JVS( 392 )
  W( 59 ) = W( 59 ) + a*JVS( 393 )
  a = -W( 51 ) / JVS(          427  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 428 )
  W( 53 ) = W( 53 ) + a*JVS( 429 )
  W( 54 ) = W( 54 ) + a*JVS( 430 )
  W( 55 ) = W( 55 ) + a*JVS( 431 )
  W( 56 ) = W( 56 ) + a*JVS( 432 )
  W( 57 ) = W( 57 ) + a*JVS( 433 )
  W( 58 ) = W( 58 ) + a*JVS( 434 )
  W( 59 ) = W( 59 ) + a*JVS( 435 )
  a = -W( 52 ) / JVS(          457  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 458 )
  W( 54 ) = W( 54 ) + a*JVS( 459 )
  W( 55 ) = W( 55 ) + a*JVS( 460 )
  W( 56 ) = W( 56 ) + a*JVS( 461 )
  W( 57 ) = W( 57 ) + a*JVS( 462 )
  W( 58 ) = W( 58 ) + a*JVS( 463 )
  W( 59 ) = W( 59 ) + a*JVS( 464 )
  a = -W( 53 ) / JVS(          505  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 506 )
  W( 55 ) = W( 55 ) + a*JVS( 507 )
  W( 56 ) = W( 56 ) + a*JVS( 508 )
  W( 57 ) = W( 57 ) + a*JVS( 509 )
  W( 58 ) = W( 58 ) + a*JVS( 510 )
  W( 59 ) = W( 59 ) + a*JVS( 511 )
  a = -W( 54 ) / JVS(          536  )
  W( 54 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 537 )
  W( 56 ) = W( 56 ) + a*JVS( 538 )
  W( 57 ) = W( 57 ) + a*JVS( 539 )
  W( 58 ) = W( 58 ) + a*JVS( 540 )
  W( 59 ) = W( 59 ) + a*JVS( 541 )
  a = -W( 55 ) / JVS(          565  )
  W( 55 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 566 )
  W( 57 ) = W( 57 ) + a*JVS( 567 )
  W( 58 ) = W( 58 ) + a*JVS( 568 )
  W( 59 ) = W( 59 ) + a*JVS( 569 )
  a = -W( 56 ) / JVS(          583  )
  W( 56 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 584 )
  W( 58 ) = W( 58 ) + a*JVS( 585 )
  W( 59 ) = W( 59 ) + a*JVS( 586 )
  a = -W( 57 ) / JVS(          600  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 601 )
  W( 59 ) = W( 59 ) + a*JVS( 602 )
  JVS( 603) = W( 13 )
  JVS( 604) = W( 19 )
  JVS( 605) = W( 35 )
  JVS( 606) = W( 36 )
  JVS( 607) = W( 37 )
  JVS( 608) = W( 38 )
  JVS( 609) = W( 39 )
  JVS( 610) = W( 40 )
  JVS( 611) = W( 41 )
  JVS( 612) = W( 44 )
  JVS( 613) = W( 45 )
  JVS( 614) = W( 46 )
  JVS( 615) = W( 49 )
  JVS( 616) = W( 50 )
  JVS( 617) = W( 51 )
  JVS( 618) = W( 52 )
  JVS( 619) = W( 53 )
  JVS( 620) = W( 54 )
  JVS( 621) = W( 55 )
  JVS( 622) = W( 56 )
  JVS( 623) = W( 57 )
  JVS( 624) = W( 58 )
  JVS( 625) = W( 59 )
  IF ( ABS(  JVS( 659 )) < TINY(a) ) THEN
         IER = 59                                      
         RETURN
  END IF
   W( 12 ) = JVS( 626 )
   W( 15 ) = JVS( 627 )
   W( 19 ) = JVS( 628 )
   W( 20 ) = JVS( 629 )
   W( 23 ) = JVS( 630 )
   W( 25 ) = JVS( 631 )
   W( 30 ) = JVS( 632 )
   W( 31 ) = JVS( 633 )
   W( 32 ) = JVS( 634 )
   W( 34 ) = JVS( 635 )
   W( 35 ) = JVS( 636 )
   W( 36 ) = JVS( 637 )
   W( 37 ) = JVS( 638 )
   W( 38 ) = JVS( 639 )
   W( 39 ) = JVS( 640 )
   W( 40 ) = JVS( 641 )
   W( 41 ) = JVS( 642 )
   W( 43 ) = JVS( 643 )
   W( 44 ) = JVS( 644 )
   W( 45 ) = JVS( 645 )
   W( 46 ) = JVS( 646 )
   W( 47 ) = JVS( 647 )
   W( 48 ) = JVS( 648 )
   W( 49 ) = JVS( 649 )
   W( 50 ) = JVS( 650 )
   W( 51 ) = JVS( 651 )
   W( 52 ) = JVS( 652 )
   W( 53 ) = JVS( 653 )
   W( 54 ) = JVS( 654 )
   W( 55 ) = JVS( 655 )
   W( 56 ) = JVS( 656 )
   W( 57 ) = JVS( 657 )
   W( 58 ) = JVS( 658 )
   W( 59 ) = JVS( 659 )
  a = -W( 12 ) / JVS(           44  )
  W( 12 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 45 )
  W( 59 ) = W( 59 ) + a*JVS( 46 )
  a = -W( 15 ) / JVS(           53  )
  W( 15 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 54 )
  W( 59 ) = W( 59 ) + a*JVS( 55 )
  a = -W( 19 ) / JVS(           68  )
  W( 19 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 69 )
  W( 54 ) = W( 54 ) + a*JVS( 70 )
  W( 59 ) = W( 59 ) + a*JVS( 71 )
  a = -W( 20 ) / JVS(           72  )
  W( 20 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 73 )
  W( 53 ) = W( 53 ) + a*JVS( 74 )
  W( 59 ) = W( 59 ) + a*JVS( 75 )
  a = -W( 23 ) / JVS(           85  )
  W( 23 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 86 )
  W( 55 ) = W( 55 ) + a*JVS( 87 )
  W( 59 ) = W( 59 ) + a*JVS( 88 )
  a = -W( 25 ) / JVS(           95  )
  W( 25 ) = -a
  W( 30 ) = W( 30 ) + a*JVS( 96 )
  W( 31 ) = W( 31 ) + a*JVS( 97 )
  W( 34 ) = W( 34 ) + a*JVS( 98 )
  W( 43 ) = W( 43 ) + a*JVS( 99 )
  W( 48 ) = W( 48 ) + a*JVS( 100 )
  W( 51 ) = W( 51 ) + a*JVS( 101 )
  W( 53 ) = W( 53 ) + a*JVS( 102 )
  W( 54 ) = W( 54 ) + a*JVS( 103 )
  W( 59 ) = W( 59 ) + a*JVS( 104 )
  a = -W( 30 ) / JVS(          133  )
  W( 30 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 134 )
  W( 36 ) = W( 36 ) + a*JVS( 135 )
  W( 52 ) = W( 52 ) + a*JVS( 136 )
  W( 53 ) = W( 53 ) + a*JVS( 137 )
  W( 54 ) = W( 54 ) + a*JVS( 138 )
  W( 55 ) = W( 55 ) + a*JVS( 139 )
  W( 58 ) = W( 58 ) + a*JVS( 140 )
  a = -W( 31 ) / JVS(          141  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 142 )
  W( 45 ) = W( 45 ) + a*JVS( 143 )
  W( 52 ) = W( 52 ) + a*JVS( 144 )
  W( 53 ) = W( 53 ) + a*JVS( 145 )
  W( 54 ) = W( 54 ) + a*JVS( 146 )
  W( 55 ) = W( 55 ) + a*JVS( 147 )
  W( 58 ) = W( 58 ) + a*JVS( 148 )
  a = -W( 32 ) / JVS(          150  )
  W( 32 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 151 )
  W( 52 ) = W( 52 ) + a*JVS( 152 )
  W( 53 ) = W( 53 ) + a*JVS( 153 )
  W( 54 ) = W( 54 ) + a*JVS( 154 )
  W( 55 ) = W( 55 ) + a*JVS( 155 )
  W( 59 ) = W( 59 ) + a*JVS( 156 )
  a = -W( 34 ) / JVS(          171  )
  W( 34 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 172 )
  W( 36 ) = W( 36 ) + a*JVS( 173 )
  W( 45 ) = W( 45 ) + a*JVS( 174 )
  W( 50 ) = W( 50 ) + a*JVS( 175 )
  W( 52 ) = W( 52 ) + a*JVS( 176 )
  W( 53 ) = W( 53 ) + a*JVS( 177 )
  W( 54 ) = W( 54 ) + a*JVS( 178 )
  W( 55 ) = W( 55 ) + a*JVS( 179 )
  W( 58 ) = W( 58 ) + a*JVS( 180 )
  a = -W( 35 ) / JVS(          182  )
  W( 35 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 183 )
  W( 52 ) = W( 52 ) + a*JVS( 184 )
  W( 53 ) = W( 53 ) + a*JVS( 185 )
  W( 55 ) = W( 55 ) + a*JVS( 186 )
  W( 58 ) = W( 58 ) + a*JVS( 187 )
  a = -W( 36 ) / JVS(          189  )
  W( 36 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 190 )
  W( 52 ) = W( 52 ) + a*JVS( 191 )
  W( 53 ) = W( 53 ) + a*JVS( 192 )
  W( 55 ) = W( 55 ) + a*JVS( 193 )
  W( 58 ) = W( 58 ) + a*JVS( 194 )
  a = -W( 37 ) / JVS(          197  )
  W( 37 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 198 )
  W( 51 ) = W( 51 ) + a*JVS( 199 )
  W( 52 ) = W( 52 ) + a*JVS( 200 )
  W( 53 ) = W( 53 ) + a*JVS( 201 )
  W( 54 ) = W( 54 ) + a*JVS( 202 )
  W( 55 ) = W( 55 ) + a*JVS( 203 )
  W( 58 ) = W( 58 ) + a*JVS( 204 )
  a = -W( 38 ) / JVS(          209  )
  W( 38 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 210 )
  W( 51 ) = W( 51 ) + a*JVS( 211 )
  W( 52 ) = W( 52 ) + a*JVS( 212 )
  W( 53 ) = W( 53 ) + a*JVS( 213 )
  W( 54 ) = W( 54 ) + a*JVS( 214 )
  W( 55 ) = W( 55 ) + a*JVS( 215 )
  W( 58 ) = W( 58 ) + a*JVS( 216 )
  a = -W( 39 ) / JVS(          221  )
  W( 39 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 222 )
  W( 51 ) = W( 51 ) + a*JVS( 223 )
  W( 52 ) = W( 52 ) + a*JVS( 224 )
  W( 53 ) = W( 53 ) + a*JVS( 225 )
  W( 54 ) = W( 54 ) + a*JVS( 226 )
  W( 55 ) = W( 55 ) + a*JVS( 227 )
  W( 58 ) = W( 58 ) + a*JVS( 228 )
  W( 59 ) = W( 59 ) + a*JVS( 229 )
  a = -W( 40 ) / JVS(          231  )
  W( 40 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 52 ) = W( 52 ) + a*JVS( 234 )
  W( 53 ) = W( 53 ) + a*JVS( 235 )
  W( 54 ) = W( 54 ) + a*JVS( 236 )
  W( 55 ) = W( 55 ) + a*JVS( 237 )
  W( 58 ) = W( 58 ) + a*JVS( 238 )
  a = -W( 41 ) / JVS(          240  )
  W( 41 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 241 )
  W( 52 ) = W( 52 ) + a*JVS( 242 )
  W( 53 ) = W( 53 ) + a*JVS( 243 )
  W( 55 ) = W( 55 ) + a*JVS( 244 )
  W( 58 ) = W( 58 ) + a*JVS( 245 )
  a = -W( 43 ) / JVS(          284  )
  W( 43 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 285 )
  W( 45 ) = W( 45 ) + a*JVS( 286 )
  W( 46 ) = W( 46 ) + a*JVS( 287 )
  W( 49 ) = W( 49 ) + a*JVS( 288 )
  W( 50 ) = W( 50 ) + a*JVS( 289 )
  W( 51 ) = W( 51 ) + a*JVS( 290 )
  W( 52 ) = W( 52 ) + a*JVS( 291 )
  W( 53 ) = W( 53 ) + a*JVS( 292 )
  W( 54 ) = W( 54 ) + a*JVS( 293 )
  W( 55 ) = W( 55 ) + a*JVS( 294 )
  W( 56 ) = W( 56 ) + a*JVS( 295 )
  W( 57 ) = W( 57 ) + a*JVS( 296 )
  W( 58 ) = W( 58 ) + a*JVS( 297 )
  W( 59 ) = W( 59 ) + a*JVS( 298 )
  a = -W( 44 ) / JVS(          300  )
  W( 44 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 301 )
  W( 52 ) = W( 52 ) + a*JVS( 302 )
  W( 53 ) = W( 53 ) + a*JVS( 303 )
  W( 55 ) = W( 55 ) + a*JVS( 304 )
  W( 58 ) = W( 58 ) + a*JVS( 305 )
  a = -W( 45 ) / JVS(          311  )
  W( 45 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 312 )
  W( 52 ) = W( 52 ) + a*JVS( 313 )
  W( 53 ) = W( 53 ) + a*JVS( 314 )
  W( 54 ) = W( 54 ) + a*JVS( 315 )
  W( 55 ) = W( 55 ) + a*JVS( 316 )
  W( 58 ) = W( 58 ) + a*JVS( 317 )
  W( 59 ) = W( 59 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          324  )
  W( 46 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 325 )
  W( 53 ) = W( 53 ) + a*JVS( 326 )
  W( 54 ) = W( 54 ) + a*JVS( 327 )
  W( 58 ) = W( 58 ) + a*JVS( 328 )
  W( 59 ) = W( 59 ) + a*JVS( 329 )
  a = -W( 47 ) / JVS(          335  )
  W( 47 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 336 )
  W( 52 ) = W( 52 ) + a*JVS( 337 )
  W( 53 ) = W( 53 ) + a*JVS( 338 )
  W( 54 ) = W( 54 ) + a*JVS( 339 )
  W( 55 ) = W( 55 ) + a*JVS( 340 )
  W( 56 ) = W( 56 ) + a*JVS( 341 )
  W( 58 ) = W( 58 ) + a*JVS( 342 )
  W( 59 ) = W( 59 ) + a*JVS( 343 )
  a = -W( 48 ) / JVS(          357  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 358 )
  W( 50 ) = W( 50 ) + a*JVS( 359 )
  W( 51 ) = W( 51 ) + a*JVS( 360 )
  W( 52 ) = W( 52 ) + a*JVS( 361 )
  W( 53 ) = W( 53 ) + a*JVS( 362 )
  W( 54 ) = W( 54 ) + a*JVS( 363 )
  W( 55 ) = W( 55 ) + a*JVS( 364 )
  W( 56 ) = W( 56 ) + a*JVS( 365 )
  W( 57 ) = W( 57 ) + a*JVS( 366 )
  W( 58 ) = W( 58 ) + a*JVS( 367 )
  W( 59 ) = W( 59 ) + a*JVS( 368 )
  a = -W( 49 ) / JVS(          371  )
  W( 49 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 372 )
  W( 52 ) = W( 52 ) + a*JVS( 373 )
  W( 53 ) = W( 53 ) + a*JVS( 374 )
  W( 54 ) = W( 54 ) + a*JVS( 375 )
  W( 55 ) = W( 55 ) + a*JVS( 376 )
  W( 58 ) = W( 58 ) + a*JVS( 377 )
  W( 59 ) = W( 59 ) + a*JVS( 378 )
  a = -W( 50 ) / JVS(          385  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 386 )
  W( 52 ) = W( 52 ) + a*JVS( 387 )
  W( 53 ) = W( 53 ) + a*JVS( 388 )
  W( 54 ) = W( 54 ) + a*JVS( 389 )
  W( 55 ) = W( 55 ) + a*JVS( 390 )
  W( 56 ) = W( 56 ) + a*JVS( 391 )
  W( 58 ) = W( 58 ) + a*JVS( 392 )
  W( 59 ) = W( 59 ) + a*JVS( 393 )
  a = -W( 51 ) / JVS(          427  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 428 )
  W( 53 ) = W( 53 ) + a*JVS( 429 )
  W( 54 ) = W( 54 ) + a*JVS( 430 )
  W( 55 ) = W( 55 ) + a*JVS( 431 )
  W( 56 ) = W( 56 ) + a*JVS( 432 )
  W( 57 ) = W( 57 ) + a*JVS( 433 )
  W( 58 ) = W( 58 ) + a*JVS( 434 )
  W( 59 ) = W( 59 ) + a*JVS( 435 )
  a = -W( 52 ) / JVS(          457  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 458 )
  W( 54 ) = W( 54 ) + a*JVS( 459 )
  W( 55 ) = W( 55 ) + a*JVS( 460 )
  W( 56 ) = W( 56 ) + a*JVS( 461 )
  W( 57 ) = W( 57 ) + a*JVS( 462 )
  W( 58 ) = W( 58 ) + a*JVS( 463 )
  W( 59 ) = W( 59 ) + a*JVS( 464 )
  a = -W( 53 ) / JVS(          505  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 506 )
  W( 55 ) = W( 55 ) + a*JVS( 507 )
  W( 56 ) = W( 56 ) + a*JVS( 508 )
  W( 57 ) = W( 57 ) + a*JVS( 509 )
  W( 58 ) = W( 58 ) + a*JVS( 510 )
  W( 59 ) = W( 59 ) + a*JVS( 511 )
  a = -W( 54 ) / JVS(          536  )
  W( 54 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 537 )
  W( 56 ) = W( 56 ) + a*JVS( 538 )
  W( 57 ) = W( 57 ) + a*JVS( 539 )
  W( 58 ) = W( 58 ) + a*JVS( 540 )
  W( 59 ) = W( 59 ) + a*JVS( 541 )
  a = -W( 55 ) / JVS(          565  )
  W( 55 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 566 )
  W( 57 ) = W( 57 ) + a*JVS( 567 )
  W( 58 ) = W( 58 ) + a*JVS( 568 )
  W( 59 ) = W( 59 ) + a*JVS( 569 )
  a = -W( 56 ) / JVS(          583  )
  W( 56 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 584 )
  W( 58 ) = W( 58 ) + a*JVS( 585 )
  W( 59 ) = W( 59 ) + a*JVS( 586 )
  a = -W( 57 ) / JVS(          600  )
  W( 57 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 601 )
  W( 59 ) = W( 59 ) + a*JVS( 602 )
  a = -W( 58 ) / JVS(          624  )
  W( 58 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 625 )
  JVS( 626) = W( 12 )
  JVS( 627) = W( 15 )
  JVS( 628) = W( 19 )
  JVS( 629) = W( 20 )
  JVS( 630) = W( 23 )
  JVS( 631) = W( 25 )
  JVS( 632) = W( 30 )
  JVS( 633) = W( 31 )
  JVS( 634) = W( 32 )
  JVS( 635) = W( 34 )
  JVS( 636) = W( 35 )
  JVS( 637) = W( 36 )
  JVS( 638) = W( 37 )
  JVS( 639) = W( 38 )
  JVS( 640) = W( 39 )
  JVS( 641) = W( 40 )
  JVS( 642) = W( 41 )
  JVS( 643) = W( 43 )
  JVS( 644) = W( 44 )
  JVS( 645) = W( 45 )
  JVS( 646) = W( 46 )
  JVS( 647) = W( 47 )
  JVS( 648) = W( 48 )
  JVS( 649) = W( 49 )
  JVS( 650) = W( 50 )
  JVS( 651) = W( 51 )
  JVS( 652) = W( 52 )
  JVS( 653) = W( 53 )
  JVS( 654) = W( 54 )
  JVS( 655) = W( 55 )
  JVS( 656) = W( 56 )
  JVS( 657) = W( 57 )
  JVS( 658) = W( 58 )
  JVS( 659) = W( 59 )
   
   END SUBROUTINE decomp_radm2sorg
 


END MODULE radm2sorg_Integrator
