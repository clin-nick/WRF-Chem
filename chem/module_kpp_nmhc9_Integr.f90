
































MODULE nmhc9_Integrator

 USE nmhc9_Parameters
 USE nmhc9_Precision
 USE nmhc9_JacobianSP

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

SUBROUTINE  nmhc9_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE nmhc9_Parameters

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

   CALL nmhc9_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  nmhc9_INTEGRATE


SUBROUTINE  nmhc9_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE nmhc9_Parameters

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
      CALL nmhc9_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL nmhc9_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = nmhc9_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL nmhc9_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL nmhc9_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL nmhc9_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL nmhc9_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL nmhc9_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL nmhc9_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL nmhc9_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL nmhc9_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL nmhc9_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL nmhc9_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL nmhc9_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL nmhc9_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL nmhc9_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL nmhc9_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL nmhc9_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  nmhc9_ros_ErrorMsg(Code,T,H,IERR)



   USE nmhc9_Precision

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

 END SUBROUTINE  nmhc9_ros_ErrorMsg


 SUBROUTINE  nmhc9_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL nmhc9_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL nmhc9_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL nmhc9_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL nmhc9_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL nmhc9_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL nmhc9_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL nmhc9_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL nmhc9_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL nmhc9_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL nmhc9_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL nmhc9_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL nmhc9_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL nmhc9_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL nmhc9_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL nmhc9_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL nmhc9_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL nmhc9_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL nmhc9_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL nmhc9_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL nmhc9_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = nmhc9_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL nmhc9_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  nmhc9_ros_Integrator



  REAL(kind=dp) FUNCTION  nmhc9_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    nmhc9_ros_ErrorNorm = Err

  END FUNCTION  nmhc9_ros_ErrorNorm



  SUBROUTINE nmhc9_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL nmhc9_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL nmhc9_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL nmhc9_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  nmhc9_ros_FunTimeDeriv



  SUBROUTINE  nmhc9_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL nmhc9_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL nmhc9_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL nmhc9_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  nmhc9_ros_PrepareMatrix



  SUBROUTINE  nmhc9_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 



CALL decomp_nmhc9 ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  nmhc9_ros_Decomp



  SUBROUTINE  nmhc9_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL nmhc9_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  nmhc9_ros_Solve




  SUBROUTINE  nmhc9_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  nmhc9_Ros2



  SUBROUTINE  nmhc9_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  nmhc9_Ros3





  SUBROUTINE  nmhc9_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  nmhc9_Ros4


  SUBROUTINE  nmhc9_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  nmhc9_Rodas3


  SUBROUTINE  nmhc9_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  nmhc9_Rodas4




END SUBROUTINE  nmhc9_Rosenbrock




SUBROUTINE  nmhc9_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE nmhc9_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL nmhc9_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  nmhc9_FunTemplate



SUBROUTINE  nmhc9_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE nmhc9_Parameters
 
 USE nmhc9_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL nmhc9_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  nmhc9_JacTemplate

















SUBROUTINE nmhc9_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(53)
  A(2) = RCT(2)*V(21)*F(2)
  A(3) = RCT(3)*V(21)*F(1)
  A(4) = RCT(4)*F(2)
  A(5) = RCT(5)*V(53)*V(55)
  A(6) = RCT(6)*V(53)*V(57)
  A(7) = RCT(7)*V(55)*V(57)
  A(8) = RCT(8)*V(57)*V(57)
  A(9) = RCT(9)*V(8)
  A(10) = RCT(10)*V(8)*V(55)
  A(11) = RCT(11)*V(31)*V(55)
  A(12) = RCT(12)*V(24)*V(55)
  A(13) = RCT(13)*V(21)*V(24)
  A(14) = RCT(14)*V(50)*V(57)
  A(15) = RCT(15)*V(50)*V(57)
  A(16) = RCT(16)*V(50)*V(56)
  A(17) = RCT(17)*V(50)*V(50)
  A(18) = RCT(18)*V(50)*V(50)
  A(19) = RCT(19)*V(50)*V(52)
  A(20) = RCT(20)*V(12)
  A(21) = RCT(21)*V(12)*V(55)
  A(22) = RCT(22)*V(45)
  A(23) = RCT(23)*V(45)
  A(24) = RCT(24)*V(45)*V(55)
  A(25) = RCT(25)*V(45)*V(52)
  A(26) = RCT(26)*V(53)*V(56)
  A(27) = RCT(27)*V(56)*V(57)
  A(28) = RCT(28)*V(54)
  A(29) = RCT(29)*V(53)*V(54)
  A(30) = RCT(30)*V(54)*V(55)
  A(31) = RCT(31)*V(54)*V(57)
  A(32) = RCT(32)*V(25)
  A(33) = RCT(33)*V(25)*V(55)
  A(34) = RCT(34)*V(52)
  A(35) = RCT(35)*V(52)
  A(36) = RCT(36)*V(52)*V(56)
  A(37) = RCT(37)*V(52)*V(54)
  A(38) = RCT(38)*V(52)*V(57)
  A(39) = RCT(39)*V(5)
  A(40) = RCT(40)*V(5)
  A(41) = RCT(41)*V(5)
  A(42) = RCT(42)*V(5)*F(1)
  A(43) = RCT(43)*V(22)
  A(44) = RCT(44)*V(22)
  A(45) = RCT(45)*V(22)*V(55)
  A(46) = RCT(46)*V(55)*F(2)
  A(47) = RCT(47)*V(32)*V(55)
  A(48) = RCT(48)*V(50)*V(54)
  A(49) = RCT(49)*V(6)
  A(50) = RCT(50)*V(6)
  A(51) = RCT(51)*V(38)*V(55)
  A(52) = RCT(52)*V(38)*V(52)
  A(53) = RCT(53)*V(38)*V(53)
  A(54) = RCT(54)*V(48)*V(57)
  A(55) = RCT(55)*V(48)*V(56)
  A(56) = RCT(56)*V(11)
  A(57) = RCT(57)*V(11)*V(55)
  A(58) = RCT(58)*V(41)*V(55)
  A(59) = RCT(59)*V(41)
  A(60) = RCT(60)*V(41)*V(53)
  A(61) = RCT(61)*V(49)*V(56)
  A(62) = RCT(62)*V(49)*V(57)
  A(63) = RCT(63)*V(29)
  A(64) = RCT(64)*V(29)*V(55)
  A(65) = RCT(65)*V(33)*V(55)
  A(66) = RCT(66)*V(33)
  A(67) = RCT(67)*V(49)*V(54)
  A(68) = RCT(68)*V(48)*V(48)
  A(69) = RCT(69)*V(49)*V(49)
  A(70) = RCT(70)*V(48)*V(50)
  A(71) = RCT(71)*V(49)*V(50)
  A(72) = RCT(72)*V(35)*V(55)
  A(73) = RCT(73)*V(35)
  A(74) = RCT(74)*V(39)*V(55)
  A(75) = RCT(75)*V(39)
  A(76) = RCT(76)*V(51)*V(57)
  A(77) = RCT(77)*V(51)*V(57)
  A(78) = RCT(78)*V(51)*V(56)
  A(79) = RCT(79)*V(51)*V(54)
  A(80) = RCT(80)*V(50)*V(51)
  A(81) = RCT(81)*V(50)*V(51)
  A(82) = RCT(82)*V(51)*V(51)
  A(83) = RCT(83)*V(51)*V(52)
  A(84) = RCT(84)*V(9)
  A(85) = RCT(85)*V(9)*V(55)
  A(86) = RCT(86)*V(13)*V(55)
  A(87) = RCT(87)*V(13)
  A(88) = RCT(88)*V(13)
  A(89) = RCT(89)*V(14)*V(55)
  A(90) = RCT(90)*V(14)
  A(91) = RCT(91)*V(14)
  A(92) = RCT(92)*V(10)*V(55)
  A(93) = RCT(93)*V(15)*V(55)
  A(94) = RCT(94)*V(7)*V(55)
  A(95) = RCT(95)*V(7)
  A(96) = RCT(96)*V(1)*V(55)
  A(97) = RCT(97)*V(47)*V(57)
  A(98) = RCT(98)*V(47)*V(50)
  A(99) = RCT(99)*V(47)*V(51)
  A(100) = RCT(100)*V(47)*V(56)
  A(101) = RCT(101)*V(47)*V(52)
  A(102) = RCT(102)*V(46)*V(55)
  A(103) = RCT(103)*V(46)
  A(104) = RCT(104)*V(46)*V(52)
  A(105) = RCT(105)*V(16)
  A(106) = RCT(106)*V(16)*V(55)
  A(107) = RCT(107)*V(4)*V(55)
  A(108) = RCT(108)*V(34)*V(56)
  A(109) = RCT(109)*V(34)*V(57)
  A(110) = RCT(110)*V(34)*V(50)
  A(111) = RCT(111)*V(17)
  A(112) = RCT(112)*V(17)*V(55)
  A(113) = RCT(113)*V(18)
  A(114) = RCT(114)*V(18)*V(55)
  A(115) = RCT(115)*V(28)*V(55)
  A(116) = RCT(116)*V(28)
  A(117) = RCT(117)*V(36)*V(56)
  A(118) = RCT(118)*V(36)*V(57)
  A(119) = RCT(119)*V(36)*V(50)
  A(120) = RCT(120)*V(23)*V(55)
  A(121) = RCT(121)*V(23)
  A(122) = RCT(122)*V(37)*V(55)
  A(123) = RCT(123)*V(37)*V(53)
  A(124) = RCT(124)*V(37)*V(52)
  A(125) = RCT(125)*V(40)*V(56)
  A(126) = RCT(126)*V(40)*V(57)
  A(127) = RCT(127)*V(26)*V(55)
  A(128) = RCT(128)*V(26)
  A(129) = RCT(129)*V(19)*V(55)
  A(130) = RCT(130)*V(19)*V(53)
  A(131) = RCT(131)*V(2)*V(55)
  A(132) = RCT(132)*V(43)*V(56)
  A(133) = RCT(133)*V(43)*V(57)
  A(134) = RCT(134)*V(43)*V(50)
  A(135) = RCT(135)*V(27)*V(55)
  A(136) = RCT(136)*V(27)
  A(137) = RCT(137)*V(30)
  A(138) = RCT(138)*V(30)*V(55)
  A(139) = RCT(139)*V(42)*V(56)
  A(140) = RCT(140)*V(42)*V(57)
  A(141) = RCT(141)*V(20)*V(55)
  A(142) = RCT(142)*V(20)
  A(143) = RCT(143)*V(3)
  A(144) = RCT(144)*V(44)*V(55)
  A(145) = RCT(145)*V(44)


  Vdot(1) = -A(96)
  Vdot(2) = -A(131)
  Vdot(3) = 0.8*A(141)-A(143)
  Vdot(4) = -A(107)
  Vdot(5) = A(37)-A(39)-A(40)-A(41)-A(42)
  Vdot(6) = A(48)-A(49)-A(50)
  Vdot(7) = A(65)-A(94)-A(95)
  Vdot(8) = A(8)-A(9)-A(10)+0.09*A(53)
  Vdot(9) = A(76)-A(84)-A(85)
  Vdot(10) = A(77)+A(81)-A(92)+0.18*A(99)
  Vdot(11) = A(54)-A(56)-A(57)
  Vdot(12) = A(14)-A(20)-A(21)
  Vdot(13) = A(79)-A(86)-A(87)-A(88)
  Vdot(14) = A(67)-A(89)-A(90)-A(91)
  Vdot(15) = 0.28*A(53)+0.45*A(60)-A(93)+0.22*A(123)+0.54*A(130)
  Vdot(16) = A(97)-A(105)-A(106)
  Vdot(17) = A(109)-A(111)-A(112)
  Vdot(18) = 0.04*A(108)-A(113)-A(114)
  Vdot(19) = -A(129)-A(130)
  Vdot(20) = A(140)-A(141)-A(142)
  Vdot(21) = A(1)-A(2)-A(3)-A(13)
  Vdot(22) = A(31)-A(43)-A(44)-A(45)
  Vdot(23) = A(118)-A(120)-A(121)
  Vdot(24) = -A(12)-A(13)+0.06*A(123)
  Vdot(25) = A(25)+A(30)-A(32)-A(33)+2*A(41)+2*A(42)+A(104)
  Vdot(26) = A(126)-A(127)-A(128)
  Vdot(27) = A(133)-A(135)-A(136)
  Vdot(28) = 0.96*A(108)+A(110)+A(111)+0.7*A(112)+A(113)+A(114)-A(115)-A(116)
  Vdot(29) = A(62)-A(63)-A(64)
  Vdot(30) = 0.56*A(132)+0.88*A(134)+0.85*A(135)+0.67*A(136)-A(137)-A(138)+A(144)+0.67*A(145)
  Vdot(31) = -A(11)+A(22)+A(23)+A(24)+A(25)+0.14*A(53)+A(59)+0.22*A(60)+0.25*A(61)+0.25*A(63)+0.5*A(69)+0.125*A(71)&
               &+A(72)+A(73)+A(94)+A(95)+A(103)+0.31*A(123)+0.23*A(130)
  Vdot(32) = A(18)-A(47)+0.25*A(70)+0.25*A(71)+0.25*A(98)+0.2*A(110)+0.5*A(119)+0.03*A(123)+0.18*A(134)
  Vdot(33) = A(52)+0.12*A(55)-A(65)-A(66)
  Vdot(34) = 0.82*A(107)-A(108)-A(109)-A(110)+0.3*A(112)
  Vdot(35) = 0.9*A(60)+0.5*A(61)+0.5*A(63)+A(69)+0.25*A(70)+0.5*A(71)-A(72)-A(73)+A(74)+0.5*A(119)+0.7*A(120)+0.04&
               &*A(123)
  Vdot(36) = A(115)-A(117)-A(118)-A(119)+0.3*A(120)
  Vdot(37) = -A(122)-A(123)-A(124)
  Vdot(38) = -A(51)-A(52)-A(53)
  Vdot(39) = 0.25*A(61)+0.25*A(63)+A(65)+A(69)+0.25*A(70)+0.375*A(71)-A(74)-A(75)+A(89)+A(91)+0.2*A(119)+0.5*A(127)
  Vdot(40) = A(122)-A(125)-A(126)+0.5*A(127)+0.666667*A(129)
  Vdot(41) = 0.65*A(53)+0.88*A(55)+A(56)+A(57)-A(58)-A(59)-A(60)+A(66)+2*A(68)+0.5*A(70)
  Vdot(42) = A(138)-A(139)-A(140)+0.2*A(141)
  Vdot(43) = A(131)-A(132)-A(133)-A(134)+0.15*A(135)
  Vdot(44) = A(124)+0.02*A(125)+0.16*A(132)+0.015*A(139)-A(144)-A(145)
  Vdot(45) = 0.25*A(13)+A(15)+A(16)+2*A(17)+A(18)+A(19)+A(20)+0.3*A(21)-A(22)-A(23)-A(24)-A(25)+A(47)+0.39*A(50)+0.58&
               &*A(53)+0.88*A(55)+A(56)+A(59)+0.75*A(61)+0.75*A(63)+A(66)+A(68)+0.5*A(69)+1.25*A(70)+1.125*A(71)+A(75)+A(80)&
               &+A(81)+A(86)+A(94)+A(95)+0.75*A(98)+0.8*A(110)+A(117)+0.8*A(119)+0.57*A(123)+0.98*A(125)+A(128)+A(130)+0.68&
               &*A(134)
  Vdot(46) = 0.75*A(98)+A(99)+A(100)+A(101)-A(102)-A(103)-A(104)+A(105)+0.7*A(106)+0.47*A(123)+0.98*A(125)+A(128)+0.84&
               &*A(132)+0.12*A(134)+0.33*A(136)+0.985*A(139)+A(142)+0.33*A(145)
  Vdot(47) = A(96)-A(97)-A(98)-A(99)-A(100)-A(101)+0.3*A(106)+0.18*A(107)+0.06*A(123)+0.28*A(132)+0.12*A(134)+0.33&
               &*A(136)+A(137)+0.33*A(145)
  Vdot(48) = A(51)-A(54)-A(55)-2*A(68)-A(70)
  Vdot(49) = 0.1*A(53)+A(58)-A(61)-A(62)+A(64)-A(67)-2*A(69)-A(71)+A(90)
  Vdot(50) = A(12)+0.75*A(13)-A(14)-A(15)-A(16)-2*A(17)-2*A(18)-A(19)+0.7*A(21)-A(48)+A(49)+0.61*A(50)+0.08*A(53)-A(70)&
               &-A(71)+A(78)-A(81)+2*A(82)+A(83)+A(84)+A(92)-A(98)+0.82*A(99)+A(103)-A(110)+A(116)-A(119)+0.07*A(123)-A(134)
  Vdot(51) = 0.1*A(53)+A(59)+0.1*A(60)+0.25*A(61)+0.25*A(63)+0.125*A(71)+A(72)+A(73)+A(75)-A(76)-A(77)-A(78)-A(79)-A(80)&
               &-A(81)-2*A(82)-A(83)+A(85)+A(87)+A(88)-A(99)+A(102)+A(104)+A(116)+A(117)+0.3*A(119)+A(121)+0.23*A(123)&
               &+A(137)+0.985*A(139)+A(142)+2*A(143)
  Vdot(52) = -A(19)-A(25)+A(29)+A(33)-A(34)-A(35)-A(36)-A(37)-A(38)+A(39)+A(40)+0.39*A(43)+0.39*A(50)-A(52)-A(83)-A(101)&
               &-A(104)-A(124)
  Vdot(53) = -A(1)+A(2)+2*A(4)-A(5)-A(6)-A(26)+A(28)-A(29)+A(34)-A(53)-A(60)+A(77)-A(123)-A(130)
  Vdot(54) = A(16)+A(19)+A(26)+A(27)-A(28)-A(29)-A(30)-A(31)+A(32)+A(34)+2*A(36)-A(37)+A(38)+A(39)+A(40)+0.61*A(43)&
               &+A(44)+A(45)-A(48)+A(49)+0.61*A(50)+0.88*A(55)+A(61)+A(66)-A(67)+A(78)-A(79)+A(83)+A(86)+A(87)+A(88)+A(89)&
               &+A(90)+A(91)+A(94)+A(95)+A(100)+A(101)+0.96*A(108)+A(113)+A(114)+A(117)+0.98*A(125)+0.84*A(132)+0.985*A(139)&
               &+A(144)+A(145)
  Vdot(55) = 2*A(3)-A(5)+A(6)-A(7)+2*A(9)-A(10)-A(11)-A(12)+0.75*A(13)+A(20)-0.7*A(21)-A(24)+A(27)-A(30)+A(32)-A(33)&
               &+A(38)+0.39*A(43)-A(45)-A(46)-A(47)-A(51)+0.25*A(53)+A(56)-A(58)+0.19*A(60)+A(63)-A(64)-A(65)-A(72)-A(74)&
               &+A(84)-A(85)-A(86)-A(89)-A(92)-A(93)-A(94)-A(96)-A(102)+A(105)-0.3*A(106)-A(107)+A(111)-0.3*A(112)-A(114)&
               &-A(115)-0.3*A(120)+A(121)-A(122)+0.33*A(123)-0.5*A(127)+A(128)-A(129)+0.12*A(130)-A(131)-0.15*A(135)+A(136)&
               &-A(138)-0.2*A(141)+A(142)-A(144)
  Vdot(56) = -A(16)-A(26)-A(27)+A(28)+A(35)-A(36)-A(55)-A(61)-A(78)-A(100)-A(108)-A(117)-A(125)-A(132)-A(139)
  Vdot(57) = A(5)-A(6)-A(7)-2*A(8)+A(10)+A(11)+0.4*A(13)-A(14)-A(15)+A(16)+2*A(17)+A(19)+A(20)+2*A(22)+A(24)+A(25)-A(27)&
               &-A(31)-A(38)+0.61*A(43)+A(44)+A(46)+A(47)+0.25*A(53)-A(54)+0.88*A(55)+A(56)+A(59)+0.32*A(60)+0.75*A(61)&
               &-A(62)+0.75*A(63)+A(66)+A(68)+A(69)+A(70)+0.875*A(71)+A(73)+A(74)+A(75)-A(76)-A(77)+A(80)+A(93)-A(97)+A(98)&
               &+0.82*A(99)+A(100)+A(101)+A(103)+A(105)+0.96*A(108)-A(109)+0.8*A(110)+A(111)+A(113)-A(118)+0.3*A(119)+A(121)&
               &+0.26*A(123)+0.98*A(125)-A(126)+A(128)+0.22*A(130)+0.56*A(132)-A(133)+1.23*A(134)+0.67*A(136)-A(140)+0.67&
               &*A(145)
      
END SUBROUTINE nmhc9_Fun
















SUBROUTINE nmhc9_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(53)
  IRR(2) = RCT(2)*V(21)*F(2)
  IRR(3) = RCT(3)*V(21)*F(1)
  IRR(4) = RCT(4)*F(2)
  IRR(5) = RCT(5)*V(53)*V(55)
  IRR(6) = RCT(6)*V(53)*V(57)
  IRR(7) = RCT(7)*V(55)*V(57)
  IRR(8) = RCT(8)*V(57)*V(57)
  IRR(9) = RCT(9)*V(8)
  IRR(10) = RCT(10)*V(8)*V(55)
  IRR(11) = RCT(11)*V(31)*V(55)
  IRR(12) = RCT(12)*V(24)*V(55)
  IRR(13) = RCT(13)*V(21)*V(24)
  IRR(14) = RCT(14)*V(50)*V(57)
  IRR(15) = RCT(15)*V(50)*V(57)
  IRR(16) = RCT(16)*V(50)*V(56)
  IRR(17) = RCT(17)*V(50)*V(50)
  IRR(18) = RCT(18)*V(50)*V(50)
  IRR(19) = RCT(19)*V(50)*V(52)
  IRR(20) = RCT(20)*V(12)
  IRR(21) = RCT(21)*V(12)*V(55)
  IRR(22) = RCT(22)*V(45)
  IRR(23) = RCT(23)*V(45)
  IRR(24) = RCT(24)*V(45)*V(55)
  IRR(25) = RCT(25)*V(45)*V(52)
  IRR(26) = RCT(26)*V(53)*V(56)
  IRR(27) = RCT(27)*V(56)*V(57)
  IRR(28) = RCT(28)*V(54)
  IRR(29) = RCT(29)*V(53)*V(54)
  IRR(30) = RCT(30)*V(54)*V(55)
  IRR(31) = RCT(31)*V(54)*V(57)
  IRR(32) = RCT(32)*V(25)
  IRR(33) = RCT(33)*V(25)*V(55)
  IRR(34) = RCT(34)*V(52)
  IRR(35) = RCT(35)*V(52)
  IRR(36) = RCT(36)*V(52)*V(56)
  IRR(37) = RCT(37)*V(52)*V(54)
  IRR(38) = RCT(38)*V(52)*V(57)
  IRR(39) = RCT(39)*V(5)
  IRR(40) = RCT(40)*V(5)
  IRR(41) = RCT(41)*V(5)
  IRR(42) = RCT(42)*V(5)*F(1)
  IRR(43) = RCT(43)*V(22)
  IRR(44) = RCT(44)*V(22)
  IRR(45) = RCT(45)*V(22)*V(55)
  IRR(46) = RCT(46)*V(55)*F(2)
  IRR(47) = RCT(47)*V(32)*V(55)
  IRR(48) = RCT(48)*V(50)*V(54)
  IRR(49) = RCT(49)*V(6)
  IRR(50) = RCT(50)*V(6)
  IRR(51) = RCT(51)*V(38)*V(55)
  IRR(52) = RCT(52)*V(38)*V(52)
  IRR(53) = RCT(53)*V(38)*V(53)
  IRR(54) = RCT(54)*V(48)*V(57)
  IRR(55) = RCT(55)*V(48)*V(56)
  IRR(56) = RCT(56)*V(11)
  IRR(57) = RCT(57)*V(11)*V(55)
  IRR(58) = RCT(58)*V(41)*V(55)
  IRR(59) = RCT(59)*V(41)
  IRR(60) = RCT(60)*V(41)*V(53)
  IRR(61) = RCT(61)*V(49)*V(56)
  IRR(62) = RCT(62)*V(49)*V(57)
  IRR(63) = RCT(63)*V(29)
  IRR(64) = RCT(64)*V(29)*V(55)
  IRR(65) = RCT(65)*V(33)*V(55)
  IRR(66) = RCT(66)*V(33)
  IRR(67) = RCT(67)*V(49)*V(54)
  IRR(68) = RCT(68)*V(48)*V(48)
  IRR(69) = RCT(69)*V(49)*V(49)
  IRR(70) = RCT(70)*V(48)*V(50)
  IRR(71) = RCT(71)*V(49)*V(50)
  IRR(72) = RCT(72)*V(35)*V(55)
  IRR(73) = RCT(73)*V(35)
  IRR(74) = RCT(74)*V(39)*V(55)
  IRR(75) = RCT(75)*V(39)
  IRR(76) = RCT(76)*V(51)*V(57)
  IRR(77) = RCT(77)*V(51)*V(57)
  IRR(78) = RCT(78)*V(51)*V(56)
  IRR(79) = RCT(79)*V(51)*V(54)
  IRR(80) = RCT(80)*V(50)*V(51)
  IRR(81) = RCT(81)*V(50)*V(51)
  IRR(82) = RCT(82)*V(51)*V(51)
  IRR(83) = RCT(83)*V(51)*V(52)
  IRR(84) = RCT(84)*V(9)
  IRR(85) = RCT(85)*V(9)*V(55)
  IRR(86) = RCT(86)*V(13)*V(55)
  IRR(87) = RCT(87)*V(13)
  IRR(88) = RCT(88)*V(13)
  IRR(89) = RCT(89)*V(14)*V(55)
  IRR(90) = RCT(90)*V(14)
  IRR(91) = RCT(91)*V(14)
  IRR(92) = RCT(92)*V(10)*V(55)
  IRR(93) = RCT(93)*V(15)*V(55)
  IRR(94) = RCT(94)*V(7)*V(55)
  IRR(95) = RCT(95)*V(7)
  IRR(96) = RCT(96)*V(1)*V(55)
  IRR(97) = RCT(97)*V(47)*V(57)
  IRR(98) = RCT(98)*V(47)*V(50)
  IRR(99) = RCT(99)*V(47)*V(51)
  IRR(100) = RCT(100)*V(47)*V(56)
  IRR(101) = RCT(101)*V(47)*V(52)
  IRR(102) = RCT(102)*V(46)*V(55)
  IRR(103) = RCT(103)*V(46)
  IRR(104) = RCT(104)*V(46)*V(52)
  IRR(105) = RCT(105)*V(16)
  IRR(106) = RCT(106)*V(16)*V(55)
  IRR(107) = RCT(107)*V(4)*V(55)
  IRR(108) = RCT(108)*V(34)*V(56)
  IRR(109) = RCT(109)*V(34)*V(57)
  IRR(110) = RCT(110)*V(34)*V(50)
  IRR(111) = RCT(111)*V(17)
  IRR(112) = RCT(112)*V(17)*V(55)
  IRR(113) = RCT(113)*V(18)
  IRR(114) = RCT(114)*V(18)*V(55)
  IRR(115) = RCT(115)*V(28)*V(55)
  IRR(116) = RCT(116)*V(28)
  IRR(117) = RCT(117)*V(36)*V(56)
  IRR(118) = RCT(118)*V(36)*V(57)
  IRR(119) = RCT(119)*V(36)*V(50)
  IRR(120) = RCT(120)*V(23)*V(55)
  IRR(121) = RCT(121)*V(23)
  IRR(122) = RCT(122)*V(37)*V(55)
  IRR(123) = RCT(123)*V(37)*V(53)
  IRR(124) = RCT(124)*V(37)*V(52)
  IRR(125) = RCT(125)*V(40)*V(56)
  IRR(126) = RCT(126)*V(40)*V(57)
  IRR(127) = RCT(127)*V(26)*V(55)
  IRR(128) = RCT(128)*V(26)
  IRR(129) = RCT(129)*V(19)*V(55)
  IRR(130) = RCT(130)*V(19)*V(53)
  IRR(131) = RCT(131)*V(2)*V(55)
  IRR(132) = RCT(132)*V(43)*V(56)
  IRR(133) = RCT(133)*V(43)*V(57)
  IRR(134) = RCT(134)*V(43)*V(50)
  IRR(135) = RCT(135)*V(27)*V(55)
  IRR(136) = RCT(136)*V(27)
  IRR(137) = RCT(137)*V(30)
  IRR(138) = RCT(138)*V(30)*V(55)
  IRR(139) = RCT(139)*V(42)*V(56)
  IRR(140) = RCT(140)*V(42)*V(57)
  IRR(141) = RCT(141)*V(20)*V(55)
  IRR(142) = RCT(142)*V(20)
  IRR(143) = RCT(143)*V(3)
  IRR(144) = RCT(144)*V(44)*V(55)
  IRR(145) = RCT(145)*V(44)
      
END SUBROUTINE nmhc9_IRRFun
















SUBROUTINE nmhc9_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(243)


  B(1) = RCT(1)

  B(2) = RCT(2)*F(2)

  B(4) = RCT(3)*F(1)

  B(7) = RCT(5)*V(55)

  B(8) = RCT(5)*V(53)

  B(9) = RCT(6)*V(57)

  B(10) = RCT(6)*V(53)

  B(11) = RCT(7)*V(57)

  B(12) = RCT(7)*V(55)

  B(13) = RCT(8)*2*V(57)

  B(14) = RCT(9)

  B(15) = RCT(10)*V(55)

  B(16) = RCT(10)*V(8)

  B(17) = RCT(11)*V(55)

  B(18) = RCT(11)*V(31)

  B(19) = RCT(12)*V(55)

  B(20) = RCT(12)*V(24)

  B(21) = RCT(13)*V(24)

  B(22) = RCT(13)*V(21)

  B(23) = RCT(14)*V(57)

  B(24) = RCT(14)*V(50)

  B(25) = RCT(15)*V(57)

  B(26) = RCT(15)*V(50)

  B(27) = RCT(16)*V(56)

  B(28) = RCT(16)*V(50)

  B(29) = RCT(17)*2*V(50)

  B(30) = RCT(18)*2*V(50)

  B(31) = RCT(19)*V(52)

  B(32) = RCT(19)*V(50)

  B(33) = RCT(20)

  B(34) = RCT(21)*V(55)

  B(35) = RCT(21)*V(12)

  B(36) = RCT(22)

  B(37) = RCT(23)

  B(38) = RCT(24)*V(55)

  B(39) = RCT(24)*V(45)

  B(40) = RCT(25)*V(52)

  B(41) = RCT(25)*V(45)

  B(42) = RCT(26)*V(56)

  B(43) = RCT(26)*V(53)

  B(44) = RCT(27)*V(57)

  B(45) = RCT(27)*V(56)

  B(46) = RCT(28)

  B(47) = RCT(29)*V(54)

  B(48) = RCT(29)*V(53)

  B(49) = RCT(30)*V(55)

  B(50) = RCT(30)*V(54)

  B(51) = RCT(31)*V(57)

  B(52) = RCT(31)*V(54)

  B(53) = RCT(32)

  B(54) = RCT(33)*V(55)

  B(55) = RCT(33)*V(25)

  B(56) = RCT(34)

  B(57) = RCT(35)

  B(58) = RCT(36)*V(56)

  B(59) = RCT(36)*V(52)

  B(60) = RCT(37)*V(54)

  B(61) = RCT(37)*V(52)

  B(62) = RCT(38)*V(57)

  B(63) = RCT(38)*V(52)

  B(64) = RCT(39)

  B(65) = RCT(40)

  B(66) = RCT(41)

  B(67) = RCT(42)*F(1)

  B(69) = RCT(43)

  B(70) = RCT(44)

  B(71) = RCT(45)*V(55)

  B(72) = RCT(45)*V(22)

  B(73) = RCT(46)*F(2)

  B(75) = RCT(47)*V(55)

  B(76) = RCT(47)*V(32)

  B(77) = RCT(48)*V(54)

  B(78) = RCT(48)*V(50)

  B(79) = RCT(49)

  B(80) = RCT(50)

  B(81) = RCT(51)*V(55)

  B(82) = RCT(51)*V(38)

  B(83) = RCT(52)*V(52)

  B(84) = RCT(52)*V(38)

  B(85) = RCT(53)*V(53)

  B(86) = RCT(53)*V(38)

  B(87) = RCT(54)*V(57)

  B(88) = RCT(54)*V(48)

  B(89) = RCT(55)*V(56)

  B(90) = RCT(55)*V(48)

  B(91) = RCT(56)

  B(92) = RCT(57)*V(55)

  B(93) = RCT(57)*V(11)

  B(94) = RCT(58)*V(55)

  B(95) = RCT(58)*V(41)

  B(96) = RCT(59)

  B(97) = RCT(60)*V(53)

  B(98) = RCT(60)*V(41)

  B(99) = RCT(61)*V(56)

  B(100) = RCT(61)*V(49)

  B(101) = RCT(62)*V(57)

  B(102) = RCT(62)*V(49)

  B(103) = RCT(63)

  B(104) = RCT(64)*V(55)

  B(105) = RCT(64)*V(29)

  B(106) = RCT(65)*V(55)

  B(107) = RCT(65)*V(33)

  B(108) = RCT(66)

  B(109) = RCT(67)*V(54)

  B(110) = RCT(67)*V(49)

  B(111) = RCT(68)*2*V(48)

  B(112) = RCT(69)*2*V(49)

  B(113) = RCT(70)*V(50)

  B(114) = RCT(70)*V(48)

  B(115) = RCT(71)*V(50)

  B(116) = RCT(71)*V(49)

  B(117) = RCT(72)*V(55)

  B(118) = RCT(72)*V(35)

  B(119) = RCT(73)

  B(120) = RCT(74)*V(55)

  B(121) = RCT(74)*V(39)

  B(122) = RCT(75)

  B(123) = RCT(76)*V(57)

  B(124) = RCT(76)*V(51)

  B(125) = RCT(77)*V(57)

  B(126) = RCT(77)*V(51)

  B(127) = RCT(78)*V(56)

  B(128) = RCT(78)*V(51)

  B(129) = RCT(79)*V(54)

  B(130) = RCT(79)*V(51)

  B(131) = RCT(80)*V(51)

  B(132) = RCT(80)*V(50)

  B(133) = RCT(81)*V(51)

  B(134) = RCT(81)*V(50)

  B(135) = RCT(82)*2*V(51)

  B(136) = RCT(83)*V(52)

  B(137) = RCT(83)*V(51)

  B(138) = RCT(84)

  B(139) = RCT(85)*V(55)

  B(140) = RCT(85)*V(9)

  B(141) = RCT(86)*V(55)

  B(142) = RCT(86)*V(13)

  B(143) = RCT(87)

  B(144) = RCT(88)

  B(145) = RCT(89)*V(55)

  B(146) = RCT(89)*V(14)

  B(147) = RCT(90)

  B(148) = RCT(91)

  B(149) = RCT(92)*V(55)

  B(150) = RCT(92)*V(10)

  B(151) = RCT(93)*V(55)

  B(152) = RCT(93)*V(15)

  B(153) = RCT(94)*V(55)

  B(154) = RCT(94)*V(7)

  B(155) = RCT(95)

  B(156) = RCT(96)*V(55)

  B(157) = RCT(96)*V(1)

  B(158) = RCT(97)*V(57)

  B(159) = RCT(97)*V(47)

  B(160) = RCT(98)*V(50)

  B(161) = RCT(98)*V(47)

  B(162) = RCT(99)*V(51)

  B(163) = RCT(99)*V(47)

  B(164) = RCT(100)*V(56)

  B(165) = RCT(100)*V(47)

  B(166) = RCT(101)*V(52)

  B(167) = RCT(101)*V(47)

  B(168) = RCT(102)*V(55)

  B(169) = RCT(102)*V(46)

  B(170) = RCT(103)

  B(171) = RCT(104)*V(52)

  B(172) = RCT(104)*V(46)

  B(173) = RCT(105)

  B(174) = RCT(106)*V(55)

  B(175) = RCT(106)*V(16)

  B(176) = RCT(107)*V(55)

  B(177) = RCT(107)*V(4)

  B(178) = RCT(108)*V(56)

  B(179) = RCT(108)*V(34)

  B(180) = RCT(109)*V(57)

  B(181) = RCT(109)*V(34)

  B(182) = RCT(110)*V(50)

  B(183) = RCT(110)*V(34)

  B(184) = RCT(111)

  B(185) = RCT(112)*V(55)

  B(186) = RCT(112)*V(17)

  B(187) = RCT(113)

  B(188) = RCT(114)*V(55)

  B(189) = RCT(114)*V(18)

  B(190) = RCT(115)*V(55)

  B(191) = RCT(115)*V(28)

  B(192) = RCT(116)

  B(193) = RCT(117)*V(56)

  B(194) = RCT(117)*V(36)

  B(195) = RCT(118)*V(57)

  B(196) = RCT(118)*V(36)

  B(197) = RCT(119)*V(50)

  B(198) = RCT(119)*V(36)

  B(199) = RCT(120)*V(55)

  B(200) = RCT(120)*V(23)

  B(201) = RCT(121)

  B(202) = RCT(122)*V(55)

  B(203) = RCT(122)*V(37)

  B(204) = RCT(123)*V(53)

  B(205) = RCT(123)*V(37)

  B(206) = RCT(124)*V(52)

  B(207) = RCT(124)*V(37)

  B(208) = RCT(125)*V(56)

  B(209) = RCT(125)*V(40)

  B(210) = RCT(126)*V(57)

  B(211) = RCT(126)*V(40)

  B(212) = RCT(127)*V(55)

  B(213) = RCT(127)*V(26)

  B(214) = RCT(128)

  B(215) = RCT(129)*V(55)

  B(216) = RCT(129)*V(19)

  B(217) = RCT(130)*V(53)

  B(218) = RCT(130)*V(19)

  B(219) = RCT(131)*V(55)

  B(220) = RCT(131)*V(2)

  B(221) = RCT(132)*V(56)

  B(222) = RCT(132)*V(43)

  B(223) = RCT(133)*V(57)

  B(224) = RCT(133)*V(43)

  B(225) = RCT(134)*V(50)

  B(226) = RCT(134)*V(43)

  B(227) = RCT(135)*V(55)

  B(228) = RCT(135)*V(27)

  B(229) = RCT(136)

  B(230) = RCT(137)

  B(231) = RCT(138)*V(55)

  B(232) = RCT(138)*V(30)

  B(233) = RCT(139)*V(56)

  B(234) = RCT(139)*V(42)

  B(235) = RCT(140)*V(57)

  B(236) = RCT(140)*V(42)

  B(237) = RCT(141)*V(55)

  B(238) = RCT(141)*V(20)

  B(239) = RCT(142)

  B(240) = RCT(143)

  B(241) = RCT(144)*V(55)

  B(242) = RCT(144)*V(44)

  B(243) = RCT(145)



  JVS(1) = -B(156)

  JVS(2) = -B(157)

  JVS(3) = -B(219)

  JVS(4) = -B(220)

  JVS(5) = -B(240)

  JVS(6) = 0.8*B(237)

  JVS(7) = 0.8*B(238)

  JVS(8) = -B(176)

  JVS(9) = -B(177)

  JVS(10) = -B(64)-B(65)-B(66)-B(67)

  JVS(11) = B(60)

  JVS(12) = B(61)

  JVS(13) = -B(79)-B(80)

  JVS(14) = B(77)

  JVS(15) = B(78)

  JVS(16) = -B(153)-B(155)

  JVS(17) = B(106)

  JVS(18) = B(107)-B(154)

  JVS(19) = -B(14)-B(15)

  JVS(20) = 0.09*B(85)

  JVS(21) = 0.09*B(86)

  JVS(22) = -B(16)

  JVS(23) = B(13)

  JVS(24) = -B(138)-B(139)

  JVS(25) = B(123)

  JVS(26) = -B(140)

  JVS(27) = B(124)

  JVS(28) = -B(149)

  JVS(29) = 0.18*B(162)

  JVS(30) = B(133)

  JVS(31) = B(125)+B(134)+0.18*B(163)

  JVS(32) = -B(150)

  JVS(33) = B(126)

  JVS(34) = -B(91)-B(92)

  JVS(35) = B(87)

  JVS(36) = -B(93)

  JVS(37) = B(88)

  JVS(38) = -B(33)-B(34)

  JVS(39) = B(23)

  JVS(40) = -B(35)

  JVS(41) = B(24)

  JVS(42) = -B(141)-B(143)-B(144)

  JVS(43) = B(129)

  JVS(44) = B(130)

  JVS(45) = -B(142)

  JVS(46) = -B(145)-B(147)-B(148)

  JVS(47) = B(109)

  JVS(48) = B(110)

  JVS(49) = -B(146)

  JVS(50) = -B(151)

  JVS(51) = 0.54*B(217)

  JVS(52) = 0.22*B(204)

  JVS(53) = 0.28*B(85)

  JVS(54) = 0.45*B(97)

  JVS(55) = 0.28*B(86)+0.45*B(98)+0.22*B(205)+0.54*B(218)

  JVS(56) = -B(152)

  JVS(57) = -B(173)-B(174)

  JVS(58) = B(158)

  JVS(59) = -B(175)

  JVS(60) = B(159)

  JVS(61) = -B(184)-B(185)

  JVS(62) = B(180)

  JVS(63) = -B(186)

  JVS(64) = B(181)

  JVS(65) = -B(187)-B(188)

  JVS(66) = 0.04*B(178)

  JVS(67) = -B(189)

  JVS(68) = 0.04*B(179)

  JVS(69) = -B(215)-B(217)

  JVS(70) = -B(218)

  JVS(71) = -B(216)

  JVS(72) = -B(237)-B(239)

  JVS(73) = B(235)

  JVS(74) = -B(238)

  JVS(75) = B(236)

  JVS(76) = -B(2)-B(4)-B(21)

  JVS(77) = -B(22)

  JVS(78) = B(1)

  JVS(79) = -B(69)-B(70)-B(71)

  JVS(80) = B(51)

  JVS(81) = -B(72)

  JVS(82) = B(52)

  JVS(83) = -B(199)-B(201)

  JVS(84) = B(195)

  JVS(85) = -B(200)

  JVS(86) = B(196)

  JVS(87) = -B(21)

  JVS(88) = -B(19)-B(22)

  JVS(89) = 0.06*B(204)

  JVS(90) = 0.06*B(205)

  JVS(91) = -B(20)

  JVS(92) = 2*B(66)+2*B(67)

  JVS(93) = -B(53)-B(54)

  JVS(94) = B(40)

  JVS(95) = B(171)

  JVS(96) = B(41)+B(172)

  JVS(97) = B(49)

  JVS(98) = B(50)-B(55)

  JVS(99) = -B(212)-B(214)

  JVS(100) = B(210)

  JVS(101) = -B(213)

  JVS(102) = B(211)

  JVS(103) = -B(227)-B(229)

  JVS(104) = B(223)

  JVS(105) = -B(228)

  JVS(106) = B(224)

  JVS(107) = B(184)+0.7*B(185)

  JVS(108) = B(187)+B(188)

  JVS(109) = -B(190)-B(192)

  JVS(110) = 0.96*B(178)+B(182)

  JVS(111) = B(183)

  JVS(112) = 0.7*B(186)+B(189)-B(191)

  JVS(113) = 0.96*B(179)

  JVS(114) = 0

  JVS(115) = -B(103)-B(104)

  JVS(116) = B(101)

  JVS(117) = -B(105)

  JVS(118) = B(102)

  JVS(119) = 0.85*B(227)+0.67*B(229)

  JVS(120) = -B(230)-B(231)

  JVS(121) = 0.56*B(221)+0.88*B(225)

  JVS(122) = B(241)+0.67*B(243)

  JVS(123) = 0.88*B(226)

  JVS(124) = 0.85*B(228)-B(232)+B(242)

  JVS(125) = 0.56*B(222)

  JVS(126) = 0

  JVS(127) = B(153)+B(155)

  JVS(128) = 0.23*B(217)

  JVS(129) = 0.25*B(103)

  JVS(130) = -B(17)

  JVS(131) = 0

  JVS(132) = B(117)+B(119)

  JVS(133) = 0.31*B(204)

  JVS(134) = 0.14*B(85)

  JVS(135) = B(96)+0.22*B(97)

  JVS(136) = B(36)+B(37)+B(38)+B(40)

  JVS(137) = B(170)

  JVS(138) = 0.25*B(99)+0.5*B(112)+0.125*B(115)

  JVS(139) = 0.125*B(116)

  JVS(140) = B(41)

  JVS(141) = 0.14*B(86)+0.22*B(98)+0.31*B(205)+0.23*B(218)

  JVS(142) = -B(18)+B(39)+B(118)+B(154)

  JVS(143) = 0.25*B(100)

  JVS(144) = 0

  JVS(145) = -B(75)

  JVS(146) = 0.2*B(182)

  JVS(147) = 0.5*B(197)

  JVS(148) = 0.03*B(204)

  JVS(149) = 0.18*B(225)

  JVS(150) = 0.25*B(160)

  JVS(151) = 0.25*B(113)

  JVS(152) = 0.25*B(115)

  JVS(153) = B(30)+0.25*B(114)+0.25*B(116)+0.25*B(161)+0.2*B(183)+0.5*B(198)+0.18*B(226)

  JVS(154) = 0.03*B(205)

  JVS(155) = -B(76)

  JVS(156) = -B(106)-B(108)

  JVS(157) = B(83)

  JVS(158) = 0.12*B(89)

  JVS(159) = B(84)

  JVS(160) = -B(107)

  JVS(161) = 0.12*B(90)

  JVS(162) = 0.82*B(176)

  JVS(163) = 0.3*B(185)

  JVS(164) = -B(178)-B(180)-B(182)

  JVS(165) = -B(183)

  JVS(166) = 0.82*B(177)+0.3*B(186)

  JVS(167) = -B(179)

  JVS(168) = -B(181)

  JVS(169) = 0.7*B(199)

  JVS(170) = 0.5*B(103)

  JVS(171) = -B(117)-B(119)

  JVS(172) = 0.5*B(197)

  JVS(173) = 0.04*B(204)

  JVS(174) = B(120)

  JVS(175) = 0.9*B(97)

  JVS(176) = 0.25*B(113)

  JVS(177) = 0.5*B(99)+B(112)+0.5*B(115)

  JVS(178) = 0.25*B(114)+0.5*B(116)+0.5*B(198)

  JVS(179) = 0.9*B(98)+0.04*B(205)

  JVS(180) = -B(118)+B(121)+0.7*B(200)

  JVS(181) = 0.5*B(100)

  JVS(182) = 0

  JVS(183) = 0.3*B(199)

  JVS(184) = B(190)

  JVS(185) = 0

  JVS(186) = -B(193)-B(195)-B(197)

  JVS(187) = -B(198)

  JVS(188) = B(191)+0.3*B(200)

  JVS(189) = -B(194)

  JVS(190) = -B(196)

  JVS(191) = -B(202)-B(204)-B(206)

  JVS(192) = -B(207)

  JVS(193) = -B(205)

  JVS(194) = -B(203)

  JVS(195) = -B(81)-B(83)-B(85)

  JVS(196) = -B(84)

  JVS(197) = -B(86)

  JVS(198) = -B(82)

  JVS(199) = B(145)+B(148)

  JVS(200) = 0.5*B(212)

  JVS(201) = 0.25*B(103)

  JVS(202) = B(106)

  JVS(203) = 0.2*B(197)

  JVS(204) = 0

  JVS(205) = -B(120)-B(122)

  JVS(206) = 0

  JVS(207) = 0.25*B(113)

  JVS(208) = 0.25*B(99)+B(112)+0.375*B(115)

  JVS(209) = 0.25*B(114)+0.375*B(116)+0.2*B(198)

  JVS(210) = 0

  JVS(211) = 0

  JVS(212) = 0

  JVS(213) = B(107)-B(121)+B(146)+0.5*B(213)

  JVS(214) = 0.25*B(100)

  JVS(215) = 0

  JVS(216) = 0.666667*B(215)

  JVS(217) = 0.5*B(212)

  JVS(218) = B(202)

  JVS(219) = -B(208)-B(210)

  JVS(220) = 0

  JVS(221) = 0

  JVS(222) = B(203)+0.5*B(213)+0.666667*B(216)

  JVS(223) = -B(209)

  JVS(224) = -B(211)

  JVS(225) = B(91)+B(92)

  JVS(226) = B(108)

  JVS(227) = 0.65*B(85)

  JVS(228) = -B(94)-B(96)-B(97)

  JVS(229) = 0.88*B(89)+2*B(111)+0.5*B(113)

  JVS(230) = 0.5*B(114)

  JVS(231) = 0

  JVS(232) = 0.65*B(86)-B(98)

  JVS(233) = B(93)-B(95)

  JVS(234) = 0.88*B(90)

  JVS(235) = 0

  JVS(236) = 0.2*B(237)

  JVS(237) = B(231)

  JVS(238) = -B(233)-B(235)

  JVS(239) = 0

  JVS(240) = 0

  JVS(241) = 0

  JVS(242) = B(232)+0.2*B(238)

  JVS(243) = -B(234)

  JVS(244) = -B(236)

  JVS(245) = B(219)

  JVS(246) = 0.15*B(227)

  JVS(247) = -B(221)-B(223)-B(225)

  JVS(248) = -B(226)

  JVS(249) = B(220)+0.15*B(228)

  JVS(250) = -B(222)

  JVS(251) = -B(224)

  JVS(252) = B(206)

  JVS(253) = 0.02*B(208)

  JVS(254) = 0.015*B(233)

  JVS(255) = 0.16*B(221)

  JVS(256) = -B(241)-B(243)

  JVS(257) = 0

  JVS(258) = B(207)

  JVS(259) = 0

  JVS(260) = -B(242)

  JVS(261) = 0.02*B(209)+0.16*B(222)+0.015*B(234)

  JVS(262) = 0

  JVS(263) = 0.39*B(80)

  JVS(264) = B(153)+B(155)

  JVS(265) = B(91)

  JVS(266) = B(33)+0.3*B(34)

  JVS(267) = B(141)

  JVS(268) = B(217)

  JVS(269) = 0.25*B(21)

  JVS(270) = 0.25*B(22)

  JVS(271) = B(214)

  JVS(272) = 0.75*B(103)

  JVS(273) = B(75)

  JVS(274) = B(108)

  JVS(275) = 0.8*B(182)

  JVS(276) = B(193)+0.8*B(197)

  JVS(277) = 0.57*B(204)

  JVS(278) = 0.58*B(85)

  JVS(279) = B(122)

  JVS(280) = 0.98*B(208)

  JVS(281) = B(96)

  JVS(282) = 0.68*B(225)

  JVS(283) = -B(36)-B(37)-B(38)-B(40)

  JVS(284) = 0.75*B(160)

  JVS(285) = 0.88*B(89)+B(111)+1.25*B(113)

  JVS(286) = 0.75*B(99)+0.5*B(112)+1.125*B(115)

  JVS(287) = B(25)+B(27)+2*B(29)+B(30)+B(31)+1.25*B(114)+1.125*B(116)+B(131)+B(133)+0.75*B(161)+0.8*B(183)+0.8*B(198)&
               &+0.68*B(226)

  JVS(288) = B(132)+B(134)

  JVS(289) = B(32)-B(41)

  JVS(290) = 0.58*B(86)+0.57*B(205)+B(218)

  JVS(291) = 0

  JVS(292) = 0.3*B(35)-B(39)+B(76)+B(142)+B(154)

  JVS(293) = B(28)+0.88*B(90)+0.75*B(100)+B(194)+0.98*B(209)

  JVS(294) = B(26)

  JVS(295) = B(173)+0.7*B(174)

  JVS(296) = B(239)

  JVS(297) = B(214)

  JVS(298) = 0.33*B(229)

  JVS(299) = 0.47*B(204)

  JVS(300) = 0.98*B(208)

  JVS(301) = 0.985*B(233)

  JVS(302) = 0.84*B(221)+0.12*B(225)

  JVS(303) = 0.33*B(243)

  JVS(304) = -B(168)-B(170)-B(171)

  JVS(305) = 0.75*B(160)+B(162)+B(164)+B(166)

  JVS(306) = 0.75*B(161)+0.12*B(226)

  JVS(307) = B(163)

  JVS(308) = B(167)-B(172)

  JVS(309) = 0.47*B(205)

  JVS(310) = -B(169)+0.7*B(175)

  JVS(311) = B(165)+0.98*B(209)+0.84*B(222)+0.985*B(234)

  JVS(312) = 0

  JVS(313) = B(156)

  JVS(314) = 0.18*B(176)

  JVS(315) = 0.3*B(174)

  JVS(316) = 0.33*B(229)

  JVS(317) = B(230)

  JVS(318) = 0.06*B(204)

  JVS(319) = 0.28*B(221)+0.12*B(225)

  JVS(320) = 0.33*B(243)

  JVS(321) = -B(158)-B(160)-B(162)-B(164)-B(166)

  JVS(322) = -B(161)+0.12*B(226)

  JVS(323) = -B(163)

  JVS(324) = -B(167)

  JVS(325) = 0.06*B(205)

  JVS(326) = B(157)+0.3*B(175)+0.18*B(177)

  JVS(327) = -B(165)+0.28*B(222)

  JVS(328) = -B(159)

  JVS(329) = B(81)

  JVS(330) = -B(87)-B(89)-2*B(111)-B(113)

  JVS(331) = -B(114)

  JVS(332) = 0

  JVS(333) = 0

  JVS(334) = B(82)

  JVS(335) = -B(90)

  JVS(336) = -B(88)

  JVS(337) = B(147)

  JVS(338) = B(104)

  JVS(339) = 0.1*B(85)

  JVS(340) = B(94)

  JVS(341) = 0

  JVS(342) = -B(99)-B(101)-B(109)-2*B(112)-B(115)

  JVS(343) = -B(116)

  JVS(344) = 0

  JVS(345) = 0.1*B(86)

  JVS(346) = -B(110)

  JVS(347) = B(95)+B(105)

  JVS(348) = -B(100)

  JVS(349) = -B(102)

  JVS(350) = B(79)+0.61*B(80)

  JVS(351) = B(138)

  JVS(352) = B(149)

  JVS(353) = 0.7*B(34)

  JVS(354) = 0.75*B(21)

  JVS(355) = B(19)+0.75*B(22)

  JVS(356) = B(192)

  JVS(357) = -B(182)

  JVS(358) = -B(197)

  JVS(359) = 0.07*B(204)

  JVS(360) = 0.08*B(85)

  JVS(361) = -B(225)

  JVS(362) = B(170)

  JVS(363) = -B(160)+0.82*B(162)

  JVS(364) = -B(113)

  JVS(365) = -B(115)

  JVS(366) = -B(23)-B(25)-B(27)-2*B(29)-2*B(30)-B(31)-B(77)-B(114)-B(116)-B(133)-B(161)-B(183)-B(198)-B(226)

  JVS(367) = B(127)-B(134)+2*B(135)+B(136)+0.82*B(163)

  JVS(368) = -B(32)+B(137)

  JVS(369) = 0.08*B(86)+0.07*B(205)

  JVS(370) = -B(78)

  JVS(371) = B(20)+0.7*B(35)+B(150)

  JVS(372) = -B(28)+B(128)

  JVS(373) = -B(24)-B(26)

  JVS(374) = 2*B(240)

  JVS(375) = B(139)

  JVS(376) = B(143)+B(144)

  JVS(377) = B(239)

  JVS(378) = B(201)

  JVS(379) = B(192)

  JVS(380) = 0.25*B(103)

  JVS(381) = B(230)

  JVS(382) = 0

  JVS(383) = B(117)+B(119)

  JVS(384) = B(193)+0.3*B(197)

  JVS(385) = 0.23*B(204)

  JVS(386) = 0.1*B(85)

  JVS(387) = B(122)

  JVS(388) = 0

  JVS(389) = B(96)+0.1*B(97)

  JVS(390) = 0.985*B(233)

  JVS(391) = 0

  JVS(392) = 0

  JVS(393) = B(168)+B(171)

  JVS(394) = -B(162)

  JVS(395) = 0

  JVS(396) = 0.25*B(99)+0.125*B(115)

  JVS(397) = 0.125*B(116)-B(131)-B(133)+0.3*B(198)

  JVS(398) = -B(123)-B(125)-B(127)-B(129)-B(132)-B(134)-2*B(135)-B(136)-B(163)

  JVS(399) = -B(137)+B(172)

  JVS(400) = 0.1*B(86)+0.1*B(98)+0.23*B(205)

  JVS(401) = -B(130)

  JVS(402) = B(118)+B(140)+B(169)

  JVS(403) = 0.25*B(100)-B(128)+B(194)+0.985*B(234)

  JVS(404) = -B(124)-B(126)

  JVS(405) = B(64)+B(65)

  JVS(406) = 0.39*B(80)

  JVS(407) = 0.39*B(69)

  JVS(408) = B(54)

  JVS(409) = -B(206)

  JVS(410) = -B(83)

  JVS(411) = -B(40)

  JVS(412) = -B(171)

  JVS(413) = -B(166)

  JVS(414) = 0

  JVS(415) = 0

  JVS(416) = -B(31)

  JVS(417) = -B(136)

  JVS(418) = -B(32)-B(41)-B(56)-B(57)-B(58)-B(60)-B(62)-B(84)-B(137)-B(167)-B(172)-B(207)

  JVS(419) = B(47)

  JVS(420) = B(48)-B(61)

  JVS(421) = B(55)

  JVS(422) = -B(59)

  JVS(423) = -B(63)

  JVS(424) = -B(217)

  JVS(425) = B(2)

  JVS(426) = 0

  JVS(427) = -B(204)

  JVS(428) = -B(85)

  JVS(429) = -B(97)

  JVS(430) = 0

  JVS(431) = 0

  JVS(432) = B(125)

  JVS(433) = B(56)

  JVS(434) = -B(1)-B(7)-B(9)-B(42)-B(47)-B(86)-B(98)-B(205)-B(218)

  JVS(435) = B(46)-B(48)

  JVS(436) = -B(8)

  JVS(437) = -B(43)

  JVS(438) = -B(10)+B(126)

  JVS(439) = B(64)+B(65)

  JVS(440) = B(79)+0.61*B(80)

  JVS(441) = B(153)+B(155)

  JVS(442) = B(141)+B(143)+B(144)

  JVS(443) = B(145)+B(147)+B(148)

  JVS(444) = B(187)+B(188)

  JVS(445) = 0.61*B(69)+B(70)+B(71)

  JVS(446) = B(53)

  JVS(447) = B(108)

  JVS(448) = 0.96*B(178)

  JVS(449) = B(193)

  JVS(450) = 0

  JVS(451) = 0.98*B(208)

  JVS(452) = 0.985*B(233)

  JVS(453) = 0.84*B(221)

  JVS(454) = B(241)+B(243)

  JVS(455) = 0

  JVS(456) = 0

  JVS(457) = B(164)+B(166)

  JVS(458) = 0.88*B(89)

  JVS(459) = B(99)-B(109)

  JVS(460) = B(27)+B(31)-B(77)

  JVS(461) = B(127)-B(129)+B(136)

  JVS(462) = B(32)+B(56)+2*B(58)-B(60)+B(62)+B(137)+B(167)

  JVS(463) = B(42)-B(47)

  JVS(464) = -B(46)-B(48)-B(49)-B(51)-B(61)-B(78)-B(110)-B(130)

  JVS(465) = -B(50)+B(72)+B(142)+B(146)+B(154)+B(189)+B(242)

  JVS(466) = B(28)+B(43)+B(44)+2*B(59)+0.88*B(90)+B(100)+B(128)+B(165)+0.96*B(179)+B(194)+0.98*B(209)+0.84*B(222)+0.985&
               &*B(234)

  JVS(467) = B(45)-B(52)+B(63)

  JVS(468) = -B(156)

  JVS(469) = -B(219)

  JVS(470) = -B(176)

  JVS(471) = -B(153)

  JVS(472) = 2*B(14)-B(15)

  JVS(473) = B(138)-B(139)

  JVS(474) = -B(149)

  JVS(475) = B(91)

  JVS(476) = B(33)-0.7*B(34)

  JVS(477) = -B(141)

  JVS(478) = -B(145)

  JVS(479) = -B(151)

  JVS(480) = B(173)-0.3*B(174)

  JVS(481) = B(184)-0.3*B(185)

  JVS(482) = -B(188)

  JVS(483) = -B(215)+0.12*B(217)

  JVS(484) = -0.2*B(237)+B(239)

  JVS(485) = 2*B(4)+0.75*B(21)

  JVS(486) = 0.39*B(69)-B(71)

  JVS(487) = -0.3*B(199)+B(201)

  JVS(488) = -B(19)+0.75*B(22)

  JVS(489) = B(53)-B(54)

  JVS(490) = -0.5*B(212)+B(214)

  JVS(491) = -0.15*B(227)+B(229)

  JVS(492) = -B(190)

  JVS(493) = B(103)-B(104)

  JVS(494) = -B(231)

  JVS(495) = -B(17)

  JVS(496) = -B(75)

  JVS(497) = -B(106)

  JVS(498) = 0

  JVS(499) = -B(117)

  JVS(500) = 0

  JVS(501) = -B(202)+0.33*B(204)

  JVS(502) = -B(81)+0.25*B(85)

  JVS(503) = -B(120)

  JVS(504) = 0

  JVS(505) = -B(94)+0.19*B(97)

  JVS(506) = 0

  JVS(507) = 0

  JVS(508) = -B(241)

  JVS(509) = -B(38)

  JVS(510) = -B(168)

  JVS(511) = 0

  JVS(512) = 0

  JVS(513) = 0

  JVS(514) = 0

  JVS(515) = 0

  JVS(516) = B(62)

  JVS(517) = -B(7)+B(9)+0.25*B(86)+0.19*B(98)+0.33*B(205)+0.12*B(218)

  JVS(518) = -B(49)

  JVS(519) = -B(8)-B(11)-B(16)-B(18)-B(20)-0.7*B(35)-B(39)-B(50)-B(55)-B(72)-B(73)-B(76)-B(82)-B(95)-B(105)-B(107)&
               &-B(118)-B(121)-B(140)-B(142)-B(146)-B(150)-B(152)-B(154)-B(157)-B(169)-0.3*B(175)-B(177)-0.3*B(186)-B(189)&
               &-B(191)-0.3*B(200)-B(203)-0.5*B(213)-B(216)-B(220)-0.15*B(228)-B(232)-0.2*B(238)-B(242)

  JVS(520) = B(44)

  JVS(521) = B(10)-B(12)+B(45)+B(63)

  JVS(522) = -B(178)

  JVS(523) = -B(193)

  JVS(524) = -B(208)

  JVS(525) = -B(233)

  JVS(526) = -B(221)

  JVS(527) = 0

  JVS(528) = -B(164)

  JVS(529) = -B(89)

  JVS(530) = -B(99)

  JVS(531) = -B(27)

  JVS(532) = -B(127)

  JVS(533) = B(57)-B(58)

  JVS(534) = -B(42)

  JVS(535) = B(46)

  JVS(536) = 0

  JVS(537) = -B(28)-B(43)-B(44)-B(59)-B(90)-B(100)-B(128)-B(165)-B(179)-B(194)-B(209)-B(222)-B(234)

  JVS(538) = -B(45)

  JVS(539) = B(15)

  JVS(540) = B(91)

  JVS(541) = B(33)

  JVS(542) = B(151)

  JVS(543) = B(173)

  JVS(544) = B(184)

  JVS(545) = B(187)

  JVS(546) = 0.22*B(217)

  JVS(547) = 0.4*B(21)

  JVS(548) = 0.61*B(69)+B(70)

  JVS(549) = B(201)

  JVS(550) = 0.4*B(22)

  JVS(551) = B(214)

  JVS(552) = 0.67*B(229)

  JVS(553) = 0.75*B(103)

  JVS(554) = B(17)

  JVS(555) = B(75)

  JVS(556) = B(108)

  JVS(557) = 0.96*B(178)-B(180)+0.8*B(182)

  JVS(558) = B(119)

  JVS(559) = -B(195)+0.3*B(197)

  JVS(560) = 0.26*B(204)

  JVS(561) = 0.25*B(85)

  JVS(562) = B(120)+B(122)

  JVS(563) = 0.98*B(208)-B(210)

  JVS(564) = B(96)+0.32*B(97)

  JVS(565) = -B(235)

  JVS(566) = 0.56*B(221)-B(223)+1.23*B(225)

  JVS(567) = 0.67*B(243)

  JVS(568) = 2*B(36)+B(38)+B(40)

  JVS(569) = B(170)

  JVS(570) = -B(158)+B(160)+0.82*B(162)+B(164)+B(166)

  JVS(571) = -B(87)+0.88*B(89)+B(111)+B(113)

  JVS(572) = 0.75*B(99)-B(101)+B(112)+0.875*B(115)

  JVS(573) = -B(23)-B(25)+B(27)+2*B(29)+B(31)+B(114)+0.875*B(116)+B(131)+B(161)+0.8*B(183)+0.3*B(198)+1.23*B(226)

  JVS(574) = -B(123)-B(125)+B(132)+0.82*B(163)

  JVS(575) = B(32)+B(41)-B(62)+B(167)

  JVS(576) = B(7)-B(9)+0.25*B(86)+0.32*B(98)+0.26*B(205)+0.22*B(218)

  JVS(577) = -B(51)

  JVS(578) = B(8)-B(11)+B(16)+B(18)+B(39)+B(73)+B(76)+B(121)+B(152)

  JVS(579) = B(28)-B(44)+0.88*B(90)+0.75*B(100)+B(165)+0.96*B(179)+0.98*B(209)+0.56*B(222)

  JVS(580) = -B(10)-B(12)-2*B(13)-B(24)-B(26)-B(45)-B(52)-B(63)-B(88)-B(102)-B(124)-B(126)-B(159)-B(181)-B(196)-B(211)&
               &-B(224)-B(236)
      
END SUBROUTINE nmhc9_Jac_SP














SUBROUTINE nmhc9_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(580), W(57), a
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
      
END SUBROUTINE nmhc9_KppDecomp



SUBROUTINE nmhc9_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(580), W(57), a
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
      
END SUBROUTINE nmhc9_KppDecompCmplx


SUBROUTINE nmhc9_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(580), X(57), sum

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
      
END SUBROUTINE nmhc9_KppSolveIndirect


SUBROUTINE nmhc9_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(580), X(57), sum

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
      
END SUBROUTINE nmhc9_KppSolveCmplx













SUBROUTINE nmhc9_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(24) = X(24)-JVS(87)*X(21)
  X(25) = X(25)-JVS(92)*X(5)
  X(28) = X(28)-JVS(107)*X(17)-JVS(108)*X(18)
  X(30) = X(30)-JVS(119)*X(27)
  X(31) = X(31)-JVS(127)*X(7)-JVS(128)*X(19)-JVS(129)*X(29)
  X(34) = X(34)-JVS(162)*X(4)-JVS(163)*X(17)
  X(35) = X(35)-JVS(169)*X(23)-JVS(170)*X(29)
  X(36) = X(36)-JVS(183)*X(23)-JVS(184)*X(28)-JVS(185)*X(34)
  X(39) = X(39)-JVS(199)*X(14)-JVS(200)*X(26)-JVS(201)*X(29)-JVS(202)*X(33)-JVS(203)*X(36)-JVS(204)*X(38)
  X(40) = X(40)-JVS(216)*X(19)-JVS(217)*X(26)-JVS(218)*X(37)
  X(41) = X(41)-JVS(225)*X(11)-JVS(226)*X(33)-JVS(227)*X(38)
  X(42) = X(42)-JVS(236)*X(20)-JVS(237)*X(30)
  X(43) = X(43)-JVS(245)*X(2)-JVS(246)*X(27)
  X(44) = X(44)-JVS(252)*X(37)-JVS(253)*X(40)-JVS(254)*X(42)-JVS(255)*X(43)
  X(45) = X(45)-JVS(263)*X(6)-JVS(264)*X(7)-JVS(265)*X(11)-JVS(266)*X(12)-JVS(267)*X(13)-JVS(268)*X(19)-JVS(269)*X(21)&
            &-JVS(270)*X(24)-JVS(271)*X(26)-JVS(272)*X(29)-JVS(273)*X(32)-JVS(274)*X(33)-JVS(275)*X(34)-JVS(276)*X(36)&
            &-JVS(277)*X(37)-JVS(278)*X(38)-JVS(279)*X(39)-JVS(280)*X(40)-JVS(281)*X(41)-JVS(282)*X(43)
  X(46) = X(46)-JVS(295)*X(16)-JVS(296)*X(20)-JVS(297)*X(26)-JVS(298)*X(27)-JVS(299)*X(37)-JVS(300)*X(40)-JVS(301)*X(42)&
            &-JVS(302)*X(43)-JVS(303)*X(44)
  X(47) = X(47)-JVS(313)*X(1)-JVS(314)*X(4)-JVS(315)*X(16)-JVS(316)*X(27)-JVS(317)*X(30)-JVS(318)*X(37)-JVS(319)*X(43)&
            &-JVS(320)*X(44)
  X(48) = X(48)-JVS(329)*X(38)
  X(49) = X(49)-JVS(337)*X(14)-JVS(338)*X(29)-JVS(339)*X(38)-JVS(340)*X(41)-JVS(341)*X(48)
  X(50) = X(50)-JVS(350)*X(6)-JVS(351)*X(9)-JVS(352)*X(10)-JVS(353)*X(12)-JVS(354)*X(21)-JVS(355)*X(24)-JVS(356)*X(28)&
            &-JVS(357)*X(34)-JVS(358)*X(36)-JVS(359)*X(37)-JVS(360)*X(38)-JVS(361)*X(43)-JVS(362)*X(46)-JVS(363)*X(47)&
            &-JVS(364)*X(48)-JVS(365)*X(49)
  X(51) = X(51)-JVS(374)*X(3)-JVS(375)*X(9)-JVS(376)*X(13)-JVS(377)*X(20)-JVS(378)*X(23)-JVS(379)*X(28)-JVS(380)*X(29)&
            &-JVS(381)*X(30)-JVS(382)*X(34)-JVS(383)*X(35)-JVS(384)*X(36)-JVS(385)*X(37)-JVS(386)*X(38)-JVS(387)*X(39)&
            &-JVS(388)*X(40)-JVS(389)*X(41)-JVS(390)*X(42)-JVS(391)*X(43)-JVS(392)*X(44)-JVS(393)*X(46)-JVS(394)*X(47)&
            &-JVS(395)*X(48)-JVS(396)*X(49)-JVS(397)*X(50)
  X(52) = X(52)-JVS(405)*X(5)-JVS(406)*X(6)-JVS(407)*X(22)-JVS(408)*X(25)-JVS(409)*X(37)-JVS(410)*X(38)-JVS(411)*X(45)&
            &-JVS(412)*X(46)-JVS(413)*X(47)-JVS(414)*X(48)-JVS(415)*X(49)-JVS(416)*X(50)-JVS(417)*X(51)
  X(53) = X(53)-JVS(424)*X(19)-JVS(425)*X(21)-JVS(426)*X(24)-JVS(427)*X(37)-JVS(428)*X(38)-JVS(429)*X(41)-JVS(430)*X(48)&
            &-JVS(431)*X(50)-JVS(432)*X(51)-JVS(433)*X(52)
  X(54) = X(54)-JVS(439)*X(5)-JVS(440)*X(6)-JVS(441)*X(7)-JVS(442)*X(13)-JVS(443)*X(14)-JVS(444)*X(18)-JVS(445)*X(22)&
            &-JVS(446)*X(25)-JVS(447)*X(33)-JVS(448)*X(34)-JVS(449)*X(36)-JVS(450)*X(38)-JVS(451)*X(40)-JVS(452)*X(42)&
            &-JVS(453)*X(43)-JVS(454)*X(44)-JVS(455)*X(45)-JVS(456)*X(46)-JVS(457)*X(47)-JVS(458)*X(48)-JVS(459)*X(49)&
            &-JVS(460)*X(50)-JVS(461)*X(51)-JVS(462)*X(52)-JVS(463)*X(53)
  X(55) = X(55)-JVS(468)*X(1)-JVS(469)*X(2)-JVS(470)*X(4)-JVS(471)*X(7)-JVS(472)*X(8)-JVS(473)*X(9)-JVS(474)*X(10)&
            &-JVS(475)*X(11)-JVS(476)*X(12)-JVS(477)*X(13)-JVS(478)*X(14)-JVS(479)*X(15)-JVS(480)*X(16)-JVS(481)*X(17)&
            &-JVS(482)*X(18)-JVS(483)*X(19)-JVS(484)*X(20)-JVS(485)*X(21)-JVS(486)*X(22)-JVS(487)*X(23)-JVS(488)*X(24)&
            &-JVS(489)*X(25)-JVS(490)*X(26)-JVS(491)*X(27)-JVS(492)*X(28)-JVS(493)*X(29)-JVS(494)*X(30)-JVS(495)*X(31)&
            &-JVS(496)*X(32)-JVS(497)*X(33)-JVS(498)*X(34)-JVS(499)*X(35)-JVS(500)*X(36)-JVS(501)*X(37)-JVS(502)*X(38)&
            &-JVS(503)*X(39)-JVS(504)*X(40)-JVS(505)*X(41)-JVS(506)*X(42)-JVS(507)*X(43)-JVS(508)*X(44)-JVS(509)*X(45)&
            &-JVS(510)*X(46)-JVS(511)*X(47)-JVS(512)*X(48)-JVS(513)*X(49)-JVS(514)*X(50)-JVS(515)*X(51)-JVS(516)*X(52)&
            &-JVS(517)*X(53)-JVS(518)*X(54)
  X(56) = X(56)-JVS(522)*X(34)-JVS(523)*X(36)-JVS(524)*X(40)-JVS(525)*X(42)-JVS(526)*X(43)-JVS(527)*X(44)-JVS(528)*X(47)&
            &-JVS(529)*X(48)-JVS(530)*X(49)-JVS(531)*X(50)-JVS(532)*X(51)-JVS(533)*X(52)-JVS(534)*X(53)-JVS(535)*X(54)&
            &-JVS(536)*X(55)
  X(57) = X(57)-JVS(539)*X(8)-JVS(540)*X(11)-JVS(541)*X(12)-JVS(542)*X(15)-JVS(543)*X(16)-JVS(544)*X(17)-JVS(545)*X(18)&
            &-JVS(546)*X(19)-JVS(547)*X(21)-JVS(548)*X(22)-JVS(549)*X(23)-JVS(550)*X(24)-JVS(551)*X(26)-JVS(552)*X(27)&
            &-JVS(553)*X(29)-JVS(554)*X(31)-JVS(555)*X(32)-JVS(556)*X(33)-JVS(557)*X(34)-JVS(558)*X(35)-JVS(559)*X(36)&
            &-JVS(560)*X(37)-JVS(561)*X(38)-JVS(562)*X(39)-JVS(563)*X(40)-JVS(564)*X(41)-JVS(565)*X(42)-JVS(566)*X(43)&
            &-JVS(567)*X(44)-JVS(568)*X(45)-JVS(569)*X(46)-JVS(570)*X(47)-JVS(571)*X(48)-JVS(572)*X(49)-JVS(573)*X(50)&
            &-JVS(574)*X(51)-JVS(575)*X(52)-JVS(576)*X(53)-JVS(577)*X(54)-JVS(578)*X(55)-JVS(579)*X(56)
  X(57) = X(57)/JVS(580)
  X(56) = (X(56)-JVS(538)*X(57))/(JVS(537))
  X(55) = (X(55)-JVS(520)*X(56)-JVS(521)*X(57))/(JVS(519))
  X(54) = (X(54)-JVS(465)*X(55)-JVS(466)*X(56)-JVS(467)*X(57))/(JVS(464))
  X(53) = (X(53)-JVS(435)*X(54)-JVS(436)*X(55)-JVS(437)*X(56)-JVS(438)*X(57))/(JVS(434))
  X(52) = (X(52)-JVS(419)*X(53)-JVS(420)*X(54)-JVS(421)*X(55)-JVS(422)*X(56)-JVS(423)*X(57))/(JVS(418))
  X(51) = (X(51)-JVS(399)*X(52)-JVS(400)*X(53)-JVS(401)*X(54)-JVS(402)*X(55)-JVS(403)*X(56)-JVS(404)*X(57))/(JVS(398))
  X(50) = (X(50)-JVS(367)*X(51)-JVS(368)*X(52)-JVS(369)*X(53)-JVS(370)*X(54)-JVS(371)*X(55)-JVS(372)*X(56)-JVS(373)&
            &*X(57))/(JVS(366))
  X(49) = (X(49)-JVS(343)*X(50)-JVS(344)*X(52)-JVS(345)*X(53)-JVS(346)*X(54)-JVS(347)*X(55)-JVS(348)*X(56)-JVS(349)&
            &*X(57))/(JVS(342))
  X(48) = (X(48)-JVS(331)*X(50)-JVS(332)*X(52)-JVS(333)*X(53)-JVS(334)*X(55)-JVS(335)*X(56)-JVS(336)*X(57))/(JVS(330))
  X(47) = (X(47)-JVS(322)*X(50)-JVS(323)*X(51)-JVS(324)*X(52)-JVS(325)*X(53)-JVS(326)*X(55)-JVS(327)*X(56)-JVS(328)&
            &*X(57))/(JVS(321))
  X(46) = (X(46)-JVS(305)*X(47)-JVS(306)*X(50)-JVS(307)*X(51)-JVS(308)*X(52)-JVS(309)*X(53)-JVS(310)*X(55)-JVS(311)&
            &*X(56)-JVS(312)*X(57))/(JVS(304))
  X(45) = (X(45)-JVS(284)*X(47)-JVS(285)*X(48)-JVS(286)*X(49)-JVS(287)*X(50)-JVS(288)*X(51)-JVS(289)*X(52)-JVS(290)&
            &*X(53)-JVS(291)*X(54)-JVS(292)*X(55)-JVS(293)*X(56)-JVS(294)*X(57))/(JVS(283))
  X(44) = (X(44)-JVS(257)*X(50)-JVS(258)*X(52)-JVS(259)*X(53)-JVS(260)*X(55)-JVS(261)*X(56)-JVS(262)*X(57))/(JVS(256))
  X(43) = (X(43)-JVS(248)*X(50)-JVS(249)*X(55)-JVS(250)*X(56)-JVS(251)*X(57))/(JVS(247))
  X(42) = (X(42)-JVS(239)*X(43)-JVS(240)*X(44)-JVS(241)*X(50)-JVS(242)*X(55)-JVS(243)*X(56)-JVS(244)*X(57))/(JVS(238))
  X(41) = (X(41)-JVS(229)*X(48)-JVS(230)*X(50)-JVS(231)*X(52)-JVS(232)*X(53)-JVS(233)*X(55)-JVS(234)*X(56)-JVS(235)&
            &*X(57))/(JVS(228))
  X(40) = (X(40)-JVS(220)*X(52)-JVS(221)*X(53)-JVS(222)*X(55)-JVS(223)*X(56)-JVS(224)*X(57))/(JVS(219))
  X(39) = (X(39)-JVS(206)*X(40)-JVS(207)*X(48)-JVS(208)*X(49)-JVS(209)*X(50)-JVS(210)*X(52)-JVS(211)*X(53)-JVS(212)&
            &*X(54)-JVS(213)*X(55)-JVS(214)*X(56)-JVS(215)*X(57))/(JVS(205))
  X(38) = (X(38)-JVS(196)*X(52)-JVS(197)*X(53)-JVS(198)*X(55))/(JVS(195))
  X(37) = (X(37)-JVS(192)*X(52)-JVS(193)*X(53)-JVS(194)*X(55))/(JVS(191))
  X(36) = (X(36)-JVS(187)*X(50)-JVS(188)*X(55)-JVS(189)*X(56)-JVS(190)*X(57))/(JVS(186))
  X(35) = (X(35)-JVS(172)*X(36)-JVS(173)*X(37)-JVS(174)*X(39)-JVS(175)*X(41)-JVS(176)*X(48)-JVS(177)*X(49)-JVS(178)&
            &*X(50)-JVS(179)*X(53)-JVS(180)*X(55)-JVS(181)*X(56)-JVS(182)*X(57))/(JVS(171))
  X(34) = (X(34)-JVS(165)*X(50)-JVS(166)*X(55)-JVS(167)*X(56)-JVS(168)*X(57))/(JVS(164))
  X(33) = (X(33)-JVS(157)*X(38)-JVS(158)*X(48)-JVS(159)*X(52)-JVS(160)*X(55)-JVS(161)*X(56))/(JVS(156))
  X(32) = (X(32)-JVS(146)*X(34)-JVS(147)*X(36)-JVS(148)*X(37)-JVS(149)*X(43)-JVS(150)*X(47)-JVS(151)*X(48)-JVS(152)&
            &*X(49)-JVS(153)*X(50)-JVS(154)*X(53)-JVS(155)*X(55))/(JVS(145))
  X(31) = (X(31)-JVS(131)*X(33)-JVS(132)*X(35)-JVS(133)*X(37)-JVS(134)*X(38)-JVS(135)*X(41)-JVS(136)*X(45)-JVS(137)&
            &*X(46)-JVS(138)*X(49)-JVS(139)*X(50)-JVS(140)*X(52)-JVS(141)*X(53)-JVS(142)*X(55)-JVS(143)*X(56)-JVS(144)&
            &*X(57))/(JVS(130))
  X(30) = (X(30)-JVS(121)*X(43)-JVS(122)*X(44)-JVS(123)*X(50)-JVS(124)*X(55)-JVS(125)*X(56)-JVS(126)*X(57))/(JVS(120))
  X(29) = (X(29)-JVS(116)*X(49)-JVS(117)*X(55)-JVS(118)*X(57))/(JVS(115))
  X(28) = (X(28)-JVS(110)*X(34)-JVS(111)*X(50)-JVS(112)*X(55)-JVS(113)*X(56)-JVS(114)*X(57))/(JVS(109))
  X(27) = (X(27)-JVS(104)*X(43)-JVS(105)*X(55)-JVS(106)*X(57))/(JVS(103))
  X(26) = (X(26)-JVS(100)*X(40)-JVS(101)*X(55)-JVS(102)*X(57))/(JVS(99))
  X(25) = (X(25)-JVS(94)*X(45)-JVS(95)*X(46)-JVS(96)*X(52)-JVS(97)*X(54)-JVS(98)*X(55))/(JVS(93))
  X(24) = (X(24)-JVS(89)*X(37)-JVS(90)*X(53)-JVS(91)*X(55))/(JVS(88))
  X(23) = (X(23)-JVS(84)*X(36)-JVS(85)*X(55)-JVS(86)*X(57))/(JVS(83))
  X(22) = (X(22)-JVS(80)*X(54)-JVS(81)*X(55)-JVS(82)*X(57))/(JVS(79))
  X(21) = (X(21)-JVS(77)*X(24)-JVS(78)*X(53))/(JVS(76))
  X(20) = (X(20)-JVS(73)*X(42)-JVS(74)*X(55)-JVS(75)*X(57))/(JVS(72))
  X(19) = (X(19)-JVS(70)*X(53)-JVS(71)*X(55))/(JVS(69))
  X(18) = (X(18)-JVS(66)*X(34)-JVS(67)*X(55)-JVS(68)*X(56))/(JVS(65))
  X(17) = (X(17)-JVS(62)*X(34)-JVS(63)*X(55)-JVS(64)*X(57))/(JVS(61))
  X(16) = (X(16)-JVS(58)*X(47)-JVS(59)*X(55)-JVS(60)*X(57))/(JVS(57))
  X(15) = (X(15)-JVS(51)*X(19)-JVS(52)*X(37)-JVS(53)*X(38)-JVS(54)*X(41)-JVS(55)*X(53)-JVS(56)*X(55))/(JVS(50))
  X(14) = (X(14)-JVS(47)*X(49)-JVS(48)*X(54)-JVS(49)*X(55))/(JVS(46))
  X(13) = (X(13)-JVS(43)*X(51)-JVS(44)*X(54)-JVS(45)*X(55))/(JVS(42))
  X(12) = (X(12)-JVS(39)*X(50)-JVS(40)*X(55)-JVS(41)*X(57))/(JVS(38))
  X(11) = (X(11)-JVS(35)*X(48)-JVS(36)*X(55)-JVS(37)*X(57))/(JVS(34))
  X(10) = (X(10)-JVS(29)*X(47)-JVS(30)*X(50)-JVS(31)*X(51)-JVS(32)*X(55)-JVS(33)*X(57))/(JVS(28))
  X(9) = (X(9)-JVS(25)*X(51)-JVS(26)*X(55)-JVS(27)*X(57))/(JVS(24))
  X(8) = (X(8)-JVS(20)*X(38)-JVS(21)*X(53)-JVS(22)*X(55)-JVS(23)*X(57))/(JVS(19))
  X(7) = (X(7)-JVS(17)*X(33)-JVS(18)*X(55))/(JVS(16))
  X(6) = (X(6)-JVS(14)*X(50)-JVS(15)*X(54))/(JVS(13))
  X(5) = (X(5)-JVS(11)*X(52)-JVS(12)*X(54))/(JVS(10))
  X(4) = (X(4)-JVS(9)*X(55))/(JVS(8))
  X(3) = (X(3)-JVS(6)*X(20)-JVS(7)*X(55))/(JVS(5))
  X(2) = (X(2)-JVS(4)*X(55))/(JVS(3))
  X(1) = (X(1)-JVS(2)*X(55))/(JVS(1))
      
END SUBROUTINE nmhc9_KppSolve
























      SUBROUTINE nmhc9_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE nmhc9_WCOPY



      SUBROUTINE nmhc9_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE nmhc9_WAXPY




      SUBROUTINE nmhc9_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE nmhc9_WSCAL


      REAL(kind=dp) FUNCTION nmhc9_WLAMCH( C )








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
          CALL nmhc9_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      nmhc9_WLAMCH = Eps

      END FUNCTION nmhc9_WLAMCH
     
      SUBROUTINE nmhc9_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE nmhc9_WLAMCH_ADD




      SUBROUTINE nmhc9_SET2ZERO(N,Y)




      
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

      END SUBROUTINE nmhc9_SET2ZERO



      REAL(kind=dp) FUNCTION nmhc9_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      nmhc9_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        nmhc9_WDOT = nmhc9_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         nmhc9_WDOT = nmhc9_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          nmhc9_WDOT = nmhc9_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        nmhc9_WDOT = nmhc9_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION nmhc9_WDOT                                          




   SUBROUTINE decomp_nmhc9( JVS, IER )
   
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
   W( 55 ) = JVS( 2 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 55 )
  IF ( ABS(  JVS( 3 )) < TINY(a) ) THEN
         IER = 2                                       
         RETURN
  END IF
   W( 2 ) = JVS( 3 )
   W( 55 ) = JVS( 4 )
  JVS( 3) = W( 2 )
  JVS( 4) = W( 55 )
  IF ( ABS(  JVS( 5 )) < TINY(a) ) THEN
         IER = 3                                       
         RETURN
  END IF
   W( 3 ) = JVS( 5 )
   W( 20 ) = JVS( 6 )
   W( 55 ) = JVS( 7 )
  JVS( 5) = W( 3 )
  JVS( 6) = W( 20 )
  JVS( 7) = W( 55 )
  IF ( ABS(  JVS( 8 )) < TINY(a) ) THEN
         IER = 4                                       
         RETURN
  END IF
   W( 4 ) = JVS( 8 )
   W( 55 ) = JVS( 9 )
  JVS( 8) = W( 4 )
  JVS( 9) = W( 55 )
  IF ( ABS(  JVS( 10 )) < TINY(a) ) THEN
         IER = 5                                       
         RETURN
  END IF
   W( 5 ) = JVS( 10 )
   W( 52 ) = JVS( 11 )
   W( 54 ) = JVS( 12 )
  JVS( 10) = W( 5 )
  JVS( 11) = W( 52 )
  JVS( 12) = W( 54 )
  IF ( ABS(  JVS( 13 )) < TINY(a) ) THEN
         IER = 6                                       
         RETURN
  END IF
   W( 6 ) = JVS( 13 )
   W( 50 ) = JVS( 14 )
   W( 54 ) = JVS( 15 )
  JVS( 13) = W( 6 )
  JVS( 14) = W( 50 )
  JVS( 15) = W( 54 )
  IF ( ABS(  JVS( 16 )) < TINY(a) ) THEN
         IER = 7                                       
         RETURN
  END IF
   W( 7 ) = JVS( 16 )
   W( 33 ) = JVS( 17 )
   W( 55 ) = JVS( 18 )
  JVS( 16) = W( 7 )
  JVS( 17) = W( 33 )
  JVS( 18) = W( 55 )
  IF ( ABS(  JVS( 19 )) < TINY(a) ) THEN
         IER = 8                                       
         RETURN
  END IF
   W( 8 ) = JVS( 19 )
   W( 38 ) = JVS( 20 )
   W( 53 ) = JVS( 21 )
   W( 55 ) = JVS( 22 )
   W( 57 ) = JVS( 23 )
  JVS( 19) = W( 8 )
  JVS( 20) = W( 38 )
  JVS( 21) = W( 53 )
  JVS( 22) = W( 55 )
  JVS( 23) = W( 57 )
  IF ( ABS(  JVS( 24 )) < TINY(a) ) THEN
         IER = 9                                       
         RETURN
  END IF
   W( 9 ) = JVS( 24 )
   W( 51 ) = JVS( 25 )
   W( 55 ) = JVS( 26 )
   W( 57 ) = JVS( 27 )
  JVS( 24) = W( 9 )
  JVS( 25) = W( 51 )
  JVS( 26) = W( 55 )
  JVS( 27) = W( 57 )
  IF ( ABS(  JVS( 28 )) < TINY(a) ) THEN
         IER = 10                                      
         RETURN
  END IF
   W( 10 ) = JVS( 28 )
   W( 47 ) = JVS( 29 )
   W( 50 ) = JVS( 30 )
   W( 51 ) = JVS( 31 )
   W( 55 ) = JVS( 32 )
   W( 57 ) = JVS( 33 )
  JVS( 28) = W( 10 )
  JVS( 29) = W( 47 )
  JVS( 30) = W( 50 )
  JVS( 31) = W( 51 )
  JVS( 32) = W( 55 )
  JVS( 33) = W( 57 )
  IF ( ABS(  JVS( 34 )) < TINY(a) ) THEN
         IER = 11                                      
         RETURN
  END IF
   W( 11 ) = JVS( 34 )
   W( 48 ) = JVS( 35 )
   W( 55 ) = JVS( 36 )
   W( 57 ) = JVS( 37 )
  JVS( 34) = W( 11 )
  JVS( 35) = W( 48 )
  JVS( 36) = W( 55 )
  JVS( 37) = W( 57 )
  IF ( ABS(  JVS( 38 )) < TINY(a) ) THEN
         IER = 12                                      
         RETURN
  END IF
   W( 12 ) = JVS( 38 )
   W( 50 ) = JVS( 39 )
   W( 55 ) = JVS( 40 )
   W( 57 ) = JVS( 41 )
  JVS( 38) = W( 12 )
  JVS( 39) = W( 50 )
  JVS( 40) = W( 55 )
  JVS( 41) = W( 57 )
  IF ( ABS(  JVS( 42 )) < TINY(a) ) THEN
         IER = 13                                      
         RETURN
  END IF
   W( 13 ) = JVS( 42 )
   W( 51 ) = JVS( 43 )
   W( 54 ) = JVS( 44 )
   W( 55 ) = JVS( 45 )
  JVS( 42) = W( 13 )
  JVS( 43) = W( 51 )
  JVS( 44) = W( 54 )
  JVS( 45) = W( 55 )
  IF ( ABS(  JVS( 46 )) < TINY(a) ) THEN
         IER = 14                                      
         RETURN
  END IF
   W( 14 ) = JVS( 46 )
   W( 49 ) = JVS( 47 )
   W( 54 ) = JVS( 48 )
   W( 55 ) = JVS( 49 )
  JVS( 46) = W( 14 )
  JVS( 47) = W( 49 )
  JVS( 48) = W( 54 )
  JVS( 49) = W( 55 )
  IF ( ABS(  JVS( 50 )) < TINY(a) ) THEN
         IER = 15                                      
         RETURN
  END IF
   W( 15 ) = JVS( 50 )
   W( 19 ) = JVS( 51 )
   W( 37 ) = JVS( 52 )
   W( 38 ) = JVS( 53 )
   W( 41 ) = JVS( 54 )
   W( 53 ) = JVS( 55 )
   W( 55 ) = JVS( 56 )
  JVS( 50) = W( 15 )
  JVS( 51) = W( 19 )
  JVS( 52) = W( 37 )
  JVS( 53) = W( 38 )
  JVS( 54) = W( 41 )
  JVS( 55) = W( 53 )
  JVS( 56) = W( 55 )
  IF ( ABS(  JVS( 57 )) < TINY(a) ) THEN
         IER = 16                                      
         RETURN
  END IF
   W( 16 ) = JVS( 57 )
   W( 47 ) = JVS( 58 )
   W( 55 ) = JVS( 59 )
   W( 57 ) = JVS( 60 )
  JVS( 57) = W( 16 )
  JVS( 58) = W( 47 )
  JVS( 59) = W( 55 )
  JVS( 60) = W( 57 )
  IF ( ABS(  JVS( 61 )) < TINY(a) ) THEN
         IER = 17                                      
         RETURN
  END IF
   W( 17 ) = JVS( 61 )
   W( 34 ) = JVS( 62 )
   W( 55 ) = JVS( 63 )
   W( 57 ) = JVS( 64 )
  JVS( 61) = W( 17 )
  JVS( 62) = W( 34 )
  JVS( 63) = W( 55 )
  JVS( 64) = W( 57 )
  IF ( ABS(  JVS( 65 )) < TINY(a) ) THEN
         IER = 18                                      
         RETURN
  END IF
   W( 18 ) = JVS( 65 )
   W( 34 ) = JVS( 66 )
   W( 55 ) = JVS( 67 )
   W( 56 ) = JVS( 68 )
  JVS( 65) = W( 18 )
  JVS( 66) = W( 34 )
  JVS( 67) = W( 55 )
  JVS( 68) = W( 56 )
  IF ( ABS(  JVS( 69 )) < TINY(a) ) THEN
         IER = 19                                      
         RETURN
  END IF
   W( 19 ) = JVS( 69 )
   W( 53 ) = JVS( 70 )
   W( 55 ) = JVS( 71 )
  JVS( 69) = W( 19 )
  JVS( 70) = W( 53 )
  JVS( 71) = W( 55 )
  IF ( ABS(  JVS( 72 )) < TINY(a) ) THEN
         IER = 20                                      
         RETURN
  END IF
   W( 20 ) = JVS( 72 )
   W( 42 ) = JVS( 73 )
   W( 55 ) = JVS( 74 )
   W( 57 ) = JVS( 75 )
  JVS( 72) = W( 20 )
  JVS( 73) = W( 42 )
  JVS( 74) = W( 55 )
  JVS( 75) = W( 57 )
  IF ( ABS(  JVS( 76 )) < TINY(a) ) THEN
         IER = 21                                      
         RETURN
  END IF
   W( 21 ) = JVS( 76 )
   W( 24 ) = JVS( 77 )
   W( 53 ) = JVS( 78 )
  JVS( 76) = W( 21 )
  JVS( 77) = W( 24 )
  JVS( 78) = W( 53 )
  IF ( ABS(  JVS( 79 )) < TINY(a) ) THEN
         IER = 22                                      
         RETURN
  END IF
   W( 22 ) = JVS( 79 )
   W( 54 ) = JVS( 80 )
   W( 55 ) = JVS( 81 )
   W( 57 ) = JVS( 82 )
  JVS( 79) = W( 22 )
  JVS( 80) = W( 54 )
  JVS( 81) = W( 55 )
  JVS( 82) = W( 57 )
  IF ( ABS(  JVS( 83 )) < TINY(a) ) THEN
         IER = 23                                      
         RETURN
  END IF
   W( 23 ) = JVS( 83 )
   W( 36 ) = JVS( 84 )
   W( 55 ) = JVS( 85 )
   W( 57 ) = JVS( 86 )
  JVS( 83) = W( 23 )
  JVS( 84) = W( 36 )
  JVS( 85) = W( 55 )
  JVS( 86) = W( 57 )
  IF ( ABS(  JVS( 88 )) < TINY(a) ) THEN
         IER = 24                                      
         RETURN
  END IF
   W( 21 ) = JVS( 87 )
   W( 24 ) = JVS( 88 )
   W( 37 ) = JVS( 89 )
   W( 53 ) = JVS( 90 )
   W( 55 ) = JVS( 91 )
  a = -W( 21 ) / JVS(           76  )
  W( 21 ) = -a
  W( 24 ) = W( 24 ) + a*JVS( 77 )
  W( 53 ) = W( 53 ) + a*JVS( 78 )
  JVS( 87) = W( 21 )
  JVS( 88) = W( 24 )
  JVS( 89) = W( 37 )
  JVS( 90) = W( 53 )
  JVS( 91) = W( 55 )
  IF ( ABS(  JVS( 93 )) < TINY(a) ) THEN
         IER = 25                                      
         RETURN
  END IF
   W( 5 ) = JVS( 92 )
   W( 25 ) = JVS( 93 )
   W( 45 ) = JVS( 94 )
   W( 46 ) = JVS( 95 )
   W( 52 ) = JVS( 96 )
   W( 54 ) = JVS( 97 )
   W( 55 ) = JVS( 98 )
  a = -W( 5 ) / JVS(           10  )
  W( 5 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 11 )
  W( 54 ) = W( 54 ) + a*JVS( 12 )
  JVS( 92) = W( 5 )
  JVS( 93) = W( 25 )
  JVS( 94) = W( 45 )
  JVS( 95) = W( 46 )
  JVS( 96) = W( 52 )
  JVS( 97) = W( 54 )
  JVS( 98) = W( 55 )
  IF ( ABS(  JVS( 99 )) < TINY(a) ) THEN
         IER = 26                                      
         RETURN
  END IF
   W( 26 ) = JVS( 99 )
   W( 40 ) = JVS( 100 )
   W( 55 ) = JVS( 101 )
   W( 57 ) = JVS( 102 )
  JVS( 99) = W( 26 )
  JVS( 100) = W( 40 )
  JVS( 101) = W( 55 )
  JVS( 102) = W( 57 )
  IF ( ABS(  JVS( 103 )) < TINY(a) ) THEN
         IER = 27                                      
         RETURN
  END IF
   W( 27 ) = JVS( 103 )
   W( 43 ) = JVS( 104 )
   W( 55 ) = JVS( 105 )
   W( 57 ) = JVS( 106 )
  JVS( 103) = W( 27 )
  JVS( 104) = W( 43 )
  JVS( 105) = W( 55 )
  JVS( 106) = W( 57 )
  IF ( ABS(  JVS( 109 )) < TINY(a) ) THEN
         IER = 28                                      
         RETURN
  END IF
   W( 17 ) = JVS( 107 )
   W( 18 ) = JVS( 108 )
   W( 28 ) = JVS( 109 )
   W( 34 ) = JVS( 110 )
   W( 50 ) = JVS( 111 )
   W( 55 ) = JVS( 112 )
   W( 56 ) = JVS( 113 )
   W( 57 ) = JVS( 114 )
  a = -W( 17 ) / JVS(           61  )
  W( 17 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 62 )
  W( 55 ) = W( 55 ) + a*JVS( 63 )
  W( 57 ) = W( 57 ) + a*JVS( 64 )
  a = -W( 18 ) / JVS(           65  )
  W( 18 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 66 )
  W( 55 ) = W( 55 ) + a*JVS( 67 )
  W( 56 ) = W( 56 ) + a*JVS( 68 )
  JVS( 107) = W( 17 )
  JVS( 108) = W( 18 )
  JVS( 109) = W( 28 )
  JVS( 110) = W( 34 )
  JVS( 111) = W( 50 )
  JVS( 112) = W( 55 )
  JVS( 113) = W( 56 )
  JVS( 114) = W( 57 )
  IF ( ABS(  JVS( 115 )) < TINY(a) ) THEN
         IER = 29                                      
         RETURN
  END IF
   W( 29 ) = JVS( 115 )
   W( 49 ) = JVS( 116 )
   W( 55 ) = JVS( 117 )
   W( 57 ) = JVS( 118 )
  JVS( 115) = W( 29 )
  JVS( 116) = W( 49 )
  JVS( 117) = W( 55 )
  JVS( 118) = W( 57 )
  IF ( ABS(  JVS( 120 )) < TINY(a) ) THEN
         IER = 30                                      
         RETURN
  END IF
   W( 27 ) = JVS( 119 )
   W( 30 ) = JVS( 120 )
   W( 43 ) = JVS( 121 )
   W( 44 ) = JVS( 122 )
   W( 50 ) = JVS( 123 )
   W( 55 ) = JVS( 124 )
   W( 56 ) = JVS( 125 )
   W( 57 ) = JVS( 126 )
  a = -W( 27 ) / JVS(          103  )
  W( 27 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 104 )
  W( 55 ) = W( 55 ) + a*JVS( 105 )
  W( 57 ) = W( 57 ) + a*JVS( 106 )
  JVS( 119) = W( 27 )
  JVS( 120) = W( 30 )
  JVS( 121) = W( 43 )
  JVS( 122) = W( 44 )
  JVS( 123) = W( 50 )
  JVS( 124) = W( 55 )
  JVS( 125) = W( 56 )
  JVS( 126) = W( 57 )
  IF ( ABS(  JVS( 130 )) < TINY(a) ) THEN
         IER = 31                                      
         RETURN
  END IF
   W( 7 ) = JVS( 127 )
   W( 19 ) = JVS( 128 )
   W( 29 ) = JVS( 129 )
   W( 31 ) = JVS( 130 )
   W( 33 ) = JVS( 131 )
   W( 35 ) = JVS( 132 )
   W( 37 ) = JVS( 133 )
   W( 38 ) = JVS( 134 )
   W( 41 ) = JVS( 135 )
   W( 45 ) = JVS( 136 )
   W( 46 ) = JVS( 137 )
   W( 49 ) = JVS( 138 )
   W( 50 ) = JVS( 139 )
   W( 52 ) = JVS( 140 )
   W( 53 ) = JVS( 141 )
   W( 55 ) = JVS( 142 )
   W( 56 ) = JVS( 143 )
   W( 57 ) = JVS( 144 )
  a = -W( 7 ) / JVS(           16  )
  W( 7 ) = -a
  W( 33 ) = W( 33 ) + a*JVS( 17 )
  W( 55 ) = W( 55 ) + a*JVS( 18 )
  a = -W( 19 ) / JVS(           69  )
  W( 19 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 70 )
  W( 55 ) = W( 55 ) + a*JVS( 71 )
  a = -W( 29 ) / JVS(          115  )
  W( 29 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 116 )
  W( 55 ) = W( 55 ) + a*JVS( 117 )
  W( 57 ) = W( 57 ) + a*JVS( 118 )
  JVS( 127) = W( 7 )
  JVS( 128) = W( 19 )
  JVS( 129) = W( 29 )
  JVS( 130) = W( 31 )
  JVS( 131) = W( 33 )
  JVS( 132) = W( 35 )
  JVS( 133) = W( 37 )
  JVS( 134) = W( 38 )
  JVS( 135) = W( 41 )
  JVS( 136) = W( 45 )
  JVS( 137) = W( 46 )
  JVS( 138) = W( 49 )
  JVS( 139) = W( 50 )
  JVS( 140) = W( 52 )
  JVS( 141) = W( 53 )
  JVS( 142) = W( 55 )
  JVS( 143) = W( 56 )
  JVS( 144) = W( 57 )
  IF ( ABS(  JVS( 145 )) < TINY(a) ) THEN
         IER = 32                                      
         RETURN
  END IF
   W( 32 ) = JVS( 145 )
   W( 34 ) = JVS( 146 )
   W( 36 ) = JVS( 147 )
   W( 37 ) = JVS( 148 )
   W( 43 ) = JVS( 149 )
   W( 47 ) = JVS( 150 )
   W( 48 ) = JVS( 151 )
   W( 49 ) = JVS( 152 )
   W( 50 ) = JVS( 153 )
   W( 53 ) = JVS( 154 )
   W( 55 ) = JVS( 155 )
  JVS( 145) = W( 32 )
  JVS( 146) = W( 34 )
  JVS( 147) = W( 36 )
  JVS( 148) = W( 37 )
  JVS( 149) = W( 43 )
  JVS( 150) = W( 47 )
  JVS( 151) = W( 48 )
  JVS( 152) = W( 49 )
  JVS( 153) = W( 50 )
  JVS( 154) = W( 53 )
  JVS( 155) = W( 55 )
  IF ( ABS(  JVS( 156 )) < TINY(a) ) THEN
         IER = 33                                      
         RETURN
  END IF
   W( 33 ) = JVS( 156 )
   W( 38 ) = JVS( 157 )
   W( 48 ) = JVS( 158 )
   W( 52 ) = JVS( 159 )
   W( 55 ) = JVS( 160 )
   W( 56 ) = JVS( 161 )
  JVS( 156) = W( 33 )
  JVS( 157) = W( 38 )
  JVS( 158) = W( 48 )
  JVS( 159) = W( 52 )
  JVS( 160) = W( 55 )
  JVS( 161) = W( 56 )
  IF ( ABS(  JVS( 164 )) < TINY(a) ) THEN
         IER = 34                                      
         RETURN
  END IF
   W( 4 ) = JVS( 162 )
   W( 17 ) = JVS( 163 )
   W( 34 ) = JVS( 164 )
   W( 50 ) = JVS( 165 )
   W( 55 ) = JVS( 166 )
   W( 56 ) = JVS( 167 )
   W( 57 ) = JVS( 168 )
  a = -W( 4 ) / JVS(            8  )
  W( 4 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 9 )
  a = -W( 17 ) / JVS(           61  )
  W( 17 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 62 )
  W( 55 ) = W( 55 ) + a*JVS( 63 )
  W( 57 ) = W( 57 ) + a*JVS( 64 )
  JVS( 162) = W( 4 )
  JVS( 163) = W( 17 )
  JVS( 164) = W( 34 )
  JVS( 165) = W( 50 )
  JVS( 166) = W( 55 )
  JVS( 167) = W( 56 )
  JVS( 168) = W( 57 )
  IF ( ABS(  JVS( 171 )) < TINY(a) ) THEN
         IER = 35                                      
         RETURN
  END IF
   W( 23 ) = JVS( 169 )
   W( 29 ) = JVS( 170 )
   W( 35 ) = JVS( 171 )
   W( 36 ) = JVS( 172 )
   W( 37 ) = JVS( 173 )
   W( 39 ) = JVS( 174 )
   W( 41 ) = JVS( 175 )
   W( 48 ) = JVS( 176 )
   W( 49 ) = JVS( 177 )
   W( 50 ) = JVS( 178 )
   W( 53 ) = JVS( 179 )
   W( 55 ) = JVS( 180 )
   W( 56 ) = JVS( 181 )
   W( 57 ) = JVS( 182 )
  a = -W( 23 ) / JVS(           83  )
  W( 23 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 84 )
  W( 55 ) = W( 55 ) + a*JVS( 85 )
  W( 57 ) = W( 57 ) + a*JVS( 86 )
  a = -W( 29 ) / JVS(          115  )
  W( 29 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 116 )
  W( 55 ) = W( 55 ) + a*JVS( 117 )
  W( 57 ) = W( 57 ) + a*JVS( 118 )
  JVS( 169) = W( 23 )
  JVS( 170) = W( 29 )
  JVS( 171) = W( 35 )
  JVS( 172) = W( 36 )
  JVS( 173) = W( 37 )
  JVS( 174) = W( 39 )
  JVS( 175) = W( 41 )
  JVS( 176) = W( 48 )
  JVS( 177) = W( 49 )
  JVS( 178) = W( 50 )
  JVS( 179) = W( 53 )
  JVS( 180) = W( 55 )
  JVS( 181) = W( 56 )
  JVS( 182) = W( 57 )
  IF ( ABS(  JVS( 186 )) < TINY(a) ) THEN
         IER = 36                                      
         RETURN
  END IF
   W( 23 ) = JVS( 183 )
   W( 28 ) = JVS( 184 )
   W( 34 ) = JVS( 185 )
   W( 36 ) = JVS( 186 )
   W( 50 ) = JVS( 187 )
   W( 55 ) = JVS( 188 )
   W( 56 ) = JVS( 189 )
   W( 57 ) = JVS( 190 )
  a = -W( 23 ) / JVS(           83  )
  W( 23 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 84 )
  W( 55 ) = W( 55 ) + a*JVS( 85 )
  W( 57 ) = W( 57 ) + a*JVS( 86 )
  a = -W( 28 ) / JVS(          109  )
  W( 28 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 110 )
  W( 50 ) = W( 50 ) + a*JVS( 111 )
  W( 55 ) = W( 55 ) + a*JVS( 112 )
  W( 56 ) = W( 56 ) + a*JVS( 113 )
  W( 57 ) = W( 57 ) + a*JVS( 114 )
  a = -W( 34 ) / JVS(          164  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 165 )
  W( 55 ) = W( 55 ) + a*JVS( 166 )
  W( 56 ) = W( 56 ) + a*JVS( 167 )
  W( 57 ) = W( 57 ) + a*JVS( 168 )
  JVS( 183) = W( 23 )
  JVS( 184) = W( 28 )
  JVS( 185) = W( 34 )
  JVS( 186) = W( 36 )
  JVS( 187) = W( 50 )
  JVS( 188) = W( 55 )
  JVS( 189) = W( 56 )
  JVS( 190) = W( 57 )
  IF ( ABS(  JVS( 191 )) < TINY(a) ) THEN
         IER = 37                                      
         RETURN
  END IF
   W( 37 ) = JVS( 191 )
   W( 52 ) = JVS( 192 )
   W( 53 ) = JVS( 193 )
   W( 55 ) = JVS( 194 )
  JVS( 191) = W( 37 )
  JVS( 192) = W( 52 )
  JVS( 193) = W( 53 )
  JVS( 194) = W( 55 )
  IF ( ABS(  JVS( 195 )) < TINY(a) ) THEN
         IER = 38                                      
         RETURN
  END IF
   W( 38 ) = JVS( 195 )
   W( 52 ) = JVS( 196 )
   W( 53 ) = JVS( 197 )
   W( 55 ) = JVS( 198 )
  JVS( 195) = W( 38 )
  JVS( 196) = W( 52 )
  JVS( 197) = W( 53 )
  JVS( 198) = W( 55 )
  IF ( ABS(  JVS( 205 )) < TINY(a) ) THEN
         IER = 39                                      
         RETURN
  END IF
   W( 14 ) = JVS( 199 )
   W( 26 ) = JVS( 200 )
   W( 29 ) = JVS( 201 )
   W( 33 ) = JVS( 202 )
   W( 36 ) = JVS( 203 )
   W( 38 ) = JVS( 204 )
   W( 39 ) = JVS( 205 )
   W( 40 ) = JVS( 206 )
   W( 48 ) = JVS( 207 )
   W( 49 ) = JVS( 208 )
   W( 50 ) = JVS( 209 )
   W( 52 ) = JVS( 210 )
   W( 53 ) = JVS( 211 )
   W( 54 ) = JVS( 212 )
   W( 55 ) = JVS( 213 )
   W( 56 ) = JVS( 214 )
   W( 57 ) = JVS( 215 )
  a = -W( 14 ) / JVS(           46  )
  W( 14 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 47 )
  W( 54 ) = W( 54 ) + a*JVS( 48 )
  W( 55 ) = W( 55 ) + a*JVS( 49 )
  a = -W( 26 ) / JVS(           99  )
  W( 26 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 100 )
  W( 55 ) = W( 55 ) + a*JVS( 101 )
  W( 57 ) = W( 57 ) + a*JVS( 102 )
  a = -W( 29 ) / JVS(          115  )
  W( 29 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 116 )
  W( 55 ) = W( 55 ) + a*JVS( 117 )
  W( 57 ) = W( 57 ) + a*JVS( 118 )
  a = -W( 33 ) / JVS(          156  )
  W( 33 ) = -a
  W( 38 ) = W( 38 ) + a*JVS( 157 )
  W( 48 ) = W( 48 ) + a*JVS( 158 )
  W( 52 ) = W( 52 ) + a*JVS( 159 )
  W( 55 ) = W( 55 ) + a*JVS( 160 )
  W( 56 ) = W( 56 ) + a*JVS( 161 )
  a = -W( 36 ) / JVS(          186  )
  W( 36 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 187 )
  W( 55 ) = W( 55 ) + a*JVS( 188 )
  W( 56 ) = W( 56 ) + a*JVS( 189 )
  W( 57 ) = W( 57 ) + a*JVS( 190 )
  a = -W( 38 ) / JVS(          195  )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 196 )
  W( 53 ) = W( 53 ) + a*JVS( 197 )
  W( 55 ) = W( 55 ) + a*JVS( 198 )
  JVS( 199) = W( 14 )
  JVS( 200) = W( 26 )
  JVS( 201) = W( 29 )
  JVS( 202) = W( 33 )
  JVS( 203) = W( 36 )
  JVS( 204) = W( 38 )
  JVS( 205) = W( 39 )
  JVS( 206) = W( 40 )
  JVS( 207) = W( 48 )
  JVS( 208) = W( 49 )
  JVS( 209) = W( 50 )
  JVS( 210) = W( 52 )
  JVS( 211) = W( 53 )
  JVS( 212) = W( 54 )
  JVS( 213) = W( 55 )
  JVS( 214) = W( 56 )
  JVS( 215) = W( 57 )
  IF ( ABS(  JVS( 219 )) < TINY(a) ) THEN
         IER = 40                                      
         RETURN
  END IF
   W( 19 ) = JVS( 216 )
   W( 26 ) = JVS( 217 )
   W( 37 ) = JVS( 218 )
   W( 40 ) = JVS( 219 )
   W( 52 ) = JVS( 220 )
   W( 53 ) = JVS( 221 )
   W( 55 ) = JVS( 222 )
   W( 56 ) = JVS( 223 )
   W( 57 ) = JVS( 224 )
  a = -W( 19 ) / JVS(           69  )
  W( 19 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 70 )
  W( 55 ) = W( 55 ) + a*JVS( 71 )
  a = -W( 26 ) / JVS(           99  )
  W( 26 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 100 )
  W( 55 ) = W( 55 ) + a*JVS( 101 )
  W( 57 ) = W( 57 ) + a*JVS( 102 )
  a = -W( 37 ) / JVS(          191  )
  W( 37 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 192 )
  W( 53 ) = W( 53 ) + a*JVS( 193 )
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  JVS( 216) = W( 19 )
  JVS( 217) = W( 26 )
  JVS( 218) = W( 37 )
  JVS( 219) = W( 40 )
  JVS( 220) = W( 52 )
  JVS( 221) = W( 53 )
  JVS( 222) = W( 55 )
  JVS( 223) = W( 56 )
  JVS( 224) = W( 57 )
  IF ( ABS(  JVS( 228 )) < TINY(a) ) THEN
         IER = 41                                      
         RETURN
  END IF
   W( 11 ) = JVS( 225 )
   W( 33 ) = JVS( 226 )
   W( 38 ) = JVS( 227 )
   W( 41 ) = JVS( 228 )
   W( 48 ) = JVS( 229 )
   W( 50 ) = JVS( 230 )
   W( 52 ) = JVS( 231 )
   W( 53 ) = JVS( 232 )
   W( 55 ) = JVS( 233 )
   W( 56 ) = JVS( 234 )
   W( 57 ) = JVS( 235 )
  a = -W( 11 ) / JVS(           34  )
  W( 11 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 35 )
  W( 55 ) = W( 55 ) + a*JVS( 36 )
  W( 57 ) = W( 57 ) + a*JVS( 37 )
  a = -W( 33 ) / JVS(          156  )
  W( 33 ) = -a
  W( 38 ) = W( 38 ) + a*JVS( 157 )
  W( 48 ) = W( 48 ) + a*JVS( 158 )
  W( 52 ) = W( 52 ) + a*JVS( 159 )
  W( 55 ) = W( 55 ) + a*JVS( 160 )
  W( 56 ) = W( 56 ) + a*JVS( 161 )
  a = -W( 38 ) / JVS(          195  )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 196 )
  W( 53 ) = W( 53 ) + a*JVS( 197 )
  W( 55 ) = W( 55 ) + a*JVS( 198 )
  JVS( 225) = W( 11 )
  JVS( 226) = W( 33 )
  JVS( 227) = W( 38 )
  JVS( 228) = W( 41 )
  JVS( 229) = W( 48 )
  JVS( 230) = W( 50 )
  JVS( 231) = W( 52 )
  JVS( 232) = W( 53 )
  JVS( 233) = W( 55 )
  JVS( 234) = W( 56 )
  JVS( 235) = W( 57 )
  IF ( ABS(  JVS( 238 )) < TINY(a) ) THEN
         IER = 42                                      
         RETURN
  END IF
   W( 20 ) = JVS( 236 )
   W( 30 ) = JVS( 237 )
   W( 42 ) = JVS( 238 )
   W( 43 ) = JVS( 239 )
   W( 44 ) = JVS( 240 )
   W( 50 ) = JVS( 241 )
   W( 55 ) = JVS( 242 )
   W( 56 ) = JVS( 243 )
   W( 57 ) = JVS( 244 )
  a = -W( 20 ) / JVS(           72  )
  W( 20 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 73 )
  W( 55 ) = W( 55 ) + a*JVS( 74 )
  W( 57 ) = W( 57 ) + a*JVS( 75 )
  a = -W( 30 ) / JVS(          120  )
  W( 30 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 121 )
  W( 44 ) = W( 44 ) + a*JVS( 122 )
  W( 50 ) = W( 50 ) + a*JVS( 123 )
  W( 55 ) = W( 55 ) + a*JVS( 124 )
  W( 56 ) = W( 56 ) + a*JVS( 125 )
  W( 57 ) = W( 57 ) + a*JVS( 126 )
  JVS( 236) = W( 20 )
  JVS( 237) = W( 30 )
  JVS( 238) = W( 42 )
  JVS( 239) = W( 43 )
  JVS( 240) = W( 44 )
  JVS( 241) = W( 50 )
  JVS( 242) = W( 55 )
  JVS( 243) = W( 56 )
  JVS( 244) = W( 57 )
  IF ( ABS(  JVS( 247 )) < TINY(a) ) THEN
         IER = 43                                      
         RETURN
  END IF
   W( 2 ) = JVS( 245 )
   W( 27 ) = JVS( 246 )
   W( 43 ) = JVS( 247 )
   W( 50 ) = JVS( 248 )
   W( 55 ) = JVS( 249 )
   W( 56 ) = JVS( 250 )
   W( 57 ) = JVS( 251 )
  a = -W( 2 ) / JVS(            3  )
  W( 2 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 4 )
  a = -W( 27 ) / JVS(          103  )
  W( 27 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 104 )
  W( 55 ) = W( 55 ) + a*JVS( 105 )
  W( 57 ) = W( 57 ) + a*JVS( 106 )
  JVS( 245) = W( 2 )
  JVS( 246) = W( 27 )
  JVS( 247) = W( 43 )
  JVS( 248) = W( 50 )
  JVS( 249) = W( 55 )
  JVS( 250) = W( 56 )
  JVS( 251) = W( 57 )
  IF ( ABS(  JVS( 256 )) < TINY(a) ) THEN
         IER = 44                                      
         RETURN
  END IF
   W( 37 ) = JVS( 252 )
   W( 40 ) = JVS( 253 )
   W( 42 ) = JVS( 254 )
   W( 43 ) = JVS( 255 )
   W( 44 ) = JVS( 256 )
   W( 50 ) = JVS( 257 )
   W( 52 ) = JVS( 258 )
   W( 53 ) = JVS( 259 )
   W( 55 ) = JVS( 260 )
   W( 56 ) = JVS( 261 )
   W( 57 ) = JVS( 262 )
  a = -W( 37 ) / JVS(          191  )
  W( 37 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 192 )
  W( 53 ) = W( 53 ) + a*JVS( 193 )
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  a = -W( 40 ) / JVS(          219  )
  W( 40 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 220 )
  W( 53 ) = W( 53 ) + a*JVS( 221 )
  W( 55 ) = W( 55 ) + a*JVS( 222 )
  W( 56 ) = W( 56 ) + a*JVS( 223 )
  W( 57 ) = W( 57 ) + a*JVS( 224 )
  a = -W( 42 ) / JVS(          238  )
  W( 42 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 239 )
  W( 44 ) = W( 44 ) + a*JVS( 240 )
  W( 50 ) = W( 50 ) + a*JVS( 241 )
  W( 55 ) = W( 55 ) + a*JVS( 242 )
  W( 56 ) = W( 56 ) + a*JVS( 243 )
  W( 57 ) = W( 57 ) + a*JVS( 244 )
  a = -W( 43 ) / JVS(          247  )
  W( 43 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 248 )
  W( 55 ) = W( 55 ) + a*JVS( 249 )
  W( 56 ) = W( 56 ) + a*JVS( 250 )
  W( 57 ) = W( 57 ) + a*JVS( 251 )
  JVS( 252) = W( 37 )
  JVS( 253) = W( 40 )
  JVS( 254) = W( 42 )
  JVS( 255) = W( 43 )
  JVS( 256) = W( 44 )
  JVS( 257) = W( 50 )
  JVS( 258) = W( 52 )
  JVS( 259) = W( 53 )
  JVS( 260) = W( 55 )
  JVS( 261) = W( 56 )
  JVS( 262) = W( 57 )
  IF ( ABS(  JVS( 283 )) < TINY(a) ) THEN
         IER = 45                                      
         RETURN
  END IF
   W( 6 ) = JVS( 263 )
   W( 7 ) = JVS( 264 )
   W( 11 ) = JVS( 265 )
   W( 12 ) = JVS( 266 )
   W( 13 ) = JVS( 267 )
   W( 19 ) = JVS( 268 )
   W( 21 ) = JVS( 269 )
   W( 24 ) = JVS( 270 )
   W( 26 ) = JVS( 271 )
   W( 29 ) = JVS( 272 )
   W( 32 ) = JVS( 273 )
   W( 33 ) = JVS( 274 )
   W( 34 ) = JVS( 275 )
   W( 36 ) = JVS( 276 )
   W( 37 ) = JVS( 277 )
   W( 38 ) = JVS( 278 )
   W( 39 ) = JVS( 279 )
   W( 40 ) = JVS( 280 )
   W( 41 ) = JVS( 281 )
   W( 43 ) = JVS( 282 )
   W( 45 ) = JVS( 283 )
   W( 47 ) = JVS( 284 )
   W( 48 ) = JVS( 285 )
   W( 49 ) = JVS( 286 )
   W( 50 ) = JVS( 287 )
   W( 51 ) = JVS( 288 )
   W( 52 ) = JVS( 289 )
   W( 53 ) = JVS( 290 )
   W( 54 ) = JVS( 291 )
   W( 55 ) = JVS( 292 )
   W( 56 ) = JVS( 293 )
   W( 57 ) = JVS( 294 )
  a = -W( 6 ) / JVS(           13  )
  W( 6 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 14 )
  W( 54 ) = W( 54 ) + a*JVS( 15 )
  a = -W( 7 ) / JVS(           16  )
  W( 7 ) = -a
  W( 33 ) = W( 33 ) + a*JVS( 17 )
  W( 55 ) = W( 55 ) + a*JVS( 18 )
  a = -W( 11 ) / JVS(           34  )
  W( 11 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 35 )
  W( 55 ) = W( 55 ) + a*JVS( 36 )
  W( 57 ) = W( 57 ) + a*JVS( 37 )
  a = -W( 12 ) / JVS(           38  )
  W( 12 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 39 )
  W( 55 ) = W( 55 ) + a*JVS( 40 )
  W( 57 ) = W( 57 ) + a*JVS( 41 )
  a = -W( 13 ) / JVS(           42  )
  W( 13 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 43 )
  W( 54 ) = W( 54 ) + a*JVS( 44 )
  W( 55 ) = W( 55 ) + a*JVS( 45 )
  a = -W( 19 ) / JVS(           69  )
  W( 19 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 70 )
  W( 55 ) = W( 55 ) + a*JVS( 71 )
  a = -W( 21 ) / JVS(           76  )
  W( 21 ) = -a
  W( 24 ) = W( 24 ) + a*JVS( 77 )
  W( 53 ) = W( 53 ) + a*JVS( 78 )
  a = -W( 24 ) / JVS(           88  )
  W( 24 ) = -a
  W( 37 ) = W( 37 ) + a*JVS( 89 )
  W( 53 ) = W( 53 ) + a*JVS( 90 )
  W( 55 ) = W( 55 ) + a*JVS( 91 )
  a = -W( 26 ) / JVS(           99  )
  W( 26 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 100 )
  W( 55 ) = W( 55 ) + a*JVS( 101 )
  W( 57 ) = W( 57 ) + a*JVS( 102 )
  a = -W( 29 ) / JVS(          115  )
  W( 29 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 116 )
  W( 55 ) = W( 55 ) + a*JVS( 117 )
  W( 57 ) = W( 57 ) + a*JVS( 118 )
  a = -W( 32 ) / JVS(          145  )
  W( 32 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 146 )
  W( 36 ) = W( 36 ) + a*JVS( 147 )
  W( 37 ) = W( 37 ) + a*JVS( 148 )
  W( 43 ) = W( 43 ) + a*JVS( 149 )
  W( 47 ) = W( 47 ) + a*JVS( 150 )
  W( 48 ) = W( 48 ) + a*JVS( 151 )
  W( 49 ) = W( 49 ) + a*JVS( 152 )
  W( 50 ) = W( 50 ) + a*JVS( 153 )
  W( 53 ) = W( 53 ) + a*JVS( 154 )
  W( 55 ) = W( 55 ) + a*JVS( 155 )
  a = -W( 33 ) / JVS(          156  )
  W( 33 ) = -a
  W( 38 ) = W( 38 ) + a*JVS( 157 )
  W( 48 ) = W( 48 ) + a*JVS( 158 )
  W( 52 ) = W( 52 ) + a*JVS( 159 )
  W( 55 ) = W( 55 ) + a*JVS( 160 )
  W( 56 ) = W( 56 ) + a*JVS( 161 )
  a = -W( 34 ) / JVS(          164  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 165 )
  W( 55 ) = W( 55 ) + a*JVS( 166 )
  W( 56 ) = W( 56 ) + a*JVS( 167 )
  W( 57 ) = W( 57 ) + a*JVS( 168 )
  a = -W( 36 ) / JVS(          186  )
  W( 36 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 187 )
  W( 55 ) = W( 55 ) + a*JVS( 188 )
  W( 56 ) = W( 56 ) + a*JVS( 189 )
  W( 57 ) = W( 57 ) + a*JVS( 190 )
  a = -W( 37 ) / JVS(          191  )
  W( 37 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 192 )
  W( 53 ) = W( 53 ) + a*JVS( 193 )
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  a = -W( 38 ) / JVS(          195  )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 196 )
  W( 53 ) = W( 53 ) + a*JVS( 197 )
  W( 55 ) = W( 55 ) + a*JVS( 198 )
  a = -W( 39 ) / JVS(          205  )
  W( 39 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 206 )
  W( 48 ) = W( 48 ) + a*JVS( 207 )
  W( 49 ) = W( 49 ) + a*JVS( 208 )
  W( 50 ) = W( 50 ) + a*JVS( 209 )
  W( 52 ) = W( 52 ) + a*JVS( 210 )
  W( 53 ) = W( 53 ) + a*JVS( 211 )
  W( 54 ) = W( 54 ) + a*JVS( 212 )
  W( 55 ) = W( 55 ) + a*JVS( 213 )
  W( 56 ) = W( 56 ) + a*JVS( 214 )
  W( 57 ) = W( 57 ) + a*JVS( 215 )
  a = -W( 40 ) / JVS(          219  )
  W( 40 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 220 )
  W( 53 ) = W( 53 ) + a*JVS( 221 )
  W( 55 ) = W( 55 ) + a*JVS( 222 )
  W( 56 ) = W( 56 ) + a*JVS( 223 )
  W( 57 ) = W( 57 ) + a*JVS( 224 )
  a = -W( 41 ) / JVS(          228  )
  W( 41 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 229 )
  W( 50 ) = W( 50 ) + a*JVS( 230 )
  W( 52 ) = W( 52 ) + a*JVS( 231 )
  W( 53 ) = W( 53 ) + a*JVS( 232 )
  W( 55 ) = W( 55 ) + a*JVS( 233 )
  W( 56 ) = W( 56 ) + a*JVS( 234 )
  W( 57 ) = W( 57 ) + a*JVS( 235 )
  a = -W( 43 ) / JVS(          247  )
  W( 43 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 248 )
  W( 55 ) = W( 55 ) + a*JVS( 249 )
  W( 56 ) = W( 56 ) + a*JVS( 250 )
  W( 57 ) = W( 57 ) + a*JVS( 251 )
  JVS( 263) = W( 6 )
  JVS( 264) = W( 7 )
  JVS( 265) = W( 11 )
  JVS( 266) = W( 12 )
  JVS( 267) = W( 13 )
  JVS( 268) = W( 19 )
  JVS( 269) = W( 21 )
  JVS( 270) = W( 24 )
  JVS( 271) = W( 26 )
  JVS( 272) = W( 29 )
  JVS( 273) = W( 32 )
  JVS( 274) = W( 33 )
  JVS( 275) = W( 34 )
  JVS( 276) = W( 36 )
  JVS( 277) = W( 37 )
  JVS( 278) = W( 38 )
  JVS( 279) = W( 39 )
  JVS( 280) = W( 40 )
  JVS( 281) = W( 41 )
  JVS( 282) = W( 43 )
  JVS( 283) = W( 45 )
  JVS( 284) = W( 47 )
  JVS( 285) = W( 48 )
  JVS( 286) = W( 49 )
  JVS( 287) = W( 50 )
  JVS( 288) = W( 51 )
  JVS( 289) = W( 52 )
  JVS( 290) = W( 53 )
  JVS( 291) = W( 54 )
  JVS( 292) = W( 55 )
  JVS( 293) = W( 56 )
  JVS( 294) = W( 57 )
  IF ( ABS(  JVS( 304 )) < TINY(a) ) THEN
         IER = 46                                      
         RETURN
  END IF
   W( 16 ) = JVS( 295 )
   W( 20 ) = JVS( 296 )
   W( 26 ) = JVS( 297 )
   W( 27 ) = JVS( 298 )
   W( 37 ) = JVS( 299 )
   W( 40 ) = JVS( 300 )
   W( 42 ) = JVS( 301 )
   W( 43 ) = JVS( 302 )
   W( 44 ) = JVS( 303 )
   W( 46 ) = JVS( 304 )
   W( 47 ) = JVS( 305 )
   W( 50 ) = JVS( 306 )
   W( 51 ) = JVS( 307 )
   W( 52 ) = JVS( 308 )
   W( 53 ) = JVS( 309 )
   W( 55 ) = JVS( 310 )
   W( 56 ) = JVS( 311 )
   W( 57 ) = JVS( 312 )
  a = -W( 16 ) / JVS(           57  )
  W( 16 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 58 )
  W( 55 ) = W( 55 ) + a*JVS( 59 )
  W( 57 ) = W( 57 ) + a*JVS( 60 )
  a = -W( 20 ) / JVS(           72  )
  W( 20 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 73 )
  W( 55 ) = W( 55 ) + a*JVS( 74 )
  W( 57 ) = W( 57 ) + a*JVS( 75 )
  a = -W( 26 ) / JVS(           99  )
  W( 26 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 100 )
  W( 55 ) = W( 55 ) + a*JVS( 101 )
  W( 57 ) = W( 57 ) + a*JVS( 102 )
  a = -W( 27 ) / JVS(          103  )
  W( 27 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 104 )
  W( 55 ) = W( 55 ) + a*JVS( 105 )
  W( 57 ) = W( 57 ) + a*JVS( 106 )
  a = -W( 37 ) / JVS(          191  )
  W( 37 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 192 )
  W( 53 ) = W( 53 ) + a*JVS( 193 )
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  a = -W( 40 ) / JVS(          219  )
  W( 40 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 220 )
  W( 53 ) = W( 53 ) + a*JVS( 221 )
  W( 55 ) = W( 55 ) + a*JVS( 222 )
  W( 56 ) = W( 56 ) + a*JVS( 223 )
  W( 57 ) = W( 57 ) + a*JVS( 224 )
  a = -W( 42 ) / JVS(          238  )
  W( 42 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 239 )
  W( 44 ) = W( 44 ) + a*JVS( 240 )
  W( 50 ) = W( 50 ) + a*JVS( 241 )
  W( 55 ) = W( 55 ) + a*JVS( 242 )
  W( 56 ) = W( 56 ) + a*JVS( 243 )
  W( 57 ) = W( 57 ) + a*JVS( 244 )
  a = -W( 43 ) / JVS(          247  )
  W( 43 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 248 )
  W( 55 ) = W( 55 ) + a*JVS( 249 )
  W( 56 ) = W( 56 ) + a*JVS( 250 )
  W( 57 ) = W( 57 ) + a*JVS( 251 )
  a = -W( 44 ) / JVS(          256  )
  W( 44 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 257 )
  W( 52 ) = W( 52 ) + a*JVS( 258 )
  W( 53 ) = W( 53 ) + a*JVS( 259 )
  W( 55 ) = W( 55 ) + a*JVS( 260 )
  W( 56 ) = W( 56 ) + a*JVS( 261 )
  W( 57 ) = W( 57 ) + a*JVS( 262 )
  JVS( 295) = W( 16 )
  JVS( 296) = W( 20 )
  JVS( 297) = W( 26 )
  JVS( 298) = W( 27 )
  JVS( 299) = W( 37 )
  JVS( 300) = W( 40 )
  JVS( 301) = W( 42 )
  JVS( 302) = W( 43 )
  JVS( 303) = W( 44 )
  JVS( 304) = W( 46 )
  JVS( 305) = W( 47 )
  JVS( 306) = W( 50 )
  JVS( 307) = W( 51 )
  JVS( 308) = W( 52 )
  JVS( 309) = W( 53 )
  JVS( 310) = W( 55 )
  JVS( 311) = W( 56 )
  JVS( 312) = W( 57 )
  IF ( ABS(  JVS( 321 )) < TINY(a) ) THEN
         IER = 47                                      
         RETURN
  END IF
   W( 1 ) = JVS( 313 )
   W( 4 ) = JVS( 314 )
   W( 16 ) = JVS( 315 )
   W( 27 ) = JVS( 316 )
   W( 30 ) = JVS( 317 )
   W( 37 ) = JVS( 318 )
   W( 43 ) = JVS( 319 )
   W( 44 ) = JVS( 320 )
   W( 47 ) = JVS( 321 )
   W( 50 ) = JVS( 322 )
   W( 51 ) = JVS( 323 )
   W( 52 ) = JVS( 324 )
   W( 53 ) = JVS( 325 )
   W( 55 ) = JVS( 326 )
   W( 56 ) = JVS( 327 )
   W( 57 ) = JVS( 328 )
  a = -W( 1 ) / JVS(            1  )
  W( 1 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 2 )
  a = -W( 4 ) / JVS(            8  )
  W( 4 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 9 )
  a = -W( 16 ) / JVS(           57  )
  W( 16 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 58 )
  W( 55 ) = W( 55 ) + a*JVS( 59 )
  W( 57 ) = W( 57 ) + a*JVS( 60 )
  a = -W( 27 ) / JVS(          103  )
  W( 27 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 104 )
  W( 55 ) = W( 55 ) + a*JVS( 105 )
  W( 57 ) = W( 57 ) + a*JVS( 106 )
  a = -W( 30 ) / JVS(          120  )
  W( 30 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 121 )
  W( 44 ) = W( 44 ) + a*JVS( 122 )
  W( 50 ) = W( 50 ) + a*JVS( 123 )
  W( 55 ) = W( 55 ) + a*JVS( 124 )
  W( 56 ) = W( 56 ) + a*JVS( 125 )
  W( 57 ) = W( 57 ) + a*JVS( 126 )
  a = -W( 37 ) / JVS(          191  )
  W( 37 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 192 )
  W( 53 ) = W( 53 ) + a*JVS( 193 )
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  a = -W( 43 ) / JVS(          247  )
  W( 43 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 248 )
  W( 55 ) = W( 55 ) + a*JVS( 249 )
  W( 56 ) = W( 56 ) + a*JVS( 250 )
  W( 57 ) = W( 57 ) + a*JVS( 251 )
  a = -W( 44 ) / JVS(          256  )
  W( 44 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 257 )
  W( 52 ) = W( 52 ) + a*JVS( 258 )
  W( 53 ) = W( 53 ) + a*JVS( 259 )
  W( 55 ) = W( 55 ) + a*JVS( 260 )
  W( 56 ) = W( 56 ) + a*JVS( 261 )
  W( 57 ) = W( 57 ) + a*JVS( 262 )
  JVS( 313) = W( 1 )
  JVS( 314) = W( 4 )
  JVS( 315) = W( 16 )
  JVS( 316) = W( 27 )
  JVS( 317) = W( 30 )
  JVS( 318) = W( 37 )
  JVS( 319) = W( 43 )
  JVS( 320) = W( 44 )
  JVS( 321) = W( 47 )
  JVS( 322) = W( 50 )
  JVS( 323) = W( 51 )
  JVS( 324) = W( 52 )
  JVS( 325) = W( 53 )
  JVS( 326) = W( 55 )
  JVS( 327) = W( 56 )
  JVS( 328) = W( 57 )
  IF ( ABS(  JVS( 330 )) < TINY(a) ) THEN
         IER = 48                                      
         RETURN
  END IF
   W( 38 ) = JVS( 329 )
   W( 48 ) = JVS( 330 )
   W( 50 ) = JVS( 331 )
   W( 52 ) = JVS( 332 )
   W( 53 ) = JVS( 333 )
   W( 55 ) = JVS( 334 )
   W( 56 ) = JVS( 335 )
   W( 57 ) = JVS( 336 )
  a = -W( 38 ) / JVS(          195  )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 196 )
  W( 53 ) = W( 53 ) + a*JVS( 197 )
  W( 55 ) = W( 55 ) + a*JVS( 198 )
  JVS( 329) = W( 38 )
  JVS( 330) = W( 48 )
  JVS( 331) = W( 50 )
  JVS( 332) = W( 52 )
  JVS( 333) = W( 53 )
  JVS( 334) = W( 55 )
  JVS( 335) = W( 56 )
  JVS( 336) = W( 57 )
  IF ( ABS(  JVS( 342 )) < TINY(a) ) THEN
         IER = 49                                      
         RETURN
  END IF
   W( 14 ) = JVS( 337 )
   W( 29 ) = JVS( 338 )
   W( 38 ) = JVS( 339 )
   W( 41 ) = JVS( 340 )
   W( 48 ) = JVS( 341 )
   W( 49 ) = JVS( 342 )
   W( 50 ) = JVS( 343 )
   W( 52 ) = JVS( 344 )
   W( 53 ) = JVS( 345 )
   W( 54 ) = JVS( 346 )
   W( 55 ) = JVS( 347 )
   W( 56 ) = JVS( 348 )
   W( 57 ) = JVS( 349 )
  a = -W( 14 ) / JVS(           46  )
  W( 14 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 47 )
  W( 54 ) = W( 54 ) + a*JVS( 48 )
  W( 55 ) = W( 55 ) + a*JVS( 49 )
  a = -W( 29 ) / JVS(          115  )
  W( 29 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 116 )
  W( 55 ) = W( 55 ) + a*JVS( 117 )
  W( 57 ) = W( 57 ) + a*JVS( 118 )
  a = -W( 38 ) / JVS(          195  )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 196 )
  W( 53 ) = W( 53 ) + a*JVS( 197 )
  W( 55 ) = W( 55 ) + a*JVS( 198 )
  a = -W( 41 ) / JVS(          228  )
  W( 41 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 229 )
  W( 50 ) = W( 50 ) + a*JVS( 230 )
  W( 52 ) = W( 52 ) + a*JVS( 231 )
  W( 53 ) = W( 53 ) + a*JVS( 232 )
  W( 55 ) = W( 55 ) + a*JVS( 233 )
  W( 56 ) = W( 56 ) + a*JVS( 234 )
  W( 57 ) = W( 57 ) + a*JVS( 235 )
  a = -W( 48 ) / JVS(          330  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 331 )
  W( 52 ) = W( 52 ) + a*JVS( 332 )
  W( 53 ) = W( 53 ) + a*JVS( 333 )
  W( 55 ) = W( 55 ) + a*JVS( 334 )
  W( 56 ) = W( 56 ) + a*JVS( 335 )
  W( 57 ) = W( 57 ) + a*JVS( 336 )
  JVS( 337) = W( 14 )
  JVS( 338) = W( 29 )
  JVS( 339) = W( 38 )
  JVS( 340) = W( 41 )
  JVS( 341) = W( 48 )
  JVS( 342) = W( 49 )
  JVS( 343) = W( 50 )
  JVS( 344) = W( 52 )
  JVS( 345) = W( 53 )
  JVS( 346) = W( 54 )
  JVS( 347) = W( 55 )
  JVS( 348) = W( 56 )
  JVS( 349) = W( 57 )
  IF ( ABS(  JVS( 366 )) < TINY(a) ) THEN
         IER = 50                                      
         RETURN
  END IF
   W( 6 ) = JVS( 350 )
   W( 9 ) = JVS( 351 )
   W( 10 ) = JVS( 352 )
   W( 12 ) = JVS( 353 )
   W( 21 ) = JVS( 354 )
   W( 24 ) = JVS( 355 )
   W( 28 ) = JVS( 356 )
   W( 34 ) = JVS( 357 )
   W( 36 ) = JVS( 358 )
   W( 37 ) = JVS( 359 )
   W( 38 ) = JVS( 360 )
   W( 43 ) = JVS( 361 )
   W( 46 ) = JVS( 362 )
   W( 47 ) = JVS( 363 )
   W( 48 ) = JVS( 364 )
   W( 49 ) = JVS( 365 )
   W( 50 ) = JVS( 366 )
   W( 51 ) = JVS( 367 )
   W( 52 ) = JVS( 368 )
   W( 53 ) = JVS( 369 )
   W( 54 ) = JVS( 370 )
   W( 55 ) = JVS( 371 )
   W( 56 ) = JVS( 372 )
   W( 57 ) = JVS( 373 )
  a = -W( 6 ) / JVS(           13  )
  W( 6 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 14 )
  W( 54 ) = W( 54 ) + a*JVS( 15 )
  a = -W( 9 ) / JVS(           24  )
  W( 9 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 25 )
  W( 55 ) = W( 55 ) + a*JVS( 26 )
  W( 57 ) = W( 57 ) + a*JVS( 27 )
  a = -W( 10 ) / JVS(           28  )
  W( 10 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 29 )
  W( 50 ) = W( 50 ) + a*JVS( 30 )
  W( 51 ) = W( 51 ) + a*JVS( 31 )
  W( 55 ) = W( 55 ) + a*JVS( 32 )
  W( 57 ) = W( 57 ) + a*JVS( 33 )
  a = -W( 12 ) / JVS(           38  )
  W( 12 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 39 )
  W( 55 ) = W( 55 ) + a*JVS( 40 )
  W( 57 ) = W( 57 ) + a*JVS( 41 )
  a = -W( 21 ) / JVS(           76  )
  W( 21 ) = -a
  W( 24 ) = W( 24 ) + a*JVS( 77 )
  W( 53 ) = W( 53 ) + a*JVS( 78 )
  a = -W( 24 ) / JVS(           88  )
  W( 24 ) = -a
  W( 37 ) = W( 37 ) + a*JVS( 89 )
  W( 53 ) = W( 53 ) + a*JVS( 90 )
  W( 55 ) = W( 55 ) + a*JVS( 91 )
  a = -W( 28 ) / JVS(          109  )
  W( 28 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 110 )
  W( 50 ) = W( 50 ) + a*JVS( 111 )
  W( 55 ) = W( 55 ) + a*JVS( 112 )
  W( 56 ) = W( 56 ) + a*JVS( 113 )
  W( 57 ) = W( 57 ) + a*JVS( 114 )
  a = -W( 34 ) / JVS(          164  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 165 )
  W( 55 ) = W( 55 ) + a*JVS( 166 )
  W( 56 ) = W( 56 ) + a*JVS( 167 )
  W( 57 ) = W( 57 ) + a*JVS( 168 )
  a = -W( 36 ) / JVS(          186  )
  W( 36 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 187 )
  W( 55 ) = W( 55 ) + a*JVS( 188 )
  W( 56 ) = W( 56 ) + a*JVS( 189 )
  W( 57 ) = W( 57 ) + a*JVS( 190 )
  a = -W( 37 ) / JVS(          191  )
  W( 37 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 192 )
  W( 53 ) = W( 53 ) + a*JVS( 193 )
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  a = -W( 38 ) / JVS(          195  )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 196 )
  W( 53 ) = W( 53 ) + a*JVS( 197 )
  W( 55 ) = W( 55 ) + a*JVS( 198 )
  a = -W( 43 ) / JVS(          247  )
  W( 43 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 248 )
  W( 55 ) = W( 55 ) + a*JVS( 249 )
  W( 56 ) = W( 56 ) + a*JVS( 250 )
  W( 57 ) = W( 57 ) + a*JVS( 251 )
  a = -W( 46 ) / JVS(          304  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 305 )
  W( 50 ) = W( 50 ) + a*JVS( 306 )
  W( 51 ) = W( 51 ) + a*JVS( 307 )
  W( 52 ) = W( 52 ) + a*JVS( 308 )
  W( 53 ) = W( 53 ) + a*JVS( 309 )
  W( 55 ) = W( 55 ) + a*JVS( 310 )
  W( 56 ) = W( 56 ) + a*JVS( 311 )
  W( 57 ) = W( 57 ) + a*JVS( 312 )
  a = -W( 47 ) / JVS(          321  )
  W( 47 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 322 )
  W( 51 ) = W( 51 ) + a*JVS( 323 )
  W( 52 ) = W( 52 ) + a*JVS( 324 )
  W( 53 ) = W( 53 ) + a*JVS( 325 )
  W( 55 ) = W( 55 ) + a*JVS( 326 )
  W( 56 ) = W( 56 ) + a*JVS( 327 )
  W( 57 ) = W( 57 ) + a*JVS( 328 )
  a = -W( 48 ) / JVS(          330  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 331 )
  W( 52 ) = W( 52 ) + a*JVS( 332 )
  W( 53 ) = W( 53 ) + a*JVS( 333 )
  W( 55 ) = W( 55 ) + a*JVS( 334 )
  W( 56 ) = W( 56 ) + a*JVS( 335 )
  W( 57 ) = W( 57 ) + a*JVS( 336 )
  a = -W( 49 ) / JVS(          342  )
  W( 49 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 343 )
  W( 52 ) = W( 52 ) + a*JVS( 344 )
  W( 53 ) = W( 53 ) + a*JVS( 345 )
  W( 54 ) = W( 54 ) + a*JVS( 346 )
  W( 55 ) = W( 55 ) + a*JVS( 347 )
  W( 56 ) = W( 56 ) + a*JVS( 348 )
  W( 57 ) = W( 57 ) + a*JVS( 349 )
  JVS( 350) = W( 6 )
  JVS( 351) = W( 9 )
  JVS( 352) = W( 10 )
  JVS( 353) = W( 12 )
  JVS( 354) = W( 21 )
  JVS( 355) = W( 24 )
  JVS( 356) = W( 28 )
  JVS( 357) = W( 34 )
  JVS( 358) = W( 36 )
  JVS( 359) = W( 37 )
  JVS( 360) = W( 38 )
  JVS( 361) = W( 43 )
  JVS( 362) = W( 46 )
  JVS( 363) = W( 47 )
  JVS( 364) = W( 48 )
  JVS( 365) = W( 49 )
  JVS( 366) = W( 50 )
  JVS( 367) = W( 51 )
  JVS( 368) = W( 52 )
  JVS( 369) = W( 53 )
  JVS( 370) = W( 54 )
  JVS( 371) = W( 55 )
  JVS( 372) = W( 56 )
  JVS( 373) = W( 57 )
  IF ( ABS(  JVS( 398 )) < TINY(a) ) THEN
         IER = 51                                      
         RETURN
  END IF
   W( 3 ) = JVS( 374 )
   W( 9 ) = JVS( 375 )
   W( 13 ) = JVS( 376 )
   W( 20 ) = JVS( 377 )
   W( 23 ) = JVS( 378 )
   W( 28 ) = JVS( 379 )
   W( 29 ) = JVS( 380 )
   W( 30 ) = JVS( 381 )
   W( 34 ) = JVS( 382 )
   W( 35 ) = JVS( 383 )
   W( 36 ) = JVS( 384 )
   W( 37 ) = JVS( 385 )
   W( 38 ) = JVS( 386 )
   W( 39 ) = JVS( 387 )
   W( 40 ) = JVS( 388 )
   W( 41 ) = JVS( 389 )
   W( 42 ) = JVS( 390 )
   W( 43 ) = JVS( 391 )
   W( 44 ) = JVS( 392 )
   W( 46 ) = JVS( 393 )
   W( 47 ) = JVS( 394 )
   W( 48 ) = JVS( 395 )
   W( 49 ) = JVS( 396 )
   W( 50 ) = JVS( 397 )
   W( 51 ) = JVS( 398 )
   W( 52 ) = JVS( 399 )
   W( 53 ) = JVS( 400 )
   W( 54 ) = JVS( 401 )
   W( 55 ) = JVS( 402 )
   W( 56 ) = JVS( 403 )
   W( 57 ) = JVS( 404 )
  a = -W( 3 ) / JVS(            5  )
  W( 3 ) = -a
  W( 20 ) = W( 20 ) + a*JVS( 6 )
  W( 55 ) = W( 55 ) + a*JVS( 7 )
  a = -W( 9 ) / JVS(           24  )
  W( 9 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 25 )
  W( 55 ) = W( 55 ) + a*JVS( 26 )
  W( 57 ) = W( 57 ) + a*JVS( 27 )
  a = -W( 13 ) / JVS(           42  )
  W( 13 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 43 )
  W( 54 ) = W( 54 ) + a*JVS( 44 )
  W( 55 ) = W( 55 ) + a*JVS( 45 )
  a = -W( 20 ) / JVS(           72  )
  W( 20 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 73 )
  W( 55 ) = W( 55 ) + a*JVS( 74 )
  W( 57 ) = W( 57 ) + a*JVS( 75 )
  a = -W( 23 ) / JVS(           83  )
  W( 23 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 84 )
  W( 55 ) = W( 55 ) + a*JVS( 85 )
  W( 57 ) = W( 57 ) + a*JVS( 86 )
  a = -W( 28 ) / JVS(          109  )
  W( 28 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 110 )
  W( 50 ) = W( 50 ) + a*JVS( 111 )
  W( 55 ) = W( 55 ) + a*JVS( 112 )
  W( 56 ) = W( 56 ) + a*JVS( 113 )
  W( 57 ) = W( 57 ) + a*JVS( 114 )
  a = -W( 29 ) / JVS(          115  )
  W( 29 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 116 )
  W( 55 ) = W( 55 ) + a*JVS( 117 )
  W( 57 ) = W( 57 ) + a*JVS( 118 )
  a = -W( 30 ) / JVS(          120  )
  W( 30 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 121 )
  W( 44 ) = W( 44 ) + a*JVS( 122 )
  W( 50 ) = W( 50 ) + a*JVS( 123 )
  W( 55 ) = W( 55 ) + a*JVS( 124 )
  W( 56 ) = W( 56 ) + a*JVS( 125 )
  W( 57 ) = W( 57 ) + a*JVS( 126 )
  a = -W( 34 ) / JVS(          164  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 165 )
  W( 55 ) = W( 55 ) + a*JVS( 166 )
  W( 56 ) = W( 56 ) + a*JVS( 167 )
  W( 57 ) = W( 57 ) + a*JVS( 168 )
  a = -W( 35 ) / JVS(          171  )
  W( 35 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 172 )
  W( 37 ) = W( 37 ) + a*JVS( 173 )
  W( 39 ) = W( 39 ) + a*JVS( 174 )
  W( 41 ) = W( 41 ) + a*JVS( 175 )
  W( 48 ) = W( 48 ) + a*JVS( 176 )
  W( 49 ) = W( 49 ) + a*JVS( 177 )
  W( 50 ) = W( 50 ) + a*JVS( 178 )
  W( 53 ) = W( 53 ) + a*JVS( 179 )
  W( 55 ) = W( 55 ) + a*JVS( 180 )
  W( 56 ) = W( 56 ) + a*JVS( 181 )
  W( 57 ) = W( 57 ) + a*JVS( 182 )
  a = -W( 36 ) / JVS(          186  )
  W( 36 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 187 )
  W( 55 ) = W( 55 ) + a*JVS( 188 )
  W( 56 ) = W( 56 ) + a*JVS( 189 )
  W( 57 ) = W( 57 ) + a*JVS( 190 )
  a = -W( 37 ) / JVS(          191  )
  W( 37 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 192 )
  W( 53 ) = W( 53 ) + a*JVS( 193 )
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  a = -W( 38 ) / JVS(          195  )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 196 )
  W( 53 ) = W( 53 ) + a*JVS( 197 )
  W( 55 ) = W( 55 ) + a*JVS( 198 )
  a = -W( 39 ) / JVS(          205  )
  W( 39 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 206 )
  W( 48 ) = W( 48 ) + a*JVS( 207 )
  W( 49 ) = W( 49 ) + a*JVS( 208 )
  W( 50 ) = W( 50 ) + a*JVS( 209 )
  W( 52 ) = W( 52 ) + a*JVS( 210 )
  W( 53 ) = W( 53 ) + a*JVS( 211 )
  W( 54 ) = W( 54 ) + a*JVS( 212 )
  W( 55 ) = W( 55 ) + a*JVS( 213 )
  W( 56 ) = W( 56 ) + a*JVS( 214 )
  W( 57 ) = W( 57 ) + a*JVS( 215 )
  a = -W( 40 ) / JVS(          219  )
  W( 40 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 220 )
  W( 53 ) = W( 53 ) + a*JVS( 221 )
  W( 55 ) = W( 55 ) + a*JVS( 222 )
  W( 56 ) = W( 56 ) + a*JVS( 223 )
  W( 57 ) = W( 57 ) + a*JVS( 224 )
  a = -W( 41 ) / JVS(          228  )
  W( 41 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 229 )
  W( 50 ) = W( 50 ) + a*JVS( 230 )
  W( 52 ) = W( 52 ) + a*JVS( 231 )
  W( 53 ) = W( 53 ) + a*JVS( 232 )
  W( 55 ) = W( 55 ) + a*JVS( 233 )
  W( 56 ) = W( 56 ) + a*JVS( 234 )
  W( 57 ) = W( 57 ) + a*JVS( 235 )
  a = -W( 42 ) / JVS(          238  )
  W( 42 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 239 )
  W( 44 ) = W( 44 ) + a*JVS( 240 )
  W( 50 ) = W( 50 ) + a*JVS( 241 )
  W( 55 ) = W( 55 ) + a*JVS( 242 )
  W( 56 ) = W( 56 ) + a*JVS( 243 )
  W( 57 ) = W( 57 ) + a*JVS( 244 )
  a = -W( 43 ) / JVS(          247  )
  W( 43 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 248 )
  W( 55 ) = W( 55 ) + a*JVS( 249 )
  W( 56 ) = W( 56 ) + a*JVS( 250 )
  W( 57 ) = W( 57 ) + a*JVS( 251 )
  a = -W( 44 ) / JVS(          256  )
  W( 44 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 257 )
  W( 52 ) = W( 52 ) + a*JVS( 258 )
  W( 53 ) = W( 53 ) + a*JVS( 259 )
  W( 55 ) = W( 55 ) + a*JVS( 260 )
  W( 56 ) = W( 56 ) + a*JVS( 261 )
  W( 57 ) = W( 57 ) + a*JVS( 262 )
  a = -W( 46 ) / JVS(          304  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 305 )
  W( 50 ) = W( 50 ) + a*JVS( 306 )
  W( 51 ) = W( 51 ) + a*JVS( 307 )
  W( 52 ) = W( 52 ) + a*JVS( 308 )
  W( 53 ) = W( 53 ) + a*JVS( 309 )
  W( 55 ) = W( 55 ) + a*JVS( 310 )
  W( 56 ) = W( 56 ) + a*JVS( 311 )
  W( 57 ) = W( 57 ) + a*JVS( 312 )
  a = -W( 47 ) / JVS(          321  )
  W( 47 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 322 )
  W( 51 ) = W( 51 ) + a*JVS( 323 )
  W( 52 ) = W( 52 ) + a*JVS( 324 )
  W( 53 ) = W( 53 ) + a*JVS( 325 )
  W( 55 ) = W( 55 ) + a*JVS( 326 )
  W( 56 ) = W( 56 ) + a*JVS( 327 )
  W( 57 ) = W( 57 ) + a*JVS( 328 )
  a = -W( 48 ) / JVS(          330  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 331 )
  W( 52 ) = W( 52 ) + a*JVS( 332 )
  W( 53 ) = W( 53 ) + a*JVS( 333 )
  W( 55 ) = W( 55 ) + a*JVS( 334 )
  W( 56 ) = W( 56 ) + a*JVS( 335 )
  W( 57 ) = W( 57 ) + a*JVS( 336 )
  a = -W( 49 ) / JVS(          342  )
  W( 49 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 343 )
  W( 52 ) = W( 52 ) + a*JVS( 344 )
  W( 53 ) = W( 53 ) + a*JVS( 345 )
  W( 54 ) = W( 54 ) + a*JVS( 346 )
  W( 55 ) = W( 55 ) + a*JVS( 347 )
  W( 56 ) = W( 56 ) + a*JVS( 348 )
  W( 57 ) = W( 57 ) + a*JVS( 349 )
  a = -W( 50 ) / JVS(          366  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 367 )
  W( 52 ) = W( 52 ) + a*JVS( 368 )
  W( 53 ) = W( 53 ) + a*JVS( 369 )
  W( 54 ) = W( 54 ) + a*JVS( 370 )
  W( 55 ) = W( 55 ) + a*JVS( 371 )
  W( 56 ) = W( 56 ) + a*JVS( 372 )
  W( 57 ) = W( 57 ) + a*JVS( 373 )
  JVS( 374) = W( 3 )
  JVS( 375) = W( 9 )
  JVS( 376) = W( 13 )
  JVS( 377) = W( 20 )
  JVS( 378) = W( 23 )
  JVS( 379) = W( 28 )
  JVS( 380) = W( 29 )
  JVS( 381) = W( 30 )
  JVS( 382) = W( 34 )
  JVS( 383) = W( 35 )
  JVS( 384) = W( 36 )
  JVS( 385) = W( 37 )
  JVS( 386) = W( 38 )
  JVS( 387) = W( 39 )
  JVS( 388) = W( 40 )
  JVS( 389) = W( 41 )
  JVS( 390) = W( 42 )
  JVS( 391) = W( 43 )
  JVS( 392) = W( 44 )
  JVS( 393) = W( 46 )
  JVS( 394) = W( 47 )
  JVS( 395) = W( 48 )
  JVS( 396) = W( 49 )
  JVS( 397) = W( 50 )
  JVS( 398) = W( 51 )
  JVS( 399) = W( 52 )
  JVS( 400) = W( 53 )
  JVS( 401) = W( 54 )
  JVS( 402) = W( 55 )
  JVS( 403) = W( 56 )
  JVS( 404) = W( 57 )
  IF ( ABS(  JVS( 418 )) < TINY(a) ) THEN
         IER = 52                                      
         RETURN
  END IF
   W( 5 ) = JVS( 405 )
   W( 6 ) = JVS( 406 )
   W( 22 ) = JVS( 407 )
   W( 25 ) = JVS( 408 )
   W( 37 ) = JVS( 409 )
   W( 38 ) = JVS( 410 )
   W( 45 ) = JVS( 411 )
   W( 46 ) = JVS( 412 )
   W( 47 ) = JVS( 413 )
   W( 48 ) = JVS( 414 )
   W( 49 ) = JVS( 415 )
   W( 50 ) = JVS( 416 )
   W( 51 ) = JVS( 417 )
   W( 52 ) = JVS( 418 )
   W( 53 ) = JVS( 419 )
   W( 54 ) = JVS( 420 )
   W( 55 ) = JVS( 421 )
   W( 56 ) = JVS( 422 )
   W( 57 ) = JVS( 423 )
  a = -W( 5 ) / JVS(           10  )
  W( 5 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 11 )
  W( 54 ) = W( 54 ) + a*JVS( 12 )
  a = -W( 6 ) / JVS(           13  )
  W( 6 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 14 )
  W( 54 ) = W( 54 ) + a*JVS( 15 )
  a = -W( 22 ) / JVS(           79  )
  W( 22 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 80 )
  W( 55 ) = W( 55 ) + a*JVS( 81 )
  W( 57 ) = W( 57 ) + a*JVS( 82 )
  a = -W( 25 ) / JVS(           93  )
  W( 25 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 94 )
  W( 46 ) = W( 46 ) + a*JVS( 95 )
  W( 52 ) = W( 52 ) + a*JVS( 96 )
  W( 54 ) = W( 54 ) + a*JVS( 97 )
  W( 55 ) = W( 55 ) + a*JVS( 98 )
  a = -W( 37 ) / JVS(          191  )
  W( 37 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 192 )
  W( 53 ) = W( 53 ) + a*JVS( 193 )
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  a = -W( 38 ) / JVS(          195  )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 196 )
  W( 53 ) = W( 53 ) + a*JVS( 197 )
  W( 55 ) = W( 55 ) + a*JVS( 198 )
  a = -W( 45 ) / JVS(          283  )
  W( 45 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 284 )
  W( 48 ) = W( 48 ) + a*JVS( 285 )
  W( 49 ) = W( 49 ) + a*JVS( 286 )
  W( 50 ) = W( 50 ) + a*JVS( 287 )
  W( 51 ) = W( 51 ) + a*JVS( 288 )
  W( 52 ) = W( 52 ) + a*JVS( 289 )
  W( 53 ) = W( 53 ) + a*JVS( 290 )
  W( 54 ) = W( 54 ) + a*JVS( 291 )
  W( 55 ) = W( 55 ) + a*JVS( 292 )
  W( 56 ) = W( 56 ) + a*JVS( 293 )
  W( 57 ) = W( 57 ) + a*JVS( 294 )
  a = -W( 46 ) / JVS(          304  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 305 )
  W( 50 ) = W( 50 ) + a*JVS( 306 )
  W( 51 ) = W( 51 ) + a*JVS( 307 )
  W( 52 ) = W( 52 ) + a*JVS( 308 )
  W( 53 ) = W( 53 ) + a*JVS( 309 )
  W( 55 ) = W( 55 ) + a*JVS( 310 )
  W( 56 ) = W( 56 ) + a*JVS( 311 )
  W( 57 ) = W( 57 ) + a*JVS( 312 )
  a = -W( 47 ) / JVS(          321  )
  W( 47 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 322 )
  W( 51 ) = W( 51 ) + a*JVS( 323 )
  W( 52 ) = W( 52 ) + a*JVS( 324 )
  W( 53 ) = W( 53 ) + a*JVS( 325 )
  W( 55 ) = W( 55 ) + a*JVS( 326 )
  W( 56 ) = W( 56 ) + a*JVS( 327 )
  W( 57 ) = W( 57 ) + a*JVS( 328 )
  a = -W( 48 ) / JVS(          330  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 331 )
  W( 52 ) = W( 52 ) + a*JVS( 332 )
  W( 53 ) = W( 53 ) + a*JVS( 333 )
  W( 55 ) = W( 55 ) + a*JVS( 334 )
  W( 56 ) = W( 56 ) + a*JVS( 335 )
  W( 57 ) = W( 57 ) + a*JVS( 336 )
  a = -W( 49 ) / JVS(          342  )
  W( 49 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 343 )
  W( 52 ) = W( 52 ) + a*JVS( 344 )
  W( 53 ) = W( 53 ) + a*JVS( 345 )
  W( 54 ) = W( 54 ) + a*JVS( 346 )
  W( 55 ) = W( 55 ) + a*JVS( 347 )
  W( 56 ) = W( 56 ) + a*JVS( 348 )
  W( 57 ) = W( 57 ) + a*JVS( 349 )
  a = -W( 50 ) / JVS(          366  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 367 )
  W( 52 ) = W( 52 ) + a*JVS( 368 )
  W( 53 ) = W( 53 ) + a*JVS( 369 )
  W( 54 ) = W( 54 ) + a*JVS( 370 )
  W( 55 ) = W( 55 ) + a*JVS( 371 )
  W( 56 ) = W( 56 ) + a*JVS( 372 )
  W( 57 ) = W( 57 ) + a*JVS( 373 )
  a = -W( 51 ) / JVS(          398  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 399 )
  W( 53 ) = W( 53 ) + a*JVS( 400 )
  W( 54 ) = W( 54 ) + a*JVS( 401 )
  W( 55 ) = W( 55 ) + a*JVS( 402 )
  W( 56 ) = W( 56 ) + a*JVS( 403 )
  W( 57 ) = W( 57 ) + a*JVS( 404 )
  JVS( 405) = W( 5 )
  JVS( 406) = W( 6 )
  JVS( 407) = W( 22 )
  JVS( 408) = W( 25 )
  JVS( 409) = W( 37 )
  JVS( 410) = W( 38 )
  JVS( 411) = W( 45 )
  JVS( 412) = W( 46 )
  JVS( 413) = W( 47 )
  JVS( 414) = W( 48 )
  JVS( 415) = W( 49 )
  JVS( 416) = W( 50 )
  JVS( 417) = W( 51 )
  JVS( 418) = W( 52 )
  JVS( 419) = W( 53 )
  JVS( 420) = W( 54 )
  JVS( 421) = W( 55 )
  JVS( 422) = W( 56 )
  JVS( 423) = W( 57 )
  IF ( ABS(  JVS( 434 )) < TINY(a) ) THEN
         IER = 53                                      
         RETURN
  END IF
   W( 19 ) = JVS( 424 )
   W( 21 ) = JVS( 425 )
   W( 24 ) = JVS( 426 )
   W( 37 ) = JVS( 427 )
   W( 38 ) = JVS( 428 )
   W( 41 ) = JVS( 429 )
   W( 48 ) = JVS( 430 )
   W( 50 ) = JVS( 431 )
   W( 51 ) = JVS( 432 )
   W( 52 ) = JVS( 433 )
   W( 53 ) = JVS( 434 )
   W( 54 ) = JVS( 435 )
   W( 55 ) = JVS( 436 )
   W( 56 ) = JVS( 437 )
   W( 57 ) = JVS( 438 )
  a = -W( 19 ) / JVS(           69  )
  W( 19 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 70 )
  W( 55 ) = W( 55 ) + a*JVS( 71 )
  a = -W( 21 ) / JVS(           76  )
  W( 21 ) = -a
  W( 24 ) = W( 24 ) + a*JVS( 77 )
  W( 53 ) = W( 53 ) + a*JVS( 78 )
  a = -W( 24 ) / JVS(           88  )
  W( 24 ) = -a
  W( 37 ) = W( 37 ) + a*JVS( 89 )
  W( 53 ) = W( 53 ) + a*JVS( 90 )
  W( 55 ) = W( 55 ) + a*JVS( 91 )
  a = -W( 37 ) / JVS(          191  )
  W( 37 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 192 )
  W( 53 ) = W( 53 ) + a*JVS( 193 )
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  a = -W( 38 ) / JVS(          195  )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 196 )
  W( 53 ) = W( 53 ) + a*JVS( 197 )
  W( 55 ) = W( 55 ) + a*JVS( 198 )
  a = -W( 41 ) / JVS(          228  )
  W( 41 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 229 )
  W( 50 ) = W( 50 ) + a*JVS( 230 )
  W( 52 ) = W( 52 ) + a*JVS( 231 )
  W( 53 ) = W( 53 ) + a*JVS( 232 )
  W( 55 ) = W( 55 ) + a*JVS( 233 )
  W( 56 ) = W( 56 ) + a*JVS( 234 )
  W( 57 ) = W( 57 ) + a*JVS( 235 )
  a = -W( 48 ) / JVS(          330  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 331 )
  W( 52 ) = W( 52 ) + a*JVS( 332 )
  W( 53 ) = W( 53 ) + a*JVS( 333 )
  W( 55 ) = W( 55 ) + a*JVS( 334 )
  W( 56 ) = W( 56 ) + a*JVS( 335 )
  W( 57 ) = W( 57 ) + a*JVS( 336 )
  a = -W( 50 ) / JVS(          366  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 367 )
  W( 52 ) = W( 52 ) + a*JVS( 368 )
  W( 53 ) = W( 53 ) + a*JVS( 369 )
  W( 54 ) = W( 54 ) + a*JVS( 370 )
  W( 55 ) = W( 55 ) + a*JVS( 371 )
  W( 56 ) = W( 56 ) + a*JVS( 372 )
  W( 57 ) = W( 57 ) + a*JVS( 373 )
  a = -W( 51 ) / JVS(          398  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 399 )
  W( 53 ) = W( 53 ) + a*JVS( 400 )
  W( 54 ) = W( 54 ) + a*JVS( 401 )
  W( 55 ) = W( 55 ) + a*JVS( 402 )
  W( 56 ) = W( 56 ) + a*JVS( 403 )
  W( 57 ) = W( 57 ) + a*JVS( 404 )
  a = -W( 52 ) / JVS(          418  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 419 )
  W( 54 ) = W( 54 ) + a*JVS( 420 )
  W( 55 ) = W( 55 ) + a*JVS( 421 )
  W( 56 ) = W( 56 ) + a*JVS( 422 )
  W( 57 ) = W( 57 ) + a*JVS( 423 )
  JVS( 424) = W( 19 )
  JVS( 425) = W( 21 )
  JVS( 426) = W( 24 )
  JVS( 427) = W( 37 )
  JVS( 428) = W( 38 )
  JVS( 429) = W( 41 )
  JVS( 430) = W( 48 )
  JVS( 431) = W( 50 )
  JVS( 432) = W( 51 )
  JVS( 433) = W( 52 )
  JVS( 434) = W( 53 )
  JVS( 435) = W( 54 )
  JVS( 436) = W( 55 )
  JVS( 437) = W( 56 )
  JVS( 438) = W( 57 )
  IF ( ABS(  JVS( 464 )) < TINY(a) ) THEN
         IER = 54                                      
         RETURN
  END IF
   W( 5 ) = JVS( 439 )
   W( 6 ) = JVS( 440 )
   W( 7 ) = JVS( 441 )
   W( 13 ) = JVS( 442 )
   W( 14 ) = JVS( 443 )
   W( 18 ) = JVS( 444 )
   W( 22 ) = JVS( 445 )
   W( 25 ) = JVS( 446 )
   W( 33 ) = JVS( 447 )
   W( 34 ) = JVS( 448 )
   W( 36 ) = JVS( 449 )
   W( 38 ) = JVS( 450 )
   W( 40 ) = JVS( 451 )
   W( 42 ) = JVS( 452 )
   W( 43 ) = JVS( 453 )
   W( 44 ) = JVS( 454 )
   W( 45 ) = JVS( 455 )
   W( 46 ) = JVS( 456 )
   W( 47 ) = JVS( 457 )
   W( 48 ) = JVS( 458 )
   W( 49 ) = JVS( 459 )
   W( 50 ) = JVS( 460 )
   W( 51 ) = JVS( 461 )
   W( 52 ) = JVS( 462 )
   W( 53 ) = JVS( 463 )
   W( 54 ) = JVS( 464 )
   W( 55 ) = JVS( 465 )
   W( 56 ) = JVS( 466 )
   W( 57 ) = JVS( 467 )
  a = -W( 5 ) / JVS(           10  )
  W( 5 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 11 )
  W( 54 ) = W( 54 ) + a*JVS( 12 )
  a = -W( 6 ) / JVS(           13  )
  W( 6 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 14 )
  W( 54 ) = W( 54 ) + a*JVS( 15 )
  a = -W( 7 ) / JVS(           16  )
  W( 7 ) = -a
  W( 33 ) = W( 33 ) + a*JVS( 17 )
  W( 55 ) = W( 55 ) + a*JVS( 18 )
  a = -W( 13 ) / JVS(           42  )
  W( 13 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 43 )
  W( 54 ) = W( 54 ) + a*JVS( 44 )
  W( 55 ) = W( 55 ) + a*JVS( 45 )
  a = -W( 14 ) / JVS(           46  )
  W( 14 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 47 )
  W( 54 ) = W( 54 ) + a*JVS( 48 )
  W( 55 ) = W( 55 ) + a*JVS( 49 )
  a = -W( 18 ) / JVS(           65  )
  W( 18 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 66 )
  W( 55 ) = W( 55 ) + a*JVS( 67 )
  W( 56 ) = W( 56 ) + a*JVS( 68 )
  a = -W( 22 ) / JVS(           79  )
  W( 22 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 80 )
  W( 55 ) = W( 55 ) + a*JVS( 81 )
  W( 57 ) = W( 57 ) + a*JVS( 82 )
  a = -W( 25 ) / JVS(           93  )
  W( 25 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 94 )
  W( 46 ) = W( 46 ) + a*JVS( 95 )
  W( 52 ) = W( 52 ) + a*JVS( 96 )
  W( 54 ) = W( 54 ) + a*JVS( 97 )
  W( 55 ) = W( 55 ) + a*JVS( 98 )
  a = -W( 33 ) / JVS(          156  )
  W( 33 ) = -a
  W( 38 ) = W( 38 ) + a*JVS( 157 )
  W( 48 ) = W( 48 ) + a*JVS( 158 )
  W( 52 ) = W( 52 ) + a*JVS( 159 )
  W( 55 ) = W( 55 ) + a*JVS( 160 )
  W( 56 ) = W( 56 ) + a*JVS( 161 )
  a = -W( 34 ) / JVS(          164  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 165 )
  W( 55 ) = W( 55 ) + a*JVS( 166 )
  W( 56 ) = W( 56 ) + a*JVS( 167 )
  W( 57 ) = W( 57 ) + a*JVS( 168 )
  a = -W( 36 ) / JVS(          186  )
  W( 36 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 187 )
  W( 55 ) = W( 55 ) + a*JVS( 188 )
  W( 56 ) = W( 56 ) + a*JVS( 189 )
  W( 57 ) = W( 57 ) + a*JVS( 190 )
  a = -W( 38 ) / JVS(          195  )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 196 )
  W( 53 ) = W( 53 ) + a*JVS( 197 )
  W( 55 ) = W( 55 ) + a*JVS( 198 )
  a = -W( 40 ) / JVS(          219  )
  W( 40 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 220 )
  W( 53 ) = W( 53 ) + a*JVS( 221 )
  W( 55 ) = W( 55 ) + a*JVS( 222 )
  W( 56 ) = W( 56 ) + a*JVS( 223 )
  W( 57 ) = W( 57 ) + a*JVS( 224 )
  a = -W( 42 ) / JVS(          238  )
  W( 42 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 239 )
  W( 44 ) = W( 44 ) + a*JVS( 240 )
  W( 50 ) = W( 50 ) + a*JVS( 241 )
  W( 55 ) = W( 55 ) + a*JVS( 242 )
  W( 56 ) = W( 56 ) + a*JVS( 243 )
  W( 57 ) = W( 57 ) + a*JVS( 244 )
  a = -W( 43 ) / JVS(          247  )
  W( 43 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 248 )
  W( 55 ) = W( 55 ) + a*JVS( 249 )
  W( 56 ) = W( 56 ) + a*JVS( 250 )
  W( 57 ) = W( 57 ) + a*JVS( 251 )
  a = -W( 44 ) / JVS(          256  )
  W( 44 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 257 )
  W( 52 ) = W( 52 ) + a*JVS( 258 )
  W( 53 ) = W( 53 ) + a*JVS( 259 )
  W( 55 ) = W( 55 ) + a*JVS( 260 )
  W( 56 ) = W( 56 ) + a*JVS( 261 )
  W( 57 ) = W( 57 ) + a*JVS( 262 )
  a = -W( 45 ) / JVS(          283  )
  W( 45 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 284 )
  W( 48 ) = W( 48 ) + a*JVS( 285 )
  W( 49 ) = W( 49 ) + a*JVS( 286 )
  W( 50 ) = W( 50 ) + a*JVS( 287 )
  W( 51 ) = W( 51 ) + a*JVS( 288 )
  W( 52 ) = W( 52 ) + a*JVS( 289 )
  W( 53 ) = W( 53 ) + a*JVS( 290 )
  W( 54 ) = W( 54 ) + a*JVS( 291 )
  W( 55 ) = W( 55 ) + a*JVS( 292 )
  W( 56 ) = W( 56 ) + a*JVS( 293 )
  W( 57 ) = W( 57 ) + a*JVS( 294 )
  a = -W( 46 ) / JVS(          304  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 305 )
  W( 50 ) = W( 50 ) + a*JVS( 306 )
  W( 51 ) = W( 51 ) + a*JVS( 307 )
  W( 52 ) = W( 52 ) + a*JVS( 308 )
  W( 53 ) = W( 53 ) + a*JVS( 309 )
  W( 55 ) = W( 55 ) + a*JVS( 310 )
  W( 56 ) = W( 56 ) + a*JVS( 311 )
  W( 57 ) = W( 57 ) + a*JVS( 312 )
  a = -W( 47 ) / JVS(          321  )
  W( 47 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 322 )
  W( 51 ) = W( 51 ) + a*JVS( 323 )
  W( 52 ) = W( 52 ) + a*JVS( 324 )
  W( 53 ) = W( 53 ) + a*JVS( 325 )
  W( 55 ) = W( 55 ) + a*JVS( 326 )
  W( 56 ) = W( 56 ) + a*JVS( 327 )
  W( 57 ) = W( 57 ) + a*JVS( 328 )
  a = -W( 48 ) / JVS(          330  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 331 )
  W( 52 ) = W( 52 ) + a*JVS( 332 )
  W( 53 ) = W( 53 ) + a*JVS( 333 )
  W( 55 ) = W( 55 ) + a*JVS( 334 )
  W( 56 ) = W( 56 ) + a*JVS( 335 )
  W( 57 ) = W( 57 ) + a*JVS( 336 )
  a = -W( 49 ) / JVS(          342  )
  W( 49 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 343 )
  W( 52 ) = W( 52 ) + a*JVS( 344 )
  W( 53 ) = W( 53 ) + a*JVS( 345 )
  W( 54 ) = W( 54 ) + a*JVS( 346 )
  W( 55 ) = W( 55 ) + a*JVS( 347 )
  W( 56 ) = W( 56 ) + a*JVS( 348 )
  W( 57 ) = W( 57 ) + a*JVS( 349 )
  a = -W( 50 ) / JVS(          366  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 367 )
  W( 52 ) = W( 52 ) + a*JVS( 368 )
  W( 53 ) = W( 53 ) + a*JVS( 369 )
  W( 54 ) = W( 54 ) + a*JVS( 370 )
  W( 55 ) = W( 55 ) + a*JVS( 371 )
  W( 56 ) = W( 56 ) + a*JVS( 372 )
  W( 57 ) = W( 57 ) + a*JVS( 373 )
  a = -W( 51 ) / JVS(          398  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 399 )
  W( 53 ) = W( 53 ) + a*JVS( 400 )
  W( 54 ) = W( 54 ) + a*JVS( 401 )
  W( 55 ) = W( 55 ) + a*JVS( 402 )
  W( 56 ) = W( 56 ) + a*JVS( 403 )
  W( 57 ) = W( 57 ) + a*JVS( 404 )
  a = -W( 52 ) / JVS(          418  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 419 )
  W( 54 ) = W( 54 ) + a*JVS( 420 )
  W( 55 ) = W( 55 ) + a*JVS( 421 )
  W( 56 ) = W( 56 ) + a*JVS( 422 )
  W( 57 ) = W( 57 ) + a*JVS( 423 )
  a = -W( 53 ) / JVS(          434  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 435 )
  W( 55 ) = W( 55 ) + a*JVS( 436 )
  W( 56 ) = W( 56 ) + a*JVS( 437 )
  W( 57 ) = W( 57 ) + a*JVS( 438 )
  JVS( 439) = W( 5 )
  JVS( 440) = W( 6 )
  JVS( 441) = W( 7 )
  JVS( 442) = W( 13 )
  JVS( 443) = W( 14 )
  JVS( 444) = W( 18 )
  JVS( 445) = W( 22 )
  JVS( 446) = W( 25 )
  JVS( 447) = W( 33 )
  JVS( 448) = W( 34 )
  JVS( 449) = W( 36 )
  JVS( 450) = W( 38 )
  JVS( 451) = W( 40 )
  JVS( 452) = W( 42 )
  JVS( 453) = W( 43 )
  JVS( 454) = W( 44 )
  JVS( 455) = W( 45 )
  JVS( 456) = W( 46 )
  JVS( 457) = W( 47 )
  JVS( 458) = W( 48 )
  JVS( 459) = W( 49 )
  JVS( 460) = W( 50 )
  JVS( 461) = W( 51 )
  JVS( 462) = W( 52 )
  JVS( 463) = W( 53 )
  JVS( 464) = W( 54 )
  JVS( 465) = W( 55 )
  JVS( 466) = W( 56 )
  JVS( 467) = W( 57 )
  IF ( ABS(  JVS( 519 )) < TINY(a) ) THEN
         IER = 55                                      
         RETURN
  END IF
   W( 1 ) = JVS( 468 )
   W( 2 ) = JVS( 469 )
   W( 4 ) = JVS( 470 )
   W( 7 ) = JVS( 471 )
   W( 8 ) = JVS( 472 )
   W( 9 ) = JVS( 473 )
   W( 10 ) = JVS( 474 )
   W( 11 ) = JVS( 475 )
   W( 12 ) = JVS( 476 )
   W( 13 ) = JVS( 477 )
   W( 14 ) = JVS( 478 )
   W( 15 ) = JVS( 479 )
   W( 16 ) = JVS( 480 )
   W( 17 ) = JVS( 481 )
   W( 18 ) = JVS( 482 )
   W( 19 ) = JVS( 483 )
   W( 20 ) = JVS( 484 )
   W( 21 ) = JVS( 485 )
   W( 22 ) = JVS( 486 )
   W( 23 ) = JVS( 487 )
   W( 24 ) = JVS( 488 )
   W( 25 ) = JVS( 489 )
   W( 26 ) = JVS( 490 )
   W( 27 ) = JVS( 491 )
   W( 28 ) = JVS( 492 )
   W( 29 ) = JVS( 493 )
   W( 30 ) = JVS( 494 )
   W( 31 ) = JVS( 495 )
   W( 32 ) = JVS( 496 )
   W( 33 ) = JVS( 497 )
   W( 34 ) = JVS( 498 )
   W( 35 ) = JVS( 499 )
   W( 36 ) = JVS( 500 )
   W( 37 ) = JVS( 501 )
   W( 38 ) = JVS( 502 )
   W( 39 ) = JVS( 503 )
   W( 40 ) = JVS( 504 )
   W( 41 ) = JVS( 505 )
   W( 42 ) = JVS( 506 )
   W( 43 ) = JVS( 507 )
   W( 44 ) = JVS( 508 )
   W( 45 ) = JVS( 509 )
   W( 46 ) = JVS( 510 )
   W( 47 ) = JVS( 511 )
   W( 48 ) = JVS( 512 )
   W( 49 ) = JVS( 513 )
   W( 50 ) = JVS( 514 )
   W( 51 ) = JVS( 515 )
   W( 52 ) = JVS( 516 )
   W( 53 ) = JVS( 517 )
   W( 54 ) = JVS( 518 )
   W( 55 ) = JVS( 519 )
   W( 56 ) = JVS( 520 )
   W( 57 ) = JVS( 521 )
  a = -W( 1 ) / JVS(            1  )
  W( 1 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 2 )
  a = -W( 2 ) / JVS(            3  )
  W( 2 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 4 )
  a = -W( 4 ) / JVS(            8  )
  W( 4 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 9 )
  a = -W( 7 ) / JVS(           16  )
  W( 7 ) = -a
  W( 33 ) = W( 33 ) + a*JVS( 17 )
  W( 55 ) = W( 55 ) + a*JVS( 18 )
  a = -W( 8 ) / JVS(           19  )
  W( 8 ) = -a
  W( 38 ) = W( 38 ) + a*JVS( 20 )
  W( 53 ) = W( 53 ) + a*JVS( 21 )
  W( 55 ) = W( 55 ) + a*JVS( 22 )
  W( 57 ) = W( 57 ) + a*JVS( 23 )
  a = -W( 9 ) / JVS(           24  )
  W( 9 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 25 )
  W( 55 ) = W( 55 ) + a*JVS( 26 )
  W( 57 ) = W( 57 ) + a*JVS( 27 )
  a = -W( 10 ) / JVS(           28  )
  W( 10 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 29 )
  W( 50 ) = W( 50 ) + a*JVS( 30 )
  W( 51 ) = W( 51 ) + a*JVS( 31 )
  W( 55 ) = W( 55 ) + a*JVS( 32 )
  W( 57 ) = W( 57 ) + a*JVS( 33 )
  a = -W( 11 ) / JVS(           34  )
  W( 11 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 35 )
  W( 55 ) = W( 55 ) + a*JVS( 36 )
  W( 57 ) = W( 57 ) + a*JVS( 37 )
  a = -W( 12 ) / JVS(           38  )
  W( 12 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 39 )
  W( 55 ) = W( 55 ) + a*JVS( 40 )
  W( 57 ) = W( 57 ) + a*JVS( 41 )
  a = -W( 13 ) / JVS(           42  )
  W( 13 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 43 )
  W( 54 ) = W( 54 ) + a*JVS( 44 )
  W( 55 ) = W( 55 ) + a*JVS( 45 )
  a = -W( 14 ) / JVS(           46  )
  W( 14 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 47 )
  W( 54 ) = W( 54 ) + a*JVS( 48 )
  W( 55 ) = W( 55 ) + a*JVS( 49 )
  a = -W( 15 ) / JVS(           50  )
  W( 15 ) = -a
  W( 19 ) = W( 19 ) + a*JVS( 51 )
  W( 37 ) = W( 37 ) + a*JVS( 52 )
  W( 38 ) = W( 38 ) + a*JVS( 53 )
  W( 41 ) = W( 41 ) + a*JVS( 54 )
  W( 53 ) = W( 53 ) + a*JVS( 55 )
  W( 55 ) = W( 55 ) + a*JVS( 56 )
  a = -W( 16 ) / JVS(           57  )
  W( 16 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 58 )
  W( 55 ) = W( 55 ) + a*JVS( 59 )
  W( 57 ) = W( 57 ) + a*JVS( 60 )
  a = -W( 17 ) / JVS(           61  )
  W( 17 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 62 )
  W( 55 ) = W( 55 ) + a*JVS( 63 )
  W( 57 ) = W( 57 ) + a*JVS( 64 )
  a = -W( 18 ) / JVS(           65  )
  W( 18 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 66 )
  W( 55 ) = W( 55 ) + a*JVS( 67 )
  W( 56 ) = W( 56 ) + a*JVS( 68 )
  a = -W( 19 ) / JVS(           69  )
  W( 19 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 70 )
  W( 55 ) = W( 55 ) + a*JVS( 71 )
  a = -W( 20 ) / JVS(           72  )
  W( 20 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 73 )
  W( 55 ) = W( 55 ) + a*JVS( 74 )
  W( 57 ) = W( 57 ) + a*JVS( 75 )
  a = -W( 21 ) / JVS(           76  )
  W( 21 ) = -a
  W( 24 ) = W( 24 ) + a*JVS( 77 )
  W( 53 ) = W( 53 ) + a*JVS( 78 )
  a = -W( 22 ) / JVS(           79  )
  W( 22 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 80 )
  W( 55 ) = W( 55 ) + a*JVS( 81 )
  W( 57 ) = W( 57 ) + a*JVS( 82 )
  a = -W( 23 ) / JVS(           83  )
  W( 23 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 84 )
  W( 55 ) = W( 55 ) + a*JVS( 85 )
  W( 57 ) = W( 57 ) + a*JVS( 86 )
  a = -W( 24 ) / JVS(           88  )
  W( 24 ) = -a
  W( 37 ) = W( 37 ) + a*JVS( 89 )
  W( 53 ) = W( 53 ) + a*JVS( 90 )
  W( 55 ) = W( 55 ) + a*JVS( 91 )
  a = -W( 25 ) / JVS(           93  )
  W( 25 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 94 )
  W( 46 ) = W( 46 ) + a*JVS( 95 )
  W( 52 ) = W( 52 ) + a*JVS( 96 )
  W( 54 ) = W( 54 ) + a*JVS( 97 )
  W( 55 ) = W( 55 ) + a*JVS( 98 )
  a = -W( 26 ) / JVS(           99  )
  W( 26 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 100 )
  W( 55 ) = W( 55 ) + a*JVS( 101 )
  W( 57 ) = W( 57 ) + a*JVS( 102 )
  a = -W( 27 ) / JVS(          103  )
  W( 27 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 104 )
  W( 55 ) = W( 55 ) + a*JVS( 105 )
  W( 57 ) = W( 57 ) + a*JVS( 106 )
  a = -W( 28 ) / JVS(          109  )
  W( 28 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 110 )
  W( 50 ) = W( 50 ) + a*JVS( 111 )
  W( 55 ) = W( 55 ) + a*JVS( 112 )
  W( 56 ) = W( 56 ) + a*JVS( 113 )
  W( 57 ) = W( 57 ) + a*JVS( 114 )
  a = -W( 29 ) / JVS(          115  )
  W( 29 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 116 )
  W( 55 ) = W( 55 ) + a*JVS( 117 )
  W( 57 ) = W( 57 ) + a*JVS( 118 )
  a = -W( 30 ) / JVS(          120  )
  W( 30 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 121 )
  W( 44 ) = W( 44 ) + a*JVS( 122 )
  W( 50 ) = W( 50 ) + a*JVS( 123 )
  W( 55 ) = W( 55 ) + a*JVS( 124 )
  W( 56 ) = W( 56 ) + a*JVS( 125 )
  W( 57 ) = W( 57 ) + a*JVS( 126 )
  a = -W( 31 ) / JVS(          130  )
  W( 31 ) = -a
  W( 33 ) = W( 33 ) + a*JVS( 131 )
  W( 35 ) = W( 35 ) + a*JVS( 132 )
  W( 37 ) = W( 37 ) + a*JVS( 133 )
  W( 38 ) = W( 38 ) + a*JVS( 134 )
  W( 41 ) = W( 41 ) + a*JVS( 135 )
  W( 45 ) = W( 45 ) + a*JVS( 136 )
  W( 46 ) = W( 46 ) + a*JVS( 137 )
  W( 49 ) = W( 49 ) + a*JVS( 138 )
  W( 50 ) = W( 50 ) + a*JVS( 139 )
  W( 52 ) = W( 52 ) + a*JVS( 140 )
  W( 53 ) = W( 53 ) + a*JVS( 141 )
  W( 55 ) = W( 55 ) + a*JVS( 142 )
  W( 56 ) = W( 56 ) + a*JVS( 143 )
  W( 57 ) = W( 57 ) + a*JVS( 144 )
  a = -W( 32 ) / JVS(          145  )
  W( 32 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 146 )
  W( 36 ) = W( 36 ) + a*JVS( 147 )
  W( 37 ) = W( 37 ) + a*JVS( 148 )
  W( 43 ) = W( 43 ) + a*JVS( 149 )
  W( 47 ) = W( 47 ) + a*JVS( 150 )
  W( 48 ) = W( 48 ) + a*JVS( 151 )
  W( 49 ) = W( 49 ) + a*JVS( 152 )
  W( 50 ) = W( 50 ) + a*JVS( 153 )
  W( 53 ) = W( 53 ) + a*JVS( 154 )
  W( 55 ) = W( 55 ) + a*JVS( 155 )
  a = -W( 33 ) / JVS(          156  )
  W( 33 ) = -a
  W( 38 ) = W( 38 ) + a*JVS( 157 )
  W( 48 ) = W( 48 ) + a*JVS( 158 )
  W( 52 ) = W( 52 ) + a*JVS( 159 )
  W( 55 ) = W( 55 ) + a*JVS( 160 )
  W( 56 ) = W( 56 ) + a*JVS( 161 )
  a = -W( 34 ) / JVS(          164  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 165 )
  W( 55 ) = W( 55 ) + a*JVS( 166 )
  W( 56 ) = W( 56 ) + a*JVS( 167 )
  W( 57 ) = W( 57 ) + a*JVS( 168 )
  a = -W( 35 ) / JVS(          171  )
  W( 35 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 172 )
  W( 37 ) = W( 37 ) + a*JVS( 173 )
  W( 39 ) = W( 39 ) + a*JVS( 174 )
  W( 41 ) = W( 41 ) + a*JVS( 175 )
  W( 48 ) = W( 48 ) + a*JVS( 176 )
  W( 49 ) = W( 49 ) + a*JVS( 177 )
  W( 50 ) = W( 50 ) + a*JVS( 178 )
  W( 53 ) = W( 53 ) + a*JVS( 179 )
  W( 55 ) = W( 55 ) + a*JVS( 180 )
  W( 56 ) = W( 56 ) + a*JVS( 181 )
  W( 57 ) = W( 57 ) + a*JVS( 182 )
  a = -W( 36 ) / JVS(          186  )
  W( 36 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 187 )
  W( 55 ) = W( 55 ) + a*JVS( 188 )
  W( 56 ) = W( 56 ) + a*JVS( 189 )
  W( 57 ) = W( 57 ) + a*JVS( 190 )
  a = -W( 37 ) / JVS(          191  )
  W( 37 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 192 )
  W( 53 ) = W( 53 ) + a*JVS( 193 )
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  a = -W( 38 ) / JVS(          195  )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 196 )
  W( 53 ) = W( 53 ) + a*JVS( 197 )
  W( 55 ) = W( 55 ) + a*JVS( 198 )
  a = -W( 39 ) / JVS(          205  )
  W( 39 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 206 )
  W( 48 ) = W( 48 ) + a*JVS( 207 )
  W( 49 ) = W( 49 ) + a*JVS( 208 )
  W( 50 ) = W( 50 ) + a*JVS( 209 )
  W( 52 ) = W( 52 ) + a*JVS( 210 )
  W( 53 ) = W( 53 ) + a*JVS( 211 )
  W( 54 ) = W( 54 ) + a*JVS( 212 )
  W( 55 ) = W( 55 ) + a*JVS( 213 )
  W( 56 ) = W( 56 ) + a*JVS( 214 )
  W( 57 ) = W( 57 ) + a*JVS( 215 )
  a = -W( 40 ) / JVS(          219  )
  W( 40 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 220 )
  W( 53 ) = W( 53 ) + a*JVS( 221 )
  W( 55 ) = W( 55 ) + a*JVS( 222 )
  W( 56 ) = W( 56 ) + a*JVS( 223 )
  W( 57 ) = W( 57 ) + a*JVS( 224 )
  a = -W( 41 ) / JVS(          228  )
  W( 41 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 229 )
  W( 50 ) = W( 50 ) + a*JVS( 230 )
  W( 52 ) = W( 52 ) + a*JVS( 231 )
  W( 53 ) = W( 53 ) + a*JVS( 232 )
  W( 55 ) = W( 55 ) + a*JVS( 233 )
  W( 56 ) = W( 56 ) + a*JVS( 234 )
  W( 57 ) = W( 57 ) + a*JVS( 235 )
  a = -W( 42 ) / JVS(          238  )
  W( 42 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 239 )
  W( 44 ) = W( 44 ) + a*JVS( 240 )
  W( 50 ) = W( 50 ) + a*JVS( 241 )
  W( 55 ) = W( 55 ) + a*JVS( 242 )
  W( 56 ) = W( 56 ) + a*JVS( 243 )
  W( 57 ) = W( 57 ) + a*JVS( 244 )
  a = -W( 43 ) / JVS(          247  )
  W( 43 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 248 )
  W( 55 ) = W( 55 ) + a*JVS( 249 )
  W( 56 ) = W( 56 ) + a*JVS( 250 )
  W( 57 ) = W( 57 ) + a*JVS( 251 )
  a = -W( 44 ) / JVS(          256  )
  W( 44 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 257 )
  W( 52 ) = W( 52 ) + a*JVS( 258 )
  W( 53 ) = W( 53 ) + a*JVS( 259 )
  W( 55 ) = W( 55 ) + a*JVS( 260 )
  W( 56 ) = W( 56 ) + a*JVS( 261 )
  W( 57 ) = W( 57 ) + a*JVS( 262 )
  a = -W( 45 ) / JVS(          283  )
  W( 45 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 284 )
  W( 48 ) = W( 48 ) + a*JVS( 285 )
  W( 49 ) = W( 49 ) + a*JVS( 286 )
  W( 50 ) = W( 50 ) + a*JVS( 287 )
  W( 51 ) = W( 51 ) + a*JVS( 288 )
  W( 52 ) = W( 52 ) + a*JVS( 289 )
  W( 53 ) = W( 53 ) + a*JVS( 290 )
  W( 54 ) = W( 54 ) + a*JVS( 291 )
  W( 55 ) = W( 55 ) + a*JVS( 292 )
  W( 56 ) = W( 56 ) + a*JVS( 293 )
  W( 57 ) = W( 57 ) + a*JVS( 294 )
  a = -W( 46 ) / JVS(          304  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 305 )
  W( 50 ) = W( 50 ) + a*JVS( 306 )
  W( 51 ) = W( 51 ) + a*JVS( 307 )
  W( 52 ) = W( 52 ) + a*JVS( 308 )
  W( 53 ) = W( 53 ) + a*JVS( 309 )
  W( 55 ) = W( 55 ) + a*JVS( 310 )
  W( 56 ) = W( 56 ) + a*JVS( 311 )
  W( 57 ) = W( 57 ) + a*JVS( 312 )
  a = -W( 47 ) / JVS(          321  )
  W( 47 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 322 )
  W( 51 ) = W( 51 ) + a*JVS( 323 )
  W( 52 ) = W( 52 ) + a*JVS( 324 )
  W( 53 ) = W( 53 ) + a*JVS( 325 )
  W( 55 ) = W( 55 ) + a*JVS( 326 )
  W( 56 ) = W( 56 ) + a*JVS( 327 )
  W( 57 ) = W( 57 ) + a*JVS( 328 )
  a = -W( 48 ) / JVS(          330  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 331 )
  W( 52 ) = W( 52 ) + a*JVS( 332 )
  W( 53 ) = W( 53 ) + a*JVS( 333 )
  W( 55 ) = W( 55 ) + a*JVS( 334 )
  W( 56 ) = W( 56 ) + a*JVS( 335 )
  W( 57 ) = W( 57 ) + a*JVS( 336 )
  a = -W( 49 ) / JVS(          342  )
  W( 49 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 343 )
  W( 52 ) = W( 52 ) + a*JVS( 344 )
  W( 53 ) = W( 53 ) + a*JVS( 345 )
  W( 54 ) = W( 54 ) + a*JVS( 346 )
  W( 55 ) = W( 55 ) + a*JVS( 347 )
  W( 56 ) = W( 56 ) + a*JVS( 348 )
  W( 57 ) = W( 57 ) + a*JVS( 349 )
  a = -W( 50 ) / JVS(          366  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 367 )
  W( 52 ) = W( 52 ) + a*JVS( 368 )
  W( 53 ) = W( 53 ) + a*JVS( 369 )
  W( 54 ) = W( 54 ) + a*JVS( 370 )
  W( 55 ) = W( 55 ) + a*JVS( 371 )
  W( 56 ) = W( 56 ) + a*JVS( 372 )
  W( 57 ) = W( 57 ) + a*JVS( 373 )
  a = -W( 51 ) / JVS(          398  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 399 )
  W( 53 ) = W( 53 ) + a*JVS( 400 )
  W( 54 ) = W( 54 ) + a*JVS( 401 )
  W( 55 ) = W( 55 ) + a*JVS( 402 )
  W( 56 ) = W( 56 ) + a*JVS( 403 )
  W( 57 ) = W( 57 ) + a*JVS( 404 )
  a = -W( 52 ) / JVS(          418  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 419 )
  W( 54 ) = W( 54 ) + a*JVS( 420 )
  W( 55 ) = W( 55 ) + a*JVS( 421 )
  W( 56 ) = W( 56 ) + a*JVS( 422 )
  W( 57 ) = W( 57 ) + a*JVS( 423 )
  a = -W( 53 ) / JVS(          434  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 435 )
  W( 55 ) = W( 55 ) + a*JVS( 436 )
  W( 56 ) = W( 56 ) + a*JVS( 437 )
  W( 57 ) = W( 57 ) + a*JVS( 438 )
  a = -W( 54 ) / JVS(          464  )
  W( 54 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 465 )
  W( 56 ) = W( 56 ) + a*JVS( 466 )
  W( 57 ) = W( 57 ) + a*JVS( 467 )
  JVS( 468) = W( 1 )
  JVS( 469) = W( 2 )
  JVS( 470) = W( 4 )
  JVS( 471) = W( 7 )
  JVS( 472) = W( 8 )
  JVS( 473) = W( 9 )
  JVS( 474) = W( 10 )
  JVS( 475) = W( 11 )
  JVS( 476) = W( 12 )
  JVS( 477) = W( 13 )
  JVS( 478) = W( 14 )
  JVS( 479) = W( 15 )
  JVS( 480) = W( 16 )
  JVS( 481) = W( 17 )
  JVS( 482) = W( 18 )
  JVS( 483) = W( 19 )
  JVS( 484) = W( 20 )
  JVS( 485) = W( 21 )
  JVS( 486) = W( 22 )
  JVS( 487) = W( 23 )
  JVS( 488) = W( 24 )
  JVS( 489) = W( 25 )
  JVS( 490) = W( 26 )
  JVS( 491) = W( 27 )
  JVS( 492) = W( 28 )
  JVS( 493) = W( 29 )
  JVS( 494) = W( 30 )
  JVS( 495) = W( 31 )
  JVS( 496) = W( 32 )
  JVS( 497) = W( 33 )
  JVS( 498) = W( 34 )
  JVS( 499) = W( 35 )
  JVS( 500) = W( 36 )
  JVS( 501) = W( 37 )
  JVS( 502) = W( 38 )
  JVS( 503) = W( 39 )
  JVS( 504) = W( 40 )
  JVS( 505) = W( 41 )
  JVS( 506) = W( 42 )
  JVS( 507) = W( 43 )
  JVS( 508) = W( 44 )
  JVS( 509) = W( 45 )
  JVS( 510) = W( 46 )
  JVS( 511) = W( 47 )
  JVS( 512) = W( 48 )
  JVS( 513) = W( 49 )
  JVS( 514) = W( 50 )
  JVS( 515) = W( 51 )
  JVS( 516) = W( 52 )
  JVS( 517) = W( 53 )
  JVS( 518) = W( 54 )
  JVS( 519) = W( 55 )
  JVS( 520) = W( 56 )
  JVS( 521) = W( 57 )
  IF ( ABS(  JVS( 537 )) < TINY(a) ) THEN
         IER = 56                                      
         RETURN
  END IF
   W( 34 ) = JVS( 522 )
   W( 36 ) = JVS( 523 )
   W( 40 ) = JVS( 524 )
   W( 42 ) = JVS( 525 )
   W( 43 ) = JVS( 526 )
   W( 44 ) = JVS( 527 )
   W( 47 ) = JVS( 528 )
   W( 48 ) = JVS( 529 )
   W( 49 ) = JVS( 530 )
   W( 50 ) = JVS( 531 )
   W( 51 ) = JVS( 532 )
   W( 52 ) = JVS( 533 )
   W( 53 ) = JVS( 534 )
   W( 54 ) = JVS( 535 )
   W( 55 ) = JVS( 536 )
   W( 56 ) = JVS( 537 )
   W( 57 ) = JVS( 538 )
  a = -W( 34 ) / JVS(          164  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 165 )
  W( 55 ) = W( 55 ) + a*JVS( 166 )
  W( 56 ) = W( 56 ) + a*JVS( 167 )
  W( 57 ) = W( 57 ) + a*JVS( 168 )
  a = -W( 36 ) / JVS(          186  )
  W( 36 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 187 )
  W( 55 ) = W( 55 ) + a*JVS( 188 )
  W( 56 ) = W( 56 ) + a*JVS( 189 )
  W( 57 ) = W( 57 ) + a*JVS( 190 )
  a = -W( 40 ) / JVS(          219  )
  W( 40 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 220 )
  W( 53 ) = W( 53 ) + a*JVS( 221 )
  W( 55 ) = W( 55 ) + a*JVS( 222 )
  W( 56 ) = W( 56 ) + a*JVS( 223 )
  W( 57 ) = W( 57 ) + a*JVS( 224 )
  a = -W( 42 ) / JVS(          238  )
  W( 42 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 239 )
  W( 44 ) = W( 44 ) + a*JVS( 240 )
  W( 50 ) = W( 50 ) + a*JVS( 241 )
  W( 55 ) = W( 55 ) + a*JVS( 242 )
  W( 56 ) = W( 56 ) + a*JVS( 243 )
  W( 57 ) = W( 57 ) + a*JVS( 244 )
  a = -W( 43 ) / JVS(          247  )
  W( 43 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 248 )
  W( 55 ) = W( 55 ) + a*JVS( 249 )
  W( 56 ) = W( 56 ) + a*JVS( 250 )
  W( 57 ) = W( 57 ) + a*JVS( 251 )
  a = -W( 44 ) / JVS(          256  )
  W( 44 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 257 )
  W( 52 ) = W( 52 ) + a*JVS( 258 )
  W( 53 ) = W( 53 ) + a*JVS( 259 )
  W( 55 ) = W( 55 ) + a*JVS( 260 )
  W( 56 ) = W( 56 ) + a*JVS( 261 )
  W( 57 ) = W( 57 ) + a*JVS( 262 )
  a = -W( 47 ) / JVS(          321  )
  W( 47 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 322 )
  W( 51 ) = W( 51 ) + a*JVS( 323 )
  W( 52 ) = W( 52 ) + a*JVS( 324 )
  W( 53 ) = W( 53 ) + a*JVS( 325 )
  W( 55 ) = W( 55 ) + a*JVS( 326 )
  W( 56 ) = W( 56 ) + a*JVS( 327 )
  W( 57 ) = W( 57 ) + a*JVS( 328 )
  a = -W( 48 ) / JVS(          330  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 331 )
  W( 52 ) = W( 52 ) + a*JVS( 332 )
  W( 53 ) = W( 53 ) + a*JVS( 333 )
  W( 55 ) = W( 55 ) + a*JVS( 334 )
  W( 56 ) = W( 56 ) + a*JVS( 335 )
  W( 57 ) = W( 57 ) + a*JVS( 336 )
  a = -W( 49 ) / JVS(          342  )
  W( 49 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 343 )
  W( 52 ) = W( 52 ) + a*JVS( 344 )
  W( 53 ) = W( 53 ) + a*JVS( 345 )
  W( 54 ) = W( 54 ) + a*JVS( 346 )
  W( 55 ) = W( 55 ) + a*JVS( 347 )
  W( 56 ) = W( 56 ) + a*JVS( 348 )
  W( 57 ) = W( 57 ) + a*JVS( 349 )
  a = -W( 50 ) / JVS(          366  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 367 )
  W( 52 ) = W( 52 ) + a*JVS( 368 )
  W( 53 ) = W( 53 ) + a*JVS( 369 )
  W( 54 ) = W( 54 ) + a*JVS( 370 )
  W( 55 ) = W( 55 ) + a*JVS( 371 )
  W( 56 ) = W( 56 ) + a*JVS( 372 )
  W( 57 ) = W( 57 ) + a*JVS( 373 )
  a = -W( 51 ) / JVS(          398  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 399 )
  W( 53 ) = W( 53 ) + a*JVS( 400 )
  W( 54 ) = W( 54 ) + a*JVS( 401 )
  W( 55 ) = W( 55 ) + a*JVS( 402 )
  W( 56 ) = W( 56 ) + a*JVS( 403 )
  W( 57 ) = W( 57 ) + a*JVS( 404 )
  a = -W( 52 ) / JVS(          418  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 419 )
  W( 54 ) = W( 54 ) + a*JVS( 420 )
  W( 55 ) = W( 55 ) + a*JVS( 421 )
  W( 56 ) = W( 56 ) + a*JVS( 422 )
  W( 57 ) = W( 57 ) + a*JVS( 423 )
  a = -W( 53 ) / JVS(          434  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 435 )
  W( 55 ) = W( 55 ) + a*JVS( 436 )
  W( 56 ) = W( 56 ) + a*JVS( 437 )
  W( 57 ) = W( 57 ) + a*JVS( 438 )
  a = -W( 54 ) / JVS(          464  )
  W( 54 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 465 )
  W( 56 ) = W( 56 ) + a*JVS( 466 )
  W( 57 ) = W( 57 ) + a*JVS( 467 )
  a = -W( 55 ) / JVS(          519  )
  W( 55 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 520 )
  W( 57 ) = W( 57 ) + a*JVS( 521 )
  JVS( 522) = W( 34 )
  JVS( 523) = W( 36 )
  JVS( 524) = W( 40 )
  JVS( 525) = W( 42 )
  JVS( 526) = W( 43 )
  JVS( 527) = W( 44 )
  JVS( 528) = W( 47 )
  JVS( 529) = W( 48 )
  JVS( 530) = W( 49 )
  JVS( 531) = W( 50 )
  JVS( 532) = W( 51 )
  JVS( 533) = W( 52 )
  JVS( 534) = W( 53 )
  JVS( 535) = W( 54 )
  JVS( 536) = W( 55 )
  JVS( 537) = W( 56 )
  JVS( 538) = W( 57 )
  IF ( ABS(  JVS( 580 )) < TINY(a) ) THEN
         IER = 57                                      
         RETURN
  END IF
   W( 8 ) = JVS( 539 )
   W( 11 ) = JVS( 540 )
   W( 12 ) = JVS( 541 )
   W( 15 ) = JVS( 542 )
   W( 16 ) = JVS( 543 )
   W( 17 ) = JVS( 544 )
   W( 18 ) = JVS( 545 )
   W( 19 ) = JVS( 546 )
   W( 21 ) = JVS( 547 )
   W( 22 ) = JVS( 548 )
   W( 23 ) = JVS( 549 )
   W( 24 ) = JVS( 550 )
   W( 26 ) = JVS( 551 )
   W( 27 ) = JVS( 552 )
   W( 29 ) = JVS( 553 )
   W( 31 ) = JVS( 554 )
   W( 32 ) = JVS( 555 )
   W( 33 ) = JVS( 556 )
   W( 34 ) = JVS( 557 )
   W( 35 ) = JVS( 558 )
   W( 36 ) = JVS( 559 )
   W( 37 ) = JVS( 560 )
   W( 38 ) = JVS( 561 )
   W( 39 ) = JVS( 562 )
   W( 40 ) = JVS( 563 )
   W( 41 ) = JVS( 564 )
   W( 42 ) = JVS( 565 )
   W( 43 ) = JVS( 566 )
   W( 44 ) = JVS( 567 )
   W( 45 ) = JVS( 568 )
   W( 46 ) = JVS( 569 )
   W( 47 ) = JVS( 570 )
   W( 48 ) = JVS( 571 )
   W( 49 ) = JVS( 572 )
   W( 50 ) = JVS( 573 )
   W( 51 ) = JVS( 574 )
   W( 52 ) = JVS( 575 )
   W( 53 ) = JVS( 576 )
   W( 54 ) = JVS( 577 )
   W( 55 ) = JVS( 578 )
   W( 56 ) = JVS( 579 )
   W( 57 ) = JVS( 580 )
  a = -W( 8 ) / JVS(           19  )
  W( 8 ) = -a
  W( 38 ) = W( 38 ) + a*JVS( 20 )
  W( 53 ) = W( 53 ) + a*JVS( 21 )
  W( 55 ) = W( 55 ) + a*JVS( 22 )
  W( 57 ) = W( 57 ) + a*JVS( 23 )
  a = -W( 11 ) / JVS(           34  )
  W( 11 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 35 )
  W( 55 ) = W( 55 ) + a*JVS( 36 )
  W( 57 ) = W( 57 ) + a*JVS( 37 )
  a = -W( 12 ) / JVS(           38  )
  W( 12 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 39 )
  W( 55 ) = W( 55 ) + a*JVS( 40 )
  W( 57 ) = W( 57 ) + a*JVS( 41 )
  a = -W( 15 ) / JVS(           50  )
  W( 15 ) = -a
  W( 19 ) = W( 19 ) + a*JVS( 51 )
  W( 37 ) = W( 37 ) + a*JVS( 52 )
  W( 38 ) = W( 38 ) + a*JVS( 53 )
  W( 41 ) = W( 41 ) + a*JVS( 54 )
  W( 53 ) = W( 53 ) + a*JVS( 55 )
  W( 55 ) = W( 55 ) + a*JVS( 56 )
  a = -W( 16 ) / JVS(           57  )
  W( 16 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 58 )
  W( 55 ) = W( 55 ) + a*JVS( 59 )
  W( 57 ) = W( 57 ) + a*JVS( 60 )
  a = -W( 17 ) / JVS(           61  )
  W( 17 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 62 )
  W( 55 ) = W( 55 ) + a*JVS( 63 )
  W( 57 ) = W( 57 ) + a*JVS( 64 )
  a = -W( 18 ) / JVS(           65  )
  W( 18 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 66 )
  W( 55 ) = W( 55 ) + a*JVS( 67 )
  W( 56 ) = W( 56 ) + a*JVS( 68 )
  a = -W( 19 ) / JVS(           69  )
  W( 19 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 70 )
  W( 55 ) = W( 55 ) + a*JVS( 71 )
  a = -W( 21 ) / JVS(           76  )
  W( 21 ) = -a
  W( 24 ) = W( 24 ) + a*JVS( 77 )
  W( 53 ) = W( 53 ) + a*JVS( 78 )
  a = -W( 22 ) / JVS(           79  )
  W( 22 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 80 )
  W( 55 ) = W( 55 ) + a*JVS( 81 )
  W( 57 ) = W( 57 ) + a*JVS( 82 )
  a = -W( 23 ) / JVS(           83  )
  W( 23 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 84 )
  W( 55 ) = W( 55 ) + a*JVS( 85 )
  W( 57 ) = W( 57 ) + a*JVS( 86 )
  a = -W( 24 ) / JVS(           88  )
  W( 24 ) = -a
  W( 37 ) = W( 37 ) + a*JVS( 89 )
  W( 53 ) = W( 53 ) + a*JVS( 90 )
  W( 55 ) = W( 55 ) + a*JVS( 91 )
  a = -W( 26 ) / JVS(           99  )
  W( 26 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 100 )
  W( 55 ) = W( 55 ) + a*JVS( 101 )
  W( 57 ) = W( 57 ) + a*JVS( 102 )
  a = -W( 27 ) / JVS(          103  )
  W( 27 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 104 )
  W( 55 ) = W( 55 ) + a*JVS( 105 )
  W( 57 ) = W( 57 ) + a*JVS( 106 )
  a = -W( 29 ) / JVS(          115  )
  W( 29 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 116 )
  W( 55 ) = W( 55 ) + a*JVS( 117 )
  W( 57 ) = W( 57 ) + a*JVS( 118 )
  a = -W( 31 ) / JVS(          130  )
  W( 31 ) = -a
  W( 33 ) = W( 33 ) + a*JVS( 131 )
  W( 35 ) = W( 35 ) + a*JVS( 132 )
  W( 37 ) = W( 37 ) + a*JVS( 133 )
  W( 38 ) = W( 38 ) + a*JVS( 134 )
  W( 41 ) = W( 41 ) + a*JVS( 135 )
  W( 45 ) = W( 45 ) + a*JVS( 136 )
  W( 46 ) = W( 46 ) + a*JVS( 137 )
  W( 49 ) = W( 49 ) + a*JVS( 138 )
  W( 50 ) = W( 50 ) + a*JVS( 139 )
  W( 52 ) = W( 52 ) + a*JVS( 140 )
  W( 53 ) = W( 53 ) + a*JVS( 141 )
  W( 55 ) = W( 55 ) + a*JVS( 142 )
  W( 56 ) = W( 56 ) + a*JVS( 143 )
  W( 57 ) = W( 57 ) + a*JVS( 144 )
  a = -W( 32 ) / JVS(          145  )
  W( 32 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 146 )
  W( 36 ) = W( 36 ) + a*JVS( 147 )
  W( 37 ) = W( 37 ) + a*JVS( 148 )
  W( 43 ) = W( 43 ) + a*JVS( 149 )
  W( 47 ) = W( 47 ) + a*JVS( 150 )
  W( 48 ) = W( 48 ) + a*JVS( 151 )
  W( 49 ) = W( 49 ) + a*JVS( 152 )
  W( 50 ) = W( 50 ) + a*JVS( 153 )
  W( 53 ) = W( 53 ) + a*JVS( 154 )
  W( 55 ) = W( 55 ) + a*JVS( 155 )
  a = -W( 33 ) / JVS(          156  )
  W( 33 ) = -a
  W( 38 ) = W( 38 ) + a*JVS( 157 )
  W( 48 ) = W( 48 ) + a*JVS( 158 )
  W( 52 ) = W( 52 ) + a*JVS( 159 )
  W( 55 ) = W( 55 ) + a*JVS( 160 )
  W( 56 ) = W( 56 ) + a*JVS( 161 )
  a = -W( 34 ) / JVS(          164  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 165 )
  W( 55 ) = W( 55 ) + a*JVS( 166 )
  W( 56 ) = W( 56 ) + a*JVS( 167 )
  W( 57 ) = W( 57 ) + a*JVS( 168 )
  a = -W( 35 ) / JVS(          171  )
  W( 35 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 172 )
  W( 37 ) = W( 37 ) + a*JVS( 173 )
  W( 39 ) = W( 39 ) + a*JVS( 174 )
  W( 41 ) = W( 41 ) + a*JVS( 175 )
  W( 48 ) = W( 48 ) + a*JVS( 176 )
  W( 49 ) = W( 49 ) + a*JVS( 177 )
  W( 50 ) = W( 50 ) + a*JVS( 178 )
  W( 53 ) = W( 53 ) + a*JVS( 179 )
  W( 55 ) = W( 55 ) + a*JVS( 180 )
  W( 56 ) = W( 56 ) + a*JVS( 181 )
  W( 57 ) = W( 57 ) + a*JVS( 182 )
  a = -W( 36 ) / JVS(          186  )
  W( 36 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 187 )
  W( 55 ) = W( 55 ) + a*JVS( 188 )
  W( 56 ) = W( 56 ) + a*JVS( 189 )
  W( 57 ) = W( 57 ) + a*JVS( 190 )
  a = -W( 37 ) / JVS(          191  )
  W( 37 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 192 )
  W( 53 ) = W( 53 ) + a*JVS( 193 )
  W( 55 ) = W( 55 ) + a*JVS( 194 )
  a = -W( 38 ) / JVS(          195  )
  W( 38 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 196 )
  W( 53 ) = W( 53 ) + a*JVS( 197 )
  W( 55 ) = W( 55 ) + a*JVS( 198 )
  a = -W( 39 ) / JVS(          205  )
  W( 39 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 206 )
  W( 48 ) = W( 48 ) + a*JVS( 207 )
  W( 49 ) = W( 49 ) + a*JVS( 208 )
  W( 50 ) = W( 50 ) + a*JVS( 209 )
  W( 52 ) = W( 52 ) + a*JVS( 210 )
  W( 53 ) = W( 53 ) + a*JVS( 211 )
  W( 54 ) = W( 54 ) + a*JVS( 212 )
  W( 55 ) = W( 55 ) + a*JVS( 213 )
  W( 56 ) = W( 56 ) + a*JVS( 214 )
  W( 57 ) = W( 57 ) + a*JVS( 215 )
  a = -W( 40 ) / JVS(          219  )
  W( 40 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 220 )
  W( 53 ) = W( 53 ) + a*JVS( 221 )
  W( 55 ) = W( 55 ) + a*JVS( 222 )
  W( 56 ) = W( 56 ) + a*JVS( 223 )
  W( 57 ) = W( 57 ) + a*JVS( 224 )
  a = -W( 41 ) / JVS(          228  )
  W( 41 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 229 )
  W( 50 ) = W( 50 ) + a*JVS( 230 )
  W( 52 ) = W( 52 ) + a*JVS( 231 )
  W( 53 ) = W( 53 ) + a*JVS( 232 )
  W( 55 ) = W( 55 ) + a*JVS( 233 )
  W( 56 ) = W( 56 ) + a*JVS( 234 )
  W( 57 ) = W( 57 ) + a*JVS( 235 )
  a = -W( 42 ) / JVS(          238  )
  W( 42 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 239 )
  W( 44 ) = W( 44 ) + a*JVS( 240 )
  W( 50 ) = W( 50 ) + a*JVS( 241 )
  W( 55 ) = W( 55 ) + a*JVS( 242 )
  W( 56 ) = W( 56 ) + a*JVS( 243 )
  W( 57 ) = W( 57 ) + a*JVS( 244 )
  a = -W( 43 ) / JVS(          247  )
  W( 43 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 248 )
  W( 55 ) = W( 55 ) + a*JVS( 249 )
  W( 56 ) = W( 56 ) + a*JVS( 250 )
  W( 57 ) = W( 57 ) + a*JVS( 251 )
  a = -W( 44 ) / JVS(          256  )
  W( 44 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 257 )
  W( 52 ) = W( 52 ) + a*JVS( 258 )
  W( 53 ) = W( 53 ) + a*JVS( 259 )
  W( 55 ) = W( 55 ) + a*JVS( 260 )
  W( 56 ) = W( 56 ) + a*JVS( 261 )
  W( 57 ) = W( 57 ) + a*JVS( 262 )
  a = -W( 45 ) / JVS(          283  )
  W( 45 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 284 )
  W( 48 ) = W( 48 ) + a*JVS( 285 )
  W( 49 ) = W( 49 ) + a*JVS( 286 )
  W( 50 ) = W( 50 ) + a*JVS( 287 )
  W( 51 ) = W( 51 ) + a*JVS( 288 )
  W( 52 ) = W( 52 ) + a*JVS( 289 )
  W( 53 ) = W( 53 ) + a*JVS( 290 )
  W( 54 ) = W( 54 ) + a*JVS( 291 )
  W( 55 ) = W( 55 ) + a*JVS( 292 )
  W( 56 ) = W( 56 ) + a*JVS( 293 )
  W( 57 ) = W( 57 ) + a*JVS( 294 )
  a = -W( 46 ) / JVS(          304  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 305 )
  W( 50 ) = W( 50 ) + a*JVS( 306 )
  W( 51 ) = W( 51 ) + a*JVS( 307 )
  W( 52 ) = W( 52 ) + a*JVS( 308 )
  W( 53 ) = W( 53 ) + a*JVS( 309 )
  W( 55 ) = W( 55 ) + a*JVS( 310 )
  W( 56 ) = W( 56 ) + a*JVS( 311 )
  W( 57 ) = W( 57 ) + a*JVS( 312 )
  a = -W( 47 ) / JVS(          321  )
  W( 47 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 322 )
  W( 51 ) = W( 51 ) + a*JVS( 323 )
  W( 52 ) = W( 52 ) + a*JVS( 324 )
  W( 53 ) = W( 53 ) + a*JVS( 325 )
  W( 55 ) = W( 55 ) + a*JVS( 326 )
  W( 56 ) = W( 56 ) + a*JVS( 327 )
  W( 57 ) = W( 57 ) + a*JVS( 328 )
  a = -W( 48 ) / JVS(          330  )
  W( 48 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 331 )
  W( 52 ) = W( 52 ) + a*JVS( 332 )
  W( 53 ) = W( 53 ) + a*JVS( 333 )
  W( 55 ) = W( 55 ) + a*JVS( 334 )
  W( 56 ) = W( 56 ) + a*JVS( 335 )
  W( 57 ) = W( 57 ) + a*JVS( 336 )
  a = -W( 49 ) / JVS(          342  )
  W( 49 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 343 )
  W( 52 ) = W( 52 ) + a*JVS( 344 )
  W( 53 ) = W( 53 ) + a*JVS( 345 )
  W( 54 ) = W( 54 ) + a*JVS( 346 )
  W( 55 ) = W( 55 ) + a*JVS( 347 )
  W( 56 ) = W( 56 ) + a*JVS( 348 )
  W( 57 ) = W( 57 ) + a*JVS( 349 )
  a = -W( 50 ) / JVS(          366  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 367 )
  W( 52 ) = W( 52 ) + a*JVS( 368 )
  W( 53 ) = W( 53 ) + a*JVS( 369 )
  W( 54 ) = W( 54 ) + a*JVS( 370 )
  W( 55 ) = W( 55 ) + a*JVS( 371 )
  W( 56 ) = W( 56 ) + a*JVS( 372 )
  W( 57 ) = W( 57 ) + a*JVS( 373 )
  a = -W( 51 ) / JVS(          398  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 399 )
  W( 53 ) = W( 53 ) + a*JVS( 400 )
  W( 54 ) = W( 54 ) + a*JVS( 401 )
  W( 55 ) = W( 55 ) + a*JVS( 402 )
  W( 56 ) = W( 56 ) + a*JVS( 403 )
  W( 57 ) = W( 57 ) + a*JVS( 404 )
  a = -W( 52 ) / JVS(          418  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 419 )
  W( 54 ) = W( 54 ) + a*JVS( 420 )
  W( 55 ) = W( 55 ) + a*JVS( 421 )
  W( 56 ) = W( 56 ) + a*JVS( 422 )
  W( 57 ) = W( 57 ) + a*JVS( 423 )
  a = -W( 53 ) / JVS(          434  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 435 )
  W( 55 ) = W( 55 ) + a*JVS( 436 )
  W( 56 ) = W( 56 ) + a*JVS( 437 )
  W( 57 ) = W( 57 ) + a*JVS( 438 )
  a = -W( 54 ) / JVS(          464  )
  W( 54 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 465 )
  W( 56 ) = W( 56 ) + a*JVS( 466 )
  W( 57 ) = W( 57 ) + a*JVS( 467 )
  a = -W( 55 ) / JVS(          519  )
  W( 55 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 520 )
  W( 57 ) = W( 57 ) + a*JVS( 521 )
  a = -W( 56 ) / JVS(          537  )
  W( 56 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 538 )
  JVS( 539) = W( 8 )
  JVS( 540) = W( 11 )
  JVS( 541) = W( 12 )
  JVS( 542) = W( 15 )
  JVS( 543) = W( 16 )
  JVS( 544) = W( 17 )
  JVS( 545) = W( 18 )
  JVS( 546) = W( 19 )
  JVS( 547) = W( 21 )
  JVS( 548) = W( 22 )
  JVS( 549) = W( 23 )
  JVS( 550) = W( 24 )
  JVS( 551) = W( 26 )
  JVS( 552) = W( 27 )
  JVS( 553) = W( 29 )
  JVS( 554) = W( 31 )
  JVS( 555) = W( 32 )
  JVS( 556) = W( 33 )
  JVS( 557) = W( 34 )
  JVS( 558) = W( 35 )
  JVS( 559) = W( 36 )
  JVS( 560) = W( 37 )
  JVS( 561) = W( 38 )
  JVS( 562) = W( 39 )
  JVS( 563) = W( 40 )
  JVS( 564) = W( 41 )
  JVS( 565) = W( 42 )
  JVS( 566) = W( 43 )
  JVS( 567) = W( 44 )
  JVS( 568) = W( 45 )
  JVS( 569) = W( 46 )
  JVS( 570) = W( 47 )
  JVS( 571) = W( 48 )
  JVS( 572) = W( 49 )
  JVS( 573) = W( 50 )
  JVS( 574) = W( 51 )
  JVS( 575) = W( 52 )
  JVS( 576) = W( 53 )
  JVS( 577) = W( 54 )
  JVS( 578) = W( 55 )
  JVS( 579) = W( 56 )
  JVS( 580) = W( 57 )
   
   END SUBROUTINE decomp_nmhc9
 


END MODULE nmhc9_Integrator
