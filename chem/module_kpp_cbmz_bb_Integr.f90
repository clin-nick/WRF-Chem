
































MODULE cbmz_bb_Integrator

 USE cbmz_bb_Parameters
 USE cbmz_bb_Precision
 USE cbmz_bb_JacobianSP

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

SUBROUTINE  cbmz_bb_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE cbmz_bb_Parameters

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

   CALL cbmz_bb_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  cbmz_bb_INTEGRATE


SUBROUTINE  cbmz_bb_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE cbmz_bb_Parameters

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
      CALL cbmz_bb_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL cbmz_bb_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = cbmz_bb_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL cbmz_bb_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL cbmz_bb_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL cbmz_bb_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL cbmz_bb_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL cbmz_bb_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL cbmz_bb_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL cbmz_bb_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL cbmz_bb_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL cbmz_bb_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL cbmz_bb_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL cbmz_bb_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL cbmz_bb_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL cbmz_bb_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL cbmz_bb_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL cbmz_bb_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  cbmz_bb_ros_ErrorMsg(Code,T,H,IERR)



   USE cbmz_bb_Precision

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

 END SUBROUTINE  cbmz_bb_ros_ErrorMsg


 SUBROUTINE  cbmz_bb_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL cbmz_bb_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL cbmz_bb_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL cbmz_bb_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL cbmz_bb_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL cbmz_bb_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL cbmz_bb_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL cbmz_bb_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL cbmz_bb_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL cbmz_bb_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL cbmz_bb_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL cbmz_bb_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL cbmz_bb_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL cbmz_bb_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL cbmz_bb_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL cbmz_bb_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL cbmz_bb_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL cbmz_bb_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL cbmz_bb_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL cbmz_bb_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL cbmz_bb_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = cbmz_bb_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL cbmz_bb_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  cbmz_bb_ros_Integrator



  REAL(kind=dp) FUNCTION  cbmz_bb_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    cbmz_bb_ros_ErrorNorm = Err

  END FUNCTION  cbmz_bb_ros_ErrorNorm



  SUBROUTINE cbmz_bb_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL cbmz_bb_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL cbmz_bb_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL cbmz_bb_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  cbmz_bb_ros_FunTimeDeriv



  SUBROUTINE  cbmz_bb_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL cbmz_bb_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL cbmz_bb_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL cbmz_bb_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  cbmz_bb_ros_PrepareMatrix



  SUBROUTINE  cbmz_bb_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 



CALL decomp_cbmz_bb ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  cbmz_bb_ros_Decomp



  SUBROUTINE  cbmz_bb_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL cbmz_bb_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  cbmz_bb_ros_Solve




  SUBROUTINE  cbmz_bb_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  cbmz_bb_Ros2



  SUBROUTINE  cbmz_bb_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cbmz_bb_Ros3





  SUBROUTINE  cbmz_bb_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cbmz_bb_Ros4


  SUBROUTINE  cbmz_bb_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cbmz_bb_Rodas3


  SUBROUTINE  cbmz_bb_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cbmz_bb_Rodas4




END SUBROUTINE  cbmz_bb_Rosenbrock




SUBROUTINE  cbmz_bb_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE cbmz_bb_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL cbmz_bb_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  cbmz_bb_FunTemplate



SUBROUTINE  cbmz_bb_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE cbmz_bb_Parameters
 
 USE cbmz_bb_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL cbmz_bb_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  cbmz_bb_JacTemplate

















SUBROUTINE cbmz_bb_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(52)
  A(2) = RCT(2)*V(54)
  A(3) = RCT(3)*V(22)
  A(4) = RCT(4)*V(27)
  A(5) = RCT(5)*V(16)
  A(6) = RCT(6)*V(12)
  A(7) = RCT(7)*V(51)
  A(8) = RCT(8)*V(51)
  A(9) = RCT(9)*V(9)
  A(10) = RCT(10)*V(36)
  A(11) = RCT(11)*V(36)
  A(12) = RCT(12)*V(20)
  A(13) = RCT(13)*V(21)
  A(14) = RCT(14)*V(43)
  A(15) = RCT(15)*V(38)
  A(16) = RCT(16)*V(42)
  A(17) = RCT(17)*V(30)
  A(18) = RCT(18)*V(47)
  A(19) = RCT(19)*V(49)
  A(20) = RCT(20)*V(45)
  A(21) = RCT(21)*V(6)*F(2)
  A(22) = RCT(22)*V(6)*F(1)
  A(23) = RCT(23)*V(25)*F(2)
  A(24) = RCT(24)*V(25)*V(51)
  A(25) = RCT(25)*V(25)*V(52)
  A(26) = RCT(26)*V(25)*V(52)
  A(27) = RCT(27)*V(25)*V(56)
  A(28) = RCT(28)*V(51)*V(56)
  A(29) = RCT(29)*V(51)*V(52)
  A(30) = RCT(30)*V(50)*V(51)
  A(31) = RCT(31)*V(51)*V(55)
  A(32) = RCT(32)*V(50)*F(2)
  A(33) = RCT(33)*V(50)*V(56)
  A(34) = RCT(34)*V(50)*V(52)
  A(35) = RCT(35)*V(50)*V(54)
  A(36) = RCT(36)*V(22)*V(50)
  A(37) = RCT(37)*V(27)*V(50)
  A(38) = RCT(38)*V(16)*V(50)
  A(39) = RCT(39)*V(50)*V(55)
  A(40) = RCT(40)*V(9)*V(50)
  A(41) = RCT(41)*V(55)*V(55)
  A(42) = RCT(42)*V(55)*V(55)*F(1)
  A(43) = RCT(43)*V(55)*V(56)
  A(44) = RCT(44)*V(52)*V(55)
  A(45) = RCT(45)*V(52)*V(55)
  A(46) = RCT(46)*V(16)
  A(47) = RCT(47)*V(54)*V(56)
  A(48) = RCT(48)*V(52)*V(54)
  A(49) = RCT(49)*V(52)*V(54)
  A(50) = RCT(50)*V(54)*V(54)
  A(51) = RCT(51)*V(54)*V(55)
  A(52) = RCT(52)*V(12)*F(1)
  A(53) = RCT(53)*V(12)
  A(54) = RCT(54)*V(28)*V(50)
  A(55) = RCT(55)*V(7)*V(50)
  A(56) = RCT(56)*V(14)*V(50)
  A(57) = RCT(57)*V(17)*V(50)
  A(58) = RCT(58)*V(29)*V(50)
  A(59) = RCT(59)*V(36)*V(50)
  A(60) = RCT(60)*V(36)*V(54)
  A(61) = RCT(61)*V(43)*V(50)
  A(62) = RCT(62)*V(43)*V(54)
  A(63) = RCT(63)*V(38)*V(50)
  A(64) = RCT(64)*V(42)*V(50)
  A(65) = RCT(65)*V(42)*V(54)
  A(66) = RCT(66)*V(23)*V(51)
  A(67) = RCT(67)*V(23)*V(50)
  A(68) = RCT(68)*V(34)*V(51)
  A(69) = RCT(69)*V(39)*V(51)
  A(70) = RCT(70)*V(34)*V(50)
  A(71) = RCT(71)*V(39)*V(50)
  A(72) = RCT(72)*V(34)*V(54)
  A(73) = RCT(73)*V(39)*V(54)
  A(74) = RCT(74)*V(11)*V(50)
  A(75) = RCT(75)*V(13)*V(50)
  A(76) = RCT(76)*V(18)*V(56)
  A(77) = RCT(77)*V(26)*V(50)
  A(78) = RCT(78)*V(26)*V(54)
  A(79) = RCT(79)*V(15)*V(52)
  A(80) = RCT(80)*V(30)*V(50)
  A(81) = RCT(81)*V(30)*V(51)
  A(82) = RCT(82)*V(35)*V(50)
  A(83) = RCT(83)*V(35)*V(51)
  A(84) = RCT(84)*V(35)*V(54)
  A(85) = RCT(85)*V(45)*V(50)
  A(86) = RCT(86)*V(45)*V(51)
  A(87) = RCT(87)*V(45)*V(54)
  A(88) = RCT(88)*V(32)*V(56)
  A(89) = RCT(89)*V(31)*V(56)
  A(90) = RCT(90)*V(33)*V(56)
  A(91) = RCT(91)*V(32)*V(55)
  A(92) = RCT(92)*V(31)*V(55)
  A(93) = RCT(93)*V(33)*V(55)
  A(94) = RCT(94)*V(20)*V(50)
  A(95) = RCT(95)*V(21)*V(50)
  A(96) = RCT(96)*V(47)*V(50)
  A(97) = RCT(97)*V(49)*V(50)
  A(98) = RCT(98)*V(52)*V(53)
  A(99) = RCT(99)*V(10)
  A(100) = RCT(100)*V(24)*V(50)
  A(101) = RCT(101)*V(44)*V(56)
  A(102) = RCT(102)*V(40)*V(56)
  A(103) = RCT(103)*V(48)*V(56)
  A(104) = RCT(104)*V(53)*V(56)
  A(105) = RCT(105)*V(46)*V(56)
  A(106) = RCT(106)*V(41)*V(56)
  A(107) = RCT(107)*V(37)*V(56)
  A(108) = RCT(108)*V(44)*V(54)
  A(109) = RCT(109)*V(40)*V(54)
  A(110) = RCT(110)*V(48)*V(54)
  A(111) = RCT(111)*V(53)*V(54)
  A(112) = RCT(112)*V(46)*V(54)
  A(113) = RCT(113)*V(41)*V(54)
  A(114) = RCT(114)*V(37)*V(54)
  A(115) = RCT(115)*V(44)*V(55)
  A(116) = RCT(116)*V(40)*V(55)
  A(117) = RCT(117)*V(48)*V(55)
  A(118) = RCT(118)*V(53)*V(55)
  A(119) = RCT(119)*V(46)*V(55)
  A(120) = RCT(120)*V(41)*V(55)
  A(121) = RCT(121)*V(37)*V(55)
  A(122) = RCT(122)*V(8)*V(50)
  A(123) = RCT(123)*V(19)*V(29)
  A(127) = RCT(127)*V(44)
  A(128) = RCT(128)*V(40)
  A(129) = RCT(129)*V(48)
  A(130) = RCT(130)*V(53)
  A(131) = RCT(131)*V(46)
  A(132) = RCT(132)*V(41)
  A(133) = RCT(133)*V(32)
  A(134) = RCT(134)*V(31)
  A(135) = RCT(135)*V(33)
  A(136) = RCT(136)*V(37)


  Vdot(1) = A(55)
  Vdot(2) = 0
  Vdot(3) = 0
  Vdot(4) = 0.52*A(66)+0.22*A(68)
  Vdot(5) = 0.09*A(68)+0.16*A(69)+0.39*A(83)+0.46*A(86)+0.4*A(118)
  Vdot(6) = A(8)-A(21)-A(22)
  Vdot(7) = -A(55)
  Vdot(8) = -A(122)
  Vdot(9) = -A(9)-A(40)+A(41)+A(42)
  Vdot(10) = A(98)-A(99)
  Vdot(11) = -A(74)
  Vdot(12) = -A(6)+A(49)-A(52)-A(53)
  Vdot(13) = -A(75)
  Vdot(14) = -A(56)+0.06*A(68)+0.08*A(69)
  Vdot(15) = 0.4*A(77)+A(78)-A(79)
  Vdot(16) = -A(5)-A(38)+A(44)-A(46)
  Vdot(17) = -A(57)+0.01*A(68)+0.01*A(69)+0.2*A(128)
  Vdot(18) = 0.8*A(74)+0.45*A(75)-A(76)
  Vdot(19) = 1.98*A(18)+1.98*A(19)+1.06*A(68)+2.26*A(69)+A(70)+2.23*A(71)+0.42*A(96)+1.68*A(103)+A(106)+1.98*A(110)&
               &+A(113)-A(123)+1.25*A(129)+A(132)
  Vdot(20) = -A(12)-A(94)+A(115)
  Vdot(21) = -A(13)-A(95)+A(116)
  Vdot(22) = -A(3)+A(33)-A(36)+A(45)
  Vdot(23) = -A(66)-A(67)
  Vdot(24) = 0.03*A(68)+0.04*A(69)-A(100)+0.34*A(127)
  Vdot(25) = A(1)+0.89*A(2)+A(7)+A(21)-A(23)-A(24)-A(25)-A(26)-A(27)
  Vdot(26) = 0.12*A(74)+0.05*A(75)-A(77)-A(78)
  Vdot(27) = -A(4)+A(34)-A(37)+0.3*A(51)+2*A(52)+A(60)+A(62)+A(65)+A(78)+0.07*A(87)
  Vdot(28) = A(10)+A(11)+A(14)+A(16)+A(17)+0.33*A(20)-A(54)+A(59)+A(60)+A(65)+0.24*A(66)+0.31*A(68)+0.3*A(69)+2*A(80)&
               &+0.69*A(81)+0.07*A(83)+0.16*A(86)+0.64*A(87)+0.59*A(90)
  Vdot(29) = -A(58)+1.1*A(75)+1.86*A(87)+0.18*A(88)+1.6*A(89)+2*A(92)-A(123)+2*A(134)
  Vdot(30) = -A(17)+0.95*A(76)+0.3*A(77)-A(80)-A(81)
  Vdot(31) = A(84)-A(89)-A(92)-A(134)
  Vdot(32) = A(82)-A(88)-A(91)-A(133)
  Vdot(33) = 0.5*A(85)-A(90)-A(93)-A(135)
  Vdot(34) = -A(68)-A(70)-A(72)
  Vdot(35) = -A(82)-A(83)-A(84)
  Vdot(36) = -A(10)-A(11)+A(12)+0.2*A(20)-A(59)-A(60)+A(66)+1.56*A(67)+0.57*A(68)+A(70)+A(80)+0.7*A(81)+0.6*A(83)+0.15&
               &*A(86)+0.28*A(87)+0.63*A(88)+0.25*A(90)+0.3*A(94)+A(100)+A(101)+A(105)+0.5*A(106)+A(108)+A(112)+0.5*A(113)&
               &+0.66*A(127)+0.7*A(131)+0.5*A(132)
  Vdot(37) = 0.4*A(18)+0.41*A(19)+A(64)+A(67)+A(70)+A(71)+0.08*A(74)+0.5*A(75)+0.6*A(77)+A(80)+0.03*A(81)+0.08*A(82)+0.2&
               &*A(83)+0.2*A(85)+0.07*A(86)+0.93*A(87)+0.34*A(103)-A(107)+0.4*A(110)-A(114)-A(121)+0.24*A(129)-A(136)
  Vdot(38) = -A(15)+0.74*A(18)+0.74*A(19)+0.03*A(20)-A(63)+0.07*A(69)+0.23*A(71)+0.09*A(86)+0.63*A(90)+0.62*A(103)+0.74&
               &*A(110)+0.57*A(129)+0.15*A(131)+0.5*A(135)
  Vdot(39) = -A(69)-A(71)-A(73)
  Vdot(40) = 0.1*A(18)+0.1*A(19)+A(57)+0.06*A(68)+0.05*A(69)+0.7*A(95)-A(102)+0.08*A(103)-A(109)+0.1*A(110)-A(116)&
               &-A(128)+0.06*A(129)
  Vdot(41) = A(72)+A(73)+A(97)-A(106)-A(113)-A(120)-A(132)
  Vdot(42) = -A(16)-A(64)-A(65)+0.04*A(68)+0.07*A(69)+0.8*A(75)+0.2*A(81)+0.85*A(86)+0.34*A(90)+0.19*A(96)+0.15*A(131)
  Vdot(43) = A(13)-A(14)+0.3*A(18)+0.3*A(19)+0.07*A(20)-A(61)-A(62)+0.22*A(67)+0.47*A(68)+1.03*A(69)+A(70)+1.77*A(71)&
               &+0.03*A(81)+0.15*A(83)+0.02*A(86)+0.28*A(87)+0.8*A(89)+0.55*A(90)+0.3*A(95)+0.04*A(96)+A(102)+0.25*A(103)&
               &+0.5*A(106)+A(109)+0.3*A(110)+0.5*A(113)+A(122)+0.8*A(128)+0.21*A(129)+0.5*A(132)+A(134)+0.5*A(135)
  Vdot(44) = A(14)+A(15)+0.7*A(20)+A(56)+0.07*A(68)+0.1*A(69)+0.05*A(86)+0.7*A(94)-A(101)+A(104)-A(108)+A(111)-A(115)&
               &-A(127)+A(130)
  Vdot(45) = -A(20)+0.65*A(83)-A(85)-A(86)-A(87)+0.91*A(88)+0.2*A(89)+A(133)
  Vdot(46) = A(63)+0.11*A(69)-A(105)-A(112)-A(119)-A(131)
  Vdot(47) = -A(18)+A(91)+A(93)-A(96)+A(117)+A(119)
  Vdot(48) = A(58)+0.03*A(68)+0.09*A(69)+0.77*A(96)-A(103)-A(110)-A(117)-A(129)
  Vdot(49) = -A(19)+0.05*A(76)+A(79)+0.93*A(87)+0.09*A(88)+0.8*A(89)+A(92)-A(97)+0.16*A(103)+0.5*A(106)+0.5*A(113)&
               &+A(120)+0.5*A(132)+A(134)
  Vdot(50) = A(3)+A(4)+2*A(9)+A(12)+A(13)+A(18)+2*A(22)-A(30)+A(31)-A(32)-A(33)-A(34)-A(35)-A(36)-A(37)-A(38)-A(39)&
               &-A(40)+A(43)+0.7*A(51)-A(54)-A(55)-A(56)-A(57)-A(58)-A(59)-A(61)-A(63)-A(64)+0.12*A(66)-A(67)+0.33*A(68)+0.6&
               &*A(69)-A(70)-A(71)-A(74)-A(75)-A(77)-A(80)+0.08*A(81)-A(82)+0.27*A(83)-A(85)+0.27*A(86)-0.7*A(94)-0.7*A(95)&
               &-0.77*A(96)-A(97)-A(100)-A(122)
  Vdot(51) = -A(7)-A(8)+A(23)-A(24)-A(28)-A(29)-A(30)-A(31)-A(66)-A(68)-A(69)-A(81)-A(83)-A(86)+0.4*A(118)
  Vdot(52) = -A(1)+0.89*A(2)+A(4)+A(5)+A(6)+A(19)-A(25)-A(26)+A(27)+A(28)-A(29)-A(34)+A(35)+A(36)+A(38)+A(43)-A(44)&
               &-A(45)+A(46)+2*A(47)-A(49)+2*A(50)+0.7*A(51)+A(53)+0.95*A(76)-A(79)+0.91*A(88)+1.2*A(89)+A(90)-A(98)+A(99)&
               &+A(101)+A(102)+0.84*A(103)+A(104)+A(105)+1.5*A(106)+A(107)+A(108)+A(109)+A(110)+A(111)+A(112)+1.5*A(113)&
               &+A(114)+0.5*A(132)
  Vdot(53) = A(15)+A(16)+A(17)+0.97*A(20)+A(61)+A(62)+A(64)+A(65)+0.13*A(68)+0.19*A(69)+A(80)+0.62*A(81)+0.2*A(83)+0.5&
               &*A(85)+0.11*A(86)+0.07*A(87)-A(98)+A(99)-A(104)+A(105)-A(111)+A(112)-A(118)-A(130)+0.7*A(131)
  Vdot(54) = -A(2)+A(6)+A(26)+A(29)-A(35)+A(37)-A(47)-A(48)-A(49)-2*A(50)-A(51)+A(53)-A(60)-A(62)-A(65)-A(72)-A(73)&
               &-A(78)-A(84)-A(87)-A(108)-A(109)-A(110)-A(111)-A(112)-A(113)-A(114)
  Vdot(55) = A(5)+2*A(10)+A(12)+A(13)+A(14)+A(16)+A(17)+0.9*A(18)+0.9*A(19)+0.33*A(20)+A(30)-A(31)+A(32)+A(35)-A(39)&
               &+A(40)-2*A(41)-2*A(42)-A(43)-A(44)-A(45)+A(46)-A(51)+A(54)+A(55)+A(59)+A(60)+0.22*A(66)+A(67)+0.26*A(68)&
               &+0.22*A(69)+A(70)+A(71)+0.2*A(74)+0.55*A(75)+0.95*A(76)+0.6*A(77)+2*A(80)+0.76*A(81)+0.07*A(83)+0.1*A(86)&
               &+0.93*A(87)+0.91*A(88)+0.8*A(89)+A(90)-A(91)-A(92)-A(93)+A(100)+A(101)+A(102)+0.76*A(103)+0.5*A(106)+A(108)&
               &+A(109)+0.9*A(110)+0.5*A(113)-A(115)-A(116)-A(117)-A(118)-A(119)-A(120)-A(121)+A(122)+0.32*A(127)+0.6*A(128)&
               &+0.54*A(129)
  Vdot(56) = A(1)+0.11*A(2)+A(3)+A(25)-A(27)-A(28)-A(33)-A(43)-A(47)+A(48)-A(76)-A(88)-A(89)-A(90)-A(101)-A(102)-A(103)&
               &-A(104)-A(105)-A(106)-A(107)
      
END SUBROUTINE cbmz_bb_Fun
















SUBROUTINE cbmz_bb_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(52)
  IRR(2) = RCT(2)*V(54)
  IRR(3) = RCT(3)*V(22)
  IRR(4) = RCT(4)*V(27)
  IRR(5) = RCT(5)*V(16)
  IRR(6) = RCT(6)*V(12)
  IRR(7) = RCT(7)*V(51)
  IRR(8) = RCT(8)*V(51)
  IRR(9) = RCT(9)*V(9)
  IRR(10) = RCT(10)*V(36)
  IRR(11) = RCT(11)*V(36)
  IRR(12) = RCT(12)*V(20)
  IRR(13) = RCT(13)*V(21)
  IRR(14) = RCT(14)*V(43)
  IRR(15) = RCT(15)*V(38)
  IRR(16) = RCT(16)*V(42)
  IRR(17) = RCT(17)*V(30)
  IRR(18) = RCT(18)*V(47)
  IRR(19) = RCT(19)*V(49)
  IRR(20) = RCT(20)*V(45)
  IRR(21) = RCT(21)*V(6)*F(2)
  IRR(22) = RCT(22)*V(6)*F(1)
  IRR(23) = RCT(23)*V(25)*F(2)
  IRR(24) = RCT(24)*V(25)*V(51)
  IRR(25) = RCT(25)*V(25)*V(52)
  IRR(26) = RCT(26)*V(25)*V(52)
  IRR(27) = RCT(27)*V(25)*V(56)
  IRR(28) = RCT(28)*V(51)*V(56)
  IRR(29) = RCT(29)*V(51)*V(52)
  IRR(30) = RCT(30)*V(50)*V(51)
  IRR(31) = RCT(31)*V(51)*V(55)
  IRR(32) = RCT(32)*V(50)*F(2)
  IRR(33) = RCT(33)*V(50)*V(56)
  IRR(34) = RCT(34)*V(50)*V(52)
  IRR(35) = RCT(35)*V(50)*V(54)
  IRR(36) = RCT(36)*V(22)*V(50)
  IRR(37) = RCT(37)*V(27)*V(50)
  IRR(38) = RCT(38)*V(16)*V(50)
  IRR(39) = RCT(39)*V(50)*V(55)
  IRR(40) = RCT(40)*V(9)*V(50)
  IRR(41) = RCT(41)*V(55)*V(55)
  IRR(42) = RCT(42)*V(55)*V(55)*F(1)
  IRR(43) = RCT(43)*V(55)*V(56)
  IRR(44) = RCT(44)*V(52)*V(55)
  IRR(45) = RCT(45)*V(52)*V(55)
  IRR(46) = RCT(46)*V(16)
  IRR(47) = RCT(47)*V(54)*V(56)
  IRR(48) = RCT(48)*V(52)*V(54)
  IRR(49) = RCT(49)*V(52)*V(54)
  IRR(50) = RCT(50)*V(54)*V(54)
  IRR(51) = RCT(51)*V(54)*V(55)
  IRR(52) = RCT(52)*V(12)*F(1)
  IRR(53) = RCT(53)*V(12)
  IRR(54) = RCT(54)*V(28)*V(50)
  IRR(55) = RCT(55)*V(7)*V(50)
  IRR(56) = RCT(56)*V(14)*V(50)
  IRR(57) = RCT(57)*V(17)*V(50)
  IRR(58) = RCT(58)*V(29)*V(50)
  IRR(59) = RCT(59)*V(36)*V(50)
  IRR(60) = RCT(60)*V(36)*V(54)
  IRR(61) = RCT(61)*V(43)*V(50)
  IRR(62) = RCT(62)*V(43)*V(54)
  IRR(63) = RCT(63)*V(38)*V(50)
  IRR(64) = RCT(64)*V(42)*V(50)
  IRR(65) = RCT(65)*V(42)*V(54)
  IRR(66) = RCT(66)*V(23)*V(51)
  IRR(67) = RCT(67)*V(23)*V(50)
  IRR(68) = RCT(68)*V(34)*V(51)
  IRR(69) = RCT(69)*V(39)*V(51)
  IRR(70) = RCT(70)*V(34)*V(50)
  IRR(71) = RCT(71)*V(39)*V(50)
  IRR(72) = RCT(72)*V(34)*V(54)
  IRR(73) = RCT(73)*V(39)*V(54)
  IRR(74) = RCT(74)*V(11)*V(50)
  IRR(75) = RCT(75)*V(13)*V(50)
  IRR(76) = RCT(76)*V(18)*V(56)
  IRR(77) = RCT(77)*V(26)*V(50)
  IRR(78) = RCT(78)*V(26)*V(54)
  IRR(79) = RCT(79)*V(15)*V(52)
  IRR(80) = RCT(80)*V(30)*V(50)
  IRR(81) = RCT(81)*V(30)*V(51)
  IRR(82) = RCT(82)*V(35)*V(50)
  IRR(83) = RCT(83)*V(35)*V(51)
  IRR(84) = RCT(84)*V(35)*V(54)
  IRR(85) = RCT(85)*V(45)*V(50)
  IRR(86) = RCT(86)*V(45)*V(51)
  IRR(87) = RCT(87)*V(45)*V(54)
  IRR(88) = RCT(88)*V(32)*V(56)
  IRR(89) = RCT(89)*V(31)*V(56)
  IRR(90) = RCT(90)*V(33)*V(56)
  IRR(91) = RCT(91)*V(32)*V(55)
  IRR(92) = RCT(92)*V(31)*V(55)
  IRR(93) = RCT(93)*V(33)*V(55)
  IRR(94) = RCT(94)*V(20)*V(50)
  IRR(95) = RCT(95)*V(21)*V(50)
  IRR(96) = RCT(96)*V(47)*V(50)
  IRR(97) = RCT(97)*V(49)*V(50)
  IRR(98) = RCT(98)*V(52)*V(53)
  IRR(99) = RCT(99)*V(10)
  IRR(100) = RCT(100)*V(24)*V(50)
  IRR(101) = RCT(101)*V(44)*V(56)
  IRR(102) = RCT(102)*V(40)*V(56)
  IRR(103) = RCT(103)*V(48)*V(56)
  IRR(104) = RCT(104)*V(53)*V(56)
  IRR(105) = RCT(105)*V(46)*V(56)
  IRR(106) = RCT(106)*V(41)*V(56)
  IRR(107) = RCT(107)*V(37)*V(56)
  IRR(108) = RCT(108)*V(44)*V(54)
  IRR(109) = RCT(109)*V(40)*V(54)
  IRR(110) = RCT(110)*V(48)*V(54)
  IRR(111) = RCT(111)*V(53)*V(54)
  IRR(112) = RCT(112)*V(46)*V(54)
  IRR(113) = RCT(113)*V(41)*V(54)
  IRR(114) = RCT(114)*V(37)*V(54)
  IRR(115) = RCT(115)*V(44)*V(55)
  IRR(116) = RCT(116)*V(40)*V(55)
  IRR(117) = RCT(117)*V(48)*V(55)
  IRR(118) = RCT(118)*V(53)*V(55)
  IRR(119) = RCT(119)*V(46)*V(55)
  IRR(120) = RCT(120)*V(41)*V(55)
  IRR(121) = RCT(121)*V(37)*V(55)
  IRR(122) = RCT(122)*V(8)*V(50)
  IRR(123) = RCT(123)*V(19)*V(29)
  IRR(127) = RCT(127)*V(44)
  IRR(128) = RCT(128)*V(40)
  IRR(129) = RCT(129)*V(48)
  IRR(130) = RCT(130)*V(53)
  IRR(131) = RCT(131)*V(46)
  IRR(132) = RCT(132)*V(41)
  IRR(133) = RCT(133)*V(32)
  IRR(134) = RCT(134)*V(31)
  IRR(135) = RCT(135)*V(33)
  IRR(136) = RCT(136)*V(37)
      
END SUBROUTINE cbmz_bb_IRRFun
















SUBROUTINE cbmz_bb_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(234)


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

  B(21) = RCT(21)*F(2)

  B(23) = RCT(22)*F(1)

  B(25) = RCT(23)*F(2)

  B(27) = RCT(24)*V(51)

  B(28) = RCT(24)*V(25)

  B(29) = RCT(25)*V(52)

  B(30) = RCT(25)*V(25)

  B(31) = RCT(26)*V(52)

  B(32) = RCT(26)*V(25)

  B(33) = RCT(27)*V(56)

  B(34) = RCT(27)*V(25)

  B(35) = RCT(28)*V(56)

  B(36) = RCT(28)*V(51)

  B(37) = RCT(29)*V(52)

  B(38) = RCT(29)*V(51)

  B(39) = RCT(30)*V(51)

  B(40) = RCT(30)*V(50)

  B(41) = RCT(31)*V(55)

  B(42) = RCT(31)*V(51)

  B(43) = RCT(32)*F(2)

  B(45) = RCT(33)*V(56)

  B(46) = RCT(33)*V(50)

  B(47) = RCT(34)*V(52)

  B(48) = RCT(34)*V(50)

  B(49) = RCT(35)*V(54)

  B(50) = RCT(35)*V(50)

  B(51) = RCT(36)*V(50)

  B(52) = RCT(36)*V(22)

  B(53) = RCT(37)*V(50)

  B(54) = RCT(37)*V(27)

  B(55) = RCT(38)*V(50)

  B(56) = RCT(38)*V(16)

  B(57) = RCT(39)*V(55)

  B(58) = RCT(39)*V(50)

  B(59) = RCT(40)*V(50)

  B(60) = RCT(40)*V(9)

  B(61) = RCT(41)*2*V(55)

  B(62) = RCT(42)*2*V(55)*F(1)

  B(64) = RCT(43)*V(56)

  B(65) = RCT(43)*V(55)

  B(66) = RCT(44)*V(55)

  B(67) = RCT(44)*V(52)

  B(68) = RCT(45)*V(55)

  B(69) = RCT(45)*V(52)

  B(70) = RCT(46)

  B(71) = RCT(47)*V(56)

  B(72) = RCT(47)*V(54)

  B(73) = RCT(48)*V(54)

  B(74) = RCT(48)*V(52)

  B(75) = RCT(49)*V(54)

  B(76) = RCT(49)*V(52)

  B(77) = RCT(50)*2*V(54)

  B(78) = RCT(51)*V(55)

  B(79) = RCT(51)*V(54)

  B(80) = RCT(52)*F(1)

  B(82) = RCT(53)

  B(83) = RCT(54)*V(50)

  B(84) = RCT(54)*V(28)

  B(85) = RCT(55)*V(50)

  B(86) = RCT(55)*V(7)

  B(87) = RCT(56)*V(50)

  B(88) = RCT(56)*V(14)

  B(89) = RCT(57)*V(50)

  B(90) = RCT(57)*V(17)

  B(91) = RCT(58)*V(50)

  B(92) = RCT(58)*V(29)

  B(93) = RCT(59)*V(50)

  B(94) = RCT(59)*V(36)

  B(95) = RCT(60)*V(54)

  B(96) = RCT(60)*V(36)

  B(97) = RCT(61)*V(50)

  B(98) = RCT(61)*V(43)

  B(99) = RCT(62)*V(54)

  B(100) = RCT(62)*V(43)

  B(101) = RCT(63)*V(50)

  B(102) = RCT(63)*V(38)

  B(103) = RCT(64)*V(50)

  B(104) = RCT(64)*V(42)

  B(105) = RCT(65)*V(54)

  B(106) = RCT(65)*V(42)

  B(107) = RCT(66)*V(51)

  B(108) = RCT(66)*V(23)

  B(109) = RCT(67)*V(50)

  B(110) = RCT(67)*V(23)

  B(111) = RCT(68)*V(51)

  B(112) = RCT(68)*V(34)

  B(113) = RCT(69)*V(51)

  B(114) = RCT(69)*V(39)

  B(115) = RCT(70)*V(50)

  B(116) = RCT(70)*V(34)

  B(117) = RCT(71)*V(50)

  B(118) = RCT(71)*V(39)

  B(119) = RCT(72)*V(54)

  B(120) = RCT(72)*V(34)

  B(121) = RCT(73)*V(54)

  B(122) = RCT(73)*V(39)

  B(123) = RCT(74)*V(50)

  B(124) = RCT(74)*V(11)

  B(125) = RCT(75)*V(50)

  B(126) = RCT(75)*V(13)

  B(127) = RCT(76)*V(56)

  B(128) = RCT(76)*V(18)

  B(129) = RCT(77)*V(50)

  B(130) = RCT(77)*V(26)

  B(131) = RCT(78)*V(54)

  B(132) = RCT(78)*V(26)

  B(133) = RCT(79)*V(52)

  B(134) = RCT(79)*V(15)

  B(135) = RCT(80)*V(50)

  B(136) = RCT(80)*V(30)

  B(137) = RCT(81)*V(51)

  B(138) = RCT(81)*V(30)

  B(139) = RCT(82)*V(50)

  B(140) = RCT(82)*V(35)

  B(141) = RCT(83)*V(51)

  B(142) = RCT(83)*V(35)

  B(143) = RCT(84)*V(54)

  B(144) = RCT(84)*V(35)

  B(145) = RCT(85)*V(50)

  B(146) = RCT(85)*V(45)

  B(147) = RCT(86)*V(51)

  B(148) = RCT(86)*V(45)

  B(149) = RCT(87)*V(54)

  B(150) = RCT(87)*V(45)

  B(151) = RCT(88)*V(56)

  B(152) = RCT(88)*V(32)

  B(153) = RCT(89)*V(56)

  B(154) = RCT(89)*V(31)

  B(155) = RCT(90)*V(56)

  B(156) = RCT(90)*V(33)

  B(157) = RCT(91)*V(55)

  B(158) = RCT(91)*V(32)

  B(159) = RCT(92)*V(55)

  B(160) = RCT(92)*V(31)

  B(161) = RCT(93)*V(55)

  B(162) = RCT(93)*V(33)

  B(163) = RCT(94)*V(50)

  B(164) = RCT(94)*V(20)

  B(165) = RCT(95)*V(50)

  B(166) = RCT(95)*V(21)

  B(167) = RCT(96)*V(50)

  B(168) = RCT(96)*V(47)

  B(169) = RCT(97)*V(50)

  B(170) = RCT(97)*V(49)

  B(171) = RCT(98)*V(53)

  B(172) = RCT(98)*V(52)

  B(173) = RCT(99)

  B(174) = RCT(100)*V(50)

  B(175) = RCT(100)*V(24)

  B(176) = RCT(101)*V(56)

  B(177) = RCT(101)*V(44)

  B(178) = RCT(102)*V(56)

  B(179) = RCT(102)*V(40)

  B(180) = RCT(103)*V(56)

  B(181) = RCT(103)*V(48)

  B(182) = RCT(104)*V(56)

  B(183) = RCT(104)*V(53)

  B(184) = RCT(105)*V(56)

  B(185) = RCT(105)*V(46)

  B(186) = RCT(106)*V(56)

  B(187) = RCT(106)*V(41)

  B(188) = RCT(107)*V(56)

  B(189) = RCT(107)*V(37)

  B(190) = RCT(108)*V(54)

  B(191) = RCT(108)*V(44)

  B(192) = RCT(109)*V(54)

  B(193) = RCT(109)*V(40)

  B(194) = RCT(110)*V(54)

  B(195) = RCT(110)*V(48)

  B(196) = RCT(111)*V(54)

  B(197) = RCT(111)*V(53)

  B(198) = RCT(112)*V(54)

  B(199) = RCT(112)*V(46)

  B(200) = RCT(113)*V(54)

  B(201) = RCT(113)*V(41)

  B(202) = RCT(114)*V(54)

  B(203) = RCT(114)*V(37)

  B(204) = RCT(115)*V(55)

  B(205) = RCT(115)*V(44)

  B(206) = RCT(116)*V(55)

  B(207) = RCT(116)*V(40)

  B(208) = RCT(117)*V(55)

  B(209) = RCT(117)*V(48)

  B(210) = RCT(118)*V(55)

  B(211) = RCT(118)*V(53)

  B(212) = RCT(119)*V(55)

  B(213) = RCT(119)*V(46)

  B(214) = RCT(120)*V(55)

  B(215) = RCT(120)*V(41)

  B(216) = RCT(121)*V(55)

  B(217) = RCT(121)*V(37)

  B(218) = RCT(122)*V(50)

  B(219) = RCT(122)*V(8)

  B(220) = RCT(123)*V(29)

  B(221) = RCT(123)*V(19)

  B(222) = 1

  B(223) = 1

  B(224) = 1

  B(225) = RCT(127)

  B(226) = RCT(128)

  B(227) = RCT(129)

  B(228) = RCT(130)

  B(229) = RCT(131)

  B(230) = RCT(132)

  B(231) = RCT(133)

  B(232) = RCT(134)

  B(233) = RCT(135)

  B(234) = RCT(136)



  JVS(1) = 0

  JVS(2) = B(85)

  JVS(3) = B(86)

  JVS(4) = 0

  JVS(5) = 0

  JVS(6) = 0

  JVS(7) = 0.52*B(107)

  JVS(8) = 0.22*B(111)

  JVS(9) = 0.52*B(108)+0.22*B(112)

  JVS(10) = 0

  JVS(11) = 0.09*B(111)

  JVS(12) = 0.39*B(141)

  JVS(13) = 0.16*B(113)

  JVS(14) = 0.46*B(147)

  JVS(15) = 0.09*B(112)+0.16*B(114)+0.39*B(142)+0.46*B(148)

  JVS(16) = 0.4*B(210)

  JVS(17) = 0.4*B(211)

  JVS(18) = -B(21)-B(23)

  JVS(19) = B(8)

  JVS(20) = -B(85)

  JVS(21) = -B(86)

  JVS(22) = -B(218)

  JVS(23) = -B(219)

  JVS(24) = -B(9)-B(59)

  JVS(25) = -B(60)

  JVS(26) = B(61)+B(62)

  JVS(27) = -B(173)

  JVS(28) = B(171)

  JVS(29) = B(172)

  JVS(30) = -B(123)

  JVS(31) = -B(124)

  JVS(32) = -B(6)-B(80)-B(82)

  JVS(33) = B(75)

  JVS(34) = B(76)

  JVS(35) = -B(125)

  JVS(36) = -B(126)

  JVS(37) = -B(87)

  JVS(38) = 0.06*B(111)

  JVS(39) = 0.08*B(113)

  JVS(40) = -B(88)

  JVS(41) = 0.06*B(112)+0.08*B(114)

  JVS(42) = -B(133)

  JVS(43) = 0.4*B(129)+B(131)

  JVS(44) = 0.4*B(130)

  JVS(45) = -B(134)

  JVS(46) = B(132)

  JVS(47) = -B(5)-B(55)-B(70)

  JVS(48) = -B(56)

  JVS(49) = B(66)

  JVS(50) = B(67)

  JVS(51) = -B(89)

  JVS(52) = 0.01*B(111)

  JVS(53) = 0.01*B(113)

  JVS(54) = 0.2*B(226)

  JVS(55) = -B(90)

  JVS(56) = 0.01*B(112)+0.01*B(114)

  JVS(57) = 0.8*B(123)

  JVS(58) = 0.45*B(125)

  JVS(59) = -B(127)

  JVS(60) = 0.8*B(124)+0.45*B(126)

  JVS(61) = -B(128)

  JVS(62) = -B(220)

  JVS(63) = -B(221)

  JVS(64) = 1.06*B(111)+B(115)

  JVS(65) = 2.26*B(113)+2.23*B(117)

  JVS(66) = B(186)+B(200)+B(230)

  JVS(67) = 1.98*B(18)+0.42*B(167)

  JVS(68) = 1.68*B(180)+1.98*B(194)+1.25*B(227)

  JVS(69) = 1.98*B(19)

  JVS(70) = B(116)+2.23*B(118)+0.42*B(168)

  JVS(71) = 1.06*B(112)+2.26*B(114)

  JVS(72) = 1.98*B(195)+B(201)

  JVS(73) = 1.68*B(181)+B(187)

  JVS(74) = -B(12)-B(163)

  JVS(75) = B(204)

  JVS(76) = -B(164)

  JVS(77) = B(205)

  JVS(78) = -B(13)-B(165)

  JVS(79) = B(206)

  JVS(80) = -B(166)

  JVS(81) = B(207)

  JVS(82) = -B(3)-B(51)

  JVS(83) = B(45)-B(52)

  JVS(84) = B(68)

  JVS(85) = B(69)

  JVS(86) = B(46)

  JVS(87) = -B(107)-B(109)

  JVS(88) = -B(110)

  JVS(89) = -B(108)

  JVS(90) = -B(174)

  JVS(91) = 0.03*B(111)

  JVS(92) = 0.04*B(113)

  JVS(93) = 0.34*B(225)

  JVS(94) = -B(175)

  JVS(95) = 0.03*B(112)+0.04*B(114)

  JVS(96) = B(21)

  JVS(97) = -B(25)-B(27)-B(29)-B(31)-B(33)

  JVS(98) = B(7)-B(28)

  JVS(99) = B(1)-B(30)-B(32)

  JVS(100) = 0.89*B(2)

  JVS(101) = -B(34)

  JVS(102) = 0.12*B(123)

  JVS(103) = 0.05*B(125)

  JVS(104) = -B(129)-B(131)

  JVS(105) = 0.12*B(124)+0.05*B(126)-B(130)

  JVS(106) = -B(132)

  JVS(107) = 2*B(80)

  JVS(108) = B(131)

  JVS(109) = -B(4)-B(53)

  JVS(110) = B(95)

  JVS(111) = B(105)

  JVS(112) = B(99)

  JVS(113) = 0.07*B(149)

  JVS(114) = B(47)-B(54)

  JVS(115) = B(48)

  JVS(116) = 0.3*B(78)+B(96)+B(100)+B(106)+B(132)+0.07*B(150)

  JVS(117) = 0.3*B(79)

  JVS(118) = 0.24*B(107)

  JVS(119) = -B(83)

  JVS(120) = B(17)+2*B(135)+0.69*B(137)

  JVS(121) = 0.59*B(155)

  JVS(122) = 0.31*B(111)

  JVS(123) = 0.07*B(141)

  JVS(124) = B(10)+B(11)+B(93)+B(95)

  JVS(125) = 0.3*B(113)

  JVS(126) = B(16)+B(105)

  JVS(127) = B(14)

  JVS(128) = 0.33*B(20)+0.16*B(147)+0.64*B(149)

  JVS(129) = -B(84)+B(94)+2*B(136)

  JVS(130) = 0.24*B(108)+0.31*B(112)+0.3*B(114)+0.69*B(138)+0.07*B(142)+0.16*B(148)

  JVS(131) = B(96)+B(106)+0.64*B(150)

  JVS(132) = 0.59*B(156)

  JVS(133) = 1.1*B(125)

  JVS(134) = -B(220)

  JVS(135) = -B(91)-B(221)

  JVS(136) = 1.6*B(153)+2*B(159)+2*B(232)

  JVS(137) = 0.18*B(151)

  JVS(138) = 0

  JVS(139) = 0

  JVS(140) = 0

  JVS(141) = 1.86*B(149)

  JVS(142) = 0

  JVS(143) = 0

  JVS(144) = 0

  JVS(145) = -B(92)+1.1*B(126)

  JVS(146) = 0

  JVS(147) = 1.86*B(150)

  JVS(148) = 2*B(160)

  JVS(149) = 0.18*B(152)+1.6*B(154)

  JVS(150) = 0.95*B(127)

  JVS(151) = 0.3*B(129)

  JVS(152) = -B(17)-B(135)-B(137)

  JVS(153) = 0.3*B(130)-B(136)

  JVS(154) = -B(138)

  JVS(155) = 0

  JVS(156) = 0.95*B(128)

  JVS(157) = -B(153)-B(159)-B(232)

  JVS(158) = B(143)

  JVS(159) = B(144)

  JVS(160) = -B(160)

  JVS(161) = -B(154)

  JVS(162) = -B(151)-B(157)-B(231)

  JVS(163) = B(139)

  JVS(164) = B(140)

  JVS(165) = -B(158)

  JVS(166) = -B(152)

  JVS(167) = -B(155)-B(161)-B(233)

  JVS(168) = 0.5*B(145)

  JVS(169) = 0.5*B(146)

  JVS(170) = -B(162)

  JVS(171) = -B(156)

  JVS(172) = -B(111)-B(115)-B(119)

  JVS(173) = -B(116)

  JVS(174) = -B(112)

  JVS(175) = -B(120)

  JVS(176) = -B(139)-B(141)-B(143)

  JVS(177) = -B(140)

  JVS(178) = -B(142)

  JVS(179) = -B(144)

  JVS(180) = B(12)+0.3*B(163)

  JVS(181) = B(107)+1.56*B(109)

  JVS(182) = B(174)

  JVS(183) = B(135)+0.7*B(137)

  JVS(184) = 0.63*B(151)

  JVS(185) = 0.25*B(155)

  JVS(186) = 0.57*B(111)+B(115)

  JVS(187) = 0.6*B(141)

  JVS(188) = -B(10)-B(11)-B(93)-B(95)

  JVS(189) = 0

  JVS(190) = 0.5*B(186)+0.5*B(200)+0.5*B(230)

  JVS(191) = B(176)+B(190)+0.66*B(225)

  JVS(192) = 0.2*B(20)+0.15*B(147)+0.28*B(149)

  JVS(193) = B(184)+B(198)+0.7*B(229)

  JVS(194) = -B(94)+1.56*B(110)+B(116)+B(136)+0.3*B(164)+B(175)

  JVS(195) = B(108)+0.57*B(112)+0.7*B(138)+0.6*B(142)+0.15*B(148)

  JVS(196) = -B(96)+0.28*B(150)+B(191)+B(199)+0.5*B(201)

  JVS(197) = 0

  JVS(198) = 0.63*B(152)+0.25*B(156)+B(177)+B(185)+0.5*B(187)

  JVS(199) = 0.08*B(123)

  JVS(200) = 0.5*B(125)

  JVS(201) = B(109)

  JVS(202) = 0.6*B(129)

  JVS(203) = B(135)+0.03*B(137)

  JVS(204) = B(115)

  JVS(205) = 0.08*B(139)+0.2*B(141)

  JVS(206) = -B(188)-B(202)-B(216)-B(234)

  JVS(207) = B(117)

  JVS(208) = B(103)

  JVS(209) = 0.2*B(145)+0.07*B(147)+0.93*B(149)

  JVS(210) = 0.4*B(18)

  JVS(211) = 0.34*B(180)+0.4*B(194)+0.24*B(227)

  JVS(212) = 0.41*B(19)

  JVS(213) = B(104)+B(110)+B(116)+B(118)+0.08*B(124)+0.5*B(126)+0.6*B(130)+B(136)+0.08*B(140)+0.2*B(146)

  JVS(214) = 0.03*B(138)+0.2*B(142)+0.07*B(148)

  JVS(215) = 0.93*B(150)+0.4*B(195)-B(203)

  JVS(216) = -B(217)

  JVS(217) = 0.34*B(181)-B(189)

  JVS(218) = 0.63*B(155)+0.5*B(233)

  JVS(219) = -B(15)-B(101)

  JVS(220) = 0.07*B(113)+0.23*B(117)

  JVS(221) = 0.03*B(20)+0.09*B(147)

  JVS(222) = 0.15*B(229)

  JVS(223) = 0.74*B(18)

  JVS(224) = 0.62*B(180)+0.74*B(194)+0.57*B(227)

  JVS(225) = 0.74*B(19)

  JVS(226) = -B(102)+0.23*B(118)

  JVS(227) = 0.07*B(114)+0.09*B(148)

  JVS(228) = 0.74*B(195)

  JVS(229) = 0

  JVS(230) = 0.63*B(156)+0.62*B(181)

  JVS(231) = -B(113)-B(117)-B(121)

  JVS(232) = -B(118)

  JVS(233) = -B(114)

  JVS(234) = -B(122)

  JVS(235) = B(89)

  JVS(236) = 0.7*B(165)

  JVS(237) = 0.06*B(111)

  JVS(238) = 0.05*B(113)

  JVS(239) = -B(178)-B(192)-B(206)-B(226)

  JVS(240) = 0.1*B(18)

  JVS(241) = 0.08*B(180)+0.1*B(194)+0.06*B(227)

  JVS(242) = 0.1*B(19)

  JVS(243) = B(90)+0.7*B(166)

  JVS(244) = 0.06*B(112)+0.05*B(114)

  JVS(245) = -B(193)+0.1*B(195)

  JVS(246) = -B(207)

  JVS(247) = -B(179)+0.08*B(181)

  JVS(248) = B(119)

  JVS(249) = B(121)

  JVS(250) = -B(186)-B(200)-B(214)-B(230)

  JVS(251) = B(169)

  JVS(252) = B(170)

  JVS(253) = 0

  JVS(254) = B(120)+B(122)-B(201)

  JVS(255) = -B(215)

  JVS(256) = -B(187)

  JVS(257) = 0.8*B(125)

  JVS(258) = 0.2*B(137)

  JVS(259) = 0.34*B(155)

  JVS(260) = 0.04*B(111)

  JVS(261) = 0.07*B(113)

  JVS(262) = -B(16)-B(103)-B(105)

  JVS(263) = 0.85*B(147)

  JVS(264) = 0.15*B(229)

  JVS(265) = 0.19*B(167)

  JVS(266) = -B(104)+0.8*B(126)+0.19*B(168)

  JVS(267) = 0.04*B(112)+0.07*B(114)+0.2*B(138)+0.85*B(148)

  JVS(268) = -B(106)

  JVS(269) = 0

  JVS(270) = 0.34*B(156)

  JVS(271) = B(218)

  JVS(272) = B(13)+0.3*B(165)

  JVS(273) = 0.22*B(109)

  JVS(274) = 0.03*B(137)

  JVS(275) = 0.8*B(153)+B(232)

  JVS(276) = 0.55*B(155)+0.5*B(233)

  JVS(277) = 0.47*B(111)+B(115)

  JVS(278) = 0.15*B(141)

  JVS(279) = 1.03*B(113)+1.77*B(117)

  JVS(280) = B(178)+B(192)+0.8*B(226)

  JVS(281) = 0.5*B(186)+0.5*B(200)+0.5*B(230)

  JVS(282) = -B(14)-B(97)-B(99)

  JVS(283) = 0.07*B(20)+0.02*B(147)+0.28*B(149)

  JVS(284) = 0.3*B(18)+0.04*B(167)

  JVS(285) = 0.25*B(180)+0.3*B(194)+0.21*B(227)

  JVS(286) = 0.3*B(19)

  JVS(287) = -B(98)+0.22*B(110)+B(116)+1.77*B(118)+0.3*B(166)+0.04*B(168)+B(219)

  JVS(288) = 0.47*B(112)+1.03*B(114)+0.03*B(138)+0.15*B(142)+0.02*B(148)

  JVS(289) = -B(100)+0.28*B(150)+B(193)+0.3*B(195)+0.5*B(201)

  JVS(290) = 0

  JVS(291) = 0.8*B(154)+0.55*B(156)+B(179)+0.25*B(181)+0.5*B(187)

  JVS(292) = B(87)

  JVS(293) = 0.7*B(163)

  JVS(294) = 0.07*B(111)

  JVS(295) = B(15)

  JVS(296) = 0.1*B(113)

  JVS(297) = B(14)

  JVS(298) = -B(176)-B(190)-B(204)-B(225)

  JVS(299) = 0.7*B(20)+0.05*B(147)

  JVS(300) = 0

  JVS(301) = 0

  JVS(302) = 0

  JVS(303) = 0

  JVS(304) = B(88)+0.7*B(164)

  JVS(305) = 0.07*B(112)+0.1*B(114)+0.05*B(148)

  JVS(306) = B(182)+B(196)+B(228)

  JVS(307) = -B(191)+B(197)

  JVS(308) = -B(205)

  JVS(309) = -B(177)+B(183)

  JVS(310) = 0.2*B(153)

  JVS(311) = 0.91*B(151)+B(231)

  JVS(312) = 0.65*B(141)

  JVS(313) = -B(20)-B(145)-B(147)-B(149)

  JVS(314) = -B(146)

  JVS(315) = 0.65*B(142)-B(148)

  JVS(316) = -B(150)

  JVS(317) = 0

  JVS(318) = 0.91*B(152)+0.2*B(154)

  JVS(319) = B(101)

  JVS(320) = 0.11*B(113)

  JVS(321) = 0

  JVS(322) = -B(184)-B(198)-B(212)-B(229)

  JVS(323) = 0

  JVS(324) = 0

  JVS(325) = 0

  JVS(326) = B(102)

  JVS(327) = 0.11*B(114)

  JVS(328) = -B(199)

  JVS(329) = -B(213)

  JVS(330) = -B(185)

  JVS(331) = B(157)

  JVS(332) = B(161)

  JVS(333) = 0

  JVS(334) = 0

  JVS(335) = B(212)

  JVS(336) = -B(18)-B(167)

  JVS(337) = B(208)

  JVS(338) = 0

  JVS(339) = -B(168)

  JVS(340) = 0

  JVS(341) = 0

  JVS(342) = B(158)+B(162)+B(209)+B(213)

  JVS(343) = 0

  JVS(344) = B(91)

  JVS(345) = 0

  JVS(346) = 0

  JVS(347) = 0.03*B(111)

  JVS(348) = 0

  JVS(349) = 0.09*B(113)

  JVS(350) = 0

  JVS(351) = 0

  JVS(352) = 0.77*B(167)

  JVS(353) = -B(180)-B(194)-B(208)-B(227)

  JVS(354) = 0

  JVS(355) = B(92)+0.77*B(168)

  JVS(356) = 0.03*B(112)+0.09*B(114)

  JVS(357) = -B(195)

  JVS(358) = -B(209)

  JVS(359) = -B(181)

  JVS(360) = B(133)

  JVS(361) = 0.05*B(127)

  JVS(362) = 0

  JVS(363) = 0.8*B(153)+B(159)+B(232)

  JVS(364) = 0.09*B(151)

  JVS(365) = 0

  JVS(366) = 0.5*B(186)+0.5*B(200)+B(214)+0.5*B(230)

  JVS(367) = 0.93*B(149)

  JVS(368) = 0.16*B(180)

  JVS(369) = -B(19)-B(169)

  JVS(370) = -B(170)

  JVS(371) = 0

  JVS(372) = B(134)

  JVS(373) = 0.93*B(150)+0.5*B(201)

  JVS(374) = B(160)+B(215)

  JVS(375) = 0.05*B(128)+0.09*B(152)+0.8*B(154)+0.16*B(181)+0.5*B(187)

  JVS(376) = 2*B(23)

  JVS(377) = -B(85)

  JVS(378) = -B(218)

  JVS(379) = 2*B(9)-B(59)

  JVS(380) = -B(123)

  JVS(381) = -B(125)

  JVS(382) = -B(87)

  JVS(383) = -B(55)

  JVS(384) = -B(89)

  JVS(385) = B(12)-0.7*B(163)

  JVS(386) = B(13)-0.7*B(165)

  JVS(387) = B(3)-B(51)

  JVS(388) = 0.12*B(107)-B(109)

  JVS(389) = -B(174)

  JVS(390) = -B(129)

  JVS(391) = B(4)-B(53)

  JVS(392) = -B(83)

  JVS(393) = -B(91)

  JVS(394) = -B(135)+0.08*B(137)

  JVS(395) = 0

  JVS(396) = 0

  JVS(397) = 0

  JVS(398) = 0.33*B(111)-B(115)

  JVS(399) = -B(139)+0.27*B(141)

  JVS(400) = -B(93)

  JVS(401) = -B(101)

  JVS(402) = 0.6*B(113)-B(117)

  JVS(403) = 0

  JVS(404) = 0

  JVS(405) = -B(103)

  JVS(406) = -B(97)

  JVS(407) = 0

  JVS(408) = -B(145)+0.27*B(147)

  JVS(409) = 0

  JVS(410) = B(18)-0.77*B(167)

  JVS(411) = 0

  JVS(412) = -B(169)

  JVS(413) = -B(39)-B(43)-B(45)-B(47)-B(49)-B(52)-B(54)-B(56)-B(57)-B(60)-B(84)-B(86)-B(88)-B(90)-B(92)-B(94)-B(98)&
               &-B(102)-B(104)-B(110)-B(116)-B(118)-B(124)-B(126)-B(130)-B(136)-B(140)-B(146)-0.7*B(164)-0.7*B(166)-0.77&
               &*B(168)-B(170)-B(175)-B(219)

  JVS(414) = -B(40)+B(41)+0.12*B(108)+0.33*B(112)+0.6*B(114)+0.08*B(138)+0.27*B(142)+0.27*B(148)

  JVS(415) = -B(48)

  JVS(416) = 0

  JVS(417) = -B(50)+0.7*B(78)

  JVS(418) = B(42)-B(58)+B(64)+0.7*B(79)

  JVS(419) = -B(46)+B(65)

  JVS(420) = -B(107)

  JVS(421) = B(25)-B(27)

  JVS(422) = -B(137)

  JVS(423) = -B(111)

  JVS(424) = -B(141)

  JVS(425) = -B(113)

  JVS(426) = -B(147)

  JVS(427) = -B(39)

  JVS(428) = -B(7)-B(8)-B(28)-B(35)-B(37)-B(40)-B(41)-B(108)-B(112)-B(114)-B(138)-B(142)-B(148)

  JVS(429) = -B(38)

  JVS(430) = 0.4*B(210)

  JVS(431) = 0

  JVS(432) = -B(42)+0.4*B(211)

  JVS(433) = -B(36)

  JVS(434) = B(173)

  JVS(435) = B(6)+B(82)

  JVS(436) = -B(133)

  JVS(437) = B(5)+B(55)+B(70)

  JVS(438) = 0.95*B(127)

  JVS(439) = B(51)

  JVS(440) = -B(29)-B(31)+B(33)

  JVS(441) = 0

  JVS(442) = B(4)

  JVS(443) = 1.2*B(153)

  JVS(444) = 0.91*B(151)

  JVS(445) = B(155)

  JVS(446) = 0

  JVS(447) = 0

  JVS(448) = B(188)+B(202)

  JVS(449) = 0

  JVS(450) = B(178)+B(192)

  JVS(451) = 1.5*B(186)+1.5*B(200)+0.5*B(230)

  JVS(452) = 0

  JVS(453) = 0

  JVS(454) = B(176)+B(190)

  JVS(455) = 0

  JVS(456) = B(184)+B(198)

  JVS(457) = 0

  JVS(458) = 0.84*B(180)+B(194)

  JVS(459) = B(19)

  JVS(460) = -B(47)+B(49)+B(52)+B(56)

  JVS(461) = B(35)-B(37)

  JVS(462) = -B(1)-B(30)-B(32)-B(38)-B(48)-B(66)-B(68)-B(75)-B(134)-B(171)

  JVS(463) = -B(172)+B(182)+B(196)

  JVS(464) = 0.89*B(2)+B(50)+2*B(71)-B(76)+2*B(77)+0.7*B(78)+B(191)+B(193)+B(195)+B(197)+B(199)+1.5*B(201)+B(203)

  JVS(465) = B(64)-B(67)-B(69)+0.7*B(79)

  JVS(466) = B(34)+B(36)+B(65)+2*B(72)+0.95*B(128)+0.91*B(152)+1.2*B(154)+B(156)+B(177)+B(179)+0.84*B(181)+B(183)+B(185)&
               &+1.5*B(187)+B(189)

  JVS(467) = B(173)

  JVS(468) = B(17)+B(135)+0.62*B(137)

  JVS(469) = 0.13*B(111)

  JVS(470) = 0.2*B(141)

  JVS(471) = B(15)

  JVS(472) = 0.19*B(113)

  JVS(473) = B(16)+B(103)+B(105)

  JVS(474) = B(97)+B(99)

  JVS(475) = 0.97*B(20)+0.5*B(145)+0.11*B(147)+0.07*B(149)

  JVS(476) = B(184)+B(198)+0.7*B(229)

  JVS(477) = 0

  JVS(478) = 0

  JVS(479) = 0

  JVS(480) = B(98)+B(104)+B(136)+0.5*B(146)

  JVS(481) = 0.13*B(112)+0.19*B(114)+0.62*B(138)+0.2*B(142)+0.11*B(148)

  JVS(482) = -B(171)

  JVS(483) = -B(172)-B(182)-B(196)-B(210)-B(228)

  JVS(484) = B(100)+B(106)+0.07*B(150)-B(197)+B(199)

  JVS(485) = -B(211)

  JVS(486) = -B(183)+B(185)

  JVS(487) = B(6)+B(82)

  JVS(488) = B(31)

  JVS(489) = -B(131)

  JVS(490) = B(53)

  JVS(491) = -B(119)

  JVS(492) = -B(143)

  JVS(493) = -B(95)

  JVS(494) = -B(202)

  JVS(495) = -B(121)

  JVS(496) = -B(192)

  JVS(497) = -B(200)

  JVS(498) = -B(105)

  JVS(499) = -B(99)

  JVS(500) = -B(190)

  JVS(501) = -B(149)

  JVS(502) = -B(198)

  JVS(503) = 0

  JVS(504) = -B(194)

  JVS(505) = 0

  JVS(506) = -B(49)+B(54)

  JVS(507) = B(37)

  JVS(508) = B(32)+B(38)-B(73)-B(75)

  JVS(509) = -B(196)

  JVS(510) = -B(2)-B(50)-B(71)-B(74)-B(76)-2*B(77)-B(78)-B(96)-B(100)-B(106)-B(120)-B(122)-B(132)-B(144)-B(150)-B(191)&
               &-B(193)-B(195)-B(197)-B(199)-B(201)-B(203)

  JVS(511) = -B(79)

  JVS(512) = -B(72)

  JVS(513) = B(85)

  JVS(514) = B(218)

  JVS(515) = B(59)

  JVS(516) = 0.2*B(123)

  JVS(517) = 0.55*B(125)

  JVS(518) = B(5)+B(70)

  JVS(519) = 0.95*B(127)

  JVS(520) = B(12)

  JVS(521) = B(13)

  JVS(522) = 0.22*B(107)+B(109)

  JVS(523) = B(174)

  JVS(524) = 0.6*B(129)

  JVS(525) = B(83)

  JVS(526) = B(17)+2*B(135)+0.76*B(137)

  JVS(527) = 0.8*B(153)-B(159)

  JVS(528) = 0.91*B(151)-B(157)

  JVS(529) = B(155)-B(161)

  JVS(530) = 0.26*B(111)+B(115)

  JVS(531) = 0.07*B(141)

  JVS(532) = 2*B(10)+B(93)+B(95)

  JVS(533) = -B(216)

  JVS(534) = 0.22*B(113)+B(117)

  JVS(535) = B(178)+B(192)-B(206)+0.6*B(226)

  JVS(536) = 0.5*B(186)+0.5*B(200)-B(214)

  JVS(537) = B(16)

  JVS(538) = B(14)

  JVS(539) = B(176)+B(190)-B(204)+0.32*B(225)

  JVS(540) = 0.33*B(20)+0.1*B(147)+0.93*B(149)

  JVS(541) = -B(212)

  JVS(542) = 0.9*B(18)

  JVS(543) = 0.76*B(180)+0.9*B(194)-B(208)+0.54*B(227)

  JVS(544) = 0.9*B(19)

  JVS(545) = B(39)+B(43)+B(49)-B(57)+B(60)+B(84)+B(86)+B(94)+B(110)+B(116)+B(118)+0.2*B(124)+0.55*B(126)+0.6*B(130)+2&
               &*B(136)+B(175)+B(219)

  JVS(546) = B(40)-B(41)+0.22*B(108)+0.26*B(112)+0.22*B(114)+0.76*B(138)+0.07*B(142)+0.1*B(148)

  JVS(547) = -B(66)-B(68)

  JVS(548) = -B(210)

  JVS(549) = B(50)-B(78)+B(96)+0.93*B(150)+B(191)+B(193)+0.9*B(195)+0.5*B(201)

  JVS(550) = -B(42)-B(58)-2*B(61)-2*B(62)-B(64)-B(67)-B(69)-B(79)-B(158)-B(160)-B(162)-B(205)-B(207)-B(209)-B(211)&
               &-B(213)-B(215)-B(217)

  JVS(551) = -B(65)+0.95*B(128)+0.91*B(152)+0.8*B(154)+B(156)+B(177)+B(179)+0.76*B(181)+0.5*B(187)

  JVS(552) = -B(127)

  JVS(553) = B(3)

  JVS(554) = B(29)-B(33)

  JVS(555) = -B(153)

  JVS(556) = -B(151)

  JVS(557) = -B(155)

  JVS(558) = 0

  JVS(559) = -B(188)

  JVS(560) = 0

  JVS(561) = -B(178)

  JVS(562) = -B(186)

  JVS(563) = 0

  JVS(564) = -B(176)

  JVS(565) = 0

  JVS(566) = -B(184)

  JVS(567) = 0

  JVS(568) = -B(180)

  JVS(569) = 0

  JVS(570) = -B(45)

  JVS(571) = -B(35)

  JVS(572) = B(1)+B(30)+B(73)

  JVS(573) = -B(182)

  JVS(574) = 0.11*B(2)-B(71)+B(74)

  JVS(575) = -B(64)

  JVS(576) = -B(34)-B(36)-B(46)-B(65)-B(72)-B(128)-B(152)-B(154)-B(156)-B(177)-B(179)-B(181)-B(183)-B(185)-B(187)-B(189)
      
END SUBROUTINE cbmz_bb_Jac_SP














SUBROUTINE cbmz_bb_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(576), W(56), a
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
      
END SUBROUTINE cbmz_bb_KppDecomp



SUBROUTINE cbmz_bb_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(576), W(56), a
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
      
END SUBROUTINE cbmz_bb_KppDecompCmplx


SUBROUTINE cbmz_bb_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(576), X(56), sum

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
      
END SUBROUTINE cbmz_bb_KppSolveIndirect


SUBROUTINE cbmz_bb_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(576), X(56), sum

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
      
END SUBROUTINE cbmz_bb_KppSolveCmplx













SUBROUTINE cbmz_bb_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(18) = X(18)-JVS(57)*X(11)-JVS(58)*X(13)
  X(25) = X(25)-JVS(96)*X(6)
  X(26) = X(26)-JVS(102)*X(11)-JVS(103)*X(13)
  X(27) = X(27)-JVS(107)*X(12)-JVS(108)*X(26)
  X(28) = X(28)-JVS(118)*X(23)
  X(29) = X(29)-JVS(133)*X(13)-JVS(134)*X(19)
  X(30) = X(30)-JVS(150)*X(18)-JVS(151)*X(26)
  X(36) = X(36)-JVS(180)*X(20)-JVS(181)*X(23)-JVS(182)*X(24)-JVS(183)*X(30)-JVS(184)*X(32)-JVS(185)*X(33)-JVS(186)*X(34)&
            &-JVS(187)*X(35)
  X(37) = X(37)-JVS(199)*X(11)-JVS(200)*X(13)-JVS(201)*X(23)-JVS(202)*X(26)-JVS(203)*X(30)-JVS(204)*X(34)-JVS(205)*X(35)
  X(38) = X(38)-JVS(218)*X(33)
  X(40) = X(40)-JVS(235)*X(17)-JVS(236)*X(21)-JVS(237)*X(34)-JVS(238)*X(39)
  X(41) = X(41)-JVS(248)*X(34)-JVS(249)*X(39)
  X(42) = X(42)-JVS(257)*X(13)-JVS(258)*X(30)-JVS(259)*X(33)-JVS(260)*X(34)-JVS(261)*X(39)
  X(43) = X(43)-JVS(271)*X(8)-JVS(272)*X(21)-JVS(273)*X(23)-JVS(274)*X(30)-JVS(275)*X(31)-JVS(276)*X(33)-JVS(277)*X(34)&
            &-JVS(278)*X(35)-JVS(279)*X(39)-JVS(280)*X(40)-JVS(281)*X(41)
  X(44) = X(44)-JVS(292)*X(14)-JVS(293)*X(20)-JVS(294)*X(34)-JVS(295)*X(38)-JVS(296)*X(39)-JVS(297)*X(43)
  X(45) = X(45)-JVS(310)*X(31)-JVS(311)*X(32)-JVS(312)*X(35)
  X(46) = X(46)-JVS(319)*X(38)-JVS(320)*X(39)-JVS(321)*X(45)
  X(47) = X(47)-JVS(331)*X(32)-JVS(332)*X(33)-JVS(333)*X(35)-JVS(334)*X(45)-JVS(335)*X(46)
  X(48) = X(48)-JVS(344)*X(29)-JVS(345)*X(31)-JVS(346)*X(32)-JVS(347)*X(34)-JVS(348)*X(35)-JVS(349)*X(39)-JVS(350)*X(41)&
            &-JVS(351)*X(45)-JVS(352)*X(47)
  X(49) = X(49)-JVS(360)*X(15)-JVS(361)*X(18)-JVS(362)*X(26)-JVS(363)*X(31)-JVS(364)*X(32)-JVS(365)*X(35)-JVS(366)*X(41)&
            &-JVS(367)*X(45)-JVS(368)*X(48)
  X(50) = X(50)-JVS(376)*X(6)-JVS(377)*X(7)-JVS(378)*X(8)-JVS(379)*X(9)-JVS(380)*X(11)-JVS(381)*X(13)-JVS(382)*X(14)&
            &-JVS(383)*X(16)-JVS(384)*X(17)-JVS(385)*X(20)-JVS(386)*X(21)-JVS(387)*X(22)-JVS(388)*X(23)-JVS(389)*X(24)&
            &-JVS(390)*X(26)-JVS(391)*X(27)-JVS(392)*X(28)-JVS(393)*X(29)-JVS(394)*X(30)-JVS(395)*X(31)-JVS(396)*X(32)&
            &-JVS(397)*X(33)-JVS(398)*X(34)-JVS(399)*X(35)-JVS(400)*X(36)-JVS(401)*X(38)-JVS(402)*X(39)-JVS(403)*X(40)&
            &-JVS(404)*X(41)-JVS(405)*X(42)-JVS(406)*X(43)-JVS(407)*X(44)-JVS(408)*X(45)-JVS(409)*X(46)-JVS(410)*X(47)&
            &-JVS(411)*X(48)-JVS(412)*X(49)
  X(51) = X(51)-JVS(420)*X(23)-JVS(421)*X(25)-JVS(422)*X(30)-JVS(423)*X(34)-JVS(424)*X(35)-JVS(425)*X(39)-JVS(426)*X(45)&
            &-JVS(427)*X(50)
  X(52) = X(52)-JVS(434)*X(10)-JVS(435)*X(12)-JVS(436)*X(15)-JVS(437)*X(16)-JVS(438)*X(18)-JVS(439)*X(22)-JVS(440)*X(25)&
            &-JVS(441)*X(26)-JVS(442)*X(27)-JVS(443)*X(31)-JVS(444)*X(32)-JVS(445)*X(33)-JVS(446)*X(35)-JVS(447)*X(36)&
            &-JVS(448)*X(37)-JVS(449)*X(39)-JVS(450)*X(40)-JVS(451)*X(41)-JVS(452)*X(42)-JVS(453)*X(43)-JVS(454)*X(44)&
            &-JVS(455)*X(45)-JVS(456)*X(46)-JVS(457)*X(47)-JVS(458)*X(48)-JVS(459)*X(49)-JVS(460)*X(50)-JVS(461)*X(51)
  X(53) = X(53)-JVS(467)*X(10)-JVS(468)*X(30)-JVS(469)*X(34)-JVS(470)*X(35)-JVS(471)*X(38)-JVS(472)*X(39)-JVS(473)*X(42)&
            &-JVS(474)*X(43)-JVS(475)*X(45)-JVS(476)*X(46)-JVS(477)*X(47)-JVS(478)*X(48)-JVS(479)*X(49)-JVS(480)*X(50)&
            &-JVS(481)*X(51)-JVS(482)*X(52)
  X(54) = X(54)-JVS(487)*X(12)-JVS(488)*X(25)-JVS(489)*X(26)-JVS(490)*X(27)-JVS(491)*X(34)-JVS(492)*X(35)-JVS(493)*X(36)&
            &-JVS(494)*X(37)-JVS(495)*X(39)-JVS(496)*X(40)-JVS(497)*X(41)-JVS(498)*X(42)-JVS(499)*X(43)-JVS(500)*X(44)&
            &-JVS(501)*X(45)-JVS(502)*X(46)-JVS(503)*X(47)-JVS(504)*X(48)-JVS(505)*X(49)-JVS(506)*X(50)-JVS(507)*X(51)&
            &-JVS(508)*X(52)-JVS(509)*X(53)
  X(55) = X(55)-JVS(513)*X(7)-JVS(514)*X(8)-JVS(515)*X(9)-JVS(516)*X(11)-JVS(517)*X(13)-JVS(518)*X(16)-JVS(519)*X(18)&
            &-JVS(520)*X(20)-JVS(521)*X(21)-JVS(522)*X(23)-JVS(523)*X(24)-JVS(524)*X(26)-JVS(525)*X(28)-JVS(526)*X(30)&
            &-JVS(527)*X(31)-JVS(528)*X(32)-JVS(529)*X(33)-JVS(530)*X(34)-JVS(531)*X(35)-JVS(532)*X(36)-JVS(533)*X(37)&
            &-JVS(534)*X(39)-JVS(535)*X(40)-JVS(536)*X(41)-JVS(537)*X(42)-JVS(538)*X(43)-JVS(539)*X(44)-JVS(540)*X(45)&
            &-JVS(541)*X(46)-JVS(542)*X(47)-JVS(543)*X(48)-JVS(544)*X(49)-JVS(545)*X(50)-JVS(546)*X(51)-JVS(547)*X(52)&
            &-JVS(548)*X(53)-JVS(549)*X(54)
  X(56) = X(56)-JVS(552)*X(18)-JVS(553)*X(22)-JVS(554)*X(25)-JVS(555)*X(31)-JVS(556)*X(32)-JVS(557)*X(33)-JVS(558)*X(35)&
            &-JVS(559)*X(37)-JVS(560)*X(39)-JVS(561)*X(40)-JVS(562)*X(41)-JVS(563)*X(42)-JVS(564)*X(44)-JVS(565)*X(45)&
            &-JVS(566)*X(46)-JVS(567)*X(47)-JVS(568)*X(48)-JVS(569)*X(49)-JVS(570)*X(50)-JVS(571)*X(51)-JVS(572)*X(52)&
            &-JVS(573)*X(53)-JVS(574)*X(54)-JVS(575)*X(55)
  X(56) = X(56)/JVS(576)
  X(55) = (X(55)-JVS(551)*X(56))/(JVS(550))
  X(54) = (X(54)-JVS(511)*X(55)-JVS(512)*X(56))/(JVS(510))
  X(53) = (X(53)-JVS(484)*X(54)-JVS(485)*X(55)-JVS(486)*X(56))/(JVS(483))
  X(52) = (X(52)-JVS(463)*X(53)-JVS(464)*X(54)-JVS(465)*X(55)-JVS(466)*X(56))/(JVS(462))
  X(51) = (X(51)-JVS(429)*X(52)-JVS(430)*X(53)-JVS(431)*X(54)-JVS(432)*X(55)-JVS(433)*X(56))/(JVS(428))
  X(50) = (X(50)-JVS(414)*X(51)-JVS(415)*X(52)-JVS(416)*X(53)-JVS(417)*X(54)-JVS(418)*X(55)-JVS(419)*X(56))/(JVS(413))
  X(49) = (X(49)-JVS(370)*X(50)-JVS(371)*X(51)-JVS(372)*X(52)-JVS(373)*X(54)-JVS(374)*X(55)-JVS(375)*X(56))/(JVS(369))
  X(48) = (X(48)-JVS(354)*X(49)-JVS(355)*X(50)-JVS(356)*X(51)-JVS(357)*X(54)-JVS(358)*X(55)-JVS(359)*X(56))/(JVS(353))
  X(47) = (X(47)-JVS(337)*X(48)-JVS(338)*X(49)-JVS(339)*X(50)-JVS(340)*X(51)-JVS(341)*X(54)-JVS(342)*X(55)-JVS(343)&
            &*X(56))/(JVS(336))
  X(46) = (X(46)-JVS(323)*X(47)-JVS(324)*X(48)-JVS(325)*X(49)-JVS(326)*X(50)-JVS(327)*X(51)-JVS(328)*X(54)-JVS(329)&
            &*X(55)-JVS(330)*X(56))/(JVS(322))
  X(45) = (X(45)-JVS(314)*X(50)-JVS(315)*X(51)-JVS(316)*X(54)-JVS(317)*X(55)-JVS(318)*X(56))/(JVS(313))
  X(44) = (X(44)-JVS(299)*X(45)-JVS(300)*X(46)-JVS(301)*X(47)-JVS(302)*X(48)-JVS(303)*X(49)-JVS(304)*X(50)-JVS(305)&
            &*X(51)-JVS(306)*X(53)-JVS(307)*X(54)-JVS(308)*X(55)-JVS(309)*X(56))/(JVS(298))
  X(43) = (X(43)-JVS(283)*X(45)-JVS(284)*X(47)-JVS(285)*X(48)-JVS(286)*X(49)-JVS(287)*X(50)-JVS(288)*X(51)-JVS(289)&
            &*X(54)-JVS(290)*X(55)-JVS(291)*X(56))/(JVS(282))
  X(42) = (X(42)-JVS(263)*X(45)-JVS(264)*X(46)-JVS(265)*X(47)-JVS(266)*X(50)-JVS(267)*X(51)-JVS(268)*X(54)-JVS(269)&
            &*X(55)-JVS(270)*X(56))/(JVS(262))
  X(41) = (X(41)-JVS(251)*X(49)-JVS(252)*X(50)-JVS(253)*X(51)-JVS(254)*X(54)-JVS(255)*X(55)-JVS(256)*X(56))/(JVS(250))
  X(40) = (X(40)-JVS(240)*X(47)-JVS(241)*X(48)-JVS(242)*X(49)-JVS(243)*X(50)-JVS(244)*X(51)-JVS(245)*X(54)-JVS(246)&
            &*X(55)-JVS(247)*X(56))/(JVS(239))
  X(39) = (X(39)-JVS(232)*X(50)-JVS(233)*X(51)-JVS(234)*X(54))/(JVS(231))
  X(38) = (X(38)-JVS(220)*X(39)-JVS(221)*X(45)-JVS(222)*X(46)-JVS(223)*X(47)-JVS(224)*X(48)-JVS(225)*X(49)-JVS(226)&
            &*X(50)-JVS(227)*X(51)-JVS(228)*X(54)-JVS(229)*X(55)-JVS(230)*X(56))/(JVS(219))
  X(37) = (X(37)-JVS(207)*X(39)-JVS(208)*X(42)-JVS(209)*X(45)-JVS(210)*X(47)-JVS(211)*X(48)-JVS(212)*X(49)-JVS(213)&
            &*X(50)-JVS(214)*X(51)-JVS(215)*X(54)-JVS(216)*X(55)-JVS(217)*X(56))/(JVS(206))
  X(36) = (X(36)-JVS(189)*X(39)-JVS(190)*X(41)-JVS(191)*X(44)-JVS(192)*X(45)-JVS(193)*X(46)-JVS(194)*X(50)-JVS(195)&
            &*X(51)-JVS(196)*X(54)-JVS(197)*X(55)-JVS(198)*X(56))/(JVS(188))
  X(35) = (X(35)-JVS(177)*X(50)-JVS(178)*X(51)-JVS(179)*X(54))/(JVS(176))
  X(34) = (X(34)-JVS(173)*X(50)-JVS(174)*X(51)-JVS(175)*X(54))/(JVS(172))
  X(33) = (X(33)-JVS(168)*X(45)-JVS(169)*X(50)-JVS(170)*X(55)-JVS(171)*X(56))/(JVS(167))
  X(32) = (X(32)-JVS(163)*X(35)-JVS(164)*X(50)-JVS(165)*X(55)-JVS(166)*X(56))/(JVS(162))
  X(31) = (X(31)-JVS(158)*X(35)-JVS(159)*X(54)-JVS(160)*X(55)-JVS(161)*X(56))/(JVS(157))
  X(30) = (X(30)-JVS(153)*X(50)-JVS(154)*X(51)-JVS(155)*X(54)-JVS(156)*X(56))/(JVS(152))
  X(29) = (X(29)-JVS(136)*X(31)-JVS(137)*X(32)-JVS(138)*X(34)-JVS(139)*X(39)-JVS(140)*X(41)-JVS(141)*X(45)-JVS(142)&
            &*X(47)-JVS(143)*X(48)-JVS(144)*X(49)-JVS(145)*X(50)-JVS(146)*X(51)-JVS(147)*X(54)-JVS(148)*X(55)-JVS(149)&
            &*X(56))/(JVS(135))
  X(28) = (X(28)-JVS(120)*X(30)-JVS(121)*X(33)-JVS(122)*X(34)-JVS(123)*X(35)-JVS(124)*X(36)-JVS(125)*X(39)-JVS(126)&
            &*X(42)-JVS(127)*X(43)-JVS(128)*X(45)-JVS(129)*X(50)-JVS(130)*X(51)-JVS(131)*X(54)-JVS(132)*X(56))/(JVS(119))
  X(27) = (X(27)-JVS(110)*X(36)-JVS(111)*X(42)-JVS(112)*X(43)-JVS(113)*X(45)-JVS(114)*X(50)-JVS(115)*X(52)-JVS(116)&
            &*X(54)-JVS(117)*X(55))/(JVS(109))
  X(26) = (X(26)-JVS(105)*X(50)-JVS(106)*X(54))/(JVS(104))
  X(25) = (X(25)-JVS(98)*X(51)-JVS(99)*X(52)-JVS(100)*X(54)-JVS(101)*X(56))/(JVS(97))
  X(24) = (X(24)-JVS(91)*X(34)-JVS(92)*X(39)-JVS(93)*X(44)-JVS(94)*X(50)-JVS(95)*X(51))/(JVS(90))
  X(23) = (X(23)-JVS(88)*X(50)-JVS(89)*X(51))/(JVS(87))
  X(22) = (X(22)-JVS(83)*X(50)-JVS(84)*X(52)-JVS(85)*X(55)-JVS(86)*X(56))/(JVS(82))
  X(21) = (X(21)-JVS(79)*X(40)-JVS(80)*X(50)-JVS(81)*X(55))/(JVS(78))
  X(20) = (X(20)-JVS(75)*X(44)-JVS(76)*X(50)-JVS(77)*X(55))/(JVS(74))
  X(19) = (X(19)-JVS(63)*X(29)-JVS(64)*X(34)-JVS(65)*X(39)-JVS(66)*X(41)-JVS(67)*X(47)-JVS(68)*X(48)-JVS(69)*X(49)&
            &-JVS(70)*X(50)-JVS(71)*X(51)-JVS(72)*X(54)-JVS(73)*X(56))/(JVS(62))
  X(18) = (X(18)-JVS(60)*X(50)-JVS(61)*X(56))/(JVS(59))
  X(17) = (X(17)-JVS(52)*X(34)-JVS(53)*X(39)-JVS(54)*X(40)-JVS(55)*X(50)-JVS(56)*X(51))/(JVS(51))
  X(16) = (X(16)-JVS(48)*X(50)-JVS(49)*X(52)-JVS(50)*X(55))/(JVS(47))
  X(15) = (X(15)-JVS(43)*X(26)-JVS(44)*X(50)-JVS(45)*X(52)-JVS(46)*X(54))/(JVS(42))
  X(14) = (X(14)-JVS(38)*X(34)-JVS(39)*X(39)-JVS(40)*X(50)-JVS(41)*X(51))/(JVS(37))
  X(13) = (X(13)-JVS(36)*X(50))/(JVS(35))
  X(12) = (X(12)-JVS(33)*X(52)-JVS(34)*X(54))/(JVS(32))
  X(11) = (X(11)-JVS(31)*X(50))/(JVS(30))
  X(10) = (X(10)-JVS(28)*X(52)-JVS(29)*X(53))/(JVS(27))
  X(9) = (X(9)-JVS(25)*X(50)-JVS(26)*X(55))/(JVS(24))
  X(8) = (X(8)-JVS(23)*X(50))/(JVS(22))
  X(7) = (X(7)-JVS(21)*X(50))/(JVS(20))
  X(6) = (X(6)-JVS(19)*X(51))/(JVS(18))
  X(5) = (X(5)-JVS(11)*X(34)-JVS(12)*X(35)-JVS(13)*X(39)-JVS(14)*X(45)-JVS(15)*X(51)-JVS(16)*X(53)-JVS(17)*X(55))&
           &/(JVS(10))
  X(4) = (X(4)-JVS(7)*X(23)-JVS(8)*X(34)-JVS(9)*X(51))/(JVS(6))
  X(3) = X(3)/JVS(5)
  X(2) = X(2)/JVS(4)
  X(1) = (X(1)-JVS(2)*X(7)-JVS(3)*X(50))/(JVS(1))
      
END SUBROUTINE cbmz_bb_KppSolve
























      SUBROUTINE cbmz_bb_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE cbmz_bb_WCOPY



      SUBROUTINE cbmz_bb_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE cbmz_bb_WAXPY




      SUBROUTINE cbmz_bb_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE cbmz_bb_WSCAL


      REAL(kind=dp) FUNCTION cbmz_bb_WLAMCH( C )








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
          CALL cbmz_bb_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      cbmz_bb_WLAMCH = Eps

      END FUNCTION cbmz_bb_WLAMCH
     
      SUBROUTINE cbmz_bb_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE cbmz_bb_WLAMCH_ADD




      SUBROUTINE cbmz_bb_SET2ZERO(N,Y)




      
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

      END SUBROUTINE cbmz_bb_SET2ZERO



      REAL(kind=dp) FUNCTION cbmz_bb_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      cbmz_bb_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        cbmz_bb_WDOT = cbmz_bb_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         cbmz_bb_WDOT = cbmz_bb_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          cbmz_bb_WDOT = cbmz_bb_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        cbmz_bb_WDOT = cbmz_bb_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION cbmz_bb_WDOT                                          




   SUBROUTINE decomp_cbmz_bb( JVS, IER )
   
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
   W( 7 ) = JVS( 2 )
   W( 50 ) = JVS( 3 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 7 )
  JVS( 3) = W( 50 )
  IF ( ABS(  JVS( 4 )) < TINY(a) ) THEN
         IER = 2                                       
         RETURN
  END IF
   W( 2 ) = JVS( 4 )
  JVS( 4) = W( 2 )
  IF ( ABS(  JVS( 5 )) < TINY(a) ) THEN
         IER = 3                                       
         RETURN
  END IF
   W( 3 ) = JVS( 5 )
  JVS( 5) = W( 3 )
  IF ( ABS(  JVS( 6 )) < TINY(a) ) THEN
         IER = 4                                       
         RETURN
  END IF
   W( 4 ) = JVS( 6 )
   W( 23 ) = JVS( 7 )
   W( 34 ) = JVS( 8 )
   W( 51 ) = JVS( 9 )
  JVS( 6) = W( 4 )
  JVS( 7) = W( 23 )
  JVS( 8) = W( 34 )
  JVS( 9) = W( 51 )
  IF ( ABS(  JVS( 10 )) < TINY(a) ) THEN
         IER = 5                                       
         RETURN
  END IF
   W( 5 ) = JVS( 10 )
   W( 34 ) = JVS( 11 )
   W( 35 ) = JVS( 12 )
   W( 39 ) = JVS( 13 )
   W( 45 ) = JVS( 14 )
   W( 51 ) = JVS( 15 )
   W( 53 ) = JVS( 16 )
   W( 55 ) = JVS( 17 )
  JVS( 10) = W( 5 )
  JVS( 11) = W( 34 )
  JVS( 12) = W( 35 )
  JVS( 13) = W( 39 )
  JVS( 14) = W( 45 )
  JVS( 15) = W( 51 )
  JVS( 16) = W( 53 )
  JVS( 17) = W( 55 )
  IF ( ABS(  JVS( 18 )) < TINY(a) ) THEN
         IER = 6                                       
         RETURN
  END IF
   W( 6 ) = JVS( 18 )
   W( 51 ) = JVS( 19 )
  JVS( 18) = W( 6 )
  JVS( 19) = W( 51 )
  IF ( ABS(  JVS( 20 )) < TINY(a) ) THEN
         IER = 7                                       
         RETURN
  END IF
   W( 7 ) = JVS( 20 )
   W( 50 ) = JVS( 21 )
  JVS( 20) = W( 7 )
  JVS( 21) = W( 50 )
  IF ( ABS(  JVS( 22 )) < TINY(a) ) THEN
         IER = 8                                       
         RETURN
  END IF
   W( 8 ) = JVS( 22 )
   W( 50 ) = JVS( 23 )
  JVS( 22) = W( 8 )
  JVS( 23) = W( 50 )
  IF ( ABS(  JVS( 24 )) < TINY(a) ) THEN
         IER = 9                                       
         RETURN
  END IF
   W( 9 ) = JVS( 24 )
   W( 50 ) = JVS( 25 )
   W( 55 ) = JVS( 26 )
  JVS( 24) = W( 9 )
  JVS( 25) = W( 50 )
  JVS( 26) = W( 55 )
  IF ( ABS(  JVS( 27 )) < TINY(a) ) THEN
         IER = 10                                      
         RETURN
  END IF
   W( 10 ) = JVS( 27 )
   W( 52 ) = JVS( 28 )
   W( 53 ) = JVS( 29 )
  JVS( 27) = W( 10 )
  JVS( 28) = W( 52 )
  JVS( 29) = W( 53 )
  IF ( ABS(  JVS( 30 )) < TINY(a) ) THEN
         IER = 11                                      
         RETURN
  END IF
   W( 11 ) = JVS( 30 )
   W( 50 ) = JVS( 31 )
  JVS( 30) = W( 11 )
  JVS( 31) = W( 50 )
  IF ( ABS(  JVS( 32 )) < TINY(a) ) THEN
         IER = 12                                      
         RETURN
  END IF
   W( 12 ) = JVS( 32 )
   W( 52 ) = JVS( 33 )
   W( 54 ) = JVS( 34 )
  JVS( 32) = W( 12 )
  JVS( 33) = W( 52 )
  JVS( 34) = W( 54 )
  IF ( ABS(  JVS( 35 )) < TINY(a) ) THEN
         IER = 13                                      
         RETURN
  END IF
   W( 13 ) = JVS( 35 )
   W( 50 ) = JVS( 36 )
  JVS( 35) = W( 13 )
  JVS( 36) = W( 50 )
  IF ( ABS(  JVS( 37 )) < TINY(a) ) THEN
         IER = 14                                      
         RETURN
  END IF
   W( 14 ) = JVS( 37 )
   W( 34 ) = JVS( 38 )
   W( 39 ) = JVS( 39 )
   W( 50 ) = JVS( 40 )
   W( 51 ) = JVS( 41 )
  JVS( 37) = W( 14 )
  JVS( 38) = W( 34 )
  JVS( 39) = W( 39 )
  JVS( 40) = W( 50 )
  JVS( 41) = W( 51 )
  IF ( ABS(  JVS( 42 )) < TINY(a) ) THEN
         IER = 15                                      
         RETURN
  END IF
   W( 15 ) = JVS( 42 )
   W( 26 ) = JVS( 43 )
   W( 50 ) = JVS( 44 )
   W( 52 ) = JVS( 45 )
   W( 54 ) = JVS( 46 )
  JVS( 42) = W( 15 )
  JVS( 43) = W( 26 )
  JVS( 44) = W( 50 )
  JVS( 45) = W( 52 )
  JVS( 46) = W( 54 )
  IF ( ABS(  JVS( 47 )) < TINY(a) ) THEN
         IER = 16                                      
         RETURN
  END IF
   W( 16 ) = JVS( 47 )
   W( 50 ) = JVS( 48 )
   W( 52 ) = JVS( 49 )
   W( 55 ) = JVS( 50 )
  JVS( 47) = W( 16 )
  JVS( 48) = W( 50 )
  JVS( 49) = W( 52 )
  JVS( 50) = W( 55 )
  IF ( ABS(  JVS( 51 )) < TINY(a) ) THEN
         IER = 17                                      
         RETURN
  END IF
   W( 17 ) = JVS( 51 )
   W( 34 ) = JVS( 52 )
   W( 39 ) = JVS( 53 )
   W( 40 ) = JVS( 54 )
   W( 50 ) = JVS( 55 )
   W( 51 ) = JVS( 56 )
  JVS( 51) = W( 17 )
  JVS( 52) = W( 34 )
  JVS( 53) = W( 39 )
  JVS( 54) = W( 40 )
  JVS( 55) = W( 50 )
  JVS( 56) = W( 51 )
  IF ( ABS(  JVS( 59 )) < TINY(a) ) THEN
         IER = 18                                      
         RETURN
  END IF
   W( 11 ) = JVS( 57 )
   W( 13 ) = JVS( 58 )
   W( 18 ) = JVS( 59 )
   W( 50 ) = JVS( 60 )
   W( 56 ) = JVS( 61 )
  a = -W( 11 ) / JVS(           30  )
  W( 11 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 31 )
  a = -W( 13 ) / JVS(           35  )
  W( 13 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 36 )
  JVS( 57) = W( 11 )
  JVS( 58) = W( 13 )
  JVS( 59) = W( 18 )
  JVS( 60) = W( 50 )
  JVS( 61) = W( 56 )
  IF ( ABS(  JVS( 62 )) < TINY(a) ) THEN
         IER = 19                                      
         RETURN
  END IF
   W( 19 ) = JVS( 62 )
   W( 29 ) = JVS( 63 )
   W( 34 ) = JVS( 64 )
   W( 39 ) = JVS( 65 )
   W( 41 ) = JVS( 66 )
   W( 47 ) = JVS( 67 )
   W( 48 ) = JVS( 68 )
   W( 49 ) = JVS( 69 )
   W( 50 ) = JVS( 70 )
   W( 51 ) = JVS( 71 )
   W( 54 ) = JVS( 72 )
   W( 56 ) = JVS( 73 )
  JVS( 62) = W( 19 )
  JVS( 63) = W( 29 )
  JVS( 64) = W( 34 )
  JVS( 65) = W( 39 )
  JVS( 66) = W( 41 )
  JVS( 67) = W( 47 )
  JVS( 68) = W( 48 )
  JVS( 69) = W( 49 )
  JVS( 70) = W( 50 )
  JVS( 71) = W( 51 )
  JVS( 72) = W( 54 )
  JVS( 73) = W( 56 )
  IF ( ABS(  JVS( 74 )) < TINY(a) ) THEN
         IER = 20                                      
         RETURN
  END IF
   W( 20 ) = JVS( 74 )
   W( 44 ) = JVS( 75 )
   W( 50 ) = JVS( 76 )
   W( 55 ) = JVS( 77 )
  JVS( 74) = W( 20 )
  JVS( 75) = W( 44 )
  JVS( 76) = W( 50 )
  JVS( 77) = W( 55 )
  IF ( ABS(  JVS( 78 )) < TINY(a) ) THEN
         IER = 21                                      
         RETURN
  END IF
   W( 21 ) = JVS( 78 )
   W( 40 ) = JVS( 79 )
   W( 50 ) = JVS( 80 )
   W( 55 ) = JVS( 81 )
  JVS( 78) = W( 21 )
  JVS( 79) = W( 40 )
  JVS( 80) = W( 50 )
  JVS( 81) = W( 55 )
  IF ( ABS(  JVS( 82 )) < TINY(a) ) THEN
         IER = 22                                      
         RETURN
  END IF
   W( 22 ) = JVS( 82 )
   W( 50 ) = JVS( 83 )
   W( 52 ) = JVS( 84 )
   W( 55 ) = JVS( 85 )
   W( 56 ) = JVS( 86 )
  JVS( 82) = W( 22 )
  JVS( 83) = W( 50 )
  JVS( 84) = W( 52 )
  JVS( 85) = W( 55 )
  JVS( 86) = W( 56 )
  IF ( ABS(  JVS( 87 )) < TINY(a) ) THEN
         IER = 23                                      
         RETURN
  END IF
   W( 23 ) = JVS( 87 )
   W( 50 ) = JVS( 88 )
   W( 51 ) = JVS( 89 )
  JVS( 87) = W( 23 )
  JVS( 88) = W( 50 )
  JVS( 89) = W( 51 )
  IF ( ABS(  JVS( 90 )) < TINY(a) ) THEN
         IER = 24                                      
         RETURN
  END IF
   W( 24 ) = JVS( 90 )
   W( 34 ) = JVS( 91 )
   W( 39 ) = JVS( 92 )
   W( 44 ) = JVS( 93 )
   W( 50 ) = JVS( 94 )
   W( 51 ) = JVS( 95 )
  JVS( 90) = W( 24 )
  JVS( 91) = W( 34 )
  JVS( 92) = W( 39 )
  JVS( 93) = W( 44 )
  JVS( 94) = W( 50 )
  JVS( 95) = W( 51 )
  IF ( ABS(  JVS( 97 )) < TINY(a) ) THEN
         IER = 25                                      
         RETURN
  END IF
   W( 6 ) = JVS( 96 )
   W( 25 ) = JVS( 97 )
   W( 51 ) = JVS( 98 )
   W( 52 ) = JVS( 99 )
   W( 54 ) = JVS( 100 )
   W( 56 ) = JVS( 101 )
  a = -W( 6 ) / JVS(           18  )
  W( 6 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 19 )
  JVS( 96) = W( 6 )
  JVS( 97) = W( 25 )
  JVS( 98) = W( 51 )
  JVS( 99) = W( 52 )
  JVS( 100) = W( 54 )
  JVS( 101) = W( 56 )
  IF ( ABS(  JVS( 104 )) < TINY(a) ) THEN
         IER = 26                                      
         RETURN
  END IF
   W( 11 ) = JVS( 102 )
   W( 13 ) = JVS( 103 )
   W( 26 ) = JVS( 104 )
   W( 50 ) = JVS( 105 )
   W( 54 ) = JVS( 106 )
  a = -W( 11 ) / JVS(           30  )
  W( 11 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 31 )
  a = -W( 13 ) / JVS(           35  )
  W( 13 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 36 )
  JVS( 102) = W( 11 )
  JVS( 103) = W( 13 )
  JVS( 104) = W( 26 )
  JVS( 105) = W( 50 )
  JVS( 106) = W( 54 )
  IF ( ABS(  JVS( 109 )) < TINY(a) ) THEN
         IER = 27                                      
         RETURN
  END IF
   W( 12 ) = JVS( 107 )
   W( 26 ) = JVS( 108 )
   W( 27 ) = JVS( 109 )
   W( 36 ) = JVS( 110 )
   W( 42 ) = JVS( 111 )
   W( 43 ) = JVS( 112 )
   W( 45 ) = JVS( 113 )
   W( 50 ) = JVS( 114 )
   W( 52 ) = JVS( 115 )
   W( 54 ) = JVS( 116 )
   W( 55 ) = JVS( 117 )
  a = -W( 12 ) / JVS(           32  )
  W( 12 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 33 )
  W( 54 ) = W( 54 ) + a*JVS( 34 )
  a = -W( 26 ) / JVS(          104  )
  W( 26 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 105 )
  W( 54 ) = W( 54 ) + a*JVS( 106 )
  JVS( 107) = W( 12 )
  JVS( 108) = W( 26 )
  JVS( 109) = W( 27 )
  JVS( 110) = W( 36 )
  JVS( 111) = W( 42 )
  JVS( 112) = W( 43 )
  JVS( 113) = W( 45 )
  JVS( 114) = W( 50 )
  JVS( 115) = W( 52 )
  JVS( 116) = W( 54 )
  JVS( 117) = W( 55 )
  IF ( ABS(  JVS( 119 )) < TINY(a) ) THEN
         IER = 28                                      
         RETURN
  END IF
   W( 23 ) = JVS( 118 )
   W( 28 ) = JVS( 119 )
   W( 30 ) = JVS( 120 )
   W( 33 ) = JVS( 121 )
   W( 34 ) = JVS( 122 )
   W( 35 ) = JVS( 123 )
   W( 36 ) = JVS( 124 )
   W( 39 ) = JVS( 125 )
   W( 42 ) = JVS( 126 )
   W( 43 ) = JVS( 127 )
   W( 45 ) = JVS( 128 )
   W( 50 ) = JVS( 129 )
   W( 51 ) = JVS( 130 )
   W( 54 ) = JVS( 131 )
   W( 56 ) = JVS( 132 )
  a = -W( 23 ) / JVS(           87  )
  W( 23 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 88 )
  W( 51 ) = W( 51 ) + a*JVS( 89 )
  JVS( 118) = W( 23 )
  JVS( 119) = W( 28 )
  JVS( 120) = W( 30 )
  JVS( 121) = W( 33 )
  JVS( 122) = W( 34 )
  JVS( 123) = W( 35 )
  JVS( 124) = W( 36 )
  JVS( 125) = W( 39 )
  JVS( 126) = W( 42 )
  JVS( 127) = W( 43 )
  JVS( 128) = W( 45 )
  JVS( 129) = W( 50 )
  JVS( 130) = W( 51 )
  JVS( 131) = W( 54 )
  JVS( 132) = W( 56 )
  IF ( ABS(  JVS( 135 )) < TINY(a) ) THEN
         IER = 29                                      
         RETURN
  END IF
   W( 13 ) = JVS( 133 )
   W( 19 ) = JVS( 134 )
   W( 29 ) = JVS( 135 )
   W( 31 ) = JVS( 136 )
   W( 32 ) = JVS( 137 )
   W( 34 ) = JVS( 138 )
   W( 39 ) = JVS( 139 )
   W( 41 ) = JVS( 140 )
   W( 45 ) = JVS( 141 )
   W( 47 ) = JVS( 142 )
   W( 48 ) = JVS( 143 )
   W( 49 ) = JVS( 144 )
   W( 50 ) = JVS( 145 )
   W( 51 ) = JVS( 146 )
   W( 54 ) = JVS( 147 )
   W( 55 ) = JVS( 148 )
   W( 56 ) = JVS( 149 )
  a = -W( 13 ) / JVS(           35  )
  W( 13 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 36 )
  a = -W( 19 ) / JVS(           62  )
  W( 19 ) = -a
  W( 29 ) = W( 29 ) + a*JVS( 63 )
  W( 34 ) = W( 34 ) + a*JVS( 64 )
  W( 39 ) = W( 39 ) + a*JVS( 65 )
  W( 41 ) = W( 41 ) + a*JVS( 66 )
  W( 47 ) = W( 47 ) + a*JVS( 67 )
  W( 48 ) = W( 48 ) + a*JVS( 68 )
  W( 49 ) = W( 49 ) + a*JVS( 69 )
  W( 50 ) = W( 50 ) + a*JVS( 70 )
  W( 51 ) = W( 51 ) + a*JVS( 71 )
  W( 54 ) = W( 54 ) + a*JVS( 72 )
  W( 56 ) = W( 56 ) + a*JVS( 73 )
  JVS( 133) = W( 13 )
  JVS( 134) = W( 19 )
  JVS( 135) = W( 29 )
  JVS( 136) = W( 31 )
  JVS( 137) = W( 32 )
  JVS( 138) = W( 34 )
  JVS( 139) = W( 39 )
  JVS( 140) = W( 41 )
  JVS( 141) = W( 45 )
  JVS( 142) = W( 47 )
  JVS( 143) = W( 48 )
  JVS( 144) = W( 49 )
  JVS( 145) = W( 50 )
  JVS( 146) = W( 51 )
  JVS( 147) = W( 54 )
  JVS( 148) = W( 55 )
  JVS( 149) = W( 56 )
  IF ( ABS(  JVS( 152 )) < TINY(a) ) THEN
         IER = 30                                      
         RETURN
  END IF
   W( 18 ) = JVS( 150 )
   W( 26 ) = JVS( 151 )
   W( 30 ) = JVS( 152 )
   W( 50 ) = JVS( 153 )
   W( 51 ) = JVS( 154 )
   W( 54 ) = JVS( 155 )
   W( 56 ) = JVS( 156 )
  a = -W( 18 ) / JVS(           59  )
  W( 18 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 60 )
  W( 56 ) = W( 56 ) + a*JVS( 61 )
  a = -W( 26 ) / JVS(          104  )
  W( 26 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 105 )
  W( 54 ) = W( 54 ) + a*JVS( 106 )
  JVS( 150) = W( 18 )
  JVS( 151) = W( 26 )
  JVS( 152) = W( 30 )
  JVS( 153) = W( 50 )
  JVS( 154) = W( 51 )
  JVS( 155) = W( 54 )
  JVS( 156) = W( 56 )
  IF ( ABS(  JVS( 157 )) < TINY(a) ) THEN
         IER = 31                                      
         RETURN
  END IF
   W( 31 ) = JVS( 157 )
   W( 35 ) = JVS( 158 )
   W( 54 ) = JVS( 159 )
   W( 55 ) = JVS( 160 )
   W( 56 ) = JVS( 161 )
  JVS( 157) = W( 31 )
  JVS( 158) = W( 35 )
  JVS( 159) = W( 54 )
  JVS( 160) = W( 55 )
  JVS( 161) = W( 56 )
  IF ( ABS(  JVS( 162 )) < TINY(a) ) THEN
         IER = 32                                      
         RETURN
  END IF
   W( 32 ) = JVS( 162 )
   W( 35 ) = JVS( 163 )
   W( 50 ) = JVS( 164 )
   W( 55 ) = JVS( 165 )
   W( 56 ) = JVS( 166 )
  JVS( 162) = W( 32 )
  JVS( 163) = W( 35 )
  JVS( 164) = W( 50 )
  JVS( 165) = W( 55 )
  JVS( 166) = W( 56 )
  IF ( ABS(  JVS( 167 )) < TINY(a) ) THEN
         IER = 33                                      
         RETURN
  END IF
   W( 33 ) = JVS( 167 )
   W( 45 ) = JVS( 168 )
   W( 50 ) = JVS( 169 )
   W( 55 ) = JVS( 170 )
   W( 56 ) = JVS( 171 )
  JVS( 167) = W( 33 )
  JVS( 168) = W( 45 )
  JVS( 169) = W( 50 )
  JVS( 170) = W( 55 )
  JVS( 171) = W( 56 )
  IF ( ABS(  JVS( 172 )) < TINY(a) ) THEN
         IER = 34                                      
         RETURN
  END IF
   W( 34 ) = JVS( 172 )
   W( 50 ) = JVS( 173 )
   W( 51 ) = JVS( 174 )
   W( 54 ) = JVS( 175 )
  JVS( 172) = W( 34 )
  JVS( 173) = W( 50 )
  JVS( 174) = W( 51 )
  JVS( 175) = W( 54 )
  IF ( ABS(  JVS( 176 )) < TINY(a) ) THEN
         IER = 35                                      
         RETURN
  END IF
   W( 35 ) = JVS( 176 )
   W( 50 ) = JVS( 177 )
   W( 51 ) = JVS( 178 )
   W( 54 ) = JVS( 179 )
  JVS( 176) = W( 35 )
  JVS( 177) = W( 50 )
  JVS( 178) = W( 51 )
  JVS( 179) = W( 54 )
  IF ( ABS(  JVS( 188 )) < TINY(a) ) THEN
         IER = 36                                      
         RETURN
  END IF
   W( 20 ) = JVS( 180 )
   W( 23 ) = JVS( 181 )
   W( 24 ) = JVS( 182 )
   W( 30 ) = JVS( 183 )
   W( 32 ) = JVS( 184 )
   W( 33 ) = JVS( 185 )
   W( 34 ) = JVS( 186 )
   W( 35 ) = JVS( 187 )
   W( 36 ) = JVS( 188 )
   W( 39 ) = JVS( 189 )
   W( 41 ) = JVS( 190 )
   W( 44 ) = JVS( 191 )
   W( 45 ) = JVS( 192 )
   W( 46 ) = JVS( 193 )
   W( 50 ) = JVS( 194 )
   W( 51 ) = JVS( 195 )
   W( 54 ) = JVS( 196 )
   W( 55 ) = JVS( 197 )
   W( 56 ) = JVS( 198 )
  a = -W( 20 ) / JVS(           74  )
  W( 20 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 75 )
  W( 50 ) = W( 50 ) + a*JVS( 76 )
  W( 55 ) = W( 55 ) + a*JVS( 77 )
  a = -W( 23 ) / JVS(           87  )
  W( 23 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 88 )
  W( 51 ) = W( 51 ) + a*JVS( 89 )
  a = -W( 24 ) / JVS(           90  )
  W( 24 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 91 )
  W( 39 ) = W( 39 ) + a*JVS( 92 )
  W( 44 ) = W( 44 ) + a*JVS( 93 )
  W( 50 ) = W( 50 ) + a*JVS( 94 )
  W( 51 ) = W( 51 ) + a*JVS( 95 )
  a = -W( 30 ) / JVS(          152  )
  W( 30 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 153 )
  W( 51 ) = W( 51 ) + a*JVS( 154 )
  W( 54 ) = W( 54 ) + a*JVS( 155 )
  W( 56 ) = W( 56 ) + a*JVS( 156 )
  a = -W( 32 ) / JVS(          162  )
  W( 32 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 163 )
  W( 50 ) = W( 50 ) + a*JVS( 164 )
  W( 55 ) = W( 55 ) + a*JVS( 165 )
  W( 56 ) = W( 56 ) + a*JVS( 166 )
  a = -W( 33 ) / JVS(          167  )
  W( 33 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 168 )
  W( 50 ) = W( 50 ) + a*JVS( 169 )
  W( 55 ) = W( 55 ) + a*JVS( 170 )
  W( 56 ) = W( 56 ) + a*JVS( 171 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 173 )
  W( 51 ) = W( 51 ) + a*JVS( 174 )
  W( 54 ) = W( 54 ) + a*JVS( 175 )
  a = -W( 35 ) / JVS(          176  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 177 )
  W( 51 ) = W( 51 ) + a*JVS( 178 )
  W( 54 ) = W( 54 ) + a*JVS( 179 )
  JVS( 180) = W( 20 )
  JVS( 181) = W( 23 )
  JVS( 182) = W( 24 )
  JVS( 183) = W( 30 )
  JVS( 184) = W( 32 )
  JVS( 185) = W( 33 )
  JVS( 186) = W( 34 )
  JVS( 187) = W( 35 )
  JVS( 188) = W( 36 )
  JVS( 189) = W( 39 )
  JVS( 190) = W( 41 )
  JVS( 191) = W( 44 )
  JVS( 192) = W( 45 )
  JVS( 193) = W( 46 )
  JVS( 194) = W( 50 )
  JVS( 195) = W( 51 )
  JVS( 196) = W( 54 )
  JVS( 197) = W( 55 )
  JVS( 198) = W( 56 )
  IF ( ABS(  JVS( 206 )) < TINY(a) ) THEN
         IER = 37                                      
         RETURN
  END IF
   W( 11 ) = JVS( 199 )
   W( 13 ) = JVS( 200 )
   W( 23 ) = JVS( 201 )
   W( 26 ) = JVS( 202 )
   W( 30 ) = JVS( 203 )
   W( 34 ) = JVS( 204 )
   W( 35 ) = JVS( 205 )
   W( 37 ) = JVS( 206 )
   W( 39 ) = JVS( 207 )
   W( 42 ) = JVS( 208 )
   W( 45 ) = JVS( 209 )
   W( 47 ) = JVS( 210 )
   W( 48 ) = JVS( 211 )
   W( 49 ) = JVS( 212 )
   W( 50 ) = JVS( 213 )
   W( 51 ) = JVS( 214 )
   W( 54 ) = JVS( 215 )
   W( 55 ) = JVS( 216 )
   W( 56 ) = JVS( 217 )
  a = -W( 11 ) / JVS(           30  )
  W( 11 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 31 )
  a = -W( 13 ) / JVS(           35  )
  W( 13 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 36 )
  a = -W( 23 ) / JVS(           87  )
  W( 23 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 88 )
  W( 51 ) = W( 51 ) + a*JVS( 89 )
  a = -W( 26 ) / JVS(          104  )
  W( 26 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 105 )
  W( 54 ) = W( 54 ) + a*JVS( 106 )
  a = -W( 30 ) / JVS(          152  )
  W( 30 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 153 )
  W( 51 ) = W( 51 ) + a*JVS( 154 )
  W( 54 ) = W( 54 ) + a*JVS( 155 )
  W( 56 ) = W( 56 ) + a*JVS( 156 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 173 )
  W( 51 ) = W( 51 ) + a*JVS( 174 )
  W( 54 ) = W( 54 ) + a*JVS( 175 )
  a = -W( 35 ) / JVS(          176  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 177 )
  W( 51 ) = W( 51 ) + a*JVS( 178 )
  W( 54 ) = W( 54 ) + a*JVS( 179 )
  JVS( 199) = W( 11 )
  JVS( 200) = W( 13 )
  JVS( 201) = W( 23 )
  JVS( 202) = W( 26 )
  JVS( 203) = W( 30 )
  JVS( 204) = W( 34 )
  JVS( 205) = W( 35 )
  JVS( 206) = W( 37 )
  JVS( 207) = W( 39 )
  JVS( 208) = W( 42 )
  JVS( 209) = W( 45 )
  JVS( 210) = W( 47 )
  JVS( 211) = W( 48 )
  JVS( 212) = W( 49 )
  JVS( 213) = W( 50 )
  JVS( 214) = W( 51 )
  JVS( 215) = W( 54 )
  JVS( 216) = W( 55 )
  JVS( 217) = W( 56 )
  IF ( ABS(  JVS( 219 )) < TINY(a) ) THEN
         IER = 38                                      
         RETURN
  END IF
   W( 33 ) = JVS( 218 )
   W( 38 ) = JVS( 219 )
   W( 39 ) = JVS( 220 )
   W( 45 ) = JVS( 221 )
   W( 46 ) = JVS( 222 )
   W( 47 ) = JVS( 223 )
   W( 48 ) = JVS( 224 )
   W( 49 ) = JVS( 225 )
   W( 50 ) = JVS( 226 )
   W( 51 ) = JVS( 227 )
   W( 54 ) = JVS( 228 )
   W( 55 ) = JVS( 229 )
   W( 56 ) = JVS( 230 )
  a = -W( 33 ) / JVS(          167  )
  W( 33 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 168 )
  W( 50 ) = W( 50 ) + a*JVS( 169 )
  W( 55 ) = W( 55 ) + a*JVS( 170 )
  W( 56 ) = W( 56 ) + a*JVS( 171 )
  JVS( 218) = W( 33 )
  JVS( 219) = W( 38 )
  JVS( 220) = W( 39 )
  JVS( 221) = W( 45 )
  JVS( 222) = W( 46 )
  JVS( 223) = W( 47 )
  JVS( 224) = W( 48 )
  JVS( 225) = W( 49 )
  JVS( 226) = W( 50 )
  JVS( 227) = W( 51 )
  JVS( 228) = W( 54 )
  JVS( 229) = W( 55 )
  JVS( 230) = W( 56 )
  IF ( ABS(  JVS( 231 )) < TINY(a) ) THEN
         IER = 39                                      
         RETURN
  END IF
   W( 39 ) = JVS( 231 )
   W( 50 ) = JVS( 232 )
   W( 51 ) = JVS( 233 )
   W( 54 ) = JVS( 234 )
  JVS( 231) = W( 39 )
  JVS( 232) = W( 50 )
  JVS( 233) = W( 51 )
  JVS( 234) = W( 54 )
  IF ( ABS(  JVS( 239 )) < TINY(a) ) THEN
         IER = 40                                      
         RETURN
  END IF
   W( 17 ) = JVS( 235 )
   W( 21 ) = JVS( 236 )
   W( 34 ) = JVS( 237 )
   W( 39 ) = JVS( 238 )
   W( 40 ) = JVS( 239 )
   W( 47 ) = JVS( 240 )
   W( 48 ) = JVS( 241 )
   W( 49 ) = JVS( 242 )
   W( 50 ) = JVS( 243 )
   W( 51 ) = JVS( 244 )
   W( 54 ) = JVS( 245 )
   W( 55 ) = JVS( 246 )
   W( 56 ) = JVS( 247 )
  a = -W( 17 ) / JVS(           51  )
  W( 17 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 52 )
  W( 39 ) = W( 39 ) + a*JVS( 53 )
  W( 40 ) = W( 40 ) + a*JVS( 54 )
  W( 50 ) = W( 50 ) + a*JVS( 55 )
  W( 51 ) = W( 51 ) + a*JVS( 56 )
  a = -W( 21 ) / JVS(           78  )
  W( 21 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 79 )
  W( 50 ) = W( 50 ) + a*JVS( 80 )
  W( 55 ) = W( 55 ) + a*JVS( 81 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 173 )
  W( 51 ) = W( 51 ) + a*JVS( 174 )
  W( 54 ) = W( 54 ) + a*JVS( 175 )
  a = -W( 39 ) / JVS(          231  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 54 ) = W( 54 ) + a*JVS( 234 )
  JVS( 235) = W( 17 )
  JVS( 236) = W( 21 )
  JVS( 237) = W( 34 )
  JVS( 238) = W( 39 )
  JVS( 239) = W( 40 )
  JVS( 240) = W( 47 )
  JVS( 241) = W( 48 )
  JVS( 242) = W( 49 )
  JVS( 243) = W( 50 )
  JVS( 244) = W( 51 )
  JVS( 245) = W( 54 )
  JVS( 246) = W( 55 )
  JVS( 247) = W( 56 )
  IF ( ABS(  JVS( 250 )) < TINY(a) ) THEN
         IER = 41                                      
         RETURN
  END IF
   W( 34 ) = JVS( 248 )
   W( 39 ) = JVS( 249 )
   W( 41 ) = JVS( 250 )
   W( 49 ) = JVS( 251 )
   W( 50 ) = JVS( 252 )
   W( 51 ) = JVS( 253 )
   W( 54 ) = JVS( 254 )
   W( 55 ) = JVS( 255 )
   W( 56 ) = JVS( 256 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 173 )
  W( 51 ) = W( 51 ) + a*JVS( 174 )
  W( 54 ) = W( 54 ) + a*JVS( 175 )
  a = -W( 39 ) / JVS(          231  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 54 ) = W( 54 ) + a*JVS( 234 )
  JVS( 248) = W( 34 )
  JVS( 249) = W( 39 )
  JVS( 250) = W( 41 )
  JVS( 251) = W( 49 )
  JVS( 252) = W( 50 )
  JVS( 253) = W( 51 )
  JVS( 254) = W( 54 )
  JVS( 255) = W( 55 )
  JVS( 256) = W( 56 )
  IF ( ABS(  JVS( 262 )) < TINY(a) ) THEN
         IER = 42                                      
         RETURN
  END IF
   W( 13 ) = JVS( 257 )
   W( 30 ) = JVS( 258 )
   W( 33 ) = JVS( 259 )
   W( 34 ) = JVS( 260 )
   W( 39 ) = JVS( 261 )
   W( 42 ) = JVS( 262 )
   W( 45 ) = JVS( 263 )
   W( 46 ) = JVS( 264 )
   W( 47 ) = JVS( 265 )
   W( 50 ) = JVS( 266 )
   W( 51 ) = JVS( 267 )
   W( 54 ) = JVS( 268 )
   W( 55 ) = JVS( 269 )
   W( 56 ) = JVS( 270 )
  a = -W( 13 ) / JVS(           35  )
  W( 13 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 36 )
  a = -W( 30 ) / JVS(          152  )
  W( 30 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 153 )
  W( 51 ) = W( 51 ) + a*JVS( 154 )
  W( 54 ) = W( 54 ) + a*JVS( 155 )
  W( 56 ) = W( 56 ) + a*JVS( 156 )
  a = -W( 33 ) / JVS(          167  )
  W( 33 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 168 )
  W( 50 ) = W( 50 ) + a*JVS( 169 )
  W( 55 ) = W( 55 ) + a*JVS( 170 )
  W( 56 ) = W( 56 ) + a*JVS( 171 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 173 )
  W( 51 ) = W( 51 ) + a*JVS( 174 )
  W( 54 ) = W( 54 ) + a*JVS( 175 )
  a = -W( 39 ) / JVS(          231  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 54 ) = W( 54 ) + a*JVS( 234 )
  JVS( 257) = W( 13 )
  JVS( 258) = W( 30 )
  JVS( 259) = W( 33 )
  JVS( 260) = W( 34 )
  JVS( 261) = W( 39 )
  JVS( 262) = W( 42 )
  JVS( 263) = W( 45 )
  JVS( 264) = W( 46 )
  JVS( 265) = W( 47 )
  JVS( 266) = W( 50 )
  JVS( 267) = W( 51 )
  JVS( 268) = W( 54 )
  JVS( 269) = W( 55 )
  JVS( 270) = W( 56 )
  IF ( ABS(  JVS( 282 )) < TINY(a) ) THEN
         IER = 43                                      
         RETURN
  END IF
   W( 8 ) = JVS( 271 )
   W( 21 ) = JVS( 272 )
   W( 23 ) = JVS( 273 )
   W( 30 ) = JVS( 274 )
   W( 31 ) = JVS( 275 )
   W( 33 ) = JVS( 276 )
   W( 34 ) = JVS( 277 )
   W( 35 ) = JVS( 278 )
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
   W( 54 ) = JVS( 289 )
   W( 55 ) = JVS( 290 )
   W( 56 ) = JVS( 291 )
  a = -W( 8 ) / JVS(           22  )
  W( 8 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 23 )
  a = -W( 21 ) / JVS(           78  )
  W( 21 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 79 )
  W( 50 ) = W( 50 ) + a*JVS( 80 )
  W( 55 ) = W( 55 ) + a*JVS( 81 )
  a = -W( 23 ) / JVS(           87  )
  W( 23 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 88 )
  W( 51 ) = W( 51 ) + a*JVS( 89 )
  a = -W( 30 ) / JVS(          152  )
  W( 30 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 153 )
  W( 51 ) = W( 51 ) + a*JVS( 154 )
  W( 54 ) = W( 54 ) + a*JVS( 155 )
  W( 56 ) = W( 56 ) + a*JVS( 156 )
  a = -W( 31 ) / JVS(          157  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 158 )
  W( 54 ) = W( 54 ) + a*JVS( 159 )
  W( 55 ) = W( 55 ) + a*JVS( 160 )
  W( 56 ) = W( 56 ) + a*JVS( 161 )
  a = -W( 33 ) / JVS(          167  )
  W( 33 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 168 )
  W( 50 ) = W( 50 ) + a*JVS( 169 )
  W( 55 ) = W( 55 ) + a*JVS( 170 )
  W( 56 ) = W( 56 ) + a*JVS( 171 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 173 )
  W( 51 ) = W( 51 ) + a*JVS( 174 )
  W( 54 ) = W( 54 ) + a*JVS( 175 )
  a = -W( 35 ) / JVS(          176  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 177 )
  W( 51 ) = W( 51 ) + a*JVS( 178 )
  W( 54 ) = W( 54 ) + a*JVS( 179 )
  a = -W( 39 ) / JVS(          231  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 54 ) = W( 54 ) + a*JVS( 234 )
  a = -W( 40 ) / JVS(          239  )
  W( 40 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 240 )
  W( 48 ) = W( 48 ) + a*JVS( 241 )
  W( 49 ) = W( 49 ) + a*JVS( 242 )
  W( 50 ) = W( 50 ) + a*JVS( 243 )
  W( 51 ) = W( 51 ) + a*JVS( 244 )
  W( 54 ) = W( 54 ) + a*JVS( 245 )
  W( 55 ) = W( 55 ) + a*JVS( 246 )
  W( 56 ) = W( 56 ) + a*JVS( 247 )
  a = -W( 41 ) / JVS(          250  )
  W( 41 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 251 )
  W( 50 ) = W( 50 ) + a*JVS( 252 )
  W( 51 ) = W( 51 ) + a*JVS( 253 )
  W( 54 ) = W( 54 ) + a*JVS( 254 )
  W( 55 ) = W( 55 ) + a*JVS( 255 )
  W( 56 ) = W( 56 ) + a*JVS( 256 )
  JVS( 271) = W( 8 )
  JVS( 272) = W( 21 )
  JVS( 273) = W( 23 )
  JVS( 274) = W( 30 )
  JVS( 275) = W( 31 )
  JVS( 276) = W( 33 )
  JVS( 277) = W( 34 )
  JVS( 278) = W( 35 )
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
  JVS( 289) = W( 54 )
  JVS( 290) = W( 55 )
  JVS( 291) = W( 56 )
  IF ( ABS(  JVS( 298 )) < TINY(a) ) THEN
         IER = 44                                      
         RETURN
  END IF
   W( 14 ) = JVS( 292 )
   W( 20 ) = JVS( 293 )
   W( 34 ) = JVS( 294 )
   W( 38 ) = JVS( 295 )
   W( 39 ) = JVS( 296 )
   W( 43 ) = JVS( 297 )
   W( 44 ) = JVS( 298 )
   W( 45 ) = JVS( 299 )
   W( 46 ) = JVS( 300 )
   W( 47 ) = JVS( 301 )
   W( 48 ) = JVS( 302 )
   W( 49 ) = JVS( 303 )
   W( 50 ) = JVS( 304 )
   W( 51 ) = JVS( 305 )
   W( 53 ) = JVS( 306 )
   W( 54 ) = JVS( 307 )
   W( 55 ) = JVS( 308 )
   W( 56 ) = JVS( 309 )
  a = -W( 14 ) / JVS(           37  )
  W( 14 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 38 )
  W( 39 ) = W( 39 ) + a*JVS( 39 )
  W( 50 ) = W( 50 ) + a*JVS( 40 )
  W( 51 ) = W( 51 ) + a*JVS( 41 )
  a = -W( 20 ) / JVS(           74  )
  W( 20 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 75 )
  W( 50 ) = W( 50 ) + a*JVS( 76 )
  W( 55 ) = W( 55 ) + a*JVS( 77 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 173 )
  W( 51 ) = W( 51 ) + a*JVS( 174 )
  W( 54 ) = W( 54 ) + a*JVS( 175 )
  a = -W( 38 ) / JVS(          219  )
  W( 38 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 220 )
  W( 45 ) = W( 45 ) + a*JVS( 221 )
  W( 46 ) = W( 46 ) + a*JVS( 222 )
  W( 47 ) = W( 47 ) + a*JVS( 223 )
  W( 48 ) = W( 48 ) + a*JVS( 224 )
  W( 49 ) = W( 49 ) + a*JVS( 225 )
  W( 50 ) = W( 50 ) + a*JVS( 226 )
  W( 51 ) = W( 51 ) + a*JVS( 227 )
  W( 54 ) = W( 54 ) + a*JVS( 228 )
  W( 55 ) = W( 55 ) + a*JVS( 229 )
  W( 56 ) = W( 56 ) + a*JVS( 230 )
  a = -W( 39 ) / JVS(          231  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 54 ) = W( 54 ) + a*JVS( 234 )
  a = -W( 43 ) / JVS(          282  )
  W( 43 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 283 )
  W( 47 ) = W( 47 ) + a*JVS( 284 )
  W( 48 ) = W( 48 ) + a*JVS( 285 )
  W( 49 ) = W( 49 ) + a*JVS( 286 )
  W( 50 ) = W( 50 ) + a*JVS( 287 )
  W( 51 ) = W( 51 ) + a*JVS( 288 )
  W( 54 ) = W( 54 ) + a*JVS( 289 )
  W( 55 ) = W( 55 ) + a*JVS( 290 )
  W( 56 ) = W( 56 ) + a*JVS( 291 )
  JVS( 292) = W( 14 )
  JVS( 293) = W( 20 )
  JVS( 294) = W( 34 )
  JVS( 295) = W( 38 )
  JVS( 296) = W( 39 )
  JVS( 297) = W( 43 )
  JVS( 298) = W( 44 )
  JVS( 299) = W( 45 )
  JVS( 300) = W( 46 )
  JVS( 301) = W( 47 )
  JVS( 302) = W( 48 )
  JVS( 303) = W( 49 )
  JVS( 304) = W( 50 )
  JVS( 305) = W( 51 )
  JVS( 306) = W( 53 )
  JVS( 307) = W( 54 )
  JVS( 308) = W( 55 )
  JVS( 309) = W( 56 )
  IF ( ABS(  JVS( 313 )) < TINY(a) ) THEN
         IER = 45                                      
         RETURN
  END IF
   W( 31 ) = JVS( 310 )
   W( 32 ) = JVS( 311 )
   W( 35 ) = JVS( 312 )
   W( 45 ) = JVS( 313 )
   W( 50 ) = JVS( 314 )
   W( 51 ) = JVS( 315 )
   W( 54 ) = JVS( 316 )
   W( 55 ) = JVS( 317 )
   W( 56 ) = JVS( 318 )
  a = -W( 31 ) / JVS(          157  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 158 )
  W( 54 ) = W( 54 ) + a*JVS( 159 )
  W( 55 ) = W( 55 ) + a*JVS( 160 )
  W( 56 ) = W( 56 ) + a*JVS( 161 )
  a = -W( 32 ) / JVS(          162  )
  W( 32 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 163 )
  W( 50 ) = W( 50 ) + a*JVS( 164 )
  W( 55 ) = W( 55 ) + a*JVS( 165 )
  W( 56 ) = W( 56 ) + a*JVS( 166 )
  a = -W( 35 ) / JVS(          176  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 177 )
  W( 51 ) = W( 51 ) + a*JVS( 178 )
  W( 54 ) = W( 54 ) + a*JVS( 179 )
  JVS( 310) = W( 31 )
  JVS( 311) = W( 32 )
  JVS( 312) = W( 35 )
  JVS( 313) = W( 45 )
  JVS( 314) = W( 50 )
  JVS( 315) = W( 51 )
  JVS( 316) = W( 54 )
  JVS( 317) = W( 55 )
  JVS( 318) = W( 56 )
  IF ( ABS(  JVS( 322 )) < TINY(a) ) THEN
         IER = 46                                      
         RETURN
  END IF
   W( 38 ) = JVS( 319 )
   W( 39 ) = JVS( 320 )
   W( 45 ) = JVS( 321 )
   W( 46 ) = JVS( 322 )
   W( 47 ) = JVS( 323 )
   W( 48 ) = JVS( 324 )
   W( 49 ) = JVS( 325 )
   W( 50 ) = JVS( 326 )
   W( 51 ) = JVS( 327 )
   W( 54 ) = JVS( 328 )
   W( 55 ) = JVS( 329 )
   W( 56 ) = JVS( 330 )
  a = -W( 38 ) / JVS(          219  )
  W( 38 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 220 )
  W( 45 ) = W( 45 ) + a*JVS( 221 )
  W( 46 ) = W( 46 ) + a*JVS( 222 )
  W( 47 ) = W( 47 ) + a*JVS( 223 )
  W( 48 ) = W( 48 ) + a*JVS( 224 )
  W( 49 ) = W( 49 ) + a*JVS( 225 )
  W( 50 ) = W( 50 ) + a*JVS( 226 )
  W( 51 ) = W( 51 ) + a*JVS( 227 )
  W( 54 ) = W( 54 ) + a*JVS( 228 )
  W( 55 ) = W( 55 ) + a*JVS( 229 )
  W( 56 ) = W( 56 ) + a*JVS( 230 )
  a = -W( 39 ) / JVS(          231  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 54 ) = W( 54 ) + a*JVS( 234 )
  a = -W( 45 ) / JVS(          313  )
  W( 45 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 314 )
  W( 51 ) = W( 51 ) + a*JVS( 315 )
  W( 54 ) = W( 54 ) + a*JVS( 316 )
  W( 55 ) = W( 55 ) + a*JVS( 317 )
  W( 56 ) = W( 56 ) + a*JVS( 318 )
  JVS( 319) = W( 38 )
  JVS( 320) = W( 39 )
  JVS( 321) = W( 45 )
  JVS( 322) = W( 46 )
  JVS( 323) = W( 47 )
  JVS( 324) = W( 48 )
  JVS( 325) = W( 49 )
  JVS( 326) = W( 50 )
  JVS( 327) = W( 51 )
  JVS( 328) = W( 54 )
  JVS( 329) = W( 55 )
  JVS( 330) = W( 56 )
  IF ( ABS(  JVS( 336 )) < TINY(a) ) THEN
         IER = 47                                      
         RETURN
  END IF
   W( 32 ) = JVS( 331 )
   W( 33 ) = JVS( 332 )
   W( 35 ) = JVS( 333 )
   W( 45 ) = JVS( 334 )
   W( 46 ) = JVS( 335 )
   W( 47 ) = JVS( 336 )
   W( 48 ) = JVS( 337 )
   W( 49 ) = JVS( 338 )
   W( 50 ) = JVS( 339 )
   W( 51 ) = JVS( 340 )
   W( 54 ) = JVS( 341 )
   W( 55 ) = JVS( 342 )
   W( 56 ) = JVS( 343 )
  a = -W( 32 ) / JVS(          162  )
  W( 32 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 163 )
  W( 50 ) = W( 50 ) + a*JVS( 164 )
  W( 55 ) = W( 55 ) + a*JVS( 165 )
  W( 56 ) = W( 56 ) + a*JVS( 166 )
  a = -W( 33 ) / JVS(          167  )
  W( 33 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 168 )
  W( 50 ) = W( 50 ) + a*JVS( 169 )
  W( 55 ) = W( 55 ) + a*JVS( 170 )
  W( 56 ) = W( 56 ) + a*JVS( 171 )
  a = -W( 35 ) / JVS(          176  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 177 )
  W( 51 ) = W( 51 ) + a*JVS( 178 )
  W( 54 ) = W( 54 ) + a*JVS( 179 )
  a = -W( 45 ) / JVS(          313  )
  W( 45 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 314 )
  W( 51 ) = W( 51 ) + a*JVS( 315 )
  W( 54 ) = W( 54 ) + a*JVS( 316 )
  W( 55 ) = W( 55 ) + a*JVS( 317 )
  W( 56 ) = W( 56 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          322  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 323 )
  W( 48 ) = W( 48 ) + a*JVS( 324 )
  W( 49 ) = W( 49 ) + a*JVS( 325 )
  W( 50 ) = W( 50 ) + a*JVS( 326 )
  W( 51 ) = W( 51 ) + a*JVS( 327 )
  W( 54 ) = W( 54 ) + a*JVS( 328 )
  W( 55 ) = W( 55 ) + a*JVS( 329 )
  W( 56 ) = W( 56 ) + a*JVS( 330 )
  JVS( 331) = W( 32 )
  JVS( 332) = W( 33 )
  JVS( 333) = W( 35 )
  JVS( 334) = W( 45 )
  JVS( 335) = W( 46 )
  JVS( 336) = W( 47 )
  JVS( 337) = W( 48 )
  JVS( 338) = W( 49 )
  JVS( 339) = W( 50 )
  JVS( 340) = W( 51 )
  JVS( 341) = W( 54 )
  JVS( 342) = W( 55 )
  JVS( 343) = W( 56 )
  IF ( ABS(  JVS( 353 )) < TINY(a) ) THEN
         IER = 48                                      
         RETURN
  END IF
   W( 29 ) = JVS( 344 )
   W( 31 ) = JVS( 345 )
   W( 32 ) = JVS( 346 )
   W( 34 ) = JVS( 347 )
   W( 35 ) = JVS( 348 )
   W( 39 ) = JVS( 349 )
   W( 41 ) = JVS( 350 )
   W( 45 ) = JVS( 351 )
   W( 47 ) = JVS( 352 )
   W( 48 ) = JVS( 353 )
   W( 49 ) = JVS( 354 )
   W( 50 ) = JVS( 355 )
   W( 51 ) = JVS( 356 )
   W( 54 ) = JVS( 357 )
   W( 55 ) = JVS( 358 )
   W( 56 ) = JVS( 359 )
  a = -W( 29 ) / JVS(          135  )
  W( 29 ) = -a
  W( 31 ) = W( 31 ) + a*JVS( 136 )
  W( 32 ) = W( 32 ) + a*JVS( 137 )
  W( 34 ) = W( 34 ) + a*JVS( 138 )
  W( 39 ) = W( 39 ) + a*JVS( 139 )
  W( 41 ) = W( 41 ) + a*JVS( 140 )
  W( 45 ) = W( 45 ) + a*JVS( 141 )
  W( 47 ) = W( 47 ) + a*JVS( 142 )
  W( 48 ) = W( 48 ) + a*JVS( 143 )
  W( 49 ) = W( 49 ) + a*JVS( 144 )
  W( 50 ) = W( 50 ) + a*JVS( 145 )
  W( 51 ) = W( 51 ) + a*JVS( 146 )
  W( 54 ) = W( 54 ) + a*JVS( 147 )
  W( 55 ) = W( 55 ) + a*JVS( 148 )
  W( 56 ) = W( 56 ) + a*JVS( 149 )
  a = -W( 31 ) / JVS(          157  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 158 )
  W( 54 ) = W( 54 ) + a*JVS( 159 )
  W( 55 ) = W( 55 ) + a*JVS( 160 )
  W( 56 ) = W( 56 ) + a*JVS( 161 )
  a = -W( 32 ) / JVS(          162  )
  W( 32 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 163 )
  W( 50 ) = W( 50 ) + a*JVS( 164 )
  W( 55 ) = W( 55 ) + a*JVS( 165 )
  W( 56 ) = W( 56 ) + a*JVS( 166 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 173 )
  W( 51 ) = W( 51 ) + a*JVS( 174 )
  W( 54 ) = W( 54 ) + a*JVS( 175 )
  a = -W( 35 ) / JVS(          176  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 177 )
  W( 51 ) = W( 51 ) + a*JVS( 178 )
  W( 54 ) = W( 54 ) + a*JVS( 179 )
  a = -W( 39 ) / JVS(          231  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 54 ) = W( 54 ) + a*JVS( 234 )
  a = -W( 41 ) / JVS(          250  )
  W( 41 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 251 )
  W( 50 ) = W( 50 ) + a*JVS( 252 )
  W( 51 ) = W( 51 ) + a*JVS( 253 )
  W( 54 ) = W( 54 ) + a*JVS( 254 )
  W( 55 ) = W( 55 ) + a*JVS( 255 )
  W( 56 ) = W( 56 ) + a*JVS( 256 )
  a = -W( 45 ) / JVS(          313  )
  W( 45 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 314 )
  W( 51 ) = W( 51 ) + a*JVS( 315 )
  W( 54 ) = W( 54 ) + a*JVS( 316 )
  W( 55 ) = W( 55 ) + a*JVS( 317 )
  W( 56 ) = W( 56 ) + a*JVS( 318 )
  a = -W( 47 ) / JVS(          336  )
  W( 47 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 337 )
  W( 49 ) = W( 49 ) + a*JVS( 338 )
  W( 50 ) = W( 50 ) + a*JVS( 339 )
  W( 51 ) = W( 51 ) + a*JVS( 340 )
  W( 54 ) = W( 54 ) + a*JVS( 341 )
  W( 55 ) = W( 55 ) + a*JVS( 342 )
  W( 56 ) = W( 56 ) + a*JVS( 343 )
  JVS( 344) = W( 29 )
  JVS( 345) = W( 31 )
  JVS( 346) = W( 32 )
  JVS( 347) = W( 34 )
  JVS( 348) = W( 35 )
  JVS( 349) = W( 39 )
  JVS( 350) = W( 41 )
  JVS( 351) = W( 45 )
  JVS( 352) = W( 47 )
  JVS( 353) = W( 48 )
  JVS( 354) = W( 49 )
  JVS( 355) = W( 50 )
  JVS( 356) = W( 51 )
  JVS( 357) = W( 54 )
  JVS( 358) = W( 55 )
  JVS( 359) = W( 56 )
  IF ( ABS(  JVS( 369 )) < TINY(a) ) THEN
         IER = 49                                      
         RETURN
  END IF
   W( 15 ) = JVS( 360 )
   W( 18 ) = JVS( 361 )
   W( 26 ) = JVS( 362 )
   W( 31 ) = JVS( 363 )
   W( 32 ) = JVS( 364 )
   W( 35 ) = JVS( 365 )
   W( 41 ) = JVS( 366 )
   W( 45 ) = JVS( 367 )
   W( 48 ) = JVS( 368 )
   W( 49 ) = JVS( 369 )
   W( 50 ) = JVS( 370 )
   W( 51 ) = JVS( 371 )
   W( 52 ) = JVS( 372 )
   W( 54 ) = JVS( 373 )
   W( 55 ) = JVS( 374 )
   W( 56 ) = JVS( 375 )
  a = -W( 15 ) / JVS(           42  )
  W( 15 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 43 )
  W( 50 ) = W( 50 ) + a*JVS( 44 )
  W( 52 ) = W( 52 ) + a*JVS( 45 )
  W( 54 ) = W( 54 ) + a*JVS( 46 )
  a = -W( 18 ) / JVS(           59  )
  W( 18 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 60 )
  W( 56 ) = W( 56 ) + a*JVS( 61 )
  a = -W( 26 ) / JVS(          104  )
  W( 26 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 105 )
  W( 54 ) = W( 54 ) + a*JVS( 106 )
  a = -W( 31 ) / JVS(          157  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 158 )
  W( 54 ) = W( 54 ) + a*JVS( 159 )
  W( 55 ) = W( 55 ) + a*JVS( 160 )
  W( 56 ) = W( 56 ) + a*JVS( 161 )
  a = -W( 32 ) / JVS(          162  )
  W( 32 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 163 )
  W( 50 ) = W( 50 ) + a*JVS( 164 )
  W( 55 ) = W( 55 ) + a*JVS( 165 )
  W( 56 ) = W( 56 ) + a*JVS( 166 )
  a = -W( 35 ) / JVS(          176  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 177 )
  W( 51 ) = W( 51 ) + a*JVS( 178 )
  W( 54 ) = W( 54 ) + a*JVS( 179 )
  a = -W( 41 ) / JVS(          250  )
  W( 41 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 251 )
  W( 50 ) = W( 50 ) + a*JVS( 252 )
  W( 51 ) = W( 51 ) + a*JVS( 253 )
  W( 54 ) = W( 54 ) + a*JVS( 254 )
  W( 55 ) = W( 55 ) + a*JVS( 255 )
  W( 56 ) = W( 56 ) + a*JVS( 256 )
  a = -W( 45 ) / JVS(          313  )
  W( 45 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 314 )
  W( 51 ) = W( 51 ) + a*JVS( 315 )
  W( 54 ) = W( 54 ) + a*JVS( 316 )
  W( 55 ) = W( 55 ) + a*JVS( 317 )
  W( 56 ) = W( 56 ) + a*JVS( 318 )
  a = -W( 48 ) / JVS(          353  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 354 )
  W( 50 ) = W( 50 ) + a*JVS( 355 )
  W( 51 ) = W( 51 ) + a*JVS( 356 )
  W( 54 ) = W( 54 ) + a*JVS( 357 )
  W( 55 ) = W( 55 ) + a*JVS( 358 )
  W( 56 ) = W( 56 ) + a*JVS( 359 )
  JVS( 360) = W( 15 )
  JVS( 361) = W( 18 )
  JVS( 362) = W( 26 )
  JVS( 363) = W( 31 )
  JVS( 364) = W( 32 )
  JVS( 365) = W( 35 )
  JVS( 366) = W( 41 )
  JVS( 367) = W( 45 )
  JVS( 368) = W( 48 )
  JVS( 369) = W( 49 )
  JVS( 370) = W( 50 )
  JVS( 371) = W( 51 )
  JVS( 372) = W( 52 )
  JVS( 373) = W( 54 )
  JVS( 374) = W( 55 )
  JVS( 375) = W( 56 )
  IF ( ABS(  JVS( 413 )) < TINY(a) ) THEN
         IER = 50                                      
         RETURN
  END IF
   W( 6 ) = JVS( 376 )
   W( 7 ) = JVS( 377 )
   W( 8 ) = JVS( 378 )
   W( 9 ) = JVS( 379 )
   W( 11 ) = JVS( 380 )
   W( 13 ) = JVS( 381 )
   W( 14 ) = JVS( 382 )
   W( 16 ) = JVS( 383 )
   W( 17 ) = JVS( 384 )
   W( 20 ) = JVS( 385 )
   W( 21 ) = JVS( 386 )
   W( 22 ) = JVS( 387 )
   W( 23 ) = JVS( 388 )
   W( 24 ) = JVS( 389 )
   W( 26 ) = JVS( 390 )
   W( 27 ) = JVS( 391 )
   W( 28 ) = JVS( 392 )
   W( 29 ) = JVS( 393 )
   W( 30 ) = JVS( 394 )
   W( 31 ) = JVS( 395 )
   W( 32 ) = JVS( 396 )
   W( 33 ) = JVS( 397 )
   W( 34 ) = JVS( 398 )
   W( 35 ) = JVS( 399 )
   W( 36 ) = JVS( 400 )
   W( 38 ) = JVS( 401 )
   W( 39 ) = JVS( 402 )
   W( 40 ) = JVS( 403 )
   W( 41 ) = JVS( 404 )
   W( 42 ) = JVS( 405 )
   W( 43 ) = JVS( 406 )
   W( 44 ) = JVS( 407 )
   W( 45 ) = JVS( 408 )
   W( 46 ) = JVS( 409 )
   W( 47 ) = JVS( 410 )
   W( 48 ) = JVS( 411 )
   W( 49 ) = JVS( 412 )
   W( 50 ) = JVS( 413 )
   W( 51 ) = JVS( 414 )
   W( 52 ) = JVS( 415 )
   W( 53 ) = JVS( 416 )
   W( 54 ) = JVS( 417 )
   W( 55 ) = JVS( 418 )
   W( 56 ) = JVS( 419 )
  a = -W( 6 ) / JVS(           18  )
  W( 6 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 19 )
  a = -W( 7 ) / JVS(           20  )
  W( 7 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 21 )
  a = -W( 8 ) / JVS(           22  )
  W( 8 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 23 )
  a = -W( 9 ) / JVS(           24  )
  W( 9 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 25 )
  W( 55 ) = W( 55 ) + a*JVS( 26 )
  a = -W( 11 ) / JVS(           30  )
  W( 11 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 31 )
  a = -W( 13 ) / JVS(           35  )
  W( 13 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 36 )
  a = -W( 14 ) / JVS(           37  )
  W( 14 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 38 )
  W( 39 ) = W( 39 ) + a*JVS( 39 )
  W( 50 ) = W( 50 ) + a*JVS( 40 )
  W( 51 ) = W( 51 ) + a*JVS( 41 )
  a = -W( 16 ) / JVS(           47  )
  W( 16 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 48 )
  W( 52 ) = W( 52 ) + a*JVS( 49 )
  W( 55 ) = W( 55 ) + a*JVS( 50 )
  a = -W( 17 ) / JVS(           51  )
  W( 17 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 52 )
  W( 39 ) = W( 39 ) + a*JVS( 53 )
  W( 40 ) = W( 40 ) + a*JVS( 54 )
  W( 50 ) = W( 50 ) + a*JVS( 55 )
  W( 51 ) = W( 51 ) + a*JVS( 56 )
  a = -W( 20 ) / JVS(           74  )
  W( 20 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 75 )
  W( 50 ) = W( 50 ) + a*JVS( 76 )
  W( 55 ) = W( 55 ) + a*JVS( 77 )
  a = -W( 21 ) / JVS(           78  )
  W( 21 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 79 )
  W( 50 ) = W( 50 ) + a*JVS( 80 )
  W( 55 ) = W( 55 ) + a*JVS( 81 )
  a = -W( 22 ) / JVS(           82  )
  W( 22 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 83 )
  W( 52 ) = W( 52 ) + a*JVS( 84 )
  W( 55 ) = W( 55 ) + a*JVS( 85 )
  W( 56 ) = W( 56 ) + a*JVS( 86 )
  a = -W( 23 ) / JVS(           87  )
  W( 23 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 88 )
  W( 51 ) = W( 51 ) + a*JVS( 89 )
  a = -W( 24 ) / JVS(           90  )
  W( 24 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 91 )
  W( 39 ) = W( 39 ) + a*JVS( 92 )
  W( 44 ) = W( 44 ) + a*JVS( 93 )
  W( 50 ) = W( 50 ) + a*JVS( 94 )
  W( 51 ) = W( 51 ) + a*JVS( 95 )
  a = -W( 26 ) / JVS(          104  )
  W( 26 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 105 )
  W( 54 ) = W( 54 ) + a*JVS( 106 )
  a = -W( 27 ) / JVS(          109  )
  W( 27 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 110 )
  W( 42 ) = W( 42 ) + a*JVS( 111 )
  W( 43 ) = W( 43 ) + a*JVS( 112 )
  W( 45 ) = W( 45 ) + a*JVS( 113 )
  W( 50 ) = W( 50 ) + a*JVS( 114 )
  W( 52 ) = W( 52 ) + a*JVS( 115 )
  W( 54 ) = W( 54 ) + a*JVS( 116 )
  W( 55 ) = W( 55 ) + a*JVS( 117 )
  a = -W( 28 ) / JVS(          119  )
  W( 28 ) = -a
  W( 30 ) = W( 30 ) + a*JVS( 120 )
  W( 33 ) = W( 33 ) + a*JVS( 121 )
  W( 34 ) = W( 34 ) + a*JVS( 122 )
  W( 35 ) = W( 35 ) + a*JVS( 123 )
  W( 36 ) = W( 36 ) + a*JVS( 124 )
  W( 39 ) = W( 39 ) + a*JVS( 125 )
  W( 42 ) = W( 42 ) + a*JVS( 126 )
  W( 43 ) = W( 43 ) + a*JVS( 127 )
  W( 45 ) = W( 45 ) + a*JVS( 128 )
  W( 50 ) = W( 50 ) + a*JVS( 129 )
  W( 51 ) = W( 51 ) + a*JVS( 130 )
  W( 54 ) = W( 54 ) + a*JVS( 131 )
  W( 56 ) = W( 56 ) + a*JVS( 132 )
  a = -W( 29 ) / JVS(          135  )
  W( 29 ) = -a
  W( 31 ) = W( 31 ) + a*JVS( 136 )
  W( 32 ) = W( 32 ) + a*JVS( 137 )
  W( 34 ) = W( 34 ) + a*JVS( 138 )
  W( 39 ) = W( 39 ) + a*JVS( 139 )
  W( 41 ) = W( 41 ) + a*JVS( 140 )
  W( 45 ) = W( 45 ) + a*JVS( 141 )
  W( 47 ) = W( 47 ) + a*JVS( 142 )
  W( 48 ) = W( 48 ) + a*JVS( 143 )
  W( 49 ) = W( 49 ) + a*JVS( 144 )
  W( 50 ) = W( 50 ) + a*JVS( 145 )
  W( 51 ) = W( 51 ) + a*JVS( 146 )
  W( 54 ) = W( 54 ) + a*JVS( 147 )
  W( 55 ) = W( 55 ) + a*JVS( 148 )
  W( 56 ) = W( 56 ) + a*JVS( 149 )
  a = -W( 30 ) / JVS(          152  )
  W( 30 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 153 )
  W( 51 ) = W( 51 ) + a*JVS( 154 )
  W( 54 ) = W( 54 ) + a*JVS( 155 )
  W( 56 ) = W( 56 ) + a*JVS( 156 )
  a = -W( 31 ) / JVS(          157  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 158 )
  W( 54 ) = W( 54 ) + a*JVS( 159 )
  W( 55 ) = W( 55 ) + a*JVS( 160 )
  W( 56 ) = W( 56 ) + a*JVS( 161 )
  a = -W( 32 ) / JVS(          162  )
  W( 32 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 163 )
  W( 50 ) = W( 50 ) + a*JVS( 164 )
  W( 55 ) = W( 55 ) + a*JVS( 165 )
  W( 56 ) = W( 56 ) + a*JVS( 166 )
  a = -W( 33 ) / JVS(          167  )
  W( 33 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 168 )
  W( 50 ) = W( 50 ) + a*JVS( 169 )
  W( 55 ) = W( 55 ) + a*JVS( 170 )
  W( 56 ) = W( 56 ) + a*JVS( 171 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 173 )
  W( 51 ) = W( 51 ) + a*JVS( 174 )
  W( 54 ) = W( 54 ) + a*JVS( 175 )
  a = -W( 35 ) / JVS(          176  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 177 )
  W( 51 ) = W( 51 ) + a*JVS( 178 )
  W( 54 ) = W( 54 ) + a*JVS( 179 )
  a = -W( 36 ) / JVS(          188  )
  W( 36 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 189 )
  W( 41 ) = W( 41 ) + a*JVS( 190 )
  W( 44 ) = W( 44 ) + a*JVS( 191 )
  W( 45 ) = W( 45 ) + a*JVS( 192 )
  W( 46 ) = W( 46 ) + a*JVS( 193 )
  W( 50 ) = W( 50 ) + a*JVS( 194 )
  W( 51 ) = W( 51 ) + a*JVS( 195 )
  W( 54 ) = W( 54 ) + a*JVS( 196 )
  W( 55 ) = W( 55 ) + a*JVS( 197 )
  W( 56 ) = W( 56 ) + a*JVS( 198 )
  a = -W( 38 ) / JVS(          219  )
  W( 38 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 220 )
  W( 45 ) = W( 45 ) + a*JVS( 221 )
  W( 46 ) = W( 46 ) + a*JVS( 222 )
  W( 47 ) = W( 47 ) + a*JVS( 223 )
  W( 48 ) = W( 48 ) + a*JVS( 224 )
  W( 49 ) = W( 49 ) + a*JVS( 225 )
  W( 50 ) = W( 50 ) + a*JVS( 226 )
  W( 51 ) = W( 51 ) + a*JVS( 227 )
  W( 54 ) = W( 54 ) + a*JVS( 228 )
  W( 55 ) = W( 55 ) + a*JVS( 229 )
  W( 56 ) = W( 56 ) + a*JVS( 230 )
  a = -W( 39 ) / JVS(          231  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 54 ) = W( 54 ) + a*JVS( 234 )
  a = -W( 40 ) / JVS(          239  )
  W( 40 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 240 )
  W( 48 ) = W( 48 ) + a*JVS( 241 )
  W( 49 ) = W( 49 ) + a*JVS( 242 )
  W( 50 ) = W( 50 ) + a*JVS( 243 )
  W( 51 ) = W( 51 ) + a*JVS( 244 )
  W( 54 ) = W( 54 ) + a*JVS( 245 )
  W( 55 ) = W( 55 ) + a*JVS( 246 )
  W( 56 ) = W( 56 ) + a*JVS( 247 )
  a = -W( 41 ) / JVS(          250  )
  W( 41 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 251 )
  W( 50 ) = W( 50 ) + a*JVS( 252 )
  W( 51 ) = W( 51 ) + a*JVS( 253 )
  W( 54 ) = W( 54 ) + a*JVS( 254 )
  W( 55 ) = W( 55 ) + a*JVS( 255 )
  W( 56 ) = W( 56 ) + a*JVS( 256 )
  a = -W( 42 ) / JVS(          262  )
  W( 42 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 263 )
  W( 46 ) = W( 46 ) + a*JVS( 264 )
  W( 47 ) = W( 47 ) + a*JVS( 265 )
  W( 50 ) = W( 50 ) + a*JVS( 266 )
  W( 51 ) = W( 51 ) + a*JVS( 267 )
  W( 54 ) = W( 54 ) + a*JVS( 268 )
  W( 55 ) = W( 55 ) + a*JVS( 269 )
  W( 56 ) = W( 56 ) + a*JVS( 270 )
  a = -W( 43 ) / JVS(          282  )
  W( 43 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 283 )
  W( 47 ) = W( 47 ) + a*JVS( 284 )
  W( 48 ) = W( 48 ) + a*JVS( 285 )
  W( 49 ) = W( 49 ) + a*JVS( 286 )
  W( 50 ) = W( 50 ) + a*JVS( 287 )
  W( 51 ) = W( 51 ) + a*JVS( 288 )
  W( 54 ) = W( 54 ) + a*JVS( 289 )
  W( 55 ) = W( 55 ) + a*JVS( 290 )
  W( 56 ) = W( 56 ) + a*JVS( 291 )
  a = -W( 44 ) / JVS(          298  )
  W( 44 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 299 )
  W( 46 ) = W( 46 ) + a*JVS( 300 )
  W( 47 ) = W( 47 ) + a*JVS( 301 )
  W( 48 ) = W( 48 ) + a*JVS( 302 )
  W( 49 ) = W( 49 ) + a*JVS( 303 )
  W( 50 ) = W( 50 ) + a*JVS( 304 )
  W( 51 ) = W( 51 ) + a*JVS( 305 )
  W( 53 ) = W( 53 ) + a*JVS( 306 )
  W( 54 ) = W( 54 ) + a*JVS( 307 )
  W( 55 ) = W( 55 ) + a*JVS( 308 )
  W( 56 ) = W( 56 ) + a*JVS( 309 )
  a = -W( 45 ) / JVS(          313  )
  W( 45 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 314 )
  W( 51 ) = W( 51 ) + a*JVS( 315 )
  W( 54 ) = W( 54 ) + a*JVS( 316 )
  W( 55 ) = W( 55 ) + a*JVS( 317 )
  W( 56 ) = W( 56 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          322  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 323 )
  W( 48 ) = W( 48 ) + a*JVS( 324 )
  W( 49 ) = W( 49 ) + a*JVS( 325 )
  W( 50 ) = W( 50 ) + a*JVS( 326 )
  W( 51 ) = W( 51 ) + a*JVS( 327 )
  W( 54 ) = W( 54 ) + a*JVS( 328 )
  W( 55 ) = W( 55 ) + a*JVS( 329 )
  W( 56 ) = W( 56 ) + a*JVS( 330 )
  a = -W( 47 ) / JVS(          336  )
  W( 47 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 337 )
  W( 49 ) = W( 49 ) + a*JVS( 338 )
  W( 50 ) = W( 50 ) + a*JVS( 339 )
  W( 51 ) = W( 51 ) + a*JVS( 340 )
  W( 54 ) = W( 54 ) + a*JVS( 341 )
  W( 55 ) = W( 55 ) + a*JVS( 342 )
  W( 56 ) = W( 56 ) + a*JVS( 343 )
  a = -W( 48 ) / JVS(          353  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 354 )
  W( 50 ) = W( 50 ) + a*JVS( 355 )
  W( 51 ) = W( 51 ) + a*JVS( 356 )
  W( 54 ) = W( 54 ) + a*JVS( 357 )
  W( 55 ) = W( 55 ) + a*JVS( 358 )
  W( 56 ) = W( 56 ) + a*JVS( 359 )
  a = -W( 49 ) / JVS(          369  )
  W( 49 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 370 )
  W( 51 ) = W( 51 ) + a*JVS( 371 )
  W( 52 ) = W( 52 ) + a*JVS( 372 )
  W( 54 ) = W( 54 ) + a*JVS( 373 )
  W( 55 ) = W( 55 ) + a*JVS( 374 )
  W( 56 ) = W( 56 ) + a*JVS( 375 )
  JVS( 376) = W( 6 )
  JVS( 377) = W( 7 )
  JVS( 378) = W( 8 )
  JVS( 379) = W( 9 )
  JVS( 380) = W( 11 )
  JVS( 381) = W( 13 )
  JVS( 382) = W( 14 )
  JVS( 383) = W( 16 )
  JVS( 384) = W( 17 )
  JVS( 385) = W( 20 )
  JVS( 386) = W( 21 )
  JVS( 387) = W( 22 )
  JVS( 388) = W( 23 )
  JVS( 389) = W( 24 )
  JVS( 390) = W( 26 )
  JVS( 391) = W( 27 )
  JVS( 392) = W( 28 )
  JVS( 393) = W( 29 )
  JVS( 394) = W( 30 )
  JVS( 395) = W( 31 )
  JVS( 396) = W( 32 )
  JVS( 397) = W( 33 )
  JVS( 398) = W( 34 )
  JVS( 399) = W( 35 )
  JVS( 400) = W( 36 )
  JVS( 401) = W( 38 )
  JVS( 402) = W( 39 )
  JVS( 403) = W( 40 )
  JVS( 404) = W( 41 )
  JVS( 405) = W( 42 )
  JVS( 406) = W( 43 )
  JVS( 407) = W( 44 )
  JVS( 408) = W( 45 )
  JVS( 409) = W( 46 )
  JVS( 410) = W( 47 )
  JVS( 411) = W( 48 )
  JVS( 412) = W( 49 )
  JVS( 413) = W( 50 )
  JVS( 414) = W( 51 )
  JVS( 415) = W( 52 )
  JVS( 416) = W( 53 )
  JVS( 417) = W( 54 )
  JVS( 418) = W( 55 )
  JVS( 419) = W( 56 )
  IF ( ABS(  JVS( 428 )) < TINY(a) ) THEN
         IER = 51                                      
         RETURN
  END IF
   W( 23 ) = JVS( 420 )
   W( 25 ) = JVS( 421 )
   W( 30 ) = JVS( 422 )
   W( 34 ) = JVS( 423 )
   W( 35 ) = JVS( 424 )
   W( 39 ) = JVS( 425 )
   W( 45 ) = JVS( 426 )
   W( 50 ) = JVS( 427 )
   W( 51 ) = JVS( 428 )
   W( 52 ) = JVS( 429 )
   W( 53 ) = JVS( 430 )
   W( 54 ) = JVS( 431 )
   W( 55 ) = JVS( 432 )
   W( 56 ) = JVS( 433 )
  a = -W( 23 ) / JVS(           87  )
  W( 23 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 88 )
  W( 51 ) = W( 51 ) + a*JVS( 89 )
  a = -W( 25 ) / JVS(           97  )
  W( 25 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 98 )
  W( 52 ) = W( 52 ) + a*JVS( 99 )
  W( 54 ) = W( 54 ) + a*JVS( 100 )
  W( 56 ) = W( 56 ) + a*JVS( 101 )
  a = -W( 30 ) / JVS(          152  )
  W( 30 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 153 )
  W( 51 ) = W( 51 ) + a*JVS( 154 )
  W( 54 ) = W( 54 ) + a*JVS( 155 )
  W( 56 ) = W( 56 ) + a*JVS( 156 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 173 )
  W( 51 ) = W( 51 ) + a*JVS( 174 )
  W( 54 ) = W( 54 ) + a*JVS( 175 )
  a = -W( 35 ) / JVS(          176  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 177 )
  W( 51 ) = W( 51 ) + a*JVS( 178 )
  W( 54 ) = W( 54 ) + a*JVS( 179 )
  a = -W( 39 ) / JVS(          231  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 54 ) = W( 54 ) + a*JVS( 234 )
  a = -W( 45 ) / JVS(          313  )
  W( 45 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 314 )
  W( 51 ) = W( 51 ) + a*JVS( 315 )
  W( 54 ) = W( 54 ) + a*JVS( 316 )
  W( 55 ) = W( 55 ) + a*JVS( 317 )
  W( 56 ) = W( 56 ) + a*JVS( 318 )
  a = -W( 50 ) / JVS(          413  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 414 )
  W( 52 ) = W( 52 ) + a*JVS( 415 )
  W( 53 ) = W( 53 ) + a*JVS( 416 )
  W( 54 ) = W( 54 ) + a*JVS( 417 )
  W( 55 ) = W( 55 ) + a*JVS( 418 )
  W( 56 ) = W( 56 ) + a*JVS( 419 )
  JVS( 420) = W( 23 )
  JVS( 421) = W( 25 )
  JVS( 422) = W( 30 )
  JVS( 423) = W( 34 )
  JVS( 424) = W( 35 )
  JVS( 425) = W( 39 )
  JVS( 426) = W( 45 )
  JVS( 427) = W( 50 )
  JVS( 428) = W( 51 )
  JVS( 429) = W( 52 )
  JVS( 430) = W( 53 )
  JVS( 431) = W( 54 )
  JVS( 432) = W( 55 )
  JVS( 433) = W( 56 )
  IF ( ABS(  JVS( 462 )) < TINY(a) ) THEN
         IER = 52                                      
         RETURN
  END IF
   W( 10 ) = JVS( 434 )
   W( 12 ) = JVS( 435 )
   W( 15 ) = JVS( 436 )
   W( 16 ) = JVS( 437 )
   W( 18 ) = JVS( 438 )
   W( 22 ) = JVS( 439 )
   W( 25 ) = JVS( 440 )
   W( 26 ) = JVS( 441 )
   W( 27 ) = JVS( 442 )
   W( 31 ) = JVS( 443 )
   W( 32 ) = JVS( 444 )
   W( 33 ) = JVS( 445 )
   W( 35 ) = JVS( 446 )
   W( 36 ) = JVS( 447 )
   W( 37 ) = JVS( 448 )
   W( 39 ) = JVS( 449 )
   W( 40 ) = JVS( 450 )
   W( 41 ) = JVS( 451 )
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
  a = -W( 10 ) / JVS(           27  )
  W( 10 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 28 )
  W( 53 ) = W( 53 ) + a*JVS( 29 )
  a = -W( 12 ) / JVS(           32  )
  W( 12 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 33 )
  W( 54 ) = W( 54 ) + a*JVS( 34 )
  a = -W( 15 ) / JVS(           42  )
  W( 15 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 43 )
  W( 50 ) = W( 50 ) + a*JVS( 44 )
  W( 52 ) = W( 52 ) + a*JVS( 45 )
  W( 54 ) = W( 54 ) + a*JVS( 46 )
  a = -W( 16 ) / JVS(           47  )
  W( 16 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 48 )
  W( 52 ) = W( 52 ) + a*JVS( 49 )
  W( 55 ) = W( 55 ) + a*JVS( 50 )
  a = -W( 18 ) / JVS(           59  )
  W( 18 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 60 )
  W( 56 ) = W( 56 ) + a*JVS( 61 )
  a = -W( 22 ) / JVS(           82  )
  W( 22 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 83 )
  W( 52 ) = W( 52 ) + a*JVS( 84 )
  W( 55 ) = W( 55 ) + a*JVS( 85 )
  W( 56 ) = W( 56 ) + a*JVS( 86 )
  a = -W( 25 ) / JVS(           97  )
  W( 25 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 98 )
  W( 52 ) = W( 52 ) + a*JVS( 99 )
  W( 54 ) = W( 54 ) + a*JVS( 100 )
  W( 56 ) = W( 56 ) + a*JVS( 101 )
  a = -W( 26 ) / JVS(          104  )
  W( 26 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 105 )
  W( 54 ) = W( 54 ) + a*JVS( 106 )
  a = -W( 27 ) / JVS(          109  )
  W( 27 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 110 )
  W( 42 ) = W( 42 ) + a*JVS( 111 )
  W( 43 ) = W( 43 ) + a*JVS( 112 )
  W( 45 ) = W( 45 ) + a*JVS( 113 )
  W( 50 ) = W( 50 ) + a*JVS( 114 )
  W( 52 ) = W( 52 ) + a*JVS( 115 )
  W( 54 ) = W( 54 ) + a*JVS( 116 )
  W( 55 ) = W( 55 ) + a*JVS( 117 )
  a = -W( 31 ) / JVS(          157  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 158 )
  W( 54 ) = W( 54 ) + a*JVS( 159 )
  W( 55 ) = W( 55 ) + a*JVS( 160 )
  W( 56 ) = W( 56 ) + a*JVS( 161 )
  a = -W( 32 ) / JVS(          162  )
  W( 32 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 163 )
  W( 50 ) = W( 50 ) + a*JVS( 164 )
  W( 55 ) = W( 55 ) + a*JVS( 165 )
  W( 56 ) = W( 56 ) + a*JVS( 166 )
  a = -W( 33 ) / JVS(          167  )
  W( 33 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 168 )
  W( 50 ) = W( 50 ) + a*JVS( 169 )
  W( 55 ) = W( 55 ) + a*JVS( 170 )
  W( 56 ) = W( 56 ) + a*JVS( 171 )
  a = -W( 35 ) / JVS(          176  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 177 )
  W( 51 ) = W( 51 ) + a*JVS( 178 )
  W( 54 ) = W( 54 ) + a*JVS( 179 )
  a = -W( 36 ) / JVS(          188  )
  W( 36 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 189 )
  W( 41 ) = W( 41 ) + a*JVS( 190 )
  W( 44 ) = W( 44 ) + a*JVS( 191 )
  W( 45 ) = W( 45 ) + a*JVS( 192 )
  W( 46 ) = W( 46 ) + a*JVS( 193 )
  W( 50 ) = W( 50 ) + a*JVS( 194 )
  W( 51 ) = W( 51 ) + a*JVS( 195 )
  W( 54 ) = W( 54 ) + a*JVS( 196 )
  W( 55 ) = W( 55 ) + a*JVS( 197 )
  W( 56 ) = W( 56 ) + a*JVS( 198 )
  a = -W( 37 ) / JVS(          206  )
  W( 37 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 207 )
  W( 42 ) = W( 42 ) + a*JVS( 208 )
  W( 45 ) = W( 45 ) + a*JVS( 209 )
  W( 47 ) = W( 47 ) + a*JVS( 210 )
  W( 48 ) = W( 48 ) + a*JVS( 211 )
  W( 49 ) = W( 49 ) + a*JVS( 212 )
  W( 50 ) = W( 50 ) + a*JVS( 213 )
  W( 51 ) = W( 51 ) + a*JVS( 214 )
  W( 54 ) = W( 54 ) + a*JVS( 215 )
  W( 55 ) = W( 55 ) + a*JVS( 216 )
  W( 56 ) = W( 56 ) + a*JVS( 217 )
  a = -W( 39 ) / JVS(          231  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 54 ) = W( 54 ) + a*JVS( 234 )
  a = -W( 40 ) / JVS(          239  )
  W( 40 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 240 )
  W( 48 ) = W( 48 ) + a*JVS( 241 )
  W( 49 ) = W( 49 ) + a*JVS( 242 )
  W( 50 ) = W( 50 ) + a*JVS( 243 )
  W( 51 ) = W( 51 ) + a*JVS( 244 )
  W( 54 ) = W( 54 ) + a*JVS( 245 )
  W( 55 ) = W( 55 ) + a*JVS( 246 )
  W( 56 ) = W( 56 ) + a*JVS( 247 )
  a = -W( 41 ) / JVS(          250  )
  W( 41 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 251 )
  W( 50 ) = W( 50 ) + a*JVS( 252 )
  W( 51 ) = W( 51 ) + a*JVS( 253 )
  W( 54 ) = W( 54 ) + a*JVS( 254 )
  W( 55 ) = W( 55 ) + a*JVS( 255 )
  W( 56 ) = W( 56 ) + a*JVS( 256 )
  a = -W( 42 ) / JVS(          262  )
  W( 42 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 263 )
  W( 46 ) = W( 46 ) + a*JVS( 264 )
  W( 47 ) = W( 47 ) + a*JVS( 265 )
  W( 50 ) = W( 50 ) + a*JVS( 266 )
  W( 51 ) = W( 51 ) + a*JVS( 267 )
  W( 54 ) = W( 54 ) + a*JVS( 268 )
  W( 55 ) = W( 55 ) + a*JVS( 269 )
  W( 56 ) = W( 56 ) + a*JVS( 270 )
  a = -W( 43 ) / JVS(          282  )
  W( 43 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 283 )
  W( 47 ) = W( 47 ) + a*JVS( 284 )
  W( 48 ) = W( 48 ) + a*JVS( 285 )
  W( 49 ) = W( 49 ) + a*JVS( 286 )
  W( 50 ) = W( 50 ) + a*JVS( 287 )
  W( 51 ) = W( 51 ) + a*JVS( 288 )
  W( 54 ) = W( 54 ) + a*JVS( 289 )
  W( 55 ) = W( 55 ) + a*JVS( 290 )
  W( 56 ) = W( 56 ) + a*JVS( 291 )
  a = -W( 44 ) / JVS(          298  )
  W( 44 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 299 )
  W( 46 ) = W( 46 ) + a*JVS( 300 )
  W( 47 ) = W( 47 ) + a*JVS( 301 )
  W( 48 ) = W( 48 ) + a*JVS( 302 )
  W( 49 ) = W( 49 ) + a*JVS( 303 )
  W( 50 ) = W( 50 ) + a*JVS( 304 )
  W( 51 ) = W( 51 ) + a*JVS( 305 )
  W( 53 ) = W( 53 ) + a*JVS( 306 )
  W( 54 ) = W( 54 ) + a*JVS( 307 )
  W( 55 ) = W( 55 ) + a*JVS( 308 )
  W( 56 ) = W( 56 ) + a*JVS( 309 )
  a = -W( 45 ) / JVS(          313  )
  W( 45 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 314 )
  W( 51 ) = W( 51 ) + a*JVS( 315 )
  W( 54 ) = W( 54 ) + a*JVS( 316 )
  W( 55 ) = W( 55 ) + a*JVS( 317 )
  W( 56 ) = W( 56 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          322  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 323 )
  W( 48 ) = W( 48 ) + a*JVS( 324 )
  W( 49 ) = W( 49 ) + a*JVS( 325 )
  W( 50 ) = W( 50 ) + a*JVS( 326 )
  W( 51 ) = W( 51 ) + a*JVS( 327 )
  W( 54 ) = W( 54 ) + a*JVS( 328 )
  W( 55 ) = W( 55 ) + a*JVS( 329 )
  W( 56 ) = W( 56 ) + a*JVS( 330 )
  a = -W( 47 ) / JVS(          336  )
  W( 47 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 337 )
  W( 49 ) = W( 49 ) + a*JVS( 338 )
  W( 50 ) = W( 50 ) + a*JVS( 339 )
  W( 51 ) = W( 51 ) + a*JVS( 340 )
  W( 54 ) = W( 54 ) + a*JVS( 341 )
  W( 55 ) = W( 55 ) + a*JVS( 342 )
  W( 56 ) = W( 56 ) + a*JVS( 343 )
  a = -W( 48 ) / JVS(          353  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 354 )
  W( 50 ) = W( 50 ) + a*JVS( 355 )
  W( 51 ) = W( 51 ) + a*JVS( 356 )
  W( 54 ) = W( 54 ) + a*JVS( 357 )
  W( 55 ) = W( 55 ) + a*JVS( 358 )
  W( 56 ) = W( 56 ) + a*JVS( 359 )
  a = -W( 49 ) / JVS(          369  )
  W( 49 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 370 )
  W( 51 ) = W( 51 ) + a*JVS( 371 )
  W( 52 ) = W( 52 ) + a*JVS( 372 )
  W( 54 ) = W( 54 ) + a*JVS( 373 )
  W( 55 ) = W( 55 ) + a*JVS( 374 )
  W( 56 ) = W( 56 ) + a*JVS( 375 )
  a = -W( 50 ) / JVS(          413  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 414 )
  W( 52 ) = W( 52 ) + a*JVS( 415 )
  W( 53 ) = W( 53 ) + a*JVS( 416 )
  W( 54 ) = W( 54 ) + a*JVS( 417 )
  W( 55 ) = W( 55 ) + a*JVS( 418 )
  W( 56 ) = W( 56 ) + a*JVS( 419 )
  a = -W( 51 ) / JVS(          428  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 429 )
  W( 53 ) = W( 53 ) + a*JVS( 430 )
  W( 54 ) = W( 54 ) + a*JVS( 431 )
  W( 55 ) = W( 55 ) + a*JVS( 432 )
  W( 56 ) = W( 56 ) + a*JVS( 433 )
  JVS( 434) = W( 10 )
  JVS( 435) = W( 12 )
  JVS( 436) = W( 15 )
  JVS( 437) = W( 16 )
  JVS( 438) = W( 18 )
  JVS( 439) = W( 22 )
  JVS( 440) = W( 25 )
  JVS( 441) = W( 26 )
  JVS( 442) = W( 27 )
  JVS( 443) = W( 31 )
  JVS( 444) = W( 32 )
  JVS( 445) = W( 33 )
  JVS( 446) = W( 35 )
  JVS( 447) = W( 36 )
  JVS( 448) = W( 37 )
  JVS( 449) = W( 39 )
  JVS( 450) = W( 40 )
  JVS( 451) = W( 41 )
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
  IF ( ABS(  JVS( 483 )) < TINY(a) ) THEN
         IER = 53                                      
         RETURN
  END IF
   W( 10 ) = JVS( 467 )
   W( 30 ) = JVS( 468 )
   W( 34 ) = JVS( 469 )
   W( 35 ) = JVS( 470 )
   W( 38 ) = JVS( 471 )
   W( 39 ) = JVS( 472 )
   W( 42 ) = JVS( 473 )
   W( 43 ) = JVS( 474 )
   W( 45 ) = JVS( 475 )
   W( 46 ) = JVS( 476 )
   W( 47 ) = JVS( 477 )
   W( 48 ) = JVS( 478 )
   W( 49 ) = JVS( 479 )
   W( 50 ) = JVS( 480 )
   W( 51 ) = JVS( 481 )
   W( 52 ) = JVS( 482 )
   W( 53 ) = JVS( 483 )
   W( 54 ) = JVS( 484 )
   W( 55 ) = JVS( 485 )
   W( 56 ) = JVS( 486 )
  a = -W( 10 ) / JVS(           27  )
  W( 10 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 28 )
  W( 53 ) = W( 53 ) + a*JVS( 29 )
  a = -W( 30 ) / JVS(          152  )
  W( 30 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 153 )
  W( 51 ) = W( 51 ) + a*JVS( 154 )
  W( 54 ) = W( 54 ) + a*JVS( 155 )
  W( 56 ) = W( 56 ) + a*JVS( 156 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 173 )
  W( 51 ) = W( 51 ) + a*JVS( 174 )
  W( 54 ) = W( 54 ) + a*JVS( 175 )
  a = -W( 35 ) / JVS(          176  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 177 )
  W( 51 ) = W( 51 ) + a*JVS( 178 )
  W( 54 ) = W( 54 ) + a*JVS( 179 )
  a = -W( 38 ) / JVS(          219  )
  W( 38 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 220 )
  W( 45 ) = W( 45 ) + a*JVS( 221 )
  W( 46 ) = W( 46 ) + a*JVS( 222 )
  W( 47 ) = W( 47 ) + a*JVS( 223 )
  W( 48 ) = W( 48 ) + a*JVS( 224 )
  W( 49 ) = W( 49 ) + a*JVS( 225 )
  W( 50 ) = W( 50 ) + a*JVS( 226 )
  W( 51 ) = W( 51 ) + a*JVS( 227 )
  W( 54 ) = W( 54 ) + a*JVS( 228 )
  W( 55 ) = W( 55 ) + a*JVS( 229 )
  W( 56 ) = W( 56 ) + a*JVS( 230 )
  a = -W( 39 ) / JVS(          231  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 54 ) = W( 54 ) + a*JVS( 234 )
  a = -W( 42 ) / JVS(          262  )
  W( 42 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 263 )
  W( 46 ) = W( 46 ) + a*JVS( 264 )
  W( 47 ) = W( 47 ) + a*JVS( 265 )
  W( 50 ) = W( 50 ) + a*JVS( 266 )
  W( 51 ) = W( 51 ) + a*JVS( 267 )
  W( 54 ) = W( 54 ) + a*JVS( 268 )
  W( 55 ) = W( 55 ) + a*JVS( 269 )
  W( 56 ) = W( 56 ) + a*JVS( 270 )
  a = -W( 43 ) / JVS(          282  )
  W( 43 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 283 )
  W( 47 ) = W( 47 ) + a*JVS( 284 )
  W( 48 ) = W( 48 ) + a*JVS( 285 )
  W( 49 ) = W( 49 ) + a*JVS( 286 )
  W( 50 ) = W( 50 ) + a*JVS( 287 )
  W( 51 ) = W( 51 ) + a*JVS( 288 )
  W( 54 ) = W( 54 ) + a*JVS( 289 )
  W( 55 ) = W( 55 ) + a*JVS( 290 )
  W( 56 ) = W( 56 ) + a*JVS( 291 )
  a = -W( 45 ) / JVS(          313  )
  W( 45 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 314 )
  W( 51 ) = W( 51 ) + a*JVS( 315 )
  W( 54 ) = W( 54 ) + a*JVS( 316 )
  W( 55 ) = W( 55 ) + a*JVS( 317 )
  W( 56 ) = W( 56 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          322  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 323 )
  W( 48 ) = W( 48 ) + a*JVS( 324 )
  W( 49 ) = W( 49 ) + a*JVS( 325 )
  W( 50 ) = W( 50 ) + a*JVS( 326 )
  W( 51 ) = W( 51 ) + a*JVS( 327 )
  W( 54 ) = W( 54 ) + a*JVS( 328 )
  W( 55 ) = W( 55 ) + a*JVS( 329 )
  W( 56 ) = W( 56 ) + a*JVS( 330 )
  a = -W( 47 ) / JVS(          336  )
  W( 47 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 337 )
  W( 49 ) = W( 49 ) + a*JVS( 338 )
  W( 50 ) = W( 50 ) + a*JVS( 339 )
  W( 51 ) = W( 51 ) + a*JVS( 340 )
  W( 54 ) = W( 54 ) + a*JVS( 341 )
  W( 55 ) = W( 55 ) + a*JVS( 342 )
  W( 56 ) = W( 56 ) + a*JVS( 343 )
  a = -W( 48 ) / JVS(          353  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 354 )
  W( 50 ) = W( 50 ) + a*JVS( 355 )
  W( 51 ) = W( 51 ) + a*JVS( 356 )
  W( 54 ) = W( 54 ) + a*JVS( 357 )
  W( 55 ) = W( 55 ) + a*JVS( 358 )
  W( 56 ) = W( 56 ) + a*JVS( 359 )
  a = -W( 49 ) / JVS(          369  )
  W( 49 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 370 )
  W( 51 ) = W( 51 ) + a*JVS( 371 )
  W( 52 ) = W( 52 ) + a*JVS( 372 )
  W( 54 ) = W( 54 ) + a*JVS( 373 )
  W( 55 ) = W( 55 ) + a*JVS( 374 )
  W( 56 ) = W( 56 ) + a*JVS( 375 )
  a = -W( 50 ) / JVS(          413  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 414 )
  W( 52 ) = W( 52 ) + a*JVS( 415 )
  W( 53 ) = W( 53 ) + a*JVS( 416 )
  W( 54 ) = W( 54 ) + a*JVS( 417 )
  W( 55 ) = W( 55 ) + a*JVS( 418 )
  W( 56 ) = W( 56 ) + a*JVS( 419 )
  a = -W( 51 ) / JVS(          428  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 429 )
  W( 53 ) = W( 53 ) + a*JVS( 430 )
  W( 54 ) = W( 54 ) + a*JVS( 431 )
  W( 55 ) = W( 55 ) + a*JVS( 432 )
  W( 56 ) = W( 56 ) + a*JVS( 433 )
  a = -W( 52 ) / JVS(          462  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 463 )
  W( 54 ) = W( 54 ) + a*JVS( 464 )
  W( 55 ) = W( 55 ) + a*JVS( 465 )
  W( 56 ) = W( 56 ) + a*JVS( 466 )
  JVS( 467) = W( 10 )
  JVS( 468) = W( 30 )
  JVS( 469) = W( 34 )
  JVS( 470) = W( 35 )
  JVS( 471) = W( 38 )
  JVS( 472) = W( 39 )
  JVS( 473) = W( 42 )
  JVS( 474) = W( 43 )
  JVS( 475) = W( 45 )
  JVS( 476) = W( 46 )
  JVS( 477) = W( 47 )
  JVS( 478) = W( 48 )
  JVS( 479) = W( 49 )
  JVS( 480) = W( 50 )
  JVS( 481) = W( 51 )
  JVS( 482) = W( 52 )
  JVS( 483) = W( 53 )
  JVS( 484) = W( 54 )
  JVS( 485) = W( 55 )
  JVS( 486) = W( 56 )
  IF ( ABS(  JVS( 510 )) < TINY(a) ) THEN
         IER = 54                                      
         RETURN
  END IF
   W( 12 ) = JVS( 487 )
   W( 25 ) = JVS( 488 )
   W( 26 ) = JVS( 489 )
   W( 27 ) = JVS( 490 )
   W( 34 ) = JVS( 491 )
   W( 35 ) = JVS( 492 )
   W( 36 ) = JVS( 493 )
   W( 37 ) = JVS( 494 )
   W( 39 ) = JVS( 495 )
   W( 40 ) = JVS( 496 )
   W( 41 ) = JVS( 497 )
   W( 42 ) = JVS( 498 )
   W( 43 ) = JVS( 499 )
   W( 44 ) = JVS( 500 )
   W( 45 ) = JVS( 501 )
   W( 46 ) = JVS( 502 )
   W( 47 ) = JVS( 503 )
   W( 48 ) = JVS( 504 )
   W( 49 ) = JVS( 505 )
   W( 50 ) = JVS( 506 )
   W( 51 ) = JVS( 507 )
   W( 52 ) = JVS( 508 )
   W( 53 ) = JVS( 509 )
   W( 54 ) = JVS( 510 )
   W( 55 ) = JVS( 511 )
   W( 56 ) = JVS( 512 )
  a = -W( 12 ) / JVS(           32  )
  W( 12 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 33 )
  W( 54 ) = W( 54 ) + a*JVS( 34 )
  a = -W( 25 ) / JVS(           97  )
  W( 25 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 98 )
  W( 52 ) = W( 52 ) + a*JVS( 99 )
  W( 54 ) = W( 54 ) + a*JVS( 100 )
  W( 56 ) = W( 56 ) + a*JVS( 101 )
  a = -W( 26 ) / JVS(          104  )
  W( 26 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 105 )
  W( 54 ) = W( 54 ) + a*JVS( 106 )
  a = -W( 27 ) / JVS(          109  )
  W( 27 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 110 )
  W( 42 ) = W( 42 ) + a*JVS( 111 )
  W( 43 ) = W( 43 ) + a*JVS( 112 )
  W( 45 ) = W( 45 ) + a*JVS( 113 )
  W( 50 ) = W( 50 ) + a*JVS( 114 )
  W( 52 ) = W( 52 ) + a*JVS( 115 )
  W( 54 ) = W( 54 ) + a*JVS( 116 )
  W( 55 ) = W( 55 ) + a*JVS( 117 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 173 )
  W( 51 ) = W( 51 ) + a*JVS( 174 )
  W( 54 ) = W( 54 ) + a*JVS( 175 )
  a = -W( 35 ) / JVS(          176  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 177 )
  W( 51 ) = W( 51 ) + a*JVS( 178 )
  W( 54 ) = W( 54 ) + a*JVS( 179 )
  a = -W( 36 ) / JVS(          188  )
  W( 36 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 189 )
  W( 41 ) = W( 41 ) + a*JVS( 190 )
  W( 44 ) = W( 44 ) + a*JVS( 191 )
  W( 45 ) = W( 45 ) + a*JVS( 192 )
  W( 46 ) = W( 46 ) + a*JVS( 193 )
  W( 50 ) = W( 50 ) + a*JVS( 194 )
  W( 51 ) = W( 51 ) + a*JVS( 195 )
  W( 54 ) = W( 54 ) + a*JVS( 196 )
  W( 55 ) = W( 55 ) + a*JVS( 197 )
  W( 56 ) = W( 56 ) + a*JVS( 198 )
  a = -W( 37 ) / JVS(          206  )
  W( 37 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 207 )
  W( 42 ) = W( 42 ) + a*JVS( 208 )
  W( 45 ) = W( 45 ) + a*JVS( 209 )
  W( 47 ) = W( 47 ) + a*JVS( 210 )
  W( 48 ) = W( 48 ) + a*JVS( 211 )
  W( 49 ) = W( 49 ) + a*JVS( 212 )
  W( 50 ) = W( 50 ) + a*JVS( 213 )
  W( 51 ) = W( 51 ) + a*JVS( 214 )
  W( 54 ) = W( 54 ) + a*JVS( 215 )
  W( 55 ) = W( 55 ) + a*JVS( 216 )
  W( 56 ) = W( 56 ) + a*JVS( 217 )
  a = -W( 39 ) / JVS(          231  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 54 ) = W( 54 ) + a*JVS( 234 )
  a = -W( 40 ) / JVS(          239  )
  W( 40 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 240 )
  W( 48 ) = W( 48 ) + a*JVS( 241 )
  W( 49 ) = W( 49 ) + a*JVS( 242 )
  W( 50 ) = W( 50 ) + a*JVS( 243 )
  W( 51 ) = W( 51 ) + a*JVS( 244 )
  W( 54 ) = W( 54 ) + a*JVS( 245 )
  W( 55 ) = W( 55 ) + a*JVS( 246 )
  W( 56 ) = W( 56 ) + a*JVS( 247 )
  a = -W( 41 ) / JVS(          250  )
  W( 41 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 251 )
  W( 50 ) = W( 50 ) + a*JVS( 252 )
  W( 51 ) = W( 51 ) + a*JVS( 253 )
  W( 54 ) = W( 54 ) + a*JVS( 254 )
  W( 55 ) = W( 55 ) + a*JVS( 255 )
  W( 56 ) = W( 56 ) + a*JVS( 256 )
  a = -W( 42 ) / JVS(          262  )
  W( 42 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 263 )
  W( 46 ) = W( 46 ) + a*JVS( 264 )
  W( 47 ) = W( 47 ) + a*JVS( 265 )
  W( 50 ) = W( 50 ) + a*JVS( 266 )
  W( 51 ) = W( 51 ) + a*JVS( 267 )
  W( 54 ) = W( 54 ) + a*JVS( 268 )
  W( 55 ) = W( 55 ) + a*JVS( 269 )
  W( 56 ) = W( 56 ) + a*JVS( 270 )
  a = -W( 43 ) / JVS(          282  )
  W( 43 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 283 )
  W( 47 ) = W( 47 ) + a*JVS( 284 )
  W( 48 ) = W( 48 ) + a*JVS( 285 )
  W( 49 ) = W( 49 ) + a*JVS( 286 )
  W( 50 ) = W( 50 ) + a*JVS( 287 )
  W( 51 ) = W( 51 ) + a*JVS( 288 )
  W( 54 ) = W( 54 ) + a*JVS( 289 )
  W( 55 ) = W( 55 ) + a*JVS( 290 )
  W( 56 ) = W( 56 ) + a*JVS( 291 )
  a = -W( 44 ) / JVS(          298  )
  W( 44 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 299 )
  W( 46 ) = W( 46 ) + a*JVS( 300 )
  W( 47 ) = W( 47 ) + a*JVS( 301 )
  W( 48 ) = W( 48 ) + a*JVS( 302 )
  W( 49 ) = W( 49 ) + a*JVS( 303 )
  W( 50 ) = W( 50 ) + a*JVS( 304 )
  W( 51 ) = W( 51 ) + a*JVS( 305 )
  W( 53 ) = W( 53 ) + a*JVS( 306 )
  W( 54 ) = W( 54 ) + a*JVS( 307 )
  W( 55 ) = W( 55 ) + a*JVS( 308 )
  W( 56 ) = W( 56 ) + a*JVS( 309 )
  a = -W( 45 ) / JVS(          313  )
  W( 45 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 314 )
  W( 51 ) = W( 51 ) + a*JVS( 315 )
  W( 54 ) = W( 54 ) + a*JVS( 316 )
  W( 55 ) = W( 55 ) + a*JVS( 317 )
  W( 56 ) = W( 56 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          322  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 323 )
  W( 48 ) = W( 48 ) + a*JVS( 324 )
  W( 49 ) = W( 49 ) + a*JVS( 325 )
  W( 50 ) = W( 50 ) + a*JVS( 326 )
  W( 51 ) = W( 51 ) + a*JVS( 327 )
  W( 54 ) = W( 54 ) + a*JVS( 328 )
  W( 55 ) = W( 55 ) + a*JVS( 329 )
  W( 56 ) = W( 56 ) + a*JVS( 330 )
  a = -W( 47 ) / JVS(          336  )
  W( 47 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 337 )
  W( 49 ) = W( 49 ) + a*JVS( 338 )
  W( 50 ) = W( 50 ) + a*JVS( 339 )
  W( 51 ) = W( 51 ) + a*JVS( 340 )
  W( 54 ) = W( 54 ) + a*JVS( 341 )
  W( 55 ) = W( 55 ) + a*JVS( 342 )
  W( 56 ) = W( 56 ) + a*JVS( 343 )
  a = -W( 48 ) / JVS(          353  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 354 )
  W( 50 ) = W( 50 ) + a*JVS( 355 )
  W( 51 ) = W( 51 ) + a*JVS( 356 )
  W( 54 ) = W( 54 ) + a*JVS( 357 )
  W( 55 ) = W( 55 ) + a*JVS( 358 )
  W( 56 ) = W( 56 ) + a*JVS( 359 )
  a = -W( 49 ) / JVS(          369  )
  W( 49 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 370 )
  W( 51 ) = W( 51 ) + a*JVS( 371 )
  W( 52 ) = W( 52 ) + a*JVS( 372 )
  W( 54 ) = W( 54 ) + a*JVS( 373 )
  W( 55 ) = W( 55 ) + a*JVS( 374 )
  W( 56 ) = W( 56 ) + a*JVS( 375 )
  a = -W( 50 ) / JVS(          413  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 414 )
  W( 52 ) = W( 52 ) + a*JVS( 415 )
  W( 53 ) = W( 53 ) + a*JVS( 416 )
  W( 54 ) = W( 54 ) + a*JVS( 417 )
  W( 55 ) = W( 55 ) + a*JVS( 418 )
  W( 56 ) = W( 56 ) + a*JVS( 419 )
  a = -W( 51 ) / JVS(          428  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 429 )
  W( 53 ) = W( 53 ) + a*JVS( 430 )
  W( 54 ) = W( 54 ) + a*JVS( 431 )
  W( 55 ) = W( 55 ) + a*JVS( 432 )
  W( 56 ) = W( 56 ) + a*JVS( 433 )
  a = -W( 52 ) / JVS(          462  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 463 )
  W( 54 ) = W( 54 ) + a*JVS( 464 )
  W( 55 ) = W( 55 ) + a*JVS( 465 )
  W( 56 ) = W( 56 ) + a*JVS( 466 )
  a = -W( 53 ) / JVS(          483  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 484 )
  W( 55 ) = W( 55 ) + a*JVS( 485 )
  W( 56 ) = W( 56 ) + a*JVS( 486 )
  JVS( 487) = W( 12 )
  JVS( 488) = W( 25 )
  JVS( 489) = W( 26 )
  JVS( 490) = W( 27 )
  JVS( 491) = W( 34 )
  JVS( 492) = W( 35 )
  JVS( 493) = W( 36 )
  JVS( 494) = W( 37 )
  JVS( 495) = W( 39 )
  JVS( 496) = W( 40 )
  JVS( 497) = W( 41 )
  JVS( 498) = W( 42 )
  JVS( 499) = W( 43 )
  JVS( 500) = W( 44 )
  JVS( 501) = W( 45 )
  JVS( 502) = W( 46 )
  JVS( 503) = W( 47 )
  JVS( 504) = W( 48 )
  JVS( 505) = W( 49 )
  JVS( 506) = W( 50 )
  JVS( 507) = W( 51 )
  JVS( 508) = W( 52 )
  JVS( 509) = W( 53 )
  JVS( 510) = W( 54 )
  JVS( 511) = W( 55 )
  JVS( 512) = W( 56 )
  IF ( ABS(  JVS( 550 )) < TINY(a) ) THEN
         IER = 55                                      
         RETURN
  END IF
   W( 7 ) = JVS( 513 )
   W( 8 ) = JVS( 514 )
   W( 9 ) = JVS( 515 )
   W( 11 ) = JVS( 516 )
   W( 13 ) = JVS( 517 )
   W( 16 ) = JVS( 518 )
   W( 18 ) = JVS( 519 )
   W( 20 ) = JVS( 520 )
   W( 21 ) = JVS( 521 )
   W( 23 ) = JVS( 522 )
   W( 24 ) = JVS( 523 )
   W( 26 ) = JVS( 524 )
   W( 28 ) = JVS( 525 )
   W( 30 ) = JVS( 526 )
   W( 31 ) = JVS( 527 )
   W( 32 ) = JVS( 528 )
   W( 33 ) = JVS( 529 )
   W( 34 ) = JVS( 530 )
   W( 35 ) = JVS( 531 )
   W( 36 ) = JVS( 532 )
   W( 37 ) = JVS( 533 )
   W( 39 ) = JVS( 534 )
   W( 40 ) = JVS( 535 )
   W( 41 ) = JVS( 536 )
   W( 42 ) = JVS( 537 )
   W( 43 ) = JVS( 538 )
   W( 44 ) = JVS( 539 )
   W( 45 ) = JVS( 540 )
   W( 46 ) = JVS( 541 )
   W( 47 ) = JVS( 542 )
   W( 48 ) = JVS( 543 )
   W( 49 ) = JVS( 544 )
   W( 50 ) = JVS( 545 )
   W( 51 ) = JVS( 546 )
   W( 52 ) = JVS( 547 )
   W( 53 ) = JVS( 548 )
   W( 54 ) = JVS( 549 )
   W( 55 ) = JVS( 550 )
   W( 56 ) = JVS( 551 )
  a = -W( 7 ) / JVS(           20  )
  W( 7 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 21 )
  a = -W( 8 ) / JVS(           22  )
  W( 8 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 23 )
  a = -W( 9 ) / JVS(           24  )
  W( 9 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 25 )
  W( 55 ) = W( 55 ) + a*JVS( 26 )
  a = -W( 11 ) / JVS(           30  )
  W( 11 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 31 )
  a = -W( 13 ) / JVS(           35  )
  W( 13 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 36 )
  a = -W( 16 ) / JVS(           47  )
  W( 16 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 48 )
  W( 52 ) = W( 52 ) + a*JVS( 49 )
  W( 55 ) = W( 55 ) + a*JVS( 50 )
  a = -W( 18 ) / JVS(           59  )
  W( 18 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 60 )
  W( 56 ) = W( 56 ) + a*JVS( 61 )
  a = -W( 20 ) / JVS(           74  )
  W( 20 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 75 )
  W( 50 ) = W( 50 ) + a*JVS( 76 )
  W( 55 ) = W( 55 ) + a*JVS( 77 )
  a = -W( 21 ) / JVS(           78  )
  W( 21 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 79 )
  W( 50 ) = W( 50 ) + a*JVS( 80 )
  W( 55 ) = W( 55 ) + a*JVS( 81 )
  a = -W( 23 ) / JVS(           87  )
  W( 23 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 88 )
  W( 51 ) = W( 51 ) + a*JVS( 89 )
  a = -W( 24 ) / JVS(           90  )
  W( 24 ) = -a
  W( 34 ) = W( 34 ) + a*JVS( 91 )
  W( 39 ) = W( 39 ) + a*JVS( 92 )
  W( 44 ) = W( 44 ) + a*JVS( 93 )
  W( 50 ) = W( 50 ) + a*JVS( 94 )
  W( 51 ) = W( 51 ) + a*JVS( 95 )
  a = -W( 26 ) / JVS(          104  )
  W( 26 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 105 )
  W( 54 ) = W( 54 ) + a*JVS( 106 )
  a = -W( 28 ) / JVS(          119  )
  W( 28 ) = -a
  W( 30 ) = W( 30 ) + a*JVS( 120 )
  W( 33 ) = W( 33 ) + a*JVS( 121 )
  W( 34 ) = W( 34 ) + a*JVS( 122 )
  W( 35 ) = W( 35 ) + a*JVS( 123 )
  W( 36 ) = W( 36 ) + a*JVS( 124 )
  W( 39 ) = W( 39 ) + a*JVS( 125 )
  W( 42 ) = W( 42 ) + a*JVS( 126 )
  W( 43 ) = W( 43 ) + a*JVS( 127 )
  W( 45 ) = W( 45 ) + a*JVS( 128 )
  W( 50 ) = W( 50 ) + a*JVS( 129 )
  W( 51 ) = W( 51 ) + a*JVS( 130 )
  W( 54 ) = W( 54 ) + a*JVS( 131 )
  W( 56 ) = W( 56 ) + a*JVS( 132 )
  a = -W( 30 ) / JVS(          152  )
  W( 30 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 153 )
  W( 51 ) = W( 51 ) + a*JVS( 154 )
  W( 54 ) = W( 54 ) + a*JVS( 155 )
  W( 56 ) = W( 56 ) + a*JVS( 156 )
  a = -W( 31 ) / JVS(          157  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 158 )
  W( 54 ) = W( 54 ) + a*JVS( 159 )
  W( 55 ) = W( 55 ) + a*JVS( 160 )
  W( 56 ) = W( 56 ) + a*JVS( 161 )
  a = -W( 32 ) / JVS(          162  )
  W( 32 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 163 )
  W( 50 ) = W( 50 ) + a*JVS( 164 )
  W( 55 ) = W( 55 ) + a*JVS( 165 )
  W( 56 ) = W( 56 ) + a*JVS( 166 )
  a = -W( 33 ) / JVS(          167  )
  W( 33 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 168 )
  W( 50 ) = W( 50 ) + a*JVS( 169 )
  W( 55 ) = W( 55 ) + a*JVS( 170 )
  W( 56 ) = W( 56 ) + a*JVS( 171 )
  a = -W( 34 ) / JVS(          172  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 173 )
  W( 51 ) = W( 51 ) + a*JVS( 174 )
  W( 54 ) = W( 54 ) + a*JVS( 175 )
  a = -W( 35 ) / JVS(          176  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 177 )
  W( 51 ) = W( 51 ) + a*JVS( 178 )
  W( 54 ) = W( 54 ) + a*JVS( 179 )
  a = -W( 36 ) / JVS(          188  )
  W( 36 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 189 )
  W( 41 ) = W( 41 ) + a*JVS( 190 )
  W( 44 ) = W( 44 ) + a*JVS( 191 )
  W( 45 ) = W( 45 ) + a*JVS( 192 )
  W( 46 ) = W( 46 ) + a*JVS( 193 )
  W( 50 ) = W( 50 ) + a*JVS( 194 )
  W( 51 ) = W( 51 ) + a*JVS( 195 )
  W( 54 ) = W( 54 ) + a*JVS( 196 )
  W( 55 ) = W( 55 ) + a*JVS( 197 )
  W( 56 ) = W( 56 ) + a*JVS( 198 )
  a = -W( 37 ) / JVS(          206  )
  W( 37 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 207 )
  W( 42 ) = W( 42 ) + a*JVS( 208 )
  W( 45 ) = W( 45 ) + a*JVS( 209 )
  W( 47 ) = W( 47 ) + a*JVS( 210 )
  W( 48 ) = W( 48 ) + a*JVS( 211 )
  W( 49 ) = W( 49 ) + a*JVS( 212 )
  W( 50 ) = W( 50 ) + a*JVS( 213 )
  W( 51 ) = W( 51 ) + a*JVS( 214 )
  W( 54 ) = W( 54 ) + a*JVS( 215 )
  W( 55 ) = W( 55 ) + a*JVS( 216 )
  W( 56 ) = W( 56 ) + a*JVS( 217 )
  a = -W( 39 ) / JVS(          231  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 54 ) = W( 54 ) + a*JVS( 234 )
  a = -W( 40 ) / JVS(          239  )
  W( 40 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 240 )
  W( 48 ) = W( 48 ) + a*JVS( 241 )
  W( 49 ) = W( 49 ) + a*JVS( 242 )
  W( 50 ) = W( 50 ) + a*JVS( 243 )
  W( 51 ) = W( 51 ) + a*JVS( 244 )
  W( 54 ) = W( 54 ) + a*JVS( 245 )
  W( 55 ) = W( 55 ) + a*JVS( 246 )
  W( 56 ) = W( 56 ) + a*JVS( 247 )
  a = -W( 41 ) / JVS(          250  )
  W( 41 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 251 )
  W( 50 ) = W( 50 ) + a*JVS( 252 )
  W( 51 ) = W( 51 ) + a*JVS( 253 )
  W( 54 ) = W( 54 ) + a*JVS( 254 )
  W( 55 ) = W( 55 ) + a*JVS( 255 )
  W( 56 ) = W( 56 ) + a*JVS( 256 )
  a = -W( 42 ) / JVS(          262  )
  W( 42 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 263 )
  W( 46 ) = W( 46 ) + a*JVS( 264 )
  W( 47 ) = W( 47 ) + a*JVS( 265 )
  W( 50 ) = W( 50 ) + a*JVS( 266 )
  W( 51 ) = W( 51 ) + a*JVS( 267 )
  W( 54 ) = W( 54 ) + a*JVS( 268 )
  W( 55 ) = W( 55 ) + a*JVS( 269 )
  W( 56 ) = W( 56 ) + a*JVS( 270 )
  a = -W( 43 ) / JVS(          282  )
  W( 43 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 283 )
  W( 47 ) = W( 47 ) + a*JVS( 284 )
  W( 48 ) = W( 48 ) + a*JVS( 285 )
  W( 49 ) = W( 49 ) + a*JVS( 286 )
  W( 50 ) = W( 50 ) + a*JVS( 287 )
  W( 51 ) = W( 51 ) + a*JVS( 288 )
  W( 54 ) = W( 54 ) + a*JVS( 289 )
  W( 55 ) = W( 55 ) + a*JVS( 290 )
  W( 56 ) = W( 56 ) + a*JVS( 291 )
  a = -W( 44 ) / JVS(          298  )
  W( 44 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 299 )
  W( 46 ) = W( 46 ) + a*JVS( 300 )
  W( 47 ) = W( 47 ) + a*JVS( 301 )
  W( 48 ) = W( 48 ) + a*JVS( 302 )
  W( 49 ) = W( 49 ) + a*JVS( 303 )
  W( 50 ) = W( 50 ) + a*JVS( 304 )
  W( 51 ) = W( 51 ) + a*JVS( 305 )
  W( 53 ) = W( 53 ) + a*JVS( 306 )
  W( 54 ) = W( 54 ) + a*JVS( 307 )
  W( 55 ) = W( 55 ) + a*JVS( 308 )
  W( 56 ) = W( 56 ) + a*JVS( 309 )
  a = -W( 45 ) / JVS(          313  )
  W( 45 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 314 )
  W( 51 ) = W( 51 ) + a*JVS( 315 )
  W( 54 ) = W( 54 ) + a*JVS( 316 )
  W( 55 ) = W( 55 ) + a*JVS( 317 )
  W( 56 ) = W( 56 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          322  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 323 )
  W( 48 ) = W( 48 ) + a*JVS( 324 )
  W( 49 ) = W( 49 ) + a*JVS( 325 )
  W( 50 ) = W( 50 ) + a*JVS( 326 )
  W( 51 ) = W( 51 ) + a*JVS( 327 )
  W( 54 ) = W( 54 ) + a*JVS( 328 )
  W( 55 ) = W( 55 ) + a*JVS( 329 )
  W( 56 ) = W( 56 ) + a*JVS( 330 )
  a = -W( 47 ) / JVS(          336  )
  W( 47 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 337 )
  W( 49 ) = W( 49 ) + a*JVS( 338 )
  W( 50 ) = W( 50 ) + a*JVS( 339 )
  W( 51 ) = W( 51 ) + a*JVS( 340 )
  W( 54 ) = W( 54 ) + a*JVS( 341 )
  W( 55 ) = W( 55 ) + a*JVS( 342 )
  W( 56 ) = W( 56 ) + a*JVS( 343 )
  a = -W( 48 ) / JVS(          353  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 354 )
  W( 50 ) = W( 50 ) + a*JVS( 355 )
  W( 51 ) = W( 51 ) + a*JVS( 356 )
  W( 54 ) = W( 54 ) + a*JVS( 357 )
  W( 55 ) = W( 55 ) + a*JVS( 358 )
  W( 56 ) = W( 56 ) + a*JVS( 359 )
  a = -W( 49 ) / JVS(          369  )
  W( 49 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 370 )
  W( 51 ) = W( 51 ) + a*JVS( 371 )
  W( 52 ) = W( 52 ) + a*JVS( 372 )
  W( 54 ) = W( 54 ) + a*JVS( 373 )
  W( 55 ) = W( 55 ) + a*JVS( 374 )
  W( 56 ) = W( 56 ) + a*JVS( 375 )
  a = -W( 50 ) / JVS(          413  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 414 )
  W( 52 ) = W( 52 ) + a*JVS( 415 )
  W( 53 ) = W( 53 ) + a*JVS( 416 )
  W( 54 ) = W( 54 ) + a*JVS( 417 )
  W( 55 ) = W( 55 ) + a*JVS( 418 )
  W( 56 ) = W( 56 ) + a*JVS( 419 )
  a = -W( 51 ) / JVS(          428  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 429 )
  W( 53 ) = W( 53 ) + a*JVS( 430 )
  W( 54 ) = W( 54 ) + a*JVS( 431 )
  W( 55 ) = W( 55 ) + a*JVS( 432 )
  W( 56 ) = W( 56 ) + a*JVS( 433 )
  a = -W( 52 ) / JVS(          462  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 463 )
  W( 54 ) = W( 54 ) + a*JVS( 464 )
  W( 55 ) = W( 55 ) + a*JVS( 465 )
  W( 56 ) = W( 56 ) + a*JVS( 466 )
  a = -W( 53 ) / JVS(          483  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 484 )
  W( 55 ) = W( 55 ) + a*JVS( 485 )
  W( 56 ) = W( 56 ) + a*JVS( 486 )
  a = -W( 54 ) / JVS(          510  )
  W( 54 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 511 )
  W( 56 ) = W( 56 ) + a*JVS( 512 )
  JVS( 513) = W( 7 )
  JVS( 514) = W( 8 )
  JVS( 515) = W( 9 )
  JVS( 516) = W( 11 )
  JVS( 517) = W( 13 )
  JVS( 518) = W( 16 )
  JVS( 519) = W( 18 )
  JVS( 520) = W( 20 )
  JVS( 521) = W( 21 )
  JVS( 522) = W( 23 )
  JVS( 523) = W( 24 )
  JVS( 524) = W( 26 )
  JVS( 525) = W( 28 )
  JVS( 526) = W( 30 )
  JVS( 527) = W( 31 )
  JVS( 528) = W( 32 )
  JVS( 529) = W( 33 )
  JVS( 530) = W( 34 )
  JVS( 531) = W( 35 )
  JVS( 532) = W( 36 )
  JVS( 533) = W( 37 )
  JVS( 534) = W( 39 )
  JVS( 535) = W( 40 )
  JVS( 536) = W( 41 )
  JVS( 537) = W( 42 )
  JVS( 538) = W( 43 )
  JVS( 539) = W( 44 )
  JVS( 540) = W( 45 )
  JVS( 541) = W( 46 )
  JVS( 542) = W( 47 )
  JVS( 543) = W( 48 )
  JVS( 544) = W( 49 )
  JVS( 545) = W( 50 )
  JVS( 546) = W( 51 )
  JVS( 547) = W( 52 )
  JVS( 548) = W( 53 )
  JVS( 549) = W( 54 )
  JVS( 550) = W( 55 )
  JVS( 551) = W( 56 )
  IF ( ABS(  JVS( 576 )) < TINY(a) ) THEN
         IER = 56                                      
         RETURN
  END IF
   W( 18 ) = JVS( 552 )
   W( 22 ) = JVS( 553 )
   W( 25 ) = JVS( 554 )
   W( 31 ) = JVS( 555 )
   W( 32 ) = JVS( 556 )
   W( 33 ) = JVS( 557 )
   W( 35 ) = JVS( 558 )
   W( 37 ) = JVS( 559 )
   W( 39 ) = JVS( 560 )
   W( 40 ) = JVS( 561 )
   W( 41 ) = JVS( 562 )
   W( 42 ) = JVS( 563 )
   W( 44 ) = JVS( 564 )
   W( 45 ) = JVS( 565 )
   W( 46 ) = JVS( 566 )
   W( 47 ) = JVS( 567 )
   W( 48 ) = JVS( 568 )
   W( 49 ) = JVS( 569 )
   W( 50 ) = JVS( 570 )
   W( 51 ) = JVS( 571 )
   W( 52 ) = JVS( 572 )
   W( 53 ) = JVS( 573 )
   W( 54 ) = JVS( 574 )
   W( 55 ) = JVS( 575 )
   W( 56 ) = JVS( 576 )
  a = -W( 18 ) / JVS(           59  )
  W( 18 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 60 )
  W( 56 ) = W( 56 ) + a*JVS( 61 )
  a = -W( 22 ) / JVS(           82  )
  W( 22 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 83 )
  W( 52 ) = W( 52 ) + a*JVS( 84 )
  W( 55 ) = W( 55 ) + a*JVS( 85 )
  W( 56 ) = W( 56 ) + a*JVS( 86 )
  a = -W( 25 ) / JVS(           97  )
  W( 25 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 98 )
  W( 52 ) = W( 52 ) + a*JVS( 99 )
  W( 54 ) = W( 54 ) + a*JVS( 100 )
  W( 56 ) = W( 56 ) + a*JVS( 101 )
  a = -W( 31 ) / JVS(          157  )
  W( 31 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 158 )
  W( 54 ) = W( 54 ) + a*JVS( 159 )
  W( 55 ) = W( 55 ) + a*JVS( 160 )
  W( 56 ) = W( 56 ) + a*JVS( 161 )
  a = -W( 32 ) / JVS(          162  )
  W( 32 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 163 )
  W( 50 ) = W( 50 ) + a*JVS( 164 )
  W( 55 ) = W( 55 ) + a*JVS( 165 )
  W( 56 ) = W( 56 ) + a*JVS( 166 )
  a = -W( 33 ) / JVS(          167  )
  W( 33 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 168 )
  W( 50 ) = W( 50 ) + a*JVS( 169 )
  W( 55 ) = W( 55 ) + a*JVS( 170 )
  W( 56 ) = W( 56 ) + a*JVS( 171 )
  a = -W( 35 ) / JVS(          176  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 177 )
  W( 51 ) = W( 51 ) + a*JVS( 178 )
  W( 54 ) = W( 54 ) + a*JVS( 179 )
  a = -W( 37 ) / JVS(          206  )
  W( 37 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 207 )
  W( 42 ) = W( 42 ) + a*JVS( 208 )
  W( 45 ) = W( 45 ) + a*JVS( 209 )
  W( 47 ) = W( 47 ) + a*JVS( 210 )
  W( 48 ) = W( 48 ) + a*JVS( 211 )
  W( 49 ) = W( 49 ) + a*JVS( 212 )
  W( 50 ) = W( 50 ) + a*JVS( 213 )
  W( 51 ) = W( 51 ) + a*JVS( 214 )
  W( 54 ) = W( 54 ) + a*JVS( 215 )
  W( 55 ) = W( 55 ) + a*JVS( 216 )
  W( 56 ) = W( 56 ) + a*JVS( 217 )
  a = -W( 39 ) / JVS(          231  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 232 )
  W( 51 ) = W( 51 ) + a*JVS( 233 )
  W( 54 ) = W( 54 ) + a*JVS( 234 )
  a = -W( 40 ) / JVS(          239  )
  W( 40 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 240 )
  W( 48 ) = W( 48 ) + a*JVS( 241 )
  W( 49 ) = W( 49 ) + a*JVS( 242 )
  W( 50 ) = W( 50 ) + a*JVS( 243 )
  W( 51 ) = W( 51 ) + a*JVS( 244 )
  W( 54 ) = W( 54 ) + a*JVS( 245 )
  W( 55 ) = W( 55 ) + a*JVS( 246 )
  W( 56 ) = W( 56 ) + a*JVS( 247 )
  a = -W( 41 ) / JVS(          250  )
  W( 41 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 251 )
  W( 50 ) = W( 50 ) + a*JVS( 252 )
  W( 51 ) = W( 51 ) + a*JVS( 253 )
  W( 54 ) = W( 54 ) + a*JVS( 254 )
  W( 55 ) = W( 55 ) + a*JVS( 255 )
  W( 56 ) = W( 56 ) + a*JVS( 256 )
  a = -W( 42 ) / JVS(          262  )
  W( 42 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 263 )
  W( 46 ) = W( 46 ) + a*JVS( 264 )
  W( 47 ) = W( 47 ) + a*JVS( 265 )
  W( 50 ) = W( 50 ) + a*JVS( 266 )
  W( 51 ) = W( 51 ) + a*JVS( 267 )
  W( 54 ) = W( 54 ) + a*JVS( 268 )
  W( 55 ) = W( 55 ) + a*JVS( 269 )
  W( 56 ) = W( 56 ) + a*JVS( 270 )
  a = -W( 44 ) / JVS(          298  )
  W( 44 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 299 )
  W( 46 ) = W( 46 ) + a*JVS( 300 )
  W( 47 ) = W( 47 ) + a*JVS( 301 )
  W( 48 ) = W( 48 ) + a*JVS( 302 )
  W( 49 ) = W( 49 ) + a*JVS( 303 )
  W( 50 ) = W( 50 ) + a*JVS( 304 )
  W( 51 ) = W( 51 ) + a*JVS( 305 )
  W( 53 ) = W( 53 ) + a*JVS( 306 )
  W( 54 ) = W( 54 ) + a*JVS( 307 )
  W( 55 ) = W( 55 ) + a*JVS( 308 )
  W( 56 ) = W( 56 ) + a*JVS( 309 )
  a = -W( 45 ) / JVS(          313  )
  W( 45 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 314 )
  W( 51 ) = W( 51 ) + a*JVS( 315 )
  W( 54 ) = W( 54 ) + a*JVS( 316 )
  W( 55 ) = W( 55 ) + a*JVS( 317 )
  W( 56 ) = W( 56 ) + a*JVS( 318 )
  a = -W( 46 ) / JVS(          322  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 323 )
  W( 48 ) = W( 48 ) + a*JVS( 324 )
  W( 49 ) = W( 49 ) + a*JVS( 325 )
  W( 50 ) = W( 50 ) + a*JVS( 326 )
  W( 51 ) = W( 51 ) + a*JVS( 327 )
  W( 54 ) = W( 54 ) + a*JVS( 328 )
  W( 55 ) = W( 55 ) + a*JVS( 329 )
  W( 56 ) = W( 56 ) + a*JVS( 330 )
  a = -W( 47 ) / JVS(          336  )
  W( 47 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 337 )
  W( 49 ) = W( 49 ) + a*JVS( 338 )
  W( 50 ) = W( 50 ) + a*JVS( 339 )
  W( 51 ) = W( 51 ) + a*JVS( 340 )
  W( 54 ) = W( 54 ) + a*JVS( 341 )
  W( 55 ) = W( 55 ) + a*JVS( 342 )
  W( 56 ) = W( 56 ) + a*JVS( 343 )
  a = -W( 48 ) / JVS(          353  )
  W( 48 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 354 )
  W( 50 ) = W( 50 ) + a*JVS( 355 )
  W( 51 ) = W( 51 ) + a*JVS( 356 )
  W( 54 ) = W( 54 ) + a*JVS( 357 )
  W( 55 ) = W( 55 ) + a*JVS( 358 )
  W( 56 ) = W( 56 ) + a*JVS( 359 )
  a = -W( 49 ) / JVS(          369  )
  W( 49 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 370 )
  W( 51 ) = W( 51 ) + a*JVS( 371 )
  W( 52 ) = W( 52 ) + a*JVS( 372 )
  W( 54 ) = W( 54 ) + a*JVS( 373 )
  W( 55 ) = W( 55 ) + a*JVS( 374 )
  W( 56 ) = W( 56 ) + a*JVS( 375 )
  a = -W( 50 ) / JVS(          413  )
  W( 50 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 414 )
  W( 52 ) = W( 52 ) + a*JVS( 415 )
  W( 53 ) = W( 53 ) + a*JVS( 416 )
  W( 54 ) = W( 54 ) + a*JVS( 417 )
  W( 55 ) = W( 55 ) + a*JVS( 418 )
  W( 56 ) = W( 56 ) + a*JVS( 419 )
  a = -W( 51 ) / JVS(          428  )
  W( 51 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 429 )
  W( 53 ) = W( 53 ) + a*JVS( 430 )
  W( 54 ) = W( 54 ) + a*JVS( 431 )
  W( 55 ) = W( 55 ) + a*JVS( 432 )
  W( 56 ) = W( 56 ) + a*JVS( 433 )
  a = -W( 52 ) / JVS(          462  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 463 )
  W( 54 ) = W( 54 ) + a*JVS( 464 )
  W( 55 ) = W( 55 ) + a*JVS( 465 )
  W( 56 ) = W( 56 ) + a*JVS( 466 )
  a = -W( 53 ) / JVS(          483  )
  W( 53 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 484 )
  W( 55 ) = W( 55 ) + a*JVS( 485 )
  W( 56 ) = W( 56 ) + a*JVS( 486 )
  a = -W( 54 ) / JVS(          510  )
  W( 54 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 511 )
  W( 56 ) = W( 56 ) + a*JVS( 512 )
  a = -W( 55 ) / JVS(          550  )
  W( 55 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 551 )
  JVS( 552) = W( 18 )
  JVS( 553) = W( 22 )
  JVS( 554) = W( 25 )
  JVS( 555) = W( 31 )
  JVS( 556) = W( 32 )
  JVS( 557) = W( 33 )
  JVS( 558) = W( 35 )
  JVS( 559) = W( 37 )
  JVS( 560) = W( 39 )
  JVS( 561) = W( 40 )
  JVS( 562) = W( 41 )
  JVS( 563) = W( 42 )
  JVS( 564) = W( 44 )
  JVS( 565) = W( 45 )
  JVS( 566) = W( 46 )
  JVS( 567) = W( 47 )
  JVS( 568) = W( 48 )
  JVS( 569) = W( 49 )
  JVS( 570) = W( 50 )
  JVS( 571) = W( 51 )
  JVS( 572) = W( 52 )
  JVS( 573) = W( 53 )
  JVS( 574) = W( 54 )
  JVS( 575) = W( 55 )
  JVS( 576) = W( 56 )
   
   END SUBROUTINE decomp_cbmz_bb
 


END MODULE cbmz_bb_Integrator
