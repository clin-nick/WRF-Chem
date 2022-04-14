
































MODULE cbm4_Integrator

 USE cbm4_Parameters
 USE cbm4_Precision
 USE cbm4_JacobianSP

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

SUBROUTINE  cbm4_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE cbm4_Parameters

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

   CALL cbm4_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  cbm4_INTEGRATE


SUBROUTINE  cbm4_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE cbm4_Parameters

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
      CALL cbm4_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL cbm4_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = cbm4_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL cbm4_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL cbm4_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL cbm4_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL cbm4_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL cbm4_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL cbm4_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL cbm4_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL cbm4_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL cbm4_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL cbm4_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL cbm4_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL cbm4_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL cbm4_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL cbm4_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL cbm4_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  cbm4_ros_ErrorMsg(Code,T,H,IERR)



   USE cbm4_Precision

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

 END SUBROUTINE  cbm4_ros_ErrorMsg


 SUBROUTINE  cbm4_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL cbm4_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL cbm4_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL cbm4_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL cbm4_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL cbm4_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL cbm4_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL cbm4_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL cbm4_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL cbm4_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL cbm4_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL cbm4_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL cbm4_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL cbm4_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL cbm4_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL cbm4_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL cbm4_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL cbm4_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL cbm4_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL cbm4_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL cbm4_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = cbm4_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL cbm4_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  cbm4_ros_Integrator



  REAL(kind=dp) FUNCTION  cbm4_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    cbm4_ros_ErrorNorm = Err

  END FUNCTION  cbm4_ros_ErrorNorm



  SUBROUTINE cbm4_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL cbm4_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL cbm4_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL cbm4_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  cbm4_ros_FunTimeDeriv



  SUBROUTINE  cbm4_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL cbm4_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL cbm4_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL cbm4_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  cbm4_ros_PrepareMatrix



  SUBROUTINE  cbm4_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 



CALL decomp_cbm4 ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  cbm4_ros_Decomp



  SUBROUTINE  cbm4_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL cbm4_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  cbm4_ros_Solve




  SUBROUTINE  cbm4_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  cbm4_Ros2



  SUBROUTINE  cbm4_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cbm4_Ros3





  SUBROUTINE  cbm4_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cbm4_Ros4


  SUBROUTINE  cbm4_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cbm4_Rodas3


  SUBROUTINE  cbm4_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cbm4_Rodas4




END SUBROUTINE  cbm4_Rosenbrock




SUBROUTINE  cbm4_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE cbm4_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL cbm4_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  cbm4_FunTemplate



SUBROUTINE  cbm4_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE cbm4_Parameters
 
 USE cbm4_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL cbm4_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  cbm4_JacTemplate

















SUBROUTINE cbm4_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(26)
  A(2) = RCT(2)*V(25)
  A(3) = RCT(3)*V(25)
  A(4) = RCT(4)*V(30)
  A(5) = RCT(5)*V(9)
  A(6) = RCT(6)*V(2)
  A(7) = RCT(7)*V(21)
  A(8) = RCT(8)*V(21)
  A(9) = RCT(9)*V(24)
  A(10) = RCT(10)*V(19)
  A(11) = RCT(11)*V(15)
  A(12) = RCT(12)*V(29)
  A(13) = RCT(13)*V(25)*V(31)
  A(14) = RCT(14)*V(26)*V(29)
  A(15) = RCT(15)*V(26)*V(29)
  A(16) = RCT(16)*V(29)*V(31)
  A(17) = RCT(17)*V(25)*V(26)
  A(18) = RCT(18)*V(1)
  A(19) = RCT(19)*V(1)*F(1)
  A(20) = RCT(20)*V(25)*V(27)
  A(21) = RCT(21)*V(25)*V(28)
  A(22) = RCT(22)*V(30)*V(31)
  A(23) = RCT(23)*V(26)*V(30)
  A(24) = RCT(24)*V(26)*V(30)
  A(25) = RCT(25)*V(6)*F(1)
  A(26) = RCT(26)*V(6)
  A(27) = RCT(27)*V(31)*V(31)
  A(28) = RCT(28)*V(26)*V(31)*F(1)
  A(29) = RCT(29)*V(27)*V(31)
  A(30) = RCT(30)*V(9)*V(27)
  A(31) = RCT(31)*V(9)*V(9)
  A(32) = RCT(32)*V(26)*V(27)
  A(33) = RCT(33)*V(12)*V(27)
  A(34) = RCT(34)*V(28)*V(31)
  A(35) = RCT(35)*V(26)*V(28)
  A(36) = RCT(36)*V(10)
  A(37) = RCT(37)*V(10)*V(27)
  A(38) = RCT(38)*V(28)*V(28)
  A(39) = RCT(39)*V(28)*V(28)*F(1)
  A(40) = RCT(40)*V(2)*V(27)
  A(41) = RCT(41)*V(16)*V(27)
  A(42) = RCT(42)*V(21)*V(27)
  A(43) = RCT(43)*V(21)*V(29)
  A(44) = RCT(44)*V(21)*V(30)
  A(45) = RCT(45)*V(24)*V(29)
  A(46) = RCT(46)*V(24)*V(27)
  A(47) = RCT(47)*V(24)*V(30)
  A(48) = RCT(48)*V(31)*V(32)
  A(49) = RCT(49)*V(26)*V(32)
  A(50) = RCT(50)*V(3)
  A(51) = RCT(51)*V(32)*V(32)
  A(52) = RCT(52)*V(28)*V(32)
  A(53) = RCT(53)*V(27)
  A(54) = RCT(54)*V(20)*V(27)
  A(55) = RCT(55)*V(13)
  A(56) = RCT(56)*V(13)
  A(57) = RCT(57)*V(13)*V(26)
  A(58) = RCT(58)*V(23)*V(29)
  A(59) = RCT(59)*V(23)*V(27)
  A(60) = RCT(60)*V(23)*V(25)
  A(61) = RCT(61)*V(23)*V(30)
  A(62) = RCT(62)*V(17)*V(29)
  A(63) = RCT(63)*V(17)*V(27)
  A(64) = RCT(64)*V(17)*V(25)
  A(65) = RCT(65)*V(5)*V(27)
  A(66) = RCT(66)*V(11)*V(31)
  A(67) = 4.2*V(11)
  A(68) = RCT(68)*V(14)*V(27)
  A(69) = RCT(69)*V(14)*V(30)
  A(70) = RCT(70)*V(4)*V(26)
  A(71) = RCT(71)*V(7)*V(27)
  A(72) = RCT(72)*V(19)*V(27)
  A(73) = RCT(73)*V(19)*V(25)
  A(74) = RCT(74)*V(15)*V(27)
  A(75) = RCT(75)*V(22)*V(29)
  A(76) = RCT(76)*V(22)*V(27)
  A(77) = RCT(77)*V(22)*V(25)
  A(78) = RCT(78)*V(22)*V(30)
  A(79) = RCT(79)*V(18)*V(31)
  A(80) = RCT(80)*V(18)*V(18)
  A(81) = RCT(81)*V(8)*V(31)


  Vdot(1) = A(3)-A(18)-A(19)
  Vdot(2) = -A(6)+A(38)+A(39)-A(40)
  Vdot(3) = A(49)-A(50)
  Vdot(4) = 0.4*A(68)+A(69)-A(70)
  Vdot(5) = -A(65)
  Vdot(6) = A(24)-A(25)-A(26)
  Vdot(7) = -A(71)
  Vdot(8) = 0.13*A(54)+0.04*A(55)+0.02*A(58)+0.09*A(61)+0.13*A(76)+A(78)-A(81)
  Vdot(9) = -A(5)+2*A(28)+A(29)-A(30)-2*A(31)
  Vdot(10) = A(35)-A(36)-A(37)
  Vdot(11) = 0.56*A(65)-A(66)-A(67)+0.3*A(71)
  Vdot(12) = 2*A(25)+A(32)-A(33)+A(44)+A(47)+A(69)
  Vdot(13) = 0.76*A(54)-0.98*A(55)-A(56)-A(57)
  Vdot(14) = 0.36*A(65)+A(67)-A(68)-A(69)+0.2*A(71)
  Vdot(15) = -A(11)+0.8*A(71)+0.2*A(73)-A(74)+0.4*A(76)+0.2*A(77)
  Vdot(16) = A(7)+A(8)+A(9)+A(10)+A(11)-A(41)+A(42)+A(43)+A(44)+0.3*A(58)+0.33*A(60)+A(62)+0.42*A(64)+2*A(72)+0.69*A(73)&
               &+0.5*A(75)+0.06*A(77)
  Vdot(17) = -A(62)-A(63)-A(64)+0.45*A(75)+A(76)+0.55*A(77)
  Vdot(18) = A(9)+A(48)+2*A(51)+0.79*A(52)+A(53)+0.87*A(54)+0.96*A(55)+0.28*A(58)+A(59)+0.22*A(60)+0.91*A(61)+0.7*A(62)&
               &+A(63)+0.08*A(65)+0.6*A(68)+0.5*A(71)+A(72)+0.03*A(73)+A(74)+0.5*A(75)+A(76)-A(79)-2*A(80)
  Vdot(19) = -A(10)+0.9*A(66)+0.3*A(68)-A(72)-A(73)
  Vdot(20) = -1.11*A(54)-2.1*A(55)+0.22*A(58)-A(59)-A(60)-A(61)+1.1*A(71)+0.9*A(75)+0.1*A(77)
  Vdot(21) = -A(7)-A(8)+A(9)-A(42)-A(43)-A(44)+A(48)+2*A(51)+0.79*A(52)+A(53)+0.2*A(58)+A(59)+0.74*A(60)+A(61)+A(62)&
               &+1.56*A(63)+A(64)+A(72)+0.7*A(73)+A(76)+A(77)
  Vdot(22) = -A(75)-A(76)-A(77)-A(78)
  Vdot(23) = -A(58)-A(59)-A(60)-A(61)+0.55*A(75)
  Vdot(24) = -A(9)-A(45)-A(46)-A(47)+0.11*A(54)+1.1*A(55)+0.63*A(58)+A(59)+0.5*A(60)+A(61)+0.22*A(63)+0.03*A(73)+0.8&
               &*A(75)+0.2*A(76)+0.4*A(77)
  Vdot(25) = -A(2)-A(3)+A(12)-A(13)-A(17)-A(20)-A(21)-A(60)-A(64)-A(73)-A(77)
  Vdot(26) = -A(1)+0.89*A(4)+A(13)-A(14)-A(15)+A(16)-A(17)+2*A(22)-A(24)+A(26)+2*A(27)-A(28)+A(30)+A(31)-A(32)+A(34)&
               &-A(35)+A(36)+A(37)+A(48)-A(49)+A(50)-A(57)+A(61)+0.9*A(66)-A(70)+A(79)
  Vdot(27) = A(5)+2*A(6)+2*A(19)-A(20)+A(21)-A(29)-A(30)-A(32)-A(33)+A(34)-A(37)-A(40)-A(41)-A(42)+A(43)+A(45)-A(46)&
               &+0.79*A(52)-A(53)-A(54)+0.2*A(58)-A(59)+0.1*A(60)+0.3*A(62)-A(63)-A(65)-A(68)-A(71)-A(72)+0.08*A(73)-A(74)&
               &-A(76)+0.1*A(77)
  Vdot(28) = 2*A(7)+2*A(9)+A(10)+A(11)+A(20)-A(21)-A(34)-A(35)+A(36)-2*A(38)-2*A(39)+A(40)+A(41)+A(42)+A(43)+A(44)+A(48)&
               &+2*A(51)-0.21*A(52)+A(53)+0.11*A(54)+0.94*A(55)+A(56)+0.38*A(58)+A(59)+0.44*A(60)+1.7*A(62)+A(63)+0.12*A(64)&
               &+0.44*A(65)+0.9*A(66)+A(67)+0.6*A(68)+0.7*A(71)+2*A(72)+0.76*A(73)+0.6*A(75)+0.67*A(76)+0.44*A(77)
  Vdot(29) = A(1)+A(2)+0.89*A(4)-A(12)-A(14)-A(15)-A(16)+A(18)-A(43)-A(45)-A(58)-A(62)-A(75)
  Vdot(30) = -A(4)+A(15)+A(17)-A(22)-A(23)-A(24)+A(26)+A(33)-A(44)-A(47)-A(61)-A(69)-A(78)
  Vdot(31) = A(1)+0.11*A(4)+A(5)-A(13)+A(14)-A(16)-A(22)+A(23)-2*A(27)-A(28)-A(29)+A(31)-A(34)-A(48)-A(66)-A(79)-A(81)
  Vdot(32) = A(10)+A(11)+A(45)+A(46)+A(47)-A(48)-A(49)+A(50)-2*A(51)-A(52)+A(72)+0.62*A(73)+A(74)+0.2*A(76)
      
END SUBROUTINE cbm4_Fun
















SUBROUTINE cbm4_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(26)
  IRR(2) = RCT(2)*V(25)
  IRR(3) = RCT(3)*V(25)
  IRR(4) = RCT(4)*V(30)
  IRR(5) = RCT(5)*V(9)
  IRR(6) = RCT(6)*V(2)
  IRR(7) = RCT(7)*V(21)
  IRR(8) = RCT(8)*V(21)
  IRR(9) = RCT(9)*V(24)
  IRR(10) = RCT(10)*V(19)
  IRR(11) = RCT(11)*V(15)
  IRR(12) = RCT(12)*V(29)
  IRR(13) = RCT(13)*V(25)*V(31)
  IRR(14) = RCT(14)*V(26)*V(29)
  IRR(15) = RCT(15)*V(26)*V(29)
  IRR(16) = RCT(16)*V(29)*V(31)
  IRR(17) = RCT(17)*V(25)*V(26)
  IRR(18) = RCT(18)*V(1)
  IRR(19) = RCT(19)*V(1)*F(1)
  IRR(20) = RCT(20)*V(25)*V(27)
  IRR(21) = RCT(21)*V(25)*V(28)
  IRR(22) = RCT(22)*V(30)*V(31)
  IRR(23) = RCT(23)*V(26)*V(30)
  IRR(24) = RCT(24)*V(26)*V(30)
  IRR(25) = RCT(25)*V(6)*F(1)
  IRR(26) = RCT(26)*V(6)
  IRR(27) = RCT(27)*V(31)*V(31)
  IRR(28) = RCT(28)*V(26)*V(31)*F(1)
  IRR(29) = RCT(29)*V(27)*V(31)
  IRR(30) = RCT(30)*V(9)*V(27)
  IRR(31) = RCT(31)*V(9)*V(9)
  IRR(32) = RCT(32)*V(26)*V(27)
  IRR(33) = RCT(33)*V(12)*V(27)
  IRR(34) = RCT(34)*V(28)*V(31)
  IRR(35) = RCT(35)*V(26)*V(28)
  IRR(36) = RCT(36)*V(10)
  IRR(37) = RCT(37)*V(10)*V(27)
  IRR(38) = RCT(38)*V(28)*V(28)
  IRR(39) = RCT(39)*V(28)*V(28)*F(1)
  IRR(40) = RCT(40)*V(2)*V(27)
  IRR(41) = RCT(41)*V(16)*V(27)
  IRR(42) = RCT(42)*V(21)*V(27)
  IRR(43) = RCT(43)*V(21)*V(29)
  IRR(44) = RCT(44)*V(21)*V(30)
  IRR(45) = RCT(45)*V(24)*V(29)
  IRR(46) = RCT(46)*V(24)*V(27)
  IRR(47) = RCT(47)*V(24)*V(30)
  IRR(48) = RCT(48)*V(31)*V(32)
  IRR(49) = RCT(49)*V(26)*V(32)
  IRR(50) = RCT(50)*V(3)
  IRR(51) = RCT(51)*V(32)*V(32)
  IRR(52) = RCT(52)*V(28)*V(32)
  IRR(53) = RCT(53)*V(27)
  IRR(54) = RCT(54)*V(20)*V(27)
  IRR(55) = RCT(55)*V(13)
  IRR(56) = RCT(56)*V(13)
  IRR(57) = RCT(57)*V(13)*V(26)
  IRR(58) = RCT(58)*V(23)*V(29)
  IRR(59) = RCT(59)*V(23)*V(27)
  IRR(60) = RCT(60)*V(23)*V(25)
  IRR(61) = RCT(61)*V(23)*V(30)
  IRR(62) = RCT(62)*V(17)*V(29)
  IRR(63) = RCT(63)*V(17)*V(27)
  IRR(64) = RCT(64)*V(17)*V(25)
  IRR(65) = RCT(65)*V(5)*V(27)
  IRR(66) = RCT(66)*V(11)*V(31)
  IRR(67) = 4.2*V(11)
  IRR(68) = RCT(68)*V(14)*V(27)
  IRR(69) = RCT(69)*V(14)*V(30)
  IRR(70) = RCT(70)*V(4)*V(26)
  IRR(71) = RCT(71)*V(7)*V(27)
  IRR(72) = RCT(72)*V(19)*V(27)
  IRR(73) = RCT(73)*V(19)*V(25)
  IRR(74) = RCT(74)*V(15)*V(27)
  IRR(75) = RCT(75)*V(22)*V(29)
  IRR(76) = RCT(76)*V(22)*V(27)
  IRR(77) = RCT(77)*V(22)*V(25)
  IRR(78) = RCT(78)*V(22)*V(30)
  IRR(79) = RCT(79)*V(18)*V(31)
  IRR(80) = RCT(80)*V(18)*V(18)
  IRR(81) = RCT(81)*V(8)*V(31)
      
END SUBROUTINE cbm4_IRRFun
















SUBROUTINE cbm4_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(138)


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

  B(13) = RCT(13)*V(31)

  B(14) = RCT(13)*V(25)

  B(15) = RCT(14)*V(29)

  B(16) = RCT(14)*V(26)

  B(17) = RCT(15)*V(29)

  B(18) = RCT(15)*V(26)

  B(19) = RCT(16)*V(31)

  B(20) = RCT(16)*V(29)

  B(21) = RCT(17)*V(26)

  B(22) = RCT(17)*V(25)

  B(23) = RCT(18)

  B(24) = RCT(19)*F(1)

  B(26) = RCT(20)*V(27)

  B(27) = RCT(20)*V(25)

  B(28) = RCT(21)*V(28)

  B(29) = RCT(21)*V(25)

  B(30) = RCT(22)*V(31)

  B(31) = RCT(22)*V(30)

  B(32) = RCT(23)*V(30)

  B(33) = RCT(23)*V(26)

  B(34) = RCT(24)*V(30)

  B(35) = RCT(24)*V(26)

  B(36) = RCT(25)*F(1)

  B(38) = RCT(26)

  B(39) = RCT(27)*2*V(31)

  B(40) = RCT(28)*V(31)*F(1)

  B(41) = RCT(28)*V(26)*F(1)

  B(43) = RCT(29)*V(31)

  B(44) = RCT(29)*V(27)

  B(45) = RCT(30)*V(27)

  B(46) = RCT(30)*V(9)

  B(47) = RCT(31)*2*V(9)

  B(48) = RCT(32)*V(27)

  B(49) = RCT(32)*V(26)

  B(50) = RCT(33)*V(27)

  B(51) = RCT(33)*V(12)

  B(52) = RCT(34)*V(31)

  B(53) = RCT(34)*V(28)

  B(54) = RCT(35)*V(28)

  B(55) = RCT(35)*V(26)

  B(56) = RCT(36)

  B(57) = RCT(37)*V(27)

  B(58) = RCT(37)*V(10)

  B(59) = RCT(38)*2*V(28)

  B(60) = RCT(39)*2*V(28)*F(1)

  B(62) = RCT(40)*V(27)

  B(63) = RCT(40)*V(2)

  B(64) = RCT(41)*V(27)

  B(65) = RCT(41)*V(16)

  B(66) = RCT(42)*V(27)

  B(67) = RCT(42)*V(21)

  B(68) = RCT(43)*V(29)

  B(69) = RCT(43)*V(21)

  B(70) = RCT(44)*V(30)

  B(71) = RCT(44)*V(21)

  B(72) = RCT(45)*V(29)

  B(73) = RCT(45)*V(24)

  B(74) = RCT(46)*V(27)

  B(75) = RCT(46)*V(24)

  B(76) = RCT(47)*V(30)

  B(77) = RCT(47)*V(24)

  B(78) = RCT(48)*V(32)

  B(79) = RCT(48)*V(31)

  B(80) = RCT(49)*V(32)

  B(81) = RCT(49)*V(26)

  B(82) = RCT(50)

  B(83) = RCT(51)*2*V(32)

  B(84) = RCT(52)*V(32)

  B(85) = RCT(52)*V(28)

  B(86) = RCT(53)

  B(87) = RCT(54)*V(27)

  B(88) = RCT(54)*V(20)

  B(89) = RCT(55)

  B(90) = RCT(56)

  B(91) = RCT(57)*V(26)

  B(92) = RCT(57)*V(13)

  B(93) = RCT(58)*V(29)

  B(94) = RCT(58)*V(23)

  B(95) = RCT(59)*V(27)

  B(96) = RCT(59)*V(23)

  B(97) = RCT(60)*V(25)

  B(98) = RCT(60)*V(23)

  B(99) = RCT(61)*V(30)

  B(100) = RCT(61)*V(23)

  B(101) = RCT(62)*V(29)

  B(102) = RCT(62)*V(17)

  B(103) = RCT(63)*V(27)

  B(104) = RCT(63)*V(17)

  B(105) = RCT(64)*V(25)

  B(106) = RCT(64)*V(17)

  B(107) = RCT(65)*V(27)

  B(108) = RCT(65)*V(5)

  B(109) = RCT(66)*V(31)

  B(110) = RCT(66)*V(11)

  B(111) = 4.2

  B(112) = RCT(68)*V(27)

  B(113) = RCT(68)*V(14)

  B(114) = RCT(69)*V(30)

  B(115) = RCT(69)*V(14)

  B(116) = RCT(70)*V(26)

  B(117) = RCT(70)*V(4)

  B(118) = RCT(71)*V(27)

  B(119) = RCT(71)*V(7)

  B(120) = RCT(72)*V(27)

  B(121) = RCT(72)*V(19)

  B(122) = RCT(73)*V(25)

  B(123) = RCT(73)*V(19)

  B(124) = RCT(74)*V(27)

  B(125) = RCT(74)*V(15)

  B(126) = RCT(75)*V(29)

  B(127) = RCT(75)*V(22)

  B(128) = RCT(76)*V(27)

  B(129) = RCT(76)*V(22)

  B(130) = RCT(77)*V(25)

  B(131) = RCT(77)*V(22)

  B(132) = RCT(78)*V(30)

  B(133) = RCT(78)*V(22)

  B(134) = RCT(79)*V(31)

  B(135) = RCT(79)*V(18)

  B(136) = RCT(80)*2*V(18)

  B(137) = RCT(81)*V(31)

  B(138) = RCT(81)*V(8)



  JVS(1) = -B(23)-B(24)

  JVS(2) = B(3)

  JVS(3) = -B(6)-B(62)

  JVS(4) = -B(63)

  JVS(5) = B(59)+B(60)

  JVS(6) = -B(82)

  JVS(7) = B(80)

  JVS(8) = B(81)

  JVS(9) = -B(116)

  JVS(10) = 0.4*B(112)+B(114)

  JVS(11) = -B(117)

  JVS(12) = 0.4*B(113)

  JVS(13) = B(115)

  JVS(14) = -B(107)

  JVS(15) = -B(108)

  JVS(16) = -B(36)-B(38)

  JVS(17) = B(34)

  JVS(18) = B(35)

  JVS(19) = -B(118)

  JVS(20) = -B(119)

  JVS(21) = -B(137)

  JVS(22) = 0.04*B(89)

  JVS(23) = 0.13*B(87)

  JVS(24) = 0.13*B(128)+B(132)

  JVS(25) = 0.02*B(93)+0.09*B(99)

  JVS(26) = 0.13*B(88)+0.13*B(129)

  JVS(27) = 0.02*B(94)

  JVS(28) = 0.09*B(100)+B(133)

  JVS(29) = -B(138)

  JVS(30) = -B(5)-B(45)-2*B(47)

  JVS(31) = 2*B(40)

  JVS(32) = B(43)-B(46)

  JVS(33) = 2*B(41)+B(44)

  JVS(34) = -B(56)-B(57)

  JVS(35) = B(54)

  JVS(36) = -B(58)

  JVS(37) = B(55)

  JVS(38) = 0.56*B(107)

  JVS(39) = 0.3*B(118)

  JVS(40) = -B(109)-B(111)

  JVS(41) = 0.56*B(108)+0.3*B(119)

  JVS(42) = -B(110)

  JVS(43) = 2*B(36)

  JVS(44) = -B(50)

  JVS(45) = B(114)

  JVS(46) = B(70)

  JVS(47) = B(76)

  JVS(48) = B(48)

  JVS(49) = B(49)-B(51)

  JVS(50) = B(71)+B(77)+B(115)

  JVS(51) = -0.98*B(89)-B(90)-B(91)

  JVS(52) = 0.76*B(87)

  JVS(53) = -B(92)

  JVS(54) = 0.76*B(88)

  JVS(55) = 0.36*B(107)

  JVS(56) = 0.2*B(118)

  JVS(57) = B(111)

  JVS(58) = -B(112)-B(114)

  JVS(59) = 0.36*B(108)-B(113)+0.2*B(119)

  JVS(60) = -B(115)

  JVS(61) = 0

  JVS(62) = 0.8*B(118)

  JVS(63) = -B(11)-B(124)

  JVS(64) = 0.2*B(122)

  JVS(65) = 0.4*B(128)+0.2*B(130)

  JVS(66) = 0.2*B(123)+0.2*B(131)

  JVS(67) = 0.8*B(119)-B(125)+0.4*B(129)

  JVS(68) = B(11)

  JVS(69) = -B(64)

  JVS(70) = B(101)+0.42*B(105)

  JVS(71) = B(10)+2*B(120)+0.69*B(122)

  JVS(72) = B(7)+B(8)+B(66)+B(68)+B(70)

  JVS(73) = 0.5*B(126)+0.06*B(130)

  JVS(74) = 0.3*B(93)+0.33*B(97)

  JVS(75) = B(9)

  JVS(76) = 0.33*B(98)+0.42*B(106)+0.69*B(123)+0.06*B(131)

  JVS(77) = -B(65)+B(67)+2*B(121)

  JVS(78) = B(69)+0.3*B(94)+B(102)+0.5*B(127)

  JVS(79) = B(71)

  JVS(80) = -B(101)-B(103)-B(105)

  JVS(81) = 0.45*B(126)+B(128)+0.55*B(130)

  JVS(82) = -B(106)+0.55*B(131)

  JVS(83) = -B(104)+B(129)

  JVS(84) = -B(102)+0.45*B(127)

  JVS(85) = 0.08*B(107)

  JVS(86) = 0.5*B(118)

  JVS(87) = 0.96*B(89)

  JVS(88) = 0.6*B(112)

  JVS(89) = B(124)

  JVS(90) = 0.7*B(101)+B(103)

  JVS(91) = -B(134)-2*B(136)

  JVS(92) = B(120)+0.03*B(122)

  JVS(93) = 0.87*B(87)

  JVS(94) = 0.5*B(126)+B(128)

  JVS(95) = 0.28*B(93)+B(95)+0.22*B(97)+0.91*B(99)

  JVS(96) = B(9)

  JVS(97) = 0.22*B(98)+0.03*B(123)

  JVS(98) = 0

  JVS(99) = B(86)+0.87*B(88)+B(96)+B(104)+0.08*B(108)+0.6*B(113)+0.5*B(119)+B(121)+B(125)+B(129)

  JVS(100) = 0.79*B(84)

  JVS(101) = 0.28*B(94)+0.7*B(102)+0.5*B(127)

  JVS(102) = 0.91*B(100)

  JVS(103) = B(78)-B(135)

  JVS(104) = B(79)+2*B(83)+0.79*B(85)

  JVS(105) = 0.9*B(109)

  JVS(106) = 0.3*B(112)

  JVS(107) = -B(10)-B(120)-B(122)

  JVS(108) = -B(123)

  JVS(109) = 0.3*B(113)-B(121)

  JVS(110) = 0

  JVS(111) = 0.9*B(110)

  JVS(112) = 1.1*B(118)

  JVS(113) = -2.1*B(89)

  JVS(114) = -1.11*B(87)

  JVS(115) = 0.9*B(126)+0.1*B(130)

  JVS(116) = 0.22*B(93)-B(95)-B(97)-B(99)

  JVS(117) = -B(98)+0.1*B(131)

  JVS(118) = 0

  JVS(119) = -1.11*B(88)-B(96)+1.1*B(119)

  JVS(120) = 0.22*B(94)+0.9*B(127)

  JVS(121) = -B(100)

  JVS(122) = B(101)+1.56*B(103)+B(105)

  JVS(123) = B(120)+0.7*B(122)

  JVS(124) = -B(7)-B(8)-B(66)-B(68)-B(70)

  JVS(125) = B(128)+B(130)

  JVS(126) = 0.2*B(93)+B(95)+0.74*B(97)+B(99)

  JVS(127) = B(9)

  JVS(128) = 0.74*B(98)+B(106)+0.7*B(123)+B(131)

  JVS(129) = -B(67)+B(86)+B(96)+1.56*B(104)+B(121)+B(129)

  JVS(130) = 0.79*B(84)

  JVS(131) = -B(69)+0.2*B(94)+B(102)

  JVS(132) = -B(71)+B(100)

  JVS(133) = B(78)

  JVS(134) = B(79)+2*B(83)+0.79*B(85)

  JVS(135) = -B(126)-B(128)-B(130)-B(132)

  JVS(136) = -B(131)

  JVS(137) = -B(129)

  JVS(138) = -B(127)

  JVS(139) = -B(133)

  JVS(140) = 0.55*B(126)

  JVS(141) = -B(93)-B(95)-B(97)-B(99)

  JVS(142) = -B(98)

  JVS(143) = -B(96)

  JVS(144) = -B(94)+0.55*B(127)

  JVS(145) = -B(100)

  JVS(146) = 1.1*B(89)

  JVS(147) = 0.22*B(103)

  JVS(148) = 0.03*B(122)

  JVS(149) = 0.11*B(87)

  JVS(150) = 0.8*B(126)+0.2*B(128)+0.4*B(130)

  JVS(151) = 0.63*B(93)+B(95)+0.5*B(97)+B(99)

  JVS(152) = -B(9)-B(72)-B(74)-B(76)

  JVS(153) = 0.5*B(98)+0.03*B(123)+0.4*B(131)

  JVS(154) = 0

  JVS(155) = -B(75)+0.11*B(88)+B(96)+0.22*B(104)+0.2*B(129)

  JVS(156) = -B(73)+0.63*B(94)+0.8*B(127)

  JVS(157) = -B(77)+B(100)

  JVS(158) = 0

  JVS(159) = -B(105)

  JVS(160) = -B(122)

  JVS(161) = -B(130)

  JVS(162) = -B(97)

  JVS(163) = -B(2)-B(3)-B(13)-B(21)-B(26)-B(28)-B(98)-B(106)-B(123)-B(131)

  JVS(164) = -B(22)

  JVS(165) = -B(27)

  JVS(166) = -B(29)

  JVS(167) = B(12)

  JVS(168) = 0

  JVS(169) = -B(14)

  JVS(170) = B(82)

  JVS(171) = -B(116)

  JVS(172) = B(38)

  JVS(173) = B(45)+B(47)

  JVS(174) = B(56)+B(57)

  JVS(175) = 0.9*B(109)

  JVS(176) = -B(91)

  JVS(177) = 0

  JVS(178) = B(134)

  JVS(179) = 0

  JVS(180) = 0

  JVS(181) = 0

  JVS(182) = B(99)

  JVS(183) = 0

  JVS(184) = B(13)-B(21)

  JVS(185) = -B(1)-B(15)-B(17)-B(22)-B(34)-B(40)-B(48)-B(54)-B(80)-B(92)-B(117)

  JVS(186) = B(46)-B(49)+B(58)

  JVS(187) = B(52)-B(55)

  JVS(188) = -B(16)-B(18)+B(19)

  JVS(189) = 0.89*B(4)+2*B(30)-B(35)+B(100)

  JVS(190) = B(14)+B(20)+2*B(31)+2*B(39)-B(41)+B(53)+B(78)+0.9*B(110)+B(135)

  JVS(191) = B(79)-B(81)

  JVS(192) = 2*B(24)

  JVS(193) = 2*B(6)-B(62)

  JVS(194) = -B(107)

  JVS(195) = -B(118)

  JVS(196) = B(5)-B(45)

  JVS(197) = -B(57)

  JVS(198) = -B(50)

  JVS(199) = -B(112)

  JVS(200) = -B(124)

  JVS(201) = -B(64)

  JVS(202) = 0.3*B(101)-B(103)

  JVS(203) = -B(120)+0.08*B(122)

  JVS(204) = -B(87)

  JVS(205) = -B(66)+B(68)

  JVS(206) = -B(128)+0.1*B(130)

  JVS(207) = 0.2*B(93)-B(95)+0.1*B(97)

  JVS(208) = B(72)-B(74)

  JVS(209) = -B(26)+B(28)+0.1*B(98)+0.08*B(123)+0.1*B(131)

  JVS(210) = -B(48)

  JVS(211) = -B(27)-B(43)-B(46)-B(49)-B(51)-B(58)-B(63)-B(65)-B(67)-B(75)-B(86)-B(88)-B(96)-B(104)-B(108)-B(113)-B(119)&
               &-B(121)-B(125)-B(129)

  JVS(212) = B(29)+B(52)+0.79*B(84)

  JVS(213) = B(69)+B(73)+0.2*B(94)+0.3*B(102)

  JVS(214) = 0

  JVS(215) = -B(44)+B(53)

  JVS(216) = 0.79*B(85)

  JVS(217) = B(62)

  JVS(218) = 0.44*B(107)

  JVS(219) = 0.7*B(118)

  JVS(220) = B(56)

  JVS(221) = 0.9*B(109)+B(111)

  JVS(222) = 0.94*B(89)+B(90)

  JVS(223) = 0.6*B(112)

  JVS(224) = B(11)

  JVS(225) = B(64)

  JVS(226) = 1.7*B(101)+B(103)+0.12*B(105)

  JVS(227) = B(10)+2*B(120)+0.76*B(122)

  JVS(228) = 0.11*B(87)

  JVS(229) = 2*B(7)+B(66)+B(68)+B(70)

  JVS(230) = 0.6*B(126)+0.67*B(128)+0.44*B(130)

  JVS(231) = 0.38*B(93)+B(95)+0.44*B(97)

  JVS(232) = 2*B(9)

  JVS(233) = B(26)-B(28)+0.44*B(98)+0.12*B(106)+0.76*B(123)+0.44*B(131)

  JVS(234) = -B(54)

  JVS(235) = B(27)+B(63)+B(65)+B(67)+B(86)+0.11*B(88)+B(96)+B(104)+0.44*B(108)+0.6*B(113)+0.7*B(119)+2*B(121)+0.67&
               &*B(129)

  JVS(236) = -B(29)-B(52)-B(55)-2*B(59)-2*B(60)-0.21*B(84)

  JVS(237) = B(69)+0.38*B(94)+1.7*B(102)+0.6*B(127)

  JVS(238) = B(71)

  JVS(239) = -B(53)+B(78)+0.9*B(110)

  JVS(240) = B(79)+2*B(83)-0.21*B(85)

  JVS(241) = B(23)

  JVS(242) = -B(101)

  JVS(243) = -B(68)

  JVS(244) = -B(126)

  JVS(245) = -B(93)

  JVS(246) = -B(72)

  JVS(247) = B(2)

  JVS(248) = B(1)-B(15)-B(17)

  JVS(249) = 0

  JVS(250) = 0

  JVS(251) = -B(12)-B(16)-B(18)-B(19)-B(69)-B(73)-B(94)-B(102)-B(127)

  JVS(252) = 0.89*B(4)

  JVS(253) = -B(20)

  JVS(254) = 0

  JVS(255) = B(38)

  JVS(256) = B(50)

  JVS(257) = -B(114)

  JVS(258) = -B(70)

  JVS(259) = -B(132)

  JVS(260) = -B(99)

  JVS(261) = -B(76)

  JVS(262) = B(21)

  JVS(263) = B(17)+B(22)-B(32)-B(34)

  JVS(264) = B(51)

  JVS(265) = 0

  JVS(266) = B(18)

  JVS(267) = -B(4)-B(30)-B(33)-B(35)-B(71)-B(77)-B(100)-B(115)-B(133)

  JVS(268) = -B(31)

  JVS(269) = 0

  JVS(270) = -B(137)

  JVS(271) = B(5)+B(47)

  JVS(272) = -B(109)

  JVS(273) = 0

  JVS(274) = -B(134)

  JVS(275) = 0

  JVS(276) = 0

  JVS(277) = 0

  JVS(278) = 0

  JVS(279) = 0

  JVS(280) = -B(13)

  JVS(281) = B(1)+B(15)+B(32)-B(40)

  JVS(282) = -B(43)

  JVS(283) = -B(52)

  JVS(284) = B(16)-B(19)

  JVS(285) = 0.11*B(4)-B(30)+B(33)

  JVS(286) = -B(14)-B(20)-B(31)-2*B(39)-B(41)-B(44)-B(53)-B(78)-B(110)-B(135)-B(138)

  JVS(287) = -B(79)

  JVS(288) = B(82)

  JVS(289) = B(11)+B(124)

  JVS(290) = B(10)+B(120)+0.62*B(122)

  JVS(291) = 0.2*B(128)

  JVS(292) = B(72)+B(74)+B(76)

  JVS(293) = 0.62*B(123)

  JVS(294) = -B(80)

  JVS(295) = B(75)+B(121)+B(125)+0.2*B(129)

  JVS(296) = -B(84)

  JVS(297) = B(73)

  JVS(298) = B(77)

  JVS(299) = -B(78)

  JVS(300) = -B(79)-B(81)-2*B(83)-B(85)
      
END SUBROUTINE cbm4_Jac_SP














SUBROUTINE cbm4_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(300), W(32), a
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
      
END SUBROUTINE cbm4_KppDecomp



SUBROUTINE cbm4_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(300), W(32), a
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
      
END SUBROUTINE cbm4_KppDecompCmplx


SUBROUTINE cbm4_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(300), X(32), sum

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
      
END SUBROUTINE cbm4_KppSolveIndirect


SUBROUTINE cbm4_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(300), X(32), sum

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
      
END SUBROUTINE cbm4_KppSolveCmplx













SUBROUTINE cbm4_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(11) = X(11)-JVS(38)*X(5)-JVS(39)*X(7)
  X(12) = X(12)-JVS(43)*X(6)
  X(14) = X(14)-JVS(55)*X(5)-JVS(56)*X(7)-JVS(57)*X(11)
  X(15) = X(15)-JVS(62)*X(7)
  X(16) = X(16)-JVS(68)*X(15)
  X(18) = X(18)-JVS(85)*X(5)-JVS(86)*X(7)-JVS(87)*X(13)-JVS(88)*X(14)-JVS(89)*X(15)-JVS(90)*X(17)
  X(19) = X(19)-JVS(105)*X(11)-JVS(106)*X(14)
  X(20) = X(20)-JVS(112)*X(7)-JVS(113)*X(13)
  X(21) = X(21)-JVS(122)*X(17)-JVS(123)*X(19)
  X(23) = X(23)-JVS(140)*X(22)
  X(24) = X(24)-JVS(146)*X(13)-JVS(147)*X(17)-JVS(148)*X(19)-JVS(149)*X(20)-JVS(150)*X(22)-JVS(151)*X(23)
  X(25) = X(25)-JVS(159)*X(17)-JVS(160)*X(19)-JVS(161)*X(22)-JVS(162)*X(23)
  X(26) = X(26)-JVS(170)*X(3)-JVS(171)*X(4)-JVS(172)*X(6)-JVS(173)*X(9)-JVS(174)*X(10)-JVS(175)*X(11)-JVS(176)*X(13)&
            &-JVS(177)*X(14)-JVS(178)*X(18)-JVS(179)*X(19)-JVS(180)*X(20)-JVS(181)*X(22)-JVS(182)*X(23)-JVS(183)*X(24)&
            &-JVS(184)*X(25)
  X(27) = X(27)-JVS(192)*X(1)-JVS(193)*X(2)-JVS(194)*X(5)-JVS(195)*X(7)-JVS(196)*X(9)-JVS(197)*X(10)-JVS(198)*X(12)&
            &-JVS(199)*X(14)-JVS(200)*X(15)-JVS(201)*X(16)-JVS(202)*X(17)-JVS(203)*X(19)-JVS(204)*X(20)-JVS(205)*X(21)&
            &-JVS(206)*X(22)-JVS(207)*X(23)-JVS(208)*X(24)-JVS(209)*X(25)-JVS(210)*X(26)
  X(28) = X(28)-JVS(217)*X(2)-JVS(218)*X(5)-JVS(219)*X(7)-JVS(220)*X(10)-JVS(221)*X(11)-JVS(222)*X(13)-JVS(223)*X(14)&
            &-JVS(224)*X(15)-JVS(225)*X(16)-JVS(226)*X(17)-JVS(227)*X(19)-JVS(228)*X(20)-JVS(229)*X(21)-JVS(230)*X(22)&
            &-JVS(231)*X(23)-JVS(232)*X(24)-JVS(233)*X(25)-JVS(234)*X(26)-JVS(235)*X(27)
  X(29) = X(29)-JVS(241)*X(1)-JVS(242)*X(17)-JVS(243)*X(21)-JVS(244)*X(22)-JVS(245)*X(23)-JVS(246)*X(24)-JVS(247)*X(25)&
            &-JVS(248)*X(26)-JVS(249)*X(27)-JVS(250)*X(28)
  X(30) = X(30)-JVS(255)*X(6)-JVS(256)*X(12)-JVS(257)*X(14)-JVS(258)*X(21)-JVS(259)*X(22)-JVS(260)*X(23)-JVS(261)*X(24)&
            &-JVS(262)*X(25)-JVS(263)*X(26)-JVS(264)*X(27)-JVS(265)*X(28)-JVS(266)*X(29)
  X(31) = X(31)-JVS(270)*X(8)-JVS(271)*X(9)-JVS(272)*X(11)-JVS(273)*X(13)-JVS(274)*X(18)-JVS(275)*X(19)-JVS(276)*X(20)&
            &-JVS(277)*X(22)-JVS(278)*X(23)-JVS(279)*X(24)-JVS(280)*X(25)-JVS(281)*X(26)-JVS(282)*X(27)-JVS(283)*X(28)&
            &-JVS(284)*X(29)-JVS(285)*X(30)
  X(32) = X(32)-JVS(288)*X(3)-JVS(289)*X(15)-JVS(290)*X(19)-JVS(291)*X(22)-JVS(292)*X(24)-JVS(293)*X(25)-JVS(294)*X(26)&
            &-JVS(295)*X(27)-JVS(296)*X(28)-JVS(297)*X(29)-JVS(298)*X(30)-JVS(299)*X(31)
  X(32) = X(32)/JVS(300)
  X(31) = (X(31)-JVS(287)*X(32))/(JVS(286))
  X(30) = (X(30)-JVS(268)*X(31)-JVS(269)*X(32))/(JVS(267))
  X(29) = (X(29)-JVS(252)*X(30)-JVS(253)*X(31)-JVS(254)*X(32))/(JVS(251))
  X(28) = (X(28)-JVS(237)*X(29)-JVS(238)*X(30)-JVS(239)*X(31)-JVS(240)*X(32))/(JVS(236))
  X(27) = (X(27)-JVS(212)*X(28)-JVS(213)*X(29)-JVS(214)*X(30)-JVS(215)*X(31)-JVS(216)*X(32))/(JVS(211))
  X(26) = (X(26)-JVS(186)*X(27)-JVS(187)*X(28)-JVS(188)*X(29)-JVS(189)*X(30)-JVS(190)*X(31)-JVS(191)*X(32))/(JVS(185))
  X(25) = (X(25)-JVS(164)*X(26)-JVS(165)*X(27)-JVS(166)*X(28)-JVS(167)*X(29)-JVS(168)*X(30)-JVS(169)*X(31))/(JVS(163))
  X(24) = (X(24)-JVS(153)*X(25)-JVS(154)*X(26)-JVS(155)*X(27)-JVS(156)*X(29)-JVS(157)*X(30)-JVS(158)*X(31))/(JVS(152))
  X(23) = (X(23)-JVS(142)*X(25)-JVS(143)*X(27)-JVS(144)*X(29)-JVS(145)*X(30))/(JVS(141))
  X(22) = (X(22)-JVS(136)*X(25)-JVS(137)*X(27)-JVS(138)*X(29)-JVS(139)*X(30))/(JVS(135))
  X(21) = (X(21)-JVS(125)*X(22)-JVS(126)*X(23)-JVS(127)*X(24)-JVS(128)*X(25)-JVS(129)*X(27)-JVS(130)*X(28)-JVS(131)&
            &*X(29)-JVS(132)*X(30)-JVS(133)*X(31)-JVS(134)*X(32))/(JVS(124))
  X(20) = (X(20)-JVS(115)*X(22)-JVS(116)*X(23)-JVS(117)*X(25)-JVS(118)*X(26)-JVS(119)*X(27)-JVS(120)*X(29)-JVS(121)&
            &*X(30))/(JVS(114))
  X(19) = (X(19)-JVS(108)*X(25)-JVS(109)*X(27)-JVS(110)*X(30)-JVS(111)*X(31))/(JVS(107))
  X(18) = (X(18)-JVS(92)*X(19)-JVS(93)*X(20)-JVS(94)*X(22)-JVS(95)*X(23)-JVS(96)*X(24)-JVS(97)*X(25)-JVS(98)*X(26)&
            &-JVS(99)*X(27)-JVS(100)*X(28)-JVS(101)*X(29)-JVS(102)*X(30)-JVS(103)*X(31)-JVS(104)*X(32))/(JVS(91))
  X(17) = (X(17)-JVS(81)*X(22)-JVS(82)*X(25)-JVS(83)*X(27)-JVS(84)*X(29))/(JVS(80))
  X(16) = (X(16)-JVS(70)*X(17)-JVS(71)*X(19)-JVS(72)*X(21)-JVS(73)*X(22)-JVS(74)*X(23)-JVS(75)*X(24)-JVS(76)*X(25)&
            &-JVS(77)*X(27)-JVS(78)*X(29)-JVS(79)*X(30))/(JVS(69))
  X(15) = (X(15)-JVS(64)*X(19)-JVS(65)*X(22)-JVS(66)*X(25)-JVS(67)*X(27))/(JVS(63))
  X(14) = (X(14)-JVS(59)*X(27)-JVS(60)*X(30)-JVS(61)*X(31))/(JVS(58))
  X(13) = (X(13)-JVS(52)*X(20)-JVS(53)*X(26)-JVS(54)*X(27))/(JVS(51))
  X(12) = (X(12)-JVS(45)*X(14)-JVS(46)*X(21)-JVS(47)*X(24)-JVS(48)*X(26)-JVS(49)*X(27)-JVS(50)*X(30))/(JVS(44))
  X(11) = (X(11)-JVS(41)*X(27)-JVS(42)*X(31))/(JVS(40))
  X(10) = (X(10)-JVS(35)*X(26)-JVS(36)*X(27)-JVS(37)*X(28))/(JVS(34))
  X(9) = (X(9)-JVS(31)*X(26)-JVS(32)*X(27)-JVS(33)*X(31))/(JVS(30))
  X(8) = (X(8)-JVS(22)*X(13)-JVS(23)*X(20)-JVS(24)*X(22)-JVS(25)*X(23)-JVS(26)*X(27)-JVS(27)*X(29)-JVS(28)*X(30)-JVS(29)&
           &*X(31))/(JVS(21))
  X(7) = (X(7)-JVS(20)*X(27))/(JVS(19))
  X(6) = (X(6)-JVS(17)*X(26)-JVS(18)*X(30))/(JVS(16))
  X(5) = (X(5)-JVS(15)*X(27))/(JVS(14))
  X(4) = (X(4)-JVS(10)*X(14)-JVS(11)*X(26)-JVS(12)*X(27)-JVS(13)*X(30))/(JVS(9))
  X(3) = (X(3)-JVS(7)*X(26)-JVS(8)*X(32))/(JVS(6))
  X(2) = (X(2)-JVS(4)*X(27)-JVS(5)*X(28))/(JVS(3))
  X(1) = (X(1)-JVS(2)*X(25))/(JVS(1))
      
END SUBROUTINE cbm4_KppSolve
























      SUBROUTINE cbm4_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE cbm4_WCOPY



      SUBROUTINE cbm4_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE cbm4_WAXPY




      SUBROUTINE cbm4_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE cbm4_WSCAL


      REAL(kind=dp) FUNCTION cbm4_WLAMCH( C )








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
          CALL cbm4_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      cbm4_WLAMCH = Eps

      END FUNCTION cbm4_WLAMCH
     
      SUBROUTINE cbm4_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE cbm4_WLAMCH_ADD




      SUBROUTINE cbm4_SET2ZERO(N,Y)




      
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

      END SUBROUTINE cbm4_SET2ZERO



      REAL(kind=dp) FUNCTION cbm4_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      cbm4_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        cbm4_WDOT = cbm4_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         cbm4_WDOT = cbm4_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          cbm4_WDOT = cbm4_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        cbm4_WDOT = cbm4_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION cbm4_WDOT                                          




   SUBROUTINE decomp_cbm4( JVS, IER )
   
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
   W( 25 ) = JVS( 2 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 25 )
  IF ( ABS(  JVS( 3 )) < TINY(a) ) THEN
         IER = 2                                       
         RETURN
  END IF
   W( 2 ) = JVS( 3 )
   W( 27 ) = JVS( 4 )
   W( 28 ) = JVS( 5 )
  JVS( 3) = W( 2 )
  JVS( 4) = W( 27 )
  JVS( 5) = W( 28 )
  IF ( ABS(  JVS( 6 )) < TINY(a) ) THEN
         IER = 3                                       
         RETURN
  END IF
   W( 3 ) = JVS( 6 )
   W( 26 ) = JVS( 7 )
   W( 32 ) = JVS( 8 )
  JVS( 6) = W( 3 )
  JVS( 7) = W( 26 )
  JVS( 8) = W( 32 )
  IF ( ABS(  JVS( 9 )) < TINY(a) ) THEN
         IER = 4                                       
         RETURN
  END IF
   W( 4 ) = JVS( 9 )
   W( 14 ) = JVS( 10 )
   W( 26 ) = JVS( 11 )
   W( 27 ) = JVS( 12 )
   W( 30 ) = JVS( 13 )
  JVS( 9) = W( 4 )
  JVS( 10) = W( 14 )
  JVS( 11) = W( 26 )
  JVS( 12) = W( 27 )
  JVS( 13) = W( 30 )
  IF ( ABS(  JVS( 14 )) < TINY(a) ) THEN
         IER = 5                                       
         RETURN
  END IF
   W( 5 ) = JVS( 14 )
   W( 27 ) = JVS( 15 )
  JVS( 14) = W( 5 )
  JVS( 15) = W( 27 )
  IF ( ABS(  JVS( 16 )) < TINY(a) ) THEN
         IER = 6                                       
         RETURN
  END IF
   W( 6 ) = JVS( 16 )
   W( 26 ) = JVS( 17 )
   W( 30 ) = JVS( 18 )
  JVS( 16) = W( 6 )
  JVS( 17) = W( 26 )
  JVS( 18) = W( 30 )
  IF ( ABS(  JVS( 19 )) < TINY(a) ) THEN
         IER = 7                                       
         RETURN
  END IF
   W( 7 ) = JVS( 19 )
   W( 27 ) = JVS( 20 )
  JVS( 19) = W( 7 )
  JVS( 20) = W( 27 )
  IF ( ABS(  JVS( 21 )) < TINY(a) ) THEN
         IER = 8                                       
         RETURN
  END IF
   W( 8 ) = JVS( 21 )
   W( 13 ) = JVS( 22 )
   W( 20 ) = JVS( 23 )
   W( 22 ) = JVS( 24 )
   W( 23 ) = JVS( 25 )
   W( 27 ) = JVS( 26 )
   W( 29 ) = JVS( 27 )
   W( 30 ) = JVS( 28 )
   W( 31 ) = JVS( 29 )
  JVS( 21) = W( 8 )
  JVS( 22) = W( 13 )
  JVS( 23) = W( 20 )
  JVS( 24) = W( 22 )
  JVS( 25) = W( 23 )
  JVS( 26) = W( 27 )
  JVS( 27) = W( 29 )
  JVS( 28) = W( 30 )
  JVS( 29) = W( 31 )
  IF ( ABS(  JVS( 30 )) < TINY(a) ) THEN
         IER = 9                                       
         RETURN
  END IF
   W( 9 ) = JVS( 30 )
   W( 26 ) = JVS( 31 )
   W( 27 ) = JVS( 32 )
   W( 31 ) = JVS( 33 )
  JVS( 30) = W( 9 )
  JVS( 31) = W( 26 )
  JVS( 32) = W( 27 )
  JVS( 33) = W( 31 )
  IF ( ABS(  JVS( 34 )) < TINY(a) ) THEN
         IER = 10                                      
         RETURN
  END IF
   W( 10 ) = JVS( 34 )
   W( 26 ) = JVS( 35 )
   W( 27 ) = JVS( 36 )
   W( 28 ) = JVS( 37 )
  JVS( 34) = W( 10 )
  JVS( 35) = W( 26 )
  JVS( 36) = W( 27 )
  JVS( 37) = W( 28 )
  IF ( ABS(  JVS( 40 )) < TINY(a) ) THEN
         IER = 11                                      
         RETURN
  END IF
   W( 5 ) = JVS( 38 )
   W( 7 ) = JVS( 39 )
   W( 11 ) = JVS( 40 )
   W( 27 ) = JVS( 41 )
   W( 31 ) = JVS( 42 )
  a = -W( 5 ) / JVS(           14  )
  W( 5 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 15 )
  a = -W( 7 ) / JVS(           19  )
  W( 7 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 20 )
  JVS( 38) = W( 5 )
  JVS( 39) = W( 7 )
  JVS( 40) = W( 11 )
  JVS( 41) = W( 27 )
  JVS( 42) = W( 31 )
  IF ( ABS(  JVS( 44 )) < TINY(a) ) THEN
         IER = 12                                      
         RETURN
  END IF
   W( 6 ) = JVS( 43 )
   W( 12 ) = JVS( 44 )
   W( 14 ) = JVS( 45 )
   W( 21 ) = JVS( 46 )
   W( 24 ) = JVS( 47 )
   W( 26 ) = JVS( 48 )
   W( 27 ) = JVS( 49 )
   W( 30 ) = JVS( 50 )
  a = -W( 6 ) / JVS(           16  )
  W( 6 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 17 )
  W( 30 ) = W( 30 ) + a*JVS( 18 )
  JVS( 43) = W( 6 )
  JVS( 44) = W( 12 )
  JVS( 45) = W( 14 )
  JVS( 46) = W( 21 )
  JVS( 47) = W( 24 )
  JVS( 48) = W( 26 )
  JVS( 49) = W( 27 )
  JVS( 50) = W( 30 )
  IF ( ABS(  JVS( 51 )) < TINY(a) ) THEN
         IER = 13                                      
         RETURN
  END IF
   W( 13 ) = JVS( 51 )
   W( 20 ) = JVS( 52 )
   W( 26 ) = JVS( 53 )
   W( 27 ) = JVS( 54 )
  JVS( 51) = W( 13 )
  JVS( 52) = W( 20 )
  JVS( 53) = W( 26 )
  JVS( 54) = W( 27 )
  IF ( ABS(  JVS( 58 )) < TINY(a) ) THEN
         IER = 14                                      
         RETURN
  END IF
   W( 5 ) = JVS( 55 )
   W( 7 ) = JVS( 56 )
   W( 11 ) = JVS( 57 )
   W( 14 ) = JVS( 58 )
   W( 27 ) = JVS( 59 )
   W( 30 ) = JVS( 60 )
   W( 31 ) = JVS( 61 )
  a = -W( 5 ) / JVS(           14  )
  W( 5 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 15 )
  a = -W( 7 ) / JVS(           19  )
  W( 7 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 20 )
  a = -W( 11 ) / JVS(           40  )
  W( 11 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 41 )
  W( 31 ) = W( 31 ) + a*JVS( 42 )
  JVS( 55) = W( 5 )
  JVS( 56) = W( 7 )
  JVS( 57) = W( 11 )
  JVS( 58) = W( 14 )
  JVS( 59) = W( 27 )
  JVS( 60) = W( 30 )
  JVS( 61) = W( 31 )
  IF ( ABS(  JVS( 63 )) < TINY(a) ) THEN
         IER = 15                                      
         RETURN
  END IF
   W( 7 ) = JVS( 62 )
   W( 15 ) = JVS( 63 )
   W( 19 ) = JVS( 64 )
   W( 22 ) = JVS( 65 )
   W( 25 ) = JVS( 66 )
   W( 27 ) = JVS( 67 )
  a = -W( 7 ) / JVS(           19  )
  W( 7 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 20 )
  JVS( 62) = W( 7 )
  JVS( 63) = W( 15 )
  JVS( 64) = W( 19 )
  JVS( 65) = W( 22 )
  JVS( 66) = W( 25 )
  JVS( 67) = W( 27 )
  IF ( ABS(  JVS( 69 )) < TINY(a) ) THEN
         IER = 16                                      
         RETURN
  END IF
   W( 15 ) = JVS( 68 )
   W( 16 ) = JVS( 69 )
   W( 17 ) = JVS( 70 )
   W( 19 ) = JVS( 71 )
   W( 21 ) = JVS( 72 )
   W( 22 ) = JVS( 73 )
   W( 23 ) = JVS( 74 )
   W( 24 ) = JVS( 75 )
   W( 25 ) = JVS( 76 )
   W( 27 ) = JVS( 77 )
   W( 29 ) = JVS( 78 )
   W( 30 ) = JVS( 79 )
  a = -W( 15 ) / JVS(           63  )
  W( 15 ) = -a
  W( 19 ) = W( 19 ) + a*JVS( 64 )
  W( 22 ) = W( 22 ) + a*JVS( 65 )
  W( 25 ) = W( 25 ) + a*JVS( 66 )
  W( 27 ) = W( 27 ) + a*JVS( 67 )
  JVS( 68) = W( 15 )
  JVS( 69) = W( 16 )
  JVS( 70) = W( 17 )
  JVS( 71) = W( 19 )
  JVS( 72) = W( 21 )
  JVS( 73) = W( 22 )
  JVS( 74) = W( 23 )
  JVS( 75) = W( 24 )
  JVS( 76) = W( 25 )
  JVS( 77) = W( 27 )
  JVS( 78) = W( 29 )
  JVS( 79) = W( 30 )
  IF ( ABS(  JVS( 80 )) < TINY(a) ) THEN
         IER = 17                                      
         RETURN
  END IF
   W( 17 ) = JVS( 80 )
   W( 22 ) = JVS( 81 )
   W( 25 ) = JVS( 82 )
   W( 27 ) = JVS( 83 )
   W( 29 ) = JVS( 84 )
  JVS( 80) = W( 17 )
  JVS( 81) = W( 22 )
  JVS( 82) = W( 25 )
  JVS( 83) = W( 27 )
  JVS( 84) = W( 29 )
  IF ( ABS(  JVS( 91 )) < TINY(a) ) THEN
         IER = 18                                      
         RETURN
  END IF
   W( 5 ) = JVS( 85 )
   W( 7 ) = JVS( 86 )
   W( 13 ) = JVS( 87 )
   W( 14 ) = JVS( 88 )
   W( 15 ) = JVS( 89 )
   W( 17 ) = JVS( 90 )
   W( 18 ) = JVS( 91 )
   W( 19 ) = JVS( 92 )
   W( 20 ) = JVS( 93 )
   W( 22 ) = JVS( 94 )
   W( 23 ) = JVS( 95 )
   W( 24 ) = JVS( 96 )
   W( 25 ) = JVS( 97 )
   W( 26 ) = JVS( 98 )
   W( 27 ) = JVS( 99 )
   W( 28 ) = JVS( 100 )
   W( 29 ) = JVS( 101 )
   W( 30 ) = JVS( 102 )
   W( 31 ) = JVS( 103 )
   W( 32 ) = JVS( 104 )
  a = -W( 5 ) / JVS(           14  )
  W( 5 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 15 )
  a = -W( 7 ) / JVS(           19  )
  W( 7 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 20 )
  a = -W( 13 ) / JVS(           51  )
  W( 13 ) = -a
  W( 20 ) = W( 20 ) + a*JVS( 52 )
  W( 26 ) = W( 26 ) + a*JVS( 53 )
  W( 27 ) = W( 27 ) + a*JVS( 54 )
  a = -W( 14 ) / JVS(           58  )
  W( 14 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 59 )
  W( 30 ) = W( 30 ) + a*JVS( 60 )
  W( 31 ) = W( 31 ) + a*JVS( 61 )
  a = -W( 15 ) / JVS(           63  )
  W( 15 ) = -a
  W( 19 ) = W( 19 ) + a*JVS( 64 )
  W( 22 ) = W( 22 ) + a*JVS( 65 )
  W( 25 ) = W( 25 ) + a*JVS( 66 )
  W( 27 ) = W( 27 ) + a*JVS( 67 )
  a = -W( 17 ) / JVS(           80  )
  W( 17 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 81 )
  W( 25 ) = W( 25 ) + a*JVS( 82 )
  W( 27 ) = W( 27 ) + a*JVS( 83 )
  W( 29 ) = W( 29 ) + a*JVS( 84 )
  JVS( 85) = W( 5 )
  JVS( 86) = W( 7 )
  JVS( 87) = W( 13 )
  JVS( 88) = W( 14 )
  JVS( 89) = W( 15 )
  JVS( 90) = W( 17 )
  JVS( 91) = W( 18 )
  JVS( 92) = W( 19 )
  JVS( 93) = W( 20 )
  JVS( 94) = W( 22 )
  JVS( 95) = W( 23 )
  JVS( 96) = W( 24 )
  JVS( 97) = W( 25 )
  JVS( 98) = W( 26 )
  JVS( 99) = W( 27 )
  JVS( 100) = W( 28 )
  JVS( 101) = W( 29 )
  JVS( 102) = W( 30 )
  JVS( 103) = W( 31 )
  JVS( 104) = W( 32 )
  IF ( ABS(  JVS( 107 )) < TINY(a) ) THEN
         IER = 19                                      
         RETURN
  END IF
   W( 11 ) = JVS( 105 )
   W( 14 ) = JVS( 106 )
   W( 19 ) = JVS( 107 )
   W( 25 ) = JVS( 108 )
   W( 27 ) = JVS( 109 )
   W( 30 ) = JVS( 110 )
   W( 31 ) = JVS( 111 )
  a = -W( 11 ) / JVS(           40  )
  W( 11 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 41 )
  W( 31 ) = W( 31 ) + a*JVS( 42 )
  a = -W( 14 ) / JVS(           58  )
  W( 14 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 59 )
  W( 30 ) = W( 30 ) + a*JVS( 60 )
  W( 31 ) = W( 31 ) + a*JVS( 61 )
  JVS( 105) = W( 11 )
  JVS( 106) = W( 14 )
  JVS( 107) = W( 19 )
  JVS( 108) = W( 25 )
  JVS( 109) = W( 27 )
  JVS( 110) = W( 30 )
  JVS( 111) = W( 31 )
  IF ( ABS(  JVS( 114 )) < TINY(a) ) THEN
         IER = 20                                      
         RETURN
  END IF
   W( 7 ) = JVS( 112 )
   W( 13 ) = JVS( 113 )
   W( 20 ) = JVS( 114 )
   W( 22 ) = JVS( 115 )
   W( 23 ) = JVS( 116 )
   W( 25 ) = JVS( 117 )
   W( 26 ) = JVS( 118 )
   W( 27 ) = JVS( 119 )
   W( 29 ) = JVS( 120 )
   W( 30 ) = JVS( 121 )
  a = -W( 7 ) / JVS(           19  )
  W( 7 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 20 )
  a = -W( 13 ) / JVS(           51  )
  W( 13 ) = -a
  W( 20 ) = W( 20 ) + a*JVS( 52 )
  W( 26 ) = W( 26 ) + a*JVS( 53 )
  W( 27 ) = W( 27 ) + a*JVS( 54 )
  JVS( 112) = W( 7 )
  JVS( 113) = W( 13 )
  JVS( 114) = W( 20 )
  JVS( 115) = W( 22 )
  JVS( 116) = W( 23 )
  JVS( 117) = W( 25 )
  JVS( 118) = W( 26 )
  JVS( 119) = W( 27 )
  JVS( 120) = W( 29 )
  JVS( 121) = W( 30 )
  IF ( ABS(  JVS( 124 )) < TINY(a) ) THEN
         IER = 21                                      
         RETURN
  END IF
   W( 17 ) = JVS( 122 )
   W( 19 ) = JVS( 123 )
   W( 21 ) = JVS( 124 )
   W( 22 ) = JVS( 125 )
   W( 23 ) = JVS( 126 )
   W( 24 ) = JVS( 127 )
   W( 25 ) = JVS( 128 )
   W( 27 ) = JVS( 129 )
   W( 28 ) = JVS( 130 )
   W( 29 ) = JVS( 131 )
   W( 30 ) = JVS( 132 )
   W( 31 ) = JVS( 133 )
   W( 32 ) = JVS( 134 )
  a = -W( 17 ) / JVS(           80  )
  W( 17 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 81 )
  W( 25 ) = W( 25 ) + a*JVS( 82 )
  W( 27 ) = W( 27 ) + a*JVS( 83 )
  W( 29 ) = W( 29 ) + a*JVS( 84 )
  a = -W( 19 ) / JVS(          107  )
  W( 19 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 108 )
  W( 27 ) = W( 27 ) + a*JVS( 109 )
  W( 30 ) = W( 30 ) + a*JVS( 110 )
  W( 31 ) = W( 31 ) + a*JVS( 111 )
  JVS( 122) = W( 17 )
  JVS( 123) = W( 19 )
  JVS( 124) = W( 21 )
  JVS( 125) = W( 22 )
  JVS( 126) = W( 23 )
  JVS( 127) = W( 24 )
  JVS( 128) = W( 25 )
  JVS( 129) = W( 27 )
  JVS( 130) = W( 28 )
  JVS( 131) = W( 29 )
  JVS( 132) = W( 30 )
  JVS( 133) = W( 31 )
  JVS( 134) = W( 32 )
  IF ( ABS(  JVS( 135 )) < TINY(a) ) THEN
         IER = 22                                      
         RETURN
  END IF
   W( 22 ) = JVS( 135 )
   W( 25 ) = JVS( 136 )
   W( 27 ) = JVS( 137 )
   W( 29 ) = JVS( 138 )
   W( 30 ) = JVS( 139 )
  JVS( 135) = W( 22 )
  JVS( 136) = W( 25 )
  JVS( 137) = W( 27 )
  JVS( 138) = W( 29 )
  JVS( 139) = W( 30 )
  IF ( ABS(  JVS( 141 )) < TINY(a) ) THEN
         IER = 23                                      
         RETURN
  END IF
   W( 22 ) = JVS( 140 )
   W( 23 ) = JVS( 141 )
   W( 25 ) = JVS( 142 )
   W( 27 ) = JVS( 143 )
   W( 29 ) = JVS( 144 )
   W( 30 ) = JVS( 145 )
  a = -W( 22 ) / JVS(          135  )
  W( 22 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 136 )
  W( 27 ) = W( 27 ) + a*JVS( 137 )
  W( 29 ) = W( 29 ) + a*JVS( 138 )
  W( 30 ) = W( 30 ) + a*JVS( 139 )
  JVS( 140) = W( 22 )
  JVS( 141) = W( 23 )
  JVS( 142) = W( 25 )
  JVS( 143) = W( 27 )
  JVS( 144) = W( 29 )
  JVS( 145) = W( 30 )
  IF ( ABS(  JVS( 152 )) < TINY(a) ) THEN
         IER = 24                                      
         RETURN
  END IF
   W( 13 ) = JVS( 146 )
   W( 17 ) = JVS( 147 )
   W( 19 ) = JVS( 148 )
   W( 20 ) = JVS( 149 )
   W( 22 ) = JVS( 150 )
   W( 23 ) = JVS( 151 )
   W( 24 ) = JVS( 152 )
   W( 25 ) = JVS( 153 )
   W( 26 ) = JVS( 154 )
   W( 27 ) = JVS( 155 )
   W( 29 ) = JVS( 156 )
   W( 30 ) = JVS( 157 )
   W( 31 ) = JVS( 158 )
  a = -W( 13 ) / JVS(           51  )
  W( 13 ) = -a
  W( 20 ) = W( 20 ) + a*JVS( 52 )
  W( 26 ) = W( 26 ) + a*JVS( 53 )
  W( 27 ) = W( 27 ) + a*JVS( 54 )
  a = -W( 17 ) / JVS(           80  )
  W( 17 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 81 )
  W( 25 ) = W( 25 ) + a*JVS( 82 )
  W( 27 ) = W( 27 ) + a*JVS( 83 )
  W( 29 ) = W( 29 ) + a*JVS( 84 )
  a = -W( 19 ) / JVS(          107  )
  W( 19 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 108 )
  W( 27 ) = W( 27 ) + a*JVS( 109 )
  W( 30 ) = W( 30 ) + a*JVS( 110 )
  W( 31 ) = W( 31 ) + a*JVS( 111 )
  a = -W( 20 ) / JVS(          114  )
  W( 20 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 115 )
  W( 23 ) = W( 23 ) + a*JVS( 116 )
  W( 25 ) = W( 25 ) + a*JVS( 117 )
  W( 26 ) = W( 26 ) + a*JVS( 118 )
  W( 27 ) = W( 27 ) + a*JVS( 119 )
  W( 29 ) = W( 29 ) + a*JVS( 120 )
  W( 30 ) = W( 30 ) + a*JVS( 121 )
  a = -W( 22 ) / JVS(          135  )
  W( 22 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 136 )
  W( 27 ) = W( 27 ) + a*JVS( 137 )
  W( 29 ) = W( 29 ) + a*JVS( 138 )
  W( 30 ) = W( 30 ) + a*JVS( 139 )
  a = -W( 23 ) / JVS(          141  )
  W( 23 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 142 )
  W( 27 ) = W( 27 ) + a*JVS( 143 )
  W( 29 ) = W( 29 ) + a*JVS( 144 )
  W( 30 ) = W( 30 ) + a*JVS( 145 )
  JVS( 146) = W( 13 )
  JVS( 147) = W( 17 )
  JVS( 148) = W( 19 )
  JVS( 149) = W( 20 )
  JVS( 150) = W( 22 )
  JVS( 151) = W( 23 )
  JVS( 152) = W( 24 )
  JVS( 153) = W( 25 )
  JVS( 154) = W( 26 )
  JVS( 155) = W( 27 )
  JVS( 156) = W( 29 )
  JVS( 157) = W( 30 )
  JVS( 158) = W( 31 )
  IF ( ABS(  JVS( 163 )) < TINY(a) ) THEN
         IER = 25                                      
         RETURN
  END IF
   W( 17 ) = JVS( 159 )
   W( 19 ) = JVS( 160 )
   W( 22 ) = JVS( 161 )
   W( 23 ) = JVS( 162 )
   W( 25 ) = JVS( 163 )
   W( 26 ) = JVS( 164 )
   W( 27 ) = JVS( 165 )
   W( 28 ) = JVS( 166 )
   W( 29 ) = JVS( 167 )
   W( 30 ) = JVS( 168 )
   W( 31 ) = JVS( 169 )
  a = -W( 17 ) / JVS(           80  )
  W( 17 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 81 )
  W( 25 ) = W( 25 ) + a*JVS( 82 )
  W( 27 ) = W( 27 ) + a*JVS( 83 )
  W( 29 ) = W( 29 ) + a*JVS( 84 )
  a = -W( 19 ) / JVS(          107  )
  W( 19 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 108 )
  W( 27 ) = W( 27 ) + a*JVS( 109 )
  W( 30 ) = W( 30 ) + a*JVS( 110 )
  W( 31 ) = W( 31 ) + a*JVS( 111 )
  a = -W( 22 ) / JVS(          135  )
  W( 22 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 136 )
  W( 27 ) = W( 27 ) + a*JVS( 137 )
  W( 29 ) = W( 29 ) + a*JVS( 138 )
  W( 30 ) = W( 30 ) + a*JVS( 139 )
  a = -W( 23 ) / JVS(          141  )
  W( 23 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 142 )
  W( 27 ) = W( 27 ) + a*JVS( 143 )
  W( 29 ) = W( 29 ) + a*JVS( 144 )
  W( 30 ) = W( 30 ) + a*JVS( 145 )
  JVS( 159) = W( 17 )
  JVS( 160) = W( 19 )
  JVS( 161) = W( 22 )
  JVS( 162) = W( 23 )
  JVS( 163) = W( 25 )
  JVS( 164) = W( 26 )
  JVS( 165) = W( 27 )
  JVS( 166) = W( 28 )
  JVS( 167) = W( 29 )
  JVS( 168) = W( 30 )
  JVS( 169) = W( 31 )
  IF ( ABS(  JVS( 185 )) < TINY(a) ) THEN
         IER = 26                                      
         RETURN
  END IF
   W( 3 ) = JVS( 170 )
   W( 4 ) = JVS( 171 )
   W( 6 ) = JVS( 172 )
   W( 9 ) = JVS( 173 )
   W( 10 ) = JVS( 174 )
   W( 11 ) = JVS( 175 )
   W( 13 ) = JVS( 176 )
   W( 14 ) = JVS( 177 )
   W( 18 ) = JVS( 178 )
   W( 19 ) = JVS( 179 )
   W( 20 ) = JVS( 180 )
   W( 22 ) = JVS( 181 )
   W( 23 ) = JVS( 182 )
   W( 24 ) = JVS( 183 )
   W( 25 ) = JVS( 184 )
   W( 26 ) = JVS( 185 )
   W( 27 ) = JVS( 186 )
   W( 28 ) = JVS( 187 )
   W( 29 ) = JVS( 188 )
   W( 30 ) = JVS( 189 )
   W( 31 ) = JVS( 190 )
   W( 32 ) = JVS( 191 )
  a = -W( 3 ) / JVS(            6  )
  W( 3 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 7 )
  W( 32 ) = W( 32 ) + a*JVS( 8 )
  a = -W( 4 ) / JVS(            9  )
  W( 4 ) = -a
  W( 14 ) = W( 14 ) + a*JVS( 10 )
  W( 26 ) = W( 26 ) + a*JVS( 11 )
  W( 27 ) = W( 27 ) + a*JVS( 12 )
  W( 30 ) = W( 30 ) + a*JVS( 13 )
  a = -W( 6 ) / JVS(           16  )
  W( 6 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 17 )
  W( 30 ) = W( 30 ) + a*JVS( 18 )
  a = -W( 9 ) / JVS(           30  )
  W( 9 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 31 )
  W( 27 ) = W( 27 ) + a*JVS( 32 )
  W( 31 ) = W( 31 ) + a*JVS( 33 )
  a = -W( 10 ) / JVS(           34  )
  W( 10 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 35 )
  W( 27 ) = W( 27 ) + a*JVS( 36 )
  W( 28 ) = W( 28 ) + a*JVS( 37 )
  a = -W( 11 ) / JVS(           40  )
  W( 11 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 41 )
  W( 31 ) = W( 31 ) + a*JVS( 42 )
  a = -W( 13 ) / JVS(           51  )
  W( 13 ) = -a
  W( 20 ) = W( 20 ) + a*JVS( 52 )
  W( 26 ) = W( 26 ) + a*JVS( 53 )
  W( 27 ) = W( 27 ) + a*JVS( 54 )
  a = -W( 14 ) / JVS(           58  )
  W( 14 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 59 )
  W( 30 ) = W( 30 ) + a*JVS( 60 )
  W( 31 ) = W( 31 ) + a*JVS( 61 )
  a = -W( 18 ) / JVS(           91  )
  W( 18 ) = -a
  W( 19 ) = W( 19 ) + a*JVS( 92 )
  W( 20 ) = W( 20 ) + a*JVS( 93 )
  W( 22 ) = W( 22 ) + a*JVS( 94 )
  W( 23 ) = W( 23 ) + a*JVS( 95 )
  W( 24 ) = W( 24 ) + a*JVS( 96 )
  W( 25 ) = W( 25 ) + a*JVS( 97 )
  W( 26 ) = W( 26 ) + a*JVS( 98 )
  W( 27 ) = W( 27 ) + a*JVS( 99 )
  W( 28 ) = W( 28 ) + a*JVS( 100 )
  W( 29 ) = W( 29 ) + a*JVS( 101 )
  W( 30 ) = W( 30 ) + a*JVS( 102 )
  W( 31 ) = W( 31 ) + a*JVS( 103 )
  W( 32 ) = W( 32 ) + a*JVS( 104 )
  a = -W( 19 ) / JVS(          107  )
  W( 19 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 108 )
  W( 27 ) = W( 27 ) + a*JVS( 109 )
  W( 30 ) = W( 30 ) + a*JVS( 110 )
  W( 31 ) = W( 31 ) + a*JVS( 111 )
  a = -W( 20 ) / JVS(          114  )
  W( 20 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 115 )
  W( 23 ) = W( 23 ) + a*JVS( 116 )
  W( 25 ) = W( 25 ) + a*JVS( 117 )
  W( 26 ) = W( 26 ) + a*JVS( 118 )
  W( 27 ) = W( 27 ) + a*JVS( 119 )
  W( 29 ) = W( 29 ) + a*JVS( 120 )
  W( 30 ) = W( 30 ) + a*JVS( 121 )
  a = -W( 22 ) / JVS(          135  )
  W( 22 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 136 )
  W( 27 ) = W( 27 ) + a*JVS( 137 )
  W( 29 ) = W( 29 ) + a*JVS( 138 )
  W( 30 ) = W( 30 ) + a*JVS( 139 )
  a = -W( 23 ) / JVS(          141  )
  W( 23 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 142 )
  W( 27 ) = W( 27 ) + a*JVS( 143 )
  W( 29 ) = W( 29 ) + a*JVS( 144 )
  W( 30 ) = W( 30 ) + a*JVS( 145 )
  a = -W( 24 ) / JVS(          152  )
  W( 24 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 153 )
  W( 26 ) = W( 26 ) + a*JVS( 154 )
  W( 27 ) = W( 27 ) + a*JVS( 155 )
  W( 29 ) = W( 29 ) + a*JVS( 156 )
  W( 30 ) = W( 30 ) + a*JVS( 157 )
  W( 31 ) = W( 31 ) + a*JVS( 158 )
  a = -W( 25 ) / JVS(          163  )
  W( 25 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 164 )
  W( 27 ) = W( 27 ) + a*JVS( 165 )
  W( 28 ) = W( 28 ) + a*JVS( 166 )
  W( 29 ) = W( 29 ) + a*JVS( 167 )
  W( 30 ) = W( 30 ) + a*JVS( 168 )
  W( 31 ) = W( 31 ) + a*JVS( 169 )
  JVS( 170) = W( 3 )
  JVS( 171) = W( 4 )
  JVS( 172) = W( 6 )
  JVS( 173) = W( 9 )
  JVS( 174) = W( 10 )
  JVS( 175) = W( 11 )
  JVS( 176) = W( 13 )
  JVS( 177) = W( 14 )
  JVS( 178) = W( 18 )
  JVS( 179) = W( 19 )
  JVS( 180) = W( 20 )
  JVS( 181) = W( 22 )
  JVS( 182) = W( 23 )
  JVS( 183) = W( 24 )
  JVS( 184) = W( 25 )
  JVS( 185) = W( 26 )
  JVS( 186) = W( 27 )
  JVS( 187) = W( 28 )
  JVS( 188) = W( 29 )
  JVS( 189) = W( 30 )
  JVS( 190) = W( 31 )
  JVS( 191) = W( 32 )
  IF ( ABS(  JVS( 211 )) < TINY(a) ) THEN
         IER = 27                                      
         RETURN
  END IF
   W( 1 ) = JVS( 192 )
   W( 2 ) = JVS( 193 )
   W( 5 ) = JVS( 194 )
   W( 7 ) = JVS( 195 )
   W( 9 ) = JVS( 196 )
   W( 10 ) = JVS( 197 )
   W( 12 ) = JVS( 198 )
   W( 14 ) = JVS( 199 )
   W( 15 ) = JVS( 200 )
   W( 16 ) = JVS( 201 )
   W( 17 ) = JVS( 202 )
   W( 19 ) = JVS( 203 )
   W( 20 ) = JVS( 204 )
   W( 21 ) = JVS( 205 )
   W( 22 ) = JVS( 206 )
   W( 23 ) = JVS( 207 )
   W( 24 ) = JVS( 208 )
   W( 25 ) = JVS( 209 )
   W( 26 ) = JVS( 210 )
   W( 27 ) = JVS( 211 )
   W( 28 ) = JVS( 212 )
   W( 29 ) = JVS( 213 )
   W( 30 ) = JVS( 214 )
   W( 31 ) = JVS( 215 )
   W( 32 ) = JVS( 216 )
  a = -W( 1 ) / JVS(            1  )
  W( 1 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 2 )
  a = -W( 2 ) / JVS(            3  )
  W( 2 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 4 )
  W( 28 ) = W( 28 ) + a*JVS( 5 )
  a = -W( 5 ) / JVS(           14  )
  W( 5 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 15 )
  a = -W( 7 ) / JVS(           19  )
  W( 7 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 20 )
  a = -W( 9 ) / JVS(           30  )
  W( 9 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 31 )
  W( 27 ) = W( 27 ) + a*JVS( 32 )
  W( 31 ) = W( 31 ) + a*JVS( 33 )
  a = -W( 10 ) / JVS(           34  )
  W( 10 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 35 )
  W( 27 ) = W( 27 ) + a*JVS( 36 )
  W( 28 ) = W( 28 ) + a*JVS( 37 )
  a = -W( 12 ) / JVS(           44  )
  W( 12 ) = -a
  W( 14 ) = W( 14 ) + a*JVS( 45 )
  W( 21 ) = W( 21 ) + a*JVS( 46 )
  W( 24 ) = W( 24 ) + a*JVS( 47 )
  W( 26 ) = W( 26 ) + a*JVS( 48 )
  W( 27 ) = W( 27 ) + a*JVS( 49 )
  W( 30 ) = W( 30 ) + a*JVS( 50 )
  a = -W( 14 ) / JVS(           58  )
  W( 14 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 59 )
  W( 30 ) = W( 30 ) + a*JVS( 60 )
  W( 31 ) = W( 31 ) + a*JVS( 61 )
  a = -W( 15 ) / JVS(           63  )
  W( 15 ) = -a
  W( 19 ) = W( 19 ) + a*JVS( 64 )
  W( 22 ) = W( 22 ) + a*JVS( 65 )
  W( 25 ) = W( 25 ) + a*JVS( 66 )
  W( 27 ) = W( 27 ) + a*JVS( 67 )
  a = -W( 16 ) / JVS(           69  )
  W( 16 ) = -a
  W( 17 ) = W( 17 ) + a*JVS( 70 )
  W( 19 ) = W( 19 ) + a*JVS( 71 )
  W( 21 ) = W( 21 ) + a*JVS( 72 )
  W( 22 ) = W( 22 ) + a*JVS( 73 )
  W( 23 ) = W( 23 ) + a*JVS( 74 )
  W( 24 ) = W( 24 ) + a*JVS( 75 )
  W( 25 ) = W( 25 ) + a*JVS( 76 )
  W( 27 ) = W( 27 ) + a*JVS( 77 )
  W( 29 ) = W( 29 ) + a*JVS( 78 )
  W( 30 ) = W( 30 ) + a*JVS( 79 )
  a = -W( 17 ) / JVS(           80  )
  W( 17 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 81 )
  W( 25 ) = W( 25 ) + a*JVS( 82 )
  W( 27 ) = W( 27 ) + a*JVS( 83 )
  W( 29 ) = W( 29 ) + a*JVS( 84 )
  a = -W( 19 ) / JVS(          107  )
  W( 19 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 108 )
  W( 27 ) = W( 27 ) + a*JVS( 109 )
  W( 30 ) = W( 30 ) + a*JVS( 110 )
  W( 31 ) = W( 31 ) + a*JVS( 111 )
  a = -W( 20 ) / JVS(          114  )
  W( 20 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 115 )
  W( 23 ) = W( 23 ) + a*JVS( 116 )
  W( 25 ) = W( 25 ) + a*JVS( 117 )
  W( 26 ) = W( 26 ) + a*JVS( 118 )
  W( 27 ) = W( 27 ) + a*JVS( 119 )
  W( 29 ) = W( 29 ) + a*JVS( 120 )
  W( 30 ) = W( 30 ) + a*JVS( 121 )
  a = -W( 21 ) / JVS(          124  )
  W( 21 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 125 )
  W( 23 ) = W( 23 ) + a*JVS( 126 )
  W( 24 ) = W( 24 ) + a*JVS( 127 )
  W( 25 ) = W( 25 ) + a*JVS( 128 )
  W( 27 ) = W( 27 ) + a*JVS( 129 )
  W( 28 ) = W( 28 ) + a*JVS( 130 )
  W( 29 ) = W( 29 ) + a*JVS( 131 )
  W( 30 ) = W( 30 ) + a*JVS( 132 )
  W( 31 ) = W( 31 ) + a*JVS( 133 )
  W( 32 ) = W( 32 ) + a*JVS( 134 )
  a = -W( 22 ) / JVS(          135  )
  W( 22 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 136 )
  W( 27 ) = W( 27 ) + a*JVS( 137 )
  W( 29 ) = W( 29 ) + a*JVS( 138 )
  W( 30 ) = W( 30 ) + a*JVS( 139 )
  a = -W( 23 ) / JVS(          141  )
  W( 23 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 142 )
  W( 27 ) = W( 27 ) + a*JVS( 143 )
  W( 29 ) = W( 29 ) + a*JVS( 144 )
  W( 30 ) = W( 30 ) + a*JVS( 145 )
  a = -W( 24 ) / JVS(          152  )
  W( 24 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 153 )
  W( 26 ) = W( 26 ) + a*JVS( 154 )
  W( 27 ) = W( 27 ) + a*JVS( 155 )
  W( 29 ) = W( 29 ) + a*JVS( 156 )
  W( 30 ) = W( 30 ) + a*JVS( 157 )
  W( 31 ) = W( 31 ) + a*JVS( 158 )
  a = -W( 25 ) / JVS(          163  )
  W( 25 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 164 )
  W( 27 ) = W( 27 ) + a*JVS( 165 )
  W( 28 ) = W( 28 ) + a*JVS( 166 )
  W( 29 ) = W( 29 ) + a*JVS( 167 )
  W( 30 ) = W( 30 ) + a*JVS( 168 )
  W( 31 ) = W( 31 ) + a*JVS( 169 )
  a = -W( 26 ) / JVS(          185  )
  W( 26 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 186 )
  W( 28 ) = W( 28 ) + a*JVS( 187 )
  W( 29 ) = W( 29 ) + a*JVS( 188 )
  W( 30 ) = W( 30 ) + a*JVS( 189 )
  W( 31 ) = W( 31 ) + a*JVS( 190 )
  W( 32 ) = W( 32 ) + a*JVS( 191 )
  JVS( 192) = W( 1 )
  JVS( 193) = W( 2 )
  JVS( 194) = W( 5 )
  JVS( 195) = W( 7 )
  JVS( 196) = W( 9 )
  JVS( 197) = W( 10 )
  JVS( 198) = W( 12 )
  JVS( 199) = W( 14 )
  JVS( 200) = W( 15 )
  JVS( 201) = W( 16 )
  JVS( 202) = W( 17 )
  JVS( 203) = W( 19 )
  JVS( 204) = W( 20 )
  JVS( 205) = W( 21 )
  JVS( 206) = W( 22 )
  JVS( 207) = W( 23 )
  JVS( 208) = W( 24 )
  JVS( 209) = W( 25 )
  JVS( 210) = W( 26 )
  JVS( 211) = W( 27 )
  JVS( 212) = W( 28 )
  JVS( 213) = W( 29 )
  JVS( 214) = W( 30 )
  JVS( 215) = W( 31 )
  JVS( 216) = W( 32 )
  IF ( ABS(  JVS( 236 )) < TINY(a) ) THEN
         IER = 28                                      
         RETURN
  END IF
   W( 2 ) = JVS( 217 )
   W( 5 ) = JVS( 218 )
   W( 7 ) = JVS( 219 )
   W( 10 ) = JVS( 220 )
   W( 11 ) = JVS( 221 )
   W( 13 ) = JVS( 222 )
   W( 14 ) = JVS( 223 )
   W( 15 ) = JVS( 224 )
   W( 16 ) = JVS( 225 )
   W( 17 ) = JVS( 226 )
   W( 19 ) = JVS( 227 )
   W( 20 ) = JVS( 228 )
   W( 21 ) = JVS( 229 )
   W( 22 ) = JVS( 230 )
   W( 23 ) = JVS( 231 )
   W( 24 ) = JVS( 232 )
   W( 25 ) = JVS( 233 )
   W( 26 ) = JVS( 234 )
   W( 27 ) = JVS( 235 )
   W( 28 ) = JVS( 236 )
   W( 29 ) = JVS( 237 )
   W( 30 ) = JVS( 238 )
   W( 31 ) = JVS( 239 )
   W( 32 ) = JVS( 240 )
  a = -W( 2 ) / JVS(            3  )
  W( 2 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 4 )
  W( 28 ) = W( 28 ) + a*JVS( 5 )
  a = -W( 5 ) / JVS(           14  )
  W( 5 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 15 )
  a = -W( 7 ) / JVS(           19  )
  W( 7 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 20 )
  a = -W( 10 ) / JVS(           34  )
  W( 10 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 35 )
  W( 27 ) = W( 27 ) + a*JVS( 36 )
  W( 28 ) = W( 28 ) + a*JVS( 37 )
  a = -W( 11 ) / JVS(           40  )
  W( 11 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 41 )
  W( 31 ) = W( 31 ) + a*JVS( 42 )
  a = -W( 13 ) / JVS(           51  )
  W( 13 ) = -a
  W( 20 ) = W( 20 ) + a*JVS( 52 )
  W( 26 ) = W( 26 ) + a*JVS( 53 )
  W( 27 ) = W( 27 ) + a*JVS( 54 )
  a = -W( 14 ) / JVS(           58  )
  W( 14 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 59 )
  W( 30 ) = W( 30 ) + a*JVS( 60 )
  W( 31 ) = W( 31 ) + a*JVS( 61 )
  a = -W( 15 ) / JVS(           63  )
  W( 15 ) = -a
  W( 19 ) = W( 19 ) + a*JVS( 64 )
  W( 22 ) = W( 22 ) + a*JVS( 65 )
  W( 25 ) = W( 25 ) + a*JVS( 66 )
  W( 27 ) = W( 27 ) + a*JVS( 67 )
  a = -W( 16 ) / JVS(           69  )
  W( 16 ) = -a
  W( 17 ) = W( 17 ) + a*JVS( 70 )
  W( 19 ) = W( 19 ) + a*JVS( 71 )
  W( 21 ) = W( 21 ) + a*JVS( 72 )
  W( 22 ) = W( 22 ) + a*JVS( 73 )
  W( 23 ) = W( 23 ) + a*JVS( 74 )
  W( 24 ) = W( 24 ) + a*JVS( 75 )
  W( 25 ) = W( 25 ) + a*JVS( 76 )
  W( 27 ) = W( 27 ) + a*JVS( 77 )
  W( 29 ) = W( 29 ) + a*JVS( 78 )
  W( 30 ) = W( 30 ) + a*JVS( 79 )
  a = -W( 17 ) / JVS(           80  )
  W( 17 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 81 )
  W( 25 ) = W( 25 ) + a*JVS( 82 )
  W( 27 ) = W( 27 ) + a*JVS( 83 )
  W( 29 ) = W( 29 ) + a*JVS( 84 )
  a = -W( 19 ) / JVS(          107  )
  W( 19 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 108 )
  W( 27 ) = W( 27 ) + a*JVS( 109 )
  W( 30 ) = W( 30 ) + a*JVS( 110 )
  W( 31 ) = W( 31 ) + a*JVS( 111 )
  a = -W( 20 ) / JVS(          114  )
  W( 20 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 115 )
  W( 23 ) = W( 23 ) + a*JVS( 116 )
  W( 25 ) = W( 25 ) + a*JVS( 117 )
  W( 26 ) = W( 26 ) + a*JVS( 118 )
  W( 27 ) = W( 27 ) + a*JVS( 119 )
  W( 29 ) = W( 29 ) + a*JVS( 120 )
  W( 30 ) = W( 30 ) + a*JVS( 121 )
  a = -W( 21 ) / JVS(          124  )
  W( 21 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 125 )
  W( 23 ) = W( 23 ) + a*JVS( 126 )
  W( 24 ) = W( 24 ) + a*JVS( 127 )
  W( 25 ) = W( 25 ) + a*JVS( 128 )
  W( 27 ) = W( 27 ) + a*JVS( 129 )
  W( 28 ) = W( 28 ) + a*JVS( 130 )
  W( 29 ) = W( 29 ) + a*JVS( 131 )
  W( 30 ) = W( 30 ) + a*JVS( 132 )
  W( 31 ) = W( 31 ) + a*JVS( 133 )
  W( 32 ) = W( 32 ) + a*JVS( 134 )
  a = -W( 22 ) / JVS(          135  )
  W( 22 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 136 )
  W( 27 ) = W( 27 ) + a*JVS( 137 )
  W( 29 ) = W( 29 ) + a*JVS( 138 )
  W( 30 ) = W( 30 ) + a*JVS( 139 )
  a = -W( 23 ) / JVS(          141  )
  W( 23 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 142 )
  W( 27 ) = W( 27 ) + a*JVS( 143 )
  W( 29 ) = W( 29 ) + a*JVS( 144 )
  W( 30 ) = W( 30 ) + a*JVS( 145 )
  a = -W( 24 ) / JVS(          152  )
  W( 24 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 153 )
  W( 26 ) = W( 26 ) + a*JVS( 154 )
  W( 27 ) = W( 27 ) + a*JVS( 155 )
  W( 29 ) = W( 29 ) + a*JVS( 156 )
  W( 30 ) = W( 30 ) + a*JVS( 157 )
  W( 31 ) = W( 31 ) + a*JVS( 158 )
  a = -W( 25 ) / JVS(          163  )
  W( 25 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 164 )
  W( 27 ) = W( 27 ) + a*JVS( 165 )
  W( 28 ) = W( 28 ) + a*JVS( 166 )
  W( 29 ) = W( 29 ) + a*JVS( 167 )
  W( 30 ) = W( 30 ) + a*JVS( 168 )
  W( 31 ) = W( 31 ) + a*JVS( 169 )
  a = -W( 26 ) / JVS(          185  )
  W( 26 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 186 )
  W( 28 ) = W( 28 ) + a*JVS( 187 )
  W( 29 ) = W( 29 ) + a*JVS( 188 )
  W( 30 ) = W( 30 ) + a*JVS( 189 )
  W( 31 ) = W( 31 ) + a*JVS( 190 )
  W( 32 ) = W( 32 ) + a*JVS( 191 )
  a = -W( 27 ) / JVS(          211  )
  W( 27 ) = -a
  W( 28 ) = W( 28 ) + a*JVS( 212 )
  W( 29 ) = W( 29 ) + a*JVS( 213 )
  W( 30 ) = W( 30 ) + a*JVS( 214 )
  W( 31 ) = W( 31 ) + a*JVS( 215 )
  W( 32 ) = W( 32 ) + a*JVS( 216 )
  JVS( 217) = W( 2 )
  JVS( 218) = W( 5 )
  JVS( 219) = W( 7 )
  JVS( 220) = W( 10 )
  JVS( 221) = W( 11 )
  JVS( 222) = W( 13 )
  JVS( 223) = W( 14 )
  JVS( 224) = W( 15 )
  JVS( 225) = W( 16 )
  JVS( 226) = W( 17 )
  JVS( 227) = W( 19 )
  JVS( 228) = W( 20 )
  JVS( 229) = W( 21 )
  JVS( 230) = W( 22 )
  JVS( 231) = W( 23 )
  JVS( 232) = W( 24 )
  JVS( 233) = W( 25 )
  JVS( 234) = W( 26 )
  JVS( 235) = W( 27 )
  JVS( 236) = W( 28 )
  JVS( 237) = W( 29 )
  JVS( 238) = W( 30 )
  JVS( 239) = W( 31 )
  JVS( 240) = W( 32 )
  IF ( ABS(  JVS( 251 )) < TINY(a) ) THEN
         IER = 29                                      
         RETURN
  END IF
   W( 1 ) = JVS( 241 )
   W( 17 ) = JVS( 242 )
   W( 21 ) = JVS( 243 )
   W( 22 ) = JVS( 244 )
   W( 23 ) = JVS( 245 )
   W( 24 ) = JVS( 246 )
   W( 25 ) = JVS( 247 )
   W( 26 ) = JVS( 248 )
   W( 27 ) = JVS( 249 )
   W( 28 ) = JVS( 250 )
   W( 29 ) = JVS( 251 )
   W( 30 ) = JVS( 252 )
   W( 31 ) = JVS( 253 )
   W( 32 ) = JVS( 254 )
  a = -W( 1 ) / JVS(            1  )
  W( 1 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 2 )
  a = -W( 17 ) / JVS(           80  )
  W( 17 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 81 )
  W( 25 ) = W( 25 ) + a*JVS( 82 )
  W( 27 ) = W( 27 ) + a*JVS( 83 )
  W( 29 ) = W( 29 ) + a*JVS( 84 )
  a = -W( 21 ) / JVS(          124  )
  W( 21 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 125 )
  W( 23 ) = W( 23 ) + a*JVS( 126 )
  W( 24 ) = W( 24 ) + a*JVS( 127 )
  W( 25 ) = W( 25 ) + a*JVS( 128 )
  W( 27 ) = W( 27 ) + a*JVS( 129 )
  W( 28 ) = W( 28 ) + a*JVS( 130 )
  W( 29 ) = W( 29 ) + a*JVS( 131 )
  W( 30 ) = W( 30 ) + a*JVS( 132 )
  W( 31 ) = W( 31 ) + a*JVS( 133 )
  W( 32 ) = W( 32 ) + a*JVS( 134 )
  a = -W( 22 ) / JVS(          135  )
  W( 22 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 136 )
  W( 27 ) = W( 27 ) + a*JVS( 137 )
  W( 29 ) = W( 29 ) + a*JVS( 138 )
  W( 30 ) = W( 30 ) + a*JVS( 139 )
  a = -W( 23 ) / JVS(          141  )
  W( 23 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 142 )
  W( 27 ) = W( 27 ) + a*JVS( 143 )
  W( 29 ) = W( 29 ) + a*JVS( 144 )
  W( 30 ) = W( 30 ) + a*JVS( 145 )
  a = -W( 24 ) / JVS(          152  )
  W( 24 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 153 )
  W( 26 ) = W( 26 ) + a*JVS( 154 )
  W( 27 ) = W( 27 ) + a*JVS( 155 )
  W( 29 ) = W( 29 ) + a*JVS( 156 )
  W( 30 ) = W( 30 ) + a*JVS( 157 )
  W( 31 ) = W( 31 ) + a*JVS( 158 )
  a = -W( 25 ) / JVS(          163  )
  W( 25 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 164 )
  W( 27 ) = W( 27 ) + a*JVS( 165 )
  W( 28 ) = W( 28 ) + a*JVS( 166 )
  W( 29 ) = W( 29 ) + a*JVS( 167 )
  W( 30 ) = W( 30 ) + a*JVS( 168 )
  W( 31 ) = W( 31 ) + a*JVS( 169 )
  a = -W( 26 ) / JVS(          185  )
  W( 26 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 186 )
  W( 28 ) = W( 28 ) + a*JVS( 187 )
  W( 29 ) = W( 29 ) + a*JVS( 188 )
  W( 30 ) = W( 30 ) + a*JVS( 189 )
  W( 31 ) = W( 31 ) + a*JVS( 190 )
  W( 32 ) = W( 32 ) + a*JVS( 191 )
  a = -W( 27 ) / JVS(          211  )
  W( 27 ) = -a
  W( 28 ) = W( 28 ) + a*JVS( 212 )
  W( 29 ) = W( 29 ) + a*JVS( 213 )
  W( 30 ) = W( 30 ) + a*JVS( 214 )
  W( 31 ) = W( 31 ) + a*JVS( 215 )
  W( 32 ) = W( 32 ) + a*JVS( 216 )
  a = -W( 28 ) / JVS(          236  )
  W( 28 ) = -a
  W( 29 ) = W( 29 ) + a*JVS( 237 )
  W( 30 ) = W( 30 ) + a*JVS( 238 )
  W( 31 ) = W( 31 ) + a*JVS( 239 )
  W( 32 ) = W( 32 ) + a*JVS( 240 )
  JVS( 241) = W( 1 )
  JVS( 242) = W( 17 )
  JVS( 243) = W( 21 )
  JVS( 244) = W( 22 )
  JVS( 245) = W( 23 )
  JVS( 246) = W( 24 )
  JVS( 247) = W( 25 )
  JVS( 248) = W( 26 )
  JVS( 249) = W( 27 )
  JVS( 250) = W( 28 )
  JVS( 251) = W( 29 )
  JVS( 252) = W( 30 )
  JVS( 253) = W( 31 )
  JVS( 254) = W( 32 )
  IF ( ABS(  JVS( 267 )) < TINY(a) ) THEN
         IER = 30                                      
         RETURN
  END IF
   W( 6 ) = JVS( 255 )
   W( 12 ) = JVS( 256 )
   W( 14 ) = JVS( 257 )
   W( 21 ) = JVS( 258 )
   W( 22 ) = JVS( 259 )
   W( 23 ) = JVS( 260 )
   W( 24 ) = JVS( 261 )
   W( 25 ) = JVS( 262 )
   W( 26 ) = JVS( 263 )
   W( 27 ) = JVS( 264 )
   W( 28 ) = JVS( 265 )
   W( 29 ) = JVS( 266 )
   W( 30 ) = JVS( 267 )
   W( 31 ) = JVS( 268 )
   W( 32 ) = JVS( 269 )
  a = -W( 6 ) / JVS(           16  )
  W( 6 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 17 )
  W( 30 ) = W( 30 ) + a*JVS( 18 )
  a = -W( 12 ) / JVS(           44  )
  W( 12 ) = -a
  W( 14 ) = W( 14 ) + a*JVS( 45 )
  W( 21 ) = W( 21 ) + a*JVS( 46 )
  W( 24 ) = W( 24 ) + a*JVS( 47 )
  W( 26 ) = W( 26 ) + a*JVS( 48 )
  W( 27 ) = W( 27 ) + a*JVS( 49 )
  W( 30 ) = W( 30 ) + a*JVS( 50 )
  a = -W( 14 ) / JVS(           58  )
  W( 14 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 59 )
  W( 30 ) = W( 30 ) + a*JVS( 60 )
  W( 31 ) = W( 31 ) + a*JVS( 61 )
  a = -W( 21 ) / JVS(          124  )
  W( 21 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 125 )
  W( 23 ) = W( 23 ) + a*JVS( 126 )
  W( 24 ) = W( 24 ) + a*JVS( 127 )
  W( 25 ) = W( 25 ) + a*JVS( 128 )
  W( 27 ) = W( 27 ) + a*JVS( 129 )
  W( 28 ) = W( 28 ) + a*JVS( 130 )
  W( 29 ) = W( 29 ) + a*JVS( 131 )
  W( 30 ) = W( 30 ) + a*JVS( 132 )
  W( 31 ) = W( 31 ) + a*JVS( 133 )
  W( 32 ) = W( 32 ) + a*JVS( 134 )
  a = -W( 22 ) / JVS(          135  )
  W( 22 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 136 )
  W( 27 ) = W( 27 ) + a*JVS( 137 )
  W( 29 ) = W( 29 ) + a*JVS( 138 )
  W( 30 ) = W( 30 ) + a*JVS( 139 )
  a = -W( 23 ) / JVS(          141  )
  W( 23 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 142 )
  W( 27 ) = W( 27 ) + a*JVS( 143 )
  W( 29 ) = W( 29 ) + a*JVS( 144 )
  W( 30 ) = W( 30 ) + a*JVS( 145 )
  a = -W( 24 ) / JVS(          152  )
  W( 24 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 153 )
  W( 26 ) = W( 26 ) + a*JVS( 154 )
  W( 27 ) = W( 27 ) + a*JVS( 155 )
  W( 29 ) = W( 29 ) + a*JVS( 156 )
  W( 30 ) = W( 30 ) + a*JVS( 157 )
  W( 31 ) = W( 31 ) + a*JVS( 158 )
  a = -W( 25 ) / JVS(          163  )
  W( 25 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 164 )
  W( 27 ) = W( 27 ) + a*JVS( 165 )
  W( 28 ) = W( 28 ) + a*JVS( 166 )
  W( 29 ) = W( 29 ) + a*JVS( 167 )
  W( 30 ) = W( 30 ) + a*JVS( 168 )
  W( 31 ) = W( 31 ) + a*JVS( 169 )
  a = -W( 26 ) / JVS(          185  )
  W( 26 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 186 )
  W( 28 ) = W( 28 ) + a*JVS( 187 )
  W( 29 ) = W( 29 ) + a*JVS( 188 )
  W( 30 ) = W( 30 ) + a*JVS( 189 )
  W( 31 ) = W( 31 ) + a*JVS( 190 )
  W( 32 ) = W( 32 ) + a*JVS( 191 )
  a = -W( 27 ) / JVS(          211  )
  W( 27 ) = -a
  W( 28 ) = W( 28 ) + a*JVS( 212 )
  W( 29 ) = W( 29 ) + a*JVS( 213 )
  W( 30 ) = W( 30 ) + a*JVS( 214 )
  W( 31 ) = W( 31 ) + a*JVS( 215 )
  W( 32 ) = W( 32 ) + a*JVS( 216 )
  a = -W( 28 ) / JVS(          236  )
  W( 28 ) = -a
  W( 29 ) = W( 29 ) + a*JVS( 237 )
  W( 30 ) = W( 30 ) + a*JVS( 238 )
  W( 31 ) = W( 31 ) + a*JVS( 239 )
  W( 32 ) = W( 32 ) + a*JVS( 240 )
  a = -W( 29 ) / JVS(          251  )
  W( 29 ) = -a
  W( 30 ) = W( 30 ) + a*JVS( 252 )
  W( 31 ) = W( 31 ) + a*JVS( 253 )
  W( 32 ) = W( 32 ) + a*JVS( 254 )
  JVS( 255) = W( 6 )
  JVS( 256) = W( 12 )
  JVS( 257) = W( 14 )
  JVS( 258) = W( 21 )
  JVS( 259) = W( 22 )
  JVS( 260) = W( 23 )
  JVS( 261) = W( 24 )
  JVS( 262) = W( 25 )
  JVS( 263) = W( 26 )
  JVS( 264) = W( 27 )
  JVS( 265) = W( 28 )
  JVS( 266) = W( 29 )
  JVS( 267) = W( 30 )
  JVS( 268) = W( 31 )
  JVS( 269) = W( 32 )
  IF ( ABS(  JVS( 286 )) < TINY(a) ) THEN
         IER = 31                                      
         RETURN
  END IF
   W( 8 ) = JVS( 270 )
   W( 9 ) = JVS( 271 )
   W( 11 ) = JVS( 272 )
   W( 13 ) = JVS( 273 )
   W( 18 ) = JVS( 274 )
   W( 19 ) = JVS( 275 )
   W( 20 ) = JVS( 276 )
   W( 22 ) = JVS( 277 )
   W( 23 ) = JVS( 278 )
   W( 24 ) = JVS( 279 )
   W( 25 ) = JVS( 280 )
   W( 26 ) = JVS( 281 )
   W( 27 ) = JVS( 282 )
   W( 28 ) = JVS( 283 )
   W( 29 ) = JVS( 284 )
   W( 30 ) = JVS( 285 )
   W( 31 ) = JVS( 286 )
   W( 32 ) = JVS( 287 )
  a = -W( 8 ) / JVS(           21  )
  W( 8 ) = -a
  W( 13 ) = W( 13 ) + a*JVS( 22 )
  W( 20 ) = W( 20 ) + a*JVS( 23 )
  W( 22 ) = W( 22 ) + a*JVS( 24 )
  W( 23 ) = W( 23 ) + a*JVS( 25 )
  W( 27 ) = W( 27 ) + a*JVS( 26 )
  W( 29 ) = W( 29 ) + a*JVS( 27 )
  W( 30 ) = W( 30 ) + a*JVS( 28 )
  W( 31 ) = W( 31 ) + a*JVS( 29 )
  a = -W( 9 ) / JVS(           30  )
  W( 9 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 31 )
  W( 27 ) = W( 27 ) + a*JVS( 32 )
  W( 31 ) = W( 31 ) + a*JVS( 33 )
  a = -W( 11 ) / JVS(           40  )
  W( 11 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 41 )
  W( 31 ) = W( 31 ) + a*JVS( 42 )
  a = -W( 13 ) / JVS(           51  )
  W( 13 ) = -a
  W( 20 ) = W( 20 ) + a*JVS( 52 )
  W( 26 ) = W( 26 ) + a*JVS( 53 )
  W( 27 ) = W( 27 ) + a*JVS( 54 )
  a = -W( 18 ) / JVS(           91  )
  W( 18 ) = -a
  W( 19 ) = W( 19 ) + a*JVS( 92 )
  W( 20 ) = W( 20 ) + a*JVS( 93 )
  W( 22 ) = W( 22 ) + a*JVS( 94 )
  W( 23 ) = W( 23 ) + a*JVS( 95 )
  W( 24 ) = W( 24 ) + a*JVS( 96 )
  W( 25 ) = W( 25 ) + a*JVS( 97 )
  W( 26 ) = W( 26 ) + a*JVS( 98 )
  W( 27 ) = W( 27 ) + a*JVS( 99 )
  W( 28 ) = W( 28 ) + a*JVS( 100 )
  W( 29 ) = W( 29 ) + a*JVS( 101 )
  W( 30 ) = W( 30 ) + a*JVS( 102 )
  W( 31 ) = W( 31 ) + a*JVS( 103 )
  W( 32 ) = W( 32 ) + a*JVS( 104 )
  a = -W( 19 ) / JVS(          107  )
  W( 19 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 108 )
  W( 27 ) = W( 27 ) + a*JVS( 109 )
  W( 30 ) = W( 30 ) + a*JVS( 110 )
  W( 31 ) = W( 31 ) + a*JVS( 111 )
  a = -W( 20 ) / JVS(          114  )
  W( 20 ) = -a
  W( 22 ) = W( 22 ) + a*JVS( 115 )
  W( 23 ) = W( 23 ) + a*JVS( 116 )
  W( 25 ) = W( 25 ) + a*JVS( 117 )
  W( 26 ) = W( 26 ) + a*JVS( 118 )
  W( 27 ) = W( 27 ) + a*JVS( 119 )
  W( 29 ) = W( 29 ) + a*JVS( 120 )
  W( 30 ) = W( 30 ) + a*JVS( 121 )
  a = -W( 22 ) / JVS(          135  )
  W( 22 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 136 )
  W( 27 ) = W( 27 ) + a*JVS( 137 )
  W( 29 ) = W( 29 ) + a*JVS( 138 )
  W( 30 ) = W( 30 ) + a*JVS( 139 )
  a = -W( 23 ) / JVS(          141  )
  W( 23 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 142 )
  W( 27 ) = W( 27 ) + a*JVS( 143 )
  W( 29 ) = W( 29 ) + a*JVS( 144 )
  W( 30 ) = W( 30 ) + a*JVS( 145 )
  a = -W( 24 ) / JVS(          152  )
  W( 24 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 153 )
  W( 26 ) = W( 26 ) + a*JVS( 154 )
  W( 27 ) = W( 27 ) + a*JVS( 155 )
  W( 29 ) = W( 29 ) + a*JVS( 156 )
  W( 30 ) = W( 30 ) + a*JVS( 157 )
  W( 31 ) = W( 31 ) + a*JVS( 158 )
  a = -W( 25 ) / JVS(          163  )
  W( 25 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 164 )
  W( 27 ) = W( 27 ) + a*JVS( 165 )
  W( 28 ) = W( 28 ) + a*JVS( 166 )
  W( 29 ) = W( 29 ) + a*JVS( 167 )
  W( 30 ) = W( 30 ) + a*JVS( 168 )
  W( 31 ) = W( 31 ) + a*JVS( 169 )
  a = -W( 26 ) / JVS(          185  )
  W( 26 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 186 )
  W( 28 ) = W( 28 ) + a*JVS( 187 )
  W( 29 ) = W( 29 ) + a*JVS( 188 )
  W( 30 ) = W( 30 ) + a*JVS( 189 )
  W( 31 ) = W( 31 ) + a*JVS( 190 )
  W( 32 ) = W( 32 ) + a*JVS( 191 )
  a = -W( 27 ) / JVS(          211  )
  W( 27 ) = -a
  W( 28 ) = W( 28 ) + a*JVS( 212 )
  W( 29 ) = W( 29 ) + a*JVS( 213 )
  W( 30 ) = W( 30 ) + a*JVS( 214 )
  W( 31 ) = W( 31 ) + a*JVS( 215 )
  W( 32 ) = W( 32 ) + a*JVS( 216 )
  a = -W( 28 ) / JVS(          236  )
  W( 28 ) = -a
  W( 29 ) = W( 29 ) + a*JVS( 237 )
  W( 30 ) = W( 30 ) + a*JVS( 238 )
  W( 31 ) = W( 31 ) + a*JVS( 239 )
  W( 32 ) = W( 32 ) + a*JVS( 240 )
  a = -W( 29 ) / JVS(          251  )
  W( 29 ) = -a
  W( 30 ) = W( 30 ) + a*JVS( 252 )
  W( 31 ) = W( 31 ) + a*JVS( 253 )
  W( 32 ) = W( 32 ) + a*JVS( 254 )
  a = -W( 30 ) / JVS(          267  )
  W( 30 ) = -a
  W( 31 ) = W( 31 ) + a*JVS( 268 )
  W( 32 ) = W( 32 ) + a*JVS( 269 )
  JVS( 270) = W( 8 )
  JVS( 271) = W( 9 )
  JVS( 272) = W( 11 )
  JVS( 273) = W( 13 )
  JVS( 274) = W( 18 )
  JVS( 275) = W( 19 )
  JVS( 276) = W( 20 )
  JVS( 277) = W( 22 )
  JVS( 278) = W( 23 )
  JVS( 279) = W( 24 )
  JVS( 280) = W( 25 )
  JVS( 281) = W( 26 )
  JVS( 282) = W( 27 )
  JVS( 283) = W( 28 )
  JVS( 284) = W( 29 )
  JVS( 285) = W( 30 )
  JVS( 286) = W( 31 )
  JVS( 287) = W( 32 )
  IF ( ABS(  JVS( 300 )) < TINY(a) ) THEN
         IER = 32                                      
         RETURN
  END IF
   W( 3 ) = JVS( 288 )
   W( 15 ) = JVS( 289 )
   W( 19 ) = JVS( 290 )
   W( 22 ) = JVS( 291 )
   W( 24 ) = JVS( 292 )
   W( 25 ) = JVS( 293 )
   W( 26 ) = JVS( 294 )
   W( 27 ) = JVS( 295 )
   W( 28 ) = JVS( 296 )
   W( 29 ) = JVS( 297 )
   W( 30 ) = JVS( 298 )
   W( 31 ) = JVS( 299 )
   W( 32 ) = JVS( 300 )
  a = -W( 3 ) / JVS(            6  )
  W( 3 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 7 )
  W( 32 ) = W( 32 ) + a*JVS( 8 )
  a = -W( 15 ) / JVS(           63  )
  W( 15 ) = -a
  W( 19 ) = W( 19 ) + a*JVS( 64 )
  W( 22 ) = W( 22 ) + a*JVS( 65 )
  W( 25 ) = W( 25 ) + a*JVS( 66 )
  W( 27 ) = W( 27 ) + a*JVS( 67 )
  a = -W( 19 ) / JVS(          107  )
  W( 19 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 108 )
  W( 27 ) = W( 27 ) + a*JVS( 109 )
  W( 30 ) = W( 30 ) + a*JVS( 110 )
  W( 31 ) = W( 31 ) + a*JVS( 111 )
  a = -W( 22 ) / JVS(          135  )
  W( 22 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 136 )
  W( 27 ) = W( 27 ) + a*JVS( 137 )
  W( 29 ) = W( 29 ) + a*JVS( 138 )
  W( 30 ) = W( 30 ) + a*JVS( 139 )
  a = -W( 24 ) / JVS(          152  )
  W( 24 ) = -a
  W( 25 ) = W( 25 ) + a*JVS( 153 )
  W( 26 ) = W( 26 ) + a*JVS( 154 )
  W( 27 ) = W( 27 ) + a*JVS( 155 )
  W( 29 ) = W( 29 ) + a*JVS( 156 )
  W( 30 ) = W( 30 ) + a*JVS( 157 )
  W( 31 ) = W( 31 ) + a*JVS( 158 )
  a = -W( 25 ) / JVS(          163  )
  W( 25 ) = -a
  W( 26 ) = W( 26 ) + a*JVS( 164 )
  W( 27 ) = W( 27 ) + a*JVS( 165 )
  W( 28 ) = W( 28 ) + a*JVS( 166 )
  W( 29 ) = W( 29 ) + a*JVS( 167 )
  W( 30 ) = W( 30 ) + a*JVS( 168 )
  W( 31 ) = W( 31 ) + a*JVS( 169 )
  a = -W( 26 ) / JVS(          185  )
  W( 26 ) = -a
  W( 27 ) = W( 27 ) + a*JVS( 186 )
  W( 28 ) = W( 28 ) + a*JVS( 187 )
  W( 29 ) = W( 29 ) + a*JVS( 188 )
  W( 30 ) = W( 30 ) + a*JVS( 189 )
  W( 31 ) = W( 31 ) + a*JVS( 190 )
  W( 32 ) = W( 32 ) + a*JVS( 191 )
  a = -W( 27 ) / JVS(          211  )
  W( 27 ) = -a
  W( 28 ) = W( 28 ) + a*JVS( 212 )
  W( 29 ) = W( 29 ) + a*JVS( 213 )
  W( 30 ) = W( 30 ) + a*JVS( 214 )
  W( 31 ) = W( 31 ) + a*JVS( 215 )
  W( 32 ) = W( 32 ) + a*JVS( 216 )
  a = -W( 28 ) / JVS(          236  )
  W( 28 ) = -a
  W( 29 ) = W( 29 ) + a*JVS( 237 )
  W( 30 ) = W( 30 ) + a*JVS( 238 )
  W( 31 ) = W( 31 ) + a*JVS( 239 )
  W( 32 ) = W( 32 ) + a*JVS( 240 )
  a = -W( 29 ) / JVS(          251  )
  W( 29 ) = -a
  W( 30 ) = W( 30 ) + a*JVS( 252 )
  W( 31 ) = W( 31 ) + a*JVS( 253 )
  W( 32 ) = W( 32 ) + a*JVS( 254 )
  a = -W( 30 ) / JVS(          267  )
  W( 30 ) = -a
  W( 31 ) = W( 31 ) + a*JVS( 268 )
  W( 32 ) = W( 32 ) + a*JVS( 269 )
  a = -W( 31 ) / JVS(          286  )
  W( 31 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 287 )
  JVS( 288) = W( 3 )
  JVS( 289) = W( 15 )
  JVS( 290) = W( 19 )
  JVS( 291) = W( 22 )
  JVS( 292) = W( 24 )
  JVS( 293) = W( 25 )
  JVS( 294) = W( 26 )
  JVS( 295) = W( 27 )
  JVS( 296) = W( 28 )
  JVS( 297) = W( 29 )
  JVS( 298) = W( 30 )
  JVS( 299) = W( 31 )
  JVS( 300) = W( 32 )
   
   END SUBROUTINE decomp_cbm4
 


END MODULE cbm4_Integrator
