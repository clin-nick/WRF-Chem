
































MODULE mozart_Integrator

 USE mozart_Parameters
 USE mozart_Precision
 USE mozart_JacobianSP

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

SUBROUTINE  mozart_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE mozart_Parameters

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

   CALL mozart_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  mozart_INTEGRATE


SUBROUTINE  mozart_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE mozart_Parameters

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
      CALL mozart_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL mozart_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = mozart_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL mozart_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL mozart_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL mozart_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL mozart_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL mozart_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL mozart_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL mozart_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL mozart_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL mozart_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL mozart_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL mozart_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL mozart_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL mozart_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL mozart_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL mozart_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  mozart_ros_ErrorMsg(Code,T,H,IERR)



   USE mozart_Precision

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

 END SUBROUTINE  mozart_ros_ErrorMsg


 SUBROUTINE  mozart_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL mozart_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL mozart_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL mozart_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL mozart_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL mozart_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL mozart_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL mozart_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL mozart_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL mozart_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL mozart_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL mozart_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL mozart_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL mozart_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL mozart_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL mozart_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL mozart_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL mozart_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL mozart_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL mozart_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL mozart_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = mozart_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL mozart_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  mozart_ros_Integrator



  REAL(kind=dp) FUNCTION  mozart_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    mozart_ros_ErrorNorm = Err

  END FUNCTION  mozart_ros_ErrorNorm



  SUBROUTINE mozart_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL mozart_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL mozart_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL mozart_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  mozart_ros_FunTimeDeriv



  SUBROUTINE  mozart_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL mozart_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL mozart_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL mozart_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  mozart_ros_PrepareMatrix



  SUBROUTINE  mozart_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 



CALL decomp_mozart ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  mozart_ros_Decomp



  SUBROUTINE  mozart_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL mozart_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  mozart_ros_Solve




  SUBROUTINE  mozart_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  mozart_Ros2



  SUBROUTINE  mozart_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  mozart_Ros3





  SUBROUTINE  mozart_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  mozart_Ros4


  SUBROUTINE  mozart_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  mozart_Rodas3


  SUBROUTINE  mozart_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  mozart_Rodas4




END SUBROUTINE  mozart_Rosenbrock




SUBROUTINE  mozart_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE mozart_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL mozart_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  mozart_FunTemplate



SUBROUTINE  mozart_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE mozart_Parameters
 
 USE mozart_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL mozart_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  mozart_JacTemplate

















SUBROUTINE mozart_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*F(2)
  A(2) = RCT(2)*V(81)
  A(3) = RCT(3)*V(81)
  A(4) = RCT(4)*V(3)
  A(5) = RCT(5)*V(73)
  A(6) = RCT(6)*V(13)
  A(7) = RCT(7)*V(39)
  A(8) = RCT(8)*V(80)
  A(9) = RCT(9)*V(26)
  A(10) = RCT(10)*V(18)
  A(11) = RCT(11)*V(72)
  A(12) = RCT(12)*V(72)
  A(13) = RCT(13)*V(9)
  A(14) = RCT(14)*V(62)
  A(15) = RCT(15)*V(35)
  A(16) = RCT(16)*V(32)
  A(17) = RCT(17)*V(36)
  A(18) = RCT(18)*V(43)
  A(19) = RCT(19)*V(66)
  A(20) = RCT(20)*V(69)
  A(21) = RCT(21)*V(21)
  A(22) = RCT(22)*V(22)
  A(23) = RCT(23)*V(23)
  A(24) = RCT(24)*V(51)
  A(25) = RCT(25)*V(67)
  A(26) = RCT(26)*V(15)
  A(27) = RCT(27)*V(60)
  A(28) = RCT(28)*V(42)
  A(29) = RCT(29)*V(65)
  A(30) = RCT(30)*V(53)
  A(31) = RCT(31)*V(33)
  A(32) = RCT(32)*V(45)
  A(33) = RCT(33)*V(37)
  A(34) = RCT(34)*V(48)
  A(35) = RCT(35)*V(25)
  A(36) = RCT(36)*V(31)
  A(37) = RCT(37)*V(38)
  A(38) = RCT(38)*V(52)*F(2)
  A(39) = RCT(39)*V(52)*V(81)
  A(40) = RCT(40)*V(47)*F(2)
  A(41) = RCT(41)*V(47)*F(1)
  A(42) = RCT(42)*V(19)*V(47)
  A(43) = RCT(43)*V(19)*V(75)
  A(44) = RCT(44)*V(52)*V(75)
  A(45) = RCT(45)*V(52)*V(79)
  A(46) = RCT(46)*V(75)*V(81)
  A(47) = RCT(47)*V(79)*V(81)
  A(48) = RCT(48)*V(79)*V(79)*F(1)
  A(49) = RCT(49)*V(9)*V(75)
  A(50) = RCT(50)*V(75)*V(79)
  A(51) = RCT(51)*V(75)*V(75)
  A(52) = RCT(52)*V(75)*V(75)
  A(53) = RCT(53)*V(3)*V(47)
  A(54) = RCT(54)*V(3)*V(47)
  A(55) = RCT(55)*V(76)*V(79)
  A(56) = RCT(56)*V(76)*V(81)
  A(57) = RCT(57)*V(52)*V(73)
  A(58) = RCT(58)*V(73)*V(81)
  A(59) = RCT(59)*V(79)*V(80)
  A(60) = RCT(60)*V(73)*V(80)
  A(61) = RCT(61)*V(13)
  A(62) = RCT(62)*V(73)*V(75)
  A(63) = RCT(63)*V(39)*V(75)
  A(64) = RCT(64)*V(76)*V(80)
  A(65) = RCT(65)*V(73)*V(79)
  A(66) = RCT(66)*V(26)*V(75)
  A(67) = RCT(67)*V(26)
  A(68) = RCT(68)*V(13)*F(2)
  A(69) = RCT(69)*V(80)
  A(70) = RCT(70)*V(73)
  A(71) = RCT(71)*V(40)*V(75)
  A(72) = RCT(72)*V(40)*V(47)
  A(73) = RCT(73)*V(76)*V(77)
  A(74) = RCT(74)*V(77)*V(77)
  A(75) = RCT(75)*V(77)*V(77)
  A(76) = RCT(76)*V(77)*V(79)
  A(77) = RCT(77)*V(18)*V(75)
  A(78) = RCT(78)*V(72)*V(80)
  A(79) = RCT(79)*V(72)*V(75)
  A(80) = RCT(80)*V(58)*V(75)
  A(81) = RCT(81)*V(27)*V(75)
  A(82) = RCT(82)*V(27)*V(81)
  A(83) = RCT(83)*V(8)*V(75)
  A(84) = RCT(84)*V(37)*V(75)
  A(85) = RCT(85)*V(30)*V(76)
  A(86) = RCT(86)*V(12)*F(2)
  A(87) = RCT(87)*V(12)
  A(88) = RCT(88)*V(4)*V(75)
  A(89) = RCT(89)*V(55)*V(76)
  A(90) = RCT(90)*V(55)*V(79)
  A(91) = RCT(91)*V(55)*V(77)
  A(92) = RCT(92)*V(21)*V(75)
  A(93) = RCT(93)*V(61)*V(75)
  A(94) = RCT(94)*V(61)*V(81)
  A(95) = RCT(95)*V(61)*V(80)
  A(96) = RCT(96)*V(50)*V(76)
  A(97) = RCT(97)*V(50)*V(79)
  A(98) = RCT(98)*V(35)*V(75)
  A(99) = RCT(99)*V(62)*V(75)
  A(100) = RCT(100)*V(62)*V(80)
  A(101) = RCT(101)*V(76)*V(78)
  A(102) = RCT(102)*V(73)*V(78)
  A(103) = RCT(103)*V(78)*V(79)
  A(104) = RCT(104)*V(77)*V(78)
  A(105) = RCT(105)*V(32)*V(75)
  A(106) = RCT(106)*V(36)
  A(107) = RCT(107)*V(78)*V(78)
  A(108) = RCT(108)*V(5)*V(75)
  A(109) = RCT(109)*V(57)*V(76)
  A(110) = RCT(110)*V(57)*V(79)
  A(111) = RCT(111)*V(57)*V(77)
  A(112) = RCT(112)*V(22)*V(75)
  A(113) = RCT(113)*V(51)*V(75)
  A(114) = RCT(114)*V(63)*V(76)
  A(115) = RCT(115)*V(63)*V(79)
  A(116) = RCT(116)*V(63)*V(77)
  A(117) = RCT(117)*V(23)*V(75)
  A(118) = RCT(118)*V(6)*V(75)
  A(119) = RCT(119)*V(24)*V(76)
  A(120) = RCT(120)*V(7)*V(75)
  A(121) = RCT(121)*V(54)*V(76)
  A(122) = RCT(122)*V(54)*V(79)
  A(123) = RCT(123)*V(48)*V(75)
  A(124) = RCT(124)*V(29)*V(75)
  A(125) = RCT(125)*V(33)*V(75)
  A(126) = RCT(126)*V(49)*V(76)
  A(127) = RCT(127)*V(49)*V(79)
  A(128) = RCT(128)*V(25)*V(75)
  A(129) = RCT(129)*V(10)*V(75)
  A(130) = RCT(130)*V(11)*V(75)
  A(131) = RCT(131)*V(14)*V(73)
  A(132) = RCT(132)*V(44)*V(76)
  A(133) = RCT(133)*V(44)*V(79)
  A(134) = RCT(134)*V(31)*V(75)
  A(135) = RCT(135)*V(59)*V(75)
  A(136) = RCT(136)*V(59)*V(81)
  A(137) = RCT(137)*V(71)*V(76)
  A(138) = RCT(138)*V(71)*V(80)
  A(139) = RCT(139)*V(71)*V(79)
  A(140) = RCT(140)*V(42)*V(75)
  A(141) = RCT(141)*V(71)*V(77)
  A(142) = RCT(142)*V(71)*V(78)
  A(143) = RCT(143)*V(69)*V(75)
  A(144) = RCT(144)*V(69)*V(81)
  A(145) = RCT(145)*V(66)*V(75)
  A(146) = RCT(146)*V(66)*V(81)
  A(147) = RCT(147)*V(70)*V(76)
  A(148) = RCT(148)*V(70)*V(76)
  A(149) = RCT(149)*V(70)*V(80)
  A(150) = RCT(150)*V(70)*V(79)
  A(151) = RCT(151)*V(70)*V(77)
  A(152) = RCT(152)*V(70)*V(78)
  A(153) = RCT(153)*V(20)*V(75)
  A(154) = RCT(154)*V(74)*V(76)
  A(155) = RCT(155)*V(74)*V(80)
  A(156) = RCT(156)*V(74)*V(79)
  A(157) = RCT(157)*V(74)*V(77)
  A(158) = RCT(158)*V(74)*V(78)
  A(159) = RCT(159)*V(74)*V(74)
  A(160) = RCT(160)*V(73)*V(74)*F(2)
  A(161) = RCT(161)*V(43)*F(2)
  A(162) = RCT(162)*V(46)*V(75)
  A(163) = RCT(163)*V(46)*V(81)
  A(164) = RCT(164)*V(46)*V(80)
  A(165) = RCT(165)*V(64)*V(76)
  A(166) = RCT(166)*V(64)*V(79)
  A(167) = RCT(167)*V(38)*V(75)
  A(168) = RCT(168)*V(34)*V(75)
  A(169) = RCT(169)*V(59)*V(80)
  A(170) = RCT(170)*V(56)*V(76)
  A(171) = RCT(171)*V(56)*V(80)
  A(172) = RCT(172)*V(56)*V(79)
  A(173) = RCT(173)*V(67)*V(75)
  A(174) = RCT(174)*V(67)*V(80)
  A(175) = RCT(175)*V(60)*V(75)
  A(176) = RCT(176)*V(60)*V(80)
  A(177) = RCT(177)*V(28)*V(75)
  A(178) = RCT(178)*V(68)*V(76)
  A(179) = RCT(179)*V(68)*V(80)
  A(180) = RCT(180)*V(68)*V(79)
  A(181) = RCT(181)*V(68)*V(77)
  A(182) = RCT(182)*V(68)*V(78)
  A(183) = RCT(183)*V(15)*V(75)
  A(184) = RCT(184)*V(15)*V(75)
  A(185) = RCT(185)*V(41)*V(75)
  A(186) = RCT(186)*V(17)*V(75)
  A(187) = RCT(187)*V(43)*V(75)
  A(188) = RCT(188)*V(36)*V(75)
  A(189) = RCT(189)*V(65)*V(75)
  A(190) = RCT(190)*V(53)*V(75)
  A(191) = RCT(191)*V(16)*V(75)
  A(192) = RCT(192)*V(16)*V(75)
  A(193) = RCT(193)*V(16)*V(80)
  A(194) = RCT(194)*V(2)*V(75)
  A(195) = RCT(195)*V(79)
  A(196) = RCT(196)*V(55)*V(55)


  Vdot(1) = A(83)
  Vdot(2) = -A(194)
  Vdot(3) = -A(4)-A(53)-A(54)
  Vdot(4) = -A(88)
  Vdot(5) = -A(108)
  Vdot(6) = -A(118)
  Vdot(7) = -A(120)
  Vdot(8) = -A(83)+A(191)+0.5*A(192)+A(193)
  Vdot(9) = -A(13)+A(48)-A(49)+A(52)+0.5*A(195)
  Vdot(10) = -A(129)
  Vdot(11) = 0.25*A(129)-A(130)
  Vdot(12) = A(85)-A(86)-A(87)
  Vdot(13) = -A(6)+A(60)-A(61)-A(68)
  Vdot(14) = A(130)-A(131)
  Vdot(15) = -A(26)+A(180)-A(183)-A(184)
  Vdot(16) = -A(191)-A(192)-A(193)
  Vdot(17) = 0.2*A(91)-A(186)+0.4*A(196)
  Vdot(18) = -A(10)+A(76)-A(77)
  Vdot(19) = A(12)-A(42)-A(43)+0.05*A(72)
  Vdot(20) = A(150)-A(153)
  Vdot(21) = -A(21)+A(90)-A(92)
  Vdot(22) = -A(22)+A(110)-A(112)
  Vdot(23) = -A(23)+A(115)-A(117)
  Vdot(24) = A(118)-A(119)
  Vdot(25) = -A(35)+A(127)-A(128)
  Vdot(26) = -A(9)+A(65)-A(66)-A(67)
  Vdot(27) = -A(81)-A(82)
  Vdot(28) = 0.37*A(137)+0.4*A(138)+0.3*A(141)+0.4*A(142)+A(175)+A(176)-A(177)
  Vdot(29) = A(95)+0.1*A(121)-A(124)
  Vdot(30) = 0.75*A(81)-A(85)
  Vdot(31) = -A(36)+A(133)-A(134)
  Vdot(32) = -A(16)+0.75*A(103)-A(105)+0.75*A(156)
  Vdot(33) = -A(31)+0.8*A(34)+0.75*A(121)-A(125)
  Vdot(34) = 0.25*A(82)+0.25*A(94)+0.25*A(103)+0.1*A(104)+0.2*A(136)+0.25*A(156)-A(168)
  Vdot(35) = -A(15)+A(97)-A(98)
  Vdot(36) = -A(17)+A(102)-A(106)-A(188)
  Vdot(37) = 0.13*A(32)-A(33)+0.45*A(36)-A(84)+0.45*A(132)+0.2*A(190)
  Vdot(38) = -A(37)+A(166)-A(167)
  Vdot(39) = -A(7)+A(62)-A(63)+2*A(68)+A(69)+0.5*A(70)+A(78)+A(100)+A(174)+A(193)
  Vdot(40) = -A(71)-A(72)+0.08*A(94)
  Vdot(41) = A(75)+0.3*A(91)+0.5*A(116)+0.25*A(141)+0.25*A(151)+0.3*A(181)-A(185)
  Vdot(42) = -A(28)+A(139)-A(140)
  Vdot(43) = -A(18)+A(160)-A(161)-A(187)
  Vdot(44) = 0.7*A(129)-A(132)-A(133)+A(134)
  Vdot(45) = -A(32)+0.9*A(36)+0.7*A(131)+0.9*A(132)
  Vdot(46) = -A(162)-A(163)-A(164)
  Vdot(47) = A(2)+A(4)-A(40)-A(41)-A(42)-A(53)-A(54)-A(72)
  Vdot(48) = -A(34)+A(122)-A(123)
  Vdot(49) = A(125)-A(126)-A(127)+A(128)
  Vdot(50) = A(93)-A(96)-A(97)+0.5*A(98)
  Vdot(51) = 0.82*A(22)-A(24)+0.25*A(34)+0.1*A(37)+0.82*A(109)+0.82*A(111)-A(113)+0.5*A(119)+0.25*A(121)+0.1*A(165)
  Vdot(52) = 2*A(1)+A(3)+A(5)-A(38)-A(39)+A(40)-A(44)-A(45)+A(51)-A(57)
  Vdot(53) = -A(30)+A(86)+0.53*A(147)+0.53*A(149)+0.26*A(151)+0.53*A(152)+0.25*A(178)+0.25*A(179)+0.1*A(181)+0.25*A(182)&
               &-A(190)
  Vdot(54) = A(120)-A(121)-A(122)+A(123)
  Vdot(55) = A(31)+A(88)-A(89)-A(90)-A(91)+0.5*A(92)-2*A(196)
  Vdot(56) = A(169)-A(170)-A(171)-A(172)
  Vdot(57) = A(108)-A(109)-A(110)-A(111)+A(112)
  Vdot(58) = A(11)+A(12)+A(14)+0.67*A(19)+0.7*A(20)+A(25)+A(27)+A(30)+0.45*A(32)+2*A(33)+A(78)+A(79)-A(80)+0.5*A(82)&
               &+A(84)+0.56*A(94)+0.3*A(136)+0.05*A(144)+0.2*A(146)+0.22*A(147)+0.22*A(149)+0.11*A(151)+0.22*A(152)+A(173)&
               &+A(174)+A(178)+A(179)+0.4*A(181)+A(182)
  Vdot(59) = -A(135)-A(136)-A(169)
  Vdot(60) = -A(27)+0.08*A(137)+0.8*A(148)+0.794*A(170)+0.794*A(171)+0.794*A(172)-A(175)-A(176)
  Vdot(61) = 0.7*A(20)-A(93)-A(94)-A(95)+0.07*A(136)
  Vdot(62) = -A(14)+A(15)+A(21)+0.4*A(34)+A(35)+A(89)+0.8*A(91)+0.5*A(92)+0.5*A(94)+A(96)-A(99)-A(100)+0.27*A(109)&
               &+A(119)+0.4*A(121)+A(126)+0.04*A(144)+A(186)+1.6*A(196)
  Vdot(63) = A(113)-A(114)-A(115)-A(116)+A(117)
  Vdot(64) = A(162)+A(164)-A(165)-A(166)+A(167)
  Vdot(65) = -A(29)+0.5*A(98)+0.2*A(116)+0.22*A(147)+0.22*A(149)+0.23*A(151)+0.22*A(152)+0.25*A(178)+0.25*A(179)+0.1&
               &*A(181)+0.25*A(182)+0.5*A(187)-A(189)
  Vdot(66) = -A(19)+0.288*A(28)+A(37)+0.4*A(136)+0.23*A(137)+0.25*A(138)+0.19*A(141)+0.25*A(142)-A(145)-A(146)+A(163)&
               &+A(165)+0.167*A(170)+0.167*A(171)+0.167*A(172)
  Vdot(67) = -A(25)+0.18*A(32)+0.45*A(36)+0.5*A(116)+A(124)+0.45*A(132)+0.95*A(144)+0.8*A(146)+0.25*A(147)+0.25*A(149)&
               &+0.24*A(151)+0.25*A(152)-A(173)-A(174)+0.25*A(178)+0.25*A(179)+0.1*A(181)+0.25*A(182)+A(189)
  Vdot(68) = 0.5*A(140)+A(177)-A(178)-A(179)-A(180)-A(181)-A(182)+A(183)
  Vdot(69) = -A(20)+0.402*A(28)+A(37)+0.2*A(136)+0.32*A(137)+0.35*A(138)+0.26*A(141)+0.35*A(142)-A(143)-A(144)+A(163)&
               &+A(165)+0.039*A(170)+0.039*A(171)+0.039*A(172)
  Vdot(70) = A(143)+0.5*A(145)-A(147)-A(148)-A(149)-A(150)-A(151)-A(152)+0.2*A(153)
  Vdot(71) = A(135)-A(137)-A(138)-A(139)+0.5*A(140)-A(141)-A(142)
  Vdot(72) = A(10)-A(11)-A(12)+A(15)+0.67*A(19)+A(23)+A(27)+0.69*A(28)+A(29)+A(30)+0.1*A(34)+0.25*A(72)+A(73)+2*A(74)&
               &+A(75)+0.3*A(77)-A(78)-A(79)+0.5*A(81)+A(82)+2*A(87)+0.7*A(91)+0.54*A(94)+A(96)+A(104)+0.5*A(105)+A(111)&
               &+A(114)+0.8*A(116)+0.5*A(119)+0.1*A(121)+0.6*A(136)+0.55*A(137)+0.6*A(138)+1.2*A(141)+0.6*A(142)+0.8*A(144)&
               &+0.7*A(146)+0.25*A(147)+0.25*A(149)+0.88*A(151)+0.25*A(152)+A(154)+A(155)+2*A(157)+A(158)+2*A(159)+0.072&
               &*A(170)+0.072*A(171)+0.206*A(172)+0.7*A(181)+A(185)+0.5*A(187)+A(188)+0.8*A(190)
  Vdot(73) = -A(5)+A(6)+A(7)+0.89*A(8)+0.66*A(9)+0.6*A(17)+A(18)+A(27)+A(55)+A(56)-A(57)-A(58)+A(59)-A(60)+A(61)-A(62)+2&
               &*A(64)-A(65)+A(66)+A(67)-A(70)+A(73)+A(85)+A(89)+A(96)+A(101)-A(102)+A(106)+A(109)+A(114)+A(119)+0.9*A(121)&
               &+A(124)+A(126)-0.3*A(131)+0.9*A(132)+0.92*A(137)+A(138)+A(147)+A(149)+A(154)+A(155)-A(160)+A(161)+A(164)&
               &+A(165)+1.206*A(170)+1.206*A(171)+0.206*A(172)+0.4*A(175)+A(176)+A(178)+A(179)
  Vdot(74) = A(18)+0.33*A(19)+0.2*A(136)+0.5*A(145)+0.5*A(153)-A(154)-A(155)-A(156)-A(157)-A(158)-2*A(159)-A(160)+A(161)
  Vdot(75) = A(7)+0.33*A(9)+A(10)+2*A(13)+A(15)+A(16)+0.33*A(19)+A(21)+A(22)+A(23)+A(26)+A(34)+A(35)+A(36)+A(37)+2*A(41)&
               &+A(42)-A(43)-A(44)+A(45)-A(46)+A(47)-A(49)-A(50)-2*A(51)-2*A(52)+A(55)+A(59)-A(62)-A(63)-A(66)+0.5*A(70)&
               &-A(71)+0.75*A(72)-0.7*A(77)-A(79)-A(80)-A(81)+0.12*A(82)-A(83)-A(84)-A(88)-0.5*A(92)-A(93)+0.33*A(94)-0.5&
               &*A(98)-A(99)-A(105)-A(108)-A(112)-A(113)-A(117)-A(118)-A(120)-A(123)-A(124)-A(125)-A(128)-A(129)-A(130)&
               &-A(134)-A(135)+0.27*A(136)-A(140)-A(143)+0.08*A(144)-A(145)+0.215*A(146)-0.9*A(153)-A(162)+0.7*A(163)-A(167)&
               &-A(168)+0.206*A(172)-A(173)-A(175)-A(177)-A(183)-A(185)-A(186)-A(187)-A(188)-A(189)-A(190)-A(191)-A(192)&
               &-A(194)
  Vdot(76) = A(5)+0.11*A(8)+2*A(53)-A(55)-A(56)+A(57)-A(64)+0.5*A(70)-A(73)-A(85)-A(89)-A(96)-A(101)-A(109)-A(114)&
               &-A(119)-A(121)-A(126)-A(132)-A(137)-A(147)-A(148)-A(154)-A(165)-A(170)-A(178)
  Vdot(77) = A(14)+A(16)+0.4*A(17)+0.3*A(20)+A(24)+A(71)+0.75*A(72)-A(73)-2*A(74)-2*A(75)-A(76)+0.7*A(77)-A(91)+0.31&
               &*A(94)+A(101)-0.1*A(104)+2*A(107)-A(111)-A(116)-A(141)+A(142)-A(151)+A(152)-A(157)+A(158)+A(168)-A(181)&
               &+A(182)
  Vdot(78) = 0.6*A(17)+0.67*A(19)+0.3*A(20)+A(23)+A(24)+A(25)+A(29)+A(31)+0.13*A(32)+A(35)+A(99)+A(100)-A(101)-A(102)&
               &-A(103)-A(104)+0.5*A(105)+A(106)-2*A(107)+A(114)+0.3*A(116)+A(126)-A(142)+0.53*A(147)+0.53*A(149)+0.26&
               &*A(151)-0.47*A(152)+A(154)+A(155)+A(157)+2*A(159)+A(173)+A(174)-A(182)
  Vdot(79) = 0.66*A(9)+A(10)+2*A(11)+A(14)+A(15)+0.67*A(19)+A(21)+A(22)+A(25)+A(27)+A(28)+A(29)+2*A(30)+0.56*A(32)+2&
               &*A(33)+0.9*A(34)+A(37)+A(42)+A(43)+A(44)-A(45)+A(46)-A(47)-2*A(48)+A(49)-A(50)-A(55)-A(59)-A(65)+A(67)+0.4&
               &*A(72)+A(73)+2*A(74)-A(76)+A(78)+A(79)+A(80)+0.25*A(81)+0.12*A(82)+A(84)+A(86)+A(87)+A(89)-A(90)+A(91)+0.19&
               &*A(94)+A(96)-A(97)-A(103)+0.9*A(104)+A(109)-A(110)+A(111)-A(115)+0.3*A(116)+A(119)+0.9*A(121)-A(122)-A(127)&
               &+0.25*A(129)+0.7*A(131)+0.9*A(132)-A(133)+0.06*A(136)+A(137)+A(138)-A(139)+A(141)+A(142)+0.06*A(144)+0.275&
               &*A(146)+0.47*A(147)+0.47*A(149)-A(150)+0.73*A(151)+0.47*A(152)+0.2*A(153)-A(156)+A(157)+A(163)+A(165)-A(166)&
               &+0.794*A(170)+0.794*A(171)-A(172)+A(175)+A(176)+1.5*A(178)+1.5*A(179)-A(180)+A(181)+1.5*A(182)+A(185)+A(186)&
               &+0.5*A(187)+A(189)+A(190)+0.5*A(192)-A(195)+1.2*A(196)
  Vdot(80) = A(6)-A(8)+0.33*A(9)+0.4*A(17)+A(58)-A(59)-A(60)+A(61)+A(63)-A(64)-A(69)-A(78)-A(95)-A(100)-A(138)-A(149)&
               &-A(155)-A(164)-A(169)-A(171)-A(174)-A(176)-A(179)+0.5*A(187)+A(188)-A(193)
  Vdot(81) = -A(2)-A(3)+0.89*A(8)+A(38)-A(39)-A(46)-A(47)-A(56)-A(58)-A(82)-A(94)+0.25*A(103)-0.9*A(136)-0.8*A(144)-0.8&
               &*A(146)+0.25*A(156)-A(163)
      
END SUBROUTINE mozart_Fun
















SUBROUTINE mozart_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*F(2)
  IRR(2) = RCT(2)*V(81)
  IRR(3) = RCT(3)*V(81)
  IRR(4) = RCT(4)*V(3)
  IRR(5) = RCT(5)*V(73)
  IRR(6) = RCT(6)*V(13)
  IRR(7) = RCT(7)*V(39)
  IRR(8) = RCT(8)*V(80)
  IRR(9) = RCT(9)*V(26)
  IRR(10) = RCT(10)*V(18)
  IRR(11) = RCT(11)*V(72)
  IRR(12) = RCT(12)*V(72)
  IRR(13) = RCT(13)*V(9)
  IRR(14) = RCT(14)*V(62)
  IRR(15) = RCT(15)*V(35)
  IRR(16) = RCT(16)*V(32)
  IRR(17) = RCT(17)*V(36)
  IRR(18) = RCT(18)*V(43)
  IRR(19) = RCT(19)*V(66)
  IRR(20) = RCT(20)*V(69)
  IRR(21) = RCT(21)*V(21)
  IRR(22) = RCT(22)*V(22)
  IRR(23) = RCT(23)*V(23)
  IRR(24) = RCT(24)*V(51)
  IRR(25) = RCT(25)*V(67)
  IRR(26) = RCT(26)*V(15)
  IRR(27) = RCT(27)*V(60)
  IRR(28) = RCT(28)*V(42)
  IRR(29) = RCT(29)*V(65)
  IRR(30) = RCT(30)*V(53)
  IRR(31) = RCT(31)*V(33)
  IRR(32) = RCT(32)*V(45)
  IRR(33) = RCT(33)*V(37)
  IRR(34) = RCT(34)*V(48)
  IRR(35) = RCT(35)*V(25)
  IRR(36) = RCT(36)*V(31)
  IRR(37) = RCT(37)*V(38)
  IRR(38) = RCT(38)*V(52)*F(2)
  IRR(39) = RCT(39)*V(52)*V(81)
  IRR(40) = RCT(40)*V(47)*F(2)
  IRR(41) = RCT(41)*V(47)*F(1)
  IRR(42) = RCT(42)*V(19)*V(47)
  IRR(43) = RCT(43)*V(19)*V(75)
  IRR(44) = RCT(44)*V(52)*V(75)
  IRR(45) = RCT(45)*V(52)*V(79)
  IRR(46) = RCT(46)*V(75)*V(81)
  IRR(47) = RCT(47)*V(79)*V(81)
  IRR(48) = RCT(48)*V(79)*V(79)*F(1)
  IRR(49) = RCT(49)*V(9)*V(75)
  IRR(50) = RCT(50)*V(75)*V(79)
  IRR(51) = RCT(51)*V(75)*V(75)
  IRR(52) = RCT(52)*V(75)*V(75)
  IRR(53) = RCT(53)*V(3)*V(47)
  IRR(54) = RCT(54)*V(3)*V(47)
  IRR(55) = RCT(55)*V(76)*V(79)
  IRR(56) = RCT(56)*V(76)*V(81)
  IRR(57) = RCT(57)*V(52)*V(73)
  IRR(58) = RCT(58)*V(73)*V(81)
  IRR(59) = RCT(59)*V(79)*V(80)
  IRR(60) = RCT(60)*V(73)*V(80)
  IRR(61) = RCT(61)*V(13)
  IRR(62) = RCT(62)*V(73)*V(75)
  IRR(63) = RCT(63)*V(39)*V(75)
  IRR(64) = RCT(64)*V(76)*V(80)
  IRR(65) = RCT(65)*V(73)*V(79)
  IRR(66) = RCT(66)*V(26)*V(75)
  IRR(67) = RCT(67)*V(26)
  IRR(68) = RCT(68)*V(13)*F(2)
  IRR(69) = RCT(69)*V(80)
  IRR(70) = RCT(70)*V(73)
  IRR(71) = RCT(71)*V(40)*V(75)
  IRR(72) = RCT(72)*V(40)*V(47)
  IRR(73) = RCT(73)*V(76)*V(77)
  IRR(74) = RCT(74)*V(77)*V(77)
  IRR(75) = RCT(75)*V(77)*V(77)
  IRR(76) = RCT(76)*V(77)*V(79)
  IRR(77) = RCT(77)*V(18)*V(75)
  IRR(78) = RCT(78)*V(72)*V(80)
  IRR(79) = RCT(79)*V(72)*V(75)
  IRR(80) = RCT(80)*V(58)*V(75)
  IRR(81) = RCT(81)*V(27)*V(75)
  IRR(82) = RCT(82)*V(27)*V(81)
  IRR(83) = RCT(83)*V(8)*V(75)
  IRR(84) = RCT(84)*V(37)*V(75)
  IRR(85) = RCT(85)*V(30)*V(76)
  IRR(86) = RCT(86)*V(12)*F(2)
  IRR(87) = RCT(87)*V(12)
  IRR(88) = RCT(88)*V(4)*V(75)
  IRR(89) = RCT(89)*V(55)*V(76)
  IRR(90) = RCT(90)*V(55)*V(79)
  IRR(91) = RCT(91)*V(55)*V(77)
  IRR(92) = RCT(92)*V(21)*V(75)
  IRR(93) = RCT(93)*V(61)*V(75)
  IRR(94) = RCT(94)*V(61)*V(81)
  IRR(95) = RCT(95)*V(61)*V(80)
  IRR(96) = RCT(96)*V(50)*V(76)
  IRR(97) = RCT(97)*V(50)*V(79)
  IRR(98) = RCT(98)*V(35)*V(75)
  IRR(99) = RCT(99)*V(62)*V(75)
  IRR(100) = RCT(100)*V(62)*V(80)
  IRR(101) = RCT(101)*V(76)*V(78)
  IRR(102) = RCT(102)*V(73)*V(78)
  IRR(103) = RCT(103)*V(78)*V(79)
  IRR(104) = RCT(104)*V(77)*V(78)
  IRR(105) = RCT(105)*V(32)*V(75)
  IRR(106) = RCT(106)*V(36)
  IRR(107) = RCT(107)*V(78)*V(78)
  IRR(108) = RCT(108)*V(5)*V(75)
  IRR(109) = RCT(109)*V(57)*V(76)
  IRR(110) = RCT(110)*V(57)*V(79)
  IRR(111) = RCT(111)*V(57)*V(77)
  IRR(112) = RCT(112)*V(22)*V(75)
  IRR(113) = RCT(113)*V(51)*V(75)
  IRR(114) = RCT(114)*V(63)*V(76)
  IRR(115) = RCT(115)*V(63)*V(79)
  IRR(116) = RCT(116)*V(63)*V(77)
  IRR(117) = RCT(117)*V(23)*V(75)
  IRR(118) = RCT(118)*V(6)*V(75)
  IRR(119) = RCT(119)*V(24)*V(76)
  IRR(120) = RCT(120)*V(7)*V(75)
  IRR(121) = RCT(121)*V(54)*V(76)
  IRR(122) = RCT(122)*V(54)*V(79)
  IRR(123) = RCT(123)*V(48)*V(75)
  IRR(124) = RCT(124)*V(29)*V(75)
  IRR(125) = RCT(125)*V(33)*V(75)
  IRR(126) = RCT(126)*V(49)*V(76)
  IRR(127) = RCT(127)*V(49)*V(79)
  IRR(128) = RCT(128)*V(25)*V(75)
  IRR(129) = RCT(129)*V(10)*V(75)
  IRR(130) = RCT(130)*V(11)*V(75)
  IRR(131) = RCT(131)*V(14)*V(73)
  IRR(132) = RCT(132)*V(44)*V(76)
  IRR(133) = RCT(133)*V(44)*V(79)
  IRR(134) = RCT(134)*V(31)*V(75)
  IRR(135) = RCT(135)*V(59)*V(75)
  IRR(136) = RCT(136)*V(59)*V(81)
  IRR(137) = RCT(137)*V(71)*V(76)
  IRR(138) = RCT(138)*V(71)*V(80)
  IRR(139) = RCT(139)*V(71)*V(79)
  IRR(140) = RCT(140)*V(42)*V(75)
  IRR(141) = RCT(141)*V(71)*V(77)
  IRR(142) = RCT(142)*V(71)*V(78)
  IRR(143) = RCT(143)*V(69)*V(75)
  IRR(144) = RCT(144)*V(69)*V(81)
  IRR(145) = RCT(145)*V(66)*V(75)
  IRR(146) = RCT(146)*V(66)*V(81)
  IRR(147) = RCT(147)*V(70)*V(76)
  IRR(148) = RCT(148)*V(70)*V(76)
  IRR(149) = RCT(149)*V(70)*V(80)
  IRR(150) = RCT(150)*V(70)*V(79)
  IRR(151) = RCT(151)*V(70)*V(77)
  IRR(152) = RCT(152)*V(70)*V(78)
  IRR(153) = RCT(153)*V(20)*V(75)
  IRR(154) = RCT(154)*V(74)*V(76)
  IRR(155) = RCT(155)*V(74)*V(80)
  IRR(156) = RCT(156)*V(74)*V(79)
  IRR(157) = RCT(157)*V(74)*V(77)
  IRR(158) = RCT(158)*V(74)*V(78)
  IRR(159) = RCT(159)*V(74)*V(74)
  IRR(160) = RCT(160)*V(73)*V(74)*F(2)
  IRR(161) = RCT(161)*V(43)*F(2)
  IRR(162) = RCT(162)*V(46)*V(75)
  IRR(163) = RCT(163)*V(46)*V(81)
  IRR(164) = RCT(164)*V(46)*V(80)
  IRR(165) = RCT(165)*V(64)*V(76)
  IRR(166) = RCT(166)*V(64)*V(79)
  IRR(167) = RCT(167)*V(38)*V(75)
  IRR(168) = RCT(168)*V(34)*V(75)
  IRR(169) = RCT(169)*V(59)*V(80)
  IRR(170) = RCT(170)*V(56)*V(76)
  IRR(171) = RCT(171)*V(56)*V(80)
  IRR(172) = RCT(172)*V(56)*V(79)
  IRR(173) = RCT(173)*V(67)*V(75)
  IRR(174) = RCT(174)*V(67)*V(80)
  IRR(175) = RCT(175)*V(60)*V(75)
  IRR(176) = RCT(176)*V(60)*V(80)
  IRR(177) = RCT(177)*V(28)*V(75)
  IRR(178) = RCT(178)*V(68)*V(76)
  IRR(179) = RCT(179)*V(68)*V(80)
  IRR(180) = RCT(180)*V(68)*V(79)
  IRR(181) = RCT(181)*V(68)*V(77)
  IRR(182) = RCT(182)*V(68)*V(78)
  IRR(183) = RCT(183)*V(15)*V(75)
  IRR(184) = RCT(184)*V(15)*V(75)
  IRR(185) = RCT(185)*V(41)*V(75)
  IRR(186) = RCT(186)*V(17)*V(75)
  IRR(187) = RCT(187)*V(43)*V(75)
  IRR(188) = RCT(188)*V(36)*V(75)
  IRR(189) = RCT(189)*V(65)*V(75)
  IRR(190) = RCT(190)*V(53)*V(75)
  IRR(191) = RCT(191)*V(16)*V(75)
  IRR(192) = RCT(192)*V(16)*V(75)
  IRR(193) = RCT(193)*V(16)*V(80)
  IRR(194) = RCT(194)*V(2)*V(75)
  IRR(195) = RCT(195)*V(79)
  IRR(196) = RCT(196)*V(55)*V(55)
      
END SUBROUTINE mozart_IRRFun
















SUBROUTINE mozart_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(342)


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

  B(24) = RCT(24)

  B(25) = RCT(25)

  B(26) = RCT(26)

  B(27) = RCT(27)

  B(28) = RCT(28)

  B(29) = RCT(29)

  B(30) = RCT(30)

  B(31) = RCT(31)

  B(32) = RCT(32)

  B(33) = RCT(33)

  B(34) = RCT(34)

  B(35) = RCT(35)

  B(36) = RCT(36)

  B(37) = RCT(37)

  B(38) = RCT(38)*F(2)

  B(40) = RCT(39)*V(81)

  B(41) = RCT(39)*V(52)

  B(42) = RCT(40)*F(2)

  B(44) = RCT(41)*F(1)

  B(46) = RCT(42)*V(47)

  B(47) = RCT(42)*V(19)

  B(48) = RCT(43)*V(75)

  B(49) = RCT(43)*V(19)

  B(50) = RCT(44)*V(75)

  B(51) = RCT(44)*V(52)

  B(52) = RCT(45)*V(79)

  B(53) = RCT(45)*V(52)

  B(54) = RCT(46)*V(81)

  B(55) = RCT(46)*V(75)

  B(56) = RCT(47)*V(81)

  B(57) = RCT(47)*V(79)

  B(58) = RCT(48)*2*V(79)*F(1)

  B(60) = RCT(49)*V(75)

  B(61) = RCT(49)*V(9)

  B(62) = RCT(50)*V(79)

  B(63) = RCT(50)*V(75)

  B(64) = RCT(51)*2*V(75)

  B(65) = RCT(52)*2*V(75)

  B(66) = RCT(53)*V(47)

  B(67) = RCT(53)*V(3)

  B(68) = RCT(54)*V(47)

  B(69) = RCT(54)*V(3)

  B(70) = RCT(55)*V(79)

  B(71) = RCT(55)*V(76)

  B(72) = RCT(56)*V(81)

  B(73) = RCT(56)*V(76)

  B(74) = RCT(57)*V(73)

  B(75) = RCT(57)*V(52)

  B(76) = RCT(58)*V(81)

  B(77) = RCT(58)*V(73)

  B(78) = RCT(59)*V(80)

  B(79) = RCT(59)*V(79)

  B(80) = RCT(60)*V(80)

  B(81) = RCT(60)*V(73)

  B(82) = RCT(61)

  B(83) = RCT(62)*V(75)

  B(84) = RCT(62)*V(73)

  B(85) = RCT(63)*V(75)

  B(86) = RCT(63)*V(39)

  B(87) = RCT(64)*V(80)

  B(88) = RCT(64)*V(76)

  B(89) = RCT(65)*V(79)

  B(90) = RCT(65)*V(73)

  B(91) = RCT(66)*V(75)

  B(92) = RCT(66)*V(26)

  B(93) = RCT(67)

  B(94) = RCT(68)*F(2)

  B(96) = RCT(69)

  B(97) = RCT(70)

  B(98) = RCT(71)*V(75)

  B(99) = RCT(71)*V(40)

  B(100) = RCT(72)*V(47)

  B(101) = RCT(72)*V(40)

  B(102) = RCT(73)*V(77)

  B(103) = RCT(73)*V(76)

  B(104) = RCT(74)*2*V(77)

  B(105) = RCT(75)*2*V(77)

  B(106) = RCT(76)*V(79)

  B(107) = RCT(76)*V(77)

  B(108) = RCT(77)*V(75)

  B(109) = RCT(77)*V(18)

  B(110) = RCT(78)*V(80)

  B(111) = RCT(78)*V(72)

  B(112) = RCT(79)*V(75)

  B(113) = RCT(79)*V(72)

  B(114) = RCT(80)*V(75)

  B(115) = RCT(80)*V(58)

  B(116) = RCT(81)*V(75)

  B(117) = RCT(81)*V(27)

  B(118) = RCT(82)*V(81)

  B(119) = RCT(82)*V(27)

  B(120) = RCT(83)*V(75)

  B(121) = RCT(83)*V(8)

  B(122) = RCT(84)*V(75)

  B(123) = RCT(84)*V(37)

  B(124) = RCT(85)*V(76)

  B(125) = RCT(85)*V(30)

  B(126) = RCT(86)*F(2)

  B(128) = RCT(87)

  B(129) = RCT(88)*V(75)

  B(130) = RCT(88)*V(4)

  B(131) = RCT(89)*V(76)

  B(132) = RCT(89)*V(55)

  B(133) = RCT(90)*V(79)

  B(134) = RCT(90)*V(55)

  B(135) = RCT(91)*V(77)

  B(136) = RCT(91)*V(55)

  B(137) = RCT(92)*V(75)

  B(138) = RCT(92)*V(21)

  B(139) = RCT(93)*V(75)

  B(140) = RCT(93)*V(61)

  B(141) = RCT(94)*V(81)

  B(142) = RCT(94)*V(61)

  B(143) = RCT(95)*V(80)

  B(144) = RCT(95)*V(61)

  B(145) = RCT(96)*V(76)

  B(146) = RCT(96)*V(50)

  B(147) = RCT(97)*V(79)

  B(148) = RCT(97)*V(50)

  B(149) = RCT(98)*V(75)

  B(150) = RCT(98)*V(35)

  B(151) = RCT(99)*V(75)

  B(152) = RCT(99)*V(62)

  B(153) = RCT(100)*V(80)

  B(154) = RCT(100)*V(62)

  B(155) = RCT(101)*V(78)

  B(156) = RCT(101)*V(76)

  B(157) = RCT(102)*V(78)

  B(158) = RCT(102)*V(73)

  B(159) = RCT(103)*V(79)

  B(160) = RCT(103)*V(78)

  B(161) = RCT(104)*V(78)

  B(162) = RCT(104)*V(77)

  B(163) = RCT(105)*V(75)

  B(164) = RCT(105)*V(32)

  B(165) = RCT(106)

  B(166) = RCT(107)*2*V(78)

  B(167) = RCT(108)*V(75)

  B(168) = RCT(108)*V(5)

  B(169) = RCT(109)*V(76)

  B(170) = RCT(109)*V(57)

  B(171) = RCT(110)*V(79)

  B(172) = RCT(110)*V(57)

  B(173) = RCT(111)*V(77)

  B(174) = RCT(111)*V(57)

  B(175) = RCT(112)*V(75)

  B(176) = RCT(112)*V(22)

  B(177) = RCT(113)*V(75)

  B(178) = RCT(113)*V(51)

  B(179) = RCT(114)*V(76)

  B(180) = RCT(114)*V(63)

  B(181) = RCT(115)*V(79)

  B(182) = RCT(115)*V(63)

  B(183) = RCT(116)*V(77)

  B(184) = RCT(116)*V(63)

  B(185) = RCT(117)*V(75)

  B(186) = RCT(117)*V(23)

  B(187) = RCT(118)*V(75)

  B(188) = RCT(118)*V(6)

  B(189) = RCT(119)*V(76)

  B(190) = RCT(119)*V(24)

  B(191) = RCT(120)*V(75)

  B(192) = RCT(120)*V(7)

  B(193) = RCT(121)*V(76)

  B(194) = RCT(121)*V(54)

  B(195) = RCT(122)*V(79)

  B(196) = RCT(122)*V(54)

  B(197) = RCT(123)*V(75)

  B(198) = RCT(123)*V(48)

  B(199) = RCT(124)*V(75)

  B(200) = RCT(124)*V(29)

  B(201) = RCT(125)*V(75)

  B(202) = RCT(125)*V(33)

  B(203) = RCT(126)*V(76)

  B(204) = RCT(126)*V(49)

  B(205) = RCT(127)*V(79)

  B(206) = RCT(127)*V(49)

  B(207) = RCT(128)*V(75)

  B(208) = RCT(128)*V(25)

  B(209) = RCT(129)*V(75)

  B(210) = RCT(129)*V(10)

  B(211) = RCT(130)*V(75)

  B(212) = RCT(130)*V(11)

  B(213) = RCT(131)*V(73)

  B(214) = RCT(131)*V(14)

  B(215) = RCT(132)*V(76)

  B(216) = RCT(132)*V(44)

  B(217) = RCT(133)*V(79)

  B(218) = RCT(133)*V(44)

  B(219) = RCT(134)*V(75)

  B(220) = RCT(134)*V(31)

  B(221) = RCT(135)*V(75)

  B(222) = RCT(135)*V(59)

  B(223) = RCT(136)*V(81)

  B(224) = RCT(136)*V(59)

  B(225) = RCT(137)*V(76)

  B(226) = RCT(137)*V(71)

  B(227) = RCT(138)*V(80)

  B(228) = RCT(138)*V(71)

  B(229) = RCT(139)*V(79)

  B(230) = RCT(139)*V(71)

  B(231) = RCT(140)*V(75)

  B(232) = RCT(140)*V(42)

  B(233) = RCT(141)*V(77)

  B(234) = RCT(141)*V(71)

  B(235) = RCT(142)*V(78)

  B(236) = RCT(142)*V(71)

  B(237) = RCT(143)*V(75)

  B(238) = RCT(143)*V(69)

  B(239) = RCT(144)*V(81)

  B(240) = RCT(144)*V(69)

  B(241) = RCT(145)*V(75)

  B(242) = RCT(145)*V(66)

  B(243) = RCT(146)*V(81)

  B(244) = RCT(146)*V(66)

  B(245) = RCT(147)*V(76)

  B(246) = RCT(147)*V(70)

  B(247) = RCT(148)*V(76)

  B(248) = RCT(148)*V(70)

  B(249) = RCT(149)*V(80)

  B(250) = RCT(149)*V(70)

  B(251) = RCT(150)*V(79)

  B(252) = RCT(150)*V(70)

  B(253) = RCT(151)*V(77)

  B(254) = RCT(151)*V(70)

  B(255) = RCT(152)*V(78)

  B(256) = RCT(152)*V(70)

  B(257) = RCT(153)*V(75)

  B(258) = RCT(153)*V(20)

  B(259) = RCT(154)*V(76)

  B(260) = RCT(154)*V(74)

  B(261) = RCT(155)*V(80)

  B(262) = RCT(155)*V(74)

  B(263) = RCT(156)*V(79)

  B(264) = RCT(156)*V(74)

  B(265) = RCT(157)*V(77)

  B(266) = RCT(157)*V(74)

  B(267) = RCT(158)*V(78)

  B(268) = RCT(158)*V(74)

  B(269) = RCT(159)*2*V(74)

  B(270) = RCT(160)*V(74)*F(2)

  B(271) = RCT(160)*V(73)*F(2)

  B(273) = RCT(161)*F(2)

  B(275) = RCT(162)*V(75)

  B(276) = RCT(162)*V(46)

  B(277) = RCT(163)*V(81)

  B(278) = RCT(163)*V(46)

  B(279) = RCT(164)*V(80)

  B(280) = RCT(164)*V(46)

  B(281) = RCT(165)*V(76)

  B(282) = RCT(165)*V(64)

  B(283) = RCT(166)*V(79)

  B(284) = RCT(166)*V(64)

  B(285) = RCT(167)*V(75)

  B(286) = RCT(167)*V(38)

  B(287) = RCT(168)*V(75)

  B(288) = RCT(168)*V(34)

  B(289) = RCT(169)*V(80)

  B(290) = RCT(169)*V(59)

  B(291) = RCT(170)*V(76)

  B(292) = RCT(170)*V(56)

  B(293) = RCT(171)*V(80)

  B(294) = RCT(171)*V(56)

  B(295) = RCT(172)*V(79)

  B(296) = RCT(172)*V(56)

  B(297) = RCT(173)*V(75)

  B(298) = RCT(173)*V(67)

  B(299) = RCT(174)*V(80)

  B(300) = RCT(174)*V(67)

  B(301) = RCT(175)*V(75)

  B(302) = RCT(175)*V(60)

  B(303) = RCT(176)*V(80)

  B(304) = RCT(176)*V(60)

  B(305) = RCT(177)*V(75)

  B(306) = RCT(177)*V(28)

  B(307) = RCT(178)*V(76)

  B(308) = RCT(178)*V(68)

  B(309) = RCT(179)*V(80)

  B(310) = RCT(179)*V(68)

  B(311) = RCT(180)*V(79)

  B(312) = RCT(180)*V(68)

  B(313) = RCT(181)*V(77)

  B(314) = RCT(181)*V(68)

  B(315) = RCT(182)*V(78)

  B(316) = RCT(182)*V(68)

  B(317) = RCT(183)*V(75)

  B(318) = RCT(183)*V(15)

  B(319) = RCT(184)*V(75)

  B(320) = RCT(184)*V(15)

  B(321) = RCT(185)*V(75)

  B(322) = RCT(185)*V(41)

  B(323) = RCT(186)*V(75)

  B(324) = RCT(186)*V(17)

  B(325) = RCT(187)*V(75)

  B(326) = RCT(187)*V(43)

  B(327) = RCT(188)*V(75)

  B(328) = RCT(188)*V(36)

  B(329) = RCT(189)*V(75)

  B(330) = RCT(189)*V(65)

  B(331) = RCT(190)*V(75)

  B(332) = RCT(190)*V(53)

  B(333) = RCT(191)*V(75)

  B(334) = RCT(191)*V(16)

  B(335) = RCT(192)*V(75)

  B(336) = RCT(192)*V(16)

  B(337) = RCT(193)*V(80)

  B(338) = RCT(193)*V(16)

  B(339) = RCT(194)*V(75)

  B(340) = RCT(194)*V(2)

  B(341) = RCT(195)

  B(342) = RCT(196)*2*V(55)



  JVS(1) = 0

  JVS(2) = B(120)

  JVS(3) = B(121)

  JVS(4) = -B(339)

  JVS(5) = -B(340)

  JVS(6) = -B(4)-B(66)-B(68)

  JVS(7) = -B(67)-B(69)

  JVS(8) = -B(129)

  JVS(9) = -B(130)

  JVS(10) = -B(167)

  JVS(11) = -B(168)

  JVS(12) = -B(187)

  JVS(13) = -B(188)

  JVS(14) = -B(191)

  JVS(15) = -B(192)

  JVS(16) = -B(120)

  JVS(17) = B(333)+0.5*B(335)+B(337)

  JVS(18) = -B(121)+B(334)+0.5*B(336)

  JVS(19) = B(338)

  JVS(20) = -B(13)-B(60)

  JVS(21) = -B(61)+B(65)

  JVS(22) = B(58)+0.5*B(341)

  JVS(23) = -B(209)

  JVS(24) = -B(210)

  JVS(25) = 0.25*B(209)

  JVS(26) = -B(211)

  JVS(27) = 0.25*B(210)-B(212)

  JVS(28) = -B(126)-B(128)

  JVS(29) = B(124)

  JVS(30) = B(125)

  JVS(31) = -B(6)-B(82)-B(94)

  JVS(32) = B(80)

  JVS(33) = B(81)

  JVS(34) = B(211)

  JVS(35) = -B(213)

  JVS(36) = -B(214)

  JVS(37) = B(212)

  JVS(38) = -B(26)-B(317)-B(319)

  JVS(39) = B(311)

  JVS(40) = -B(318)-B(320)

  JVS(41) = B(312)

  JVS(42) = -B(333)-B(335)-B(337)

  JVS(43) = -B(334)-B(336)

  JVS(44) = -B(338)

  JVS(45) = -B(323)

  JVS(46) = 0.2*B(135)+0.4*B(342)

  JVS(47) = -B(324)

  JVS(48) = 0.2*B(136)

  JVS(49) = -B(10)-B(108)

  JVS(50) = -B(109)

  JVS(51) = B(106)

  JVS(52) = B(107)

  JVS(53) = -B(46)-B(48)

  JVS(54) = 0.05*B(100)

  JVS(55) = -B(47)+0.05*B(101)

  JVS(56) = B(12)

  JVS(57) = -B(49)

  JVS(58) = -B(257)

  JVS(59) = B(251)

  JVS(60) = -B(258)

  JVS(61) = B(252)

  JVS(62) = -B(21)-B(137)

  JVS(63) = B(133)

  JVS(64) = -B(138)

  JVS(65) = B(134)

  JVS(66) = -B(22)-B(175)

  JVS(67) = B(171)

  JVS(68) = -B(176)

  JVS(69) = B(172)

  JVS(70) = -B(23)-B(185)

  JVS(71) = B(181)

  JVS(72) = -B(186)

  JVS(73) = B(182)

  JVS(74) = B(187)

  JVS(75) = -B(189)

  JVS(76) = B(188)

  JVS(77) = -B(190)

  JVS(78) = -B(35)-B(207)

  JVS(79) = B(205)

  JVS(80) = -B(208)

  JVS(81) = B(206)

  JVS(82) = -B(9)-B(91)-B(93)

  JVS(83) = B(89)

  JVS(84) = -B(92)

  JVS(85) = B(90)

  JVS(86) = -B(116)-B(118)

  JVS(87) = -B(117)

  JVS(88) = -B(119)

  JVS(89) = -B(305)

  JVS(90) = B(301)+B(303)

  JVS(91) = 0.37*B(225)+0.4*B(227)+0.3*B(233)+0.4*B(235)

  JVS(92) = B(302)-B(306)

  JVS(93) = 0.37*B(226)

  JVS(94) = 0.3*B(234)

  JVS(95) = 0.4*B(236)

  JVS(96) = 0.4*B(228)+B(304)

  JVS(97) = -B(199)

  JVS(98) = 0.1*B(193)

  JVS(99) = B(143)

  JVS(100) = -B(200)

  JVS(101) = 0.1*B(194)

  JVS(102) = B(144)

  JVS(103) = 0.75*B(116)

  JVS(104) = -B(124)

  JVS(105) = 0.75*B(117)

  JVS(106) = -B(125)

  JVS(107) = 0

  JVS(108) = -B(36)-B(219)

  JVS(109) = B(217)

  JVS(110) = -B(220)

  JVS(111) = B(218)

  JVS(112) = -B(16)-B(163)

  JVS(113) = 0.75*B(263)

  JVS(114) = -B(164)

  JVS(115) = 0.75*B(159)

  JVS(116) = 0.75*B(160)+0.75*B(264)

  JVS(117) = -B(31)-B(201)

  JVS(118) = 0.8*B(34)

  JVS(119) = 0.75*B(193)

  JVS(120) = -B(202)

  JVS(121) = 0.75*B(194)

  JVS(122) = 0.25*B(118)

  JVS(123) = -B(287)

  JVS(124) = 0.2*B(223)

  JVS(125) = 0.25*B(141)

  JVS(126) = 0.25*B(263)

  JVS(127) = -B(288)

  JVS(128) = 0.1*B(161)

  JVS(129) = 0.25*B(159)+0.1*B(162)

  JVS(130) = 0.25*B(160)+0.25*B(264)

  JVS(131) = 0.25*B(119)+0.25*B(142)+0.2*B(224)

  JVS(132) = -B(15)-B(149)

  JVS(133) = B(147)

  JVS(134) = -B(150)

  JVS(135) = B(148)

  JVS(136) = -B(17)-B(165)-B(327)

  JVS(137) = B(157)

  JVS(138) = -B(328)

  JVS(139) = B(158)

  JVS(140) = 0.45*B(36)

  JVS(141) = -B(33)-B(122)

  JVS(142) = 0.45*B(215)

  JVS(143) = 0.13*B(32)

  JVS(144) = 0.2*B(331)

  JVS(145) = -B(123)+0.2*B(332)

  JVS(146) = 0.45*B(216)

  JVS(147) = 0

  JVS(148) = -B(37)-B(285)

  JVS(149) = B(283)

  JVS(150) = -B(286)

  JVS(151) = B(284)

  JVS(152) = 2*B(94)

  JVS(153) = B(337)

  JVS(154) = -B(7)-B(85)

  JVS(155) = B(153)

  JVS(156) = B(299)

  JVS(157) = B(110)

  JVS(158) = B(83)+0.5*B(97)

  JVS(159) = B(84)-B(86)

  JVS(160) = B(96)+B(111)+B(154)+B(300)+B(338)

  JVS(161) = -B(98)-B(100)

  JVS(162) = -B(101)

  JVS(163) = 0.08*B(141)

  JVS(164) = -B(99)

  JVS(165) = 0.08*B(142)

  JVS(166) = -B(321)

  JVS(167) = 0.3*B(135)

  JVS(168) = 0.5*B(183)

  JVS(169) = 0.3*B(313)

  JVS(170) = 0.25*B(253)

  JVS(171) = 0.25*B(233)

  JVS(172) = -B(322)

  JVS(173) = B(105)+0.3*B(136)+0.5*B(184)+0.25*B(234)+0.25*B(254)+0.3*B(314)

  JVS(174) = -B(28)-B(231)

  JVS(175) = B(229)

  JVS(176) = -B(232)

  JVS(177) = B(230)

  JVS(178) = -B(18)-B(273)-B(325)

  JVS(179) = B(270)

  JVS(180) = B(271)

  JVS(181) = -B(326)

  JVS(182) = 0.7*B(209)

  JVS(183) = B(219)

  JVS(184) = -B(215)-B(217)

  JVS(185) = 0.7*B(210)+B(220)

  JVS(186) = -B(216)

  JVS(187) = -B(218)

  JVS(188) = 0.7*B(213)

  JVS(189) = 0.9*B(36)

  JVS(190) = 0.9*B(215)

  JVS(191) = -B(32)

  JVS(192) = 0.7*B(214)

  JVS(193) = 0

  JVS(194) = 0.9*B(216)

  JVS(195) = 0

  JVS(196) = -B(275)-B(277)-B(279)

  JVS(197) = -B(276)

  JVS(198) = -B(280)

  JVS(199) = -B(278)

  JVS(200) = B(4)-B(66)-B(68)

  JVS(201) = -B(46)

  JVS(202) = -B(100)

  JVS(203) = -B(42)-B(44)-B(47)-B(67)-B(69)-B(101)

  JVS(204) = 0

  JVS(205) = 0

  JVS(206) = 0

  JVS(207) = B(2)

  JVS(208) = -B(34)-B(197)

  JVS(209) = B(195)

  JVS(210) = -B(198)

  JVS(211) = B(196)

  JVS(212) = B(207)

  JVS(213) = B(201)

  JVS(214) = 0

  JVS(215) = -B(203)-B(205)

  JVS(216) = 0

  JVS(217) = B(202)+B(208)

  JVS(218) = -B(204)

  JVS(219) = -B(206)

  JVS(220) = 0.5*B(149)

  JVS(221) = -B(145)-B(147)

  JVS(222) = B(139)

  JVS(223) = B(140)+0.5*B(150)

  JVS(224) = -B(146)

  JVS(225) = -B(148)

  JVS(226) = 0.82*B(22)

  JVS(227) = 0.5*B(189)

  JVS(228) = 0.1*B(37)

  JVS(229) = 0.25*B(34)

  JVS(230) = -B(24)-B(177)

  JVS(231) = 0.25*B(193)

  JVS(232) = 0.82*B(169)+0.82*B(173)

  JVS(233) = 0.1*B(281)

  JVS(234) = -B(178)

  JVS(235) = 0.82*B(170)+0.5*B(190)+0.25*B(194)+0.1*B(282)

  JVS(236) = 0.82*B(174)

  JVS(237) = 0

  JVS(238) = B(42)

  JVS(239) = -B(38)-B(40)-B(50)-B(52)-B(74)

  JVS(240) = 0

  JVS(241) = 0

  JVS(242) = B(5)-B(75)

  JVS(243) = -B(51)+B(64)

  JVS(244) = -B(53)

  JVS(245) = B(3)-B(41)

  JVS(246) = B(126)

  JVS(247) = 0

  JVS(248) = -B(30)-B(331)

  JVS(249) = 0.25*B(307)+0.25*B(309)+0.1*B(313)+0.25*B(315)

  JVS(250) = 0.53*B(245)+0.53*B(249)+0.26*B(253)+0.53*B(255)

  JVS(251) = -B(332)

  JVS(252) = 0.53*B(246)+0.25*B(308)

  JVS(253) = 0.26*B(254)+0.1*B(314)

  JVS(254) = 0.53*B(256)+0.25*B(316)

  JVS(255) = 0.53*B(250)+0.25*B(310)

  JVS(256) = 0

  JVS(257) = B(191)

  JVS(258) = B(197)

  JVS(259) = -B(193)-B(195)

  JVS(260) = B(192)+B(198)

  JVS(261) = -B(194)

  JVS(262) = -B(196)

  JVS(263) = B(129)

  JVS(264) = 0.5*B(137)

  JVS(265) = B(31)

  JVS(266) = 0

  JVS(267) = 0

  JVS(268) = -B(131)-B(133)-B(135)-2*B(342)

  JVS(269) = B(130)+0.5*B(138)

  JVS(270) = -B(132)

  JVS(271) = -B(136)

  JVS(272) = -B(134)

  JVS(273) = -B(291)-B(293)-B(295)

  JVS(274) = B(289)

  JVS(275) = -B(292)

  JVS(276) = -B(296)

  JVS(277) = B(290)-B(294)

  JVS(278) = B(167)

  JVS(279) = B(175)

  JVS(280) = -B(169)-B(171)-B(173)

  JVS(281) = B(168)+B(176)

  JVS(282) = -B(170)

  JVS(283) = -B(174)

  JVS(284) = -B(172)

  JVS(285) = 0.5*B(118)

  JVS(286) = 2*B(33)+B(122)

  JVS(287) = 0

  JVS(288) = 0.45*B(32)

  JVS(289) = B(30)

  JVS(290) = -B(114)

  JVS(291) = 0.3*B(223)

  JVS(292) = B(27)

  JVS(293) = 0.56*B(141)

  JVS(294) = B(14)

  JVS(295) = 0.67*B(19)+0.2*B(243)

  JVS(296) = B(25)+B(297)+B(299)

  JVS(297) = B(307)+B(309)+0.4*B(313)+B(315)

  JVS(298) = 0.7*B(20)+0.05*B(239)

  JVS(299) = 0.22*B(245)+0.22*B(249)+0.11*B(253)+0.22*B(255)

  JVS(300) = B(11)+B(12)+B(110)+B(112)

  JVS(301) = 0

  JVS(302) = B(113)-B(115)+B(123)+B(298)

  JVS(303) = 0.22*B(246)+B(308)

  JVS(304) = 0.11*B(254)+0.4*B(314)

  JVS(305) = 0.22*B(256)+B(316)

  JVS(306) = 0

  JVS(307) = B(111)+0.22*B(250)+B(300)+B(310)

  JVS(308) = 0.5*B(119)+0.56*B(142)+0.3*B(224)+0.05*B(240)+0.2*B(244)

  JVS(309) = -B(221)-B(223)-B(289)

  JVS(310) = -B(222)

  JVS(311) = -B(290)

  JVS(312) = -B(224)

  JVS(313) = 0.794*B(291)+0.794*B(293)+0.794*B(295)

  JVS(314) = 0

  JVS(315) = -B(27)-B(301)-B(303)

  JVS(316) = 0.8*B(247)

  JVS(317) = 0.08*B(225)

  JVS(318) = -B(302)

  JVS(319) = 0.08*B(226)+0.8*B(248)+0.794*B(292)

  JVS(320) = 0.794*B(296)

  JVS(321) = 0.794*B(294)-B(304)

  JVS(322) = 0

  JVS(323) = 0.07*B(223)

  JVS(324) = -B(139)-B(141)-B(143)

  JVS(325) = 0.7*B(20)

  JVS(326) = -B(140)

  JVS(327) = -B(144)

  JVS(328) = -B(142)+0.07*B(224)

  JVS(329) = B(323)

  JVS(330) = B(21)+0.5*B(137)

  JVS(331) = B(189)

  JVS(332) = B(35)

  JVS(333) = B(15)

  JVS(334) = 0.4*B(34)

  JVS(335) = B(203)

  JVS(336) = B(145)

  JVS(337) = 0.4*B(193)

  JVS(338) = B(131)+0.8*B(135)+1.6*B(342)

  JVS(339) = 0.27*B(169)

  JVS(340) = 0.5*B(141)

  JVS(341) = -B(14)-B(151)-B(153)

  JVS(342) = 0.04*B(239)

  JVS(343) = 0.5*B(138)-B(152)+B(324)

  JVS(344) = B(132)+B(146)+0.27*B(170)+B(190)+0.4*B(194)+B(204)

  JVS(345) = 0.8*B(136)

  JVS(346) = 0

  JVS(347) = -B(154)

  JVS(348) = 0.5*B(142)+0.04*B(240)

  JVS(349) = B(185)

  JVS(350) = B(177)

  JVS(351) = 0

  JVS(352) = 0

  JVS(353) = -B(179)-B(181)-B(183)

  JVS(354) = 0

  JVS(355) = B(178)+B(186)

  JVS(356) = -B(180)

  JVS(357) = -B(184)

  JVS(358) = -B(182)

  JVS(359) = B(285)

  JVS(360) = B(275)+B(279)

  JVS(361) = -B(281)-B(283)

  JVS(362) = B(276)+B(286)

  JVS(363) = -B(282)

  JVS(364) = -B(284)

  JVS(365) = B(280)

  JVS(366) = 0

  JVS(367) = 0.5*B(149)

  JVS(368) = 0.5*B(325)

  JVS(369) = 0

  JVS(370) = 0

  JVS(371) = 0.2*B(183)

  JVS(372) = 0

  JVS(373) = -B(29)-B(329)

  JVS(374) = 0.25*B(307)+0.25*B(309)+0.1*B(313)+0.25*B(315)

  JVS(375) = 0

  JVS(376) = 0.22*B(245)+0.22*B(249)+0.23*B(253)+0.22*B(255)

  JVS(377) = 0

  JVS(378) = 0

  JVS(379) = 0.5*B(150)+0.5*B(326)-B(330)

  JVS(380) = 0.22*B(246)+0.25*B(308)

  JVS(381) = 0.2*B(184)+0.23*B(254)+0.1*B(314)

  JVS(382) = 0.22*B(256)+0.25*B(316)

  JVS(383) = 0

  JVS(384) = 0.22*B(250)+0.25*B(310)

  JVS(385) = 0

  JVS(386) = B(37)

  JVS(387) = 0.288*B(28)

  JVS(388) = B(277)

  JVS(389) = 0.167*B(291)+0.167*B(293)+0.167*B(295)

  JVS(390) = 0.4*B(223)

  JVS(391) = B(281)

  JVS(392) = -B(19)-B(241)-B(243)

  JVS(393) = 0.23*B(225)+0.25*B(227)+0.19*B(233)+0.25*B(235)

  JVS(394) = -B(242)

  JVS(395) = 0.23*B(226)+B(282)+0.167*B(292)

  JVS(396) = 0.19*B(234)

  JVS(397) = 0.25*B(236)

  JVS(398) = 0.167*B(296)

  JVS(399) = 0.25*B(228)+0.167*B(294)

  JVS(400) = 0.4*B(224)-B(244)+B(278)

  JVS(401) = B(199)

  JVS(402) = 0.45*B(36)

  JVS(403) = 0.45*B(215)

  JVS(404) = 0.18*B(32)

  JVS(405) = 0

  JVS(406) = 0

  JVS(407) = 0.5*B(183)

  JVS(408) = 0

  JVS(409) = B(329)

  JVS(410) = 0.8*B(243)

  JVS(411) = -B(25)-B(297)-B(299)

  JVS(412) = 0.25*B(307)+0.25*B(309)+0.1*B(313)+0.25*B(315)

  JVS(413) = 0.95*B(239)

  JVS(414) = 0.25*B(245)+0.25*B(249)+0.24*B(253)+0.25*B(255)

  JVS(415) = 0

  JVS(416) = 0

  JVS(417) = 0

  JVS(418) = B(200)-B(298)+B(330)

  JVS(419) = 0.45*B(216)+0.25*B(246)+0.25*B(308)

  JVS(420) = 0.5*B(184)+0.24*B(254)+0.1*B(314)

  JVS(421) = 0.25*B(256)+0.25*B(316)

  JVS(422) = 0

  JVS(423) = 0.25*B(250)-B(300)+0.25*B(310)

  JVS(424) = 0.95*B(240)+0.8*B(244)

  JVS(425) = B(317)

  JVS(426) = B(305)

  JVS(427) = 0.5*B(231)

  JVS(428) = 0

  JVS(429) = -B(307)-B(309)-B(311)-B(313)-B(315)

  JVS(430) = 0

  JVS(431) = 0

  JVS(432) = 0.5*B(232)+B(306)+B(318)

  JVS(433) = -B(308)

  JVS(434) = -B(314)

  JVS(435) = -B(316)

  JVS(436) = -B(312)

  JVS(437) = -B(310)

  JVS(438) = 0

  JVS(439) = B(37)

  JVS(440) = 0.402*B(28)

  JVS(441) = B(277)

  JVS(442) = 0.039*B(291)+0.039*B(293)+0.039*B(295)

  JVS(443) = 0.2*B(223)

  JVS(444) = B(281)

  JVS(445) = -B(20)-B(237)-B(239)

  JVS(446) = 0.32*B(225)+0.35*B(227)+0.26*B(233)+0.35*B(235)

  JVS(447) = -B(238)

  JVS(448) = 0.32*B(226)+B(282)+0.039*B(292)

  JVS(449) = 0.26*B(234)

  JVS(450) = 0.35*B(236)

  JVS(451) = 0.039*B(296)

  JVS(452) = 0.35*B(228)+0.039*B(294)

  JVS(453) = 0.2*B(224)-B(240)+B(278)

  JVS(454) = 0.2*B(257)

  JVS(455) = 0.5*B(241)

  JVS(456) = B(237)

  JVS(457) = -B(245)-B(247)-B(249)-B(251)-B(253)-B(255)

  JVS(458) = 0

  JVS(459) = B(238)+0.5*B(242)+0.2*B(258)

  JVS(460) = -B(246)-B(248)

  JVS(461) = -B(254)

  JVS(462) = -B(256)

  JVS(463) = -B(252)

  JVS(464) = -B(250)

  JVS(465) = 0

  JVS(466) = 0.5*B(231)

  JVS(467) = B(221)

  JVS(468) = -B(225)-B(227)-B(229)-B(233)-B(235)

  JVS(469) = B(222)+0.5*B(232)

  JVS(470) = -B(226)

  JVS(471) = -B(234)

  JVS(472) = -B(236)

  JVS(473) = -B(230)

  JVS(474) = -B(228)

  JVS(475) = 0

  JVS(476) = 2*B(128)

  JVS(477) = B(10)+0.3*B(108)

  JVS(478) = B(23)

  JVS(479) = 0.5*B(189)

  JVS(480) = 0.5*B(116)+B(118)

  JVS(481) = 0

  JVS(482) = 0.5*B(163)

  JVS(483) = B(15)

  JVS(484) = B(327)

  JVS(485) = 0.25*B(100)

  JVS(486) = B(321)

  JVS(487) = 0.69*B(28)

  JVS(488) = 0.5*B(325)

  JVS(489) = 0.25*B(101)

  JVS(490) = 0.1*B(34)

  JVS(491) = B(145)

  JVS(492) = B(30)+0.8*B(331)

  JVS(493) = 0.1*B(193)

  JVS(494) = 0.7*B(135)

  JVS(495) = 0.072*B(291)+0.072*B(293)+0.206*B(295)

  JVS(496) = B(173)

  JVS(497) = 0.6*B(223)

  JVS(498) = B(27)

  JVS(499) = 0.54*B(141)

  JVS(500) = B(179)+0.8*B(183)

  JVS(501) = 0

  JVS(502) = B(29)

  JVS(503) = 0.67*B(19)+0.7*B(243)

  JVS(504) = 0.7*B(313)

  JVS(505) = 0.8*B(239)

  JVS(506) = 0.25*B(245)+0.25*B(249)+0.88*B(253)+0.25*B(255)

  JVS(507) = 0.55*B(225)+0.6*B(227)+1.2*B(233)+0.6*B(235)

  JVS(508) = -B(11)-B(12)-B(110)-B(112)

  JVS(509) = 0

  JVS(510) = B(259)+B(261)+2*B(265)+B(267)+2*B(269)

  JVS(511) = 0.3*B(109)-B(113)+0.5*B(117)+0.5*B(164)+B(322)+0.5*B(326)+B(328)+0.8*B(332)

  JVS(512) = B(102)+B(146)+B(180)+0.5*B(190)+0.1*B(194)+0.55*B(226)+0.25*B(246)+B(260)+0.072*B(292)

  JVS(513) = B(103)+2*B(104)+B(105)+0.7*B(136)+B(161)+B(174)+0.8*B(184)+1.2*B(234)+0.88*B(254)+2*B(266)+0.7*B(314)

  JVS(514) = B(162)+0.6*B(236)+0.25*B(256)+B(268)

  JVS(515) = 0.206*B(296)

  JVS(516) = -B(111)+0.6*B(228)+0.25*B(250)+B(262)+0.072*B(294)

  JVS(517) = B(119)+0.54*B(142)+0.6*B(224)+0.8*B(240)+0.7*B(244)

  JVS(518) = B(6)+B(82)

  JVS(519) = -0.3*B(213)

  JVS(520) = B(189)

  JVS(521) = 0.66*B(9)+B(91)+B(93)

  JVS(522) = B(199)

  JVS(523) = B(124)

  JVS(524) = 0.6*B(17)+B(165)

  JVS(525) = B(7)

  JVS(526) = B(18)+B(273)

  JVS(527) = 0.9*B(215)

  JVS(528) = B(279)

  JVS(529) = B(203)

  JVS(530) = B(145)

  JVS(531) = -B(74)

  JVS(532) = 0.9*B(193)

  JVS(533) = B(131)

  JVS(534) = 1.206*B(291)+1.206*B(293)+0.206*B(295)

  JVS(535) = B(169)

  JVS(536) = 0

  JVS(537) = B(27)+0.4*B(301)+B(303)

  JVS(538) = 0

  JVS(539) = 0

  JVS(540) = B(179)

  JVS(541) = B(281)

  JVS(542) = 0

  JVS(543) = B(307)+B(309)

  JVS(544) = 0

  JVS(545) = B(245)+B(249)

  JVS(546) = 0.92*B(225)+B(227)

  JVS(547) = 0

  JVS(548) = -B(5)-B(75)-B(76)-B(80)-B(83)-B(89)-B(97)-B(157)-0.3*B(214)-B(270)

  JVS(549) = B(259)+B(261)-B(271)

  JVS(550) = -B(84)+B(92)+B(200)+0.4*B(302)

  JVS(551) = B(70)+B(72)+2*B(87)+B(102)+B(125)+B(132)+B(146)+B(155)+B(170)+B(180)+B(190)+0.9*B(194)+B(204)+0.9*B(216)&
               &+0.92*B(226)+B(246)+B(260)+B(282)+1.206*B(292)+B(308)

  JVS(552) = B(103)

  JVS(553) = B(156)-B(158)

  JVS(554) = B(71)+B(78)-B(90)+0.206*B(296)

  JVS(555) = 0.89*B(8)+B(79)-B(81)+2*B(88)+B(228)+B(250)+B(262)+B(280)+1.206*B(294)+B(304)+B(310)

  JVS(556) = B(73)-B(77)

  JVS(557) = 0.5*B(257)

  JVS(558) = B(18)+B(273)

  JVS(559) = 0.2*B(223)

  JVS(560) = 0.33*B(19)+0.5*B(241)

  JVS(561) = 0

  JVS(562) = 0

  JVS(563) = -B(270)

  JVS(564) = -B(259)-B(261)-B(263)-B(265)-B(267)-2*B(269)-B(271)

  JVS(565) = 0.5*B(242)+0.5*B(258)

  JVS(566) = -B(260)

  JVS(567) = -B(266)

  JVS(568) = -B(268)

  JVS(569) = -B(264)

  JVS(570) = -B(262)

  JVS(571) = 0.2*B(224)

  JVS(572) = -B(339)

  JVS(573) = -B(129)

  JVS(574) = -B(167)

  JVS(575) = -B(187)

  JVS(576) = -B(191)

  JVS(577) = -B(120)

  JVS(578) = 2*B(13)-B(60)

  JVS(579) = -B(209)

  JVS(580) = -B(211)

  JVS(581) = B(26)-B(317)

  JVS(582) = -B(333)-B(335)

  JVS(583) = -B(323)

  JVS(584) = B(10)-0.7*B(108)

  JVS(585) = B(46)-B(48)

  JVS(586) = -0.9*B(257)

  JVS(587) = B(21)-0.5*B(137)

  JVS(588) = B(22)-B(175)

  JVS(589) = B(23)-B(185)

  JVS(590) = B(35)-B(207)

  JVS(591) = 0.33*B(9)-B(91)

  JVS(592) = -B(116)+0.12*B(118)

  JVS(593) = -B(305)

  JVS(594) = -B(199)

  JVS(595) = B(36)-B(219)

  JVS(596) = B(16)-B(163)

  JVS(597) = -B(201)

  JVS(598) = -B(287)

  JVS(599) = B(15)-0.5*B(149)

  JVS(600) = -B(327)

  JVS(601) = -B(122)

  JVS(602) = B(37)-B(285)

  JVS(603) = B(7)-B(85)

  JVS(604) = -B(98)+0.75*B(100)

  JVS(605) = -B(321)

  JVS(606) = -B(231)

  JVS(607) = -B(325)

  JVS(608) = 0

  JVS(609) = 0

  JVS(610) = -B(275)+0.7*B(277)

  JVS(611) = 2*B(44)+B(47)+0.75*B(101)

  JVS(612) = B(34)-B(197)

  JVS(613) = 0

  JVS(614) = 0

  JVS(615) = -B(177)

  JVS(616) = -B(50)+B(52)

  JVS(617) = -B(331)

  JVS(618) = 0

  JVS(619) = 0

  JVS(620) = 0.206*B(295)

  JVS(621) = 0

  JVS(622) = -B(114)

  JVS(623) = -B(221)+0.27*B(223)

  JVS(624) = -B(301)

  JVS(625) = -B(139)+0.33*B(141)

  JVS(626) = -B(151)

  JVS(627) = 0

  JVS(628) = 0

  JVS(629) = -B(329)

  JVS(630) = 0.33*B(19)-B(241)+0.215*B(243)

  JVS(631) = -B(297)

  JVS(632) = 0

  JVS(633) = -B(237)+0.08*B(239)

  JVS(634) = 0

  JVS(635) = 0

  JVS(636) = -B(112)

  JVS(637) = -B(83)+0.5*B(97)

  JVS(638) = 0

  JVS(639) = -B(49)-B(51)-B(54)-B(61)-B(62)-2*B(64)-2*B(65)-B(84)-B(86)-B(92)-B(99)-0.7*B(109)-B(113)-B(115)-B(117)&
               &-B(121)-B(123)-B(130)-0.5*B(138)-B(140)-0.5*B(150)-B(152)-B(164)-B(168)-B(176)-B(178)-B(186)-B(188)-B(192)&
               &-B(198)-B(200)-B(202)-B(208)-B(210)-B(212)-B(220)-B(222)-B(232)-B(238)-B(242)-0.9*B(258)-B(276)-B(286)&
               &-B(288)-B(298)-B(302)-B(306)-B(318)-B(322)-B(324)-B(326)-B(328)-B(330)-B(332)-B(334)-B(336)-B(340)

  JVS(640) = B(70)

  JVS(641) = 0

  JVS(642) = 0

  JVS(643) = B(53)+B(56)-B(63)+B(71)+B(78)+0.206*B(296)

  JVS(644) = B(79)

  JVS(645) = -B(55)+B(57)+0.12*B(119)+0.33*B(142)+0.27*B(224)+0.08*B(240)+0.215*B(244)+0.7*B(278)

  JVS(646) = 2*B(66)

  JVS(647) = -B(189)

  JVS(648) = -B(124)

  JVS(649) = -B(215)

  JVS(650) = 2*B(67)

  JVS(651) = -B(203)

  JVS(652) = -B(145)

  JVS(653) = B(74)

  JVS(654) = -B(193)

  JVS(655) = -B(131)

  JVS(656) = -B(291)

  JVS(657) = -B(169)

  JVS(658) = 0

  JVS(659) = 0

  JVS(660) = -B(179)

  JVS(661) = -B(281)

  JVS(662) = -B(307)

  JVS(663) = 0

  JVS(664) = -B(245)-B(247)

  JVS(665) = -B(225)

  JVS(666) = 0

  JVS(667) = B(5)+B(75)+0.5*B(97)

  JVS(668) = -B(259)

  JVS(669) = 0

  JVS(670) = -B(70)-B(72)-B(87)-B(102)-B(125)-B(132)-B(146)-B(155)-B(170)-B(180)-B(190)-B(194)-B(204)-B(216)-B(226)&
               &-B(246)-B(248)-B(260)-B(282)-B(292)-B(308)

  JVS(671) = -B(103)

  JVS(672) = -B(156)

  JVS(673) = -B(71)

  JVS(674) = 0.11*B(8)-B(88)

  JVS(675) = -B(73)

  JVS(676) = 0.7*B(108)

  JVS(677) = B(16)

  JVS(678) = B(287)

  JVS(679) = 0.4*B(17)

  JVS(680) = B(98)+0.75*B(100)

  JVS(681) = 0.75*B(101)

  JVS(682) = B(24)

  JVS(683) = 0

  JVS(684) = -B(135)

  JVS(685) = -B(173)

  JVS(686) = 0

  JVS(687) = 0.31*B(141)

  JVS(688) = B(14)

  JVS(689) = -B(183)

  JVS(690) = 0

  JVS(691) = -B(313)+B(315)

  JVS(692) = 0.3*B(20)

  JVS(693) = -B(253)+B(255)

  JVS(694) = -B(233)+B(235)

  JVS(695) = 0

  JVS(696) = 0

  JVS(697) = -B(265)+B(267)

  JVS(698) = B(99)+0.7*B(109)+B(288)

  JVS(699) = -B(102)+B(155)

  JVS(700) = -B(103)-2*B(104)-2*B(105)-B(106)-B(136)-0.1*B(161)-B(174)-B(184)-B(234)-B(254)-B(266)-B(314)

  JVS(701) = B(156)-0.1*B(162)+2*B(166)+B(236)+B(256)+B(268)+B(316)

  JVS(702) = -B(107)

  JVS(703) = 0

  JVS(704) = 0.31*B(142)

  JVS(705) = B(23)

  JVS(706) = B(35)

  JVS(707) = 0.5*B(163)

  JVS(708) = B(31)

  JVS(709) = 0.6*B(17)+B(165)

  JVS(710) = 0.13*B(32)

  JVS(711) = 0

  JVS(712) = B(203)

  JVS(713) = B(24)

  JVS(714) = 0

  JVS(715) = 0

  JVS(716) = B(151)+B(153)

  JVS(717) = B(179)+0.3*B(183)

  JVS(718) = 0

  JVS(719) = B(29)

  JVS(720) = 0.67*B(19)

  JVS(721) = B(25)+B(297)+B(299)

  JVS(722) = -B(315)

  JVS(723) = 0.3*B(20)

  JVS(724) = 0.53*B(245)+0.53*B(249)+0.26*B(253)-0.47*B(255)

  JVS(725) = -B(235)

  JVS(726) = -B(157)

  JVS(727) = B(259)+B(261)+B(265)+2*B(269)

  JVS(728) = B(152)+0.5*B(164)+B(298)

  JVS(729) = -B(155)+B(180)+B(204)+0.53*B(246)+B(260)

  JVS(730) = -B(161)+0.3*B(184)+0.26*B(254)+B(266)

  JVS(731) = -B(156)-B(158)-B(159)-B(162)-2*B(166)-B(236)-0.47*B(256)-B(316)

  JVS(732) = -B(160)

  JVS(733) = B(154)+0.53*B(250)+B(262)+B(300)

  JVS(734) = 0

  JVS(735) = B(60)

  JVS(736) = 0.25*B(209)

  JVS(737) = B(126)+B(128)

  JVS(738) = 0.7*B(213)

  JVS(739) = 0.5*B(335)

  JVS(740) = B(323)

  JVS(741) = B(10)

  JVS(742) = B(46)+B(48)

  JVS(743) = 0.2*B(257)

  JVS(744) = B(21)

  JVS(745) = B(22)

  JVS(746) = B(189)

  JVS(747) = 0.66*B(9)+B(93)

  JVS(748) = 0.25*B(116)+0.12*B(118)

  JVS(749) = 0

  JVS(750) = B(15)

  JVS(751) = 2*B(33)+B(122)

  JVS(752) = B(37)

  JVS(753) = 0.4*B(100)

  JVS(754) = B(321)

  JVS(755) = B(28)

  JVS(756) = 0.5*B(325)

  JVS(757) = 0.9*B(215)-B(217)

  JVS(758) = 0.56*B(32)

  JVS(759) = B(277)

  JVS(760) = B(47)+0.4*B(101)

  JVS(761) = 0.9*B(34)

  JVS(762) = -B(205)

  JVS(763) = B(145)-B(147)

  JVS(764) = B(50)-B(52)

  JVS(765) = 2*B(30)+B(331)

  JVS(766) = 0.9*B(193)-B(195)

  JVS(767) = B(131)-B(133)+B(135)+1.2*B(342)

  JVS(768) = 0.794*B(291)+0.794*B(293)-B(295)

  JVS(769) = B(169)-B(171)+B(173)

  JVS(770) = B(114)

  JVS(771) = 0.06*B(223)

  JVS(772) = B(27)+B(301)+B(303)

  JVS(773) = 0.19*B(141)

  JVS(774) = B(14)

  JVS(775) = -B(181)+0.3*B(183)

  JVS(776) = B(281)-B(283)

  JVS(777) = B(29)+B(329)

  JVS(778) = 0.67*B(19)+0.275*B(243)

  JVS(779) = B(25)

  JVS(780) = 1.5*B(307)+1.5*B(309)-B(311)+B(313)+1.5*B(315)

  JVS(781) = 0.06*B(239)

  JVS(782) = 0.47*B(245)+0.47*B(249)-B(251)+0.73*B(253)+0.47*B(255)

  JVS(783) = B(225)+B(227)-B(229)+B(233)+B(235)

  JVS(784) = 2*B(11)+B(110)+B(112)

  JVS(785) = -B(89)+0.7*B(214)

  JVS(786) = -B(263)+B(265)

  JVS(787) = B(49)+B(51)+B(54)+B(61)-B(62)+B(113)+B(115)+0.25*B(117)+B(123)+0.25*B(210)+0.2*B(258)+B(302)+B(322)+B(324)&
               &+0.5*B(326)+B(330)+B(332)+0.5*B(336)

  JVS(788) = -B(70)+B(102)+B(132)+B(146)+B(170)+B(190)+0.9*B(194)+0.9*B(216)+B(226)+0.47*B(246)+B(282)+0.794*B(292)+1.5&
               &*B(308)

  JVS(789) = B(103)+2*B(104)-B(106)+B(136)+0.9*B(161)+B(174)+0.3*B(184)+B(234)+0.73*B(254)+B(266)+B(314)

  JVS(790) = -B(159)+0.9*B(162)+B(236)+0.47*B(256)+1.5*B(316)

  JVS(791) = -B(53)-B(56)-2*B(58)-B(63)-B(71)-B(78)-B(90)-B(107)-B(134)-B(148)-B(160)-B(172)-B(182)-B(196)-B(206)-B(218)&
               &-B(230)-B(252)-B(264)-B(284)-B(296)-B(312)-B(341)

  JVS(792) = -B(79)+B(111)+B(228)+0.47*B(250)+0.794*B(294)+B(304)+1.5*B(310)

  JVS(793) = B(55)-B(57)+0.12*B(119)+0.19*B(142)+0.06*B(224)+0.06*B(240)+0.275*B(244)+B(278)

  JVS(794) = B(6)+B(82)

  JVS(795) = -B(337)

  JVS(796) = 0.33*B(9)

  JVS(797) = 0.4*B(17)+B(327)

  JVS(798) = B(85)

  JVS(799) = 0.5*B(325)

  JVS(800) = -B(279)

  JVS(801) = -B(293)

  JVS(802) = -B(289)

  JVS(803) = -B(303)

  JVS(804) = -B(143)

  JVS(805) = -B(153)

  JVS(806) = -B(299)

  JVS(807) = -B(309)

  JVS(808) = 0

  JVS(809) = -B(249)

  JVS(810) = -B(227)

  JVS(811) = -B(110)

  JVS(812) = B(76)-B(80)

  JVS(813) = -B(261)

  JVS(814) = B(86)+0.5*B(326)+B(328)

  JVS(815) = -B(87)

  JVS(816) = 0

  JVS(817) = 0

  JVS(818) = -B(78)

  JVS(819) = -B(8)-B(79)-B(81)-B(88)-B(96)-B(111)-B(144)-B(154)-B(228)-B(250)-B(262)-B(280)-B(290)-B(294)-B(300)-B(304)&
               &-B(310)-B(338)

  JVS(820) = B(77)

  JVS(821) = -B(118)

  JVS(822) = -B(277)

  JVS(823) = B(38)-B(40)

  JVS(824) = -0.9*B(223)

  JVS(825) = -B(141)

  JVS(826) = -0.8*B(243)

  JVS(827) = -0.8*B(239)

  JVS(828) = 0

  JVS(829) = 0

  JVS(830) = -B(76)

  JVS(831) = 0.25*B(263)

  JVS(832) = -B(54)

  JVS(833) = -B(72)

  JVS(834) = 0

  JVS(835) = 0.25*B(159)

  JVS(836) = -B(56)+0.25*B(160)+0.25*B(264)

  JVS(837) = 0.89*B(8)

  JVS(838) = -B(2)-B(3)-B(41)-B(55)-B(57)-B(73)-B(77)-B(119)-B(142)-0.9*B(224)-0.8*B(240)-0.8*B(244)-B(278)
      
END SUBROUTINE mozart_Jac_SP














SUBROUTINE mozart_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(838), W(81), a
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
      
END SUBROUTINE mozart_KppDecomp



SUBROUTINE mozart_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(838), W(81), a
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
      
END SUBROUTINE mozart_KppDecompCmplx


SUBROUTINE mozart_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(838), X(81), sum

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
      
END SUBROUTINE mozart_KppSolveIndirect


SUBROUTINE mozart_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(838), X(81), sum

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
      
END SUBROUTINE mozart_KppSolveCmplx













SUBROUTINE mozart_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(11) = X(11)-JVS(25)*X(10)
  X(14) = X(14)-JVS(34)*X(11)
  X(24) = X(24)-JVS(74)*X(6)
  X(30) = X(30)-JVS(103)*X(27)
  X(34) = X(34)-JVS(122)*X(27)
  X(37) = X(37)-JVS(140)*X(31)
  X(39) = X(39)-JVS(152)*X(13)-JVS(153)*X(16)
  X(44) = X(44)-JVS(182)*X(10)-JVS(183)*X(31)
  X(45) = X(45)-JVS(188)*X(14)-JVS(189)*X(31)-JVS(190)*X(44)
  X(47) = X(47)-JVS(200)*X(3)-JVS(201)*X(19)-JVS(202)*X(40)
  X(49) = X(49)-JVS(212)*X(25)-JVS(213)*X(33)-JVS(214)*X(48)
  X(50) = X(50)-JVS(220)*X(35)
  X(51) = X(51)-JVS(226)*X(22)-JVS(227)*X(24)-JVS(228)*X(38)-JVS(229)*X(48)
  X(52) = X(52)-JVS(238)*X(47)
  X(53) = X(53)-JVS(246)*X(12)-JVS(247)*X(30)
  X(54) = X(54)-JVS(257)*X(7)-JVS(258)*X(48)
  X(55) = X(55)-JVS(263)*X(4)-JVS(264)*X(21)-JVS(265)*X(33)-JVS(266)*X(48)-JVS(267)*X(54)
  X(57) = X(57)-JVS(278)*X(5)-JVS(279)*X(22)
  X(58) = X(58)-JVS(285)*X(27)-JVS(286)*X(37)-JVS(287)*X(44)-JVS(288)*X(45)-JVS(289)*X(53)
  X(60) = X(60)-JVS(313)*X(56)-JVS(314)*X(59)
  X(61) = X(61)-JVS(323)*X(59)
  X(62) = X(62)-JVS(329)*X(17)-JVS(330)*X(21)-JVS(331)*X(24)-JVS(332)*X(25)-JVS(333)*X(35)-JVS(334)*X(48)-JVS(335)*X(49)&
            &-JVS(336)*X(50)-JVS(337)*X(54)-JVS(338)*X(55)-JVS(339)*X(57)-JVS(340)*X(61)
  X(63) = X(63)-JVS(349)*X(23)-JVS(350)*X(51)-JVS(351)*X(54)-JVS(352)*X(57)
  X(64) = X(64)-JVS(359)*X(38)-JVS(360)*X(46)
  X(65) = X(65)-JVS(367)*X(35)-JVS(368)*X(43)-JVS(369)*X(50)-JVS(370)*X(61)-JVS(371)*X(63)-JVS(372)*X(64)
  X(66) = X(66)-JVS(386)*X(38)-JVS(387)*X(42)-JVS(388)*X(46)-JVS(389)*X(56)-JVS(390)*X(59)-JVS(391)*X(64)
  X(67) = X(67)-JVS(401)*X(29)-JVS(402)*X(31)-JVS(403)*X(44)-JVS(404)*X(45)-JVS(405)*X(54)-JVS(406)*X(61)-JVS(407)*X(63)&
            &-JVS(408)*X(64)-JVS(409)*X(65)-JVS(410)*X(66)
  X(68) = X(68)-JVS(425)*X(15)-JVS(426)*X(28)-JVS(427)*X(42)-JVS(428)*X(60)
  X(69) = X(69)-JVS(439)*X(38)-JVS(440)*X(42)-JVS(441)*X(46)-JVS(442)*X(56)-JVS(443)*X(59)-JVS(444)*X(64)
  X(70) = X(70)-JVS(454)*X(20)-JVS(455)*X(66)-JVS(456)*X(69)
  X(71) = X(71)-JVS(466)*X(42)-JVS(467)*X(59)
  X(72) = X(72)-JVS(476)*X(12)-JVS(477)*X(18)-JVS(478)*X(23)-JVS(479)*X(24)-JVS(480)*X(27)-JVS(481)*X(30)-JVS(482)*X(32)&
            &-JVS(483)*X(35)-JVS(484)*X(36)-JVS(485)*X(40)-JVS(486)*X(41)-JVS(487)*X(42)-JVS(488)*X(43)-JVS(489)*X(47)&
            &-JVS(490)*X(48)-JVS(491)*X(50)-JVS(492)*X(53)-JVS(493)*X(54)-JVS(494)*X(55)-JVS(495)*X(56)-JVS(496)*X(57)&
            &-JVS(497)*X(59)-JVS(498)*X(60)-JVS(499)*X(61)-JVS(500)*X(63)-JVS(501)*X(64)-JVS(502)*X(65)-JVS(503)*X(66)&
            &-JVS(504)*X(68)-JVS(505)*X(69)-JVS(506)*X(70)-JVS(507)*X(71)
  X(73) = X(73)-JVS(518)*X(13)-JVS(519)*X(14)-JVS(520)*X(24)-JVS(521)*X(26)-JVS(522)*X(29)-JVS(523)*X(30)-JVS(524)*X(36)&
            &-JVS(525)*X(39)-JVS(526)*X(43)-JVS(527)*X(44)-JVS(528)*X(46)-JVS(529)*X(49)-JVS(530)*X(50)-JVS(531)*X(52)&
            &-JVS(532)*X(54)-JVS(533)*X(55)-JVS(534)*X(56)-JVS(535)*X(57)-JVS(536)*X(59)-JVS(537)*X(60)-JVS(538)*X(61)&
            &-JVS(539)*X(62)-JVS(540)*X(63)-JVS(541)*X(64)-JVS(542)*X(67)-JVS(543)*X(68)-JVS(544)*X(69)-JVS(545)*X(70)&
            &-JVS(546)*X(71)-JVS(547)*X(72)
  X(74) = X(74)-JVS(557)*X(20)-JVS(558)*X(43)-JVS(559)*X(59)-JVS(560)*X(66)-JVS(561)*X(70)-JVS(562)*X(71)-JVS(563)*X(73)
  X(75) = X(75)-JVS(572)*X(2)-JVS(573)*X(4)-JVS(574)*X(5)-JVS(575)*X(6)-JVS(576)*X(7)-JVS(577)*X(8)-JVS(578)*X(9)&
            &-JVS(579)*X(10)-JVS(580)*X(11)-JVS(581)*X(15)-JVS(582)*X(16)-JVS(583)*X(17)-JVS(584)*X(18)-JVS(585)*X(19)&
            &-JVS(586)*X(20)-JVS(587)*X(21)-JVS(588)*X(22)-JVS(589)*X(23)-JVS(590)*X(25)-JVS(591)*X(26)-JVS(592)*X(27)&
            &-JVS(593)*X(28)-JVS(594)*X(29)-JVS(595)*X(31)-JVS(596)*X(32)-JVS(597)*X(33)-JVS(598)*X(34)-JVS(599)*X(35)&
            &-JVS(600)*X(36)-JVS(601)*X(37)-JVS(602)*X(38)-JVS(603)*X(39)-JVS(604)*X(40)-JVS(605)*X(41)-JVS(606)*X(42)&
            &-JVS(607)*X(43)-JVS(608)*X(44)-JVS(609)*X(45)-JVS(610)*X(46)-JVS(611)*X(47)-JVS(612)*X(48)-JVS(613)*X(49)&
            &-JVS(614)*X(50)-JVS(615)*X(51)-JVS(616)*X(52)-JVS(617)*X(53)-JVS(618)*X(54)-JVS(619)*X(55)-JVS(620)*X(56)&
            &-JVS(621)*X(57)-JVS(622)*X(58)-JVS(623)*X(59)-JVS(624)*X(60)-JVS(625)*X(61)-JVS(626)*X(62)-JVS(627)*X(63)&
            &-JVS(628)*X(64)-JVS(629)*X(65)-JVS(630)*X(66)-JVS(631)*X(67)-JVS(632)*X(68)-JVS(633)*X(69)-JVS(634)*X(70)&
            &-JVS(635)*X(71)-JVS(636)*X(72)-JVS(637)*X(73)-JVS(638)*X(74)
  X(76) = X(76)-JVS(646)*X(3)-JVS(647)*X(24)-JVS(648)*X(30)-JVS(649)*X(44)-JVS(650)*X(47)-JVS(651)*X(49)-JVS(652)*X(50)&
            &-JVS(653)*X(52)-JVS(654)*X(54)-JVS(655)*X(55)-JVS(656)*X(56)-JVS(657)*X(57)-JVS(658)*X(59)-JVS(659)*X(61)&
            &-JVS(660)*X(63)-JVS(661)*X(64)-JVS(662)*X(68)-JVS(663)*X(69)-JVS(664)*X(70)-JVS(665)*X(71)-JVS(666)*X(72)&
            &-JVS(667)*X(73)-JVS(668)*X(74)-JVS(669)*X(75)
  X(77) = X(77)-JVS(676)*X(18)-JVS(677)*X(32)-JVS(678)*X(34)-JVS(679)*X(36)-JVS(680)*X(40)-JVS(681)*X(47)-JVS(682)*X(51)&
            &-JVS(683)*X(54)-JVS(684)*X(55)-JVS(685)*X(57)-JVS(686)*X(59)-JVS(687)*X(61)-JVS(688)*X(62)-JVS(689)*X(63)&
            &-JVS(690)*X(64)-JVS(691)*X(68)-JVS(692)*X(69)-JVS(693)*X(70)-JVS(694)*X(71)-JVS(695)*X(72)-JVS(696)*X(73)&
            &-JVS(697)*X(74)-JVS(698)*X(75)-JVS(699)*X(76)
  X(78) = X(78)-JVS(705)*X(23)-JVS(706)*X(25)-JVS(707)*X(32)-JVS(708)*X(33)-JVS(709)*X(36)-JVS(710)*X(45)-JVS(711)*X(48)&
            &-JVS(712)*X(49)-JVS(713)*X(51)-JVS(714)*X(54)-JVS(715)*X(57)-JVS(716)*X(62)-JVS(717)*X(63)-JVS(718)*X(64)&
            &-JVS(719)*X(65)-JVS(720)*X(66)-JVS(721)*X(67)-JVS(722)*X(68)-JVS(723)*X(69)-JVS(724)*X(70)-JVS(725)*X(71)&
            &-JVS(726)*X(73)-JVS(727)*X(74)-JVS(728)*X(75)-JVS(729)*X(76)-JVS(730)*X(77)
  X(79) = X(79)-JVS(735)*X(9)-JVS(736)*X(10)-JVS(737)*X(12)-JVS(738)*X(14)-JVS(739)*X(16)-JVS(740)*X(17)-JVS(741)*X(18)&
            &-JVS(742)*X(19)-JVS(743)*X(20)-JVS(744)*X(21)-JVS(745)*X(22)-JVS(746)*X(24)-JVS(747)*X(26)-JVS(748)*X(27)&
            &-JVS(749)*X(30)-JVS(750)*X(35)-JVS(751)*X(37)-JVS(752)*X(38)-JVS(753)*X(40)-JVS(754)*X(41)-JVS(755)*X(42)&
            &-JVS(756)*X(43)-JVS(757)*X(44)-JVS(758)*X(45)-JVS(759)*X(46)-JVS(760)*X(47)-JVS(761)*X(48)-JVS(762)*X(49)&
            &-JVS(763)*X(50)-JVS(764)*X(52)-JVS(765)*X(53)-JVS(766)*X(54)-JVS(767)*X(55)-JVS(768)*X(56)-JVS(769)*X(57)&
            &-JVS(770)*X(58)-JVS(771)*X(59)-JVS(772)*X(60)-JVS(773)*X(61)-JVS(774)*X(62)-JVS(775)*X(63)-JVS(776)*X(64)&
            &-JVS(777)*X(65)-JVS(778)*X(66)-JVS(779)*X(67)-JVS(780)*X(68)-JVS(781)*X(69)-JVS(782)*X(70)-JVS(783)*X(71)&
            &-JVS(784)*X(72)-JVS(785)*X(73)-JVS(786)*X(74)-JVS(787)*X(75)-JVS(788)*X(76)-JVS(789)*X(77)-JVS(790)*X(78)
  X(80) = X(80)-JVS(794)*X(13)-JVS(795)*X(16)-JVS(796)*X(26)-JVS(797)*X(36)-JVS(798)*X(39)-JVS(799)*X(43)-JVS(800)*X(46)&
            &-JVS(801)*X(56)-JVS(802)*X(59)-JVS(803)*X(60)-JVS(804)*X(61)-JVS(805)*X(62)-JVS(806)*X(67)-JVS(807)*X(68)&
            &-JVS(808)*X(69)-JVS(809)*X(70)-JVS(810)*X(71)-JVS(811)*X(72)-JVS(812)*X(73)-JVS(813)*X(74)-JVS(814)*X(75)&
            &-JVS(815)*X(76)-JVS(816)*X(77)-JVS(817)*X(78)-JVS(818)*X(79)
  X(81) = X(81)-JVS(821)*X(27)-JVS(822)*X(46)-JVS(823)*X(52)-JVS(824)*X(59)-JVS(825)*X(61)-JVS(826)*X(66)-JVS(827)*X(69)&
            &-JVS(828)*X(71)-JVS(829)*X(72)-JVS(830)*X(73)-JVS(831)*X(74)-JVS(832)*X(75)-JVS(833)*X(76)-JVS(834)*X(77)&
            &-JVS(835)*X(78)-JVS(836)*X(79)-JVS(837)*X(80)
  X(81) = X(81)/JVS(838)
  X(80) = (X(80)-JVS(820)*X(81))/(JVS(819))
  X(79) = (X(79)-JVS(792)*X(80)-JVS(793)*X(81))/(JVS(791))
  X(78) = (X(78)-JVS(732)*X(79)-JVS(733)*X(80)-JVS(734)*X(81))/(JVS(731))
  X(77) = (X(77)-JVS(701)*X(78)-JVS(702)*X(79)-JVS(703)*X(80)-JVS(704)*X(81))/(JVS(700))
  X(76) = (X(76)-JVS(671)*X(77)-JVS(672)*X(78)-JVS(673)*X(79)-JVS(674)*X(80)-JVS(675)*X(81))/(JVS(670))
  X(75) = (X(75)-JVS(640)*X(76)-JVS(641)*X(77)-JVS(642)*X(78)-JVS(643)*X(79)-JVS(644)*X(80)-JVS(645)*X(81))/(JVS(639))
  X(74) = (X(74)-JVS(565)*X(75)-JVS(566)*X(76)-JVS(567)*X(77)-JVS(568)*X(78)-JVS(569)*X(79)-JVS(570)*X(80)-JVS(571)&
            &*X(81))/(JVS(564))
  X(73) = (X(73)-JVS(549)*X(74)-JVS(550)*X(75)-JVS(551)*X(76)-JVS(552)*X(77)-JVS(553)*X(78)-JVS(554)*X(79)-JVS(555)&
            &*X(80)-JVS(556)*X(81))/(JVS(548))
  X(72) = (X(72)-JVS(509)*X(73)-JVS(510)*X(74)-JVS(511)*X(75)-JVS(512)*X(76)-JVS(513)*X(77)-JVS(514)*X(78)-JVS(515)&
            &*X(79)-JVS(516)*X(80)-JVS(517)*X(81))/(JVS(508))
  X(71) = (X(71)-JVS(469)*X(75)-JVS(470)*X(76)-JVS(471)*X(77)-JVS(472)*X(78)-JVS(473)*X(79)-JVS(474)*X(80)-JVS(475)&
            &*X(81))/(JVS(468))
  X(70) = (X(70)-JVS(458)*X(71)-JVS(459)*X(75)-JVS(460)*X(76)-JVS(461)*X(77)-JVS(462)*X(78)-JVS(463)*X(79)-JVS(464)&
            &*X(80)-JVS(465)*X(81))/(JVS(457))
  X(69) = (X(69)-JVS(446)*X(71)-JVS(447)*X(75)-JVS(448)*X(76)-JVS(449)*X(77)-JVS(450)*X(78)-JVS(451)*X(79)-JVS(452)&
            &*X(80)-JVS(453)*X(81))/(JVS(445))
  X(68) = (X(68)-JVS(430)*X(70)-JVS(431)*X(71)-JVS(432)*X(75)-JVS(433)*X(76)-JVS(434)*X(77)-JVS(435)*X(78)-JVS(436)&
            &*X(79)-JVS(437)*X(80)-JVS(438)*X(81))/(JVS(429))
  X(67) = (X(67)-JVS(412)*X(68)-JVS(413)*X(69)-JVS(414)*X(70)-JVS(415)*X(71)-JVS(416)*X(73)-JVS(417)*X(74)-JVS(418)&
            &*X(75)-JVS(419)*X(76)-JVS(420)*X(77)-JVS(421)*X(78)-JVS(422)*X(79)-JVS(423)*X(80)-JVS(424)*X(81))/(JVS(411))
  X(66) = (X(66)-JVS(393)*X(71)-JVS(394)*X(75)-JVS(395)*X(76)-JVS(396)*X(77)-JVS(397)*X(78)-JVS(398)*X(79)-JVS(399)&
            &*X(80)-JVS(400)*X(81))/(JVS(392))
  X(65) = (X(65)-JVS(374)*X(68)-JVS(375)*X(69)-JVS(376)*X(70)-JVS(377)*X(73)-JVS(378)*X(74)-JVS(379)*X(75)-JVS(380)&
            &*X(76)-JVS(381)*X(77)-JVS(382)*X(78)-JVS(383)*X(79)-JVS(384)*X(80)-JVS(385)*X(81))/(JVS(373))
  X(64) = (X(64)-JVS(362)*X(75)-JVS(363)*X(76)-JVS(364)*X(79)-JVS(365)*X(80)-JVS(366)*X(81))/(JVS(361))
  X(63) = (X(63)-JVS(354)*X(64)-JVS(355)*X(75)-JVS(356)*X(76)-JVS(357)*X(77)-JVS(358)*X(79))/(JVS(353))
  X(62) = (X(62)-JVS(342)*X(69)-JVS(343)*X(75)-JVS(344)*X(76)-JVS(345)*X(77)-JVS(346)*X(79)-JVS(347)*X(80)-JVS(348)&
            &*X(81))/(JVS(341))
  X(61) = (X(61)-JVS(325)*X(69)-JVS(326)*X(75)-JVS(327)*X(80)-JVS(328)*X(81))/(JVS(324))
  X(60) = (X(60)-JVS(316)*X(70)-JVS(317)*X(71)-JVS(318)*X(75)-JVS(319)*X(76)-JVS(320)*X(79)-JVS(321)*X(80)-JVS(322)&
            &*X(81))/(JVS(315))
  X(59) = (X(59)-JVS(310)*X(75)-JVS(311)*X(80)-JVS(312)*X(81))/(JVS(309))
  X(58) = (X(58)-JVS(291)*X(59)-JVS(292)*X(60)-JVS(293)*X(61)-JVS(294)*X(62)-JVS(295)*X(66)-JVS(296)*X(67)-JVS(297)&
            &*X(68)-JVS(298)*X(69)-JVS(299)*X(70)-JVS(300)*X(72)-JVS(301)*X(73)-JVS(302)*X(75)-JVS(303)*X(76)-JVS(304)*X(77)&
            &-JVS(305)*X(78)-JVS(306)*X(79)-JVS(307)*X(80)-JVS(308)*X(81))/(JVS(290))
  X(57) = (X(57)-JVS(281)*X(75)-JVS(282)*X(76)-JVS(283)*X(77)-JVS(284)*X(79))/(JVS(280))
  X(56) = (X(56)-JVS(274)*X(59)-JVS(275)*X(76)-JVS(276)*X(79)-JVS(277)*X(80))/(JVS(273))
  X(55) = (X(55)-JVS(269)*X(75)-JVS(270)*X(76)-JVS(271)*X(77)-JVS(272)*X(79))/(JVS(268))
  X(54) = (X(54)-JVS(260)*X(75)-JVS(261)*X(76)-JVS(262)*X(79))/(JVS(259))
  X(53) = (X(53)-JVS(249)*X(68)-JVS(250)*X(70)-JVS(251)*X(75)-JVS(252)*X(76)-JVS(253)*X(77)-JVS(254)*X(78)-JVS(255)&
            &*X(80)-JVS(256)*X(81))/(JVS(248))
  X(52) = (X(52)-JVS(240)*X(61)-JVS(241)*X(72)-JVS(242)*X(73)-JVS(243)*X(75)-JVS(244)*X(79)-JVS(245)*X(81))/(JVS(239))
  X(51) = (X(51)-JVS(231)*X(54)-JVS(232)*X(57)-JVS(233)*X(64)-JVS(234)*X(75)-JVS(235)*X(76)-JVS(236)*X(77)-JVS(237)&
            &*X(79))/(JVS(230))
  X(50) = (X(50)-JVS(222)*X(61)-JVS(223)*X(75)-JVS(224)*X(76)-JVS(225)*X(79))/(JVS(221))
  X(49) = (X(49)-JVS(216)*X(54)-JVS(217)*X(75)-JVS(218)*X(76)-JVS(219)*X(79))/(JVS(215))
  X(48) = (X(48)-JVS(209)*X(54)-JVS(210)*X(75)-JVS(211)*X(79))/(JVS(208))
  X(47) = (X(47)-JVS(204)*X(61)-JVS(205)*X(72)-JVS(206)*X(75)-JVS(207)*X(81))/(JVS(203))
  X(46) = (X(46)-JVS(197)*X(75)-JVS(198)*X(80)-JVS(199)*X(81))/(JVS(196))
  X(45) = (X(45)-JVS(192)*X(73)-JVS(193)*X(75)-JVS(194)*X(76)-JVS(195)*X(79))/(JVS(191))
  X(44) = (X(44)-JVS(185)*X(75)-JVS(186)*X(76)-JVS(187)*X(79))/(JVS(184))
  X(43) = (X(43)-JVS(179)*X(73)-JVS(180)*X(74)-JVS(181)*X(75))/(JVS(178))
  X(42) = (X(42)-JVS(175)*X(71)-JVS(176)*X(75)-JVS(177)*X(79))/(JVS(174))
  X(41) = (X(41)-JVS(167)*X(55)-JVS(168)*X(63)-JVS(169)*X(68)-JVS(170)*X(70)-JVS(171)*X(71)-JVS(172)*X(75)-JVS(173)&
            &*X(77))/(JVS(166))
  X(40) = (X(40)-JVS(162)*X(47)-JVS(163)*X(61)-JVS(164)*X(75)-JVS(165)*X(81))/(JVS(161))
  X(39) = (X(39)-JVS(155)*X(62)-JVS(156)*X(67)-JVS(157)*X(72)-JVS(158)*X(73)-JVS(159)*X(75)-JVS(160)*X(80))/(JVS(154))
  X(38) = (X(38)-JVS(149)*X(64)-JVS(150)*X(75)-JVS(151)*X(79))/(JVS(148))
  X(37) = (X(37)-JVS(142)*X(44)-JVS(143)*X(45)-JVS(144)*X(53)-JVS(145)*X(75)-JVS(146)*X(76)-JVS(147)*X(79))/(JVS(141))
  X(36) = (X(36)-JVS(137)*X(73)-JVS(138)*X(75)-JVS(139)*X(78))/(JVS(136))
  X(35) = (X(35)-JVS(133)*X(50)-JVS(134)*X(75)-JVS(135)*X(79))/(JVS(132))
  X(34) = (X(34)-JVS(124)*X(59)-JVS(125)*X(61)-JVS(126)*X(74)-JVS(127)*X(75)-JVS(128)*X(77)-JVS(129)*X(78)-JVS(130)&
            &*X(79)-JVS(131)*X(81))/(JVS(123))
  X(33) = (X(33)-JVS(118)*X(48)-JVS(119)*X(54)-JVS(120)*X(75)-JVS(121)*X(76))/(JVS(117))
  X(32) = (X(32)-JVS(113)*X(74)-JVS(114)*X(75)-JVS(115)*X(78)-JVS(116)*X(79))/(JVS(112))
  X(31) = (X(31)-JVS(109)*X(44)-JVS(110)*X(75)-JVS(111)*X(79))/(JVS(108))
  X(30) = (X(30)-JVS(105)*X(75)-JVS(106)*X(76)-JVS(107)*X(81))/(JVS(104))
  X(29) = (X(29)-JVS(98)*X(54)-JVS(99)*X(61)-JVS(100)*X(75)-JVS(101)*X(76)-JVS(102)*X(80))/(JVS(97))
  X(28) = (X(28)-JVS(90)*X(60)-JVS(91)*X(71)-JVS(92)*X(75)-JVS(93)*X(76)-JVS(94)*X(77)-JVS(95)*X(78)-JVS(96)*X(80))&
            &/(JVS(89))
  X(27) = (X(27)-JVS(87)*X(75)-JVS(88)*X(81))/(JVS(86))
  X(26) = (X(26)-JVS(83)*X(73)-JVS(84)*X(75)-JVS(85)*X(79))/(JVS(82))
  X(25) = (X(25)-JVS(79)*X(49)-JVS(80)*X(75)-JVS(81)*X(79))/(JVS(78))
  X(24) = (X(24)-JVS(76)*X(75)-JVS(77)*X(76))/(JVS(75))
  X(23) = (X(23)-JVS(71)*X(63)-JVS(72)*X(75)-JVS(73)*X(79))/(JVS(70))
  X(22) = (X(22)-JVS(67)*X(57)-JVS(68)*X(75)-JVS(69)*X(79))/(JVS(66))
  X(21) = (X(21)-JVS(63)*X(55)-JVS(64)*X(75)-JVS(65)*X(79))/(JVS(62))
  X(20) = (X(20)-JVS(59)*X(70)-JVS(60)*X(75)-JVS(61)*X(79))/(JVS(58))
  X(19) = (X(19)-JVS(54)*X(40)-JVS(55)*X(47)-JVS(56)*X(72)-JVS(57)*X(75))/(JVS(53))
  X(18) = (X(18)-JVS(50)*X(75)-JVS(51)*X(77)-JVS(52)*X(79))/(JVS(49))
  X(17) = (X(17)-JVS(46)*X(55)-JVS(47)*X(75)-JVS(48)*X(77))/(JVS(45))
  X(16) = (X(16)-JVS(43)*X(75)-JVS(44)*X(80))/(JVS(42))
  X(15) = (X(15)-JVS(39)*X(68)-JVS(40)*X(75)-JVS(41)*X(79))/(JVS(38))
  X(14) = (X(14)-JVS(36)*X(73)-JVS(37)*X(75))/(JVS(35))
  X(13) = (X(13)-JVS(32)*X(73)-JVS(33)*X(80))/(JVS(31))
  X(12) = (X(12)-JVS(29)*X(30)-JVS(30)*X(76))/(JVS(28))
  X(11) = (X(11)-JVS(27)*X(75))/(JVS(26))
  X(10) = (X(10)-JVS(24)*X(75))/(JVS(23))
  X(9) = (X(9)-JVS(21)*X(75)-JVS(22)*X(79))/(JVS(20))
  X(8) = (X(8)-JVS(17)*X(16)-JVS(18)*X(75)-JVS(19)*X(80))/(JVS(16))
  X(7) = (X(7)-JVS(15)*X(75))/(JVS(14))
  X(6) = (X(6)-JVS(13)*X(75))/(JVS(12))
  X(5) = (X(5)-JVS(11)*X(75))/(JVS(10))
  X(4) = (X(4)-JVS(9)*X(75))/(JVS(8))
  X(3) = (X(3)-JVS(7)*X(47))/(JVS(6))
  X(2) = (X(2)-JVS(5)*X(75))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(8)-JVS(3)*X(75))/(JVS(1))
      
END SUBROUTINE mozart_KppSolve
























      SUBROUTINE mozart_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE mozart_WCOPY



      SUBROUTINE mozart_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE mozart_WAXPY




      SUBROUTINE mozart_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE mozart_WSCAL


      REAL(kind=dp) FUNCTION mozart_WLAMCH( C )








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
          CALL mozart_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      mozart_WLAMCH = Eps

      END FUNCTION mozart_WLAMCH
     
      SUBROUTINE mozart_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE mozart_WLAMCH_ADD




      SUBROUTINE mozart_SET2ZERO(N,Y)




      
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

      END SUBROUTINE mozart_SET2ZERO



      REAL(kind=dp) FUNCTION mozart_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      mozart_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        mozart_WDOT = mozart_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         mozart_WDOT = mozart_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          mozart_WDOT = mozart_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        mozart_WDOT = mozart_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION mozart_WDOT                                          




   SUBROUTINE decomp_mozart( JVS, IER )
   
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
   W( 8 ) = JVS( 2 )
   W( 75 ) = JVS( 3 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 8 )
  JVS( 3) = W( 75 )
  IF ( ABS(  JVS( 4 )) < TINY(a) ) THEN
         IER = 2                                       
         RETURN
  END IF
   W( 2 ) = JVS( 4 )
   W( 75 ) = JVS( 5 )
  JVS( 4) = W( 2 )
  JVS( 5) = W( 75 )
  IF ( ABS(  JVS( 6 )) < TINY(a) ) THEN
         IER = 3                                       
         RETURN
  END IF
   W( 3 ) = JVS( 6 )
   W( 47 ) = JVS( 7 )
  JVS( 6) = W( 3 )
  JVS( 7) = W( 47 )
  IF ( ABS(  JVS( 8 )) < TINY(a) ) THEN
         IER = 4                                       
         RETURN
  END IF
   W( 4 ) = JVS( 8 )
   W( 75 ) = JVS( 9 )
  JVS( 8) = W( 4 )
  JVS( 9) = W( 75 )
  IF ( ABS(  JVS( 10 )) < TINY(a) ) THEN
         IER = 5                                       
         RETURN
  END IF
   W( 5 ) = JVS( 10 )
   W( 75 ) = JVS( 11 )
  JVS( 10) = W( 5 )
  JVS( 11) = W( 75 )
  IF ( ABS(  JVS( 12 )) < TINY(a) ) THEN
         IER = 6                                       
         RETURN
  END IF
   W( 6 ) = JVS( 12 )
   W( 75 ) = JVS( 13 )
  JVS( 12) = W( 6 )
  JVS( 13) = W( 75 )
  IF ( ABS(  JVS( 14 )) < TINY(a) ) THEN
         IER = 7                                       
         RETURN
  END IF
   W( 7 ) = JVS( 14 )
   W( 75 ) = JVS( 15 )
  JVS( 14) = W( 7 )
  JVS( 15) = W( 75 )
  IF ( ABS(  JVS( 16 )) < TINY(a) ) THEN
         IER = 8                                       
         RETURN
  END IF
   W( 8 ) = JVS( 16 )
   W( 16 ) = JVS( 17 )
   W( 75 ) = JVS( 18 )
   W( 80 ) = JVS( 19 )
  JVS( 16) = W( 8 )
  JVS( 17) = W( 16 )
  JVS( 18) = W( 75 )
  JVS( 19) = W( 80 )
  IF ( ABS(  JVS( 20 )) < TINY(a) ) THEN
         IER = 9                                       
         RETURN
  END IF
   W( 9 ) = JVS( 20 )
   W( 75 ) = JVS( 21 )
   W( 79 ) = JVS( 22 )
  JVS( 20) = W( 9 )
  JVS( 21) = W( 75 )
  JVS( 22) = W( 79 )
  IF ( ABS(  JVS( 23 )) < TINY(a) ) THEN
         IER = 10                                      
         RETURN
  END IF
   W( 10 ) = JVS( 23 )
   W( 75 ) = JVS( 24 )
  JVS( 23) = W( 10 )
  JVS( 24) = W( 75 )
  IF ( ABS(  JVS( 26 )) < TINY(a) ) THEN
         IER = 11                                      
         RETURN
  END IF
   W( 10 ) = JVS( 25 )
   W( 11 ) = JVS( 26 )
   W( 75 ) = JVS( 27 )
  a = -W( 10 ) / JVS(           23  )
  W( 10 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 24 )
  JVS( 25) = W( 10 )
  JVS( 26) = W( 11 )
  JVS( 27) = W( 75 )
  IF ( ABS(  JVS( 28 )) < TINY(a) ) THEN
         IER = 12                                      
         RETURN
  END IF
   W( 12 ) = JVS( 28 )
   W( 30 ) = JVS( 29 )
   W( 76 ) = JVS( 30 )
  JVS( 28) = W( 12 )
  JVS( 29) = W( 30 )
  JVS( 30) = W( 76 )
  IF ( ABS(  JVS( 31 )) < TINY(a) ) THEN
         IER = 13                                      
         RETURN
  END IF
   W( 13 ) = JVS( 31 )
   W( 73 ) = JVS( 32 )
   W( 80 ) = JVS( 33 )
  JVS( 31) = W( 13 )
  JVS( 32) = W( 73 )
  JVS( 33) = W( 80 )
  IF ( ABS(  JVS( 35 )) < TINY(a) ) THEN
         IER = 14                                      
         RETURN
  END IF
   W( 11 ) = JVS( 34 )
   W( 14 ) = JVS( 35 )
   W( 73 ) = JVS( 36 )
   W( 75 ) = JVS( 37 )
  a = -W( 11 ) / JVS(           26  )
  W( 11 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 27 )
  JVS( 34) = W( 11 )
  JVS( 35) = W( 14 )
  JVS( 36) = W( 73 )
  JVS( 37) = W( 75 )
  IF ( ABS(  JVS( 38 )) < TINY(a) ) THEN
         IER = 15                                      
         RETURN
  END IF
   W( 15 ) = JVS( 38 )
   W( 68 ) = JVS( 39 )
   W( 75 ) = JVS( 40 )
   W( 79 ) = JVS( 41 )
  JVS( 38) = W( 15 )
  JVS( 39) = W( 68 )
  JVS( 40) = W( 75 )
  JVS( 41) = W( 79 )
  IF ( ABS(  JVS( 42 )) < TINY(a) ) THEN
         IER = 16                                      
         RETURN
  END IF
   W( 16 ) = JVS( 42 )
   W( 75 ) = JVS( 43 )
   W( 80 ) = JVS( 44 )
  JVS( 42) = W( 16 )
  JVS( 43) = W( 75 )
  JVS( 44) = W( 80 )
  IF ( ABS(  JVS( 45 )) < TINY(a) ) THEN
         IER = 17                                      
         RETURN
  END IF
   W( 17 ) = JVS( 45 )
   W( 55 ) = JVS( 46 )
   W( 75 ) = JVS( 47 )
   W( 77 ) = JVS( 48 )
  JVS( 45) = W( 17 )
  JVS( 46) = W( 55 )
  JVS( 47) = W( 75 )
  JVS( 48) = W( 77 )
  IF ( ABS(  JVS( 49 )) < TINY(a) ) THEN
         IER = 18                                      
         RETURN
  END IF
   W( 18 ) = JVS( 49 )
   W( 75 ) = JVS( 50 )
   W( 77 ) = JVS( 51 )
   W( 79 ) = JVS( 52 )
  JVS( 49) = W( 18 )
  JVS( 50) = W( 75 )
  JVS( 51) = W( 77 )
  JVS( 52) = W( 79 )
  IF ( ABS(  JVS( 53 )) < TINY(a) ) THEN
         IER = 19                                      
         RETURN
  END IF
   W( 19 ) = JVS( 53 )
   W( 40 ) = JVS( 54 )
   W( 47 ) = JVS( 55 )
   W( 72 ) = JVS( 56 )
   W( 75 ) = JVS( 57 )
  JVS( 53) = W( 19 )
  JVS( 54) = W( 40 )
  JVS( 55) = W( 47 )
  JVS( 56) = W( 72 )
  JVS( 57) = W( 75 )
  IF ( ABS(  JVS( 58 )) < TINY(a) ) THEN
         IER = 20                                      
         RETURN
  END IF
   W( 20 ) = JVS( 58 )
   W( 70 ) = JVS( 59 )
   W( 75 ) = JVS( 60 )
   W( 79 ) = JVS( 61 )
  JVS( 58) = W( 20 )
  JVS( 59) = W( 70 )
  JVS( 60) = W( 75 )
  JVS( 61) = W( 79 )
  IF ( ABS(  JVS( 62 )) < TINY(a) ) THEN
         IER = 21                                      
         RETURN
  END IF
   W( 21 ) = JVS( 62 )
   W( 55 ) = JVS( 63 )
   W( 75 ) = JVS( 64 )
   W( 79 ) = JVS( 65 )
  JVS( 62) = W( 21 )
  JVS( 63) = W( 55 )
  JVS( 64) = W( 75 )
  JVS( 65) = W( 79 )
  IF ( ABS(  JVS( 66 )) < TINY(a) ) THEN
         IER = 22                                      
         RETURN
  END IF
   W( 22 ) = JVS( 66 )
   W( 57 ) = JVS( 67 )
   W( 75 ) = JVS( 68 )
   W( 79 ) = JVS( 69 )
  JVS( 66) = W( 22 )
  JVS( 67) = W( 57 )
  JVS( 68) = W( 75 )
  JVS( 69) = W( 79 )
  IF ( ABS(  JVS( 70 )) < TINY(a) ) THEN
         IER = 23                                      
         RETURN
  END IF
   W( 23 ) = JVS( 70 )
   W( 63 ) = JVS( 71 )
   W( 75 ) = JVS( 72 )
   W( 79 ) = JVS( 73 )
  JVS( 70) = W( 23 )
  JVS( 71) = W( 63 )
  JVS( 72) = W( 75 )
  JVS( 73) = W( 79 )
  IF ( ABS(  JVS( 75 )) < TINY(a) ) THEN
         IER = 24                                      
         RETURN
  END IF
   W( 6 ) = JVS( 74 )
   W( 24 ) = JVS( 75 )
   W( 75 ) = JVS( 76 )
   W( 76 ) = JVS( 77 )
  a = -W( 6 ) / JVS(           12  )
  W( 6 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 13 )
  JVS( 74) = W( 6 )
  JVS( 75) = W( 24 )
  JVS( 76) = W( 75 )
  JVS( 77) = W( 76 )
  IF ( ABS(  JVS( 78 )) < TINY(a) ) THEN
         IER = 25                                      
         RETURN
  END IF
   W( 25 ) = JVS( 78 )
   W( 49 ) = JVS( 79 )
   W( 75 ) = JVS( 80 )
   W( 79 ) = JVS( 81 )
  JVS( 78) = W( 25 )
  JVS( 79) = W( 49 )
  JVS( 80) = W( 75 )
  JVS( 81) = W( 79 )
  IF ( ABS(  JVS( 82 )) < TINY(a) ) THEN
         IER = 26                                      
         RETURN
  END IF
   W( 26 ) = JVS( 82 )
   W( 73 ) = JVS( 83 )
   W( 75 ) = JVS( 84 )
   W( 79 ) = JVS( 85 )
  JVS( 82) = W( 26 )
  JVS( 83) = W( 73 )
  JVS( 84) = W( 75 )
  JVS( 85) = W( 79 )
  IF ( ABS(  JVS( 86 )) < TINY(a) ) THEN
         IER = 27                                      
         RETURN
  END IF
   W( 27 ) = JVS( 86 )
   W( 75 ) = JVS( 87 )
   W( 81 ) = JVS( 88 )
  JVS( 86) = W( 27 )
  JVS( 87) = W( 75 )
  JVS( 88) = W( 81 )
  IF ( ABS(  JVS( 89 )) < TINY(a) ) THEN
         IER = 28                                      
         RETURN
  END IF
   W( 28 ) = JVS( 89 )
   W( 60 ) = JVS( 90 )
   W( 71 ) = JVS( 91 )
   W( 75 ) = JVS( 92 )
   W( 76 ) = JVS( 93 )
   W( 77 ) = JVS( 94 )
   W( 78 ) = JVS( 95 )
   W( 80 ) = JVS( 96 )
  JVS( 89) = W( 28 )
  JVS( 90) = W( 60 )
  JVS( 91) = W( 71 )
  JVS( 92) = W( 75 )
  JVS( 93) = W( 76 )
  JVS( 94) = W( 77 )
  JVS( 95) = W( 78 )
  JVS( 96) = W( 80 )
  IF ( ABS(  JVS( 97 )) < TINY(a) ) THEN
         IER = 29                                      
         RETURN
  END IF
   W( 29 ) = JVS( 97 )
   W( 54 ) = JVS( 98 )
   W( 61 ) = JVS( 99 )
   W( 75 ) = JVS( 100 )
   W( 76 ) = JVS( 101 )
   W( 80 ) = JVS( 102 )
  JVS( 97) = W( 29 )
  JVS( 98) = W( 54 )
  JVS( 99) = W( 61 )
  JVS( 100) = W( 75 )
  JVS( 101) = W( 76 )
  JVS( 102) = W( 80 )
  IF ( ABS(  JVS( 104 )) < TINY(a) ) THEN
         IER = 30                                      
         RETURN
  END IF
   W( 27 ) = JVS( 103 )
   W( 30 ) = JVS( 104 )
   W( 75 ) = JVS( 105 )
   W( 76 ) = JVS( 106 )
   W( 81 ) = JVS( 107 )
  a = -W( 27 ) / JVS(           86  )
  W( 27 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 87 )
  W( 81 ) = W( 81 ) + a*JVS( 88 )
  JVS( 103) = W( 27 )
  JVS( 104) = W( 30 )
  JVS( 105) = W( 75 )
  JVS( 106) = W( 76 )
  JVS( 107) = W( 81 )
  IF ( ABS(  JVS( 108 )) < TINY(a) ) THEN
         IER = 31                                      
         RETURN
  END IF
   W( 31 ) = JVS( 108 )
   W( 44 ) = JVS( 109 )
   W( 75 ) = JVS( 110 )
   W( 79 ) = JVS( 111 )
  JVS( 108) = W( 31 )
  JVS( 109) = W( 44 )
  JVS( 110) = W( 75 )
  JVS( 111) = W( 79 )
  IF ( ABS(  JVS( 112 )) < TINY(a) ) THEN
         IER = 32                                      
         RETURN
  END IF
   W( 32 ) = JVS( 112 )
   W( 74 ) = JVS( 113 )
   W( 75 ) = JVS( 114 )
   W( 78 ) = JVS( 115 )
   W( 79 ) = JVS( 116 )
  JVS( 112) = W( 32 )
  JVS( 113) = W( 74 )
  JVS( 114) = W( 75 )
  JVS( 115) = W( 78 )
  JVS( 116) = W( 79 )
  IF ( ABS(  JVS( 117 )) < TINY(a) ) THEN
         IER = 33                                      
         RETURN
  END IF
   W( 33 ) = JVS( 117 )
   W( 48 ) = JVS( 118 )
   W( 54 ) = JVS( 119 )
   W( 75 ) = JVS( 120 )
   W( 76 ) = JVS( 121 )
  JVS( 117) = W( 33 )
  JVS( 118) = W( 48 )
  JVS( 119) = W( 54 )
  JVS( 120) = W( 75 )
  JVS( 121) = W( 76 )
  IF ( ABS(  JVS( 123 )) < TINY(a) ) THEN
         IER = 34                                      
         RETURN
  END IF
   W( 27 ) = JVS( 122 )
   W( 34 ) = JVS( 123 )
   W( 59 ) = JVS( 124 )
   W( 61 ) = JVS( 125 )
   W( 74 ) = JVS( 126 )
   W( 75 ) = JVS( 127 )
   W( 77 ) = JVS( 128 )
   W( 78 ) = JVS( 129 )
   W( 79 ) = JVS( 130 )
   W( 81 ) = JVS( 131 )
  a = -W( 27 ) / JVS(           86  )
  W( 27 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 87 )
  W( 81 ) = W( 81 ) + a*JVS( 88 )
  JVS( 122) = W( 27 )
  JVS( 123) = W( 34 )
  JVS( 124) = W( 59 )
  JVS( 125) = W( 61 )
  JVS( 126) = W( 74 )
  JVS( 127) = W( 75 )
  JVS( 128) = W( 77 )
  JVS( 129) = W( 78 )
  JVS( 130) = W( 79 )
  JVS( 131) = W( 81 )
  IF ( ABS(  JVS( 132 )) < TINY(a) ) THEN
         IER = 35                                      
         RETURN
  END IF
   W( 35 ) = JVS( 132 )
   W( 50 ) = JVS( 133 )
   W( 75 ) = JVS( 134 )
   W( 79 ) = JVS( 135 )
  JVS( 132) = W( 35 )
  JVS( 133) = W( 50 )
  JVS( 134) = W( 75 )
  JVS( 135) = W( 79 )
  IF ( ABS(  JVS( 136 )) < TINY(a) ) THEN
         IER = 36                                      
         RETURN
  END IF
   W( 36 ) = JVS( 136 )
   W( 73 ) = JVS( 137 )
   W( 75 ) = JVS( 138 )
   W( 78 ) = JVS( 139 )
  JVS( 136) = W( 36 )
  JVS( 137) = W( 73 )
  JVS( 138) = W( 75 )
  JVS( 139) = W( 78 )
  IF ( ABS(  JVS( 141 )) < TINY(a) ) THEN
         IER = 37                                      
         RETURN
  END IF
   W( 31 ) = JVS( 140 )
   W( 37 ) = JVS( 141 )
   W( 44 ) = JVS( 142 )
   W( 45 ) = JVS( 143 )
   W( 53 ) = JVS( 144 )
   W( 75 ) = JVS( 145 )
   W( 76 ) = JVS( 146 )
   W( 79 ) = JVS( 147 )
  a = -W( 31 ) / JVS(          108  )
  W( 31 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 109 )
  W( 75 ) = W( 75 ) + a*JVS( 110 )
  W( 79 ) = W( 79 ) + a*JVS( 111 )
  JVS( 140) = W( 31 )
  JVS( 141) = W( 37 )
  JVS( 142) = W( 44 )
  JVS( 143) = W( 45 )
  JVS( 144) = W( 53 )
  JVS( 145) = W( 75 )
  JVS( 146) = W( 76 )
  JVS( 147) = W( 79 )
  IF ( ABS(  JVS( 148 )) < TINY(a) ) THEN
         IER = 38                                      
         RETURN
  END IF
   W( 38 ) = JVS( 148 )
   W( 64 ) = JVS( 149 )
   W( 75 ) = JVS( 150 )
   W( 79 ) = JVS( 151 )
  JVS( 148) = W( 38 )
  JVS( 149) = W( 64 )
  JVS( 150) = W( 75 )
  JVS( 151) = W( 79 )
  IF ( ABS(  JVS( 154 )) < TINY(a) ) THEN
         IER = 39                                      
         RETURN
  END IF
   W( 13 ) = JVS( 152 )
   W( 16 ) = JVS( 153 )
   W( 39 ) = JVS( 154 )
   W( 62 ) = JVS( 155 )
   W( 67 ) = JVS( 156 )
   W( 72 ) = JVS( 157 )
   W( 73 ) = JVS( 158 )
   W( 75 ) = JVS( 159 )
   W( 80 ) = JVS( 160 )
  a = -W( 13 ) / JVS(           31  )
  W( 13 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 32 )
  W( 80 ) = W( 80 ) + a*JVS( 33 )
  a = -W( 16 ) / JVS(           42  )
  W( 16 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 43 )
  W( 80 ) = W( 80 ) + a*JVS( 44 )
  JVS( 152) = W( 13 )
  JVS( 153) = W( 16 )
  JVS( 154) = W( 39 )
  JVS( 155) = W( 62 )
  JVS( 156) = W( 67 )
  JVS( 157) = W( 72 )
  JVS( 158) = W( 73 )
  JVS( 159) = W( 75 )
  JVS( 160) = W( 80 )
  IF ( ABS(  JVS( 161 )) < TINY(a) ) THEN
         IER = 40                                      
         RETURN
  END IF
   W( 40 ) = JVS( 161 )
   W( 47 ) = JVS( 162 )
   W( 61 ) = JVS( 163 )
   W( 75 ) = JVS( 164 )
   W( 81 ) = JVS( 165 )
  JVS( 161) = W( 40 )
  JVS( 162) = W( 47 )
  JVS( 163) = W( 61 )
  JVS( 164) = W( 75 )
  JVS( 165) = W( 81 )
  IF ( ABS(  JVS( 166 )) < TINY(a) ) THEN
         IER = 41                                      
         RETURN
  END IF
   W( 41 ) = JVS( 166 )
   W( 55 ) = JVS( 167 )
   W( 63 ) = JVS( 168 )
   W( 68 ) = JVS( 169 )
   W( 70 ) = JVS( 170 )
   W( 71 ) = JVS( 171 )
   W( 75 ) = JVS( 172 )
   W( 77 ) = JVS( 173 )
  JVS( 166) = W( 41 )
  JVS( 167) = W( 55 )
  JVS( 168) = W( 63 )
  JVS( 169) = W( 68 )
  JVS( 170) = W( 70 )
  JVS( 171) = W( 71 )
  JVS( 172) = W( 75 )
  JVS( 173) = W( 77 )
  IF ( ABS(  JVS( 174 )) < TINY(a) ) THEN
         IER = 42                                      
         RETURN
  END IF
   W( 42 ) = JVS( 174 )
   W( 71 ) = JVS( 175 )
   W( 75 ) = JVS( 176 )
   W( 79 ) = JVS( 177 )
  JVS( 174) = W( 42 )
  JVS( 175) = W( 71 )
  JVS( 176) = W( 75 )
  JVS( 177) = W( 79 )
  IF ( ABS(  JVS( 178 )) < TINY(a) ) THEN
         IER = 43                                      
         RETURN
  END IF
   W( 43 ) = JVS( 178 )
   W( 73 ) = JVS( 179 )
   W( 74 ) = JVS( 180 )
   W( 75 ) = JVS( 181 )
  JVS( 178) = W( 43 )
  JVS( 179) = W( 73 )
  JVS( 180) = W( 74 )
  JVS( 181) = W( 75 )
  IF ( ABS(  JVS( 184 )) < TINY(a) ) THEN
         IER = 44                                      
         RETURN
  END IF
   W( 10 ) = JVS( 182 )
   W( 31 ) = JVS( 183 )
   W( 44 ) = JVS( 184 )
   W( 75 ) = JVS( 185 )
   W( 76 ) = JVS( 186 )
   W( 79 ) = JVS( 187 )
  a = -W( 10 ) / JVS(           23  )
  W( 10 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 24 )
  a = -W( 31 ) / JVS(          108  )
  W( 31 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 109 )
  W( 75 ) = W( 75 ) + a*JVS( 110 )
  W( 79 ) = W( 79 ) + a*JVS( 111 )
  JVS( 182) = W( 10 )
  JVS( 183) = W( 31 )
  JVS( 184) = W( 44 )
  JVS( 185) = W( 75 )
  JVS( 186) = W( 76 )
  JVS( 187) = W( 79 )
  IF ( ABS(  JVS( 191 )) < TINY(a) ) THEN
         IER = 45                                      
         RETURN
  END IF
   W( 14 ) = JVS( 188 )
   W( 31 ) = JVS( 189 )
   W( 44 ) = JVS( 190 )
   W( 45 ) = JVS( 191 )
   W( 73 ) = JVS( 192 )
   W( 75 ) = JVS( 193 )
   W( 76 ) = JVS( 194 )
   W( 79 ) = JVS( 195 )
  a = -W( 14 ) / JVS(           35  )
  W( 14 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 36 )
  W( 75 ) = W( 75 ) + a*JVS( 37 )
  a = -W( 31 ) / JVS(          108  )
  W( 31 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 109 )
  W( 75 ) = W( 75 ) + a*JVS( 110 )
  W( 79 ) = W( 79 ) + a*JVS( 111 )
  a = -W( 44 ) / JVS(          184  )
  W( 44 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 185 )
  W( 76 ) = W( 76 ) + a*JVS( 186 )
  W( 79 ) = W( 79 ) + a*JVS( 187 )
  JVS( 188) = W( 14 )
  JVS( 189) = W( 31 )
  JVS( 190) = W( 44 )
  JVS( 191) = W( 45 )
  JVS( 192) = W( 73 )
  JVS( 193) = W( 75 )
  JVS( 194) = W( 76 )
  JVS( 195) = W( 79 )
  IF ( ABS(  JVS( 196 )) < TINY(a) ) THEN
         IER = 46                                      
         RETURN
  END IF
   W( 46 ) = JVS( 196 )
   W( 75 ) = JVS( 197 )
   W( 80 ) = JVS( 198 )
   W( 81 ) = JVS( 199 )
  JVS( 196) = W( 46 )
  JVS( 197) = W( 75 )
  JVS( 198) = W( 80 )
  JVS( 199) = W( 81 )
  IF ( ABS(  JVS( 203 )) < TINY(a) ) THEN
         IER = 47                                      
         RETURN
  END IF
   W( 3 ) = JVS( 200 )
   W( 19 ) = JVS( 201 )
   W( 40 ) = JVS( 202 )
   W( 47 ) = JVS( 203 )
   W( 61 ) = JVS( 204 )
   W( 72 ) = JVS( 205 )
   W( 75 ) = JVS( 206 )
   W( 81 ) = JVS( 207 )
  a = -W( 3 ) / JVS(            6  )
  W( 3 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 7 )
  a = -W( 19 ) / JVS(           53  )
  W( 19 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 54 )
  W( 47 ) = W( 47 ) + a*JVS( 55 )
  W( 72 ) = W( 72 ) + a*JVS( 56 )
  W( 75 ) = W( 75 ) + a*JVS( 57 )
  a = -W( 40 ) / JVS(          161  )
  W( 40 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 162 )
  W( 61 ) = W( 61 ) + a*JVS( 163 )
  W( 75 ) = W( 75 ) + a*JVS( 164 )
  W( 81 ) = W( 81 ) + a*JVS( 165 )
  JVS( 200) = W( 3 )
  JVS( 201) = W( 19 )
  JVS( 202) = W( 40 )
  JVS( 203) = W( 47 )
  JVS( 204) = W( 61 )
  JVS( 205) = W( 72 )
  JVS( 206) = W( 75 )
  JVS( 207) = W( 81 )
  IF ( ABS(  JVS( 208 )) < TINY(a) ) THEN
         IER = 48                                      
         RETURN
  END IF
   W( 48 ) = JVS( 208 )
   W( 54 ) = JVS( 209 )
   W( 75 ) = JVS( 210 )
   W( 79 ) = JVS( 211 )
  JVS( 208) = W( 48 )
  JVS( 209) = W( 54 )
  JVS( 210) = W( 75 )
  JVS( 211) = W( 79 )
  IF ( ABS(  JVS( 215 )) < TINY(a) ) THEN
         IER = 49                                      
         RETURN
  END IF
   W( 25 ) = JVS( 212 )
   W( 33 ) = JVS( 213 )
   W( 48 ) = JVS( 214 )
   W( 49 ) = JVS( 215 )
   W( 54 ) = JVS( 216 )
   W( 75 ) = JVS( 217 )
   W( 76 ) = JVS( 218 )
   W( 79 ) = JVS( 219 )
  a = -W( 25 ) / JVS(           78  )
  W( 25 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 79 )
  W( 75 ) = W( 75 ) + a*JVS( 80 )
  W( 79 ) = W( 79 ) + a*JVS( 81 )
  a = -W( 33 ) / JVS(          117  )
  W( 33 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 118 )
  W( 54 ) = W( 54 ) + a*JVS( 119 )
  W( 75 ) = W( 75 ) + a*JVS( 120 )
  W( 76 ) = W( 76 ) + a*JVS( 121 )
  a = -W( 48 ) / JVS(          208  )
  W( 48 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 209 )
  W( 75 ) = W( 75 ) + a*JVS( 210 )
  W( 79 ) = W( 79 ) + a*JVS( 211 )
  JVS( 212) = W( 25 )
  JVS( 213) = W( 33 )
  JVS( 214) = W( 48 )
  JVS( 215) = W( 49 )
  JVS( 216) = W( 54 )
  JVS( 217) = W( 75 )
  JVS( 218) = W( 76 )
  JVS( 219) = W( 79 )
  IF ( ABS(  JVS( 221 )) < TINY(a) ) THEN
         IER = 50                                      
         RETURN
  END IF
   W( 35 ) = JVS( 220 )
   W( 50 ) = JVS( 221 )
   W( 61 ) = JVS( 222 )
   W( 75 ) = JVS( 223 )
   W( 76 ) = JVS( 224 )
   W( 79 ) = JVS( 225 )
  a = -W( 35 ) / JVS(          132  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 133 )
  W( 75 ) = W( 75 ) + a*JVS( 134 )
  W( 79 ) = W( 79 ) + a*JVS( 135 )
  JVS( 220) = W( 35 )
  JVS( 221) = W( 50 )
  JVS( 222) = W( 61 )
  JVS( 223) = W( 75 )
  JVS( 224) = W( 76 )
  JVS( 225) = W( 79 )
  IF ( ABS(  JVS( 230 )) < TINY(a) ) THEN
         IER = 51                                      
         RETURN
  END IF
   W( 22 ) = JVS( 226 )
   W( 24 ) = JVS( 227 )
   W( 38 ) = JVS( 228 )
   W( 48 ) = JVS( 229 )
   W( 51 ) = JVS( 230 )
   W( 54 ) = JVS( 231 )
   W( 57 ) = JVS( 232 )
   W( 64 ) = JVS( 233 )
   W( 75 ) = JVS( 234 )
   W( 76 ) = JVS( 235 )
   W( 77 ) = JVS( 236 )
   W( 79 ) = JVS( 237 )
  a = -W( 22 ) / JVS(           66  )
  W( 22 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 67 )
  W( 75 ) = W( 75 ) + a*JVS( 68 )
  W( 79 ) = W( 79 ) + a*JVS( 69 )
  a = -W( 24 ) / JVS(           75  )
  W( 24 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 76 )
  W( 76 ) = W( 76 ) + a*JVS( 77 )
  a = -W( 38 ) / JVS(          148  )
  W( 38 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 149 )
  W( 75 ) = W( 75 ) + a*JVS( 150 )
  W( 79 ) = W( 79 ) + a*JVS( 151 )
  a = -W( 48 ) / JVS(          208  )
  W( 48 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 209 )
  W( 75 ) = W( 75 ) + a*JVS( 210 )
  W( 79 ) = W( 79 ) + a*JVS( 211 )
  JVS( 226) = W( 22 )
  JVS( 227) = W( 24 )
  JVS( 228) = W( 38 )
  JVS( 229) = W( 48 )
  JVS( 230) = W( 51 )
  JVS( 231) = W( 54 )
  JVS( 232) = W( 57 )
  JVS( 233) = W( 64 )
  JVS( 234) = W( 75 )
  JVS( 235) = W( 76 )
  JVS( 236) = W( 77 )
  JVS( 237) = W( 79 )
  IF ( ABS(  JVS( 239 )) < TINY(a) ) THEN
         IER = 52                                      
         RETURN
  END IF
   W( 47 ) = JVS( 238 )
   W( 52 ) = JVS( 239 )
   W( 61 ) = JVS( 240 )
   W( 72 ) = JVS( 241 )
   W( 73 ) = JVS( 242 )
   W( 75 ) = JVS( 243 )
   W( 79 ) = JVS( 244 )
   W( 81 ) = JVS( 245 )
  a = -W( 47 ) / JVS(          203  )
  W( 47 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 204 )
  W( 72 ) = W( 72 ) + a*JVS( 205 )
  W( 75 ) = W( 75 ) + a*JVS( 206 )
  W( 81 ) = W( 81 ) + a*JVS( 207 )
  JVS( 238) = W( 47 )
  JVS( 239) = W( 52 )
  JVS( 240) = W( 61 )
  JVS( 241) = W( 72 )
  JVS( 242) = W( 73 )
  JVS( 243) = W( 75 )
  JVS( 244) = W( 79 )
  JVS( 245) = W( 81 )
  IF ( ABS(  JVS( 248 )) < TINY(a) ) THEN
         IER = 53                                      
         RETURN
  END IF
   W( 12 ) = JVS( 246 )
   W( 30 ) = JVS( 247 )
   W( 53 ) = JVS( 248 )
   W( 68 ) = JVS( 249 )
   W( 70 ) = JVS( 250 )
   W( 75 ) = JVS( 251 )
   W( 76 ) = JVS( 252 )
   W( 77 ) = JVS( 253 )
   W( 78 ) = JVS( 254 )
   W( 80 ) = JVS( 255 )
   W( 81 ) = JVS( 256 )
  a = -W( 12 ) / JVS(           28  )
  W( 12 ) = -a
  W( 30 ) = W( 30 ) + a*JVS( 29 )
  W( 76 ) = W( 76 ) + a*JVS( 30 )
  a = -W( 30 ) / JVS(          104  )
  W( 30 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 105 )
  W( 76 ) = W( 76 ) + a*JVS( 106 )
  W( 81 ) = W( 81 ) + a*JVS( 107 )
  JVS( 246) = W( 12 )
  JVS( 247) = W( 30 )
  JVS( 248) = W( 53 )
  JVS( 249) = W( 68 )
  JVS( 250) = W( 70 )
  JVS( 251) = W( 75 )
  JVS( 252) = W( 76 )
  JVS( 253) = W( 77 )
  JVS( 254) = W( 78 )
  JVS( 255) = W( 80 )
  JVS( 256) = W( 81 )
  IF ( ABS(  JVS( 259 )) < TINY(a) ) THEN
         IER = 54                                      
         RETURN
  END IF
   W( 7 ) = JVS( 257 )
   W( 48 ) = JVS( 258 )
   W( 54 ) = JVS( 259 )
   W( 75 ) = JVS( 260 )
   W( 76 ) = JVS( 261 )
   W( 79 ) = JVS( 262 )
  a = -W( 7 ) / JVS(           14  )
  W( 7 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 15 )
  a = -W( 48 ) / JVS(          208  )
  W( 48 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 209 )
  W( 75 ) = W( 75 ) + a*JVS( 210 )
  W( 79 ) = W( 79 ) + a*JVS( 211 )
  JVS( 257) = W( 7 )
  JVS( 258) = W( 48 )
  JVS( 259) = W( 54 )
  JVS( 260) = W( 75 )
  JVS( 261) = W( 76 )
  JVS( 262) = W( 79 )
  IF ( ABS(  JVS( 268 )) < TINY(a) ) THEN
         IER = 55                                      
         RETURN
  END IF
   W( 4 ) = JVS( 263 )
   W( 21 ) = JVS( 264 )
   W( 33 ) = JVS( 265 )
   W( 48 ) = JVS( 266 )
   W( 54 ) = JVS( 267 )
   W( 55 ) = JVS( 268 )
   W( 75 ) = JVS( 269 )
   W( 76 ) = JVS( 270 )
   W( 77 ) = JVS( 271 )
   W( 79 ) = JVS( 272 )
  a = -W( 4 ) / JVS(            8  )
  W( 4 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 9 )
  a = -W( 21 ) / JVS(           62  )
  W( 21 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 63 )
  W( 75 ) = W( 75 ) + a*JVS( 64 )
  W( 79 ) = W( 79 ) + a*JVS( 65 )
  a = -W( 33 ) / JVS(          117  )
  W( 33 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 118 )
  W( 54 ) = W( 54 ) + a*JVS( 119 )
  W( 75 ) = W( 75 ) + a*JVS( 120 )
  W( 76 ) = W( 76 ) + a*JVS( 121 )
  a = -W( 48 ) / JVS(          208  )
  W( 48 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 209 )
  W( 75 ) = W( 75 ) + a*JVS( 210 )
  W( 79 ) = W( 79 ) + a*JVS( 211 )
  a = -W( 54 ) / JVS(          259  )
  W( 54 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 260 )
  W( 76 ) = W( 76 ) + a*JVS( 261 )
  W( 79 ) = W( 79 ) + a*JVS( 262 )
  JVS( 263) = W( 4 )
  JVS( 264) = W( 21 )
  JVS( 265) = W( 33 )
  JVS( 266) = W( 48 )
  JVS( 267) = W( 54 )
  JVS( 268) = W( 55 )
  JVS( 269) = W( 75 )
  JVS( 270) = W( 76 )
  JVS( 271) = W( 77 )
  JVS( 272) = W( 79 )
  IF ( ABS(  JVS( 273 )) < TINY(a) ) THEN
         IER = 56                                      
         RETURN
  END IF
   W( 56 ) = JVS( 273 )
   W( 59 ) = JVS( 274 )
   W( 76 ) = JVS( 275 )
   W( 79 ) = JVS( 276 )
   W( 80 ) = JVS( 277 )
  JVS( 273) = W( 56 )
  JVS( 274) = W( 59 )
  JVS( 275) = W( 76 )
  JVS( 276) = W( 79 )
  JVS( 277) = W( 80 )
  IF ( ABS(  JVS( 280 )) < TINY(a) ) THEN
         IER = 57                                      
         RETURN
  END IF
   W( 5 ) = JVS( 278 )
   W( 22 ) = JVS( 279 )
   W( 57 ) = JVS( 280 )
   W( 75 ) = JVS( 281 )
   W( 76 ) = JVS( 282 )
   W( 77 ) = JVS( 283 )
   W( 79 ) = JVS( 284 )
  a = -W( 5 ) / JVS(           10  )
  W( 5 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 11 )
  a = -W( 22 ) / JVS(           66  )
  W( 22 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 67 )
  W( 75 ) = W( 75 ) + a*JVS( 68 )
  W( 79 ) = W( 79 ) + a*JVS( 69 )
  JVS( 278) = W( 5 )
  JVS( 279) = W( 22 )
  JVS( 280) = W( 57 )
  JVS( 281) = W( 75 )
  JVS( 282) = W( 76 )
  JVS( 283) = W( 77 )
  JVS( 284) = W( 79 )
  IF ( ABS(  JVS( 290 )) < TINY(a) ) THEN
         IER = 58                                      
         RETURN
  END IF
   W( 27 ) = JVS( 285 )
   W( 37 ) = JVS( 286 )
   W( 44 ) = JVS( 287 )
   W( 45 ) = JVS( 288 )
   W( 53 ) = JVS( 289 )
   W( 58 ) = JVS( 290 )
   W( 59 ) = JVS( 291 )
   W( 60 ) = JVS( 292 )
   W( 61 ) = JVS( 293 )
   W( 62 ) = JVS( 294 )
   W( 66 ) = JVS( 295 )
   W( 67 ) = JVS( 296 )
   W( 68 ) = JVS( 297 )
   W( 69 ) = JVS( 298 )
   W( 70 ) = JVS( 299 )
   W( 72 ) = JVS( 300 )
   W( 73 ) = JVS( 301 )
   W( 75 ) = JVS( 302 )
   W( 76 ) = JVS( 303 )
   W( 77 ) = JVS( 304 )
   W( 78 ) = JVS( 305 )
   W( 79 ) = JVS( 306 )
   W( 80 ) = JVS( 307 )
   W( 81 ) = JVS( 308 )
  a = -W( 27 ) / JVS(           86  )
  W( 27 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 87 )
  W( 81 ) = W( 81 ) + a*JVS( 88 )
  a = -W( 37 ) / JVS(          141  )
  W( 37 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 142 )
  W( 45 ) = W( 45 ) + a*JVS( 143 )
  W( 53 ) = W( 53 ) + a*JVS( 144 )
  W( 75 ) = W( 75 ) + a*JVS( 145 )
  W( 76 ) = W( 76 ) + a*JVS( 146 )
  W( 79 ) = W( 79 ) + a*JVS( 147 )
  a = -W( 44 ) / JVS(          184  )
  W( 44 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 185 )
  W( 76 ) = W( 76 ) + a*JVS( 186 )
  W( 79 ) = W( 79 ) + a*JVS( 187 )
  a = -W( 45 ) / JVS(          191  )
  W( 45 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 192 )
  W( 75 ) = W( 75 ) + a*JVS( 193 )
  W( 76 ) = W( 76 ) + a*JVS( 194 )
  W( 79 ) = W( 79 ) + a*JVS( 195 )
  a = -W( 53 ) / JVS(          248  )
  W( 53 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 249 )
  W( 70 ) = W( 70 ) + a*JVS( 250 )
  W( 75 ) = W( 75 ) + a*JVS( 251 )
  W( 76 ) = W( 76 ) + a*JVS( 252 )
  W( 77 ) = W( 77 ) + a*JVS( 253 )
  W( 78 ) = W( 78 ) + a*JVS( 254 )
  W( 80 ) = W( 80 ) + a*JVS( 255 )
  W( 81 ) = W( 81 ) + a*JVS( 256 )
  JVS( 285) = W( 27 )
  JVS( 286) = W( 37 )
  JVS( 287) = W( 44 )
  JVS( 288) = W( 45 )
  JVS( 289) = W( 53 )
  JVS( 290) = W( 58 )
  JVS( 291) = W( 59 )
  JVS( 292) = W( 60 )
  JVS( 293) = W( 61 )
  JVS( 294) = W( 62 )
  JVS( 295) = W( 66 )
  JVS( 296) = W( 67 )
  JVS( 297) = W( 68 )
  JVS( 298) = W( 69 )
  JVS( 299) = W( 70 )
  JVS( 300) = W( 72 )
  JVS( 301) = W( 73 )
  JVS( 302) = W( 75 )
  JVS( 303) = W( 76 )
  JVS( 304) = W( 77 )
  JVS( 305) = W( 78 )
  JVS( 306) = W( 79 )
  JVS( 307) = W( 80 )
  JVS( 308) = W( 81 )
  IF ( ABS(  JVS( 309 )) < TINY(a) ) THEN
         IER = 59                                      
         RETURN
  END IF
   W( 59 ) = JVS( 309 )
   W( 75 ) = JVS( 310 )
   W( 80 ) = JVS( 311 )
   W( 81 ) = JVS( 312 )
  JVS( 309) = W( 59 )
  JVS( 310) = W( 75 )
  JVS( 311) = W( 80 )
  JVS( 312) = W( 81 )
  IF ( ABS(  JVS( 315 )) < TINY(a) ) THEN
         IER = 60                                      
         RETURN
  END IF
   W( 56 ) = JVS( 313 )
   W( 59 ) = JVS( 314 )
   W( 60 ) = JVS( 315 )
   W( 70 ) = JVS( 316 )
   W( 71 ) = JVS( 317 )
   W( 75 ) = JVS( 318 )
   W( 76 ) = JVS( 319 )
   W( 79 ) = JVS( 320 )
   W( 80 ) = JVS( 321 )
   W( 81 ) = JVS( 322 )
  a = -W( 56 ) / JVS(          273  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 274 )
  W( 76 ) = W( 76 ) + a*JVS( 275 )
  W( 79 ) = W( 79 ) + a*JVS( 276 )
  W( 80 ) = W( 80 ) + a*JVS( 277 )
  a = -W( 59 ) / JVS(          309  )
  W( 59 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 310 )
  W( 80 ) = W( 80 ) + a*JVS( 311 )
  W( 81 ) = W( 81 ) + a*JVS( 312 )
  JVS( 313) = W( 56 )
  JVS( 314) = W( 59 )
  JVS( 315) = W( 60 )
  JVS( 316) = W( 70 )
  JVS( 317) = W( 71 )
  JVS( 318) = W( 75 )
  JVS( 319) = W( 76 )
  JVS( 320) = W( 79 )
  JVS( 321) = W( 80 )
  JVS( 322) = W( 81 )
  IF ( ABS(  JVS( 324 )) < TINY(a) ) THEN
         IER = 61                                      
         RETURN
  END IF
   W( 59 ) = JVS( 323 )
   W( 61 ) = JVS( 324 )
   W( 69 ) = JVS( 325 )
   W( 75 ) = JVS( 326 )
   W( 80 ) = JVS( 327 )
   W( 81 ) = JVS( 328 )
  a = -W( 59 ) / JVS(          309  )
  W( 59 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 310 )
  W( 80 ) = W( 80 ) + a*JVS( 311 )
  W( 81 ) = W( 81 ) + a*JVS( 312 )
  JVS( 323) = W( 59 )
  JVS( 324) = W( 61 )
  JVS( 325) = W( 69 )
  JVS( 326) = W( 75 )
  JVS( 327) = W( 80 )
  JVS( 328) = W( 81 )
  IF ( ABS(  JVS( 341 )) < TINY(a) ) THEN
         IER = 62                                      
         RETURN
  END IF
   W( 17 ) = JVS( 329 )
   W( 21 ) = JVS( 330 )
   W( 24 ) = JVS( 331 )
   W( 25 ) = JVS( 332 )
   W( 35 ) = JVS( 333 )
   W( 48 ) = JVS( 334 )
   W( 49 ) = JVS( 335 )
   W( 50 ) = JVS( 336 )
   W( 54 ) = JVS( 337 )
   W( 55 ) = JVS( 338 )
   W( 57 ) = JVS( 339 )
   W( 61 ) = JVS( 340 )
   W( 62 ) = JVS( 341 )
   W( 69 ) = JVS( 342 )
   W( 75 ) = JVS( 343 )
   W( 76 ) = JVS( 344 )
   W( 77 ) = JVS( 345 )
   W( 79 ) = JVS( 346 )
   W( 80 ) = JVS( 347 )
   W( 81 ) = JVS( 348 )
  a = -W( 17 ) / JVS(           45  )
  W( 17 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 46 )
  W( 75 ) = W( 75 ) + a*JVS( 47 )
  W( 77 ) = W( 77 ) + a*JVS( 48 )
  a = -W( 21 ) / JVS(           62  )
  W( 21 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 63 )
  W( 75 ) = W( 75 ) + a*JVS( 64 )
  W( 79 ) = W( 79 ) + a*JVS( 65 )
  a = -W( 24 ) / JVS(           75  )
  W( 24 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 76 )
  W( 76 ) = W( 76 ) + a*JVS( 77 )
  a = -W( 25 ) / JVS(           78  )
  W( 25 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 79 )
  W( 75 ) = W( 75 ) + a*JVS( 80 )
  W( 79 ) = W( 79 ) + a*JVS( 81 )
  a = -W( 35 ) / JVS(          132  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 133 )
  W( 75 ) = W( 75 ) + a*JVS( 134 )
  W( 79 ) = W( 79 ) + a*JVS( 135 )
  a = -W( 48 ) / JVS(          208  )
  W( 48 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 209 )
  W( 75 ) = W( 75 ) + a*JVS( 210 )
  W( 79 ) = W( 79 ) + a*JVS( 211 )
  a = -W( 49 ) / JVS(          215  )
  W( 49 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 216 )
  W( 75 ) = W( 75 ) + a*JVS( 217 )
  W( 76 ) = W( 76 ) + a*JVS( 218 )
  W( 79 ) = W( 79 ) + a*JVS( 219 )
  a = -W( 50 ) / JVS(          221  )
  W( 50 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 222 )
  W( 75 ) = W( 75 ) + a*JVS( 223 )
  W( 76 ) = W( 76 ) + a*JVS( 224 )
  W( 79 ) = W( 79 ) + a*JVS( 225 )
  a = -W( 54 ) / JVS(          259  )
  W( 54 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 260 )
  W( 76 ) = W( 76 ) + a*JVS( 261 )
  W( 79 ) = W( 79 ) + a*JVS( 262 )
  a = -W( 55 ) / JVS(          268  )
  W( 55 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 269 )
  W( 76 ) = W( 76 ) + a*JVS( 270 )
  W( 77 ) = W( 77 ) + a*JVS( 271 )
  W( 79 ) = W( 79 ) + a*JVS( 272 )
  a = -W( 57 ) / JVS(          280  )
  W( 57 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 281 )
  W( 76 ) = W( 76 ) + a*JVS( 282 )
  W( 77 ) = W( 77 ) + a*JVS( 283 )
  W( 79 ) = W( 79 ) + a*JVS( 284 )
  a = -W( 61 ) / JVS(          324  )
  W( 61 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 325 )
  W( 75 ) = W( 75 ) + a*JVS( 326 )
  W( 80 ) = W( 80 ) + a*JVS( 327 )
  W( 81 ) = W( 81 ) + a*JVS( 328 )
  JVS( 329) = W( 17 )
  JVS( 330) = W( 21 )
  JVS( 331) = W( 24 )
  JVS( 332) = W( 25 )
  JVS( 333) = W( 35 )
  JVS( 334) = W( 48 )
  JVS( 335) = W( 49 )
  JVS( 336) = W( 50 )
  JVS( 337) = W( 54 )
  JVS( 338) = W( 55 )
  JVS( 339) = W( 57 )
  JVS( 340) = W( 61 )
  JVS( 341) = W( 62 )
  JVS( 342) = W( 69 )
  JVS( 343) = W( 75 )
  JVS( 344) = W( 76 )
  JVS( 345) = W( 77 )
  JVS( 346) = W( 79 )
  JVS( 347) = W( 80 )
  JVS( 348) = W( 81 )
  IF ( ABS(  JVS( 353 )) < TINY(a) ) THEN
         IER = 63                                      
         RETURN
  END IF
   W( 23 ) = JVS( 349 )
   W( 51 ) = JVS( 350 )
   W( 54 ) = JVS( 351 )
   W( 57 ) = JVS( 352 )
   W( 63 ) = JVS( 353 )
   W( 64 ) = JVS( 354 )
   W( 75 ) = JVS( 355 )
   W( 76 ) = JVS( 356 )
   W( 77 ) = JVS( 357 )
   W( 79 ) = JVS( 358 )
  a = -W( 23 ) / JVS(           70  )
  W( 23 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 71 )
  W( 75 ) = W( 75 ) + a*JVS( 72 )
  W( 79 ) = W( 79 ) + a*JVS( 73 )
  a = -W( 51 ) / JVS(          230  )
  W( 51 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 231 )
  W( 57 ) = W( 57 ) + a*JVS( 232 )
  W( 64 ) = W( 64 ) + a*JVS( 233 )
  W( 75 ) = W( 75 ) + a*JVS( 234 )
  W( 76 ) = W( 76 ) + a*JVS( 235 )
  W( 77 ) = W( 77 ) + a*JVS( 236 )
  W( 79 ) = W( 79 ) + a*JVS( 237 )
  a = -W( 54 ) / JVS(          259  )
  W( 54 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 260 )
  W( 76 ) = W( 76 ) + a*JVS( 261 )
  W( 79 ) = W( 79 ) + a*JVS( 262 )
  a = -W( 57 ) / JVS(          280  )
  W( 57 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 281 )
  W( 76 ) = W( 76 ) + a*JVS( 282 )
  W( 77 ) = W( 77 ) + a*JVS( 283 )
  W( 79 ) = W( 79 ) + a*JVS( 284 )
  JVS( 349) = W( 23 )
  JVS( 350) = W( 51 )
  JVS( 351) = W( 54 )
  JVS( 352) = W( 57 )
  JVS( 353) = W( 63 )
  JVS( 354) = W( 64 )
  JVS( 355) = W( 75 )
  JVS( 356) = W( 76 )
  JVS( 357) = W( 77 )
  JVS( 358) = W( 79 )
  IF ( ABS(  JVS( 361 )) < TINY(a) ) THEN
         IER = 64                                      
         RETURN
  END IF
   W( 38 ) = JVS( 359 )
   W( 46 ) = JVS( 360 )
   W( 64 ) = JVS( 361 )
   W( 75 ) = JVS( 362 )
   W( 76 ) = JVS( 363 )
   W( 79 ) = JVS( 364 )
   W( 80 ) = JVS( 365 )
   W( 81 ) = JVS( 366 )
  a = -W( 38 ) / JVS(          148  )
  W( 38 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 149 )
  W( 75 ) = W( 75 ) + a*JVS( 150 )
  W( 79 ) = W( 79 ) + a*JVS( 151 )
  a = -W( 46 ) / JVS(          196  )
  W( 46 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 197 )
  W( 80 ) = W( 80 ) + a*JVS( 198 )
  W( 81 ) = W( 81 ) + a*JVS( 199 )
  JVS( 359) = W( 38 )
  JVS( 360) = W( 46 )
  JVS( 361) = W( 64 )
  JVS( 362) = W( 75 )
  JVS( 363) = W( 76 )
  JVS( 364) = W( 79 )
  JVS( 365) = W( 80 )
  JVS( 366) = W( 81 )
  IF ( ABS(  JVS( 373 )) < TINY(a) ) THEN
         IER = 65                                      
         RETURN
  END IF
   W( 35 ) = JVS( 367 )
   W( 43 ) = JVS( 368 )
   W( 50 ) = JVS( 369 )
   W( 61 ) = JVS( 370 )
   W( 63 ) = JVS( 371 )
   W( 64 ) = JVS( 372 )
   W( 65 ) = JVS( 373 )
   W( 68 ) = JVS( 374 )
   W( 69 ) = JVS( 375 )
   W( 70 ) = JVS( 376 )
   W( 73 ) = JVS( 377 )
   W( 74 ) = JVS( 378 )
   W( 75 ) = JVS( 379 )
   W( 76 ) = JVS( 380 )
   W( 77 ) = JVS( 381 )
   W( 78 ) = JVS( 382 )
   W( 79 ) = JVS( 383 )
   W( 80 ) = JVS( 384 )
   W( 81 ) = JVS( 385 )
  a = -W( 35 ) / JVS(          132  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 133 )
  W( 75 ) = W( 75 ) + a*JVS( 134 )
  W( 79 ) = W( 79 ) + a*JVS( 135 )
  a = -W( 43 ) / JVS(          178  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 179 )
  W( 74 ) = W( 74 ) + a*JVS( 180 )
  W( 75 ) = W( 75 ) + a*JVS( 181 )
  a = -W( 50 ) / JVS(          221  )
  W( 50 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 222 )
  W( 75 ) = W( 75 ) + a*JVS( 223 )
  W( 76 ) = W( 76 ) + a*JVS( 224 )
  W( 79 ) = W( 79 ) + a*JVS( 225 )
  a = -W( 61 ) / JVS(          324  )
  W( 61 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 325 )
  W( 75 ) = W( 75 ) + a*JVS( 326 )
  W( 80 ) = W( 80 ) + a*JVS( 327 )
  W( 81 ) = W( 81 ) + a*JVS( 328 )
  a = -W( 63 ) / JVS(          353  )
  W( 63 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  W( 76 ) = W( 76 ) + a*JVS( 356 )
  W( 77 ) = W( 77 ) + a*JVS( 357 )
  W( 79 ) = W( 79 ) + a*JVS( 358 )
  a = -W( 64 ) / JVS(          361  )
  W( 64 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 362 )
  W( 76 ) = W( 76 ) + a*JVS( 363 )
  W( 79 ) = W( 79 ) + a*JVS( 364 )
  W( 80 ) = W( 80 ) + a*JVS( 365 )
  W( 81 ) = W( 81 ) + a*JVS( 366 )
  JVS( 367) = W( 35 )
  JVS( 368) = W( 43 )
  JVS( 369) = W( 50 )
  JVS( 370) = W( 61 )
  JVS( 371) = W( 63 )
  JVS( 372) = W( 64 )
  JVS( 373) = W( 65 )
  JVS( 374) = W( 68 )
  JVS( 375) = W( 69 )
  JVS( 376) = W( 70 )
  JVS( 377) = W( 73 )
  JVS( 378) = W( 74 )
  JVS( 379) = W( 75 )
  JVS( 380) = W( 76 )
  JVS( 381) = W( 77 )
  JVS( 382) = W( 78 )
  JVS( 383) = W( 79 )
  JVS( 384) = W( 80 )
  JVS( 385) = W( 81 )
  IF ( ABS(  JVS( 392 )) < TINY(a) ) THEN
         IER = 66                                      
         RETURN
  END IF
   W( 38 ) = JVS( 386 )
   W( 42 ) = JVS( 387 )
   W( 46 ) = JVS( 388 )
   W( 56 ) = JVS( 389 )
   W( 59 ) = JVS( 390 )
   W( 64 ) = JVS( 391 )
   W( 66 ) = JVS( 392 )
   W( 71 ) = JVS( 393 )
   W( 75 ) = JVS( 394 )
   W( 76 ) = JVS( 395 )
   W( 77 ) = JVS( 396 )
   W( 78 ) = JVS( 397 )
   W( 79 ) = JVS( 398 )
   W( 80 ) = JVS( 399 )
   W( 81 ) = JVS( 400 )
  a = -W( 38 ) / JVS(          148  )
  W( 38 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 149 )
  W( 75 ) = W( 75 ) + a*JVS( 150 )
  W( 79 ) = W( 79 ) + a*JVS( 151 )
  a = -W( 42 ) / JVS(          174  )
  W( 42 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 175 )
  W( 75 ) = W( 75 ) + a*JVS( 176 )
  W( 79 ) = W( 79 ) + a*JVS( 177 )
  a = -W( 46 ) / JVS(          196  )
  W( 46 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 197 )
  W( 80 ) = W( 80 ) + a*JVS( 198 )
  W( 81 ) = W( 81 ) + a*JVS( 199 )
  a = -W( 56 ) / JVS(          273  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 274 )
  W( 76 ) = W( 76 ) + a*JVS( 275 )
  W( 79 ) = W( 79 ) + a*JVS( 276 )
  W( 80 ) = W( 80 ) + a*JVS( 277 )
  a = -W( 59 ) / JVS(          309  )
  W( 59 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 310 )
  W( 80 ) = W( 80 ) + a*JVS( 311 )
  W( 81 ) = W( 81 ) + a*JVS( 312 )
  a = -W( 64 ) / JVS(          361  )
  W( 64 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 362 )
  W( 76 ) = W( 76 ) + a*JVS( 363 )
  W( 79 ) = W( 79 ) + a*JVS( 364 )
  W( 80 ) = W( 80 ) + a*JVS( 365 )
  W( 81 ) = W( 81 ) + a*JVS( 366 )
  JVS( 386) = W( 38 )
  JVS( 387) = W( 42 )
  JVS( 388) = W( 46 )
  JVS( 389) = W( 56 )
  JVS( 390) = W( 59 )
  JVS( 391) = W( 64 )
  JVS( 392) = W( 66 )
  JVS( 393) = W( 71 )
  JVS( 394) = W( 75 )
  JVS( 395) = W( 76 )
  JVS( 396) = W( 77 )
  JVS( 397) = W( 78 )
  JVS( 398) = W( 79 )
  JVS( 399) = W( 80 )
  JVS( 400) = W( 81 )
  IF ( ABS(  JVS( 411 )) < TINY(a) ) THEN
         IER = 67                                      
         RETURN
  END IF
   W( 29 ) = JVS( 401 )
   W( 31 ) = JVS( 402 )
   W( 44 ) = JVS( 403 )
   W( 45 ) = JVS( 404 )
   W( 54 ) = JVS( 405 )
   W( 61 ) = JVS( 406 )
   W( 63 ) = JVS( 407 )
   W( 64 ) = JVS( 408 )
   W( 65 ) = JVS( 409 )
   W( 66 ) = JVS( 410 )
   W( 67 ) = JVS( 411 )
   W( 68 ) = JVS( 412 )
   W( 69 ) = JVS( 413 )
   W( 70 ) = JVS( 414 )
   W( 71 ) = JVS( 415 )
   W( 73 ) = JVS( 416 )
   W( 74 ) = JVS( 417 )
   W( 75 ) = JVS( 418 )
   W( 76 ) = JVS( 419 )
   W( 77 ) = JVS( 420 )
   W( 78 ) = JVS( 421 )
   W( 79 ) = JVS( 422 )
   W( 80 ) = JVS( 423 )
   W( 81 ) = JVS( 424 )
  a = -W( 29 ) / JVS(           97  )
  W( 29 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 98 )
  W( 61 ) = W( 61 ) + a*JVS( 99 )
  W( 75 ) = W( 75 ) + a*JVS( 100 )
  W( 76 ) = W( 76 ) + a*JVS( 101 )
  W( 80 ) = W( 80 ) + a*JVS( 102 )
  a = -W( 31 ) / JVS(          108  )
  W( 31 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 109 )
  W( 75 ) = W( 75 ) + a*JVS( 110 )
  W( 79 ) = W( 79 ) + a*JVS( 111 )
  a = -W( 44 ) / JVS(          184  )
  W( 44 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 185 )
  W( 76 ) = W( 76 ) + a*JVS( 186 )
  W( 79 ) = W( 79 ) + a*JVS( 187 )
  a = -W( 45 ) / JVS(          191  )
  W( 45 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 192 )
  W( 75 ) = W( 75 ) + a*JVS( 193 )
  W( 76 ) = W( 76 ) + a*JVS( 194 )
  W( 79 ) = W( 79 ) + a*JVS( 195 )
  a = -W( 54 ) / JVS(          259  )
  W( 54 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 260 )
  W( 76 ) = W( 76 ) + a*JVS( 261 )
  W( 79 ) = W( 79 ) + a*JVS( 262 )
  a = -W( 61 ) / JVS(          324  )
  W( 61 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 325 )
  W( 75 ) = W( 75 ) + a*JVS( 326 )
  W( 80 ) = W( 80 ) + a*JVS( 327 )
  W( 81 ) = W( 81 ) + a*JVS( 328 )
  a = -W( 63 ) / JVS(          353  )
  W( 63 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  W( 76 ) = W( 76 ) + a*JVS( 356 )
  W( 77 ) = W( 77 ) + a*JVS( 357 )
  W( 79 ) = W( 79 ) + a*JVS( 358 )
  a = -W( 64 ) / JVS(          361  )
  W( 64 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 362 )
  W( 76 ) = W( 76 ) + a*JVS( 363 )
  W( 79 ) = W( 79 ) + a*JVS( 364 )
  W( 80 ) = W( 80 ) + a*JVS( 365 )
  W( 81 ) = W( 81 ) + a*JVS( 366 )
  a = -W( 65 ) / JVS(          373  )
  W( 65 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 374 )
  W( 69 ) = W( 69 ) + a*JVS( 375 )
  W( 70 ) = W( 70 ) + a*JVS( 376 )
  W( 73 ) = W( 73 ) + a*JVS( 377 )
  W( 74 ) = W( 74 ) + a*JVS( 378 )
  W( 75 ) = W( 75 ) + a*JVS( 379 )
  W( 76 ) = W( 76 ) + a*JVS( 380 )
  W( 77 ) = W( 77 ) + a*JVS( 381 )
  W( 78 ) = W( 78 ) + a*JVS( 382 )
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 80 ) = W( 80 ) + a*JVS( 384 )
  W( 81 ) = W( 81 ) + a*JVS( 385 )
  a = -W( 66 ) / JVS(          392  )
  W( 66 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 393 )
  W( 75 ) = W( 75 ) + a*JVS( 394 )
  W( 76 ) = W( 76 ) + a*JVS( 395 )
  W( 77 ) = W( 77 ) + a*JVS( 396 )
  W( 78 ) = W( 78 ) + a*JVS( 397 )
  W( 79 ) = W( 79 ) + a*JVS( 398 )
  W( 80 ) = W( 80 ) + a*JVS( 399 )
  W( 81 ) = W( 81 ) + a*JVS( 400 )
  JVS( 401) = W( 29 )
  JVS( 402) = W( 31 )
  JVS( 403) = W( 44 )
  JVS( 404) = W( 45 )
  JVS( 405) = W( 54 )
  JVS( 406) = W( 61 )
  JVS( 407) = W( 63 )
  JVS( 408) = W( 64 )
  JVS( 409) = W( 65 )
  JVS( 410) = W( 66 )
  JVS( 411) = W( 67 )
  JVS( 412) = W( 68 )
  JVS( 413) = W( 69 )
  JVS( 414) = W( 70 )
  JVS( 415) = W( 71 )
  JVS( 416) = W( 73 )
  JVS( 417) = W( 74 )
  JVS( 418) = W( 75 )
  JVS( 419) = W( 76 )
  JVS( 420) = W( 77 )
  JVS( 421) = W( 78 )
  JVS( 422) = W( 79 )
  JVS( 423) = W( 80 )
  JVS( 424) = W( 81 )
  IF ( ABS(  JVS( 429 )) < TINY(a) ) THEN
         IER = 68                                      
         RETURN
  END IF
   W( 15 ) = JVS( 425 )
   W( 28 ) = JVS( 426 )
   W( 42 ) = JVS( 427 )
   W( 60 ) = JVS( 428 )
   W( 68 ) = JVS( 429 )
   W( 70 ) = JVS( 430 )
   W( 71 ) = JVS( 431 )
   W( 75 ) = JVS( 432 )
   W( 76 ) = JVS( 433 )
   W( 77 ) = JVS( 434 )
   W( 78 ) = JVS( 435 )
   W( 79 ) = JVS( 436 )
   W( 80 ) = JVS( 437 )
   W( 81 ) = JVS( 438 )
  a = -W( 15 ) / JVS(           38  )
  W( 15 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 39 )
  W( 75 ) = W( 75 ) + a*JVS( 40 )
  W( 79 ) = W( 79 ) + a*JVS( 41 )
  a = -W( 28 ) / JVS(           89  )
  W( 28 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 90 )
  W( 71 ) = W( 71 ) + a*JVS( 91 )
  W( 75 ) = W( 75 ) + a*JVS( 92 )
  W( 76 ) = W( 76 ) + a*JVS( 93 )
  W( 77 ) = W( 77 ) + a*JVS( 94 )
  W( 78 ) = W( 78 ) + a*JVS( 95 )
  W( 80 ) = W( 80 ) + a*JVS( 96 )
  a = -W( 42 ) / JVS(          174  )
  W( 42 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 175 )
  W( 75 ) = W( 75 ) + a*JVS( 176 )
  W( 79 ) = W( 79 ) + a*JVS( 177 )
  a = -W( 60 ) / JVS(          315  )
  W( 60 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 316 )
  W( 71 ) = W( 71 ) + a*JVS( 317 )
  W( 75 ) = W( 75 ) + a*JVS( 318 )
  W( 76 ) = W( 76 ) + a*JVS( 319 )
  W( 79 ) = W( 79 ) + a*JVS( 320 )
  W( 80 ) = W( 80 ) + a*JVS( 321 )
  W( 81 ) = W( 81 ) + a*JVS( 322 )
  JVS( 425) = W( 15 )
  JVS( 426) = W( 28 )
  JVS( 427) = W( 42 )
  JVS( 428) = W( 60 )
  JVS( 429) = W( 68 )
  JVS( 430) = W( 70 )
  JVS( 431) = W( 71 )
  JVS( 432) = W( 75 )
  JVS( 433) = W( 76 )
  JVS( 434) = W( 77 )
  JVS( 435) = W( 78 )
  JVS( 436) = W( 79 )
  JVS( 437) = W( 80 )
  JVS( 438) = W( 81 )
  IF ( ABS(  JVS( 445 )) < TINY(a) ) THEN
         IER = 69                                      
         RETURN
  END IF
   W( 38 ) = JVS( 439 )
   W( 42 ) = JVS( 440 )
   W( 46 ) = JVS( 441 )
   W( 56 ) = JVS( 442 )
   W( 59 ) = JVS( 443 )
   W( 64 ) = JVS( 444 )
   W( 69 ) = JVS( 445 )
   W( 71 ) = JVS( 446 )
   W( 75 ) = JVS( 447 )
   W( 76 ) = JVS( 448 )
   W( 77 ) = JVS( 449 )
   W( 78 ) = JVS( 450 )
   W( 79 ) = JVS( 451 )
   W( 80 ) = JVS( 452 )
   W( 81 ) = JVS( 453 )
  a = -W( 38 ) / JVS(          148  )
  W( 38 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 149 )
  W( 75 ) = W( 75 ) + a*JVS( 150 )
  W( 79 ) = W( 79 ) + a*JVS( 151 )
  a = -W( 42 ) / JVS(          174  )
  W( 42 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 175 )
  W( 75 ) = W( 75 ) + a*JVS( 176 )
  W( 79 ) = W( 79 ) + a*JVS( 177 )
  a = -W( 46 ) / JVS(          196  )
  W( 46 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 197 )
  W( 80 ) = W( 80 ) + a*JVS( 198 )
  W( 81 ) = W( 81 ) + a*JVS( 199 )
  a = -W( 56 ) / JVS(          273  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 274 )
  W( 76 ) = W( 76 ) + a*JVS( 275 )
  W( 79 ) = W( 79 ) + a*JVS( 276 )
  W( 80 ) = W( 80 ) + a*JVS( 277 )
  a = -W( 59 ) / JVS(          309  )
  W( 59 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 310 )
  W( 80 ) = W( 80 ) + a*JVS( 311 )
  W( 81 ) = W( 81 ) + a*JVS( 312 )
  a = -W( 64 ) / JVS(          361  )
  W( 64 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 362 )
  W( 76 ) = W( 76 ) + a*JVS( 363 )
  W( 79 ) = W( 79 ) + a*JVS( 364 )
  W( 80 ) = W( 80 ) + a*JVS( 365 )
  W( 81 ) = W( 81 ) + a*JVS( 366 )
  JVS( 439) = W( 38 )
  JVS( 440) = W( 42 )
  JVS( 441) = W( 46 )
  JVS( 442) = W( 56 )
  JVS( 443) = W( 59 )
  JVS( 444) = W( 64 )
  JVS( 445) = W( 69 )
  JVS( 446) = W( 71 )
  JVS( 447) = W( 75 )
  JVS( 448) = W( 76 )
  JVS( 449) = W( 77 )
  JVS( 450) = W( 78 )
  JVS( 451) = W( 79 )
  JVS( 452) = W( 80 )
  JVS( 453) = W( 81 )
  IF ( ABS(  JVS( 457 )) < TINY(a) ) THEN
         IER = 70                                      
         RETURN
  END IF
   W( 20 ) = JVS( 454 )
   W( 66 ) = JVS( 455 )
   W( 69 ) = JVS( 456 )
   W( 70 ) = JVS( 457 )
   W( 71 ) = JVS( 458 )
   W( 75 ) = JVS( 459 )
   W( 76 ) = JVS( 460 )
   W( 77 ) = JVS( 461 )
   W( 78 ) = JVS( 462 )
   W( 79 ) = JVS( 463 )
   W( 80 ) = JVS( 464 )
   W( 81 ) = JVS( 465 )
  a = -W( 20 ) / JVS(           58  )
  W( 20 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 59 )
  W( 75 ) = W( 75 ) + a*JVS( 60 )
  W( 79 ) = W( 79 ) + a*JVS( 61 )
  a = -W( 66 ) / JVS(          392  )
  W( 66 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 393 )
  W( 75 ) = W( 75 ) + a*JVS( 394 )
  W( 76 ) = W( 76 ) + a*JVS( 395 )
  W( 77 ) = W( 77 ) + a*JVS( 396 )
  W( 78 ) = W( 78 ) + a*JVS( 397 )
  W( 79 ) = W( 79 ) + a*JVS( 398 )
  W( 80 ) = W( 80 ) + a*JVS( 399 )
  W( 81 ) = W( 81 ) + a*JVS( 400 )
  a = -W( 69 ) / JVS(          445  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 446 )
  W( 75 ) = W( 75 ) + a*JVS( 447 )
  W( 76 ) = W( 76 ) + a*JVS( 448 )
  W( 77 ) = W( 77 ) + a*JVS( 449 )
  W( 78 ) = W( 78 ) + a*JVS( 450 )
  W( 79 ) = W( 79 ) + a*JVS( 451 )
  W( 80 ) = W( 80 ) + a*JVS( 452 )
  W( 81 ) = W( 81 ) + a*JVS( 453 )
  JVS( 454) = W( 20 )
  JVS( 455) = W( 66 )
  JVS( 456) = W( 69 )
  JVS( 457) = W( 70 )
  JVS( 458) = W( 71 )
  JVS( 459) = W( 75 )
  JVS( 460) = W( 76 )
  JVS( 461) = W( 77 )
  JVS( 462) = W( 78 )
  JVS( 463) = W( 79 )
  JVS( 464) = W( 80 )
  JVS( 465) = W( 81 )
  IF ( ABS(  JVS( 468 )) < TINY(a) ) THEN
         IER = 71                                      
         RETURN
  END IF
   W( 42 ) = JVS( 466 )
   W( 59 ) = JVS( 467 )
   W( 71 ) = JVS( 468 )
   W( 75 ) = JVS( 469 )
   W( 76 ) = JVS( 470 )
   W( 77 ) = JVS( 471 )
   W( 78 ) = JVS( 472 )
   W( 79 ) = JVS( 473 )
   W( 80 ) = JVS( 474 )
   W( 81 ) = JVS( 475 )
  a = -W( 42 ) / JVS(          174  )
  W( 42 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 175 )
  W( 75 ) = W( 75 ) + a*JVS( 176 )
  W( 79 ) = W( 79 ) + a*JVS( 177 )
  a = -W( 59 ) / JVS(          309  )
  W( 59 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 310 )
  W( 80 ) = W( 80 ) + a*JVS( 311 )
  W( 81 ) = W( 81 ) + a*JVS( 312 )
  JVS( 466) = W( 42 )
  JVS( 467) = W( 59 )
  JVS( 468) = W( 71 )
  JVS( 469) = W( 75 )
  JVS( 470) = W( 76 )
  JVS( 471) = W( 77 )
  JVS( 472) = W( 78 )
  JVS( 473) = W( 79 )
  JVS( 474) = W( 80 )
  JVS( 475) = W( 81 )
  IF ( ABS(  JVS( 508 )) < TINY(a) ) THEN
         IER = 72                                      
         RETURN
  END IF
   W( 12 ) = JVS( 476 )
   W( 18 ) = JVS( 477 )
   W( 23 ) = JVS( 478 )
   W( 24 ) = JVS( 479 )
   W( 27 ) = JVS( 480 )
   W( 30 ) = JVS( 481 )
   W( 32 ) = JVS( 482 )
   W( 35 ) = JVS( 483 )
   W( 36 ) = JVS( 484 )
   W( 40 ) = JVS( 485 )
   W( 41 ) = JVS( 486 )
   W( 42 ) = JVS( 487 )
   W( 43 ) = JVS( 488 )
   W( 47 ) = JVS( 489 )
   W( 48 ) = JVS( 490 )
   W( 50 ) = JVS( 491 )
   W( 53 ) = JVS( 492 )
   W( 54 ) = JVS( 493 )
   W( 55 ) = JVS( 494 )
   W( 56 ) = JVS( 495 )
   W( 57 ) = JVS( 496 )
   W( 59 ) = JVS( 497 )
   W( 60 ) = JVS( 498 )
   W( 61 ) = JVS( 499 )
   W( 63 ) = JVS( 500 )
   W( 64 ) = JVS( 501 )
   W( 65 ) = JVS( 502 )
   W( 66 ) = JVS( 503 )
   W( 68 ) = JVS( 504 )
   W( 69 ) = JVS( 505 )
   W( 70 ) = JVS( 506 )
   W( 71 ) = JVS( 507 )
   W( 72 ) = JVS( 508 )
   W( 73 ) = JVS( 509 )
   W( 74 ) = JVS( 510 )
   W( 75 ) = JVS( 511 )
   W( 76 ) = JVS( 512 )
   W( 77 ) = JVS( 513 )
   W( 78 ) = JVS( 514 )
   W( 79 ) = JVS( 515 )
   W( 80 ) = JVS( 516 )
   W( 81 ) = JVS( 517 )
  a = -W( 12 ) / JVS(           28  )
  W( 12 ) = -a
  W( 30 ) = W( 30 ) + a*JVS( 29 )
  W( 76 ) = W( 76 ) + a*JVS( 30 )
  a = -W( 18 ) / JVS(           49  )
  W( 18 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 50 )
  W( 77 ) = W( 77 ) + a*JVS( 51 )
  W( 79 ) = W( 79 ) + a*JVS( 52 )
  a = -W( 23 ) / JVS(           70  )
  W( 23 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 71 )
  W( 75 ) = W( 75 ) + a*JVS( 72 )
  W( 79 ) = W( 79 ) + a*JVS( 73 )
  a = -W( 24 ) / JVS(           75  )
  W( 24 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 76 )
  W( 76 ) = W( 76 ) + a*JVS( 77 )
  a = -W( 27 ) / JVS(           86  )
  W( 27 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 87 )
  W( 81 ) = W( 81 ) + a*JVS( 88 )
  a = -W( 30 ) / JVS(          104  )
  W( 30 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 105 )
  W( 76 ) = W( 76 ) + a*JVS( 106 )
  W( 81 ) = W( 81 ) + a*JVS( 107 )
  a = -W( 32 ) / JVS(          112  )
  W( 32 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 113 )
  W( 75 ) = W( 75 ) + a*JVS( 114 )
  W( 78 ) = W( 78 ) + a*JVS( 115 )
  W( 79 ) = W( 79 ) + a*JVS( 116 )
  a = -W( 35 ) / JVS(          132  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 133 )
  W( 75 ) = W( 75 ) + a*JVS( 134 )
  W( 79 ) = W( 79 ) + a*JVS( 135 )
  a = -W( 36 ) / JVS(          136  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 137 )
  W( 75 ) = W( 75 ) + a*JVS( 138 )
  W( 78 ) = W( 78 ) + a*JVS( 139 )
  a = -W( 40 ) / JVS(          161  )
  W( 40 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 162 )
  W( 61 ) = W( 61 ) + a*JVS( 163 )
  W( 75 ) = W( 75 ) + a*JVS( 164 )
  W( 81 ) = W( 81 ) + a*JVS( 165 )
  a = -W( 41 ) / JVS(          166  )
  W( 41 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 167 )
  W( 63 ) = W( 63 ) + a*JVS( 168 )
  W( 68 ) = W( 68 ) + a*JVS( 169 )
  W( 70 ) = W( 70 ) + a*JVS( 170 )
  W( 71 ) = W( 71 ) + a*JVS( 171 )
  W( 75 ) = W( 75 ) + a*JVS( 172 )
  W( 77 ) = W( 77 ) + a*JVS( 173 )
  a = -W( 42 ) / JVS(          174  )
  W( 42 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 175 )
  W( 75 ) = W( 75 ) + a*JVS( 176 )
  W( 79 ) = W( 79 ) + a*JVS( 177 )
  a = -W( 43 ) / JVS(          178  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 179 )
  W( 74 ) = W( 74 ) + a*JVS( 180 )
  W( 75 ) = W( 75 ) + a*JVS( 181 )
  a = -W( 47 ) / JVS(          203  )
  W( 47 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 204 )
  W( 72 ) = W( 72 ) + a*JVS( 205 )
  W( 75 ) = W( 75 ) + a*JVS( 206 )
  W( 81 ) = W( 81 ) + a*JVS( 207 )
  a = -W( 48 ) / JVS(          208  )
  W( 48 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 209 )
  W( 75 ) = W( 75 ) + a*JVS( 210 )
  W( 79 ) = W( 79 ) + a*JVS( 211 )
  a = -W( 50 ) / JVS(          221  )
  W( 50 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 222 )
  W( 75 ) = W( 75 ) + a*JVS( 223 )
  W( 76 ) = W( 76 ) + a*JVS( 224 )
  W( 79 ) = W( 79 ) + a*JVS( 225 )
  a = -W( 53 ) / JVS(          248  )
  W( 53 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 249 )
  W( 70 ) = W( 70 ) + a*JVS( 250 )
  W( 75 ) = W( 75 ) + a*JVS( 251 )
  W( 76 ) = W( 76 ) + a*JVS( 252 )
  W( 77 ) = W( 77 ) + a*JVS( 253 )
  W( 78 ) = W( 78 ) + a*JVS( 254 )
  W( 80 ) = W( 80 ) + a*JVS( 255 )
  W( 81 ) = W( 81 ) + a*JVS( 256 )
  a = -W( 54 ) / JVS(          259  )
  W( 54 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 260 )
  W( 76 ) = W( 76 ) + a*JVS( 261 )
  W( 79 ) = W( 79 ) + a*JVS( 262 )
  a = -W( 55 ) / JVS(          268  )
  W( 55 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 269 )
  W( 76 ) = W( 76 ) + a*JVS( 270 )
  W( 77 ) = W( 77 ) + a*JVS( 271 )
  W( 79 ) = W( 79 ) + a*JVS( 272 )
  a = -W( 56 ) / JVS(          273  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 274 )
  W( 76 ) = W( 76 ) + a*JVS( 275 )
  W( 79 ) = W( 79 ) + a*JVS( 276 )
  W( 80 ) = W( 80 ) + a*JVS( 277 )
  a = -W( 57 ) / JVS(          280  )
  W( 57 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 281 )
  W( 76 ) = W( 76 ) + a*JVS( 282 )
  W( 77 ) = W( 77 ) + a*JVS( 283 )
  W( 79 ) = W( 79 ) + a*JVS( 284 )
  a = -W( 59 ) / JVS(          309  )
  W( 59 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 310 )
  W( 80 ) = W( 80 ) + a*JVS( 311 )
  W( 81 ) = W( 81 ) + a*JVS( 312 )
  a = -W( 60 ) / JVS(          315  )
  W( 60 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 316 )
  W( 71 ) = W( 71 ) + a*JVS( 317 )
  W( 75 ) = W( 75 ) + a*JVS( 318 )
  W( 76 ) = W( 76 ) + a*JVS( 319 )
  W( 79 ) = W( 79 ) + a*JVS( 320 )
  W( 80 ) = W( 80 ) + a*JVS( 321 )
  W( 81 ) = W( 81 ) + a*JVS( 322 )
  a = -W( 61 ) / JVS(          324  )
  W( 61 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 325 )
  W( 75 ) = W( 75 ) + a*JVS( 326 )
  W( 80 ) = W( 80 ) + a*JVS( 327 )
  W( 81 ) = W( 81 ) + a*JVS( 328 )
  a = -W( 63 ) / JVS(          353  )
  W( 63 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  W( 76 ) = W( 76 ) + a*JVS( 356 )
  W( 77 ) = W( 77 ) + a*JVS( 357 )
  W( 79 ) = W( 79 ) + a*JVS( 358 )
  a = -W( 64 ) / JVS(          361  )
  W( 64 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 362 )
  W( 76 ) = W( 76 ) + a*JVS( 363 )
  W( 79 ) = W( 79 ) + a*JVS( 364 )
  W( 80 ) = W( 80 ) + a*JVS( 365 )
  W( 81 ) = W( 81 ) + a*JVS( 366 )
  a = -W( 65 ) / JVS(          373  )
  W( 65 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 374 )
  W( 69 ) = W( 69 ) + a*JVS( 375 )
  W( 70 ) = W( 70 ) + a*JVS( 376 )
  W( 73 ) = W( 73 ) + a*JVS( 377 )
  W( 74 ) = W( 74 ) + a*JVS( 378 )
  W( 75 ) = W( 75 ) + a*JVS( 379 )
  W( 76 ) = W( 76 ) + a*JVS( 380 )
  W( 77 ) = W( 77 ) + a*JVS( 381 )
  W( 78 ) = W( 78 ) + a*JVS( 382 )
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 80 ) = W( 80 ) + a*JVS( 384 )
  W( 81 ) = W( 81 ) + a*JVS( 385 )
  a = -W( 66 ) / JVS(          392  )
  W( 66 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 393 )
  W( 75 ) = W( 75 ) + a*JVS( 394 )
  W( 76 ) = W( 76 ) + a*JVS( 395 )
  W( 77 ) = W( 77 ) + a*JVS( 396 )
  W( 78 ) = W( 78 ) + a*JVS( 397 )
  W( 79 ) = W( 79 ) + a*JVS( 398 )
  W( 80 ) = W( 80 ) + a*JVS( 399 )
  W( 81 ) = W( 81 ) + a*JVS( 400 )
  a = -W( 68 ) / JVS(          429  )
  W( 68 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 430 )
  W( 71 ) = W( 71 ) + a*JVS( 431 )
  W( 75 ) = W( 75 ) + a*JVS( 432 )
  W( 76 ) = W( 76 ) + a*JVS( 433 )
  W( 77 ) = W( 77 ) + a*JVS( 434 )
  W( 78 ) = W( 78 ) + a*JVS( 435 )
  W( 79 ) = W( 79 ) + a*JVS( 436 )
  W( 80 ) = W( 80 ) + a*JVS( 437 )
  W( 81 ) = W( 81 ) + a*JVS( 438 )
  a = -W( 69 ) / JVS(          445  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 446 )
  W( 75 ) = W( 75 ) + a*JVS( 447 )
  W( 76 ) = W( 76 ) + a*JVS( 448 )
  W( 77 ) = W( 77 ) + a*JVS( 449 )
  W( 78 ) = W( 78 ) + a*JVS( 450 )
  W( 79 ) = W( 79 ) + a*JVS( 451 )
  W( 80 ) = W( 80 ) + a*JVS( 452 )
  W( 81 ) = W( 81 ) + a*JVS( 453 )
  a = -W( 70 ) / JVS(          457  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 75 ) = W( 75 ) + a*JVS( 459 )
  W( 76 ) = W( 76 ) + a*JVS( 460 )
  W( 77 ) = W( 77 ) + a*JVS( 461 )
  W( 78 ) = W( 78 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  W( 81 ) = W( 81 ) + a*JVS( 465 )
  a = -W( 71 ) / JVS(          468  )
  W( 71 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 469 )
  W( 76 ) = W( 76 ) + a*JVS( 470 )
  W( 77 ) = W( 77 ) + a*JVS( 471 )
  W( 78 ) = W( 78 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  W( 80 ) = W( 80 ) + a*JVS( 474 )
  W( 81 ) = W( 81 ) + a*JVS( 475 )
  JVS( 476) = W( 12 )
  JVS( 477) = W( 18 )
  JVS( 478) = W( 23 )
  JVS( 479) = W( 24 )
  JVS( 480) = W( 27 )
  JVS( 481) = W( 30 )
  JVS( 482) = W( 32 )
  JVS( 483) = W( 35 )
  JVS( 484) = W( 36 )
  JVS( 485) = W( 40 )
  JVS( 486) = W( 41 )
  JVS( 487) = W( 42 )
  JVS( 488) = W( 43 )
  JVS( 489) = W( 47 )
  JVS( 490) = W( 48 )
  JVS( 491) = W( 50 )
  JVS( 492) = W( 53 )
  JVS( 493) = W( 54 )
  JVS( 494) = W( 55 )
  JVS( 495) = W( 56 )
  JVS( 496) = W( 57 )
  JVS( 497) = W( 59 )
  JVS( 498) = W( 60 )
  JVS( 499) = W( 61 )
  JVS( 500) = W( 63 )
  JVS( 501) = W( 64 )
  JVS( 502) = W( 65 )
  JVS( 503) = W( 66 )
  JVS( 504) = W( 68 )
  JVS( 505) = W( 69 )
  JVS( 506) = W( 70 )
  JVS( 507) = W( 71 )
  JVS( 508) = W( 72 )
  JVS( 509) = W( 73 )
  JVS( 510) = W( 74 )
  JVS( 511) = W( 75 )
  JVS( 512) = W( 76 )
  JVS( 513) = W( 77 )
  JVS( 514) = W( 78 )
  JVS( 515) = W( 79 )
  JVS( 516) = W( 80 )
  JVS( 517) = W( 81 )
  IF ( ABS(  JVS( 548 )) < TINY(a) ) THEN
         IER = 73                                      
         RETURN
  END IF
   W( 13 ) = JVS( 518 )
   W( 14 ) = JVS( 519 )
   W( 24 ) = JVS( 520 )
   W( 26 ) = JVS( 521 )
   W( 29 ) = JVS( 522 )
   W( 30 ) = JVS( 523 )
   W( 36 ) = JVS( 524 )
   W( 39 ) = JVS( 525 )
   W( 43 ) = JVS( 526 )
   W( 44 ) = JVS( 527 )
   W( 46 ) = JVS( 528 )
   W( 49 ) = JVS( 529 )
   W( 50 ) = JVS( 530 )
   W( 52 ) = JVS( 531 )
   W( 54 ) = JVS( 532 )
   W( 55 ) = JVS( 533 )
   W( 56 ) = JVS( 534 )
   W( 57 ) = JVS( 535 )
   W( 59 ) = JVS( 536 )
   W( 60 ) = JVS( 537 )
   W( 61 ) = JVS( 538 )
   W( 62 ) = JVS( 539 )
   W( 63 ) = JVS( 540 )
   W( 64 ) = JVS( 541 )
   W( 67 ) = JVS( 542 )
   W( 68 ) = JVS( 543 )
   W( 69 ) = JVS( 544 )
   W( 70 ) = JVS( 545 )
   W( 71 ) = JVS( 546 )
   W( 72 ) = JVS( 547 )
   W( 73 ) = JVS( 548 )
   W( 74 ) = JVS( 549 )
   W( 75 ) = JVS( 550 )
   W( 76 ) = JVS( 551 )
   W( 77 ) = JVS( 552 )
   W( 78 ) = JVS( 553 )
   W( 79 ) = JVS( 554 )
   W( 80 ) = JVS( 555 )
   W( 81 ) = JVS( 556 )
  a = -W( 13 ) / JVS(           31  )
  W( 13 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 32 )
  W( 80 ) = W( 80 ) + a*JVS( 33 )
  a = -W( 14 ) / JVS(           35  )
  W( 14 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 36 )
  W( 75 ) = W( 75 ) + a*JVS( 37 )
  a = -W( 24 ) / JVS(           75  )
  W( 24 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 76 )
  W( 76 ) = W( 76 ) + a*JVS( 77 )
  a = -W( 26 ) / JVS(           82  )
  W( 26 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 83 )
  W( 75 ) = W( 75 ) + a*JVS( 84 )
  W( 79 ) = W( 79 ) + a*JVS( 85 )
  a = -W( 29 ) / JVS(           97  )
  W( 29 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 98 )
  W( 61 ) = W( 61 ) + a*JVS( 99 )
  W( 75 ) = W( 75 ) + a*JVS( 100 )
  W( 76 ) = W( 76 ) + a*JVS( 101 )
  W( 80 ) = W( 80 ) + a*JVS( 102 )
  a = -W( 30 ) / JVS(          104  )
  W( 30 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 105 )
  W( 76 ) = W( 76 ) + a*JVS( 106 )
  W( 81 ) = W( 81 ) + a*JVS( 107 )
  a = -W( 36 ) / JVS(          136  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 137 )
  W( 75 ) = W( 75 ) + a*JVS( 138 )
  W( 78 ) = W( 78 ) + a*JVS( 139 )
  a = -W( 39 ) / JVS(          154  )
  W( 39 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 155 )
  W( 67 ) = W( 67 ) + a*JVS( 156 )
  W( 72 ) = W( 72 ) + a*JVS( 157 )
  W( 73 ) = W( 73 ) + a*JVS( 158 )
  W( 75 ) = W( 75 ) + a*JVS( 159 )
  W( 80 ) = W( 80 ) + a*JVS( 160 )
  a = -W( 43 ) / JVS(          178  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 179 )
  W( 74 ) = W( 74 ) + a*JVS( 180 )
  W( 75 ) = W( 75 ) + a*JVS( 181 )
  a = -W( 44 ) / JVS(          184  )
  W( 44 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 185 )
  W( 76 ) = W( 76 ) + a*JVS( 186 )
  W( 79 ) = W( 79 ) + a*JVS( 187 )
  a = -W( 46 ) / JVS(          196  )
  W( 46 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 197 )
  W( 80 ) = W( 80 ) + a*JVS( 198 )
  W( 81 ) = W( 81 ) + a*JVS( 199 )
  a = -W( 49 ) / JVS(          215  )
  W( 49 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 216 )
  W( 75 ) = W( 75 ) + a*JVS( 217 )
  W( 76 ) = W( 76 ) + a*JVS( 218 )
  W( 79 ) = W( 79 ) + a*JVS( 219 )
  a = -W( 50 ) / JVS(          221  )
  W( 50 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 222 )
  W( 75 ) = W( 75 ) + a*JVS( 223 )
  W( 76 ) = W( 76 ) + a*JVS( 224 )
  W( 79 ) = W( 79 ) + a*JVS( 225 )
  a = -W( 52 ) / JVS(          239  )
  W( 52 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 240 )
  W( 72 ) = W( 72 ) + a*JVS( 241 )
  W( 73 ) = W( 73 ) + a*JVS( 242 )
  W( 75 ) = W( 75 ) + a*JVS( 243 )
  W( 79 ) = W( 79 ) + a*JVS( 244 )
  W( 81 ) = W( 81 ) + a*JVS( 245 )
  a = -W( 54 ) / JVS(          259  )
  W( 54 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 260 )
  W( 76 ) = W( 76 ) + a*JVS( 261 )
  W( 79 ) = W( 79 ) + a*JVS( 262 )
  a = -W( 55 ) / JVS(          268  )
  W( 55 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 269 )
  W( 76 ) = W( 76 ) + a*JVS( 270 )
  W( 77 ) = W( 77 ) + a*JVS( 271 )
  W( 79 ) = W( 79 ) + a*JVS( 272 )
  a = -W( 56 ) / JVS(          273  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 274 )
  W( 76 ) = W( 76 ) + a*JVS( 275 )
  W( 79 ) = W( 79 ) + a*JVS( 276 )
  W( 80 ) = W( 80 ) + a*JVS( 277 )
  a = -W( 57 ) / JVS(          280  )
  W( 57 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 281 )
  W( 76 ) = W( 76 ) + a*JVS( 282 )
  W( 77 ) = W( 77 ) + a*JVS( 283 )
  W( 79 ) = W( 79 ) + a*JVS( 284 )
  a = -W( 59 ) / JVS(          309  )
  W( 59 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 310 )
  W( 80 ) = W( 80 ) + a*JVS( 311 )
  W( 81 ) = W( 81 ) + a*JVS( 312 )
  a = -W( 60 ) / JVS(          315  )
  W( 60 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 316 )
  W( 71 ) = W( 71 ) + a*JVS( 317 )
  W( 75 ) = W( 75 ) + a*JVS( 318 )
  W( 76 ) = W( 76 ) + a*JVS( 319 )
  W( 79 ) = W( 79 ) + a*JVS( 320 )
  W( 80 ) = W( 80 ) + a*JVS( 321 )
  W( 81 ) = W( 81 ) + a*JVS( 322 )
  a = -W( 61 ) / JVS(          324  )
  W( 61 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 325 )
  W( 75 ) = W( 75 ) + a*JVS( 326 )
  W( 80 ) = W( 80 ) + a*JVS( 327 )
  W( 81 ) = W( 81 ) + a*JVS( 328 )
  a = -W( 62 ) / JVS(          341  )
  W( 62 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  W( 76 ) = W( 76 ) + a*JVS( 344 )
  W( 77 ) = W( 77 ) + a*JVS( 345 )
  W( 79 ) = W( 79 ) + a*JVS( 346 )
  W( 80 ) = W( 80 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  a = -W( 63 ) / JVS(          353  )
  W( 63 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  W( 76 ) = W( 76 ) + a*JVS( 356 )
  W( 77 ) = W( 77 ) + a*JVS( 357 )
  W( 79 ) = W( 79 ) + a*JVS( 358 )
  a = -W( 64 ) / JVS(          361  )
  W( 64 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 362 )
  W( 76 ) = W( 76 ) + a*JVS( 363 )
  W( 79 ) = W( 79 ) + a*JVS( 364 )
  W( 80 ) = W( 80 ) + a*JVS( 365 )
  W( 81 ) = W( 81 ) + a*JVS( 366 )
  a = -W( 67 ) / JVS(          411  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 412 )
  W( 69 ) = W( 69 ) + a*JVS( 413 )
  W( 70 ) = W( 70 ) + a*JVS( 414 )
  W( 71 ) = W( 71 ) + a*JVS( 415 )
  W( 73 ) = W( 73 ) + a*JVS( 416 )
  W( 74 ) = W( 74 ) + a*JVS( 417 )
  W( 75 ) = W( 75 ) + a*JVS( 418 )
  W( 76 ) = W( 76 ) + a*JVS( 419 )
  W( 77 ) = W( 77 ) + a*JVS( 420 )
  W( 78 ) = W( 78 ) + a*JVS( 421 )
  W( 79 ) = W( 79 ) + a*JVS( 422 )
  W( 80 ) = W( 80 ) + a*JVS( 423 )
  W( 81 ) = W( 81 ) + a*JVS( 424 )
  a = -W( 68 ) / JVS(          429  )
  W( 68 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 430 )
  W( 71 ) = W( 71 ) + a*JVS( 431 )
  W( 75 ) = W( 75 ) + a*JVS( 432 )
  W( 76 ) = W( 76 ) + a*JVS( 433 )
  W( 77 ) = W( 77 ) + a*JVS( 434 )
  W( 78 ) = W( 78 ) + a*JVS( 435 )
  W( 79 ) = W( 79 ) + a*JVS( 436 )
  W( 80 ) = W( 80 ) + a*JVS( 437 )
  W( 81 ) = W( 81 ) + a*JVS( 438 )
  a = -W( 69 ) / JVS(          445  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 446 )
  W( 75 ) = W( 75 ) + a*JVS( 447 )
  W( 76 ) = W( 76 ) + a*JVS( 448 )
  W( 77 ) = W( 77 ) + a*JVS( 449 )
  W( 78 ) = W( 78 ) + a*JVS( 450 )
  W( 79 ) = W( 79 ) + a*JVS( 451 )
  W( 80 ) = W( 80 ) + a*JVS( 452 )
  W( 81 ) = W( 81 ) + a*JVS( 453 )
  a = -W( 70 ) / JVS(          457  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 75 ) = W( 75 ) + a*JVS( 459 )
  W( 76 ) = W( 76 ) + a*JVS( 460 )
  W( 77 ) = W( 77 ) + a*JVS( 461 )
  W( 78 ) = W( 78 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  W( 81 ) = W( 81 ) + a*JVS( 465 )
  a = -W( 71 ) / JVS(          468  )
  W( 71 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 469 )
  W( 76 ) = W( 76 ) + a*JVS( 470 )
  W( 77 ) = W( 77 ) + a*JVS( 471 )
  W( 78 ) = W( 78 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  W( 80 ) = W( 80 ) + a*JVS( 474 )
  W( 81 ) = W( 81 ) + a*JVS( 475 )
  a = -W( 72 ) / JVS(          508  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 509 )
  W( 74 ) = W( 74 ) + a*JVS( 510 )
  W( 75 ) = W( 75 ) + a*JVS( 511 )
  W( 76 ) = W( 76 ) + a*JVS( 512 )
  W( 77 ) = W( 77 ) + a*JVS( 513 )
  W( 78 ) = W( 78 ) + a*JVS( 514 )
  W( 79 ) = W( 79 ) + a*JVS( 515 )
  W( 80 ) = W( 80 ) + a*JVS( 516 )
  W( 81 ) = W( 81 ) + a*JVS( 517 )
  JVS( 518) = W( 13 )
  JVS( 519) = W( 14 )
  JVS( 520) = W( 24 )
  JVS( 521) = W( 26 )
  JVS( 522) = W( 29 )
  JVS( 523) = W( 30 )
  JVS( 524) = W( 36 )
  JVS( 525) = W( 39 )
  JVS( 526) = W( 43 )
  JVS( 527) = W( 44 )
  JVS( 528) = W( 46 )
  JVS( 529) = W( 49 )
  JVS( 530) = W( 50 )
  JVS( 531) = W( 52 )
  JVS( 532) = W( 54 )
  JVS( 533) = W( 55 )
  JVS( 534) = W( 56 )
  JVS( 535) = W( 57 )
  JVS( 536) = W( 59 )
  JVS( 537) = W( 60 )
  JVS( 538) = W( 61 )
  JVS( 539) = W( 62 )
  JVS( 540) = W( 63 )
  JVS( 541) = W( 64 )
  JVS( 542) = W( 67 )
  JVS( 543) = W( 68 )
  JVS( 544) = W( 69 )
  JVS( 545) = W( 70 )
  JVS( 546) = W( 71 )
  JVS( 547) = W( 72 )
  JVS( 548) = W( 73 )
  JVS( 549) = W( 74 )
  JVS( 550) = W( 75 )
  JVS( 551) = W( 76 )
  JVS( 552) = W( 77 )
  JVS( 553) = W( 78 )
  JVS( 554) = W( 79 )
  JVS( 555) = W( 80 )
  JVS( 556) = W( 81 )
  IF ( ABS(  JVS( 564 )) < TINY(a) ) THEN
         IER = 74                                      
         RETURN
  END IF
   W( 20 ) = JVS( 557 )
   W( 43 ) = JVS( 558 )
   W( 59 ) = JVS( 559 )
   W( 66 ) = JVS( 560 )
   W( 70 ) = JVS( 561 )
   W( 71 ) = JVS( 562 )
   W( 73 ) = JVS( 563 )
   W( 74 ) = JVS( 564 )
   W( 75 ) = JVS( 565 )
   W( 76 ) = JVS( 566 )
   W( 77 ) = JVS( 567 )
   W( 78 ) = JVS( 568 )
   W( 79 ) = JVS( 569 )
   W( 80 ) = JVS( 570 )
   W( 81 ) = JVS( 571 )
  a = -W( 20 ) / JVS(           58  )
  W( 20 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 59 )
  W( 75 ) = W( 75 ) + a*JVS( 60 )
  W( 79 ) = W( 79 ) + a*JVS( 61 )
  a = -W( 43 ) / JVS(          178  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 179 )
  W( 74 ) = W( 74 ) + a*JVS( 180 )
  W( 75 ) = W( 75 ) + a*JVS( 181 )
  a = -W( 59 ) / JVS(          309  )
  W( 59 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 310 )
  W( 80 ) = W( 80 ) + a*JVS( 311 )
  W( 81 ) = W( 81 ) + a*JVS( 312 )
  a = -W( 66 ) / JVS(          392  )
  W( 66 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 393 )
  W( 75 ) = W( 75 ) + a*JVS( 394 )
  W( 76 ) = W( 76 ) + a*JVS( 395 )
  W( 77 ) = W( 77 ) + a*JVS( 396 )
  W( 78 ) = W( 78 ) + a*JVS( 397 )
  W( 79 ) = W( 79 ) + a*JVS( 398 )
  W( 80 ) = W( 80 ) + a*JVS( 399 )
  W( 81 ) = W( 81 ) + a*JVS( 400 )
  a = -W( 70 ) / JVS(          457  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 75 ) = W( 75 ) + a*JVS( 459 )
  W( 76 ) = W( 76 ) + a*JVS( 460 )
  W( 77 ) = W( 77 ) + a*JVS( 461 )
  W( 78 ) = W( 78 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  W( 81 ) = W( 81 ) + a*JVS( 465 )
  a = -W( 71 ) / JVS(          468  )
  W( 71 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 469 )
  W( 76 ) = W( 76 ) + a*JVS( 470 )
  W( 77 ) = W( 77 ) + a*JVS( 471 )
  W( 78 ) = W( 78 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  W( 80 ) = W( 80 ) + a*JVS( 474 )
  W( 81 ) = W( 81 ) + a*JVS( 475 )
  a = -W( 73 ) / JVS(          548  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 76 ) = W( 76 ) + a*JVS( 551 )
  W( 77 ) = W( 77 ) + a*JVS( 552 )
  W( 78 ) = W( 78 ) + a*JVS( 553 )
  W( 79 ) = W( 79 ) + a*JVS( 554 )
  W( 80 ) = W( 80 ) + a*JVS( 555 )
  W( 81 ) = W( 81 ) + a*JVS( 556 )
  JVS( 557) = W( 20 )
  JVS( 558) = W( 43 )
  JVS( 559) = W( 59 )
  JVS( 560) = W( 66 )
  JVS( 561) = W( 70 )
  JVS( 562) = W( 71 )
  JVS( 563) = W( 73 )
  JVS( 564) = W( 74 )
  JVS( 565) = W( 75 )
  JVS( 566) = W( 76 )
  JVS( 567) = W( 77 )
  JVS( 568) = W( 78 )
  JVS( 569) = W( 79 )
  JVS( 570) = W( 80 )
  JVS( 571) = W( 81 )
  IF ( ABS(  JVS( 639 )) < TINY(a) ) THEN
         IER = 75                                      
         RETURN
  END IF
   W( 2 ) = JVS( 572 )
   W( 4 ) = JVS( 573 )
   W( 5 ) = JVS( 574 )
   W( 6 ) = JVS( 575 )
   W( 7 ) = JVS( 576 )
   W( 8 ) = JVS( 577 )
   W( 9 ) = JVS( 578 )
   W( 10 ) = JVS( 579 )
   W( 11 ) = JVS( 580 )
   W( 15 ) = JVS( 581 )
   W( 16 ) = JVS( 582 )
   W( 17 ) = JVS( 583 )
   W( 18 ) = JVS( 584 )
   W( 19 ) = JVS( 585 )
   W( 20 ) = JVS( 586 )
   W( 21 ) = JVS( 587 )
   W( 22 ) = JVS( 588 )
   W( 23 ) = JVS( 589 )
   W( 25 ) = JVS( 590 )
   W( 26 ) = JVS( 591 )
   W( 27 ) = JVS( 592 )
   W( 28 ) = JVS( 593 )
   W( 29 ) = JVS( 594 )
   W( 31 ) = JVS( 595 )
   W( 32 ) = JVS( 596 )
   W( 33 ) = JVS( 597 )
   W( 34 ) = JVS( 598 )
   W( 35 ) = JVS( 599 )
   W( 36 ) = JVS( 600 )
   W( 37 ) = JVS( 601 )
   W( 38 ) = JVS( 602 )
   W( 39 ) = JVS( 603 )
   W( 40 ) = JVS( 604 )
   W( 41 ) = JVS( 605 )
   W( 42 ) = JVS( 606 )
   W( 43 ) = JVS( 607 )
   W( 44 ) = JVS( 608 )
   W( 45 ) = JVS( 609 )
   W( 46 ) = JVS( 610 )
   W( 47 ) = JVS( 611 )
   W( 48 ) = JVS( 612 )
   W( 49 ) = JVS( 613 )
   W( 50 ) = JVS( 614 )
   W( 51 ) = JVS( 615 )
   W( 52 ) = JVS( 616 )
   W( 53 ) = JVS( 617 )
   W( 54 ) = JVS( 618 )
   W( 55 ) = JVS( 619 )
   W( 56 ) = JVS( 620 )
   W( 57 ) = JVS( 621 )
   W( 58 ) = JVS( 622 )
   W( 59 ) = JVS( 623 )
   W( 60 ) = JVS( 624 )
   W( 61 ) = JVS( 625 )
   W( 62 ) = JVS( 626 )
   W( 63 ) = JVS( 627 )
   W( 64 ) = JVS( 628 )
   W( 65 ) = JVS( 629 )
   W( 66 ) = JVS( 630 )
   W( 67 ) = JVS( 631 )
   W( 68 ) = JVS( 632 )
   W( 69 ) = JVS( 633 )
   W( 70 ) = JVS( 634 )
   W( 71 ) = JVS( 635 )
   W( 72 ) = JVS( 636 )
   W( 73 ) = JVS( 637 )
   W( 74 ) = JVS( 638 )
   W( 75 ) = JVS( 639 )
   W( 76 ) = JVS( 640 )
   W( 77 ) = JVS( 641 )
   W( 78 ) = JVS( 642 )
   W( 79 ) = JVS( 643 )
   W( 80 ) = JVS( 644 )
   W( 81 ) = JVS( 645 )
  a = -W( 2 ) / JVS(            4  )
  W( 2 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 5 )
  a = -W( 4 ) / JVS(            8  )
  W( 4 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 9 )
  a = -W( 5 ) / JVS(           10  )
  W( 5 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 11 )
  a = -W( 6 ) / JVS(           12  )
  W( 6 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 13 )
  a = -W( 7 ) / JVS(           14  )
  W( 7 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 15 )
  a = -W( 8 ) / JVS(           16  )
  W( 8 ) = -a
  W( 16 ) = W( 16 ) + a*JVS( 17 )
  W( 75 ) = W( 75 ) + a*JVS( 18 )
  W( 80 ) = W( 80 ) + a*JVS( 19 )
  a = -W( 9 ) / JVS(           20  )
  W( 9 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 21 )
  W( 79 ) = W( 79 ) + a*JVS( 22 )
  a = -W( 10 ) / JVS(           23  )
  W( 10 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 24 )
  a = -W( 11 ) / JVS(           26  )
  W( 11 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 27 )
  a = -W( 15 ) / JVS(           38  )
  W( 15 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 39 )
  W( 75 ) = W( 75 ) + a*JVS( 40 )
  W( 79 ) = W( 79 ) + a*JVS( 41 )
  a = -W( 16 ) / JVS(           42  )
  W( 16 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 43 )
  W( 80 ) = W( 80 ) + a*JVS( 44 )
  a = -W( 17 ) / JVS(           45  )
  W( 17 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 46 )
  W( 75 ) = W( 75 ) + a*JVS( 47 )
  W( 77 ) = W( 77 ) + a*JVS( 48 )
  a = -W( 18 ) / JVS(           49  )
  W( 18 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 50 )
  W( 77 ) = W( 77 ) + a*JVS( 51 )
  W( 79 ) = W( 79 ) + a*JVS( 52 )
  a = -W( 19 ) / JVS(           53  )
  W( 19 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 54 )
  W( 47 ) = W( 47 ) + a*JVS( 55 )
  W( 72 ) = W( 72 ) + a*JVS( 56 )
  W( 75 ) = W( 75 ) + a*JVS( 57 )
  a = -W( 20 ) / JVS(           58  )
  W( 20 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 59 )
  W( 75 ) = W( 75 ) + a*JVS( 60 )
  W( 79 ) = W( 79 ) + a*JVS( 61 )
  a = -W( 21 ) / JVS(           62  )
  W( 21 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 63 )
  W( 75 ) = W( 75 ) + a*JVS( 64 )
  W( 79 ) = W( 79 ) + a*JVS( 65 )
  a = -W( 22 ) / JVS(           66  )
  W( 22 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 67 )
  W( 75 ) = W( 75 ) + a*JVS( 68 )
  W( 79 ) = W( 79 ) + a*JVS( 69 )
  a = -W( 23 ) / JVS(           70  )
  W( 23 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 71 )
  W( 75 ) = W( 75 ) + a*JVS( 72 )
  W( 79 ) = W( 79 ) + a*JVS( 73 )
  a = -W( 25 ) / JVS(           78  )
  W( 25 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 79 )
  W( 75 ) = W( 75 ) + a*JVS( 80 )
  W( 79 ) = W( 79 ) + a*JVS( 81 )
  a = -W( 26 ) / JVS(           82  )
  W( 26 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 83 )
  W( 75 ) = W( 75 ) + a*JVS( 84 )
  W( 79 ) = W( 79 ) + a*JVS( 85 )
  a = -W( 27 ) / JVS(           86  )
  W( 27 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 87 )
  W( 81 ) = W( 81 ) + a*JVS( 88 )
  a = -W( 28 ) / JVS(           89  )
  W( 28 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 90 )
  W( 71 ) = W( 71 ) + a*JVS( 91 )
  W( 75 ) = W( 75 ) + a*JVS( 92 )
  W( 76 ) = W( 76 ) + a*JVS( 93 )
  W( 77 ) = W( 77 ) + a*JVS( 94 )
  W( 78 ) = W( 78 ) + a*JVS( 95 )
  W( 80 ) = W( 80 ) + a*JVS( 96 )
  a = -W( 29 ) / JVS(           97  )
  W( 29 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 98 )
  W( 61 ) = W( 61 ) + a*JVS( 99 )
  W( 75 ) = W( 75 ) + a*JVS( 100 )
  W( 76 ) = W( 76 ) + a*JVS( 101 )
  W( 80 ) = W( 80 ) + a*JVS( 102 )
  a = -W( 31 ) / JVS(          108  )
  W( 31 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 109 )
  W( 75 ) = W( 75 ) + a*JVS( 110 )
  W( 79 ) = W( 79 ) + a*JVS( 111 )
  a = -W( 32 ) / JVS(          112  )
  W( 32 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 113 )
  W( 75 ) = W( 75 ) + a*JVS( 114 )
  W( 78 ) = W( 78 ) + a*JVS( 115 )
  W( 79 ) = W( 79 ) + a*JVS( 116 )
  a = -W( 33 ) / JVS(          117  )
  W( 33 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 118 )
  W( 54 ) = W( 54 ) + a*JVS( 119 )
  W( 75 ) = W( 75 ) + a*JVS( 120 )
  W( 76 ) = W( 76 ) + a*JVS( 121 )
  a = -W( 34 ) / JVS(          123  )
  W( 34 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 124 )
  W( 61 ) = W( 61 ) + a*JVS( 125 )
  W( 74 ) = W( 74 ) + a*JVS( 126 )
  W( 75 ) = W( 75 ) + a*JVS( 127 )
  W( 77 ) = W( 77 ) + a*JVS( 128 )
  W( 78 ) = W( 78 ) + a*JVS( 129 )
  W( 79 ) = W( 79 ) + a*JVS( 130 )
  W( 81 ) = W( 81 ) + a*JVS( 131 )
  a = -W( 35 ) / JVS(          132  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 133 )
  W( 75 ) = W( 75 ) + a*JVS( 134 )
  W( 79 ) = W( 79 ) + a*JVS( 135 )
  a = -W( 36 ) / JVS(          136  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 137 )
  W( 75 ) = W( 75 ) + a*JVS( 138 )
  W( 78 ) = W( 78 ) + a*JVS( 139 )
  a = -W( 37 ) / JVS(          141  )
  W( 37 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 142 )
  W( 45 ) = W( 45 ) + a*JVS( 143 )
  W( 53 ) = W( 53 ) + a*JVS( 144 )
  W( 75 ) = W( 75 ) + a*JVS( 145 )
  W( 76 ) = W( 76 ) + a*JVS( 146 )
  W( 79 ) = W( 79 ) + a*JVS( 147 )
  a = -W( 38 ) / JVS(          148  )
  W( 38 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 149 )
  W( 75 ) = W( 75 ) + a*JVS( 150 )
  W( 79 ) = W( 79 ) + a*JVS( 151 )
  a = -W( 39 ) / JVS(          154  )
  W( 39 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 155 )
  W( 67 ) = W( 67 ) + a*JVS( 156 )
  W( 72 ) = W( 72 ) + a*JVS( 157 )
  W( 73 ) = W( 73 ) + a*JVS( 158 )
  W( 75 ) = W( 75 ) + a*JVS( 159 )
  W( 80 ) = W( 80 ) + a*JVS( 160 )
  a = -W( 40 ) / JVS(          161  )
  W( 40 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 162 )
  W( 61 ) = W( 61 ) + a*JVS( 163 )
  W( 75 ) = W( 75 ) + a*JVS( 164 )
  W( 81 ) = W( 81 ) + a*JVS( 165 )
  a = -W( 41 ) / JVS(          166  )
  W( 41 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 167 )
  W( 63 ) = W( 63 ) + a*JVS( 168 )
  W( 68 ) = W( 68 ) + a*JVS( 169 )
  W( 70 ) = W( 70 ) + a*JVS( 170 )
  W( 71 ) = W( 71 ) + a*JVS( 171 )
  W( 75 ) = W( 75 ) + a*JVS( 172 )
  W( 77 ) = W( 77 ) + a*JVS( 173 )
  a = -W( 42 ) / JVS(          174  )
  W( 42 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 175 )
  W( 75 ) = W( 75 ) + a*JVS( 176 )
  W( 79 ) = W( 79 ) + a*JVS( 177 )
  a = -W( 43 ) / JVS(          178  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 179 )
  W( 74 ) = W( 74 ) + a*JVS( 180 )
  W( 75 ) = W( 75 ) + a*JVS( 181 )
  a = -W( 44 ) / JVS(          184  )
  W( 44 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 185 )
  W( 76 ) = W( 76 ) + a*JVS( 186 )
  W( 79 ) = W( 79 ) + a*JVS( 187 )
  a = -W( 45 ) / JVS(          191  )
  W( 45 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 192 )
  W( 75 ) = W( 75 ) + a*JVS( 193 )
  W( 76 ) = W( 76 ) + a*JVS( 194 )
  W( 79 ) = W( 79 ) + a*JVS( 195 )
  a = -W( 46 ) / JVS(          196  )
  W( 46 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 197 )
  W( 80 ) = W( 80 ) + a*JVS( 198 )
  W( 81 ) = W( 81 ) + a*JVS( 199 )
  a = -W( 47 ) / JVS(          203  )
  W( 47 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 204 )
  W( 72 ) = W( 72 ) + a*JVS( 205 )
  W( 75 ) = W( 75 ) + a*JVS( 206 )
  W( 81 ) = W( 81 ) + a*JVS( 207 )
  a = -W( 48 ) / JVS(          208  )
  W( 48 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 209 )
  W( 75 ) = W( 75 ) + a*JVS( 210 )
  W( 79 ) = W( 79 ) + a*JVS( 211 )
  a = -W( 49 ) / JVS(          215  )
  W( 49 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 216 )
  W( 75 ) = W( 75 ) + a*JVS( 217 )
  W( 76 ) = W( 76 ) + a*JVS( 218 )
  W( 79 ) = W( 79 ) + a*JVS( 219 )
  a = -W( 50 ) / JVS(          221  )
  W( 50 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 222 )
  W( 75 ) = W( 75 ) + a*JVS( 223 )
  W( 76 ) = W( 76 ) + a*JVS( 224 )
  W( 79 ) = W( 79 ) + a*JVS( 225 )
  a = -W( 51 ) / JVS(          230  )
  W( 51 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 231 )
  W( 57 ) = W( 57 ) + a*JVS( 232 )
  W( 64 ) = W( 64 ) + a*JVS( 233 )
  W( 75 ) = W( 75 ) + a*JVS( 234 )
  W( 76 ) = W( 76 ) + a*JVS( 235 )
  W( 77 ) = W( 77 ) + a*JVS( 236 )
  W( 79 ) = W( 79 ) + a*JVS( 237 )
  a = -W( 52 ) / JVS(          239  )
  W( 52 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 240 )
  W( 72 ) = W( 72 ) + a*JVS( 241 )
  W( 73 ) = W( 73 ) + a*JVS( 242 )
  W( 75 ) = W( 75 ) + a*JVS( 243 )
  W( 79 ) = W( 79 ) + a*JVS( 244 )
  W( 81 ) = W( 81 ) + a*JVS( 245 )
  a = -W( 53 ) / JVS(          248  )
  W( 53 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 249 )
  W( 70 ) = W( 70 ) + a*JVS( 250 )
  W( 75 ) = W( 75 ) + a*JVS( 251 )
  W( 76 ) = W( 76 ) + a*JVS( 252 )
  W( 77 ) = W( 77 ) + a*JVS( 253 )
  W( 78 ) = W( 78 ) + a*JVS( 254 )
  W( 80 ) = W( 80 ) + a*JVS( 255 )
  W( 81 ) = W( 81 ) + a*JVS( 256 )
  a = -W( 54 ) / JVS(          259  )
  W( 54 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 260 )
  W( 76 ) = W( 76 ) + a*JVS( 261 )
  W( 79 ) = W( 79 ) + a*JVS( 262 )
  a = -W( 55 ) / JVS(          268  )
  W( 55 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 269 )
  W( 76 ) = W( 76 ) + a*JVS( 270 )
  W( 77 ) = W( 77 ) + a*JVS( 271 )
  W( 79 ) = W( 79 ) + a*JVS( 272 )
  a = -W( 56 ) / JVS(          273  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 274 )
  W( 76 ) = W( 76 ) + a*JVS( 275 )
  W( 79 ) = W( 79 ) + a*JVS( 276 )
  W( 80 ) = W( 80 ) + a*JVS( 277 )
  a = -W( 57 ) / JVS(          280  )
  W( 57 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 281 )
  W( 76 ) = W( 76 ) + a*JVS( 282 )
  W( 77 ) = W( 77 ) + a*JVS( 283 )
  W( 79 ) = W( 79 ) + a*JVS( 284 )
  a = -W( 58 ) / JVS(          290  )
  W( 58 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 291 )
  W( 60 ) = W( 60 ) + a*JVS( 292 )
  W( 61 ) = W( 61 ) + a*JVS( 293 )
  W( 62 ) = W( 62 ) + a*JVS( 294 )
  W( 66 ) = W( 66 ) + a*JVS( 295 )
  W( 67 ) = W( 67 ) + a*JVS( 296 )
  W( 68 ) = W( 68 ) + a*JVS( 297 )
  W( 69 ) = W( 69 ) + a*JVS( 298 )
  W( 70 ) = W( 70 ) + a*JVS( 299 )
  W( 72 ) = W( 72 ) + a*JVS( 300 )
  W( 73 ) = W( 73 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  W( 76 ) = W( 76 ) + a*JVS( 303 )
  W( 77 ) = W( 77 ) + a*JVS( 304 )
  W( 78 ) = W( 78 ) + a*JVS( 305 )
  W( 79 ) = W( 79 ) + a*JVS( 306 )
  W( 80 ) = W( 80 ) + a*JVS( 307 )
  W( 81 ) = W( 81 ) + a*JVS( 308 )
  a = -W( 59 ) / JVS(          309  )
  W( 59 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 310 )
  W( 80 ) = W( 80 ) + a*JVS( 311 )
  W( 81 ) = W( 81 ) + a*JVS( 312 )
  a = -W( 60 ) / JVS(          315  )
  W( 60 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 316 )
  W( 71 ) = W( 71 ) + a*JVS( 317 )
  W( 75 ) = W( 75 ) + a*JVS( 318 )
  W( 76 ) = W( 76 ) + a*JVS( 319 )
  W( 79 ) = W( 79 ) + a*JVS( 320 )
  W( 80 ) = W( 80 ) + a*JVS( 321 )
  W( 81 ) = W( 81 ) + a*JVS( 322 )
  a = -W( 61 ) / JVS(          324  )
  W( 61 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 325 )
  W( 75 ) = W( 75 ) + a*JVS( 326 )
  W( 80 ) = W( 80 ) + a*JVS( 327 )
  W( 81 ) = W( 81 ) + a*JVS( 328 )
  a = -W( 62 ) / JVS(          341  )
  W( 62 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  W( 76 ) = W( 76 ) + a*JVS( 344 )
  W( 77 ) = W( 77 ) + a*JVS( 345 )
  W( 79 ) = W( 79 ) + a*JVS( 346 )
  W( 80 ) = W( 80 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  a = -W( 63 ) / JVS(          353  )
  W( 63 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  W( 76 ) = W( 76 ) + a*JVS( 356 )
  W( 77 ) = W( 77 ) + a*JVS( 357 )
  W( 79 ) = W( 79 ) + a*JVS( 358 )
  a = -W( 64 ) / JVS(          361  )
  W( 64 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 362 )
  W( 76 ) = W( 76 ) + a*JVS( 363 )
  W( 79 ) = W( 79 ) + a*JVS( 364 )
  W( 80 ) = W( 80 ) + a*JVS( 365 )
  W( 81 ) = W( 81 ) + a*JVS( 366 )
  a = -W( 65 ) / JVS(          373  )
  W( 65 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 374 )
  W( 69 ) = W( 69 ) + a*JVS( 375 )
  W( 70 ) = W( 70 ) + a*JVS( 376 )
  W( 73 ) = W( 73 ) + a*JVS( 377 )
  W( 74 ) = W( 74 ) + a*JVS( 378 )
  W( 75 ) = W( 75 ) + a*JVS( 379 )
  W( 76 ) = W( 76 ) + a*JVS( 380 )
  W( 77 ) = W( 77 ) + a*JVS( 381 )
  W( 78 ) = W( 78 ) + a*JVS( 382 )
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 80 ) = W( 80 ) + a*JVS( 384 )
  W( 81 ) = W( 81 ) + a*JVS( 385 )
  a = -W( 66 ) / JVS(          392  )
  W( 66 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 393 )
  W( 75 ) = W( 75 ) + a*JVS( 394 )
  W( 76 ) = W( 76 ) + a*JVS( 395 )
  W( 77 ) = W( 77 ) + a*JVS( 396 )
  W( 78 ) = W( 78 ) + a*JVS( 397 )
  W( 79 ) = W( 79 ) + a*JVS( 398 )
  W( 80 ) = W( 80 ) + a*JVS( 399 )
  W( 81 ) = W( 81 ) + a*JVS( 400 )
  a = -W( 67 ) / JVS(          411  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 412 )
  W( 69 ) = W( 69 ) + a*JVS( 413 )
  W( 70 ) = W( 70 ) + a*JVS( 414 )
  W( 71 ) = W( 71 ) + a*JVS( 415 )
  W( 73 ) = W( 73 ) + a*JVS( 416 )
  W( 74 ) = W( 74 ) + a*JVS( 417 )
  W( 75 ) = W( 75 ) + a*JVS( 418 )
  W( 76 ) = W( 76 ) + a*JVS( 419 )
  W( 77 ) = W( 77 ) + a*JVS( 420 )
  W( 78 ) = W( 78 ) + a*JVS( 421 )
  W( 79 ) = W( 79 ) + a*JVS( 422 )
  W( 80 ) = W( 80 ) + a*JVS( 423 )
  W( 81 ) = W( 81 ) + a*JVS( 424 )
  a = -W( 68 ) / JVS(          429  )
  W( 68 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 430 )
  W( 71 ) = W( 71 ) + a*JVS( 431 )
  W( 75 ) = W( 75 ) + a*JVS( 432 )
  W( 76 ) = W( 76 ) + a*JVS( 433 )
  W( 77 ) = W( 77 ) + a*JVS( 434 )
  W( 78 ) = W( 78 ) + a*JVS( 435 )
  W( 79 ) = W( 79 ) + a*JVS( 436 )
  W( 80 ) = W( 80 ) + a*JVS( 437 )
  W( 81 ) = W( 81 ) + a*JVS( 438 )
  a = -W( 69 ) / JVS(          445  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 446 )
  W( 75 ) = W( 75 ) + a*JVS( 447 )
  W( 76 ) = W( 76 ) + a*JVS( 448 )
  W( 77 ) = W( 77 ) + a*JVS( 449 )
  W( 78 ) = W( 78 ) + a*JVS( 450 )
  W( 79 ) = W( 79 ) + a*JVS( 451 )
  W( 80 ) = W( 80 ) + a*JVS( 452 )
  W( 81 ) = W( 81 ) + a*JVS( 453 )
  a = -W( 70 ) / JVS(          457  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 75 ) = W( 75 ) + a*JVS( 459 )
  W( 76 ) = W( 76 ) + a*JVS( 460 )
  W( 77 ) = W( 77 ) + a*JVS( 461 )
  W( 78 ) = W( 78 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  W( 81 ) = W( 81 ) + a*JVS( 465 )
  a = -W( 71 ) / JVS(          468  )
  W( 71 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 469 )
  W( 76 ) = W( 76 ) + a*JVS( 470 )
  W( 77 ) = W( 77 ) + a*JVS( 471 )
  W( 78 ) = W( 78 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  W( 80 ) = W( 80 ) + a*JVS( 474 )
  W( 81 ) = W( 81 ) + a*JVS( 475 )
  a = -W( 72 ) / JVS(          508  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 509 )
  W( 74 ) = W( 74 ) + a*JVS( 510 )
  W( 75 ) = W( 75 ) + a*JVS( 511 )
  W( 76 ) = W( 76 ) + a*JVS( 512 )
  W( 77 ) = W( 77 ) + a*JVS( 513 )
  W( 78 ) = W( 78 ) + a*JVS( 514 )
  W( 79 ) = W( 79 ) + a*JVS( 515 )
  W( 80 ) = W( 80 ) + a*JVS( 516 )
  W( 81 ) = W( 81 ) + a*JVS( 517 )
  a = -W( 73 ) / JVS(          548  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 76 ) = W( 76 ) + a*JVS( 551 )
  W( 77 ) = W( 77 ) + a*JVS( 552 )
  W( 78 ) = W( 78 ) + a*JVS( 553 )
  W( 79 ) = W( 79 ) + a*JVS( 554 )
  W( 80 ) = W( 80 ) + a*JVS( 555 )
  W( 81 ) = W( 81 ) + a*JVS( 556 )
  a = -W( 74 ) / JVS(          564  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 565 )
  W( 76 ) = W( 76 ) + a*JVS( 566 )
  W( 77 ) = W( 77 ) + a*JVS( 567 )
  W( 78 ) = W( 78 ) + a*JVS( 568 )
  W( 79 ) = W( 79 ) + a*JVS( 569 )
  W( 80 ) = W( 80 ) + a*JVS( 570 )
  W( 81 ) = W( 81 ) + a*JVS( 571 )
  JVS( 572) = W( 2 )
  JVS( 573) = W( 4 )
  JVS( 574) = W( 5 )
  JVS( 575) = W( 6 )
  JVS( 576) = W( 7 )
  JVS( 577) = W( 8 )
  JVS( 578) = W( 9 )
  JVS( 579) = W( 10 )
  JVS( 580) = W( 11 )
  JVS( 581) = W( 15 )
  JVS( 582) = W( 16 )
  JVS( 583) = W( 17 )
  JVS( 584) = W( 18 )
  JVS( 585) = W( 19 )
  JVS( 586) = W( 20 )
  JVS( 587) = W( 21 )
  JVS( 588) = W( 22 )
  JVS( 589) = W( 23 )
  JVS( 590) = W( 25 )
  JVS( 591) = W( 26 )
  JVS( 592) = W( 27 )
  JVS( 593) = W( 28 )
  JVS( 594) = W( 29 )
  JVS( 595) = W( 31 )
  JVS( 596) = W( 32 )
  JVS( 597) = W( 33 )
  JVS( 598) = W( 34 )
  JVS( 599) = W( 35 )
  JVS( 600) = W( 36 )
  JVS( 601) = W( 37 )
  JVS( 602) = W( 38 )
  JVS( 603) = W( 39 )
  JVS( 604) = W( 40 )
  JVS( 605) = W( 41 )
  JVS( 606) = W( 42 )
  JVS( 607) = W( 43 )
  JVS( 608) = W( 44 )
  JVS( 609) = W( 45 )
  JVS( 610) = W( 46 )
  JVS( 611) = W( 47 )
  JVS( 612) = W( 48 )
  JVS( 613) = W( 49 )
  JVS( 614) = W( 50 )
  JVS( 615) = W( 51 )
  JVS( 616) = W( 52 )
  JVS( 617) = W( 53 )
  JVS( 618) = W( 54 )
  JVS( 619) = W( 55 )
  JVS( 620) = W( 56 )
  JVS( 621) = W( 57 )
  JVS( 622) = W( 58 )
  JVS( 623) = W( 59 )
  JVS( 624) = W( 60 )
  JVS( 625) = W( 61 )
  JVS( 626) = W( 62 )
  JVS( 627) = W( 63 )
  JVS( 628) = W( 64 )
  JVS( 629) = W( 65 )
  JVS( 630) = W( 66 )
  JVS( 631) = W( 67 )
  JVS( 632) = W( 68 )
  JVS( 633) = W( 69 )
  JVS( 634) = W( 70 )
  JVS( 635) = W( 71 )
  JVS( 636) = W( 72 )
  JVS( 637) = W( 73 )
  JVS( 638) = W( 74 )
  JVS( 639) = W( 75 )
  JVS( 640) = W( 76 )
  JVS( 641) = W( 77 )
  JVS( 642) = W( 78 )
  JVS( 643) = W( 79 )
  JVS( 644) = W( 80 )
  JVS( 645) = W( 81 )
  IF ( ABS(  JVS( 670 )) < TINY(a) ) THEN
         IER = 76                                      
         RETURN
  END IF
   W( 3 ) = JVS( 646 )
   W( 24 ) = JVS( 647 )
   W( 30 ) = JVS( 648 )
   W( 44 ) = JVS( 649 )
   W( 47 ) = JVS( 650 )
   W( 49 ) = JVS( 651 )
   W( 50 ) = JVS( 652 )
   W( 52 ) = JVS( 653 )
   W( 54 ) = JVS( 654 )
   W( 55 ) = JVS( 655 )
   W( 56 ) = JVS( 656 )
   W( 57 ) = JVS( 657 )
   W( 59 ) = JVS( 658 )
   W( 61 ) = JVS( 659 )
   W( 63 ) = JVS( 660 )
   W( 64 ) = JVS( 661 )
   W( 68 ) = JVS( 662 )
   W( 69 ) = JVS( 663 )
   W( 70 ) = JVS( 664 )
   W( 71 ) = JVS( 665 )
   W( 72 ) = JVS( 666 )
   W( 73 ) = JVS( 667 )
   W( 74 ) = JVS( 668 )
   W( 75 ) = JVS( 669 )
   W( 76 ) = JVS( 670 )
   W( 77 ) = JVS( 671 )
   W( 78 ) = JVS( 672 )
   W( 79 ) = JVS( 673 )
   W( 80 ) = JVS( 674 )
   W( 81 ) = JVS( 675 )
  a = -W( 3 ) / JVS(            6  )
  W( 3 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 7 )
  a = -W( 24 ) / JVS(           75  )
  W( 24 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 76 )
  W( 76 ) = W( 76 ) + a*JVS( 77 )
  a = -W( 30 ) / JVS(          104  )
  W( 30 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 105 )
  W( 76 ) = W( 76 ) + a*JVS( 106 )
  W( 81 ) = W( 81 ) + a*JVS( 107 )
  a = -W( 44 ) / JVS(          184  )
  W( 44 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 185 )
  W( 76 ) = W( 76 ) + a*JVS( 186 )
  W( 79 ) = W( 79 ) + a*JVS( 187 )
  a = -W( 47 ) / JVS(          203  )
  W( 47 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 204 )
  W( 72 ) = W( 72 ) + a*JVS( 205 )
  W( 75 ) = W( 75 ) + a*JVS( 206 )
  W( 81 ) = W( 81 ) + a*JVS( 207 )
  a = -W( 49 ) / JVS(          215  )
  W( 49 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 216 )
  W( 75 ) = W( 75 ) + a*JVS( 217 )
  W( 76 ) = W( 76 ) + a*JVS( 218 )
  W( 79 ) = W( 79 ) + a*JVS( 219 )
  a = -W( 50 ) / JVS(          221  )
  W( 50 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 222 )
  W( 75 ) = W( 75 ) + a*JVS( 223 )
  W( 76 ) = W( 76 ) + a*JVS( 224 )
  W( 79 ) = W( 79 ) + a*JVS( 225 )
  a = -W( 52 ) / JVS(          239  )
  W( 52 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 240 )
  W( 72 ) = W( 72 ) + a*JVS( 241 )
  W( 73 ) = W( 73 ) + a*JVS( 242 )
  W( 75 ) = W( 75 ) + a*JVS( 243 )
  W( 79 ) = W( 79 ) + a*JVS( 244 )
  W( 81 ) = W( 81 ) + a*JVS( 245 )
  a = -W( 54 ) / JVS(          259  )
  W( 54 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 260 )
  W( 76 ) = W( 76 ) + a*JVS( 261 )
  W( 79 ) = W( 79 ) + a*JVS( 262 )
  a = -W( 55 ) / JVS(          268  )
  W( 55 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 269 )
  W( 76 ) = W( 76 ) + a*JVS( 270 )
  W( 77 ) = W( 77 ) + a*JVS( 271 )
  W( 79 ) = W( 79 ) + a*JVS( 272 )
  a = -W( 56 ) / JVS(          273  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 274 )
  W( 76 ) = W( 76 ) + a*JVS( 275 )
  W( 79 ) = W( 79 ) + a*JVS( 276 )
  W( 80 ) = W( 80 ) + a*JVS( 277 )
  a = -W( 57 ) / JVS(          280  )
  W( 57 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 281 )
  W( 76 ) = W( 76 ) + a*JVS( 282 )
  W( 77 ) = W( 77 ) + a*JVS( 283 )
  W( 79 ) = W( 79 ) + a*JVS( 284 )
  a = -W( 59 ) / JVS(          309  )
  W( 59 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 310 )
  W( 80 ) = W( 80 ) + a*JVS( 311 )
  W( 81 ) = W( 81 ) + a*JVS( 312 )
  a = -W( 61 ) / JVS(          324  )
  W( 61 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 325 )
  W( 75 ) = W( 75 ) + a*JVS( 326 )
  W( 80 ) = W( 80 ) + a*JVS( 327 )
  W( 81 ) = W( 81 ) + a*JVS( 328 )
  a = -W( 63 ) / JVS(          353  )
  W( 63 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  W( 76 ) = W( 76 ) + a*JVS( 356 )
  W( 77 ) = W( 77 ) + a*JVS( 357 )
  W( 79 ) = W( 79 ) + a*JVS( 358 )
  a = -W( 64 ) / JVS(          361  )
  W( 64 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 362 )
  W( 76 ) = W( 76 ) + a*JVS( 363 )
  W( 79 ) = W( 79 ) + a*JVS( 364 )
  W( 80 ) = W( 80 ) + a*JVS( 365 )
  W( 81 ) = W( 81 ) + a*JVS( 366 )
  a = -W( 68 ) / JVS(          429  )
  W( 68 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 430 )
  W( 71 ) = W( 71 ) + a*JVS( 431 )
  W( 75 ) = W( 75 ) + a*JVS( 432 )
  W( 76 ) = W( 76 ) + a*JVS( 433 )
  W( 77 ) = W( 77 ) + a*JVS( 434 )
  W( 78 ) = W( 78 ) + a*JVS( 435 )
  W( 79 ) = W( 79 ) + a*JVS( 436 )
  W( 80 ) = W( 80 ) + a*JVS( 437 )
  W( 81 ) = W( 81 ) + a*JVS( 438 )
  a = -W( 69 ) / JVS(          445  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 446 )
  W( 75 ) = W( 75 ) + a*JVS( 447 )
  W( 76 ) = W( 76 ) + a*JVS( 448 )
  W( 77 ) = W( 77 ) + a*JVS( 449 )
  W( 78 ) = W( 78 ) + a*JVS( 450 )
  W( 79 ) = W( 79 ) + a*JVS( 451 )
  W( 80 ) = W( 80 ) + a*JVS( 452 )
  W( 81 ) = W( 81 ) + a*JVS( 453 )
  a = -W( 70 ) / JVS(          457  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 75 ) = W( 75 ) + a*JVS( 459 )
  W( 76 ) = W( 76 ) + a*JVS( 460 )
  W( 77 ) = W( 77 ) + a*JVS( 461 )
  W( 78 ) = W( 78 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  W( 81 ) = W( 81 ) + a*JVS( 465 )
  a = -W( 71 ) / JVS(          468  )
  W( 71 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 469 )
  W( 76 ) = W( 76 ) + a*JVS( 470 )
  W( 77 ) = W( 77 ) + a*JVS( 471 )
  W( 78 ) = W( 78 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  W( 80 ) = W( 80 ) + a*JVS( 474 )
  W( 81 ) = W( 81 ) + a*JVS( 475 )
  a = -W( 72 ) / JVS(          508  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 509 )
  W( 74 ) = W( 74 ) + a*JVS( 510 )
  W( 75 ) = W( 75 ) + a*JVS( 511 )
  W( 76 ) = W( 76 ) + a*JVS( 512 )
  W( 77 ) = W( 77 ) + a*JVS( 513 )
  W( 78 ) = W( 78 ) + a*JVS( 514 )
  W( 79 ) = W( 79 ) + a*JVS( 515 )
  W( 80 ) = W( 80 ) + a*JVS( 516 )
  W( 81 ) = W( 81 ) + a*JVS( 517 )
  a = -W( 73 ) / JVS(          548  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 76 ) = W( 76 ) + a*JVS( 551 )
  W( 77 ) = W( 77 ) + a*JVS( 552 )
  W( 78 ) = W( 78 ) + a*JVS( 553 )
  W( 79 ) = W( 79 ) + a*JVS( 554 )
  W( 80 ) = W( 80 ) + a*JVS( 555 )
  W( 81 ) = W( 81 ) + a*JVS( 556 )
  a = -W( 74 ) / JVS(          564  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 565 )
  W( 76 ) = W( 76 ) + a*JVS( 566 )
  W( 77 ) = W( 77 ) + a*JVS( 567 )
  W( 78 ) = W( 78 ) + a*JVS( 568 )
  W( 79 ) = W( 79 ) + a*JVS( 569 )
  W( 80 ) = W( 80 ) + a*JVS( 570 )
  W( 81 ) = W( 81 ) + a*JVS( 571 )
  a = -W( 75 ) / JVS(          639  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 640 )
  W( 77 ) = W( 77 ) + a*JVS( 641 )
  W( 78 ) = W( 78 ) + a*JVS( 642 )
  W( 79 ) = W( 79 ) + a*JVS( 643 )
  W( 80 ) = W( 80 ) + a*JVS( 644 )
  W( 81 ) = W( 81 ) + a*JVS( 645 )
  JVS( 646) = W( 3 )
  JVS( 647) = W( 24 )
  JVS( 648) = W( 30 )
  JVS( 649) = W( 44 )
  JVS( 650) = W( 47 )
  JVS( 651) = W( 49 )
  JVS( 652) = W( 50 )
  JVS( 653) = W( 52 )
  JVS( 654) = W( 54 )
  JVS( 655) = W( 55 )
  JVS( 656) = W( 56 )
  JVS( 657) = W( 57 )
  JVS( 658) = W( 59 )
  JVS( 659) = W( 61 )
  JVS( 660) = W( 63 )
  JVS( 661) = W( 64 )
  JVS( 662) = W( 68 )
  JVS( 663) = W( 69 )
  JVS( 664) = W( 70 )
  JVS( 665) = W( 71 )
  JVS( 666) = W( 72 )
  JVS( 667) = W( 73 )
  JVS( 668) = W( 74 )
  JVS( 669) = W( 75 )
  JVS( 670) = W( 76 )
  JVS( 671) = W( 77 )
  JVS( 672) = W( 78 )
  JVS( 673) = W( 79 )
  JVS( 674) = W( 80 )
  JVS( 675) = W( 81 )
  IF ( ABS(  JVS( 700 )) < TINY(a) ) THEN
         IER = 77                                      
         RETURN
  END IF
   W( 18 ) = JVS( 676 )
   W( 32 ) = JVS( 677 )
   W( 34 ) = JVS( 678 )
   W( 36 ) = JVS( 679 )
   W( 40 ) = JVS( 680 )
   W( 47 ) = JVS( 681 )
   W( 51 ) = JVS( 682 )
   W( 54 ) = JVS( 683 )
   W( 55 ) = JVS( 684 )
   W( 57 ) = JVS( 685 )
   W( 59 ) = JVS( 686 )
   W( 61 ) = JVS( 687 )
   W( 62 ) = JVS( 688 )
   W( 63 ) = JVS( 689 )
   W( 64 ) = JVS( 690 )
   W( 68 ) = JVS( 691 )
   W( 69 ) = JVS( 692 )
   W( 70 ) = JVS( 693 )
   W( 71 ) = JVS( 694 )
   W( 72 ) = JVS( 695 )
   W( 73 ) = JVS( 696 )
   W( 74 ) = JVS( 697 )
   W( 75 ) = JVS( 698 )
   W( 76 ) = JVS( 699 )
   W( 77 ) = JVS( 700 )
   W( 78 ) = JVS( 701 )
   W( 79 ) = JVS( 702 )
   W( 80 ) = JVS( 703 )
   W( 81 ) = JVS( 704 )
  a = -W( 18 ) / JVS(           49  )
  W( 18 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 50 )
  W( 77 ) = W( 77 ) + a*JVS( 51 )
  W( 79 ) = W( 79 ) + a*JVS( 52 )
  a = -W( 32 ) / JVS(          112  )
  W( 32 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 113 )
  W( 75 ) = W( 75 ) + a*JVS( 114 )
  W( 78 ) = W( 78 ) + a*JVS( 115 )
  W( 79 ) = W( 79 ) + a*JVS( 116 )
  a = -W( 34 ) / JVS(          123  )
  W( 34 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 124 )
  W( 61 ) = W( 61 ) + a*JVS( 125 )
  W( 74 ) = W( 74 ) + a*JVS( 126 )
  W( 75 ) = W( 75 ) + a*JVS( 127 )
  W( 77 ) = W( 77 ) + a*JVS( 128 )
  W( 78 ) = W( 78 ) + a*JVS( 129 )
  W( 79 ) = W( 79 ) + a*JVS( 130 )
  W( 81 ) = W( 81 ) + a*JVS( 131 )
  a = -W( 36 ) / JVS(          136  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 137 )
  W( 75 ) = W( 75 ) + a*JVS( 138 )
  W( 78 ) = W( 78 ) + a*JVS( 139 )
  a = -W( 40 ) / JVS(          161  )
  W( 40 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 162 )
  W( 61 ) = W( 61 ) + a*JVS( 163 )
  W( 75 ) = W( 75 ) + a*JVS( 164 )
  W( 81 ) = W( 81 ) + a*JVS( 165 )
  a = -W( 47 ) / JVS(          203  )
  W( 47 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 204 )
  W( 72 ) = W( 72 ) + a*JVS( 205 )
  W( 75 ) = W( 75 ) + a*JVS( 206 )
  W( 81 ) = W( 81 ) + a*JVS( 207 )
  a = -W( 51 ) / JVS(          230  )
  W( 51 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 231 )
  W( 57 ) = W( 57 ) + a*JVS( 232 )
  W( 64 ) = W( 64 ) + a*JVS( 233 )
  W( 75 ) = W( 75 ) + a*JVS( 234 )
  W( 76 ) = W( 76 ) + a*JVS( 235 )
  W( 77 ) = W( 77 ) + a*JVS( 236 )
  W( 79 ) = W( 79 ) + a*JVS( 237 )
  a = -W( 54 ) / JVS(          259  )
  W( 54 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 260 )
  W( 76 ) = W( 76 ) + a*JVS( 261 )
  W( 79 ) = W( 79 ) + a*JVS( 262 )
  a = -W( 55 ) / JVS(          268  )
  W( 55 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 269 )
  W( 76 ) = W( 76 ) + a*JVS( 270 )
  W( 77 ) = W( 77 ) + a*JVS( 271 )
  W( 79 ) = W( 79 ) + a*JVS( 272 )
  a = -W( 57 ) / JVS(          280  )
  W( 57 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 281 )
  W( 76 ) = W( 76 ) + a*JVS( 282 )
  W( 77 ) = W( 77 ) + a*JVS( 283 )
  W( 79 ) = W( 79 ) + a*JVS( 284 )
  a = -W( 59 ) / JVS(          309  )
  W( 59 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 310 )
  W( 80 ) = W( 80 ) + a*JVS( 311 )
  W( 81 ) = W( 81 ) + a*JVS( 312 )
  a = -W( 61 ) / JVS(          324  )
  W( 61 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 325 )
  W( 75 ) = W( 75 ) + a*JVS( 326 )
  W( 80 ) = W( 80 ) + a*JVS( 327 )
  W( 81 ) = W( 81 ) + a*JVS( 328 )
  a = -W( 62 ) / JVS(          341  )
  W( 62 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  W( 76 ) = W( 76 ) + a*JVS( 344 )
  W( 77 ) = W( 77 ) + a*JVS( 345 )
  W( 79 ) = W( 79 ) + a*JVS( 346 )
  W( 80 ) = W( 80 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  a = -W( 63 ) / JVS(          353  )
  W( 63 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  W( 76 ) = W( 76 ) + a*JVS( 356 )
  W( 77 ) = W( 77 ) + a*JVS( 357 )
  W( 79 ) = W( 79 ) + a*JVS( 358 )
  a = -W( 64 ) / JVS(          361  )
  W( 64 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 362 )
  W( 76 ) = W( 76 ) + a*JVS( 363 )
  W( 79 ) = W( 79 ) + a*JVS( 364 )
  W( 80 ) = W( 80 ) + a*JVS( 365 )
  W( 81 ) = W( 81 ) + a*JVS( 366 )
  a = -W( 68 ) / JVS(          429  )
  W( 68 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 430 )
  W( 71 ) = W( 71 ) + a*JVS( 431 )
  W( 75 ) = W( 75 ) + a*JVS( 432 )
  W( 76 ) = W( 76 ) + a*JVS( 433 )
  W( 77 ) = W( 77 ) + a*JVS( 434 )
  W( 78 ) = W( 78 ) + a*JVS( 435 )
  W( 79 ) = W( 79 ) + a*JVS( 436 )
  W( 80 ) = W( 80 ) + a*JVS( 437 )
  W( 81 ) = W( 81 ) + a*JVS( 438 )
  a = -W( 69 ) / JVS(          445  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 446 )
  W( 75 ) = W( 75 ) + a*JVS( 447 )
  W( 76 ) = W( 76 ) + a*JVS( 448 )
  W( 77 ) = W( 77 ) + a*JVS( 449 )
  W( 78 ) = W( 78 ) + a*JVS( 450 )
  W( 79 ) = W( 79 ) + a*JVS( 451 )
  W( 80 ) = W( 80 ) + a*JVS( 452 )
  W( 81 ) = W( 81 ) + a*JVS( 453 )
  a = -W( 70 ) / JVS(          457  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 75 ) = W( 75 ) + a*JVS( 459 )
  W( 76 ) = W( 76 ) + a*JVS( 460 )
  W( 77 ) = W( 77 ) + a*JVS( 461 )
  W( 78 ) = W( 78 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  W( 81 ) = W( 81 ) + a*JVS( 465 )
  a = -W( 71 ) / JVS(          468  )
  W( 71 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 469 )
  W( 76 ) = W( 76 ) + a*JVS( 470 )
  W( 77 ) = W( 77 ) + a*JVS( 471 )
  W( 78 ) = W( 78 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  W( 80 ) = W( 80 ) + a*JVS( 474 )
  W( 81 ) = W( 81 ) + a*JVS( 475 )
  a = -W( 72 ) / JVS(          508  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 509 )
  W( 74 ) = W( 74 ) + a*JVS( 510 )
  W( 75 ) = W( 75 ) + a*JVS( 511 )
  W( 76 ) = W( 76 ) + a*JVS( 512 )
  W( 77 ) = W( 77 ) + a*JVS( 513 )
  W( 78 ) = W( 78 ) + a*JVS( 514 )
  W( 79 ) = W( 79 ) + a*JVS( 515 )
  W( 80 ) = W( 80 ) + a*JVS( 516 )
  W( 81 ) = W( 81 ) + a*JVS( 517 )
  a = -W( 73 ) / JVS(          548  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 76 ) = W( 76 ) + a*JVS( 551 )
  W( 77 ) = W( 77 ) + a*JVS( 552 )
  W( 78 ) = W( 78 ) + a*JVS( 553 )
  W( 79 ) = W( 79 ) + a*JVS( 554 )
  W( 80 ) = W( 80 ) + a*JVS( 555 )
  W( 81 ) = W( 81 ) + a*JVS( 556 )
  a = -W( 74 ) / JVS(          564  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 565 )
  W( 76 ) = W( 76 ) + a*JVS( 566 )
  W( 77 ) = W( 77 ) + a*JVS( 567 )
  W( 78 ) = W( 78 ) + a*JVS( 568 )
  W( 79 ) = W( 79 ) + a*JVS( 569 )
  W( 80 ) = W( 80 ) + a*JVS( 570 )
  W( 81 ) = W( 81 ) + a*JVS( 571 )
  a = -W( 75 ) / JVS(          639  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 640 )
  W( 77 ) = W( 77 ) + a*JVS( 641 )
  W( 78 ) = W( 78 ) + a*JVS( 642 )
  W( 79 ) = W( 79 ) + a*JVS( 643 )
  W( 80 ) = W( 80 ) + a*JVS( 644 )
  W( 81 ) = W( 81 ) + a*JVS( 645 )
  a = -W( 76 ) / JVS(          670  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 671 )
  W( 78 ) = W( 78 ) + a*JVS( 672 )
  W( 79 ) = W( 79 ) + a*JVS( 673 )
  W( 80 ) = W( 80 ) + a*JVS( 674 )
  W( 81 ) = W( 81 ) + a*JVS( 675 )
  JVS( 676) = W( 18 )
  JVS( 677) = W( 32 )
  JVS( 678) = W( 34 )
  JVS( 679) = W( 36 )
  JVS( 680) = W( 40 )
  JVS( 681) = W( 47 )
  JVS( 682) = W( 51 )
  JVS( 683) = W( 54 )
  JVS( 684) = W( 55 )
  JVS( 685) = W( 57 )
  JVS( 686) = W( 59 )
  JVS( 687) = W( 61 )
  JVS( 688) = W( 62 )
  JVS( 689) = W( 63 )
  JVS( 690) = W( 64 )
  JVS( 691) = W( 68 )
  JVS( 692) = W( 69 )
  JVS( 693) = W( 70 )
  JVS( 694) = W( 71 )
  JVS( 695) = W( 72 )
  JVS( 696) = W( 73 )
  JVS( 697) = W( 74 )
  JVS( 698) = W( 75 )
  JVS( 699) = W( 76 )
  JVS( 700) = W( 77 )
  JVS( 701) = W( 78 )
  JVS( 702) = W( 79 )
  JVS( 703) = W( 80 )
  JVS( 704) = W( 81 )
  IF ( ABS(  JVS( 731 )) < TINY(a) ) THEN
         IER = 78                                      
         RETURN
  END IF
   W( 23 ) = JVS( 705 )
   W( 25 ) = JVS( 706 )
   W( 32 ) = JVS( 707 )
   W( 33 ) = JVS( 708 )
   W( 36 ) = JVS( 709 )
   W( 45 ) = JVS( 710 )
   W( 48 ) = JVS( 711 )
   W( 49 ) = JVS( 712 )
   W( 51 ) = JVS( 713 )
   W( 54 ) = JVS( 714 )
   W( 57 ) = JVS( 715 )
   W( 62 ) = JVS( 716 )
   W( 63 ) = JVS( 717 )
   W( 64 ) = JVS( 718 )
   W( 65 ) = JVS( 719 )
   W( 66 ) = JVS( 720 )
   W( 67 ) = JVS( 721 )
   W( 68 ) = JVS( 722 )
   W( 69 ) = JVS( 723 )
   W( 70 ) = JVS( 724 )
   W( 71 ) = JVS( 725 )
   W( 73 ) = JVS( 726 )
   W( 74 ) = JVS( 727 )
   W( 75 ) = JVS( 728 )
   W( 76 ) = JVS( 729 )
   W( 77 ) = JVS( 730 )
   W( 78 ) = JVS( 731 )
   W( 79 ) = JVS( 732 )
   W( 80 ) = JVS( 733 )
   W( 81 ) = JVS( 734 )
  a = -W( 23 ) / JVS(           70  )
  W( 23 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 71 )
  W( 75 ) = W( 75 ) + a*JVS( 72 )
  W( 79 ) = W( 79 ) + a*JVS( 73 )
  a = -W( 25 ) / JVS(           78  )
  W( 25 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 79 )
  W( 75 ) = W( 75 ) + a*JVS( 80 )
  W( 79 ) = W( 79 ) + a*JVS( 81 )
  a = -W( 32 ) / JVS(          112  )
  W( 32 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 113 )
  W( 75 ) = W( 75 ) + a*JVS( 114 )
  W( 78 ) = W( 78 ) + a*JVS( 115 )
  W( 79 ) = W( 79 ) + a*JVS( 116 )
  a = -W( 33 ) / JVS(          117  )
  W( 33 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 118 )
  W( 54 ) = W( 54 ) + a*JVS( 119 )
  W( 75 ) = W( 75 ) + a*JVS( 120 )
  W( 76 ) = W( 76 ) + a*JVS( 121 )
  a = -W( 36 ) / JVS(          136  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 137 )
  W( 75 ) = W( 75 ) + a*JVS( 138 )
  W( 78 ) = W( 78 ) + a*JVS( 139 )
  a = -W( 45 ) / JVS(          191  )
  W( 45 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 192 )
  W( 75 ) = W( 75 ) + a*JVS( 193 )
  W( 76 ) = W( 76 ) + a*JVS( 194 )
  W( 79 ) = W( 79 ) + a*JVS( 195 )
  a = -W( 48 ) / JVS(          208  )
  W( 48 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 209 )
  W( 75 ) = W( 75 ) + a*JVS( 210 )
  W( 79 ) = W( 79 ) + a*JVS( 211 )
  a = -W( 49 ) / JVS(          215  )
  W( 49 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 216 )
  W( 75 ) = W( 75 ) + a*JVS( 217 )
  W( 76 ) = W( 76 ) + a*JVS( 218 )
  W( 79 ) = W( 79 ) + a*JVS( 219 )
  a = -W( 51 ) / JVS(          230  )
  W( 51 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 231 )
  W( 57 ) = W( 57 ) + a*JVS( 232 )
  W( 64 ) = W( 64 ) + a*JVS( 233 )
  W( 75 ) = W( 75 ) + a*JVS( 234 )
  W( 76 ) = W( 76 ) + a*JVS( 235 )
  W( 77 ) = W( 77 ) + a*JVS( 236 )
  W( 79 ) = W( 79 ) + a*JVS( 237 )
  a = -W( 54 ) / JVS(          259  )
  W( 54 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 260 )
  W( 76 ) = W( 76 ) + a*JVS( 261 )
  W( 79 ) = W( 79 ) + a*JVS( 262 )
  a = -W( 57 ) / JVS(          280  )
  W( 57 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 281 )
  W( 76 ) = W( 76 ) + a*JVS( 282 )
  W( 77 ) = W( 77 ) + a*JVS( 283 )
  W( 79 ) = W( 79 ) + a*JVS( 284 )
  a = -W( 62 ) / JVS(          341  )
  W( 62 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  W( 76 ) = W( 76 ) + a*JVS( 344 )
  W( 77 ) = W( 77 ) + a*JVS( 345 )
  W( 79 ) = W( 79 ) + a*JVS( 346 )
  W( 80 ) = W( 80 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  a = -W( 63 ) / JVS(          353  )
  W( 63 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  W( 76 ) = W( 76 ) + a*JVS( 356 )
  W( 77 ) = W( 77 ) + a*JVS( 357 )
  W( 79 ) = W( 79 ) + a*JVS( 358 )
  a = -W( 64 ) / JVS(          361  )
  W( 64 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 362 )
  W( 76 ) = W( 76 ) + a*JVS( 363 )
  W( 79 ) = W( 79 ) + a*JVS( 364 )
  W( 80 ) = W( 80 ) + a*JVS( 365 )
  W( 81 ) = W( 81 ) + a*JVS( 366 )
  a = -W( 65 ) / JVS(          373  )
  W( 65 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 374 )
  W( 69 ) = W( 69 ) + a*JVS( 375 )
  W( 70 ) = W( 70 ) + a*JVS( 376 )
  W( 73 ) = W( 73 ) + a*JVS( 377 )
  W( 74 ) = W( 74 ) + a*JVS( 378 )
  W( 75 ) = W( 75 ) + a*JVS( 379 )
  W( 76 ) = W( 76 ) + a*JVS( 380 )
  W( 77 ) = W( 77 ) + a*JVS( 381 )
  W( 78 ) = W( 78 ) + a*JVS( 382 )
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 80 ) = W( 80 ) + a*JVS( 384 )
  W( 81 ) = W( 81 ) + a*JVS( 385 )
  a = -W( 66 ) / JVS(          392  )
  W( 66 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 393 )
  W( 75 ) = W( 75 ) + a*JVS( 394 )
  W( 76 ) = W( 76 ) + a*JVS( 395 )
  W( 77 ) = W( 77 ) + a*JVS( 396 )
  W( 78 ) = W( 78 ) + a*JVS( 397 )
  W( 79 ) = W( 79 ) + a*JVS( 398 )
  W( 80 ) = W( 80 ) + a*JVS( 399 )
  W( 81 ) = W( 81 ) + a*JVS( 400 )
  a = -W( 67 ) / JVS(          411  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 412 )
  W( 69 ) = W( 69 ) + a*JVS( 413 )
  W( 70 ) = W( 70 ) + a*JVS( 414 )
  W( 71 ) = W( 71 ) + a*JVS( 415 )
  W( 73 ) = W( 73 ) + a*JVS( 416 )
  W( 74 ) = W( 74 ) + a*JVS( 417 )
  W( 75 ) = W( 75 ) + a*JVS( 418 )
  W( 76 ) = W( 76 ) + a*JVS( 419 )
  W( 77 ) = W( 77 ) + a*JVS( 420 )
  W( 78 ) = W( 78 ) + a*JVS( 421 )
  W( 79 ) = W( 79 ) + a*JVS( 422 )
  W( 80 ) = W( 80 ) + a*JVS( 423 )
  W( 81 ) = W( 81 ) + a*JVS( 424 )
  a = -W( 68 ) / JVS(          429  )
  W( 68 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 430 )
  W( 71 ) = W( 71 ) + a*JVS( 431 )
  W( 75 ) = W( 75 ) + a*JVS( 432 )
  W( 76 ) = W( 76 ) + a*JVS( 433 )
  W( 77 ) = W( 77 ) + a*JVS( 434 )
  W( 78 ) = W( 78 ) + a*JVS( 435 )
  W( 79 ) = W( 79 ) + a*JVS( 436 )
  W( 80 ) = W( 80 ) + a*JVS( 437 )
  W( 81 ) = W( 81 ) + a*JVS( 438 )
  a = -W( 69 ) / JVS(          445  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 446 )
  W( 75 ) = W( 75 ) + a*JVS( 447 )
  W( 76 ) = W( 76 ) + a*JVS( 448 )
  W( 77 ) = W( 77 ) + a*JVS( 449 )
  W( 78 ) = W( 78 ) + a*JVS( 450 )
  W( 79 ) = W( 79 ) + a*JVS( 451 )
  W( 80 ) = W( 80 ) + a*JVS( 452 )
  W( 81 ) = W( 81 ) + a*JVS( 453 )
  a = -W( 70 ) / JVS(          457  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 75 ) = W( 75 ) + a*JVS( 459 )
  W( 76 ) = W( 76 ) + a*JVS( 460 )
  W( 77 ) = W( 77 ) + a*JVS( 461 )
  W( 78 ) = W( 78 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  W( 81 ) = W( 81 ) + a*JVS( 465 )
  a = -W( 71 ) / JVS(          468  )
  W( 71 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 469 )
  W( 76 ) = W( 76 ) + a*JVS( 470 )
  W( 77 ) = W( 77 ) + a*JVS( 471 )
  W( 78 ) = W( 78 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  W( 80 ) = W( 80 ) + a*JVS( 474 )
  W( 81 ) = W( 81 ) + a*JVS( 475 )
  a = -W( 73 ) / JVS(          548  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 76 ) = W( 76 ) + a*JVS( 551 )
  W( 77 ) = W( 77 ) + a*JVS( 552 )
  W( 78 ) = W( 78 ) + a*JVS( 553 )
  W( 79 ) = W( 79 ) + a*JVS( 554 )
  W( 80 ) = W( 80 ) + a*JVS( 555 )
  W( 81 ) = W( 81 ) + a*JVS( 556 )
  a = -W( 74 ) / JVS(          564  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 565 )
  W( 76 ) = W( 76 ) + a*JVS( 566 )
  W( 77 ) = W( 77 ) + a*JVS( 567 )
  W( 78 ) = W( 78 ) + a*JVS( 568 )
  W( 79 ) = W( 79 ) + a*JVS( 569 )
  W( 80 ) = W( 80 ) + a*JVS( 570 )
  W( 81 ) = W( 81 ) + a*JVS( 571 )
  a = -W( 75 ) / JVS(          639  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 640 )
  W( 77 ) = W( 77 ) + a*JVS( 641 )
  W( 78 ) = W( 78 ) + a*JVS( 642 )
  W( 79 ) = W( 79 ) + a*JVS( 643 )
  W( 80 ) = W( 80 ) + a*JVS( 644 )
  W( 81 ) = W( 81 ) + a*JVS( 645 )
  a = -W( 76 ) / JVS(          670  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 671 )
  W( 78 ) = W( 78 ) + a*JVS( 672 )
  W( 79 ) = W( 79 ) + a*JVS( 673 )
  W( 80 ) = W( 80 ) + a*JVS( 674 )
  W( 81 ) = W( 81 ) + a*JVS( 675 )
  a = -W( 77 ) / JVS(          700  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 701 )
  W( 79 ) = W( 79 ) + a*JVS( 702 )
  W( 80 ) = W( 80 ) + a*JVS( 703 )
  W( 81 ) = W( 81 ) + a*JVS( 704 )
  JVS( 705) = W( 23 )
  JVS( 706) = W( 25 )
  JVS( 707) = W( 32 )
  JVS( 708) = W( 33 )
  JVS( 709) = W( 36 )
  JVS( 710) = W( 45 )
  JVS( 711) = W( 48 )
  JVS( 712) = W( 49 )
  JVS( 713) = W( 51 )
  JVS( 714) = W( 54 )
  JVS( 715) = W( 57 )
  JVS( 716) = W( 62 )
  JVS( 717) = W( 63 )
  JVS( 718) = W( 64 )
  JVS( 719) = W( 65 )
  JVS( 720) = W( 66 )
  JVS( 721) = W( 67 )
  JVS( 722) = W( 68 )
  JVS( 723) = W( 69 )
  JVS( 724) = W( 70 )
  JVS( 725) = W( 71 )
  JVS( 726) = W( 73 )
  JVS( 727) = W( 74 )
  JVS( 728) = W( 75 )
  JVS( 729) = W( 76 )
  JVS( 730) = W( 77 )
  JVS( 731) = W( 78 )
  JVS( 732) = W( 79 )
  JVS( 733) = W( 80 )
  JVS( 734) = W( 81 )
  IF ( ABS(  JVS( 791 )) < TINY(a) ) THEN
         IER = 79                                      
         RETURN
  END IF
   W( 9 ) = JVS( 735 )
   W( 10 ) = JVS( 736 )
   W( 12 ) = JVS( 737 )
   W( 14 ) = JVS( 738 )
   W( 16 ) = JVS( 739 )
   W( 17 ) = JVS( 740 )
   W( 18 ) = JVS( 741 )
   W( 19 ) = JVS( 742 )
   W( 20 ) = JVS( 743 )
   W( 21 ) = JVS( 744 )
   W( 22 ) = JVS( 745 )
   W( 24 ) = JVS( 746 )
   W( 26 ) = JVS( 747 )
   W( 27 ) = JVS( 748 )
   W( 30 ) = JVS( 749 )
   W( 35 ) = JVS( 750 )
   W( 37 ) = JVS( 751 )
   W( 38 ) = JVS( 752 )
   W( 40 ) = JVS( 753 )
   W( 41 ) = JVS( 754 )
   W( 42 ) = JVS( 755 )
   W( 43 ) = JVS( 756 )
   W( 44 ) = JVS( 757 )
   W( 45 ) = JVS( 758 )
   W( 46 ) = JVS( 759 )
   W( 47 ) = JVS( 760 )
   W( 48 ) = JVS( 761 )
   W( 49 ) = JVS( 762 )
   W( 50 ) = JVS( 763 )
   W( 52 ) = JVS( 764 )
   W( 53 ) = JVS( 765 )
   W( 54 ) = JVS( 766 )
   W( 55 ) = JVS( 767 )
   W( 56 ) = JVS( 768 )
   W( 57 ) = JVS( 769 )
   W( 58 ) = JVS( 770 )
   W( 59 ) = JVS( 771 )
   W( 60 ) = JVS( 772 )
   W( 61 ) = JVS( 773 )
   W( 62 ) = JVS( 774 )
   W( 63 ) = JVS( 775 )
   W( 64 ) = JVS( 776 )
   W( 65 ) = JVS( 777 )
   W( 66 ) = JVS( 778 )
   W( 67 ) = JVS( 779 )
   W( 68 ) = JVS( 780 )
   W( 69 ) = JVS( 781 )
   W( 70 ) = JVS( 782 )
   W( 71 ) = JVS( 783 )
   W( 72 ) = JVS( 784 )
   W( 73 ) = JVS( 785 )
   W( 74 ) = JVS( 786 )
   W( 75 ) = JVS( 787 )
   W( 76 ) = JVS( 788 )
   W( 77 ) = JVS( 789 )
   W( 78 ) = JVS( 790 )
   W( 79 ) = JVS( 791 )
   W( 80 ) = JVS( 792 )
   W( 81 ) = JVS( 793 )
  a = -W( 9 ) / JVS(           20  )
  W( 9 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 21 )
  W( 79 ) = W( 79 ) + a*JVS( 22 )
  a = -W( 10 ) / JVS(           23  )
  W( 10 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 24 )
  a = -W( 12 ) / JVS(           28  )
  W( 12 ) = -a
  W( 30 ) = W( 30 ) + a*JVS( 29 )
  W( 76 ) = W( 76 ) + a*JVS( 30 )
  a = -W( 14 ) / JVS(           35  )
  W( 14 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 36 )
  W( 75 ) = W( 75 ) + a*JVS( 37 )
  a = -W( 16 ) / JVS(           42  )
  W( 16 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 43 )
  W( 80 ) = W( 80 ) + a*JVS( 44 )
  a = -W( 17 ) / JVS(           45  )
  W( 17 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 46 )
  W( 75 ) = W( 75 ) + a*JVS( 47 )
  W( 77 ) = W( 77 ) + a*JVS( 48 )
  a = -W( 18 ) / JVS(           49  )
  W( 18 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 50 )
  W( 77 ) = W( 77 ) + a*JVS( 51 )
  W( 79 ) = W( 79 ) + a*JVS( 52 )
  a = -W( 19 ) / JVS(           53  )
  W( 19 ) = -a
  W( 40 ) = W( 40 ) + a*JVS( 54 )
  W( 47 ) = W( 47 ) + a*JVS( 55 )
  W( 72 ) = W( 72 ) + a*JVS( 56 )
  W( 75 ) = W( 75 ) + a*JVS( 57 )
  a = -W( 20 ) / JVS(           58  )
  W( 20 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 59 )
  W( 75 ) = W( 75 ) + a*JVS( 60 )
  W( 79 ) = W( 79 ) + a*JVS( 61 )
  a = -W( 21 ) / JVS(           62  )
  W( 21 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 63 )
  W( 75 ) = W( 75 ) + a*JVS( 64 )
  W( 79 ) = W( 79 ) + a*JVS( 65 )
  a = -W( 22 ) / JVS(           66  )
  W( 22 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 67 )
  W( 75 ) = W( 75 ) + a*JVS( 68 )
  W( 79 ) = W( 79 ) + a*JVS( 69 )
  a = -W( 24 ) / JVS(           75  )
  W( 24 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 76 )
  W( 76 ) = W( 76 ) + a*JVS( 77 )
  a = -W( 26 ) / JVS(           82  )
  W( 26 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 83 )
  W( 75 ) = W( 75 ) + a*JVS( 84 )
  W( 79 ) = W( 79 ) + a*JVS( 85 )
  a = -W( 27 ) / JVS(           86  )
  W( 27 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 87 )
  W( 81 ) = W( 81 ) + a*JVS( 88 )
  a = -W( 30 ) / JVS(          104  )
  W( 30 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 105 )
  W( 76 ) = W( 76 ) + a*JVS( 106 )
  W( 81 ) = W( 81 ) + a*JVS( 107 )
  a = -W( 35 ) / JVS(          132  )
  W( 35 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 133 )
  W( 75 ) = W( 75 ) + a*JVS( 134 )
  W( 79 ) = W( 79 ) + a*JVS( 135 )
  a = -W( 37 ) / JVS(          141  )
  W( 37 ) = -a
  W( 44 ) = W( 44 ) + a*JVS( 142 )
  W( 45 ) = W( 45 ) + a*JVS( 143 )
  W( 53 ) = W( 53 ) + a*JVS( 144 )
  W( 75 ) = W( 75 ) + a*JVS( 145 )
  W( 76 ) = W( 76 ) + a*JVS( 146 )
  W( 79 ) = W( 79 ) + a*JVS( 147 )
  a = -W( 38 ) / JVS(          148  )
  W( 38 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 149 )
  W( 75 ) = W( 75 ) + a*JVS( 150 )
  W( 79 ) = W( 79 ) + a*JVS( 151 )
  a = -W( 40 ) / JVS(          161  )
  W( 40 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 162 )
  W( 61 ) = W( 61 ) + a*JVS( 163 )
  W( 75 ) = W( 75 ) + a*JVS( 164 )
  W( 81 ) = W( 81 ) + a*JVS( 165 )
  a = -W( 41 ) / JVS(          166  )
  W( 41 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 167 )
  W( 63 ) = W( 63 ) + a*JVS( 168 )
  W( 68 ) = W( 68 ) + a*JVS( 169 )
  W( 70 ) = W( 70 ) + a*JVS( 170 )
  W( 71 ) = W( 71 ) + a*JVS( 171 )
  W( 75 ) = W( 75 ) + a*JVS( 172 )
  W( 77 ) = W( 77 ) + a*JVS( 173 )
  a = -W( 42 ) / JVS(          174  )
  W( 42 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 175 )
  W( 75 ) = W( 75 ) + a*JVS( 176 )
  W( 79 ) = W( 79 ) + a*JVS( 177 )
  a = -W( 43 ) / JVS(          178  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 179 )
  W( 74 ) = W( 74 ) + a*JVS( 180 )
  W( 75 ) = W( 75 ) + a*JVS( 181 )
  a = -W( 44 ) / JVS(          184  )
  W( 44 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 185 )
  W( 76 ) = W( 76 ) + a*JVS( 186 )
  W( 79 ) = W( 79 ) + a*JVS( 187 )
  a = -W( 45 ) / JVS(          191  )
  W( 45 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 192 )
  W( 75 ) = W( 75 ) + a*JVS( 193 )
  W( 76 ) = W( 76 ) + a*JVS( 194 )
  W( 79 ) = W( 79 ) + a*JVS( 195 )
  a = -W( 46 ) / JVS(          196  )
  W( 46 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 197 )
  W( 80 ) = W( 80 ) + a*JVS( 198 )
  W( 81 ) = W( 81 ) + a*JVS( 199 )
  a = -W( 47 ) / JVS(          203  )
  W( 47 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 204 )
  W( 72 ) = W( 72 ) + a*JVS( 205 )
  W( 75 ) = W( 75 ) + a*JVS( 206 )
  W( 81 ) = W( 81 ) + a*JVS( 207 )
  a = -W( 48 ) / JVS(          208  )
  W( 48 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 209 )
  W( 75 ) = W( 75 ) + a*JVS( 210 )
  W( 79 ) = W( 79 ) + a*JVS( 211 )
  a = -W( 49 ) / JVS(          215  )
  W( 49 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 216 )
  W( 75 ) = W( 75 ) + a*JVS( 217 )
  W( 76 ) = W( 76 ) + a*JVS( 218 )
  W( 79 ) = W( 79 ) + a*JVS( 219 )
  a = -W( 50 ) / JVS(          221  )
  W( 50 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 222 )
  W( 75 ) = W( 75 ) + a*JVS( 223 )
  W( 76 ) = W( 76 ) + a*JVS( 224 )
  W( 79 ) = W( 79 ) + a*JVS( 225 )
  a = -W( 52 ) / JVS(          239  )
  W( 52 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 240 )
  W( 72 ) = W( 72 ) + a*JVS( 241 )
  W( 73 ) = W( 73 ) + a*JVS( 242 )
  W( 75 ) = W( 75 ) + a*JVS( 243 )
  W( 79 ) = W( 79 ) + a*JVS( 244 )
  W( 81 ) = W( 81 ) + a*JVS( 245 )
  a = -W( 53 ) / JVS(          248  )
  W( 53 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 249 )
  W( 70 ) = W( 70 ) + a*JVS( 250 )
  W( 75 ) = W( 75 ) + a*JVS( 251 )
  W( 76 ) = W( 76 ) + a*JVS( 252 )
  W( 77 ) = W( 77 ) + a*JVS( 253 )
  W( 78 ) = W( 78 ) + a*JVS( 254 )
  W( 80 ) = W( 80 ) + a*JVS( 255 )
  W( 81 ) = W( 81 ) + a*JVS( 256 )
  a = -W( 54 ) / JVS(          259  )
  W( 54 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 260 )
  W( 76 ) = W( 76 ) + a*JVS( 261 )
  W( 79 ) = W( 79 ) + a*JVS( 262 )
  a = -W( 55 ) / JVS(          268  )
  W( 55 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 269 )
  W( 76 ) = W( 76 ) + a*JVS( 270 )
  W( 77 ) = W( 77 ) + a*JVS( 271 )
  W( 79 ) = W( 79 ) + a*JVS( 272 )
  a = -W( 56 ) / JVS(          273  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 274 )
  W( 76 ) = W( 76 ) + a*JVS( 275 )
  W( 79 ) = W( 79 ) + a*JVS( 276 )
  W( 80 ) = W( 80 ) + a*JVS( 277 )
  a = -W( 57 ) / JVS(          280  )
  W( 57 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 281 )
  W( 76 ) = W( 76 ) + a*JVS( 282 )
  W( 77 ) = W( 77 ) + a*JVS( 283 )
  W( 79 ) = W( 79 ) + a*JVS( 284 )
  a = -W( 58 ) / JVS(          290  )
  W( 58 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 291 )
  W( 60 ) = W( 60 ) + a*JVS( 292 )
  W( 61 ) = W( 61 ) + a*JVS( 293 )
  W( 62 ) = W( 62 ) + a*JVS( 294 )
  W( 66 ) = W( 66 ) + a*JVS( 295 )
  W( 67 ) = W( 67 ) + a*JVS( 296 )
  W( 68 ) = W( 68 ) + a*JVS( 297 )
  W( 69 ) = W( 69 ) + a*JVS( 298 )
  W( 70 ) = W( 70 ) + a*JVS( 299 )
  W( 72 ) = W( 72 ) + a*JVS( 300 )
  W( 73 ) = W( 73 ) + a*JVS( 301 )
  W( 75 ) = W( 75 ) + a*JVS( 302 )
  W( 76 ) = W( 76 ) + a*JVS( 303 )
  W( 77 ) = W( 77 ) + a*JVS( 304 )
  W( 78 ) = W( 78 ) + a*JVS( 305 )
  W( 79 ) = W( 79 ) + a*JVS( 306 )
  W( 80 ) = W( 80 ) + a*JVS( 307 )
  W( 81 ) = W( 81 ) + a*JVS( 308 )
  a = -W( 59 ) / JVS(          309  )
  W( 59 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 310 )
  W( 80 ) = W( 80 ) + a*JVS( 311 )
  W( 81 ) = W( 81 ) + a*JVS( 312 )
  a = -W( 60 ) / JVS(          315  )
  W( 60 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 316 )
  W( 71 ) = W( 71 ) + a*JVS( 317 )
  W( 75 ) = W( 75 ) + a*JVS( 318 )
  W( 76 ) = W( 76 ) + a*JVS( 319 )
  W( 79 ) = W( 79 ) + a*JVS( 320 )
  W( 80 ) = W( 80 ) + a*JVS( 321 )
  W( 81 ) = W( 81 ) + a*JVS( 322 )
  a = -W( 61 ) / JVS(          324  )
  W( 61 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 325 )
  W( 75 ) = W( 75 ) + a*JVS( 326 )
  W( 80 ) = W( 80 ) + a*JVS( 327 )
  W( 81 ) = W( 81 ) + a*JVS( 328 )
  a = -W( 62 ) / JVS(          341  )
  W( 62 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  W( 76 ) = W( 76 ) + a*JVS( 344 )
  W( 77 ) = W( 77 ) + a*JVS( 345 )
  W( 79 ) = W( 79 ) + a*JVS( 346 )
  W( 80 ) = W( 80 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  a = -W( 63 ) / JVS(          353  )
  W( 63 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 354 )
  W( 75 ) = W( 75 ) + a*JVS( 355 )
  W( 76 ) = W( 76 ) + a*JVS( 356 )
  W( 77 ) = W( 77 ) + a*JVS( 357 )
  W( 79 ) = W( 79 ) + a*JVS( 358 )
  a = -W( 64 ) / JVS(          361  )
  W( 64 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 362 )
  W( 76 ) = W( 76 ) + a*JVS( 363 )
  W( 79 ) = W( 79 ) + a*JVS( 364 )
  W( 80 ) = W( 80 ) + a*JVS( 365 )
  W( 81 ) = W( 81 ) + a*JVS( 366 )
  a = -W( 65 ) / JVS(          373  )
  W( 65 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 374 )
  W( 69 ) = W( 69 ) + a*JVS( 375 )
  W( 70 ) = W( 70 ) + a*JVS( 376 )
  W( 73 ) = W( 73 ) + a*JVS( 377 )
  W( 74 ) = W( 74 ) + a*JVS( 378 )
  W( 75 ) = W( 75 ) + a*JVS( 379 )
  W( 76 ) = W( 76 ) + a*JVS( 380 )
  W( 77 ) = W( 77 ) + a*JVS( 381 )
  W( 78 ) = W( 78 ) + a*JVS( 382 )
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 80 ) = W( 80 ) + a*JVS( 384 )
  W( 81 ) = W( 81 ) + a*JVS( 385 )
  a = -W( 66 ) / JVS(          392  )
  W( 66 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 393 )
  W( 75 ) = W( 75 ) + a*JVS( 394 )
  W( 76 ) = W( 76 ) + a*JVS( 395 )
  W( 77 ) = W( 77 ) + a*JVS( 396 )
  W( 78 ) = W( 78 ) + a*JVS( 397 )
  W( 79 ) = W( 79 ) + a*JVS( 398 )
  W( 80 ) = W( 80 ) + a*JVS( 399 )
  W( 81 ) = W( 81 ) + a*JVS( 400 )
  a = -W( 67 ) / JVS(          411  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 412 )
  W( 69 ) = W( 69 ) + a*JVS( 413 )
  W( 70 ) = W( 70 ) + a*JVS( 414 )
  W( 71 ) = W( 71 ) + a*JVS( 415 )
  W( 73 ) = W( 73 ) + a*JVS( 416 )
  W( 74 ) = W( 74 ) + a*JVS( 417 )
  W( 75 ) = W( 75 ) + a*JVS( 418 )
  W( 76 ) = W( 76 ) + a*JVS( 419 )
  W( 77 ) = W( 77 ) + a*JVS( 420 )
  W( 78 ) = W( 78 ) + a*JVS( 421 )
  W( 79 ) = W( 79 ) + a*JVS( 422 )
  W( 80 ) = W( 80 ) + a*JVS( 423 )
  W( 81 ) = W( 81 ) + a*JVS( 424 )
  a = -W( 68 ) / JVS(          429  )
  W( 68 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 430 )
  W( 71 ) = W( 71 ) + a*JVS( 431 )
  W( 75 ) = W( 75 ) + a*JVS( 432 )
  W( 76 ) = W( 76 ) + a*JVS( 433 )
  W( 77 ) = W( 77 ) + a*JVS( 434 )
  W( 78 ) = W( 78 ) + a*JVS( 435 )
  W( 79 ) = W( 79 ) + a*JVS( 436 )
  W( 80 ) = W( 80 ) + a*JVS( 437 )
  W( 81 ) = W( 81 ) + a*JVS( 438 )
  a = -W( 69 ) / JVS(          445  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 446 )
  W( 75 ) = W( 75 ) + a*JVS( 447 )
  W( 76 ) = W( 76 ) + a*JVS( 448 )
  W( 77 ) = W( 77 ) + a*JVS( 449 )
  W( 78 ) = W( 78 ) + a*JVS( 450 )
  W( 79 ) = W( 79 ) + a*JVS( 451 )
  W( 80 ) = W( 80 ) + a*JVS( 452 )
  W( 81 ) = W( 81 ) + a*JVS( 453 )
  a = -W( 70 ) / JVS(          457  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 75 ) = W( 75 ) + a*JVS( 459 )
  W( 76 ) = W( 76 ) + a*JVS( 460 )
  W( 77 ) = W( 77 ) + a*JVS( 461 )
  W( 78 ) = W( 78 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  W( 81 ) = W( 81 ) + a*JVS( 465 )
  a = -W( 71 ) / JVS(          468  )
  W( 71 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 469 )
  W( 76 ) = W( 76 ) + a*JVS( 470 )
  W( 77 ) = W( 77 ) + a*JVS( 471 )
  W( 78 ) = W( 78 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  W( 80 ) = W( 80 ) + a*JVS( 474 )
  W( 81 ) = W( 81 ) + a*JVS( 475 )
  a = -W( 72 ) / JVS(          508  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 509 )
  W( 74 ) = W( 74 ) + a*JVS( 510 )
  W( 75 ) = W( 75 ) + a*JVS( 511 )
  W( 76 ) = W( 76 ) + a*JVS( 512 )
  W( 77 ) = W( 77 ) + a*JVS( 513 )
  W( 78 ) = W( 78 ) + a*JVS( 514 )
  W( 79 ) = W( 79 ) + a*JVS( 515 )
  W( 80 ) = W( 80 ) + a*JVS( 516 )
  W( 81 ) = W( 81 ) + a*JVS( 517 )
  a = -W( 73 ) / JVS(          548  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 76 ) = W( 76 ) + a*JVS( 551 )
  W( 77 ) = W( 77 ) + a*JVS( 552 )
  W( 78 ) = W( 78 ) + a*JVS( 553 )
  W( 79 ) = W( 79 ) + a*JVS( 554 )
  W( 80 ) = W( 80 ) + a*JVS( 555 )
  W( 81 ) = W( 81 ) + a*JVS( 556 )
  a = -W( 74 ) / JVS(          564  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 565 )
  W( 76 ) = W( 76 ) + a*JVS( 566 )
  W( 77 ) = W( 77 ) + a*JVS( 567 )
  W( 78 ) = W( 78 ) + a*JVS( 568 )
  W( 79 ) = W( 79 ) + a*JVS( 569 )
  W( 80 ) = W( 80 ) + a*JVS( 570 )
  W( 81 ) = W( 81 ) + a*JVS( 571 )
  a = -W( 75 ) / JVS(          639  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 640 )
  W( 77 ) = W( 77 ) + a*JVS( 641 )
  W( 78 ) = W( 78 ) + a*JVS( 642 )
  W( 79 ) = W( 79 ) + a*JVS( 643 )
  W( 80 ) = W( 80 ) + a*JVS( 644 )
  W( 81 ) = W( 81 ) + a*JVS( 645 )
  a = -W( 76 ) / JVS(          670  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 671 )
  W( 78 ) = W( 78 ) + a*JVS( 672 )
  W( 79 ) = W( 79 ) + a*JVS( 673 )
  W( 80 ) = W( 80 ) + a*JVS( 674 )
  W( 81 ) = W( 81 ) + a*JVS( 675 )
  a = -W( 77 ) / JVS(          700  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 701 )
  W( 79 ) = W( 79 ) + a*JVS( 702 )
  W( 80 ) = W( 80 ) + a*JVS( 703 )
  W( 81 ) = W( 81 ) + a*JVS( 704 )
  a = -W( 78 ) / JVS(          731  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 732 )
  W( 80 ) = W( 80 ) + a*JVS( 733 )
  W( 81 ) = W( 81 ) + a*JVS( 734 )
  JVS( 735) = W( 9 )
  JVS( 736) = W( 10 )
  JVS( 737) = W( 12 )
  JVS( 738) = W( 14 )
  JVS( 739) = W( 16 )
  JVS( 740) = W( 17 )
  JVS( 741) = W( 18 )
  JVS( 742) = W( 19 )
  JVS( 743) = W( 20 )
  JVS( 744) = W( 21 )
  JVS( 745) = W( 22 )
  JVS( 746) = W( 24 )
  JVS( 747) = W( 26 )
  JVS( 748) = W( 27 )
  JVS( 749) = W( 30 )
  JVS( 750) = W( 35 )
  JVS( 751) = W( 37 )
  JVS( 752) = W( 38 )
  JVS( 753) = W( 40 )
  JVS( 754) = W( 41 )
  JVS( 755) = W( 42 )
  JVS( 756) = W( 43 )
  JVS( 757) = W( 44 )
  JVS( 758) = W( 45 )
  JVS( 759) = W( 46 )
  JVS( 760) = W( 47 )
  JVS( 761) = W( 48 )
  JVS( 762) = W( 49 )
  JVS( 763) = W( 50 )
  JVS( 764) = W( 52 )
  JVS( 765) = W( 53 )
  JVS( 766) = W( 54 )
  JVS( 767) = W( 55 )
  JVS( 768) = W( 56 )
  JVS( 769) = W( 57 )
  JVS( 770) = W( 58 )
  JVS( 771) = W( 59 )
  JVS( 772) = W( 60 )
  JVS( 773) = W( 61 )
  JVS( 774) = W( 62 )
  JVS( 775) = W( 63 )
  JVS( 776) = W( 64 )
  JVS( 777) = W( 65 )
  JVS( 778) = W( 66 )
  JVS( 779) = W( 67 )
  JVS( 780) = W( 68 )
  JVS( 781) = W( 69 )
  JVS( 782) = W( 70 )
  JVS( 783) = W( 71 )
  JVS( 784) = W( 72 )
  JVS( 785) = W( 73 )
  JVS( 786) = W( 74 )
  JVS( 787) = W( 75 )
  JVS( 788) = W( 76 )
  JVS( 789) = W( 77 )
  JVS( 790) = W( 78 )
  JVS( 791) = W( 79 )
  JVS( 792) = W( 80 )
  JVS( 793) = W( 81 )
  IF ( ABS(  JVS( 819 )) < TINY(a) ) THEN
         IER = 80                                      
         RETURN
  END IF
   W( 13 ) = JVS( 794 )
   W( 16 ) = JVS( 795 )
   W( 26 ) = JVS( 796 )
   W( 36 ) = JVS( 797 )
   W( 39 ) = JVS( 798 )
   W( 43 ) = JVS( 799 )
   W( 46 ) = JVS( 800 )
   W( 56 ) = JVS( 801 )
   W( 59 ) = JVS( 802 )
   W( 60 ) = JVS( 803 )
   W( 61 ) = JVS( 804 )
   W( 62 ) = JVS( 805 )
   W( 67 ) = JVS( 806 )
   W( 68 ) = JVS( 807 )
   W( 69 ) = JVS( 808 )
   W( 70 ) = JVS( 809 )
   W( 71 ) = JVS( 810 )
   W( 72 ) = JVS( 811 )
   W( 73 ) = JVS( 812 )
   W( 74 ) = JVS( 813 )
   W( 75 ) = JVS( 814 )
   W( 76 ) = JVS( 815 )
   W( 77 ) = JVS( 816 )
   W( 78 ) = JVS( 817 )
   W( 79 ) = JVS( 818 )
   W( 80 ) = JVS( 819 )
   W( 81 ) = JVS( 820 )
  a = -W( 13 ) / JVS(           31  )
  W( 13 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 32 )
  W( 80 ) = W( 80 ) + a*JVS( 33 )
  a = -W( 16 ) / JVS(           42  )
  W( 16 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 43 )
  W( 80 ) = W( 80 ) + a*JVS( 44 )
  a = -W( 26 ) / JVS(           82  )
  W( 26 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 83 )
  W( 75 ) = W( 75 ) + a*JVS( 84 )
  W( 79 ) = W( 79 ) + a*JVS( 85 )
  a = -W( 36 ) / JVS(          136  )
  W( 36 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 137 )
  W( 75 ) = W( 75 ) + a*JVS( 138 )
  W( 78 ) = W( 78 ) + a*JVS( 139 )
  a = -W( 39 ) / JVS(          154  )
  W( 39 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 155 )
  W( 67 ) = W( 67 ) + a*JVS( 156 )
  W( 72 ) = W( 72 ) + a*JVS( 157 )
  W( 73 ) = W( 73 ) + a*JVS( 158 )
  W( 75 ) = W( 75 ) + a*JVS( 159 )
  W( 80 ) = W( 80 ) + a*JVS( 160 )
  a = -W( 43 ) / JVS(          178  )
  W( 43 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 179 )
  W( 74 ) = W( 74 ) + a*JVS( 180 )
  W( 75 ) = W( 75 ) + a*JVS( 181 )
  a = -W( 46 ) / JVS(          196  )
  W( 46 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 197 )
  W( 80 ) = W( 80 ) + a*JVS( 198 )
  W( 81 ) = W( 81 ) + a*JVS( 199 )
  a = -W( 56 ) / JVS(          273  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 274 )
  W( 76 ) = W( 76 ) + a*JVS( 275 )
  W( 79 ) = W( 79 ) + a*JVS( 276 )
  W( 80 ) = W( 80 ) + a*JVS( 277 )
  a = -W( 59 ) / JVS(          309  )
  W( 59 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 310 )
  W( 80 ) = W( 80 ) + a*JVS( 311 )
  W( 81 ) = W( 81 ) + a*JVS( 312 )
  a = -W( 60 ) / JVS(          315  )
  W( 60 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 316 )
  W( 71 ) = W( 71 ) + a*JVS( 317 )
  W( 75 ) = W( 75 ) + a*JVS( 318 )
  W( 76 ) = W( 76 ) + a*JVS( 319 )
  W( 79 ) = W( 79 ) + a*JVS( 320 )
  W( 80 ) = W( 80 ) + a*JVS( 321 )
  W( 81 ) = W( 81 ) + a*JVS( 322 )
  a = -W( 61 ) / JVS(          324  )
  W( 61 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 325 )
  W( 75 ) = W( 75 ) + a*JVS( 326 )
  W( 80 ) = W( 80 ) + a*JVS( 327 )
  W( 81 ) = W( 81 ) + a*JVS( 328 )
  a = -W( 62 ) / JVS(          341  )
  W( 62 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 342 )
  W( 75 ) = W( 75 ) + a*JVS( 343 )
  W( 76 ) = W( 76 ) + a*JVS( 344 )
  W( 77 ) = W( 77 ) + a*JVS( 345 )
  W( 79 ) = W( 79 ) + a*JVS( 346 )
  W( 80 ) = W( 80 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  a = -W( 67 ) / JVS(          411  )
  W( 67 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 412 )
  W( 69 ) = W( 69 ) + a*JVS( 413 )
  W( 70 ) = W( 70 ) + a*JVS( 414 )
  W( 71 ) = W( 71 ) + a*JVS( 415 )
  W( 73 ) = W( 73 ) + a*JVS( 416 )
  W( 74 ) = W( 74 ) + a*JVS( 417 )
  W( 75 ) = W( 75 ) + a*JVS( 418 )
  W( 76 ) = W( 76 ) + a*JVS( 419 )
  W( 77 ) = W( 77 ) + a*JVS( 420 )
  W( 78 ) = W( 78 ) + a*JVS( 421 )
  W( 79 ) = W( 79 ) + a*JVS( 422 )
  W( 80 ) = W( 80 ) + a*JVS( 423 )
  W( 81 ) = W( 81 ) + a*JVS( 424 )
  a = -W( 68 ) / JVS(          429  )
  W( 68 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 430 )
  W( 71 ) = W( 71 ) + a*JVS( 431 )
  W( 75 ) = W( 75 ) + a*JVS( 432 )
  W( 76 ) = W( 76 ) + a*JVS( 433 )
  W( 77 ) = W( 77 ) + a*JVS( 434 )
  W( 78 ) = W( 78 ) + a*JVS( 435 )
  W( 79 ) = W( 79 ) + a*JVS( 436 )
  W( 80 ) = W( 80 ) + a*JVS( 437 )
  W( 81 ) = W( 81 ) + a*JVS( 438 )
  a = -W( 69 ) / JVS(          445  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 446 )
  W( 75 ) = W( 75 ) + a*JVS( 447 )
  W( 76 ) = W( 76 ) + a*JVS( 448 )
  W( 77 ) = W( 77 ) + a*JVS( 449 )
  W( 78 ) = W( 78 ) + a*JVS( 450 )
  W( 79 ) = W( 79 ) + a*JVS( 451 )
  W( 80 ) = W( 80 ) + a*JVS( 452 )
  W( 81 ) = W( 81 ) + a*JVS( 453 )
  a = -W( 70 ) / JVS(          457  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 458 )
  W( 75 ) = W( 75 ) + a*JVS( 459 )
  W( 76 ) = W( 76 ) + a*JVS( 460 )
  W( 77 ) = W( 77 ) + a*JVS( 461 )
  W( 78 ) = W( 78 ) + a*JVS( 462 )
  W( 79 ) = W( 79 ) + a*JVS( 463 )
  W( 80 ) = W( 80 ) + a*JVS( 464 )
  W( 81 ) = W( 81 ) + a*JVS( 465 )
  a = -W( 71 ) / JVS(          468  )
  W( 71 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 469 )
  W( 76 ) = W( 76 ) + a*JVS( 470 )
  W( 77 ) = W( 77 ) + a*JVS( 471 )
  W( 78 ) = W( 78 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  W( 80 ) = W( 80 ) + a*JVS( 474 )
  W( 81 ) = W( 81 ) + a*JVS( 475 )
  a = -W( 72 ) / JVS(          508  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 509 )
  W( 74 ) = W( 74 ) + a*JVS( 510 )
  W( 75 ) = W( 75 ) + a*JVS( 511 )
  W( 76 ) = W( 76 ) + a*JVS( 512 )
  W( 77 ) = W( 77 ) + a*JVS( 513 )
  W( 78 ) = W( 78 ) + a*JVS( 514 )
  W( 79 ) = W( 79 ) + a*JVS( 515 )
  W( 80 ) = W( 80 ) + a*JVS( 516 )
  W( 81 ) = W( 81 ) + a*JVS( 517 )
  a = -W( 73 ) / JVS(          548  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 76 ) = W( 76 ) + a*JVS( 551 )
  W( 77 ) = W( 77 ) + a*JVS( 552 )
  W( 78 ) = W( 78 ) + a*JVS( 553 )
  W( 79 ) = W( 79 ) + a*JVS( 554 )
  W( 80 ) = W( 80 ) + a*JVS( 555 )
  W( 81 ) = W( 81 ) + a*JVS( 556 )
  a = -W( 74 ) / JVS(          564  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 565 )
  W( 76 ) = W( 76 ) + a*JVS( 566 )
  W( 77 ) = W( 77 ) + a*JVS( 567 )
  W( 78 ) = W( 78 ) + a*JVS( 568 )
  W( 79 ) = W( 79 ) + a*JVS( 569 )
  W( 80 ) = W( 80 ) + a*JVS( 570 )
  W( 81 ) = W( 81 ) + a*JVS( 571 )
  a = -W( 75 ) / JVS(          639  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 640 )
  W( 77 ) = W( 77 ) + a*JVS( 641 )
  W( 78 ) = W( 78 ) + a*JVS( 642 )
  W( 79 ) = W( 79 ) + a*JVS( 643 )
  W( 80 ) = W( 80 ) + a*JVS( 644 )
  W( 81 ) = W( 81 ) + a*JVS( 645 )
  a = -W( 76 ) / JVS(          670  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 671 )
  W( 78 ) = W( 78 ) + a*JVS( 672 )
  W( 79 ) = W( 79 ) + a*JVS( 673 )
  W( 80 ) = W( 80 ) + a*JVS( 674 )
  W( 81 ) = W( 81 ) + a*JVS( 675 )
  a = -W( 77 ) / JVS(          700  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 701 )
  W( 79 ) = W( 79 ) + a*JVS( 702 )
  W( 80 ) = W( 80 ) + a*JVS( 703 )
  W( 81 ) = W( 81 ) + a*JVS( 704 )
  a = -W( 78 ) / JVS(          731  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 732 )
  W( 80 ) = W( 80 ) + a*JVS( 733 )
  W( 81 ) = W( 81 ) + a*JVS( 734 )
  a = -W( 79 ) / JVS(          791  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 792 )
  W( 81 ) = W( 81 ) + a*JVS( 793 )
  JVS( 794) = W( 13 )
  JVS( 795) = W( 16 )
  JVS( 796) = W( 26 )
  JVS( 797) = W( 36 )
  JVS( 798) = W( 39 )
  JVS( 799) = W( 43 )
  JVS( 800) = W( 46 )
  JVS( 801) = W( 56 )
  JVS( 802) = W( 59 )
  JVS( 803) = W( 60 )
  JVS( 804) = W( 61 )
  JVS( 805) = W( 62 )
  JVS( 806) = W( 67 )
  JVS( 807) = W( 68 )
  JVS( 808) = W( 69 )
  JVS( 809) = W( 70 )
  JVS( 810) = W( 71 )
  JVS( 811) = W( 72 )
  JVS( 812) = W( 73 )
  JVS( 813) = W( 74 )
  JVS( 814) = W( 75 )
  JVS( 815) = W( 76 )
  JVS( 816) = W( 77 )
  JVS( 817) = W( 78 )
  JVS( 818) = W( 79 )
  JVS( 819) = W( 80 )
  JVS( 820) = W( 81 )
  IF ( ABS(  JVS( 838 )) < TINY(a) ) THEN
         IER = 81                                      
         RETURN
  END IF
   W( 27 ) = JVS( 821 )
   W( 46 ) = JVS( 822 )
   W( 52 ) = JVS( 823 )
   W( 59 ) = JVS( 824 )
   W( 61 ) = JVS( 825 )
   W( 66 ) = JVS( 826 )
   W( 69 ) = JVS( 827 )
   W( 71 ) = JVS( 828 )
   W( 72 ) = JVS( 829 )
   W( 73 ) = JVS( 830 )
   W( 74 ) = JVS( 831 )
   W( 75 ) = JVS( 832 )
   W( 76 ) = JVS( 833 )
   W( 77 ) = JVS( 834 )
   W( 78 ) = JVS( 835 )
   W( 79 ) = JVS( 836 )
   W( 80 ) = JVS( 837 )
   W( 81 ) = JVS( 838 )
  a = -W( 27 ) / JVS(           86  )
  W( 27 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 87 )
  W( 81 ) = W( 81 ) + a*JVS( 88 )
  a = -W( 46 ) / JVS(          196  )
  W( 46 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 197 )
  W( 80 ) = W( 80 ) + a*JVS( 198 )
  W( 81 ) = W( 81 ) + a*JVS( 199 )
  a = -W( 52 ) / JVS(          239  )
  W( 52 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 240 )
  W( 72 ) = W( 72 ) + a*JVS( 241 )
  W( 73 ) = W( 73 ) + a*JVS( 242 )
  W( 75 ) = W( 75 ) + a*JVS( 243 )
  W( 79 ) = W( 79 ) + a*JVS( 244 )
  W( 81 ) = W( 81 ) + a*JVS( 245 )
  a = -W( 59 ) / JVS(          309  )
  W( 59 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 310 )
  W( 80 ) = W( 80 ) + a*JVS( 311 )
  W( 81 ) = W( 81 ) + a*JVS( 312 )
  a = -W( 61 ) / JVS(          324  )
  W( 61 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 325 )
  W( 75 ) = W( 75 ) + a*JVS( 326 )
  W( 80 ) = W( 80 ) + a*JVS( 327 )
  W( 81 ) = W( 81 ) + a*JVS( 328 )
  a = -W( 66 ) / JVS(          392  )
  W( 66 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 393 )
  W( 75 ) = W( 75 ) + a*JVS( 394 )
  W( 76 ) = W( 76 ) + a*JVS( 395 )
  W( 77 ) = W( 77 ) + a*JVS( 396 )
  W( 78 ) = W( 78 ) + a*JVS( 397 )
  W( 79 ) = W( 79 ) + a*JVS( 398 )
  W( 80 ) = W( 80 ) + a*JVS( 399 )
  W( 81 ) = W( 81 ) + a*JVS( 400 )
  a = -W( 69 ) / JVS(          445  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 446 )
  W( 75 ) = W( 75 ) + a*JVS( 447 )
  W( 76 ) = W( 76 ) + a*JVS( 448 )
  W( 77 ) = W( 77 ) + a*JVS( 449 )
  W( 78 ) = W( 78 ) + a*JVS( 450 )
  W( 79 ) = W( 79 ) + a*JVS( 451 )
  W( 80 ) = W( 80 ) + a*JVS( 452 )
  W( 81 ) = W( 81 ) + a*JVS( 453 )
  a = -W( 71 ) / JVS(          468  )
  W( 71 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 469 )
  W( 76 ) = W( 76 ) + a*JVS( 470 )
  W( 77 ) = W( 77 ) + a*JVS( 471 )
  W( 78 ) = W( 78 ) + a*JVS( 472 )
  W( 79 ) = W( 79 ) + a*JVS( 473 )
  W( 80 ) = W( 80 ) + a*JVS( 474 )
  W( 81 ) = W( 81 ) + a*JVS( 475 )
  a = -W( 72 ) / JVS(          508  )
  W( 72 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 509 )
  W( 74 ) = W( 74 ) + a*JVS( 510 )
  W( 75 ) = W( 75 ) + a*JVS( 511 )
  W( 76 ) = W( 76 ) + a*JVS( 512 )
  W( 77 ) = W( 77 ) + a*JVS( 513 )
  W( 78 ) = W( 78 ) + a*JVS( 514 )
  W( 79 ) = W( 79 ) + a*JVS( 515 )
  W( 80 ) = W( 80 ) + a*JVS( 516 )
  W( 81 ) = W( 81 ) + a*JVS( 517 )
  a = -W( 73 ) / JVS(          548  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 549 )
  W( 75 ) = W( 75 ) + a*JVS( 550 )
  W( 76 ) = W( 76 ) + a*JVS( 551 )
  W( 77 ) = W( 77 ) + a*JVS( 552 )
  W( 78 ) = W( 78 ) + a*JVS( 553 )
  W( 79 ) = W( 79 ) + a*JVS( 554 )
  W( 80 ) = W( 80 ) + a*JVS( 555 )
  W( 81 ) = W( 81 ) + a*JVS( 556 )
  a = -W( 74 ) / JVS(          564  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 565 )
  W( 76 ) = W( 76 ) + a*JVS( 566 )
  W( 77 ) = W( 77 ) + a*JVS( 567 )
  W( 78 ) = W( 78 ) + a*JVS( 568 )
  W( 79 ) = W( 79 ) + a*JVS( 569 )
  W( 80 ) = W( 80 ) + a*JVS( 570 )
  W( 81 ) = W( 81 ) + a*JVS( 571 )
  a = -W( 75 ) / JVS(          639  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 640 )
  W( 77 ) = W( 77 ) + a*JVS( 641 )
  W( 78 ) = W( 78 ) + a*JVS( 642 )
  W( 79 ) = W( 79 ) + a*JVS( 643 )
  W( 80 ) = W( 80 ) + a*JVS( 644 )
  W( 81 ) = W( 81 ) + a*JVS( 645 )
  a = -W( 76 ) / JVS(          670  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 671 )
  W( 78 ) = W( 78 ) + a*JVS( 672 )
  W( 79 ) = W( 79 ) + a*JVS( 673 )
  W( 80 ) = W( 80 ) + a*JVS( 674 )
  W( 81 ) = W( 81 ) + a*JVS( 675 )
  a = -W( 77 ) / JVS(          700  )
  W( 77 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 701 )
  W( 79 ) = W( 79 ) + a*JVS( 702 )
  W( 80 ) = W( 80 ) + a*JVS( 703 )
  W( 81 ) = W( 81 ) + a*JVS( 704 )
  a = -W( 78 ) / JVS(          731  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 732 )
  W( 80 ) = W( 80 ) + a*JVS( 733 )
  W( 81 ) = W( 81 ) + a*JVS( 734 )
  a = -W( 79 ) / JVS(          791  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 792 )
  W( 81 ) = W( 81 ) + a*JVS( 793 )
  a = -W( 80 ) / JVS(          819  )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 820 )
  JVS( 821) = W( 27 )
  JVS( 822) = W( 46 )
  JVS( 823) = W( 52 )
  JVS( 824) = W( 59 )
  JVS( 825) = W( 61 )
  JVS( 826) = W( 66 )
  JVS( 827) = W( 69 )
  JVS( 828) = W( 71 )
  JVS( 829) = W( 72 )
  JVS( 830) = W( 73 )
  JVS( 831) = W( 74 )
  JVS( 832) = W( 75 )
  JVS( 833) = W( 76 )
  JVS( 834) = W( 77 )
  JVS( 835) = W( 78 )
  JVS( 836) = W( 79 )
  JVS( 837) = W( 80 )
  JVS( 838) = W( 81 )
   
   END SUBROUTINE decomp_mozart
 


END MODULE mozart_Integrator
