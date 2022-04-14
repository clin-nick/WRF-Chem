
































MODULE racm_soa_vbs_aqchem_Integrator

 USE racm_soa_vbs_aqchem_Parameters
 USE racm_soa_vbs_aqchem_Precision
 USE racm_soa_vbs_aqchem_JacobianSP

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

SUBROUTINE  racm_soa_vbs_aqchem_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE racm_soa_vbs_aqchem_Parameters

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

   CALL racm_soa_vbs_aqchem_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  racm_soa_vbs_aqchem_INTEGRATE


SUBROUTINE  racm_soa_vbs_aqchem_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE racm_soa_vbs_aqchem_Parameters

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
      CALL racm_soa_vbs_aqchem_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL racm_soa_vbs_aqchem_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = racm_soa_vbs_aqchem_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL racm_soa_vbs_aqchem_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL racm_soa_vbs_aqchem_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL racm_soa_vbs_aqchem_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL racm_soa_vbs_aqchem_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL racm_soa_vbs_aqchem_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL racm_soa_vbs_aqchem_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL racm_soa_vbs_aqchem_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL racm_soa_vbs_aqchem_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL racm_soa_vbs_aqchem_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL racm_soa_vbs_aqchem_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL racm_soa_vbs_aqchem_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL racm_soa_vbs_aqchem_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL racm_soa_vbs_aqchem_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL racm_soa_vbs_aqchem_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL racm_soa_vbs_aqchem_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  racm_soa_vbs_aqchem_ros_ErrorMsg(Code,T,H,IERR)



   USE racm_soa_vbs_aqchem_Precision

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

 END SUBROUTINE  racm_soa_vbs_aqchem_ros_ErrorMsg


 SUBROUTINE  racm_soa_vbs_aqchem_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL racm_soa_vbs_aqchem_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL racm_soa_vbs_aqchem_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL racm_soa_vbs_aqchem_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL racm_soa_vbs_aqchem_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL racm_soa_vbs_aqchem_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL racm_soa_vbs_aqchem_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL racm_soa_vbs_aqchem_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL racm_soa_vbs_aqchem_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL racm_soa_vbs_aqchem_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL racm_soa_vbs_aqchem_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL racm_soa_vbs_aqchem_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL racm_soa_vbs_aqchem_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL racm_soa_vbs_aqchem_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL racm_soa_vbs_aqchem_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL racm_soa_vbs_aqchem_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL racm_soa_vbs_aqchem_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL racm_soa_vbs_aqchem_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL racm_soa_vbs_aqchem_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL racm_soa_vbs_aqchem_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL racm_soa_vbs_aqchem_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = racm_soa_vbs_aqchem_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL racm_soa_vbs_aqchem_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  racm_soa_vbs_aqchem_ros_Integrator



  REAL(kind=dp) FUNCTION  racm_soa_vbs_aqchem_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    racm_soa_vbs_aqchem_ros_ErrorNorm = Err

  END FUNCTION  racm_soa_vbs_aqchem_ros_ErrorNorm



  SUBROUTINE racm_soa_vbs_aqchem_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL racm_soa_vbs_aqchem_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL racm_soa_vbs_aqchem_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL racm_soa_vbs_aqchem_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  racm_soa_vbs_aqchem_ros_FunTimeDeriv



  SUBROUTINE  racm_soa_vbs_aqchem_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL racm_soa_vbs_aqchem_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL racm_soa_vbs_aqchem_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL racm_soa_vbs_aqchem_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  racm_soa_vbs_aqchem_ros_PrepareMatrix



  SUBROUTINE  racm_soa_vbs_aqchem_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 



CALL decomp_racm_soa_vbs_aqchem ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  racm_soa_vbs_aqchem_ros_Decomp



  SUBROUTINE  racm_soa_vbs_aqchem_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL racm_soa_vbs_aqchem_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  racm_soa_vbs_aqchem_ros_Solve




  SUBROUTINE  racm_soa_vbs_aqchem_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  racm_soa_vbs_aqchem_Ros2



  SUBROUTINE  racm_soa_vbs_aqchem_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racm_soa_vbs_aqchem_Ros3





  SUBROUTINE  racm_soa_vbs_aqchem_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racm_soa_vbs_aqchem_Ros4


  SUBROUTINE  racm_soa_vbs_aqchem_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racm_soa_vbs_aqchem_Rodas3


  SUBROUTINE  racm_soa_vbs_aqchem_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racm_soa_vbs_aqchem_Rodas4




END SUBROUTINE  racm_soa_vbs_aqchem_Rosenbrock




SUBROUTINE  racm_soa_vbs_aqchem_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE racm_soa_vbs_aqchem_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL racm_soa_vbs_aqchem_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  racm_soa_vbs_aqchem_FunTemplate



SUBROUTINE  racm_soa_vbs_aqchem_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE racm_soa_vbs_aqchem_Parameters
 
 USE racm_soa_vbs_aqchem_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL racm_soa_vbs_aqchem_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  racm_soa_vbs_aqchem_JacTemplate

















SUBROUTINE racm_soa_vbs_aqchem_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(88)
  A(2) = RCT(2)*V(79)
  A(3) = RCT(3)*V(79)
  A(4) = RCT(4)*V(31)
  A(5) = RCT(5)*V(43)
  A(6) = RCT(6)*V(29)
  A(7) = RCT(7)*V(89)
  A(8) = RCT(8)*V(89)
  A(9) = RCT(9)*V(35)
  A(10) = RCT(10)*V(66)
  A(11) = RCT(11)*V(66)
  A(12) = RCT(12)*V(74)
  A(13) = RCT(13)*V(27)
  A(14) = RCT(14)*V(86)
  A(15) = RCT(15)*V(40)
  A(16) = RCT(16)*V(59)
  A(17) = RCT(17)*V(53)
  A(18) = RCT(18)*V(53)
  A(19) = RCT(19)*V(62)
  A(20) = RCT(20)*V(71)
  A(21) = RCT(21)*V(70)
  A(22) = RCT(22)*V(57)
  A(23) = RCT(23)*V(32)
  A(24) = RCT(24)*V(33)*F(2)
  A(25) = RCT(25)*V(33)*V(79)
  A(26) = RCT(26)*V(13)*F(2)
  A(27) = RCT(27)*V(13)*F(1)
  A(28) = RCT(28)*V(79)*V(84)
  A(29) = RCT(29)*V(79)*V(85)
  A(30) = RCT(30)*V(84)*V(85)
  A(31) = RCT(31)*V(35)*V(84)
  A(32) = RCT(32)*V(85)*V(85)
  A(33) = RCT(33)*V(85)*V(85)*F(1)
  A(34) = RCT(34)*V(33)*V(87)
  A(35) = RCT(35)*V(33)*V(88)
  A(36) = RCT(36)*V(33)*V(88)
  A(37) = RCT(37)*V(84)*V(87)
  A(38) = RCT(38)*V(84)*V(88)
  A(39) = RCT(39)*V(84)*V(89)
  A(40) = RCT(40)*V(85)*V(87)
  A(41) = RCT(41)*V(85)*V(88)
  A(42) = RCT(42)*V(29)
  A(43) = RCT(43)*V(85)*V(89)
  A(44) = RCT(44)*V(31)*V(84)
  A(45) = RCT(45)*V(43)*V(84)
  A(46) = RCT(46)*V(29)*V(84)
  A(47) = RCT(47)*V(79)*V(87)
  A(48) = RCT(48)*V(79)*V(88)
  A(49) = RCT(49)*V(87)*V(87)*F(2)
  A(50) = RCT(50)*V(87)*V(89)
  A(51) = RCT(51)*V(88)*V(89)
  A(52) = RCT(52)*V(88)*V(89)
  A(53) = RCT(53)*V(19)
  A(54) = RCT(54)*V(89)*V(89)
  A(55) = RCT(55)*V(84)*F(2)
  A(56) = RCT(56)*V(14)*V(84)
  A(57) = RCT(57)*V(46)*V(84)
  A(58) = RCT(58)*V(24)*V(84)
  A(59) = RCT(59)*V(30)*V(84)
  A(60) = RCT(60)*V(26)*V(84)
  A(61) = RCT(61)*V(25)*V(84)
  A(62) = RCT(62)*V(22)*V(84)
  A(63) = RCT(63)*V(18)*V(84)
  A(64) = RCT(64)*V(20)*V(84)
  A(65) = RCT(65)*V(38)*V(84)
  A(66) = RCT(66)*V(75)*V(84)
  A(67) = RCT(67)*V(73)*V(84)
  A(68) = RCT(68)*V(51)*V(84)
  A(69) = RCT(69)*V(48)*V(84)
  A(70) = RCT(70)*V(45)*V(84)
  A(71) = RCT(71)*V(47)*V(84)
  A(72) = RCT(72)*V(16)*V(84)
  A(73) = RCT(73)*V(17)*V(84)
  A(74) = RCT(74)*V(50)*V(84)
  A(75) = RCT(75)*V(66)*V(84)
  A(76) = RCT(76)*V(74)*V(84)
  A(77) = RCT(77)*V(59)*V(84)
  A(78) = RCT(78)*V(32)*V(84)
  A(79) = RCT(79)*V(53)*V(84)
  A(80) = RCT(80)*V(62)*V(84)
  A(81) = RCT(81)*V(57)*V(84)
  A(82) = RCT(82)*V(71)*V(84)
  A(83) = RCT(83)*V(23)*V(84)
  A(84) = RCT(84)*V(27)*V(84)
  A(85) = RCT(85)*V(84)*V(86)
  A(86) = RCT(86)*V(40)*V(84)
  A(87) = RCT(87)*V(44)*V(84)
  A(88) = RCT(88)*V(54)*V(84)
  A(89) = RCT(89)*V(70)*V(84)
  A(90) = RCT(90)*V(66)*V(89)
  A(91) = RCT(91)*V(74)*V(89)
  A(92) = RCT(92)*V(53)*V(89)
  A(93) = RCT(93)*V(62)*V(89)
  A(94) = RCT(94)*V(21)*V(84)
  A(95) = RCT(95)*V(71)*V(89)
  A(96) = RCT(96)*V(50)*V(89)
  A(97) = RCT(97)*V(38)*V(89)
  A(98) = RCT(98)*V(75)*V(89)
  A(99) = RCT(99)*V(73)*V(89)
  A(100) = RCT(100)*V(51)*V(89)
  A(101) = RCT(101)*V(48)*V(89)
  A(102) = RCT(102)*V(45)*V(89)
  A(103) = RCT(103)*V(47)*V(89)
  A(104) = RCT(104)*V(54)*V(89)
  A(105) = RCT(105)*V(38)*V(79)
  A(106) = RCT(106)*V(75)*V(79)
  A(107) = RCT(107)*V(73)*V(79)
  A(108) = RCT(108)*V(51)*V(79)
  A(109) = RCT(109)*V(48)*V(79)
  A(110) = RCT(110)*V(45)*V(79)
  A(111) = RCT(111)*V(47)*V(79)
  A(112) = RCT(112)*V(57)*V(79)
  A(113) = RCT(113)*V(71)*V(79)
  A(114) = RCT(114)*V(54)*V(79)
  A(115) = RCT(115)*V(34)*V(88)
  A(116) = RCT(116)*V(34)*V(85)
  A(117) = RCT(117)*V(36)*V(88)
  A(118) = RCT(118)*V(36)*F(2)
  A(119) = RCT(119)*V(36)*V(79)
  A(120) = RCT(120)*V(37)*V(88)
  A(121) = RCT(121)*V(37)*F(2)
  A(122) = RCT(122)*V(37)*V(79)
  A(123) = RCT(123)*V(39)*V(88)
  A(124) = RCT(124)*V(39)*F(2)
  A(125) = RCT(125)*V(39)*V(79)
  A(126) = RCT(126)*V(82)*V(88)
  A(127) = RCT(127)*V(44)
  A(128) = RCT(128)*V(78)*V(88)
  A(129) = RCT(129)*V(54)
  A(130) = RCT(130)*V(83)*V(87)
  A(131) = RCT(131)*V(87)*V(90)
  A(132) = RCT(132)*V(80)*V(87)
  A(133) = RCT(133)*V(64)*V(87)
  A(134) = RCT(134)*V(65)*V(87)
  A(135) = RCT(135)*V(55)*V(87)
  A(136) = RCT(136)*V(56)*V(87)
  A(137) = RCT(137)*V(69)*V(87)
  A(138) = RCT(138)*V(60)*V(87)
  A(139) = RCT(139)*V(58)*V(87)
  A(140) = RCT(140)*V(63)*V(87)
  A(141) = RCT(141)*V(67)*V(87)
  A(142) = RCT(142)*V(68)*V(87)
  A(143) = RCT(143)*V(61)*V(87)
  A(144) = RCT(144)*V(82)*V(87)
  A(145) = RCT(145)*V(78)*V(87)
  A(146) = RCT(146)*V(81)*V(87)
  A(147) = RCT(147)*V(77)*V(87)
  A(148) = RCT(148)*V(76)*V(87)
  A(149) = RCT(149)*V(83)*V(85)
  A(150) = RCT(150)*V(85)*V(90)
  A(151) = RCT(151)*V(80)*V(85)
  A(152) = RCT(152)*V(64)*V(85)
  A(153) = RCT(153)*V(65)*V(85)
  A(154) = RCT(154)*V(55)*V(85)
  A(155) = RCT(155)*V(56)*V(85)
  A(156) = RCT(156)*V(69)*V(85)
  A(157) = RCT(157)*V(60)*V(85)
  A(158) = RCT(158)*V(58)*V(85)
  A(159) = RCT(159)*V(63)*V(85)
  A(160) = RCT(160)*V(67)*V(85)
  A(161) = RCT(161)*V(68)*V(85)
  A(162) = RCT(162)*V(61)*V(85)
  A(163) = RCT(163)*V(82)*V(85)
  A(164) = RCT(164)*V(82)*V(85)
  A(165) = RCT(165)*V(78)*V(85)
  A(166) = RCT(166)*V(78)*V(85)
  A(167) = RCT(167)*V(81)*V(85)
  A(168) = RCT(168)*V(77)*V(85)
  A(169) = RCT(169)*V(76)*V(85)
  A(170) = RCT(170)*V(83)*V(83)
  A(171) = RCT(171)*V(83)*V(90)
  A(172) = RCT(172)*V(80)*V(83)
  A(173) = RCT(173)*V(64)*V(83)
  A(174) = RCT(174)*V(65)*V(83)
  A(175) = RCT(175)*V(55)*V(83)
  A(176) = RCT(176)*V(56)*V(83)
  A(177) = RCT(177)*V(69)*V(83)
  A(178) = RCT(178)*V(60)*V(83)
  A(179) = RCT(179)*V(58)*V(83)
  A(180) = RCT(180)*V(63)*V(83)
  A(181) = RCT(181)*V(67)*V(83)
  A(182) = RCT(182)*V(68)*V(83)
  A(183) = RCT(183)*V(61)*V(83)
  A(184) = RCT(184)*V(82)*V(83)
  A(185) = RCT(185)*V(82)*V(83)
  A(186) = RCT(186)*V(78)*V(83)
  A(187) = RCT(187)*V(78)*V(83)
  A(188) = RCT(188)*V(81)*V(83)
  A(189) = RCT(189)*V(77)*V(83)
  A(190) = RCT(190)*V(76)*V(83)
  A(191) = RCT(191)*V(82)*V(90)
  A(192) = RCT(192)*V(80)*V(82)
  A(193) = RCT(193)*V(64)*V(82)
  A(194) = RCT(194)*V(65)*V(82)
  A(195) = RCT(195)*V(55)*V(82)
  A(196) = RCT(196)*V(56)*V(82)
  A(197) = RCT(197)*V(69)*V(82)
  A(198) = RCT(198)*V(60)*V(82)
  A(199) = RCT(199)*V(58)*V(82)
  A(200) = RCT(200)*V(63)*V(82)
  A(201) = RCT(201)*V(67)*V(82)
  A(202) = RCT(202)*V(68)*V(82)
  A(203) = RCT(203)*V(61)*V(82)
  A(204) = RCT(204)*V(82)*V(82)
  A(205) = RCT(205)*V(78)*V(82)
  A(206) = RCT(206)*V(81)*V(82)
  A(207) = RCT(207)*V(77)*V(82)
  A(208) = RCT(208)*V(76)*V(82)
  A(209) = RCT(209)*V(77)*V(77)
  A(210) = RCT(210)*V(76)*V(77)
  A(211) = RCT(211)*V(76)*V(76)
  A(212) = RCT(212)*V(83)*V(89)
  A(213) = RCT(213)*V(89)*V(90)
  A(214) = RCT(214)*V(80)*V(89)
  A(215) = RCT(215)*V(64)*V(89)
  A(216) = RCT(216)*V(65)*V(89)
  A(217) = RCT(217)*V(55)*V(89)
  A(218) = RCT(218)*V(56)*V(89)
  A(219) = RCT(219)*V(69)*V(89)
  A(220) = RCT(220)*V(28)*V(84)
  A(221) = RCT(221)*V(58)*V(89)
  A(222) = RCT(222)*V(63)*V(89)
  A(223) = RCT(223)*V(67)*V(89)
  A(224) = RCT(224)*V(68)*V(89)
  A(225) = RCT(225)*V(61)*V(89)
  A(226) = RCT(226)*V(82)*V(89)
  A(227) = RCT(227)*V(78)*V(89)
  A(228) = RCT(228)*V(81)*V(89)
  A(229) = RCT(229)*V(77)*V(89)
  A(230) = RCT(230)*V(76)*V(89)
  A(231) = RCT(231)*V(72)*V(85)
  A(232) = RCT(232)*V(72)*V(83)
  A(233) = RCT(233)*V(72)*V(82)
  A(234) = RCT(234)*V(72)*V(72)
  A(235) = RCT(235)*V(72)*V(87)
  A(236) = RCT(236)*V(72)*V(89)
  A(237) = RCT(237)*V(60)*V(60)
  A(238) = RCT(238)*V(15)*V(84)
  A(239) = RCT(239)*V(41)*V(84)
  A(240) = RCT(240)*V(52)*V(87)
  A(241) = RCT(241)*V(52)*V(85)
  A(242) = RCT(242)*V(52)*V(52)
  A(243) = RCT(243)*V(52)*V(88)
  A(244) = RCT(244)*V(28)
  A(245) = RCT(245)*V(42)*V(84)
  A(246) = RCT(246)*V(42)*V(79)
  A(247) = RCT(247)*V(42)*V(89)
  A(248) = RCT(248)*V(49)*V(84)
  A(249) = RCT(249)*V(49)*V(89)
  A(250) = RCT(250)*V(49)*V(79)
  A(251) = RCT(251)*V(8)*V(84)
  A(252) = RCT(252)*V(7)*V(84)
  A(253) = RCT(253)*V(6)*V(84)
  A(254) = RCT(254)*V(12)*V(84)
  A(255) = RCT(255)*V(11)*V(84)
  A(256) = RCT(256)*V(10)*V(84)


  Vdot(1) = A(56)
  Vdot(2) = A(57)
  Vdot(3) = 0.036*A(62)+0.37*A(105)+0.14*A(106)+0.15*A(108)+0.28*A(109)+0.01*A(111)+0.45*A(112)+0.11*A(113)+0.11*A(114)&
              &+0.05*A(245)+0.039*A(246)
  Vdot(4) = 0.1*A(106)+0.14*A(107)+0.07*A(111)+0.21*A(113)+A(164)+A(166)+A(185)+A(187)+0.5*A(191)+0.499*A(192)+0.495&
              &*A(193)+0.495*A(194)+0.5*A(195)+0.499*A(196)+0.49*A(197)+0.494*A(198)+0.5*A(206)+0.5*A(207)+0.484*A(208)&
              &+0.053*A(246)+0.14*A(250)
  Vdot(5) = 1.075*A(253)
  Vdot(6) = 1.075*A(252)-A(253)
  Vdot(7) = 1.075*A(251)-A(252)
  Vdot(8) = -A(251)
  Vdot(9) = 1.075*A(256)
  Vdot(10) = 1.075*A(255)-A(256)
  Vdot(11) = 1.075*A(254)-A(255)
  Vdot(12) = -A(254)
  Vdot(13) = A(2)-A(26)-A(27)
  Vdot(14) = -A(56)
  Vdot(15) = A(157)-A(238)
  Vdot(16) = -A(72)
  Vdot(17) = -A(73)
  Vdot(18) = -A(63)
  Vdot(19) = A(52)-A(53)
  Vdot(20) = -A(64)
  Vdot(21) = -A(94)+A(241)
  Vdot(22) = -A(62)
  Vdot(23) = 0.35*A(82)-A(83)
  Vdot(24) = -A(58)+A(239)
  Vdot(25) = -A(61)+0.03*A(106)+0.06*A(107)+0.06*A(250)
  Vdot(26) = -A(60)+0.06*A(106)+0.07*A(107)+0.07*A(250)
  Vdot(27) = -A(13)-A(84)+A(149)
  Vdot(28) = -A(220)+A(243)-A(244)
  Vdot(29) = -A(6)+A(41)-A(42)-A(46)
  Vdot(30) = -A(59)+A(220)+A(239)+0.25*A(240)+A(242)
  Vdot(31) = -A(4)+A(37)-A(44)+A(117)+A(120)+A(123)
  Vdot(32) = -A(23)+0.024*A(64)-A(78)+0.6*A(88)+0.3*A(188)
  Vdot(33) = A(1)+A(3)+A(8)-A(24)-A(25)+A(26)-A(34)-A(35)-A(36)+0.09*A(108)
  Vdot(34) = 0.1*A(74)+A(96)-A(115)-A(116)
  Vdot(35) = -A(9)-A(31)+A(32)+A(33)+0.006*A(106)+0.011*A(107)+0.001*A(108)+0.09*A(109)+0.02*A(110)+0.02*A(111)+0.011&
               &*A(250)
  Vdot(36) = 0.9*A(72)-A(117)-A(118)-A(119)
  Vdot(37) = 0.9*A(73)-A(120)-A(121)-A(122)
  Vdot(38) = -A(65)-A(97)-A(105)
  Vdot(39) = 0.85*A(74)-A(123)-A(124)-A(125)
  Vdot(40) = -A(15)-A(86)+0.11*A(113)+A(163)
  Vdot(41) = A(101)+0.046*A(138)-A(239)
  Vdot(42) = -A(245)-A(246)-A(247)
  Vdot(43) = -A(5)+A(38)+0.3*A(43)-A(45)+A(90)+A(91)+A(92)+A(93)+0.5*A(95)+A(96)
  Vdot(44) = -A(87)+0.4*A(88)+0.4*A(104)+0.3*A(114)+A(126)-A(127)
  Vdot(45) = -A(70)-A(102)-A(110)
  Vdot(46) = A(10)+A(11)+A(12)+1.87*A(17)+1.55*A(18)+A(19)+A(22)-A(57)+A(58)+0.036*A(62)+A(75)+2*A(79)+A(80)+A(90)+2&
               &*A(92)+A(93)+0.43*A(105)+0.37*A(106)+0.3*A(107)+0.36*A(108)+0.14*A(109)+0.14*A(110)+0.14*A(111)+0.22*A(112)&
               &+0.66*A(113)+0.13*A(114)+0.25*A(240)+0.5*A(242)+0.3*A(250)
  Vdot(47) = -A(71)-A(103)-A(111)
  Vdot(48) = -A(69)-A(101)-A(109)
  Vdot(49) = -A(248)-A(249)-A(250)
  Vdot(50) = -A(74)-A(96)+0.1*A(115)+A(116)+A(117)+0.02*A(118)+A(119)+A(120)+0.02*A(121)+A(122)+A(123)+0.02*A(124)&
               &+A(125)
  Vdot(51) = -A(68)-A(100)-A(108)
  Vdot(52) = A(81)+A(94)+0.1*A(109)-A(240)-A(241)-2*A(242)-A(243)+A(244)
  Vdot(53) = -A(17)-A(18)+0.036*A(62)-A(79)+0.15*A(82)-A(92)+0.25*A(95)+0.5*A(113)+0.063*A(132)+1.2*A(141)+0.35*A(142)&
               &+A(143)+0.119*A(172)+0.65*A(181)+0.37*A(182)+A(183)+0.1*A(192)+0.65*A(201)+0.37*A(202)+A(203)+0.063*A(214)&
               &+1.3*A(223)+0.74*A(224)+A(225)
  Vdot(54) = -A(88)-A(104)-A(114)+A(128)-A(129)
  Vdot(55) = A(65)-A(135)-A(154)-A(175)-A(195)-A(217)
  Vdot(56) = A(66)-A(136)-A(155)-A(176)-A(196)-A(218)
  Vdot(57) = -A(22)-A(81)+0.9*A(100)+0.39*A(108)+0.65*A(109)+0.79*A(111)-A(112)+A(138)+0.4*A(140)+0.55*A(178)+0.6*A(180)&
               &+0.771*A(198)+0.6*A(200)+0.6*A(222)+2*A(237)+A(238)+0.9*A(247)
  Vdot(58) = A(70)-A(139)-A(158)-A(179)-A(199)-A(221)
  Vdot(59) = -A(16)+0.8*A(21)+0.25*A(63)-A(77)+0.12*A(83)+0.41*A(85)+0.03*A(95)+0.03*A(106)+0.16*A(107)+0.53*A(110)&
               &+0.623*A(132)+0.722*A(133)+0.642*A(134)+0.06*A(136)+0.29*A(137)+0.8*A(139)+0.464*A(148)+0.018*A(172)+0.24&
               &*A(173)+0.419*A(174)+0.081*A(176)+0.313*A(177)+A(179)+0.149*A(190)+0.127*A(192)+0.33*A(193)+0.581*A(194)&
               &+0.141*A(196)+0.569*A(197)+A(199)+0.11*A(206)+0.167*A(208)+0.149*A(210)+0.285*A(211)+0.67*A(214)+0.828&
               &*A(215)+0.88*A(216)+0.06*A(218)+0.29*A(219)+A(221)+0.469*A(230)+0.36*A(245)+0.23*A(246)+0.16*A(250)
  Vdot(60) = A(68)+A(69)-A(138)-A(157)-A(178)-A(198)-2*A(237)
  Vdot(61) = 0.98*A(124)-A(143)-A(162)-A(183)-A(203)-A(225)
  Vdot(62) = -A(19)+A(59)+A(78)-A(80)+0.15*A(82)-A(93)+0.25*A(95)+0.9*A(112)+0.62*A(113)+0.65*A(141)+0.6*A(142)+A(143)&
               &+0.54*A(146)+0.005*A(172)+0.35*A(181)+0.63*A(182)+A(183)+0.4*A(188)+0.004*A(192)+0.35*A(201)+0.63*A(202)&
               &+A(203)+0.54*A(206)+0.7*A(223)+1.26*A(224)+A(225)+0.54*A(228)+0.5*A(240)+A(242)
  Vdot(63) = A(71)-A(140)-A(159)-A(180)-A(200)-A(222)
  Vdot(64) = 0.75*A(63)-A(133)-A(152)-A(173)-A(193)-A(215)
  Vdot(65) = 0.9511*A(64)-A(134)-A(153)-A(174)-A(194)-A(216)
  Vdot(66) = -A(10)-A(11)+A(13)+0.13*A(17)+0.45*A(18)+A(22)+A(23)+A(58)+0.01*A(62)-A(75)+0.35*A(84)+0.35*A(86)+A(87)+0.4&
               &*A(88)-A(90)+0.4*A(104)+A(105)+0.64*A(106)+0.02*A(107)+0.9*A(108)+0.58*A(109)+0.04*A(111)+0.7*A(114)+A(130)&
               &+0.047*A(132)+0.021*A(133)+1.6*A(135)+A(136)+A(138)+0.25*A(140)+A(145)+0.287*A(148)+1.33*A(170)+0.75*A(171)&
               &+0.81*A(172)+0.829*A(173)+0.753*A(174)+1.55*A(175)+1.25*A(176)+0.755*A(177)+1.09*A(178)+A(179)+1.4*A(180)&
               &+A(181)+A(182)+A(183)+A(184)+A(185)+2*A(186)+A(187)+0.75*A(188)+0.75*A(189)+0.96*A(190)+0.091*A(192)+0.076&
               &*A(193)+0.8*A(195)+0.501*A(196)+0.34*A(198)+0.4*A(200)+A(205)+0.207*A(208)+0.202*A(210)+0.504*A(211)+A(212)&
               &+0.048*A(214)+0.021*A(215)+1.6*A(217)+A(218)+0.4*A(222)+A(227)+0.28*A(230)+A(232)+A(237)+0.75*A(240)+0.5&
               &*A(242)+0.3*A(245)+0.51*A(246)+0.02*A(250)
  Vdot(67) = 0.98*A(118)-A(141)-A(160)-A(181)-A(201)-A(223)
  Vdot(68) = 0.98*A(121)-A(142)-A(161)-A(182)-A(202)-A(224)
  Vdot(69) = A(67)-A(137)-A(156)-A(177)-A(197)-A(219)+0.19*A(245)+A(248)
  Vdot(70) = -A(21)-A(89)+0.6*A(104)+A(115)+0.059*A(132)+0.124*A(133)+0.261*A(134)+0.2*A(139)+0.35*A(140)+0.05*A(141)&
               &+0.05*A(142)+A(147)+A(168)+A(169)+A(189)+0.5*A(190)+A(207)+0.484*A(208)+2*A(209)+1.5*A(210)+A(211)+A(229)
  Vdot(71) = -A(20)-A(82)-A(95)-A(113)+0.5*A(141)+0.95*A(142)+A(181)+A(182)+A(201)+A(202)+0.5*A(223)+A(224)
  Vdot(72) = 0.1*A(72)+0.1*A(73)+0.05*A(74)+0.5*A(82)+0.07*A(85)+0.35*A(86)+A(87)+A(88)+0.5*A(95)+A(104)+0.13*A(108)&
               &+0.048*A(132)+0.334*A(133)+0.416*A(134)+0.16*A(146)+0.085*A(172)+0.245*A(173)+0.322*A(174)+0.08*A(188)+0.071&
               &*A(192)+0.237*A(193)+0.318*A(194)+0.08*A(206)+0.051*A(214)+0.391*A(215)+0.587*A(216)+0.16*A(228)-A(231)&
               &-A(232)-A(233)-2*A(234)-A(235)-A(236)
  Vdot(73) = -A(67)-A(99)-A(107)+0.25*A(140)+0.08*A(178)+0.4*A(180)+0.4*A(200)+0.4*A(222)
  Vdot(74) = -A(12)+A(14)+0.2*A(21)+0.335*A(62)+0.025*A(64)-A(76)+0.88*A(83)+0.08*A(85)-A(91)+0.25*A(95)+0.44*A(106)&
               &+0.99*A(107)+0.65*A(110)+0.16*A(113)+A(131)+0.233*A(132)+0.211*A(133)+0.15*A(134)+0.2*A(135)+0.94*A(136)&
               &+1.71*A(137)+0.8*A(139)+0.46*A(146)+1.24*A(148)+0.75*A(171)+0.58*A(172)+0.523*A(173)+0.411*A(174)+0.35&
               &*A(175)+0.669*A(176)+0.932*A(177)+A(179)+0.3*A(188)+0.64*A(190)+A(191)+0.724*A(192)+0.677*A(193)+0.497&
               &*A(194)+0.6*A(195)+0.859*A(196)+0.941*A(197)+A(199)+0.35*A(206)+0.65*A(208)+0.64*A(210)+1.21*A(211)+A(213)&
               &+0.243*A(214)+0.239*A(215)+0.187*A(216)+0.2*A(217)+0.94*A(218)+1.71*A(219)+A(221)+0.46*A(228)+1.24*A(230)&
               &+0.85*A(246)+0.99*A(250)
  Vdot(75) = -A(66)-A(98)-A(106)+0.35*A(108)+0.46*A(111)+0.37*A(178)+0.229*A(198)
  Vdot(76) = 0.2*A(97)+0.57*A(98)+0.89*A(99)+0.1*A(100)+0.9*A(102)+0.87*A(103)-A(148)-A(169)-A(190)-A(208)-A(210)-2&
               &*A(211)-A(230)+0.1*A(247)+0.89*A(249)
  Vdot(77) = 0.8*A(97)+0.43*A(98)+0.11*A(99)+0.9*A(100)+0.1*A(102)+0.13*A(103)-A(147)-A(168)-A(189)-A(207)-2*A(209)&
               &-A(210)-A(229)+0.9*A(247)+0.11*A(249)
  Vdot(78) = A(20)+0.5*A(82)+0.5*A(95)-A(128)+A(129)-A(145)-A(165)-A(166)-A(186)-A(187)-A(205)-A(227)
  Vdot(79) = -A(2)-A(3)+A(24)-A(25)-A(28)-A(29)-A(47)-A(48)-A(105)-A(106)-A(107)-A(108)-A(109)-A(110)-A(111)-A(112)&
               &-A(113)-A(114)-A(119)-A(122)-A(125)+A(164)+A(166)-A(246)-A(250)
  Vdot(80) = 0.583*A(62)+0.44*A(85)+A(89)-A(132)-A(151)-A(172)-A(192)-A(214)
  Vdot(81) = A(77)+0.03*A(106)+0.12*A(107)+0.02*A(108)+0.42*A(110)+0.42*A(111)-A(146)-A(167)-A(188)-A(206)-A(228)+0.12&
               &*A(250)
  Vdot(82) = A(16)+A(19)+A(22)+A(23)+A(76)+A(80)+0.65*A(86)+A(91)+A(93)+0.15*A(108)+0.1*A(109)+0.1*A(112)+0.28*A(113)&
               &+0.7*A(114)-A(126)+A(127)-A(144)+A(145)+0.23*A(146)-A(163)-A(164)-A(184)-A(185)+A(186)+0.12*A(188)-A(191)&
               &-A(192)-A(193)-A(194)-A(195)-A(196)-A(197)-A(198)-A(199)-A(200)-A(201)-A(202)-A(203)-2*A(204)-0.88*A(206)&
               &-A(207)-A(208)-A(226)+A(227)+0.23*A(228)-A(233)+0.25*A(240)
  Vdot(83) = A(12)+A(15)+A(60)+0.65*A(84)+0.19*A(106)+0.23*A(107)+0.03*A(108)+0.08*A(109)-A(130)+0.15*A(132)+0.031&
               &*A(133)+A(144)-A(149)-2*A(170)-A(171)-0.993*A(172)-0.951*A(173)-A(174)-A(175)-A(176)-A(177)-A(178)-A(179)&
               &-A(180)-A(181)-A(182)-A(183)-A(185)-A(186)-A(187)-A(188)-A(189)-A(190)+0.5*A(191)+0.508*A(192)+0.554*A(193)&
               &+0.507*A(194)+0.5*A(195)+0.501*A(196)+0.51*A(197)+0.506*A(198)+A(199)+A(200)+A(201)+A(202)+A(203)+2*A(204)&
               &+A(205)+0.5*A(206)+0.5*A(207)+0.516*A(208)-A(212)+0.155*A(214)+0.04*A(215)+A(226)-A(232)+A(233)+0.23*A(250)
  Vdot(84) = A(4)+A(5)+0.35*A(6)+2*A(9)+A(13)+A(14)+A(15)+2*A(27)-A(28)+A(29)-A(30)-A(31)-A(37)-A(38)-A(39)+A(40)+0.7&
               &*A(43)-A(44)-A(45)-A(46)-A(55)-A(56)-A(57)-A(58)-A(59)-A(60)-A(61)-0.964*A(62)-A(63)-A(64)-A(65)-A(66)-A(67)&
               &-A(68)-A(69)-A(70)-A(71)-A(72)-A(73)-A(74)-A(75)-A(76)-A(77)-A(78)-A(79)-A(80)-A(81)-A(82)-A(83)-0.65*A(84)&
               &-0.51*A(85)-A(86)-A(87)-A(88)-A(89)-A(94)+0.12*A(105)+0.4*A(106)+0.63*A(107)+0.28*A(108)+0.25*A(109)+0.85&
               &*A(110)+0.85*A(111)+0.19*A(112)+0.21*A(113)+0.036*A(114)+A(119)+A(122)+A(125)-A(220)-A(239)-A(245)+0.63&
               &*A(246)-A(248)+0.63*A(250)
  Vdot(85) = 0.65*A(6)+2*A(11)+A(12)+A(13)+A(14)+0.8*A(18)+A(19)+A(20)+A(21)+A(22)+A(23)+A(28)-A(29)-A(30)+A(31)-2*A(32)&
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
               &+A(232)+A(237)+0.75*A(240)-A(241)+A(242)+0.22*A(250)
  Vdot(86) = -A(14)-A(85)+A(150)+A(151)+A(152)+A(153)+A(154)+A(155)+A(156)+A(158)+A(159)+A(160)+A(161)+A(162)+A(165)&
               &+A(167)+A(231)
  Vdot(87) = A(1)+A(4)+A(7)-A(34)+A(35)-A(37)-A(40)-A(47)-2*A(49)-A(50)+A(51)-A(130)-A(131)-A(132)-A(133)-A(134)-A(135)&
               &-A(136)-A(137)-A(138)-A(139)-A(140)-A(141)-A(142)-A(143)-A(144)-A(145)-A(146)-A(147)-A(148)-A(235)-A(240)
  Vdot(88) = -A(1)+A(5)+0.65*A(6)+A(8)+A(21)+A(34)-A(35)-A(36)-A(38)+A(39)+A(40)-A(41)+A(42)+0.7*A(43)+A(44)+A(46)+A(47)&
               &-A(48)+2*A(49)+2*A(50)-A(52)+A(53)+2*A(54)+A(58)+A(89)+0.5*A(95)+0.4*A(104)+0.7*A(114)-A(115)-A(117)-A(120)&
               &-A(123)-A(126)+A(127)-A(128)+A(129)+A(130)+A(131)+0.941*A(132)+0.876*A(133)+0.739*A(134)+A(135)+A(136)&
               &+A(137)+A(138)+0.8*A(139)+0.65*A(140)+0.95*A(141)+0.95*A(142)+A(143)+A(144)+A(145)+A(146)+A(147)+2*A(148)&
               &+0.5*A(190)+0.516*A(208)+0.5*A(210)+A(211)+A(212)+A(213)+A(214)+A(215)+A(216)+A(217)+A(218)+A(219)+A(220)&
               &+A(221)+A(222)+A(223)+A(224)+A(225)+A(226)+A(227)+A(228)+A(229)+2*A(230)+A(235)+A(236)+A(240)-A(243)+A(244)
  Vdot(89) = 0.35*A(6)-A(7)-A(8)+A(36)-A(39)-A(43)+A(45)+A(48)-A(50)-A(51)-A(52)+A(53)-2*A(54)+A(87)+0.6*A(88)-A(90)&
               &-A(91)-A(92)-A(93)-A(95)-A(96)-A(97)-A(98)-A(99)-A(100)-A(101)-A(102)-A(103)-0.4*A(104)-A(212)-A(213)-A(214)&
               &-A(215)-A(216)-A(217)-A(218)-A(219)-A(221)-A(222)-A(223)-A(224)-A(225)-A(226)-A(227)-A(228)-A(229)-A(230)&
               &-A(236)-A(247)-A(249)
  Vdot(90) = A(16)+A(61)+0.1*A(106)+0.18*A(107)+0.2*A(110)+0.16*A(111)-A(131)+0.048*A(132)+0.245*A(133)+0.133*A(134)&
               &-A(150)-A(171)+0.014*A(173)+0.013*A(174)-A(191)+0.006*A(192)+0.018*A(193)+0.015*A(194)-A(213)+0.053*A(214)&
               &+0.262*A(215)+0.155*A(216)+0.18*A(250)
      
END SUBROUTINE racm_soa_vbs_aqchem_Fun
















SUBROUTINE racm_soa_vbs_aqchem_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(88)
  IRR(2) = RCT(2)*V(79)
  IRR(3) = RCT(3)*V(79)
  IRR(4) = RCT(4)*V(31)
  IRR(5) = RCT(5)*V(43)
  IRR(6) = RCT(6)*V(29)
  IRR(7) = RCT(7)*V(89)
  IRR(8) = RCT(8)*V(89)
  IRR(9) = RCT(9)*V(35)
  IRR(10) = RCT(10)*V(66)
  IRR(11) = RCT(11)*V(66)
  IRR(12) = RCT(12)*V(74)
  IRR(13) = RCT(13)*V(27)
  IRR(14) = RCT(14)*V(86)
  IRR(15) = RCT(15)*V(40)
  IRR(16) = RCT(16)*V(59)
  IRR(17) = RCT(17)*V(53)
  IRR(18) = RCT(18)*V(53)
  IRR(19) = RCT(19)*V(62)
  IRR(20) = RCT(20)*V(71)
  IRR(21) = RCT(21)*V(70)
  IRR(22) = RCT(22)*V(57)
  IRR(23) = RCT(23)*V(32)
  IRR(24) = RCT(24)*V(33)*F(2)
  IRR(25) = RCT(25)*V(33)*V(79)
  IRR(26) = RCT(26)*V(13)*F(2)
  IRR(27) = RCT(27)*V(13)*F(1)
  IRR(28) = RCT(28)*V(79)*V(84)
  IRR(29) = RCT(29)*V(79)*V(85)
  IRR(30) = RCT(30)*V(84)*V(85)
  IRR(31) = RCT(31)*V(35)*V(84)
  IRR(32) = RCT(32)*V(85)*V(85)
  IRR(33) = RCT(33)*V(85)*V(85)*F(1)
  IRR(34) = RCT(34)*V(33)*V(87)
  IRR(35) = RCT(35)*V(33)*V(88)
  IRR(36) = RCT(36)*V(33)*V(88)
  IRR(37) = RCT(37)*V(84)*V(87)
  IRR(38) = RCT(38)*V(84)*V(88)
  IRR(39) = RCT(39)*V(84)*V(89)
  IRR(40) = RCT(40)*V(85)*V(87)
  IRR(41) = RCT(41)*V(85)*V(88)
  IRR(42) = RCT(42)*V(29)
  IRR(43) = RCT(43)*V(85)*V(89)
  IRR(44) = RCT(44)*V(31)*V(84)
  IRR(45) = RCT(45)*V(43)*V(84)
  IRR(46) = RCT(46)*V(29)*V(84)
  IRR(47) = RCT(47)*V(79)*V(87)
  IRR(48) = RCT(48)*V(79)*V(88)
  IRR(49) = RCT(49)*V(87)*V(87)*F(2)
  IRR(50) = RCT(50)*V(87)*V(89)
  IRR(51) = RCT(51)*V(88)*V(89)
  IRR(52) = RCT(52)*V(88)*V(89)
  IRR(53) = RCT(53)*V(19)
  IRR(54) = RCT(54)*V(89)*V(89)
  IRR(55) = RCT(55)*V(84)*F(2)
  IRR(56) = RCT(56)*V(14)*V(84)
  IRR(57) = RCT(57)*V(46)*V(84)
  IRR(58) = RCT(58)*V(24)*V(84)
  IRR(59) = RCT(59)*V(30)*V(84)
  IRR(60) = RCT(60)*V(26)*V(84)
  IRR(61) = RCT(61)*V(25)*V(84)
  IRR(62) = RCT(62)*V(22)*V(84)
  IRR(63) = RCT(63)*V(18)*V(84)
  IRR(64) = RCT(64)*V(20)*V(84)
  IRR(65) = RCT(65)*V(38)*V(84)
  IRR(66) = RCT(66)*V(75)*V(84)
  IRR(67) = RCT(67)*V(73)*V(84)
  IRR(68) = RCT(68)*V(51)*V(84)
  IRR(69) = RCT(69)*V(48)*V(84)
  IRR(70) = RCT(70)*V(45)*V(84)
  IRR(71) = RCT(71)*V(47)*V(84)
  IRR(72) = RCT(72)*V(16)*V(84)
  IRR(73) = RCT(73)*V(17)*V(84)
  IRR(74) = RCT(74)*V(50)*V(84)
  IRR(75) = RCT(75)*V(66)*V(84)
  IRR(76) = RCT(76)*V(74)*V(84)
  IRR(77) = RCT(77)*V(59)*V(84)
  IRR(78) = RCT(78)*V(32)*V(84)
  IRR(79) = RCT(79)*V(53)*V(84)
  IRR(80) = RCT(80)*V(62)*V(84)
  IRR(81) = RCT(81)*V(57)*V(84)
  IRR(82) = RCT(82)*V(71)*V(84)
  IRR(83) = RCT(83)*V(23)*V(84)
  IRR(84) = RCT(84)*V(27)*V(84)
  IRR(85) = RCT(85)*V(84)*V(86)
  IRR(86) = RCT(86)*V(40)*V(84)
  IRR(87) = RCT(87)*V(44)*V(84)
  IRR(88) = RCT(88)*V(54)*V(84)
  IRR(89) = RCT(89)*V(70)*V(84)
  IRR(90) = RCT(90)*V(66)*V(89)
  IRR(91) = RCT(91)*V(74)*V(89)
  IRR(92) = RCT(92)*V(53)*V(89)
  IRR(93) = RCT(93)*V(62)*V(89)
  IRR(94) = RCT(94)*V(21)*V(84)
  IRR(95) = RCT(95)*V(71)*V(89)
  IRR(96) = RCT(96)*V(50)*V(89)
  IRR(97) = RCT(97)*V(38)*V(89)
  IRR(98) = RCT(98)*V(75)*V(89)
  IRR(99) = RCT(99)*V(73)*V(89)
  IRR(100) = RCT(100)*V(51)*V(89)
  IRR(101) = RCT(101)*V(48)*V(89)
  IRR(102) = RCT(102)*V(45)*V(89)
  IRR(103) = RCT(103)*V(47)*V(89)
  IRR(104) = RCT(104)*V(54)*V(89)
  IRR(105) = RCT(105)*V(38)*V(79)
  IRR(106) = RCT(106)*V(75)*V(79)
  IRR(107) = RCT(107)*V(73)*V(79)
  IRR(108) = RCT(108)*V(51)*V(79)
  IRR(109) = RCT(109)*V(48)*V(79)
  IRR(110) = RCT(110)*V(45)*V(79)
  IRR(111) = RCT(111)*V(47)*V(79)
  IRR(112) = RCT(112)*V(57)*V(79)
  IRR(113) = RCT(113)*V(71)*V(79)
  IRR(114) = RCT(114)*V(54)*V(79)
  IRR(115) = RCT(115)*V(34)*V(88)
  IRR(116) = RCT(116)*V(34)*V(85)
  IRR(117) = RCT(117)*V(36)*V(88)
  IRR(118) = RCT(118)*V(36)*F(2)
  IRR(119) = RCT(119)*V(36)*V(79)
  IRR(120) = RCT(120)*V(37)*V(88)
  IRR(121) = RCT(121)*V(37)*F(2)
  IRR(122) = RCT(122)*V(37)*V(79)
  IRR(123) = RCT(123)*V(39)*V(88)
  IRR(124) = RCT(124)*V(39)*F(2)
  IRR(125) = RCT(125)*V(39)*V(79)
  IRR(126) = RCT(126)*V(82)*V(88)
  IRR(127) = RCT(127)*V(44)
  IRR(128) = RCT(128)*V(78)*V(88)
  IRR(129) = RCT(129)*V(54)
  IRR(130) = RCT(130)*V(83)*V(87)
  IRR(131) = RCT(131)*V(87)*V(90)
  IRR(132) = RCT(132)*V(80)*V(87)
  IRR(133) = RCT(133)*V(64)*V(87)
  IRR(134) = RCT(134)*V(65)*V(87)
  IRR(135) = RCT(135)*V(55)*V(87)
  IRR(136) = RCT(136)*V(56)*V(87)
  IRR(137) = RCT(137)*V(69)*V(87)
  IRR(138) = RCT(138)*V(60)*V(87)
  IRR(139) = RCT(139)*V(58)*V(87)
  IRR(140) = RCT(140)*V(63)*V(87)
  IRR(141) = RCT(141)*V(67)*V(87)
  IRR(142) = RCT(142)*V(68)*V(87)
  IRR(143) = RCT(143)*V(61)*V(87)
  IRR(144) = RCT(144)*V(82)*V(87)
  IRR(145) = RCT(145)*V(78)*V(87)
  IRR(146) = RCT(146)*V(81)*V(87)
  IRR(147) = RCT(147)*V(77)*V(87)
  IRR(148) = RCT(148)*V(76)*V(87)
  IRR(149) = RCT(149)*V(83)*V(85)
  IRR(150) = RCT(150)*V(85)*V(90)
  IRR(151) = RCT(151)*V(80)*V(85)
  IRR(152) = RCT(152)*V(64)*V(85)
  IRR(153) = RCT(153)*V(65)*V(85)
  IRR(154) = RCT(154)*V(55)*V(85)
  IRR(155) = RCT(155)*V(56)*V(85)
  IRR(156) = RCT(156)*V(69)*V(85)
  IRR(157) = RCT(157)*V(60)*V(85)
  IRR(158) = RCT(158)*V(58)*V(85)
  IRR(159) = RCT(159)*V(63)*V(85)
  IRR(160) = RCT(160)*V(67)*V(85)
  IRR(161) = RCT(161)*V(68)*V(85)
  IRR(162) = RCT(162)*V(61)*V(85)
  IRR(163) = RCT(163)*V(82)*V(85)
  IRR(164) = RCT(164)*V(82)*V(85)
  IRR(165) = RCT(165)*V(78)*V(85)
  IRR(166) = RCT(166)*V(78)*V(85)
  IRR(167) = RCT(167)*V(81)*V(85)
  IRR(168) = RCT(168)*V(77)*V(85)
  IRR(169) = RCT(169)*V(76)*V(85)
  IRR(170) = RCT(170)*V(83)*V(83)
  IRR(171) = RCT(171)*V(83)*V(90)
  IRR(172) = RCT(172)*V(80)*V(83)
  IRR(173) = RCT(173)*V(64)*V(83)
  IRR(174) = RCT(174)*V(65)*V(83)
  IRR(175) = RCT(175)*V(55)*V(83)
  IRR(176) = RCT(176)*V(56)*V(83)
  IRR(177) = RCT(177)*V(69)*V(83)
  IRR(178) = RCT(178)*V(60)*V(83)
  IRR(179) = RCT(179)*V(58)*V(83)
  IRR(180) = RCT(180)*V(63)*V(83)
  IRR(181) = RCT(181)*V(67)*V(83)
  IRR(182) = RCT(182)*V(68)*V(83)
  IRR(183) = RCT(183)*V(61)*V(83)
  IRR(184) = RCT(184)*V(82)*V(83)
  IRR(185) = RCT(185)*V(82)*V(83)
  IRR(186) = RCT(186)*V(78)*V(83)
  IRR(187) = RCT(187)*V(78)*V(83)
  IRR(188) = RCT(188)*V(81)*V(83)
  IRR(189) = RCT(189)*V(77)*V(83)
  IRR(190) = RCT(190)*V(76)*V(83)
  IRR(191) = RCT(191)*V(82)*V(90)
  IRR(192) = RCT(192)*V(80)*V(82)
  IRR(193) = RCT(193)*V(64)*V(82)
  IRR(194) = RCT(194)*V(65)*V(82)
  IRR(195) = RCT(195)*V(55)*V(82)
  IRR(196) = RCT(196)*V(56)*V(82)
  IRR(197) = RCT(197)*V(69)*V(82)
  IRR(198) = RCT(198)*V(60)*V(82)
  IRR(199) = RCT(199)*V(58)*V(82)
  IRR(200) = RCT(200)*V(63)*V(82)
  IRR(201) = RCT(201)*V(67)*V(82)
  IRR(202) = RCT(202)*V(68)*V(82)
  IRR(203) = RCT(203)*V(61)*V(82)
  IRR(204) = RCT(204)*V(82)*V(82)
  IRR(205) = RCT(205)*V(78)*V(82)
  IRR(206) = RCT(206)*V(81)*V(82)
  IRR(207) = RCT(207)*V(77)*V(82)
  IRR(208) = RCT(208)*V(76)*V(82)
  IRR(209) = RCT(209)*V(77)*V(77)
  IRR(210) = RCT(210)*V(76)*V(77)
  IRR(211) = RCT(211)*V(76)*V(76)
  IRR(212) = RCT(212)*V(83)*V(89)
  IRR(213) = RCT(213)*V(89)*V(90)
  IRR(214) = RCT(214)*V(80)*V(89)
  IRR(215) = RCT(215)*V(64)*V(89)
  IRR(216) = RCT(216)*V(65)*V(89)
  IRR(217) = RCT(217)*V(55)*V(89)
  IRR(218) = RCT(218)*V(56)*V(89)
  IRR(219) = RCT(219)*V(69)*V(89)
  IRR(220) = RCT(220)*V(28)*V(84)
  IRR(221) = RCT(221)*V(58)*V(89)
  IRR(222) = RCT(222)*V(63)*V(89)
  IRR(223) = RCT(223)*V(67)*V(89)
  IRR(224) = RCT(224)*V(68)*V(89)
  IRR(225) = RCT(225)*V(61)*V(89)
  IRR(226) = RCT(226)*V(82)*V(89)
  IRR(227) = RCT(227)*V(78)*V(89)
  IRR(228) = RCT(228)*V(81)*V(89)
  IRR(229) = RCT(229)*V(77)*V(89)
  IRR(230) = RCT(230)*V(76)*V(89)
  IRR(231) = RCT(231)*V(72)*V(85)
  IRR(232) = RCT(232)*V(72)*V(83)
  IRR(233) = RCT(233)*V(72)*V(82)
  IRR(234) = RCT(234)*V(72)*V(72)
  IRR(235) = RCT(235)*V(72)*V(87)
  IRR(236) = RCT(236)*V(72)*V(89)
  IRR(237) = RCT(237)*V(60)*V(60)
  IRR(238) = RCT(238)*V(15)*V(84)
  IRR(239) = RCT(239)*V(41)*V(84)
  IRR(240) = RCT(240)*V(52)*V(87)
  IRR(241) = RCT(241)*V(52)*V(85)
  IRR(242) = RCT(242)*V(52)*V(52)
  IRR(243) = RCT(243)*V(52)*V(88)
  IRR(244) = RCT(244)*V(28)
  IRR(245) = RCT(245)*V(42)*V(84)
  IRR(246) = RCT(246)*V(42)*V(79)
  IRR(247) = RCT(247)*V(42)*V(89)
  IRR(248) = RCT(248)*V(49)*V(84)
  IRR(249) = RCT(249)*V(49)*V(89)
  IRR(250) = RCT(250)*V(49)*V(79)
  IRR(251) = RCT(251)*V(8)*V(84)
  IRR(252) = RCT(252)*V(7)*V(84)
  IRR(253) = RCT(253)*V(6)*V(84)
  IRR(254) = RCT(254)*V(12)*V(84)
  IRR(255) = RCT(255)*V(11)*V(84)
  IRR(256) = RCT(256)*V(10)*V(84)
      
END SUBROUTINE racm_soa_vbs_aqchem_IRRFun
















SUBROUTINE racm_soa_vbs_aqchem_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(475)


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

  B(26) = RCT(25)*V(79)

  B(27) = RCT(25)*V(33)

  B(28) = RCT(26)*F(2)

  B(30) = RCT(27)*F(1)

  B(32) = RCT(28)*V(84)

  B(33) = RCT(28)*V(79)

  B(34) = RCT(29)*V(85)

  B(35) = RCT(29)*V(79)

  B(36) = RCT(30)*V(85)

  B(37) = RCT(30)*V(84)

  B(38) = RCT(31)*V(84)

  B(39) = RCT(31)*V(35)

  B(40) = RCT(32)*2*V(85)

  B(41) = RCT(33)*2*V(85)*F(1)

  B(43) = RCT(34)*V(87)

  B(44) = RCT(34)*V(33)

  B(45) = RCT(35)*V(88)

  B(46) = RCT(35)*V(33)

  B(47) = RCT(36)*V(88)

  B(48) = RCT(36)*V(33)

  B(49) = RCT(37)*V(87)

  B(50) = RCT(37)*V(84)

  B(51) = RCT(38)*V(88)

  B(52) = RCT(38)*V(84)

  B(53) = RCT(39)*V(89)

  B(54) = RCT(39)*V(84)

  B(55) = RCT(40)*V(87)

  B(56) = RCT(40)*V(85)

  B(57) = RCT(41)*V(88)

  B(58) = RCT(41)*V(85)

  B(59) = RCT(42)

  B(60) = RCT(43)*V(89)

  B(61) = RCT(43)*V(85)

  B(62) = RCT(44)*V(84)

  B(63) = RCT(44)*V(31)

  B(64) = RCT(45)*V(84)

  B(65) = RCT(45)*V(43)

  B(66) = RCT(46)*V(84)

  B(67) = RCT(46)*V(29)

  B(68) = RCT(47)*V(87)

  B(69) = RCT(47)*V(79)

  B(70) = RCT(48)*V(88)

  B(71) = RCT(48)*V(79)

  B(72) = RCT(49)*2*V(87)*F(2)

  B(74) = RCT(50)*V(89)

  B(75) = RCT(50)*V(87)

  B(76) = RCT(51)*V(89)

  B(77) = RCT(51)*V(88)

  B(78) = RCT(52)*V(89)

  B(79) = RCT(52)*V(88)

  B(80) = RCT(53)

  B(81) = RCT(54)*2*V(89)

  B(82) = RCT(55)*F(2)

  B(84) = RCT(56)*V(84)

  B(85) = RCT(56)*V(14)

  B(86) = RCT(57)*V(84)

  B(87) = RCT(57)*V(46)

  B(88) = RCT(58)*V(84)

  B(89) = RCT(58)*V(24)

  B(90) = RCT(59)*V(84)

  B(91) = RCT(59)*V(30)

  B(92) = RCT(60)*V(84)

  B(93) = RCT(60)*V(26)

  B(94) = RCT(61)*V(84)

  B(95) = RCT(61)*V(25)

  B(96) = RCT(62)*V(84)

  B(97) = RCT(62)*V(22)

  B(98) = RCT(63)*V(84)

  B(99) = RCT(63)*V(18)

  B(100) = RCT(64)*V(84)

  B(101) = RCT(64)*V(20)

  B(102) = RCT(65)*V(84)

  B(103) = RCT(65)*V(38)

  B(104) = RCT(66)*V(84)

  B(105) = RCT(66)*V(75)

  B(106) = RCT(67)*V(84)

  B(107) = RCT(67)*V(73)

  B(108) = RCT(68)*V(84)

  B(109) = RCT(68)*V(51)

  B(110) = RCT(69)*V(84)

  B(111) = RCT(69)*V(48)

  B(112) = RCT(70)*V(84)

  B(113) = RCT(70)*V(45)

  B(114) = RCT(71)*V(84)

  B(115) = RCT(71)*V(47)

  B(116) = RCT(72)*V(84)

  B(117) = RCT(72)*V(16)

  B(118) = RCT(73)*V(84)

  B(119) = RCT(73)*V(17)

  B(120) = RCT(74)*V(84)

  B(121) = RCT(74)*V(50)

  B(122) = RCT(75)*V(84)

  B(123) = RCT(75)*V(66)

  B(124) = RCT(76)*V(84)

  B(125) = RCT(76)*V(74)

  B(126) = RCT(77)*V(84)

  B(127) = RCT(77)*V(59)

  B(128) = RCT(78)*V(84)

  B(129) = RCT(78)*V(32)

  B(130) = RCT(79)*V(84)

  B(131) = RCT(79)*V(53)

  B(132) = RCT(80)*V(84)

  B(133) = RCT(80)*V(62)

  B(134) = RCT(81)*V(84)

  B(135) = RCT(81)*V(57)

  B(136) = RCT(82)*V(84)

  B(137) = RCT(82)*V(71)

  B(138) = RCT(83)*V(84)

  B(139) = RCT(83)*V(23)

  B(140) = RCT(84)*V(84)

  B(141) = RCT(84)*V(27)

  B(142) = RCT(85)*V(86)

  B(143) = RCT(85)*V(84)

  B(144) = RCT(86)*V(84)

  B(145) = RCT(86)*V(40)

  B(146) = RCT(87)*V(84)

  B(147) = RCT(87)*V(44)

  B(148) = RCT(88)*V(84)

  B(149) = RCT(88)*V(54)

  B(150) = RCT(89)*V(84)

  B(151) = RCT(89)*V(70)

  B(152) = RCT(90)*V(89)

  B(153) = RCT(90)*V(66)

  B(154) = RCT(91)*V(89)

  B(155) = RCT(91)*V(74)

  B(156) = RCT(92)*V(89)

  B(157) = RCT(92)*V(53)

  B(158) = RCT(93)*V(89)

  B(159) = RCT(93)*V(62)

  B(160) = RCT(94)*V(84)

  B(161) = RCT(94)*V(21)

  B(162) = RCT(95)*V(89)

  B(163) = RCT(95)*V(71)

  B(164) = RCT(96)*V(89)

  B(165) = RCT(96)*V(50)

  B(166) = RCT(97)*V(89)

  B(167) = RCT(97)*V(38)

  B(168) = RCT(98)*V(89)

  B(169) = RCT(98)*V(75)

  B(170) = RCT(99)*V(89)

  B(171) = RCT(99)*V(73)

  B(172) = RCT(100)*V(89)

  B(173) = RCT(100)*V(51)

  B(174) = RCT(101)*V(89)

  B(175) = RCT(101)*V(48)

  B(176) = RCT(102)*V(89)

  B(177) = RCT(102)*V(45)

  B(178) = RCT(103)*V(89)

  B(179) = RCT(103)*V(47)

  B(180) = RCT(104)*V(89)

  B(181) = RCT(104)*V(54)

  B(182) = RCT(105)*V(79)

  B(183) = RCT(105)*V(38)

  B(184) = RCT(106)*V(79)

  B(185) = RCT(106)*V(75)

  B(186) = RCT(107)*V(79)

  B(187) = RCT(107)*V(73)

  B(188) = RCT(108)*V(79)

  B(189) = RCT(108)*V(51)

  B(190) = RCT(109)*V(79)

  B(191) = RCT(109)*V(48)

  B(192) = RCT(110)*V(79)

  B(193) = RCT(110)*V(45)

  B(194) = RCT(111)*V(79)

  B(195) = RCT(111)*V(47)

  B(196) = RCT(112)*V(79)

  B(197) = RCT(112)*V(57)

  B(198) = RCT(113)*V(79)

  B(199) = RCT(113)*V(71)

  B(200) = RCT(114)*V(79)

  B(201) = RCT(114)*V(54)

  B(202) = RCT(115)*V(88)

  B(203) = RCT(115)*V(34)

  B(204) = RCT(116)*V(85)

  B(205) = RCT(116)*V(34)

  B(206) = RCT(117)*V(88)

  B(207) = RCT(117)*V(36)

  B(208) = RCT(118)*F(2)

  B(210) = RCT(119)*V(79)

  B(211) = RCT(119)*V(36)

  B(212) = RCT(120)*V(88)

  B(213) = RCT(120)*V(37)

  B(214) = RCT(121)*F(2)

  B(216) = RCT(122)*V(79)

  B(217) = RCT(122)*V(37)

  B(218) = RCT(123)*V(88)

  B(219) = RCT(123)*V(39)

  B(220) = RCT(124)*F(2)

  B(222) = RCT(125)*V(79)

  B(223) = RCT(125)*V(39)

  B(224) = RCT(126)*V(88)

  B(225) = RCT(126)*V(82)

  B(226) = RCT(127)

  B(227) = RCT(128)*V(88)

  B(228) = RCT(128)*V(78)

  B(229) = RCT(129)

  B(230) = RCT(130)*V(87)

  B(231) = RCT(130)*V(83)

  B(232) = RCT(131)*V(90)

  B(233) = RCT(131)*V(87)

  B(234) = RCT(132)*V(87)

  B(235) = RCT(132)*V(80)

  B(236) = RCT(133)*V(87)

  B(237) = RCT(133)*V(64)

  B(238) = RCT(134)*V(87)

  B(239) = RCT(134)*V(65)

  B(240) = RCT(135)*V(87)

  B(241) = RCT(135)*V(55)

  B(242) = RCT(136)*V(87)

  B(243) = RCT(136)*V(56)

  B(244) = RCT(137)*V(87)

  B(245) = RCT(137)*V(69)

  B(246) = RCT(138)*V(87)

  B(247) = RCT(138)*V(60)

  B(248) = RCT(139)*V(87)

  B(249) = RCT(139)*V(58)

  B(250) = RCT(140)*V(87)

  B(251) = RCT(140)*V(63)

  B(252) = RCT(141)*V(87)

  B(253) = RCT(141)*V(67)

  B(254) = RCT(142)*V(87)

  B(255) = RCT(142)*V(68)

  B(256) = RCT(143)*V(87)

  B(257) = RCT(143)*V(61)

  B(258) = RCT(144)*V(87)

  B(259) = RCT(144)*V(82)

  B(260) = RCT(145)*V(87)

  B(261) = RCT(145)*V(78)

  B(262) = RCT(146)*V(87)

  B(263) = RCT(146)*V(81)

  B(264) = RCT(147)*V(87)

  B(265) = RCT(147)*V(77)

  B(266) = RCT(148)*V(87)

  B(267) = RCT(148)*V(76)

  B(268) = RCT(149)*V(85)

  B(269) = RCT(149)*V(83)

  B(270) = RCT(150)*V(90)

  B(271) = RCT(150)*V(85)

  B(272) = RCT(151)*V(85)

  B(273) = RCT(151)*V(80)

  B(274) = RCT(152)*V(85)

  B(275) = RCT(152)*V(64)

  B(276) = RCT(153)*V(85)

  B(277) = RCT(153)*V(65)

  B(278) = RCT(154)*V(85)

  B(279) = RCT(154)*V(55)

  B(280) = RCT(155)*V(85)

  B(281) = RCT(155)*V(56)

  B(282) = RCT(156)*V(85)

  B(283) = RCT(156)*V(69)

  B(284) = RCT(157)*V(85)

  B(285) = RCT(157)*V(60)

  B(286) = RCT(158)*V(85)

  B(287) = RCT(158)*V(58)

  B(288) = RCT(159)*V(85)

  B(289) = RCT(159)*V(63)

  B(290) = RCT(160)*V(85)

  B(291) = RCT(160)*V(67)

  B(292) = RCT(161)*V(85)

  B(293) = RCT(161)*V(68)

  B(294) = RCT(162)*V(85)

  B(295) = RCT(162)*V(61)

  B(296) = RCT(163)*V(85)

  B(297) = RCT(163)*V(82)

  B(298) = RCT(164)*V(85)

  B(299) = RCT(164)*V(82)

  B(300) = RCT(165)*V(85)

  B(301) = RCT(165)*V(78)

  B(302) = RCT(166)*V(85)

  B(303) = RCT(166)*V(78)

  B(304) = RCT(167)*V(85)

  B(305) = RCT(167)*V(81)

  B(306) = RCT(168)*V(85)

  B(307) = RCT(168)*V(77)

  B(308) = RCT(169)*V(85)

  B(309) = RCT(169)*V(76)

  B(310) = RCT(170)*2*V(83)

  B(311) = RCT(171)*V(90)

  B(312) = RCT(171)*V(83)

  B(313) = RCT(172)*V(83)

  B(314) = RCT(172)*V(80)

  B(315) = RCT(173)*V(83)

  B(316) = RCT(173)*V(64)

  B(317) = RCT(174)*V(83)

  B(318) = RCT(174)*V(65)

  B(319) = RCT(175)*V(83)

  B(320) = RCT(175)*V(55)

  B(321) = RCT(176)*V(83)

  B(322) = RCT(176)*V(56)

  B(323) = RCT(177)*V(83)

  B(324) = RCT(177)*V(69)

  B(325) = RCT(178)*V(83)

  B(326) = RCT(178)*V(60)

  B(327) = RCT(179)*V(83)

  B(328) = RCT(179)*V(58)

  B(329) = RCT(180)*V(83)

  B(330) = RCT(180)*V(63)

  B(331) = RCT(181)*V(83)

  B(332) = RCT(181)*V(67)

  B(333) = RCT(182)*V(83)

  B(334) = RCT(182)*V(68)

  B(335) = RCT(183)*V(83)

  B(336) = RCT(183)*V(61)

  B(337) = RCT(184)*V(83)

  B(338) = RCT(184)*V(82)

  B(339) = RCT(185)*V(83)

  B(340) = RCT(185)*V(82)

  B(341) = RCT(186)*V(83)

  B(342) = RCT(186)*V(78)

  B(343) = RCT(187)*V(83)

  B(344) = RCT(187)*V(78)

  B(345) = RCT(188)*V(83)

  B(346) = RCT(188)*V(81)

  B(347) = RCT(189)*V(83)

  B(348) = RCT(189)*V(77)

  B(349) = RCT(190)*V(83)

  B(350) = RCT(190)*V(76)

  B(351) = RCT(191)*V(90)

  B(352) = RCT(191)*V(82)

  B(353) = RCT(192)*V(82)

  B(354) = RCT(192)*V(80)

  B(355) = RCT(193)*V(82)

  B(356) = RCT(193)*V(64)

  B(357) = RCT(194)*V(82)

  B(358) = RCT(194)*V(65)

  B(359) = RCT(195)*V(82)

  B(360) = RCT(195)*V(55)

  B(361) = RCT(196)*V(82)

  B(362) = RCT(196)*V(56)

  B(363) = RCT(197)*V(82)

  B(364) = RCT(197)*V(69)

  B(365) = RCT(198)*V(82)

  B(366) = RCT(198)*V(60)

  B(367) = RCT(199)*V(82)

  B(368) = RCT(199)*V(58)

  B(369) = RCT(200)*V(82)

  B(370) = RCT(200)*V(63)

  B(371) = RCT(201)*V(82)

  B(372) = RCT(201)*V(67)

  B(373) = RCT(202)*V(82)

  B(374) = RCT(202)*V(68)

  B(375) = RCT(203)*V(82)

  B(376) = RCT(203)*V(61)

  B(377) = RCT(204)*2*V(82)

  B(378) = RCT(205)*V(82)

  B(379) = RCT(205)*V(78)

  B(380) = RCT(206)*V(82)

  B(381) = RCT(206)*V(81)

  B(382) = RCT(207)*V(82)

  B(383) = RCT(207)*V(77)

  B(384) = RCT(208)*V(82)

  B(385) = RCT(208)*V(76)

  B(386) = RCT(209)*2*V(77)

  B(387) = RCT(210)*V(77)

  B(388) = RCT(210)*V(76)

  B(389) = RCT(211)*2*V(76)

  B(390) = RCT(212)*V(89)

  B(391) = RCT(212)*V(83)

  B(392) = RCT(213)*V(90)

  B(393) = RCT(213)*V(89)

  B(394) = RCT(214)*V(89)

  B(395) = RCT(214)*V(80)

  B(396) = RCT(215)*V(89)

  B(397) = RCT(215)*V(64)

  B(398) = RCT(216)*V(89)

  B(399) = RCT(216)*V(65)

  B(400) = RCT(217)*V(89)

  B(401) = RCT(217)*V(55)

  B(402) = RCT(218)*V(89)

  B(403) = RCT(218)*V(56)

  B(404) = RCT(219)*V(89)

  B(405) = RCT(219)*V(69)

  B(406) = RCT(220)*V(84)

  B(407) = RCT(220)*V(28)

  B(408) = RCT(221)*V(89)

  B(409) = RCT(221)*V(58)

  B(410) = RCT(222)*V(89)

  B(411) = RCT(222)*V(63)

  B(412) = RCT(223)*V(89)

  B(413) = RCT(223)*V(67)

  B(414) = RCT(224)*V(89)

  B(415) = RCT(224)*V(68)

  B(416) = RCT(225)*V(89)

  B(417) = RCT(225)*V(61)

  B(418) = RCT(226)*V(89)

  B(419) = RCT(226)*V(82)

  B(420) = RCT(227)*V(89)

  B(421) = RCT(227)*V(78)

  B(422) = RCT(228)*V(89)

  B(423) = RCT(228)*V(81)

  B(424) = RCT(229)*V(89)

  B(425) = RCT(229)*V(77)

  B(426) = RCT(230)*V(89)

  B(427) = RCT(230)*V(76)

  B(428) = RCT(231)*V(85)

  B(429) = RCT(231)*V(72)

  B(430) = RCT(232)*V(83)

  B(431) = RCT(232)*V(72)

  B(432) = RCT(233)*V(82)

  B(433) = RCT(233)*V(72)

  B(434) = RCT(234)*2*V(72)

  B(435) = RCT(235)*V(87)

  B(436) = RCT(235)*V(72)

  B(437) = RCT(236)*V(89)

  B(438) = RCT(236)*V(72)

  B(439) = RCT(237)*2*V(60)

  B(440) = RCT(238)*V(84)

  B(441) = RCT(238)*V(15)

  B(442) = RCT(239)*V(84)

  B(443) = RCT(239)*V(41)

  B(444) = RCT(240)*V(87)

  B(445) = RCT(240)*V(52)

  B(446) = RCT(241)*V(85)

  B(447) = RCT(241)*V(52)

  B(448) = RCT(242)*2*V(52)

  B(449) = RCT(243)*V(88)

  B(450) = RCT(243)*V(52)

  B(451) = RCT(244)

  B(452) = RCT(245)*V(84)

  B(453) = RCT(245)*V(42)

  B(454) = RCT(246)*V(79)

  B(455) = RCT(246)*V(42)

  B(456) = RCT(247)*V(89)

  B(457) = RCT(247)*V(42)

  B(458) = RCT(248)*V(84)

  B(459) = RCT(248)*V(49)

  B(460) = RCT(249)*V(89)

  B(461) = RCT(249)*V(49)

  B(462) = RCT(250)*V(79)

  B(463) = RCT(250)*V(49)

  B(464) = RCT(251)*V(84)

  B(465) = RCT(251)*V(8)

  B(466) = RCT(252)*V(84)

  B(467) = RCT(252)*V(7)

  B(468) = RCT(253)*V(84)

  B(469) = RCT(253)*V(6)

  B(470) = RCT(254)*V(84)

  B(471) = RCT(254)*V(12)

  B(472) = RCT(255)*V(84)

  B(473) = RCT(255)*V(11)

  B(474) = RCT(256)*V(84)

  B(475) = RCT(256)*V(10)



  JVS(1) = 0

  JVS(2) = B(84)

  JVS(3) = B(85)

  JVS(4) = 0

  JVS(5) = B(86)

  JVS(6) = B(87)

  JVS(7) = 0

  JVS(8) = 0.036*B(96)

  JVS(9) = 0.37*B(182)

  JVS(10) = 0.05*B(452)+0.039*B(454)

  JVS(11) = 0.01*B(194)

  JVS(12) = 0.28*B(190)

  JVS(13) = 0.15*B(188)

  JVS(14) = 0.11*B(200)

  JVS(15) = 0.45*B(196)

  JVS(16) = 0.11*B(198)

  JVS(17) = 0.14*B(184)

  JVS(18) = 0.37*B(183)+0.14*B(185)+0.15*B(189)+0.28*B(191)+0.01*B(195)+0.45*B(197)+0.11*B(199)+0.11*B(201)+0.039*B(455)

  JVS(19) = 0.036*B(97)+0.05*B(453)

  JVS(20) = 0

  JVS(21) = 0.053*B(454)

  JVS(22) = 0.07*B(194)

  JVS(23) = 0.14*B(462)

  JVS(24) = 0.5*B(359)

  JVS(25) = 0.499*B(361)

  JVS(26) = 0.494*B(365)

  JVS(27) = 0.495*B(355)

  JVS(28) = 0.495*B(357)

  JVS(29) = 0.49*B(363)

  JVS(30) = 0.21*B(198)

  JVS(31) = 0.14*B(186)

  JVS(32) = 0.1*B(184)

  JVS(33) = 0.484*B(384)

  JVS(34) = 0.5*B(382)

  JVS(35) = B(302)+B(343)

  JVS(36) = 0.1*B(185)+0.14*B(187)+0.07*B(195)+0.21*B(199)+0.053*B(455)+0.14*B(463)

  JVS(37) = 0.499*B(353)

  JVS(38) = 0.5*B(380)

  JVS(39) = B(298)+B(339)+0.5*B(351)+0.499*B(354)+0.495*B(356)+0.495*B(358)+0.5*B(360)+0.499*B(362)+0.49*B(364)+0.494&
              &*B(366)+0.5*B(381)+0.5*B(383)+0.484*B(385)

  JVS(40) = B(340)+B(344)

  JVS(41) = B(299)+B(303)

  JVS(42) = 0.5*B(352)

  JVS(43) = 0

  JVS(44) = 1.075*B(468)

  JVS(45) = 1.075*B(469)

  JVS(46) = -B(468)

  JVS(47) = 1.075*B(466)

  JVS(48) = 1.075*B(467)-B(469)

  JVS(49) = -B(466)

  JVS(50) = 1.075*B(464)

  JVS(51) = 1.075*B(465)-B(467)

  JVS(52) = -B(464)

  JVS(53) = -B(465)

  JVS(54) = 0

  JVS(55) = 1.075*B(474)

  JVS(56) = 1.075*B(475)

  JVS(57) = -B(474)

  JVS(58) = 1.075*B(472)

  JVS(59) = 1.075*B(473)-B(475)

  JVS(60) = -B(472)

  JVS(61) = 1.075*B(470)

  JVS(62) = 1.075*B(471)-B(473)

  JVS(63) = -B(470)

  JVS(64) = -B(471)

  JVS(65) = -B(28)-B(30)

  JVS(66) = B(2)

  JVS(67) = -B(84)

  JVS(68) = -B(85)

  JVS(69) = -B(440)

  JVS(70) = B(284)

  JVS(71) = -B(441)

  JVS(72) = B(285)

  JVS(73) = -B(116)

  JVS(74) = -B(117)

  JVS(75) = -B(118)

  JVS(76) = -B(119)

  JVS(77) = -B(98)

  JVS(78) = -B(99)

  JVS(79) = -B(80)

  JVS(80) = B(78)

  JVS(81) = B(79)

  JVS(82) = -B(100)

  JVS(83) = -B(101)

  JVS(84) = -B(160)

  JVS(85) = B(446)

  JVS(86) = -B(161)

  JVS(87) = B(447)

  JVS(88) = -B(96)

  JVS(89) = -B(97)

  JVS(90) = -B(138)

  JVS(91) = 0.35*B(136)

  JVS(92) = 0.35*B(137)-B(139)

  JVS(93) = -B(88)

  JVS(94) = B(442)

  JVS(95) = -B(89)+B(443)

  JVS(96) = -B(94)

  JVS(97) = 0.06*B(462)

  JVS(98) = 0.06*B(186)

  JVS(99) = 0.03*B(184)

  JVS(100) = 0.03*B(185)+0.06*B(187)+0.06*B(463)

  JVS(101) = -B(95)

  JVS(102) = -B(92)

  JVS(103) = 0.07*B(462)

  JVS(104) = 0.07*B(186)

  JVS(105) = 0.06*B(184)

  JVS(106) = 0.06*B(185)+0.07*B(187)+0.07*B(463)

  JVS(107) = -B(93)

  JVS(108) = -B(13)-B(140)

  JVS(109) = B(268)

  JVS(110) = -B(141)

  JVS(111) = B(269)

  JVS(112) = -B(406)-B(451)

  JVS(113) = B(449)

  JVS(114) = -B(407)

  JVS(115) = B(450)

  JVS(116) = -B(6)-B(59)-B(66)

  JVS(117) = -B(67)

  JVS(118) = B(57)

  JVS(119) = B(58)

  JVS(120) = B(406)

  JVS(121) = -B(90)

  JVS(122) = B(442)

  JVS(123) = 0.25*B(444)+B(448)

  JVS(124) = -B(91)+B(407)+B(443)

  JVS(125) = 0.25*B(445)

  JVS(126) = 0

  JVS(127) = -B(4)-B(62)

  JVS(128) = B(206)

  JVS(129) = B(212)

  JVS(130) = B(218)

  JVS(131) = B(49)-B(63)

  JVS(132) = B(50)

  JVS(133) = B(207)+B(213)+B(219)

  JVS(134) = 0.024*B(100)

  JVS(135) = -B(23)-B(128)

  JVS(136) = 0.6*B(148)

  JVS(137) = 0.3*B(345)

  JVS(138) = 0.3*B(346)

  JVS(139) = 0.024*B(101)-B(129)+0.6*B(149)

  JVS(140) = B(28)

  JVS(141) = -B(24)-B(26)-B(43)-B(45)-B(47)

  JVS(142) = 0.09*B(188)

  JVS(143) = B(3)-B(27)+0.09*B(189)

  JVS(144) = -B(44)

  JVS(145) = B(1)-B(46)-B(48)

  JVS(146) = B(8)

  JVS(147) = -B(202)-B(204)

  JVS(148) = 0.1*B(120)+B(164)

  JVS(149) = 0.1*B(121)

  JVS(150) = -B(205)

  JVS(151) = -B(203)

  JVS(152) = B(165)

  JVS(153) = -B(9)-B(38)

  JVS(154) = 0.02*B(192)

  JVS(155) = 0.02*B(194)

  JVS(156) = 0.09*B(190)

  JVS(157) = 0.011*B(462)

  JVS(158) = 0.001*B(188)

  JVS(159) = 0.011*B(186)

  JVS(160) = 0.006*B(184)

  JVS(161) = 0.006*B(185)+0.011*B(187)+0.001*B(189)+0.09*B(191)+0.02*B(193)+0.02*B(195)+0.011*B(463)

  JVS(162) = -B(39)

  JVS(163) = B(40)+B(41)

  JVS(164) = 0.9*B(116)

  JVS(165) = -B(206)-B(208)-B(210)

  JVS(166) = -B(211)

  JVS(167) = 0.9*B(117)

  JVS(168) = -B(207)

  JVS(169) = 0.9*B(118)

  JVS(170) = -B(212)-B(214)-B(216)

  JVS(171) = -B(217)

  JVS(172) = 0.9*B(119)

  JVS(173) = -B(213)

  JVS(174) = -B(102)-B(166)-B(182)

  JVS(175) = -B(183)

  JVS(176) = -B(103)

  JVS(177) = -B(167)

  JVS(178) = -B(218)-B(220)-B(222)

  JVS(179) = 0.85*B(120)

  JVS(180) = -B(223)

  JVS(181) = 0.85*B(121)

  JVS(182) = -B(219)

  JVS(183) = -B(15)-B(144)

  JVS(184) = 0.11*B(198)

  JVS(185) = 0.11*B(199)

  JVS(186) = B(296)

  JVS(187) = -B(145)

  JVS(188) = B(297)

  JVS(189) = -B(442)

  JVS(190) = B(174)

  JVS(191) = 0.046*B(246)

  JVS(192) = -B(443)

  JVS(193) = 0.046*B(247)

  JVS(194) = B(175)

  JVS(195) = -B(452)-B(454)-B(456)

  JVS(196) = -B(455)

  JVS(197) = -B(453)

  JVS(198) = -B(457)

  JVS(199) = -B(5)-B(64)

  JVS(200) = B(164)

  JVS(201) = B(156)

  JVS(202) = B(158)

  JVS(203) = B(152)

  JVS(204) = 0.5*B(162)

  JVS(205) = B(154)

  JVS(206) = B(51)-B(65)

  JVS(207) = 0.3*B(60)

  JVS(208) = B(52)

  JVS(209) = 0.3*B(61)+B(153)+B(155)+B(157)+B(159)+0.5*B(163)+B(165)

  JVS(210) = -B(146)-B(226)

  JVS(211) = 0.4*B(148)+0.4*B(180)+0.3*B(200)

  JVS(212) = 0.3*B(201)

  JVS(213) = B(224)

  JVS(214) = -B(147)+0.4*B(149)

  JVS(215) = B(225)

  JVS(216) = 0.4*B(181)

  JVS(217) = -B(112)-B(176)-B(192)

  JVS(218) = -B(193)

  JVS(219) = -B(113)

  JVS(220) = -B(177)

  JVS(221) = 0.036*B(96)

  JVS(222) = B(88)

  JVS(223) = 0.43*B(182)

  JVS(224) = 0

  JVS(225) = 0.14*B(192)

  JVS(226) = -B(86)

  JVS(227) = 0.14*B(194)

  JVS(228) = 0.14*B(190)

  JVS(229) = 0.3*B(462)

  JVS(230) = 0.36*B(188)

  JVS(231) = 0.25*B(444)+0.5*B(448)

  JVS(232) = 1.87*B(17)+1.55*B(18)+2*B(130)+2*B(156)

  JVS(233) = 0.13*B(200)

  JVS(234) = B(22)+0.22*B(196)

  JVS(235) = 0

  JVS(236) = B(19)+B(132)+B(158)

  JVS(237) = B(10)+B(11)+B(122)+B(152)

  JVS(238) = 0.66*B(198)

  JVS(239) = 0.3*B(186)

  JVS(240) = B(12)

  JVS(241) = 0.37*B(184)

  JVS(242) = 0.43*B(183)+0.37*B(185)+0.3*B(187)+0.36*B(189)+0.14*B(191)+0.14*B(193)+0.14*B(195)+0.22*B(197)+0.66*B(199)&
               &+0.13*B(201)+0.3*B(463)

  JVS(243) = -B(87)+B(89)+0.036*B(97)+B(123)+2*B(131)+B(133)

  JVS(244) = 0.25*B(445)

  JVS(245) = B(153)+2*B(157)+B(159)

  JVS(246) = -B(114)-B(178)-B(194)

  JVS(247) = -B(195)

  JVS(248) = -B(115)

  JVS(249) = -B(179)

  JVS(250) = -B(110)-B(174)-B(190)

  JVS(251) = -B(191)

  JVS(252) = -B(111)

  JVS(253) = -B(175)

  JVS(254) = -B(458)-B(460)-B(462)

  JVS(255) = -B(463)

  JVS(256) = -B(459)

  JVS(257) = -B(461)

  JVS(258) = 0.1*B(202)+B(204)

  JVS(259) = B(206)+0.02*B(208)+B(210)

  JVS(260) = B(212)+0.02*B(214)+B(216)

  JVS(261) = B(218)+0.02*B(220)+B(222)

  JVS(262) = -B(120)-B(164)

  JVS(263) = B(211)+B(217)+B(223)

  JVS(264) = -B(121)

  JVS(265) = B(205)

  JVS(266) = 0.1*B(203)+B(207)+B(213)+B(219)

  JVS(267) = -B(165)

  JVS(268) = -B(108)-B(172)-B(188)

  JVS(269) = -B(189)

  JVS(270) = -B(109)

  JVS(271) = -B(173)

  JVS(272) = B(160)

  JVS(273) = B(451)

  JVS(274) = 0.1*B(190)

  JVS(275) = -B(444)-B(446)-2*B(448)-B(449)

  JVS(276) = B(134)

  JVS(277) = 0.1*B(191)

  JVS(278) = B(135)+B(161)

  JVS(279) = -B(447)

  JVS(280) = -B(445)

  JVS(281) = -B(450)

  JVS(282) = 0

  JVS(283) = 0.036*B(96)

  JVS(284) = -B(17)-B(18)-B(130)-B(156)

  JVS(285) = B(256)+B(335)+B(375)+B(416)

  JVS(286) = 1.2*B(252)+0.65*B(331)+0.65*B(371)+1.3*B(412)

  JVS(287) = 0.35*B(254)+0.37*B(333)+0.37*B(373)+0.74*B(414)

  JVS(288) = 0.15*B(136)+0.25*B(162)+0.5*B(198)

  JVS(289) = 0.5*B(199)

  JVS(290) = 0.063*B(234)+0.119*B(313)+0.1*B(353)+0.063*B(394)

  JVS(291) = 0.1*B(354)+0.65*B(372)+0.37*B(374)+B(376)

  JVS(292) = 0.119*B(314)+0.65*B(332)+0.37*B(334)+B(336)

  JVS(293) = 0.036*B(97)-B(131)+0.15*B(137)

  JVS(294) = 0.063*B(235)+1.2*B(253)+0.35*B(255)+B(257)

  JVS(295) = -B(157)+0.25*B(163)+0.063*B(395)+1.3*B(413)+0.74*B(415)+B(417)

  JVS(296) = -B(148)-B(180)-B(200)-B(229)

  JVS(297) = B(227)

  JVS(298) = -B(201)

  JVS(299) = -B(149)

  JVS(300) = B(228)

  JVS(301) = -B(181)

  JVS(302) = B(102)

  JVS(303) = -B(240)-B(278)-B(319)-B(359)-B(400)

  JVS(304) = 0

  JVS(305) = -B(360)

  JVS(306) = -B(320)

  JVS(307) = B(103)

  JVS(308) = -B(279)

  JVS(309) = -B(241)

  JVS(310) = -B(401)

  JVS(311) = -B(242)-B(280)-B(321)-B(361)-B(402)

  JVS(312) = B(104)

  JVS(313) = -B(362)

  JVS(314) = -B(322)

  JVS(315) = B(105)

  JVS(316) = -B(281)

  JVS(317) = -B(243)

  JVS(318) = -B(403)

  JVS(319) = B(440)

  JVS(320) = 0.9*B(456)

  JVS(321) = 0.79*B(194)

  JVS(322) = 0.65*B(190)

  JVS(323) = 0.9*B(172)+0.39*B(188)

  JVS(324) = -B(22)-B(134)-B(196)

  JVS(325) = B(246)+0.55*B(325)+0.771*B(365)+2*B(439)

  JVS(326) = 0.4*B(250)+0.6*B(329)+0.6*B(369)+0.6*B(410)

  JVS(327) = 0.39*B(189)+0.65*B(191)+0.79*B(195)-B(197)

  JVS(328) = 0.771*B(366)+0.6*B(370)

  JVS(329) = 0.55*B(326)+0.6*B(330)

  JVS(330) = -B(135)+B(441)

  JVS(331) = 0

  JVS(332) = B(247)+0.4*B(251)

  JVS(333) = 0.9*B(173)+0.6*B(411)+0.9*B(457)

  JVS(334) = B(112)

  JVS(335) = -B(248)-B(286)-B(327)-B(367)-B(408)

  JVS(336) = 0

  JVS(337) = -B(368)

  JVS(338) = -B(328)

  JVS(339) = B(113)

  JVS(340) = -B(287)

  JVS(341) = -B(249)

  JVS(342) = -B(409)

  JVS(343) = 0.25*B(98)

  JVS(344) = 0.12*B(138)

  JVS(345) = 0.36*B(452)+0.23*B(454)

  JVS(346) = 0.53*B(192)

  JVS(347) = 0.16*B(462)

  JVS(348) = 0.06*B(242)+0.081*B(321)+0.141*B(361)+0.06*B(402)

  JVS(349) = 0.8*B(248)+B(327)+B(367)+B(408)

  JVS(350) = -B(16)-B(126)

  JVS(351) = 0.722*B(236)+0.24*B(315)+0.33*B(355)+0.828*B(396)

  JVS(352) = 0.642*B(238)+0.419*B(317)+0.581*B(357)+0.88*B(398)

  JVS(353) = 0.29*B(244)+0.313*B(323)+0.569*B(363)+0.29*B(404)

  JVS(354) = 0.8*B(21)

  JVS(355) = 0.03*B(162)

  JVS(356) = 0.16*B(186)

  JVS(357) = 0.03*B(184)

  JVS(358) = 0.464*B(266)+0.149*B(349)+0.167*B(384)+0.149*B(387)+0.285*B(389)+0.469*B(426)

  JVS(359) = 0.149*B(388)

  JVS(360) = 0.03*B(185)+0.16*B(187)+0.53*B(193)+0.23*B(455)+0.16*B(463)

  JVS(361) = 0.623*B(234)+0.018*B(313)+0.127*B(353)+0.67*B(394)

  JVS(362) = 0.11*B(380)

  JVS(363) = 0.127*B(354)+0.33*B(356)+0.581*B(358)+0.141*B(362)+0.569*B(364)+B(368)+0.11*B(381)+0.167*B(385)

  JVS(364) = 0.018*B(314)+0.24*B(316)+0.419*B(318)+0.081*B(322)+0.313*B(324)+B(328)+0.149*B(350)

  JVS(365) = 0.25*B(99)-B(127)+0.12*B(139)+0.41*B(142)+0.36*B(453)

  JVS(366) = 0

  JVS(367) = 0.41*B(143)

  JVS(368) = 0.623*B(235)+0.722*B(237)+0.642*B(239)+0.06*B(243)+0.29*B(245)+0.8*B(249)+0.464*B(267)

  JVS(369) = 0.03*B(163)+0.67*B(395)+0.828*B(397)+0.88*B(399)+0.06*B(403)+0.29*B(405)+B(409)+0.469*B(427)

  JVS(370) = B(110)

  JVS(371) = B(108)

  JVS(372) = -B(246)-B(284)-B(325)-B(365)-2*B(439)

  JVS(373) = 0

  JVS(374) = -B(366)

  JVS(375) = -B(326)

  JVS(376) = B(109)+B(111)

  JVS(377) = -B(285)

  JVS(378) = -B(247)

  JVS(379) = 0

  JVS(380) = 0.98*B(220)

  JVS(381) = 0

  JVS(382) = -B(256)-B(294)-B(335)-B(375)-B(416)

  JVS(383) = 0

  JVS(384) = -B(376)

  JVS(385) = -B(336)

  JVS(386) = 0

  JVS(387) = -B(295)

  JVS(388) = -B(257)

  JVS(389) = 0

  JVS(390) = -B(417)

  JVS(391) = B(90)

  JVS(392) = B(128)

  JVS(393) = 0

  JVS(394) = 0

  JVS(395) = 0.5*B(444)+B(448)

  JVS(396) = 0

  JVS(397) = 0.9*B(196)

  JVS(398) = 0

  JVS(399) = B(256)+B(335)+B(375)+B(416)

  JVS(400) = -B(19)-B(132)-B(158)

  JVS(401) = 0

  JVS(402) = 0.65*B(252)+0.35*B(331)+0.35*B(371)+0.7*B(412)

  JVS(403) = 0.6*B(254)+0.63*B(333)+0.63*B(373)+1.26*B(414)

  JVS(404) = 0.15*B(136)+0.25*B(162)+0.62*B(198)

  JVS(405) = 0

  JVS(406) = 0.9*B(197)+0.62*B(199)

  JVS(407) = 0.005*B(313)+0.004*B(353)

  JVS(408) = 0.54*B(262)+0.4*B(345)+0.54*B(380)+0.54*B(422)

  JVS(409) = 0.004*B(354)+0.35*B(372)+0.63*B(374)+B(376)+0.54*B(381)

  JVS(410) = 0.005*B(314)+0.35*B(332)+0.63*B(334)+B(336)+0.4*B(346)

  JVS(411) = B(91)+B(129)-B(133)+0.15*B(137)

  JVS(412) = 0

  JVS(413) = 0.65*B(253)+0.6*B(255)+B(257)+0.54*B(263)+0.5*B(445)

  JVS(414) = 0

  JVS(415) = -B(159)+0.25*B(163)+0.7*B(413)+1.26*B(415)+B(417)+0.54*B(423)

  JVS(416) = B(114)

  JVS(417) = -B(250)-B(288)-B(329)-B(369)-B(410)

  JVS(418) = 0

  JVS(419) = -B(370)

  JVS(420) = -B(330)

  JVS(421) = B(115)

  JVS(422) = -B(289)

  JVS(423) = -B(251)

  JVS(424) = -B(411)

  JVS(425) = 0.75*B(98)

  JVS(426) = -B(236)-B(274)-B(315)-B(355)-B(396)

  JVS(427) = -B(356)

  JVS(428) = -B(316)

  JVS(429) = 0.75*B(99)

  JVS(430) = -B(275)

  JVS(431) = -B(237)

  JVS(432) = -B(397)

  JVS(433) = 0.9511*B(100)

  JVS(434) = -B(238)-B(276)-B(317)-B(357)-B(398)

  JVS(435) = -B(358)

  JVS(436) = -B(318)

  JVS(437) = 0.9511*B(101)

  JVS(438) = -B(277)

  JVS(439) = -B(239)

  JVS(440) = -B(399)

  JVS(441) = 0.01*B(96)

  JVS(442) = B(88)

  JVS(443) = B(13)+0.35*B(140)

  JVS(444) = B(23)

  JVS(445) = B(182)

  JVS(446) = 0.35*B(144)

  JVS(447) = 0

  JVS(448) = 0.3*B(452)+0.51*B(454)

  JVS(449) = B(146)

  JVS(450) = 0.04*B(194)

  JVS(451) = 0.58*B(190)

  JVS(452) = 0.02*B(462)

  JVS(453) = 0.9*B(188)

  JVS(454) = 0.75*B(444)+0.5*B(448)

  JVS(455) = 0.13*B(17)+0.45*B(18)

  JVS(456) = 0.4*B(148)+0.4*B(180)+0.7*B(200)

  JVS(457) = 1.6*B(240)+1.55*B(319)+0.8*B(359)+1.6*B(400)

  JVS(458) = B(242)+1.25*B(321)+0.501*B(361)+B(402)

  JVS(459) = B(22)

  JVS(460) = B(327)

  JVS(461) = B(246)+1.09*B(325)+0.34*B(365)+B(439)

  JVS(462) = B(335)

  JVS(463) = 0.25*B(250)+1.4*B(329)+0.4*B(369)+0.4*B(410)

  JVS(464) = 0.021*B(236)+0.829*B(315)+0.076*B(355)+0.021*B(396)

  JVS(465) = 0.753*B(317)

  JVS(466) = -B(10)-B(11)-B(122)-B(152)

  JVS(467) = B(331)

  JVS(468) = B(333)

  JVS(469) = 0.755*B(323)

  JVS(470) = 0

  JVS(471) = B(430)

  JVS(472) = 0.02*B(186)

  JVS(473) = 0.64*B(184)

  JVS(474) = 0.287*B(266)+0.96*B(349)+0.207*B(384)+0.202*B(387)+0.504*B(389)+0.28*B(426)

  JVS(475) = 0.75*B(347)+0.202*B(388)

  JVS(476) = B(260)+2*B(341)+B(343)+B(378)+B(420)

  JVS(477) = B(183)+0.64*B(185)+0.02*B(187)+0.9*B(189)+0.58*B(191)+0.04*B(195)+0.7*B(201)+0.51*B(455)+0.02*B(463)

  JVS(478) = 0.047*B(234)+0.81*B(313)+0.091*B(353)+0.048*B(394)

  JVS(479) = 0.75*B(345)

  JVS(480) = B(337)+B(339)+0.091*B(354)+0.076*B(356)+0.8*B(360)+0.501*B(362)+0.34*B(366)+0.4*B(370)+B(379)+0.207*B(385)

  JVS(481) = B(230)+1.33*B(310)+0.75*B(311)+0.81*B(314)+0.829*B(316)+0.753*B(318)+1.55*B(320)+1.25*B(322)+0.755*B(324)&
               &+1.09*B(326)+B(328)+1.4*B(330)+B(332)+B(334)+B(336)+B(338)+B(340)+2*B(342)+B(344)+0.75*B(346)+0.75*B(348)&
               &+0.96*B(350)+B(390)+B(431)

  JVS(482) = B(89)+0.01*B(97)-B(123)+0.35*B(141)+0.35*B(145)+B(147)+0.4*B(149)+0.3*B(453)

  JVS(483) = 0

  JVS(484) = B(231)+0.047*B(235)+0.021*B(237)+1.6*B(241)+B(243)+B(247)+0.25*B(251)+B(261)+0.287*B(267)+0.75*B(445)

  JVS(485) = 0

  JVS(486) = -B(153)+0.4*B(181)+B(391)+0.048*B(395)+0.021*B(397)+1.6*B(401)+B(403)+0.4*B(411)+B(421)+0.28*B(427)

  JVS(487) = 0.75*B(312)

  JVS(488) = 0.98*B(208)

  JVS(489) = -B(252)-B(290)-B(331)-B(371)-B(412)

  JVS(490) = 0

  JVS(491) = -B(372)

  JVS(492) = -B(332)

  JVS(493) = 0

  JVS(494) = -B(291)

  JVS(495) = -B(253)

  JVS(496) = 0

  JVS(497) = -B(413)

  JVS(498) = 0.98*B(214)

  JVS(499) = -B(254)-B(292)-B(333)-B(373)-B(414)

  JVS(500) = 0

  JVS(501) = -B(374)

  JVS(502) = -B(334)

  JVS(503) = 0

  JVS(504) = -B(293)

  JVS(505) = -B(255)

  JVS(506) = 0

  JVS(507) = -B(415)

  JVS(508) = 0.19*B(452)

  JVS(509) = B(458)

  JVS(510) = -B(244)-B(282)-B(323)-B(363)-B(404)

  JVS(511) = B(106)

  JVS(512) = 0

  JVS(513) = -B(364)

  JVS(514) = -B(324)

  JVS(515) = B(107)+0.19*B(453)+B(459)

  JVS(516) = -B(283)

  JVS(517) = -B(245)

  JVS(518) = -B(405)

  JVS(519) = B(202)

  JVS(520) = 0

  JVS(521) = 0.6*B(180)

  JVS(522) = 0.2*B(248)

  JVS(523) = 0.35*B(250)

  JVS(524) = 0.124*B(236)

  JVS(525) = 0.261*B(238)

  JVS(526) = 0.05*B(252)

  JVS(527) = 0.05*B(254)

  JVS(528) = -B(21)-B(150)

  JVS(529) = B(308)+0.5*B(349)+0.484*B(384)+1.5*B(387)+B(389)

  JVS(530) = B(264)+B(306)+B(347)+B(382)+2*B(386)+1.5*B(388)+B(424)

  JVS(531) = 0

  JVS(532) = 0

  JVS(533) = 0.059*B(234)

  JVS(534) = B(383)+0.484*B(385)

  JVS(535) = B(348)+0.5*B(350)

  JVS(536) = -B(151)

  JVS(537) = B(307)+B(309)

  JVS(538) = 0.059*B(235)+0.124*B(237)+0.261*B(239)+0.2*B(249)+0.35*B(251)+0.05*B(253)+0.05*B(255)+B(265)

  JVS(539) = B(203)

  JVS(540) = 0.6*B(181)+B(425)

  JVS(541) = 0.5*B(252)+B(331)+B(371)+0.5*B(412)

  JVS(542) = 0.95*B(254)+B(333)+B(373)+B(414)

  JVS(543) = -B(20)-B(136)-B(162)-B(198)

  JVS(544) = -B(199)

  JVS(545) = B(372)+B(374)

  JVS(546) = B(332)+B(334)

  JVS(547) = -B(137)

  JVS(548) = 0

  JVS(549) = 0.5*B(253)+0.95*B(255)

  JVS(550) = 0

  JVS(551) = -B(163)+0.5*B(413)+B(415)

  JVS(552) = 0.1*B(116)

  JVS(553) = 0.1*B(118)

  JVS(554) = 0.35*B(144)

  JVS(555) = B(146)

  JVS(556) = 0.05*B(120)

  JVS(557) = 0.13*B(188)

  JVS(558) = B(148)+B(180)

  JVS(559) = 0.334*B(236)+0.245*B(315)+0.237*B(355)+0.391*B(396)

  JVS(560) = 0.416*B(238)+0.322*B(317)+0.318*B(357)+0.587*B(398)

  JVS(561) = 0.5*B(136)+0.5*B(162)

  JVS(562) = -B(428)-B(430)-B(432)-2*B(434)-B(435)-B(437)

  JVS(563) = 0

  JVS(564) = 0.13*B(189)

  JVS(565) = 0.048*B(234)+0.085*B(313)+0.071*B(353)+0.051*B(394)

  JVS(566) = 0.16*B(262)+0.08*B(345)+0.08*B(380)+0.16*B(422)

  JVS(567) = 0.071*B(354)+0.237*B(356)+0.318*B(358)+0.08*B(381)-B(433)

  JVS(568) = 0.085*B(314)+0.245*B(316)+0.322*B(318)+0.08*B(346)-B(431)

  JVS(569) = 0.1*B(117)+0.1*B(119)+0.05*B(121)+0.5*B(137)+0.07*B(142)+0.35*B(145)+B(147)+B(149)

  JVS(570) = -B(429)

  JVS(571) = 0.07*B(143)

  JVS(572) = 0.048*B(235)+0.334*B(237)+0.416*B(239)+0.16*B(263)-B(436)

  JVS(573) = 0

  JVS(574) = 0.5*B(163)+B(181)+0.051*B(395)+0.391*B(397)+0.587*B(399)+0.16*B(423)-B(438)

  JVS(575) = 0.08*B(325)

  JVS(576) = 0.25*B(250)+0.4*B(329)+0.4*B(369)+0.4*B(410)

  JVS(577) = -B(106)-B(170)-B(186)

  JVS(578) = -B(187)

  JVS(579) = 0.4*B(370)

  JVS(580) = 0.08*B(326)+0.4*B(330)

  JVS(581) = -B(107)

  JVS(582) = 0

  JVS(583) = 0.25*B(251)

  JVS(584) = -B(171)+0.4*B(411)

  JVS(585) = 0.025*B(100)

  JVS(586) = 0.335*B(96)

  JVS(587) = 0.88*B(138)

  JVS(588) = 0.85*B(454)

  JVS(589) = 0.65*B(192)

  JVS(590) = 0.99*B(462)

  JVS(591) = 0.2*B(240)+0.35*B(319)+0.6*B(359)+0.2*B(400)

  JVS(592) = 0.94*B(242)+0.669*B(321)+0.859*B(361)+0.94*B(402)

  JVS(593) = 0.8*B(248)+B(327)+B(367)+B(408)

  JVS(594) = 0.211*B(236)+0.523*B(315)+0.677*B(355)+0.239*B(396)

  JVS(595) = 0.15*B(238)+0.411*B(317)+0.497*B(357)+0.187*B(398)

  JVS(596) = 1.71*B(244)+0.932*B(323)+0.941*B(363)+1.71*B(404)

  JVS(597) = 0.2*B(21)

  JVS(598) = 0.25*B(162)+0.16*B(198)

  JVS(599) = 0.99*B(186)

  JVS(600) = -B(12)-B(124)-B(154)

  JVS(601) = 0.44*B(184)

  JVS(602) = 1.24*B(266)+0.64*B(349)+0.65*B(384)+0.64*B(387)+1.21*B(389)+1.24*B(426)

  JVS(603) = 0.64*B(388)

  JVS(604) = 0

  JVS(605) = 0.44*B(185)+0.99*B(187)+0.65*B(193)+0.16*B(199)+0.85*B(455)+0.99*B(463)

  JVS(606) = 0.233*B(234)+0.58*B(313)+0.724*B(353)+0.243*B(394)

  JVS(607) = 0.46*B(262)+0.3*B(345)+0.35*B(380)+0.46*B(422)

  JVS(608) = B(351)+0.724*B(354)+0.677*B(356)+0.497*B(358)+0.6*B(360)+0.859*B(362)+0.941*B(364)+B(368)+0.35*B(381)+0.65&
               &*B(385)

  JVS(609) = 0.75*B(311)+0.58*B(314)+0.523*B(316)+0.411*B(318)+0.35*B(320)+0.669*B(322)+0.932*B(324)+B(328)+0.3*B(346)&
               &+0.64*B(350)

  JVS(610) = 0.335*B(97)+0.025*B(101)-B(125)+0.88*B(139)+0.08*B(142)

  JVS(611) = 0

  JVS(612) = B(14)+0.08*B(143)

  JVS(613) = B(232)+0.233*B(235)+0.211*B(237)+0.15*B(239)+0.2*B(241)+0.94*B(243)+1.71*B(245)+0.8*B(249)+0.46*B(263)+1.24&
               &*B(267)

  JVS(614) = 0

  JVS(615) = -B(155)+0.25*B(163)+B(392)+0.243*B(395)+0.239*B(397)+0.187*B(399)+0.2*B(401)+0.94*B(403)+1.71*B(405)+B(409)&
               &+0.46*B(423)+1.24*B(427)

  JVS(616) = B(233)+0.75*B(312)+B(352)+B(393)

  JVS(617) = 0.46*B(194)

  JVS(618) = 0.35*B(188)

  JVS(619) = 0.37*B(325)+0.229*B(365)

  JVS(620) = -B(104)-B(168)-B(184)

  JVS(621) = -B(185)+0.35*B(189)+0.46*B(195)

  JVS(622) = 0.229*B(366)

  JVS(623) = 0.37*B(326)

  JVS(624) = -B(105)

  JVS(625) = 0

  JVS(626) = 0

  JVS(627) = -B(169)

  JVS(628) = 0.2*B(166)

  JVS(629) = 0.1*B(456)

  JVS(630) = 0.9*B(176)

  JVS(631) = 0.87*B(178)

  JVS(632) = 0.89*B(460)

  JVS(633) = 0.1*B(172)

  JVS(634) = 0.89*B(170)

  JVS(635) = 0.57*B(168)

  JVS(636) = -B(266)-B(308)-B(349)-B(384)-B(387)-2*B(389)-B(426)

  JVS(637) = -B(388)

  JVS(638) = 0

  JVS(639) = -B(385)

  JVS(640) = -B(350)

  JVS(641) = 0

  JVS(642) = -B(309)

  JVS(643) = -B(267)

  JVS(644) = 0.2*B(167)+0.57*B(169)+0.89*B(171)+0.1*B(173)+0.9*B(177)+0.87*B(179)-B(427)+0.1*B(457)+0.89*B(461)

  JVS(645) = 0.8*B(166)

  JVS(646) = 0.9*B(456)

  JVS(647) = 0.1*B(176)

  JVS(648) = 0.13*B(178)

  JVS(649) = 0.11*B(460)

  JVS(650) = 0.9*B(172)

  JVS(651) = 0.11*B(170)

  JVS(652) = 0.43*B(168)

  JVS(653) = -B(387)

  JVS(654) = -B(264)-B(306)-B(347)-B(382)-2*B(386)-B(388)-B(424)

  JVS(655) = 0

  JVS(656) = -B(383)

  JVS(657) = -B(348)

  JVS(658) = 0

  JVS(659) = -B(307)

  JVS(660) = -B(265)

  JVS(661) = 0.8*B(167)+0.43*B(169)+0.11*B(171)+0.9*B(173)+0.1*B(177)+0.13*B(179)-B(425)+0.9*B(457)+0.11*B(461)

  JVS(662) = B(229)

  JVS(663) = B(20)+0.5*B(136)+0.5*B(162)

  JVS(664) = -B(227)-B(260)-B(300)-B(302)-B(341)-B(343)-B(378)-B(420)

  JVS(665) = 0

  JVS(666) = -B(379)

  JVS(667) = -B(342)-B(344)

  JVS(668) = 0.5*B(137)

  JVS(669) = -B(301)-B(303)

  JVS(670) = -B(261)

  JVS(671) = -B(228)

  JVS(672) = 0.5*B(163)-B(421)

  JVS(673) = B(24)-B(26)

  JVS(674) = -B(210)

  JVS(675) = -B(216)

  JVS(676) = -B(182)

  JVS(677) = -B(222)

  JVS(678) = -B(454)

  JVS(679) = -B(192)

  JVS(680) = -B(194)

  JVS(681) = -B(190)

  JVS(682) = -B(462)

  JVS(683) = 0

  JVS(684) = -B(188)

  JVS(685) = -B(200)

  JVS(686) = -B(196)

  JVS(687) = 0

  JVS(688) = 0

  JVS(689) = -B(198)

  JVS(690) = -B(186)

  JVS(691) = -B(184)

  JVS(692) = B(302)

  JVS(693) = -B(2)-B(3)-B(27)-B(32)-B(34)-B(68)-B(70)-B(183)-B(185)-B(187)-B(189)-B(191)-B(193)-B(195)-B(197)-B(199)&
               &-B(201)-B(211)-B(217)-B(223)-B(455)-B(463)

  JVS(694) = B(298)

  JVS(695) = 0

  JVS(696) = -B(33)

  JVS(697) = -B(35)+B(299)+B(303)

  JVS(698) = -B(69)

  JVS(699) = -B(71)

  JVS(700) = 0

  JVS(701) = 0.583*B(96)

  JVS(702) = B(150)

  JVS(703) = 0

  JVS(704) = 0

  JVS(705) = 0

  JVS(706) = 0

  JVS(707) = -B(234)-B(272)-B(313)-B(353)-B(394)

  JVS(708) = -B(354)

  JVS(709) = -B(314)

  JVS(710) = 0.583*B(97)+0.44*B(142)+B(151)

  JVS(711) = -B(273)

  JVS(712) = 0.44*B(143)

  JVS(713) = -B(235)

  JVS(714) = 0

  JVS(715) = -B(395)

  JVS(716) = 0.42*B(192)

  JVS(717) = 0.42*B(194)

  JVS(718) = 0.12*B(462)

  JVS(719) = 0.02*B(188)

  JVS(720) = B(126)

  JVS(721) = 0

  JVS(722) = 0

  JVS(723) = 0

  JVS(724) = 0

  JVS(725) = 0

  JVS(726) = 0.12*B(186)

  JVS(727) = 0.03*B(184)

  JVS(728) = 0

  JVS(729) = 0

  JVS(730) = 0

  JVS(731) = 0.03*B(185)+0.12*B(187)+0.02*B(189)+0.42*B(193)+0.42*B(195)+0.12*B(463)

  JVS(732) = 0

  JVS(733) = -B(262)-B(304)-B(345)-B(380)-B(422)

  JVS(734) = -B(381)

  JVS(735) = -B(346)

  JVS(736) = B(127)

  JVS(737) = -B(305)

  JVS(738) = 0

  JVS(739) = -B(263)

  JVS(740) = 0

  JVS(741) = -B(423)

  JVS(742) = B(23)

  JVS(743) = 0.65*B(144)

  JVS(744) = B(226)

  JVS(745) = 0.1*B(190)

  JVS(746) = 0.15*B(188)

  JVS(747) = 0.25*B(444)

  JVS(748) = 0.7*B(200)

  JVS(749) = -B(359)

  JVS(750) = -B(361)

  JVS(751) = B(22)+0.1*B(196)

  JVS(752) = -B(367)

  JVS(753) = B(16)

  JVS(754) = -B(365)

  JVS(755) = -B(375)

  JVS(756) = B(19)+B(132)+B(158)

  JVS(757) = -B(369)

  JVS(758) = -B(355)

  JVS(759) = -B(357)

  JVS(760) = -B(371)

  JVS(761) = -B(373)

  JVS(762) = -B(363)

  JVS(763) = 0

  JVS(764) = 0.28*B(198)

  JVS(765) = -B(432)

  JVS(766) = 0

  JVS(767) = B(124)+B(154)

  JVS(768) = 0

  JVS(769) = -B(384)

  JVS(770) = -B(382)

  JVS(771) = B(260)+B(341)+B(420)

  JVS(772) = 0.15*B(189)+0.1*B(191)+0.1*B(197)+0.28*B(199)+0.7*B(201)

  JVS(773) = -B(353)

  JVS(774) = 0.23*B(262)+0.12*B(345)-0.88*B(380)+0.23*B(422)

  JVS(775) = -B(224)-B(258)-B(296)-B(298)-B(337)-B(339)-B(351)-B(354)-B(356)-B(358)-B(360)-B(362)-B(364)-B(366)-B(368)&
               &-B(370)-B(372)-B(374)-B(376)-2*B(377)-0.88*B(381)-B(383)-B(385)-B(418)-B(433)

  JVS(776) = -B(338)-B(340)+B(342)+0.12*B(346)

  JVS(777) = B(125)+B(133)+0.65*B(145)

  JVS(778) = -B(297)-B(299)

  JVS(779) = 0

  JVS(780) = -B(259)+B(261)+0.23*B(263)+0.25*B(445)

  JVS(781) = -B(225)

  JVS(782) = B(155)+B(159)-B(419)+B(421)+0.23*B(423)

  JVS(783) = -B(352)

  JVS(784) = B(92)

  JVS(785) = 0.65*B(140)

  JVS(786) = B(15)

  JVS(787) = 0.08*B(190)

  JVS(788) = 0.23*B(462)

  JVS(789) = 0.03*B(188)

  JVS(790) = -B(319)+0.5*B(359)

  JVS(791) = -B(321)+0.501*B(361)

  JVS(792) = -B(327)+B(367)

  JVS(793) = -B(325)+0.506*B(365)

  JVS(794) = -B(335)+B(375)

  JVS(795) = -B(329)+B(369)

  JVS(796) = 0.031*B(236)-0.951*B(315)+0.554*B(355)+0.04*B(396)

  JVS(797) = -B(317)+0.507*B(357)

  JVS(798) = -B(331)+B(371)

  JVS(799) = -B(333)+B(373)

  JVS(800) = -B(323)+0.51*B(363)

  JVS(801) = 0

  JVS(802) = -B(430)+B(432)

  JVS(803) = 0.23*B(186)

  JVS(804) = B(12)

  JVS(805) = 0.19*B(184)

  JVS(806) = -B(349)+0.516*B(384)

  JVS(807) = -B(347)+0.5*B(382)

  JVS(808) = -B(341)-B(343)+B(378)

  JVS(809) = 0.19*B(185)+0.23*B(187)+0.03*B(189)+0.08*B(191)+0.23*B(463)

  JVS(810) = 0.15*B(234)-0.993*B(313)+0.508*B(353)+0.155*B(394)

  JVS(811) = -B(345)+0.5*B(380)

  JVS(812) = B(258)-B(339)+0.5*B(351)+0.508*B(354)+0.554*B(356)+0.507*B(358)+0.5*B(360)+0.501*B(362)+0.51*B(364)+0.506&
               &*B(366)+B(368)+B(370)+B(372)+B(374)+B(376)+2*B(377)+B(379)+0.5*B(381)+0.5*B(383)+0.516*B(385)+B(418)+B(433)

  JVS(813) = -B(230)-B(268)-2*B(310)-B(311)-0.993*B(314)-0.951*B(316)-B(318)-B(320)-B(322)-B(324)-B(326)-B(328)-B(330)&
               &-B(332)-B(334)-B(336)-B(340)-B(342)-B(344)-B(346)-B(348)-B(350)-B(390)-B(431)

  JVS(814) = B(93)+0.65*B(141)

  JVS(815) = -B(269)

  JVS(816) = 0

  JVS(817) = -B(231)+0.15*B(235)+0.031*B(237)+B(259)

  JVS(818) = 0

  JVS(819) = -B(391)+0.155*B(395)+0.04*B(397)+B(419)

  JVS(820) = -B(312)+0.5*B(352)

  JVS(821) = 2*B(30)

  JVS(822) = -B(84)

  JVS(823) = -B(116)

  JVS(824) = -B(118)

  JVS(825) = -B(98)

  JVS(826) = -B(100)

  JVS(827) = -B(160)

  JVS(828) = -0.964*B(96)

  JVS(829) = -B(138)

  JVS(830) = -B(88)

  JVS(831) = -B(94)

  JVS(832) = -B(92)

  JVS(833) = B(13)-0.65*B(140)

  JVS(834) = -B(406)

  JVS(835) = 0.35*B(6)-B(66)

  JVS(836) = -B(90)

  JVS(837) = B(4)-B(62)

  JVS(838) = -B(128)

  JVS(839) = 2*B(9)-B(38)

  JVS(840) = B(210)

  JVS(841) = B(216)

  JVS(842) = -B(102)+0.12*B(182)

  JVS(843) = B(222)

  JVS(844) = B(15)-B(144)

  JVS(845) = -B(442)

  JVS(846) = -B(452)+0.63*B(454)

  JVS(847) = B(5)-B(64)

  JVS(848) = -B(146)

  JVS(849) = -B(112)+0.85*B(192)

  JVS(850) = -B(86)

  JVS(851) = -B(114)+0.85*B(194)

  JVS(852) = -B(110)+0.25*B(190)

  JVS(853) = -B(458)+0.63*B(462)

  JVS(854) = -B(120)

  JVS(855) = -B(108)+0.28*B(188)

  JVS(856) = 0

  JVS(857) = -B(130)

  JVS(858) = -B(148)+0.036*B(200)

  JVS(859) = -B(134)+0.19*B(196)

  JVS(860) = -B(126)

  JVS(861) = 0

  JVS(862) = 0

  JVS(863) = -B(132)

  JVS(864) = 0

  JVS(865) = 0

  JVS(866) = 0

  JVS(867) = -B(122)

  JVS(868) = 0

  JVS(869) = 0

  JVS(870) = 0

  JVS(871) = -B(150)

  JVS(872) = -B(136)+0.21*B(198)

  JVS(873) = 0

  JVS(874) = -B(106)+0.63*B(186)

  JVS(875) = -B(124)

  JVS(876) = -B(104)+0.4*B(184)

  JVS(877) = 0

  JVS(878) = 0

  JVS(879) = 0

  JVS(880) = -B(32)+B(34)+0.12*B(183)+0.4*B(185)+0.63*B(187)+0.28*B(189)+0.25*B(191)+0.85*B(193)+0.85*B(195)+0.19*B(197)&
               &+0.21*B(199)+0.036*B(201)+B(211)+B(217)+B(223)+0.63*B(455)+0.63*B(463)

  JVS(881) = 0

  JVS(882) = 0

  JVS(883) = 0

  JVS(884) = 0

  JVS(885) = -B(33)-B(36)-B(39)-B(49)-B(51)-B(53)-B(63)-B(65)-B(67)-B(82)-B(85)-B(87)-B(89)-B(91)-B(93)-B(95)-0.964&
               &*B(97)-B(99)-B(101)-B(103)-B(105)-B(107)-B(109)-B(111)-B(113)-B(115)-B(117)-B(119)-B(121)-B(123)-B(125)&
               &-B(127)-B(129)-B(131)-B(133)-B(135)-B(137)-B(139)-0.65*B(141)-0.51*B(142)-B(145)-B(147)-B(149)-B(151)-B(161)&
               &-B(407)-B(443)-B(453)-B(459)

  JVS(886) = B(35)-B(37)+B(55)+0.7*B(60)

  JVS(887) = B(14)-0.51*B(143)

  JVS(888) = -B(50)+B(56)

  JVS(889) = -B(52)

  JVS(890) = -B(54)+0.7*B(61)

  JVS(891) = 0

  JVS(892) = B(84)

  JVS(893) = 0.1*B(116)

  JVS(894) = 0.1*B(118)

  JVS(895) = 0.25*B(98)

  JVS(896) = 0.049*B(100)

  JVS(897) = 0.381*B(96)

  JVS(898) = B(138)

  JVS(899) = B(13)

  JVS(900) = 0.65*B(6)+B(59)

  JVS(901) = B(90)

  JVS(902) = B(23)+B(128)

  JVS(903) = -B(204)

  JVS(904) = B(38)

  JVS(905) = 0.02*B(208)

  JVS(906) = 0.02*B(214)

  JVS(907) = 0.26*B(182)

  JVS(908) = 0.02*B(220)

  JVS(909) = 0.35*B(144)

  JVS(910) = 0

  JVS(911) = 0.1*B(192)

  JVS(912) = B(86)

  JVS(913) = 0.1*B(194)

  JVS(914) = 0.25*B(190)

  JVS(915) = 0.22*B(462)

  JVS(916) = 0.05*B(120)

  JVS(917) = 0.3*B(188)

  JVS(918) = 0.75*B(444)-B(446)+B(448)

  JVS(919) = 0.8*B(18)+B(130)+B(156)

  JVS(920) = 0.4*B(148)+0.08*B(200)

  JVS(921) = B(240)-B(278)+B(319)+0.5*B(359)+B(400)

  JVS(922) = B(242)-B(280)+B(321)+0.501*B(361)+B(402)

  JVS(923) = B(22)+0.32*B(196)

  JVS(924) = 0.8*B(248)-B(286)+2*B(327)+B(367)+B(408)

  JVS(925) = B(246)-B(284)+B(325)+0.506*B(365)+B(439)

  JVS(926) = B(256)-B(294)+2*B(335)+B(375)+B(416)

  JVS(927) = B(19)

  JVS(928) = 0.65*B(250)-B(288)+2*B(329)+B(369)+B(410)

  JVS(929) = 0.599*B(236)-B(274)+0.946*B(315)+0.438*B(355)+0.699*B(396)

  JVS(930) = 0.606*B(238)-B(276)+0.993*B(317)+0.489*B(357)+0.845*B(398)

  JVS(931) = 2*B(11)+B(122)+B(152)

  JVS(932) = 0.95*B(252)-B(290)+B(331)+B(371)+B(412)

  JVS(933) = 0.95*B(254)-B(292)+B(333)+B(373)+B(414)

  JVS(934) = B(244)-B(282)+B(323)+0.51*B(363)+B(404)

  JVS(935) = B(21)

  JVS(936) = B(20)+0.5*B(136)+0.5*B(162)+0.29*B(198)

  JVS(937) = -B(428)+B(430)

  JVS(938) = 0.22*B(186)

  JVS(939) = B(12)

  JVS(940) = 0.25*B(184)

  JVS(941) = -B(308)+0.5*B(349)+0.5*B(387)

  JVS(942) = B(264)-B(306)+B(347)+0.5*B(382)+B(386)+0.5*B(388)+B(424)

  JVS(943) = -B(300)-B(302)+B(341)

  JVS(944) = B(32)-B(34)+0.26*B(183)+0.25*B(185)+0.22*B(187)+0.3*B(189)+0.25*B(191)+0.1*B(193)+0.1*B(195)+0.32*B(197)&
               &+0.29*B(199)+0.08*B(201)+0.22*B(463)

  JVS(945) = 0.742*B(234)-B(272)+0.992*B(313)+0.488*B(353)+0.792*B(394)

  JVS(946) = 0.77*B(262)-B(304)+0.88*B(345)+0.38*B(380)+0.77*B(422)

  JVS(947) = -B(296)-B(298)+B(337)+0.5*B(351)+0.488*B(354)+0.438*B(356)+0.489*B(358)+0.5*B(360)+0.501*B(362)+0.51*B(364)&
               &+0.506*B(366)+B(368)+B(370)+B(372)+B(374)+B(376)+0.38*B(381)+0.5*B(383)

  JVS(948) = B(230)-B(268)+0.66*B(310)+B(311)+0.992*B(314)+0.946*B(316)+0.993*B(318)+B(320)+B(322)+B(324)+B(326)+2&
               &*B(328)+2*B(330)+B(332)+B(334)+2*B(336)+B(338)+B(342)+0.88*B(346)+B(348)+0.5*B(350)+B(390)+B(431)

  JVS(949) = B(33)-B(36)+B(39)+B(53)+B(82)+B(85)+B(87)+B(91)+0.381*B(97)+0.25*B(99)+0.049*B(101)+0.1*B(117)+0.1*B(119)&
               &+0.05*B(121)+B(123)+B(129)+B(131)+0.5*B(137)+B(139)+0.35*B(145)+0.4*B(149)

  JVS(950) = -B(35)-B(37)-2*B(40)-2*B(41)-B(55)-B(57)-B(60)-B(205)-B(269)-B(270)-B(273)-B(275)-B(277)-B(279)-B(281)&
               &-B(283)-B(285)-B(287)-B(289)-B(291)-B(293)-B(295)-B(297)-B(299)-B(301)-B(303)-B(305)-B(307)-B(309)-B(429)&
               &-B(447)

  JVS(951) = B(14)

  JVS(952) = -B(56)+B(231)+B(232)+0.742*B(235)+0.599*B(237)+0.606*B(239)+B(241)+B(243)+B(245)+B(247)+0.8*B(249)+0.65&
               &*B(251)+0.95*B(253)+0.95*B(255)+B(257)+0.77*B(263)+B(265)+0.75*B(445)

  JVS(953) = -B(58)

  JVS(954) = B(54)-B(61)+B(153)+B(157)+0.5*B(163)+B(391)+B(392)+0.792*B(395)+0.699*B(397)+0.845*B(399)+B(401)+B(403)&
               &+B(405)+B(409)+B(411)+B(413)+B(415)+B(417)+0.77*B(423)+B(425)

  JVS(955) = B(233)-B(271)+B(312)+0.5*B(352)+B(393)

  JVS(956) = B(278)

  JVS(957) = B(280)

  JVS(958) = B(286)

  JVS(959) = B(294)

  JVS(960) = B(288)

  JVS(961) = B(274)

  JVS(962) = B(276)

  JVS(963) = B(290)

  JVS(964) = B(292)

  JVS(965) = B(282)

  JVS(966) = B(428)

  JVS(967) = 0

  JVS(968) = 0

  JVS(969) = B(300)

  JVS(970) = 0

  JVS(971) = B(272)

  JVS(972) = B(304)

  JVS(973) = 0

  JVS(974) = 0

  JVS(975) = -B(142)

  JVS(976) = B(270)+B(273)+B(275)+B(277)+B(279)+B(281)+B(283)+B(287)+B(289)+B(291)+B(293)+B(295)+B(301)+B(305)+B(429)

  JVS(977) = -B(14)-B(143)

  JVS(978) = 0

  JVS(979) = 0

  JVS(980) = 0

  JVS(981) = B(271)

  JVS(982) = B(4)

  JVS(983) = -B(43)+B(45)

  JVS(984) = 0

  JVS(985) = 0

  JVS(986) = 0

  JVS(987) = 0

  JVS(988) = 0

  JVS(989) = -B(444)

  JVS(990) = -B(240)

  JVS(991) = -B(242)

  JVS(992) = 0

  JVS(993) = -B(248)

  JVS(994) = -B(246)

  JVS(995) = -B(256)

  JVS(996) = -B(250)

  JVS(997) = -B(236)

  JVS(998) = -B(238)

  JVS(999) = -B(252)

  JVS(1000) = -B(254)

  JVS(1001) = -B(244)

  JVS(1002) = -B(435)

  JVS(1003) = 0

  JVS(1004) = 0

  JVS(1005) = -B(266)

  JVS(1006) = -B(264)

  JVS(1007) = -B(260)

  JVS(1008) = -B(68)

  JVS(1009) = -B(234)

  JVS(1010) = -B(262)

  JVS(1011) = -B(258)

  JVS(1012) = -B(230)

  JVS(1013) = -B(49)

  JVS(1014) = -B(55)

  JVS(1015) = 0

  JVS(1016) = -B(44)-B(50)-B(56)-B(69)-2*B(72)-B(74)-B(231)-B(232)-B(235)-B(237)-B(239)-B(241)-B(243)-B(245)-B(247)&
                &-B(249)-B(251)-B(253)-B(255)-B(257)-B(259)-B(261)-B(263)-B(265)-B(267)-B(436)-B(445)

  JVS(1017) = B(1)+B(46)+B(76)

  JVS(1018) = B(7)-B(75)+B(77)

  JVS(1019) = -B(233)

  JVS(1020) = B(80)

  JVS(1021) = B(88)

  JVS(1022) = B(406)+B(451)

  JVS(1023) = 0.65*B(6)+B(59)+B(66)

  JVS(1024) = B(62)

  JVS(1025) = B(43)-B(45)-B(47)

  JVS(1026) = -B(202)

  JVS(1027) = -B(206)

  JVS(1028) = -B(212)

  JVS(1029) = -B(218)

  JVS(1030) = 0

  JVS(1031) = B(5)

  JVS(1032) = B(226)

  JVS(1033) = 0

  JVS(1034) = 0

  JVS(1035) = 0

  JVS(1036) = B(444)-B(449)

  JVS(1037) = 0

  JVS(1038) = 0.4*B(180)+0.7*B(200)+B(229)

  JVS(1039) = B(240)+B(400)

  JVS(1040) = B(242)+B(402)

  JVS(1041) = 0

  JVS(1042) = 0.8*B(248)+B(408)

  JVS(1043) = B(246)

  JVS(1044) = B(256)+B(416)

  JVS(1045) = 0

  JVS(1046) = 0.65*B(250)+B(410)

  JVS(1047) = 0.876*B(236)+B(396)

  JVS(1048) = 0.739*B(238)+B(398)

  JVS(1049) = 0

  JVS(1050) = 0.95*B(252)+B(412)

  JVS(1051) = 0.95*B(254)+B(414)

  JVS(1052) = B(244)+B(404)

  JVS(1053) = B(21)+B(150)

  JVS(1054) = 0.5*B(162)

  JVS(1055) = B(435)+B(437)

  JVS(1056) = 0

  JVS(1057) = 0

  JVS(1058) = 0

  JVS(1059) = 2*B(266)+0.5*B(349)+0.516*B(384)+0.5*B(387)+B(389)+2*B(426)

  JVS(1060) = B(264)+0.5*B(388)+B(424)

  JVS(1061) = -B(227)+B(260)+B(420)

  JVS(1062) = B(68)-B(70)+0.7*B(201)

  JVS(1063) = 0.941*B(234)+B(394)

  JVS(1064) = B(262)+B(422)

  JVS(1065) = -B(224)+B(258)+0.516*B(385)+B(418)

  JVS(1066) = B(230)+0.5*B(350)+B(390)

  JVS(1067) = -B(51)+B(53)+B(63)+B(67)+B(89)+B(151)+B(407)

  JVS(1068) = B(55)-B(57)+0.7*B(60)

  JVS(1069) = 0

  JVS(1070) = B(44)+B(56)+B(69)+2*B(72)+2*B(74)+B(231)+B(232)+0.941*B(235)+0.876*B(237)+0.739*B(239)+B(241)+B(243)&
                &+B(245)+B(247)+0.8*B(249)+0.65*B(251)+0.95*B(253)+0.95*B(255)+B(257)+B(259)+B(261)+B(263)+B(265)+2*B(267)&
                &+B(436)+B(445)

  JVS(1071) = -B(1)-B(46)-B(48)-B(52)-B(58)-B(71)-B(78)-B(203)-B(207)-B(213)-B(219)-B(225)-B(228)-B(450)

  JVS(1072) = B(8)+B(54)+0.7*B(61)+2*B(75)-B(79)+2*B(81)+0.5*B(163)+0.4*B(181)+B(391)+B(392)+B(395)+B(397)+B(399)+B(401)&
                &+B(403)+B(405)+B(409)+B(411)+B(413)+B(415)+B(417)+B(419)+B(421)+B(423)+B(425)+2*B(427)+B(438)

  JVS(1073) = B(233)+B(393)

  JVS(1074) = B(80)

  JVS(1075) = 0.35*B(6)

  JVS(1076) = B(47)

  JVS(1077) = -B(166)

  JVS(1078) = -B(456)

  JVS(1079) = B(64)

  JVS(1080) = B(146)

  JVS(1081) = -B(176)

  JVS(1082) = -B(178)

  JVS(1083) = -B(174)

  JVS(1084) = -B(460)

  JVS(1085) = -B(164)

  JVS(1086) = -B(172)

  JVS(1087) = -B(156)

  JVS(1088) = 0.6*B(148)-0.4*B(180)

  JVS(1089) = -B(400)

  JVS(1090) = -B(402)

  JVS(1091) = -B(408)

  JVS(1092) = -B(416)

  JVS(1093) = -B(158)

  JVS(1094) = -B(410)

  JVS(1095) = -B(396)

  JVS(1096) = -B(398)

  JVS(1097) = -B(152)

  JVS(1098) = -B(412)

  JVS(1099) = -B(414)

  JVS(1100) = -B(404)

  JVS(1101) = -B(162)

  JVS(1102) = -B(437)

  JVS(1103) = -B(170)

  JVS(1104) = -B(154)

  JVS(1105) = -B(168)

  JVS(1106) = -B(426)

  JVS(1107) = -B(424)

  JVS(1108) = -B(420)

  JVS(1109) = B(70)

  JVS(1110) = -B(394)

  JVS(1111) = -B(422)

  JVS(1112) = -B(418)

  JVS(1113) = -B(390)

  JVS(1114) = -B(53)+B(65)+B(147)+0.6*B(149)

  JVS(1115) = -B(60)

  JVS(1116) = 0

  JVS(1117) = -B(74)

  JVS(1118) = B(48)+B(71)-B(76)-B(78)

  JVS(1119) = -B(7)-B(8)-B(54)-B(61)-B(75)-B(77)-B(79)-2*B(81)-B(153)-B(155)-B(157)-B(159)-B(163)-B(165)-B(167)-B(169)&
                &-B(171)-B(173)-B(175)-B(177)-B(179)-0.4*B(181)-B(391)-B(392)-B(395)-B(397)-B(399)-B(401)-B(403)-B(405)&
                &-B(409)-B(411)-B(413)-B(415)-B(417)-B(419)-B(421)-B(423)-B(425)-B(427)-B(438)-B(457)-B(461)

  JVS(1120) = -B(393)

  JVS(1121) = B(94)

  JVS(1122) = 0.2*B(192)

  JVS(1123) = 0.16*B(194)

  JVS(1124) = 0.18*B(462)

  JVS(1125) = B(16)

  JVS(1126) = 0.245*B(236)+0.014*B(315)+0.018*B(355)+0.262*B(396)

  JVS(1127) = 0.133*B(238)+0.013*B(317)+0.015*B(357)+0.155*B(398)

  JVS(1128) = 0

  JVS(1129) = 0

  JVS(1130) = 0

  JVS(1131) = 0.18*B(186)

  JVS(1132) = 0.1*B(184)

  JVS(1133) = 0

  JVS(1134) = 0

  JVS(1135) = 0

  JVS(1136) = 0.1*B(185)+0.18*B(187)+0.2*B(193)+0.16*B(195)+0.18*B(463)

  JVS(1137) = 0.048*B(234)+0.006*B(353)+0.053*B(394)

  JVS(1138) = 0

  JVS(1139) = -B(351)+0.006*B(354)+0.018*B(356)+0.015*B(358)

  JVS(1140) = -B(311)+0.014*B(316)+0.013*B(318)

  JVS(1141) = B(95)

  JVS(1142) = -B(270)

  JVS(1143) = 0

  JVS(1144) = -B(232)+0.048*B(235)+0.245*B(237)+0.133*B(239)

  JVS(1145) = 0

  JVS(1146) = -B(392)+0.053*B(395)+0.262*B(397)+0.155*B(399)

  JVS(1147) = -B(233)-B(271)-B(312)-B(352)-B(393)
      
END SUBROUTINE racm_soa_vbs_aqchem_Jac_SP














SUBROUTINE racm_soa_vbs_aqchem_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(1147), W(90), a
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
      
END SUBROUTINE racm_soa_vbs_aqchem_KppDecomp



SUBROUTINE racm_soa_vbs_aqchem_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(1147), W(90), a
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
      
END SUBROUTINE racm_soa_vbs_aqchem_KppDecompCmplx


SUBROUTINE racm_soa_vbs_aqchem_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(1147), X(90), sum

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
      
END SUBROUTINE racm_soa_vbs_aqchem_KppSolveIndirect


SUBROUTINE racm_soa_vbs_aqchem_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(1147), X(90), sum

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
      
END SUBROUTINE racm_soa_vbs_aqchem_KppSolveCmplx













SUBROUTINE racm_soa_vbs_aqchem_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(30) = X(30)-JVS(120)*X(28)
  X(32) = X(32)-JVS(134)*X(20)
  X(33) = X(33)-JVS(140)*X(13)
  X(36) = X(36)-JVS(164)*X(16)
  X(37) = X(37)-JVS(169)*X(17)
  X(46) = X(46)-JVS(221)*X(22)-JVS(222)*X(24)-JVS(223)*X(38)-JVS(224)*X(41)-JVS(225)*X(45)
  X(50) = X(50)-JVS(258)*X(34)-JVS(259)*X(36)-JVS(260)*X(37)-JVS(261)*X(39)
  X(52) = X(52)-JVS(272)*X(21)-JVS(273)*X(28)-JVS(274)*X(48)
  X(53) = X(53)-JVS(283)*X(22)
  X(55) = X(55)-JVS(302)*X(38)
  X(57) = X(57)-JVS(319)*X(15)-JVS(320)*X(42)-JVS(321)*X(47)-JVS(322)*X(48)-JVS(323)*X(51)
  X(58) = X(58)-JVS(334)*X(45)
  X(59) = X(59)-JVS(343)*X(18)-JVS(344)*X(23)-JVS(345)*X(42)-JVS(346)*X(45)-JVS(347)*X(49)-JVS(348)*X(56)-JVS(349)*X(58)
  X(60) = X(60)-JVS(370)*X(48)-JVS(371)*X(51)
  X(61) = X(61)-JVS(380)*X(39)-JVS(381)*X(50)
  X(62) = X(62)-JVS(391)*X(30)-JVS(392)*X(32)-JVS(393)*X(41)-JVS(394)*X(48)-JVS(395)*X(52)-JVS(396)*X(54)-JVS(397)*X(57)&
            &-JVS(398)*X(60)-JVS(399)*X(61)
  X(63) = X(63)-JVS(416)*X(47)
  X(64) = X(64)-JVS(425)*X(18)
  X(65) = X(65)-JVS(433)*X(20)
  X(66) = X(66)-JVS(441)*X(22)-JVS(442)*X(24)-JVS(443)*X(27)-JVS(444)*X(32)-JVS(445)*X(38)-JVS(446)*X(40)-JVS(447)*X(41)&
            &-JVS(448)*X(42)-JVS(449)*X(44)-JVS(450)*X(47)-JVS(451)*X(48)-JVS(452)*X(49)-JVS(453)*X(51)-JVS(454)*X(52)&
            &-JVS(455)*X(53)-JVS(456)*X(54)-JVS(457)*X(55)-JVS(458)*X(56)-JVS(459)*X(57)-JVS(460)*X(58)-JVS(461)*X(60)&
            &-JVS(462)*X(61)-JVS(463)*X(63)-JVS(464)*X(64)-JVS(465)*X(65)
  X(67) = X(67)-JVS(488)*X(36)
  X(68) = X(68)-JVS(498)*X(37)
  X(69) = X(69)-JVS(508)*X(42)-JVS(509)*X(49)
  X(70) = X(70)-JVS(519)*X(34)-JVS(520)*X(50)-JVS(521)*X(54)-JVS(522)*X(58)-JVS(523)*X(63)-JVS(524)*X(64)-JVS(525)*X(65)&
            &-JVS(526)*X(67)-JVS(527)*X(68)
  X(71) = X(71)-JVS(541)*X(67)-JVS(542)*X(68)
  X(72) = X(72)-JVS(552)*X(16)-JVS(553)*X(17)-JVS(554)*X(40)-JVS(555)*X(44)-JVS(556)*X(50)-JVS(557)*X(51)-JVS(558)*X(54)&
            &-JVS(559)*X(64)-JVS(560)*X(65)-JVS(561)*X(71)
  X(73) = X(73)-JVS(575)*X(60)-JVS(576)*X(63)
  X(74) = X(74)-JVS(585)*X(20)-JVS(586)*X(22)-JVS(587)*X(23)-JVS(588)*X(42)-JVS(589)*X(45)-JVS(590)*X(49)-JVS(591)*X(55)&
            &-JVS(592)*X(56)-JVS(593)*X(58)-JVS(594)*X(64)-JVS(595)*X(65)-JVS(596)*X(69)-JVS(597)*X(70)-JVS(598)*X(71)&
            &-JVS(599)*X(73)
  X(75) = X(75)-JVS(617)*X(47)-JVS(618)*X(51)-JVS(619)*X(60)
  X(76) = X(76)-JVS(628)*X(38)-JVS(629)*X(42)-JVS(630)*X(45)-JVS(631)*X(47)-JVS(632)*X(49)-JVS(633)*X(51)-JVS(634)*X(73)&
            &-JVS(635)*X(75)
  X(77) = X(77)-JVS(645)*X(38)-JVS(646)*X(42)-JVS(647)*X(45)-JVS(648)*X(47)-JVS(649)*X(49)-JVS(650)*X(51)-JVS(651)*X(73)&
            &-JVS(652)*X(75)-JVS(653)*X(76)
  X(78) = X(78)-JVS(662)*X(54)-JVS(663)*X(71)
  X(79) = X(79)-JVS(673)*X(33)-JVS(674)*X(36)-JVS(675)*X(37)-JVS(676)*X(38)-JVS(677)*X(39)-JVS(678)*X(42)-JVS(679)*X(45)&
            &-JVS(680)*X(47)-JVS(681)*X(48)-JVS(682)*X(49)-JVS(683)*X(50)-JVS(684)*X(51)-JVS(685)*X(54)-JVS(686)*X(57)&
            &-JVS(687)*X(60)-JVS(688)*X(63)-JVS(689)*X(71)-JVS(690)*X(73)-JVS(691)*X(75)-JVS(692)*X(78)
  X(80) = X(80)-JVS(701)*X(22)-JVS(702)*X(70)-JVS(703)*X(76)-JVS(704)*X(77)-JVS(705)*X(78)-JVS(706)*X(79)
  X(81) = X(81)-JVS(716)*X(45)-JVS(717)*X(47)-JVS(718)*X(49)-JVS(719)*X(51)-JVS(720)*X(59)-JVS(721)*X(64)-JVS(722)*X(65)&
            &-JVS(723)*X(69)-JVS(724)*X(70)-JVS(725)*X(71)-JVS(726)*X(73)-JVS(727)*X(75)-JVS(728)*X(76)-JVS(729)*X(77)&
            &-JVS(730)*X(78)-JVS(731)*X(79)-JVS(732)*X(80)
  X(82) = X(82)-JVS(742)*X(32)-JVS(743)*X(40)-JVS(744)*X(44)-JVS(745)*X(48)-JVS(746)*X(51)-JVS(747)*X(52)-JVS(748)*X(54)&
            &-JVS(749)*X(55)-JVS(750)*X(56)-JVS(751)*X(57)-JVS(752)*X(58)-JVS(753)*X(59)-JVS(754)*X(60)-JVS(755)*X(61)&
            &-JVS(756)*X(62)-JVS(757)*X(63)-JVS(758)*X(64)-JVS(759)*X(65)-JVS(760)*X(67)-JVS(761)*X(68)-JVS(762)*X(69)&
            &-JVS(763)*X(70)-JVS(764)*X(71)-JVS(765)*X(72)-JVS(766)*X(73)-JVS(767)*X(74)-JVS(768)*X(75)-JVS(769)*X(76)&
            &-JVS(770)*X(77)-JVS(771)*X(78)-JVS(772)*X(79)-JVS(773)*X(80)-JVS(774)*X(81)
  X(83) = X(83)-JVS(784)*X(26)-JVS(785)*X(27)-JVS(786)*X(40)-JVS(787)*X(48)-JVS(788)*X(49)-JVS(789)*X(51)-JVS(790)*X(55)&
            &-JVS(791)*X(56)-JVS(792)*X(58)-JVS(793)*X(60)-JVS(794)*X(61)-JVS(795)*X(63)-JVS(796)*X(64)-JVS(797)*X(65)&
            &-JVS(798)*X(67)-JVS(799)*X(68)-JVS(800)*X(69)-JVS(801)*X(71)-JVS(802)*X(72)-JVS(803)*X(73)-JVS(804)*X(74)&
            &-JVS(805)*X(75)-JVS(806)*X(76)-JVS(807)*X(77)-JVS(808)*X(78)-JVS(809)*X(79)-JVS(810)*X(80)-JVS(811)*X(81)&
            &-JVS(812)*X(82)
  X(84) = X(84)-JVS(821)*X(13)-JVS(822)*X(14)-JVS(823)*X(16)-JVS(824)*X(17)-JVS(825)*X(18)-JVS(826)*X(20)-JVS(827)*X(21)&
            &-JVS(828)*X(22)-JVS(829)*X(23)-JVS(830)*X(24)-JVS(831)*X(25)-JVS(832)*X(26)-JVS(833)*X(27)-JVS(834)*X(28)&
            &-JVS(835)*X(29)-JVS(836)*X(30)-JVS(837)*X(31)-JVS(838)*X(32)-JVS(839)*X(35)-JVS(840)*X(36)-JVS(841)*X(37)&
            &-JVS(842)*X(38)-JVS(843)*X(39)-JVS(844)*X(40)-JVS(845)*X(41)-JVS(846)*X(42)-JVS(847)*X(43)-JVS(848)*X(44)&
            &-JVS(849)*X(45)-JVS(850)*X(46)-JVS(851)*X(47)-JVS(852)*X(48)-JVS(853)*X(49)-JVS(854)*X(50)-JVS(855)*X(51)&
            &-JVS(856)*X(52)-JVS(857)*X(53)-JVS(858)*X(54)-JVS(859)*X(57)-JVS(860)*X(59)-JVS(861)*X(60)-JVS(862)*X(61)&
            &-JVS(863)*X(62)-JVS(864)*X(63)-JVS(865)*X(64)-JVS(866)*X(65)-JVS(867)*X(66)-JVS(868)*X(67)-JVS(869)*X(68)&
            &-JVS(870)*X(69)-JVS(871)*X(70)-JVS(872)*X(71)-JVS(873)*X(72)-JVS(874)*X(73)-JVS(875)*X(74)-JVS(876)*X(75)&
            &-JVS(877)*X(76)-JVS(878)*X(77)-JVS(879)*X(78)-JVS(880)*X(79)-JVS(881)*X(80)-JVS(882)*X(81)-JVS(883)*X(82)&
            &-JVS(884)*X(83)
  X(85) = X(85)-JVS(892)*X(14)-JVS(893)*X(16)-JVS(894)*X(17)-JVS(895)*X(18)-JVS(896)*X(20)-JVS(897)*X(22)-JVS(898)*X(23)&
            &-JVS(899)*X(27)-JVS(900)*X(29)-JVS(901)*X(30)-JVS(902)*X(32)-JVS(903)*X(34)-JVS(904)*X(35)-JVS(905)*X(36)&
            &-JVS(906)*X(37)-JVS(907)*X(38)-JVS(908)*X(39)-JVS(909)*X(40)-JVS(910)*X(41)-JVS(911)*X(45)-JVS(912)*X(46)&
            &-JVS(913)*X(47)-JVS(914)*X(48)-JVS(915)*X(49)-JVS(916)*X(50)-JVS(917)*X(51)-JVS(918)*X(52)-JVS(919)*X(53)&
            &-JVS(920)*X(54)-JVS(921)*X(55)-JVS(922)*X(56)-JVS(923)*X(57)-JVS(924)*X(58)-JVS(925)*X(60)-JVS(926)*X(61)&
            &-JVS(927)*X(62)-JVS(928)*X(63)-JVS(929)*X(64)-JVS(930)*X(65)-JVS(931)*X(66)-JVS(932)*X(67)-JVS(933)*X(68)&
            &-JVS(934)*X(69)-JVS(935)*X(70)-JVS(936)*X(71)-JVS(937)*X(72)-JVS(938)*X(73)-JVS(939)*X(74)-JVS(940)*X(75)&
            &-JVS(941)*X(76)-JVS(942)*X(77)-JVS(943)*X(78)-JVS(944)*X(79)-JVS(945)*X(80)-JVS(946)*X(81)-JVS(947)*X(82)&
            &-JVS(948)*X(83)-JVS(949)*X(84)
  X(86) = X(86)-JVS(956)*X(55)-JVS(957)*X(56)-JVS(958)*X(58)-JVS(959)*X(61)-JVS(960)*X(63)-JVS(961)*X(64)-JVS(962)*X(65)&
            &-JVS(963)*X(67)-JVS(964)*X(68)-JVS(965)*X(69)-JVS(966)*X(72)-JVS(967)*X(73)-JVS(968)*X(75)-JVS(969)*X(78)&
            &-JVS(970)*X(79)-JVS(971)*X(80)-JVS(972)*X(81)-JVS(973)*X(82)-JVS(974)*X(83)-JVS(975)*X(84)-JVS(976)*X(85)
  X(87) = X(87)-JVS(982)*X(31)-JVS(983)*X(33)-JVS(984)*X(36)-JVS(985)*X(37)-JVS(986)*X(39)-JVS(987)*X(50)-JVS(988)*X(51)&
            &-JVS(989)*X(52)-JVS(990)*X(55)-JVS(991)*X(56)-JVS(992)*X(57)-JVS(993)*X(58)-JVS(994)*X(60)-JVS(995)*X(61)&
            &-JVS(996)*X(63)-JVS(997)*X(64)-JVS(998)*X(65)-JVS(999)*X(67)-JVS(1000)*X(68)-JVS(1001)*X(69)-JVS(1002)*X(72)&
            &-JVS(1003)*X(73)-JVS(1004)*X(75)-JVS(1005)*X(76)-JVS(1006)*X(77)-JVS(1007)*X(78)-JVS(1008)*X(79)-JVS(1009)&
            &*X(80)-JVS(1010)*X(81)-JVS(1011)*X(82)-JVS(1012)*X(83)-JVS(1013)*X(84)-JVS(1014)*X(85)-JVS(1015)*X(86)
  X(88) = X(88)-JVS(1020)*X(19)-JVS(1021)*X(24)-JVS(1022)*X(28)-JVS(1023)*X(29)-JVS(1024)*X(31)-JVS(1025)*X(33)&
            &-JVS(1026)*X(34)-JVS(1027)*X(36)-JVS(1028)*X(37)-JVS(1029)*X(39)-JVS(1030)*X(41)-JVS(1031)*X(43)-JVS(1032)&
            &*X(44)-JVS(1033)*X(48)-JVS(1034)*X(50)-JVS(1035)*X(51)-JVS(1036)*X(52)-JVS(1037)*X(53)-JVS(1038)*X(54)&
            &-JVS(1039)*X(55)-JVS(1040)*X(56)-JVS(1041)*X(57)-JVS(1042)*X(58)-JVS(1043)*X(60)-JVS(1044)*X(61)-JVS(1045)&
            &*X(62)-JVS(1046)*X(63)-JVS(1047)*X(64)-JVS(1048)*X(65)-JVS(1049)*X(66)-JVS(1050)*X(67)-JVS(1051)*X(68)&
            &-JVS(1052)*X(69)-JVS(1053)*X(70)-JVS(1054)*X(71)-JVS(1055)*X(72)-JVS(1056)*X(73)-JVS(1057)*X(74)-JVS(1058)&
            &*X(75)-JVS(1059)*X(76)-JVS(1060)*X(77)-JVS(1061)*X(78)-JVS(1062)*X(79)-JVS(1063)*X(80)-JVS(1064)*X(81)&
            &-JVS(1065)*X(82)-JVS(1066)*X(83)-JVS(1067)*X(84)-JVS(1068)*X(85)-JVS(1069)*X(86)-JVS(1070)*X(87)
  X(89) = X(89)-JVS(1074)*X(19)-JVS(1075)*X(29)-JVS(1076)*X(33)-JVS(1077)*X(38)-JVS(1078)*X(42)-JVS(1079)*X(43)&
            &-JVS(1080)*X(44)-JVS(1081)*X(45)-JVS(1082)*X(47)-JVS(1083)*X(48)-JVS(1084)*X(49)-JVS(1085)*X(50)-JVS(1086)&
            &*X(51)-JVS(1087)*X(53)-JVS(1088)*X(54)-JVS(1089)*X(55)-JVS(1090)*X(56)-JVS(1091)*X(58)-JVS(1092)*X(61)&
            &-JVS(1093)*X(62)-JVS(1094)*X(63)-JVS(1095)*X(64)-JVS(1096)*X(65)-JVS(1097)*X(66)-JVS(1098)*X(67)-JVS(1099)&
            &*X(68)-JVS(1100)*X(69)-JVS(1101)*X(71)-JVS(1102)*X(72)-JVS(1103)*X(73)-JVS(1104)*X(74)-JVS(1105)*X(75)&
            &-JVS(1106)*X(76)-JVS(1107)*X(77)-JVS(1108)*X(78)-JVS(1109)*X(79)-JVS(1110)*X(80)-JVS(1111)*X(81)-JVS(1112)&
            &*X(82)-JVS(1113)*X(83)-JVS(1114)*X(84)-JVS(1115)*X(85)-JVS(1116)*X(86)-JVS(1117)*X(87)-JVS(1118)*X(88)
  X(90) = X(90)-JVS(1121)*X(25)-JVS(1122)*X(45)-JVS(1123)*X(47)-JVS(1124)*X(49)-JVS(1125)*X(59)-JVS(1126)*X(64)&
            &-JVS(1127)*X(65)-JVS(1128)*X(69)-JVS(1129)*X(70)-JVS(1130)*X(71)-JVS(1131)*X(73)-JVS(1132)*X(75)-JVS(1133)&
            &*X(76)-JVS(1134)*X(77)-JVS(1135)*X(78)-JVS(1136)*X(79)-JVS(1137)*X(80)-JVS(1138)*X(81)-JVS(1139)*X(82)&
            &-JVS(1140)*X(83)-JVS(1141)*X(84)-JVS(1142)*X(85)-JVS(1143)*X(86)-JVS(1144)*X(87)-JVS(1145)*X(88)-JVS(1146)&
            &*X(89)
  X(90) = X(90)/JVS(1147)
  X(89) = (X(89)-JVS(1120)*X(90))/(JVS(1119))
  X(88) = (X(88)-JVS(1072)*X(89)-JVS(1073)*X(90))/(JVS(1071))
  X(87) = (X(87)-JVS(1017)*X(88)-JVS(1018)*X(89)-JVS(1019)*X(90))/(JVS(1016))
  X(86) = (X(86)-JVS(978)*X(87)-JVS(979)*X(88)-JVS(980)*X(89)-JVS(981)*X(90))/(JVS(977))
  X(85) = (X(85)-JVS(951)*X(86)-JVS(952)*X(87)-JVS(953)*X(88)-JVS(954)*X(89)-JVS(955)*X(90))/(JVS(950))
  X(84) = (X(84)-JVS(886)*X(85)-JVS(887)*X(86)-JVS(888)*X(87)-JVS(889)*X(88)-JVS(890)*X(89)-JVS(891)*X(90))/(JVS(885))
  X(83) = (X(83)-JVS(814)*X(84)-JVS(815)*X(85)-JVS(816)*X(86)-JVS(817)*X(87)-JVS(818)*X(88)-JVS(819)*X(89)-JVS(820)&
            &*X(90))/(JVS(813))
  X(82) = (X(82)-JVS(776)*X(83)-JVS(777)*X(84)-JVS(778)*X(85)-JVS(779)*X(86)-JVS(780)*X(87)-JVS(781)*X(88)-JVS(782)&
            &*X(89)-JVS(783)*X(90))/(JVS(775))
  X(81) = (X(81)-JVS(734)*X(82)-JVS(735)*X(83)-JVS(736)*X(84)-JVS(737)*X(85)-JVS(738)*X(86)-JVS(739)*X(87)-JVS(740)&
            &*X(88)-JVS(741)*X(89))/(JVS(733))
  X(80) = (X(80)-JVS(708)*X(82)-JVS(709)*X(83)-JVS(710)*X(84)-JVS(711)*X(85)-JVS(712)*X(86)-JVS(713)*X(87)-JVS(714)&
            &*X(88)-JVS(715)*X(89))/(JVS(707))
  X(79) = (X(79)-JVS(694)*X(82)-JVS(695)*X(83)-JVS(696)*X(84)-JVS(697)*X(85)-JVS(698)*X(87)-JVS(699)*X(88)-JVS(700)&
            &*X(89))/(JVS(693))
  X(78) = (X(78)-JVS(665)*X(79)-JVS(666)*X(82)-JVS(667)*X(83)-JVS(668)*X(84)-JVS(669)*X(85)-JVS(670)*X(87)-JVS(671)&
            &*X(88)-JVS(672)*X(89))/(JVS(664))
  X(77) = (X(77)-JVS(655)*X(79)-JVS(656)*X(82)-JVS(657)*X(83)-JVS(658)*X(84)-JVS(659)*X(85)-JVS(660)*X(87)-JVS(661)&
            &*X(89))/(JVS(654))
  X(76) = (X(76)-JVS(637)*X(77)-JVS(638)*X(79)-JVS(639)*X(82)-JVS(640)*X(83)-JVS(641)*X(84)-JVS(642)*X(85)-JVS(643)&
            &*X(87)-JVS(644)*X(89))/(JVS(636))
  X(75) = (X(75)-JVS(621)*X(79)-JVS(622)*X(82)-JVS(623)*X(83)-JVS(624)*X(84)-JVS(625)*X(85)-JVS(626)*X(87)-JVS(627)&
            &*X(89))/(JVS(620))
  X(74) = (X(74)-JVS(601)*X(75)-JVS(602)*X(76)-JVS(603)*X(77)-JVS(604)*X(78)-JVS(605)*X(79)-JVS(606)*X(80)-JVS(607)&
            &*X(81)-JVS(608)*X(82)-JVS(609)*X(83)-JVS(610)*X(84)-JVS(611)*X(85)-JVS(612)*X(86)-JVS(613)*X(87)-JVS(614)*X(88)&
            &-JVS(615)*X(89)-JVS(616)*X(90))/(JVS(600))
  X(73) = (X(73)-JVS(578)*X(79)-JVS(579)*X(82)-JVS(580)*X(83)-JVS(581)*X(84)-JVS(582)*X(85)-JVS(583)*X(87)-JVS(584)&
            &*X(89))/(JVS(577))
  X(72) = (X(72)-JVS(563)*X(78)-JVS(564)*X(79)-JVS(565)*X(80)-JVS(566)*X(81)-JVS(567)*X(82)-JVS(568)*X(83)-JVS(569)&
            &*X(84)-JVS(570)*X(85)-JVS(571)*X(86)-JVS(572)*X(87)-JVS(573)*X(88)-JVS(574)*X(89))/(JVS(562))
  X(71) = (X(71)-JVS(544)*X(79)-JVS(545)*X(82)-JVS(546)*X(83)-JVS(547)*X(84)-JVS(548)*X(85)-JVS(549)*X(87)-JVS(550)&
            &*X(88)-JVS(551)*X(89))/(JVS(543))
  X(70) = (X(70)-JVS(529)*X(76)-JVS(530)*X(77)-JVS(531)*X(78)-JVS(532)*X(79)-JVS(533)*X(80)-JVS(534)*X(82)-JVS(535)&
            &*X(83)-JVS(536)*X(84)-JVS(537)*X(85)-JVS(538)*X(87)-JVS(539)*X(88)-JVS(540)*X(89))/(JVS(528))
  X(69) = (X(69)-JVS(511)*X(73)-JVS(512)*X(79)-JVS(513)*X(82)-JVS(514)*X(83)-JVS(515)*X(84)-JVS(516)*X(85)-JVS(517)&
            &*X(87)-JVS(518)*X(89))/(JVS(510))
  X(68) = (X(68)-JVS(500)*X(79)-JVS(501)*X(82)-JVS(502)*X(83)-JVS(503)*X(84)-JVS(504)*X(85)-JVS(505)*X(87)-JVS(506)&
            &*X(88)-JVS(507)*X(89))/(JVS(499))
  X(67) = (X(67)-JVS(490)*X(79)-JVS(491)*X(82)-JVS(492)*X(83)-JVS(493)*X(84)-JVS(494)*X(85)-JVS(495)*X(87)-JVS(496)&
            &*X(88)-JVS(497)*X(89))/(JVS(489))
  X(66) = (X(66)-JVS(467)*X(67)-JVS(468)*X(68)-JVS(469)*X(69)-JVS(470)*X(71)-JVS(471)*X(72)-JVS(472)*X(73)-JVS(473)&
            &*X(75)-JVS(474)*X(76)-JVS(475)*X(77)-JVS(476)*X(78)-JVS(477)*X(79)-JVS(478)*X(80)-JVS(479)*X(81)-JVS(480)*X(82)&
            &-JVS(481)*X(83)-JVS(482)*X(84)-JVS(483)*X(85)-JVS(484)*X(87)-JVS(485)*X(88)-JVS(486)*X(89)-JVS(487)*X(90))&
            &/(JVS(466))
  X(65) = (X(65)-JVS(435)*X(82)-JVS(436)*X(83)-JVS(437)*X(84)-JVS(438)*X(85)-JVS(439)*X(87)-JVS(440)*X(89))/(JVS(434))
  X(64) = (X(64)-JVS(427)*X(82)-JVS(428)*X(83)-JVS(429)*X(84)-JVS(430)*X(85)-JVS(431)*X(87)-JVS(432)*X(89))/(JVS(426))
  X(63) = (X(63)-JVS(418)*X(79)-JVS(419)*X(82)-JVS(420)*X(83)-JVS(421)*X(84)-JVS(422)*X(85)-JVS(423)*X(87)-JVS(424)&
            &*X(89))/(JVS(417))
  X(62) = (X(62)-JVS(401)*X(63)-JVS(402)*X(67)-JVS(403)*X(68)-JVS(404)*X(71)-JVS(405)*X(78)-JVS(406)*X(79)-JVS(407)&
            &*X(80)-JVS(408)*X(81)-JVS(409)*X(82)-JVS(410)*X(83)-JVS(411)*X(84)-JVS(412)*X(85)-JVS(413)*X(87)-JVS(414)*X(88)&
            &-JVS(415)*X(89))/(JVS(400))
  X(61) = (X(61)-JVS(383)*X(79)-JVS(384)*X(82)-JVS(385)*X(83)-JVS(386)*X(84)-JVS(387)*X(85)-JVS(388)*X(87)-JVS(389)&
            &*X(88)-JVS(390)*X(89))/(JVS(382))
  X(60) = (X(60)-JVS(373)*X(79)-JVS(374)*X(82)-JVS(375)*X(83)-JVS(376)*X(84)-JVS(377)*X(85)-JVS(378)*X(87)-JVS(379)&
            &*X(89))/(JVS(372))
  X(59) = (X(59)-JVS(351)*X(64)-JVS(352)*X(65)-JVS(353)*X(69)-JVS(354)*X(70)-JVS(355)*X(71)-JVS(356)*X(73)-JVS(357)&
            &*X(75)-JVS(358)*X(76)-JVS(359)*X(77)-JVS(360)*X(79)-JVS(361)*X(80)-JVS(362)*X(81)-JVS(363)*X(82)-JVS(364)*X(83)&
            &-JVS(365)*X(84)-JVS(366)*X(85)-JVS(367)*X(86)-JVS(368)*X(87)-JVS(369)*X(89))/(JVS(350))
  X(58) = (X(58)-JVS(336)*X(79)-JVS(337)*X(82)-JVS(338)*X(83)-JVS(339)*X(84)-JVS(340)*X(85)-JVS(341)*X(87)-JVS(342)&
            &*X(89))/(JVS(335))
  X(57) = (X(57)-JVS(325)*X(60)-JVS(326)*X(63)-JVS(327)*X(79)-JVS(328)*X(82)-JVS(329)*X(83)-JVS(330)*X(84)-JVS(331)&
            &*X(85)-JVS(332)*X(87)-JVS(333)*X(89))/(JVS(324))
  X(56) = (X(56)-JVS(312)*X(75)-JVS(313)*X(82)-JVS(314)*X(83)-JVS(315)*X(84)-JVS(316)*X(85)-JVS(317)*X(87)-JVS(318)&
            &*X(89))/(JVS(311))
  X(55) = (X(55)-JVS(304)*X(79)-JVS(305)*X(82)-JVS(306)*X(83)-JVS(307)*X(84)-JVS(308)*X(85)-JVS(309)*X(87)-JVS(310)&
            &*X(89))/(JVS(303))
  X(54) = (X(54)-JVS(297)*X(78)-JVS(298)*X(79)-JVS(299)*X(84)-JVS(300)*X(88)-JVS(301)*X(89))/(JVS(296))
  X(53) = (X(53)-JVS(285)*X(61)-JVS(286)*X(67)-JVS(287)*X(68)-JVS(288)*X(71)-JVS(289)*X(79)-JVS(290)*X(80)-JVS(291)&
            &*X(82)-JVS(292)*X(83)-JVS(293)*X(84)-JVS(294)*X(87)-JVS(295)*X(89))/(JVS(284))
  X(52) = (X(52)-JVS(276)*X(57)-JVS(277)*X(79)-JVS(278)*X(84)-JVS(279)*X(85)-JVS(280)*X(87)-JVS(281)*X(88)-JVS(282)&
            &*X(89))/(JVS(275))
  X(51) = (X(51)-JVS(269)*X(79)-JVS(270)*X(84)-JVS(271)*X(89))/(JVS(268))
  X(50) = (X(50)-JVS(263)*X(79)-JVS(264)*X(84)-JVS(265)*X(85)-JVS(266)*X(88)-JVS(267)*X(89))/(JVS(262))
  X(49) = (X(49)-JVS(255)*X(79)-JVS(256)*X(84)-JVS(257)*X(89))/(JVS(254))
  X(48) = (X(48)-JVS(251)*X(79)-JVS(252)*X(84)-JVS(253)*X(89))/(JVS(250))
  X(47) = (X(47)-JVS(247)*X(79)-JVS(248)*X(84)-JVS(249)*X(89))/(JVS(246))
  X(46) = (X(46)-JVS(227)*X(47)-JVS(228)*X(48)-JVS(229)*X(49)-JVS(230)*X(51)-JVS(231)*X(52)-JVS(232)*X(53)-JVS(233)&
            &*X(54)-JVS(234)*X(57)-JVS(235)*X(60)-JVS(236)*X(62)-JVS(237)*X(66)-JVS(238)*X(71)-JVS(239)*X(73)-JVS(240)*X(74)&
            &-JVS(241)*X(75)-JVS(242)*X(79)-JVS(243)*X(84)-JVS(244)*X(87)-JVS(245)*X(89))/(JVS(226))
  X(45) = (X(45)-JVS(218)*X(79)-JVS(219)*X(84)-JVS(220)*X(89))/(JVS(217))
  X(44) = (X(44)-JVS(211)*X(54)-JVS(212)*X(79)-JVS(213)*X(82)-JVS(214)*X(84)-JVS(215)*X(88)-JVS(216)*X(89))/(JVS(210))
  X(43) = (X(43)-JVS(200)*X(50)-JVS(201)*X(53)-JVS(202)*X(62)-JVS(203)*X(66)-JVS(204)*X(71)-JVS(205)*X(74)-JVS(206)&
            &*X(84)-JVS(207)*X(85)-JVS(208)*X(88)-JVS(209)*X(89))/(JVS(199))
  X(42) = (X(42)-JVS(196)*X(79)-JVS(197)*X(84)-JVS(198)*X(89))/(JVS(195))
  X(41) = (X(41)-JVS(190)*X(48)-JVS(191)*X(60)-JVS(192)*X(84)-JVS(193)*X(87)-JVS(194)*X(89))/(JVS(189))
  X(40) = (X(40)-JVS(184)*X(71)-JVS(185)*X(79)-JVS(186)*X(82)-JVS(187)*X(84)-JVS(188)*X(85))/(JVS(183))
  X(39) = (X(39)-JVS(179)*X(50)-JVS(180)*X(79)-JVS(181)*X(84)-JVS(182)*X(88))/(JVS(178))
  X(38) = (X(38)-JVS(175)*X(79)-JVS(176)*X(84)-JVS(177)*X(89))/(JVS(174))
  X(37) = (X(37)-JVS(171)*X(79)-JVS(172)*X(84)-JVS(173)*X(88))/(JVS(170))
  X(36) = (X(36)-JVS(166)*X(79)-JVS(167)*X(84)-JVS(168)*X(88))/(JVS(165))
  X(35) = (X(35)-JVS(154)*X(45)-JVS(155)*X(47)-JVS(156)*X(48)-JVS(157)*X(49)-JVS(158)*X(51)-JVS(159)*X(73)-JVS(160)&
            &*X(75)-JVS(161)*X(79)-JVS(162)*X(84)-JVS(163)*X(85))/(JVS(153))
  X(34) = (X(34)-JVS(148)*X(50)-JVS(149)*X(84)-JVS(150)*X(85)-JVS(151)*X(88)-JVS(152)*X(89))/(JVS(147))
  X(33) = (X(33)-JVS(142)*X(51)-JVS(143)*X(79)-JVS(144)*X(87)-JVS(145)*X(88)-JVS(146)*X(89))/(JVS(141))
  X(32) = (X(32)-JVS(136)*X(54)-JVS(137)*X(81)-JVS(138)*X(83)-JVS(139)*X(84))/(JVS(135))
  X(31) = (X(31)-JVS(128)*X(36)-JVS(129)*X(37)-JVS(130)*X(39)-JVS(131)*X(84)-JVS(132)*X(87)-JVS(133)*X(88))/(JVS(127))
  X(30) = (X(30)-JVS(122)*X(41)-JVS(123)*X(52)-JVS(124)*X(84)-JVS(125)*X(87)-JVS(126)*X(88))/(JVS(121))
  X(29) = (X(29)-JVS(117)*X(84)-JVS(118)*X(85)-JVS(119)*X(88))/(JVS(116))
  X(28) = (X(28)-JVS(113)*X(52)-JVS(114)*X(84)-JVS(115)*X(88))/(JVS(112))
  X(27) = (X(27)-JVS(109)*X(83)-JVS(110)*X(84)-JVS(111)*X(85))/(JVS(108))
  X(26) = (X(26)-JVS(103)*X(49)-JVS(104)*X(73)-JVS(105)*X(75)-JVS(106)*X(79)-JVS(107)*X(84))/(JVS(102))
  X(25) = (X(25)-JVS(97)*X(49)-JVS(98)*X(73)-JVS(99)*X(75)-JVS(100)*X(79)-JVS(101)*X(84))/(JVS(96))
  X(24) = (X(24)-JVS(94)*X(41)-JVS(95)*X(84))/(JVS(93))
  X(23) = (X(23)-JVS(91)*X(71)-JVS(92)*X(84))/(JVS(90))
  X(22) = (X(22)-JVS(89)*X(84))/(JVS(88))
  X(21) = (X(21)-JVS(85)*X(52)-JVS(86)*X(84)-JVS(87)*X(85))/(JVS(84))
  X(20) = (X(20)-JVS(83)*X(84))/(JVS(82))
  X(19) = (X(19)-JVS(80)*X(88)-JVS(81)*X(89))/(JVS(79))
  X(18) = (X(18)-JVS(78)*X(84))/(JVS(77))
  X(17) = (X(17)-JVS(76)*X(84))/(JVS(75))
  X(16) = (X(16)-JVS(74)*X(84))/(JVS(73))
  X(15) = (X(15)-JVS(70)*X(60)-JVS(71)*X(84)-JVS(72)*X(85))/(JVS(69))
  X(14) = (X(14)-JVS(68)*X(84))/(JVS(67))
  X(13) = (X(13)-JVS(66)*X(79))/(JVS(65))
  X(12) = (X(12)-JVS(64)*X(84))/(JVS(63))
  X(11) = (X(11)-JVS(61)*X(12)-JVS(62)*X(84))/(JVS(60))
  X(10) = (X(10)-JVS(58)*X(11)-JVS(59)*X(84))/(JVS(57))
  X(9) = (X(9)-JVS(55)*X(10)-JVS(56)*X(84))/(JVS(54))
  X(8) = (X(8)-JVS(53)*X(84))/(JVS(52))
  X(7) = (X(7)-JVS(50)*X(8)-JVS(51)*X(84))/(JVS(49))
  X(6) = (X(6)-JVS(47)*X(7)-JVS(48)*X(84))/(JVS(46))
  X(5) = (X(5)-JVS(44)*X(6)-JVS(45)*X(84))/(JVS(43))
  X(4) = (X(4)-JVS(21)*X(42)-JVS(22)*X(47)-JVS(23)*X(49)-JVS(24)*X(55)-JVS(25)*X(56)-JVS(26)*X(60)-JVS(27)*X(64)-JVS(28)&
           &*X(65)-JVS(29)*X(69)-JVS(30)*X(71)-JVS(31)*X(73)-JVS(32)*X(75)-JVS(33)*X(76)-JVS(34)*X(77)-JVS(35)*X(78)-JVS(36)&
           &*X(79)-JVS(37)*X(80)-JVS(38)*X(81)-JVS(39)*X(82)-JVS(40)*X(83)-JVS(41)*X(85)-JVS(42)*X(90))/(JVS(20))
  X(3) = (X(3)-JVS(8)*X(22)-JVS(9)*X(38)-JVS(10)*X(42)-JVS(11)*X(47)-JVS(12)*X(48)-JVS(13)*X(51)-JVS(14)*X(54)-JVS(15)&
           &*X(57)-JVS(16)*X(71)-JVS(17)*X(75)-JVS(18)*X(79)-JVS(19)*X(84))/(JVS(7))
  X(2) = (X(2)-JVS(5)*X(46)-JVS(6)*X(84))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(14)-JVS(3)*X(84))/(JVS(1))
      
END SUBROUTINE racm_soa_vbs_aqchem_KppSolve
























      SUBROUTINE racm_soa_vbs_aqchem_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE racm_soa_vbs_aqchem_WCOPY



      SUBROUTINE racm_soa_vbs_aqchem_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE racm_soa_vbs_aqchem_WAXPY




      SUBROUTINE racm_soa_vbs_aqchem_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE racm_soa_vbs_aqchem_WSCAL


      REAL(kind=dp) FUNCTION racm_soa_vbs_aqchem_WLAMCH( C )








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
          CALL racm_soa_vbs_aqchem_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      racm_soa_vbs_aqchem_WLAMCH = Eps

      END FUNCTION racm_soa_vbs_aqchem_WLAMCH
     
      SUBROUTINE racm_soa_vbs_aqchem_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE racm_soa_vbs_aqchem_WLAMCH_ADD




      SUBROUTINE racm_soa_vbs_aqchem_SET2ZERO(N,Y)




      
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

      END SUBROUTINE racm_soa_vbs_aqchem_SET2ZERO



      REAL(kind=dp) FUNCTION racm_soa_vbs_aqchem_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      racm_soa_vbs_aqchem_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        racm_soa_vbs_aqchem_WDOT = racm_soa_vbs_aqchem_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         racm_soa_vbs_aqchem_WDOT = racm_soa_vbs_aqchem_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          racm_soa_vbs_aqchem_WDOT = racm_soa_vbs_aqchem_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        racm_soa_vbs_aqchem_WDOT = racm_soa_vbs_aqchem_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION racm_soa_vbs_aqchem_WDOT                                          




   SUBROUTINE decomp_racm_soa_vbs_aqchem( JVS, IER )
   
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
   W( 14 ) = JVS( 2 )
   W( 84 ) = JVS( 3 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 14 )
  JVS( 3) = W( 84 )
  IF ( ABS(  JVS( 4 )) < TINY(a) ) THEN
         IER = 2                                       
         RETURN
  END IF
   W( 2 ) = JVS( 4 )
   W( 46 ) = JVS( 5 )
   W( 84 ) = JVS( 6 )
  JVS( 4) = W( 2 )
  JVS( 5) = W( 46 )
  JVS( 6) = W( 84 )
  IF ( ABS(  JVS( 7 )) < TINY(a) ) THEN
         IER = 3                                       
         RETURN
  END IF
   W( 3 ) = JVS( 7 )
   W( 22 ) = JVS( 8 )
   W( 38 ) = JVS( 9 )
   W( 42 ) = JVS( 10 )
   W( 47 ) = JVS( 11 )
   W( 48 ) = JVS( 12 )
   W( 51 ) = JVS( 13 )
   W( 54 ) = JVS( 14 )
   W( 57 ) = JVS( 15 )
   W( 71 ) = JVS( 16 )
   W( 75 ) = JVS( 17 )
   W( 79 ) = JVS( 18 )
   W( 84 ) = JVS( 19 )
  JVS( 7) = W( 3 )
  JVS( 8) = W( 22 )
  JVS( 9) = W( 38 )
  JVS( 10) = W( 42 )
  JVS( 11) = W( 47 )
  JVS( 12) = W( 48 )
  JVS( 13) = W( 51 )
  JVS( 14) = W( 54 )
  JVS( 15) = W( 57 )
  JVS( 16) = W( 71 )
  JVS( 17) = W( 75 )
  JVS( 18) = W( 79 )
  JVS( 19) = W( 84 )
  IF ( ABS(  JVS( 20 )) < TINY(a) ) THEN
         IER = 4                                       
         RETURN
  END IF
   W( 4 ) = JVS( 20 )
   W( 42 ) = JVS( 21 )
   W( 47 ) = JVS( 22 )
   W( 49 ) = JVS( 23 )
   W( 55 ) = JVS( 24 )
   W( 56 ) = JVS( 25 )
   W( 60 ) = JVS( 26 )
   W( 64 ) = JVS( 27 )
   W( 65 ) = JVS( 28 )
   W( 69 ) = JVS( 29 )
   W( 71 ) = JVS( 30 )
   W( 73 ) = JVS( 31 )
   W( 75 ) = JVS( 32 )
   W( 76 ) = JVS( 33 )
   W( 77 ) = JVS( 34 )
   W( 78 ) = JVS( 35 )
   W( 79 ) = JVS( 36 )
   W( 80 ) = JVS( 37 )
   W( 81 ) = JVS( 38 )
   W( 82 ) = JVS( 39 )
   W( 83 ) = JVS( 40 )
   W( 85 ) = JVS( 41 )
   W( 90 ) = JVS( 42 )
  JVS( 20) = W( 4 )
  JVS( 21) = W( 42 )
  JVS( 22) = W( 47 )
  JVS( 23) = W( 49 )
  JVS( 24) = W( 55 )
  JVS( 25) = W( 56 )
  JVS( 26) = W( 60 )
  JVS( 27) = W( 64 )
  JVS( 28) = W( 65 )
  JVS( 29) = W( 69 )
  JVS( 30) = W( 71 )
  JVS( 31) = W( 73 )
  JVS( 32) = W( 75 )
  JVS( 33) = W( 76 )
  JVS( 34) = W( 77 )
  JVS( 35) = W( 78 )
  JVS( 36) = W( 79 )
  JVS( 37) = W( 80 )
  JVS( 38) = W( 81 )
  JVS( 39) = W( 82 )
  JVS( 40) = W( 83 )
  JVS( 41) = W( 85 )
  JVS( 42) = W( 90 )
  IF ( ABS(  JVS( 43 )) < TINY(a) ) THEN
         IER = 5                                       
         RETURN
  END IF
   W( 5 ) = JVS( 43 )
   W( 6 ) = JVS( 44 )
   W( 84 ) = JVS( 45 )
  JVS( 43) = W( 5 )
  JVS( 44) = W( 6 )
  JVS( 45) = W( 84 )
  IF ( ABS(  JVS( 46 )) < TINY(a) ) THEN
         IER = 6                                       
         RETURN
  END IF
   W( 6 ) = JVS( 46 )
   W( 7 ) = JVS( 47 )
   W( 84 ) = JVS( 48 )
  JVS( 46) = W( 6 )
  JVS( 47) = W( 7 )
  JVS( 48) = W( 84 )
  IF ( ABS(  JVS( 49 )) < TINY(a) ) THEN
         IER = 7                                       
         RETURN
  END IF
   W( 7 ) = JVS( 49 )
   W( 8 ) = JVS( 50 )
   W( 84 ) = JVS( 51 )
  JVS( 49) = W( 7 )
  JVS( 50) = W( 8 )
  JVS( 51) = W( 84 )
  IF ( ABS(  JVS( 52 )) < TINY(a) ) THEN
         IER = 8                                       
         RETURN
  END IF
   W( 8 ) = JVS( 52 )
   W( 84 ) = JVS( 53 )
  JVS( 52) = W( 8 )
  JVS( 53) = W( 84 )
  IF ( ABS(  JVS( 54 )) < TINY(a) ) THEN
         IER = 9                                       
         RETURN
  END IF
   W( 9 ) = JVS( 54 )
   W( 10 ) = JVS( 55 )
   W( 84 ) = JVS( 56 )
  JVS( 54) = W( 9 )
  JVS( 55) = W( 10 )
  JVS( 56) = W( 84 )
  IF ( ABS(  JVS( 57 )) < TINY(a) ) THEN
         IER = 10                                      
         RETURN
  END IF
   W( 10 ) = JVS( 57 )
   W( 11 ) = JVS( 58 )
   W( 84 ) = JVS( 59 )
  JVS( 57) = W( 10 )
  JVS( 58) = W( 11 )
  JVS( 59) = W( 84 )
  IF ( ABS(  JVS( 60 )) < TINY(a) ) THEN
         IER = 11                                      
         RETURN
  END IF
   W( 11 ) = JVS( 60 )
   W( 12 ) = JVS( 61 )
   W( 84 ) = JVS( 62 )
  JVS( 60) = W( 11 )
  JVS( 61) = W( 12 )
  JVS( 62) = W( 84 )
  IF ( ABS(  JVS( 63 )) < TINY(a) ) THEN
         IER = 12                                      
         RETURN
  END IF
   W( 12 ) = JVS( 63 )
   W( 84 ) = JVS( 64 )
  JVS( 63) = W( 12 )
  JVS( 64) = W( 84 )
  IF ( ABS(  JVS( 65 )) < TINY(a) ) THEN
         IER = 13                                      
         RETURN
  END IF
   W( 13 ) = JVS( 65 )
   W( 79 ) = JVS( 66 )
  JVS( 65) = W( 13 )
  JVS( 66) = W( 79 )
  IF ( ABS(  JVS( 67 )) < TINY(a) ) THEN
         IER = 14                                      
         RETURN
  END IF
   W( 14 ) = JVS( 67 )
   W( 84 ) = JVS( 68 )
  JVS( 67) = W( 14 )
  JVS( 68) = W( 84 )
  IF ( ABS(  JVS( 69 )) < TINY(a) ) THEN
         IER = 15                                      
         RETURN
  END IF
   W( 15 ) = JVS( 69 )
   W( 60 ) = JVS( 70 )
   W( 84 ) = JVS( 71 )
   W( 85 ) = JVS( 72 )
  JVS( 69) = W( 15 )
  JVS( 70) = W( 60 )
  JVS( 71) = W( 84 )
  JVS( 72) = W( 85 )
  IF ( ABS(  JVS( 73 )) < TINY(a) ) THEN
         IER = 16                                      
         RETURN
  END IF
   W( 16 ) = JVS( 73 )
   W( 84 ) = JVS( 74 )
  JVS( 73) = W( 16 )
  JVS( 74) = W( 84 )
  IF ( ABS(  JVS( 75 )) < TINY(a) ) THEN
         IER = 17                                      
         RETURN
  END IF
   W( 17 ) = JVS( 75 )
   W( 84 ) = JVS( 76 )
  JVS( 75) = W( 17 )
  JVS( 76) = W( 84 )
  IF ( ABS(  JVS( 77 )) < TINY(a) ) THEN
         IER = 18                                      
         RETURN
  END IF
   W( 18 ) = JVS( 77 )
   W( 84 ) = JVS( 78 )
  JVS( 77) = W( 18 )
  JVS( 78) = W( 84 )
  IF ( ABS(  JVS( 79 )) < TINY(a) ) THEN
         IER = 19                                      
         RETURN
  END IF
   W( 19 ) = JVS( 79 )
   W( 88 ) = JVS( 80 )
   W( 89 ) = JVS( 81 )
  JVS( 79) = W( 19 )
  JVS( 80) = W( 88 )
  JVS( 81) = W( 89 )
  IF ( ABS(  JVS( 82 )) < TINY(a) ) THEN
         IER = 20                                      
         RETURN
  END IF
   W( 20 ) = JVS( 82 )
   W( 84 ) = JVS( 83 )
  JVS( 82) = W( 20 )
  JVS( 83) = W( 84 )
  IF ( ABS(  JVS( 84 )) < TINY(a) ) THEN
         IER = 21                                      
         RETURN
  END IF
   W( 21 ) = JVS( 84 )
   W( 52 ) = JVS( 85 )
   W( 84 ) = JVS( 86 )
   W( 85 ) = JVS( 87 )
  JVS( 84) = W( 21 )
  JVS( 85) = W( 52 )
  JVS( 86) = W( 84 )
  JVS( 87) = W( 85 )
  IF ( ABS(  JVS( 88 )) < TINY(a) ) THEN
         IER = 22                                      
         RETURN
  END IF
   W( 22 ) = JVS( 88 )
   W( 84 ) = JVS( 89 )
  JVS( 88) = W( 22 )
  JVS( 89) = W( 84 )
  IF ( ABS(  JVS( 90 )) < TINY(a) ) THEN
         IER = 23                                      
         RETURN
  END IF
   W( 23 ) = JVS( 90 )
   W( 71 ) = JVS( 91 )
   W( 84 ) = JVS( 92 )
  JVS( 90) = W( 23 )
  JVS( 91) = W( 71 )
  JVS( 92) = W( 84 )
  IF ( ABS(  JVS( 93 )) < TINY(a) ) THEN
         IER = 24                                      
         RETURN
  END IF
   W( 24 ) = JVS( 93 )
   W( 41 ) = JVS( 94 )
   W( 84 ) = JVS( 95 )
  JVS( 93) = W( 24 )
  JVS( 94) = W( 41 )
  JVS( 95) = W( 84 )
  IF ( ABS(  JVS( 96 )) < TINY(a) ) THEN
         IER = 25                                      
         RETURN
  END IF
   W( 25 ) = JVS( 96 )
   W( 49 ) = JVS( 97 )
   W( 73 ) = JVS( 98 )
   W( 75 ) = JVS( 99 )
   W( 79 ) = JVS( 100 )
   W( 84 ) = JVS( 101 )
  JVS( 96) = W( 25 )
  JVS( 97) = W( 49 )
  JVS( 98) = W( 73 )
  JVS( 99) = W( 75 )
  JVS( 100) = W( 79 )
  JVS( 101) = W( 84 )
  IF ( ABS(  JVS( 102 )) < TINY(a) ) THEN
         IER = 26                                      
         RETURN
  END IF
   W( 26 ) = JVS( 102 )
   W( 49 ) = JVS( 103 )
   W( 73 ) = JVS( 104 )
   W( 75 ) = JVS( 105 )
   W( 79 ) = JVS( 106 )
   W( 84 ) = JVS( 107 )
  JVS( 102) = W( 26 )
  JVS( 103) = W( 49 )
  JVS( 104) = W( 73 )
  JVS( 105) = W( 75 )
  JVS( 106) = W( 79 )
  JVS( 107) = W( 84 )
  IF ( ABS(  JVS( 108 )) < TINY(a) ) THEN
         IER = 27                                      
         RETURN
  END IF
   W( 27 ) = JVS( 108 )
   W( 83 ) = JVS( 109 )
   W( 84 ) = JVS( 110 )
   W( 85 ) = JVS( 111 )
  JVS( 108) = W( 27 )
  JVS( 109) = W( 83 )
  JVS( 110) = W( 84 )
  JVS( 111) = W( 85 )
  IF ( ABS(  JVS( 112 )) < TINY(a) ) THEN
         IER = 28                                      
         RETURN
  END IF
   W( 28 ) = JVS( 112 )
   W( 52 ) = JVS( 113 )
   W( 84 ) = JVS( 114 )
   W( 88 ) = JVS( 115 )
  JVS( 112) = W( 28 )
  JVS( 113) = W( 52 )
  JVS( 114) = W( 84 )
  JVS( 115) = W( 88 )
  IF ( ABS(  JVS( 116 )) < TINY(a) ) THEN
         IER = 29                                      
         RETURN
  END IF
   W( 29 ) = JVS( 116 )
   W( 84 ) = JVS( 117 )
   W( 85 ) = JVS( 118 )
   W( 88 ) = JVS( 119 )
  JVS( 116) = W( 29 )
  JVS( 117) = W( 84 )
  JVS( 118) = W( 85 )
  JVS( 119) = W( 88 )
  IF ( ABS(  JVS( 121 )) < TINY(a) ) THEN
         IER = 30                                      
         RETURN
  END IF
   W( 28 ) = JVS( 120 )
   W( 30 ) = JVS( 121 )
   W( 41 ) = JVS( 122 )
   W( 52 ) = JVS( 123 )
   W( 84 ) = JVS( 124 )
   W( 87 ) = JVS( 125 )
   W( 88 ) = JVS( 126 )
  a = -W( 28 ) / JVS(          112  )
  W( 28 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 113 )
  W( 84 ) = W( 84 ) + a*JVS( 114 )
  W( 88 ) = W( 88 ) + a*JVS( 115 )
  JVS( 120) = W( 28 )
  JVS( 121) = W( 30 )
  JVS( 122) = W( 41 )
  JVS( 123) = W( 52 )
  JVS( 124) = W( 84 )
  JVS( 125) = W( 87 )
  JVS( 126) = W( 88 )
  IF ( ABS(  JVS( 127 )) < TINY(a) ) THEN
         IER = 31                                      
         RETURN
  END IF
   W( 31 ) = JVS( 127 )
   W( 36 ) = JVS( 128 )
   W( 37 ) = JVS( 129 )
   W( 39 ) = JVS( 130 )
   W( 84 ) = JVS( 131 )
   W( 87 ) = JVS( 132 )
   W( 88 ) = JVS( 133 )
  JVS( 127) = W( 31 )
  JVS( 128) = W( 36 )
  JVS( 129) = W( 37 )
  JVS( 130) = W( 39 )
  JVS( 131) = W( 84 )
  JVS( 132) = W( 87 )
  JVS( 133) = W( 88 )
  IF ( ABS(  JVS( 135 )) < TINY(a) ) THEN
         IER = 32                                      
         RETURN
  END IF
   W( 20 ) = JVS( 134 )
   W( 32 ) = JVS( 135 )
   W( 54 ) = JVS( 136 )
   W( 81 ) = JVS( 137 )
   W( 83 ) = JVS( 138 )
   W( 84 ) = JVS( 139 )
  a = -W( 20 ) / JVS(           82  )
  W( 20 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 83 )
  JVS( 134) = W( 20 )
  JVS( 135) = W( 32 )
  JVS( 136) = W( 54 )
  JVS( 137) = W( 81 )
  JVS( 138) = W( 83 )
  JVS( 139) = W( 84 )
  IF ( ABS(  JVS( 141 )) < TINY(a) ) THEN
         IER = 33                                      
         RETURN
  END IF
   W( 13 ) = JVS( 140 )
   W( 33 ) = JVS( 141 )
   W( 51 ) = JVS( 142 )
   W( 79 ) = JVS( 143 )
   W( 87 ) = JVS( 144 )
   W( 88 ) = JVS( 145 )
   W( 89 ) = JVS( 146 )
  a = -W( 13 ) / JVS(           65  )
  W( 13 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 66 )
  JVS( 140) = W( 13 )
  JVS( 141) = W( 33 )
  JVS( 142) = W( 51 )
  JVS( 143) = W( 79 )
  JVS( 144) = W( 87 )
  JVS( 145) = W( 88 )
  JVS( 146) = W( 89 )
  IF ( ABS(  JVS( 147 )) < TINY(a) ) THEN
         IER = 34                                      
         RETURN
  END IF
   W( 34 ) = JVS( 147 )
   W( 50 ) = JVS( 148 )
   W( 84 ) = JVS( 149 )
   W( 85 ) = JVS( 150 )
   W( 88 ) = JVS( 151 )
   W( 89 ) = JVS( 152 )
  JVS( 147) = W( 34 )
  JVS( 148) = W( 50 )
  JVS( 149) = W( 84 )
  JVS( 150) = W( 85 )
  JVS( 151) = W( 88 )
  JVS( 152) = W( 89 )
  IF ( ABS(  JVS( 153 )) < TINY(a) ) THEN
         IER = 35                                      
         RETURN
  END IF
   W( 35 ) = JVS( 153 )
   W( 45 ) = JVS( 154 )
   W( 47 ) = JVS( 155 )
   W( 48 ) = JVS( 156 )
   W( 49 ) = JVS( 157 )
   W( 51 ) = JVS( 158 )
   W( 73 ) = JVS( 159 )
   W( 75 ) = JVS( 160 )
   W( 79 ) = JVS( 161 )
   W( 84 ) = JVS( 162 )
   W( 85 ) = JVS( 163 )
  JVS( 153) = W( 35 )
  JVS( 154) = W( 45 )
  JVS( 155) = W( 47 )
  JVS( 156) = W( 48 )
  JVS( 157) = W( 49 )
  JVS( 158) = W( 51 )
  JVS( 159) = W( 73 )
  JVS( 160) = W( 75 )
  JVS( 161) = W( 79 )
  JVS( 162) = W( 84 )
  JVS( 163) = W( 85 )
  IF ( ABS(  JVS( 165 )) < TINY(a) ) THEN
         IER = 36                                      
         RETURN
  END IF
   W( 16 ) = JVS( 164 )
   W( 36 ) = JVS( 165 )
   W( 79 ) = JVS( 166 )
   W( 84 ) = JVS( 167 )
   W( 88 ) = JVS( 168 )
  a = -W( 16 ) / JVS(           73  )
  W( 16 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 74 )
  JVS( 164) = W( 16 )
  JVS( 165) = W( 36 )
  JVS( 166) = W( 79 )
  JVS( 167) = W( 84 )
  JVS( 168) = W( 88 )
  IF ( ABS(  JVS( 170 )) < TINY(a) ) THEN
         IER = 37                                      
         RETURN
  END IF
   W( 17 ) = JVS( 169 )
   W( 37 ) = JVS( 170 )
   W( 79 ) = JVS( 171 )
   W( 84 ) = JVS( 172 )
   W( 88 ) = JVS( 173 )
  a = -W( 17 ) / JVS(           75  )
  W( 17 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 76 )
  JVS( 169) = W( 17 )
  JVS( 170) = W( 37 )
  JVS( 171) = W( 79 )
  JVS( 172) = W( 84 )
  JVS( 173) = W( 88 )
  IF ( ABS(  JVS( 174 )) < TINY(a) ) THEN
         IER = 38                                      
         RETURN
  END IF
   W( 38 ) = JVS( 174 )
   W( 79 ) = JVS( 175 )
   W( 84 ) = JVS( 176 )
   W( 89 ) = JVS( 177 )
  JVS( 174) = W( 38 )
  JVS( 175) = W( 79 )
  JVS( 176) = W( 84 )
  JVS( 177) = W( 89 )
  IF ( ABS(  JVS( 178 )) < TINY(a) ) THEN
         IER = 39                                      
         RETURN
  END IF
   W( 39 ) = JVS( 178 )
   W( 50 ) = JVS( 179 )
   W( 79 ) = JVS( 180 )
   W( 84 ) = JVS( 181 )
   W( 88 ) = JVS( 182 )
  JVS( 178) = W( 39 )
  JVS( 179) = W( 50 )
  JVS( 180) = W( 79 )
  JVS( 181) = W( 84 )
  JVS( 182) = W( 88 )
  IF ( ABS(  JVS( 183 )) < TINY(a) ) THEN
         IER = 40                                      
         RETURN
  END IF
   W( 40 ) = JVS( 183 )
   W( 71 ) = JVS( 184 )
   W( 79 ) = JVS( 185 )
   W( 82 ) = JVS( 186 )
   W( 84 ) = JVS( 187 )
   W( 85 ) = JVS( 188 )
  JVS( 183) = W( 40 )
  JVS( 184) = W( 71 )
  JVS( 185) = W( 79 )
  JVS( 186) = W( 82 )
  JVS( 187) = W( 84 )
  JVS( 188) = W( 85 )
  IF ( ABS(  JVS( 189 )) < TINY(a) ) THEN
         IER = 41                                      
         RETURN
  END IF
   W( 41 ) = JVS( 189 )
   W( 48 ) = JVS( 190 )
   W( 60 ) = JVS( 191 )
   W( 84 ) = JVS( 192 )
   W( 87 ) = JVS( 193 )
   W( 89 ) = JVS( 194 )
  JVS( 189) = W( 41 )
  JVS( 190) = W( 48 )
  JVS( 191) = W( 60 )
  JVS( 192) = W( 84 )
  JVS( 193) = W( 87 )
  JVS( 194) = W( 89 )
  IF ( ABS(  JVS( 195 )) < TINY(a) ) THEN
         IER = 42                                      
         RETURN
  END IF
   W( 42 ) = JVS( 195 )
   W( 79 ) = JVS( 196 )
   W( 84 ) = JVS( 197 )
   W( 89 ) = JVS( 198 )
  JVS( 195) = W( 42 )
  JVS( 196) = W( 79 )
  JVS( 197) = W( 84 )
  JVS( 198) = W( 89 )
  IF ( ABS(  JVS( 199 )) < TINY(a) ) THEN
         IER = 43                                      
         RETURN
  END IF
   W( 43 ) = JVS( 199 )
   W( 50 ) = JVS( 200 )
   W( 53 ) = JVS( 201 )
   W( 62 ) = JVS( 202 )
   W( 66 ) = JVS( 203 )
   W( 71 ) = JVS( 204 )
   W( 74 ) = JVS( 205 )
   W( 84 ) = JVS( 206 )
   W( 85 ) = JVS( 207 )
   W( 88 ) = JVS( 208 )
   W( 89 ) = JVS( 209 )
  JVS( 199) = W( 43 )
  JVS( 200) = W( 50 )
  JVS( 201) = W( 53 )
  JVS( 202) = W( 62 )
  JVS( 203) = W( 66 )
  JVS( 204) = W( 71 )
  JVS( 205) = W( 74 )
  JVS( 206) = W( 84 )
  JVS( 207) = W( 85 )
  JVS( 208) = W( 88 )
  JVS( 209) = W( 89 )
  IF ( ABS(  JVS( 210 )) < TINY(a) ) THEN
         IER = 44                                      
         RETURN
  END IF
   W( 44 ) = JVS( 210 )
   W( 54 ) = JVS( 211 )
   W( 79 ) = JVS( 212 )
   W( 82 ) = JVS( 213 )
   W( 84 ) = JVS( 214 )
   W( 88 ) = JVS( 215 )
   W( 89 ) = JVS( 216 )
  JVS( 210) = W( 44 )
  JVS( 211) = W( 54 )
  JVS( 212) = W( 79 )
  JVS( 213) = W( 82 )
  JVS( 214) = W( 84 )
  JVS( 215) = W( 88 )
  JVS( 216) = W( 89 )
  IF ( ABS(  JVS( 217 )) < TINY(a) ) THEN
         IER = 45                                      
         RETURN
  END IF
   W( 45 ) = JVS( 217 )
   W( 79 ) = JVS( 218 )
   W( 84 ) = JVS( 219 )
   W( 89 ) = JVS( 220 )
  JVS( 217) = W( 45 )
  JVS( 218) = W( 79 )
  JVS( 219) = W( 84 )
  JVS( 220) = W( 89 )
  IF ( ABS(  JVS( 226 )) < TINY(a) ) THEN
         IER = 46                                      
         RETURN
  END IF
   W( 22 ) = JVS( 221 )
   W( 24 ) = JVS( 222 )
   W( 38 ) = JVS( 223 )
   W( 41 ) = JVS( 224 )
   W( 45 ) = JVS( 225 )
   W( 46 ) = JVS( 226 )
   W( 47 ) = JVS( 227 )
   W( 48 ) = JVS( 228 )
   W( 49 ) = JVS( 229 )
   W( 51 ) = JVS( 230 )
   W( 52 ) = JVS( 231 )
   W( 53 ) = JVS( 232 )
   W( 54 ) = JVS( 233 )
   W( 57 ) = JVS( 234 )
   W( 60 ) = JVS( 235 )
   W( 62 ) = JVS( 236 )
   W( 66 ) = JVS( 237 )
   W( 71 ) = JVS( 238 )
   W( 73 ) = JVS( 239 )
   W( 74 ) = JVS( 240 )
   W( 75 ) = JVS( 241 )
   W( 79 ) = JVS( 242 )
   W( 84 ) = JVS( 243 )
   W( 87 ) = JVS( 244 )
   W( 89 ) = JVS( 245 )
  a = -W( 22 ) / JVS(           88  )
  W( 22 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 89 )
  a = -W( 24 ) / JVS(           93  )
  W( 24 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 94 )
  W( 84 ) = W( 84 ) + a*JVS( 95 )
  a = -W( 38 ) / JVS(          174  )
  W( 38 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 175 )
  W( 84 ) = W( 84 ) + a*JVS( 176 )
  W( 89 ) = W( 89 ) + a*JVS( 177 )
  a = -W( 41 ) / JVS(          189  )
  W( 41 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 190 )
  W( 60 ) = W( 60 ) + a*JVS( 191 )
  W( 84 ) = W( 84 ) + a*JVS( 192 )
  W( 87 ) = W( 87 ) + a*JVS( 193 )
  W( 89 ) = W( 89 ) + a*JVS( 194 )
  a = -W( 45 ) / JVS(          217  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 218 )
  W( 84 ) = W( 84 ) + a*JVS( 219 )
  W( 89 ) = W( 89 ) + a*JVS( 220 )
  JVS( 221) = W( 22 )
  JVS( 222) = W( 24 )
  JVS( 223) = W( 38 )
  JVS( 224) = W( 41 )
  JVS( 225) = W( 45 )
  JVS( 226) = W( 46 )
  JVS( 227) = W( 47 )
  JVS( 228) = W( 48 )
  JVS( 229) = W( 49 )
  JVS( 230) = W( 51 )
  JVS( 231) = W( 52 )
  JVS( 232) = W( 53 )
  JVS( 233) = W( 54 )
  JVS( 234) = W( 57 )
  JVS( 235) = W( 60 )
  JVS( 236) = W( 62 )
  JVS( 237) = W( 66 )
  JVS( 238) = W( 71 )
  JVS( 239) = W( 73 )
  JVS( 240) = W( 74 )
  JVS( 241) = W( 75 )
  JVS( 242) = W( 79 )
  JVS( 243) = W( 84 )
  JVS( 244) = W( 87 )
  JVS( 245) = W( 89 )
  IF ( ABS(  JVS( 246 )) < TINY(a) ) THEN
         IER = 47                                      
         RETURN
  END IF
   W( 47 ) = JVS( 246 )
   W( 79 ) = JVS( 247 )
   W( 84 ) = JVS( 248 )
   W( 89 ) = JVS( 249 )
  JVS( 246) = W( 47 )
  JVS( 247) = W( 79 )
  JVS( 248) = W( 84 )
  JVS( 249) = W( 89 )
  IF ( ABS(  JVS( 250 )) < TINY(a) ) THEN
         IER = 48                                      
         RETURN
  END IF
   W( 48 ) = JVS( 250 )
   W( 79 ) = JVS( 251 )
   W( 84 ) = JVS( 252 )
   W( 89 ) = JVS( 253 )
  JVS( 250) = W( 48 )
  JVS( 251) = W( 79 )
  JVS( 252) = W( 84 )
  JVS( 253) = W( 89 )
  IF ( ABS(  JVS( 254 )) < TINY(a) ) THEN
         IER = 49                                      
         RETURN
  END IF
   W( 49 ) = JVS( 254 )
   W( 79 ) = JVS( 255 )
   W( 84 ) = JVS( 256 )
   W( 89 ) = JVS( 257 )
  JVS( 254) = W( 49 )
  JVS( 255) = W( 79 )
  JVS( 256) = W( 84 )
  JVS( 257) = W( 89 )
  IF ( ABS(  JVS( 262 )) < TINY(a) ) THEN
         IER = 50                                      
         RETURN
  END IF
   W( 34 ) = JVS( 258 )
   W( 36 ) = JVS( 259 )
   W( 37 ) = JVS( 260 )
   W( 39 ) = JVS( 261 )
   W( 50 ) = JVS( 262 )
   W( 79 ) = JVS( 263 )
   W( 84 ) = JVS( 264 )
   W( 85 ) = JVS( 265 )
   W( 88 ) = JVS( 266 )
   W( 89 ) = JVS( 267 )
  a = -W( 34 ) / JVS(          147  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 148 )
  W( 84 ) = W( 84 ) + a*JVS( 149 )
  W( 85 ) = W( 85 ) + a*JVS( 150 )
  W( 88 ) = W( 88 ) + a*JVS( 151 )
  W( 89 ) = W( 89 ) + a*JVS( 152 )
  a = -W( 36 ) / JVS(          165  )
  W( 36 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 166 )
  W( 84 ) = W( 84 ) + a*JVS( 167 )
  W( 88 ) = W( 88 ) + a*JVS( 168 )
  a = -W( 37 ) / JVS(          170  )
  W( 37 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 171 )
  W( 84 ) = W( 84 ) + a*JVS( 172 )
  W( 88 ) = W( 88 ) + a*JVS( 173 )
  a = -W( 39 ) / JVS(          178  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 179 )
  W( 79 ) = W( 79 ) + a*JVS( 180 )
  W( 84 ) = W( 84 ) + a*JVS( 181 )
  W( 88 ) = W( 88 ) + a*JVS( 182 )
  JVS( 258) = W( 34 )
  JVS( 259) = W( 36 )
  JVS( 260) = W( 37 )
  JVS( 261) = W( 39 )
  JVS( 262) = W( 50 )
  JVS( 263) = W( 79 )
  JVS( 264) = W( 84 )
  JVS( 265) = W( 85 )
  JVS( 266) = W( 88 )
  JVS( 267) = W( 89 )
  IF ( ABS(  JVS( 268 )) < TINY(a) ) THEN
         IER = 51                                      
         RETURN
  END IF
   W( 51 ) = JVS( 268 )
   W( 79 ) = JVS( 269 )
   W( 84 ) = JVS( 270 )
   W( 89 ) = JVS( 271 )
  JVS( 268) = W( 51 )
  JVS( 269) = W( 79 )
  JVS( 270) = W( 84 )
  JVS( 271) = W( 89 )
  IF ( ABS(  JVS( 275 )) < TINY(a) ) THEN
         IER = 52                                      
         RETURN
  END IF
   W( 21 ) = JVS( 272 )
   W( 28 ) = JVS( 273 )
   W( 48 ) = JVS( 274 )
   W( 52 ) = JVS( 275 )
   W( 57 ) = JVS( 276 )
   W( 79 ) = JVS( 277 )
   W( 84 ) = JVS( 278 )
   W( 85 ) = JVS( 279 )
   W( 87 ) = JVS( 280 )
   W( 88 ) = JVS( 281 )
   W( 89 ) = JVS( 282 )
  a = -W( 21 ) / JVS(           84  )
  W( 21 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 85 )
  W( 84 ) = W( 84 ) + a*JVS( 86 )
  W( 85 ) = W( 85 ) + a*JVS( 87 )
  a = -W( 28 ) / JVS(          112  )
  W( 28 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 113 )
  W( 84 ) = W( 84 ) + a*JVS( 114 )
  W( 88 ) = W( 88 ) + a*JVS( 115 )
  a = -W( 48 ) / JVS(          250  )
  W( 48 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 251 )
  W( 84 ) = W( 84 ) + a*JVS( 252 )
  W( 89 ) = W( 89 ) + a*JVS( 253 )
  JVS( 272) = W( 21 )
  JVS( 273) = W( 28 )
  JVS( 274) = W( 48 )
  JVS( 275) = W( 52 )
  JVS( 276) = W( 57 )
  JVS( 277) = W( 79 )
  JVS( 278) = W( 84 )
  JVS( 279) = W( 85 )
  JVS( 280) = W( 87 )
  JVS( 281) = W( 88 )
  JVS( 282) = W( 89 )
  IF ( ABS(  JVS( 284 )) < TINY(a) ) THEN
         IER = 53                                      
         RETURN
  END IF
   W( 22 ) = JVS( 283 )
   W( 53 ) = JVS( 284 )
   W( 61 ) = JVS( 285 )
   W( 67 ) = JVS( 286 )
   W( 68 ) = JVS( 287 )
   W( 71 ) = JVS( 288 )
   W( 79 ) = JVS( 289 )
   W( 80 ) = JVS( 290 )
   W( 82 ) = JVS( 291 )
   W( 83 ) = JVS( 292 )
   W( 84 ) = JVS( 293 )
   W( 87 ) = JVS( 294 )
   W( 89 ) = JVS( 295 )
  a = -W( 22 ) / JVS(           88  )
  W( 22 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 89 )
  JVS( 283) = W( 22 )
  JVS( 284) = W( 53 )
  JVS( 285) = W( 61 )
  JVS( 286) = W( 67 )
  JVS( 287) = W( 68 )
  JVS( 288) = W( 71 )
  JVS( 289) = W( 79 )
  JVS( 290) = W( 80 )
  JVS( 291) = W( 82 )
  JVS( 292) = W( 83 )
  JVS( 293) = W( 84 )
  JVS( 294) = W( 87 )
  JVS( 295) = W( 89 )
  IF ( ABS(  JVS( 296 )) < TINY(a) ) THEN
         IER = 54                                      
         RETURN
  END IF
   W( 54 ) = JVS( 296 )
   W( 78 ) = JVS( 297 )
   W( 79 ) = JVS( 298 )
   W( 84 ) = JVS( 299 )
   W( 88 ) = JVS( 300 )
   W( 89 ) = JVS( 301 )
  JVS( 296) = W( 54 )
  JVS( 297) = W( 78 )
  JVS( 298) = W( 79 )
  JVS( 299) = W( 84 )
  JVS( 300) = W( 88 )
  JVS( 301) = W( 89 )
  IF ( ABS(  JVS( 303 )) < TINY(a) ) THEN
         IER = 55                                      
         RETURN
  END IF
   W( 38 ) = JVS( 302 )
   W( 55 ) = JVS( 303 )
   W( 79 ) = JVS( 304 )
   W( 82 ) = JVS( 305 )
   W( 83 ) = JVS( 306 )
   W( 84 ) = JVS( 307 )
   W( 85 ) = JVS( 308 )
   W( 87 ) = JVS( 309 )
   W( 89 ) = JVS( 310 )
  a = -W( 38 ) / JVS(          174  )
  W( 38 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 175 )
  W( 84 ) = W( 84 ) + a*JVS( 176 )
  W( 89 ) = W( 89 ) + a*JVS( 177 )
  JVS( 302) = W( 38 )
  JVS( 303) = W( 55 )
  JVS( 304) = W( 79 )
  JVS( 305) = W( 82 )
  JVS( 306) = W( 83 )
  JVS( 307) = W( 84 )
  JVS( 308) = W( 85 )
  JVS( 309) = W( 87 )
  JVS( 310) = W( 89 )
  IF ( ABS(  JVS( 311 )) < TINY(a) ) THEN
         IER = 56                                      
         RETURN
  END IF
   W( 56 ) = JVS( 311 )
   W( 75 ) = JVS( 312 )
   W( 82 ) = JVS( 313 )
   W( 83 ) = JVS( 314 )
   W( 84 ) = JVS( 315 )
   W( 85 ) = JVS( 316 )
   W( 87 ) = JVS( 317 )
   W( 89 ) = JVS( 318 )
  JVS( 311) = W( 56 )
  JVS( 312) = W( 75 )
  JVS( 313) = W( 82 )
  JVS( 314) = W( 83 )
  JVS( 315) = W( 84 )
  JVS( 316) = W( 85 )
  JVS( 317) = W( 87 )
  JVS( 318) = W( 89 )
  IF ( ABS(  JVS( 324 )) < TINY(a) ) THEN
         IER = 57                                      
         RETURN
  END IF
   W( 15 ) = JVS( 319 )
   W( 42 ) = JVS( 320 )
   W( 47 ) = JVS( 321 )
   W( 48 ) = JVS( 322 )
   W( 51 ) = JVS( 323 )
   W( 57 ) = JVS( 324 )
   W( 60 ) = JVS( 325 )
   W( 63 ) = JVS( 326 )
   W( 79 ) = JVS( 327 )
   W( 82 ) = JVS( 328 )
   W( 83 ) = JVS( 329 )
   W( 84 ) = JVS( 330 )
   W( 85 ) = JVS( 331 )
   W( 87 ) = JVS( 332 )
   W( 89 ) = JVS( 333 )
  a = -W( 15 ) / JVS(           69  )
  W( 15 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 70 )
  W( 84 ) = W( 84 ) + a*JVS( 71 )
  W( 85 ) = W( 85 ) + a*JVS( 72 )
  a = -W( 42 ) / JVS(          195  )
  W( 42 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 196 )
  W( 84 ) = W( 84 ) + a*JVS( 197 )
  W( 89 ) = W( 89 ) + a*JVS( 198 )
  a = -W( 47 ) / JVS(          246  )
  W( 47 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 247 )
  W( 84 ) = W( 84 ) + a*JVS( 248 )
  W( 89 ) = W( 89 ) + a*JVS( 249 )
  a = -W( 48 ) / JVS(          250  )
  W( 48 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 251 )
  W( 84 ) = W( 84 ) + a*JVS( 252 )
  W( 89 ) = W( 89 ) + a*JVS( 253 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  JVS( 319) = W( 15 )
  JVS( 320) = W( 42 )
  JVS( 321) = W( 47 )
  JVS( 322) = W( 48 )
  JVS( 323) = W( 51 )
  JVS( 324) = W( 57 )
  JVS( 325) = W( 60 )
  JVS( 326) = W( 63 )
  JVS( 327) = W( 79 )
  JVS( 328) = W( 82 )
  JVS( 329) = W( 83 )
  JVS( 330) = W( 84 )
  JVS( 331) = W( 85 )
  JVS( 332) = W( 87 )
  JVS( 333) = W( 89 )
  IF ( ABS(  JVS( 335 )) < TINY(a) ) THEN
         IER = 58                                      
         RETURN
  END IF
   W( 45 ) = JVS( 334 )
   W( 58 ) = JVS( 335 )
   W( 79 ) = JVS( 336 )
   W( 82 ) = JVS( 337 )
   W( 83 ) = JVS( 338 )
   W( 84 ) = JVS( 339 )
   W( 85 ) = JVS( 340 )
   W( 87 ) = JVS( 341 )
   W( 89 ) = JVS( 342 )
  a = -W( 45 ) / JVS(          217  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 218 )
  W( 84 ) = W( 84 ) + a*JVS( 219 )
  W( 89 ) = W( 89 ) + a*JVS( 220 )
  JVS( 334) = W( 45 )
  JVS( 335) = W( 58 )
  JVS( 336) = W( 79 )
  JVS( 337) = W( 82 )
  JVS( 338) = W( 83 )
  JVS( 339) = W( 84 )
  JVS( 340) = W( 85 )
  JVS( 341) = W( 87 )
  JVS( 342) = W( 89 )
  IF ( ABS(  JVS( 350 )) < TINY(a) ) THEN
         IER = 59                                      
         RETURN
  END IF
   W( 18 ) = JVS( 343 )
   W( 23 ) = JVS( 344 )
   W( 42 ) = JVS( 345 )
   W( 45 ) = JVS( 346 )
   W( 49 ) = JVS( 347 )
   W( 56 ) = JVS( 348 )
   W( 58 ) = JVS( 349 )
   W( 59 ) = JVS( 350 )
   W( 64 ) = JVS( 351 )
   W( 65 ) = JVS( 352 )
   W( 69 ) = JVS( 353 )
   W( 70 ) = JVS( 354 )
   W( 71 ) = JVS( 355 )
   W( 73 ) = JVS( 356 )
   W( 75 ) = JVS( 357 )
   W( 76 ) = JVS( 358 )
   W( 77 ) = JVS( 359 )
   W( 79 ) = JVS( 360 )
   W( 80 ) = JVS( 361 )
   W( 81 ) = JVS( 362 )
   W( 82 ) = JVS( 363 )
   W( 83 ) = JVS( 364 )
   W( 84 ) = JVS( 365 )
   W( 85 ) = JVS( 366 )
   W( 86 ) = JVS( 367 )
   W( 87 ) = JVS( 368 )
   W( 89 ) = JVS( 369 )
  a = -W( 18 ) / JVS(           77  )
  W( 18 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 78 )
  a = -W( 23 ) / JVS(           90  )
  W( 23 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 91 )
  W( 84 ) = W( 84 ) + a*JVS( 92 )
  a = -W( 42 ) / JVS(          195  )
  W( 42 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 196 )
  W( 84 ) = W( 84 ) + a*JVS( 197 )
  W( 89 ) = W( 89 ) + a*JVS( 198 )
  a = -W( 45 ) / JVS(          217  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 218 )
  W( 84 ) = W( 84 ) + a*JVS( 219 )
  W( 89 ) = W( 89 ) + a*JVS( 220 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 255 )
  W( 84 ) = W( 84 ) + a*JVS( 256 )
  W( 89 ) = W( 89 ) + a*JVS( 257 )
  a = -W( 56 ) / JVS(          311  )
  W( 56 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 312 )
  W( 82 ) = W( 82 ) + a*JVS( 313 )
  W( 83 ) = W( 83 ) + a*JVS( 314 )
  W( 84 ) = W( 84 ) + a*JVS( 315 )
  W( 85 ) = W( 85 ) + a*JVS( 316 )
  W( 87 ) = W( 87 ) + a*JVS( 317 )
  W( 89 ) = W( 89 ) + a*JVS( 318 )
  a = -W( 58 ) / JVS(          335  )
  W( 58 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 336 )
  W( 82 ) = W( 82 ) + a*JVS( 337 )
  W( 83 ) = W( 83 ) + a*JVS( 338 )
  W( 84 ) = W( 84 ) + a*JVS( 339 )
  W( 85 ) = W( 85 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 89 ) = W( 89 ) + a*JVS( 342 )
  JVS( 343) = W( 18 )
  JVS( 344) = W( 23 )
  JVS( 345) = W( 42 )
  JVS( 346) = W( 45 )
  JVS( 347) = W( 49 )
  JVS( 348) = W( 56 )
  JVS( 349) = W( 58 )
  JVS( 350) = W( 59 )
  JVS( 351) = W( 64 )
  JVS( 352) = W( 65 )
  JVS( 353) = W( 69 )
  JVS( 354) = W( 70 )
  JVS( 355) = W( 71 )
  JVS( 356) = W( 73 )
  JVS( 357) = W( 75 )
  JVS( 358) = W( 76 )
  JVS( 359) = W( 77 )
  JVS( 360) = W( 79 )
  JVS( 361) = W( 80 )
  JVS( 362) = W( 81 )
  JVS( 363) = W( 82 )
  JVS( 364) = W( 83 )
  JVS( 365) = W( 84 )
  JVS( 366) = W( 85 )
  JVS( 367) = W( 86 )
  JVS( 368) = W( 87 )
  JVS( 369) = W( 89 )
  IF ( ABS(  JVS( 372 )) < TINY(a) ) THEN
         IER = 60                                      
         RETURN
  END IF
   W( 48 ) = JVS( 370 )
   W( 51 ) = JVS( 371 )
   W( 60 ) = JVS( 372 )
   W( 79 ) = JVS( 373 )
   W( 82 ) = JVS( 374 )
   W( 83 ) = JVS( 375 )
   W( 84 ) = JVS( 376 )
   W( 85 ) = JVS( 377 )
   W( 87 ) = JVS( 378 )
   W( 89 ) = JVS( 379 )
  a = -W( 48 ) / JVS(          250  )
  W( 48 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 251 )
  W( 84 ) = W( 84 ) + a*JVS( 252 )
  W( 89 ) = W( 89 ) + a*JVS( 253 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  JVS( 370) = W( 48 )
  JVS( 371) = W( 51 )
  JVS( 372) = W( 60 )
  JVS( 373) = W( 79 )
  JVS( 374) = W( 82 )
  JVS( 375) = W( 83 )
  JVS( 376) = W( 84 )
  JVS( 377) = W( 85 )
  JVS( 378) = W( 87 )
  JVS( 379) = W( 89 )
  IF ( ABS(  JVS( 382 )) < TINY(a) ) THEN
         IER = 61                                      
         RETURN
  END IF
   W( 39 ) = JVS( 380 )
   W( 50 ) = JVS( 381 )
   W( 61 ) = JVS( 382 )
   W( 79 ) = JVS( 383 )
   W( 82 ) = JVS( 384 )
   W( 83 ) = JVS( 385 )
   W( 84 ) = JVS( 386 )
   W( 85 ) = JVS( 387 )
   W( 87 ) = JVS( 388 )
   W( 88 ) = JVS( 389 )
   W( 89 ) = JVS( 390 )
  a = -W( 39 ) / JVS(          178  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 179 )
  W( 79 ) = W( 79 ) + a*JVS( 180 )
  W( 84 ) = W( 84 ) + a*JVS( 181 )
  W( 88 ) = W( 88 ) + a*JVS( 182 )
  a = -W( 50 ) / JVS(          262  )
  W( 50 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 263 )
  W( 84 ) = W( 84 ) + a*JVS( 264 )
  W( 85 ) = W( 85 ) + a*JVS( 265 )
  W( 88 ) = W( 88 ) + a*JVS( 266 )
  W( 89 ) = W( 89 ) + a*JVS( 267 )
  JVS( 380) = W( 39 )
  JVS( 381) = W( 50 )
  JVS( 382) = W( 61 )
  JVS( 383) = W( 79 )
  JVS( 384) = W( 82 )
  JVS( 385) = W( 83 )
  JVS( 386) = W( 84 )
  JVS( 387) = W( 85 )
  JVS( 388) = W( 87 )
  JVS( 389) = W( 88 )
  JVS( 390) = W( 89 )
  IF ( ABS(  JVS( 400 )) < TINY(a) ) THEN
         IER = 62                                      
         RETURN
  END IF
   W( 30 ) = JVS( 391 )
   W( 32 ) = JVS( 392 )
   W( 41 ) = JVS( 393 )
   W( 48 ) = JVS( 394 )
   W( 52 ) = JVS( 395 )
   W( 54 ) = JVS( 396 )
   W( 57 ) = JVS( 397 )
   W( 60 ) = JVS( 398 )
   W( 61 ) = JVS( 399 )
   W( 62 ) = JVS( 400 )
   W( 63 ) = JVS( 401 )
   W( 67 ) = JVS( 402 )
   W( 68 ) = JVS( 403 )
   W( 71 ) = JVS( 404 )
   W( 78 ) = JVS( 405 )
   W( 79 ) = JVS( 406 )
   W( 80 ) = JVS( 407 )
   W( 81 ) = JVS( 408 )
   W( 82 ) = JVS( 409 )
   W( 83 ) = JVS( 410 )
   W( 84 ) = JVS( 411 )
   W( 85 ) = JVS( 412 )
   W( 87 ) = JVS( 413 )
   W( 88 ) = JVS( 414 )
   W( 89 ) = JVS( 415 )
  a = -W( 30 ) / JVS(          121  )
  W( 30 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 122 )
  W( 52 ) = W( 52 ) + a*JVS( 123 )
  W( 84 ) = W( 84 ) + a*JVS( 124 )
  W( 87 ) = W( 87 ) + a*JVS( 125 )
  W( 88 ) = W( 88 ) + a*JVS( 126 )
  a = -W( 32 ) / JVS(          135  )
  W( 32 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 136 )
  W( 81 ) = W( 81 ) + a*JVS( 137 )
  W( 83 ) = W( 83 ) + a*JVS( 138 )
  W( 84 ) = W( 84 ) + a*JVS( 139 )
  a = -W( 41 ) / JVS(          189  )
  W( 41 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 190 )
  W( 60 ) = W( 60 ) + a*JVS( 191 )
  W( 84 ) = W( 84 ) + a*JVS( 192 )
  W( 87 ) = W( 87 ) + a*JVS( 193 )
  W( 89 ) = W( 89 ) + a*JVS( 194 )
  a = -W( 48 ) / JVS(          250  )
  W( 48 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 251 )
  W( 84 ) = W( 84 ) + a*JVS( 252 )
  W( 89 ) = W( 89 ) + a*JVS( 253 )
  a = -W( 52 ) / JVS(          275  )
  W( 52 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 276 )
  W( 79 ) = W( 79 ) + a*JVS( 277 )
  W( 84 ) = W( 84 ) + a*JVS( 278 )
  W( 85 ) = W( 85 ) + a*JVS( 279 )
  W( 87 ) = W( 87 ) + a*JVS( 280 )
  W( 88 ) = W( 88 ) + a*JVS( 281 )
  W( 89 ) = W( 89 ) + a*JVS( 282 )
  a = -W( 54 ) / JVS(          296  )
  W( 54 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 297 )
  W( 79 ) = W( 79 ) + a*JVS( 298 )
  W( 84 ) = W( 84 ) + a*JVS( 299 )
  W( 88 ) = W( 88 ) + a*JVS( 300 )
  W( 89 ) = W( 89 ) + a*JVS( 301 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 325 )
  W( 63 ) = W( 63 ) + a*JVS( 326 )
  W( 79 ) = W( 79 ) + a*JVS( 327 )
  W( 82 ) = W( 82 ) + a*JVS( 328 )
  W( 83 ) = W( 83 ) + a*JVS( 329 )
  W( 84 ) = W( 84 ) + a*JVS( 330 )
  W( 85 ) = W( 85 ) + a*JVS( 331 )
  W( 87 ) = W( 87 ) + a*JVS( 332 )
  W( 89 ) = W( 89 ) + a*JVS( 333 )
  a = -W( 60 ) / JVS(          372  )
  W( 60 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 82 ) = W( 82 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 85 ) = W( 85 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  a = -W( 61 ) / JVS(          382  )
  W( 61 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 82 ) = W( 82 ) + a*JVS( 384 )
  W( 83 ) = W( 83 ) + a*JVS( 385 )
  W( 84 ) = W( 84 ) + a*JVS( 386 )
  W( 85 ) = W( 85 ) + a*JVS( 387 )
  W( 87 ) = W( 87 ) + a*JVS( 388 )
  W( 88 ) = W( 88 ) + a*JVS( 389 )
  W( 89 ) = W( 89 ) + a*JVS( 390 )
  JVS( 391) = W( 30 )
  JVS( 392) = W( 32 )
  JVS( 393) = W( 41 )
  JVS( 394) = W( 48 )
  JVS( 395) = W( 52 )
  JVS( 396) = W( 54 )
  JVS( 397) = W( 57 )
  JVS( 398) = W( 60 )
  JVS( 399) = W( 61 )
  JVS( 400) = W( 62 )
  JVS( 401) = W( 63 )
  JVS( 402) = W( 67 )
  JVS( 403) = W( 68 )
  JVS( 404) = W( 71 )
  JVS( 405) = W( 78 )
  JVS( 406) = W( 79 )
  JVS( 407) = W( 80 )
  JVS( 408) = W( 81 )
  JVS( 409) = W( 82 )
  JVS( 410) = W( 83 )
  JVS( 411) = W( 84 )
  JVS( 412) = W( 85 )
  JVS( 413) = W( 87 )
  JVS( 414) = W( 88 )
  JVS( 415) = W( 89 )
  IF ( ABS(  JVS( 417 )) < TINY(a) ) THEN
         IER = 63                                      
         RETURN
  END IF
   W( 47 ) = JVS( 416 )
   W( 63 ) = JVS( 417 )
   W( 79 ) = JVS( 418 )
   W( 82 ) = JVS( 419 )
   W( 83 ) = JVS( 420 )
   W( 84 ) = JVS( 421 )
   W( 85 ) = JVS( 422 )
   W( 87 ) = JVS( 423 )
   W( 89 ) = JVS( 424 )
  a = -W( 47 ) / JVS(          246  )
  W( 47 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 247 )
  W( 84 ) = W( 84 ) + a*JVS( 248 )
  W( 89 ) = W( 89 ) + a*JVS( 249 )
  JVS( 416) = W( 47 )
  JVS( 417) = W( 63 )
  JVS( 418) = W( 79 )
  JVS( 419) = W( 82 )
  JVS( 420) = W( 83 )
  JVS( 421) = W( 84 )
  JVS( 422) = W( 85 )
  JVS( 423) = W( 87 )
  JVS( 424) = W( 89 )
  IF ( ABS(  JVS( 426 )) < TINY(a) ) THEN
         IER = 64                                      
         RETURN
  END IF
   W( 18 ) = JVS( 425 )
   W( 64 ) = JVS( 426 )
   W( 82 ) = JVS( 427 )
   W( 83 ) = JVS( 428 )
   W( 84 ) = JVS( 429 )
   W( 85 ) = JVS( 430 )
   W( 87 ) = JVS( 431 )
   W( 89 ) = JVS( 432 )
  a = -W( 18 ) / JVS(           77  )
  W( 18 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 78 )
  JVS( 425) = W( 18 )
  JVS( 426) = W( 64 )
  JVS( 427) = W( 82 )
  JVS( 428) = W( 83 )
  JVS( 429) = W( 84 )
  JVS( 430) = W( 85 )
  JVS( 431) = W( 87 )
  JVS( 432) = W( 89 )
  IF ( ABS(  JVS( 434 )) < TINY(a) ) THEN
         IER = 65                                      
         RETURN
  END IF
   W( 20 ) = JVS( 433 )
   W( 65 ) = JVS( 434 )
   W( 82 ) = JVS( 435 )
   W( 83 ) = JVS( 436 )
   W( 84 ) = JVS( 437 )
   W( 85 ) = JVS( 438 )
   W( 87 ) = JVS( 439 )
   W( 89 ) = JVS( 440 )
  a = -W( 20 ) / JVS(           82  )
  W( 20 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 83 )
  JVS( 433) = W( 20 )
  JVS( 434) = W( 65 )
  JVS( 435) = W( 82 )
  JVS( 436) = W( 83 )
  JVS( 437) = W( 84 )
  JVS( 438) = W( 85 )
  JVS( 439) = W( 87 )
  JVS( 440) = W( 89 )
  IF ( ABS(  JVS( 466 )) < TINY(a) ) THEN
         IER = 66                                      
         RETURN
  END IF
   W( 22 ) = JVS( 441 )
   W( 24 ) = JVS( 442 )
   W( 27 ) = JVS( 443 )
   W( 32 ) = JVS( 444 )
   W( 38 ) = JVS( 445 )
   W( 40 ) = JVS( 446 )
   W( 41 ) = JVS( 447 )
   W( 42 ) = JVS( 448 )
   W( 44 ) = JVS( 449 )
   W( 47 ) = JVS( 450 )
   W( 48 ) = JVS( 451 )
   W( 49 ) = JVS( 452 )
   W( 51 ) = JVS( 453 )
   W( 52 ) = JVS( 454 )
   W( 53 ) = JVS( 455 )
   W( 54 ) = JVS( 456 )
   W( 55 ) = JVS( 457 )
   W( 56 ) = JVS( 458 )
   W( 57 ) = JVS( 459 )
   W( 58 ) = JVS( 460 )
   W( 60 ) = JVS( 461 )
   W( 61 ) = JVS( 462 )
   W( 63 ) = JVS( 463 )
   W( 64 ) = JVS( 464 )
   W( 65 ) = JVS( 465 )
   W( 66 ) = JVS( 466 )
   W( 67 ) = JVS( 467 )
   W( 68 ) = JVS( 468 )
   W( 69 ) = JVS( 469 )
   W( 71 ) = JVS( 470 )
   W( 72 ) = JVS( 471 )
   W( 73 ) = JVS( 472 )
   W( 75 ) = JVS( 473 )
   W( 76 ) = JVS( 474 )
   W( 77 ) = JVS( 475 )
   W( 78 ) = JVS( 476 )
   W( 79 ) = JVS( 477 )
   W( 80 ) = JVS( 478 )
   W( 81 ) = JVS( 479 )
   W( 82 ) = JVS( 480 )
   W( 83 ) = JVS( 481 )
   W( 84 ) = JVS( 482 )
   W( 85 ) = JVS( 483 )
   W( 87 ) = JVS( 484 )
   W( 88 ) = JVS( 485 )
   W( 89 ) = JVS( 486 )
   W( 90 ) = JVS( 487 )
  a = -W( 22 ) / JVS(           88  )
  W( 22 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 89 )
  a = -W( 24 ) / JVS(           93  )
  W( 24 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 94 )
  W( 84 ) = W( 84 ) + a*JVS( 95 )
  a = -W( 27 ) / JVS(          108  )
  W( 27 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 109 )
  W( 84 ) = W( 84 ) + a*JVS( 110 )
  W( 85 ) = W( 85 ) + a*JVS( 111 )
  a = -W( 32 ) / JVS(          135  )
  W( 32 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 136 )
  W( 81 ) = W( 81 ) + a*JVS( 137 )
  W( 83 ) = W( 83 ) + a*JVS( 138 )
  W( 84 ) = W( 84 ) + a*JVS( 139 )
  a = -W( 38 ) / JVS(          174  )
  W( 38 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 175 )
  W( 84 ) = W( 84 ) + a*JVS( 176 )
  W( 89 ) = W( 89 ) + a*JVS( 177 )
  a = -W( 40 ) / JVS(          183  )
  W( 40 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 184 )
  W( 79 ) = W( 79 ) + a*JVS( 185 )
  W( 82 ) = W( 82 ) + a*JVS( 186 )
  W( 84 ) = W( 84 ) + a*JVS( 187 )
  W( 85 ) = W( 85 ) + a*JVS( 188 )
  a = -W( 41 ) / JVS(          189  )
  W( 41 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 190 )
  W( 60 ) = W( 60 ) + a*JVS( 191 )
  W( 84 ) = W( 84 ) + a*JVS( 192 )
  W( 87 ) = W( 87 ) + a*JVS( 193 )
  W( 89 ) = W( 89 ) + a*JVS( 194 )
  a = -W( 42 ) / JVS(          195  )
  W( 42 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 196 )
  W( 84 ) = W( 84 ) + a*JVS( 197 )
  W( 89 ) = W( 89 ) + a*JVS( 198 )
  a = -W( 44 ) / JVS(          210  )
  W( 44 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 211 )
  W( 79 ) = W( 79 ) + a*JVS( 212 )
  W( 82 ) = W( 82 ) + a*JVS( 213 )
  W( 84 ) = W( 84 ) + a*JVS( 214 )
  W( 88 ) = W( 88 ) + a*JVS( 215 )
  W( 89 ) = W( 89 ) + a*JVS( 216 )
  a = -W( 47 ) / JVS(          246  )
  W( 47 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 247 )
  W( 84 ) = W( 84 ) + a*JVS( 248 )
  W( 89 ) = W( 89 ) + a*JVS( 249 )
  a = -W( 48 ) / JVS(          250  )
  W( 48 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 251 )
  W( 84 ) = W( 84 ) + a*JVS( 252 )
  W( 89 ) = W( 89 ) + a*JVS( 253 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 255 )
  W( 84 ) = W( 84 ) + a*JVS( 256 )
  W( 89 ) = W( 89 ) + a*JVS( 257 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  a = -W( 52 ) / JVS(          275  )
  W( 52 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 276 )
  W( 79 ) = W( 79 ) + a*JVS( 277 )
  W( 84 ) = W( 84 ) + a*JVS( 278 )
  W( 85 ) = W( 85 ) + a*JVS( 279 )
  W( 87 ) = W( 87 ) + a*JVS( 280 )
  W( 88 ) = W( 88 ) + a*JVS( 281 )
  W( 89 ) = W( 89 ) + a*JVS( 282 )
  a = -W( 53 ) / JVS(          284  )
  W( 53 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 285 )
  W( 67 ) = W( 67 ) + a*JVS( 286 )
  W( 68 ) = W( 68 ) + a*JVS( 287 )
  W( 71 ) = W( 71 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  W( 80 ) = W( 80 ) + a*JVS( 290 )
  W( 82 ) = W( 82 ) + a*JVS( 291 )
  W( 83 ) = W( 83 ) + a*JVS( 292 )
  W( 84 ) = W( 84 ) + a*JVS( 293 )
  W( 87 ) = W( 87 ) + a*JVS( 294 )
  W( 89 ) = W( 89 ) + a*JVS( 295 )
  a = -W( 54 ) / JVS(          296  )
  W( 54 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 297 )
  W( 79 ) = W( 79 ) + a*JVS( 298 )
  W( 84 ) = W( 84 ) + a*JVS( 299 )
  W( 88 ) = W( 88 ) + a*JVS( 300 )
  W( 89 ) = W( 89 ) + a*JVS( 301 )
  a = -W( 55 ) / JVS(          303  )
  W( 55 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 304 )
  W( 82 ) = W( 82 ) + a*JVS( 305 )
  W( 83 ) = W( 83 ) + a*JVS( 306 )
  W( 84 ) = W( 84 ) + a*JVS( 307 )
  W( 85 ) = W( 85 ) + a*JVS( 308 )
  W( 87 ) = W( 87 ) + a*JVS( 309 )
  W( 89 ) = W( 89 ) + a*JVS( 310 )
  a = -W( 56 ) / JVS(          311  )
  W( 56 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 312 )
  W( 82 ) = W( 82 ) + a*JVS( 313 )
  W( 83 ) = W( 83 ) + a*JVS( 314 )
  W( 84 ) = W( 84 ) + a*JVS( 315 )
  W( 85 ) = W( 85 ) + a*JVS( 316 )
  W( 87 ) = W( 87 ) + a*JVS( 317 )
  W( 89 ) = W( 89 ) + a*JVS( 318 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 325 )
  W( 63 ) = W( 63 ) + a*JVS( 326 )
  W( 79 ) = W( 79 ) + a*JVS( 327 )
  W( 82 ) = W( 82 ) + a*JVS( 328 )
  W( 83 ) = W( 83 ) + a*JVS( 329 )
  W( 84 ) = W( 84 ) + a*JVS( 330 )
  W( 85 ) = W( 85 ) + a*JVS( 331 )
  W( 87 ) = W( 87 ) + a*JVS( 332 )
  W( 89 ) = W( 89 ) + a*JVS( 333 )
  a = -W( 58 ) / JVS(          335  )
  W( 58 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 336 )
  W( 82 ) = W( 82 ) + a*JVS( 337 )
  W( 83 ) = W( 83 ) + a*JVS( 338 )
  W( 84 ) = W( 84 ) + a*JVS( 339 )
  W( 85 ) = W( 85 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 89 ) = W( 89 ) + a*JVS( 342 )
  a = -W( 60 ) / JVS(          372  )
  W( 60 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 82 ) = W( 82 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 85 ) = W( 85 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  a = -W( 61 ) / JVS(          382  )
  W( 61 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 82 ) = W( 82 ) + a*JVS( 384 )
  W( 83 ) = W( 83 ) + a*JVS( 385 )
  W( 84 ) = W( 84 ) + a*JVS( 386 )
  W( 85 ) = W( 85 ) + a*JVS( 387 )
  W( 87 ) = W( 87 ) + a*JVS( 388 )
  W( 88 ) = W( 88 ) + a*JVS( 389 )
  W( 89 ) = W( 89 ) + a*JVS( 390 )
  a = -W( 63 ) / JVS(          417  )
  W( 63 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 418 )
  W( 82 ) = W( 82 ) + a*JVS( 419 )
  W( 83 ) = W( 83 ) + a*JVS( 420 )
  W( 84 ) = W( 84 ) + a*JVS( 421 )
  W( 85 ) = W( 85 ) + a*JVS( 422 )
  W( 87 ) = W( 87 ) + a*JVS( 423 )
  W( 89 ) = W( 89 ) + a*JVS( 424 )
  a = -W( 64 ) / JVS(          426  )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 85 ) = W( 85 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 89 ) = W( 89 ) + a*JVS( 432 )
  a = -W( 65 ) / JVS(          434  )
  W( 65 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 435 )
  W( 83 ) = W( 83 ) + a*JVS( 436 )
  W( 84 ) = W( 84 ) + a*JVS( 437 )
  W( 85 ) = W( 85 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 89 ) = W( 89 ) + a*JVS( 440 )
  JVS( 441) = W( 22 )
  JVS( 442) = W( 24 )
  JVS( 443) = W( 27 )
  JVS( 444) = W( 32 )
  JVS( 445) = W( 38 )
  JVS( 446) = W( 40 )
  JVS( 447) = W( 41 )
  JVS( 448) = W( 42 )
  JVS( 449) = W( 44 )
  JVS( 450) = W( 47 )
  JVS( 451) = W( 48 )
  JVS( 452) = W( 49 )
  JVS( 453) = W( 51 )
  JVS( 454) = W( 52 )
  JVS( 455) = W( 53 )
  JVS( 456) = W( 54 )
  JVS( 457) = W( 55 )
  JVS( 458) = W( 56 )
  JVS( 459) = W( 57 )
  JVS( 460) = W( 58 )
  JVS( 461) = W( 60 )
  JVS( 462) = W( 61 )
  JVS( 463) = W( 63 )
  JVS( 464) = W( 64 )
  JVS( 465) = W( 65 )
  JVS( 466) = W( 66 )
  JVS( 467) = W( 67 )
  JVS( 468) = W( 68 )
  JVS( 469) = W( 69 )
  JVS( 470) = W( 71 )
  JVS( 471) = W( 72 )
  JVS( 472) = W( 73 )
  JVS( 473) = W( 75 )
  JVS( 474) = W( 76 )
  JVS( 475) = W( 77 )
  JVS( 476) = W( 78 )
  JVS( 477) = W( 79 )
  JVS( 478) = W( 80 )
  JVS( 479) = W( 81 )
  JVS( 480) = W( 82 )
  JVS( 481) = W( 83 )
  JVS( 482) = W( 84 )
  JVS( 483) = W( 85 )
  JVS( 484) = W( 87 )
  JVS( 485) = W( 88 )
  JVS( 486) = W( 89 )
  JVS( 487) = W( 90 )
  IF ( ABS(  JVS( 489 )) < TINY(a) ) THEN
         IER = 67                                      
         RETURN
  END IF
   W( 36 ) = JVS( 488 )
   W( 67 ) = JVS( 489 )
   W( 79 ) = JVS( 490 )
   W( 82 ) = JVS( 491 )
   W( 83 ) = JVS( 492 )
   W( 84 ) = JVS( 493 )
   W( 85 ) = JVS( 494 )
   W( 87 ) = JVS( 495 )
   W( 88 ) = JVS( 496 )
   W( 89 ) = JVS( 497 )
  a = -W( 36 ) / JVS(          165  )
  W( 36 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 166 )
  W( 84 ) = W( 84 ) + a*JVS( 167 )
  W( 88 ) = W( 88 ) + a*JVS( 168 )
  JVS( 488) = W( 36 )
  JVS( 489) = W( 67 )
  JVS( 490) = W( 79 )
  JVS( 491) = W( 82 )
  JVS( 492) = W( 83 )
  JVS( 493) = W( 84 )
  JVS( 494) = W( 85 )
  JVS( 495) = W( 87 )
  JVS( 496) = W( 88 )
  JVS( 497) = W( 89 )
  IF ( ABS(  JVS( 499 )) < TINY(a) ) THEN
         IER = 68                                      
         RETURN
  END IF
   W( 37 ) = JVS( 498 )
   W( 68 ) = JVS( 499 )
   W( 79 ) = JVS( 500 )
   W( 82 ) = JVS( 501 )
   W( 83 ) = JVS( 502 )
   W( 84 ) = JVS( 503 )
   W( 85 ) = JVS( 504 )
   W( 87 ) = JVS( 505 )
   W( 88 ) = JVS( 506 )
   W( 89 ) = JVS( 507 )
  a = -W( 37 ) / JVS(          170  )
  W( 37 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 171 )
  W( 84 ) = W( 84 ) + a*JVS( 172 )
  W( 88 ) = W( 88 ) + a*JVS( 173 )
  JVS( 498) = W( 37 )
  JVS( 499) = W( 68 )
  JVS( 500) = W( 79 )
  JVS( 501) = W( 82 )
  JVS( 502) = W( 83 )
  JVS( 503) = W( 84 )
  JVS( 504) = W( 85 )
  JVS( 505) = W( 87 )
  JVS( 506) = W( 88 )
  JVS( 507) = W( 89 )
  IF ( ABS(  JVS( 510 )) < TINY(a) ) THEN
         IER = 69                                      
         RETURN
  END IF
   W( 42 ) = JVS( 508 )
   W( 49 ) = JVS( 509 )
   W( 69 ) = JVS( 510 )
   W( 73 ) = JVS( 511 )
   W( 79 ) = JVS( 512 )
   W( 82 ) = JVS( 513 )
   W( 83 ) = JVS( 514 )
   W( 84 ) = JVS( 515 )
   W( 85 ) = JVS( 516 )
   W( 87 ) = JVS( 517 )
   W( 89 ) = JVS( 518 )
  a = -W( 42 ) / JVS(          195  )
  W( 42 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 196 )
  W( 84 ) = W( 84 ) + a*JVS( 197 )
  W( 89 ) = W( 89 ) + a*JVS( 198 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 255 )
  W( 84 ) = W( 84 ) + a*JVS( 256 )
  W( 89 ) = W( 89 ) + a*JVS( 257 )
  JVS( 508) = W( 42 )
  JVS( 509) = W( 49 )
  JVS( 510) = W( 69 )
  JVS( 511) = W( 73 )
  JVS( 512) = W( 79 )
  JVS( 513) = W( 82 )
  JVS( 514) = W( 83 )
  JVS( 515) = W( 84 )
  JVS( 516) = W( 85 )
  JVS( 517) = W( 87 )
  JVS( 518) = W( 89 )
  IF ( ABS(  JVS( 528 )) < TINY(a) ) THEN
         IER = 70                                      
         RETURN
  END IF
   W( 34 ) = JVS( 519 )
   W( 50 ) = JVS( 520 )
   W( 54 ) = JVS( 521 )
   W( 58 ) = JVS( 522 )
   W( 63 ) = JVS( 523 )
   W( 64 ) = JVS( 524 )
   W( 65 ) = JVS( 525 )
   W( 67 ) = JVS( 526 )
   W( 68 ) = JVS( 527 )
   W( 70 ) = JVS( 528 )
   W( 76 ) = JVS( 529 )
   W( 77 ) = JVS( 530 )
   W( 78 ) = JVS( 531 )
   W( 79 ) = JVS( 532 )
   W( 80 ) = JVS( 533 )
   W( 82 ) = JVS( 534 )
   W( 83 ) = JVS( 535 )
   W( 84 ) = JVS( 536 )
   W( 85 ) = JVS( 537 )
   W( 87 ) = JVS( 538 )
   W( 88 ) = JVS( 539 )
   W( 89 ) = JVS( 540 )
  a = -W( 34 ) / JVS(          147  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 148 )
  W( 84 ) = W( 84 ) + a*JVS( 149 )
  W( 85 ) = W( 85 ) + a*JVS( 150 )
  W( 88 ) = W( 88 ) + a*JVS( 151 )
  W( 89 ) = W( 89 ) + a*JVS( 152 )
  a = -W( 50 ) / JVS(          262  )
  W( 50 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 263 )
  W( 84 ) = W( 84 ) + a*JVS( 264 )
  W( 85 ) = W( 85 ) + a*JVS( 265 )
  W( 88 ) = W( 88 ) + a*JVS( 266 )
  W( 89 ) = W( 89 ) + a*JVS( 267 )
  a = -W( 54 ) / JVS(          296  )
  W( 54 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 297 )
  W( 79 ) = W( 79 ) + a*JVS( 298 )
  W( 84 ) = W( 84 ) + a*JVS( 299 )
  W( 88 ) = W( 88 ) + a*JVS( 300 )
  W( 89 ) = W( 89 ) + a*JVS( 301 )
  a = -W( 58 ) / JVS(          335  )
  W( 58 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 336 )
  W( 82 ) = W( 82 ) + a*JVS( 337 )
  W( 83 ) = W( 83 ) + a*JVS( 338 )
  W( 84 ) = W( 84 ) + a*JVS( 339 )
  W( 85 ) = W( 85 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 89 ) = W( 89 ) + a*JVS( 342 )
  a = -W( 63 ) / JVS(          417  )
  W( 63 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 418 )
  W( 82 ) = W( 82 ) + a*JVS( 419 )
  W( 83 ) = W( 83 ) + a*JVS( 420 )
  W( 84 ) = W( 84 ) + a*JVS( 421 )
  W( 85 ) = W( 85 ) + a*JVS( 422 )
  W( 87 ) = W( 87 ) + a*JVS( 423 )
  W( 89 ) = W( 89 ) + a*JVS( 424 )
  a = -W( 64 ) / JVS(          426  )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 85 ) = W( 85 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 89 ) = W( 89 ) + a*JVS( 432 )
  a = -W( 65 ) / JVS(          434  )
  W( 65 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 435 )
  W( 83 ) = W( 83 ) + a*JVS( 436 )
  W( 84 ) = W( 84 ) + a*JVS( 437 )
  W( 85 ) = W( 85 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 89 ) = W( 89 ) + a*JVS( 440 )
  a = -W( 67 ) / JVS(          489  )
  W( 67 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 490 )
  W( 82 ) = W( 82 ) + a*JVS( 491 )
  W( 83 ) = W( 83 ) + a*JVS( 492 )
  W( 84 ) = W( 84 ) + a*JVS( 493 )
  W( 85 ) = W( 85 ) + a*JVS( 494 )
  W( 87 ) = W( 87 ) + a*JVS( 495 )
  W( 88 ) = W( 88 ) + a*JVS( 496 )
  W( 89 ) = W( 89 ) + a*JVS( 497 )
  a = -W( 68 ) / JVS(          499  )
  W( 68 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 500 )
  W( 82 ) = W( 82 ) + a*JVS( 501 )
  W( 83 ) = W( 83 ) + a*JVS( 502 )
  W( 84 ) = W( 84 ) + a*JVS( 503 )
  W( 85 ) = W( 85 ) + a*JVS( 504 )
  W( 87 ) = W( 87 ) + a*JVS( 505 )
  W( 88 ) = W( 88 ) + a*JVS( 506 )
  W( 89 ) = W( 89 ) + a*JVS( 507 )
  JVS( 519) = W( 34 )
  JVS( 520) = W( 50 )
  JVS( 521) = W( 54 )
  JVS( 522) = W( 58 )
  JVS( 523) = W( 63 )
  JVS( 524) = W( 64 )
  JVS( 525) = W( 65 )
  JVS( 526) = W( 67 )
  JVS( 527) = W( 68 )
  JVS( 528) = W( 70 )
  JVS( 529) = W( 76 )
  JVS( 530) = W( 77 )
  JVS( 531) = W( 78 )
  JVS( 532) = W( 79 )
  JVS( 533) = W( 80 )
  JVS( 534) = W( 82 )
  JVS( 535) = W( 83 )
  JVS( 536) = W( 84 )
  JVS( 537) = W( 85 )
  JVS( 538) = W( 87 )
  JVS( 539) = W( 88 )
  JVS( 540) = W( 89 )
  IF ( ABS(  JVS( 543 )) < TINY(a) ) THEN
         IER = 71                                      
         RETURN
  END IF
   W( 67 ) = JVS( 541 )
   W( 68 ) = JVS( 542 )
   W( 71 ) = JVS( 543 )
   W( 79 ) = JVS( 544 )
   W( 82 ) = JVS( 545 )
   W( 83 ) = JVS( 546 )
   W( 84 ) = JVS( 547 )
   W( 85 ) = JVS( 548 )
   W( 87 ) = JVS( 549 )
   W( 88 ) = JVS( 550 )
   W( 89 ) = JVS( 551 )
  a = -W( 67 ) / JVS(          489  )
  W( 67 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 490 )
  W( 82 ) = W( 82 ) + a*JVS( 491 )
  W( 83 ) = W( 83 ) + a*JVS( 492 )
  W( 84 ) = W( 84 ) + a*JVS( 493 )
  W( 85 ) = W( 85 ) + a*JVS( 494 )
  W( 87 ) = W( 87 ) + a*JVS( 495 )
  W( 88 ) = W( 88 ) + a*JVS( 496 )
  W( 89 ) = W( 89 ) + a*JVS( 497 )
  a = -W( 68 ) / JVS(          499  )
  W( 68 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 500 )
  W( 82 ) = W( 82 ) + a*JVS( 501 )
  W( 83 ) = W( 83 ) + a*JVS( 502 )
  W( 84 ) = W( 84 ) + a*JVS( 503 )
  W( 85 ) = W( 85 ) + a*JVS( 504 )
  W( 87 ) = W( 87 ) + a*JVS( 505 )
  W( 88 ) = W( 88 ) + a*JVS( 506 )
  W( 89 ) = W( 89 ) + a*JVS( 507 )
  JVS( 541) = W( 67 )
  JVS( 542) = W( 68 )
  JVS( 543) = W( 71 )
  JVS( 544) = W( 79 )
  JVS( 545) = W( 82 )
  JVS( 546) = W( 83 )
  JVS( 547) = W( 84 )
  JVS( 548) = W( 85 )
  JVS( 549) = W( 87 )
  JVS( 550) = W( 88 )
  JVS( 551) = W( 89 )
  IF ( ABS(  JVS( 562 )) < TINY(a) ) THEN
         IER = 72                                      
         RETURN
  END IF
   W( 16 ) = JVS( 552 )
   W( 17 ) = JVS( 553 )
   W( 40 ) = JVS( 554 )
   W( 44 ) = JVS( 555 )
   W( 50 ) = JVS( 556 )
   W( 51 ) = JVS( 557 )
   W( 54 ) = JVS( 558 )
   W( 64 ) = JVS( 559 )
   W( 65 ) = JVS( 560 )
   W( 71 ) = JVS( 561 )
   W( 72 ) = JVS( 562 )
   W( 78 ) = JVS( 563 )
   W( 79 ) = JVS( 564 )
   W( 80 ) = JVS( 565 )
   W( 81 ) = JVS( 566 )
   W( 82 ) = JVS( 567 )
   W( 83 ) = JVS( 568 )
   W( 84 ) = JVS( 569 )
   W( 85 ) = JVS( 570 )
   W( 86 ) = JVS( 571 )
   W( 87 ) = JVS( 572 )
   W( 88 ) = JVS( 573 )
   W( 89 ) = JVS( 574 )
  a = -W( 16 ) / JVS(           73  )
  W( 16 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 74 )
  a = -W( 17 ) / JVS(           75  )
  W( 17 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 76 )
  a = -W( 40 ) / JVS(          183  )
  W( 40 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 184 )
  W( 79 ) = W( 79 ) + a*JVS( 185 )
  W( 82 ) = W( 82 ) + a*JVS( 186 )
  W( 84 ) = W( 84 ) + a*JVS( 187 )
  W( 85 ) = W( 85 ) + a*JVS( 188 )
  a = -W( 44 ) / JVS(          210  )
  W( 44 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 211 )
  W( 79 ) = W( 79 ) + a*JVS( 212 )
  W( 82 ) = W( 82 ) + a*JVS( 213 )
  W( 84 ) = W( 84 ) + a*JVS( 214 )
  W( 88 ) = W( 88 ) + a*JVS( 215 )
  W( 89 ) = W( 89 ) + a*JVS( 216 )
  a = -W( 50 ) / JVS(          262  )
  W( 50 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 263 )
  W( 84 ) = W( 84 ) + a*JVS( 264 )
  W( 85 ) = W( 85 ) + a*JVS( 265 )
  W( 88 ) = W( 88 ) + a*JVS( 266 )
  W( 89 ) = W( 89 ) + a*JVS( 267 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  a = -W( 54 ) / JVS(          296  )
  W( 54 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 297 )
  W( 79 ) = W( 79 ) + a*JVS( 298 )
  W( 84 ) = W( 84 ) + a*JVS( 299 )
  W( 88 ) = W( 88 ) + a*JVS( 300 )
  W( 89 ) = W( 89 ) + a*JVS( 301 )
  a = -W( 64 ) / JVS(          426  )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 85 ) = W( 85 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 89 ) = W( 89 ) + a*JVS( 432 )
  a = -W( 65 ) / JVS(          434  )
  W( 65 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 435 )
  W( 83 ) = W( 83 ) + a*JVS( 436 )
  W( 84 ) = W( 84 ) + a*JVS( 437 )
  W( 85 ) = W( 85 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 89 ) = W( 89 ) + a*JVS( 440 )
  a = -W( 71 ) / JVS(          543  )
  W( 71 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 544 )
  W( 82 ) = W( 82 ) + a*JVS( 545 )
  W( 83 ) = W( 83 ) + a*JVS( 546 )
  W( 84 ) = W( 84 ) + a*JVS( 547 )
  W( 85 ) = W( 85 ) + a*JVS( 548 )
  W( 87 ) = W( 87 ) + a*JVS( 549 )
  W( 88 ) = W( 88 ) + a*JVS( 550 )
  W( 89 ) = W( 89 ) + a*JVS( 551 )
  JVS( 552) = W( 16 )
  JVS( 553) = W( 17 )
  JVS( 554) = W( 40 )
  JVS( 555) = W( 44 )
  JVS( 556) = W( 50 )
  JVS( 557) = W( 51 )
  JVS( 558) = W( 54 )
  JVS( 559) = W( 64 )
  JVS( 560) = W( 65 )
  JVS( 561) = W( 71 )
  JVS( 562) = W( 72 )
  JVS( 563) = W( 78 )
  JVS( 564) = W( 79 )
  JVS( 565) = W( 80 )
  JVS( 566) = W( 81 )
  JVS( 567) = W( 82 )
  JVS( 568) = W( 83 )
  JVS( 569) = W( 84 )
  JVS( 570) = W( 85 )
  JVS( 571) = W( 86 )
  JVS( 572) = W( 87 )
  JVS( 573) = W( 88 )
  JVS( 574) = W( 89 )
  IF ( ABS(  JVS( 577 )) < TINY(a) ) THEN
         IER = 73                                      
         RETURN
  END IF
   W( 60 ) = JVS( 575 )
   W( 63 ) = JVS( 576 )
   W( 73 ) = JVS( 577 )
   W( 79 ) = JVS( 578 )
   W( 82 ) = JVS( 579 )
   W( 83 ) = JVS( 580 )
   W( 84 ) = JVS( 581 )
   W( 85 ) = JVS( 582 )
   W( 87 ) = JVS( 583 )
   W( 89 ) = JVS( 584 )
  a = -W( 60 ) / JVS(          372  )
  W( 60 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 82 ) = W( 82 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 85 ) = W( 85 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  a = -W( 63 ) / JVS(          417  )
  W( 63 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 418 )
  W( 82 ) = W( 82 ) + a*JVS( 419 )
  W( 83 ) = W( 83 ) + a*JVS( 420 )
  W( 84 ) = W( 84 ) + a*JVS( 421 )
  W( 85 ) = W( 85 ) + a*JVS( 422 )
  W( 87 ) = W( 87 ) + a*JVS( 423 )
  W( 89 ) = W( 89 ) + a*JVS( 424 )
  JVS( 575) = W( 60 )
  JVS( 576) = W( 63 )
  JVS( 577) = W( 73 )
  JVS( 578) = W( 79 )
  JVS( 579) = W( 82 )
  JVS( 580) = W( 83 )
  JVS( 581) = W( 84 )
  JVS( 582) = W( 85 )
  JVS( 583) = W( 87 )
  JVS( 584) = W( 89 )
  IF ( ABS(  JVS( 600 )) < TINY(a) ) THEN
         IER = 74                                      
         RETURN
  END IF
   W( 20 ) = JVS( 585 )
   W( 22 ) = JVS( 586 )
   W( 23 ) = JVS( 587 )
   W( 42 ) = JVS( 588 )
   W( 45 ) = JVS( 589 )
   W( 49 ) = JVS( 590 )
   W( 55 ) = JVS( 591 )
   W( 56 ) = JVS( 592 )
   W( 58 ) = JVS( 593 )
   W( 64 ) = JVS( 594 )
   W( 65 ) = JVS( 595 )
   W( 69 ) = JVS( 596 )
   W( 70 ) = JVS( 597 )
   W( 71 ) = JVS( 598 )
   W( 73 ) = JVS( 599 )
   W( 74 ) = JVS( 600 )
   W( 75 ) = JVS( 601 )
   W( 76 ) = JVS( 602 )
   W( 77 ) = JVS( 603 )
   W( 78 ) = JVS( 604 )
   W( 79 ) = JVS( 605 )
   W( 80 ) = JVS( 606 )
   W( 81 ) = JVS( 607 )
   W( 82 ) = JVS( 608 )
   W( 83 ) = JVS( 609 )
   W( 84 ) = JVS( 610 )
   W( 85 ) = JVS( 611 )
   W( 86 ) = JVS( 612 )
   W( 87 ) = JVS( 613 )
   W( 88 ) = JVS( 614 )
   W( 89 ) = JVS( 615 )
   W( 90 ) = JVS( 616 )
  a = -W( 20 ) / JVS(           82  )
  W( 20 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 83 )
  a = -W( 22 ) / JVS(           88  )
  W( 22 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 89 )
  a = -W( 23 ) / JVS(           90  )
  W( 23 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 91 )
  W( 84 ) = W( 84 ) + a*JVS( 92 )
  a = -W( 42 ) / JVS(          195  )
  W( 42 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 196 )
  W( 84 ) = W( 84 ) + a*JVS( 197 )
  W( 89 ) = W( 89 ) + a*JVS( 198 )
  a = -W( 45 ) / JVS(          217  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 218 )
  W( 84 ) = W( 84 ) + a*JVS( 219 )
  W( 89 ) = W( 89 ) + a*JVS( 220 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 255 )
  W( 84 ) = W( 84 ) + a*JVS( 256 )
  W( 89 ) = W( 89 ) + a*JVS( 257 )
  a = -W( 55 ) / JVS(          303  )
  W( 55 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 304 )
  W( 82 ) = W( 82 ) + a*JVS( 305 )
  W( 83 ) = W( 83 ) + a*JVS( 306 )
  W( 84 ) = W( 84 ) + a*JVS( 307 )
  W( 85 ) = W( 85 ) + a*JVS( 308 )
  W( 87 ) = W( 87 ) + a*JVS( 309 )
  W( 89 ) = W( 89 ) + a*JVS( 310 )
  a = -W( 56 ) / JVS(          311  )
  W( 56 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 312 )
  W( 82 ) = W( 82 ) + a*JVS( 313 )
  W( 83 ) = W( 83 ) + a*JVS( 314 )
  W( 84 ) = W( 84 ) + a*JVS( 315 )
  W( 85 ) = W( 85 ) + a*JVS( 316 )
  W( 87 ) = W( 87 ) + a*JVS( 317 )
  W( 89 ) = W( 89 ) + a*JVS( 318 )
  a = -W( 58 ) / JVS(          335  )
  W( 58 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 336 )
  W( 82 ) = W( 82 ) + a*JVS( 337 )
  W( 83 ) = W( 83 ) + a*JVS( 338 )
  W( 84 ) = W( 84 ) + a*JVS( 339 )
  W( 85 ) = W( 85 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 89 ) = W( 89 ) + a*JVS( 342 )
  a = -W( 64 ) / JVS(          426  )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 85 ) = W( 85 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 89 ) = W( 89 ) + a*JVS( 432 )
  a = -W( 65 ) / JVS(          434  )
  W( 65 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 435 )
  W( 83 ) = W( 83 ) + a*JVS( 436 )
  W( 84 ) = W( 84 ) + a*JVS( 437 )
  W( 85 ) = W( 85 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 89 ) = W( 89 ) + a*JVS( 440 )
  a = -W( 69 ) / JVS(          510  )
  W( 69 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 511 )
  W( 79 ) = W( 79 ) + a*JVS( 512 )
  W( 82 ) = W( 82 ) + a*JVS( 513 )
  W( 83 ) = W( 83 ) + a*JVS( 514 )
  W( 84 ) = W( 84 ) + a*JVS( 515 )
  W( 85 ) = W( 85 ) + a*JVS( 516 )
  W( 87 ) = W( 87 ) + a*JVS( 517 )
  W( 89 ) = W( 89 ) + a*JVS( 518 )
  a = -W( 70 ) / JVS(          528  )
  W( 70 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 529 )
  W( 77 ) = W( 77 ) + a*JVS( 530 )
  W( 78 ) = W( 78 ) + a*JVS( 531 )
  W( 79 ) = W( 79 ) + a*JVS( 532 )
  W( 80 ) = W( 80 ) + a*JVS( 533 )
  W( 82 ) = W( 82 ) + a*JVS( 534 )
  W( 83 ) = W( 83 ) + a*JVS( 535 )
  W( 84 ) = W( 84 ) + a*JVS( 536 )
  W( 85 ) = W( 85 ) + a*JVS( 537 )
  W( 87 ) = W( 87 ) + a*JVS( 538 )
  W( 88 ) = W( 88 ) + a*JVS( 539 )
  W( 89 ) = W( 89 ) + a*JVS( 540 )
  a = -W( 71 ) / JVS(          543  )
  W( 71 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 544 )
  W( 82 ) = W( 82 ) + a*JVS( 545 )
  W( 83 ) = W( 83 ) + a*JVS( 546 )
  W( 84 ) = W( 84 ) + a*JVS( 547 )
  W( 85 ) = W( 85 ) + a*JVS( 548 )
  W( 87 ) = W( 87 ) + a*JVS( 549 )
  W( 88 ) = W( 88 ) + a*JVS( 550 )
  W( 89 ) = W( 89 ) + a*JVS( 551 )
  a = -W( 73 ) / JVS(          577  )
  W( 73 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 578 )
  W( 82 ) = W( 82 ) + a*JVS( 579 )
  W( 83 ) = W( 83 ) + a*JVS( 580 )
  W( 84 ) = W( 84 ) + a*JVS( 581 )
  W( 85 ) = W( 85 ) + a*JVS( 582 )
  W( 87 ) = W( 87 ) + a*JVS( 583 )
  W( 89 ) = W( 89 ) + a*JVS( 584 )
  JVS( 585) = W( 20 )
  JVS( 586) = W( 22 )
  JVS( 587) = W( 23 )
  JVS( 588) = W( 42 )
  JVS( 589) = W( 45 )
  JVS( 590) = W( 49 )
  JVS( 591) = W( 55 )
  JVS( 592) = W( 56 )
  JVS( 593) = W( 58 )
  JVS( 594) = W( 64 )
  JVS( 595) = W( 65 )
  JVS( 596) = W( 69 )
  JVS( 597) = W( 70 )
  JVS( 598) = W( 71 )
  JVS( 599) = W( 73 )
  JVS( 600) = W( 74 )
  JVS( 601) = W( 75 )
  JVS( 602) = W( 76 )
  JVS( 603) = W( 77 )
  JVS( 604) = W( 78 )
  JVS( 605) = W( 79 )
  JVS( 606) = W( 80 )
  JVS( 607) = W( 81 )
  JVS( 608) = W( 82 )
  JVS( 609) = W( 83 )
  JVS( 610) = W( 84 )
  JVS( 611) = W( 85 )
  JVS( 612) = W( 86 )
  JVS( 613) = W( 87 )
  JVS( 614) = W( 88 )
  JVS( 615) = W( 89 )
  JVS( 616) = W( 90 )
  IF ( ABS(  JVS( 620 )) < TINY(a) ) THEN
         IER = 75                                      
         RETURN
  END IF
   W( 47 ) = JVS( 617 )
   W( 51 ) = JVS( 618 )
   W( 60 ) = JVS( 619 )
   W( 75 ) = JVS( 620 )
   W( 79 ) = JVS( 621 )
   W( 82 ) = JVS( 622 )
   W( 83 ) = JVS( 623 )
   W( 84 ) = JVS( 624 )
   W( 85 ) = JVS( 625 )
   W( 87 ) = JVS( 626 )
   W( 89 ) = JVS( 627 )
  a = -W( 47 ) / JVS(          246  )
  W( 47 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 247 )
  W( 84 ) = W( 84 ) + a*JVS( 248 )
  W( 89 ) = W( 89 ) + a*JVS( 249 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  a = -W( 60 ) / JVS(          372  )
  W( 60 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 82 ) = W( 82 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 85 ) = W( 85 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  JVS( 617) = W( 47 )
  JVS( 618) = W( 51 )
  JVS( 619) = W( 60 )
  JVS( 620) = W( 75 )
  JVS( 621) = W( 79 )
  JVS( 622) = W( 82 )
  JVS( 623) = W( 83 )
  JVS( 624) = W( 84 )
  JVS( 625) = W( 85 )
  JVS( 626) = W( 87 )
  JVS( 627) = W( 89 )
  IF ( ABS(  JVS( 636 )) < TINY(a) ) THEN
         IER = 76                                      
         RETURN
  END IF
   W( 38 ) = JVS( 628 )
   W( 42 ) = JVS( 629 )
   W( 45 ) = JVS( 630 )
   W( 47 ) = JVS( 631 )
   W( 49 ) = JVS( 632 )
   W( 51 ) = JVS( 633 )
   W( 73 ) = JVS( 634 )
   W( 75 ) = JVS( 635 )
   W( 76 ) = JVS( 636 )
   W( 77 ) = JVS( 637 )
   W( 79 ) = JVS( 638 )
   W( 82 ) = JVS( 639 )
   W( 83 ) = JVS( 640 )
   W( 84 ) = JVS( 641 )
   W( 85 ) = JVS( 642 )
   W( 87 ) = JVS( 643 )
   W( 89 ) = JVS( 644 )
  a = -W( 38 ) / JVS(          174  )
  W( 38 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 175 )
  W( 84 ) = W( 84 ) + a*JVS( 176 )
  W( 89 ) = W( 89 ) + a*JVS( 177 )
  a = -W( 42 ) / JVS(          195  )
  W( 42 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 196 )
  W( 84 ) = W( 84 ) + a*JVS( 197 )
  W( 89 ) = W( 89 ) + a*JVS( 198 )
  a = -W( 45 ) / JVS(          217  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 218 )
  W( 84 ) = W( 84 ) + a*JVS( 219 )
  W( 89 ) = W( 89 ) + a*JVS( 220 )
  a = -W( 47 ) / JVS(          246  )
  W( 47 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 247 )
  W( 84 ) = W( 84 ) + a*JVS( 248 )
  W( 89 ) = W( 89 ) + a*JVS( 249 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 255 )
  W( 84 ) = W( 84 ) + a*JVS( 256 )
  W( 89 ) = W( 89 ) + a*JVS( 257 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  a = -W( 73 ) / JVS(          577  )
  W( 73 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 578 )
  W( 82 ) = W( 82 ) + a*JVS( 579 )
  W( 83 ) = W( 83 ) + a*JVS( 580 )
  W( 84 ) = W( 84 ) + a*JVS( 581 )
  W( 85 ) = W( 85 ) + a*JVS( 582 )
  W( 87 ) = W( 87 ) + a*JVS( 583 )
  W( 89 ) = W( 89 ) + a*JVS( 584 )
  a = -W( 75 ) / JVS(          620  )
  W( 75 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 621 )
  W( 82 ) = W( 82 ) + a*JVS( 622 )
  W( 83 ) = W( 83 ) + a*JVS( 623 )
  W( 84 ) = W( 84 ) + a*JVS( 624 )
  W( 85 ) = W( 85 ) + a*JVS( 625 )
  W( 87 ) = W( 87 ) + a*JVS( 626 )
  W( 89 ) = W( 89 ) + a*JVS( 627 )
  JVS( 628) = W( 38 )
  JVS( 629) = W( 42 )
  JVS( 630) = W( 45 )
  JVS( 631) = W( 47 )
  JVS( 632) = W( 49 )
  JVS( 633) = W( 51 )
  JVS( 634) = W( 73 )
  JVS( 635) = W( 75 )
  JVS( 636) = W( 76 )
  JVS( 637) = W( 77 )
  JVS( 638) = W( 79 )
  JVS( 639) = W( 82 )
  JVS( 640) = W( 83 )
  JVS( 641) = W( 84 )
  JVS( 642) = W( 85 )
  JVS( 643) = W( 87 )
  JVS( 644) = W( 89 )
  IF ( ABS(  JVS( 654 )) < TINY(a) ) THEN
         IER = 77                                      
         RETURN
  END IF
   W( 38 ) = JVS( 645 )
   W( 42 ) = JVS( 646 )
   W( 45 ) = JVS( 647 )
   W( 47 ) = JVS( 648 )
   W( 49 ) = JVS( 649 )
   W( 51 ) = JVS( 650 )
   W( 73 ) = JVS( 651 )
   W( 75 ) = JVS( 652 )
   W( 76 ) = JVS( 653 )
   W( 77 ) = JVS( 654 )
   W( 79 ) = JVS( 655 )
   W( 82 ) = JVS( 656 )
   W( 83 ) = JVS( 657 )
   W( 84 ) = JVS( 658 )
   W( 85 ) = JVS( 659 )
   W( 87 ) = JVS( 660 )
   W( 89 ) = JVS( 661 )
  a = -W( 38 ) / JVS(          174  )
  W( 38 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 175 )
  W( 84 ) = W( 84 ) + a*JVS( 176 )
  W( 89 ) = W( 89 ) + a*JVS( 177 )
  a = -W( 42 ) / JVS(          195  )
  W( 42 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 196 )
  W( 84 ) = W( 84 ) + a*JVS( 197 )
  W( 89 ) = W( 89 ) + a*JVS( 198 )
  a = -W( 45 ) / JVS(          217  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 218 )
  W( 84 ) = W( 84 ) + a*JVS( 219 )
  W( 89 ) = W( 89 ) + a*JVS( 220 )
  a = -W( 47 ) / JVS(          246  )
  W( 47 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 247 )
  W( 84 ) = W( 84 ) + a*JVS( 248 )
  W( 89 ) = W( 89 ) + a*JVS( 249 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 255 )
  W( 84 ) = W( 84 ) + a*JVS( 256 )
  W( 89 ) = W( 89 ) + a*JVS( 257 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  a = -W( 73 ) / JVS(          577  )
  W( 73 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 578 )
  W( 82 ) = W( 82 ) + a*JVS( 579 )
  W( 83 ) = W( 83 ) + a*JVS( 580 )
  W( 84 ) = W( 84 ) + a*JVS( 581 )
  W( 85 ) = W( 85 ) + a*JVS( 582 )
  W( 87 ) = W( 87 ) + a*JVS( 583 )
  W( 89 ) = W( 89 ) + a*JVS( 584 )
  a = -W( 75 ) / JVS(          620  )
  W( 75 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 621 )
  W( 82 ) = W( 82 ) + a*JVS( 622 )
  W( 83 ) = W( 83 ) + a*JVS( 623 )
  W( 84 ) = W( 84 ) + a*JVS( 624 )
  W( 85 ) = W( 85 ) + a*JVS( 625 )
  W( 87 ) = W( 87 ) + a*JVS( 626 )
  W( 89 ) = W( 89 ) + a*JVS( 627 )
  a = -W( 76 ) / JVS(          636  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 637 )
  W( 79 ) = W( 79 ) + a*JVS( 638 )
  W( 82 ) = W( 82 ) + a*JVS( 639 )
  W( 83 ) = W( 83 ) + a*JVS( 640 )
  W( 84 ) = W( 84 ) + a*JVS( 641 )
  W( 85 ) = W( 85 ) + a*JVS( 642 )
  W( 87 ) = W( 87 ) + a*JVS( 643 )
  W( 89 ) = W( 89 ) + a*JVS( 644 )
  JVS( 645) = W( 38 )
  JVS( 646) = W( 42 )
  JVS( 647) = W( 45 )
  JVS( 648) = W( 47 )
  JVS( 649) = W( 49 )
  JVS( 650) = W( 51 )
  JVS( 651) = W( 73 )
  JVS( 652) = W( 75 )
  JVS( 653) = W( 76 )
  JVS( 654) = W( 77 )
  JVS( 655) = W( 79 )
  JVS( 656) = W( 82 )
  JVS( 657) = W( 83 )
  JVS( 658) = W( 84 )
  JVS( 659) = W( 85 )
  JVS( 660) = W( 87 )
  JVS( 661) = W( 89 )
  IF ( ABS(  JVS( 664 )) < TINY(a) ) THEN
         IER = 78                                      
         RETURN
  END IF
   W( 54 ) = JVS( 662 )
   W( 71 ) = JVS( 663 )
   W( 78 ) = JVS( 664 )
   W( 79 ) = JVS( 665 )
   W( 82 ) = JVS( 666 )
   W( 83 ) = JVS( 667 )
   W( 84 ) = JVS( 668 )
   W( 85 ) = JVS( 669 )
   W( 87 ) = JVS( 670 )
   W( 88 ) = JVS( 671 )
   W( 89 ) = JVS( 672 )
  a = -W( 54 ) / JVS(          296  )
  W( 54 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 297 )
  W( 79 ) = W( 79 ) + a*JVS( 298 )
  W( 84 ) = W( 84 ) + a*JVS( 299 )
  W( 88 ) = W( 88 ) + a*JVS( 300 )
  W( 89 ) = W( 89 ) + a*JVS( 301 )
  a = -W( 71 ) / JVS(          543  )
  W( 71 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 544 )
  W( 82 ) = W( 82 ) + a*JVS( 545 )
  W( 83 ) = W( 83 ) + a*JVS( 546 )
  W( 84 ) = W( 84 ) + a*JVS( 547 )
  W( 85 ) = W( 85 ) + a*JVS( 548 )
  W( 87 ) = W( 87 ) + a*JVS( 549 )
  W( 88 ) = W( 88 ) + a*JVS( 550 )
  W( 89 ) = W( 89 ) + a*JVS( 551 )
  JVS( 662) = W( 54 )
  JVS( 663) = W( 71 )
  JVS( 664) = W( 78 )
  JVS( 665) = W( 79 )
  JVS( 666) = W( 82 )
  JVS( 667) = W( 83 )
  JVS( 668) = W( 84 )
  JVS( 669) = W( 85 )
  JVS( 670) = W( 87 )
  JVS( 671) = W( 88 )
  JVS( 672) = W( 89 )
  IF ( ABS(  JVS( 693 )) < TINY(a) ) THEN
         IER = 79                                      
         RETURN
  END IF
   W( 33 ) = JVS( 673 )
   W( 36 ) = JVS( 674 )
   W( 37 ) = JVS( 675 )
   W( 38 ) = JVS( 676 )
   W( 39 ) = JVS( 677 )
   W( 42 ) = JVS( 678 )
   W( 45 ) = JVS( 679 )
   W( 47 ) = JVS( 680 )
   W( 48 ) = JVS( 681 )
   W( 49 ) = JVS( 682 )
   W( 50 ) = JVS( 683 )
   W( 51 ) = JVS( 684 )
   W( 54 ) = JVS( 685 )
   W( 57 ) = JVS( 686 )
   W( 60 ) = JVS( 687 )
   W( 63 ) = JVS( 688 )
   W( 71 ) = JVS( 689 )
   W( 73 ) = JVS( 690 )
   W( 75 ) = JVS( 691 )
   W( 78 ) = JVS( 692 )
   W( 79 ) = JVS( 693 )
   W( 82 ) = JVS( 694 )
   W( 83 ) = JVS( 695 )
   W( 84 ) = JVS( 696 )
   W( 85 ) = JVS( 697 )
   W( 87 ) = JVS( 698 )
   W( 88 ) = JVS( 699 )
   W( 89 ) = JVS( 700 )
  a = -W( 33 ) / JVS(          141  )
  W( 33 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 142 )
  W( 79 ) = W( 79 ) + a*JVS( 143 )
  W( 87 ) = W( 87 ) + a*JVS( 144 )
  W( 88 ) = W( 88 ) + a*JVS( 145 )
  W( 89 ) = W( 89 ) + a*JVS( 146 )
  a = -W( 36 ) / JVS(          165  )
  W( 36 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 166 )
  W( 84 ) = W( 84 ) + a*JVS( 167 )
  W( 88 ) = W( 88 ) + a*JVS( 168 )
  a = -W( 37 ) / JVS(          170  )
  W( 37 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 171 )
  W( 84 ) = W( 84 ) + a*JVS( 172 )
  W( 88 ) = W( 88 ) + a*JVS( 173 )
  a = -W( 38 ) / JVS(          174  )
  W( 38 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 175 )
  W( 84 ) = W( 84 ) + a*JVS( 176 )
  W( 89 ) = W( 89 ) + a*JVS( 177 )
  a = -W( 39 ) / JVS(          178  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 179 )
  W( 79 ) = W( 79 ) + a*JVS( 180 )
  W( 84 ) = W( 84 ) + a*JVS( 181 )
  W( 88 ) = W( 88 ) + a*JVS( 182 )
  a = -W( 42 ) / JVS(          195  )
  W( 42 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 196 )
  W( 84 ) = W( 84 ) + a*JVS( 197 )
  W( 89 ) = W( 89 ) + a*JVS( 198 )
  a = -W( 45 ) / JVS(          217  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 218 )
  W( 84 ) = W( 84 ) + a*JVS( 219 )
  W( 89 ) = W( 89 ) + a*JVS( 220 )
  a = -W( 47 ) / JVS(          246  )
  W( 47 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 247 )
  W( 84 ) = W( 84 ) + a*JVS( 248 )
  W( 89 ) = W( 89 ) + a*JVS( 249 )
  a = -W( 48 ) / JVS(          250  )
  W( 48 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 251 )
  W( 84 ) = W( 84 ) + a*JVS( 252 )
  W( 89 ) = W( 89 ) + a*JVS( 253 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 255 )
  W( 84 ) = W( 84 ) + a*JVS( 256 )
  W( 89 ) = W( 89 ) + a*JVS( 257 )
  a = -W( 50 ) / JVS(          262  )
  W( 50 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 263 )
  W( 84 ) = W( 84 ) + a*JVS( 264 )
  W( 85 ) = W( 85 ) + a*JVS( 265 )
  W( 88 ) = W( 88 ) + a*JVS( 266 )
  W( 89 ) = W( 89 ) + a*JVS( 267 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  a = -W( 54 ) / JVS(          296  )
  W( 54 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 297 )
  W( 79 ) = W( 79 ) + a*JVS( 298 )
  W( 84 ) = W( 84 ) + a*JVS( 299 )
  W( 88 ) = W( 88 ) + a*JVS( 300 )
  W( 89 ) = W( 89 ) + a*JVS( 301 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 325 )
  W( 63 ) = W( 63 ) + a*JVS( 326 )
  W( 79 ) = W( 79 ) + a*JVS( 327 )
  W( 82 ) = W( 82 ) + a*JVS( 328 )
  W( 83 ) = W( 83 ) + a*JVS( 329 )
  W( 84 ) = W( 84 ) + a*JVS( 330 )
  W( 85 ) = W( 85 ) + a*JVS( 331 )
  W( 87 ) = W( 87 ) + a*JVS( 332 )
  W( 89 ) = W( 89 ) + a*JVS( 333 )
  a = -W( 60 ) / JVS(          372  )
  W( 60 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 82 ) = W( 82 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 85 ) = W( 85 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  a = -W( 63 ) / JVS(          417  )
  W( 63 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 418 )
  W( 82 ) = W( 82 ) + a*JVS( 419 )
  W( 83 ) = W( 83 ) + a*JVS( 420 )
  W( 84 ) = W( 84 ) + a*JVS( 421 )
  W( 85 ) = W( 85 ) + a*JVS( 422 )
  W( 87 ) = W( 87 ) + a*JVS( 423 )
  W( 89 ) = W( 89 ) + a*JVS( 424 )
  a = -W( 71 ) / JVS(          543  )
  W( 71 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 544 )
  W( 82 ) = W( 82 ) + a*JVS( 545 )
  W( 83 ) = W( 83 ) + a*JVS( 546 )
  W( 84 ) = W( 84 ) + a*JVS( 547 )
  W( 85 ) = W( 85 ) + a*JVS( 548 )
  W( 87 ) = W( 87 ) + a*JVS( 549 )
  W( 88 ) = W( 88 ) + a*JVS( 550 )
  W( 89 ) = W( 89 ) + a*JVS( 551 )
  a = -W( 73 ) / JVS(          577  )
  W( 73 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 578 )
  W( 82 ) = W( 82 ) + a*JVS( 579 )
  W( 83 ) = W( 83 ) + a*JVS( 580 )
  W( 84 ) = W( 84 ) + a*JVS( 581 )
  W( 85 ) = W( 85 ) + a*JVS( 582 )
  W( 87 ) = W( 87 ) + a*JVS( 583 )
  W( 89 ) = W( 89 ) + a*JVS( 584 )
  a = -W( 75 ) / JVS(          620  )
  W( 75 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 621 )
  W( 82 ) = W( 82 ) + a*JVS( 622 )
  W( 83 ) = W( 83 ) + a*JVS( 623 )
  W( 84 ) = W( 84 ) + a*JVS( 624 )
  W( 85 ) = W( 85 ) + a*JVS( 625 )
  W( 87 ) = W( 87 ) + a*JVS( 626 )
  W( 89 ) = W( 89 ) + a*JVS( 627 )
  a = -W( 78 ) / JVS(          664  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 665 )
  W( 82 ) = W( 82 ) + a*JVS( 666 )
  W( 83 ) = W( 83 ) + a*JVS( 667 )
  W( 84 ) = W( 84 ) + a*JVS( 668 )
  W( 85 ) = W( 85 ) + a*JVS( 669 )
  W( 87 ) = W( 87 ) + a*JVS( 670 )
  W( 88 ) = W( 88 ) + a*JVS( 671 )
  W( 89 ) = W( 89 ) + a*JVS( 672 )
  JVS( 673) = W( 33 )
  JVS( 674) = W( 36 )
  JVS( 675) = W( 37 )
  JVS( 676) = W( 38 )
  JVS( 677) = W( 39 )
  JVS( 678) = W( 42 )
  JVS( 679) = W( 45 )
  JVS( 680) = W( 47 )
  JVS( 681) = W( 48 )
  JVS( 682) = W( 49 )
  JVS( 683) = W( 50 )
  JVS( 684) = W( 51 )
  JVS( 685) = W( 54 )
  JVS( 686) = W( 57 )
  JVS( 687) = W( 60 )
  JVS( 688) = W( 63 )
  JVS( 689) = W( 71 )
  JVS( 690) = W( 73 )
  JVS( 691) = W( 75 )
  JVS( 692) = W( 78 )
  JVS( 693) = W( 79 )
  JVS( 694) = W( 82 )
  JVS( 695) = W( 83 )
  JVS( 696) = W( 84 )
  JVS( 697) = W( 85 )
  JVS( 698) = W( 87 )
  JVS( 699) = W( 88 )
  JVS( 700) = W( 89 )
  IF ( ABS(  JVS( 707 )) < TINY(a) ) THEN
         IER = 80                                      
         RETURN
  END IF
   W( 22 ) = JVS( 701 )
   W( 70 ) = JVS( 702 )
   W( 76 ) = JVS( 703 )
   W( 77 ) = JVS( 704 )
   W( 78 ) = JVS( 705 )
   W( 79 ) = JVS( 706 )
   W( 80 ) = JVS( 707 )
   W( 82 ) = JVS( 708 )
   W( 83 ) = JVS( 709 )
   W( 84 ) = JVS( 710 )
   W( 85 ) = JVS( 711 )
   W( 86 ) = JVS( 712 )
   W( 87 ) = JVS( 713 )
   W( 88 ) = JVS( 714 )
   W( 89 ) = JVS( 715 )
  a = -W( 22 ) / JVS(           88  )
  W( 22 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 89 )
  a = -W( 70 ) / JVS(          528  )
  W( 70 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 529 )
  W( 77 ) = W( 77 ) + a*JVS( 530 )
  W( 78 ) = W( 78 ) + a*JVS( 531 )
  W( 79 ) = W( 79 ) + a*JVS( 532 )
  W( 80 ) = W( 80 ) + a*JVS( 533 )
  W( 82 ) = W( 82 ) + a*JVS( 534 )
  W( 83 ) = W( 83 ) + a*JVS( 535 )
  W( 84 ) = W( 84 ) + a*JVS( 536 )
  W( 85 ) = W( 85 ) + a*JVS( 537 )
  W( 87 ) = W( 87 ) + a*JVS( 538 )
  W( 88 ) = W( 88 ) + a*JVS( 539 )
  W( 89 ) = W( 89 ) + a*JVS( 540 )
  a = -W( 76 ) / JVS(          636  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 637 )
  W( 79 ) = W( 79 ) + a*JVS( 638 )
  W( 82 ) = W( 82 ) + a*JVS( 639 )
  W( 83 ) = W( 83 ) + a*JVS( 640 )
  W( 84 ) = W( 84 ) + a*JVS( 641 )
  W( 85 ) = W( 85 ) + a*JVS( 642 )
  W( 87 ) = W( 87 ) + a*JVS( 643 )
  W( 89 ) = W( 89 ) + a*JVS( 644 )
  a = -W( 77 ) / JVS(          654  )
  W( 77 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 655 )
  W( 82 ) = W( 82 ) + a*JVS( 656 )
  W( 83 ) = W( 83 ) + a*JVS( 657 )
  W( 84 ) = W( 84 ) + a*JVS( 658 )
  W( 85 ) = W( 85 ) + a*JVS( 659 )
  W( 87 ) = W( 87 ) + a*JVS( 660 )
  W( 89 ) = W( 89 ) + a*JVS( 661 )
  a = -W( 78 ) / JVS(          664  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 665 )
  W( 82 ) = W( 82 ) + a*JVS( 666 )
  W( 83 ) = W( 83 ) + a*JVS( 667 )
  W( 84 ) = W( 84 ) + a*JVS( 668 )
  W( 85 ) = W( 85 ) + a*JVS( 669 )
  W( 87 ) = W( 87 ) + a*JVS( 670 )
  W( 88 ) = W( 88 ) + a*JVS( 671 )
  W( 89 ) = W( 89 ) + a*JVS( 672 )
  a = -W( 79 ) / JVS(          693  )
  W( 79 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 83 ) = W( 83 ) + a*JVS( 695 )
  W( 84 ) = W( 84 ) + a*JVS( 696 )
  W( 85 ) = W( 85 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  JVS( 701) = W( 22 )
  JVS( 702) = W( 70 )
  JVS( 703) = W( 76 )
  JVS( 704) = W( 77 )
  JVS( 705) = W( 78 )
  JVS( 706) = W( 79 )
  JVS( 707) = W( 80 )
  JVS( 708) = W( 82 )
  JVS( 709) = W( 83 )
  JVS( 710) = W( 84 )
  JVS( 711) = W( 85 )
  JVS( 712) = W( 86 )
  JVS( 713) = W( 87 )
  JVS( 714) = W( 88 )
  JVS( 715) = W( 89 )
  IF ( ABS(  JVS( 733 )) < TINY(a) ) THEN
         IER = 81                                      
         RETURN
  END IF
   W( 45 ) = JVS( 716 )
   W( 47 ) = JVS( 717 )
   W( 49 ) = JVS( 718 )
   W( 51 ) = JVS( 719 )
   W( 59 ) = JVS( 720 )
   W( 64 ) = JVS( 721 )
   W( 65 ) = JVS( 722 )
   W( 69 ) = JVS( 723 )
   W( 70 ) = JVS( 724 )
   W( 71 ) = JVS( 725 )
   W( 73 ) = JVS( 726 )
   W( 75 ) = JVS( 727 )
   W( 76 ) = JVS( 728 )
   W( 77 ) = JVS( 729 )
   W( 78 ) = JVS( 730 )
   W( 79 ) = JVS( 731 )
   W( 80 ) = JVS( 732 )
   W( 81 ) = JVS( 733 )
   W( 82 ) = JVS( 734 )
   W( 83 ) = JVS( 735 )
   W( 84 ) = JVS( 736 )
   W( 85 ) = JVS( 737 )
   W( 86 ) = JVS( 738 )
   W( 87 ) = JVS( 739 )
   W( 88 ) = JVS( 740 )
   W( 89 ) = JVS( 741 )
  a = -W( 45 ) / JVS(          217  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 218 )
  W( 84 ) = W( 84 ) + a*JVS( 219 )
  W( 89 ) = W( 89 ) + a*JVS( 220 )
  a = -W( 47 ) / JVS(          246  )
  W( 47 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 247 )
  W( 84 ) = W( 84 ) + a*JVS( 248 )
  W( 89 ) = W( 89 ) + a*JVS( 249 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 255 )
  W( 84 ) = W( 84 ) + a*JVS( 256 )
  W( 89 ) = W( 89 ) + a*JVS( 257 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  a = -W( 59 ) / JVS(          350  )
  W( 59 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 351 )
  W( 65 ) = W( 65 ) + a*JVS( 352 )
  W( 69 ) = W( 69 ) + a*JVS( 353 )
  W( 70 ) = W( 70 ) + a*JVS( 354 )
  W( 71 ) = W( 71 ) + a*JVS( 355 )
  W( 73 ) = W( 73 ) + a*JVS( 356 )
  W( 75 ) = W( 75 ) + a*JVS( 357 )
  W( 76 ) = W( 76 ) + a*JVS( 358 )
  W( 77 ) = W( 77 ) + a*JVS( 359 )
  W( 79 ) = W( 79 ) + a*JVS( 360 )
  W( 80 ) = W( 80 ) + a*JVS( 361 )
  W( 81 ) = W( 81 ) + a*JVS( 362 )
  W( 82 ) = W( 82 ) + a*JVS( 363 )
  W( 83 ) = W( 83 ) + a*JVS( 364 )
  W( 84 ) = W( 84 ) + a*JVS( 365 )
  W( 85 ) = W( 85 ) + a*JVS( 366 )
  W( 86 ) = W( 86 ) + a*JVS( 367 )
  W( 87 ) = W( 87 ) + a*JVS( 368 )
  W( 89 ) = W( 89 ) + a*JVS( 369 )
  a = -W( 64 ) / JVS(          426  )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 85 ) = W( 85 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 89 ) = W( 89 ) + a*JVS( 432 )
  a = -W( 65 ) / JVS(          434  )
  W( 65 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 435 )
  W( 83 ) = W( 83 ) + a*JVS( 436 )
  W( 84 ) = W( 84 ) + a*JVS( 437 )
  W( 85 ) = W( 85 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 89 ) = W( 89 ) + a*JVS( 440 )
  a = -W( 69 ) / JVS(          510  )
  W( 69 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 511 )
  W( 79 ) = W( 79 ) + a*JVS( 512 )
  W( 82 ) = W( 82 ) + a*JVS( 513 )
  W( 83 ) = W( 83 ) + a*JVS( 514 )
  W( 84 ) = W( 84 ) + a*JVS( 515 )
  W( 85 ) = W( 85 ) + a*JVS( 516 )
  W( 87 ) = W( 87 ) + a*JVS( 517 )
  W( 89 ) = W( 89 ) + a*JVS( 518 )
  a = -W( 70 ) / JVS(          528  )
  W( 70 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 529 )
  W( 77 ) = W( 77 ) + a*JVS( 530 )
  W( 78 ) = W( 78 ) + a*JVS( 531 )
  W( 79 ) = W( 79 ) + a*JVS( 532 )
  W( 80 ) = W( 80 ) + a*JVS( 533 )
  W( 82 ) = W( 82 ) + a*JVS( 534 )
  W( 83 ) = W( 83 ) + a*JVS( 535 )
  W( 84 ) = W( 84 ) + a*JVS( 536 )
  W( 85 ) = W( 85 ) + a*JVS( 537 )
  W( 87 ) = W( 87 ) + a*JVS( 538 )
  W( 88 ) = W( 88 ) + a*JVS( 539 )
  W( 89 ) = W( 89 ) + a*JVS( 540 )
  a = -W( 71 ) / JVS(          543  )
  W( 71 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 544 )
  W( 82 ) = W( 82 ) + a*JVS( 545 )
  W( 83 ) = W( 83 ) + a*JVS( 546 )
  W( 84 ) = W( 84 ) + a*JVS( 547 )
  W( 85 ) = W( 85 ) + a*JVS( 548 )
  W( 87 ) = W( 87 ) + a*JVS( 549 )
  W( 88 ) = W( 88 ) + a*JVS( 550 )
  W( 89 ) = W( 89 ) + a*JVS( 551 )
  a = -W( 73 ) / JVS(          577  )
  W( 73 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 578 )
  W( 82 ) = W( 82 ) + a*JVS( 579 )
  W( 83 ) = W( 83 ) + a*JVS( 580 )
  W( 84 ) = W( 84 ) + a*JVS( 581 )
  W( 85 ) = W( 85 ) + a*JVS( 582 )
  W( 87 ) = W( 87 ) + a*JVS( 583 )
  W( 89 ) = W( 89 ) + a*JVS( 584 )
  a = -W( 75 ) / JVS(          620  )
  W( 75 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 621 )
  W( 82 ) = W( 82 ) + a*JVS( 622 )
  W( 83 ) = W( 83 ) + a*JVS( 623 )
  W( 84 ) = W( 84 ) + a*JVS( 624 )
  W( 85 ) = W( 85 ) + a*JVS( 625 )
  W( 87 ) = W( 87 ) + a*JVS( 626 )
  W( 89 ) = W( 89 ) + a*JVS( 627 )
  a = -W( 76 ) / JVS(          636  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 637 )
  W( 79 ) = W( 79 ) + a*JVS( 638 )
  W( 82 ) = W( 82 ) + a*JVS( 639 )
  W( 83 ) = W( 83 ) + a*JVS( 640 )
  W( 84 ) = W( 84 ) + a*JVS( 641 )
  W( 85 ) = W( 85 ) + a*JVS( 642 )
  W( 87 ) = W( 87 ) + a*JVS( 643 )
  W( 89 ) = W( 89 ) + a*JVS( 644 )
  a = -W( 77 ) / JVS(          654  )
  W( 77 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 655 )
  W( 82 ) = W( 82 ) + a*JVS( 656 )
  W( 83 ) = W( 83 ) + a*JVS( 657 )
  W( 84 ) = W( 84 ) + a*JVS( 658 )
  W( 85 ) = W( 85 ) + a*JVS( 659 )
  W( 87 ) = W( 87 ) + a*JVS( 660 )
  W( 89 ) = W( 89 ) + a*JVS( 661 )
  a = -W( 78 ) / JVS(          664  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 665 )
  W( 82 ) = W( 82 ) + a*JVS( 666 )
  W( 83 ) = W( 83 ) + a*JVS( 667 )
  W( 84 ) = W( 84 ) + a*JVS( 668 )
  W( 85 ) = W( 85 ) + a*JVS( 669 )
  W( 87 ) = W( 87 ) + a*JVS( 670 )
  W( 88 ) = W( 88 ) + a*JVS( 671 )
  W( 89 ) = W( 89 ) + a*JVS( 672 )
  a = -W( 79 ) / JVS(          693  )
  W( 79 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 83 ) = W( 83 ) + a*JVS( 695 )
  W( 84 ) = W( 84 ) + a*JVS( 696 )
  W( 85 ) = W( 85 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  a = -W( 80 ) / JVS(          707  )
  W( 80 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 708 )
  W( 83 ) = W( 83 ) + a*JVS( 709 )
  W( 84 ) = W( 84 ) + a*JVS( 710 )
  W( 85 ) = W( 85 ) + a*JVS( 711 )
  W( 86 ) = W( 86 ) + a*JVS( 712 )
  W( 87 ) = W( 87 ) + a*JVS( 713 )
  W( 88 ) = W( 88 ) + a*JVS( 714 )
  W( 89 ) = W( 89 ) + a*JVS( 715 )
  JVS( 716) = W( 45 )
  JVS( 717) = W( 47 )
  JVS( 718) = W( 49 )
  JVS( 719) = W( 51 )
  JVS( 720) = W( 59 )
  JVS( 721) = W( 64 )
  JVS( 722) = W( 65 )
  JVS( 723) = W( 69 )
  JVS( 724) = W( 70 )
  JVS( 725) = W( 71 )
  JVS( 726) = W( 73 )
  JVS( 727) = W( 75 )
  JVS( 728) = W( 76 )
  JVS( 729) = W( 77 )
  JVS( 730) = W( 78 )
  JVS( 731) = W( 79 )
  JVS( 732) = W( 80 )
  JVS( 733) = W( 81 )
  JVS( 734) = W( 82 )
  JVS( 735) = W( 83 )
  JVS( 736) = W( 84 )
  JVS( 737) = W( 85 )
  JVS( 738) = W( 86 )
  JVS( 739) = W( 87 )
  JVS( 740) = W( 88 )
  JVS( 741) = W( 89 )
  IF ( ABS(  JVS( 775 )) < TINY(a) ) THEN
         IER = 82                                      
         RETURN
  END IF
   W( 32 ) = JVS( 742 )
   W( 40 ) = JVS( 743 )
   W( 44 ) = JVS( 744 )
   W( 48 ) = JVS( 745 )
   W( 51 ) = JVS( 746 )
   W( 52 ) = JVS( 747 )
   W( 54 ) = JVS( 748 )
   W( 55 ) = JVS( 749 )
   W( 56 ) = JVS( 750 )
   W( 57 ) = JVS( 751 )
   W( 58 ) = JVS( 752 )
   W( 59 ) = JVS( 753 )
   W( 60 ) = JVS( 754 )
   W( 61 ) = JVS( 755 )
   W( 62 ) = JVS( 756 )
   W( 63 ) = JVS( 757 )
   W( 64 ) = JVS( 758 )
   W( 65 ) = JVS( 759 )
   W( 67 ) = JVS( 760 )
   W( 68 ) = JVS( 761 )
   W( 69 ) = JVS( 762 )
   W( 70 ) = JVS( 763 )
   W( 71 ) = JVS( 764 )
   W( 72 ) = JVS( 765 )
   W( 73 ) = JVS( 766 )
   W( 74 ) = JVS( 767 )
   W( 75 ) = JVS( 768 )
   W( 76 ) = JVS( 769 )
   W( 77 ) = JVS( 770 )
   W( 78 ) = JVS( 771 )
   W( 79 ) = JVS( 772 )
   W( 80 ) = JVS( 773 )
   W( 81 ) = JVS( 774 )
   W( 82 ) = JVS( 775 )
   W( 83 ) = JVS( 776 )
   W( 84 ) = JVS( 777 )
   W( 85 ) = JVS( 778 )
   W( 86 ) = JVS( 779 )
   W( 87 ) = JVS( 780 )
   W( 88 ) = JVS( 781 )
   W( 89 ) = JVS( 782 )
   W( 90 ) = JVS( 783 )
  a = -W( 32 ) / JVS(          135  )
  W( 32 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 136 )
  W( 81 ) = W( 81 ) + a*JVS( 137 )
  W( 83 ) = W( 83 ) + a*JVS( 138 )
  W( 84 ) = W( 84 ) + a*JVS( 139 )
  a = -W( 40 ) / JVS(          183  )
  W( 40 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 184 )
  W( 79 ) = W( 79 ) + a*JVS( 185 )
  W( 82 ) = W( 82 ) + a*JVS( 186 )
  W( 84 ) = W( 84 ) + a*JVS( 187 )
  W( 85 ) = W( 85 ) + a*JVS( 188 )
  a = -W( 44 ) / JVS(          210  )
  W( 44 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 211 )
  W( 79 ) = W( 79 ) + a*JVS( 212 )
  W( 82 ) = W( 82 ) + a*JVS( 213 )
  W( 84 ) = W( 84 ) + a*JVS( 214 )
  W( 88 ) = W( 88 ) + a*JVS( 215 )
  W( 89 ) = W( 89 ) + a*JVS( 216 )
  a = -W( 48 ) / JVS(          250  )
  W( 48 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 251 )
  W( 84 ) = W( 84 ) + a*JVS( 252 )
  W( 89 ) = W( 89 ) + a*JVS( 253 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  a = -W( 52 ) / JVS(          275  )
  W( 52 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 276 )
  W( 79 ) = W( 79 ) + a*JVS( 277 )
  W( 84 ) = W( 84 ) + a*JVS( 278 )
  W( 85 ) = W( 85 ) + a*JVS( 279 )
  W( 87 ) = W( 87 ) + a*JVS( 280 )
  W( 88 ) = W( 88 ) + a*JVS( 281 )
  W( 89 ) = W( 89 ) + a*JVS( 282 )
  a = -W( 54 ) / JVS(          296  )
  W( 54 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 297 )
  W( 79 ) = W( 79 ) + a*JVS( 298 )
  W( 84 ) = W( 84 ) + a*JVS( 299 )
  W( 88 ) = W( 88 ) + a*JVS( 300 )
  W( 89 ) = W( 89 ) + a*JVS( 301 )
  a = -W( 55 ) / JVS(          303  )
  W( 55 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 304 )
  W( 82 ) = W( 82 ) + a*JVS( 305 )
  W( 83 ) = W( 83 ) + a*JVS( 306 )
  W( 84 ) = W( 84 ) + a*JVS( 307 )
  W( 85 ) = W( 85 ) + a*JVS( 308 )
  W( 87 ) = W( 87 ) + a*JVS( 309 )
  W( 89 ) = W( 89 ) + a*JVS( 310 )
  a = -W( 56 ) / JVS(          311  )
  W( 56 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 312 )
  W( 82 ) = W( 82 ) + a*JVS( 313 )
  W( 83 ) = W( 83 ) + a*JVS( 314 )
  W( 84 ) = W( 84 ) + a*JVS( 315 )
  W( 85 ) = W( 85 ) + a*JVS( 316 )
  W( 87 ) = W( 87 ) + a*JVS( 317 )
  W( 89 ) = W( 89 ) + a*JVS( 318 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 325 )
  W( 63 ) = W( 63 ) + a*JVS( 326 )
  W( 79 ) = W( 79 ) + a*JVS( 327 )
  W( 82 ) = W( 82 ) + a*JVS( 328 )
  W( 83 ) = W( 83 ) + a*JVS( 329 )
  W( 84 ) = W( 84 ) + a*JVS( 330 )
  W( 85 ) = W( 85 ) + a*JVS( 331 )
  W( 87 ) = W( 87 ) + a*JVS( 332 )
  W( 89 ) = W( 89 ) + a*JVS( 333 )
  a = -W( 58 ) / JVS(          335  )
  W( 58 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 336 )
  W( 82 ) = W( 82 ) + a*JVS( 337 )
  W( 83 ) = W( 83 ) + a*JVS( 338 )
  W( 84 ) = W( 84 ) + a*JVS( 339 )
  W( 85 ) = W( 85 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 89 ) = W( 89 ) + a*JVS( 342 )
  a = -W( 59 ) / JVS(          350  )
  W( 59 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 351 )
  W( 65 ) = W( 65 ) + a*JVS( 352 )
  W( 69 ) = W( 69 ) + a*JVS( 353 )
  W( 70 ) = W( 70 ) + a*JVS( 354 )
  W( 71 ) = W( 71 ) + a*JVS( 355 )
  W( 73 ) = W( 73 ) + a*JVS( 356 )
  W( 75 ) = W( 75 ) + a*JVS( 357 )
  W( 76 ) = W( 76 ) + a*JVS( 358 )
  W( 77 ) = W( 77 ) + a*JVS( 359 )
  W( 79 ) = W( 79 ) + a*JVS( 360 )
  W( 80 ) = W( 80 ) + a*JVS( 361 )
  W( 81 ) = W( 81 ) + a*JVS( 362 )
  W( 82 ) = W( 82 ) + a*JVS( 363 )
  W( 83 ) = W( 83 ) + a*JVS( 364 )
  W( 84 ) = W( 84 ) + a*JVS( 365 )
  W( 85 ) = W( 85 ) + a*JVS( 366 )
  W( 86 ) = W( 86 ) + a*JVS( 367 )
  W( 87 ) = W( 87 ) + a*JVS( 368 )
  W( 89 ) = W( 89 ) + a*JVS( 369 )
  a = -W( 60 ) / JVS(          372  )
  W( 60 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 82 ) = W( 82 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 85 ) = W( 85 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  a = -W( 61 ) / JVS(          382  )
  W( 61 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 82 ) = W( 82 ) + a*JVS( 384 )
  W( 83 ) = W( 83 ) + a*JVS( 385 )
  W( 84 ) = W( 84 ) + a*JVS( 386 )
  W( 85 ) = W( 85 ) + a*JVS( 387 )
  W( 87 ) = W( 87 ) + a*JVS( 388 )
  W( 88 ) = W( 88 ) + a*JVS( 389 )
  W( 89 ) = W( 89 ) + a*JVS( 390 )
  a = -W( 62 ) / JVS(          400  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 71 ) = W( 71 ) + a*JVS( 404 )
  W( 78 ) = W( 78 ) + a*JVS( 405 )
  W( 79 ) = W( 79 ) + a*JVS( 406 )
  W( 80 ) = W( 80 ) + a*JVS( 407 )
  W( 81 ) = W( 81 ) + a*JVS( 408 )
  W( 82 ) = W( 82 ) + a*JVS( 409 )
  W( 83 ) = W( 83 ) + a*JVS( 410 )
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 85 ) = W( 85 ) + a*JVS( 412 )
  W( 87 ) = W( 87 ) + a*JVS( 413 )
  W( 88 ) = W( 88 ) + a*JVS( 414 )
  W( 89 ) = W( 89 ) + a*JVS( 415 )
  a = -W( 63 ) / JVS(          417  )
  W( 63 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 418 )
  W( 82 ) = W( 82 ) + a*JVS( 419 )
  W( 83 ) = W( 83 ) + a*JVS( 420 )
  W( 84 ) = W( 84 ) + a*JVS( 421 )
  W( 85 ) = W( 85 ) + a*JVS( 422 )
  W( 87 ) = W( 87 ) + a*JVS( 423 )
  W( 89 ) = W( 89 ) + a*JVS( 424 )
  a = -W( 64 ) / JVS(          426  )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 85 ) = W( 85 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 89 ) = W( 89 ) + a*JVS( 432 )
  a = -W( 65 ) / JVS(          434  )
  W( 65 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 435 )
  W( 83 ) = W( 83 ) + a*JVS( 436 )
  W( 84 ) = W( 84 ) + a*JVS( 437 )
  W( 85 ) = W( 85 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 89 ) = W( 89 ) + a*JVS( 440 )
  a = -W( 67 ) / JVS(          489  )
  W( 67 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 490 )
  W( 82 ) = W( 82 ) + a*JVS( 491 )
  W( 83 ) = W( 83 ) + a*JVS( 492 )
  W( 84 ) = W( 84 ) + a*JVS( 493 )
  W( 85 ) = W( 85 ) + a*JVS( 494 )
  W( 87 ) = W( 87 ) + a*JVS( 495 )
  W( 88 ) = W( 88 ) + a*JVS( 496 )
  W( 89 ) = W( 89 ) + a*JVS( 497 )
  a = -W( 68 ) / JVS(          499  )
  W( 68 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 500 )
  W( 82 ) = W( 82 ) + a*JVS( 501 )
  W( 83 ) = W( 83 ) + a*JVS( 502 )
  W( 84 ) = W( 84 ) + a*JVS( 503 )
  W( 85 ) = W( 85 ) + a*JVS( 504 )
  W( 87 ) = W( 87 ) + a*JVS( 505 )
  W( 88 ) = W( 88 ) + a*JVS( 506 )
  W( 89 ) = W( 89 ) + a*JVS( 507 )
  a = -W( 69 ) / JVS(          510  )
  W( 69 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 511 )
  W( 79 ) = W( 79 ) + a*JVS( 512 )
  W( 82 ) = W( 82 ) + a*JVS( 513 )
  W( 83 ) = W( 83 ) + a*JVS( 514 )
  W( 84 ) = W( 84 ) + a*JVS( 515 )
  W( 85 ) = W( 85 ) + a*JVS( 516 )
  W( 87 ) = W( 87 ) + a*JVS( 517 )
  W( 89 ) = W( 89 ) + a*JVS( 518 )
  a = -W( 70 ) / JVS(          528  )
  W( 70 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 529 )
  W( 77 ) = W( 77 ) + a*JVS( 530 )
  W( 78 ) = W( 78 ) + a*JVS( 531 )
  W( 79 ) = W( 79 ) + a*JVS( 532 )
  W( 80 ) = W( 80 ) + a*JVS( 533 )
  W( 82 ) = W( 82 ) + a*JVS( 534 )
  W( 83 ) = W( 83 ) + a*JVS( 535 )
  W( 84 ) = W( 84 ) + a*JVS( 536 )
  W( 85 ) = W( 85 ) + a*JVS( 537 )
  W( 87 ) = W( 87 ) + a*JVS( 538 )
  W( 88 ) = W( 88 ) + a*JVS( 539 )
  W( 89 ) = W( 89 ) + a*JVS( 540 )
  a = -W( 71 ) / JVS(          543  )
  W( 71 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 544 )
  W( 82 ) = W( 82 ) + a*JVS( 545 )
  W( 83 ) = W( 83 ) + a*JVS( 546 )
  W( 84 ) = W( 84 ) + a*JVS( 547 )
  W( 85 ) = W( 85 ) + a*JVS( 548 )
  W( 87 ) = W( 87 ) + a*JVS( 549 )
  W( 88 ) = W( 88 ) + a*JVS( 550 )
  W( 89 ) = W( 89 ) + a*JVS( 551 )
  a = -W( 72 ) / JVS(          562  )
  W( 72 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 563 )
  W( 79 ) = W( 79 ) + a*JVS( 564 )
  W( 80 ) = W( 80 ) + a*JVS( 565 )
  W( 81 ) = W( 81 ) + a*JVS( 566 )
  W( 82 ) = W( 82 ) + a*JVS( 567 )
  W( 83 ) = W( 83 ) + a*JVS( 568 )
  W( 84 ) = W( 84 ) + a*JVS( 569 )
  W( 85 ) = W( 85 ) + a*JVS( 570 )
  W( 86 ) = W( 86 ) + a*JVS( 571 )
  W( 87 ) = W( 87 ) + a*JVS( 572 )
  W( 88 ) = W( 88 ) + a*JVS( 573 )
  W( 89 ) = W( 89 ) + a*JVS( 574 )
  a = -W( 73 ) / JVS(          577  )
  W( 73 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 578 )
  W( 82 ) = W( 82 ) + a*JVS( 579 )
  W( 83 ) = W( 83 ) + a*JVS( 580 )
  W( 84 ) = W( 84 ) + a*JVS( 581 )
  W( 85 ) = W( 85 ) + a*JVS( 582 )
  W( 87 ) = W( 87 ) + a*JVS( 583 )
  W( 89 ) = W( 89 ) + a*JVS( 584 )
  a = -W( 74 ) / JVS(          600  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 601 )
  W( 76 ) = W( 76 ) + a*JVS( 602 )
  W( 77 ) = W( 77 ) + a*JVS( 603 )
  W( 78 ) = W( 78 ) + a*JVS( 604 )
  W( 79 ) = W( 79 ) + a*JVS( 605 )
  W( 80 ) = W( 80 ) + a*JVS( 606 )
  W( 81 ) = W( 81 ) + a*JVS( 607 )
  W( 82 ) = W( 82 ) + a*JVS( 608 )
  W( 83 ) = W( 83 ) + a*JVS( 609 )
  W( 84 ) = W( 84 ) + a*JVS( 610 )
  W( 85 ) = W( 85 ) + a*JVS( 611 )
  W( 86 ) = W( 86 ) + a*JVS( 612 )
  W( 87 ) = W( 87 ) + a*JVS( 613 )
  W( 88 ) = W( 88 ) + a*JVS( 614 )
  W( 89 ) = W( 89 ) + a*JVS( 615 )
  W( 90 ) = W( 90 ) + a*JVS( 616 )
  a = -W( 75 ) / JVS(          620  )
  W( 75 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 621 )
  W( 82 ) = W( 82 ) + a*JVS( 622 )
  W( 83 ) = W( 83 ) + a*JVS( 623 )
  W( 84 ) = W( 84 ) + a*JVS( 624 )
  W( 85 ) = W( 85 ) + a*JVS( 625 )
  W( 87 ) = W( 87 ) + a*JVS( 626 )
  W( 89 ) = W( 89 ) + a*JVS( 627 )
  a = -W( 76 ) / JVS(          636  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 637 )
  W( 79 ) = W( 79 ) + a*JVS( 638 )
  W( 82 ) = W( 82 ) + a*JVS( 639 )
  W( 83 ) = W( 83 ) + a*JVS( 640 )
  W( 84 ) = W( 84 ) + a*JVS( 641 )
  W( 85 ) = W( 85 ) + a*JVS( 642 )
  W( 87 ) = W( 87 ) + a*JVS( 643 )
  W( 89 ) = W( 89 ) + a*JVS( 644 )
  a = -W( 77 ) / JVS(          654  )
  W( 77 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 655 )
  W( 82 ) = W( 82 ) + a*JVS( 656 )
  W( 83 ) = W( 83 ) + a*JVS( 657 )
  W( 84 ) = W( 84 ) + a*JVS( 658 )
  W( 85 ) = W( 85 ) + a*JVS( 659 )
  W( 87 ) = W( 87 ) + a*JVS( 660 )
  W( 89 ) = W( 89 ) + a*JVS( 661 )
  a = -W( 78 ) / JVS(          664  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 665 )
  W( 82 ) = W( 82 ) + a*JVS( 666 )
  W( 83 ) = W( 83 ) + a*JVS( 667 )
  W( 84 ) = W( 84 ) + a*JVS( 668 )
  W( 85 ) = W( 85 ) + a*JVS( 669 )
  W( 87 ) = W( 87 ) + a*JVS( 670 )
  W( 88 ) = W( 88 ) + a*JVS( 671 )
  W( 89 ) = W( 89 ) + a*JVS( 672 )
  a = -W( 79 ) / JVS(          693  )
  W( 79 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 83 ) = W( 83 ) + a*JVS( 695 )
  W( 84 ) = W( 84 ) + a*JVS( 696 )
  W( 85 ) = W( 85 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  a = -W( 80 ) / JVS(          707  )
  W( 80 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 708 )
  W( 83 ) = W( 83 ) + a*JVS( 709 )
  W( 84 ) = W( 84 ) + a*JVS( 710 )
  W( 85 ) = W( 85 ) + a*JVS( 711 )
  W( 86 ) = W( 86 ) + a*JVS( 712 )
  W( 87 ) = W( 87 ) + a*JVS( 713 )
  W( 88 ) = W( 88 ) + a*JVS( 714 )
  W( 89 ) = W( 89 ) + a*JVS( 715 )
  a = -W( 81 ) / JVS(          733  )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 734 )
  W( 83 ) = W( 83 ) + a*JVS( 735 )
  W( 84 ) = W( 84 ) + a*JVS( 736 )
  W( 85 ) = W( 85 ) + a*JVS( 737 )
  W( 86 ) = W( 86 ) + a*JVS( 738 )
  W( 87 ) = W( 87 ) + a*JVS( 739 )
  W( 88 ) = W( 88 ) + a*JVS( 740 )
  W( 89 ) = W( 89 ) + a*JVS( 741 )
  JVS( 742) = W( 32 )
  JVS( 743) = W( 40 )
  JVS( 744) = W( 44 )
  JVS( 745) = W( 48 )
  JVS( 746) = W( 51 )
  JVS( 747) = W( 52 )
  JVS( 748) = W( 54 )
  JVS( 749) = W( 55 )
  JVS( 750) = W( 56 )
  JVS( 751) = W( 57 )
  JVS( 752) = W( 58 )
  JVS( 753) = W( 59 )
  JVS( 754) = W( 60 )
  JVS( 755) = W( 61 )
  JVS( 756) = W( 62 )
  JVS( 757) = W( 63 )
  JVS( 758) = W( 64 )
  JVS( 759) = W( 65 )
  JVS( 760) = W( 67 )
  JVS( 761) = W( 68 )
  JVS( 762) = W( 69 )
  JVS( 763) = W( 70 )
  JVS( 764) = W( 71 )
  JVS( 765) = W( 72 )
  JVS( 766) = W( 73 )
  JVS( 767) = W( 74 )
  JVS( 768) = W( 75 )
  JVS( 769) = W( 76 )
  JVS( 770) = W( 77 )
  JVS( 771) = W( 78 )
  JVS( 772) = W( 79 )
  JVS( 773) = W( 80 )
  JVS( 774) = W( 81 )
  JVS( 775) = W( 82 )
  JVS( 776) = W( 83 )
  JVS( 777) = W( 84 )
  JVS( 778) = W( 85 )
  JVS( 779) = W( 86 )
  JVS( 780) = W( 87 )
  JVS( 781) = W( 88 )
  JVS( 782) = W( 89 )
  JVS( 783) = W( 90 )
  IF ( ABS(  JVS( 813 )) < TINY(a) ) THEN
         IER = 83                                      
         RETURN
  END IF
   W( 26 ) = JVS( 784 )
   W( 27 ) = JVS( 785 )
   W( 40 ) = JVS( 786 )
   W( 48 ) = JVS( 787 )
   W( 49 ) = JVS( 788 )
   W( 51 ) = JVS( 789 )
   W( 55 ) = JVS( 790 )
   W( 56 ) = JVS( 791 )
   W( 58 ) = JVS( 792 )
   W( 60 ) = JVS( 793 )
   W( 61 ) = JVS( 794 )
   W( 63 ) = JVS( 795 )
   W( 64 ) = JVS( 796 )
   W( 65 ) = JVS( 797 )
   W( 67 ) = JVS( 798 )
   W( 68 ) = JVS( 799 )
   W( 69 ) = JVS( 800 )
   W( 71 ) = JVS( 801 )
   W( 72 ) = JVS( 802 )
   W( 73 ) = JVS( 803 )
   W( 74 ) = JVS( 804 )
   W( 75 ) = JVS( 805 )
   W( 76 ) = JVS( 806 )
   W( 77 ) = JVS( 807 )
   W( 78 ) = JVS( 808 )
   W( 79 ) = JVS( 809 )
   W( 80 ) = JVS( 810 )
   W( 81 ) = JVS( 811 )
   W( 82 ) = JVS( 812 )
   W( 83 ) = JVS( 813 )
   W( 84 ) = JVS( 814 )
   W( 85 ) = JVS( 815 )
   W( 86 ) = JVS( 816 )
   W( 87 ) = JVS( 817 )
   W( 88 ) = JVS( 818 )
   W( 89 ) = JVS( 819 )
   W( 90 ) = JVS( 820 )
  a = -W( 26 ) / JVS(          102  )
  W( 26 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 103 )
  W( 73 ) = W( 73 ) + a*JVS( 104 )
  W( 75 ) = W( 75 ) + a*JVS( 105 )
  W( 79 ) = W( 79 ) + a*JVS( 106 )
  W( 84 ) = W( 84 ) + a*JVS( 107 )
  a = -W( 27 ) / JVS(          108  )
  W( 27 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 109 )
  W( 84 ) = W( 84 ) + a*JVS( 110 )
  W( 85 ) = W( 85 ) + a*JVS( 111 )
  a = -W( 40 ) / JVS(          183  )
  W( 40 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 184 )
  W( 79 ) = W( 79 ) + a*JVS( 185 )
  W( 82 ) = W( 82 ) + a*JVS( 186 )
  W( 84 ) = W( 84 ) + a*JVS( 187 )
  W( 85 ) = W( 85 ) + a*JVS( 188 )
  a = -W( 48 ) / JVS(          250  )
  W( 48 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 251 )
  W( 84 ) = W( 84 ) + a*JVS( 252 )
  W( 89 ) = W( 89 ) + a*JVS( 253 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 255 )
  W( 84 ) = W( 84 ) + a*JVS( 256 )
  W( 89 ) = W( 89 ) + a*JVS( 257 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  a = -W( 55 ) / JVS(          303  )
  W( 55 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 304 )
  W( 82 ) = W( 82 ) + a*JVS( 305 )
  W( 83 ) = W( 83 ) + a*JVS( 306 )
  W( 84 ) = W( 84 ) + a*JVS( 307 )
  W( 85 ) = W( 85 ) + a*JVS( 308 )
  W( 87 ) = W( 87 ) + a*JVS( 309 )
  W( 89 ) = W( 89 ) + a*JVS( 310 )
  a = -W( 56 ) / JVS(          311  )
  W( 56 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 312 )
  W( 82 ) = W( 82 ) + a*JVS( 313 )
  W( 83 ) = W( 83 ) + a*JVS( 314 )
  W( 84 ) = W( 84 ) + a*JVS( 315 )
  W( 85 ) = W( 85 ) + a*JVS( 316 )
  W( 87 ) = W( 87 ) + a*JVS( 317 )
  W( 89 ) = W( 89 ) + a*JVS( 318 )
  a = -W( 58 ) / JVS(          335  )
  W( 58 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 336 )
  W( 82 ) = W( 82 ) + a*JVS( 337 )
  W( 83 ) = W( 83 ) + a*JVS( 338 )
  W( 84 ) = W( 84 ) + a*JVS( 339 )
  W( 85 ) = W( 85 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 89 ) = W( 89 ) + a*JVS( 342 )
  a = -W( 60 ) / JVS(          372  )
  W( 60 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 82 ) = W( 82 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 85 ) = W( 85 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  a = -W( 61 ) / JVS(          382  )
  W( 61 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 82 ) = W( 82 ) + a*JVS( 384 )
  W( 83 ) = W( 83 ) + a*JVS( 385 )
  W( 84 ) = W( 84 ) + a*JVS( 386 )
  W( 85 ) = W( 85 ) + a*JVS( 387 )
  W( 87 ) = W( 87 ) + a*JVS( 388 )
  W( 88 ) = W( 88 ) + a*JVS( 389 )
  W( 89 ) = W( 89 ) + a*JVS( 390 )
  a = -W( 63 ) / JVS(          417  )
  W( 63 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 418 )
  W( 82 ) = W( 82 ) + a*JVS( 419 )
  W( 83 ) = W( 83 ) + a*JVS( 420 )
  W( 84 ) = W( 84 ) + a*JVS( 421 )
  W( 85 ) = W( 85 ) + a*JVS( 422 )
  W( 87 ) = W( 87 ) + a*JVS( 423 )
  W( 89 ) = W( 89 ) + a*JVS( 424 )
  a = -W( 64 ) / JVS(          426  )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 85 ) = W( 85 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 89 ) = W( 89 ) + a*JVS( 432 )
  a = -W( 65 ) / JVS(          434  )
  W( 65 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 435 )
  W( 83 ) = W( 83 ) + a*JVS( 436 )
  W( 84 ) = W( 84 ) + a*JVS( 437 )
  W( 85 ) = W( 85 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 89 ) = W( 89 ) + a*JVS( 440 )
  a = -W( 67 ) / JVS(          489  )
  W( 67 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 490 )
  W( 82 ) = W( 82 ) + a*JVS( 491 )
  W( 83 ) = W( 83 ) + a*JVS( 492 )
  W( 84 ) = W( 84 ) + a*JVS( 493 )
  W( 85 ) = W( 85 ) + a*JVS( 494 )
  W( 87 ) = W( 87 ) + a*JVS( 495 )
  W( 88 ) = W( 88 ) + a*JVS( 496 )
  W( 89 ) = W( 89 ) + a*JVS( 497 )
  a = -W( 68 ) / JVS(          499  )
  W( 68 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 500 )
  W( 82 ) = W( 82 ) + a*JVS( 501 )
  W( 83 ) = W( 83 ) + a*JVS( 502 )
  W( 84 ) = W( 84 ) + a*JVS( 503 )
  W( 85 ) = W( 85 ) + a*JVS( 504 )
  W( 87 ) = W( 87 ) + a*JVS( 505 )
  W( 88 ) = W( 88 ) + a*JVS( 506 )
  W( 89 ) = W( 89 ) + a*JVS( 507 )
  a = -W( 69 ) / JVS(          510  )
  W( 69 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 511 )
  W( 79 ) = W( 79 ) + a*JVS( 512 )
  W( 82 ) = W( 82 ) + a*JVS( 513 )
  W( 83 ) = W( 83 ) + a*JVS( 514 )
  W( 84 ) = W( 84 ) + a*JVS( 515 )
  W( 85 ) = W( 85 ) + a*JVS( 516 )
  W( 87 ) = W( 87 ) + a*JVS( 517 )
  W( 89 ) = W( 89 ) + a*JVS( 518 )
  a = -W( 71 ) / JVS(          543  )
  W( 71 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 544 )
  W( 82 ) = W( 82 ) + a*JVS( 545 )
  W( 83 ) = W( 83 ) + a*JVS( 546 )
  W( 84 ) = W( 84 ) + a*JVS( 547 )
  W( 85 ) = W( 85 ) + a*JVS( 548 )
  W( 87 ) = W( 87 ) + a*JVS( 549 )
  W( 88 ) = W( 88 ) + a*JVS( 550 )
  W( 89 ) = W( 89 ) + a*JVS( 551 )
  a = -W( 72 ) / JVS(          562  )
  W( 72 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 563 )
  W( 79 ) = W( 79 ) + a*JVS( 564 )
  W( 80 ) = W( 80 ) + a*JVS( 565 )
  W( 81 ) = W( 81 ) + a*JVS( 566 )
  W( 82 ) = W( 82 ) + a*JVS( 567 )
  W( 83 ) = W( 83 ) + a*JVS( 568 )
  W( 84 ) = W( 84 ) + a*JVS( 569 )
  W( 85 ) = W( 85 ) + a*JVS( 570 )
  W( 86 ) = W( 86 ) + a*JVS( 571 )
  W( 87 ) = W( 87 ) + a*JVS( 572 )
  W( 88 ) = W( 88 ) + a*JVS( 573 )
  W( 89 ) = W( 89 ) + a*JVS( 574 )
  a = -W( 73 ) / JVS(          577  )
  W( 73 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 578 )
  W( 82 ) = W( 82 ) + a*JVS( 579 )
  W( 83 ) = W( 83 ) + a*JVS( 580 )
  W( 84 ) = W( 84 ) + a*JVS( 581 )
  W( 85 ) = W( 85 ) + a*JVS( 582 )
  W( 87 ) = W( 87 ) + a*JVS( 583 )
  W( 89 ) = W( 89 ) + a*JVS( 584 )
  a = -W( 74 ) / JVS(          600  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 601 )
  W( 76 ) = W( 76 ) + a*JVS( 602 )
  W( 77 ) = W( 77 ) + a*JVS( 603 )
  W( 78 ) = W( 78 ) + a*JVS( 604 )
  W( 79 ) = W( 79 ) + a*JVS( 605 )
  W( 80 ) = W( 80 ) + a*JVS( 606 )
  W( 81 ) = W( 81 ) + a*JVS( 607 )
  W( 82 ) = W( 82 ) + a*JVS( 608 )
  W( 83 ) = W( 83 ) + a*JVS( 609 )
  W( 84 ) = W( 84 ) + a*JVS( 610 )
  W( 85 ) = W( 85 ) + a*JVS( 611 )
  W( 86 ) = W( 86 ) + a*JVS( 612 )
  W( 87 ) = W( 87 ) + a*JVS( 613 )
  W( 88 ) = W( 88 ) + a*JVS( 614 )
  W( 89 ) = W( 89 ) + a*JVS( 615 )
  W( 90 ) = W( 90 ) + a*JVS( 616 )
  a = -W( 75 ) / JVS(          620  )
  W( 75 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 621 )
  W( 82 ) = W( 82 ) + a*JVS( 622 )
  W( 83 ) = W( 83 ) + a*JVS( 623 )
  W( 84 ) = W( 84 ) + a*JVS( 624 )
  W( 85 ) = W( 85 ) + a*JVS( 625 )
  W( 87 ) = W( 87 ) + a*JVS( 626 )
  W( 89 ) = W( 89 ) + a*JVS( 627 )
  a = -W( 76 ) / JVS(          636  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 637 )
  W( 79 ) = W( 79 ) + a*JVS( 638 )
  W( 82 ) = W( 82 ) + a*JVS( 639 )
  W( 83 ) = W( 83 ) + a*JVS( 640 )
  W( 84 ) = W( 84 ) + a*JVS( 641 )
  W( 85 ) = W( 85 ) + a*JVS( 642 )
  W( 87 ) = W( 87 ) + a*JVS( 643 )
  W( 89 ) = W( 89 ) + a*JVS( 644 )
  a = -W( 77 ) / JVS(          654  )
  W( 77 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 655 )
  W( 82 ) = W( 82 ) + a*JVS( 656 )
  W( 83 ) = W( 83 ) + a*JVS( 657 )
  W( 84 ) = W( 84 ) + a*JVS( 658 )
  W( 85 ) = W( 85 ) + a*JVS( 659 )
  W( 87 ) = W( 87 ) + a*JVS( 660 )
  W( 89 ) = W( 89 ) + a*JVS( 661 )
  a = -W( 78 ) / JVS(          664  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 665 )
  W( 82 ) = W( 82 ) + a*JVS( 666 )
  W( 83 ) = W( 83 ) + a*JVS( 667 )
  W( 84 ) = W( 84 ) + a*JVS( 668 )
  W( 85 ) = W( 85 ) + a*JVS( 669 )
  W( 87 ) = W( 87 ) + a*JVS( 670 )
  W( 88 ) = W( 88 ) + a*JVS( 671 )
  W( 89 ) = W( 89 ) + a*JVS( 672 )
  a = -W( 79 ) / JVS(          693  )
  W( 79 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 83 ) = W( 83 ) + a*JVS( 695 )
  W( 84 ) = W( 84 ) + a*JVS( 696 )
  W( 85 ) = W( 85 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  a = -W( 80 ) / JVS(          707  )
  W( 80 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 708 )
  W( 83 ) = W( 83 ) + a*JVS( 709 )
  W( 84 ) = W( 84 ) + a*JVS( 710 )
  W( 85 ) = W( 85 ) + a*JVS( 711 )
  W( 86 ) = W( 86 ) + a*JVS( 712 )
  W( 87 ) = W( 87 ) + a*JVS( 713 )
  W( 88 ) = W( 88 ) + a*JVS( 714 )
  W( 89 ) = W( 89 ) + a*JVS( 715 )
  a = -W( 81 ) / JVS(          733  )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 734 )
  W( 83 ) = W( 83 ) + a*JVS( 735 )
  W( 84 ) = W( 84 ) + a*JVS( 736 )
  W( 85 ) = W( 85 ) + a*JVS( 737 )
  W( 86 ) = W( 86 ) + a*JVS( 738 )
  W( 87 ) = W( 87 ) + a*JVS( 739 )
  W( 88 ) = W( 88 ) + a*JVS( 740 )
  W( 89 ) = W( 89 ) + a*JVS( 741 )
  a = -W( 82 ) / JVS(          775  )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 776 )
  W( 84 ) = W( 84 ) + a*JVS( 777 )
  W( 85 ) = W( 85 ) + a*JVS( 778 )
  W( 86 ) = W( 86 ) + a*JVS( 779 )
  W( 87 ) = W( 87 ) + a*JVS( 780 )
  W( 88 ) = W( 88 ) + a*JVS( 781 )
  W( 89 ) = W( 89 ) + a*JVS( 782 )
  W( 90 ) = W( 90 ) + a*JVS( 783 )
  JVS( 784) = W( 26 )
  JVS( 785) = W( 27 )
  JVS( 786) = W( 40 )
  JVS( 787) = W( 48 )
  JVS( 788) = W( 49 )
  JVS( 789) = W( 51 )
  JVS( 790) = W( 55 )
  JVS( 791) = W( 56 )
  JVS( 792) = W( 58 )
  JVS( 793) = W( 60 )
  JVS( 794) = W( 61 )
  JVS( 795) = W( 63 )
  JVS( 796) = W( 64 )
  JVS( 797) = W( 65 )
  JVS( 798) = W( 67 )
  JVS( 799) = W( 68 )
  JVS( 800) = W( 69 )
  JVS( 801) = W( 71 )
  JVS( 802) = W( 72 )
  JVS( 803) = W( 73 )
  JVS( 804) = W( 74 )
  JVS( 805) = W( 75 )
  JVS( 806) = W( 76 )
  JVS( 807) = W( 77 )
  JVS( 808) = W( 78 )
  JVS( 809) = W( 79 )
  JVS( 810) = W( 80 )
  JVS( 811) = W( 81 )
  JVS( 812) = W( 82 )
  JVS( 813) = W( 83 )
  JVS( 814) = W( 84 )
  JVS( 815) = W( 85 )
  JVS( 816) = W( 86 )
  JVS( 817) = W( 87 )
  JVS( 818) = W( 88 )
  JVS( 819) = W( 89 )
  JVS( 820) = W( 90 )
  IF ( ABS(  JVS( 885 )) < TINY(a) ) THEN
         IER = 84                                      
         RETURN
  END IF
   W( 13 ) = JVS( 821 )
   W( 14 ) = JVS( 822 )
   W( 16 ) = JVS( 823 )
   W( 17 ) = JVS( 824 )
   W( 18 ) = JVS( 825 )
   W( 20 ) = JVS( 826 )
   W( 21 ) = JVS( 827 )
   W( 22 ) = JVS( 828 )
   W( 23 ) = JVS( 829 )
   W( 24 ) = JVS( 830 )
   W( 25 ) = JVS( 831 )
   W( 26 ) = JVS( 832 )
   W( 27 ) = JVS( 833 )
   W( 28 ) = JVS( 834 )
   W( 29 ) = JVS( 835 )
   W( 30 ) = JVS( 836 )
   W( 31 ) = JVS( 837 )
   W( 32 ) = JVS( 838 )
   W( 35 ) = JVS( 839 )
   W( 36 ) = JVS( 840 )
   W( 37 ) = JVS( 841 )
   W( 38 ) = JVS( 842 )
   W( 39 ) = JVS( 843 )
   W( 40 ) = JVS( 844 )
   W( 41 ) = JVS( 845 )
   W( 42 ) = JVS( 846 )
   W( 43 ) = JVS( 847 )
   W( 44 ) = JVS( 848 )
   W( 45 ) = JVS( 849 )
   W( 46 ) = JVS( 850 )
   W( 47 ) = JVS( 851 )
   W( 48 ) = JVS( 852 )
   W( 49 ) = JVS( 853 )
   W( 50 ) = JVS( 854 )
   W( 51 ) = JVS( 855 )
   W( 52 ) = JVS( 856 )
   W( 53 ) = JVS( 857 )
   W( 54 ) = JVS( 858 )
   W( 57 ) = JVS( 859 )
   W( 59 ) = JVS( 860 )
   W( 60 ) = JVS( 861 )
   W( 61 ) = JVS( 862 )
   W( 62 ) = JVS( 863 )
   W( 63 ) = JVS( 864 )
   W( 64 ) = JVS( 865 )
   W( 65 ) = JVS( 866 )
   W( 66 ) = JVS( 867 )
   W( 67 ) = JVS( 868 )
   W( 68 ) = JVS( 869 )
   W( 69 ) = JVS( 870 )
   W( 70 ) = JVS( 871 )
   W( 71 ) = JVS( 872 )
   W( 72 ) = JVS( 873 )
   W( 73 ) = JVS( 874 )
   W( 74 ) = JVS( 875 )
   W( 75 ) = JVS( 876 )
   W( 76 ) = JVS( 877 )
   W( 77 ) = JVS( 878 )
   W( 78 ) = JVS( 879 )
   W( 79 ) = JVS( 880 )
   W( 80 ) = JVS( 881 )
   W( 81 ) = JVS( 882 )
   W( 82 ) = JVS( 883 )
   W( 83 ) = JVS( 884 )
   W( 84 ) = JVS( 885 )
   W( 85 ) = JVS( 886 )
   W( 86 ) = JVS( 887 )
   W( 87 ) = JVS( 888 )
   W( 88 ) = JVS( 889 )
   W( 89 ) = JVS( 890 )
   W( 90 ) = JVS( 891 )
  a = -W( 13 ) / JVS(           65  )
  W( 13 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 66 )
  a = -W( 14 ) / JVS(           67  )
  W( 14 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 68 )
  a = -W( 16 ) / JVS(           73  )
  W( 16 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 74 )
  a = -W( 17 ) / JVS(           75  )
  W( 17 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 76 )
  a = -W( 18 ) / JVS(           77  )
  W( 18 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 78 )
  a = -W( 20 ) / JVS(           82  )
  W( 20 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 83 )
  a = -W( 21 ) / JVS(           84  )
  W( 21 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 85 )
  W( 84 ) = W( 84 ) + a*JVS( 86 )
  W( 85 ) = W( 85 ) + a*JVS( 87 )
  a = -W( 22 ) / JVS(           88  )
  W( 22 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 89 )
  a = -W( 23 ) / JVS(           90  )
  W( 23 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 91 )
  W( 84 ) = W( 84 ) + a*JVS( 92 )
  a = -W( 24 ) / JVS(           93  )
  W( 24 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 94 )
  W( 84 ) = W( 84 ) + a*JVS( 95 )
  a = -W( 25 ) / JVS(           96  )
  W( 25 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 97 )
  W( 73 ) = W( 73 ) + a*JVS( 98 )
  W( 75 ) = W( 75 ) + a*JVS( 99 )
  W( 79 ) = W( 79 ) + a*JVS( 100 )
  W( 84 ) = W( 84 ) + a*JVS( 101 )
  a = -W( 26 ) / JVS(          102  )
  W( 26 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 103 )
  W( 73 ) = W( 73 ) + a*JVS( 104 )
  W( 75 ) = W( 75 ) + a*JVS( 105 )
  W( 79 ) = W( 79 ) + a*JVS( 106 )
  W( 84 ) = W( 84 ) + a*JVS( 107 )
  a = -W( 27 ) / JVS(          108  )
  W( 27 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 109 )
  W( 84 ) = W( 84 ) + a*JVS( 110 )
  W( 85 ) = W( 85 ) + a*JVS( 111 )
  a = -W( 28 ) / JVS(          112  )
  W( 28 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 113 )
  W( 84 ) = W( 84 ) + a*JVS( 114 )
  W( 88 ) = W( 88 ) + a*JVS( 115 )
  a = -W( 29 ) / JVS(          116  )
  W( 29 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 117 )
  W( 85 ) = W( 85 ) + a*JVS( 118 )
  W( 88 ) = W( 88 ) + a*JVS( 119 )
  a = -W( 30 ) / JVS(          121  )
  W( 30 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 122 )
  W( 52 ) = W( 52 ) + a*JVS( 123 )
  W( 84 ) = W( 84 ) + a*JVS( 124 )
  W( 87 ) = W( 87 ) + a*JVS( 125 )
  W( 88 ) = W( 88 ) + a*JVS( 126 )
  a = -W( 31 ) / JVS(          127  )
  W( 31 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 128 )
  W( 37 ) = W( 37 ) + a*JVS( 129 )
  W( 39 ) = W( 39 ) + a*JVS( 130 )
  W( 84 ) = W( 84 ) + a*JVS( 131 )
  W( 87 ) = W( 87 ) + a*JVS( 132 )
  W( 88 ) = W( 88 ) + a*JVS( 133 )
  a = -W( 32 ) / JVS(          135  )
  W( 32 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 136 )
  W( 81 ) = W( 81 ) + a*JVS( 137 )
  W( 83 ) = W( 83 ) + a*JVS( 138 )
  W( 84 ) = W( 84 ) + a*JVS( 139 )
  a = -W( 35 ) / JVS(          153  )
  W( 35 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 154 )
  W( 47 ) = W( 47 ) + a*JVS( 155 )
  W( 48 ) = W( 48 ) + a*JVS( 156 )
  W( 49 ) = W( 49 ) + a*JVS( 157 )
  W( 51 ) = W( 51 ) + a*JVS( 158 )
  W( 73 ) = W( 73 ) + a*JVS( 159 )
  W( 75 ) = W( 75 ) + a*JVS( 160 )
  W( 79 ) = W( 79 ) + a*JVS( 161 )
  W( 84 ) = W( 84 ) + a*JVS( 162 )
  W( 85 ) = W( 85 ) + a*JVS( 163 )
  a = -W( 36 ) / JVS(          165  )
  W( 36 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 166 )
  W( 84 ) = W( 84 ) + a*JVS( 167 )
  W( 88 ) = W( 88 ) + a*JVS( 168 )
  a = -W( 37 ) / JVS(          170  )
  W( 37 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 171 )
  W( 84 ) = W( 84 ) + a*JVS( 172 )
  W( 88 ) = W( 88 ) + a*JVS( 173 )
  a = -W( 38 ) / JVS(          174  )
  W( 38 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 175 )
  W( 84 ) = W( 84 ) + a*JVS( 176 )
  W( 89 ) = W( 89 ) + a*JVS( 177 )
  a = -W( 39 ) / JVS(          178  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 179 )
  W( 79 ) = W( 79 ) + a*JVS( 180 )
  W( 84 ) = W( 84 ) + a*JVS( 181 )
  W( 88 ) = W( 88 ) + a*JVS( 182 )
  a = -W( 40 ) / JVS(          183  )
  W( 40 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 184 )
  W( 79 ) = W( 79 ) + a*JVS( 185 )
  W( 82 ) = W( 82 ) + a*JVS( 186 )
  W( 84 ) = W( 84 ) + a*JVS( 187 )
  W( 85 ) = W( 85 ) + a*JVS( 188 )
  a = -W( 41 ) / JVS(          189  )
  W( 41 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 190 )
  W( 60 ) = W( 60 ) + a*JVS( 191 )
  W( 84 ) = W( 84 ) + a*JVS( 192 )
  W( 87 ) = W( 87 ) + a*JVS( 193 )
  W( 89 ) = W( 89 ) + a*JVS( 194 )
  a = -W( 42 ) / JVS(          195  )
  W( 42 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 196 )
  W( 84 ) = W( 84 ) + a*JVS( 197 )
  W( 89 ) = W( 89 ) + a*JVS( 198 )
  a = -W( 43 ) / JVS(          199  )
  W( 43 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 200 )
  W( 53 ) = W( 53 ) + a*JVS( 201 )
  W( 62 ) = W( 62 ) + a*JVS( 202 )
  W( 66 ) = W( 66 ) + a*JVS( 203 )
  W( 71 ) = W( 71 ) + a*JVS( 204 )
  W( 74 ) = W( 74 ) + a*JVS( 205 )
  W( 84 ) = W( 84 ) + a*JVS( 206 )
  W( 85 ) = W( 85 ) + a*JVS( 207 )
  W( 88 ) = W( 88 ) + a*JVS( 208 )
  W( 89 ) = W( 89 ) + a*JVS( 209 )
  a = -W( 44 ) / JVS(          210  )
  W( 44 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 211 )
  W( 79 ) = W( 79 ) + a*JVS( 212 )
  W( 82 ) = W( 82 ) + a*JVS( 213 )
  W( 84 ) = W( 84 ) + a*JVS( 214 )
  W( 88 ) = W( 88 ) + a*JVS( 215 )
  W( 89 ) = W( 89 ) + a*JVS( 216 )
  a = -W( 45 ) / JVS(          217  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 218 )
  W( 84 ) = W( 84 ) + a*JVS( 219 )
  W( 89 ) = W( 89 ) + a*JVS( 220 )
  a = -W( 46 ) / JVS(          226  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 227 )
  W( 48 ) = W( 48 ) + a*JVS( 228 )
  W( 49 ) = W( 49 ) + a*JVS( 229 )
  W( 51 ) = W( 51 ) + a*JVS( 230 )
  W( 52 ) = W( 52 ) + a*JVS( 231 )
  W( 53 ) = W( 53 ) + a*JVS( 232 )
  W( 54 ) = W( 54 ) + a*JVS( 233 )
  W( 57 ) = W( 57 ) + a*JVS( 234 )
  W( 60 ) = W( 60 ) + a*JVS( 235 )
  W( 62 ) = W( 62 ) + a*JVS( 236 )
  W( 66 ) = W( 66 ) + a*JVS( 237 )
  W( 71 ) = W( 71 ) + a*JVS( 238 )
  W( 73 ) = W( 73 ) + a*JVS( 239 )
  W( 74 ) = W( 74 ) + a*JVS( 240 )
  W( 75 ) = W( 75 ) + a*JVS( 241 )
  W( 79 ) = W( 79 ) + a*JVS( 242 )
  W( 84 ) = W( 84 ) + a*JVS( 243 )
  W( 87 ) = W( 87 ) + a*JVS( 244 )
  W( 89 ) = W( 89 ) + a*JVS( 245 )
  a = -W( 47 ) / JVS(          246  )
  W( 47 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 247 )
  W( 84 ) = W( 84 ) + a*JVS( 248 )
  W( 89 ) = W( 89 ) + a*JVS( 249 )
  a = -W( 48 ) / JVS(          250  )
  W( 48 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 251 )
  W( 84 ) = W( 84 ) + a*JVS( 252 )
  W( 89 ) = W( 89 ) + a*JVS( 253 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 255 )
  W( 84 ) = W( 84 ) + a*JVS( 256 )
  W( 89 ) = W( 89 ) + a*JVS( 257 )
  a = -W( 50 ) / JVS(          262  )
  W( 50 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 263 )
  W( 84 ) = W( 84 ) + a*JVS( 264 )
  W( 85 ) = W( 85 ) + a*JVS( 265 )
  W( 88 ) = W( 88 ) + a*JVS( 266 )
  W( 89 ) = W( 89 ) + a*JVS( 267 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  a = -W( 52 ) / JVS(          275  )
  W( 52 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 276 )
  W( 79 ) = W( 79 ) + a*JVS( 277 )
  W( 84 ) = W( 84 ) + a*JVS( 278 )
  W( 85 ) = W( 85 ) + a*JVS( 279 )
  W( 87 ) = W( 87 ) + a*JVS( 280 )
  W( 88 ) = W( 88 ) + a*JVS( 281 )
  W( 89 ) = W( 89 ) + a*JVS( 282 )
  a = -W( 53 ) / JVS(          284  )
  W( 53 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 285 )
  W( 67 ) = W( 67 ) + a*JVS( 286 )
  W( 68 ) = W( 68 ) + a*JVS( 287 )
  W( 71 ) = W( 71 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  W( 80 ) = W( 80 ) + a*JVS( 290 )
  W( 82 ) = W( 82 ) + a*JVS( 291 )
  W( 83 ) = W( 83 ) + a*JVS( 292 )
  W( 84 ) = W( 84 ) + a*JVS( 293 )
  W( 87 ) = W( 87 ) + a*JVS( 294 )
  W( 89 ) = W( 89 ) + a*JVS( 295 )
  a = -W( 54 ) / JVS(          296  )
  W( 54 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 297 )
  W( 79 ) = W( 79 ) + a*JVS( 298 )
  W( 84 ) = W( 84 ) + a*JVS( 299 )
  W( 88 ) = W( 88 ) + a*JVS( 300 )
  W( 89 ) = W( 89 ) + a*JVS( 301 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 325 )
  W( 63 ) = W( 63 ) + a*JVS( 326 )
  W( 79 ) = W( 79 ) + a*JVS( 327 )
  W( 82 ) = W( 82 ) + a*JVS( 328 )
  W( 83 ) = W( 83 ) + a*JVS( 329 )
  W( 84 ) = W( 84 ) + a*JVS( 330 )
  W( 85 ) = W( 85 ) + a*JVS( 331 )
  W( 87 ) = W( 87 ) + a*JVS( 332 )
  W( 89 ) = W( 89 ) + a*JVS( 333 )
  a = -W( 59 ) / JVS(          350  )
  W( 59 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 351 )
  W( 65 ) = W( 65 ) + a*JVS( 352 )
  W( 69 ) = W( 69 ) + a*JVS( 353 )
  W( 70 ) = W( 70 ) + a*JVS( 354 )
  W( 71 ) = W( 71 ) + a*JVS( 355 )
  W( 73 ) = W( 73 ) + a*JVS( 356 )
  W( 75 ) = W( 75 ) + a*JVS( 357 )
  W( 76 ) = W( 76 ) + a*JVS( 358 )
  W( 77 ) = W( 77 ) + a*JVS( 359 )
  W( 79 ) = W( 79 ) + a*JVS( 360 )
  W( 80 ) = W( 80 ) + a*JVS( 361 )
  W( 81 ) = W( 81 ) + a*JVS( 362 )
  W( 82 ) = W( 82 ) + a*JVS( 363 )
  W( 83 ) = W( 83 ) + a*JVS( 364 )
  W( 84 ) = W( 84 ) + a*JVS( 365 )
  W( 85 ) = W( 85 ) + a*JVS( 366 )
  W( 86 ) = W( 86 ) + a*JVS( 367 )
  W( 87 ) = W( 87 ) + a*JVS( 368 )
  W( 89 ) = W( 89 ) + a*JVS( 369 )
  a = -W( 60 ) / JVS(          372  )
  W( 60 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 82 ) = W( 82 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 85 ) = W( 85 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  a = -W( 61 ) / JVS(          382  )
  W( 61 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 82 ) = W( 82 ) + a*JVS( 384 )
  W( 83 ) = W( 83 ) + a*JVS( 385 )
  W( 84 ) = W( 84 ) + a*JVS( 386 )
  W( 85 ) = W( 85 ) + a*JVS( 387 )
  W( 87 ) = W( 87 ) + a*JVS( 388 )
  W( 88 ) = W( 88 ) + a*JVS( 389 )
  W( 89 ) = W( 89 ) + a*JVS( 390 )
  a = -W( 62 ) / JVS(          400  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 71 ) = W( 71 ) + a*JVS( 404 )
  W( 78 ) = W( 78 ) + a*JVS( 405 )
  W( 79 ) = W( 79 ) + a*JVS( 406 )
  W( 80 ) = W( 80 ) + a*JVS( 407 )
  W( 81 ) = W( 81 ) + a*JVS( 408 )
  W( 82 ) = W( 82 ) + a*JVS( 409 )
  W( 83 ) = W( 83 ) + a*JVS( 410 )
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 85 ) = W( 85 ) + a*JVS( 412 )
  W( 87 ) = W( 87 ) + a*JVS( 413 )
  W( 88 ) = W( 88 ) + a*JVS( 414 )
  W( 89 ) = W( 89 ) + a*JVS( 415 )
  a = -W( 63 ) / JVS(          417  )
  W( 63 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 418 )
  W( 82 ) = W( 82 ) + a*JVS( 419 )
  W( 83 ) = W( 83 ) + a*JVS( 420 )
  W( 84 ) = W( 84 ) + a*JVS( 421 )
  W( 85 ) = W( 85 ) + a*JVS( 422 )
  W( 87 ) = W( 87 ) + a*JVS( 423 )
  W( 89 ) = W( 89 ) + a*JVS( 424 )
  a = -W( 64 ) / JVS(          426  )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 85 ) = W( 85 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 89 ) = W( 89 ) + a*JVS( 432 )
  a = -W( 65 ) / JVS(          434  )
  W( 65 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 435 )
  W( 83 ) = W( 83 ) + a*JVS( 436 )
  W( 84 ) = W( 84 ) + a*JVS( 437 )
  W( 85 ) = W( 85 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 89 ) = W( 89 ) + a*JVS( 440 )
  a = -W( 66 ) / JVS(          466  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 467 )
  W( 68 ) = W( 68 ) + a*JVS( 468 )
  W( 69 ) = W( 69 ) + a*JVS( 469 )
  W( 71 ) = W( 71 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  W( 75 ) = W( 75 ) + a*JVS( 473 )
  W( 76 ) = W( 76 ) + a*JVS( 474 )
  W( 77 ) = W( 77 ) + a*JVS( 475 )
  W( 78 ) = W( 78 ) + a*JVS( 476 )
  W( 79 ) = W( 79 ) + a*JVS( 477 )
  W( 80 ) = W( 80 ) + a*JVS( 478 )
  W( 81 ) = W( 81 ) + a*JVS( 479 )
  W( 82 ) = W( 82 ) + a*JVS( 480 )
  W( 83 ) = W( 83 ) + a*JVS( 481 )
  W( 84 ) = W( 84 ) + a*JVS( 482 )
  W( 85 ) = W( 85 ) + a*JVS( 483 )
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 89 ) = W( 89 ) + a*JVS( 486 )
  W( 90 ) = W( 90 ) + a*JVS( 487 )
  a = -W( 67 ) / JVS(          489  )
  W( 67 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 490 )
  W( 82 ) = W( 82 ) + a*JVS( 491 )
  W( 83 ) = W( 83 ) + a*JVS( 492 )
  W( 84 ) = W( 84 ) + a*JVS( 493 )
  W( 85 ) = W( 85 ) + a*JVS( 494 )
  W( 87 ) = W( 87 ) + a*JVS( 495 )
  W( 88 ) = W( 88 ) + a*JVS( 496 )
  W( 89 ) = W( 89 ) + a*JVS( 497 )
  a = -W( 68 ) / JVS(          499  )
  W( 68 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 500 )
  W( 82 ) = W( 82 ) + a*JVS( 501 )
  W( 83 ) = W( 83 ) + a*JVS( 502 )
  W( 84 ) = W( 84 ) + a*JVS( 503 )
  W( 85 ) = W( 85 ) + a*JVS( 504 )
  W( 87 ) = W( 87 ) + a*JVS( 505 )
  W( 88 ) = W( 88 ) + a*JVS( 506 )
  W( 89 ) = W( 89 ) + a*JVS( 507 )
  a = -W( 69 ) / JVS(          510  )
  W( 69 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 511 )
  W( 79 ) = W( 79 ) + a*JVS( 512 )
  W( 82 ) = W( 82 ) + a*JVS( 513 )
  W( 83 ) = W( 83 ) + a*JVS( 514 )
  W( 84 ) = W( 84 ) + a*JVS( 515 )
  W( 85 ) = W( 85 ) + a*JVS( 516 )
  W( 87 ) = W( 87 ) + a*JVS( 517 )
  W( 89 ) = W( 89 ) + a*JVS( 518 )
  a = -W( 70 ) / JVS(          528  )
  W( 70 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 529 )
  W( 77 ) = W( 77 ) + a*JVS( 530 )
  W( 78 ) = W( 78 ) + a*JVS( 531 )
  W( 79 ) = W( 79 ) + a*JVS( 532 )
  W( 80 ) = W( 80 ) + a*JVS( 533 )
  W( 82 ) = W( 82 ) + a*JVS( 534 )
  W( 83 ) = W( 83 ) + a*JVS( 535 )
  W( 84 ) = W( 84 ) + a*JVS( 536 )
  W( 85 ) = W( 85 ) + a*JVS( 537 )
  W( 87 ) = W( 87 ) + a*JVS( 538 )
  W( 88 ) = W( 88 ) + a*JVS( 539 )
  W( 89 ) = W( 89 ) + a*JVS( 540 )
  a = -W( 71 ) / JVS(          543  )
  W( 71 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 544 )
  W( 82 ) = W( 82 ) + a*JVS( 545 )
  W( 83 ) = W( 83 ) + a*JVS( 546 )
  W( 84 ) = W( 84 ) + a*JVS( 547 )
  W( 85 ) = W( 85 ) + a*JVS( 548 )
  W( 87 ) = W( 87 ) + a*JVS( 549 )
  W( 88 ) = W( 88 ) + a*JVS( 550 )
  W( 89 ) = W( 89 ) + a*JVS( 551 )
  a = -W( 72 ) / JVS(          562  )
  W( 72 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 563 )
  W( 79 ) = W( 79 ) + a*JVS( 564 )
  W( 80 ) = W( 80 ) + a*JVS( 565 )
  W( 81 ) = W( 81 ) + a*JVS( 566 )
  W( 82 ) = W( 82 ) + a*JVS( 567 )
  W( 83 ) = W( 83 ) + a*JVS( 568 )
  W( 84 ) = W( 84 ) + a*JVS( 569 )
  W( 85 ) = W( 85 ) + a*JVS( 570 )
  W( 86 ) = W( 86 ) + a*JVS( 571 )
  W( 87 ) = W( 87 ) + a*JVS( 572 )
  W( 88 ) = W( 88 ) + a*JVS( 573 )
  W( 89 ) = W( 89 ) + a*JVS( 574 )
  a = -W( 73 ) / JVS(          577  )
  W( 73 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 578 )
  W( 82 ) = W( 82 ) + a*JVS( 579 )
  W( 83 ) = W( 83 ) + a*JVS( 580 )
  W( 84 ) = W( 84 ) + a*JVS( 581 )
  W( 85 ) = W( 85 ) + a*JVS( 582 )
  W( 87 ) = W( 87 ) + a*JVS( 583 )
  W( 89 ) = W( 89 ) + a*JVS( 584 )
  a = -W( 74 ) / JVS(          600  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 601 )
  W( 76 ) = W( 76 ) + a*JVS( 602 )
  W( 77 ) = W( 77 ) + a*JVS( 603 )
  W( 78 ) = W( 78 ) + a*JVS( 604 )
  W( 79 ) = W( 79 ) + a*JVS( 605 )
  W( 80 ) = W( 80 ) + a*JVS( 606 )
  W( 81 ) = W( 81 ) + a*JVS( 607 )
  W( 82 ) = W( 82 ) + a*JVS( 608 )
  W( 83 ) = W( 83 ) + a*JVS( 609 )
  W( 84 ) = W( 84 ) + a*JVS( 610 )
  W( 85 ) = W( 85 ) + a*JVS( 611 )
  W( 86 ) = W( 86 ) + a*JVS( 612 )
  W( 87 ) = W( 87 ) + a*JVS( 613 )
  W( 88 ) = W( 88 ) + a*JVS( 614 )
  W( 89 ) = W( 89 ) + a*JVS( 615 )
  W( 90 ) = W( 90 ) + a*JVS( 616 )
  a = -W( 75 ) / JVS(          620  )
  W( 75 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 621 )
  W( 82 ) = W( 82 ) + a*JVS( 622 )
  W( 83 ) = W( 83 ) + a*JVS( 623 )
  W( 84 ) = W( 84 ) + a*JVS( 624 )
  W( 85 ) = W( 85 ) + a*JVS( 625 )
  W( 87 ) = W( 87 ) + a*JVS( 626 )
  W( 89 ) = W( 89 ) + a*JVS( 627 )
  a = -W( 76 ) / JVS(          636  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 637 )
  W( 79 ) = W( 79 ) + a*JVS( 638 )
  W( 82 ) = W( 82 ) + a*JVS( 639 )
  W( 83 ) = W( 83 ) + a*JVS( 640 )
  W( 84 ) = W( 84 ) + a*JVS( 641 )
  W( 85 ) = W( 85 ) + a*JVS( 642 )
  W( 87 ) = W( 87 ) + a*JVS( 643 )
  W( 89 ) = W( 89 ) + a*JVS( 644 )
  a = -W( 77 ) / JVS(          654  )
  W( 77 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 655 )
  W( 82 ) = W( 82 ) + a*JVS( 656 )
  W( 83 ) = W( 83 ) + a*JVS( 657 )
  W( 84 ) = W( 84 ) + a*JVS( 658 )
  W( 85 ) = W( 85 ) + a*JVS( 659 )
  W( 87 ) = W( 87 ) + a*JVS( 660 )
  W( 89 ) = W( 89 ) + a*JVS( 661 )
  a = -W( 78 ) / JVS(          664  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 665 )
  W( 82 ) = W( 82 ) + a*JVS( 666 )
  W( 83 ) = W( 83 ) + a*JVS( 667 )
  W( 84 ) = W( 84 ) + a*JVS( 668 )
  W( 85 ) = W( 85 ) + a*JVS( 669 )
  W( 87 ) = W( 87 ) + a*JVS( 670 )
  W( 88 ) = W( 88 ) + a*JVS( 671 )
  W( 89 ) = W( 89 ) + a*JVS( 672 )
  a = -W( 79 ) / JVS(          693  )
  W( 79 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 83 ) = W( 83 ) + a*JVS( 695 )
  W( 84 ) = W( 84 ) + a*JVS( 696 )
  W( 85 ) = W( 85 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  a = -W( 80 ) / JVS(          707  )
  W( 80 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 708 )
  W( 83 ) = W( 83 ) + a*JVS( 709 )
  W( 84 ) = W( 84 ) + a*JVS( 710 )
  W( 85 ) = W( 85 ) + a*JVS( 711 )
  W( 86 ) = W( 86 ) + a*JVS( 712 )
  W( 87 ) = W( 87 ) + a*JVS( 713 )
  W( 88 ) = W( 88 ) + a*JVS( 714 )
  W( 89 ) = W( 89 ) + a*JVS( 715 )
  a = -W( 81 ) / JVS(          733  )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 734 )
  W( 83 ) = W( 83 ) + a*JVS( 735 )
  W( 84 ) = W( 84 ) + a*JVS( 736 )
  W( 85 ) = W( 85 ) + a*JVS( 737 )
  W( 86 ) = W( 86 ) + a*JVS( 738 )
  W( 87 ) = W( 87 ) + a*JVS( 739 )
  W( 88 ) = W( 88 ) + a*JVS( 740 )
  W( 89 ) = W( 89 ) + a*JVS( 741 )
  a = -W( 82 ) / JVS(          775  )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 776 )
  W( 84 ) = W( 84 ) + a*JVS( 777 )
  W( 85 ) = W( 85 ) + a*JVS( 778 )
  W( 86 ) = W( 86 ) + a*JVS( 779 )
  W( 87 ) = W( 87 ) + a*JVS( 780 )
  W( 88 ) = W( 88 ) + a*JVS( 781 )
  W( 89 ) = W( 89 ) + a*JVS( 782 )
  W( 90 ) = W( 90 ) + a*JVS( 783 )
  a = -W( 83 ) / JVS(          813  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 814 )
  W( 85 ) = W( 85 ) + a*JVS( 815 )
  W( 86 ) = W( 86 ) + a*JVS( 816 )
  W( 87 ) = W( 87 ) + a*JVS( 817 )
  W( 88 ) = W( 88 ) + a*JVS( 818 )
  W( 89 ) = W( 89 ) + a*JVS( 819 )
  W( 90 ) = W( 90 ) + a*JVS( 820 )
  JVS( 821) = W( 13 )
  JVS( 822) = W( 14 )
  JVS( 823) = W( 16 )
  JVS( 824) = W( 17 )
  JVS( 825) = W( 18 )
  JVS( 826) = W( 20 )
  JVS( 827) = W( 21 )
  JVS( 828) = W( 22 )
  JVS( 829) = W( 23 )
  JVS( 830) = W( 24 )
  JVS( 831) = W( 25 )
  JVS( 832) = W( 26 )
  JVS( 833) = W( 27 )
  JVS( 834) = W( 28 )
  JVS( 835) = W( 29 )
  JVS( 836) = W( 30 )
  JVS( 837) = W( 31 )
  JVS( 838) = W( 32 )
  JVS( 839) = W( 35 )
  JVS( 840) = W( 36 )
  JVS( 841) = W( 37 )
  JVS( 842) = W( 38 )
  JVS( 843) = W( 39 )
  JVS( 844) = W( 40 )
  JVS( 845) = W( 41 )
  JVS( 846) = W( 42 )
  JVS( 847) = W( 43 )
  JVS( 848) = W( 44 )
  JVS( 849) = W( 45 )
  JVS( 850) = W( 46 )
  JVS( 851) = W( 47 )
  JVS( 852) = W( 48 )
  JVS( 853) = W( 49 )
  JVS( 854) = W( 50 )
  JVS( 855) = W( 51 )
  JVS( 856) = W( 52 )
  JVS( 857) = W( 53 )
  JVS( 858) = W( 54 )
  JVS( 859) = W( 57 )
  JVS( 860) = W( 59 )
  JVS( 861) = W( 60 )
  JVS( 862) = W( 61 )
  JVS( 863) = W( 62 )
  JVS( 864) = W( 63 )
  JVS( 865) = W( 64 )
  JVS( 866) = W( 65 )
  JVS( 867) = W( 66 )
  JVS( 868) = W( 67 )
  JVS( 869) = W( 68 )
  JVS( 870) = W( 69 )
  JVS( 871) = W( 70 )
  JVS( 872) = W( 71 )
  JVS( 873) = W( 72 )
  JVS( 874) = W( 73 )
  JVS( 875) = W( 74 )
  JVS( 876) = W( 75 )
  JVS( 877) = W( 76 )
  JVS( 878) = W( 77 )
  JVS( 879) = W( 78 )
  JVS( 880) = W( 79 )
  JVS( 881) = W( 80 )
  JVS( 882) = W( 81 )
  JVS( 883) = W( 82 )
  JVS( 884) = W( 83 )
  JVS( 885) = W( 84 )
  JVS( 886) = W( 85 )
  JVS( 887) = W( 86 )
  JVS( 888) = W( 87 )
  JVS( 889) = W( 88 )
  JVS( 890) = W( 89 )
  JVS( 891) = W( 90 )
  IF ( ABS(  JVS( 950 )) < TINY(a) ) THEN
         IER = 85                                      
         RETURN
  END IF
   W( 14 ) = JVS( 892 )
   W( 16 ) = JVS( 893 )
   W( 17 ) = JVS( 894 )
   W( 18 ) = JVS( 895 )
   W( 20 ) = JVS( 896 )
   W( 22 ) = JVS( 897 )
   W( 23 ) = JVS( 898 )
   W( 27 ) = JVS( 899 )
   W( 29 ) = JVS( 900 )
   W( 30 ) = JVS( 901 )
   W( 32 ) = JVS( 902 )
   W( 34 ) = JVS( 903 )
   W( 35 ) = JVS( 904 )
   W( 36 ) = JVS( 905 )
   W( 37 ) = JVS( 906 )
   W( 38 ) = JVS( 907 )
   W( 39 ) = JVS( 908 )
   W( 40 ) = JVS( 909 )
   W( 41 ) = JVS( 910 )
   W( 45 ) = JVS( 911 )
   W( 46 ) = JVS( 912 )
   W( 47 ) = JVS( 913 )
   W( 48 ) = JVS( 914 )
   W( 49 ) = JVS( 915 )
   W( 50 ) = JVS( 916 )
   W( 51 ) = JVS( 917 )
   W( 52 ) = JVS( 918 )
   W( 53 ) = JVS( 919 )
   W( 54 ) = JVS( 920 )
   W( 55 ) = JVS( 921 )
   W( 56 ) = JVS( 922 )
   W( 57 ) = JVS( 923 )
   W( 58 ) = JVS( 924 )
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
   W( 74 ) = JVS( 939 )
   W( 75 ) = JVS( 940 )
   W( 76 ) = JVS( 941 )
   W( 77 ) = JVS( 942 )
   W( 78 ) = JVS( 943 )
   W( 79 ) = JVS( 944 )
   W( 80 ) = JVS( 945 )
   W( 81 ) = JVS( 946 )
   W( 82 ) = JVS( 947 )
   W( 83 ) = JVS( 948 )
   W( 84 ) = JVS( 949 )
   W( 85 ) = JVS( 950 )
   W( 86 ) = JVS( 951 )
   W( 87 ) = JVS( 952 )
   W( 88 ) = JVS( 953 )
   W( 89 ) = JVS( 954 )
   W( 90 ) = JVS( 955 )
  a = -W( 14 ) / JVS(           67  )
  W( 14 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 68 )
  a = -W( 16 ) / JVS(           73  )
  W( 16 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 74 )
  a = -W( 17 ) / JVS(           75  )
  W( 17 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 76 )
  a = -W( 18 ) / JVS(           77  )
  W( 18 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 78 )
  a = -W( 20 ) / JVS(           82  )
  W( 20 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 83 )
  a = -W( 22 ) / JVS(           88  )
  W( 22 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 89 )
  a = -W( 23 ) / JVS(           90  )
  W( 23 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 91 )
  W( 84 ) = W( 84 ) + a*JVS( 92 )
  a = -W( 27 ) / JVS(          108  )
  W( 27 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 109 )
  W( 84 ) = W( 84 ) + a*JVS( 110 )
  W( 85 ) = W( 85 ) + a*JVS( 111 )
  a = -W( 29 ) / JVS(          116  )
  W( 29 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 117 )
  W( 85 ) = W( 85 ) + a*JVS( 118 )
  W( 88 ) = W( 88 ) + a*JVS( 119 )
  a = -W( 30 ) / JVS(          121  )
  W( 30 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 122 )
  W( 52 ) = W( 52 ) + a*JVS( 123 )
  W( 84 ) = W( 84 ) + a*JVS( 124 )
  W( 87 ) = W( 87 ) + a*JVS( 125 )
  W( 88 ) = W( 88 ) + a*JVS( 126 )
  a = -W( 32 ) / JVS(          135  )
  W( 32 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 136 )
  W( 81 ) = W( 81 ) + a*JVS( 137 )
  W( 83 ) = W( 83 ) + a*JVS( 138 )
  W( 84 ) = W( 84 ) + a*JVS( 139 )
  a = -W( 34 ) / JVS(          147  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 148 )
  W( 84 ) = W( 84 ) + a*JVS( 149 )
  W( 85 ) = W( 85 ) + a*JVS( 150 )
  W( 88 ) = W( 88 ) + a*JVS( 151 )
  W( 89 ) = W( 89 ) + a*JVS( 152 )
  a = -W( 35 ) / JVS(          153  )
  W( 35 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 154 )
  W( 47 ) = W( 47 ) + a*JVS( 155 )
  W( 48 ) = W( 48 ) + a*JVS( 156 )
  W( 49 ) = W( 49 ) + a*JVS( 157 )
  W( 51 ) = W( 51 ) + a*JVS( 158 )
  W( 73 ) = W( 73 ) + a*JVS( 159 )
  W( 75 ) = W( 75 ) + a*JVS( 160 )
  W( 79 ) = W( 79 ) + a*JVS( 161 )
  W( 84 ) = W( 84 ) + a*JVS( 162 )
  W( 85 ) = W( 85 ) + a*JVS( 163 )
  a = -W( 36 ) / JVS(          165  )
  W( 36 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 166 )
  W( 84 ) = W( 84 ) + a*JVS( 167 )
  W( 88 ) = W( 88 ) + a*JVS( 168 )
  a = -W( 37 ) / JVS(          170  )
  W( 37 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 171 )
  W( 84 ) = W( 84 ) + a*JVS( 172 )
  W( 88 ) = W( 88 ) + a*JVS( 173 )
  a = -W( 38 ) / JVS(          174  )
  W( 38 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 175 )
  W( 84 ) = W( 84 ) + a*JVS( 176 )
  W( 89 ) = W( 89 ) + a*JVS( 177 )
  a = -W( 39 ) / JVS(          178  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 179 )
  W( 79 ) = W( 79 ) + a*JVS( 180 )
  W( 84 ) = W( 84 ) + a*JVS( 181 )
  W( 88 ) = W( 88 ) + a*JVS( 182 )
  a = -W( 40 ) / JVS(          183  )
  W( 40 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 184 )
  W( 79 ) = W( 79 ) + a*JVS( 185 )
  W( 82 ) = W( 82 ) + a*JVS( 186 )
  W( 84 ) = W( 84 ) + a*JVS( 187 )
  W( 85 ) = W( 85 ) + a*JVS( 188 )
  a = -W( 41 ) / JVS(          189  )
  W( 41 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 190 )
  W( 60 ) = W( 60 ) + a*JVS( 191 )
  W( 84 ) = W( 84 ) + a*JVS( 192 )
  W( 87 ) = W( 87 ) + a*JVS( 193 )
  W( 89 ) = W( 89 ) + a*JVS( 194 )
  a = -W( 45 ) / JVS(          217  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 218 )
  W( 84 ) = W( 84 ) + a*JVS( 219 )
  W( 89 ) = W( 89 ) + a*JVS( 220 )
  a = -W( 46 ) / JVS(          226  )
  W( 46 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 227 )
  W( 48 ) = W( 48 ) + a*JVS( 228 )
  W( 49 ) = W( 49 ) + a*JVS( 229 )
  W( 51 ) = W( 51 ) + a*JVS( 230 )
  W( 52 ) = W( 52 ) + a*JVS( 231 )
  W( 53 ) = W( 53 ) + a*JVS( 232 )
  W( 54 ) = W( 54 ) + a*JVS( 233 )
  W( 57 ) = W( 57 ) + a*JVS( 234 )
  W( 60 ) = W( 60 ) + a*JVS( 235 )
  W( 62 ) = W( 62 ) + a*JVS( 236 )
  W( 66 ) = W( 66 ) + a*JVS( 237 )
  W( 71 ) = W( 71 ) + a*JVS( 238 )
  W( 73 ) = W( 73 ) + a*JVS( 239 )
  W( 74 ) = W( 74 ) + a*JVS( 240 )
  W( 75 ) = W( 75 ) + a*JVS( 241 )
  W( 79 ) = W( 79 ) + a*JVS( 242 )
  W( 84 ) = W( 84 ) + a*JVS( 243 )
  W( 87 ) = W( 87 ) + a*JVS( 244 )
  W( 89 ) = W( 89 ) + a*JVS( 245 )
  a = -W( 47 ) / JVS(          246  )
  W( 47 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 247 )
  W( 84 ) = W( 84 ) + a*JVS( 248 )
  W( 89 ) = W( 89 ) + a*JVS( 249 )
  a = -W( 48 ) / JVS(          250  )
  W( 48 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 251 )
  W( 84 ) = W( 84 ) + a*JVS( 252 )
  W( 89 ) = W( 89 ) + a*JVS( 253 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 255 )
  W( 84 ) = W( 84 ) + a*JVS( 256 )
  W( 89 ) = W( 89 ) + a*JVS( 257 )
  a = -W( 50 ) / JVS(          262  )
  W( 50 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 263 )
  W( 84 ) = W( 84 ) + a*JVS( 264 )
  W( 85 ) = W( 85 ) + a*JVS( 265 )
  W( 88 ) = W( 88 ) + a*JVS( 266 )
  W( 89 ) = W( 89 ) + a*JVS( 267 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  a = -W( 52 ) / JVS(          275  )
  W( 52 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 276 )
  W( 79 ) = W( 79 ) + a*JVS( 277 )
  W( 84 ) = W( 84 ) + a*JVS( 278 )
  W( 85 ) = W( 85 ) + a*JVS( 279 )
  W( 87 ) = W( 87 ) + a*JVS( 280 )
  W( 88 ) = W( 88 ) + a*JVS( 281 )
  W( 89 ) = W( 89 ) + a*JVS( 282 )
  a = -W( 53 ) / JVS(          284  )
  W( 53 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 285 )
  W( 67 ) = W( 67 ) + a*JVS( 286 )
  W( 68 ) = W( 68 ) + a*JVS( 287 )
  W( 71 ) = W( 71 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  W( 80 ) = W( 80 ) + a*JVS( 290 )
  W( 82 ) = W( 82 ) + a*JVS( 291 )
  W( 83 ) = W( 83 ) + a*JVS( 292 )
  W( 84 ) = W( 84 ) + a*JVS( 293 )
  W( 87 ) = W( 87 ) + a*JVS( 294 )
  W( 89 ) = W( 89 ) + a*JVS( 295 )
  a = -W( 54 ) / JVS(          296  )
  W( 54 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 297 )
  W( 79 ) = W( 79 ) + a*JVS( 298 )
  W( 84 ) = W( 84 ) + a*JVS( 299 )
  W( 88 ) = W( 88 ) + a*JVS( 300 )
  W( 89 ) = W( 89 ) + a*JVS( 301 )
  a = -W( 55 ) / JVS(          303  )
  W( 55 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 304 )
  W( 82 ) = W( 82 ) + a*JVS( 305 )
  W( 83 ) = W( 83 ) + a*JVS( 306 )
  W( 84 ) = W( 84 ) + a*JVS( 307 )
  W( 85 ) = W( 85 ) + a*JVS( 308 )
  W( 87 ) = W( 87 ) + a*JVS( 309 )
  W( 89 ) = W( 89 ) + a*JVS( 310 )
  a = -W( 56 ) / JVS(          311  )
  W( 56 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 312 )
  W( 82 ) = W( 82 ) + a*JVS( 313 )
  W( 83 ) = W( 83 ) + a*JVS( 314 )
  W( 84 ) = W( 84 ) + a*JVS( 315 )
  W( 85 ) = W( 85 ) + a*JVS( 316 )
  W( 87 ) = W( 87 ) + a*JVS( 317 )
  W( 89 ) = W( 89 ) + a*JVS( 318 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 325 )
  W( 63 ) = W( 63 ) + a*JVS( 326 )
  W( 79 ) = W( 79 ) + a*JVS( 327 )
  W( 82 ) = W( 82 ) + a*JVS( 328 )
  W( 83 ) = W( 83 ) + a*JVS( 329 )
  W( 84 ) = W( 84 ) + a*JVS( 330 )
  W( 85 ) = W( 85 ) + a*JVS( 331 )
  W( 87 ) = W( 87 ) + a*JVS( 332 )
  W( 89 ) = W( 89 ) + a*JVS( 333 )
  a = -W( 58 ) / JVS(          335  )
  W( 58 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 336 )
  W( 82 ) = W( 82 ) + a*JVS( 337 )
  W( 83 ) = W( 83 ) + a*JVS( 338 )
  W( 84 ) = W( 84 ) + a*JVS( 339 )
  W( 85 ) = W( 85 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 89 ) = W( 89 ) + a*JVS( 342 )
  a = -W( 60 ) / JVS(          372  )
  W( 60 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 82 ) = W( 82 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 85 ) = W( 85 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  a = -W( 61 ) / JVS(          382  )
  W( 61 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 82 ) = W( 82 ) + a*JVS( 384 )
  W( 83 ) = W( 83 ) + a*JVS( 385 )
  W( 84 ) = W( 84 ) + a*JVS( 386 )
  W( 85 ) = W( 85 ) + a*JVS( 387 )
  W( 87 ) = W( 87 ) + a*JVS( 388 )
  W( 88 ) = W( 88 ) + a*JVS( 389 )
  W( 89 ) = W( 89 ) + a*JVS( 390 )
  a = -W( 62 ) / JVS(          400  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 71 ) = W( 71 ) + a*JVS( 404 )
  W( 78 ) = W( 78 ) + a*JVS( 405 )
  W( 79 ) = W( 79 ) + a*JVS( 406 )
  W( 80 ) = W( 80 ) + a*JVS( 407 )
  W( 81 ) = W( 81 ) + a*JVS( 408 )
  W( 82 ) = W( 82 ) + a*JVS( 409 )
  W( 83 ) = W( 83 ) + a*JVS( 410 )
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 85 ) = W( 85 ) + a*JVS( 412 )
  W( 87 ) = W( 87 ) + a*JVS( 413 )
  W( 88 ) = W( 88 ) + a*JVS( 414 )
  W( 89 ) = W( 89 ) + a*JVS( 415 )
  a = -W( 63 ) / JVS(          417  )
  W( 63 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 418 )
  W( 82 ) = W( 82 ) + a*JVS( 419 )
  W( 83 ) = W( 83 ) + a*JVS( 420 )
  W( 84 ) = W( 84 ) + a*JVS( 421 )
  W( 85 ) = W( 85 ) + a*JVS( 422 )
  W( 87 ) = W( 87 ) + a*JVS( 423 )
  W( 89 ) = W( 89 ) + a*JVS( 424 )
  a = -W( 64 ) / JVS(          426  )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 85 ) = W( 85 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 89 ) = W( 89 ) + a*JVS( 432 )
  a = -W( 65 ) / JVS(          434  )
  W( 65 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 435 )
  W( 83 ) = W( 83 ) + a*JVS( 436 )
  W( 84 ) = W( 84 ) + a*JVS( 437 )
  W( 85 ) = W( 85 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 89 ) = W( 89 ) + a*JVS( 440 )
  a = -W( 66 ) / JVS(          466  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 467 )
  W( 68 ) = W( 68 ) + a*JVS( 468 )
  W( 69 ) = W( 69 ) + a*JVS( 469 )
  W( 71 ) = W( 71 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  W( 75 ) = W( 75 ) + a*JVS( 473 )
  W( 76 ) = W( 76 ) + a*JVS( 474 )
  W( 77 ) = W( 77 ) + a*JVS( 475 )
  W( 78 ) = W( 78 ) + a*JVS( 476 )
  W( 79 ) = W( 79 ) + a*JVS( 477 )
  W( 80 ) = W( 80 ) + a*JVS( 478 )
  W( 81 ) = W( 81 ) + a*JVS( 479 )
  W( 82 ) = W( 82 ) + a*JVS( 480 )
  W( 83 ) = W( 83 ) + a*JVS( 481 )
  W( 84 ) = W( 84 ) + a*JVS( 482 )
  W( 85 ) = W( 85 ) + a*JVS( 483 )
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 89 ) = W( 89 ) + a*JVS( 486 )
  W( 90 ) = W( 90 ) + a*JVS( 487 )
  a = -W( 67 ) / JVS(          489  )
  W( 67 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 490 )
  W( 82 ) = W( 82 ) + a*JVS( 491 )
  W( 83 ) = W( 83 ) + a*JVS( 492 )
  W( 84 ) = W( 84 ) + a*JVS( 493 )
  W( 85 ) = W( 85 ) + a*JVS( 494 )
  W( 87 ) = W( 87 ) + a*JVS( 495 )
  W( 88 ) = W( 88 ) + a*JVS( 496 )
  W( 89 ) = W( 89 ) + a*JVS( 497 )
  a = -W( 68 ) / JVS(          499  )
  W( 68 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 500 )
  W( 82 ) = W( 82 ) + a*JVS( 501 )
  W( 83 ) = W( 83 ) + a*JVS( 502 )
  W( 84 ) = W( 84 ) + a*JVS( 503 )
  W( 85 ) = W( 85 ) + a*JVS( 504 )
  W( 87 ) = W( 87 ) + a*JVS( 505 )
  W( 88 ) = W( 88 ) + a*JVS( 506 )
  W( 89 ) = W( 89 ) + a*JVS( 507 )
  a = -W( 69 ) / JVS(          510  )
  W( 69 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 511 )
  W( 79 ) = W( 79 ) + a*JVS( 512 )
  W( 82 ) = W( 82 ) + a*JVS( 513 )
  W( 83 ) = W( 83 ) + a*JVS( 514 )
  W( 84 ) = W( 84 ) + a*JVS( 515 )
  W( 85 ) = W( 85 ) + a*JVS( 516 )
  W( 87 ) = W( 87 ) + a*JVS( 517 )
  W( 89 ) = W( 89 ) + a*JVS( 518 )
  a = -W( 70 ) / JVS(          528  )
  W( 70 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 529 )
  W( 77 ) = W( 77 ) + a*JVS( 530 )
  W( 78 ) = W( 78 ) + a*JVS( 531 )
  W( 79 ) = W( 79 ) + a*JVS( 532 )
  W( 80 ) = W( 80 ) + a*JVS( 533 )
  W( 82 ) = W( 82 ) + a*JVS( 534 )
  W( 83 ) = W( 83 ) + a*JVS( 535 )
  W( 84 ) = W( 84 ) + a*JVS( 536 )
  W( 85 ) = W( 85 ) + a*JVS( 537 )
  W( 87 ) = W( 87 ) + a*JVS( 538 )
  W( 88 ) = W( 88 ) + a*JVS( 539 )
  W( 89 ) = W( 89 ) + a*JVS( 540 )
  a = -W( 71 ) / JVS(          543  )
  W( 71 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 544 )
  W( 82 ) = W( 82 ) + a*JVS( 545 )
  W( 83 ) = W( 83 ) + a*JVS( 546 )
  W( 84 ) = W( 84 ) + a*JVS( 547 )
  W( 85 ) = W( 85 ) + a*JVS( 548 )
  W( 87 ) = W( 87 ) + a*JVS( 549 )
  W( 88 ) = W( 88 ) + a*JVS( 550 )
  W( 89 ) = W( 89 ) + a*JVS( 551 )
  a = -W( 72 ) / JVS(          562  )
  W( 72 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 563 )
  W( 79 ) = W( 79 ) + a*JVS( 564 )
  W( 80 ) = W( 80 ) + a*JVS( 565 )
  W( 81 ) = W( 81 ) + a*JVS( 566 )
  W( 82 ) = W( 82 ) + a*JVS( 567 )
  W( 83 ) = W( 83 ) + a*JVS( 568 )
  W( 84 ) = W( 84 ) + a*JVS( 569 )
  W( 85 ) = W( 85 ) + a*JVS( 570 )
  W( 86 ) = W( 86 ) + a*JVS( 571 )
  W( 87 ) = W( 87 ) + a*JVS( 572 )
  W( 88 ) = W( 88 ) + a*JVS( 573 )
  W( 89 ) = W( 89 ) + a*JVS( 574 )
  a = -W( 73 ) / JVS(          577  )
  W( 73 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 578 )
  W( 82 ) = W( 82 ) + a*JVS( 579 )
  W( 83 ) = W( 83 ) + a*JVS( 580 )
  W( 84 ) = W( 84 ) + a*JVS( 581 )
  W( 85 ) = W( 85 ) + a*JVS( 582 )
  W( 87 ) = W( 87 ) + a*JVS( 583 )
  W( 89 ) = W( 89 ) + a*JVS( 584 )
  a = -W( 74 ) / JVS(          600  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 601 )
  W( 76 ) = W( 76 ) + a*JVS( 602 )
  W( 77 ) = W( 77 ) + a*JVS( 603 )
  W( 78 ) = W( 78 ) + a*JVS( 604 )
  W( 79 ) = W( 79 ) + a*JVS( 605 )
  W( 80 ) = W( 80 ) + a*JVS( 606 )
  W( 81 ) = W( 81 ) + a*JVS( 607 )
  W( 82 ) = W( 82 ) + a*JVS( 608 )
  W( 83 ) = W( 83 ) + a*JVS( 609 )
  W( 84 ) = W( 84 ) + a*JVS( 610 )
  W( 85 ) = W( 85 ) + a*JVS( 611 )
  W( 86 ) = W( 86 ) + a*JVS( 612 )
  W( 87 ) = W( 87 ) + a*JVS( 613 )
  W( 88 ) = W( 88 ) + a*JVS( 614 )
  W( 89 ) = W( 89 ) + a*JVS( 615 )
  W( 90 ) = W( 90 ) + a*JVS( 616 )
  a = -W( 75 ) / JVS(          620  )
  W( 75 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 621 )
  W( 82 ) = W( 82 ) + a*JVS( 622 )
  W( 83 ) = W( 83 ) + a*JVS( 623 )
  W( 84 ) = W( 84 ) + a*JVS( 624 )
  W( 85 ) = W( 85 ) + a*JVS( 625 )
  W( 87 ) = W( 87 ) + a*JVS( 626 )
  W( 89 ) = W( 89 ) + a*JVS( 627 )
  a = -W( 76 ) / JVS(          636  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 637 )
  W( 79 ) = W( 79 ) + a*JVS( 638 )
  W( 82 ) = W( 82 ) + a*JVS( 639 )
  W( 83 ) = W( 83 ) + a*JVS( 640 )
  W( 84 ) = W( 84 ) + a*JVS( 641 )
  W( 85 ) = W( 85 ) + a*JVS( 642 )
  W( 87 ) = W( 87 ) + a*JVS( 643 )
  W( 89 ) = W( 89 ) + a*JVS( 644 )
  a = -W( 77 ) / JVS(          654  )
  W( 77 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 655 )
  W( 82 ) = W( 82 ) + a*JVS( 656 )
  W( 83 ) = W( 83 ) + a*JVS( 657 )
  W( 84 ) = W( 84 ) + a*JVS( 658 )
  W( 85 ) = W( 85 ) + a*JVS( 659 )
  W( 87 ) = W( 87 ) + a*JVS( 660 )
  W( 89 ) = W( 89 ) + a*JVS( 661 )
  a = -W( 78 ) / JVS(          664  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 665 )
  W( 82 ) = W( 82 ) + a*JVS( 666 )
  W( 83 ) = W( 83 ) + a*JVS( 667 )
  W( 84 ) = W( 84 ) + a*JVS( 668 )
  W( 85 ) = W( 85 ) + a*JVS( 669 )
  W( 87 ) = W( 87 ) + a*JVS( 670 )
  W( 88 ) = W( 88 ) + a*JVS( 671 )
  W( 89 ) = W( 89 ) + a*JVS( 672 )
  a = -W( 79 ) / JVS(          693  )
  W( 79 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 83 ) = W( 83 ) + a*JVS( 695 )
  W( 84 ) = W( 84 ) + a*JVS( 696 )
  W( 85 ) = W( 85 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  a = -W( 80 ) / JVS(          707  )
  W( 80 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 708 )
  W( 83 ) = W( 83 ) + a*JVS( 709 )
  W( 84 ) = W( 84 ) + a*JVS( 710 )
  W( 85 ) = W( 85 ) + a*JVS( 711 )
  W( 86 ) = W( 86 ) + a*JVS( 712 )
  W( 87 ) = W( 87 ) + a*JVS( 713 )
  W( 88 ) = W( 88 ) + a*JVS( 714 )
  W( 89 ) = W( 89 ) + a*JVS( 715 )
  a = -W( 81 ) / JVS(          733  )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 734 )
  W( 83 ) = W( 83 ) + a*JVS( 735 )
  W( 84 ) = W( 84 ) + a*JVS( 736 )
  W( 85 ) = W( 85 ) + a*JVS( 737 )
  W( 86 ) = W( 86 ) + a*JVS( 738 )
  W( 87 ) = W( 87 ) + a*JVS( 739 )
  W( 88 ) = W( 88 ) + a*JVS( 740 )
  W( 89 ) = W( 89 ) + a*JVS( 741 )
  a = -W( 82 ) / JVS(          775  )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 776 )
  W( 84 ) = W( 84 ) + a*JVS( 777 )
  W( 85 ) = W( 85 ) + a*JVS( 778 )
  W( 86 ) = W( 86 ) + a*JVS( 779 )
  W( 87 ) = W( 87 ) + a*JVS( 780 )
  W( 88 ) = W( 88 ) + a*JVS( 781 )
  W( 89 ) = W( 89 ) + a*JVS( 782 )
  W( 90 ) = W( 90 ) + a*JVS( 783 )
  a = -W( 83 ) / JVS(          813  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 814 )
  W( 85 ) = W( 85 ) + a*JVS( 815 )
  W( 86 ) = W( 86 ) + a*JVS( 816 )
  W( 87 ) = W( 87 ) + a*JVS( 817 )
  W( 88 ) = W( 88 ) + a*JVS( 818 )
  W( 89 ) = W( 89 ) + a*JVS( 819 )
  W( 90 ) = W( 90 ) + a*JVS( 820 )
  a = -W( 84 ) / JVS(          885  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 886 )
  W( 86 ) = W( 86 ) + a*JVS( 887 )
  W( 87 ) = W( 87 ) + a*JVS( 888 )
  W( 88 ) = W( 88 ) + a*JVS( 889 )
  W( 89 ) = W( 89 ) + a*JVS( 890 )
  W( 90 ) = W( 90 ) + a*JVS( 891 )
  JVS( 892) = W( 14 )
  JVS( 893) = W( 16 )
  JVS( 894) = W( 17 )
  JVS( 895) = W( 18 )
  JVS( 896) = W( 20 )
  JVS( 897) = W( 22 )
  JVS( 898) = W( 23 )
  JVS( 899) = W( 27 )
  JVS( 900) = W( 29 )
  JVS( 901) = W( 30 )
  JVS( 902) = W( 32 )
  JVS( 903) = W( 34 )
  JVS( 904) = W( 35 )
  JVS( 905) = W( 36 )
  JVS( 906) = W( 37 )
  JVS( 907) = W( 38 )
  JVS( 908) = W( 39 )
  JVS( 909) = W( 40 )
  JVS( 910) = W( 41 )
  JVS( 911) = W( 45 )
  JVS( 912) = W( 46 )
  JVS( 913) = W( 47 )
  JVS( 914) = W( 48 )
  JVS( 915) = W( 49 )
  JVS( 916) = W( 50 )
  JVS( 917) = W( 51 )
  JVS( 918) = W( 52 )
  JVS( 919) = W( 53 )
  JVS( 920) = W( 54 )
  JVS( 921) = W( 55 )
  JVS( 922) = W( 56 )
  JVS( 923) = W( 57 )
  JVS( 924) = W( 58 )
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
  JVS( 939) = W( 74 )
  JVS( 940) = W( 75 )
  JVS( 941) = W( 76 )
  JVS( 942) = W( 77 )
  JVS( 943) = W( 78 )
  JVS( 944) = W( 79 )
  JVS( 945) = W( 80 )
  JVS( 946) = W( 81 )
  JVS( 947) = W( 82 )
  JVS( 948) = W( 83 )
  JVS( 949) = W( 84 )
  JVS( 950) = W( 85 )
  JVS( 951) = W( 86 )
  JVS( 952) = W( 87 )
  JVS( 953) = W( 88 )
  JVS( 954) = W( 89 )
  JVS( 955) = W( 90 )
  IF ( ABS(  JVS( 977 )) < TINY(a) ) THEN
         IER = 86                                      
         RETURN
  END IF
   W( 55 ) = JVS( 956 )
   W( 56 ) = JVS( 957 )
   W( 58 ) = JVS( 958 )
   W( 61 ) = JVS( 959 )
   W( 63 ) = JVS( 960 )
   W( 64 ) = JVS( 961 )
   W( 65 ) = JVS( 962 )
   W( 67 ) = JVS( 963 )
   W( 68 ) = JVS( 964 )
   W( 69 ) = JVS( 965 )
   W( 72 ) = JVS( 966 )
   W( 73 ) = JVS( 967 )
   W( 75 ) = JVS( 968 )
   W( 78 ) = JVS( 969 )
   W( 79 ) = JVS( 970 )
   W( 80 ) = JVS( 971 )
   W( 81 ) = JVS( 972 )
   W( 82 ) = JVS( 973 )
   W( 83 ) = JVS( 974 )
   W( 84 ) = JVS( 975 )
   W( 85 ) = JVS( 976 )
   W( 86 ) = JVS( 977 )
   W( 87 ) = JVS( 978 )
   W( 88 ) = JVS( 979 )
   W( 89 ) = JVS( 980 )
   W( 90 ) = JVS( 981 )
  a = -W( 55 ) / JVS(          303  )
  W( 55 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 304 )
  W( 82 ) = W( 82 ) + a*JVS( 305 )
  W( 83 ) = W( 83 ) + a*JVS( 306 )
  W( 84 ) = W( 84 ) + a*JVS( 307 )
  W( 85 ) = W( 85 ) + a*JVS( 308 )
  W( 87 ) = W( 87 ) + a*JVS( 309 )
  W( 89 ) = W( 89 ) + a*JVS( 310 )
  a = -W( 56 ) / JVS(          311  )
  W( 56 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 312 )
  W( 82 ) = W( 82 ) + a*JVS( 313 )
  W( 83 ) = W( 83 ) + a*JVS( 314 )
  W( 84 ) = W( 84 ) + a*JVS( 315 )
  W( 85 ) = W( 85 ) + a*JVS( 316 )
  W( 87 ) = W( 87 ) + a*JVS( 317 )
  W( 89 ) = W( 89 ) + a*JVS( 318 )
  a = -W( 58 ) / JVS(          335  )
  W( 58 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 336 )
  W( 82 ) = W( 82 ) + a*JVS( 337 )
  W( 83 ) = W( 83 ) + a*JVS( 338 )
  W( 84 ) = W( 84 ) + a*JVS( 339 )
  W( 85 ) = W( 85 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 89 ) = W( 89 ) + a*JVS( 342 )
  a = -W( 61 ) / JVS(          382  )
  W( 61 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 82 ) = W( 82 ) + a*JVS( 384 )
  W( 83 ) = W( 83 ) + a*JVS( 385 )
  W( 84 ) = W( 84 ) + a*JVS( 386 )
  W( 85 ) = W( 85 ) + a*JVS( 387 )
  W( 87 ) = W( 87 ) + a*JVS( 388 )
  W( 88 ) = W( 88 ) + a*JVS( 389 )
  W( 89 ) = W( 89 ) + a*JVS( 390 )
  a = -W( 63 ) / JVS(          417  )
  W( 63 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 418 )
  W( 82 ) = W( 82 ) + a*JVS( 419 )
  W( 83 ) = W( 83 ) + a*JVS( 420 )
  W( 84 ) = W( 84 ) + a*JVS( 421 )
  W( 85 ) = W( 85 ) + a*JVS( 422 )
  W( 87 ) = W( 87 ) + a*JVS( 423 )
  W( 89 ) = W( 89 ) + a*JVS( 424 )
  a = -W( 64 ) / JVS(          426  )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 85 ) = W( 85 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 89 ) = W( 89 ) + a*JVS( 432 )
  a = -W( 65 ) / JVS(          434  )
  W( 65 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 435 )
  W( 83 ) = W( 83 ) + a*JVS( 436 )
  W( 84 ) = W( 84 ) + a*JVS( 437 )
  W( 85 ) = W( 85 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 89 ) = W( 89 ) + a*JVS( 440 )
  a = -W( 67 ) / JVS(          489  )
  W( 67 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 490 )
  W( 82 ) = W( 82 ) + a*JVS( 491 )
  W( 83 ) = W( 83 ) + a*JVS( 492 )
  W( 84 ) = W( 84 ) + a*JVS( 493 )
  W( 85 ) = W( 85 ) + a*JVS( 494 )
  W( 87 ) = W( 87 ) + a*JVS( 495 )
  W( 88 ) = W( 88 ) + a*JVS( 496 )
  W( 89 ) = W( 89 ) + a*JVS( 497 )
  a = -W( 68 ) / JVS(          499  )
  W( 68 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 500 )
  W( 82 ) = W( 82 ) + a*JVS( 501 )
  W( 83 ) = W( 83 ) + a*JVS( 502 )
  W( 84 ) = W( 84 ) + a*JVS( 503 )
  W( 85 ) = W( 85 ) + a*JVS( 504 )
  W( 87 ) = W( 87 ) + a*JVS( 505 )
  W( 88 ) = W( 88 ) + a*JVS( 506 )
  W( 89 ) = W( 89 ) + a*JVS( 507 )
  a = -W( 69 ) / JVS(          510  )
  W( 69 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 511 )
  W( 79 ) = W( 79 ) + a*JVS( 512 )
  W( 82 ) = W( 82 ) + a*JVS( 513 )
  W( 83 ) = W( 83 ) + a*JVS( 514 )
  W( 84 ) = W( 84 ) + a*JVS( 515 )
  W( 85 ) = W( 85 ) + a*JVS( 516 )
  W( 87 ) = W( 87 ) + a*JVS( 517 )
  W( 89 ) = W( 89 ) + a*JVS( 518 )
  a = -W( 72 ) / JVS(          562  )
  W( 72 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 563 )
  W( 79 ) = W( 79 ) + a*JVS( 564 )
  W( 80 ) = W( 80 ) + a*JVS( 565 )
  W( 81 ) = W( 81 ) + a*JVS( 566 )
  W( 82 ) = W( 82 ) + a*JVS( 567 )
  W( 83 ) = W( 83 ) + a*JVS( 568 )
  W( 84 ) = W( 84 ) + a*JVS( 569 )
  W( 85 ) = W( 85 ) + a*JVS( 570 )
  W( 86 ) = W( 86 ) + a*JVS( 571 )
  W( 87 ) = W( 87 ) + a*JVS( 572 )
  W( 88 ) = W( 88 ) + a*JVS( 573 )
  W( 89 ) = W( 89 ) + a*JVS( 574 )
  a = -W( 73 ) / JVS(          577  )
  W( 73 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 578 )
  W( 82 ) = W( 82 ) + a*JVS( 579 )
  W( 83 ) = W( 83 ) + a*JVS( 580 )
  W( 84 ) = W( 84 ) + a*JVS( 581 )
  W( 85 ) = W( 85 ) + a*JVS( 582 )
  W( 87 ) = W( 87 ) + a*JVS( 583 )
  W( 89 ) = W( 89 ) + a*JVS( 584 )
  a = -W( 75 ) / JVS(          620  )
  W( 75 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 621 )
  W( 82 ) = W( 82 ) + a*JVS( 622 )
  W( 83 ) = W( 83 ) + a*JVS( 623 )
  W( 84 ) = W( 84 ) + a*JVS( 624 )
  W( 85 ) = W( 85 ) + a*JVS( 625 )
  W( 87 ) = W( 87 ) + a*JVS( 626 )
  W( 89 ) = W( 89 ) + a*JVS( 627 )
  a = -W( 78 ) / JVS(          664  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 665 )
  W( 82 ) = W( 82 ) + a*JVS( 666 )
  W( 83 ) = W( 83 ) + a*JVS( 667 )
  W( 84 ) = W( 84 ) + a*JVS( 668 )
  W( 85 ) = W( 85 ) + a*JVS( 669 )
  W( 87 ) = W( 87 ) + a*JVS( 670 )
  W( 88 ) = W( 88 ) + a*JVS( 671 )
  W( 89 ) = W( 89 ) + a*JVS( 672 )
  a = -W( 79 ) / JVS(          693  )
  W( 79 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 83 ) = W( 83 ) + a*JVS( 695 )
  W( 84 ) = W( 84 ) + a*JVS( 696 )
  W( 85 ) = W( 85 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  a = -W( 80 ) / JVS(          707  )
  W( 80 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 708 )
  W( 83 ) = W( 83 ) + a*JVS( 709 )
  W( 84 ) = W( 84 ) + a*JVS( 710 )
  W( 85 ) = W( 85 ) + a*JVS( 711 )
  W( 86 ) = W( 86 ) + a*JVS( 712 )
  W( 87 ) = W( 87 ) + a*JVS( 713 )
  W( 88 ) = W( 88 ) + a*JVS( 714 )
  W( 89 ) = W( 89 ) + a*JVS( 715 )
  a = -W( 81 ) / JVS(          733  )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 734 )
  W( 83 ) = W( 83 ) + a*JVS( 735 )
  W( 84 ) = W( 84 ) + a*JVS( 736 )
  W( 85 ) = W( 85 ) + a*JVS( 737 )
  W( 86 ) = W( 86 ) + a*JVS( 738 )
  W( 87 ) = W( 87 ) + a*JVS( 739 )
  W( 88 ) = W( 88 ) + a*JVS( 740 )
  W( 89 ) = W( 89 ) + a*JVS( 741 )
  a = -W( 82 ) / JVS(          775  )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 776 )
  W( 84 ) = W( 84 ) + a*JVS( 777 )
  W( 85 ) = W( 85 ) + a*JVS( 778 )
  W( 86 ) = W( 86 ) + a*JVS( 779 )
  W( 87 ) = W( 87 ) + a*JVS( 780 )
  W( 88 ) = W( 88 ) + a*JVS( 781 )
  W( 89 ) = W( 89 ) + a*JVS( 782 )
  W( 90 ) = W( 90 ) + a*JVS( 783 )
  a = -W( 83 ) / JVS(          813  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 814 )
  W( 85 ) = W( 85 ) + a*JVS( 815 )
  W( 86 ) = W( 86 ) + a*JVS( 816 )
  W( 87 ) = W( 87 ) + a*JVS( 817 )
  W( 88 ) = W( 88 ) + a*JVS( 818 )
  W( 89 ) = W( 89 ) + a*JVS( 819 )
  W( 90 ) = W( 90 ) + a*JVS( 820 )
  a = -W( 84 ) / JVS(          885  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 886 )
  W( 86 ) = W( 86 ) + a*JVS( 887 )
  W( 87 ) = W( 87 ) + a*JVS( 888 )
  W( 88 ) = W( 88 ) + a*JVS( 889 )
  W( 89 ) = W( 89 ) + a*JVS( 890 )
  W( 90 ) = W( 90 ) + a*JVS( 891 )
  a = -W( 85 ) / JVS(          950  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 951 )
  W( 87 ) = W( 87 ) + a*JVS( 952 )
  W( 88 ) = W( 88 ) + a*JVS( 953 )
  W( 89 ) = W( 89 ) + a*JVS( 954 )
  W( 90 ) = W( 90 ) + a*JVS( 955 )
  JVS( 956) = W( 55 )
  JVS( 957) = W( 56 )
  JVS( 958) = W( 58 )
  JVS( 959) = W( 61 )
  JVS( 960) = W( 63 )
  JVS( 961) = W( 64 )
  JVS( 962) = W( 65 )
  JVS( 963) = W( 67 )
  JVS( 964) = W( 68 )
  JVS( 965) = W( 69 )
  JVS( 966) = W( 72 )
  JVS( 967) = W( 73 )
  JVS( 968) = W( 75 )
  JVS( 969) = W( 78 )
  JVS( 970) = W( 79 )
  JVS( 971) = W( 80 )
  JVS( 972) = W( 81 )
  JVS( 973) = W( 82 )
  JVS( 974) = W( 83 )
  JVS( 975) = W( 84 )
  JVS( 976) = W( 85 )
  JVS( 977) = W( 86 )
  JVS( 978) = W( 87 )
  JVS( 979) = W( 88 )
  JVS( 980) = W( 89 )
  JVS( 981) = W( 90 )
  IF ( ABS(  JVS( 1016 )) < TINY(a) ) THEN
         IER = 87                                      
         RETURN
  END IF
   W( 31 ) = JVS( 982 )
   W( 33 ) = JVS( 983 )
   W( 36 ) = JVS( 984 )
   W( 37 ) = JVS( 985 )
   W( 39 ) = JVS( 986 )
   W( 50 ) = JVS( 987 )
   W( 51 ) = JVS( 988 )
   W( 52 ) = JVS( 989 )
   W( 55 ) = JVS( 990 )
   W( 56 ) = JVS( 991 )
   W( 57 ) = JVS( 992 )
   W( 58 ) = JVS( 993 )
   W( 60 ) = JVS( 994 )
   W( 61 ) = JVS( 995 )
   W( 63 ) = JVS( 996 )
   W( 64 ) = JVS( 997 )
   W( 65 ) = JVS( 998 )
   W( 67 ) = JVS( 999 )
   W( 68 ) = JVS( 1000 )
   W( 69 ) = JVS( 1001 )
   W( 72 ) = JVS( 1002 )
   W( 73 ) = JVS( 1003 )
   W( 75 ) = JVS( 1004 )
   W( 76 ) = JVS( 1005 )
   W( 77 ) = JVS( 1006 )
   W( 78 ) = JVS( 1007 )
   W( 79 ) = JVS( 1008 )
   W( 80 ) = JVS( 1009 )
   W( 81 ) = JVS( 1010 )
   W( 82 ) = JVS( 1011 )
   W( 83 ) = JVS( 1012 )
   W( 84 ) = JVS( 1013 )
   W( 85 ) = JVS( 1014 )
   W( 86 ) = JVS( 1015 )
   W( 87 ) = JVS( 1016 )
   W( 88 ) = JVS( 1017 )
   W( 89 ) = JVS( 1018 )
   W( 90 ) = JVS( 1019 )
  a = -W( 31 ) / JVS(          127  )
  W( 31 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 128 )
  W( 37 ) = W( 37 ) + a*JVS( 129 )
  W( 39 ) = W( 39 ) + a*JVS( 130 )
  W( 84 ) = W( 84 ) + a*JVS( 131 )
  W( 87 ) = W( 87 ) + a*JVS( 132 )
  W( 88 ) = W( 88 ) + a*JVS( 133 )
  a = -W( 33 ) / JVS(          141  )
  W( 33 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 142 )
  W( 79 ) = W( 79 ) + a*JVS( 143 )
  W( 87 ) = W( 87 ) + a*JVS( 144 )
  W( 88 ) = W( 88 ) + a*JVS( 145 )
  W( 89 ) = W( 89 ) + a*JVS( 146 )
  a = -W( 36 ) / JVS(          165  )
  W( 36 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 166 )
  W( 84 ) = W( 84 ) + a*JVS( 167 )
  W( 88 ) = W( 88 ) + a*JVS( 168 )
  a = -W( 37 ) / JVS(          170  )
  W( 37 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 171 )
  W( 84 ) = W( 84 ) + a*JVS( 172 )
  W( 88 ) = W( 88 ) + a*JVS( 173 )
  a = -W( 39 ) / JVS(          178  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 179 )
  W( 79 ) = W( 79 ) + a*JVS( 180 )
  W( 84 ) = W( 84 ) + a*JVS( 181 )
  W( 88 ) = W( 88 ) + a*JVS( 182 )
  a = -W( 50 ) / JVS(          262  )
  W( 50 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 263 )
  W( 84 ) = W( 84 ) + a*JVS( 264 )
  W( 85 ) = W( 85 ) + a*JVS( 265 )
  W( 88 ) = W( 88 ) + a*JVS( 266 )
  W( 89 ) = W( 89 ) + a*JVS( 267 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  a = -W( 52 ) / JVS(          275  )
  W( 52 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 276 )
  W( 79 ) = W( 79 ) + a*JVS( 277 )
  W( 84 ) = W( 84 ) + a*JVS( 278 )
  W( 85 ) = W( 85 ) + a*JVS( 279 )
  W( 87 ) = W( 87 ) + a*JVS( 280 )
  W( 88 ) = W( 88 ) + a*JVS( 281 )
  W( 89 ) = W( 89 ) + a*JVS( 282 )
  a = -W( 55 ) / JVS(          303  )
  W( 55 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 304 )
  W( 82 ) = W( 82 ) + a*JVS( 305 )
  W( 83 ) = W( 83 ) + a*JVS( 306 )
  W( 84 ) = W( 84 ) + a*JVS( 307 )
  W( 85 ) = W( 85 ) + a*JVS( 308 )
  W( 87 ) = W( 87 ) + a*JVS( 309 )
  W( 89 ) = W( 89 ) + a*JVS( 310 )
  a = -W( 56 ) / JVS(          311  )
  W( 56 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 312 )
  W( 82 ) = W( 82 ) + a*JVS( 313 )
  W( 83 ) = W( 83 ) + a*JVS( 314 )
  W( 84 ) = W( 84 ) + a*JVS( 315 )
  W( 85 ) = W( 85 ) + a*JVS( 316 )
  W( 87 ) = W( 87 ) + a*JVS( 317 )
  W( 89 ) = W( 89 ) + a*JVS( 318 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 325 )
  W( 63 ) = W( 63 ) + a*JVS( 326 )
  W( 79 ) = W( 79 ) + a*JVS( 327 )
  W( 82 ) = W( 82 ) + a*JVS( 328 )
  W( 83 ) = W( 83 ) + a*JVS( 329 )
  W( 84 ) = W( 84 ) + a*JVS( 330 )
  W( 85 ) = W( 85 ) + a*JVS( 331 )
  W( 87 ) = W( 87 ) + a*JVS( 332 )
  W( 89 ) = W( 89 ) + a*JVS( 333 )
  a = -W( 58 ) / JVS(          335  )
  W( 58 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 336 )
  W( 82 ) = W( 82 ) + a*JVS( 337 )
  W( 83 ) = W( 83 ) + a*JVS( 338 )
  W( 84 ) = W( 84 ) + a*JVS( 339 )
  W( 85 ) = W( 85 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 89 ) = W( 89 ) + a*JVS( 342 )
  a = -W( 60 ) / JVS(          372  )
  W( 60 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 82 ) = W( 82 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 85 ) = W( 85 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  a = -W( 61 ) / JVS(          382  )
  W( 61 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 82 ) = W( 82 ) + a*JVS( 384 )
  W( 83 ) = W( 83 ) + a*JVS( 385 )
  W( 84 ) = W( 84 ) + a*JVS( 386 )
  W( 85 ) = W( 85 ) + a*JVS( 387 )
  W( 87 ) = W( 87 ) + a*JVS( 388 )
  W( 88 ) = W( 88 ) + a*JVS( 389 )
  W( 89 ) = W( 89 ) + a*JVS( 390 )
  a = -W( 63 ) / JVS(          417  )
  W( 63 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 418 )
  W( 82 ) = W( 82 ) + a*JVS( 419 )
  W( 83 ) = W( 83 ) + a*JVS( 420 )
  W( 84 ) = W( 84 ) + a*JVS( 421 )
  W( 85 ) = W( 85 ) + a*JVS( 422 )
  W( 87 ) = W( 87 ) + a*JVS( 423 )
  W( 89 ) = W( 89 ) + a*JVS( 424 )
  a = -W( 64 ) / JVS(          426  )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 85 ) = W( 85 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 89 ) = W( 89 ) + a*JVS( 432 )
  a = -W( 65 ) / JVS(          434  )
  W( 65 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 435 )
  W( 83 ) = W( 83 ) + a*JVS( 436 )
  W( 84 ) = W( 84 ) + a*JVS( 437 )
  W( 85 ) = W( 85 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 89 ) = W( 89 ) + a*JVS( 440 )
  a = -W( 67 ) / JVS(          489  )
  W( 67 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 490 )
  W( 82 ) = W( 82 ) + a*JVS( 491 )
  W( 83 ) = W( 83 ) + a*JVS( 492 )
  W( 84 ) = W( 84 ) + a*JVS( 493 )
  W( 85 ) = W( 85 ) + a*JVS( 494 )
  W( 87 ) = W( 87 ) + a*JVS( 495 )
  W( 88 ) = W( 88 ) + a*JVS( 496 )
  W( 89 ) = W( 89 ) + a*JVS( 497 )
  a = -W( 68 ) / JVS(          499  )
  W( 68 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 500 )
  W( 82 ) = W( 82 ) + a*JVS( 501 )
  W( 83 ) = W( 83 ) + a*JVS( 502 )
  W( 84 ) = W( 84 ) + a*JVS( 503 )
  W( 85 ) = W( 85 ) + a*JVS( 504 )
  W( 87 ) = W( 87 ) + a*JVS( 505 )
  W( 88 ) = W( 88 ) + a*JVS( 506 )
  W( 89 ) = W( 89 ) + a*JVS( 507 )
  a = -W( 69 ) / JVS(          510  )
  W( 69 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 511 )
  W( 79 ) = W( 79 ) + a*JVS( 512 )
  W( 82 ) = W( 82 ) + a*JVS( 513 )
  W( 83 ) = W( 83 ) + a*JVS( 514 )
  W( 84 ) = W( 84 ) + a*JVS( 515 )
  W( 85 ) = W( 85 ) + a*JVS( 516 )
  W( 87 ) = W( 87 ) + a*JVS( 517 )
  W( 89 ) = W( 89 ) + a*JVS( 518 )
  a = -W( 72 ) / JVS(          562  )
  W( 72 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 563 )
  W( 79 ) = W( 79 ) + a*JVS( 564 )
  W( 80 ) = W( 80 ) + a*JVS( 565 )
  W( 81 ) = W( 81 ) + a*JVS( 566 )
  W( 82 ) = W( 82 ) + a*JVS( 567 )
  W( 83 ) = W( 83 ) + a*JVS( 568 )
  W( 84 ) = W( 84 ) + a*JVS( 569 )
  W( 85 ) = W( 85 ) + a*JVS( 570 )
  W( 86 ) = W( 86 ) + a*JVS( 571 )
  W( 87 ) = W( 87 ) + a*JVS( 572 )
  W( 88 ) = W( 88 ) + a*JVS( 573 )
  W( 89 ) = W( 89 ) + a*JVS( 574 )
  a = -W( 73 ) / JVS(          577  )
  W( 73 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 578 )
  W( 82 ) = W( 82 ) + a*JVS( 579 )
  W( 83 ) = W( 83 ) + a*JVS( 580 )
  W( 84 ) = W( 84 ) + a*JVS( 581 )
  W( 85 ) = W( 85 ) + a*JVS( 582 )
  W( 87 ) = W( 87 ) + a*JVS( 583 )
  W( 89 ) = W( 89 ) + a*JVS( 584 )
  a = -W( 75 ) / JVS(          620  )
  W( 75 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 621 )
  W( 82 ) = W( 82 ) + a*JVS( 622 )
  W( 83 ) = W( 83 ) + a*JVS( 623 )
  W( 84 ) = W( 84 ) + a*JVS( 624 )
  W( 85 ) = W( 85 ) + a*JVS( 625 )
  W( 87 ) = W( 87 ) + a*JVS( 626 )
  W( 89 ) = W( 89 ) + a*JVS( 627 )
  a = -W( 76 ) / JVS(          636  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 637 )
  W( 79 ) = W( 79 ) + a*JVS( 638 )
  W( 82 ) = W( 82 ) + a*JVS( 639 )
  W( 83 ) = W( 83 ) + a*JVS( 640 )
  W( 84 ) = W( 84 ) + a*JVS( 641 )
  W( 85 ) = W( 85 ) + a*JVS( 642 )
  W( 87 ) = W( 87 ) + a*JVS( 643 )
  W( 89 ) = W( 89 ) + a*JVS( 644 )
  a = -W( 77 ) / JVS(          654  )
  W( 77 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 655 )
  W( 82 ) = W( 82 ) + a*JVS( 656 )
  W( 83 ) = W( 83 ) + a*JVS( 657 )
  W( 84 ) = W( 84 ) + a*JVS( 658 )
  W( 85 ) = W( 85 ) + a*JVS( 659 )
  W( 87 ) = W( 87 ) + a*JVS( 660 )
  W( 89 ) = W( 89 ) + a*JVS( 661 )
  a = -W( 78 ) / JVS(          664  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 665 )
  W( 82 ) = W( 82 ) + a*JVS( 666 )
  W( 83 ) = W( 83 ) + a*JVS( 667 )
  W( 84 ) = W( 84 ) + a*JVS( 668 )
  W( 85 ) = W( 85 ) + a*JVS( 669 )
  W( 87 ) = W( 87 ) + a*JVS( 670 )
  W( 88 ) = W( 88 ) + a*JVS( 671 )
  W( 89 ) = W( 89 ) + a*JVS( 672 )
  a = -W( 79 ) / JVS(          693  )
  W( 79 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 83 ) = W( 83 ) + a*JVS( 695 )
  W( 84 ) = W( 84 ) + a*JVS( 696 )
  W( 85 ) = W( 85 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  a = -W( 80 ) / JVS(          707  )
  W( 80 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 708 )
  W( 83 ) = W( 83 ) + a*JVS( 709 )
  W( 84 ) = W( 84 ) + a*JVS( 710 )
  W( 85 ) = W( 85 ) + a*JVS( 711 )
  W( 86 ) = W( 86 ) + a*JVS( 712 )
  W( 87 ) = W( 87 ) + a*JVS( 713 )
  W( 88 ) = W( 88 ) + a*JVS( 714 )
  W( 89 ) = W( 89 ) + a*JVS( 715 )
  a = -W( 81 ) / JVS(          733  )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 734 )
  W( 83 ) = W( 83 ) + a*JVS( 735 )
  W( 84 ) = W( 84 ) + a*JVS( 736 )
  W( 85 ) = W( 85 ) + a*JVS( 737 )
  W( 86 ) = W( 86 ) + a*JVS( 738 )
  W( 87 ) = W( 87 ) + a*JVS( 739 )
  W( 88 ) = W( 88 ) + a*JVS( 740 )
  W( 89 ) = W( 89 ) + a*JVS( 741 )
  a = -W( 82 ) / JVS(          775  )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 776 )
  W( 84 ) = W( 84 ) + a*JVS( 777 )
  W( 85 ) = W( 85 ) + a*JVS( 778 )
  W( 86 ) = W( 86 ) + a*JVS( 779 )
  W( 87 ) = W( 87 ) + a*JVS( 780 )
  W( 88 ) = W( 88 ) + a*JVS( 781 )
  W( 89 ) = W( 89 ) + a*JVS( 782 )
  W( 90 ) = W( 90 ) + a*JVS( 783 )
  a = -W( 83 ) / JVS(          813  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 814 )
  W( 85 ) = W( 85 ) + a*JVS( 815 )
  W( 86 ) = W( 86 ) + a*JVS( 816 )
  W( 87 ) = W( 87 ) + a*JVS( 817 )
  W( 88 ) = W( 88 ) + a*JVS( 818 )
  W( 89 ) = W( 89 ) + a*JVS( 819 )
  W( 90 ) = W( 90 ) + a*JVS( 820 )
  a = -W( 84 ) / JVS(          885  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 886 )
  W( 86 ) = W( 86 ) + a*JVS( 887 )
  W( 87 ) = W( 87 ) + a*JVS( 888 )
  W( 88 ) = W( 88 ) + a*JVS( 889 )
  W( 89 ) = W( 89 ) + a*JVS( 890 )
  W( 90 ) = W( 90 ) + a*JVS( 891 )
  a = -W( 85 ) / JVS(          950  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 951 )
  W( 87 ) = W( 87 ) + a*JVS( 952 )
  W( 88 ) = W( 88 ) + a*JVS( 953 )
  W( 89 ) = W( 89 ) + a*JVS( 954 )
  W( 90 ) = W( 90 ) + a*JVS( 955 )
  a = -W( 86 ) / JVS(          977  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 978 )
  W( 88 ) = W( 88 ) + a*JVS( 979 )
  W( 89 ) = W( 89 ) + a*JVS( 980 )
  W( 90 ) = W( 90 ) + a*JVS( 981 )
  JVS( 982) = W( 31 )
  JVS( 983) = W( 33 )
  JVS( 984) = W( 36 )
  JVS( 985) = W( 37 )
  JVS( 986) = W( 39 )
  JVS( 987) = W( 50 )
  JVS( 988) = W( 51 )
  JVS( 989) = W( 52 )
  JVS( 990) = W( 55 )
  JVS( 991) = W( 56 )
  JVS( 992) = W( 57 )
  JVS( 993) = W( 58 )
  JVS( 994) = W( 60 )
  JVS( 995) = W( 61 )
  JVS( 996) = W( 63 )
  JVS( 997) = W( 64 )
  JVS( 998) = W( 65 )
  JVS( 999) = W( 67 )
  JVS( 1000) = W( 68 )
  JVS( 1001) = W( 69 )
  JVS( 1002) = W( 72 )
  JVS( 1003) = W( 73 )
  JVS( 1004) = W( 75 )
  JVS( 1005) = W( 76 )
  JVS( 1006) = W( 77 )
  JVS( 1007) = W( 78 )
  JVS( 1008) = W( 79 )
  JVS( 1009) = W( 80 )
  JVS( 1010) = W( 81 )
  JVS( 1011) = W( 82 )
  JVS( 1012) = W( 83 )
  JVS( 1013) = W( 84 )
  JVS( 1014) = W( 85 )
  JVS( 1015) = W( 86 )
  JVS( 1016) = W( 87 )
  JVS( 1017) = W( 88 )
  JVS( 1018) = W( 89 )
  JVS( 1019) = W( 90 )
  IF ( ABS(  JVS( 1071 )) < TINY(a) ) THEN
         IER = 88                                      
         RETURN
  END IF
   W( 19 ) = JVS( 1020 )
   W( 24 ) = JVS( 1021 )
   W( 28 ) = JVS( 1022 )
   W( 29 ) = JVS( 1023 )
   W( 31 ) = JVS( 1024 )
   W( 33 ) = JVS( 1025 )
   W( 34 ) = JVS( 1026 )
   W( 36 ) = JVS( 1027 )
   W( 37 ) = JVS( 1028 )
   W( 39 ) = JVS( 1029 )
   W( 41 ) = JVS( 1030 )
   W( 43 ) = JVS( 1031 )
   W( 44 ) = JVS( 1032 )
   W( 48 ) = JVS( 1033 )
   W( 50 ) = JVS( 1034 )
   W( 51 ) = JVS( 1035 )
   W( 52 ) = JVS( 1036 )
   W( 53 ) = JVS( 1037 )
   W( 54 ) = JVS( 1038 )
   W( 55 ) = JVS( 1039 )
   W( 56 ) = JVS( 1040 )
   W( 57 ) = JVS( 1041 )
   W( 58 ) = JVS( 1042 )
   W( 60 ) = JVS( 1043 )
   W( 61 ) = JVS( 1044 )
   W( 62 ) = JVS( 1045 )
   W( 63 ) = JVS( 1046 )
   W( 64 ) = JVS( 1047 )
   W( 65 ) = JVS( 1048 )
   W( 66 ) = JVS( 1049 )
   W( 67 ) = JVS( 1050 )
   W( 68 ) = JVS( 1051 )
   W( 69 ) = JVS( 1052 )
   W( 70 ) = JVS( 1053 )
   W( 71 ) = JVS( 1054 )
   W( 72 ) = JVS( 1055 )
   W( 73 ) = JVS( 1056 )
   W( 74 ) = JVS( 1057 )
   W( 75 ) = JVS( 1058 )
   W( 76 ) = JVS( 1059 )
   W( 77 ) = JVS( 1060 )
   W( 78 ) = JVS( 1061 )
   W( 79 ) = JVS( 1062 )
   W( 80 ) = JVS( 1063 )
   W( 81 ) = JVS( 1064 )
   W( 82 ) = JVS( 1065 )
   W( 83 ) = JVS( 1066 )
   W( 84 ) = JVS( 1067 )
   W( 85 ) = JVS( 1068 )
   W( 86 ) = JVS( 1069 )
   W( 87 ) = JVS( 1070 )
   W( 88 ) = JVS( 1071 )
   W( 89 ) = JVS( 1072 )
   W( 90 ) = JVS( 1073 )
  a = -W( 19 ) / JVS(           79  )
  W( 19 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 80 )
  W( 89 ) = W( 89 ) + a*JVS( 81 )
  a = -W( 24 ) / JVS(           93  )
  W( 24 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 94 )
  W( 84 ) = W( 84 ) + a*JVS( 95 )
  a = -W( 28 ) / JVS(          112  )
  W( 28 ) = -a
  W( 52 ) = W( 52 ) + a*JVS( 113 )
  W( 84 ) = W( 84 ) + a*JVS( 114 )
  W( 88 ) = W( 88 ) + a*JVS( 115 )
  a = -W( 29 ) / JVS(          116  )
  W( 29 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 117 )
  W( 85 ) = W( 85 ) + a*JVS( 118 )
  W( 88 ) = W( 88 ) + a*JVS( 119 )
  a = -W( 31 ) / JVS(          127  )
  W( 31 ) = -a
  W( 36 ) = W( 36 ) + a*JVS( 128 )
  W( 37 ) = W( 37 ) + a*JVS( 129 )
  W( 39 ) = W( 39 ) + a*JVS( 130 )
  W( 84 ) = W( 84 ) + a*JVS( 131 )
  W( 87 ) = W( 87 ) + a*JVS( 132 )
  W( 88 ) = W( 88 ) + a*JVS( 133 )
  a = -W( 33 ) / JVS(          141  )
  W( 33 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 142 )
  W( 79 ) = W( 79 ) + a*JVS( 143 )
  W( 87 ) = W( 87 ) + a*JVS( 144 )
  W( 88 ) = W( 88 ) + a*JVS( 145 )
  W( 89 ) = W( 89 ) + a*JVS( 146 )
  a = -W( 34 ) / JVS(          147  )
  W( 34 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 148 )
  W( 84 ) = W( 84 ) + a*JVS( 149 )
  W( 85 ) = W( 85 ) + a*JVS( 150 )
  W( 88 ) = W( 88 ) + a*JVS( 151 )
  W( 89 ) = W( 89 ) + a*JVS( 152 )
  a = -W( 36 ) / JVS(          165  )
  W( 36 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 166 )
  W( 84 ) = W( 84 ) + a*JVS( 167 )
  W( 88 ) = W( 88 ) + a*JVS( 168 )
  a = -W( 37 ) / JVS(          170  )
  W( 37 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 171 )
  W( 84 ) = W( 84 ) + a*JVS( 172 )
  W( 88 ) = W( 88 ) + a*JVS( 173 )
  a = -W( 39 ) / JVS(          178  )
  W( 39 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 179 )
  W( 79 ) = W( 79 ) + a*JVS( 180 )
  W( 84 ) = W( 84 ) + a*JVS( 181 )
  W( 88 ) = W( 88 ) + a*JVS( 182 )
  a = -W( 41 ) / JVS(          189  )
  W( 41 ) = -a
  W( 48 ) = W( 48 ) + a*JVS( 190 )
  W( 60 ) = W( 60 ) + a*JVS( 191 )
  W( 84 ) = W( 84 ) + a*JVS( 192 )
  W( 87 ) = W( 87 ) + a*JVS( 193 )
  W( 89 ) = W( 89 ) + a*JVS( 194 )
  a = -W( 43 ) / JVS(          199  )
  W( 43 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 200 )
  W( 53 ) = W( 53 ) + a*JVS( 201 )
  W( 62 ) = W( 62 ) + a*JVS( 202 )
  W( 66 ) = W( 66 ) + a*JVS( 203 )
  W( 71 ) = W( 71 ) + a*JVS( 204 )
  W( 74 ) = W( 74 ) + a*JVS( 205 )
  W( 84 ) = W( 84 ) + a*JVS( 206 )
  W( 85 ) = W( 85 ) + a*JVS( 207 )
  W( 88 ) = W( 88 ) + a*JVS( 208 )
  W( 89 ) = W( 89 ) + a*JVS( 209 )
  a = -W( 44 ) / JVS(          210  )
  W( 44 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 211 )
  W( 79 ) = W( 79 ) + a*JVS( 212 )
  W( 82 ) = W( 82 ) + a*JVS( 213 )
  W( 84 ) = W( 84 ) + a*JVS( 214 )
  W( 88 ) = W( 88 ) + a*JVS( 215 )
  W( 89 ) = W( 89 ) + a*JVS( 216 )
  a = -W( 48 ) / JVS(          250  )
  W( 48 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 251 )
  W( 84 ) = W( 84 ) + a*JVS( 252 )
  W( 89 ) = W( 89 ) + a*JVS( 253 )
  a = -W( 50 ) / JVS(          262  )
  W( 50 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 263 )
  W( 84 ) = W( 84 ) + a*JVS( 264 )
  W( 85 ) = W( 85 ) + a*JVS( 265 )
  W( 88 ) = W( 88 ) + a*JVS( 266 )
  W( 89 ) = W( 89 ) + a*JVS( 267 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  a = -W( 52 ) / JVS(          275  )
  W( 52 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 276 )
  W( 79 ) = W( 79 ) + a*JVS( 277 )
  W( 84 ) = W( 84 ) + a*JVS( 278 )
  W( 85 ) = W( 85 ) + a*JVS( 279 )
  W( 87 ) = W( 87 ) + a*JVS( 280 )
  W( 88 ) = W( 88 ) + a*JVS( 281 )
  W( 89 ) = W( 89 ) + a*JVS( 282 )
  a = -W( 53 ) / JVS(          284  )
  W( 53 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 285 )
  W( 67 ) = W( 67 ) + a*JVS( 286 )
  W( 68 ) = W( 68 ) + a*JVS( 287 )
  W( 71 ) = W( 71 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  W( 80 ) = W( 80 ) + a*JVS( 290 )
  W( 82 ) = W( 82 ) + a*JVS( 291 )
  W( 83 ) = W( 83 ) + a*JVS( 292 )
  W( 84 ) = W( 84 ) + a*JVS( 293 )
  W( 87 ) = W( 87 ) + a*JVS( 294 )
  W( 89 ) = W( 89 ) + a*JVS( 295 )
  a = -W( 54 ) / JVS(          296  )
  W( 54 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 297 )
  W( 79 ) = W( 79 ) + a*JVS( 298 )
  W( 84 ) = W( 84 ) + a*JVS( 299 )
  W( 88 ) = W( 88 ) + a*JVS( 300 )
  W( 89 ) = W( 89 ) + a*JVS( 301 )
  a = -W( 55 ) / JVS(          303  )
  W( 55 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 304 )
  W( 82 ) = W( 82 ) + a*JVS( 305 )
  W( 83 ) = W( 83 ) + a*JVS( 306 )
  W( 84 ) = W( 84 ) + a*JVS( 307 )
  W( 85 ) = W( 85 ) + a*JVS( 308 )
  W( 87 ) = W( 87 ) + a*JVS( 309 )
  W( 89 ) = W( 89 ) + a*JVS( 310 )
  a = -W( 56 ) / JVS(          311  )
  W( 56 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 312 )
  W( 82 ) = W( 82 ) + a*JVS( 313 )
  W( 83 ) = W( 83 ) + a*JVS( 314 )
  W( 84 ) = W( 84 ) + a*JVS( 315 )
  W( 85 ) = W( 85 ) + a*JVS( 316 )
  W( 87 ) = W( 87 ) + a*JVS( 317 )
  W( 89 ) = W( 89 ) + a*JVS( 318 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 325 )
  W( 63 ) = W( 63 ) + a*JVS( 326 )
  W( 79 ) = W( 79 ) + a*JVS( 327 )
  W( 82 ) = W( 82 ) + a*JVS( 328 )
  W( 83 ) = W( 83 ) + a*JVS( 329 )
  W( 84 ) = W( 84 ) + a*JVS( 330 )
  W( 85 ) = W( 85 ) + a*JVS( 331 )
  W( 87 ) = W( 87 ) + a*JVS( 332 )
  W( 89 ) = W( 89 ) + a*JVS( 333 )
  a = -W( 58 ) / JVS(          335  )
  W( 58 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 336 )
  W( 82 ) = W( 82 ) + a*JVS( 337 )
  W( 83 ) = W( 83 ) + a*JVS( 338 )
  W( 84 ) = W( 84 ) + a*JVS( 339 )
  W( 85 ) = W( 85 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 89 ) = W( 89 ) + a*JVS( 342 )
  a = -W( 60 ) / JVS(          372  )
  W( 60 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 373 )
  W( 82 ) = W( 82 ) + a*JVS( 374 )
  W( 83 ) = W( 83 ) + a*JVS( 375 )
  W( 84 ) = W( 84 ) + a*JVS( 376 )
  W( 85 ) = W( 85 ) + a*JVS( 377 )
  W( 87 ) = W( 87 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  a = -W( 61 ) / JVS(          382  )
  W( 61 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 82 ) = W( 82 ) + a*JVS( 384 )
  W( 83 ) = W( 83 ) + a*JVS( 385 )
  W( 84 ) = W( 84 ) + a*JVS( 386 )
  W( 85 ) = W( 85 ) + a*JVS( 387 )
  W( 87 ) = W( 87 ) + a*JVS( 388 )
  W( 88 ) = W( 88 ) + a*JVS( 389 )
  W( 89 ) = W( 89 ) + a*JVS( 390 )
  a = -W( 62 ) / JVS(          400  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 71 ) = W( 71 ) + a*JVS( 404 )
  W( 78 ) = W( 78 ) + a*JVS( 405 )
  W( 79 ) = W( 79 ) + a*JVS( 406 )
  W( 80 ) = W( 80 ) + a*JVS( 407 )
  W( 81 ) = W( 81 ) + a*JVS( 408 )
  W( 82 ) = W( 82 ) + a*JVS( 409 )
  W( 83 ) = W( 83 ) + a*JVS( 410 )
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 85 ) = W( 85 ) + a*JVS( 412 )
  W( 87 ) = W( 87 ) + a*JVS( 413 )
  W( 88 ) = W( 88 ) + a*JVS( 414 )
  W( 89 ) = W( 89 ) + a*JVS( 415 )
  a = -W( 63 ) / JVS(          417  )
  W( 63 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 418 )
  W( 82 ) = W( 82 ) + a*JVS( 419 )
  W( 83 ) = W( 83 ) + a*JVS( 420 )
  W( 84 ) = W( 84 ) + a*JVS( 421 )
  W( 85 ) = W( 85 ) + a*JVS( 422 )
  W( 87 ) = W( 87 ) + a*JVS( 423 )
  W( 89 ) = W( 89 ) + a*JVS( 424 )
  a = -W( 64 ) / JVS(          426  )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 85 ) = W( 85 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 89 ) = W( 89 ) + a*JVS( 432 )
  a = -W( 65 ) / JVS(          434  )
  W( 65 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 435 )
  W( 83 ) = W( 83 ) + a*JVS( 436 )
  W( 84 ) = W( 84 ) + a*JVS( 437 )
  W( 85 ) = W( 85 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 89 ) = W( 89 ) + a*JVS( 440 )
  a = -W( 66 ) / JVS(          466  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 467 )
  W( 68 ) = W( 68 ) + a*JVS( 468 )
  W( 69 ) = W( 69 ) + a*JVS( 469 )
  W( 71 ) = W( 71 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  W( 75 ) = W( 75 ) + a*JVS( 473 )
  W( 76 ) = W( 76 ) + a*JVS( 474 )
  W( 77 ) = W( 77 ) + a*JVS( 475 )
  W( 78 ) = W( 78 ) + a*JVS( 476 )
  W( 79 ) = W( 79 ) + a*JVS( 477 )
  W( 80 ) = W( 80 ) + a*JVS( 478 )
  W( 81 ) = W( 81 ) + a*JVS( 479 )
  W( 82 ) = W( 82 ) + a*JVS( 480 )
  W( 83 ) = W( 83 ) + a*JVS( 481 )
  W( 84 ) = W( 84 ) + a*JVS( 482 )
  W( 85 ) = W( 85 ) + a*JVS( 483 )
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 89 ) = W( 89 ) + a*JVS( 486 )
  W( 90 ) = W( 90 ) + a*JVS( 487 )
  a = -W( 67 ) / JVS(          489  )
  W( 67 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 490 )
  W( 82 ) = W( 82 ) + a*JVS( 491 )
  W( 83 ) = W( 83 ) + a*JVS( 492 )
  W( 84 ) = W( 84 ) + a*JVS( 493 )
  W( 85 ) = W( 85 ) + a*JVS( 494 )
  W( 87 ) = W( 87 ) + a*JVS( 495 )
  W( 88 ) = W( 88 ) + a*JVS( 496 )
  W( 89 ) = W( 89 ) + a*JVS( 497 )
  a = -W( 68 ) / JVS(          499  )
  W( 68 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 500 )
  W( 82 ) = W( 82 ) + a*JVS( 501 )
  W( 83 ) = W( 83 ) + a*JVS( 502 )
  W( 84 ) = W( 84 ) + a*JVS( 503 )
  W( 85 ) = W( 85 ) + a*JVS( 504 )
  W( 87 ) = W( 87 ) + a*JVS( 505 )
  W( 88 ) = W( 88 ) + a*JVS( 506 )
  W( 89 ) = W( 89 ) + a*JVS( 507 )
  a = -W( 69 ) / JVS(          510  )
  W( 69 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 511 )
  W( 79 ) = W( 79 ) + a*JVS( 512 )
  W( 82 ) = W( 82 ) + a*JVS( 513 )
  W( 83 ) = W( 83 ) + a*JVS( 514 )
  W( 84 ) = W( 84 ) + a*JVS( 515 )
  W( 85 ) = W( 85 ) + a*JVS( 516 )
  W( 87 ) = W( 87 ) + a*JVS( 517 )
  W( 89 ) = W( 89 ) + a*JVS( 518 )
  a = -W( 70 ) / JVS(          528  )
  W( 70 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 529 )
  W( 77 ) = W( 77 ) + a*JVS( 530 )
  W( 78 ) = W( 78 ) + a*JVS( 531 )
  W( 79 ) = W( 79 ) + a*JVS( 532 )
  W( 80 ) = W( 80 ) + a*JVS( 533 )
  W( 82 ) = W( 82 ) + a*JVS( 534 )
  W( 83 ) = W( 83 ) + a*JVS( 535 )
  W( 84 ) = W( 84 ) + a*JVS( 536 )
  W( 85 ) = W( 85 ) + a*JVS( 537 )
  W( 87 ) = W( 87 ) + a*JVS( 538 )
  W( 88 ) = W( 88 ) + a*JVS( 539 )
  W( 89 ) = W( 89 ) + a*JVS( 540 )
  a = -W( 71 ) / JVS(          543  )
  W( 71 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 544 )
  W( 82 ) = W( 82 ) + a*JVS( 545 )
  W( 83 ) = W( 83 ) + a*JVS( 546 )
  W( 84 ) = W( 84 ) + a*JVS( 547 )
  W( 85 ) = W( 85 ) + a*JVS( 548 )
  W( 87 ) = W( 87 ) + a*JVS( 549 )
  W( 88 ) = W( 88 ) + a*JVS( 550 )
  W( 89 ) = W( 89 ) + a*JVS( 551 )
  a = -W( 72 ) / JVS(          562  )
  W( 72 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 563 )
  W( 79 ) = W( 79 ) + a*JVS( 564 )
  W( 80 ) = W( 80 ) + a*JVS( 565 )
  W( 81 ) = W( 81 ) + a*JVS( 566 )
  W( 82 ) = W( 82 ) + a*JVS( 567 )
  W( 83 ) = W( 83 ) + a*JVS( 568 )
  W( 84 ) = W( 84 ) + a*JVS( 569 )
  W( 85 ) = W( 85 ) + a*JVS( 570 )
  W( 86 ) = W( 86 ) + a*JVS( 571 )
  W( 87 ) = W( 87 ) + a*JVS( 572 )
  W( 88 ) = W( 88 ) + a*JVS( 573 )
  W( 89 ) = W( 89 ) + a*JVS( 574 )
  a = -W( 73 ) / JVS(          577  )
  W( 73 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 578 )
  W( 82 ) = W( 82 ) + a*JVS( 579 )
  W( 83 ) = W( 83 ) + a*JVS( 580 )
  W( 84 ) = W( 84 ) + a*JVS( 581 )
  W( 85 ) = W( 85 ) + a*JVS( 582 )
  W( 87 ) = W( 87 ) + a*JVS( 583 )
  W( 89 ) = W( 89 ) + a*JVS( 584 )
  a = -W( 74 ) / JVS(          600  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 601 )
  W( 76 ) = W( 76 ) + a*JVS( 602 )
  W( 77 ) = W( 77 ) + a*JVS( 603 )
  W( 78 ) = W( 78 ) + a*JVS( 604 )
  W( 79 ) = W( 79 ) + a*JVS( 605 )
  W( 80 ) = W( 80 ) + a*JVS( 606 )
  W( 81 ) = W( 81 ) + a*JVS( 607 )
  W( 82 ) = W( 82 ) + a*JVS( 608 )
  W( 83 ) = W( 83 ) + a*JVS( 609 )
  W( 84 ) = W( 84 ) + a*JVS( 610 )
  W( 85 ) = W( 85 ) + a*JVS( 611 )
  W( 86 ) = W( 86 ) + a*JVS( 612 )
  W( 87 ) = W( 87 ) + a*JVS( 613 )
  W( 88 ) = W( 88 ) + a*JVS( 614 )
  W( 89 ) = W( 89 ) + a*JVS( 615 )
  W( 90 ) = W( 90 ) + a*JVS( 616 )
  a = -W( 75 ) / JVS(          620  )
  W( 75 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 621 )
  W( 82 ) = W( 82 ) + a*JVS( 622 )
  W( 83 ) = W( 83 ) + a*JVS( 623 )
  W( 84 ) = W( 84 ) + a*JVS( 624 )
  W( 85 ) = W( 85 ) + a*JVS( 625 )
  W( 87 ) = W( 87 ) + a*JVS( 626 )
  W( 89 ) = W( 89 ) + a*JVS( 627 )
  a = -W( 76 ) / JVS(          636  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 637 )
  W( 79 ) = W( 79 ) + a*JVS( 638 )
  W( 82 ) = W( 82 ) + a*JVS( 639 )
  W( 83 ) = W( 83 ) + a*JVS( 640 )
  W( 84 ) = W( 84 ) + a*JVS( 641 )
  W( 85 ) = W( 85 ) + a*JVS( 642 )
  W( 87 ) = W( 87 ) + a*JVS( 643 )
  W( 89 ) = W( 89 ) + a*JVS( 644 )
  a = -W( 77 ) / JVS(          654  )
  W( 77 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 655 )
  W( 82 ) = W( 82 ) + a*JVS( 656 )
  W( 83 ) = W( 83 ) + a*JVS( 657 )
  W( 84 ) = W( 84 ) + a*JVS( 658 )
  W( 85 ) = W( 85 ) + a*JVS( 659 )
  W( 87 ) = W( 87 ) + a*JVS( 660 )
  W( 89 ) = W( 89 ) + a*JVS( 661 )
  a = -W( 78 ) / JVS(          664  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 665 )
  W( 82 ) = W( 82 ) + a*JVS( 666 )
  W( 83 ) = W( 83 ) + a*JVS( 667 )
  W( 84 ) = W( 84 ) + a*JVS( 668 )
  W( 85 ) = W( 85 ) + a*JVS( 669 )
  W( 87 ) = W( 87 ) + a*JVS( 670 )
  W( 88 ) = W( 88 ) + a*JVS( 671 )
  W( 89 ) = W( 89 ) + a*JVS( 672 )
  a = -W( 79 ) / JVS(          693  )
  W( 79 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 83 ) = W( 83 ) + a*JVS( 695 )
  W( 84 ) = W( 84 ) + a*JVS( 696 )
  W( 85 ) = W( 85 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  a = -W( 80 ) / JVS(          707  )
  W( 80 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 708 )
  W( 83 ) = W( 83 ) + a*JVS( 709 )
  W( 84 ) = W( 84 ) + a*JVS( 710 )
  W( 85 ) = W( 85 ) + a*JVS( 711 )
  W( 86 ) = W( 86 ) + a*JVS( 712 )
  W( 87 ) = W( 87 ) + a*JVS( 713 )
  W( 88 ) = W( 88 ) + a*JVS( 714 )
  W( 89 ) = W( 89 ) + a*JVS( 715 )
  a = -W( 81 ) / JVS(          733  )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 734 )
  W( 83 ) = W( 83 ) + a*JVS( 735 )
  W( 84 ) = W( 84 ) + a*JVS( 736 )
  W( 85 ) = W( 85 ) + a*JVS( 737 )
  W( 86 ) = W( 86 ) + a*JVS( 738 )
  W( 87 ) = W( 87 ) + a*JVS( 739 )
  W( 88 ) = W( 88 ) + a*JVS( 740 )
  W( 89 ) = W( 89 ) + a*JVS( 741 )
  a = -W( 82 ) / JVS(          775  )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 776 )
  W( 84 ) = W( 84 ) + a*JVS( 777 )
  W( 85 ) = W( 85 ) + a*JVS( 778 )
  W( 86 ) = W( 86 ) + a*JVS( 779 )
  W( 87 ) = W( 87 ) + a*JVS( 780 )
  W( 88 ) = W( 88 ) + a*JVS( 781 )
  W( 89 ) = W( 89 ) + a*JVS( 782 )
  W( 90 ) = W( 90 ) + a*JVS( 783 )
  a = -W( 83 ) / JVS(          813  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 814 )
  W( 85 ) = W( 85 ) + a*JVS( 815 )
  W( 86 ) = W( 86 ) + a*JVS( 816 )
  W( 87 ) = W( 87 ) + a*JVS( 817 )
  W( 88 ) = W( 88 ) + a*JVS( 818 )
  W( 89 ) = W( 89 ) + a*JVS( 819 )
  W( 90 ) = W( 90 ) + a*JVS( 820 )
  a = -W( 84 ) / JVS(          885  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 886 )
  W( 86 ) = W( 86 ) + a*JVS( 887 )
  W( 87 ) = W( 87 ) + a*JVS( 888 )
  W( 88 ) = W( 88 ) + a*JVS( 889 )
  W( 89 ) = W( 89 ) + a*JVS( 890 )
  W( 90 ) = W( 90 ) + a*JVS( 891 )
  a = -W( 85 ) / JVS(          950  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 951 )
  W( 87 ) = W( 87 ) + a*JVS( 952 )
  W( 88 ) = W( 88 ) + a*JVS( 953 )
  W( 89 ) = W( 89 ) + a*JVS( 954 )
  W( 90 ) = W( 90 ) + a*JVS( 955 )
  a = -W( 86 ) / JVS(          977  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 978 )
  W( 88 ) = W( 88 ) + a*JVS( 979 )
  W( 89 ) = W( 89 ) + a*JVS( 980 )
  W( 90 ) = W( 90 ) + a*JVS( 981 )
  a = -W( 87 ) / JVS(         1016  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 1017 )
  W( 89 ) = W( 89 ) + a*JVS( 1018 )
  W( 90 ) = W( 90 ) + a*JVS( 1019 )
  JVS( 1020) = W( 19 )
  JVS( 1021) = W( 24 )
  JVS( 1022) = W( 28 )
  JVS( 1023) = W( 29 )
  JVS( 1024) = W( 31 )
  JVS( 1025) = W( 33 )
  JVS( 1026) = W( 34 )
  JVS( 1027) = W( 36 )
  JVS( 1028) = W( 37 )
  JVS( 1029) = W( 39 )
  JVS( 1030) = W( 41 )
  JVS( 1031) = W( 43 )
  JVS( 1032) = W( 44 )
  JVS( 1033) = W( 48 )
  JVS( 1034) = W( 50 )
  JVS( 1035) = W( 51 )
  JVS( 1036) = W( 52 )
  JVS( 1037) = W( 53 )
  JVS( 1038) = W( 54 )
  JVS( 1039) = W( 55 )
  JVS( 1040) = W( 56 )
  JVS( 1041) = W( 57 )
  JVS( 1042) = W( 58 )
  JVS( 1043) = W( 60 )
  JVS( 1044) = W( 61 )
  JVS( 1045) = W( 62 )
  JVS( 1046) = W( 63 )
  JVS( 1047) = W( 64 )
  JVS( 1048) = W( 65 )
  JVS( 1049) = W( 66 )
  JVS( 1050) = W( 67 )
  JVS( 1051) = W( 68 )
  JVS( 1052) = W( 69 )
  JVS( 1053) = W( 70 )
  JVS( 1054) = W( 71 )
  JVS( 1055) = W( 72 )
  JVS( 1056) = W( 73 )
  JVS( 1057) = W( 74 )
  JVS( 1058) = W( 75 )
  JVS( 1059) = W( 76 )
  JVS( 1060) = W( 77 )
  JVS( 1061) = W( 78 )
  JVS( 1062) = W( 79 )
  JVS( 1063) = W( 80 )
  JVS( 1064) = W( 81 )
  JVS( 1065) = W( 82 )
  JVS( 1066) = W( 83 )
  JVS( 1067) = W( 84 )
  JVS( 1068) = W( 85 )
  JVS( 1069) = W( 86 )
  JVS( 1070) = W( 87 )
  JVS( 1071) = W( 88 )
  JVS( 1072) = W( 89 )
  JVS( 1073) = W( 90 )
  IF ( ABS(  JVS( 1119 )) < TINY(a) ) THEN
         IER = 89                                      
         RETURN
  END IF
   W( 19 ) = JVS( 1074 )
   W( 29 ) = JVS( 1075 )
   W( 33 ) = JVS( 1076 )
   W( 38 ) = JVS( 1077 )
   W( 42 ) = JVS( 1078 )
   W( 43 ) = JVS( 1079 )
   W( 44 ) = JVS( 1080 )
   W( 45 ) = JVS( 1081 )
   W( 47 ) = JVS( 1082 )
   W( 48 ) = JVS( 1083 )
   W( 49 ) = JVS( 1084 )
   W( 50 ) = JVS( 1085 )
   W( 51 ) = JVS( 1086 )
   W( 53 ) = JVS( 1087 )
   W( 54 ) = JVS( 1088 )
   W( 55 ) = JVS( 1089 )
   W( 56 ) = JVS( 1090 )
   W( 58 ) = JVS( 1091 )
   W( 61 ) = JVS( 1092 )
   W( 62 ) = JVS( 1093 )
   W( 63 ) = JVS( 1094 )
   W( 64 ) = JVS( 1095 )
   W( 65 ) = JVS( 1096 )
   W( 66 ) = JVS( 1097 )
   W( 67 ) = JVS( 1098 )
   W( 68 ) = JVS( 1099 )
   W( 69 ) = JVS( 1100 )
   W( 71 ) = JVS( 1101 )
   W( 72 ) = JVS( 1102 )
   W( 73 ) = JVS( 1103 )
   W( 74 ) = JVS( 1104 )
   W( 75 ) = JVS( 1105 )
   W( 76 ) = JVS( 1106 )
   W( 77 ) = JVS( 1107 )
   W( 78 ) = JVS( 1108 )
   W( 79 ) = JVS( 1109 )
   W( 80 ) = JVS( 1110 )
   W( 81 ) = JVS( 1111 )
   W( 82 ) = JVS( 1112 )
   W( 83 ) = JVS( 1113 )
   W( 84 ) = JVS( 1114 )
   W( 85 ) = JVS( 1115 )
   W( 86 ) = JVS( 1116 )
   W( 87 ) = JVS( 1117 )
   W( 88 ) = JVS( 1118 )
   W( 89 ) = JVS( 1119 )
   W( 90 ) = JVS( 1120 )
  a = -W( 19 ) / JVS(           79  )
  W( 19 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 80 )
  W( 89 ) = W( 89 ) + a*JVS( 81 )
  a = -W( 29 ) / JVS(          116  )
  W( 29 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 117 )
  W( 85 ) = W( 85 ) + a*JVS( 118 )
  W( 88 ) = W( 88 ) + a*JVS( 119 )
  a = -W( 33 ) / JVS(          141  )
  W( 33 ) = -a
  W( 51 ) = W( 51 ) + a*JVS( 142 )
  W( 79 ) = W( 79 ) + a*JVS( 143 )
  W( 87 ) = W( 87 ) + a*JVS( 144 )
  W( 88 ) = W( 88 ) + a*JVS( 145 )
  W( 89 ) = W( 89 ) + a*JVS( 146 )
  a = -W( 38 ) / JVS(          174  )
  W( 38 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 175 )
  W( 84 ) = W( 84 ) + a*JVS( 176 )
  W( 89 ) = W( 89 ) + a*JVS( 177 )
  a = -W( 42 ) / JVS(          195  )
  W( 42 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 196 )
  W( 84 ) = W( 84 ) + a*JVS( 197 )
  W( 89 ) = W( 89 ) + a*JVS( 198 )
  a = -W( 43 ) / JVS(          199  )
  W( 43 ) = -a
  W( 50 ) = W( 50 ) + a*JVS( 200 )
  W( 53 ) = W( 53 ) + a*JVS( 201 )
  W( 62 ) = W( 62 ) + a*JVS( 202 )
  W( 66 ) = W( 66 ) + a*JVS( 203 )
  W( 71 ) = W( 71 ) + a*JVS( 204 )
  W( 74 ) = W( 74 ) + a*JVS( 205 )
  W( 84 ) = W( 84 ) + a*JVS( 206 )
  W( 85 ) = W( 85 ) + a*JVS( 207 )
  W( 88 ) = W( 88 ) + a*JVS( 208 )
  W( 89 ) = W( 89 ) + a*JVS( 209 )
  a = -W( 44 ) / JVS(          210  )
  W( 44 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 211 )
  W( 79 ) = W( 79 ) + a*JVS( 212 )
  W( 82 ) = W( 82 ) + a*JVS( 213 )
  W( 84 ) = W( 84 ) + a*JVS( 214 )
  W( 88 ) = W( 88 ) + a*JVS( 215 )
  W( 89 ) = W( 89 ) + a*JVS( 216 )
  a = -W( 45 ) / JVS(          217  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 218 )
  W( 84 ) = W( 84 ) + a*JVS( 219 )
  W( 89 ) = W( 89 ) + a*JVS( 220 )
  a = -W( 47 ) / JVS(          246  )
  W( 47 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 247 )
  W( 84 ) = W( 84 ) + a*JVS( 248 )
  W( 89 ) = W( 89 ) + a*JVS( 249 )
  a = -W( 48 ) / JVS(          250  )
  W( 48 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 251 )
  W( 84 ) = W( 84 ) + a*JVS( 252 )
  W( 89 ) = W( 89 ) + a*JVS( 253 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 255 )
  W( 84 ) = W( 84 ) + a*JVS( 256 )
  W( 89 ) = W( 89 ) + a*JVS( 257 )
  a = -W( 50 ) / JVS(          262  )
  W( 50 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 263 )
  W( 84 ) = W( 84 ) + a*JVS( 264 )
  W( 85 ) = W( 85 ) + a*JVS( 265 )
  W( 88 ) = W( 88 ) + a*JVS( 266 )
  W( 89 ) = W( 89 ) + a*JVS( 267 )
  a = -W( 51 ) / JVS(          268  )
  W( 51 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 269 )
  W( 84 ) = W( 84 ) + a*JVS( 270 )
  W( 89 ) = W( 89 ) + a*JVS( 271 )
  a = -W( 53 ) / JVS(          284  )
  W( 53 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 285 )
  W( 67 ) = W( 67 ) + a*JVS( 286 )
  W( 68 ) = W( 68 ) + a*JVS( 287 )
  W( 71 ) = W( 71 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  W( 80 ) = W( 80 ) + a*JVS( 290 )
  W( 82 ) = W( 82 ) + a*JVS( 291 )
  W( 83 ) = W( 83 ) + a*JVS( 292 )
  W( 84 ) = W( 84 ) + a*JVS( 293 )
  W( 87 ) = W( 87 ) + a*JVS( 294 )
  W( 89 ) = W( 89 ) + a*JVS( 295 )
  a = -W( 54 ) / JVS(          296  )
  W( 54 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 297 )
  W( 79 ) = W( 79 ) + a*JVS( 298 )
  W( 84 ) = W( 84 ) + a*JVS( 299 )
  W( 88 ) = W( 88 ) + a*JVS( 300 )
  W( 89 ) = W( 89 ) + a*JVS( 301 )
  a = -W( 55 ) / JVS(          303  )
  W( 55 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 304 )
  W( 82 ) = W( 82 ) + a*JVS( 305 )
  W( 83 ) = W( 83 ) + a*JVS( 306 )
  W( 84 ) = W( 84 ) + a*JVS( 307 )
  W( 85 ) = W( 85 ) + a*JVS( 308 )
  W( 87 ) = W( 87 ) + a*JVS( 309 )
  W( 89 ) = W( 89 ) + a*JVS( 310 )
  a = -W( 56 ) / JVS(          311  )
  W( 56 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 312 )
  W( 82 ) = W( 82 ) + a*JVS( 313 )
  W( 83 ) = W( 83 ) + a*JVS( 314 )
  W( 84 ) = W( 84 ) + a*JVS( 315 )
  W( 85 ) = W( 85 ) + a*JVS( 316 )
  W( 87 ) = W( 87 ) + a*JVS( 317 )
  W( 89 ) = W( 89 ) + a*JVS( 318 )
  a = -W( 58 ) / JVS(          335  )
  W( 58 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 336 )
  W( 82 ) = W( 82 ) + a*JVS( 337 )
  W( 83 ) = W( 83 ) + a*JVS( 338 )
  W( 84 ) = W( 84 ) + a*JVS( 339 )
  W( 85 ) = W( 85 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 89 ) = W( 89 ) + a*JVS( 342 )
  a = -W( 61 ) / JVS(          382  )
  W( 61 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 383 )
  W( 82 ) = W( 82 ) + a*JVS( 384 )
  W( 83 ) = W( 83 ) + a*JVS( 385 )
  W( 84 ) = W( 84 ) + a*JVS( 386 )
  W( 85 ) = W( 85 ) + a*JVS( 387 )
  W( 87 ) = W( 87 ) + a*JVS( 388 )
  W( 88 ) = W( 88 ) + a*JVS( 389 )
  W( 89 ) = W( 89 ) + a*JVS( 390 )
  a = -W( 62 ) / JVS(          400  )
  W( 62 ) = -a
  W( 63 ) = W( 63 ) + a*JVS( 401 )
  W( 67 ) = W( 67 ) + a*JVS( 402 )
  W( 68 ) = W( 68 ) + a*JVS( 403 )
  W( 71 ) = W( 71 ) + a*JVS( 404 )
  W( 78 ) = W( 78 ) + a*JVS( 405 )
  W( 79 ) = W( 79 ) + a*JVS( 406 )
  W( 80 ) = W( 80 ) + a*JVS( 407 )
  W( 81 ) = W( 81 ) + a*JVS( 408 )
  W( 82 ) = W( 82 ) + a*JVS( 409 )
  W( 83 ) = W( 83 ) + a*JVS( 410 )
  W( 84 ) = W( 84 ) + a*JVS( 411 )
  W( 85 ) = W( 85 ) + a*JVS( 412 )
  W( 87 ) = W( 87 ) + a*JVS( 413 )
  W( 88 ) = W( 88 ) + a*JVS( 414 )
  W( 89 ) = W( 89 ) + a*JVS( 415 )
  a = -W( 63 ) / JVS(          417  )
  W( 63 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 418 )
  W( 82 ) = W( 82 ) + a*JVS( 419 )
  W( 83 ) = W( 83 ) + a*JVS( 420 )
  W( 84 ) = W( 84 ) + a*JVS( 421 )
  W( 85 ) = W( 85 ) + a*JVS( 422 )
  W( 87 ) = W( 87 ) + a*JVS( 423 )
  W( 89 ) = W( 89 ) + a*JVS( 424 )
  a = -W( 64 ) / JVS(          426  )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 85 ) = W( 85 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 89 ) = W( 89 ) + a*JVS( 432 )
  a = -W( 65 ) / JVS(          434  )
  W( 65 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 435 )
  W( 83 ) = W( 83 ) + a*JVS( 436 )
  W( 84 ) = W( 84 ) + a*JVS( 437 )
  W( 85 ) = W( 85 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 89 ) = W( 89 ) + a*JVS( 440 )
  a = -W( 66 ) / JVS(          466  )
  W( 66 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 467 )
  W( 68 ) = W( 68 ) + a*JVS( 468 )
  W( 69 ) = W( 69 ) + a*JVS( 469 )
  W( 71 ) = W( 71 ) + a*JVS( 470 )
  W( 72 ) = W( 72 ) + a*JVS( 471 )
  W( 73 ) = W( 73 ) + a*JVS( 472 )
  W( 75 ) = W( 75 ) + a*JVS( 473 )
  W( 76 ) = W( 76 ) + a*JVS( 474 )
  W( 77 ) = W( 77 ) + a*JVS( 475 )
  W( 78 ) = W( 78 ) + a*JVS( 476 )
  W( 79 ) = W( 79 ) + a*JVS( 477 )
  W( 80 ) = W( 80 ) + a*JVS( 478 )
  W( 81 ) = W( 81 ) + a*JVS( 479 )
  W( 82 ) = W( 82 ) + a*JVS( 480 )
  W( 83 ) = W( 83 ) + a*JVS( 481 )
  W( 84 ) = W( 84 ) + a*JVS( 482 )
  W( 85 ) = W( 85 ) + a*JVS( 483 )
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 89 ) = W( 89 ) + a*JVS( 486 )
  W( 90 ) = W( 90 ) + a*JVS( 487 )
  a = -W( 67 ) / JVS(          489  )
  W( 67 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 490 )
  W( 82 ) = W( 82 ) + a*JVS( 491 )
  W( 83 ) = W( 83 ) + a*JVS( 492 )
  W( 84 ) = W( 84 ) + a*JVS( 493 )
  W( 85 ) = W( 85 ) + a*JVS( 494 )
  W( 87 ) = W( 87 ) + a*JVS( 495 )
  W( 88 ) = W( 88 ) + a*JVS( 496 )
  W( 89 ) = W( 89 ) + a*JVS( 497 )
  a = -W( 68 ) / JVS(          499  )
  W( 68 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 500 )
  W( 82 ) = W( 82 ) + a*JVS( 501 )
  W( 83 ) = W( 83 ) + a*JVS( 502 )
  W( 84 ) = W( 84 ) + a*JVS( 503 )
  W( 85 ) = W( 85 ) + a*JVS( 504 )
  W( 87 ) = W( 87 ) + a*JVS( 505 )
  W( 88 ) = W( 88 ) + a*JVS( 506 )
  W( 89 ) = W( 89 ) + a*JVS( 507 )
  a = -W( 69 ) / JVS(          510  )
  W( 69 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 511 )
  W( 79 ) = W( 79 ) + a*JVS( 512 )
  W( 82 ) = W( 82 ) + a*JVS( 513 )
  W( 83 ) = W( 83 ) + a*JVS( 514 )
  W( 84 ) = W( 84 ) + a*JVS( 515 )
  W( 85 ) = W( 85 ) + a*JVS( 516 )
  W( 87 ) = W( 87 ) + a*JVS( 517 )
  W( 89 ) = W( 89 ) + a*JVS( 518 )
  a = -W( 71 ) / JVS(          543  )
  W( 71 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 544 )
  W( 82 ) = W( 82 ) + a*JVS( 545 )
  W( 83 ) = W( 83 ) + a*JVS( 546 )
  W( 84 ) = W( 84 ) + a*JVS( 547 )
  W( 85 ) = W( 85 ) + a*JVS( 548 )
  W( 87 ) = W( 87 ) + a*JVS( 549 )
  W( 88 ) = W( 88 ) + a*JVS( 550 )
  W( 89 ) = W( 89 ) + a*JVS( 551 )
  a = -W( 72 ) / JVS(          562  )
  W( 72 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 563 )
  W( 79 ) = W( 79 ) + a*JVS( 564 )
  W( 80 ) = W( 80 ) + a*JVS( 565 )
  W( 81 ) = W( 81 ) + a*JVS( 566 )
  W( 82 ) = W( 82 ) + a*JVS( 567 )
  W( 83 ) = W( 83 ) + a*JVS( 568 )
  W( 84 ) = W( 84 ) + a*JVS( 569 )
  W( 85 ) = W( 85 ) + a*JVS( 570 )
  W( 86 ) = W( 86 ) + a*JVS( 571 )
  W( 87 ) = W( 87 ) + a*JVS( 572 )
  W( 88 ) = W( 88 ) + a*JVS( 573 )
  W( 89 ) = W( 89 ) + a*JVS( 574 )
  a = -W( 73 ) / JVS(          577  )
  W( 73 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 578 )
  W( 82 ) = W( 82 ) + a*JVS( 579 )
  W( 83 ) = W( 83 ) + a*JVS( 580 )
  W( 84 ) = W( 84 ) + a*JVS( 581 )
  W( 85 ) = W( 85 ) + a*JVS( 582 )
  W( 87 ) = W( 87 ) + a*JVS( 583 )
  W( 89 ) = W( 89 ) + a*JVS( 584 )
  a = -W( 74 ) / JVS(          600  )
  W( 74 ) = -a
  W( 75 ) = W( 75 ) + a*JVS( 601 )
  W( 76 ) = W( 76 ) + a*JVS( 602 )
  W( 77 ) = W( 77 ) + a*JVS( 603 )
  W( 78 ) = W( 78 ) + a*JVS( 604 )
  W( 79 ) = W( 79 ) + a*JVS( 605 )
  W( 80 ) = W( 80 ) + a*JVS( 606 )
  W( 81 ) = W( 81 ) + a*JVS( 607 )
  W( 82 ) = W( 82 ) + a*JVS( 608 )
  W( 83 ) = W( 83 ) + a*JVS( 609 )
  W( 84 ) = W( 84 ) + a*JVS( 610 )
  W( 85 ) = W( 85 ) + a*JVS( 611 )
  W( 86 ) = W( 86 ) + a*JVS( 612 )
  W( 87 ) = W( 87 ) + a*JVS( 613 )
  W( 88 ) = W( 88 ) + a*JVS( 614 )
  W( 89 ) = W( 89 ) + a*JVS( 615 )
  W( 90 ) = W( 90 ) + a*JVS( 616 )
  a = -W( 75 ) / JVS(          620  )
  W( 75 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 621 )
  W( 82 ) = W( 82 ) + a*JVS( 622 )
  W( 83 ) = W( 83 ) + a*JVS( 623 )
  W( 84 ) = W( 84 ) + a*JVS( 624 )
  W( 85 ) = W( 85 ) + a*JVS( 625 )
  W( 87 ) = W( 87 ) + a*JVS( 626 )
  W( 89 ) = W( 89 ) + a*JVS( 627 )
  a = -W( 76 ) / JVS(          636  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 637 )
  W( 79 ) = W( 79 ) + a*JVS( 638 )
  W( 82 ) = W( 82 ) + a*JVS( 639 )
  W( 83 ) = W( 83 ) + a*JVS( 640 )
  W( 84 ) = W( 84 ) + a*JVS( 641 )
  W( 85 ) = W( 85 ) + a*JVS( 642 )
  W( 87 ) = W( 87 ) + a*JVS( 643 )
  W( 89 ) = W( 89 ) + a*JVS( 644 )
  a = -W( 77 ) / JVS(          654  )
  W( 77 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 655 )
  W( 82 ) = W( 82 ) + a*JVS( 656 )
  W( 83 ) = W( 83 ) + a*JVS( 657 )
  W( 84 ) = W( 84 ) + a*JVS( 658 )
  W( 85 ) = W( 85 ) + a*JVS( 659 )
  W( 87 ) = W( 87 ) + a*JVS( 660 )
  W( 89 ) = W( 89 ) + a*JVS( 661 )
  a = -W( 78 ) / JVS(          664  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 665 )
  W( 82 ) = W( 82 ) + a*JVS( 666 )
  W( 83 ) = W( 83 ) + a*JVS( 667 )
  W( 84 ) = W( 84 ) + a*JVS( 668 )
  W( 85 ) = W( 85 ) + a*JVS( 669 )
  W( 87 ) = W( 87 ) + a*JVS( 670 )
  W( 88 ) = W( 88 ) + a*JVS( 671 )
  W( 89 ) = W( 89 ) + a*JVS( 672 )
  a = -W( 79 ) / JVS(          693  )
  W( 79 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 83 ) = W( 83 ) + a*JVS( 695 )
  W( 84 ) = W( 84 ) + a*JVS( 696 )
  W( 85 ) = W( 85 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  a = -W( 80 ) / JVS(          707  )
  W( 80 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 708 )
  W( 83 ) = W( 83 ) + a*JVS( 709 )
  W( 84 ) = W( 84 ) + a*JVS( 710 )
  W( 85 ) = W( 85 ) + a*JVS( 711 )
  W( 86 ) = W( 86 ) + a*JVS( 712 )
  W( 87 ) = W( 87 ) + a*JVS( 713 )
  W( 88 ) = W( 88 ) + a*JVS( 714 )
  W( 89 ) = W( 89 ) + a*JVS( 715 )
  a = -W( 81 ) / JVS(          733  )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 734 )
  W( 83 ) = W( 83 ) + a*JVS( 735 )
  W( 84 ) = W( 84 ) + a*JVS( 736 )
  W( 85 ) = W( 85 ) + a*JVS( 737 )
  W( 86 ) = W( 86 ) + a*JVS( 738 )
  W( 87 ) = W( 87 ) + a*JVS( 739 )
  W( 88 ) = W( 88 ) + a*JVS( 740 )
  W( 89 ) = W( 89 ) + a*JVS( 741 )
  a = -W( 82 ) / JVS(          775  )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 776 )
  W( 84 ) = W( 84 ) + a*JVS( 777 )
  W( 85 ) = W( 85 ) + a*JVS( 778 )
  W( 86 ) = W( 86 ) + a*JVS( 779 )
  W( 87 ) = W( 87 ) + a*JVS( 780 )
  W( 88 ) = W( 88 ) + a*JVS( 781 )
  W( 89 ) = W( 89 ) + a*JVS( 782 )
  W( 90 ) = W( 90 ) + a*JVS( 783 )
  a = -W( 83 ) / JVS(          813  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 814 )
  W( 85 ) = W( 85 ) + a*JVS( 815 )
  W( 86 ) = W( 86 ) + a*JVS( 816 )
  W( 87 ) = W( 87 ) + a*JVS( 817 )
  W( 88 ) = W( 88 ) + a*JVS( 818 )
  W( 89 ) = W( 89 ) + a*JVS( 819 )
  W( 90 ) = W( 90 ) + a*JVS( 820 )
  a = -W( 84 ) / JVS(          885  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 886 )
  W( 86 ) = W( 86 ) + a*JVS( 887 )
  W( 87 ) = W( 87 ) + a*JVS( 888 )
  W( 88 ) = W( 88 ) + a*JVS( 889 )
  W( 89 ) = W( 89 ) + a*JVS( 890 )
  W( 90 ) = W( 90 ) + a*JVS( 891 )
  a = -W( 85 ) / JVS(          950  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 951 )
  W( 87 ) = W( 87 ) + a*JVS( 952 )
  W( 88 ) = W( 88 ) + a*JVS( 953 )
  W( 89 ) = W( 89 ) + a*JVS( 954 )
  W( 90 ) = W( 90 ) + a*JVS( 955 )
  a = -W( 86 ) / JVS(          977  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 978 )
  W( 88 ) = W( 88 ) + a*JVS( 979 )
  W( 89 ) = W( 89 ) + a*JVS( 980 )
  W( 90 ) = W( 90 ) + a*JVS( 981 )
  a = -W( 87 ) / JVS(         1016  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 1017 )
  W( 89 ) = W( 89 ) + a*JVS( 1018 )
  W( 90 ) = W( 90 ) + a*JVS( 1019 )
  a = -W( 88 ) / JVS(         1071  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 1072 )
  W( 90 ) = W( 90 ) + a*JVS( 1073 )
  JVS( 1074) = W( 19 )
  JVS( 1075) = W( 29 )
  JVS( 1076) = W( 33 )
  JVS( 1077) = W( 38 )
  JVS( 1078) = W( 42 )
  JVS( 1079) = W( 43 )
  JVS( 1080) = W( 44 )
  JVS( 1081) = W( 45 )
  JVS( 1082) = W( 47 )
  JVS( 1083) = W( 48 )
  JVS( 1084) = W( 49 )
  JVS( 1085) = W( 50 )
  JVS( 1086) = W( 51 )
  JVS( 1087) = W( 53 )
  JVS( 1088) = W( 54 )
  JVS( 1089) = W( 55 )
  JVS( 1090) = W( 56 )
  JVS( 1091) = W( 58 )
  JVS( 1092) = W( 61 )
  JVS( 1093) = W( 62 )
  JVS( 1094) = W( 63 )
  JVS( 1095) = W( 64 )
  JVS( 1096) = W( 65 )
  JVS( 1097) = W( 66 )
  JVS( 1098) = W( 67 )
  JVS( 1099) = W( 68 )
  JVS( 1100) = W( 69 )
  JVS( 1101) = W( 71 )
  JVS( 1102) = W( 72 )
  JVS( 1103) = W( 73 )
  JVS( 1104) = W( 74 )
  JVS( 1105) = W( 75 )
  JVS( 1106) = W( 76 )
  JVS( 1107) = W( 77 )
  JVS( 1108) = W( 78 )
  JVS( 1109) = W( 79 )
  JVS( 1110) = W( 80 )
  JVS( 1111) = W( 81 )
  JVS( 1112) = W( 82 )
  JVS( 1113) = W( 83 )
  JVS( 1114) = W( 84 )
  JVS( 1115) = W( 85 )
  JVS( 1116) = W( 86 )
  JVS( 1117) = W( 87 )
  JVS( 1118) = W( 88 )
  JVS( 1119) = W( 89 )
  JVS( 1120) = W( 90 )
  IF ( ABS(  JVS( 1147 )) < TINY(a) ) THEN
         IER = 90                                      
         RETURN
  END IF
   W( 25 ) = JVS( 1121 )
   W( 45 ) = JVS( 1122 )
   W( 47 ) = JVS( 1123 )
   W( 49 ) = JVS( 1124 )
   W( 59 ) = JVS( 1125 )
   W( 64 ) = JVS( 1126 )
   W( 65 ) = JVS( 1127 )
   W( 69 ) = JVS( 1128 )
   W( 70 ) = JVS( 1129 )
   W( 71 ) = JVS( 1130 )
   W( 73 ) = JVS( 1131 )
   W( 75 ) = JVS( 1132 )
   W( 76 ) = JVS( 1133 )
   W( 77 ) = JVS( 1134 )
   W( 78 ) = JVS( 1135 )
   W( 79 ) = JVS( 1136 )
   W( 80 ) = JVS( 1137 )
   W( 81 ) = JVS( 1138 )
   W( 82 ) = JVS( 1139 )
   W( 83 ) = JVS( 1140 )
   W( 84 ) = JVS( 1141 )
   W( 85 ) = JVS( 1142 )
   W( 86 ) = JVS( 1143 )
   W( 87 ) = JVS( 1144 )
   W( 88 ) = JVS( 1145 )
   W( 89 ) = JVS( 1146 )
   W( 90 ) = JVS( 1147 )
  a = -W( 25 ) / JVS(           96  )
  W( 25 ) = -a
  W( 49 ) = W( 49 ) + a*JVS( 97 )
  W( 73 ) = W( 73 ) + a*JVS( 98 )
  W( 75 ) = W( 75 ) + a*JVS( 99 )
  W( 79 ) = W( 79 ) + a*JVS( 100 )
  W( 84 ) = W( 84 ) + a*JVS( 101 )
  a = -W( 45 ) / JVS(          217  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 218 )
  W( 84 ) = W( 84 ) + a*JVS( 219 )
  W( 89 ) = W( 89 ) + a*JVS( 220 )
  a = -W( 47 ) / JVS(          246  )
  W( 47 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 247 )
  W( 84 ) = W( 84 ) + a*JVS( 248 )
  W( 89 ) = W( 89 ) + a*JVS( 249 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 255 )
  W( 84 ) = W( 84 ) + a*JVS( 256 )
  W( 89 ) = W( 89 ) + a*JVS( 257 )
  a = -W( 59 ) / JVS(          350  )
  W( 59 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 351 )
  W( 65 ) = W( 65 ) + a*JVS( 352 )
  W( 69 ) = W( 69 ) + a*JVS( 353 )
  W( 70 ) = W( 70 ) + a*JVS( 354 )
  W( 71 ) = W( 71 ) + a*JVS( 355 )
  W( 73 ) = W( 73 ) + a*JVS( 356 )
  W( 75 ) = W( 75 ) + a*JVS( 357 )
  W( 76 ) = W( 76 ) + a*JVS( 358 )
  W( 77 ) = W( 77 ) + a*JVS( 359 )
  W( 79 ) = W( 79 ) + a*JVS( 360 )
  W( 80 ) = W( 80 ) + a*JVS( 361 )
  W( 81 ) = W( 81 ) + a*JVS( 362 )
  W( 82 ) = W( 82 ) + a*JVS( 363 )
  W( 83 ) = W( 83 ) + a*JVS( 364 )
  W( 84 ) = W( 84 ) + a*JVS( 365 )
  W( 85 ) = W( 85 ) + a*JVS( 366 )
  W( 86 ) = W( 86 ) + a*JVS( 367 )
  W( 87 ) = W( 87 ) + a*JVS( 368 )
  W( 89 ) = W( 89 ) + a*JVS( 369 )
  a = -W( 64 ) / JVS(          426  )
  W( 64 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 427 )
  W( 83 ) = W( 83 ) + a*JVS( 428 )
  W( 84 ) = W( 84 ) + a*JVS( 429 )
  W( 85 ) = W( 85 ) + a*JVS( 430 )
  W( 87 ) = W( 87 ) + a*JVS( 431 )
  W( 89 ) = W( 89 ) + a*JVS( 432 )
  a = -W( 65 ) / JVS(          434  )
  W( 65 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 435 )
  W( 83 ) = W( 83 ) + a*JVS( 436 )
  W( 84 ) = W( 84 ) + a*JVS( 437 )
  W( 85 ) = W( 85 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 89 ) = W( 89 ) + a*JVS( 440 )
  a = -W( 69 ) / JVS(          510  )
  W( 69 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 511 )
  W( 79 ) = W( 79 ) + a*JVS( 512 )
  W( 82 ) = W( 82 ) + a*JVS( 513 )
  W( 83 ) = W( 83 ) + a*JVS( 514 )
  W( 84 ) = W( 84 ) + a*JVS( 515 )
  W( 85 ) = W( 85 ) + a*JVS( 516 )
  W( 87 ) = W( 87 ) + a*JVS( 517 )
  W( 89 ) = W( 89 ) + a*JVS( 518 )
  a = -W( 70 ) / JVS(          528  )
  W( 70 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 529 )
  W( 77 ) = W( 77 ) + a*JVS( 530 )
  W( 78 ) = W( 78 ) + a*JVS( 531 )
  W( 79 ) = W( 79 ) + a*JVS( 532 )
  W( 80 ) = W( 80 ) + a*JVS( 533 )
  W( 82 ) = W( 82 ) + a*JVS( 534 )
  W( 83 ) = W( 83 ) + a*JVS( 535 )
  W( 84 ) = W( 84 ) + a*JVS( 536 )
  W( 85 ) = W( 85 ) + a*JVS( 537 )
  W( 87 ) = W( 87 ) + a*JVS( 538 )
  W( 88 ) = W( 88 ) + a*JVS( 539 )
  W( 89 ) = W( 89 ) + a*JVS( 540 )
  a = -W( 71 ) / JVS(          543  )
  W( 71 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 544 )
  W( 82 ) = W( 82 ) + a*JVS( 545 )
  W( 83 ) = W( 83 ) + a*JVS( 546 )
  W( 84 ) = W( 84 ) + a*JVS( 547 )
  W( 85 ) = W( 85 ) + a*JVS( 548 )
  W( 87 ) = W( 87 ) + a*JVS( 549 )
  W( 88 ) = W( 88 ) + a*JVS( 550 )
  W( 89 ) = W( 89 ) + a*JVS( 551 )
  a = -W( 73 ) / JVS(          577  )
  W( 73 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 578 )
  W( 82 ) = W( 82 ) + a*JVS( 579 )
  W( 83 ) = W( 83 ) + a*JVS( 580 )
  W( 84 ) = W( 84 ) + a*JVS( 581 )
  W( 85 ) = W( 85 ) + a*JVS( 582 )
  W( 87 ) = W( 87 ) + a*JVS( 583 )
  W( 89 ) = W( 89 ) + a*JVS( 584 )
  a = -W( 75 ) / JVS(          620  )
  W( 75 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 621 )
  W( 82 ) = W( 82 ) + a*JVS( 622 )
  W( 83 ) = W( 83 ) + a*JVS( 623 )
  W( 84 ) = W( 84 ) + a*JVS( 624 )
  W( 85 ) = W( 85 ) + a*JVS( 625 )
  W( 87 ) = W( 87 ) + a*JVS( 626 )
  W( 89 ) = W( 89 ) + a*JVS( 627 )
  a = -W( 76 ) / JVS(          636  )
  W( 76 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 637 )
  W( 79 ) = W( 79 ) + a*JVS( 638 )
  W( 82 ) = W( 82 ) + a*JVS( 639 )
  W( 83 ) = W( 83 ) + a*JVS( 640 )
  W( 84 ) = W( 84 ) + a*JVS( 641 )
  W( 85 ) = W( 85 ) + a*JVS( 642 )
  W( 87 ) = W( 87 ) + a*JVS( 643 )
  W( 89 ) = W( 89 ) + a*JVS( 644 )
  a = -W( 77 ) / JVS(          654  )
  W( 77 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 655 )
  W( 82 ) = W( 82 ) + a*JVS( 656 )
  W( 83 ) = W( 83 ) + a*JVS( 657 )
  W( 84 ) = W( 84 ) + a*JVS( 658 )
  W( 85 ) = W( 85 ) + a*JVS( 659 )
  W( 87 ) = W( 87 ) + a*JVS( 660 )
  W( 89 ) = W( 89 ) + a*JVS( 661 )
  a = -W( 78 ) / JVS(          664  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 665 )
  W( 82 ) = W( 82 ) + a*JVS( 666 )
  W( 83 ) = W( 83 ) + a*JVS( 667 )
  W( 84 ) = W( 84 ) + a*JVS( 668 )
  W( 85 ) = W( 85 ) + a*JVS( 669 )
  W( 87 ) = W( 87 ) + a*JVS( 670 )
  W( 88 ) = W( 88 ) + a*JVS( 671 )
  W( 89 ) = W( 89 ) + a*JVS( 672 )
  a = -W( 79 ) / JVS(          693  )
  W( 79 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 83 ) = W( 83 ) + a*JVS( 695 )
  W( 84 ) = W( 84 ) + a*JVS( 696 )
  W( 85 ) = W( 85 ) + a*JVS( 697 )
  W( 87 ) = W( 87 ) + a*JVS( 698 )
  W( 88 ) = W( 88 ) + a*JVS( 699 )
  W( 89 ) = W( 89 ) + a*JVS( 700 )
  a = -W( 80 ) / JVS(          707  )
  W( 80 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 708 )
  W( 83 ) = W( 83 ) + a*JVS( 709 )
  W( 84 ) = W( 84 ) + a*JVS( 710 )
  W( 85 ) = W( 85 ) + a*JVS( 711 )
  W( 86 ) = W( 86 ) + a*JVS( 712 )
  W( 87 ) = W( 87 ) + a*JVS( 713 )
  W( 88 ) = W( 88 ) + a*JVS( 714 )
  W( 89 ) = W( 89 ) + a*JVS( 715 )
  a = -W( 81 ) / JVS(          733  )
  W( 81 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 734 )
  W( 83 ) = W( 83 ) + a*JVS( 735 )
  W( 84 ) = W( 84 ) + a*JVS( 736 )
  W( 85 ) = W( 85 ) + a*JVS( 737 )
  W( 86 ) = W( 86 ) + a*JVS( 738 )
  W( 87 ) = W( 87 ) + a*JVS( 739 )
  W( 88 ) = W( 88 ) + a*JVS( 740 )
  W( 89 ) = W( 89 ) + a*JVS( 741 )
  a = -W( 82 ) / JVS(          775  )
  W( 82 ) = -a
  W( 83 ) = W( 83 ) + a*JVS( 776 )
  W( 84 ) = W( 84 ) + a*JVS( 777 )
  W( 85 ) = W( 85 ) + a*JVS( 778 )
  W( 86 ) = W( 86 ) + a*JVS( 779 )
  W( 87 ) = W( 87 ) + a*JVS( 780 )
  W( 88 ) = W( 88 ) + a*JVS( 781 )
  W( 89 ) = W( 89 ) + a*JVS( 782 )
  W( 90 ) = W( 90 ) + a*JVS( 783 )
  a = -W( 83 ) / JVS(          813  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 814 )
  W( 85 ) = W( 85 ) + a*JVS( 815 )
  W( 86 ) = W( 86 ) + a*JVS( 816 )
  W( 87 ) = W( 87 ) + a*JVS( 817 )
  W( 88 ) = W( 88 ) + a*JVS( 818 )
  W( 89 ) = W( 89 ) + a*JVS( 819 )
  W( 90 ) = W( 90 ) + a*JVS( 820 )
  a = -W( 84 ) / JVS(          885  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 886 )
  W( 86 ) = W( 86 ) + a*JVS( 887 )
  W( 87 ) = W( 87 ) + a*JVS( 888 )
  W( 88 ) = W( 88 ) + a*JVS( 889 )
  W( 89 ) = W( 89 ) + a*JVS( 890 )
  W( 90 ) = W( 90 ) + a*JVS( 891 )
  a = -W( 85 ) / JVS(          950  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 951 )
  W( 87 ) = W( 87 ) + a*JVS( 952 )
  W( 88 ) = W( 88 ) + a*JVS( 953 )
  W( 89 ) = W( 89 ) + a*JVS( 954 )
  W( 90 ) = W( 90 ) + a*JVS( 955 )
  a = -W( 86 ) / JVS(          977  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 978 )
  W( 88 ) = W( 88 ) + a*JVS( 979 )
  W( 89 ) = W( 89 ) + a*JVS( 980 )
  W( 90 ) = W( 90 ) + a*JVS( 981 )
  a = -W( 87 ) / JVS(         1016  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 1017 )
  W( 89 ) = W( 89 ) + a*JVS( 1018 )
  W( 90 ) = W( 90 ) + a*JVS( 1019 )
  a = -W( 88 ) / JVS(         1071  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 1072 )
  W( 90 ) = W( 90 ) + a*JVS( 1073 )
  a = -W( 89 ) / JVS(         1119  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 1120 )
  JVS( 1121) = W( 25 )
  JVS( 1122) = W( 45 )
  JVS( 1123) = W( 47 )
  JVS( 1124) = W( 49 )
  JVS( 1125) = W( 59 )
  JVS( 1126) = W( 64 )
  JVS( 1127) = W( 65 )
  JVS( 1128) = W( 69 )
  JVS( 1129) = W( 70 )
  JVS( 1130) = W( 71 )
  JVS( 1131) = W( 73 )
  JVS( 1132) = W( 75 )
  JVS( 1133) = W( 76 )
  JVS( 1134) = W( 77 )
  JVS( 1135) = W( 78 )
  JVS( 1136) = W( 79 )
  JVS( 1137) = W( 80 )
  JVS( 1138) = W( 81 )
  JVS( 1139) = W( 82 )
  JVS( 1140) = W( 83 )
  JVS( 1141) = W( 84 )
  JVS( 1142) = W( 85 )
  JVS( 1143) = W( 86 )
  JVS( 1144) = W( 87 )
  JVS( 1145) = W( 88 )
  JVS( 1146) = W( 89 )
  JVS( 1147) = W( 90 )
   
   END SUBROUTINE decomp_racm_soa_vbs_aqchem
 


END MODULE racm_soa_vbs_aqchem_Integrator
