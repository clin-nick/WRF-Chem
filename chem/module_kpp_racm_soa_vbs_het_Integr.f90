
































MODULE racm_soa_vbs_het_Integrator

 USE racm_soa_vbs_het_Parameters
 USE racm_soa_vbs_het_Precision
 USE racm_soa_vbs_het_JacobianSP

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

SUBROUTINE  racm_soa_vbs_het_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE racm_soa_vbs_het_Parameters

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

   CALL racm_soa_vbs_het_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  racm_soa_vbs_het_INTEGRATE


SUBROUTINE  racm_soa_vbs_het_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE racm_soa_vbs_het_Parameters

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
      CALL racm_soa_vbs_het_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL racm_soa_vbs_het_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = racm_soa_vbs_het_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL racm_soa_vbs_het_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL racm_soa_vbs_het_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL racm_soa_vbs_het_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL racm_soa_vbs_het_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL racm_soa_vbs_het_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL racm_soa_vbs_het_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL racm_soa_vbs_het_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL racm_soa_vbs_het_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL racm_soa_vbs_het_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL racm_soa_vbs_het_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL racm_soa_vbs_het_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL racm_soa_vbs_het_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL racm_soa_vbs_het_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL racm_soa_vbs_het_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL racm_soa_vbs_het_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  racm_soa_vbs_het_ros_ErrorMsg(Code,T,H,IERR)



   USE racm_soa_vbs_het_Precision

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

 END SUBROUTINE  racm_soa_vbs_het_ros_ErrorMsg


 SUBROUTINE  racm_soa_vbs_het_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL racm_soa_vbs_het_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL racm_soa_vbs_het_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL racm_soa_vbs_het_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL racm_soa_vbs_het_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL racm_soa_vbs_het_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL racm_soa_vbs_het_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL racm_soa_vbs_het_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL racm_soa_vbs_het_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL racm_soa_vbs_het_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL racm_soa_vbs_het_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL racm_soa_vbs_het_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL racm_soa_vbs_het_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL racm_soa_vbs_het_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL racm_soa_vbs_het_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL racm_soa_vbs_het_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL racm_soa_vbs_het_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL racm_soa_vbs_het_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL racm_soa_vbs_het_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL racm_soa_vbs_het_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL racm_soa_vbs_het_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = racm_soa_vbs_het_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL racm_soa_vbs_het_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  racm_soa_vbs_het_ros_Integrator



  REAL(kind=dp) FUNCTION  racm_soa_vbs_het_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    racm_soa_vbs_het_ros_ErrorNorm = Err

  END FUNCTION  racm_soa_vbs_het_ros_ErrorNorm



  SUBROUTINE racm_soa_vbs_het_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL racm_soa_vbs_het_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL racm_soa_vbs_het_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL racm_soa_vbs_het_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  racm_soa_vbs_het_ros_FunTimeDeriv



  SUBROUTINE  racm_soa_vbs_het_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL racm_soa_vbs_het_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL racm_soa_vbs_het_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL racm_soa_vbs_het_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  racm_soa_vbs_het_ros_PrepareMatrix



  SUBROUTINE  racm_soa_vbs_het_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 



CALL decomp_racm_soa_vbs_het ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  racm_soa_vbs_het_ros_Decomp



  SUBROUTINE  racm_soa_vbs_het_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL racm_soa_vbs_het_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  racm_soa_vbs_het_ros_Solve




  SUBROUTINE  racm_soa_vbs_het_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  racm_soa_vbs_het_Ros2



  SUBROUTINE  racm_soa_vbs_het_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racm_soa_vbs_het_Ros3





  SUBROUTINE  racm_soa_vbs_het_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racm_soa_vbs_het_Ros4


  SUBROUTINE  racm_soa_vbs_het_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racm_soa_vbs_het_Rodas3


  SUBROUTINE  racm_soa_vbs_het_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  racm_soa_vbs_het_Rodas4




END SUBROUTINE  racm_soa_vbs_het_Rosenbrock




SUBROUTINE  racm_soa_vbs_het_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE racm_soa_vbs_het_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL racm_soa_vbs_het_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  racm_soa_vbs_het_FunTemplate



SUBROUTINE  racm_soa_vbs_het_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE racm_soa_vbs_het_Parameters
 
 USE racm_soa_vbs_het_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL racm_soa_vbs_het_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  racm_soa_vbs_het_JacTemplate

















SUBROUTINE racm_soa_vbs_het_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(95)
  A(2) = RCT(2)*V(85)
  A(3) = RCT(3)*V(85)
  A(4) = RCT(4)*V(29)
  A(5) = RCT(5)*V(44)
  A(6) = RCT(6)*V(27)
  A(7) = RCT(7)*V(96)
  A(8) = RCT(8)*V(96)
  A(9) = RCT(9)*V(33)
  A(10) = RCT(10)*V(73)
  A(11) = RCT(11)*V(73)
  A(12) = RCT(12)*V(78)
  A(13) = RCT(13)*V(25)
  A(14) = RCT(14)*V(92)
  A(15) = RCT(15)*V(45)
  A(16) = RCT(16)*V(66)
  A(17) = RCT(17)*V(60)
  A(18) = RCT(18)*V(60)
  A(19) = RCT(19)*V(69)
  A(20) = RCT(20)*V(79)
  A(21) = RCT(21)*V(75)
  A(22) = RCT(22)*V(64)
  A(23) = RCT(23)*V(37)
  A(24) = RCT(24)*V(13)
  A(25) = RCT(25)*V(19)
  A(26) = RCT(26)*V(18)
  A(27) = RCT(27)*V(36)
  A(28) = RCT(28)*V(30)*F(2)
  A(29) = RCT(29)*V(30)*V(85)
  A(30) = RCT(30)*V(14)*F(2)
  A(31) = RCT(31)*V(14)*F(1)
  A(32) = RCT(32)*V(85)*V(91)
  A(33) = RCT(33)*V(85)*V(89)
  A(34) = RCT(34)*V(89)*V(91)
  A(35) = RCT(35)*V(33)*V(91)
  A(36) = RCT(36)*V(89)*V(89)
  A(37) = RCT(37)*V(89)*V(89)*F(1)
  A(38) = RCT(38)*V(30)*V(94)
  A(39) = RCT(39)*V(30)*V(95)
  A(40) = RCT(40)*V(30)*V(95)
  A(41) = RCT(41)*V(91)*V(94)
  A(42) = RCT(42)*V(91)*V(95)
  A(43) = RCT(43)*V(91)*V(96)
  A(44) = RCT(44)*V(89)*V(94)
  A(45) = RCT(45)*V(89)*V(95)
  A(46) = RCT(46)*V(27)
  A(47) = RCT(47)*V(89)*V(96)
  A(48) = RCT(48)*V(29)*V(91)
  A(49) = RCT(49)*V(44)*V(91)
  A(50) = RCT(50)*V(27)*V(91)
  A(51) = RCT(51)*V(85)*V(94)
  A(52) = RCT(52)*V(85)*V(95)
  A(53) = RCT(53)*V(94)*V(94)*F(2)
  A(54) = RCT(54)*V(94)*V(96)
  A(55) = RCT(55)*V(95)*V(96)
  A(56) = RCT(56)*V(95)*V(96)
  A(57) = RCT(57)*V(17)
  A(58) = RCT(58)*V(96)*V(96)
  A(59) = RCT(59)*V(91)*F(2)
  A(60) = RCT(60)*V(15)*V(91)
  A(61) = RCT(61)*V(52)*V(91)
  A(62) = RCT(62)*V(22)*V(91)
  A(63) = RCT(63)*V(28)*V(91)
  A(64) = RCT(64)*V(35)*V(91)
  A(65) = RCT(65)*V(50)*V(91)
  A(66) = RCT(66)*V(49)*V(91)
  A(67) = RCT(67)*V(41)*V(91)
  A(68) = RCT(68)*V(46)*V(91)
  A(69) = RCT(69)*V(53)*V(91)
  A(70) = RCT(70)*V(82)*V(91)
  A(71) = RCT(71)*V(81)*V(91)
  A(72) = RCT(72)*V(55)*V(91)
  A(73) = RCT(73)*V(59)*V(91)
  A(74) = RCT(74)*V(47)*V(91)
  A(75) = RCT(75)*V(51)*V(91)
  A(76) = RCT(76)*V(23)*V(91)
  A(77) = RCT(77)*V(24)*V(91)
  A(78) = RCT(78)*V(57)*V(91)
  A(79) = RCT(79)*V(73)*V(91)
  A(80) = RCT(80)*V(78)*V(91)
  A(81) = RCT(81)*V(66)*V(91)
  A(82) = RCT(82)*V(37)*V(91)
  A(83) = RCT(83)*V(60)*V(91)
  A(84) = RCT(84)*V(69)*V(91)
  A(85) = RCT(85)*V(64)*V(91)
  A(86) = RCT(86)*V(79)*V(91)
  A(87) = RCT(87)*V(21)*V(91)
  A(88) = RCT(88)*V(25)*V(91)
  A(89) = RCT(89)*V(91)*V(92)
  A(90) = RCT(90)*V(45)*V(91)
  A(91) = RCT(91)*V(48)*V(91)
  A(92) = RCT(92)*V(58)*V(91)
  A(93) = RCT(93)*V(75)*V(91)
  A(94) = RCT(94)*V(73)*V(96)
  A(95) = RCT(95)*V(78)*V(96)
  A(96) = RCT(96)*V(60)*V(96)
  A(97) = RCT(97)*V(69)*V(96)
  A(98) = RCT(98)*V(20)*V(91)
  A(99) = RCT(99)*V(79)*V(96)
  A(100) = RCT(100)*V(57)*V(96)
  A(101) = RCT(101)*V(53)*V(96)
  A(102) = RCT(102)*V(82)*V(96)
  A(103) = RCT(103)*V(81)*V(96)
  A(104) = RCT(104)*V(55)*V(96)
  A(105) = RCT(105)*V(59)*V(96)
  A(106) = RCT(106)*V(47)*V(96)
  A(107) = RCT(107)*V(51)*V(96)
  A(108) = RCT(108)*V(58)*V(96)
  A(109) = RCT(109)*V(53)*V(85)
  A(110) = RCT(110)*V(82)*V(85)
  A(111) = RCT(111)*V(81)*V(85)
  A(112) = RCT(112)*V(55)*V(85)
  A(113) = RCT(113)*V(59)*V(85)
  A(114) = RCT(114)*V(47)*V(85)
  A(115) = RCT(115)*V(51)*V(85)
  A(116) = RCT(116)*V(64)*V(85)
  A(117) = RCT(117)*V(79)*V(85)
  A(118) = RCT(118)*V(58)*V(85)
  A(119) = RCT(119)*V(31)*V(95)
  A(120) = RCT(120)*V(31)*V(89)
  A(121) = RCT(121)*V(38)*V(95)
  A(122) = RCT(122)*V(38)*F(2)
  A(123) = RCT(123)*V(38)*V(85)
  A(124) = RCT(124)*V(39)*V(95)
  A(125) = RCT(125)*V(39)*F(2)
  A(126) = RCT(126)*V(39)*V(85)
  A(127) = RCT(127)*V(40)*V(95)
  A(128) = RCT(128)*V(40)*F(2)
  A(129) = RCT(129)*V(40)*V(85)
  A(130) = RCT(130)*V(93)*V(95)
  A(131) = RCT(131)*V(48)
  A(132) = RCT(132)*V(74)*V(95)
  A(133) = RCT(133)*V(58)
  A(134) = RCT(134)*V(90)*V(94)
  A(135) = RCT(135)*V(94)*V(97)
  A(136) = RCT(136)*V(86)*V(94)
  A(137) = RCT(137)*V(76)*V(94)
  A(138) = RCT(138)*V(77)*V(94)
  A(139) = RCT(139)*V(62)*V(94)
  A(140) = RCT(140)*V(61)*V(94)
  A(141) = RCT(141)*V(70)*V(94)
  A(142) = RCT(142)*V(67)*V(94)
  A(143) = RCT(143)*V(65)*V(94)
  A(144) = RCT(144)*V(63)*V(94)
  A(145) = RCT(145)*V(71)*V(94)
  A(146) = RCT(146)*V(72)*V(94)
  A(147) = RCT(147)*V(68)*V(94)
  A(148) = RCT(148)*V(93)*V(94)
  A(149) = RCT(149)*V(74)*V(94)
  A(150) = RCT(150)*V(87)*V(94)
  A(151) = RCT(151)*V(83)*V(94)
  A(152) = RCT(152)*V(84)*V(94)
  A(153) = RCT(153)*V(89)*V(90)
  A(154) = RCT(154)*V(89)*V(97)
  A(155) = RCT(155)*V(86)*V(89)
  A(156) = RCT(156)*V(76)*V(89)
  A(157) = RCT(157)*V(77)*V(89)
  A(158) = RCT(158)*V(62)*V(89)
  A(159) = RCT(159)*V(61)*V(89)
  A(160) = RCT(160)*V(70)*V(89)
  A(161) = RCT(161)*V(67)*V(89)
  A(162) = RCT(162)*V(65)*V(89)
  A(163) = RCT(163)*V(63)*V(89)
  A(164) = RCT(164)*V(71)*V(89)
  A(165) = RCT(165)*V(72)*V(89)
  A(166) = RCT(166)*V(68)*V(89)
  A(167) = RCT(167)*V(89)*V(93)
  A(168) = RCT(168)*V(89)*V(93)
  A(169) = RCT(169)*V(74)*V(89)
  A(170) = RCT(170)*V(74)*V(89)
  A(171) = RCT(171)*V(87)*V(89)
  A(172) = RCT(172)*V(83)*V(89)
  A(173) = RCT(173)*V(84)*V(89)
  A(174) = RCT(174)*V(90)*V(90)
  A(175) = RCT(175)*V(90)*V(97)
  A(176) = RCT(176)*V(86)*V(90)
  A(177) = RCT(177)*V(76)*V(90)
  A(178) = RCT(178)*V(77)*V(90)
  A(179) = RCT(179)*V(62)*V(90)
  A(180) = RCT(180)*V(61)*V(90)
  A(181) = RCT(181)*V(70)*V(90)
  A(182) = RCT(182)*V(67)*V(90)
  A(183) = RCT(183)*V(65)*V(90)
  A(184) = RCT(184)*V(63)*V(90)
  A(185) = RCT(185)*V(71)*V(90)
  A(186) = RCT(186)*V(72)*V(90)
  A(187) = RCT(187)*V(68)*V(90)
  A(188) = RCT(188)*V(90)*V(93)
  A(189) = RCT(189)*V(90)*V(93)
  A(190) = RCT(190)*V(74)*V(90)
  A(191) = RCT(191)*V(74)*V(90)
  A(192) = RCT(192)*V(87)*V(90)
  A(193) = RCT(193)*V(83)*V(90)
  A(194) = RCT(194)*V(84)*V(90)
  A(195) = RCT(195)*V(93)*V(97)
  A(196) = RCT(196)*V(86)*V(93)
  A(197) = RCT(197)*V(76)*V(93)
  A(198) = RCT(198)*V(77)*V(93)
  A(199) = RCT(199)*V(62)*V(93)
  A(200) = RCT(200)*V(61)*V(93)
  A(201) = RCT(201)*V(70)*V(93)
  A(202) = RCT(202)*V(67)*V(93)
  A(203) = RCT(203)*V(65)*V(93)
  A(204) = RCT(204)*V(63)*V(93)
  A(205) = RCT(205)*V(71)*V(93)
  A(206) = RCT(206)*V(72)*V(93)
  A(207) = RCT(207)*V(68)*V(93)
  A(208) = RCT(208)*V(93)*V(93)
  A(209) = RCT(209)*V(74)*V(93)
  A(210) = RCT(210)*V(87)*V(93)
  A(211) = RCT(211)*V(83)*V(93)
  A(212) = RCT(212)*V(84)*V(93)
  A(213) = RCT(213)*V(83)*V(83)
  A(214) = RCT(214)*V(83)*V(84)
  A(215) = RCT(215)*V(84)*V(84)
  A(216) = RCT(216)*V(90)*V(96)
  A(217) = RCT(217)*V(96)*V(97)
  A(218) = RCT(218)*V(86)*V(96)
  A(219) = RCT(219)*V(76)*V(96)
  A(220) = RCT(220)*V(77)*V(96)
  A(221) = RCT(221)*V(62)*V(96)
  A(222) = RCT(222)*V(61)*V(96)
  A(223) = RCT(223)*V(70)*V(96)
  A(224) = RCT(224)*V(26)*V(91)
  A(225) = RCT(225)*V(65)*V(96)
  A(226) = RCT(226)*V(63)*V(96)
  A(227) = RCT(227)*V(71)*V(96)
  A(228) = RCT(228)*V(72)*V(96)
  A(229) = RCT(229)*V(68)*V(96)
  A(230) = RCT(230)*V(93)*V(96)
  A(231) = RCT(231)*V(74)*V(96)
  A(232) = RCT(232)*V(87)*V(96)
  A(233) = RCT(233)*V(83)*V(96)
  A(234) = RCT(234)*V(84)*V(96)
  A(235) = RCT(235)*V(80)*V(89)
  A(236) = RCT(236)*V(80)*V(90)
  A(237) = RCT(237)*V(80)*V(93)
  A(238) = RCT(238)*V(80)*V(80)
  A(239) = RCT(239)*V(80)*V(94)
  A(240) = RCT(240)*V(80)*V(96)
  A(241) = RCT(241)*V(67)*V(67)
  A(242) = RCT(242)*V(16)*V(91)
  A(243) = RCT(243)*V(42)*V(91)
  A(244) = RCT(244)*V(56)*V(94)
  A(245) = RCT(245)*V(56)*V(89)
  A(246) = RCT(246)*V(56)*V(56)
  A(247) = RCT(247)*V(56)*V(95)
  A(248) = RCT(248)*V(26)
  A(249) = RCT(249)*V(43)*V(91)
  A(250) = RCT(250)*V(43)*V(85)
  A(251) = RCT(251)*V(43)*V(96)
  A(252) = RCT(252)*V(54)*V(91)
  A(253) = RCT(253)*V(54)*V(96)
  A(254) = RCT(254)*V(54)*V(85)
  A(255) = RCT(255)*V(8)*V(91)
  A(256) = RCT(256)*V(7)*V(91)
  A(257) = RCT(257)*V(6)*V(91)
  A(258) = RCT(258)*V(12)*V(91)
  A(259) = RCT(259)*V(11)*V(91)
  A(260) = RCT(260)*V(10)*V(91)
  A(261) = RCT(261)*V(34)*V(91)
  A(262) = RCT(262)*V(85)*V(88)
  A(263) = RCT(263)*V(32)*V(32)
  A(264) = RCT(264)*V(32)*V(94)
  A(265) = RCT(265)*V(32)*V(89)
  A(266) = RCT(266)*V(88)*V(95)
  A(267) = RCT(267)*V(35)*V(88)
  A(268) = RCT(268)*V(50)*V(88)
  A(269) = RCT(269)*V(49)*V(88)
  A(270) = RCT(270)*V(41)*V(88)
  A(271) = RCT(271)*V(46)*V(88)
  A(272) = RCT(272)*V(53)*V(88)
  A(273) = RCT(273)*V(82)*V(88)
  A(274) = RCT(274)*V(81)*V(88)
  A(275) = RCT(275)*V(59)*V(88)
  A(276) = RCT(276)*V(36)*V(91)
  A(277) = RCT(277)*V(73)*V(88)
  A(278) = RCT(278)*V(78)*V(88)
  A(279) = RCT(279)*V(23)*V(88)
  A(280) = RCT(280)*V(24)*V(88)


  Vdot(1) = A(60)
  Vdot(2) = A(61)
  Vdot(3) = 0.036*A(66)+0.37*A(109)+0.14*A(110)+0.15*A(112)+0.28*A(113)+0.01*A(115)+0.45*A(116)+0.11*A(117)+0.11*A(118)&
              &+0.05*A(249)+0.039*A(250)
  Vdot(4) = 0.1*A(110)+0.14*A(111)+0.07*A(115)+0.21*A(117)+A(168)+A(170)+A(189)+A(191)+0.5*A(195)+0.499*A(196)+0.495&
              &*A(197)+0.495*A(198)+0.5*A(199)+0.499*A(200)+0.49*A(201)+0.494*A(202)+0.5*A(210)+0.5*A(211)+0.484*A(212)&
              &+0.053*A(250)+0.14*A(254)
  Vdot(5) = 1.075*A(257)
  Vdot(6) = 1.075*A(256)-A(257)
  Vdot(7) = 1.075*A(255)-A(256)
  Vdot(8) = -A(255)
  Vdot(9) = 1.075*A(260)
  Vdot(10) = 1.075*A(259)-A(260)
  Vdot(11) = 1.075*A(258)-A(259)
  Vdot(12) = -A(258)
  Vdot(13) = -A(24)+0.3*A(263)
  Vdot(14) = A(2)-A(30)-A(31)
  Vdot(15) = -A(60)
  Vdot(16) = A(161)-A(242)
  Vdot(17) = A(56)-A(57)
  Vdot(18) = -A(26)+A(266)
  Vdot(19) = -A(25)+A(265)
  Vdot(20) = -A(98)+A(245)
  Vdot(21) = 0.35*A(86)-A(87)
  Vdot(22) = -A(62)+A(243)
  Vdot(23) = -A(76)-A(279)
  Vdot(24) = -A(77)-A(280)
  Vdot(25) = -A(13)-A(88)+A(153)
  Vdot(26) = -A(224)+A(247)-A(248)
  Vdot(27) = -A(6)+A(45)-A(46)-A(50)
  Vdot(28) = -A(63)+A(224)+A(243)+0.25*A(244)+A(246)
  Vdot(29) = -A(4)+A(41)-A(48)+A(121)+A(124)+A(127)
  Vdot(30) = A(1)+A(3)+A(8)-A(28)-A(29)+A(30)-A(38)-A(39)-A(40)+0.09*A(112)
  Vdot(31) = 0.1*A(78)+A(100)-A(119)-A(120)
  Vdot(32) = A(262)-2*A(263)-A(264)-A(265)
  Vdot(33) = -A(9)-A(35)+A(36)+A(37)+0.006*A(110)+0.011*A(111)+0.001*A(112)+0.09*A(113)+0.02*A(114)+0.02*A(115)+0.011&
               &*A(254)
  Vdot(34) = -A(261)+A(267)+A(268)+A(269)+A(270)+A(271)+0.3*A(274)+0.15*A(275)+A(277)+A(278)+A(279)+A(280)
  Vdot(35) = -A(64)+0.06*A(110)+0.07*A(111)+0.07*A(254)-A(267)
  Vdot(36) = -A(27)+A(272)+A(273)+0.7*A(274)+0.85*A(275)-A(276)
  Vdot(37) = -A(23)+0.024*A(68)-A(82)+0.6*A(92)+0.3*A(192)
  Vdot(38) = 0.9*A(76)-A(121)-A(122)-A(123)
  Vdot(39) = 0.9*A(77)-A(124)-A(125)-A(126)
  Vdot(40) = 0.85*A(78)-A(127)-A(128)-A(129)
  Vdot(41) = -A(67)-1.11*A(270)-0.33*A(273)+0.1*A(274)
  Vdot(42) = A(105)+0.046*A(142)-A(243)
  Vdot(43) = -A(249)-A(250)-A(251)
  Vdot(44) = -A(5)+A(42)+0.3*A(47)-A(49)+A(94)+A(95)+A(96)+A(97)+0.5*A(99)+A(100)
  Vdot(45) = -A(15)-A(90)+0.11*A(117)+A(167)
  Vdot(46) = -A(68)-1.11*A(271)-0.33*A(273)+0.1*A(274)
  Vdot(47) = -A(74)-A(106)-A(114)
  Vdot(48) = -A(91)+0.4*A(92)+0.4*A(108)+0.3*A(118)+A(130)-A(131)
  Vdot(49) = -A(66)-1.11*A(269)-0.33*A(273)+0.1*A(274)
  Vdot(50) = -A(65)+0.03*A(110)+0.06*A(111)+0.06*A(254)-A(268)
  Vdot(51) = -A(75)-A(107)-A(115)
  Vdot(52) = A(10)+A(11)+A(12)+1.87*A(17)+1.55*A(18)+A(19)+A(22)+A(27)-A(61)+A(62)+0.036*A(66)+A(79)+2*A(83)+A(84)+A(94)&
               &+2*A(96)+A(97)+0.43*A(109)+0.37*A(110)+0.3*A(111)+0.36*A(112)+0.14*A(113)+0.14*A(114)+0.14*A(115)+0.22&
               &*A(116)+0.66*A(117)+0.13*A(118)+0.25*A(244)+0.5*A(246)+0.3*A(254)+A(276)+A(277)
  Vdot(53) = -A(69)-A(101)-A(109)-A(272)
  Vdot(54) = -A(252)-A(253)-A(254)
  Vdot(55) = -A(72)-A(104)-A(112)
  Vdot(56) = A(85)+A(98)+0.1*A(113)-A(244)-A(245)-2*A(246)-A(247)+A(248)
  Vdot(57) = -A(78)-A(100)+0.1*A(119)+A(120)+A(121)+0.02*A(122)+A(123)+A(124)+0.02*A(125)+A(126)+A(127)+0.02*A(128)&
               &+A(129)
  Vdot(58) = -A(92)-A(108)-A(118)+A(132)-A(133)
  Vdot(59) = -A(73)-A(105)-A(113)-A(275)
  Vdot(60) = -A(17)-A(18)+0.036*A(66)-A(83)+0.15*A(86)-A(96)+0.25*A(99)+0.5*A(117)+0.063*A(136)+1.2*A(145)+0.35*A(146)&
               &+A(147)+0.119*A(176)+0.65*A(185)+0.37*A(186)+A(187)+0.1*A(196)+0.65*A(205)+0.37*A(206)+A(207)+0.063*A(218)&
               &+1.3*A(227)+0.74*A(228)+A(229)
  Vdot(61) = A(70)-A(140)-A(159)-A(180)-A(200)-A(222)
  Vdot(62) = A(69)-A(139)-A(158)-A(179)-A(199)-A(221)
  Vdot(63) = A(75)-A(144)-A(163)-A(184)-A(204)-A(226)
  Vdot(64) = -A(22)-A(85)+0.9*A(104)+0.39*A(112)+0.65*A(113)+0.79*A(115)-A(116)+A(142)+0.4*A(144)+0.55*A(182)+0.6*A(184)&
               &+0.771*A(202)+0.6*A(204)+0.6*A(226)+2*A(241)+A(242)+0.9*A(251)
  Vdot(65) = A(74)-A(143)-A(162)-A(183)-A(203)-A(225)
  Vdot(66) = -A(16)+0.8*A(21)+0.25*A(67)-A(81)+0.12*A(87)+0.41*A(89)+0.03*A(99)+0.03*A(110)+0.16*A(111)+0.53*A(114)&
               &+0.623*A(136)+0.722*A(137)+0.642*A(138)+0.06*A(140)+0.29*A(141)+0.8*A(143)+0.464*A(152)+0.018*A(176)+0.24&
               &*A(177)+0.419*A(178)+0.081*A(180)+0.313*A(181)+A(183)+0.149*A(194)+0.127*A(196)+0.33*A(197)+0.581*A(198)&
               &+0.141*A(200)+0.569*A(201)+A(203)+0.11*A(210)+0.167*A(212)+0.149*A(214)+0.285*A(215)+0.67*A(218)+0.828&
               &*A(219)+0.88*A(220)+0.06*A(222)+0.29*A(223)+A(225)+0.469*A(234)+0.36*A(249)+0.23*A(250)+0.16*A(254)
  Vdot(67) = A(72)+A(73)-A(142)-A(161)-A(182)-A(202)-2*A(241)+A(275)
  Vdot(68) = 0.98*A(128)-A(147)-A(166)-A(187)-A(207)-A(229)
  Vdot(69) = -A(19)+A(63)+A(82)-A(84)+0.15*A(86)-A(97)+0.25*A(99)+0.9*A(116)+0.62*A(117)+0.65*A(145)+0.6*A(146)+A(147)&
               &+0.54*A(150)+0.005*A(176)+0.35*A(185)+0.63*A(186)+A(187)+0.4*A(192)+0.004*A(196)+0.35*A(205)+0.63*A(206)&
               &+A(207)+0.54*A(210)+0.7*A(227)+1.26*A(228)+A(229)+0.54*A(232)+0.5*A(244)+A(246)
  Vdot(70) = A(71)-A(141)-A(160)-A(181)-A(201)-A(223)+0.19*A(249)+A(252)
  Vdot(71) = 0.98*A(122)-A(145)-A(164)-A(185)-A(205)-A(227)
  Vdot(72) = 0.98*A(125)-A(146)-A(165)-A(186)-A(206)-A(228)
  Vdot(73) = -A(10)-A(11)+A(13)+0.13*A(17)+0.45*A(18)+A(22)+A(23)+A(62)+0.01*A(66)-A(79)+0.35*A(88)+0.35*A(90)+A(91)+0.4&
               &*A(92)-A(94)+0.4*A(108)+A(109)+0.64*A(110)+0.02*A(111)+0.9*A(112)+0.58*A(113)+0.04*A(115)+0.7*A(118)+A(134)&
               &+0.047*A(136)+0.021*A(137)+1.6*A(139)+A(140)+A(142)+0.25*A(144)+A(149)+0.287*A(152)+1.33*A(174)+0.75*A(175)&
               &+0.81*A(176)+0.829*A(177)+0.753*A(178)+1.55*A(179)+1.25*A(180)+0.755*A(181)+1.09*A(182)+A(183)+1.4*A(184)&
               &+A(185)+A(186)+A(187)+A(188)+A(189)+2*A(190)+A(191)+0.75*A(192)+0.75*A(193)+0.96*A(194)+0.091*A(196)+0.076&
               &*A(197)+0.8*A(199)+0.501*A(200)+0.34*A(202)+0.4*A(204)+A(209)+0.207*A(212)+0.202*A(214)+0.504*A(215)+A(216)&
               &+0.048*A(218)+0.021*A(219)+1.6*A(221)+A(222)+0.4*A(226)+A(231)+0.28*A(234)+A(236)+A(241)+0.75*A(244)+0.5&
               &*A(246)+0.3*A(249)+0.51*A(250)+0.02*A(254)+A(272)-A(277)
  Vdot(74) = A(20)+0.5*A(86)+0.5*A(99)-A(132)+A(133)-A(149)-A(169)-A(170)-A(190)-A(191)-A(209)-A(231)
  Vdot(75) = -A(21)-A(93)+0.6*A(108)+A(119)+0.059*A(136)+0.124*A(137)+0.261*A(138)+0.2*A(143)+0.35*A(144)+0.05*A(145)&
               &+0.05*A(146)+A(151)+A(172)+A(173)+A(193)+0.5*A(194)+A(211)+0.484*A(212)+2*A(213)+1.5*A(214)+A(215)+A(233)
  Vdot(76) = 0.75*A(67)-A(137)-A(156)-A(177)-A(197)-A(219)
  Vdot(77) = 0.9511*A(68)-A(138)-A(157)-A(178)-A(198)-A(220)
  Vdot(78) = -A(12)+A(14)+0.2*A(21)+0.335*A(66)+0.025*A(68)-A(80)+0.88*A(87)+0.08*A(89)-A(95)+0.25*A(99)+0.44*A(110)&
               &+0.99*A(111)+0.65*A(114)+0.16*A(117)+A(135)+0.233*A(136)+0.211*A(137)+0.15*A(138)+0.2*A(139)+0.94*A(140)&
               &+1.71*A(141)+0.8*A(143)+0.46*A(150)+1.24*A(152)+0.75*A(175)+0.58*A(176)+0.523*A(177)+0.411*A(178)+0.35&
               &*A(179)+0.669*A(180)+0.932*A(181)+A(183)+0.3*A(192)+0.64*A(194)+A(195)+0.724*A(196)+0.677*A(197)+0.497&
               &*A(198)+0.6*A(199)+0.859*A(200)+0.941*A(201)+A(203)+0.35*A(210)+0.65*A(212)+0.64*A(214)+1.21*A(215)+A(217)&
               &+0.243*A(218)+0.239*A(219)+0.187*A(220)+0.2*A(221)+0.94*A(222)+1.71*A(223)+A(225)+0.46*A(232)+1.24*A(234)&
               &+0.85*A(250)+0.99*A(254)+0.991*A(268)+0.11*A(269)+0.11*A(270)+0.11*A(271)+A(273)+A(274)-A(278)
  Vdot(79) = -A(20)-A(86)-A(99)-A(117)+0.5*A(145)+0.95*A(146)+A(185)+A(186)+A(205)+A(206)+0.5*A(227)+A(228)
  Vdot(80) = 0.1*A(76)+0.1*A(77)+0.05*A(78)+0.5*A(86)+0.07*A(89)+0.35*A(90)+A(91)+A(92)+0.5*A(99)+A(108)+0.13*A(112)&
               &+0.048*A(136)+0.334*A(137)+0.416*A(138)+0.16*A(150)+0.085*A(176)+0.245*A(177)+0.322*A(178)+0.08*A(192)+0.071&
               &*A(196)+0.237*A(197)+0.318*A(198)+0.08*A(210)+0.051*A(218)+0.391*A(219)+0.587*A(220)+0.16*A(232)-A(235)&
               &-A(236)-A(237)-2*A(238)-A(239)-A(240)+A(268)+A(269)+A(270)+A(271)+2*A(272)+2*A(273)+1.7*A(274)+A(275)+A(279)&
               &+A(280)
  Vdot(81) = -A(71)-A(103)-A(111)+0.25*A(144)+0.08*A(182)+0.4*A(184)+0.4*A(204)+0.4*A(226)-A(274)
  Vdot(82) = -A(70)-A(102)-A(110)+0.35*A(112)+0.46*A(115)+0.37*A(182)+0.229*A(202)-A(273)+0.3*A(274)
  Vdot(83) = 0.8*A(101)+0.43*A(102)+0.11*A(103)+0.9*A(104)+0.1*A(106)+0.13*A(107)-A(151)-A(172)-A(193)-A(211)-2*A(213)&
               &-A(214)-A(233)+0.9*A(251)+0.11*A(253)
  Vdot(84) = 0.2*A(101)+0.57*A(102)+0.89*A(103)+0.1*A(104)+0.9*A(106)+0.87*A(107)-A(152)-A(173)-A(194)-A(212)-A(214)-2&
               &*A(215)-A(234)+0.1*A(251)+0.89*A(253)
  Vdot(85) = -A(2)-A(3)+A(28)-A(29)-A(32)-A(33)-A(51)-A(52)-A(109)-A(110)-A(111)-A(112)-A(113)-A(114)-A(115)-A(116)&
               &-A(117)-A(118)-A(123)-A(126)-A(129)+A(168)+A(170)-A(250)-A(254)-A(262)
  Vdot(86) = 0.583*A(66)+0.44*A(89)+A(93)-A(136)-A(155)-A(176)-A(196)-A(218)
  Vdot(87) = A(81)+0.03*A(110)+0.12*A(111)+0.02*A(112)+0.42*A(114)+0.42*A(115)-A(150)-A(171)-A(192)-A(210)-A(232)+0.12&
               &*A(254)
  Vdot(88) = 2*A(24)+A(25)+A(26)+A(27)+A(261)-A(262)+1.4*A(263)+A(264)-A(266)-A(267)-A(268)-A(269)-A(270)-A(271)-A(272)&
               &-A(273)-A(274)-A(275)+A(276)-A(277)-A(278)-A(279)-A(280)
  Vdot(89) = 0.65*A(6)+2*A(11)+A(12)+A(13)+A(14)+0.8*A(18)+A(19)+A(20)+A(21)+A(22)+A(23)+A(27)+A(32)-A(33)-A(34)+A(35)-2&
               &*A(36)-2*A(37)+A(43)-A(44)-A(45)+A(46)-A(47)+A(59)+A(60)+A(61)+A(63)+0.381*A(66)+0.25*A(67)+0.049*A(68)+0.1&
               &*A(76)+0.1*A(77)+0.05*A(78)+A(79)+A(82)+A(83)+0.5*A(86)+A(87)+0.35*A(90)+0.4*A(92)+A(94)+A(96)+0.5*A(99)&
               &+0.26*A(109)+0.25*A(110)+0.22*A(111)+0.3*A(112)+0.25*A(113)+0.1*A(114)+0.1*A(115)+0.32*A(116)+0.29*A(117)&
               &+0.08*A(118)-A(120)+0.02*A(122)+0.02*A(125)+0.02*A(128)+A(134)+A(135)+0.742*A(136)+0.599*A(137)+0.606*A(138)&
               &+A(139)+A(140)+A(141)+A(142)+0.8*A(143)+0.65*A(144)+0.95*A(145)+0.95*A(146)+A(147)+0.77*A(150)+A(151)-A(153)&
               &-A(154)-A(155)-A(156)-A(157)-A(158)-A(159)-A(160)-A(161)-A(162)-A(163)-A(164)-A(165)-A(166)-A(167)-A(168)&
               &-A(169)-A(170)-A(171)-A(172)-A(173)+0.66*A(174)+A(175)+0.992*A(176)+0.946*A(177)+0.993*A(178)+A(179)+A(180)&
               &+A(181)+A(182)+2*A(183)+2*A(184)+A(185)+A(186)+2*A(187)+A(188)+A(190)+0.88*A(192)+A(193)+0.5*A(194)+0.5&
               &*A(195)+0.488*A(196)+0.438*A(197)+0.489*A(198)+0.5*A(199)+0.501*A(200)+0.51*A(201)+0.506*A(202)+A(203)&
               &+A(204)+A(205)+A(206)+A(207)+0.38*A(210)+0.5*A(211)+A(213)+0.5*A(214)+A(216)+A(217)+0.792*A(218)+0.699&
               &*A(219)+0.845*A(220)+A(221)+A(222)+A(223)+A(225)+A(226)+A(227)+A(228)+A(229)+0.77*A(232)+A(233)-A(235)&
               &+A(236)+A(241)+0.75*A(244)-A(245)+A(246)+0.22*A(254)-A(265)+A(268)+0.11*A(269)+0.11*A(270)+0.11*A(271)&
               &+A(272)+A(273)+A(274)+A(275)+A(277)+0.88*A(279)+0.84*A(280)
  Vdot(90) = A(12)+A(15)+A(64)+0.65*A(88)+0.19*A(110)+0.23*A(111)+0.03*A(112)+0.08*A(113)-A(134)+0.15*A(136)+0.031&
               &*A(137)+A(148)-A(153)-2*A(174)-A(175)-0.993*A(176)-0.951*A(177)-A(178)-A(179)-A(180)-A(181)-A(182)-A(183)&
               &-A(184)-A(185)-A(186)-A(187)-A(189)-A(190)-A(191)-A(192)-A(193)-A(194)+0.5*A(195)+0.508*A(196)+0.554*A(197)&
               &+0.507*A(198)+0.5*A(199)+0.501*A(200)+0.51*A(201)+0.506*A(202)+A(203)+A(204)+A(205)+A(206)+A(207)+2*A(208)&
               &+A(209)+0.5*A(210)+0.5*A(211)+0.516*A(212)-A(216)+0.155*A(218)+0.04*A(219)+A(230)-A(236)+A(237)+0.23*A(254)&
               &+A(267)
  Vdot(91) = A(4)+A(5)+0.35*A(6)+2*A(9)+A(13)+A(14)+A(15)+A(25)+2*A(31)-A(32)+A(33)-A(34)-A(35)-A(41)-A(42)-A(43)+A(44)&
               &+0.7*A(47)-A(48)-A(49)-A(50)-A(59)-A(60)-A(61)-A(62)-A(63)-A(64)-A(65)-0.964*A(66)-A(67)-A(68)-A(69)-A(70)&
               &-A(71)-A(72)-A(73)-A(74)-A(75)-A(76)-A(77)-A(78)-A(79)-A(80)-A(81)-A(82)-A(83)-A(84)-A(85)-A(86)-A(87)-0.65&
               &*A(88)-0.51*A(89)-A(90)-A(91)-A(92)-A(93)-A(98)+0.12*A(109)+0.4*A(110)+0.63*A(111)+0.28*A(112)+0.25*A(113)&
               &+0.85*A(114)+0.85*A(115)+0.19*A(116)+0.21*A(117)+0.036*A(118)+A(123)+A(126)+A(129)-A(224)-A(243)-A(249)+0.63&
               &*A(250)-A(252)+0.63*A(254)-A(261)-A(276)
  Vdot(92) = -A(14)-A(89)+A(154)+A(155)+A(156)+A(157)+A(158)+A(159)+A(160)+A(162)+A(163)+A(164)+A(165)+A(166)+A(169)&
               &+A(171)+A(235)
  Vdot(93) = A(16)+A(19)+A(22)+A(23)+A(80)+A(84)+0.65*A(90)+A(95)+A(97)+0.15*A(112)+0.1*A(113)+0.1*A(116)+0.28*A(117)&
               &+0.7*A(118)-A(130)+A(131)-A(148)+A(149)+0.23*A(150)-A(167)-A(168)-A(188)-A(189)+A(190)+0.12*A(192)-A(195)&
               &-A(196)-A(197)-A(198)-A(199)-A(200)-A(201)-A(202)-A(203)-A(204)-A(205)-A(206)-A(207)-2*A(208)-0.88*A(210)&
               &-A(211)-A(212)-A(230)+A(231)+0.23*A(232)-A(237)+0.25*A(244)+A(278)
  Vdot(94) = A(1)+A(4)+A(7)-A(38)+A(39)-A(41)-A(44)-A(51)-2*A(53)-A(54)+A(55)-A(134)-A(135)-A(136)-A(137)-A(138)-A(139)&
               &-A(140)-A(141)-A(142)-A(143)-A(144)-A(145)-A(146)-A(147)-A(148)-A(149)-A(150)-A(151)-A(152)-A(239)-A(244)&
               &-A(264)
  Vdot(95) = -A(1)+A(5)+0.65*A(6)+A(8)+A(21)+A(26)+A(38)-A(39)-A(40)-A(42)+A(43)+A(44)-A(45)+A(46)+0.7*A(47)+A(48)+A(50)&
               &+A(51)-A(52)+2*A(53)+2*A(54)-A(56)+A(57)+2*A(58)+A(62)+A(93)+0.5*A(99)+0.4*A(108)+0.7*A(118)-A(119)-A(121)&
               &-A(124)-A(127)-A(130)+A(131)-A(132)+A(133)+A(134)+A(135)+0.941*A(136)+0.876*A(137)+0.739*A(138)+A(139)&
               &+A(140)+A(141)+A(142)+0.8*A(143)+0.65*A(144)+0.95*A(145)+0.95*A(146)+A(147)+A(148)+A(149)+A(150)+A(151)+2&
               &*A(152)+0.5*A(194)+0.516*A(212)+0.5*A(214)+A(215)+A(216)+A(217)+A(218)+A(219)+A(220)+A(221)+A(222)+A(223)&
               &+A(224)+A(225)+A(226)+A(227)+A(228)+A(229)+A(230)+A(231)+A(232)+A(233)+2*A(234)+A(239)+A(240)+A(244)-A(247)&
               &+A(248)+A(264)-A(266)
  Vdot(96) = 0.35*A(6)-A(7)-A(8)+A(40)-A(43)-A(47)+A(49)+A(52)-A(54)-A(55)-A(56)+A(57)-2*A(58)+A(91)+0.6*A(92)-A(94)&
               &-A(95)-A(96)-A(97)-A(99)-A(100)-A(101)-A(102)-A(103)-A(104)-A(105)-A(106)-A(107)-0.4*A(108)-A(216)-A(217)&
               &-A(218)-A(219)-A(220)-A(221)-A(222)-A(223)-A(225)-A(226)-A(227)-A(228)-A(229)-A(230)-A(231)-A(232)-A(233)&
               &-A(234)-A(240)-A(251)-A(253)
  Vdot(97) = A(16)+A(65)+0.1*A(110)+0.18*A(111)+0.2*A(114)+0.16*A(115)-A(135)+0.048*A(136)+0.245*A(137)+0.133*A(138)&
               &-A(154)-A(175)+0.014*A(177)+0.013*A(178)-A(195)+0.006*A(196)+0.018*A(197)+0.015*A(198)-A(217)+0.053*A(218)&
               &+0.262*A(219)+0.155*A(220)+0.18*A(254)
      
END SUBROUTINE racm_soa_vbs_het_Fun
















SUBROUTINE racm_soa_vbs_het_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(95)
  IRR(2) = RCT(2)*V(85)
  IRR(3) = RCT(3)*V(85)
  IRR(4) = RCT(4)*V(29)
  IRR(5) = RCT(5)*V(44)
  IRR(6) = RCT(6)*V(27)
  IRR(7) = RCT(7)*V(96)
  IRR(8) = RCT(8)*V(96)
  IRR(9) = RCT(9)*V(33)
  IRR(10) = RCT(10)*V(73)
  IRR(11) = RCT(11)*V(73)
  IRR(12) = RCT(12)*V(78)
  IRR(13) = RCT(13)*V(25)
  IRR(14) = RCT(14)*V(92)
  IRR(15) = RCT(15)*V(45)
  IRR(16) = RCT(16)*V(66)
  IRR(17) = RCT(17)*V(60)
  IRR(18) = RCT(18)*V(60)
  IRR(19) = RCT(19)*V(69)
  IRR(20) = RCT(20)*V(79)
  IRR(21) = RCT(21)*V(75)
  IRR(22) = RCT(22)*V(64)
  IRR(23) = RCT(23)*V(37)
  IRR(24) = RCT(24)*V(13)
  IRR(25) = RCT(25)*V(19)
  IRR(26) = RCT(26)*V(18)
  IRR(27) = RCT(27)*V(36)
  IRR(28) = RCT(28)*V(30)*F(2)
  IRR(29) = RCT(29)*V(30)*V(85)
  IRR(30) = RCT(30)*V(14)*F(2)
  IRR(31) = RCT(31)*V(14)*F(1)
  IRR(32) = RCT(32)*V(85)*V(91)
  IRR(33) = RCT(33)*V(85)*V(89)
  IRR(34) = RCT(34)*V(89)*V(91)
  IRR(35) = RCT(35)*V(33)*V(91)
  IRR(36) = RCT(36)*V(89)*V(89)
  IRR(37) = RCT(37)*V(89)*V(89)*F(1)
  IRR(38) = RCT(38)*V(30)*V(94)
  IRR(39) = RCT(39)*V(30)*V(95)
  IRR(40) = RCT(40)*V(30)*V(95)
  IRR(41) = RCT(41)*V(91)*V(94)
  IRR(42) = RCT(42)*V(91)*V(95)
  IRR(43) = RCT(43)*V(91)*V(96)
  IRR(44) = RCT(44)*V(89)*V(94)
  IRR(45) = RCT(45)*V(89)*V(95)
  IRR(46) = RCT(46)*V(27)
  IRR(47) = RCT(47)*V(89)*V(96)
  IRR(48) = RCT(48)*V(29)*V(91)
  IRR(49) = RCT(49)*V(44)*V(91)
  IRR(50) = RCT(50)*V(27)*V(91)
  IRR(51) = RCT(51)*V(85)*V(94)
  IRR(52) = RCT(52)*V(85)*V(95)
  IRR(53) = RCT(53)*V(94)*V(94)*F(2)
  IRR(54) = RCT(54)*V(94)*V(96)
  IRR(55) = RCT(55)*V(95)*V(96)
  IRR(56) = RCT(56)*V(95)*V(96)
  IRR(57) = RCT(57)*V(17)
  IRR(58) = RCT(58)*V(96)*V(96)
  IRR(59) = RCT(59)*V(91)*F(2)
  IRR(60) = RCT(60)*V(15)*V(91)
  IRR(61) = RCT(61)*V(52)*V(91)
  IRR(62) = RCT(62)*V(22)*V(91)
  IRR(63) = RCT(63)*V(28)*V(91)
  IRR(64) = RCT(64)*V(35)*V(91)
  IRR(65) = RCT(65)*V(50)*V(91)
  IRR(66) = RCT(66)*V(49)*V(91)
  IRR(67) = RCT(67)*V(41)*V(91)
  IRR(68) = RCT(68)*V(46)*V(91)
  IRR(69) = RCT(69)*V(53)*V(91)
  IRR(70) = RCT(70)*V(82)*V(91)
  IRR(71) = RCT(71)*V(81)*V(91)
  IRR(72) = RCT(72)*V(55)*V(91)
  IRR(73) = RCT(73)*V(59)*V(91)
  IRR(74) = RCT(74)*V(47)*V(91)
  IRR(75) = RCT(75)*V(51)*V(91)
  IRR(76) = RCT(76)*V(23)*V(91)
  IRR(77) = RCT(77)*V(24)*V(91)
  IRR(78) = RCT(78)*V(57)*V(91)
  IRR(79) = RCT(79)*V(73)*V(91)
  IRR(80) = RCT(80)*V(78)*V(91)
  IRR(81) = RCT(81)*V(66)*V(91)
  IRR(82) = RCT(82)*V(37)*V(91)
  IRR(83) = RCT(83)*V(60)*V(91)
  IRR(84) = RCT(84)*V(69)*V(91)
  IRR(85) = RCT(85)*V(64)*V(91)
  IRR(86) = RCT(86)*V(79)*V(91)
  IRR(87) = RCT(87)*V(21)*V(91)
  IRR(88) = RCT(88)*V(25)*V(91)
  IRR(89) = RCT(89)*V(91)*V(92)
  IRR(90) = RCT(90)*V(45)*V(91)
  IRR(91) = RCT(91)*V(48)*V(91)
  IRR(92) = RCT(92)*V(58)*V(91)
  IRR(93) = RCT(93)*V(75)*V(91)
  IRR(94) = RCT(94)*V(73)*V(96)
  IRR(95) = RCT(95)*V(78)*V(96)
  IRR(96) = RCT(96)*V(60)*V(96)
  IRR(97) = RCT(97)*V(69)*V(96)
  IRR(98) = RCT(98)*V(20)*V(91)
  IRR(99) = RCT(99)*V(79)*V(96)
  IRR(100) = RCT(100)*V(57)*V(96)
  IRR(101) = RCT(101)*V(53)*V(96)
  IRR(102) = RCT(102)*V(82)*V(96)
  IRR(103) = RCT(103)*V(81)*V(96)
  IRR(104) = RCT(104)*V(55)*V(96)
  IRR(105) = RCT(105)*V(59)*V(96)
  IRR(106) = RCT(106)*V(47)*V(96)
  IRR(107) = RCT(107)*V(51)*V(96)
  IRR(108) = RCT(108)*V(58)*V(96)
  IRR(109) = RCT(109)*V(53)*V(85)
  IRR(110) = RCT(110)*V(82)*V(85)
  IRR(111) = RCT(111)*V(81)*V(85)
  IRR(112) = RCT(112)*V(55)*V(85)
  IRR(113) = RCT(113)*V(59)*V(85)
  IRR(114) = RCT(114)*V(47)*V(85)
  IRR(115) = RCT(115)*V(51)*V(85)
  IRR(116) = RCT(116)*V(64)*V(85)
  IRR(117) = RCT(117)*V(79)*V(85)
  IRR(118) = RCT(118)*V(58)*V(85)
  IRR(119) = RCT(119)*V(31)*V(95)
  IRR(120) = RCT(120)*V(31)*V(89)
  IRR(121) = RCT(121)*V(38)*V(95)
  IRR(122) = RCT(122)*V(38)*F(2)
  IRR(123) = RCT(123)*V(38)*V(85)
  IRR(124) = RCT(124)*V(39)*V(95)
  IRR(125) = RCT(125)*V(39)*F(2)
  IRR(126) = RCT(126)*V(39)*V(85)
  IRR(127) = RCT(127)*V(40)*V(95)
  IRR(128) = RCT(128)*V(40)*F(2)
  IRR(129) = RCT(129)*V(40)*V(85)
  IRR(130) = RCT(130)*V(93)*V(95)
  IRR(131) = RCT(131)*V(48)
  IRR(132) = RCT(132)*V(74)*V(95)
  IRR(133) = RCT(133)*V(58)
  IRR(134) = RCT(134)*V(90)*V(94)
  IRR(135) = RCT(135)*V(94)*V(97)
  IRR(136) = RCT(136)*V(86)*V(94)
  IRR(137) = RCT(137)*V(76)*V(94)
  IRR(138) = RCT(138)*V(77)*V(94)
  IRR(139) = RCT(139)*V(62)*V(94)
  IRR(140) = RCT(140)*V(61)*V(94)
  IRR(141) = RCT(141)*V(70)*V(94)
  IRR(142) = RCT(142)*V(67)*V(94)
  IRR(143) = RCT(143)*V(65)*V(94)
  IRR(144) = RCT(144)*V(63)*V(94)
  IRR(145) = RCT(145)*V(71)*V(94)
  IRR(146) = RCT(146)*V(72)*V(94)
  IRR(147) = RCT(147)*V(68)*V(94)
  IRR(148) = RCT(148)*V(93)*V(94)
  IRR(149) = RCT(149)*V(74)*V(94)
  IRR(150) = RCT(150)*V(87)*V(94)
  IRR(151) = RCT(151)*V(83)*V(94)
  IRR(152) = RCT(152)*V(84)*V(94)
  IRR(153) = RCT(153)*V(89)*V(90)
  IRR(154) = RCT(154)*V(89)*V(97)
  IRR(155) = RCT(155)*V(86)*V(89)
  IRR(156) = RCT(156)*V(76)*V(89)
  IRR(157) = RCT(157)*V(77)*V(89)
  IRR(158) = RCT(158)*V(62)*V(89)
  IRR(159) = RCT(159)*V(61)*V(89)
  IRR(160) = RCT(160)*V(70)*V(89)
  IRR(161) = RCT(161)*V(67)*V(89)
  IRR(162) = RCT(162)*V(65)*V(89)
  IRR(163) = RCT(163)*V(63)*V(89)
  IRR(164) = RCT(164)*V(71)*V(89)
  IRR(165) = RCT(165)*V(72)*V(89)
  IRR(166) = RCT(166)*V(68)*V(89)
  IRR(167) = RCT(167)*V(89)*V(93)
  IRR(168) = RCT(168)*V(89)*V(93)
  IRR(169) = RCT(169)*V(74)*V(89)
  IRR(170) = RCT(170)*V(74)*V(89)
  IRR(171) = RCT(171)*V(87)*V(89)
  IRR(172) = RCT(172)*V(83)*V(89)
  IRR(173) = RCT(173)*V(84)*V(89)
  IRR(174) = RCT(174)*V(90)*V(90)
  IRR(175) = RCT(175)*V(90)*V(97)
  IRR(176) = RCT(176)*V(86)*V(90)
  IRR(177) = RCT(177)*V(76)*V(90)
  IRR(178) = RCT(178)*V(77)*V(90)
  IRR(179) = RCT(179)*V(62)*V(90)
  IRR(180) = RCT(180)*V(61)*V(90)
  IRR(181) = RCT(181)*V(70)*V(90)
  IRR(182) = RCT(182)*V(67)*V(90)
  IRR(183) = RCT(183)*V(65)*V(90)
  IRR(184) = RCT(184)*V(63)*V(90)
  IRR(185) = RCT(185)*V(71)*V(90)
  IRR(186) = RCT(186)*V(72)*V(90)
  IRR(187) = RCT(187)*V(68)*V(90)
  IRR(188) = RCT(188)*V(90)*V(93)
  IRR(189) = RCT(189)*V(90)*V(93)
  IRR(190) = RCT(190)*V(74)*V(90)
  IRR(191) = RCT(191)*V(74)*V(90)
  IRR(192) = RCT(192)*V(87)*V(90)
  IRR(193) = RCT(193)*V(83)*V(90)
  IRR(194) = RCT(194)*V(84)*V(90)
  IRR(195) = RCT(195)*V(93)*V(97)
  IRR(196) = RCT(196)*V(86)*V(93)
  IRR(197) = RCT(197)*V(76)*V(93)
  IRR(198) = RCT(198)*V(77)*V(93)
  IRR(199) = RCT(199)*V(62)*V(93)
  IRR(200) = RCT(200)*V(61)*V(93)
  IRR(201) = RCT(201)*V(70)*V(93)
  IRR(202) = RCT(202)*V(67)*V(93)
  IRR(203) = RCT(203)*V(65)*V(93)
  IRR(204) = RCT(204)*V(63)*V(93)
  IRR(205) = RCT(205)*V(71)*V(93)
  IRR(206) = RCT(206)*V(72)*V(93)
  IRR(207) = RCT(207)*V(68)*V(93)
  IRR(208) = RCT(208)*V(93)*V(93)
  IRR(209) = RCT(209)*V(74)*V(93)
  IRR(210) = RCT(210)*V(87)*V(93)
  IRR(211) = RCT(211)*V(83)*V(93)
  IRR(212) = RCT(212)*V(84)*V(93)
  IRR(213) = RCT(213)*V(83)*V(83)
  IRR(214) = RCT(214)*V(83)*V(84)
  IRR(215) = RCT(215)*V(84)*V(84)
  IRR(216) = RCT(216)*V(90)*V(96)
  IRR(217) = RCT(217)*V(96)*V(97)
  IRR(218) = RCT(218)*V(86)*V(96)
  IRR(219) = RCT(219)*V(76)*V(96)
  IRR(220) = RCT(220)*V(77)*V(96)
  IRR(221) = RCT(221)*V(62)*V(96)
  IRR(222) = RCT(222)*V(61)*V(96)
  IRR(223) = RCT(223)*V(70)*V(96)
  IRR(224) = RCT(224)*V(26)*V(91)
  IRR(225) = RCT(225)*V(65)*V(96)
  IRR(226) = RCT(226)*V(63)*V(96)
  IRR(227) = RCT(227)*V(71)*V(96)
  IRR(228) = RCT(228)*V(72)*V(96)
  IRR(229) = RCT(229)*V(68)*V(96)
  IRR(230) = RCT(230)*V(93)*V(96)
  IRR(231) = RCT(231)*V(74)*V(96)
  IRR(232) = RCT(232)*V(87)*V(96)
  IRR(233) = RCT(233)*V(83)*V(96)
  IRR(234) = RCT(234)*V(84)*V(96)
  IRR(235) = RCT(235)*V(80)*V(89)
  IRR(236) = RCT(236)*V(80)*V(90)
  IRR(237) = RCT(237)*V(80)*V(93)
  IRR(238) = RCT(238)*V(80)*V(80)
  IRR(239) = RCT(239)*V(80)*V(94)
  IRR(240) = RCT(240)*V(80)*V(96)
  IRR(241) = RCT(241)*V(67)*V(67)
  IRR(242) = RCT(242)*V(16)*V(91)
  IRR(243) = RCT(243)*V(42)*V(91)
  IRR(244) = RCT(244)*V(56)*V(94)
  IRR(245) = RCT(245)*V(56)*V(89)
  IRR(246) = RCT(246)*V(56)*V(56)
  IRR(247) = RCT(247)*V(56)*V(95)
  IRR(248) = RCT(248)*V(26)
  IRR(249) = RCT(249)*V(43)*V(91)
  IRR(250) = RCT(250)*V(43)*V(85)
  IRR(251) = RCT(251)*V(43)*V(96)
  IRR(252) = RCT(252)*V(54)*V(91)
  IRR(253) = RCT(253)*V(54)*V(96)
  IRR(254) = RCT(254)*V(54)*V(85)
  IRR(255) = RCT(255)*V(8)*V(91)
  IRR(256) = RCT(256)*V(7)*V(91)
  IRR(257) = RCT(257)*V(6)*V(91)
  IRR(258) = RCT(258)*V(12)*V(91)
  IRR(259) = RCT(259)*V(11)*V(91)
  IRR(260) = RCT(260)*V(10)*V(91)
  IRR(261) = RCT(261)*V(34)*V(91)
  IRR(262) = RCT(262)*V(85)*V(88)
  IRR(263) = RCT(263)*V(32)*V(32)
  IRR(264) = RCT(264)*V(32)*V(94)
  IRR(265) = RCT(265)*V(32)*V(89)
  IRR(266) = RCT(266)*V(88)*V(95)
  IRR(267) = RCT(267)*V(35)*V(88)
  IRR(268) = RCT(268)*V(50)*V(88)
  IRR(269) = RCT(269)*V(49)*V(88)
  IRR(270) = RCT(270)*V(41)*V(88)
  IRR(271) = RCT(271)*V(46)*V(88)
  IRR(272) = RCT(272)*V(53)*V(88)
  IRR(273) = RCT(273)*V(82)*V(88)
  IRR(274) = RCT(274)*V(81)*V(88)
  IRR(275) = RCT(275)*V(59)*V(88)
  IRR(276) = RCT(276)*V(36)*V(91)
  IRR(277) = RCT(277)*V(73)*V(88)
  IRR(278) = RCT(278)*V(78)*V(88)
  IRR(279) = RCT(279)*V(23)*V(88)
  IRR(280) = RCT(280)*V(24)*V(88)
      
END SUBROUTINE racm_soa_vbs_het_IRRFun
















SUBROUTINE racm_soa_vbs_het_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(518)


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

  B(24) = RCT(24)

  B(25) = RCT(25)

  B(26) = RCT(26)

  B(27) = RCT(27)

  B(28) = RCT(28)*F(2)

  B(30) = RCT(29)*V(85)

  B(31) = RCT(29)*V(30)

  B(32) = RCT(30)*F(2)

  B(34) = RCT(31)*F(1)

  B(36) = RCT(32)*V(91)

  B(37) = RCT(32)*V(85)

  B(38) = RCT(33)*V(89)

  B(39) = RCT(33)*V(85)

  B(40) = RCT(34)*V(91)

  B(41) = RCT(34)*V(89)

  B(42) = RCT(35)*V(91)

  B(43) = RCT(35)*V(33)

  B(44) = RCT(36)*2*V(89)

  B(45) = RCT(37)*2*V(89)*F(1)

  B(47) = RCT(38)*V(94)

  B(48) = RCT(38)*V(30)

  B(49) = RCT(39)*V(95)

  B(50) = RCT(39)*V(30)

  B(51) = RCT(40)*V(95)

  B(52) = RCT(40)*V(30)

  B(53) = RCT(41)*V(94)

  B(54) = RCT(41)*V(91)

  B(55) = RCT(42)*V(95)

  B(56) = RCT(42)*V(91)

  B(57) = RCT(43)*V(96)

  B(58) = RCT(43)*V(91)

  B(59) = RCT(44)*V(94)

  B(60) = RCT(44)*V(89)

  B(61) = RCT(45)*V(95)

  B(62) = RCT(45)*V(89)

  B(63) = RCT(46)

  B(64) = RCT(47)*V(96)

  B(65) = RCT(47)*V(89)

  B(66) = RCT(48)*V(91)

  B(67) = RCT(48)*V(29)

  B(68) = RCT(49)*V(91)

  B(69) = RCT(49)*V(44)

  B(70) = RCT(50)*V(91)

  B(71) = RCT(50)*V(27)

  B(72) = RCT(51)*V(94)

  B(73) = RCT(51)*V(85)

  B(74) = RCT(52)*V(95)

  B(75) = RCT(52)*V(85)

  B(76) = RCT(53)*2*V(94)*F(2)

  B(78) = RCT(54)*V(96)

  B(79) = RCT(54)*V(94)

  B(80) = RCT(55)*V(96)

  B(81) = RCT(55)*V(95)

  B(82) = RCT(56)*V(96)

  B(83) = RCT(56)*V(95)

  B(84) = RCT(57)

  B(85) = RCT(58)*2*V(96)

  B(86) = RCT(59)*F(2)

  B(88) = RCT(60)*V(91)

  B(89) = RCT(60)*V(15)

  B(90) = RCT(61)*V(91)

  B(91) = RCT(61)*V(52)

  B(92) = RCT(62)*V(91)

  B(93) = RCT(62)*V(22)

  B(94) = RCT(63)*V(91)

  B(95) = RCT(63)*V(28)

  B(96) = RCT(64)*V(91)

  B(97) = RCT(64)*V(35)

  B(98) = RCT(65)*V(91)

  B(99) = RCT(65)*V(50)

  B(100) = RCT(66)*V(91)

  B(101) = RCT(66)*V(49)

  B(102) = RCT(67)*V(91)

  B(103) = RCT(67)*V(41)

  B(104) = RCT(68)*V(91)

  B(105) = RCT(68)*V(46)

  B(106) = RCT(69)*V(91)

  B(107) = RCT(69)*V(53)

  B(108) = RCT(70)*V(91)

  B(109) = RCT(70)*V(82)

  B(110) = RCT(71)*V(91)

  B(111) = RCT(71)*V(81)

  B(112) = RCT(72)*V(91)

  B(113) = RCT(72)*V(55)

  B(114) = RCT(73)*V(91)

  B(115) = RCT(73)*V(59)

  B(116) = RCT(74)*V(91)

  B(117) = RCT(74)*V(47)

  B(118) = RCT(75)*V(91)

  B(119) = RCT(75)*V(51)

  B(120) = RCT(76)*V(91)

  B(121) = RCT(76)*V(23)

  B(122) = RCT(77)*V(91)

  B(123) = RCT(77)*V(24)

  B(124) = RCT(78)*V(91)

  B(125) = RCT(78)*V(57)

  B(126) = RCT(79)*V(91)

  B(127) = RCT(79)*V(73)

  B(128) = RCT(80)*V(91)

  B(129) = RCT(80)*V(78)

  B(130) = RCT(81)*V(91)

  B(131) = RCT(81)*V(66)

  B(132) = RCT(82)*V(91)

  B(133) = RCT(82)*V(37)

  B(134) = RCT(83)*V(91)

  B(135) = RCT(83)*V(60)

  B(136) = RCT(84)*V(91)

  B(137) = RCT(84)*V(69)

  B(138) = RCT(85)*V(91)

  B(139) = RCT(85)*V(64)

  B(140) = RCT(86)*V(91)

  B(141) = RCT(86)*V(79)

  B(142) = RCT(87)*V(91)

  B(143) = RCT(87)*V(21)

  B(144) = RCT(88)*V(91)

  B(145) = RCT(88)*V(25)

  B(146) = RCT(89)*V(92)

  B(147) = RCT(89)*V(91)

  B(148) = RCT(90)*V(91)

  B(149) = RCT(90)*V(45)

  B(150) = RCT(91)*V(91)

  B(151) = RCT(91)*V(48)

  B(152) = RCT(92)*V(91)

  B(153) = RCT(92)*V(58)

  B(154) = RCT(93)*V(91)

  B(155) = RCT(93)*V(75)

  B(156) = RCT(94)*V(96)

  B(157) = RCT(94)*V(73)

  B(158) = RCT(95)*V(96)

  B(159) = RCT(95)*V(78)

  B(160) = RCT(96)*V(96)

  B(161) = RCT(96)*V(60)

  B(162) = RCT(97)*V(96)

  B(163) = RCT(97)*V(69)

  B(164) = RCT(98)*V(91)

  B(165) = RCT(98)*V(20)

  B(166) = RCT(99)*V(96)

  B(167) = RCT(99)*V(79)

  B(168) = RCT(100)*V(96)

  B(169) = RCT(100)*V(57)

  B(170) = RCT(101)*V(96)

  B(171) = RCT(101)*V(53)

  B(172) = RCT(102)*V(96)

  B(173) = RCT(102)*V(82)

  B(174) = RCT(103)*V(96)

  B(175) = RCT(103)*V(81)

  B(176) = RCT(104)*V(96)

  B(177) = RCT(104)*V(55)

  B(178) = RCT(105)*V(96)

  B(179) = RCT(105)*V(59)

  B(180) = RCT(106)*V(96)

  B(181) = RCT(106)*V(47)

  B(182) = RCT(107)*V(96)

  B(183) = RCT(107)*V(51)

  B(184) = RCT(108)*V(96)

  B(185) = RCT(108)*V(58)

  B(186) = RCT(109)*V(85)

  B(187) = RCT(109)*V(53)

  B(188) = RCT(110)*V(85)

  B(189) = RCT(110)*V(82)

  B(190) = RCT(111)*V(85)

  B(191) = RCT(111)*V(81)

  B(192) = RCT(112)*V(85)

  B(193) = RCT(112)*V(55)

  B(194) = RCT(113)*V(85)

  B(195) = RCT(113)*V(59)

  B(196) = RCT(114)*V(85)

  B(197) = RCT(114)*V(47)

  B(198) = RCT(115)*V(85)

  B(199) = RCT(115)*V(51)

  B(200) = RCT(116)*V(85)

  B(201) = RCT(116)*V(64)

  B(202) = RCT(117)*V(85)

  B(203) = RCT(117)*V(79)

  B(204) = RCT(118)*V(85)

  B(205) = RCT(118)*V(58)

  B(206) = RCT(119)*V(95)

  B(207) = RCT(119)*V(31)

  B(208) = RCT(120)*V(89)

  B(209) = RCT(120)*V(31)

  B(210) = RCT(121)*V(95)

  B(211) = RCT(121)*V(38)

  B(212) = RCT(122)*F(2)

  B(214) = RCT(123)*V(85)

  B(215) = RCT(123)*V(38)

  B(216) = RCT(124)*V(95)

  B(217) = RCT(124)*V(39)

  B(218) = RCT(125)*F(2)

  B(220) = RCT(126)*V(85)

  B(221) = RCT(126)*V(39)

  B(222) = RCT(127)*V(95)

  B(223) = RCT(127)*V(40)

  B(224) = RCT(128)*F(2)

  B(226) = RCT(129)*V(85)

  B(227) = RCT(129)*V(40)

  B(228) = RCT(130)*V(95)

  B(229) = RCT(130)*V(93)

  B(230) = RCT(131)

  B(231) = RCT(132)*V(95)

  B(232) = RCT(132)*V(74)

  B(233) = RCT(133)

  B(234) = RCT(134)*V(94)

  B(235) = RCT(134)*V(90)

  B(236) = RCT(135)*V(97)

  B(237) = RCT(135)*V(94)

  B(238) = RCT(136)*V(94)

  B(239) = RCT(136)*V(86)

  B(240) = RCT(137)*V(94)

  B(241) = RCT(137)*V(76)

  B(242) = RCT(138)*V(94)

  B(243) = RCT(138)*V(77)

  B(244) = RCT(139)*V(94)

  B(245) = RCT(139)*V(62)

  B(246) = RCT(140)*V(94)

  B(247) = RCT(140)*V(61)

  B(248) = RCT(141)*V(94)

  B(249) = RCT(141)*V(70)

  B(250) = RCT(142)*V(94)

  B(251) = RCT(142)*V(67)

  B(252) = RCT(143)*V(94)

  B(253) = RCT(143)*V(65)

  B(254) = RCT(144)*V(94)

  B(255) = RCT(144)*V(63)

  B(256) = RCT(145)*V(94)

  B(257) = RCT(145)*V(71)

  B(258) = RCT(146)*V(94)

  B(259) = RCT(146)*V(72)

  B(260) = RCT(147)*V(94)

  B(261) = RCT(147)*V(68)

  B(262) = RCT(148)*V(94)

  B(263) = RCT(148)*V(93)

  B(264) = RCT(149)*V(94)

  B(265) = RCT(149)*V(74)

  B(266) = RCT(150)*V(94)

  B(267) = RCT(150)*V(87)

  B(268) = RCT(151)*V(94)

  B(269) = RCT(151)*V(83)

  B(270) = RCT(152)*V(94)

  B(271) = RCT(152)*V(84)

  B(272) = RCT(153)*V(90)

  B(273) = RCT(153)*V(89)

  B(274) = RCT(154)*V(97)

  B(275) = RCT(154)*V(89)

  B(276) = RCT(155)*V(89)

  B(277) = RCT(155)*V(86)

  B(278) = RCT(156)*V(89)

  B(279) = RCT(156)*V(76)

  B(280) = RCT(157)*V(89)

  B(281) = RCT(157)*V(77)

  B(282) = RCT(158)*V(89)

  B(283) = RCT(158)*V(62)

  B(284) = RCT(159)*V(89)

  B(285) = RCT(159)*V(61)

  B(286) = RCT(160)*V(89)

  B(287) = RCT(160)*V(70)

  B(288) = RCT(161)*V(89)

  B(289) = RCT(161)*V(67)

  B(290) = RCT(162)*V(89)

  B(291) = RCT(162)*V(65)

  B(292) = RCT(163)*V(89)

  B(293) = RCT(163)*V(63)

  B(294) = RCT(164)*V(89)

  B(295) = RCT(164)*V(71)

  B(296) = RCT(165)*V(89)

  B(297) = RCT(165)*V(72)

  B(298) = RCT(166)*V(89)

  B(299) = RCT(166)*V(68)

  B(300) = RCT(167)*V(93)

  B(301) = RCT(167)*V(89)

  B(302) = RCT(168)*V(93)

  B(303) = RCT(168)*V(89)

  B(304) = RCT(169)*V(89)

  B(305) = RCT(169)*V(74)

  B(306) = RCT(170)*V(89)

  B(307) = RCT(170)*V(74)

  B(308) = RCT(171)*V(89)

  B(309) = RCT(171)*V(87)

  B(310) = RCT(172)*V(89)

  B(311) = RCT(172)*V(83)

  B(312) = RCT(173)*V(89)

  B(313) = RCT(173)*V(84)

  B(314) = RCT(174)*2*V(90)

  B(315) = RCT(175)*V(97)

  B(316) = RCT(175)*V(90)

  B(317) = RCT(176)*V(90)

  B(318) = RCT(176)*V(86)

  B(319) = RCT(177)*V(90)

  B(320) = RCT(177)*V(76)

  B(321) = RCT(178)*V(90)

  B(322) = RCT(178)*V(77)

  B(323) = RCT(179)*V(90)

  B(324) = RCT(179)*V(62)

  B(325) = RCT(180)*V(90)

  B(326) = RCT(180)*V(61)

  B(327) = RCT(181)*V(90)

  B(328) = RCT(181)*V(70)

  B(329) = RCT(182)*V(90)

  B(330) = RCT(182)*V(67)

  B(331) = RCT(183)*V(90)

  B(332) = RCT(183)*V(65)

  B(333) = RCT(184)*V(90)

  B(334) = RCT(184)*V(63)

  B(335) = RCT(185)*V(90)

  B(336) = RCT(185)*V(71)

  B(337) = RCT(186)*V(90)

  B(338) = RCT(186)*V(72)

  B(339) = RCT(187)*V(90)

  B(340) = RCT(187)*V(68)

  B(341) = RCT(188)*V(93)

  B(342) = RCT(188)*V(90)

  B(343) = RCT(189)*V(93)

  B(344) = RCT(189)*V(90)

  B(345) = RCT(190)*V(90)

  B(346) = RCT(190)*V(74)

  B(347) = RCT(191)*V(90)

  B(348) = RCT(191)*V(74)

  B(349) = RCT(192)*V(90)

  B(350) = RCT(192)*V(87)

  B(351) = RCT(193)*V(90)

  B(352) = RCT(193)*V(83)

  B(353) = RCT(194)*V(90)

  B(354) = RCT(194)*V(84)

  B(355) = RCT(195)*V(97)

  B(356) = RCT(195)*V(93)

  B(357) = RCT(196)*V(93)

  B(358) = RCT(196)*V(86)

  B(359) = RCT(197)*V(93)

  B(360) = RCT(197)*V(76)

  B(361) = RCT(198)*V(93)

  B(362) = RCT(198)*V(77)

  B(363) = RCT(199)*V(93)

  B(364) = RCT(199)*V(62)

  B(365) = RCT(200)*V(93)

  B(366) = RCT(200)*V(61)

  B(367) = RCT(201)*V(93)

  B(368) = RCT(201)*V(70)

  B(369) = RCT(202)*V(93)

  B(370) = RCT(202)*V(67)

  B(371) = RCT(203)*V(93)

  B(372) = RCT(203)*V(65)

  B(373) = RCT(204)*V(93)

  B(374) = RCT(204)*V(63)

  B(375) = RCT(205)*V(93)

  B(376) = RCT(205)*V(71)

  B(377) = RCT(206)*V(93)

  B(378) = RCT(206)*V(72)

  B(379) = RCT(207)*V(93)

  B(380) = RCT(207)*V(68)

  B(381) = RCT(208)*2*V(93)

  B(382) = RCT(209)*V(93)

  B(383) = RCT(209)*V(74)

  B(384) = RCT(210)*V(93)

  B(385) = RCT(210)*V(87)

  B(386) = RCT(211)*V(93)

  B(387) = RCT(211)*V(83)

  B(388) = RCT(212)*V(93)

  B(389) = RCT(212)*V(84)

  B(390) = RCT(213)*2*V(83)

  B(391) = RCT(214)*V(84)

  B(392) = RCT(214)*V(83)

  B(393) = RCT(215)*2*V(84)

  B(394) = RCT(216)*V(96)

  B(395) = RCT(216)*V(90)

  B(396) = RCT(217)*V(97)

  B(397) = RCT(217)*V(96)

  B(398) = RCT(218)*V(96)

  B(399) = RCT(218)*V(86)

  B(400) = RCT(219)*V(96)

  B(401) = RCT(219)*V(76)

  B(402) = RCT(220)*V(96)

  B(403) = RCT(220)*V(77)

  B(404) = RCT(221)*V(96)

  B(405) = RCT(221)*V(62)

  B(406) = RCT(222)*V(96)

  B(407) = RCT(222)*V(61)

  B(408) = RCT(223)*V(96)

  B(409) = RCT(223)*V(70)

  B(410) = RCT(224)*V(91)

  B(411) = RCT(224)*V(26)

  B(412) = RCT(225)*V(96)

  B(413) = RCT(225)*V(65)

  B(414) = RCT(226)*V(96)

  B(415) = RCT(226)*V(63)

  B(416) = RCT(227)*V(96)

  B(417) = RCT(227)*V(71)

  B(418) = RCT(228)*V(96)

  B(419) = RCT(228)*V(72)

  B(420) = RCT(229)*V(96)

  B(421) = RCT(229)*V(68)

  B(422) = RCT(230)*V(96)

  B(423) = RCT(230)*V(93)

  B(424) = RCT(231)*V(96)

  B(425) = RCT(231)*V(74)

  B(426) = RCT(232)*V(96)

  B(427) = RCT(232)*V(87)

  B(428) = RCT(233)*V(96)

  B(429) = RCT(233)*V(83)

  B(430) = RCT(234)*V(96)

  B(431) = RCT(234)*V(84)

  B(432) = RCT(235)*V(89)

  B(433) = RCT(235)*V(80)

  B(434) = RCT(236)*V(90)

  B(435) = RCT(236)*V(80)

  B(436) = RCT(237)*V(93)

  B(437) = RCT(237)*V(80)

  B(438) = RCT(238)*2*V(80)

  B(439) = RCT(239)*V(94)

  B(440) = RCT(239)*V(80)

  B(441) = RCT(240)*V(96)

  B(442) = RCT(240)*V(80)

  B(443) = RCT(241)*2*V(67)

  B(444) = RCT(242)*V(91)

  B(445) = RCT(242)*V(16)

  B(446) = RCT(243)*V(91)

  B(447) = RCT(243)*V(42)

  B(448) = RCT(244)*V(94)

  B(449) = RCT(244)*V(56)

  B(450) = RCT(245)*V(89)

  B(451) = RCT(245)*V(56)

  B(452) = RCT(246)*2*V(56)

  B(453) = RCT(247)*V(95)

  B(454) = RCT(247)*V(56)

  B(455) = RCT(248)

  B(456) = RCT(249)*V(91)

  B(457) = RCT(249)*V(43)

  B(458) = RCT(250)*V(85)

  B(459) = RCT(250)*V(43)

  B(460) = RCT(251)*V(96)

  B(461) = RCT(251)*V(43)

  B(462) = RCT(252)*V(91)

  B(463) = RCT(252)*V(54)

  B(464) = RCT(253)*V(96)

  B(465) = RCT(253)*V(54)

  B(466) = RCT(254)*V(85)

  B(467) = RCT(254)*V(54)

  B(468) = RCT(255)*V(91)

  B(469) = RCT(255)*V(8)

  B(470) = RCT(256)*V(91)

  B(471) = RCT(256)*V(7)

  B(472) = RCT(257)*V(91)

  B(473) = RCT(257)*V(6)

  B(474) = RCT(258)*V(91)

  B(475) = RCT(258)*V(12)

  B(476) = RCT(259)*V(91)

  B(477) = RCT(259)*V(11)

  B(478) = RCT(260)*V(91)

  B(479) = RCT(260)*V(10)

  B(480) = RCT(261)*V(91)

  B(481) = RCT(261)*V(34)

  B(482) = RCT(262)*V(88)

  B(483) = RCT(262)*V(85)

  B(484) = RCT(263)*2*V(32)

  B(485) = RCT(264)*V(94)

  B(486) = RCT(264)*V(32)

  B(487) = RCT(265)*V(89)

  B(488) = RCT(265)*V(32)

  B(489) = RCT(266)*V(95)

  B(490) = RCT(266)*V(88)

  B(491) = RCT(267)*V(88)

  B(492) = RCT(267)*V(35)

  B(493) = RCT(268)*V(88)

  B(494) = RCT(268)*V(50)

  B(495) = RCT(269)*V(88)

  B(496) = RCT(269)*V(49)

  B(497) = RCT(270)*V(88)

  B(498) = RCT(270)*V(41)

  B(499) = RCT(271)*V(88)

  B(500) = RCT(271)*V(46)

  B(501) = RCT(272)*V(88)

  B(502) = RCT(272)*V(53)

  B(503) = RCT(273)*V(88)

  B(504) = RCT(273)*V(82)

  B(505) = RCT(274)*V(88)

  B(506) = RCT(274)*V(81)

  B(507) = RCT(275)*V(88)

  B(508) = RCT(275)*V(59)

  B(509) = RCT(276)*V(91)

  B(510) = RCT(276)*V(36)

  B(511) = RCT(277)*V(88)

  B(512) = RCT(277)*V(73)

  B(513) = RCT(278)*V(88)

  B(514) = RCT(278)*V(78)

  B(515) = RCT(279)*V(88)

  B(516) = RCT(279)*V(23)

  B(517) = RCT(280)*V(88)

  B(518) = RCT(280)*V(24)



  JVS(1) = 0

  JVS(2) = B(88)

  JVS(3) = B(89)

  JVS(4) = 0

  JVS(5) = B(90)

  JVS(6) = B(91)

  JVS(7) = 0

  JVS(8) = 0.05*B(456)+0.039*B(458)

  JVS(9) = 0.036*B(100)

  JVS(10) = 0.01*B(198)

  JVS(11) = 0.37*B(186)

  JVS(12) = 0.15*B(192)

  JVS(13) = 0.11*B(204)

  JVS(14) = 0.28*B(194)

  JVS(15) = 0.45*B(200)

  JVS(16) = 0.11*B(202)

  JVS(17) = 0.14*B(188)

  JVS(18) = 0.37*B(187)+0.14*B(189)+0.15*B(193)+0.28*B(195)+0.01*B(199)+0.45*B(201)+0.11*B(203)+0.11*B(205)+0.039*B(459)

  JVS(19) = 0.036*B(101)+0.05*B(457)

  JVS(20) = 0

  JVS(21) = 0.053*B(458)

  JVS(22) = 0.07*B(198)

  JVS(23) = 0.14*B(466)

  JVS(24) = 0.499*B(365)

  JVS(25) = 0.5*B(363)

  JVS(26) = 0.494*B(369)

  JVS(27) = 0.49*B(367)

  JVS(28) = B(306)+B(347)

  JVS(29) = 0.495*B(359)

  JVS(30) = 0.495*B(361)

  JVS(31) = 0.21*B(202)

  JVS(32) = 0.14*B(190)

  JVS(33) = 0.1*B(188)

  JVS(34) = 0.5*B(386)

  JVS(35) = 0.484*B(388)

  JVS(36) = 0.1*B(189)+0.14*B(191)+0.07*B(199)+0.21*B(203)+0.053*B(459)+0.14*B(467)

  JVS(37) = 0.499*B(357)

  JVS(38) = 0.5*B(384)

  JVS(39) = B(302)+B(307)

  JVS(40) = B(343)+B(348)

  JVS(41) = B(303)+B(344)+0.5*B(355)+0.499*B(358)+0.495*B(360)+0.495*B(362)+0.5*B(364)+0.499*B(366)+0.49*B(368)+0.494&
              &*B(370)+0.5*B(385)+0.5*B(387)+0.484*B(389)

  JVS(42) = 0.5*B(356)

  JVS(43) = 0

  JVS(44) = 1.075*B(472)

  JVS(45) = 1.075*B(473)

  JVS(46) = -B(472)

  JVS(47) = 1.075*B(470)

  JVS(48) = 1.075*B(471)-B(473)

  JVS(49) = -B(470)

  JVS(50) = 1.075*B(468)

  JVS(51) = 1.075*B(469)-B(471)

  JVS(52) = -B(468)

  JVS(53) = -B(469)

  JVS(54) = 0

  JVS(55) = 1.075*B(478)

  JVS(56) = 1.075*B(479)

  JVS(57) = -B(478)

  JVS(58) = 1.075*B(476)

  JVS(59) = 1.075*B(477)-B(479)

  JVS(60) = -B(476)

  JVS(61) = 1.075*B(474)

  JVS(62) = 1.075*B(475)-B(477)

  JVS(63) = -B(474)

  JVS(64) = -B(475)

  JVS(65) = -B(24)

  JVS(66) = 0.3*B(484)

  JVS(67) = -B(32)-B(34)

  JVS(68) = B(2)

  JVS(69) = -B(88)

  JVS(70) = -B(89)

  JVS(71) = -B(444)

  JVS(72) = B(288)

  JVS(73) = B(289)

  JVS(74) = -B(445)

  JVS(75) = -B(84)

  JVS(76) = B(82)

  JVS(77) = B(83)

  JVS(78) = -B(26)

  JVS(79) = B(489)

  JVS(80) = B(490)

  JVS(81) = -B(25)

  JVS(82) = B(487)

  JVS(83) = B(488)

  JVS(84) = -B(164)

  JVS(85) = B(450)

  JVS(86) = B(451)

  JVS(87) = -B(165)

  JVS(88) = -B(142)

  JVS(89) = 0.35*B(140)

  JVS(90) = 0.35*B(141)-B(143)

  JVS(91) = -B(92)

  JVS(92) = B(446)

  JVS(93) = -B(93)+B(447)

  JVS(94) = -B(120)-B(515)

  JVS(95) = -B(516)

  JVS(96) = -B(121)

  JVS(97) = -B(122)-B(517)

  JVS(98) = -B(518)

  JVS(99) = -B(123)

  JVS(100) = -B(13)-B(144)

  JVS(101) = B(272)

  JVS(102) = B(273)

  JVS(103) = -B(145)

  JVS(104) = -B(410)-B(455)

  JVS(105) = B(453)

  JVS(106) = -B(411)

  JVS(107) = B(454)

  JVS(108) = -B(6)-B(63)-B(70)

  JVS(109) = B(61)

  JVS(110) = -B(71)

  JVS(111) = B(62)

  JVS(112) = B(410)

  JVS(113) = -B(94)

  JVS(114) = B(446)

  JVS(115) = 0.25*B(448)+B(452)

  JVS(116) = -B(95)+B(411)+B(447)

  JVS(117) = 0.25*B(449)

  JVS(118) = 0

  JVS(119) = -B(4)-B(66)

  JVS(120) = B(210)

  JVS(121) = B(216)

  JVS(122) = B(222)

  JVS(123) = B(53)-B(67)

  JVS(124) = B(54)

  JVS(125) = B(211)+B(217)+B(223)

  JVS(126) = B(32)

  JVS(127) = -B(28)-B(30)-B(47)-B(49)-B(51)

  JVS(128) = 0.09*B(192)

  JVS(129) = B(3)-B(31)+0.09*B(193)

  JVS(130) = -B(48)

  JVS(131) = B(1)-B(50)-B(52)

  JVS(132) = B(8)

  JVS(133) = -B(206)-B(208)

  JVS(134) = 0.1*B(124)+B(168)

  JVS(135) = -B(209)

  JVS(136) = 0.1*B(125)

  JVS(137) = -B(207)

  JVS(138) = B(169)

  JVS(139) = -2*B(484)-B(485)-B(487)

  JVS(140) = B(482)

  JVS(141) = B(483)

  JVS(142) = -B(488)

  JVS(143) = -B(486)

  JVS(144) = -B(9)-B(42)

  JVS(145) = 0.02*B(196)

  JVS(146) = 0.02*B(198)

  JVS(147) = 0.011*B(466)

  JVS(148) = 0.001*B(192)

  JVS(149) = 0.09*B(194)

  JVS(150) = 0.011*B(190)

  JVS(151) = 0.006*B(188)

  JVS(152) = 0.006*B(189)+0.011*B(191)+0.001*B(193)+0.09*B(195)+0.02*B(197)+0.02*B(199)+0.011*B(467)

  JVS(153) = B(44)+B(45)

  JVS(154) = -B(43)

  JVS(155) = B(515)

  JVS(156) = B(517)

  JVS(157) = -B(480)

  JVS(158) = B(491)

  JVS(159) = B(497)

  JVS(160) = B(499)

  JVS(161) = B(495)

  JVS(162) = B(493)

  JVS(163) = 0.15*B(507)

  JVS(164) = B(511)

  JVS(165) = B(513)

  JVS(166) = 0.3*B(505)

  JVS(167) = B(492)+B(494)+B(496)+B(498)+B(500)+0.3*B(506)+0.15*B(508)+B(512)+B(514)+B(516)+B(518)

  JVS(168) = -B(481)

  JVS(169) = -B(96)-B(491)

  JVS(170) = 0.07*B(466)

  JVS(171) = 0.07*B(190)

  JVS(172) = 0.06*B(188)

  JVS(173) = 0.06*B(189)+0.07*B(191)+0.07*B(467)

  JVS(174) = -B(492)

  JVS(175) = -B(97)

  JVS(176) = -B(27)-B(509)

  JVS(177) = B(501)

  JVS(178) = 0.85*B(507)

  JVS(179) = 0.7*B(505)

  JVS(180) = B(503)

  JVS(181) = B(502)+B(504)+0.7*B(506)+0.85*B(508)

  JVS(182) = -B(510)

  JVS(183) = -B(23)-B(132)

  JVS(184) = 0.024*B(104)

  JVS(185) = 0.6*B(152)

  JVS(186) = 0.3*B(349)

  JVS(187) = 0.3*B(350)

  JVS(188) = 0.024*B(105)-B(133)+0.6*B(153)

  JVS(189) = 0.9*B(120)

  JVS(190) = -B(210)-B(212)-B(214)

  JVS(191) = -B(215)

  JVS(192) = 0

  JVS(193) = 0.9*B(121)

  JVS(194) = -B(211)

  JVS(195) = 0.9*B(122)

  JVS(196) = -B(216)-B(218)-B(220)

  JVS(197) = -B(221)

  JVS(198) = 0

  JVS(199) = 0.9*B(123)

  JVS(200) = -B(217)

  JVS(201) = -B(222)-B(224)-B(226)

  JVS(202) = 0.85*B(124)

  JVS(203) = -B(227)

  JVS(204) = 0.85*B(125)

  JVS(205) = -B(223)

  JVS(206) = -B(102)-1.11*B(497)

  JVS(207) = 0.1*B(505)

  JVS(208) = -0.33*B(503)

  JVS(209) = -1.11*B(498)-0.33*B(504)+0.1*B(506)

  JVS(210) = -B(103)

  JVS(211) = -B(446)

  JVS(212) = B(178)

  JVS(213) = 0.046*B(250)

  JVS(214) = -B(447)

  JVS(215) = 0.046*B(251)

  JVS(216) = B(179)

  JVS(217) = -B(456)-B(458)-B(460)

  JVS(218) = -B(459)

  JVS(219) = -B(457)

  JVS(220) = -B(461)

  JVS(221) = -B(5)-B(68)

  JVS(222) = B(168)

  JVS(223) = B(160)

  JVS(224) = B(162)

  JVS(225) = B(156)

  JVS(226) = B(158)

  JVS(227) = 0.5*B(166)

  JVS(228) = 0.3*B(64)

  JVS(229) = B(55)-B(69)

  JVS(230) = B(56)

  JVS(231) = 0.3*B(65)+B(157)+B(159)+B(161)+B(163)+0.5*B(167)+B(169)

  JVS(232) = -B(15)-B(148)

  JVS(233) = 0.11*B(202)

  JVS(234) = 0.11*B(203)

  JVS(235) = B(300)

  JVS(236) = -B(149)

  JVS(237) = B(301)

  JVS(238) = -B(104)-1.11*B(499)

  JVS(239) = 0.1*B(505)

  JVS(240) = -0.33*B(503)

  JVS(241) = -1.11*B(500)-0.33*B(504)+0.1*B(506)

  JVS(242) = -B(105)

  JVS(243) = -B(116)-B(180)-B(196)

  JVS(244) = -B(197)

  JVS(245) = -B(117)

  JVS(246) = -B(181)

  JVS(247) = -B(150)-B(230)

  JVS(248) = 0.4*B(152)+0.4*B(184)+0.3*B(204)

  JVS(249) = 0.3*B(205)

  JVS(250) = -B(151)+0.4*B(153)

  JVS(251) = B(228)

  JVS(252) = B(229)

  JVS(253) = 0.4*B(185)

  JVS(254) = -B(100)-1.11*B(495)

  JVS(255) = 0.1*B(505)

  JVS(256) = -0.33*B(503)

  JVS(257) = -1.11*B(496)-0.33*B(504)+0.1*B(506)

  JVS(258) = -B(101)

  JVS(259) = -B(98)-B(493)

  JVS(260) = 0.06*B(466)

  JVS(261) = 0.06*B(190)

  JVS(262) = 0.03*B(188)

  JVS(263) = 0.03*B(189)+0.06*B(191)+0.06*B(467)

  JVS(264) = -B(494)

  JVS(265) = -B(99)

  JVS(266) = -B(118)-B(182)-B(198)

  JVS(267) = -B(199)

  JVS(268) = -B(119)

  JVS(269) = -B(183)

  JVS(270) = B(92)

  JVS(271) = B(27)+B(509)

  JVS(272) = 0

  JVS(273) = 0.14*B(196)

  JVS(274) = 0.036*B(100)

  JVS(275) = 0.14*B(198)

  JVS(276) = -B(90)

  JVS(277) = 0.43*B(186)

  JVS(278) = 0.3*B(466)

  JVS(279) = 0.36*B(192)

  JVS(280) = 0.25*B(448)+0.5*B(452)

  JVS(281) = 0.13*B(204)

  JVS(282) = 0.14*B(194)

  JVS(283) = 1.87*B(17)+1.55*B(18)+2*B(134)+2*B(160)

  JVS(284) = B(22)+0.22*B(200)

  JVS(285) = 0

  JVS(286) = B(19)+B(136)+B(162)

  JVS(287) = B(10)+B(11)+B(126)+B(156)+B(511)

  JVS(288) = B(12)

  JVS(289) = 0.66*B(202)

  JVS(290) = 0.3*B(190)

  JVS(291) = 0.37*B(188)

  JVS(292) = 0.43*B(187)+0.37*B(189)+0.3*B(191)+0.36*B(193)+0.14*B(195)+0.14*B(197)+0.14*B(199)+0.22*B(201)+0.66*B(203)&
               &+0.13*B(205)+0.3*B(467)

  JVS(293) = B(512)

  JVS(294) = -B(91)+B(93)+0.036*B(101)+B(127)+2*B(135)+B(137)+B(510)

  JVS(295) = 0.25*B(449)

  JVS(296) = B(157)+2*B(161)+B(163)

  JVS(297) = -B(106)-B(170)-B(186)-B(501)

  JVS(298) = -B(187)

  JVS(299) = -B(502)

  JVS(300) = -B(107)

  JVS(301) = -B(171)

  JVS(302) = -B(462)-B(464)-B(466)

  JVS(303) = -B(467)

  JVS(304) = -B(463)

  JVS(305) = -B(465)

  JVS(306) = -B(112)-B(176)-B(192)

  JVS(307) = -B(193)

  JVS(308) = -B(113)

  JVS(309) = -B(177)

  JVS(310) = B(164)

  JVS(311) = B(455)

  JVS(312) = -B(448)-B(450)-2*B(452)-B(453)

  JVS(313) = 0.1*B(194)

  JVS(314) = B(138)

  JVS(315) = 0.1*B(195)

  JVS(316) = -B(451)

  JVS(317) = B(139)+B(165)

  JVS(318) = -B(449)

  JVS(319) = -B(454)

  JVS(320) = 0.1*B(206)+B(208)

  JVS(321) = B(210)+0.02*B(212)+B(214)

  JVS(322) = B(216)+0.02*B(218)+B(220)

  JVS(323) = B(222)+0.02*B(224)+B(226)

  JVS(324) = -B(124)-B(168)

  JVS(325) = B(215)+B(221)+B(227)

  JVS(326) = 0

  JVS(327) = B(209)

  JVS(328) = -B(125)

  JVS(329) = 0.1*B(207)+B(211)+B(217)+B(223)

  JVS(330) = -B(169)

  JVS(331) = -B(152)-B(184)-B(204)-B(233)

  JVS(332) = B(231)

  JVS(333) = -B(205)

  JVS(334) = -B(153)

  JVS(335) = B(232)

  JVS(336) = -B(185)

  JVS(337) = -B(114)-B(178)-B(194)-B(507)

  JVS(338) = -B(195)

  JVS(339) = -B(508)

  JVS(340) = -B(115)

  JVS(341) = -B(179)

  JVS(342) = 0.036*B(100)

  JVS(343) = -B(17)-B(18)-B(134)-B(160)

  JVS(344) = B(260)+B(339)+B(379)+B(420)

  JVS(345) = 1.2*B(256)+0.65*B(335)+0.65*B(375)+1.3*B(416)

  JVS(346) = 0.35*B(258)+0.37*B(337)+0.37*B(377)+0.74*B(418)

  JVS(347) = 0.15*B(140)+0.25*B(166)+0.5*B(202)

  JVS(348) = 0

  JVS(349) = 0

  JVS(350) = 0.5*B(203)

  JVS(351) = 0.063*B(238)+0.119*B(317)+0.1*B(357)+0.063*B(398)

  JVS(352) = 0

  JVS(353) = 0.119*B(318)+0.65*B(336)+0.37*B(338)+B(340)

  JVS(354) = 0.036*B(101)-B(135)+0.15*B(141)

  JVS(355) = 0.1*B(358)+0.65*B(376)+0.37*B(378)+B(380)

  JVS(356) = 0.063*B(239)+1.2*B(257)+0.35*B(259)+B(261)

  JVS(357) = -B(161)+0.25*B(167)+0.063*B(399)+1.3*B(417)+0.74*B(419)+B(421)

  JVS(358) = -B(246)-B(284)-B(325)-B(365)-B(406)

  JVS(359) = B(108)

  JVS(360) = -B(285)

  JVS(361) = -B(326)

  JVS(362) = B(109)

  JVS(363) = -B(366)

  JVS(364) = -B(247)

  JVS(365) = -B(407)

  JVS(366) = B(106)

  JVS(367) = -B(244)-B(282)-B(323)-B(363)-B(404)

  JVS(368) = 0

  JVS(369) = 0

  JVS(370) = -B(283)

  JVS(371) = -B(324)

  JVS(372) = B(107)

  JVS(373) = -B(364)

  JVS(374) = -B(245)

  JVS(375) = -B(405)

  JVS(376) = B(118)

  JVS(377) = -B(254)-B(292)-B(333)-B(373)-B(414)

  JVS(378) = 0

  JVS(379) = -B(293)

  JVS(380) = -B(334)

  JVS(381) = B(119)

  JVS(382) = -B(374)

  JVS(383) = -B(255)

  JVS(384) = -B(415)

  JVS(385) = B(444)

  JVS(386) = 0.9*B(460)

  JVS(387) = 0.79*B(198)

  JVS(388) = 0.9*B(176)+0.39*B(192)

  JVS(389) = 0.65*B(194)

  JVS(390) = 0.4*B(254)+0.6*B(333)+0.6*B(373)+0.6*B(414)

  JVS(391) = -B(22)-B(138)-B(200)

  JVS(392) = B(250)+0.55*B(329)+0.771*B(369)+2*B(443)

  JVS(393) = 0.39*B(193)+0.65*B(195)+0.79*B(199)-B(201)

  JVS(394) = 0

  JVS(395) = 0

  JVS(396) = 0.55*B(330)+0.6*B(334)

  JVS(397) = -B(139)+B(445)

  JVS(398) = 0.771*B(370)+0.6*B(374)

  JVS(399) = B(251)+0.4*B(255)

  JVS(400) = 0.9*B(177)+0.6*B(415)+0.9*B(461)

  JVS(401) = B(116)

  JVS(402) = -B(252)-B(290)-B(331)-B(371)-B(412)

  JVS(403) = 0

  JVS(404) = -B(291)

  JVS(405) = -B(332)

  JVS(406) = B(117)

  JVS(407) = -B(372)

  JVS(408) = -B(253)

  JVS(409) = -B(413)

  JVS(410) = 0.12*B(142)

  JVS(411) = 0.25*B(102)

  JVS(412) = 0.36*B(456)+0.23*B(458)

  JVS(413) = 0.53*B(196)

  JVS(414) = 0.16*B(466)

  JVS(415) = 0.06*B(246)+0.081*B(325)+0.141*B(365)+0.06*B(406)

  JVS(416) = 0.8*B(252)+B(331)+B(371)+B(412)

  JVS(417) = -B(16)-B(130)

  JVS(418) = 0.29*B(248)+0.313*B(327)+0.569*B(367)+0.29*B(408)

  JVS(419) = 0.8*B(21)

  JVS(420) = 0.722*B(240)+0.24*B(319)+0.33*B(359)+0.828*B(400)

  JVS(421) = 0.642*B(242)+0.419*B(321)+0.581*B(361)+0.88*B(402)

  JVS(422) = 0.03*B(166)

  JVS(423) = 0.16*B(190)

  JVS(424) = 0.03*B(188)

  JVS(425) = 0.149*B(391)

  JVS(426) = 0.464*B(270)+0.149*B(353)+0.167*B(388)+0.149*B(392)+0.285*B(393)+0.469*B(430)

  JVS(427) = 0.03*B(189)+0.16*B(191)+0.53*B(197)+0.23*B(459)+0.16*B(467)

  JVS(428) = 0.623*B(238)+0.018*B(317)+0.127*B(357)+0.67*B(398)

  JVS(429) = 0.11*B(384)

  JVS(430) = 0

  JVS(431) = 0

  JVS(432) = 0.018*B(318)+0.24*B(320)+0.419*B(322)+0.081*B(326)+0.313*B(328)+B(332)+0.149*B(354)

  JVS(433) = 0.25*B(103)-B(131)+0.12*B(143)+0.41*B(146)+0.36*B(457)

  JVS(434) = 0.41*B(147)

  JVS(435) = 0.127*B(358)+0.33*B(360)+0.581*B(362)+0.141*B(366)+0.569*B(368)+B(372)+0.11*B(385)+0.167*B(389)

  JVS(436) = 0.623*B(239)+0.722*B(241)+0.642*B(243)+0.06*B(247)+0.29*B(249)+0.8*B(253)+0.464*B(271)

  JVS(437) = 0.03*B(167)+0.67*B(399)+0.828*B(401)+0.88*B(403)+0.06*B(407)+0.29*B(409)+B(413)+0.469*B(431)

  JVS(438) = B(112)

  JVS(439) = B(114)+B(507)

  JVS(440) = -B(250)-B(288)-B(329)-B(369)-2*B(443)

  JVS(441) = 0

  JVS(442) = B(508)

  JVS(443) = -B(289)

  JVS(444) = -B(330)

  JVS(445) = B(113)+B(115)

  JVS(446) = -B(370)

  JVS(447) = -B(251)

  JVS(448) = 0

  JVS(449) = 0.98*B(224)

  JVS(450) = 0

  JVS(451) = -B(260)-B(298)-B(339)-B(379)-B(420)

  JVS(452) = 0

  JVS(453) = 0

  JVS(454) = -B(299)

  JVS(455) = -B(340)

  JVS(456) = 0

  JVS(457) = -B(380)

  JVS(458) = -B(261)

  JVS(459) = 0

  JVS(460) = -B(421)

  JVS(461) = B(94)

  JVS(462) = B(132)

  JVS(463) = 0

  JVS(464) = 0

  JVS(465) = 0.5*B(448)+B(452)

  JVS(466) = 0

  JVS(467) = 0

  JVS(468) = 0.9*B(200)

  JVS(469) = 0

  JVS(470) = B(260)+B(339)+B(379)+B(420)

  JVS(471) = -B(19)-B(136)-B(162)

  JVS(472) = 0.65*B(256)+0.35*B(335)+0.35*B(375)+0.7*B(416)

  JVS(473) = 0.6*B(258)+0.63*B(337)+0.63*B(377)+1.26*B(418)

  JVS(474) = 0

  JVS(475) = 0.15*B(140)+0.25*B(166)+0.62*B(202)

  JVS(476) = 0

  JVS(477) = 0

  JVS(478) = 0.9*B(201)+0.62*B(203)

  JVS(479) = 0.005*B(317)+0.004*B(357)

  JVS(480) = 0.54*B(266)+0.4*B(349)+0.54*B(384)+0.54*B(426)

  JVS(481) = 0

  JVS(482) = 0

  JVS(483) = 0.005*B(318)+0.35*B(336)+0.63*B(338)+B(340)+0.4*B(350)

  JVS(484) = B(95)+B(133)-B(137)+0.15*B(141)

  JVS(485) = 0.004*B(358)+0.35*B(376)+0.63*B(378)+B(380)+0.54*B(385)

  JVS(486) = 0.65*B(257)+0.6*B(259)+B(261)+0.54*B(267)+0.5*B(449)

  JVS(487) = 0

  JVS(488) = -B(163)+0.25*B(167)+0.7*B(417)+1.26*B(419)+B(421)+0.54*B(427)

  JVS(489) = 0.19*B(456)

  JVS(490) = B(462)

  JVS(491) = -B(248)-B(286)-B(327)-B(367)-B(408)

  JVS(492) = B(110)

  JVS(493) = 0

  JVS(494) = -B(287)

  JVS(495) = -B(328)

  JVS(496) = B(111)+0.19*B(457)+B(463)

  JVS(497) = -B(368)

  JVS(498) = -B(249)

  JVS(499) = -B(409)

  JVS(500) = 0.98*B(212)

  JVS(501) = -B(256)-B(294)-B(335)-B(375)-B(416)

  JVS(502) = 0

  JVS(503) = 0

  JVS(504) = -B(295)

  JVS(505) = -B(336)

  JVS(506) = 0

  JVS(507) = -B(376)

  JVS(508) = -B(257)

  JVS(509) = 0

  JVS(510) = -B(417)

  JVS(511) = 0.98*B(218)

  JVS(512) = -B(258)-B(296)-B(337)-B(377)-B(418)

  JVS(513) = 0

  JVS(514) = 0

  JVS(515) = -B(297)

  JVS(516) = -B(338)

  JVS(517) = 0

  JVS(518) = -B(378)

  JVS(519) = -B(259)

  JVS(520) = 0

  JVS(521) = -B(419)

  JVS(522) = B(92)

  JVS(523) = B(13)+0.35*B(144)

  JVS(524) = B(23)

  JVS(525) = 0

  JVS(526) = 0.3*B(456)+0.51*B(458)

  JVS(527) = 0.35*B(148)

  JVS(528) = 0

  JVS(529) = B(150)

  JVS(530) = 0.01*B(100)

  JVS(531) = 0.04*B(198)

  JVS(532) = B(186)+B(501)

  JVS(533) = 0.02*B(466)

  JVS(534) = 0.9*B(192)

  JVS(535) = 0.75*B(448)+0.5*B(452)

  JVS(536) = 0.4*B(152)+0.4*B(184)+0.7*B(204)

  JVS(537) = 0.58*B(194)

  JVS(538) = 0.13*B(17)+0.45*B(18)

  JVS(539) = B(246)+1.25*B(325)+0.501*B(365)+B(406)

  JVS(540) = 1.6*B(244)+1.55*B(323)+0.8*B(363)+1.6*B(404)

  JVS(541) = 0.25*B(254)+1.4*B(333)+0.4*B(373)+0.4*B(414)

  JVS(542) = B(22)

  JVS(543) = B(331)

  JVS(544) = B(250)+1.09*B(329)+0.34*B(369)+B(443)

  JVS(545) = B(339)

  JVS(546) = 0.755*B(327)

  JVS(547) = B(335)

  JVS(548) = B(337)

  JVS(549) = -B(10)-B(11)-B(126)-B(156)-B(511)

  JVS(550) = B(264)+2*B(345)+B(347)+B(382)+B(424)

  JVS(551) = 0.021*B(240)+0.829*B(319)+0.076*B(359)+0.021*B(400)

  JVS(552) = 0.753*B(321)

  JVS(553) = 0

  JVS(554) = B(434)

  JVS(555) = 0.02*B(190)

  JVS(556) = 0.64*B(188)

  JVS(557) = 0.75*B(351)+0.202*B(391)

  JVS(558) = 0.287*B(270)+0.96*B(353)+0.207*B(388)+0.202*B(392)+0.504*B(393)+0.28*B(430)

  JVS(559) = B(187)+0.64*B(189)+0.02*B(191)+0.9*B(193)+0.58*B(195)+0.04*B(199)+0.7*B(205)+0.51*B(459)+0.02*B(467)

  JVS(560) = 0.047*B(238)+0.81*B(317)+0.091*B(357)+0.048*B(398)

  JVS(561) = 0.75*B(349)

  JVS(562) = B(502)-B(512)

  JVS(563) = 0

  JVS(564) = B(234)+1.33*B(314)+0.75*B(315)+0.81*B(318)+0.829*B(320)+0.753*B(322)+1.55*B(324)+1.25*B(326)+0.755*B(328)&
               &+1.09*B(330)+B(332)+1.4*B(334)+B(336)+B(338)+B(340)+B(341)+B(343)+2*B(346)+B(348)+0.75*B(350)+0.75*B(352)&
               &+0.96*B(354)+B(394)+B(435)

  JVS(565) = B(93)+0.01*B(101)-B(127)+0.35*B(145)+0.35*B(149)+B(151)+0.4*B(153)+0.3*B(457)

  JVS(566) = B(342)+B(344)+0.091*B(358)+0.076*B(360)+0.8*B(364)+0.501*B(366)+0.34*B(370)+0.4*B(374)+B(383)+0.207*B(389)

  JVS(567) = B(235)+0.047*B(239)+0.021*B(241)+1.6*B(245)+B(247)+B(251)+0.25*B(255)+B(265)+0.287*B(271)+0.75*B(449)

  JVS(568) = 0

  JVS(569) = -B(157)+0.4*B(185)+B(395)+0.048*B(399)+0.021*B(401)+1.6*B(405)+B(407)+0.4*B(415)+B(425)+0.28*B(431)

  JVS(570) = 0.75*B(316)

  JVS(571) = B(233)

  JVS(572) = -B(231)-B(264)-B(304)-B(306)-B(345)-B(347)-B(382)-B(424)

  JVS(573) = B(20)+0.5*B(140)+0.5*B(166)

  JVS(574) = 0

  JVS(575) = -B(305)-B(307)

  JVS(576) = -B(346)-B(348)

  JVS(577) = 0.5*B(141)

  JVS(578) = -B(383)

  JVS(579) = -B(265)

  JVS(580) = -B(232)

  JVS(581) = 0.5*B(167)-B(425)

  JVS(582) = B(206)

  JVS(583) = 0

  JVS(584) = 0.6*B(184)

  JVS(585) = 0.35*B(254)

  JVS(586) = 0.2*B(252)

  JVS(587) = 0.05*B(256)

  JVS(588) = 0.05*B(258)

  JVS(589) = 0

  JVS(590) = -B(21)-B(154)

  JVS(591) = 0.124*B(240)

  JVS(592) = 0.261*B(242)

  JVS(593) = 0

  JVS(594) = B(268)+B(310)+B(351)+B(386)+2*B(390)+1.5*B(391)+B(428)

  JVS(595) = B(312)+0.5*B(353)+0.484*B(388)+1.5*B(392)+B(393)

  JVS(596) = 0

  JVS(597) = 0.059*B(238)

  JVS(598) = 0

  JVS(599) = B(311)+B(313)

  JVS(600) = B(352)+0.5*B(354)

  JVS(601) = -B(155)

  JVS(602) = B(387)+0.484*B(389)

  JVS(603) = 0.059*B(239)+0.124*B(241)+0.261*B(243)+0.2*B(253)+0.35*B(255)+0.05*B(257)+0.05*B(259)+B(269)

  JVS(604) = B(207)

  JVS(605) = 0.6*B(185)+B(429)

  JVS(606) = 0.75*B(102)

  JVS(607) = -B(240)-B(278)-B(319)-B(359)-B(400)

  JVS(608) = 0

  JVS(609) = 0

  JVS(610) = 0

  JVS(611) = -B(279)

  JVS(612) = -B(320)

  JVS(613) = 0.75*B(103)

  JVS(614) = -B(360)

  JVS(615) = -B(241)

  JVS(616) = -B(401)

  JVS(617) = 0.9511*B(104)

  JVS(618) = -B(242)-B(280)-B(321)-B(361)-B(402)

  JVS(619) = 0

  JVS(620) = 0

  JVS(621) = 0

  JVS(622) = -B(281)

  JVS(623) = -B(322)

  JVS(624) = 0.9511*B(105)

  JVS(625) = -B(362)

  JVS(626) = -B(243)

  JVS(627) = -B(403)

  JVS(628) = 0.88*B(142)

  JVS(629) = 0.11*B(497)

  JVS(630) = 0.85*B(458)

  JVS(631) = 0.025*B(104)+0.11*B(499)

  JVS(632) = 0.65*B(196)

  JVS(633) = 0.335*B(100)+0.11*B(495)

  JVS(634) = 0.991*B(493)

  JVS(635) = 0.99*B(466)

  JVS(636) = 0.94*B(246)+0.669*B(325)+0.859*B(365)+0.94*B(406)

  JVS(637) = 0.2*B(244)+0.35*B(323)+0.6*B(363)+0.2*B(404)

  JVS(638) = 0.8*B(252)+B(331)+B(371)+B(412)

  JVS(639) = 1.71*B(248)+0.932*B(327)+0.941*B(367)+1.71*B(408)

  JVS(640) = 0.2*B(21)

  JVS(641) = 0.211*B(240)+0.523*B(319)+0.677*B(359)+0.239*B(400)

  JVS(642) = 0.15*B(242)+0.411*B(321)+0.497*B(361)+0.187*B(402)

  JVS(643) = -B(12)-B(128)-B(158)-B(513)

  JVS(644) = 0.25*B(166)+0.16*B(202)

  JVS(645) = 0.99*B(190)+B(505)

  JVS(646) = 0.44*B(188)+B(503)

  JVS(647) = 0.64*B(391)

  JVS(648) = 1.24*B(270)+0.64*B(353)+0.65*B(388)+0.64*B(392)+1.21*B(393)+1.24*B(430)

  JVS(649) = 0.44*B(189)+0.99*B(191)+0.65*B(197)+0.16*B(203)+0.85*B(459)+0.99*B(467)

  JVS(650) = 0.233*B(238)+0.58*B(317)+0.724*B(357)+0.243*B(398)

  JVS(651) = 0.46*B(266)+0.3*B(349)+0.35*B(384)+0.46*B(426)

  JVS(652) = 0.991*B(494)+0.11*B(496)+0.11*B(498)+0.11*B(500)+B(504)+B(506)-B(514)

  JVS(653) = 0

  JVS(654) = 0.75*B(315)+0.58*B(318)+0.523*B(320)+0.411*B(322)+0.35*B(324)+0.669*B(326)+0.932*B(328)+B(332)+0.3*B(350)&
               &+0.64*B(354)

  JVS(655) = 0.335*B(101)+0.025*B(105)-B(129)+0.88*B(143)+0.08*B(146)

  JVS(656) = B(14)+0.08*B(147)

  JVS(657) = B(355)+0.724*B(358)+0.677*B(360)+0.497*B(362)+0.6*B(364)+0.859*B(366)+0.941*B(368)+B(372)+0.35*B(385)+0.65&
               &*B(389)

  JVS(658) = B(236)+0.233*B(239)+0.211*B(241)+0.15*B(243)+0.2*B(245)+0.94*B(247)+1.71*B(249)+0.8*B(253)+0.46*B(267)+1.24&
               &*B(271)

  JVS(659) = 0

  JVS(660) = -B(159)+0.25*B(167)+B(396)+0.243*B(399)+0.239*B(401)+0.187*B(403)+0.2*B(405)+0.94*B(407)+1.71*B(409)+B(413)&
               &+0.46*B(427)+1.24*B(431)

  JVS(661) = B(237)+0.75*B(316)+B(356)+B(397)

  JVS(662) = 0.5*B(256)+B(335)+B(375)+0.5*B(416)

  JVS(663) = 0.95*B(258)+B(337)+B(377)+B(418)

  JVS(664) = -B(20)-B(140)-B(166)-B(202)

  JVS(665) = -B(203)

  JVS(666) = 0

  JVS(667) = 0

  JVS(668) = B(336)+B(338)

  JVS(669) = -B(141)

  JVS(670) = B(376)+B(378)

  JVS(671) = 0.5*B(257)+0.95*B(259)

  JVS(672) = 0

  JVS(673) = -B(167)+0.5*B(417)+B(419)

  JVS(674) = 0.1*B(120)+B(515)

  JVS(675) = 0.1*B(122)+B(517)

  JVS(676) = B(497)

  JVS(677) = 0.35*B(148)

  JVS(678) = B(499)

  JVS(679) = B(150)

  JVS(680) = B(495)

  JVS(681) = B(493)

  JVS(682) = 2*B(501)

  JVS(683) = 0

  JVS(684) = 0.13*B(192)

  JVS(685) = 0.05*B(124)

  JVS(686) = B(152)+B(184)

  JVS(687) = B(507)

  JVS(688) = 0

  JVS(689) = 0.334*B(240)+0.245*B(319)+0.237*B(359)+0.391*B(400)

  JVS(690) = 0.416*B(242)+0.322*B(321)+0.318*B(361)+0.587*B(402)

  JVS(691) = 0.5*B(140)+0.5*B(166)

  JVS(692) = -B(432)-B(434)-B(436)-2*B(438)-B(439)-B(441)

  JVS(693) = 1.7*B(505)

  JVS(694) = 2*B(503)

  JVS(695) = 0.13*B(193)

  JVS(696) = 0.048*B(238)+0.085*B(317)+0.071*B(357)+0.051*B(398)

  JVS(697) = 0.16*B(266)+0.08*B(349)+0.08*B(384)+0.16*B(426)

  JVS(698) = B(494)+B(496)+B(498)+B(500)+2*B(502)+2*B(504)+1.7*B(506)+B(508)+B(516)+B(518)

  JVS(699) = -B(433)

  JVS(700) = 0.085*B(318)+0.245*B(320)+0.322*B(322)+0.08*B(350)-B(435)

  JVS(701) = 0.1*B(121)+0.1*B(123)+0.05*B(125)+0.5*B(141)+0.07*B(146)+0.35*B(149)+B(151)+B(153)

  JVS(702) = 0.07*B(147)

  JVS(703) = 0.071*B(358)+0.237*B(360)+0.318*B(362)+0.08*B(385)-B(437)

  JVS(704) = 0.048*B(239)+0.334*B(241)+0.416*B(243)+0.16*B(267)-B(440)

  JVS(705) = 0

  JVS(706) = 0.5*B(167)+B(185)+0.051*B(399)+0.391*B(401)+0.587*B(403)+0.16*B(427)-B(442)

  JVS(707) = 0.25*B(254)+0.4*B(333)+0.4*B(373)+0.4*B(414)

  JVS(708) = 0.08*B(329)

  JVS(709) = -B(110)-B(174)-B(190)-B(505)

  JVS(710) = -B(191)

  JVS(711) = -B(506)

  JVS(712) = 0

  JVS(713) = 0.08*B(330)+0.4*B(334)

  JVS(714) = -B(111)

  JVS(715) = 0.4*B(374)

  JVS(716) = 0.25*B(255)

  JVS(717) = -B(175)+0.4*B(415)

  JVS(718) = 0.46*B(198)

  JVS(719) = 0.35*B(192)

  JVS(720) = 0.37*B(329)+0.229*B(369)

  JVS(721) = 0.3*B(505)

  JVS(722) = -B(108)-B(172)-B(188)-B(503)

  JVS(723) = -B(189)+0.35*B(193)+0.46*B(199)

  JVS(724) = -B(504)+0.3*B(506)

  JVS(725) = 0

  JVS(726) = 0.37*B(330)

  JVS(727) = -B(109)

  JVS(728) = 0.229*B(370)

  JVS(729) = 0

  JVS(730) = -B(173)

  JVS(731) = 0.9*B(460)

  JVS(732) = 0.1*B(180)

  JVS(733) = 0.13*B(182)

  JVS(734) = 0.8*B(170)

  JVS(735) = 0.11*B(464)

  JVS(736) = 0.9*B(176)

  JVS(737) = 0.11*B(174)

  JVS(738) = 0.43*B(172)

  JVS(739) = -B(268)-B(310)-B(351)-B(386)-2*B(390)-B(391)-B(428)

  JVS(740) = -B(392)

  JVS(741) = 0

  JVS(742) = 0

  JVS(743) = -B(311)

  JVS(744) = -B(352)

  JVS(745) = 0

  JVS(746) = -B(387)

  JVS(747) = -B(269)

  JVS(748) = 0.8*B(171)+0.43*B(173)+0.11*B(175)+0.9*B(177)+0.1*B(181)+0.13*B(183)-B(429)+0.9*B(461)+0.11*B(465)

  JVS(749) = 0.1*B(460)

  JVS(750) = 0.9*B(180)

  JVS(751) = 0.87*B(182)

  JVS(752) = 0.2*B(170)

  JVS(753) = 0.89*B(464)

  JVS(754) = 0.1*B(176)

  JVS(755) = 0.89*B(174)

  JVS(756) = 0.57*B(172)

  JVS(757) = -B(391)

  JVS(758) = -B(270)-B(312)-B(353)-B(388)-B(392)-2*B(393)-B(430)

  JVS(759) = 0

  JVS(760) = 0

  JVS(761) = -B(313)

  JVS(762) = -B(354)

  JVS(763) = 0

  JVS(764) = -B(389)

  JVS(765) = -B(271)

  JVS(766) = 0.2*B(171)+0.57*B(173)+0.89*B(175)+0.1*B(177)+0.9*B(181)+0.87*B(183)-B(431)+0.1*B(461)+0.89*B(465)

  JVS(767) = B(28)-B(30)

  JVS(768) = -B(214)

  JVS(769) = -B(220)

  JVS(770) = -B(226)

  JVS(771) = -B(458)

  JVS(772) = -B(196)

  JVS(773) = -B(198)

  JVS(774) = -B(186)

  JVS(775) = -B(466)

  JVS(776) = -B(192)

  JVS(777) = 0

  JVS(778) = -B(204)

  JVS(779) = -B(194)

  JVS(780) = -B(200)

  JVS(781) = 0

  JVS(782) = B(306)

  JVS(783) = -B(202)

  JVS(784) = -B(190)

  JVS(785) = -B(188)

  JVS(786) = -B(2)-B(3)-B(31)-B(36)-B(38)-B(72)-B(74)-B(187)-B(189)-B(191)-B(193)-B(195)-B(197)-B(199)-B(201)-B(203)&
               &-B(205)-B(215)-B(221)-B(227)-B(459)-B(467)-B(482)

  JVS(787) = -B(483)

  JVS(788) = -B(39)+B(302)+B(307)

  JVS(789) = 0

  JVS(790) = -B(37)

  JVS(791) = B(303)

  JVS(792) = -B(73)

  JVS(793) = -B(75)

  JVS(794) = 0

  JVS(795) = 0.583*B(100)

  JVS(796) = B(154)

  JVS(797) = 0

  JVS(798) = 0

  JVS(799) = 0

  JVS(800) = 0

  JVS(801) = 0

  JVS(802) = 0

  JVS(803) = 0

  JVS(804) = 0

  JVS(805) = -B(238)-B(276)-B(317)-B(357)-B(398)

  JVS(806) = 0

  JVS(807) = -B(277)

  JVS(808) = -B(318)

  JVS(809) = 0.583*B(101)+0.44*B(146)+B(155)

  JVS(810) = 0.44*B(147)

  JVS(811) = -B(358)

  JVS(812) = -B(239)

  JVS(813) = 0

  JVS(814) = -B(399)

  JVS(815) = 0.42*B(196)

  JVS(816) = 0.42*B(198)

  JVS(817) = 0.12*B(466)

  JVS(818) = 0.02*B(192)

  JVS(819) = B(130)

  JVS(820) = 0

  JVS(821) = 0

  JVS(822) = 0

  JVS(823) = 0

  JVS(824) = 0

  JVS(825) = 0.12*B(190)

  JVS(826) = 0.03*B(188)

  JVS(827) = 0

  JVS(828) = 0

  JVS(829) = 0.03*B(189)+0.12*B(191)+0.02*B(193)+0.42*B(197)+0.42*B(199)+0.12*B(467)

  JVS(830) = 0

  JVS(831) = -B(266)-B(308)-B(349)-B(384)-B(426)

  JVS(832) = 0

  JVS(833) = -B(309)

  JVS(834) = -B(350)

  JVS(835) = B(131)

  JVS(836) = 0

  JVS(837) = -B(385)

  JVS(838) = -B(267)

  JVS(839) = 0

  JVS(840) = -B(427)

  JVS(841) = 2*B(24)

  JVS(842) = B(26)

  JVS(843) = B(25)

  JVS(844) = -B(515)

  JVS(845) = -B(517)

  JVS(846) = 1.4*B(484)+B(485)

  JVS(847) = B(480)

  JVS(848) = -B(491)

  JVS(849) = B(27)+B(509)

  JVS(850) = -B(497)

  JVS(851) = -B(499)

  JVS(852) = -B(495)

  JVS(853) = -B(493)

  JVS(854) = -B(501)

  JVS(855) = 0

  JVS(856) = -B(507)

  JVS(857) = -B(511)

  JVS(858) = 0

  JVS(859) = 0

  JVS(860) = 0

  JVS(861) = -B(513)

  JVS(862) = 0

  JVS(863) = 0

  JVS(864) = -B(505)

  JVS(865) = -B(503)

  JVS(866) = 0

  JVS(867) = 0

  JVS(868) = -B(482)

  JVS(869) = 0

  JVS(870) = 0

  JVS(871) = -B(483)-B(489)-B(492)-B(494)-B(496)-B(498)-B(500)-B(502)-B(504)-B(506)-B(508)-B(512)-B(514)-B(516)-B(518)

  JVS(872) = 0

  JVS(873) = 0

  JVS(874) = B(481)+B(510)

  JVS(875) = 0

  JVS(876) = 0

  JVS(877) = B(486)

  JVS(878) = -B(490)

  JVS(879) = 0

  JVS(880) = 0

  JVS(881) = B(88)

  JVS(882) = B(142)

  JVS(883) = 0.1*B(120)+0.88*B(515)

  JVS(884) = 0.1*B(122)+0.84*B(517)

  JVS(885) = B(13)

  JVS(886) = 0.65*B(6)+B(63)

  JVS(887) = B(94)

  JVS(888) = -B(208)

  JVS(889) = -B(487)

  JVS(890) = B(42)

  JVS(891) = B(27)

  JVS(892) = B(23)+B(132)

  JVS(893) = 0.02*B(212)

  JVS(894) = 0.02*B(218)

  JVS(895) = 0.02*B(224)

  JVS(896) = 0.25*B(102)+0.11*B(497)

  JVS(897) = 0

  JVS(898) = 0.35*B(148)

  JVS(899) = 0.049*B(104)+0.11*B(499)

  JVS(900) = 0.1*B(196)

  JVS(901) = 0.381*B(100)+0.11*B(495)

  JVS(902) = B(493)

  JVS(903) = 0.1*B(198)

  JVS(904) = B(90)

  JVS(905) = 0.26*B(186)+B(501)

  JVS(906) = 0.22*B(466)

  JVS(907) = 0.3*B(192)

  JVS(908) = 0.75*B(448)-B(450)+B(452)

  JVS(909) = 0.05*B(124)

  JVS(910) = 0.4*B(152)+0.08*B(204)

  JVS(911) = 0.25*B(194)+B(507)

  JVS(912) = 0.8*B(18)+B(134)+B(160)

  JVS(913) = B(246)-B(284)+B(325)+0.501*B(365)+B(406)

  JVS(914) = B(244)-B(282)+B(323)+0.5*B(363)+B(404)

  JVS(915) = 0.65*B(254)-B(292)+2*B(333)+B(373)+B(414)

  JVS(916) = B(22)+0.32*B(200)

  JVS(917) = 0.8*B(252)-B(290)+2*B(331)+B(371)+B(412)

  JVS(918) = B(250)-B(288)+B(329)+0.506*B(369)+B(443)

  JVS(919) = B(260)-B(298)+2*B(339)+B(379)+B(420)

  JVS(920) = B(19)

  JVS(921) = B(248)-B(286)+B(327)+0.51*B(367)+B(408)

  JVS(922) = 0.95*B(256)-B(294)+B(335)+B(375)+B(416)

  JVS(923) = 0.95*B(258)-B(296)+B(337)+B(377)+B(418)

  JVS(924) = 2*B(11)+B(126)+B(156)+B(511)

  JVS(925) = -B(304)-B(306)+B(345)

  JVS(926) = B(21)

  JVS(927) = 0.599*B(240)-B(278)+0.946*B(319)+0.438*B(359)+0.699*B(400)

  JVS(928) = 0.606*B(242)-B(280)+0.993*B(321)+0.489*B(361)+0.845*B(402)

  JVS(929) = B(12)

  JVS(930) = B(20)+0.5*B(140)+0.5*B(166)+0.29*B(202)

  JVS(931) = -B(432)+B(434)

  JVS(932) = 0.22*B(190)+B(505)

  JVS(933) = 0.25*B(188)+B(503)

  JVS(934) = B(268)-B(310)+B(351)+0.5*B(386)+B(390)+0.5*B(391)+B(428)

  JVS(935) = -B(312)+0.5*B(353)+0.5*B(392)

  JVS(936) = B(36)-B(38)+0.26*B(187)+0.25*B(189)+0.22*B(191)+0.3*B(193)+0.25*B(195)+0.1*B(197)+0.1*B(199)+0.32*B(201)&
               &+0.29*B(203)+0.08*B(205)+0.22*B(467)

  JVS(937) = 0.742*B(238)-B(276)+0.992*B(317)+0.488*B(357)+0.792*B(398)

  JVS(938) = 0.77*B(266)-B(308)+0.88*B(349)+0.38*B(384)+0.77*B(426)

  JVS(939) = B(494)+0.11*B(496)+0.11*B(498)+0.11*B(500)+B(502)+B(504)+B(506)+B(508)+B(512)+0.88*B(516)+0.84*B(518)

  JVS(940) = -B(39)-B(40)-2*B(44)-2*B(45)-B(59)-B(61)-B(64)-B(209)-B(272)-B(274)-B(277)-B(279)-B(281)-B(283)-B(285)&
               &-B(287)-B(289)-B(291)-B(293)-B(295)-B(297)-B(299)-B(300)-B(302)-B(305)-B(307)-B(309)-B(311)-B(313)-B(433)&
               &-B(451)-B(488)

  JVS(941) = B(234)-B(273)+0.66*B(314)+B(315)+0.992*B(318)+0.946*B(320)+0.993*B(322)+B(324)+B(326)+B(328)+B(330)+2&
               &*B(332)+2*B(334)+B(336)+B(338)+2*B(340)+B(341)+B(346)+0.88*B(350)+B(352)+0.5*B(354)+B(394)+B(435)

  JVS(942) = B(37)-B(41)+B(43)+B(57)+B(86)+B(89)+B(91)+B(95)+0.381*B(101)+0.25*B(103)+0.049*B(105)+0.1*B(121)+0.1*B(123)&
               &+0.05*B(125)+B(127)+B(133)+B(135)+0.5*B(141)+B(143)+0.35*B(149)+0.4*B(153)

  JVS(943) = B(14)

  JVS(944) = -B(301)-B(303)+B(342)+0.5*B(355)+0.488*B(358)+0.438*B(360)+0.489*B(362)+0.5*B(364)+0.501*B(366)+0.51*B(368)&
               &+0.506*B(370)+B(372)+B(374)+B(376)+B(378)+B(380)+0.38*B(385)+0.5*B(387)

  JVS(945) = -B(60)+B(235)+B(236)+0.742*B(239)+0.599*B(241)+0.606*B(243)+B(245)+B(247)+B(249)+B(251)+0.8*B(253)+0.65&
               &*B(255)+0.95*B(257)+0.95*B(259)+B(261)+0.77*B(267)+B(269)+0.75*B(449)

  JVS(946) = -B(62)

  JVS(947) = B(58)-B(65)+B(157)+B(161)+0.5*B(167)+B(395)+B(396)+0.792*B(399)+0.699*B(401)+0.845*B(403)+B(405)+B(407)&
               &+B(409)+B(413)+B(415)+B(417)+B(419)+B(421)+0.77*B(427)+B(429)

  JVS(948) = B(237)-B(275)+B(316)+0.5*B(356)+B(397)

  JVS(949) = 0.65*B(144)

  JVS(950) = B(96)+B(491)

  JVS(951) = B(15)

  JVS(952) = 0.23*B(466)

  JVS(953) = 0.03*B(192)

  JVS(954) = 0.08*B(194)

  JVS(955) = -B(325)+0.501*B(365)

  JVS(956) = -B(323)+0.5*B(363)

  JVS(957) = -B(333)+B(373)

  JVS(958) = -B(331)+B(371)

  JVS(959) = -B(329)+0.506*B(369)

  JVS(960) = -B(339)+B(379)

  JVS(961) = -B(327)+0.51*B(367)

  JVS(962) = -B(335)+B(375)

  JVS(963) = -B(337)+B(377)

  JVS(964) = -B(345)-B(347)+B(382)

  JVS(965) = 0.031*B(240)-0.951*B(319)+0.554*B(359)+0.04*B(400)

  JVS(966) = -B(321)+0.507*B(361)

  JVS(967) = B(12)

  JVS(968) = 0

  JVS(969) = -B(434)+B(436)

  JVS(970) = 0.23*B(190)

  JVS(971) = 0.19*B(188)

  JVS(972) = -B(351)+0.5*B(386)

  JVS(973) = -B(353)+0.516*B(388)

  JVS(974) = 0.19*B(189)+0.23*B(191)+0.03*B(193)+0.08*B(195)+0.23*B(467)

  JVS(975) = 0.15*B(238)-0.993*B(317)+0.508*B(357)+0.155*B(398)

  JVS(976) = -B(349)+0.5*B(384)

  JVS(977) = B(492)

  JVS(978) = -B(272)

  JVS(979) = -B(234)-B(273)-2*B(314)-B(315)-0.993*B(318)-0.951*B(320)-B(322)-B(324)-B(326)-B(328)-B(330)-B(332)-B(334)&
               &-B(336)-B(338)-B(340)-B(343)-B(346)-B(348)-B(350)-B(352)-B(354)-B(394)-B(435)

  JVS(980) = B(97)+0.65*B(145)

  JVS(981) = 0

  JVS(982) = B(262)-B(344)+0.5*B(355)+0.508*B(358)+0.554*B(360)+0.507*B(362)+0.5*B(364)+0.501*B(366)+0.51*B(368)+0.506&
               &*B(370)+B(372)+B(374)+B(376)+B(378)+B(380)+2*B(381)+B(383)+0.5*B(385)+0.5*B(387)+0.516*B(389)+B(422)+B(437)

  JVS(983) = -B(235)+0.15*B(239)+0.031*B(241)+B(263)

  JVS(984) = 0

  JVS(985) = -B(395)+0.155*B(399)+0.04*B(401)+B(423)

  JVS(986) = -B(316)+0.5*B(356)

  JVS(987) = 2*B(34)

  JVS(988) = -B(88)

  JVS(989) = B(25)

  JVS(990) = -B(164)

  JVS(991) = -B(142)

  JVS(992) = -B(92)

  JVS(993) = -B(120)

  JVS(994) = -B(122)

  JVS(995) = B(13)-0.65*B(144)

  JVS(996) = -B(410)

  JVS(997) = 0.35*B(6)-B(70)

  JVS(998) = -B(94)

  JVS(999) = B(4)-B(66)

  JVS(1000) = 0

  JVS(1001) = 2*B(9)-B(42)

  JVS(1002) = -B(480)

  JVS(1003) = -B(96)

  JVS(1004) = -B(509)

  JVS(1005) = -B(132)

  JVS(1006) = B(214)

  JVS(1007) = B(220)

  JVS(1008) = B(226)

  JVS(1009) = -B(102)

  JVS(1010) = -B(446)

  JVS(1011) = -B(456)+0.63*B(458)

  JVS(1012) = B(5)-B(68)

  JVS(1013) = B(15)-B(148)

  JVS(1014) = -B(104)

  JVS(1015) = -B(116)+0.85*B(196)

  JVS(1016) = -B(150)

  JVS(1017) = -0.964*B(100)

  JVS(1018) = -B(98)

  JVS(1019) = -B(118)+0.85*B(198)

  JVS(1020) = -B(90)

  JVS(1021) = -B(106)+0.12*B(186)

  JVS(1022) = -B(462)+0.63*B(466)

  JVS(1023) = -B(112)+0.28*B(192)

  JVS(1024) = 0

  JVS(1025) = -B(124)

  JVS(1026) = -B(152)+0.036*B(204)

  JVS(1027) = -B(114)+0.25*B(194)

  JVS(1028) = -B(134)

  JVS(1029) = -B(138)+0.19*B(200)

  JVS(1030) = -B(130)

  JVS(1031) = 0

  JVS(1032) = 0

  JVS(1033) = -B(136)

  JVS(1034) = 0

  JVS(1035) = 0

  JVS(1036) = 0

  JVS(1037) = -B(126)

  JVS(1038) = 0

  JVS(1039) = -B(154)

  JVS(1040) = 0

  JVS(1041) = 0

  JVS(1042) = -B(128)

  JVS(1043) = -B(140)+0.21*B(202)

  JVS(1044) = 0

  JVS(1045) = -B(110)+0.63*B(190)

  JVS(1046) = -B(108)+0.4*B(188)

  JVS(1047) = 0

  JVS(1048) = 0

  JVS(1049) = -B(36)+B(38)+0.12*B(187)+0.4*B(189)+0.63*B(191)+0.28*B(193)+0.25*B(195)+0.85*B(197)+0.85*B(199)+0.19&
                &*B(201)+0.21*B(203)+0.036*B(205)+B(215)+B(221)+B(227)+0.63*B(459)+0.63*B(467)

  JVS(1050) = 0

  JVS(1051) = 0

  JVS(1052) = 0

  JVS(1053) = B(39)-B(40)+B(59)+0.7*B(64)

  JVS(1054) = 0

  JVS(1055) = -B(37)-B(41)-B(43)-B(53)-B(55)-B(57)-B(67)-B(69)-B(71)-B(86)-B(89)-B(91)-B(93)-B(95)-B(97)-B(99)-0.964&
                &*B(101)-B(103)-B(105)-B(107)-B(109)-B(111)-B(113)-B(115)-B(117)-B(119)-B(121)-B(123)-B(125)-B(127)-B(129)&
                &-B(131)-B(133)-B(135)-B(137)-B(139)-B(141)-B(143)-0.65*B(145)-0.51*B(146)-B(149)-B(151)-B(153)-B(155)&
                &-B(165)-B(411)-B(447)-B(457)-B(463)-B(481)-B(510)

  JVS(1056) = B(14)-0.51*B(147)

  JVS(1057) = 0

  JVS(1058) = -B(54)+B(60)

  JVS(1059) = -B(56)

  JVS(1060) = -B(58)+0.7*B(65)

  JVS(1061) = 0

  JVS(1062) = B(284)

  JVS(1063) = B(282)

  JVS(1064) = B(292)

  JVS(1065) = B(290)

  JVS(1066) = B(298)

  JVS(1067) = B(286)

  JVS(1068) = B(294)

  JVS(1069) = B(296)

  JVS(1070) = B(304)

  JVS(1071) = B(278)

  JVS(1072) = B(280)

  JVS(1073) = 0

  JVS(1074) = B(432)

  JVS(1075) = 0

  JVS(1076) = 0

  JVS(1077) = 0

  JVS(1078) = B(276)

  JVS(1079) = B(308)

  JVS(1080) = 0

  JVS(1081) = B(274)+B(277)+B(279)+B(281)+B(283)+B(285)+B(287)+B(291)+B(293)+B(295)+B(297)+B(299)+B(305)+B(309)+B(433)

  JVS(1082) = 0

  JVS(1083) = -B(146)

  JVS(1084) = -B(14)-B(147)

  JVS(1085) = 0

  JVS(1086) = 0

  JVS(1087) = 0

  JVS(1088) = 0

  JVS(1089) = B(275)

  JVS(1090) = B(23)

  JVS(1091) = 0.65*B(148)

  JVS(1092) = 0

  JVS(1093) = B(230)

  JVS(1094) = 0.15*B(192)

  JVS(1095) = 0.25*B(448)

  JVS(1096) = 0.7*B(204)

  JVS(1097) = 0.1*B(194)

  JVS(1098) = -B(365)

  JVS(1099) = -B(363)

  JVS(1100) = -B(373)

  JVS(1101) = B(22)+0.1*B(200)

  JVS(1102) = -B(371)

  JVS(1103) = B(16)

  JVS(1104) = -B(369)

  JVS(1105) = -B(379)

  JVS(1106) = B(19)+B(136)+B(162)

  JVS(1107) = -B(367)

  JVS(1108) = -B(375)

  JVS(1109) = -B(377)

  JVS(1110) = B(264)+B(345)+B(424)

  JVS(1111) = 0

  JVS(1112) = -B(359)

  JVS(1113) = -B(361)

  JVS(1114) = B(128)+B(158)+B(513)

  JVS(1115) = 0.28*B(202)

  JVS(1116) = -B(436)

  JVS(1117) = 0

  JVS(1118) = 0

  JVS(1119) = -B(386)

  JVS(1120) = -B(388)

  JVS(1121) = 0.15*B(193)+0.1*B(195)+0.1*B(201)+0.28*B(203)+0.7*B(205)

  JVS(1122) = -B(357)

  JVS(1123) = 0.23*B(266)+0.12*B(349)-0.88*B(384)+0.23*B(426)

  JVS(1124) = B(514)

  JVS(1125) = -B(300)-B(302)

  JVS(1126) = -B(341)-B(343)+B(346)+0.12*B(350)

  JVS(1127) = B(129)+B(137)+0.65*B(149)

  JVS(1128) = 0

  JVS(1129) = -B(228)-B(262)-B(301)-B(303)-B(342)-B(344)-B(355)-B(358)-B(360)-B(362)-B(364)-B(366)-B(368)-B(370)-B(372)&
                &-B(374)-B(376)-B(378)-B(380)-2*B(381)-0.88*B(385)-B(387)-B(389)-B(422)-B(437)

  JVS(1130) = -B(263)+B(265)+0.23*B(267)+0.25*B(449)

  JVS(1131) = -B(229)

  JVS(1132) = B(159)+B(163)-B(423)+B(425)+0.23*B(427)

  JVS(1133) = -B(356)

  JVS(1134) = B(4)

  JVS(1135) = -B(47)+B(49)

  JVS(1136) = -B(485)

  JVS(1137) = 0

  JVS(1138) = 0

  JVS(1139) = 0

  JVS(1140) = 0

  JVS(1141) = -B(448)

  JVS(1142) = 0

  JVS(1143) = 0

  JVS(1144) = -B(246)

  JVS(1145) = -B(244)

  JVS(1146) = -B(254)

  JVS(1147) = 0

  JVS(1148) = -B(252)

  JVS(1149) = -B(250)

  JVS(1150) = -B(260)

  JVS(1151) = -B(248)

  JVS(1152) = -B(256)

  JVS(1153) = -B(258)

  JVS(1154) = -B(264)

  JVS(1155) = -B(240)

  JVS(1156) = -B(242)

  JVS(1157) = 0

  JVS(1158) = -B(439)

  JVS(1159) = 0

  JVS(1160) = 0

  JVS(1161) = -B(268)

  JVS(1162) = -B(270)

  JVS(1163) = -B(72)

  JVS(1164) = -B(238)

  JVS(1165) = -B(266)

  JVS(1166) = 0

  JVS(1167) = -B(59)

  JVS(1168) = -B(234)

  JVS(1169) = -B(53)

  JVS(1170) = 0

  JVS(1171) = -B(262)

  JVS(1172) = -B(48)-B(54)-B(60)-B(73)-2*B(76)-B(78)-B(235)-B(236)-B(239)-B(241)-B(243)-B(245)-B(247)-B(249)-B(251)&
                &-B(253)-B(255)-B(257)-B(259)-B(261)-B(263)-B(265)-B(267)-B(269)-B(271)-B(440)-B(449)-B(486)

  JVS(1173) = B(1)+B(50)+B(80)

  JVS(1174) = B(7)-B(79)+B(81)

  JVS(1175) = -B(237)

  JVS(1176) = B(84)

  JVS(1177) = B(26)

  JVS(1178) = B(92)

  JVS(1179) = B(410)+B(455)

  JVS(1180) = 0.65*B(6)+B(63)+B(70)

  JVS(1181) = B(66)

  JVS(1182) = B(47)-B(49)-B(51)

  JVS(1183) = -B(206)

  JVS(1184) = B(485)

  JVS(1185) = -B(210)

  JVS(1186) = -B(216)

  JVS(1187) = -B(222)

  JVS(1188) = 0

  JVS(1189) = B(5)

  JVS(1190) = B(230)

  JVS(1191) = 0

  JVS(1192) = B(448)-B(453)

  JVS(1193) = 0

  JVS(1194) = 0.4*B(184)+0.7*B(204)+B(233)

  JVS(1195) = 0

  JVS(1196) = 0

  JVS(1197) = B(246)+B(406)

  JVS(1198) = B(244)+B(404)

  JVS(1199) = 0.65*B(254)+B(414)

  JVS(1200) = 0

  JVS(1201) = 0.8*B(252)+B(412)

  JVS(1202) = B(250)

  JVS(1203) = B(260)+B(420)

  JVS(1204) = 0

  JVS(1205) = B(248)+B(408)

  JVS(1206) = 0.95*B(256)+B(416)

  JVS(1207) = 0.95*B(258)+B(418)

  JVS(1208) = 0

  JVS(1209) = -B(231)+B(264)+B(424)

  JVS(1210) = B(21)+B(154)

  JVS(1211) = 0.876*B(240)+B(400)

  JVS(1212) = 0.739*B(242)+B(402)

  JVS(1213) = 0

  JVS(1214) = 0.5*B(166)

  JVS(1215) = B(439)+B(441)

  JVS(1216) = 0

  JVS(1217) = 0

  JVS(1218) = B(268)+0.5*B(391)+B(428)

  JVS(1219) = 2*B(270)+0.5*B(353)+0.516*B(388)+0.5*B(392)+B(393)+2*B(430)

  JVS(1220) = B(72)-B(74)+0.7*B(205)

  JVS(1221) = 0.941*B(238)+B(398)

  JVS(1222) = B(266)+B(426)

  JVS(1223) = -B(489)

  JVS(1224) = B(59)-B(61)+0.7*B(64)

  JVS(1225) = B(234)+0.5*B(354)+B(394)

  JVS(1226) = -B(55)+B(57)+B(67)+B(71)+B(93)+B(155)+B(411)

  JVS(1227) = 0

  JVS(1228) = -B(228)+B(262)+0.516*B(389)+B(422)

  JVS(1229) = B(48)+B(60)+B(73)+2*B(76)+2*B(78)+B(235)+B(236)+0.941*B(239)+0.876*B(241)+0.739*B(243)+B(245)+B(247)&
                &+B(249)+B(251)+0.8*B(253)+0.65*B(255)+0.95*B(257)+0.95*B(259)+B(261)+B(263)+B(265)+B(267)+B(269)+2*B(271)&
                &+B(440)+B(449)+B(486)

  JVS(1230) = -B(1)-B(50)-B(52)-B(56)-B(62)-B(75)-B(82)-B(207)-B(211)-B(217)-B(223)-B(229)-B(232)-B(454)-B(490)

  JVS(1231) = B(8)+B(58)+0.7*B(65)+2*B(79)-B(83)+2*B(85)+0.5*B(167)+0.4*B(185)+B(395)+B(396)+B(399)+B(401)+B(403)+B(405)&
                &+B(407)+B(409)+B(413)+B(415)+B(417)+B(419)+B(421)+B(423)+B(425)+B(427)+B(429)+2*B(431)+B(442)

  JVS(1232) = B(237)+B(397)

  JVS(1233) = B(84)

  JVS(1234) = 0.35*B(6)

  JVS(1235) = B(51)

  JVS(1236) = -B(460)

  JVS(1237) = B(68)

  JVS(1238) = -B(180)

  JVS(1239) = B(150)

  JVS(1240) = -B(182)

  JVS(1241) = -B(170)

  JVS(1242) = -B(464)

  JVS(1243) = -B(176)

  JVS(1244) = -B(168)

  JVS(1245) = 0.6*B(152)-0.4*B(184)

  JVS(1246) = -B(178)

  JVS(1247) = -B(160)

  JVS(1248) = -B(406)

  JVS(1249) = -B(404)

  JVS(1250) = -B(414)

  JVS(1251) = -B(412)

  JVS(1252) = -B(420)

  JVS(1253) = -B(162)

  JVS(1254) = -B(408)

  JVS(1255) = -B(416)

  JVS(1256) = -B(418)

  JVS(1257) = -B(156)

  JVS(1258) = -B(424)

  JVS(1259) = -B(400)

  JVS(1260) = -B(402)

  JVS(1261) = -B(158)

  JVS(1262) = -B(166)

  JVS(1263) = -B(441)

  JVS(1264) = -B(174)

  JVS(1265) = -B(172)

  JVS(1266) = -B(428)

  JVS(1267) = -B(430)

  JVS(1268) = B(74)

  JVS(1269) = -B(398)

  JVS(1270) = -B(426)

  JVS(1271) = 0

  JVS(1272) = -B(64)

  JVS(1273) = -B(394)

  JVS(1274) = -B(57)+B(69)+B(151)+0.6*B(153)

  JVS(1275) = 0

  JVS(1276) = -B(422)

  JVS(1277) = -B(78)

  JVS(1278) = B(52)+B(75)-B(80)-B(82)

  JVS(1279) = -B(7)-B(8)-B(58)-B(65)-B(79)-B(81)-B(83)-2*B(85)-B(157)-B(159)-B(161)-B(163)-B(167)-B(169)-B(171)-B(173)&
                &-B(175)-B(177)-B(179)-B(181)-B(183)-0.4*B(185)-B(395)-B(396)-B(399)-B(401)-B(403)-B(405)-B(407)-B(409)&
                &-B(413)-B(415)-B(417)-B(419)-B(421)-B(423)-B(425)-B(427)-B(429)-B(431)-B(442)-B(461)-B(465)

  JVS(1280) = -B(397)

  JVS(1281) = 0.2*B(196)

  JVS(1282) = B(98)

  JVS(1283) = 0.16*B(198)

  JVS(1284) = 0.18*B(466)

  JVS(1285) = B(16)

  JVS(1286) = 0

  JVS(1287) = 0

  JVS(1288) = 0.245*B(240)+0.014*B(319)+0.018*B(359)+0.262*B(400)

  JVS(1289) = 0.133*B(242)+0.013*B(321)+0.015*B(361)+0.155*B(402)

  JVS(1290) = 0

  JVS(1291) = 0.18*B(190)

  JVS(1292) = 0.1*B(188)

  JVS(1293) = 0

  JVS(1294) = 0

  JVS(1295) = 0.1*B(189)+0.18*B(191)+0.2*B(197)+0.16*B(199)+0.18*B(467)

  JVS(1296) = 0.048*B(238)+0.006*B(357)+0.053*B(398)

  JVS(1297) = 0

  JVS(1298) = 0

  JVS(1299) = -B(274)

  JVS(1300) = -B(315)+0.014*B(320)+0.013*B(322)

  JVS(1301) = B(99)

  JVS(1302) = 0

  JVS(1303) = -B(355)+0.006*B(358)+0.018*B(360)+0.015*B(362)

  JVS(1304) = -B(236)+0.048*B(239)+0.245*B(241)+0.133*B(243)

  JVS(1305) = 0

  JVS(1306) = -B(396)+0.053*B(399)+0.262*B(401)+0.155*B(403)

  JVS(1307) = -B(237)-B(275)-B(316)-B(356)-B(397)
      
END SUBROUTINE racm_soa_vbs_het_Jac_SP














SUBROUTINE racm_soa_vbs_het_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(1307), W(97), a
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
      
END SUBROUTINE racm_soa_vbs_het_KppDecomp



SUBROUTINE racm_soa_vbs_het_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(1307), W(97), a
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
      
END SUBROUTINE racm_soa_vbs_het_KppDecompCmplx


SUBROUTINE racm_soa_vbs_het_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(1307), X(97), sum

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
      
END SUBROUTINE racm_soa_vbs_het_KppSolveIndirect


SUBROUTINE racm_soa_vbs_het_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(1307), X(97), sum

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
      
END SUBROUTINE racm_soa_vbs_het_KppSolveCmplx













SUBROUTINE racm_soa_vbs_het_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(28) = X(28)-JVS(112)*X(26)
  X(30) = X(30)-JVS(126)*X(14)
  X(34) = X(34)-JVS(155)*X(23)-JVS(156)*X(24)
  X(38) = X(38)-JVS(189)*X(23)
  X(39) = X(39)-JVS(195)*X(24)
  X(52) = X(52)-JVS(270)*X(22)-JVS(271)*X(36)-JVS(272)*X(42)-JVS(273)*X(47)-JVS(274)*X(49)-JVS(275)*X(51)
  X(56) = X(56)-JVS(310)*X(20)-JVS(311)*X(26)
  X(57) = X(57)-JVS(320)*X(31)-JVS(321)*X(38)-JVS(322)*X(39)-JVS(323)*X(40)
  X(60) = X(60)-JVS(342)*X(49)
  X(62) = X(62)-JVS(366)*X(53)
  X(63) = X(63)-JVS(376)*X(51)
  X(64) = X(64)-JVS(385)*X(16)-JVS(386)*X(43)-JVS(387)*X(51)-JVS(388)*X(55)-JVS(389)*X(59)-JVS(390)*X(63)
  X(65) = X(65)-JVS(401)*X(47)
  X(66) = X(66)-JVS(410)*X(21)-JVS(411)*X(41)-JVS(412)*X(43)-JVS(413)*X(47)-JVS(414)*X(54)-JVS(415)*X(61)-JVS(416)*X(65)
  X(67) = X(67)-JVS(438)*X(55)-JVS(439)*X(59)
  X(68) = X(68)-JVS(449)*X(40)-JVS(450)*X(57)
  X(69) = X(69)-JVS(461)*X(28)-JVS(462)*X(37)-JVS(463)*X(42)-JVS(464)*X(46)-JVS(465)*X(56)-JVS(466)*X(58)-JVS(467)*X(59)&
            &-JVS(468)*X(64)-JVS(469)*X(67)-JVS(470)*X(68)
  X(70) = X(70)-JVS(489)*X(43)-JVS(490)*X(54)
  X(71) = X(71)-JVS(500)*X(38)
  X(72) = X(72)-JVS(511)*X(39)
  X(73) = X(73)-JVS(522)*X(22)-JVS(523)*X(25)-JVS(524)*X(37)-JVS(525)*X(42)-JVS(526)*X(43)-JVS(527)*X(45)-JVS(528)*X(46)&
            &-JVS(529)*X(48)-JVS(530)*X(49)-JVS(531)*X(51)-JVS(532)*X(53)-JVS(533)*X(54)-JVS(534)*X(55)-JVS(535)*X(56)&
            &-JVS(536)*X(58)-JVS(537)*X(59)-JVS(538)*X(60)-JVS(539)*X(61)-JVS(540)*X(62)-JVS(541)*X(63)-JVS(542)*X(64)&
            &-JVS(543)*X(65)-JVS(544)*X(67)-JVS(545)*X(68)-JVS(546)*X(70)-JVS(547)*X(71)-JVS(548)*X(72)
  X(74) = X(74)-JVS(571)*X(58)
  X(75) = X(75)-JVS(582)*X(31)-JVS(583)*X(57)-JVS(584)*X(58)-JVS(585)*X(63)-JVS(586)*X(65)-JVS(587)*X(71)-JVS(588)*X(72)&
            &-JVS(589)*X(74)
  X(76) = X(76)-JVS(606)*X(41)
  X(77) = X(77)-JVS(617)*X(46)
  X(78) = X(78)-JVS(628)*X(21)-JVS(629)*X(41)-JVS(630)*X(43)-JVS(631)*X(46)-JVS(632)*X(47)-JVS(633)*X(49)-JVS(634)*X(50)&
            &-JVS(635)*X(54)-JVS(636)*X(61)-JVS(637)*X(62)-JVS(638)*X(65)-JVS(639)*X(70)-JVS(640)*X(75)-JVS(641)*X(76)&
            &-JVS(642)*X(77)
  X(79) = X(79)-JVS(662)*X(71)-JVS(663)*X(72)
  X(80) = X(80)-JVS(674)*X(23)-JVS(675)*X(24)-JVS(676)*X(41)-JVS(677)*X(45)-JVS(678)*X(46)-JVS(679)*X(48)-JVS(680)*X(49)&
            &-JVS(681)*X(50)-JVS(682)*X(53)-JVS(683)*X(54)-JVS(684)*X(55)-JVS(685)*X(57)-JVS(686)*X(58)-JVS(687)*X(59)&
            &-JVS(688)*X(74)-JVS(689)*X(76)-JVS(690)*X(77)-JVS(691)*X(79)
  X(81) = X(81)-JVS(707)*X(63)-JVS(708)*X(67)
  X(82) = X(82)-JVS(718)*X(51)-JVS(719)*X(55)-JVS(720)*X(67)-JVS(721)*X(81)
  X(83) = X(83)-JVS(731)*X(43)-JVS(732)*X(47)-JVS(733)*X(51)-JVS(734)*X(53)-JVS(735)*X(54)-JVS(736)*X(55)-JVS(737)*X(81)&
            &-JVS(738)*X(82)
  X(84) = X(84)-JVS(749)*X(43)-JVS(750)*X(47)-JVS(751)*X(51)-JVS(752)*X(53)-JVS(753)*X(54)-JVS(754)*X(55)-JVS(755)*X(81)&
            &-JVS(756)*X(82)-JVS(757)*X(83)
  X(85) = X(85)-JVS(767)*X(30)-JVS(768)*X(38)-JVS(769)*X(39)-JVS(770)*X(40)-JVS(771)*X(43)-JVS(772)*X(47)-JVS(773)*X(51)&
            &-JVS(774)*X(53)-JVS(775)*X(54)-JVS(776)*X(55)-JVS(777)*X(57)-JVS(778)*X(58)-JVS(779)*X(59)-JVS(780)*X(64)&
            &-JVS(781)*X(67)-JVS(782)*X(74)-JVS(783)*X(79)-JVS(784)*X(81)-JVS(785)*X(82)
  X(86) = X(86)-JVS(795)*X(49)-JVS(796)*X(75)-JVS(797)*X(76)-JVS(798)*X(77)-JVS(799)*X(79)-JVS(800)*X(81)-JVS(801)*X(82)&
            &-JVS(802)*X(83)-JVS(803)*X(84)-JVS(804)*X(85)
  X(87) = X(87)-JVS(815)*X(47)-JVS(816)*X(51)-JVS(817)*X(54)-JVS(818)*X(55)-JVS(819)*X(66)-JVS(820)*X(70)-JVS(821)*X(75)&
            &-JVS(822)*X(76)-JVS(823)*X(77)-JVS(824)*X(79)-JVS(825)*X(81)-JVS(826)*X(82)-JVS(827)*X(83)-JVS(828)*X(84)&
            &-JVS(829)*X(85)-JVS(830)*X(86)
  X(88) = X(88)-JVS(841)*X(13)-JVS(842)*X(18)-JVS(843)*X(19)-JVS(844)*X(23)-JVS(845)*X(24)-JVS(846)*X(32)-JVS(847)*X(34)&
            &-JVS(848)*X(35)-JVS(849)*X(36)-JVS(850)*X(41)-JVS(851)*X(46)-JVS(852)*X(49)-JVS(853)*X(50)-JVS(854)*X(53)&
            &-JVS(855)*X(54)-JVS(856)*X(59)-JVS(857)*X(73)-JVS(858)*X(74)-JVS(859)*X(76)-JVS(860)*X(77)-JVS(861)*X(78)&
            &-JVS(862)*X(79)-JVS(863)*X(80)-JVS(864)*X(81)-JVS(865)*X(82)-JVS(866)*X(83)-JVS(867)*X(84)-JVS(868)*X(85)&
            &-JVS(869)*X(86)-JVS(870)*X(87)
  X(89) = X(89)-JVS(881)*X(15)-JVS(882)*X(21)-JVS(883)*X(23)-JVS(884)*X(24)-JVS(885)*X(25)-JVS(886)*X(27)-JVS(887)*X(28)&
            &-JVS(888)*X(31)-JVS(889)*X(32)-JVS(890)*X(33)-JVS(891)*X(36)-JVS(892)*X(37)-JVS(893)*X(38)-JVS(894)*X(39)&
            &-JVS(895)*X(40)-JVS(896)*X(41)-JVS(897)*X(42)-JVS(898)*X(45)-JVS(899)*X(46)-JVS(900)*X(47)-JVS(901)*X(49)&
            &-JVS(902)*X(50)-JVS(903)*X(51)-JVS(904)*X(52)-JVS(905)*X(53)-JVS(906)*X(54)-JVS(907)*X(55)-JVS(908)*X(56)&
            &-JVS(909)*X(57)-JVS(910)*X(58)-JVS(911)*X(59)-JVS(912)*X(60)-JVS(913)*X(61)-JVS(914)*X(62)-JVS(915)*X(63)&
            &-JVS(916)*X(64)-JVS(917)*X(65)-JVS(918)*X(67)-JVS(919)*X(68)-JVS(920)*X(69)-JVS(921)*X(70)-JVS(922)*X(71)&
            &-JVS(923)*X(72)-JVS(924)*X(73)-JVS(925)*X(74)-JVS(926)*X(75)-JVS(927)*X(76)-JVS(928)*X(77)-JVS(929)*X(78)&
            &-JVS(930)*X(79)-JVS(931)*X(80)-JVS(932)*X(81)-JVS(933)*X(82)-JVS(934)*X(83)-JVS(935)*X(84)-JVS(936)*X(85)&
            &-JVS(937)*X(86)-JVS(938)*X(87)-JVS(939)*X(88)
  X(90) = X(90)-JVS(949)*X(25)-JVS(950)*X(35)-JVS(951)*X(45)-JVS(952)*X(54)-JVS(953)*X(55)-JVS(954)*X(59)-JVS(955)*X(61)&
            &-JVS(956)*X(62)-JVS(957)*X(63)-JVS(958)*X(65)-JVS(959)*X(67)-JVS(960)*X(68)-JVS(961)*X(70)-JVS(962)*X(71)&
            &-JVS(963)*X(72)-JVS(964)*X(74)-JVS(965)*X(76)-JVS(966)*X(77)-JVS(967)*X(78)-JVS(968)*X(79)-JVS(969)*X(80)&
            &-JVS(970)*X(81)-JVS(971)*X(82)-JVS(972)*X(83)-JVS(973)*X(84)-JVS(974)*X(85)-JVS(975)*X(86)-JVS(976)*X(87)&
            &-JVS(977)*X(88)-JVS(978)*X(89)
  X(91) = X(91)-JVS(987)*X(14)-JVS(988)*X(15)-JVS(989)*X(19)-JVS(990)*X(20)-JVS(991)*X(21)-JVS(992)*X(22)-JVS(993)*X(23)&
            &-JVS(994)*X(24)-JVS(995)*X(25)-JVS(996)*X(26)-JVS(997)*X(27)-JVS(998)*X(28)-JVS(999)*X(29)-JVS(1000)*X(32)&
            &-JVS(1001)*X(33)-JVS(1002)*X(34)-JVS(1003)*X(35)-JVS(1004)*X(36)-JVS(1005)*X(37)-JVS(1006)*X(38)-JVS(1007)&
            &*X(39)-JVS(1008)*X(40)-JVS(1009)*X(41)-JVS(1010)*X(42)-JVS(1011)*X(43)-JVS(1012)*X(44)-JVS(1013)*X(45)&
            &-JVS(1014)*X(46)-JVS(1015)*X(47)-JVS(1016)*X(48)-JVS(1017)*X(49)-JVS(1018)*X(50)-JVS(1019)*X(51)-JVS(1020)&
            &*X(52)-JVS(1021)*X(53)-JVS(1022)*X(54)-JVS(1023)*X(55)-JVS(1024)*X(56)-JVS(1025)*X(57)-JVS(1026)*X(58)&
            &-JVS(1027)*X(59)-JVS(1028)*X(60)-JVS(1029)*X(64)-JVS(1030)*X(66)-JVS(1031)*X(67)-JVS(1032)*X(68)-JVS(1033)&
            &*X(69)-JVS(1034)*X(70)-JVS(1035)*X(71)-JVS(1036)*X(72)-JVS(1037)*X(73)-JVS(1038)*X(74)-JVS(1039)*X(75)&
            &-JVS(1040)*X(76)-JVS(1041)*X(77)-JVS(1042)*X(78)-JVS(1043)*X(79)-JVS(1044)*X(80)-JVS(1045)*X(81)-JVS(1046)&
            &*X(82)-JVS(1047)*X(83)-JVS(1048)*X(84)-JVS(1049)*X(85)-JVS(1050)*X(86)-JVS(1051)*X(87)-JVS(1052)*X(88)&
            &-JVS(1053)*X(89)-JVS(1054)*X(90)
  X(92) = X(92)-JVS(1062)*X(61)-JVS(1063)*X(62)-JVS(1064)*X(63)-JVS(1065)*X(65)-JVS(1066)*X(68)-JVS(1067)*X(70)&
            &-JVS(1068)*X(71)-JVS(1069)*X(72)-JVS(1070)*X(74)-JVS(1071)*X(76)-JVS(1072)*X(77)-JVS(1073)*X(79)-JVS(1074)&
            &*X(80)-JVS(1075)*X(81)-JVS(1076)*X(82)-JVS(1077)*X(85)-JVS(1078)*X(86)-JVS(1079)*X(87)-JVS(1080)*X(88)&
            &-JVS(1081)*X(89)-JVS(1082)*X(90)-JVS(1083)*X(91)
  X(93) = X(93)-JVS(1090)*X(37)-JVS(1091)*X(45)-JVS(1092)*X(46)-JVS(1093)*X(48)-JVS(1094)*X(55)-JVS(1095)*X(56)&
            &-JVS(1096)*X(58)-JVS(1097)*X(59)-JVS(1098)*X(61)-JVS(1099)*X(62)-JVS(1100)*X(63)-JVS(1101)*X(64)-JVS(1102)&
            &*X(65)-JVS(1103)*X(66)-JVS(1104)*X(67)-JVS(1105)*X(68)-JVS(1106)*X(69)-JVS(1107)*X(70)-JVS(1108)*X(71)&
            &-JVS(1109)*X(72)-JVS(1110)*X(74)-JVS(1111)*X(75)-JVS(1112)*X(76)-JVS(1113)*X(77)-JVS(1114)*X(78)-JVS(1115)&
            &*X(79)-JVS(1116)*X(80)-JVS(1117)*X(81)-JVS(1118)*X(82)-JVS(1119)*X(83)-JVS(1120)*X(84)-JVS(1121)*X(85)&
            &-JVS(1122)*X(86)-JVS(1123)*X(87)-JVS(1124)*X(88)-JVS(1125)*X(89)-JVS(1126)*X(90)-JVS(1127)*X(91)-JVS(1128)&
            &*X(92)
  X(94) = X(94)-JVS(1134)*X(29)-JVS(1135)*X(30)-JVS(1136)*X(32)-JVS(1137)*X(38)-JVS(1138)*X(39)-JVS(1139)*X(40)&
            &-JVS(1140)*X(55)-JVS(1141)*X(56)-JVS(1142)*X(57)-JVS(1143)*X(59)-JVS(1144)*X(61)-JVS(1145)*X(62)-JVS(1146)&
            &*X(63)-JVS(1147)*X(64)-JVS(1148)*X(65)-JVS(1149)*X(67)-JVS(1150)*X(68)-JVS(1151)*X(70)-JVS(1152)*X(71)&
            &-JVS(1153)*X(72)-JVS(1154)*X(74)-JVS(1155)*X(76)-JVS(1156)*X(77)-JVS(1157)*X(79)-JVS(1158)*X(80)-JVS(1159)&
            &*X(81)-JVS(1160)*X(82)-JVS(1161)*X(83)-JVS(1162)*X(84)-JVS(1163)*X(85)-JVS(1164)*X(86)-JVS(1165)*X(87)&
            &-JVS(1166)*X(88)-JVS(1167)*X(89)-JVS(1168)*X(90)-JVS(1169)*X(91)-JVS(1170)*X(92)-JVS(1171)*X(93)
  X(95) = X(95)-JVS(1176)*X(17)-JVS(1177)*X(18)-JVS(1178)*X(22)-JVS(1179)*X(26)-JVS(1180)*X(27)-JVS(1181)*X(29)&
            &-JVS(1182)*X(30)-JVS(1183)*X(31)-JVS(1184)*X(32)-JVS(1185)*X(38)-JVS(1186)*X(39)-JVS(1187)*X(40)-JVS(1188)&
            &*X(42)-JVS(1189)*X(44)-JVS(1190)*X(48)-JVS(1191)*X(55)-JVS(1192)*X(56)-JVS(1193)*X(57)-JVS(1194)*X(58)&
            &-JVS(1195)*X(59)-JVS(1196)*X(60)-JVS(1197)*X(61)-JVS(1198)*X(62)-JVS(1199)*X(63)-JVS(1200)*X(64)-JVS(1201)&
            &*X(65)-JVS(1202)*X(67)-JVS(1203)*X(68)-JVS(1204)*X(69)-JVS(1205)*X(70)-JVS(1206)*X(71)-JVS(1207)*X(72)&
            &-JVS(1208)*X(73)-JVS(1209)*X(74)-JVS(1210)*X(75)-JVS(1211)*X(76)-JVS(1212)*X(77)-JVS(1213)*X(78)-JVS(1214)&
            &*X(79)-JVS(1215)*X(80)-JVS(1216)*X(81)-JVS(1217)*X(82)-JVS(1218)*X(83)-JVS(1219)*X(84)-JVS(1220)*X(85)&
            &-JVS(1221)*X(86)-JVS(1222)*X(87)-JVS(1223)*X(88)-JVS(1224)*X(89)-JVS(1225)*X(90)-JVS(1226)*X(91)-JVS(1227)&
            &*X(92)-JVS(1228)*X(93)-JVS(1229)*X(94)
  X(96) = X(96)-JVS(1233)*X(17)-JVS(1234)*X(27)-JVS(1235)*X(30)-JVS(1236)*X(43)-JVS(1237)*X(44)-JVS(1238)*X(47)&
            &-JVS(1239)*X(48)-JVS(1240)*X(51)-JVS(1241)*X(53)-JVS(1242)*X(54)-JVS(1243)*X(55)-JVS(1244)*X(57)-JVS(1245)&
            &*X(58)-JVS(1246)*X(59)-JVS(1247)*X(60)-JVS(1248)*X(61)-JVS(1249)*X(62)-JVS(1250)*X(63)-JVS(1251)*X(65)&
            &-JVS(1252)*X(68)-JVS(1253)*X(69)-JVS(1254)*X(70)-JVS(1255)*X(71)-JVS(1256)*X(72)-JVS(1257)*X(73)-JVS(1258)&
            &*X(74)-JVS(1259)*X(76)-JVS(1260)*X(77)-JVS(1261)*X(78)-JVS(1262)*X(79)-JVS(1263)*X(80)-JVS(1264)*X(81)&
            &-JVS(1265)*X(82)-JVS(1266)*X(83)-JVS(1267)*X(84)-JVS(1268)*X(85)-JVS(1269)*X(86)-JVS(1270)*X(87)-JVS(1271)&
            &*X(88)-JVS(1272)*X(89)-JVS(1273)*X(90)-JVS(1274)*X(91)-JVS(1275)*X(92)-JVS(1276)*X(93)-JVS(1277)*X(94)&
            &-JVS(1278)*X(95)
  X(97) = X(97)-JVS(1281)*X(47)-JVS(1282)*X(50)-JVS(1283)*X(51)-JVS(1284)*X(54)-JVS(1285)*X(66)-JVS(1286)*X(70)&
            &-JVS(1287)*X(75)-JVS(1288)*X(76)-JVS(1289)*X(77)-JVS(1290)*X(79)-JVS(1291)*X(81)-JVS(1292)*X(82)-JVS(1293)&
            &*X(83)-JVS(1294)*X(84)-JVS(1295)*X(85)-JVS(1296)*X(86)-JVS(1297)*X(87)-JVS(1298)*X(88)-JVS(1299)*X(89)&
            &-JVS(1300)*X(90)-JVS(1301)*X(91)-JVS(1302)*X(92)-JVS(1303)*X(93)-JVS(1304)*X(94)-JVS(1305)*X(95)-JVS(1306)&
            &*X(96)
  X(97) = X(97)/JVS(1307)
  X(96) = (X(96)-JVS(1280)*X(97))/(JVS(1279))
  X(95) = (X(95)-JVS(1231)*X(96)-JVS(1232)*X(97))/(JVS(1230))
  X(94) = (X(94)-JVS(1173)*X(95)-JVS(1174)*X(96)-JVS(1175)*X(97))/(JVS(1172))
  X(93) = (X(93)-JVS(1130)*X(94)-JVS(1131)*X(95)-JVS(1132)*X(96)-JVS(1133)*X(97))/(JVS(1129))
  X(92) = (X(92)-JVS(1085)*X(93)-JVS(1086)*X(94)-JVS(1087)*X(95)-JVS(1088)*X(96)-JVS(1089)*X(97))/(JVS(1084))
  X(91) = (X(91)-JVS(1056)*X(92)-JVS(1057)*X(93)-JVS(1058)*X(94)-JVS(1059)*X(95)-JVS(1060)*X(96)-JVS(1061)*X(97))&
            &/(JVS(1055))
  X(90) = (X(90)-JVS(980)*X(91)-JVS(981)*X(92)-JVS(982)*X(93)-JVS(983)*X(94)-JVS(984)*X(95)-JVS(985)*X(96)-JVS(986)&
            &*X(97))/(JVS(979))
  X(89) = (X(89)-JVS(941)*X(90)-JVS(942)*X(91)-JVS(943)*X(92)-JVS(944)*X(93)-JVS(945)*X(94)-JVS(946)*X(95)-JVS(947)&
            &*X(96)-JVS(948)*X(97))/(JVS(940))
  X(88) = (X(88)-JVS(872)*X(89)-JVS(873)*X(90)-JVS(874)*X(91)-JVS(875)*X(92)-JVS(876)*X(93)-JVS(877)*X(94)-JVS(878)&
            &*X(95)-JVS(879)*X(96)-JVS(880)*X(97))/(JVS(871))
  X(87) = (X(87)-JVS(832)*X(88)-JVS(833)*X(89)-JVS(834)*X(90)-JVS(835)*X(91)-JVS(836)*X(92)-JVS(837)*X(93)-JVS(838)&
            &*X(94)-JVS(839)*X(95)-JVS(840)*X(96))/(JVS(831))
  X(86) = (X(86)-JVS(806)*X(88)-JVS(807)*X(89)-JVS(808)*X(90)-JVS(809)*X(91)-JVS(810)*X(92)-JVS(811)*X(93)-JVS(812)&
            &*X(94)-JVS(813)*X(95)-JVS(814)*X(96))/(JVS(805))
  X(85) = (X(85)-JVS(787)*X(88)-JVS(788)*X(89)-JVS(789)*X(90)-JVS(790)*X(91)-JVS(791)*X(93)-JVS(792)*X(94)-JVS(793)&
            &*X(95)-JVS(794)*X(96))/(JVS(786))
  X(84) = (X(84)-JVS(759)*X(85)-JVS(760)*X(88)-JVS(761)*X(89)-JVS(762)*X(90)-JVS(763)*X(91)-JVS(764)*X(93)-JVS(765)&
            &*X(94)-JVS(766)*X(96))/(JVS(758))
  X(83) = (X(83)-JVS(740)*X(84)-JVS(741)*X(85)-JVS(742)*X(88)-JVS(743)*X(89)-JVS(744)*X(90)-JVS(745)*X(91)-JVS(746)&
            &*X(93)-JVS(747)*X(94)-JVS(748)*X(96))/(JVS(739))
  X(82) = (X(82)-JVS(723)*X(85)-JVS(724)*X(88)-JVS(725)*X(89)-JVS(726)*X(90)-JVS(727)*X(91)-JVS(728)*X(93)-JVS(729)&
            &*X(94)-JVS(730)*X(96))/(JVS(722))
  X(81) = (X(81)-JVS(710)*X(85)-JVS(711)*X(88)-JVS(712)*X(89)-JVS(713)*X(90)-JVS(714)*X(91)-JVS(715)*X(93)-JVS(716)&
            &*X(94)-JVS(717)*X(96))/(JVS(709))
  X(80) = (X(80)-JVS(693)*X(81)-JVS(694)*X(82)-JVS(695)*X(85)-JVS(696)*X(86)-JVS(697)*X(87)-JVS(698)*X(88)-JVS(699)&
            &*X(89)-JVS(700)*X(90)-JVS(701)*X(91)-JVS(702)*X(92)-JVS(703)*X(93)-JVS(704)*X(94)-JVS(705)*X(95)-JVS(706)&
            &*X(96))/(JVS(692))
  X(79) = (X(79)-JVS(665)*X(85)-JVS(666)*X(88)-JVS(667)*X(89)-JVS(668)*X(90)-JVS(669)*X(91)-JVS(670)*X(93)-JVS(671)&
            &*X(94)-JVS(672)*X(95)-JVS(673)*X(96))/(JVS(664))
  X(78) = (X(78)-JVS(644)*X(79)-JVS(645)*X(81)-JVS(646)*X(82)-JVS(647)*X(83)-JVS(648)*X(84)-JVS(649)*X(85)-JVS(650)&
            &*X(86)-JVS(651)*X(87)-JVS(652)*X(88)-JVS(653)*X(89)-JVS(654)*X(90)-JVS(655)*X(91)-JVS(656)*X(92)-JVS(657)*X(93)&
            &-JVS(658)*X(94)-JVS(659)*X(95)-JVS(660)*X(96)-JVS(661)*X(97))/(JVS(643))
  X(77) = (X(77)-JVS(619)*X(81)-JVS(620)*X(82)-JVS(621)*X(88)-JVS(622)*X(89)-JVS(623)*X(90)-JVS(624)*X(91)-JVS(625)&
            &*X(93)-JVS(626)*X(94)-JVS(627)*X(96))/(JVS(618))
  X(76) = (X(76)-JVS(608)*X(81)-JVS(609)*X(82)-JVS(610)*X(88)-JVS(611)*X(89)-JVS(612)*X(90)-JVS(613)*X(91)-JVS(614)&
            &*X(93)-JVS(615)*X(94)-JVS(616)*X(96))/(JVS(607))
  X(75) = (X(75)-JVS(591)*X(76)-JVS(592)*X(77)-JVS(593)*X(79)-JVS(594)*X(83)-JVS(595)*X(84)-JVS(596)*X(85)-JVS(597)&
            &*X(86)-JVS(598)*X(88)-JVS(599)*X(89)-JVS(600)*X(90)-JVS(601)*X(91)-JVS(602)*X(93)-JVS(603)*X(94)-JVS(604)*X(95)&
            &-JVS(605)*X(96))/(JVS(590))
  X(74) = (X(74)-JVS(573)*X(79)-JVS(574)*X(85)-JVS(575)*X(89)-JVS(576)*X(90)-JVS(577)*X(91)-JVS(578)*X(93)-JVS(579)&
            &*X(94)-JVS(580)*X(95)-JVS(581)*X(96))/(JVS(572))
  X(73) = (X(73)-JVS(550)*X(74)-JVS(551)*X(76)-JVS(552)*X(77)-JVS(553)*X(79)-JVS(554)*X(80)-JVS(555)*X(81)-JVS(556)&
            &*X(82)-JVS(557)*X(83)-JVS(558)*X(84)-JVS(559)*X(85)-JVS(560)*X(86)-JVS(561)*X(87)-JVS(562)*X(88)-JVS(563)*X(89)&
            &-JVS(564)*X(90)-JVS(565)*X(91)-JVS(566)*X(93)-JVS(567)*X(94)-JVS(568)*X(95)-JVS(569)*X(96)-JVS(570)*X(97))&
            &/(JVS(549))
  X(72) = (X(72)-JVS(513)*X(85)-JVS(514)*X(88)-JVS(515)*X(89)-JVS(516)*X(90)-JVS(517)*X(91)-JVS(518)*X(93)-JVS(519)&
            &*X(94)-JVS(520)*X(95)-JVS(521)*X(96))/(JVS(512))
  X(71) = (X(71)-JVS(502)*X(85)-JVS(503)*X(88)-JVS(504)*X(89)-JVS(505)*X(90)-JVS(506)*X(91)-JVS(507)*X(93)-JVS(508)&
            &*X(94)-JVS(509)*X(95)-JVS(510)*X(96))/(JVS(501))
  X(70) = (X(70)-JVS(492)*X(81)-JVS(493)*X(85)-JVS(494)*X(89)-JVS(495)*X(90)-JVS(496)*X(91)-JVS(497)*X(93)-JVS(498)&
            &*X(94)-JVS(499)*X(96))/(JVS(491))
  X(69) = (X(69)-JVS(472)*X(71)-JVS(473)*X(72)-JVS(474)*X(74)-JVS(475)*X(79)-JVS(476)*X(81)-JVS(477)*X(82)-JVS(478)&
            &*X(85)-JVS(479)*X(86)-JVS(480)*X(87)-JVS(481)*X(88)-JVS(482)*X(89)-JVS(483)*X(90)-JVS(484)*X(91)-JVS(485)*X(93)&
            &-JVS(486)*X(94)-JVS(487)*X(95)-JVS(488)*X(96))/(JVS(471))
  X(68) = (X(68)-JVS(452)*X(85)-JVS(453)*X(88)-JVS(454)*X(89)-JVS(455)*X(90)-JVS(456)*X(91)-JVS(457)*X(93)-JVS(458)&
            &*X(94)-JVS(459)*X(95)-JVS(460)*X(96))/(JVS(451))
  X(67) = (X(67)-JVS(441)*X(85)-JVS(442)*X(88)-JVS(443)*X(89)-JVS(444)*X(90)-JVS(445)*X(91)-JVS(446)*X(93)-JVS(447)&
            &*X(94)-JVS(448)*X(96))/(JVS(440))
  X(66) = (X(66)-JVS(418)*X(70)-JVS(419)*X(75)-JVS(420)*X(76)-JVS(421)*X(77)-JVS(422)*X(79)-JVS(423)*X(81)-JVS(424)&
            &*X(82)-JVS(425)*X(83)-JVS(426)*X(84)-JVS(427)*X(85)-JVS(428)*X(86)-JVS(429)*X(87)-JVS(430)*X(88)-JVS(431)*X(89)&
            &-JVS(432)*X(90)-JVS(433)*X(91)-JVS(434)*X(92)-JVS(435)*X(93)-JVS(436)*X(94)-JVS(437)*X(96))/(JVS(417))
  X(65) = (X(65)-JVS(403)*X(85)-JVS(404)*X(89)-JVS(405)*X(90)-JVS(406)*X(91)-JVS(407)*X(93)-JVS(408)*X(94)-JVS(409)&
            &*X(96))/(JVS(402))
  X(64) = (X(64)-JVS(392)*X(67)-JVS(393)*X(85)-JVS(394)*X(88)-JVS(395)*X(89)-JVS(396)*X(90)-JVS(397)*X(91)-JVS(398)&
            &*X(93)-JVS(399)*X(94)-JVS(400)*X(96))/(JVS(391))
  X(63) = (X(63)-JVS(378)*X(85)-JVS(379)*X(89)-JVS(380)*X(90)-JVS(381)*X(91)-JVS(382)*X(93)-JVS(383)*X(94)-JVS(384)&
            &*X(96))/(JVS(377))
  X(62) = (X(62)-JVS(368)*X(85)-JVS(369)*X(88)-JVS(370)*X(89)-JVS(371)*X(90)-JVS(372)*X(91)-JVS(373)*X(93)-JVS(374)&
            &*X(94)-JVS(375)*X(96))/(JVS(367))
  X(61) = (X(61)-JVS(359)*X(82)-JVS(360)*X(89)-JVS(361)*X(90)-JVS(362)*X(91)-JVS(363)*X(93)-JVS(364)*X(94)-JVS(365)&
            &*X(96))/(JVS(358))
  X(60) = (X(60)-JVS(344)*X(68)-JVS(345)*X(71)-JVS(346)*X(72)-JVS(347)*X(79)-JVS(348)*X(81)-JVS(349)*X(82)-JVS(350)&
            &*X(85)-JVS(351)*X(86)-JVS(352)*X(88)-JVS(353)*X(90)-JVS(354)*X(91)-JVS(355)*X(93)-JVS(356)*X(94)-JVS(357)&
            &*X(96))/(JVS(343))
  X(59) = (X(59)-JVS(338)*X(85)-JVS(339)*X(88)-JVS(340)*X(91)-JVS(341)*X(96))/(JVS(337))
  X(58) = (X(58)-JVS(332)*X(74)-JVS(333)*X(85)-JVS(334)*X(91)-JVS(335)*X(95)-JVS(336)*X(96))/(JVS(331))
  X(57) = (X(57)-JVS(325)*X(85)-JVS(326)*X(88)-JVS(327)*X(89)-JVS(328)*X(91)-JVS(329)*X(95)-JVS(330)*X(96))/(JVS(324))
  X(56) = (X(56)-JVS(313)*X(59)-JVS(314)*X(64)-JVS(315)*X(85)-JVS(316)*X(89)-JVS(317)*X(91)-JVS(318)*X(94)-JVS(319)&
            &*X(95))/(JVS(312))
  X(55) = (X(55)-JVS(307)*X(85)-JVS(308)*X(91)-JVS(309)*X(96))/(JVS(306))
  X(54) = (X(54)-JVS(303)*X(85)-JVS(304)*X(91)-JVS(305)*X(96))/(JVS(302))
  X(53) = (X(53)-JVS(298)*X(85)-JVS(299)*X(88)-JVS(300)*X(91)-JVS(301)*X(96))/(JVS(297))
  X(52) = (X(52)-JVS(277)*X(53)-JVS(278)*X(54)-JVS(279)*X(55)-JVS(280)*X(56)-JVS(281)*X(58)-JVS(282)*X(59)-JVS(283)&
            &*X(60)-JVS(284)*X(64)-JVS(285)*X(67)-JVS(286)*X(69)-JVS(287)*X(73)-JVS(288)*X(78)-JVS(289)*X(79)-JVS(290)*X(81)&
            &-JVS(291)*X(82)-JVS(292)*X(85)-JVS(293)*X(88)-JVS(294)*X(91)-JVS(295)*X(94)-JVS(296)*X(96))/(JVS(276))
  X(51) = (X(51)-JVS(267)*X(85)-JVS(268)*X(91)-JVS(269)*X(96))/(JVS(266))
  X(50) = (X(50)-JVS(260)*X(54)-JVS(261)*X(81)-JVS(262)*X(82)-JVS(263)*X(85)-JVS(264)*X(88)-JVS(265)*X(91))/(JVS(259))
  X(49) = (X(49)-JVS(255)*X(81)-JVS(256)*X(82)-JVS(257)*X(88)-JVS(258)*X(91))/(JVS(254))
  X(48) = (X(48)-JVS(248)*X(58)-JVS(249)*X(85)-JVS(250)*X(91)-JVS(251)*X(93)-JVS(252)*X(95)-JVS(253)*X(96))/(JVS(247))
  X(47) = (X(47)-JVS(244)*X(85)-JVS(245)*X(91)-JVS(246)*X(96))/(JVS(243))
  X(46) = (X(46)-JVS(239)*X(81)-JVS(240)*X(82)-JVS(241)*X(88)-JVS(242)*X(91))/(JVS(238))
  X(45) = (X(45)-JVS(233)*X(79)-JVS(234)*X(85)-JVS(235)*X(89)-JVS(236)*X(91)-JVS(237)*X(93))/(JVS(232))
  X(44) = (X(44)-JVS(222)*X(57)-JVS(223)*X(60)-JVS(224)*X(69)-JVS(225)*X(73)-JVS(226)*X(78)-JVS(227)*X(79)-JVS(228)&
            &*X(89)-JVS(229)*X(91)-JVS(230)*X(95)-JVS(231)*X(96))/(JVS(221))
  X(43) = (X(43)-JVS(218)*X(85)-JVS(219)*X(91)-JVS(220)*X(96))/(JVS(217))
  X(42) = (X(42)-JVS(212)*X(59)-JVS(213)*X(67)-JVS(214)*X(91)-JVS(215)*X(94)-JVS(216)*X(96))/(JVS(211))
  X(41) = (X(41)-JVS(207)*X(81)-JVS(208)*X(82)-JVS(209)*X(88)-JVS(210)*X(91))/(JVS(206))
  X(40) = (X(40)-JVS(202)*X(57)-JVS(203)*X(85)-JVS(204)*X(91)-JVS(205)*X(95))/(JVS(201))
  X(39) = (X(39)-JVS(197)*X(85)-JVS(198)*X(88)-JVS(199)*X(91)-JVS(200)*X(95))/(JVS(196))
  X(38) = (X(38)-JVS(191)*X(85)-JVS(192)*X(88)-JVS(193)*X(91)-JVS(194)*X(95))/(JVS(190))
  X(37) = (X(37)-JVS(184)*X(46)-JVS(185)*X(58)-JVS(186)*X(87)-JVS(187)*X(90)-JVS(188)*X(91))/(JVS(183))
  X(36) = (X(36)-JVS(177)*X(53)-JVS(178)*X(59)-JVS(179)*X(81)-JVS(180)*X(82)-JVS(181)*X(88)-JVS(182)*X(91))/(JVS(176))
  X(35) = (X(35)-JVS(170)*X(54)-JVS(171)*X(81)-JVS(172)*X(82)-JVS(173)*X(85)-JVS(174)*X(88)-JVS(175)*X(91))/(JVS(169))
  X(34) = (X(34)-JVS(158)*X(35)-JVS(159)*X(41)-JVS(160)*X(46)-JVS(161)*X(49)-JVS(162)*X(50)-JVS(163)*X(59)-JVS(164)&
            &*X(73)-JVS(165)*X(78)-JVS(166)*X(81)-JVS(167)*X(88)-JVS(168)*X(91))/(JVS(157))
  X(33) = (X(33)-JVS(145)*X(47)-JVS(146)*X(51)-JVS(147)*X(54)-JVS(148)*X(55)-JVS(149)*X(59)-JVS(150)*X(81)-JVS(151)&
            &*X(82)-JVS(152)*X(85)-JVS(153)*X(89)-JVS(154)*X(91))/(JVS(144))
  X(32) = (X(32)-JVS(140)*X(85)-JVS(141)*X(88)-JVS(142)*X(89)-JVS(143)*X(94))/(JVS(139))
  X(31) = (X(31)-JVS(134)*X(57)-JVS(135)*X(89)-JVS(136)*X(91)-JVS(137)*X(95)-JVS(138)*X(96))/(JVS(133))
  X(30) = (X(30)-JVS(128)*X(55)-JVS(129)*X(85)-JVS(130)*X(94)-JVS(131)*X(95)-JVS(132)*X(96))/(JVS(127))
  X(29) = (X(29)-JVS(120)*X(38)-JVS(121)*X(39)-JVS(122)*X(40)-JVS(123)*X(91)-JVS(124)*X(94)-JVS(125)*X(95))/(JVS(119))
  X(28) = (X(28)-JVS(114)*X(42)-JVS(115)*X(56)-JVS(116)*X(91)-JVS(117)*X(94)-JVS(118)*X(95))/(JVS(113))
  X(27) = (X(27)-JVS(109)*X(89)-JVS(110)*X(91)-JVS(111)*X(95))/(JVS(108))
  X(26) = (X(26)-JVS(105)*X(56)-JVS(106)*X(91)-JVS(107)*X(95))/(JVS(104))
  X(25) = (X(25)-JVS(101)*X(89)-JVS(102)*X(90)-JVS(103)*X(91))/(JVS(100))
  X(24) = (X(24)-JVS(98)*X(88)-JVS(99)*X(91))/(JVS(97))
  X(23) = (X(23)-JVS(95)*X(88)-JVS(96)*X(91))/(JVS(94))
  X(22) = (X(22)-JVS(92)*X(42)-JVS(93)*X(91))/(JVS(91))
  X(21) = (X(21)-JVS(89)*X(79)-JVS(90)*X(91))/(JVS(88))
  X(20) = (X(20)-JVS(85)*X(56)-JVS(86)*X(89)-JVS(87)*X(91))/(JVS(84))
  X(19) = (X(19)-JVS(82)*X(32)-JVS(83)*X(89))/(JVS(81))
  X(18) = (X(18)-JVS(79)*X(88)-JVS(80)*X(95))/(JVS(78))
  X(17) = (X(17)-JVS(76)*X(95)-JVS(77)*X(96))/(JVS(75))
  X(16) = (X(16)-JVS(72)*X(67)-JVS(73)*X(89)-JVS(74)*X(91))/(JVS(71))
  X(15) = (X(15)-JVS(70)*X(91))/(JVS(69))
  X(14) = (X(14)-JVS(68)*X(85))/(JVS(67))
  X(13) = (X(13)-JVS(66)*X(32))/(JVS(65))
  X(12) = (X(12)-JVS(64)*X(91))/(JVS(63))
  X(11) = (X(11)-JVS(61)*X(12)-JVS(62)*X(91))/(JVS(60))
  X(10) = (X(10)-JVS(58)*X(11)-JVS(59)*X(91))/(JVS(57))
  X(9) = (X(9)-JVS(55)*X(10)-JVS(56)*X(91))/(JVS(54))
  X(8) = (X(8)-JVS(53)*X(91))/(JVS(52))
  X(7) = (X(7)-JVS(50)*X(8)-JVS(51)*X(91))/(JVS(49))
  X(6) = (X(6)-JVS(47)*X(7)-JVS(48)*X(91))/(JVS(46))
  X(5) = (X(5)-JVS(44)*X(6)-JVS(45)*X(91))/(JVS(43))
  X(4) = (X(4)-JVS(21)*X(43)-JVS(22)*X(51)-JVS(23)*X(54)-JVS(24)*X(61)-JVS(25)*X(62)-JVS(26)*X(67)-JVS(27)*X(70)-JVS(28)&
           &*X(74)-JVS(29)*X(76)-JVS(30)*X(77)-JVS(31)*X(79)-JVS(32)*X(81)-JVS(33)*X(82)-JVS(34)*X(83)-JVS(35)*X(84)-JVS(36)&
           &*X(85)-JVS(37)*X(86)-JVS(38)*X(87)-JVS(39)*X(89)-JVS(40)*X(90)-JVS(41)*X(93)-JVS(42)*X(97))/(JVS(20))
  X(3) = (X(3)-JVS(8)*X(43)-JVS(9)*X(49)-JVS(10)*X(51)-JVS(11)*X(53)-JVS(12)*X(55)-JVS(13)*X(58)-JVS(14)*X(59)-JVS(15)&
           &*X(64)-JVS(16)*X(79)-JVS(17)*X(82)-JVS(18)*X(85)-JVS(19)*X(91))/(JVS(7))
  X(2) = (X(2)-JVS(5)*X(52)-JVS(6)*X(91))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(15)-JVS(3)*X(91))/(JVS(1))
      
END SUBROUTINE racm_soa_vbs_het_KppSolve
























      SUBROUTINE racm_soa_vbs_het_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE racm_soa_vbs_het_WCOPY



      SUBROUTINE racm_soa_vbs_het_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE racm_soa_vbs_het_WAXPY




      SUBROUTINE racm_soa_vbs_het_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE racm_soa_vbs_het_WSCAL


      REAL(kind=dp) FUNCTION racm_soa_vbs_het_WLAMCH( C )








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
          CALL racm_soa_vbs_het_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      racm_soa_vbs_het_WLAMCH = Eps

      END FUNCTION racm_soa_vbs_het_WLAMCH
     
      SUBROUTINE racm_soa_vbs_het_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE racm_soa_vbs_het_WLAMCH_ADD




      SUBROUTINE racm_soa_vbs_het_SET2ZERO(N,Y)




      
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

      END SUBROUTINE racm_soa_vbs_het_SET2ZERO



      REAL(kind=dp) FUNCTION racm_soa_vbs_het_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      racm_soa_vbs_het_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        racm_soa_vbs_het_WDOT = racm_soa_vbs_het_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         racm_soa_vbs_het_WDOT = racm_soa_vbs_het_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          racm_soa_vbs_het_WDOT = racm_soa_vbs_het_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        racm_soa_vbs_het_WDOT = racm_soa_vbs_het_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION racm_soa_vbs_het_WDOT                                          




   SUBROUTINE decomp_racm_soa_vbs_het( JVS, IER )
   
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
   W( 15 ) = JVS( 2 )
   W( 91 ) = JVS( 3 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 15 )
  JVS( 3) = W( 91 )
  IF ( ABS(  JVS( 4 )) < TINY(a) ) THEN
         IER = 2                                       
         RETURN
  END IF
   W( 2 ) = JVS( 4 )
   W( 52 ) = JVS( 5 )
   W( 91 ) = JVS( 6 )
  JVS( 4) = W( 2 )
  JVS( 5) = W( 52 )
  JVS( 6) = W( 91 )
  IF ( ABS(  JVS( 7 )) < TINY(a) ) THEN
         IER = 3                                       
         RETURN
  END IF
   W( 3 ) = JVS( 7 )
   W( 43 ) = JVS( 8 )
   W( 49 ) = JVS( 9 )
   W( 51 ) = JVS( 10 )
   W( 53 ) = JVS( 11 )
   W( 55 ) = JVS( 12 )
   W( 58 ) = JVS( 13 )
   W( 59 ) = JVS( 14 )
   W( 64 ) = JVS( 15 )
   W( 79 ) = JVS( 16 )
   W( 82 ) = JVS( 17 )
   W( 85 ) = JVS( 18 )
   W( 91 ) = JVS( 19 )
  JVS( 7) = W( 3 )
  JVS( 8) = W( 43 )
  JVS( 9) = W( 49 )
  JVS( 10) = W( 51 )
  JVS( 11) = W( 53 )
  JVS( 12) = W( 55 )
  JVS( 13) = W( 58 )
  JVS( 14) = W( 59 )
  JVS( 15) = W( 64 )
  JVS( 16) = W( 79 )
  JVS( 17) = W( 82 )
  JVS( 18) = W( 85 )
  JVS( 19) = W( 91 )
  IF ( ABS(  JVS( 20 )) < TINY(a) ) THEN
         IER = 4                                       
         RETURN
  END IF
   W( 4 ) = JVS( 20 )
   W( 43 ) = JVS( 21 )
   W( 51 ) = JVS( 22 )
   W( 54 ) = JVS( 23 )
   W( 61 ) = JVS( 24 )
   W( 62 ) = JVS( 25 )
   W( 67 ) = JVS( 26 )
   W( 70 ) = JVS( 27 )
   W( 74 ) = JVS( 28 )
   W( 76 ) = JVS( 29 )
   W( 77 ) = JVS( 30 )
   W( 79 ) = JVS( 31 )
   W( 81 ) = JVS( 32 )
   W( 82 ) = JVS( 33 )
   W( 83 ) = JVS( 34 )
   W( 84 ) = JVS( 35 )
   W( 85 ) = JVS( 36 )
   W( 86 ) = JVS( 37 )
   W( 87 ) = JVS( 38 )
   W( 89 ) = JVS( 39 )
   W( 90 ) = JVS( 40 )
   W( 93 ) = JVS( 41 )
   W( 97 ) = JVS( 42 )
  JVS( 20) = W( 4 )
  JVS( 21) = W( 43 )
  JVS( 22) = W( 51 )
  JVS( 23) = W( 54 )
  JVS( 24) = W( 61 )
  JVS( 25) = W( 62 )
  JVS( 26) = W( 67 )
  JVS( 27) = W( 70 )
  JVS( 28) = W( 74 )
  JVS( 29) = W( 76 )
  JVS( 30) = W( 77 )
  JVS( 31) = W( 79 )
  JVS( 32) = W( 81 )
  JVS( 33) = W( 82 )
  JVS( 34) = W( 83 )
  JVS( 35) = W( 84 )
  JVS( 36) = W( 85 )
  JVS( 37) = W( 86 )
  JVS( 38) = W( 87 )
  JVS( 39) = W( 89 )
  JVS( 40) = W( 90 )
  JVS( 41) = W( 93 )
  JVS( 42) = W( 97 )
  IF ( ABS(  JVS( 43 )) < TINY(a) ) THEN
         IER = 5                                       
         RETURN
  END IF
   W( 5 ) = JVS( 43 )
   W( 6 ) = JVS( 44 )
   W( 91 ) = JVS( 45 )
  JVS( 43) = W( 5 )
  JVS( 44) = W( 6 )
  JVS( 45) = W( 91 )
  IF ( ABS(  JVS( 46 )) < TINY(a) ) THEN
         IER = 6                                       
         RETURN
  END IF
   W( 6 ) = JVS( 46 )
   W( 7 ) = JVS( 47 )
   W( 91 ) = JVS( 48 )
  JVS( 46) = W( 6 )
  JVS( 47) = W( 7 )
  JVS( 48) = W( 91 )
  IF ( ABS(  JVS( 49 )) < TINY(a) ) THEN
         IER = 7                                       
         RETURN
  END IF
   W( 7 ) = JVS( 49 )
   W( 8 ) = JVS( 50 )
   W( 91 ) = JVS( 51 )
  JVS( 49) = W( 7 )
  JVS( 50) = W( 8 )
  JVS( 51) = W( 91 )
  IF ( ABS(  JVS( 52 )) < TINY(a) ) THEN
         IER = 8                                       
         RETURN
  END IF
   W( 8 ) = JVS( 52 )
   W( 91 ) = JVS( 53 )
  JVS( 52) = W( 8 )
  JVS( 53) = W( 91 )
  IF ( ABS(  JVS( 54 )) < TINY(a) ) THEN
         IER = 9                                       
         RETURN
  END IF
   W( 9 ) = JVS( 54 )
   W( 10 ) = JVS( 55 )
   W( 91 ) = JVS( 56 )
  JVS( 54) = W( 9 )
  JVS( 55) = W( 10 )
  JVS( 56) = W( 91 )
  IF ( ABS(  JVS( 57 )) < TINY(a) ) THEN
         IER = 10                                      
         RETURN
  END IF
   W( 10 ) = JVS( 57 )
   W( 11 ) = JVS( 58 )
   W( 91 ) = JVS( 59 )
  JVS( 57) = W( 10 )
  JVS( 58) = W( 11 )
  JVS( 59) = W( 91 )
  IF ( ABS(  JVS( 60 )) < TINY(a) ) THEN
         IER = 11                                      
         RETURN
  END IF
   W( 11 ) = JVS( 60 )
   W( 12 ) = JVS( 61 )
   W( 91 ) = JVS( 62 )
  JVS( 60) = W( 11 )
  JVS( 61) = W( 12 )
  JVS( 62) = W( 91 )
  IF ( ABS(  JVS( 63 )) < TINY(a) ) THEN
         IER = 12                                      
         RETURN
  END IF
   W( 12 ) = JVS( 63 )
   W( 91 ) = JVS( 64 )
  JVS( 63) = W( 12 )
  JVS( 64) = W( 91 )
  IF ( ABS(  JVS( 65 )) < TINY(a) ) THEN
         IER = 13                                      
         RETURN
  END IF
   W( 13 ) = JVS( 65 )
   W( 32 ) = JVS( 66 )
  JVS( 65) = W( 13 )
  JVS( 66) = W( 32 )
  IF ( ABS(  JVS( 67 )) < TINY(a) ) THEN
         IER = 14                                      
         RETURN
  END IF
   W( 14 ) = JVS( 67 )
   W( 85 ) = JVS( 68 )
  JVS( 67) = W( 14 )
  JVS( 68) = W( 85 )
  IF ( ABS(  JVS( 69 )) < TINY(a) ) THEN
         IER = 15                                      
         RETURN
  END IF
   W( 15 ) = JVS( 69 )
   W( 91 ) = JVS( 70 )
  JVS( 69) = W( 15 )
  JVS( 70) = W( 91 )
  IF ( ABS(  JVS( 71 )) < TINY(a) ) THEN
         IER = 16                                      
         RETURN
  END IF
   W( 16 ) = JVS( 71 )
   W( 67 ) = JVS( 72 )
   W( 89 ) = JVS( 73 )
   W( 91 ) = JVS( 74 )
  JVS( 71) = W( 16 )
  JVS( 72) = W( 67 )
  JVS( 73) = W( 89 )
  JVS( 74) = W( 91 )
  IF ( ABS(  JVS( 75 )) < TINY(a) ) THEN
         IER = 17                                      
         RETURN
  END IF
   W( 17 ) = JVS( 75 )
   W( 95 ) = JVS( 76 )
   W( 96 ) = JVS( 77 )
  JVS( 75) = W( 17 )
  JVS( 76) = W( 95 )
  JVS( 77) = W( 96 )
  IF ( ABS(  JVS( 78 )) < TINY(a) ) THEN
         IER = 18                                      
         RETURN
  END IF
   W( 18 ) = JVS( 78 )
   W( 88 ) = JVS( 79 )
   W( 95 ) = JVS( 80 )
  JVS( 78) = W( 18 )
  JVS( 79) = W( 88 )
  JVS( 80) = W( 95 )
  IF ( ABS(  JVS( 81 )) < TINY(a) ) THEN
         IER = 19                                      
         RETURN
  END IF
   W( 19 ) = JVS( 81 )
   W( 32 ) = JVS( 82 )
   W( 89 ) = JVS( 83 )
  JVS( 81) = W( 19 )
  JVS( 82) = W( 32 )
  JVS( 83) = W( 89 )
  IF ( ABS(  JVS( 84 )) < TINY(a) ) THEN
         IER = 20                                      
         RETURN
  END IF
   W( 20 ) = JVS( 84 )
   W( 56 ) = JVS( 85 )
   W( 89 ) = JVS( 86 )
   W( 91 ) = JVS( 87 )
  JVS( 84) = W( 20 )
  JVS( 85) = W( 56 )
  JVS( 86) = W( 89 )
  JVS( 87) = W( 91 )
  IF ( ABS(  JVS( 88 )) < TINY(a) ) THEN
         IER = 21                                      
         RETURN
  END IF
   W( 21 ) = JVS( 88 )
   W( 79 ) = JVS( 89 )
   W( 91 ) = JVS( 90 )
  JVS( 88) = W( 21 )
  JVS( 89) = W( 79 )
  JVS( 90) = W( 91 )
  IF ( ABS(  JVS( 91 )) < TINY(a) ) THEN
         IER = 22                                      
         RETURN
  END IF
   W( 22 ) = JVS( 91 )
   W( 42 ) = JVS( 92 )
   W( 91 ) = JVS( 93 )
  JVS( 91) = W( 22 )
  JVS( 92) = W( 42 )
  JVS( 93) = W( 91 )
  IF ( ABS(  JVS( 94 )) < TINY(a) ) THEN
         IER = 23                                      
         RETURN
  END IF
   W( 23 ) = JVS( 94 )
   W( 88 ) = JVS( 95 )
   W( 91 ) = JVS( 96 )
  JVS( 94) = W( 23 )
  JVS( 95) = W( 88 )
  JVS( 96) = W( 91 )
  IF ( ABS(  JVS( 97 )) < TINY(a) ) THEN
         IER = 24                                      
         RETURN
  END IF
   W( 24 ) = JVS( 97 )
   W( 88 ) = JVS( 98 )
   W( 91 ) = JVS( 99 )
  JVS( 97) = W( 24 )
  JVS( 98) = W( 88 )
  JVS( 99) = W( 91 )
  IF ( ABS(  JVS( 100 )) < TINY(a) ) THEN
         IER = 25                                      
         RETURN
  END IF
   W( 25 ) = JVS( 100 )
   W( 89 ) = JVS( 101 )
   W( 90 ) = JVS( 102 )
   W( 91 ) = JVS( 103 )
  JVS( 100) = W( 25 )
  JVS( 101) = W( 89 )
  JVS( 102) = W( 90 )
  JVS( 103) = W( 91 )
  IF ( ABS(  JVS( 104 )) < TINY(a) ) THEN
         IER = 26                                      
         RETURN
  END IF
   W( 26 ) = JVS( 104 )
   W( 56 ) = JVS( 105 )
   W( 91 ) = JVS( 106 )
   W( 95 ) = JVS( 107 )
  JVS( 104) = W( 26 )
  JVS( 105) = W( 56 )
  JVS( 106) = W( 91 )
  JVS( 107) = W( 95 )
  IF ( ABS(  JVS( 108 )) < TINY(a) ) THEN
         IER = 27                                      
         RETURN
  END IF
   W( 27 ) = JVS( 108 )
   W( 89 ) = JVS( 109 )
   W( 91 ) = JVS( 110 )
   W( 95 ) = JVS( 111 )
  JVS( 108) = W( 27 )
  JVS( 109) = W( 89 )
  JVS( 110) = W( 91 )
  JVS( 111) = W( 95 )
  IF ( ABS(  JVS( 113 )) < TINY(a) ) THEN
         IER = 28                                      
         RETURN
  END IF
   W( 26 ) = JVS( 112 )
   W( 28 ) = JVS( 113 )
   W( 42 ) = JVS( 114 )
   W( 56 ) = JVS( 115 )
   W( 91 ) = JVS( 116 )
   W( 94 ) = JVS( 117 )
   W( 95 ) = JVS( 118 )
  a = -W( 26 ) / JVS(          104  )
  W( 26 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 105 )
  W( 91 ) = W( 91 ) + a*JVS( 106 )
  W( 95 ) = W( 95 ) + a*JVS( 107 )
  JVS( 112) = W( 26 )
  JVS( 113) = W( 28 )
  JVS( 114) = W( 42 )
  JVS( 115) = W( 56 )
  JVS( 116) = W( 91 )
  JVS( 117) = W( 94 )
  JVS( 118) = W( 95 )
  IF ( ABS(  JVS( 119 )) < TINY(a) ) THEN
         IER = 29                                      
         RETURN
  END IF
   W( 29 ) = JVS( 119 )
   W( 38 ) = JVS( 120 )
   W( 39 ) = JVS( 121 )
   W( 40 ) = JVS( 122 )
   W( 91 ) = JVS( 123 )
   W( 94 ) = JVS( 124 )
   W( 95 ) = JVS( 125 )
  JVS( 119) = W( 29 )
  JVS( 120) = W( 38 )
  JVS( 121) = W( 39 )
  JVS( 122) = W( 40 )
  JVS( 123) = W( 91 )
  JVS( 124) = W( 94 )
  JVS( 125) = W( 95 )
  IF ( ABS(  JVS( 127 )) < TINY(a) ) THEN
         IER = 30                                      
         RETURN
  END IF
   W( 14 ) = JVS( 126 )
   W( 30 ) = JVS( 127 )
   W( 55 ) = JVS( 128 )
   W( 85 ) = JVS( 129 )
   W( 94 ) = JVS( 130 )
   W( 95 ) = JVS( 131 )
   W( 96 ) = JVS( 132 )
  a = -W( 14 ) / JVS(           67  )
  W( 14 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 68 )
  JVS( 126) = W( 14 )
  JVS( 127) = W( 30 )
  JVS( 128) = W( 55 )
  JVS( 129) = W( 85 )
  JVS( 130) = W( 94 )
  JVS( 131) = W( 95 )
  JVS( 132) = W( 96 )
  IF ( ABS(  JVS( 133 )) < TINY(a) ) THEN
         IER = 31                                      
         RETURN
  END IF
   W( 31 ) = JVS( 133 )
   W( 57 ) = JVS( 134 )
   W( 89 ) = JVS( 135 )
   W( 91 ) = JVS( 136 )
   W( 95 ) = JVS( 137 )
   W( 96 ) = JVS( 138 )
  JVS( 133) = W( 31 )
  JVS( 134) = W( 57 )
  JVS( 135) = W( 89 )
  JVS( 136) = W( 91 )
  JVS( 137) = W( 95 )
  JVS( 138) = W( 96 )
  IF ( ABS(  JVS( 139 )) < TINY(a) ) THEN
         IER = 32                                      
         RETURN
  END IF
   W( 32 ) = JVS( 139 )
   W( 85 ) = JVS( 140 )
   W( 88 ) = JVS( 141 )
   W( 89 ) = JVS( 142 )
   W( 94 ) = JVS( 143 )
  JVS( 139) = W( 32 )
  JVS( 140) = W( 85 )
  JVS( 141) = W( 88 )
  JVS( 142) = W( 89 )
  JVS( 143) = W( 94 )
  IF ( ABS(  JVS( 144 )) < TINY(a) ) THEN
         IER = 33                                      
         RETURN
  END IF
   W( 33 ) = JVS( 144 )
   W( 47 ) = JVS( 145 )
   W( 51 ) = JVS( 146 )
   W( 54 ) = JVS( 147 )
   W( 55 ) = JVS( 148 )
   W( 59 ) = JVS( 149 )
   W( 81 ) = JVS( 150 )
   W( 82 ) = JVS( 151 )
   W( 85 ) = JVS( 152 )
   W( 89 ) = JVS( 153 )
   W( 91 ) = JVS( 154 )
  JVS( 144) = W( 33 )
  JVS( 145) = W( 47 )
  JVS( 146) = W( 51 )
  JVS( 147) = W( 54 )
  JVS( 148) = W( 55 )
  JVS( 149) = W( 59 )
  JVS( 150) = W( 81 )
  JVS( 151) = W( 82 )
  JVS( 152) = W( 85 )
  JVS( 153) = W( 89 )
  JVS( 154) = W( 91 )
  IF ( ABS(  JVS( 157 )) < TINY(a) ) THEN
         IER = 34                                      
         RETURN
  END IF
   W( 23 ) = JVS( 155 )
   W( 24 ) = JVS( 156 )
   W( 34 ) = JVS( 157 )
   W( 35 ) = JVS( 158 )
   W( 41 ) = JVS( 159 )
   W( 46 ) = JVS( 160 )
   W( 49 ) = JVS( 161 )
   W( 50 ) = JVS( 162 )
   W( 59 ) = JVS( 163 )
   W( 73 ) = JVS( 164 )
   W( 78 ) = JVS( 165 )
   W( 81 ) = JVS( 166 )
   W( 88 ) = JVS( 167 )
   W( 91 ) = JVS( 168 )
  a = -W( 23 ) / JVS(           94  )
  W( 23 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 95 )
  W( 91 ) = W( 91 ) + a*JVS( 96 )
  a = -W( 24 ) / JVS(           97  )
  W( 24 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 98 )
  W( 91 ) = W( 91 ) + a*JVS( 99 )
  JVS( 155) = W( 23 )
  JVS( 156) = W( 24 )
  JVS( 157) = W( 34 )
  JVS( 158) = W( 35 )
  JVS( 159) = W( 41 )
  JVS( 160) = W( 46 )
  JVS( 161) = W( 49 )
  JVS( 162) = W( 50 )
  JVS( 163) = W( 59 )
  JVS( 164) = W( 73 )
  JVS( 165) = W( 78 )
  JVS( 166) = W( 81 )
  JVS( 167) = W( 88 )
  JVS( 168) = W( 91 )
  IF ( ABS(  JVS( 169 )) < TINY(a) ) THEN
         IER = 35                                      
         RETURN
  END IF
   W( 35 ) = JVS( 169 )
   W( 54 ) = JVS( 170 )
   W( 81 ) = JVS( 171 )
   W( 82 ) = JVS( 172 )
   W( 85 ) = JVS( 173 )
   W( 88 ) = JVS( 174 )
   W( 91 ) = JVS( 175 )
  JVS( 169) = W( 35 )
  JVS( 170) = W( 54 )
  JVS( 171) = W( 81 )
  JVS( 172) = W( 82 )
  JVS( 173) = W( 85 )
  JVS( 174) = W( 88 )
  JVS( 175) = W( 91 )
  IF ( ABS(  JVS( 176 )) < TINY(a) ) THEN
         IER = 36                                      
         RETURN
  END IF
   W( 36 ) = JVS( 176 )
   W( 53 ) = JVS( 177 )
   W( 59 ) = JVS( 178 )
   W( 81 ) = JVS( 179 )
   W( 82 ) = JVS( 180 )
   W( 88 ) = JVS( 181 )
   W( 91 ) = JVS( 182 )
  JVS( 176) = W( 36 )
  JVS( 177) = W( 53 )
  JVS( 178) = W( 59 )
  JVS( 179) = W( 81 )
  JVS( 180) = W( 82 )
  JVS( 181) = W( 88 )
  JVS( 182) = W( 91 )
  IF ( ABS(  JVS( 183 )) < TINY(a) ) THEN
         IER = 37                                      
         RETURN
  END IF
   W( 37 ) = JVS( 183 )
   W( 46 ) = JVS( 184 )
   W( 58 ) = JVS( 185 )
   W( 87 ) = JVS( 186 )
   W( 90 ) = JVS( 187 )
   W( 91 ) = JVS( 188 )
  JVS( 183) = W( 37 )
  JVS( 184) = W( 46 )
  JVS( 185) = W( 58 )
  JVS( 186) = W( 87 )
  JVS( 187) = W( 90 )
  JVS( 188) = W( 91 )
  IF ( ABS(  JVS( 190 )) < TINY(a) ) THEN
         IER = 38                                      
         RETURN
  END IF
   W( 23 ) = JVS( 189 )
   W( 38 ) = JVS( 190 )
   W( 85 ) = JVS( 191 )
   W( 88 ) = JVS( 192 )
   W( 91 ) = JVS( 193 )
   W( 95 ) = JVS( 194 )
  a = -W( 23 ) / JVS(           94  )
  W( 23 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 95 )
  W( 91 ) = W( 91 ) + a*JVS( 96 )
  JVS( 189) = W( 23 )
  JVS( 190) = W( 38 )
  JVS( 191) = W( 85 )
  JVS( 192) = W( 88 )
  JVS( 193) = W( 91 )
  JVS( 194) = W( 95 )
  IF ( ABS(  JVS( 196 )) < TINY(a) ) THEN
         IER = 39                                      
         RETURN
  END IF
   W( 24 ) = JVS( 195 )
   W( 39 ) = JVS( 196 )
   W( 85 ) = JVS( 197 )
   W( 88 ) = JVS( 198 )
   W( 91 ) = JVS( 199 )
   W( 95 ) = JVS( 200 )
  a = -W( 24 ) / JVS(           97  )
  W( 24 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 98 )
  W( 91 ) = W( 91 ) + a*JVS( 99 )
  JVS( 195) = W( 24 )
  JVS( 196) = W( 39 )
  JVS( 197) = W( 85 )
  JVS( 198) = W( 88 )
  JVS( 199) = W( 91 )
  JVS( 200) = W( 95 )
  IF ( ABS(  JVS( 201 )) < TINY(a) ) THEN
         IER = 40                                      
         RETURN
  END IF
   W( 40 ) = JVS( 201 )
   W( 57 ) = JVS( 202 )
   W( 85 ) = JVS( 203 )
   W( 91 ) = JVS( 204 )
   W( 95 ) = JVS( 205 )
  JVS( 201) = W( 40 )
  JVS( 202) = W( 57 )
  JVS( 203) = W( 85 )
  JVS( 204) = W( 91 )
  JVS( 205) = W( 95 )
  IF ( ABS(  JVS( 206 )) < TINY(a) ) THEN
         IER = 41                                      
         RETURN
  END IF
   W( 41 ) = JVS( 206 )
   W( 81 ) = JVS( 207 )
   W( 82 ) = JVS( 208 )
   W( 88 ) = JVS( 209 )
   W( 91 ) = JVS( 210 )
  JVS( 206) = W( 41 )
  JVS( 207) = W( 81 )
  JVS( 208) = W( 82 )
  JVS( 209) = W( 88 )
  JVS( 210) = W( 91 )
  IF ( ABS(  JVS( 211 )) < TINY(a) ) THEN
         IER = 42                                      
         RETURN
  END IF
   W( 42 ) = JVS( 211 )
   W( 59 ) = JVS( 212 )
   W( 67 ) = JVS( 213 )
   W( 91 ) = JVS( 214 )
   W( 94 ) = JVS( 215 )
   W( 96 ) = JVS( 216 )
  JVS( 211) = W( 42 )
  JVS( 212) = W( 59 )
  JVS( 213) = W( 67 )
  JVS( 214) = W( 91 )
  JVS( 215) = W( 94 )
  JVS( 216) = W( 96 )
  IF ( ABS(  JVS( 217 )) < TINY(a) ) THEN
         IER = 43                                      
         RETURN
  END IF
   W( 43 ) = JVS( 217 )
   W( 85 ) = JVS( 218 )
   W( 91 ) = JVS( 219 )
   W( 96 ) = JVS( 220 )
  JVS( 217) = W( 43 )
  JVS( 218) = W( 85 )
  JVS( 219) = W( 91 )
  JVS( 220) = W( 96 )
  IF ( ABS(  JVS( 221 )) < TINY(a) ) THEN
         IER = 44                                      
         RETURN
  END IF
   W( 44 ) = JVS( 221 )
   W( 57 ) = JVS( 222 )
   W( 60 ) = JVS( 223 )
   W( 69 ) = JVS( 224 )
   W( 73 ) = JVS( 225 )
   W( 78 ) = JVS( 226 )
   W( 79 ) = JVS( 227 )
   W( 89 ) = JVS( 228 )
   W( 91 ) = JVS( 229 )
   W( 95 ) = JVS( 230 )
   W( 96 ) = JVS( 231 )
  JVS( 221) = W( 44 )
  JVS( 222) = W( 57 )
  JVS( 223) = W( 60 )
  JVS( 224) = W( 69 )
  JVS( 225) = W( 73 )
  JVS( 226) = W( 78 )
  JVS( 227) = W( 79 )
  JVS( 228) = W( 89 )
  JVS( 229) = W( 91 )
  JVS( 230) = W( 95 )
  JVS( 231) = W( 96 )
  IF ( ABS(  JVS( 232 )) < TINY(a) ) THEN
         IER = 45                                      
         RETURN
  END IF
   W( 45 ) = JVS( 232 )
   W( 79 ) = JVS( 233 )
   W( 85 ) = JVS( 234 )
   W( 89 ) = JVS( 235 )
   W( 91 ) = JVS( 236 )
   W( 93 ) = JVS( 237 )
  JVS( 232) = W( 45 )
  JVS( 233) = W( 79 )
  JVS( 234) = W( 85 )
  JVS( 235) = W( 89 )
  JVS( 236) = W( 91 )
  JVS( 237) = W( 93 )
  IF ( ABS(  JVS( 238 )) < TINY(a) ) THEN
         IER = 46                                      
         RETURN
  END IF
   W( 46 ) = JVS( 238 )
   W( 81 ) = JVS( 239 )
   W( 82 ) = JVS( 240 )
   W( 88 ) = JVS( 241 )
   W( 91 ) = JVS( 242 )
  JVS( 238) = W( 46 )
  JVS( 239) = W( 81 )
  JVS( 240) = W( 82 )
  JVS( 241) = W( 88 )
  JVS( 242) = W( 91 )
  IF ( ABS(  JVS( 243 )) < TINY(a) ) THEN
         IER = 47                                      
         RETURN
  END IF
   W( 47 ) = JVS( 243 )
   W( 85 ) = JVS( 244 )
   W( 91 ) = JVS( 245 )
   W( 96 ) = JVS( 246 )
  JVS( 243) = W( 47 )
  JVS( 244) = W( 85 )
  JVS( 245) = W( 91 )
  JVS( 246) = W( 96 )
  IF ( ABS(  JVS( 247 )) < TINY(a) ) THEN
         IER = 48                                      
         RETURN
  END IF
   W( 48 ) = JVS( 247 )
   W( 58 ) = JVS( 248 )
   W( 85 ) = JVS( 249 )
   W( 91 ) = JVS( 250 )
   W( 93 ) = JVS( 251 )
   W( 95 ) = JVS( 252 )
   W( 96 ) = JVS( 253 )
  JVS( 247) = W( 48 )
  JVS( 248) = W( 58 )
  JVS( 249) = W( 85 )
  JVS( 250) = W( 91 )
  JVS( 251) = W( 93 )
  JVS( 252) = W( 95 )
  JVS( 253) = W( 96 )
  IF ( ABS(  JVS( 254 )) < TINY(a) ) THEN
         IER = 49                                      
         RETURN
  END IF
   W( 49 ) = JVS( 254 )
   W( 81 ) = JVS( 255 )
   W( 82 ) = JVS( 256 )
   W( 88 ) = JVS( 257 )
   W( 91 ) = JVS( 258 )
  JVS( 254) = W( 49 )
  JVS( 255) = W( 81 )
  JVS( 256) = W( 82 )
  JVS( 257) = W( 88 )
  JVS( 258) = W( 91 )
  IF ( ABS(  JVS( 259 )) < TINY(a) ) THEN
         IER = 50                                      
         RETURN
  END IF
   W( 50 ) = JVS( 259 )
   W( 54 ) = JVS( 260 )
   W( 81 ) = JVS( 261 )
   W( 82 ) = JVS( 262 )
   W( 85 ) = JVS( 263 )
   W( 88 ) = JVS( 264 )
   W( 91 ) = JVS( 265 )
  JVS( 259) = W( 50 )
  JVS( 260) = W( 54 )
  JVS( 261) = W( 81 )
  JVS( 262) = W( 82 )
  JVS( 263) = W( 85 )
  JVS( 264) = W( 88 )
  JVS( 265) = W( 91 )
  IF ( ABS(  JVS( 266 )) < TINY(a) ) THEN
         IER = 51                                      
         RETURN
  END IF
   W( 51 ) = JVS( 266 )
   W( 85 ) = JVS( 267 )
   W( 91 ) = JVS( 268 )
   W( 96 ) = JVS( 269 )
  JVS( 266) = W( 51 )
  JVS( 267) = W( 85 )
  JVS( 268) = W( 91 )
  JVS( 269) = W( 96 )
  IF ( ABS(  JVS( 276 )) < TINY(a) ) THEN
         IER = 52                                      
         RETURN
  END IF
   W( 22 ) = JVS( 270 )
   W( 36 ) = JVS( 271 )
   W( 42 ) = JVS( 272 )
   W( 47 ) = JVS( 273 )
   W( 49 ) = JVS( 274 )
   W( 51 ) = JVS( 275 )
   W( 52 ) = JVS( 276 )
   W( 53 ) = JVS( 277 )
   W( 54 ) = JVS( 278 )
   W( 55 ) = JVS( 279 )
   W( 56 ) = JVS( 280 )
   W( 58 ) = JVS( 281 )
   W( 59 ) = JVS( 282 )
   W( 60 ) = JVS( 283 )
   W( 64 ) = JVS( 284 )
   W( 67 ) = JVS( 285 )
   W( 69 ) = JVS( 286 )
   W( 73 ) = JVS( 287 )
   W( 78 ) = JVS( 288 )
   W( 79 ) = JVS( 289 )
   W( 81 ) = JVS( 290 )
   W( 82 ) = JVS( 291 )
   W( 85 ) = JVS( 292 )
   W( 88 ) = JVS( 293 )
   W( 91 ) = JVS( 294 )
   W( 94 ) = JVS( 295 )
   W( 96 ) = JVS( 296 )
  a = -W( 22 ) / JVS(           91  )
  W( 22 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 92 )
  W( 91 ) = W( 91 ) + a*JVS( 93 )
  a = -W( 36 ) / JVS(          176  )
  W( 36 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 177 )
  W( 59 ) = W( 59 ) + a*JVS( 178 )
  W( 81 ) = W( 81 ) + a*JVS( 179 )
  W( 82 ) = W( 82 ) + a*JVS( 180 )
  W( 88 ) = W( 88 ) + a*JVS( 181 )
  W( 91 ) = W( 91 ) + a*JVS( 182 )
  a = -W( 42 ) / JVS(          211  )
  W( 42 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 212 )
  W( 67 ) = W( 67 ) + a*JVS( 213 )
  W( 91 ) = W( 91 ) + a*JVS( 214 )
  W( 94 ) = W( 94 ) + a*JVS( 215 )
  W( 96 ) = W( 96 ) + a*JVS( 216 )
  a = -W( 47 ) / JVS(          243  )
  W( 47 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 244 )
  W( 91 ) = W( 91 ) + a*JVS( 245 )
  W( 96 ) = W( 96 ) + a*JVS( 246 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 255 )
  W( 82 ) = W( 82 ) + a*JVS( 256 )
  W( 88 ) = W( 88 ) + a*JVS( 257 )
  W( 91 ) = W( 91 ) + a*JVS( 258 )
  a = -W( 51 ) / JVS(          266  )
  W( 51 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 267 )
  W( 91 ) = W( 91 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  JVS( 270) = W( 22 )
  JVS( 271) = W( 36 )
  JVS( 272) = W( 42 )
  JVS( 273) = W( 47 )
  JVS( 274) = W( 49 )
  JVS( 275) = W( 51 )
  JVS( 276) = W( 52 )
  JVS( 277) = W( 53 )
  JVS( 278) = W( 54 )
  JVS( 279) = W( 55 )
  JVS( 280) = W( 56 )
  JVS( 281) = W( 58 )
  JVS( 282) = W( 59 )
  JVS( 283) = W( 60 )
  JVS( 284) = W( 64 )
  JVS( 285) = W( 67 )
  JVS( 286) = W( 69 )
  JVS( 287) = W( 73 )
  JVS( 288) = W( 78 )
  JVS( 289) = W( 79 )
  JVS( 290) = W( 81 )
  JVS( 291) = W( 82 )
  JVS( 292) = W( 85 )
  JVS( 293) = W( 88 )
  JVS( 294) = W( 91 )
  JVS( 295) = W( 94 )
  JVS( 296) = W( 96 )
  IF ( ABS(  JVS( 297 )) < TINY(a) ) THEN
         IER = 53                                      
         RETURN
  END IF
   W( 53 ) = JVS( 297 )
   W( 85 ) = JVS( 298 )
   W( 88 ) = JVS( 299 )
   W( 91 ) = JVS( 300 )
   W( 96 ) = JVS( 301 )
  JVS( 297) = W( 53 )
  JVS( 298) = W( 85 )
  JVS( 299) = W( 88 )
  JVS( 300) = W( 91 )
  JVS( 301) = W( 96 )
  IF ( ABS(  JVS( 302 )) < TINY(a) ) THEN
         IER = 54                                      
         RETURN
  END IF
   W( 54 ) = JVS( 302 )
   W( 85 ) = JVS( 303 )
   W( 91 ) = JVS( 304 )
   W( 96 ) = JVS( 305 )
  JVS( 302) = W( 54 )
  JVS( 303) = W( 85 )
  JVS( 304) = W( 91 )
  JVS( 305) = W( 96 )
  IF ( ABS(  JVS( 306 )) < TINY(a) ) THEN
         IER = 55                                      
         RETURN
  END IF
   W( 55 ) = JVS( 306 )
   W( 85 ) = JVS( 307 )
   W( 91 ) = JVS( 308 )
   W( 96 ) = JVS( 309 )
  JVS( 306) = W( 55 )
  JVS( 307) = W( 85 )
  JVS( 308) = W( 91 )
  JVS( 309) = W( 96 )
  IF ( ABS(  JVS( 312 )) < TINY(a) ) THEN
         IER = 56                                      
         RETURN
  END IF
   W( 20 ) = JVS( 310 )
   W( 26 ) = JVS( 311 )
   W( 56 ) = JVS( 312 )
   W( 59 ) = JVS( 313 )
   W( 64 ) = JVS( 314 )
   W( 85 ) = JVS( 315 )
   W( 89 ) = JVS( 316 )
   W( 91 ) = JVS( 317 )
   W( 94 ) = JVS( 318 )
   W( 95 ) = JVS( 319 )
  a = -W( 20 ) / JVS(           84  )
  W( 20 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 85 )
  W( 89 ) = W( 89 ) + a*JVS( 86 )
  W( 91 ) = W( 91 ) + a*JVS( 87 )
  a = -W( 26 ) / JVS(          104  )
  W( 26 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 105 )
  W( 91 ) = W( 91 ) + a*JVS( 106 )
  W( 95 ) = W( 95 ) + a*JVS( 107 )
  JVS( 310) = W( 20 )
  JVS( 311) = W( 26 )
  JVS( 312) = W( 56 )
  JVS( 313) = W( 59 )
  JVS( 314) = W( 64 )
  JVS( 315) = W( 85 )
  JVS( 316) = W( 89 )
  JVS( 317) = W( 91 )
  JVS( 318) = W( 94 )
  JVS( 319) = W( 95 )
  IF ( ABS(  JVS( 324 )) < TINY(a) ) THEN
         IER = 57                                      
         RETURN
  END IF
   W( 31 ) = JVS( 320 )
   W( 38 ) = JVS( 321 )
   W( 39 ) = JVS( 322 )
   W( 40 ) = JVS( 323 )
   W( 57 ) = JVS( 324 )
   W( 85 ) = JVS( 325 )
   W( 88 ) = JVS( 326 )
   W( 89 ) = JVS( 327 )
   W( 91 ) = JVS( 328 )
   W( 95 ) = JVS( 329 )
   W( 96 ) = JVS( 330 )
  a = -W( 31 ) / JVS(          133  )
  W( 31 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 134 )
  W( 89 ) = W( 89 ) + a*JVS( 135 )
  W( 91 ) = W( 91 ) + a*JVS( 136 )
  W( 95 ) = W( 95 ) + a*JVS( 137 )
  W( 96 ) = W( 96 ) + a*JVS( 138 )
  a = -W( 38 ) / JVS(          190  )
  W( 38 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 191 )
  W( 88 ) = W( 88 ) + a*JVS( 192 )
  W( 91 ) = W( 91 ) + a*JVS( 193 )
  W( 95 ) = W( 95 ) + a*JVS( 194 )
  a = -W( 39 ) / JVS(          196  )
  W( 39 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 197 )
  W( 88 ) = W( 88 ) + a*JVS( 198 )
  W( 91 ) = W( 91 ) + a*JVS( 199 )
  W( 95 ) = W( 95 ) + a*JVS( 200 )
  a = -W( 40 ) / JVS(          201  )
  W( 40 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 202 )
  W( 85 ) = W( 85 ) + a*JVS( 203 )
  W( 91 ) = W( 91 ) + a*JVS( 204 )
  W( 95 ) = W( 95 ) + a*JVS( 205 )
  JVS( 320) = W( 31 )
  JVS( 321) = W( 38 )
  JVS( 322) = W( 39 )
  JVS( 323) = W( 40 )
  JVS( 324) = W( 57 )
  JVS( 325) = W( 85 )
  JVS( 326) = W( 88 )
  JVS( 327) = W( 89 )
  JVS( 328) = W( 91 )
  JVS( 329) = W( 95 )
  JVS( 330) = W( 96 )
  IF ( ABS(  JVS( 331 )) < TINY(a) ) THEN
         IER = 58                                      
         RETURN
  END IF
   W( 58 ) = JVS( 331 )
   W( 74 ) = JVS( 332 )
   W( 85 ) = JVS( 333 )
   W( 91 ) = JVS( 334 )
   W( 95 ) = JVS( 335 )
   W( 96 ) = JVS( 336 )
  JVS( 331) = W( 58 )
  JVS( 332) = W( 74 )
  JVS( 333) = W( 85 )
  JVS( 334) = W( 91 )
  JVS( 335) = W( 95 )
  JVS( 336) = W( 96 )
  IF ( ABS(  JVS( 337 )) < TINY(a) ) THEN
         IER = 59                                      
         RETURN
  END IF
   W( 59 ) = JVS( 337 )
   W( 85 ) = JVS( 338 )
   W( 88 ) = JVS( 339 )
   W( 91 ) = JVS( 340 )
   W( 96 ) = JVS( 341 )
  JVS( 337) = W( 59 )
  JVS( 338) = W( 85 )
  JVS( 339) = W( 88 )
  JVS( 340) = W( 91 )
  JVS( 341) = W( 96 )
  IF ( ABS(  JVS( 343 )) < TINY(a) ) THEN
         IER = 60                                      
         RETURN
  END IF
   W( 49 ) = JVS( 342 )
   W( 60 ) = JVS( 343 )
   W( 68 ) = JVS( 344 )
   W( 71 ) = JVS( 345 )
   W( 72 ) = JVS( 346 )
   W( 79 ) = JVS( 347 )
   W( 81 ) = JVS( 348 )
   W( 82 ) = JVS( 349 )
   W( 85 ) = JVS( 350 )
   W( 86 ) = JVS( 351 )
   W( 88 ) = JVS( 352 )
   W( 90 ) = JVS( 353 )
   W( 91 ) = JVS( 354 )
   W( 93 ) = JVS( 355 )
   W( 94 ) = JVS( 356 )
   W( 96 ) = JVS( 357 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 255 )
  W( 82 ) = W( 82 ) + a*JVS( 256 )
  W( 88 ) = W( 88 ) + a*JVS( 257 )
  W( 91 ) = W( 91 ) + a*JVS( 258 )
  JVS( 342) = W( 49 )
  JVS( 343) = W( 60 )
  JVS( 344) = W( 68 )
  JVS( 345) = W( 71 )
  JVS( 346) = W( 72 )
  JVS( 347) = W( 79 )
  JVS( 348) = W( 81 )
  JVS( 349) = W( 82 )
  JVS( 350) = W( 85 )
  JVS( 351) = W( 86 )
  JVS( 352) = W( 88 )
  JVS( 353) = W( 90 )
  JVS( 354) = W( 91 )
  JVS( 355) = W( 93 )
  JVS( 356) = W( 94 )
  JVS( 357) = W( 96 )
  IF ( ABS(  JVS( 358 )) < TINY(a) ) THEN
         IER = 61                                      
         RETURN
  END IF
   W( 61 ) = JVS( 358 )
   W( 82 ) = JVS( 359 )
   W( 89 ) = JVS( 360 )
   W( 90 ) = JVS( 361 )
   W( 91 ) = JVS( 362 )
   W( 93 ) = JVS( 363 )
   W( 94 ) = JVS( 364 )
   W( 96 ) = JVS( 365 )
  JVS( 358) = W( 61 )
  JVS( 359) = W( 82 )
  JVS( 360) = W( 89 )
  JVS( 361) = W( 90 )
  JVS( 362) = W( 91 )
  JVS( 363) = W( 93 )
  JVS( 364) = W( 94 )
  JVS( 365) = W( 96 )
  IF ( ABS(  JVS( 367 )) < TINY(a) ) THEN
         IER = 62                                      
         RETURN
  END IF
   W( 53 ) = JVS( 366 )
   W( 62 ) = JVS( 367 )
   W( 85 ) = JVS( 368 )
   W( 88 ) = JVS( 369 )
   W( 89 ) = JVS( 370 )
   W( 90 ) = JVS( 371 )
   W( 91 ) = JVS( 372 )
   W( 93 ) = JVS( 373 )
   W( 94 ) = JVS( 374 )
   W( 96 ) = JVS( 375 )
  a = -W( 53 ) / JVS(          297  )
  W( 53 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 298 )
  W( 88 ) = W( 88 ) + a*JVS( 299 )
  W( 91 ) = W( 91 ) + a*JVS( 300 )
  W( 96 ) = W( 96 ) + a*JVS( 301 )
  JVS( 366) = W( 53 )
  JVS( 367) = W( 62 )
  JVS( 368) = W( 85 )
  JVS( 369) = W( 88 )
  JVS( 370) = W( 89 )
  JVS( 371) = W( 90 )
  JVS( 372) = W( 91 )
  JVS( 373) = W( 93 )
  JVS( 374) = W( 94 )
  JVS( 375) = W( 96 )
  IF ( ABS(  JVS( 377 )) < TINY(a) ) THEN
         IER = 63                                      
         RETURN
  END IF
   W( 51 ) = JVS( 376 )
   W( 63 ) = JVS( 377 )
   W( 85 ) = JVS( 378 )
   W( 89 ) = JVS( 379 )
   W( 90 ) = JVS( 380 )
   W( 91 ) = JVS( 381 )
   W( 93 ) = JVS( 382 )
   W( 94 ) = JVS( 383 )
   W( 96 ) = JVS( 384 )
  a = -W( 51 ) / JVS(          266  )
  W( 51 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 267 )
  W( 91 ) = W( 91 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  JVS( 376) = W( 51 )
  JVS( 377) = W( 63 )
  JVS( 378) = W( 85 )
  JVS( 379) = W( 89 )
  JVS( 380) = W( 90 )
  JVS( 381) = W( 91 )
  JVS( 382) = W( 93 )
  JVS( 383) = W( 94 )
  JVS( 384) = W( 96 )
  IF ( ABS(  JVS( 391 )) < TINY(a) ) THEN
         IER = 64                                      
         RETURN
  END IF
   W( 16 ) = JVS( 385 )
   W( 43 ) = JVS( 386 )
   W( 51 ) = JVS( 387 )
   W( 55 ) = JVS( 388 )
   W( 59 ) = JVS( 389 )
   W( 63 ) = JVS( 390 )
   W( 64 ) = JVS( 391 )
   W( 67 ) = JVS( 392 )
   W( 85 ) = JVS( 393 )
   W( 88 ) = JVS( 394 )
   W( 89 ) = JVS( 395 )
   W( 90 ) = JVS( 396 )
   W( 91 ) = JVS( 397 )
   W( 93 ) = JVS( 398 )
   W( 94 ) = JVS( 399 )
   W( 96 ) = JVS( 400 )
  a = -W( 16 ) / JVS(           71  )
  W( 16 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 72 )
  W( 89 ) = W( 89 ) + a*JVS( 73 )
  W( 91 ) = W( 91 ) + a*JVS( 74 )
  a = -W( 43 ) / JVS(          217  )
  W( 43 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 218 )
  W( 91 ) = W( 91 ) + a*JVS( 219 )
  W( 96 ) = W( 96 ) + a*JVS( 220 )
  a = -W( 51 ) / JVS(          266  )
  W( 51 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 267 )
  W( 91 ) = W( 91 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 59 ) / JVS(          337  )
  W( 59 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 338 )
  W( 88 ) = W( 88 ) + a*JVS( 339 )
  W( 91 ) = W( 91 ) + a*JVS( 340 )
  W( 96 ) = W( 96 ) + a*JVS( 341 )
  a = -W( 63 ) / JVS(          377  )
  W( 63 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  W( 90 ) = W( 90 ) + a*JVS( 380 )
  W( 91 ) = W( 91 ) + a*JVS( 381 )
  W( 93 ) = W( 93 ) + a*JVS( 382 )
  W( 94 ) = W( 94 ) + a*JVS( 383 )
  W( 96 ) = W( 96 ) + a*JVS( 384 )
  JVS( 385) = W( 16 )
  JVS( 386) = W( 43 )
  JVS( 387) = W( 51 )
  JVS( 388) = W( 55 )
  JVS( 389) = W( 59 )
  JVS( 390) = W( 63 )
  JVS( 391) = W( 64 )
  JVS( 392) = W( 67 )
  JVS( 393) = W( 85 )
  JVS( 394) = W( 88 )
  JVS( 395) = W( 89 )
  JVS( 396) = W( 90 )
  JVS( 397) = W( 91 )
  JVS( 398) = W( 93 )
  JVS( 399) = W( 94 )
  JVS( 400) = W( 96 )
  IF ( ABS(  JVS( 402 )) < TINY(a) ) THEN
         IER = 65                                      
         RETURN
  END IF
   W( 47 ) = JVS( 401 )
   W( 65 ) = JVS( 402 )
   W( 85 ) = JVS( 403 )
   W( 89 ) = JVS( 404 )
   W( 90 ) = JVS( 405 )
   W( 91 ) = JVS( 406 )
   W( 93 ) = JVS( 407 )
   W( 94 ) = JVS( 408 )
   W( 96 ) = JVS( 409 )
  a = -W( 47 ) / JVS(          243  )
  W( 47 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 244 )
  W( 91 ) = W( 91 ) + a*JVS( 245 )
  W( 96 ) = W( 96 ) + a*JVS( 246 )
  JVS( 401) = W( 47 )
  JVS( 402) = W( 65 )
  JVS( 403) = W( 85 )
  JVS( 404) = W( 89 )
  JVS( 405) = W( 90 )
  JVS( 406) = W( 91 )
  JVS( 407) = W( 93 )
  JVS( 408) = W( 94 )
  JVS( 409) = W( 96 )
  IF ( ABS(  JVS( 417 )) < TINY(a) ) THEN
         IER = 66                                      
         RETURN
  END IF
   W( 21 ) = JVS( 410 )
   W( 41 ) = JVS( 411 )
   W( 43 ) = JVS( 412 )
   W( 47 ) = JVS( 413 )
   W( 54 ) = JVS( 414 )
   W( 61 ) = JVS( 415 )
   W( 65 ) = JVS( 416 )
   W( 66 ) = JVS( 417 )
   W( 70 ) = JVS( 418 )
   W( 75 ) = JVS( 419 )
   W( 76 ) = JVS( 420 )
   W( 77 ) = JVS( 421 )
   W( 79 ) = JVS( 422 )
   W( 81 ) = JVS( 423 )
   W( 82 ) = JVS( 424 )
   W( 83 ) = JVS( 425 )
   W( 84 ) = JVS( 426 )
   W( 85 ) = JVS( 427 )
   W( 86 ) = JVS( 428 )
   W( 87 ) = JVS( 429 )
   W( 88 ) = JVS( 430 )
   W( 89 ) = JVS( 431 )
   W( 90 ) = JVS( 432 )
   W( 91 ) = JVS( 433 )
   W( 92 ) = JVS( 434 )
   W( 93 ) = JVS( 435 )
   W( 94 ) = JVS( 436 )
   W( 96 ) = JVS( 437 )
  a = -W( 21 ) / JVS(           88  )
  W( 21 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 89 )
  W( 91 ) = W( 91 ) + a*JVS( 90 )
  a = -W( 41 ) / JVS(          206  )
  W( 41 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 207 )
  W( 82 ) = W( 82 ) + a*JVS( 208 )
  W( 88 ) = W( 88 ) + a*JVS( 209 )
  W( 91 ) = W( 91 ) + a*JVS( 210 )
  a = -W( 43 ) / JVS(          217  )
  W( 43 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 218 )
  W( 91 ) = W( 91 ) + a*JVS( 219 )
  W( 96 ) = W( 96 ) + a*JVS( 220 )
  a = -W( 47 ) / JVS(          243  )
  W( 47 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 244 )
  W( 91 ) = W( 91 ) + a*JVS( 245 )
  W( 96 ) = W( 96 ) + a*JVS( 246 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  a = -W( 61 ) / JVS(          358  )
  W( 61 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 359 )
  W( 89 ) = W( 89 ) + a*JVS( 360 )
  W( 90 ) = W( 90 ) + a*JVS( 361 )
  W( 91 ) = W( 91 ) + a*JVS( 362 )
  W( 93 ) = W( 93 ) + a*JVS( 363 )
  W( 94 ) = W( 94 ) + a*JVS( 364 )
  W( 96 ) = W( 96 ) + a*JVS( 365 )
  a = -W( 65 ) / JVS(          402  )
  W( 65 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 403 )
  W( 89 ) = W( 89 ) + a*JVS( 404 )
  W( 90 ) = W( 90 ) + a*JVS( 405 )
  W( 91 ) = W( 91 ) + a*JVS( 406 )
  W( 93 ) = W( 93 ) + a*JVS( 407 )
  W( 94 ) = W( 94 ) + a*JVS( 408 )
  W( 96 ) = W( 96 ) + a*JVS( 409 )
  JVS( 410) = W( 21 )
  JVS( 411) = W( 41 )
  JVS( 412) = W( 43 )
  JVS( 413) = W( 47 )
  JVS( 414) = W( 54 )
  JVS( 415) = W( 61 )
  JVS( 416) = W( 65 )
  JVS( 417) = W( 66 )
  JVS( 418) = W( 70 )
  JVS( 419) = W( 75 )
  JVS( 420) = W( 76 )
  JVS( 421) = W( 77 )
  JVS( 422) = W( 79 )
  JVS( 423) = W( 81 )
  JVS( 424) = W( 82 )
  JVS( 425) = W( 83 )
  JVS( 426) = W( 84 )
  JVS( 427) = W( 85 )
  JVS( 428) = W( 86 )
  JVS( 429) = W( 87 )
  JVS( 430) = W( 88 )
  JVS( 431) = W( 89 )
  JVS( 432) = W( 90 )
  JVS( 433) = W( 91 )
  JVS( 434) = W( 92 )
  JVS( 435) = W( 93 )
  JVS( 436) = W( 94 )
  JVS( 437) = W( 96 )
  IF ( ABS(  JVS( 440 )) < TINY(a) ) THEN
         IER = 67                                      
         RETURN
  END IF
   W( 55 ) = JVS( 438 )
   W( 59 ) = JVS( 439 )
   W( 67 ) = JVS( 440 )
   W( 85 ) = JVS( 441 )
   W( 88 ) = JVS( 442 )
   W( 89 ) = JVS( 443 )
   W( 90 ) = JVS( 444 )
   W( 91 ) = JVS( 445 )
   W( 93 ) = JVS( 446 )
   W( 94 ) = JVS( 447 )
   W( 96 ) = JVS( 448 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 59 ) / JVS(          337  )
  W( 59 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 338 )
  W( 88 ) = W( 88 ) + a*JVS( 339 )
  W( 91 ) = W( 91 ) + a*JVS( 340 )
  W( 96 ) = W( 96 ) + a*JVS( 341 )
  JVS( 438) = W( 55 )
  JVS( 439) = W( 59 )
  JVS( 440) = W( 67 )
  JVS( 441) = W( 85 )
  JVS( 442) = W( 88 )
  JVS( 443) = W( 89 )
  JVS( 444) = W( 90 )
  JVS( 445) = W( 91 )
  JVS( 446) = W( 93 )
  JVS( 447) = W( 94 )
  JVS( 448) = W( 96 )
  IF ( ABS(  JVS( 451 )) < TINY(a) ) THEN
         IER = 68                                      
         RETURN
  END IF
   W( 40 ) = JVS( 449 )
   W( 57 ) = JVS( 450 )
   W( 68 ) = JVS( 451 )
   W( 85 ) = JVS( 452 )
   W( 88 ) = JVS( 453 )
   W( 89 ) = JVS( 454 )
   W( 90 ) = JVS( 455 )
   W( 91 ) = JVS( 456 )
   W( 93 ) = JVS( 457 )
   W( 94 ) = JVS( 458 )
   W( 95 ) = JVS( 459 )
   W( 96 ) = JVS( 460 )
  a = -W( 40 ) / JVS(          201  )
  W( 40 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 202 )
  W( 85 ) = W( 85 ) + a*JVS( 203 )
  W( 91 ) = W( 91 ) + a*JVS( 204 )
  W( 95 ) = W( 95 ) + a*JVS( 205 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 325 )
  W( 88 ) = W( 88 ) + a*JVS( 326 )
  W( 89 ) = W( 89 ) + a*JVS( 327 )
  W( 91 ) = W( 91 ) + a*JVS( 328 )
  W( 95 ) = W( 95 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  JVS( 449) = W( 40 )
  JVS( 450) = W( 57 )
  JVS( 451) = W( 68 )
  JVS( 452) = W( 85 )
  JVS( 453) = W( 88 )
  JVS( 454) = W( 89 )
  JVS( 455) = W( 90 )
  JVS( 456) = W( 91 )
  JVS( 457) = W( 93 )
  JVS( 458) = W( 94 )
  JVS( 459) = W( 95 )
  JVS( 460) = W( 96 )
  IF ( ABS(  JVS( 471 )) < TINY(a) ) THEN
         IER = 69                                      
         RETURN
  END IF
   W( 28 ) = JVS( 461 )
   W( 37 ) = JVS( 462 )
   W( 42 ) = JVS( 463 )
   W( 46 ) = JVS( 464 )
   W( 56 ) = JVS( 465 )
   W( 58 ) = JVS( 466 )
   W( 59 ) = JVS( 467 )
   W( 64 ) = JVS( 468 )
   W( 67 ) = JVS( 469 )
   W( 68 ) = JVS( 470 )
   W( 69 ) = JVS( 471 )
   W( 71 ) = JVS( 472 )
   W( 72 ) = JVS( 473 )
   W( 74 ) = JVS( 474 )
   W( 79 ) = JVS( 475 )
   W( 81 ) = JVS( 476 )
   W( 82 ) = JVS( 477 )
   W( 85 ) = JVS( 478 )
   W( 86 ) = JVS( 479 )
   W( 87 ) = JVS( 480 )
   W( 88 ) = JVS( 481 )
   W( 89 ) = JVS( 482 )
   W( 90 ) = JVS( 483 )
   W( 91 ) = JVS( 484 )
   W( 93 ) = JVS( 485 )
   W( 94 ) = JVS( 486 )
   W( 95 ) = JVS( 487 )
   W( 96 ) = JVS( 488 )
  a = -W( 28 ) / JVS(          113  )
  W( 28 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 114 )
  W( 56 ) = W( 56 ) + a*JVS( 115 )
  W( 91 ) = W( 91 ) + a*JVS( 116 )
  W( 94 ) = W( 94 ) + a*JVS( 117 )
  W( 95 ) = W( 95 ) + a*JVS( 118 )
  a = -W( 37 ) / JVS(          183  )
  W( 37 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 184 )
  W( 58 ) = W( 58 ) + a*JVS( 185 )
  W( 87 ) = W( 87 ) + a*JVS( 186 )
  W( 90 ) = W( 90 ) + a*JVS( 187 )
  W( 91 ) = W( 91 ) + a*JVS( 188 )
  a = -W( 42 ) / JVS(          211  )
  W( 42 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 212 )
  W( 67 ) = W( 67 ) + a*JVS( 213 )
  W( 91 ) = W( 91 ) + a*JVS( 214 )
  W( 94 ) = W( 94 ) + a*JVS( 215 )
  W( 96 ) = W( 96 ) + a*JVS( 216 )
  a = -W( 46 ) / JVS(          238  )
  W( 46 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 239 )
  W( 82 ) = W( 82 ) + a*JVS( 240 )
  W( 88 ) = W( 88 ) + a*JVS( 241 )
  W( 91 ) = W( 91 ) + a*JVS( 242 )
  a = -W( 56 ) / JVS(          312  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 313 )
  W( 64 ) = W( 64 ) + a*JVS( 314 )
  W( 85 ) = W( 85 ) + a*JVS( 315 )
  W( 89 ) = W( 89 ) + a*JVS( 316 )
  W( 91 ) = W( 91 ) + a*JVS( 317 )
  W( 94 ) = W( 94 ) + a*JVS( 318 )
  W( 95 ) = W( 95 ) + a*JVS( 319 )
  a = -W( 58 ) / JVS(          331  )
  W( 58 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 332 )
  W( 85 ) = W( 85 ) + a*JVS( 333 )
  W( 91 ) = W( 91 ) + a*JVS( 334 )
  W( 95 ) = W( 95 ) + a*JVS( 335 )
  W( 96 ) = W( 96 ) + a*JVS( 336 )
  a = -W( 59 ) / JVS(          337  )
  W( 59 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 338 )
  W( 88 ) = W( 88 ) + a*JVS( 339 )
  W( 91 ) = W( 91 ) + a*JVS( 340 )
  W( 96 ) = W( 96 ) + a*JVS( 341 )
  a = -W( 64 ) / JVS(          391  )
  W( 64 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 392 )
  W( 85 ) = W( 85 ) + a*JVS( 393 )
  W( 88 ) = W( 88 ) + a*JVS( 394 )
  W( 89 ) = W( 89 ) + a*JVS( 395 )
  W( 90 ) = W( 90 ) + a*JVS( 396 )
  W( 91 ) = W( 91 ) + a*JVS( 397 )
  W( 93 ) = W( 93 ) + a*JVS( 398 )
  W( 94 ) = W( 94 ) + a*JVS( 399 )
  W( 96 ) = W( 96 ) + a*JVS( 400 )
  a = -W( 67 ) / JVS(          440  )
  W( 67 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 441 )
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 89 ) = W( 89 ) + a*JVS( 443 )
  W( 90 ) = W( 90 ) + a*JVS( 444 )
  W( 91 ) = W( 91 ) + a*JVS( 445 )
  W( 93 ) = W( 93 ) + a*JVS( 446 )
  W( 94 ) = W( 94 ) + a*JVS( 447 )
  W( 96 ) = W( 96 ) + a*JVS( 448 )
  a = -W( 68 ) / JVS(          451  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 452 )
  W( 88 ) = W( 88 ) + a*JVS( 453 )
  W( 89 ) = W( 89 ) + a*JVS( 454 )
  W( 90 ) = W( 90 ) + a*JVS( 455 )
  W( 91 ) = W( 91 ) + a*JVS( 456 )
  W( 93 ) = W( 93 ) + a*JVS( 457 )
  W( 94 ) = W( 94 ) + a*JVS( 458 )
  W( 95 ) = W( 95 ) + a*JVS( 459 )
  W( 96 ) = W( 96 ) + a*JVS( 460 )
  JVS( 461) = W( 28 )
  JVS( 462) = W( 37 )
  JVS( 463) = W( 42 )
  JVS( 464) = W( 46 )
  JVS( 465) = W( 56 )
  JVS( 466) = W( 58 )
  JVS( 467) = W( 59 )
  JVS( 468) = W( 64 )
  JVS( 469) = W( 67 )
  JVS( 470) = W( 68 )
  JVS( 471) = W( 69 )
  JVS( 472) = W( 71 )
  JVS( 473) = W( 72 )
  JVS( 474) = W( 74 )
  JVS( 475) = W( 79 )
  JVS( 476) = W( 81 )
  JVS( 477) = W( 82 )
  JVS( 478) = W( 85 )
  JVS( 479) = W( 86 )
  JVS( 480) = W( 87 )
  JVS( 481) = W( 88 )
  JVS( 482) = W( 89 )
  JVS( 483) = W( 90 )
  JVS( 484) = W( 91 )
  JVS( 485) = W( 93 )
  JVS( 486) = W( 94 )
  JVS( 487) = W( 95 )
  JVS( 488) = W( 96 )
  IF ( ABS(  JVS( 491 )) < TINY(a) ) THEN
         IER = 70                                      
         RETURN
  END IF
   W( 43 ) = JVS( 489 )
   W( 54 ) = JVS( 490 )
   W( 70 ) = JVS( 491 )
   W( 81 ) = JVS( 492 )
   W( 85 ) = JVS( 493 )
   W( 89 ) = JVS( 494 )
   W( 90 ) = JVS( 495 )
   W( 91 ) = JVS( 496 )
   W( 93 ) = JVS( 497 )
   W( 94 ) = JVS( 498 )
   W( 96 ) = JVS( 499 )
  a = -W( 43 ) / JVS(          217  )
  W( 43 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 218 )
  W( 91 ) = W( 91 ) + a*JVS( 219 )
  W( 96 ) = W( 96 ) + a*JVS( 220 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  JVS( 489) = W( 43 )
  JVS( 490) = W( 54 )
  JVS( 491) = W( 70 )
  JVS( 492) = W( 81 )
  JVS( 493) = W( 85 )
  JVS( 494) = W( 89 )
  JVS( 495) = W( 90 )
  JVS( 496) = W( 91 )
  JVS( 497) = W( 93 )
  JVS( 498) = W( 94 )
  JVS( 499) = W( 96 )
  IF ( ABS(  JVS( 501 )) < TINY(a) ) THEN
         IER = 71                                      
         RETURN
  END IF
   W( 38 ) = JVS( 500 )
   W( 71 ) = JVS( 501 )
   W( 85 ) = JVS( 502 )
   W( 88 ) = JVS( 503 )
   W( 89 ) = JVS( 504 )
   W( 90 ) = JVS( 505 )
   W( 91 ) = JVS( 506 )
   W( 93 ) = JVS( 507 )
   W( 94 ) = JVS( 508 )
   W( 95 ) = JVS( 509 )
   W( 96 ) = JVS( 510 )
  a = -W( 38 ) / JVS(          190  )
  W( 38 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 191 )
  W( 88 ) = W( 88 ) + a*JVS( 192 )
  W( 91 ) = W( 91 ) + a*JVS( 193 )
  W( 95 ) = W( 95 ) + a*JVS( 194 )
  JVS( 500) = W( 38 )
  JVS( 501) = W( 71 )
  JVS( 502) = W( 85 )
  JVS( 503) = W( 88 )
  JVS( 504) = W( 89 )
  JVS( 505) = W( 90 )
  JVS( 506) = W( 91 )
  JVS( 507) = W( 93 )
  JVS( 508) = W( 94 )
  JVS( 509) = W( 95 )
  JVS( 510) = W( 96 )
  IF ( ABS(  JVS( 512 )) < TINY(a) ) THEN
         IER = 72                                      
         RETURN
  END IF
   W( 39 ) = JVS( 511 )
   W( 72 ) = JVS( 512 )
   W( 85 ) = JVS( 513 )
   W( 88 ) = JVS( 514 )
   W( 89 ) = JVS( 515 )
   W( 90 ) = JVS( 516 )
   W( 91 ) = JVS( 517 )
   W( 93 ) = JVS( 518 )
   W( 94 ) = JVS( 519 )
   W( 95 ) = JVS( 520 )
   W( 96 ) = JVS( 521 )
  a = -W( 39 ) / JVS(          196  )
  W( 39 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 197 )
  W( 88 ) = W( 88 ) + a*JVS( 198 )
  W( 91 ) = W( 91 ) + a*JVS( 199 )
  W( 95 ) = W( 95 ) + a*JVS( 200 )
  JVS( 511) = W( 39 )
  JVS( 512) = W( 72 )
  JVS( 513) = W( 85 )
  JVS( 514) = W( 88 )
  JVS( 515) = W( 89 )
  JVS( 516) = W( 90 )
  JVS( 517) = W( 91 )
  JVS( 518) = W( 93 )
  JVS( 519) = W( 94 )
  JVS( 520) = W( 95 )
  JVS( 521) = W( 96 )
  IF ( ABS(  JVS( 549 )) < TINY(a) ) THEN
         IER = 73                                      
         RETURN
  END IF
   W( 22 ) = JVS( 522 )
   W( 25 ) = JVS( 523 )
   W( 37 ) = JVS( 524 )
   W( 42 ) = JVS( 525 )
   W( 43 ) = JVS( 526 )
   W( 45 ) = JVS( 527 )
   W( 46 ) = JVS( 528 )
   W( 48 ) = JVS( 529 )
   W( 49 ) = JVS( 530 )
   W( 51 ) = JVS( 531 )
   W( 53 ) = JVS( 532 )
   W( 54 ) = JVS( 533 )
   W( 55 ) = JVS( 534 )
   W( 56 ) = JVS( 535 )
   W( 58 ) = JVS( 536 )
   W( 59 ) = JVS( 537 )
   W( 60 ) = JVS( 538 )
   W( 61 ) = JVS( 539 )
   W( 62 ) = JVS( 540 )
   W( 63 ) = JVS( 541 )
   W( 64 ) = JVS( 542 )
   W( 65 ) = JVS( 543 )
   W( 67 ) = JVS( 544 )
   W( 68 ) = JVS( 545 )
   W( 70 ) = JVS( 546 )
   W( 71 ) = JVS( 547 )
   W( 72 ) = JVS( 548 )
   W( 73 ) = JVS( 549 )
   W( 74 ) = JVS( 550 )
   W( 76 ) = JVS( 551 )
   W( 77 ) = JVS( 552 )
   W( 79 ) = JVS( 553 )
   W( 80 ) = JVS( 554 )
   W( 81 ) = JVS( 555 )
   W( 82 ) = JVS( 556 )
   W( 83 ) = JVS( 557 )
   W( 84 ) = JVS( 558 )
   W( 85 ) = JVS( 559 )
   W( 86 ) = JVS( 560 )
   W( 87 ) = JVS( 561 )
   W( 88 ) = JVS( 562 )
   W( 89 ) = JVS( 563 )
   W( 90 ) = JVS( 564 )
   W( 91 ) = JVS( 565 )
   W( 93 ) = JVS( 566 )
   W( 94 ) = JVS( 567 )
   W( 95 ) = JVS( 568 )
   W( 96 ) = JVS( 569 )
   W( 97 ) = JVS( 570 )
  a = -W( 22 ) / JVS(           91  )
  W( 22 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 92 )
  W( 91 ) = W( 91 ) + a*JVS( 93 )
  a = -W( 25 ) / JVS(          100  )
  W( 25 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 101 )
  W( 90 ) = W( 90 ) + a*JVS( 102 )
  W( 91 ) = W( 91 ) + a*JVS( 103 )
  a = -W( 37 ) / JVS(          183  )
  W( 37 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 184 )
  W( 58 ) = W( 58 ) + a*JVS( 185 )
  W( 87 ) = W( 87 ) + a*JVS( 186 )
  W( 90 ) = W( 90 ) + a*JVS( 187 )
  W( 91 ) = W( 91 ) + a*JVS( 188 )
  a = -W( 42 ) / JVS(          211  )
  W( 42 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 212 )
  W( 67 ) = W( 67 ) + a*JVS( 213 )
  W( 91 ) = W( 91 ) + a*JVS( 214 )
  W( 94 ) = W( 94 ) + a*JVS( 215 )
  W( 96 ) = W( 96 ) + a*JVS( 216 )
  a = -W( 43 ) / JVS(          217  )
  W( 43 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 218 )
  W( 91 ) = W( 91 ) + a*JVS( 219 )
  W( 96 ) = W( 96 ) + a*JVS( 220 )
  a = -W( 45 ) / JVS(          232  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 233 )
  W( 85 ) = W( 85 ) + a*JVS( 234 )
  W( 89 ) = W( 89 ) + a*JVS( 235 )
  W( 91 ) = W( 91 ) + a*JVS( 236 )
  W( 93 ) = W( 93 ) + a*JVS( 237 )
  a = -W( 46 ) / JVS(          238  )
  W( 46 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 239 )
  W( 82 ) = W( 82 ) + a*JVS( 240 )
  W( 88 ) = W( 88 ) + a*JVS( 241 )
  W( 91 ) = W( 91 ) + a*JVS( 242 )
  a = -W( 48 ) / JVS(          247  )
  W( 48 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 248 )
  W( 85 ) = W( 85 ) + a*JVS( 249 )
  W( 91 ) = W( 91 ) + a*JVS( 250 )
  W( 93 ) = W( 93 ) + a*JVS( 251 )
  W( 95 ) = W( 95 ) + a*JVS( 252 )
  W( 96 ) = W( 96 ) + a*JVS( 253 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 255 )
  W( 82 ) = W( 82 ) + a*JVS( 256 )
  W( 88 ) = W( 88 ) + a*JVS( 257 )
  W( 91 ) = W( 91 ) + a*JVS( 258 )
  a = -W( 51 ) / JVS(          266  )
  W( 51 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 267 )
  W( 91 ) = W( 91 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 53 ) / JVS(          297  )
  W( 53 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 298 )
  W( 88 ) = W( 88 ) + a*JVS( 299 )
  W( 91 ) = W( 91 ) + a*JVS( 300 )
  W( 96 ) = W( 96 ) + a*JVS( 301 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 56 ) / JVS(          312  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 313 )
  W( 64 ) = W( 64 ) + a*JVS( 314 )
  W( 85 ) = W( 85 ) + a*JVS( 315 )
  W( 89 ) = W( 89 ) + a*JVS( 316 )
  W( 91 ) = W( 91 ) + a*JVS( 317 )
  W( 94 ) = W( 94 ) + a*JVS( 318 )
  W( 95 ) = W( 95 ) + a*JVS( 319 )
  a = -W( 58 ) / JVS(          331  )
  W( 58 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 332 )
  W( 85 ) = W( 85 ) + a*JVS( 333 )
  W( 91 ) = W( 91 ) + a*JVS( 334 )
  W( 95 ) = W( 95 ) + a*JVS( 335 )
  W( 96 ) = W( 96 ) + a*JVS( 336 )
  a = -W( 59 ) / JVS(          337  )
  W( 59 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 338 )
  W( 88 ) = W( 88 ) + a*JVS( 339 )
  W( 91 ) = W( 91 ) + a*JVS( 340 )
  W( 96 ) = W( 96 ) + a*JVS( 341 )
  a = -W( 60 ) / JVS(          343  )
  W( 60 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 344 )
  W( 71 ) = W( 71 ) + a*JVS( 345 )
  W( 72 ) = W( 72 ) + a*JVS( 346 )
  W( 79 ) = W( 79 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 82 ) = W( 82 ) + a*JVS( 349 )
  W( 85 ) = W( 85 ) + a*JVS( 350 )
  W( 86 ) = W( 86 ) + a*JVS( 351 )
  W( 88 ) = W( 88 ) + a*JVS( 352 )
  W( 90 ) = W( 90 ) + a*JVS( 353 )
  W( 91 ) = W( 91 ) + a*JVS( 354 )
  W( 93 ) = W( 93 ) + a*JVS( 355 )
  W( 94 ) = W( 94 ) + a*JVS( 356 )
  W( 96 ) = W( 96 ) + a*JVS( 357 )
  a = -W( 61 ) / JVS(          358  )
  W( 61 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 359 )
  W( 89 ) = W( 89 ) + a*JVS( 360 )
  W( 90 ) = W( 90 ) + a*JVS( 361 )
  W( 91 ) = W( 91 ) + a*JVS( 362 )
  W( 93 ) = W( 93 ) + a*JVS( 363 )
  W( 94 ) = W( 94 ) + a*JVS( 364 )
  W( 96 ) = W( 96 ) + a*JVS( 365 )
  a = -W( 62 ) / JVS(          367  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 368 )
  W( 88 ) = W( 88 ) + a*JVS( 369 )
  W( 89 ) = W( 89 ) + a*JVS( 370 )
  W( 90 ) = W( 90 ) + a*JVS( 371 )
  W( 91 ) = W( 91 ) + a*JVS( 372 )
  W( 93 ) = W( 93 ) + a*JVS( 373 )
  W( 94 ) = W( 94 ) + a*JVS( 374 )
  W( 96 ) = W( 96 ) + a*JVS( 375 )
  a = -W( 63 ) / JVS(          377  )
  W( 63 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  W( 90 ) = W( 90 ) + a*JVS( 380 )
  W( 91 ) = W( 91 ) + a*JVS( 381 )
  W( 93 ) = W( 93 ) + a*JVS( 382 )
  W( 94 ) = W( 94 ) + a*JVS( 383 )
  W( 96 ) = W( 96 ) + a*JVS( 384 )
  a = -W( 64 ) / JVS(          391  )
  W( 64 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 392 )
  W( 85 ) = W( 85 ) + a*JVS( 393 )
  W( 88 ) = W( 88 ) + a*JVS( 394 )
  W( 89 ) = W( 89 ) + a*JVS( 395 )
  W( 90 ) = W( 90 ) + a*JVS( 396 )
  W( 91 ) = W( 91 ) + a*JVS( 397 )
  W( 93 ) = W( 93 ) + a*JVS( 398 )
  W( 94 ) = W( 94 ) + a*JVS( 399 )
  W( 96 ) = W( 96 ) + a*JVS( 400 )
  a = -W( 65 ) / JVS(          402  )
  W( 65 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 403 )
  W( 89 ) = W( 89 ) + a*JVS( 404 )
  W( 90 ) = W( 90 ) + a*JVS( 405 )
  W( 91 ) = W( 91 ) + a*JVS( 406 )
  W( 93 ) = W( 93 ) + a*JVS( 407 )
  W( 94 ) = W( 94 ) + a*JVS( 408 )
  W( 96 ) = W( 96 ) + a*JVS( 409 )
  a = -W( 67 ) / JVS(          440  )
  W( 67 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 441 )
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 89 ) = W( 89 ) + a*JVS( 443 )
  W( 90 ) = W( 90 ) + a*JVS( 444 )
  W( 91 ) = W( 91 ) + a*JVS( 445 )
  W( 93 ) = W( 93 ) + a*JVS( 446 )
  W( 94 ) = W( 94 ) + a*JVS( 447 )
  W( 96 ) = W( 96 ) + a*JVS( 448 )
  a = -W( 68 ) / JVS(          451  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 452 )
  W( 88 ) = W( 88 ) + a*JVS( 453 )
  W( 89 ) = W( 89 ) + a*JVS( 454 )
  W( 90 ) = W( 90 ) + a*JVS( 455 )
  W( 91 ) = W( 91 ) + a*JVS( 456 )
  W( 93 ) = W( 93 ) + a*JVS( 457 )
  W( 94 ) = W( 94 ) + a*JVS( 458 )
  W( 95 ) = W( 95 ) + a*JVS( 459 )
  W( 96 ) = W( 96 ) + a*JVS( 460 )
  a = -W( 70 ) / JVS(          491  )
  W( 70 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 492 )
  W( 85 ) = W( 85 ) + a*JVS( 493 )
  W( 89 ) = W( 89 ) + a*JVS( 494 )
  W( 90 ) = W( 90 ) + a*JVS( 495 )
  W( 91 ) = W( 91 ) + a*JVS( 496 )
  W( 93 ) = W( 93 ) + a*JVS( 497 )
  W( 94 ) = W( 94 ) + a*JVS( 498 )
  W( 96 ) = W( 96 ) + a*JVS( 499 )
  a = -W( 71 ) / JVS(          501  )
  W( 71 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 502 )
  W( 88 ) = W( 88 ) + a*JVS( 503 )
  W( 89 ) = W( 89 ) + a*JVS( 504 )
  W( 90 ) = W( 90 ) + a*JVS( 505 )
  W( 91 ) = W( 91 ) + a*JVS( 506 )
  W( 93 ) = W( 93 ) + a*JVS( 507 )
  W( 94 ) = W( 94 ) + a*JVS( 508 )
  W( 95 ) = W( 95 ) + a*JVS( 509 )
  W( 96 ) = W( 96 ) + a*JVS( 510 )
  a = -W( 72 ) / JVS(          512  )
  W( 72 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 513 )
  W( 88 ) = W( 88 ) + a*JVS( 514 )
  W( 89 ) = W( 89 ) + a*JVS( 515 )
  W( 90 ) = W( 90 ) + a*JVS( 516 )
  W( 91 ) = W( 91 ) + a*JVS( 517 )
  W( 93 ) = W( 93 ) + a*JVS( 518 )
  W( 94 ) = W( 94 ) + a*JVS( 519 )
  W( 95 ) = W( 95 ) + a*JVS( 520 )
  W( 96 ) = W( 96 ) + a*JVS( 521 )
  JVS( 522) = W( 22 )
  JVS( 523) = W( 25 )
  JVS( 524) = W( 37 )
  JVS( 525) = W( 42 )
  JVS( 526) = W( 43 )
  JVS( 527) = W( 45 )
  JVS( 528) = W( 46 )
  JVS( 529) = W( 48 )
  JVS( 530) = W( 49 )
  JVS( 531) = W( 51 )
  JVS( 532) = W( 53 )
  JVS( 533) = W( 54 )
  JVS( 534) = W( 55 )
  JVS( 535) = W( 56 )
  JVS( 536) = W( 58 )
  JVS( 537) = W( 59 )
  JVS( 538) = W( 60 )
  JVS( 539) = W( 61 )
  JVS( 540) = W( 62 )
  JVS( 541) = W( 63 )
  JVS( 542) = W( 64 )
  JVS( 543) = W( 65 )
  JVS( 544) = W( 67 )
  JVS( 545) = W( 68 )
  JVS( 546) = W( 70 )
  JVS( 547) = W( 71 )
  JVS( 548) = W( 72 )
  JVS( 549) = W( 73 )
  JVS( 550) = W( 74 )
  JVS( 551) = W( 76 )
  JVS( 552) = W( 77 )
  JVS( 553) = W( 79 )
  JVS( 554) = W( 80 )
  JVS( 555) = W( 81 )
  JVS( 556) = W( 82 )
  JVS( 557) = W( 83 )
  JVS( 558) = W( 84 )
  JVS( 559) = W( 85 )
  JVS( 560) = W( 86 )
  JVS( 561) = W( 87 )
  JVS( 562) = W( 88 )
  JVS( 563) = W( 89 )
  JVS( 564) = W( 90 )
  JVS( 565) = W( 91 )
  JVS( 566) = W( 93 )
  JVS( 567) = W( 94 )
  JVS( 568) = W( 95 )
  JVS( 569) = W( 96 )
  JVS( 570) = W( 97 )
  IF ( ABS(  JVS( 572 )) < TINY(a) ) THEN
         IER = 74                                      
         RETURN
  END IF
   W( 58 ) = JVS( 571 )
   W( 74 ) = JVS( 572 )
   W( 79 ) = JVS( 573 )
   W( 85 ) = JVS( 574 )
   W( 89 ) = JVS( 575 )
   W( 90 ) = JVS( 576 )
   W( 91 ) = JVS( 577 )
   W( 93 ) = JVS( 578 )
   W( 94 ) = JVS( 579 )
   W( 95 ) = JVS( 580 )
   W( 96 ) = JVS( 581 )
  a = -W( 58 ) / JVS(          331  )
  W( 58 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 332 )
  W( 85 ) = W( 85 ) + a*JVS( 333 )
  W( 91 ) = W( 91 ) + a*JVS( 334 )
  W( 95 ) = W( 95 ) + a*JVS( 335 )
  W( 96 ) = W( 96 ) + a*JVS( 336 )
  JVS( 571) = W( 58 )
  JVS( 572) = W( 74 )
  JVS( 573) = W( 79 )
  JVS( 574) = W( 85 )
  JVS( 575) = W( 89 )
  JVS( 576) = W( 90 )
  JVS( 577) = W( 91 )
  JVS( 578) = W( 93 )
  JVS( 579) = W( 94 )
  JVS( 580) = W( 95 )
  JVS( 581) = W( 96 )
  IF ( ABS(  JVS( 590 )) < TINY(a) ) THEN
         IER = 75                                      
         RETURN
  END IF
   W( 31 ) = JVS( 582 )
   W( 57 ) = JVS( 583 )
   W( 58 ) = JVS( 584 )
   W( 63 ) = JVS( 585 )
   W( 65 ) = JVS( 586 )
   W( 71 ) = JVS( 587 )
   W( 72 ) = JVS( 588 )
   W( 74 ) = JVS( 589 )
   W( 75 ) = JVS( 590 )
   W( 76 ) = JVS( 591 )
   W( 77 ) = JVS( 592 )
   W( 79 ) = JVS( 593 )
   W( 83 ) = JVS( 594 )
   W( 84 ) = JVS( 595 )
   W( 85 ) = JVS( 596 )
   W( 86 ) = JVS( 597 )
   W( 88 ) = JVS( 598 )
   W( 89 ) = JVS( 599 )
   W( 90 ) = JVS( 600 )
   W( 91 ) = JVS( 601 )
   W( 93 ) = JVS( 602 )
   W( 94 ) = JVS( 603 )
   W( 95 ) = JVS( 604 )
   W( 96 ) = JVS( 605 )
  a = -W( 31 ) / JVS(          133  )
  W( 31 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 134 )
  W( 89 ) = W( 89 ) + a*JVS( 135 )
  W( 91 ) = W( 91 ) + a*JVS( 136 )
  W( 95 ) = W( 95 ) + a*JVS( 137 )
  W( 96 ) = W( 96 ) + a*JVS( 138 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 325 )
  W( 88 ) = W( 88 ) + a*JVS( 326 )
  W( 89 ) = W( 89 ) + a*JVS( 327 )
  W( 91 ) = W( 91 ) + a*JVS( 328 )
  W( 95 ) = W( 95 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  a = -W( 58 ) / JVS(          331  )
  W( 58 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 332 )
  W( 85 ) = W( 85 ) + a*JVS( 333 )
  W( 91 ) = W( 91 ) + a*JVS( 334 )
  W( 95 ) = W( 95 ) + a*JVS( 335 )
  W( 96 ) = W( 96 ) + a*JVS( 336 )
  a = -W( 63 ) / JVS(          377  )
  W( 63 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  W( 90 ) = W( 90 ) + a*JVS( 380 )
  W( 91 ) = W( 91 ) + a*JVS( 381 )
  W( 93 ) = W( 93 ) + a*JVS( 382 )
  W( 94 ) = W( 94 ) + a*JVS( 383 )
  W( 96 ) = W( 96 ) + a*JVS( 384 )
  a = -W( 65 ) / JVS(          402  )
  W( 65 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 403 )
  W( 89 ) = W( 89 ) + a*JVS( 404 )
  W( 90 ) = W( 90 ) + a*JVS( 405 )
  W( 91 ) = W( 91 ) + a*JVS( 406 )
  W( 93 ) = W( 93 ) + a*JVS( 407 )
  W( 94 ) = W( 94 ) + a*JVS( 408 )
  W( 96 ) = W( 96 ) + a*JVS( 409 )
  a = -W( 71 ) / JVS(          501  )
  W( 71 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 502 )
  W( 88 ) = W( 88 ) + a*JVS( 503 )
  W( 89 ) = W( 89 ) + a*JVS( 504 )
  W( 90 ) = W( 90 ) + a*JVS( 505 )
  W( 91 ) = W( 91 ) + a*JVS( 506 )
  W( 93 ) = W( 93 ) + a*JVS( 507 )
  W( 94 ) = W( 94 ) + a*JVS( 508 )
  W( 95 ) = W( 95 ) + a*JVS( 509 )
  W( 96 ) = W( 96 ) + a*JVS( 510 )
  a = -W( 72 ) / JVS(          512  )
  W( 72 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 513 )
  W( 88 ) = W( 88 ) + a*JVS( 514 )
  W( 89 ) = W( 89 ) + a*JVS( 515 )
  W( 90 ) = W( 90 ) + a*JVS( 516 )
  W( 91 ) = W( 91 ) + a*JVS( 517 )
  W( 93 ) = W( 93 ) + a*JVS( 518 )
  W( 94 ) = W( 94 ) + a*JVS( 519 )
  W( 95 ) = W( 95 ) + a*JVS( 520 )
  W( 96 ) = W( 96 ) + a*JVS( 521 )
  a = -W( 74 ) / JVS(          572  )
  W( 74 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 573 )
  W( 85 ) = W( 85 ) + a*JVS( 574 )
  W( 89 ) = W( 89 ) + a*JVS( 575 )
  W( 90 ) = W( 90 ) + a*JVS( 576 )
  W( 91 ) = W( 91 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  JVS( 582) = W( 31 )
  JVS( 583) = W( 57 )
  JVS( 584) = W( 58 )
  JVS( 585) = W( 63 )
  JVS( 586) = W( 65 )
  JVS( 587) = W( 71 )
  JVS( 588) = W( 72 )
  JVS( 589) = W( 74 )
  JVS( 590) = W( 75 )
  JVS( 591) = W( 76 )
  JVS( 592) = W( 77 )
  JVS( 593) = W( 79 )
  JVS( 594) = W( 83 )
  JVS( 595) = W( 84 )
  JVS( 596) = W( 85 )
  JVS( 597) = W( 86 )
  JVS( 598) = W( 88 )
  JVS( 599) = W( 89 )
  JVS( 600) = W( 90 )
  JVS( 601) = W( 91 )
  JVS( 602) = W( 93 )
  JVS( 603) = W( 94 )
  JVS( 604) = W( 95 )
  JVS( 605) = W( 96 )
  IF ( ABS(  JVS( 607 )) < TINY(a) ) THEN
         IER = 76                                      
         RETURN
  END IF
   W( 41 ) = JVS( 606 )
   W( 76 ) = JVS( 607 )
   W( 81 ) = JVS( 608 )
   W( 82 ) = JVS( 609 )
   W( 88 ) = JVS( 610 )
   W( 89 ) = JVS( 611 )
   W( 90 ) = JVS( 612 )
   W( 91 ) = JVS( 613 )
   W( 93 ) = JVS( 614 )
   W( 94 ) = JVS( 615 )
   W( 96 ) = JVS( 616 )
  a = -W( 41 ) / JVS(          206  )
  W( 41 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 207 )
  W( 82 ) = W( 82 ) + a*JVS( 208 )
  W( 88 ) = W( 88 ) + a*JVS( 209 )
  W( 91 ) = W( 91 ) + a*JVS( 210 )
  JVS( 606) = W( 41 )
  JVS( 607) = W( 76 )
  JVS( 608) = W( 81 )
  JVS( 609) = W( 82 )
  JVS( 610) = W( 88 )
  JVS( 611) = W( 89 )
  JVS( 612) = W( 90 )
  JVS( 613) = W( 91 )
  JVS( 614) = W( 93 )
  JVS( 615) = W( 94 )
  JVS( 616) = W( 96 )
  IF ( ABS(  JVS( 618 )) < TINY(a) ) THEN
         IER = 77                                      
         RETURN
  END IF
   W( 46 ) = JVS( 617 )
   W( 77 ) = JVS( 618 )
   W( 81 ) = JVS( 619 )
   W( 82 ) = JVS( 620 )
   W( 88 ) = JVS( 621 )
   W( 89 ) = JVS( 622 )
   W( 90 ) = JVS( 623 )
   W( 91 ) = JVS( 624 )
   W( 93 ) = JVS( 625 )
   W( 94 ) = JVS( 626 )
   W( 96 ) = JVS( 627 )
  a = -W( 46 ) / JVS(          238  )
  W( 46 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 239 )
  W( 82 ) = W( 82 ) + a*JVS( 240 )
  W( 88 ) = W( 88 ) + a*JVS( 241 )
  W( 91 ) = W( 91 ) + a*JVS( 242 )
  JVS( 617) = W( 46 )
  JVS( 618) = W( 77 )
  JVS( 619) = W( 81 )
  JVS( 620) = W( 82 )
  JVS( 621) = W( 88 )
  JVS( 622) = W( 89 )
  JVS( 623) = W( 90 )
  JVS( 624) = W( 91 )
  JVS( 625) = W( 93 )
  JVS( 626) = W( 94 )
  JVS( 627) = W( 96 )
  IF ( ABS(  JVS( 643 )) < TINY(a) ) THEN
         IER = 78                                      
         RETURN
  END IF
   W( 21 ) = JVS( 628 )
   W( 41 ) = JVS( 629 )
   W( 43 ) = JVS( 630 )
   W( 46 ) = JVS( 631 )
   W( 47 ) = JVS( 632 )
   W( 49 ) = JVS( 633 )
   W( 50 ) = JVS( 634 )
   W( 54 ) = JVS( 635 )
   W( 61 ) = JVS( 636 )
   W( 62 ) = JVS( 637 )
   W( 65 ) = JVS( 638 )
   W( 70 ) = JVS( 639 )
   W( 75 ) = JVS( 640 )
   W( 76 ) = JVS( 641 )
   W( 77 ) = JVS( 642 )
   W( 78 ) = JVS( 643 )
   W( 79 ) = JVS( 644 )
   W( 81 ) = JVS( 645 )
   W( 82 ) = JVS( 646 )
   W( 83 ) = JVS( 647 )
   W( 84 ) = JVS( 648 )
   W( 85 ) = JVS( 649 )
   W( 86 ) = JVS( 650 )
   W( 87 ) = JVS( 651 )
   W( 88 ) = JVS( 652 )
   W( 89 ) = JVS( 653 )
   W( 90 ) = JVS( 654 )
   W( 91 ) = JVS( 655 )
   W( 92 ) = JVS( 656 )
   W( 93 ) = JVS( 657 )
   W( 94 ) = JVS( 658 )
   W( 95 ) = JVS( 659 )
   W( 96 ) = JVS( 660 )
   W( 97 ) = JVS( 661 )
  a = -W( 21 ) / JVS(           88  )
  W( 21 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 89 )
  W( 91 ) = W( 91 ) + a*JVS( 90 )
  a = -W( 41 ) / JVS(          206  )
  W( 41 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 207 )
  W( 82 ) = W( 82 ) + a*JVS( 208 )
  W( 88 ) = W( 88 ) + a*JVS( 209 )
  W( 91 ) = W( 91 ) + a*JVS( 210 )
  a = -W( 43 ) / JVS(          217  )
  W( 43 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 218 )
  W( 91 ) = W( 91 ) + a*JVS( 219 )
  W( 96 ) = W( 96 ) + a*JVS( 220 )
  a = -W( 46 ) / JVS(          238  )
  W( 46 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 239 )
  W( 82 ) = W( 82 ) + a*JVS( 240 )
  W( 88 ) = W( 88 ) + a*JVS( 241 )
  W( 91 ) = W( 91 ) + a*JVS( 242 )
  a = -W( 47 ) / JVS(          243  )
  W( 47 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 244 )
  W( 91 ) = W( 91 ) + a*JVS( 245 )
  W( 96 ) = W( 96 ) + a*JVS( 246 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 255 )
  W( 82 ) = W( 82 ) + a*JVS( 256 )
  W( 88 ) = W( 88 ) + a*JVS( 257 )
  W( 91 ) = W( 91 ) + a*JVS( 258 )
  a = -W( 50 ) / JVS(          259  )
  W( 50 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 260 )
  W( 81 ) = W( 81 ) + a*JVS( 261 )
  W( 82 ) = W( 82 ) + a*JVS( 262 )
  W( 85 ) = W( 85 ) + a*JVS( 263 )
  W( 88 ) = W( 88 ) + a*JVS( 264 )
  W( 91 ) = W( 91 ) + a*JVS( 265 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  a = -W( 61 ) / JVS(          358  )
  W( 61 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 359 )
  W( 89 ) = W( 89 ) + a*JVS( 360 )
  W( 90 ) = W( 90 ) + a*JVS( 361 )
  W( 91 ) = W( 91 ) + a*JVS( 362 )
  W( 93 ) = W( 93 ) + a*JVS( 363 )
  W( 94 ) = W( 94 ) + a*JVS( 364 )
  W( 96 ) = W( 96 ) + a*JVS( 365 )
  a = -W( 62 ) / JVS(          367  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 368 )
  W( 88 ) = W( 88 ) + a*JVS( 369 )
  W( 89 ) = W( 89 ) + a*JVS( 370 )
  W( 90 ) = W( 90 ) + a*JVS( 371 )
  W( 91 ) = W( 91 ) + a*JVS( 372 )
  W( 93 ) = W( 93 ) + a*JVS( 373 )
  W( 94 ) = W( 94 ) + a*JVS( 374 )
  W( 96 ) = W( 96 ) + a*JVS( 375 )
  a = -W( 65 ) / JVS(          402  )
  W( 65 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 403 )
  W( 89 ) = W( 89 ) + a*JVS( 404 )
  W( 90 ) = W( 90 ) + a*JVS( 405 )
  W( 91 ) = W( 91 ) + a*JVS( 406 )
  W( 93 ) = W( 93 ) + a*JVS( 407 )
  W( 94 ) = W( 94 ) + a*JVS( 408 )
  W( 96 ) = W( 96 ) + a*JVS( 409 )
  a = -W( 70 ) / JVS(          491  )
  W( 70 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 492 )
  W( 85 ) = W( 85 ) + a*JVS( 493 )
  W( 89 ) = W( 89 ) + a*JVS( 494 )
  W( 90 ) = W( 90 ) + a*JVS( 495 )
  W( 91 ) = W( 91 ) + a*JVS( 496 )
  W( 93 ) = W( 93 ) + a*JVS( 497 )
  W( 94 ) = W( 94 ) + a*JVS( 498 )
  W( 96 ) = W( 96 ) + a*JVS( 499 )
  a = -W( 75 ) / JVS(          590  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 591 )
  W( 77 ) = W( 77 ) + a*JVS( 592 )
  W( 79 ) = W( 79 ) + a*JVS( 593 )
  W( 83 ) = W( 83 ) + a*JVS( 594 )
  W( 84 ) = W( 84 ) + a*JVS( 595 )
  W( 85 ) = W( 85 ) + a*JVS( 596 )
  W( 86 ) = W( 86 ) + a*JVS( 597 )
  W( 88 ) = W( 88 ) + a*JVS( 598 )
  W( 89 ) = W( 89 ) + a*JVS( 599 )
  W( 90 ) = W( 90 ) + a*JVS( 600 )
  W( 91 ) = W( 91 ) + a*JVS( 601 )
  W( 93 ) = W( 93 ) + a*JVS( 602 )
  W( 94 ) = W( 94 ) + a*JVS( 603 )
  W( 95 ) = W( 95 ) + a*JVS( 604 )
  W( 96 ) = W( 96 ) + a*JVS( 605 )
  a = -W( 76 ) / JVS(          607  )
  W( 76 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 608 )
  W( 82 ) = W( 82 ) + a*JVS( 609 )
  W( 88 ) = W( 88 ) + a*JVS( 610 )
  W( 89 ) = W( 89 ) + a*JVS( 611 )
  W( 90 ) = W( 90 ) + a*JVS( 612 )
  W( 91 ) = W( 91 ) + a*JVS( 613 )
  W( 93 ) = W( 93 ) + a*JVS( 614 )
  W( 94 ) = W( 94 ) + a*JVS( 615 )
  W( 96 ) = W( 96 ) + a*JVS( 616 )
  a = -W( 77 ) / JVS(          618  )
  W( 77 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 619 )
  W( 82 ) = W( 82 ) + a*JVS( 620 )
  W( 88 ) = W( 88 ) + a*JVS( 621 )
  W( 89 ) = W( 89 ) + a*JVS( 622 )
  W( 90 ) = W( 90 ) + a*JVS( 623 )
  W( 91 ) = W( 91 ) + a*JVS( 624 )
  W( 93 ) = W( 93 ) + a*JVS( 625 )
  W( 94 ) = W( 94 ) + a*JVS( 626 )
  W( 96 ) = W( 96 ) + a*JVS( 627 )
  JVS( 628) = W( 21 )
  JVS( 629) = W( 41 )
  JVS( 630) = W( 43 )
  JVS( 631) = W( 46 )
  JVS( 632) = W( 47 )
  JVS( 633) = W( 49 )
  JVS( 634) = W( 50 )
  JVS( 635) = W( 54 )
  JVS( 636) = W( 61 )
  JVS( 637) = W( 62 )
  JVS( 638) = W( 65 )
  JVS( 639) = W( 70 )
  JVS( 640) = W( 75 )
  JVS( 641) = W( 76 )
  JVS( 642) = W( 77 )
  JVS( 643) = W( 78 )
  JVS( 644) = W( 79 )
  JVS( 645) = W( 81 )
  JVS( 646) = W( 82 )
  JVS( 647) = W( 83 )
  JVS( 648) = W( 84 )
  JVS( 649) = W( 85 )
  JVS( 650) = W( 86 )
  JVS( 651) = W( 87 )
  JVS( 652) = W( 88 )
  JVS( 653) = W( 89 )
  JVS( 654) = W( 90 )
  JVS( 655) = W( 91 )
  JVS( 656) = W( 92 )
  JVS( 657) = W( 93 )
  JVS( 658) = W( 94 )
  JVS( 659) = W( 95 )
  JVS( 660) = W( 96 )
  JVS( 661) = W( 97 )
  IF ( ABS(  JVS( 664 )) < TINY(a) ) THEN
         IER = 79                                      
         RETURN
  END IF
   W( 71 ) = JVS( 662 )
   W( 72 ) = JVS( 663 )
   W( 79 ) = JVS( 664 )
   W( 85 ) = JVS( 665 )
   W( 88 ) = JVS( 666 )
   W( 89 ) = JVS( 667 )
   W( 90 ) = JVS( 668 )
   W( 91 ) = JVS( 669 )
   W( 93 ) = JVS( 670 )
   W( 94 ) = JVS( 671 )
   W( 95 ) = JVS( 672 )
   W( 96 ) = JVS( 673 )
  a = -W( 71 ) / JVS(          501  )
  W( 71 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 502 )
  W( 88 ) = W( 88 ) + a*JVS( 503 )
  W( 89 ) = W( 89 ) + a*JVS( 504 )
  W( 90 ) = W( 90 ) + a*JVS( 505 )
  W( 91 ) = W( 91 ) + a*JVS( 506 )
  W( 93 ) = W( 93 ) + a*JVS( 507 )
  W( 94 ) = W( 94 ) + a*JVS( 508 )
  W( 95 ) = W( 95 ) + a*JVS( 509 )
  W( 96 ) = W( 96 ) + a*JVS( 510 )
  a = -W( 72 ) / JVS(          512  )
  W( 72 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 513 )
  W( 88 ) = W( 88 ) + a*JVS( 514 )
  W( 89 ) = W( 89 ) + a*JVS( 515 )
  W( 90 ) = W( 90 ) + a*JVS( 516 )
  W( 91 ) = W( 91 ) + a*JVS( 517 )
  W( 93 ) = W( 93 ) + a*JVS( 518 )
  W( 94 ) = W( 94 ) + a*JVS( 519 )
  W( 95 ) = W( 95 ) + a*JVS( 520 )
  W( 96 ) = W( 96 ) + a*JVS( 521 )
  JVS( 662) = W( 71 )
  JVS( 663) = W( 72 )
  JVS( 664) = W( 79 )
  JVS( 665) = W( 85 )
  JVS( 666) = W( 88 )
  JVS( 667) = W( 89 )
  JVS( 668) = W( 90 )
  JVS( 669) = W( 91 )
  JVS( 670) = W( 93 )
  JVS( 671) = W( 94 )
  JVS( 672) = W( 95 )
  JVS( 673) = W( 96 )
  IF ( ABS(  JVS( 692 )) < TINY(a) ) THEN
         IER = 80                                      
         RETURN
  END IF
   W( 23 ) = JVS( 674 )
   W( 24 ) = JVS( 675 )
   W( 41 ) = JVS( 676 )
   W( 45 ) = JVS( 677 )
   W( 46 ) = JVS( 678 )
   W( 48 ) = JVS( 679 )
   W( 49 ) = JVS( 680 )
   W( 50 ) = JVS( 681 )
   W( 53 ) = JVS( 682 )
   W( 54 ) = JVS( 683 )
   W( 55 ) = JVS( 684 )
   W( 57 ) = JVS( 685 )
   W( 58 ) = JVS( 686 )
   W( 59 ) = JVS( 687 )
   W( 74 ) = JVS( 688 )
   W( 76 ) = JVS( 689 )
   W( 77 ) = JVS( 690 )
   W( 79 ) = JVS( 691 )
   W( 80 ) = JVS( 692 )
   W( 81 ) = JVS( 693 )
   W( 82 ) = JVS( 694 )
   W( 85 ) = JVS( 695 )
   W( 86 ) = JVS( 696 )
   W( 87 ) = JVS( 697 )
   W( 88 ) = JVS( 698 )
   W( 89 ) = JVS( 699 )
   W( 90 ) = JVS( 700 )
   W( 91 ) = JVS( 701 )
   W( 92 ) = JVS( 702 )
   W( 93 ) = JVS( 703 )
   W( 94 ) = JVS( 704 )
   W( 95 ) = JVS( 705 )
   W( 96 ) = JVS( 706 )
  a = -W( 23 ) / JVS(           94  )
  W( 23 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 95 )
  W( 91 ) = W( 91 ) + a*JVS( 96 )
  a = -W( 24 ) / JVS(           97  )
  W( 24 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 98 )
  W( 91 ) = W( 91 ) + a*JVS( 99 )
  a = -W( 41 ) / JVS(          206  )
  W( 41 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 207 )
  W( 82 ) = W( 82 ) + a*JVS( 208 )
  W( 88 ) = W( 88 ) + a*JVS( 209 )
  W( 91 ) = W( 91 ) + a*JVS( 210 )
  a = -W( 45 ) / JVS(          232  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 233 )
  W( 85 ) = W( 85 ) + a*JVS( 234 )
  W( 89 ) = W( 89 ) + a*JVS( 235 )
  W( 91 ) = W( 91 ) + a*JVS( 236 )
  W( 93 ) = W( 93 ) + a*JVS( 237 )
  a = -W( 46 ) / JVS(          238  )
  W( 46 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 239 )
  W( 82 ) = W( 82 ) + a*JVS( 240 )
  W( 88 ) = W( 88 ) + a*JVS( 241 )
  W( 91 ) = W( 91 ) + a*JVS( 242 )
  a = -W( 48 ) / JVS(          247  )
  W( 48 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 248 )
  W( 85 ) = W( 85 ) + a*JVS( 249 )
  W( 91 ) = W( 91 ) + a*JVS( 250 )
  W( 93 ) = W( 93 ) + a*JVS( 251 )
  W( 95 ) = W( 95 ) + a*JVS( 252 )
  W( 96 ) = W( 96 ) + a*JVS( 253 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 255 )
  W( 82 ) = W( 82 ) + a*JVS( 256 )
  W( 88 ) = W( 88 ) + a*JVS( 257 )
  W( 91 ) = W( 91 ) + a*JVS( 258 )
  a = -W( 50 ) / JVS(          259  )
  W( 50 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 260 )
  W( 81 ) = W( 81 ) + a*JVS( 261 )
  W( 82 ) = W( 82 ) + a*JVS( 262 )
  W( 85 ) = W( 85 ) + a*JVS( 263 )
  W( 88 ) = W( 88 ) + a*JVS( 264 )
  W( 91 ) = W( 91 ) + a*JVS( 265 )
  a = -W( 53 ) / JVS(          297  )
  W( 53 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 298 )
  W( 88 ) = W( 88 ) + a*JVS( 299 )
  W( 91 ) = W( 91 ) + a*JVS( 300 )
  W( 96 ) = W( 96 ) + a*JVS( 301 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 325 )
  W( 88 ) = W( 88 ) + a*JVS( 326 )
  W( 89 ) = W( 89 ) + a*JVS( 327 )
  W( 91 ) = W( 91 ) + a*JVS( 328 )
  W( 95 ) = W( 95 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  a = -W( 58 ) / JVS(          331  )
  W( 58 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 332 )
  W( 85 ) = W( 85 ) + a*JVS( 333 )
  W( 91 ) = W( 91 ) + a*JVS( 334 )
  W( 95 ) = W( 95 ) + a*JVS( 335 )
  W( 96 ) = W( 96 ) + a*JVS( 336 )
  a = -W( 59 ) / JVS(          337  )
  W( 59 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 338 )
  W( 88 ) = W( 88 ) + a*JVS( 339 )
  W( 91 ) = W( 91 ) + a*JVS( 340 )
  W( 96 ) = W( 96 ) + a*JVS( 341 )
  a = -W( 74 ) / JVS(          572  )
  W( 74 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 573 )
  W( 85 ) = W( 85 ) + a*JVS( 574 )
  W( 89 ) = W( 89 ) + a*JVS( 575 )
  W( 90 ) = W( 90 ) + a*JVS( 576 )
  W( 91 ) = W( 91 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  a = -W( 76 ) / JVS(          607  )
  W( 76 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 608 )
  W( 82 ) = W( 82 ) + a*JVS( 609 )
  W( 88 ) = W( 88 ) + a*JVS( 610 )
  W( 89 ) = W( 89 ) + a*JVS( 611 )
  W( 90 ) = W( 90 ) + a*JVS( 612 )
  W( 91 ) = W( 91 ) + a*JVS( 613 )
  W( 93 ) = W( 93 ) + a*JVS( 614 )
  W( 94 ) = W( 94 ) + a*JVS( 615 )
  W( 96 ) = W( 96 ) + a*JVS( 616 )
  a = -W( 77 ) / JVS(          618  )
  W( 77 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 619 )
  W( 82 ) = W( 82 ) + a*JVS( 620 )
  W( 88 ) = W( 88 ) + a*JVS( 621 )
  W( 89 ) = W( 89 ) + a*JVS( 622 )
  W( 90 ) = W( 90 ) + a*JVS( 623 )
  W( 91 ) = W( 91 ) + a*JVS( 624 )
  W( 93 ) = W( 93 ) + a*JVS( 625 )
  W( 94 ) = W( 94 ) + a*JVS( 626 )
  W( 96 ) = W( 96 ) + a*JVS( 627 )
  a = -W( 79 ) / JVS(          664  )
  W( 79 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 665 )
  W( 88 ) = W( 88 ) + a*JVS( 666 )
  W( 89 ) = W( 89 ) + a*JVS( 667 )
  W( 90 ) = W( 90 ) + a*JVS( 668 )
  W( 91 ) = W( 91 ) + a*JVS( 669 )
  W( 93 ) = W( 93 ) + a*JVS( 670 )
  W( 94 ) = W( 94 ) + a*JVS( 671 )
  W( 95 ) = W( 95 ) + a*JVS( 672 )
  W( 96 ) = W( 96 ) + a*JVS( 673 )
  JVS( 674) = W( 23 )
  JVS( 675) = W( 24 )
  JVS( 676) = W( 41 )
  JVS( 677) = W( 45 )
  JVS( 678) = W( 46 )
  JVS( 679) = W( 48 )
  JVS( 680) = W( 49 )
  JVS( 681) = W( 50 )
  JVS( 682) = W( 53 )
  JVS( 683) = W( 54 )
  JVS( 684) = W( 55 )
  JVS( 685) = W( 57 )
  JVS( 686) = W( 58 )
  JVS( 687) = W( 59 )
  JVS( 688) = W( 74 )
  JVS( 689) = W( 76 )
  JVS( 690) = W( 77 )
  JVS( 691) = W( 79 )
  JVS( 692) = W( 80 )
  JVS( 693) = W( 81 )
  JVS( 694) = W( 82 )
  JVS( 695) = W( 85 )
  JVS( 696) = W( 86 )
  JVS( 697) = W( 87 )
  JVS( 698) = W( 88 )
  JVS( 699) = W( 89 )
  JVS( 700) = W( 90 )
  JVS( 701) = W( 91 )
  JVS( 702) = W( 92 )
  JVS( 703) = W( 93 )
  JVS( 704) = W( 94 )
  JVS( 705) = W( 95 )
  JVS( 706) = W( 96 )
  IF ( ABS(  JVS( 709 )) < TINY(a) ) THEN
         IER = 81                                      
         RETURN
  END IF
   W( 63 ) = JVS( 707 )
   W( 67 ) = JVS( 708 )
   W( 81 ) = JVS( 709 )
   W( 85 ) = JVS( 710 )
   W( 88 ) = JVS( 711 )
   W( 89 ) = JVS( 712 )
   W( 90 ) = JVS( 713 )
   W( 91 ) = JVS( 714 )
   W( 93 ) = JVS( 715 )
   W( 94 ) = JVS( 716 )
   W( 96 ) = JVS( 717 )
  a = -W( 63 ) / JVS(          377  )
  W( 63 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  W( 90 ) = W( 90 ) + a*JVS( 380 )
  W( 91 ) = W( 91 ) + a*JVS( 381 )
  W( 93 ) = W( 93 ) + a*JVS( 382 )
  W( 94 ) = W( 94 ) + a*JVS( 383 )
  W( 96 ) = W( 96 ) + a*JVS( 384 )
  a = -W( 67 ) / JVS(          440  )
  W( 67 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 441 )
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 89 ) = W( 89 ) + a*JVS( 443 )
  W( 90 ) = W( 90 ) + a*JVS( 444 )
  W( 91 ) = W( 91 ) + a*JVS( 445 )
  W( 93 ) = W( 93 ) + a*JVS( 446 )
  W( 94 ) = W( 94 ) + a*JVS( 447 )
  W( 96 ) = W( 96 ) + a*JVS( 448 )
  JVS( 707) = W( 63 )
  JVS( 708) = W( 67 )
  JVS( 709) = W( 81 )
  JVS( 710) = W( 85 )
  JVS( 711) = W( 88 )
  JVS( 712) = W( 89 )
  JVS( 713) = W( 90 )
  JVS( 714) = W( 91 )
  JVS( 715) = W( 93 )
  JVS( 716) = W( 94 )
  JVS( 717) = W( 96 )
  IF ( ABS(  JVS( 722 )) < TINY(a) ) THEN
         IER = 82                                      
         RETURN
  END IF
   W( 51 ) = JVS( 718 )
   W( 55 ) = JVS( 719 )
   W( 67 ) = JVS( 720 )
   W( 81 ) = JVS( 721 )
   W( 82 ) = JVS( 722 )
   W( 85 ) = JVS( 723 )
   W( 88 ) = JVS( 724 )
   W( 89 ) = JVS( 725 )
   W( 90 ) = JVS( 726 )
   W( 91 ) = JVS( 727 )
   W( 93 ) = JVS( 728 )
   W( 94 ) = JVS( 729 )
   W( 96 ) = JVS( 730 )
  a = -W( 51 ) / JVS(          266  )
  W( 51 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 267 )
  W( 91 ) = W( 91 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 67 ) / JVS(          440  )
  W( 67 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 441 )
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 89 ) = W( 89 ) + a*JVS( 443 )
  W( 90 ) = W( 90 ) + a*JVS( 444 )
  W( 91 ) = W( 91 ) + a*JVS( 445 )
  W( 93 ) = W( 93 ) + a*JVS( 446 )
  W( 94 ) = W( 94 ) + a*JVS( 447 )
  W( 96 ) = W( 96 ) + a*JVS( 448 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  JVS( 718) = W( 51 )
  JVS( 719) = W( 55 )
  JVS( 720) = W( 67 )
  JVS( 721) = W( 81 )
  JVS( 722) = W( 82 )
  JVS( 723) = W( 85 )
  JVS( 724) = W( 88 )
  JVS( 725) = W( 89 )
  JVS( 726) = W( 90 )
  JVS( 727) = W( 91 )
  JVS( 728) = W( 93 )
  JVS( 729) = W( 94 )
  JVS( 730) = W( 96 )
  IF ( ABS(  JVS( 739 )) < TINY(a) ) THEN
         IER = 83                                      
         RETURN
  END IF
   W( 43 ) = JVS( 731 )
   W( 47 ) = JVS( 732 )
   W( 51 ) = JVS( 733 )
   W( 53 ) = JVS( 734 )
   W( 54 ) = JVS( 735 )
   W( 55 ) = JVS( 736 )
   W( 81 ) = JVS( 737 )
   W( 82 ) = JVS( 738 )
   W( 83 ) = JVS( 739 )
   W( 84 ) = JVS( 740 )
   W( 85 ) = JVS( 741 )
   W( 88 ) = JVS( 742 )
   W( 89 ) = JVS( 743 )
   W( 90 ) = JVS( 744 )
   W( 91 ) = JVS( 745 )
   W( 93 ) = JVS( 746 )
   W( 94 ) = JVS( 747 )
   W( 96 ) = JVS( 748 )
  a = -W( 43 ) / JVS(          217  )
  W( 43 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 218 )
  W( 91 ) = W( 91 ) + a*JVS( 219 )
  W( 96 ) = W( 96 ) + a*JVS( 220 )
  a = -W( 47 ) / JVS(          243  )
  W( 47 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 244 )
  W( 91 ) = W( 91 ) + a*JVS( 245 )
  W( 96 ) = W( 96 ) + a*JVS( 246 )
  a = -W( 51 ) / JVS(          266  )
  W( 51 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 267 )
  W( 91 ) = W( 91 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 53 ) / JVS(          297  )
  W( 53 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 298 )
  W( 88 ) = W( 88 ) + a*JVS( 299 )
  W( 91 ) = W( 91 ) + a*JVS( 300 )
  W( 96 ) = W( 96 ) + a*JVS( 301 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  JVS( 731) = W( 43 )
  JVS( 732) = W( 47 )
  JVS( 733) = W( 51 )
  JVS( 734) = W( 53 )
  JVS( 735) = W( 54 )
  JVS( 736) = W( 55 )
  JVS( 737) = W( 81 )
  JVS( 738) = W( 82 )
  JVS( 739) = W( 83 )
  JVS( 740) = W( 84 )
  JVS( 741) = W( 85 )
  JVS( 742) = W( 88 )
  JVS( 743) = W( 89 )
  JVS( 744) = W( 90 )
  JVS( 745) = W( 91 )
  JVS( 746) = W( 93 )
  JVS( 747) = W( 94 )
  JVS( 748) = W( 96 )
  IF ( ABS(  JVS( 758 )) < TINY(a) ) THEN
         IER = 84                                      
         RETURN
  END IF
   W( 43 ) = JVS( 749 )
   W( 47 ) = JVS( 750 )
   W( 51 ) = JVS( 751 )
   W( 53 ) = JVS( 752 )
   W( 54 ) = JVS( 753 )
   W( 55 ) = JVS( 754 )
   W( 81 ) = JVS( 755 )
   W( 82 ) = JVS( 756 )
   W( 83 ) = JVS( 757 )
   W( 84 ) = JVS( 758 )
   W( 85 ) = JVS( 759 )
   W( 88 ) = JVS( 760 )
   W( 89 ) = JVS( 761 )
   W( 90 ) = JVS( 762 )
   W( 91 ) = JVS( 763 )
   W( 93 ) = JVS( 764 )
   W( 94 ) = JVS( 765 )
   W( 96 ) = JVS( 766 )
  a = -W( 43 ) / JVS(          217  )
  W( 43 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 218 )
  W( 91 ) = W( 91 ) + a*JVS( 219 )
  W( 96 ) = W( 96 ) + a*JVS( 220 )
  a = -W( 47 ) / JVS(          243  )
  W( 47 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 244 )
  W( 91 ) = W( 91 ) + a*JVS( 245 )
  W( 96 ) = W( 96 ) + a*JVS( 246 )
  a = -W( 51 ) / JVS(          266  )
  W( 51 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 267 )
  W( 91 ) = W( 91 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 53 ) / JVS(          297  )
  W( 53 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 298 )
  W( 88 ) = W( 88 ) + a*JVS( 299 )
  W( 91 ) = W( 91 ) + a*JVS( 300 )
  W( 96 ) = W( 96 ) + a*JVS( 301 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  a = -W( 83 ) / JVS(          739  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 740 )
  W( 85 ) = W( 85 ) + a*JVS( 741 )
  W( 88 ) = W( 88 ) + a*JVS( 742 )
  W( 89 ) = W( 89 ) + a*JVS( 743 )
  W( 90 ) = W( 90 ) + a*JVS( 744 )
  W( 91 ) = W( 91 ) + a*JVS( 745 )
  W( 93 ) = W( 93 ) + a*JVS( 746 )
  W( 94 ) = W( 94 ) + a*JVS( 747 )
  W( 96 ) = W( 96 ) + a*JVS( 748 )
  JVS( 749) = W( 43 )
  JVS( 750) = W( 47 )
  JVS( 751) = W( 51 )
  JVS( 752) = W( 53 )
  JVS( 753) = W( 54 )
  JVS( 754) = W( 55 )
  JVS( 755) = W( 81 )
  JVS( 756) = W( 82 )
  JVS( 757) = W( 83 )
  JVS( 758) = W( 84 )
  JVS( 759) = W( 85 )
  JVS( 760) = W( 88 )
  JVS( 761) = W( 89 )
  JVS( 762) = W( 90 )
  JVS( 763) = W( 91 )
  JVS( 764) = W( 93 )
  JVS( 765) = W( 94 )
  JVS( 766) = W( 96 )
  IF ( ABS(  JVS( 786 )) < TINY(a) ) THEN
         IER = 85                                      
         RETURN
  END IF
   W( 30 ) = JVS( 767 )
   W( 38 ) = JVS( 768 )
   W( 39 ) = JVS( 769 )
   W( 40 ) = JVS( 770 )
   W( 43 ) = JVS( 771 )
   W( 47 ) = JVS( 772 )
   W( 51 ) = JVS( 773 )
   W( 53 ) = JVS( 774 )
   W( 54 ) = JVS( 775 )
   W( 55 ) = JVS( 776 )
   W( 57 ) = JVS( 777 )
   W( 58 ) = JVS( 778 )
   W( 59 ) = JVS( 779 )
   W( 64 ) = JVS( 780 )
   W( 67 ) = JVS( 781 )
   W( 74 ) = JVS( 782 )
   W( 79 ) = JVS( 783 )
   W( 81 ) = JVS( 784 )
   W( 82 ) = JVS( 785 )
   W( 85 ) = JVS( 786 )
   W( 88 ) = JVS( 787 )
   W( 89 ) = JVS( 788 )
   W( 90 ) = JVS( 789 )
   W( 91 ) = JVS( 790 )
   W( 93 ) = JVS( 791 )
   W( 94 ) = JVS( 792 )
   W( 95 ) = JVS( 793 )
   W( 96 ) = JVS( 794 )
  a = -W( 30 ) / JVS(          127  )
  W( 30 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 128 )
  W( 85 ) = W( 85 ) + a*JVS( 129 )
  W( 94 ) = W( 94 ) + a*JVS( 130 )
  W( 95 ) = W( 95 ) + a*JVS( 131 )
  W( 96 ) = W( 96 ) + a*JVS( 132 )
  a = -W( 38 ) / JVS(          190  )
  W( 38 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 191 )
  W( 88 ) = W( 88 ) + a*JVS( 192 )
  W( 91 ) = W( 91 ) + a*JVS( 193 )
  W( 95 ) = W( 95 ) + a*JVS( 194 )
  a = -W( 39 ) / JVS(          196  )
  W( 39 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 197 )
  W( 88 ) = W( 88 ) + a*JVS( 198 )
  W( 91 ) = W( 91 ) + a*JVS( 199 )
  W( 95 ) = W( 95 ) + a*JVS( 200 )
  a = -W( 40 ) / JVS(          201  )
  W( 40 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 202 )
  W( 85 ) = W( 85 ) + a*JVS( 203 )
  W( 91 ) = W( 91 ) + a*JVS( 204 )
  W( 95 ) = W( 95 ) + a*JVS( 205 )
  a = -W( 43 ) / JVS(          217  )
  W( 43 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 218 )
  W( 91 ) = W( 91 ) + a*JVS( 219 )
  W( 96 ) = W( 96 ) + a*JVS( 220 )
  a = -W( 47 ) / JVS(          243  )
  W( 47 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 244 )
  W( 91 ) = W( 91 ) + a*JVS( 245 )
  W( 96 ) = W( 96 ) + a*JVS( 246 )
  a = -W( 51 ) / JVS(          266  )
  W( 51 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 267 )
  W( 91 ) = W( 91 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 53 ) / JVS(          297  )
  W( 53 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 298 )
  W( 88 ) = W( 88 ) + a*JVS( 299 )
  W( 91 ) = W( 91 ) + a*JVS( 300 )
  W( 96 ) = W( 96 ) + a*JVS( 301 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 325 )
  W( 88 ) = W( 88 ) + a*JVS( 326 )
  W( 89 ) = W( 89 ) + a*JVS( 327 )
  W( 91 ) = W( 91 ) + a*JVS( 328 )
  W( 95 ) = W( 95 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  a = -W( 58 ) / JVS(          331  )
  W( 58 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 332 )
  W( 85 ) = W( 85 ) + a*JVS( 333 )
  W( 91 ) = W( 91 ) + a*JVS( 334 )
  W( 95 ) = W( 95 ) + a*JVS( 335 )
  W( 96 ) = W( 96 ) + a*JVS( 336 )
  a = -W( 59 ) / JVS(          337  )
  W( 59 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 338 )
  W( 88 ) = W( 88 ) + a*JVS( 339 )
  W( 91 ) = W( 91 ) + a*JVS( 340 )
  W( 96 ) = W( 96 ) + a*JVS( 341 )
  a = -W( 64 ) / JVS(          391  )
  W( 64 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 392 )
  W( 85 ) = W( 85 ) + a*JVS( 393 )
  W( 88 ) = W( 88 ) + a*JVS( 394 )
  W( 89 ) = W( 89 ) + a*JVS( 395 )
  W( 90 ) = W( 90 ) + a*JVS( 396 )
  W( 91 ) = W( 91 ) + a*JVS( 397 )
  W( 93 ) = W( 93 ) + a*JVS( 398 )
  W( 94 ) = W( 94 ) + a*JVS( 399 )
  W( 96 ) = W( 96 ) + a*JVS( 400 )
  a = -W( 67 ) / JVS(          440  )
  W( 67 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 441 )
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 89 ) = W( 89 ) + a*JVS( 443 )
  W( 90 ) = W( 90 ) + a*JVS( 444 )
  W( 91 ) = W( 91 ) + a*JVS( 445 )
  W( 93 ) = W( 93 ) + a*JVS( 446 )
  W( 94 ) = W( 94 ) + a*JVS( 447 )
  W( 96 ) = W( 96 ) + a*JVS( 448 )
  a = -W( 74 ) / JVS(          572  )
  W( 74 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 573 )
  W( 85 ) = W( 85 ) + a*JVS( 574 )
  W( 89 ) = W( 89 ) + a*JVS( 575 )
  W( 90 ) = W( 90 ) + a*JVS( 576 )
  W( 91 ) = W( 91 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  a = -W( 79 ) / JVS(          664  )
  W( 79 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 665 )
  W( 88 ) = W( 88 ) + a*JVS( 666 )
  W( 89 ) = W( 89 ) + a*JVS( 667 )
  W( 90 ) = W( 90 ) + a*JVS( 668 )
  W( 91 ) = W( 91 ) + a*JVS( 669 )
  W( 93 ) = W( 93 ) + a*JVS( 670 )
  W( 94 ) = W( 94 ) + a*JVS( 671 )
  W( 95 ) = W( 95 ) + a*JVS( 672 )
  W( 96 ) = W( 96 ) + a*JVS( 673 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  JVS( 767) = W( 30 )
  JVS( 768) = W( 38 )
  JVS( 769) = W( 39 )
  JVS( 770) = W( 40 )
  JVS( 771) = W( 43 )
  JVS( 772) = W( 47 )
  JVS( 773) = W( 51 )
  JVS( 774) = W( 53 )
  JVS( 775) = W( 54 )
  JVS( 776) = W( 55 )
  JVS( 777) = W( 57 )
  JVS( 778) = W( 58 )
  JVS( 779) = W( 59 )
  JVS( 780) = W( 64 )
  JVS( 781) = W( 67 )
  JVS( 782) = W( 74 )
  JVS( 783) = W( 79 )
  JVS( 784) = W( 81 )
  JVS( 785) = W( 82 )
  JVS( 786) = W( 85 )
  JVS( 787) = W( 88 )
  JVS( 788) = W( 89 )
  JVS( 789) = W( 90 )
  JVS( 790) = W( 91 )
  JVS( 791) = W( 93 )
  JVS( 792) = W( 94 )
  JVS( 793) = W( 95 )
  JVS( 794) = W( 96 )
  IF ( ABS(  JVS( 805 )) < TINY(a) ) THEN
         IER = 86                                      
         RETURN
  END IF
   W( 49 ) = JVS( 795 )
   W( 75 ) = JVS( 796 )
   W( 76 ) = JVS( 797 )
   W( 77 ) = JVS( 798 )
   W( 79 ) = JVS( 799 )
   W( 81 ) = JVS( 800 )
   W( 82 ) = JVS( 801 )
   W( 83 ) = JVS( 802 )
   W( 84 ) = JVS( 803 )
   W( 85 ) = JVS( 804 )
   W( 86 ) = JVS( 805 )
   W( 88 ) = JVS( 806 )
   W( 89 ) = JVS( 807 )
   W( 90 ) = JVS( 808 )
   W( 91 ) = JVS( 809 )
   W( 92 ) = JVS( 810 )
   W( 93 ) = JVS( 811 )
   W( 94 ) = JVS( 812 )
   W( 95 ) = JVS( 813 )
   W( 96 ) = JVS( 814 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 255 )
  W( 82 ) = W( 82 ) + a*JVS( 256 )
  W( 88 ) = W( 88 ) + a*JVS( 257 )
  W( 91 ) = W( 91 ) + a*JVS( 258 )
  a = -W( 75 ) / JVS(          590  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 591 )
  W( 77 ) = W( 77 ) + a*JVS( 592 )
  W( 79 ) = W( 79 ) + a*JVS( 593 )
  W( 83 ) = W( 83 ) + a*JVS( 594 )
  W( 84 ) = W( 84 ) + a*JVS( 595 )
  W( 85 ) = W( 85 ) + a*JVS( 596 )
  W( 86 ) = W( 86 ) + a*JVS( 597 )
  W( 88 ) = W( 88 ) + a*JVS( 598 )
  W( 89 ) = W( 89 ) + a*JVS( 599 )
  W( 90 ) = W( 90 ) + a*JVS( 600 )
  W( 91 ) = W( 91 ) + a*JVS( 601 )
  W( 93 ) = W( 93 ) + a*JVS( 602 )
  W( 94 ) = W( 94 ) + a*JVS( 603 )
  W( 95 ) = W( 95 ) + a*JVS( 604 )
  W( 96 ) = W( 96 ) + a*JVS( 605 )
  a = -W( 76 ) / JVS(          607  )
  W( 76 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 608 )
  W( 82 ) = W( 82 ) + a*JVS( 609 )
  W( 88 ) = W( 88 ) + a*JVS( 610 )
  W( 89 ) = W( 89 ) + a*JVS( 611 )
  W( 90 ) = W( 90 ) + a*JVS( 612 )
  W( 91 ) = W( 91 ) + a*JVS( 613 )
  W( 93 ) = W( 93 ) + a*JVS( 614 )
  W( 94 ) = W( 94 ) + a*JVS( 615 )
  W( 96 ) = W( 96 ) + a*JVS( 616 )
  a = -W( 77 ) / JVS(          618  )
  W( 77 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 619 )
  W( 82 ) = W( 82 ) + a*JVS( 620 )
  W( 88 ) = W( 88 ) + a*JVS( 621 )
  W( 89 ) = W( 89 ) + a*JVS( 622 )
  W( 90 ) = W( 90 ) + a*JVS( 623 )
  W( 91 ) = W( 91 ) + a*JVS( 624 )
  W( 93 ) = W( 93 ) + a*JVS( 625 )
  W( 94 ) = W( 94 ) + a*JVS( 626 )
  W( 96 ) = W( 96 ) + a*JVS( 627 )
  a = -W( 79 ) / JVS(          664  )
  W( 79 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 665 )
  W( 88 ) = W( 88 ) + a*JVS( 666 )
  W( 89 ) = W( 89 ) + a*JVS( 667 )
  W( 90 ) = W( 90 ) + a*JVS( 668 )
  W( 91 ) = W( 91 ) + a*JVS( 669 )
  W( 93 ) = W( 93 ) + a*JVS( 670 )
  W( 94 ) = W( 94 ) + a*JVS( 671 )
  W( 95 ) = W( 95 ) + a*JVS( 672 )
  W( 96 ) = W( 96 ) + a*JVS( 673 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  a = -W( 83 ) / JVS(          739  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 740 )
  W( 85 ) = W( 85 ) + a*JVS( 741 )
  W( 88 ) = W( 88 ) + a*JVS( 742 )
  W( 89 ) = W( 89 ) + a*JVS( 743 )
  W( 90 ) = W( 90 ) + a*JVS( 744 )
  W( 91 ) = W( 91 ) + a*JVS( 745 )
  W( 93 ) = W( 93 ) + a*JVS( 746 )
  W( 94 ) = W( 94 ) + a*JVS( 747 )
  W( 96 ) = W( 96 ) + a*JVS( 748 )
  a = -W( 84 ) / JVS(          758  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 759 )
  W( 88 ) = W( 88 ) + a*JVS( 760 )
  W( 89 ) = W( 89 ) + a*JVS( 761 )
  W( 90 ) = W( 90 ) + a*JVS( 762 )
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 93 ) = W( 93 ) + a*JVS( 764 )
  W( 94 ) = W( 94 ) + a*JVS( 765 )
  W( 96 ) = W( 96 ) + a*JVS( 766 )
  a = -W( 85 ) / JVS(          786  )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 787 )
  W( 89 ) = W( 89 ) + a*JVS( 788 )
  W( 90 ) = W( 90 ) + a*JVS( 789 )
  W( 91 ) = W( 91 ) + a*JVS( 790 )
  W( 93 ) = W( 93 ) + a*JVS( 791 )
  W( 94 ) = W( 94 ) + a*JVS( 792 )
  W( 95 ) = W( 95 ) + a*JVS( 793 )
  W( 96 ) = W( 96 ) + a*JVS( 794 )
  JVS( 795) = W( 49 )
  JVS( 796) = W( 75 )
  JVS( 797) = W( 76 )
  JVS( 798) = W( 77 )
  JVS( 799) = W( 79 )
  JVS( 800) = W( 81 )
  JVS( 801) = W( 82 )
  JVS( 802) = W( 83 )
  JVS( 803) = W( 84 )
  JVS( 804) = W( 85 )
  JVS( 805) = W( 86 )
  JVS( 806) = W( 88 )
  JVS( 807) = W( 89 )
  JVS( 808) = W( 90 )
  JVS( 809) = W( 91 )
  JVS( 810) = W( 92 )
  JVS( 811) = W( 93 )
  JVS( 812) = W( 94 )
  JVS( 813) = W( 95 )
  JVS( 814) = W( 96 )
  IF ( ABS(  JVS( 831 )) < TINY(a) ) THEN
         IER = 87                                      
         RETURN
  END IF
   W( 47 ) = JVS( 815 )
   W( 51 ) = JVS( 816 )
   W( 54 ) = JVS( 817 )
   W( 55 ) = JVS( 818 )
   W( 66 ) = JVS( 819 )
   W( 70 ) = JVS( 820 )
   W( 75 ) = JVS( 821 )
   W( 76 ) = JVS( 822 )
   W( 77 ) = JVS( 823 )
   W( 79 ) = JVS( 824 )
   W( 81 ) = JVS( 825 )
   W( 82 ) = JVS( 826 )
   W( 83 ) = JVS( 827 )
   W( 84 ) = JVS( 828 )
   W( 85 ) = JVS( 829 )
   W( 86 ) = JVS( 830 )
   W( 87 ) = JVS( 831 )
   W( 88 ) = JVS( 832 )
   W( 89 ) = JVS( 833 )
   W( 90 ) = JVS( 834 )
   W( 91 ) = JVS( 835 )
   W( 92 ) = JVS( 836 )
   W( 93 ) = JVS( 837 )
   W( 94 ) = JVS( 838 )
   W( 95 ) = JVS( 839 )
   W( 96 ) = JVS( 840 )
  a = -W( 47 ) / JVS(          243  )
  W( 47 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 244 )
  W( 91 ) = W( 91 ) + a*JVS( 245 )
  W( 96 ) = W( 96 ) + a*JVS( 246 )
  a = -W( 51 ) / JVS(          266  )
  W( 51 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 267 )
  W( 91 ) = W( 91 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 66 ) / JVS(          417  )
  W( 66 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 75 ) = W( 75 ) + a*JVS( 419 )
  W( 76 ) = W( 76 ) + a*JVS( 420 )
  W( 77 ) = W( 77 ) + a*JVS( 421 )
  W( 79 ) = W( 79 ) + a*JVS( 422 )
  W( 81 ) = W( 81 ) + a*JVS( 423 )
  W( 82 ) = W( 82 ) + a*JVS( 424 )
  W( 83 ) = W( 83 ) + a*JVS( 425 )
  W( 84 ) = W( 84 ) + a*JVS( 426 )
  W( 85 ) = W( 85 ) + a*JVS( 427 )
  W( 86 ) = W( 86 ) + a*JVS( 428 )
  W( 87 ) = W( 87 ) + a*JVS( 429 )
  W( 88 ) = W( 88 ) + a*JVS( 430 )
  W( 89 ) = W( 89 ) + a*JVS( 431 )
  W( 90 ) = W( 90 ) + a*JVS( 432 )
  W( 91 ) = W( 91 ) + a*JVS( 433 )
  W( 92 ) = W( 92 ) + a*JVS( 434 )
  W( 93 ) = W( 93 ) + a*JVS( 435 )
  W( 94 ) = W( 94 ) + a*JVS( 436 )
  W( 96 ) = W( 96 ) + a*JVS( 437 )
  a = -W( 70 ) / JVS(          491  )
  W( 70 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 492 )
  W( 85 ) = W( 85 ) + a*JVS( 493 )
  W( 89 ) = W( 89 ) + a*JVS( 494 )
  W( 90 ) = W( 90 ) + a*JVS( 495 )
  W( 91 ) = W( 91 ) + a*JVS( 496 )
  W( 93 ) = W( 93 ) + a*JVS( 497 )
  W( 94 ) = W( 94 ) + a*JVS( 498 )
  W( 96 ) = W( 96 ) + a*JVS( 499 )
  a = -W( 75 ) / JVS(          590  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 591 )
  W( 77 ) = W( 77 ) + a*JVS( 592 )
  W( 79 ) = W( 79 ) + a*JVS( 593 )
  W( 83 ) = W( 83 ) + a*JVS( 594 )
  W( 84 ) = W( 84 ) + a*JVS( 595 )
  W( 85 ) = W( 85 ) + a*JVS( 596 )
  W( 86 ) = W( 86 ) + a*JVS( 597 )
  W( 88 ) = W( 88 ) + a*JVS( 598 )
  W( 89 ) = W( 89 ) + a*JVS( 599 )
  W( 90 ) = W( 90 ) + a*JVS( 600 )
  W( 91 ) = W( 91 ) + a*JVS( 601 )
  W( 93 ) = W( 93 ) + a*JVS( 602 )
  W( 94 ) = W( 94 ) + a*JVS( 603 )
  W( 95 ) = W( 95 ) + a*JVS( 604 )
  W( 96 ) = W( 96 ) + a*JVS( 605 )
  a = -W( 76 ) / JVS(          607  )
  W( 76 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 608 )
  W( 82 ) = W( 82 ) + a*JVS( 609 )
  W( 88 ) = W( 88 ) + a*JVS( 610 )
  W( 89 ) = W( 89 ) + a*JVS( 611 )
  W( 90 ) = W( 90 ) + a*JVS( 612 )
  W( 91 ) = W( 91 ) + a*JVS( 613 )
  W( 93 ) = W( 93 ) + a*JVS( 614 )
  W( 94 ) = W( 94 ) + a*JVS( 615 )
  W( 96 ) = W( 96 ) + a*JVS( 616 )
  a = -W( 77 ) / JVS(          618  )
  W( 77 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 619 )
  W( 82 ) = W( 82 ) + a*JVS( 620 )
  W( 88 ) = W( 88 ) + a*JVS( 621 )
  W( 89 ) = W( 89 ) + a*JVS( 622 )
  W( 90 ) = W( 90 ) + a*JVS( 623 )
  W( 91 ) = W( 91 ) + a*JVS( 624 )
  W( 93 ) = W( 93 ) + a*JVS( 625 )
  W( 94 ) = W( 94 ) + a*JVS( 626 )
  W( 96 ) = W( 96 ) + a*JVS( 627 )
  a = -W( 79 ) / JVS(          664  )
  W( 79 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 665 )
  W( 88 ) = W( 88 ) + a*JVS( 666 )
  W( 89 ) = W( 89 ) + a*JVS( 667 )
  W( 90 ) = W( 90 ) + a*JVS( 668 )
  W( 91 ) = W( 91 ) + a*JVS( 669 )
  W( 93 ) = W( 93 ) + a*JVS( 670 )
  W( 94 ) = W( 94 ) + a*JVS( 671 )
  W( 95 ) = W( 95 ) + a*JVS( 672 )
  W( 96 ) = W( 96 ) + a*JVS( 673 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  a = -W( 83 ) / JVS(          739  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 740 )
  W( 85 ) = W( 85 ) + a*JVS( 741 )
  W( 88 ) = W( 88 ) + a*JVS( 742 )
  W( 89 ) = W( 89 ) + a*JVS( 743 )
  W( 90 ) = W( 90 ) + a*JVS( 744 )
  W( 91 ) = W( 91 ) + a*JVS( 745 )
  W( 93 ) = W( 93 ) + a*JVS( 746 )
  W( 94 ) = W( 94 ) + a*JVS( 747 )
  W( 96 ) = W( 96 ) + a*JVS( 748 )
  a = -W( 84 ) / JVS(          758  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 759 )
  W( 88 ) = W( 88 ) + a*JVS( 760 )
  W( 89 ) = W( 89 ) + a*JVS( 761 )
  W( 90 ) = W( 90 ) + a*JVS( 762 )
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 93 ) = W( 93 ) + a*JVS( 764 )
  W( 94 ) = W( 94 ) + a*JVS( 765 )
  W( 96 ) = W( 96 ) + a*JVS( 766 )
  a = -W( 85 ) / JVS(          786  )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 787 )
  W( 89 ) = W( 89 ) + a*JVS( 788 )
  W( 90 ) = W( 90 ) + a*JVS( 789 )
  W( 91 ) = W( 91 ) + a*JVS( 790 )
  W( 93 ) = W( 93 ) + a*JVS( 791 )
  W( 94 ) = W( 94 ) + a*JVS( 792 )
  W( 95 ) = W( 95 ) + a*JVS( 793 )
  W( 96 ) = W( 96 ) + a*JVS( 794 )
  a = -W( 86 ) / JVS(          805  )
  W( 86 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 806 )
  W( 89 ) = W( 89 ) + a*JVS( 807 )
  W( 90 ) = W( 90 ) + a*JVS( 808 )
  W( 91 ) = W( 91 ) + a*JVS( 809 )
  W( 92 ) = W( 92 ) + a*JVS( 810 )
  W( 93 ) = W( 93 ) + a*JVS( 811 )
  W( 94 ) = W( 94 ) + a*JVS( 812 )
  W( 95 ) = W( 95 ) + a*JVS( 813 )
  W( 96 ) = W( 96 ) + a*JVS( 814 )
  JVS( 815) = W( 47 )
  JVS( 816) = W( 51 )
  JVS( 817) = W( 54 )
  JVS( 818) = W( 55 )
  JVS( 819) = W( 66 )
  JVS( 820) = W( 70 )
  JVS( 821) = W( 75 )
  JVS( 822) = W( 76 )
  JVS( 823) = W( 77 )
  JVS( 824) = W( 79 )
  JVS( 825) = W( 81 )
  JVS( 826) = W( 82 )
  JVS( 827) = W( 83 )
  JVS( 828) = W( 84 )
  JVS( 829) = W( 85 )
  JVS( 830) = W( 86 )
  JVS( 831) = W( 87 )
  JVS( 832) = W( 88 )
  JVS( 833) = W( 89 )
  JVS( 834) = W( 90 )
  JVS( 835) = W( 91 )
  JVS( 836) = W( 92 )
  JVS( 837) = W( 93 )
  JVS( 838) = W( 94 )
  JVS( 839) = W( 95 )
  JVS( 840) = W( 96 )
  IF ( ABS(  JVS( 871 )) < TINY(a) ) THEN
         IER = 88                                      
         RETURN
  END IF
   W( 13 ) = JVS( 841 )
   W( 18 ) = JVS( 842 )
   W( 19 ) = JVS( 843 )
   W( 23 ) = JVS( 844 )
   W( 24 ) = JVS( 845 )
   W( 32 ) = JVS( 846 )
   W( 34 ) = JVS( 847 )
   W( 35 ) = JVS( 848 )
   W( 36 ) = JVS( 849 )
   W( 41 ) = JVS( 850 )
   W( 46 ) = JVS( 851 )
   W( 49 ) = JVS( 852 )
   W( 50 ) = JVS( 853 )
   W( 53 ) = JVS( 854 )
   W( 54 ) = JVS( 855 )
   W( 59 ) = JVS( 856 )
   W( 73 ) = JVS( 857 )
   W( 74 ) = JVS( 858 )
   W( 76 ) = JVS( 859 )
   W( 77 ) = JVS( 860 )
   W( 78 ) = JVS( 861 )
   W( 79 ) = JVS( 862 )
   W( 80 ) = JVS( 863 )
   W( 81 ) = JVS( 864 )
   W( 82 ) = JVS( 865 )
   W( 83 ) = JVS( 866 )
   W( 84 ) = JVS( 867 )
   W( 85 ) = JVS( 868 )
   W( 86 ) = JVS( 869 )
   W( 87 ) = JVS( 870 )
   W( 88 ) = JVS( 871 )
   W( 89 ) = JVS( 872 )
   W( 90 ) = JVS( 873 )
   W( 91 ) = JVS( 874 )
   W( 92 ) = JVS( 875 )
   W( 93 ) = JVS( 876 )
   W( 94 ) = JVS( 877 )
   W( 95 ) = JVS( 878 )
   W( 96 ) = JVS( 879 )
   W( 97 ) = JVS( 880 )
  a = -W( 13 ) / JVS(           65  )
  W( 13 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 66 )
  a = -W( 18 ) / JVS(           78  )
  W( 18 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 79 )
  W( 95 ) = W( 95 ) + a*JVS( 80 )
  a = -W( 19 ) / JVS(           81  )
  W( 19 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 82 )
  W( 89 ) = W( 89 ) + a*JVS( 83 )
  a = -W( 23 ) / JVS(           94  )
  W( 23 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 95 )
  W( 91 ) = W( 91 ) + a*JVS( 96 )
  a = -W( 24 ) / JVS(           97  )
  W( 24 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 98 )
  W( 91 ) = W( 91 ) + a*JVS( 99 )
  a = -W( 32 ) / JVS(          139  )
  W( 32 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 140 )
  W( 88 ) = W( 88 ) + a*JVS( 141 )
  W( 89 ) = W( 89 ) + a*JVS( 142 )
  W( 94 ) = W( 94 ) + a*JVS( 143 )
  a = -W( 34 ) / JVS(          157  )
  W( 34 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 158 )
  W( 41 ) = W( 41 ) + a*JVS( 159 )
  W( 46 ) = W( 46 ) + a*JVS( 160 )
  W( 49 ) = W( 49 ) + a*JVS( 161 )
  W( 50 ) = W( 50 ) + a*JVS( 162 )
  W( 59 ) = W( 59 ) + a*JVS( 163 )
  W( 73 ) = W( 73 ) + a*JVS( 164 )
  W( 78 ) = W( 78 ) + a*JVS( 165 )
  W( 81 ) = W( 81 ) + a*JVS( 166 )
  W( 88 ) = W( 88 ) + a*JVS( 167 )
  W( 91 ) = W( 91 ) + a*JVS( 168 )
  a = -W( 35 ) / JVS(          169  )
  W( 35 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 170 )
  W( 81 ) = W( 81 ) + a*JVS( 171 )
  W( 82 ) = W( 82 ) + a*JVS( 172 )
  W( 85 ) = W( 85 ) + a*JVS( 173 )
  W( 88 ) = W( 88 ) + a*JVS( 174 )
  W( 91 ) = W( 91 ) + a*JVS( 175 )
  a = -W( 36 ) / JVS(          176  )
  W( 36 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 177 )
  W( 59 ) = W( 59 ) + a*JVS( 178 )
  W( 81 ) = W( 81 ) + a*JVS( 179 )
  W( 82 ) = W( 82 ) + a*JVS( 180 )
  W( 88 ) = W( 88 ) + a*JVS( 181 )
  W( 91 ) = W( 91 ) + a*JVS( 182 )
  a = -W( 41 ) / JVS(          206  )
  W( 41 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 207 )
  W( 82 ) = W( 82 ) + a*JVS( 208 )
  W( 88 ) = W( 88 ) + a*JVS( 209 )
  W( 91 ) = W( 91 ) + a*JVS( 210 )
  a = -W( 46 ) / JVS(          238  )
  W( 46 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 239 )
  W( 82 ) = W( 82 ) + a*JVS( 240 )
  W( 88 ) = W( 88 ) + a*JVS( 241 )
  W( 91 ) = W( 91 ) + a*JVS( 242 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 255 )
  W( 82 ) = W( 82 ) + a*JVS( 256 )
  W( 88 ) = W( 88 ) + a*JVS( 257 )
  W( 91 ) = W( 91 ) + a*JVS( 258 )
  a = -W( 50 ) / JVS(          259  )
  W( 50 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 260 )
  W( 81 ) = W( 81 ) + a*JVS( 261 )
  W( 82 ) = W( 82 ) + a*JVS( 262 )
  W( 85 ) = W( 85 ) + a*JVS( 263 )
  W( 88 ) = W( 88 ) + a*JVS( 264 )
  W( 91 ) = W( 91 ) + a*JVS( 265 )
  a = -W( 53 ) / JVS(          297  )
  W( 53 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 298 )
  W( 88 ) = W( 88 ) + a*JVS( 299 )
  W( 91 ) = W( 91 ) + a*JVS( 300 )
  W( 96 ) = W( 96 ) + a*JVS( 301 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  a = -W( 59 ) / JVS(          337  )
  W( 59 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 338 )
  W( 88 ) = W( 88 ) + a*JVS( 339 )
  W( 91 ) = W( 91 ) + a*JVS( 340 )
  W( 96 ) = W( 96 ) + a*JVS( 341 )
  a = -W( 73 ) / JVS(          549  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 550 )
  W( 76 ) = W( 76 ) + a*JVS( 551 )
  W( 77 ) = W( 77 ) + a*JVS( 552 )
  W( 79 ) = W( 79 ) + a*JVS( 553 )
  W( 80 ) = W( 80 ) + a*JVS( 554 )
  W( 81 ) = W( 81 ) + a*JVS( 555 )
  W( 82 ) = W( 82 ) + a*JVS( 556 )
  W( 83 ) = W( 83 ) + a*JVS( 557 )
  W( 84 ) = W( 84 ) + a*JVS( 558 )
  W( 85 ) = W( 85 ) + a*JVS( 559 )
  W( 86 ) = W( 86 ) + a*JVS( 560 )
  W( 87 ) = W( 87 ) + a*JVS( 561 )
  W( 88 ) = W( 88 ) + a*JVS( 562 )
  W( 89 ) = W( 89 ) + a*JVS( 563 )
  W( 90 ) = W( 90 ) + a*JVS( 564 )
  W( 91 ) = W( 91 ) + a*JVS( 565 )
  W( 93 ) = W( 93 ) + a*JVS( 566 )
  W( 94 ) = W( 94 ) + a*JVS( 567 )
  W( 95 ) = W( 95 ) + a*JVS( 568 )
  W( 96 ) = W( 96 ) + a*JVS( 569 )
  W( 97 ) = W( 97 ) + a*JVS( 570 )
  a = -W( 74 ) / JVS(          572  )
  W( 74 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 573 )
  W( 85 ) = W( 85 ) + a*JVS( 574 )
  W( 89 ) = W( 89 ) + a*JVS( 575 )
  W( 90 ) = W( 90 ) + a*JVS( 576 )
  W( 91 ) = W( 91 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  a = -W( 76 ) / JVS(          607  )
  W( 76 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 608 )
  W( 82 ) = W( 82 ) + a*JVS( 609 )
  W( 88 ) = W( 88 ) + a*JVS( 610 )
  W( 89 ) = W( 89 ) + a*JVS( 611 )
  W( 90 ) = W( 90 ) + a*JVS( 612 )
  W( 91 ) = W( 91 ) + a*JVS( 613 )
  W( 93 ) = W( 93 ) + a*JVS( 614 )
  W( 94 ) = W( 94 ) + a*JVS( 615 )
  W( 96 ) = W( 96 ) + a*JVS( 616 )
  a = -W( 77 ) / JVS(          618  )
  W( 77 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 619 )
  W( 82 ) = W( 82 ) + a*JVS( 620 )
  W( 88 ) = W( 88 ) + a*JVS( 621 )
  W( 89 ) = W( 89 ) + a*JVS( 622 )
  W( 90 ) = W( 90 ) + a*JVS( 623 )
  W( 91 ) = W( 91 ) + a*JVS( 624 )
  W( 93 ) = W( 93 ) + a*JVS( 625 )
  W( 94 ) = W( 94 ) + a*JVS( 626 )
  W( 96 ) = W( 96 ) + a*JVS( 627 )
  a = -W( 78 ) / JVS(          643  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 644 )
  W( 81 ) = W( 81 ) + a*JVS( 645 )
  W( 82 ) = W( 82 ) + a*JVS( 646 )
  W( 83 ) = W( 83 ) + a*JVS( 647 )
  W( 84 ) = W( 84 ) + a*JVS( 648 )
  W( 85 ) = W( 85 ) + a*JVS( 649 )
  W( 86 ) = W( 86 ) + a*JVS( 650 )
  W( 87 ) = W( 87 ) + a*JVS( 651 )
  W( 88 ) = W( 88 ) + a*JVS( 652 )
  W( 89 ) = W( 89 ) + a*JVS( 653 )
  W( 90 ) = W( 90 ) + a*JVS( 654 )
  W( 91 ) = W( 91 ) + a*JVS( 655 )
  W( 92 ) = W( 92 ) + a*JVS( 656 )
  W( 93 ) = W( 93 ) + a*JVS( 657 )
  W( 94 ) = W( 94 ) + a*JVS( 658 )
  W( 95 ) = W( 95 ) + a*JVS( 659 )
  W( 96 ) = W( 96 ) + a*JVS( 660 )
  W( 97 ) = W( 97 ) + a*JVS( 661 )
  a = -W( 79 ) / JVS(          664  )
  W( 79 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 665 )
  W( 88 ) = W( 88 ) + a*JVS( 666 )
  W( 89 ) = W( 89 ) + a*JVS( 667 )
  W( 90 ) = W( 90 ) + a*JVS( 668 )
  W( 91 ) = W( 91 ) + a*JVS( 669 )
  W( 93 ) = W( 93 ) + a*JVS( 670 )
  W( 94 ) = W( 94 ) + a*JVS( 671 )
  W( 95 ) = W( 95 ) + a*JVS( 672 )
  W( 96 ) = W( 96 ) + a*JVS( 673 )
  a = -W( 80 ) / JVS(          692  )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 693 )
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 85 ) = W( 85 ) + a*JVS( 695 )
  W( 86 ) = W( 86 ) + a*JVS( 696 )
  W( 87 ) = W( 87 ) + a*JVS( 697 )
  W( 88 ) = W( 88 ) + a*JVS( 698 )
  W( 89 ) = W( 89 ) + a*JVS( 699 )
  W( 90 ) = W( 90 ) + a*JVS( 700 )
  W( 91 ) = W( 91 ) + a*JVS( 701 )
  W( 92 ) = W( 92 ) + a*JVS( 702 )
  W( 93 ) = W( 93 ) + a*JVS( 703 )
  W( 94 ) = W( 94 ) + a*JVS( 704 )
  W( 95 ) = W( 95 ) + a*JVS( 705 )
  W( 96 ) = W( 96 ) + a*JVS( 706 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  a = -W( 83 ) / JVS(          739  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 740 )
  W( 85 ) = W( 85 ) + a*JVS( 741 )
  W( 88 ) = W( 88 ) + a*JVS( 742 )
  W( 89 ) = W( 89 ) + a*JVS( 743 )
  W( 90 ) = W( 90 ) + a*JVS( 744 )
  W( 91 ) = W( 91 ) + a*JVS( 745 )
  W( 93 ) = W( 93 ) + a*JVS( 746 )
  W( 94 ) = W( 94 ) + a*JVS( 747 )
  W( 96 ) = W( 96 ) + a*JVS( 748 )
  a = -W( 84 ) / JVS(          758  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 759 )
  W( 88 ) = W( 88 ) + a*JVS( 760 )
  W( 89 ) = W( 89 ) + a*JVS( 761 )
  W( 90 ) = W( 90 ) + a*JVS( 762 )
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 93 ) = W( 93 ) + a*JVS( 764 )
  W( 94 ) = W( 94 ) + a*JVS( 765 )
  W( 96 ) = W( 96 ) + a*JVS( 766 )
  a = -W( 85 ) / JVS(          786  )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 787 )
  W( 89 ) = W( 89 ) + a*JVS( 788 )
  W( 90 ) = W( 90 ) + a*JVS( 789 )
  W( 91 ) = W( 91 ) + a*JVS( 790 )
  W( 93 ) = W( 93 ) + a*JVS( 791 )
  W( 94 ) = W( 94 ) + a*JVS( 792 )
  W( 95 ) = W( 95 ) + a*JVS( 793 )
  W( 96 ) = W( 96 ) + a*JVS( 794 )
  a = -W( 86 ) / JVS(          805  )
  W( 86 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 806 )
  W( 89 ) = W( 89 ) + a*JVS( 807 )
  W( 90 ) = W( 90 ) + a*JVS( 808 )
  W( 91 ) = W( 91 ) + a*JVS( 809 )
  W( 92 ) = W( 92 ) + a*JVS( 810 )
  W( 93 ) = W( 93 ) + a*JVS( 811 )
  W( 94 ) = W( 94 ) + a*JVS( 812 )
  W( 95 ) = W( 95 ) + a*JVS( 813 )
  W( 96 ) = W( 96 ) + a*JVS( 814 )
  a = -W( 87 ) / JVS(          831  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 832 )
  W( 89 ) = W( 89 ) + a*JVS( 833 )
  W( 90 ) = W( 90 ) + a*JVS( 834 )
  W( 91 ) = W( 91 ) + a*JVS( 835 )
  W( 92 ) = W( 92 ) + a*JVS( 836 )
  W( 93 ) = W( 93 ) + a*JVS( 837 )
  W( 94 ) = W( 94 ) + a*JVS( 838 )
  W( 95 ) = W( 95 ) + a*JVS( 839 )
  W( 96 ) = W( 96 ) + a*JVS( 840 )
  JVS( 841) = W( 13 )
  JVS( 842) = W( 18 )
  JVS( 843) = W( 19 )
  JVS( 844) = W( 23 )
  JVS( 845) = W( 24 )
  JVS( 846) = W( 32 )
  JVS( 847) = W( 34 )
  JVS( 848) = W( 35 )
  JVS( 849) = W( 36 )
  JVS( 850) = W( 41 )
  JVS( 851) = W( 46 )
  JVS( 852) = W( 49 )
  JVS( 853) = W( 50 )
  JVS( 854) = W( 53 )
  JVS( 855) = W( 54 )
  JVS( 856) = W( 59 )
  JVS( 857) = W( 73 )
  JVS( 858) = W( 74 )
  JVS( 859) = W( 76 )
  JVS( 860) = W( 77 )
  JVS( 861) = W( 78 )
  JVS( 862) = W( 79 )
  JVS( 863) = W( 80 )
  JVS( 864) = W( 81 )
  JVS( 865) = W( 82 )
  JVS( 866) = W( 83 )
  JVS( 867) = W( 84 )
  JVS( 868) = W( 85 )
  JVS( 869) = W( 86 )
  JVS( 870) = W( 87 )
  JVS( 871) = W( 88 )
  JVS( 872) = W( 89 )
  JVS( 873) = W( 90 )
  JVS( 874) = W( 91 )
  JVS( 875) = W( 92 )
  JVS( 876) = W( 93 )
  JVS( 877) = W( 94 )
  JVS( 878) = W( 95 )
  JVS( 879) = W( 96 )
  JVS( 880) = W( 97 )
  IF ( ABS(  JVS( 940 )) < TINY(a) ) THEN
         IER = 89                                      
         RETURN
  END IF
   W( 15 ) = JVS( 881 )
   W( 21 ) = JVS( 882 )
   W( 23 ) = JVS( 883 )
   W( 24 ) = JVS( 884 )
   W( 25 ) = JVS( 885 )
   W( 27 ) = JVS( 886 )
   W( 28 ) = JVS( 887 )
   W( 31 ) = JVS( 888 )
   W( 32 ) = JVS( 889 )
   W( 33 ) = JVS( 890 )
   W( 36 ) = JVS( 891 )
   W( 37 ) = JVS( 892 )
   W( 38 ) = JVS( 893 )
   W( 39 ) = JVS( 894 )
   W( 40 ) = JVS( 895 )
   W( 41 ) = JVS( 896 )
   W( 42 ) = JVS( 897 )
   W( 45 ) = JVS( 898 )
   W( 46 ) = JVS( 899 )
   W( 47 ) = JVS( 900 )
   W( 49 ) = JVS( 901 )
   W( 50 ) = JVS( 902 )
   W( 51 ) = JVS( 903 )
   W( 52 ) = JVS( 904 )
   W( 53 ) = JVS( 905 )
   W( 54 ) = JVS( 906 )
   W( 55 ) = JVS( 907 )
   W( 56 ) = JVS( 908 )
   W( 57 ) = JVS( 909 )
   W( 58 ) = JVS( 910 )
   W( 59 ) = JVS( 911 )
   W( 60 ) = JVS( 912 )
   W( 61 ) = JVS( 913 )
   W( 62 ) = JVS( 914 )
   W( 63 ) = JVS( 915 )
   W( 64 ) = JVS( 916 )
   W( 65 ) = JVS( 917 )
   W( 67 ) = JVS( 918 )
   W( 68 ) = JVS( 919 )
   W( 69 ) = JVS( 920 )
   W( 70 ) = JVS( 921 )
   W( 71 ) = JVS( 922 )
   W( 72 ) = JVS( 923 )
   W( 73 ) = JVS( 924 )
   W( 74 ) = JVS( 925 )
   W( 75 ) = JVS( 926 )
   W( 76 ) = JVS( 927 )
   W( 77 ) = JVS( 928 )
   W( 78 ) = JVS( 929 )
   W( 79 ) = JVS( 930 )
   W( 80 ) = JVS( 931 )
   W( 81 ) = JVS( 932 )
   W( 82 ) = JVS( 933 )
   W( 83 ) = JVS( 934 )
   W( 84 ) = JVS( 935 )
   W( 85 ) = JVS( 936 )
   W( 86 ) = JVS( 937 )
   W( 87 ) = JVS( 938 )
   W( 88 ) = JVS( 939 )
   W( 89 ) = JVS( 940 )
   W( 90 ) = JVS( 941 )
   W( 91 ) = JVS( 942 )
   W( 92 ) = JVS( 943 )
   W( 93 ) = JVS( 944 )
   W( 94 ) = JVS( 945 )
   W( 95 ) = JVS( 946 )
   W( 96 ) = JVS( 947 )
   W( 97 ) = JVS( 948 )
  a = -W( 15 ) / JVS(           69  )
  W( 15 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 70 )
  a = -W( 21 ) / JVS(           88  )
  W( 21 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 89 )
  W( 91 ) = W( 91 ) + a*JVS( 90 )
  a = -W( 23 ) / JVS(           94  )
  W( 23 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 95 )
  W( 91 ) = W( 91 ) + a*JVS( 96 )
  a = -W( 24 ) / JVS(           97  )
  W( 24 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 98 )
  W( 91 ) = W( 91 ) + a*JVS( 99 )
  a = -W( 25 ) / JVS(          100  )
  W( 25 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 101 )
  W( 90 ) = W( 90 ) + a*JVS( 102 )
  W( 91 ) = W( 91 ) + a*JVS( 103 )
  a = -W( 27 ) / JVS(          108  )
  W( 27 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 109 )
  W( 91 ) = W( 91 ) + a*JVS( 110 )
  W( 95 ) = W( 95 ) + a*JVS( 111 )
  a = -W( 28 ) / JVS(          113  )
  W( 28 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 114 )
  W( 56 ) = W( 56 ) + a*JVS( 115 )
  W( 91 ) = W( 91 ) + a*JVS( 116 )
  W( 94 ) = W( 94 ) + a*JVS( 117 )
  W( 95 ) = W( 95 ) + a*JVS( 118 )
  a = -W( 31 ) / JVS(          133  )
  W( 31 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 134 )
  W( 89 ) = W( 89 ) + a*JVS( 135 )
  W( 91 ) = W( 91 ) + a*JVS( 136 )
  W( 95 ) = W( 95 ) + a*JVS( 137 )
  W( 96 ) = W( 96 ) + a*JVS( 138 )
  a = -W( 32 ) / JVS(          139  )
  W( 32 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 140 )
  W( 88 ) = W( 88 ) + a*JVS( 141 )
  W( 89 ) = W( 89 ) + a*JVS( 142 )
  W( 94 ) = W( 94 ) + a*JVS( 143 )
  a = -W( 33 ) / JVS(          144  )
  W( 33 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 145 )
  W( 51 ) = W( 51 ) + a*JVS( 146 )
  W( 54 ) = W( 54 ) + a*JVS( 147 )
  W( 55 ) = W( 55 ) + a*JVS( 148 )
  W( 59 ) = W( 59 ) + a*JVS( 149 )
  W( 81 ) = W( 81 ) + a*JVS( 150 )
  W( 82 ) = W( 82 ) + a*JVS( 151 )
  W( 85 ) = W( 85 ) + a*JVS( 152 )
  W( 89 ) = W( 89 ) + a*JVS( 153 )
  W( 91 ) = W( 91 ) + a*JVS( 154 )
  a = -W( 36 ) / JVS(          176  )
  W( 36 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 177 )
  W( 59 ) = W( 59 ) + a*JVS( 178 )
  W( 81 ) = W( 81 ) + a*JVS( 179 )
  W( 82 ) = W( 82 ) + a*JVS( 180 )
  W( 88 ) = W( 88 ) + a*JVS( 181 )
  W( 91 ) = W( 91 ) + a*JVS( 182 )
  a = -W( 37 ) / JVS(          183  )
  W( 37 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 184 )
  W( 58 ) = W( 58 ) + a*JVS( 185 )
  W( 87 ) = W( 87 ) + a*JVS( 186 )
  W( 90 ) = W( 90 ) + a*JVS( 187 )
  W( 91 ) = W( 91 ) + a*JVS( 188 )
  a = -W( 38 ) / JVS(          190  )
  W( 38 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 191 )
  W( 88 ) = W( 88 ) + a*JVS( 192 )
  W( 91 ) = W( 91 ) + a*JVS( 193 )
  W( 95 ) = W( 95 ) + a*JVS( 194 )
  a = -W( 39 ) / JVS(          196  )
  W( 39 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 197 )
  W( 88 ) = W( 88 ) + a*JVS( 198 )
  W( 91 ) = W( 91 ) + a*JVS( 199 )
  W( 95 ) = W( 95 ) + a*JVS( 200 )
  a = -W( 40 ) / JVS(          201  )
  W( 40 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 202 )
  W( 85 ) = W( 85 ) + a*JVS( 203 )
  W( 91 ) = W( 91 ) + a*JVS( 204 )
  W( 95 ) = W( 95 ) + a*JVS( 205 )
  a = -W( 41 ) / JVS(          206  )
  W( 41 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 207 )
  W( 82 ) = W( 82 ) + a*JVS( 208 )
  W( 88 ) = W( 88 ) + a*JVS( 209 )
  W( 91 ) = W( 91 ) + a*JVS( 210 )
  a = -W( 42 ) / JVS(          211  )
  W( 42 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 212 )
  W( 67 ) = W( 67 ) + a*JVS( 213 )
  W( 91 ) = W( 91 ) + a*JVS( 214 )
  W( 94 ) = W( 94 ) + a*JVS( 215 )
  W( 96 ) = W( 96 ) + a*JVS( 216 )
  a = -W( 45 ) / JVS(          232  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 233 )
  W( 85 ) = W( 85 ) + a*JVS( 234 )
  W( 89 ) = W( 89 ) + a*JVS( 235 )
  W( 91 ) = W( 91 ) + a*JVS( 236 )
  W( 93 ) = W( 93 ) + a*JVS( 237 )
  a = -W( 46 ) / JVS(          238  )
  W( 46 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 239 )
  W( 82 ) = W( 82 ) + a*JVS( 240 )
  W( 88 ) = W( 88 ) + a*JVS( 241 )
  W( 91 ) = W( 91 ) + a*JVS( 242 )
  a = -W( 47 ) / JVS(          243  )
  W( 47 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 244 )
  W( 91 ) = W( 91 ) + a*JVS( 245 )
  W( 96 ) = W( 96 ) + a*JVS( 246 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 255 )
  W( 82 ) = W( 82 ) + a*JVS( 256 )
  W( 88 ) = W( 88 ) + a*JVS( 257 )
  W( 91 ) = W( 91 ) + a*JVS( 258 )
  a = -W( 50 ) / JVS(          259  )
  W( 50 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 260 )
  W( 81 ) = W( 81 ) + a*JVS( 261 )
  W( 82 ) = W( 82 ) + a*JVS( 262 )
  W( 85 ) = W( 85 ) + a*JVS( 263 )
  W( 88 ) = W( 88 ) + a*JVS( 264 )
  W( 91 ) = W( 91 ) + a*JVS( 265 )
  a = -W( 51 ) / JVS(          266  )
  W( 51 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 267 )
  W( 91 ) = W( 91 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 52 ) / JVS(          276  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 277 )
  W( 54 ) = W( 54 ) + a*JVS( 278 )
  W( 55 ) = W( 55 ) + a*JVS( 279 )
  W( 56 ) = W( 56 ) + a*JVS( 280 )
  W( 58 ) = W( 58 ) + a*JVS( 281 )
  W( 59 ) = W( 59 ) + a*JVS( 282 )
  W( 60 ) = W( 60 ) + a*JVS( 283 )
  W( 64 ) = W( 64 ) + a*JVS( 284 )
  W( 67 ) = W( 67 ) + a*JVS( 285 )
  W( 69 ) = W( 69 ) + a*JVS( 286 )
  W( 73 ) = W( 73 ) + a*JVS( 287 )
  W( 78 ) = W( 78 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  W( 81 ) = W( 81 ) + a*JVS( 290 )
  W( 82 ) = W( 82 ) + a*JVS( 291 )
  W( 85 ) = W( 85 ) + a*JVS( 292 )
  W( 88 ) = W( 88 ) + a*JVS( 293 )
  W( 91 ) = W( 91 ) + a*JVS( 294 )
  W( 94 ) = W( 94 ) + a*JVS( 295 )
  W( 96 ) = W( 96 ) + a*JVS( 296 )
  a = -W( 53 ) / JVS(          297  )
  W( 53 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 298 )
  W( 88 ) = W( 88 ) + a*JVS( 299 )
  W( 91 ) = W( 91 ) + a*JVS( 300 )
  W( 96 ) = W( 96 ) + a*JVS( 301 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 56 ) / JVS(          312  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 313 )
  W( 64 ) = W( 64 ) + a*JVS( 314 )
  W( 85 ) = W( 85 ) + a*JVS( 315 )
  W( 89 ) = W( 89 ) + a*JVS( 316 )
  W( 91 ) = W( 91 ) + a*JVS( 317 )
  W( 94 ) = W( 94 ) + a*JVS( 318 )
  W( 95 ) = W( 95 ) + a*JVS( 319 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 325 )
  W( 88 ) = W( 88 ) + a*JVS( 326 )
  W( 89 ) = W( 89 ) + a*JVS( 327 )
  W( 91 ) = W( 91 ) + a*JVS( 328 )
  W( 95 ) = W( 95 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  a = -W( 58 ) / JVS(          331  )
  W( 58 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 332 )
  W( 85 ) = W( 85 ) + a*JVS( 333 )
  W( 91 ) = W( 91 ) + a*JVS( 334 )
  W( 95 ) = W( 95 ) + a*JVS( 335 )
  W( 96 ) = W( 96 ) + a*JVS( 336 )
  a = -W( 59 ) / JVS(          337  )
  W( 59 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 338 )
  W( 88 ) = W( 88 ) + a*JVS( 339 )
  W( 91 ) = W( 91 ) + a*JVS( 340 )
  W( 96 ) = W( 96 ) + a*JVS( 341 )
  a = -W( 60 ) / JVS(          343  )
  W( 60 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 344 )
  W( 71 ) = W( 71 ) + a*JVS( 345 )
  W( 72 ) = W( 72 ) + a*JVS( 346 )
  W( 79 ) = W( 79 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 82 ) = W( 82 ) + a*JVS( 349 )
  W( 85 ) = W( 85 ) + a*JVS( 350 )
  W( 86 ) = W( 86 ) + a*JVS( 351 )
  W( 88 ) = W( 88 ) + a*JVS( 352 )
  W( 90 ) = W( 90 ) + a*JVS( 353 )
  W( 91 ) = W( 91 ) + a*JVS( 354 )
  W( 93 ) = W( 93 ) + a*JVS( 355 )
  W( 94 ) = W( 94 ) + a*JVS( 356 )
  W( 96 ) = W( 96 ) + a*JVS( 357 )
  a = -W( 61 ) / JVS(          358  )
  W( 61 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 359 )
  W( 89 ) = W( 89 ) + a*JVS( 360 )
  W( 90 ) = W( 90 ) + a*JVS( 361 )
  W( 91 ) = W( 91 ) + a*JVS( 362 )
  W( 93 ) = W( 93 ) + a*JVS( 363 )
  W( 94 ) = W( 94 ) + a*JVS( 364 )
  W( 96 ) = W( 96 ) + a*JVS( 365 )
  a = -W( 62 ) / JVS(          367  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 368 )
  W( 88 ) = W( 88 ) + a*JVS( 369 )
  W( 89 ) = W( 89 ) + a*JVS( 370 )
  W( 90 ) = W( 90 ) + a*JVS( 371 )
  W( 91 ) = W( 91 ) + a*JVS( 372 )
  W( 93 ) = W( 93 ) + a*JVS( 373 )
  W( 94 ) = W( 94 ) + a*JVS( 374 )
  W( 96 ) = W( 96 ) + a*JVS( 375 )
  a = -W( 63 ) / JVS(          377  )
  W( 63 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  W( 90 ) = W( 90 ) + a*JVS( 380 )
  W( 91 ) = W( 91 ) + a*JVS( 381 )
  W( 93 ) = W( 93 ) + a*JVS( 382 )
  W( 94 ) = W( 94 ) + a*JVS( 383 )
  W( 96 ) = W( 96 ) + a*JVS( 384 )
  a = -W( 64 ) / JVS(          391  )
  W( 64 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 392 )
  W( 85 ) = W( 85 ) + a*JVS( 393 )
  W( 88 ) = W( 88 ) + a*JVS( 394 )
  W( 89 ) = W( 89 ) + a*JVS( 395 )
  W( 90 ) = W( 90 ) + a*JVS( 396 )
  W( 91 ) = W( 91 ) + a*JVS( 397 )
  W( 93 ) = W( 93 ) + a*JVS( 398 )
  W( 94 ) = W( 94 ) + a*JVS( 399 )
  W( 96 ) = W( 96 ) + a*JVS( 400 )
  a = -W( 65 ) / JVS(          402  )
  W( 65 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 403 )
  W( 89 ) = W( 89 ) + a*JVS( 404 )
  W( 90 ) = W( 90 ) + a*JVS( 405 )
  W( 91 ) = W( 91 ) + a*JVS( 406 )
  W( 93 ) = W( 93 ) + a*JVS( 407 )
  W( 94 ) = W( 94 ) + a*JVS( 408 )
  W( 96 ) = W( 96 ) + a*JVS( 409 )
  a = -W( 67 ) / JVS(          440  )
  W( 67 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 441 )
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 89 ) = W( 89 ) + a*JVS( 443 )
  W( 90 ) = W( 90 ) + a*JVS( 444 )
  W( 91 ) = W( 91 ) + a*JVS( 445 )
  W( 93 ) = W( 93 ) + a*JVS( 446 )
  W( 94 ) = W( 94 ) + a*JVS( 447 )
  W( 96 ) = W( 96 ) + a*JVS( 448 )
  a = -W( 68 ) / JVS(          451  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 452 )
  W( 88 ) = W( 88 ) + a*JVS( 453 )
  W( 89 ) = W( 89 ) + a*JVS( 454 )
  W( 90 ) = W( 90 ) + a*JVS( 455 )
  W( 91 ) = W( 91 ) + a*JVS( 456 )
  W( 93 ) = W( 93 ) + a*JVS( 457 )
  W( 94 ) = W( 94 ) + a*JVS( 458 )
  W( 95 ) = W( 95 ) + a*JVS( 459 )
  W( 96 ) = W( 96 ) + a*JVS( 460 )
  a = -W( 69 ) / JVS(          471  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 472 )
  W( 72 ) = W( 72 ) + a*JVS( 473 )
  W( 74 ) = W( 74 ) + a*JVS( 474 )
  W( 79 ) = W( 79 ) + a*JVS( 475 )
  W( 81 ) = W( 81 ) + a*JVS( 476 )
  W( 82 ) = W( 82 ) + a*JVS( 477 )
  W( 85 ) = W( 85 ) + a*JVS( 478 )
  W( 86 ) = W( 86 ) + a*JVS( 479 )
  W( 87 ) = W( 87 ) + a*JVS( 480 )
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 89 ) = W( 89 ) + a*JVS( 482 )
  W( 90 ) = W( 90 ) + a*JVS( 483 )
  W( 91 ) = W( 91 ) + a*JVS( 484 )
  W( 93 ) = W( 93 ) + a*JVS( 485 )
  W( 94 ) = W( 94 ) + a*JVS( 486 )
  W( 95 ) = W( 95 ) + a*JVS( 487 )
  W( 96 ) = W( 96 ) + a*JVS( 488 )
  a = -W( 70 ) / JVS(          491  )
  W( 70 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 492 )
  W( 85 ) = W( 85 ) + a*JVS( 493 )
  W( 89 ) = W( 89 ) + a*JVS( 494 )
  W( 90 ) = W( 90 ) + a*JVS( 495 )
  W( 91 ) = W( 91 ) + a*JVS( 496 )
  W( 93 ) = W( 93 ) + a*JVS( 497 )
  W( 94 ) = W( 94 ) + a*JVS( 498 )
  W( 96 ) = W( 96 ) + a*JVS( 499 )
  a = -W( 71 ) / JVS(          501  )
  W( 71 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 502 )
  W( 88 ) = W( 88 ) + a*JVS( 503 )
  W( 89 ) = W( 89 ) + a*JVS( 504 )
  W( 90 ) = W( 90 ) + a*JVS( 505 )
  W( 91 ) = W( 91 ) + a*JVS( 506 )
  W( 93 ) = W( 93 ) + a*JVS( 507 )
  W( 94 ) = W( 94 ) + a*JVS( 508 )
  W( 95 ) = W( 95 ) + a*JVS( 509 )
  W( 96 ) = W( 96 ) + a*JVS( 510 )
  a = -W( 72 ) / JVS(          512  )
  W( 72 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 513 )
  W( 88 ) = W( 88 ) + a*JVS( 514 )
  W( 89 ) = W( 89 ) + a*JVS( 515 )
  W( 90 ) = W( 90 ) + a*JVS( 516 )
  W( 91 ) = W( 91 ) + a*JVS( 517 )
  W( 93 ) = W( 93 ) + a*JVS( 518 )
  W( 94 ) = W( 94 ) + a*JVS( 519 )
  W( 95 ) = W( 95 ) + a*JVS( 520 )
  W( 96 ) = W( 96 ) + a*JVS( 521 )
  a = -W( 73 ) / JVS(          549  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 550 )
  W( 76 ) = W( 76 ) + a*JVS( 551 )
  W( 77 ) = W( 77 ) + a*JVS( 552 )
  W( 79 ) = W( 79 ) + a*JVS( 553 )
  W( 80 ) = W( 80 ) + a*JVS( 554 )
  W( 81 ) = W( 81 ) + a*JVS( 555 )
  W( 82 ) = W( 82 ) + a*JVS( 556 )
  W( 83 ) = W( 83 ) + a*JVS( 557 )
  W( 84 ) = W( 84 ) + a*JVS( 558 )
  W( 85 ) = W( 85 ) + a*JVS( 559 )
  W( 86 ) = W( 86 ) + a*JVS( 560 )
  W( 87 ) = W( 87 ) + a*JVS( 561 )
  W( 88 ) = W( 88 ) + a*JVS( 562 )
  W( 89 ) = W( 89 ) + a*JVS( 563 )
  W( 90 ) = W( 90 ) + a*JVS( 564 )
  W( 91 ) = W( 91 ) + a*JVS( 565 )
  W( 93 ) = W( 93 ) + a*JVS( 566 )
  W( 94 ) = W( 94 ) + a*JVS( 567 )
  W( 95 ) = W( 95 ) + a*JVS( 568 )
  W( 96 ) = W( 96 ) + a*JVS( 569 )
  W( 97 ) = W( 97 ) + a*JVS( 570 )
  a = -W( 74 ) / JVS(          572  )
  W( 74 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 573 )
  W( 85 ) = W( 85 ) + a*JVS( 574 )
  W( 89 ) = W( 89 ) + a*JVS( 575 )
  W( 90 ) = W( 90 ) + a*JVS( 576 )
  W( 91 ) = W( 91 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  a = -W( 75 ) / JVS(          590  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 591 )
  W( 77 ) = W( 77 ) + a*JVS( 592 )
  W( 79 ) = W( 79 ) + a*JVS( 593 )
  W( 83 ) = W( 83 ) + a*JVS( 594 )
  W( 84 ) = W( 84 ) + a*JVS( 595 )
  W( 85 ) = W( 85 ) + a*JVS( 596 )
  W( 86 ) = W( 86 ) + a*JVS( 597 )
  W( 88 ) = W( 88 ) + a*JVS( 598 )
  W( 89 ) = W( 89 ) + a*JVS( 599 )
  W( 90 ) = W( 90 ) + a*JVS( 600 )
  W( 91 ) = W( 91 ) + a*JVS( 601 )
  W( 93 ) = W( 93 ) + a*JVS( 602 )
  W( 94 ) = W( 94 ) + a*JVS( 603 )
  W( 95 ) = W( 95 ) + a*JVS( 604 )
  W( 96 ) = W( 96 ) + a*JVS( 605 )
  a = -W( 76 ) / JVS(          607  )
  W( 76 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 608 )
  W( 82 ) = W( 82 ) + a*JVS( 609 )
  W( 88 ) = W( 88 ) + a*JVS( 610 )
  W( 89 ) = W( 89 ) + a*JVS( 611 )
  W( 90 ) = W( 90 ) + a*JVS( 612 )
  W( 91 ) = W( 91 ) + a*JVS( 613 )
  W( 93 ) = W( 93 ) + a*JVS( 614 )
  W( 94 ) = W( 94 ) + a*JVS( 615 )
  W( 96 ) = W( 96 ) + a*JVS( 616 )
  a = -W( 77 ) / JVS(          618  )
  W( 77 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 619 )
  W( 82 ) = W( 82 ) + a*JVS( 620 )
  W( 88 ) = W( 88 ) + a*JVS( 621 )
  W( 89 ) = W( 89 ) + a*JVS( 622 )
  W( 90 ) = W( 90 ) + a*JVS( 623 )
  W( 91 ) = W( 91 ) + a*JVS( 624 )
  W( 93 ) = W( 93 ) + a*JVS( 625 )
  W( 94 ) = W( 94 ) + a*JVS( 626 )
  W( 96 ) = W( 96 ) + a*JVS( 627 )
  a = -W( 78 ) / JVS(          643  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 644 )
  W( 81 ) = W( 81 ) + a*JVS( 645 )
  W( 82 ) = W( 82 ) + a*JVS( 646 )
  W( 83 ) = W( 83 ) + a*JVS( 647 )
  W( 84 ) = W( 84 ) + a*JVS( 648 )
  W( 85 ) = W( 85 ) + a*JVS( 649 )
  W( 86 ) = W( 86 ) + a*JVS( 650 )
  W( 87 ) = W( 87 ) + a*JVS( 651 )
  W( 88 ) = W( 88 ) + a*JVS( 652 )
  W( 89 ) = W( 89 ) + a*JVS( 653 )
  W( 90 ) = W( 90 ) + a*JVS( 654 )
  W( 91 ) = W( 91 ) + a*JVS( 655 )
  W( 92 ) = W( 92 ) + a*JVS( 656 )
  W( 93 ) = W( 93 ) + a*JVS( 657 )
  W( 94 ) = W( 94 ) + a*JVS( 658 )
  W( 95 ) = W( 95 ) + a*JVS( 659 )
  W( 96 ) = W( 96 ) + a*JVS( 660 )
  W( 97 ) = W( 97 ) + a*JVS( 661 )
  a = -W( 79 ) / JVS(          664  )
  W( 79 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 665 )
  W( 88 ) = W( 88 ) + a*JVS( 666 )
  W( 89 ) = W( 89 ) + a*JVS( 667 )
  W( 90 ) = W( 90 ) + a*JVS( 668 )
  W( 91 ) = W( 91 ) + a*JVS( 669 )
  W( 93 ) = W( 93 ) + a*JVS( 670 )
  W( 94 ) = W( 94 ) + a*JVS( 671 )
  W( 95 ) = W( 95 ) + a*JVS( 672 )
  W( 96 ) = W( 96 ) + a*JVS( 673 )
  a = -W( 80 ) / JVS(          692  )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 693 )
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 85 ) = W( 85 ) + a*JVS( 695 )
  W( 86 ) = W( 86 ) + a*JVS( 696 )
  W( 87 ) = W( 87 ) + a*JVS( 697 )
  W( 88 ) = W( 88 ) + a*JVS( 698 )
  W( 89 ) = W( 89 ) + a*JVS( 699 )
  W( 90 ) = W( 90 ) + a*JVS( 700 )
  W( 91 ) = W( 91 ) + a*JVS( 701 )
  W( 92 ) = W( 92 ) + a*JVS( 702 )
  W( 93 ) = W( 93 ) + a*JVS( 703 )
  W( 94 ) = W( 94 ) + a*JVS( 704 )
  W( 95 ) = W( 95 ) + a*JVS( 705 )
  W( 96 ) = W( 96 ) + a*JVS( 706 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  a = -W( 83 ) / JVS(          739  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 740 )
  W( 85 ) = W( 85 ) + a*JVS( 741 )
  W( 88 ) = W( 88 ) + a*JVS( 742 )
  W( 89 ) = W( 89 ) + a*JVS( 743 )
  W( 90 ) = W( 90 ) + a*JVS( 744 )
  W( 91 ) = W( 91 ) + a*JVS( 745 )
  W( 93 ) = W( 93 ) + a*JVS( 746 )
  W( 94 ) = W( 94 ) + a*JVS( 747 )
  W( 96 ) = W( 96 ) + a*JVS( 748 )
  a = -W( 84 ) / JVS(          758  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 759 )
  W( 88 ) = W( 88 ) + a*JVS( 760 )
  W( 89 ) = W( 89 ) + a*JVS( 761 )
  W( 90 ) = W( 90 ) + a*JVS( 762 )
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 93 ) = W( 93 ) + a*JVS( 764 )
  W( 94 ) = W( 94 ) + a*JVS( 765 )
  W( 96 ) = W( 96 ) + a*JVS( 766 )
  a = -W( 85 ) / JVS(          786  )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 787 )
  W( 89 ) = W( 89 ) + a*JVS( 788 )
  W( 90 ) = W( 90 ) + a*JVS( 789 )
  W( 91 ) = W( 91 ) + a*JVS( 790 )
  W( 93 ) = W( 93 ) + a*JVS( 791 )
  W( 94 ) = W( 94 ) + a*JVS( 792 )
  W( 95 ) = W( 95 ) + a*JVS( 793 )
  W( 96 ) = W( 96 ) + a*JVS( 794 )
  a = -W( 86 ) / JVS(          805  )
  W( 86 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 806 )
  W( 89 ) = W( 89 ) + a*JVS( 807 )
  W( 90 ) = W( 90 ) + a*JVS( 808 )
  W( 91 ) = W( 91 ) + a*JVS( 809 )
  W( 92 ) = W( 92 ) + a*JVS( 810 )
  W( 93 ) = W( 93 ) + a*JVS( 811 )
  W( 94 ) = W( 94 ) + a*JVS( 812 )
  W( 95 ) = W( 95 ) + a*JVS( 813 )
  W( 96 ) = W( 96 ) + a*JVS( 814 )
  a = -W( 87 ) / JVS(          831  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 832 )
  W( 89 ) = W( 89 ) + a*JVS( 833 )
  W( 90 ) = W( 90 ) + a*JVS( 834 )
  W( 91 ) = W( 91 ) + a*JVS( 835 )
  W( 92 ) = W( 92 ) + a*JVS( 836 )
  W( 93 ) = W( 93 ) + a*JVS( 837 )
  W( 94 ) = W( 94 ) + a*JVS( 838 )
  W( 95 ) = W( 95 ) + a*JVS( 839 )
  W( 96 ) = W( 96 ) + a*JVS( 840 )
  a = -W( 88 ) / JVS(          871  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 872 )
  W( 90 ) = W( 90 ) + a*JVS( 873 )
  W( 91 ) = W( 91 ) + a*JVS( 874 )
  W( 92 ) = W( 92 ) + a*JVS( 875 )
  W( 93 ) = W( 93 ) + a*JVS( 876 )
  W( 94 ) = W( 94 ) + a*JVS( 877 )
  W( 95 ) = W( 95 ) + a*JVS( 878 )
  W( 96 ) = W( 96 ) + a*JVS( 879 )
  W( 97 ) = W( 97 ) + a*JVS( 880 )
  JVS( 881) = W( 15 )
  JVS( 882) = W( 21 )
  JVS( 883) = W( 23 )
  JVS( 884) = W( 24 )
  JVS( 885) = W( 25 )
  JVS( 886) = W( 27 )
  JVS( 887) = W( 28 )
  JVS( 888) = W( 31 )
  JVS( 889) = W( 32 )
  JVS( 890) = W( 33 )
  JVS( 891) = W( 36 )
  JVS( 892) = W( 37 )
  JVS( 893) = W( 38 )
  JVS( 894) = W( 39 )
  JVS( 895) = W( 40 )
  JVS( 896) = W( 41 )
  JVS( 897) = W( 42 )
  JVS( 898) = W( 45 )
  JVS( 899) = W( 46 )
  JVS( 900) = W( 47 )
  JVS( 901) = W( 49 )
  JVS( 902) = W( 50 )
  JVS( 903) = W( 51 )
  JVS( 904) = W( 52 )
  JVS( 905) = W( 53 )
  JVS( 906) = W( 54 )
  JVS( 907) = W( 55 )
  JVS( 908) = W( 56 )
  JVS( 909) = W( 57 )
  JVS( 910) = W( 58 )
  JVS( 911) = W( 59 )
  JVS( 912) = W( 60 )
  JVS( 913) = W( 61 )
  JVS( 914) = W( 62 )
  JVS( 915) = W( 63 )
  JVS( 916) = W( 64 )
  JVS( 917) = W( 65 )
  JVS( 918) = W( 67 )
  JVS( 919) = W( 68 )
  JVS( 920) = W( 69 )
  JVS( 921) = W( 70 )
  JVS( 922) = W( 71 )
  JVS( 923) = W( 72 )
  JVS( 924) = W( 73 )
  JVS( 925) = W( 74 )
  JVS( 926) = W( 75 )
  JVS( 927) = W( 76 )
  JVS( 928) = W( 77 )
  JVS( 929) = W( 78 )
  JVS( 930) = W( 79 )
  JVS( 931) = W( 80 )
  JVS( 932) = W( 81 )
  JVS( 933) = W( 82 )
  JVS( 934) = W( 83 )
  JVS( 935) = W( 84 )
  JVS( 936) = W( 85 )
  JVS( 937) = W( 86 )
  JVS( 938) = W( 87 )
  JVS( 939) = W( 88 )
  JVS( 940) = W( 89 )
  JVS( 941) = W( 90 )
  JVS( 942) = W( 91 )
  JVS( 943) = W( 92 )
  JVS( 944) = W( 93 )
  JVS( 945) = W( 94 )
  JVS( 946) = W( 95 )
  JVS( 947) = W( 96 )
  JVS( 948) = W( 97 )
  IF ( ABS(  JVS( 979 )) < TINY(a) ) THEN
         IER = 90                                      
         RETURN
  END IF
   W( 25 ) = JVS( 949 )
   W( 35 ) = JVS( 950 )
   W( 45 ) = JVS( 951 )
   W( 54 ) = JVS( 952 )
   W( 55 ) = JVS( 953 )
   W( 59 ) = JVS( 954 )
   W( 61 ) = JVS( 955 )
   W( 62 ) = JVS( 956 )
   W( 63 ) = JVS( 957 )
   W( 65 ) = JVS( 958 )
   W( 67 ) = JVS( 959 )
   W( 68 ) = JVS( 960 )
   W( 70 ) = JVS( 961 )
   W( 71 ) = JVS( 962 )
   W( 72 ) = JVS( 963 )
   W( 74 ) = JVS( 964 )
   W( 76 ) = JVS( 965 )
   W( 77 ) = JVS( 966 )
   W( 78 ) = JVS( 967 )
   W( 79 ) = JVS( 968 )
   W( 80 ) = JVS( 969 )
   W( 81 ) = JVS( 970 )
   W( 82 ) = JVS( 971 )
   W( 83 ) = JVS( 972 )
   W( 84 ) = JVS( 973 )
   W( 85 ) = JVS( 974 )
   W( 86 ) = JVS( 975 )
   W( 87 ) = JVS( 976 )
   W( 88 ) = JVS( 977 )
   W( 89 ) = JVS( 978 )
   W( 90 ) = JVS( 979 )
   W( 91 ) = JVS( 980 )
   W( 92 ) = JVS( 981 )
   W( 93 ) = JVS( 982 )
   W( 94 ) = JVS( 983 )
   W( 95 ) = JVS( 984 )
   W( 96 ) = JVS( 985 )
   W( 97 ) = JVS( 986 )
  a = -W( 25 ) / JVS(          100  )
  W( 25 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 101 )
  W( 90 ) = W( 90 ) + a*JVS( 102 )
  W( 91 ) = W( 91 ) + a*JVS( 103 )
  a = -W( 35 ) / JVS(          169  )
  W( 35 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 170 )
  W( 81 ) = W( 81 ) + a*JVS( 171 )
  W( 82 ) = W( 82 ) + a*JVS( 172 )
  W( 85 ) = W( 85 ) + a*JVS( 173 )
  W( 88 ) = W( 88 ) + a*JVS( 174 )
  W( 91 ) = W( 91 ) + a*JVS( 175 )
  a = -W( 45 ) / JVS(          232  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 233 )
  W( 85 ) = W( 85 ) + a*JVS( 234 )
  W( 89 ) = W( 89 ) + a*JVS( 235 )
  W( 91 ) = W( 91 ) + a*JVS( 236 )
  W( 93 ) = W( 93 ) + a*JVS( 237 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 59 ) / JVS(          337  )
  W( 59 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 338 )
  W( 88 ) = W( 88 ) + a*JVS( 339 )
  W( 91 ) = W( 91 ) + a*JVS( 340 )
  W( 96 ) = W( 96 ) + a*JVS( 341 )
  a = -W( 61 ) / JVS(          358  )
  W( 61 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 359 )
  W( 89 ) = W( 89 ) + a*JVS( 360 )
  W( 90 ) = W( 90 ) + a*JVS( 361 )
  W( 91 ) = W( 91 ) + a*JVS( 362 )
  W( 93 ) = W( 93 ) + a*JVS( 363 )
  W( 94 ) = W( 94 ) + a*JVS( 364 )
  W( 96 ) = W( 96 ) + a*JVS( 365 )
  a = -W( 62 ) / JVS(          367  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 368 )
  W( 88 ) = W( 88 ) + a*JVS( 369 )
  W( 89 ) = W( 89 ) + a*JVS( 370 )
  W( 90 ) = W( 90 ) + a*JVS( 371 )
  W( 91 ) = W( 91 ) + a*JVS( 372 )
  W( 93 ) = W( 93 ) + a*JVS( 373 )
  W( 94 ) = W( 94 ) + a*JVS( 374 )
  W( 96 ) = W( 96 ) + a*JVS( 375 )
  a = -W( 63 ) / JVS(          377  )
  W( 63 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  W( 90 ) = W( 90 ) + a*JVS( 380 )
  W( 91 ) = W( 91 ) + a*JVS( 381 )
  W( 93 ) = W( 93 ) + a*JVS( 382 )
  W( 94 ) = W( 94 ) + a*JVS( 383 )
  W( 96 ) = W( 96 ) + a*JVS( 384 )
  a = -W( 65 ) / JVS(          402  )
  W( 65 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 403 )
  W( 89 ) = W( 89 ) + a*JVS( 404 )
  W( 90 ) = W( 90 ) + a*JVS( 405 )
  W( 91 ) = W( 91 ) + a*JVS( 406 )
  W( 93 ) = W( 93 ) + a*JVS( 407 )
  W( 94 ) = W( 94 ) + a*JVS( 408 )
  W( 96 ) = W( 96 ) + a*JVS( 409 )
  a = -W( 67 ) / JVS(          440  )
  W( 67 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 441 )
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 89 ) = W( 89 ) + a*JVS( 443 )
  W( 90 ) = W( 90 ) + a*JVS( 444 )
  W( 91 ) = W( 91 ) + a*JVS( 445 )
  W( 93 ) = W( 93 ) + a*JVS( 446 )
  W( 94 ) = W( 94 ) + a*JVS( 447 )
  W( 96 ) = W( 96 ) + a*JVS( 448 )
  a = -W( 68 ) / JVS(          451  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 452 )
  W( 88 ) = W( 88 ) + a*JVS( 453 )
  W( 89 ) = W( 89 ) + a*JVS( 454 )
  W( 90 ) = W( 90 ) + a*JVS( 455 )
  W( 91 ) = W( 91 ) + a*JVS( 456 )
  W( 93 ) = W( 93 ) + a*JVS( 457 )
  W( 94 ) = W( 94 ) + a*JVS( 458 )
  W( 95 ) = W( 95 ) + a*JVS( 459 )
  W( 96 ) = W( 96 ) + a*JVS( 460 )
  a = -W( 70 ) / JVS(          491  )
  W( 70 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 492 )
  W( 85 ) = W( 85 ) + a*JVS( 493 )
  W( 89 ) = W( 89 ) + a*JVS( 494 )
  W( 90 ) = W( 90 ) + a*JVS( 495 )
  W( 91 ) = W( 91 ) + a*JVS( 496 )
  W( 93 ) = W( 93 ) + a*JVS( 497 )
  W( 94 ) = W( 94 ) + a*JVS( 498 )
  W( 96 ) = W( 96 ) + a*JVS( 499 )
  a = -W( 71 ) / JVS(          501  )
  W( 71 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 502 )
  W( 88 ) = W( 88 ) + a*JVS( 503 )
  W( 89 ) = W( 89 ) + a*JVS( 504 )
  W( 90 ) = W( 90 ) + a*JVS( 505 )
  W( 91 ) = W( 91 ) + a*JVS( 506 )
  W( 93 ) = W( 93 ) + a*JVS( 507 )
  W( 94 ) = W( 94 ) + a*JVS( 508 )
  W( 95 ) = W( 95 ) + a*JVS( 509 )
  W( 96 ) = W( 96 ) + a*JVS( 510 )
  a = -W( 72 ) / JVS(          512  )
  W( 72 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 513 )
  W( 88 ) = W( 88 ) + a*JVS( 514 )
  W( 89 ) = W( 89 ) + a*JVS( 515 )
  W( 90 ) = W( 90 ) + a*JVS( 516 )
  W( 91 ) = W( 91 ) + a*JVS( 517 )
  W( 93 ) = W( 93 ) + a*JVS( 518 )
  W( 94 ) = W( 94 ) + a*JVS( 519 )
  W( 95 ) = W( 95 ) + a*JVS( 520 )
  W( 96 ) = W( 96 ) + a*JVS( 521 )
  a = -W( 74 ) / JVS(          572  )
  W( 74 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 573 )
  W( 85 ) = W( 85 ) + a*JVS( 574 )
  W( 89 ) = W( 89 ) + a*JVS( 575 )
  W( 90 ) = W( 90 ) + a*JVS( 576 )
  W( 91 ) = W( 91 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  a = -W( 76 ) / JVS(          607  )
  W( 76 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 608 )
  W( 82 ) = W( 82 ) + a*JVS( 609 )
  W( 88 ) = W( 88 ) + a*JVS( 610 )
  W( 89 ) = W( 89 ) + a*JVS( 611 )
  W( 90 ) = W( 90 ) + a*JVS( 612 )
  W( 91 ) = W( 91 ) + a*JVS( 613 )
  W( 93 ) = W( 93 ) + a*JVS( 614 )
  W( 94 ) = W( 94 ) + a*JVS( 615 )
  W( 96 ) = W( 96 ) + a*JVS( 616 )
  a = -W( 77 ) / JVS(          618  )
  W( 77 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 619 )
  W( 82 ) = W( 82 ) + a*JVS( 620 )
  W( 88 ) = W( 88 ) + a*JVS( 621 )
  W( 89 ) = W( 89 ) + a*JVS( 622 )
  W( 90 ) = W( 90 ) + a*JVS( 623 )
  W( 91 ) = W( 91 ) + a*JVS( 624 )
  W( 93 ) = W( 93 ) + a*JVS( 625 )
  W( 94 ) = W( 94 ) + a*JVS( 626 )
  W( 96 ) = W( 96 ) + a*JVS( 627 )
  a = -W( 78 ) / JVS(          643  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 644 )
  W( 81 ) = W( 81 ) + a*JVS( 645 )
  W( 82 ) = W( 82 ) + a*JVS( 646 )
  W( 83 ) = W( 83 ) + a*JVS( 647 )
  W( 84 ) = W( 84 ) + a*JVS( 648 )
  W( 85 ) = W( 85 ) + a*JVS( 649 )
  W( 86 ) = W( 86 ) + a*JVS( 650 )
  W( 87 ) = W( 87 ) + a*JVS( 651 )
  W( 88 ) = W( 88 ) + a*JVS( 652 )
  W( 89 ) = W( 89 ) + a*JVS( 653 )
  W( 90 ) = W( 90 ) + a*JVS( 654 )
  W( 91 ) = W( 91 ) + a*JVS( 655 )
  W( 92 ) = W( 92 ) + a*JVS( 656 )
  W( 93 ) = W( 93 ) + a*JVS( 657 )
  W( 94 ) = W( 94 ) + a*JVS( 658 )
  W( 95 ) = W( 95 ) + a*JVS( 659 )
  W( 96 ) = W( 96 ) + a*JVS( 660 )
  W( 97 ) = W( 97 ) + a*JVS( 661 )
  a = -W( 79 ) / JVS(          664  )
  W( 79 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 665 )
  W( 88 ) = W( 88 ) + a*JVS( 666 )
  W( 89 ) = W( 89 ) + a*JVS( 667 )
  W( 90 ) = W( 90 ) + a*JVS( 668 )
  W( 91 ) = W( 91 ) + a*JVS( 669 )
  W( 93 ) = W( 93 ) + a*JVS( 670 )
  W( 94 ) = W( 94 ) + a*JVS( 671 )
  W( 95 ) = W( 95 ) + a*JVS( 672 )
  W( 96 ) = W( 96 ) + a*JVS( 673 )
  a = -W( 80 ) / JVS(          692  )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 693 )
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 85 ) = W( 85 ) + a*JVS( 695 )
  W( 86 ) = W( 86 ) + a*JVS( 696 )
  W( 87 ) = W( 87 ) + a*JVS( 697 )
  W( 88 ) = W( 88 ) + a*JVS( 698 )
  W( 89 ) = W( 89 ) + a*JVS( 699 )
  W( 90 ) = W( 90 ) + a*JVS( 700 )
  W( 91 ) = W( 91 ) + a*JVS( 701 )
  W( 92 ) = W( 92 ) + a*JVS( 702 )
  W( 93 ) = W( 93 ) + a*JVS( 703 )
  W( 94 ) = W( 94 ) + a*JVS( 704 )
  W( 95 ) = W( 95 ) + a*JVS( 705 )
  W( 96 ) = W( 96 ) + a*JVS( 706 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  a = -W( 83 ) / JVS(          739  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 740 )
  W( 85 ) = W( 85 ) + a*JVS( 741 )
  W( 88 ) = W( 88 ) + a*JVS( 742 )
  W( 89 ) = W( 89 ) + a*JVS( 743 )
  W( 90 ) = W( 90 ) + a*JVS( 744 )
  W( 91 ) = W( 91 ) + a*JVS( 745 )
  W( 93 ) = W( 93 ) + a*JVS( 746 )
  W( 94 ) = W( 94 ) + a*JVS( 747 )
  W( 96 ) = W( 96 ) + a*JVS( 748 )
  a = -W( 84 ) / JVS(          758  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 759 )
  W( 88 ) = W( 88 ) + a*JVS( 760 )
  W( 89 ) = W( 89 ) + a*JVS( 761 )
  W( 90 ) = W( 90 ) + a*JVS( 762 )
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 93 ) = W( 93 ) + a*JVS( 764 )
  W( 94 ) = W( 94 ) + a*JVS( 765 )
  W( 96 ) = W( 96 ) + a*JVS( 766 )
  a = -W( 85 ) / JVS(          786  )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 787 )
  W( 89 ) = W( 89 ) + a*JVS( 788 )
  W( 90 ) = W( 90 ) + a*JVS( 789 )
  W( 91 ) = W( 91 ) + a*JVS( 790 )
  W( 93 ) = W( 93 ) + a*JVS( 791 )
  W( 94 ) = W( 94 ) + a*JVS( 792 )
  W( 95 ) = W( 95 ) + a*JVS( 793 )
  W( 96 ) = W( 96 ) + a*JVS( 794 )
  a = -W( 86 ) / JVS(          805  )
  W( 86 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 806 )
  W( 89 ) = W( 89 ) + a*JVS( 807 )
  W( 90 ) = W( 90 ) + a*JVS( 808 )
  W( 91 ) = W( 91 ) + a*JVS( 809 )
  W( 92 ) = W( 92 ) + a*JVS( 810 )
  W( 93 ) = W( 93 ) + a*JVS( 811 )
  W( 94 ) = W( 94 ) + a*JVS( 812 )
  W( 95 ) = W( 95 ) + a*JVS( 813 )
  W( 96 ) = W( 96 ) + a*JVS( 814 )
  a = -W( 87 ) / JVS(          831  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 832 )
  W( 89 ) = W( 89 ) + a*JVS( 833 )
  W( 90 ) = W( 90 ) + a*JVS( 834 )
  W( 91 ) = W( 91 ) + a*JVS( 835 )
  W( 92 ) = W( 92 ) + a*JVS( 836 )
  W( 93 ) = W( 93 ) + a*JVS( 837 )
  W( 94 ) = W( 94 ) + a*JVS( 838 )
  W( 95 ) = W( 95 ) + a*JVS( 839 )
  W( 96 ) = W( 96 ) + a*JVS( 840 )
  a = -W( 88 ) / JVS(          871  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 872 )
  W( 90 ) = W( 90 ) + a*JVS( 873 )
  W( 91 ) = W( 91 ) + a*JVS( 874 )
  W( 92 ) = W( 92 ) + a*JVS( 875 )
  W( 93 ) = W( 93 ) + a*JVS( 876 )
  W( 94 ) = W( 94 ) + a*JVS( 877 )
  W( 95 ) = W( 95 ) + a*JVS( 878 )
  W( 96 ) = W( 96 ) + a*JVS( 879 )
  W( 97 ) = W( 97 ) + a*JVS( 880 )
  a = -W( 89 ) / JVS(          940  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 941 )
  W( 91 ) = W( 91 ) + a*JVS( 942 )
  W( 92 ) = W( 92 ) + a*JVS( 943 )
  W( 93 ) = W( 93 ) + a*JVS( 944 )
  W( 94 ) = W( 94 ) + a*JVS( 945 )
  W( 95 ) = W( 95 ) + a*JVS( 946 )
  W( 96 ) = W( 96 ) + a*JVS( 947 )
  W( 97 ) = W( 97 ) + a*JVS( 948 )
  JVS( 949) = W( 25 )
  JVS( 950) = W( 35 )
  JVS( 951) = W( 45 )
  JVS( 952) = W( 54 )
  JVS( 953) = W( 55 )
  JVS( 954) = W( 59 )
  JVS( 955) = W( 61 )
  JVS( 956) = W( 62 )
  JVS( 957) = W( 63 )
  JVS( 958) = W( 65 )
  JVS( 959) = W( 67 )
  JVS( 960) = W( 68 )
  JVS( 961) = W( 70 )
  JVS( 962) = W( 71 )
  JVS( 963) = W( 72 )
  JVS( 964) = W( 74 )
  JVS( 965) = W( 76 )
  JVS( 966) = W( 77 )
  JVS( 967) = W( 78 )
  JVS( 968) = W( 79 )
  JVS( 969) = W( 80 )
  JVS( 970) = W( 81 )
  JVS( 971) = W( 82 )
  JVS( 972) = W( 83 )
  JVS( 973) = W( 84 )
  JVS( 974) = W( 85 )
  JVS( 975) = W( 86 )
  JVS( 976) = W( 87 )
  JVS( 977) = W( 88 )
  JVS( 978) = W( 89 )
  JVS( 979) = W( 90 )
  JVS( 980) = W( 91 )
  JVS( 981) = W( 92 )
  JVS( 982) = W( 93 )
  JVS( 983) = W( 94 )
  JVS( 984) = W( 95 )
  JVS( 985) = W( 96 )
  JVS( 986) = W( 97 )
  IF ( ABS(  JVS( 1055 )) < TINY(a) ) THEN
         IER = 91                                      
         RETURN
  END IF
   W( 14 ) = JVS( 987 )
   W( 15 ) = JVS( 988 )
   W( 19 ) = JVS( 989 )
   W( 20 ) = JVS( 990 )
   W( 21 ) = JVS( 991 )
   W( 22 ) = JVS( 992 )
   W( 23 ) = JVS( 993 )
   W( 24 ) = JVS( 994 )
   W( 25 ) = JVS( 995 )
   W( 26 ) = JVS( 996 )
   W( 27 ) = JVS( 997 )
   W( 28 ) = JVS( 998 )
   W( 29 ) = JVS( 999 )
   W( 32 ) = JVS( 1000 )
   W( 33 ) = JVS( 1001 )
   W( 34 ) = JVS( 1002 )
   W( 35 ) = JVS( 1003 )
   W( 36 ) = JVS( 1004 )
   W( 37 ) = JVS( 1005 )
   W( 38 ) = JVS( 1006 )
   W( 39 ) = JVS( 1007 )
   W( 40 ) = JVS( 1008 )
   W( 41 ) = JVS( 1009 )
   W( 42 ) = JVS( 1010 )
   W( 43 ) = JVS( 1011 )
   W( 44 ) = JVS( 1012 )
   W( 45 ) = JVS( 1013 )
   W( 46 ) = JVS( 1014 )
   W( 47 ) = JVS( 1015 )
   W( 48 ) = JVS( 1016 )
   W( 49 ) = JVS( 1017 )
   W( 50 ) = JVS( 1018 )
   W( 51 ) = JVS( 1019 )
   W( 52 ) = JVS( 1020 )
   W( 53 ) = JVS( 1021 )
   W( 54 ) = JVS( 1022 )
   W( 55 ) = JVS( 1023 )
   W( 56 ) = JVS( 1024 )
   W( 57 ) = JVS( 1025 )
   W( 58 ) = JVS( 1026 )
   W( 59 ) = JVS( 1027 )
   W( 60 ) = JVS( 1028 )
   W( 64 ) = JVS( 1029 )
   W( 66 ) = JVS( 1030 )
   W( 67 ) = JVS( 1031 )
   W( 68 ) = JVS( 1032 )
   W( 69 ) = JVS( 1033 )
   W( 70 ) = JVS( 1034 )
   W( 71 ) = JVS( 1035 )
   W( 72 ) = JVS( 1036 )
   W( 73 ) = JVS( 1037 )
   W( 74 ) = JVS( 1038 )
   W( 75 ) = JVS( 1039 )
   W( 76 ) = JVS( 1040 )
   W( 77 ) = JVS( 1041 )
   W( 78 ) = JVS( 1042 )
   W( 79 ) = JVS( 1043 )
   W( 80 ) = JVS( 1044 )
   W( 81 ) = JVS( 1045 )
   W( 82 ) = JVS( 1046 )
   W( 83 ) = JVS( 1047 )
   W( 84 ) = JVS( 1048 )
   W( 85 ) = JVS( 1049 )
   W( 86 ) = JVS( 1050 )
   W( 87 ) = JVS( 1051 )
   W( 88 ) = JVS( 1052 )
   W( 89 ) = JVS( 1053 )
   W( 90 ) = JVS( 1054 )
   W( 91 ) = JVS( 1055 )
   W( 92 ) = JVS( 1056 )
   W( 93 ) = JVS( 1057 )
   W( 94 ) = JVS( 1058 )
   W( 95 ) = JVS( 1059 )
   W( 96 ) = JVS( 1060 )
   W( 97 ) = JVS( 1061 )
  a = -W( 14 ) / JVS(           67  )
  W( 14 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 68 )
  a = -W( 15 ) / JVS(           69  )
  W( 15 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 70 )
  a = -W( 19 ) / JVS(           81  )
  W( 19 ) = -a
  W( 32 ) = W( 32 ) + a*JVS( 82 )
  W( 89 ) = W( 89 ) + a*JVS( 83 )
  a = -W( 20 ) / JVS(           84  )
  W( 20 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 85 )
  W( 89 ) = W( 89 ) + a*JVS( 86 )
  W( 91 ) = W( 91 ) + a*JVS( 87 )
  a = -W( 21 ) / JVS(           88  )
  W( 21 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 89 )
  W( 91 ) = W( 91 ) + a*JVS( 90 )
  a = -W( 22 ) / JVS(           91  )
  W( 22 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 92 )
  W( 91 ) = W( 91 ) + a*JVS( 93 )
  a = -W( 23 ) / JVS(           94  )
  W( 23 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 95 )
  W( 91 ) = W( 91 ) + a*JVS( 96 )
  a = -W( 24 ) / JVS(           97  )
  W( 24 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 98 )
  W( 91 ) = W( 91 ) + a*JVS( 99 )
  a = -W( 25 ) / JVS(          100  )
  W( 25 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 101 )
  W( 90 ) = W( 90 ) + a*JVS( 102 )
  W( 91 ) = W( 91 ) + a*JVS( 103 )
  a = -W( 26 ) / JVS(          104  )
  W( 26 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 105 )
  W( 91 ) = W( 91 ) + a*JVS( 106 )
  W( 95 ) = W( 95 ) + a*JVS( 107 )
  a = -W( 27 ) / JVS(          108  )
  W( 27 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 109 )
  W( 91 ) = W( 91 ) + a*JVS( 110 )
  W( 95 ) = W( 95 ) + a*JVS( 111 )
  a = -W( 28 ) / JVS(          113  )
  W( 28 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 114 )
  W( 56 ) = W( 56 ) + a*JVS( 115 )
  W( 91 ) = W( 91 ) + a*JVS( 116 )
  W( 94 ) = W( 94 ) + a*JVS( 117 )
  W( 95 ) = W( 95 ) + a*JVS( 118 )
  a = -W( 29 ) / JVS(          119  )
  W( 29 ) = -a
  W( 38 ) = W( 38 ) + a*JVS( 120 )
  W( 39 ) = W( 39 ) + a*JVS( 121 )
  W( 40 ) = W( 40 ) + a*JVS( 122 )
  W( 91 ) = W( 91 ) + a*JVS( 123 )
  W( 94 ) = W( 94 ) + a*JVS( 124 )
  W( 95 ) = W( 95 ) + a*JVS( 125 )
  a = -W( 32 ) / JVS(          139  )
  W( 32 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 140 )
  W( 88 ) = W( 88 ) + a*JVS( 141 )
  W( 89 ) = W( 89 ) + a*JVS( 142 )
  W( 94 ) = W( 94 ) + a*JVS( 143 )
  a = -W( 33 ) / JVS(          144  )
  W( 33 ) = -a
  W( 47 ) = W( 47 ) + a*JVS( 145 )
  W( 51 ) = W( 51 ) + a*JVS( 146 )
  W( 54 ) = W( 54 ) + a*JVS( 147 )
  W( 55 ) = W( 55 ) + a*JVS( 148 )
  W( 59 ) = W( 59 ) + a*JVS( 149 )
  W( 81 ) = W( 81 ) + a*JVS( 150 )
  W( 82 ) = W( 82 ) + a*JVS( 151 )
  W( 85 ) = W( 85 ) + a*JVS( 152 )
  W( 89 ) = W( 89 ) + a*JVS( 153 )
  W( 91 ) = W( 91 ) + a*JVS( 154 )
  a = -W( 34 ) / JVS(          157  )
  W( 34 ) = -a
  W( 35 ) = W( 35 ) + a*JVS( 158 )
  W( 41 ) = W( 41 ) + a*JVS( 159 )
  W( 46 ) = W( 46 ) + a*JVS( 160 )
  W( 49 ) = W( 49 ) + a*JVS( 161 )
  W( 50 ) = W( 50 ) + a*JVS( 162 )
  W( 59 ) = W( 59 ) + a*JVS( 163 )
  W( 73 ) = W( 73 ) + a*JVS( 164 )
  W( 78 ) = W( 78 ) + a*JVS( 165 )
  W( 81 ) = W( 81 ) + a*JVS( 166 )
  W( 88 ) = W( 88 ) + a*JVS( 167 )
  W( 91 ) = W( 91 ) + a*JVS( 168 )
  a = -W( 35 ) / JVS(          169  )
  W( 35 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 170 )
  W( 81 ) = W( 81 ) + a*JVS( 171 )
  W( 82 ) = W( 82 ) + a*JVS( 172 )
  W( 85 ) = W( 85 ) + a*JVS( 173 )
  W( 88 ) = W( 88 ) + a*JVS( 174 )
  W( 91 ) = W( 91 ) + a*JVS( 175 )
  a = -W( 36 ) / JVS(          176  )
  W( 36 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 177 )
  W( 59 ) = W( 59 ) + a*JVS( 178 )
  W( 81 ) = W( 81 ) + a*JVS( 179 )
  W( 82 ) = W( 82 ) + a*JVS( 180 )
  W( 88 ) = W( 88 ) + a*JVS( 181 )
  W( 91 ) = W( 91 ) + a*JVS( 182 )
  a = -W( 37 ) / JVS(          183  )
  W( 37 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 184 )
  W( 58 ) = W( 58 ) + a*JVS( 185 )
  W( 87 ) = W( 87 ) + a*JVS( 186 )
  W( 90 ) = W( 90 ) + a*JVS( 187 )
  W( 91 ) = W( 91 ) + a*JVS( 188 )
  a = -W( 38 ) / JVS(          190  )
  W( 38 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 191 )
  W( 88 ) = W( 88 ) + a*JVS( 192 )
  W( 91 ) = W( 91 ) + a*JVS( 193 )
  W( 95 ) = W( 95 ) + a*JVS( 194 )
  a = -W( 39 ) / JVS(          196  )
  W( 39 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 197 )
  W( 88 ) = W( 88 ) + a*JVS( 198 )
  W( 91 ) = W( 91 ) + a*JVS( 199 )
  W( 95 ) = W( 95 ) + a*JVS( 200 )
  a = -W( 40 ) / JVS(          201  )
  W( 40 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 202 )
  W( 85 ) = W( 85 ) + a*JVS( 203 )
  W( 91 ) = W( 91 ) + a*JVS( 204 )
  W( 95 ) = W( 95 ) + a*JVS( 205 )
  a = -W( 41 ) / JVS(          206  )
  W( 41 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 207 )
  W( 82 ) = W( 82 ) + a*JVS( 208 )
  W( 88 ) = W( 88 ) + a*JVS( 209 )
  W( 91 ) = W( 91 ) + a*JVS( 210 )
  a = -W( 42 ) / JVS(          211  )
  W( 42 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 212 )
  W( 67 ) = W( 67 ) + a*JVS( 213 )
  W( 91 ) = W( 91 ) + a*JVS( 214 )
  W( 94 ) = W( 94 ) + a*JVS( 215 )
  W( 96 ) = W( 96 ) + a*JVS( 216 )
  a = -W( 43 ) / JVS(          217  )
  W( 43 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 218 )
  W( 91 ) = W( 91 ) + a*JVS( 219 )
  W( 96 ) = W( 96 ) + a*JVS( 220 )
  a = -W( 44 ) / JVS(          221  )
  W( 44 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 222 )
  W( 60 ) = W( 60 ) + a*JVS( 223 )
  W( 69 ) = W( 69 ) + a*JVS( 224 )
  W( 73 ) = W( 73 ) + a*JVS( 225 )
  W( 78 ) = W( 78 ) + a*JVS( 226 )
  W( 79 ) = W( 79 ) + a*JVS( 227 )
  W( 89 ) = W( 89 ) + a*JVS( 228 )
  W( 91 ) = W( 91 ) + a*JVS( 229 )
  W( 95 ) = W( 95 ) + a*JVS( 230 )
  W( 96 ) = W( 96 ) + a*JVS( 231 )
  a = -W( 45 ) / JVS(          232  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 233 )
  W( 85 ) = W( 85 ) + a*JVS( 234 )
  W( 89 ) = W( 89 ) + a*JVS( 235 )
  W( 91 ) = W( 91 ) + a*JVS( 236 )
  W( 93 ) = W( 93 ) + a*JVS( 237 )
  a = -W( 46 ) / JVS(          238  )
  W( 46 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 239 )
  W( 82 ) = W( 82 ) + a*JVS( 240 )
  W( 88 ) = W( 88 ) + a*JVS( 241 )
  W( 91 ) = W( 91 ) + a*JVS( 242 )
  a = -W( 47 ) / JVS(          243  )
  W( 47 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 244 )
  W( 91 ) = W( 91 ) + a*JVS( 245 )
  W( 96 ) = W( 96 ) + a*JVS( 246 )
  a = -W( 48 ) / JVS(          247  )
  W( 48 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 248 )
  W( 85 ) = W( 85 ) + a*JVS( 249 )
  W( 91 ) = W( 91 ) + a*JVS( 250 )
  W( 93 ) = W( 93 ) + a*JVS( 251 )
  W( 95 ) = W( 95 ) + a*JVS( 252 )
  W( 96 ) = W( 96 ) + a*JVS( 253 )
  a = -W( 49 ) / JVS(          254  )
  W( 49 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 255 )
  W( 82 ) = W( 82 ) + a*JVS( 256 )
  W( 88 ) = W( 88 ) + a*JVS( 257 )
  W( 91 ) = W( 91 ) + a*JVS( 258 )
  a = -W( 50 ) / JVS(          259  )
  W( 50 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 260 )
  W( 81 ) = W( 81 ) + a*JVS( 261 )
  W( 82 ) = W( 82 ) + a*JVS( 262 )
  W( 85 ) = W( 85 ) + a*JVS( 263 )
  W( 88 ) = W( 88 ) + a*JVS( 264 )
  W( 91 ) = W( 91 ) + a*JVS( 265 )
  a = -W( 51 ) / JVS(          266  )
  W( 51 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 267 )
  W( 91 ) = W( 91 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 52 ) / JVS(          276  )
  W( 52 ) = -a
  W( 53 ) = W( 53 ) + a*JVS( 277 )
  W( 54 ) = W( 54 ) + a*JVS( 278 )
  W( 55 ) = W( 55 ) + a*JVS( 279 )
  W( 56 ) = W( 56 ) + a*JVS( 280 )
  W( 58 ) = W( 58 ) + a*JVS( 281 )
  W( 59 ) = W( 59 ) + a*JVS( 282 )
  W( 60 ) = W( 60 ) + a*JVS( 283 )
  W( 64 ) = W( 64 ) + a*JVS( 284 )
  W( 67 ) = W( 67 ) + a*JVS( 285 )
  W( 69 ) = W( 69 ) + a*JVS( 286 )
  W( 73 ) = W( 73 ) + a*JVS( 287 )
  W( 78 ) = W( 78 ) + a*JVS( 288 )
  W( 79 ) = W( 79 ) + a*JVS( 289 )
  W( 81 ) = W( 81 ) + a*JVS( 290 )
  W( 82 ) = W( 82 ) + a*JVS( 291 )
  W( 85 ) = W( 85 ) + a*JVS( 292 )
  W( 88 ) = W( 88 ) + a*JVS( 293 )
  W( 91 ) = W( 91 ) + a*JVS( 294 )
  W( 94 ) = W( 94 ) + a*JVS( 295 )
  W( 96 ) = W( 96 ) + a*JVS( 296 )
  a = -W( 53 ) / JVS(          297  )
  W( 53 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 298 )
  W( 88 ) = W( 88 ) + a*JVS( 299 )
  W( 91 ) = W( 91 ) + a*JVS( 300 )
  W( 96 ) = W( 96 ) + a*JVS( 301 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 56 ) / JVS(          312  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 313 )
  W( 64 ) = W( 64 ) + a*JVS( 314 )
  W( 85 ) = W( 85 ) + a*JVS( 315 )
  W( 89 ) = W( 89 ) + a*JVS( 316 )
  W( 91 ) = W( 91 ) + a*JVS( 317 )
  W( 94 ) = W( 94 ) + a*JVS( 318 )
  W( 95 ) = W( 95 ) + a*JVS( 319 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 325 )
  W( 88 ) = W( 88 ) + a*JVS( 326 )
  W( 89 ) = W( 89 ) + a*JVS( 327 )
  W( 91 ) = W( 91 ) + a*JVS( 328 )
  W( 95 ) = W( 95 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  a = -W( 58 ) / JVS(          331  )
  W( 58 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 332 )
  W( 85 ) = W( 85 ) + a*JVS( 333 )
  W( 91 ) = W( 91 ) + a*JVS( 334 )
  W( 95 ) = W( 95 ) + a*JVS( 335 )
  W( 96 ) = W( 96 ) + a*JVS( 336 )
  a = -W( 59 ) / JVS(          337  )
  W( 59 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 338 )
  W( 88 ) = W( 88 ) + a*JVS( 339 )
  W( 91 ) = W( 91 ) + a*JVS( 340 )
  W( 96 ) = W( 96 ) + a*JVS( 341 )
  a = -W( 60 ) / JVS(          343  )
  W( 60 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 344 )
  W( 71 ) = W( 71 ) + a*JVS( 345 )
  W( 72 ) = W( 72 ) + a*JVS( 346 )
  W( 79 ) = W( 79 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 82 ) = W( 82 ) + a*JVS( 349 )
  W( 85 ) = W( 85 ) + a*JVS( 350 )
  W( 86 ) = W( 86 ) + a*JVS( 351 )
  W( 88 ) = W( 88 ) + a*JVS( 352 )
  W( 90 ) = W( 90 ) + a*JVS( 353 )
  W( 91 ) = W( 91 ) + a*JVS( 354 )
  W( 93 ) = W( 93 ) + a*JVS( 355 )
  W( 94 ) = W( 94 ) + a*JVS( 356 )
  W( 96 ) = W( 96 ) + a*JVS( 357 )
  a = -W( 64 ) / JVS(          391  )
  W( 64 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 392 )
  W( 85 ) = W( 85 ) + a*JVS( 393 )
  W( 88 ) = W( 88 ) + a*JVS( 394 )
  W( 89 ) = W( 89 ) + a*JVS( 395 )
  W( 90 ) = W( 90 ) + a*JVS( 396 )
  W( 91 ) = W( 91 ) + a*JVS( 397 )
  W( 93 ) = W( 93 ) + a*JVS( 398 )
  W( 94 ) = W( 94 ) + a*JVS( 399 )
  W( 96 ) = W( 96 ) + a*JVS( 400 )
  a = -W( 66 ) / JVS(          417  )
  W( 66 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 75 ) = W( 75 ) + a*JVS( 419 )
  W( 76 ) = W( 76 ) + a*JVS( 420 )
  W( 77 ) = W( 77 ) + a*JVS( 421 )
  W( 79 ) = W( 79 ) + a*JVS( 422 )
  W( 81 ) = W( 81 ) + a*JVS( 423 )
  W( 82 ) = W( 82 ) + a*JVS( 424 )
  W( 83 ) = W( 83 ) + a*JVS( 425 )
  W( 84 ) = W( 84 ) + a*JVS( 426 )
  W( 85 ) = W( 85 ) + a*JVS( 427 )
  W( 86 ) = W( 86 ) + a*JVS( 428 )
  W( 87 ) = W( 87 ) + a*JVS( 429 )
  W( 88 ) = W( 88 ) + a*JVS( 430 )
  W( 89 ) = W( 89 ) + a*JVS( 431 )
  W( 90 ) = W( 90 ) + a*JVS( 432 )
  W( 91 ) = W( 91 ) + a*JVS( 433 )
  W( 92 ) = W( 92 ) + a*JVS( 434 )
  W( 93 ) = W( 93 ) + a*JVS( 435 )
  W( 94 ) = W( 94 ) + a*JVS( 436 )
  W( 96 ) = W( 96 ) + a*JVS( 437 )
  a = -W( 67 ) / JVS(          440  )
  W( 67 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 441 )
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 89 ) = W( 89 ) + a*JVS( 443 )
  W( 90 ) = W( 90 ) + a*JVS( 444 )
  W( 91 ) = W( 91 ) + a*JVS( 445 )
  W( 93 ) = W( 93 ) + a*JVS( 446 )
  W( 94 ) = W( 94 ) + a*JVS( 447 )
  W( 96 ) = W( 96 ) + a*JVS( 448 )
  a = -W( 68 ) / JVS(          451  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 452 )
  W( 88 ) = W( 88 ) + a*JVS( 453 )
  W( 89 ) = W( 89 ) + a*JVS( 454 )
  W( 90 ) = W( 90 ) + a*JVS( 455 )
  W( 91 ) = W( 91 ) + a*JVS( 456 )
  W( 93 ) = W( 93 ) + a*JVS( 457 )
  W( 94 ) = W( 94 ) + a*JVS( 458 )
  W( 95 ) = W( 95 ) + a*JVS( 459 )
  W( 96 ) = W( 96 ) + a*JVS( 460 )
  a = -W( 69 ) / JVS(          471  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 472 )
  W( 72 ) = W( 72 ) + a*JVS( 473 )
  W( 74 ) = W( 74 ) + a*JVS( 474 )
  W( 79 ) = W( 79 ) + a*JVS( 475 )
  W( 81 ) = W( 81 ) + a*JVS( 476 )
  W( 82 ) = W( 82 ) + a*JVS( 477 )
  W( 85 ) = W( 85 ) + a*JVS( 478 )
  W( 86 ) = W( 86 ) + a*JVS( 479 )
  W( 87 ) = W( 87 ) + a*JVS( 480 )
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 89 ) = W( 89 ) + a*JVS( 482 )
  W( 90 ) = W( 90 ) + a*JVS( 483 )
  W( 91 ) = W( 91 ) + a*JVS( 484 )
  W( 93 ) = W( 93 ) + a*JVS( 485 )
  W( 94 ) = W( 94 ) + a*JVS( 486 )
  W( 95 ) = W( 95 ) + a*JVS( 487 )
  W( 96 ) = W( 96 ) + a*JVS( 488 )
  a = -W( 70 ) / JVS(          491  )
  W( 70 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 492 )
  W( 85 ) = W( 85 ) + a*JVS( 493 )
  W( 89 ) = W( 89 ) + a*JVS( 494 )
  W( 90 ) = W( 90 ) + a*JVS( 495 )
  W( 91 ) = W( 91 ) + a*JVS( 496 )
  W( 93 ) = W( 93 ) + a*JVS( 497 )
  W( 94 ) = W( 94 ) + a*JVS( 498 )
  W( 96 ) = W( 96 ) + a*JVS( 499 )
  a = -W( 71 ) / JVS(          501  )
  W( 71 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 502 )
  W( 88 ) = W( 88 ) + a*JVS( 503 )
  W( 89 ) = W( 89 ) + a*JVS( 504 )
  W( 90 ) = W( 90 ) + a*JVS( 505 )
  W( 91 ) = W( 91 ) + a*JVS( 506 )
  W( 93 ) = W( 93 ) + a*JVS( 507 )
  W( 94 ) = W( 94 ) + a*JVS( 508 )
  W( 95 ) = W( 95 ) + a*JVS( 509 )
  W( 96 ) = W( 96 ) + a*JVS( 510 )
  a = -W( 72 ) / JVS(          512  )
  W( 72 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 513 )
  W( 88 ) = W( 88 ) + a*JVS( 514 )
  W( 89 ) = W( 89 ) + a*JVS( 515 )
  W( 90 ) = W( 90 ) + a*JVS( 516 )
  W( 91 ) = W( 91 ) + a*JVS( 517 )
  W( 93 ) = W( 93 ) + a*JVS( 518 )
  W( 94 ) = W( 94 ) + a*JVS( 519 )
  W( 95 ) = W( 95 ) + a*JVS( 520 )
  W( 96 ) = W( 96 ) + a*JVS( 521 )
  a = -W( 73 ) / JVS(          549  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 550 )
  W( 76 ) = W( 76 ) + a*JVS( 551 )
  W( 77 ) = W( 77 ) + a*JVS( 552 )
  W( 79 ) = W( 79 ) + a*JVS( 553 )
  W( 80 ) = W( 80 ) + a*JVS( 554 )
  W( 81 ) = W( 81 ) + a*JVS( 555 )
  W( 82 ) = W( 82 ) + a*JVS( 556 )
  W( 83 ) = W( 83 ) + a*JVS( 557 )
  W( 84 ) = W( 84 ) + a*JVS( 558 )
  W( 85 ) = W( 85 ) + a*JVS( 559 )
  W( 86 ) = W( 86 ) + a*JVS( 560 )
  W( 87 ) = W( 87 ) + a*JVS( 561 )
  W( 88 ) = W( 88 ) + a*JVS( 562 )
  W( 89 ) = W( 89 ) + a*JVS( 563 )
  W( 90 ) = W( 90 ) + a*JVS( 564 )
  W( 91 ) = W( 91 ) + a*JVS( 565 )
  W( 93 ) = W( 93 ) + a*JVS( 566 )
  W( 94 ) = W( 94 ) + a*JVS( 567 )
  W( 95 ) = W( 95 ) + a*JVS( 568 )
  W( 96 ) = W( 96 ) + a*JVS( 569 )
  W( 97 ) = W( 97 ) + a*JVS( 570 )
  a = -W( 74 ) / JVS(          572  )
  W( 74 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 573 )
  W( 85 ) = W( 85 ) + a*JVS( 574 )
  W( 89 ) = W( 89 ) + a*JVS( 575 )
  W( 90 ) = W( 90 ) + a*JVS( 576 )
  W( 91 ) = W( 91 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  a = -W( 75 ) / JVS(          590  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 591 )
  W( 77 ) = W( 77 ) + a*JVS( 592 )
  W( 79 ) = W( 79 ) + a*JVS( 593 )
  W( 83 ) = W( 83 ) + a*JVS( 594 )
  W( 84 ) = W( 84 ) + a*JVS( 595 )
  W( 85 ) = W( 85 ) + a*JVS( 596 )
  W( 86 ) = W( 86 ) + a*JVS( 597 )
  W( 88 ) = W( 88 ) + a*JVS( 598 )
  W( 89 ) = W( 89 ) + a*JVS( 599 )
  W( 90 ) = W( 90 ) + a*JVS( 600 )
  W( 91 ) = W( 91 ) + a*JVS( 601 )
  W( 93 ) = W( 93 ) + a*JVS( 602 )
  W( 94 ) = W( 94 ) + a*JVS( 603 )
  W( 95 ) = W( 95 ) + a*JVS( 604 )
  W( 96 ) = W( 96 ) + a*JVS( 605 )
  a = -W( 76 ) / JVS(          607  )
  W( 76 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 608 )
  W( 82 ) = W( 82 ) + a*JVS( 609 )
  W( 88 ) = W( 88 ) + a*JVS( 610 )
  W( 89 ) = W( 89 ) + a*JVS( 611 )
  W( 90 ) = W( 90 ) + a*JVS( 612 )
  W( 91 ) = W( 91 ) + a*JVS( 613 )
  W( 93 ) = W( 93 ) + a*JVS( 614 )
  W( 94 ) = W( 94 ) + a*JVS( 615 )
  W( 96 ) = W( 96 ) + a*JVS( 616 )
  a = -W( 77 ) / JVS(          618  )
  W( 77 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 619 )
  W( 82 ) = W( 82 ) + a*JVS( 620 )
  W( 88 ) = W( 88 ) + a*JVS( 621 )
  W( 89 ) = W( 89 ) + a*JVS( 622 )
  W( 90 ) = W( 90 ) + a*JVS( 623 )
  W( 91 ) = W( 91 ) + a*JVS( 624 )
  W( 93 ) = W( 93 ) + a*JVS( 625 )
  W( 94 ) = W( 94 ) + a*JVS( 626 )
  W( 96 ) = W( 96 ) + a*JVS( 627 )
  a = -W( 78 ) / JVS(          643  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 644 )
  W( 81 ) = W( 81 ) + a*JVS( 645 )
  W( 82 ) = W( 82 ) + a*JVS( 646 )
  W( 83 ) = W( 83 ) + a*JVS( 647 )
  W( 84 ) = W( 84 ) + a*JVS( 648 )
  W( 85 ) = W( 85 ) + a*JVS( 649 )
  W( 86 ) = W( 86 ) + a*JVS( 650 )
  W( 87 ) = W( 87 ) + a*JVS( 651 )
  W( 88 ) = W( 88 ) + a*JVS( 652 )
  W( 89 ) = W( 89 ) + a*JVS( 653 )
  W( 90 ) = W( 90 ) + a*JVS( 654 )
  W( 91 ) = W( 91 ) + a*JVS( 655 )
  W( 92 ) = W( 92 ) + a*JVS( 656 )
  W( 93 ) = W( 93 ) + a*JVS( 657 )
  W( 94 ) = W( 94 ) + a*JVS( 658 )
  W( 95 ) = W( 95 ) + a*JVS( 659 )
  W( 96 ) = W( 96 ) + a*JVS( 660 )
  W( 97 ) = W( 97 ) + a*JVS( 661 )
  a = -W( 79 ) / JVS(          664  )
  W( 79 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 665 )
  W( 88 ) = W( 88 ) + a*JVS( 666 )
  W( 89 ) = W( 89 ) + a*JVS( 667 )
  W( 90 ) = W( 90 ) + a*JVS( 668 )
  W( 91 ) = W( 91 ) + a*JVS( 669 )
  W( 93 ) = W( 93 ) + a*JVS( 670 )
  W( 94 ) = W( 94 ) + a*JVS( 671 )
  W( 95 ) = W( 95 ) + a*JVS( 672 )
  W( 96 ) = W( 96 ) + a*JVS( 673 )
  a = -W( 80 ) / JVS(          692  )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 693 )
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 85 ) = W( 85 ) + a*JVS( 695 )
  W( 86 ) = W( 86 ) + a*JVS( 696 )
  W( 87 ) = W( 87 ) + a*JVS( 697 )
  W( 88 ) = W( 88 ) + a*JVS( 698 )
  W( 89 ) = W( 89 ) + a*JVS( 699 )
  W( 90 ) = W( 90 ) + a*JVS( 700 )
  W( 91 ) = W( 91 ) + a*JVS( 701 )
  W( 92 ) = W( 92 ) + a*JVS( 702 )
  W( 93 ) = W( 93 ) + a*JVS( 703 )
  W( 94 ) = W( 94 ) + a*JVS( 704 )
  W( 95 ) = W( 95 ) + a*JVS( 705 )
  W( 96 ) = W( 96 ) + a*JVS( 706 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  a = -W( 83 ) / JVS(          739  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 740 )
  W( 85 ) = W( 85 ) + a*JVS( 741 )
  W( 88 ) = W( 88 ) + a*JVS( 742 )
  W( 89 ) = W( 89 ) + a*JVS( 743 )
  W( 90 ) = W( 90 ) + a*JVS( 744 )
  W( 91 ) = W( 91 ) + a*JVS( 745 )
  W( 93 ) = W( 93 ) + a*JVS( 746 )
  W( 94 ) = W( 94 ) + a*JVS( 747 )
  W( 96 ) = W( 96 ) + a*JVS( 748 )
  a = -W( 84 ) / JVS(          758  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 759 )
  W( 88 ) = W( 88 ) + a*JVS( 760 )
  W( 89 ) = W( 89 ) + a*JVS( 761 )
  W( 90 ) = W( 90 ) + a*JVS( 762 )
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 93 ) = W( 93 ) + a*JVS( 764 )
  W( 94 ) = W( 94 ) + a*JVS( 765 )
  W( 96 ) = W( 96 ) + a*JVS( 766 )
  a = -W( 85 ) / JVS(          786  )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 787 )
  W( 89 ) = W( 89 ) + a*JVS( 788 )
  W( 90 ) = W( 90 ) + a*JVS( 789 )
  W( 91 ) = W( 91 ) + a*JVS( 790 )
  W( 93 ) = W( 93 ) + a*JVS( 791 )
  W( 94 ) = W( 94 ) + a*JVS( 792 )
  W( 95 ) = W( 95 ) + a*JVS( 793 )
  W( 96 ) = W( 96 ) + a*JVS( 794 )
  a = -W( 86 ) / JVS(          805  )
  W( 86 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 806 )
  W( 89 ) = W( 89 ) + a*JVS( 807 )
  W( 90 ) = W( 90 ) + a*JVS( 808 )
  W( 91 ) = W( 91 ) + a*JVS( 809 )
  W( 92 ) = W( 92 ) + a*JVS( 810 )
  W( 93 ) = W( 93 ) + a*JVS( 811 )
  W( 94 ) = W( 94 ) + a*JVS( 812 )
  W( 95 ) = W( 95 ) + a*JVS( 813 )
  W( 96 ) = W( 96 ) + a*JVS( 814 )
  a = -W( 87 ) / JVS(          831  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 832 )
  W( 89 ) = W( 89 ) + a*JVS( 833 )
  W( 90 ) = W( 90 ) + a*JVS( 834 )
  W( 91 ) = W( 91 ) + a*JVS( 835 )
  W( 92 ) = W( 92 ) + a*JVS( 836 )
  W( 93 ) = W( 93 ) + a*JVS( 837 )
  W( 94 ) = W( 94 ) + a*JVS( 838 )
  W( 95 ) = W( 95 ) + a*JVS( 839 )
  W( 96 ) = W( 96 ) + a*JVS( 840 )
  a = -W( 88 ) / JVS(          871  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 872 )
  W( 90 ) = W( 90 ) + a*JVS( 873 )
  W( 91 ) = W( 91 ) + a*JVS( 874 )
  W( 92 ) = W( 92 ) + a*JVS( 875 )
  W( 93 ) = W( 93 ) + a*JVS( 876 )
  W( 94 ) = W( 94 ) + a*JVS( 877 )
  W( 95 ) = W( 95 ) + a*JVS( 878 )
  W( 96 ) = W( 96 ) + a*JVS( 879 )
  W( 97 ) = W( 97 ) + a*JVS( 880 )
  a = -W( 89 ) / JVS(          940  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 941 )
  W( 91 ) = W( 91 ) + a*JVS( 942 )
  W( 92 ) = W( 92 ) + a*JVS( 943 )
  W( 93 ) = W( 93 ) + a*JVS( 944 )
  W( 94 ) = W( 94 ) + a*JVS( 945 )
  W( 95 ) = W( 95 ) + a*JVS( 946 )
  W( 96 ) = W( 96 ) + a*JVS( 947 )
  W( 97 ) = W( 97 ) + a*JVS( 948 )
  a = -W( 90 ) / JVS(          979  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 980 )
  W( 92 ) = W( 92 ) + a*JVS( 981 )
  W( 93 ) = W( 93 ) + a*JVS( 982 )
  W( 94 ) = W( 94 ) + a*JVS( 983 )
  W( 95 ) = W( 95 ) + a*JVS( 984 )
  W( 96 ) = W( 96 ) + a*JVS( 985 )
  W( 97 ) = W( 97 ) + a*JVS( 986 )
  JVS( 987) = W( 14 )
  JVS( 988) = W( 15 )
  JVS( 989) = W( 19 )
  JVS( 990) = W( 20 )
  JVS( 991) = W( 21 )
  JVS( 992) = W( 22 )
  JVS( 993) = W( 23 )
  JVS( 994) = W( 24 )
  JVS( 995) = W( 25 )
  JVS( 996) = W( 26 )
  JVS( 997) = W( 27 )
  JVS( 998) = W( 28 )
  JVS( 999) = W( 29 )
  JVS( 1000) = W( 32 )
  JVS( 1001) = W( 33 )
  JVS( 1002) = W( 34 )
  JVS( 1003) = W( 35 )
  JVS( 1004) = W( 36 )
  JVS( 1005) = W( 37 )
  JVS( 1006) = W( 38 )
  JVS( 1007) = W( 39 )
  JVS( 1008) = W( 40 )
  JVS( 1009) = W( 41 )
  JVS( 1010) = W( 42 )
  JVS( 1011) = W( 43 )
  JVS( 1012) = W( 44 )
  JVS( 1013) = W( 45 )
  JVS( 1014) = W( 46 )
  JVS( 1015) = W( 47 )
  JVS( 1016) = W( 48 )
  JVS( 1017) = W( 49 )
  JVS( 1018) = W( 50 )
  JVS( 1019) = W( 51 )
  JVS( 1020) = W( 52 )
  JVS( 1021) = W( 53 )
  JVS( 1022) = W( 54 )
  JVS( 1023) = W( 55 )
  JVS( 1024) = W( 56 )
  JVS( 1025) = W( 57 )
  JVS( 1026) = W( 58 )
  JVS( 1027) = W( 59 )
  JVS( 1028) = W( 60 )
  JVS( 1029) = W( 64 )
  JVS( 1030) = W( 66 )
  JVS( 1031) = W( 67 )
  JVS( 1032) = W( 68 )
  JVS( 1033) = W( 69 )
  JVS( 1034) = W( 70 )
  JVS( 1035) = W( 71 )
  JVS( 1036) = W( 72 )
  JVS( 1037) = W( 73 )
  JVS( 1038) = W( 74 )
  JVS( 1039) = W( 75 )
  JVS( 1040) = W( 76 )
  JVS( 1041) = W( 77 )
  JVS( 1042) = W( 78 )
  JVS( 1043) = W( 79 )
  JVS( 1044) = W( 80 )
  JVS( 1045) = W( 81 )
  JVS( 1046) = W( 82 )
  JVS( 1047) = W( 83 )
  JVS( 1048) = W( 84 )
  JVS( 1049) = W( 85 )
  JVS( 1050) = W( 86 )
  JVS( 1051) = W( 87 )
  JVS( 1052) = W( 88 )
  JVS( 1053) = W( 89 )
  JVS( 1054) = W( 90 )
  JVS( 1055) = W( 91 )
  JVS( 1056) = W( 92 )
  JVS( 1057) = W( 93 )
  JVS( 1058) = W( 94 )
  JVS( 1059) = W( 95 )
  JVS( 1060) = W( 96 )
  JVS( 1061) = W( 97 )
  IF ( ABS(  JVS( 1084 )) < TINY(a) ) THEN
         IER = 92                                      
         RETURN
  END IF
   W( 61 ) = JVS( 1062 )
   W( 62 ) = JVS( 1063 )
   W( 63 ) = JVS( 1064 )
   W( 65 ) = JVS( 1065 )
   W( 68 ) = JVS( 1066 )
   W( 70 ) = JVS( 1067 )
   W( 71 ) = JVS( 1068 )
   W( 72 ) = JVS( 1069 )
   W( 74 ) = JVS( 1070 )
   W( 76 ) = JVS( 1071 )
   W( 77 ) = JVS( 1072 )
   W( 79 ) = JVS( 1073 )
   W( 80 ) = JVS( 1074 )
   W( 81 ) = JVS( 1075 )
   W( 82 ) = JVS( 1076 )
   W( 85 ) = JVS( 1077 )
   W( 86 ) = JVS( 1078 )
   W( 87 ) = JVS( 1079 )
   W( 88 ) = JVS( 1080 )
   W( 89 ) = JVS( 1081 )
   W( 90 ) = JVS( 1082 )
   W( 91 ) = JVS( 1083 )
   W( 92 ) = JVS( 1084 )
   W( 93 ) = JVS( 1085 )
   W( 94 ) = JVS( 1086 )
   W( 95 ) = JVS( 1087 )
   W( 96 ) = JVS( 1088 )
   W( 97 ) = JVS( 1089 )
  a = -W( 61 ) / JVS(          358  )
  W( 61 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 359 )
  W( 89 ) = W( 89 ) + a*JVS( 360 )
  W( 90 ) = W( 90 ) + a*JVS( 361 )
  W( 91 ) = W( 91 ) + a*JVS( 362 )
  W( 93 ) = W( 93 ) + a*JVS( 363 )
  W( 94 ) = W( 94 ) + a*JVS( 364 )
  W( 96 ) = W( 96 ) + a*JVS( 365 )
  a = -W( 62 ) / JVS(          367  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 368 )
  W( 88 ) = W( 88 ) + a*JVS( 369 )
  W( 89 ) = W( 89 ) + a*JVS( 370 )
  W( 90 ) = W( 90 ) + a*JVS( 371 )
  W( 91 ) = W( 91 ) + a*JVS( 372 )
  W( 93 ) = W( 93 ) + a*JVS( 373 )
  W( 94 ) = W( 94 ) + a*JVS( 374 )
  W( 96 ) = W( 96 ) + a*JVS( 375 )
  a = -W( 63 ) / JVS(          377  )
  W( 63 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  W( 90 ) = W( 90 ) + a*JVS( 380 )
  W( 91 ) = W( 91 ) + a*JVS( 381 )
  W( 93 ) = W( 93 ) + a*JVS( 382 )
  W( 94 ) = W( 94 ) + a*JVS( 383 )
  W( 96 ) = W( 96 ) + a*JVS( 384 )
  a = -W( 65 ) / JVS(          402  )
  W( 65 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 403 )
  W( 89 ) = W( 89 ) + a*JVS( 404 )
  W( 90 ) = W( 90 ) + a*JVS( 405 )
  W( 91 ) = W( 91 ) + a*JVS( 406 )
  W( 93 ) = W( 93 ) + a*JVS( 407 )
  W( 94 ) = W( 94 ) + a*JVS( 408 )
  W( 96 ) = W( 96 ) + a*JVS( 409 )
  a = -W( 68 ) / JVS(          451  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 452 )
  W( 88 ) = W( 88 ) + a*JVS( 453 )
  W( 89 ) = W( 89 ) + a*JVS( 454 )
  W( 90 ) = W( 90 ) + a*JVS( 455 )
  W( 91 ) = W( 91 ) + a*JVS( 456 )
  W( 93 ) = W( 93 ) + a*JVS( 457 )
  W( 94 ) = W( 94 ) + a*JVS( 458 )
  W( 95 ) = W( 95 ) + a*JVS( 459 )
  W( 96 ) = W( 96 ) + a*JVS( 460 )
  a = -W( 70 ) / JVS(          491  )
  W( 70 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 492 )
  W( 85 ) = W( 85 ) + a*JVS( 493 )
  W( 89 ) = W( 89 ) + a*JVS( 494 )
  W( 90 ) = W( 90 ) + a*JVS( 495 )
  W( 91 ) = W( 91 ) + a*JVS( 496 )
  W( 93 ) = W( 93 ) + a*JVS( 497 )
  W( 94 ) = W( 94 ) + a*JVS( 498 )
  W( 96 ) = W( 96 ) + a*JVS( 499 )
  a = -W( 71 ) / JVS(          501  )
  W( 71 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 502 )
  W( 88 ) = W( 88 ) + a*JVS( 503 )
  W( 89 ) = W( 89 ) + a*JVS( 504 )
  W( 90 ) = W( 90 ) + a*JVS( 505 )
  W( 91 ) = W( 91 ) + a*JVS( 506 )
  W( 93 ) = W( 93 ) + a*JVS( 507 )
  W( 94 ) = W( 94 ) + a*JVS( 508 )
  W( 95 ) = W( 95 ) + a*JVS( 509 )
  W( 96 ) = W( 96 ) + a*JVS( 510 )
  a = -W( 72 ) / JVS(          512  )
  W( 72 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 513 )
  W( 88 ) = W( 88 ) + a*JVS( 514 )
  W( 89 ) = W( 89 ) + a*JVS( 515 )
  W( 90 ) = W( 90 ) + a*JVS( 516 )
  W( 91 ) = W( 91 ) + a*JVS( 517 )
  W( 93 ) = W( 93 ) + a*JVS( 518 )
  W( 94 ) = W( 94 ) + a*JVS( 519 )
  W( 95 ) = W( 95 ) + a*JVS( 520 )
  W( 96 ) = W( 96 ) + a*JVS( 521 )
  a = -W( 74 ) / JVS(          572  )
  W( 74 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 573 )
  W( 85 ) = W( 85 ) + a*JVS( 574 )
  W( 89 ) = W( 89 ) + a*JVS( 575 )
  W( 90 ) = W( 90 ) + a*JVS( 576 )
  W( 91 ) = W( 91 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  a = -W( 76 ) / JVS(          607  )
  W( 76 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 608 )
  W( 82 ) = W( 82 ) + a*JVS( 609 )
  W( 88 ) = W( 88 ) + a*JVS( 610 )
  W( 89 ) = W( 89 ) + a*JVS( 611 )
  W( 90 ) = W( 90 ) + a*JVS( 612 )
  W( 91 ) = W( 91 ) + a*JVS( 613 )
  W( 93 ) = W( 93 ) + a*JVS( 614 )
  W( 94 ) = W( 94 ) + a*JVS( 615 )
  W( 96 ) = W( 96 ) + a*JVS( 616 )
  a = -W( 77 ) / JVS(          618  )
  W( 77 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 619 )
  W( 82 ) = W( 82 ) + a*JVS( 620 )
  W( 88 ) = W( 88 ) + a*JVS( 621 )
  W( 89 ) = W( 89 ) + a*JVS( 622 )
  W( 90 ) = W( 90 ) + a*JVS( 623 )
  W( 91 ) = W( 91 ) + a*JVS( 624 )
  W( 93 ) = W( 93 ) + a*JVS( 625 )
  W( 94 ) = W( 94 ) + a*JVS( 626 )
  W( 96 ) = W( 96 ) + a*JVS( 627 )
  a = -W( 79 ) / JVS(          664  )
  W( 79 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 665 )
  W( 88 ) = W( 88 ) + a*JVS( 666 )
  W( 89 ) = W( 89 ) + a*JVS( 667 )
  W( 90 ) = W( 90 ) + a*JVS( 668 )
  W( 91 ) = W( 91 ) + a*JVS( 669 )
  W( 93 ) = W( 93 ) + a*JVS( 670 )
  W( 94 ) = W( 94 ) + a*JVS( 671 )
  W( 95 ) = W( 95 ) + a*JVS( 672 )
  W( 96 ) = W( 96 ) + a*JVS( 673 )
  a = -W( 80 ) / JVS(          692  )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 693 )
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 85 ) = W( 85 ) + a*JVS( 695 )
  W( 86 ) = W( 86 ) + a*JVS( 696 )
  W( 87 ) = W( 87 ) + a*JVS( 697 )
  W( 88 ) = W( 88 ) + a*JVS( 698 )
  W( 89 ) = W( 89 ) + a*JVS( 699 )
  W( 90 ) = W( 90 ) + a*JVS( 700 )
  W( 91 ) = W( 91 ) + a*JVS( 701 )
  W( 92 ) = W( 92 ) + a*JVS( 702 )
  W( 93 ) = W( 93 ) + a*JVS( 703 )
  W( 94 ) = W( 94 ) + a*JVS( 704 )
  W( 95 ) = W( 95 ) + a*JVS( 705 )
  W( 96 ) = W( 96 ) + a*JVS( 706 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  a = -W( 85 ) / JVS(          786  )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 787 )
  W( 89 ) = W( 89 ) + a*JVS( 788 )
  W( 90 ) = W( 90 ) + a*JVS( 789 )
  W( 91 ) = W( 91 ) + a*JVS( 790 )
  W( 93 ) = W( 93 ) + a*JVS( 791 )
  W( 94 ) = W( 94 ) + a*JVS( 792 )
  W( 95 ) = W( 95 ) + a*JVS( 793 )
  W( 96 ) = W( 96 ) + a*JVS( 794 )
  a = -W( 86 ) / JVS(          805  )
  W( 86 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 806 )
  W( 89 ) = W( 89 ) + a*JVS( 807 )
  W( 90 ) = W( 90 ) + a*JVS( 808 )
  W( 91 ) = W( 91 ) + a*JVS( 809 )
  W( 92 ) = W( 92 ) + a*JVS( 810 )
  W( 93 ) = W( 93 ) + a*JVS( 811 )
  W( 94 ) = W( 94 ) + a*JVS( 812 )
  W( 95 ) = W( 95 ) + a*JVS( 813 )
  W( 96 ) = W( 96 ) + a*JVS( 814 )
  a = -W( 87 ) / JVS(          831  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 832 )
  W( 89 ) = W( 89 ) + a*JVS( 833 )
  W( 90 ) = W( 90 ) + a*JVS( 834 )
  W( 91 ) = W( 91 ) + a*JVS( 835 )
  W( 92 ) = W( 92 ) + a*JVS( 836 )
  W( 93 ) = W( 93 ) + a*JVS( 837 )
  W( 94 ) = W( 94 ) + a*JVS( 838 )
  W( 95 ) = W( 95 ) + a*JVS( 839 )
  W( 96 ) = W( 96 ) + a*JVS( 840 )
  a = -W( 88 ) / JVS(          871  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 872 )
  W( 90 ) = W( 90 ) + a*JVS( 873 )
  W( 91 ) = W( 91 ) + a*JVS( 874 )
  W( 92 ) = W( 92 ) + a*JVS( 875 )
  W( 93 ) = W( 93 ) + a*JVS( 876 )
  W( 94 ) = W( 94 ) + a*JVS( 877 )
  W( 95 ) = W( 95 ) + a*JVS( 878 )
  W( 96 ) = W( 96 ) + a*JVS( 879 )
  W( 97 ) = W( 97 ) + a*JVS( 880 )
  a = -W( 89 ) / JVS(          940  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 941 )
  W( 91 ) = W( 91 ) + a*JVS( 942 )
  W( 92 ) = W( 92 ) + a*JVS( 943 )
  W( 93 ) = W( 93 ) + a*JVS( 944 )
  W( 94 ) = W( 94 ) + a*JVS( 945 )
  W( 95 ) = W( 95 ) + a*JVS( 946 )
  W( 96 ) = W( 96 ) + a*JVS( 947 )
  W( 97 ) = W( 97 ) + a*JVS( 948 )
  a = -W( 90 ) / JVS(          979  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 980 )
  W( 92 ) = W( 92 ) + a*JVS( 981 )
  W( 93 ) = W( 93 ) + a*JVS( 982 )
  W( 94 ) = W( 94 ) + a*JVS( 983 )
  W( 95 ) = W( 95 ) + a*JVS( 984 )
  W( 96 ) = W( 96 ) + a*JVS( 985 )
  W( 97 ) = W( 97 ) + a*JVS( 986 )
  a = -W( 91 ) / JVS(         1055  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 1056 )
  W( 93 ) = W( 93 ) + a*JVS( 1057 )
  W( 94 ) = W( 94 ) + a*JVS( 1058 )
  W( 95 ) = W( 95 ) + a*JVS( 1059 )
  W( 96 ) = W( 96 ) + a*JVS( 1060 )
  W( 97 ) = W( 97 ) + a*JVS( 1061 )
  JVS( 1062) = W( 61 )
  JVS( 1063) = W( 62 )
  JVS( 1064) = W( 63 )
  JVS( 1065) = W( 65 )
  JVS( 1066) = W( 68 )
  JVS( 1067) = W( 70 )
  JVS( 1068) = W( 71 )
  JVS( 1069) = W( 72 )
  JVS( 1070) = W( 74 )
  JVS( 1071) = W( 76 )
  JVS( 1072) = W( 77 )
  JVS( 1073) = W( 79 )
  JVS( 1074) = W( 80 )
  JVS( 1075) = W( 81 )
  JVS( 1076) = W( 82 )
  JVS( 1077) = W( 85 )
  JVS( 1078) = W( 86 )
  JVS( 1079) = W( 87 )
  JVS( 1080) = W( 88 )
  JVS( 1081) = W( 89 )
  JVS( 1082) = W( 90 )
  JVS( 1083) = W( 91 )
  JVS( 1084) = W( 92 )
  JVS( 1085) = W( 93 )
  JVS( 1086) = W( 94 )
  JVS( 1087) = W( 95 )
  JVS( 1088) = W( 96 )
  JVS( 1089) = W( 97 )
  IF ( ABS(  JVS( 1129 )) < TINY(a) ) THEN
         IER = 93                                      
         RETURN
  END IF
   W( 37 ) = JVS( 1090 )
   W( 45 ) = JVS( 1091 )
   W( 46 ) = JVS( 1092 )
   W( 48 ) = JVS( 1093 )
   W( 55 ) = JVS( 1094 )
   W( 56 ) = JVS( 1095 )
   W( 58 ) = JVS( 1096 )
   W( 59 ) = JVS( 1097 )
   W( 61 ) = JVS( 1098 )
   W( 62 ) = JVS( 1099 )
   W( 63 ) = JVS( 1100 )
   W( 64 ) = JVS( 1101 )
   W( 65 ) = JVS( 1102 )
   W( 66 ) = JVS( 1103 )
   W( 67 ) = JVS( 1104 )
   W( 68 ) = JVS( 1105 )
   W( 69 ) = JVS( 1106 )
   W( 70 ) = JVS( 1107 )
   W( 71 ) = JVS( 1108 )
   W( 72 ) = JVS( 1109 )
   W( 74 ) = JVS( 1110 )
   W( 75 ) = JVS( 1111 )
   W( 76 ) = JVS( 1112 )
   W( 77 ) = JVS( 1113 )
   W( 78 ) = JVS( 1114 )
   W( 79 ) = JVS( 1115 )
   W( 80 ) = JVS( 1116 )
   W( 81 ) = JVS( 1117 )
   W( 82 ) = JVS( 1118 )
   W( 83 ) = JVS( 1119 )
   W( 84 ) = JVS( 1120 )
   W( 85 ) = JVS( 1121 )
   W( 86 ) = JVS( 1122 )
   W( 87 ) = JVS( 1123 )
   W( 88 ) = JVS( 1124 )
   W( 89 ) = JVS( 1125 )
   W( 90 ) = JVS( 1126 )
   W( 91 ) = JVS( 1127 )
   W( 92 ) = JVS( 1128 )
   W( 93 ) = JVS( 1129 )
   W( 94 ) = JVS( 1130 )
   W( 95 ) = JVS( 1131 )
   W( 96 ) = JVS( 1132 )
   W( 97 ) = JVS( 1133 )
  a = -W( 37 ) / JVS(          183  )
  W( 37 ) = -a
  W( 46 ) = W( 46 ) + a*JVS( 184 )
  W( 58 ) = W( 58 ) + a*JVS( 185 )
  W( 87 ) = W( 87 ) + a*JVS( 186 )
  W( 90 ) = W( 90 ) + a*JVS( 187 )
  W( 91 ) = W( 91 ) + a*JVS( 188 )
  a = -W( 45 ) / JVS(          232  )
  W( 45 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 233 )
  W( 85 ) = W( 85 ) + a*JVS( 234 )
  W( 89 ) = W( 89 ) + a*JVS( 235 )
  W( 91 ) = W( 91 ) + a*JVS( 236 )
  W( 93 ) = W( 93 ) + a*JVS( 237 )
  a = -W( 46 ) / JVS(          238  )
  W( 46 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 239 )
  W( 82 ) = W( 82 ) + a*JVS( 240 )
  W( 88 ) = W( 88 ) + a*JVS( 241 )
  W( 91 ) = W( 91 ) + a*JVS( 242 )
  a = -W( 48 ) / JVS(          247  )
  W( 48 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 248 )
  W( 85 ) = W( 85 ) + a*JVS( 249 )
  W( 91 ) = W( 91 ) + a*JVS( 250 )
  W( 93 ) = W( 93 ) + a*JVS( 251 )
  W( 95 ) = W( 95 ) + a*JVS( 252 )
  W( 96 ) = W( 96 ) + a*JVS( 253 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 56 ) / JVS(          312  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 313 )
  W( 64 ) = W( 64 ) + a*JVS( 314 )
  W( 85 ) = W( 85 ) + a*JVS( 315 )
  W( 89 ) = W( 89 ) + a*JVS( 316 )
  W( 91 ) = W( 91 ) + a*JVS( 317 )
  W( 94 ) = W( 94 ) + a*JVS( 318 )
  W( 95 ) = W( 95 ) + a*JVS( 319 )
  a = -W( 58 ) / JVS(          331  )
  W( 58 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 332 )
  W( 85 ) = W( 85 ) + a*JVS( 333 )
  W( 91 ) = W( 91 ) + a*JVS( 334 )
  W( 95 ) = W( 95 ) + a*JVS( 335 )
  W( 96 ) = W( 96 ) + a*JVS( 336 )
  a = -W( 59 ) / JVS(          337  )
  W( 59 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 338 )
  W( 88 ) = W( 88 ) + a*JVS( 339 )
  W( 91 ) = W( 91 ) + a*JVS( 340 )
  W( 96 ) = W( 96 ) + a*JVS( 341 )
  a = -W( 61 ) / JVS(          358  )
  W( 61 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 359 )
  W( 89 ) = W( 89 ) + a*JVS( 360 )
  W( 90 ) = W( 90 ) + a*JVS( 361 )
  W( 91 ) = W( 91 ) + a*JVS( 362 )
  W( 93 ) = W( 93 ) + a*JVS( 363 )
  W( 94 ) = W( 94 ) + a*JVS( 364 )
  W( 96 ) = W( 96 ) + a*JVS( 365 )
  a = -W( 62 ) / JVS(          367  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 368 )
  W( 88 ) = W( 88 ) + a*JVS( 369 )
  W( 89 ) = W( 89 ) + a*JVS( 370 )
  W( 90 ) = W( 90 ) + a*JVS( 371 )
  W( 91 ) = W( 91 ) + a*JVS( 372 )
  W( 93 ) = W( 93 ) + a*JVS( 373 )
  W( 94 ) = W( 94 ) + a*JVS( 374 )
  W( 96 ) = W( 96 ) + a*JVS( 375 )
  a = -W( 63 ) / JVS(          377  )
  W( 63 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  W( 90 ) = W( 90 ) + a*JVS( 380 )
  W( 91 ) = W( 91 ) + a*JVS( 381 )
  W( 93 ) = W( 93 ) + a*JVS( 382 )
  W( 94 ) = W( 94 ) + a*JVS( 383 )
  W( 96 ) = W( 96 ) + a*JVS( 384 )
  a = -W( 64 ) / JVS(          391  )
  W( 64 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 392 )
  W( 85 ) = W( 85 ) + a*JVS( 393 )
  W( 88 ) = W( 88 ) + a*JVS( 394 )
  W( 89 ) = W( 89 ) + a*JVS( 395 )
  W( 90 ) = W( 90 ) + a*JVS( 396 )
  W( 91 ) = W( 91 ) + a*JVS( 397 )
  W( 93 ) = W( 93 ) + a*JVS( 398 )
  W( 94 ) = W( 94 ) + a*JVS( 399 )
  W( 96 ) = W( 96 ) + a*JVS( 400 )
  a = -W( 65 ) / JVS(          402  )
  W( 65 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 403 )
  W( 89 ) = W( 89 ) + a*JVS( 404 )
  W( 90 ) = W( 90 ) + a*JVS( 405 )
  W( 91 ) = W( 91 ) + a*JVS( 406 )
  W( 93 ) = W( 93 ) + a*JVS( 407 )
  W( 94 ) = W( 94 ) + a*JVS( 408 )
  W( 96 ) = W( 96 ) + a*JVS( 409 )
  a = -W( 66 ) / JVS(          417  )
  W( 66 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 75 ) = W( 75 ) + a*JVS( 419 )
  W( 76 ) = W( 76 ) + a*JVS( 420 )
  W( 77 ) = W( 77 ) + a*JVS( 421 )
  W( 79 ) = W( 79 ) + a*JVS( 422 )
  W( 81 ) = W( 81 ) + a*JVS( 423 )
  W( 82 ) = W( 82 ) + a*JVS( 424 )
  W( 83 ) = W( 83 ) + a*JVS( 425 )
  W( 84 ) = W( 84 ) + a*JVS( 426 )
  W( 85 ) = W( 85 ) + a*JVS( 427 )
  W( 86 ) = W( 86 ) + a*JVS( 428 )
  W( 87 ) = W( 87 ) + a*JVS( 429 )
  W( 88 ) = W( 88 ) + a*JVS( 430 )
  W( 89 ) = W( 89 ) + a*JVS( 431 )
  W( 90 ) = W( 90 ) + a*JVS( 432 )
  W( 91 ) = W( 91 ) + a*JVS( 433 )
  W( 92 ) = W( 92 ) + a*JVS( 434 )
  W( 93 ) = W( 93 ) + a*JVS( 435 )
  W( 94 ) = W( 94 ) + a*JVS( 436 )
  W( 96 ) = W( 96 ) + a*JVS( 437 )
  a = -W( 67 ) / JVS(          440  )
  W( 67 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 441 )
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 89 ) = W( 89 ) + a*JVS( 443 )
  W( 90 ) = W( 90 ) + a*JVS( 444 )
  W( 91 ) = W( 91 ) + a*JVS( 445 )
  W( 93 ) = W( 93 ) + a*JVS( 446 )
  W( 94 ) = W( 94 ) + a*JVS( 447 )
  W( 96 ) = W( 96 ) + a*JVS( 448 )
  a = -W( 68 ) / JVS(          451  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 452 )
  W( 88 ) = W( 88 ) + a*JVS( 453 )
  W( 89 ) = W( 89 ) + a*JVS( 454 )
  W( 90 ) = W( 90 ) + a*JVS( 455 )
  W( 91 ) = W( 91 ) + a*JVS( 456 )
  W( 93 ) = W( 93 ) + a*JVS( 457 )
  W( 94 ) = W( 94 ) + a*JVS( 458 )
  W( 95 ) = W( 95 ) + a*JVS( 459 )
  W( 96 ) = W( 96 ) + a*JVS( 460 )
  a = -W( 69 ) / JVS(          471  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 472 )
  W( 72 ) = W( 72 ) + a*JVS( 473 )
  W( 74 ) = W( 74 ) + a*JVS( 474 )
  W( 79 ) = W( 79 ) + a*JVS( 475 )
  W( 81 ) = W( 81 ) + a*JVS( 476 )
  W( 82 ) = W( 82 ) + a*JVS( 477 )
  W( 85 ) = W( 85 ) + a*JVS( 478 )
  W( 86 ) = W( 86 ) + a*JVS( 479 )
  W( 87 ) = W( 87 ) + a*JVS( 480 )
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 89 ) = W( 89 ) + a*JVS( 482 )
  W( 90 ) = W( 90 ) + a*JVS( 483 )
  W( 91 ) = W( 91 ) + a*JVS( 484 )
  W( 93 ) = W( 93 ) + a*JVS( 485 )
  W( 94 ) = W( 94 ) + a*JVS( 486 )
  W( 95 ) = W( 95 ) + a*JVS( 487 )
  W( 96 ) = W( 96 ) + a*JVS( 488 )
  a = -W( 70 ) / JVS(          491  )
  W( 70 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 492 )
  W( 85 ) = W( 85 ) + a*JVS( 493 )
  W( 89 ) = W( 89 ) + a*JVS( 494 )
  W( 90 ) = W( 90 ) + a*JVS( 495 )
  W( 91 ) = W( 91 ) + a*JVS( 496 )
  W( 93 ) = W( 93 ) + a*JVS( 497 )
  W( 94 ) = W( 94 ) + a*JVS( 498 )
  W( 96 ) = W( 96 ) + a*JVS( 499 )
  a = -W( 71 ) / JVS(          501  )
  W( 71 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 502 )
  W( 88 ) = W( 88 ) + a*JVS( 503 )
  W( 89 ) = W( 89 ) + a*JVS( 504 )
  W( 90 ) = W( 90 ) + a*JVS( 505 )
  W( 91 ) = W( 91 ) + a*JVS( 506 )
  W( 93 ) = W( 93 ) + a*JVS( 507 )
  W( 94 ) = W( 94 ) + a*JVS( 508 )
  W( 95 ) = W( 95 ) + a*JVS( 509 )
  W( 96 ) = W( 96 ) + a*JVS( 510 )
  a = -W( 72 ) / JVS(          512  )
  W( 72 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 513 )
  W( 88 ) = W( 88 ) + a*JVS( 514 )
  W( 89 ) = W( 89 ) + a*JVS( 515 )
  W( 90 ) = W( 90 ) + a*JVS( 516 )
  W( 91 ) = W( 91 ) + a*JVS( 517 )
  W( 93 ) = W( 93 ) + a*JVS( 518 )
  W( 94 ) = W( 94 ) + a*JVS( 519 )
  W( 95 ) = W( 95 ) + a*JVS( 520 )
  W( 96 ) = W( 96 ) + a*JVS( 521 )
  a = -W( 74 ) / JVS(          572  )
  W( 74 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 573 )
  W( 85 ) = W( 85 ) + a*JVS( 574 )
  W( 89 ) = W( 89 ) + a*JVS( 575 )
  W( 90 ) = W( 90 ) + a*JVS( 576 )
  W( 91 ) = W( 91 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  a = -W( 75 ) / JVS(          590  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 591 )
  W( 77 ) = W( 77 ) + a*JVS( 592 )
  W( 79 ) = W( 79 ) + a*JVS( 593 )
  W( 83 ) = W( 83 ) + a*JVS( 594 )
  W( 84 ) = W( 84 ) + a*JVS( 595 )
  W( 85 ) = W( 85 ) + a*JVS( 596 )
  W( 86 ) = W( 86 ) + a*JVS( 597 )
  W( 88 ) = W( 88 ) + a*JVS( 598 )
  W( 89 ) = W( 89 ) + a*JVS( 599 )
  W( 90 ) = W( 90 ) + a*JVS( 600 )
  W( 91 ) = W( 91 ) + a*JVS( 601 )
  W( 93 ) = W( 93 ) + a*JVS( 602 )
  W( 94 ) = W( 94 ) + a*JVS( 603 )
  W( 95 ) = W( 95 ) + a*JVS( 604 )
  W( 96 ) = W( 96 ) + a*JVS( 605 )
  a = -W( 76 ) / JVS(          607  )
  W( 76 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 608 )
  W( 82 ) = W( 82 ) + a*JVS( 609 )
  W( 88 ) = W( 88 ) + a*JVS( 610 )
  W( 89 ) = W( 89 ) + a*JVS( 611 )
  W( 90 ) = W( 90 ) + a*JVS( 612 )
  W( 91 ) = W( 91 ) + a*JVS( 613 )
  W( 93 ) = W( 93 ) + a*JVS( 614 )
  W( 94 ) = W( 94 ) + a*JVS( 615 )
  W( 96 ) = W( 96 ) + a*JVS( 616 )
  a = -W( 77 ) / JVS(          618  )
  W( 77 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 619 )
  W( 82 ) = W( 82 ) + a*JVS( 620 )
  W( 88 ) = W( 88 ) + a*JVS( 621 )
  W( 89 ) = W( 89 ) + a*JVS( 622 )
  W( 90 ) = W( 90 ) + a*JVS( 623 )
  W( 91 ) = W( 91 ) + a*JVS( 624 )
  W( 93 ) = W( 93 ) + a*JVS( 625 )
  W( 94 ) = W( 94 ) + a*JVS( 626 )
  W( 96 ) = W( 96 ) + a*JVS( 627 )
  a = -W( 78 ) / JVS(          643  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 644 )
  W( 81 ) = W( 81 ) + a*JVS( 645 )
  W( 82 ) = W( 82 ) + a*JVS( 646 )
  W( 83 ) = W( 83 ) + a*JVS( 647 )
  W( 84 ) = W( 84 ) + a*JVS( 648 )
  W( 85 ) = W( 85 ) + a*JVS( 649 )
  W( 86 ) = W( 86 ) + a*JVS( 650 )
  W( 87 ) = W( 87 ) + a*JVS( 651 )
  W( 88 ) = W( 88 ) + a*JVS( 652 )
  W( 89 ) = W( 89 ) + a*JVS( 653 )
  W( 90 ) = W( 90 ) + a*JVS( 654 )
  W( 91 ) = W( 91 ) + a*JVS( 655 )
  W( 92 ) = W( 92 ) + a*JVS( 656 )
  W( 93 ) = W( 93 ) + a*JVS( 657 )
  W( 94 ) = W( 94 ) + a*JVS( 658 )
  W( 95 ) = W( 95 ) + a*JVS( 659 )
  W( 96 ) = W( 96 ) + a*JVS( 660 )
  W( 97 ) = W( 97 ) + a*JVS( 661 )
  a = -W( 79 ) / JVS(          664  )
  W( 79 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 665 )
  W( 88 ) = W( 88 ) + a*JVS( 666 )
  W( 89 ) = W( 89 ) + a*JVS( 667 )
  W( 90 ) = W( 90 ) + a*JVS( 668 )
  W( 91 ) = W( 91 ) + a*JVS( 669 )
  W( 93 ) = W( 93 ) + a*JVS( 670 )
  W( 94 ) = W( 94 ) + a*JVS( 671 )
  W( 95 ) = W( 95 ) + a*JVS( 672 )
  W( 96 ) = W( 96 ) + a*JVS( 673 )
  a = -W( 80 ) / JVS(          692  )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 693 )
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 85 ) = W( 85 ) + a*JVS( 695 )
  W( 86 ) = W( 86 ) + a*JVS( 696 )
  W( 87 ) = W( 87 ) + a*JVS( 697 )
  W( 88 ) = W( 88 ) + a*JVS( 698 )
  W( 89 ) = W( 89 ) + a*JVS( 699 )
  W( 90 ) = W( 90 ) + a*JVS( 700 )
  W( 91 ) = W( 91 ) + a*JVS( 701 )
  W( 92 ) = W( 92 ) + a*JVS( 702 )
  W( 93 ) = W( 93 ) + a*JVS( 703 )
  W( 94 ) = W( 94 ) + a*JVS( 704 )
  W( 95 ) = W( 95 ) + a*JVS( 705 )
  W( 96 ) = W( 96 ) + a*JVS( 706 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  a = -W( 83 ) / JVS(          739  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 740 )
  W( 85 ) = W( 85 ) + a*JVS( 741 )
  W( 88 ) = W( 88 ) + a*JVS( 742 )
  W( 89 ) = W( 89 ) + a*JVS( 743 )
  W( 90 ) = W( 90 ) + a*JVS( 744 )
  W( 91 ) = W( 91 ) + a*JVS( 745 )
  W( 93 ) = W( 93 ) + a*JVS( 746 )
  W( 94 ) = W( 94 ) + a*JVS( 747 )
  W( 96 ) = W( 96 ) + a*JVS( 748 )
  a = -W( 84 ) / JVS(          758  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 759 )
  W( 88 ) = W( 88 ) + a*JVS( 760 )
  W( 89 ) = W( 89 ) + a*JVS( 761 )
  W( 90 ) = W( 90 ) + a*JVS( 762 )
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 93 ) = W( 93 ) + a*JVS( 764 )
  W( 94 ) = W( 94 ) + a*JVS( 765 )
  W( 96 ) = W( 96 ) + a*JVS( 766 )
  a = -W( 85 ) / JVS(          786  )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 787 )
  W( 89 ) = W( 89 ) + a*JVS( 788 )
  W( 90 ) = W( 90 ) + a*JVS( 789 )
  W( 91 ) = W( 91 ) + a*JVS( 790 )
  W( 93 ) = W( 93 ) + a*JVS( 791 )
  W( 94 ) = W( 94 ) + a*JVS( 792 )
  W( 95 ) = W( 95 ) + a*JVS( 793 )
  W( 96 ) = W( 96 ) + a*JVS( 794 )
  a = -W( 86 ) / JVS(          805  )
  W( 86 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 806 )
  W( 89 ) = W( 89 ) + a*JVS( 807 )
  W( 90 ) = W( 90 ) + a*JVS( 808 )
  W( 91 ) = W( 91 ) + a*JVS( 809 )
  W( 92 ) = W( 92 ) + a*JVS( 810 )
  W( 93 ) = W( 93 ) + a*JVS( 811 )
  W( 94 ) = W( 94 ) + a*JVS( 812 )
  W( 95 ) = W( 95 ) + a*JVS( 813 )
  W( 96 ) = W( 96 ) + a*JVS( 814 )
  a = -W( 87 ) / JVS(          831  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 832 )
  W( 89 ) = W( 89 ) + a*JVS( 833 )
  W( 90 ) = W( 90 ) + a*JVS( 834 )
  W( 91 ) = W( 91 ) + a*JVS( 835 )
  W( 92 ) = W( 92 ) + a*JVS( 836 )
  W( 93 ) = W( 93 ) + a*JVS( 837 )
  W( 94 ) = W( 94 ) + a*JVS( 838 )
  W( 95 ) = W( 95 ) + a*JVS( 839 )
  W( 96 ) = W( 96 ) + a*JVS( 840 )
  a = -W( 88 ) / JVS(          871  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 872 )
  W( 90 ) = W( 90 ) + a*JVS( 873 )
  W( 91 ) = W( 91 ) + a*JVS( 874 )
  W( 92 ) = W( 92 ) + a*JVS( 875 )
  W( 93 ) = W( 93 ) + a*JVS( 876 )
  W( 94 ) = W( 94 ) + a*JVS( 877 )
  W( 95 ) = W( 95 ) + a*JVS( 878 )
  W( 96 ) = W( 96 ) + a*JVS( 879 )
  W( 97 ) = W( 97 ) + a*JVS( 880 )
  a = -W( 89 ) / JVS(          940  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 941 )
  W( 91 ) = W( 91 ) + a*JVS( 942 )
  W( 92 ) = W( 92 ) + a*JVS( 943 )
  W( 93 ) = W( 93 ) + a*JVS( 944 )
  W( 94 ) = W( 94 ) + a*JVS( 945 )
  W( 95 ) = W( 95 ) + a*JVS( 946 )
  W( 96 ) = W( 96 ) + a*JVS( 947 )
  W( 97 ) = W( 97 ) + a*JVS( 948 )
  a = -W( 90 ) / JVS(          979  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 980 )
  W( 92 ) = W( 92 ) + a*JVS( 981 )
  W( 93 ) = W( 93 ) + a*JVS( 982 )
  W( 94 ) = W( 94 ) + a*JVS( 983 )
  W( 95 ) = W( 95 ) + a*JVS( 984 )
  W( 96 ) = W( 96 ) + a*JVS( 985 )
  W( 97 ) = W( 97 ) + a*JVS( 986 )
  a = -W( 91 ) / JVS(         1055  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 1056 )
  W( 93 ) = W( 93 ) + a*JVS( 1057 )
  W( 94 ) = W( 94 ) + a*JVS( 1058 )
  W( 95 ) = W( 95 ) + a*JVS( 1059 )
  W( 96 ) = W( 96 ) + a*JVS( 1060 )
  W( 97 ) = W( 97 ) + a*JVS( 1061 )
  a = -W( 92 ) / JVS(         1084  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 1085 )
  W( 94 ) = W( 94 ) + a*JVS( 1086 )
  W( 95 ) = W( 95 ) + a*JVS( 1087 )
  W( 96 ) = W( 96 ) + a*JVS( 1088 )
  W( 97 ) = W( 97 ) + a*JVS( 1089 )
  JVS( 1090) = W( 37 )
  JVS( 1091) = W( 45 )
  JVS( 1092) = W( 46 )
  JVS( 1093) = W( 48 )
  JVS( 1094) = W( 55 )
  JVS( 1095) = W( 56 )
  JVS( 1096) = W( 58 )
  JVS( 1097) = W( 59 )
  JVS( 1098) = W( 61 )
  JVS( 1099) = W( 62 )
  JVS( 1100) = W( 63 )
  JVS( 1101) = W( 64 )
  JVS( 1102) = W( 65 )
  JVS( 1103) = W( 66 )
  JVS( 1104) = W( 67 )
  JVS( 1105) = W( 68 )
  JVS( 1106) = W( 69 )
  JVS( 1107) = W( 70 )
  JVS( 1108) = W( 71 )
  JVS( 1109) = W( 72 )
  JVS( 1110) = W( 74 )
  JVS( 1111) = W( 75 )
  JVS( 1112) = W( 76 )
  JVS( 1113) = W( 77 )
  JVS( 1114) = W( 78 )
  JVS( 1115) = W( 79 )
  JVS( 1116) = W( 80 )
  JVS( 1117) = W( 81 )
  JVS( 1118) = W( 82 )
  JVS( 1119) = W( 83 )
  JVS( 1120) = W( 84 )
  JVS( 1121) = W( 85 )
  JVS( 1122) = W( 86 )
  JVS( 1123) = W( 87 )
  JVS( 1124) = W( 88 )
  JVS( 1125) = W( 89 )
  JVS( 1126) = W( 90 )
  JVS( 1127) = W( 91 )
  JVS( 1128) = W( 92 )
  JVS( 1129) = W( 93 )
  JVS( 1130) = W( 94 )
  JVS( 1131) = W( 95 )
  JVS( 1132) = W( 96 )
  JVS( 1133) = W( 97 )
  IF ( ABS(  JVS( 1172 )) < TINY(a) ) THEN
         IER = 94                                      
         RETURN
  END IF
   W( 29 ) = JVS( 1134 )
   W( 30 ) = JVS( 1135 )
   W( 32 ) = JVS( 1136 )
   W( 38 ) = JVS( 1137 )
   W( 39 ) = JVS( 1138 )
   W( 40 ) = JVS( 1139 )
   W( 55 ) = JVS( 1140 )
   W( 56 ) = JVS( 1141 )
   W( 57 ) = JVS( 1142 )
   W( 59 ) = JVS( 1143 )
   W( 61 ) = JVS( 1144 )
   W( 62 ) = JVS( 1145 )
   W( 63 ) = JVS( 1146 )
   W( 64 ) = JVS( 1147 )
   W( 65 ) = JVS( 1148 )
   W( 67 ) = JVS( 1149 )
   W( 68 ) = JVS( 1150 )
   W( 70 ) = JVS( 1151 )
   W( 71 ) = JVS( 1152 )
   W( 72 ) = JVS( 1153 )
   W( 74 ) = JVS( 1154 )
   W( 76 ) = JVS( 1155 )
   W( 77 ) = JVS( 1156 )
   W( 79 ) = JVS( 1157 )
   W( 80 ) = JVS( 1158 )
   W( 81 ) = JVS( 1159 )
   W( 82 ) = JVS( 1160 )
   W( 83 ) = JVS( 1161 )
   W( 84 ) = JVS( 1162 )
   W( 85 ) = JVS( 1163 )
   W( 86 ) = JVS( 1164 )
   W( 87 ) = JVS( 1165 )
   W( 88 ) = JVS( 1166 )
   W( 89 ) = JVS( 1167 )
   W( 90 ) = JVS( 1168 )
   W( 91 ) = JVS( 1169 )
   W( 92 ) = JVS( 1170 )
   W( 93 ) = JVS( 1171 )
   W( 94 ) = JVS( 1172 )
   W( 95 ) = JVS( 1173 )
   W( 96 ) = JVS( 1174 )
   W( 97 ) = JVS( 1175 )
  a = -W( 29 ) / JVS(          119  )
  W( 29 ) = -a
  W( 38 ) = W( 38 ) + a*JVS( 120 )
  W( 39 ) = W( 39 ) + a*JVS( 121 )
  W( 40 ) = W( 40 ) + a*JVS( 122 )
  W( 91 ) = W( 91 ) + a*JVS( 123 )
  W( 94 ) = W( 94 ) + a*JVS( 124 )
  W( 95 ) = W( 95 ) + a*JVS( 125 )
  a = -W( 30 ) / JVS(          127  )
  W( 30 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 128 )
  W( 85 ) = W( 85 ) + a*JVS( 129 )
  W( 94 ) = W( 94 ) + a*JVS( 130 )
  W( 95 ) = W( 95 ) + a*JVS( 131 )
  W( 96 ) = W( 96 ) + a*JVS( 132 )
  a = -W( 32 ) / JVS(          139  )
  W( 32 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 140 )
  W( 88 ) = W( 88 ) + a*JVS( 141 )
  W( 89 ) = W( 89 ) + a*JVS( 142 )
  W( 94 ) = W( 94 ) + a*JVS( 143 )
  a = -W( 38 ) / JVS(          190  )
  W( 38 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 191 )
  W( 88 ) = W( 88 ) + a*JVS( 192 )
  W( 91 ) = W( 91 ) + a*JVS( 193 )
  W( 95 ) = W( 95 ) + a*JVS( 194 )
  a = -W( 39 ) / JVS(          196  )
  W( 39 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 197 )
  W( 88 ) = W( 88 ) + a*JVS( 198 )
  W( 91 ) = W( 91 ) + a*JVS( 199 )
  W( 95 ) = W( 95 ) + a*JVS( 200 )
  a = -W( 40 ) / JVS(          201  )
  W( 40 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 202 )
  W( 85 ) = W( 85 ) + a*JVS( 203 )
  W( 91 ) = W( 91 ) + a*JVS( 204 )
  W( 95 ) = W( 95 ) + a*JVS( 205 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 56 ) / JVS(          312  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 313 )
  W( 64 ) = W( 64 ) + a*JVS( 314 )
  W( 85 ) = W( 85 ) + a*JVS( 315 )
  W( 89 ) = W( 89 ) + a*JVS( 316 )
  W( 91 ) = W( 91 ) + a*JVS( 317 )
  W( 94 ) = W( 94 ) + a*JVS( 318 )
  W( 95 ) = W( 95 ) + a*JVS( 319 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 325 )
  W( 88 ) = W( 88 ) + a*JVS( 326 )
  W( 89 ) = W( 89 ) + a*JVS( 327 )
  W( 91 ) = W( 91 ) + a*JVS( 328 )
  W( 95 ) = W( 95 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  a = -W( 59 ) / JVS(          337  )
  W( 59 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 338 )
  W( 88 ) = W( 88 ) + a*JVS( 339 )
  W( 91 ) = W( 91 ) + a*JVS( 340 )
  W( 96 ) = W( 96 ) + a*JVS( 341 )
  a = -W( 61 ) / JVS(          358  )
  W( 61 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 359 )
  W( 89 ) = W( 89 ) + a*JVS( 360 )
  W( 90 ) = W( 90 ) + a*JVS( 361 )
  W( 91 ) = W( 91 ) + a*JVS( 362 )
  W( 93 ) = W( 93 ) + a*JVS( 363 )
  W( 94 ) = W( 94 ) + a*JVS( 364 )
  W( 96 ) = W( 96 ) + a*JVS( 365 )
  a = -W( 62 ) / JVS(          367  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 368 )
  W( 88 ) = W( 88 ) + a*JVS( 369 )
  W( 89 ) = W( 89 ) + a*JVS( 370 )
  W( 90 ) = W( 90 ) + a*JVS( 371 )
  W( 91 ) = W( 91 ) + a*JVS( 372 )
  W( 93 ) = W( 93 ) + a*JVS( 373 )
  W( 94 ) = W( 94 ) + a*JVS( 374 )
  W( 96 ) = W( 96 ) + a*JVS( 375 )
  a = -W( 63 ) / JVS(          377  )
  W( 63 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  W( 90 ) = W( 90 ) + a*JVS( 380 )
  W( 91 ) = W( 91 ) + a*JVS( 381 )
  W( 93 ) = W( 93 ) + a*JVS( 382 )
  W( 94 ) = W( 94 ) + a*JVS( 383 )
  W( 96 ) = W( 96 ) + a*JVS( 384 )
  a = -W( 64 ) / JVS(          391  )
  W( 64 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 392 )
  W( 85 ) = W( 85 ) + a*JVS( 393 )
  W( 88 ) = W( 88 ) + a*JVS( 394 )
  W( 89 ) = W( 89 ) + a*JVS( 395 )
  W( 90 ) = W( 90 ) + a*JVS( 396 )
  W( 91 ) = W( 91 ) + a*JVS( 397 )
  W( 93 ) = W( 93 ) + a*JVS( 398 )
  W( 94 ) = W( 94 ) + a*JVS( 399 )
  W( 96 ) = W( 96 ) + a*JVS( 400 )
  a = -W( 65 ) / JVS(          402  )
  W( 65 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 403 )
  W( 89 ) = W( 89 ) + a*JVS( 404 )
  W( 90 ) = W( 90 ) + a*JVS( 405 )
  W( 91 ) = W( 91 ) + a*JVS( 406 )
  W( 93 ) = W( 93 ) + a*JVS( 407 )
  W( 94 ) = W( 94 ) + a*JVS( 408 )
  W( 96 ) = W( 96 ) + a*JVS( 409 )
  a = -W( 67 ) / JVS(          440  )
  W( 67 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 441 )
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 89 ) = W( 89 ) + a*JVS( 443 )
  W( 90 ) = W( 90 ) + a*JVS( 444 )
  W( 91 ) = W( 91 ) + a*JVS( 445 )
  W( 93 ) = W( 93 ) + a*JVS( 446 )
  W( 94 ) = W( 94 ) + a*JVS( 447 )
  W( 96 ) = W( 96 ) + a*JVS( 448 )
  a = -W( 68 ) / JVS(          451  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 452 )
  W( 88 ) = W( 88 ) + a*JVS( 453 )
  W( 89 ) = W( 89 ) + a*JVS( 454 )
  W( 90 ) = W( 90 ) + a*JVS( 455 )
  W( 91 ) = W( 91 ) + a*JVS( 456 )
  W( 93 ) = W( 93 ) + a*JVS( 457 )
  W( 94 ) = W( 94 ) + a*JVS( 458 )
  W( 95 ) = W( 95 ) + a*JVS( 459 )
  W( 96 ) = W( 96 ) + a*JVS( 460 )
  a = -W( 70 ) / JVS(          491  )
  W( 70 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 492 )
  W( 85 ) = W( 85 ) + a*JVS( 493 )
  W( 89 ) = W( 89 ) + a*JVS( 494 )
  W( 90 ) = W( 90 ) + a*JVS( 495 )
  W( 91 ) = W( 91 ) + a*JVS( 496 )
  W( 93 ) = W( 93 ) + a*JVS( 497 )
  W( 94 ) = W( 94 ) + a*JVS( 498 )
  W( 96 ) = W( 96 ) + a*JVS( 499 )
  a = -W( 71 ) / JVS(          501  )
  W( 71 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 502 )
  W( 88 ) = W( 88 ) + a*JVS( 503 )
  W( 89 ) = W( 89 ) + a*JVS( 504 )
  W( 90 ) = W( 90 ) + a*JVS( 505 )
  W( 91 ) = W( 91 ) + a*JVS( 506 )
  W( 93 ) = W( 93 ) + a*JVS( 507 )
  W( 94 ) = W( 94 ) + a*JVS( 508 )
  W( 95 ) = W( 95 ) + a*JVS( 509 )
  W( 96 ) = W( 96 ) + a*JVS( 510 )
  a = -W( 72 ) / JVS(          512  )
  W( 72 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 513 )
  W( 88 ) = W( 88 ) + a*JVS( 514 )
  W( 89 ) = W( 89 ) + a*JVS( 515 )
  W( 90 ) = W( 90 ) + a*JVS( 516 )
  W( 91 ) = W( 91 ) + a*JVS( 517 )
  W( 93 ) = W( 93 ) + a*JVS( 518 )
  W( 94 ) = W( 94 ) + a*JVS( 519 )
  W( 95 ) = W( 95 ) + a*JVS( 520 )
  W( 96 ) = W( 96 ) + a*JVS( 521 )
  a = -W( 74 ) / JVS(          572  )
  W( 74 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 573 )
  W( 85 ) = W( 85 ) + a*JVS( 574 )
  W( 89 ) = W( 89 ) + a*JVS( 575 )
  W( 90 ) = W( 90 ) + a*JVS( 576 )
  W( 91 ) = W( 91 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  a = -W( 76 ) / JVS(          607  )
  W( 76 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 608 )
  W( 82 ) = W( 82 ) + a*JVS( 609 )
  W( 88 ) = W( 88 ) + a*JVS( 610 )
  W( 89 ) = W( 89 ) + a*JVS( 611 )
  W( 90 ) = W( 90 ) + a*JVS( 612 )
  W( 91 ) = W( 91 ) + a*JVS( 613 )
  W( 93 ) = W( 93 ) + a*JVS( 614 )
  W( 94 ) = W( 94 ) + a*JVS( 615 )
  W( 96 ) = W( 96 ) + a*JVS( 616 )
  a = -W( 77 ) / JVS(          618  )
  W( 77 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 619 )
  W( 82 ) = W( 82 ) + a*JVS( 620 )
  W( 88 ) = W( 88 ) + a*JVS( 621 )
  W( 89 ) = W( 89 ) + a*JVS( 622 )
  W( 90 ) = W( 90 ) + a*JVS( 623 )
  W( 91 ) = W( 91 ) + a*JVS( 624 )
  W( 93 ) = W( 93 ) + a*JVS( 625 )
  W( 94 ) = W( 94 ) + a*JVS( 626 )
  W( 96 ) = W( 96 ) + a*JVS( 627 )
  a = -W( 79 ) / JVS(          664  )
  W( 79 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 665 )
  W( 88 ) = W( 88 ) + a*JVS( 666 )
  W( 89 ) = W( 89 ) + a*JVS( 667 )
  W( 90 ) = W( 90 ) + a*JVS( 668 )
  W( 91 ) = W( 91 ) + a*JVS( 669 )
  W( 93 ) = W( 93 ) + a*JVS( 670 )
  W( 94 ) = W( 94 ) + a*JVS( 671 )
  W( 95 ) = W( 95 ) + a*JVS( 672 )
  W( 96 ) = W( 96 ) + a*JVS( 673 )
  a = -W( 80 ) / JVS(          692  )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 693 )
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 85 ) = W( 85 ) + a*JVS( 695 )
  W( 86 ) = W( 86 ) + a*JVS( 696 )
  W( 87 ) = W( 87 ) + a*JVS( 697 )
  W( 88 ) = W( 88 ) + a*JVS( 698 )
  W( 89 ) = W( 89 ) + a*JVS( 699 )
  W( 90 ) = W( 90 ) + a*JVS( 700 )
  W( 91 ) = W( 91 ) + a*JVS( 701 )
  W( 92 ) = W( 92 ) + a*JVS( 702 )
  W( 93 ) = W( 93 ) + a*JVS( 703 )
  W( 94 ) = W( 94 ) + a*JVS( 704 )
  W( 95 ) = W( 95 ) + a*JVS( 705 )
  W( 96 ) = W( 96 ) + a*JVS( 706 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  a = -W( 83 ) / JVS(          739  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 740 )
  W( 85 ) = W( 85 ) + a*JVS( 741 )
  W( 88 ) = W( 88 ) + a*JVS( 742 )
  W( 89 ) = W( 89 ) + a*JVS( 743 )
  W( 90 ) = W( 90 ) + a*JVS( 744 )
  W( 91 ) = W( 91 ) + a*JVS( 745 )
  W( 93 ) = W( 93 ) + a*JVS( 746 )
  W( 94 ) = W( 94 ) + a*JVS( 747 )
  W( 96 ) = W( 96 ) + a*JVS( 748 )
  a = -W( 84 ) / JVS(          758  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 759 )
  W( 88 ) = W( 88 ) + a*JVS( 760 )
  W( 89 ) = W( 89 ) + a*JVS( 761 )
  W( 90 ) = W( 90 ) + a*JVS( 762 )
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 93 ) = W( 93 ) + a*JVS( 764 )
  W( 94 ) = W( 94 ) + a*JVS( 765 )
  W( 96 ) = W( 96 ) + a*JVS( 766 )
  a = -W( 85 ) / JVS(          786  )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 787 )
  W( 89 ) = W( 89 ) + a*JVS( 788 )
  W( 90 ) = W( 90 ) + a*JVS( 789 )
  W( 91 ) = W( 91 ) + a*JVS( 790 )
  W( 93 ) = W( 93 ) + a*JVS( 791 )
  W( 94 ) = W( 94 ) + a*JVS( 792 )
  W( 95 ) = W( 95 ) + a*JVS( 793 )
  W( 96 ) = W( 96 ) + a*JVS( 794 )
  a = -W( 86 ) / JVS(          805  )
  W( 86 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 806 )
  W( 89 ) = W( 89 ) + a*JVS( 807 )
  W( 90 ) = W( 90 ) + a*JVS( 808 )
  W( 91 ) = W( 91 ) + a*JVS( 809 )
  W( 92 ) = W( 92 ) + a*JVS( 810 )
  W( 93 ) = W( 93 ) + a*JVS( 811 )
  W( 94 ) = W( 94 ) + a*JVS( 812 )
  W( 95 ) = W( 95 ) + a*JVS( 813 )
  W( 96 ) = W( 96 ) + a*JVS( 814 )
  a = -W( 87 ) / JVS(          831  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 832 )
  W( 89 ) = W( 89 ) + a*JVS( 833 )
  W( 90 ) = W( 90 ) + a*JVS( 834 )
  W( 91 ) = W( 91 ) + a*JVS( 835 )
  W( 92 ) = W( 92 ) + a*JVS( 836 )
  W( 93 ) = W( 93 ) + a*JVS( 837 )
  W( 94 ) = W( 94 ) + a*JVS( 838 )
  W( 95 ) = W( 95 ) + a*JVS( 839 )
  W( 96 ) = W( 96 ) + a*JVS( 840 )
  a = -W( 88 ) / JVS(          871  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 872 )
  W( 90 ) = W( 90 ) + a*JVS( 873 )
  W( 91 ) = W( 91 ) + a*JVS( 874 )
  W( 92 ) = W( 92 ) + a*JVS( 875 )
  W( 93 ) = W( 93 ) + a*JVS( 876 )
  W( 94 ) = W( 94 ) + a*JVS( 877 )
  W( 95 ) = W( 95 ) + a*JVS( 878 )
  W( 96 ) = W( 96 ) + a*JVS( 879 )
  W( 97 ) = W( 97 ) + a*JVS( 880 )
  a = -W( 89 ) / JVS(          940  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 941 )
  W( 91 ) = W( 91 ) + a*JVS( 942 )
  W( 92 ) = W( 92 ) + a*JVS( 943 )
  W( 93 ) = W( 93 ) + a*JVS( 944 )
  W( 94 ) = W( 94 ) + a*JVS( 945 )
  W( 95 ) = W( 95 ) + a*JVS( 946 )
  W( 96 ) = W( 96 ) + a*JVS( 947 )
  W( 97 ) = W( 97 ) + a*JVS( 948 )
  a = -W( 90 ) / JVS(          979  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 980 )
  W( 92 ) = W( 92 ) + a*JVS( 981 )
  W( 93 ) = W( 93 ) + a*JVS( 982 )
  W( 94 ) = W( 94 ) + a*JVS( 983 )
  W( 95 ) = W( 95 ) + a*JVS( 984 )
  W( 96 ) = W( 96 ) + a*JVS( 985 )
  W( 97 ) = W( 97 ) + a*JVS( 986 )
  a = -W( 91 ) / JVS(         1055  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 1056 )
  W( 93 ) = W( 93 ) + a*JVS( 1057 )
  W( 94 ) = W( 94 ) + a*JVS( 1058 )
  W( 95 ) = W( 95 ) + a*JVS( 1059 )
  W( 96 ) = W( 96 ) + a*JVS( 1060 )
  W( 97 ) = W( 97 ) + a*JVS( 1061 )
  a = -W( 92 ) / JVS(         1084  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 1085 )
  W( 94 ) = W( 94 ) + a*JVS( 1086 )
  W( 95 ) = W( 95 ) + a*JVS( 1087 )
  W( 96 ) = W( 96 ) + a*JVS( 1088 )
  W( 97 ) = W( 97 ) + a*JVS( 1089 )
  a = -W( 93 ) / JVS(         1129  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 1130 )
  W( 95 ) = W( 95 ) + a*JVS( 1131 )
  W( 96 ) = W( 96 ) + a*JVS( 1132 )
  W( 97 ) = W( 97 ) + a*JVS( 1133 )
  JVS( 1134) = W( 29 )
  JVS( 1135) = W( 30 )
  JVS( 1136) = W( 32 )
  JVS( 1137) = W( 38 )
  JVS( 1138) = W( 39 )
  JVS( 1139) = W( 40 )
  JVS( 1140) = W( 55 )
  JVS( 1141) = W( 56 )
  JVS( 1142) = W( 57 )
  JVS( 1143) = W( 59 )
  JVS( 1144) = W( 61 )
  JVS( 1145) = W( 62 )
  JVS( 1146) = W( 63 )
  JVS( 1147) = W( 64 )
  JVS( 1148) = W( 65 )
  JVS( 1149) = W( 67 )
  JVS( 1150) = W( 68 )
  JVS( 1151) = W( 70 )
  JVS( 1152) = W( 71 )
  JVS( 1153) = W( 72 )
  JVS( 1154) = W( 74 )
  JVS( 1155) = W( 76 )
  JVS( 1156) = W( 77 )
  JVS( 1157) = W( 79 )
  JVS( 1158) = W( 80 )
  JVS( 1159) = W( 81 )
  JVS( 1160) = W( 82 )
  JVS( 1161) = W( 83 )
  JVS( 1162) = W( 84 )
  JVS( 1163) = W( 85 )
  JVS( 1164) = W( 86 )
  JVS( 1165) = W( 87 )
  JVS( 1166) = W( 88 )
  JVS( 1167) = W( 89 )
  JVS( 1168) = W( 90 )
  JVS( 1169) = W( 91 )
  JVS( 1170) = W( 92 )
  JVS( 1171) = W( 93 )
  JVS( 1172) = W( 94 )
  JVS( 1173) = W( 95 )
  JVS( 1174) = W( 96 )
  JVS( 1175) = W( 97 )
  IF ( ABS(  JVS( 1230 )) < TINY(a) ) THEN
         IER = 95                                      
         RETURN
  END IF
   W( 17 ) = JVS( 1176 )
   W( 18 ) = JVS( 1177 )
   W( 22 ) = JVS( 1178 )
   W( 26 ) = JVS( 1179 )
   W( 27 ) = JVS( 1180 )
   W( 29 ) = JVS( 1181 )
   W( 30 ) = JVS( 1182 )
   W( 31 ) = JVS( 1183 )
   W( 32 ) = JVS( 1184 )
   W( 38 ) = JVS( 1185 )
   W( 39 ) = JVS( 1186 )
   W( 40 ) = JVS( 1187 )
   W( 42 ) = JVS( 1188 )
   W( 44 ) = JVS( 1189 )
   W( 48 ) = JVS( 1190 )
   W( 55 ) = JVS( 1191 )
   W( 56 ) = JVS( 1192 )
   W( 57 ) = JVS( 1193 )
   W( 58 ) = JVS( 1194 )
   W( 59 ) = JVS( 1195 )
   W( 60 ) = JVS( 1196 )
   W( 61 ) = JVS( 1197 )
   W( 62 ) = JVS( 1198 )
   W( 63 ) = JVS( 1199 )
   W( 64 ) = JVS( 1200 )
   W( 65 ) = JVS( 1201 )
   W( 67 ) = JVS( 1202 )
   W( 68 ) = JVS( 1203 )
   W( 69 ) = JVS( 1204 )
   W( 70 ) = JVS( 1205 )
   W( 71 ) = JVS( 1206 )
   W( 72 ) = JVS( 1207 )
   W( 73 ) = JVS( 1208 )
   W( 74 ) = JVS( 1209 )
   W( 75 ) = JVS( 1210 )
   W( 76 ) = JVS( 1211 )
   W( 77 ) = JVS( 1212 )
   W( 78 ) = JVS( 1213 )
   W( 79 ) = JVS( 1214 )
   W( 80 ) = JVS( 1215 )
   W( 81 ) = JVS( 1216 )
   W( 82 ) = JVS( 1217 )
   W( 83 ) = JVS( 1218 )
   W( 84 ) = JVS( 1219 )
   W( 85 ) = JVS( 1220 )
   W( 86 ) = JVS( 1221 )
   W( 87 ) = JVS( 1222 )
   W( 88 ) = JVS( 1223 )
   W( 89 ) = JVS( 1224 )
   W( 90 ) = JVS( 1225 )
   W( 91 ) = JVS( 1226 )
   W( 92 ) = JVS( 1227 )
   W( 93 ) = JVS( 1228 )
   W( 94 ) = JVS( 1229 )
   W( 95 ) = JVS( 1230 )
   W( 96 ) = JVS( 1231 )
   W( 97 ) = JVS( 1232 )
  a = -W( 17 ) / JVS(           75  )
  W( 17 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 76 )
  W( 96 ) = W( 96 ) + a*JVS( 77 )
  a = -W( 18 ) / JVS(           78  )
  W( 18 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 79 )
  W( 95 ) = W( 95 ) + a*JVS( 80 )
  a = -W( 22 ) / JVS(           91  )
  W( 22 ) = -a
  W( 42 ) = W( 42 ) + a*JVS( 92 )
  W( 91 ) = W( 91 ) + a*JVS( 93 )
  a = -W( 26 ) / JVS(          104  )
  W( 26 ) = -a
  W( 56 ) = W( 56 ) + a*JVS( 105 )
  W( 91 ) = W( 91 ) + a*JVS( 106 )
  W( 95 ) = W( 95 ) + a*JVS( 107 )
  a = -W( 27 ) / JVS(          108  )
  W( 27 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 109 )
  W( 91 ) = W( 91 ) + a*JVS( 110 )
  W( 95 ) = W( 95 ) + a*JVS( 111 )
  a = -W( 29 ) / JVS(          119  )
  W( 29 ) = -a
  W( 38 ) = W( 38 ) + a*JVS( 120 )
  W( 39 ) = W( 39 ) + a*JVS( 121 )
  W( 40 ) = W( 40 ) + a*JVS( 122 )
  W( 91 ) = W( 91 ) + a*JVS( 123 )
  W( 94 ) = W( 94 ) + a*JVS( 124 )
  W( 95 ) = W( 95 ) + a*JVS( 125 )
  a = -W( 30 ) / JVS(          127  )
  W( 30 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 128 )
  W( 85 ) = W( 85 ) + a*JVS( 129 )
  W( 94 ) = W( 94 ) + a*JVS( 130 )
  W( 95 ) = W( 95 ) + a*JVS( 131 )
  W( 96 ) = W( 96 ) + a*JVS( 132 )
  a = -W( 31 ) / JVS(          133  )
  W( 31 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 134 )
  W( 89 ) = W( 89 ) + a*JVS( 135 )
  W( 91 ) = W( 91 ) + a*JVS( 136 )
  W( 95 ) = W( 95 ) + a*JVS( 137 )
  W( 96 ) = W( 96 ) + a*JVS( 138 )
  a = -W( 32 ) / JVS(          139  )
  W( 32 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 140 )
  W( 88 ) = W( 88 ) + a*JVS( 141 )
  W( 89 ) = W( 89 ) + a*JVS( 142 )
  W( 94 ) = W( 94 ) + a*JVS( 143 )
  a = -W( 38 ) / JVS(          190  )
  W( 38 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 191 )
  W( 88 ) = W( 88 ) + a*JVS( 192 )
  W( 91 ) = W( 91 ) + a*JVS( 193 )
  W( 95 ) = W( 95 ) + a*JVS( 194 )
  a = -W( 39 ) / JVS(          196  )
  W( 39 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 197 )
  W( 88 ) = W( 88 ) + a*JVS( 198 )
  W( 91 ) = W( 91 ) + a*JVS( 199 )
  W( 95 ) = W( 95 ) + a*JVS( 200 )
  a = -W( 40 ) / JVS(          201  )
  W( 40 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 202 )
  W( 85 ) = W( 85 ) + a*JVS( 203 )
  W( 91 ) = W( 91 ) + a*JVS( 204 )
  W( 95 ) = W( 95 ) + a*JVS( 205 )
  a = -W( 42 ) / JVS(          211  )
  W( 42 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 212 )
  W( 67 ) = W( 67 ) + a*JVS( 213 )
  W( 91 ) = W( 91 ) + a*JVS( 214 )
  W( 94 ) = W( 94 ) + a*JVS( 215 )
  W( 96 ) = W( 96 ) + a*JVS( 216 )
  a = -W( 44 ) / JVS(          221  )
  W( 44 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 222 )
  W( 60 ) = W( 60 ) + a*JVS( 223 )
  W( 69 ) = W( 69 ) + a*JVS( 224 )
  W( 73 ) = W( 73 ) + a*JVS( 225 )
  W( 78 ) = W( 78 ) + a*JVS( 226 )
  W( 79 ) = W( 79 ) + a*JVS( 227 )
  W( 89 ) = W( 89 ) + a*JVS( 228 )
  W( 91 ) = W( 91 ) + a*JVS( 229 )
  W( 95 ) = W( 95 ) + a*JVS( 230 )
  W( 96 ) = W( 96 ) + a*JVS( 231 )
  a = -W( 48 ) / JVS(          247  )
  W( 48 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 248 )
  W( 85 ) = W( 85 ) + a*JVS( 249 )
  W( 91 ) = W( 91 ) + a*JVS( 250 )
  W( 93 ) = W( 93 ) + a*JVS( 251 )
  W( 95 ) = W( 95 ) + a*JVS( 252 )
  W( 96 ) = W( 96 ) + a*JVS( 253 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 56 ) / JVS(          312  )
  W( 56 ) = -a
  W( 59 ) = W( 59 ) + a*JVS( 313 )
  W( 64 ) = W( 64 ) + a*JVS( 314 )
  W( 85 ) = W( 85 ) + a*JVS( 315 )
  W( 89 ) = W( 89 ) + a*JVS( 316 )
  W( 91 ) = W( 91 ) + a*JVS( 317 )
  W( 94 ) = W( 94 ) + a*JVS( 318 )
  W( 95 ) = W( 95 ) + a*JVS( 319 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 325 )
  W( 88 ) = W( 88 ) + a*JVS( 326 )
  W( 89 ) = W( 89 ) + a*JVS( 327 )
  W( 91 ) = W( 91 ) + a*JVS( 328 )
  W( 95 ) = W( 95 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  a = -W( 58 ) / JVS(          331  )
  W( 58 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 332 )
  W( 85 ) = W( 85 ) + a*JVS( 333 )
  W( 91 ) = W( 91 ) + a*JVS( 334 )
  W( 95 ) = W( 95 ) + a*JVS( 335 )
  W( 96 ) = W( 96 ) + a*JVS( 336 )
  a = -W( 59 ) / JVS(          337  )
  W( 59 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 338 )
  W( 88 ) = W( 88 ) + a*JVS( 339 )
  W( 91 ) = W( 91 ) + a*JVS( 340 )
  W( 96 ) = W( 96 ) + a*JVS( 341 )
  a = -W( 60 ) / JVS(          343  )
  W( 60 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 344 )
  W( 71 ) = W( 71 ) + a*JVS( 345 )
  W( 72 ) = W( 72 ) + a*JVS( 346 )
  W( 79 ) = W( 79 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 82 ) = W( 82 ) + a*JVS( 349 )
  W( 85 ) = W( 85 ) + a*JVS( 350 )
  W( 86 ) = W( 86 ) + a*JVS( 351 )
  W( 88 ) = W( 88 ) + a*JVS( 352 )
  W( 90 ) = W( 90 ) + a*JVS( 353 )
  W( 91 ) = W( 91 ) + a*JVS( 354 )
  W( 93 ) = W( 93 ) + a*JVS( 355 )
  W( 94 ) = W( 94 ) + a*JVS( 356 )
  W( 96 ) = W( 96 ) + a*JVS( 357 )
  a = -W( 61 ) / JVS(          358  )
  W( 61 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 359 )
  W( 89 ) = W( 89 ) + a*JVS( 360 )
  W( 90 ) = W( 90 ) + a*JVS( 361 )
  W( 91 ) = W( 91 ) + a*JVS( 362 )
  W( 93 ) = W( 93 ) + a*JVS( 363 )
  W( 94 ) = W( 94 ) + a*JVS( 364 )
  W( 96 ) = W( 96 ) + a*JVS( 365 )
  a = -W( 62 ) / JVS(          367  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 368 )
  W( 88 ) = W( 88 ) + a*JVS( 369 )
  W( 89 ) = W( 89 ) + a*JVS( 370 )
  W( 90 ) = W( 90 ) + a*JVS( 371 )
  W( 91 ) = W( 91 ) + a*JVS( 372 )
  W( 93 ) = W( 93 ) + a*JVS( 373 )
  W( 94 ) = W( 94 ) + a*JVS( 374 )
  W( 96 ) = W( 96 ) + a*JVS( 375 )
  a = -W( 63 ) / JVS(          377  )
  W( 63 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  W( 90 ) = W( 90 ) + a*JVS( 380 )
  W( 91 ) = W( 91 ) + a*JVS( 381 )
  W( 93 ) = W( 93 ) + a*JVS( 382 )
  W( 94 ) = W( 94 ) + a*JVS( 383 )
  W( 96 ) = W( 96 ) + a*JVS( 384 )
  a = -W( 64 ) / JVS(          391  )
  W( 64 ) = -a
  W( 67 ) = W( 67 ) + a*JVS( 392 )
  W( 85 ) = W( 85 ) + a*JVS( 393 )
  W( 88 ) = W( 88 ) + a*JVS( 394 )
  W( 89 ) = W( 89 ) + a*JVS( 395 )
  W( 90 ) = W( 90 ) + a*JVS( 396 )
  W( 91 ) = W( 91 ) + a*JVS( 397 )
  W( 93 ) = W( 93 ) + a*JVS( 398 )
  W( 94 ) = W( 94 ) + a*JVS( 399 )
  W( 96 ) = W( 96 ) + a*JVS( 400 )
  a = -W( 65 ) / JVS(          402  )
  W( 65 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 403 )
  W( 89 ) = W( 89 ) + a*JVS( 404 )
  W( 90 ) = W( 90 ) + a*JVS( 405 )
  W( 91 ) = W( 91 ) + a*JVS( 406 )
  W( 93 ) = W( 93 ) + a*JVS( 407 )
  W( 94 ) = W( 94 ) + a*JVS( 408 )
  W( 96 ) = W( 96 ) + a*JVS( 409 )
  a = -W( 67 ) / JVS(          440  )
  W( 67 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 441 )
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 89 ) = W( 89 ) + a*JVS( 443 )
  W( 90 ) = W( 90 ) + a*JVS( 444 )
  W( 91 ) = W( 91 ) + a*JVS( 445 )
  W( 93 ) = W( 93 ) + a*JVS( 446 )
  W( 94 ) = W( 94 ) + a*JVS( 447 )
  W( 96 ) = W( 96 ) + a*JVS( 448 )
  a = -W( 68 ) / JVS(          451  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 452 )
  W( 88 ) = W( 88 ) + a*JVS( 453 )
  W( 89 ) = W( 89 ) + a*JVS( 454 )
  W( 90 ) = W( 90 ) + a*JVS( 455 )
  W( 91 ) = W( 91 ) + a*JVS( 456 )
  W( 93 ) = W( 93 ) + a*JVS( 457 )
  W( 94 ) = W( 94 ) + a*JVS( 458 )
  W( 95 ) = W( 95 ) + a*JVS( 459 )
  W( 96 ) = W( 96 ) + a*JVS( 460 )
  a = -W( 69 ) / JVS(          471  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 472 )
  W( 72 ) = W( 72 ) + a*JVS( 473 )
  W( 74 ) = W( 74 ) + a*JVS( 474 )
  W( 79 ) = W( 79 ) + a*JVS( 475 )
  W( 81 ) = W( 81 ) + a*JVS( 476 )
  W( 82 ) = W( 82 ) + a*JVS( 477 )
  W( 85 ) = W( 85 ) + a*JVS( 478 )
  W( 86 ) = W( 86 ) + a*JVS( 479 )
  W( 87 ) = W( 87 ) + a*JVS( 480 )
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 89 ) = W( 89 ) + a*JVS( 482 )
  W( 90 ) = W( 90 ) + a*JVS( 483 )
  W( 91 ) = W( 91 ) + a*JVS( 484 )
  W( 93 ) = W( 93 ) + a*JVS( 485 )
  W( 94 ) = W( 94 ) + a*JVS( 486 )
  W( 95 ) = W( 95 ) + a*JVS( 487 )
  W( 96 ) = W( 96 ) + a*JVS( 488 )
  a = -W( 70 ) / JVS(          491  )
  W( 70 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 492 )
  W( 85 ) = W( 85 ) + a*JVS( 493 )
  W( 89 ) = W( 89 ) + a*JVS( 494 )
  W( 90 ) = W( 90 ) + a*JVS( 495 )
  W( 91 ) = W( 91 ) + a*JVS( 496 )
  W( 93 ) = W( 93 ) + a*JVS( 497 )
  W( 94 ) = W( 94 ) + a*JVS( 498 )
  W( 96 ) = W( 96 ) + a*JVS( 499 )
  a = -W( 71 ) / JVS(          501  )
  W( 71 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 502 )
  W( 88 ) = W( 88 ) + a*JVS( 503 )
  W( 89 ) = W( 89 ) + a*JVS( 504 )
  W( 90 ) = W( 90 ) + a*JVS( 505 )
  W( 91 ) = W( 91 ) + a*JVS( 506 )
  W( 93 ) = W( 93 ) + a*JVS( 507 )
  W( 94 ) = W( 94 ) + a*JVS( 508 )
  W( 95 ) = W( 95 ) + a*JVS( 509 )
  W( 96 ) = W( 96 ) + a*JVS( 510 )
  a = -W( 72 ) / JVS(          512  )
  W( 72 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 513 )
  W( 88 ) = W( 88 ) + a*JVS( 514 )
  W( 89 ) = W( 89 ) + a*JVS( 515 )
  W( 90 ) = W( 90 ) + a*JVS( 516 )
  W( 91 ) = W( 91 ) + a*JVS( 517 )
  W( 93 ) = W( 93 ) + a*JVS( 518 )
  W( 94 ) = W( 94 ) + a*JVS( 519 )
  W( 95 ) = W( 95 ) + a*JVS( 520 )
  W( 96 ) = W( 96 ) + a*JVS( 521 )
  a = -W( 73 ) / JVS(          549  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 550 )
  W( 76 ) = W( 76 ) + a*JVS( 551 )
  W( 77 ) = W( 77 ) + a*JVS( 552 )
  W( 79 ) = W( 79 ) + a*JVS( 553 )
  W( 80 ) = W( 80 ) + a*JVS( 554 )
  W( 81 ) = W( 81 ) + a*JVS( 555 )
  W( 82 ) = W( 82 ) + a*JVS( 556 )
  W( 83 ) = W( 83 ) + a*JVS( 557 )
  W( 84 ) = W( 84 ) + a*JVS( 558 )
  W( 85 ) = W( 85 ) + a*JVS( 559 )
  W( 86 ) = W( 86 ) + a*JVS( 560 )
  W( 87 ) = W( 87 ) + a*JVS( 561 )
  W( 88 ) = W( 88 ) + a*JVS( 562 )
  W( 89 ) = W( 89 ) + a*JVS( 563 )
  W( 90 ) = W( 90 ) + a*JVS( 564 )
  W( 91 ) = W( 91 ) + a*JVS( 565 )
  W( 93 ) = W( 93 ) + a*JVS( 566 )
  W( 94 ) = W( 94 ) + a*JVS( 567 )
  W( 95 ) = W( 95 ) + a*JVS( 568 )
  W( 96 ) = W( 96 ) + a*JVS( 569 )
  W( 97 ) = W( 97 ) + a*JVS( 570 )
  a = -W( 74 ) / JVS(          572  )
  W( 74 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 573 )
  W( 85 ) = W( 85 ) + a*JVS( 574 )
  W( 89 ) = W( 89 ) + a*JVS( 575 )
  W( 90 ) = W( 90 ) + a*JVS( 576 )
  W( 91 ) = W( 91 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  a = -W( 75 ) / JVS(          590  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 591 )
  W( 77 ) = W( 77 ) + a*JVS( 592 )
  W( 79 ) = W( 79 ) + a*JVS( 593 )
  W( 83 ) = W( 83 ) + a*JVS( 594 )
  W( 84 ) = W( 84 ) + a*JVS( 595 )
  W( 85 ) = W( 85 ) + a*JVS( 596 )
  W( 86 ) = W( 86 ) + a*JVS( 597 )
  W( 88 ) = W( 88 ) + a*JVS( 598 )
  W( 89 ) = W( 89 ) + a*JVS( 599 )
  W( 90 ) = W( 90 ) + a*JVS( 600 )
  W( 91 ) = W( 91 ) + a*JVS( 601 )
  W( 93 ) = W( 93 ) + a*JVS( 602 )
  W( 94 ) = W( 94 ) + a*JVS( 603 )
  W( 95 ) = W( 95 ) + a*JVS( 604 )
  W( 96 ) = W( 96 ) + a*JVS( 605 )
  a = -W( 76 ) / JVS(          607  )
  W( 76 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 608 )
  W( 82 ) = W( 82 ) + a*JVS( 609 )
  W( 88 ) = W( 88 ) + a*JVS( 610 )
  W( 89 ) = W( 89 ) + a*JVS( 611 )
  W( 90 ) = W( 90 ) + a*JVS( 612 )
  W( 91 ) = W( 91 ) + a*JVS( 613 )
  W( 93 ) = W( 93 ) + a*JVS( 614 )
  W( 94 ) = W( 94 ) + a*JVS( 615 )
  W( 96 ) = W( 96 ) + a*JVS( 616 )
  a = -W( 77 ) / JVS(          618  )
  W( 77 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 619 )
  W( 82 ) = W( 82 ) + a*JVS( 620 )
  W( 88 ) = W( 88 ) + a*JVS( 621 )
  W( 89 ) = W( 89 ) + a*JVS( 622 )
  W( 90 ) = W( 90 ) + a*JVS( 623 )
  W( 91 ) = W( 91 ) + a*JVS( 624 )
  W( 93 ) = W( 93 ) + a*JVS( 625 )
  W( 94 ) = W( 94 ) + a*JVS( 626 )
  W( 96 ) = W( 96 ) + a*JVS( 627 )
  a = -W( 78 ) / JVS(          643  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 644 )
  W( 81 ) = W( 81 ) + a*JVS( 645 )
  W( 82 ) = W( 82 ) + a*JVS( 646 )
  W( 83 ) = W( 83 ) + a*JVS( 647 )
  W( 84 ) = W( 84 ) + a*JVS( 648 )
  W( 85 ) = W( 85 ) + a*JVS( 649 )
  W( 86 ) = W( 86 ) + a*JVS( 650 )
  W( 87 ) = W( 87 ) + a*JVS( 651 )
  W( 88 ) = W( 88 ) + a*JVS( 652 )
  W( 89 ) = W( 89 ) + a*JVS( 653 )
  W( 90 ) = W( 90 ) + a*JVS( 654 )
  W( 91 ) = W( 91 ) + a*JVS( 655 )
  W( 92 ) = W( 92 ) + a*JVS( 656 )
  W( 93 ) = W( 93 ) + a*JVS( 657 )
  W( 94 ) = W( 94 ) + a*JVS( 658 )
  W( 95 ) = W( 95 ) + a*JVS( 659 )
  W( 96 ) = W( 96 ) + a*JVS( 660 )
  W( 97 ) = W( 97 ) + a*JVS( 661 )
  a = -W( 79 ) / JVS(          664  )
  W( 79 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 665 )
  W( 88 ) = W( 88 ) + a*JVS( 666 )
  W( 89 ) = W( 89 ) + a*JVS( 667 )
  W( 90 ) = W( 90 ) + a*JVS( 668 )
  W( 91 ) = W( 91 ) + a*JVS( 669 )
  W( 93 ) = W( 93 ) + a*JVS( 670 )
  W( 94 ) = W( 94 ) + a*JVS( 671 )
  W( 95 ) = W( 95 ) + a*JVS( 672 )
  W( 96 ) = W( 96 ) + a*JVS( 673 )
  a = -W( 80 ) / JVS(          692  )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 693 )
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 85 ) = W( 85 ) + a*JVS( 695 )
  W( 86 ) = W( 86 ) + a*JVS( 696 )
  W( 87 ) = W( 87 ) + a*JVS( 697 )
  W( 88 ) = W( 88 ) + a*JVS( 698 )
  W( 89 ) = W( 89 ) + a*JVS( 699 )
  W( 90 ) = W( 90 ) + a*JVS( 700 )
  W( 91 ) = W( 91 ) + a*JVS( 701 )
  W( 92 ) = W( 92 ) + a*JVS( 702 )
  W( 93 ) = W( 93 ) + a*JVS( 703 )
  W( 94 ) = W( 94 ) + a*JVS( 704 )
  W( 95 ) = W( 95 ) + a*JVS( 705 )
  W( 96 ) = W( 96 ) + a*JVS( 706 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  a = -W( 83 ) / JVS(          739  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 740 )
  W( 85 ) = W( 85 ) + a*JVS( 741 )
  W( 88 ) = W( 88 ) + a*JVS( 742 )
  W( 89 ) = W( 89 ) + a*JVS( 743 )
  W( 90 ) = W( 90 ) + a*JVS( 744 )
  W( 91 ) = W( 91 ) + a*JVS( 745 )
  W( 93 ) = W( 93 ) + a*JVS( 746 )
  W( 94 ) = W( 94 ) + a*JVS( 747 )
  W( 96 ) = W( 96 ) + a*JVS( 748 )
  a = -W( 84 ) / JVS(          758  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 759 )
  W( 88 ) = W( 88 ) + a*JVS( 760 )
  W( 89 ) = W( 89 ) + a*JVS( 761 )
  W( 90 ) = W( 90 ) + a*JVS( 762 )
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 93 ) = W( 93 ) + a*JVS( 764 )
  W( 94 ) = W( 94 ) + a*JVS( 765 )
  W( 96 ) = W( 96 ) + a*JVS( 766 )
  a = -W( 85 ) / JVS(          786  )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 787 )
  W( 89 ) = W( 89 ) + a*JVS( 788 )
  W( 90 ) = W( 90 ) + a*JVS( 789 )
  W( 91 ) = W( 91 ) + a*JVS( 790 )
  W( 93 ) = W( 93 ) + a*JVS( 791 )
  W( 94 ) = W( 94 ) + a*JVS( 792 )
  W( 95 ) = W( 95 ) + a*JVS( 793 )
  W( 96 ) = W( 96 ) + a*JVS( 794 )
  a = -W( 86 ) / JVS(          805  )
  W( 86 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 806 )
  W( 89 ) = W( 89 ) + a*JVS( 807 )
  W( 90 ) = W( 90 ) + a*JVS( 808 )
  W( 91 ) = W( 91 ) + a*JVS( 809 )
  W( 92 ) = W( 92 ) + a*JVS( 810 )
  W( 93 ) = W( 93 ) + a*JVS( 811 )
  W( 94 ) = W( 94 ) + a*JVS( 812 )
  W( 95 ) = W( 95 ) + a*JVS( 813 )
  W( 96 ) = W( 96 ) + a*JVS( 814 )
  a = -W( 87 ) / JVS(          831  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 832 )
  W( 89 ) = W( 89 ) + a*JVS( 833 )
  W( 90 ) = W( 90 ) + a*JVS( 834 )
  W( 91 ) = W( 91 ) + a*JVS( 835 )
  W( 92 ) = W( 92 ) + a*JVS( 836 )
  W( 93 ) = W( 93 ) + a*JVS( 837 )
  W( 94 ) = W( 94 ) + a*JVS( 838 )
  W( 95 ) = W( 95 ) + a*JVS( 839 )
  W( 96 ) = W( 96 ) + a*JVS( 840 )
  a = -W( 88 ) / JVS(          871  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 872 )
  W( 90 ) = W( 90 ) + a*JVS( 873 )
  W( 91 ) = W( 91 ) + a*JVS( 874 )
  W( 92 ) = W( 92 ) + a*JVS( 875 )
  W( 93 ) = W( 93 ) + a*JVS( 876 )
  W( 94 ) = W( 94 ) + a*JVS( 877 )
  W( 95 ) = W( 95 ) + a*JVS( 878 )
  W( 96 ) = W( 96 ) + a*JVS( 879 )
  W( 97 ) = W( 97 ) + a*JVS( 880 )
  a = -W( 89 ) / JVS(          940  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 941 )
  W( 91 ) = W( 91 ) + a*JVS( 942 )
  W( 92 ) = W( 92 ) + a*JVS( 943 )
  W( 93 ) = W( 93 ) + a*JVS( 944 )
  W( 94 ) = W( 94 ) + a*JVS( 945 )
  W( 95 ) = W( 95 ) + a*JVS( 946 )
  W( 96 ) = W( 96 ) + a*JVS( 947 )
  W( 97 ) = W( 97 ) + a*JVS( 948 )
  a = -W( 90 ) / JVS(          979  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 980 )
  W( 92 ) = W( 92 ) + a*JVS( 981 )
  W( 93 ) = W( 93 ) + a*JVS( 982 )
  W( 94 ) = W( 94 ) + a*JVS( 983 )
  W( 95 ) = W( 95 ) + a*JVS( 984 )
  W( 96 ) = W( 96 ) + a*JVS( 985 )
  W( 97 ) = W( 97 ) + a*JVS( 986 )
  a = -W( 91 ) / JVS(         1055  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 1056 )
  W( 93 ) = W( 93 ) + a*JVS( 1057 )
  W( 94 ) = W( 94 ) + a*JVS( 1058 )
  W( 95 ) = W( 95 ) + a*JVS( 1059 )
  W( 96 ) = W( 96 ) + a*JVS( 1060 )
  W( 97 ) = W( 97 ) + a*JVS( 1061 )
  a = -W( 92 ) / JVS(         1084  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 1085 )
  W( 94 ) = W( 94 ) + a*JVS( 1086 )
  W( 95 ) = W( 95 ) + a*JVS( 1087 )
  W( 96 ) = W( 96 ) + a*JVS( 1088 )
  W( 97 ) = W( 97 ) + a*JVS( 1089 )
  a = -W( 93 ) / JVS(         1129  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 1130 )
  W( 95 ) = W( 95 ) + a*JVS( 1131 )
  W( 96 ) = W( 96 ) + a*JVS( 1132 )
  W( 97 ) = W( 97 ) + a*JVS( 1133 )
  a = -W( 94 ) / JVS(         1172  )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 1173 )
  W( 96 ) = W( 96 ) + a*JVS( 1174 )
  W( 97 ) = W( 97 ) + a*JVS( 1175 )
  JVS( 1176) = W( 17 )
  JVS( 1177) = W( 18 )
  JVS( 1178) = W( 22 )
  JVS( 1179) = W( 26 )
  JVS( 1180) = W( 27 )
  JVS( 1181) = W( 29 )
  JVS( 1182) = W( 30 )
  JVS( 1183) = W( 31 )
  JVS( 1184) = W( 32 )
  JVS( 1185) = W( 38 )
  JVS( 1186) = W( 39 )
  JVS( 1187) = W( 40 )
  JVS( 1188) = W( 42 )
  JVS( 1189) = W( 44 )
  JVS( 1190) = W( 48 )
  JVS( 1191) = W( 55 )
  JVS( 1192) = W( 56 )
  JVS( 1193) = W( 57 )
  JVS( 1194) = W( 58 )
  JVS( 1195) = W( 59 )
  JVS( 1196) = W( 60 )
  JVS( 1197) = W( 61 )
  JVS( 1198) = W( 62 )
  JVS( 1199) = W( 63 )
  JVS( 1200) = W( 64 )
  JVS( 1201) = W( 65 )
  JVS( 1202) = W( 67 )
  JVS( 1203) = W( 68 )
  JVS( 1204) = W( 69 )
  JVS( 1205) = W( 70 )
  JVS( 1206) = W( 71 )
  JVS( 1207) = W( 72 )
  JVS( 1208) = W( 73 )
  JVS( 1209) = W( 74 )
  JVS( 1210) = W( 75 )
  JVS( 1211) = W( 76 )
  JVS( 1212) = W( 77 )
  JVS( 1213) = W( 78 )
  JVS( 1214) = W( 79 )
  JVS( 1215) = W( 80 )
  JVS( 1216) = W( 81 )
  JVS( 1217) = W( 82 )
  JVS( 1218) = W( 83 )
  JVS( 1219) = W( 84 )
  JVS( 1220) = W( 85 )
  JVS( 1221) = W( 86 )
  JVS( 1222) = W( 87 )
  JVS( 1223) = W( 88 )
  JVS( 1224) = W( 89 )
  JVS( 1225) = W( 90 )
  JVS( 1226) = W( 91 )
  JVS( 1227) = W( 92 )
  JVS( 1228) = W( 93 )
  JVS( 1229) = W( 94 )
  JVS( 1230) = W( 95 )
  JVS( 1231) = W( 96 )
  JVS( 1232) = W( 97 )
  IF ( ABS(  JVS( 1279 )) < TINY(a) ) THEN
         IER = 96                                      
         RETURN
  END IF
   W( 17 ) = JVS( 1233 )
   W( 27 ) = JVS( 1234 )
   W( 30 ) = JVS( 1235 )
   W( 43 ) = JVS( 1236 )
   W( 44 ) = JVS( 1237 )
   W( 47 ) = JVS( 1238 )
   W( 48 ) = JVS( 1239 )
   W( 51 ) = JVS( 1240 )
   W( 53 ) = JVS( 1241 )
   W( 54 ) = JVS( 1242 )
   W( 55 ) = JVS( 1243 )
   W( 57 ) = JVS( 1244 )
   W( 58 ) = JVS( 1245 )
   W( 59 ) = JVS( 1246 )
   W( 60 ) = JVS( 1247 )
   W( 61 ) = JVS( 1248 )
   W( 62 ) = JVS( 1249 )
   W( 63 ) = JVS( 1250 )
   W( 65 ) = JVS( 1251 )
   W( 68 ) = JVS( 1252 )
   W( 69 ) = JVS( 1253 )
   W( 70 ) = JVS( 1254 )
   W( 71 ) = JVS( 1255 )
   W( 72 ) = JVS( 1256 )
   W( 73 ) = JVS( 1257 )
   W( 74 ) = JVS( 1258 )
   W( 76 ) = JVS( 1259 )
   W( 77 ) = JVS( 1260 )
   W( 78 ) = JVS( 1261 )
   W( 79 ) = JVS( 1262 )
   W( 80 ) = JVS( 1263 )
   W( 81 ) = JVS( 1264 )
   W( 82 ) = JVS( 1265 )
   W( 83 ) = JVS( 1266 )
   W( 84 ) = JVS( 1267 )
   W( 85 ) = JVS( 1268 )
   W( 86 ) = JVS( 1269 )
   W( 87 ) = JVS( 1270 )
   W( 88 ) = JVS( 1271 )
   W( 89 ) = JVS( 1272 )
   W( 90 ) = JVS( 1273 )
   W( 91 ) = JVS( 1274 )
   W( 92 ) = JVS( 1275 )
   W( 93 ) = JVS( 1276 )
   W( 94 ) = JVS( 1277 )
   W( 95 ) = JVS( 1278 )
   W( 96 ) = JVS( 1279 )
   W( 97 ) = JVS( 1280 )
  a = -W( 17 ) / JVS(           75  )
  W( 17 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 76 )
  W( 96 ) = W( 96 ) + a*JVS( 77 )
  a = -W( 27 ) / JVS(          108  )
  W( 27 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 109 )
  W( 91 ) = W( 91 ) + a*JVS( 110 )
  W( 95 ) = W( 95 ) + a*JVS( 111 )
  a = -W( 30 ) / JVS(          127  )
  W( 30 ) = -a
  W( 55 ) = W( 55 ) + a*JVS( 128 )
  W( 85 ) = W( 85 ) + a*JVS( 129 )
  W( 94 ) = W( 94 ) + a*JVS( 130 )
  W( 95 ) = W( 95 ) + a*JVS( 131 )
  W( 96 ) = W( 96 ) + a*JVS( 132 )
  a = -W( 43 ) / JVS(          217  )
  W( 43 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 218 )
  W( 91 ) = W( 91 ) + a*JVS( 219 )
  W( 96 ) = W( 96 ) + a*JVS( 220 )
  a = -W( 44 ) / JVS(          221  )
  W( 44 ) = -a
  W( 57 ) = W( 57 ) + a*JVS( 222 )
  W( 60 ) = W( 60 ) + a*JVS( 223 )
  W( 69 ) = W( 69 ) + a*JVS( 224 )
  W( 73 ) = W( 73 ) + a*JVS( 225 )
  W( 78 ) = W( 78 ) + a*JVS( 226 )
  W( 79 ) = W( 79 ) + a*JVS( 227 )
  W( 89 ) = W( 89 ) + a*JVS( 228 )
  W( 91 ) = W( 91 ) + a*JVS( 229 )
  W( 95 ) = W( 95 ) + a*JVS( 230 )
  W( 96 ) = W( 96 ) + a*JVS( 231 )
  a = -W( 47 ) / JVS(          243  )
  W( 47 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 244 )
  W( 91 ) = W( 91 ) + a*JVS( 245 )
  W( 96 ) = W( 96 ) + a*JVS( 246 )
  a = -W( 48 ) / JVS(          247  )
  W( 48 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 248 )
  W( 85 ) = W( 85 ) + a*JVS( 249 )
  W( 91 ) = W( 91 ) + a*JVS( 250 )
  W( 93 ) = W( 93 ) + a*JVS( 251 )
  W( 95 ) = W( 95 ) + a*JVS( 252 )
  W( 96 ) = W( 96 ) + a*JVS( 253 )
  a = -W( 51 ) / JVS(          266  )
  W( 51 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 267 )
  W( 91 ) = W( 91 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 53 ) / JVS(          297  )
  W( 53 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 298 )
  W( 88 ) = W( 88 ) + a*JVS( 299 )
  W( 91 ) = W( 91 ) + a*JVS( 300 )
  W( 96 ) = W( 96 ) + a*JVS( 301 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  a = -W( 55 ) / JVS(          306  )
  W( 55 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 307 )
  W( 91 ) = W( 91 ) + a*JVS( 308 )
  W( 96 ) = W( 96 ) + a*JVS( 309 )
  a = -W( 57 ) / JVS(          324  )
  W( 57 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 325 )
  W( 88 ) = W( 88 ) + a*JVS( 326 )
  W( 89 ) = W( 89 ) + a*JVS( 327 )
  W( 91 ) = W( 91 ) + a*JVS( 328 )
  W( 95 ) = W( 95 ) + a*JVS( 329 )
  W( 96 ) = W( 96 ) + a*JVS( 330 )
  a = -W( 58 ) / JVS(          331  )
  W( 58 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 332 )
  W( 85 ) = W( 85 ) + a*JVS( 333 )
  W( 91 ) = W( 91 ) + a*JVS( 334 )
  W( 95 ) = W( 95 ) + a*JVS( 335 )
  W( 96 ) = W( 96 ) + a*JVS( 336 )
  a = -W( 59 ) / JVS(          337  )
  W( 59 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 338 )
  W( 88 ) = W( 88 ) + a*JVS( 339 )
  W( 91 ) = W( 91 ) + a*JVS( 340 )
  W( 96 ) = W( 96 ) + a*JVS( 341 )
  a = -W( 60 ) / JVS(          343  )
  W( 60 ) = -a
  W( 68 ) = W( 68 ) + a*JVS( 344 )
  W( 71 ) = W( 71 ) + a*JVS( 345 )
  W( 72 ) = W( 72 ) + a*JVS( 346 )
  W( 79 ) = W( 79 ) + a*JVS( 347 )
  W( 81 ) = W( 81 ) + a*JVS( 348 )
  W( 82 ) = W( 82 ) + a*JVS( 349 )
  W( 85 ) = W( 85 ) + a*JVS( 350 )
  W( 86 ) = W( 86 ) + a*JVS( 351 )
  W( 88 ) = W( 88 ) + a*JVS( 352 )
  W( 90 ) = W( 90 ) + a*JVS( 353 )
  W( 91 ) = W( 91 ) + a*JVS( 354 )
  W( 93 ) = W( 93 ) + a*JVS( 355 )
  W( 94 ) = W( 94 ) + a*JVS( 356 )
  W( 96 ) = W( 96 ) + a*JVS( 357 )
  a = -W( 61 ) / JVS(          358  )
  W( 61 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 359 )
  W( 89 ) = W( 89 ) + a*JVS( 360 )
  W( 90 ) = W( 90 ) + a*JVS( 361 )
  W( 91 ) = W( 91 ) + a*JVS( 362 )
  W( 93 ) = W( 93 ) + a*JVS( 363 )
  W( 94 ) = W( 94 ) + a*JVS( 364 )
  W( 96 ) = W( 96 ) + a*JVS( 365 )
  a = -W( 62 ) / JVS(          367  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 368 )
  W( 88 ) = W( 88 ) + a*JVS( 369 )
  W( 89 ) = W( 89 ) + a*JVS( 370 )
  W( 90 ) = W( 90 ) + a*JVS( 371 )
  W( 91 ) = W( 91 ) + a*JVS( 372 )
  W( 93 ) = W( 93 ) + a*JVS( 373 )
  W( 94 ) = W( 94 ) + a*JVS( 374 )
  W( 96 ) = W( 96 ) + a*JVS( 375 )
  a = -W( 63 ) / JVS(          377  )
  W( 63 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 378 )
  W( 89 ) = W( 89 ) + a*JVS( 379 )
  W( 90 ) = W( 90 ) + a*JVS( 380 )
  W( 91 ) = W( 91 ) + a*JVS( 381 )
  W( 93 ) = W( 93 ) + a*JVS( 382 )
  W( 94 ) = W( 94 ) + a*JVS( 383 )
  W( 96 ) = W( 96 ) + a*JVS( 384 )
  a = -W( 65 ) / JVS(          402  )
  W( 65 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 403 )
  W( 89 ) = W( 89 ) + a*JVS( 404 )
  W( 90 ) = W( 90 ) + a*JVS( 405 )
  W( 91 ) = W( 91 ) + a*JVS( 406 )
  W( 93 ) = W( 93 ) + a*JVS( 407 )
  W( 94 ) = W( 94 ) + a*JVS( 408 )
  W( 96 ) = W( 96 ) + a*JVS( 409 )
  a = -W( 68 ) / JVS(          451  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 452 )
  W( 88 ) = W( 88 ) + a*JVS( 453 )
  W( 89 ) = W( 89 ) + a*JVS( 454 )
  W( 90 ) = W( 90 ) + a*JVS( 455 )
  W( 91 ) = W( 91 ) + a*JVS( 456 )
  W( 93 ) = W( 93 ) + a*JVS( 457 )
  W( 94 ) = W( 94 ) + a*JVS( 458 )
  W( 95 ) = W( 95 ) + a*JVS( 459 )
  W( 96 ) = W( 96 ) + a*JVS( 460 )
  a = -W( 69 ) / JVS(          471  )
  W( 69 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 472 )
  W( 72 ) = W( 72 ) + a*JVS( 473 )
  W( 74 ) = W( 74 ) + a*JVS( 474 )
  W( 79 ) = W( 79 ) + a*JVS( 475 )
  W( 81 ) = W( 81 ) + a*JVS( 476 )
  W( 82 ) = W( 82 ) + a*JVS( 477 )
  W( 85 ) = W( 85 ) + a*JVS( 478 )
  W( 86 ) = W( 86 ) + a*JVS( 479 )
  W( 87 ) = W( 87 ) + a*JVS( 480 )
  W( 88 ) = W( 88 ) + a*JVS( 481 )
  W( 89 ) = W( 89 ) + a*JVS( 482 )
  W( 90 ) = W( 90 ) + a*JVS( 483 )
  W( 91 ) = W( 91 ) + a*JVS( 484 )
  W( 93 ) = W( 93 ) + a*JVS( 485 )
  W( 94 ) = W( 94 ) + a*JVS( 486 )
  W( 95 ) = W( 95 ) + a*JVS( 487 )
  W( 96 ) = W( 96 ) + a*JVS( 488 )
  a = -W( 70 ) / JVS(          491  )
  W( 70 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 492 )
  W( 85 ) = W( 85 ) + a*JVS( 493 )
  W( 89 ) = W( 89 ) + a*JVS( 494 )
  W( 90 ) = W( 90 ) + a*JVS( 495 )
  W( 91 ) = W( 91 ) + a*JVS( 496 )
  W( 93 ) = W( 93 ) + a*JVS( 497 )
  W( 94 ) = W( 94 ) + a*JVS( 498 )
  W( 96 ) = W( 96 ) + a*JVS( 499 )
  a = -W( 71 ) / JVS(          501  )
  W( 71 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 502 )
  W( 88 ) = W( 88 ) + a*JVS( 503 )
  W( 89 ) = W( 89 ) + a*JVS( 504 )
  W( 90 ) = W( 90 ) + a*JVS( 505 )
  W( 91 ) = W( 91 ) + a*JVS( 506 )
  W( 93 ) = W( 93 ) + a*JVS( 507 )
  W( 94 ) = W( 94 ) + a*JVS( 508 )
  W( 95 ) = W( 95 ) + a*JVS( 509 )
  W( 96 ) = W( 96 ) + a*JVS( 510 )
  a = -W( 72 ) / JVS(          512  )
  W( 72 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 513 )
  W( 88 ) = W( 88 ) + a*JVS( 514 )
  W( 89 ) = W( 89 ) + a*JVS( 515 )
  W( 90 ) = W( 90 ) + a*JVS( 516 )
  W( 91 ) = W( 91 ) + a*JVS( 517 )
  W( 93 ) = W( 93 ) + a*JVS( 518 )
  W( 94 ) = W( 94 ) + a*JVS( 519 )
  W( 95 ) = W( 95 ) + a*JVS( 520 )
  W( 96 ) = W( 96 ) + a*JVS( 521 )
  a = -W( 73 ) / JVS(          549  )
  W( 73 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 550 )
  W( 76 ) = W( 76 ) + a*JVS( 551 )
  W( 77 ) = W( 77 ) + a*JVS( 552 )
  W( 79 ) = W( 79 ) + a*JVS( 553 )
  W( 80 ) = W( 80 ) + a*JVS( 554 )
  W( 81 ) = W( 81 ) + a*JVS( 555 )
  W( 82 ) = W( 82 ) + a*JVS( 556 )
  W( 83 ) = W( 83 ) + a*JVS( 557 )
  W( 84 ) = W( 84 ) + a*JVS( 558 )
  W( 85 ) = W( 85 ) + a*JVS( 559 )
  W( 86 ) = W( 86 ) + a*JVS( 560 )
  W( 87 ) = W( 87 ) + a*JVS( 561 )
  W( 88 ) = W( 88 ) + a*JVS( 562 )
  W( 89 ) = W( 89 ) + a*JVS( 563 )
  W( 90 ) = W( 90 ) + a*JVS( 564 )
  W( 91 ) = W( 91 ) + a*JVS( 565 )
  W( 93 ) = W( 93 ) + a*JVS( 566 )
  W( 94 ) = W( 94 ) + a*JVS( 567 )
  W( 95 ) = W( 95 ) + a*JVS( 568 )
  W( 96 ) = W( 96 ) + a*JVS( 569 )
  W( 97 ) = W( 97 ) + a*JVS( 570 )
  a = -W( 74 ) / JVS(          572  )
  W( 74 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 573 )
  W( 85 ) = W( 85 ) + a*JVS( 574 )
  W( 89 ) = W( 89 ) + a*JVS( 575 )
  W( 90 ) = W( 90 ) + a*JVS( 576 )
  W( 91 ) = W( 91 ) + a*JVS( 577 )
  W( 93 ) = W( 93 ) + a*JVS( 578 )
  W( 94 ) = W( 94 ) + a*JVS( 579 )
  W( 95 ) = W( 95 ) + a*JVS( 580 )
  W( 96 ) = W( 96 ) + a*JVS( 581 )
  a = -W( 76 ) / JVS(          607  )
  W( 76 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 608 )
  W( 82 ) = W( 82 ) + a*JVS( 609 )
  W( 88 ) = W( 88 ) + a*JVS( 610 )
  W( 89 ) = W( 89 ) + a*JVS( 611 )
  W( 90 ) = W( 90 ) + a*JVS( 612 )
  W( 91 ) = W( 91 ) + a*JVS( 613 )
  W( 93 ) = W( 93 ) + a*JVS( 614 )
  W( 94 ) = W( 94 ) + a*JVS( 615 )
  W( 96 ) = W( 96 ) + a*JVS( 616 )
  a = -W( 77 ) / JVS(          618  )
  W( 77 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 619 )
  W( 82 ) = W( 82 ) + a*JVS( 620 )
  W( 88 ) = W( 88 ) + a*JVS( 621 )
  W( 89 ) = W( 89 ) + a*JVS( 622 )
  W( 90 ) = W( 90 ) + a*JVS( 623 )
  W( 91 ) = W( 91 ) + a*JVS( 624 )
  W( 93 ) = W( 93 ) + a*JVS( 625 )
  W( 94 ) = W( 94 ) + a*JVS( 626 )
  W( 96 ) = W( 96 ) + a*JVS( 627 )
  a = -W( 78 ) / JVS(          643  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 644 )
  W( 81 ) = W( 81 ) + a*JVS( 645 )
  W( 82 ) = W( 82 ) + a*JVS( 646 )
  W( 83 ) = W( 83 ) + a*JVS( 647 )
  W( 84 ) = W( 84 ) + a*JVS( 648 )
  W( 85 ) = W( 85 ) + a*JVS( 649 )
  W( 86 ) = W( 86 ) + a*JVS( 650 )
  W( 87 ) = W( 87 ) + a*JVS( 651 )
  W( 88 ) = W( 88 ) + a*JVS( 652 )
  W( 89 ) = W( 89 ) + a*JVS( 653 )
  W( 90 ) = W( 90 ) + a*JVS( 654 )
  W( 91 ) = W( 91 ) + a*JVS( 655 )
  W( 92 ) = W( 92 ) + a*JVS( 656 )
  W( 93 ) = W( 93 ) + a*JVS( 657 )
  W( 94 ) = W( 94 ) + a*JVS( 658 )
  W( 95 ) = W( 95 ) + a*JVS( 659 )
  W( 96 ) = W( 96 ) + a*JVS( 660 )
  W( 97 ) = W( 97 ) + a*JVS( 661 )
  a = -W( 79 ) / JVS(          664  )
  W( 79 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 665 )
  W( 88 ) = W( 88 ) + a*JVS( 666 )
  W( 89 ) = W( 89 ) + a*JVS( 667 )
  W( 90 ) = W( 90 ) + a*JVS( 668 )
  W( 91 ) = W( 91 ) + a*JVS( 669 )
  W( 93 ) = W( 93 ) + a*JVS( 670 )
  W( 94 ) = W( 94 ) + a*JVS( 671 )
  W( 95 ) = W( 95 ) + a*JVS( 672 )
  W( 96 ) = W( 96 ) + a*JVS( 673 )
  a = -W( 80 ) / JVS(          692  )
  W( 80 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 693 )
  W( 82 ) = W( 82 ) + a*JVS( 694 )
  W( 85 ) = W( 85 ) + a*JVS( 695 )
  W( 86 ) = W( 86 ) + a*JVS( 696 )
  W( 87 ) = W( 87 ) + a*JVS( 697 )
  W( 88 ) = W( 88 ) + a*JVS( 698 )
  W( 89 ) = W( 89 ) + a*JVS( 699 )
  W( 90 ) = W( 90 ) + a*JVS( 700 )
  W( 91 ) = W( 91 ) + a*JVS( 701 )
  W( 92 ) = W( 92 ) + a*JVS( 702 )
  W( 93 ) = W( 93 ) + a*JVS( 703 )
  W( 94 ) = W( 94 ) + a*JVS( 704 )
  W( 95 ) = W( 95 ) + a*JVS( 705 )
  W( 96 ) = W( 96 ) + a*JVS( 706 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  a = -W( 83 ) / JVS(          739  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 740 )
  W( 85 ) = W( 85 ) + a*JVS( 741 )
  W( 88 ) = W( 88 ) + a*JVS( 742 )
  W( 89 ) = W( 89 ) + a*JVS( 743 )
  W( 90 ) = W( 90 ) + a*JVS( 744 )
  W( 91 ) = W( 91 ) + a*JVS( 745 )
  W( 93 ) = W( 93 ) + a*JVS( 746 )
  W( 94 ) = W( 94 ) + a*JVS( 747 )
  W( 96 ) = W( 96 ) + a*JVS( 748 )
  a = -W( 84 ) / JVS(          758  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 759 )
  W( 88 ) = W( 88 ) + a*JVS( 760 )
  W( 89 ) = W( 89 ) + a*JVS( 761 )
  W( 90 ) = W( 90 ) + a*JVS( 762 )
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 93 ) = W( 93 ) + a*JVS( 764 )
  W( 94 ) = W( 94 ) + a*JVS( 765 )
  W( 96 ) = W( 96 ) + a*JVS( 766 )
  a = -W( 85 ) / JVS(          786  )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 787 )
  W( 89 ) = W( 89 ) + a*JVS( 788 )
  W( 90 ) = W( 90 ) + a*JVS( 789 )
  W( 91 ) = W( 91 ) + a*JVS( 790 )
  W( 93 ) = W( 93 ) + a*JVS( 791 )
  W( 94 ) = W( 94 ) + a*JVS( 792 )
  W( 95 ) = W( 95 ) + a*JVS( 793 )
  W( 96 ) = W( 96 ) + a*JVS( 794 )
  a = -W( 86 ) / JVS(          805  )
  W( 86 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 806 )
  W( 89 ) = W( 89 ) + a*JVS( 807 )
  W( 90 ) = W( 90 ) + a*JVS( 808 )
  W( 91 ) = W( 91 ) + a*JVS( 809 )
  W( 92 ) = W( 92 ) + a*JVS( 810 )
  W( 93 ) = W( 93 ) + a*JVS( 811 )
  W( 94 ) = W( 94 ) + a*JVS( 812 )
  W( 95 ) = W( 95 ) + a*JVS( 813 )
  W( 96 ) = W( 96 ) + a*JVS( 814 )
  a = -W( 87 ) / JVS(          831  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 832 )
  W( 89 ) = W( 89 ) + a*JVS( 833 )
  W( 90 ) = W( 90 ) + a*JVS( 834 )
  W( 91 ) = W( 91 ) + a*JVS( 835 )
  W( 92 ) = W( 92 ) + a*JVS( 836 )
  W( 93 ) = W( 93 ) + a*JVS( 837 )
  W( 94 ) = W( 94 ) + a*JVS( 838 )
  W( 95 ) = W( 95 ) + a*JVS( 839 )
  W( 96 ) = W( 96 ) + a*JVS( 840 )
  a = -W( 88 ) / JVS(          871  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 872 )
  W( 90 ) = W( 90 ) + a*JVS( 873 )
  W( 91 ) = W( 91 ) + a*JVS( 874 )
  W( 92 ) = W( 92 ) + a*JVS( 875 )
  W( 93 ) = W( 93 ) + a*JVS( 876 )
  W( 94 ) = W( 94 ) + a*JVS( 877 )
  W( 95 ) = W( 95 ) + a*JVS( 878 )
  W( 96 ) = W( 96 ) + a*JVS( 879 )
  W( 97 ) = W( 97 ) + a*JVS( 880 )
  a = -W( 89 ) / JVS(          940  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 941 )
  W( 91 ) = W( 91 ) + a*JVS( 942 )
  W( 92 ) = W( 92 ) + a*JVS( 943 )
  W( 93 ) = W( 93 ) + a*JVS( 944 )
  W( 94 ) = W( 94 ) + a*JVS( 945 )
  W( 95 ) = W( 95 ) + a*JVS( 946 )
  W( 96 ) = W( 96 ) + a*JVS( 947 )
  W( 97 ) = W( 97 ) + a*JVS( 948 )
  a = -W( 90 ) / JVS(          979  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 980 )
  W( 92 ) = W( 92 ) + a*JVS( 981 )
  W( 93 ) = W( 93 ) + a*JVS( 982 )
  W( 94 ) = W( 94 ) + a*JVS( 983 )
  W( 95 ) = W( 95 ) + a*JVS( 984 )
  W( 96 ) = W( 96 ) + a*JVS( 985 )
  W( 97 ) = W( 97 ) + a*JVS( 986 )
  a = -W( 91 ) / JVS(         1055  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 1056 )
  W( 93 ) = W( 93 ) + a*JVS( 1057 )
  W( 94 ) = W( 94 ) + a*JVS( 1058 )
  W( 95 ) = W( 95 ) + a*JVS( 1059 )
  W( 96 ) = W( 96 ) + a*JVS( 1060 )
  W( 97 ) = W( 97 ) + a*JVS( 1061 )
  a = -W( 92 ) / JVS(         1084  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 1085 )
  W( 94 ) = W( 94 ) + a*JVS( 1086 )
  W( 95 ) = W( 95 ) + a*JVS( 1087 )
  W( 96 ) = W( 96 ) + a*JVS( 1088 )
  W( 97 ) = W( 97 ) + a*JVS( 1089 )
  a = -W( 93 ) / JVS(         1129  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 1130 )
  W( 95 ) = W( 95 ) + a*JVS( 1131 )
  W( 96 ) = W( 96 ) + a*JVS( 1132 )
  W( 97 ) = W( 97 ) + a*JVS( 1133 )
  a = -W( 94 ) / JVS(         1172  )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 1173 )
  W( 96 ) = W( 96 ) + a*JVS( 1174 )
  W( 97 ) = W( 97 ) + a*JVS( 1175 )
  a = -W( 95 ) / JVS(         1230  )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 1231 )
  W( 97 ) = W( 97 ) + a*JVS( 1232 )
  JVS( 1233) = W( 17 )
  JVS( 1234) = W( 27 )
  JVS( 1235) = W( 30 )
  JVS( 1236) = W( 43 )
  JVS( 1237) = W( 44 )
  JVS( 1238) = W( 47 )
  JVS( 1239) = W( 48 )
  JVS( 1240) = W( 51 )
  JVS( 1241) = W( 53 )
  JVS( 1242) = W( 54 )
  JVS( 1243) = W( 55 )
  JVS( 1244) = W( 57 )
  JVS( 1245) = W( 58 )
  JVS( 1246) = W( 59 )
  JVS( 1247) = W( 60 )
  JVS( 1248) = W( 61 )
  JVS( 1249) = W( 62 )
  JVS( 1250) = W( 63 )
  JVS( 1251) = W( 65 )
  JVS( 1252) = W( 68 )
  JVS( 1253) = W( 69 )
  JVS( 1254) = W( 70 )
  JVS( 1255) = W( 71 )
  JVS( 1256) = W( 72 )
  JVS( 1257) = W( 73 )
  JVS( 1258) = W( 74 )
  JVS( 1259) = W( 76 )
  JVS( 1260) = W( 77 )
  JVS( 1261) = W( 78 )
  JVS( 1262) = W( 79 )
  JVS( 1263) = W( 80 )
  JVS( 1264) = W( 81 )
  JVS( 1265) = W( 82 )
  JVS( 1266) = W( 83 )
  JVS( 1267) = W( 84 )
  JVS( 1268) = W( 85 )
  JVS( 1269) = W( 86 )
  JVS( 1270) = W( 87 )
  JVS( 1271) = W( 88 )
  JVS( 1272) = W( 89 )
  JVS( 1273) = W( 90 )
  JVS( 1274) = W( 91 )
  JVS( 1275) = W( 92 )
  JVS( 1276) = W( 93 )
  JVS( 1277) = W( 94 )
  JVS( 1278) = W( 95 )
  JVS( 1279) = W( 96 )
  JVS( 1280) = W( 97 )
  IF ( ABS(  JVS( 1307 )) < TINY(a) ) THEN
         IER = 97                                      
         RETURN
  END IF
   W( 47 ) = JVS( 1281 )
   W( 50 ) = JVS( 1282 )
   W( 51 ) = JVS( 1283 )
   W( 54 ) = JVS( 1284 )
   W( 66 ) = JVS( 1285 )
   W( 70 ) = JVS( 1286 )
   W( 75 ) = JVS( 1287 )
   W( 76 ) = JVS( 1288 )
   W( 77 ) = JVS( 1289 )
   W( 79 ) = JVS( 1290 )
   W( 81 ) = JVS( 1291 )
   W( 82 ) = JVS( 1292 )
   W( 83 ) = JVS( 1293 )
   W( 84 ) = JVS( 1294 )
   W( 85 ) = JVS( 1295 )
   W( 86 ) = JVS( 1296 )
   W( 87 ) = JVS( 1297 )
   W( 88 ) = JVS( 1298 )
   W( 89 ) = JVS( 1299 )
   W( 90 ) = JVS( 1300 )
   W( 91 ) = JVS( 1301 )
   W( 92 ) = JVS( 1302 )
   W( 93 ) = JVS( 1303 )
   W( 94 ) = JVS( 1304 )
   W( 95 ) = JVS( 1305 )
   W( 96 ) = JVS( 1306 )
   W( 97 ) = JVS( 1307 )
  a = -W( 47 ) / JVS(          243  )
  W( 47 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 244 )
  W( 91 ) = W( 91 ) + a*JVS( 245 )
  W( 96 ) = W( 96 ) + a*JVS( 246 )
  a = -W( 50 ) / JVS(          259  )
  W( 50 ) = -a
  W( 54 ) = W( 54 ) + a*JVS( 260 )
  W( 81 ) = W( 81 ) + a*JVS( 261 )
  W( 82 ) = W( 82 ) + a*JVS( 262 )
  W( 85 ) = W( 85 ) + a*JVS( 263 )
  W( 88 ) = W( 88 ) + a*JVS( 264 )
  W( 91 ) = W( 91 ) + a*JVS( 265 )
  a = -W( 51 ) / JVS(          266  )
  W( 51 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 267 )
  W( 91 ) = W( 91 ) + a*JVS( 268 )
  W( 96 ) = W( 96 ) + a*JVS( 269 )
  a = -W( 54 ) / JVS(          302  )
  W( 54 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 303 )
  W( 91 ) = W( 91 ) + a*JVS( 304 )
  W( 96 ) = W( 96 ) + a*JVS( 305 )
  a = -W( 66 ) / JVS(          417  )
  W( 66 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 418 )
  W( 75 ) = W( 75 ) + a*JVS( 419 )
  W( 76 ) = W( 76 ) + a*JVS( 420 )
  W( 77 ) = W( 77 ) + a*JVS( 421 )
  W( 79 ) = W( 79 ) + a*JVS( 422 )
  W( 81 ) = W( 81 ) + a*JVS( 423 )
  W( 82 ) = W( 82 ) + a*JVS( 424 )
  W( 83 ) = W( 83 ) + a*JVS( 425 )
  W( 84 ) = W( 84 ) + a*JVS( 426 )
  W( 85 ) = W( 85 ) + a*JVS( 427 )
  W( 86 ) = W( 86 ) + a*JVS( 428 )
  W( 87 ) = W( 87 ) + a*JVS( 429 )
  W( 88 ) = W( 88 ) + a*JVS( 430 )
  W( 89 ) = W( 89 ) + a*JVS( 431 )
  W( 90 ) = W( 90 ) + a*JVS( 432 )
  W( 91 ) = W( 91 ) + a*JVS( 433 )
  W( 92 ) = W( 92 ) + a*JVS( 434 )
  W( 93 ) = W( 93 ) + a*JVS( 435 )
  W( 94 ) = W( 94 ) + a*JVS( 436 )
  W( 96 ) = W( 96 ) + a*JVS( 437 )
  a = -W( 70 ) / JVS(          491  )
  W( 70 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 492 )
  W( 85 ) = W( 85 ) + a*JVS( 493 )
  W( 89 ) = W( 89 ) + a*JVS( 494 )
  W( 90 ) = W( 90 ) + a*JVS( 495 )
  W( 91 ) = W( 91 ) + a*JVS( 496 )
  W( 93 ) = W( 93 ) + a*JVS( 497 )
  W( 94 ) = W( 94 ) + a*JVS( 498 )
  W( 96 ) = W( 96 ) + a*JVS( 499 )
  a = -W( 75 ) / JVS(          590  )
  W( 75 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 591 )
  W( 77 ) = W( 77 ) + a*JVS( 592 )
  W( 79 ) = W( 79 ) + a*JVS( 593 )
  W( 83 ) = W( 83 ) + a*JVS( 594 )
  W( 84 ) = W( 84 ) + a*JVS( 595 )
  W( 85 ) = W( 85 ) + a*JVS( 596 )
  W( 86 ) = W( 86 ) + a*JVS( 597 )
  W( 88 ) = W( 88 ) + a*JVS( 598 )
  W( 89 ) = W( 89 ) + a*JVS( 599 )
  W( 90 ) = W( 90 ) + a*JVS( 600 )
  W( 91 ) = W( 91 ) + a*JVS( 601 )
  W( 93 ) = W( 93 ) + a*JVS( 602 )
  W( 94 ) = W( 94 ) + a*JVS( 603 )
  W( 95 ) = W( 95 ) + a*JVS( 604 )
  W( 96 ) = W( 96 ) + a*JVS( 605 )
  a = -W( 76 ) / JVS(          607  )
  W( 76 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 608 )
  W( 82 ) = W( 82 ) + a*JVS( 609 )
  W( 88 ) = W( 88 ) + a*JVS( 610 )
  W( 89 ) = W( 89 ) + a*JVS( 611 )
  W( 90 ) = W( 90 ) + a*JVS( 612 )
  W( 91 ) = W( 91 ) + a*JVS( 613 )
  W( 93 ) = W( 93 ) + a*JVS( 614 )
  W( 94 ) = W( 94 ) + a*JVS( 615 )
  W( 96 ) = W( 96 ) + a*JVS( 616 )
  a = -W( 77 ) / JVS(          618  )
  W( 77 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 619 )
  W( 82 ) = W( 82 ) + a*JVS( 620 )
  W( 88 ) = W( 88 ) + a*JVS( 621 )
  W( 89 ) = W( 89 ) + a*JVS( 622 )
  W( 90 ) = W( 90 ) + a*JVS( 623 )
  W( 91 ) = W( 91 ) + a*JVS( 624 )
  W( 93 ) = W( 93 ) + a*JVS( 625 )
  W( 94 ) = W( 94 ) + a*JVS( 626 )
  W( 96 ) = W( 96 ) + a*JVS( 627 )
  a = -W( 79 ) / JVS(          664  )
  W( 79 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 665 )
  W( 88 ) = W( 88 ) + a*JVS( 666 )
  W( 89 ) = W( 89 ) + a*JVS( 667 )
  W( 90 ) = W( 90 ) + a*JVS( 668 )
  W( 91 ) = W( 91 ) + a*JVS( 669 )
  W( 93 ) = W( 93 ) + a*JVS( 670 )
  W( 94 ) = W( 94 ) + a*JVS( 671 )
  W( 95 ) = W( 95 ) + a*JVS( 672 )
  W( 96 ) = W( 96 ) + a*JVS( 673 )
  a = -W( 81 ) / JVS(          709  )
  W( 81 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 710 )
  W( 88 ) = W( 88 ) + a*JVS( 711 )
  W( 89 ) = W( 89 ) + a*JVS( 712 )
  W( 90 ) = W( 90 ) + a*JVS( 713 )
  W( 91 ) = W( 91 ) + a*JVS( 714 )
  W( 93 ) = W( 93 ) + a*JVS( 715 )
  W( 94 ) = W( 94 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  a = -W( 82 ) / JVS(          722  )
  W( 82 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 723 )
  W( 88 ) = W( 88 ) + a*JVS( 724 )
  W( 89 ) = W( 89 ) + a*JVS( 725 )
  W( 90 ) = W( 90 ) + a*JVS( 726 )
  W( 91 ) = W( 91 ) + a*JVS( 727 )
  W( 93 ) = W( 93 ) + a*JVS( 728 )
  W( 94 ) = W( 94 ) + a*JVS( 729 )
  W( 96 ) = W( 96 ) + a*JVS( 730 )
  a = -W( 83 ) / JVS(          739  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 740 )
  W( 85 ) = W( 85 ) + a*JVS( 741 )
  W( 88 ) = W( 88 ) + a*JVS( 742 )
  W( 89 ) = W( 89 ) + a*JVS( 743 )
  W( 90 ) = W( 90 ) + a*JVS( 744 )
  W( 91 ) = W( 91 ) + a*JVS( 745 )
  W( 93 ) = W( 93 ) + a*JVS( 746 )
  W( 94 ) = W( 94 ) + a*JVS( 747 )
  W( 96 ) = W( 96 ) + a*JVS( 748 )
  a = -W( 84 ) / JVS(          758  )
  W( 84 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 759 )
  W( 88 ) = W( 88 ) + a*JVS( 760 )
  W( 89 ) = W( 89 ) + a*JVS( 761 )
  W( 90 ) = W( 90 ) + a*JVS( 762 )
  W( 91 ) = W( 91 ) + a*JVS( 763 )
  W( 93 ) = W( 93 ) + a*JVS( 764 )
  W( 94 ) = W( 94 ) + a*JVS( 765 )
  W( 96 ) = W( 96 ) + a*JVS( 766 )
  a = -W( 85 ) / JVS(          786  )
  W( 85 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 787 )
  W( 89 ) = W( 89 ) + a*JVS( 788 )
  W( 90 ) = W( 90 ) + a*JVS( 789 )
  W( 91 ) = W( 91 ) + a*JVS( 790 )
  W( 93 ) = W( 93 ) + a*JVS( 791 )
  W( 94 ) = W( 94 ) + a*JVS( 792 )
  W( 95 ) = W( 95 ) + a*JVS( 793 )
  W( 96 ) = W( 96 ) + a*JVS( 794 )
  a = -W( 86 ) / JVS(          805  )
  W( 86 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 806 )
  W( 89 ) = W( 89 ) + a*JVS( 807 )
  W( 90 ) = W( 90 ) + a*JVS( 808 )
  W( 91 ) = W( 91 ) + a*JVS( 809 )
  W( 92 ) = W( 92 ) + a*JVS( 810 )
  W( 93 ) = W( 93 ) + a*JVS( 811 )
  W( 94 ) = W( 94 ) + a*JVS( 812 )
  W( 95 ) = W( 95 ) + a*JVS( 813 )
  W( 96 ) = W( 96 ) + a*JVS( 814 )
  a = -W( 87 ) / JVS(          831  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 832 )
  W( 89 ) = W( 89 ) + a*JVS( 833 )
  W( 90 ) = W( 90 ) + a*JVS( 834 )
  W( 91 ) = W( 91 ) + a*JVS( 835 )
  W( 92 ) = W( 92 ) + a*JVS( 836 )
  W( 93 ) = W( 93 ) + a*JVS( 837 )
  W( 94 ) = W( 94 ) + a*JVS( 838 )
  W( 95 ) = W( 95 ) + a*JVS( 839 )
  W( 96 ) = W( 96 ) + a*JVS( 840 )
  a = -W( 88 ) / JVS(          871  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 872 )
  W( 90 ) = W( 90 ) + a*JVS( 873 )
  W( 91 ) = W( 91 ) + a*JVS( 874 )
  W( 92 ) = W( 92 ) + a*JVS( 875 )
  W( 93 ) = W( 93 ) + a*JVS( 876 )
  W( 94 ) = W( 94 ) + a*JVS( 877 )
  W( 95 ) = W( 95 ) + a*JVS( 878 )
  W( 96 ) = W( 96 ) + a*JVS( 879 )
  W( 97 ) = W( 97 ) + a*JVS( 880 )
  a = -W( 89 ) / JVS(          940  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 941 )
  W( 91 ) = W( 91 ) + a*JVS( 942 )
  W( 92 ) = W( 92 ) + a*JVS( 943 )
  W( 93 ) = W( 93 ) + a*JVS( 944 )
  W( 94 ) = W( 94 ) + a*JVS( 945 )
  W( 95 ) = W( 95 ) + a*JVS( 946 )
  W( 96 ) = W( 96 ) + a*JVS( 947 )
  W( 97 ) = W( 97 ) + a*JVS( 948 )
  a = -W( 90 ) / JVS(          979  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 980 )
  W( 92 ) = W( 92 ) + a*JVS( 981 )
  W( 93 ) = W( 93 ) + a*JVS( 982 )
  W( 94 ) = W( 94 ) + a*JVS( 983 )
  W( 95 ) = W( 95 ) + a*JVS( 984 )
  W( 96 ) = W( 96 ) + a*JVS( 985 )
  W( 97 ) = W( 97 ) + a*JVS( 986 )
  a = -W( 91 ) / JVS(         1055  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 1056 )
  W( 93 ) = W( 93 ) + a*JVS( 1057 )
  W( 94 ) = W( 94 ) + a*JVS( 1058 )
  W( 95 ) = W( 95 ) + a*JVS( 1059 )
  W( 96 ) = W( 96 ) + a*JVS( 1060 )
  W( 97 ) = W( 97 ) + a*JVS( 1061 )
  a = -W( 92 ) / JVS(         1084  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 1085 )
  W( 94 ) = W( 94 ) + a*JVS( 1086 )
  W( 95 ) = W( 95 ) + a*JVS( 1087 )
  W( 96 ) = W( 96 ) + a*JVS( 1088 )
  W( 97 ) = W( 97 ) + a*JVS( 1089 )
  a = -W( 93 ) / JVS(         1129  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 1130 )
  W( 95 ) = W( 95 ) + a*JVS( 1131 )
  W( 96 ) = W( 96 ) + a*JVS( 1132 )
  W( 97 ) = W( 97 ) + a*JVS( 1133 )
  a = -W( 94 ) / JVS(         1172  )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 1173 )
  W( 96 ) = W( 96 ) + a*JVS( 1174 )
  W( 97 ) = W( 97 ) + a*JVS( 1175 )
  a = -W( 95 ) / JVS(         1230  )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 1231 )
  W( 97 ) = W( 97 ) + a*JVS( 1232 )
  a = -W( 96 ) / JVS(         1279  )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 1280 )
  JVS( 1281) = W( 47 )
  JVS( 1282) = W( 50 )
  JVS( 1283) = W( 51 )
  JVS( 1284) = W( 54 )
  JVS( 1285) = W( 66 )
  JVS( 1286) = W( 70 )
  JVS( 1287) = W( 75 )
  JVS( 1288) = W( 76 )
  JVS( 1289) = W( 77 )
  JVS( 1290) = W( 79 )
  JVS( 1291) = W( 81 )
  JVS( 1292) = W( 82 )
  JVS( 1293) = W( 83 )
  JVS( 1294) = W( 84 )
  JVS( 1295) = W( 85 )
  JVS( 1296) = W( 86 )
  JVS( 1297) = W( 87 )
  JVS( 1298) = W( 88 )
  JVS( 1299) = W( 89 )
  JVS( 1300) = W( 90 )
  JVS( 1301) = W( 91 )
  JVS( 1302) = W( 92 )
  JVS( 1303) = W( 93 )
  JVS( 1304) = W( 94 )
  JVS( 1305) = W( 95 )
  JVS( 1306) = W( 96 )
  JVS( 1307) = W( 97 )
   
   END SUBROUTINE decomp_racm_soa_vbs_het
 


END MODULE racm_soa_vbs_het_Integrator
