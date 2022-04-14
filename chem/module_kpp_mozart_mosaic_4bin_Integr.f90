
































MODULE mozart_mosaic_4bin_Integrator

 USE mozart_mosaic_4bin_Parameters
 USE mozart_mosaic_4bin_Precision
 USE mozart_mosaic_4bin_JacobianSP

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

SUBROUTINE  mozart_mosaic_4bin_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE mozart_mosaic_4bin_Parameters

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

   CALL mozart_mosaic_4bin_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  mozart_mosaic_4bin_INTEGRATE


SUBROUTINE  mozart_mosaic_4bin_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE mozart_mosaic_4bin_Parameters

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
      CALL mozart_mosaic_4bin_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL mozart_mosaic_4bin_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = mozart_mosaic_4bin_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL mozart_mosaic_4bin_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL mozart_mosaic_4bin_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL mozart_mosaic_4bin_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL mozart_mosaic_4bin_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL mozart_mosaic_4bin_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL mozart_mosaic_4bin_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL mozart_mosaic_4bin_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL mozart_mosaic_4bin_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL mozart_mosaic_4bin_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL mozart_mosaic_4bin_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL mozart_mosaic_4bin_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL mozart_mosaic_4bin_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL mozart_mosaic_4bin_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL mozart_mosaic_4bin_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL mozart_mosaic_4bin_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  mozart_mosaic_4bin_ros_ErrorMsg(Code,T,H,IERR)



   USE mozart_mosaic_4bin_Precision

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

 END SUBROUTINE  mozart_mosaic_4bin_ros_ErrorMsg


 SUBROUTINE  mozart_mosaic_4bin_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL mozart_mosaic_4bin_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL mozart_mosaic_4bin_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL mozart_mosaic_4bin_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL mozart_mosaic_4bin_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL mozart_mosaic_4bin_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL mozart_mosaic_4bin_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL mozart_mosaic_4bin_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL mozart_mosaic_4bin_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL mozart_mosaic_4bin_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL mozart_mosaic_4bin_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL mozart_mosaic_4bin_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL mozart_mosaic_4bin_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL mozart_mosaic_4bin_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL mozart_mosaic_4bin_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL mozart_mosaic_4bin_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL mozart_mosaic_4bin_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL mozart_mosaic_4bin_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL mozart_mosaic_4bin_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL mozart_mosaic_4bin_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL mozart_mosaic_4bin_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = mozart_mosaic_4bin_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL mozart_mosaic_4bin_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  mozart_mosaic_4bin_ros_Integrator



  REAL(kind=dp) FUNCTION  mozart_mosaic_4bin_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    mozart_mosaic_4bin_ros_ErrorNorm = Err

  END FUNCTION  mozart_mosaic_4bin_ros_ErrorNorm



  SUBROUTINE mozart_mosaic_4bin_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL mozart_mosaic_4bin_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL mozart_mosaic_4bin_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL mozart_mosaic_4bin_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  mozart_mosaic_4bin_ros_FunTimeDeriv



  SUBROUTINE  mozart_mosaic_4bin_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL mozart_mosaic_4bin_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL mozart_mosaic_4bin_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL mozart_mosaic_4bin_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  mozart_mosaic_4bin_ros_PrepareMatrix



  SUBROUTINE  mozart_mosaic_4bin_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 



CALL decomp_mozart_mosaic_4bin ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  mozart_mosaic_4bin_ros_Decomp



  SUBROUTINE  mozart_mosaic_4bin_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL mozart_mosaic_4bin_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  mozart_mosaic_4bin_ros_Solve




  SUBROUTINE  mozart_mosaic_4bin_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  mozart_mosaic_4bin_Ros2



  SUBROUTINE  mozart_mosaic_4bin_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  mozart_mosaic_4bin_Ros3





  SUBROUTINE  mozart_mosaic_4bin_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  mozart_mosaic_4bin_Ros4


  SUBROUTINE  mozart_mosaic_4bin_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  mozart_mosaic_4bin_Rodas3


  SUBROUTINE  mozart_mosaic_4bin_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  mozart_mosaic_4bin_Rodas4




END SUBROUTINE  mozart_mosaic_4bin_Rosenbrock




SUBROUTINE  mozart_mosaic_4bin_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE mozart_mosaic_4bin_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL mozart_mosaic_4bin_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  mozart_mosaic_4bin_FunTemplate



SUBROUTINE  mozart_mosaic_4bin_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE mozart_mosaic_4bin_Parameters
 
 USE mozart_mosaic_4bin_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL mozart_mosaic_4bin_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  mozart_mosaic_4bin_JacTemplate

















SUBROUTINE mozart_mosaic_4bin_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*F(2)
  A(2) = RCT(2)*V(132)
  A(3) = RCT(3)*V(132)
  A(4) = RCT(4)*V(15)
  A(5) = RCT(5)*V(134)
  A(6) = RCT(6)*V(36)
  A(7) = RCT(7)*V(71)
  A(8) = RCT(8)*V(135)
  A(9) = RCT(9)*V(52)
  A(10) = RCT(10)*V(50)
  A(11) = RCT(11)*V(126)
  A(12) = RCT(12)*V(126)
  A(13) = RCT(13)*V(18)
  A(14) = RCT(14)*V(112)
  A(15) = RCT(15)*V(70)
  A(16) = RCT(16)*V(64)
  A(17) = RCT(17)*V(68)
  A(18) = RCT(18)*V(77)
  A(19) = RCT(19)*V(115)
  A(20) = RCT(20)*V(124)
  A(21) = RCT(21)*V(45)
  A(22) = RCT(22)*V(46)
  A(23) = RCT(23)*V(47)
  A(24) = RCT(24)*V(105)
  A(25) = RCT(25)*V(117)
  A(26) = RCT(26)*V(26)
  A(27) = RCT(27)*V(116)
  A(28) = RCT(28)*V(75)
  A(29) = RCT(29)*V(122)
  A(30) = RCT(30)*V(108)
  A(31) = RCT(31)*V(62)
  A(32) = RCT(32)*V(98)
  A(33) = RCT(33)*V(96)
  A(34) = RCT(34)*V(84)
  A(35) = RCT(35)*V(49)
  A(36) = RCT(36)*V(57)
  A(37) = RCT(37)*V(43)
  A(38) = RCT(38)*V(21)
  A(39) = RCT(39)*V(31)
  A(40) = RCT(40)*V(30)
  A(41) = RCT(41)*V(42)
  A(42) = RCT(42)*V(40)
  A(43) = RCT(43)*V(76)
  A(44) = RCT(44)*V(103)
  A(45) = RCT(45)*V(110)
  A(46) = RCT(46)*V(114)
  A(47) = RCT(47)*V(65)
  A(48) = RCT(48)*V(79)
  A(49) = RCT(49)*V(1)
  A(50) = RCT(50)*V(92)*F(2)
  A(51) = RCT(51)*V(92)*V(132)
  A(52) = RCT(52)*V(80)*F(2)
  A(53) = RCT(53)*V(80)*F(1)
  A(54) = RCT(54)*V(44)*V(80)
  A(55) = RCT(55)*V(44)*V(136)
  A(56) = RCT(56)*V(92)*V(136)
  A(57) = RCT(57)*V(92)*V(131)
  A(58) = RCT(58)*V(132)*V(136)
  A(59) = RCT(59)*V(131)*V(132)
  A(60) = RCT(60)*V(131)*V(131)*F(1)
  A(61) = RCT(61)*V(18)*V(136)
  A(62) = RCT(62)*V(131)*V(136)
  A(63) = RCT(63)*V(136)*V(136)
  A(64) = RCT(64)*V(136)*V(136)
  A(65) = RCT(65)*V(15)*V(80)
  A(66) = RCT(66)*V(15)*V(80)
  A(67) = RCT(67)*V(131)*V(133)
  A(68) = RCT(68)*V(132)*V(133)
  A(69) = RCT(69)*V(92)*V(134)
  A(70) = RCT(70)*V(132)*V(134)
  A(71) = RCT(71)*V(131)*V(135)
  A(72) = RCT(72)*V(134)*V(135)
  A(73) = RCT(73)*V(36)
  A(74) = RCT(74)*V(134)*V(136)
  A(75) = RCT(75)*V(71)*V(136)
  A(76) = RCT(76)*V(133)*V(135)
  A(77) = RCT(77)*V(131)*V(134)
  A(78) = RCT(78)*V(52)*V(136)
  A(79) = RCT(79)*V(52)
  A(80) = RCT(80)*V(36)*F(2)
  A(81) = RCT(81)*V(135)
  A(82) = RCT(82)*V(134)
  A(83) = RCT(83)*V(72)*V(136)
  A(84) = RCT(84)*V(72)*V(80)
  A(85) = RCT(85)*V(129)*V(133)
  A(86) = RCT(86)*V(129)*V(129)
  A(87) = RCT(87)*V(129)*V(129)
  A(88) = RCT(88)*V(129)*V(131)
  A(89) = RCT(89)*V(50)*V(136)
  A(90) = RCT(90)*V(126)*V(135)
  A(91) = RCT(91)*V(126)*V(136)
  A(92) = RCT(92)*V(107)*V(136)
  A(93) = RCT(93)*V(55)*V(136)
  A(94) = RCT(94)*V(55)*V(132)
  A(95) = RCT(95)*V(16)*V(136)
  A(96) = RCT(96)*V(96)*V(136)
  A(97) = RCT(97)*V(60)*V(133)
  A(98) = RCT(98)*V(23)*F(2)
  A(99) = RCT(99)*V(23)
  A(100) = RCT(100)*V(12)*V(136)
  A(101) = RCT(101)*V(97)*V(133)
  A(102) = RCT(102)*V(97)*V(131)
  A(103) = RCT(103)*V(97)*V(129)
  A(104) = RCT(104)*V(45)*V(136)
  A(105) = RCT(105)*V(109)*V(136)
  A(106) = RCT(106)*V(109)*V(132)
  A(107) = RCT(107)*V(109)*V(135)
  A(108) = RCT(108)*V(86)*V(133)
  A(109) = RCT(109)*V(86)*V(131)
  A(110) = RCT(110)*V(70)*V(136)
  A(111) = RCT(111)*V(112)*V(136)
  A(112) = RCT(112)*V(112)*V(135)
  A(113) = RCT(113)*V(128)*V(133)
  A(114) = RCT(114)*V(128)*V(134)
  A(115) = RCT(115)*V(128)*V(131)
  A(116) = RCT(116)*V(128)*V(129)
  A(117) = RCT(117)*V(64)*V(136)
  A(118) = RCT(118)*V(68)
  A(119) = RCT(119)*V(128)*V(128)
  A(120) = RCT(120)*V(13)*V(136)
  A(121) = RCT(121)*V(93)*V(133)
  A(122) = RCT(122)*V(93)*V(131)
  A(123) = RCT(123)*V(93)*V(129)
  A(124) = RCT(124)*V(46)*V(136)
  A(125) = RCT(125)*V(105)*V(136)
  A(126) = RCT(126)*V(121)*V(133)
  A(127) = RCT(127)*V(121)*V(131)
  A(128) = RCT(128)*V(121)*V(129)
  A(129) = RCT(129)*V(47)*V(136)
  A(130) = RCT(130)*V(14)*V(136)
  A(131) = RCT(131)*V(48)*V(133)
  A(132) = RCT(132)*V(54)*V(136)
  A(133) = RCT(133)*V(113)*V(133)
  A(134) = RCT(134)*V(113)*V(131)
  A(135) = RCT(135)*V(84)*V(136)
  A(136) = RCT(136)*V(67)*V(136)
  A(137) = RCT(137)*V(62)*V(136)
  A(138) = RCT(138)*V(85)*V(133)
  A(139) = RCT(139)*V(85)*V(131)
  A(140) = RCT(140)*V(49)*V(136)
  A(141) = RCT(141)*V(88)*V(131)
  A(142) = RCT(142)*V(57)*V(136)
  A(143) = RCT(143)*V(106)*V(136)
  A(144) = RCT(144)*V(106)*V(132)
  A(145) = RCT(145)*V(127)*V(133)
  A(146) = RCT(146)*V(127)*V(135)
  A(147) = RCT(147)*V(127)*V(131)
  A(148) = RCT(148)*V(75)*V(136)
  A(149) = RCT(149)*V(127)*V(129)
  A(150) = RCT(150)*V(127)*V(128)
  A(151) = RCT(151)*V(124)*V(136)
  A(152) = RCT(152)*V(124)*V(132)
  A(153) = RCT(153)*V(115)*V(136)
  A(154) = RCT(154)*V(115)*V(132)
  A(155) = RCT(155)*V(125)*V(133)
  A(156) = RCT(156)*V(125)*V(133)
  A(157) = RCT(157)*V(125)*V(135)
  A(158) = RCT(158)*V(125)*V(131)
  A(159) = RCT(159)*V(125)*V(129)
  A(160) = RCT(160)*V(125)*V(128)
  A(161) = RCT(161)*V(51)*V(136)
  A(162) = RCT(162)*V(130)*V(133)
  A(163) = RCT(163)*V(130)*V(135)
  A(164) = RCT(164)*V(130)*V(131)
  A(165) = RCT(165)*V(129)*V(130)
  A(166) = RCT(166)*V(128)*V(130)
  A(167) = RCT(167)*V(130)*V(130)
  A(168) = RCT(168)*V(130)*V(134)*F(2)
  A(169) = RCT(169)*V(77)*F(2)
  A(170) = RCT(170)*V(20)*V(136)
  A(171) = RCT(171)*V(22)*V(136)
  A(172) = RCT(172)*V(58)*V(133)
  A(173) = RCT(173)*V(58)*V(131)
  A(174) = RCT(174)*V(27)*V(136)
  A(175) = RCT(175)*V(59)*V(134)
  A(176) = RCT(176)*V(59)*V(132)
  A(177) = RCT(177)*V(81)*V(133)
  A(178) = RCT(178)*V(81)*V(131)
  A(179) = RCT(179)*V(28)*V(136)
  A(180) = RCT(180)*V(66)*V(133)
  A(181) = RCT(181)*V(66)*V(131)
  A(182) = RCT(182)*V(29)*V(136)
  A(183) = RCT(183)*V(87)*V(134)
  A(184) = RCT(184)*V(87)*V(133)
  A(185) = RCT(185)*V(87)*V(131)
  A(186) = RCT(186)*V(24)*V(136)
  A(187) = RCT(187)*V(25)*V(136)
  A(188) = RCT(188)*V(61)*V(131)
  A(189) = RCT(189)*V(32)*V(136)
  A(190) = RCT(190)*V(61)*V(133)
  A(191) = RCT(191)*V(33)*V(136)
  A(192) = RCT(192)*V(74)*V(134)
  A(193) = RCT(193)*V(17)
  A(194) = RCT(194)*V(74)*V(133)
  A(195) = RCT(195)*V(74)*V(131)
  A(196) = RCT(196)*V(88)*V(133)
  A(197) = RCT(197)*V(91)*V(131)
  A(198) = RCT(198)*V(91)*V(133)
  A(199) = RCT(199)*V(94)*V(131)
  A(200) = RCT(200)*V(94)*V(133)
  A(201) = RCT(201)*V(91)*V(134)
  A(202) = RCT(202)*V(94)*V(134)
  A(203) = RCT(203)*V(37)*V(136)
  A(204) = RCT(204)*V(38)*V(136)
  A(205) = RCT(205)*V(69)*V(133)
  A(206) = RCT(206)*V(69)*V(131)
  A(207) = RCT(207)*V(34)*V(136)
  A(208) = RCT(208)*V(90)*V(131)
  A(209) = RCT(209)*V(35)*V(136)
  A(210) = RCT(210)*V(90)*V(133)
  A(211) = RCT(211)*V(82)*V(136)
  A(212) = RCT(212)*V(83)*V(136)
  A(213) = RCT(213)*V(101)*V(136)
  A(214) = RCT(214)*V(104)*V(136)
  A(215) = RCT(215)*V(100)*V(136)
  A(216) = RCT(216)*V(82)*V(132)
  A(217) = RCT(217)*V(83)*V(132)
  A(218) = RCT(218)*V(101)*V(132)
  A(219) = RCT(219)*V(104)*V(132)
  A(220) = RCT(220)*V(100)*V(132)
  A(221) = RCT(221)*V(101)*V(135)
  A(222) = RCT(222)*V(104)*V(135)
  A(223) = RCT(223)*V(100)*V(135)
  A(224) = RCT(224)*V(118)*V(133)
  A(225) = RCT(225)*V(118)*V(131)
  A(226) = RCT(226)*V(118)*V(129)
  A(227) = RCT(227)*V(65)*V(136)
  A(228) = RCT(228)*V(79)*V(136)
  A(229) = RCT(229)*V(110)*V(136)
  A(230) = RCT(230)*V(110)*V(135)
  A(231) = RCT(231)*V(120)*V(133)
  A(232) = RCT(232)*V(120)*V(131)
  A(233) = RCT(233)*V(120)*V(129)
  A(234) = RCT(234)*V(114)*V(136)
  A(235) = RCT(235)*V(119)*V(133)
  A(236) = RCT(236)*V(119)*V(131)
  A(237) = RCT(237)*V(119)*V(129)
  A(238) = RCT(238)*V(119)*V(135)
  A(239) = RCT(239)*V(63)*V(136)
  A(240) = RCT(240)*V(106)*V(135)
  A(241) = RCT(241)*V(99)*V(133)
  A(242) = RCT(242)*V(99)*V(135)
  A(243) = RCT(243)*V(99)*V(131)
  A(244) = RCT(244)*V(117)*V(136)
  A(245) = RCT(245)*V(117)*V(135)
  A(246) = RCT(246)*V(116)*V(136)
  A(247) = RCT(247)*V(116)*V(135)
  A(248) = RCT(248)*V(53)*V(136)
  A(249) = RCT(249)*V(123)*V(133)
  A(250) = RCT(250)*V(123)*V(135)
  A(251) = RCT(251)*V(123)*V(131)
  A(252) = RCT(252)*V(123)*V(129)
  A(253) = RCT(253)*V(123)*V(128)
  A(254) = RCT(254)*V(26)*V(136)
  A(255) = RCT(255)*V(26)*V(136)
  A(256) = RCT(256)*V(95)*V(136)
  A(257) = RCT(257)*V(41)*V(136)
  A(258) = RCT(258)*V(77)*V(136)
  A(259) = RCT(259)*V(68)*V(136)
  A(260) = RCT(260)*V(122)*V(136)
  A(261) = RCT(261)*V(108)*V(136)
  A(262) = RCT(262)*V(39)*V(136)
  A(263) = RCT(263)*V(39)*V(136)
  A(264) = RCT(264)*V(39)*V(135)
  A(265) = RCT(265)*V(11)*V(136)
  A(266) = RCT(266)*V(131)
  A(267) = RCT(267)*V(97)*V(97)
  A(268) = RCT(268)*V(89)*V(136)
  A(269) = RCT(269)*V(111)*V(133)
  A(270) = RCT(270)*V(111)*V(129)
  A(271) = RCT(271)*V(103)*V(136)
  A(272) = RCT(272)*V(73)*V(133)
  A(273) = RCT(273)*V(73)*V(131)
  A(274) = RCT(274)*V(111)*V(131)
  A(275) = RCT(275)*V(76)*V(136)
  A(276) = RCT(276)*V(89)*V(132)
  A(277) = RCT(277)*V(89)*V(135)
  A(278) = RCT(278)*V(102)*V(131)
  A(279) = RCT(279)*V(102)*V(133)
  A(280) = RCT(280)*V(102)*V(135)
  A(281) = RCT(281)*V(19)*V(136)
  A(282) = 4.5e-13*V(78)*V(136)
  A(283) = RCT(283)*V(126)*V(131)
  A(284) = RCT(284)*V(56)
  A(285) = RCT(285)*V(56)*V(133)
  A(286) = RCT(286)*V(56)*V(131)
  A(287) = RCT(287)*V(8)*V(136)
  A(288) = RCT(288)*V(10)*V(136)
  A(289) = RCT(289)*V(106)*V(136)
  A(290) = RCT(290)*V(82)*V(136)
  A(291) = RCT(291)*V(83)*V(136)
  A(292) = RCT(292)*V(101)*V(136)


  Vdot(1) = -A(49)
  Vdot(2) = A(95)
  Vdot(3) = A(85)+A(97)+A(101)+A(108)+A(113)+A(121)+A(126)+A(131)+A(133)+A(138)+A(145)+A(155)+A(156)+A(162)+A(172)&
              &+A(177)+A(180)+A(184)+A(190)+A(194)+A(196)+A(198)+A(200)+A(205)+A(210)+A(224)+A(231)+A(235)+A(241)+A(249)&
              &+A(269)+A(272)+A(279)+A(285)
  Vdot(4) = A(86)+A(87)+A(88)+A(102)+A(103)+A(109)+A(115)+A(116)+A(119)+A(122)+A(123)+A(127)+A(128)+A(134)+A(139)+A(141)&
              &+A(147)+A(149)+A(150)+A(158)+A(159)+A(160)+A(164)+A(165)+A(166)+A(167)+A(173)+A(178)+A(181)+A(185)+A(188)&
              &+A(195)+A(197)+A(199)+A(206)+A(208)+A(225)+A(226)+A(232)+A(233)+A(236)+A(237)+A(243)+A(251)+A(252)+A(253)&
              &+A(267)+A(270)+A(273)+A(274)+A(278)+A(286)
  Vdot(5) = A(289)
  Vdot(6) = A(290)+A(291)+A(292)
  Vdot(7) = A(287)
  Vdot(8) = -A(287)
  Vdot(9) = A(288)
  Vdot(10) = -A(288)
  Vdot(11) = -A(265)
  Vdot(12) = -A(100)
  Vdot(13) = -A(120)
  Vdot(14) = -A(130)
  Vdot(15) = -A(4)-A(65)-A(66)
  Vdot(16) = -A(95)+A(262)+0.5*A(263)+A(264)
  Vdot(17) = A(192)-A(193)
  Vdot(18) = -A(13)+A(60)-A(61)+A(64)+0.5*A(266)
  Vdot(19) = -A(281)
  Vdot(20) = -A(170)
  Vdot(21) = -A(38)+0.12*A(170)
  Vdot(22) = 0.53*A(170)-A(171)
  Vdot(23) = A(97)-A(98)-A(99)
  Vdot(24) = -A(186)
  Vdot(25) = 0.18*A(186)-A(187)
  Vdot(26) = -A(26)+A(251)-A(254)-A(255)
  Vdot(27) = A(173)-A(174)
  Vdot(28) = A(178)-A(179)
  Vdot(29) = A(181)-A(182)
  Vdot(30) = -A(40)+0.2*A(196)+0.2*A(210)
  Vdot(31) = -A(39)+0.1*A(186)+0.23*A(203)
  Vdot(32) = A(188)-A(189)
  Vdot(33) = A(190)-A(191)
  Vdot(34) = A(206)-A(207)
  Vdot(35) = A(208)-A(209)
  Vdot(36) = -A(6)+A(72)-A(73)-A(80)
  Vdot(37) = -A(203)
  Vdot(38) = 0.15*A(203)-A(204)
  Vdot(39) = -A(262)-A(263)-A(264)
  Vdot(40) = -A(42)+0.21*A(210)
  Vdot(41) = 0.2*A(103)-A(257)+0.4*A(267)
  Vdot(42) = -A(41)+0.2*A(196)+0.15*A(210)
  Vdot(43) = -A(37)+A(38)+0.5*A(180)+0.2*A(196)+0.06*A(210)
  Vdot(44) = A(12)-A(54)-A(55)+0.05*A(84)
  Vdot(45) = -A(21)+A(102)-A(104)
  Vdot(46) = -A(22)+A(122)-A(124)
  Vdot(47) = -A(23)+A(127)-A(129)
  Vdot(48) = A(130)-A(131)
  Vdot(49) = -A(35)+A(139)-A(140)
  Vdot(50) = -A(10)+A(88)-A(89)
  Vdot(51) = A(158)-A(161)
  Vdot(52) = -A(9)+A(77)-A(78)-A(79)
  Vdot(53) = 0.33*A(145)+0.36*A(146)+0.3*A(149)+0.4*A(150)+A(246)+A(247)-A(248)
  Vdot(54) = -A(132)+0.05*A(216)+0.05*A(217)+0.05*A(218)+0.05*A(219)+0.05*A(220)
  Vdot(55) = -A(93)-A(94)
  Vdot(56) = A(283)-A(284)-A(285)-A(286)
  Vdot(57) = -A(36)+A(141)-A(142)
  Vdot(58) = 0.14*A(171)-A(172)-A(173)+A(174)+0.2*A(187)
  Vdot(59) = 0.06*A(171)-A(175)-A(176)+A(177)+0.07*A(187)+0.07*A(204)
  Vdot(60) = 0.75*A(93)-A(97)
  Vdot(61) = 0.07*A(186)-A(188)+A(189)-A(190)+0.06*A(203)
  Vdot(62) = -A(31)+0.8*A(34)+0.75*A(133)-A(137)
  Vdot(63) = 0.25*A(106)+0.25*A(115)+0.1*A(116)+0.2*A(144)+0.25*A(164)-A(239)
  Vdot(64) = -A(16)+0.75*A(115)-A(117)+0.75*A(164)
  Vdot(65) = -A(47)+A(225)-A(227)
  Vdot(66) = 0.35*A(170)-A(180)-A(181)+A(182)
  Vdot(67) = A(107)+0.1*A(133)-A(136)+0.75*A(279)+0.75*A(280)
  Vdot(68) = -A(17)+A(114)-A(118)-A(259)
  Vdot(69) = 0.3*A(204)-A(205)-A(206)+A(207)
  Vdot(70) = -A(15)+A(109)-A(110)
  Vdot(71) = -A(7)+A(74)-A(75)+2*A(80)+A(81)+0.5*A(82)+A(90)+A(112)+A(245)+A(264)
  Vdot(72) = -A(83)-A(84)+0.08*A(106)
  Vdot(73) = A(271)-A(272)-A(273)
  Vdot(74) = A(191)-A(192)+A(193)-A(194)-A(195)
  Vdot(75) = -A(28)+A(147)-A(148)
  Vdot(76) = -A(43)+A(274)-A(275)
  Vdot(77) = -A(18)+A(168)-A(169)-A(258)
  Vdot(78) = 0.5*A(94)+0.05*A(216)+0.05*A(217)+0.05*A(218)+0.05*A(219)+0.05*A(220)+0.25*A(276)+0.35*A(281)-A(282)+A(285)&
               &+A(286)
  Vdot(79) = -A(48)-A(228)+A(232)
  Vdot(80) = A(2)+A(4)-A(52)-A(53)-A(54)-A(65)-A(66)-A(84)
  Vdot(81) = A(176)-A(177)-A(178)+A(179)+A(194)+0.4*A(195)
  Vdot(82) = -A(211)-A(216)
  Vdot(83) = -A(212)-A(217)
  Vdot(84) = -A(34)+A(134)-A(135)
  Vdot(85) = A(137)-A(138)-A(139)+A(140)
  Vdot(86) = A(105)-A(108)-A(109)+0.5*A(110)
  Vdot(87) = 0.6*A(37)-A(183)-A(184)-A(185)
  Vdot(88) = -A(141)+A(142)+0.65*A(186)-A(196)
  Vdot(89) = -A(268)-A(276)-A(277)
  Vdot(90) = 0.56*A(203)-A(208)+A(209)-A(210)
  Vdot(91) = 0.6*A(40)-A(197)-A(198)-A(201)
  Vdot(92) = 2*A(1)+A(3)+A(5)-A(50)-A(51)+A(52)-A(56)-A(57)+A(63)-A(69)
  Vdot(93) = A(120)-A(121)-A(122)-A(123)+A(124)
  Vdot(94) = 0.6*A(41)-A(199)-A(200)-A(202)
  Vdot(95) = A(87)+0.3*A(103)+0.5*A(128)+0.25*A(149)+0.25*A(159)+0.25*A(226)+0.25*A(233)+0.25*A(237)+0.3*A(252)-A(256)&
               &+0.25*A(270)
  Vdot(96) = 0.13*A(32)-A(33)+0.45*A(36)-A(96)+0.02*A(145)+0.02*A(146)+0.7*A(172)+A(180)+0.4*A(184)+0.16*A(185)+0.6&
               &*A(196)+0.07*A(199)+0.17*A(200)+0.17*A(205)+0.34*A(210)+0.25*A(249)+0.25*A(250)+0.1*A(252)+0.25*A(253)+0.2&
               &*A(261)+0.65*A(281)
  Vdot(97) = A(31)+A(100)-A(101)-A(102)-A(103)+0.5*A(104)-2*A(267)
  Vdot(98) = -A(32)+0.9*A(36)+0.1*A(216)+0.1*A(217)+0.1*A(218)+0.1*A(219)+0.1*A(220)
  Vdot(99) = A(240)-A(241)-A(242)-A(243)
  Vdot(100) = -A(215)-A(220)-A(223)
  Vdot(101) = -A(213)-A(218)-A(221)
  Vdot(102) = A(277)-A(278)-A(279)-A(280)
  Vdot(103) = 0.33*A(43)-A(44)+0.33*A(269)+0.167*A(270)-A(271)+0.9*A(276)+0.25*A(279)+0.25*A(280)
  Vdot(104) = -A(214)-A(219)-A(222)
  Vdot(105) = 0.82*A(22)-A(24)+0.25*A(34)+0.67*A(43)+A(44)+0.5*A(46)+0.05*A(47)+0.3*A(48)+0.82*A(121)+0.82*A(123)-A(125)&
                &+0.5*A(131)+0.25*A(133)+0.52*A(216)+0.52*A(217)+0.52*A(218)+0.52*A(219)+0.52*A(220)+0.045*A(224)+0.025&
                &*A(226)+0.27*A(231)+0.15*A(233)+0.5*A(234)+0.67*A(269)+0.333*A(270)+A(272)+0.4*A(273)+0.1*A(276)+0.75&
                &*A(279)+0.75*A(280)
  Vdot(106) = -A(143)-A(144)-A(240)
  Vdot(107) = A(11)+A(12)+A(14)+0.67*A(19)+0.7*A(20)+A(25)+A(27)+A(30)+0.45*A(32)+2*A(33)+1.5*A(38)+1.5*A(39)+0.6*A(41)&
                &+A(42)+A(44)+A(45)+1.7*A(46)+0.25*A(48)+A(90)+A(91)-A(92)+0.5*A(94)+A(96)+0.56*A(106)+0.3*A(144)+0.05&
                &*A(152)+0.2*A(154)+0.22*A(155)+0.22*A(157)+0.11*A(159)+0.22*A(160)+0.4*A(184)+0.16*A(185)+0.07*A(197)+0.17&
                &*A(198)+0.14*A(199)+0.35*A(200)+0.23*A(216)+0.23*A(217)+0.23*A(218)+0.23*A(219)+0.23*A(220)+0.225*A(231)&
                &+0.125*A(233)+0.7*A(234)+A(244)+A(245)+0.25*A(249)+0.25*A(250)+0.2*A(252)+0.25*A(253)+0.35*A(276)+0.35&
                &*A(281)
  Vdot(108) = -A(30)+0.67*A(43)+0.25*A(48)+A(98)+0.02*A(145)+0.02*A(146)+0.53*A(155)+0.53*A(157)+0.26*A(159)+0.53*A(160)&
                &+0.225*A(231)+0.125*A(233)+0.25*A(249)+0.25*A(250)+0.1*A(252)+0.25*A(253)-A(261)+0.67*A(269)+0.333*A(270)
  Vdot(109) = 0.7*A(20)-A(105)-A(106)-A(107)+0.07*A(144)
  Vdot(110) = -A(45)+A(47)+0.33*A(216)+0.33*A(217)+0.33*A(218)+0.33*A(219)+0.33*A(220)+0.9*A(224)+A(226)-A(229)-A(230)&
                &+0.9*A(235)+0.5*A(237)+A(238)
  Vdot(111) = A(268)-A(269)-A(270)-A(274)+0.5*A(275)
  Vdot(112) = -A(14)+A(15)+A(21)+0.4*A(34)+A(35)+A(101)+0.8*A(103)+0.5*A(104)+0.5*A(106)+A(108)-A(111)-A(112)+0.27&
                &*A(121)+A(131)+0.4*A(133)+A(138)+0.04*A(152)+A(257)+1.6*A(267)
  Vdot(113) = A(132)-A(133)-A(134)+A(135)
  Vdot(114) = A(45)-A(46)+A(48)+0.3*A(216)+0.3*A(217)+0.3*A(218)+0.3*A(219)+0.3*A(220)+0.9*A(231)+A(233)-A(234)
  Vdot(115) = -A(19)+0.288*A(28)+0.4*A(144)+0.23*A(145)+0.25*A(146)+0.19*A(149)+0.25*A(150)-A(153)-A(154)+0.167*A(241)&
                &+0.167*A(242)+0.167*A(243)
  Vdot(116) = -A(27)+0.08*A(145)+0.8*A(156)+0.1*A(224)+0.1*A(231)+0.1*A(235)+A(236)+0.5*A(237)+0.794*A(241)+0.794*A(242)&
                &+0.794*A(243)-A(246)-A(247)
  Vdot(117) = -A(25)+0.18*A(32)+0.45*A(36)+A(42)+0.5*A(128)+A(136)+0.02*A(145)+0.02*A(146)+0.95*A(152)+0.8*A(154)+0.25&
                &*A(155)+0.25*A(157)+0.24*A(159)+0.25*A(160)+0.4*A(196)+0.07*A(197)+0.17*A(198)+0.07*A(199)+0.17*A(200)+0.51&
                &*A(205)+0.54*A(210)-A(244)-A(245)+0.25*A(249)+0.25*A(250)+0.1*A(252)+0.25*A(253)+A(260)
  Vdot(118) = A(211)+A(212)+A(213)+A(214)+A(215)-A(224)-A(225)-A(226)+A(227)
  Vdot(119) = A(221)+A(222)+A(223)+0.5*A(230)-A(235)-A(236)-A(237)-A(238)
  Vdot(120) = A(228)+A(229)+0.5*A(230)-A(231)-A(232)-A(233)
  Vdot(121) = 0.15*A(46)+A(125)-A(126)-A(127)-A(128)+A(129)+0.06*A(216)+0.06*A(217)+0.06*A(218)+0.06*A(219)+0.06*A(220)&
                &+0.15*A(234)
  Vdot(122) = -A(29)+0.5*A(110)+0.2*A(128)+0.02*A(145)+0.02*A(146)+0.22*A(155)+0.22*A(157)+0.23*A(159)+0.22*A(160)+0.25&
                &*A(249)+0.25*A(250)+0.1*A(252)+0.25*A(253)+0.5*A(258)-A(260)
  Vdot(123) = 0.5*A(148)+A(248)-A(249)-A(250)-A(251)-A(252)-A(253)+A(254)
  Vdot(124) = -A(20)+0.402*A(28)+0.2*A(144)+0.32*A(145)+0.35*A(146)+0.26*A(149)+0.35*A(150)-A(151)-A(152)+0.039*A(241)&
                &+0.039*A(242)+0.039*A(243)
  Vdot(125) = A(151)+0.5*A(153)-A(155)-A(156)-A(157)-A(158)-A(159)-A(160)+0.2*A(161)
  Vdot(126) = A(10)-A(11)-A(12)+A(15)+0.67*A(19)+A(23)+A(27)+0.69*A(28)+A(29)+A(30)+0.1*A(34)+0.33*A(43)+0.68*A(46)+0.4&
                &*A(47)+0.375*A(48)+0.25*A(84)+A(85)+2*A(86)+A(87)+0.3*A(89)-A(90)-A(91)+0.5*A(93)+A(94)+2*A(99)+0.7*A(103)&
                &+0.54*A(106)+A(108)+A(116)+0.5*A(117)+A(123)+A(126)+0.8*A(128)+0.5*A(131)+0.1*A(133)+0.6*A(144)+0.55*A(145)&
                &+0.6*A(146)+1.2*A(149)+0.6*A(150)+0.8*A(152)+0.7*A(154)+0.25*A(155)+0.25*A(157)+0.88*A(159)+0.25*A(160)&
                &+A(162)+A(163)+2*A(165)+A(166)+2*A(167)+0.34*A(216)+0.34*A(217)+0.34*A(218)+0.34*A(219)+0.34*A(220)+0.36&
                &*A(224)+0.95*A(226)+0.34*A(231)+0.93*A(233)+0.68*A(234)+0.75*A(237)+0.072*A(241)+0.072*A(242)+0.206*A(243)&
                &+0.25*A(249)+0.25*A(250)+0.8*A(252)+0.25*A(253)+A(256)+0.5*A(258)+A(259)+0.8*A(261)+0.33*A(269)+0.917&
                &*A(270)+0.5*A(276)+0.25*A(279)+0.25*A(280)-A(283)+A(284)
  Vdot(127) = A(143)-A(145)-A(146)-A(147)+0.5*A(148)-A(149)-A(150)
  Vdot(128) = 0.6*A(17)+0.67*A(19)+0.3*A(20)+A(23)+A(24)+A(25)+A(29)+A(31)+0.13*A(32)+A(35)+0.5*A(39)+A(42)+0.65*A(46)&
                &+A(111)+A(112)-A(113)-A(114)-A(115)-A(116)+0.5*A(117)+A(118)-2*A(119)+A(126)+0.3*A(128)+A(138)-A(150)+0.53&
                &*A(155)+0.53*A(157)+0.26*A(159)-0.47*A(160)+A(162)+A(163)+A(165)+2*A(167)+0.06*A(216)+0.06*A(217)+0.06&
                &*A(218)+0.06*A(219)+0.06*A(220)+0.65*A(234)+A(244)+A(245)-A(253)
  Vdot(129) = A(14)+A(16)+0.4*A(17)+0.3*A(20)+A(24)+A(83)+0.75*A(84)-A(85)-2*A(86)-2*A(87)-A(88)+0.7*A(89)-A(103)+0.31&
                &*A(106)+A(113)-0.1*A(116)+2*A(119)-A(123)-A(128)-A(149)+A(150)-A(159)+A(160)-A(165)+A(166)+0.33*A(197)+0.83&
                &*A(198)+0.07*A(199)+0.17*A(200)-A(226)-A(233)-A(237)+A(239)-A(252)+A(253)-A(270)
  Vdot(130) = A(18)+0.33*A(19)+0.2*A(144)+0.5*A(153)+0.5*A(161)-A(162)-A(163)-A(164)-A(165)-A(166)-2*A(167)-A(168)&
                &+A(169)
  Vdot(131) = 0.66*A(9)+A(10)+2*A(11)+A(14)+A(15)+0.67*A(19)+A(21)+A(22)+A(25)+A(27)+A(28)+A(29)+2*A(30)+0.56*A(32)+2&
                &*A(33)+0.9*A(34)+A(37)+1.5*A(38)+A(39)+0.6*A(40)+0.6*A(41)+A(42)+A(43)+2*A(44)+A(45)+1.2*A(46)+A(47)+A(48)&
                &+A(54)+A(55)+A(56)-A(57)+A(58)-A(59)-2*A(60)+A(61)-A(62)-A(67)-A(71)-A(77)+A(79)+0.4*A(84)+A(85)+2*A(86)&
                &-A(88)+A(90)+A(91)+A(92)+0.25*A(93)+0.12*A(94)+A(96)+A(98)+A(99)+A(101)-A(102)+A(103)+0.19*A(106)+A(108)&
                &-A(109)-A(115)+0.9*A(116)+A(121)-A(122)+A(123)-A(127)+0.3*A(128)+A(131)+0.9*A(133)-A(134)-A(139)-A(141)&
                &+0.06*A(144)+0.92*A(145)+A(146)-A(147)+A(149)+A(150)+0.06*A(152)+0.275*A(154)+0.47*A(155)+0.47*A(157)&
                &-A(158)+0.73*A(159)+0.47*A(160)+0.2*A(161)-A(164)+A(165)+0.65*A(170)+0.8*A(171)+A(172)-A(173)-A(178)+A(180)&
                &-A(181)+0.4*A(184)-0.84*A(185)+0.28*A(186)+0.73*A(187)-A(188)+A(190)-A(195)+A(196)-0.93*A(197)+0.17*A(198)&
                &-0.67*A(199)+0.83*A(200)+0.38*A(203)+0.63*A(204)+A(205)-A(206)-A(208)+A(210)+0.57*A(216)+0.57*A(217)+0.57&
                &*A(218)+0.57*A(219)+0.57*A(220)+0.9*A(224)-A(225)+A(226)+0.9*A(231)-A(232)+A(233)+0.2*A(234)-A(236)+0.5&
                &*A(237)+0.794*A(241)+0.794*A(242)-A(243)+A(246)+A(247)+A(249)+A(250)-A(251)+0.8*A(252)+A(253)+A(256)+A(257)&
                &+0.5*A(258)+A(260)+A(261)+0.5*A(263)-A(266)+1.2*A(267)+A(269)+A(270)+A(272)-0.6*A(273)-A(274)+0.06*A(276)&
                &-A(278)+0.75*A(279)+0.75*A(280)+0.35*A(281)+A(282)-A(283)+A(284)+A(285)-A(286)
  Vdot(132) = -A(2)-A(3)+0.89*A(8)+A(50)-A(51)-A(58)-A(59)-A(68)-A(70)-A(94)-A(106)+0.25*A(115)-0.9*A(144)-0.8*A(152)&
                &-0.8*A(154)+0.25*A(164)-A(176)-A(216)-A(217)-A(218)-A(219)-A(220)-A(276)
  Vdot(133) = A(5)+0.11*A(8)+A(49)+2*A(65)-A(67)-A(68)+A(69)-A(76)+0.5*A(82)-A(85)-A(97)-A(101)-A(108)-A(113)-A(121)&
                &-A(126)-A(131)-A(133)-A(138)-A(145)-A(155)-A(156)-A(162)-A(172)-A(177)-A(180)-A(184)-A(190)-A(194)-A(196)&
                &-A(198)-A(200)-A(205)-A(210)-A(224)-A(231)-A(235)-A(241)-A(249)-A(269)-A(272)-A(279)-A(285)
  Vdot(134) = -A(5)+A(6)+A(7)+0.89*A(8)+0.66*A(9)+0.6*A(17)+A(18)+A(27)+A(67)+A(68)-A(69)-A(70)+A(71)-A(72)+A(73)-A(74)&
                &+2*A(76)-A(77)+A(78)+A(79)-A(82)+A(85)+A(97)+A(101)+A(108)+A(113)-A(114)+A(118)+A(121)+A(126)+A(131)+0.9&
                &*A(133)+A(136)+A(138)+0.92*A(145)+A(146)+A(155)+A(157)+A(162)+A(163)-A(168)+A(169)+A(172)-A(175)+A(177)&
                &+A(180)-A(183)+A(190)-A(192)+A(193)+A(194)+A(196)+A(198)+A(200)-A(201)-A(202)+A(205)+A(210)+0.9*A(224)+0.9&
                &*A(231)+1.9*A(235)+0.5*A(237)+2*A(238)+1.206*A(241)+1.206*A(242)+0.206*A(243)+0.4*A(246)+A(247)+A(249)&
                &+A(250)+A(269)+A(272)+1.25*A(279)+1.25*A(280)+A(285)
  Vdot(135) = A(6)-A(8)+0.33*A(9)+0.4*A(17)+A(70)-A(71)-A(72)+A(73)+A(75)-A(76)-A(81)-A(90)-A(107)-A(112)-A(146)-A(157)&
                &-A(163)-A(221)-A(222)-A(223)-A(230)-A(238)-A(240)-A(242)-A(245)-A(247)-A(250)+0.5*A(258)+A(259)-A(264)&
                &-A(277)-A(280)
  Vdot(136) = A(7)+0.33*A(9)+A(10)+2*A(13)+A(15)+A(16)+0.33*A(19)+A(21)+A(22)+A(23)+A(26)+A(34)+A(35)+A(36)+A(43)+A(47)&
                &+A(48)+A(49)+2*A(53)+A(54)-A(55)-A(56)+A(57)-A(58)+A(59)-A(61)-A(62)-2*A(63)-2*A(64)+A(67)+A(71)-A(74)&
                &-A(75)-A(78)+0.5*A(82)-A(83)+0.75*A(84)-0.7*A(89)-A(91)-A(92)-A(93)+0.12*A(94)-A(95)-A(96)-A(100)-0.5&
                &*A(104)-A(105)+0.33*A(106)-0.5*A(110)-A(111)-A(117)-A(120)-A(124)-A(125)-A(129)-A(130)-A(132)-A(135)-A(136)&
                &-A(137)-A(140)-A(142)-A(143)+0.27*A(144)-A(148)-A(151)+0.08*A(152)-A(153)+0.215*A(154)-0.9*A(161)-A(170)&
                &-A(171)-A(174)-A(179)-A(182)-A(186)-A(187)-A(189)-A(191)+0.4*A(195)+0.4*A(197)+0.4*A(199)-A(203)-A(204)&
                &-A(207)-A(209)-A(211)-A(212)-A(213)-A(214)-A(215)+0.63*A(216)+0.63*A(217)+0.63*A(218)+0.63*A(219)+0.63&
                &*A(220)-A(227)-A(228)-A(229)-A(234)-A(239)+0.206*A(243)-A(244)-A(246)-A(248)-A(254)-A(256)-A(257)-A(258)&
                &-A(259)-A(260)-A(261)-A(262)-A(263)-A(265)-A(268)-A(271)+0.4*A(273)-0.5*A(275)+0.06*A(276)-0.35*A(281)&
                &-A(282)
      
END SUBROUTINE mozart_mosaic_4bin_Fun
















SUBROUTINE mozart_mosaic_4bin_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*F(2)
  IRR(2) = RCT(2)*V(132)
  IRR(3) = RCT(3)*V(132)
  IRR(4) = RCT(4)*V(15)
  IRR(5) = RCT(5)*V(134)
  IRR(6) = RCT(6)*V(36)
  IRR(7) = RCT(7)*V(71)
  IRR(8) = RCT(8)*V(135)
  IRR(9) = RCT(9)*V(52)
  IRR(10) = RCT(10)*V(50)
  IRR(11) = RCT(11)*V(126)
  IRR(12) = RCT(12)*V(126)
  IRR(13) = RCT(13)*V(18)
  IRR(14) = RCT(14)*V(112)
  IRR(15) = RCT(15)*V(70)
  IRR(16) = RCT(16)*V(64)
  IRR(17) = RCT(17)*V(68)
  IRR(18) = RCT(18)*V(77)
  IRR(19) = RCT(19)*V(115)
  IRR(20) = RCT(20)*V(124)
  IRR(21) = RCT(21)*V(45)
  IRR(22) = RCT(22)*V(46)
  IRR(23) = RCT(23)*V(47)
  IRR(24) = RCT(24)*V(105)
  IRR(25) = RCT(25)*V(117)
  IRR(26) = RCT(26)*V(26)
  IRR(27) = RCT(27)*V(116)
  IRR(28) = RCT(28)*V(75)
  IRR(29) = RCT(29)*V(122)
  IRR(30) = RCT(30)*V(108)
  IRR(31) = RCT(31)*V(62)
  IRR(32) = RCT(32)*V(98)
  IRR(33) = RCT(33)*V(96)
  IRR(34) = RCT(34)*V(84)
  IRR(35) = RCT(35)*V(49)
  IRR(36) = RCT(36)*V(57)
  IRR(37) = RCT(37)*V(43)
  IRR(38) = RCT(38)*V(21)
  IRR(39) = RCT(39)*V(31)
  IRR(40) = RCT(40)*V(30)
  IRR(41) = RCT(41)*V(42)
  IRR(42) = RCT(42)*V(40)
  IRR(43) = RCT(43)*V(76)
  IRR(44) = RCT(44)*V(103)
  IRR(45) = RCT(45)*V(110)
  IRR(46) = RCT(46)*V(114)
  IRR(47) = RCT(47)*V(65)
  IRR(48) = RCT(48)*V(79)
  IRR(49) = RCT(49)*V(1)
  IRR(50) = RCT(50)*V(92)*F(2)
  IRR(51) = RCT(51)*V(92)*V(132)
  IRR(52) = RCT(52)*V(80)*F(2)
  IRR(53) = RCT(53)*V(80)*F(1)
  IRR(54) = RCT(54)*V(44)*V(80)
  IRR(55) = RCT(55)*V(44)*V(136)
  IRR(56) = RCT(56)*V(92)*V(136)
  IRR(57) = RCT(57)*V(92)*V(131)
  IRR(58) = RCT(58)*V(132)*V(136)
  IRR(59) = RCT(59)*V(131)*V(132)
  IRR(60) = RCT(60)*V(131)*V(131)*F(1)
  IRR(61) = RCT(61)*V(18)*V(136)
  IRR(62) = RCT(62)*V(131)*V(136)
  IRR(63) = RCT(63)*V(136)*V(136)
  IRR(64) = RCT(64)*V(136)*V(136)
  IRR(65) = RCT(65)*V(15)*V(80)
  IRR(66) = RCT(66)*V(15)*V(80)
  IRR(67) = RCT(67)*V(131)*V(133)
  IRR(68) = RCT(68)*V(132)*V(133)
  IRR(69) = RCT(69)*V(92)*V(134)
  IRR(70) = RCT(70)*V(132)*V(134)
  IRR(71) = RCT(71)*V(131)*V(135)
  IRR(72) = RCT(72)*V(134)*V(135)
  IRR(73) = RCT(73)*V(36)
  IRR(74) = RCT(74)*V(134)*V(136)
  IRR(75) = RCT(75)*V(71)*V(136)
  IRR(76) = RCT(76)*V(133)*V(135)
  IRR(77) = RCT(77)*V(131)*V(134)
  IRR(78) = RCT(78)*V(52)*V(136)
  IRR(79) = RCT(79)*V(52)
  IRR(80) = RCT(80)*V(36)*F(2)
  IRR(81) = RCT(81)*V(135)
  IRR(82) = RCT(82)*V(134)
  IRR(83) = RCT(83)*V(72)*V(136)
  IRR(84) = RCT(84)*V(72)*V(80)
  IRR(85) = RCT(85)*V(129)*V(133)
  IRR(86) = RCT(86)*V(129)*V(129)
  IRR(87) = RCT(87)*V(129)*V(129)
  IRR(88) = RCT(88)*V(129)*V(131)
  IRR(89) = RCT(89)*V(50)*V(136)
  IRR(90) = RCT(90)*V(126)*V(135)
  IRR(91) = RCT(91)*V(126)*V(136)
  IRR(92) = RCT(92)*V(107)*V(136)
  IRR(93) = RCT(93)*V(55)*V(136)
  IRR(94) = RCT(94)*V(55)*V(132)
  IRR(95) = RCT(95)*V(16)*V(136)
  IRR(96) = RCT(96)*V(96)*V(136)
  IRR(97) = RCT(97)*V(60)*V(133)
  IRR(98) = RCT(98)*V(23)*F(2)
  IRR(99) = RCT(99)*V(23)
  IRR(100) = RCT(100)*V(12)*V(136)
  IRR(101) = RCT(101)*V(97)*V(133)
  IRR(102) = RCT(102)*V(97)*V(131)
  IRR(103) = RCT(103)*V(97)*V(129)
  IRR(104) = RCT(104)*V(45)*V(136)
  IRR(105) = RCT(105)*V(109)*V(136)
  IRR(106) = RCT(106)*V(109)*V(132)
  IRR(107) = RCT(107)*V(109)*V(135)
  IRR(108) = RCT(108)*V(86)*V(133)
  IRR(109) = RCT(109)*V(86)*V(131)
  IRR(110) = RCT(110)*V(70)*V(136)
  IRR(111) = RCT(111)*V(112)*V(136)
  IRR(112) = RCT(112)*V(112)*V(135)
  IRR(113) = RCT(113)*V(128)*V(133)
  IRR(114) = RCT(114)*V(128)*V(134)
  IRR(115) = RCT(115)*V(128)*V(131)
  IRR(116) = RCT(116)*V(128)*V(129)
  IRR(117) = RCT(117)*V(64)*V(136)
  IRR(118) = RCT(118)*V(68)
  IRR(119) = RCT(119)*V(128)*V(128)
  IRR(120) = RCT(120)*V(13)*V(136)
  IRR(121) = RCT(121)*V(93)*V(133)
  IRR(122) = RCT(122)*V(93)*V(131)
  IRR(123) = RCT(123)*V(93)*V(129)
  IRR(124) = RCT(124)*V(46)*V(136)
  IRR(125) = RCT(125)*V(105)*V(136)
  IRR(126) = RCT(126)*V(121)*V(133)
  IRR(127) = RCT(127)*V(121)*V(131)
  IRR(128) = RCT(128)*V(121)*V(129)
  IRR(129) = RCT(129)*V(47)*V(136)
  IRR(130) = RCT(130)*V(14)*V(136)
  IRR(131) = RCT(131)*V(48)*V(133)
  IRR(132) = RCT(132)*V(54)*V(136)
  IRR(133) = RCT(133)*V(113)*V(133)
  IRR(134) = RCT(134)*V(113)*V(131)
  IRR(135) = RCT(135)*V(84)*V(136)
  IRR(136) = RCT(136)*V(67)*V(136)
  IRR(137) = RCT(137)*V(62)*V(136)
  IRR(138) = RCT(138)*V(85)*V(133)
  IRR(139) = RCT(139)*V(85)*V(131)
  IRR(140) = RCT(140)*V(49)*V(136)
  IRR(141) = RCT(141)*V(88)*V(131)
  IRR(142) = RCT(142)*V(57)*V(136)
  IRR(143) = RCT(143)*V(106)*V(136)
  IRR(144) = RCT(144)*V(106)*V(132)
  IRR(145) = RCT(145)*V(127)*V(133)
  IRR(146) = RCT(146)*V(127)*V(135)
  IRR(147) = RCT(147)*V(127)*V(131)
  IRR(148) = RCT(148)*V(75)*V(136)
  IRR(149) = RCT(149)*V(127)*V(129)
  IRR(150) = RCT(150)*V(127)*V(128)
  IRR(151) = RCT(151)*V(124)*V(136)
  IRR(152) = RCT(152)*V(124)*V(132)
  IRR(153) = RCT(153)*V(115)*V(136)
  IRR(154) = RCT(154)*V(115)*V(132)
  IRR(155) = RCT(155)*V(125)*V(133)
  IRR(156) = RCT(156)*V(125)*V(133)
  IRR(157) = RCT(157)*V(125)*V(135)
  IRR(158) = RCT(158)*V(125)*V(131)
  IRR(159) = RCT(159)*V(125)*V(129)
  IRR(160) = RCT(160)*V(125)*V(128)
  IRR(161) = RCT(161)*V(51)*V(136)
  IRR(162) = RCT(162)*V(130)*V(133)
  IRR(163) = RCT(163)*V(130)*V(135)
  IRR(164) = RCT(164)*V(130)*V(131)
  IRR(165) = RCT(165)*V(129)*V(130)
  IRR(166) = RCT(166)*V(128)*V(130)
  IRR(167) = RCT(167)*V(130)*V(130)
  IRR(168) = RCT(168)*V(130)*V(134)*F(2)
  IRR(169) = RCT(169)*V(77)*F(2)
  IRR(170) = RCT(170)*V(20)*V(136)
  IRR(171) = RCT(171)*V(22)*V(136)
  IRR(172) = RCT(172)*V(58)*V(133)
  IRR(173) = RCT(173)*V(58)*V(131)
  IRR(174) = RCT(174)*V(27)*V(136)
  IRR(175) = RCT(175)*V(59)*V(134)
  IRR(176) = RCT(176)*V(59)*V(132)
  IRR(177) = RCT(177)*V(81)*V(133)
  IRR(178) = RCT(178)*V(81)*V(131)
  IRR(179) = RCT(179)*V(28)*V(136)
  IRR(180) = RCT(180)*V(66)*V(133)
  IRR(181) = RCT(181)*V(66)*V(131)
  IRR(182) = RCT(182)*V(29)*V(136)
  IRR(183) = RCT(183)*V(87)*V(134)
  IRR(184) = RCT(184)*V(87)*V(133)
  IRR(185) = RCT(185)*V(87)*V(131)
  IRR(186) = RCT(186)*V(24)*V(136)
  IRR(187) = RCT(187)*V(25)*V(136)
  IRR(188) = RCT(188)*V(61)*V(131)
  IRR(189) = RCT(189)*V(32)*V(136)
  IRR(190) = RCT(190)*V(61)*V(133)
  IRR(191) = RCT(191)*V(33)*V(136)
  IRR(192) = RCT(192)*V(74)*V(134)
  IRR(193) = RCT(193)*V(17)
  IRR(194) = RCT(194)*V(74)*V(133)
  IRR(195) = RCT(195)*V(74)*V(131)
  IRR(196) = RCT(196)*V(88)*V(133)
  IRR(197) = RCT(197)*V(91)*V(131)
  IRR(198) = RCT(198)*V(91)*V(133)
  IRR(199) = RCT(199)*V(94)*V(131)
  IRR(200) = RCT(200)*V(94)*V(133)
  IRR(201) = RCT(201)*V(91)*V(134)
  IRR(202) = RCT(202)*V(94)*V(134)
  IRR(203) = RCT(203)*V(37)*V(136)
  IRR(204) = RCT(204)*V(38)*V(136)
  IRR(205) = RCT(205)*V(69)*V(133)
  IRR(206) = RCT(206)*V(69)*V(131)
  IRR(207) = RCT(207)*V(34)*V(136)
  IRR(208) = RCT(208)*V(90)*V(131)
  IRR(209) = RCT(209)*V(35)*V(136)
  IRR(210) = RCT(210)*V(90)*V(133)
  IRR(211) = RCT(211)*V(82)*V(136)
  IRR(212) = RCT(212)*V(83)*V(136)
  IRR(213) = RCT(213)*V(101)*V(136)
  IRR(214) = RCT(214)*V(104)*V(136)
  IRR(215) = RCT(215)*V(100)*V(136)
  IRR(216) = RCT(216)*V(82)*V(132)
  IRR(217) = RCT(217)*V(83)*V(132)
  IRR(218) = RCT(218)*V(101)*V(132)
  IRR(219) = RCT(219)*V(104)*V(132)
  IRR(220) = RCT(220)*V(100)*V(132)
  IRR(221) = RCT(221)*V(101)*V(135)
  IRR(222) = RCT(222)*V(104)*V(135)
  IRR(223) = RCT(223)*V(100)*V(135)
  IRR(224) = RCT(224)*V(118)*V(133)
  IRR(225) = RCT(225)*V(118)*V(131)
  IRR(226) = RCT(226)*V(118)*V(129)
  IRR(227) = RCT(227)*V(65)*V(136)
  IRR(228) = RCT(228)*V(79)*V(136)
  IRR(229) = RCT(229)*V(110)*V(136)
  IRR(230) = RCT(230)*V(110)*V(135)
  IRR(231) = RCT(231)*V(120)*V(133)
  IRR(232) = RCT(232)*V(120)*V(131)
  IRR(233) = RCT(233)*V(120)*V(129)
  IRR(234) = RCT(234)*V(114)*V(136)
  IRR(235) = RCT(235)*V(119)*V(133)
  IRR(236) = RCT(236)*V(119)*V(131)
  IRR(237) = RCT(237)*V(119)*V(129)
  IRR(238) = RCT(238)*V(119)*V(135)
  IRR(239) = RCT(239)*V(63)*V(136)
  IRR(240) = RCT(240)*V(106)*V(135)
  IRR(241) = RCT(241)*V(99)*V(133)
  IRR(242) = RCT(242)*V(99)*V(135)
  IRR(243) = RCT(243)*V(99)*V(131)
  IRR(244) = RCT(244)*V(117)*V(136)
  IRR(245) = RCT(245)*V(117)*V(135)
  IRR(246) = RCT(246)*V(116)*V(136)
  IRR(247) = RCT(247)*V(116)*V(135)
  IRR(248) = RCT(248)*V(53)*V(136)
  IRR(249) = RCT(249)*V(123)*V(133)
  IRR(250) = RCT(250)*V(123)*V(135)
  IRR(251) = RCT(251)*V(123)*V(131)
  IRR(252) = RCT(252)*V(123)*V(129)
  IRR(253) = RCT(253)*V(123)*V(128)
  IRR(254) = RCT(254)*V(26)*V(136)
  IRR(255) = RCT(255)*V(26)*V(136)
  IRR(256) = RCT(256)*V(95)*V(136)
  IRR(257) = RCT(257)*V(41)*V(136)
  IRR(258) = RCT(258)*V(77)*V(136)
  IRR(259) = RCT(259)*V(68)*V(136)
  IRR(260) = RCT(260)*V(122)*V(136)
  IRR(261) = RCT(261)*V(108)*V(136)
  IRR(262) = RCT(262)*V(39)*V(136)
  IRR(263) = RCT(263)*V(39)*V(136)
  IRR(264) = RCT(264)*V(39)*V(135)
  IRR(265) = RCT(265)*V(11)*V(136)
  IRR(266) = RCT(266)*V(131)
  IRR(267) = RCT(267)*V(97)*V(97)
  IRR(268) = RCT(268)*V(89)*V(136)
  IRR(269) = RCT(269)*V(111)*V(133)
  IRR(270) = RCT(270)*V(111)*V(129)
  IRR(271) = RCT(271)*V(103)*V(136)
  IRR(272) = RCT(272)*V(73)*V(133)
  IRR(273) = RCT(273)*V(73)*V(131)
  IRR(274) = RCT(274)*V(111)*V(131)
  IRR(275) = RCT(275)*V(76)*V(136)
  IRR(276) = RCT(276)*V(89)*V(132)
  IRR(277) = RCT(277)*V(89)*V(135)
  IRR(278) = RCT(278)*V(102)*V(131)
  IRR(279) = RCT(279)*V(102)*V(133)
  IRR(280) = RCT(280)*V(102)*V(135)
  IRR(281) = RCT(281)*V(19)*V(136)
  IRR(282) = 4.5e-13*V(78)*V(136)
  IRR(283) = RCT(283)*V(126)*V(131)
  IRR(284) = RCT(284)*V(56)
  IRR(285) = RCT(285)*V(56)*V(133)
  IRR(286) = RCT(286)*V(56)*V(131)
  IRR(287) = RCT(287)*V(8)*V(136)
  IRR(288) = RCT(288)*V(10)*V(136)
  IRR(289) = RCT(289)*V(106)*V(136)
  IRR(290) = RCT(290)*V(82)*V(136)
  IRR(291) = RCT(291)*V(83)*V(136)
  IRR(292) = RCT(292)*V(101)*V(136)
      
END SUBROUTINE mozart_mosaic_4bin_IRRFun
















SUBROUTINE mozart_mosaic_4bin_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(520)


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

  B(38) = RCT(38)

  B(39) = RCT(39)

  B(40) = RCT(40)

  B(41) = RCT(41)

  B(42) = RCT(42)

  B(43) = RCT(43)

  B(44) = RCT(44)

  B(45) = RCT(45)

  B(46) = RCT(46)

  B(47) = RCT(47)

  B(48) = RCT(48)

  B(49) = RCT(49)

  B(50) = RCT(50)*F(2)

  B(52) = RCT(51)*V(132)

  B(53) = RCT(51)*V(92)

  B(54) = RCT(52)*F(2)

  B(56) = RCT(53)*F(1)

  B(58) = RCT(54)*V(80)

  B(59) = RCT(54)*V(44)

  B(60) = RCT(55)*V(136)

  B(61) = RCT(55)*V(44)

  B(62) = RCT(56)*V(136)

  B(63) = RCT(56)*V(92)

  B(64) = RCT(57)*V(131)

  B(65) = RCT(57)*V(92)

  B(66) = RCT(58)*V(136)

  B(67) = RCT(58)*V(132)

  B(68) = RCT(59)*V(132)

  B(69) = RCT(59)*V(131)

  B(70) = RCT(60)*2*V(131)*F(1)

  B(72) = RCT(61)*V(136)

  B(73) = RCT(61)*V(18)

  B(74) = RCT(62)*V(136)

  B(75) = RCT(62)*V(131)

  B(76) = RCT(63)*2*V(136)

  B(77) = RCT(64)*2*V(136)

  B(78) = RCT(65)*V(80)

  B(79) = RCT(65)*V(15)

  B(80) = RCT(66)*V(80)

  B(81) = RCT(66)*V(15)

  B(82) = RCT(67)*V(133)

  B(83) = RCT(67)*V(131)

  B(84) = RCT(68)*V(133)

  B(85) = RCT(68)*V(132)

  B(86) = RCT(69)*V(134)

  B(87) = RCT(69)*V(92)

  B(88) = RCT(70)*V(134)

  B(89) = RCT(70)*V(132)

  B(90) = RCT(71)*V(135)

  B(91) = RCT(71)*V(131)

  B(92) = RCT(72)*V(135)

  B(93) = RCT(72)*V(134)

  B(94) = RCT(73)

  B(95) = RCT(74)*V(136)

  B(96) = RCT(74)*V(134)

  B(97) = RCT(75)*V(136)

  B(98) = RCT(75)*V(71)

  B(99) = RCT(76)*V(135)

  B(100) = RCT(76)*V(133)

  B(101) = RCT(77)*V(134)

  B(102) = RCT(77)*V(131)

  B(103) = RCT(78)*V(136)

  B(104) = RCT(78)*V(52)

  B(105) = RCT(79)

  B(106) = RCT(80)*F(2)

  B(108) = RCT(81)

  B(109) = RCT(82)

  B(110) = RCT(83)*V(136)

  B(111) = RCT(83)*V(72)

  B(112) = RCT(84)*V(80)

  B(113) = RCT(84)*V(72)

  B(114) = RCT(85)*V(133)

  B(115) = RCT(85)*V(129)

  B(116) = RCT(86)*2*V(129)

  B(117) = RCT(87)*2*V(129)

  B(118) = RCT(88)*V(131)

  B(119) = RCT(88)*V(129)

  B(120) = RCT(89)*V(136)

  B(121) = RCT(89)*V(50)

  B(122) = RCT(90)*V(135)

  B(123) = RCT(90)*V(126)

  B(124) = RCT(91)*V(136)

  B(125) = RCT(91)*V(126)

  B(126) = RCT(92)*V(136)

  B(127) = RCT(92)*V(107)

  B(128) = RCT(93)*V(136)

  B(129) = RCT(93)*V(55)

  B(130) = RCT(94)*V(132)

  B(131) = RCT(94)*V(55)

  B(132) = RCT(95)*V(136)

  B(133) = RCT(95)*V(16)

  B(134) = RCT(96)*V(136)

  B(135) = RCT(96)*V(96)

  B(136) = RCT(97)*V(133)

  B(137) = RCT(97)*V(60)

  B(138) = RCT(98)*F(2)

  B(140) = RCT(99)

  B(141) = RCT(100)*V(136)

  B(142) = RCT(100)*V(12)

  B(143) = RCT(101)*V(133)

  B(144) = RCT(101)*V(97)

  B(145) = RCT(102)*V(131)

  B(146) = RCT(102)*V(97)

  B(147) = RCT(103)*V(129)

  B(148) = RCT(103)*V(97)

  B(149) = RCT(104)*V(136)

  B(150) = RCT(104)*V(45)

  B(151) = RCT(105)*V(136)

  B(152) = RCT(105)*V(109)

  B(153) = RCT(106)*V(132)

  B(154) = RCT(106)*V(109)

  B(155) = RCT(107)*V(135)

  B(156) = RCT(107)*V(109)

  B(157) = RCT(108)*V(133)

  B(158) = RCT(108)*V(86)

  B(159) = RCT(109)*V(131)

  B(160) = RCT(109)*V(86)

  B(161) = RCT(110)*V(136)

  B(162) = RCT(110)*V(70)

  B(163) = RCT(111)*V(136)

  B(164) = RCT(111)*V(112)

  B(165) = RCT(112)*V(135)

  B(166) = RCT(112)*V(112)

  B(167) = RCT(113)*V(133)

  B(168) = RCT(113)*V(128)

  B(169) = RCT(114)*V(134)

  B(170) = RCT(114)*V(128)

  B(171) = RCT(115)*V(131)

  B(172) = RCT(115)*V(128)

  B(173) = RCT(116)*V(129)

  B(174) = RCT(116)*V(128)

  B(175) = RCT(117)*V(136)

  B(176) = RCT(117)*V(64)

  B(177) = RCT(118)

  B(178) = RCT(119)*2*V(128)

  B(179) = RCT(120)*V(136)

  B(180) = RCT(120)*V(13)

  B(181) = RCT(121)*V(133)

  B(182) = RCT(121)*V(93)

  B(183) = RCT(122)*V(131)

  B(184) = RCT(122)*V(93)

  B(185) = RCT(123)*V(129)

  B(186) = RCT(123)*V(93)

  B(187) = RCT(124)*V(136)

  B(188) = RCT(124)*V(46)

  B(189) = RCT(125)*V(136)

  B(190) = RCT(125)*V(105)

  B(191) = RCT(126)*V(133)

  B(192) = RCT(126)*V(121)

  B(193) = RCT(127)*V(131)

  B(194) = RCT(127)*V(121)

  B(195) = RCT(128)*V(129)

  B(196) = RCT(128)*V(121)

  B(197) = RCT(129)*V(136)

  B(198) = RCT(129)*V(47)

  B(199) = RCT(130)*V(136)

  B(200) = RCT(130)*V(14)

  B(201) = RCT(131)*V(133)

  B(202) = RCT(131)*V(48)

  B(203) = RCT(132)*V(136)

  B(204) = RCT(132)*V(54)

  B(205) = RCT(133)*V(133)

  B(206) = RCT(133)*V(113)

  B(207) = RCT(134)*V(131)

  B(208) = RCT(134)*V(113)

  B(209) = RCT(135)*V(136)

  B(210) = RCT(135)*V(84)

  B(211) = RCT(136)*V(136)

  B(212) = RCT(136)*V(67)

  B(213) = RCT(137)*V(136)

  B(214) = RCT(137)*V(62)

  B(215) = RCT(138)*V(133)

  B(216) = RCT(138)*V(85)

  B(217) = RCT(139)*V(131)

  B(218) = RCT(139)*V(85)

  B(219) = RCT(140)*V(136)

  B(220) = RCT(140)*V(49)

  B(221) = RCT(141)*V(131)

  B(222) = RCT(141)*V(88)

  B(223) = RCT(142)*V(136)

  B(224) = RCT(142)*V(57)

  B(225) = RCT(143)*V(136)

  B(226) = RCT(143)*V(106)

  B(227) = RCT(144)*V(132)

  B(228) = RCT(144)*V(106)

  B(229) = RCT(145)*V(133)

  B(230) = RCT(145)*V(127)

  B(231) = RCT(146)*V(135)

  B(232) = RCT(146)*V(127)

  B(233) = RCT(147)*V(131)

  B(234) = RCT(147)*V(127)

  B(235) = RCT(148)*V(136)

  B(236) = RCT(148)*V(75)

  B(237) = RCT(149)*V(129)

  B(238) = RCT(149)*V(127)

  B(239) = RCT(150)*V(128)

  B(240) = RCT(150)*V(127)

  B(241) = RCT(151)*V(136)

  B(242) = RCT(151)*V(124)

  B(243) = RCT(152)*V(132)

  B(244) = RCT(152)*V(124)

  B(245) = RCT(153)*V(136)

  B(246) = RCT(153)*V(115)

  B(247) = RCT(154)*V(132)

  B(248) = RCT(154)*V(115)

  B(249) = RCT(155)*V(133)

  B(250) = RCT(155)*V(125)

  B(251) = RCT(156)*V(133)

  B(252) = RCT(156)*V(125)

  B(253) = RCT(157)*V(135)

  B(254) = RCT(157)*V(125)

  B(255) = RCT(158)*V(131)

  B(256) = RCT(158)*V(125)

  B(257) = RCT(159)*V(129)

  B(258) = RCT(159)*V(125)

  B(259) = RCT(160)*V(128)

  B(260) = RCT(160)*V(125)

  B(261) = RCT(161)*V(136)

  B(262) = RCT(161)*V(51)

  B(263) = RCT(162)*V(133)

  B(264) = RCT(162)*V(130)

  B(265) = RCT(163)*V(135)

  B(266) = RCT(163)*V(130)

  B(267) = RCT(164)*V(131)

  B(268) = RCT(164)*V(130)

  B(269) = RCT(165)*V(130)

  B(270) = RCT(165)*V(129)

  B(271) = RCT(166)*V(130)

  B(272) = RCT(166)*V(128)

  B(273) = RCT(167)*2*V(130)

  B(274) = RCT(168)*V(134)*F(2)

  B(275) = RCT(168)*V(130)*F(2)

  B(277) = RCT(169)*F(2)

  B(279) = RCT(170)*V(136)

  B(280) = RCT(170)*V(20)

  B(281) = RCT(171)*V(136)

  B(282) = RCT(171)*V(22)

  B(283) = RCT(172)*V(133)

  B(284) = RCT(172)*V(58)

  B(285) = RCT(173)*V(131)

  B(286) = RCT(173)*V(58)

  B(287) = RCT(174)*V(136)

  B(288) = RCT(174)*V(27)

  B(289) = RCT(175)*V(134)

  B(290) = RCT(175)*V(59)

  B(291) = RCT(176)*V(132)

  B(292) = RCT(176)*V(59)

  B(293) = RCT(177)*V(133)

  B(294) = RCT(177)*V(81)

  B(295) = RCT(178)*V(131)

  B(296) = RCT(178)*V(81)

  B(297) = RCT(179)*V(136)

  B(298) = RCT(179)*V(28)

  B(299) = RCT(180)*V(133)

  B(300) = RCT(180)*V(66)

  B(301) = RCT(181)*V(131)

  B(302) = RCT(181)*V(66)

  B(303) = RCT(182)*V(136)

  B(304) = RCT(182)*V(29)

  B(305) = RCT(183)*V(134)

  B(306) = RCT(183)*V(87)

  B(307) = RCT(184)*V(133)

  B(308) = RCT(184)*V(87)

  B(309) = RCT(185)*V(131)

  B(310) = RCT(185)*V(87)

  B(311) = RCT(186)*V(136)

  B(312) = RCT(186)*V(24)

  B(313) = RCT(187)*V(136)

  B(314) = RCT(187)*V(25)

  B(315) = RCT(188)*V(131)

  B(316) = RCT(188)*V(61)

  B(317) = RCT(189)*V(136)

  B(318) = RCT(189)*V(32)

  B(319) = RCT(190)*V(133)

  B(320) = RCT(190)*V(61)

  B(321) = RCT(191)*V(136)

  B(322) = RCT(191)*V(33)

  B(323) = RCT(192)*V(134)

  B(324) = RCT(192)*V(74)

  B(325) = RCT(193)

  B(326) = RCT(194)*V(133)

  B(327) = RCT(194)*V(74)

  B(328) = RCT(195)*V(131)

  B(329) = RCT(195)*V(74)

  B(330) = RCT(196)*V(133)

  B(331) = RCT(196)*V(88)

  B(332) = RCT(197)*V(131)

  B(333) = RCT(197)*V(91)

  B(334) = RCT(198)*V(133)

  B(335) = RCT(198)*V(91)

  B(336) = RCT(199)*V(131)

  B(337) = RCT(199)*V(94)

  B(338) = RCT(200)*V(133)

  B(339) = RCT(200)*V(94)

  B(340) = RCT(201)*V(134)

  B(341) = RCT(201)*V(91)

  B(342) = RCT(202)*V(134)

  B(343) = RCT(202)*V(94)

  B(344) = RCT(203)*V(136)

  B(345) = RCT(203)*V(37)

  B(346) = RCT(204)*V(136)

  B(347) = RCT(204)*V(38)

  B(348) = RCT(205)*V(133)

  B(349) = RCT(205)*V(69)

  B(350) = RCT(206)*V(131)

  B(351) = RCT(206)*V(69)

  B(352) = RCT(207)*V(136)

  B(353) = RCT(207)*V(34)

  B(354) = RCT(208)*V(131)

  B(355) = RCT(208)*V(90)

  B(356) = RCT(209)*V(136)

  B(357) = RCT(209)*V(35)

  B(358) = RCT(210)*V(133)

  B(359) = RCT(210)*V(90)

  B(360) = RCT(211)*V(136)

  B(361) = RCT(211)*V(82)

  B(362) = RCT(212)*V(136)

  B(363) = RCT(212)*V(83)

  B(364) = RCT(213)*V(136)

  B(365) = RCT(213)*V(101)

  B(366) = RCT(214)*V(136)

  B(367) = RCT(214)*V(104)

  B(368) = RCT(215)*V(136)

  B(369) = RCT(215)*V(100)

  B(370) = RCT(216)*V(132)

  B(371) = RCT(216)*V(82)

  B(372) = RCT(217)*V(132)

  B(373) = RCT(217)*V(83)

  B(374) = RCT(218)*V(132)

  B(375) = RCT(218)*V(101)

  B(376) = RCT(219)*V(132)

  B(377) = RCT(219)*V(104)

  B(378) = RCT(220)*V(132)

  B(379) = RCT(220)*V(100)

  B(380) = RCT(221)*V(135)

  B(381) = RCT(221)*V(101)

  B(382) = RCT(222)*V(135)

  B(383) = RCT(222)*V(104)

  B(384) = RCT(223)*V(135)

  B(385) = RCT(223)*V(100)

  B(386) = RCT(224)*V(133)

  B(387) = RCT(224)*V(118)

  B(388) = RCT(225)*V(131)

  B(389) = RCT(225)*V(118)

  B(390) = RCT(226)*V(129)

  B(391) = RCT(226)*V(118)

  B(392) = RCT(227)*V(136)

  B(393) = RCT(227)*V(65)

  B(394) = RCT(228)*V(136)

  B(395) = RCT(228)*V(79)

  B(396) = RCT(229)*V(136)

  B(397) = RCT(229)*V(110)

  B(398) = RCT(230)*V(135)

  B(399) = RCT(230)*V(110)

  B(400) = RCT(231)*V(133)

  B(401) = RCT(231)*V(120)

  B(402) = RCT(232)*V(131)

  B(403) = RCT(232)*V(120)

  B(404) = RCT(233)*V(129)

  B(405) = RCT(233)*V(120)

  B(406) = RCT(234)*V(136)

  B(407) = RCT(234)*V(114)

  B(408) = RCT(235)*V(133)

  B(409) = RCT(235)*V(119)

  B(410) = RCT(236)*V(131)

  B(411) = RCT(236)*V(119)

  B(412) = RCT(237)*V(129)

  B(413) = RCT(237)*V(119)

  B(414) = RCT(238)*V(135)

  B(415) = RCT(238)*V(119)

  B(416) = RCT(239)*V(136)

  B(417) = RCT(239)*V(63)

  B(418) = RCT(240)*V(135)

  B(419) = RCT(240)*V(106)

  B(420) = RCT(241)*V(133)

  B(421) = RCT(241)*V(99)

  B(422) = RCT(242)*V(135)

  B(423) = RCT(242)*V(99)

  B(424) = RCT(243)*V(131)

  B(425) = RCT(243)*V(99)

  B(426) = RCT(244)*V(136)

  B(427) = RCT(244)*V(117)

  B(428) = RCT(245)*V(135)

  B(429) = RCT(245)*V(117)

  B(430) = RCT(246)*V(136)

  B(431) = RCT(246)*V(116)

  B(432) = RCT(247)*V(135)

  B(433) = RCT(247)*V(116)

  B(434) = RCT(248)*V(136)

  B(435) = RCT(248)*V(53)

  B(436) = RCT(249)*V(133)

  B(437) = RCT(249)*V(123)

  B(438) = RCT(250)*V(135)

  B(439) = RCT(250)*V(123)

  B(440) = RCT(251)*V(131)

  B(441) = RCT(251)*V(123)

  B(442) = RCT(252)*V(129)

  B(443) = RCT(252)*V(123)

  B(444) = RCT(253)*V(128)

  B(445) = RCT(253)*V(123)

  B(446) = RCT(254)*V(136)

  B(447) = RCT(254)*V(26)

  B(448) = RCT(255)*V(136)

  B(449) = RCT(255)*V(26)

  B(450) = RCT(256)*V(136)

  B(451) = RCT(256)*V(95)

  B(452) = RCT(257)*V(136)

  B(453) = RCT(257)*V(41)

  B(454) = RCT(258)*V(136)

  B(455) = RCT(258)*V(77)

  B(456) = RCT(259)*V(136)

  B(457) = RCT(259)*V(68)

  B(458) = RCT(260)*V(136)

  B(459) = RCT(260)*V(122)

  B(460) = RCT(261)*V(136)

  B(461) = RCT(261)*V(108)

  B(462) = RCT(262)*V(136)

  B(463) = RCT(262)*V(39)

  B(464) = RCT(263)*V(136)

  B(465) = RCT(263)*V(39)

  B(466) = RCT(264)*V(135)

  B(467) = RCT(264)*V(39)

  B(468) = RCT(265)*V(136)

  B(469) = RCT(265)*V(11)

  B(470) = RCT(266)

  B(471) = RCT(267)*2*V(97)

  B(472) = RCT(268)*V(136)

  B(473) = RCT(268)*V(89)

  B(474) = RCT(269)*V(133)

  B(475) = RCT(269)*V(111)

  B(476) = RCT(270)*V(129)

  B(477) = RCT(270)*V(111)

  B(478) = RCT(271)*V(136)

  B(479) = RCT(271)*V(103)

  B(480) = RCT(272)*V(133)

  B(481) = RCT(272)*V(73)

  B(482) = RCT(273)*V(131)

  B(483) = RCT(273)*V(73)

  B(484) = RCT(274)*V(131)

  B(485) = RCT(274)*V(111)

  B(486) = RCT(275)*V(136)

  B(487) = RCT(275)*V(76)

  B(488) = RCT(276)*V(132)

  B(489) = RCT(276)*V(89)

  B(490) = RCT(277)*V(135)

  B(491) = RCT(277)*V(89)

  B(492) = RCT(278)*V(131)

  B(493) = RCT(278)*V(102)

  B(494) = RCT(279)*V(133)

  B(495) = RCT(279)*V(102)

  B(496) = RCT(280)*V(135)

  B(497) = RCT(280)*V(102)

  B(498) = RCT(281)*V(136)

  B(499) = RCT(281)*V(19)

  B(500) = 4.5e-13*V(136)

  B(501) = 4.5e-13*V(78)

  B(502) = RCT(283)*V(131)

  B(503) = RCT(283)*V(126)

  B(504) = RCT(284)

  B(505) = RCT(285)*V(133)

  B(506) = RCT(285)*V(56)

  B(507) = RCT(286)*V(131)

  B(508) = RCT(286)*V(56)

  B(509) = RCT(287)*V(136)

  B(510) = RCT(287)*V(8)

  B(511) = RCT(288)*V(136)

  B(512) = RCT(288)*V(10)

  B(513) = RCT(289)*V(136)

  B(514) = RCT(289)*V(106)

  B(515) = RCT(290)*V(136)

  B(516) = RCT(290)*V(82)

  B(517) = RCT(291)*V(136)

  B(518) = RCT(291)*V(83)

  B(519) = RCT(292)*V(136)

  B(520) = RCT(292)*V(101)



  JVS(1) = -B(49)

  JVS(2) = 0

  JVS(3) = B(132)

  JVS(4) = B(133)

  JVS(5) = 0

  JVS(6) = B(201)

  JVS(7) = B(505)

  JVS(8) = B(283)

  JVS(9) = B(136)

  JVS(10) = B(319)

  JVS(11) = B(299)

  JVS(12) = B(348)

  JVS(13) = B(480)

  JVS(14) = B(326)

  JVS(15) = B(293)

  JVS(16) = B(215)

  JVS(17) = B(157)

  JVS(18) = B(307)

  JVS(19) = B(330)

  JVS(20) = B(358)

  JVS(21) = B(334)

  JVS(22) = B(181)

  JVS(23) = B(338)

  JVS(24) = B(143)

  JVS(25) = B(420)

  JVS(26) = B(494)

  JVS(27) = B(474)

  JVS(28) = B(205)

  JVS(29) = B(386)

  JVS(30) = B(408)

  JVS(31) = B(400)

  JVS(32) = B(191)

  JVS(33) = B(436)

  JVS(34) = B(249)+B(251)

  JVS(35) = B(229)

  JVS(36) = B(167)

  JVS(37) = B(114)

  JVS(38) = B(263)

  JVS(39) = B(115)+B(137)+B(144)+B(158)+B(168)+B(182)+B(192)+B(202)+B(206)+B(216)+B(230)+B(250)+B(252)+B(264)+B(284)&
              &+B(294)+B(300)+B(308)+B(320)+B(327)+B(331)+B(335)+B(339)+B(349)+B(359)+B(387)+B(401)+B(409)+B(421)+B(437)&
              &+B(475)+B(481)+B(495)+B(506)

  JVS(40) = 0

  JVS(41) = B(507)

  JVS(42) = B(285)

  JVS(43) = B(315)

  JVS(44) = B(301)

  JVS(45) = B(350)

  JVS(46) = B(482)

  JVS(47) = B(328)

  JVS(48) = B(295)

  JVS(49) = B(217)

  JVS(50) = B(159)

  JVS(51) = B(309)

  JVS(52) = B(221)

  JVS(53) = B(354)

  JVS(54) = B(332)

  JVS(55) = B(183)+B(185)

  JVS(56) = B(336)

  JVS(57) = B(145)+B(147)+B(471)

  JVS(58) = B(424)

  JVS(59) = B(492)

  JVS(60) = B(476)+B(484)

  JVS(61) = B(207)

  JVS(62) = B(388)+B(390)

  JVS(63) = B(410)+B(412)

  JVS(64) = B(402)+B(404)

  JVS(65) = B(193)+B(195)

  JVS(66) = B(440)+B(442)+B(444)

  JVS(67) = B(255)+B(257)+B(259)

  JVS(68) = B(233)+B(237)+B(239)

  JVS(69) = B(171)+B(173)+B(178)+B(240)+B(260)+B(271)+B(445)

  JVS(70) = B(116)+B(117)+B(118)+B(148)+B(174)+B(186)+B(196)+B(238)+B(258)+B(269)+B(391)+B(405)+B(413)+B(443)+B(477)

  JVS(71) = B(267)+B(270)+B(272)+B(273)

  JVS(72) = B(119)+B(146)+B(160)+B(172)+B(184)+B(194)+B(208)+B(218)+B(222)+B(234)+B(256)+B(268)+B(286)+B(296)+B(302)&
              &+B(310)+B(316)+B(329)+B(333)+B(337)+B(351)+B(355)+B(389)+B(403)+B(411)+B(425)+B(441)+B(483)+B(485)+B(493)&
              &+B(508)

  JVS(73) = 0

  JVS(74) = B(513)

  JVS(75) = B(514)

  JVS(76) = 0

  JVS(77) = B(515)

  JVS(78) = B(517)

  JVS(79) = B(519)

  JVS(80) = B(516)+B(518)+B(520)

  JVS(81) = 0

  JVS(82) = B(509)

  JVS(83) = B(510)

  JVS(84) = -B(509)

  JVS(85) = -B(510)

  JVS(86) = 0

  JVS(87) = B(511)

  JVS(88) = B(512)

  JVS(89) = -B(511)

  JVS(90) = -B(512)

  JVS(91) = -B(468)

  JVS(92) = -B(469)

  JVS(93) = -B(141)

  JVS(94) = -B(142)

  JVS(95) = -B(179)

  JVS(96) = -B(180)

  JVS(97) = -B(199)

  JVS(98) = -B(200)

  JVS(99) = -B(4)-B(78)-B(80)

  JVS(100) = -B(79)-B(81)

  JVS(101) = -B(132)

  JVS(102) = B(462)+0.5*B(464)+B(466)

  JVS(103) = B(467)

  JVS(104) = -B(133)+B(463)+0.5*B(465)

  JVS(105) = -B(325)

  JVS(106) = B(323)

  JVS(107) = B(324)

  JVS(108) = -B(13)-B(72)

  JVS(109) = B(70)+0.5*B(470)

  JVS(110) = -B(73)+B(77)

  JVS(111) = -B(498)

  JVS(112) = -B(499)

  JVS(113) = -B(279)

  JVS(114) = -B(280)

  JVS(115) = 0.12*B(279)

  JVS(116) = -B(38)

  JVS(117) = 0.12*B(280)

  JVS(118) = 0.53*B(279)

  JVS(119) = -B(281)

  JVS(120) = 0.53*B(280)-B(282)

  JVS(121) = -B(138)-B(140)

  JVS(122) = B(136)

  JVS(123) = B(137)

  JVS(124) = -B(311)

  JVS(125) = -B(312)

  JVS(126) = 0.18*B(311)

  JVS(127) = -B(313)

  JVS(128) = 0.18*B(312)-B(314)

  JVS(129) = -B(26)-B(446)-B(448)

  JVS(130) = B(440)

  JVS(131) = B(441)

  JVS(132) = -B(447)-B(449)

  JVS(133) = -B(287)

  JVS(134) = B(285)

  JVS(135) = B(286)

  JVS(136) = -B(288)

  JVS(137) = -B(297)

  JVS(138) = B(295)

  JVS(139) = B(296)

  JVS(140) = -B(298)

  JVS(141) = -B(303)

  JVS(142) = B(301)

  JVS(143) = B(302)

  JVS(144) = -B(304)

  JVS(145) = -B(40)

  JVS(146) = 0.2*B(330)

  JVS(147) = 0.2*B(358)

  JVS(148) = 0.2*B(331)+0.2*B(359)

  JVS(149) = 0.1*B(311)

  JVS(150) = -B(39)

  JVS(151) = 0.23*B(344)

  JVS(152) = 0.1*B(312)+0.23*B(345)

  JVS(153) = -B(317)

  JVS(154) = B(315)

  JVS(155) = B(316)

  JVS(156) = -B(318)

  JVS(157) = -B(321)

  JVS(158) = B(319)

  JVS(159) = B(320)

  JVS(160) = -B(322)

  JVS(161) = -B(352)

  JVS(162) = B(350)

  JVS(163) = B(351)

  JVS(164) = -B(353)

  JVS(165) = -B(356)

  JVS(166) = B(354)

  JVS(167) = B(355)

  JVS(168) = -B(357)

  JVS(169) = -B(6)-B(94)-B(106)

  JVS(170) = B(92)

  JVS(171) = B(93)

  JVS(172) = -B(344)

  JVS(173) = -B(345)

  JVS(174) = 0.15*B(344)

  JVS(175) = -B(346)

  JVS(176) = 0.15*B(345)-B(347)

  JVS(177) = -B(462)-B(464)-B(466)

  JVS(178) = -B(467)

  JVS(179) = -B(463)-B(465)

  JVS(180) = -B(42)

  JVS(181) = 0.21*B(358)

  JVS(182) = 0.21*B(359)

  JVS(183) = -B(452)

  JVS(184) = 0.2*B(147)+0.4*B(471)

  JVS(185) = 0.2*B(148)

  JVS(186) = -B(453)

  JVS(187) = -B(41)

  JVS(188) = 0.2*B(330)

  JVS(189) = 0.15*B(358)

  JVS(190) = 0.2*B(331)+0.15*B(359)

  JVS(191) = B(38)

  JVS(192) = -B(37)

  JVS(193) = 0.5*B(299)

  JVS(194) = 0.2*B(330)

  JVS(195) = 0.06*B(358)

  JVS(196) = 0.5*B(300)+0.2*B(331)+0.06*B(359)

  JVS(197) = 0

  JVS(198) = -B(58)-B(60)

  JVS(199) = 0.05*B(112)

  JVS(200) = -B(59)+0.05*B(113)

  JVS(201) = B(12)

  JVS(202) = -B(61)

  JVS(203) = -B(21)-B(149)

  JVS(204) = B(145)

  JVS(205) = B(146)

  JVS(206) = -B(150)

  JVS(207) = -B(22)-B(187)

  JVS(208) = B(183)

  JVS(209) = B(184)

  JVS(210) = -B(188)

  JVS(211) = -B(23)-B(197)

  JVS(212) = B(193)

  JVS(213) = B(194)

  JVS(214) = -B(198)

  JVS(215) = B(199)

  JVS(216) = -B(201)

  JVS(217) = -B(202)

  JVS(218) = B(200)

  JVS(219) = -B(35)-B(219)

  JVS(220) = B(217)

  JVS(221) = B(218)

  JVS(222) = -B(220)

  JVS(223) = -B(10)-B(120)

  JVS(224) = B(118)

  JVS(225) = B(119)

  JVS(226) = -B(121)

  JVS(227) = -B(261)

  JVS(228) = B(255)

  JVS(229) = B(256)

  JVS(230) = -B(262)

  JVS(231) = -B(9)-B(103)-B(105)

  JVS(232) = B(101)

  JVS(233) = B(102)

  JVS(234) = -B(104)

  JVS(235) = -B(434)

  JVS(236) = B(430)+B(432)

  JVS(237) = 0.33*B(229)+0.36*B(231)+0.3*B(237)+0.4*B(239)

  JVS(238) = 0.4*B(240)

  JVS(239) = 0.3*B(238)

  JVS(240) = 0.33*B(230)

  JVS(241) = 0.36*B(232)+B(433)

  JVS(242) = B(431)-B(435)

  JVS(243) = -B(203)

  JVS(244) = 0.05*B(370)

  JVS(245) = 0.05*B(372)

  JVS(246) = 0.05*B(378)

  JVS(247) = 0.05*B(374)

  JVS(248) = 0.05*B(376)

  JVS(249) = 0.05*B(371)+0.05*B(373)+0.05*B(375)+0.05*B(377)+0.05*B(379)

  JVS(250) = -B(204)

  JVS(251) = -B(128)-B(130)

  JVS(252) = -B(131)

  JVS(253) = -B(129)

  JVS(254) = -B(504)-B(505)-B(507)

  JVS(255) = B(502)

  JVS(256) = B(503)-B(508)

  JVS(257) = -B(506)

  JVS(258) = -B(36)-B(223)

  JVS(259) = B(221)

  JVS(260) = B(222)

  JVS(261) = -B(224)

  JVS(262) = 0.14*B(281)

  JVS(263) = 0.2*B(313)

  JVS(264) = B(287)

  JVS(265) = -B(283)-B(285)

  JVS(266) = -B(286)

  JVS(267) = -B(284)

  JVS(268) = 0.14*B(282)+B(288)+0.2*B(314)

  JVS(269) = 0.06*B(281)

  JVS(270) = 0.07*B(313)

  JVS(271) = 0.07*B(346)

  JVS(272) = -B(289)-B(291)

  JVS(273) = B(293)

  JVS(274) = -B(292)

  JVS(275) = B(294)

  JVS(276) = -B(290)

  JVS(277) = 0.06*B(282)+0.07*B(314)+0.07*B(347)

  JVS(278) = 0.75*B(128)

  JVS(279) = -B(136)

  JVS(280) = 0

  JVS(281) = -B(137)

  JVS(282) = 0.75*B(129)

  JVS(283) = 0.07*B(311)

  JVS(284) = B(317)

  JVS(285) = 0.06*B(344)

  JVS(286) = -B(315)-B(319)

  JVS(287) = -B(316)

  JVS(288) = -B(320)

  JVS(289) = 0.07*B(312)+B(318)+0.06*B(345)

  JVS(290) = -B(31)-B(213)

  JVS(291) = 0.8*B(34)

  JVS(292) = 0.75*B(205)

  JVS(293) = 0.75*B(206)

  JVS(294) = -B(214)

  JVS(295) = -B(416)

  JVS(296) = 0.2*B(227)

  JVS(297) = 0.25*B(153)

  JVS(298) = 0.25*B(171)+0.1*B(173)

  JVS(299) = 0.1*B(174)

  JVS(300) = 0.25*B(267)

  JVS(301) = 0.25*B(172)+0.25*B(268)

  JVS(302) = 0.25*B(154)+0.2*B(228)

  JVS(303) = -B(417)

  JVS(304) = -B(16)-B(175)

  JVS(305) = 0.75*B(171)

  JVS(306) = 0.75*B(267)

  JVS(307) = 0.75*B(172)+0.75*B(268)

  JVS(308) = -B(176)

  JVS(309) = -B(47)-B(392)

  JVS(310) = B(388)

  JVS(311) = B(389)

  JVS(312) = -B(393)

  JVS(313) = 0.35*B(279)

  JVS(314) = B(303)

  JVS(315) = -B(299)-B(301)

  JVS(316) = -B(302)

  JVS(317) = -B(300)

  JVS(318) = 0.35*B(280)+B(304)

  JVS(319) = -B(211)

  JVS(320) = 0.75*B(494)+0.75*B(496)

  JVS(321) = B(155)

  JVS(322) = 0.1*B(205)

  JVS(323) = 0.1*B(206)+0.75*B(495)

  JVS(324) = B(156)+0.75*B(497)

  JVS(325) = -B(212)

  JVS(326) = -B(17)-B(177)-B(456)

  JVS(327) = B(169)

  JVS(328) = B(170)

  JVS(329) = -B(457)

  JVS(330) = B(352)

  JVS(331) = 0.3*B(346)

  JVS(332) = -B(348)-B(350)

  JVS(333) = -B(351)

  JVS(334) = -B(349)

  JVS(335) = 0.3*B(347)+B(353)

  JVS(336) = -B(15)-B(161)

  JVS(337) = B(159)

  JVS(338) = B(160)

  JVS(339) = -B(162)

  JVS(340) = 2*B(106)

  JVS(341) = B(466)

  JVS(342) = -B(7)-B(97)

  JVS(343) = B(165)

  JVS(344) = B(428)

  JVS(345) = B(122)

  JVS(346) = B(95)+0.5*B(109)

  JVS(347) = B(108)+B(123)+B(166)+B(429)+B(467)

  JVS(348) = B(96)-B(98)

  JVS(349) = -B(110)-B(112)

  JVS(350) = -B(113)

  JVS(351) = 0.08*B(153)

  JVS(352) = 0.08*B(154)

  JVS(353) = -B(111)

  JVS(354) = -B(480)-B(482)

  JVS(355) = B(478)

  JVS(356) = -B(483)

  JVS(357) = -B(481)

  JVS(358) = B(479)

  JVS(359) = B(325)

  JVS(360) = B(321)

  JVS(361) = 0

  JVS(362) = -B(323)-B(326)-B(328)

  JVS(363) = -B(329)

  JVS(364) = -B(327)

  JVS(365) = -B(324)

  JVS(366) = B(322)

  JVS(367) = -B(28)-B(235)

  JVS(368) = B(233)

  JVS(369) = B(234)

  JVS(370) = -B(236)

  JVS(371) = -B(43)-B(486)

  JVS(372) = B(484)

  JVS(373) = B(485)

  JVS(374) = -B(487)

  JVS(375) = -B(18)-B(277)-B(454)

  JVS(376) = B(274)

  JVS(377) = B(275)

  JVS(378) = -B(455)

  JVS(379) = 0.35*B(498)

  JVS(380) = 0.5*B(130)

  JVS(381) = B(505)+B(507)

  JVS(382) = -B(500)

  JVS(383) = 0.05*B(370)

  JVS(384) = 0.05*B(372)

  JVS(385) = 0.25*B(488)

  JVS(386) = 0.05*B(378)

  JVS(387) = 0.05*B(374)

  JVS(388) = 0.05*B(376)

  JVS(389) = 0

  JVS(390) = B(508)

  JVS(391) = 0.5*B(131)+0.05*B(371)+0.05*B(373)+0.05*B(375)+0.05*B(377)+0.05*B(379)+0.25*B(489)

  JVS(392) = B(506)

  JVS(393) = 0.35*B(499)-B(501)

  JVS(394) = -B(48)-B(394)

  JVS(395) = B(402)

  JVS(396) = B(403)

  JVS(397) = -B(395)

  JVS(398) = B(4)-B(78)-B(80)

  JVS(399) = -B(58)

  JVS(400) = -B(112)

  JVS(401) = -B(54)-B(56)-B(59)-B(79)-B(81)-B(113)

  JVS(402) = 0

  JVS(403) = 0

  JVS(404) = B(2)

  JVS(405) = 0

  JVS(406) = B(297)

  JVS(407) = B(291)

  JVS(408) = B(326)+0.4*B(328)

  JVS(409) = -B(293)-B(295)

  JVS(410) = -B(296)+0.4*B(329)

  JVS(411) = B(292)

  JVS(412) = -B(294)+B(327)

  JVS(413) = 0

  JVS(414) = B(298)

  JVS(415) = -B(360)-B(370)

  JVS(416) = -B(371)

  JVS(417) = -B(361)

  JVS(418) = -B(362)-B(372)

  JVS(419) = -B(373)

  JVS(420) = -B(363)

  JVS(421) = -B(34)-B(209)

  JVS(422) = B(207)

  JVS(423) = B(208)

  JVS(424) = -B(210)

  JVS(425) = B(219)

  JVS(426) = B(213)

  JVS(427) = 0

  JVS(428) = -B(215)-B(217)

  JVS(429) = 0

  JVS(430) = -B(218)

  JVS(431) = -B(216)

  JVS(432) = B(214)+B(220)

  JVS(433) = 0.5*B(161)

  JVS(434) = -B(157)-B(159)

  JVS(435) = B(151)

  JVS(436) = -B(160)

  JVS(437) = -B(158)

  JVS(438) = B(152)+0.5*B(162)

  JVS(439) = 0.6*B(37)

  JVS(440) = 0

  JVS(441) = -B(305)-B(307)-B(309)

  JVS(442) = 0

  JVS(443) = 0

  JVS(444) = -B(310)

  JVS(445) = -B(308)

  JVS(446) = -B(306)

  JVS(447) = 0

  JVS(448) = 0.65*B(311)

  JVS(449) = B(223)

  JVS(450) = -B(221)-B(330)

  JVS(451) = -B(222)

  JVS(452) = -B(331)

  JVS(453) = B(224)+0.65*B(312)

  JVS(454) = -B(472)-B(488)-B(490)

  JVS(455) = -B(489)

  JVS(456) = -B(491)

  JVS(457) = -B(473)

  JVS(458) = B(356)

  JVS(459) = 0.56*B(344)

  JVS(460) = -B(354)-B(358)

  JVS(461) = -B(355)

  JVS(462) = -B(359)

  JVS(463) = 0.56*B(345)+B(357)

  JVS(464) = 0.6*B(40)

  JVS(465) = 0

  JVS(466) = 0

  JVS(467) = -B(332)-B(334)-B(340)

  JVS(468) = -B(333)

  JVS(469) = -B(335)

  JVS(470) = -B(341)

  JVS(471) = 0

  JVS(472) = B(54)

  JVS(473) = -B(50)-B(52)-B(62)-B(64)-B(86)

  JVS(474) = 0

  JVS(475) = 0

  JVS(476) = -B(65)

  JVS(477) = B(3)-B(53)

  JVS(478) = B(5)-B(87)

  JVS(479) = -B(63)+B(76)

  JVS(480) = B(179)

  JVS(481) = B(187)

  JVS(482) = -B(181)-B(183)-B(185)

  JVS(483) = -B(186)

  JVS(484) = -B(184)

  JVS(485) = -B(182)

  JVS(486) = B(180)+B(188)

  JVS(487) = 0.6*B(41)

  JVS(488) = 0

  JVS(489) = 0

  JVS(490) = -B(336)-B(338)-B(342)

  JVS(491) = -B(337)

  JVS(492) = -B(339)

  JVS(493) = -B(343)

  JVS(494) = 0

  JVS(495) = -B(450)

  JVS(496) = 0.3*B(147)

  JVS(497) = 0.25*B(476)

  JVS(498) = 0.25*B(390)

  JVS(499) = 0.25*B(412)

  JVS(500) = 0.25*B(404)

  JVS(501) = 0.5*B(195)

  JVS(502) = 0.3*B(442)

  JVS(503) = 0.25*B(257)

  JVS(504) = 0.25*B(237)

  JVS(505) = B(117)+0.3*B(148)+0.5*B(196)+0.25*B(238)+0.25*B(258)+0.25*B(391)+0.25*B(405)+0.25*B(413)+0.3*B(443)+0.25&
               &*B(477)

  JVS(506) = -B(451)

  JVS(507) = 0.65*B(498)

  JVS(508) = 0.45*B(36)

  JVS(509) = 0.7*B(283)

  JVS(510) = B(299)

  JVS(511) = 0.17*B(348)

  JVS(512) = 0.4*B(307)+0.16*B(309)

  JVS(513) = 0.6*B(330)

  JVS(514) = 0.34*B(358)

  JVS(515) = 0.07*B(336)+0.17*B(338)

  JVS(516) = -B(33)-B(134)

  JVS(517) = 0.13*B(32)

  JVS(518) = 0.2*B(460)

  JVS(519) = 0.25*B(436)+0.25*B(438)+0.1*B(442)+0.25*B(444)

  JVS(520) = 0.02*B(229)+0.02*B(231)

  JVS(521) = 0.25*B(445)

  JVS(522) = 0.1*B(443)

  JVS(523) = 0.16*B(310)+0.07*B(337)

  JVS(524) = 0.02*B(230)+0.7*B(284)+B(300)+0.4*B(308)+0.6*B(331)+0.17*B(339)+0.17*B(349)+0.34*B(359)+0.25*B(437)

  JVS(525) = 0

  JVS(526) = 0.02*B(232)+0.25*B(439)

  JVS(527) = -B(135)+0.2*B(461)+0.65*B(499)

  JVS(528) = B(141)

  JVS(529) = 0.5*B(149)

  JVS(530) = B(31)

  JVS(531) = 0

  JVS(532) = -B(143)-B(145)-B(147)-2*B(471)

  JVS(533) = 0

  JVS(534) = -B(148)

  JVS(535) = -B(146)

  JVS(536) = -B(144)

  JVS(537) = B(142)+0.5*B(150)

  JVS(538) = 0.9*B(36)

  JVS(539) = 0.1*B(370)

  JVS(540) = 0.1*B(372)

  JVS(541) = 0

  JVS(542) = -B(32)

  JVS(543) = 0.1*B(378)

  JVS(544) = 0.1*B(374)

  JVS(545) = 0.1*B(376)

  JVS(546) = 0

  JVS(547) = 0.1*B(371)+0.1*B(373)+0.1*B(375)+0.1*B(377)+0.1*B(379)

  JVS(548) = 0

  JVS(549) = 0

  JVS(550) = -B(420)-B(422)-B(424)

  JVS(551) = B(418)

  JVS(552) = -B(425)

  JVS(553) = -B(421)

  JVS(554) = B(419)-B(423)

  JVS(555) = -B(368)-B(378)-B(384)

  JVS(556) = -B(379)

  JVS(557) = -B(385)

  JVS(558) = -B(369)

  JVS(559) = -B(364)-B(374)-B(380)

  JVS(560) = -B(375)

  JVS(561) = -B(381)

  JVS(562) = -B(365)

  JVS(563) = B(490)

  JVS(564) = -B(492)-B(494)-B(496)

  JVS(565) = -B(493)

  JVS(566) = 0

  JVS(567) = -B(495)

  JVS(568) = B(491)-B(497)

  JVS(569) = 0

  JVS(570) = 0.33*B(43)

  JVS(571) = 0.9*B(488)

  JVS(572) = 0.25*B(494)+0.25*B(496)

  JVS(573) = -B(44)-B(478)

  JVS(574) = 0.33*B(474)+0.167*B(476)

  JVS(575) = 0.167*B(477)

  JVS(576) = 0

  JVS(577) = 0.9*B(489)

  JVS(578) = 0.33*B(475)+0.25*B(495)

  JVS(579) = 0.25*B(497)

  JVS(580) = -B(479)

  JVS(581) = -B(366)-B(376)-B(382)

  JVS(582) = -B(377)

  JVS(583) = -B(383)

  JVS(584) = -B(367)

  JVS(585) = 0.82*B(22)

  JVS(586) = 0.5*B(201)

  JVS(587) = 0.05*B(47)

  JVS(588) = B(480)+0.4*B(482)

  JVS(589) = 0.67*B(43)

  JVS(590) = 0.3*B(48)

  JVS(591) = 0.52*B(370)

  JVS(592) = 0.52*B(372)

  JVS(593) = 0.25*B(34)

  JVS(594) = 0.1*B(488)

  JVS(595) = 0.82*B(181)+0.82*B(185)

  JVS(596) = 0.52*B(378)

  JVS(597) = 0.52*B(374)

  JVS(598) = 0.75*B(494)+0.75*B(496)

  JVS(599) = B(44)

  JVS(600) = 0.52*B(376)

  JVS(601) = -B(24)-B(189)

  JVS(602) = 0.67*B(474)+0.333*B(476)

  JVS(603) = 0.25*B(205)

  JVS(604) = 0.5*B(46)+0.5*B(406)

  JVS(605) = 0.045*B(386)+0.025*B(390)

  JVS(606) = 0.27*B(400)+0.15*B(404)

  JVS(607) = 0.82*B(186)+0.025*B(391)+0.15*B(405)+0.333*B(477)

  JVS(608) = 0.4*B(483)

  JVS(609) = 0.52*B(371)+0.52*B(373)+0.52*B(375)+0.52*B(377)+0.52*B(379)+0.1*B(489)

  JVS(610) = 0.82*B(182)+0.5*B(202)+0.25*B(206)+0.045*B(387)+0.27*B(401)+0.67*B(475)+B(481)+0.75*B(495)

  JVS(611) = 0.75*B(497)

  JVS(612) = -B(190)+0.5*B(407)

  JVS(613) = -B(225)-B(227)-B(418)

  JVS(614) = -B(228)

  JVS(615) = -B(419)

  JVS(616) = -B(226)

  JVS(617) = 0.35*B(498)

  JVS(618) = 1.5*B(38)

  JVS(619) = 1.5*B(39)

  JVS(620) = 0

  JVS(621) = B(42)

  JVS(622) = 0.6*B(41)

  JVS(623) = 0.5*B(130)

  JVS(624) = 0.25*B(48)

  JVS(625) = 0.23*B(370)

  JVS(626) = 0.23*B(372)

  JVS(627) = 0.4*B(307)+0.16*B(309)

  JVS(628) = 0

  JVS(629) = 0.35*B(488)

  JVS(630) = 0

  JVS(631) = 0.07*B(332)+0.17*B(334)

  JVS(632) = 0.14*B(336)+0.35*B(338)

  JVS(633) = 2*B(33)+B(134)

  JVS(634) = 0.45*B(32)

  JVS(635) = 0.23*B(378)

  JVS(636) = 0.23*B(374)

  JVS(637) = B(44)

  JVS(638) = 0.23*B(376)

  JVS(639) = 0.3*B(227)

  JVS(640) = -B(126)

  JVS(641) = B(30)

  JVS(642) = 0.56*B(153)

  JVS(643) = B(45)

  JVS(644) = 0

  JVS(645) = B(14)

  JVS(646) = 1.7*B(46)+0.7*B(406)

  JVS(647) = 0.67*B(19)+0.2*B(247)

  JVS(648) = B(27)

  JVS(649) = B(25)+B(426)+B(428)

  JVS(650) = 0.225*B(400)+0.125*B(404)

  JVS(651) = 0.25*B(436)+0.25*B(438)+0.2*B(442)+0.25*B(444)

  JVS(652) = 0.7*B(20)+0.05*B(243)

  JVS(653) = 0.22*B(249)+0.22*B(253)+0.11*B(257)+0.22*B(259)

  JVS(654) = B(11)+B(12)+B(122)+B(124)

  JVS(655) = 0

  JVS(656) = 0.22*B(260)+0.25*B(445)

  JVS(657) = 0.11*B(258)+0.125*B(405)+0.2*B(443)

  JVS(658) = 0.16*B(310)+0.07*B(333)+0.14*B(337)

  JVS(659) = 0.5*B(131)+0.56*B(154)+0.3*B(228)+0.05*B(244)+0.2*B(248)+0.23*B(371)+0.23*B(373)+0.23*B(375)+0.23*B(377)&
               &+0.23*B(379)+0.35*B(489)

  JVS(660) = 0.22*B(250)+0.4*B(308)+0.17*B(335)+0.35*B(339)+0.225*B(401)+0.25*B(437)

  JVS(661) = 0

  JVS(662) = B(123)+0.22*B(254)+B(429)+0.25*B(439)

  JVS(663) = B(125)-B(127)+B(135)+0.7*B(407)+B(427)+0.35*B(499)

  JVS(664) = B(138)

  JVS(665) = 0

  JVS(666) = 0.67*B(43)

  JVS(667) = 0.25*B(48)

  JVS(668) = -B(30)-B(460)

  JVS(669) = 0.67*B(474)+0.333*B(476)

  JVS(670) = 0.225*B(400)+0.125*B(404)

  JVS(671) = 0.25*B(436)+0.25*B(438)+0.1*B(442)+0.25*B(444)

  JVS(672) = 0.53*B(249)+0.53*B(253)+0.26*B(257)+0.53*B(259)

  JVS(673) = 0.02*B(229)+0.02*B(231)

  JVS(674) = 0.53*B(260)+0.25*B(445)

  JVS(675) = 0.26*B(258)+0.125*B(405)+0.1*B(443)+0.333*B(477)

  JVS(676) = 0

  JVS(677) = 0

  JVS(678) = 0.02*B(230)+0.53*B(250)+0.225*B(401)+0.25*B(437)+0.67*B(475)

  JVS(679) = 0.02*B(232)+0.53*B(254)+0.25*B(439)

  JVS(680) = -B(461)

  JVS(681) = 0.07*B(227)

  JVS(682) = -B(151)-B(153)-B(155)

  JVS(683) = 0.7*B(20)

  JVS(684) = -B(154)+0.07*B(228)

  JVS(685) = -B(156)

  JVS(686) = -B(152)

  JVS(687) = B(47)

  JVS(688) = 0.33*B(370)

  JVS(689) = 0.33*B(372)

  JVS(690) = 0.33*B(378)

  JVS(691) = 0.33*B(374)

  JVS(692) = 0.33*B(376)

  JVS(693) = -B(45)-B(396)-B(398)

  JVS(694) = 0.9*B(386)+B(390)

  JVS(695) = 0.9*B(408)+0.5*B(412)+B(414)

  JVS(696) = B(391)+0.5*B(413)

  JVS(697) = 0

  JVS(698) = 0.33*B(371)+0.33*B(373)+0.33*B(375)+0.33*B(377)+0.33*B(379)

  JVS(699) = 0.9*B(387)+0.9*B(409)

  JVS(700) = -B(399)+B(415)

  JVS(701) = -B(397)

  JVS(702) = 0.5*B(486)

  JVS(703) = B(472)

  JVS(704) = -B(474)-B(476)-B(484)

  JVS(705) = -B(477)

  JVS(706) = -B(485)

  JVS(707) = 0

  JVS(708) = -B(475)

  JVS(709) = 0

  JVS(710) = B(473)+0.5*B(487)

  JVS(711) = B(452)

  JVS(712) = B(21)+0.5*B(149)

  JVS(713) = B(201)

  JVS(714) = B(35)

  JVS(715) = B(15)

  JVS(716) = 0.4*B(34)

  JVS(717) = B(215)

  JVS(718) = B(157)

  JVS(719) = 0.27*B(181)

  JVS(720) = B(143)+0.8*B(147)+1.6*B(471)

  JVS(721) = 0.5*B(153)

  JVS(722) = -B(14)-B(163)-B(165)

  JVS(723) = 0.4*B(205)

  JVS(724) = 0.04*B(243)

  JVS(725) = 0.8*B(148)

  JVS(726) = 0

  JVS(727) = 0.5*B(154)+0.04*B(244)

  JVS(728) = B(144)+B(158)+0.27*B(182)+B(202)+0.4*B(206)+B(216)

  JVS(729) = -B(166)

  JVS(730) = 0.5*B(150)-B(164)+B(453)

  JVS(731) = B(203)

  JVS(732) = 0

  JVS(733) = 0

  JVS(734) = B(209)

  JVS(735) = 0

  JVS(736) = 0

  JVS(737) = 0

  JVS(738) = -B(205)-B(207)

  JVS(739) = -B(208)

  JVS(740) = 0

  JVS(741) = -B(206)

  JVS(742) = 0

  JVS(743) = B(204)+B(210)

  JVS(744) = B(48)

  JVS(745) = 0.3*B(370)

  JVS(746) = 0.3*B(372)

  JVS(747) = 0.3*B(378)

  JVS(748) = 0.3*B(374)

  JVS(749) = 0.3*B(376)

  JVS(750) = B(45)

  JVS(751) = -B(46)-B(406)

  JVS(752) = 0

  JVS(753) = 0

  JVS(754) = 0.9*B(400)+B(404)

  JVS(755) = B(405)

  JVS(756) = 0

  JVS(757) = 0.3*B(371)+0.3*B(373)+0.3*B(375)+0.3*B(377)+0.3*B(379)

  JVS(758) = 0.9*B(401)

  JVS(759) = 0

  JVS(760) = -B(407)

  JVS(761) = 0.288*B(28)

  JVS(762) = 0.167*B(420)+0.167*B(422)+0.167*B(424)

  JVS(763) = 0.4*B(227)

  JVS(764) = -B(19)-B(245)-B(247)

  JVS(765) = 0.23*B(229)+0.25*B(231)+0.19*B(237)+0.25*B(239)

  JVS(766) = 0.25*B(240)

  JVS(767) = 0.19*B(238)

  JVS(768) = 0.167*B(425)

  JVS(769) = 0.4*B(228)-B(248)

  JVS(770) = 0.23*B(230)+0.167*B(421)

  JVS(771) = 0.25*B(232)+0.167*B(423)

  JVS(772) = -B(246)

  JVS(773) = 0.794*B(420)+0.794*B(422)+0.794*B(424)

  JVS(774) = 0

  JVS(775) = -B(27)-B(430)-B(432)

  JVS(776) = 0.1*B(386)

  JVS(777) = 0.1*B(408)+B(410)+0.5*B(412)

  JVS(778) = 0.1*B(400)

  JVS(779) = 0.8*B(251)

  JVS(780) = 0.08*B(229)

  JVS(781) = 0.5*B(413)

  JVS(782) = B(411)+0.794*B(425)

  JVS(783) = 0

  JVS(784) = 0.08*B(230)+0.8*B(252)+0.1*B(387)+0.1*B(401)+0.1*B(409)+0.794*B(421)

  JVS(785) = 0.794*B(423)-B(433)

  JVS(786) = -B(431)

  JVS(787) = B(42)

  JVS(788) = 0.45*B(36)

  JVS(789) = B(211)

  JVS(790) = 0.51*B(348)

  JVS(791) = 0.4*B(330)

  JVS(792) = 0.54*B(358)

  JVS(793) = 0.07*B(332)+0.17*B(334)

  JVS(794) = 0.07*B(336)+0.17*B(338)

  JVS(795) = 0.18*B(32)

  JVS(796) = 0

  JVS(797) = 0

  JVS(798) = 0

  JVS(799) = 0

  JVS(800) = 0

  JVS(801) = 0

  JVS(802) = 0.8*B(247)

  JVS(803) = -B(25)-B(426)-B(428)

  JVS(804) = 0.5*B(195)

  JVS(805) = B(458)

  JVS(806) = 0.25*B(436)+0.25*B(438)+0.1*B(442)+0.25*B(444)

  JVS(807) = 0.95*B(243)

  JVS(808) = 0.25*B(249)+0.25*B(253)+0.24*B(257)+0.25*B(259)

  JVS(809) = 0.02*B(229)+0.02*B(231)

  JVS(810) = 0.25*B(260)+0.25*B(445)

  JVS(811) = 0.5*B(196)+0.24*B(258)+0.1*B(443)

  JVS(812) = 0.07*B(333)+0.07*B(337)

  JVS(813) = 0.95*B(244)+0.8*B(248)

  JVS(814) = 0.02*B(230)+0.25*B(250)+0.4*B(331)+0.17*B(335)+0.17*B(339)+0.51*B(349)+0.54*B(359)+0.25*B(437)

  JVS(815) = 0

  JVS(816) = 0.02*B(232)+0.25*B(254)-B(429)+0.25*B(439)

  JVS(817) = B(212)-B(427)+B(459)

  JVS(818) = B(392)

  JVS(819) = B(360)

  JVS(820) = B(362)

  JVS(821) = B(368)

  JVS(822) = B(364)

  JVS(823) = B(366)

  JVS(824) = -B(386)-B(388)-B(390)

  JVS(825) = -B(391)

  JVS(826) = -B(389)

  JVS(827) = 0

  JVS(828) = -B(387)

  JVS(829) = 0

  JVS(830) = B(361)+B(363)+B(365)+B(367)+B(369)+B(393)

  JVS(831) = B(384)

  JVS(832) = B(380)

  JVS(833) = B(382)

  JVS(834) = 0.5*B(398)

  JVS(835) = 0

  JVS(836) = -B(408)-B(410)-B(412)-B(414)

  JVS(837) = -B(413)

  JVS(838) = -B(411)

  JVS(839) = 0

  JVS(840) = -B(409)

  JVS(841) = B(381)+B(383)+B(385)+0.5*B(399)-B(415)

  JVS(842) = 0

  JVS(843) = B(394)

  JVS(844) = B(396)+0.5*B(398)

  JVS(845) = 0

  JVS(846) = 0

  JVS(847) = -B(400)-B(402)-B(404)

  JVS(848) = -B(405)

  JVS(849) = -B(403)

  JVS(850) = 0

  JVS(851) = -B(401)

  JVS(852) = 0.5*B(399)

  JVS(853) = B(395)+B(397)

  JVS(854) = B(197)

  JVS(855) = 0.06*B(370)

  JVS(856) = 0.06*B(372)

  JVS(857) = 0.06*B(378)

  JVS(858) = 0.06*B(374)

  JVS(859) = 0.06*B(376)

  JVS(860) = B(189)

  JVS(861) = 0

  JVS(862) = 0

  JVS(863) = 0.15*B(46)+0.15*B(406)

  JVS(864) = 0

  JVS(865) = 0

  JVS(866) = 0

  JVS(867) = -B(191)-B(193)-B(195)

  JVS(868) = -B(196)

  JVS(869) = -B(194)

  JVS(870) = 0.06*B(371)+0.06*B(373)+0.06*B(375)+0.06*B(377)+0.06*B(379)

  JVS(871) = -B(192)

  JVS(872) = 0

  JVS(873) = B(190)+B(198)+0.15*B(407)

  JVS(874) = 0.5*B(161)

  JVS(875) = 0.5*B(454)

  JVS(876) = 0

  JVS(877) = 0

  JVS(878) = 0.2*B(195)

  JVS(879) = -B(29)-B(458)

  JVS(880) = 0.25*B(436)+0.25*B(438)+0.1*B(442)+0.25*B(444)

  JVS(881) = 0

  JVS(882) = 0.22*B(249)+0.22*B(253)+0.23*B(257)+0.22*B(259)

  JVS(883) = 0.02*B(229)+0.02*B(231)

  JVS(884) = 0.22*B(260)+0.25*B(445)

  JVS(885) = 0.2*B(196)+0.23*B(258)+0.1*B(443)

  JVS(886) = 0

  JVS(887) = 0

  JVS(888) = 0

  JVS(889) = 0.02*B(230)+0.22*B(250)+0.25*B(437)

  JVS(890) = 0

  JVS(891) = 0.02*B(232)+0.22*B(254)+0.25*B(439)

  JVS(892) = 0.5*B(162)+0.5*B(455)-B(459)

  JVS(893) = B(446)

  JVS(894) = B(434)

  JVS(895) = 0.5*B(235)

  JVS(896) = 0

  JVS(897) = 0

  JVS(898) = 0

  JVS(899) = 0

  JVS(900) = -B(436)-B(438)-B(440)-B(442)-B(444)

  JVS(901) = 0

  JVS(902) = 0

  JVS(903) = -B(445)

  JVS(904) = -B(443)

  JVS(905) = -B(441)

  JVS(906) = 0

  JVS(907) = -B(437)

  JVS(908) = -B(439)

  JVS(909) = 0.5*B(236)+B(435)+B(447)

  JVS(910) = 0.402*B(28)

  JVS(911) = 0.039*B(420)+0.039*B(422)+0.039*B(424)

  JVS(912) = 0.2*B(227)

  JVS(913) = -B(20)-B(241)-B(243)

  JVS(914) = 0.32*B(229)+0.35*B(231)+0.26*B(237)+0.35*B(239)

  JVS(915) = 0.35*B(240)

  JVS(916) = 0.26*B(238)

  JVS(917) = 0.039*B(425)

  JVS(918) = 0.2*B(228)-B(244)

  JVS(919) = 0.32*B(230)+0.039*B(421)

  JVS(920) = 0.35*B(232)+0.039*B(423)

  JVS(921) = -B(242)

  JVS(922) = 0.2*B(261)

  JVS(923) = 0.5*B(245)

  JVS(924) = B(241)

  JVS(925) = -B(249)-B(251)-B(253)-B(255)-B(257)-B(259)

  JVS(926) = 0

  JVS(927) = -B(260)

  JVS(928) = -B(258)

  JVS(929) = -B(256)

  JVS(930) = 0

  JVS(931) = -B(250)-B(252)

  JVS(932) = -B(254)

  JVS(933) = B(242)+0.5*B(246)+0.2*B(262)

  JVS(934) = 2*B(140)

  JVS(935) = B(23)

  JVS(936) = 0.5*B(201)

  JVS(937) = B(10)+0.3*B(120)

  JVS(938) = 0.5*B(128)+B(130)

  JVS(939) = B(504)

  JVS(940) = 0

  JVS(941) = 0.5*B(175)

  JVS(942) = 0.4*B(47)

  JVS(943) = B(456)

  JVS(944) = B(15)

  JVS(945) = 0.25*B(112)

  JVS(946) = 0.69*B(28)

  JVS(947) = 0.33*B(43)

  JVS(948) = 0.5*B(454)

  JVS(949) = 0.375*B(48)

  JVS(950) = 0.25*B(113)

  JVS(951) = 0.34*B(370)

  JVS(952) = 0.34*B(372)

  JVS(953) = 0.1*B(34)

  JVS(954) = B(157)

  JVS(955) = 0.5*B(488)

  JVS(956) = B(185)

  JVS(957) = B(450)

  JVS(958) = 0.7*B(147)

  JVS(959) = 0.072*B(420)+0.072*B(422)+0.206*B(424)

  JVS(960) = 0.34*B(378)

  JVS(961) = 0.34*B(374)

  JVS(962) = 0.25*B(494)+0.25*B(496)

  JVS(963) = 0.34*B(376)

  JVS(964) = 0.6*B(227)

  JVS(965) = B(30)+0.8*B(460)

  JVS(966) = 0.54*B(153)

  JVS(967) = 0.33*B(474)+0.917*B(476)

  JVS(968) = 0.1*B(205)

  JVS(969) = 0.68*B(46)+0.68*B(406)

  JVS(970) = 0.67*B(19)+0.7*B(247)

  JVS(971) = B(27)

  JVS(972) = 0.36*B(386)+0.95*B(390)

  JVS(973) = 0.75*B(412)

  JVS(974) = 0.34*B(400)+0.93*B(404)

  JVS(975) = B(191)+0.8*B(195)

  JVS(976) = B(29)

  JVS(977) = 0.25*B(436)+0.25*B(438)+0.8*B(442)+0.25*B(444)

  JVS(978) = 0.8*B(243)

  JVS(979) = 0.25*B(249)+0.25*B(253)+0.88*B(257)+0.25*B(259)

  JVS(980) = -B(11)-B(12)-B(122)-B(124)-B(502)

  JVS(981) = 0.55*B(229)+0.6*B(231)+1.2*B(237)+0.6*B(239)

  JVS(982) = B(173)+0.6*B(240)+0.25*B(260)+B(271)+0.25*B(445)

  JVS(983) = B(114)+2*B(116)+B(117)+0.7*B(148)+B(174)+B(186)+0.8*B(196)+1.2*B(238)+0.88*B(258)+2*B(269)+0.95*B(391)+0.93&
               &*B(405)+0.75*B(413)+0.8*B(443)+0.917*B(477)

  JVS(984) = B(263)+B(265)+2*B(270)+B(272)+2*B(273)

  JVS(985) = 0.206*B(425)-B(503)

  JVS(986) = B(131)+0.54*B(154)+0.6*B(228)+0.8*B(244)+0.7*B(248)+0.34*B(371)+0.34*B(373)+0.34*B(375)+0.34*B(377)+0.34&
               &*B(379)+0.5*B(489)

  JVS(987) = B(115)+B(158)+B(192)+0.5*B(202)+0.1*B(206)+0.55*B(230)+0.25*B(250)+B(264)+0.36*B(387)+0.34*B(401)+0.072&
               &*B(421)+0.25*B(437)+0.33*B(475)+0.25*B(495)

  JVS(988) = 0

  JVS(989) = -B(123)+0.6*B(232)+0.25*B(254)+B(266)+0.072*B(423)+0.25*B(439)+0.25*B(497)

  JVS(990) = 0.3*B(121)-B(125)+0.5*B(129)+0.5*B(176)+0.68*B(407)+B(451)+0.5*B(455)+B(457)+0.8*B(461)

  JVS(991) = 0.5*B(235)

  JVS(992) = B(225)

  JVS(993) = -B(229)-B(231)-B(233)-B(237)-B(239)

  JVS(994) = -B(240)

  JVS(995) = -B(238)

  JVS(996) = -B(234)

  JVS(997) = 0

  JVS(998) = -B(230)

  JVS(999) = -B(232)

  JVS(1000) = B(226)+0.5*B(236)

  JVS(1001) = 0.5*B(39)

  JVS(1002) = 0

  JVS(1003) = B(42)

  JVS(1004) = B(23)

  JVS(1005) = B(35)

  JVS(1006) = B(31)

  JVS(1007) = 0.5*B(175)

  JVS(1008) = 0.6*B(17)+B(177)

  JVS(1009) = 0.06*B(370)

  JVS(1010) = 0.06*B(372)

  JVS(1011) = 0

  JVS(1012) = B(215)

  JVS(1013) = 0

  JVS(1014) = 0.13*B(32)

  JVS(1015) = 0.06*B(378)

  JVS(1016) = 0.06*B(374)

  JVS(1017) = 0.06*B(376)

  JVS(1018) = B(24)

  JVS(1019) = 0

  JVS(1020) = B(163)+B(165)

  JVS(1021) = 0

  JVS(1022) = 0.65*B(46)+0.65*B(406)

  JVS(1023) = 0.67*B(19)

  JVS(1024) = B(25)+B(426)+B(428)

  JVS(1025) = 0

  JVS(1026) = 0

  JVS(1027) = 0

  JVS(1028) = B(191)+0.3*B(195)

  JVS(1029) = B(29)

  JVS(1030) = -B(444)

  JVS(1031) = 0.3*B(20)

  JVS(1032) = 0.53*B(249)+0.53*B(253)+0.26*B(257)-0.47*B(259)

  JVS(1033) = -B(239)

  JVS(1034) = -B(167)-B(169)-B(171)-B(173)-2*B(178)-B(240)-0.47*B(260)-B(445)

  JVS(1035) = -B(174)+0.3*B(196)+0.26*B(258)+B(269)

  JVS(1036) = B(263)+B(265)+B(270)+2*B(273)

  JVS(1037) = -B(172)

  JVS(1038) = 0.06*B(371)+0.06*B(373)+0.06*B(375)+0.06*B(377)+0.06*B(379)

  JVS(1039) = -B(168)+B(192)+B(216)+0.53*B(250)+B(264)

  JVS(1040) = -B(170)

  JVS(1041) = B(166)+0.53*B(254)+B(266)+B(429)

  JVS(1042) = B(164)+0.5*B(176)+0.65*B(407)+B(427)

  JVS(1043) = 0.7*B(120)

  JVS(1044) = B(416)

  JVS(1045) = B(16)

  JVS(1046) = 0.4*B(17)

  JVS(1047) = B(110)+0.75*B(112)

  JVS(1048) = 0.75*B(113)

  JVS(1049) = 0.33*B(332)+0.83*B(334)

  JVS(1050) = -B(185)

  JVS(1051) = 0.07*B(336)+0.17*B(338)

  JVS(1052) = -B(147)

  JVS(1053) = B(24)

  JVS(1054) = 0

  JVS(1055) = 0.31*B(153)

  JVS(1056) = -B(476)

  JVS(1057) = B(14)

  JVS(1058) = 0

  JVS(1059) = 0

  JVS(1060) = -B(390)

  JVS(1061) = -B(412)

  JVS(1062) = -B(404)

  JVS(1063) = -B(195)

  JVS(1064) = -B(442)+B(444)

  JVS(1065) = 0.3*B(20)

  JVS(1066) = -B(257)+B(259)

  JVS(1067) = 0

  JVS(1068) = -B(237)+B(239)

  JVS(1069) = B(167)-0.1*B(173)+2*B(178)+B(240)+B(260)+B(271)+B(445)

  JVS(1070) = -B(114)-2*B(116)-2*B(117)-B(118)-B(148)-0.1*B(174)-B(186)-B(196)-B(238)-B(258)-B(269)-B(391)-B(405)-B(413)&
                &-B(443)-B(477)

  JVS(1071) = -B(270)+B(272)

  JVS(1072) = -B(119)+0.33*B(333)+0.07*B(337)

  JVS(1073) = 0.31*B(154)

  JVS(1074) = -B(115)+B(168)+0.83*B(335)+0.17*B(339)

  JVS(1075) = 0

  JVS(1076) = 0

  JVS(1077) = B(111)+0.7*B(121)+B(417)

  JVS(1078) = 0.5*B(261)

  JVS(1079) = B(18)+B(277)

  JVS(1080) = 0.2*B(227)

  JVS(1081) = 0.33*B(19)+0.5*B(245)

  JVS(1082) = 0

  JVS(1083) = 0

  JVS(1084) = -B(271)

  JVS(1085) = -B(269)

  JVS(1086) = -B(263)-B(265)-B(267)-B(270)-B(272)-2*B(273)-B(274)

  JVS(1087) = -B(268)

  JVS(1088) = 0.2*B(228)

  JVS(1089) = -B(264)

  JVS(1090) = -B(275)

  JVS(1091) = -B(266)

  JVS(1092) = 0.5*B(246)+0.5*B(262)

  JVS(1093) = B(72)

  JVS(1094) = 0.35*B(498)

  JVS(1095) = 0.65*B(279)

  JVS(1096) = 1.5*B(38)

  JVS(1097) = 0.8*B(281)

  JVS(1098) = B(138)+B(140)

  JVS(1099) = 0.28*B(311)

  JVS(1100) = 0.73*B(313)

  JVS(1101) = 0.6*B(40)

  JVS(1102) = B(39)

  JVS(1103) = 0.38*B(344)

  JVS(1104) = 0.63*B(346)

  JVS(1105) = 0.5*B(464)

  JVS(1106) = B(42)

  JVS(1107) = B(452)

  JVS(1108) = 0.6*B(41)

  JVS(1109) = B(37)

  JVS(1110) = B(58)+B(60)

  JVS(1111) = B(21)

  JVS(1112) = B(22)

  JVS(1113) = B(201)

  JVS(1114) = B(10)

  JVS(1115) = 0.2*B(261)

  JVS(1116) = 0.66*B(9)+B(105)

  JVS(1117) = 0.25*B(128)+0.12*B(130)

  JVS(1118) = B(504)+B(505)-B(507)

  JVS(1119) = B(283)-B(285)

  JVS(1120) = 0

  JVS(1121) = -B(315)+B(319)

  JVS(1122) = B(47)

  JVS(1123) = B(299)-B(301)

  JVS(1124) = B(348)-B(350)

  JVS(1125) = B(15)

  JVS(1126) = 0.4*B(112)

  JVS(1127) = B(480)-0.6*B(482)

  JVS(1128) = -B(328)

  JVS(1129) = B(28)

  JVS(1130) = B(43)

  JVS(1131) = 0.5*B(454)

  JVS(1132) = B(500)

  JVS(1133) = B(48)

  JVS(1134) = B(59)+0.4*B(113)

  JVS(1135) = -B(295)

  JVS(1136) = 0.57*B(370)

  JVS(1137) = 0.57*B(372)

  JVS(1138) = 0.9*B(34)

  JVS(1139) = -B(217)

  JVS(1140) = B(157)-B(159)

  JVS(1141) = 0.4*B(307)-0.84*B(309)

  JVS(1142) = -B(221)+B(330)

  JVS(1143) = 0.06*B(488)

  JVS(1144) = -B(354)+B(358)

  JVS(1145) = -0.93*B(332)+0.17*B(334)

  JVS(1146) = B(62)-B(64)

  JVS(1147) = B(181)-B(183)+B(185)

  JVS(1148) = -0.67*B(336)+0.83*B(338)

  JVS(1149) = B(450)

  JVS(1150) = 2*B(33)+B(134)

  JVS(1151) = B(143)-B(145)+B(147)+1.2*B(471)

  JVS(1152) = 0.56*B(32)

  JVS(1153) = 0.794*B(420)+0.794*B(422)-B(424)

  JVS(1154) = 0.57*B(378)

  JVS(1155) = 0.57*B(374)

  JVS(1156) = -B(492)+0.75*B(494)+0.75*B(496)

  JVS(1157) = 2*B(44)

  JVS(1158) = 0.57*B(376)

  JVS(1159) = 0.06*B(227)

  JVS(1160) = B(126)

  JVS(1161) = 2*B(30)+B(460)

  JVS(1162) = 0.19*B(153)

  JVS(1163) = B(45)

  JVS(1164) = B(474)+B(476)-B(484)

  JVS(1165) = B(14)

  JVS(1166) = 0.9*B(205)-B(207)

  JVS(1167) = 1.2*B(46)+0.2*B(406)

  JVS(1168) = 0.67*B(19)+0.275*B(247)

  JVS(1169) = B(27)+B(430)+B(432)

  JVS(1170) = B(25)

  JVS(1171) = 0.9*B(386)-B(388)+B(390)

  JVS(1172) = -B(410)+0.5*B(412)

  JVS(1173) = 0.9*B(400)-B(402)+B(404)

  JVS(1174) = -B(193)+0.3*B(195)

  JVS(1175) = B(29)+B(458)

  JVS(1176) = B(436)+B(438)-B(440)+0.8*B(442)+B(444)

  JVS(1177) = 0.06*B(243)

  JVS(1178) = 0.47*B(249)+0.47*B(253)-B(255)+0.73*B(257)+0.47*B(259)

  JVS(1179) = 2*B(11)+B(122)+B(124)-B(502)

  JVS(1180) = 0.92*B(229)+B(231)-B(233)+B(237)+B(239)

  JVS(1181) = -B(171)+0.9*B(173)+B(240)+0.47*B(260)+B(445)

  JVS(1182) = B(114)+2*B(116)-B(118)+B(148)+0.9*B(174)+B(186)+0.3*B(196)+B(238)+0.73*B(258)+B(269)+B(391)+B(405)+0.5&
                &*B(413)+0.8*B(443)+B(477)

  JVS(1183) = -B(267)+B(270)

  JVS(1184) = -B(65)-B(68)-2*B(70)-B(74)-B(82)-B(90)-B(101)-B(119)-B(146)-B(160)-B(172)-B(184)-B(194)-B(208)-B(218)&
                &-B(222)-B(234)-B(256)-B(268)-B(286)-B(296)-B(302)-0.84*B(310)-B(316)-B(329)-0.93*B(333)-0.67*B(337)-B(351)&
                &-B(355)-B(389)-B(403)-B(411)-B(425)-B(441)-B(470)-0.6*B(483)-B(485)-B(493)-B(503)-B(508)

  JVS(1185) = B(66)-B(69)+0.12*B(131)+0.19*B(154)+0.06*B(228)+0.06*B(244)+0.275*B(248)+0.57*B(371)+0.57*B(373)+0.57&
                &*B(375)+0.57*B(377)+0.57*B(379)+0.06*B(489)

  JVS(1186) = -B(83)+B(115)+B(144)+B(158)+B(182)+B(202)+0.9*B(206)+0.92*B(230)+0.47*B(250)+B(284)+B(300)+0.4*B(308)&
                &+B(320)+B(331)+0.17*B(335)+0.83*B(339)+B(349)+B(359)+0.9*B(387)+0.9*B(401)+0.794*B(421)+B(437)+B(475)&
                &+B(481)+0.75*B(495)+B(506)

  JVS(1187) = -B(102)

  JVS(1188) = -B(91)+B(123)+B(232)+0.47*B(254)+0.794*B(423)+B(433)+B(439)+0.75*B(497)

  JVS(1189) = B(61)+B(63)+B(67)+B(73)-B(75)+B(125)+B(127)+0.25*B(129)+B(135)+0.2*B(262)+0.65*B(280)+0.8*B(282)+0.28&
                &*B(312)+0.73*B(314)+0.38*B(345)+0.63*B(347)+0.2*B(407)+B(431)+B(451)+B(453)+0.5*B(455)+B(459)+B(461)+0.5&
                &*B(465)+0.35*B(499)+B(501)

  JVS(1190) = -B(130)

  JVS(1191) = -B(291)

  JVS(1192) = 0

  JVS(1193) = -B(370)

  JVS(1194) = -B(372)

  JVS(1195) = -B(488)

  JVS(1196) = B(50)-B(52)

  JVS(1197) = -B(378)

  JVS(1198) = -B(374)

  JVS(1199) = -B(376)

  JVS(1200) = -0.9*B(227)

  JVS(1201) = -B(153)

  JVS(1202) = -0.8*B(247)

  JVS(1203) = -0.8*B(243)

  JVS(1204) = 0

  JVS(1205) = 0

  JVS(1206) = 0.25*B(171)

  JVS(1207) = 0

  JVS(1208) = 0.25*B(267)

  JVS(1209) = -B(68)+0.25*B(172)+0.25*B(268)

  JVS(1210) = -B(2)-B(3)-B(53)-B(66)-B(69)-B(84)-B(88)-B(131)-B(154)-0.9*B(228)-0.8*B(244)-0.8*B(248)-B(292)-B(371)&
                &-B(373)-B(375)-B(377)-B(379)-B(489)

  JVS(1211) = -B(85)

  JVS(1212) = -B(89)

  JVS(1213) = 0.89*B(8)

  JVS(1214) = -B(67)

  JVS(1215) = B(49)

  JVS(1216) = 2*B(78)

  JVS(1217) = -B(201)

  JVS(1218) = -B(505)

  JVS(1219) = -B(283)

  JVS(1220) = -B(136)

  JVS(1221) = -B(319)

  JVS(1222) = -B(299)

  JVS(1223) = -B(348)

  JVS(1224) = -B(480)

  JVS(1225) = -B(326)

  JVS(1226) = 2*B(79)

  JVS(1227) = -B(293)

  JVS(1228) = -B(215)

  JVS(1229) = -B(157)

  JVS(1230) = -B(307)

  JVS(1231) = -B(330)

  JVS(1232) = -B(358)

  JVS(1233) = -B(334)

  JVS(1234) = B(86)

  JVS(1235) = -B(181)

  JVS(1236) = -B(338)

  JVS(1237) = -B(143)

  JVS(1238) = -B(420)

  JVS(1239) = -B(494)

  JVS(1240) = 0

  JVS(1241) = 0

  JVS(1242) = 0

  JVS(1243) = -B(474)

  JVS(1244) = -B(205)

  JVS(1245) = -B(386)

  JVS(1246) = -B(408)

  JVS(1247) = -B(400)

  JVS(1248) = -B(191)

  JVS(1249) = -B(436)

  JVS(1250) = 0

  JVS(1251) = -B(249)-B(251)

  JVS(1252) = 0

  JVS(1253) = -B(229)

  JVS(1254) = -B(167)

  JVS(1255) = -B(114)

  JVS(1256) = -B(263)

  JVS(1257) = -B(82)

  JVS(1258) = -B(84)

  JVS(1259) = -B(83)-B(85)-B(99)-B(115)-B(137)-B(144)-B(158)-B(168)-B(182)-B(192)-B(202)-B(206)-B(216)-B(230)-B(250)&
                &-B(252)-B(264)-B(284)-B(294)-B(300)-B(308)-B(320)-B(327)-B(331)-B(335)-B(339)-B(349)-B(359)-B(387)-B(401)&
                &-B(409)-B(421)-B(437)-B(475)-B(481)-B(495)-B(506)

  JVS(1260) = B(5)+B(87)+0.5*B(109)

  JVS(1261) = 0.11*B(8)-B(100)

  JVS(1262) = 0

  JVS(1263) = B(325)

  JVS(1264) = B(6)+B(94)

  JVS(1265) = B(201)

  JVS(1266) = 0.66*B(9)+B(103)+B(105)

  JVS(1267) = B(505)

  JVS(1268) = B(283)

  JVS(1269) = -B(289)

  JVS(1270) = B(136)

  JVS(1271) = B(319)

  JVS(1272) = B(299)

  JVS(1273) = B(211)

  JVS(1274) = 0.6*B(17)+B(177)

  JVS(1275) = B(348)

  JVS(1276) = B(7)

  JVS(1277) = B(480)

  JVS(1278) = -B(323)+B(326)

  JVS(1279) = B(18)+B(277)

  JVS(1280) = B(293)

  JVS(1281) = B(215)

  JVS(1282) = B(157)

  JVS(1283) = -B(305)

  JVS(1284) = B(330)

  JVS(1285) = B(358)

  JVS(1286) = B(334)-B(340)

  JVS(1287) = -B(86)

  JVS(1288) = B(181)

  JVS(1289) = B(338)-B(342)

  JVS(1290) = B(143)

  JVS(1291) = 1.206*B(420)+1.206*B(422)+0.206*B(424)

  JVS(1292) = 1.25*B(494)+1.25*B(496)

  JVS(1293) = 0

  JVS(1294) = 0

  JVS(1295) = 0

  JVS(1296) = B(474)

  JVS(1297) = 0

  JVS(1298) = 0.9*B(205)

  JVS(1299) = B(27)+0.4*B(430)+B(432)

  JVS(1300) = 0

  JVS(1301) = 0.9*B(386)

  JVS(1302) = 1.9*B(408)+0.5*B(412)+2*B(414)

  JVS(1303) = 0.9*B(400)

  JVS(1304) = B(191)

  JVS(1305) = 0

  JVS(1306) = B(436)+B(438)

  JVS(1307) = 0

  JVS(1308) = B(249)+B(253)

  JVS(1309) = 0

  JVS(1310) = 0.92*B(229)+B(231)

  JVS(1311) = B(167)-B(169)

  JVS(1312) = B(114)+0.5*B(413)

  JVS(1313) = B(263)+B(265)-B(274)

  JVS(1314) = B(82)+B(90)-B(101)+0.206*B(425)

  JVS(1315) = B(84)-B(88)

  JVS(1316) = B(83)+B(85)+2*B(99)+B(115)+B(137)+B(144)+B(158)+B(168)+B(182)+B(192)+B(202)+0.9*B(206)+B(216)+0.92*B(230)&
                &+B(250)+B(264)+B(284)+B(294)+B(300)+B(320)+B(327)+B(331)+B(335)+B(339)+B(349)+B(359)+0.9*B(387)+0.9*B(401)&
                &+1.9*B(409)+1.206*B(421)+B(437)+B(475)+B(481)+1.25*B(495)+B(506)

  JVS(1317) = -B(5)-B(87)-B(89)-B(92)-B(95)-B(102)-B(109)-B(170)-B(275)-B(290)-B(306)-B(324)-B(341)-B(343)

  JVS(1318) = 0.89*B(8)+B(91)-B(93)+2*B(100)+B(232)+B(254)+B(266)+2*B(415)+1.206*B(423)+B(433)+B(439)+1.25*B(497)

  JVS(1319) = -B(96)+B(104)+B(212)+0.4*B(431)

  JVS(1320) = B(6)+B(94)

  JVS(1321) = -B(466)

  JVS(1322) = 0.33*B(9)

  JVS(1323) = 0.4*B(17)+B(456)

  JVS(1324) = B(97)

  JVS(1325) = 0.5*B(454)

  JVS(1326) = -B(490)

  JVS(1327) = -B(422)

  JVS(1328) = -B(384)

  JVS(1329) = -B(380)

  JVS(1330) = -B(496)

  JVS(1331) = -B(382)

  JVS(1332) = -B(418)

  JVS(1333) = -B(155)

  JVS(1334) = -B(398)

  JVS(1335) = -B(165)

  JVS(1336) = 0

  JVS(1337) = -B(432)

  JVS(1338) = -B(428)

  JVS(1339) = 0

  JVS(1340) = -B(414)

  JVS(1341) = 0

  JVS(1342) = 0

  JVS(1343) = 0

  JVS(1344) = -B(438)

  JVS(1345) = 0

  JVS(1346) = -B(253)

  JVS(1347) = -B(122)

  JVS(1348) = -B(231)

  JVS(1349) = 0

  JVS(1350) = 0

  JVS(1351) = -B(265)

  JVS(1352) = -B(90)

  JVS(1353) = B(88)

  JVS(1354) = -B(99)

  JVS(1355) = B(89)-B(92)

  JVS(1356) = -B(8)-B(91)-B(93)-B(100)-B(108)-B(123)-B(156)-B(166)-B(232)-B(254)-B(266)-B(381)-B(383)-B(385)-B(399)&
                &-B(415)-B(419)-B(423)-B(429)-B(433)-B(439)-B(467)-B(491)-B(497)

  JVS(1357) = B(98)+0.5*B(455)+B(457)

  JVS(1358) = B(49)

  JVS(1359) = -B(468)

  JVS(1360) = -B(141)

  JVS(1361) = -B(179)

  JVS(1362) = -B(199)

  JVS(1363) = -B(132)

  JVS(1364) = 2*B(13)-B(72)

  JVS(1365) = -0.35*B(498)

  JVS(1366) = -B(279)

  JVS(1367) = -B(281)

  JVS(1368) = -B(311)

  JVS(1369) = -B(313)

  JVS(1370) = B(26)-B(446)

  JVS(1371) = -B(287)

  JVS(1372) = -B(297)

  JVS(1373) = -B(303)

  JVS(1374) = -B(317)

  JVS(1375) = -B(321)

  JVS(1376) = -B(352)

  JVS(1377) = -B(356)

  JVS(1378) = -B(344)

  JVS(1379) = -B(346)

  JVS(1380) = -B(462)-B(464)

  JVS(1381) = -B(452)

  JVS(1382) = B(58)-B(60)

  JVS(1383) = B(21)-0.5*B(149)

  JVS(1384) = B(22)-B(187)

  JVS(1385) = B(23)-B(197)

  JVS(1386) = B(35)-B(219)

  JVS(1387) = B(10)-0.7*B(120)

  JVS(1388) = -0.9*B(261)

  JVS(1389) = 0.33*B(9)-B(103)

  JVS(1390) = -B(434)

  JVS(1391) = -B(203)

  JVS(1392) = -B(128)+0.12*B(130)

  JVS(1393) = B(36)-B(223)

  JVS(1394) = 0

  JVS(1395) = 0

  JVS(1396) = -B(213)

  JVS(1397) = -B(416)

  JVS(1398) = B(16)-B(175)

  JVS(1399) = B(47)-B(392)

  JVS(1400) = 0

  JVS(1401) = -B(211)

  JVS(1402) = -B(456)

  JVS(1403) = 0

  JVS(1404) = B(15)-0.5*B(161)

  JVS(1405) = B(7)-B(97)

  JVS(1406) = -B(110)+0.75*B(112)

  JVS(1407) = 0.4*B(482)

  JVS(1408) = 0.4*B(328)

  JVS(1409) = -B(235)

  JVS(1410) = B(43)-0.5*B(486)

  JVS(1411) = -B(454)

  JVS(1412) = -B(500)

  JVS(1413) = B(48)-B(394)

  JVS(1414) = 2*B(56)+B(59)+0.75*B(113)

  JVS(1415) = 0

  JVS(1416) = -B(360)+0.63*B(370)

  JVS(1417) = -B(362)+0.63*B(372)

  JVS(1418) = B(34)-B(209)

  JVS(1419) = 0

  JVS(1420) = 0

  JVS(1421) = 0

  JVS(1422) = -B(472)+0.06*B(488)

  JVS(1423) = 0

  JVS(1424) = 0.4*B(332)

  JVS(1425) = -B(62)+B(64)

  JVS(1426) = 0

  JVS(1427) = 0.4*B(336)

  JVS(1428) = -B(450)

  JVS(1429) = -B(134)

  JVS(1430) = 0

  JVS(1431) = 0

  JVS(1432) = 0.206*B(424)

  JVS(1433) = -B(368)+0.63*B(378)

  JVS(1434) = -B(364)+0.63*B(374)

  JVS(1435) = 0

  JVS(1436) = -B(478)

  JVS(1437) = -B(366)+0.63*B(376)

  JVS(1438) = -B(189)

  JVS(1439) = -B(225)+0.27*B(227)

  JVS(1440) = -B(126)

  JVS(1441) = -B(460)

  JVS(1442) = -B(151)+0.33*B(153)

  JVS(1443) = -B(396)

  JVS(1444) = 0

  JVS(1445) = -B(163)

  JVS(1446) = 0

  JVS(1447) = -B(406)

  JVS(1448) = 0.33*B(19)-B(245)+0.215*B(247)

  JVS(1449) = -B(430)

  JVS(1450) = -B(426)

  JVS(1451) = 0

  JVS(1452) = 0

  JVS(1453) = 0

  JVS(1454) = 0

  JVS(1455) = -B(458)

  JVS(1456) = 0

  JVS(1457) = -B(241)+0.08*B(243)

  JVS(1458) = 0

  JVS(1459) = -B(124)

  JVS(1460) = 0

  JVS(1461) = 0

  JVS(1462) = 0

  JVS(1463) = 0

  JVS(1464) = B(65)+B(68)-B(74)+B(82)+B(90)+0.4*B(329)+0.4*B(333)+0.4*B(337)+0.206*B(425)+0.4*B(483)

  JVS(1465) = -B(66)+B(69)+0.12*B(131)+0.33*B(154)+0.27*B(228)+0.08*B(244)+0.215*B(248)+0.63*B(371)+0.63*B(373)+0.63&
                &*B(375)+0.63*B(377)+0.63*B(379)+0.06*B(489)

  JVS(1466) = B(83)

  JVS(1467) = -B(95)+0.5*B(109)

  JVS(1468) = B(91)

  JVS(1469) = -B(61)-B(63)-B(67)-B(73)-B(75)-2*B(76)-2*B(77)-B(96)-B(98)-B(104)-B(111)-0.7*B(121)-B(125)-B(127)-B(129)&
                &-B(133)-B(135)-B(142)-0.5*B(150)-B(152)-0.5*B(162)-B(164)-B(176)-B(180)-B(188)-B(190)-B(198)-B(200)-B(204)&
                &-B(210)-B(212)-B(214)-B(220)-B(224)-B(226)-B(236)-B(242)-B(246)-0.9*B(262)-B(280)-B(282)-B(288)-B(298)&
                &-B(304)-B(312)-B(314)-B(318)-B(322)-B(345)-B(347)-B(353)-B(357)-B(361)-B(363)-B(365)-B(367)-B(369)-B(393)&
                &-B(395)-B(397)-B(407)-B(417)-B(427)-B(431)-B(435)-B(447)-B(451)-B(453)-B(455)-B(457)-B(459)-B(461)-B(463)&
                &-B(465)-B(469)-B(473)-B(479)-0.5*B(487)-0.35*B(499)-B(501)
      
END SUBROUTINE mozart_mosaic_4bin_Jac_SP














SUBROUTINE mozart_mosaic_4bin_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(1469), W(136), a
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
      
END SUBROUTINE mozart_mosaic_4bin_KppDecomp



SUBROUTINE mozart_mosaic_4bin_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(1469), W(136), a
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
      
END SUBROUTINE mozart_mosaic_4bin_KppDecompCmplx


SUBROUTINE mozart_mosaic_4bin_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(1469), X(136), sum

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
      
END SUBROUTINE mozart_mosaic_4bin_KppSolveIndirect


SUBROUTINE mozart_mosaic_4bin_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(1469), X(136), sum

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
      
END SUBROUTINE mozart_mosaic_4bin_KppSolveCmplx













SUBROUTINE mozart_mosaic_4bin_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(21) = X(21)-JVS(115)*X(20)
  X(22) = X(22)-JVS(118)*X(20)
  X(25) = X(25)-JVS(126)*X(24)
  X(31) = X(31)-JVS(149)*X(24)
  X(38) = X(38)-JVS(174)*X(37)
  X(43) = X(43)-JVS(191)*X(21)
  X(48) = X(48)-JVS(215)*X(14)
  X(58) = X(58)-JVS(262)*X(22)-JVS(263)*X(25)-JVS(264)*X(27)
  X(59) = X(59)-JVS(269)*X(22)-JVS(270)*X(25)-JVS(271)*X(38)
  X(60) = X(60)-JVS(278)*X(55)
  X(61) = X(61)-JVS(283)*X(24)-JVS(284)*X(32)-JVS(285)*X(37)
  X(66) = X(66)-JVS(313)*X(20)-JVS(314)*X(29)
  X(69) = X(69)-JVS(330)*X(34)-JVS(331)*X(38)
  X(71) = X(71)-JVS(340)*X(36)-JVS(341)*X(39)
  X(74) = X(74)-JVS(359)*X(17)-JVS(360)*X(33)-JVS(361)*X(61)
  X(78) = X(78)-JVS(379)*X(19)-JVS(380)*X(55)-JVS(381)*X(56)
  X(80) = X(80)-JVS(398)*X(15)-JVS(399)*X(44)-JVS(400)*X(72)
  X(81) = X(81)-JVS(406)*X(28)-JVS(407)*X(59)-JVS(408)*X(74)
  X(85) = X(85)-JVS(425)*X(49)-JVS(426)*X(62)-JVS(427)*X(84)
  X(86) = X(86)-JVS(433)*X(70)
  X(87) = X(87)-JVS(439)*X(43)-JVS(440)*X(66)
  X(88) = X(88)-JVS(448)*X(24)-JVS(449)*X(57)
  X(90) = X(90)-JVS(458)*X(35)-JVS(459)*X(37)
  X(91) = X(91)-JVS(464)*X(30)-JVS(465)*X(88)-JVS(466)*X(90)
  X(92) = X(92)-JVS(472)*X(80)
  X(93) = X(93)-JVS(480)*X(13)-JVS(481)*X(46)
  X(94) = X(94)-JVS(487)*X(42)-JVS(488)*X(88)-JVS(489)*X(90)
  X(96) = X(96)-JVS(507)*X(19)-JVS(508)*X(57)-JVS(509)*X(58)-JVS(510)*X(66)-JVS(511)*X(69)-JVS(512)*X(87)-JVS(513)*X(88)&
            &-JVS(514)*X(90)-JVS(515)*X(94)
  X(97) = X(97)-JVS(528)*X(12)-JVS(529)*X(45)-JVS(530)*X(62)-JVS(531)*X(84)
  X(98) = X(98)-JVS(538)*X(57)-JVS(539)*X(82)-JVS(540)*X(83)-JVS(541)*X(88)
  X(102) = X(102)-JVS(563)*X(89)
  X(103) = X(103)-JVS(570)*X(76)-JVS(571)*X(89)-JVS(572)*X(102)
  X(105) = X(105)-JVS(585)*X(46)-JVS(586)*X(48)-JVS(587)*X(65)-JVS(588)*X(73)-JVS(589)*X(76)-JVS(590)*X(79)-JVS(591)&
             &*X(82)-JVS(592)*X(83)-JVS(593)*X(84)-JVS(594)*X(89)-JVS(595)*X(93)-JVS(596)*X(100)-JVS(597)*X(101)-JVS(598)&
             &*X(102)-JVS(599)*X(103)-JVS(600)*X(104)
  X(107) = X(107)-JVS(617)*X(19)-JVS(618)*X(21)-JVS(619)*X(31)-JVS(620)*X(37)-JVS(621)*X(40)-JVS(622)*X(42)-JVS(623)&
             &*X(55)-JVS(624)*X(79)-JVS(625)*X(82)-JVS(626)*X(83)-JVS(627)*X(87)-JVS(628)*X(88)-JVS(629)*X(89)-JVS(630)&
             &*X(90)-JVS(631)*X(91)-JVS(632)*X(94)-JVS(633)*X(96)-JVS(634)*X(98)-JVS(635)*X(100)-JVS(636)*X(101)-JVS(637)&
             &*X(103)-JVS(638)*X(104)-JVS(639)*X(106)
  X(108) = X(108)-JVS(664)*X(23)-JVS(665)*X(60)-JVS(666)*X(76)-JVS(667)*X(79)
  X(109) = X(109)-JVS(681)*X(106)
  X(110) = X(110)-JVS(687)*X(65)-JVS(688)*X(82)-JVS(689)*X(83)-JVS(690)*X(100)-JVS(691)*X(101)-JVS(692)*X(104)
  X(111) = X(111)-JVS(702)*X(76)-JVS(703)*X(89)
  X(112) = X(112)-JVS(711)*X(41)-JVS(712)*X(45)-JVS(713)*X(48)-JVS(714)*X(49)-JVS(715)*X(70)-JVS(716)*X(84)-JVS(717)&
             &*X(85)-JVS(718)*X(86)-JVS(719)*X(93)-JVS(720)*X(97)-JVS(721)*X(109)
  X(113) = X(113)-JVS(731)*X(54)-JVS(732)*X(82)-JVS(733)*X(83)-JVS(734)*X(84)-JVS(735)*X(100)-JVS(736)*X(101)-JVS(737)&
             &*X(104)
  X(114) = X(114)-JVS(744)*X(79)-JVS(745)*X(82)-JVS(746)*X(83)-JVS(747)*X(100)-JVS(748)*X(101)-JVS(749)*X(104)-JVS(750)&
             &*X(110)
  X(115) = X(115)-JVS(761)*X(75)-JVS(762)*X(99)-JVS(763)*X(106)
  X(116) = X(116)-JVS(773)*X(99)-JVS(774)*X(106)
  X(117) = X(117)-JVS(787)*X(40)-JVS(788)*X(57)-JVS(789)*X(67)-JVS(790)*X(69)-JVS(791)*X(88)-JVS(792)*X(90)-JVS(793)&
             &*X(91)-JVS(794)*X(94)-JVS(795)*X(98)-JVS(796)*X(100)-JVS(797)*X(101)-JVS(798)*X(102)-JVS(799)*X(104)-JVS(800)&
             &*X(109)-JVS(801)*X(113)-JVS(802)*X(115)
  X(118) = X(118)-JVS(818)*X(65)-JVS(819)*X(82)-JVS(820)*X(83)-JVS(821)*X(100)-JVS(822)*X(101)-JVS(823)*X(104)
  X(119) = X(119)-JVS(831)*X(100)-JVS(832)*X(101)-JVS(833)*X(104)-JVS(834)*X(110)-JVS(835)*X(118)
  X(120) = X(120)-JVS(843)*X(79)-JVS(844)*X(110)-JVS(845)*X(118)-JVS(846)*X(119)
  X(121) = X(121)-JVS(854)*X(47)-JVS(855)*X(82)-JVS(856)*X(83)-JVS(857)*X(100)-JVS(858)*X(101)-JVS(859)*X(104)-JVS(860)&
             &*X(105)-JVS(861)*X(111)-JVS(862)*X(113)-JVS(863)*X(114)-JVS(864)*X(118)-JVS(865)*X(119)-JVS(866)*X(120)
  X(122) = X(122)-JVS(874)*X(70)-JVS(875)*X(77)-JVS(876)*X(86)-JVS(877)*X(109)-JVS(878)*X(121)
  X(123) = X(123)-JVS(893)*X(26)-JVS(894)*X(53)-JVS(895)*X(75)-JVS(896)*X(116)-JVS(897)*X(118)-JVS(898)*X(119)-JVS(899)&
             &*X(120)
  X(124) = X(124)-JVS(910)*X(75)-JVS(911)*X(99)-JVS(912)*X(106)
  X(125) = X(125)-JVS(922)*X(51)-JVS(923)*X(115)-JVS(924)*X(124)
  X(126) = X(126)-JVS(934)*X(23)-JVS(935)*X(47)-JVS(936)*X(48)-JVS(937)*X(50)-JVS(938)*X(55)-JVS(939)*X(56)-JVS(940)&
             &*X(60)-JVS(941)*X(64)-JVS(942)*X(65)-JVS(943)*X(68)-JVS(944)*X(70)-JVS(945)*X(72)-JVS(946)*X(75)-JVS(947)&
             &*X(76)-JVS(948)*X(77)-JVS(949)*X(79)-JVS(950)*X(80)-JVS(951)*X(82)-JVS(952)*X(83)-JVS(953)*X(84)-JVS(954)&
             &*X(86)-JVS(955)*X(89)-JVS(956)*X(93)-JVS(957)*X(95)-JVS(958)*X(97)-JVS(959)*X(99)-JVS(960)*X(100)-JVS(961)&
             &*X(101)-JVS(962)*X(102)-JVS(963)*X(104)-JVS(964)*X(106)-JVS(965)*X(108)-JVS(966)*X(109)-JVS(967)*X(111)&
             &-JVS(968)*X(113)-JVS(969)*X(114)-JVS(970)*X(115)-JVS(971)*X(116)-JVS(972)*X(118)-JVS(973)*X(119)-JVS(974)&
             &*X(120)-JVS(975)*X(121)-JVS(976)*X(122)-JVS(977)*X(123)-JVS(978)*X(124)-JVS(979)*X(125)
  X(127) = X(127)-JVS(991)*X(75)-JVS(992)*X(106)
  X(128) = X(128)-JVS(1001)*X(31)-JVS(1002)*X(37)-JVS(1003)*X(40)-JVS(1004)*X(47)-JVS(1005)*X(49)-JVS(1006)*X(62)&
             &-JVS(1007)*X(64)-JVS(1008)*X(68)-JVS(1009)*X(82)-JVS(1010)*X(83)-JVS(1011)*X(84)-JVS(1012)*X(85)-JVS(1013)&
             &*X(90)-JVS(1014)*X(98)-JVS(1015)*X(100)-JVS(1016)*X(101)-JVS(1017)*X(104)-JVS(1018)*X(105)-JVS(1019)*X(111)&
             &-JVS(1020)*X(112)-JVS(1021)*X(113)-JVS(1022)*X(114)-JVS(1023)*X(115)-JVS(1024)*X(117)-JVS(1025)*X(118)&
             &-JVS(1026)*X(119)-JVS(1027)*X(120)-JVS(1028)*X(121)-JVS(1029)*X(122)-JVS(1030)*X(123)-JVS(1031)*X(124)&
             &-JVS(1032)*X(125)-JVS(1033)*X(127)
  X(129) = X(129)-JVS(1043)*X(50)-JVS(1044)*X(63)-JVS(1045)*X(64)-JVS(1046)*X(68)-JVS(1047)*X(72)-JVS(1048)*X(80)&
             &-JVS(1049)*X(91)-JVS(1050)*X(93)-JVS(1051)*X(94)-JVS(1052)*X(97)-JVS(1053)*X(105)-JVS(1054)*X(106)-JVS(1055)&
             &*X(109)-JVS(1056)*X(111)-JVS(1057)*X(112)-JVS(1058)*X(113)-JVS(1059)*X(114)-JVS(1060)*X(118)-JVS(1061)*X(119)&
             &-JVS(1062)*X(120)-JVS(1063)*X(121)-JVS(1064)*X(123)-JVS(1065)*X(124)-JVS(1066)*X(125)-JVS(1067)*X(126)&
             &-JVS(1068)*X(127)-JVS(1069)*X(128)
  X(130) = X(130)-JVS(1078)*X(51)-JVS(1079)*X(77)-JVS(1080)*X(106)-JVS(1081)*X(115)-JVS(1082)*X(125)-JVS(1083)*X(127)&
             &-JVS(1084)*X(128)-JVS(1085)*X(129)
  X(131) = X(131)-JVS(1093)*X(18)-JVS(1094)*X(19)-JVS(1095)*X(20)-JVS(1096)*X(21)-JVS(1097)*X(22)-JVS(1098)*X(23)&
             &-JVS(1099)*X(24)-JVS(1100)*X(25)-JVS(1101)*X(30)-JVS(1102)*X(31)-JVS(1103)*X(37)-JVS(1104)*X(38)-JVS(1105)&
             &*X(39)-JVS(1106)*X(40)-JVS(1107)*X(41)-JVS(1108)*X(42)-JVS(1109)*X(43)-JVS(1110)*X(44)-JVS(1111)*X(45)&
             &-JVS(1112)*X(46)-JVS(1113)*X(48)-JVS(1114)*X(50)-JVS(1115)*X(51)-JVS(1116)*X(52)-JVS(1117)*X(55)-JVS(1118)&
             &*X(56)-JVS(1119)*X(58)-JVS(1120)*X(60)-JVS(1121)*X(61)-JVS(1122)*X(65)-JVS(1123)*X(66)-JVS(1124)*X(69)&
             &-JVS(1125)*X(70)-JVS(1126)*X(72)-JVS(1127)*X(73)-JVS(1128)*X(74)-JVS(1129)*X(75)-JVS(1130)*X(76)-JVS(1131)&
             &*X(77)-JVS(1132)*X(78)-JVS(1133)*X(79)-JVS(1134)*X(80)-JVS(1135)*X(81)-JVS(1136)*X(82)-JVS(1137)*X(83)&
             &-JVS(1138)*X(84)-JVS(1139)*X(85)-JVS(1140)*X(86)-JVS(1141)*X(87)-JVS(1142)*X(88)-JVS(1143)*X(89)-JVS(1144)&
             &*X(90)-JVS(1145)*X(91)-JVS(1146)*X(92)-JVS(1147)*X(93)-JVS(1148)*X(94)-JVS(1149)*X(95)-JVS(1150)*X(96)&
             &-JVS(1151)*X(97)-JVS(1152)*X(98)-JVS(1153)*X(99)-JVS(1154)*X(100)-JVS(1155)*X(101)-JVS(1156)*X(102)-JVS(1157)&
             &*X(103)-JVS(1158)*X(104)-JVS(1159)*X(106)-JVS(1160)*X(107)-JVS(1161)*X(108)-JVS(1162)*X(109)-JVS(1163)*X(110)&
             &-JVS(1164)*X(111)-JVS(1165)*X(112)-JVS(1166)*X(113)-JVS(1167)*X(114)-JVS(1168)*X(115)-JVS(1169)*X(116)&
             &-JVS(1170)*X(117)-JVS(1171)*X(118)-JVS(1172)*X(119)-JVS(1173)*X(120)-JVS(1174)*X(121)-JVS(1175)*X(122)&
             &-JVS(1176)*X(123)-JVS(1177)*X(124)-JVS(1178)*X(125)-JVS(1179)*X(126)-JVS(1180)*X(127)-JVS(1181)*X(128)&
             &-JVS(1182)*X(129)-JVS(1183)*X(130)
  X(132) = X(132)-JVS(1190)*X(55)-JVS(1191)*X(59)-JVS(1192)*X(81)-JVS(1193)*X(82)-JVS(1194)*X(83)-JVS(1195)*X(89)&
             &-JVS(1196)*X(92)-JVS(1197)*X(100)-JVS(1198)*X(101)-JVS(1199)*X(104)-JVS(1200)*X(106)-JVS(1201)*X(109)&
             &-JVS(1202)*X(115)-JVS(1203)*X(124)-JVS(1204)*X(126)-JVS(1205)*X(127)-JVS(1206)*X(128)-JVS(1207)*X(129)&
             &-JVS(1208)*X(130)-JVS(1209)*X(131)
  X(133) = X(133)-JVS(1215)*X(1)-JVS(1216)*X(15)-JVS(1217)*X(48)-JVS(1218)*X(56)-JVS(1219)*X(58)-JVS(1220)*X(60)&
             &-JVS(1221)*X(61)-JVS(1222)*X(66)-JVS(1223)*X(69)-JVS(1224)*X(73)-JVS(1225)*X(74)-JVS(1226)*X(80)-JVS(1227)&
             &*X(81)-JVS(1228)*X(85)-JVS(1229)*X(86)-JVS(1230)*X(87)-JVS(1231)*X(88)-JVS(1232)*X(90)-JVS(1233)*X(91)&
             &-JVS(1234)*X(92)-JVS(1235)*X(93)-JVS(1236)*X(94)-JVS(1237)*X(97)-JVS(1238)*X(99)-JVS(1239)*X(102)-JVS(1240)&
             &*X(103)-JVS(1241)*X(106)-JVS(1242)*X(109)-JVS(1243)*X(111)-JVS(1244)*X(113)-JVS(1245)*X(118)-JVS(1246)*X(119)&
             &-JVS(1247)*X(120)-JVS(1248)*X(121)-JVS(1249)*X(123)-JVS(1250)*X(124)-JVS(1251)*X(125)-JVS(1252)*X(126)&
             &-JVS(1253)*X(127)-JVS(1254)*X(128)-JVS(1255)*X(129)-JVS(1256)*X(130)-JVS(1257)*X(131)-JVS(1258)*X(132)
  X(134) = X(134)-JVS(1263)*X(17)-JVS(1264)*X(36)-JVS(1265)*X(48)-JVS(1266)*X(52)-JVS(1267)*X(56)-JVS(1268)*X(58)&
             &-JVS(1269)*X(59)-JVS(1270)*X(60)-JVS(1271)*X(61)-JVS(1272)*X(66)-JVS(1273)*X(67)-JVS(1274)*X(68)-JVS(1275)&
             &*X(69)-JVS(1276)*X(71)-JVS(1277)*X(73)-JVS(1278)*X(74)-JVS(1279)*X(77)-JVS(1280)*X(81)-JVS(1281)*X(85)&
             &-JVS(1282)*X(86)-JVS(1283)*X(87)-JVS(1284)*X(88)-JVS(1285)*X(90)-JVS(1286)*X(91)-JVS(1287)*X(92)-JVS(1288)&
             &*X(93)-JVS(1289)*X(94)-JVS(1290)*X(97)-JVS(1291)*X(99)-JVS(1292)*X(102)-JVS(1293)*X(103)-JVS(1294)*X(106)&
             &-JVS(1295)*X(109)-JVS(1296)*X(111)-JVS(1297)*X(112)-JVS(1298)*X(113)-JVS(1299)*X(116)-JVS(1300)*X(117)&
             &-JVS(1301)*X(118)-JVS(1302)*X(119)-JVS(1303)*X(120)-JVS(1304)*X(121)-JVS(1305)*X(122)-JVS(1306)*X(123)&
             &-JVS(1307)*X(124)-JVS(1308)*X(125)-JVS(1309)*X(126)-JVS(1310)*X(127)-JVS(1311)*X(128)-JVS(1312)*X(129)&
             &-JVS(1313)*X(130)-JVS(1314)*X(131)-JVS(1315)*X(132)-JVS(1316)*X(133)
  X(135) = X(135)-JVS(1320)*X(36)-JVS(1321)*X(39)-JVS(1322)*X(52)-JVS(1323)*X(68)-JVS(1324)*X(71)-JVS(1325)*X(77)&
             &-JVS(1326)*X(89)-JVS(1327)*X(99)-JVS(1328)*X(100)-JVS(1329)*X(101)-JVS(1330)*X(102)-JVS(1331)*X(104)-JVS(1332)&
             &*X(106)-JVS(1333)*X(109)-JVS(1334)*X(110)-JVS(1335)*X(112)-JVS(1336)*X(113)-JVS(1337)*X(116)-JVS(1338)*X(117)&
             &-JVS(1339)*X(118)-JVS(1340)*X(119)-JVS(1341)*X(120)-JVS(1342)*X(121)-JVS(1343)*X(122)-JVS(1344)*X(123)&
             &-JVS(1345)*X(124)-JVS(1346)*X(125)-JVS(1347)*X(126)-JVS(1348)*X(127)-JVS(1349)*X(128)-JVS(1350)*X(129)&
             &-JVS(1351)*X(130)-JVS(1352)*X(131)-JVS(1353)*X(132)-JVS(1354)*X(133)-JVS(1355)*X(134)
  X(136) = X(136)-JVS(1358)*X(1)-JVS(1359)*X(11)-JVS(1360)*X(12)-JVS(1361)*X(13)-JVS(1362)*X(14)-JVS(1363)*X(16)&
             &-JVS(1364)*X(18)-JVS(1365)*X(19)-JVS(1366)*X(20)-JVS(1367)*X(22)-JVS(1368)*X(24)-JVS(1369)*X(25)-JVS(1370)&
             &*X(26)-JVS(1371)*X(27)-JVS(1372)*X(28)-JVS(1373)*X(29)-JVS(1374)*X(32)-JVS(1375)*X(33)-JVS(1376)*X(34)&
             &-JVS(1377)*X(35)-JVS(1378)*X(37)-JVS(1379)*X(38)-JVS(1380)*X(39)-JVS(1381)*X(41)-JVS(1382)*X(44)-JVS(1383)&
             &*X(45)-JVS(1384)*X(46)-JVS(1385)*X(47)-JVS(1386)*X(49)-JVS(1387)*X(50)-JVS(1388)*X(51)-JVS(1389)*X(52)&
             &-JVS(1390)*X(53)-JVS(1391)*X(54)-JVS(1392)*X(55)-JVS(1393)*X(57)-JVS(1394)*X(58)-JVS(1395)*X(61)-JVS(1396)&
             &*X(62)-JVS(1397)*X(63)-JVS(1398)*X(64)-JVS(1399)*X(65)-JVS(1400)*X(66)-JVS(1401)*X(67)-JVS(1402)*X(68)&
             &-JVS(1403)*X(69)-JVS(1404)*X(70)-JVS(1405)*X(71)-JVS(1406)*X(72)-JVS(1407)*X(73)-JVS(1408)*X(74)-JVS(1409)&
             &*X(75)-JVS(1410)*X(76)-JVS(1411)*X(77)-JVS(1412)*X(78)-JVS(1413)*X(79)-JVS(1414)*X(80)-JVS(1415)*X(81)&
             &-JVS(1416)*X(82)-JVS(1417)*X(83)-JVS(1418)*X(84)-JVS(1419)*X(85)-JVS(1420)*X(86)-JVS(1421)*X(88)-JVS(1422)&
             &*X(89)-JVS(1423)*X(90)-JVS(1424)*X(91)-JVS(1425)*X(92)-JVS(1426)*X(93)-JVS(1427)*X(94)-JVS(1428)*X(95)&
             &-JVS(1429)*X(96)-JVS(1430)*X(97)-JVS(1431)*X(98)-JVS(1432)*X(99)-JVS(1433)*X(100)-JVS(1434)*X(101)-JVS(1435)&
             &*X(102)-JVS(1436)*X(103)-JVS(1437)*X(104)-JVS(1438)*X(105)-JVS(1439)*X(106)-JVS(1440)*X(107)-JVS(1441)*X(108)&
             &-JVS(1442)*X(109)-JVS(1443)*X(110)-JVS(1444)*X(111)-JVS(1445)*X(112)-JVS(1446)*X(113)-JVS(1447)*X(114)&
             &-JVS(1448)*X(115)-JVS(1449)*X(116)-JVS(1450)*X(117)-JVS(1451)*X(118)-JVS(1452)*X(119)-JVS(1453)*X(120)&
             &-JVS(1454)*X(121)-JVS(1455)*X(122)-JVS(1456)*X(123)-JVS(1457)*X(124)-JVS(1458)*X(125)-JVS(1459)*X(126)&
             &-JVS(1460)*X(127)-JVS(1461)*X(128)-JVS(1462)*X(129)-JVS(1463)*X(130)-JVS(1464)*X(131)-JVS(1465)*X(132)&
             &-JVS(1466)*X(133)-JVS(1467)*X(134)-JVS(1468)*X(135)
  X(136) = X(136)/JVS(1469)
  X(135) = (X(135)-JVS(1357)*X(136))/(JVS(1356))
  X(134) = (X(134)-JVS(1318)*X(135)-JVS(1319)*X(136))/(JVS(1317))
  X(133) = (X(133)-JVS(1260)*X(134)-JVS(1261)*X(135)-JVS(1262)*X(136))/(JVS(1259))
  X(132) = (X(132)-JVS(1211)*X(133)-JVS(1212)*X(134)-JVS(1213)*X(135)-JVS(1214)*X(136))/(JVS(1210))
  X(131) = (X(131)-JVS(1185)*X(132)-JVS(1186)*X(133)-JVS(1187)*X(134)-JVS(1188)*X(135)-JVS(1189)*X(136))/(JVS(1184))
  X(130) = (X(130)-JVS(1087)*X(131)-JVS(1088)*X(132)-JVS(1089)*X(133)-JVS(1090)*X(134)-JVS(1091)*X(135)-JVS(1092)&
             &*X(136))/(JVS(1086))
  X(129) = (X(129)-JVS(1071)*X(130)-JVS(1072)*X(131)-JVS(1073)*X(132)-JVS(1074)*X(133)-JVS(1075)*X(134)-JVS(1076)*X(135)&
             &-JVS(1077)*X(136))/(JVS(1070))
  X(128) = (X(128)-JVS(1035)*X(129)-JVS(1036)*X(130)-JVS(1037)*X(131)-JVS(1038)*X(132)-JVS(1039)*X(133)-JVS(1040)*X(134)&
             &-JVS(1041)*X(135)-JVS(1042)*X(136))/(JVS(1034))
  X(127) = (X(127)-JVS(994)*X(128)-JVS(995)*X(129)-JVS(996)*X(131)-JVS(997)*X(132)-JVS(998)*X(133)-JVS(999)*X(135)&
             &-JVS(1000)*X(136))/(JVS(993))
  X(126) = (X(126)-JVS(981)*X(127)-JVS(982)*X(128)-JVS(983)*X(129)-JVS(984)*X(130)-JVS(985)*X(131)-JVS(986)*X(132)&
             &-JVS(987)*X(133)-JVS(988)*X(134)-JVS(989)*X(135)-JVS(990)*X(136))/(JVS(980))
  X(125) = (X(125)-JVS(926)*X(127)-JVS(927)*X(128)-JVS(928)*X(129)-JVS(929)*X(131)-JVS(930)*X(132)-JVS(931)*X(133)&
             &-JVS(932)*X(135)-JVS(933)*X(136))/(JVS(925))
  X(124) = (X(124)-JVS(914)*X(127)-JVS(915)*X(128)-JVS(916)*X(129)-JVS(917)*X(131)-JVS(918)*X(132)-JVS(919)*X(133)&
             &-JVS(920)*X(135)-JVS(921)*X(136))/(JVS(913))
  X(123) = (X(123)-JVS(901)*X(125)-JVS(902)*X(127)-JVS(903)*X(128)-JVS(904)*X(129)-JVS(905)*X(131)-JVS(906)*X(132)&
             &-JVS(907)*X(133)-JVS(908)*X(135)-JVS(909)*X(136))/(JVS(900))
  X(122) = (X(122)-JVS(880)*X(123)-JVS(881)*X(124)-JVS(882)*X(125)-JVS(883)*X(127)-JVS(884)*X(128)-JVS(885)*X(129)&
             &-JVS(886)*X(130)-JVS(887)*X(131)-JVS(888)*X(132)-JVS(889)*X(133)-JVS(890)*X(134)-JVS(891)*X(135)-JVS(892)&
             &*X(136))/(JVS(879))
  X(121) = (X(121)-JVS(868)*X(129)-JVS(869)*X(131)-JVS(870)*X(132)-JVS(871)*X(133)-JVS(872)*X(135)-JVS(873)*X(136))&
             &/(JVS(867))
  X(120) = (X(120)-JVS(848)*X(129)-JVS(849)*X(131)-JVS(850)*X(132)-JVS(851)*X(133)-JVS(852)*X(135)-JVS(853)*X(136))&
             &/(JVS(847))
  X(119) = (X(119)-JVS(837)*X(129)-JVS(838)*X(131)-JVS(839)*X(132)-JVS(840)*X(133)-JVS(841)*X(135)-JVS(842)*X(136))&
             &/(JVS(836))
  X(118) = (X(118)-JVS(825)*X(129)-JVS(826)*X(131)-JVS(827)*X(132)-JVS(828)*X(133)-JVS(829)*X(135)-JVS(830)*X(136))&
             &/(JVS(824))
  X(117) = (X(117)-JVS(804)*X(121)-JVS(805)*X(122)-JVS(806)*X(123)-JVS(807)*X(124)-JVS(808)*X(125)-JVS(809)*X(127)&
             &-JVS(810)*X(128)-JVS(811)*X(129)-JVS(812)*X(131)-JVS(813)*X(132)-JVS(814)*X(133)-JVS(815)*X(134)-JVS(816)&
             &*X(135)-JVS(817)*X(136))/(JVS(803))
  X(116) = (X(116)-JVS(776)*X(118)-JVS(777)*X(119)-JVS(778)*X(120)-JVS(779)*X(125)-JVS(780)*X(127)-JVS(781)*X(129)&
             &-JVS(782)*X(131)-JVS(783)*X(132)-JVS(784)*X(133)-JVS(785)*X(135)-JVS(786)*X(136))/(JVS(775))
  X(115) = (X(115)-JVS(765)*X(127)-JVS(766)*X(128)-JVS(767)*X(129)-JVS(768)*X(131)-JVS(769)*X(132)-JVS(770)*X(133)&
             &-JVS(771)*X(135)-JVS(772)*X(136))/(JVS(764))
  X(114) = (X(114)-JVS(752)*X(118)-JVS(753)*X(119)-JVS(754)*X(120)-JVS(755)*X(129)-JVS(756)*X(131)-JVS(757)*X(132)&
             &-JVS(758)*X(133)-JVS(759)*X(135)-JVS(760)*X(136))/(JVS(751))
  X(113) = (X(113)-JVS(739)*X(131)-JVS(740)*X(132)-JVS(741)*X(133)-JVS(742)*X(135)-JVS(743)*X(136))/(JVS(738))
  X(112) = (X(112)-JVS(723)*X(113)-JVS(724)*X(124)-JVS(725)*X(129)-JVS(726)*X(131)-JVS(727)*X(132)-JVS(728)*X(133)&
             &-JVS(729)*X(135)-JVS(730)*X(136))/(JVS(722))
  X(111) = (X(111)-JVS(705)*X(129)-JVS(706)*X(131)-JVS(707)*X(132)-JVS(708)*X(133)-JVS(709)*X(135)-JVS(710)*X(136))&
             &/(JVS(704))
  X(110) = (X(110)-JVS(694)*X(118)-JVS(695)*X(119)-JVS(696)*X(129)-JVS(697)*X(131)-JVS(698)*X(132)-JVS(699)*X(133)&
             &-JVS(700)*X(135)-JVS(701)*X(136))/(JVS(693))
  X(109) = (X(109)-JVS(683)*X(124)-JVS(684)*X(132)-JVS(685)*X(135)-JVS(686)*X(136))/(JVS(682))
  X(108) = (X(108)-JVS(669)*X(111)-JVS(670)*X(120)-JVS(671)*X(123)-JVS(672)*X(125)-JVS(673)*X(127)-JVS(674)*X(128)&
             &-JVS(675)*X(129)-JVS(676)*X(131)-JVS(677)*X(132)-JVS(678)*X(133)-JVS(679)*X(135)-JVS(680)*X(136))/(JVS(668))
  X(107) = (X(107)-JVS(641)*X(108)-JVS(642)*X(109)-JVS(643)*X(110)-JVS(644)*X(111)-JVS(645)*X(112)-JVS(646)*X(114)&
             &-JVS(647)*X(115)-JVS(648)*X(116)-JVS(649)*X(117)-JVS(650)*X(120)-JVS(651)*X(123)-JVS(652)*X(124)-JVS(653)&
             &*X(125)-JVS(654)*X(126)-JVS(655)*X(127)-JVS(656)*X(128)-JVS(657)*X(129)-JVS(658)*X(131)-JVS(659)*X(132)&
             &-JVS(660)*X(133)-JVS(661)*X(134)-JVS(662)*X(135)-JVS(663)*X(136))/(JVS(640))
  X(106) = (X(106)-JVS(614)*X(132)-JVS(615)*X(135)-JVS(616)*X(136))/(JVS(613))
  X(105) = (X(105)-JVS(602)*X(111)-JVS(603)*X(113)-JVS(604)*X(114)-JVS(605)*X(118)-JVS(606)*X(120)-JVS(607)*X(129)&
             &-JVS(608)*X(131)-JVS(609)*X(132)-JVS(610)*X(133)-JVS(611)*X(135)-JVS(612)*X(136))/(JVS(601))
  X(104) = (X(104)-JVS(582)*X(132)-JVS(583)*X(135)-JVS(584)*X(136))/(JVS(581))
  X(103) = (X(103)-JVS(574)*X(111)-JVS(575)*X(129)-JVS(576)*X(131)-JVS(577)*X(132)-JVS(578)*X(133)-JVS(579)*X(135)&
             &-JVS(580)*X(136))/(JVS(573))
  X(102) = (X(102)-JVS(565)*X(131)-JVS(566)*X(132)-JVS(567)*X(133)-JVS(568)*X(135)-JVS(569)*X(136))/(JVS(564))
  X(101) = (X(101)-JVS(560)*X(132)-JVS(561)*X(135)-JVS(562)*X(136))/(JVS(559))
  X(100) = (X(100)-JVS(556)*X(132)-JVS(557)*X(135)-JVS(558)*X(136))/(JVS(555))
  X(99) = (X(99)-JVS(551)*X(106)-JVS(552)*X(131)-JVS(553)*X(133)-JVS(554)*X(135))/(JVS(550))
  X(98) = (X(98)-JVS(543)*X(100)-JVS(544)*X(101)-JVS(545)*X(104)-JVS(546)*X(131)-JVS(547)*X(132)-JVS(548)*X(133)&
            &-JVS(549)*X(136))/(JVS(542))
  X(97) = (X(97)-JVS(533)*X(113)-JVS(534)*X(129)-JVS(535)*X(131)-JVS(536)*X(133)-JVS(537)*X(136))/(JVS(532))
  X(96) = (X(96)-JVS(517)*X(98)-JVS(518)*X(108)-JVS(519)*X(123)-JVS(520)*X(127)-JVS(521)*X(128)-JVS(522)*X(129)-JVS(523)&
            &*X(131)-JVS(524)*X(133)-JVS(525)*X(134)-JVS(526)*X(135)-JVS(527)*X(136))/(JVS(516))
  X(95) = (X(95)-JVS(496)*X(97)-JVS(497)*X(111)-JVS(498)*X(118)-JVS(499)*X(119)-JVS(500)*X(120)-JVS(501)*X(121)-JVS(502)&
            &*X(123)-JVS(503)*X(125)-JVS(504)*X(127)-JVS(505)*X(129)-JVS(506)*X(136))/(JVS(495))
  X(94) = (X(94)-JVS(491)*X(131)-JVS(492)*X(133)-JVS(493)*X(134)-JVS(494)*X(136))/(JVS(490))
  X(93) = (X(93)-JVS(483)*X(129)-JVS(484)*X(131)-JVS(485)*X(133)-JVS(486)*X(136))/(JVS(482))
  X(92) = (X(92)-JVS(474)*X(109)-JVS(475)*X(126)-JVS(476)*X(131)-JVS(477)*X(132)-JVS(478)*X(134)-JVS(479)*X(136))&
            &/(JVS(473))
  X(91) = (X(91)-JVS(468)*X(131)-JVS(469)*X(133)-JVS(470)*X(134)-JVS(471)*X(136))/(JVS(467))
  X(90) = (X(90)-JVS(461)*X(131)-JVS(462)*X(133)-JVS(463)*X(136))/(JVS(460))
  X(89) = (X(89)-JVS(455)*X(132)-JVS(456)*X(135)-JVS(457)*X(136))/(JVS(454))
  X(88) = (X(88)-JVS(451)*X(131)-JVS(452)*X(133)-JVS(453)*X(136))/(JVS(450))
  X(87) = (X(87)-JVS(442)*X(88)-JVS(443)*X(90)-JVS(444)*X(131)-JVS(445)*X(133)-JVS(446)*X(134)-JVS(447)*X(136))&
            &/(JVS(441))
  X(86) = (X(86)-JVS(435)*X(109)-JVS(436)*X(131)-JVS(437)*X(133)-JVS(438)*X(136))/(JVS(434))
  X(85) = (X(85)-JVS(429)*X(113)-JVS(430)*X(131)-JVS(431)*X(133)-JVS(432)*X(136))/(JVS(428))
  X(84) = (X(84)-JVS(422)*X(113)-JVS(423)*X(131)-JVS(424)*X(136))/(JVS(421))
  X(83) = (X(83)-JVS(419)*X(132)-JVS(420)*X(136))/(JVS(418))
  X(82) = (X(82)-JVS(416)*X(132)-JVS(417)*X(136))/(JVS(415))
  X(81) = (X(81)-JVS(410)*X(131)-JVS(411)*X(132)-JVS(412)*X(133)-JVS(413)*X(134)-JVS(414)*X(136))/(JVS(409))
  X(80) = (X(80)-JVS(402)*X(109)-JVS(403)*X(126)-JVS(404)*X(132)-JVS(405)*X(136))/(JVS(401))
  X(79) = (X(79)-JVS(395)*X(120)-JVS(396)*X(131)-JVS(397)*X(136))/(JVS(394))
  X(78) = (X(78)-JVS(383)*X(82)-JVS(384)*X(83)-JVS(385)*X(89)-JVS(386)*X(100)-JVS(387)*X(101)-JVS(388)*X(104)-JVS(389)&
            &*X(126)-JVS(390)*X(131)-JVS(391)*X(132)-JVS(392)*X(133)-JVS(393)*X(136))/(JVS(382))
  X(77) = (X(77)-JVS(376)*X(130)-JVS(377)*X(134)-JVS(378)*X(136))/(JVS(375))
  X(76) = (X(76)-JVS(372)*X(111)-JVS(373)*X(131)-JVS(374)*X(136))/(JVS(371))
  X(75) = (X(75)-JVS(368)*X(127)-JVS(369)*X(131)-JVS(370)*X(136))/(JVS(367))
  X(74) = (X(74)-JVS(363)*X(131)-JVS(364)*X(133)-JVS(365)*X(134)-JVS(366)*X(136))/(JVS(362))
  X(73) = (X(73)-JVS(355)*X(103)-JVS(356)*X(131)-JVS(357)*X(133)-JVS(358)*X(136))/(JVS(354))
  X(72) = (X(72)-JVS(350)*X(80)-JVS(351)*X(109)-JVS(352)*X(132)-JVS(353)*X(136))/(JVS(349))
  X(71) = (X(71)-JVS(343)*X(112)-JVS(344)*X(117)-JVS(345)*X(126)-JVS(346)*X(134)-JVS(347)*X(135)-JVS(348)*X(136))&
            &/(JVS(342))
  X(70) = (X(70)-JVS(337)*X(86)-JVS(338)*X(131)-JVS(339)*X(136))/(JVS(336))
  X(69) = (X(69)-JVS(333)*X(131)-JVS(334)*X(133)-JVS(335)*X(136))/(JVS(332))
  X(68) = (X(68)-JVS(327)*X(128)-JVS(328)*X(134)-JVS(329)*X(136))/(JVS(326))
  X(67) = (X(67)-JVS(320)*X(102)-JVS(321)*X(109)-JVS(322)*X(113)-JVS(323)*X(133)-JVS(324)*X(135)-JVS(325)*X(136))&
            &/(JVS(319))
  X(66) = (X(66)-JVS(316)*X(131)-JVS(317)*X(133)-JVS(318)*X(136))/(JVS(315))
  X(65) = (X(65)-JVS(310)*X(118)-JVS(311)*X(131)-JVS(312)*X(136))/(JVS(309))
  X(64) = (X(64)-JVS(305)*X(128)-JVS(306)*X(130)-JVS(307)*X(131)-JVS(308)*X(136))/(JVS(304))
  X(63) = (X(63)-JVS(296)*X(106)-JVS(297)*X(109)-JVS(298)*X(128)-JVS(299)*X(129)-JVS(300)*X(130)-JVS(301)*X(131)&
            &-JVS(302)*X(132)-JVS(303)*X(136))/(JVS(295))
  X(62) = (X(62)-JVS(291)*X(84)-JVS(292)*X(113)-JVS(293)*X(133)-JVS(294)*X(136))/(JVS(290))
  X(61) = (X(61)-JVS(287)*X(131)-JVS(288)*X(133)-JVS(289)*X(136))/(JVS(286))
  X(60) = (X(60)-JVS(280)*X(132)-JVS(281)*X(133)-JVS(282)*X(136))/(JVS(279))
  X(59) = (X(59)-JVS(273)*X(81)-JVS(274)*X(132)-JVS(275)*X(133)-JVS(276)*X(134)-JVS(277)*X(136))/(JVS(272))
  X(58) = (X(58)-JVS(266)*X(131)-JVS(267)*X(133)-JVS(268)*X(136))/(JVS(265))
  X(57) = (X(57)-JVS(259)*X(88)-JVS(260)*X(131)-JVS(261)*X(136))/(JVS(258))
  X(56) = (X(56)-JVS(255)*X(126)-JVS(256)*X(131)-JVS(257)*X(133))/(JVS(254))
  X(55) = (X(55)-JVS(252)*X(132)-JVS(253)*X(136))/(JVS(251))
  X(54) = (X(54)-JVS(244)*X(82)-JVS(245)*X(83)-JVS(246)*X(100)-JVS(247)*X(101)-JVS(248)*X(104)-JVS(249)*X(132)-JVS(250)&
            &*X(136))/(JVS(243))
  X(53) = (X(53)-JVS(236)*X(116)-JVS(237)*X(127)-JVS(238)*X(128)-JVS(239)*X(129)-JVS(240)*X(133)-JVS(241)*X(135)&
            &-JVS(242)*X(136))/(JVS(235))
  X(52) = (X(52)-JVS(232)*X(131)-JVS(233)*X(134)-JVS(234)*X(136))/(JVS(231))
  X(51) = (X(51)-JVS(228)*X(125)-JVS(229)*X(131)-JVS(230)*X(136))/(JVS(227))
  X(50) = (X(50)-JVS(224)*X(129)-JVS(225)*X(131)-JVS(226)*X(136))/(JVS(223))
  X(49) = (X(49)-JVS(220)*X(85)-JVS(221)*X(131)-JVS(222)*X(136))/(JVS(219))
  X(48) = (X(48)-JVS(217)*X(133)-JVS(218)*X(136))/(JVS(216))
  X(47) = (X(47)-JVS(212)*X(121)-JVS(213)*X(131)-JVS(214)*X(136))/(JVS(211))
  X(46) = (X(46)-JVS(208)*X(93)-JVS(209)*X(131)-JVS(210)*X(136))/(JVS(207))
  X(45) = (X(45)-JVS(204)*X(97)-JVS(205)*X(131)-JVS(206)*X(136))/(JVS(203))
  X(44) = (X(44)-JVS(199)*X(72)-JVS(200)*X(80)-JVS(201)*X(126)-JVS(202)*X(136))/(JVS(198))
  X(43) = (X(43)-JVS(193)*X(66)-JVS(194)*X(88)-JVS(195)*X(90)-JVS(196)*X(133)-JVS(197)*X(136))/(JVS(192))
  X(42) = (X(42)-JVS(188)*X(88)-JVS(189)*X(90)-JVS(190)*X(133))/(JVS(187))
  X(41) = (X(41)-JVS(184)*X(97)-JVS(185)*X(129)-JVS(186)*X(136))/(JVS(183))
  X(40) = (X(40)-JVS(181)*X(90)-JVS(182)*X(133))/(JVS(180))
  X(39) = (X(39)-JVS(178)*X(135)-JVS(179)*X(136))/(JVS(177))
  X(38) = (X(38)-JVS(176)*X(136))/(JVS(175))
  X(37) = (X(37)-JVS(173)*X(136))/(JVS(172))
  X(36) = (X(36)-JVS(170)*X(134)-JVS(171)*X(135))/(JVS(169))
  X(35) = (X(35)-JVS(166)*X(90)-JVS(167)*X(131)-JVS(168)*X(136))/(JVS(165))
  X(34) = (X(34)-JVS(162)*X(69)-JVS(163)*X(131)-JVS(164)*X(136))/(JVS(161))
  X(33) = (X(33)-JVS(158)*X(61)-JVS(159)*X(133)-JVS(160)*X(136))/(JVS(157))
  X(32) = (X(32)-JVS(154)*X(61)-JVS(155)*X(131)-JVS(156)*X(136))/(JVS(153))
  X(31) = (X(31)-JVS(151)*X(37)-JVS(152)*X(136))/(JVS(150))
  X(30) = (X(30)-JVS(146)*X(88)-JVS(147)*X(90)-JVS(148)*X(133))/(JVS(145))
  X(29) = (X(29)-JVS(142)*X(66)-JVS(143)*X(131)-JVS(144)*X(136))/(JVS(141))
  X(28) = (X(28)-JVS(138)*X(81)-JVS(139)*X(131)-JVS(140)*X(136))/(JVS(137))
  X(27) = (X(27)-JVS(134)*X(58)-JVS(135)*X(131)-JVS(136)*X(136))/(JVS(133))
  X(26) = (X(26)-JVS(130)*X(123)-JVS(131)*X(131)-JVS(132)*X(136))/(JVS(129))
  X(25) = (X(25)-JVS(128)*X(136))/(JVS(127))
  X(24) = (X(24)-JVS(125)*X(136))/(JVS(124))
  X(23) = (X(23)-JVS(122)*X(60)-JVS(123)*X(133))/(JVS(121))
  X(22) = (X(22)-JVS(120)*X(136))/(JVS(119))
  X(21) = (X(21)-JVS(117)*X(136))/(JVS(116))
  X(20) = (X(20)-JVS(114)*X(136))/(JVS(113))
  X(19) = (X(19)-JVS(112)*X(136))/(JVS(111))
  X(18) = (X(18)-JVS(109)*X(131)-JVS(110)*X(136))/(JVS(108))
  X(17) = (X(17)-JVS(106)*X(74)-JVS(107)*X(134))/(JVS(105))
  X(16) = (X(16)-JVS(102)*X(39)-JVS(103)*X(135)-JVS(104)*X(136))/(JVS(101))
  X(15) = (X(15)-JVS(100)*X(80))/(JVS(99))
  X(14) = (X(14)-JVS(98)*X(136))/(JVS(97))
  X(13) = (X(13)-JVS(96)*X(136))/(JVS(95))
  X(12) = (X(12)-JVS(94)*X(136))/(JVS(93))
  X(11) = (X(11)-JVS(92)*X(136))/(JVS(91))
  X(10) = (X(10)-JVS(90)*X(136))/(JVS(89))
  X(9) = (X(9)-JVS(87)*X(10)-JVS(88)*X(136))/(JVS(86))
  X(8) = (X(8)-JVS(85)*X(136))/(JVS(84))
  X(7) = (X(7)-JVS(82)*X(8)-JVS(83)*X(136))/(JVS(81))
  X(6) = (X(6)-JVS(77)*X(82)-JVS(78)*X(83)-JVS(79)*X(101)-JVS(80)*X(136))/(JVS(76))
  X(5) = (X(5)-JVS(74)*X(106)-JVS(75)*X(136))/(JVS(73))
  X(4) = (X(4)-JVS(41)*X(56)-JVS(42)*X(58)-JVS(43)*X(61)-JVS(44)*X(66)-JVS(45)*X(69)-JVS(46)*X(73)-JVS(47)*X(74)-JVS(48)&
           &*X(81)-JVS(49)*X(85)-JVS(50)*X(86)-JVS(51)*X(87)-JVS(52)*X(88)-JVS(53)*X(90)-JVS(54)*X(91)-JVS(55)*X(93)-JVS(56)&
           &*X(94)-JVS(57)*X(97)-JVS(58)*X(99)-JVS(59)*X(102)-JVS(60)*X(111)-JVS(61)*X(113)-JVS(62)*X(118)-JVS(63)*X(119)&
           &-JVS(64)*X(120)-JVS(65)*X(121)-JVS(66)*X(123)-JVS(67)*X(125)-JVS(68)*X(127)-JVS(69)*X(128)-JVS(70)*X(129)&
           &-JVS(71)*X(130)-JVS(72)*X(131))/(JVS(40))
  X(3) = (X(3)-JVS(6)*X(48)-JVS(7)*X(56)-JVS(8)*X(58)-JVS(9)*X(60)-JVS(10)*X(61)-JVS(11)*X(66)-JVS(12)*X(69)-JVS(13)&
           &*X(73)-JVS(14)*X(74)-JVS(15)*X(81)-JVS(16)*X(85)-JVS(17)*X(86)-JVS(18)*X(87)-JVS(19)*X(88)-JVS(20)*X(90)-JVS(21)&
           &*X(91)-JVS(22)*X(93)-JVS(23)*X(94)-JVS(24)*X(97)-JVS(25)*X(99)-JVS(26)*X(102)-JVS(27)*X(111)-JVS(28)*X(113)&
           &-JVS(29)*X(118)-JVS(30)*X(119)-JVS(31)*X(120)-JVS(32)*X(121)-JVS(33)*X(123)-JVS(34)*X(125)-JVS(35)*X(127)&
           &-JVS(36)*X(128)-JVS(37)*X(129)-JVS(38)*X(130)-JVS(39)*X(133))/(JVS(5))
  X(2) = (X(2)-JVS(3)*X(16)-JVS(4)*X(136))/(JVS(2))
  X(1) = X(1)/JVS(1)
      
END SUBROUTINE mozart_mosaic_4bin_KppSolve
























      SUBROUTINE mozart_mosaic_4bin_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE mozart_mosaic_4bin_WCOPY



      SUBROUTINE mozart_mosaic_4bin_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE mozart_mosaic_4bin_WAXPY




      SUBROUTINE mozart_mosaic_4bin_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE mozart_mosaic_4bin_WSCAL


      REAL(kind=dp) FUNCTION mozart_mosaic_4bin_WLAMCH( C )








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
          CALL mozart_mosaic_4bin_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      mozart_mosaic_4bin_WLAMCH = Eps

      END FUNCTION mozart_mosaic_4bin_WLAMCH
     
      SUBROUTINE mozart_mosaic_4bin_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE mozart_mosaic_4bin_WLAMCH_ADD




      SUBROUTINE mozart_mosaic_4bin_SET2ZERO(N,Y)




      
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

      END SUBROUTINE mozart_mosaic_4bin_SET2ZERO



      REAL(kind=dp) FUNCTION mozart_mosaic_4bin_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      mozart_mosaic_4bin_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        mozart_mosaic_4bin_WDOT = mozart_mosaic_4bin_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         mozart_mosaic_4bin_WDOT = mozart_mosaic_4bin_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          mozart_mosaic_4bin_WDOT = mozart_mosaic_4bin_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        mozart_mosaic_4bin_WDOT = mozart_mosaic_4bin_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION mozart_mosaic_4bin_WDOT                                          




   SUBROUTINE decomp_mozart_mosaic_4bin( JVS, IER )
   
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
  JVS( 1) = W( 1 )
  IF ( ABS(  JVS( 2 )) < TINY(a) ) THEN
         IER = 2                                       
         RETURN
  END IF
   W( 2 ) = JVS( 2 )
   W( 16 ) = JVS( 3 )
   W( 136 ) = JVS( 4 )
  JVS( 2) = W( 2 )
  JVS( 3) = W( 16 )
  JVS( 4) = W( 136 )
  IF ( ABS(  JVS( 5 )) < TINY(a) ) THEN
         IER = 3                                       
         RETURN
  END IF
   W( 3 ) = JVS( 5 )
   W( 48 ) = JVS( 6 )
   W( 56 ) = JVS( 7 )
   W( 58 ) = JVS( 8 )
   W( 60 ) = JVS( 9 )
   W( 61 ) = JVS( 10 )
   W( 66 ) = JVS( 11 )
   W( 69 ) = JVS( 12 )
   W( 73 ) = JVS( 13 )
   W( 74 ) = JVS( 14 )
   W( 81 ) = JVS( 15 )
   W( 85 ) = JVS( 16 )
   W( 86 ) = JVS( 17 )
   W( 87 ) = JVS( 18 )
   W( 88 ) = JVS( 19 )
   W( 90 ) = JVS( 20 )
   W( 91 ) = JVS( 21 )
   W( 93 ) = JVS( 22 )
   W( 94 ) = JVS( 23 )
   W( 97 ) = JVS( 24 )
   W( 99 ) = JVS( 25 )
   W( 102 ) = JVS( 26 )
   W( 111 ) = JVS( 27 )
   W( 113 ) = JVS( 28 )
   W( 118 ) = JVS( 29 )
   W( 119 ) = JVS( 30 )
   W( 120 ) = JVS( 31 )
   W( 121 ) = JVS( 32 )
   W( 123 ) = JVS( 33 )
   W( 125 ) = JVS( 34 )
   W( 127 ) = JVS( 35 )
   W( 128 ) = JVS( 36 )
   W( 129 ) = JVS( 37 )
   W( 130 ) = JVS( 38 )
   W( 133 ) = JVS( 39 )
  JVS( 5) = W( 3 )
  JVS( 6) = W( 48 )
  JVS( 7) = W( 56 )
  JVS( 8) = W( 58 )
  JVS( 9) = W( 60 )
  JVS( 10) = W( 61 )
  JVS( 11) = W( 66 )
  JVS( 12) = W( 69 )
  JVS( 13) = W( 73 )
  JVS( 14) = W( 74 )
  JVS( 15) = W( 81 )
  JVS( 16) = W( 85 )
  JVS( 17) = W( 86 )
  JVS( 18) = W( 87 )
  JVS( 19) = W( 88 )
  JVS( 20) = W( 90 )
  JVS( 21) = W( 91 )
  JVS( 22) = W( 93 )
  JVS( 23) = W( 94 )
  JVS( 24) = W( 97 )
  JVS( 25) = W( 99 )
  JVS( 26) = W( 102 )
  JVS( 27) = W( 111 )
  JVS( 28) = W( 113 )
  JVS( 29) = W( 118 )
  JVS( 30) = W( 119 )
  JVS( 31) = W( 120 )
  JVS( 32) = W( 121 )
  JVS( 33) = W( 123 )
  JVS( 34) = W( 125 )
  JVS( 35) = W( 127 )
  JVS( 36) = W( 128 )
  JVS( 37) = W( 129 )
  JVS( 38) = W( 130 )
  JVS( 39) = W( 133 )
  IF ( ABS(  JVS( 40 )) < TINY(a) ) THEN
         IER = 4                                       
         RETURN
  END IF
   W( 4 ) = JVS( 40 )
   W( 56 ) = JVS( 41 )
   W( 58 ) = JVS( 42 )
   W( 61 ) = JVS( 43 )
   W( 66 ) = JVS( 44 )
   W( 69 ) = JVS( 45 )
   W( 73 ) = JVS( 46 )
   W( 74 ) = JVS( 47 )
   W( 81 ) = JVS( 48 )
   W( 85 ) = JVS( 49 )
   W( 86 ) = JVS( 50 )
   W( 87 ) = JVS( 51 )
   W( 88 ) = JVS( 52 )
   W( 90 ) = JVS( 53 )
   W( 91 ) = JVS( 54 )
   W( 93 ) = JVS( 55 )
   W( 94 ) = JVS( 56 )
   W( 97 ) = JVS( 57 )
   W( 99 ) = JVS( 58 )
   W( 102 ) = JVS( 59 )
   W( 111 ) = JVS( 60 )
   W( 113 ) = JVS( 61 )
   W( 118 ) = JVS( 62 )
   W( 119 ) = JVS( 63 )
   W( 120 ) = JVS( 64 )
   W( 121 ) = JVS( 65 )
   W( 123 ) = JVS( 66 )
   W( 125 ) = JVS( 67 )
   W( 127 ) = JVS( 68 )
   W( 128 ) = JVS( 69 )
   W( 129 ) = JVS( 70 )
   W( 130 ) = JVS( 71 )
   W( 131 ) = JVS( 72 )
  JVS( 40) = W( 4 )
  JVS( 41) = W( 56 )
  JVS( 42) = W( 58 )
  JVS( 43) = W( 61 )
  JVS( 44) = W( 66 )
  JVS( 45) = W( 69 )
  JVS( 46) = W( 73 )
  JVS( 47) = W( 74 )
  JVS( 48) = W( 81 )
  JVS( 49) = W( 85 )
  JVS( 50) = W( 86 )
  JVS( 51) = W( 87 )
  JVS( 52) = W( 88 )
  JVS( 53) = W( 90 )
  JVS( 54) = W( 91 )
  JVS( 55) = W( 93 )
  JVS( 56) = W( 94 )
  JVS( 57) = W( 97 )
  JVS( 58) = W( 99 )
  JVS( 59) = W( 102 )
  JVS( 60) = W( 111 )
  JVS( 61) = W( 113 )
  JVS( 62) = W( 118 )
  JVS( 63) = W( 119 )
  JVS( 64) = W( 120 )
  JVS( 65) = W( 121 )
  JVS( 66) = W( 123 )
  JVS( 67) = W( 125 )
  JVS( 68) = W( 127 )
  JVS( 69) = W( 128 )
  JVS( 70) = W( 129 )
  JVS( 71) = W( 130 )
  JVS( 72) = W( 131 )
  IF ( ABS(  JVS( 73 )) < TINY(a) ) THEN
         IER = 5                                       
         RETURN
  END IF
   W( 5 ) = JVS( 73 )
   W( 106 ) = JVS( 74 )
   W( 136 ) = JVS( 75 )
  JVS( 73) = W( 5 )
  JVS( 74) = W( 106 )
  JVS( 75) = W( 136 )
  IF ( ABS(  JVS( 76 )) < TINY(a) ) THEN
         IER = 6                                       
         RETURN
  END IF
   W( 6 ) = JVS( 76 )
   W( 82 ) = JVS( 77 )
   W( 83 ) = JVS( 78 )
   W( 101 ) = JVS( 79 )
   W( 136 ) = JVS( 80 )
  JVS( 76) = W( 6 )
  JVS( 77) = W( 82 )
  JVS( 78) = W( 83 )
  JVS( 79) = W( 101 )
  JVS( 80) = W( 136 )
  IF ( ABS(  JVS( 81 )) < TINY(a) ) THEN
         IER = 7                                       
         RETURN
  END IF
   W( 7 ) = JVS( 81 )
   W( 8 ) = JVS( 82 )
   W( 136 ) = JVS( 83 )
  JVS( 81) = W( 7 )
  JVS( 82) = W( 8 )
  JVS( 83) = W( 136 )
  IF ( ABS(  JVS( 84 )) < TINY(a) ) THEN
         IER = 8                                       
         RETURN
  END IF
   W( 8 ) = JVS( 84 )
   W( 136 ) = JVS( 85 )
  JVS( 84) = W( 8 )
  JVS( 85) = W( 136 )
  IF ( ABS(  JVS( 86 )) < TINY(a) ) THEN
         IER = 9                                       
         RETURN
  END IF
   W( 9 ) = JVS( 86 )
   W( 10 ) = JVS( 87 )
   W( 136 ) = JVS( 88 )
  JVS( 86) = W( 9 )
  JVS( 87) = W( 10 )
  JVS( 88) = W( 136 )
  IF ( ABS(  JVS( 89 )) < TINY(a) ) THEN
         IER = 10                                      
         RETURN
  END IF
   W( 10 ) = JVS( 89 )
   W( 136 ) = JVS( 90 )
  JVS( 89) = W( 10 )
  JVS( 90) = W( 136 )
  IF ( ABS(  JVS( 91 )) < TINY(a) ) THEN
         IER = 11                                      
         RETURN
  END IF
   W( 11 ) = JVS( 91 )
   W( 136 ) = JVS( 92 )
  JVS( 91) = W( 11 )
  JVS( 92) = W( 136 )
  IF ( ABS(  JVS( 93 )) < TINY(a) ) THEN
         IER = 12                                      
         RETURN
  END IF
   W( 12 ) = JVS( 93 )
   W( 136 ) = JVS( 94 )
  JVS( 93) = W( 12 )
  JVS( 94) = W( 136 )
  IF ( ABS(  JVS( 95 )) < TINY(a) ) THEN
         IER = 13                                      
         RETURN
  END IF
   W( 13 ) = JVS( 95 )
   W( 136 ) = JVS( 96 )
  JVS( 95) = W( 13 )
  JVS( 96) = W( 136 )
  IF ( ABS(  JVS( 97 )) < TINY(a) ) THEN
         IER = 14                                      
         RETURN
  END IF
   W( 14 ) = JVS( 97 )
   W( 136 ) = JVS( 98 )
  JVS( 97) = W( 14 )
  JVS( 98) = W( 136 )
  IF ( ABS(  JVS( 99 )) < TINY(a) ) THEN
         IER = 15                                      
         RETURN
  END IF
   W( 15 ) = JVS( 99 )
   W( 80 ) = JVS( 100 )
  JVS( 99) = W( 15 )
  JVS( 100) = W( 80 )
  IF ( ABS(  JVS( 101 )) < TINY(a) ) THEN
         IER = 16                                      
         RETURN
  END IF
   W( 16 ) = JVS( 101 )
   W( 39 ) = JVS( 102 )
   W( 135 ) = JVS( 103 )
   W( 136 ) = JVS( 104 )
  JVS( 101) = W( 16 )
  JVS( 102) = W( 39 )
  JVS( 103) = W( 135 )
  JVS( 104) = W( 136 )
  IF ( ABS(  JVS( 105 )) < TINY(a) ) THEN
         IER = 17                                      
         RETURN
  END IF
   W( 17 ) = JVS( 105 )
   W( 74 ) = JVS( 106 )
   W( 134 ) = JVS( 107 )
  JVS( 105) = W( 17 )
  JVS( 106) = W( 74 )
  JVS( 107) = W( 134 )
  IF ( ABS(  JVS( 108 )) < TINY(a) ) THEN
         IER = 18                                      
         RETURN
  END IF
   W( 18 ) = JVS( 108 )
   W( 131 ) = JVS( 109 )
   W( 136 ) = JVS( 110 )
  JVS( 108) = W( 18 )
  JVS( 109) = W( 131 )
  JVS( 110) = W( 136 )
  IF ( ABS(  JVS( 111 )) < TINY(a) ) THEN
         IER = 19                                      
         RETURN
  END IF
   W( 19 ) = JVS( 111 )
   W( 136 ) = JVS( 112 )
  JVS( 111) = W( 19 )
  JVS( 112) = W( 136 )
  IF ( ABS(  JVS( 113 )) < TINY(a) ) THEN
         IER = 20                                      
         RETURN
  END IF
   W( 20 ) = JVS( 113 )
   W( 136 ) = JVS( 114 )
  JVS( 113) = W( 20 )
  JVS( 114) = W( 136 )
  IF ( ABS(  JVS( 116 )) < TINY(a) ) THEN
         IER = 21                                      
         RETURN
  END IF
   W( 20 ) = JVS( 115 )
   W( 21 ) = JVS( 116 )
   W( 136 ) = JVS( 117 )
  a = -W( 20 ) / JVS(          113  )
  W( 20 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 114 )
  JVS( 115) = W( 20 )
  JVS( 116) = W( 21 )
  JVS( 117) = W( 136 )
  IF ( ABS(  JVS( 119 )) < TINY(a) ) THEN
         IER = 22                                      
         RETURN
  END IF
   W( 20 ) = JVS( 118 )
   W( 22 ) = JVS( 119 )
   W( 136 ) = JVS( 120 )
  a = -W( 20 ) / JVS(          113  )
  W( 20 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 114 )
  JVS( 118) = W( 20 )
  JVS( 119) = W( 22 )
  JVS( 120) = W( 136 )
  IF ( ABS(  JVS( 121 )) < TINY(a) ) THEN
         IER = 23                                      
         RETURN
  END IF
   W( 23 ) = JVS( 121 )
   W( 60 ) = JVS( 122 )
   W( 133 ) = JVS( 123 )
  JVS( 121) = W( 23 )
  JVS( 122) = W( 60 )
  JVS( 123) = W( 133 )
  IF ( ABS(  JVS( 124 )) < TINY(a) ) THEN
         IER = 24                                      
         RETURN
  END IF
   W( 24 ) = JVS( 124 )
   W( 136 ) = JVS( 125 )
  JVS( 124) = W( 24 )
  JVS( 125) = W( 136 )
  IF ( ABS(  JVS( 127 )) < TINY(a) ) THEN
         IER = 25                                      
         RETURN
  END IF
   W( 24 ) = JVS( 126 )
   W( 25 ) = JVS( 127 )
   W( 136 ) = JVS( 128 )
  a = -W( 24 ) / JVS(          124  )
  W( 24 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 125 )
  JVS( 126) = W( 24 )
  JVS( 127) = W( 25 )
  JVS( 128) = W( 136 )
  IF ( ABS(  JVS( 129 )) < TINY(a) ) THEN
         IER = 26                                      
         RETURN
  END IF
   W( 26 ) = JVS( 129 )
   W( 123 ) = JVS( 130 )
   W( 131 ) = JVS( 131 )
   W( 136 ) = JVS( 132 )
  JVS( 129) = W( 26 )
  JVS( 130) = W( 123 )
  JVS( 131) = W( 131 )
  JVS( 132) = W( 136 )
  IF ( ABS(  JVS( 133 )) < TINY(a) ) THEN
         IER = 27                                      
         RETURN
  END IF
   W( 27 ) = JVS( 133 )
   W( 58 ) = JVS( 134 )
   W( 131 ) = JVS( 135 )
   W( 136 ) = JVS( 136 )
  JVS( 133) = W( 27 )
  JVS( 134) = W( 58 )
  JVS( 135) = W( 131 )
  JVS( 136) = W( 136 )
  IF ( ABS(  JVS( 137 )) < TINY(a) ) THEN
         IER = 28                                      
         RETURN
  END IF
   W( 28 ) = JVS( 137 )
   W( 81 ) = JVS( 138 )
   W( 131 ) = JVS( 139 )
   W( 136 ) = JVS( 140 )
  JVS( 137) = W( 28 )
  JVS( 138) = W( 81 )
  JVS( 139) = W( 131 )
  JVS( 140) = W( 136 )
  IF ( ABS(  JVS( 141 )) < TINY(a) ) THEN
         IER = 29                                      
         RETURN
  END IF
   W( 29 ) = JVS( 141 )
   W( 66 ) = JVS( 142 )
   W( 131 ) = JVS( 143 )
   W( 136 ) = JVS( 144 )
  JVS( 141) = W( 29 )
  JVS( 142) = W( 66 )
  JVS( 143) = W( 131 )
  JVS( 144) = W( 136 )
  IF ( ABS(  JVS( 145 )) < TINY(a) ) THEN
         IER = 30                                      
         RETURN
  END IF
   W( 30 ) = JVS( 145 )
   W( 88 ) = JVS( 146 )
   W( 90 ) = JVS( 147 )
   W( 133 ) = JVS( 148 )
  JVS( 145) = W( 30 )
  JVS( 146) = W( 88 )
  JVS( 147) = W( 90 )
  JVS( 148) = W( 133 )
  IF ( ABS(  JVS( 150 )) < TINY(a) ) THEN
         IER = 31                                      
         RETURN
  END IF
   W( 24 ) = JVS( 149 )
   W( 31 ) = JVS( 150 )
   W( 37 ) = JVS( 151 )
   W( 136 ) = JVS( 152 )
  a = -W( 24 ) / JVS(          124  )
  W( 24 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 125 )
  JVS( 149) = W( 24 )
  JVS( 150) = W( 31 )
  JVS( 151) = W( 37 )
  JVS( 152) = W( 136 )
  IF ( ABS(  JVS( 153 )) < TINY(a) ) THEN
         IER = 32                                      
         RETURN
  END IF
   W( 32 ) = JVS( 153 )
   W( 61 ) = JVS( 154 )
   W( 131 ) = JVS( 155 )
   W( 136 ) = JVS( 156 )
  JVS( 153) = W( 32 )
  JVS( 154) = W( 61 )
  JVS( 155) = W( 131 )
  JVS( 156) = W( 136 )
  IF ( ABS(  JVS( 157 )) < TINY(a) ) THEN
         IER = 33                                      
         RETURN
  END IF
   W( 33 ) = JVS( 157 )
   W( 61 ) = JVS( 158 )
   W( 133 ) = JVS( 159 )
   W( 136 ) = JVS( 160 )
  JVS( 157) = W( 33 )
  JVS( 158) = W( 61 )
  JVS( 159) = W( 133 )
  JVS( 160) = W( 136 )
  IF ( ABS(  JVS( 161 )) < TINY(a) ) THEN
         IER = 34                                      
         RETURN
  END IF
   W( 34 ) = JVS( 161 )
   W( 69 ) = JVS( 162 )
   W( 131 ) = JVS( 163 )
   W( 136 ) = JVS( 164 )
  JVS( 161) = W( 34 )
  JVS( 162) = W( 69 )
  JVS( 163) = W( 131 )
  JVS( 164) = W( 136 )
  IF ( ABS(  JVS( 165 )) < TINY(a) ) THEN
         IER = 35                                      
         RETURN
  END IF
   W( 35 ) = JVS( 165 )
   W( 90 ) = JVS( 166 )
   W( 131 ) = JVS( 167 )
   W( 136 ) = JVS( 168 )
  JVS( 165) = W( 35 )
  JVS( 166) = W( 90 )
  JVS( 167) = W( 131 )
  JVS( 168) = W( 136 )
  IF ( ABS(  JVS( 169 )) < TINY(a) ) THEN
         IER = 36                                      
         RETURN
  END IF
   W( 36 ) = JVS( 169 )
   W( 134 ) = JVS( 170 )
   W( 135 ) = JVS( 171 )
  JVS( 169) = W( 36 )
  JVS( 170) = W( 134 )
  JVS( 171) = W( 135 )
  IF ( ABS(  JVS( 172 )) < TINY(a) ) THEN
         IER = 37                                      
         RETURN
  END IF
   W( 37 ) = JVS( 172 )
   W( 136 ) = JVS( 173 )
  JVS( 172) = W( 37 )
  JVS( 173) = W( 136 )
  IF ( ABS(  JVS( 175 )) < TINY(a) ) THEN
         IER = 38                                      
         RETURN
  END IF
   W( 37 ) = JVS( 174 )
   W( 38 ) = JVS( 175 )
   W( 136 ) = JVS( 176 )
  a = -W( 37 ) / JVS(          172  )
  W( 37 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 173 )
  JVS( 174) = W( 37 )
  JVS( 175) = W( 38 )
  JVS( 176) = W( 136 )
  IF ( ABS(  JVS( 177 )) < TINY(a) ) THEN
         IER = 39                                      
         RETURN
  END IF
   W( 39 ) = JVS( 177 )
   W( 135 ) = JVS( 178 )
   W( 136 ) = JVS( 179 )
  JVS( 177) = W( 39 )
  JVS( 178) = W( 135 )
  JVS( 179) = W( 136 )
  IF ( ABS(  JVS( 180 )) < TINY(a) ) THEN
         IER = 40                                      
         RETURN
  END IF
   W( 40 ) = JVS( 180 )
   W( 90 ) = JVS( 181 )
   W( 133 ) = JVS( 182 )
  JVS( 180) = W( 40 )
  JVS( 181) = W( 90 )
  JVS( 182) = W( 133 )
  IF ( ABS(  JVS( 183 )) < TINY(a) ) THEN
         IER = 41                                      
         RETURN
  END IF
   W( 41 ) = JVS( 183 )
   W( 97 ) = JVS( 184 )
   W( 129 ) = JVS( 185 )
   W( 136 ) = JVS( 186 )
  JVS( 183) = W( 41 )
  JVS( 184) = W( 97 )
  JVS( 185) = W( 129 )
  JVS( 186) = W( 136 )
  IF ( ABS(  JVS( 187 )) < TINY(a) ) THEN
         IER = 42                                      
         RETURN
  END IF
   W( 42 ) = JVS( 187 )
   W( 88 ) = JVS( 188 )
   W( 90 ) = JVS( 189 )
   W( 133 ) = JVS( 190 )
  JVS( 187) = W( 42 )
  JVS( 188) = W( 88 )
  JVS( 189) = W( 90 )
  JVS( 190) = W( 133 )
  IF ( ABS(  JVS( 192 )) < TINY(a) ) THEN
         IER = 43                                      
         RETURN
  END IF
   W( 21 ) = JVS( 191 )
   W( 43 ) = JVS( 192 )
   W( 66 ) = JVS( 193 )
   W( 88 ) = JVS( 194 )
   W( 90 ) = JVS( 195 )
   W( 133 ) = JVS( 196 )
   W( 136 ) = JVS( 197 )
  a = -W( 21 ) / JVS(          116  )
  W( 21 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 117 )
  JVS( 191) = W( 21 )
  JVS( 192) = W( 43 )
  JVS( 193) = W( 66 )
  JVS( 194) = W( 88 )
  JVS( 195) = W( 90 )
  JVS( 196) = W( 133 )
  JVS( 197) = W( 136 )
  IF ( ABS(  JVS( 198 )) < TINY(a) ) THEN
         IER = 44                                      
         RETURN
  END IF
   W( 44 ) = JVS( 198 )
   W( 72 ) = JVS( 199 )
   W( 80 ) = JVS( 200 )
   W( 126 ) = JVS( 201 )
   W( 136 ) = JVS( 202 )
  JVS( 198) = W( 44 )
  JVS( 199) = W( 72 )
  JVS( 200) = W( 80 )
  JVS( 201) = W( 126 )
  JVS( 202) = W( 136 )
  IF ( ABS(  JVS( 203 )) < TINY(a) ) THEN
         IER = 45                                      
         RETURN
  END IF
   W( 45 ) = JVS( 203 )
   W( 97 ) = JVS( 204 )
   W( 131 ) = JVS( 205 )
   W( 136 ) = JVS( 206 )
  JVS( 203) = W( 45 )
  JVS( 204) = W( 97 )
  JVS( 205) = W( 131 )
  JVS( 206) = W( 136 )
  IF ( ABS(  JVS( 207 )) < TINY(a) ) THEN
         IER = 46                                      
         RETURN
  END IF
   W( 46 ) = JVS( 207 )
   W( 93 ) = JVS( 208 )
   W( 131 ) = JVS( 209 )
   W( 136 ) = JVS( 210 )
  JVS( 207) = W( 46 )
  JVS( 208) = W( 93 )
  JVS( 209) = W( 131 )
  JVS( 210) = W( 136 )
  IF ( ABS(  JVS( 211 )) < TINY(a) ) THEN
         IER = 47                                      
         RETURN
  END IF
   W( 47 ) = JVS( 211 )
   W( 121 ) = JVS( 212 )
   W( 131 ) = JVS( 213 )
   W( 136 ) = JVS( 214 )
  JVS( 211) = W( 47 )
  JVS( 212) = W( 121 )
  JVS( 213) = W( 131 )
  JVS( 214) = W( 136 )
  IF ( ABS(  JVS( 216 )) < TINY(a) ) THEN
         IER = 48                                      
         RETURN
  END IF
   W( 14 ) = JVS( 215 )
   W( 48 ) = JVS( 216 )
   W( 133 ) = JVS( 217 )
   W( 136 ) = JVS( 218 )
  a = -W( 14 ) / JVS(           97  )
  W( 14 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 98 )
  JVS( 215) = W( 14 )
  JVS( 216) = W( 48 )
  JVS( 217) = W( 133 )
  JVS( 218) = W( 136 )
  IF ( ABS(  JVS( 219 )) < TINY(a) ) THEN
         IER = 49                                      
         RETURN
  END IF
   W( 49 ) = JVS( 219 )
   W( 85 ) = JVS( 220 )
   W( 131 ) = JVS( 221 )
   W( 136 ) = JVS( 222 )
  JVS( 219) = W( 49 )
  JVS( 220) = W( 85 )
  JVS( 221) = W( 131 )
  JVS( 222) = W( 136 )
  IF ( ABS(  JVS( 223 )) < TINY(a) ) THEN
         IER = 50                                      
         RETURN
  END IF
   W( 50 ) = JVS( 223 )
   W( 129 ) = JVS( 224 )
   W( 131 ) = JVS( 225 )
   W( 136 ) = JVS( 226 )
  JVS( 223) = W( 50 )
  JVS( 224) = W( 129 )
  JVS( 225) = W( 131 )
  JVS( 226) = W( 136 )
  IF ( ABS(  JVS( 227 )) < TINY(a) ) THEN
         IER = 51                                      
         RETURN
  END IF
   W( 51 ) = JVS( 227 )
   W( 125 ) = JVS( 228 )
   W( 131 ) = JVS( 229 )
   W( 136 ) = JVS( 230 )
  JVS( 227) = W( 51 )
  JVS( 228) = W( 125 )
  JVS( 229) = W( 131 )
  JVS( 230) = W( 136 )
  IF ( ABS(  JVS( 231 )) < TINY(a) ) THEN
         IER = 52                                      
         RETURN
  END IF
   W( 52 ) = JVS( 231 )
   W( 131 ) = JVS( 232 )
   W( 134 ) = JVS( 233 )
   W( 136 ) = JVS( 234 )
  JVS( 231) = W( 52 )
  JVS( 232) = W( 131 )
  JVS( 233) = W( 134 )
  JVS( 234) = W( 136 )
  IF ( ABS(  JVS( 235 )) < TINY(a) ) THEN
         IER = 53                                      
         RETURN
  END IF
   W( 53 ) = JVS( 235 )
   W( 116 ) = JVS( 236 )
   W( 127 ) = JVS( 237 )
   W( 128 ) = JVS( 238 )
   W( 129 ) = JVS( 239 )
   W( 133 ) = JVS( 240 )
   W( 135 ) = JVS( 241 )
   W( 136 ) = JVS( 242 )
  JVS( 235) = W( 53 )
  JVS( 236) = W( 116 )
  JVS( 237) = W( 127 )
  JVS( 238) = W( 128 )
  JVS( 239) = W( 129 )
  JVS( 240) = W( 133 )
  JVS( 241) = W( 135 )
  JVS( 242) = W( 136 )
  IF ( ABS(  JVS( 243 )) < TINY(a) ) THEN
         IER = 54                                      
         RETURN
  END IF
   W( 54 ) = JVS( 243 )
   W( 82 ) = JVS( 244 )
   W( 83 ) = JVS( 245 )
   W( 100 ) = JVS( 246 )
   W( 101 ) = JVS( 247 )
   W( 104 ) = JVS( 248 )
   W( 132 ) = JVS( 249 )
   W( 136 ) = JVS( 250 )
  JVS( 243) = W( 54 )
  JVS( 244) = W( 82 )
  JVS( 245) = W( 83 )
  JVS( 246) = W( 100 )
  JVS( 247) = W( 101 )
  JVS( 248) = W( 104 )
  JVS( 249) = W( 132 )
  JVS( 250) = W( 136 )
  IF ( ABS(  JVS( 251 )) < TINY(a) ) THEN
         IER = 55                                      
         RETURN
  END IF
   W( 55 ) = JVS( 251 )
   W( 132 ) = JVS( 252 )
   W( 136 ) = JVS( 253 )
  JVS( 251) = W( 55 )
  JVS( 252) = W( 132 )
  JVS( 253) = W( 136 )
  IF ( ABS(  JVS( 254 )) < TINY(a) ) THEN
         IER = 56                                      
         RETURN
  END IF
   W( 56 ) = JVS( 254 )
   W( 126 ) = JVS( 255 )
   W( 131 ) = JVS( 256 )
   W( 133 ) = JVS( 257 )
  JVS( 254) = W( 56 )
  JVS( 255) = W( 126 )
  JVS( 256) = W( 131 )
  JVS( 257) = W( 133 )
  IF ( ABS(  JVS( 258 )) < TINY(a) ) THEN
         IER = 57                                      
         RETURN
  END IF
   W( 57 ) = JVS( 258 )
   W( 88 ) = JVS( 259 )
   W( 131 ) = JVS( 260 )
   W( 136 ) = JVS( 261 )
  JVS( 258) = W( 57 )
  JVS( 259) = W( 88 )
  JVS( 260) = W( 131 )
  JVS( 261) = W( 136 )
  IF ( ABS(  JVS( 265 )) < TINY(a) ) THEN
         IER = 58                                      
         RETURN
  END IF
   W( 22 ) = JVS( 262 )
   W( 25 ) = JVS( 263 )
   W( 27 ) = JVS( 264 )
   W( 58 ) = JVS( 265 )
   W( 131 ) = JVS( 266 )
   W( 133 ) = JVS( 267 )
   W( 136 ) = JVS( 268 )
  a = -W( 22 ) / JVS(          119  )
  W( 22 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 120 )
  a = -W( 25 ) / JVS(          127  )
  W( 25 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 128 )
  a = -W( 27 ) / JVS(          133  )
  W( 27 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 134 )
  W( 131 ) = W( 131 ) + a*JVS( 135 )
  W( 136 ) = W( 136 ) + a*JVS( 136 )
  JVS( 262) = W( 22 )
  JVS( 263) = W( 25 )
  JVS( 264) = W( 27 )
  JVS( 265) = W( 58 )
  JVS( 266) = W( 131 )
  JVS( 267) = W( 133 )
  JVS( 268) = W( 136 )
  IF ( ABS(  JVS( 272 )) < TINY(a) ) THEN
         IER = 59                                      
         RETURN
  END IF
   W( 22 ) = JVS( 269 )
   W( 25 ) = JVS( 270 )
   W( 38 ) = JVS( 271 )
   W( 59 ) = JVS( 272 )
   W( 81 ) = JVS( 273 )
   W( 132 ) = JVS( 274 )
   W( 133 ) = JVS( 275 )
   W( 134 ) = JVS( 276 )
   W( 136 ) = JVS( 277 )
  a = -W( 22 ) / JVS(          119  )
  W( 22 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 120 )
  a = -W( 25 ) / JVS(          127  )
  W( 25 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 128 )
  a = -W( 38 ) / JVS(          175  )
  W( 38 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 176 )
  JVS( 269) = W( 22 )
  JVS( 270) = W( 25 )
  JVS( 271) = W( 38 )
  JVS( 272) = W( 59 )
  JVS( 273) = W( 81 )
  JVS( 274) = W( 132 )
  JVS( 275) = W( 133 )
  JVS( 276) = W( 134 )
  JVS( 277) = W( 136 )
  IF ( ABS(  JVS( 279 )) < TINY(a) ) THEN
         IER = 60                                      
         RETURN
  END IF
   W( 55 ) = JVS( 278 )
   W( 60 ) = JVS( 279 )
   W( 132 ) = JVS( 280 )
   W( 133 ) = JVS( 281 )
   W( 136 ) = JVS( 282 )
  a = -W( 55 ) / JVS(          251  )
  W( 55 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 252 )
  W( 136 ) = W( 136 ) + a*JVS( 253 )
  JVS( 278) = W( 55 )
  JVS( 279) = W( 60 )
  JVS( 280) = W( 132 )
  JVS( 281) = W( 133 )
  JVS( 282) = W( 136 )
  IF ( ABS(  JVS( 286 )) < TINY(a) ) THEN
         IER = 61                                      
         RETURN
  END IF
   W( 24 ) = JVS( 283 )
   W( 32 ) = JVS( 284 )
   W( 37 ) = JVS( 285 )
   W( 61 ) = JVS( 286 )
   W( 131 ) = JVS( 287 )
   W( 133 ) = JVS( 288 )
   W( 136 ) = JVS( 289 )
  a = -W( 24 ) / JVS(          124  )
  W( 24 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 125 )
  a = -W( 32 ) / JVS(          153  )
  W( 32 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 154 )
  W( 131 ) = W( 131 ) + a*JVS( 155 )
  W( 136 ) = W( 136 ) + a*JVS( 156 )
  a = -W( 37 ) / JVS(          172  )
  W( 37 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 173 )
  JVS( 283) = W( 24 )
  JVS( 284) = W( 32 )
  JVS( 285) = W( 37 )
  JVS( 286) = W( 61 )
  JVS( 287) = W( 131 )
  JVS( 288) = W( 133 )
  JVS( 289) = W( 136 )
  IF ( ABS(  JVS( 290 )) < TINY(a) ) THEN
         IER = 62                                      
         RETURN
  END IF
   W( 62 ) = JVS( 290 )
   W( 84 ) = JVS( 291 )
   W( 113 ) = JVS( 292 )
   W( 133 ) = JVS( 293 )
   W( 136 ) = JVS( 294 )
  JVS( 290) = W( 62 )
  JVS( 291) = W( 84 )
  JVS( 292) = W( 113 )
  JVS( 293) = W( 133 )
  JVS( 294) = W( 136 )
  IF ( ABS(  JVS( 295 )) < TINY(a) ) THEN
         IER = 63                                      
         RETURN
  END IF
   W( 63 ) = JVS( 295 )
   W( 106 ) = JVS( 296 )
   W( 109 ) = JVS( 297 )
   W( 128 ) = JVS( 298 )
   W( 129 ) = JVS( 299 )
   W( 130 ) = JVS( 300 )
   W( 131 ) = JVS( 301 )
   W( 132 ) = JVS( 302 )
   W( 136 ) = JVS( 303 )
  JVS( 295) = W( 63 )
  JVS( 296) = W( 106 )
  JVS( 297) = W( 109 )
  JVS( 298) = W( 128 )
  JVS( 299) = W( 129 )
  JVS( 300) = W( 130 )
  JVS( 301) = W( 131 )
  JVS( 302) = W( 132 )
  JVS( 303) = W( 136 )
  IF ( ABS(  JVS( 304 )) < TINY(a) ) THEN
         IER = 64                                      
         RETURN
  END IF
   W( 64 ) = JVS( 304 )
   W( 128 ) = JVS( 305 )
   W( 130 ) = JVS( 306 )
   W( 131 ) = JVS( 307 )
   W( 136 ) = JVS( 308 )
  JVS( 304) = W( 64 )
  JVS( 305) = W( 128 )
  JVS( 306) = W( 130 )
  JVS( 307) = W( 131 )
  JVS( 308) = W( 136 )
  IF ( ABS(  JVS( 309 )) < TINY(a) ) THEN
         IER = 65                                      
         RETURN
  END IF
   W( 65 ) = JVS( 309 )
   W( 118 ) = JVS( 310 )
   W( 131 ) = JVS( 311 )
   W( 136 ) = JVS( 312 )
  JVS( 309) = W( 65 )
  JVS( 310) = W( 118 )
  JVS( 311) = W( 131 )
  JVS( 312) = W( 136 )
  IF ( ABS(  JVS( 315 )) < TINY(a) ) THEN
         IER = 66                                      
         RETURN
  END IF
   W( 20 ) = JVS( 313 )
   W( 29 ) = JVS( 314 )
   W( 66 ) = JVS( 315 )
   W( 131 ) = JVS( 316 )
   W( 133 ) = JVS( 317 )
   W( 136 ) = JVS( 318 )
  a = -W( 20 ) / JVS(          113  )
  W( 20 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 114 )
  a = -W( 29 ) / JVS(          141  )
  W( 29 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 142 )
  W( 131 ) = W( 131 ) + a*JVS( 143 )
  W( 136 ) = W( 136 ) + a*JVS( 144 )
  JVS( 313) = W( 20 )
  JVS( 314) = W( 29 )
  JVS( 315) = W( 66 )
  JVS( 316) = W( 131 )
  JVS( 317) = W( 133 )
  JVS( 318) = W( 136 )
  IF ( ABS(  JVS( 319 )) < TINY(a) ) THEN
         IER = 67                                      
         RETURN
  END IF
   W( 67 ) = JVS( 319 )
   W( 102 ) = JVS( 320 )
   W( 109 ) = JVS( 321 )
   W( 113 ) = JVS( 322 )
   W( 133 ) = JVS( 323 )
   W( 135 ) = JVS( 324 )
   W( 136 ) = JVS( 325 )
  JVS( 319) = W( 67 )
  JVS( 320) = W( 102 )
  JVS( 321) = W( 109 )
  JVS( 322) = W( 113 )
  JVS( 323) = W( 133 )
  JVS( 324) = W( 135 )
  JVS( 325) = W( 136 )
  IF ( ABS(  JVS( 326 )) < TINY(a) ) THEN
         IER = 68                                      
         RETURN
  END IF
   W( 68 ) = JVS( 326 )
   W( 128 ) = JVS( 327 )
   W( 134 ) = JVS( 328 )
   W( 136 ) = JVS( 329 )
  JVS( 326) = W( 68 )
  JVS( 327) = W( 128 )
  JVS( 328) = W( 134 )
  JVS( 329) = W( 136 )
  IF ( ABS(  JVS( 332 )) < TINY(a) ) THEN
         IER = 69                                      
         RETURN
  END IF
   W( 34 ) = JVS( 330 )
   W( 38 ) = JVS( 331 )
   W( 69 ) = JVS( 332 )
   W( 131 ) = JVS( 333 )
   W( 133 ) = JVS( 334 )
   W( 136 ) = JVS( 335 )
  a = -W( 34 ) / JVS(          161  )
  W( 34 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 162 )
  W( 131 ) = W( 131 ) + a*JVS( 163 )
  W( 136 ) = W( 136 ) + a*JVS( 164 )
  a = -W( 38 ) / JVS(          175  )
  W( 38 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 176 )
  JVS( 330) = W( 34 )
  JVS( 331) = W( 38 )
  JVS( 332) = W( 69 )
  JVS( 333) = W( 131 )
  JVS( 334) = W( 133 )
  JVS( 335) = W( 136 )
  IF ( ABS(  JVS( 336 )) < TINY(a) ) THEN
         IER = 70                                      
         RETURN
  END IF
   W( 70 ) = JVS( 336 )
   W( 86 ) = JVS( 337 )
   W( 131 ) = JVS( 338 )
   W( 136 ) = JVS( 339 )
  JVS( 336) = W( 70 )
  JVS( 337) = W( 86 )
  JVS( 338) = W( 131 )
  JVS( 339) = W( 136 )
  IF ( ABS(  JVS( 342 )) < TINY(a) ) THEN
         IER = 71                                      
         RETURN
  END IF
   W( 36 ) = JVS( 340 )
   W( 39 ) = JVS( 341 )
   W( 71 ) = JVS( 342 )
   W( 112 ) = JVS( 343 )
   W( 117 ) = JVS( 344 )
   W( 126 ) = JVS( 345 )
   W( 134 ) = JVS( 346 )
   W( 135 ) = JVS( 347 )
   W( 136 ) = JVS( 348 )
  a = -W( 36 ) / JVS(          169  )
  W( 36 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 170 )
  W( 135 ) = W( 135 ) + a*JVS( 171 )
  a = -W( 39 ) / JVS(          177  )
  W( 39 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 178 )
  W( 136 ) = W( 136 ) + a*JVS( 179 )
  JVS( 340) = W( 36 )
  JVS( 341) = W( 39 )
  JVS( 342) = W( 71 )
  JVS( 343) = W( 112 )
  JVS( 344) = W( 117 )
  JVS( 345) = W( 126 )
  JVS( 346) = W( 134 )
  JVS( 347) = W( 135 )
  JVS( 348) = W( 136 )
  IF ( ABS(  JVS( 349 )) < TINY(a) ) THEN
         IER = 72                                      
         RETURN
  END IF
   W( 72 ) = JVS( 349 )
   W( 80 ) = JVS( 350 )
   W( 109 ) = JVS( 351 )
   W( 132 ) = JVS( 352 )
   W( 136 ) = JVS( 353 )
  JVS( 349) = W( 72 )
  JVS( 350) = W( 80 )
  JVS( 351) = W( 109 )
  JVS( 352) = W( 132 )
  JVS( 353) = W( 136 )
  IF ( ABS(  JVS( 354 )) < TINY(a) ) THEN
         IER = 73                                      
         RETURN
  END IF
   W( 73 ) = JVS( 354 )
   W( 103 ) = JVS( 355 )
   W( 131 ) = JVS( 356 )
   W( 133 ) = JVS( 357 )
   W( 136 ) = JVS( 358 )
  JVS( 354) = W( 73 )
  JVS( 355) = W( 103 )
  JVS( 356) = W( 131 )
  JVS( 357) = W( 133 )
  JVS( 358) = W( 136 )
  IF ( ABS(  JVS( 362 )) < TINY(a) ) THEN
         IER = 74                                      
         RETURN
  END IF
   W( 17 ) = JVS( 359 )
   W( 33 ) = JVS( 360 )
   W( 61 ) = JVS( 361 )
   W( 74 ) = JVS( 362 )
   W( 131 ) = JVS( 363 )
   W( 133 ) = JVS( 364 )
   W( 134 ) = JVS( 365 )
   W( 136 ) = JVS( 366 )
  a = -W( 17 ) / JVS(          105  )
  W( 17 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 106 )
  W( 134 ) = W( 134 ) + a*JVS( 107 )
  a = -W( 33 ) / JVS(          157  )
  W( 33 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 158 )
  W( 133 ) = W( 133 ) + a*JVS( 159 )
  W( 136 ) = W( 136 ) + a*JVS( 160 )
  a = -W( 61 ) / JVS(          286  )
  W( 61 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 287 )
  W( 133 ) = W( 133 ) + a*JVS( 288 )
  W( 136 ) = W( 136 ) + a*JVS( 289 )
  JVS( 359) = W( 17 )
  JVS( 360) = W( 33 )
  JVS( 361) = W( 61 )
  JVS( 362) = W( 74 )
  JVS( 363) = W( 131 )
  JVS( 364) = W( 133 )
  JVS( 365) = W( 134 )
  JVS( 366) = W( 136 )
  IF ( ABS(  JVS( 367 )) < TINY(a) ) THEN
         IER = 75                                      
         RETURN
  END IF
   W( 75 ) = JVS( 367 )
   W( 127 ) = JVS( 368 )
   W( 131 ) = JVS( 369 )
   W( 136 ) = JVS( 370 )
  JVS( 367) = W( 75 )
  JVS( 368) = W( 127 )
  JVS( 369) = W( 131 )
  JVS( 370) = W( 136 )
  IF ( ABS(  JVS( 371 )) < TINY(a) ) THEN
         IER = 76                                      
         RETURN
  END IF
   W( 76 ) = JVS( 371 )
   W( 111 ) = JVS( 372 )
   W( 131 ) = JVS( 373 )
   W( 136 ) = JVS( 374 )
  JVS( 371) = W( 76 )
  JVS( 372) = W( 111 )
  JVS( 373) = W( 131 )
  JVS( 374) = W( 136 )
  IF ( ABS(  JVS( 375 )) < TINY(a) ) THEN
         IER = 77                                      
         RETURN
  END IF
   W( 77 ) = JVS( 375 )
   W( 130 ) = JVS( 376 )
   W( 134 ) = JVS( 377 )
   W( 136 ) = JVS( 378 )
  JVS( 375) = W( 77 )
  JVS( 376) = W( 130 )
  JVS( 377) = W( 134 )
  JVS( 378) = W( 136 )
  IF ( ABS(  JVS( 382 )) < TINY(a) ) THEN
         IER = 78                                      
         RETURN
  END IF
   W( 19 ) = JVS( 379 )
   W( 55 ) = JVS( 380 )
   W( 56 ) = JVS( 381 )
   W( 78 ) = JVS( 382 )
   W( 82 ) = JVS( 383 )
   W( 83 ) = JVS( 384 )
   W( 89 ) = JVS( 385 )
   W( 100 ) = JVS( 386 )
   W( 101 ) = JVS( 387 )
   W( 104 ) = JVS( 388 )
   W( 126 ) = JVS( 389 )
   W( 131 ) = JVS( 390 )
   W( 132 ) = JVS( 391 )
   W( 133 ) = JVS( 392 )
   W( 136 ) = JVS( 393 )
  a = -W( 19 ) / JVS(          111  )
  W( 19 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 112 )
  a = -W( 55 ) / JVS(          251  )
  W( 55 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 252 )
  W( 136 ) = W( 136 ) + a*JVS( 253 )
  a = -W( 56 ) / JVS(          254  )
  W( 56 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 255 )
  W( 131 ) = W( 131 ) + a*JVS( 256 )
  W( 133 ) = W( 133 ) + a*JVS( 257 )
  JVS( 379) = W( 19 )
  JVS( 380) = W( 55 )
  JVS( 381) = W( 56 )
  JVS( 382) = W( 78 )
  JVS( 383) = W( 82 )
  JVS( 384) = W( 83 )
  JVS( 385) = W( 89 )
  JVS( 386) = W( 100 )
  JVS( 387) = W( 101 )
  JVS( 388) = W( 104 )
  JVS( 389) = W( 126 )
  JVS( 390) = W( 131 )
  JVS( 391) = W( 132 )
  JVS( 392) = W( 133 )
  JVS( 393) = W( 136 )
  IF ( ABS(  JVS( 394 )) < TINY(a) ) THEN
         IER = 79                                      
         RETURN
  END IF
   W( 79 ) = JVS( 394 )
   W( 120 ) = JVS( 395 )
   W( 131 ) = JVS( 396 )
   W( 136 ) = JVS( 397 )
  JVS( 394) = W( 79 )
  JVS( 395) = W( 120 )
  JVS( 396) = W( 131 )
  JVS( 397) = W( 136 )
  IF ( ABS(  JVS( 401 )) < TINY(a) ) THEN
         IER = 80                                      
         RETURN
  END IF
   W( 15 ) = JVS( 398 )
   W( 44 ) = JVS( 399 )
   W( 72 ) = JVS( 400 )
   W( 80 ) = JVS( 401 )
   W( 109 ) = JVS( 402 )
   W( 126 ) = JVS( 403 )
   W( 132 ) = JVS( 404 )
   W( 136 ) = JVS( 405 )
  a = -W( 15 ) / JVS(           99  )
  W( 15 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 100 )
  a = -W( 44 ) / JVS(          198  )
  W( 44 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 199 )
  W( 80 ) = W( 80 ) + a*JVS( 200 )
  W( 126 ) = W( 126 ) + a*JVS( 201 )
  W( 136 ) = W( 136 ) + a*JVS( 202 )
  a = -W( 72 ) / JVS(          349  )
  W( 72 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 350 )
  W( 109 ) = W( 109 ) + a*JVS( 351 )
  W( 132 ) = W( 132 ) + a*JVS( 352 )
  W( 136 ) = W( 136 ) + a*JVS( 353 )
  JVS( 398) = W( 15 )
  JVS( 399) = W( 44 )
  JVS( 400) = W( 72 )
  JVS( 401) = W( 80 )
  JVS( 402) = W( 109 )
  JVS( 403) = W( 126 )
  JVS( 404) = W( 132 )
  JVS( 405) = W( 136 )
  IF ( ABS(  JVS( 409 )) < TINY(a) ) THEN
         IER = 81                                      
         RETURN
  END IF
   W( 28 ) = JVS( 406 )
   W( 59 ) = JVS( 407 )
   W( 74 ) = JVS( 408 )
   W( 81 ) = JVS( 409 )
   W( 131 ) = JVS( 410 )
   W( 132 ) = JVS( 411 )
   W( 133 ) = JVS( 412 )
   W( 134 ) = JVS( 413 )
   W( 136 ) = JVS( 414 )
  a = -W( 28 ) / JVS(          137  )
  W( 28 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 138 )
  W( 131 ) = W( 131 ) + a*JVS( 139 )
  W( 136 ) = W( 136 ) + a*JVS( 140 )
  a = -W( 59 ) / JVS(          272  )
  W( 59 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 273 )
  W( 132 ) = W( 132 ) + a*JVS( 274 )
  W( 133 ) = W( 133 ) + a*JVS( 275 )
  W( 134 ) = W( 134 ) + a*JVS( 276 )
  W( 136 ) = W( 136 ) + a*JVS( 277 )
  a = -W( 74 ) / JVS(          362  )
  W( 74 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 363 )
  W( 133 ) = W( 133 ) + a*JVS( 364 )
  W( 134 ) = W( 134 ) + a*JVS( 365 )
  W( 136 ) = W( 136 ) + a*JVS( 366 )
  JVS( 406) = W( 28 )
  JVS( 407) = W( 59 )
  JVS( 408) = W( 74 )
  JVS( 409) = W( 81 )
  JVS( 410) = W( 131 )
  JVS( 411) = W( 132 )
  JVS( 412) = W( 133 )
  JVS( 413) = W( 134 )
  JVS( 414) = W( 136 )
  IF ( ABS(  JVS( 415 )) < TINY(a) ) THEN
         IER = 82                                      
         RETURN
  END IF
   W( 82 ) = JVS( 415 )
   W( 132 ) = JVS( 416 )
   W( 136 ) = JVS( 417 )
  JVS( 415) = W( 82 )
  JVS( 416) = W( 132 )
  JVS( 417) = W( 136 )
  IF ( ABS(  JVS( 418 )) < TINY(a) ) THEN
         IER = 83                                      
         RETURN
  END IF
   W( 83 ) = JVS( 418 )
   W( 132 ) = JVS( 419 )
   W( 136 ) = JVS( 420 )
  JVS( 418) = W( 83 )
  JVS( 419) = W( 132 )
  JVS( 420) = W( 136 )
  IF ( ABS(  JVS( 421 )) < TINY(a) ) THEN
         IER = 84                                      
         RETURN
  END IF
   W( 84 ) = JVS( 421 )
   W( 113 ) = JVS( 422 )
   W( 131 ) = JVS( 423 )
   W( 136 ) = JVS( 424 )
  JVS( 421) = W( 84 )
  JVS( 422) = W( 113 )
  JVS( 423) = W( 131 )
  JVS( 424) = W( 136 )
  IF ( ABS(  JVS( 428 )) < TINY(a) ) THEN
         IER = 85                                      
         RETURN
  END IF
   W( 49 ) = JVS( 425 )
   W( 62 ) = JVS( 426 )
   W( 84 ) = JVS( 427 )
   W( 85 ) = JVS( 428 )
   W( 113 ) = JVS( 429 )
   W( 131 ) = JVS( 430 )
   W( 133 ) = JVS( 431 )
   W( 136 ) = JVS( 432 )
  a = -W( 49 ) / JVS(          219  )
  W( 49 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 220 )
  W( 131 ) = W( 131 ) + a*JVS( 221 )
  W( 136 ) = W( 136 ) + a*JVS( 222 )
  a = -W( 62 ) / JVS(          290  )
  W( 62 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 291 )
  W( 113 ) = W( 113 ) + a*JVS( 292 )
  W( 133 ) = W( 133 ) + a*JVS( 293 )
  W( 136 ) = W( 136 ) + a*JVS( 294 )
  a = -W( 84 ) / JVS(          421  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 422 )
  W( 131 ) = W( 131 ) + a*JVS( 423 )
  W( 136 ) = W( 136 ) + a*JVS( 424 )
  JVS( 425) = W( 49 )
  JVS( 426) = W( 62 )
  JVS( 427) = W( 84 )
  JVS( 428) = W( 85 )
  JVS( 429) = W( 113 )
  JVS( 430) = W( 131 )
  JVS( 431) = W( 133 )
  JVS( 432) = W( 136 )
  IF ( ABS(  JVS( 434 )) < TINY(a) ) THEN
         IER = 86                                      
         RETURN
  END IF
   W( 70 ) = JVS( 433 )
   W( 86 ) = JVS( 434 )
   W( 109 ) = JVS( 435 )
   W( 131 ) = JVS( 436 )
   W( 133 ) = JVS( 437 )
   W( 136 ) = JVS( 438 )
  a = -W( 70 ) / JVS(          336  )
  W( 70 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 337 )
  W( 131 ) = W( 131 ) + a*JVS( 338 )
  W( 136 ) = W( 136 ) + a*JVS( 339 )
  JVS( 433) = W( 70 )
  JVS( 434) = W( 86 )
  JVS( 435) = W( 109 )
  JVS( 436) = W( 131 )
  JVS( 437) = W( 133 )
  JVS( 438) = W( 136 )
  IF ( ABS(  JVS( 441 )) < TINY(a) ) THEN
         IER = 87                                      
         RETURN
  END IF
   W( 43 ) = JVS( 439 )
   W( 66 ) = JVS( 440 )
   W( 87 ) = JVS( 441 )
   W( 88 ) = JVS( 442 )
   W( 90 ) = JVS( 443 )
   W( 131 ) = JVS( 444 )
   W( 133 ) = JVS( 445 )
   W( 134 ) = JVS( 446 )
   W( 136 ) = JVS( 447 )
  a = -W( 43 ) / JVS(          192  )
  W( 43 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 193 )
  W( 88 ) = W( 88 ) + a*JVS( 194 )
  W( 90 ) = W( 90 ) + a*JVS( 195 )
  W( 133 ) = W( 133 ) + a*JVS( 196 )
  W( 136 ) = W( 136 ) + a*JVS( 197 )
  a = -W( 66 ) / JVS(          315  )
  W( 66 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 316 )
  W( 133 ) = W( 133 ) + a*JVS( 317 )
  W( 136 ) = W( 136 ) + a*JVS( 318 )
  JVS( 439) = W( 43 )
  JVS( 440) = W( 66 )
  JVS( 441) = W( 87 )
  JVS( 442) = W( 88 )
  JVS( 443) = W( 90 )
  JVS( 444) = W( 131 )
  JVS( 445) = W( 133 )
  JVS( 446) = W( 134 )
  JVS( 447) = W( 136 )
  IF ( ABS(  JVS( 450 )) < TINY(a) ) THEN
         IER = 88                                      
         RETURN
  END IF
   W( 24 ) = JVS( 448 )
   W( 57 ) = JVS( 449 )
   W( 88 ) = JVS( 450 )
   W( 131 ) = JVS( 451 )
   W( 133 ) = JVS( 452 )
   W( 136 ) = JVS( 453 )
  a = -W( 24 ) / JVS(          124  )
  W( 24 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 125 )
  a = -W( 57 ) / JVS(          258  )
  W( 57 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 259 )
  W( 131 ) = W( 131 ) + a*JVS( 260 )
  W( 136 ) = W( 136 ) + a*JVS( 261 )
  JVS( 448) = W( 24 )
  JVS( 449) = W( 57 )
  JVS( 450) = W( 88 )
  JVS( 451) = W( 131 )
  JVS( 452) = W( 133 )
  JVS( 453) = W( 136 )
  IF ( ABS(  JVS( 454 )) < TINY(a) ) THEN
         IER = 89                                      
         RETURN
  END IF
   W( 89 ) = JVS( 454 )
   W( 132 ) = JVS( 455 )
   W( 135 ) = JVS( 456 )
   W( 136 ) = JVS( 457 )
  JVS( 454) = W( 89 )
  JVS( 455) = W( 132 )
  JVS( 456) = W( 135 )
  JVS( 457) = W( 136 )
  IF ( ABS(  JVS( 460 )) < TINY(a) ) THEN
         IER = 90                                      
         RETURN
  END IF
   W( 35 ) = JVS( 458 )
   W( 37 ) = JVS( 459 )
   W( 90 ) = JVS( 460 )
   W( 131 ) = JVS( 461 )
   W( 133 ) = JVS( 462 )
   W( 136 ) = JVS( 463 )
  a = -W( 35 ) / JVS(          165  )
  W( 35 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 166 )
  W( 131 ) = W( 131 ) + a*JVS( 167 )
  W( 136 ) = W( 136 ) + a*JVS( 168 )
  a = -W( 37 ) / JVS(          172  )
  W( 37 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 173 )
  JVS( 458) = W( 35 )
  JVS( 459) = W( 37 )
  JVS( 460) = W( 90 )
  JVS( 461) = W( 131 )
  JVS( 462) = W( 133 )
  JVS( 463) = W( 136 )
  IF ( ABS(  JVS( 467 )) < TINY(a) ) THEN
         IER = 91                                      
         RETURN
  END IF
   W( 30 ) = JVS( 464 )
   W( 88 ) = JVS( 465 )
   W( 90 ) = JVS( 466 )
   W( 91 ) = JVS( 467 )
   W( 131 ) = JVS( 468 )
   W( 133 ) = JVS( 469 )
   W( 134 ) = JVS( 470 )
   W( 136 ) = JVS( 471 )
  a = -W( 30 ) / JVS(          145  )
  W( 30 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 146 )
  W( 90 ) = W( 90 ) + a*JVS( 147 )
  W( 133 ) = W( 133 ) + a*JVS( 148 )
  a = -W( 88 ) / JVS(          450  )
  W( 88 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 451 )
  W( 133 ) = W( 133 ) + a*JVS( 452 )
  W( 136 ) = W( 136 ) + a*JVS( 453 )
  a = -W( 90 ) / JVS(          460  )
  W( 90 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 461 )
  W( 133 ) = W( 133 ) + a*JVS( 462 )
  W( 136 ) = W( 136 ) + a*JVS( 463 )
  JVS( 464) = W( 30 )
  JVS( 465) = W( 88 )
  JVS( 466) = W( 90 )
  JVS( 467) = W( 91 )
  JVS( 468) = W( 131 )
  JVS( 469) = W( 133 )
  JVS( 470) = W( 134 )
  JVS( 471) = W( 136 )
  IF ( ABS(  JVS( 473 )) < TINY(a) ) THEN
         IER = 92                                      
         RETURN
  END IF
   W( 80 ) = JVS( 472 )
   W( 92 ) = JVS( 473 )
   W( 109 ) = JVS( 474 )
   W( 126 ) = JVS( 475 )
   W( 131 ) = JVS( 476 )
   W( 132 ) = JVS( 477 )
   W( 134 ) = JVS( 478 )
   W( 136 ) = JVS( 479 )
  a = -W( 80 ) / JVS(          401  )
  W( 80 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 402 )
  W( 126 ) = W( 126 ) + a*JVS( 403 )
  W( 132 ) = W( 132 ) + a*JVS( 404 )
  W( 136 ) = W( 136 ) + a*JVS( 405 )
  JVS( 472) = W( 80 )
  JVS( 473) = W( 92 )
  JVS( 474) = W( 109 )
  JVS( 475) = W( 126 )
  JVS( 476) = W( 131 )
  JVS( 477) = W( 132 )
  JVS( 478) = W( 134 )
  JVS( 479) = W( 136 )
  IF ( ABS(  JVS( 482 )) < TINY(a) ) THEN
         IER = 93                                      
         RETURN
  END IF
   W( 13 ) = JVS( 480 )
   W( 46 ) = JVS( 481 )
   W( 93 ) = JVS( 482 )
   W( 129 ) = JVS( 483 )
   W( 131 ) = JVS( 484 )
   W( 133 ) = JVS( 485 )
   W( 136 ) = JVS( 486 )
  a = -W( 13 ) / JVS(           95  )
  W( 13 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 96 )
  a = -W( 46 ) / JVS(          207  )
  W( 46 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 208 )
  W( 131 ) = W( 131 ) + a*JVS( 209 )
  W( 136 ) = W( 136 ) + a*JVS( 210 )
  JVS( 480) = W( 13 )
  JVS( 481) = W( 46 )
  JVS( 482) = W( 93 )
  JVS( 483) = W( 129 )
  JVS( 484) = W( 131 )
  JVS( 485) = W( 133 )
  JVS( 486) = W( 136 )
  IF ( ABS(  JVS( 490 )) < TINY(a) ) THEN
         IER = 94                                      
         RETURN
  END IF
   W( 42 ) = JVS( 487 )
   W( 88 ) = JVS( 488 )
   W( 90 ) = JVS( 489 )
   W( 94 ) = JVS( 490 )
   W( 131 ) = JVS( 491 )
   W( 133 ) = JVS( 492 )
   W( 134 ) = JVS( 493 )
   W( 136 ) = JVS( 494 )
  a = -W( 42 ) / JVS(          187  )
  W( 42 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 188 )
  W( 90 ) = W( 90 ) + a*JVS( 189 )
  W( 133 ) = W( 133 ) + a*JVS( 190 )
  a = -W( 88 ) / JVS(          450  )
  W( 88 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 451 )
  W( 133 ) = W( 133 ) + a*JVS( 452 )
  W( 136 ) = W( 136 ) + a*JVS( 453 )
  a = -W( 90 ) / JVS(          460  )
  W( 90 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 461 )
  W( 133 ) = W( 133 ) + a*JVS( 462 )
  W( 136 ) = W( 136 ) + a*JVS( 463 )
  JVS( 487) = W( 42 )
  JVS( 488) = W( 88 )
  JVS( 489) = W( 90 )
  JVS( 490) = W( 94 )
  JVS( 491) = W( 131 )
  JVS( 492) = W( 133 )
  JVS( 493) = W( 134 )
  JVS( 494) = W( 136 )
  IF ( ABS(  JVS( 495 )) < TINY(a) ) THEN
         IER = 95                                      
         RETURN
  END IF
   W( 95 ) = JVS( 495 )
   W( 97 ) = JVS( 496 )
   W( 111 ) = JVS( 497 )
   W( 118 ) = JVS( 498 )
   W( 119 ) = JVS( 499 )
   W( 120 ) = JVS( 500 )
   W( 121 ) = JVS( 501 )
   W( 123 ) = JVS( 502 )
   W( 125 ) = JVS( 503 )
   W( 127 ) = JVS( 504 )
   W( 129 ) = JVS( 505 )
   W( 136 ) = JVS( 506 )
  JVS( 495) = W( 95 )
  JVS( 496) = W( 97 )
  JVS( 497) = W( 111 )
  JVS( 498) = W( 118 )
  JVS( 499) = W( 119 )
  JVS( 500) = W( 120 )
  JVS( 501) = W( 121 )
  JVS( 502) = W( 123 )
  JVS( 503) = W( 125 )
  JVS( 504) = W( 127 )
  JVS( 505) = W( 129 )
  JVS( 506) = W( 136 )
  IF ( ABS(  JVS( 516 )) < TINY(a) ) THEN
         IER = 96                                      
         RETURN
  END IF
   W( 19 ) = JVS( 507 )
   W( 57 ) = JVS( 508 )
   W( 58 ) = JVS( 509 )
   W( 66 ) = JVS( 510 )
   W( 69 ) = JVS( 511 )
   W( 87 ) = JVS( 512 )
   W( 88 ) = JVS( 513 )
   W( 90 ) = JVS( 514 )
   W( 94 ) = JVS( 515 )
   W( 96 ) = JVS( 516 )
   W( 98 ) = JVS( 517 )
   W( 108 ) = JVS( 518 )
   W( 123 ) = JVS( 519 )
   W( 127 ) = JVS( 520 )
   W( 128 ) = JVS( 521 )
   W( 129 ) = JVS( 522 )
   W( 131 ) = JVS( 523 )
   W( 133 ) = JVS( 524 )
   W( 134 ) = JVS( 525 )
   W( 135 ) = JVS( 526 )
   W( 136 ) = JVS( 527 )
  a = -W( 19 ) / JVS(          111  )
  W( 19 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 112 )
  a = -W( 57 ) / JVS(          258  )
  W( 57 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 259 )
  W( 131 ) = W( 131 ) + a*JVS( 260 )
  W( 136 ) = W( 136 ) + a*JVS( 261 )
  a = -W( 58 ) / JVS(          265  )
  W( 58 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 266 )
  W( 133 ) = W( 133 ) + a*JVS( 267 )
  W( 136 ) = W( 136 ) + a*JVS( 268 )
  a = -W( 66 ) / JVS(          315  )
  W( 66 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 316 )
  W( 133 ) = W( 133 ) + a*JVS( 317 )
  W( 136 ) = W( 136 ) + a*JVS( 318 )
  a = -W( 69 ) / JVS(          332  )
  W( 69 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 333 )
  W( 133 ) = W( 133 ) + a*JVS( 334 )
  W( 136 ) = W( 136 ) + a*JVS( 335 )
  a = -W( 87 ) / JVS(          441  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 90 ) = W( 90 ) + a*JVS( 443 )
  W( 131 ) = W( 131 ) + a*JVS( 444 )
  W( 133 ) = W( 133 ) + a*JVS( 445 )
  W( 134 ) = W( 134 ) + a*JVS( 446 )
  W( 136 ) = W( 136 ) + a*JVS( 447 )
  a = -W( 88 ) / JVS(          450  )
  W( 88 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 451 )
  W( 133 ) = W( 133 ) + a*JVS( 452 )
  W( 136 ) = W( 136 ) + a*JVS( 453 )
  a = -W( 90 ) / JVS(          460  )
  W( 90 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 461 )
  W( 133 ) = W( 133 ) + a*JVS( 462 )
  W( 136 ) = W( 136 ) + a*JVS( 463 )
  a = -W( 94 ) / JVS(          490  )
  W( 94 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 491 )
  W( 133 ) = W( 133 ) + a*JVS( 492 )
  W( 134 ) = W( 134 ) + a*JVS( 493 )
  W( 136 ) = W( 136 ) + a*JVS( 494 )
  JVS( 507) = W( 19 )
  JVS( 508) = W( 57 )
  JVS( 509) = W( 58 )
  JVS( 510) = W( 66 )
  JVS( 511) = W( 69 )
  JVS( 512) = W( 87 )
  JVS( 513) = W( 88 )
  JVS( 514) = W( 90 )
  JVS( 515) = W( 94 )
  JVS( 516) = W( 96 )
  JVS( 517) = W( 98 )
  JVS( 518) = W( 108 )
  JVS( 519) = W( 123 )
  JVS( 520) = W( 127 )
  JVS( 521) = W( 128 )
  JVS( 522) = W( 129 )
  JVS( 523) = W( 131 )
  JVS( 524) = W( 133 )
  JVS( 525) = W( 134 )
  JVS( 526) = W( 135 )
  JVS( 527) = W( 136 )
  IF ( ABS(  JVS( 532 )) < TINY(a) ) THEN
         IER = 97                                      
         RETURN
  END IF
   W( 12 ) = JVS( 528 )
   W( 45 ) = JVS( 529 )
   W( 62 ) = JVS( 530 )
   W( 84 ) = JVS( 531 )
   W( 97 ) = JVS( 532 )
   W( 113 ) = JVS( 533 )
   W( 129 ) = JVS( 534 )
   W( 131 ) = JVS( 535 )
   W( 133 ) = JVS( 536 )
   W( 136 ) = JVS( 537 )
  a = -W( 12 ) / JVS(           93  )
  W( 12 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 94 )
  a = -W( 45 ) / JVS(          203  )
  W( 45 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 204 )
  W( 131 ) = W( 131 ) + a*JVS( 205 )
  W( 136 ) = W( 136 ) + a*JVS( 206 )
  a = -W( 62 ) / JVS(          290  )
  W( 62 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 291 )
  W( 113 ) = W( 113 ) + a*JVS( 292 )
  W( 133 ) = W( 133 ) + a*JVS( 293 )
  W( 136 ) = W( 136 ) + a*JVS( 294 )
  a = -W( 84 ) / JVS(          421  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 422 )
  W( 131 ) = W( 131 ) + a*JVS( 423 )
  W( 136 ) = W( 136 ) + a*JVS( 424 )
  JVS( 528) = W( 12 )
  JVS( 529) = W( 45 )
  JVS( 530) = W( 62 )
  JVS( 531) = W( 84 )
  JVS( 532) = W( 97 )
  JVS( 533) = W( 113 )
  JVS( 534) = W( 129 )
  JVS( 535) = W( 131 )
  JVS( 536) = W( 133 )
  JVS( 537) = W( 136 )
  IF ( ABS(  JVS( 542 )) < TINY(a) ) THEN
         IER = 98                                      
         RETURN
  END IF
   W( 57 ) = JVS( 538 )
   W( 82 ) = JVS( 539 )
   W( 83 ) = JVS( 540 )
   W( 88 ) = JVS( 541 )
   W( 98 ) = JVS( 542 )
   W( 100 ) = JVS( 543 )
   W( 101 ) = JVS( 544 )
   W( 104 ) = JVS( 545 )
   W( 131 ) = JVS( 546 )
   W( 132 ) = JVS( 547 )
   W( 133 ) = JVS( 548 )
   W( 136 ) = JVS( 549 )
  a = -W( 57 ) / JVS(          258  )
  W( 57 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 259 )
  W( 131 ) = W( 131 ) + a*JVS( 260 )
  W( 136 ) = W( 136 ) + a*JVS( 261 )
  a = -W( 82 ) / JVS(          415  )
  W( 82 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 416 )
  W( 136 ) = W( 136 ) + a*JVS( 417 )
  a = -W( 83 ) / JVS(          418  )
  W( 83 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 88 ) / JVS(          450  )
  W( 88 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 451 )
  W( 133 ) = W( 133 ) + a*JVS( 452 )
  W( 136 ) = W( 136 ) + a*JVS( 453 )
  JVS( 538) = W( 57 )
  JVS( 539) = W( 82 )
  JVS( 540) = W( 83 )
  JVS( 541) = W( 88 )
  JVS( 542) = W( 98 )
  JVS( 543) = W( 100 )
  JVS( 544) = W( 101 )
  JVS( 545) = W( 104 )
  JVS( 546) = W( 131 )
  JVS( 547) = W( 132 )
  JVS( 548) = W( 133 )
  JVS( 549) = W( 136 )
  IF ( ABS(  JVS( 550 )) < TINY(a) ) THEN
         IER = 99                                      
         RETURN
  END IF
   W( 99 ) = JVS( 550 )
   W( 106 ) = JVS( 551 )
   W( 131 ) = JVS( 552 )
   W( 133 ) = JVS( 553 )
   W( 135 ) = JVS( 554 )
  JVS( 550) = W( 99 )
  JVS( 551) = W( 106 )
  JVS( 552) = W( 131 )
  JVS( 553) = W( 133 )
  JVS( 554) = W( 135 )
  IF ( ABS(  JVS( 555 )) < TINY(a) ) THEN
         IER = 100                                     
         RETURN
  END IF
   W( 100 ) = JVS( 555 )
   W( 132 ) = JVS( 556 )
   W( 135 ) = JVS( 557 )
   W( 136 ) = JVS( 558 )
  JVS( 555) = W( 100 )
  JVS( 556) = W( 132 )
  JVS( 557) = W( 135 )
  JVS( 558) = W( 136 )
  IF ( ABS(  JVS( 559 )) < TINY(a) ) THEN
         IER = 101                                     
         RETURN
  END IF
   W( 101 ) = JVS( 559 )
   W( 132 ) = JVS( 560 )
   W( 135 ) = JVS( 561 )
   W( 136 ) = JVS( 562 )
  JVS( 559) = W( 101 )
  JVS( 560) = W( 132 )
  JVS( 561) = W( 135 )
  JVS( 562) = W( 136 )
  IF ( ABS(  JVS( 564 )) < TINY(a) ) THEN
         IER = 102                                     
         RETURN
  END IF
   W( 89 ) = JVS( 563 )
   W( 102 ) = JVS( 564 )
   W( 131 ) = JVS( 565 )
   W( 132 ) = JVS( 566 )
   W( 133 ) = JVS( 567 )
   W( 135 ) = JVS( 568 )
   W( 136 ) = JVS( 569 )
  a = -W( 89 ) / JVS(          454  )
  W( 89 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 455 )
  W( 135 ) = W( 135 ) + a*JVS( 456 )
  W( 136 ) = W( 136 ) + a*JVS( 457 )
  JVS( 563) = W( 89 )
  JVS( 564) = W( 102 )
  JVS( 565) = W( 131 )
  JVS( 566) = W( 132 )
  JVS( 567) = W( 133 )
  JVS( 568) = W( 135 )
  JVS( 569) = W( 136 )
  IF ( ABS(  JVS( 573 )) < TINY(a) ) THEN
         IER = 103                                     
         RETURN
  END IF
   W( 76 ) = JVS( 570 )
   W( 89 ) = JVS( 571 )
   W( 102 ) = JVS( 572 )
   W( 103 ) = JVS( 573 )
   W( 111 ) = JVS( 574 )
   W( 129 ) = JVS( 575 )
   W( 131 ) = JVS( 576 )
   W( 132 ) = JVS( 577 )
   W( 133 ) = JVS( 578 )
   W( 135 ) = JVS( 579 )
   W( 136 ) = JVS( 580 )
  a = -W( 76 ) / JVS(          371  )
  W( 76 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 372 )
  W( 131 ) = W( 131 ) + a*JVS( 373 )
  W( 136 ) = W( 136 ) + a*JVS( 374 )
  a = -W( 89 ) / JVS(          454  )
  W( 89 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 455 )
  W( 135 ) = W( 135 ) + a*JVS( 456 )
  W( 136 ) = W( 136 ) + a*JVS( 457 )
  a = -W( 102 ) / JVS(          564  )
  W( 102 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 565 )
  W( 132 ) = W( 132 ) + a*JVS( 566 )
  W( 133 ) = W( 133 ) + a*JVS( 567 )
  W( 135 ) = W( 135 ) + a*JVS( 568 )
  W( 136 ) = W( 136 ) + a*JVS( 569 )
  JVS( 570) = W( 76 )
  JVS( 571) = W( 89 )
  JVS( 572) = W( 102 )
  JVS( 573) = W( 103 )
  JVS( 574) = W( 111 )
  JVS( 575) = W( 129 )
  JVS( 576) = W( 131 )
  JVS( 577) = W( 132 )
  JVS( 578) = W( 133 )
  JVS( 579) = W( 135 )
  JVS( 580) = W( 136 )
  IF ( ABS(  JVS( 581 )) < TINY(a) ) THEN
         IER = 104                                     
         RETURN
  END IF
   W( 104 ) = JVS( 581 )
   W( 132 ) = JVS( 582 )
   W( 135 ) = JVS( 583 )
   W( 136 ) = JVS( 584 )
  JVS( 581) = W( 104 )
  JVS( 582) = W( 132 )
  JVS( 583) = W( 135 )
  JVS( 584) = W( 136 )
  IF ( ABS(  JVS( 601 )) < TINY(a) ) THEN
         IER = 105                                     
         RETURN
  END IF
   W( 46 ) = JVS( 585 )
   W( 48 ) = JVS( 586 )
   W( 65 ) = JVS( 587 )
   W( 73 ) = JVS( 588 )
   W( 76 ) = JVS( 589 )
   W( 79 ) = JVS( 590 )
   W( 82 ) = JVS( 591 )
   W( 83 ) = JVS( 592 )
   W( 84 ) = JVS( 593 )
   W( 89 ) = JVS( 594 )
   W( 93 ) = JVS( 595 )
   W( 100 ) = JVS( 596 )
   W( 101 ) = JVS( 597 )
   W( 102 ) = JVS( 598 )
   W( 103 ) = JVS( 599 )
   W( 104 ) = JVS( 600 )
   W( 105 ) = JVS( 601 )
   W( 111 ) = JVS( 602 )
   W( 113 ) = JVS( 603 )
   W( 114 ) = JVS( 604 )
   W( 118 ) = JVS( 605 )
   W( 120 ) = JVS( 606 )
   W( 129 ) = JVS( 607 )
   W( 131 ) = JVS( 608 )
   W( 132 ) = JVS( 609 )
   W( 133 ) = JVS( 610 )
   W( 135 ) = JVS( 611 )
   W( 136 ) = JVS( 612 )
  a = -W( 46 ) / JVS(          207  )
  W( 46 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 208 )
  W( 131 ) = W( 131 ) + a*JVS( 209 )
  W( 136 ) = W( 136 ) + a*JVS( 210 )
  a = -W( 48 ) / JVS(          216  )
  W( 48 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 217 )
  W( 136 ) = W( 136 ) + a*JVS( 218 )
  a = -W( 65 ) / JVS(          309  )
  W( 65 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 310 )
  W( 131 ) = W( 131 ) + a*JVS( 311 )
  W( 136 ) = W( 136 ) + a*JVS( 312 )
  a = -W( 73 ) / JVS(          354  )
  W( 73 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 355 )
  W( 131 ) = W( 131 ) + a*JVS( 356 )
  W( 133 ) = W( 133 ) + a*JVS( 357 )
  W( 136 ) = W( 136 ) + a*JVS( 358 )
  a = -W( 76 ) / JVS(          371  )
  W( 76 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 372 )
  W( 131 ) = W( 131 ) + a*JVS( 373 )
  W( 136 ) = W( 136 ) + a*JVS( 374 )
  a = -W( 79 ) / JVS(          394  )
  W( 79 ) = -a
  W( 120 ) = W( 120 ) + a*JVS( 395 )
  W( 131 ) = W( 131 ) + a*JVS( 396 )
  W( 136 ) = W( 136 ) + a*JVS( 397 )
  a = -W( 82 ) / JVS(          415  )
  W( 82 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 416 )
  W( 136 ) = W( 136 ) + a*JVS( 417 )
  a = -W( 83 ) / JVS(          418  )
  W( 83 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 84 ) / JVS(          421  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 422 )
  W( 131 ) = W( 131 ) + a*JVS( 423 )
  W( 136 ) = W( 136 ) + a*JVS( 424 )
  a = -W( 89 ) / JVS(          454  )
  W( 89 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 455 )
  W( 135 ) = W( 135 ) + a*JVS( 456 )
  W( 136 ) = W( 136 ) + a*JVS( 457 )
  a = -W( 93 ) / JVS(          482  )
  W( 93 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 483 )
  W( 131 ) = W( 131 ) + a*JVS( 484 )
  W( 133 ) = W( 133 ) + a*JVS( 485 )
  W( 136 ) = W( 136 ) + a*JVS( 486 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 102 ) / JVS(          564  )
  W( 102 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 565 )
  W( 132 ) = W( 132 ) + a*JVS( 566 )
  W( 133 ) = W( 133 ) + a*JVS( 567 )
  W( 135 ) = W( 135 ) + a*JVS( 568 )
  W( 136 ) = W( 136 ) + a*JVS( 569 )
  a = -W( 103 ) / JVS(          573  )
  W( 103 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 574 )
  W( 129 ) = W( 129 ) + a*JVS( 575 )
  W( 131 ) = W( 131 ) + a*JVS( 576 )
  W( 132 ) = W( 132 ) + a*JVS( 577 )
  W( 133 ) = W( 133 ) + a*JVS( 578 )
  W( 135 ) = W( 135 ) + a*JVS( 579 )
  W( 136 ) = W( 136 ) + a*JVS( 580 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  JVS( 585) = W( 46 )
  JVS( 586) = W( 48 )
  JVS( 587) = W( 65 )
  JVS( 588) = W( 73 )
  JVS( 589) = W( 76 )
  JVS( 590) = W( 79 )
  JVS( 591) = W( 82 )
  JVS( 592) = W( 83 )
  JVS( 593) = W( 84 )
  JVS( 594) = W( 89 )
  JVS( 595) = W( 93 )
  JVS( 596) = W( 100 )
  JVS( 597) = W( 101 )
  JVS( 598) = W( 102 )
  JVS( 599) = W( 103 )
  JVS( 600) = W( 104 )
  JVS( 601) = W( 105 )
  JVS( 602) = W( 111 )
  JVS( 603) = W( 113 )
  JVS( 604) = W( 114 )
  JVS( 605) = W( 118 )
  JVS( 606) = W( 120 )
  JVS( 607) = W( 129 )
  JVS( 608) = W( 131 )
  JVS( 609) = W( 132 )
  JVS( 610) = W( 133 )
  JVS( 611) = W( 135 )
  JVS( 612) = W( 136 )
  IF ( ABS(  JVS( 613 )) < TINY(a) ) THEN
         IER = 106                                     
         RETURN
  END IF
   W( 106 ) = JVS( 613 )
   W( 132 ) = JVS( 614 )
   W( 135 ) = JVS( 615 )
   W( 136 ) = JVS( 616 )
  JVS( 613) = W( 106 )
  JVS( 614) = W( 132 )
  JVS( 615) = W( 135 )
  JVS( 616) = W( 136 )
  IF ( ABS(  JVS( 640 )) < TINY(a) ) THEN
         IER = 107                                     
         RETURN
  END IF
   W( 19 ) = JVS( 617 )
   W( 21 ) = JVS( 618 )
   W( 31 ) = JVS( 619 )
   W( 37 ) = JVS( 620 )
   W( 40 ) = JVS( 621 )
   W( 42 ) = JVS( 622 )
   W( 55 ) = JVS( 623 )
   W( 79 ) = JVS( 624 )
   W( 82 ) = JVS( 625 )
   W( 83 ) = JVS( 626 )
   W( 87 ) = JVS( 627 )
   W( 88 ) = JVS( 628 )
   W( 89 ) = JVS( 629 )
   W( 90 ) = JVS( 630 )
   W( 91 ) = JVS( 631 )
   W( 94 ) = JVS( 632 )
   W( 96 ) = JVS( 633 )
   W( 98 ) = JVS( 634 )
   W( 100 ) = JVS( 635 )
   W( 101 ) = JVS( 636 )
   W( 103 ) = JVS( 637 )
   W( 104 ) = JVS( 638 )
   W( 106 ) = JVS( 639 )
   W( 107 ) = JVS( 640 )
   W( 108 ) = JVS( 641 )
   W( 109 ) = JVS( 642 )
   W( 110 ) = JVS( 643 )
   W( 111 ) = JVS( 644 )
   W( 112 ) = JVS( 645 )
   W( 114 ) = JVS( 646 )
   W( 115 ) = JVS( 647 )
   W( 116 ) = JVS( 648 )
   W( 117 ) = JVS( 649 )
   W( 120 ) = JVS( 650 )
   W( 123 ) = JVS( 651 )
   W( 124 ) = JVS( 652 )
   W( 125 ) = JVS( 653 )
   W( 126 ) = JVS( 654 )
   W( 127 ) = JVS( 655 )
   W( 128 ) = JVS( 656 )
   W( 129 ) = JVS( 657 )
   W( 131 ) = JVS( 658 )
   W( 132 ) = JVS( 659 )
   W( 133 ) = JVS( 660 )
   W( 134 ) = JVS( 661 )
   W( 135 ) = JVS( 662 )
   W( 136 ) = JVS( 663 )
  a = -W( 19 ) / JVS(          111  )
  W( 19 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 112 )
  a = -W( 21 ) / JVS(          116  )
  W( 21 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 117 )
  a = -W( 31 ) / JVS(          150  )
  W( 31 ) = -a
  W( 37 ) = W( 37 ) + a*JVS( 151 )
  W( 136 ) = W( 136 ) + a*JVS( 152 )
  a = -W( 37 ) / JVS(          172  )
  W( 37 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 173 )
  a = -W( 40 ) / JVS(          180  )
  W( 40 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 181 )
  W( 133 ) = W( 133 ) + a*JVS( 182 )
  a = -W( 42 ) / JVS(          187  )
  W( 42 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 188 )
  W( 90 ) = W( 90 ) + a*JVS( 189 )
  W( 133 ) = W( 133 ) + a*JVS( 190 )
  a = -W( 55 ) / JVS(          251  )
  W( 55 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 252 )
  W( 136 ) = W( 136 ) + a*JVS( 253 )
  a = -W( 79 ) / JVS(          394  )
  W( 79 ) = -a
  W( 120 ) = W( 120 ) + a*JVS( 395 )
  W( 131 ) = W( 131 ) + a*JVS( 396 )
  W( 136 ) = W( 136 ) + a*JVS( 397 )
  a = -W( 82 ) / JVS(          415  )
  W( 82 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 416 )
  W( 136 ) = W( 136 ) + a*JVS( 417 )
  a = -W( 83 ) / JVS(          418  )
  W( 83 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 87 ) / JVS(          441  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 90 ) = W( 90 ) + a*JVS( 443 )
  W( 131 ) = W( 131 ) + a*JVS( 444 )
  W( 133 ) = W( 133 ) + a*JVS( 445 )
  W( 134 ) = W( 134 ) + a*JVS( 446 )
  W( 136 ) = W( 136 ) + a*JVS( 447 )
  a = -W( 88 ) / JVS(          450  )
  W( 88 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 451 )
  W( 133 ) = W( 133 ) + a*JVS( 452 )
  W( 136 ) = W( 136 ) + a*JVS( 453 )
  a = -W( 89 ) / JVS(          454  )
  W( 89 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 455 )
  W( 135 ) = W( 135 ) + a*JVS( 456 )
  W( 136 ) = W( 136 ) + a*JVS( 457 )
  a = -W( 90 ) / JVS(          460  )
  W( 90 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 461 )
  W( 133 ) = W( 133 ) + a*JVS( 462 )
  W( 136 ) = W( 136 ) + a*JVS( 463 )
  a = -W( 91 ) / JVS(          467  )
  W( 91 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 468 )
  W( 133 ) = W( 133 ) + a*JVS( 469 )
  W( 134 ) = W( 134 ) + a*JVS( 470 )
  W( 136 ) = W( 136 ) + a*JVS( 471 )
  a = -W( 94 ) / JVS(          490  )
  W( 94 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 491 )
  W( 133 ) = W( 133 ) + a*JVS( 492 )
  W( 134 ) = W( 134 ) + a*JVS( 493 )
  W( 136 ) = W( 136 ) + a*JVS( 494 )
  a = -W( 96 ) / JVS(          516  )
  W( 96 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 517 )
  W( 108 ) = W( 108 ) + a*JVS( 518 )
  W( 123 ) = W( 123 ) + a*JVS( 519 )
  W( 127 ) = W( 127 ) + a*JVS( 520 )
  W( 128 ) = W( 128 ) + a*JVS( 521 )
  W( 129 ) = W( 129 ) + a*JVS( 522 )
  W( 131 ) = W( 131 ) + a*JVS( 523 )
  W( 133 ) = W( 133 ) + a*JVS( 524 )
  W( 134 ) = W( 134 ) + a*JVS( 525 )
  W( 135 ) = W( 135 ) + a*JVS( 526 )
  W( 136 ) = W( 136 ) + a*JVS( 527 )
  a = -W( 98 ) / JVS(          542  )
  W( 98 ) = -a
  W( 100 ) = W( 100 ) + a*JVS( 543 )
  W( 101 ) = W( 101 ) + a*JVS( 544 )
  W( 104 ) = W( 104 ) + a*JVS( 545 )
  W( 131 ) = W( 131 ) + a*JVS( 546 )
  W( 132 ) = W( 132 ) + a*JVS( 547 )
  W( 133 ) = W( 133 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 103 ) / JVS(          573  )
  W( 103 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 574 )
  W( 129 ) = W( 129 ) + a*JVS( 575 )
  W( 131 ) = W( 131 ) + a*JVS( 576 )
  W( 132 ) = W( 132 ) + a*JVS( 577 )
  W( 133 ) = W( 133 ) + a*JVS( 578 )
  W( 135 ) = W( 135 ) + a*JVS( 579 )
  W( 136 ) = W( 136 ) + a*JVS( 580 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  JVS( 617) = W( 19 )
  JVS( 618) = W( 21 )
  JVS( 619) = W( 31 )
  JVS( 620) = W( 37 )
  JVS( 621) = W( 40 )
  JVS( 622) = W( 42 )
  JVS( 623) = W( 55 )
  JVS( 624) = W( 79 )
  JVS( 625) = W( 82 )
  JVS( 626) = W( 83 )
  JVS( 627) = W( 87 )
  JVS( 628) = W( 88 )
  JVS( 629) = W( 89 )
  JVS( 630) = W( 90 )
  JVS( 631) = W( 91 )
  JVS( 632) = W( 94 )
  JVS( 633) = W( 96 )
  JVS( 634) = W( 98 )
  JVS( 635) = W( 100 )
  JVS( 636) = W( 101 )
  JVS( 637) = W( 103 )
  JVS( 638) = W( 104 )
  JVS( 639) = W( 106 )
  JVS( 640) = W( 107 )
  JVS( 641) = W( 108 )
  JVS( 642) = W( 109 )
  JVS( 643) = W( 110 )
  JVS( 644) = W( 111 )
  JVS( 645) = W( 112 )
  JVS( 646) = W( 114 )
  JVS( 647) = W( 115 )
  JVS( 648) = W( 116 )
  JVS( 649) = W( 117 )
  JVS( 650) = W( 120 )
  JVS( 651) = W( 123 )
  JVS( 652) = W( 124 )
  JVS( 653) = W( 125 )
  JVS( 654) = W( 126 )
  JVS( 655) = W( 127 )
  JVS( 656) = W( 128 )
  JVS( 657) = W( 129 )
  JVS( 658) = W( 131 )
  JVS( 659) = W( 132 )
  JVS( 660) = W( 133 )
  JVS( 661) = W( 134 )
  JVS( 662) = W( 135 )
  JVS( 663) = W( 136 )
  IF ( ABS(  JVS( 668 )) < TINY(a) ) THEN
         IER = 108                                     
         RETURN
  END IF
   W( 23 ) = JVS( 664 )
   W( 60 ) = JVS( 665 )
   W( 76 ) = JVS( 666 )
   W( 79 ) = JVS( 667 )
   W( 108 ) = JVS( 668 )
   W( 111 ) = JVS( 669 )
   W( 120 ) = JVS( 670 )
   W( 123 ) = JVS( 671 )
   W( 125 ) = JVS( 672 )
   W( 127 ) = JVS( 673 )
   W( 128 ) = JVS( 674 )
   W( 129 ) = JVS( 675 )
   W( 131 ) = JVS( 676 )
   W( 132 ) = JVS( 677 )
   W( 133 ) = JVS( 678 )
   W( 135 ) = JVS( 679 )
   W( 136 ) = JVS( 680 )
  a = -W( 23 ) / JVS(          121  )
  W( 23 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 122 )
  W( 133 ) = W( 133 ) + a*JVS( 123 )
  a = -W( 60 ) / JVS(          279  )
  W( 60 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 280 )
  W( 133 ) = W( 133 ) + a*JVS( 281 )
  W( 136 ) = W( 136 ) + a*JVS( 282 )
  a = -W( 76 ) / JVS(          371  )
  W( 76 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 372 )
  W( 131 ) = W( 131 ) + a*JVS( 373 )
  W( 136 ) = W( 136 ) + a*JVS( 374 )
  a = -W( 79 ) / JVS(          394  )
  W( 79 ) = -a
  W( 120 ) = W( 120 ) + a*JVS( 395 )
  W( 131 ) = W( 131 ) + a*JVS( 396 )
  W( 136 ) = W( 136 ) + a*JVS( 397 )
  JVS( 664) = W( 23 )
  JVS( 665) = W( 60 )
  JVS( 666) = W( 76 )
  JVS( 667) = W( 79 )
  JVS( 668) = W( 108 )
  JVS( 669) = W( 111 )
  JVS( 670) = W( 120 )
  JVS( 671) = W( 123 )
  JVS( 672) = W( 125 )
  JVS( 673) = W( 127 )
  JVS( 674) = W( 128 )
  JVS( 675) = W( 129 )
  JVS( 676) = W( 131 )
  JVS( 677) = W( 132 )
  JVS( 678) = W( 133 )
  JVS( 679) = W( 135 )
  JVS( 680) = W( 136 )
  IF ( ABS(  JVS( 682 )) < TINY(a) ) THEN
         IER = 109                                     
         RETURN
  END IF
   W( 106 ) = JVS( 681 )
   W( 109 ) = JVS( 682 )
   W( 124 ) = JVS( 683 )
   W( 132 ) = JVS( 684 )
   W( 135 ) = JVS( 685 )
   W( 136 ) = JVS( 686 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  JVS( 681) = W( 106 )
  JVS( 682) = W( 109 )
  JVS( 683) = W( 124 )
  JVS( 684) = W( 132 )
  JVS( 685) = W( 135 )
  JVS( 686) = W( 136 )
  IF ( ABS(  JVS( 693 )) < TINY(a) ) THEN
         IER = 110                                     
         RETURN
  END IF
   W( 65 ) = JVS( 687 )
   W( 82 ) = JVS( 688 )
   W( 83 ) = JVS( 689 )
   W( 100 ) = JVS( 690 )
   W( 101 ) = JVS( 691 )
   W( 104 ) = JVS( 692 )
   W( 110 ) = JVS( 693 )
   W( 118 ) = JVS( 694 )
   W( 119 ) = JVS( 695 )
   W( 129 ) = JVS( 696 )
   W( 131 ) = JVS( 697 )
   W( 132 ) = JVS( 698 )
   W( 133 ) = JVS( 699 )
   W( 135 ) = JVS( 700 )
   W( 136 ) = JVS( 701 )
  a = -W( 65 ) / JVS(          309  )
  W( 65 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 310 )
  W( 131 ) = W( 131 ) + a*JVS( 311 )
  W( 136 ) = W( 136 ) + a*JVS( 312 )
  a = -W( 82 ) / JVS(          415  )
  W( 82 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 416 )
  W( 136 ) = W( 136 ) + a*JVS( 417 )
  a = -W( 83 ) / JVS(          418  )
  W( 83 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  JVS( 687) = W( 65 )
  JVS( 688) = W( 82 )
  JVS( 689) = W( 83 )
  JVS( 690) = W( 100 )
  JVS( 691) = W( 101 )
  JVS( 692) = W( 104 )
  JVS( 693) = W( 110 )
  JVS( 694) = W( 118 )
  JVS( 695) = W( 119 )
  JVS( 696) = W( 129 )
  JVS( 697) = W( 131 )
  JVS( 698) = W( 132 )
  JVS( 699) = W( 133 )
  JVS( 700) = W( 135 )
  JVS( 701) = W( 136 )
  IF ( ABS(  JVS( 704 )) < TINY(a) ) THEN
         IER = 111                                     
         RETURN
  END IF
   W( 76 ) = JVS( 702 )
   W( 89 ) = JVS( 703 )
   W( 111 ) = JVS( 704 )
   W( 129 ) = JVS( 705 )
   W( 131 ) = JVS( 706 )
   W( 132 ) = JVS( 707 )
   W( 133 ) = JVS( 708 )
   W( 135 ) = JVS( 709 )
   W( 136 ) = JVS( 710 )
  a = -W( 76 ) / JVS(          371  )
  W( 76 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 372 )
  W( 131 ) = W( 131 ) + a*JVS( 373 )
  W( 136 ) = W( 136 ) + a*JVS( 374 )
  a = -W( 89 ) / JVS(          454  )
  W( 89 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 455 )
  W( 135 ) = W( 135 ) + a*JVS( 456 )
  W( 136 ) = W( 136 ) + a*JVS( 457 )
  JVS( 702) = W( 76 )
  JVS( 703) = W( 89 )
  JVS( 704) = W( 111 )
  JVS( 705) = W( 129 )
  JVS( 706) = W( 131 )
  JVS( 707) = W( 132 )
  JVS( 708) = W( 133 )
  JVS( 709) = W( 135 )
  JVS( 710) = W( 136 )
  IF ( ABS(  JVS( 722 )) < TINY(a) ) THEN
         IER = 112                                     
         RETURN
  END IF
   W( 41 ) = JVS( 711 )
   W( 45 ) = JVS( 712 )
   W( 48 ) = JVS( 713 )
   W( 49 ) = JVS( 714 )
   W( 70 ) = JVS( 715 )
   W( 84 ) = JVS( 716 )
   W( 85 ) = JVS( 717 )
   W( 86 ) = JVS( 718 )
   W( 93 ) = JVS( 719 )
   W( 97 ) = JVS( 720 )
   W( 109 ) = JVS( 721 )
   W( 112 ) = JVS( 722 )
   W( 113 ) = JVS( 723 )
   W( 124 ) = JVS( 724 )
   W( 129 ) = JVS( 725 )
   W( 131 ) = JVS( 726 )
   W( 132 ) = JVS( 727 )
   W( 133 ) = JVS( 728 )
   W( 135 ) = JVS( 729 )
   W( 136 ) = JVS( 730 )
  a = -W( 41 ) / JVS(          183  )
  W( 41 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 184 )
  W( 129 ) = W( 129 ) + a*JVS( 185 )
  W( 136 ) = W( 136 ) + a*JVS( 186 )
  a = -W( 45 ) / JVS(          203  )
  W( 45 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 204 )
  W( 131 ) = W( 131 ) + a*JVS( 205 )
  W( 136 ) = W( 136 ) + a*JVS( 206 )
  a = -W( 48 ) / JVS(          216  )
  W( 48 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 217 )
  W( 136 ) = W( 136 ) + a*JVS( 218 )
  a = -W( 49 ) / JVS(          219  )
  W( 49 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 220 )
  W( 131 ) = W( 131 ) + a*JVS( 221 )
  W( 136 ) = W( 136 ) + a*JVS( 222 )
  a = -W( 70 ) / JVS(          336  )
  W( 70 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 337 )
  W( 131 ) = W( 131 ) + a*JVS( 338 )
  W( 136 ) = W( 136 ) + a*JVS( 339 )
  a = -W( 84 ) / JVS(          421  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 422 )
  W( 131 ) = W( 131 ) + a*JVS( 423 )
  W( 136 ) = W( 136 ) + a*JVS( 424 )
  a = -W( 85 ) / JVS(          428  )
  W( 85 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 429 )
  W( 131 ) = W( 131 ) + a*JVS( 430 )
  W( 133 ) = W( 133 ) + a*JVS( 431 )
  W( 136 ) = W( 136 ) + a*JVS( 432 )
  a = -W( 86 ) / JVS(          434  )
  W( 86 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 435 )
  W( 131 ) = W( 131 ) + a*JVS( 436 )
  W( 133 ) = W( 133 ) + a*JVS( 437 )
  W( 136 ) = W( 136 ) + a*JVS( 438 )
  a = -W( 93 ) / JVS(          482  )
  W( 93 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 483 )
  W( 131 ) = W( 131 ) + a*JVS( 484 )
  W( 133 ) = W( 133 ) + a*JVS( 485 )
  W( 136 ) = W( 136 ) + a*JVS( 486 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 533 )
  W( 129 ) = W( 129 ) + a*JVS( 534 )
  W( 131 ) = W( 131 ) + a*JVS( 535 )
  W( 133 ) = W( 133 ) + a*JVS( 536 )
  W( 136 ) = W( 136 ) + a*JVS( 537 )
  a = -W( 109 ) / JVS(          682  )
  W( 109 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 683 )
  W( 132 ) = W( 132 ) + a*JVS( 684 )
  W( 135 ) = W( 135 ) + a*JVS( 685 )
  W( 136 ) = W( 136 ) + a*JVS( 686 )
  JVS( 711) = W( 41 )
  JVS( 712) = W( 45 )
  JVS( 713) = W( 48 )
  JVS( 714) = W( 49 )
  JVS( 715) = W( 70 )
  JVS( 716) = W( 84 )
  JVS( 717) = W( 85 )
  JVS( 718) = W( 86 )
  JVS( 719) = W( 93 )
  JVS( 720) = W( 97 )
  JVS( 721) = W( 109 )
  JVS( 722) = W( 112 )
  JVS( 723) = W( 113 )
  JVS( 724) = W( 124 )
  JVS( 725) = W( 129 )
  JVS( 726) = W( 131 )
  JVS( 727) = W( 132 )
  JVS( 728) = W( 133 )
  JVS( 729) = W( 135 )
  JVS( 730) = W( 136 )
  IF ( ABS(  JVS( 738 )) < TINY(a) ) THEN
         IER = 113                                     
         RETURN
  END IF
   W( 54 ) = JVS( 731 )
   W( 82 ) = JVS( 732 )
   W( 83 ) = JVS( 733 )
   W( 84 ) = JVS( 734 )
   W( 100 ) = JVS( 735 )
   W( 101 ) = JVS( 736 )
   W( 104 ) = JVS( 737 )
   W( 113 ) = JVS( 738 )
   W( 131 ) = JVS( 739 )
   W( 132 ) = JVS( 740 )
   W( 133 ) = JVS( 741 )
   W( 135 ) = JVS( 742 )
   W( 136 ) = JVS( 743 )
  a = -W( 54 ) / JVS(          243  )
  W( 54 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 244 )
  W( 83 ) = W( 83 ) + a*JVS( 245 )
  W( 100 ) = W( 100 ) + a*JVS( 246 )
  W( 101 ) = W( 101 ) + a*JVS( 247 )
  W( 104 ) = W( 104 ) + a*JVS( 248 )
  W( 132 ) = W( 132 ) + a*JVS( 249 )
  W( 136 ) = W( 136 ) + a*JVS( 250 )
  a = -W( 82 ) / JVS(          415  )
  W( 82 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 416 )
  W( 136 ) = W( 136 ) + a*JVS( 417 )
  a = -W( 83 ) / JVS(          418  )
  W( 83 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 84 ) / JVS(          421  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 422 )
  W( 131 ) = W( 131 ) + a*JVS( 423 )
  W( 136 ) = W( 136 ) + a*JVS( 424 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  JVS( 731) = W( 54 )
  JVS( 732) = W( 82 )
  JVS( 733) = W( 83 )
  JVS( 734) = W( 84 )
  JVS( 735) = W( 100 )
  JVS( 736) = W( 101 )
  JVS( 737) = W( 104 )
  JVS( 738) = W( 113 )
  JVS( 739) = W( 131 )
  JVS( 740) = W( 132 )
  JVS( 741) = W( 133 )
  JVS( 742) = W( 135 )
  JVS( 743) = W( 136 )
  IF ( ABS(  JVS( 751 )) < TINY(a) ) THEN
         IER = 114                                     
         RETURN
  END IF
   W( 79 ) = JVS( 744 )
   W( 82 ) = JVS( 745 )
   W( 83 ) = JVS( 746 )
   W( 100 ) = JVS( 747 )
   W( 101 ) = JVS( 748 )
   W( 104 ) = JVS( 749 )
   W( 110 ) = JVS( 750 )
   W( 114 ) = JVS( 751 )
   W( 118 ) = JVS( 752 )
   W( 119 ) = JVS( 753 )
   W( 120 ) = JVS( 754 )
   W( 129 ) = JVS( 755 )
   W( 131 ) = JVS( 756 )
   W( 132 ) = JVS( 757 )
   W( 133 ) = JVS( 758 )
   W( 135 ) = JVS( 759 )
   W( 136 ) = JVS( 760 )
  a = -W( 79 ) / JVS(          394  )
  W( 79 ) = -a
  W( 120 ) = W( 120 ) + a*JVS( 395 )
  W( 131 ) = W( 131 ) + a*JVS( 396 )
  W( 136 ) = W( 136 ) + a*JVS( 397 )
  a = -W( 82 ) / JVS(          415  )
  W( 82 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 416 )
  W( 136 ) = W( 136 ) + a*JVS( 417 )
  a = -W( 83 ) / JVS(          418  )
  W( 83 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  a = -W( 110 ) / JVS(          693  )
  W( 110 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 694 )
  W( 119 ) = W( 119 ) + a*JVS( 695 )
  W( 129 ) = W( 129 ) + a*JVS( 696 )
  W( 131 ) = W( 131 ) + a*JVS( 697 )
  W( 132 ) = W( 132 ) + a*JVS( 698 )
  W( 133 ) = W( 133 ) + a*JVS( 699 )
  W( 135 ) = W( 135 ) + a*JVS( 700 )
  W( 136 ) = W( 136 ) + a*JVS( 701 )
  JVS( 744) = W( 79 )
  JVS( 745) = W( 82 )
  JVS( 746) = W( 83 )
  JVS( 747) = W( 100 )
  JVS( 748) = W( 101 )
  JVS( 749) = W( 104 )
  JVS( 750) = W( 110 )
  JVS( 751) = W( 114 )
  JVS( 752) = W( 118 )
  JVS( 753) = W( 119 )
  JVS( 754) = W( 120 )
  JVS( 755) = W( 129 )
  JVS( 756) = W( 131 )
  JVS( 757) = W( 132 )
  JVS( 758) = W( 133 )
  JVS( 759) = W( 135 )
  JVS( 760) = W( 136 )
  IF ( ABS(  JVS( 764 )) < TINY(a) ) THEN
         IER = 115                                     
         RETURN
  END IF
   W( 75 ) = JVS( 761 )
   W( 99 ) = JVS( 762 )
   W( 106 ) = JVS( 763 )
   W( 115 ) = JVS( 764 )
   W( 127 ) = JVS( 765 )
   W( 128 ) = JVS( 766 )
   W( 129 ) = JVS( 767 )
   W( 131 ) = JVS( 768 )
   W( 132 ) = JVS( 769 )
   W( 133 ) = JVS( 770 )
   W( 135 ) = JVS( 771 )
   W( 136 ) = JVS( 772 )
  a = -W( 75 ) / JVS(          367  )
  W( 75 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 368 )
  W( 131 ) = W( 131 ) + a*JVS( 369 )
  W( 136 ) = W( 136 ) + a*JVS( 370 )
  a = -W( 99 ) / JVS(          550  )
  W( 99 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 133 ) = W( 133 ) + a*JVS( 553 )
  W( 135 ) = W( 135 ) + a*JVS( 554 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  JVS( 761) = W( 75 )
  JVS( 762) = W( 99 )
  JVS( 763) = W( 106 )
  JVS( 764) = W( 115 )
  JVS( 765) = W( 127 )
  JVS( 766) = W( 128 )
  JVS( 767) = W( 129 )
  JVS( 768) = W( 131 )
  JVS( 769) = W( 132 )
  JVS( 770) = W( 133 )
  JVS( 771) = W( 135 )
  JVS( 772) = W( 136 )
  IF ( ABS(  JVS( 775 )) < TINY(a) ) THEN
         IER = 116                                     
         RETURN
  END IF
   W( 99 ) = JVS( 773 )
   W( 106 ) = JVS( 774 )
   W( 116 ) = JVS( 775 )
   W( 118 ) = JVS( 776 )
   W( 119 ) = JVS( 777 )
   W( 120 ) = JVS( 778 )
   W( 125 ) = JVS( 779 )
   W( 127 ) = JVS( 780 )
   W( 129 ) = JVS( 781 )
   W( 131 ) = JVS( 782 )
   W( 132 ) = JVS( 783 )
   W( 133 ) = JVS( 784 )
   W( 135 ) = JVS( 785 )
   W( 136 ) = JVS( 786 )
  a = -W( 99 ) / JVS(          550  )
  W( 99 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 133 ) = W( 133 ) + a*JVS( 553 )
  W( 135 ) = W( 135 ) + a*JVS( 554 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  JVS( 773) = W( 99 )
  JVS( 774) = W( 106 )
  JVS( 775) = W( 116 )
  JVS( 776) = W( 118 )
  JVS( 777) = W( 119 )
  JVS( 778) = W( 120 )
  JVS( 779) = W( 125 )
  JVS( 780) = W( 127 )
  JVS( 781) = W( 129 )
  JVS( 782) = W( 131 )
  JVS( 783) = W( 132 )
  JVS( 784) = W( 133 )
  JVS( 785) = W( 135 )
  JVS( 786) = W( 136 )
  IF ( ABS(  JVS( 803 )) < TINY(a) ) THEN
         IER = 117                                     
         RETURN
  END IF
   W( 40 ) = JVS( 787 )
   W( 57 ) = JVS( 788 )
   W( 67 ) = JVS( 789 )
   W( 69 ) = JVS( 790 )
   W( 88 ) = JVS( 791 )
   W( 90 ) = JVS( 792 )
   W( 91 ) = JVS( 793 )
   W( 94 ) = JVS( 794 )
   W( 98 ) = JVS( 795 )
   W( 100 ) = JVS( 796 )
   W( 101 ) = JVS( 797 )
   W( 102 ) = JVS( 798 )
   W( 104 ) = JVS( 799 )
   W( 109 ) = JVS( 800 )
   W( 113 ) = JVS( 801 )
   W( 115 ) = JVS( 802 )
   W( 117 ) = JVS( 803 )
   W( 121 ) = JVS( 804 )
   W( 122 ) = JVS( 805 )
   W( 123 ) = JVS( 806 )
   W( 124 ) = JVS( 807 )
   W( 125 ) = JVS( 808 )
   W( 127 ) = JVS( 809 )
   W( 128 ) = JVS( 810 )
   W( 129 ) = JVS( 811 )
   W( 131 ) = JVS( 812 )
   W( 132 ) = JVS( 813 )
   W( 133 ) = JVS( 814 )
   W( 134 ) = JVS( 815 )
   W( 135 ) = JVS( 816 )
   W( 136 ) = JVS( 817 )
  a = -W( 40 ) / JVS(          180  )
  W( 40 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 181 )
  W( 133 ) = W( 133 ) + a*JVS( 182 )
  a = -W( 57 ) / JVS(          258  )
  W( 57 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 259 )
  W( 131 ) = W( 131 ) + a*JVS( 260 )
  W( 136 ) = W( 136 ) + a*JVS( 261 )
  a = -W( 67 ) / JVS(          319  )
  W( 67 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 320 )
  W( 109 ) = W( 109 ) + a*JVS( 321 )
  W( 113 ) = W( 113 ) + a*JVS( 322 )
  W( 133 ) = W( 133 ) + a*JVS( 323 )
  W( 135 ) = W( 135 ) + a*JVS( 324 )
  W( 136 ) = W( 136 ) + a*JVS( 325 )
  a = -W( 69 ) / JVS(          332  )
  W( 69 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 333 )
  W( 133 ) = W( 133 ) + a*JVS( 334 )
  W( 136 ) = W( 136 ) + a*JVS( 335 )
  a = -W( 88 ) / JVS(          450  )
  W( 88 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 451 )
  W( 133 ) = W( 133 ) + a*JVS( 452 )
  W( 136 ) = W( 136 ) + a*JVS( 453 )
  a = -W( 90 ) / JVS(          460  )
  W( 90 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 461 )
  W( 133 ) = W( 133 ) + a*JVS( 462 )
  W( 136 ) = W( 136 ) + a*JVS( 463 )
  a = -W( 91 ) / JVS(          467  )
  W( 91 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 468 )
  W( 133 ) = W( 133 ) + a*JVS( 469 )
  W( 134 ) = W( 134 ) + a*JVS( 470 )
  W( 136 ) = W( 136 ) + a*JVS( 471 )
  a = -W( 94 ) / JVS(          490  )
  W( 94 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 491 )
  W( 133 ) = W( 133 ) + a*JVS( 492 )
  W( 134 ) = W( 134 ) + a*JVS( 493 )
  W( 136 ) = W( 136 ) + a*JVS( 494 )
  a = -W( 98 ) / JVS(          542  )
  W( 98 ) = -a
  W( 100 ) = W( 100 ) + a*JVS( 543 )
  W( 101 ) = W( 101 ) + a*JVS( 544 )
  W( 104 ) = W( 104 ) + a*JVS( 545 )
  W( 131 ) = W( 131 ) + a*JVS( 546 )
  W( 132 ) = W( 132 ) + a*JVS( 547 )
  W( 133 ) = W( 133 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 102 ) / JVS(          564  )
  W( 102 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 565 )
  W( 132 ) = W( 132 ) + a*JVS( 566 )
  W( 133 ) = W( 133 ) + a*JVS( 567 )
  W( 135 ) = W( 135 ) + a*JVS( 568 )
  W( 136 ) = W( 136 ) + a*JVS( 569 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  a = -W( 109 ) / JVS(          682  )
  W( 109 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 683 )
  W( 132 ) = W( 132 ) + a*JVS( 684 )
  W( 135 ) = W( 135 ) + a*JVS( 685 )
  W( 136 ) = W( 136 ) + a*JVS( 686 )
  a = -W( 113 ) / JVS(          738  )
  W( 113 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 739 )
  W( 132 ) = W( 132 ) + a*JVS( 740 )
  W( 133 ) = W( 133 ) + a*JVS( 741 )
  W( 135 ) = W( 135 ) + a*JVS( 742 )
  W( 136 ) = W( 136 ) + a*JVS( 743 )
  a = -W( 115 ) / JVS(          764  )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 765 )
  W( 128 ) = W( 128 ) + a*JVS( 766 )
  W( 129 ) = W( 129 ) + a*JVS( 767 )
  W( 131 ) = W( 131 ) + a*JVS( 768 )
  W( 132 ) = W( 132 ) + a*JVS( 769 )
  W( 133 ) = W( 133 ) + a*JVS( 770 )
  W( 135 ) = W( 135 ) + a*JVS( 771 )
  W( 136 ) = W( 136 ) + a*JVS( 772 )
  JVS( 787) = W( 40 )
  JVS( 788) = W( 57 )
  JVS( 789) = W( 67 )
  JVS( 790) = W( 69 )
  JVS( 791) = W( 88 )
  JVS( 792) = W( 90 )
  JVS( 793) = W( 91 )
  JVS( 794) = W( 94 )
  JVS( 795) = W( 98 )
  JVS( 796) = W( 100 )
  JVS( 797) = W( 101 )
  JVS( 798) = W( 102 )
  JVS( 799) = W( 104 )
  JVS( 800) = W( 109 )
  JVS( 801) = W( 113 )
  JVS( 802) = W( 115 )
  JVS( 803) = W( 117 )
  JVS( 804) = W( 121 )
  JVS( 805) = W( 122 )
  JVS( 806) = W( 123 )
  JVS( 807) = W( 124 )
  JVS( 808) = W( 125 )
  JVS( 809) = W( 127 )
  JVS( 810) = W( 128 )
  JVS( 811) = W( 129 )
  JVS( 812) = W( 131 )
  JVS( 813) = W( 132 )
  JVS( 814) = W( 133 )
  JVS( 815) = W( 134 )
  JVS( 816) = W( 135 )
  JVS( 817) = W( 136 )
  IF ( ABS(  JVS( 824 )) < TINY(a) ) THEN
         IER = 118                                     
         RETURN
  END IF
   W( 65 ) = JVS( 818 )
   W( 82 ) = JVS( 819 )
   W( 83 ) = JVS( 820 )
   W( 100 ) = JVS( 821 )
   W( 101 ) = JVS( 822 )
   W( 104 ) = JVS( 823 )
   W( 118 ) = JVS( 824 )
   W( 129 ) = JVS( 825 )
   W( 131 ) = JVS( 826 )
   W( 132 ) = JVS( 827 )
   W( 133 ) = JVS( 828 )
   W( 135 ) = JVS( 829 )
   W( 136 ) = JVS( 830 )
  a = -W( 65 ) / JVS(          309  )
  W( 65 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 310 )
  W( 131 ) = W( 131 ) + a*JVS( 311 )
  W( 136 ) = W( 136 ) + a*JVS( 312 )
  a = -W( 82 ) / JVS(          415  )
  W( 82 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 416 )
  W( 136 ) = W( 136 ) + a*JVS( 417 )
  a = -W( 83 ) / JVS(          418  )
  W( 83 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  JVS( 818) = W( 65 )
  JVS( 819) = W( 82 )
  JVS( 820) = W( 83 )
  JVS( 821) = W( 100 )
  JVS( 822) = W( 101 )
  JVS( 823) = W( 104 )
  JVS( 824) = W( 118 )
  JVS( 825) = W( 129 )
  JVS( 826) = W( 131 )
  JVS( 827) = W( 132 )
  JVS( 828) = W( 133 )
  JVS( 829) = W( 135 )
  JVS( 830) = W( 136 )
  IF ( ABS(  JVS( 836 )) < TINY(a) ) THEN
         IER = 119                                     
         RETURN
  END IF
   W( 100 ) = JVS( 831 )
   W( 101 ) = JVS( 832 )
   W( 104 ) = JVS( 833 )
   W( 110 ) = JVS( 834 )
   W( 118 ) = JVS( 835 )
   W( 119 ) = JVS( 836 )
   W( 129 ) = JVS( 837 )
   W( 131 ) = JVS( 838 )
   W( 132 ) = JVS( 839 )
   W( 133 ) = JVS( 840 )
   W( 135 ) = JVS( 841 )
   W( 136 ) = JVS( 842 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  a = -W( 110 ) / JVS(          693  )
  W( 110 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 694 )
  W( 119 ) = W( 119 ) + a*JVS( 695 )
  W( 129 ) = W( 129 ) + a*JVS( 696 )
  W( 131 ) = W( 131 ) + a*JVS( 697 )
  W( 132 ) = W( 132 ) + a*JVS( 698 )
  W( 133 ) = W( 133 ) + a*JVS( 699 )
  W( 135 ) = W( 135 ) + a*JVS( 700 )
  W( 136 ) = W( 136 ) + a*JVS( 701 )
  a = -W( 118 ) / JVS(          824  )
  W( 118 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 825 )
  W( 131 ) = W( 131 ) + a*JVS( 826 )
  W( 132 ) = W( 132 ) + a*JVS( 827 )
  W( 133 ) = W( 133 ) + a*JVS( 828 )
  W( 135 ) = W( 135 ) + a*JVS( 829 )
  W( 136 ) = W( 136 ) + a*JVS( 830 )
  JVS( 831) = W( 100 )
  JVS( 832) = W( 101 )
  JVS( 833) = W( 104 )
  JVS( 834) = W( 110 )
  JVS( 835) = W( 118 )
  JVS( 836) = W( 119 )
  JVS( 837) = W( 129 )
  JVS( 838) = W( 131 )
  JVS( 839) = W( 132 )
  JVS( 840) = W( 133 )
  JVS( 841) = W( 135 )
  JVS( 842) = W( 136 )
  IF ( ABS(  JVS( 847 )) < TINY(a) ) THEN
         IER = 120                                     
         RETURN
  END IF
   W( 79 ) = JVS( 843 )
   W( 110 ) = JVS( 844 )
   W( 118 ) = JVS( 845 )
   W( 119 ) = JVS( 846 )
   W( 120 ) = JVS( 847 )
   W( 129 ) = JVS( 848 )
   W( 131 ) = JVS( 849 )
   W( 132 ) = JVS( 850 )
   W( 133 ) = JVS( 851 )
   W( 135 ) = JVS( 852 )
   W( 136 ) = JVS( 853 )
  a = -W( 79 ) / JVS(          394  )
  W( 79 ) = -a
  W( 120 ) = W( 120 ) + a*JVS( 395 )
  W( 131 ) = W( 131 ) + a*JVS( 396 )
  W( 136 ) = W( 136 ) + a*JVS( 397 )
  a = -W( 110 ) / JVS(          693  )
  W( 110 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 694 )
  W( 119 ) = W( 119 ) + a*JVS( 695 )
  W( 129 ) = W( 129 ) + a*JVS( 696 )
  W( 131 ) = W( 131 ) + a*JVS( 697 )
  W( 132 ) = W( 132 ) + a*JVS( 698 )
  W( 133 ) = W( 133 ) + a*JVS( 699 )
  W( 135 ) = W( 135 ) + a*JVS( 700 )
  W( 136 ) = W( 136 ) + a*JVS( 701 )
  a = -W( 118 ) / JVS(          824  )
  W( 118 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 825 )
  W( 131 ) = W( 131 ) + a*JVS( 826 )
  W( 132 ) = W( 132 ) + a*JVS( 827 )
  W( 133 ) = W( 133 ) + a*JVS( 828 )
  W( 135 ) = W( 135 ) + a*JVS( 829 )
  W( 136 ) = W( 136 ) + a*JVS( 830 )
  a = -W( 119 ) / JVS(          836  )
  W( 119 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 837 )
  W( 131 ) = W( 131 ) + a*JVS( 838 )
  W( 132 ) = W( 132 ) + a*JVS( 839 )
  W( 133 ) = W( 133 ) + a*JVS( 840 )
  W( 135 ) = W( 135 ) + a*JVS( 841 )
  W( 136 ) = W( 136 ) + a*JVS( 842 )
  JVS( 843) = W( 79 )
  JVS( 844) = W( 110 )
  JVS( 845) = W( 118 )
  JVS( 846) = W( 119 )
  JVS( 847) = W( 120 )
  JVS( 848) = W( 129 )
  JVS( 849) = W( 131 )
  JVS( 850) = W( 132 )
  JVS( 851) = W( 133 )
  JVS( 852) = W( 135 )
  JVS( 853) = W( 136 )
  IF ( ABS(  JVS( 867 )) < TINY(a) ) THEN
         IER = 121                                     
         RETURN
  END IF
   W( 47 ) = JVS( 854 )
   W( 82 ) = JVS( 855 )
   W( 83 ) = JVS( 856 )
   W( 100 ) = JVS( 857 )
   W( 101 ) = JVS( 858 )
   W( 104 ) = JVS( 859 )
   W( 105 ) = JVS( 860 )
   W( 111 ) = JVS( 861 )
   W( 113 ) = JVS( 862 )
   W( 114 ) = JVS( 863 )
   W( 118 ) = JVS( 864 )
   W( 119 ) = JVS( 865 )
   W( 120 ) = JVS( 866 )
   W( 121 ) = JVS( 867 )
   W( 129 ) = JVS( 868 )
   W( 131 ) = JVS( 869 )
   W( 132 ) = JVS( 870 )
   W( 133 ) = JVS( 871 )
   W( 135 ) = JVS( 872 )
   W( 136 ) = JVS( 873 )
  a = -W( 47 ) / JVS(          211  )
  W( 47 ) = -a
  W( 121 ) = W( 121 ) + a*JVS( 212 )
  W( 131 ) = W( 131 ) + a*JVS( 213 )
  W( 136 ) = W( 136 ) + a*JVS( 214 )
  a = -W( 82 ) / JVS(          415  )
  W( 82 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 416 )
  W( 136 ) = W( 136 ) + a*JVS( 417 )
  a = -W( 83 ) / JVS(          418  )
  W( 83 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  a = -W( 105 ) / JVS(          601  )
  W( 105 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 602 )
  W( 113 ) = W( 113 ) + a*JVS( 603 )
  W( 114 ) = W( 114 ) + a*JVS( 604 )
  W( 118 ) = W( 118 ) + a*JVS( 605 )
  W( 120 ) = W( 120 ) + a*JVS( 606 )
  W( 129 ) = W( 129 ) + a*JVS( 607 )
  W( 131 ) = W( 131 ) + a*JVS( 608 )
  W( 132 ) = W( 132 ) + a*JVS( 609 )
  W( 133 ) = W( 133 ) + a*JVS( 610 )
  W( 135 ) = W( 135 ) + a*JVS( 611 )
  W( 136 ) = W( 136 ) + a*JVS( 612 )
  a = -W( 111 ) / JVS(          704  )
  W( 111 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 705 )
  W( 131 ) = W( 131 ) + a*JVS( 706 )
  W( 132 ) = W( 132 ) + a*JVS( 707 )
  W( 133 ) = W( 133 ) + a*JVS( 708 )
  W( 135 ) = W( 135 ) + a*JVS( 709 )
  W( 136 ) = W( 136 ) + a*JVS( 710 )
  a = -W( 113 ) / JVS(          738  )
  W( 113 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 739 )
  W( 132 ) = W( 132 ) + a*JVS( 740 )
  W( 133 ) = W( 133 ) + a*JVS( 741 )
  W( 135 ) = W( 135 ) + a*JVS( 742 )
  W( 136 ) = W( 136 ) + a*JVS( 743 )
  a = -W( 114 ) / JVS(          751  )
  W( 114 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 752 )
  W( 119 ) = W( 119 ) + a*JVS( 753 )
  W( 120 ) = W( 120 ) + a*JVS( 754 )
  W( 129 ) = W( 129 ) + a*JVS( 755 )
  W( 131 ) = W( 131 ) + a*JVS( 756 )
  W( 132 ) = W( 132 ) + a*JVS( 757 )
  W( 133 ) = W( 133 ) + a*JVS( 758 )
  W( 135 ) = W( 135 ) + a*JVS( 759 )
  W( 136 ) = W( 136 ) + a*JVS( 760 )
  a = -W( 118 ) / JVS(          824  )
  W( 118 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 825 )
  W( 131 ) = W( 131 ) + a*JVS( 826 )
  W( 132 ) = W( 132 ) + a*JVS( 827 )
  W( 133 ) = W( 133 ) + a*JVS( 828 )
  W( 135 ) = W( 135 ) + a*JVS( 829 )
  W( 136 ) = W( 136 ) + a*JVS( 830 )
  a = -W( 119 ) / JVS(          836  )
  W( 119 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 837 )
  W( 131 ) = W( 131 ) + a*JVS( 838 )
  W( 132 ) = W( 132 ) + a*JVS( 839 )
  W( 133 ) = W( 133 ) + a*JVS( 840 )
  W( 135 ) = W( 135 ) + a*JVS( 841 )
  W( 136 ) = W( 136 ) + a*JVS( 842 )
  a = -W( 120 ) / JVS(          847  )
  W( 120 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 848 )
  W( 131 ) = W( 131 ) + a*JVS( 849 )
  W( 132 ) = W( 132 ) + a*JVS( 850 )
  W( 133 ) = W( 133 ) + a*JVS( 851 )
  W( 135 ) = W( 135 ) + a*JVS( 852 )
  W( 136 ) = W( 136 ) + a*JVS( 853 )
  JVS( 854) = W( 47 )
  JVS( 855) = W( 82 )
  JVS( 856) = W( 83 )
  JVS( 857) = W( 100 )
  JVS( 858) = W( 101 )
  JVS( 859) = W( 104 )
  JVS( 860) = W( 105 )
  JVS( 861) = W( 111 )
  JVS( 862) = W( 113 )
  JVS( 863) = W( 114 )
  JVS( 864) = W( 118 )
  JVS( 865) = W( 119 )
  JVS( 866) = W( 120 )
  JVS( 867) = W( 121 )
  JVS( 868) = W( 129 )
  JVS( 869) = W( 131 )
  JVS( 870) = W( 132 )
  JVS( 871) = W( 133 )
  JVS( 872) = W( 135 )
  JVS( 873) = W( 136 )
  IF ( ABS(  JVS( 879 )) < TINY(a) ) THEN
         IER = 122                                     
         RETURN
  END IF
   W( 70 ) = JVS( 874 )
   W( 77 ) = JVS( 875 )
   W( 86 ) = JVS( 876 )
   W( 109 ) = JVS( 877 )
   W( 121 ) = JVS( 878 )
   W( 122 ) = JVS( 879 )
   W( 123 ) = JVS( 880 )
   W( 124 ) = JVS( 881 )
   W( 125 ) = JVS( 882 )
   W( 127 ) = JVS( 883 )
   W( 128 ) = JVS( 884 )
   W( 129 ) = JVS( 885 )
   W( 130 ) = JVS( 886 )
   W( 131 ) = JVS( 887 )
   W( 132 ) = JVS( 888 )
   W( 133 ) = JVS( 889 )
   W( 134 ) = JVS( 890 )
   W( 135 ) = JVS( 891 )
   W( 136 ) = JVS( 892 )
  a = -W( 70 ) / JVS(          336  )
  W( 70 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 337 )
  W( 131 ) = W( 131 ) + a*JVS( 338 )
  W( 136 ) = W( 136 ) + a*JVS( 339 )
  a = -W( 77 ) / JVS(          375  )
  W( 77 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 376 )
  W( 134 ) = W( 134 ) + a*JVS( 377 )
  W( 136 ) = W( 136 ) + a*JVS( 378 )
  a = -W( 86 ) / JVS(          434  )
  W( 86 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 435 )
  W( 131 ) = W( 131 ) + a*JVS( 436 )
  W( 133 ) = W( 133 ) + a*JVS( 437 )
  W( 136 ) = W( 136 ) + a*JVS( 438 )
  a = -W( 109 ) / JVS(          682  )
  W( 109 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 683 )
  W( 132 ) = W( 132 ) + a*JVS( 684 )
  W( 135 ) = W( 135 ) + a*JVS( 685 )
  W( 136 ) = W( 136 ) + a*JVS( 686 )
  a = -W( 121 ) / JVS(          867  )
  W( 121 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 868 )
  W( 131 ) = W( 131 ) + a*JVS( 869 )
  W( 132 ) = W( 132 ) + a*JVS( 870 )
  W( 133 ) = W( 133 ) + a*JVS( 871 )
  W( 135 ) = W( 135 ) + a*JVS( 872 )
  W( 136 ) = W( 136 ) + a*JVS( 873 )
  JVS( 874) = W( 70 )
  JVS( 875) = W( 77 )
  JVS( 876) = W( 86 )
  JVS( 877) = W( 109 )
  JVS( 878) = W( 121 )
  JVS( 879) = W( 122 )
  JVS( 880) = W( 123 )
  JVS( 881) = W( 124 )
  JVS( 882) = W( 125 )
  JVS( 883) = W( 127 )
  JVS( 884) = W( 128 )
  JVS( 885) = W( 129 )
  JVS( 886) = W( 130 )
  JVS( 887) = W( 131 )
  JVS( 888) = W( 132 )
  JVS( 889) = W( 133 )
  JVS( 890) = W( 134 )
  JVS( 891) = W( 135 )
  JVS( 892) = W( 136 )
  IF ( ABS(  JVS( 900 )) < TINY(a) ) THEN
         IER = 123                                     
         RETURN
  END IF
   W( 26 ) = JVS( 893 )
   W( 53 ) = JVS( 894 )
   W( 75 ) = JVS( 895 )
   W( 116 ) = JVS( 896 )
   W( 118 ) = JVS( 897 )
   W( 119 ) = JVS( 898 )
   W( 120 ) = JVS( 899 )
   W( 123 ) = JVS( 900 )
   W( 125 ) = JVS( 901 )
   W( 127 ) = JVS( 902 )
   W( 128 ) = JVS( 903 )
   W( 129 ) = JVS( 904 )
   W( 131 ) = JVS( 905 )
   W( 132 ) = JVS( 906 )
   W( 133 ) = JVS( 907 )
   W( 135 ) = JVS( 908 )
   W( 136 ) = JVS( 909 )
  a = -W( 26 ) / JVS(          129  )
  W( 26 ) = -a
  W( 123 ) = W( 123 ) + a*JVS( 130 )
  W( 131 ) = W( 131 ) + a*JVS( 131 )
  W( 136 ) = W( 136 ) + a*JVS( 132 )
  a = -W( 53 ) / JVS(          235  )
  W( 53 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 236 )
  W( 127 ) = W( 127 ) + a*JVS( 237 )
  W( 128 ) = W( 128 ) + a*JVS( 238 )
  W( 129 ) = W( 129 ) + a*JVS( 239 )
  W( 133 ) = W( 133 ) + a*JVS( 240 )
  W( 135 ) = W( 135 ) + a*JVS( 241 )
  W( 136 ) = W( 136 ) + a*JVS( 242 )
  a = -W( 75 ) / JVS(          367  )
  W( 75 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 368 )
  W( 131 ) = W( 131 ) + a*JVS( 369 )
  W( 136 ) = W( 136 ) + a*JVS( 370 )
  a = -W( 116 ) / JVS(          775  )
  W( 116 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 776 )
  W( 119 ) = W( 119 ) + a*JVS( 777 )
  W( 120 ) = W( 120 ) + a*JVS( 778 )
  W( 125 ) = W( 125 ) + a*JVS( 779 )
  W( 127 ) = W( 127 ) + a*JVS( 780 )
  W( 129 ) = W( 129 ) + a*JVS( 781 )
  W( 131 ) = W( 131 ) + a*JVS( 782 )
  W( 132 ) = W( 132 ) + a*JVS( 783 )
  W( 133 ) = W( 133 ) + a*JVS( 784 )
  W( 135 ) = W( 135 ) + a*JVS( 785 )
  W( 136 ) = W( 136 ) + a*JVS( 786 )
  a = -W( 118 ) / JVS(          824  )
  W( 118 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 825 )
  W( 131 ) = W( 131 ) + a*JVS( 826 )
  W( 132 ) = W( 132 ) + a*JVS( 827 )
  W( 133 ) = W( 133 ) + a*JVS( 828 )
  W( 135 ) = W( 135 ) + a*JVS( 829 )
  W( 136 ) = W( 136 ) + a*JVS( 830 )
  a = -W( 119 ) / JVS(          836  )
  W( 119 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 837 )
  W( 131 ) = W( 131 ) + a*JVS( 838 )
  W( 132 ) = W( 132 ) + a*JVS( 839 )
  W( 133 ) = W( 133 ) + a*JVS( 840 )
  W( 135 ) = W( 135 ) + a*JVS( 841 )
  W( 136 ) = W( 136 ) + a*JVS( 842 )
  a = -W( 120 ) / JVS(          847  )
  W( 120 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 848 )
  W( 131 ) = W( 131 ) + a*JVS( 849 )
  W( 132 ) = W( 132 ) + a*JVS( 850 )
  W( 133 ) = W( 133 ) + a*JVS( 851 )
  W( 135 ) = W( 135 ) + a*JVS( 852 )
  W( 136 ) = W( 136 ) + a*JVS( 853 )
  JVS( 893) = W( 26 )
  JVS( 894) = W( 53 )
  JVS( 895) = W( 75 )
  JVS( 896) = W( 116 )
  JVS( 897) = W( 118 )
  JVS( 898) = W( 119 )
  JVS( 899) = W( 120 )
  JVS( 900) = W( 123 )
  JVS( 901) = W( 125 )
  JVS( 902) = W( 127 )
  JVS( 903) = W( 128 )
  JVS( 904) = W( 129 )
  JVS( 905) = W( 131 )
  JVS( 906) = W( 132 )
  JVS( 907) = W( 133 )
  JVS( 908) = W( 135 )
  JVS( 909) = W( 136 )
  IF ( ABS(  JVS( 913 )) < TINY(a) ) THEN
         IER = 124                                     
         RETURN
  END IF
   W( 75 ) = JVS( 910 )
   W( 99 ) = JVS( 911 )
   W( 106 ) = JVS( 912 )
   W( 124 ) = JVS( 913 )
   W( 127 ) = JVS( 914 )
   W( 128 ) = JVS( 915 )
   W( 129 ) = JVS( 916 )
   W( 131 ) = JVS( 917 )
   W( 132 ) = JVS( 918 )
   W( 133 ) = JVS( 919 )
   W( 135 ) = JVS( 920 )
   W( 136 ) = JVS( 921 )
  a = -W( 75 ) / JVS(          367  )
  W( 75 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 368 )
  W( 131 ) = W( 131 ) + a*JVS( 369 )
  W( 136 ) = W( 136 ) + a*JVS( 370 )
  a = -W( 99 ) / JVS(          550  )
  W( 99 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 133 ) = W( 133 ) + a*JVS( 553 )
  W( 135 ) = W( 135 ) + a*JVS( 554 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  JVS( 910) = W( 75 )
  JVS( 911) = W( 99 )
  JVS( 912) = W( 106 )
  JVS( 913) = W( 124 )
  JVS( 914) = W( 127 )
  JVS( 915) = W( 128 )
  JVS( 916) = W( 129 )
  JVS( 917) = W( 131 )
  JVS( 918) = W( 132 )
  JVS( 919) = W( 133 )
  JVS( 920) = W( 135 )
  JVS( 921) = W( 136 )
  IF ( ABS(  JVS( 925 )) < TINY(a) ) THEN
         IER = 125                                     
         RETURN
  END IF
   W( 51 ) = JVS( 922 )
   W( 115 ) = JVS( 923 )
   W( 124 ) = JVS( 924 )
   W( 125 ) = JVS( 925 )
   W( 127 ) = JVS( 926 )
   W( 128 ) = JVS( 927 )
   W( 129 ) = JVS( 928 )
   W( 131 ) = JVS( 929 )
   W( 132 ) = JVS( 930 )
   W( 133 ) = JVS( 931 )
   W( 135 ) = JVS( 932 )
   W( 136 ) = JVS( 933 )
  a = -W( 51 ) / JVS(          227  )
  W( 51 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 228 )
  W( 131 ) = W( 131 ) + a*JVS( 229 )
  W( 136 ) = W( 136 ) + a*JVS( 230 )
  a = -W( 115 ) / JVS(          764  )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 765 )
  W( 128 ) = W( 128 ) + a*JVS( 766 )
  W( 129 ) = W( 129 ) + a*JVS( 767 )
  W( 131 ) = W( 131 ) + a*JVS( 768 )
  W( 132 ) = W( 132 ) + a*JVS( 769 )
  W( 133 ) = W( 133 ) + a*JVS( 770 )
  W( 135 ) = W( 135 ) + a*JVS( 771 )
  W( 136 ) = W( 136 ) + a*JVS( 772 )
  a = -W( 124 ) / JVS(          913  )
  W( 124 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 914 )
  W( 128 ) = W( 128 ) + a*JVS( 915 )
  W( 129 ) = W( 129 ) + a*JVS( 916 )
  W( 131 ) = W( 131 ) + a*JVS( 917 )
  W( 132 ) = W( 132 ) + a*JVS( 918 )
  W( 133 ) = W( 133 ) + a*JVS( 919 )
  W( 135 ) = W( 135 ) + a*JVS( 920 )
  W( 136 ) = W( 136 ) + a*JVS( 921 )
  JVS( 922) = W( 51 )
  JVS( 923) = W( 115 )
  JVS( 924) = W( 124 )
  JVS( 925) = W( 125 )
  JVS( 926) = W( 127 )
  JVS( 927) = W( 128 )
  JVS( 928) = W( 129 )
  JVS( 929) = W( 131 )
  JVS( 930) = W( 132 )
  JVS( 931) = W( 133 )
  JVS( 932) = W( 135 )
  JVS( 933) = W( 136 )
  IF ( ABS(  JVS( 980 )) < TINY(a) ) THEN
         IER = 126                                     
         RETURN
  END IF
   W( 23 ) = JVS( 934 )
   W( 47 ) = JVS( 935 )
   W( 48 ) = JVS( 936 )
   W( 50 ) = JVS( 937 )
   W( 55 ) = JVS( 938 )
   W( 56 ) = JVS( 939 )
   W( 60 ) = JVS( 940 )
   W( 64 ) = JVS( 941 )
   W( 65 ) = JVS( 942 )
   W( 68 ) = JVS( 943 )
   W( 70 ) = JVS( 944 )
   W( 72 ) = JVS( 945 )
   W( 75 ) = JVS( 946 )
   W( 76 ) = JVS( 947 )
   W( 77 ) = JVS( 948 )
   W( 79 ) = JVS( 949 )
   W( 80 ) = JVS( 950 )
   W( 82 ) = JVS( 951 )
   W( 83 ) = JVS( 952 )
   W( 84 ) = JVS( 953 )
   W( 86 ) = JVS( 954 )
   W( 89 ) = JVS( 955 )
   W( 93 ) = JVS( 956 )
   W( 95 ) = JVS( 957 )
   W( 97 ) = JVS( 958 )
   W( 99 ) = JVS( 959 )
   W( 100 ) = JVS( 960 )
   W( 101 ) = JVS( 961 )
   W( 102 ) = JVS( 962 )
   W( 104 ) = JVS( 963 )
   W( 106 ) = JVS( 964 )
   W( 108 ) = JVS( 965 )
   W( 109 ) = JVS( 966 )
   W( 111 ) = JVS( 967 )
   W( 113 ) = JVS( 968 )
   W( 114 ) = JVS( 969 )
   W( 115 ) = JVS( 970 )
   W( 116 ) = JVS( 971 )
   W( 118 ) = JVS( 972 )
   W( 119 ) = JVS( 973 )
   W( 120 ) = JVS( 974 )
   W( 121 ) = JVS( 975 )
   W( 122 ) = JVS( 976 )
   W( 123 ) = JVS( 977 )
   W( 124 ) = JVS( 978 )
   W( 125 ) = JVS( 979 )
   W( 126 ) = JVS( 980 )
   W( 127 ) = JVS( 981 )
   W( 128 ) = JVS( 982 )
   W( 129 ) = JVS( 983 )
   W( 130 ) = JVS( 984 )
   W( 131 ) = JVS( 985 )
   W( 132 ) = JVS( 986 )
   W( 133 ) = JVS( 987 )
   W( 134 ) = JVS( 988 )
   W( 135 ) = JVS( 989 )
   W( 136 ) = JVS( 990 )
  a = -W( 23 ) / JVS(          121  )
  W( 23 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 122 )
  W( 133 ) = W( 133 ) + a*JVS( 123 )
  a = -W( 47 ) / JVS(          211  )
  W( 47 ) = -a
  W( 121 ) = W( 121 ) + a*JVS( 212 )
  W( 131 ) = W( 131 ) + a*JVS( 213 )
  W( 136 ) = W( 136 ) + a*JVS( 214 )
  a = -W( 48 ) / JVS(          216  )
  W( 48 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 217 )
  W( 136 ) = W( 136 ) + a*JVS( 218 )
  a = -W( 50 ) / JVS(          223  )
  W( 50 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 224 )
  W( 131 ) = W( 131 ) + a*JVS( 225 )
  W( 136 ) = W( 136 ) + a*JVS( 226 )
  a = -W( 55 ) / JVS(          251  )
  W( 55 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 252 )
  W( 136 ) = W( 136 ) + a*JVS( 253 )
  a = -W( 56 ) / JVS(          254  )
  W( 56 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 255 )
  W( 131 ) = W( 131 ) + a*JVS( 256 )
  W( 133 ) = W( 133 ) + a*JVS( 257 )
  a = -W( 60 ) / JVS(          279  )
  W( 60 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 280 )
  W( 133 ) = W( 133 ) + a*JVS( 281 )
  W( 136 ) = W( 136 ) + a*JVS( 282 )
  a = -W( 64 ) / JVS(          304  )
  W( 64 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 305 )
  W( 130 ) = W( 130 ) + a*JVS( 306 )
  W( 131 ) = W( 131 ) + a*JVS( 307 )
  W( 136 ) = W( 136 ) + a*JVS( 308 )
  a = -W( 65 ) / JVS(          309  )
  W( 65 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 310 )
  W( 131 ) = W( 131 ) + a*JVS( 311 )
  W( 136 ) = W( 136 ) + a*JVS( 312 )
  a = -W( 68 ) / JVS(          326  )
  W( 68 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 327 )
  W( 134 ) = W( 134 ) + a*JVS( 328 )
  W( 136 ) = W( 136 ) + a*JVS( 329 )
  a = -W( 70 ) / JVS(          336  )
  W( 70 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 337 )
  W( 131 ) = W( 131 ) + a*JVS( 338 )
  W( 136 ) = W( 136 ) + a*JVS( 339 )
  a = -W( 72 ) / JVS(          349  )
  W( 72 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 350 )
  W( 109 ) = W( 109 ) + a*JVS( 351 )
  W( 132 ) = W( 132 ) + a*JVS( 352 )
  W( 136 ) = W( 136 ) + a*JVS( 353 )
  a = -W( 75 ) / JVS(          367  )
  W( 75 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 368 )
  W( 131 ) = W( 131 ) + a*JVS( 369 )
  W( 136 ) = W( 136 ) + a*JVS( 370 )
  a = -W( 76 ) / JVS(          371  )
  W( 76 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 372 )
  W( 131 ) = W( 131 ) + a*JVS( 373 )
  W( 136 ) = W( 136 ) + a*JVS( 374 )
  a = -W( 77 ) / JVS(          375  )
  W( 77 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 376 )
  W( 134 ) = W( 134 ) + a*JVS( 377 )
  W( 136 ) = W( 136 ) + a*JVS( 378 )
  a = -W( 79 ) / JVS(          394  )
  W( 79 ) = -a
  W( 120 ) = W( 120 ) + a*JVS( 395 )
  W( 131 ) = W( 131 ) + a*JVS( 396 )
  W( 136 ) = W( 136 ) + a*JVS( 397 )
  a = -W( 80 ) / JVS(          401  )
  W( 80 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 402 )
  W( 126 ) = W( 126 ) + a*JVS( 403 )
  W( 132 ) = W( 132 ) + a*JVS( 404 )
  W( 136 ) = W( 136 ) + a*JVS( 405 )
  a = -W( 82 ) / JVS(          415  )
  W( 82 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 416 )
  W( 136 ) = W( 136 ) + a*JVS( 417 )
  a = -W( 83 ) / JVS(          418  )
  W( 83 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 84 ) / JVS(          421  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 422 )
  W( 131 ) = W( 131 ) + a*JVS( 423 )
  W( 136 ) = W( 136 ) + a*JVS( 424 )
  a = -W( 86 ) / JVS(          434  )
  W( 86 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 435 )
  W( 131 ) = W( 131 ) + a*JVS( 436 )
  W( 133 ) = W( 133 ) + a*JVS( 437 )
  W( 136 ) = W( 136 ) + a*JVS( 438 )
  a = -W( 89 ) / JVS(          454  )
  W( 89 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 455 )
  W( 135 ) = W( 135 ) + a*JVS( 456 )
  W( 136 ) = W( 136 ) + a*JVS( 457 )
  a = -W( 93 ) / JVS(          482  )
  W( 93 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 483 )
  W( 131 ) = W( 131 ) + a*JVS( 484 )
  W( 133 ) = W( 133 ) + a*JVS( 485 )
  W( 136 ) = W( 136 ) + a*JVS( 486 )
  a = -W( 95 ) / JVS(          495  )
  W( 95 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 496 )
  W( 111 ) = W( 111 ) + a*JVS( 497 )
  W( 118 ) = W( 118 ) + a*JVS( 498 )
  W( 119 ) = W( 119 ) + a*JVS( 499 )
  W( 120 ) = W( 120 ) + a*JVS( 500 )
  W( 121 ) = W( 121 ) + a*JVS( 501 )
  W( 123 ) = W( 123 ) + a*JVS( 502 )
  W( 125 ) = W( 125 ) + a*JVS( 503 )
  W( 127 ) = W( 127 ) + a*JVS( 504 )
  W( 129 ) = W( 129 ) + a*JVS( 505 )
  W( 136 ) = W( 136 ) + a*JVS( 506 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 533 )
  W( 129 ) = W( 129 ) + a*JVS( 534 )
  W( 131 ) = W( 131 ) + a*JVS( 535 )
  W( 133 ) = W( 133 ) + a*JVS( 536 )
  W( 136 ) = W( 136 ) + a*JVS( 537 )
  a = -W( 99 ) / JVS(          550  )
  W( 99 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 133 ) = W( 133 ) + a*JVS( 553 )
  W( 135 ) = W( 135 ) + a*JVS( 554 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 102 ) / JVS(          564  )
  W( 102 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 565 )
  W( 132 ) = W( 132 ) + a*JVS( 566 )
  W( 133 ) = W( 133 ) + a*JVS( 567 )
  W( 135 ) = W( 135 ) + a*JVS( 568 )
  W( 136 ) = W( 136 ) + a*JVS( 569 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  a = -W( 108 ) / JVS(          668  )
  W( 108 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 669 )
  W( 120 ) = W( 120 ) + a*JVS( 670 )
  W( 123 ) = W( 123 ) + a*JVS( 671 )
  W( 125 ) = W( 125 ) + a*JVS( 672 )
  W( 127 ) = W( 127 ) + a*JVS( 673 )
  W( 128 ) = W( 128 ) + a*JVS( 674 )
  W( 129 ) = W( 129 ) + a*JVS( 675 )
  W( 131 ) = W( 131 ) + a*JVS( 676 )
  W( 132 ) = W( 132 ) + a*JVS( 677 )
  W( 133 ) = W( 133 ) + a*JVS( 678 )
  W( 135 ) = W( 135 ) + a*JVS( 679 )
  W( 136 ) = W( 136 ) + a*JVS( 680 )
  a = -W( 109 ) / JVS(          682  )
  W( 109 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 683 )
  W( 132 ) = W( 132 ) + a*JVS( 684 )
  W( 135 ) = W( 135 ) + a*JVS( 685 )
  W( 136 ) = W( 136 ) + a*JVS( 686 )
  a = -W( 111 ) / JVS(          704  )
  W( 111 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 705 )
  W( 131 ) = W( 131 ) + a*JVS( 706 )
  W( 132 ) = W( 132 ) + a*JVS( 707 )
  W( 133 ) = W( 133 ) + a*JVS( 708 )
  W( 135 ) = W( 135 ) + a*JVS( 709 )
  W( 136 ) = W( 136 ) + a*JVS( 710 )
  a = -W( 113 ) / JVS(          738  )
  W( 113 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 739 )
  W( 132 ) = W( 132 ) + a*JVS( 740 )
  W( 133 ) = W( 133 ) + a*JVS( 741 )
  W( 135 ) = W( 135 ) + a*JVS( 742 )
  W( 136 ) = W( 136 ) + a*JVS( 743 )
  a = -W( 114 ) / JVS(          751  )
  W( 114 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 752 )
  W( 119 ) = W( 119 ) + a*JVS( 753 )
  W( 120 ) = W( 120 ) + a*JVS( 754 )
  W( 129 ) = W( 129 ) + a*JVS( 755 )
  W( 131 ) = W( 131 ) + a*JVS( 756 )
  W( 132 ) = W( 132 ) + a*JVS( 757 )
  W( 133 ) = W( 133 ) + a*JVS( 758 )
  W( 135 ) = W( 135 ) + a*JVS( 759 )
  W( 136 ) = W( 136 ) + a*JVS( 760 )
  a = -W( 115 ) / JVS(          764  )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 765 )
  W( 128 ) = W( 128 ) + a*JVS( 766 )
  W( 129 ) = W( 129 ) + a*JVS( 767 )
  W( 131 ) = W( 131 ) + a*JVS( 768 )
  W( 132 ) = W( 132 ) + a*JVS( 769 )
  W( 133 ) = W( 133 ) + a*JVS( 770 )
  W( 135 ) = W( 135 ) + a*JVS( 771 )
  W( 136 ) = W( 136 ) + a*JVS( 772 )
  a = -W( 116 ) / JVS(          775  )
  W( 116 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 776 )
  W( 119 ) = W( 119 ) + a*JVS( 777 )
  W( 120 ) = W( 120 ) + a*JVS( 778 )
  W( 125 ) = W( 125 ) + a*JVS( 779 )
  W( 127 ) = W( 127 ) + a*JVS( 780 )
  W( 129 ) = W( 129 ) + a*JVS( 781 )
  W( 131 ) = W( 131 ) + a*JVS( 782 )
  W( 132 ) = W( 132 ) + a*JVS( 783 )
  W( 133 ) = W( 133 ) + a*JVS( 784 )
  W( 135 ) = W( 135 ) + a*JVS( 785 )
  W( 136 ) = W( 136 ) + a*JVS( 786 )
  a = -W( 118 ) / JVS(          824  )
  W( 118 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 825 )
  W( 131 ) = W( 131 ) + a*JVS( 826 )
  W( 132 ) = W( 132 ) + a*JVS( 827 )
  W( 133 ) = W( 133 ) + a*JVS( 828 )
  W( 135 ) = W( 135 ) + a*JVS( 829 )
  W( 136 ) = W( 136 ) + a*JVS( 830 )
  a = -W( 119 ) / JVS(          836  )
  W( 119 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 837 )
  W( 131 ) = W( 131 ) + a*JVS( 838 )
  W( 132 ) = W( 132 ) + a*JVS( 839 )
  W( 133 ) = W( 133 ) + a*JVS( 840 )
  W( 135 ) = W( 135 ) + a*JVS( 841 )
  W( 136 ) = W( 136 ) + a*JVS( 842 )
  a = -W( 120 ) / JVS(          847  )
  W( 120 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 848 )
  W( 131 ) = W( 131 ) + a*JVS( 849 )
  W( 132 ) = W( 132 ) + a*JVS( 850 )
  W( 133 ) = W( 133 ) + a*JVS( 851 )
  W( 135 ) = W( 135 ) + a*JVS( 852 )
  W( 136 ) = W( 136 ) + a*JVS( 853 )
  a = -W( 121 ) / JVS(          867  )
  W( 121 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 868 )
  W( 131 ) = W( 131 ) + a*JVS( 869 )
  W( 132 ) = W( 132 ) + a*JVS( 870 )
  W( 133 ) = W( 133 ) + a*JVS( 871 )
  W( 135 ) = W( 135 ) + a*JVS( 872 )
  W( 136 ) = W( 136 ) + a*JVS( 873 )
  a = -W( 122 ) / JVS(          879  )
  W( 122 ) = -a
  W( 123 ) = W( 123 ) + a*JVS( 880 )
  W( 124 ) = W( 124 ) + a*JVS( 881 )
  W( 125 ) = W( 125 ) + a*JVS( 882 )
  W( 127 ) = W( 127 ) + a*JVS( 883 )
  W( 128 ) = W( 128 ) + a*JVS( 884 )
  W( 129 ) = W( 129 ) + a*JVS( 885 )
  W( 130 ) = W( 130 ) + a*JVS( 886 )
  W( 131 ) = W( 131 ) + a*JVS( 887 )
  W( 132 ) = W( 132 ) + a*JVS( 888 )
  W( 133 ) = W( 133 ) + a*JVS( 889 )
  W( 134 ) = W( 134 ) + a*JVS( 890 )
  W( 135 ) = W( 135 ) + a*JVS( 891 )
  W( 136 ) = W( 136 ) + a*JVS( 892 )
  a = -W( 123 ) / JVS(          900  )
  W( 123 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 901 )
  W( 127 ) = W( 127 ) + a*JVS( 902 )
  W( 128 ) = W( 128 ) + a*JVS( 903 )
  W( 129 ) = W( 129 ) + a*JVS( 904 )
  W( 131 ) = W( 131 ) + a*JVS( 905 )
  W( 132 ) = W( 132 ) + a*JVS( 906 )
  W( 133 ) = W( 133 ) + a*JVS( 907 )
  W( 135 ) = W( 135 ) + a*JVS( 908 )
  W( 136 ) = W( 136 ) + a*JVS( 909 )
  a = -W( 124 ) / JVS(          913  )
  W( 124 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 914 )
  W( 128 ) = W( 128 ) + a*JVS( 915 )
  W( 129 ) = W( 129 ) + a*JVS( 916 )
  W( 131 ) = W( 131 ) + a*JVS( 917 )
  W( 132 ) = W( 132 ) + a*JVS( 918 )
  W( 133 ) = W( 133 ) + a*JVS( 919 )
  W( 135 ) = W( 135 ) + a*JVS( 920 )
  W( 136 ) = W( 136 ) + a*JVS( 921 )
  a = -W( 125 ) / JVS(          925  )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 926 )
  W( 128 ) = W( 128 ) + a*JVS( 927 )
  W( 129 ) = W( 129 ) + a*JVS( 928 )
  W( 131 ) = W( 131 ) + a*JVS( 929 )
  W( 132 ) = W( 132 ) + a*JVS( 930 )
  W( 133 ) = W( 133 ) + a*JVS( 931 )
  W( 135 ) = W( 135 ) + a*JVS( 932 )
  W( 136 ) = W( 136 ) + a*JVS( 933 )
  JVS( 934) = W( 23 )
  JVS( 935) = W( 47 )
  JVS( 936) = W( 48 )
  JVS( 937) = W( 50 )
  JVS( 938) = W( 55 )
  JVS( 939) = W( 56 )
  JVS( 940) = W( 60 )
  JVS( 941) = W( 64 )
  JVS( 942) = W( 65 )
  JVS( 943) = W( 68 )
  JVS( 944) = W( 70 )
  JVS( 945) = W( 72 )
  JVS( 946) = W( 75 )
  JVS( 947) = W( 76 )
  JVS( 948) = W( 77 )
  JVS( 949) = W( 79 )
  JVS( 950) = W( 80 )
  JVS( 951) = W( 82 )
  JVS( 952) = W( 83 )
  JVS( 953) = W( 84 )
  JVS( 954) = W( 86 )
  JVS( 955) = W( 89 )
  JVS( 956) = W( 93 )
  JVS( 957) = W( 95 )
  JVS( 958) = W( 97 )
  JVS( 959) = W( 99 )
  JVS( 960) = W( 100 )
  JVS( 961) = W( 101 )
  JVS( 962) = W( 102 )
  JVS( 963) = W( 104 )
  JVS( 964) = W( 106 )
  JVS( 965) = W( 108 )
  JVS( 966) = W( 109 )
  JVS( 967) = W( 111 )
  JVS( 968) = W( 113 )
  JVS( 969) = W( 114 )
  JVS( 970) = W( 115 )
  JVS( 971) = W( 116 )
  JVS( 972) = W( 118 )
  JVS( 973) = W( 119 )
  JVS( 974) = W( 120 )
  JVS( 975) = W( 121 )
  JVS( 976) = W( 122 )
  JVS( 977) = W( 123 )
  JVS( 978) = W( 124 )
  JVS( 979) = W( 125 )
  JVS( 980) = W( 126 )
  JVS( 981) = W( 127 )
  JVS( 982) = W( 128 )
  JVS( 983) = W( 129 )
  JVS( 984) = W( 130 )
  JVS( 985) = W( 131 )
  JVS( 986) = W( 132 )
  JVS( 987) = W( 133 )
  JVS( 988) = W( 134 )
  JVS( 989) = W( 135 )
  JVS( 990) = W( 136 )
  IF ( ABS(  JVS( 993 )) < TINY(a) ) THEN
         IER = 127                                     
         RETURN
  END IF
   W( 75 ) = JVS( 991 )
   W( 106 ) = JVS( 992 )
   W( 127 ) = JVS( 993 )
   W( 128 ) = JVS( 994 )
   W( 129 ) = JVS( 995 )
   W( 131 ) = JVS( 996 )
   W( 132 ) = JVS( 997 )
   W( 133 ) = JVS( 998 )
   W( 135 ) = JVS( 999 )
   W( 136 ) = JVS( 1000 )
  a = -W( 75 ) / JVS(          367  )
  W( 75 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 368 )
  W( 131 ) = W( 131 ) + a*JVS( 369 )
  W( 136 ) = W( 136 ) + a*JVS( 370 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  JVS( 991) = W( 75 )
  JVS( 992) = W( 106 )
  JVS( 993) = W( 127 )
  JVS( 994) = W( 128 )
  JVS( 995) = W( 129 )
  JVS( 996) = W( 131 )
  JVS( 997) = W( 132 )
  JVS( 998) = W( 133 )
  JVS( 999) = W( 135 )
  JVS( 1000) = W( 136 )
  IF ( ABS(  JVS( 1034 )) < TINY(a) ) THEN
         IER = 128                                     
         RETURN
  END IF
   W( 31 ) = JVS( 1001 )
   W( 37 ) = JVS( 1002 )
   W( 40 ) = JVS( 1003 )
   W( 47 ) = JVS( 1004 )
   W( 49 ) = JVS( 1005 )
   W( 62 ) = JVS( 1006 )
   W( 64 ) = JVS( 1007 )
   W( 68 ) = JVS( 1008 )
   W( 82 ) = JVS( 1009 )
   W( 83 ) = JVS( 1010 )
   W( 84 ) = JVS( 1011 )
   W( 85 ) = JVS( 1012 )
   W( 90 ) = JVS( 1013 )
   W( 98 ) = JVS( 1014 )
   W( 100 ) = JVS( 1015 )
   W( 101 ) = JVS( 1016 )
   W( 104 ) = JVS( 1017 )
   W( 105 ) = JVS( 1018 )
   W( 111 ) = JVS( 1019 )
   W( 112 ) = JVS( 1020 )
   W( 113 ) = JVS( 1021 )
   W( 114 ) = JVS( 1022 )
   W( 115 ) = JVS( 1023 )
   W( 117 ) = JVS( 1024 )
   W( 118 ) = JVS( 1025 )
   W( 119 ) = JVS( 1026 )
   W( 120 ) = JVS( 1027 )
   W( 121 ) = JVS( 1028 )
   W( 122 ) = JVS( 1029 )
   W( 123 ) = JVS( 1030 )
   W( 124 ) = JVS( 1031 )
   W( 125 ) = JVS( 1032 )
   W( 127 ) = JVS( 1033 )
   W( 128 ) = JVS( 1034 )
   W( 129 ) = JVS( 1035 )
   W( 130 ) = JVS( 1036 )
   W( 131 ) = JVS( 1037 )
   W( 132 ) = JVS( 1038 )
   W( 133 ) = JVS( 1039 )
   W( 134 ) = JVS( 1040 )
   W( 135 ) = JVS( 1041 )
   W( 136 ) = JVS( 1042 )
  a = -W( 31 ) / JVS(          150  )
  W( 31 ) = -a
  W( 37 ) = W( 37 ) + a*JVS( 151 )
  W( 136 ) = W( 136 ) + a*JVS( 152 )
  a = -W( 37 ) / JVS(          172  )
  W( 37 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 173 )
  a = -W( 40 ) / JVS(          180  )
  W( 40 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 181 )
  W( 133 ) = W( 133 ) + a*JVS( 182 )
  a = -W( 47 ) / JVS(          211  )
  W( 47 ) = -a
  W( 121 ) = W( 121 ) + a*JVS( 212 )
  W( 131 ) = W( 131 ) + a*JVS( 213 )
  W( 136 ) = W( 136 ) + a*JVS( 214 )
  a = -W( 49 ) / JVS(          219  )
  W( 49 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 220 )
  W( 131 ) = W( 131 ) + a*JVS( 221 )
  W( 136 ) = W( 136 ) + a*JVS( 222 )
  a = -W( 62 ) / JVS(          290  )
  W( 62 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 291 )
  W( 113 ) = W( 113 ) + a*JVS( 292 )
  W( 133 ) = W( 133 ) + a*JVS( 293 )
  W( 136 ) = W( 136 ) + a*JVS( 294 )
  a = -W( 64 ) / JVS(          304  )
  W( 64 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 305 )
  W( 130 ) = W( 130 ) + a*JVS( 306 )
  W( 131 ) = W( 131 ) + a*JVS( 307 )
  W( 136 ) = W( 136 ) + a*JVS( 308 )
  a = -W( 68 ) / JVS(          326  )
  W( 68 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 327 )
  W( 134 ) = W( 134 ) + a*JVS( 328 )
  W( 136 ) = W( 136 ) + a*JVS( 329 )
  a = -W( 82 ) / JVS(          415  )
  W( 82 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 416 )
  W( 136 ) = W( 136 ) + a*JVS( 417 )
  a = -W( 83 ) / JVS(          418  )
  W( 83 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 84 ) / JVS(          421  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 422 )
  W( 131 ) = W( 131 ) + a*JVS( 423 )
  W( 136 ) = W( 136 ) + a*JVS( 424 )
  a = -W( 85 ) / JVS(          428  )
  W( 85 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 429 )
  W( 131 ) = W( 131 ) + a*JVS( 430 )
  W( 133 ) = W( 133 ) + a*JVS( 431 )
  W( 136 ) = W( 136 ) + a*JVS( 432 )
  a = -W( 90 ) / JVS(          460  )
  W( 90 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 461 )
  W( 133 ) = W( 133 ) + a*JVS( 462 )
  W( 136 ) = W( 136 ) + a*JVS( 463 )
  a = -W( 98 ) / JVS(          542  )
  W( 98 ) = -a
  W( 100 ) = W( 100 ) + a*JVS( 543 )
  W( 101 ) = W( 101 ) + a*JVS( 544 )
  W( 104 ) = W( 104 ) + a*JVS( 545 )
  W( 131 ) = W( 131 ) + a*JVS( 546 )
  W( 132 ) = W( 132 ) + a*JVS( 547 )
  W( 133 ) = W( 133 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  a = -W( 105 ) / JVS(          601  )
  W( 105 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 602 )
  W( 113 ) = W( 113 ) + a*JVS( 603 )
  W( 114 ) = W( 114 ) + a*JVS( 604 )
  W( 118 ) = W( 118 ) + a*JVS( 605 )
  W( 120 ) = W( 120 ) + a*JVS( 606 )
  W( 129 ) = W( 129 ) + a*JVS( 607 )
  W( 131 ) = W( 131 ) + a*JVS( 608 )
  W( 132 ) = W( 132 ) + a*JVS( 609 )
  W( 133 ) = W( 133 ) + a*JVS( 610 )
  W( 135 ) = W( 135 ) + a*JVS( 611 )
  W( 136 ) = W( 136 ) + a*JVS( 612 )
  a = -W( 111 ) / JVS(          704  )
  W( 111 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 705 )
  W( 131 ) = W( 131 ) + a*JVS( 706 )
  W( 132 ) = W( 132 ) + a*JVS( 707 )
  W( 133 ) = W( 133 ) + a*JVS( 708 )
  W( 135 ) = W( 135 ) + a*JVS( 709 )
  W( 136 ) = W( 136 ) + a*JVS( 710 )
  a = -W( 112 ) / JVS(          722  )
  W( 112 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 723 )
  W( 124 ) = W( 124 ) + a*JVS( 724 )
  W( 129 ) = W( 129 ) + a*JVS( 725 )
  W( 131 ) = W( 131 ) + a*JVS( 726 )
  W( 132 ) = W( 132 ) + a*JVS( 727 )
  W( 133 ) = W( 133 ) + a*JVS( 728 )
  W( 135 ) = W( 135 ) + a*JVS( 729 )
  W( 136 ) = W( 136 ) + a*JVS( 730 )
  a = -W( 113 ) / JVS(          738  )
  W( 113 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 739 )
  W( 132 ) = W( 132 ) + a*JVS( 740 )
  W( 133 ) = W( 133 ) + a*JVS( 741 )
  W( 135 ) = W( 135 ) + a*JVS( 742 )
  W( 136 ) = W( 136 ) + a*JVS( 743 )
  a = -W( 114 ) / JVS(          751  )
  W( 114 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 752 )
  W( 119 ) = W( 119 ) + a*JVS( 753 )
  W( 120 ) = W( 120 ) + a*JVS( 754 )
  W( 129 ) = W( 129 ) + a*JVS( 755 )
  W( 131 ) = W( 131 ) + a*JVS( 756 )
  W( 132 ) = W( 132 ) + a*JVS( 757 )
  W( 133 ) = W( 133 ) + a*JVS( 758 )
  W( 135 ) = W( 135 ) + a*JVS( 759 )
  W( 136 ) = W( 136 ) + a*JVS( 760 )
  a = -W( 115 ) / JVS(          764  )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 765 )
  W( 128 ) = W( 128 ) + a*JVS( 766 )
  W( 129 ) = W( 129 ) + a*JVS( 767 )
  W( 131 ) = W( 131 ) + a*JVS( 768 )
  W( 132 ) = W( 132 ) + a*JVS( 769 )
  W( 133 ) = W( 133 ) + a*JVS( 770 )
  W( 135 ) = W( 135 ) + a*JVS( 771 )
  W( 136 ) = W( 136 ) + a*JVS( 772 )
  a = -W( 117 ) / JVS(          803  )
  W( 117 ) = -a
  W( 121 ) = W( 121 ) + a*JVS( 804 )
  W( 122 ) = W( 122 ) + a*JVS( 805 )
  W( 123 ) = W( 123 ) + a*JVS( 806 )
  W( 124 ) = W( 124 ) + a*JVS( 807 )
  W( 125 ) = W( 125 ) + a*JVS( 808 )
  W( 127 ) = W( 127 ) + a*JVS( 809 )
  W( 128 ) = W( 128 ) + a*JVS( 810 )
  W( 129 ) = W( 129 ) + a*JVS( 811 )
  W( 131 ) = W( 131 ) + a*JVS( 812 )
  W( 132 ) = W( 132 ) + a*JVS( 813 )
  W( 133 ) = W( 133 ) + a*JVS( 814 )
  W( 134 ) = W( 134 ) + a*JVS( 815 )
  W( 135 ) = W( 135 ) + a*JVS( 816 )
  W( 136 ) = W( 136 ) + a*JVS( 817 )
  a = -W( 118 ) / JVS(          824  )
  W( 118 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 825 )
  W( 131 ) = W( 131 ) + a*JVS( 826 )
  W( 132 ) = W( 132 ) + a*JVS( 827 )
  W( 133 ) = W( 133 ) + a*JVS( 828 )
  W( 135 ) = W( 135 ) + a*JVS( 829 )
  W( 136 ) = W( 136 ) + a*JVS( 830 )
  a = -W( 119 ) / JVS(          836  )
  W( 119 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 837 )
  W( 131 ) = W( 131 ) + a*JVS( 838 )
  W( 132 ) = W( 132 ) + a*JVS( 839 )
  W( 133 ) = W( 133 ) + a*JVS( 840 )
  W( 135 ) = W( 135 ) + a*JVS( 841 )
  W( 136 ) = W( 136 ) + a*JVS( 842 )
  a = -W( 120 ) / JVS(          847  )
  W( 120 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 848 )
  W( 131 ) = W( 131 ) + a*JVS( 849 )
  W( 132 ) = W( 132 ) + a*JVS( 850 )
  W( 133 ) = W( 133 ) + a*JVS( 851 )
  W( 135 ) = W( 135 ) + a*JVS( 852 )
  W( 136 ) = W( 136 ) + a*JVS( 853 )
  a = -W( 121 ) / JVS(          867  )
  W( 121 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 868 )
  W( 131 ) = W( 131 ) + a*JVS( 869 )
  W( 132 ) = W( 132 ) + a*JVS( 870 )
  W( 133 ) = W( 133 ) + a*JVS( 871 )
  W( 135 ) = W( 135 ) + a*JVS( 872 )
  W( 136 ) = W( 136 ) + a*JVS( 873 )
  a = -W( 122 ) / JVS(          879  )
  W( 122 ) = -a
  W( 123 ) = W( 123 ) + a*JVS( 880 )
  W( 124 ) = W( 124 ) + a*JVS( 881 )
  W( 125 ) = W( 125 ) + a*JVS( 882 )
  W( 127 ) = W( 127 ) + a*JVS( 883 )
  W( 128 ) = W( 128 ) + a*JVS( 884 )
  W( 129 ) = W( 129 ) + a*JVS( 885 )
  W( 130 ) = W( 130 ) + a*JVS( 886 )
  W( 131 ) = W( 131 ) + a*JVS( 887 )
  W( 132 ) = W( 132 ) + a*JVS( 888 )
  W( 133 ) = W( 133 ) + a*JVS( 889 )
  W( 134 ) = W( 134 ) + a*JVS( 890 )
  W( 135 ) = W( 135 ) + a*JVS( 891 )
  W( 136 ) = W( 136 ) + a*JVS( 892 )
  a = -W( 123 ) / JVS(          900  )
  W( 123 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 901 )
  W( 127 ) = W( 127 ) + a*JVS( 902 )
  W( 128 ) = W( 128 ) + a*JVS( 903 )
  W( 129 ) = W( 129 ) + a*JVS( 904 )
  W( 131 ) = W( 131 ) + a*JVS( 905 )
  W( 132 ) = W( 132 ) + a*JVS( 906 )
  W( 133 ) = W( 133 ) + a*JVS( 907 )
  W( 135 ) = W( 135 ) + a*JVS( 908 )
  W( 136 ) = W( 136 ) + a*JVS( 909 )
  a = -W( 124 ) / JVS(          913  )
  W( 124 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 914 )
  W( 128 ) = W( 128 ) + a*JVS( 915 )
  W( 129 ) = W( 129 ) + a*JVS( 916 )
  W( 131 ) = W( 131 ) + a*JVS( 917 )
  W( 132 ) = W( 132 ) + a*JVS( 918 )
  W( 133 ) = W( 133 ) + a*JVS( 919 )
  W( 135 ) = W( 135 ) + a*JVS( 920 )
  W( 136 ) = W( 136 ) + a*JVS( 921 )
  a = -W( 125 ) / JVS(          925  )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 926 )
  W( 128 ) = W( 128 ) + a*JVS( 927 )
  W( 129 ) = W( 129 ) + a*JVS( 928 )
  W( 131 ) = W( 131 ) + a*JVS( 929 )
  W( 132 ) = W( 132 ) + a*JVS( 930 )
  W( 133 ) = W( 133 ) + a*JVS( 931 )
  W( 135 ) = W( 135 ) + a*JVS( 932 )
  W( 136 ) = W( 136 ) + a*JVS( 933 )
  a = -W( 127 ) / JVS(          993  )
  W( 127 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 994 )
  W( 129 ) = W( 129 ) + a*JVS( 995 )
  W( 131 ) = W( 131 ) + a*JVS( 996 )
  W( 132 ) = W( 132 ) + a*JVS( 997 )
  W( 133 ) = W( 133 ) + a*JVS( 998 )
  W( 135 ) = W( 135 ) + a*JVS( 999 )
  W( 136 ) = W( 136 ) + a*JVS( 1000 )
  JVS( 1001) = W( 31 )
  JVS( 1002) = W( 37 )
  JVS( 1003) = W( 40 )
  JVS( 1004) = W( 47 )
  JVS( 1005) = W( 49 )
  JVS( 1006) = W( 62 )
  JVS( 1007) = W( 64 )
  JVS( 1008) = W( 68 )
  JVS( 1009) = W( 82 )
  JVS( 1010) = W( 83 )
  JVS( 1011) = W( 84 )
  JVS( 1012) = W( 85 )
  JVS( 1013) = W( 90 )
  JVS( 1014) = W( 98 )
  JVS( 1015) = W( 100 )
  JVS( 1016) = W( 101 )
  JVS( 1017) = W( 104 )
  JVS( 1018) = W( 105 )
  JVS( 1019) = W( 111 )
  JVS( 1020) = W( 112 )
  JVS( 1021) = W( 113 )
  JVS( 1022) = W( 114 )
  JVS( 1023) = W( 115 )
  JVS( 1024) = W( 117 )
  JVS( 1025) = W( 118 )
  JVS( 1026) = W( 119 )
  JVS( 1027) = W( 120 )
  JVS( 1028) = W( 121 )
  JVS( 1029) = W( 122 )
  JVS( 1030) = W( 123 )
  JVS( 1031) = W( 124 )
  JVS( 1032) = W( 125 )
  JVS( 1033) = W( 127 )
  JVS( 1034) = W( 128 )
  JVS( 1035) = W( 129 )
  JVS( 1036) = W( 130 )
  JVS( 1037) = W( 131 )
  JVS( 1038) = W( 132 )
  JVS( 1039) = W( 133 )
  JVS( 1040) = W( 134 )
  JVS( 1041) = W( 135 )
  JVS( 1042) = W( 136 )
  IF ( ABS(  JVS( 1070 )) < TINY(a) ) THEN
         IER = 129                                     
         RETURN
  END IF
   W( 50 ) = JVS( 1043 )
   W( 63 ) = JVS( 1044 )
   W( 64 ) = JVS( 1045 )
   W( 68 ) = JVS( 1046 )
   W( 72 ) = JVS( 1047 )
   W( 80 ) = JVS( 1048 )
   W( 91 ) = JVS( 1049 )
   W( 93 ) = JVS( 1050 )
   W( 94 ) = JVS( 1051 )
   W( 97 ) = JVS( 1052 )
   W( 105 ) = JVS( 1053 )
   W( 106 ) = JVS( 1054 )
   W( 109 ) = JVS( 1055 )
   W( 111 ) = JVS( 1056 )
   W( 112 ) = JVS( 1057 )
   W( 113 ) = JVS( 1058 )
   W( 114 ) = JVS( 1059 )
   W( 118 ) = JVS( 1060 )
   W( 119 ) = JVS( 1061 )
   W( 120 ) = JVS( 1062 )
   W( 121 ) = JVS( 1063 )
   W( 123 ) = JVS( 1064 )
   W( 124 ) = JVS( 1065 )
   W( 125 ) = JVS( 1066 )
   W( 126 ) = JVS( 1067 )
   W( 127 ) = JVS( 1068 )
   W( 128 ) = JVS( 1069 )
   W( 129 ) = JVS( 1070 )
   W( 130 ) = JVS( 1071 )
   W( 131 ) = JVS( 1072 )
   W( 132 ) = JVS( 1073 )
   W( 133 ) = JVS( 1074 )
   W( 134 ) = JVS( 1075 )
   W( 135 ) = JVS( 1076 )
   W( 136 ) = JVS( 1077 )
  a = -W( 50 ) / JVS(          223  )
  W( 50 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 224 )
  W( 131 ) = W( 131 ) + a*JVS( 225 )
  W( 136 ) = W( 136 ) + a*JVS( 226 )
  a = -W( 63 ) / JVS(          295  )
  W( 63 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 296 )
  W( 109 ) = W( 109 ) + a*JVS( 297 )
  W( 128 ) = W( 128 ) + a*JVS( 298 )
  W( 129 ) = W( 129 ) + a*JVS( 299 )
  W( 130 ) = W( 130 ) + a*JVS( 300 )
  W( 131 ) = W( 131 ) + a*JVS( 301 )
  W( 132 ) = W( 132 ) + a*JVS( 302 )
  W( 136 ) = W( 136 ) + a*JVS( 303 )
  a = -W( 64 ) / JVS(          304  )
  W( 64 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 305 )
  W( 130 ) = W( 130 ) + a*JVS( 306 )
  W( 131 ) = W( 131 ) + a*JVS( 307 )
  W( 136 ) = W( 136 ) + a*JVS( 308 )
  a = -W( 68 ) / JVS(          326  )
  W( 68 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 327 )
  W( 134 ) = W( 134 ) + a*JVS( 328 )
  W( 136 ) = W( 136 ) + a*JVS( 329 )
  a = -W( 72 ) / JVS(          349  )
  W( 72 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 350 )
  W( 109 ) = W( 109 ) + a*JVS( 351 )
  W( 132 ) = W( 132 ) + a*JVS( 352 )
  W( 136 ) = W( 136 ) + a*JVS( 353 )
  a = -W( 80 ) / JVS(          401  )
  W( 80 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 402 )
  W( 126 ) = W( 126 ) + a*JVS( 403 )
  W( 132 ) = W( 132 ) + a*JVS( 404 )
  W( 136 ) = W( 136 ) + a*JVS( 405 )
  a = -W( 91 ) / JVS(          467  )
  W( 91 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 468 )
  W( 133 ) = W( 133 ) + a*JVS( 469 )
  W( 134 ) = W( 134 ) + a*JVS( 470 )
  W( 136 ) = W( 136 ) + a*JVS( 471 )
  a = -W( 93 ) / JVS(          482  )
  W( 93 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 483 )
  W( 131 ) = W( 131 ) + a*JVS( 484 )
  W( 133 ) = W( 133 ) + a*JVS( 485 )
  W( 136 ) = W( 136 ) + a*JVS( 486 )
  a = -W( 94 ) / JVS(          490  )
  W( 94 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 491 )
  W( 133 ) = W( 133 ) + a*JVS( 492 )
  W( 134 ) = W( 134 ) + a*JVS( 493 )
  W( 136 ) = W( 136 ) + a*JVS( 494 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 533 )
  W( 129 ) = W( 129 ) + a*JVS( 534 )
  W( 131 ) = W( 131 ) + a*JVS( 535 )
  W( 133 ) = W( 133 ) + a*JVS( 536 )
  W( 136 ) = W( 136 ) + a*JVS( 537 )
  a = -W( 105 ) / JVS(          601  )
  W( 105 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 602 )
  W( 113 ) = W( 113 ) + a*JVS( 603 )
  W( 114 ) = W( 114 ) + a*JVS( 604 )
  W( 118 ) = W( 118 ) + a*JVS( 605 )
  W( 120 ) = W( 120 ) + a*JVS( 606 )
  W( 129 ) = W( 129 ) + a*JVS( 607 )
  W( 131 ) = W( 131 ) + a*JVS( 608 )
  W( 132 ) = W( 132 ) + a*JVS( 609 )
  W( 133 ) = W( 133 ) + a*JVS( 610 )
  W( 135 ) = W( 135 ) + a*JVS( 611 )
  W( 136 ) = W( 136 ) + a*JVS( 612 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  a = -W( 109 ) / JVS(          682  )
  W( 109 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 683 )
  W( 132 ) = W( 132 ) + a*JVS( 684 )
  W( 135 ) = W( 135 ) + a*JVS( 685 )
  W( 136 ) = W( 136 ) + a*JVS( 686 )
  a = -W( 111 ) / JVS(          704  )
  W( 111 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 705 )
  W( 131 ) = W( 131 ) + a*JVS( 706 )
  W( 132 ) = W( 132 ) + a*JVS( 707 )
  W( 133 ) = W( 133 ) + a*JVS( 708 )
  W( 135 ) = W( 135 ) + a*JVS( 709 )
  W( 136 ) = W( 136 ) + a*JVS( 710 )
  a = -W( 112 ) / JVS(          722  )
  W( 112 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 723 )
  W( 124 ) = W( 124 ) + a*JVS( 724 )
  W( 129 ) = W( 129 ) + a*JVS( 725 )
  W( 131 ) = W( 131 ) + a*JVS( 726 )
  W( 132 ) = W( 132 ) + a*JVS( 727 )
  W( 133 ) = W( 133 ) + a*JVS( 728 )
  W( 135 ) = W( 135 ) + a*JVS( 729 )
  W( 136 ) = W( 136 ) + a*JVS( 730 )
  a = -W( 113 ) / JVS(          738  )
  W( 113 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 739 )
  W( 132 ) = W( 132 ) + a*JVS( 740 )
  W( 133 ) = W( 133 ) + a*JVS( 741 )
  W( 135 ) = W( 135 ) + a*JVS( 742 )
  W( 136 ) = W( 136 ) + a*JVS( 743 )
  a = -W( 114 ) / JVS(          751  )
  W( 114 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 752 )
  W( 119 ) = W( 119 ) + a*JVS( 753 )
  W( 120 ) = W( 120 ) + a*JVS( 754 )
  W( 129 ) = W( 129 ) + a*JVS( 755 )
  W( 131 ) = W( 131 ) + a*JVS( 756 )
  W( 132 ) = W( 132 ) + a*JVS( 757 )
  W( 133 ) = W( 133 ) + a*JVS( 758 )
  W( 135 ) = W( 135 ) + a*JVS( 759 )
  W( 136 ) = W( 136 ) + a*JVS( 760 )
  a = -W( 118 ) / JVS(          824  )
  W( 118 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 825 )
  W( 131 ) = W( 131 ) + a*JVS( 826 )
  W( 132 ) = W( 132 ) + a*JVS( 827 )
  W( 133 ) = W( 133 ) + a*JVS( 828 )
  W( 135 ) = W( 135 ) + a*JVS( 829 )
  W( 136 ) = W( 136 ) + a*JVS( 830 )
  a = -W( 119 ) / JVS(          836  )
  W( 119 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 837 )
  W( 131 ) = W( 131 ) + a*JVS( 838 )
  W( 132 ) = W( 132 ) + a*JVS( 839 )
  W( 133 ) = W( 133 ) + a*JVS( 840 )
  W( 135 ) = W( 135 ) + a*JVS( 841 )
  W( 136 ) = W( 136 ) + a*JVS( 842 )
  a = -W( 120 ) / JVS(          847  )
  W( 120 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 848 )
  W( 131 ) = W( 131 ) + a*JVS( 849 )
  W( 132 ) = W( 132 ) + a*JVS( 850 )
  W( 133 ) = W( 133 ) + a*JVS( 851 )
  W( 135 ) = W( 135 ) + a*JVS( 852 )
  W( 136 ) = W( 136 ) + a*JVS( 853 )
  a = -W( 121 ) / JVS(          867  )
  W( 121 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 868 )
  W( 131 ) = W( 131 ) + a*JVS( 869 )
  W( 132 ) = W( 132 ) + a*JVS( 870 )
  W( 133 ) = W( 133 ) + a*JVS( 871 )
  W( 135 ) = W( 135 ) + a*JVS( 872 )
  W( 136 ) = W( 136 ) + a*JVS( 873 )
  a = -W( 123 ) / JVS(          900  )
  W( 123 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 901 )
  W( 127 ) = W( 127 ) + a*JVS( 902 )
  W( 128 ) = W( 128 ) + a*JVS( 903 )
  W( 129 ) = W( 129 ) + a*JVS( 904 )
  W( 131 ) = W( 131 ) + a*JVS( 905 )
  W( 132 ) = W( 132 ) + a*JVS( 906 )
  W( 133 ) = W( 133 ) + a*JVS( 907 )
  W( 135 ) = W( 135 ) + a*JVS( 908 )
  W( 136 ) = W( 136 ) + a*JVS( 909 )
  a = -W( 124 ) / JVS(          913  )
  W( 124 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 914 )
  W( 128 ) = W( 128 ) + a*JVS( 915 )
  W( 129 ) = W( 129 ) + a*JVS( 916 )
  W( 131 ) = W( 131 ) + a*JVS( 917 )
  W( 132 ) = W( 132 ) + a*JVS( 918 )
  W( 133 ) = W( 133 ) + a*JVS( 919 )
  W( 135 ) = W( 135 ) + a*JVS( 920 )
  W( 136 ) = W( 136 ) + a*JVS( 921 )
  a = -W( 125 ) / JVS(          925  )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 926 )
  W( 128 ) = W( 128 ) + a*JVS( 927 )
  W( 129 ) = W( 129 ) + a*JVS( 928 )
  W( 131 ) = W( 131 ) + a*JVS( 929 )
  W( 132 ) = W( 132 ) + a*JVS( 930 )
  W( 133 ) = W( 133 ) + a*JVS( 931 )
  W( 135 ) = W( 135 ) + a*JVS( 932 )
  W( 136 ) = W( 136 ) + a*JVS( 933 )
  a = -W( 126 ) / JVS(          980  )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 981 )
  W( 128 ) = W( 128 ) + a*JVS( 982 )
  W( 129 ) = W( 129 ) + a*JVS( 983 )
  W( 130 ) = W( 130 ) + a*JVS( 984 )
  W( 131 ) = W( 131 ) + a*JVS( 985 )
  W( 132 ) = W( 132 ) + a*JVS( 986 )
  W( 133 ) = W( 133 ) + a*JVS( 987 )
  W( 134 ) = W( 134 ) + a*JVS( 988 )
  W( 135 ) = W( 135 ) + a*JVS( 989 )
  W( 136 ) = W( 136 ) + a*JVS( 990 )
  a = -W( 127 ) / JVS(          993  )
  W( 127 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 994 )
  W( 129 ) = W( 129 ) + a*JVS( 995 )
  W( 131 ) = W( 131 ) + a*JVS( 996 )
  W( 132 ) = W( 132 ) + a*JVS( 997 )
  W( 133 ) = W( 133 ) + a*JVS( 998 )
  W( 135 ) = W( 135 ) + a*JVS( 999 )
  W( 136 ) = W( 136 ) + a*JVS( 1000 )
  a = -W( 128 ) / JVS(         1034  )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 1035 )
  W( 130 ) = W( 130 ) + a*JVS( 1036 )
  W( 131 ) = W( 131 ) + a*JVS( 1037 )
  W( 132 ) = W( 132 ) + a*JVS( 1038 )
  W( 133 ) = W( 133 ) + a*JVS( 1039 )
  W( 134 ) = W( 134 ) + a*JVS( 1040 )
  W( 135 ) = W( 135 ) + a*JVS( 1041 )
  W( 136 ) = W( 136 ) + a*JVS( 1042 )
  JVS( 1043) = W( 50 )
  JVS( 1044) = W( 63 )
  JVS( 1045) = W( 64 )
  JVS( 1046) = W( 68 )
  JVS( 1047) = W( 72 )
  JVS( 1048) = W( 80 )
  JVS( 1049) = W( 91 )
  JVS( 1050) = W( 93 )
  JVS( 1051) = W( 94 )
  JVS( 1052) = W( 97 )
  JVS( 1053) = W( 105 )
  JVS( 1054) = W( 106 )
  JVS( 1055) = W( 109 )
  JVS( 1056) = W( 111 )
  JVS( 1057) = W( 112 )
  JVS( 1058) = W( 113 )
  JVS( 1059) = W( 114 )
  JVS( 1060) = W( 118 )
  JVS( 1061) = W( 119 )
  JVS( 1062) = W( 120 )
  JVS( 1063) = W( 121 )
  JVS( 1064) = W( 123 )
  JVS( 1065) = W( 124 )
  JVS( 1066) = W( 125 )
  JVS( 1067) = W( 126 )
  JVS( 1068) = W( 127 )
  JVS( 1069) = W( 128 )
  JVS( 1070) = W( 129 )
  JVS( 1071) = W( 130 )
  JVS( 1072) = W( 131 )
  JVS( 1073) = W( 132 )
  JVS( 1074) = W( 133 )
  JVS( 1075) = W( 134 )
  JVS( 1076) = W( 135 )
  JVS( 1077) = W( 136 )
  IF ( ABS(  JVS( 1086 )) < TINY(a) ) THEN
         IER = 130                                     
         RETURN
  END IF
   W( 51 ) = JVS( 1078 )
   W( 77 ) = JVS( 1079 )
   W( 106 ) = JVS( 1080 )
   W( 115 ) = JVS( 1081 )
   W( 125 ) = JVS( 1082 )
   W( 127 ) = JVS( 1083 )
   W( 128 ) = JVS( 1084 )
   W( 129 ) = JVS( 1085 )
   W( 130 ) = JVS( 1086 )
   W( 131 ) = JVS( 1087 )
   W( 132 ) = JVS( 1088 )
   W( 133 ) = JVS( 1089 )
   W( 134 ) = JVS( 1090 )
   W( 135 ) = JVS( 1091 )
   W( 136 ) = JVS( 1092 )
  a = -W( 51 ) / JVS(          227  )
  W( 51 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 228 )
  W( 131 ) = W( 131 ) + a*JVS( 229 )
  W( 136 ) = W( 136 ) + a*JVS( 230 )
  a = -W( 77 ) / JVS(          375  )
  W( 77 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 376 )
  W( 134 ) = W( 134 ) + a*JVS( 377 )
  W( 136 ) = W( 136 ) + a*JVS( 378 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  a = -W( 115 ) / JVS(          764  )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 765 )
  W( 128 ) = W( 128 ) + a*JVS( 766 )
  W( 129 ) = W( 129 ) + a*JVS( 767 )
  W( 131 ) = W( 131 ) + a*JVS( 768 )
  W( 132 ) = W( 132 ) + a*JVS( 769 )
  W( 133 ) = W( 133 ) + a*JVS( 770 )
  W( 135 ) = W( 135 ) + a*JVS( 771 )
  W( 136 ) = W( 136 ) + a*JVS( 772 )
  a = -W( 125 ) / JVS(          925  )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 926 )
  W( 128 ) = W( 128 ) + a*JVS( 927 )
  W( 129 ) = W( 129 ) + a*JVS( 928 )
  W( 131 ) = W( 131 ) + a*JVS( 929 )
  W( 132 ) = W( 132 ) + a*JVS( 930 )
  W( 133 ) = W( 133 ) + a*JVS( 931 )
  W( 135 ) = W( 135 ) + a*JVS( 932 )
  W( 136 ) = W( 136 ) + a*JVS( 933 )
  a = -W( 127 ) / JVS(          993  )
  W( 127 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 994 )
  W( 129 ) = W( 129 ) + a*JVS( 995 )
  W( 131 ) = W( 131 ) + a*JVS( 996 )
  W( 132 ) = W( 132 ) + a*JVS( 997 )
  W( 133 ) = W( 133 ) + a*JVS( 998 )
  W( 135 ) = W( 135 ) + a*JVS( 999 )
  W( 136 ) = W( 136 ) + a*JVS( 1000 )
  a = -W( 128 ) / JVS(         1034  )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 1035 )
  W( 130 ) = W( 130 ) + a*JVS( 1036 )
  W( 131 ) = W( 131 ) + a*JVS( 1037 )
  W( 132 ) = W( 132 ) + a*JVS( 1038 )
  W( 133 ) = W( 133 ) + a*JVS( 1039 )
  W( 134 ) = W( 134 ) + a*JVS( 1040 )
  W( 135 ) = W( 135 ) + a*JVS( 1041 )
  W( 136 ) = W( 136 ) + a*JVS( 1042 )
  a = -W( 129 ) / JVS(         1070  )
  W( 129 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 1071 )
  W( 131 ) = W( 131 ) + a*JVS( 1072 )
  W( 132 ) = W( 132 ) + a*JVS( 1073 )
  W( 133 ) = W( 133 ) + a*JVS( 1074 )
  W( 134 ) = W( 134 ) + a*JVS( 1075 )
  W( 135 ) = W( 135 ) + a*JVS( 1076 )
  W( 136 ) = W( 136 ) + a*JVS( 1077 )
  JVS( 1078) = W( 51 )
  JVS( 1079) = W( 77 )
  JVS( 1080) = W( 106 )
  JVS( 1081) = W( 115 )
  JVS( 1082) = W( 125 )
  JVS( 1083) = W( 127 )
  JVS( 1084) = W( 128 )
  JVS( 1085) = W( 129 )
  JVS( 1086) = W( 130 )
  JVS( 1087) = W( 131 )
  JVS( 1088) = W( 132 )
  JVS( 1089) = W( 133 )
  JVS( 1090) = W( 134 )
  JVS( 1091) = W( 135 )
  JVS( 1092) = W( 136 )
  IF ( ABS(  JVS( 1184 )) < TINY(a) ) THEN
         IER = 131                                     
         RETURN
  END IF
   W( 18 ) = JVS( 1093 )
   W( 19 ) = JVS( 1094 )
   W( 20 ) = JVS( 1095 )
   W( 21 ) = JVS( 1096 )
   W( 22 ) = JVS( 1097 )
   W( 23 ) = JVS( 1098 )
   W( 24 ) = JVS( 1099 )
   W( 25 ) = JVS( 1100 )
   W( 30 ) = JVS( 1101 )
   W( 31 ) = JVS( 1102 )
   W( 37 ) = JVS( 1103 )
   W( 38 ) = JVS( 1104 )
   W( 39 ) = JVS( 1105 )
   W( 40 ) = JVS( 1106 )
   W( 41 ) = JVS( 1107 )
   W( 42 ) = JVS( 1108 )
   W( 43 ) = JVS( 1109 )
   W( 44 ) = JVS( 1110 )
   W( 45 ) = JVS( 1111 )
   W( 46 ) = JVS( 1112 )
   W( 48 ) = JVS( 1113 )
   W( 50 ) = JVS( 1114 )
   W( 51 ) = JVS( 1115 )
   W( 52 ) = JVS( 1116 )
   W( 55 ) = JVS( 1117 )
   W( 56 ) = JVS( 1118 )
   W( 58 ) = JVS( 1119 )
   W( 60 ) = JVS( 1120 )
   W( 61 ) = JVS( 1121 )
   W( 65 ) = JVS( 1122 )
   W( 66 ) = JVS( 1123 )
   W( 69 ) = JVS( 1124 )
   W( 70 ) = JVS( 1125 )
   W( 72 ) = JVS( 1126 )
   W( 73 ) = JVS( 1127 )
   W( 74 ) = JVS( 1128 )
   W( 75 ) = JVS( 1129 )
   W( 76 ) = JVS( 1130 )
   W( 77 ) = JVS( 1131 )
   W( 78 ) = JVS( 1132 )
   W( 79 ) = JVS( 1133 )
   W( 80 ) = JVS( 1134 )
   W( 81 ) = JVS( 1135 )
   W( 82 ) = JVS( 1136 )
   W( 83 ) = JVS( 1137 )
   W( 84 ) = JVS( 1138 )
   W( 85 ) = JVS( 1139 )
   W( 86 ) = JVS( 1140 )
   W( 87 ) = JVS( 1141 )
   W( 88 ) = JVS( 1142 )
   W( 89 ) = JVS( 1143 )
   W( 90 ) = JVS( 1144 )
   W( 91 ) = JVS( 1145 )
   W( 92 ) = JVS( 1146 )
   W( 93 ) = JVS( 1147 )
   W( 94 ) = JVS( 1148 )
   W( 95 ) = JVS( 1149 )
   W( 96 ) = JVS( 1150 )
   W( 97 ) = JVS( 1151 )
   W( 98 ) = JVS( 1152 )
   W( 99 ) = JVS( 1153 )
   W( 100 ) = JVS( 1154 )
   W( 101 ) = JVS( 1155 )
   W( 102 ) = JVS( 1156 )
   W( 103 ) = JVS( 1157 )
   W( 104 ) = JVS( 1158 )
   W( 106 ) = JVS( 1159 )
   W( 107 ) = JVS( 1160 )
   W( 108 ) = JVS( 1161 )
   W( 109 ) = JVS( 1162 )
   W( 110 ) = JVS( 1163 )
   W( 111 ) = JVS( 1164 )
   W( 112 ) = JVS( 1165 )
   W( 113 ) = JVS( 1166 )
   W( 114 ) = JVS( 1167 )
   W( 115 ) = JVS( 1168 )
   W( 116 ) = JVS( 1169 )
   W( 117 ) = JVS( 1170 )
   W( 118 ) = JVS( 1171 )
   W( 119 ) = JVS( 1172 )
   W( 120 ) = JVS( 1173 )
   W( 121 ) = JVS( 1174 )
   W( 122 ) = JVS( 1175 )
   W( 123 ) = JVS( 1176 )
   W( 124 ) = JVS( 1177 )
   W( 125 ) = JVS( 1178 )
   W( 126 ) = JVS( 1179 )
   W( 127 ) = JVS( 1180 )
   W( 128 ) = JVS( 1181 )
   W( 129 ) = JVS( 1182 )
   W( 130 ) = JVS( 1183 )
   W( 131 ) = JVS( 1184 )
   W( 132 ) = JVS( 1185 )
   W( 133 ) = JVS( 1186 )
   W( 134 ) = JVS( 1187 )
   W( 135 ) = JVS( 1188 )
   W( 136 ) = JVS( 1189 )
  a = -W( 18 ) / JVS(          108  )
  W( 18 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 109 )
  W( 136 ) = W( 136 ) + a*JVS( 110 )
  a = -W( 19 ) / JVS(          111  )
  W( 19 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 112 )
  a = -W( 20 ) / JVS(          113  )
  W( 20 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 114 )
  a = -W( 21 ) / JVS(          116  )
  W( 21 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 117 )
  a = -W( 22 ) / JVS(          119  )
  W( 22 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 120 )
  a = -W( 23 ) / JVS(          121  )
  W( 23 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 122 )
  W( 133 ) = W( 133 ) + a*JVS( 123 )
  a = -W( 24 ) / JVS(          124  )
  W( 24 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 125 )
  a = -W( 25 ) / JVS(          127  )
  W( 25 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 128 )
  a = -W( 30 ) / JVS(          145  )
  W( 30 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 146 )
  W( 90 ) = W( 90 ) + a*JVS( 147 )
  W( 133 ) = W( 133 ) + a*JVS( 148 )
  a = -W( 31 ) / JVS(          150  )
  W( 31 ) = -a
  W( 37 ) = W( 37 ) + a*JVS( 151 )
  W( 136 ) = W( 136 ) + a*JVS( 152 )
  a = -W( 37 ) / JVS(          172  )
  W( 37 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 173 )
  a = -W( 38 ) / JVS(          175  )
  W( 38 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 176 )
  a = -W( 39 ) / JVS(          177  )
  W( 39 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 178 )
  W( 136 ) = W( 136 ) + a*JVS( 179 )
  a = -W( 40 ) / JVS(          180  )
  W( 40 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 181 )
  W( 133 ) = W( 133 ) + a*JVS( 182 )
  a = -W( 41 ) / JVS(          183  )
  W( 41 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 184 )
  W( 129 ) = W( 129 ) + a*JVS( 185 )
  W( 136 ) = W( 136 ) + a*JVS( 186 )
  a = -W( 42 ) / JVS(          187  )
  W( 42 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 188 )
  W( 90 ) = W( 90 ) + a*JVS( 189 )
  W( 133 ) = W( 133 ) + a*JVS( 190 )
  a = -W( 43 ) / JVS(          192  )
  W( 43 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 193 )
  W( 88 ) = W( 88 ) + a*JVS( 194 )
  W( 90 ) = W( 90 ) + a*JVS( 195 )
  W( 133 ) = W( 133 ) + a*JVS( 196 )
  W( 136 ) = W( 136 ) + a*JVS( 197 )
  a = -W( 44 ) / JVS(          198  )
  W( 44 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 199 )
  W( 80 ) = W( 80 ) + a*JVS( 200 )
  W( 126 ) = W( 126 ) + a*JVS( 201 )
  W( 136 ) = W( 136 ) + a*JVS( 202 )
  a = -W( 45 ) / JVS(          203  )
  W( 45 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 204 )
  W( 131 ) = W( 131 ) + a*JVS( 205 )
  W( 136 ) = W( 136 ) + a*JVS( 206 )
  a = -W( 46 ) / JVS(          207  )
  W( 46 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 208 )
  W( 131 ) = W( 131 ) + a*JVS( 209 )
  W( 136 ) = W( 136 ) + a*JVS( 210 )
  a = -W( 48 ) / JVS(          216  )
  W( 48 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 217 )
  W( 136 ) = W( 136 ) + a*JVS( 218 )
  a = -W( 50 ) / JVS(          223  )
  W( 50 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 224 )
  W( 131 ) = W( 131 ) + a*JVS( 225 )
  W( 136 ) = W( 136 ) + a*JVS( 226 )
  a = -W( 51 ) / JVS(          227  )
  W( 51 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 228 )
  W( 131 ) = W( 131 ) + a*JVS( 229 )
  W( 136 ) = W( 136 ) + a*JVS( 230 )
  a = -W( 52 ) / JVS(          231  )
  W( 52 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 232 )
  W( 134 ) = W( 134 ) + a*JVS( 233 )
  W( 136 ) = W( 136 ) + a*JVS( 234 )
  a = -W( 55 ) / JVS(          251  )
  W( 55 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 252 )
  W( 136 ) = W( 136 ) + a*JVS( 253 )
  a = -W( 56 ) / JVS(          254  )
  W( 56 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 255 )
  W( 131 ) = W( 131 ) + a*JVS( 256 )
  W( 133 ) = W( 133 ) + a*JVS( 257 )
  a = -W( 58 ) / JVS(          265  )
  W( 58 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 266 )
  W( 133 ) = W( 133 ) + a*JVS( 267 )
  W( 136 ) = W( 136 ) + a*JVS( 268 )
  a = -W( 60 ) / JVS(          279  )
  W( 60 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 280 )
  W( 133 ) = W( 133 ) + a*JVS( 281 )
  W( 136 ) = W( 136 ) + a*JVS( 282 )
  a = -W( 61 ) / JVS(          286  )
  W( 61 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 287 )
  W( 133 ) = W( 133 ) + a*JVS( 288 )
  W( 136 ) = W( 136 ) + a*JVS( 289 )
  a = -W( 65 ) / JVS(          309  )
  W( 65 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 310 )
  W( 131 ) = W( 131 ) + a*JVS( 311 )
  W( 136 ) = W( 136 ) + a*JVS( 312 )
  a = -W( 66 ) / JVS(          315  )
  W( 66 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 316 )
  W( 133 ) = W( 133 ) + a*JVS( 317 )
  W( 136 ) = W( 136 ) + a*JVS( 318 )
  a = -W( 69 ) / JVS(          332  )
  W( 69 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 333 )
  W( 133 ) = W( 133 ) + a*JVS( 334 )
  W( 136 ) = W( 136 ) + a*JVS( 335 )
  a = -W( 70 ) / JVS(          336  )
  W( 70 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 337 )
  W( 131 ) = W( 131 ) + a*JVS( 338 )
  W( 136 ) = W( 136 ) + a*JVS( 339 )
  a = -W( 72 ) / JVS(          349  )
  W( 72 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 350 )
  W( 109 ) = W( 109 ) + a*JVS( 351 )
  W( 132 ) = W( 132 ) + a*JVS( 352 )
  W( 136 ) = W( 136 ) + a*JVS( 353 )
  a = -W( 73 ) / JVS(          354  )
  W( 73 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 355 )
  W( 131 ) = W( 131 ) + a*JVS( 356 )
  W( 133 ) = W( 133 ) + a*JVS( 357 )
  W( 136 ) = W( 136 ) + a*JVS( 358 )
  a = -W( 74 ) / JVS(          362  )
  W( 74 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 363 )
  W( 133 ) = W( 133 ) + a*JVS( 364 )
  W( 134 ) = W( 134 ) + a*JVS( 365 )
  W( 136 ) = W( 136 ) + a*JVS( 366 )
  a = -W( 75 ) / JVS(          367  )
  W( 75 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 368 )
  W( 131 ) = W( 131 ) + a*JVS( 369 )
  W( 136 ) = W( 136 ) + a*JVS( 370 )
  a = -W( 76 ) / JVS(          371  )
  W( 76 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 372 )
  W( 131 ) = W( 131 ) + a*JVS( 373 )
  W( 136 ) = W( 136 ) + a*JVS( 374 )
  a = -W( 77 ) / JVS(          375  )
  W( 77 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 376 )
  W( 134 ) = W( 134 ) + a*JVS( 377 )
  W( 136 ) = W( 136 ) + a*JVS( 378 )
  a = -W( 78 ) / JVS(          382  )
  W( 78 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 383 )
  W( 83 ) = W( 83 ) + a*JVS( 384 )
  W( 89 ) = W( 89 ) + a*JVS( 385 )
  W( 100 ) = W( 100 ) + a*JVS( 386 )
  W( 101 ) = W( 101 ) + a*JVS( 387 )
  W( 104 ) = W( 104 ) + a*JVS( 388 )
  W( 126 ) = W( 126 ) + a*JVS( 389 )
  W( 131 ) = W( 131 ) + a*JVS( 390 )
  W( 132 ) = W( 132 ) + a*JVS( 391 )
  W( 133 ) = W( 133 ) + a*JVS( 392 )
  W( 136 ) = W( 136 ) + a*JVS( 393 )
  a = -W( 79 ) / JVS(          394  )
  W( 79 ) = -a
  W( 120 ) = W( 120 ) + a*JVS( 395 )
  W( 131 ) = W( 131 ) + a*JVS( 396 )
  W( 136 ) = W( 136 ) + a*JVS( 397 )
  a = -W( 80 ) / JVS(          401  )
  W( 80 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 402 )
  W( 126 ) = W( 126 ) + a*JVS( 403 )
  W( 132 ) = W( 132 ) + a*JVS( 404 )
  W( 136 ) = W( 136 ) + a*JVS( 405 )
  a = -W( 81 ) / JVS(          409  )
  W( 81 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 410 )
  W( 132 ) = W( 132 ) + a*JVS( 411 )
  W( 133 ) = W( 133 ) + a*JVS( 412 )
  W( 134 ) = W( 134 ) + a*JVS( 413 )
  W( 136 ) = W( 136 ) + a*JVS( 414 )
  a = -W( 82 ) / JVS(          415  )
  W( 82 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 416 )
  W( 136 ) = W( 136 ) + a*JVS( 417 )
  a = -W( 83 ) / JVS(          418  )
  W( 83 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 84 ) / JVS(          421  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 422 )
  W( 131 ) = W( 131 ) + a*JVS( 423 )
  W( 136 ) = W( 136 ) + a*JVS( 424 )
  a = -W( 85 ) / JVS(          428  )
  W( 85 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 429 )
  W( 131 ) = W( 131 ) + a*JVS( 430 )
  W( 133 ) = W( 133 ) + a*JVS( 431 )
  W( 136 ) = W( 136 ) + a*JVS( 432 )
  a = -W( 86 ) / JVS(          434  )
  W( 86 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 435 )
  W( 131 ) = W( 131 ) + a*JVS( 436 )
  W( 133 ) = W( 133 ) + a*JVS( 437 )
  W( 136 ) = W( 136 ) + a*JVS( 438 )
  a = -W( 87 ) / JVS(          441  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 90 ) = W( 90 ) + a*JVS( 443 )
  W( 131 ) = W( 131 ) + a*JVS( 444 )
  W( 133 ) = W( 133 ) + a*JVS( 445 )
  W( 134 ) = W( 134 ) + a*JVS( 446 )
  W( 136 ) = W( 136 ) + a*JVS( 447 )
  a = -W( 88 ) / JVS(          450  )
  W( 88 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 451 )
  W( 133 ) = W( 133 ) + a*JVS( 452 )
  W( 136 ) = W( 136 ) + a*JVS( 453 )
  a = -W( 89 ) / JVS(          454  )
  W( 89 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 455 )
  W( 135 ) = W( 135 ) + a*JVS( 456 )
  W( 136 ) = W( 136 ) + a*JVS( 457 )
  a = -W( 90 ) / JVS(          460  )
  W( 90 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 461 )
  W( 133 ) = W( 133 ) + a*JVS( 462 )
  W( 136 ) = W( 136 ) + a*JVS( 463 )
  a = -W( 91 ) / JVS(          467  )
  W( 91 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 468 )
  W( 133 ) = W( 133 ) + a*JVS( 469 )
  W( 134 ) = W( 134 ) + a*JVS( 470 )
  W( 136 ) = W( 136 ) + a*JVS( 471 )
  a = -W( 92 ) / JVS(          473  )
  W( 92 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 474 )
  W( 126 ) = W( 126 ) + a*JVS( 475 )
  W( 131 ) = W( 131 ) + a*JVS( 476 )
  W( 132 ) = W( 132 ) + a*JVS( 477 )
  W( 134 ) = W( 134 ) + a*JVS( 478 )
  W( 136 ) = W( 136 ) + a*JVS( 479 )
  a = -W( 93 ) / JVS(          482  )
  W( 93 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 483 )
  W( 131 ) = W( 131 ) + a*JVS( 484 )
  W( 133 ) = W( 133 ) + a*JVS( 485 )
  W( 136 ) = W( 136 ) + a*JVS( 486 )
  a = -W( 94 ) / JVS(          490  )
  W( 94 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 491 )
  W( 133 ) = W( 133 ) + a*JVS( 492 )
  W( 134 ) = W( 134 ) + a*JVS( 493 )
  W( 136 ) = W( 136 ) + a*JVS( 494 )
  a = -W( 95 ) / JVS(          495  )
  W( 95 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 496 )
  W( 111 ) = W( 111 ) + a*JVS( 497 )
  W( 118 ) = W( 118 ) + a*JVS( 498 )
  W( 119 ) = W( 119 ) + a*JVS( 499 )
  W( 120 ) = W( 120 ) + a*JVS( 500 )
  W( 121 ) = W( 121 ) + a*JVS( 501 )
  W( 123 ) = W( 123 ) + a*JVS( 502 )
  W( 125 ) = W( 125 ) + a*JVS( 503 )
  W( 127 ) = W( 127 ) + a*JVS( 504 )
  W( 129 ) = W( 129 ) + a*JVS( 505 )
  W( 136 ) = W( 136 ) + a*JVS( 506 )
  a = -W( 96 ) / JVS(          516  )
  W( 96 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 517 )
  W( 108 ) = W( 108 ) + a*JVS( 518 )
  W( 123 ) = W( 123 ) + a*JVS( 519 )
  W( 127 ) = W( 127 ) + a*JVS( 520 )
  W( 128 ) = W( 128 ) + a*JVS( 521 )
  W( 129 ) = W( 129 ) + a*JVS( 522 )
  W( 131 ) = W( 131 ) + a*JVS( 523 )
  W( 133 ) = W( 133 ) + a*JVS( 524 )
  W( 134 ) = W( 134 ) + a*JVS( 525 )
  W( 135 ) = W( 135 ) + a*JVS( 526 )
  W( 136 ) = W( 136 ) + a*JVS( 527 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 533 )
  W( 129 ) = W( 129 ) + a*JVS( 534 )
  W( 131 ) = W( 131 ) + a*JVS( 535 )
  W( 133 ) = W( 133 ) + a*JVS( 536 )
  W( 136 ) = W( 136 ) + a*JVS( 537 )
  a = -W( 98 ) / JVS(          542  )
  W( 98 ) = -a
  W( 100 ) = W( 100 ) + a*JVS( 543 )
  W( 101 ) = W( 101 ) + a*JVS( 544 )
  W( 104 ) = W( 104 ) + a*JVS( 545 )
  W( 131 ) = W( 131 ) + a*JVS( 546 )
  W( 132 ) = W( 132 ) + a*JVS( 547 )
  W( 133 ) = W( 133 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 99 ) / JVS(          550  )
  W( 99 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 133 ) = W( 133 ) + a*JVS( 553 )
  W( 135 ) = W( 135 ) + a*JVS( 554 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 102 ) / JVS(          564  )
  W( 102 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 565 )
  W( 132 ) = W( 132 ) + a*JVS( 566 )
  W( 133 ) = W( 133 ) + a*JVS( 567 )
  W( 135 ) = W( 135 ) + a*JVS( 568 )
  W( 136 ) = W( 136 ) + a*JVS( 569 )
  a = -W( 103 ) / JVS(          573  )
  W( 103 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 574 )
  W( 129 ) = W( 129 ) + a*JVS( 575 )
  W( 131 ) = W( 131 ) + a*JVS( 576 )
  W( 132 ) = W( 132 ) + a*JVS( 577 )
  W( 133 ) = W( 133 ) + a*JVS( 578 )
  W( 135 ) = W( 135 ) + a*JVS( 579 )
  W( 136 ) = W( 136 ) + a*JVS( 580 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  a = -W( 107 ) / JVS(          640  )
  W( 107 ) = -a
  W( 108 ) = W( 108 ) + a*JVS( 641 )
  W( 109 ) = W( 109 ) + a*JVS( 642 )
  W( 110 ) = W( 110 ) + a*JVS( 643 )
  W( 111 ) = W( 111 ) + a*JVS( 644 )
  W( 112 ) = W( 112 ) + a*JVS( 645 )
  W( 114 ) = W( 114 ) + a*JVS( 646 )
  W( 115 ) = W( 115 ) + a*JVS( 647 )
  W( 116 ) = W( 116 ) + a*JVS( 648 )
  W( 117 ) = W( 117 ) + a*JVS( 649 )
  W( 120 ) = W( 120 ) + a*JVS( 650 )
  W( 123 ) = W( 123 ) + a*JVS( 651 )
  W( 124 ) = W( 124 ) + a*JVS( 652 )
  W( 125 ) = W( 125 ) + a*JVS( 653 )
  W( 126 ) = W( 126 ) + a*JVS( 654 )
  W( 127 ) = W( 127 ) + a*JVS( 655 )
  W( 128 ) = W( 128 ) + a*JVS( 656 )
  W( 129 ) = W( 129 ) + a*JVS( 657 )
  W( 131 ) = W( 131 ) + a*JVS( 658 )
  W( 132 ) = W( 132 ) + a*JVS( 659 )
  W( 133 ) = W( 133 ) + a*JVS( 660 )
  W( 134 ) = W( 134 ) + a*JVS( 661 )
  W( 135 ) = W( 135 ) + a*JVS( 662 )
  W( 136 ) = W( 136 ) + a*JVS( 663 )
  a = -W( 108 ) / JVS(          668  )
  W( 108 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 669 )
  W( 120 ) = W( 120 ) + a*JVS( 670 )
  W( 123 ) = W( 123 ) + a*JVS( 671 )
  W( 125 ) = W( 125 ) + a*JVS( 672 )
  W( 127 ) = W( 127 ) + a*JVS( 673 )
  W( 128 ) = W( 128 ) + a*JVS( 674 )
  W( 129 ) = W( 129 ) + a*JVS( 675 )
  W( 131 ) = W( 131 ) + a*JVS( 676 )
  W( 132 ) = W( 132 ) + a*JVS( 677 )
  W( 133 ) = W( 133 ) + a*JVS( 678 )
  W( 135 ) = W( 135 ) + a*JVS( 679 )
  W( 136 ) = W( 136 ) + a*JVS( 680 )
  a = -W( 109 ) / JVS(          682  )
  W( 109 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 683 )
  W( 132 ) = W( 132 ) + a*JVS( 684 )
  W( 135 ) = W( 135 ) + a*JVS( 685 )
  W( 136 ) = W( 136 ) + a*JVS( 686 )
  a = -W( 110 ) / JVS(          693  )
  W( 110 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 694 )
  W( 119 ) = W( 119 ) + a*JVS( 695 )
  W( 129 ) = W( 129 ) + a*JVS( 696 )
  W( 131 ) = W( 131 ) + a*JVS( 697 )
  W( 132 ) = W( 132 ) + a*JVS( 698 )
  W( 133 ) = W( 133 ) + a*JVS( 699 )
  W( 135 ) = W( 135 ) + a*JVS( 700 )
  W( 136 ) = W( 136 ) + a*JVS( 701 )
  a = -W( 111 ) / JVS(          704  )
  W( 111 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 705 )
  W( 131 ) = W( 131 ) + a*JVS( 706 )
  W( 132 ) = W( 132 ) + a*JVS( 707 )
  W( 133 ) = W( 133 ) + a*JVS( 708 )
  W( 135 ) = W( 135 ) + a*JVS( 709 )
  W( 136 ) = W( 136 ) + a*JVS( 710 )
  a = -W( 112 ) / JVS(          722  )
  W( 112 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 723 )
  W( 124 ) = W( 124 ) + a*JVS( 724 )
  W( 129 ) = W( 129 ) + a*JVS( 725 )
  W( 131 ) = W( 131 ) + a*JVS( 726 )
  W( 132 ) = W( 132 ) + a*JVS( 727 )
  W( 133 ) = W( 133 ) + a*JVS( 728 )
  W( 135 ) = W( 135 ) + a*JVS( 729 )
  W( 136 ) = W( 136 ) + a*JVS( 730 )
  a = -W( 113 ) / JVS(          738  )
  W( 113 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 739 )
  W( 132 ) = W( 132 ) + a*JVS( 740 )
  W( 133 ) = W( 133 ) + a*JVS( 741 )
  W( 135 ) = W( 135 ) + a*JVS( 742 )
  W( 136 ) = W( 136 ) + a*JVS( 743 )
  a = -W( 114 ) / JVS(          751  )
  W( 114 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 752 )
  W( 119 ) = W( 119 ) + a*JVS( 753 )
  W( 120 ) = W( 120 ) + a*JVS( 754 )
  W( 129 ) = W( 129 ) + a*JVS( 755 )
  W( 131 ) = W( 131 ) + a*JVS( 756 )
  W( 132 ) = W( 132 ) + a*JVS( 757 )
  W( 133 ) = W( 133 ) + a*JVS( 758 )
  W( 135 ) = W( 135 ) + a*JVS( 759 )
  W( 136 ) = W( 136 ) + a*JVS( 760 )
  a = -W( 115 ) / JVS(          764  )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 765 )
  W( 128 ) = W( 128 ) + a*JVS( 766 )
  W( 129 ) = W( 129 ) + a*JVS( 767 )
  W( 131 ) = W( 131 ) + a*JVS( 768 )
  W( 132 ) = W( 132 ) + a*JVS( 769 )
  W( 133 ) = W( 133 ) + a*JVS( 770 )
  W( 135 ) = W( 135 ) + a*JVS( 771 )
  W( 136 ) = W( 136 ) + a*JVS( 772 )
  a = -W( 116 ) / JVS(          775  )
  W( 116 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 776 )
  W( 119 ) = W( 119 ) + a*JVS( 777 )
  W( 120 ) = W( 120 ) + a*JVS( 778 )
  W( 125 ) = W( 125 ) + a*JVS( 779 )
  W( 127 ) = W( 127 ) + a*JVS( 780 )
  W( 129 ) = W( 129 ) + a*JVS( 781 )
  W( 131 ) = W( 131 ) + a*JVS( 782 )
  W( 132 ) = W( 132 ) + a*JVS( 783 )
  W( 133 ) = W( 133 ) + a*JVS( 784 )
  W( 135 ) = W( 135 ) + a*JVS( 785 )
  W( 136 ) = W( 136 ) + a*JVS( 786 )
  a = -W( 117 ) / JVS(          803  )
  W( 117 ) = -a
  W( 121 ) = W( 121 ) + a*JVS( 804 )
  W( 122 ) = W( 122 ) + a*JVS( 805 )
  W( 123 ) = W( 123 ) + a*JVS( 806 )
  W( 124 ) = W( 124 ) + a*JVS( 807 )
  W( 125 ) = W( 125 ) + a*JVS( 808 )
  W( 127 ) = W( 127 ) + a*JVS( 809 )
  W( 128 ) = W( 128 ) + a*JVS( 810 )
  W( 129 ) = W( 129 ) + a*JVS( 811 )
  W( 131 ) = W( 131 ) + a*JVS( 812 )
  W( 132 ) = W( 132 ) + a*JVS( 813 )
  W( 133 ) = W( 133 ) + a*JVS( 814 )
  W( 134 ) = W( 134 ) + a*JVS( 815 )
  W( 135 ) = W( 135 ) + a*JVS( 816 )
  W( 136 ) = W( 136 ) + a*JVS( 817 )
  a = -W( 118 ) / JVS(          824  )
  W( 118 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 825 )
  W( 131 ) = W( 131 ) + a*JVS( 826 )
  W( 132 ) = W( 132 ) + a*JVS( 827 )
  W( 133 ) = W( 133 ) + a*JVS( 828 )
  W( 135 ) = W( 135 ) + a*JVS( 829 )
  W( 136 ) = W( 136 ) + a*JVS( 830 )
  a = -W( 119 ) / JVS(          836  )
  W( 119 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 837 )
  W( 131 ) = W( 131 ) + a*JVS( 838 )
  W( 132 ) = W( 132 ) + a*JVS( 839 )
  W( 133 ) = W( 133 ) + a*JVS( 840 )
  W( 135 ) = W( 135 ) + a*JVS( 841 )
  W( 136 ) = W( 136 ) + a*JVS( 842 )
  a = -W( 120 ) / JVS(          847  )
  W( 120 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 848 )
  W( 131 ) = W( 131 ) + a*JVS( 849 )
  W( 132 ) = W( 132 ) + a*JVS( 850 )
  W( 133 ) = W( 133 ) + a*JVS( 851 )
  W( 135 ) = W( 135 ) + a*JVS( 852 )
  W( 136 ) = W( 136 ) + a*JVS( 853 )
  a = -W( 121 ) / JVS(          867  )
  W( 121 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 868 )
  W( 131 ) = W( 131 ) + a*JVS( 869 )
  W( 132 ) = W( 132 ) + a*JVS( 870 )
  W( 133 ) = W( 133 ) + a*JVS( 871 )
  W( 135 ) = W( 135 ) + a*JVS( 872 )
  W( 136 ) = W( 136 ) + a*JVS( 873 )
  a = -W( 122 ) / JVS(          879  )
  W( 122 ) = -a
  W( 123 ) = W( 123 ) + a*JVS( 880 )
  W( 124 ) = W( 124 ) + a*JVS( 881 )
  W( 125 ) = W( 125 ) + a*JVS( 882 )
  W( 127 ) = W( 127 ) + a*JVS( 883 )
  W( 128 ) = W( 128 ) + a*JVS( 884 )
  W( 129 ) = W( 129 ) + a*JVS( 885 )
  W( 130 ) = W( 130 ) + a*JVS( 886 )
  W( 131 ) = W( 131 ) + a*JVS( 887 )
  W( 132 ) = W( 132 ) + a*JVS( 888 )
  W( 133 ) = W( 133 ) + a*JVS( 889 )
  W( 134 ) = W( 134 ) + a*JVS( 890 )
  W( 135 ) = W( 135 ) + a*JVS( 891 )
  W( 136 ) = W( 136 ) + a*JVS( 892 )
  a = -W( 123 ) / JVS(          900  )
  W( 123 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 901 )
  W( 127 ) = W( 127 ) + a*JVS( 902 )
  W( 128 ) = W( 128 ) + a*JVS( 903 )
  W( 129 ) = W( 129 ) + a*JVS( 904 )
  W( 131 ) = W( 131 ) + a*JVS( 905 )
  W( 132 ) = W( 132 ) + a*JVS( 906 )
  W( 133 ) = W( 133 ) + a*JVS( 907 )
  W( 135 ) = W( 135 ) + a*JVS( 908 )
  W( 136 ) = W( 136 ) + a*JVS( 909 )
  a = -W( 124 ) / JVS(          913  )
  W( 124 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 914 )
  W( 128 ) = W( 128 ) + a*JVS( 915 )
  W( 129 ) = W( 129 ) + a*JVS( 916 )
  W( 131 ) = W( 131 ) + a*JVS( 917 )
  W( 132 ) = W( 132 ) + a*JVS( 918 )
  W( 133 ) = W( 133 ) + a*JVS( 919 )
  W( 135 ) = W( 135 ) + a*JVS( 920 )
  W( 136 ) = W( 136 ) + a*JVS( 921 )
  a = -W( 125 ) / JVS(          925  )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 926 )
  W( 128 ) = W( 128 ) + a*JVS( 927 )
  W( 129 ) = W( 129 ) + a*JVS( 928 )
  W( 131 ) = W( 131 ) + a*JVS( 929 )
  W( 132 ) = W( 132 ) + a*JVS( 930 )
  W( 133 ) = W( 133 ) + a*JVS( 931 )
  W( 135 ) = W( 135 ) + a*JVS( 932 )
  W( 136 ) = W( 136 ) + a*JVS( 933 )
  a = -W( 126 ) / JVS(          980  )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 981 )
  W( 128 ) = W( 128 ) + a*JVS( 982 )
  W( 129 ) = W( 129 ) + a*JVS( 983 )
  W( 130 ) = W( 130 ) + a*JVS( 984 )
  W( 131 ) = W( 131 ) + a*JVS( 985 )
  W( 132 ) = W( 132 ) + a*JVS( 986 )
  W( 133 ) = W( 133 ) + a*JVS( 987 )
  W( 134 ) = W( 134 ) + a*JVS( 988 )
  W( 135 ) = W( 135 ) + a*JVS( 989 )
  W( 136 ) = W( 136 ) + a*JVS( 990 )
  a = -W( 127 ) / JVS(          993  )
  W( 127 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 994 )
  W( 129 ) = W( 129 ) + a*JVS( 995 )
  W( 131 ) = W( 131 ) + a*JVS( 996 )
  W( 132 ) = W( 132 ) + a*JVS( 997 )
  W( 133 ) = W( 133 ) + a*JVS( 998 )
  W( 135 ) = W( 135 ) + a*JVS( 999 )
  W( 136 ) = W( 136 ) + a*JVS( 1000 )
  a = -W( 128 ) / JVS(         1034  )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 1035 )
  W( 130 ) = W( 130 ) + a*JVS( 1036 )
  W( 131 ) = W( 131 ) + a*JVS( 1037 )
  W( 132 ) = W( 132 ) + a*JVS( 1038 )
  W( 133 ) = W( 133 ) + a*JVS( 1039 )
  W( 134 ) = W( 134 ) + a*JVS( 1040 )
  W( 135 ) = W( 135 ) + a*JVS( 1041 )
  W( 136 ) = W( 136 ) + a*JVS( 1042 )
  a = -W( 129 ) / JVS(         1070  )
  W( 129 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 1071 )
  W( 131 ) = W( 131 ) + a*JVS( 1072 )
  W( 132 ) = W( 132 ) + a*JVS( 1073 )
  W( 133 ) = W( 133 ) + a*JVS( 1074 )
  W( 134 ) = W( 134 ) + a*JVS( 1075 )
  W( 135 ) = W( 135 ) + a*JVS( 1076 )
  W( 136 ) = W( 136 ) + a*JVS( 1077 )
  a = -W( 130 ) / JVS(         1086  )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 1087 )
  W( 132 ) = W( 132 ) + a*JVS( 1088 )
  W( 133 ) = W( 133 ) + a*JVS( 1089 )
  W( 134 ) = W( 134 ) + a*JVS( 1090 )
  W( 135 ) = W( 135 ) + a*JVS( 1091 )
  W( 136 ) = W( 136 ) + a*JVS( 1092 )
  JVS( 1093) = W( 18 )
  JVS( 1094) = W( 19 )
  JVS( 1095) = W( 20 )
  JVS( 1096) = W( 21 )
  JVS( 1097) = W( 22 )
  JVS( 1098) = W( 23 )
  JVS( 1099) = W( 24 )
  JVS( 1100) = W( 25 )
  JVS( 1101) = W( 30 )
  JVS( 1102) = W( 31 )
  JVS( 1103) = W( 37 )
  JVS( 1104) = W( 38 )
  JVS( 1105) = W( 39 )
  JVS( 1106) = W( 40 )
  JVS( 1107) = W( 41 )
  JVS( 1108) = W( 42 )
  JVS( 1109) = W( 43 )
  JVS( 1110) = W( 44 )
  JVS( 1111) = W( 45 )
  JVS( 1112) = W( 46 )
  JVS( 1113) = W( 48 )
  JVS( 1114) = W( 50 )
  JVS( 1115) = W( 51 )
  JVS( 1116) = W( 52 )
  JVS( 1117) = W( 55 )
  JVS( 1118) = W( 56 )
  JVS( 1119) = W( 58 )
  JVS( 1120) = W( 60 )
  JVS( 1121) = W( 61 )
  JVS( 1122) = W( 65 )
  JVS( 1123) = W( 66 )
  JVS( 1124) = W( 69 )
  JVS( 1125) = W( 70 )
  JVS( 1126) = W( 72 )
  JVS( 1127) = W( 73 )
  JVS( 1128) = W( 74 )
  JVS( 1129) = W( 75 )
  JVS( 1130) = W( 76 )
  JVS( 1131) = W( 77 )
  JVS( 1132) = W( 78 )
  JVS( 1133) = W( 79 )
  JVS( 1134) = W( 80 )
  JVS( 1135) = W( 81 )
  JVS( 1136) = W( 82 )
  JVS( 1137) = W( 83 )
  JVS( 1138) = W( 84 )
  JVS( 1139) = W( 85 )
  JVS( 1140) = W( 86 )
  JVS( 1141) = W( 87 )
  JVS( 1142) = W( 88 )
  JVS( 1143) = W( 89 )
  JVS( 1144) = W( 90 )
  JVS( 1145) = W( 91 )
  JVS( 1146) = W( 92 )
  JVS( 1147) = W( 93 )
  JVS( 1148) = W( 94 )
  JVS( 1149) = W( 95 )
  JVS( 1150) = W( 96 )
  JVS( 1151) = W( 97 )
  JVS( 1152) = W( 98 )
  JVS( 1153) = W( 99 )
  JVS( 1154) = W( 100 )
  JVS( 1155) = W( 101 )
  JVS( 1156) = W( 102 )
  JVS( 1157) = W( 103 )
  JVS( 1158) = W( 104 )
  JVS( 1159) = W( 106 )
  JVS( 1160) = W( 107 )
  JVS( 1161) = W( 108 )
  JVS( 1162) = W( 109 )
  JVS( 1163) = W( 110 )
  JVS( 1164) = W( 111 )
  JVS( 1165) = W( 112 )
  JVS( 1166) = W( 113 )
  JVS( 1167) = W( 114 )
  JVS( 1168) = W( 115 )
  JVS( 1169) = W( 116 )
  JVS( 1170) = W( 117 )
  JVS( 1171) = W( 118 )
  JVS( 1172) = W( 119 )
  JVS( 1173) = W( 120 )
  JVS( 1174) = W( 121 )
  JVS( 1175) = W( 122 )
  JVS( 1176) = W( 123 )
  JVS( 1177) = W( 124 )
  JVS( 1178) = W( 125 )
  JVS( 1179) = W( 126 )
  JVS( 1180) = W( 127 )
  JVS( 1181) = W( 128 )
  JVS( 1182) = W( 129 )
  JVS( 1183) = W( 130 )
  JVS( 1184) = W( 131 )
  JVS( 1185) = W( 132 )
  JVS( 1186) = W( 133 )
  JVS( 1187) = W( 134 )
  JVS( 1188) = W( 135 )
  JVS( 1189) = W( 136 )
  IF ( ABS(  JVS( 1210 )) < TINY(a) ) THEN
         IER = 132                                     
         RETURN
  END IF
   W( 55 ) = JVS( 1190 )
   W( 59 ) = JVS( 1191 )
   W( 81 ) = JVS( 1192 )
   W( 82 ) = JVS( 1193 )
   W( 83 ) = JVS( 1194 )
   W( 89 ) = JVS( 1195 )
   W( 92 ) = JVS( 1196 )
   W( 100 ) = JVS( 1197 )
   W( 101 ) = JVS( 1198 )
   W( 104 ) = JVS( 1199 )
   W( 106 ) = JVS( 1200 )
   W( 109 ) = JVS( 1201 )
   W( 115 ) = JVS( 1202 )
   W( 124 ) = JVS( 1203 )
   W( 126 ) = JVS( 1204 )
   W( 127 ) = JVS( 1205 )
   W( 128 ) = JVS( 1206 )
   W( 129 ) = JVS( 1207 )
   W( 130 ) = JVS( 1208 )
   W( 131 ) = JVS( 1209 )
   W( 132 ) = JVS( 1210 )
   W( 133 ) = JVS( 1211 )
   W( 134 ) = JVS( 1212 )
   W( 135 ) = JVS( 1213 )
   W( 136 ) = JVS( 1214 )
  a = -W( 55 ) / JVS(          251  )
  W( 55 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 252 )
  W( 136 ) = W( 136 ) + a*JVS( 253 )
  a = -W( 59 ) / JVS(          272  )
  W( 59 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 273 )
  W( 132 ) = W( 132 ) + a*JVS( 274 )
  W( 133 ) = W( 133 ) + a*JVS( 275 )
  W( 134 ) = W( 134 ) + a*JVS( 276 )
  W( 136 ) = W( 136 ) + a*JVS( 277 )
  a = -W( 81 ) / JVS(          409  )
  W( 81 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 410 )
  W( 132 ) = W( 132 ) + a*JVS( 411 )
  W( 133 ) = W( 133 ) + a*JVS( 412 )
  W( 134 ) = W( 134 ) + a*JVS( 413 )
  W( 136 ) = W( 136 ) + a*JVS( 414 )
  a = -W( 82 ) / JVS(          415  )
  W( 82 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 416 )
  W( 136 ) = W( 136 ) + a*JVS( 417 )
  a = -W( 83 ) / JVS(          418  )
  W( 83 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 89 ) / JVS(          454  )
  W( 89 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 455 )
  W( 135 ) = W( 135 ) + a*JVS( 456 )
  W( 136 ) = W( 136 ) + a*JVS( 457 )
  a = -W( 92 ) / JVS(          473  )
  W( 92 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 474 )
  W( 126 ) = W( 126 ) + a*JVS( 475 )
  W( 131 ) = W( 131 ) + a*JVS( 476 )
  W( 132 ) = W( 132 ) + a*JVS( 477 )
  W( 134 ) = W( 134 ) + a*JVS( 478 )
  W( 136 ) = W( 136 ) + a*JVS( 479 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  a = -W( 109 ) / JVS(          682  )
  W( 109 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 683 )
  W( 132 ) = W( 132 ) + a*JVS( 684 )
  W( 135 ) = W( 135 ) + a*JVS( 685 )
  W( 136 ) = W( 136 ) + a*JVS( 686 )
  a = -W( 115 ) / JVS(          764  )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 765 )
  W( 128 ) = W( 128 ) + a*JVS( 766 )
  W( 129 ) = W( 129 ) + a*JVS( 767 )
  W( 131 ) = W( 131 ) + a*JVS( 768 )
  W( 132 ) = W( 132 ) + a*JVS( 769 )
  W( 133 ) = W( 133 ) + a*JVS( 770 )
  W( 135 ) = W( 135 ) + a*JVS( 771 )
  W( 136 ) = W( 136 ) + a*JVS( 772 )
  a = -W( 124 ) / JVS(          913  )
  W( 124 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 914 )
  W( 128 ) = W( 128 ) + a*JVS( 915 )
  W( 129 ) = W( 129 ) + a*JVS( 916 )
  W( 131 ) = W( 131 ) + a*JVS( 917 )
  W( 132 ) = W( 132 ) + a*JVS( 918 )
  W( 133 ) = W( 133 ) + a*JVS( 919 )
  W( 135 ) = W( 135 ) + a*JVS( 920 )
  W( 136 ) = W( 136 ) + a*JVS( 921 )
  a = -W( 126 ) / JVS(          980  )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 981 )
  W( 128 ) = W( 128 ) + a*JVS( 982 )
  W( 129 ) = W( 129 ) + a*JVS( 983 )
  W( 130 ) = W( 130 ) + a*JVS( 984 )
  W( 131 ) = W( 131 ) + a*JVS( 985 )
  W( 132 ) = W( 132 ) + a*JVS( 986 )
  W( 133 ) = W( 133 ) + a*JVS( 987 )
  W( 134 ) = W( 134 ) + a*JVS( 988 )
  W( 135 ) = W( 135 ) + a*JVS( 989 )
  W( 136 ) = W( 136 ) + a*JVS( 990 )
  a = -W( 127 ) / JVS(          993  )
  W( 127 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 994 )
  W( 129 ) = W( 129 ) + a*JVS( 995 )
  W( 131 ) = W( 131 ) + a*JVS( 996 )
  W( 132 ) = W( 132 ) + a*JVS( 997 )
  W( 133 ) = W( 133 ) + a*JVS( 998 )
  W( 135 ) = W( 135 ) + a*JVS( 999 )
  W( 136 ) = W( 136 ) + a*JVS( 1000 )
  a = -W( 128 ) / JVS(         1034  )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 1035 )
  W( 130 ) = W( 130 ) + a*JVS( 1036 )
  W( 131 ) = W( 131 ) + a*JVS( 1037 )
  W( 132 ) = W( 132 ) + a*JVS( 1038 )
  W( 133 ) = W( 133 ) + a*JVS( 1039 )
  W( 134 ) = W( 134 ) + a*JVS( 1040 )
  W( 135 ) = W( 135 ) + a*JVS( 1041 )
  W( 136 ) = W( 136 ) + a*JVS( 1042 )
  a = -W( 129 ) / JVS(         1070  )
  W( 129 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 1071 )
  W( 131 ) = W( 131 ) + a*JVS( 1072 )
  W( 132 ) = W( 132 ) + a*JVS( 1073 )
  W( 133 ) = W( 133 ) + a*JVS( 1074 )
  W( 134 ) = W( 134 ) + a*JVS( 1075 )
  W( 135 ) = W( 135 ) + a*JVS( 1076 )
  W( 136 ) = W( 136 ) + a*JVS( 1077 )
  a = -W( 130 ) / JVS(         1086  )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 1087 )
  W( 132 ) = W( 132 ) + a*JVS( 1088 )
  W( 133 ) = W( 133 ) + a*JVS( 1089 )
  W( 134 ) = W( 134 ) + a*JVS( 1090 )
  W( 135 ) = W( 135 ) + a*JVS( 1091 )
  W( 136 ) = W( 136 ) + a*JVS( 1092 )
  a = -W( 131 ) / JVS(         1184  )
  W( 131 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 1185 )
  W( 133 ) = W( 133 ) + a*JVS( 1186 )
  W( 134 ) = W( 134 ) + a*JVS( 1187 )
  W( 135 ) = W( 135 ) + a*JVS( 1188 )
  W( 136 ) = W( 136 ) + a*JVS( 1189 )
  JVS( 1190) = W( 55 )
  JVS( 1191) = W( 59 )
  JVS( 1192) = W( 81 )
  JVS( 1193) = W( 82 )
  JVS( 1194) = W( 83 )
  JVS( 1195) = W( 89 )
  JVS( 1196) = W( 92 )
  JVS( 1197) = W( 100 )
  JVS( 1198) = W( 101 )
  JVS( 1199) = W( 104 )
  JVS( 1200) = W( 106 )
  JVS( 1201) = W( 109 )
  JVS( 1202) = W( 115 )
  JVS( 1203) = W( 124 )
  JVS( 1204) = W( 126 )
  JVS( 1205) = W( 127 )
  JVS( 1206) = W( 128 )
  JVS( 1207) = W( 129 )
  JVS( 1208) = W( 130 )
  JVS( 1209) = W( 131 )
  JVS( 1210) = W( 132 )
  JVS( 1211) = W( 133 )
  JVS( 1212) = W( 134 )
  JVS( 1213) = W( 135 )
  JVS( 1214) = W( 136 )
  IF ( ABS(  JVS( 1259 )) < TINY(a) ) THEN
         IER = 133                                     
         RETURN
  END IF
   W( 1 ) = JVS( 1215 )
   W( 15 ) = JVS( 1216 )
   W( 48 ) = JVS( 1217 )
   W( 56 ) = JVS( 1218 )
   W( 58 ) = JVS( 1219 )
   W( 60 ) = JVS( 1220 )
   W( 61 ) = JVS( 1221 )
   W( 66 ) = JVS( 1222 )
   W( 69 ) = JVS( 1223 )
   W( 73 ) = JVS( 1224 )
   W( 74 ) = JVS( 1225 )
   W( 80 ) = JVS( 1226 )
   W( 81 ) = JVS( 1227 )
   W( 85 ) = JVS( 1228 )
   W( 86 ) = JVS( 1229 )
   W( 87 ) = JVS( 1230 )
   W( 88 ) = JVS( 1231 )
   W( 90 ) = JVS( 1232 )
   W( 91 ) = JVS( 1233 )
   W( 92 ) = JVS( 1234 )
   W( 93 ) = JVS( 1235 )
   W( 94 ) = JVS( 1236 )
   W( 97 ) = JVS( 1237 )
   W( 99 ) = JVS( 1238 )
   W( 102 ) = JVS( 1239 )
   W( 103 ) = JVS( 1240 )
   W( 106 ) = JVS( 1241 )
   W( 109 ) = JVS( 1242 )
   W( 111 ) = JVS( 1243 )
   W( 113 ) = JVS( 1244 )
   W( 118 ) = JVS( 1245 )
   W( 119 ) = JVS( 1246 )
   W( 120 ) = JVS( 1247 )
   W( 121 ) = JVS( 1248 )
   W( 123 ) = JVS( 1249 )
   W( 124 ) = JVS( 1250 )
   W( 125 ) = JVS( 1251 )
   W( 126 ) = JVS( 1252 )
   W( 127 ) = JVS( 1253 )
   W( 128 ) = JVS( 1254 )
   W( 129 ) = JVS( 1255 )
   W( 130 ) = JVS( 1256 )
   W( 131 ) = JVS( 1257 )
   W( 132 ) = JVS( 1258 )
   W( 133 ) = JVS( 1259 )
   W( 134 ) = JVS( 1260 )
   W( 135 ) = JVS( 1261 )
   W( 136 ) = JVS( 1262 )
  a = -W( 1 ) / JVS(            1  )
  W( 1 ) = -a
  a = -W( 15 ) / JVS(           99  )
  W( 15 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 100 )
  a = -W( 48 ) / JVS(          216  )
  W( 48 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 217 )
  W( 136 ) = W( 136 ) + a*JVS( 218 )
  a = -W( 56 ) / JVS(          254  )
  W( 56 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 255 )
  W( 131 ) = W( 131 ) + a*JVS( 256 )
  W( 133 ) = W( 133 ) + a*JVS( 257 )
  a = -W( 58 ) / JVS(          265  )
  W( 58 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 266 )
  W( 133 ) = W( 133 ) + a*JVS( 267 )
  W( 136 ) = W( 136 ) + a*JVS( 268 )
  a = -W( 60 ) / JVS(          279  )
  W( 60 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 280 )
  W( 133 ) = W( 133 ) + a*JVS( 281 )
  W( 136 ) = W( 136 ) + a*JVS( 282 )
  a = -W( 61 ) / JVS(          286  )
  W( 61 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 287 )
  W( 133 ) = W( 133 ) + a*JVS( 288 )
  W( 136 ) = W( 136 ) + a*JVS( 289 )
  a = -W( 66 ) / JVS(          315  )
  W( 66 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 316 )
  W( 133 ) = W( 133 ) + a*JVS( 317 )
  W( 136 ) = W( 136 ) + a*JVS( 318 )
  a = -W( 69 ) / JVS(          332  )
  W( 69 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 333 )
  W( 133 ) = W( 133 ) + a*JVS( 334 )
  W( 136 ) = W( 136 ) + a*JVS( 335 )
  a = -W( 73 ) / JVS(          354  )
  W( 73 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 355 )
  W( 131 ) = W( 131 ) + a*JVS( 356 )
  W( 133 ) = W( 133 ) + a*JVS( 357 )
  W( 136 ) = W( 136 ) + a*JVS( 358 )
  a = -W( 74 ) / JVS(          362  )
  W( 74 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 363 )
  W( 133 ) = W( 133 ) + a*JVS( 364 )
  W( 134 ) = W( 134 ) + a*JVS( 365 )
  W( 136 ) = W( 136 ) + a*JVS( 366 )
  a = -W( 80 ) / JVS(          401  )
  W( 80 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 402 )
  W( 126 ) = W( 126 ) + a*JVS( 403 )
  W( 132 ) = W( 132 ) + a*JVS( 404 )
  W( 136 ) = W( 136 ) + a*JVS( 405 )
  a = -W( 81 ) / JVS(          409  )
  W( 81 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 410 )
  W( 132 ) = W( 132 ) + a*JVS( 411 )
  W( 133 ) = W( 133 ) + a*JVS( 412 )
  W( 134 ) = W( 134 ) + a*JVS( 413 )
  W( 136 ) = W( 136 ) + a*JVS( 414 )
  a = -W( 85 ) / JVS(          428  )
  W( 85 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 429 )
  W( 131 ) = W( 131 ) + a*JVS( 430 )
  W( 133 ) = W( 133 ) + a*JVS( 431 )
  W( 136 ) = W( 136 ) + a*JVS( 432 )
  a = -W( 86 ) / JVS(          434  )
  W( 86 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 435 )
  W( 131 ) = W( 131 ) + a*JVS( 436 )
  W( 133 ) = W( 133 ) + a*JVS( 437 )
  W( 136 ) = W( 136 ) + a*JVS( 438 )
  a = -W( 87 ) / JVS(          441  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 90 ) = W( 90 ) + a*JVS( 443 )
  W( 131 ) = W( 131 ) + a*JVS( 444 )
  W( 133 ) = W( 133 ) + a*JVS( 445 )
  W( 134 ) = W( 134 ) + a*JVS( 446 )
  W( 136 ) = W( 136 ) + a*JVS( 447 )
  a = -W( 88 ) / JVS(          450  )
  W( 88 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 451 )
  W( 133 ) = W( 133 ) + a*JVS( 452 )
  W( 136 ) = W( 136 ) + a*JVS( 453 )
  a = -W( 90 ) / JVS(          460  )
  W( 90 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 461 )
  W( 133 ) = W( 133 ) + a*JVS( 462 )
  W( 136 ) = W( 136 ) + a*JVS( 463 )
  a = -W( 91 ) / JVS(          467  )
  W( 91 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 468 )
  W( 133 ) = W( 133 ) + a*JVS( 469 )
  W( 134 ) = W( 134 ) + a*JVS( 470 )
  W( 136 ) = W( 136 ) + a*JVS( 471 )
  a = -W( 92 ) / JVS(          473  )
  W( 92 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 474 )
  W( 126 ) = W( 126 ) + a*JVS( 475 )
  W( 131 ) = W( 131 ) + a*JVS( 476 )
  W( 132 ) = W( 132 ) + a*JVS( 477 )
  W( 134 ) = W( 134 ) + a*JVS( 478 )
  W( 136 ) = W( 136 ) + a*JVS( 479 )
  a = -W( 93 ) / JVS(          482  )
  W( 93 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 483 )
  W( 131 ) = W( 131 ) + a*JVS( 484 )
  W( 133 ) = W( 133 ) + a*JVS( 485 )
  W( 136 ) = W( 136 ) + a*JVS( 486 )
  a = -W( 94 ) / JVS(          490  )
  W( 94 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 491 )
  W( 133 ) = W( 133 ) + a*JVS( 492 )
  W( 134 ) = W( 134 ) + a*JVS( 493 )
  W( 136 ) = W( 136 ) + a*JVS( 494 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 533 )
  W( 129 ) = W( 129 ) + a*JVS( 534 )
  W( 131 ) = W( 131 ) + a*JVS( 535 )
  W( 133 ) = W( 133 ) + a*JVS( 536 )
  W( 136 ) = W( 136 ) + a*JVS( 537 )
  a = -W( 99 ) / JVS(          550  )
  W( 99 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 133 ) = W( 133 ) + a*JVS( 553 )
  W( 135 ) = W( 135 ) + a*JVS( 554 )
  a = -W( 102 ) / JVS(          564  )
  W( 102 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 565 )
  W( 132 ) = W( 132 ) + a*JVS( 566 )
  W( 133 ) = W( 133 ) + a*JVS( 567 )
  W( 135 ) = W( 135 ) + a*JVS( 568 )
  W( 136 ) = W( 136 ) + a*JVS( 569 )
  a = -W( 103 ) / JVS(          573  )
  W( 103 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 574 )
  W( 129 ) = W( 129 ) + a*JVS( 575 )
  W( 131 ) = W( 131 ) + a*JVS( 576 )
  W( 132 ) = W( 132 ) + a*JVS( 577 )
  W( 133 ) = W( 133 ) + a*JVS( 578 )
  W( 135 ) = W( 135 ) + a*JVS( 579 )
  W( 136 ) = W( 136 ) + a*JVS( 580 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  a = -W( 109 ) / JVS(          682  )
  W( 109 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 683 )
  W( 132 ) = W( 132 ) + a*JVS( 684 )
  W( 135 ) = W( 135 ) + a*JVS( 685 )
  W( 136 ) = W( 136 ) + a*JVS( 686 )
  a = -W( 111 ) / JVS(          704  )
  W( 111 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 705 )
  W( 131 ) = W( 131 ) + a*JVS( 706 )
  W( 132 ) = W( 132 ) + a*JVS( 707 )
  W( 133 ) = W( 133 ) + a*JVS( 708 )
  W( 135 ) = W( 135 ) + a*JVS( 709 )
  W( 136 ) = W( 136 ) + a*JVS( 710 )
  a = -W( 113 ) / JVS(          738  )
  W( 113 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 739 )
  W( 132 ) = W( 132 ) + a*JVS( 740 )
  W( 133 ) = W( 133 ) + a*JVS( 741 )
  W( 135 ) = W( 135 ) + a*JVS( 742 )
  W( 136 ) = W( 136 ) + a*JVS( 743 )
  a = -W( 118 ) / JVS(          824  )
  W( 118 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 825 )
  W( 131 ) = W( 131 ) + a*JVS( 826 )
  W( 132 ) = W( 132 ) + a*JVS( 827 )
  W( 133 ) = W( 133 ) + a*JVS( 828 )
  W( 135 ) = W( 135 ) + a*JVS( 829 )
  W( 136 ) = W( 136 ) + a*JVS( 830 )
  a = -W( 119 ) / JVS(          836  )
  W( 119 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 837 )
  W( 131 ) = W( 131 ) + a*JVS( 838 )
  W( 132 ) = W( 132 ) + a*JVS( 839 )
  W( 133 ) = W( 133 ) + a*JVS( 840 )
  W( 135 ) = W( 135 ) + a*JVS( 841 )
  W( 136 ) = W( 136 ) + a*JVS( 842 )
  a = -W( 120 ) / JVS(          847  )
  W( 120 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 848 )
  W( 131 ) = W( 131 ) + a*JVS( 849 )
  W( 132 ) = W( 132 ) + a*JVS( 850 )
  W( 133 ) = W( 133 ) + a*JVS( 851 )
  W( 135 ) = W( 135 ) + a*JVS( 852 )
  W( 136 ) = W( 136 ) + a*JVS( 853 )
  a = -W( 121 ) / JVS(          867  )
  W( 121 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 868 )
  W( 131 ) = W( 131 ) + a*JVS( 869 )
  W( 132 ) = W( 132 ) + a*JVS( 870 )
  W( 133 ) = W( 133 ) + a*JVS( 871 )
  W( 135 ) = W( 135 ) + a*JVS( 872 )
  W( 136 ) = W( 136 ) + a*JVS( 873 )
  a = -W( 123 ) / JVS(          900  )
  W( 123 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 901 )
  W( 127 ) = W( 127 ) + a*JVS( 902 )
  W( 128 ) = W( 128 ) + a*JVS( 903 )
  W( 129 ) = W( 129 ) + a*JVS( 904 )
  W( 131 ) = W( 131 ) + a*JVS( 905 )
  W( 132 ) = W( 132 ) + a*JVS( 906 )
  W( 133 ) = W( 133 ) + a*JVS( 907 )
  W( 135 ) = W( 135 ) + a*JVS( 908 )
  W( 136 ) = W( 136 ) + a*JVS( 909 )
  a = -W( 124 ) / JVS(          913  )
  W( 124 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 914 )
  W( 128 ) = W( 128 ) + a*JVS( 915 )
  W( 129 ) = W( 129 ) + a*JVS( 916 )
  W( 131 ) = W( 131 ) + a*JVS( 917 )
  W( 132 ) = W( 132 ) + a*JVS( 918 )
  W( 133 ) = W( 133 ) + a*JVS( 919 )
  W( 135 ) = W( 135 ) + a*JVS( 920 )
  W( 136 ) = W( 136 ) + a*JVS( 921 )
  a = -W( 125 ) / JVS(          925  )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 926 )
  W( 128 ) = W( 128 ) + a*JVS( 927 )
  W( 129 ) = W( 129 ) + a*JVS( 928 )
  W( 131 ) = W( 131 ) + a*JVS( 929 )
  W( 132 ) = W( 132 ) + a*JVS( 930 )
  W( 133 ) = W( 133 ) + a*JVS( 931 )
  W( 135 ) = W( 135 ) + a*JVS( 932 )
  W( 136 ) = W( 136 ) + a*JVS( 933 )
  a = -W( 126 ) / JVS(          980  )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 981 )
  W( 128 ) = W( 128 ) + a*JVS( 982 )
  W( 129 ) = W( 129 ) + a*JVS( 983 )
  W( 130 ) = W( 130 ) + a*JVS( 984 )
  W( 131 ) = W( 131 ) + a*JVS( 985 )
  W( 132 ) = W( 132 ) + a*JVS( 986 )
  W( 133 ) = W( 133 ) + a*JVS( 987 )
  W( 134 ) = W( 134 ) + a*JVS( 988 )
  W( 135 ) = W( 135 ) + a*JVS( 989 )
  W( 136 ) = W( 136 ) + a*JVS( 990 )
  a = -W( 127 ) / JVS(          993  )
  W( 127 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 994 )
  W( 129 ) = W( 129 ) + a*JVS( 995 )
  W( 131 ) = W( 131 ) + a*JVS( 996 )
  W( 132 ) = W( 132 ) + a*JVS( 997 )
  W( 133 ) = W( 133 ) + a*JVS( 998 )
  W( 135 ) = W( 135 ) + a*JVS( 999 )
  W( 136 ) = W( 136 ) + a*JVS( 1000 )
  a = -W( 128 ) / JVS(         1034  )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 1035 )
  W( 130 ) = W( 130 ) + a*JVS( 1036 )
  W( 131 ) = W( 131 ) + a*JVS( 1037 )
  W( 132 ) = W( 132 ) + a*JVS( 1038 )
  W( 133 ) = W( 133 ) + a*JVS( 1039 )
  W( 134 ) = W( 134 ) + a*JVS( 1040 )
  W( 135 ) = W( 135 ) + a*JVS( 1041 )
  W( 136 ) = W( 136 ) + a*JVS( 1042 )
  a = -W( 129 ) / JVS(         1070  )
  W( 129 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 1071 )
  W( 131 ) = W( 131 ) + a*JVS( 1072 )
  W( 132 ) = W( 132 ) + a*JVS( 1073 )
  W( 133 ) = W( 133 ) + a*JVS( 1074 )
  W( 134 ) = W( 134 ) + a*JVS( 1075 )
  W( 135 ) = W( 135 ) + a*JVS( 1076 )
  W( 136 ) = W( 136 ) + a*JVS( 1077 )
  a = -W( 130 ) / JVS(         1086  )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 1087 )
  W( 132 ) = W( 132 ) + a*JVS( 1088 )
  W( 133 ) = W( 133 ) + a*JVS( 1089 )
  W( 134 ) = W( 134 ) + a*JVS( 1090 )
  W( 135 ) = W( 135 ) + a*JVS( 1091 )
  W( 136 ) = W( 136 ) + a*JVS( 1092 )
  a = -W( 131 ) / JVS(         1184  )
  W( 131 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 1185 )
  W( 133 ) = W( 133 ) + a*JVS( 1186 )
  W( 134 ) = W( 134 ) + a*JVS( 1187 )
  W( 135 ) = W( 135 ) + a*JVS( 1188 )
  W( 136 ) = W( 136 ) + a*JVS( 1189 )
  a = -W( 132 ) / JVS(         1210  )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1211 )
  W( 134 ) = W( 134 ) + a*JVS( 1212 )
  W( 135 ) = W( 135 ) + a*JVS( 1213 )
  W( 136 ) = W( 136 ) + a*JVS( 1214 )
  JVS( 1215) = W( 1 )
  JVS( 1216) = W( 15 )
  JVS( 1217) = W( 48 )
  JVS( 1218) = W( 56 )
  JVS( 1219) = W( 58 )
  JVS( 1220) = W( 60 )
  JVS( 1221) = W( 61 )
  JVS( 1222) = W( 66 )
  JVS( 1223) = W( 69 )
  JVS( 1224) = W( 73 )
  JVS( 1225) = W( 74 )
  JVS( 1226) = W( 80 )
  JVS( 1227) = W( 81 )
  JVS( 1228) = W( 85 )
  JVS( 1229) = W( 86 )
  JVS( 1230) = W( 87 )
  JVS( 1231) = W( 88 )
  JVS( 1232) = W( 90 )
  JVS( 1233) = W( 91 )
  JVS( 1234) = W( 92 )
  JVS( 1235) = W( 93 )
  JVS( 1236) = W( 94 )
  JVS( 1237) = W( 97 )
  JVS( 1238) = W( 99 )
  JVS( 1239) = W( 102 )
  JVS( 1240) = W( 103 )
  JVS( 1241) = W( 106 )
  JVS( 1242) = W( 109 )
  JVS( 1243) = W( 111 )
  JVS( 1244) = W( 113 )
  JVS( 1245) = W( 118 )
  JVS( 1246) = W( 119 )
  JVS( 1247) = W( 120 )
  JVS( 1248) = W( 121 )
  JVS( 1249) = W( 123 )
  JVS( 1250) = W( 124 )
  JVS( 1251) = W( 125 )
  JVS( 1252) = W( 126 )
  JVS( 1253) = W( 127 )
  JVS( 1254) = W( 128 )
  JVS( 1255) = W( 129 )
  JVS( 1256) = W( 130 )
  JVS( 1257) = W( 131 )
  JVS( 1258) = W( 132 )
  JVS( 1259) = W( 133 )
  JVS( 1260) = W( 134 )
  JVS( 1261) = W( 135 )
  JVS( 1262) = W( 136 )
  IF ( ABS(  JVS( 1317 )) < TINY(a) ) THEN
         IER = 134                                     
         RETURN
  END IF
   W( 17 ) = JVS( 1263 )
   W( 36 ) = JVS( 1264 )
   W( 48 ) = JVS( 1265 )
   W( 52 ) = JVS( 1266 )
   W( 56 ) = JVS( 1267 )
   W( 58 ) = JVS( 1268 )
   W( 59 ) = JVS( 1269 )
   W( 60 ) = JVS( 1270 )
   W( 61 ) = JVS( 1271 )
   W( 66 ) = JVS( 1272 )
   W( 67 ) = JVS( 1273 )
   W( 68 ) = JVS( 1274 )
   W( 69 ) = JVS( 1275 )
   W( 71 ) = JVS( 1276 )
   W( 73 ) = JVS( 1277 )
   W( 74 ) = JVS( 1278 )
   W( 77 ) = JVS( 1279 )
   W( 81 ) = JVS( 1280 )
   W( 85 ) = JVS( 1281 )
   W( 86 ) = JVS( 1282 )
   W( 87 ) = JVS( 1283 )
   W( 88 ) = JVS( 1284 )
   W( 90 ) = JVS( 1285 )
   W( 91 ) = JVS( 1286 )
   W( 92 ) = JVS( 1287 )
   W( 93 ) = JVS( 1288 )
   W( 94 ) = JVS( 1289 )
   W( 97 ) = JVS( 1290 )
   W( 99 ) = JVS( 1291 )
   W( 102 ) = JVS( 1292 )
   W( 103 ) = JVS( 1293 )
   W( 106 ) = JVS( 1294 )
   W( 109 ) = JVS( 1295 )
   W( 111 ) = JVS( 1296 )
   W( 112 ) = JVS( 1297 )
   W( 113 ) = JVS( 1298 )
   W( 116 ) = JVS( 1299 )
   W( 117 ) = JVS( 1300 )
   W( 118 ) = JVS( 1301 )
   W( 119 ) = JVS( 1302 )
   W( 120 ) = JVS( 1303 )
   W( 121 ) = JVS( 1304 )
   W( 122 ) = JVS( 1305 )
   W( 123 ) = JVS( 1306 )
   W( 124 ) = JVS( 1307 )
   W( 125 ) = JVS( 1308 )
   W( 126 ) = JVS( 1309 )
   W( 127 ) = JVS( 1310 )
   W( 128 ) = JVS( 1311 )
   W( 129 ) = JVS( 1312 )
   W( 130 ) = JVS( 1313 )
   W( 131 ) = JVS( 1314 )
   W( 132 ) = JVS( 1315 )
   W( 133 ) = JVS( 1316 )
   W( 134 ) = JVS( 1317 )
   W( 135 ) = JVS( 1318 )
   W( 136 ) = JVS( 1319 )
  a = -W( 17 ) / JVS(          105  )
  W( 17 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 106 )
  W( 134 ) = W( 134 ) + a*JVS( 107 )
  a = -W( 36 ) / JVS(          169  )
  W( 36 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 170 )
  W( 135 ) = W( 135 ) + a*JVS( 171 )
  a = -W( 48 ) / JVS(          216  )
  W( 48 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 217 )
  W( 136 ) = W( 136 ) + a*JVS( 218 )
  a = -W( 52 ) / JVS(          231  )
  W( 52 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 232 )
  W( 134 ) = W( 134 ) + a*JVS( 233 )
  W( 136 ) = W( 136 ) + a*JVS( 234 )
  a = -W( 56 ) / JVS(          254  )
  W( 56 ) = -a
  W( 126 ) = W( 126 ) + a*JVS( 255 )
  W( 131 ) = W( 131 ) + a*JVS( 256 )
  W( 133 ) = W( 133 ) + a*JVS( 257 )
  a = -W( 58 ) / JVS(          265  )
  W( 58 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 266 )
  W( 133 ) = W( 133 ) + a*JVS( 267 )
  W( 136 ) = W( 136 ) + a*JVS( 268 )
  a = -W( 59 ) / JVS(          272  )
  W( 59 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 273 )
  W( 132 ) = W( 132 ) + a*JVS( 274 )
  W( 133 ) = W( 133 ) + a*JVS( 275 )
  W( 134 ) = W( 134 ) + a*JVS( 276 )
  W( 136 ) = W( 136 ) + a*JVS( 277 )
  a = -W( 60 ) / JVS(          279  )
  W( 60 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 280 )
  W( 133 ) = W( 133 ) + a*JVS( 281 )
  W( 136 ) = W( 136 ) + a*JVS( 282 )
  a = -W( 61 ) / JVS(          286  )
  W( 61 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 287 )
  W( 133 ) = W( 133 ) + a*JVS( 288 )
  W( 136 ) = W( 136 ) + a*JVS( 289 )
  a = -W( 66 ) / JVS(          315  )
  W( 66 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 316 )
  W( 133 ) = W( 133 ) + a*JVS( 317 )
  W( 136 ) = W( 136 ) + a*JVS( 318 )
  a = -W( 67 ) / JVS(          319  )
  W( 67 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 320 )
  W( 109 ) = W( 109 ) + a*JVS( 321 )
  W( 113 ) = W( 113 ) + a*JVS( 322 )
  W( 133 ) = W( 133 ) + a*JVS( 323 )
  W( 135 ) = W( 135 ) + a*JVS( 324 )
  W( 136 ) = W( 136 ) + a*JVS( 325 )
  a = -W( 68 ) / JVS(          326  )
  W( 68 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 327 )
  W( 134 ) = W( 134 ) + a*JVS( 328 )
  W( 136 ) = W( 136 ) + a*JVS( 329 )
  a = -W( 69 ) / JVS(          332  )
  W( 69 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 333 )
  W( 133 ) = W( 133 ) + a*JVS( 334 )
  W( 136 ) = W( 136 ) + a*JVS( 335 )
  a = -W( 71 ) / JVS(          342  )
  W( 71 ) = -a
  W( 112 ) = W( 112 ) + a*JVS( 343 )
  W( 117 ) = W( 117 ) + a*JVS( 344 )
  W( 126 ) = W( 126 ) + a*JVS( 345 )
  W( 134 ) = W( 134 ) + a*JVS( 346 )
  W( 135 ) = W( 135 ) + a*JVS( 347 )
  W( 136 ) = W( 136 ) + a*JVS( 348 )
  a = -W( 73 ) / JVS(          354  )
  W( 73 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 355 )
  W( 131 ) = W( 131 ) + a*JVS( 356 )
  W( 133 ) = W( 133 ) + a*JVS( 357 )
  W( 136 ) = W( 136 ) + a*JVS( 358 )
  a = -W( 74 ) / JVS(          362  )
  W( 74 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 363 )
  W( 133 ) = W( 133 ) + a*JVS( 364 )
  W( 134 ) = W( 134 ) + a*JVS( 365 )
  W( 136 ) = W( 136 ) + a*JVS( 366 )
  a = -W( 77 ) / JVS(          375  )
  W( 77 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 376 )
  W( 134 ) = W( 134 ) + a*JVS( 377 )
  W( 136 ) = W( 136 ) + a*JVS( 378 )
  a = -W( 81 ) / JVS(          409  )
  W( 81 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 410 )
  W( 132 ) = W( 132 ) + a*JVS( 411 )
  W( 133 ) = W( 133 ) + a*JVS( 412 )
  W( 134 ) = W( 134 ) + a*JVS( 413 )
  W( 136 ) = W( 136 ) + a*JVS( 414 )
  a = -W( 85 ) / JVS(          428  )
  W( 85 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 429 )
  W( 131 ) = W( 131 ) + a*JVS( 430 )
  W( 133 ) = W( 133 ) + a*JVS( 431 )
  W( 136 ) = W( 136 ) + a*JVS( 432 )
  a = -W( 86 ) / JVS(          434  )
  W( 86 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 435 )
  W( 131 ) = W( 131 ) + a*JVS( 436 )
  W( 133 ) = W( 133 ) + a*JVS( 437 )
  W( 136 ) = W( 136 ) + a*JVS( 438 )
  a = -W( 87 ) / JVS(          441  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 442 )
  W( 90 ) = W( 90 ) + a*JVS( 443 )
  W( 131 ) = W( 131 ) + a*JVS( 444 )
  W( 133 ) = W( 133 ) + a*JVS( 445 )
  W( 134 ) = W( 134 ) + a*JVS( 446 )
  W( 136 ) = W( 136 ) + a*JVS( 447 )
  a = -W( 88 ) / JVS(          450  )
  W( 88 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 451 )
  W( 133 ) = W( 133 ) + a*JVS( 452 )
  W( 136 ) = W( 136 ) + a*JVS( 453 )
  a = -W( 90 ) / JVS(          460  )
  W( 90 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 461 )
  W( 133 ) = W( 133 ) + a*JVS( 462 )
  W( 136 ) = W( 136 ) + a*JVS( 463 )
  a = -W( 91 ) / JVS(          467  )
  W( 91 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 468 )
  W( 133 ) = W( 133 ) + a*JVS( 469 )
  W( 134 ) = W( 134 ) + a*JVS( 470 )
  W( 136 ) = W( 136 ) + a*JVS( 471 )
  a = -W( 92 ) / JVS(          473  )
  W( 92 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 474 )
  W( 126 ) = W( 126 ) + a*JVS( 475 )
  W( 131 ) = W( 131 ) + a*JVS( 476 )
  W( 132 ) = W( 132 ) + a*JVS( 477 )
  W( 134 ) = W( 134 ) + a*JVS( 478 )
  W( 136 ) = W( 136 ) + a*JVS( 479 )
  a = -W( 93 ) / JVS(          482  )
  W( 93 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 483 )
  W( 131 ) = W( 131 ) + a*JVS( 484 )
  W( 133 ) = W( 133 ) + a*JVS( 485 )
  W( 136 ) = W( 136 ) + a*JVS( 486 )
  a = -W( 94 ) / JVS(          490  )
  W( 94 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 491 )
  W( 133 ) = W( 133 ) + a*JVS( 492 )
  W( 134 ) = W( 134 ) + a*JVS( 493 )
  W( 136 ) = W( 136 ) + a*JVS( 494 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 533 )
  W( 129 ) = W( 129 ) + a*JVS( 534 )
  W( 131 ) = W( 131 ) + a*JVS( 535 )
  W( 133 ) = W( 133 ) + a*JVS( 536 )
  W( 136 ) = W( 136 ) + a*JVS( 537 )
  a = -W( 99 ) / JVS(          550  )
  W( 99 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 133 ) = W( 133 ) + a*JVS( 553 )
  W( 135 ) = W( 135 ) + a*JVS( 554 )
  a = -W( 102 ) / JVS(          564  )
  W( 102 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 565 )
  W( 132 ) = W( 132 ) + a*JVS( 566 )
  W( 133 ) = W( 133 ) + a*JVS( 567 )
  W( 135 ) = W( 135 ) + a*JVS( 568 )
  W( 136 ) = W( 136 ) + a*JVS( 569 )
  a = -W( 103 ) / JVS(          573  )
  W( 103 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 574 )
  W( 129 ) = W( 129 ) + a*JVS( 575 )
  W( 131 ) = W( 131 ) + a*JVS( 576 )
  W( 132 ) = W( 132 ) + a*JVS( 577 )
  W( 133 ) = W( 133 ) + a*JVS( 578 )
  W( 135 ) = W( 135 ) + a*JVS( 579 )
  W( 136 ) = W( 136 ) + a*JVS( 580 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  a = -W( 109 ) / JVS(          682  )
  W( 109 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 683 )
  W( 132 ) = W( 132 ) + a*JVS( 684 )
  W( 135 ) = W( 135 ) + a*JVS( 685 )
  W( 136 ) = W( 136 ) + a*JVS( 686 )
  a = -W( 111 ) / JVS(          704  )
  W( 111 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 705 )
  W( 131 ) = W( 131 ) + a*JVS( 706 )
  W( 132 ) = W( 132 ) + a*JVS( 707 )
  W( 133 ) = W( 133 ) + a*JVS( 708 )
  W( 135 ) = W( 135 ) + a*JVS( 709 )
  W( 136 ) = W( 136 ) + a*JVS( 710 )
  a = -W( 112 ) / JVS(          722  )
  W( 112 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 723 )
  W( 124 ) = W( 124 ) + a*JVS( 724 )
  W( 129 ) = W( 129 ) + a*JVS( 725 )
  W( 131 ) = W( 131 ) + a*JVS( 726 )
  W( 132 ) = W( 132 ) + a*JVS( 727 )
  W( 133 ) = W( 133 ) + a*JVS( 728 )
  W( 135 ) = W( 135 ) + a*JVS( 729 )
  W( 136 ) = W( 136 ) + a*JVS( 730 )
  a = -W( 113 ) / JVS(          738  )
  W( 113 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 739 )
  W( 132 ) = W( 132 ) + a*JVS( 740 )
  W( 133 ) = W( 133 ) + a*JVS( 741 )
  W( 135 ) = W( 135 ) + a*JVS( 742 )
  W( 136 ) = W( 136 ) + a*JVS( 743 )
  a = -W( 116 ) / JVS(          775  )
  W( 116 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 776 )
  W( 119 ) = W( 119 ) + a*JVS( 777 )
  W( 120 ) = W( 120 ) + a*JVS( 778 )
  W( 125 ) = W( 125 ) + a*JVS( 779 )
  W( 127 ) = W( 127 ) + a*JVS( 780 )
  W( 129 ) = W( 129 ) + a*JVS( 781 )
  W( 131 ) = W( 131 ) + a*JVS( 782 )
  W( 132 ) = W( 132 ) + a*JVS( 783 )
  W( 133 ) = W( 133 ) + a*JVS( 784 )
  W( 135 ) = W( 135 ) + a*JVS( 785 )
  W( 136 ) = W( 136 ) + a*JVS( 786 )
  a = -W( 117 ) / JVS(          803  )
  W( 117 ) = -a
  W( 121 ) = W( 121 ) + a*JVS( 804 )
  W( 122 ) = W( 122 ) + a*JVS( 805 )
  W( 123 ) = W( 123 ) + a*JVS( 806 )
  W( 124 ) = W( 124 ) + a*JVS( 807 )
  W( 125 ) = W( 125 ) + a*JVS( 808 )
  W( 127 ) = W( 127 ) + a*JVS( 809 )
  W( 128 ) = W( 128 ) + a*JVS( 810 )
  W( 129 ) = W( 129 ) + a*JVS( 811 )
  W( 131 ) = W( 131 ) + a*JVS( 812 )
  W( 132 ) = W( 132 ) + a*JVS( 813 )
  W( 133 ) = W( 133 ) + a*JVS( 814 )
  W( 134 ) = W( 134 ) + a*JVS( 815 )
  W( 135 ) = W( 135 ) + a*JVS( 816 )
  W( 136 ) = W( 136 ) + a*JVS( 817 )
  a = -W( 118 ) / JVS(          824  )
  W( 118 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 825 )
  W( 131 ) = W( 131 ) + a*JVS( 826 )
  W( 132 ) = W( 132 ) + a*JVS( 827 )
  W( 133 ) = W( 133 ) + a*JVS( 828 )
  W( 135 ) = W( 135 ) + a*JVS( 829 )
  W( 136 ) = W( 136 ) + a*JVS( 830 )
  a = -W( 119 ) / JVS(          836  )
  W( 119 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 837 )
  W( 131 ) = W( 131 ) + a*JVS( 838 )
  W( 132 ) = W( 132 ) + a*JVS( 839 )
  W( 133 ) = W( 133 ) + a*JVS( 840 )
  W( 135 ) = W( 135 ) + a*JVS( 841 )
  W( 136 ) = W( 136 ) + a*JVS( 842 )
  a = -W( 120 ) / JVS(          847  )
  W( 120 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 848 )
  W( 131 ) = W( 131 ) + a*JVS( 849 )
  W( 132 ) = W( 132 ) + a*JVS( 850 )
  W( 133 ) = W( 133 ) + a*JVS( 851 )
  W( 135 ) = W( 135 ) + a*JVS( 852 )
  W( 136 ) = W( 136 ) + a*JVS( 853 )
  a = -W( 121 ) / JVS(          867  )
  W( 121 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 868 )
  W( 131 ) = W( 131 ) + a*JVS( 869 )
  W( 132 ) = W( 132 ) + a*JVS( 870 )
  W( 133 ) = W( 133 ) + a*JVS( 871 )
  W( 135 ) = W( 135 ) + a*JVS( 872 )
  W( 136 ) = W( 136 ) + a*JVS( 873 )
  a = -W( 122 ) / JVS(          879  )
  W( 122 ) = -a
  W( 123 ) = W( 123 ) + a*JVS( 880 )
  W( 124 ) = W( 124 ) + a*JVS( 881 )
  W( 125 ) = W( 125 ) + a*JVS( 882 )
  W( 127 ) = W( 127 ) + a*JVS( 883 )
  W( 128 ) = W( 128 ) + a*JVS( 884 )
  W( 129 ) = W( 129 ) + a*JVS( 885 )
  W( 130 ) = W( 130 ) + a*JVS( 886 )
  W( 131 ) = W( 131 ) + a*JVS( 887 )
  W( 132 ) = W( 132 ) + a*JVS( 888 )
  W( 133 ) = W( 133 ) + a*JVS( 889 )
  W( 134 ) = W( 134 ) + a*JVS( 890 )
  W( 135 ) = W( 135 ) + a*JVS( 891 )
  W( 136 ) = W( 136 ) + a*JVS( 892 )
  a = -W( 123 ) / JVS(          900  )
  W( 123 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 901 )
  W( 127 ) = W( 127 ) + a*JVS( 902 )
  W( 128 ) = W( 128 ) + a*JVS( 903 )
  W( 129 ) = W( 129 ) + a*JVS( 904 )
  W( 131 ) = W( 131 ) + a*JVS( 905 )
  W( 132 ) = W( 132 ) + a*JVS( 906 )
  W( 133 ) = W( 133 ) + a*JVS( 907 )
  W( 135 ) = W( 135 ) + a*JVS( 908 )
  W( 136 ) = W( 136 ) + a*JVS( 909 )
  a = -W( 124 ) / JVS(          913  )
  W( 124 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 914 )
  W( 128 ) = W( 128 ) + a*JVS( 915 )
  W( 129 ) = W( 129 ) + a*JVS( 916 )
  W( 131 ) = W( 131 ) + a*JVS( 917 )
  W( 132 ) = W( 132 ) + a*JVS( 918 )
  W( 133 ) = W( 133 ) + a*JVS( 919 )
  W( 135 ) = W( 135 ) + a*JVS( 920 )
  W( 136 ) = W( 136 ) + a*JVS( 921 )
  a = -W( 125 ) / JVS(          925  )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 926 )
  W( 128 ) = W( 128 ) + a*JVS( 927 )
  W( 129 ) = W( 129 ) + a*JVS( 928 )
  W( 131 ) = W( 131 ) + a*JVS( 929 )
  W( 132 ) = W( 132 ) + a*JVS( 930 )
  W( 133 ) = W( 133 ) + a*JVS( 931 )
  W( 135 ) = W( 135 ) + a*JVS( 932 )
  W( 136 ) = W( 136 ) + a*JVS( 933 )
  a = -W( 126 ) / JVS(          980  )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 981 )
  W( 128 ) = W( 128 ) + a*JVS( 982 )
  W( 129 ) = W( 129 ) + a*JVS( 983 )
  W( 130 ) = W( 130 ) + a*JVS( 984 )
  W( 131 ) = W( 131 ) + a*JVS( 985 )
  W( 132 ) = W( 132 ) + a*JVS( 986 )
  W( 133 ) = W( 133 ) + a*JVS( 987 )
  W( 134 ) = W( 134 ) + a*JVS( 988 )
  W( 135 ) = W( 135 ) + a*JVS( 989 )
  W( 136 ) = W( 136 ) + a*JVS( 990 )
  a = -W( 127 ) / JVS(          993  )
  W( 127 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 994 )
  W( 129 ) = W( 129 ) + a*JVS( 995 )
  W( 131 ) = W( 131 ) + a*JVS( 996 )
  W( 132 ) = W( 132 ) + a*JVS( 997 )
  W( 133 ) = W( 133 ) + a*JVS( 998 )
  W( 135 ) = W( 135 ) + a*JVS( 999 )
  W( 136 ) = W( 136 ) + a*JVS( 1000 )
  a = -W( 128 ) / JVS(         1034  )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 1035 )
  W( 130 ) = W( 130 ) + a*JVS( 1036 )
  W( 131 ) = W( 131 ) + a*JVS( 1037 )
  W( 132 ) = W( 132 ) + a*JVS( 1038 )
  W( 133 ) = W( 133 ) + a*JVS( 1039 )
  W( 134 ) = W( 134 ) + a*JVS( 1040 )
  W( 135 ) = W( 135 ) + a*JVS( 1041 )
  W( 136 ) = W( 136 ) + a*JVS( 1042 )
  a = -W( 129 ) / JVS(         1070  )
  W( 129 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 1071 )
  W( 131 ) = W( 131 ) + a*JVS( 1072 )
  W( 132 ) = W( 132 ) + a*JVS( 1073 )
  W( 133 ) = W( 133 ) + a*JVS( 1074 )
  W( 134 ) = W( 134 ) + a*JVS( 1075 )
  W( 135 ) = W( 135 ) + a*JVS( 1076 )
  W( 136 ) = W( 136 ) + a*JVS( 1077 )
  a = -W( 130 ) / JVS(         1086  )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 1087 )
  W( 132 ) = W( 132 ) + a*JVS( 1088 )
  W( 133 ) = W( 133 ) + a*JVS( 1089 )
  W( 134 ) = W( 134 ) + a*JVS( 1090 )
  W( 135 ) = W( 135 ) + a*JVS( 1091 )
  W( 136 ) = W( 136 ) + a*JVS( 1092 )
  a = -W( 131 ) / JVS(         1184  )
  W( 131 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 1185 )
  W( 133 ) = W( 133 ) + a*JVS( 1186 )
  W( 134 ) = W( 134 ) + a*JVS( 1187 )
  W( 135 ) = W( 135 ) + a*JVS( 1188 )
  W( 136 ) = W( 136 ) + a*JVS( 1189 )
  a = -W( 132 ) / JVS(         1210  )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1211 )
  W( 134 ) = W( 134 ) + a*JVS( 1212 )
  W( 135 ) = W( 135 ) + a*JVS( 1213 )
  W( 136 ) = W( 136 ) + a*JVS( 1214 )
  a = -W( 133 ) / JVS(         1259  )
  W( 133 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 1260 )
  W( 135 ) = W( 135 ) + a*JVS( 1261 )
  W( 136 ) = W( 136 ) + a*JVS( 1262 )
  JVS( 1263) = W( 17 )
  JVS( 1264) = W( 36 )
  JVS( 1265) = W( 48 )
  JVS( 1266) = W( 52 )
  JVS( 1267) = W( 56 )
  JVS( 1268) = W( 58 )
  JVS( 1269) = W( 59 )
  JVS( 1270) = W( 60 )
  JVS( 1271) = W( 61 )
  JVS( 1272) = W( 66 )
  JVS( 1273) = W( 67 )
  JVS( 1274) = W( 68 )
  JVS( 1275) = W( 69 )
  JVS( 1276) = W( 71 )
  JVS( 1277) = W( 73 )
  JVS( 1278) = W( 74 )
  JVS( 1279) = W( 77 )
  JVS( 1280) = W( 81 )
  JVS( 1281) = W( 85 )
  JVS( 1282) = W( 86 )
  JVS( 1283) = W( 87 )
  JVS( 1284) = W( 88 )
  JVS( 1285) = W( 90 )
  JVS( 1286) = W( 91 )
  JVS( 1287) = W( 92 )
  JVS( 1288) = W( 93 )
  JVS( 1289) = W( 94 )
  JVS( 1290) = W( 97 )
  JVS( 1291) = W( 99 )
  JVS( 1292) = W( 102 )
  JVS( 1293) = W( 103 )
  JVS( 1294) = W( 106 )
  JVS( 1295) = W( 109 )
  JVS( 1296) = W( 111 )
  JVS( 1297) = W( 112 )
  JVS( 1298) = W( 113 )
  JVS( 1299) = W( 116 )
  JVS( 1300) = W( 117 )
  JVS( 1301) = W( 118 )
  JVS( 1302) = W( 119 )
  JVS( 1303) = W( 120 )
  JVS( 1304) = W( 121 )
  JVS( 1305) = W( 122 )
  JVS( 1306) = W( 123 )
  JVS( 1307) = W( 124 )
  JVS( 1308) = W( 125 )
  JVS( 1309) = W( 126 )
  JVS( 1310) = W( 127 )
  JVS( 1311) = W( 128 )
  JVS( 1312) = W( 129 )
  JVS( 1313) = W( 130 )
  JVS( 1314) = W( 131 )
  JVS( 1315) = W( 132 )
  JVS( 1316) = W( 133 )
  JVS( 1317) = W( 134 )
  JVS( 1318) = W( 135 )
  JVS( 1319) = W( 136 )
  IF ( ABS(  JVS( 1356 )) < TINY(a) ) THEN
         IER = 135                                     
         RETURN
  END IF
   W( 36 ) = JVS( 1320 )
   W( 39 ) = JVS( 1321 )
   W( 52 ) = JVS( 1322 )
   W( 68 ) = JVS( 1323 )
   W( 71 ) = JVS( 1324 )
   W( 77 ) = JVS( 1325 )
   W( 89 ) = JVS( 1326 )
   W( 99 ) = JVS( 1327 )
   W( 100 ) = JVS( 1328 )
   W( 101 ) = JVS( 1329 )
   W( 102 ) = JVS( 1330 )
   W( 104 ) = JVS( 1331 )
   W( 106 ) = JVS( 1332 )
   W( 109 ) = JVS( 1333 )
   W( 110 ) = JVS( 1334 )
   W( 112 ) = JVS( 1335 )
   W( 113 ) = JVS( 1336 )
   W( 116 ) = JVS( 1337 )
   W( 117 ) = JVS( 1338 )
   W( 118 ) = JVS( 1339 )
   W( 119 ) = JVS( 1340 )
   W( 120 ) = JVS( 1341 )
   W( 121 ) = JVS( 1342 )
   W( 122 ) = JVS( 1343 )
   W( 123 ) = JVS( 1344 )
   W( 124 ) = JVS( 1345 )
   W( 125 ) = JVS( 1346 )
   W( 126 ) = JVS( 1347 )
   W( 127 ) = JVS( 1348 )
   W( 128 ) = JVS( 1349 )
   W( 129 ) = JVS( 1350 )
   W( 130 ) = JVS( 1351 )
   W( 131 ) = JVS( 1352 )
   W( 132 ) = JVS( 1353 )
   W( 133 ) = JVS( 1354 )
   W( 134 ) = JVS( 1355 )
   W( 135 ) = JVS( 1356 )
   W( 136 ) = JVS( 1357 )
  a = -W( 36 ) / JVS(          169  )
  W( 36 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 170 )
  W( 135 ) = W( 135 ) + a*JVS( 171 )
  a = -W( 39 ) / JVS(          177  )
  W( 39 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 178 )
  W( 136 ) = W( 136 ) + a*JVS( 179 )
  a = -W( 52 ) / JVS(          231  )
  W( 52 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 232 )
  W( 134 ) = W( 134 ) + a*JVS( 233 )
  W( 136 ) = W( 136 ) + a*JVS( 234 )
  a = -W( 68 ) / JVS(          326  )
  W( 68 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 327 )
  W( 134 ) = W( 134 ) + a*JVS( 328 )
  W( 136 ) = W( 136 ) + a*JVS( 329 )
  a = -W( 71 ) / JVS(          342  )
  W( 71 ) = -a
  W( 112 ) = W( 112 ) + a*JVS( 343 )
  W( 117 ) = W( 117 ) + a*JVS( 344 )
  W( 126 ) = W( 126 ) + a*JVS( 345 )
  W( 134 ) = W( 134 ) + a*JVS( 346 )
  W( 135 ) = W( 135 ) + a*JVS( 347 )
  W( 136 ) = W( 136 ) + a*JVS( 348 )
  a = -W( 77 ) / JVS(          375  )
  W( 77 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 376 )
  W( 134 ) = W( 134 ) + a*JVS( 377 )
  W( 136 ) = W( 136 ) + a*JVS( 378 )
  a = -W( 89 ) / JVS(          454  )
  W( 89 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 455 )
  W( 135 ) = W( 135 ) + a*JVS( 456 )
  W( 136 ) = W( 136 ) + a*JVS( 457 )
  a = -W( 99 ) / JVS(          550  )
  W( 99 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 133 ) = W( 133 ) + a*JVS( 553 )
  W( 135 ) = W( 135 ) + a*JVS( 554 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 102 ) / JVS(          564  )
  W( 102 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 565 )
  W( 132 ) = W( 132 ) + a*JVS( 566 )
  W( 133 ) = W( 133 ) + a*JVS( 567 )
  W( 135 ) = W( 135 ) + a*JVS( 568 )
  W( 136 ) = W( 136 ) + a*JVS( 569 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  a = -W( 109 ) / JVS(          682  )
  W( 109 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 683 )
  W( 132 ) = W( 132 ) + a*JVS( 684 )
  W( 135 ) = W( 135 ) + a*JVS( 685 )
  W( 136 ) = W( 136 ) + a*JVS( 686 )
  a = -W( 110 ) / JVS(          693  )
  W( 110 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 694 )
  W( 119 ) = W( 119 ) + a*JVS( 695 )
  W( 129 ) = W( 129 ) + a*JVS( 696 )
  W( 131 ) = W( 131 ) + a*JVS( 697 )
  W( 132 ) = W( 132 ) + a*JVS( 698 )
  W( 133 ) = W( 133 ) + a*JVS( 699 )
  W( 135 ) = W( 135 ) + a*JVS( 700 )
  W( 136 ) = W( 136 ) + a*JVS( 701 )
  a = -W( 112 ) / JVS(          722  )
  W( 112 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 723 )
  W( 124 ) = W( 124 ) + a*JVS( 724 )
  W( 129 ) = W( 129 ) + a*JVS( 725 )
  W( 131 ) = W( 131 ) + a*JVS( 726 )
  W( 132 ) = W( 132 ) + a*JVS( 727 )
  W( 133 ) = W( 133 ) + a*JVS( 728 )
  W( 135 ) = W( 135 ) + a*JVS( 729 )
  W( 136 ) = W( 136 ) + a*JVS( 730 )
  a = -W( 113 ) / JVS(          738  )
  W( 113 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 739 )
  W( 132 ) = W( 132 ) + a*JVS( 740 )
  W( 133 ) = W( 133 ) + a*JVS( 741 )
  W( 135 ) = W( 135 ) + a*JVS( 742 )
  W( 136 ) = W( 136 ) + a*JVS( 743 )
  a = -W( 116 ) / JVS(          775  )
  W( 116 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 776 )
  W( 119 ) = W( 119 ) + a*JVS( 777 )
  W( 120 ) = W( 120 ) + a*JVS( 778 )
  W( 125 ) = W( 125 ) + a*JVS( 779 )
  W( 127 ) = W( 127 ) + a*JVS( 780 )
  W( 129 ) = W( 129 ) + a*JVS( 781 )
  W( 131 ) = W( 131 ) + a*JVS( 782 )
  W( 132 ) = W( 132 ) + a*JVS( 783 )
  W( 133 ) = W( 133 ) + a*JVS( 784 )
  W( 135 ) = W( 135 ) + a*JVS( 785 )
  W( 136 ) = W( 136 ) + a*JVS( 786 )
  a = -W( 117 ) / JVS(          803  )
  W( 117 ) = -a
  W( 121 ) = W( 121 ) + a*JVS( 804 )
  W( 122 ) = W( 122 ) + a*JVS( 805 )
  W( 123 ) = W( 123 ) + a*JVS( 806 )
  W( 124 ) = W( 124 ) + a*JVS( 807 )
  W( 125 ) = W( 125 ) + a*JVS( 808 )
  W( 127 ) = W( 127 ) + a*JVS( 809 )
  W( 128 ) = W( 128 ) + a*JVS( 810 )
  W( 129 ) = W( 129 ) + a*JVS( 811 )
  W( 131 ) = W( 131 ) + a*JVS( 812 )
  W( 132 ) = W( 132 ) + a*JVS( 813 )
  W( 133 ) = W( 133 ) + a*JVS( 814 )
  W( 134 ) = W( 134 ) + a*JVS( 815 )
  W( 135 ) = W( 135 ) + a*JVS( 816 )
  W( 136 ) = W( 136 ) + a*JVS( 817 )
  a = -W( 118 ) / JVS(          824  )
  W( 118 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 825 )
  W( 131 ) = W( 131 ) + a*JVS( 826 )
  W( 132 ) = W( 132 ) + a*JVS( 827 )
  W( 133 ) = W( 133 ) + a*JVS( 828 )
  W( 135 ) = W( 135 ) + a*JVS( 829 )
  W( 136 ) = W( 136 ) + a*JVS( 830 )
  a = -W( 119 ) / JVS(          836  )
  W( 119 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 837 )
  W( 131 ) = W( 131 ) + a*JVS( 838 )
  W( 132 ) = W( 132 ) + a*JVS( 839 )
  W( 133 ) = W( 133 ) + a*JVS( 840 )
  W( 135 ) = W( 135 ) + a*JVS( 841 )
  W( 136 ) = W( 136 ) + a*JVS( 842 )
  a = -W( 120 ) / JVS(          847  )
  W( 120 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 848 )
  W( 131 ) = W( 131 ) + a*JVS( 849 )
  W( 132 ) = W( 132 ) + a*JVS( 850 )
  W( 133 ) = W( 133 ) + a*JVS( 851 )
  W( 135 ) = W( 135 ) + a*JVS( 852 )
  W( 136 ) = W( 136 ) + a*JVS( 853 )
  a = -W( 121 ) / JVS(          867  )
  W( 121 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 868 )
  W( 131 ) = W( 131 ) + a*JVS( 869 )
  W( 132 ) = W( 132 ) + a*JVS( 870 )
  W( 133 ) = W( 133 ) + a*JVS( 871 )
  W( 135 ) = W( 135 ) + a*JVS( 872 )
  W( 136 ) = W( 136 ) + a*JVS( 873 )
  a = -W( 122 ) / JVS(          879  )
  W( 122 ) = -a
  W( 123 ) = W( 123 ) + a*JVS( 880 )
  W( 124 ) = W( 124 ) + a*JVS( 881 )
  W( 125 ) = W( 125 ) + a*JVS( 882 )
  W( 127 ) = W( 127 ) + a*JVS( 883 )
  W( 128 ) = W( 128 ) + a*JVS( 884 )
  W( 129 ) = W( 129 ) + a*JVS( 885 )
  W( 130 ) = W( 130 ) + a*JVS( 886 )
  W( 131 ) = W( 131 ) + a*JVS( 887 )
  W( 132 ) = W( 132 ) + a*JVS( 888 )
  W( 133 ) = W( 133 ) + a*JVS( 889 )
  W( 134 ) = W( 134 ) + a*JVS( 890 )
  W( 135 ) = W( 135 ) + a*JVS( 891 )
  W( 136 ) = W( 136 ) + a*JVS( 892 )
  a = -W( 123 ) / JVS(          900  )
  W( 123 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 901 )
  W( 127 ) = W( 127 ) + a*JVS( 902 )
  W( 128 ) = W( 128 ) + a*JVS( 903 )
  W( 129 ) = W( 129 ) + a*JVS( 904 )
  W( 131 ) = W( 131 ) + a*JVS( 905 )
  W( 132 ) = W( 132 ) + a*JVS( 906 )
  W( 133 ) = W( 133 ) + a*JVS( 907 )
  W( 135 ) = W( 135 ) + a*JVS( 908 )
  W( 136 ) = W( 136 ) + a*JVS( 909 )
  a = -W( 124 ) / JVS(          913  )
  W( 124 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 914 )
  W( 128 ) = W( 128 ) + a*JVS( 915 )
  W( 129 ) = W( 129 ) + a*JVS( 916 )
  W( 131 ) = W( 131 ) + a*JVS( 917 )
  W( 132 ) = W( 132 ) + a*JVS( 918 )
  W( 133 ) = W( 133 ) + a*JVS( 919 )
  W( 135 ) = W( 135 ) + a*JVS( 920 )
  W( 136 ) = W( 136 ) + a*JVS( 921 )
  a = -W( 125 ) / JVS(          925  )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 926 )
  W( 128 ) = W( 128 ) + a*JVS( 927 )
  W( 129 ) = W( 129 ) + a*JVS( 928 )
  W( 131 ) = W( 131 ) + a*JVS( 929 )
  W( 132 ) = W( 132 ) + a*JVS( 930 )
  W( 133 ) = W( 133 ) + a*JVS( 931 )
  W( 135 ) = W( 135 ) + a*JVS( 932 )
  W( 136 ) = W( 136 ) + a*JVS( 933 )
  a = -W( 126 ) / JVS(          980  )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 981 )
  W( 128 ) = W( 128 ) + a*JVS( 982 )
  W( 129 ) = W( 129 ) + a*JVS( 983 )
  W( 130 ) = W( 130 ) + a*JVS( 984 )
  W( 131 ) = W( 131 ) + a*JVS( 985 )
  W( 132 ) = W( 132 ) + a*JVS( 986 )
  W( 133 ) = W( 133 ) + a*JVS( 987 )
  W( 134 ) = W( 134 ) + a*JVS( 988 )
  W( 135 ) = W( 135 ) + a*JVS( 989 )
  W( 136 ) = W( 136 ) + a*JVS( 990 )
  a = -W( 127 ) / JVS(          993  )
  W( 127 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 994 )
  W( 129 ) = W( 129 ) + a*JVS( 995 )
  W( 131 ) = W( 131 ) + a*JVS( 996 )
  W( 132 ) = W( 132 ) + a*JVS( 997 )
  W( 133 ) = W( 133 ) + a*JVS( 998 )
  W( 135 ) = W( 135 ) + a*JVS( 999 )
  W( 136 ) = W( 136 ) + a*JVS( 1000 )
  a = -W( 128 ) / JVS(         1034  )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 1035 )
  W( 130 ) = W( 130 ) + a*JVS( 1036 )
  W( 131 ) = W( 131 ) + a*JVS( 1037 )
  W( 132 ) = W( 132 ) + a*JVS( 1038 )
  W( 133 ) = W( 133 ) + a*JVS( 1039 )
  W( 134 ) = W( 134 ) + a*JVS( 1040 )
  W( 135 ) = W( 135 ) + a*JVS( 1041 )
  W( 136 ) = W( 136 ) + a*JVS( 1042 )
  a = -W( 129 ) / JVS(         1070  )
  W( 129 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 1071 )
  W( 131 ) = W( 131 ) + a*JVS( 1072 )
  W( 132 ) = W( 132 ) + a*JVS( 1073 )
  W( 133 ) = W( 133 ) + a*JVS( 1074 )
  W( 134 ) = W( 134 ) + a*JVS( 1075 )
  W( 135 ) = W( 135 ) + a*JVS( 1076 )
  W( 136 ) = W( 136 ) + a*JVS( 1077 )
  a = -W( 130 ) / JVS(         1086  )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 1087 )
  W( 132 ) = W( 132 ) + a*JVS( 1088 )
  W( 133 ) = W( 133 ) + a*JVS( 1089 )
  W( 134 ) = W( 134 ) + a*JVS( 1090 )
  W( 135 ) = W( 135 ) + a*JVS( 1091 )
  W( 136 ) = W( 136 ) + a*JVS( 1092 )
  a = -W( 131 ) / JVS(         1184  )
  W( 131 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 1185 )
  W( 133 ) = W( 133 ) + a*JVS( 1186 )
  W( 134 ) = W( 134 ) + a*JVS( 1187 )
  W( 135 ) = W( 135 ) + a*JVS( 1188 )
  W( 136 ) = W( 136 ) + a*JVS( 1189 )
  a = -W( 132 ) / JVS(         1210  )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1211 )
  W( 134 ) = W( 134 ) + a*JVS( 1212 )
  W( 135 ) = W( 135 ) + a*JVS( 1213 )
  W( 136 ) = W( 136 ) + a*JVS( 1214 )
  a = -W( 133 ) / JVS(         1259  )
  W( 133 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 1260 )
  W( 135 ) = W( 135 ) + a*JVS( 1261 )
  W( 136 ) = W( 136 ) + a*JVS( 1262 )
  a = -W( 134 ) / JVS(         1317  )
  W( 134 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 1318 )
  W( 136 ) = W( 136 ) + a*JVS( 1319 )
  JVS( 1320) = W( 36 )
  JVS( 1321) = W( 39 )
  JVS( 1322) = W( 52 )
  JVS( 1323) = W( 68 )
  JVS( 1324) = W( 71 )
  JVS( 1325) = W( 77 )
  JVS( 1326) = W( 89 )
  JVS( 1327) = W( 99 )
  JVS( 1328) = W( 100 )
  JVS( 1329) = W( 101 )
  JVS( 1330) = W( 102 )
  JVS( 1331) = W( 104 )
  JVS( 1332) = W( 106 )
  JVS( 1333) = W( 109 )
  JVS( 1334) = W( 110 )
  JVS( 1335) = W( 112 )
  JVS( 1336) = W( 113 )
  JVS( 1337) = W( 116 )
  JVS( 1338) = W( 117 )
  JVS( 1339) = W( 118 )
  JVS( 1340) = W( 119 )
  JVS( 1341) = W( 120 )
  JVS( 1342) = W( 121 )
  JVS( 1343) = W( 122 )
  JVS( 1344) = W( 123 )
  JVS( 1345) = W( 124 )
  JVS( 1346) = W( 125 )
  JVS( 1347) = W( 126 )
  JVS( 1348) = W( 127 )
  JVS( 1349) = W( 128 )
  JVS( 1350) = W( 129 )
  JVS( 1351) = W( 130 )
  JVS( 1352) = W( 131 )
  JVS( 1353) = W( 132 )
  JVS( 1354) = W( 133 )
  JVS( 1355) = W( 134 )
  JVS( 1356) = W( 135 )
  JVS( 1357) = W( 136 )
  IF ( ABS(  JVS( 1469 )) < TINY(a) ) THEN
         IER = 136                                     
         RETURN
  END IF
   W( 1 ) = JVS( 1358 )
   W( 11 ) = JVS( 1359 )
   W( 12 ) = JVS( 1360 )
   W( 13 ) = JVS( 1361 )
   W( 14 ) = JVS( 1362 )
   W( 16 ) = JVS( 1363 )
   W( 18 ) = JVS( 1364 )
   W( 19 ) = JVS( 1365 )
   W( 20 ) = JVS( 1366 )
   W( 22 ) = JVS( 1367 )
   W( 24 ) = JVS( 1368 )
   W( 25 ) = JVS( 1369 )
   W( 26 ) = JVS( 1370 )
   W( 27 ) = JVS( 1371 )
   W( 28 ) = JVS( 1372 )
   W( 29 ) = JVS( 1373 )
   W( 32 ) = JVS( 1374 )
   W( 33 ) = JVS( 1375 )
   W( 34 ) = JVS( 1376 )
   W( 35 ) = JVS( 1377 )
   W( 37 ) = JVS( 1378 )
   W( 38 ) = JVS( 1379 )
   W( 39 ) = JVS( 1380 )
   W( 41 ) = JVS( 1381 )
   W( 44 ) = JVS( 1382 )
   W( 45 ) = JVS( 1383 )
   W( 46 ) = JVS( 1384 )
   W( 47 ) = JVS( 1385 )
   W( 49 ) = JVS( 1386 )
   W( 50 ) = JVS( 1387 )
   W( 51 ) = JVS( 1388 )
   W( 52 ) = JVS( 1389 )
   W( 53 ) = JVS( 1390 )
   W( 54 ) = JVS( 1391 )
   W( 55 ) = JVS( 1392 )
   W( 57 ) = JVS( 1393 )
   W( 58 ) = JVS( 1394 )
   W( 61 ) = JVS( 1395 )
   W( 62 ) = JVS( 1396 )
   W( 63 ) = JVS( 1397 )
   W( 64 ) = JVS( 1398 )
   W( 65 ) = JVS( 1399 )
   W( 66 ) = JVS( 1400 )
   W( 67 ) = JVS( 1401 )
   W( 68 ) = JVS( 1402 )
   W( 69 ) = JVS( 1403 )
   W( 70 ) = JVS( 1404 )
   W( 71 ) = JVS( 1405 )
   W( 72 ) = JVS( 1406 )
   W( 73 ) = JVS( 1407 )
   W( 74 ) = JVS( 1408 )
   W( 75 ) = JVS( 1409 )
   W( 76 ) = JVS( 1410 )
   W( 77 ) = JVS( 1411 )
   W( 78 ) = JVS( 1412 )
   W( 79 ) = JVS( 1413 )
   W( 80 ) = JVS( 1414 )
   W( 81 ) = JVS( 1415 )
   W( 82 ) = JVS( 1416 )
   W( 83 ) = JVS( 1417 )
   W( 84 ) = JVS( 1418 )
   W( 85 ) = JVS( 1419 )
   W( 86 ) = JVS( 1420 )
   W( 88 ) = JVS( 1421 )
   W( 89 ) = JVS( 1422 )
   W( 90 ) = JVS( 1423 )
   W( 91 ) = JVS( 1424 )
   W( 92 ) = JVS( 1425 )
   W( 93 ) = JVS( 1426 )
   W( 94 ) = JVS( 1427 )
   W( 95 ) = JVS( 1428 )
   W( 96 ) = JVS( 1429 )
   W( 97 ) = JVS( 1430 )
   W( 98 ) = JVS( 1431 )
   W( 99 ) = JVS( 1432 )
   W( 100 ) = JVS( 1433 )
   W( 101 ) = JVS( 1434 )
   W( 102 ) = JVS( 1435 )
   W( 103 ) = JVS( 1436 )
   W( 104 ) = JVS( 1437 )
   W( 105 ) = JVS( 1438 )
   W( 106 ) = JVS( 1439 )
   W( 107 ) = JVS( 1440 )
   W( 108 ) = JVS( 1441 )
   W( 109 ) = JVS( 1442 )
   W( 110 ) = JVS( 1443 )
   W( 111 ) = JVS( 1444 )
   W( 112 ) = JVS( 1445 )
   W( 113 ) = JVS( 1446 )
   W( 114 ) = JVS( 1447 )
   W( 115 ) = JVS( 1448 )
   W( 116 ) = JVS( 1449 )
   W( 117 ) = JVS( 1450 )
   W( 118 ) = JVS( 1451 )
   W( 119 ) = JVS( 1452 )
   W( 120 ) = JVS( 1453 )
   W( 121 ) = JVS( 1454 )
   W( 122 ) = JVS( 1455 )
   W( 123 ) = JVS( 1456 )
   W( 124 ) = JVS( 1457 )
   W( 125 ) = JVS( 1458 )
   W( 126 ) = JVS( 1459 )
   W( 127 ) = JVS( 1460 )
   W( 128 ) = JVS( 1461 )
   W( 129 ) = JVS( 1462 )
   W( 130 ) = JVS( 1463 )
   W( 131 ) = JVS( 1464 )
   W( 132 ) = JVS( 1465 )
   W( 133 ) = JVS( 1466 )
   W( 134 ) = JVS( 1467 )
   W( 135 ) = JVS( 1468 )
   W( 136 ) = JVS( 1469 )
  a = -W( 1 ) / JVS(            1  )
  W( 1 ) = -a
  a = -W( 11 ) / JVS(           91  )
  W( 11 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 92 )
  a = -W( 12 ) / JVS(           93  )
  W( 12 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 94 )
  a = -W( 13 ) / JVS(           95  )
  W( 13 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 96 )
  a = -W( 14 ) / JVS(           97  )
  W( 14 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 98 )
  a = -W( 16 ) / JVS(          101  )
  W( 16 ) = -a
  W( 39 ) = W( 39 ) + a*JVS( 102 )
  W( 135 ) = W( 135 ) + a*JVS( 103 )
  W( 136 ) = W( 136 ) + a*JVS( 104 )
  a = -W( 18 ) / JVS(          108  )
  W( 18 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 109 )
  W( 136 ) = W( 136 ) + a*JVS( 110 )
  a = -W( 19 ) / JVS(          111  )
  W( 19 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 112 )
  a = -W( 20 ) / JVS(          113  )
  W( 20 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 114 )
  a = -W( 22 ) / JVS(          119  )
  W( 22 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 120 )
  a = -W( 24 ) / JVS(          124  )
  W( 24 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 125 )
  a = -W( 25 ) / JVS(          127  )
  W( 25 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 128 )
  a = -W( 26 ) / JVS(          129  )
  W( 26 ) = -a
  W( 123 ) = W( 123 ) + a*JVS( 130 )
  W( 131 ) = W( 131 ) + a*JVS( 131 )
  W( 136 ) = W( 136 ) + a*JVS( 132 )
  a = -W( 27 ) / JVS(          133  )
  W( 27 ) = -a
  W( 58 ) = W( 58 ) + a*JVS( 134 )
  W( 131 ) = W( 131 ) + a*JVS( 135 )
  W( 136 ) = W( 136 ) + a*JVS( 136 )
  a = -W( 28 ) / JVS(          137  )
  W( 28 ) = -a
  W( 81 ) = W( 81 ) + a*JVS( 138 )
  W( 131 ) = W( 131 ) + a*JVS( 139 )
  W( 136 ) = W( 136 ) + a*JVS( 140 )
  a = -W( 29 ) / JVS(          141  )
  W( 29 ) = -a
  W( 66 ) = W( 66 ) + a*JVS( 142 )
  W( 131 ) = W( 131 ) + a*JVS( 143 )
  W( 136 ) = W( 136 ) + a*JVS( 144 )
  a = -W( 32 ) / JVS(          153  )
  W( 32 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 154 )
  W( 131 ) = W( 131 ) + a*JVS( 155 )
  W( 136 ) = W( 136 ) + a*JVS( 156 )
  a = -W( 33 ) / JVS(          157  )
  W( 33 ) = -a
  W( 61 ) = W( 61 ) + a*JVS( 158 )
  W( 133 ) = W( 133 ) + a*JVS( 159 )
  W( 136 ) = W( 136 ) + a*JVS( 160 )
  a = -W( 34 ) / JVS(          161  )
  W( 34 ) = -a
  W( 69 ) = W( 69 ) + a*JVS( 162 )
  W( 131 ) = W( 131 ) + a*JVS( 163 )
  W( 136 ) = W( 136 ) + a*JVS( 164 )
  a = -W( 35 ) / JVS(          165  )
  W( 35 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 166 )
  W( 131 ) = W( 131 ) + a*JVS( 167 )
  W( 136 ) = W( 136 ) + a*JVS( 168 )
  a = -W( 37 ) / JVS(          172  )
  W( 37 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 173 )
  a = -W( 38 ) / JVS(          175  )
  W( 38 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 176 )
  a = -W( 39 ) / JVS(          177  )
  W( 39 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 178 )
  W( 136 ) = W( 136 ) + a*JVS( 179 )
  a = -W( 41 ) / JVS(          183  )
  W( 41 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 184 )
  W( 129 ) = W( 129 ) + a*JVS( 185 )
  W( 136 ) = W( 136 ) + a*JVS( 186 )
  a = -W( 44 ) / JVS(          198  )
  W( 44 ) = -a
  W( 72 ) = W( 72 ) + a*JVS( 199 )
  W( 80 ) = W( 80 ) + a*JVS( 200 )
  W( 126 ) = W( 126 ) + a*JVS( 201 )
  W( 136 ) = W( 136 ) + a*JVS( 202 )
  a = -W( 45 ) / JVS(          203  )
  W( 45 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 204 )
  W( 131 ) = W( 131 ) + a*JVS( 205 )
  W( 136 ) = W( 136 ) + a*JVS( 206 )
  a = -W( 46 ) / JVS(          207  )
  W( 46 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 208 )
  W( 131 ) = W( 131 ) + a*JVS( 209 )
  W( 136 ) = W( 136 ) + a*JVS( 210 )
  a = -W( 47 ) / JVS(          211  )
  W( 47 ) = -a
  W( 121 ) = W( 121 ) + a*JVS( 212 )
  W( 131 ) = W( 131 ) + a*JVS( 213 )
  W( 136 ) = W( 136 ) + a*JVS( 214 )
  a = -W( 49 ) / JVS(          219  )
  W( 49 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 220 )
  W( 131 ) = W( 131 ) + a*JVS( 221 )
  W( 136 ) = W( 136 ) + a*JVS( 222 )
  a = -W( 50 ) / JVS(          223  )
  W( 50 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 224 )
  W( 131 ) = W( 131 ) + a*JVS( 225 )
  W( 136 ) = W( 136 ) + a*JVS( 226 )
  a = -W( 51 ) / JVS(          227  )
  W( 51 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 228 )
  W( 131 ) = W( 131 ) + a*JVS( 229 )
  W( 136 ) = W( 136 ) + a*JVS( 230 )
  a = -W( 52 ) / JVS(          231  )
  W( 52 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 232 )
  W( 134 ) = W( 134 ) + a*JVS( 233 )
  W( 136 ) = W( 136 ) + a*JVS( 234 )
  a = -W( 53 ) / JVS(          235  )
  W( 53 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 236 )
  W( 127 ) = W( 127 ) + a*JVS( 237 )
  W( 128 ) = W( 128 ) + a*JVS( 238 )
  W( 129 ) = W( 129 ) + a*JVS( 239 )
  W( 133 ) = W( 133 ) + a*JVS( 240 )
  W( 135 ) = W( 135 ) + a*JVS( 241 )
  W( 136 ) = W( 136 ) + a*JVS( 242 )
  a = -W( 54 ) / JVS(          243  )
  W( 54 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 244 )
  W( 83 ) = W( 83 ) + a*JVS( 245 )
  W( 100 ) = W( 100 ) + a*JVS( 246 )
  W( 101 ) = W( 101 ) + a*JVS( 247 )
  W( 104 ) = W( 104 ) + a*JVS( 248 )
  W( 132 ) = W( 132 ) + a*JVS( 249 )
  W( 136 ) = W( 136 ) + a*JVS( 250 )
  a = -W( 55 ) / JVS(          251  )
  W( 55 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 252 )
  W( 136 ) = W( 136 ) + a*JVS( 253 )
  a = -W( 57 ) / JVS(          258  )
  W( 57 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 259 )
  W( 131 ) = W( 131 ) + a*JVS( 260 )
  W( 136 ) = W( 136 ) + a*JVS( 261 )
  a = -W( 58 ) / JVS(          265  )
  W( 58 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 266 )
  W( 133 ) = W( 133 ) + a*JVS( 267 )
  W( 136 ) = W( 136 ) + a*JVS( 268 )
  a = -W( 61 ) / JVS(          286  )
  W( 61 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 287 )
  W( 133 ) = W( 133 ) + a*JVS( 288 )
  W( 136 ) = W( 136 ) + a*JVS( 289 )
  a = -W( 62 ) / JVS(          290  )
  W( 62 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 291 )
  W( 113 ) = W( 113 ) + a*JVS( 292 )
  W( 133 ) = W( 133 ) + a*JVS( 293 )
  W( 136 ) = W( 136 ) + a*JVS( 294 )
  a = -W( 63 ) / JVS(          295  )
  W( 63 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 296 )
  W( 109 ) = W( 109 ) + a*JVS( 297 )
  W( 128 ) = W( 128 ) + a*JVS( 298 )
  W( 129 ) = W( 129 ) + a*JVS( 299 )
  W( 130 ) = W( 130 ) + a*JVS( 300 )
  W( 131 ) = W( 131 ) + a*JVS( 301 )
  W( 132 ) = W( 132 ) + a*JVS( 302 )
  W( 136 ) = W( 136 ) + a*JVS( 303 )
  a = -W( 64 ) / JVS(          304  )
  W( 64 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 305 )
  W( 130 ) = W( 130 ) + a*JVS( 306 )
  W( 131 ) = W( 131 ) + a*JVS( 307 )
  W( 136 ) = W( 136 ) + a*JVS( 308 )
  a = -W( 65 ) / JVS(          309  )
  W( 65 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 310 )
  W( 131 ) = W( 131 ) + a*JVS( 311 )
  W( 136 ) = W( 136 ) + a*JVS( 312 )
  a = -W( 66 ) / JVS(          315  )
  W( 66 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 316 )
  W( 133 ) = W( 133 ) + a*JVS( 317 )
  W( 136 ) = W( 136 ) + a*JVS( 318 )
  a = -W( 67 ) / JVS(          319  )
  W( 67 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 320 )
  W( 109 ) = W( 109 ) + a*JVS( 321 )
  W( 113 ) = W( 113 ) + a*JVS( 322 )
  W( 133 ) = W( 133 ) + a*JVS( 323 )
  W( 135 ) = W( 135 ) + a*JVS( 324 )
  W( 136 ) = W( 136 ) + a*JVS( 325 )
  a = -W( 68 ) / JVS(          326  )
  W( 68 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 327 )
  W( 134 ) = W( 134 ) + a*JVS( 328 )
  W( 136 ) = W( 136 ) + a*JVS( 329 )
  a = -W( 69 ) / JVS(          332  )
  W( 69 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 333 )
  W( 133 ) = W( 133 ) + a*JVS( 334 )
  W( 136 ) = W( 136 ) + a*JVS( 335 )
  a = -W( 70 ) / JVS(          336  )
  W( 70 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 337 )
  W( 131 ) = W( 131 ) + a*JVS( 338 )
  W( 136 ) = W( 136 ) + a*JVS( 339 )
  a = -W( 71 ) / JVS(          342  )
  W( 71 ) = -a
  W( 112 ) = W( 112 ) + a*JVS( 343 )
  W( 117 ) = W( 117 ) + a*JVS( 344 )
  W( 126 ) = W( 126 ) + a*JVS( 345 )
  W( 134 ) = W( 134 ) + a*JVS( 346 )
  W( 135 ) = W( 135 ) + a*JVS( 347 )
  W( 136 ) = W( 136 ) + a*JVS( 348 )
  a = -W( 72 ) / JVS(          349  )
  W( 72 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 350 )
  W( 109 ) = W( 109 ) + a*JVS( 351 )
  W( 132 ) = W( 132 ) + a*JVS( 352 )
  W( 136 ) = W( 136 ) + a*JVS( 353 )
  a = -W( 73 ) / JVS(          354  )
  W( 73 ) = -a
  W( 103 ) = W( 103 ) + a*JVS( 355 )
  W( 131 ) = W( 131 ) + a*JVS( 356 )
  W( 133 ) = W( 133 ) + a*JVS( 357 )
  W( 136 ) = W( 136 ) + a*JVS( 358 )
  a = -W( 74 ) / JVS(          362  )
  W( 74 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 363 )
  W( 133 ) = W( 133 ) + a*JVS( 364 )
  W( 134 ) = W( 134 ) + a*JVS( 365 )
  W( 136 ) = W( 136 ) + a*JVS( 366 )
  a = -W( 75 ) / JVS(          367  )
  W( 75 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 368 )
  W( 131 ) = W( 131 ) + a*JVS( 369 )
  W( 136 ) = W( 136 ) + a*JVS( 370 )
  a = -W( 76 ) / JVS(          371  )
  W( 76 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 372 )
  W( 131 ) = W( 131 ) + a*JVS( 373 )
  W( 136 ) = W( 136 ) + a*JVS( 374 )
  a = -W( 77 ) / JVS(          375  )
  W( 77 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 376 )
  W( 134 ) = W( 134 ) + a*JVS( 377 )
  W( 136 ) = W( 136 ) + a*JVS( 378 )
  a = -W( 78 ) / JVS(          382  )
  W( 78 ) = -a
  W( 82 ) = W( 82 ) + a*JVS( 383 )
  W( 83 ) = W( 83 ) + a*JVS( 384 )
  W( 89 ) = W( 89 ) + a*JVS( 385 )
  W( 100 ) = W( 100 ) + a*JVS( 386 )
  W( 101 ) = W( 101 ) + a*JVS( 387 )
  W( 104 ) = W( 104 ) + a*JVS( 388 )
  W( 126 ) = W( 126 ) + a*JVS( 389 )
  W( 131 ) = W( 131 ) + a*JVS( 390 )
  W( 132 ) = W( 132 ) + a*JVS( 391 )
  W( 133 ) = W( 133 ) + a*JVS( 392 )
  W( 136 ) = W( 136 ) + a*JVS( 393 )
  a = -W( 79 ) / JVS(          394  )
  W( 79 ) = -a
  W( 120 ) = W( 120 ) + a*JVS( 395 )
  W( 131 ) = W( 131 ) + a*JVS( 396 )
  W( 136 ) = W( 136 ) + a*JVS( 397 )
  a = -W( 80 ) / JVS(          401  )
  W( 80 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 402 )
  W( 126 ) = W( 126 ) + a*JVS( 403 )
  W( 132 ) = W( 132 ) + a*JVS( 404 )
  W( 136 ) = W( 136 ) + a*JVS( 405 )
  a = -W( 81 ) / JVS(          409  )
  W( 81 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 410 )
  W( 132 ) = W( 132 ) + a*JVS( 411 )
  W( 133 ) = W( 133 ) + a*JVS( 412 )
  W( 134 ) = W( 134 ) + a*JVS( 413 )
  W( 136 ) = W( 136 ) + a*JVS( 414 )
  a = -W( 82 ) / JVS(          415  )
  W( 82 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 416 )
  W( 136 ) = W( 136 ) + a*JVS( 417 )
  a = -W( 83 ) / JVS(          418  )
  W( 83 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 419 )
  W( 136 ) = W( 136 ) + a*JVS( 420 )
  a = -W( 84 ) / JVS(          421  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 422 )
  W( 131 ) = W( 131 ) + a*JVS( 423 )
  W( 136 ) = W( 136 ) + a*JVS( 424 )
  a = -W( 85 ) / JVS(          428  )
  W( 85 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 429 )
  W( 131 ) = W( 131 ) + a*JVS( 430 )
  W( 133 ) = W( 133 ) + a*JVS( 431 )
  W( 136 ) = W( 136 ) + a*JVS( 432 )
  a = -W( 86 ) / JVS(          434  )
  W( 86 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 435 )
  W( 131 ) = W( 131 ) + a*JVS( 436 )
  W( 133 ) = W( 133 ) + a*JVS( 437 )
  W( 136 ) = W( 136 ) + a*JVS( 438 )
  a = -W( 88 ) / JVS(          450  )
  W( 88 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 451 )
  W( 133 ) = W( 133 ) + a*JVS( 452 )
  W( 136 ) = W( 136 ) + a*JVS( 453 )
  a = -W( 89 ) / JVS(          454  )
  W( 89 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 455 )
  W( 135 ) = W( 135 ) + a*JVS( 456 )
  W( 136 ) = W( 136 ) + a*JVS( 457 )
  a = -W( 90 ) / JVS(          460  )
  W( 90 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 461 )
  W( 133 ) = W( 133 ) + a*JVS( 462 )
  W( 136 ) = W( 136 ) + a*JVS( 463 )
  a = -W( 91 ) / JVS(          467  )
  W( 91 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 468 )
  W( 133 ) = W( 133 ) + a*JVS( 469 )
  W( 134 ) = W( 134 ) + a*JVS( 470 )
  W( 136 ) = W( 136 ) + a*JVS( 471 )
  a = -W( 92 ) / JVS(          473  )
  W( 92 ) = -a
  W( 109 ) = W( 109 ) + a*JVS( 474 )
  W( 126 ) = W( 126 ) + a*JVS( 475 )
  W( 131 ) = W( 131 ) + a*JVS( 476 )
  W( 132 ) = W( 132 ) + a*JVS( 477 )
  W( 134 ) = W( 134 ) + a*JVS( 478 )
  W( 136 ) = W( 136 ) + a*JVS( 479 )
  a = -W( 93 ) / JVS(          482  )
  W( 93 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 483 )
  W( 131 ) = W( 131 ) + a*JVS( 484 )
  W( 133 ) = W( 133 ) + a*JVS( 485 )
  W( 136 ) = W( 136 ) + a*JVS( 486 )
  a = -W( 94 ) / JVS(          490  )
  W( 94 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 491 )
  W( 133 ) = W( 133 ) + a*JVS( 492 )
  W( 134 ) = W( 134 ) + a*JVS( 493 )
  W( 136 ) = W( 136 ) + a*JVS( 494 )
  a = -W( 95 ) / JVS(          495  )
  W( 95 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 496 )
  W( 111 ) = W( 111 ) + a*JVS( 497 )
  W( 118 ) = W( 118 ) + a*JVS( 498 )
  W( 119 ) = W( 119 ) + a*JVS( 499 )
  W( 120 ) = W( 120 ) + a*JVS( 500 )
  W( 121 ) = W( 121 ) + a*JVS( 501 )
  W( 123 ) = W( 123 ) + a*JVS( 502 )
  W( 125 ) = W( 125 ) + a*JVS( 503 )
  W( 127 ) = W( 127 ) + a*JVS( 504 )
  W( 129 ) = W( 129 ) + a*JVS( 505 )
  W( 136 ) = W( 136 ) + a*JVS( 506 )
  a = -W( 96 ) / JVS(          516  )
  W( 96 ) = -a
  W( 98 ) = W( 98 ) + a*JVS( 517 )
  W( 108 ) = W( 108 ) + a*JVS( 518 )
  W( 123 ) = W( 123 ) + a*JVS( 519 )
  W( 127 ) = W( 127 ) + a*JVS( 520 )
  W( 128 ) = W( 128 ) + a*JVS( 521 )
  W( 129 ) = W( 129 ) + a*JVS( 522 )
  W( 131 ) = W( 131 ) + a*JVS( 523 )
  W( 133 ) = W( 133 ) + a*JVS( 524 )
  W( 134 ) = W( 134 ) + a*JVS( 525 )
  W( 135 ) = W( 135 ) + a*JVS( 526 )
  W( 136 ) = W( 136 ) + a*JVS( 527 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 533 )
  W( 129 ) = W( 129 ) + a*JVS( 534 )
  W( 131 ) = W( 131 ) + a*JVS( 535 )
  W( 133 ) = W( 133 ) + a*JVS( 536 )
  W( 136 ) = W( 136 ) + a*JVS( 537 )
  a = -W( 98 ) / JVS(          542  )
  W( 98 ) = -a
  W( 100 ) = W( 100 ) + a*JVS( 543 )
  W( 101 ) = W( 101 ) + a*JVS( 544 )
  W( 104 ) = W( 104 ) + a*JVS( 545 )
  W( 131 ) = W( 131 ) + a*JVS( 546 )
  W( 132 ) = W( 132 ) + a*JVS( 547 )
  W( 133 ) = W( 133 ) + a*JVS( 548 )
  W( 136 ) = W( 136 ) + a*JVS( 549 )
  a = -W( 99 ) / JVS(          550  )
  W( 99 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 551 )
  W( 131 ) = W( 131 ) + a*JVS( 552 )
  W( 133 ) = W( 133 ) + a*JVS( 553 )
  W( 135 ) = W( 135 ) + a*JVS( 554 )
  a = -W( 100 ) / JVS(          555  )
  W( 100 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 556 )
  W( 135 ) = W( 135 ) + a*JVS( 557 )
  W( 136 ) = W( 136 ) + a*JVS( 558 )
  a = -W( 101 ) / JVS(          559  )
  W( 101 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 136 ) = W( 136 ) + a*JVS( 562 )
  a = -W( 102 ) / JVS(          564  )
  W( 102 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 565 )
  W( 132 ) = W( 132 ) + a*JVS( 566 )
  W( 133 ) = W( 133 ) + a*JVS( 567 )
  W( 135 ) = W( 135 ) + a*JVS( 568 )
  W( 136 ) = W( 136 ) + a*JVS( 569 )
  a = -W( 103 ) / JVS(          573  )
  W( 103 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 574 )
  W( 129 ) = W( 129 ) + a*JVS( 575 )
  W( 131 ) = W( 131 ) + a*JVS( 576 )
  W( 132 ) = W( 132 ) + a*JVS( 577 )
  W( 133 ) = W( 133 ) + a*JVS( 578 )
  W( 135 ) = W( 135 ) + a*JVS( 579 )
  W( 136 ) = W( 136 ) + a*JVS( 580 )
  a = -W( 104 ) / JVS(          581  )
  W( 104 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 582 )
  W( 135 ) = W( 135 ) + a*JVS( 583 )
  W( 136 ) = W( 136 ) + a*JVS( 584 )
  a = -W( 105 ) / JVS(          601  )
  W( 105 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 602 )
  W( 113 ) = W( 113 ) + a*JVS( 603 )
  W( 114 ) = W( 114 ) + a*JVS( 604 )
  W( 118 ) = W( 118 ) + a*JVS( 605 )
  W( 120 ) = W( 120 ) + a*JVS( 606 )
  W( 129 ) = W( 129 ) + a*JVS( 607 )
  W( 131 ) = W( 131 ) + a*JVS( 608 )
  W( 132 ) = W( 132 ) + a*JVS( 609 )
  W( 133 ) = W( 133 ) + a*JVS( 610 )
  W( 135 ) = W( 135 ) + a*JVS( 611 )
  W( 136 ) = W( 136 ) + a*JVS( 612 )
  a = -W( 106 ) / JVS(          613  )
  W( 106 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 614 )
  W( 135 ) = W( 135 ) + a*JVS( 615 )
  W( 136 ) = W( 136 ) + a*JVS( 616 )
  a = -W( 107 ) / JVS(          640  )
  W( 107 ) = -a
  W( 108 ) = W( 108 ) + a*JVS( 641 )
  W( 109 ) = W( 109 ) + a*JVS( 642 )
  W( 110 ) = W( 110 ) + a*JVS( 643 )
  W( 111 ) = W( 111 ) + a*JVS( 644 )
  W( 112 ) = W( 112 ) + a*JVS( 645 )
  W( 114 ) = W( 114 ) + a*JVS( 646 )
  W( 115 ) = W( 115 ) + a*JVS( 647 )
  W( 116 ) = W( 116 ) + a*JVS( 648 )
  W( 117 ) = W( 117 ) + a*JVS( 649 )
  W( 120 ) = W( 120 ) + a*JVS( 650 )
  W( 123 ) = W( 123 ) + a*JVS( 651 )
  W( 124 ) = W( 124 ) + a*JVS( 652 )
  W( 125 ) = W( 125 ) + a*JVS( 653 )
  W( 126 ) = W( 126 ) + a*JVS( 654 )
  W( 127 ) = W( 127 ) + a*JVS( 655 )
  W( 128 ) = W( 128 ) + a*JVS( 656 )
  W( 129 ) = W( 129 ) + a*JVS( 657 )
  W( 131 ) = W( 131 ) + a*JVS( 658 )
  W( 132 ) = W( 132 ) + a*JVS( 659 )
  W( 133 ) = W( 133 ) + a*JVS( 660 )
  W( 134 ) = W( 134 ) + a*JVS( 661 )
  W( 135 ) = W( 135 ) + a*JVS( 662 )
  W( 136 ) = W( 136 ) + a*JVS( 663 )
  a = -W( 108 ) / JVS(          668  )
  W( 108 ) = -a
  W( 111 ) = W( 111 ) + a*JVS( 669 )
  W( 120 ) = W( 120 ) + a*JVS( 670 )
  W( 123 ) = W( 123 ) + a*JVS( 671 )
  W( 125 ) = W( 125 ) + a*JVS( 672 )
  W( 127 ) = W( 127 ) + a*JVS( 673 )
  W( 128 ) = W( 128 ) + a*JVS( 674 )
  W( 129 ) = W( 129 ) + a*JVS( 675 )
  W( 131 ) = W( 131 ) + a*JVS( 676 )
  W( 132 ) = W( 132 ) + a*JVS( 677 )
  W( 133 ) = W( 133 ) + a*JVS( 678 )
  W( 135 ) = W( 135 ) + a*JVS( 679 )
  W( 136 ) = W( 136 ) + a*JVS( 680 )
  a = -W( 109 ) / JVS(          682  )
  W( 109 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 683 )
  W( 132 ) = W( 132 ) + a*JVS( 684 )
  W( 135 ) = W( 135 ) + a*JVS( 685 )
  W( 136 ) = W( 136 ) + a*JVS( 686 )
  a = -W( 110 ) / JVS(          693  )
  W( 110 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 694 )
  W( 119 ) = W( 119 ) + a*JVS( 695 )
  W( 129 ) = W( 129 ) + a*JVS( 696 )
  W( 131 ) = W( 131 ) + a*JVS( 697 )
  W( 132 ) = W( 132 ) + a*JVS( 698 )
  W( 133 ) = W( 133 ) + a*JVS( 699 )
  W( 135 ) = W( 135 ) + a*JVS( 700 )
  W( 136 ) = W( 136 ) + a*JVS( 701 )
  a = -W( 111 ) / JVS(          704  )
  W( 111 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 705 )
  W( 131 ) = W( 131 ) + a*JVS( 706 )
  W( 132 ) = W( 132 ) + a*JVS( 707 )
  W( 133 ) = W( 133 ) + a*JVS( 708 )
  W( 135 ) = W( 135 ) + a*JVS( 709 )
  W( 136 ) = W( 136 ) + a*JVS( 710 )
  a = -W( 112 ) / JVS(          722  )
  W( 112 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 723 )
  W( 124 ) = W( 124 ) + a*JVS( 724 )
  W( 129 ) = W( 129 ) + a*JVS( 725 )
  W( 131 ) = W( 131 ) + a*JVS( 726 )
  W( 132 ) = W( 132 ) + a*JVS( 727 )
  W( 133 ) = W( 133 ) + a*JVS( 728 )
  W( 135 ) = W( 135 ) + a*JVS( 729 )
  W( 136 ) = W( 136 ) + a*JVS( 730 )
  a = -W( 113 ) / JVS(          738  )
  W( 113 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 739 )
  W( 132 ) = W( 132 ) + a*JVS( 740 )
  W( 133 ) = W( 133 ) + a*JVS( 741 )
  W( 135 ) = W( 135 ) + a*JVS( 742 )
  W( 136 ) = W( 136 ) + a*JVS( 743 )
  a = -W( 114 ) / JVS(          751  )
  W( 114 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 752 )
  W( 119 ) = W( 119 ) + a*JVS( 753 )
  W( 120 ) = W( 120 ) + a*JVS( 754 )
  W( 129 ) = W( 129 ) + a*JVS( 755 )
  W( 131 ) = W( 131 ) + a*JVS( 756 )
  W( 132 ) = W( 132 ) + a*JVS( 757 )
  W( 133 ) = W( 133 ) + a*JVS( 758 )
  W( 135 ) = W( 135 ) + a*JVS( 759 )
  W( 136 ) = W( 136 ) + a*JVS( 760 )
  a = -W( 115 ) / JVS(          764  )
  W( 115 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 765 )
  W( 128 ) = W( 128 ) + a*JVS( 766 )
  W( 129 ) = W( 129 ) + a*JVS( 767 )
  W( 131 ) = W( 131 ) + a*JVS( 768 )
  W( 132 ) = W( 132 ) + a*JVS( 769 )
  W( 133 ) = W( 133 ) + a*JVS( 770 )
  W( 135 ) = W( 135 ) + a*JVS( 771 )
  W( 136 ) = W( 136 ) + a*JVS( 772 )
  a = -W( 116 ) / JVS(          775  )
  W( 116 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 776 )
  W( 119 ) = W( 119 ) + a*JVS( 777 )
  W( 120 ) = W( 120 ) + a*JVS( 778 )
  W( 125 ) = W( 125 ) + a*JVS( 779 )
  W( 127 ) = W( 127 ) + a*JVS( 780 )
  W( 129 ) = W( 129 ) + a*JVS( 781 )
  W( 131 ) = W( 131 ) + a*JVS( 782 )
  W( 132 ) = W( 132 ) + a*JVS( 783 )
  W( 133 ) = W( 133 ) + a*JVS( 784 )
  W( 135 ) = W( 135 ) + a*JVS( 785 )
  W( 136 ) = W( 136 ) + a*JVS( 786 )
  a = -W( 117 ) / JVS(          803  )
  W( 117 ) = -a
  W( 121 ) = W( 121 ) + a*JVS( 804 )
  W( 122 ) = W( 122 ) + a*JVS( 805 )
  W( 123 ) = W( 123 ) + a*JVS( 806 )
  W( 124 ) = W( 124 ) + a*JVS( 807 )
  W( 125 ) = W( 125 ) + a*JVS( 808 )
  W( 127 ) = W( 127 ) + a*JVS( 809 )
  W( 128 ) = W( 128 ) + a*JVS( 810 )
  W( 129 ) = W( 129 ) + a*JVS( 811 )
  W( 131 ) = W( 131 ) + a*JVS( 812 )
  W( 132 ) = W( 132 ) + a*JVS( 813 )
  W( 133 ) = W( 133 ) + a*JVS( 814 )
  W( 134 ) = W( 134 ) + a*JVS( 815 )
  W( 135 ) = W( 135 ) + a*JVS( 816 )
  W( 136 ) = W( 136 ) + a*JVS( 817 )
  a = -W( 118 ) / JVS(          824  )
  W( 118 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 825 )
  W( 131 ) = W( 131 ) + a*JVS( 826 )
  W( 132 ) = W( 132 ) + a*JVS( 827 )
  W( 133 ) = W( 133 ) + a*JVS( 828 )
  W( 135 ) = W( 135 ) + a*JVS( 829 )
  W( 136 ) = W( 136 ) + a*JVS( 830 )
  a = -W( 119 ) / JVS(          836  )
  W( 119 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 837 )
  W( 131 ) = W( 131 ) + a*JVS( 838 )
  W( 132 ) = W( 132 ) + a*JVS( 839 )
  W( 133 ) = W( 133 ) + a*JVS( 840 )
  W( 135 ) = W( 135 ) + a*JVS( 841 )
  W( 136 ) = W( 136 ) + a*JVS( 842 )
  a = -W( 120 ) / JVS(          847  )
  W( 120 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 848 )
  W( 131 ) = W( 131 ) + a*JVS( 849 )
  W( 132 ) = W( 132 ) + a*JVS( 850 )
  W( 133 ) = W( 133 ) + a*JVS( 851 )
  W( 135 ) = W( 135 ) + a*JVS( 852 )
  W( 136 ) = W( 136 ) + a*JVS( 853 )
  a = -W( 121 ) / JVS(          867  )
  W( 121 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 868 )
  W( 131 ) = W( 131 ) + a*JVS( 869 )
  W( 132 ) = W( 132 ) + a*JVS( 870 )
  W( 133 ) = W( 133 ) + a*JVS( 871 )
  W( 135 ) = W( 135 ) + a*JVS( 872 )
  W( 136 ) = W( 136 ) + a*JVS( 873 )
  a = -W( 122 ) / JVS(          879  )
  W( 122 ) = -a
  W( 123 ) = W( 123 ) + a*JVS( 880 )
  W( 124 ) = W( 124 ) + a*JVS( 881 )
  W( 125 ) = W( 125 ) + a*JVS( 882 )
  W( 127 ) = W( 127 ) + a*JVS( 883 )
  W( 128 ) = W( 128 ) + a*JVS( 884 )
  W( 129 ) = W( 129 ) + a*JVS( 885 )
  W( 130 ) = W( 130 ) + a*JVS( 886 )
  W( 131 ) = W( 131 ) + a*JVS( 887 )
  W( 132 ) = W( 132 ) + a*JVS( 888 )
  W( 133 ) = W( 133 ) + a*JVS( 889 )
  W( 134 ) = W( 134 ) + a*JVS( 890 )
  W( 135 ) = W( 135 ) + a*JVS( 891 )
  W( 136 ) = W( 136 ) + a*JVS( 892 )
  a = -W( 123 ) / JVS(          900  )
  W( 123 ) = -a
  W( 125 ) = W( 125 ) + a*JVS( 901 )
  W( 127 ) = W( 127 ) + a*JVS( 902 )
  W( 128 ) = W( 128 ) + a*JVS( 903 )
  W( 129 ) = W( 129 ) + a*JVS( 904 )
  W( 131 ) = W( 131 ) + a*JVS( 905 )
  W( 132 ) = W( 132 ) + a*JVS( 906 )
  W( 133 ) = W( 133 ) + a*JVS( 907 )
  W( 135 ) = W( 135 ) + a*JVS( 908 )
  W( 136 ) = W( 136 ) + a*JVS( 909 )
  a = -W( 124 ) / JVS(          913  )
  W( 124 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 914 )
  W( 128 ) = W( 128 ) + a*JVS( 915 )
  W( 129 ) = W( 129 ) + a*JVS( 916 )
  W( 131 ) = W( 131 ) + a*JVS( 917 )
  W( 132 ) = W( 132 ) + a*JVS( 918 )
  W( 133 ) = W( 133 ) + a*JVS( 919 )
  W( 135 ) = W( 135 ) + a*JVS( 920 )
  W( 136 ) = W( 136 ) + a*JVS( 921 )
  a = -W( 125 ) / JVS(          925  )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 926 )
  W( 128 ) = W( 128 ) + a*JVS( 927 )
  W( 129 ) = W( 129 ) + a*JVS( 928 )
  W( 131 ) = W( 131 ) + a*JVS( 929 )
  W( 132 ) = W( 132 ) + a*JVS( 930 )
  W( 133 ) = W( 133 ) + a*JVS( 931 )
  W( 135 ) = W( 135 ) + a*JVS( 932 )
  W( 136 ) = W( 136 ) + a*JVS( 933 )
  a = -W( 126 ) / JVS(          980  )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 981 )
  W( 128 ) = W( 128 ) + a*JVS( 982 )
  W( 129 ) = W( 129 ) + a*JVS( 983 )
  W( 130 ) = W( 130 ) + a*JVS( 984 )
  W( 131 ) = W( 131 ) + a*JVS( 985 )
  W( 132 ) = W( 132 ) + a*JVS( 986 )
  W( 133 ) = W( 133 ) + a*JVS( 987 )
  W( 134 ) = W( 134 ) + a*JVS( 988 )
  W( 135 ) = W( 135 ) + a*JVS( 989 )
  W( 136 ) = W( 136 ) + a*JVS( 990 )
  a = -W( 127 ) / JVS(          993  )
  W( 127 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 994 )
  W( 129 ) = W( 129 ) + a*JVS( 995 )
  W( 131 ) = W( 131 ) + a*JVS( 996 )
  W( 132 ) = W( 132 ) + a*JVS( 997 )
  W( 133 ) = W( 133 ) + a*JVS( 998 )
  W( 135 ) = W( 135 ) + a*JVS( 999 )
  W( 136 ) = W( 136 ) + a*JVS( 1000 )
  a = -W( 128 ) / JVS(         1034  )
  W( 128 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 1035 )
  W( 130 ) = W( 130 ) + a*JVS( 1036 )
  W( 131 ) = W( 131 ) + a*JVS( 1037 )
  W( 132 ) = W( 132 ) + a*JVS( 1038 )
  W( 133 ) = W( 133 ) + a*JVS( 1039 )
  W( 134 ) = W( 134 ) + a*JVS( 1040 )
  W( 135 ) = W( 135 ) + a*JVS( 1041 )
  W( 136 ) = W( 136 ) + a*JVS( 1042 )
  a = -W( 129 ) / JVS(         1070  )
  W( 129 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 1071 )
  W( 131 ) = W( 131 ) + a*JVS( 1072 )
  W( 132 ) = W( 132 ) + a*JVS( 1073 )
  W( 133 ) = W( 133 ) + a*JVS( 1074 )
  W( 134 ) = W( 134 ) + a*JVS( 1075 )
  W( 135 ) = W( 135 ) + a*JVS( 1076 )
  W( 136 ) = W( 136 ) + a*JVS( 1077 )
  a = -W( 130 ) / JVS(         1086  )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 1087 )
  W( 132 ) = W( 132 ) + a*JVS( 1088 )
  W( 133 ) = W( 133 ) + a*JVS( 1089 )
  W( 134 ) = W( 134 ) + a*JVS( 1090 )
  W( 135 ) = W( 135 ) + a*JVS( 1091 )
  W( 136 ) = W( 136 ) + a*JVS( 1092 )
  a = -W( 131 ) / JVS(         1184  )
  W( 131 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 1185 )
  W( 133 ) = W( 133 ) + a*JVS( 1186 )
  W( 134 ) = W( 134 ) + a*JVS( 1187 )
  W( 135 ) = W( 135 ) + a*JVS( 1188 )
  W( 136 ) = W( 136 ) + a*JVS( 1189 )
  a = -W( 132 ) / JVS(         1210  )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1211 )
  W( 134 ) = W( 134 ) + a*JVS( 1212 )
  W( 135 ) = W( 135 ) + a*JVS( 1213 )
  W( 136 ) = W( 136 ) + a*JVS( 1214 )
  a = -W( 133 ) / JVS(         1259  )
  W( 133 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 1260 )
  W( 135 ) = W( 135 ) + a*JVS( 1261 )
  W( 136 ) = W( 136 ) + a*JVS( 1262 )
  a = -W( 134 ) / JVS(         1317  )
  W( 134 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 1318 )
  W( 136 ) = W( 136 ) + a*JVS( 1319 )
  a = -W( 135 ) / JVS(         1356  )
  W( 135 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 1357 )
  JVS( 1358) = W( 1 )
  JVS( 1359) = W( 11 )
  JVS( 1360) = W( 12 )
  JVS( 1361) = W( 13 )
  JVS( 1362) = W( 14 )
  JVS( 1363) = W( 16 )
  JVS( 1364) = W( 18 )
  JVS( 1365) = W( 19 )
  JVS( 1366) = W( 20 )
  JVS( 1367) = W( 22 )
  JVS( 1368) = W( 24 )
  JVS( 1369) = W( 25 )
  JVS( 1370) = W( 26 )
  JVS( 1371) = W( 27 )
  JVS( 1372) = W( 28 )
  JVS( 1373) = W( 29 )
  JVS( 1374) = W( 32 )
  JVS( 1375) = W( 33 )
  JVS( 1376) = W( 34 )
  JVS( 1377) = W( 35 )
  JVS( 1378) = W( 37 )
  JVS( 1379) = W( 38 )
  JVS( 1380) = W( 39 )
  JVS( 1381) = W( 41 )
  JVS( 1382) = W( 44 )
  JVS( 1383) = W( 45 )
  JVS( 1384) = W( 46 )
  JVS( 1385) = W( 47 )
  JVS( 1386) = W( 49 )
  JVS( 1387) = W( 50 )
  JVS( 1388) = W( 51 )
  JVS( 1389) = W( 52 )
  JVS( 1390) = W( 53 )
  JVS( 1391) = W( 54 )
  JVS( 1392) = W( 55 )
  JVS( 1393) = W( 57 )
  JVS( 1394) = W( 58 )
  JVS( 1395) = W( 61 )
  JVS( 1396) = W( 62 )
  JVS( 1397) = W( 63 )
  JVS( 1398) = W( 64 )
  JVS( 1399) = W( 65 )
  JVS( 1400) = W( 66 )
  JVS( 1401) = W( 67 )
  JVS( 1402) = W( 68 )
  JVS( 1403) = W( 69 )
  JVS( 1404) = W( 70 )
  JVS( 1405) = W( 71 )
  JVS( 1406) = W( 72 )
  JVS( 1407) = W( 73 )
  JVS( 1408) = W( 74 )
  JVS( 1409) = W( 75 )
  JVS( 1410) = W( 76 )
  JVS( 1411) = W( 77 )
  JVS( 1412) = W( 78 )
  JVS( 1413) = W( 79 )
  JVS( 1414) = W( 80 )
  JVS( 1415) = W( 81 )
  JVS( 1416) = W( 82 )
  JVS( 1417) = W( 83 )
  JVS( 1418) = W( 84 )
  JVS( 1419) = W( 85 )
  JVS( 1420) = W( 86 )
  JVS( 1421) = W( 88 )
  JVS( 1422) = W( 89 )
  JVS( 1423) = W( 90 )
  JVS( 1424) = W( 91 )
  JVS( 1425) = W( 92 )
  JVS( 1426) = W( 93 )
  JVS( 1427) = W( 94 )
  JVS( 1428) = W( 95 )
  JVS( 1429) = W( 96 )
  JVS( 1430) = W( 97 )
  JVS( 1431) = W( 98 )
  JVS( 1432) = W( 99 )
  JVS( 1433) = W( 100 )
  JVS( 1434) = W( 101 )
  JVS( 1435) = W( 102 )
  JVS( 1436) = W( 103 )
  JVS( 1437) = W( 104 )
  JVS( 1438) = W( 105 )
  JVS( 1439) = W( 106 )
  JVS( 1440) = W( 107 )
  JVS( 1441) = W( 108 )
  JVS( 1442) = W( 109 )
  JVS( 1443) = W( 110 )
  JVS( 1444) = W( 111 )
  JVS( 1445) = W( 112 )
  JVS( 1446) = W( 113 )
  JVS( 1447) = W( 114 )
  JVS( 1448) = W( 115 )
  JVS( 1449) = W( 116 )
  JVS( 1450) = W( 117 )
  JVS( 1451) = W( 118 )
  JVS( 1452) = W( 119 )
  JVS( 1453) = W( 120 )
  JVS( 1454) = W( 121 )
  JVS( 1455) = W( 122 )
  JVS( 1456) = W( 123 )
  JVS( 1457) = W( 124 )
  JVS( 1458) = W( 125 )
  JVS( 1459) = W( 126 )
  JVS( 1460) = W( 127 )
  JVS( 1461) = W( 128 )
  JVS( 1462) = W( 129 )
  JVS( 1463) = W( 130 )
  JVS( 1464) = W( 131 )
  JVS( 1465) = W( 132 )
  JVS( 1466) = W( 133 )
  JVS( 1467) = W( 134 )
  JVS( 1468) = W( 135 )
  JVS( 1469) = W( 136 )
   
   END SUBROUTINE decomp_mozart_mosaic_4bin
 


END MODULE mozart_mosaic_4bin_Integrator
