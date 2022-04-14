
































MODULE mozart_mosaic_4bin_aq_Integrator

 USE mozart_mosaic_4bin_aq_Parameters
 USE mozart_mosaic_4bin_aq_Precision
 USE mozart_mosaic_4bin_aq_JacobianSP

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

SUBROUTINE  mozart_mosaic_4bin_aq_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE mozart_mosaic_4bin_aq_Parameters

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

   CALL mozart_mosaic_4bin_aq_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  mozart_mosaic_4bin_aq_INTEGRATE


SUBROUTINE  mozart_mosaic_4bin_aq_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE mozart_mosaic_4bin_aq_Parameters

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
      CALL mozart_mosaic_4bin_aq_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL mozart_mosaic_4bin_aq_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = mozart_mosaic_4bin_aq_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL mozart_mosaic_4bin_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL mozart_mosaic_4bin_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL mozart_mosaic_4bin_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL mozart_mosaic_4bin_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL mozart_mosaic_4bin_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL mozart_mosaic_4bin_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL mozart_mosaic_4bin_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL mozart_mosaic_4bin_aq_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL mozart_mosaic_4bin_aq_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL mozart_mosaic_4bin_aq_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL mozart_mosaic_4bin_aq_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL mozart_mosaic_4bin_aq_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL mozart_mosaic_4bin_aq_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL mozart_mosaic_4bin_aq_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL mozart_mosaic_4bin_aq_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  mozart_mosaic_4bin_aq_ros_ErrorMsg(Code,T,H,IERR)



   USE mozart_mosaic_4bin_aq_Precision

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

 END SUBROUTINE  mozart_mosaic_4bin_aq_ros_ErrorMsg


 SUBROUTINE  mozart_mosaic_4bin_aq_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL mozart_mosaic_4bin_aq_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL mozart_mosaic_4bin_aq_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL mozart_mosaic_4bin_aq_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL mozart_mosaic_4bin_aq_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL mozart_mosaic_4bin_aq_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL mozart_mosaic_4bin_aq_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL mozart_mosaic_4bin_aq_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL mozart_mosaic_4bin_aq_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL mozart_mosaic_4bin_aq_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL mozart_mosaic_4bin_aq_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL mozart_mosaic_4bin_aq_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL mozart_mosaic_4bin_aq_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL mozart_mosaic_4bin_aq_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL mozart_mosaic_4bin_aq_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL mozart_mosaic_4bin_aq_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL mozart_mosaic_4bin_aq_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL mozart_mosaic_4bin_aq_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL mozart_mosaic_4bin_aq_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL mozart_mosaic_4bin_aq_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL mozart_mosaic_4bin_aq_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = mozart_mosaic_4bin_aq_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL mozart_mosaic_4bin_aq_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  mozart_mosaic_4bin_aq_ros_Integrator



  REAL(kind=dp) FUNCTION  mozart_mosaic_4bin_aq_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    mozart_mosaic_4bin_aq_ros_ErrorNorm = Err

  END FUNCTION  mozart_mosaic_4bin_aq_ros_ErrorNorm



  SUBROUTINE mozart_mosaic_4bin_aq_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL mozart_mosaic_4bin_aq_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL mozart_mosaic_4bin_aq_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL mozart_mosaic_4bin_aq_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  mozart_mosaic_4bin_aq_ros_FunTimeDeriv



  SUBROUTINE  mozart_mosaic_4bin_aq_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL mozart_mosaic_4bin_aq_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL mozart_mosaic_4bin_aq_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL mozart_mosaic_4bin_aq_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  mozart_mosaic_4bin_aq_ros_PrepareMatrix



  SUBROUTINE  mozart_mosaic_4bin_aq_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 



CALL decomp_mozart_mosaic_4bin_aq ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  mozart_mosaic_4bin_aq_ros_Decomp



  SUBROUTINE  mozart_mosaic_4bin_aq_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL mozart_mosaic_4bin_aq_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  mozart_mosaic_4bin_aq_ros_Solve




  SUBROUTINE  mozart_mosaic_4bin_aq_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  mozart_mosaic_4bin_aq_Ros2



  SUBROUTINE  mozart_mosaic_4bin_aq_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  mozart_mosaic_4bin_aq_Ros3





  SUBROUTINE  mozart_mosaic_4bin_aq_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  mozart_mosaic_4bin_aq_Ros4


  SUBROUTINE  mozart_mosaic_4bin_aq_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  mozart_mosaic_4bin_aq_Rodas3


  SUBROUTINE  mozart_mosaic_4bin_aq_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  mozart_mosaic_4bin_aq_Rodas4




END SUBROUTINE  mozart_mosaic_4bin_aq_Rosenbrock




SUBROUTINE  mozart_mosaic_4bin_aq_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE mozart_mosaic_4bin_aq_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL mozart_mosaic_4bin_aq_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  mozart_mosaic_4bin_aq_FunTemplate



SUBROUTINE  mozart_mosaic_4bin_aq_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE mozart_mosaic_4bin_aq_Parameters
 
 USE mozart_mosaic_4bin_aq_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL mozart_mosaic_4bin_aq_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  mozart_mosaic_4bin_aq_JacTemplate

















SUBROUTINE mozart_mosaic_4bin_aq_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*F(2)
  A(2) = RCT(2)*V(136)
  A(3) = RCT(3)*V(136)
  A(4) = RCT(4)*V(19)
  A(5) = RCT(5)*V(138)
  A(6) = RCT(6)*V(40)
  A(7) = RCT(7)*V(75)
  A(8) = RCT(8)*V(133)
  A(9) = RCT(9)*V(56)
  A(10) = RCT(10)*V(53)
  A(11) = RCT(11)*V(130)
  A(12) = RCT(12)*V(130)
  A(13) = RCT(13)*V(22)
  A(14) = RCT(14)*V(114)
  A(15) = RCT(15)*V(72)
  A(16) = RCT(16)*V(68)
  A(17) = RCT(17)*V(73)
  A(18) = RCT(18)*V(81)
  A(19) = RCT(19)*V(123)
  A(20) = RCT(20)*V(128)
  A(21) = RCT(21)*V(48)
  A(22) = RCT(22)*V(49)
  A(23) = RCT(23)*V(50)
  A(24) = RCT(24)*V(110)
  A(25) = RCT(25)*V(126)
  A(26) = RCT(26)*V(29)
  A(27) = RCT(27)*V(120)
  A(28) = RCT(28)*V(80)
  A(29) = RCT(29)*V(125)
  A(30) = RCT(30)*V(112)
  A(31) = RCT(31)*V(66)
  A(32) = RCT(32)*V(101)
  A(33) = RCT(33)*V(99)
  A(34) = RCT(34)*V(88)
  A(35) = RCT(35)*V(52)
  A(36) = RCT(36)*V(61)
  A(37) = RCT(37)*V(47)
  A(38) = RCT(38)*V(25)
  A(39) = RCT(39)*V(35)
  A(40) = RCT(40)*V(34)
  A(41) = RCT(41)*V(46)
  A(42) = RCT(42)*V(44)
  A(43) = RCT(43)*V(79)
  A(44) = RCT(44)*V(106)
  A(45) = RCT(45)*V(124)
  A(46) = RCT(46)*V(115)
  A(47) = RCT(47)*V(69)
  A(48) = RCT(48)*V(83)
  A(49) = RCT(49)*V(1)
  A(50) = RCT(50)*V(96)*F(2)
  A(51) = RCT(51)*V(96)*V(136)
  A(52) = RCT(52)*V(84)*F(2)
  A(53) = RCT(53)*V(84)*F(1)
  A(54) = RCT(54)*V(55)*V(84)
  A(55) = RCT(55)*V(55)*V(134)
  A(56) = RCT(56)*V(96)*V(134)
  A(57) = RCT(57)*V(96)*V(135)
  A(58) = RCT(58)*V(134)*V(136)
  A(59) = RCT(59)*V(135)*V(136)
  A(60) = RCT(60)*V(135)*V(135)*F(1)
  A(61) = RCT(61)*V(22)*V(134)
  A(62) = RCT(62)*V(134)*V(135)
  A(63) = RCT(63)*V(134)*V(134)
  A(64) = RCT(64)*V(134)*V(134)
  A(65) = RCT(65)*V(19)*V(84)
  A(66) = RCT(66)*V(19)*V(84)
  A(67) = RCT(67)*V(135)*V(139)
  A(68) = RCT(68)*V(136)*V(139)
  A(69) = RCT(69)*V(96)*V(138)
  A(70) = RCT(70)*V(136)*V(138)
  A(71) = RCT(71)*V(133)*V(135)
  A(72) = RCT(72)*V(133)*V(138)
  A(73) = RCT(73)*V(40)
  A(74) = RCT(74)*V(134)*V(138)
  A(75) = RCT(75)*V(75)*V(134)
  A(76) = RCT(76)*V(133)*V(139)
  A(77) = RCT(77)*V(135)*V(138)
  A(78) = RCT(78)*V(56)*V(134)
  A(79) = RCT(79)*V(56)
  A(80) = RCT(80)*V(40)*F(2)
  A(81) = RCT(81)*V(133)
  A(82) = RCT(82)*V(138)
  A(83) = RCT(83)*V(77)*V(134)
  A(84) = RCT(84)*V(77)*V(84)
  A(85) = RCT(85)*V(139)*V(140)
  A(86) = RCT(86)*V(140)*V(140)
  A(87) = RCT(87)*V(140)*V(140)
  A(88) = RCT(88)*V(135)*V(140)
  A(89) = RCT(89)*V(53)*V(134)
  A(90) = RCT(90)*V(130)*V(133)
  A(91) = RCT(91)*V(130)*V(134)
  A(92) = RCT(92)*V(111)*V(134)
  A(93) = RCT(93)*V(59)*V(134)
  A(94) = RCT(94)*V(59)*V(136)
  A(95) = RCT(95)*V(20)*V(134)
  A(96) = RCT(96)*V(99)*V(134)
  A(97) = RCT(97)*V(60)*V(139)
  A(98) = RCT(98)*V(30)*F(2)
  A(99) = RCT(99)*V(30)
  A(100) = RCT(100)*V(16)*V(134)
  A(101) = RCT(101)*V(102)*V(139)
  A(102) = RCT(102)*V(102)*V(135)
  A(103) = RCT(103)*V(102)*V(140)
  A(104) = RCT(104)*V(48)*V(134)
  A(105) = RCT(105)*V(113)*V(134)
  A(106) = RCT(106)*V(113)*V(136)
  A(107) = RCT(107)*V(113)*V(133)
  A(108) = RCT(108)*V(90)*V(139)
  A(109) = RCT(109)*V(90)*V(135)
  A(110) = RCT(110)*V(72)*V(134)
  A(111) = RCT(111)*V(114)*V(134)
  A(112) = RCT(112)*V(114)*V(133)
  A(113) = RCT(113)*V(137)*V(139)
  A(114) = RCT(114)*V(137)*V(138)
  A(115) = RCT(115)*V(135)*V(137)
  A(116) = RCT(116)*V(137)*V(140)
  A(117) = RCT(117)*V(68)*V(134)
  A(118) = RCT(118)*V(73)
  A(119) = RCT(119)*V(137)*V(137)
  A(120) = RCT(120)*V(17)*V(134)
  A(121) = RCT(121)*V(97)*V(139)
  A(122) = RCT(122)*V(97)*V(135)
  A(123) = RCT(123)*V(97)*V(140)
  A(124) = RCT(124)*V(49)*V(134)
  A(125) = RCT(125)*V(110)*V(134)
  A(126) = RCT(126)*V(121)*V(139)
  A(127) = RCT(127)*V(121)*V(135)
  A(128) = RCT(128)*V(121)*V(140)
  A(129) = RCT(129)*V(50)*V(134)
  A(130) = RCT(130)*V(18)*V(134)
  A(131) = RCT(131)*V(18)*V(134)
  A(132) = RCT(132)*V(18)*V(134)
  A(133) = RCT(133)*V(18)*V(134)
  A(134) = RCT(134)*V(18)*V(134)
  A(135) = RCT(135)*V(51)*V(139)
  A(136) = RCT(136)*V(58)*V(134)
  A(137) = RCT(137)*V(58)*V(134)
  A(138) = RCT(138)*V(58)*V(134)
  A(139) = RCT(139)*V(58)*V(134)
  A(140) = RCT(140)*V(58)*V(134)
  A(141) = RCT(141)*V(117)*V(139)
  A(142) = RCT(142)*V(117)*V(135)
  A(143) = RCT(143)*V(88)*V(134)
  A(144) = RCT(144)*V(71)*V(134)
  A(145) = RCT(145)*V(66)*V(134)
  A(146) = RCT(146)*V(89)*V(139)
  A(147) = RCT(147)*V(89)*V(135)
  A(148) = RCT(148)*V(52)*V(134)
  A(149) = RCT(149)*V(93)*V(135)
  A(150) = RCT(150)*V(61)*V(134)
  A(151) = RCT(151)*V(108)*V(134)
  A(152) = RCT(152)*V(108)*V(134)
  A(153) = RCT(153)*V(108)*V(134)
  A(154) = RCT(154)*V(108)*V(134)
  A(155) = RCT(155)*V(108)*V(134)
  A(156) = RCT(156)*V(108)*V(136)
  A(157) = RCT(157)*V(108)*V(136)
  A(158) = RCT(158)*V(108)*V(136)
  A(159) = RCT(159)*V(108)*V(136)
  A(160) = RCT(160)*V(108)*V(136)
  A(161) = RCT(161)*V(131)*V(139)
  A(162) = RCT(162)*V(131)*V(133)
  A(163) = RCT(163)*V(131)*V(135)
  A(164) = RCT(164)*V(80)*V(134)
  A(165) = RCT(165)*V(131)*V(140)
  A(166) = RCT(166)*V(131)*V(137)
  A(167) = RCT(167)*V(128)*V(134)
  A(168) = RCT(168)*V(128)*V(136)
  A(169) = RCT(169)*V(123)*V(134)
  A(170) = RCT(170)*V(123)*V(136)
  A(171) = RCT(171)*V(129)*V(139)
  A(172) = RCT(172)*V(129)*V(139)
  A(173) = RCT(173)*V(129)*V(133)
  A(174) = RCT(174)*V(129)*V(135)
  A(175) = RCT(175)*V(129)*V(140)
  A(176) = RCT(176)*V(129)*V(137)
  A(177) = RCT(177)*V(54)*V(134)
  A(178) = RCT(178)*V(132)*V(139)
  A(179) = RCT(179)*V(132)*V(133)
  A(180) = RCT(180)*V(132)*V(135)
  A(181) = RCT(181)*V(132)*V(140)
  A(182) = RCT(182)*V(132)*V(137)
  A(183) = RCT(183)*V(132)*V(132)
  A(184) = RCT(184)*V(132)*V(138)*F(2)
  A(185) = RCT(185)*V(81)*F(2)
  A(186) = RCT(186)*V(24)*V(134)
  A(187) = RCT(187)*V(24)*V(134)
  A(188) = RCT(188)*V(24)*V(134)
  A(189) = RCT(189)*V(24)*V(134)
  A(190) = RCT(190)*V(24)*V(134)
  A(191) = RCT(191)*V(26)*V(134)
  A(192) = RCT(192)*V(62)*V(139)
  A(193) = RCT(193)*V(62)*V(135)
  A(194) = RCT(194)*V(31)*V(134)
  A(195) = RCT(195)*V(63)*V(138)
  A(196) = RCT(196)*V(63)*V(136)
  A(197) = RCT(197)*V(85)*V(139)
  A(198) = RCT(198)*V(85)*V(135)
  A(199) = RCT(199)*V(32)*V(134)
  A(200) = RCT(200)*V(70)*V(139)
  A(201) = RCT(201)*V(70)*V(135)
  A(202) = RCT(202)*V(33)*V(134)
  A(203) = RCT(203)*V(91)*V(138)
  A(204) = RCT(204)*V(91)*V(139)
  A(205) = RCT(205)*V(91)*V(135)
  A(206) = RCT(206)*V(27)*V(134)
  A(207) = RCT(207)*V(27)*V(134)
  A(208) = RCT(208)*V(27)*V(134)
  A(209) = RCT(209)*V(27)*V(134)
  A(210) = RCT(210)*V(27)*V(134)
  A(211) = RCT(211)*V(28)*V(134)
  A(212) = RCT(212)*V(64)*V(135)
  A(213) = RCT(213)*V(36)*V(134)
  A(214) = RCT(214)*V(64)*V(139)
  A(215) = RCT(215)*V(37)*V(134)
  A(216) = RCT(216)*V(78)*V(138)
  A(217) = RCT(217)*V(21)
  A(218) = RCT(218)*V(78)*V(139)
  A(219) = RCT(219)*V(78)*V(135)
  A(220) = RCT(220)*V(93)*V(139)
  A(221) = RCT(221)*V(95)*V(135)
  A(222) = RCT(222)*V(95)*V(139)
  A(223) = RCT(223)*V(98)*V(135)
  A(224) = RCT(224)*V(98)*V(139)
  A(225) = RCT(225)*V(95)*V(138)
  A(226) = RCT(226)*V(98)*V(138)
  A(227) = RCT(227)*V(41)*V(134)
  A(228) = RCT(228)*V(41)*V(134)
  A(229) = RCT(229)*V(41)*V(134)
  A(230) = RCT(230)*V(41)*V(134)
  A(231) = RCT(231)*V(41)*V(134)
  A(232) = RCT(232)*V(42)*V(134)
  A(233) = RCT(233)*V(74)*V(139)
  A(234) = RCT(234)*V(74)*V(135)
  A(235) = RCT(235)*V(38)*V(134)
  A(236) = RCT(236)*V(94)*V(135)
  A(237) = RCT(237)*V(39)*V(134)
  A(238) = RCT(238)*V(94)*V(139)
  A(239) = RCT(239)*V(87)*V(134)
  A(240) = RCT(240)*V(87)*V(134)
  A(241) = RCT(241)*V(87)*V(134)
  A(242) = RCT(242)*V(87)*V(134)
  A(243) = RCT(243)*V(87)*V(134)
  A(244) = RCT(244)*V(86)*V(134)
  A(245) = RCT(245)*V(86)*V(134)
  A(246) = RCT(246)*V(86)*V(134)
  A(247) = RCT(247)*V(86)*V(134)
  A(248) = RCT(248)*V(86)*V(134)
  A(249) = RCT(249)*V(104)*V(134)
  A(250) = RCT(250)*V(104)*V(134)
  A(251) = RCT(251)*V(104)*V(134)
  A(252) = RCT(252)*V(104)*V(134)
  A(253) = RCT(253)*V(104)*V(134)
  A(254) = RCT(254)*V(107)*V(134)
  A(255) = RCT(255)*V(109)*V(134)
  A(256) = RCT(256)*V(87)*V(136)
  A(257) = RCT(257)*V(87)*V(136)
  A(258) = RCT(258)*V(87)*V(136)
  A(259) = RCT(259)*V(87)*V(136)
  A(260) = RCT(260)*V(87)*V(136)
  A(261) = RCT(261)*V(86)*V(136)
  A(262) = RCT(262)*V(86)*V(136)
  A(263) = RCT(263)*V(86)*V(136)
  A(264) = RCT(264)*V(86)*V(136)
  A(265) = RCT(265)*V(86)*V(136)
  A(266) = RCT(266)*V(104)*V(136)
  A(267) = RCT(267)*V(104)*V(136)
  A(268) = RCT(268)*V(104)*V(136)
  A(269) = RCT(269)*V(104)*V(136)
  A(270) = RCT(270)*V(104)*V(136)
  A(271) = RCT(271)*V(107)*V(136)
  A(272) = RCT(272)*V(109)*V(136)
  A(273) = RCT(273)*V(104)*V(133)
  A(274) = RCT(274)*V(107)*V(133)
  A(275) = RCT(275)*V(109)*V(133)
  A(276) = RCT(276)*V(119)*V(139)
  A(277) = RCT(277)*V(119)*V(135)
  A(278) = RCT(278)*V(119)*V(140)
  A(279) = RCT(279)*V(69)*V(134)
  A(280) = RCT(280)*V(83)*V(134)
  A(281) = RCT(281)*V(124)*V(134)
  A(282) = RCT(282)*V(124)*V(133)
  A(283) = RCT(283)*V(118)*V(139)
  A(284) = RCT(284)*V(118)*V(135)
  A(285) = RCT(285)*V(118)*V(140)
  A(286) = RCT(286)*V(115)*V(134)
  A(287) = RCT(287)*V(122)*V(139)
  A(288) = RCT(288)*V(122)*V(135)
  A(289) = RCT(289)*V(122)*V(140)
  A(290) = RCT(290)*V(122)*V(133)
  A(291) = RCT(291)*V(67)*V(134)
  A(292) = RCT(292)*V(108)*V(133)
  A(293) = RCT(293)*V(103)*V(139)
  A(294) = RCT(294)*V(103)*V(133)
  A(295) = RCT(295)*V(103)*V(135)
  A(296) = RCT(296)*V(126)*V(134)
  A(297) = RCT(297)*V(126)*V(133)
  A(298) = RCT(298)*V(120)*V(134)
  A(299) = RCT(299)*V(120)*V(133)
  A(300) = RCT(300)*V(57)*V(134)
  A(301) = RCT(301)*V(127)*V(139)
  A(302) = RCT(302)*V(127)*V(133)
  A(303) = RCT(303)*V(127)*V(135)
  A(304) = RCT(304)*V(127)*V(140)
  A(305) = RCT(305)*V(127)*V(137)
  A(306) = RCT(306)*V(29)*V(134)
  A(307) = RCT(307)*V(29)*V(134)
  A(308) = RCT(308)*V(100)*V(134)
  A(309) = RCT(309)*V(45)*V(134)
  A(310) = RCT(310)*V(81)*V(134)
  A(311) = RCT(311)*V(73)*V(134)
  A(312) = RCT(312)*V(125)*V(134)
  A(313) = RCT(313)*V(112)*V(134)
  A(314) = RCT(314)*V(43)*V(134)
  A(315) = RCT(315)*V(43)*V(134)
  A(316) = RCT(316)*V(43)*V(133)
  A(317) = RCT(317)*V(15)*V(134)
  A(318) = RCT(318)*V(135)
  A(319) = RCT(319)*V(102)*V(102)
  A(320) = RCT(320)*V(92)*V(134)
  A(321) = RCT(321)*V(116)*V(139)
  A(322) = RCT(322)*V(116)*V(140)
  A(323) = RCT(323)*V(106)*V(134)
  A(324) = RCT(324)*V(76)*V(139)
  A(325) = RCT(325)*V(76)*V(135)
  A(326) = RCT(326)*V(116)*V(135)
  A(327) = RCT(327)*V(79)*V(134)
  A(328) = RCT(328)*V(92)*V(136)
  A(329) = RCT(329)*V(92)*V(133)
  A(330) = RCT(330)*V(105)*V(135)
  A(331) = RCT(331)*V(105)*V(139)
  A(332) = RCT(332)*V(105)*V(133)
  A(333) = RCT(333)*V(23)*V(134)
  A(334) = 4.5e-13*V(82)*V(134)
  A(335) = RCT(335)*V(130)*V(135)
  A(336) = RCT(336)*V(65)
  A(337) = RCT(337)*V(65)*V(139)
  A(338) = RCT(338)*V(65)*V(135)
  A(339) = RCT(339)*V(9)*V(134)
  A(340) = RCT(340)*V(8)*V(134)
  A(341) = RCT(341)*V(7)*V(134)
  A(342) = RCT(342)*V(6)*V(134)
  A(343) = RCT(343)*V(14)*V(134)
  A(344) = RCT(344)*V(13)*V(134)
  A(345) = RCT(345)*V(12)*V(134)
  A(346) = RCT(346)*V(11)*V(134)


  Vdot(1) = -A(49)
  Vdot(2) = A(95)
  Vdot(3) = A(85)+A(97)+A(101)+A(108)+A(113)+A(121)+A(126)+A(135)+A(141)+A(146)+A(161)+A(171)+A(172)+A(178)+A(192)&
              &+A(197)+A(200)+A(204)+A(214)+A(218)+A(220)+A(222)+A(224)+A(233)+A(238)+A(276)+A(283)+A(287)+A(293)+A(301)&
              &+A(321)+A(324)+A(331)+A(337)
  Vdot(4) = A(86)+A(87)+A(88)+A(102)+A(103)+A(109)+A(115)+A(116)+A(119)+A(122)+A(123)+A(127)+A(128)+A(142)+A(147)+A(149)&
              &+A(163)+A(165)+A(166)+A(174)+A(175)+A(176)+A(180)+A(181)+A(182)+A(183)+A(193)+A(198)+A(201)+A(205)+A(212)&
              &+A(219)+A(221)+A(223)+A(234)+A(236)+A(277)+A(278)+A(284)+A(285)+A(288)+A(289)+A(295)+A(303)+A(304)+A(305)&
              &+A(319)+A(322)+A(325)+A(326)+A(330)+A(338)
  Vdot(5) = 1.075*A(342)
  Vdot(6) = A(134)+A(140)+A(190)+A(210)+A(231)+1.075*A(341)-A(342)
  Vdot(7) = A(133)+A(139)+A(189)+A(209)+A(230)+1.075*A(340)-A(341)
  Vdot(8) = A(132)+A(138)+A(188)+A(208)+A(229)+1.075*A(339)-A(340)
  Vdot(9) = A(131)+A(137)+A(187)+A(207)+A(228)-A(339)
  Vdot(10) = 1.075*A(346)
  Vdot(11) = A(155)+A(160)+A(243)+A(248)+A(253)+A(260)+A(265)+A(270)+1.075*A(345)-A(346)
  Vdot(12) = A(154)+A(159)+A(242)+A(247)+A(252)+A(259)+A(264)+A(269)+1.075*A(344)-A(345)
  Vdot(13) = A(153)+A(158)+A(241)+A(246)+A(251)+A(258)+A(263)+A(268)+1.075*A(343)-A(344)
  Vdot(14) = A(152)+A(157)+A(240)+A(245)+A(250)+A(257)+A(262)+A(267)-A(343)
  Vdot(15) = -A(317)
  Vdot(16) = -A(100)
  Vdot(17) = -A(120)
  Vdot(18) = -A(130)
  Vdot(19) = -A(4)-A(65)-A(66)
  Vdot(20) = -A(95)+A(314)+0.5*A(315)+A(316)
  Vdot(21) = A(216)-A(217)
  Vdot(22) = -A(13)+A(60)-A(61)+A(64)+0.5*A(318)
  Vdot(23) = -A(333)
  Vdot(24) = -A(186)
  Vdot(25) = -A(38)+0.12*A(186)
  Vdot(26) = 0.53*A(186)-A(191)
  Vdot(27) = -A(206)
  Vdot(28) = 0.18*A(206)-A(211)
  Vdot(29) = -A(26)+A(303)-A(306)-A(307)
  Vdot(30) = A(97)-A(98)-A(99)
  Vdot(31) = A(193)-A(194)
  Vdot(32) = A(198)-A(199)
  Vdot(33) = A(201)-A(202)
  Vdot(34) = -A(40)+0.2*A(220)+0.2*A(238)
  Vdot(35) = -A(39)+0.1*A(206)+0.23*A(227)
  Vdot(36) = A(212)-A(213)
  Vdot(37) = A(214)-A(215)
  Vdot(38) = A(234)-A(235)
  Vdot(39) = A(236)-A(237)
  Vdot(40) = -A(6)+A(72)-A(73)-A(80)
  Vdot(41) = -A(227)
  Vdot(42) = 0.15*A(227)-A(232)
  Vdot(43) = -A(314)-A(315)-A(316)
  Vdot(44) = -A(42)+0.21*A(238)
  Vdot(45) = 0.2*A(103)-A(309)+0.4*A(319)
  Vdot(46) = -A(41)+0.2*A(220)+0.15*A(238)
  Vdot(47) = -A(37)+A(38)+0.5*A(200)+0.2*A(220)+0.06*A(238)
  Vdot(48) = -A(21)+A(102)-A(104)
  Vdot(49) = -A(22)+A(122)-A(124)
  Vdot(50) = -A(23)+A(127)-A(129)
  Vdot(51) = A(130)-A(135)
  Vdot(52) = -A(35)+A(147)-A(148)
  Vdot(53) = -A(10)+A(88)-A(89)
  Vdot(54) = A(174)-A(177)
  Vdot(55) = A(12)-A(54)-A(55)+0.05*A(84)
  Vdot(56) = -A(9)+A(77)-A(78)-A(79)
  Vdot(57) = 0.33*A(161)+0.36*A(162)+0.3*A(165)+0.4*A(166)+A(298)+A(299)-A(300)
  Vdot(58) = -A(136)+0.05*A(256)+0.05*A(261)+0.05*A(266)+0.05*A(271)+0.05*A(272)
  Vdot(59) = -A(93)-A(94)
  Vdot(60) = 0.75*A(93)-A(97)
  Vdot(61) = -A(36)+A(149)-A(150)
  Vdot(62) = 0.14*A(191)-A(192)-A(193)+A(194)+0.2*A(211)
  Vdot(63) = 0.06*A(191)-A(195)-A(196)+A(197)+0.07*A(211)+0.07*A(232)
  Vdot(64) = 0.07*A(206)-A(212)+A(213)-A(214)+0.06*A(227)
  Vdot(65) = A(335)-A(336)-A(337)-A(338)
  Vdot(66) = -A(31)+0.8*A(34)+0.75*A(141)-A(145)
  Vdot(67) = 0.25*A(106)+0.25*A(115)+0.1*A(116)+0.2*A(156)+0.25*A(180)-A(291)
  Vdot(68) = -A(16)+0.75*A(115)-A(117)+0.75*A(180)
  Vdot(69) = -A(47)+A(277)-A(279)
  Vdot(70) = 0.35*A(186)-A(200)-A(201)+A(202)
  Vdot(71) = A(107)+0.1*A(141)-A(144)+0.75*A(331)+0.75*A(332)
  Vdot(72) = -A(15)+A(109)-A(110)
  Vdot(73) = -A(17)+A(114)-A(118)-A(311)
  Vdot(74) = 0.3*A(232)-A(233)-A(234)+A(235)
  Vdot(75) = -A(7)+A(74)-A(75)+2*A(80)+A(81)+0.5*A(82)+A(90)+A(112)+A(297)+A(316)
  Vdot(76) = A(323)-A(324)-A(325)
  Vdot(77) = -A(83)-A(84)+0.08*A(106)
  Vdot(78) = A(215)-A(216)+A(217)-A(218)-A(219)
  Vdot(79) = -A(43)+A(326)-A(327)
  Vdot(80) = -A(28)+A(163)-A(164)
  Vdot(81) = -A(18)+A(184)-A(185)-A(310)
  Vdot(82) = 0.5*A(94)+0.05*A(256)+0.05*A(261)+0.05*A(266)+0.05*A(271)+0.05*A(272)+0.25*A(328)+0.35*A(333)-A(334)+A(337)&
               &+A(338)
  Vdot(83) = -A(48)-A(280)+A(284)
  Vdot(84) = A(2)+A(4)-A(52)-A(53)-A(54)-A(65)-A(66)-A(84)
  Vdot(85) = A(196)-A(197)-A(198)+A(199)+A(218)+0.4*A(219)
  Vdot(86) = -A(244)-A(261)
  Vdot(87) = -A(239)-A(256)
  Vdot(88) = -A(34)+A(142)-A(143)
  Vdot(89) = A(145)-A(146)-A(147)+A(148)
  Vdot(90) = A(105)-A(108)-A(109)+0.5*A(110)
  Vdot(91) = 0.6*A(37)-A(203)-A(204)-A(205)
  Vdot(92) = -A(320)-A(328)-A(329)
  Vdot(93) = -A(149)+A(150)+0.65*A(206)-A(220)
  Vdot(94) = 0.56*A(227)-A(236)+A(237)-A(238)
  Vdot(95) = 0.6*A(40)-A(221)-A(222)-A(225)
  Vdot(96) = 2*A(1)+A(3)+A(5)-A(50)-A(51)+A(52)-A(56)-A(57)+A(63)-A(69)
  Vdot(97) = A(120)-A(121)-A(122)-A(123)+A(124)
  Vdot(98) = 0.6*A(41)-A(223)-A(224)-A(226)
  Vdot(99) = 0.13*A(32)-A(33)+0.45*A(36)-A(96)+0.02*A(161)+0.02*A(162)+0.7*A(192)+A(200)+0.4*A(204)+0.16*A(205)+0.6&
               &*A(220)+0.07*A(223)+0.17*A(224)+0.17*A(233)+0.34*A(238)+0.25*A(301)+0.25*A(302)+0.1*A(304)+0.25*A(305)+0.2&
               &*A(313)+0.65*A(333)
  Vdot(100) = A(87)+0.3*A(103)+0.5*A(128)+0.25*A(165)+0.25*A(175)+0.25*A(278)+0.25*A(285)+0.25*A(289)+0.3*A(304)-A(308)&
                &+0.25*A(322)
  Vdot(101) = -A(32)+0.9*A(36)+0.1*A(256)+0.1*A(261)+0.1*A(266)+0.1*A(271)+0.1*A(272)
  Vdot(102) = A(31)+A(100)-A(101)-A(102)-A(103)+0.5*A(104)-2*A(319)
  Vdot(103) = A(292)-A(293)-A(294)-A(295)
  Vdot(104) = -A(249)-A(266)-A(273)
  Vdot(105) = A(329)-A(330)-A(331)-A(332)
  Vdot(106) = 0.33*A(43)-A(44)+0.33*A(321)+0.167*A(322)-A(323)+0.9*A(328)+0.25*A(331)+0.25*A(332)
  Vdot(107) = -A(254)-A(271)-A(274)
  Vdot(108) = -A(151)-A(156)-A(292)
  Vdot(109) = -A(255)-A(272)-A(275)
  Vdot(110) = 0.82*A(22)-A(24)+0.25*A(34)+0.67*A(43)+A(44)+0.5*A(46)+0.05*A(47)+0.3*A(48)+0.82*A(121)+0.82*A(123)-A(125)&
                &+0.5*A(135)+0.25*A(141)+0.52*A(256)+0.52*A(261)+0.52*A(266)+0.52*A(271)+0.52*A(272)+0.045*A(276)+0.025&
                &*A(278)+0.27*A(283)+0.15*A(285)+0.5*A(286)+0.67*A(321)+0.333*A(322)+A(324)+0.4*A(325)+0.1*A(328)+0.75&
                &*A(331)+0.75*A(332)
  Vdot(111) = A(11)+A(12)+A(14)+0.67*A(19)+0.7*A(20)+A(25)+A(27)+A(30)+0.45*A(32)+2*A(33)+1.5*A(38)+1.5*A(39)+0.6*A(41)&
                &+A(42)+A(44)+A(45)+1.7*A(46)+0.25*A(48)+A(90)+A(91)-A(92)+0.5*A(94)+A(96)+0.56*A(106)+0.3*A(156)+0.05&
                &*A(168)+0.2*A(170)+0.22*A(171)+0.22*A(173)+0.11*A(175)+0.22*A(176)+0.4*A(204)+0.16*A(205)+0.07*A(221)+0.17&
                &*A(222)+0.14*A(223)+0.35*A(224)+0.23*A(256)+0.23*A(261)+0.23*A(266)+0.23*A(271)+0.23*A(272)+0.225*A(283)&
                &+0.125*A(285)+0.7*A(286)+A(296)+A(297)+0.25*A(301)+0.25*A(302)+0.2*A(304)+0.25*A(305)+0.35*A(328)+0.35&
                &*A(333)
  Vdot(112) = -A(30)+0.67*A(43)+0.25*A(48)+A(98)+0.02*A(161)+0.02*A(162)+0.53*A(171)+0.53*A(173)+0.26*A(175)+0.53*A(176)&
                &+0.225*A(283)+0.125*A(285)+0.25*A(301)+0.25*A(302)+0.1*A(304)+0.25*A(305)-A(313)+0.67*A(321)+0.333*A(322)
  Vdot(113) = 0.7*A(20)-A(105)-A(106)-A(107)+0.07*A(156)
  Vdot(114) = -A(14)+A(15)+A(21)+0.4*A(34)+A(35)+A(101)+0.8*A(103)+0.5*A(104)+0.5*A(106)+A(108)-A(111)-A(112)+0.27&
                &*A(121)+A(135)+0.4*A(141)+A(146)+0.04*A(168)+A(309)+1.6*A(319)
  Vdot(115) = A(45)-A(46)+A(48)+0.3*A(256)+0.3*A(261)+0.3*A(266)+0.3*A(271)+0.3*A(272)+0.9*A(283)+A(285)-A(286)
  Vdot(116) = A(320)-A(321)-A(322)-A(326)+0.5*A(327)
  Vdot(117) = A(136)-A(141)-A(142)+A(143)
  Vdot(118) = A(280)+A(281)+0.5*A(282)-A(283)-A(284)-A(285)
  Vdot(119) = A(239)+A(244)+A(249)+A(254)+A(255)-A(276)-A(277)-A(278)+A(279)
  Vdot(120) = -A(27)+0.08*A(161)+0.8*A(172)+0.1*A(276)+0.1*A(283)+0.1*A(287)+A(288)+0.5*A(289)+0.794*A(293)+0.794*A(294)&
                &+0.794*A(295)-A(298)-A(299)
  Vdot(121) = 0.15*A(46)+A(125)-A(126)-A(127)-A(128)+A(129)+0.06*A(256)+0.06*A(261)+0.06*A(266)+0.06*A(271)+0.06*A(272)&
                &+0.15*A(286)
  Vdot(122) = A(273)+A(274)+A(275)+0.5*A(282)-A(287)-A(288)-A(289)-A(290)
  Vdot(123) = -A(19)+0.288*A(28)+0.4*A(156)+0.23*A(161)+0.25*A(162)+0.19*A(165)+0.25*A(166)-A(169)-A(170)+0.167*A(293)&
                &+0.167*A(294)+0.167*A(295)
  Vdot(124) = -A(45)+A(47)+0.33*A(256)+0.33*A(261)+0.33*A(266)+0.33*A(271)+0.33*A(272)+0.9*A(276)+A(278)-A(281)-A(282)&
                &+0.9*A(287)+0.5*A(289)+A(290)
  Vdot(125) = -A(29)+0.5*A(110)+0.2*A(128)+0.02*A(161)+0.02*A(162)+0.22*A(171)+0.22*A(173)+0.23*A(175)+0.22*A(176)+0.25&
                &*A(301)+0.25*A(302)+0.1*A(304)+0.25*A(305)+0.5*A(310)-A(312)
  Vdot(126) = -A(25)+0.18*A(32)+0.45*A(36)+A(42)+0.5*A(128)+A(144)+0.02*A(161)+0.02*A(162)+0.95*A(168)+0.8*A(170)+0.25&
                &*A(171)+0.25*A(173)+0.24*A(175)+0.25*A(176)+0.4*A(220)+0.07*A(221)+0.17*A(222)+0.07*A(223)+0.17*A(224)+0.51&
                &*A(233)+0.54*A(238)-A(296)-A(297)+0.25*A(301)+0.25*A(302)+0.1*A(304)+0.25*A(305)+A(312)
  Vdot(127) = 0.5*A(164)+A(300)-A(301)-A(302)-A(303)-A(304)-A(305)+A(306)
  Vdot(128) = -A(20)+0.402*A(28)+0.2*A(156)+0.32*A(161)+0.35*A(162)+0.26*A(165)+0.35*A(166)-A(167)-A(168)+0.039*A(293)&
                &+0.039*A(294)+0.039*A(295)
  Vdot(129) = A(167)+0.5*A(169)-A(171)-A(172)-A(173)-A(174)-A(175)-A(176)+0.2*A(177)
  Vdot(130) = A(10)-A(11)-A(12)+A(15)+0.67*A(19)+A(23)+A(27)+0.69*A(28)+A(29)+A(30)+0.1*A(34)+0.33*A(43)+0.68*A(46)+0.4&
                &*A(47)+0.375*A(48)+0.25*A(84)+A(85)+2*A(86)+A(87)+0.3*A(89)-A(90)-A(91)+0.5*A(93)+A(94)+2*A(99)+0.7*A(103)&
                &+0.54*A(106)+A(108)+A(116)+0.5*A(117)+A(123)+A(126)+0.8*A(128)+0.5*A(135)+0.1*A(141)+0.6*A(156)+0.55*A(161)&
                &+0.6*A(162)+1.2*A(165)+0.6*A(166)+0.8*A(168)+0.7*A(170)+0.25*A(171)+0.25*A(173)+0.88*A(175)+0.25*A(176)&
                &+A(178)+A(179)+2*A(181)+A(182)+2*A(183)+0.34*A(256)+0.34*A(261)+0.34*A(266)+0.34*A(271)+0.34*A(272)+0.36&
                &*A(276)+0.95*A(278)+0.34*A(283)+0.93*A(285)+0.68*A(286)+0.75*A(289)+0.072*A(293)+0.072*A(294)+0.206*A(295)&
                &+0.25*A(301)+0.25*A(302)+0.8*A(304)+0.25*A(305)+A(308)+0.5*A(310)+A(311)+0.8*A(313)+0.33*A(321)+0.917&
                &*A(322)+0.5*A(328)+0.25*A(331)+0.25*A(332)-A(335)+A(336)
  Vdot(131) = A(151)-A(161)-A(162)-A(163)+0.5*A(164)-A(165)-A(166)
  Vdot(132) = A(18)+0.33*A(19)+0.2*A(156)+0.5*A(169)+0.5*A(177)-A(178)-A(179)-A(180)-A(181)-A(182)-2*A(183)-A(184)&
                &+A(185)
  Vdot(133) = A(6)-A(8)+0.33*A(9)+0.4*A(17)+A(70)-A(71)-A(72)+A(73)+A(75)-A(76)-A(81)-A(90)-A(107)-A(112)-A(162)-A(173)&
                &-A(179)-A(273)-A(274)-A(275)-A(282)-A(290)-A(292)-A(294)-A(297)-A(299)-A(302)+0.5*A(310)+A(311)-A(316)&
                &-A(329)-A(332)
  Vdot(134) = A(7)+0.33*A(9)+A(10)+2*A(13)+A(15)+A(16)+0.33*A(19)+A(21)+A(22)+A(23)+A(26)+A(34)+A(35)+A(36)+A(43)+A(47)&
                &+A(48)+A(49)+2*A(53)+A(54)-A(55)-A(56)+A(57)-A(58)+A(59)-A(61)-A(62)-2*A(63)-2*A(64)+A(67)+A(71)-A(74)&
                &-A(75)-A(78)+0.5*A(82)-A(83)+0.75*A(84)-0.7*A(89)-A(91)-A(92)-A(93)+0.12*A(94)-A(95)-A(96)-A(100)-0.5&
                &*A(104)-A(105)+0.33*A(106)-0.5*A(110)-A(111)-A(117)-A(120)-A(124)-A(125)-A(129)-A(130)-A(136)-A(143)-A(144)&
                &-A(145)-A(148)-A(150)-A(151)+0.27*A(156)-A(164)-A(167)+0.08*A(168)-A(169)+0.215*A(170)-0.9*A(177)-A(186)&
                &-A(191)-A(194)-A(199)-A(202)-A(206)-A(211)-A(213)-A(215)+0.4*A(219)+0.4*A(221)+0.4*A(223)-A(227)-A(232)&
                &-A(235)-A(237)-A(239)-A(244)-A(249)-A(254)-A(255)+0.63*A(256)+0.63*A(261)+0.63*A(266)+0.63*A(271)+0.63&
                &*A(272)-A(279)-A(280)-A(281)-A(286)-A(291)+0.206*A(295)-A(296)-A(298)-A(300)-A(306)-A(308)-A(309)-A(310)&
                &-A(311)-A(312)-A(313)-A(314)-A(315)-A(317)-A(320)-A(323)+0.4*A(325)-0.5*A(327)+0.06*A(328)-0.35*A(333)&
                &-A(334)
  Vdot(135) = 0.66*A(9)+A(10)+2*A(11)+A(14)+A(15)+0.67*A(19)+A(21)+A(22)+A(25)+A(27)+A(28)+A(29)+2*A(30)+0.56*A(32)+2&
                &*A(33)+0.9*A(34)+A(37)+1.5*A(38)+A(39)+0.6*A(40)+0.6*A(41)+A(42)+A(43)+2*A(44)+A(45)+1.2*A(46)+A(47)+A(48)&
                &+A(54)+A(55)+A(56)-A(57)+A(58)-A(59)-2*A(60)+A(61)-A(62)-A(67)-A(71)-A(77)+A(79)+0.4*A(84)+A(85)+2*A(86)&
                &-A(88)+A(90)+A(91)+A(92)+0.25*A(93)+0.12*A(94)+A(96)+A(98)+A(99)+A(101)-A(102)+A(103)+0.19*A(106)+A(108)&
                &-A(109)-A(115)+0.9*A(116)+A(121)-A(122)+A(123)-A(127)+0.3*A(128)+A(135)+0.9*A(141)-A(142)-A(147)-A(149)&
                &+0.06*A(156)+0.92*A(161)+A(162)-A(163)+A(165)+A(166)+0.06*A(168)+0.275*A(170)+0.47*A(171)+0.47*A(173)&
                &-A(174)+0.73*A(175)+0.47*A(176)+0.2*A(177)-A(180)+A(181)+0.65*A(186)+0.8*A(191)+A(192)-A(193)-A(198)+A(200)&
                &-A(201)+0.4*A(204)-0.84*A(205)+0.28*A(206)+0.73*A(211)-A(212)+A(214)-A(219)+A(220)-0.93*A(221)+0.17*A(222)&
                &-0.67*A(223)+0.83*A(224)+0.38*A(227)+0.63*A(232)+A(233)-A(234)-A(236)+A(238)+0.57*A(256)+0.57*A(261)+0.57&
                &*A(266)+0.57*A(271)+0.57*A(272)+0.9*A(276)-A(277)+A(278)+0.9*A(283)-A(284)+A(285)+0.2*A(286)-A(288)+0.5&
                &*A(289)+0.794*A(293)+0.794*A(294)-A(295)+A(298)+A(299)+A(301)+A(302)-A(303)+0.8*A(304)+A(305)+A(308)+A(309)&
                &+0.5*A(310)+A(312)+A(313)+0.5*A(315)-A(318)+1.2*A(319)+A(321)+A(322)+A(324)-0.6*A(325)-A(326)+0.06*A(328)&
                &-A(330)+0.75*A(331)+0.75*A(332)+0.35*A(333)+A(334)-A(335)+A(336)+A(337)-A(338)
  Vdot(136) = -A(2)-A(3)+0.89*A(8)+A(50)-A(51)-A(58)-A(59)-A(68)-A(70)-A(94)-A(106)+0.25*A(115)-0.9*A(156)-0.8*A(168)&
                &-0.8*A(170)+0.25*A(180)-A(196)-A(256)-A(261)-A(266)-A(271)-A(272)-A(328)
  Vdot(137) = 0.6*A(17)+0.67*A(19)+0.3*A(20)+A(23)+A(24)+A(25)+A(29)+A(31)+0.13*A(32)+A(35)+0.5*A(39)+A(42)+0.65*A(46)&
                &+A(111)+A(112)-A(113)-A(114)-A(115)-A(116)+0.5*A(117)+A(118)-2*A(119)+A(126)+0.3*A(128)+A(146)-A(166)+0.53&
                &*A(171)+0.53*A(173)+0.26*A(175)-0.47*A(176)+A(178)+A(179)+A(181)+2*A(183)+0.06*A(256)+0.06*A(261)+0.06&
                &*A(266)+0.06*A(271)+0.06*A(272)+0.65*A(286)+A(296)+A(297)-A(305)
  Vdot(138) = -A(5)+A(6)+A(7)+0.89*A(8)+0.66*A(9)+0.6*A(17)+A(18)+A(27)+A(67)+A(68)-A(69)-A(70)+A(71)-A(72)+A(73)-A(74)&
                &+2*A(76)-A(77)+A(78)+A(79)-A(82)+A(85)+A(97)+A(101)+A(108)+A(113)-A(114)+A(118)+A(121)+A(126)+A(135)+0.9&
                &*A(141)+A(144)+A(146)+0.92*A(161)+A(162)+A(171)+A(173)+A(178)+A(179)-A(184)+A(185)+A(192)-A(195)+A(197)&
                &+A(200)-A(203)+A(214)-A(216)+A(217)+A(218)+A(220)+A(222)+A(224)-A(225)-A(226)+A(233)+A(238)+0.9*A(276)+0.9&
                &*A(283)+1.9*A(287)+0.5*A(289)+2*A(290)+1.206*A(293)+1.206*A(294)+0.206*A(295)+0.4*A(298)+A(299)+A(301)&
                &+A(302)+A(321)+A(324)+1.25*A(331)+1.25*A(332)+A(337)
  Vdot(139) = A(5)+0.11*A(8)+A(49)+2*A(65)-A(67)-A(68)+A(69)-A(76)+0.5*A(82)-A(85)-A(97)-A(101)-A(108)-A(113)-A(121)&
                &-A(126)-A(135)-A(141)-A(146)-A(161)-A(171)-A(172)-A(178)-A(192)-A(197)-A(200)-A(204)-A(214)-A(218)-A(220)&
                &-A(222)-A(224)-A(233)-A(238)-A(276)-A(283)-A(287)-A(293)-A(301)-A(321)-A(324)-A(331)-A(337)
  Vdot(140) = A(14)+A(16)+0.4*A(17)+0.3*A(20)+A(24)+A(83)+0.75*A(84)-A(85)-2*A(86)-2*A(87)-A(88)+0.7*A(89)-A(103)+0.31&
                &*A(106)+A(113)-0.1*A(116)+2*A(119)-A(123)-A(128)-A(165)+A(166)-A(175)+A(176)-A(181)+A(182)+0.33*A(221)+0.83&
                &*A(222)+0.07*A(223)+0.17*A(224)-A(278)-A(285)-A(289)+A(291)-A(304)+A(305)-A(322)
      
END SUBROUTINE mozart_mosaic_4bin_aq_Fun
















SUBROUTINE mozart_mosaic_4bin_aq_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*F(2)
  IRR(2) = RCT(2)*V(136)
  IRR(3) = RCT(3)*V(136)
  IRR(4) = RCT(4)*V(19)
  IRR(5) = RCT(5)*V(138)
  IRR(6) = RCT(6)*V(40)
  IRR(7) = RCT(7)*V(75)
  IRR(8) = RCT(8)*V(133)
  IRR(9) = RCT(9)*V(56)
  IRR(10) = RCT(10)*V(53)
  IRR(11) = RCT(11)*V(130)
  IRR(12) = RCT(12)*V(130)
  IRR(13) = RCT(13)*V(22)
  IRR(14) = RCT(14)*V(114)
  IRR(15) = RCT(15)*V(72)
  IRR(16) = RCT(16)*V(68)
  IRR(17) = RCT(17)*V(73)
  IRR(18) = RCT(18)*V(81)
  IRR(19) = RCT(19)*V(123)
  IRR(20) = RCT(20)*V(128)
  IRR(21) = RCT(21)*V(48)
  IRR(22) = RCT(22)*V(49)
  IRR(23) = RCT(23)*V(50)
  IRR(24) = RCT(24)*V(110)
  IRR(25) = RCT(25)*V(126)
  IRR(26) = RCT(26)*V(29)
  IRR(27) = RCT(27)*V(120)
  IRR(28) = RCT(28)*V(80)
  IRR(29) = RCT(29)*V(125)
  IRR(30) = RCT(30)*V(112)
  IRR(31) = RCT(31)*V(66)
  IRR(32) = RCT(32)*V(101)
  IRR(33) = RCT(33)*V(99)
  IRR(34) = RCT(34)*V(88)
  IRR(35) = RCT(35)*V(52)
  IRR(36) = RCT(36)*V(61)
  IRR(37) = RCT(37)*V(47)
  IRR(38) = RCT(38)*V(25)
  IRR(39) = RCT(39)*V(35)
  IRR(40) = RCT(40)*V(34)
  IRR(41) = RCT(41)*V(46)
  IRR(42) = RCT(42)*V(44)
  IRR(43) = RCT(43)*V(79)
  IRR(44) = RCT(44)*V(106)
  IRR(45) = RCT(45)*V(124)
  IRR(46) = RCT(46)*V(115)
  IRR(47) = RCT(47)*V(69)
  IRR(48) = RCT(48)*V(83)
  IRR(49) = RCT(49)*V(1)
  IRR(50) = RCT(50)*V(96)*F(2)
  IRR(51) = RCT(51)*V(96)*V(136)
  IRR(52) = RCT(52)*V(84)*F(2)
  IRR(53) = RCT(53)*V(84)*F(1)
  IRR(54) = RCT(54)*V(55)*V(84)
  IRR(55) = RCT(55)*V(55)*V(134)
  IRR(56) = RCT(56)*V(96)*V(134)
  IRR(57) = RCT(57)*V(96)*V(135)
  IRR(58) = RCT(58)*V(134)*V(136)
  IRR(59) = RCT(59)*V(135)*V(136)
  IRR(60) = RCT(60)*V(135)*V(135)*F(1)
  IRR(61) = RCT(61)*V(22)*V(134)
  IRR(62) = RCT(62)*V(134)*V(135)
  IRR(63) = RCT(63)*V(134)*V(134)
  IRR(64) = RCT(64)*V(134)*V(134)
  IRR(65) = RCT(65)*V(19)*V(84)
  IRR(66) = RCT(66)*V(19)*V(84)
  IRR(67) = RCT(67)*V(135)*V(139)
  IRR(68) = RCT(68)*V(136)*V(139)
  IRR(69) = RCT(69)*V(96)*V(138)
  IRR(70) = RCT(70)*V(136)*V(138)
  IRR(71) = RCT(71)*V(133)*V(135)
  IRR(72) = RCT(72)*V(133)*V(138)
  IRR(73) = RCT(73)*V(40)
  IRR(74) = RCT(74)*V(134)*V(138)
  IRR(75) = RCT(75)*V(75)*V(134)
  IRR(76) = RCT(76)*V(133)*V(139)
  IRR(77) = RCT(77)*V(135)*V(138)
  IRR(78) = RCT(78)*V(56)*V(134)
  IRR(79) = RCT(79)*V(56)
  IRR(80) = RCT(80)*V(40)*F(2)
  IRR(81) = RCT(81)*V(133)
  IRR(82) = RCT(82)*V(138)
  IRR(83) = RCT(83)*V(77)*V(134)
  IRR(84) = RCT(84)*V(77)*V(84)
  IRR(85) = RCT(85)*V(139)*V(140)
  IRR(86) = RCT(86)*V(140)*V(140)
  IRR(87) = RCT(87)*V(140)*V(140)
  IRR(88) = RCT(88)*V(135)*V(140)
  IRR(89) = RCT(89)*V(53)*V(134)
  IRR(90) = RCT(90)*V(130)*V(133)
  IRR(91) = RCT(91)*V(130)*V(134)
  IRR(92) = RCT(92)*V(111)*V(134)
  IRR(93) = RCT(93)*V(59)*V(134)
  IRR(94) = RCT(94)*V(59)*V(136)
  IRR(95) = RCT(95)*V(20)*V(134)
  IRR(96) = RCT(96)*V(99)*V(134)
  IRR(97) = RCT(97)*V(60)*V(139)
  IRR(98) = RCT(98)*V(30)*F(2)
  IRR(99) = RCT(99)*V(30)
  IRR(100) = RCT(100)*V(16)*V(134)
  IRR(101) = RCT(101)*V(102)*V(139)
  IRR(102) = RCT(102)*V(102)*V(135)
  IRR(103) = RCT(103)*V(102)*V(140)
  IRR(104) = RCT(104)*V(48)*V(134)
  IRR(105) = RCT(105)*V(113)*V(134)
  IRR(106) = RCT(106)*V(113)*V(136)
  IRR(107) = RCT(107)*V(113)*V(133)
  IRR(108) = RCT(108)*V(90)*V(139)
  IRR(109) = RCT(109)*V(90)*V(135)
  IRR(110) = RCT(110)*V(72)*V(134)
  IRR(111) = RCT(111)*V(114)*V(134)
  IRR(112) = RCT(112)*V(114)*V(133)
  IRR(113) = RCT(113)*V(137)*V(139)
  IRR(114) = RCT(114)*V(137)*V(138)
  IRR(115) = RCT(115)*V(135)*V(137)
  IRR(116) = RCT(116)*V(137)*V(140)
  IRR(117) = RCT(117)*V(68)*V(134)
  IRR(118) = RCT(118)*V(73)
  IRR(119) = RCT(119)*V(137)*V(137)
  IRR(120) = RCT(120)*V(17)*V(134)
  IRR(121) = RCT(121)*V(97)*V(139)
  IRR(122) = RCT(122)*V(97)*V(135)
  IRR(123) = RCT(123)*V(97)*V(140)
  IRR(124) = RCT(124)*V(49)*V(134)
  IRR(125) = RCT(125)*V(110)*V(134)
  IRR(126) = RCT(126)*V(121)*V(139)
  IRR(127) = RCT(127)*V(121)*V(135)
  IRR(128) = RCT(128)*V(121)*V(140)
  IRR(129) = RCT(129)*V(50)*V(134)
  IRR(130) = RCT(130)*V(18)*V(134)
  IRR(131) = RCT(131)*V(18)*V(134)
  IRR(132) = RCT(132)*V(18)*V(134)
  IRR(133) = RCT(133)*V(18)*V(134)
  IRR(134) = RCT(134)*V(18)*V(134)
  IRR(135) = RCT(135)*V(51)*V(139)
  IRR(136) = RCT(136)*V(58)*V(134)
  IRR(137) = RCT(137)*V(58)*V(134)
  IRR(138) = RCT(138)*V(58)*V(134)
  IRR(139) = RCT(139)*V(58)*V(134)
  IRR(140) = RCT(140)*V(58)*V(134)
  IRR(141) = RCT(141)*V(117)*V(139)
  IRR(142) = RCT(142)*V(117)*V(135)
  IRR(143) = RCT(143)*V(88)*V(134)
  IRR(144) = RCT(144)*V(71)*V(134)
  IRR(145) = RCT(145)*V(66)*V(134)
  IRR(146) = RCT(146)*V(89)*V(139)
  IRR(147) = RCT(147)*V(89)*V(135)
  IRR(148) = RCT(148)*V(52)*V(134)
  IRR(149) = RCT(149)*V(93)*V(135)
  IRR(150) = RCT(150)*V(61)*V(134)
  IRR(151) = RCT(151)*V(108)*V(134)
  IRR(152) = RCT(152)*V(108)*V(134)
  IRR(153) = RCT(153)*V(108)*V(134)
  IRR(154) = RCT(154)*V(108)*V(134)
  IRR(155) = RCT(155)*V(108)*V(134)
  IRR(156) = RCT(156)*V(108)*V(136)
  IRR(157) = RCT(157)*V(108)*V(136)
  IRR(158) = RCT(158)*V(108)*V(136)
  IRR(159) = RCT(159)*V(108)*V(136)
  IRR(160) = RCT(160)*V(108)*V(136)
  IRR(161) = RCT(161)*V(131)*V(139)
  IRR(162) = RCT(162)*V(131)*V(133)
  IRR(163) = RCT(163)*V(131)*V(135)
  IRR(164) = RCT(164)*V(80)*V(134)
  IRR(165) = RCT(165)*V(131)*V(140)
  IRR(166) = RCT(166)*V(131)*V(137)
  IRR(167) = RCT(167)*V(128)*V(134)
  IRR(168) = RCT(168)*V(128)*V(136)
  IRR(169) = RCT(169)*V(123)*V(134)
  IRR(170) = RCT(170)*V(123)*V(136)
  IRR(171) = RCT(171)*V(129)*V(139)
  IRR(172) = RCT(172)*V(129)*V(139)
  IRR(173) = RCT(173)*V(129)*V(133)
  IRR(174) = RCT(174)*V(129)*V(135)
  IRR(175) = RCT(175)*V(129)*V(140)
  IRR(176) = RCT(176)*V(129)*V(137)
  IRR(177) = RCT(177)*V(54)*V(134)
  IRR(178) = RCT(178)*V(132)*V(139)
  IRR(179) = RCT(179)*V(132)*V(133)
  IRR(180) = RCT(180)*V(132)*V(135)
  IRR(181) = RCT(181)*V(132)*V(140)
  IRR(182) = RCT(182)*V(132)*V(137)
  IRR(183) = RCT(183)*V(132)*V(132)
  IRR(184) = RCT(184)*V(132)*V(138)*F(2)
  IRR(185) = RCT(185)*V(81)*F(2)
  IRR(186) = RCT(186)*V(24)*V(134)
  IRR(187) = RCT(187)*V(24)*V(134)
  IRR(188) = RCT(188)*V(24)*V(134)
  IRR(189) = RCT(189)*V(24)*V(134)
  IRR(190) = RCT(190)*V(24)*V(134)
  IRR(191) = RCT(191)*V(26)*V(134)
  IRR(192) = RCT(192)*V(62)*V(139)
  IRR(193) = RCT(193)*V(62)*V(135)
  IRR(194) = RCT(194)*V(31)*V(134)
  IRR(195) = RCT(195)*V(63)*V(138)
  IRR(196) = RCT(196)*V(63)*V(136)
  IRR(197) = RCT(197)*V(85)*V(139)
  IRR(198) = RCT(198)*V(85)*V(135)
  IRR(199) = RCT(199)*V(32)*V(134)
  IRR(200) = RCT(200)*V(70)*V(139)
  IRR(201) = RCT(201)*V(70)*V(135)
  IRR(202) = RCT(202)*V(33)*V(134)
  IRR(203) = RCT(203)*V(91)*V(138)
  IRR(204) = RCT(204)*V(91)*V(139)
  IRR(205) = RCT(205)*V(91)*V(135)
  IRR(206) = RCT(206)*V(27)*V(134)
  IRR(207) = RCT(207)*V(27)*V(134)
  IRR(208) = RCT(208)*V(27)*V(134)
  IRR(209) = RCT(209)*V(27)*V(134)
  IRR(210) = RCT(210)*V(27)*V(134)
  IRR(211) = RCT(211)*V(28)*V(134)
  IRR(212) = RCT(212)*V(64)*V(135)
  IRR(213) = RCT(213)*V(36)*V(134)
  IRR(214) = RCT(214)*V(64)*V(139)
  IRR(215) = RCT(215)*V(37)*V(134)
  IRR(216) = RCT(216)*V(78)*V(138)
  IRR(217) = RCT(217)*V(21)
  IRR(218) = RCT(218)*V(78)*V(139)
  IRR(219) = RCT(219)*V(78)*V(135)
  IRR(220) = RCT(220)*V(93)*V(139)
  IRR(221) = RCT(221)*V(95)*V(135)
  IRR(222) = RCT(222)*V(95)*V(139)
  IRR(223) = RCT(223)*V(98)*V(135)
  IRR(224) = RCT(224)*V(98)*V(139)
  IRR(225) = RCT(225)*V(95)*V(138)
  IRR(226) = RCT(226)*V(98)*V(138)
  IRR(227) = RCT(227)*V(41)*V(134)
  IRR(228) = RCT(228)*V(41)*V(134)
  IRR(229) = RCT(229)*V(41)*V(134)
  IRR(230) = RCT(230)*V(41)*V(134)
  IRR(231) = RCT(231)*V(41)*V(134)
  IRR(232) = RCT(232)*V(42)*V(134)
  IRR(233) = RCT(233)*V(74)*V(139)
  IRR(234) = RCT(234)*V(74)*V(135)
  IRR(235) = RCT(235)*V(38)*V(134)
  IRR(236) = RCT(236)*V(94)*V(135)
  IRR(237) = RCT(237)*V(39)*V(134)
  IRR(238) = RCT(238)*V(94)*V(139)
  IRR(239) = RCT(239)*V(87)*V(134)
  IRR(240) = RCT(240)*V(87)*V(134)
  IRR(241) = RCT(241)*V(87)*V(134)
  IRR(242) = RCT(242)*V(87)*V(134)
  IRR(243) = RCT(243)*V(87)*V(134)
  IRR(244) = RCT(244)*V(86)*V(134)
  IRR(245) = RCT(245)*V(86)*V(134)
  IRR(246) = RCT(246)*V(86)*V(134)
  IRR(247) = RCT(247)*V(86)*V(134)
  IRR(248) = RCT(248)*V(86)*V(134)
  IRR(249) = RCT(249)*V(104)*V(134)
  IRR(250) = RCT(250)*V(104)*V(134)
  IRR(251) = RCT(251)*V(104)*V(134)
  IRR(252) = RCT(252)*V(104)*V(134)
  IRR(253) = RCT(253)*V(104)*V(134)
  IRR(254) = RCT(254)*V(107)*V(134)
  IRR(255) = RCT(255)*V(109)*V(134)
  IRR(256) = RCT(256)*V(87)*V(136)
  IRR(257) = RCT(257)*V(87)*V(136)
  IRR(258) = RCT(258)*V(87)*V(136)
  IRR(259) = RCT(259)*V(87)*V(136)
  IRR(260) = RCT(260)*V(87)*V(136)
  IRR(261) = RCT(261)*V(86)*V(136)
  IRR(262) = RCT(262)*V(86)*V(136)
  IRR(263) = RCT(263)*V(86)*V(136)
  IRR(264) = RCT(264)*V(86)*V(136)
  IRR(265) = RCT(265)*V(86)*V(136)
  IRR(266) = RCT(266)*V(104)*V(136)
  IRR(267) = RCT(267)*V(104)*V(136)
  IRR(268) = RCT(268)*V(104)*V(136)
  IRR(269) = RCT(269)*V(104)*V(136)
  IRR(270) = RCT(270)*V(104)*V(136)
  IRR(271) = RCT(271)*V(107)*V(136)
  IRR(272) = RCT(272)*V(109)*V(136)
  IRR(273) = RCT(273)*V(104)*V(133)
  IRR(274) = RCT(274)*V(107)*V(133)
  IRR(275) = RCT(275)*V(109)*V(133)
  IRR(276) = RCT(276)*V(119)*V(139)
  IRR(277) = RCT(277)*V(119)*V(135)
  IRR(278) = RCT(278)*V(119)*V(140)
  IRR(279) = RCT(279)*V(69)*V(134)
  IRR(280) = RCT(280)*V(83)*V(134)
  IRR(281) = RCT(281)*V(124)*V(134)
  IRR(282) = RCT(282)*V(124)*V(133)
  IRR(283) = RCT(283)*V(118)*V(139)
  IRR(284) = RCT(284)*V(118)*V(135)
  IRR(285) = RCT(285)*V(118)*V(140)
  IRR(286) = RCT(286)*V(115)*V(134)
  IRR(287) = RCT(287)*V(122)*V(139)
  IRR(288) = RCT(288)*V(122)*V(135)
  IRR(289) = RCT(289)*V(122)*V(140)
  IRR(290) = RCT(290)*V(122)*V(133)
  IRR(291) = RCT(291)*V(67)*V(134)
  IRR(292) = RCT(292)*V(108)*V(133)
  IRR(293) = RCT(293)*V(103)*V(139)
  IRR(294) = RCT(294)*V(103)*V(133)
  IRR(295) = RCT(295)*V(103)*V(135)
  IRR(296) = RCT(296)*V(126)*V(134)
  IRR(297) = RCT(297)*V(126)*V(133)
  IRR(298) = RCT(298)*V(120)*V(134)
  IRR(299) = RCT(299)*V(120)*V(133)
  IRR(300) = RCT(300)*V(57)*V(134)
  IRR(301) = RCT(301)*V(127)*V(139)
  IRR(302) = RCT(302)*V(127)*V(133)
  IRR(303) = RCT(303)*V(127)*V(135)
  IRR(304) = RCT(304)*V(127)*V(140)
  IRR(305) = RCT(305)*V(127)*V(137)
  IRR(306) = RCT(306)*V(29)*V(134)
  IRR(307) = RCT(307)*V(29)*V(134)
  IRR(308) = RCT(308)*V(100)*V(134)
  IRR(309) = RCT(309)*V(45)*V(134)
  IRR(310) = RCT(310)*V(81)*V(134)
  IRR(311) = RCT(311)*V(73)*V(134)
  IRR(312) = RCT(312)*V(125)*V(134)
  IRR(313) = RCT(313)*V(112)*V(134)
  IRR(314) = RCT(314)*V(43)*V(134)
  IRR(315) = RCT(315)*V(43)*V(134)
  IRR(316) = RCT(316)*V(43)*V(133)
  IRR(317) = RCT(317)*V(15)*V(134)
  IRR(318) = RCT(318)*V(135)
  IRR(319) = RCT(319)*V(102)*V(102)
  IRR(320) = RCT(320)*V(92)*V(134)
  IRR(321) = RCT(321)*V(116)*V(139)
  IRR(322) = RCT(322)*V(116)*V(140)
  IRR(323) = RCT(323)*V(106)*V(134)
  IRR(324) = RCT(324)*V(76)*V(139)
  IRR(325) = RCT(325)*V(76)*V(135)
  IRR(326) = RCT(326)*V(116)*V(135)
  IRR(327) = RCT(327)*V(79)*V(134)
  IRR(328) = RCT(328)*V(92)*V(136)
  IRR(329) = RCT(329)*V(92)*V(133)
  IRR(330) = RCT(330)*V(105)*V(135)
  IRR(331) = RCT(331)*V(105)*V(139)
  IRR(332) = RCT(332)*V(105)*V(133)
  IRR(333) = RCT(333)*V(23)*V(134)
  IRR(334) = 4.5e-13*V(82)*V(134)
  IRR(335) = RCT(335)*V(130)*V(135)
  IRR(336) = RCT(336)*V(65)
  IRR(337) = RCT(337)*V(65)*V(139)
  IRR(338) = RCT(338)*V(65)*V(135)
  IRR(339) = RCT(339)*V(9)*V(134)
  IRR(340) = RCT(340)*V(8)*V(134)
  IRR(341) = RCT(341)*V(7)*V(134)
  IRR(342) = RCT(342)*V(6)*V(134)
  IRR(343) = RCT(343)*V(14)*V(134)
  IRR(344) = RCT(344)*V(13)*V(134)
  IRR(345) = RCT(345)*V(12)*V(134)
  IRR(346) = RCT(346)*V(11)*V(134)
      
END SUBROUTINE mozart_mosaic_4bin_aq_IRRFun
















SUBROUTINE mozart_mosaic_4bin_aq_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(628)


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

  B(52) = RCT(51)*V(136)

  B(53) = RCT(51)*V(96)

  B(54) = RCT(52)*F(2)

  B(56) = RCT(53)*F(1)

  B(58) = RCT(54)*V(84)

  B(59) = RCT(54)*V(55)

  B(60) = RCT(55)*V(134)

  B(61) = RCT(55)*V(55)

  B(62) = RCT(56)*V(134)

  B(63) = RCT(56)*V(96)

  B(64) = RCT(57)*V(135)

  B(65) = RCT(57)*V(96)

  B(66) = RCT(58)*V(136)

  B(67) = RCT(58)*V(134)

  B(68) = RCT(59)*V(136)

  B(69) = RCT(59)*V(135)

  B(70) = RCT(60)*2*V(135)*F(1)

  B(72) = RCT(61)*V(134)

  B(73) = RCT(61)*V(22)

  B(74) = RCT(62)*V(135)

  B(75) = RCT(62)*V(134)

  B(76) = RCT(63)*2*V(134)

  B(77) = RCT(64)*2*V(134)

  B(78) = RCT(65)*V(84)

  B(79) = RCT(65)*V(19)

  B(80) = RCT(66)*V(84)

  B(81) = RCT(66)*V(19)

  B(82) = RCT(67)*V(139)

  B(83) = RCT(67)*V(135)

  B(84) = RCT(68)*V(139)

  B(85) = RCT(68)*V(136)

  B(86) = RCT(69)*V(138)

  B(87) = RCT(69)*V(96)

  B(88) = RCT(70)*V(138)

  B(89) = RCT(70)*V(136)

  B(90) = RCT(71)*V(135)

  B(91) = RCT(71)*V(133)

  B(92) = RCT(72)*V(138)

  B(93) = RCT(72)*V(133)

  B(94) = RCT(73)

  B(95) = RCT(74)*V(138)

  B(96) = RCT(74)*V(134)

  B(97) = RCT(75)*V(134)

  B(98) = RCT(75)*V(75)

  B(99) = RCT(76)*V(139)

  B(100) = RCT(76)*V(133)

  B(101) = RCT(77)*V(138)

  B(102) = RCT(77)*V(135)

  B(103) = RCT(78)*V(134)

  B(104) = RCT(78)*V(56)

  B(105) = RCT(79)

  B(106) = RCT(80)*F(2)

  B(108) = RCT(81)

  B(109) = RCT(82)

  B(110) = RCT(83)*V(134)

  B(111) = RCT(83)*V(77)

  B(112) = RCT(84)*V(84)

  B(113) = RCT(84)*V(77)

  B(114) = RCT(85)*V(140)

  B(115) = RCT(85)*V(139)

  B(116) = RCT(86)*2*V(140)

  B(117) = RCT(87)*2*V(140)

  B(118) = RCT(88)*V(140)

  B(119) = RCT(88)*V(135)

  B(120) = RCT(89)*V(134)

  B(121) = RCT(89)*V(53)

  B(122) = RCT(90)*V(133)

  B(123) = RCT(90)*V(130)

  B(124) = RCT(91)*V(134)

  B(125) = RCT(91)*V(130)

  B(126) = RCT(92)*V(134)

  B(127) = RCT(92)*V(111)

  B(128) = RCT(93)*V(134)

  B(129) = RCT(93)*V(59)

  B(130) = RCT(94)*V(136)

  B(131) = RCT(94)*V(59)

  B(132) = RCT(95)*V(134)

  B(133) = RCT(95)*V(20)

  B(134) = RCT(96)*V(134)

  B(135) = RCT(96)*V(99)

  B(136) = RCT(97)*V(139)

  B(137) = RCT(97)*V(60)

  B(138) = RCT(98)*F(2)

  B(140) = RCT(99)

  B(141) = RCT(100)*V(134)

  B(142) = RCT(100)*V(16)

  B(143) = RCT(101)*V(139)

  B(144) = RCT(101)*V(102)

  B(145) = RCT(102)*V(135)

  B(146) = RCT(102)*V(102)

  B(147) = RCT(103)*V(140)

  B(148) = RCT(103)*V(102)

  B(149) = RCT(104)*V(134)

  B(150) = RCT(104)*V(48)

  B(151) = RCT(105)*V(134)

  B(152) = RCT(105)*V(113)

  B(153) = RCT(106)*V(136)

  B(154) = RCT(106)*V(113)

  B(155) = RCT(107)*V(133)

  B(156) = RCT(107)*V(113)

  B(157) = RCT(108)*V(139)

  B(158) = RCT(108)*V(90)

  B(159) = RCT(109)*V(135)

  B(160) = RCT(109)*V(90)

  B(161) = RCT(110)*V(134)

  B(162) = RCT(110)*V(72)

  B(163) = RCT(111)*V(134)

  B(164) = RCT(111)*V(114)

  B(165) = RCT(112)*V(133)

  B(166) = RCT(112)*V(114)

  B(167) = RCT(113)*V(139)

  B(168) = RCT(113)*V(137)

  B(169) = RCT(114)*V(138)

  B(170) = RCT(114)*V(137)

  B(171) = RCT(115)*V(137)

  B(172) = RCT(115)*V(135)

  B(173) = RCT(116)*V(140)

  B(174) = RCT(116)*V(137)

  B(175) = RCT(117)*V(134)

  B(176) = RCT(117)*V(68)

  B(177) = RCT(118)

  B(178) = RCT(119)*2*V(137)

  B(179) = RCT(120)*V(134)

  B(180) = RCT(120)*V(17)

  B(181) = RCT(121)*V(139)

  B(182) = RCT(121)*V(97)

  B(183) = RCT(122)*V(135)

  B(184) = RCT(122)*V(97)

  B(185) = RCT(123)*V(140)

  B(186) = RCT(123)*V(97)

  B(187) = RCT(124)*V(134)

  B(188) = RCT(124)*V(49)

  B(189) = RCT(125)*V(134)

  B(190) = RCT(125)*V(110)

  B(191) = RCT(126)*V(139)

  B(192) = RCT(126)*V(121)

  B(193) = RCT(127)*V(135)

  B(194) = RCT(127)*V(121)

  B(195) = RCT(128)*V(140)

  B(196) = RCT(128)*V(121)

  B(197) = RCT(129)*V(134)

  B(198) = RCT(129)*V(50)

  B(199) = RCT(130)*V(134)

  B(200) = RCT(130)*V(18)

  B(201) = RCT(131)*V(134)

  B(202) = RCT(131)*V(18)

  B(203) = RCT(132)*V(134)

  B(204) = RCT(132)*V(18)

  B(205) = RCT(133)*V(134)

  B(206) = RCT(133)*V(18)

  B(207) = RCT(134)*V(134)

  B(208) = RCT(134)*V(18)

  B(209) = RCT(135)*V(139)

  B(210) = RCT(135)*V(51)

  B(211) = RCT(136)*V(134)

  B(212) = RCT(136)*V(58)

  B(213) = RCT(137)*V(134)

  B(214) = RCT(137)*V(58)

  B(215) = RCT(138)*V(134)

  B(216) = RCT(138)*V(58)

  B(217) = RCT(139)*V(134)

  B(218) = RCT(139)*V(58)

  B(219) = RCT(140)*V(134)

  B(220) = RCT(140)*V(58)

  B(221) = RCT(141)*V(139)

  B(222) = RCT(141)*V(117)

  B(223) = RCT(142)*V(135)

  B(224) = RCT(142)*V(117)

  B(225) = RCT(143)*V(134)

  B(226) = RCT(143)*V(88)

  B(227) = RCT(144)*V(134)

  B(228) = RCT(144)*V(71)

  B(229) = RCT(145)*V(134)

  B(230) = RCT(145)*V(66)

  B(231) = RCT(146)*V(139)

  B(232) = RCT(146)*V(89)

  B(233) = RCT(147)*V(135)

  B(234) = RCT(147)*V(89)

  B(235) = RCT(148)*V(134)

  B(236) = RCT(148)*V(52)

  B(237) = RCT(149)*V(135)

  B(238) = RCT(149)*V(93)

  B(239) = RCT(150)*V(134)

  B(240) = RCT(150)*V(61)

  B(241) = RCT(151)*V(134)

  B(242) = RCT(151)*V(108)

  B(243) = RCT(152)*V(134)

  B(244) = RCT(152)*V(108)

  B(245) = RCT(153)*V(134)

  B(246) = RCT(153)*V(108)

  B(247) = RCT(154)*V(134)

  B(248) = RCT(154)*V(108)

  B(249) = RCT(155)*V(134)

  B(250) = RCT(155)*V(108)

  B(251) = RCT(156)*V(136)

  B(252) = RCT(156)*V(108)

  B(253) = RCT(157)*V(136)

  B(254) = RCT(157)*V(108)

  B(255) = RCT(158)*V(136)

  B(256) = RCT(158)*V(108)

  B(257) = RCT(159)*V(136)

  B(258) = RCT(159)*V(108)

  B(259) = RCT(160)*V(136)

  B(260) = RCT(160)*V(108)

  B(261) = RCT(161)*V(139)

  B(262) = RCT(161)*V(131)

  B(263) = RCT(162)*V(133)

  B(264) = RCT(162)*V(131)

  B(265) = RCT(163)*V(135)

  B(266) = RCT(163)*V(131)

  B(267) = RCT(164)*V(134)

  B(268) = RCT(164)*V(80)

  B(269) = RCT(165)*V(140)

  B(270) = RCT(165)*V(131)

  B(271) = RCT(166)*V(137)

  B(272) = RCT(166)*V(131)

  B(273) = RCT(167)*V(134)

  B(274) = RCT(167)*V(128)

  B(275) = RCT(168)*V(136)

  B(276) = RCT(168)*V(128)

  B(277) = RCT(169)*V(134)

  B(278) = RCT(169)*V(123)

  B(279) = RCT(170)*V(136)

  B(280) = RCT(170)*V(123)

  B(281) = RCT(171)*V(139)

  B(282) = RCT(171)*V(129)

  B(283) = RCT(172)*V(139)

  B(284) = RCT(172)*V(129)

  B(285) = RCT(173)*V(133)

  B(286) = RCT(173)*V(129)

  B(287) = RCT(174)*V(135)

  B(288) = RCT(174)*V(129)

  B(289) = RCT(175)*V(140)

  B(290) = RCT(175)*V(129)

  B(291) = RCT(176)*V(137)

  B(292) = RCT(176)*V(129)

  B(293) = RCT(177)*V(134)

  B(294) = RCT(177)*V(54)

  B(295) = RCT(178)*V(139)

  B(296) = RCT(178)*V(132)

  B(297) = RCT(179)*V(133)

  B(298) = RCT(179)*V(132)

  B(299) = RCT(180)*V(135)

  B(300) = RCT(180)*V(132)

  B(301) = RCT(181)*V(140)

  B(302) = RCT(181)*V(132)

  B(303) = RCT(182)*V(137)

  B(304) = RCT(182)*V(132)

  B(305) = RCT(183)*2*V(132)

  B(306) = RCT(184)*V(138)*F(2)

  B(307) = RCT(184)*V(132)*F(2)

  B(309) = RCT(185)*F(2)

  B(311) = RCT(186)*V(134)

  B(312) = RCT(186)*V(24)

  B(313) = RCT(187)*V(134)

  B(314) = RCT(187)*V(24)

  B(315) = RCT(188)*V(134)

  B(316) = RCT(188)*V(24)

  B(317) = RCT(189)*V(134)

  B(318) = RCT(189)*V(24)

  B(319) = RCT(190)*V(134)

  B(320) = RCT(190)*V(24)

  B(321) = RCT(191)*V(134)

  B(322) = RCT(191)*V(26)

  B(323) = RCT(192)*V(139)

  B(324) = RCT(192)*V(62)

  B(325) = RCT(193)*V(135)

  B(326) = RCT(193)*V(62)

  B(327) = RCT(194)*V(134)

  B(328) = RCT(194)*V(31)

  B(329) = RCT(195)*V(138)

  B(330) = RCT(195)*V(63)

  B(331) = RCT(196)*V(136)

  B(332) = RCT(196)*V(63)

  B(333) = RCT(197)*V(139)

  B(334) = RCT(197)*V(85)

  B(335) = RCT(198)*V(135)

  B(336) = RCT(198)*V(85)

  B(337) = RCT(199)*V(134)

  B(338) = RCT(199)*V(32)

  B(339) = RCT(200)*V(139)

  B(340) = RCT(200)*V(70)

  B(341) = RCT(201)*V(135)

  B(342) = RCT(201)*V(70)

  B(343) = RCT(202)*V(134)

  B(344) = RCT(202)*V(33)

  B(345) = RCT(203)*V(138)

  B(346) = RCT(203)*V(91)

  B(347) = RCT(204)*V(139)

  B(348) = RCT(204)*V(91)

  B(349) = RCT(205)*V(135)

  B(350) = RCT(205)*V(91)

  B(351) = RCT(206)*V(134)

  B(352) = RCT(206)*V(27)

  B(353) = RCT(207)*V(134)

  B(354) = RCT(207)*V(27)

  B(355) = RCT(208)*V(134)

  B(356) = RCT(208)*V(27)

  B(357) = RCT(209)*V(134)

  B(358) = RCT(209)*V(27)

  B(359) = RCT(210)*V(134)

  B(360) = RCT(210)*V(27)

  B(361) = RCT(211)*V(134)

  B(362) = RCT(211)*V(28)

  B(363) = RCT(212)*V(135)

  B(364) = RCT(212)*V(64)

  B(365) = RCT(213)*V(134)

  B(366) = RCT(213)*V(36)

  B(367) = RCT(214)*V(139)

  B(368) = RCT(214)*V(64)

  B(369) = RCT(215)*V(134)

  B(370) = RCT(215)*V(37)

  B(371) = RCT(216)*V(138)

  B(372) = RCT(216)*V(78)

  B(373) = RCT(217)

  B(374) = RCT(218)*V(139)

  B(375) = RCT(218)*V(78)

  B(376) = RCT(219)*V(135)

  B(377) = RCT(219)*V(78)

  B(378) = RCT(220)*V(139)

  B(379) = RCT(220)*V(93)

  B(380) = RCT(221)*V(135)

  B(381) = RCT(221)*V(95)

  B(382) = RCT(222)*V(139)

  B(383) = RCT(222)*V(95)

  B(384) = RCT(223)*V(135)

  B(385) = RCT(223)*V(98)

  B(386) = RCT(224)*V(139)

  B(387) = RCT(224)*V(98)

  B(388) = RCT(225)*V(138)

  B(389) = RCT(225)*V(95)

  B(390) = RCT(226)*V(138)

  B(391) = RCT(226)*V(98)

  B(392) = RCT(227)*V(134)

  B(393) = RCT(227)*V(41)

  B(394) = RCT(228)*V(134)

  B(395) = RCT(228)*V(41)

  B(396) = RCT(229)*V(134)

  B(397) = RCT(229)*V(41)

  B(398) = RCT(230)*V(134)

  B(399) = RCT(230)*V(41)

  B(400) = RCT(231)*V(134)

  B(401) = RCT(231)*V(41)

  B(402) = RCT(232)*V(134)

  B(403) = RCT(232)*V(42)

  B(404) = RCT(233)*V(139)

  B(405) = RCT(233)*V(74)

  B(406) = RCT(234)*V(135)

  B(407) = RCT(234)*V(74)

  B(408) = RCT(235)*V(134)

  B(409) = RCT(235)*V(38)

  B(410) = RCT(236)*V(135)

  B(411) = RCT(236)*V(94)

  B(412) = RCT(237)*V(134)

  B(413) = RCT(237)*V(39)

  B(414) = RCT(238)*V(139)

  B(415) = RCT(238)*V(94)

  B(416) = RCT(239)*V(134)

  B(417) = RCT(239)*V(87)

  B(418) = RCT(240)*V(134)

  B(419) = RCT(240)*V(87)

  B(420) = RCT(241)*V(134)

  B(421) = RCT(241)*V(87)

  B(422) = RCT(242)*V(134)

  B(423) = RCT(242)*V(87)

  B(424) = RCT(243)*V(134)

  B(425) = RCT(243)*V(87)

  B(426) = RCT(244)*V(134)

  B(427) = RCT(244)*V(86)

  B(428) = RCT(245)*V(134)

  B(429) = RCT(245)*V(86)

  B(430) = RCT(246)*V(134)

  B(431) = RCT(246)*V(86)

  B(432) = RCT(247)*V(134)

  B(433) = RCT(247)*V(86)

  B(434) = RCT(248)*V(134)

  B(435) = RCT(248)*V(86)

  B(436) = RCT(249)*V(134)

  B(437) = RCT(249)*V(104)

  B(438) = RCT(250)*V(134)

  B(439) = RCT(250)*V(104)

  B(440) = RCT(251)*V(134)

  B(441) = RCT(251)*V(104)

  B(442) = RCT(252)*V(134)

  B(443) = RCT(252)*V(104)

  B(444) = RCT(253)*V(134)

  B(445) = RCT(253)*V(104)

  B(446) = RCT(254)*V(134)

  B(447) = RCT(254)*V(107)

  B(448) = RCT(255)*V(134)

  B(449) = RCT(255)*V(109)

  B(450) = RCT(256)*V(136)

  B(451) = RCT(256)*V(87)

  B(452) = RCT(257)*V(136)

  B(453) = RCT(257)*V(87)

  B(454) = RCT(258)*V(136)

  B(455) = RCT(258)*V(87)

  B(456) = RCT(259)*V(136)

  B(457) = RCT(259)*V(87)

  B(458) = RCT(260)*V(136)

  B(459) = RCT(260)*V(87)

  B(460) = RCT(261)*V(136)

  B(461) = RCT(261)*V(86)

  B(462) = RCT(262)*V(136)

  B(463) = RCT(262)*V(86)

  B(464) = RCT(263)*V(136)

  B(465) = RCT(263)*V(86)

  B(466) = RCT(264)*V(136)

  B(467) = RCT(264)*V(86)

  B(468) = RCT(265)*V(136)

  B(469) = RCT(265)*V(86)

  B(470) = RCT(266)*V(136)

  B(471) = RCT(266)*V(104)

  B(472) = RCT(267)*V(136)

  B(473) = RCT(267)*V(104)

  B(474) = RCT(268)*V(136)

  B(475) = RCT(268)*V(104)

  B(476) = RCT(269)*V(136)

  B(477) = RCT(269)*V(104)

  B(478) = RCT(270)*V(136)

  B(479) = RCT(270)*V(104)

  B(480) = RCT(271)*V(136)

  B(481) = RCT(271)*V(107)

  B(482) = RCT(272)*V(136)

  B(483) = RCT(272)*V(109)

  B(484) = RCT(273)*V(133)

  B(485) = RCT(273)*V(104)

  B(486) = RCT(274)*V(133)

  B(487) = RCT(274)*V(107)

  B(488) = RCT(275)*V(133)

  B(489) = RCT(275)*V(109)

  B(490) = RCT(276)*V(139)

  B(491) = RCT(276)*V(119)

  B(492) = RCT(277)*V(135)

  B(493) = RCT(277)*V(119)

  B(494) = RCT(278)*V(140)

  B(495) = RCT(278)*V(119)

  B(496) = RCT(279)*V(134)

  B(497) = RCT(279)*V(69)

  B(498) = RCT(280)*V(134)

  B(499) = RCT(280)*V(83)

  B(500) = RCT(281)*V(134)

  B(501) = RCT(281)*V(124)

  B(502) = RCT(282)*V(133)

  B(503) = RCT(282)*V(124)

  B(504) = RCT(283)*V(139)

  B(505) = RCT(283)*V(118)

  B(506) = RCT(284)*V(135)

  B(507) = RCT(284)*V(118)

  B(508) = RCT(285)*V(140)

  B(509) = RCT(285)*V(118)

  B(510) = RCT(286)*V(134)

  B(511) = RCT(286)*V(115)

  B(512) = RCT(287)*V(139)

  B(513) = RCT(287)*V(122)

  B(514) = RCT(288)*V(135)

  B(515) = RCT(288)*V(122)

  B(516) = RCT(289)*V(140)

  B(517) = RCT(289)*V(122)

  B(518) = RCT(290)*V(133)

  B(519) = RCT(290)*V(122)

  B(520) = RCT(291)*V(134)

  B(521) = RCT(291)*V(67)

  B(522) = RCT(292)*V(133)

  B(523) = RCT(292)*V(108)

  B(524) = RCT(293)*V(139)

  B(525) = RCT(293)*V(103)

  B(526) = RCT(294)*V(133)

  B(527) = RCT(294)*V(103)

  B(528) = RCT(295)*V(135)

  B(529) = RCT(295)*V(103)

  B(530) = RCT(296)*V(134)

  B(531) = RCT(296)*V(126)

  B(532) = RCT(297)*V(133)

  B(533) = RCT(297)*V(126)

  B(534) = RCT(298)*V(134)

  B(535) = RCT(298)*V(120)

  B(536) = RCT(299)*V(133)

  B(537) = RCT(299)*V(120)

  B(538) = RCT(300)*V(134)

  B(539) = RCT(300)*V(57)

  B(540) = RCT(301)*V(139)

  B(541) = RCT(301)*V(127)

  B(542) = RCT(302)*V(133)

  B(543) = RCT(302)*V(127)

  B(544) = RCT(303)*V(135)

  B(545) = RCT(303)*V(127)

  B(546) = RCT(304)*V(140)

  B(547) = RCT(304)*V(127)

  B(548) = RCT(305)*V(137)

  B(549) = RCT(305)*V(127)

  B(550) = RCT(306)*V(134)

  B(551) = RCT(306)*V(29)

  B(552) = RCT(307)*V(134)

  B(553) = RCT(307)*V(29)

  B(554) = RCT(308)*V(134)

  B(555) = RCT(308)*V(100)

  B(556) = RCT(309)*V(134)

  B(557) = RCT(309)*V(45)

  B(558) = RCT(310)*V(134)

  B(559) = RCT(310)*V(81)

  B(560) = RCT(311)*V(134)

  B(561) = RCT(311)*V(73)

  B(562) = RCT(312)*V(134)

  B(563) = RCT(312)*V(125)

  B(564) = RCT(313)*V(134)

  B(565) = RCT(313)*V(112)

  B(566) = RCT(314)*V(134)

  B(567) = RCT(314)*V(43)

  B(568) = RCT(315)*V(134)

  B(569) = RCT(315)*V(43)

  B(570) = RCT(316)*V(133)

  B(571) = RCT(316)*V(43)

  B(572) = RCT(317)*V(134)

  B(573) = RCT(317)*V(15)

  B(574) = RCT(318)

  B(575) = RCT(319)*2*V(102)

  B(576) = RCT(320)*V(134)

  B(577) = RCT(320)*V(92)

  B(578) = RCT(321)*V(139)

  B(579) = RCT(321)*V(116)

  B(580) = RCT(322)*V(140)

  B(581) = RCT(322)*V(116)

  B(582) = RCT(323)*V(134)

  B(583) = RCT(323)*V(106)

  B(584) = RCT(324)*V(139)

  B(585) = RCT(324)*V(76)

  B(586) = RCT(325)*V(135)

  B(587) = RCT(325)*V(76)

  B(588) = RCT(326)*V(135)

  B(589) = RCT(326)*V(116)

  B(590) = RCT(327)*V(134)

  B(591) = RCT(327)*V(79)

  B(592) = RCT(328)*V(136)

  B(593) = RCT(328)*V(92)

  B(594) = RCT(329)*V(133)

  B(595) = RCT(329)*V(92)

  B(596) = RCT(330)*V(135)

  B(597) = RCT(330)*V(105)

  B(598) = RCT(331)*V(139)

  B(599) = RCT(331)*V(105)

  B(600) = RCT(332)*V(133)

  B(601) = RCT(332)*V(105)

  B(602) = RCT(333)*V(134)

  B(603) = RCT(333)*V(23)

  B(604) = 4.5e-13*V(134)

  B(605) = 4.5e-13*V(82)

  B(606) = RCT(335)*V(135)

  B(607) = RCT(335)*V(130)

  B(608) = RCT(336)

  B(609) = RCT(337)*V(139)

  B(610) = RCT(337)*V(65)

  B(611) = RCT(338)*V(135)

  B(612) = RCT(338)*V(65)

  B(613) = RCT(339)*V(134)

  B(614) = RCT(339)*V(9)

  B(615) = RCT(340)*V(134)

  B(616) = RCT(340)*V(8)

  B(617) = RCT(341)*V(134)

  B(618) = RCT(341)*V(7)

  B(619) = RCT(342)*V(134)

  B(620) = RCT(342)*V(6)

  B(621) = RCT(343)*V(134)

  B(622) = RCT(343)*V(14)

  B(623) = RCT(344)*V(134)

  B(624) = RCT(344)*V(13)

  B(625) = RCT(345)*V(134)

  B(626) = RCT(345)*V(12)

  B(627) = RCT(346)*V(134)

  B(628) = RCT(346)*V(11)



  JVS(1) = -B(49)

  JVS(2) = 0

  JVS(3) = B(132)

  JVS(4) = B(133)

  JVS(5) = 0

  JVS(6) = B(209)

  JVS(7) = B(136)

  JVS(8) = B(323)

  JVS(9) = B(367)

  JVS(10) = B(609)

  JVS(11) = B(339)

  JVS(12) = B(404)

  JVS(13) = B(584)

  JVS(14) = B(374)

  JVS(15) = B(333)

  JVS(16) = B(231)

  JVS(17) = B(157)

  JVS(18) = B(347)

  JVS(19) = B(378)

  JVS(20) = B(414)

  JVS(21) = B(382)

  JVS(22) = B(181)

  JVS(23) = B(386)

  JVS(24) = B(143)

  JVS(25) = B(524)

  JVS(26) = B(598)

  JVS(27) = B(578)

  JVS(28) = B(221)

  JVS(29) = B(504)

  JVS(30) = B(490)

  JVS(31) = B(191)

  JVS(32) = B(512)

  JVS(33) = B(540)

  JVS(34) = B(281)+B(283)

  JVS(35) = B(261)

  JVS(36) = B(295)

  JVS(37) = B(167)

  JVS(38) = B(114)+B(137)+B(144)+B(158)+B(168)+B(182)+B(192)+B(210)+B(222)+B(232)+B(262)+B(282)+B(284)+B(296)+B(324)&
              &+B(334)+B(340)+B(348)+B(368)+B(375)+B(379)+B(383)+B(387)+B(405)+B(415)+B(491)+B(505)+B(513)+B(525)+B(541)&
              &+B(579)+B(585)+B(599)+B(610)

  JVS(39) = B(115)

  JVS(40) = 0

  JVS(41) = B(325)

  JVS(42) = B(363)

  JVS(43) = B(611)

  JVS(44) = B(341)

  JVS(45) = B(406)

  JVS(46) = B(586)

  JVS(47) = B(376)

  JVS(48) = B(335)

  JVS(49) = B(233)

  JVS(50) = B(159)

  JVS(51) = B(349)

  JVS(52) = B(237)

  JVS(53) = B(410)

  JVS(54) = B(380)

  JVS(55) = B(183)+B(185)

  JVS(56) = B(384)

  JVS(57) = B(145)+B(147)+B(575)

  JVS(58) = B(528)

  JVS(59) = B(596)

  JVS(60) = B(580)+B(588)

  JVS(61) = B(223)

  JVS(62) = B(506)+B(508)

  JVS(63) = B(492)+B(494)

  JVS(64) = B(193)+B(195)

  JVS(65) = B(514)+B(516)

  JVS(66) = B(544)+B(546)+B(548)

  JVS(67) = B(287)+B(289)+B(291)

  JVS(68) = B(265)+B(269)+B(271)

  JVS(69) = B(299)+B(301)+B(303)+B(305)

  JVS(70) = B(118)+B(146)+B(160)+B(171)+B(184)+B(194)+B(224)+B(234)+B(238)+B(266)+B(288)+B(300)+B(326)+B(336)+B(342)&
              &+B(350)+B(364)+B(377)+B(381)+B(385)+B(407)+B(411)+B(493)+B(507)+B(515)+B(529)+B(545)+B(587)+B(589)+B(597)&
              &+B(612)

  JVS(71) = B(172)+B(173)+B(178)+B(272)+B(292)+B(304)+B(549)

  JVS(72) = B(116)+B(117)+B(119)+B(148)+B(174)+B(186)+B(196)+B(270)+B(290)+B(302)+B(495)+B(509)+B(517)+B(547)+B(581)

  JVS(73) = 0

  JVS(74) = 1.075*B(619)

  JVS(75) = 1.075*B(620)

  JVS(76) = -B(619)

  JVS(77) = 1.075*B(617)

  JVS(78) = B(207)

  JVS(79) = B(319)

  JVS(80) = B(359)

  JVS(81) = B(400)

  JVS(82) = B(219)

  JVS(83) = B(208)+B(220)+B(320)+B(360)+B(401)+1.075*B(618)-B(620)

  JVS(84) = -B(617)

  JVS(85) = 1.075*B(615)

  JVS(86) = B(205)

  JVS(87) = B(317)

  JVS(88) = B(357)

  JVS(89) = B(398)

  JVS(90) = B(217)

  JVS(91) = B(206)+B(218)+B(318)+B(358)+B(399)+1.075*B(616)-B(618)

  JVS(92) = -B(615)

  JVS(93) = 1.075*B(613)

  JVS(94) = B(203)

  JVS(95) = B(315)

  JVS(96) = B(355)

  JVS(97) = B(396)

  JVS(98) = B(215)

  JVS(99) = B(204)+B(216)+B(316)+B(356)+B(397)+1.075*B(614)-B(616)

  JVS(100) = -B(613)

  JVS(101) = B(201)

  JVS(102) = B(313)

  JVS(103) = B(353)

  JVS(104) = B(394)

  JVS(105) = B(213)

  JVS(106) = B(202)+B(214)+B(314)+B(354)+B(395)-B(614)

  JVS(107) = 0

  JVS(108) = 1.075*B(627)

  JVS(109) = 1.075*B(628)

  JVS(110) = -B(627)

  JVS(111) = 1.075*B(625)

  JVS(112) = B(434)+B(468)

  JVS(113) = B(424)+B(458)

  JVS(114) = B(444)+B(478)

  JVS(115) = B(249)+B(259)

  JVS(116) = B(250)+B(425)+B(435)+B(445)+1.075*B(626)-B(628)

  JVS(117) = B(260)+B(459)+B(469)+B(479)

  JVS(118) = -B(625)

  JVS(119) = 1.075*B(623)

  JVS(120) = B(432)+B(466)

  JVS(121) = B(422)+B(456)

  JVS(122) = B(442)+B(476)

  JVS(123) = B(247)+B(257)

  JVS(124) = B(248)+B(423)+B(433)+B(443)+1.075*B(624)-B(626)

  JVS(125) = B(258)+B(457)+B(467)+B(477)

  JVS(126) = -B(623)

  JVS(127) = 1.075*B(621)

  JVS(128) = B(430)+B(464)

  JVS(129) = B(420)+B(454)

  JVS(130) = B(440)+B(474)

  JVS(131) = B(245)+B(255)

  JVS(132) = B(246)+B(421)+B(431)+B(441)+1.075*B(622)-B(624)

  JVS(133) = B(256)+B(455)+B(465)+B(475)

  JVS(134) = -B(621)

  JVS(135) = B(428)+B(462)

  JVS(136) = B(418)+B(452)

  JVS(137) = B(438)+B(472)

  JVS(138) = B(243)+B(253)

  JVS(139) = B(244)+B(419)+B(429)+B(439)-B(622)

  JVS(140) = B(254)+B(453)+B(463)+B(473)

  JVS(141) = -B(572)

  JVS(142) = -B(573)

  JVS(143) = -B(141)

  JVS(144) = -B(142)

  JVS(145) = -B(179)

  JVS(146) = -B(180)

  JVS(147) = -B(199)

  JVS(148) = -B(200)

  JVS(149) = -B(4)-B(78)-B(80)

  JVS(150) = -B(79)-B(81)

  JVS(151) = -B(132)

  JVS(152) = B(566)+0.5*B(568)+B(570)

  JVS(153) = B(571)

  JVS(154) = -B(133)+B(567)+0.5*B(569)

  JVS(155) = -B(373)

  JVS(156) = B(371)

  JVS(157) = B(372)

  JVS(158) = -B(13)-B(72)

  JVS(159) = -B(73)+B(77)

  JVS(160) = B(70)+0.5*B(574)

  JVS(161) = -B(602)

  JVS(162) = -B(603)

  JVS(163) = -B(311)

  JVS(164) = -B(312)

  JVS(165) = 0.12*B(311)

  JVS(166) = -B(38)

  JVS(167) = 0.12*B(312)

  JVS(168) = 0.53*B(311)

  JVS(169) = -B(321)

  JVS(170) = 0.53*B(312)-B(322)

  JVS(171) = -B(351)

  JVS(172) = -B(352)

  JVS(173) = 0.18*B(351)

  JVS(174) = -B(361)

  JVS(175) = 0.18*B(352)-B(362)

  JVS(176) = -B(26)-B(550)-B(552)

  JVS(177) = B(544)

  JVS(178) = -B(551)-B(553)

  JVS(179) = B(545)

  JVS(180) = -B(138)-B(140)

  JVS(181) = B(136)

  JVS(182) = B(137)

  JVS(183) = -B(327)

  JVS(184) = B(325)

  JVS(185) = -B(328)

  JVS(186) = B(326)

  JVS(187) = -B(337)

  JVS(188) = B(335)

  JVS(189) = -B(338)

  JVS(190) = B(336)

  JVS(191) = -B(343)

  JVS(192) = B(341)

  JVS(193) = -B(344)

  JVS(194) = B(342)

  JVS(195) = -B(40)

  JVS(196) = 0.2*B(378)

  JVS(197) = 0.2*B(414)

  JVS(198) = 0.2*B(379)+0.2*B(415)

  JVS(199) = 0.1*B(351)

  JVS(200) = -B(39)

  JVS(201) = 0.23*B(392)

  JVS(202) = 0.1*B(352)+0.23*B(393)

  JVS(203) = -B(365)

  JVS(204) = B(363)

  JVS(205) = -B(366)

  JVS(206) = B(364)

  JVS(207) = -B(369)

  JVS(208) = B(367)

  JVS(209) = -B(370)

  JVS(210) = B(368)

  JVS(211) = -B(408)

  JVS(212) = B(406)

  JVS(213) = -B(409)

  JVS(214) = B(407)

  JVS(215) = -B(412)

  JVS(216) = B(410)

  JVS(217) = -B(413)

  JVS(218) = B(411)

  JVS(219) = -B(6)-B(94)-B(106)

  JVS(220) = B(92)

  JVS(221) = B(93)

  JVS(222) = -B(392)

  JVS(223) = -B(393)

  JVS(224) = 0.15*B(392)

  JVS(225) = -B(402)

  JVS(226) = 0.15*B(393)-B(403)

  JVS(227) = -B(566)-B(568)-B(570)

  JVS(228) = -B(571)

  JVS(229) = -B(567)-B(569)

  JVS(230) = -B(42)

  JVS(231) = 0.21*B(414)

  JVS(232) = 0.21*B(415)

  JVS(233) = -B(556)

  JVS(234) = 0.2*B(147)+0.4*B(575)

  JVS(235) = -B(557)

  JVS(236) = 0.2*B(148)

  JVS(237) = -B(41)

  JVS(238) = 0.2*B(378)

  JVS(239) = 0.15*B(414)

  JVS(240) = 0.2*B(379)+0.15*B(415)

  JVS(241) = B(38)

  JVS(242) = -B(37)

  JVS(243) = 0.5*B(339)

  JVS(244) = 0.2*B(378)

  JVS(245) = 0.06*B(414)

  JVS(246) = 0

  JVS(247) = 0.5*B(340)+0.2*B(379)+0.06*B(415)

  JVS(248) = -B(21)-B(149)

  JVS(249) = B(145)

  JVS(250) = -B(150)

  JVS(251) = B(146)

  JVS(252) = -B(22)-B(187)

  JVS(253) = B(183)

  JVS(254) = -B(188)

  JVS(255) = B(184)

  JVS(256) = -B(23)-B(197)

  JVS(257) = B(193)

  JVS(258) = -B(198)

  JVS(259) = B(194)

  JVS(260) = B(199)

  JVS(261) = -B(209)

  JVS(262) = B(200)

  JVS(263) = -B(210)

  JVS(264) = -B(35)-B(235)

  JVS(265) = B(233)

  JVS(266) = -B(236)

  JVS(267) = B(234)

  JVS(268) = -B(10)-B(120)

  JVS(269) = -B(121)

  JVS(270) = B(118)

  JVS(271) = B(119)

  JVS(272) = -B(293)

  JVS(273) = B(287)

  JVS(274) = -B(294)

  JVS(275) = B(288)

  JVS(276) = -B(58)-B(60)

  JVS(277) = 0.05*B(112)

  JVS(278) = -B(59)+0.05*B(113)

  JVS(279) = B(12)

  JVS(280) = -B(61)

  JVS(281) = -B(9)-B(103)-B(105)

  JVS(282) = -B(104)

  JVS(283) = B(101)

  JVS(284) = B(102)

  JVS(285) = -B(538)

  JVS(286) = B(534)+B(536)

  JVS(287) = 0.33*B(261)+0.36*B(263)+0.3*B(269)+0.4*B(271)

  JVS(288) = 0.36*B(264)+B(537)

  JVS(289) = B(535)-B(539)

  JVS(290) = 0.4*B(272)

  JVS(291) = 0.33*B(262)

  JVS(292) = 0.3*B(270)

  JVS(293) = -B(211)

  JVS(294) = 0.05*B(460)

  JVS(295) = 0.05*B(450)

  JVS(296) = 0.05*B(470)

  JVS(297) = 0.05*B(480)

  JVS(298) = 0.05*B(482)

  JVS(299) = -B(212)

  JVS(300) = 0.05*B(451)+0.05*B(461)+0.05*B(471)+0.05*B(481)+0.05*B(483)

  JVS(301) = -B(128)-B(130)

  JVS(302) = -B(129)

  JVS(303) = -B(131)

  JVS(304) = 0.75*B(128)

  JVS(305) = -B(136)

  JVS(306) = 0.75*B(129)

  JVS(307) = 0

  JVS(308) = -B(137)

  JVS(309) = -B(36)-B(239)

  JVS(310) = B(237)

  JVS(311) = -B(240)

  JVS(312) = B(238)

  JVS(313) = 0.14*B(321)

  JVS(314) = 0.2*B(361)

  JVS(315) = B(327)

  JVS(316) = -B(323)-B(325)

  JVS(317) = 0.14*B(322)+B(328)+0.2*B(362)

  JVS(318) = -B(326)

  JVS(319) = -B(324)

  JVS(320) = 0.06*B(321)

  JVS(321) = 0.07*B(361)

  JVS(322) = 0.07*B(402)

  JVS(323) = -B(329)-B(331)

  JVS(324) = B(333)

  JVS(325) = 0.06*B(322)+0.07*B(362)+0.07*B(403)

  JVS(326) = -B(332)

  JVS(327) = -B(330)

  JVS(328) = B(334)

  JVS(329) = 0.07*B(351)

  JVS(330) = B(365)

  JVS(331) = 0.06*B(392)

  JVS(332) = -B(363)-B(367)

  JVS(333) = 0.07*B(352)+B(366)+0.06*B(393)

  JVS(334) = -B(364)

  JVS(335) = -B(368)

  JVS(336) = -B(608)-B(609)-B(611)

  JVS(337) = B(606)

  JVS(338) = B(607)-B(612)

  JVS(339) = -B(610)

  JVS(340) = -B(31)-B(229)

  JVS(341) = 0.8*B(34)

  JVS(342) = 0.75*B(221)

  JVS(343) = -B(230)

  JVS(344) = 0.75*B(222)

  JVS(345) = -B(520)

  JVS(346) = 0.2*B(251)

  JVS(347) = 0.25*B(153)

  JVS(348) = 0.25*B(299)

  JVS(349) = -B(521)

  JVS(350) = 0.25*B(171)+0.25*B(300)

  JVS(351) = 0.25*B(154)+0.2*B(252)

  JVS(352) = 0.25*B(172)+0.1*B(173)

  JVS(353) = 0.1*B(174)

  JVS(354) = -B(16)-B(175)

  JVS(355) = 0.75*B(299)

  JVS(356) = -B(176)

  JVS(357) = 0.75*B(171)+0.75*B(300)

  JVS(358) = 0.75*B(172)

  JVS(359) = -B(47)-B(496)

  JVS(360) = B(492)

  JVS(361) = -B(497)

  JVS(362) = B(493)

  JVS(363) = 0.35*B(311)

  JVS(364) = B(343)

  JVS(365) = -B(339)-B(341)

  JVS(366) = 0.35*B(312)+B(344)

  JVS(367) = -B(342)

  JVS(368) = -B(340)

  JVS(369) = -B(227)

  JVS(370) = 0.75*B(598)+0.75*B(600)

  JVS(371) = B(155)

  JVS(372) = 0.1*B(221)

  JVS(373) = B(156)+0.75*B(601)

  JVS(374) = -B(228)

  JVS(375) = 0.1*B(222)+0.75*B(599)

  JVS(376) = -B(15)-B(161)

  JVS(377) = B(159)

  JVS(378) = -B(162)

  JVS(379) = B(160)

  JVS(380) = -B(17)-B(177)-B(560)

  JVS(381) = -B(561)

  JVS(382) = B(169)

  JVS(383) = B(170)

  JVS(384) = B(408)

  JVS(385) = 0.3*B(402)

  JVS(386) = -B(404)-B(406)

  JVS(387) = 0.3*B(403)+B(409)

  JVS(388) = -B(407)

  JVS(389) = -B(405)

  JVS(390) = 2*B(106)

  JVS(391) = B(570)

  JVS(392) = -B(7)-B(97)

  JVS(393) = B(165)

  JVS(394) = B(532)

  JVS(395) = B(122)

  JVS(396) = B(108)+B(123)+B(166)+B(533)+B(571)

  JVS(397) = B(95)-B(98)

  JVS(398) = B(96)+0.5*B(109)

  JVS(399) = -B(584)-B(586)

  JVS(400) = B(582)

  JVS(401) = B(583)

  JVS(402) = -B(587)

  JVS(403) = -B(585)

  JVS(404) = -B(110)-B(112)

  JVS(405) = -B(113)

  JVS(406) = 0.08*B(153)

  JVS(407) = -B(111)

  JVS(408) = 0.08*B(154)

  JVS(409) = B(373)

  JVS(410) = B(369)

  JVS(411) = 0

  JVS(412) = -B(371)-B(374)-B(376)

  JVS(413) = B(370)

  JVS(414) = -B(377)

  JVS(415) = -B(372)

  JVS(416) = -B(375)

  JVS(417) = -B(43)-B(590)

  JVS(418) = B(588)

  JVS(419) = -B(591)

  JVS(420) = B(589)

  JVS(421) = -B(28)-B(267)

  JVS(422) = B(265)

  JVS(423) = -B(268)

  JVS(424) = B(266)

  JVS(425) = -B(18)-B(309)-B(558)

  JVS(426) = B(306)

  JVS(427) = -B(559)

  JVS(428) = B(307)

  JVS(429) = 0.35*B(602)

  JVS(430) = 0.5*B(130)

  JVS(431) = B(609)+B(611)

  JVS(432) = -B(604)

  JVS(433) = 0.05*B(460)

  JVS(434) = 0.05*B(450)

  JVS(435) = 0.25*B(592)

  JVS(436) = 0.05*B(470)

  JVS(437) = 0.05*B(480)

  JVS(438) = 0.05*B(482)

  JVS(439) = 0

  JVS(440) = 0.35*B(603)-B(605)

  JVS(441) = B(612)

  JVS(442) = 0.5*B(131)+0.05*B(451)+0.05*B(461)+0.05*B(471)+0.05*B(481)+0.05*B(483)+0.25*B(593)

  JVS(443) = B(610)

  JVS(444) = -B(48)-B(498)

  JVS(445) = B(506)

  JVS(446) = -B(499)

  JVS(447) = B(507)

  JVS(448) = B(4)-B(78)-B(80)

  JVS(449) = -B(58)

  JVS(450) = -B(112)

  JVS(451) = -B(54)-B(56)-B(59)-B(79)-B(81)-B(113)

  JVS(452) = 0

  JVS(453) = 0

  JVS(454) = 0

  JVS(455) = B(2)

  JVS(456) = B(337)

  JVS(457) = B(331)

  JVS(458) = B(374)+0.4*B(376)

  JVS(459) = -B(333)-B(335)

  JVS(460) = B(338)

  JVS(461) = -B(336)+0.4*B(377)

  JVS(462) = B(332)

  JVS(463) = 0

  JVS(464) = -B(334)+B(375)

  JVS(465) = -B(426)-B(460)

  JVS(466) = -B(427)

  JVS(467) = -B(461)

  JVS(468) = -B(416)-B(450)

  JVS(469) = -B(417)

  JVS(470) = -B(451)

  JVS(471) = -B(34)-B(225)

  JVS(472) = B(223)

  JVS(473) = -B(226)

  JVS(474) = B(224)

  JVS(475) = B(235)

  JVS(476) = B(229)

  JVS(477) = 0

  JVS(478) = -B(231)-B(233)

  JVS(479) = 0

  JVS(480) = B(230)+B(236)

  JVS(481) = -B(234)

  JVS(482) = -B(232)

  JVS(483) = 0.5*B(161)

  JVS(484) = -B(157)-B(159)

  JVS(485) = B(151)

  JVS(486) = B(152)+0.5*B(162)

  JVS(487) = -B(160)

  JVS(488) = -B(158)

  JVS(489) = 0.6*B(37)

  JVS(490) = 0

  JVS(491) = -B(345)-B(347)-B(349)

  JVS(492) = 0

  JVS(493) = 0

  JVS(494) = 0

  JVS(495) = -B(350)

  JVS(496) = -B(346)

  JVS(497) = -B(348)

  JVS(498) = -B(576)-B(592)-B(594)

  JVS(499) = -B(595)

  JVS(500) = -B(577)

  JVS(501) = -B(593)

  JVS(502) = 0.65*B(351)

  JVS(503) = B(239)

  JVS(504) = -B(237)-B(378)

  JVS(505) = B(240)+0.65*B(352)

  JVS(506) = -B(238)

  JVS(507) = -B(379)

  JVS(508) = B(412)

  JVS(509) = 0.56*B(392)

  JVS(510) = -B(410)-B(414)

  JVS(511) = 0.56*B(393)+B(413)

  JVS(512) = -B(411)

  JVS(513) = -B(415)

  JVS(514) = 0.6*B(40)

  JVS(515) = 0

  JVS(516) = 0

  JVS(517) = -B(380)-B(382)-B(388)

  JVS(518) = 0

  JVS(519) = -B(381)

  JVS(520) = -B(389)

  JVS(521) = -B(383)

  JVS(522) = B(54)

  JVS(523) = -B(50)-B(52)-B(62)-B(64)-B(86)

  JVS(524) = 0

  JVS(525) = 0

  JVS(526) = -B(63)+B(76)

  JVS(527) = -B(65)

  JVS(528) = B(3)-B(53)

  JVS(529) = B(5)-B(87)

  JVS(530) = B(179)

  JVS(531) = B(187)

  JVS(532) = -B(181)-B(183)-B(185)

  JVS(533) = B(180)+B(188)

  JVS(534) = -B(184)

  JVS(535) = -B(182)

  JVS(536) = -B(186)

  JVS(537) = 0.6*B(41)

  JVS(538) = 0

  JVS(539) = 0

  JVS(540) = -B(384)-B(386)-B(390)

  JVS(541) = 0

  JVS(542) = -B(385)

  JVS(543) = -B(391)

  JVS(544) = -B(387)

  JVS(545) = 0.65*B(602)

  JVS(546) = 0.45*B(36)

  JVS(547) = 0.7*B(323)

  JVS(548) = B(339)

  JVS(549) = 0.17*B(404)

  JVS(550) = 0.4*B(347)+0.16*B(349)

  JVS(551) = 0.6*B(378)

  JVS(552) = 0.34*B(414)

  JVS(553) = 0.07*B(384)+0.17*B(386)

  JVS(554) = -B(33)-B(134)

  JVS(555) = 0.13*B(32)

  JVS(556) = 0.2*B(564)

  JVS(557) = 0.25*B(540)+0.25*B(542)+0.1*B(546)+0.25*B(548)

  JVS(558) = 0.02*B(261)+0.02*B(263)

  JVS(559) = 0.02*B(264)+0.25*B(543)

  JVS(560) = -B(135)+0.2*B(565)+0.65*B(603)

  JVS(561) = 0.16*B(350)+0.07*B(385)

  JVS(562) = 0.25*B(549)

  JVS(563) = 0

  JVS(564) = 0.02*B(262)+0.7*B(324)+B(340)+0.4*B(348)+0.6*B(379)+0.17*B(387)+0.17*B(405)+0.34*B(415)+0.25*B(541)

  JVS(565) = 0.1*B(547)

  JVS(566) = -B(554)

  JVS(567) = 0.3*B(147)

  JVS(568) = 0.25*B(580)

  JVS(569) = 0.25*B(508)

  JVS(570) = 0.25*B(494)

  JVS(571) = 0.5*B(195)

  JVS(572) = 0.25*B(516)

  JVS(573) = 0.3*B(546)

  JVS(574) = 0.25*B(289)

  JVS(575) = 0.25*B(269)

  JVS(576) = -B(555)

  JVS(577) = B(117)+0.3*B(148)+0.5*B(196)+0.25*B(270)+0.25*B(290)+0.25*B(495)+0.25*B(509)+0.25*B(517)+0.3*B(547)+0.25&
               &*B(581)

  JVS(578) = 0.9*B(36)

  JVS(579) = 0.1*B(460)

  JVS(580) = 0.1*B(450)

  JVS(581) = 0

  JVS(582) = -B(32)

  JVS(583) = 0.1*B(470)

  JVS(584) = 0.1*B(480)

  JVS(585) = 0.1*B(482)

  JVS(586) = 0

  JVS(587) = 0

  JVS(588) = 0.1*B(451)+0.1*B(461)+0.1*B(471)+0.1*B(481)+0.1*B(483)

  JVS(589) = 0

  JVS(590) = B(141)

  JVS(591) = 0.5*B(149)

  JVS(592) = B(31)

  JVS(593) = 0

  JVS(594) = -B(143)-B(145)-B(147)-2*B(575)

  JVS(595) = 0

  JVS(596) = B(142)+0.5*B(150)

  JVS(597) = -B(146)

  JVS(598) = -B(144)

  JVS(599) = -B(148)

  JVS(600) = -B(524)-B(526)-B(528)

  JVS(601) = B(522)

  JVS(602) = B(523)-B(527)

  JVS(603) = -B(529)

  JVS(604) = -B(525)

  JVS(605) = -B(436)-B(470)-B(484)

  JVS(606) = -B(485)

  JVS(607) = -B(437)

  JVS(608) = -B(471)

  JVS(609) = B(594)

  JVS(610) = -B(596)-B(598)-B(600)

  JVS(611) = B(595)-B(601)

  JVS(612) = 0

  JVS(613) = -B(597)

  JVS(614) = 0

  JVS(615) = -B(599)

  JVS(616) = 0.33*B(43)

  JVS(617) = 0.9*B(592)

  JVS(618) = 0.25*B(598)+0.25*B(600)

  JVS(619) = -B(44)-B(582)

  JVS(620) = 0.33*B(578)+0.167*B(580)

  JVS(621) = 0.25*B(601)

  JVS(622) = -B(583)

  JVS(623) = 0

  JVS(624) = 0.9*B(593)

  JVS(625) = 0.33*B(579)+0.25*B(599)

  JVS(626) = 0.167*B(581)

  JVS(627) = -B(446)-B(480)-B(486)

  JVS(628) = -B(487)

  JVS(629) = -B(447)

  JVS(630) = -B(481)

  JVS(631) = -B(241)-B(251)-B(522)

  JVS(632) = -B(523)

  JVS(633) = -B(242)

  JVS(634) = -B(252)

  JVS(635) = -B(448)-B(482)-B(488)

  JVS(636) = -B(489)

  JVS(637) = -B(449)

  JVS(638) = -B(483)

  JVS(639) = 0.82*B(22)

  JVS(640) = 0.5*B(209)

  JVS(641) = 0.05*B(47)

  JVS(642) = B(584)+0.4*B(586)

  JVS(643) = 0.67*B(43)

  JVS(644) = 0.3*B(48)

  JVS(645) = 0.52*B(460)

  JVS(646) = 0.52*B(450)

  JVS(647) = 0.25*B(34)

  JVS(648) = 0.1*B(592)

  JVS(649) = 0.82*B(181)+0.82*B(185)

  JVS(650) = 0.52*B(470)

  JVS(651) = 0.75*B(598)+0.75*B(600)

  JVS(652) = B(44)

  JVS(653) = 0.52*B(480)

  JVS(654) = 0.52*B(482)

  JVS(655) = -B(24)-B(189)

  JVS(656) = 0.5*B(46)+0.5*B(510)

  JVS(657) = 0.67*B(578)+0.333*B(580)

  JVS(658) = 0.25*B(221)

  JVS(659) = 0.27*B(504)+0.15*B(508)

  JVS(660) = 0.045*B(490)+0.025*B(494)

  JVS(661) = 0.75*B(601)

  JVS(662) = -B(190)+0.5*B(511)

  JVS(663) = 0.4*B(587)

  JVS(664) = 0.52*B(451)+0.52*B(461)+0.52*B(471)+0.52*B(481)+0.52*B(483)+0.1*B(593)

  JVS(665) = 0.82*B(182)+0.5*B(210)+0.25*B(222)+0.045*B(491)+0.27*B(505)+0.67*B(579)+B(585)+0.75*B(599)

  JVS(666) = 0.82*B(186)+0.025*B(495)+0.15*B(509)+0.333*B(581)

  JVS(667) = 0.35*B(602)

  JVS(668) = 1.5*B(38)

  JVS(669) = 1.5*B(39)

  JVS(670) = 0

  JVS(671) = B(42)

  JVS(672) = 0.6*B(41)

  JVS(673) = 0.5*B(130)

  JVS(674) = 0.25*B(48)

  JVS(675) = 0.23*B(460)

  JVS(676) = 0.23*B(450)

  JVS(677) = 0.4*B(347)+0.16*B(349)

  JVS(678) = 0.35*B(592)

  JVS(679) = 0

  JVS(680) = 0

  JVS(681) = 0.07*B(380)+0.17*B(382)

  JVS(682) = 0.14*B(384)+0.35*B(386)

  JVS(683) = 2*B(33)+B(134)

  JVS(684) = 0.45*B(32)

  JVS(685) = 0.23*B(470)

  JVS(686) = B(44)

  JVS(687) = 0.23*B(480)

  JVS(688) = 0.3*B(251)

  JVS(689) = 0.23*B(482)

  JVS(690) = -B(126)

  JVS(691) = B(30)

  JVS(692) = 0.56*B(153)

  JVS(693) = B(14)

  JVS(694) = 1.7*B(46)+0.7*B(510)

  JVS(695) = 0

  JVS(696) = 0.225*B(504)+0.125*B(508)

  JVS(697) = B(27)

  JVS(698) = 0.67*B(19)+0.2*B(279)

  JVS(699) = B(45)

  JVS(700) = B(25)+B(530)+B(532)

  JVS(701) = 0.25*B(540)+0.25*B(542)+0.2*B(546)+0.25*B(548)

  JVS(702) = 0.7*B(20)+0.05*B(275)

  JVS(703) = 0.22*B(281)+0.22*B(285)+0.11*B(289)+0.22*B(291)

  JVS(704) = B(11)+B(12)+B(122)+B(124)

  JVS(705) = 0

  JVS(706) = B(123)+0.22*B(286)+B(533)+0.25*B(543)

  JVS(707) = B(125)-B(127)+B(135)+0.7*B(511)+B(531)+0.35*B(603)

  JVS(708) = 0.16*B(350)+0.07*B(381)+0.14*B(385)

  JVS(709) = 0.5*B(131)+0.56*B(154)+0.3*B(252)+0.05*B(276)+0.2*B(280)+0.23*B(451)+0.23*B(461)+0.23*B(471)+0.23*B(481)&
               &+0.23*B(483)+0.35*B(593)

  JVS(710) = 0.22*B(292)+0.25*B(549)

  JVS(711) = 0

  JVS(712) = 0.22*B(282)+0.4*B(348)+0.17*B(383)+0.35*B(387)+0.225*B(505)+0.25*B(541)

  JVS(713) = 0.11*B(290)+0.125*B(509)+0.2*B(547)

  JVS(714) = B(138)

  JVS(715) = 0

  JVS(716) = 0.67*B(43)

  JVS(717) = 0.25*B(48)

  JVS(718) = -B(30)-B(564)

  JVS(719) = 0.67*B(578)+0.333*B(580)

  JVS(720) = 0.225*B(504)+0.125*B(508)

  JVS(721) = 0.25*B(540)+0.25*B(542)+0.1*B(546)+0.25*B(548)

  JVS(722) = 0.53*B(281)+0.53*B(285)+0.26*B(289)+0.53*B(291)

  JVS(723) = 0.02*B(261)+0.02*B(263)

  JVS(724) = 0.02*B(264)+0.53*B(286)+0.25*B(543)

  JVS(725) = -B(565)

  JVS(726) = 0

  JVS(727) = 0

  JVS(728) = 0.53*B(292)+0.25*B(549)

  JVS(729) = 0.02*B(262)+0.53*B(282)+0.225*B(505)+0.25*B(541)+0.67*B(579)

  JVS(730) = 0.26*B(290)+0.125*B(509)+0.1*B(547)+0.333*B(581)

  JVS(731) = 0.07*B(251)

  JVS(732) = -B(151)-B(153)-B(155)

  JVS(733) = 0.7*B(20)

  JVS(734) = -B(156)

  JVS(735) = -B(152)

  JVS(736) = -B(154)+0.07*B(252)

  JVS(737) = B(556)

  JVS(738) = B(21)+0.5*B(149)

  JVS(739) = B(209)

  JVS(740) = B(35)

  JVS(741) = B(15)

  JVS(742) = 0.4*B(34)

  JVS(743) = B(231)

  JVS(744) = B(157)

  JVS(745) = 0.27*B(181)

  JVS(746) = B(143)+0.8*B(147)+1.6*B(575)

  JVS(747) = 0.5*B(153)

  JVS(748) = -B(14)-B(163)-B(165)

  JVS(749) = 0.4*B(221)

  JVS(750) = 0.04*B(275)

  JVS(751) = -B(166)

  JVS(752) = 0.5*B(150)-B(164)+B(557)

  JVS(753) = 0

  JVS(754) = 0.5*B(154)+0.04*B(276)

  JVS(755) = B(144)+B(158)+0.27*B(182)+B(210)+0.4*B(222)+B(232)

  JVS(756) = 0.8*B(148)

  JVS(757) = B(48)

  JVS(758) = 0.3*B(460)

  JVS(759) = 0.3*B(450)

  JVS(760) = 0.3*B(470)

  JVS(761) = 0.3*B(480)

  JVS(762) = 0.3*B(482)

  JVS(763) = -B(46)-B(510)

  JVS(764) = 0.9*B(504)+B(508)

  JVS(765) = B(45)

  JVS(766) = 0

  JVS(767) = -B(511)

  JVS(768) = 0

  JVS(769) = 0.3*B(451)+0.3*B(461)+0.3*B(471)+0.3*B(481)+0.3*B(483)

  JVS(770) = 0.9*B(505)

  JVS(771) = B(509)

  JVS(772) = 0.5*B(590)

  JVS(773) = B(576)

  JVS(774) = -B(578)-B(580)-B(588)

  JVS(775) = 0

  JVS(776) = B(577)+0.5*B(591)

  JVS(777) = -B(589)

  JVS(778) = 0

  JVS(779) = -B(579)

  JVS(780) = -B(581)

  JVS(781) = B(211)

  JVS(782) = 0

  JVS(783) = 0

  JVS(784) = B(225)

  JVS(785) = 0

  JVS(786) = 0

  JVS(787) = 0

  JVS(788) = -B(221)-B(223)

  JVS(789) = 0

  JVS(790) = B(212)+B(226)

  JVS(791) = -B(224)

  JVS(792) = 0

  JVS(793) = -B(222)

  JVS(794) = B(498)

  JVS(795) = -B(504)-B(506)-B(508)

  JVS(796) = B(500)+0.5*B(502)

  JVS(797) = 0.5*B(503)

  JVS(798) = B(499)+B(501)

  JVS(799) = -B(507)

  JVS(800) = -B(505)

  JVS(801) = -B(509)

  JVS(802) = B(496)

  JVS(803) = B(426)

  JVS(804) = B(416)

  JVS(805) = B(436)

  JVS(806) = B(446)

  JVS(807) = B(448)

  JVS(808) = -B(490)-B(492)-B(494)

  JVS(809) = 0

  JVS(810) = B(417)+B(427)+B(437)+B(447)+B(449)+B(497)

  JVS(811) = -B(493)

  JVS(812) = 0

  JVS(813) = -B(491)

  JVS(814) = -B(495)

  JVS(815) = 0.794*B(524)+0.794*B(526)+0.794*B(528)

  JVS(816) = 0

  JVS(817) = 0.1*B(504)

  JVS(818) = 0.1*B(490)

  JVS(819) = -B(27)-B(534)-B(536)

  JVS(820) = 0.1*B(512)+B(514)+0.5*B(516)

  JVS(821) = 0

  JVS(822) = 0.8*B(283)

  JVS(823) = 0.08*B(261)

  JVS(824) = 0.794*B(527)-B(537)

  JVS(825) = -B(535)

  JVS(826) = B(515)+0.794*B(529)

  JVS(827) = 0

  JVS(828) = 0.08*B(262)+0.8*B(284)+0.1*B(491)+0.1*B(505)+0.1*B(513)+0.794*B(525)

  JVS(829) = 0.5*B(517)

  JVS(830) = B(197)

  JVS(831) = 0.06*B(460)

  JVS(832) = 0.06*B(450)

  JVS(833) = 0.06*B(470)

  JVS(834) = 0.06*B(480)

  JVS(835) = 0.06*B(482)

  JVS(836) = B(189)

  JVS(837) = 0.15*B(46)+0.15*B(510)

  JVS(838) = 0

  JVS(839) = 0

  JVS(840) = 0

  JVS(841) = 0

  JVS(842) = -B(191)-B(193)-B(195)

  JVS(843) = 0

  JVS(844) = 0

  JVS(845) = B(190)+B(198)+0.15*B(511)

  JVS(846) = -B(194)

  JVS(847) = 0.06*B(451)+0.06*B(461)+0.06*B(471)+0.06*B(481)+0.06*B(483)

  JVS(848) = -B(192)

  JVS(849) = -B(196)

  JVS(850) = B(484)

  JVS(851) = B(486)

  JVS(852) = B(488)

  JVS(853) = -B(512)-B(514)-B(516)-B(518)

  JVS(854) = 0.5*B(502)

  JVS(855) = B(485)+B(487)+B(489)+0.5*B(503)-B(519)

  JVS(856) = 0

  JVS(857) = -B(515)

  JVS(858) = 0

  JVS(859) = -B(513)

  JVS(860) = -B(517)

  JVS(861) = 0.288*B(28)

  JVS(862) = 0.167*B(524)+0.167*B(526)+0.167*B(528)

  JVS(863) = 0.4*B(251)

  JVS(864) = -B(19)-B(277)-B(279)

  JVS(865) = 0.23*B(261)+0.25*B(263)+0.19*B(269)+0.25*B(271)

  JVS(866) = 0.25*B(264)+0.167*B(527)

  JVS(867) = -B(278)

  JVS(868) = 0.167*B(529)

  JVS(869) = 0.4*B(252)-B(280)

  JVS(870) = 0.25*B(272)

  JVS(871) = 0.23*B(262)+0.167*B(525)

  JVS(872) = 0.19*B(270)

  JVS(873) = B(47)

  JVS(874) = 0.33*B(460)

  JVS(875) = 0.33*B(450)

  JVS(876) = 0.33*B(470)

  JVS(877) = 0.33*B(480)

  JVS(878) = 0.33*B(482)

  JVS(879) = 0.9*B(490)+B(494)

  JVS(880) = 0.9*B(512)+0.5*B(516)+B(518)

  JVS(881) = -B(45)-B(500)-B(502)

  JVS(882) = -B(503)+B(519)

  JVS(883) = -B(501)

  JVS(884) = 0

  JVS(885) = 0.33*B(451)+0.33*B(461)+0.33*B(471)+0.33*B(481)+0.33*B(483)

  JVS(886) = 0.9*B(491)+0.9*B(513)

  JVS(887) = B(495)+0.5*B(517)

  JVS(888) = 0.5*B(161)

  JVS(889) = 0.5*B(558)

  JVS(890) = 0

  JVS(891) = 0

  JVS(892) = 0.2*B(195)

  JVS(893) = 0

  JVS(894) = -B(29)-B(562)

  JVS(895) = 0.25*B(540)+0.25*B(542)+0.1*B(546)+0.25*B(548)

  JVS(896) = 0

  JVS(897) = 0.22*B(281)+0.22*B(285)+0.23*B(289)+0.22*B(291)

  JVS(898) = 0.02*B(261)+0.02*B(263)

  JVS(899) = 0

  JVS(900) = 0.02*B(264)+0.22*B(286)+0.25*B(543)

  JVS(901) = 0.5*B(162)+0.5*B(559)-B(563)

  JVS(902) = 0

  JVS(903) = 0

  JVS(904) = 0.22*B(292)+0.25*B(549)

  JVS(905) = 0

  JVS(906) = 0.02*B(262)+0.22*B(282)+0.25*B(541)

  JVS(907) = 0.2*B(196)+0.23*B(290)+0.1*B(547)

  JVS(908) = B(42)

  JVS(909) = 0.45*B(36)

  JVS(910) = B(227)

  JVS(911) = 0.51*B(404)

  JVS(912) = 0.4*B(378)

  JVS(913) = 0.54*B(414)

  JVS(914) = 0.07*B(380)+0.17*B(382)

  JVS(915) = 0.07*B(384)+0.17*B(386)

  JVS(916) = 0.18*B(32)

  JVS(917) = 0

  JVS(918) = 0

  JVS(919) = 0

  JVS(920) = 0

  JVS(921) = 0

  JVS(922) = 0

  JVS(923) = 0.5*B(195)

  JVS(924) = 0.8*B(279)

  JVS(925) = 0

  JVS(926) = B(562)

  JVS(927) = -B(25)-B(530)-B(532)

  JVS(928) = 0.25*B(540)+0.25*B(542)+0.1*B(546)+0.25*B(548)

  JVS(929) = 0.95*B(275)

  JVS(930) = 0.25*B(281)+0.25*B(285)+0.24*B(289)+0.25*B(291)

  JVS(931) = 0.02*B(261)+0.02*B(263)

  JVS(932) = 0

  JVS(933) = 0.02*B(264)+0.25*B(286)-B(533)+0.25*B(543)

  JVS(934) = B(228)-B(531)+B(563)

  JVS(935) = 0.07*B(381)+0.07*B(385)

  JVS(936) = 0.95*B(276)+0.8*B(280)

  JVS(937) = 0.25*B(292)+0.25*B(549)

  JVS(938) = 0

  JVS(939) = 0.02*B(262)+0.25*B(282)+0.4*B(379)+0.17*B(383)+0.17*B(387)+0.51*B(405)+0.54*B(415)+0.25*B(541)

  JVS(940) = 0.5*B(196)+0.24*B(290)+0.1*B(547)

  JVS(941) = B(550)

  JVS(942) = B(538)

  JVS(943) = 0.5*B(267)

  JVS(944) = 0

  JVS(945) = 0

  JVS(946) = 0

  JVS(947) = -B(540)-B(542)-B(544)-B(546)-B(548)

  JVS(948) = 0

  JVS(949) = 0

  JVS(950) = -B(543)

  JVS(951) = 0.5*B(268)+B(539)+B(551)

  JVS(952) = -B(545)

  JVS(953) = 0

  JVS(954) = -B(549)

  JVS(955) = -B(541)

  JVS(956) = -B(547)

  JVS(957) = 0.402*B(28)

  JVS(958) = 0.039*B(524)+0.039*B(526)+0.039*B(528)

  JVS(959) = 0.2*B(251)

  JVS(960) = -B(20)-B(273)-B(275)

  JVS(961) = 0.32*B(261)+0.35*B(263)+0.26*B(269)+0.35*B(271)

  JVS(962) = 0.35*B(264)+0.039*B(527)

  JVS(963) = -B(274)

  JVS(964) = 0.039*B(529)

  JVS(965) = 0.2*B(252)-B(276)

  JVS(966) = 0.35*B(272)

  JVS(967) = 0.32*B(262)+0.039*B(525)

  JVS(968) = 0.26*B(270)

  JVS(969) = 0.2*B(293)

  JVS(970) = 0.5*B(277)

  JVS(971) = B(273)

  JVS(972) = -B(281)-B(283)-B(285)-B(287)-B(289)-B(291)

  JVS(973) = 0

  JVS(974) = -B(286)

  JVS(975) = B(274)+0.5*B(278)+0.2*B(294)

  JVS(976) = -B(288)

  JVS(977) = 0

  JVS(978) = -B(292)

  JVS(979) = -B(282)-B(284)

  JVS(980) = -B(290)

  JVS(981) = 2*B(140)

  JVS(982) = B(23)

  JVS(983) = 0.5*B(209)

  JVS(984) = B(10)+0.3*B(120)

  JVS(985) = 0.5*B(128)+B(130)

  JVS(986) = 0

  JVS(987) = B(608)

  JVS(988) = 0.5*B(175)

  JVS(989) = 0.4*B(47)

  JVS(990) = B(15)

  JVS(991) = B(560)

  JVS(992) = 0.25*B(112)

  JVS(993) = 0.33*B(43)

  JVS(994) = 0.69*B(28)

  JVS(995) = 0.5*B(558)

  JVS(996) = 0.375*B(48)

  JVS(997) = 0.25*B(113)

  JVS(998) = 0.34*B(460)

  JVS(999) = 0.34*B(450)

  JVS(1000) = 0.1*B(34)

  JVS(1001) = B(157)

  JVS(1002) = 0.5*B(592)

  JVS(1003) = B(185)

  JVS(1004) = B(554)

  JVS(1005) = 0.7*B(147)

  JVS(1006) = 0.072*B(524)+0.072*B(526)+0.206*B(528)

  JVS(1007) = 0.34*B(470)

  JVS(1008) = 0.25*B(598)+0.25*B(600)

  JVS(1009) = 0.34*B(480)

  JVS(1010) = 0.6*B(251)

  JVS(1011) = 0.34*B(482)

  JVS(1012) = B(30)+0.8*B(564)

  JVS(1013) = 0.54*B(153)

  JVS(1014) = 0.68*B(46)+0.68*B(510)

  JVS(1015) = 0.33*B(578)+0.917*B(580)

  JVS(1016) = 0.1*B(221)

  JVS(1017) = 0.34*B(504)+0.93*B(508)

  JVS(1018) = 0.36*B(490)+0.95*B(494)

  JVS(1019) = B(27)

  JVS(1020) = B(191)+0.8*B(195)

  JVS(1021) = 0.75*B(516)

  JVS(1022) = 0.67*B(19)+0.7*B(279)

  JVS(1023) = 0

  JVS(1024) = B(29)

  JVS(1025) = 0.25*B(540)+0.25*B(542)+0.8*B(546)+0.25*B(548)

  JVS(1026) = 0.8*B(275)

  JVS(1027) = 0.25*B(281)+0.25*B(285)+0.88*B(289)+0.25*B(291)

  JVS(1028) = -B(11)-B(12)-B(122)-B(124)-B(606)

  JVS(1029) = 0.55*B(261)+0.6*B(263)+1.2*B(269)+0.6*B(271)

  JVS(1030) = B(295)+B(297)+2*B(301)+B(303)+2*B(305)

  JVS(1031) = -B(123)+0.6*B(264)+0.25*B(286)+B(298)+0.072*B(527)+0.25*B(543)+0.25*B(601)

  JVS(1032) = 0.3*B(121)-B(125)+0.5*B(129)+0.5*B(176)+0.68*B(511)+B(555)+0.5*B(559)+B(561)+0.8*B(565)

  JVS(1033) = 0.206*B(529)-B(607)

  JVS(1034) = B(131)+0.54*B(154)+0.6*B(252)+0.8*B(276)+0.7*B(280)+0.34*B(451)+0.34*B(461)+0.34*B(471)+0.34*B(481)+0.34&
                &*B(483)+0.5*B(593)

  JVS(1035) = B(173)+0.6*B(272)+0.25*B(292)+B(304)+0.25*B(549)

  JVS(1036) = 0

  JVS(1037) = B(114)+B(158)+B(192)+0.5*B(210)+0.1*B(222)+0.55*B(262)+0.25*B(282)+B(296)+0.36*B(491)+0.34*B(505)+0.072&
                &*B(525)+0.25*B(541)+0.33*B(579)+0.25*B(599)

  JVS(1038) = B(115)+2*B(116)+B(117)+0.7*B(148)+B(174)+B(186)+0.8*B(196)+1.2*B(270)+0.88*B(290)+2*B(302)+0.95*B(495)&
                &+0.93*B(509)+0.75*B(517)+0.8*B(547)+0.917*B(581)

  JVS(1039) = 0.5*B(267)

  JVS(1040) = B(241)

  JVS(1041) = -B(261)-B(263)-B(265)-B(269)-B(271)

  JVS(1042) = -B(264)

  JVS(1043) = B(242)+0.5*B(268)

  JVS(1044) = -B(266)

  JVS(1045) = 0

  JVS(1046) = -B(272)

  JVS(1047) = -B(262)

  JVS(1048) = -B(270)

  JVS(1049) = 0.5*B(293)

  JVS(1050) = B(18)+B(309)

  JVS(1051) = 0.2*B(251)

  JVS(1052) = 0.33*B(19)+0.5*B(277)

  JVS(1053) = 0

  JVS(1054) = 0

  JVS(1055) = -B(295)-B(297)-B(299)-B(301)-B(303)-2*B(305)-B(306)

  JVS(1056) = -B(298)

  JVS(1057) = 0.5*B(278)+0.5*B(294)

  JVS(1058) = -B(300)

  JVS(1059) = 0.2*B(252)

  JVS(1060) = -B(304)

  JVS(1061) = -B(307)

  JVS(1062) = -B(296)

  JVS(1063) = -B(302)

  JVS(1064) = B(6)+B(94)

  JVS(1065) = -B(570)

  JVS(1066) = 0.33*B(9)

  JVS(1067) = 0.4*B(17)+B(560)

  JVS(1068) = B(97)

  JVS(1069) = 0.5*B(558)

  JVS(1070) = -B(594)

  JVS(1071) = -B(526)

  JVS(1072) = -B(484)

  JVS(1073) = -B(600)

  JVS(1074) = -B(486)

  JVS(1075) = -B(522)

  JVS(1076) = -B(488)

  JVS(1077) = -B(155)

  JVS(1078) = -B(165)

  JVS(1079) = 0

  JVS(1080) = -B(536)

  JVS(1081) = -B(518)

  JVS(1082) = -B(502)

  JVS(1083) = -B(532)

  JVS(1084) = -B(542)

  JVS(1085) = 0

  JVS(1086) = -B(285)

  JVS(1087) = -B(122)

  JVS(1088) = -B(263)

  JVS(1089) = -B(297)

  JVS(1090) = -B(8)-B(90)-B(92)-B(99)-B(108)-B(123)-B(156)-B(166)-B(264)-B(286)-B(298)-B(485)-B(487)-B(489)-B(503)&
                &-B(519)-B(523)-B(527)-B(533)-B(537)-B(543)-B(571)-B(595)-B(601)

  JVS(1091) = B(98)+0.5*B(559)+B(561)

  JVS(1092) = -B(91)

  JVS(1093) = B(88)

  JVS(1094) = 0

  JVS(1095) = B(89)-B(93)

  JVS(1096) = -B(100)

  JVS(1097) = 0

  JVS(1098) = B(49)

  JVS(1099) = -B(572)

  JVS(1100) = -B(141)

  JVS(1101) = -B(179)

  JVS(1102) = -B(199)

  JVS(1103) = -B(132)

  JVS(1104) = 2*B(13)-B(72)

  JVS(1105) = -0.35*B(602)

  JVS(1106) = -B(311)

  JVS(1107) = -B(321)

  JVS(1108) = -B(351)

  JVS(1109) = -B(361)

  JVS(1110) = B(26)-B(550)

  JVS(1111) = -B(327)

  JVS(1112) = -B(337)

  JVS(1113) = -B(343)

  JVS(1114) = -B(365)

  JVS(1115) = -B(369)

  JVS(1116) = -B(408)

  JVS(1117) = -B(412)

  JVS(1118) = -B(392)

  JVS(1119) = -B(402)

  JVS(1120) = -B(566)-B(568)

  JVS(1121) = -B(556)

  JVS(1122) = B(21)-0.5*B(149)

  JVS(1123) = B(22)-B(187)

  JVS(1124) = B(23)-B(197)

  JVS(1125) = B(35)-B(235)

  JVS(1126) = B(10)-0.7*B(120)

  JVS(1127) = -0.9*B(293)

  JVS(1128) = B(58)-B(60)

  JVS(1129) = 0.33*B(9)-B(103)

  JVS(1130) = -B(538)

  JVS(1131) = -B(211)

  JVS(1132) = -B(128)+0.12*B(130)

  JVS(1133) = B(36)-B(239)

  JVS(1134) = 0

  JVS(1135) = 0

  JVS(1136) = -B(229)

  JVS(1137) = -B(520)

  JVS(1138) = B(16)-B(175)

  JVS(1139) = B(47)-B(496)

  JVS(1140) = 0

  JVS(1141) = -B(227)

  JVS(1142) = B(15)-0.5*B(161)

  JVS(1143) = -B(560)

  JVS(1144) = 0

  JVS(1145) = B(7)-B(97)

  JVS(1146) = 0.4*B(586)

  JVS(1147) = -B(110)+0.75*B(112)

  JVS(1148) = 0.4*B(376)

  JVS(1149) = B(43)-0.5*B(590)

  JVS(1150) = -B(267)

  JVS(1151) = -B(558)

  JVS(1152) = -B(604)

  JVS(1153) = B(48)-B(498)

  JVS(1154) = 2*B(56)+B(59)+0.75*B(113)

  JVS(1155) = 0

  JVS(1156) = -B(426)+0.63*B(460)

  JVS(1157) = -B(416)+0.63*B(450)

  JVS(1158) = B(34)-B(225)

  JVS(1159) = 0

  JVS(1160) = 0

  JVS(1161) = -B(576)+0.06*B(592)

  JVS(1162) = 0

  JVS(1163) = 0

  JVS(1164) = 0.4*B(380)

  JVS(1165) = -B(62)+B(64)

  JVS(1166) = 0

  JVS(1167) = 0.4*B(384)

  JVS(1168) = -B(134)

  JVS(1169) = -B(554)

  JVS(1170) = 0

  JVS(1171) = 0

  JVS(1172) = 0.206*B(528)

  JVS(1173) = -B(436)+0.63*B(470)

  JVS(1174) = 0

  JVS(1175) = -B(582)

  JVS(1176) = -B(446)+0.63*B(480)

  JVS(1177) = -B(241)+0.27*B(251)

  JVS(1178) = -B(448)+0.63*B(482)

  JVS(1179) = -B(189)

  JVS(1180) = -B(126)

  JVS(1181) = -B(564)

  JVS(1182) = -B(151)+0.33*B(153)

  JVS(1183) = -B(163)

  JVS(1184) = -B(510)

  JVS(1185) = 0

  JVS(1186) = 0

  JVS(1187) = 0

  JVS(1188) = 0

  JVS(1189) = -B(534)

  JVS(1190) = 0

  JVS(1191) = 0

  JVS(1192) = 0.33*B(19)-B(277)+0.215*B(279)

  JVS(1193) = -B(500)

  JVS(1194) = -B(562)

  JVS(1195) = -B(530)

  JVS(1196) = 0

  JVS(1197) = -B(273)+0.08*B(275)

  JVS(1198) = 0

  JVS(1199) = -B(124)

  JVS(1200) = 0

  JVS(1201) = 0

  JVS(1202) = B(90)

  JVS(1203) = -B(61)-B(63)-B(66)-B(73)-B(74)-2*B(76)-2*B(77)-B(95)-B(98)-B(104)-B(111)-0.7*B(121)-B(125)-B(127)-B(129)&
                &-B(133)-B(135)-B(142)-0.5*B(150)-B(152)-0.5*B(162)-B(164)-B(176)-B(180)-B(188)-B(190)-B(198)-B(200)-B(212)&
                &-B(226)-B(228)-B(230)-B(236)-B(240)-B(242)-B(268)-B(274)-B(278)-0.9*B(294)-B(312)-B(322)-B(328)-B(338)&
                &-B(344)-B(352)-B(362)-B(366)-B(370)-B(393)-B(403)-B(409)-B(413)-B(417)-B(427)-B(437)-B(447)-B(449)-B(497)&
                &-B(499)-B(501)-B(511)-B(521)-B(531)-B(535)-B(539)-B(551)-B(555)-B(557)-B(559)-B(561)-B(563)-B(565)-B(567)&
                &-B(569)-B(573)-B(577)-B(583)-0.5*B(591)-0.35*B(603)-B(605)

  JVS(1204) = B(65)+B(68)-B(75)+B(82)+B(91)+0.4*B(377)+0.4*B(381)+0.4*B(385)+0.206*B(529)+0.4*B(587)

  JVS(1205) = -B(67)+B(69)+0.12*B(131)+0.33*B(154)+0.27*B(252)+0.08*B(276)+0.215*B(280)+0.63*B(451)+0.63*B(461)+0.63&
                &*B(471)+0.63*B(481)+0.63*B(483)+0.06*B(593)

  JVS(1206) = 0

  JVS(1207) = -B(96)+0.5*B(109)

  JVS(1208) = B(83)

  JVS(1209) = 0

  JVS(1210) = B(72)

  JVS(1211) = 0.35*B(602)

  JVS(1212) = 0.65*B(311)

  JVS(1213) = 1.5*B(38)

  JVS(1214) = 0.8*B(321)

  JVS(1215) = 0.28*B(351)

  JVS(1216) = 0.73*B(361)

  JVS(1217) = B(138)+B(140)

  JVS(1218) = 0.6*B(40)

  JVS(1219) = B(39)

  JVS(1220) = 0.38*B(392)

  JVS(1221) = 0.63*B(402)

  JVS(1222) = 0.5*B(568)

  JVS(1223) = B(42)

  JVS(1224) = B(556)

  JVS(1225) = 0.6*B(41)

  JVS(1226) = B(37)

  JVS(1227) = B(21)

  JVS(1228) = B(22)

  JVS(1229) = B(209)

  JVS(1230) = B(10)

  JVS(1231) = 0.2*B(293)

  JVS(1232) = B(58)+B(60)

  JVS(1233) = 0.66*B(9)+B(105)

  JVS(1234) = 0.25*B(128)+0.12*B(130)

  JVS(1235) = 0

  JVS(1236) = B(323)-B(325)

  JVS(1237) = -B(363)+B(367)

  JVS(1238) = B(608)+B(609)-B(611)

  JVS(1239) = B(47)

  JVS(1240) = B(339)-B(341)

  JVS(1241) = B(15)

  JVS(1242) = B(404)-B(406)

  JVS(1243) = B(584)-0.6*B(586)

  JVS(1244) = 0.4*B(112)

  JVS(1245) = -B(376)

  JVS(1246) = B(43)

  JVS(1247) = B(28)

  JVS(1248) = 0.5*B(558)

  JVS(1249) = B(604)

  JVS(1250) = B(48)

  JVS(1251) = B(59)+0.4*B(113)

  JVS(1252) = -B(335)

  JVS(1253) = 0.57*B(460)

  JVS(1254) = 0.57*B(450)

  JVS(1255) = 0.9*B(34)

  JVS(1256) = -B(233)

  JVS(1257) = B(157)-B(159)

  JVS(1258) = 0.4*B(347)-0.84*B(349)

  JVS(1259) = 0.06*B(592)

  JVS(1260) = -B(237)+B(378)

  JVS(1261) = -B(410)+B(414)

  JVS(1262) = -0.93*B(380)+0.17*B(382)

  JVS(1263) = B(62)-B(64)

  JVS(1264) = B(181)-B(183)+B(185)

  JVS(1265) = -0.67*B(384)+0.83*B(386)

  JVS(1266) = 2*B(33)+B(134)

  JVS(1267) = B(554)

  JVS(1268) = 0.56*B(32)

  JVS(1269) = B(143)-B(145)+B(147)+1.2*B(575)

  JVS(1270) = 0.794*B(524)+0.794*B(526)-B(528)

  JVS(1271) = 0.57*B(470)

  JVS(1272) = -B(596)+0.75*B(598)+0.75*B(600)

  JVS(1273) = 2*B(44)

  JVS(1274) = 0.57*B(480)

  JVS(1275) = 0.06*B(251)

  JVS(1276) = 0.57*B(482)

  JVS(1277) = B(126)

  JVS(1278) = 2*B(30)+B(564)

  JVS(1279) = 0.19*B(153)

  JVS(1280) = B(14)

  JVS(1281) = 1.2*B(46)+0.2*B(510)

  JVS(1282) = B(578)+B(580)-B(588)

  JVS(1283) = 0.9*B(221)-B(223)

  JVS(1284) = 0.9*B(504)-B(506)+B(508)

  JVS(1285) = 0.9*B(490)-B(492)+B(494)

  JVS(1286) = B(27)+B(534)+B(536)

  JVS(1287) = -B(193)+0.3*B(195)

  JVS(1288) = -B(514)+0.5*B(516)

  JVS(1289) = 0.67*B(19)+0.275*B(279)

  JVS(1290) = B(45)

  JVS(1291) = B(29)+B(562)

  JVS(1292) = B(25)

  JVS(1293) = B(540)+B(542)-B(544)+0.8*B(546)+B(548)

  JVS(1294) = 0.06*B(275)

  JVS(1295) = 0.47*B(281)+0.47*B(285)-B(287)+0.73*B(289)+0.47*B(291)

  JVS(1296) = 2*B(11)+B(122)+B(124)-B(606)

  JVS(1297) = 0.92*B(261)+B(263)-B(265)+B(269)+B(271)

  JVS(1298) = -B(299)+B(301)

  JVS(1299) = -B(90)+B(123)+B(264)+0.47*B(286)+0.794*B(527)+B(537)+B(543)+0.75*B(601)

  JVS(1300) = B(61)+B(63)+B(66)+B(73)-B(74)+B(125)+B(127)+0.25*B(129)+B(135)+0.2*B(294)+0.65*B(312)+0.8*B(322)+0.28&
                &*B(352)+0.73*B(362)+0.38*B(393)+0.63*B(403)+0.2*B(511)+B(535)+B(555)+B(557)+0.5*B(559)+B(563)+B(565)+0.5&
                &*B(569)+0.35*B(603)+B(605)

  JVS(1301) = -B(65)-B(68)-2*B(70)-B(75)-B(82)-B(91)-B(101)-B(118)-B(146)-B(160)-B(171)-B(184)-B(194)-B(224)-B(234)&
                &-B(238)-B(266)-B(288)-B(300)-B(326)-B(336)-B(342)-0.84*B(350)-B(364)-B(377)-0.93*B(381)-0.67*B(385)-B(407)&
                &-B(411)-B(493)-B(507)-B(515)-B(529)-B(545)-B(574)-0.6*B(587)-B(589)-B(597)-B(607)-B(612)

  JVS(1302) = B(67)-B(69)+0.12*B(131)+0.19*B(154)+0.06*B(252)+0.06*B(276)+0.275*B(280)+0.57*B(451)+0.57*B(461)+0.57&
                &*B(471)+0.57*B(481)+0.57*B(483)+0.06*B(593)

  JVS(1303) = -B(172)+0.9*B(173)+B(272)+0.47*B(292)+B(549)

  JVS(1304) = -B(102)

  JVS(1305) = -B(83)+B(114)+B(144)+B(158)+B(182)+B(210)+0.9*B(222)+0.92*B(262)+0.47*B(282)+B(324)+B(340)+0.4*B(348)&
                &+B(368)+B(379)+0.17*B(383)+0.83*B(387)+B(405)+B(415)+0.9*B(491)+0.9*B(505)+0.794*B(525)+B(541)+B(579)&
                &+B(585)+0.75*B(599)+B(610)

  JVS(1306) = B(115)+2*B(116)-B(119)+B(148)+0.9*B(174)+B(186)+0.3*B(196)+B(270)+0.73*B(290)+B(302)+B(495)+B(509)+0.5&
                &*B(517)+0.8*B(547)+B(581)

  JVS(1307) = -B(130)

  JVS(1308) = -B(331)

  JVS(1309) = 0

  JVS(1310) = -B(460)

  JVS(1311) = -B(450)

  JVS(1312) = -B(592)

  JVS(1313) = B(50)-B(52)

  JVS(1314) = -B(470)

  JVS(1315) = -B(480)

  JVS(1316) = -0.9*B(251)

  JVS(1317) = -B(482)

  JVS(1318) = -B(153)

  JVS(1319) = -0.8*B(279)

  JVS(1320) = -0.8*B(275)

  JVS(1321) = 0

  JVS(1322) = 0

  JVS(1323) = 0.25*B(299)

  JVS(1324) = 0.89*B(8)

  JVS(1325) = -B(66)

  JVS(1326) = -B(68)+0.25*B(171)+0.25*B(300)

  JVS(1327) = -B(2)-B(3)-B(53)-B(67)-B(69)-B(84)-B(88)-B(131)-B(154)-0.9*B(252)-0.8*B(276)-0.8*B(280)-B(332)-B(451)&
                &-B(461)-B(471)-B(481)-B(483)-B(593)

  JVS(1328) = 0.25*B(172)

  JVS(1329) = -B(89)

  JVS(1330) = -B(85)

  JVS(1331) = 0

  JVS(1332) = 0.5*B(39)

  JVS(1333) = 0

  JVS(1334) = B(42)

  JVS(1335) = B(23)

  JVS(1336) = B(35)

  JVS(1337) = B(31)

  JVS(1338) = 0.5*B(175)

  JVS(1339) = 0.6*B(17)+B(177)

  JVS(1340) = 0.06*B(460)

  JVS(1341) = 0.06*B(450)

  JVS(1342) = 0

  JVS(1343) = B(231)

  JVS(1344) = 0

  JVS(1345) = 0.13*B(32)

  JVS(1346) = 0.06*B(470)

  JVS(1347) = 0.06*B(480)

  JVS(1348) = 0.06*B(482)

  JVS(1349) = B(24)

  JVS(1350) = B(163)+B(165)

  JVS(1351) = 0.65*B(46)+0.65*B(510)

  JVS(1352) = 0

  JVS(1353) = 0

  JVS(1354) = 0

  JVS(1355) = 0

  JVS(1356) = B(191)+0.3*B(195)

  JVS(1357) = 0.67*B(19)

  JVS(1358) = 0

  JVS(1359) = B(29)

  JVS(1360) = B(25)+B(530)+B(532)

  JVS(1361) = -B(548)

  JVS(1362) = 0.3*B(20)

  JVS(1363) = 0.53*B(281)+0.53*B(285)+0.26*B(289)-0.47*B(291)

  JVS(1364) = -B(271)

  JVS(1365) = B(295)+B(297)+B(301)+2*B(305)

  JVS(1366) = B(166)+0.53*B(286)+B(298)+B(533)

  JVS(1367) = B(164)+0.5*B(176)+0.65*B(511)+B(531)

  JVS(1368) = -B(171)

  JVS(1369) = 0.06*B(451)+0.06*B(461)+0.06*B(471)+0.06*B(481)+0.06*B(483)

  JVS(1370) = -B(167)-B(169)-B(172)-B(173)-2*B(178)-B(272)-0.47*B(292)-B(549)

  JVS(1371) = -B(170)

  JVS(1372) = -B(168)+B(192)+B(232)+0.53*B(282)+B(296)

  JVS(1373) = -B(174)+0.3*B(196)+0.26*B(290)+B(302)

  JVS(1374) = B(373)

  JVS(1375) = B(6)+B(94)

  JVS(1376) = B(209)

  JVS(1377) = 0.66*B(9)+B(103)+B(105)

  JVS(1378) = B(136)

  JVS(1379) = B(323)

  JVS(1380) = -B(329)

  JVS(1381) = B(367)

  JVS(1382) = B(609)

  JVS(1383) = B(339)

  JVS(1384) = B(227)

  JVS(1385) = 0.6*B(17)+B(177)

  JVS(1386) = B(404)

  JVS(1387) = B(7)

  JVS(1388) = B(584)

  JVS(1389) = -B(371)+B(374)

  JVS(1390) = B(18)+B(309)

  JVS(1391) = B(333)

  JVS(1392) = B(231)

  JVS(1393) = B(157)

  JVS(1394) = -B(345)

  JVS(1395) = B(378)

  JVS(1396) = B(414)

  JVS(1397) = B(382)-B(388)

  JVS(1398) = -B(86)

  JVS(1399) = B(181)

  JVS(1400) = B(386)-B(390)

  JVS(1401) = B(143)

  JVS(1402) = 1.206*B(524)+1.206*B(526)+0.206*B(528)

  JVS(1403) = 1.25*B(598)+1.25*B(600)

  JVS(1404) = 0

  JVS(1405) = 0

  JVS(1406) = 0

  JVS(1407) = 0

  JVS(1408) = B(578)

  JVS(1409) = 0.9*B(221)

  JVS(1410) = 0.9*B(504)

  JVS(1411) = 0.9*B(490)

  JVS(1412) = B(27)+0.4*B(534)+B(536)

  JVS(1413) = B(191)

  JVS(1414) = 1.9*B(512)+0.5*B(516)+2*B(518)

  JVS(1415) = 0

  JVS(1416) = 0

  JVS(1417) = B(540)+B(542)

  JVS(1418) = 0

  JVS(1419) = B(281)+B(285)

  JVS(1420) = 0

  JVS(1421) = 0.92*B(261)+B(263)

  JVS(1422) = B(295)+B(297)-B(306)

  JVS(1423) = 0.89*B(8)+B(90)-B(92)+2*B(99)+B(264)+B(286)+B(298)+2*B(519)+1.206*B(527)+B(537)+B(543)+1.25*B(601)

  JVS(1424) = -B(95)+B(104)+B(228)+0.4*B(535)

  JVS(1425) = B(82)+B(91)-B(101)+0.206*B(529)

  JVS(1426) = B(84)-B(88)

  JVS(1427) = B(167)-B(169)

  JVS(1428) = -B(5)-B(87)-B(89)-B(93)-B(96)-B(102)-B(109)-B(170)-B(307)-B(330)-B(346)-B(372)-B(389)-B(391)

  JVS(1429) = B(83)+B(85)+2*B(100)+B(114)+B(137)+B(144)+B(158)+B(168)+B(182)+B(192)+B(210)+0.9*B(222)+B(232)+0.92*B(262)&
                &+B(282)+B(296)+B(324)+B(334)+B(340)+B(368)+B(375)+B(379)+B(383)+B(387)+B(405)+B(415)+0.9*B(491)+0.9*B(505)&
                &+1.9*B(513)+1.206*B(525)+B(541)+B(579)+B(585)+1.25*B(599)+B(610)

  JVS(1430) = B(115)+0.5*B(517)

  JVS(1431) = B(49)

  JVS(1432) = 2*B(78)

  JVS(1433) = -B(209)

  JVS(1434) = -B(136)

  JVS(1435) = -B(323)

  JVS(1436) = -B(367)

  JVS(1437) = -B(609)

  JVS(1438) = -B(339)

  JVS(1439) = -B(404)

  JVS(1440) = -B(584)

  JVS(1441) = -B(374)

  JVS(1442) = 2*B(79)

  JVS(1443) = -B(333)

  JVS(1444) = -B(231)

  JVS(1445) = -B(157)

  JVS(1446) = -B(347)

  JVS(1447) = -B(378)

  JVS(1448) = -B(414)

  JVS(1449) = -B(382)

  JVS(1450) = B(86)

  JVS(1451) = -B(181)

  JVS(1452) = -B(386)

  JVS(1453) = -B(143)

  JVS(1454) = -B(524)

  JVS(1455) = -B(598)

  JVS(1456) = 0

  JVS(1457) = 0

  JVS(1458) = 0

  JVS(1459) = -B(578)

  JVS(1460) = -B(221)

  JVS(1461) = -B(504)

  JVS(1462) = -B(490)

  JVS(1463) = -B(191)

  JVS(1464) = -B(512)

  JVS(1465) = 0

  JVS(1466) = -B(540)

  JVS(1467) = 0

  JVS(1468) = -B(281)-B(283)

  JVS(1469) = 0

  JVS(1470) = -B(261)

  JVS(1471) = -B(295)

  JVS(1472) = 0.11*B(8)-B(99)

  JVS(1473) = 0

  JVS(1474) = -B(82)

  JVS(1475) = -B(84)

  JVS(1476) = -B(167)

  JVS(1477) = B(5)+B(87)+0.5*B(109)

  JVS(1478) = -B(83)-B(85)-B(100)-B(114)-B(137)-B(144)-B(158)-B(168)-B(182)-B(192)-B(210)-B(222)-B(232)-B(262)-B(282)&
                &-B(284)-B(296)-B(324)-B(334)-B(340)-B(348)-B(368)-B(375)-B(379)-B(383)-B(387)-B(405)-B(415)-B(491)-B(505)&
                &-B(513)-B(525)-B(541)-B(579)-B(585)-B(599)-B(610)

  JVS(1479) = -B(115)

  JVS(1480) = 0.7*B(120)

  JVS(1481) = B(520)

  JVS(1482) = B(16)

  JVS(1483) = 0.4*B(17)

  JVS(1484) = B(110)+0.75*B(112)

  JVS(1485) = 0.75*B(113)

  JVS(1486) = 0.33*B(380)+0.83*B(382)

  JVS(1487) = -B(185)

  JVS(1488) = 0.07*B(384)+0.17*B(386)

  JVS(1489) = -B(147)

  JVS(1490) = 0

  JVS(1491) = B(24)

  JVS(1492) = 0.31*B(153)

  JVS(1493) = B(14)

  JVS(1494) = 0

  JVS(1495) = -B(580)

  JVS(1496) = 0

  JVS(1497) = -B(508)

  JVS(1498) = -B(494)

  JVS(1499) = -B(195)

  JVS(1500) = -B(516)

  JVS(1501) = 0

  JVS(1502) = -B(546)+B(548)

  JVS(1503) = 0.3*B(20)

  JVS(1504) = -B(289)+B(291)

  JVS(1505) = 0

  JVS(1506) = -B(269)+B(271)

  JVS(1507) = -B(301)+B(303)

  JVS(1508) = 0

  JVS(1509) = B(111)+0.7*B(121)+B(521)

  JVS(1510) = -B(118)+0.33*B(381)+0.07*B(385)

  JVS(1511) = 0.31*B(154)

  JVS(1512) = B(167)-0.1*B(173)+2*B(178)+B(272)+B(292)+B(304)+B(549)

  JVS(1513) = 0

  JVS(1514) = -B(114)+B(168)+0.83*B(383)+0.17*B(387)

  JVS(1515) = -B(115)-2*B(116)-2*B(117)-B(119)-B(148)-0.1*B(174)-B(186)-B(196)-B(270)-B(290)-B(302)-B(495)-B(509)-B(517)&
                &-B(547)-B(581)
      
END SUBROUTINE mozart_mosaic_4bin_aq_Jac_SP














SUBROUTINE mozart_mosaic_4bin_aq_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(1515), W(140), a
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
      
END SUBROUTINE mozart_mosaic_4bin_aq_KppDecomp



SUBROUTINE mozart_mosaic_4bin_aq_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(1515), W(140), a
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
      
END SUBROUTINE mozart_mosaic_4bin_aq_KppDecompCmplx


SUBROUTINE mozart_mosaic_4bin_aq_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(1515), X(140), sum

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
      
END SUBROUTINE mozart_mosaic_4bin_aq_KppSolveIndirect


SUBROUTINE mozart_mosaic_4bin_aq_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(1515), X(140), sum

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
      
END SUBROUTINE mozart_mosaic_4bin_aq_KppSolveCmplx













SUBROUTINE mozart_mosaic_4bin_aq_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(25) = X(25)-JVS(165)*X(24)
  X(26) = X(26)-JVS(168)*X(24)
  X(28) = X(28)-JVS(173)*X(27)
  X(35) = X(35)-JVS(199)*X(27)
  X(42) = X(42)-JVS(224)*X(41)
  X(47) = X(47)-JVS(241)*X(25)
  X(51) = X(51)-JVS(260)*X(18)
  X(60) = X(60)-JVS(304)*X(59)
  X(62) = X(62)-JVS(313)*X(26)-JVS(314)*X(28)-JVS(315)*X(31)
  X(63) = X(63)-JVS(320)*X(26)-JVS(321)*X(28)-JVS(322)*X(42)
  X(64) = X(64)-JVS(329)*X(27)-JVS(330)*X(36)-JVS(331)*X(41)
  X(70) = X(70)-JVS(363)*X(24)-JVS(364)*X(33)
  X(74) = X(74)-JVS(384)*X(38)-JVS(385)*X(42)
  X(75) = X(75)-JVS(390)*X(40)-JVS(391)*X(43)
  X(78) = X(78)-JVS(409)*X(21)-JVS(410)*X(37)-JVS(411)*X(64)
  X(82) = X(82)-JVS(429)*X(23)-JVS(430)*X(59)-JVS(431)*X(65)
  X(84) = X(84)-JVS(448)*X(19)-JVS(449)*X(55)-JVS(450)*X(77)
  X(85) = X(85)-JVS(456)*X(32)-JVS(457)*X(63)-JVS(458)*X(78)
  X(89) = X(89)-JVS(475)*X(52)-JVS(476)*X(66)-JVS(477)*X(88)
  X(90) = X(90)-JVS(483)*X(72)
  X(91) = X(91)-JVS(489)*X(47)-JVS(490)*X(70)
  X(93) = X(93)-JVS(502)*X(27)-JVS(503)*X(61)
  X(94) = X(94)-JVS(508)*X(39)-JVS(509)*X(41)
  X(95) = X(95)-JVS(514)*X(34)-JVS(515)*X(93)-JVS(516)*X(94)
  X(96) = X(96)-JVS(522)*X(84)
  X(97) = X(97)-JVS(530)*X(17)-JVS(531)*X(49)
  X(98) = X(98)-JVS(537)*X(46)-JVS(538)*X(93)-JVS(539)*X(94)
  X(99) = X(99)-JVS(545)*X(23)-JVS(546)*X(61)-JVS(547)*X(62)-JVS(548)*X(70)-JVS(549)*X(74)-JVS(550)*X(91)-JVS(551)*X(93)&
            &-JVS(552)*X(94)-JVS(553)*X(98)
  X(101) = X(101)-JVS(578)*X(61)-JVS(579)*X(86)-JVS(580)*X(87)-JVS(581)*X(93)
  X(102) = X(102)-JVS(590)*X(16)-JVS(591)*X(48)-JVS(592)*X(66)-JVS(593)*X(88)
  X(105) = X(105)-JVS(609)*X(92)
  X(106) = X(106)-JVS(616)*X(79)-JVS(617)*X(92)-JVS(618)*X(105)
  X(110) = X(110)-JVS(639)*X(49)-JVS(640)*X(51)-JVS(641)*X(69)-JVS(642)*X(76)-JVS(643)*X(79)-JVS(644)*X(83)-JVS(645)&
             &*X(86)-JVS(646)*X(87)-JVS(647)*X(88)-JVS(648)*X(92)-JVS(649)*X(97)-JVS(650)*X(104)-JVS(651)*X(105)-JVS(652)&
             &*X(106)-JVS(653)*X(107)-JVS(654)*X(109)
  X(111) = X(111)-JVS(667)*X(23)-JVS(668)*X(25)-JVS(669)*X(35)-JVS(670)*X(41)-JVS(671)*X(44)-JVS(672)*X(46)-JVS(673)&
             &*X(59)-JVS(674)*X(83)-JVS(675)*X(86)-JVS(676)*X(87)-JVS(677)*X(91)-JVS(678)*X(92)-JVS(679)*X(93)-JVS(680)&
             &*X(94)-JVS(681)*X(95)-JVS(682)*X(98)-JVS(683)*X(99)-JVS(684)*X(101)-JVS(685)*X(104)-JVS(686)*X(106)-JVS(687)&
             &*X(107)-JVS(688)*X(108)-JVS(689)*X(109)
  X(112) = X(112)-JVS(714)*X(30)-JVS(715)*X(60)-JVS(716)*X(79)-JVS(717)*X(83)
  X(113) = X(113)-JVS(731)*X(108)
  X(114) = X(114)-JVS(737)*X(45)-JVS(738)*X(48)-JVS(739)*X(51)-JVS(740)*X(52)-JVS(741)*X(72)-JVS(742)*X(88)-JVS(743)&
             &*X(89)-JVS(744)*X(90)-JVS(745)*X(97)-JVS(746)*X(102)-JVS(747)*X(113)
  X(115) = X(115)-JVS(757)*X(83)-JVS(758)*X(86)-JVS(759)*X(87)-JVS(760)*X(104)-JVS(761)*X(107)-JVS(762)*X(109)
  X(116) = X(116)-JVS(772)*X(79)-JVS(773)*X(92)
  X(117) = X(117)-JVS(781)*X(58)-JVS(782)*X(86)-JVS(783)*X(87)-JVS(784)*X(88)-JVS(785)*X(104)-JVS(786)*X(107)-JVS(787)&
             &*X(109)
  X(118) = X(118)-JVS(794)*X(83)
  X(119) = X(119)-JVS(802)*X(69)-JVS(803)*X(86)-JVS(804)*X(87)-JVS(805)*X(104)-JVS(806)*X(107)-JVS(807)*X(109)
  X(120) = X(120)-JVS(815)*X(103)-JVS(816)*X(108)-JVS(817)*X(118)-JVS(818)*X(119)
  X(121) = X(121)-JVS(830)*X(50)-JVS(831)*X(86)-JVS(832)*X(87)-JVS(833)*X(104)-JVS(834)*X(107)-JVS(835)*X(109)-JVS(836)&
             &*X(110)-JVS(837)*X(115)-JVS(838)*X(116)-JVS(839)*X(117)-JVS(840)*X(118)-JVS(841)*X(119)
  X(122) = X(122)-JVS(850)*X(104)-JVS(851)*X(107)-JVS(852)*X(109)
  X(123) = X(123)-JVS(861)*X(80)-JVS(862)*X(103)-JVS(863)*X(108)
  X(124) = X(124)-JVS(873)*X(69)-JVS(874)*X(86)-JVS(875)*X(87)-JVS(876)*X(104)-JVS(877)*X(107)-JVS(878)*X(109)-JVS(879)&
             &*X(119)-JVS(880)*X(122)
  X(125) = X(125)-JVS(888)*X(72)-JVS(889)*X(81)-JVS(890)*X(90)-JVS(891)*X(113)-JVS(892)*X(121)-JVS(893)*X(124)
  X(126) = X(126)-JVS(908)*X(44)-JVS(909)*X(61)-JVS(910)*X(71)-JVS(911)*X(74)-JVS(912)*X(93)-JVS(913)*X(94)-JVS(914)&
             &*X(95)-JVS(915)*X(98)-JVS(916)*X(101)-JVS(917)*X(104)-JVS(918)*X(105)-JVS(919)*X(107)-JVS(920)*X(109)-JVS(921)&
             &*X(113)-JVS(922)*X(117)-JVS(923)*X(121)-JVS(924)*X(123)-JVS(925)*X(124)-JVS(926)*X(125)
  X(127) = X(127)-JVS(941)*X(29)-JVS(942)*X(57)-JVS(943)*X(80)-JVS(944)*X(120)-JVS(945)*X(122)-JVS(946)*X(124)
  X(128) = X(128)-JVS(957)*X(80)-JVS(958)*X(103)-JVS(959)*X(108)
  X(129) = X(129)-JVS(969)*X(54)-JVS(970)*X(123)-JVS(971)*X(128)
  X(130) = X(130)-JVS(981)*X(30)-JVS(982)*X(50)-JVS(983)*X(51)-JVS(984)*X(53)-JVS(985)*X(59)-JVS(986)*X(60)-JVS(987)&
             &*X(65)-JVS(988)*X(68)-JVS(989)*X(69)-JVS(990)*X(72)-JVS(991)*X(73)-JVS(992)*X(77)-JVS(993)*X(79)-JVS(994)&
             &*X(80)-JVS(995)*X(81)-JVS(996)*X(83)-JVS(997)*X(84)-JVS(998)*X(86)-JVS(999)*X(87)-JVS(1000)*X(88)-JVS(1001)&
             &*X(90)-JVS(1002)*X(92)-JVS(1003)*X(97)-JVS(1004)*X(100)-JVS(1005)*X(102)-JVS(1006)*X(103)-JVS(1007)*X(104)&
             &-JVS(1008)*X(105)-JVS(1009)*X(107)-JVS(1010)*X(108)-JVS(1011)*X(109)-JVS(1012)*X(112)-JVS(1013)*X(113)&
             &-JVS(1014)*X(115)-JVS(1015)*X(116)-JVS(1016)*X(117)-JVS(1017)*X(118)-JVS(1018)*X(119)-JVS(1019)*X(120)&
             &-JVS(1020)*X(121)-JVS(1021)*X(122)-JVS(1022)*X(123)-JVS(1023)*X(124)-JVS(1024)*X(125)-JVS(1025)*X(127)&
             &-JVS(1026)*X(128)-JVS(1027)*X(129)
  X(131) = X(131)-JVS(1039)*X(80)-JVS(1040)*X(108)
  X(132) = X(132)-JVS(1049)*X(54)-JVS(1050)*X(81)-JVS(1051)*X(108)-JVS(1052)*X(123)-JVS(1053)*X(129)-JVS(1054)*X(131)
  X(133) = X(133)-JVS(1064)*X(40)-JVS(1065)*X(43)-JVS(1066)*X(56)-JVS(1067)*X(73)-JVS(1068)*X(75)-JVS(1069)*X(81)&
             &-JVS(1070)*X(92)-JVS(1071)*X(103)-JVS(1072)*X(104)-JVS(1073)*X(105)-JVS(1074)*X(107)-JVS(1075)*X(108)&
             &-JVS(1076)*X(109)-JVS(1077)*X(113)-JVS(1078)*X(114)-JVS(1079)*X(117)-JVS(1080)*X(120)-JVS(1081)*X(122)&
             &-JVS(1082)*X(124)-JVS(1083)*X(126)-JVS(1084)*X(127)-JVS(1085)*X(128)-JVS(1086)*X(129)-JVS(1087)*X(130)&
             &-JVS(1088)*X(131)-JVS(1089)*X(132)
  X(134) = X(134)-JVS(1098)*X(1)-JVS(1099)*X(15)-JVS(1100)*X(16)-JVS(1101)*X(17)-JVS(1102)*X(18)-JVS(1103)*X(20)&
             &-JVS(1104)*X(22)-JVS(1105)*X(23)-JVS(1106)*X(24)-JVS(1107)*X(26)-JVS(1108)*X(27)-JVS(1109)*X(28)-JVS(1110)&
             &*X(29)-JVS(1111)*X(31)-JVS(1112)*X(32)-JVS(1113)*X(33)-JVS(1114)*X(36)-JVS(1115)*X(37)-JVS(1116)*X(38)&
             &-JVS(1117)*X(39)-JVS(1118)*X(41)-JVS(1119)*X(42)-JVS(1120)*X(43)-JVS(1121)*X(45)-JVS(1122)*X(48)-JVS(1123)&
             &*X(49)-JVS(1124)*X(50)-JVS(1125)*X(52)-JVS(1126)*X(53)-JVS(1127)*X(54)-JVS(1128)*X(55)-JVS(1129)*X(56)&
             &-JVS(1130)*X(57)-JVS(1131)*X(58)-JVS(1132)*X(59)-JVS(1133)*X(61)-JVS(1134)*X(62)-JVS(1135)*X(64)-JVS(1136)&
             &*X(66)-JVS(1137)*X(67)-JVS(1138)*X(68)-JVS(1139)*X(69)-JVS(1140)*X(70)-JVS(1141)*X(71)-JVS(1142)*X(72)&
             &-JVS(1143)*X(73)-JVS(1144)*X(74)-JVS(1145)*X(75)-JVS(1146)*X(76)-JVS(1147)*X(77)-JVS(1148)*X(78)-JVS(1149)&
             &*X(79)-JVS(1150)*X(80)-JVS(1151)*X(81)-JVS(1152)*X(82)-JVS(1153)*X(83)-JVS(1154)*X(84)-JVS(1155)*X(85)&
             &-JVS(1156)*X(86)-JVS(1157)*X(87)-JVS(1158)*X(88)-JVS(1159)*X(89)-JVS(1160)*X(90)-JVS(1161)*X(92)-JVS(1162)&
             &*X(93)-JVS(1163)*X(94)-JVS(1164)*X(95)-JVS(1165)*X(96)-JVS(1166)*X(97)-JVS(1167)*X(98)-JVS(1168)*X(99)&
             &-JVS(1169)*X(100)-JVS(1170)*X(101)-JVS(1171)*X(102)-JVS(1172)*X(103)-JVS(1173)*X(104)-JVS(1174)*X(105)&
             &-JVS(1175)*X(106)-JVS(1176)*X(107)-JVS(1177)*X(108)-JVS(1178)*X(109)-JVS(1179)*X(110)-JVS(1180)*X(111)&
             &-JVS(1181)*X(112)-JVS(1182)*X(113)-JVS(1183)*X(114)-JVS(1184)*X(115)-JVS(1185)*X(116)-JVS(1186)*X(117)&
             &-JVS(1187)*X(118)-JVS(1188)*X(119)-JVS(1189)*X(120)-JVS(1190)*X(121)-JVS(1191)*X(122)-JVS(1192)*X(123)&
             &-JVS(1193)*X(124)-JVS(1194)*X(125)-JVS(1195)*X(126)-JVS(1196)*X(127)-JVS(1197)*X(128)-JVS(1198)*X(129)&
             &-JVS(1199)*X(130)-JVS(1200)*X(131)-JVS(1201)*X(132)-JVS(1202)*X(133)
  X(135) = X(135)-JVS(1210)*X(22)-JVS(1211)*X(23)-JVS(1212)*X(24)-JVS(1213)*X(25)-JVS(1214)*X(26)-JVS(1215)*X(27)&
             &-JVS(1216)*X(28)-JVS(1217)*X(30)-JVS(1218)*X(34)-JVS(1219)*X(35)-JVS(1220)*X(41)-JVS(1221)*X(42)-JVS(1222)&
             &*X(43)-JVS(1223)*X(44)-JVS(1224)*X(45)-JVS(1225)*X(46)-JVS(1226)*X(47)-JVS(1227)*X(48)-JVS(1228)*X(49)&
             &-JVS(1229)*X(51)-JVS(1230)*X(53)-JVS(1231)*X(54)-JVS(1232)*X(55)-JVS(1233)*X(56)-JVS(1234)*X(59)-JVS(1235)&
             &*X(60)-JVS(1236)*X(62)-JVS(1237)*X(64)-JVS(1238)*X(65)-JVS(1239)*X(69)-JVS(1240)*X(70)-JVS(1241)*X(72)&
             &-JVS(1242)*X(74)-JVS(1243)*X(76)-JVS(1244)*X(77)-JVS(1245)*X(78)-JVS(1246)*X(79)-JVS(1247)*X(80)-JVS(1248)&
             &*X(81)-JVS(1249)*X(82)-JVS(1250)*X(83)-JVS(1251)*X(84)-JVS(1252)*X(85)-JVS(1253)*X(86)-JVS(1254)*X(87)&
             &-JVS(1255)*X(88)-JVS(1256)*X(89)-JVS(1257)*X(90)-JVS(1258)*X(91)-JVS(1259)*X(92)-JVS(1260)*X(93)-JVS(1261)&
             &*X(94)-JVS(1262)*X(95)-JVS(1263)*X(96)-JVS(1264)*X(97)-JVS(1265)*X(98)-JVS(1266)*X(99)-JVS(1267)*X(100)&
             &-JVS(1268)*X(101)-JVS(1269)*X(102)-JVS(1270)*X(103)-JVS(1271)*X(104)-JVS(1272)*X(105)-JVS(1273)*X(106)&
             &-JVS(1274)*X(107)-JVS(1275)*X(108)-JVS(1276)*X(109)-JVS(1277)*X(111)-JVS(1278)*X(112)-JVS(1279)*X(113)&
             &-JVS(1280)*X(114)-JVS(1281)*X(115)-JVS(1282)*X(116)-JVS(1283)*X(117)-JVS(1284)*X(118)-JVS(1285)*X(119)&
             &-JVS(1286)*X(120)-JVS(1287)*X(121)-JVS(1288)*X(122)-JVS(1289)*X(123)-JVS(1290)*X(124)-JVS(1291)*X(125)&
             &-JVS(1292)*X(126)-JVS(1293)*X(127)-JVS(1294)*X(128)-JVS(1295)*X(129)-JVS(1296)*X(130)-JVS(1297)*X(131)&
             &-JVS(1298)*X(132)-JVS(1299)*X(133)-JVS(1300)*X(134)
  X(136) = X(136)-JVS(1307)*X(59)-JVS(1308)*X(63)-JVS(1309)*X(85)-JVS(1310)*X(86)-JVS(1311)*X(87)-JVS(1312)*X(92)&
             &-JVS(1313)*X(96)-JVS(1314)*X(104)-JVS(1315)*X(107)-JVS(1316)*X(108)-JVS(1317)*X(109)-JVS(1318)*X(113)&
             &-JVS(1319)*X(123)-JVS(1320)*X(128)-JVS(1321)*X(130)-JVS(1322)*X(131)-JVS(1323)*X(132)-JVS(1324)*X(133)&
             &-JVS(1325)*X(134)-JVS(1326)*X(135)
  X(137) = X(137)-JVS(1332)*X(35)-JVS(1333)*X(41)-JVS(1334)*X(44)-JVS(1335)*X(50)-JVS(1336)*X(52)-JVS(1337)*X(66)&
             &-JVS(1338)*X(68)-JVS(1339)*X(73)-JVS(1340)*X(86)-JVS(1341)*X(87)-JVS(1342)*X(88)-JVS(1343)*X(89)-JVS(1344)&
             &*X(94)-JVS(1345)*X(101)-JVS(1346)*X(104)-JVS(1347)*X(107)-JVS(1348)*X(109)-JVS(1349)*X(110)-JVS(1350)*X(114)&
             &-JVS(1351)*X(115)-JVS(1352)*X(116)-JVS(1353)*X(117)-JVS(1354)*X(118)-JVS(1355)*X(119)-JVS(1356)*X(121)&
             &-JVS(1357)*X(123)-JVS(1358)*X(124)-JVS(1359)*X(125)-JVS(1360)*X(126)-JVS(1361)*X(127)-JVS(1362)*X(128)&
             &-JVS(1363)*X(129)-JVS(1364)*X(131)-JVS(1365)*X(132)-JVS(1366)*X(133)-JVS(1367)*X(134)-JVS(1368)*X(135)&
             &-JVS(1369)*X(136)
  X(138) = X(138)-JVS(1374)*X(21)-JVS(1375)*X(40)-JVS(1376)*X(51)-JVS(1377)*X(56)-JVS(1378)*X(60)-JVS(1379)*X(62)&
             &-JVS(1380)*X(63)-JVS(1381)*X(64)-JVS(1382)*X(65)-JVS(1383)*X(70)-JVS(1384)*X(71)-JVS(1385)*X(73)-JVS(1386)&
             &*X(74)-JVS(1387)*X(75)-JVS(1388)*X(76)-JVS(1389)*X(78)-JVS(1390)*X(81)-JVS(1391)*X(85)-JVS(1392)*X(89)&
             &-JVS(1393)*X(90)-JVS(1394)*X(91)-JVS(1395)*X(93)-JVS(1396)*X(94)-JVS(1397)*X(95)-JVS(1398)*X(96)-JVS(1399)&
             &*X(97)-JVS(1400)*X(98)-JVS(1401)*X(102)-JVS(1402)*X(103)-JVS(1403)*X(105)-JVS(1404)*X(106)-JVS(1405)*X(108)&
             &-JVS(1406)*X(113)-JVS(1407)*X(114)-JVS(1408)*X(116)-JVS(1409)*X(117)-JVS(1410)*X(118)-JVS(1411)*X(119)&
             &-JVS(1412)*X(120)-JVS(1413)*X(121)-JVS(1414)*X(122)-JVS(1415)*X(124)-JVS(1416)*X(126)-JVS(1417)*X(127)&
             &-JVS(1418)*X(128)-JVS(1419)*X(129)-JVS(1420)*X(130)-JVS(1421)*X(131)-JVS(1422)*X(132)-JVS(1423)*X(133)&
             &-JVS(1424)*X(134)-JVS(1425)*X(135)-JVS(1426)*X(136)-JVS(1427)*X(137)
  X(139) = X(139)-JVS(1431)*X(1)-JVS(1432)*X(19)-JVS(1433)*X(51)-JVS(1434)*X(60)-JVS(1435)*X(62)-JVS(1436)*X(64)&
             &-JVS(1437)*X(65)-JVS(1438)*X(70)-JVS(1439)*X(74)-JVS(1440)*X(76)-JVS(1441)*X(78)-JVS(1442)*X(84)-JVS(1443)&
             &*X(85)-JVS(1444)*X(89)-JVS(1445)*X(90)-JVS(1446)*X(91)-JVS(1447)*X(93)-JVS(1448)*X(94)-JVS(1449)*X(95)&
             &-JVS(1450)*X(96)-JVS(1451)*X(97)-JVS(1452)*X(98)-JVS(1453)*X(102)-JVS(1454)*X(103)-JVS(1455)*X(105)-JVS(1456)&
             &*X(106)-JVS(1457)*X(108)-JVS(1458)*X(113)-JVS(1459)*X(116)-JVS(1460)*X(117)-JVS(1461)*X(118)-JVS(1462)*X(119)&
             &-JVS(1463)*X(121)-JVS(1464)*X(122)-JVS(1465)*X(124)-JVS(1466)*X(127)-JVS(1467)*X(128)-JVS(1468)*X(129)&
             &-JVS(1469)*X(130)-JVS(1470)*X(131)-JVS(1471)*X(132)-JVS(1472)*X(133)-JVS(1473)*X(134)-JVS(1474)*X(135)&
             &-JVS(1475)*X(136)-JVS(1476)*X(137)-JVS(1477)*X(138)
  X(140) = X(140)-JVS(1480)*X(53)-JVS(1481)*X(67)-JVS(1482)*X(68)-JVS(1483)*X(73)-JVS(1484)*X(77)-JVS(1485)*X(84)&
             &-JVS(1486)*X(95)-JVS(1487)*X(97)-JVS(1488)*X(98)-JVS(1489)*X(102)-JVS(1490)*X(108)-JVS(1491)*X(110)-JVS(1492)&
             &*X(113)-JVS(1493)*X(114)-JVS(1494)*X(115)-JVS(1495)*X(116)-JVS(1496)*X(117)-JVS(1497)*X(118)-JVS(1498)*X(119)&
             &-JVS(1499)*X(121)-JVS(1500)*X(122)-JVS(1501)*X(124)-JVS(1502)*X(127)-JVS(1503)*X(128)-JVS(1504)*X(129)&
             &-JVS(1505)*X(130)-JVS(1506)*X(131)-JVS(1507)*X(132)-JVS(1508)*X(133)-JVS(1509)*X(134)-JVS(1510)*X(135)&
             &-JVS(1511)*X(136)-JVS(1512)*X(137)-JVS(1513)*X(138)-JVS(1514)*X(139)
  X(140) = X(140)/JVS(1515)
  X(139) = (X(139)-JVS(1479)*X(140))/(JVS(1478))
  X(138) = (X(138)-JVS(1429)*X(139)-JVS(1430)*X(140))/(JVS(1428))
  X(137) = (X(137)-JVS(1371)*X(138)-JVS(1372)*X(139)-JVS(1373)*X(140))/(JVS(1370))
  X(136) = (X(136)-JVS(1328)*X(137)-JVS(1329)*X(138)-JVS(1330)*X(139)-JVS(1331)*X(140))/(JVS(1327))
  X(135) = (X(135)-JVS(1302)*X(136)-JVS(1303)*X(137)-JVS(1304)*X(138)-JVS(1305)*X(139)-JVS(1306)*X(140))/(JVS(1301))
  X(134) = (X(134)-JVS(1204)*X(135)-JVS(1205)*X(136)-JVS(1206)*X(137)-JVS(1207)*X(138)-JVS(1208)*X(139)-JVS(1209)&
             &*X(140))/(JVS(1203))
  X(133) = (X(133)-JVS(1091)*X(134)-JVS(1092)*X(135)-JVS(1093)*X(136)-JVS(1094)*X(137)-JVS(1095)*X(138)-JVS(1096)*X(139)&
             &-JVS(1097)*X(140))/(JVS(1090))
  X(132) = (X(132)-JVS(1056)*X(133)-JVS(1057)*X(134)-JVS(1058)*X(135)-JVS(1059)*X(136)-JVS(1060)*X(137)-JVS(1061)*X(138)&
             &-JVS(1062)*X(139)-JVS(1063)*X(140))/(JVS(1055))
  X(131) = (X(131)-JVS(1042)*X(133)-JVS(1043)*X(134)-JVS(1044)*X(135)-JVS(1045)*X(136)-JVS(1046)*X(137)-JVS(1047)*X(139)&
             &-JVS(1048)*X(140))/(JVS(1041))
  X(130) = (X(130)-JVS(1029)*X(131)-JVS(1030)*X(132)-JVS(1031)*X(133)-JVS(1032)*X(134)-JVS(1033)*X(135)-JVS(1034)*X(136)&
             &-JVS(1035)*X(137)-JVS(1036)*X(138)-JVS(1037)*X(139)-JVS(1038)*X(140))/(JVS(1028))
  X(129) = (X(129)-JVS(973)*X(131)-JVS(974)*X(133)-JVS(975)*X(134)-JVS(976)*X(135)-JVS(977)*X(136)-JVS(978)*X(137)&
             &-JVS(979)*X(139)-JVS(980)*X(140))/(JVS(972))
  X(128) = (X(128)-JVS(961)*X(131)-JVS(962)*X(133)-JVS(963)*X(134)-JVS(964)*X(135)-JVS(965)*X(136)-JVS(966)*X(137)&
             &-JVS(967)*X(139)-JVS(968)*X(140))/(JVS(960))
  X(127) = (X(127)-JVS(948)*X(129)-JVS(949)*X(131)-JVS(950)*X(133)-JVS(951)*X(134)-JVS(952)*X(135)-JVS(953)*X(136)&
             &-JVS(954)*X(137)-JVS(955)*X(139)-JVS(956)*X(140))/(JVS(947))
  X(126) = (X(126)-JVS(928)*X(127)-JVS(929)*X(128)-JVS(930)*X(129)-JVS(931)*X(131)-JVS(932)*X(132)-JVS(933)*X(133)&
             &-JVS(934)*X(134)-JVS(935)*X(135)-JVS(936)*X(136)-JVS(937)*X(137)-JVS(938)*X(138)-JVS(939)*X(139)-JVS(940)&
             &*X(140))/(JVS(927))
  X(125) = (X(125)-JVS(895)*X(127)-JVS(896)*X(128)-JVS(897)*X(129)-JVS(898)*X(131)-JVS(899)*X(132)-JVS(900)*X(133)&
             &-JVS(901)*X(134)-JVS(902)*X(135)-JVS(903)*X(136)-JVS(904)*X(137)-JVS(905)*X(138)-JVS(906)*X(139)-JVS(907)&
             &*X(140))/(JVS(894))
  X(124) = (X(124)-JVS(882)*X(133)-JVS(883)*X(134)-JVS(884)*X(135)-JVS(885)*X(136)-JVS(886)*X(139)-JVS(887)*X(140))&
             &/(JVS(881))
  X(123) = (X(123)-JVS(865)*X(131)-JVS(866)*X(133)-JVS(867)*X(134)-JVS(868)*X(135)-JVS(869)*X(136)-JVS(870)*X(137)&
             &-JVS(871)*X(139)-JVS(872)*X(140))/(JVS(864))
  X(122) = (X(122)-JVS(854)*X(124)-JVS(855)*X(133)-JVS(856)*X(134)-JVS(857)*X(135)-JVS(858)*X(136)-JVS(859)*X(139)&
             &-JVS(860)*X(140))/(JVS(853))
  X(121) = (X(121)-JVS(843)*X(124)-JVS(844)*X(133)-JVS(845)*X(134)-JVS(846)*X(135)-JVS(847)*X(136)-JVS(848)*X(139)&
             &-JVS(849)*X(140))/(JVS(842))
  X(120) = (X(120)-JVS(820)*X(122)-JVS(821)*X(124)-JVS(822)*X(129)-JVS(823)*X(131)-JVS(824)*X(133)-JVS(825)*X(134)&
             &-JVS(826)*X(135)-JVS(827)*X(136)-JVS(828)*X(139)-JVS(829)*X(140))/(JVS(819))
  X(119) = (X(119)-JVS(809)*X(133)-JVS(810)*X(134)-JVS(811)*X(135)-JVS(812)*X(136)-JVS(813)*X(139)-JVS(814)*X(140))&
             &/(JVS(808))
  X(118) = (X(118)-JVS(796)*X(124)-JVS(797)*X(133)-JVS(798)*X(134)-JVS(799)*X(135)-JVS(800)*X(139)-JVS(801)*X(140))&
             &/(JVS(795))
  X(117) = (X(117)-JVS(789)*X(133)-JVS(790)*X(134)-JVS(791)*X(135)-JVS(792)*X(136)-JVS(793)*X(139))/(JVS(788))
  X(116) = (X(116)-JVS(775)*X(133)-JVS(776)*X(134)-JVS(777)*X(135)-JVS(778)*X(136)-JVS(779)*X(139)-JVS(780)*X(140))&
             &/(JVS(774))
  X(115) = (X(115)-JVS(764)*X(118)-JVS(765)*X(124)-JVS(766)*X(133)-JVS(767)*X(134)-JVS(768)*X(135)-JVS(769)*X(136)&
             &-JVS(770)*X(139)-JVS(771)*X(140))/(JVS(763))
  X(114) = (X(114)-JVS(749)*X(117)-JVS(750)*X(128)-JVS(751)*X(133)-JVS(752)*X(134)-JVS(753)*X(135)-JVS(754)*X(136)&
             &-JVS(755)*X(139)-JVS(756)*X(140))/(JVS(748))
  X(113) = (X(113)-JVS(733)*X(128)-JVS(734)*X(133)-JVS(735)*X(134)-JVS(736)*X(136))/(JVS(732))
  X(112) = (X(112)-JVS(719)*X(116)-JVS(720)*X(118)-JVS(721)*X(127)-JVS(722)*X(129)-JVS(723)*X(131)-JVS(724)*X(133)&
             &-JVS(725)*X(134)-JVS(726)*X(135)-JVS(727)*X(136)-JVS(728)*X(137)-JVS(729)*X(139)-JVS(730)*X(140))/(JVS(718))
  X(111) = (X(111)-JVS(691)*X(112)-JVS(692)*X(113)-JVS(693)*X(114)-JVS(694)*X(115)-JVS(695)*X(116)-JVS(696)*X(118)&
             &-JVS(697)*X(120)-JVS(698)*X(123)-JVS(699)*X(124)-JVS(700)*X(126)-JVS(701)*X(127)-JVS(702)*X(128)-JVS(703)&
             &*X(129)-JVS(704)*X(130)-JVS(705)*X(131)-JVS(706)*X(133)-JVS(707)*X(134)-JVS(708)*X(135)-JVS(709)*X(136)&
             &-JVS(710)*X(137)-JVS(711)*X(138)-JVS(712)*X(139)-JVS(713)*X(140))/(JVS(690))
  X(110) = (X(110)-JVS(656)*X(115)-JVS(657)*X(116)-JVS(658)*X(117)-JVS(659)*X(118)-JVS(660)*X(119)-JVS(661)*X(133)&
             &-JVS(662)*X(134)-JVS(663)*X(135)-JVS(664)*X(136)-JVS(665)*X(139)-JVS(666)*X(140))/(JVS(655))
  X(109) = (X(109)-JVS(636)*X(133)-JVS(637)*X(134)-JVS(638)*X(136))/(JVS(635))
  X(108) = (X(108)-JVS(632)*X(133)-JVS(633)*X(134)-JVS(634)*X(136))/(JVS(631))
  X(107) = (X(107)-JVS(628)*X(133)-JVS(629)*X(134)-JVS(630)*X(136))/(JVS(627))
  X(106) = (X(106)-JVS(620)*X(116)-JVS(621)*X(133)-JVS(622)*X(134)-JVS(623)*X(135)-JVS(624)*X(136)-JVS(625)*X(139)&
             &-JVS(626)*X(140))/(JVS(619))
  X(105) = (X(105)-JVS(611)*X(133)-JVS(612)*X(134)-JVS(613)*X(135)-JVS(614)*X(136)-JVS(615)*X(139))/(JVS(610))
  X(104) = (X(104)-JVS(606)*X(133)-JVS(607)*X(134)-JVS(608)*X(136))/(JVS(605))
  X(103) = (X(103)-JVS(601)*X(108)-JVS(602)*X(133)-JVS(603)*X(135)-JVS(604)*X(139))/(JVS(600))
  X(102) = (X(102)-JVS(595)*X(117)-JVS(596)*X(134)-JVS(597)*X(135)-JVS(598)*X(139)-JVS(599)*X(140))/(JVS(594))
  X(101) = (X(101)-JVS(583)*X(104)-JVS(584)*X(107)-JVS(585)*X(109)-JVS(586)*X(134)-JVS(587)*X(135)-JVS(588)*X(136)&
             &-JVS(589)*X(139))/(JVS(582))
  X(100) = (X(100)-JVS(567)*X(102)-JVS(568)*X(116)-JVS(569)*X(118)-JVS(570)*X(119)-JVS(571)*X(121)-JVS(572)*X(122)&
             &-JVS(573)*X(127)-JVS(574)*X(129)-JVS(575)*X(131)-JVS(576)*X(134)-JVS(577)*X(140))/(JVS(566))
  X(99) = (X(99)-JVS(555)*X(101)-JVS(556)*X(112)-JVS(557)*X(127)-JVS(558)*X(131)-JVS(559)*X(133)-JVS(560)*X(134)&
            &-JVS(561)*X(135)-JVS(562)*X(137)-JVS(563)*X(138)-JVS(564)*X(139)-JVS(565)*X(140))/(JVS(554))
  X(98) = (X(98)-JVS(541)*X(134)-JVS(542)*X(135)-JVS(543)*X(138)-JVS(544)*X(139))/(JVS(540))
  X(97) = (X(97)-JVS(533)*X(134)-JVS(534)*X(135)-JVS(535)*X(139)-JVS(536)*X(140))/(JVS(532))
  X(96) = (X(96)-JVS(524)*X(113)-JVS(525)*X(130)-JVS(526)*X(134)-JVS(527)*X(135)-JVS(528)*X(136)-JVS(529)*X(138))&
            &/(JVS(523))
  X(95) = (X(95)-JVS(518)*X(134)-JVS(519)*X(135)-JVS(520)*X(138)-JVS(521)*X(139))/(JVS(517))
  X(94) = (X(94)-JVS(511)*X(134)-JVS(512)*X(135)-JVS(513)*X(139))/(JVS(510))
  X(93) = (X(93)-JVS(505)*X(134)-JVS(506)*X(135)-JVS(507)*X(139))/(JVS(504))
  X(92) = (X(92)-JVS(499)*X(133)-JVS(500)*X(134)-JVS(501)*X(136))/(JVS(498))
  X(91) = (X(91)-JVS(492)*X(93)-JVS(493)*X(94)-JVS(494)*X(134)-JVS(495)*X(135)-JVS(496)*X(138)-JVS(497)*X(139))&
            &/(JVS(491))
  X(90) = (X(90)-JVS(485)*X(113)-JVS(486)*X(134)-JVS(487)*X(135)-JVS(488)*X(139))/(JVS(484))
  X(89) = (X(89)-JVS(479)*X(117)-JVS(480)*X(134)-JVS(481)*X(135)-JVS(482)*X(139))/(JVS(478))
  X(88) = (X(88)-JVS(472)*X(117)-JVS(473)*X(134)-JVS(474)*X(135))/(JVS(471))
  X(87) = (X(87)-JVS(469)*X(134)-JVS(470)*X(136))/(JVS(468))
  X(86) = (X(86)-JVS(466)*X(134)-JVS(467)*X(136))/(JVS(465))
  X(85) = (X(85)-JVS(460)*X(134)-JVS(461)*X(135)-JVS(462)*X(136)-JVS(463)*X(138)-JVS(464)*X(139))/(JVS(459))
  X(84) = (X(84)-JVS(452)*X(113)-JVS(453)*X(130)-JVS(454)*X(134)-JVS(455)*X(136))/(JVS(451))
  X(83) = (X(83)-JVS(445)*X(118)-JVS(446)*X(134)-JVS(447)*X(135))/(JVS(444))
  X(82) = (X(82)-JVS(433)*X(86)-JVS(434)*X(87)-JVS(435)*X(92)-JVS(436)*X(104)-JVS(437)*X(107)-JVS(438)*X(109)-JVS(439)&
            &*X(130)-JVS(440)*X(134)-JVS(441)*X(135)-JVS(442)*X(136)-JVS(443)*X(139))/(JVS(432))
  X(81) = (X(81)-JVS(426)*X(132)-JVS(427)*X(134)-JVS(428)*X(138))/(JVS(425))
  X(80) = (X(80)-JVS(422)*X(131)-JVS(423)*X(134)-JVS(424)*X(135))/(JVS(421))
  X(79) = (X(79)-JVS(418)*X(116)-JVS(419)*X(134)-JVS(420)*X(135))/(JVS(417))
  X(78) = (X(78)-JVS(413)*X(134)-JVS(414)*X(135)-JVS(415)*X(138)-JVS(416)*X(139))/(JVS(412))
  X(77) = (X(77)-JVS(405)*X(84)-JVS(406)*X(113)-JVS(407)*X(134)-JVS(408)*X(136))/(JVS(404))
  X(76) = (X(76)-JVS(400)*X(106)-JVS(401)*X(134)-JVS(402)*X(135)-JVS(403)*X(139))/(JVS(399))
  X(75) = (X(75)-JVS(393)*X(114)-JVS(394)*X(126)-JVS(395)*X(130)-JVS(396)*X(133)-JVS(397)*X(134)-JVS(398)*X(138))&
            &/(JVS(392))
  X(74) = (X(74)-JVS(387)*X(134)-JVS(388)*X(135)-JVS(389)*X(139))/(JVS(386))
  X(73) = (X(73)-JVS(381)*X(134)-JVS(382)*X(137)-JVS(383)*X(138))/(JVS(380))
  X(72) = (X(72)-JVS(377)*X(90)-JVS(378)*X(134)-JVS(379)*X(135))/(JVS(376))
  X(71) = (X(71)-JVS(370)*X(105)-JVS(371)*X(113)-JVS(372)*X(117)-JVS(373)*X(133)-JVS(374)*X(134)-JVS(375)*X(139))&
            &/(JVS(369))
  X(70) = (X(70)-JVS(366)*X(134)-JVS(367)*X(135)-JVS(368)*X(139))/(JVS(365))
  X(69) = (X(69)-JVS(360)*X(119)-JVS(361)*X(134)-JVS(362)*X(135))/(JVS(359))
  X(68) = (X(68)-JVS(355)*X(132)-JVS(356)*X(134)-JVS(357)*X(135)-JVS(358)*X(137))/(JVS(354))
  X(67) = (X(67)-JVS(346)*X(108)-JVS(347)*X(113)-JVS(348)*X(132)-JVS(349)*X(134)-JVS(350)*X(135)-JVS(351)*X(136)&
            &-JVS(352)*X(137)-JVS(353)*X(140))/(JVS(345))
  X(66) = (X(66)-JVS(341)*X(88)-JVS(342)*X(117)-JVS(343)*X(134)-JVS(344)*X(139))/(JVS(340))
  X(65) = (X(65)-JVS(337)*X(130)-JVS(338)*X(135)-JVS(339)*X(139))/(JVS(336))
  X(64) = (X(64)-JVS(333)*X(134)-JVS(334)*X(135)-JVS(335)*X(139))/(JVS(332))
  X(63) = (X(63)-JVS(324)*X(85)-JVS(325)*X(134)-JVS(326)*X(136)-JVS(327)*X(138)-JVS(328)*X(139))/(JVS(323))
  X(62) = (X(62)-JVS(317)*X(134)-JVS(318)*X(135)-JVS(319)*X(139))/(JVS(316))
  X(61) = (X(61)-JVS(310)*X(93)-JVS(311)*X(134)-JVS(312)*X(135))/(JVS(309))
  X(60) = (X(60)-JVS(306)*X(134)-JVS(307)*X(136)-JVS(308)*X(139))/(JVS(305))
  X(59) = (X(59)-JVS(302)*X(134)-JVS(303)*X(136))/(JVS(301))
  X(58) = (X(58)-JVS(294)*X(86)-JVS(295)*X(87)-JVS(296)*X(104)-JVS(297)*X(107)-JVS(298)*X(109)-JVS(299)*X(134)-JVS(300)&
            &*X(136))/(JVS(293))
  X(57) = (X(57)-JVS(286)*X(120)-JVS(287)*X(131)-JVS(288)*X(133)-JVS(289)*X(134)-JVS(290)*X(137)-JVS(291)*X(139)&
            &-JVS(292)*X(140))/(JVS(285))
  X(56) = (X(56)-JVS(282)*X(134)-JVS(283)*X(135)-JVS(284)*X(138))/(JVS(281))
  X(55) = (X(55)-JVS(277)*X(77)-JVS(278)*X(84)-JVS(279)*X(130)-JVS(280)*X(134))/(JVS(276))
  X(54) = (X(54)-JVS(273)*X(129)-JVS(274)*X(134)-JVS(275)*X(135))/(JVS(272))
  X(53) = (X(53)-JVS(269)*X(134)-JVS(270)*X(135)-JVS(271)*X(140))/(JVS(268))
  X(52) = (X(52)-JVS(265)*X(89)-JVS(266)*X(134)-JVS(267)*X(135))/(JVS(264))
  X(51) = (X(51)-JVS(262)*X(134)-JVS(263)*X(139))/(JVS(261))
  X(50) = (X(50)-JVS(257)*X(121)-JVS(258)*X(134)-JVS(259)*X(135))/(JVS(256))
  X(49) = (X(49)-JVS(253)*X(97)-JVS(254)*X(134)-JVS(255)*X(135))/(JVS(252))
  X(48) = (X(48)-JVS(249)*X(102)-JVS(250)*X(134)-JVS(251)*X(135))/(JVS(248))
  X(47) = (X(47)-JVS(243)*X(70)-JVS(244)*X(93)-JVS(245)*X(94)-JVS(246)*X(134)-JVS(247)*X(139))/(JVS(242))
  X(46) = (X(46)-JVS(238)*X(93)-JVS(239)*X(94)-JVS(240)*X(139))/(JVS(237))
  X(45) = (X(45)-JVS(234)*X(102)-JVS(235)*X(134)-JVS(236)*X(140))/(JVS(233))
  X(44) = (X(44)-JVS(231)*X(94)-JVS(232)*X(139))/(JVS(230))
  X(43) = (X(43)-JVS(228)*X(133)-JVS(229)*X(134))/(JVS(227))
  X(42) = (X(42)-JVS(226)*X(134))/(JVS(225))
  X(41) = (X(41)-JVS(223)*X(134))/(JVS(222))
  X(40) = (X(40)-JVS(220)*X(133)-JVS(221)*X(138))/(JVS(219))
  X(39) = (X(39)-JVS(216)*X(94)-JVS(217)*X(134)-JVS(218)*X(135))/(JVS(215))
  X(38) = (X(38)-JVS(212)*X(74)-JVS(213)*X(134)-JVS(214)*X(135))/(JVS(211))
  X(37) = (X(37)-JVS(208)*X(64)-JVS(209)*X(134)-JVS(210)*X(139))/(JVS(207))
  X(36) = (X(36)-JVS(204)*X(64)-JVS(205)*X(134)-JVS(206)*X(135))/(JVS(203))
  X(35) = (X(35)-JVS(201)*X(41)-JVS(202)*X(134))/(JVS(200))
  X(34) = (X(34)-JVS(196)*X(93)-JVS(197)*X(94)-JVS(198)*X(139))/(JVS(195))
  X(33) = (X(33)-JVS(192)*X(70)-JVS(193)*X(134)-JVS(194)*X(135))/(JVS(191))
  X(32) = (X(32)-JVS(188)*X(85)-JVS(189)*X(134)-JVS(190)*X(135))/(JVS(187))
  X(31) = (X(31)-JVS(184)*X(62)-JVS(185)*X(134)-JVS(186)*X(135))/(JVS(183))
  X(30) = (X(30)-JVS(181)*X(60)-JVS(182)*X(139))/(JVS(180))
  X(29) = (X(29)-JVS(177)*X(127)-JVS(178)*X(134)-JVS(179)*X(135))/(JVS(176))
  X(28) = (X(28)-JVS(175)*X(134))/(JVS(174))
  X(27) = (X(27)-JVS(172)*X(134))/(JVS(171))
  X(26) = (X(26)-JVS(170)*X(134))/(JVS(169))
  X(25) = (X(25)-JVS(167)*X(134))/(JVS(166))
  X(24) = (X(24)-JVS(164)*X(134))/(JVS(163))
  X(23) = (X(23)-JVS(162)*X(134))/(JVS(161))
  X(22) = (X(22)-JVS(159)*X(134)-JVS(160)*X(135))/(JVS(158))
  X(21) = (X(21)-JVS(156)*X(78)-JVS(157)*X(138))/(JVS(155))
  X(20) = (X(20)-JVS(152)*X(43)-JVS(153)*X(133)-JVS(154)*X(134))/(JVS(151))
  X(19) = (X(19)-JVS(150)*X(84))/(JVS(149))
  X(18) = (X(18)-JVS(148)*X(134))/(JVS(147))
  X(17) = (X(17)-JVS(146)*X(134))/(JVS(145))
  X(16) = (X(16)-JVS(144)*X(134))/(JVS(143))
  X(15) = (X(15)-JVS(142)*X(134))/(JVS(141))
  X(14) = (X(14)-JVS(135)*X(86)-JVS(136)*X(87)-JVS(137)*X(104)-JVS(138)*X(108)-JVS(139)*X(134)-JVS(140)*X(136))&
            &/(JVS(134))
  X(13) = (X(13)-JVS(127)*X(14)-JVS(128)*X(86)-JVS(129)*X(87)-JVS(130)*X(104)-JVS(131)*X(108)-JVS(132)*X(134)-JVS(133)&
            &*X(136))/(JVS(126))
  X(12) = (X(12)-JVS(119)*X(13)-JVS(120)*X(86)-JVS(121)*X(87)-JVS(122)*X(104)-JVS(123)*X(108)-JVS(124)*X(134)-JVS(125)&
            &*X(136))/(JVS(118))
  X(11) = (X(11)-JVS(111)*X(12)-JVS(112)*X(86)-JVS(113)*X(87)-JVS(114)*X(104)-JVS(115)*X(108)-JVS(116)*X(134)-JVS(117)&
            &*X(136))/(JVS(110))
  X(10) = (X(10)-JVS(108)*X(11)-JVS(109)*X(134))/(JVS(107))
  X(9) = (X(9)-JVS(101)*X(18)-JVS(102)*X(24)-JVS(103)*X(27)-JVS(104)*X(41)-JVS(105)*X(58)-JVS(106)*X(134))/(JVS(100))
  X(8) = (X(8)-JVS(93)*X(9)-JVS(94)*X(18)-JVS(95)*X(24)-JVS(96)*X(27)-JVS(97)*X(41)-JVS(98)*X(58)-JVS(99)*X(134))&
           &/(JVS(92))
  X(7) = (X(7)-JVS(85)*X(8)-JVS(86)*X(18)-JVS(87)*X(24)-JVS(88)*X(27)-JVS(89)*X(41)-JVS(90)*X(58)-JVS(91)*X(134))&
           &/(JVS(84))
  X(6) = (X(6)-JVS(77)*X(7)-JVS(78)*X(18)-JVS(79)*X(24)-JVS(80)*X(27)-JVS(81)*X(41)-JVS(82)*X(58)-JVS(83)*X(134))&
           &/(JVS(76))
  X(5) = (X(5)-JVS(74)*X(6)-JVS(75)*X(134))/(JVS(73))
  X(4) = (X(4)-JVS(41)*X(62)-JVS(42)*X(64)-JVS(43)*X(65)-JVS(44)*X(70)-JVS(45)*X(74)-JVS(46)*X(76)-JVS(47)*X(78)-JVS(48)&
           &*X(85)-JVS(49)*X(89)-JVS(50)*X(90)-JVS(51)*X(91)-JVS(52)*X(93)-JVS(53)*X(94)-JVS(54)*X(95)-JVS(55)*X(97)-JVS(56)&
           &*X(98)-JVS(57)*X(102)-JVS(58)*X(103)-JVS(59)*X(105)-JVS(60)*X(116)-JVS(61)*X(117)-JVS(62)*X(118)-JVS(63)*X(119)&
           &-JVS(64)*X(121)-JVS(65)*X(122)-JVS(66)*X(127)-JVS(67)*X(129)-JVS(68)*X(131)-JVS(69)*X(132)-JVS(70)*X(135)&
           &-JVS(71)*X(137)-JVS(72)*X(140))/(JVS(40))
  X(3) = (X(3)-JVS(6)*X(51)-JVS(7)*X(60)-JVS(8)*X(62)-JVS(9)*X(64)-JVS(10)*X(65)-JVS(11)*X(70)-JVS(12)*X(74)-JVS(13)&
           &*X(76)-JVS(14)*X(78)-JVS(15)*X(85)-JVS(16)*X(89)-JVS(17)*X(90)-JVS(18)*X(91)-JVS(19)*X(93)-JVS(20)*X(94)-JVS(21)&
           &*X(95)-JVS(22)*X(97)-JVS(23)*X(98)-JVS(24)*X(102)-JVS(25)*X(103)-JVS(26)*X(105)-JVS(27)*X(116)-JVS(28)*X(117)&
           &-JVS(29)*X(118)-JVS(30)*X(119)-JVS(31)*X(121)-JVS(32)*X(122)-JVS(33)*X(127)-JVS(34)*X(129)-JVS(35)*X(131)&
           &-JVS(36)*X(132)-JVS(37)*X(137)-JVS(38)*X(139)-JVS(39)*X(140))/(JVS(5))
  X(2) = (X(2)-JVS(3)*X(20)-JVS(4)*X(134))/(JVS(2))
  X(1) = X(1)/JVS(1)
      
END SUBROUTINE mozart_mosaic_4bin_aq_KppSolve
























      SUBROUTINE mozart_mosaic_4bin_aq_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE mozart_mosaic_4bin_aq_WCOPY



      SUBROUTINE mozart_mosaic_4bin_aq_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE mozart_mosaic_4bin_aq_WAXPY




      SUBROUTINE mozart_mosaic_4bin_aq_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE mozart_mosaic_4bin_aq_WSCAL


      REAL(kind=dp) FUNCTION mozart_mosaic_4bin_aq_WLAMCH( C )








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
          CALL mozart_mosaic_4bin_aq_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      mozart_mosaic_4bin_aq_WLAMCH = Eps

      END FUNCTION mozart_mosaic_4bin_aq_WLAMCH
     
      SUBROUTINE mozart_mosaic_4bin_aq_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE mozart_mosaic_4bin_aq_WLAMCH_ADD




      SUBROUTINE mozart_mosaic_4bin_aq_SET2ZERO(N,Y)




      
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

      END SUBROUTINE mozart_mosaic_4bin_aq_SET2ZERO



      REAL(kind=dp) FUNCTION mozart_mosaic_4bin_aq_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      mozart_mosaic_4bin_aq_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        mozart_mosaic_4bin_aq_WDOT = mozart_mosaic_4bin_aq_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         mozart_mosaic_4bin_aq_WDOT = mozart_mosaic_4bin_aq_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          mozart_mosaic_4bin_aq_WDOT = mozart_mosaic_4bin_aq_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        mozart_mosaic_4bin_aq_WDOT = mozart_mosaic_4bin_aq_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION mozart_mosaic_4bin_aq_WDOT                                          




   SUBROUTINE decomp_mozart_mosaic_4bin_aq( JVS, IER )
   
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
   W( 20 ) = JVS( 3 )
   W( 134 ) = JVS( 4 )
  JVS( 2) = W( 2 )
  JVS( 3) = W( 20 )
  JVS( 4) = W( 134 )
  IF ( ABS(  JVS( 5 )) < TINY(a) ) THEN
         IER = 3                                       
         RETURN
  END IF
   W( 3 ) = JVS( 5 )
   W( 51 ) = JVS( 6 )
   W( 60 ) = JVS( 7 )
   W( 62 ) = JVS( 8 )
   W( 64 ) = JVS( 9 )
   W( 65 ) = JVS( 10 )
   W( 70 ) = JVS( 11 )
   W( 74 ) = JVS( 12 )
   W( 76 ) = JVS( 13 )
   W( 78 ) = JVS( 14 )
   W( 85 ) = JVS( 15 )
   W( 89 ) = JVS( 16 )
   W( 90 ) = JVS( 17 )
   W( 91 ) = JVS( 18 )
   W( 93 ) = JVS( 19 )
   W( 94 ) = JVS( 20 )
   W( 95 ) = JVS( 21 )
   W( 97 ) = JVS( 22 )
   W( 98 ) = JVS( 23 )
   W( 102 ) = JVS( 24 )
   W( 103 ) = JVS( 25 )
   W( 105 ) = JVS( 26 )
   W( 116 ) = JVS( 27 )
   W( 117 ) = JVS( 28 )
   W( 118 ) = JVS( 29 )
   W( 119 ) = JVS( 30 )
   W( 121 ) = JVS( 31 )
   W( 122 ) = JVS( 32 )
   W( 127 ) = JVS( 33 )
   W( 129 ) = JVS( 34 )
   W( 131 ) = JVS( 35 )
   W( 132 ) = JVS( 36 )
   W( 137 ) = JVS( 37 )
   W( 139 ) = JVS( 38 )
   W( 140 ) = JVS( 39 )
  JVS( 5) = W( 3 )
  JVS( 6) = W( 51 )
  JVS( 7) = W( 60 )
  JVS( 8) = W( 62 )
  JVS( 9) = W( 64 )
  JVS( 10) = W( 65 )
  JVS( 11) = W( 70 )
  JVS( 12) = W( 74 )
  JVS( 13) = W( 76 )
  JVS( 14) = W( 78 )
  JVS( 15) = W( 85 )
  JVS( 16) = W( 89 )
  JVS( 17) = W( 90 )
  JVS( 18) = W( 91 )
  JVS( 19) = W( 93 )
  JVS( 20) = W( 94 )
  JVS( 21) = W( 95 )
  JVS( 22) = W( 97 )
  JVS( 23) = W( 98 )
  JVS( 24) = W( 102 )
  JVS( 25) = W( 103 )
  JVS( 26) = W( 105 )
  JVS( 27) = W( 116 )
  JVS( 28) = W( 117 )
  JVS( 29) = W( 118 )
  JVS( 30) = W( 119 )
  JVS( 31) = W( 121 )
  JVS( 32) = W( 122 )
  JVS( 33) = W( 127 )
  JVS( 34) = W( 129 )
  JVS( 35) = W( 131 )
  JVS( 36) = W( 132 )
  JVS( 37) = W( 137 )
  JVS( 38) = W( 139 )
  JVS( 39) = W( 140 )
  IF ( ABS(  JVS( 40 )) < TINY(a) ) THEN
         IER = 4                                       
         RETURN
  END IF
   W( 4 ) = JVS( 40 )
   W( 62 ) = JVS( 41 )
   W( 64 ) = JVS( 42 )
   W( 65 ) = JVS( 43 )
   W( 70 ) = JVS( 44 )
   W( 74 ) = JVS( 45 )
   W( 76 ) = JVS( 46 )
   W( 78 ) = JVS( 47 )
   W( 85 ) = JVS( 48 )
   W( 89 ) = JVS( 49 )
   W( 90 ) = JVS( 50 )
   W( 91 ) = JVS( 51 )
   W( 93 ) = JVS( 52 )
   W( 94 ) = JVS( 53 )
   W( 95 ) = JVS( 54 )
   W( 97 ) = JVS( 55 )
   W( 98 ) = JVS( 56 )
   W( 102 ) = JVS( 57 )
   W( 103 ) = JVS( 58 )
   W( 105 ) = JVS( 59 )
   W( 116 ) = JVS( 60 )
   W( 117 ) = JVS( 61 )
   W( 118 ) = JVS( 62 )
   W( 119 ) = JVS( 63 )
   W( 121 ) = JVS( 64 )
   W( 122 ) = JVS( 65 )
   W( 127 ) = JVS( 66 )
   W( 129 ) = JVS( 67 )
   W( 131 ) = JVS( 68 )
   W( 132 ) = JVS( 69 )
   W( 135 ) = JVS( 70 )
   W( 137 ) = JVS( 71 )
   W( 140 ) = JVS( 72 )
  JVS( 40) = W( 4 )
  JVS( 41) = W( 62 )
  JVS( 42) = W( 64 )
  JVS( 43) = W( 65 )
  JVS( 44) = W( 70 )
  JVS( 45) = W( 74 )
  JVS( 46) = W( 76 )
  JVS( 47) = W( 78 )
  JVS( 48) = W( 85 )
  JVS( 49) = W( 89 )
  JVS( 50) = W( 90 )
  JVS( 51) = W( 91 )
  JVS( 52) = W( 93 )
  JVS( 53) = W( 94 )
  JVS( 54) = W( 95 )
  JVS( 55) = W( 97 )
  JVS( 56) = W( 98 )
  JVS( 57) = W( 102 )
  JVS( 58) = W( 103 )
  JVS( 59) = W( 105 )
  JVS( 60) = W( 116 )
  JVS( 61) = W( 117 )
  JVS( 62) = W( 118 )
  JVS( 63) = W( 119 )
  JVS( 64) = W( 121 )
  JVS( 65) = W( 122 )
  JVS( 66) = W( 127 )
  JVS( 67) = W( 129 )
  JVS( 68) = W( 131 )
  JVS( 69) = W( 132 )
  JVS( 70) = W( 135 )
  JVS( 71) = W( 137 )
  JVS( 72) = W( 140 )
  IF ( ABS(  JVS( 73 )) < TINY(a) ) THEN
         IER = 5                                       
         RETURN
  END IF
   W( 5 ) = JVS( 73 )
   W( 6 ) = JVS( 74 )
   W( 134 ) = JVS( 75 )
  JVS( 73) = W( 5 )
  JVS( 74) = W( 6 )
  JVS( 75) = W( 134 )
  IF ( ABS(  JVS( 76 )) < TINY(a) ) THEN
         IER = 6                                       
         RETURN
  END IF
   W( 6 ) = JVS( 76 )
   W( 7 ) = JVS( 77 )
   W( 18 ) = JVS( 78 )
   W( 24 ) = JVS( 79 )
   W( 27 ) = JVS( 80 )
   W( 41 ) = JVS( 81 )
   W( 58 ) = JVS( 82 )
   W( 134 ) = JVS( 83 )
  JVS( 76) = W( 6 )
  JVS( 77) = W( 7 )
  JVS( 78) = W( 18 )
  JVS( 79) = W( 24 )
  JVS( 80) = W( 27 )
  JVS( 81) = W( 41 )
  JVS( 82) = W( 58 )
  JVS( 83) = W( 134 )
  IF ( ABS(  JVS( 84 )) < TINY(a) ) THEN
         IER = 7                                       
         RETURN
  END IF
   W( 7 ) = JVS( 84 )
   W( 8 ) = JVS( 85 )
   W( 18 ) = JVS( 86 )
   W( 24 ) = JVS( 87 )
   W( 27 ) = JVS( 88 )
   W( 41 ) = JVS( 89 )
   W( 58 ) = JVS( 90 )
   W( 134 ) = JVS( 91 )
  JVS( 84) = W( 7 )
  JVS( 85) = W( 8 )
  JVS( 86) = W( 18 )
  JVS( 87) = W( 24 )
  JVS( 88) = W( 27 )
  JVS( 89) = W( 41 )
  JVS( 90) = W( 58 )
  JVS( 91) = W( 134 )
  IF ( ABS(  JVS( 92 )) < TINY(a) ) THEN
         IER = 8                                       
         RETURN
  END IF
   W( 8 ) = JVS( 92 )
   W( 9 ) = JVS( 93 )
   W( 18 ) = JVS( 94 )
   W( 24 ) = JVS( 95 )
   W( 27 ) = JVS( 96 )
   W( 41 ) = JVS( 97 )
   W( 58 ) = JVS( 98 )
   W( 134 ) = JVS( 99 )
  JVS( 92) = W( 8 )
  JVS( 93) = W( 9 )
  JVS( 94) = W( 18 )
  JVS( 95) = W( 24 )
  JVS( 96) = W( 27 )
  JVS( 97) = W( 41 )
  JVS( 98) = W( 58 )
  JVS( 99) = W( 134 )
  IF ( ABS(  JVS( 100 )) < TINY(a) ) THEN
         IER = 9                                       
         RETURN
  END IF
   W( 9 ) = JVS( 100 )
   W( 18 ) = JVS( 101 )
   W( 24 ) = JVS( 102 )
   W( 27 ) = JVS( 103 )
   W( 41 ) = JVS( 104 )
   W( 58 ) = JVS( 105 )
   W( 134 ) = JVS( 106 )
  JVS( 100) = W( 9 )
  JVS( 101) = W( 18 )
  JVS( 102) = W( 24 )
  JVS( 103) = W( 27 )
  JVS( 104) = W( 41 )
  JVS( 105) = W( 58 )
  JVS( 106) = W( 134 )
  IF ( ABS(  JVS( 107 )) < TINY(a) ) THEN
         IER = 10                                      
         RETURN
  END IF
   W( 10 ) = JVS( 107 )
   W( 11 ) = JVS( 108 )
   W( 134 ) = JVS( 109 )
  JVS( 107) = W( 10 )
  JVS( 108) = W( 11 )
  JVS( 109) = W( 134 )
  IF ( ABS(  JVS( 110 )) < TINY(a) ) THEN
         IER = 11                                      
         RETURN
  END IF
   W( 11 ) = JVS( 110 )
   W( 12 ) = JVS( 111 )
   W( 86 ) = JVS( 112 )
   W( 87 ) = JVS( 113 )
   W( 104 ) = JVS( 114 )
   W( 108 ) = JVS( 115 )
   W( 134 ) = JVS( 116 )
   W( 136 ) = JVS( 117 )
  JVS( 110) = W( 11 )
  JVS( 111) = W( 12 )
  JVS( 112) = W( 86 )
  JVS( 113) = W( 87 )
  JVS( 114) = W( 104 )
  JVS( 115) = W( 108 )
  JVS( 116) = W( 134 )
  JVS( 117) = W( 136 )
  IF ( ABS(  JVS( 118 )) < TINY(a) ) THEN
         IER = 12                                      
         RETURN
  END IF
   W( 12 ) = JVS( 118 )
   W( 13 ) = JVS( 119 )
   W( 86 ) = JVS( 120 )
   W( 87 ) = JVS( 121 )
   W( 104 ) = JVS( 122 )
   W( 108 ) = JVS( 123 )
   W( 134 ) = JVS( 124 )
   W( 136 ) = JVS( 125 )
  JVS( 118) = W( 12 )
  JVS( 119) = W( 13 )
  JVS( 120) = W( 86 )
  JVS( 121) = W( 87 )
  JVS( 122) = W( 104 )
  JVS( 123) = W( 108 )
  JVS( 124) = W( 134 )
  JVS( 125) = W( 136 )
  IF ( ABS(  JVS( 126 )) < TINY(a) ) THEN
         IER = 13                                      
         RETURN
  END IF
   W( 13 ) = JVS( 126 )
   W( 14 ) = JVS( 127 )
   W( 86 ) = JVS( 128 )
   W( 87 ) = JVS( 129 )
   W( 104 ) = JVS( 130 )
   W( 108 ) = JVS( 131 )
   W( 134 ) = JVS( 132 )
   W( 136 ) = JVS( 133 )
  JVS( 126) = W( 13 )
  JVS( 127) = W( 14 )
  JVS( 128) = W( 86 )
  JVS( 129) = W( 87 )
  JVS( 130) = W( 104 )
  JVS( 131) = W( 108 )
  JVS( 132) = W( 134 )
  JVS( 133) = W( 136 )
  IF ( ABS(  JVS( 134 )) < TINY(a) ) THEN
         IER = 14                                      
         RETURN
  END IF
   W( 14 ) = JVS( 134 )
   W( 86 ) = JVS( 135 )
   W( 87 ) = JVS( 136 )
   W( 104 ) = JVS( 137 )
   W( 108 ) = JVS( 138 )
   W( 134 ) = JVS( 139 )
   W( 136 ) = JVS( 140 )
  JVS( 134) = W( 14 )
  JVS( 135) = W( 86 )
  JVS( 136) = W( 87 )
  JVS( 137) = W( 104 )
  JVS( 138) = W( 108 )
  JVS( 139) = W( 134 )
  JVS( 140) = W( 136 )
  IF ( ABS(  JVS( 141 )) < TINY(a) ) THEN
         IER = 15                                      
         RETURN
  END IF
   W( 15 ) = JVS( 141 )
   W( 134 ) = JVS( 142 )
  JVS( 141) = W( 15 )
  JVS( 142) = W( 134 )
  IF ( ABS(  JVS( 143 )) < TINY(a) ) THEN
         IER = 16                                      
         RETURN
  END IF
   W( 16 ) = JVS( 143 )
   W( 134 ) = JVS( 144 )
  JVS( 143) = W( 16 )
  JVS( 144) = W( 134 )
  IF ( ABS(  JVS( 145 )) < TINY(a) ) THEN
         IER = 17                                      
         RETURN
  END IF
   W( 17 ) = JVS( 145 )
   W( 134 ) = JVS( 146 )
  JVS( 145) = W( 17 )
  JVS( 146) = W( 134 )
  IF ( ABS(  JVS( 147 )) < TINY(a) ) THEN
         IER = 18                                      
         RETURN
  END IF
   W( 18 ) = JVS( 147 )
   W( 134 ) = JVS( 148 )
  JVS( 147) = W( 18 )
  JVS( 148) = W( 134 )
  IF ( ABS(  JVS( 149 )) < TINY(a) ) THEN
         IER = 19                                      
         RETURN
  END IF
   W( 19 ) = JVS( 149 )
   W( 84 ) = JVS( 150 )
  JVS( 149) = W( 19 )
  JVS( 150) = W( 84 )
  IF ( ABS(  JVS( 151 )) < TINY(a) ) THEN
         IER = 20                                      
         RETURN
  END IF
   W( 20 ) = JVS( 151 )
   W( 43 ) = JVS( 152 )
   W( 133 ) = JVS( 153 )
   W( 134 ) = JVS( 154 )
  JVS( 151) = W( 20 )
  JVS( 152) = W( 43 )
  JVS( 153) = W( 133 )
  JVS( 154) = W( 134 )
  IF ( ABS(  JVS( 155 )) < TINY(a) ) THEN
         IER = 21                                      
         RETURN
  END IF
   W( 21 ) = JVS( 155 )
   W( 78 ) = JVS( 156 )
   W( 138 ) = JVS( 157 )
  JVS( 155) = W( 21 )
  JVS( 156) = W( 78 )
  JVS( 157) = W( 138 )
  IF ( ABS(  JVS( 158 )) < TINY(a) ) THEN
         IER = 22                                      
         RETURN
  END IF
   W( 22 ) = JVS( 158 )
   W( 134 ) = JVS( 159 )
   W( 135 ) = JVS( 160 )
  JVS( 158) = W( 22 )
  JVS( 159) = W( 134 )
  JVS( 160) = W( 135 )
  IF ( ABS(  JVS( 161 )) < TINY(a) ) THEN
         IER = 23                                      
         RETURN
  END IF
   W( 23 ) = JVS( 161 )
   W( 134 ) = JVS( 162 )
  JVS( 161) = W( 23 )
  JVS( 162) = W( 134 )
  IF ( ABS(  JVS( 163 )) < TINY(a) ) THEN
         IER = 24                                      
         RETURN
  END IF
   W( 24 ) = JVS( 163 )
   W( 134 ) = JVS( 164 )
  JVS( 163) = W( 24 )
  JVS( 164) = W( 134 )
  IF ( ABS(  JVS( 166 )) < TINY(a) ) THEN
         IER = 25                                      
         RETURN
  END IF
   W( 24 ) = JVS( 165 )
   W( 25 ) = JVS( 166 )
   W( 134 ) = JVS( 167 )
  a = -W( 24 ) / JVS(          163  )
  W( 24 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 164 )
  JVS( 165) = W( 24 )
  JVS( 166) = W( 25 )
  JVS( 167) = W( 134 )
  IF ( ABS(  JVS( 169 )) < TINY(a) ) THEN
         IER = 26                                      
         RETURN
  END IF
   W( 24 ) = JVS( 168 )
   W( 26 ) = JVS( 169 )
   W( 134 ) = JVS( 170 )
  a = -W( 24 ) / JVS(          163  )
  W( 24 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 164 )
  JVS( 168) = W( 24 )
  JVS( 169) = W( 26 )
  JVS( 170) = W( 134 )
  IF ( ABS(  JVS( 171 )) < TINY(a) ) THEN
         IER = 27                                      
         RETURN
  END IF
   W( 27 ) = JVS( 171 )
   W( 134 ) = JVS( 172 )
  JVS( 171) = W( 27 )
  JVS( 172) = W( 134 )
  IF ( ABS(  JVS( 174 )) < TINY(a) ) THEN
         IER = 28                                      
         RETURN
  END IF
   W( 27 ) = JVS( 173 )
   W( 28 ) = JVS( 174 )
   W( 134 ) = JVS( 175 )
  a = -W( 27 ) / JVS(          171  )
  W( 27 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 172 )
  JVS( 173) = W( 27 )
  JVS( 174) = W( 28 )
  JVS( 175) = W( 134 )
  IF ( ABS(  JVS( 176 )) < TINY(a) ) THEN
         IER = 29                                      
         RETURN
  END IF
   W( 29 ) = JVS( 176 )
   W( 127 ) = JVS( 177 )
   W( 134 ) = JVS( 178 )
   W( 135 ) = JVS( 179 )
  JVS( 176) = W( 29 )
  JVS( 177) = W( 127 )
  JVS( 178) = W( 134 )
  JVS( 179) = W( 135 )
  IF ( ABS(  JVS( 180 )) < TINY(a) ) THEN
         IER = 30                                      
         RETURN
  END IF
   W( 30 ) = JVS( 180 )
   W( 60 ) = JVS( 181 )
   W( 139 ) = JVS( 182 )
  JVS( 180) = W( 30 )
  JVS( 181) = W( 60 )
  JVS( 182) = W( 139 )
  IF ( ABS(  JVS( 183 )) < TINY(a) ) THEN
         IER = 31                                      
         RETURN
  END IF
   W( 31 ) = JVS( 183 )
   W( 62 ) = JVS( 184 )
   W( 134 ) = JVS( 185 )
   W( 135 ) = JVS( 186 )
  JVS( 183) = W( 31 )
  JVS( 184) = W( 62 )
  JVS( 185) = W( 134 )
  JVS( 186) = W( 135 )
  IF ( ABS(  JVS( 187 )) < TINY(a) ) THEN
         IER = 32                                      
         RETURN
  END IF
   W( 32 ) = JVS( 187 )
   W( 85 ) = JVS( 188 )
   W( 134 ) = JVS( 189 )
   W( 135 ) = JVS( 190 )
  JVS( 187) = W( 32 )
  JVS( 188) = W( 85 )
  JVS( 189) = W( 134 )
  JVS( 190) = W( 135 )
  IF ( ABS(  JVS( 191 )) < TINY(a) ) THEN
         IER = 33                                      
         RETURN
  END IF
   W( 33 ) = JVS( 191 )
   W( 70 ) = JVS( 192 )
   W( 134 ) = JVS( 193 )
   W( 135 ) = JVS( 194 )
  JVS( 191) = W( 33 )
  JVS( 192) = W( 70 )
  JVS( 193) = W( 134 )
  JVS( 194) = W( 135 )
  IF ( ABS(  JVS( 195 )) < TINY(a) ) THEN
         IER = 34                                      
         RETURN
  END IF
   W( 34 ) = JVS( 195 )
   W( 93 ) = JVS( 196 )
   W( 94 ) = JVS( 197 )
   W( 139 ) = JVS( 198 )
  JVS( 195) = W( 34 )
  JVS( 196) = W( 93 )
  JVS( 197) = W( 94 )
  JVS( 198) = W( 139 )
  IF ( ABS(  JVS( 200 )) < TINY(a) ) THEN
         IER = 35                                      
         RETURN
  END IF
   W( 27 ) = JVS( 199 )
   W( 35 ) = JVS( 200 )
   W( 41 ) = JVS( 201 )
   W( 134 ) = JVS( 202 )
  a = -W( 27 ) / JVS(          171  )
  W( 27 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 172 )
  JVS( 199) = W( 27 )
  JVS( 200) = W( 35 )
  JVS( 201) = W( 41 )
  JVS( 202) = W( 134 )
  IF ( ABS(  JVS( 203 )) < TINY(a) ) THEN
         IER = 36                                      
         RETURN
  END IF
   W( 36 ) = JVS( 203 )
   W( 64 ) = JVS( 204 )
   W( 134 ) = JVS( 205 )
   W( 135 ) = JVS( 206 )
  JVS( 203) = W( 36 )
  JVS( 204) = W( 64 )
  JVS( 205) = W( 134 )
  JVS( 206) = W( 135 )
  IF ( ABS(  JVS( 207 )) < TINY(a) ) THEN
         IER = 37                                      
         RETURN
  END IF
   W( 37 ) = JVS( 207 )
   W( 64 ) = JVS( 208 )
   W( 134 ) = JVS( 209 )
   W( 139 ) = JVS( 210 )
  JVS( 207) = W( 37 )
  JVS( 208) = W( 64 )
  JVS( 209) = W( 134 )
  JVS( 210) = W( 139 )
  IF ( ABS(  JVS( 211 )) < TINY(a) ) THEN
         IER = 38                                      
         RETURN
  END IF
   W( 38 ) = JVS( 211 )
   W( 74 ) = JVS( 212 )
   W( 134 ) = JVS( 213 )
   W( 135 ) = JVS( 214 )
  JVS( 211) = W( 38 )
  JVS( 212) = W( 74 )
  JVS( 213) = W( 134 )
  JVS( 214) = W( 135 )
  IF ( ABS(  JVS( 215 )) < TINY(a) ) THEN
         IER = 39                                      
         RETURN
  END IF
   W( 39 ) = JVS( 215 )
   W( 94 ) = JVS( 216 )
   W( 134 ) = JVS( 217 )
   W( 135 ) = JVS( 218 )
  JVS( 215) = W( 39 )
  JVS( 216) = W( 94 )
  JVS( 217) = W( 134 )
  JVS( 218) = W( 135 )
  IF ( ABS(  JVS( 219 )) < TINY(a) ) THEN
         IER = 40                                      
         RETURN
  END IF
   W( 40 ) = JVS( 219 )
   W( 133 ) = JVS( 220 )
   W( 138 ) = JVS( 221 )
  JVS( 219) = W( 40 )
  JVS( 220) = W( 133 )
  JVS( 221) = W( 138 )
  IF ( ABS(  JVS( 222 )) < TINY(a) ) THEN
         IER = 41                                      
         RETURN
  END IF
   W( 41 ) = JVS( 222 )
   W( 134 ) = JVS( 223 )
  JVS( 222) = W( 41 )
  JVS( 223) = W( 134 )
  IF ( ABS(  JVS( 225 )) < TINY(a) ) THEN
         IER = 42                                      
         RETURN
  END IF
   W( 41 ) = JVS( 224 )
   W( 42 ) = JVS( 225 )
   W( 134 ) = JVS( 226 )
  a = -W( 41 ) / JVS(          222  )
  W( 41 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 223 )
  JVS( 224) = W( 41 )
  JVS( 225) = W( 42 )
  JVS( 226) = W( 134 )
  IF ( ABS(  JVS( 227 )) < TINY(a) ) THEN
         IER = 43                                      
         RETURN
  END IF
   W( 43 ) = JVS( 227 )
   W( 133 ) = JVS( 228 )
   W( 134 ) = JVS( 229 )
  JVS( 227) = W( 43 )
  JVS( 228) = W( 133 )
  JVS( 229) = W( 134 )
  IF ( ABS(  JVS( 230 )) < TINY(a) ) THEN
         IER = 44                                      
         RETURN
  END IF
   W( 44 ) = JVS( 230 )
   W( 94 ) = JVS( 231 )
   W( 139 ) = JVS( 232 )
  JVS( 230) = W( 44 )
  JVS( 231) = W( 94 )
  JVS( 232) = W( 139 )
  IF ( ABS(  JVS( 233 )) < TINY(a) ) THEN
         IER = 45                                      
         RETURN
  END IF
   W( 45 ) = JVS( 233 )
   W( 102 ) = JVS( 234 )
   W( 134 ) = JVS( 235 )
   W( 140 ) = JVS( 236 )
  JVS( 233) = W( 45 )
  JVS( 234) = W( 102 )
  JVS( 235) = W( 134 )
  JVS( 236) = W( 140 )
  IF ( ABS(  JVS( 237 )) < TINY(a) ) THEN
         IER = 46                                      
         RETURN
  END IF
   W( 46 ) = JVS( 237 )
   W( 93 ) = JVS( 238 )
   W( 94 ) = JVS( 239 )
   W( 139 ) = JVS( 240 )
  JVS( 237) = W( 46 )
  JVS( 238) = W( 93 )
  JVS( 239) = W( 94 )
  JVS( 240) = W( 139 )
  IF ( ABS(  JVS( 242 )) < TINY(a) ) THEN
         IER = 47                                      
         RETURN
  END IF
   W( 25 ) = JVS( 241 )
   W( 47 ) = JVS( 242 )
   W( 70 ) = JVS( 243 )
   W( 93 ) = JVS( 244 )
   W( 94 ) = JVS( 245 )
   W( 134 ) = JVS( 246 )
   W( 139 ) = JVS( 247 )
  a = -W( 25 ) / JVS(          166  )
  W( 25 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 167 )
  JVS( 241) = W( 25 )
  JVS( 242) = W( 47 )
  JVS( 243) = W( 70 )
  JVS( 244) = W( 93 )
  JVS( 245) = W( 94 )
  JVS( 246) = W( 134 )
  JVS( 247) = W( 139 )
  IF ( ABS(  JVS( 248 )) < TINY(a) ) THEN
         IER = 48                                      
         RETURN
  END IF
   W( 48 ) = JVS( 248 )
   W( 102 ) = JVS( 249 )
   W( 134 ) = JVS( 250 )
   W( 135 ) = JVS( 251 )
  JVS( 248) = W( 48 )
  JVS( 249) = W( 102 )
  JVS( 250) = W( 134 )
  JVS( 251) = W( 135 )
  IF ( ABS(  JVS( 252 )) < TINY(a) ) THEN
         IER = 49                                      
         RETURN
  END IF
   W( 49 ) = JVS( 252 )
   W( 97 ) = JVS( 253 )
   W( 134 ) = JVS( 254 )
   W( 135 ) = JVS( 255 )
  JVS( 252) = W( 49 )
  JVS( 253) = W( 97 )
  JVS( 254) = W( 134 )
  JVS( 255) = W( 135 )
  IF ( ABS(  JVS( 256 )) < TINY(a) ) THEN
         IER = 50                                      
         RETURN
  END IF
   W( 50 ) = JVS( 256 )
   W( 121 ) = JVS( 257 )
   W( 134 ) = JVS( 258 )
   W( 135 ) = JVS( 259 )
  JVS( 256) = W( 50 )
  JVS( 257) = W( 121 )
  JVS( 258) = W( 134 )
  JVS( 259) = W( 135 )
  IF ( ABS(  JVS( 261 )) < TINY(a) ) THEN
         IER = 51                                      
         RETURN
  END IF
   W( 18 ) = JVS( 260 )
   W( 51 ) = JVS( 261 )
   W( 134 ) = JVS( 262 )
   W( 139 ) = JVS( 263 )
  a = -W( 18 ) / JVS(          147  )
  W( 18 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 148 )
  JVS( 260) = W( 18 )
  JVS( 261) = W( 51 )
  JVS( 262) = W( 134 )
  JVS( 263) = W( 139 )
  IF ( ABS(  JVS( 264 )) < TINY(a) ) THEN
         IER = 52                                      
         RETURN
  END IF
   W( 52 ) = JVS( 264 )
   W( 89 ) = JVS( 265 )
   W( 134 ) = JVS( 266 )
   W( 135 ) = JVS( 267 )
  JVS( 264) = W( 52 )
  JVS( 265) = W( 89 )
  JVS( 266) = W( 134 )
  JVS( 267) = W( 135 )
  IF ( ABS(  JVS( 268 )) < TINY(a) ) THEN
         IER = 53                                      
         RETURN
  END IF
   W( 53 ) = JVS( 268 )
   W( 134 ) = JVS( 269 )
   W( 135 ) = JVS( 270 )
   W( 140 ) = JVS( 271 )
  JVS( 268) = W( 53 )
  JVS( 269) = W( 134 )
  JVS( 270) = W( 135 )
  JVS( 271) = W( 140 )
  IF ( ABS(  JVS( 272 )) < TINY(a) ) THEN
         IER = 54                                      
         RETURN
  END IF
   W( 54 ) = JVS( 272 )
   W( 129 ) = JVS( 273 )
   W( 134 ) = JVS( 274 )
   W( 135 ) = JVS( 275 )
  JVS( 272) = W( 54 )
  JVS( 273) = W( 129 )
  JVS( 274) = W( 134 )
  JVS( 275) = W( 135 )
  IF ( ABS(  JVS( 276 )) < TINY(a) ) THEN
         IER = 55                                      
         RETURN
  END IF
   W( 55 ) = JVS( 276 )
   W( 77 ) = JVS( 277 )
   W( 84 ) = JVS( 278 )
   W( 130 ) = JVS( 279 )
   W( 134 ) = JVS( 280 )
  JVS( 276) = W( 55 )
  JVS( 277) = W( 77 )
  JVS( 278) = W( 84 )
  JVS( 279) = W( 130 )
  JVS( 280) = W( 134 )
  IF ( ABS(  JVS( 281 )) < TINY(a) ) THEN
         IER = 56                                      
         RETURN
  END IF
   W( 56 ) = JVS( 281 )
   W( 134 ) = JVS( 282 )
   W( 135 ) = JVS( 283 )
   W( 138 ) = JVS( 284 )
  JVS( 281) = W( 56 )
  JVS( 282) = W( 134 )
  JVS( 283) = W( 135 )
  JVS( 284) = W( 138 )
  IF ( ABS(  JVS( 285 )) < TINY(a) ) THEN
         IER = 57                                      
         RETURN
  END IF
   W( 57 ) = JVS( 285 )
   W( 120 ) = JVS( 286 )
   W( 131 ) = JVS( 287 )
   W( 133 ) = JVS( 288 )
   W( 134 ) = JVS( 289 )
   W( 137 ) = JVS( 290 )
   W( 139 ) = JVS( 291 )
   W( 140 ) = JVS( 292 )
  JVS( 285) = W( 57 )
  JVS( 286) = W( 120 )
  JVS( 287) = W( 131 )
  JVS( 288) = W( 133 )
  JVS( 289) = W( 134 )
  JVS( 290) = W( 137 )
  JVS( 291) = W( 139 )
  JVS( 292) = W( 140 )
  IF ( ABS(  JVS( 293 )) < TINY(a) ) THEN
         IER = 58                                      
         RETURN
  END IF
   W( 58 ) = JVS( 293 )
   W( 86 ) = JVS( 294 )
   W( 87 ) = JVS( 295 )
   W( 104 ) = JVS( 296 )
   W( 107 ) = JVS( 297 )
   W( 109 ) = JVS( 298 )
   W( 134 ) = JVS( 299 )
   W( 136 ) = JVS( 300 )
  JVS( 293) = W( 58 )
  JVS( 294) = W( 86 )
  JVS( 295) = W( 87 )
  JVS( 296) = W( 104 )
  JVS( 297) = W( 107 )
  JVS( 298) = W( 109 )
  JVS( 299) = W( 134 )
  JVS( 300) = W( 136 )
  IF ( ABS(  JVS( 301 )) < TINY(a) ) THEN
         IER = 59                                      
         RETURN
  END IF
   W( 59 ) = JVS( 301 )
   W( 134 ) = JVS( 302 )
   W( 136 ) = JVS( 303 )
  JVS( 301) = W( 59 )
  JVS( 302) = W( 134 )
  JVS( 303) = W( 136 )
  IF ( ABS(  JVS( 305 )) < TINY(a) ) THEN
         IER = 60                                      
         RETURN
  END IF
   W( 59 ) = JVS( 304 )
   W( 60 ) = JVS( 305 )
   W( 134 ) = JVS( 306 )
   W( 136 ) = JVS( 307 )
   W( 139 ) = JVS( 308 )
  a = -W( 59 ) / JVS(          301  )
  W( 59 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 302 )
  W( 136 ) = W( 136 ) + a*JVS( 303 )
  JVS( 304) = W( 59 )
  JVS( 305) = W( 60 )
  JVS( 306) = W( 134 )
  JVS( 307) = W( 136 )
  JVS( 308) = W( 139 )
  IF ( ABS(  JVS( 309 )) < TINY(a) ) THEN
         IER = 61                                      
         RETURN
  END IF
   W( 61 ) = JVS( 309 )
   W( 93 ) = JVS( 310 )
   W( 134 ) = JVS( 311 )
   W( 135 ) = JVS( 312 )
  JVS( 309) = W( 61 )
  JVS( 310) = W( 93 )
  JVS( 311) = W( 134 )
  JVS( 312) = W( 135 )
  IF ( ABS(  JVS( 316 )) < TINY(a) ) THEN
         IER = 62                                      
         RETURN
  END IF
   W( 26 ) = JVS( 313 )
   W( 28 ) = JVS( 314 )
   W( 31 ) = JVS( 315 )
   W( 62 ) = JVS( 316 )
   W( 134 ) = JVS( 317 )
   W( 135 ) = JVS( 318 )
   W( 139 ) = JVS( 319 )
  a = -W( 26 ) / JVS(          169  )
  W( 26 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 170 )
  a = -W( 28 ) / JVS(          174  )
  W( 28 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 175 )
  a = -W( 31 ) / JVS(          183  )
  W( 31 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 184 )
  W( 134 ) = W( 134 ) + a*JVS( 185 )
  W( 135 ) = W( 135 ) + a*JVS( 186 )
  JVS( 313) = W( 26 )
  JVS( 314) = W( 28 )
  JVS( 315) = W( 31 )
  JVS( 316) = W( 62 )
  JVS( 317) = W( 134 )
  JVS( 318) = W( 135 )
  JVS( 319) = W( 139 )
  IF ( ABS(  JVS( 323 )) < TINY(a) ) THEN
         IER = 63                                      
         RETURN
  END IF
   W( 26 ) = JVS( 320 )
   W( 28 ) = JVS( 321 )
   W( 42 ) = JVS( 322 )
   W( 63 ) = JVS( 323 )
   W( 85 ) = JVS( 324 )
   W( 134 ) = JVS( 325 )
   W( 136 ) = JVS( 326 )
   W( 138 ) = JVS( 327 )
   W( 139 ) = JVS( 328 )
  a = -W( 26 ) / JVS(          169  )
  W( 26 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 170 )
  a = -W( 28 ) / JVS(          174  )
  W( 28 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 175 )
  a = -W( 42 ) / JVS(          225  )
  W( 42 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 226 )
  JVS( 320) = W( 26 )
  JVS( 321) = W( 28 )
  JVS( 322) = W( 42 )
  JVS( 323) = W( 63 )
  JVS( 324) = W( 85 )
  JVS( 325) = W( 134 )
  JVS( 326) = W( 136 )
  JVS( 327) = W( 138 )
  JVS( 328) = W( 139 )
  IF ( ABS(  JVS( 332 )) < TINY(a) ) THEN
         IER = 64                                      
         RETURN
  END IF
   W( 27 ) = JVS( 329 )
   W( 36 ) = JVS( 330 )
   W( 41 ) = JVS( 331 )
   W( 64 ) = JVS( 332 )
   W( 134 ) = JVS( 333 )
   W( 135 ) = JVS( 334 )
   W( 139 ) = JVS( 335 )
  a = -W( 27 ) / JVS(          171  )
  W( 27 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 172 )
  a = -W( 36 ) / JVS(          203  )
  W( 36 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 204 )
  W( 134 ) = W( 134 ) + a*JVS( 205 )
  W( 135 ) = W( 135 ) + a*JVS( 206 )
  a = -W( 41 ) / JVS(          222  )
  W( 41 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 223 )
  JVS( 329) = W( 27 )
  JVS( 330) = W( 36 )
  JVS( 331) = W( 41 )
  JVS( 332) = W( 64 )
  JVS( 333) = W( 134 )
  JVS( 334) = W( 135 )
  JVS( 335) = W( 139 )
  IF ( ABS(  JVS( 336 )) < TINY(a) ) THEN
         IER = 65                                      
         RETURN
  END IF
   W( 65 ) = JVS( 336 )
   W( 130 ) = JVS( 337 )
   W( 135 ) = JVS( 338 )
   W( 139 ) = JVS( 339 )
  JVS( 336) = W( 65 )
  JVS( 337) = W( 130 )
  JVS( 338) = W( 135 )
  JVS( 339) = W( 139 )
  IF ( ABS(  JVS( 340 )) < TINY(a) ) THEN
         IER = 66                                      
         RETURN
  END IF
   W( 66 ) = JVS( 340 )
   W( 88 ) = JVS( 341 )
   W( 117 ) = JVS( 342 )
   W( 134 ) = JVS( 343 )
   W( 139 ) = JVS( 344 )
  JVS( 340) = W( 66 )
  JVS( 341) = W( 88 )
  JVS( 342) = W( 117 )
  JVS( 343) = W( 134 )
  JVS( 344) = W( 139 )
  IF ( ABS(  JVS( 345 )) < TINY(a) ) THEN
         IER = 67                                      
         RETURN
  END IF
   W( 67 ) = JVS( 345 )
   W( 108 ) = JVS( 346 )
   W( 113 ) = JVS( 347 )
   W( 132 ) = JVS( 348 )
   W( 134 ) = JVS( 349 )
   W( 135 ) = JVS( 350 )
   W( 136 ) = JVS( 351 )
   W( 137 ) = JVS( 352 )
   W( 140 ) = JVS( 353 )
  JVS( 345) = W( 67 )
  JVS( 346) = W( 108 )
  JVS( 347) = W( 113 )
  JVS( 348) = W( 132 )
  JVS( 349) = W( 134 )
  JVS( 350) = W( 135 )
  JVS( 351) = W( 136 )
  JVS( 352) = W( 137 )
  JVS( 353) = W( 140 )
  IF ( ABS(  JVS( 354 )) < TINY(a) ) THEN
         IER = 68                                      
         RETURN
  END IF
   W( 68 ) = JVS( 354 )
   W( 132 ) = JVS( 355 )
   W( 134 ) = JVS( 356 )
   W( 135 ) = JVS( 357 )
   W( 137 ) = JVS( 358 )
  JVS( 354) = W( 68 )
  JVS( 355) = W( 132 )
  JVS( 356) = W( 134 )
  JVS( 357) = W( 135 )
  JVS( 358) = W( 137 )
  IF ( ABS(  JVS( 359 )) < TINY(a) ) THEN
         IER = 69                                      
         RETURN
  END IF
   W( 69 ) = JVS( 359 )
   W( 119 ) = JVS( 360 )
   W( 134 ) = JVS( 361 )
   W( 135 ) = JVS( 362 )
  JVS( 359) = W( 69 )
  JVS( 360) = W( 119 )
  JVS( 361) = W( 134 )
  JVS( 362) = W( 135 )
  IF ( ABS(  JVS( 365 )) < TINY(a) ) THEN
         IER = 70                                      
         RETURN
  END IF
   W( 24 ) = JVS( 363 )
   W( 33 ) = JVS( 364 )
   W( 70 ) = JVS( 365 )
   W( 134 ) = JVS( 366 )
   W( 135 ) = JVS( 367 )
   W( 139 ) = JVS( 368 )
  a = -W( 24 ) / JVS(          163  )
  W( 24 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 164 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 192 )
  W( 134 ) = W( 134 ) + a*JVS( 193 )
  W( 135 ) = W( 135 ) + a*JVS( 194 )
  JVS( 363) = W( 24 )
  JVS( 364) = W( 33 )
  JVS( 365) = W( 70 )
  JVS( 366) = W( 134 )
  JVS( 367) = W( 135 )
  JVS( 368) = W( 139 )
  IF ( ABS(  JVS( 369 )) < TINY(a) ) THEN
         IER = 71                                      
         RETURN
  END IF
   W( 71 ) = JVS( 369 )
   W( 105 ) = JVS( 370 )
   W( 113 ) = JVS( 371 )
   W( 117 ) = JVS( 372 )
   W( 133 ) = JVS( 373 )
   W( 134 ) = JVS( 374 )
   W( 139 ) = JVS( 375 )
  JVS( 369) = W( 71 )
  JVS( 370) = W( 105 )
  JVS( 371) = W( 113 )
  JVS( 372) = W( 117 )
  JVS( 373) = W( 133 )
  JVS( 374) = W( 134 )
  JVS( 375) = W( 139 )
  IF ( ABS(  JVS( 376 )) < TINY(a) ) THEN
         IER = 72                                      
         RETURN
  END IF
   W( 72 ) = JVS( 376 )
   W( 90 ) = JVS( 377 )
   W( 134 ) = JVS( 378 )
   W( 135 ) = JVS( 379 )
  JVS( 376) = W( 72 )
  JVS( 377) = W( 90 )
  JVS( 378) = W( 134 )
  JVS( 379) = W( 135 )
  IF ( ABS(  JVS( 380 )) < TINY(a) ) THEN
         IER = 73                                      
         RETURN
  END IF
   W( 73 ) = JVS( 380 )
   W( 134 ) = JVS( 381 )
   W( 137 ) = JVS( 382 )
   W( 138 ) = JVS( 383 )
  JVS( 380) = W( 73 )
  JVS( 381) = W( 134 )
  JVS( 382) = W( 137 )
  JVS( 383) = W( 138 )
  IF ( ABS(  JVS( 386 )) < TINY(a) ) THEN
         IER = 74                                      
         RETURN
  END IF
   W( 38 ) = JVS( 384 )
   W( 42 ) = JVS( 385 )
   W( 74 ) = JVS( 386 )
   W( 134 ) = JVS( 387 )
   W( 135 ) = JVS( 388 )
   W( 139 ) = JVS( 389 )
  a = -W( 38 ) / JVS(          211  )
  W( 38 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 212 )
  W( 134 ) = W( 134 ) + a*JVS( 213 )
  W( 135 ) = W( 135 ) + a*JVS( 214 )
  a = -W( 42 ) / JVS(          225  )
  W( 42 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 226 )
  JVS( 384) = W( 38 )
  JVS( 385) = W( 42 )
  JVS( 386) = W( 74 )
  JVS( 387) = W( 134 )
  JVS( 388) = W( 135 )
  JVS( 389) = W( 139 )
  IF ( ABS(  JVS( 392 )) < TINY(a) ) THEN
         IER = 75                                      
         RETURN
  END IF
   W( 40 ) = JVS( 390 )
   W( 43 ) = JVS( 391 )
   W( 75 ) = JVS( 392 )
   W( 114 ) = JVS( 393 )
   W( 126 ) = JVS( 394 )
   W( 130 ) = JVS( 395 )
   W( 133 ) = JVS( 396 )
   W( 134 ) = JVS( 397 )
   W( 138 ) = JVS( 398 )
  a = -W( 40 ) / JVS(          219  )
  W( 40 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 220 )
  W( 138 ) = W( 138 ) + a*JVS( 221 )
  a = -W( 43 ) / JVS(          227  )
  W( 43 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 228 )
  W( 134 ) = W( 134 ) + a*JVS( 229 )
  JVS( 390) = W( 40 )
  JVS( 391) = W( 43 )
  JVS( 392) = W( 75 )
  JVS( 393) = W( 114 )
  JVS( 394) = W( 126 )
  JVS( 395) = W( 130 )
  JVS( 396) = W( 133 )
  JVS( 397) = W( 134 )
  JVS( 398) = W( 138 )
  IF ( ABS(  JVS( 399 )) < TINY(a) ) THEN
         IER = 76                                      
         RETURN
  END IF
   W( 76 ) = JVS( 399 )
   W( 106 ) = JVS( 400 )
   W( 134 ) = JVS( 401 )
   W( 135 ) = JVS( 402 )
   W( 139 ) = JVS( 403 )
  JVS( 399) = W( 76 )
  JVS( 400) = W( 106 )
  JVS( 401) = W( 134 )
  JVS( 402) = W( 135 )
  JVS( 403) = W( 139 )
  IF ( ABS(  JVS( 404 )) < TINY(a) ) THEN
         IER = 77                                      
         RETURN
  END IF
   W( 77 ) = JVS( 404 )
   W( 84 ) = JVS( 405 )
   W( 113 ) = JVS( 406 )
   W( 134 ) = JVS( 407 )
   W( 136 ) = JVS( 408 )
  JVS( 404) = W( 77 )
  JVS( 405) = W( 84 )
  JVS( 406) = W( 113 )
  JVS( 407) = W( 134 )
  JVS( 408) = W( 136 )
  IF ( ABS(  JVS( 412 )) < TINY(a) ) THEN
         IER = 78                                      
         RETURN
  END IF
   W( 21 ) = JVS( 409 )
   W( 37 ) = JVS( 410 )
   W( 64 ) = JVS( 411 )
   W( 78 ) = JVS( 412 )
   W( 134 ) = JVS( 413 )
   W( 135 ) = JVS( 414 )
   W( 138 ) = JVS( 415 )
   W( 139 ) = JVS( 416 )
  a = -W( 21 ) / JVS(          155  )
  W( 21 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 156 )
  W( 138 ) = W( 138 ) + a*JVS( 157 )
  a = -W( 37 ) / JVS(          207  )
  W( 37 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 208 )
  W( 134 ) = W( 134 ) + a*JVS( 209 )
  W( 139 ) = W( 139 ) + a*JVS( 210 )
  a = -W( 64 ) / JVS(          332  )
  W( 64 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 333 )
  W( 135 ) = W( 135 ) + a*JVS( 334 )
  W( 139 ) = W( 139 ) + a*JVS( 335 )
  JVS( 409) = W( 21 )
  JVS( 410) = W( 37 )
  JVS( 411) = W( 64 )
  JVS( 412) = W( 78 )
  JVS( 413) = W( 134 )
  JVS( 414) = W( 135 )
  JVS( 415) = W( 138 )
  JVS( 416) = W( 139 )
  IF ( ABS(  JVS( 417 )) < TINY(a) ) THEN
         IER = 79                                      
         RETURN
  END IF
   W( 79 ) = JVS( 417 )
   W( 116 ) = JVS( 418 )
   W( 134 ) = JVS( 419 )
   W( 135 ) = JVS( 420 )
  JVS( 417) = W( 79 )
  JVS( 418) = W( 116 )
  JVS( 419) = W( 134 )
  JVS( 420) = W( 135 )
  IF ( ABS(  JVS( 421 )) < TINY(a) ) THEN
         IER = 80                                      
         RETURN
  END IF
   W( 80 ) = JVS( 421 )
   W( 131 ) = JVS( 422 )
   W( 134 ) = JVS( 423 )
   W( 135 ) = JVS( 424 )
  JVS( 421) = W( 80 )
  JVS( 422) = W( 131 )
  JVS( 423) = W( 134 )
  JVS( 424) = W( 135 )
  IF ( ABS(  JVS( 425 )) < TINY(a) ) THEN
         IER = 81                                      
         RETURN
  END IF
   W( 81 ) = JVS( 425 )
   W( 132 ) = JVS( 426 )
   W( 134 ) = JVS( 427 )
   W( 138 ) = JVS( 428 )
  JVS( 425) = W( 81 )
  JVS( 426) = W( 132 )
  JVS( 427) = W( 134 )
  JVS( 428) = W( 138 )
  IF ( ABS(  JVS( 432 )) < TINY(a) ) THEN
         IER = 82                                      
         RETURN
  END IF
   W( 23 ) = JVS( 429 )
   W( 59 ) = JVS( 430 )
   W( 65 ) = JVS( 431 )
   W( 82 ) = JVS( 432 )
   W( 86 ) = JVS( 433 )
   W( 87 ) = JVS( 434 )
   W( 92 ) = JVS( 435 )
   W( 104 ) = JVS( 436 )
   W( 107 ) = JVS( 437 )
   W( 109 ) = JVS( 438 )
   W( 130 ) = JVS( 439 )
   W( 134 ) = JVS( 440 )
   W( 135 ) = JVS( 441 )
   W( 136 ) = JVS( 442 )
   W( 139 ) = JVS( 443 )
  a = -W( 23 ) / JVS(          161  )
  W( 23 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 162 )
  a = -W( 59 ) / JVS(          301  )
  W( 59 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 302 )
  W( 136 ) = W( 136 ) + a*JVS( 303 )
  a = -W( 65 ) / JVS(          336  )
  W( 65 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 337 )
  W( 135 ) = W( 135 ) + a*JVS( 338 )
  W( 139 ) = W( 139 ) + a*JVS( 339 )
  JVS( 429) = W( 23 )
  JVS( 430) = W( 59 )
  JVS( 431) = W( 65 )
  JVS( 432) = W( 82 )
  JVS( 433) = W( 86 )
  JVS( 434) = W( 87 )
  JVS( 435) = W( 92 )
  JVS( 436) = W( 104 )
  JVS( 437) = W( 107 )
  JVS( 438) = W( 109 )
  JVS( 439) = W( 130 )
  JVS( 440) = W( 134 )
  JVS( 441) = W( 135 )
  JVS( 442) = W( 136 )
  JVS( 443) = W( 139 )
  IF ( ABS(  JVS( 444 )) < TINY(a) ) THEN
         IER = 83                                      
         RETURN
  END IF
   W( 83 ) = JVS( 444 )
   W( 118 ) = JVS( 445 )
   W( 134 ) = JVS( 446 )
   W( 135 ) = JVS( 447 )
  JVS( 444) = W( 83 )
  JVS( 445) = W( 118 )
  JVS( 446) = W( 134 )
  JVS( 447) = W( 135 )
  IF ( ABS(  JVS( 451 )) < TINY(a) ) THEN
         IER = 84                                      
         RETURN
  END IF
   W( 19 ) = JVS( 448 )
   W( 55 ) = JVS( 449 )
   W( 77 ) = JVS( 450 )
   W( 84 ) = JVS( 451 )
   W( 113 ) = JVS( 452 )
   W( 130 ) = JVS( 453 )
   W( 134 ) = JVS( 454 )
   W( 136 ) = JVS( 455 )
  a = -W( 19 ) / JVS(          149  )
  W( 19 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 150 )
  a = -W( 55 ) / JVS(          276  )
  W( 55 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 277 )
  W( 84 ) = W( 84 ) + a*JVS( 278 )
  W( 130 ) = W( 130 ) + a*JVS( 279 )
  W( 134 ) = W( 134 ) + a*JVS( 280 )
  a = -W( 77 ) / JVS(          404  )
  W( 77 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 405 )
  W( 113 ) = W( 113 ) + a*JVS( 406 )
  W( 134 ) = W( 134 ) + a*JVS( 407 )
  W( 136 ) = W( 136 ) + a*JVS( 408 )
  JVS( 448) = W( 19 )
  JVS( 449) = W( 55 )
  JVS( 450) = W( 77 )
  JVS( 451) = W( 84 )
  JVS( 452) = W( 113 )
  JVS( 453) = W( 130 )
  JVS( 454) = W( 134 )
  JVS( 455) = W( 136 )
  IF ( ABS(  JVS( 459 )) < TINY(a) ) THEN
         IER = 85                                      
         RETURN
  END IF
   W( 32 ) = JVS( 456 )
   W( 63 ) = JVS( 457 )
   W( 78 ) = JVS( 458 )
   W( 85 ) = JVS( 459 )
   W( 134 ) = JVS( 460 )
   W( 135 ) = JVS( 461 )
   W( 136 ) = JVS( 462 )
   W( 138 ) = JVS( 463 )
   W( 139 ) = JVS( 464 )
  a = -W( 32 ) / JVS(          187  )
  W( 32 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 188 )
  W( 134 ) = W( 134 ) + a*JVS( 189 )
  W( 135 ) = W( 135 ) + a*JVS( 190 )
  a = -W( 63 ) / JVS(          323  )
  W( 63 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 324 )
  W( 134 ) = W( 134 ) + a*JVS( 325 )
  W( 136 ) = W( 136 ) + a*JVS( 326 )
  W( 138 ) = W( 138 ) + a*JVS( 327 )
  W( 139 ) = W( 139 ) + a*JVS( 328 )
  a = -W( 78 ) / JVS(          412  )
  W( 78 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 413 )
  W( 135 ) = W( 135 ) + a*JVS( 414 )
  W( 138 ) = W( 138 ) + a*JVS( 415 )
  W( 139 ) = W( 139 ) + a*JVS( 416 )
  JVS( 456) = W( 32 )
  JVS( 457) = W( 63 )
  JVS( 458) = W( 78 )
  JVS( 459) = W( 85 )
  JVS( 460) = W( 134 )
  JVS( 461) = W( 135 )
  JVS( 462) = W( 136 )
  JVS( 463) = W( 138 )
  JVS( 464) = W( 139 )
  IF ( ABS(  JVS( 465 )) < TINY(a) ) THEN
         IER = 86                                      
         RETURN
  END IF
   W( 86 ) = JVS( 465 )
   W( 134 ) = JVS( 466 )
   W( 136 ) = JVS( 467 )
  JVS( 465) = W( 86 )
  JVS( 466) = W( 134 )
  JVS( 467) = W( 136 )
  IF ( ABS(  JVS( 468 )) < TINY(a) ) THEN
         IER = 87                                      
         RETURN
  END IF
   W( 87 ) = JVS( 468 )
   W( 134 ) = JVS( 469 )
   W( 136 ) = JVS( 470 )
  JVS( 468) = W( 87 )
  JVS( 469) = W( 134 )
  JVS( 470) = W( 136 )
  IF ( ABS(  JVS( 471 )) < TINY(a) ) THEN
         IER = 88                                      
         RETURN
  END IF
   W( 88 ) = JVS( 471 )
   W( 117 ) = JVS( 472 )
   W( 134 ) = JVS( 473 )
   W( 135 ) = JVS( 474 )
  JVS( 471) = W( 88 )
  JVS( 472) = W( 117 )
  JVS( 473) = W( 134 )
  JVS( 474) = W( 135 )
  IF ( ABS(  JVS( 478 )) < TINY(a) ) THEN
         IER = 89                                      
         RETURN
  END IF
   W( 52 ) = JVS( 475 )
   W( 66 ) = JVS( 476 )
   W( 88 ) = JVS( 477 )
   W( 89 ) = JVS( 478 )
   W( 117 ) = JVS( 479 )
   W( 134 ) = JVS( 480 )
   W( 135 ) = JVS( 481 )
   W( 139 ) = JVS( 482 )
  a = -W( 52 ) / JVS(          264  )
  W( 52 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 265 )
  W( 134 ) = W( 134 ) + a*JVS( 266 )
  W( 135 ) = W( 135 ) + a*JVS( 267 )
  a = -W( 66 ) / JVS(          340  )
  W( 66 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 341 )
  W( 117 ) = W( 117 ) + a*JVS( 342 )
  W( 134 ) = W( 134 ) + a*JVS( 343 )
  W( 139 ) = W( 139 ) + a*JVS( 344 )
  a = -W( 88 ) / JVS(          471  )
  W( 88 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 472 )
  W( 134 ) = W( 134 ) + a*JVS( 473 )
  W( 135 ) = W( 135 ) + a*JVS( 474 )
  JVS( 475) = W( 52 )
  JVS( 476) = W( 66 )
  JVS( 477) = W( 88 )
  JVS( 478) = W( 89 )
  JVS( 479) = W( 117 )
  JVS( 480) = W( 134 )
  JVS( 481) = W( 135 )
  JVS( 482) = W( 139 )
  IF ( ABS(  JVS( 484 )) < TINY(a) ) THEN
         IER = 90                                      
         RETURN
  END IF
   W( 72 ) = JVS( 483 )
   W( 90 ) = JVS( 484 )
   W( 113 ) = JVS( 485 )
   W( 134 ) = JVS( 486 )
   W( 135 ) = JVS( 487 )
   W( 139 ) = JVS( 488 )
  a = -W( 72 ) / JVS(          376  )
  W( 72 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 377 )
  W( 134 ) = W( 134 ) + a*JVS( 378 )
  W( 135 ) = W( 135 ) + a*JVS( 379 )
  JVS( 483) = W( 72 )
  JVS( 484) = W( 90 )
  JVS( 485) = W( 113 )
  JVS( 486) = W( 134 )
  JVS( 487) = W( 135 )
  JVS( 488) = W( 139 )
  IF ( ABS(  JVS( 491 )) < TINY(a) ) THEN
         IER = 91                                      
         RETURN
  END IF
   W( 47 ) = JVS( 489 )
   W( 70 ) = JVS( 490 )
   W( 91 ) = JVS( 491 )
   W( 93 ) = JVS( 492 )
   W( 94 ) = JVS( 493 )
   W( 134 ) = JVS( 494 )
   W( 135 ) = JVS( 495 )
   W( 138 ) = JVS( 496 )
   W( 139 ) = JVS( 497 )
  a = -W( 47 ) / JVS(          242  )
  W( 47 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 243 )
  W( 93 ) = W( 93 ) + a*JVS( 244 )
  W( 94 ) = W( 94 ) + a*JVS( 245 )
  W( 134 ) = W( 134 ) + a*JVS( 246 )
  W( 139 ) = W( 139 ) + a*JVS( 247 )
  a = -W( 70 ) / JVS(          365  )
  W( 70 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 366 )
  W( 135 ) = W( 135 ) + a*JVS( 367 )
  W( 139 ) = W( 139 ) + a*JVS( 368 )
  JVS( 489) = W( 47 )
  JVS( 490) = W( 70 )
  JVS( 491) = W( 91 )
  JVS( 492) = W( 93 )
  JVS( 493) = W( 94 )
  JVS( 494) = W( 134 )
  JVS( 495) = W( 135 )
  JVS( 496) = W( 138 )
  JVS( 497) = W( 139 )
  IF ( ABS(  JVS( 498 )) < TINY(a) ) THEN
         IER = 92                                      
         RETURN
  END IF
   W( 92 ) = JVS( 498 )
   W( 133 ) = JVS( 499 )
   W( 134 ) = JVS( 500 )
   W( 136 ) = JVS( 501 )
  JVS( 498) = W( 92 )
  JVS( 499) = W( 133 )
  JVS( 500) = W( 134 )
  JVS( 501) = W( 136 )
  IF ( ABS(  JVS( 504 )) < TINY(a) ) THEN
         IER = 93                                      
         RETURN
  END IF
   W( 27 ) = JVS( 502 )
   W( 61 ) = JVS( 503 )
   W( 93 ) = JVS( 504 )
   W( 134 ) = JVS( 505 )
   W( 135 ) = JVS( 506 )
   W( 139 ) = JVS( 507 )
  a = -W( 27 ) / JVS(          171  )
  W( 27 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 172 )
  a = -W( 61 ) / JVS(          309  )
  W( 61 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 310 )
  W( 134 ) = W( 134 ) + a*JVS( 311 )
  W( 135 ) = W( 135 ) + a*JVS( 312 )
  JVS( 502) = W( 27 )
  JVS( 503) = W( 61 )
  JVS( 504) = W( 93 )
  JVS( 505) = W( 134 )
  JVS( 506) = W( 135 )
  JVS( 507) = W( 139 )
  IF ( ABS(  JVS( 510 )) < TINY(a) ) THEN
         IER = 94                                      
         RETURN
  END IF
   W( 39 ) = JVS( 508 )
   W( 41 ) = JVS( 509 )
   W( 94 ) = JVS( 510 )
   W( 134 ) = JVS( 511 )
   W( 135 ) = JVS( 512 )
   W( 139 ) = JVS( 513 )
  a = -W( 39 ) / JVS(          215  )
  W( 39 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 216 )
  W( 134 ) = W( 134 ) + a*JVS( 217 )
  W( 135 ) = W( 135 ) + a*JVS( 218 )
  a = -W( 41 ) / JVS(          222  )
  W( 41 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 223 )
  JVS( 508) = W( 39 )
  JVS( 509) = W( 41 )
  JVS( 510) = W( 94 )
  JVS( 511) = W( 134 )
  JVS( 512) = W( 135 )
  JVS( 513) = W( 139 )
  IF ( ABS(  JVS( 517 )) < TINY(a) ) THEN
         IER = 95                                      
         RETURN
  END IF
   W( 34 ) = JVS( 514 )
   W( 93 ) = JVS( 515 )
   W( 94 ) = JVS( 516 )
   W( 95 ) = JVS( 517 )
   W( 134 ) = JVS( 518 )
   W( 135 ) = JVS( 519 )
   W( 138 ) = JVS( 520 )
   W( 139 ) = JVS( 521 )
  a = -W( 34 ) / JVS(          195  )
  W( 34 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 196 )
  W( 94 ) = W( 94 ) + a*JVS( 197 )
  W( 139 ) = W( 139 ) + a*JVS( 198 )
  a = -W( 93 ) / JVS(          504  )
  W( 93 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 505 )
  W( 135 ) = W( 135 ) + a*JVS( 506 )
  W( 139 ) = W( 139 ) + a*JVS( 507 )
  a = -W( 94 ) / JVS(          510  )
  W( 94 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 511 )
  W( 135 ) = W( 135 ) + a*JVS( 512 )
  W( 139 ) = W( 139 ) + a*JVS( 513 )
  JVS( 514) = W( 34 )
  JVS( 515) = W( 93 )
  JVS( 516) = W( 94 )
  JVS( 517) = W( 95 )
  JVS( 518) = W( 134 )
  JVS( 519) = W( 135 )
  JVS( 520) = W( 138 )
  JVS( 521) = W( 139 )
  IF ( ABS(  JVS( 523 )) < TINY(a) ) THEN
         IER = 96                                      
         RETURN
  END IF
   W( 84 ) = JVS( 522 )
   W( 96 ) = JVS( 523 )
   W( 113 ) = JVS( 524 )
   W( 130 ) = JVS( 525 )
   W( 134 ) = JVS( 526 )
   W( 135 ) = JVS( 527 )
   W( 136 ) = JVS( 528 )
   W( 138 ) = JVS( 529 )
  a = -W( 84 ) / JVS(          451  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 452 )
  W( 130 ) = W( 130 ) + a*JVS( 453 )
  W( 134 ) = W( 134 ) + a*JVS( 454 )
  W( 136 ) = W( 136 ) + a*JVS( 455 )
  JVS( 522) = W( 84 )
  JVS( 523) = W( 96 )
  JVS( 524) = W( 113 )
  JVS( 525) = W( 130 )
  JVS( 526) = W( 134 )
  JVS( 527) = W( 135 )
  JVS( 528) = W( 136 )
  JVS( 529) = W( 138 )
  IF ( ABS(  JVS( 532 )) < TINY(a) ) THEN
         IER = 97                                      
         RETURN
  END IF
   W( 17 ) = JVS( 530 )
   W( 49 ) = JVS( 531 )
   W( 97 ) = JVS( 532 )
   W( 134 ) = JVS( 533 )
   W( 135 ) = JVS( 534 )
   W( 139 ) = JVS( 535 )
   W( 140 ) = JVS( 536 )
  a = -W( 17 ) / JVS(          145  )
  W( 17 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 146 )
  a = -W( 49 ) / JVS(          252  )
  W( 49 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 253 )
  W( 134 ) = W( 134 ) + a*JVS( 254 )
  W( 135 ) = W( 135 ) + a*JVS( 255 )
  JVS( 530) = W( 17 )
  JVS( 531) = W( 49 )
  JVS( 532) = W( 97 )
  JVS( 533) = W( 134 )
  JVS( 534) = W( 135 )
  JVS( 535) = W( 139 )
  JVS( 536) = W( 140 )
  IF ( ABS(  JVS( 540 )) < TINY(a) ) THEN
         IER = 98                                      
         RETURN
  END IF
   W( 46 ) = JVS( 537 )
   W( 93 ) = JVS( 538 )
   W( 94 ) = JVS( 539 )
   W( 98 ) = JVS( 540 )
   W( 134 ) = JVS( 541 )
   W( 135 ) = JVS( 542 )
   W( 138 ) = JVS( 543 )
   W( 139 ) = JVS( 544 )
  a = -W( 46 ) / JVS(          237  )
  W( 46 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 238 )
  W( 94 ) = W( 94 ) + a*JVS( 239 )
  W( 139 ) = W( 139 ) + a*JVS( 240 )
  a = -W( 93 ) / JVS(          504  )
  W( 93 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 505 )
  W( 135 ) = W( 135 ) + a*JVS( 506 )
  W( 139 ) = W( 139 ) + a*JVS( 507 )
  a = -W( 94 ) / JVS(          510  )
  W( 94 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 511 )
  W( 135 ) = W( 135 ) + a*JVS( 512 )
  W( 139 ) = W( 139 ) + a*JVS( 513 )
  JVS( 537) = W( 46 )
  JVS( 538) = W( 93 )
  JVS( 539) = W( 94 )
  JVS( 540) = W( 98 )
  JVS( 541) = W( 134 )
  JVS( 542) = W( 135 )
  JVS( 543) = W( 138 )
  JVS( 544) = W( 139 )
  IF ( ABS(  JVS( 554 )) < TINY(a) ) THEN
         IER = 99                                      
         RETURN
  END IF
   W( 23 ) = JVS( 545 )
   W( 61 ) = JVS( 546 )
   W( 62 ) = JVS( 547 )
   W( 70 ) = JVS( 548 )
   W( 74 ) = JVS( 549 )
   W( 91 ) = JVS( 550 )
   W( 93 ) = JVS( 551 )
   W( 94 ) = JVS( 552 )
   W( 98 ) = JVS( 553 )
   W( 99 ) = JVS( 554 )
   W( 101 ) = JVS( 555 )
   W( 112 ) = JVS( 556 )
   W( 127 ) = JVS( 557 )
   W( 131 ) = JVS( 558 )
   W( 133 ) = JVS( 559 )
   W( 134 ) = JVS( 560 )
   W( 135 ) = JVS( 561 )
   W( 137 ) = JVS( 562 )
   W( 138 ) = JVS( 563 )
   W( 139 ) = JVS( 564 )
   W( 140 ) = JVS( 565 )
  a = -W( 23 ) / JVS(          161  )
  W( 23 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 162 )
  a = -W( 61 ) / JVS(          309  )
  W( 61 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 310 )
  W( 134 ) = W( 134 ) + a*JVS( 311 )
  W( 135 ) = W( 135 ) + a*JVS( 312 )
  a = -W( 62 ) / JVS(          316  )
  W( 62 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 317 )
  W( 135 ) = W( 135 ) + a*JVS( 318 )
  W( 139 ) = W( 139 ) + a*JVS( 319 )
  a = -W( 70 ) / JVS(          365  )
  W( 70 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 366 )
  W( 135 ) = W( 135 ) + a*JVS( 367 )
  W( 139 ) = W( 139 ) + a*JVS( 368 )
  a = -W( 74 ) / JVS(          386  )
  W( 74 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 387 )
  W( 135 ) = W( 135 ) + a*JVS( 388 )
  W( 139 ) = W( 139 ) + a*JVS( 389 )
  a = -W( 91 ) / JVS(          491  )
  W( 91 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 492 )
  W( 94 ) = W( 94 ) + a*JVS( 493 )
  W( 134 ) = W( 134 ) + a*JVS( 494 )
  W( 135 ) = W( 135 ) + a*JVS( 495 )
  W( 138 ) = W( 138 ) + a*JVS( 496 )
  W( 139 ) = W( 139 ) + a*JVS( 497 )
  a = -W( 93 ) / JVS(          504  )
  W( 93 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 505 )
  W( 135 ) = W( 135 ) + a*JVS( 506 )
  W( 139 ) = W( 139 ) + a*JVS( 507 )
  a = -W( 94 ) / JVS(          510  )
  W( 94 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 511 )
  W( 135 ) = W( 135 ) + a*JVS( 512 )
  W( 139 ) = W( 139 ) + a*JVS( 513 )
  a = -W( 98 ) / JVS(          540  )
  W( 98 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 541 )
  W( 135 ) = W( 135 ) + a*JVS( 542 )
  W( 138 ) = W( 138 ) + a*JVS( 543 )
  W( 139 ) = W( 139 ) + a*JVS( 544 )
  JVS( 545) = W( 23 )
  JVS( 546) = W( 61 )
  JVS( 547) = W( 62 )
  JVS( 548) = W( 70 )
  JVS( 549) = W( 74 )
  JVS( 550) = W( 91 )
  JVS( 551) = W( 93 )
  JVS( 552) = W( 94 )
  JVS( 553) = W( 98 )
  JVS( 554) = W( 99 )
  JVS( 555) = W( 101 )
  JVS( 556) = W( 112 )
  JVS( 557) = W( 127 )
  JVS( 558) = W( 131 )
  JVS( 559) = W( 133 )
  JVS( 560) = W( 134 )
  JVS( 561) = W( 135 )
  JVS( 562) = W( 137 )
  JVS( 563) = W( 138 )
  JVS( 564) = W( 139 )
  JVS( 565) = W( 140 )
  IF ( ABS(  JVS( 566 )) < TINY(a) ) THEN
         IER = 100                                     
         RETURN
  END IF
   W( 100 ) = JVS( 566 )
   W( 102 ) = JVS( 567 )
   W( 116 ) = JVS( 568 )
   W( 118 ) = JVS( 569 )
   W( 119 ) = JVS( 570 )
   W( 121 ) = JVS( 571 )
   W( 122 ) = JVS( 572 )
   W( 127 ) = JVS( 573 )
   W( 129 ) = JVS( 574 )
   W( 131 ) = JVS( 575 )
   W( 134 ) = JVS( 576 )
   W( 140 ) = JVS( 577 )
  JVS( 566) = W( 100 )
  JVS( 567) = W( 102 )
  JVS( 568) = W( 116 )
  JVS( 569) = W( 118 )
  JVS( 570) = W( 119 )
  JVS( 571) = W( 121 )
  JVS( 572) = W( 122 )
  JVS( 573) = W( 127 )
  JVS( 574) = W( 129 )
  JVS( 575) = W( 131 )
  JVS( 576) = W( 134 )
  JVS( 577) = W( 140 )
  IF ( ABS(  JVS( 582 )) < TINY(a) ) THEN
         IER = 101                                     
         RETURN
  END IF
   W( 61 ) = JVS( 578 )
   W( 86 ) = JVS( 579 )
   W( 87 ) = JVS( 580 )
   W( 93 ) = JVS( 581 )
   W( 101 ) = JVS( 582 )
   W( 104 ) = JVS( 583 )
   W( 107 ) = JVS( 584 )
   W( 109 ) = JVS( 585 )
   W( 134 ) = JVS( 586 )
   W( 135 ) = JVS( 587 )
   W( 136 ) = JVS( 588 )
   W( 139 ) = JVS( 589 )
  a = -W( 61 ) / JVS(          309  )
  W( 61 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 310 )
  W( 134 ) = W( 134 ) + a*JVS( 311 )
  W( 135 ) = W( 135 ) + a*JVS( 312 )
  a = -W( 86 ) / JVS(          465  )
  W( 86 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 466 )
  W( 136 ) = W( 136 ) + a*JVS( 467 )
  a = -W( 87 ) / JVS(          468  )
  W( 87 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 469 )
  W( 136 ) = W( 136 ) + a*JVS( 470 )
  a = -W( 93 ) / JVS(          504  )
  W( 93 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 505 )
  W( 135 ) = W( 135 ) + a*JVS( 506 )
  W( 139 ) = W( 139 ) + a*JVS( 507 )
  JVS( 578) = W( 61 )
  JVS( 579) = W( 86 )
  JVS( 580) = W( 87 )
  JVS( 581) = W( 93 )
  JVS( 582) = W( 101 )
  JVS( 583) = W( 104 )
  JVS( 584) = W( 107 )
  JVS( 585) = W( 109 )
  JVS( 586) = W( 134 )
  JVS( 587) = W( 135 )
  JVS( 588) = W( 136 )
  JVS( 589) = W( 139 )
  IF ( ABS(  JVS( 594 )) < TINY(a) ) THEN
         IER = 102                                     
         RETURN
  END IF
   W( 16 ) = JVS( 590 )
   W( 48 ) = JVS( 591 )
   W( 66 ) = JVS( 592 )
   W( 88 ) = JVS( 593 )
   W( 102 ) = JVS( 594 )
   W( 117 ) = JVS( 595 )
   W( 134 ) = JVS( 596 )
   W( 135 ) = JVS( 597 )
   W( 139 ) = JVS( 598 )
   W( 140 ) = JVS( 599 )
  a = -W( 16 ) / JVS(          143  )
  W( 16 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 144 )
  a = -W( 48 ) / JVS(          248  )
  W( 48 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 249 )
  W( 134 ) = W( 134 ) + a*JVS( 250 )
  W( 135 ) = W( 135 ) + a*JVS( 251 )
  a = -W( 66 ) / JVS(          340  )
  W( 66 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 341 )
  W( 117 ) = W( 117 ) + a*JVS( 342 )
  W( 134 ) = W( 134 ) + a*JVS( 343 )
  W( 139 ) = W( 139 ) + a*JVS( 344 )
  a = -W( 88 ) / JVS(          471  )
  W( 88 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 472 )
  W( 134 ) = W( 134 ) + a*JVS( 473 )
  W( 135 ) = W( 135 ) + a*JVS( 474 )
  JVS( 590) = W( 16 )
  JVS( 591) = W( 48 )
  JVS( 592) = W( 66 )
  JVS( 593) = W( 88 )
  JVS( 594) = W( 102 )
  JVS( 595) = W( 117 )
  JVS( 596) = W( 134 )
  JVS( 597) = W( 135 )
  JVS( 598) = W( 139 )
  JVS( 599) = W( 140 )
  IF ( ABS(  JVS( 600 )) < TINY(a) ) THEN
         IER = 103                                     
         RETURN
  END IF
   W( 103 ) = JVS( 600 )
   W( 108 ) = JVS( 601 )
   W( 133 ) = JVS( 602 )
   W( 135 ) = JVS( 603 )
   W( 139 ) = JVS( 604 )
  JVS( 600) = W( 103 )
  JVS( 601) = W( 108 )
  JVS( 602) = W( 133 )
  JVS( 603) = W( 135 )
  JVS( 604) = W( 139 )
  IF ( ABS(  JVS( 605 )) < TINY(a) ) THEN
         IER = 104                                     
         RETURN
  END IF
   W( 104 ) = JVS( 605 )
   W( 133 ) = JVS( 606 )
   W( 134 ) = JVS( 607 )
   W( 136 ) = JVS( 608 )
  JVS( 605) = W( 104 )
  JVS( 606) = W( 133 )
  JVS( 607) = W( 134 )
  JVS( 608) = W( 136 )
  IF ( ABS(  JVS( 610 )) < TINY(a) ) THEN
         IER = 105                                     
         RETURN
  END IF
   W( 92 ) = JVS( 609 )
   W( 105 ) = JVS( 610 )
   W( 133 ) = JVS( 611 )
   W( 134 ) = JVS( 612 )
   W( 135 ) = JVS( 613 )
   W( 136 ) = JVS( 614 )
   W( 139 ) = JVS( 615 )
  a = -W( 92 ) / JVS(          498  )
  W( 92 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 499 )
  W( 134 ) = W( 134 ) + a*JVS( 500 )
  W( 136 ) = W( 136 ) + a*JVS( 501 )
  JVS( 609) = W( 92 )
  JVS( 610) = W( 105 )
  JVS( 611) = W( 133 )
  JVS( 612) = W( 134 )
  JVS( 613) = W( 135 )
  JVS( 614) = W( 136 )
  JVS( 615) = W( 139 )
  IF ( ABS(  JVS( 619 )) < TINY(a) ) THEN
         IER = 106                                     
         RETURN
  END IF
   W( 79 ) = JVS( 616 )
   W( 92 ) = JVS( 617 )
   W( 105 ) = JVS( 618 )
   W( 106 ) = JVS( 619 )
   W( 116 ) = JVS( 620 )
   W( 133 ) = JVS( 621 )
   W( 134 ) = JVS( 622 )
   W( 135 ) = JVS( 623 )
   W( 136 ) = JVS( 624 )
   W( 139 ) = JVS( 625 )
   W( 140 ) = JVS( 626 )
  a = -W( 79 ) / JVS(          417  )
  W( 79 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 418 )
  W( 134 ) = W( 134 ) + a*JVS( 419 )
  W( 135 ) = W( 135 ) + a*JVS( 420 )
  a = -W( 92 ) / JVS(          498  )
  W( 92 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 499 )
  W( 134 ) = W( 134 ) + a*JVS( 500 )
  W( 136 ) = W( 136 ) + a*JVS( 501 )
  a = -W( 105 ) / JVS(          610  )
  W( 105 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 611 )
  W( 134 ) = W( 134 ) + a*JVS( 612 )
  W( 135 ) = W( 135 ) + a*JVS( 613 )
  W( 136 ) = W( 136 ) + a*JVS( 614 )
  W( 139 ) = W( 139 ) + a*JVS( 615 )
  JVS( 616) = W( 79 )
  JVS( 617) = W( 92 )
  JVS( 618) = W( 105 )
  JVS( 619) = W( 106 )
  JVS( 620) = W( 116 )
  JVS( 621) = W( 133 )
  JVS( 622) = W( 134 )
  JVS( 623) = W( 135 )
  JVS( 624) = W( 136 )
  JVS( 625) = W( 139 )
  JVS( 626) = W( 140 )
  IF ( ABS(  JVS( 627 )) < TINY(a) ) THEN
         IER = 107                                     
         RETURN
  END IF
   W( 107 ) = JVS( 627 )
   W( 133 ) = JVS( 628 )
   W( 134 ) = JVS( 629 )
   W( 136 ) = JVS( 630 )
  JVS( 627) = W( 107 )
  JVS( 628) = W( 133 )
  JVS( 629) = W( 134 )
  JVS( 630) = W( 136 )
  IF ( ABS(  JVS( 631 )) < TINY(a) ) THEN
         IER = 108                                     
         RETURN
  END IF
   W( 108 ) = JVS( 631 )
   W( 133 ) = JVS( 632 )
   W( 134 ) = JVS( 633 )
   W( 136 ) = JVS( 634 )
  JVS( 631) = W( 108 )
  JVS( 632) = W( 133 )
  JVS( 633) = W( 134 )
  JVS( 634) = W( 136 )
  IF ( ABS(  JVS( 635 )) < TINY(a) ) THEN
         IER = 109                                     
         RETURN
  END IF
   W( 109 ) = JVS( 635 )
   W( 133 ) = JVS( 636 )
   W( 134 ) = JVS( 637 )
   W( 136 ) = JVS( 638 )
  JVS( 635) = W( 109 )
  JVS( 636) = W( 133 )
  JVS( 637) = W( 134 )
  JVS( 638) = W( 136 )
  IF ( ABS(  JVS( 655 )) < TINY(a) ) THEN
         IER = 110                                     
         RETURN
  END IF
   W( 49 ) = JVS( 639 )
   W( 51 ) = JVS( 640 )
   W( 69 ) = JVS( 641 )
   W( 76 ) = JVS( 642 )
   W( 79 ) = JVS( 643 )
   W( 83 ) = JVS( 644 )
   W( 86 ) = JVS( 645 )
   W( 87 ) = JVS( 646 )
   W( 88 ) = JVS( 647 )
   W( 92 ) = JVS( 648 )
   W( 97 ) = JVS( 649 )
   W( 104 ) = JVS( 650 )
   W( 105 ) = JVS( 651 )
   W( 106 ) = JVS( 652 )
   W( 107 ) = JVS( 653 )
   W( 109 ) = JVS( 654 )
   W( 110 ) = JVS( 655 )
   W( 115 ) = JVS( 656 )
   W( 116 ) = JVS( 657 )
   W( 117 ) = JVS( 658 )
   W( 118 ) = JVS( 659 )
   W( 119 ) = JVS( 660 )
   W( 133 ) = JVS( 661 )
   W( 134 ) = JVS( 662 )
   W( 135 ) = JVS( 663 )
   W( 136 ) = JVS( 664 )
   W( 139 ) = JVS( 665 )
   W( 140 ) = JVS( 666 )
  a = -W( 49 ) / JVS(          252  )
  W( 49 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 253 )
  W( 134 ) = W( 134 ) + a*JVS( 254 )
  W( 135 ) = W( 135 ) + a*JVS( 255 )
  a = -W( 51 ) / JVS(          261  )
  W( 51 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 262 )
  W( 139 ) = W( 139 ) + a*JVS( 263 )
  a = -W( 69 ) / JVS(          359  )
  W( 69 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 360 )
  W( 134 ) = W( 134 ) + a*JVS( 361 )
  W( 135 ) = W( 135 ) + a*JVS( 362 )
  a = -W( 76 ) / JVS(          399  )
  W( 76 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 400 )
  W( 134 ) = W( 134 ) + a*JVS( 401 )
  W( 135 ) = W( 135 ) + a*JVS( 402 )
  W( 139 ) = W( 139 ) + a*JVS( 403 )
  a = -W( 79 ) / JVS(          417  )
  W( 79 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 418 )
  W( 134 ) = W( 134 ) + a*JVS( 419 )
  W( 135 ) = W( 135 ) + a*JVS( 420 )
  a = -W( 83 ) / JVS(          444  )
  W( 83 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 445 )
  W( 134 ) = W( 134 ) + a*JVS( 446 )
  W( 135 ) = W( 135 ) + a*JVS( 447 )
  a = -W( 86 ) / JVS(          465  )
  W( 86 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 466 )
  W( 136 ) = W( 136 ) + a*JVS( 467 )
  a = -W( 87 ) / JVS(          468  )
  W( 87 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 469 )
  W( 136 ) = W( 136 ) + a*JVS( 470 )
  a = -W( 88 ) / JVS(          471  )
  W( 88 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 472 )
  W( 134 ) = W( 134 ) + a*JVS( 473 )
  W( 135 ) = W( 135 ) + a*JVS( 474 )
  a = -W( 92 ) / JVS(          498  )
  W( 92 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 499 )
  W( 134 ) = W( 134 ) + a*JVS( 500 )
  W( 136 ) = W( 136 ) + a*JVS( 501 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 533 )
  W( 135 ) = W( 135 ) + a*JVS( 534 )
  W( 139 ) = W( 139 ) + a*JVS( 535 )
  W( 140 ) = W( 140 ) + a*JVS( 536 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 105 ) / JVS(          610  )
  W( 105 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 611 )
  W( 134 ) = W( 134 ) + a*JVS( 612 )
  W( 135 ) = W( 135 ) + a*JVS( 613 )
  W( 136 ) = W( 136 ) + a*JVS( 614 )
  W( 139 ) = W( 139 ) + a*JVS( 615 )
  a = -W( 106 ) / JVS(          619  )
  W( 106 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 620 )
  W( 133 ) = W( 133 ) + a*JVS( 621 )
  W( 134 ) = W( 134 ) + a*JVS( 622 )
  W( 135 ) = W( 135 ) + a*JVS( 623 )
  W( 136 ) = W( 136 ) + a*JVS( 624 )
  W( 139 ) = W( 139 ) + a*JVS( 625 )
  W( 140 ) = W( 140 ) + a*JVS( 626 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  JVS( 639) = W( 49 )
  JVS( 640) = W( 51 )
  JVS( 641) = W( 69 )
  JVS( 642) = W( 76 )
  JVS( 643) = W( 79 )
  JVS( 644) = W( 83 )
  JVS( 645) = W( 86 )
  JVS( 646) = W( 87 )
  JVS( 647) = W( 88 )
  JVS( 648) = W( 92 )
  JVS( 649) = W( 97 )
  JVS( 650) = W( 104 )
  JVS( 651) = W( 105 )
  JVS( 652) = W( 106 )
  JVS( 653) = W( 107 )
  JVS( 654) = W( 109 )
  JVS( 655) = W( 110 )
  JVS( 656) = W( 115 )
  JVS( 657) = W( 116 )
  JVS( 658) = W( 117 )
  JVS( 659) = W( 118 )
  JVS( 660) = W( 119 )
  JVS( 661) = W( 133 )
  JVS( 662) = W( 134 )
  JVS( 663) = W( 135 )
  JVS( 664) = W( 136 )
  JVS( 665) = W( 139 )
  JVS( 666) = W( 140 )
  IF ( ABS(  JVS( 690 )) < TINY(a) ) THEN
         IER = 111                                     
         RETURN
  END IF
   W( 23 ) = JVS( 667 )
   W( 25 ) = JVS( 668 )
   W( 35 ) = JVS( 669 )
   W( 41 ) = JVS( 670 )
   W( 44 ) = JVS( 671 )
   W( 46 ) = JVS( 672 )
   W( 59 ) = JVS( 673 )
   W( 83 ) = JVS( 674 )
   W( 86 ) = JVS( 675 )
   W( 87 ) = JVS( 676 )
   W( 91 ) = JVS( 677 )
   W( 92 ) = JVS( 678 )
   W( 93 ) = JVS( 679 )
   W( 94 ) = JVS( 680 )
   W( 95 ) = JVS( 681 )
   W( 98 ) = JVS( 682 )
   W( 99 ) = JVS( 683 )
   W( 101 ) = JVS( 684 )
   W( 104 ) = JVS( 685 )
   W( 106 ) = JVS( 686 )
   W( 107 ) = JVS( 687 )
   W( 108 ) = JVS( 688 )
   W( 109 ) = JVS( 689 )
   W( 111 ) = JVS( 690 )
   W( 112 ) = JVS( 691 )
   W( 113 ) = JVS( 692 )
   W( 114 ) = JVS( 693 )
   W( 115 ) = JVS( 694 )
   W( 116 ) = JVS( 695 )
   W( 118 ) = JVS( 696 )
   W( 120 ) = JVS( 697 )
   W( 123 ) = JVS( 698 )
   W( 124 ) = JVS( 699 )
   W( 126 ) = JVS( 700 )
   W( 127 ) = JVS( 701 )
   W( 128 ) = JVS( 702 )
   W( 129 ) = JVS( 703 )
   W( 130 ) = JVS( 704 )
   W( 131 ) = JVS( 705 )
   W( 133 ) = JVS( 706 )
   W( 134 ) = JVS( 707 )
   W( 135 ) = JVS( 708 )
   W( 136 ) = JVS( 709 )
   W( 137 ) = JVS( 710 )
   W( 138 ) = JVS( 711 )
   W( 139 ) = JVS( 712 )
   W( 140 ) = JVS( 713 )
  a = -W( 23 ) / JVS(          161  )
  W( 23 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 162 )
  a = -W( 25 ) / JVS(          166  )
  W( 25 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 167 )
  a = -W( 35 ) / JVS(          200  )
  W( 35 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 201 )
  W( 134 ) = W( 134 ) + a*JVS( 202 )
  a = -W( 41 ) / JVS(          222  )
  W( 41 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 223 )
  a = -W( 44 ) / JVS(          230  )
  W( 44 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 231 )
  W( 139 ) = W( 139 ) + a*JVS( 232 )
  a = -W( 46 ) / JVS(          237  )
  W( 46 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 238 )
  W( 94 ) = W( 94 ) + a*JVS( 239 )
  W( 139 ) = W( 139 ) + a*JVS( 240 )
  a = -W( 59 ) / JVS(          301  )
  W( 59 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 302 )
  W( 136 ) = W( 136 ) + a*JVS( 303 )
  a = -W( 83 ) / JVS(          444  )
  W( 83 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 445 )
  W( 134 ) = W( 134 ) + a*JVS( 446 )
  W( 135 ) = W( 135 ) + a*JVS( 447 )
  a = -W( 86 ) / JVS(          465  )
  W( 86 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 466 )
  W( 136 ) = W( 136 ) + a*JVS( 467 )
  a = -W( 87 ) / JVS(          468  )
  W( 87 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 469 )
  W( 136 ) = W( 136 ) + a*JVS( 470 )
  a = -W( 91 ) / JVS(          491  )
  W( 91 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 492 )
  W( 94 ) = W( 94 ) + a*JVS( 493 )
  W( 134 ) = W( 134 ) + a*JVS( 494 )
  W( 135 ) = W( 135 ) + a*JVS( 495 )
  W( 138 ) = W( 138 ) + a*JVS( 496 )
  W( 139 ) = W( 139 ) + a*JVS( 497 )
  a = -W( 92 ) / JVS(          498  )
  W( 92 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 499 )
  W( 134 ) = W( 134 ) + a*JVS( 500 )
  W( 136 ) = W( 136 ) + a*JVS( 501 )
  a = -W( 93 ) / JVS(          504  )
  W( 93 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 505 )
  W( 135 ) = W( 135 ) + a*JVS( 506 )
  W( 139 ) = W( 139 ) + a*JVS( 507 )
  a = -W( 94 ) / JVS(          510  )
  W( 94 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 511 )
  W( 135 ) = W( 135 ) + a*JVS( 512 )
  W( 139 ) = W( 139 ) + a*JVS( 513 )
  a = -W( 95 ) / JVS(          517  )
  W( 95 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 518 )
  W( 135 ) = W( 135 ) + a*JVS( 519 )
  W( 138 ) = W( 138 ) + a*JVS( 520 )
  W( 139 ) = W( 139 ) + a*JVS( 521 )
  a = -W( 98 ) / JVS(          540  )
  W( 98 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 541 )
  W( 135 ) = W( 135 ) + a*JVS( 542 )
  W( 138 ) = W( 138 ) + a*JVS( 543 )
  W( 139 ) = W( 139 ) + a*JVS( 544 )
  a = -W( 99 ) / JVS(          554  )
  W( 99 ) = -a
  W( 101 ) = W( 101 ) + a*JVS( 555 )
  W( 112 ) = W( 112 ) + a*JVS( 556 )
  W( 127 ) = W( 127 ) + a*JVS( 557 )
  W( 131 ) = W( 131 ) + a*JVS( 558 )
  W( 133 ) = W( 133 ) + a*JVS( 559 )
  W( 134 ) = W( 134 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 137 ) = W( 137 ) + a*JVS( 562 )
  W( 138 ) = W( 138 ) + a*JVS( 563 )
  W( 139 ) = W( 139 ) + a*JVS( 564 )
  W( 140 ) = W( 140 ) + a*JVS( 565 )
  a = -W( 101 ) / JVS(          582  )
  W( 101 ) = -a
  W( 104 ) = W( 104 ) + a*JVS( 583 )
  W( 107 ) = W( 107 ) + a*JVS( 584 )
  W( 109 ) = W( 109 ) + a*JVS( 585 )
  W( 134 ) = W( 134 ) + a*JVS( 586 )
  W( 135 ) = W( 135 ) + a*JVS( 587 )
  W( 136 ) = W( 136 ) + a*JVS( 588 )
  W( 139 ) = W( 139 ) + a*JVS( 589 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 106 ) / JVS(          619  )
  W( 106 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 620 )
  W( 133 ) = W( 133 ) + a*JVS( 621 )
  W( 134 ) = W( 134 ) + a*JVS( 622 )
  W( 135 ) = W( 135 ) + a*JVS( 623 )
  W( 136 ) = W( 136 ) + a*JVS( 624 )
  W( 139 ) = W( 139 ) + a*JVS( 625 )
  W( 140 ) = W( 140 ) + a*JVS( 626 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  JVS( 667) = W( 23 )
  JVS( 668) = W( 25 )
  JVS( 669) = W( 35 )
  JVS( 670) = W( 41 )
  JVS( 671) = W( 44 )
  JVS( 672) = W( 46 )
  JVS( 673) = W( 59 )
  JVS( 674) = W( 83 )
  JVS( 675) = W( 86 )
  JVS( 676) = W( 87 )
  JVS( 677) = W( 91 )
  JVS( 678) = W( 92 )
  JVS( 679) = W( 93 )
  JVS( 680) = W( 94 )
  JVS( 681) = W( 95 )
  JVS( 682) = W( 98 )
  JVS( 683) = W( 99 )
  JVS( 684) = W( 101 )
  JVS( 685) = W( 104 )
  JVS( 686) = W( 106 )
  JVS( 687) = W( 107 )
  JVS( 688) = W( 108 )
  JVS( 689) = W( 109 )
  JVS( 690) = W( 111 )
  JVS( 691) = W( 112 )
  JVS( 692) = W( 113 )
  JVS( 693) = W( 114 )
  JVS( 694) = W( 115 )
  JVS( 695) = W( 116 )
  JVS( 696) = W( 118 )
  JVS( 697) = W( 120 )
  JVS( 698) = W( 123 )
  JVS( 699) = W( 124 )
  JVS( 700) = W( 126 )
  JVS( 701) = W( 127 )
  JVS( 702) = W( 128 )
  JVS( 703) = W( 129 )
  JVS( 704) = W( 130 )
  JVS( 705) = W( 131 )
  JVS( 706) = W( 133 )
  JVS( 707) = W( 134 )
  JVS( 708) = W( 135 )
  JVS( 709) = W( 136 )
  JVS( 710) = W( 137 )
  JVS( 711) = W( 138 )
  JVS( 712) = W( 139 )
  JVS( 713) = W( 140 )
  IF ( ABS(  JVS( 718 )) < TINY(a) ) THEN
         IER = 112                                     
         RETURN
  END IF
   W( 30 ) = JVS( 714 )
   W( 60 ) = JVS( 715 )
   W( 79 ) = JVS( 716 )
   W( 83 ) = JVS( 717 )
   W( 112 ) = JVS( 718 )
   W( 116 ) = JVS( 719 )
   W( 118 ) = JVS( 720 )
   W( 127 ) = JVS( 721 )
   W( 129 ) = JVS( 722 )
   W( 131 ) = JVS( 723 )
   W( 133 ) = JVS( 724 )
   W( 134 ) = JVS( 725 )
   W( 135 ) = JVS( 726 )
   W( 136 ) = JVS( 727 )
   W( 137 ) = JVS( 728 )
   W( 139 ) = JVS( 729 )
   W( 140 ) = JVS( 730 )
  a = -W( 30 ) / JVS(          180  )
  W( 30 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 181 )
  W( 139 ) = W( 139 ) + a*JVS( 182 )
  a = -W( 60 ) / JVS(          305  )
  W( 60 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 306 )
  W( 136 ) = W( 136 ) + a*JVS( 307 )
  W( 139 ) = W( 139 ) + a*JVS( 308 )
  a = -W( 79 ) / JVS(          417  )
  W( 79 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 418 )
  W( 134 ) = W( 134 ) + a*JVS( 419 )
  W( 135 ) = W( 135 ) + a*JVS( 420 )
  a = -W( 83 ) / JVS(          444  )
  W( 83 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 445 )
  W( 134 ) = W( 134 ) + a*JVS( 446 )
  W( 135 ) = W( 135 ) + a*JVS( 447 )
  JVS( 714) = W( 30 )
  JVS( 715) = W( 60 )
  JVS( 716) = W( 79 )
  JVS( 717) = W( 83 )
  JVS( 718) = W( 112 )
  JVS( 719) = W( 116 )
  JVS( 720) = W( 118 )
  JVS( 721) = W( 127 )
  JVS( 722) = W( 129 )
  JVS( 723) = W( 131 )
  JVS( 724) = W( 133 )
  JVS( 725) = W( 134 )
  JVS( 726) = W( 135 )
  JVS( 727) = W( 136 )
  JVS( 728) = W( 137 )
  JVS( 729) = W( 139 )
  JVS( 730) = W( 140 )
  IF ( ABS(  JVS( 732 )) < TINY(a) ) THEN
         IER = 113                                     
         RETURN
  END IF
   W( 108 ) = JVS( 731 )
   W( 113 ) = JVS( 732 )
   W( 128 ) = JVS( 733 )
   W( 133 ) = JVS( 734 )
   W( 134 ) = JVS( 735 )
   W( 136 ) = JVS( 736 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  JVS( 731) = W( 108 )
  JVS( 732) = W( 113 )
  JVS( 733) = W( 128 )
  JVS( 734) = W( 133 )
  JVS( 735) = W( 134 )
  JVS( 736) = W( 136 )
  IF ( ABS(  JVS( 748 )) < TINY(a) ) THEN
         IER = 114                                     
         RETURN
  END IF
   W( 45 ) = JVS( 737 )
   W( 48 ) = JVS( 738 )
   W( 51 ) = JVS( 739 )
   W( 52 ) = JVS( 740 )
   W( 72 ) = JVS( 741 )
   W( 88 ) = JVS( 742 )
   W( 89 ) = JVS( 743 )
   W( 90 ) = JVS( 744 )
   W( 97 ) = JVS( 745 )
   W( 102 ) = JVS( 746 )
   W( 113 ) = JVS( 747 )
   W( 114 ) = JVS( 748 )
   W( 117 ) = JVS( 749 )
   W( 128 ) = JVS( 750 )
   W( 133 ) = JVS( 751 )
   W( 134 ) = JVS( 752 )
   W( 135 ) = JVS( 753 )
   W( 136 ) = JVS( 754 )
   W( 139 ) = JVS( 755 )
   W( 140 ) = JVS( 756 )
  a = -W( 45 ) / JVS(          233  )
  W( 45 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 234 )
  W( 134 ) = W( 134 ) + a*JVS( 235 )
  W( 140 ) = W( 140 ) + a*JVS( 236 )
  a = -W( 48 ) / JVS(          248  )
  W( 48 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 249 )
  W( 134 ) = W( 134 ) + a*JVS( 250 )
  W( 135 ) = W( 135 ) + a*JVS( 251 )
  a = -W( 51 ) / JVS(          261  )
  W( 51 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 262 )
  W( 139 ) = W( 139 ) + a*JVS( 263 )
  a = -W( 52 ) / JVS(          264  )
  W( 52 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 265 )
  W( 134 ) = W( 134 ) + a*JVS( 266 )
  W( 135 ) = W( 135 ) + a*JVS( 267 )
  a = -W( 72 ) / JVS(          376  )
  W( 72 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 377 )
  W( 134 ) = W( 134 ) + a*JVS( 378 )
  W( 135 ) = W( 135 ) + a*JVS( 379 )
  a = -W( 88 ) / JVS(          471  )
  W( 88 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 472 )
  W( 134 ) = W( 134 ) + a*JVS( 473 )
  W( 135 ) = W( 135 ) + a*JVS( 474 )
  a = -W( 89 ) / JVS(          478  )
  W( 89 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 479 )
  W( 134 ) = W( 134 ) + a*JVS( 480 )
  W( 135 ) = W( 135 ) + a*JVS( 481 )
  W( 139 ) = W( 139 ) + a*JVS( 482 )
  a = -W( 90 ) / JVS(          484  )
  W( 90 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 485 )
  W( 134 ) = W( 134 ) + a*JVS( 486 )
  W( 135 ) = W( 135 ) + a*JVS( 487 )
  W( 139 ) = W( 139 ) + a*JVS( 488 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 533 )
  W( 135 ) = W( 135 ) + a*JVS( 534 )
  W( 139 ) = W( 139 ) + a*JVS( 535 )
  W( 140 ) = W( 140 ) + a*JVS( 536 )
  a = -W( 102 ) / JVS(          594  )
  W( 102 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 595 )
  W( 134 ) = W( 134 ) + a*JVS( 596 )
  W( 135 ) = W( 135 ) + a*JVS( 597 )
  W( 139 ) = W( 139 ) + a*JVS( 598 )
  W( 140 ) = W( 140 ) + a*JVS( 599 )
  a = -W( 113 ) / JVS(          732  )
  W( 113 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 733 )
  W( 133 ) = W( 133 ) + a*JVS( 734 )
  W( 134 ) = W( 134 ) + a*JVS( 735 )
  W( 136 ) = W( 136 ) + a*JVS( 736 )
  JVS( 737) = W( 45 )
  JVS( 738) = W( 48 )
  JVS( 739) = W( 51 )
  JVS( 740) = W( 52 )
  JVS( 741) = W( 72 )
  JVS( 742) = W( 88 )
  JVS( 743) = W( 89 )
  JVS( 744) = W( 90 )
  JVS( 745) = W( 97 )
  JVS( 746) = W( 102 )
  JVS( 747) = W( 113 )
  JVS( 748) = W( 114 )
  JVS( 749) = W( 117 )
  JVS( 750) = W( 128 )
  JVS( 751) = W( 133 )
  JVS( 752) = W( 134 )
  JVS( 753) = W( 135 )
  JVS( 754) = W( 136 )
  JVS( 755) = W( 139 )
  JVS( 756) = W( 140 )
  IF ( ABS(  JVS( 763 )) < TINY(a) ) THEN
         IER = 115                                     
         RETURN
  END IF
   W( 83 ) = JVS( 757 )
   W( 86 ) = JVS( 758 )
   W( 87 ) = JVS( 759 )
   W( 104 ) = JVS( 760 )
   W( 107 ) = JVS( 761 )
   W( 109 ) = JVS( 762 )
   W( 115 ) = JVS( 763 )
   W( 118 ) = JVS( 764 )
   W( 124 ) = JVS( 765 )
   W( 133 ) = JVS( 766 )
   W( 134 ) = JVS( 767 )
   W( 135 ) = JVS( 768 )
   W( 136 ) = JVS( 769 )
   W( 139 ) = JVS( 770 )
   W( 140 ) = JVS( 771 )
  a = -W( 83 ) / JVS(          444  )
  W( 83 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 445 )
  W( 134 ) = W( 134 ) + a*JVS( 446 )
  W( 135 ) = W( 135 ) + a*JVS( 447 )
  a = -W( 86 ) / JVS(          465  )
  W( 86 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 466 )
  W( 136 ) = W( 136 ) + a*JVS( 467 )
  a = -W( 87 ) / JVS(          468  )
  W( 87 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 469 )
  W( 136 ) = W( 136 ) + a*JVS( 470 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  JVS( 757) = W( 83 )
  JVS( 758) = W( 86 )
  JVS( 759) = W( 87 )
  JVS( 760) = W( 104 )
  JVS( 761) = W( 107 )
  JVS( 762) = W( 109 )
  JVS( 763) = W( 115 )
  JVS( 764) = W( 118 )
  JVS( 765) = W( 124 )
  JVS( 766) = W( 133 )
  JVS( 767) = W( 134 )
  JVS( 768) = W( 135 )
  JVS( 769) = W( 136 )
  JVS( 770) = W( 139 )
  JVS( 771) = W( 140 )
  IF ( ABS(  JVS( 774 )) < TINY(a) ) THEN
         IER = 116                                     
         RETURN
  END IF
   W( 79 ) = JVS( 772 )
   W( 92 ) = JVS( 773 )
   W( 116 ) = JVS( 774 )
   W( 133 ) = JVS( 775 )
   W( 134 ) = JVS( 776 )
   W( 135 ) = JVS( 777 )
   W( 136 ) = JVS( 778 )
   W( 139 ) = JVS( 779 )
   W( 140 ) = JVS( 780 )
  a = -W( 79 ) / JVS(          417  )
  W( 79 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 418 )
  W( 134 ) = W( 134 ) + a*JVS( 419 )
  W( 135 ) = W( 135 ) + a*JVS( 420 )
  a = -W( 92 ) / JVS(          498  )
  W( 92 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 499 )
  W( 134 ) = W( 134 ) + a*JVS( 500 )
  W( 136 ) = W( 136 ) + a*JVS( 501 )
  JVS( 772) = W( 79 )
  JVS( 773) = W( 92 )
  JVS( 774) = W( 116 )
  JVS( 775) = W( 133 )
  JVS( 776) = W( 134 )
  JVS( 777) = W( 135 )
  JVS( 778) = W( 136 )
  JVS( 779) = W( 139 )
  JVS( 780) = W( 140 )
  IF ( ABS(  JVS( 788 )) < TINY(a) ) THEN
         IER = 117                                     
         RETURN
  END IF
   W( 58 ) = JVS( 781 )
   W( 86 ) = JVS( 782 )
   W( 87 ) = JVS( 783 )
   W( 88 ) = JVS( 784 )
   W( 104 ) = JVS( 785 )
   W( 107 ) = JVS( 786 )
   W( 109 ) = JVS( 787 )
   W( 117 ) = JVS( 788 )
   W( 133 ) = JVS( 789 )
   W( 134 ) = JVS( 790 )
   W( 135 ) = JVS( 791 )
   W( 136 ) = JVS( 792 )
   W( 139 ) = JVS( 793 )
  a = -W( 58 ) / JVS(          293  )
  W( 58 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 294 )
  W( 87 ) = W( 87 ) + a*JVS( 295 )
  W( 104 ) = W( 104 ) + a*JVS( 296 )
  W( 107 ) = W( 107 ) + a*JVS( 297 )
  W( 109 ) = W( 109 ) + a*JVS( 298 )
  W( 134 ) = W( 134 ) + a*JVS( 299 )
  W( 136 ) = W( 136 ) + a*JVS( 300 )
  a = -W( 86 ) / JVS(          465  )
  W( 86 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 466 )
  W( 136 ) = W( 136 ) + a*JVS( 467 )
  a = -W( 87 ) / JVS(          468  )
  W( 87 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 469 )
  W( 136 ) = W( 136 ) + a*JVS( 470 )
  a = -W( 88 ) / JVS(          471  )
  W( 88 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 472 )
  W( 134 ) = W( 134 ) + a*JVS( 473 )
  W( 135 ) = W( 135 ) + a*JVS( 474 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  JVS( 781) = W( 58 )
  JVS( 782) = W( 86 )
  JVS( 783) = W( 87 )
  JVS( 784) = W( 88 )
  JVS( 785) = W( 104 )
  JVS( 786) = W( 107 )
  JVS( 787) = W( 109 )
  JVS( 788) = W( 117 )
  JVS( 789) = W( 133 )
  JVS( 790) = W( 134 )
  JVS( 791) = W( 135 )
  JVS( 792) = W( 136 )
  JVS( 793) = W( 139 )
  IF ( ABS(  JVS( 795 )) < TINY(a) ) THEN
         IER = 118                                     
         RETURN
  END IF
   W( 83 ) = JVS( 794 )
   W( 118 ) = JVS( 795 )
   W( 124 ) = JVS( 796 )
   W( 133 ) = JVS( 797 )
   W( 134 ) = JVS( 798 )
   W( 135 ) = JVS( 799 )
   W( 139 ) = JVS( 800 )
   W( 140 ) = JVS( 801 )
  a = -W( 83 ) / JVS(          444  )
  W( 83 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 445 )
  W( 134 ) = W( 134 ) + a*JVS( 446 )
  W( 135 ) = W( 135 ) + a*JVS( 447 )
  JVS( 794) = W( 83 )
  JVS( 795) = W( 118 )
  JVS( 796) = W( 124 )
  JVS( 797) = W( 133 )
  JVS( 798) = W( 134 )
  JVS( 799) = W( 135 )
  JVS( 800) = W( 139 )
  JVS( 801) = W( 140 )
  IF ( ABS(  JVS( 808 )) < TINY(a) ) THEN
         IER = 119                                     
         RETURN
  END IF
   W( 69 ) = JVS( 802 )
   W( 86 ) = JVS( 803 )
   W( 87 ) = JVS( 804 )
   W( 104 ) = JVS( 805 )
   W( 107 ) = JVS( 806 )
   W( 109 ) = JVS( 807 )
   W( 119 ) = JVS( 808 )
   W( 133 ) = JVS( 809 )
   W( 134 ) = JVS( 810 )
   W( 135 ) = JVS( 811 )
   W( 136 ) = JVS( 812 )
   W( 139 ) = JVS( 813 )
   W( 140 ) = JVS( 814 )
  a = -W( 69 ) / JVS(          359  )
  W( 69 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 360 )
  W( 134 ) = W( 134 ) + a*JVS( 361 )
  W( 135 ) = W( 135 ) + a*JVS( 362 )
  a = -W( 86 ) / JVS(          465  )
  W( 86 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 466 )
  W( 136 ) = W( 136 ) + a*JVS( 467 )
  a = -W( 87 ) / JVS(          468  )
  W( 87 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 469 )
  W( 136 ) = W( 136 ) + a*JVS( 470 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  JVS( 802) = W( 69 )
  JVS( 803) = W( 86 )
  JVS( 804) = W( 87 )
  JVS( 805) = W( 104 )
  JVS( 806) = W( 107 )
  JVS( 807) = W( 109 )
  JVS( 808) = W( 119 )
  JVS( 809) = W( 133 )
  JVS( 810) = W( 134 )
  JVS( 811) = W( 135 )
  JVS( 812) = W( 136 )
  JVS( 813) = W( 139 )
  JVS( 814) = W( 140 )
  IF ( ABS(  JVS( 819 )) < TINY(a) ) THEN
         IER = 120                                     
         RETURN
  END IF
   W( 103 ) = JVS( 815 )
   W( 108 ) = JVS( 816 )
   W( 118 ) = JVS( 817 )
   W( 119 ) = JVS( 818 )
   W( 120 ) = JVS( 819 )
   W( 122 ) = JVS( 820 )
   W( 124 ) = JVS( 821 )
   W( 129 ) = JVS( 822 )
   W( 131 ) = JVS( 823 )
   W( 133 ) = JVS( 824 )
   W( 134 ) = JVS( 825 )
   W( 135 ) = JVS( 826 )
   W( 136 ) = JVS( 827 )
   W( 139 ) = JVS( 828 )
   W( 140 ) = JVS( 829 )
  a = -W( 103 ) / JVS(          600  )
  W( 103 ) = -a
  W( 108 ) = W( 108 ) + a*JVS( 601 )
  W( 133 ) = W( 133 ) + a*JVS( 602 )
  W( 135 ) = W( 135 ) + a*JVS( 603 )
  W( 139 ) = W( 139 ) + a*JVS( 604 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  a = -W( 118 ) / JVS(          795  )
  W( 118 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 796 )
  W( 133 ) = W( 133 ) + a*JVS( 797 )
  W( 134 ) = W( 134 ) + a*JVS( 798 )
  W( 135 ) = W( 135 ) + a*JVS( 799 )
  W( 139 ) = W( 139 ) + a*JVS( 800 )
  W( 140 ) = W( 140 ) + a*JVS( 801 )
  a = -W( 119 ) / JVS(          808  )
  W( 119 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 809 )
  W( 134 ) = W( 134 ) + a*JVS( 810 )
  W( 135 ) = W( 135 ) + a*JVS( 811 )
  W( 136 ) = W( 136 ) + a*JVS( 812 )
  W( 139 ) = W( 139 ) + a*JVS( 813 )
  W( 140 ) = W( 140 ) + a*JVS( 814 )
  JVS( 815) = W( 103 )
  JVS( 816) = W( 108 )
  JVS( 817) = W( 118 )
  JVS( 818) = W( 119 )
  JVS( 819) = W( 120 )
  JVS( 820) = W( 122 )
  JVS( 821) = W( 124 )
  JVS( 822) = W( 129 )
  JVS( 823) = W( 131 )
  JVS( 824) = W( 133 )
  JVS( 825) = W( 134 )
  JVS( 826) = W( 135 )
  JVS( 827) = W( 136 )
  JVS( 828) = W( 139 )
  JVS( 829) = W( 140 )
  IF ( ABS(  JVS( 842 )) < TINY(a) ) THEN
         IER = 121                                     
         RETURN
  END IF
   W( 50 ) = JVS( 830 )
   W( 86 ) = JVS( 831 )
   W( 87 ) = JVS( 832 )
   W( 104 ) = JVS( 833 )
   W( 107 ) = JVS( 834 )
   W( 109 ) = JVS( 835 )
   W( 110 ) = JVS( 836 )
   W( 115 ) = JVS( 837 )
   W( 116 ) = JVS( 838 )
   W( 117 ) = JVS( 839 )
   W( 118 ) = JVS( 840 )
   W( 119 ) = JVS( 841 )
   W( 121 ) = JVS( 842 )
   W( 124 ) = JVS( 843 )
   W( 133 ) = JVS( 844 )
   W( 134 ) = JVS( 845 )
   W( 135 ) = JVS( 846 )
   W( 136 ) = JVS( 847 )
   W( 139 ) = JVS( 848 )
   W( 140 ) = JVS( 849 )
  a = -W( 50 ) / JVS(          256  )
  W( 50 ) = -a
  W( 121 ) = W( 121 ) + a*JVS( 257 )
  W( 134 ) = W( 134 ) + a*JVS( 258 )
  W( 135 ) = W( 135 ) + a*JVS( 259 )
  a = -W( 86 ) / JVS(          465  )
  W( 86 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 466 )
  W( 136 ) = W( 136 ) + a*JVS( 467 )
  a = -W( 87 ) / JVS(          468  )
  W( 87 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 469 )
  W( 136 ) = W( 136 ) + a*JVS( 470 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  a = -W( 110 ) / JVS(          655  )
  W( 110 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 656 )
  W( 116 ) = W( 116 ) + a*JVS( 657 )
  W( 117 ) = W( 117 ) + a*JVS( 658 )
  W( 118 ) = W( 118 ) + a*JVS( 659 )
  W( 119 ) = W( 119 ) + a*JVS( 660 )
  W( 133 ) = W( 133 ) + a*JVS( 661 )
  W( 134 ) = W( 134 ) + a*JVS( 662 )
  W( 135 ) = W( 135 ) + a*JVS( 663 )
  W( 136 ) = W( 136 ) + a*JVS( 664 )
  W( 139 ) = W( 139 ) + a*JVS( 665 )
  W( 140 ) = W( 140 ) + a*JVS( 666 )
  a = -W( 115 ) / JVS(          763  )
  W( 115 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 764 )
  W( 124 ) = W( 124 ) + a*JVS( 765 )
  W( 133 ) = W( 133 ) + a*JVS( 766 )
  W( 134 ) = W( 134 ) + a*JVS( 767 )
  W( 135 ) = W( 135 ) + a*JVS( 768 )
  W( 136 ) = W( 136 ) + a*JVS( 769 )
  W( 139 ) = W( 139 ) + a*JVS( 770 )
  W( 140 ) = W( 140 ) + a*JVS( 771 )
  a = -W( 116 ) / JVS(          774  )
  W( 116 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 775 )
  W( 134 ) = W( 134 ) + a*JVS( 776 )
  W( 135 ) = W( 135 ) + a*JVS( 777 )
  W( 136 ) = W( 136 ) + a*JVS( 778 )
  W( 139 ) = W( 139 ) + a*JVS( 779 )
  W( 140 ) = W( 140 ) + a*JVS( 780 )
  a = -W( 117 ) / JVS(          788  )
  W( 117 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 789 )
  W( 134 ) = W( 134 ) + a*JVS( 790 )
  W( 135 ) = W( 135 ) + a*JVS( 791 )
  W( 136 ) = W( 136 ) + a*JVS( 792 )
  W( 139 ) = W( 139 ) + a*JVS( 793 )
  a = -W( 118 ) / JVS(          795  )
  W( 118 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 796 )
  W( 133 ) = W( 133 ) + a*JVS( 797 )
  W( 134 ) = W( 134 ) + a*JVS( 798 )
  W( 135 ) = W( 135 ) + a*JVS( 799 )
  W( 139 ) = W( 139 ) + a*JVS( 800 )
  W( 140 ) = W( 140 ) + a*JVS( 801 )
  a = -W( 119 ) / JVS(          808  )
  W( 119 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 809 )
  W( 134 ) = W( 134 ) + a*JVS( 810 )
  W( 135 ) = W( 135 ) + a*JVS( 811 )
  W( 136 ) = W( 136 ) + a*JVS( 812 )
  W( 139 ) = W( 139 ) + a*JVS( 813 )
  W( 140 ) = W( 140 ) + a*JVS( 814 )
  JVS( 830) = W( 50 )
  JVS( 831) = W( 86 )
  JVS( 832) = W( 87 )
  JVS( 833) = W( 104 )
  JVS( 834) = W( 107 )
  JVS( 835) = W( 109 )
  JVS( 836) = W( 110 )
  JVS( 837) = W( 115 )
  JVS( 838) = W( 116 )
  JVS( 839) = W( 117 )
  JVS( 840) = W( 118 )
  JVS( 841) = W( 119 )
  JVS( 842) = W( 121 )
  JVS( 843) = W( 124 )
  JVS( 844) = W( 133 )
  JVS( 845) = W( 134 )
  JVS( 846) = W( 135 )
  JVS( 847) = W( 136 )
  JVS( 848) = W( 139 )
  JVS( 849) = W( 140 )
  IF ( ABS(  JVS( 853 )) < TINY(a) ) THEN
         IER = 122                                     
         RETURN
  END IF
   W( 104 ) = JVS( 850 )
   W( 107 ) = JVS( 851 )
   W( 109 ) = JVS( 852 )
   W( 122 ) = JVS( 853 )
   W( 124 ) = JVS( 854 )
   W( 133 ) = JVS( 855 )
   W( 134 ) = JVS( 856 )
   W( 135 ) = JVS( 857 )
   W( 136 ) = JVS( 858 )
   W( 139 ) = JVS( 859 )
   W( 140 ) = JVS( 860 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  JVS( 850) = W( 104 )
  JVS( 851) = W( 107 )
  JVS( 852) = W( 109 )
  JVS( 853) = W( 122 )
  JVS( 854) = W( 124 )
  JVS( 855) = W( 133 )
  JVS( 856) = W( 134 )
  JVS( 857) = W( 135 )
  JVS( 858) = W( 136 )
  JVS( 859) = W( 139 )
  JVS( 860) = W( 140 )
  IF ( ABS(  JVS( 864 )) < TINY(a) ) THEN
         IER = 123                                     
         RETURN
  END IF
   W( 80 ) = JVS( 861 )
   W( 103 ) = JVS( 862 )
   W( 108 ) = JVS( 863 )
   W( 123 ) = JVS( 864 )
   W( 131 ) = JVS( 865 )
   W( 133 ) = JVS( 866 )
   W( 134 ) = JVS( 867 )
   W( 135 ) = JVS( 868 )
   W( 136 ) = JVS( 869 )
   W( 137 ) = JVS( 870 )
   W( 139 ) = JVS( 871 )
   W( 140 ) = JVS( 872 )
  a = -W( 80 ) / JVS(          421  )
  W( 80 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 422 )
  W( 134 ) = W( 134 ) + a*JVS( 423 )
  W( 135 ) = W( 135 ) + a*JVS( 424 )
  a = -W( 103 ) / JVS(          600  )
  W( 103 ) = -a
  W( 108 ) = W( 108 ) + a*JVS( 601 )
  W( 133 ) = W( 133 ) + a*JVS( 602 )
  W( 135 ) = W( 135 ) + a*JVS( 603 )
  W( 139 ) = W( 139 ) + a*JVS( 604 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  JVS( 861) = W( 80 )
  JVS( 862) = W( 103 )
  JVS( 863) = W( 108 )
  JVS( 864) = W( 123 )
  JVS( 865) = W( 131 )
  JVS( 866) = W( 133 )
  JVS( 867) = W( 134 )
  JVS( 868) = W( 135 )
  JVS( 869) = W( 136 )
  JVS( 870) = W( 137 )
  JVS( 871) = W( 139 )
  JVS( 872) = W( 140 )
  IF ( ABS(  JVS( 881 )) < TINY(a) ) THEN
         IER = 124                                     
         RETURN
  END IF
   W( 69 ) = JVS( 873 )
   W( 86 ) = JVS( 874 )
   W( 87 ) = JVS( 875 )
   W( 104 ) = JVS( 876 )
   W( 107 ) = JVS( 877 )
   W( 109 ) = JVS( 878 )
   W( 119 ) = JVS( 879 )
   W( 122 ) = JVS( 880 )
   W( 124 ) = JVS( 881 )
   W( 133 ) = JVS( 882 )
   W( 134 ) = JVS( 883 )
   W( 135 ) = JVS( 884 )
   W( 136 ) = JVS( 885 )
   W( 139 ) = JVS( 886 )
   W( 140 ) = JVS( 887 )
  a = -W( 69 ) / JVS(          359  )
  W( 69 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 360 )
  W( 134 ) = W( 134 ) + a*JVS( 361 )
  W( 135 ) = W( 135 ) + a*JVS( 362 )
  a = -W( 86 ) / JVS(          465  )
  W( 86 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 466 )
  W( 136 ) = W( 136 ) + a*JVS( 467 )
  a = -W( 87 ) / JVS(          468  )
  W( 87 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 469 )
  W( 136 ) = W( 136 ) + a*JVS( 470 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  a = -W( 119 ) / JVS(          808  )
  W( 119 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 809 )
  W( 134 ) = W( 134 ) + a*JVS( 810 )
  W( 135 ) = W( 135 ) + a*JVS( 811 )
  W( 136 ) = W( 136 ) + a*JVS( 812 )
  W( 139 ) = W( 139 ) + a*JVS( 813 )
  W( 140 ) = W( 140 ) + a*JVS( 814 )
  a = -W( 122 ) / JVS(          853  )
  W( 122 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 854 )
  W( 133 ) = W( 133 ) + a*JVS( 855 )
  W( 134 ) = W( 134 ) + a*JVS( 856 )
  W( 135 ) = W( 135 ) + a*JVS( 857 )
  W( 136 ) = W( 136 ) + a*JVS( 858 )
  W( 139 ) = W( 139 ) + a*JVS( 859 )
  W( 140 ) = W( 140 ) + a*JVS( 860 )
  JVS( 873) = W( 69 )
  JVS( 874) = W( 86 )
  JVS( 875) = W( 87 )
  JVS( 876) = W( 104 )
  JVS( 877) = W( 107 )
  JVS( 878) = W( 109 )
  JVS( 879) = W( 119 )
  JVS( 880) = W( 122 )
  JVS( 881) = W( 124 )
  JVS( 882) = W( 133 )
  JVS( 883) = W( 134 )
  JVS( 884) = W( 135 )
  JVS( 885) = W( 136 )
  JVS( 886) = W( 139 )
  JVS( 887) = W( 140 )
  IF ( ABS(  JVS( 894 )) < TINY(a) ) THEN
         IER = 125                                     
         RETURN
  END IF
   W( 72 ) = JVS( 888 )
   W( 81 ) = JVS( 889 )
   W( 90 ) = JVS( 890 )
   W( 113 ) = JVS( 891 )
   W( 121 ) = JVS( 892 )
   W( 124 ) = JVS( 893 )
   W( 125 ) = JVS( 894 )
   W( 127 ) = JVS( 895 )
   W( 128 ) = JVS( 896 )
   W( 129 ) = JVS( 897 )
   W( 131 ) = JVS( 898 )
   W( 132 ) = JVS( 899 )
   W( 133 ) = JVS( 900 )
   W( 134 ) = JVS( 901 )
   W( 135 ) = JVS( 902 )
   W( 136 ) = JVS( 903 )
   W( 137 ) = JVS( 904 )
   W( 138 ) = JVS( 905 )
   W( 139 ) = JVS( 906 )
   W( 140 ) = JVS( 907 )
  a = -W( 72 ) / JVS(          376  )
  W( 72 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 377 )
  W( 134 ) = W( 134 ) + a*JVS( 378 )
  W( 135 ) = W( 135 ) + a*JVS( 379 )
  a = -W( 81 ) / JVS(          425  )
  W( 81 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 426 )
  W( 134 ) = W( 134 ) + a*JVS( 427 )
  W( 138 ) = W( 138 ) + a*JVS( 428 )
  a = -W( 90 ) / JVS(          484  )
  W( 90 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 485 )
  W( 134 ) = W( 134 ) + a*JVS( 486 )
  W( 135 ) = W( 135 ) + a*JVS( 487 )
  W( 139 ) = W( 139 ) + a*JVS( 488 )
  a = -W( 113 ) / JVS(          732  )
  W( 113 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 733 )
  W( 133 ) = W( 133 ) + a*JVS( 734 )
  W( 134 ) = W( 134 ) + a*JVS( 735 )
  W( 136 ) = W( 136 ) + a*JVS( 736 )
  a = -W( 121 ) / JVS(          842  )
  W( 121 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 843 )
  W( 133 ) = W( 133 ) + a*JVS( 844 )
  W( 134 ) = W( 134 ) + a*JVS( 845 )
  W( 135 ) = W( 135 ) + a*JVS( 846 )
  W( 136 ) = W( 136 ) + a*JVS( 847 )
  W( 139 ) = W( 139 ) + a*JVS( 848 )
  W( 140 ) = W( 140 ) + a*JVS( 849 )
  a = -W( 124 ) / JVS(          881  )
  W( 124 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 882 )
  W( 134 ) = W( 134 ) + a*JVS( 883 )
  W( 135 ) = W( 135 ) + a*JVS( 884 )
  W( 136 ) = W( 136 ) + a*JVS( 885 )
  W( 139 ) = W( 139 ) + a*JVS( 886 )
  W( 140 ) = W( 140 ) + a*JVS( 887 )
  JVS( 888) = W( 72 )
  JVS( 889) = W( 81 )
  JVS( 890) = W( 90 )
  JVS( 891) = W( 113 )
  JVS( 892) = W( 121 )
  JVS( 893) = W( 124 )
  JVS( 894) = W( 125 )
  JVS( 895) = W( 127 )
  JVS( 896) = W( 128 )
  JVS( 897) = W( 129 )
  JVS( 898) = W( 131 )
  JVS( 899) = W( 132 )
  JVS( 900) = W( 133 )
  JVS( 901) = W( 134 )
  JVS( 902) = W( 135 )
  JVS( 903) = W( 136 )
  JVS( 904) = W( 137 )
  JVS( 905) = W( 138 )
  JVS( 906) = W( 139 )
  JVS( 907) = W( 140 )
  IF ( ABS(  JVS( 927 )) < TINY(a) ) THEN
         IER = 126                                     
         RETURN
  END IF
   W( 44 ) = JVS( 908 )
   W( 61 ) = JVS( 909 )
   W( 71 ) = JVS( 910 )
   W( 74 ) = JVS( 911 )
   W( 93 ) = JVS( 912 )
   W( 94 ) = JVS( 913 )
   W( 95 ) = JVS( 914 )
   W( 98 ) = JVS( 915 )
   W( 101 ) = JVS( 916 )
   W( 104 ) = JVS( 917 )
   W( 105 ) = JVS( 918 )
   W( 107 ) = JVS( 919 )
   W( 109 ) = JVS( 920 )
   W( 113 ) = JVS( 921 )
   W( 117 ) = JVS( 922 )
   W( 121 ) = JVS( 923 )
   W( 123 ) = JVS( 924 )
   W( 124 ) = JVS( 925 )
   W( 125 ) = JVS( 926 )
   W( 126 ) = JVS( 927 )
   W( 127 ) = JVS( 928 )
   W( 128 ) = JVS( 929 )
   W( 129 ) = JVS( 930 )
   W( 131 ) = JVS( 931 )
   W( 132 ) = JVS( 932 )
   W( 133 ) = JVS( 933 )
   W( 134 ) = JVS( 934 )
   W( 135 ) = JVS( 935 )
   W( 136 ) = JVS( 936 )
   W( 137 ) = JVS( 937 )
   W( 138 ) = JVS( 938 )
   W( 139 ) = JVS( 939 )
   W( 140 ) = JVS( 940 )
  a = -W( 44 ) / JVS(          230  )
  W( 44 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 231 )
  W( 139 ) = W( 139 ) + a*JVS( 232 )
  a = -W( 61 ) / JVS(          309  )
  W( 61 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 310 )
  W( 134 ) = W( 134 ) + a*JVS( 311 )
  W( 135 ) = W( 135 ) + a*JVS( 312 )
  a = -W( 71 ) / JVS(          369  )
  W( 71 ) = -a
  W( 105 ) = W( 105 ) + a*JVS( 370 )
  W( 113 ) = W( 113 ) + a*JVS( 371 )
  W( 117 ) = W( 117 ) + a*JVS( 372 )
  W( 133 ) = W( 133 ) + a*JVS( 373 )
  W( 134 ) = W( 134 ) + a*JVS( 374 )
  W( 139 ) = W( 139 ) + a*JVS( 375 )
  a = -W( 74 ) / JVS(          386  )
  W( 74 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 387 )
  W( 135 ) = W( 135 ) + a*JVS( 388 )
  W( 139 ) = W( 139 ) + a*JVS( 389 )
  a = -W( 93 ) / JVS(          504  )
  W( 93 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 505 )
  W( 135 ) = W( 135 ) + a*JVS( 506 )
  W( 139 ) = W( 139 ) + a*JVS( 507 )
  a = -W( 94 ) / JVS(          510  )
  W( 94 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 511 )
  W( 135 ) = W( 135 ) + a*JVS( 512 )
  W( 139 ) = W( 139 ) + a*JVS( 513 )
  a = -W( 95 ) / JVS(          517  )
  W( 95 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 518 )
  W( 135 ) = W( 135 ) + a*JVS( 519 )
  W( 138 ) = W( 138 ) + a*JVS( 520 )
  W( 139 ) = W( 139 ) + a*JVS( 521 )
  a = -W( 98 ) / JVS(          540  )
  W( 98 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 541 )
  W( 135 ) = W( 135 ) + a*JVS( 542 )
  W( 138 ) = W( 138 ) + a*JVS( 543 )
  W( 139 ) = W( 139 ) + a*JVS( 544 )
  a = -W( 101 ) / JVS(          582  )
  W( 101 ) = -a
  W( 104 ) = W( 104 ) + a*JVS( 583 )
  W( 107 ) = W( 107 ) + a*JVS( 584 )
  W( 109 ) = W( 109 ) + a*JVS( 585 )
  W( 134 ) = W( 134 ) + a*JVS( 586 )
  W( 135 ) = W( 135 ) + a*JVS( 587 )
  W( 136 ) = W( 136 ) + a*JVS( 588 )
  W( 139 ) = W( 139 ) + a*JVS( 589 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 105 ) / JVS(          610  )
  W( 105 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 611 )
  W( 134 ) = W( 134 ) + a*JVS( 612 )
  W( 135 ) = W( 135 ) + a*JVS( 613 )
  W( 136 ) = W( 136 ) + a*JVS( 614 )
  W( 139 ) = W( 139 ) + a*JVS( 615 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  a = -W( 113 ) / JVS(          732  )
  W( 113 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 733 )
  W( 133 ) = W( 133 ) + a*JVS( 734 )
  W( 134 ) = W( 134 ) + a*JVS( 735 )
  W( 136 ) = W( 136 ) + a*JVS( 736 )
  a = -W( 117 ) / JVS(          788  )
  W( 117 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 789 )
  W( 134 ) = W( 134 ) + a*JVS( 790 )
  W( 135 ) = W( 135 ) + a*JVS( 791 )
  W( 136 ) = W( 136 ) + a*JVS( 792 )
  W( 139 ) = W( 139 ) + a*JVS( 793 )
  a = -W( 121 ) / JVS(          842  )
  W( 121 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 843 )
  W( 133 ) = W( 133 ) + a*JVS( 844 )
  W( 134 ) = W( 134 ) + a*JVS( 845 )
  W( 135 ) = W( 135 ) + a*JVS( 846 )
  W( 136 ) = W( 136 ) + a*JVS( 847 )
  W( 139 ) = W( 139 ) + a*JVS( 848 )
  W( 140 ) = W( 140 ) + a*JVS( 849 )
  a = -W( 123 ) / JVS(          864  )
  W( 123 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 865 )
  W( 133 ) = W( 133 ) + a*JVS( 866 )
  W( 134 ) = W( 134 ) + a*JVS( 867 )
  W( 135 ) = W( 135 ) + a*JVS( 868 )
  W( 136 ) = W( 136 ) + a*JVS( 869 )
  W( 137 ) = W( 137 ) + a*JVS( 870 )
  W( 139 ) = W( 139 ) + a*JVS( 871 )
  W( 140 ) = W( 140 ) + a*JVS( 872 )
  a = -W( 124 ) / JVS(          881  )
  W( 124 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 882 )
  W( 134 ) = W( 134 ) + a*JVS( 883 )
  W( 135 ) = W( 135 ) + a*JVS( 884 )
  W( 136 ) = W( 136 ) + a*JVS( 885 )
  W( 139 ) = W( 139 ) + a*JVS( 886 )
  W( 140 ) = W( 140 ) + a*JVS( 887 )
  a = -W( 125 ) / JVS(          894  )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 895 )
  W( 128 ) = W( 128 ) + a*JVS( 896 )
  W( 129 ) = W( 129 ) + a*JVS( 897 )
  W( 131 ) = W( 131 ) + a*JVS( 898 )
  W( 132 ) = W( 132 ) + a*JVS( 899 )
  W( 133 ) = W( 133 ) + a*JVS( 900 )
  W( 134 ) = W( 134 ) + a*JVS( 901 )
  W( 135 ) = W( 135 ) + a*JVS( 902 )
  W( 136 ) = W( 136 ) + a*JVS( 903 )
  W( 137 ) = W( 137 ) + a*JVS( 904 )
  W( 138 ) = W( 138 ) + a*JVS( 905 )
  W( 139 ) = W( 139 ) + a*JVS( 906 )
  W( 140 ) = W( 140 ) + a*JVS( 907 )
  JVS( 908) = W( 44 )
  JVS( 909) = W( 61 )
  JVS( 910) = W( 71 )
  JVS( 911) = W( 74 )
  JVS( 912) = W( 93 )
  JVS( 913) = W( 94 )
  JVS( 914) = W( 95 )
  JVS( 915) = W( 98 )
  JVS( 916) = W( 101 )
  JVS( 917) = W( 104 )
  JVS( 918) = W( 105 )
  JVS( 919) = W( 107 )
  JVS( 920) = W( 109 )
  JVS( 921) = W( 113 )
  JVS( 922) = W( 117 )
  JVS( 923) = W( 121 )
  JVS( 924) = W( 123 )
  JVS( 925) = W( 124 )
  JVS( 926) = W( 125 )
  JVS( 927) = W( 126 )
  JVS( 928) = W( 127 )
  JVS( 929) = W( 128 )
  JVS( 930) = W( 129 )
  JVS( 931) = W( 131 )
  JVS( 932) = W( 132 )
  JVS( 933) = W( 133 )
  JVS( 934) = W( 134 )
  JVS( 935) = W( 135 )
  JVS( 936) = W( 136 )
  JVS( 937) = W( 137 )
  JVS( 938) = W( 138 )
  JVS( 939) = W( 139 )
  JVS( 940) = W( 140 )
  IF ( ABS(  JVS( 947 )) < TINY(a) ) THEN
         IER = 127                                     
         RETURN
  END IF
   W( 29 ) = JVS( 941 )
   W( 57 ) = JVS( 942 )
   W( 80 ) = JVS( 943 )
   W( 120 ) = JVS( 944 )
   W( 122 ) = JVS( 945 )
   W( 124 ) = JVS( 946 )
   W( 127 ) = JVS( 947 )
   W( 129 ) = JVS( 948 )
   W( 131 ) = JVS( 949 )
   W( 133 ) = JVS( 950 )
   W( 134 ) = JVS( 951 )
   W( 135 ) = JVS( 952 )
   W( 136 ) = JVS( 953 )
   W( 137 ) = JVS( 954 )
   W( 139 ) = JVS( 955 )
   W( 140 ) = JVS( 956 )
  a = -W( 29 ) / JVS(          176  )
  W( 29 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 177 )
  W( 134 ) = W( 134 ) + a*JVS( 178 )
  W( 135 ) = W( 135 ) + a*JVS( 179 )
  a = -W( 57 ) / JVS(          285  )
  W( 57 ) = -a
  W( 120 ) = W( 120 ) + a*JVS( 286 )
  W( 131 ) = W( 131 ) + a*JVS( 287 )
  W( 133 ) = W( 133 ) + a*JVS( 288 )
  W( 134 ) = W( 134 ) + a*JVS( 289 )
  W( 137 ) = W( 137 ) + a*JVS( 290 )
  W( 139 ) = W( 139 ) + a*JVS( 291 )
  W( 140 ) = W( 140 ) + a*JVS( 292 )
  a = -W( 80 ) / JVS(          421  )
  W( 80 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 422 )
  W( 134 ) = W( 134 ) + a*JVS( 423 )
  W( 135 ) = W( 135 ) + a*JVS( 424 )
  a = -W( 120 ) / JVS(          819  )
  W( 120 ) = -a
  W( 122 ) = W( 122 ) + a*JVS( 820 )
  W( 124 ) = W( 124 ) + a*JVS( 821 )
  W( 129 ) = W( 129 ) + a*JVS( 822 )
  W( 131 ) = W( 131 ) + a*JVS( 823 )
  W( 133 ) = W( 133 ) + a*JVS( 824 )
  W( 134 ) = W( 134 ) + a*JVS( 825 )
  W( 135 ) = W( 135 ) + a*JVS( 826 )
  W( 136 ) = W( 136 ) + a*JVS( 827 )
  W( 139 ) = W( 139 ) + a*JVS( 828 )
  W( 140 ) = W( 140 ) + a*JVS( 829 )
  a = -W( 122 ) / JVS(          853  )
  W( 122 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 854 )
  W( 133 ) = W( 133 ) + a*JVS( 855 )
  W( 134 ) = W( 134 ) + a*JVS( 856 )
  W( 135 ) = W( 135 ) + a*JVS( 857 )
  W( 136 ) = W( 136 ) + a*JVS( 858 )
  W( 139 ) = W( 139 ) + a*JVS( 859 )
  W( 140 ) = W( 140 ) + a*JVS( 860 )
  a = -W( 124 ) / JVS(          881  )
  W( 124 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 882 )
  W( 134 ) = W( 134 ) + a*JVS( 883 )
  W( 135 ) = W( 135 ) + a*JVS( 884 )
  W( 136 ) = W( 136 ) + a*JVS( 885 )
  W( 139 ) = W( 139 ) + a*JVS( 886 )
  W( 140 ) = W( 140 ) + a*JVS( 887 )
  JVS( 941) = W( 29 )
  JVS( 942) = W( 57 )
  JVS( 943) = W( 80 )
  JVS( 944) = W( 120 )
  JVS( 945) = W( 122 )
  JVS( 946) = W( 124 )
  JVS( 947) = W( 127 )
  JVS( 948) = W( 129 )
  JVS( 949) = W( 131 )
  JVS( 950) = W( 133 )
  JVS( 951) = W( 134 )
  JVS( 952) = W( 135 )
  JVS( 953) = W( 136 )
  JVS( 954) = W( 137 )
  JVS( 955) = W( 139 )
  JVS( 956) = W( 140 )
  IF ( ABS(  JVS( 960 )) < TINY(a) ) THEN
         IER = 128                                     
         RETURN
  END IF
   W( 80 ) = JVS( 957 )
   W( 103 ) = JVS( 958 )
   W( 108 ) = JVS( 959 )
   W( 128 ) = JVS( 960 )
   W( 131 ) = JVS( 961 )
   W( 133 ) = JVS( 962 )
   W( 134 ) = JVS( 963 )
   W( 135 ) = JVS( 964 )
   W( 136 ) = JVS( 965 )
   W( 137 ) = JVS( 966 )
   W( 139 ) = JVS( 967 )
   W( 140 ) = JVS( 968 )
  a = -W( 80 ) / JVS(          421  )
  W( 80 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 422 )
  W( 134 ) = W( 134 ) + a*JVS( 423 )
  W( 135 ) = W( 135 ) + a*JVS( 424 )
  a = -W( 103 ) / JVS(          600  )
  W( 103 ) = -a
  W( 108 ) = W( 108 ) + a*JVS( 601 )
  W( 133 ) = W( 133 ) + a*JVS( 602 )
  W( 135 ) = W( 135 ) + a*JVS( 603 )
  W( 139 ) = W( 139 ) + a*JVS( 604 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  JVS( 957) = W( 80 )
  JVS( 958) = W( 103 )
  JVS( 959) = W( 108 )
  JVS( 960) = W( 128 )
  JVS( 961) = W( 131 )
  JVS( 962) = W( 133 )
  JVS( 963) = W( 134 )
  JVS( 964) = W( 135 )
  JVS( 965) = W( 136 )
  JVS( 966) = W( 137 )
  JVS( 967) = W( 139 )
  JVS( 968) = W( 140 )
  IF ( ABS(  JVS( 972 )) < TINY(a) ) THEN
         IER = 129                                     
         RETURN
  END IF
   W( 54 ) = JVS( 969 )
   W( 123 ) = JVS( 970 )
   W( 128 ) = JVS( 971 )
   W( 129 ) = JVS( 972 )
   W( 131 ) = JVS( 973 )
   W( 133 ) = JVS( 974 )
   W( 134 ) = JVS( 975 )
   W( 135 ) = JVS( 976 )
   W( 136 ) = JVS( 977 )
   W( 137 ) = JVS( 978 )
   W( 139 ) = JVS( 979 )
   W( 140 ) = JVS( 980 )
  a = -W( 54 ) / JVS(          272  )
  W( 54 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 273 )
  W( 134 ) = W( 134 ) + a*JVS( 274 )
  W( 135 ) = W( 135 ) + a*JVS( 275 )
  a = -W( 123 ) / JVS(          864  )
  W( 123 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 865 )
  W( 133 ) = W( 133 ) + a*JVS( 866 )
  W( 134 ) = W( 134 ) + a*JVS( 867 )
  W( 135 ) = W( 135 ) + a*JVS( 868 )
  W( 136 ) = W( 136 ) + a*JVS( 869 )
  W( 137 ) = W( 137 ) + a*JVS( 870 )
  W( 139 ) = W( 139 ) + a*JVS( 871 )
  W( 140 ) = W( 140 ) + a*JVS( 872 )
  a = -W( 128 ) / JVS(          960  )
  W( 128 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 961 )
  W( 133 ) = W( 133 ) + a*JVS( 962 )
  W( 134 ) = W( 134 ) + a*JVS( 963 )
  W( 135 ) = W( 135 ) + a*JVS( 964 )
  W( 136 ) = W( 136 ) + a*JVS( 965 )
  W( 137 ) = W( 137 ) + a*JVS( 966 )
  W( 139 ) = W( 139 ) + a*JVS( 967 )
  W( 140 ) = W( 140 ) + a*JVS( 968 )
  JVS( 969) = W( 54 )
  JVS( 970) = W( 123 )
  JVS( 971) = W( 128 )
  JVS( 972) = W( 129 )
  JVS( 973) = W( 131 )
  JVS( 974) = W( 133 )
  JVS( 975) = W( 134 )
  JVS( 976) = W( 135 )
  JVS( 977) = W( 136 )
  JVS( 978) = W( 137 )
  JVS( 979) = W( 139 )
  JVS( 980) = W( 140 )
  IF ( ABS(  JVS( 1028 )) < TINY(a) ) THEN
         IER = 130                                     
         RETURN
  END IF
   W( 30 ) = JVS( 981 )
   W( 50 ) = JVS( 982 )
   W( 51 ) = JVS( 983 )
   W( 53 ) = JVS( 984 )
   W( 59 ) = JVS( 985 )
   W( 60 ) = JVS( 986 )
   W( 65 ) = JVS( 987 )
   W( 68 ) = JVS( 988 )
   W( 69 ) = JVS( 989 )
   W( 72 ) = JVS( 990 )
   W( 73 ) = JVS( 991 )
   W( 77 ) = JVS( 992 )
   W( 79 ) = JVS( 993 )
   W( 80 ) = JVS( 994 )
   W( 81 ) = JVS( 995 )
   W( 83 ) = JVS( 996 )
   W( 84 ) = JVS( 997 )
   W( 86 ) = JVS( 998 )
   W( 87 ) = JVS( 999 )
   W( 88 ) = JVS( 1000 )
   W( 90 ) = JVS( 1001 )
   W( 92 ) = JVS( 1002 )
   W( 97 ) = JVS( 1003 )
   W( 100 ) = JVS( 1004 )
   W( 102 ) = JVS( 1005 )
   W( 103 ) = JVS( 1006 )
   W( 104 ) = JVS( 1007 )
   W( 105 ) = JVS( 1008 )
   W( 107 ) = JVS( 1009 )
   W( 108 ) = JVS( 1010 )
   W( 109 ) = JVS( 1011 )
   W( 112 ) = JVS( 1012 )
   W( 113 ) = JVS( 1013 )
   W( 115 ) = JVS( 1014 )
   W( 116 ) = JVS( 1015 )
   W( 117 ) = JVS( 1016 )
   W( 118 ) = JVS( 1017 )
   W( 119 ) = JVS( 1018 )
   W( 120 ) = JVS( 1019 )
   W( 121 ) = JVS( 1020 )
   W( 122 ) = JVS( 1021 )
   W( 123 ) = JVS( 1022 )
   W( 124 ) = JVS( 1023 )
   W( 125 ) = JVS( 1024 )
   W( 127 ) = JVS( 1025 )
   W( 128 ) = JVS( 1026 )
   W( 129 ) = JVS( 1027 )
   W( 130 ) = JVS( 1028 )
   W( 131 ) = JVS( 1029 )
   W( 132 ) = JVS( 1030 )
   W( 133 ) = JVS( 1031 )
   W( 134 ) = JVS( 1032 )
   W( 135 ) = JVS( 1033 )
   W( 136 ) = JVS( 1034 )
   W( 137 ) = JVS( 1035 )
   W( 138 ) = JVS( 1036 )
   W( 139 ) = JVS( 1037 )
   W( 140 ) = JVS( 1038 )
  a = -W( 30 ) / JVS(          180  )
  W( 30 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 181 )
  W( 139 ) = W( 139 ) + a*JVS( 182 )
  a = -W( 50 ) / JVS(          256  )
  W( 50 ) = -a
  W( 121 ) = W( 121 ) + a*JVS( 257 )
  W( 134 ) = W( 134 ) + a*JVS( 258 )
  W( 135 ) = W( 135 ) + a*JVS( 259 )
  a = -W( 51 ) / JVS(          261  )
  W( 51 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 262 )
  W( 139 ) = W( 139 ) + a*JVS( 263 )
  a = -W( 53 ) / JVS(          268  )
  W( 53 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 269 )
  W( 135 ) = W( 135 ) + a*JVS( 270 )
  W( 140 ) = W( 140 ) + a*JVS( 271 )
  a = -W( 59 ) / JVS(          301  )
  W( 59 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 302 )
  W( 136 ) = W( 136 ) + a*JVS( 303 )
  a = -W( 60 ) / JVS(          305  )
  W( 60 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 306 )
  W( 136 ) = W( 136 ) + a*JVS( 307 )
  W( 139 ) = W( 139 ) + a*JVS( 308 )
  a = -W( 65 ) / JVS(          336  )
  W( 65 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 337 )
  W( 135 ) = W( 135 ) + a*JVS( 338 )
  W( 139 ) = W( 139 ) + a*JVS( 339 )
  a = -W( 68 ) / JVS(          354  )
  W( 68 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 355 )
  W( 134 ) = W( 134 ) + a*JVS( 356 )
  W( 135 ) = W( 135 ) + a*JVS( 357 )
  W( 137 ) = W( 137 ) + a*JVS( 358 )
  a = -W( 69 ) / JVS(          359  )
  W( 69 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 360 )
  W( 134 ) = W( 134 ) + a*JVS( 361 )
  W( 135 ) = W( 135 ) + a*JVS( 362 )
  a = -W( 72 ) / JVS(          376  )
  W( 72 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 377 )
  W( 134 ) = W( 134 ) + a*JVS( 378 )
  W( 135 ) = W( 135 ) + a*JVS( 379 )
  a = -W( 73 ) / JVS(          380  )
  W( 73 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 381 )
  W( 137 ) = W( 137 ) + a*JVS( 382 )
  W( 138 ) = W( 138 ) + a*JVS( 383 )
  a = -W( 77 ) / JVS(          404  )
  W( 77 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 405 )
  W( 113 ) = W( 113 ) + a*JVS( 406 )
  W( 134 ) = W( 134 ) + a*JVS( 407 )
  W( 136 ) = W( 136 ) + a*JVS( 408 )
  a = -W( 79 ) / JVS(          417  )
  W( 79 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 418 )
  W( 134 ) = W( 134 ) + a*JVS( 419 )
  W( 135 ) = W( 135 ) + a*JVS( 420 )
  a = -W( 80 ) / JVS(          421  )
  W( 80 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 422 )
  W( 134 ) = W( 134 ) + a*JVS( 423 )
  W( 135 ) = W( 135 ) + a*JVS( 424 )
  a = -W( 81 ) / JVS(          425  )
  W( 81 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 426 )
  W( 134 ) = W( 134 ) + a*JVS( 427 )
  W( 138 ) = W( 138 ) + a*JVS( 428 )
  a = -W( 83 ) / JVS(          444  )
  W( 83 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 445 )
  W( 134 ) = W( 134 ) + a*JVS( 446 )
  W( 135 ) = W( 135 ) + a*JVS( 447 )
  a = -W( 84 ) / JVS(          451  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 452 )
  W( 130 ) = W( 130 ) + a*JVS( 453 )
  W( 134 ) = W( 134 ) + a*JVS( 454 )
  W( 136 ) = W( 136 ) + a*JVS( 455 )
  a = -W( 86 ) / JVS(          465  )
  W( 86 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 466 )
  W( 136 ) = W( 136 ) + a*JVS( 467 )
  a = -W( 87 ) / JVS(          468  )
  W( 87 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 469 )
  W( 136 ) = W( 136 ) + a*JVS( 470 )
  a = -W( 88 ) / JVS(          471  )
  W( 88 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 472 )
  W( 134 ) = W( 134 ) + a*JVS( 473 )
  W( 135 ) = W( 135 ) + a*JVS( 474 )
  a = -W( 90 ) / JVS(          484  )
  W( 90 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 485 )
  W( 134 ) = W( 134 ) + a*JVS( 486 )
  W( 135 ) = W( 135 ) + a*JVS( 487 )
  W( 139 ) = W( 139 ) + a*JVS( 488 )
  a = -W( 92 ) / JVS(          498  )
  W( 92 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 499 )
  W( 134 ) = W( 134 ) + a*JVS( 500 )
  W( 136 ) = W( 136 ) + a*JVS( 501 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 533 )
  W( 135 ) = W( 135 ) + a*JVS( 534 )
  W( 139 ) = W( 139 ) + a*JVS( 535 )
  W( 140 ) = W( 140 ) + a*JVS( 536 )
  a = -W( 100 ) / JVS(          566  )
  W( 100 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 567 )
  W( 116 ) = W( 116 ) + a*JVS( 568 )
  W( 118 ) = W( 118 ) + a*JVS( 569 )
  W( 119 ) = W( 119 ) + a*JVS( 570 )
  W( 121 ) = W( 121 ) + a*JVS( 571 )
  W( 122 ) = W( 122 ) + a*JVS( 572 )
  W( 127 ) = W( 127 ) + a*JVS( 573 )
  W( 129 ) = W( 129 ) + a*JVS( 574 )
  W( 131 ) = W( 131 ) + a*JVS( 575 )
  W( 134 ) = W( 134 ) + a*JVS( 576 )
  W( 140 ) = W( 140 ) + a*JVS( 577 )
  a = -W( 102 ) / JVS(          594  )
  W( 102 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 595 )
  W( 134 ) = W( 134 ) + a*JVS( 596 )
  W( 135 ) = W( 135 ) + a*JVS( 597 )
  W( 139 ) = W( 139 ) + a*JVS( 598 )
  W( 140 ) = W( 140 ) + a*JVS( 599 )
  a = -W( 103 ) / JVS(          600  )
  W( 103 ) = -a
  W( 108 ) = W( 108 ) + a*JVS( 601 )
  W( 133 ) = W( 133 ) + a*JVS( 602 )
  W( 135 ) = W( 135 ) + a*JVS( 603 )
  W( 139 ) = W( 139 ) + a*JVS( 604 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 105 ) / JVS(          610  )
  W( 105 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 611 )
  W( 134 ) = W( 134 ) + a*JVS( 612 )
  W( 135 ) = W( 135 ) + a*JVS( 613 )
  W( 136 ) = W( 136 ) + a*JVS( 614 )
  W( 139 ) = W( 139 ) + a*JVS( 615 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  a = -W( 112 ) / JVS(          718  )
  W( 112 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 719 )
  W( 118 ) = W( 118 ) + a*JVS( 720 )
  W( 127 ) = W( 127 ) + a*JVS( 721 )
  W( 129 ) = W( 129 ) + a*JVS( 722 )
  W( 131 ) = W( 131 ) + a*JVS( 723 )
  W( 133 ) = W( 133 ) + a*JVS( 724 )
  W( 134 ) = W( 134 ) + a*JVS( 725 )
  W( 135 ) = W( 135 ) + a*JVS( 726 )
  W( 136 ) = W( 136 ) + a*JVS( 727 )
  W( 137 ) = W( 137 ) + a*JVS( 728 )
  W( 139 ) = W( 139 ) + a*JVS( 729 )
  W( 140 ) = W( 140 ) + a*JVS( 730 )
  a = -W( 113 ) / JVS(          732  )
  W( 113 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 733 )
  W( 133 ) = W( 133 ) + a*JVS( 734 )
  W( 134 ) = W( 134 ) + a*JVS( 735 )
  W( 136 ) = W( 136 ) + a*JVS( 736 )
  a = -W( 115 ) / JVS(          763  )
  W( 115 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 764 )
  W( 124 ) = W( 124 ) + a*JVS( 765 )
  W( 133 ) = W( 133 ) + a*JVS( 766 )
  W( 134 ) = W( 134 ) + a*JVS( 767 )
  W( 135 ) = W( 135 ) + a*JVS( 768 )
  W( 136 ) = W( 136 ) + a*JVS( 769 )
  W( 139 ) = W( 139 ) + a*JVS( 770 )
  W( 140 ) = W( 140 ) + a*JVS( 771 )
  a = -W( 116 ) / JVS(          774  )
  W( 116 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 775 )
  W( 134 ) = W( 134 ) + a*JVS( 776 )
  W( 135 ) = W( 135 ) + a*JVS( 777 )
  W( 136 ) = W( 136 ) + a*JVS( 778 )
  W( 139 ) = W( 139 ) + a*JVS( 779 )
  W( 140 ) = W( 140 ) + a*JVS( 780 )
  a = -W( 117 ) / JVS(          788  )
  W( 117 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 789 )
  W( 134 ) = W( 134 ) + a*JVS( 790 )
  W( 135 ) = W( 135 ) + a*JVS( 791 )
  W( 136 ) = W( 136 ) + a*JVS( 792 )
  W( 139 ) = W( 139 ) + a*JVS( 793 )
  a = -W( 118 ) / JVS(          795  )
  W( 118 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 796 )
  W( 133 ) = W( 133 ) + a*JVS( 797 )
  W( 134 ) = W( 134 ) + a*JVS( 798 )
  W( 135 ) = W( 135 ) + a*JVS( 799 )
  W( 139 ) = W( 139 ) + a*JVS( 800 )
  W( 140 ) = W( 140 ) + a*JVS( 801 )
  a = -W( 119 ) / JVS(          808  )
  W( 119 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 809 )
  W( 134 ) = W( 134 ) + a*JVS( 810 )
  W( 135 ) = W( 135 ) + a*JVS( 811 )
  W( 136 ) = W( 136 ) + a*JVS( 812 )
  W( 139 ) = W( 139 ) + a*JVS( 813 )
  W( 140 ) = W( 140 ) + a*JVS( 814 )
  a = -W( 120 ) / JVS(          819  )
  W( 120 ) = -a
  W( 122 ) = W( 122 ) + a*JVS( 820 )
  W( 124 ) = W( 124 ) + a*JVS( 821 )
  W( 129 ) = W( 129 ) + a*JVS( 822 )
  W( 131 ) = W( 131 ) + a*JVS( 823 )
  W( 133 ) = W( 133 ) + a*JVS( 824 )
  W( 134 ) = W( 134 ) + a*JVS( 825 )
  W( 135 ) = W( 135 ) + a*JVS( 826 )
  W( 136 ) = W( 136 ) + a*JVS( 827 )
  W( 139 ) = W( 139 ) + a*JVS( 828 )
  W( 140 ) = W( 140 ) + a*JVS( 829 )
  a = -W( 121 ) / JVS(          842  )
  W( 121 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 843 )
  W( 133 ) = W( 133 ) + a*JVS( 844 )
  W( 134 ) = W( 134 ) + a*JVS( 845 )
  W( 135 ) = W( 135 ) + a*JVS( 846 )
  W( 136 ) = W( 136 ) + a*JVS( 847 )
  W( 139 ) = W( 139 ) + a*JVS( 848 )
  W( 140 ) = W( 140 ) + a*JVS( 849 )
  a = -W( 122 ) / JVS(          853  )
  W( 122 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 854 )
  W( 133 ) = W( 133 ) + a*JVS( 855 )
  W( 134 ) = W( 134 ) + a*JVS( 856 )
  W( 135 ) = W( 135 ) + a*JVS( 857 )
  W( 136 ) = W( 136 ) + a*JVS( 858 )
  W( 139 ) = W( 139 ) + a*JVS( 859 )
  W( 140 ) = W( 140 ) + a*JVS( 860 )
  a = -W( 123 ) / JVS(          864  )
  W( 123 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 865 )
  W( 133 ) = W( 133 ) + a*JVS( 866 )
  W( 134 ) = W( 134 ) + a*JVS( 867 )
  W( 135 ) = W( 135 ) + a*JVS( 868 )
  W( 136 ) = W( 136 ) + a*JVS( 869 )
  W( 137 ) = W( 137 ) + a*JVS( 870 )
  W( 139 ) = W( 139 ) + a*JVS( 871 )
  W( 140 ) = W( 140 ) + a*JVS( 872 )
  a = -W( 124 ) / JVS(          881  )
  W( 124 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 882 )
  W( 134 ) = W( 134 ) + a*JVS( 883 )
  W( 135 ) = W( 135 ) + a*JVS( 884 )
  W( 136 ) = W( 136 ) + a*JVS( 885 )
  W( 139 ) = W( 139 ) + a*JVS( 886 )
  W( 140 ) = W( 140 ) + a*JVS( 887 )
  a = -W( 125 ) / JVS(          894  )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 895 )
  W( 128 ) = W( 128 ) + a*JVS( 896 )
  W( 129 ) = W( 129 ) + a*JVS( 897 )
  W( 131 ) = W( 131 ) + a*JVS( 898 )
  W( 132 ) = W( 132 ) + a*JVS( 899 )
  W( 133 ) = W( 133 ) + a*JVS( 900 )
  W( 134 ) = W( 134 ) + a*JVS( 901 )
  W( 135 ) = W( 135 ) + a*JVS( 902 )
  W( 136 ) = W( 136 ) + a*JVS( 903 )
  W( 137 ) = W( 137 ) + a*JVS( 904 )
  W( 138 ) = W( 138 ) + a*JVS( 905 )
  W( 139 ) = W( 139 ) + a*JVS( 906 )
  W( 140 ) = W( 140 ) + a*JVS( 907 )
  a = -W( 127 ) / JVS(          947  )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 948 )
  W( 131 ) = W( 131 ) + a*JVS( 949 )
  W( 133 ) = W( 133 ) + a*JVS( 950 )
  W( 134 ) = W( 134 ) + a*JVS( 951 )
  W( 135 ) = W( 135 ) + a*JVS( 952 )
  W( 136 ) = W( 136 ) + a*JVS( 953 )
  W( 137 ) = W( 137 ) + a*JVS( 954 )
  W( 139 ) = W( 139 ) + a*JVS( 955 )
  W( 140 ) = W( 140 ) + a*JVS( 956 )
  a = -W( 128 ) / JVS(          960  )
  W( 128 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 961 )
  W( 133 ) = W( 133 ) + a*JVS( 962 )
  W( 134 ) = W( 134 ) + a*JVS( 963 )
  W( 135 ) = W( 135 ) + a*JVS( 964 )
  W( 136 ) = W( 136 ) + a*JVS( 965 )
  W( 137 ) = W( 137 ) + a*JVS( 966 )
  W( 139 ) = W( 139 ) + a*JVS( 967 )
  W( 140 ) = W( 140 ) + a*JVS( 968 )
  a = -W( 129 ) / JVS(          972  )
  W( 129 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 973 )
  W( 133 ) = W( 133 ) + a*JVS( 974 )
  W( 134 ) = W( 134 ) + a*JVS( 975 )
  W( 135 ) = W( 135 ) + a*JVS( 976 )
  W( 136 ) = W( 136 ) + a*JVS( 977 )
  W( 137 ) = W( 137 ) + a*JVS( 978 )
  W( 139 ) = W( 139 ) + a*JVS( 979 )
  W( 140 ) = W( 140 ) + a*JVS( 980 )
  JVS( 981) = W( 30 )
  JVS( 982) = W( 50 )
  JVS( 983) = W( 51 )
  JVS( 984) = W( 53 )
  JVS( 985) = W( 59 )
  JVS( 986) = W( 60 )
  JVS( 987) = W( 65 )
  JVS( 988) = W( 68 )
  JVS( 989) = W( 69 )
  JVS( 990) = W( 72 )
  JVS( 991) = W( 73 )
  JVS( 992) = W( 77 )
  JVS( 993) = W( 79 )
  JVS( 994) = W( 80 )
  JVS( 995) = W( 81 )
  JVS( 996) = W( 83 )
  JVS( 997) = W( 84 )
  JVS( 998) = W( 86 )
  JVS( 999) = W( 87 )
  JVS( 1000) = W( 88 )
  JVS( 1001) = W( 90 )
  JVS( 1002) = W( 92 )
  JVS( 1003) = W( 97 )
  JVS( 1004) = W( 100 )
  JVS( 1005) = W( 102 )
  JVS( 1006) = W( 103 )
  JVS( 1007) = W( 104 )
  JVS( 1008) = W( 105 )
  JVS( 1009) = W( 107 )
  JVS( 1010) = W( 108 )
  JVS( 1011) = W( 109 )
  JVS( 1012) = W( 112 )
  JVS( 1013) = W( 113 )
  JVS( 1014) = W( 115 )
  JVS( 1015) = W( 116 )
  JVS( 1016) = W( 117 )
  JVS( 1017) = W( 118 )
  JVS( 1018) = W( 119 )
  JVS( 1019) = W( 120 )
  JVS( 1020) = W( 121 )
  JVS( 1021) = W( 122 )
  JVS( 1022) = W( 123 )
  JVS( 1023) = W( 124 )
  JVS( 1024) = W( 125 )
  JVS( 1025) = W( 127 )
  JVS( 1026) = W( 128 )
  JVS( 1027) = W( 129 )
  JVS( 1028) = W( 130 )
  JVS( 1029) = W( 131 )
  JVS( 1030) = W( 132 )
  JVS( 1031) = W( 133 )
  JVS( 1032) = W( 134 )
  JVS( 1033) = W( 135 )
  JVS( 1034) = W( 136 )
  JVS( 1035) = W( 137 )
  JVS( 1036) = W( 138 )
  JVS( 1037) = W( 139 )
  JVS( 1038) = W( 140 )
  IF ( ABS(  JVS( 1041 )) < TINY(a) ) THEN
         IER = 131                                     
         RETURN
  END IF
   W( 80 ) = JVS( 1039 )
   W( 108 ) = JVS( 1040 )
   W( 131 ) = JVS( 1041 )
   W( 133 ) = JVS( 1042 )
   W( 134 ) = JVS( 1043 )
   W( 135 ) = JVS( 1044 )
   W( 136 ) = JVS( 1045 )
   W( 137 ) = JVS( 1046 )
   W( 139 ) = JVS( 1047 )
   W( 140 ) = JVS( 1048 )
  a = -W( 80 ) / JVS(          421  )
  W( 80 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 422 )
  W( 134 ) = W( 134 ) + a*JVS( 423 )
  W( 135 ) = W( 135 ) + a*JVS( 424 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  JVS( 1039) = W( 80 )
  JVS( 1040) = W( 108 )
  JVS( 1041) = W( 131 )
  JVS( 1042) = W( 133 )
  JVS( 1043) = W( 134 )
  JVS( 1044) = W( 135 )
  JVS( 1045) = W( 136 )
  JVS( 1046) = W( 137 )
  JVS( 1047) = W( 139 )
  JVS( 1048) = W( 140 )
  IF ( ABS(  JVS( 1055 )) < TINY(a) ) THEN
         IER = 132                                     
         RETURN
  END IF
   W( 54 ) = JVS( 1049 )
   W( 81 ) = JVS( 1050 )
   W( 108 ) = JVS( 1051 )
   W( 123 ) = JVS( 1052 )
   W( 129 ) = JVS( 1053 )
   W( 131 ) = JVS( 1054 )
   W( 132 ) = JVS( 1055 )
   W( 133 ) = JVS( 1056 )
   W( 134 ) = JVS( 1057 )
   W( 135 ) = JVS( 1058 )
   W( 136 ) = JVS( 1059 )
   W( 137 ) = JVS( 1060 )
   W( 138 ) = JVS( 1061 )
   W( 139 ) = JVS( 1062 )
   W( 140 ) = JVS( 1063 )
  a = -W( 54 ) / JVS(          272  )
  W( 54 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 273 )
  W( 134 ) = W( 134 ) + a*JVS( 274 )
  W( 135 ) = W( 135 ) + a*JVS( 275 )
  a = -W( 81 ) / JVS(          425  )
  W( 81 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 426 )
  W( 134 ) = W( 134 ) + a*JVS( 427 )
  W( 138 ) = W( 138 ) + a*JVS( 428 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  a = -W( 123 ) / JVS(          864  )
  W( 123 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 865 )
  W( 133 ) = W( 133 ) + a*JVS( 866 )
  W( 134 ) = W( 134 ) + a*JVS( 867 )
  W( 135 ) = W( 135 ) + a*JVS( 868 )
  W( 136 ) = W( 136 ) + a*JVS( 869 )
  W( 137 ) = W( 137 ) + a*JVS( 870 )
  W( 139 ) = W( 139 ) + a*JVS( 871 )
  W( 140 ) = W( 140 ) + a*JVS( 872 )
  a = -W( 129 ) / JVS(          972  )
  W( 129 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 973 )
  W( 133 ) = W( 133 ) + a*JVS( 974 )
  W( 134 ) = W( 134 ) + a*JVS( 975 )
  W( 135 ) = W( 135 ) + a*JVS( 976 )
  W( 136 ) = W( 136 ) + a*JVS( 977 )
  W( 137 ) = W( 137 ) + a*JVS( 978 )
  W( 139 ) = W( 139 ) + a*JVS( 979 )
  W( 140 ) = W( 140 ) + a*JVS( 980 )
  a = -W( 131 ) / JVS(         1041  )
  W( 131 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1042 )
  W( 134 ) = W( 134 ) + a*JVS( 1043 )
  W( 135 ) = W( 135 ) + a*JVS( 1044 )
  W( 136 ) = W( 136 ) + a*JVS( 1045 )
  W( 137 ) = W( 137 ) + a*JVS( 1046 )
  W( 139 ) = W( 139 ) + a*JVS( 1047 )
  W( 140 ) = W( 140 ) + a*JVS( 1048 )
  JVS( 1049) = W( 54 )
  JVS( 1050) = W( 81 )
  JVS( 1051) = W( 108 )
  JVS( 1052) = W( 123 )
  JVS( 1053) = W( 129 )
  JVS( 1054) = W( 131 )
  JVS( 1055) = W( 132 )
  JVS( 1056) = W( 133 )
  JVS( 1057) = W( 134 )
  JVS( 1058) = W( 135 )
  JVS( 1059) = W( 136 )
  JVS( 1060) = W( 137 )
  JVS( 1061) = W( 138 )
  JVS( 1062) = W( 139 )
  JVS( 1063) = W( 140 )
  IF ( ABS(  JVS( 1090 )) < TINY(a) ) THEN
         IER = 133                                     
         RETURN
  END IF
   W( 40 ) = JVS( 1064 )
   W( 43 ) = JVS( 1065 )
   W( 56 ) = JVS( 1066 )
   W( 73 ) = JVS( 1067 )
   W( 75 ) = JVS( 1068 )
   W( 81 ) = JVS( 1069 )
   W( 92 ) = JVS( 1070 )
   W( 103 ) = JVS( 1071 )
   W( 104 ) = JVS( 1072 )
   W( 105 ) = JVS( 1073 )
   W( 107 ) = JVS( 1074 )
   W( 108 ) = JVS( 1075 )
   W( 109 ) = JVS( 1076 )
   W( 113 ) = JVS( 1077 )
   W( 114 ) = JVS( 1078 )
   W( 117 ) = JVS( 1079 )
   W( 120 ) = JVS( 1080 )
   W( 122 ) = JVS( 1081 )
   W( 124 ) = JVS( 1082 )
   W( 126 ) = JVS( 1083 )
   W( 127 ) = JVS( 1084 )
   W( 128 ) = JVS( 1085 )
   W( 129 ) = JVS( 1086 )
   W( 130 ) = JVS( 1087 )
   W( 131 ) = JVS( 1088 )
   W( 132 ) = JVS( 1089 )
   W( 133 ) = JVS( 1090 )
   W( 134 ) = JVS( 1091 )
   W( 135 ) = JVS( 1092 )
   W( 136 ) = JVS( 1093 )
   W( 137 ) = JVS( 1094 )
   W( 138 ) = JVS( 1095 )
   W( 139 ) = JVS( 1096 )
   W( 140 ) = JVS( 1097 )
  a = -W( 40 ) / JVS(          219  )
  W( 40 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 220 )
  W( 138 ) = W( 138 ) + a*JVS( 221 )
  a = -W( 43 ) / JVS(          227  )
  W( 43 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 228 )
  W( 134 ) = W( 134 ) + a*JVS( 229 )
  a = -W( 56 ) / JVS(          281  )
  W( 56 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 282 )
  W( 135 ) = W( 135 ) + a*JVS( 283 )
  W( 138 ) = W( 138 ) + a*JVS( 284 )
  a = -W( 73 ) / JVS(          380  )
  W( 73 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 381 )
  W( 137 ) = W( 137 ) + a*JVS( 382 )
  W( 138 ) = W( 138 ) + a*JVS( 383 )
  a = -W( 75 ) / JVS(          392  )
  W( 75 ) = -a
  W( 114 ) = W( 114 ) + a*JVS( 393 )
  W( 126 ) = W( 126 ) + a*JVS( 394 )
  W( 130 ) = W( 130 ) + a*JVS( 395 )
  W( 133 ) = W( 133 ) + a*JVS( 396 )
  W( 134 ) = W( 134 ) + a*JVS( 397 )
  W( 138 ) = W( 138 ) + a*JVS( 398 )
  a = -W( 81 ) / JVS(          425  )
  W( 81 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 426 )
  W( 134 ) = W( 134 ) + a*JVS( 427 )
  W( 138 ) = W( 138 ) + a*JVS( 428 )
  a = -W( 92 ) / JVS(          498  )
  W( 92 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 499 )
  W( 134 ) = W( 134 ) + a*JVS( 500 )
  W( 136 ) = W( 136 ) + a*JVS( 501 )
  a = -W( 103 ) / JVS(          600  )
  W( 103 ) = -a
  W( 108 ) = W( 108 ) + a*JVS( 601 )
  W( 133 ) = W( 133 ) + a*JVS( 602 )
  W( 135 ) = W( 135 ) + a*JVS( 603 )
  W( 139 ) = W( 139 ) + a*JVS( 604 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 105 ) / JVS(          610  )
  W( 105 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 611 )
  W( 134 ) = W( 134 ) + a*JVS( 612 )
  W( 135 ) = W( 135 ) + a*JVS( 613 )
  W( 136 ) = W( 136 ) + a*JVS( 614 )
  W( 139 ) = W( 139 ) + a*JVS( 615 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  a = -W( 113 ) / JVS(          732  )
  W( 113 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 733 )
  W( 133 ) = W( 133 ) + a*JVS( 734 )
  W( 134 ) = W( 134 ) + a*JVS( 735 )
  W( 136 ) = W( 136 ) + a*JVS( 736 )
  a = -W( 114 ) / JVS(          748  )
  W( 114 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 749 )
  W( 128 ) = W( 128 ) + a*JVS( 750 )
  W( 133 ) = W( 133 ) + a*JVS( 751 )
  W( 134 ) = W( 134 ) + a*JVS( 752 )
  W( 135 ) = W( 135 ) + a*JVS( 753 )
  W( 136 ) = W( 136 ) + a*JVS( 754 )
  W( 139 ) = W( 139 ) + a*JVS( 755 )
  W( 140 ) = W( 140 ) + a*JVS( 756 )
  a = -W( 117 ) / JVS(          788  )
  W( 117 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 789 )
  W( 134 ) = W( 134 ) + a*JVS( 790 )
  W( 135 ) = W( 135 ) + a*JVS( 791 )
  W( 136 ) = W( 136 ) + a*JVS( 792 )
  W( 139 ) = W( 139 ) + a*JVS( 793 )
  a = -W( 120 ) / JVS(          819  )
  W( 120 ) = -a
  W( 122 ) = W( 122 ) + a*JVS( 820 )
  W( 124 ) = W( 124 ) + a*JVS( 821 )
  W( 129 ) = W( 129 ) + a*JVS( 822 )
  W( 131 ) = W( 131 ) + a*JVS( 823 )
  W( 133 ) = W( 133 ) + a*JVS( 824 )
  W( 134 ) = W( 134 ) + a*JVS( 825 )
  W( 135 ) = W( 135 ) + a*JVS( 826 )
  W( 136 ) = W( 136 ) + a*JVS( 827 )
  W( 139 ) = W( 139 ) + a*JVS( 828 )
  W( 140 ) = W( 140 ) + a*JVS( 829 )
  a = -W( 122 ) / JVS(          853  )
  W( 122 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 854 )
  W( 133 ) = W( 133 ) + a*JVS( 855 )
  W( 134 ) = W( 134 ) + a*JVS( 856 )
  W( 135 ) = W( 135 ) + a*JVS( 857 )
  W( 136 ) = W( 136 ) + a*JVS( 858 )
  W( 139 ) = W( 139 ) + a*JVS( 859 )
  W( 140 ) = W( 140 ) + a*JVS( 860 )
  a = -W( 124 ) / JVS(          881  )
  W( 124 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 882 )
  W( 134 ) = W( 134 ) + a*JVS( 883 )
  W( 135 ) = W( 135 ) + a*JVS( 884 )
  W( 136 ) = W( 136 ) + a*JVS( 885 )
  W( 139 ) = W( 139 ) + a*JVS( 886 )
  W( 140 ) = W( 140 ) + a*JVS( 887 )
  a = -W( 126 ) / JVS(          927  )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 928 )
  W( 128 ) = W( 128 ) + a*JVS( 929 )
  W( 129 ) = W( 129 ) + a*JVS( 930 )
  W( 131 ) = W( 131 ) + a*JVS( 931 )
  W( 132 ) = W( 132 ) + a*JVS( 932 )
  W( 133 ) = W( 133 ) + a*JVS( 933 )
  W( 134 ) = W( 134 ) + a*JVS( 934 )
  W( 135 ) = W( 135 ) + a*JVS( 935 )
  W( 136 ) = W( 136 ) + a*JVS( 936 )
  W( 137 ) = W( 137 ) + a*JVS( 937 )
  W( 138 ) = W( 138 ) + a*JVS( 938 )
  W( 139 ) = W( 139 ) + a*JVS( 939 )
  W( 140 ) = W( 140 ) + a*JVS( 940 )
  a = -W( 127 ) / JVS(          947  )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 948 )
  W( 131 ) = W( 131 ) + a*JVS( 949 )
  W( 133 ) = W( 133 ) + a*JVS( 950 )
  W( 134 ) = W( 134 ) + a*JVS( 951 )
  W( 135 ) = W( 135 ) + a*JVS( 952 )
  W( 136 ) = W( 136 ) + a*JVS( 953 )
  W( 137 ) = W( 137 ) + a*JVS( 954 )
  W( 139 ) = W( 139 ) + a*JVS( 955 )
  W( 140 ) = W( 140 ) + a*JVS( 956 )
  a = -W( 128 ) / JVS(          960  )
  W( 128 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 961 )
  W( 133 ) = W( 133 ) + a*JVS( 962 )
  W( 134 ) = W( 134 ) + a*JVS( 963 )
  W( 135 ) = W( 135 ) + a*JVS( 964 )
  W( 136 ) = W( 136 ) + a*JVS( 965 )
  W( 137 ) = W( 137 ) + a*JVS( 966 )
  W( 139 ) = W( 139 ) + a*JVS( 967 )
  W( 140 ) = W( 140 ) + a*JVS( 968 )
  a = -W( 129 ) / JVS(          972  )
  W( 129 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 973 )
  W( 133 ) = W( 133 ) + a*JVS( 974 )
  W( 134 ) = W( 134 ) + a*JVS( 975 )
  W( 135 ) = W( 135 ) + a*JVS( 976 )
  W( 136 ) = W( 136 ) + a*JVS( 977 )
  W( 137 ) = W( 137 ) + a*JVS( 978 )
  W( 139 ) = W( 139 ) + a*JVS( 979 )
  W( 140 ) = W( 140 ) + a*JVS( 980 )
  a = -W( 130 ) / JVS(         1028  )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 1029 )
  W( 132 ) = W( 132 ) + a*JVS( 1030 )
  W( 133 ) = W( 133 ) + a*JVS( 1031 )
  W( 134 ) = W( 134 ) + a*JVS( 1032 )
  W( 135 ) = W( 135 ) + a*JVS( 1033 )
  W( 136 ) = W( 136 ) + a*JVS( 1034 )
  W( 137 ) = W( 137 ) + a*JVS( 1035 )
  W( 138 ) = W( 138 ) + a*JVS( 1036 )
  W( 139 ) = W( 139 ) + a*JVS( 1037 )
  W( 140 ) = W( 140 ) + a*JVS( 1038 )
  a = -W( 131 ) / JVS(         1041  )
  W( 131 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1042 )
  W( 134 ) = W( 134 ) + a*JVS( 1043 )
  W( 135 ) = W( 135 ) + a*JVS( 1044 )
  W( 136 ) = W( 136 ) + a*JVS( 1045 )
  W( 137 ) = W( 137 ) + a*JVS( 1046 )
  W( 139 ) = W( 139 ) + a*JVS( 1047 )
  W( 140 ) = W( 140 ) + a*JVS( 1048 )
  a = -W( 132 ) / JVS(         1055  )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1056 )
  W( 134 ) = W( 134 ) + a*JVS( 1057 )
  W( 135 ) = W( 135 ) + a*JVS( 1058 )
  W( 136 ) = W( 136 ) + a*JVS( 1059 )
  W( 137 ) = W( 137 ) + a*JVS( 1060 )
  W( 138 ) = W( 138 ) + a*JVS( 1061 )
  W( 139 ) = W( 139 ) + a*JVS( 1062 )
  W( 140 ) = W( 140 ) + a*JVS( 1063 )
  JVS( 1064) = W( 40 )
  JVS( 1065) = W( 43 )
  JVS( 1066) = W( 56 )
  JVS( 1067) = W( 73 )
  JVS( 1068) = W( 75 )
  JVS( 1069) = W( 81 )
  JVS( 1070) = W( 92 )
  JVS( 1071) = W( 103 )
  JVS( 1072) = W( 104 )
  JVS( 1073) = W( 105 )
  JVS( 1074) = W( 107 )
  JVS( 1075) = W( 108 )
  JVS( 1076) = W( 109 )
  JVS( 1077) = W( 113 )
  JVS( 1078) = W( 114 )
  JVS( 1079) = W( 117 )
  JVS( 1080) = W( 120 )
  JVS( 1081) = W( 122 )
  JVS( 1082) = W( 124 )
  JVS( 1083) = W( 126 )
  JVS( 1084) = W( 127 )
  JVS( 1085) = W( 128 )
  JVS( 1086) = W( 129 )
  JVS( 1087) = W( 130 )
  JVS( 1088) = W( 131 )
  JVS( 1089) = W( 132 )
  JVS( 1090) = W( 133 )
  JVS( 1091) = W( 134 )
  JVS( 1092) = W( 135 )
  JVS( 1093) = W( 136 )
  JVS( 1094) = W( 137 )
  JVS( 1095) = W( 138 )
  JVS( 1096) = W( 139 )
  JVS( 1097) = W( 140 )
  IF ( ABS(  JVS( 1203 )) < TINY(a) ) THEN
         IER = 134                                     
         RETURN
  END IF
   W( 1 ) = JVS( 1098 )
   W( 15 ) = JVS( 1099 )
   W( 16 ) = JVS( 1100 )
   W( 17 ) = JVS( 1101 )
   W( 18 ) = JVS( 1102 )
   W( 20 ) = JVS( 1103 )
   W( 22 ) = JVS( 1104 )
   W( 23 ) = JVS( 1105 )
   W( 24 ) = JVS( 1106 )
   W( 26 ) = JVS( 1107 )
   W( 27 ) = JVS( 1108 )
   W( 28 ) = JVS( 1109 )
   W( 29 ) = JVS( 1110 )
   W( 31 ) = JVS( 1111 )
   W( 32 ) = JVS( 1112 )
   W( 33 ) = JVS( 1113 )
   W( 36 ) = JVS( 1114 )
   W( 37 ) = JVS( 1115 )
   W( 38 ) = JVS( 1116 )
   W( 39 ) = JVS( 1117 )
   W( 41 ) = JVS( 1118 )
   W( 42 ) = JVS( 1119 )
   W( 43 ) = JVS( 1120 )
   W( 45 ) = JVS( 1121 )
   W( 48 ) = JVS( 1122 )
   W( 49 ) = JVS( 1123 )
   W( 50 ) = JVS( 1124 )
   W( 52 ) = JVS( 1125 )
   W( 53 ) = JVS( 1126 )
   W( 54 ) = JVS( 1127 )
   W( 55 ) = JVS( 1128 )
   W( 56 ) = JVS( 1129 )
   W( 57 ) = JVS( 1130 )
   W( 58 ) = JVS( 1131 )
   W( 59 ) = JVS( 1132 )
   W( 61 ) = JVS( 1133 )
   W( 62 ) = JVS( 1134 )
   W( 64 ) = JVS( 1135 )
   W( 66 ) = JVS( 1136 )
   W( 67 ) = JVS( 1137 )
   W( 68 ) = JVS( 1138 )
   W( 69 ) = JVS( 1139 )
   W( 70 ) = JVS( 1140 )
   W( 71 ) = JVS( 1141 )
   W( 72 ) = JVS( 1142 )
   W( 73 ) = JVS( 1143 )
   W( 74 ) = JVS( 1144 )
   W( 75 ) = JVS( 1145 )
   W( 76 ) = JVS( 1146 )
   W( 77 ) = JVS( 1147 )
   W( 78 ) = JVS( 1148 )
   W( 79 ) = JVS( 1149 )
   W( 80 ) = JVS( 1150 )
   W( 81 ) = JVS( 1151 )
   W( 82 ) = JVS( 1152 )
   W( 83 ) = JVS( 1153 )
   W( 84 ) = JVS( 1154 )
   W( 85 ) = JVS( 1155 )
   W( 86 ) = JVS( 1156 )
   W( 87 ) = JVS( 1157 )
   W( 88 ) = JVS( 1158 )
   W( 89 ) = JVS( 1159 )
   W( 90 ) = JVS( 1160 )
   W( 92 ) = JVS( 1161 )
   W( 93 ) = JVS( 1162 )
   W( 94 ) = JVS( 1163 )
   W( 95 ) = JVS( 1164 )
   W( 96 ) = JVS( 1165 )
   W( 97 ) = JVS( 1166 )
   W( 98 ) = JVS( 1167 )
   W( 99 ) = JVS( 1168 )
   W( 100 ) = JVS( 1169 )
   W( 101 ) = JVS( 1170 )
   W( 102 ) = JVS( 1171 )
   W( 103 ) = JVS( 1172 )
   W( 104 ) = JVS( 1173 )
   W( 105 ) = JVS( 1174 )
   W( 106 ) = JVS( 1175 )
   W( 107 ) = JVS( 1176 )
   W( 108 ) = JVS( 1177 )
   W( 109 ) = JVS( 1178 )
   W( 110 ) = JVS( 1179 )
   W( 111 ) = JVS( 1180 )
   W( 112 ) = JVS( 1181 )
   W( 113 ) = JVS( 1182 )
   W( 114 ) = JVS( 1183 )
   W( 115 ) = JVS( 1184 )
   W( 116 ) = JVS( 1185 )
   W( 117 ) = JVS( 1186 )
   W( 118 ) = JVS( 1187 )
   W( 119 ) = JVS( 1188 )
   W( 120 ) = JVS( 1189 )
   W( 121 ) = JVS( 1190 )
   W( 122 ) = JVS( 1191 )
   W( 123 ) = JVS( 1192 )
   W( 124 ) = JVS( 1193 )
   W( 125 ) = JVS( 1194 )
   W( 126 ) = JVS( 1195 )
   W( 127 ) = JVS( 1196 )
   W( 128 ) = JVS( 1197 )
   W( 129 ) = JVS( 1198 )
   W( 130 ) = JVS( 1199 )
   W( 131 ) = JVS( 1200 )
   W( 132 ) = JVS( 1201 )
   W( 133 ) = JVS( 1202 )
   W( 134 ) = JVS( 1203 )
   W( 135 ) = JVS( 1204 )
   W( 136 ) = JVS( 1205 )
   W( 137 ) = JVS( 1206 )
   W( 138 ) = JVS( 1207 )
   W( 139 ) = JVS( 1208 )
   W( 140 ) = JVS( 1209 )
  a = -W( 1 ) / JVS(            1  )
  W( 1 ) = -a
  a = -W( 15 ) / JVS(          141  )
  W( 15 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 142 )
  a = -W( 16 ) / JVS(          143  )
  W( 16 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 144 )
  a = -W( 17 ) / JVS(          145  )
  W( 17 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 146 )
  a = -W( 18 ) / JVS(          147  )
  W( 18 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 148 )
  a = -W( 20 ) / JVS(          151  )
  W( 20 ) = -a
  W( 43 ) = W( 43 ) + a*JVS( 152 )
  W( 133 ) = W( 133 ) + a*JVS( 153 )
  W( 134 ) = W( 134 ) + a*JVS( 154 )
  a = -W( 22 ) / JVS(          158  )
  W( 22 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 159 )
  W( 135 ) = W( 135 ) + a*JVS( 160 )
  a = -W( 23 ) / JVS(          161  )
  W( 23 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 162 )
  a = -W( 24 ) / JVS(          163  )
  W( 24 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 164 )
  a = -W( 26 ) / JVS(          169  )
  W( 26 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 170 )
  a = -W( 27 ) / JVS(          171  )
  W( 27 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 172 )
  a = -W( 28 ) / JVS(          174  )
  W( 28 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 175 )
  a = -W( 29 ) / JVS(          176  )
  W( 29 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 177 )
  W( 134 ) = W( 134 ) + a*JVS( 178 )
  W( 135 ) = W( 135 ) + a*JVS( 179 )
  a = -W( 31 ) / JVS(          183  )
  W( 31 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 184 )
  W( 134 ) = W( 134 ) + a*JVS( 185 )
  W( 135 ) = W( 135 ) + a*JVS( 186 )
  a = -W( 32 ) / JVS(          187  )
  W( 32 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 188 )
  W( 134 ) = W( 134 ) + a*JVS( 189 )
  W( 135 ) = W( 135 ) + a*JVS( 190 )
  a = -W( 33 ) / JVS(          191  )
  W( 33 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 192 )
  W( 134 ) = W( 134 ) + a*JVS( 193 )
  W( 135 ) = W( 135 ) + a*JVS( 194 )
  a = -W( 36 ) / JVS(          203  )
  W( 36 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 204 )
  W( 134 ) = W( 134 ) + a*JVS( 205 )
  W( 135 ) = W( 135 ) + a*JVS( 206 )
  a = -W( 37 ) / JVS(          207  )
  W( 37 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 208 )
  W( 134 ) = W( 134 ) + a*JVS( 209 )
  W( 139 ) = W( 139 ) + a*JVS( 210 )
  a = -W( 38 ) / JVS(          211  )
  W( 38 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 212 )
  W( 134 ) = W( 134 ) + a*JVS( 213 )
  W( 135 ) = W( 135 ) + a*JVS( 214 )
  a = -W( 39 ) / JVS(          215  )
  W( 39 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 216 )
  W( 134 ) = W( 134 ) + a*JVS( 217 )
  W( 135 ) = W( 135 ) + a*JVS( 218 )
  a = -W( 41 ) / JVS(          222  )
  W( 41 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 223 )
  a = -W( 42 ) / JVS(          225  )
  W( 42 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 226 )
  a = -W( 43 ) / JVS(          227  )
  W( 43 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 228 )
  W( 134 ) = W( 134 ) + a*JVS( 229 )
  a = -W( 45 ) / JVS(          233  )
  W( 45 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 234 )
  W( 134 ) = W( 134 ) + a*JVS( 235 )
  W( 140 ) = W( 140 ) + a*JVS( 236 )
  a = -W( 48 ) / JVS(          248  )
  W( 48 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 249 )
  W( 134 ) = W( 134 ) + a*JVS( 250 )
  W( 135 ) = W( 135 ) + a*JVS( 251 )
  a = -W( 49 ) / JVS(          252  )
  W( 49 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 253 )
  W( 134 ) = W( 134 ) + a*JVS( 254 )
  W( 135 ) = W( 135 ) + a*JVS( 255 )
  a = -W( 50 ) / JVS(          256  )
  W( 50 ) = -a
  W( 121 ) = W( 121 ) + a*JVS( 257 )
  W( 134 ) = W( 134 ) + a*JVS( 258 )
  W( 135 ) = W( 135 ) + a*JVS( 259 )
  a = -W( 52 ) / JVS(          264  )
  W( 52 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 265 )
  W( 134 ) = W( 134 ) + a*JVS( 266 )
  W( 135 ) = W( 135 ) + a*JVS( 267 )
  a = -W( 53 ) / JVS(          268  )
  W( 53 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 269 )
  W( 135 ) = W( 135 ) + a*JVS( 270 )
  W( 140 ) = W( 140 ) + a*JVS( 271 )
  a = -W( 54 ) / JVS(          272  )
  W( 54 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 273 )
  W( 134 ) = W( 134 ) + a*JVS( 274 )
  W( 135 ) = W( 135 ) + a*JVS( 275 )
  a = -W( 55 ) / JVS(          276  )
  W( 55 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 277 )
  W( 84 ) = W( 84 ) + a*JVS( 278 )
  W( 130 ) = W( 130 ) + a*JVS( 279 )
  W( 134 ) = W( 134 ) + a*JVS( 280 )
  a = -W( 56 ) / JVS(          281  )
  W( 56 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 282 )
  W( 135 ) = W( 135 ) + a*JVS( 283 )
  W( 138 ) = W( 138 ) + a*JVS( 284 )
  a = -W( 57 ) / JVS(          285  )
  W( 57 ) = -a
  W( 120 ) = W( 120 ) + a*JVS( 286 )
  W( 131 ) = W( 131 ) + a*JVS( 287 )
  W( 133 ) = W( 133 ) + a*JVS( 288 )
  W( 134 ) = W( 134 ) + a*JVS( 289 )
  W( 137 ) = W( 137 ) + a*JVS( 290 )
  W( 139 ) = W( 139 ) + a*JVS( 291 )
  W( 140 ) = W( 140 ) + a*JVS( 292 )
  a = -W( 58 ) / JVS(          293  )
  W( 58 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 294 )
  W( 87 ) = W( 87 ) + a*JVS( 295 )
  W( 104 ) = W( 104 ) + a*JVS( 296 )
  W( 107 ) = W( 107 ) + a*JVS( 297 )
  W( 109 ) = W( 109 ) + a*JVS( 298 )
  W( 134 ) = W( 134 ) + a*JVS( 299 )
  W( 136 ) = W( 136 ) + a*JVS( 300 )
  a = -W( 59 ) / JVS(          301  )
  W( 59 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 302 )
  W( 136 ) = W( 136 ) + a*JVS( 303 )
  a = -W( 61 ) / JVS(          309  )
  W( 61 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 310 )
  W( 134 ) = W( 134 ) + a*JVS( 311 )
  W( 135 ) = W( 135 ) + a*JVS( 312 )
  a = -W( 62 ) / JVS(          316  )
  W( 62 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 317 )
  W( 135 ) = W( 135 ) + a*JVS( 318 )
  W( 139 ) = W( 139 ) + a*JVS( 319 )
  a = -W( 64 ) / JVS(          332  )
  W( 64 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 333 )
  W( 135 ) = W( 135 ) + a*JVS( 334 )
  W( 139 ) = W( 139 ) + a*JVS( 335 )
  a = -W( 66 ) / JVS(          340  )
  W( 66 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 341 )
  W( 117 ) = W( 117 ) + a*JVS( 342 )
  W( 134 ) = W( 134 ) + a*JVS( 343 )
  W( 139 ) = W( 139 ) + a*JVS( 344 )
  a = -W( 67 ) / JVS(          345  )
  W( 67 ) = -a
  W( 108 ) = W( 108 ) + a*JVS( 346 )
  W( 113 ) = W( 113 ) + a*JVS( 347 )
  W( 132 ) = W( 132 ) + a*JVS( 348 )
  W( 134 ) = W( 134 ) + a*JVS( 349 )
  W( 135 ) = W( 135 ) + a*JVS( 350 )
  W( 136 ) = W( 136 ) + a*JVS( 351 )
  W( 137 ) = W( 137 ) + a*JVS( 352 )
  W( 140 ) = W( 140 ) + a*JVS( 353 )
  a = -W( 68 ) / JVS(          354  )
  W( 68 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 355 )
  W( 134 ) = W( 134 ) + a*JVS( 356 )
  W( 135 ) = W( 135 ) + a*JVS( 357 )
  W( 137 ) = W( 137 ) + a*JVS( 358 )
  a = -W( 69 ) / JVS(          359  )
  W( 69 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 360 )
  W( 134 ) = W( 134 ) + a*JVS( 361 )
  W( 135 ) = W( 135 ) + a*JVS( 362 )
  a = -W( 70 ) / JVS(          365  )
  W( 70 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 366 )
  W( 135 ) = W( 135 ) + a*JVS( 367 )
  W( 139 ) = W( 139 ) + a*JVS( 368 )
  a = -W( 71 ) / JVS(          369  )
  W( 71 ) = -a
  W( 105 ) = W( 105 ) + a*JVS( 370 )
  W( 113 ) = W( 113 ) + a*JVS( 371 )
  W( 117 ) = W( 117 ) + a*JVS( 372 )
  W( 133 ) = W( 133 ) + a*JVS( 373 )
  W( 134 ) = W( 134 ) + a*JVS( 374 )
  W( 139 ) = W( 139 ) + a*JVS( 375 )
  a = -W( 72 ) / JVS(          376  )
  W( 72 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 377 )
  W( 134 ) = W( 134 ) + a*JVS( 378 )
  W( 135 ) = W( 135 ) + a*JVS( 379 )
  a = -W( 73 ) / JVS(          380  )
  W( 73 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 381 )
  W( 137 ) = W( 137 ) + a*JVS( 382 )
  W( 138 ) = W( 138 ) + a*JVS( 383 )
  a = -W( 74 ) / JVS(          386  )
  W( 74 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 387 )
  W( 135 ) = W( 135 ) + a*JVS( 388 )
  W( 139 ) = W( 139 ) + a*JVS( 389 )
  a = -W( 75 ) / JVS(          392  )
  W( 75 ) = -a
  W( 114 ) = W( 114 ) + a*JVS( 393 )
  W( 126 ) = W( 126 ) + a*JVS( 394 )
  W( 130 ) = W( 130 ) + a*JVS( 395 )
  W( 133 ) = W( 133 ) + a*JVS( 396 )
  W( 134 ) = W( 134 ) + a*JVS( 397 )
  W( 138 ) = W( 138 ) + a*JVS( 398 )
  a = -W( 76 ) / JVS(          399  )
  W( 76 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 400 )
  W( 134 ) = W( 134 ) + a*JVS( 401 )
  W( 135 ) = W( 135 ) + a*JVS( 402 )
  W( 139 ) = W( 139 ) + a*JVS( 403 )
  a = -W( 77 ) / JVS(          404  )
  W( 77 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 405 )
  W( 113 ) = W( 113 ) + a*JVS( 406 )
  W( 134 ) = W( 134 ) + a*JVS( 407 )
  W( 136 ) = W( 136 ) + a*JVS( 408 )
  a = -W( 78 ) / JVS(          412  )
  W( 78 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 413 )
  W( 135 ) = W( 135 ) + a*JVS( 414 )
  W( 138 ) = W( 138 ) + a*JVS( 415 )
  W( 139 ) = W( 139 ) + a*JVS( 416 )
  a = -W( 79 ) / JVS(          417  )
  W( 79 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 418 )
  W( 134 ) = W( 134 ) + a*JVS( 419 )
  W( 135 ) = W( 135 ) + a*JVS( 420 )
  a = -W( 80 ) / JVS(          421  )
  W( 80 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 422 )
  W( 134 ) = W( 134 ) + a*JVS( 423 )
  W( 135 ) = W( 135 ) + a*JVS( 424 )
  a = -W( 81 ) / JVS(          425  )
  W( 81 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 426 )
  W( 134 ) = W( 134 ) + a*JVS( 427 )
  W( 138 ) = W( 138 ) + a*JVS( 428 )
  a = -W( 82 ) / JVS(          432  )
  W( 82 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 433 )
  W( 87 ) = W( 87 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 104 ) = W( 104 ) + a*JVS( 436 )
  W( 107 ) = W( 107 ) + a*JVS( 437 )
  W( 109 ) = W( 109 ) + a*JVS( 438 )
  W( 130 ) = W( 130 ) + a*JVS( 439 )
  W( 134 ) = W( 134 ) + a*JVS( 440 )
  W( 135 ) = W( 135 ) + a*JVS( 441 )
  W( 136 ) = W( 136 ) + a*JVS( 442 )
  W( 139 ) = W( 139 ) + a*JVS( 443 )
  a = -W( 83 ) / JVS(          444  )
  W( 83 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 445 )
  W( 134 ) = W( 134 ) + a*JVS( 446 )
  W( 135 ) = W( 135 ) + a*JVS( 447 )
  a = -W( 84 ) / JVS(          451  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 452 )
  W( 130 ) = W( 130 ) + a*JVS( 453 )
  W( 134 ) = W( 134 ) + a*JVS( 454 )
  W( 136 ) = W( 136 ) + a*JVS( 455 )
  a = -W( 85 ) / JVS(          459  )
  W( 85 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 460 )
  W( 135 ) = W( 135 ) + a*JVS( 461 )
  W( 136 ) = W( 136 ) + a*JVS( 462 )
  W( 138 ) = W( 138 ) + a*JVS( 463 )
  W( 139 ) = W( 139 ) + a*JVS( 464 )
  a = -W( 86 ) / JVS(          465  )
  W( 86 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 466 )
  W( 136 ) = W( 136 ) + a*JVS( 467 )
  a = -W( 87 ) / JVS(          468  )
  W( 87 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 469 )
  W( 136 ) = W( 136 ) + a*JVS( 470 )
  a = -W( 88 ) / JVS(          471  )
  W( 88 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 472 )
  W( 134 ) = W( 134 ) + a*JVS( 473 )
  W( 135 ) = W( 135 ) + a*JVS( 474 )
  a = -W( 89 ) / JVS(          478  )
  W( 89 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 479 )
  W( 134 ) = W( 134 ) + a*JVS( 480 )
  W( 135 ) = W( 135 ) + a*JVS( 481 )
  W( 139 ) = W( 139 ) + a*JVS( 482 )
  a = -W( 90 ) / JVS(          484  )
  W( 90 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 485 )
  W( 134 ) = W( 134 ) + a*JVS( 486 )
  W( 135 ) = W( 135 ) + a*JVS( 487 )
  W( 139 ) = W( 139 ) + a*JVS( 488 )
  a = -W( 92 ) / JVS(          498  )
  W( 92 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 499 )
  W( 134 ) = W( 134 ) + a*JVS( 500 )
  W( 136 ) = W( 136 ) + a*JVS( 501 )
  a = -W( 93 ) / JVS(          504  )
  W( 93 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 505 )
  W( 135 ) = W( 135 ) + a*JVS( 506 )
  W( 139 ) = W( 139 ) + a*JVS( 507 )
  a = -W( 94 ) / JVS(          510  )
  W( 94 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 511 )
  W( 135 ) = W( 135 ) + a*JVS( 512 )
  W( 139 ) = W( 139 ) + a*JVS( 513 )
  a = -W( 95 ) / JVS(          517  )
  W( 95 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 518 )
  W( 135 ) = W( 135 ) + a*JVS( 519 )
  W( 138 ) = W( 138 ) + a*JVS( 520 )
  W( 139 ) = W( 139 ) + a*JVS( 521 )
  a = -W( 96 ) / JVS(          523  )
  W( 96 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 524 )
  W( 130 ) = W( 130 ) + a*JVS( 525 )
  W( 134 ) = W( 134 ) + a*JVS( 526 )
  W( 135 ) = W( 135 ) + a*JVS( 527 )
  W( 136 ) = W( 136 ) + a*JVS( 528 )
  W( 138 ) = W( 138 ) + a*JVS( 529 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 533 )
  W( 135 ) = W( 135 ) + a*JVS( 534 )
  W( 139 ) = W( 139 ) + a*JVS( 535 )
  W( 140 ) = W( 140 ) + a*JVS( 536 )
  a = -W( 98 ) / JVS(          540  )
  W( 98 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 541 )
  W( 135 ) = W( 135 ) + a*JVS( 542 )
  W( 138 ) = W( 138 ) + a*JVS( 543 )
  W( 139 ) = W( 139 ) + a*JVS( 544 )
  a = -W( 99 ) / JVS(          554  )
  W( 99 ) = -a
  W( 101 ) = W( 101 ) + a*JVS( 555 )
  W( 112 ) = W( 112 ) + a*JVS( 556 )
  W( 127 ) = W( 127 ) + a*JVS( 557 )
  W( 131 ) = W( 131 ) + a*JVS( 558 )
  W( 133 ) = W( 133 ) + a*JVS( 559 )
  W( 134 ) = W( 134 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 137 ) = W( 137 ) + a*JVS( 562 )
  W( 138 ) = W( 138 ) + a*JVS( 563 )
  W( 139 ) = W( 139 ) + a*JVS( 564 )
  W( 140 ) = W( 140 ) + a*JVS( 565 )
  a = -W( 100 ) / JVS(          566  )
  W( 100 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 567 )
  W( 116 ) = W( 116 ) + a*JVS( 568 )
  W( 118 ) = W( 118 ) + a*JVS( 569 )
  W( 119 ) = W( 119 ) + a*JVS( 570 )
  W( 121 ) = W( 121 ) + a*JVS( 571 )
  W( 122 ) = W( 122 ) + a*JVS( 572 )
  W( 127 ) = W( 127 ) + a*JVS( 573 )
  W( 129 ) = W( 129 ) + a*JVS( 574 )
  W( 131 ) = W( 131 ) + a*JVS( 575 )
  W( 134 ) = W( 134 ) + a*JVS( 576 )
  W( 140 ) = W( 140 ) + a*JVS( 577 )
  a = -W( 101 ) / JVS(          582  )
  W( 101 ) = -a
  W( 104 ) = W( 104 ) + a*JVS( 583 )
  W( 107 ) = W( 107 ) + a*JVS( 584 )
  W( 109 ) = W( 109 ) + a*JVS( 585 )
  W( 134 ) = W( 134 ) + a*JVS( 586 )
  W( 135 ) = W( 135 ) + a*JVS( 587 )
  W( 136 ) = W( 136 ) + a*JVS( 588 )
  W( 139 ) = W( 139 ) + a*JVS( 589 )
  a = -W( 102 ) / JVS(          594  )
  W( 102 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 595 )
  W( 134 ) = W( 134 ) + a*JVS( 596 )
  W( 135 ) = W( 135 ) + a*JVS( 597 )
  W( 139 ) = W( 139 ) + a*JVS( 598 )
  W( 140 ) = W( 140 ) + a*JVS( 599 )
  a = -W( 103 ) / JVS(          600  )
  W( 103 ) = -a
  W( 108 ) = W( 108 ) + a*JVS( 601 )
  W( 133 ) = W( 133 ) + a*JVS( 602 )
  W( 135 ) = W( 135 ) + a*JVS( 603 )
  W( 139 ) = W( 139 ) + a*JVS( 604 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 105 ) / JVS(          610  )
  W( 105 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 611 )
  W( 134 ) = W( 134 ) + a*JVS( 612 )
  W( 135 ) = W( 135 ) + a*JVS( 613 )
  W( 136 ) = W( 136 ) + a*JVS( 614 )
  W( 139 ) = W( 139 ) + a*JVS( 615 )
  a = -W( 106 ) / JVS(          619  )
  W( 106 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 620 )
  W( 133 ) = W( 133 ) + a*JVS( 621 )
  W( 134 ) = W( 134 ) + a*JVS( 622 )
  W( 135 ) = W( 135 ) + a*JVS( 623 )
  W( 136 ) = W( 136 ) + a*JVS( 624 )
  W( 139 ) = W( 139 ) + a*JVS( 625 )
  W( 140 ) = W( 140 ) + a*JVS( 626 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  a = -W( 110 ) / JVS(          655  )
  W( 110 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 656 )
  W( 116 ) = W( 116 ) + a*JVS( 657 )
  W( 117 ) = W( 117 ) + a*JVS( 658 )
  W( 118 ) = W( 118 ) + a*JVS( 659 )
  W( 119 ) = W( 119 ) + a*JVS( 660 )
  W( 133 ) = W( 133 ) + a*JVS( 661 )
  W( 134 ) = W( 134 ) + a*JVS( 662 )
  W( 135 ) = W( 135 ) + a*JVS( 663 )
  W( 136 ) = W( 136 ) + a*JVS( 664 )
  W( 139 ) = W( 139 ) + a*JVS( 665 )
  W( 140 ) = W( 140 ) + a*JVS( 666 )
  a = -W( 111 ) / JVS(          690  )
  W( 111 ) = -a
  W( 112 ) = W( 112 ) + a*JVS( 691 )
  W( 113 ) = W( 113 ) + a*JVS( 692 )
  W( 114 ) = W( 114 ) + a*JVS( 693 )
  W( 115 ) = W( 115 ) + a*JVS( 694 )
  W( 116 ) = W( 116 ) + a*JVS( 695 )
  W( 118 ) = W( 118 ) + a*JVS( 696 )
  W( 120 ) = W( 120 ) + a*JVS( 697 )
  W( 123 ) = W( 123 ) + a*JVS( 698 )
  W( 124 ) = W( 124 ) + a*JVS( 699 )
  W( 126 ) = W( 126 ) + a*JVS( 700 )
  W( 127 ) = W( 127 ) + a*JVS( 701 )
  W( 128 ) = W( 128 ) + a*JVS( 702 )
  W( 129 ) = W( 129 ) + a*JVS( 703 )
  W( 130 ) = W( 130 ) + a*JVS( 704 )
  W( 131 ) = W( 131 ) + a*JVS( 705 )
  W( 133 ) = W( 133 ) + a*JVS( 706 )
  W( 134 ) = W( 134 ) + a*JVS( 707 )
  W( 135 ) = W( 135 ) + a*JVS( 708 )
  W( 136 ) = W( 136 ) + a*JVS( 709 )
  W( 137 ) = W( 137 ) + a*JVS( 710 )
  W( 138 ) = W( 138 ) + a*JVS( 711 )
  W( 139 ) = W( 139 ) + a*JVS( 712 )
  W( 140 ) = W( 140 ) + a*JVS( 713 )
  a = -W( 112 ) / JVS(          718  )
  W( 112 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 719 )
  W( 118 ) = W( 118 ) + a*JVS( 720 )
  W( 127 ) = W( 127 ) + a*JVS( 721 )
  W( 129 ) = W( 129 ) + a*JVS( 722 )
  W( 131 ) = W( 131 ) + a*JVS( 723 )
  W( 133 ) = W( 133 ) + a*JVS( 724 )
  W( 134 ) = W( 134 ) + a*JVS( 725 )
  W( 135 ) = W( 135 ) + a*JVS( 726 )
  W( 136 ) = W( 136 ) + a*JVS( 727 )
  W( 137 ) = W( 137 ) + a*JVS( 728 )
  W( 139 ) = W( 139 ) + a*JVS( 729 )
  W( 140 ) = W( 140 ) + a*JVS( 730 )
  a = -W( 113 ) / JVS(          732  )
  W( 113 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 733 )
  W( 133 ) = W( 133 ) + a*JVS( 734 )
  W( 134 ) = W( 134 ) + a*JVS( 735 )
  W( 136 ) = W( 136 ) + a*JVS( 736 )
  a = -W( 114 ) / JVS(          748  )
  W( 114 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 749 )
  W( 128 ) = W( 128 ) + a*JVS( 750 )
  W( 133 ) = W( 133 ) + a*JVS( 751 )
  W( 134 ) = W( 134 ) + a*JVS( 752 )
  W( 135 ) = W( 135 ) + a*JVS( 753 )
  W( 136 ) = W( 136 ) + a*JVS( 754 )
  W( 139 ) = W( 139 ) + a*JVS( 755 )
  W( 140 ) = W( 140 ) + a*JVS( 756 )
  a = -W( 115 ) / JVS(          763  )
  W( 115 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 764 )
  W( 124 ) = W( 124 ) + a*JVS( 765 )
  W( 133 ) = W( 133 ) + a*JVS( 766 )
  W( 134 ) = W( 134 ) + a*JVS( 767 )
  W( 135 ) = W( 135 ) + a*JVS( 768 )
  W( 136 ) = W( 136 ) + a*JVS( 769 )
  W( 139 ) = W( 139 ) + a*JVS( 770 )
  W( 140 ) = W( 140 ) + a*JVS( 771 )
  a = -W( 116 ) / JVS(          774  )
  W( 116 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 775 )
  W( 134 ) = W( 134 ) + a*JVS( 776 )
  W( 135 ) = W( 135 ) + a*JVS( 777 )
  W( 136 ) = W( 136 ) + a*JVS( 778 )
  W( 139 ) = W( 139 ) + a*JVS( 779 )
  W( 140 ) = W( 140 ) + a*JVS( 780 )
  a = -W( 117 ) / JVS(          788  )
  W( 117 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 789 )
  W( 134 ) = W( 134 ) + a*JVS( 790 )
  W( 135 ) = W( 135 ) + a*JVS( 791 )
  W( 136 ) = W( 136 ) + a*JVS( 792 )
  W( 139 ) = W( 139 ) + a*JVS( 793 )
  a = -W( 118 ) / JVS(          795  )
  W( 118 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 796 )
  W( 133 ) = W( 133 ) + a*JVS( 797 )
  W( 134 ) = W( 134 ) + a*JVS( 798 )
  W( 135 ) = W( 135 ) + a*JVS( 799 )
  W( 139 ) = W( 139 ) + a*JVS( 800 )
  W( 140 ) = W( 140 ) + a*JVS( 801 )
  a = -W( 119 ) / JVS(          808  )
  W( 119 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 809 )
  W( 134 ) = W( 134 ) + a*JVS( 810 )
  W( 135 ) = W( 135 ) + a*JVS( 811 )
  W( 136 ) = W( 136 ) + a*JVS( 812 )
  W( 139 ) = W( 139 ) + a*JVS( 813 )
  W( 140 ) = W( 140 ) + a*JVS( 814 )
  a = -W( 120 ) / JVS(          819  )
  W( 120 ) = -a
  W( 122 ) = W( 122 ) + a*JVS( 820 )
  W( 124 ) = W( 124 ) + a*JVS( 821 )
  W( 129 ) = W( 129 ) + a*JVS( 822 )
  W( 131 ) = W( 131 ) + a*JVS( 823 )
  W( 133 ) = W( 133 ) + a*JVS( 824 )
  W( 134 ) = W( 134 ) + a*JVS( 825 )
  W( 135 ) = W( 135 ) + a*JVS( 826 )
  W( 136 ) = W( 136 ) + a*JVS( 827 )
  W( 139 ) = W( 139 ) + a*JVS( 828 )
  W( 140 ) = W( 140 ) + a*JVS( 829 )
  a = -W( 121 ) / JVS(          842  )
  W( 121 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 843 )
  W( 133 ) = W( 133 ) + a*JVS( 844 )
  W( 134 ) = W( 134 ) + a*JVS( 845 )
  W( 135 ) = W( 135 ) + a*JVS( 846 )
  W( 136 ) = W( 136 ) + a*JVS( 847 )
  W( 139 ) = W( 139 ) + a*JVS( 848 )
  W( 140 ) = W( 140 ) + a*JVS( 849 )
  a = -W( 122 ) / JVS(          853  )
  W( 122 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 854 )
  W( 133 ) = W( 133 ) + a*JVS( 855 )
  W( 134 ) = W( 134 ) + a*JVS( 856 )
  W( 135 ) = W( 135 ) + a*JVS( 857 )
  W( 136 ) = W( 136 ) + a*JVS( 858 )
  W( 139 ) = W( 139 ) + a*JVS( 859 )
  W( 140 ) = W( 140 ) + a*JVS( 860 )
  a = -W( 123 ) / JVS(          864  )
  W( 123 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 865 )
  W( 133 ) = W( 133 ) + a*JVS( 866 )
  W( 134 ) = W( 134 ) + a*JVS( 867 )
  W( 135 ) = W( 135 ) + a*JVS( 868 )
  W( 136 ) = W( 136 ) + a*JVS( 869 )
  W( 137 ) = W( 137 ) + a*JVS( 870 )
  W( 139 ) = W( 139 ) + a*JVS( 871 )
  W( 140 ) = W( 140 ) + a*JVS( 872 )
  a = -W( 124 ) / JVS(          881  )
  W( 124 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 882 )
  W( 134 ) = W( 134 ) + a*JVS( 883 )
  W( 135 ) = W( 135 ) + a*JVS( 884 )
  W( 136 ) = W( 136 ) + a*JVS( 885 )
  W( 139 ) = W( 139 ) + a*JVS( 886 )
  W( 140 ) = W( 140 ) + a*JVS( 887 )
  a = -W( 125 ) / JVS(          894  )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 895 )
  W( 128 ) = W( 128 ) + a*JVS( 896 )
  W( 129 ) = W( 129 ) + a*JVS( 897 )
  W( 131 ) = W( 131 ) + a*JVS( 898 )
  W( 132 ) = W( 132 ) + a*JVS( 899 )
  W( 133 ) = W( 133 ) + a*JVS( 900 )
  W( 134 ) = W( 134 ) + a*JVS( 901 )
  W( 135 ) = W( 135 ) + a*JVS( 902 )
  W( 136 ) = W( 136 ) + a*JVS( 903 )
  W( 137 ) = W( 137 ) + a*JVS( 904 )
  W( 138 ) = W( 138 ) + a*JVS( 905 )
  W( 139 ) = W( 139 ) + a*JVS( 906 )
  W( 140 ) = W( 140 ) + a*JVS( 907 )
  a = -W( 126 ) / JVS(          927  )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 928 )
  W( 128 ) = W( 128 ) + a*JVS( 929 )
  W( 129 ) = W( 129 ) + a*JVS( 930 )
  W( 131 ) = W( 131 ) + a*JVS( 931 )
  W( 132 ) = W( 132 ) + a*JVS( 932 )
  W( 133 ) = W( 133 ) + a*JVS( 933 )
  W( 134 ) = W( 134 ) + a*JVS( 934 )
  W( 135 ) = W( 135 ) + a*JVS( 935 )
  W( 136 ) = W( 136 ) + a*JVS( 936 )
  W( 137 ) = W( 137 ) + a*JVS( 937 )
  W( 138 ) = W( 138 ) + a*JVS( 938 )
  W( 139 ) = W( 139 ) + a*JVS( 939 )
  W( 140 ) = W( 140 ) + a*JVS( 940 )
  a = -W( 127 ) / JVS(          947  )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 948 )
  W( 131 ) = W( 131 ) + a*JVS( 949 )
  W( 133 ) = W( 133 ) + a*JVS( 950 )
  W( 134 ) = W( 134 ) + a*JVS( 951 )
  W( 135 ) = W( 135 ) + a*JVS( 952 )
  W( 136 ) = W( 136 ) + a*JVS( 953 )
  W( 137 ) = W( 137 ) + a*JVS( 954 )
  W( 139 ) = W( 139 ) + a*JVS( 955 )
  W( 140 ) = W( 140 ) + a*JVS( 956 )
  a = -W( 128 ) / JVS(          960  )
  W( 128 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 961 )
  W( 133 ) = W( 133 ) + a*JVS( 962 )
  W( 134 ) = W( 134 ) + a*JVS( 963 )
  W( 135 ) = W( 135 ) + a*JVS( 964 )
  W( 136 ) = W( 136 ) + a*JVS( 965 )
  W( 137 ) = W( 137 ) + a*JVS( 966 )
  W( 139 ) = W( 139 ) + a*JVS( 967 )
  W( 140 ) = W( 140 ) + a*JVS( 968 )
  a = -W( 129 ) / JVS(          972  )
  W( 129 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 973 )
  W( 133 ) = W( 133 ) + a*JVS( 974 )
  W( 134 ) = W( 134 ) + a*JVS( 975 )
  W( 135 ) = W( 135 ) + a*JVS( 976 )
  W( 136 ) = W( 136 ) + a*JVS( 977 )
  W( 137 ) = W( 137 ) + a*JVS( 978 )
  W( 139 ) = W( 139 ) + a*JVS( 979 )
  W( 140 ) = W( 140 ) + a*JVS( 980 )
  a = -W( 130 ) / JVS(         1028  )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 1029 )
  W( 132 ) = W( 132 ) + a*JVS( 1030 )
  W( 133 ) = W( 133 ) + a*JVS( 1031 )
  W( 134 ) = W( 134 ) + a*JVS( 1032 )
  W( 135 ) = W( 135 ) + a*JVS( 1033 )
  W( 136 ) = W( 136 ) + a*JVS( 1034 )
  W( 137 ) = W( 137 ) + a*JVS( 1035 )
  W( 138 ) = W( 138 ) + a*JVS( 1036 )
  W( 139 ) = W( 139 ) + a*JVS( 1037 )
  W( 140 ) = W( 140 ) + a*JVS( 1038 )
  a = -W( 131 ) / JVS(         1041  )
  W( 131 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1042 )
  W( 134 ) = W( 134 ) + a*JVS( 1043 )
  W( 135 ) = W( 135 ) + a*JVS( 1044 )
  W( 136 ) = W( 136 ) + a*JVS( 1045 )
  W( 137 ) = W( 137 ) + a*JVS( 1046 )
  W( 139 ) = W( 139 ) + a*JVS( 1047 )
  W( 140 ) = W( 140 ) + a*JVS( 1048 )
  a = -W( 132 ) / JVS(         1055  )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1056 )
  W( 134 ) = W( 134 ) + a*JVS( 1057 )
  W( 135 ) = W( 135 ) + a*JVS( 1058 )
  W( 136 ) = W( 136 ) + a*JVS( 1059 )
  W( 137 ) = W( 137 ) + a*JVS( 1060 )
  W( 138 ) = W( 138 ) + a*JVS( 1061 )
  W( 139 ) = W( 139 ) + a*JVS( 1062 )
  W( 140 ) = W( 140 ) + a*JVS( 1063 )
  a = -W( 133 ) / JVS(         1090  )
  W( 133 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 1091 )
  W( 135 ) = W( 135 ) + a*JVS( 1092 )
  W( 136 ) = W( 136 ) + a*JVS( 1093 )
  W( 137 ) = W( 137 ) + a*JVS( 1094 )
  W( 138 ) = W( 138 ) + a*JVS( 1095 )
  W( 139 ) = W( 139 ) + a*JVS( 1096 )
  W( 140 ) = W( 140 ) + a*JVS( 1097 )
  JVS( 1098) = W( 1 )
  JVS( 1099) = W( 15 )
  JVS( 1100) = W( 16 )
  JVS( 1101) = W( 17 )
  JVS( 1102) = W( 18 )
  JVS( 1103) = W( 20 )
  JVS( 1104) = W( 22 )
  JVS( 1105) = W( 23 )
  JVS( 1106) = W( 24 )
  JVS( 1107) = W( 26 )
  JVS( 1108) = W( 27 )
  JVS( 1109) = W( 28 )
  JVS( 1110) = W( 29 )
  JVS( 1111) = W( 31 )
  JVS( 1112) = W( 32 )
  JVS( 1113) = W( 33 )
  JVS( 1114) = W( 36 )
  JVS( 1115) = W( 37 )
  JVS( 1116) = W( 38 )
  JVS( 1117) = W( 39 )
  JVS( 1118) = W( 41 )
  JVS( 1119) = W( 42 )
  JVS( 1120) = W( 43 )
  JVS( 1121) = W( 45 )
  JVS( 1122) = W( 48 )
  JVS( 1123) = W( 49 )
  JVS( 1124) = W( 50 )
  JVS( 1125) = W( 52 )
  JVS( 1126) = W( 53 )
  JVS( 1127) = W( 54 )
  JVS( 1128) = W( 55 )
  JVS( 1129) = W( 56 )
  JVS( 1130) = W( 57 )
  JVS( 1131) = W( 58 )
  JVS( 1132) = W( 59 )
  JVS( 1133) = W( 61 )
  JVS( 1134) = W( 62 )
  JVS( 1135) = W( 64 )
  JVS( 1136) = W( 66 )
  JVS( 1137) = W( 67 )
  JVS( 1138) = W( 68 )
  JVS( 1139) = W( 69 )
  JVS( 1140) = W( 70 )
  JVS( 1141) = W( 71 )
  JVS( 1142) = W( 72 )
  JVS( 1143) = W( 73 )
  JVS( 1144) = W( 74 )
  JVS( 1145) = W( 75 )
  JVS( 1146) = W( 76 )
  JVS( 1147) = W( 77 )
  JVS( 1148) = W( 78 )
  JVS( 1149) = W( 79 )
  JVS( 1150) = W( 80 )
  JVS( 1151) = W( 81 )
  JVS( 1152) = W( 82 )
  JVS( 1153) = W( 83 )
  JVS( 1154) = W( 84 )
  JVS( 1155) = W( 85 )
  JVS( 1156) = W( 86 )
  JVS( 1157) = W( 87 )
  JVS( 1158) = W( 88 )
  JVS( 1159) = W( 89 )
  JVS( 1160) = W( 90 )
  JVS( 1161) = W( 92 )
  JVS( 1162) = W( 93 )
  JVS( 1163) = W( 94 )
  JVS( 1164) = W( 95 )
  JVS( 1165) = W( 96 )
  JVS( 1166) = W( 97 )
  JVS( 1167) = W( 98 )
  JVS( 1168) = W( 99 )
  JVS( 1169) = W( 100 )
  JVS( 1170) = W( 101 )
  JVS( 1171) = W( 102 )
  JVS( 1172) = W( 103 )
  JVS( 1173) = W( 104 )
  JVS( 1174) = W( 105 )
  JVS( 1175) = W( 106 )
  JVS( 1176) = W( 107 )
  JVS( 1177) = W( 108 )
  JVS( 1178) = W( 109 )
  JVS( 1179) = W( 110 )
  JVS( 1180) = W( 111 )
  JVS( 1181) = W( 112 )
  JVS( 1182) = W( 113 )
  JVS( 1183) = W( 114 )
  JVS( 1184) = W( 115 )
  JVS( 1185) = W( 116 )
  JVS( 1186) = W( 117 )
  JVS( 1187) = W( 118 )
  JVS( 1188) = W( 119 )
  JVS( 1189) = W( 120 )
  JVS( 1190) = W( 121 )
  JVS( 1191) = W( 122 )
  JVS( 1192) = W( 123 )
  JVS( 1193) = W( 124 )
  JVS( 1194) = W( 125 )
  JVS( 1195) = W( 126 )
  JVS( 1196) = W( 127 )
  JVS( 1197) = W( 128 )
  JVS( 1198) = W( 129 )
  JVS( 1199) = W( 130 )
  JVS( 1200) = W( 131 )
  JVS( 1201) = W( 132 )
  JVS( 1202) = W( 133 )
  JVS( 1203) = W( 134 )
  JVS( 1204) = W( 135 )
  JVS( 1205) = W( 136 )
  JVS( 1206) = W( 137 )
  JVS( 1207) = W( 138 )
  JVS( 1208) = W( 139 )
  JVS( 1209) = W( 140 )
  IF ( ABS(  JVS( 1301 )) < TINY(a) ) THEN
         IER = 135                                     
         RETURN
  END IF
   W( 22 ) = JVS( 1210 )
   W( 23 ) = JVS( 1211 )
   W( 24 ) = JVS( 1212 )
   W( 25 ) = JVS( 1213 )
   W( 26 ) = JVS( 1214 )
   W( 27 ) = JVS( 1215 )
   W( 28 ) = JVS( 1216 )
   W( 30 ) = JVS( 1217 )
   W( 34 ) = JVS( 1218 )
   W( 35 ) = JVS( 1219 )
   W( 41 ) = JVS( 1220 )
   W( 42 ) = JVS( 1221 )
   W( 43 ) = JVS( 1222 )
   W( 44 ) = JVS( 1223 )
   W( 45 ) = JVS( 1224 )
   W( 46 ) = JVS( 1225 )
   W( 47 ) = JVS( 1226 )
   W( 48 ) = JVS( 1227 )
   W( 49 ) = JVS( 1228 )
   W( 51 ) = JVS( 1229 )
   W( 53 ) = JVS( 1230 )
   W( 54 ) = JVS( 1231 )
   W( 55 ) = JVS( 1232 )
   W( 56 ) = JVS( 1233 )
   W( 59 ) = JVS( 1234 )
   W( 60 ) = JVS( 1235 )
   W( 62 ) = JVS( 1236 )
   W( 64 ) = JVS( 1237 )
   W( 65 ) = JVS( 1238 )
   W( 69 ) = JVS( 1239 )
   W( 70 ) = JVS( 1240 )
   W( 72 ) = JVS( 1241 )
   W( 74 ) = JVS( 1242 )
   W( 76 ) = JVS( 1243 )
   W( 77 ) = JVS( 1244 )
   W( 78 ) = JVS( 1245 )
   W( 79 ) = JVS( 1246 )
   W( 80 ) = JVS( 1247 )
   W( 81 ) = JVS( 1248 )
   W( 82 ) = JVS( 1249 )
   W( 83 ) = JVS( 1250 )
   W( 84 ) = JVS( 1251 )
   W( 85 ) = JVS( 1252 )
   W( 86 ) = JVS( 1253 )
   W( 87 ) = JVS( 1254 )
   W( 88 ) = JVS( 1255 )
   W( 89 ) = JVS( 1256 )
   W( 90 ) = JVS( 1257 )
   W( 91 ) = JVS( 1258 )
   W( 92 ) = JVS( 1259 )
   W( 93 ) = JVS( 1260 )
   W( 94 ) = JVS( 1261 )
   W( 95 ) = JVS( 1262 )
   W( 96 ) = JVS( 1263 )
   W( 97 ) = JVS( 1264 )
   W( 98 ) = JVS( 1265 )
   W( 99 ) = JVS( 1266 )
   W( 100 ) = JVS( 1267 )
   W( 101 ) = JVS( 1268 )
   W( 102 ) = JVS( 1269 )
   W( 103 ) = JVS( 1270 )
   W( 104 ) = JVS( 1271 )
   W( 105 ) = JVS( 1272 )
   W( 106 ) = JVS( 1273 )
   W( 107 ) = JVS( 1274 )
   W( 108 ) = JVS( 1275 )
   W( 109 ) = JVS( 1276 )
   W( 111 ) = JVS( 1277 )
   W( 112 ) = JVS( 1278 )
   W( 113 ) = JVS( 1279 )
   W( 114 ) = JVS( 1280 )
   W( 115 ) = JVS( 1281 )
   W( 116 ) = JVS( 1282 )
   W( 117 ) = JVS( 1283 )
   W( 118 ) = JVS( 1284 )
   W( 119 ) = JVS( 1285 )
   W( 120 ) = JVS( 1286 )
   W( 121 ) = JVS( 1287 )
   W( 122 ) = JVS( 1288 )
   W( 123 ) = JVS( 1289 )
   W( 124 ) = JVS( 1290 )
   W( 125 ) = JVS( 1291 )
   W( 126 ) = JVS( 1292 )
   W( 127 ) = JVS( 1293 )
   W( 128 ) = JVS( 1294 )
   W( 129 ) = JVS( 1295 )
   W( 130 ) = JVS( 1296 )
   W( 131 ) = JVS( 1297 )
   W( 132 ) = JVS( 1298 )
   W( 133 ) = JVS( 1299 )
   W( 134 ) = JVS( 1300 )
   W( 135 ) = JVS( 1301 )
   W( 136 ) = JVS( 1302 )
   W( 137 ) = JVS( 1303 )
   W( 138 ) = JVS( 1304 )
   W( 139 ) = JVS( 1305 )
   W( 140 ) = JVS( 1306 )
  a = -W( 22 ) / JVS(          158  )
  W( 22 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 159 )
  W( 135 ) = W( 135 ) + a*JVS( 160 )
  a = -W( 23 ) / JVS(          161  )
  W( 23 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 162 )
  a = -W( 24 ) / JVS(          163  )
  W( 24 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 164 )
  a = -W( 25 ) / JVS(          166  )
  W( 25 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 167 )
  a = -W( 26 ) / JVS(          169  )
  W( 26 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 170 )
  a = -W( 27 ) / JVS(          171  )
  W( 27 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 172 )
  a = -W( 28 ) / JVS(          174  )
  W( 28 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 175 )
  a = -W( 30 ) / JVS(          180  )
  W( 30 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 181 )
  W( 139 ) = W( 139 ) + a*JVS( 182 )
  a = -W( 34 ) / JVS(          195  )
  W( 34 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 196 )
  W( 94 ) = W( 94 ) + a*JVS( 197 )
  W( 139 ) = W( 139 ) + a*JVS( 198 )
  a = -W( 35 ) / JVS(          200  )
  W( 35 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 201 )
  W( 134 ) = W( 134 ) + a*JVS( 202 )
  a = -W( 41 ) / JVS(          222  )
  W( 41 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 223 )
  a = -W( 42 ) / JVS(          225  )
  W( 42 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 226 )
  a = -W( 43 ) / JVS(          227  )
  W( 43 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 228 )
  W( 134 ) = W( 134 ) + a*JVS( 229 )
  a = -W( 44 ) / JVS(          230  )
  W( 44 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 231 )
  W( 139 ) = W( 139 ) + a*JVS( 232 )
  a = -W( 45 ) / JVS(          233  )
  W( 45 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 234 )
  W( 134 ) = W( 134 ) + a*JVS( 235 )
  W( 140 ) = W( 140 ) + a*JVS( 236 )
  a = -W( 46 ) / JVS(          237  )
  W( 46 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 238 )
  W( 94 ) = W( 94 ) + a*JVS( 239 )
  W( 139 ) = W( 139 ) + a*JVS( 240 )
  a = -W( 47 ) / JVS(          242  )
  W( 47 ) = -a
  W( 70 ) = W( 70 ) + a*JVS( 243 )
  W( 93 ) = W( 93 ) + a*JVS( 244 )
  W( 94 ) = W( 94 ) + a*JVS( 245 )
  W( 134 ) = W( 134 ) + a*JVS( 246 )
  W( 139 ) = W( 139 ) + a*JVS( 247 )
  a = -W( 48 ) / JVS(          248  )
  W( 48 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 249 )
  W( 134 ) = W( 134 ) + a*JVS( 250 )
  W( 135 ) = W( 135 ) + a*JVS( 251 )
  a = -W( 49 ) / JVS(          252  )
  W( 49 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 253 )
  W( 134 ) = W( 134 ) + a*JVS( 254 )
  W( 135 ) = W( 135 ) + a*JVS( 255 )
  a = -W( 51 ) / JVS(          261  )
  W( 51 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 262 )
  W( 139 ) = W( 139 ) + a*JVS( 263 )
  a = -W( 53 ) / JVS(          268  )
  W( 53 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 269 )
  W( 135 ) = W( 135 ) + a*JVS( 270 )
  W( 140 ) = W( 140 ) + a*JVS( 271 )
  a = -W( 54 ) / JVS(          272  )
  W( 54 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 273 )
  W( 134 ) = W( 134 ) + a*JVS( 274 )
  W( 135 ) = W( 135 ) + a*JVS( 275 )
  a = -W( 55 ) / JVS(          276  )
  W( 55 ) = -a
  W( 77 ) = W( 77 ) + a*JVS( 277 )
  W( 84 ) = W( 84 ) + a*JVS( 278 )
  W( 130 ) = W( 130 ) + a*JVS( 279 )
  W( 134 ) = W( 134 ) + a*JVS( 280 )
  a = -W( 56 ) / JVS(          281  )
  W( 56 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 282 )
  W( 135 ) = W( 135 ) + a*JVS( 283 )
  W( 138 ) = W( 138 ) + a*JVS( 284 )
  a = -W( 59 ) / JVS(          301  )
  W( 59 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 302 )
  W( 136 ) = W( 136 ) + a*JVS( 303 )
  a = -W( 60 ) / JVS(          305  )
  W( 60 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 306 )
  W( 136 ) = W( 136 ) + a*JVS( 307 )
  W( 139 ) = W( 139 ) + a*JVS( 308 )
  a = -W( 62 ) / JVS(          316  )
  W( 62 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 317 )
  W( 135 ) = W( 135 ) + a*JVS( 318 )
  W( 139 ) = W( 139 ) + a*JVS( 319 )
  a = -W( 64 ) / JVS(          332  )
  W( 64 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 333 )
  W( 135 ) = W( 135 ) + a*JVS( 334 )
  W( 139 ) = W( 139 ) + a*JVS( 335 )
  a = -W( 65 ) / JVS(          336  )
  W( 65 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 337 )
  W( 135 ) = W( 135 ) + a*JVS( 338 )
  W( 139 ) = W( 139 ) + a*JVS( 339 )
  a = -W( 69 ) / JVS(          359  )
  W( 69 ) = -a
  W( 119 ) = W( 119 ) + a*JVS( 360 )
  W( 134 ) = W( 134 ) + a*JVS( 361 )
  W( 135 ) = W( 135 ) + a*JVS( 362 )
  a = -W( 70 ) / JVS(          365  )
  W( 70 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 366 )
  W( 135 ) = W( 135 ) + a*JVS( 367 )
  W( 139 ) = W( 139 ) + a*JVS( 368 )
  a = -W( 72 ) / JVS(          376  )
  W( 72 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 377 )
  W( 134 ) = W( 134 ) + a*JVS( 378 )
  W( 135 ) = W( 135 ) + a*JVS( 379 )
  a = -W( 74 ) / JVS(          386  )
  W( 74 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 387 )
  W( 135 ) = W( 135 ) + a*JVS( 388 )
  W( 139 ) = W( 139 ) + a*JVS( 389 )
  a = -W( 76 ) / JVS(          399  )
  W( 76 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 400 )
  W( 134 ) = W( 134 ) + a*JVS( 401 )
  W( 135 ) = W( 135 ) + a*JVS( 402 )
  W( 139 ) = W( 139 ) + a*JVS( 403 )
  a = -W( 77 ) / JVS(          404  )
  W( 77 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 405 )
  W( 113 ) = W( 113 ) + a*JVS( 406 )
  W( 134 ) = W( 134 ) + a*JVS( 407 )
  W( 136 ) = W( 136 ) + a*JVS( 408 )
  a = -W( 78 ) / JVS(          412  )
  W( 78 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 413 )
  W( 135 ) = W( 135 ) + a*JVS( 414 )
  W( 138 ) = W( 138 ) + a*JVS( 415 )
  W( 139 ) = W( 139 ) + a*JVS( 416 )
  a = -W( 79 ) / JVS(          417  )
  W( 79 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 418 )
  W( 134 ) = W( 134 ) + a*JVS( 419 )
  W( 135 ) = W( 135 ) + a*JVS( 420 )
  a = -W( 80 ) / JVS(          421  )
  W( 80 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 422 )
  W( 134 ) = W( 134 ) + a*JVS( 423 )
  W( 135 ) = W( 135 ) + a*JVS( 424 )
  a = -W( 81 ) / JVS(          425  )
  W( 81 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 426 )
  W( 134 ) = W( 134 ) + a*JVS( 427 )
  W( 138 ) = W( 138 ) + a*JVS( 428 )
  a = -W( 82 ) / JVS(          432  )
  W( 82 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 433 )
  W( 87 ) = W( 87 ) + a*JVS( 434 )
  W( 92 ) = W( 92 ) + a*JVS( 435 )
  W( 104 ) = W( 104 ) + a*JVS( 436 )
  W( 107 ) = W( 107 ) + a*JVS( 437 )
  W( 109 ) = W( 109 ) + a*JVS( 438 )
  W( 130 ) = W( 130 ) + a*JVS( 439 )
  W( 134 ) = W( 134 ) + a*JVS( 440 )
  W( 135 ) = W( 135 ) + a*JVS( 441 )
  W( 136 ) = W( 136 ) + a*JVS( 442 )
  W( 139 ) = W( 139 ) + a*JVS( 443 )
  a = -W( 83 ) / JVS(          444  )
  W( 83 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 445 )
  W( 134 ) = W( 134 ) + a*JVS( 446 )
  W( 135 ) = W( 135 ) + a*JVS( 447 )
  a = -W( 84 ) / JVS(          451  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 452 )
  W( 130 ) = W( 130 ) + a*JVS( 453 )
  W( 134 ) = W( 134 ) + a*JVS( 454 )
  W( 136 ) = W( 136 ) + a*JVS( 455 )
  a = -W( 85 ) / JVS(          459  )
  W( 85 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 460 )
  W( 135 ) = W( 135 ) + a*JVS( 461 )
  W( 136 ) = W( 136 ) + a*JVS( 462 )
  W( 138 ) = W( 138 ) + a*JVS( 463 )
  W( 139 ) = W( 139 ) + a*JVS( 464 )
  a = -W( 86 ) / JVS(          465  )
  W( 86 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 466 )
  W( 136 ) = W( 136 ) + a*JVS( 467 )
  a = -W( 87 ) / JVS(          468  )
  W( 87 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 469 )
  W( 136 ) = W( 136 ) + a*JVS( 470 )
  a = -W( 88 ) / JVS(          471  )
  W( 88 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 472 )
  W( 134 ) = W( 134 ) + a*JVS( 473 )
  W( 135 ) = W( 135 ) + a*JVS( 474 )
  a = -W( 89 ) / JVS(          478  )
  W( 89 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 479 )
  W( 134 ) = W( 134 ) + a*JVS( 480 )
  W( 135 ) = W( 135 ) + a*JVS( 481 )
  W( 139 ) = W( 139 ) + a*JVS( 482 )
  a = -W( 90 ) / JVS(          484  )
  W( 90 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 485 )
  W( 134 ) = W( 134 ) + a*JVS( 486 )
  W( 135 ) = W( 135 ) + a*JVS( 487 )
  W( 139 ) = W( 139 ) + a*JVS( 488 )
  a = -W( 91 ) / JVS(          491  )
  W( 91 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 492 )
  W( 94 ) = W( 94 ) + a*JVS( 493 )
  W( 134 ) = W( 134 ) + a*JVS( 494 )
  W( 135 ) = W( 135 ) + a*JVS( 495 )
  W( 138 ) = W( 138 ) + a*JVS( 496 )
  W( 139 ) = W( 139 ) + a*JVS( 497 )
  a = -W( 92 ) / JVS(          498  )
  W( 92 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 499 )
  W( 134 ) = W( 134 ) + a*JVS( 500 )
  W( 136 ) = W( 136 ) + a*JVS( 501 )
  a = -W( 93 ) / JVS(          504  )
  W( 93 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 505 )
  W( 135 ) = W( 135 ) + a*JVS( 506 )
  W( 139 ) = W( 139 ) + a*JVS( 507 )
  a = -W( 94 ) / JVS(          510  )
  W( 94 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 511 )
  W( 135 ) = W( 135 ) + a*JVS( 512 )
  W( 139 ) = W( 139 ) + a*JVS( 513 )
  a = -W( 95 ) / JVS(          517  )
  W( 95 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 518 )
  W( 135 ) = W( 135 ) + a*JVS( 519 )
  W( 138 ) = W( 138 ) + a*JVS( 520 )
  W( 139 ) = W( 139 ) + a*JVS( 521 )
  a = -W( 96 ) / JVS(          523  )
  W( 96 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 524 )
  W( 130 ) = W( 130 ) + a*JVS( 525 )
  W( 134 ) = W( 134 ) + a*JVS( 526 )
  W( 135 ) = W( 135 ) + a*JVS( 527 )
  W( 136 ) = W( 136 ) + a*JVS( 528 )
  W( 138 ) = W( 138 ) + a*JVS( 529 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 533 )
  W( 135 ) = W( 135 ) + a*JVS( 534 )
  W( 139 ) = W( 139 ) + a*JVS( 535 )
  W( 140 ) = W( 140 ) + a*JVS( 536 )
  a = -W( 98 ) / JVS(          540  )
  W( 98 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 541 )
  W( 135 ) = W( 135 ) + a*JVS( 542 )
  W( 138 ) = W( 138 ) + a*JVS( 543 )
  W( 139 ) = W( 139 ) + a*JVS( 544 )
  a = -W( 99 ) / JVS(          554  )
  W( 99 ) = -a
  W( 101 ) = W( 101 ) + a*JVS( 555 )
  W( 112 ) = W( 112 ) + a*JVS( 556 )
  W( 127 ) = W( 127 ) + a*JVS( 557 )
  W( 131 ) = W( 131 ) + a*JVS( 558 )
  W( 133 ) = W( 133 ) + a*JVS( 559 )
  W( 134 ) = W( 134 ) + a*JVS( 560 )
  W( 135 ) = W( 135 ) + a*JVS( 561 )
  W( 137 ) = W( 137 ) + a*JVS( 562 )
  W( 138 ) = W( 138 ) + a*JVS( 563 )
  W( 139 ) = W( 139 ) + a*JVS( 564 )
  W( 140 ) = W( 140 ) + a*JVS( 565 )
  a = -W( 100 ) / JVS(          566  )
  W( 100 ) = -a
  W( 102 ) = W( 102 ) + a*JVS( 567 )
  W( 116 ) = W( 116 ) + a*JVS( 568 )
  W( 118 ) = W( 118 ) + a*JVS( 569 )
  W( 119 ) = W( 119 ) + a*JVS( 570 )
  W( 121 ) = W( 121 ) + a*JVS( 571 )
  W( 122 ) = W( 122 ) + a*JVS( 572 )
  W( 127 ) = W( 127 ) + a*JVS( 573 )
  W( 129 ) = W( 129 ) + a*JVS( 574 )
  W( 131 ) = W( 131 ) + a*JVS( 575 )
  W( 134 ) = W( 134 ) + a*JVS( 576 )
  W( 140 ) = W( 140 ) + a*JVS( 577 )
  a = -W( 101 ) / JVS(          582  )
  W( 101 ) = -a
  W( 104 ) = W( 104 ) + a*JVS( 583 )
  W( 107 ) = W( 107 ) + a*JVS( 584 )
  W( 109 ) = W( 109 ) + a*JVS( 585 )
  W( 134 ) = W( 134 ) + a*JVS( 586 )
  W( 135 ) = W( 135 ) + a*JVS( 587 )
  W( 136 ) = W( 136 ) + a*JVS( 588 )
  W( 139 ) = W( 139 ) + a*JVS( 589 )
  a = -W( 102 ) / JVS(          594  )
  W( 102 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 595 )
  W( 134 ) = W( 134 ) + a*JVS( 596 )
  W( 135 ) = W( 135 ) + a*JVS( 597 )
  W( 139 ) = W( 139 ) + a*JVS( 598 )
  W( 140 ) = W( 140 ) + a*JVS( 599 )
  a = -W( 103 ) / JVS(          600  )
  W( 103 ) = -a
  W( 108 ) = W( 108 ) + a*JVS( 601 )
  W( 133 ) = W( 133 ) + a*JVS( 602 )
  W( 135 ) = W( 135 ) + a*JVS( 603 )
  W( 139 ) = W( 139 ) + a*JVS( 604 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 105 ) / JVS(          610  )
  W( 105 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 611 )
  W( 134 ) = W( 134 ) + a*JVS( 612 )
  W( 135 ) = W( 135 ) + a*JVS( 613 )
  W( 136 ) = W( 136 ) + a*JVS( 614 )
  W( 139 ) = W( 139 ) + a*JVS( 615 )
  a = -W( 106 ) / JVS(          619  )
  W( 106 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 620 )
  W( 133 ) = W( 133 ) + a*JVS( 621 )
  W( 134 ) = W( 134 ) + a*JVS( 622 )
  W( 135 ) = W( 135 ) + a*JVS( 623 )
  W( 136 ) = W( 136 ) + a*JVS( 624 )
  W( 139 ) = W( 139 ) + a*JVS( 625 )
  W( 140 ) = W( 140 ) + a*JVS( 626 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  a = -W( 111 ) / JVS(          690  )
  W( 111 ) = -a
  W( 112 ) = W( 112 ) + a*JVS( 691 )
  W( 113 ) = W( 113 ) + a*JVS( 692 )
  W( 114 ) = W( 114 ) + a*JVS( 693 )
  W( 115 ) = W( 115 ) + a*JVS( 694 )
  W( 116 ) = W( 116 ) + a*JVS( 695 )
  W( 118 ) = W( 118 ) + a*JVS( 696 )
  W( 120 ) = W( 120 ) + a*JVS( 697 )
  W( 123 ) = W( 123 ) + a*JVS( 698 )
  W( 124 ) = W( 124 ) + a*JVS( 699 )
  W( 126 ) = W( 126 ) + a*JVS( 700 )
  W( 127 ) = W( 127 ) + a*JVS( 701 )
  W( 128 ) = W( 128 ) + a*JVS( 702 )
  W( 129 ) = W( 129 ) + a*JVS( 703 )
  W( 130 ) = W( 130 ) + a*JVS( 704 )
  W( 131 ) = W( 131 ) + a*JVS( 705 )
  W( 133 ) = W( 133 ) + a*JVS( 706 )
  W( 134 ) = W( 134 ) + a*JVS( 707 )
  W( 135 ) = W( 135 ) + a*JVS( 708 )
  W( 136 ) = W( 136 ) + a*JVS( 709 )
  W( 137 ) = W( 137 ) + a*JVS( 710 )
  W( 138 ) = W( 138 ) + a*JVS( 711 )
  W( 139 ) = W( 139 ) + a*JVS( 712 )
  W( 140 ) = W( 140 ) + a*JVS( 713 )
  a = -W( 112 ) / JVS(          718  )
  W( 112 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 719 )
  W( 118 ) = W( 118 ) + a*JVS( 720 )
  W( 127 ) = W( 127 ) + a*JVS( 721 )
  W( 129 ) = W( 129 ) + a*JVS( 722 )
  W( 131 ) = W( 131 ) + a*JVS( 723 )
  W( 133 ) = W( 133 ) + a*JVS( 724 )
  W( 134 ) = W( 134 ) + a*JVS( 725 )
  W( 135 ) = W( 135 ) + a*JVS( 726 )
  W( 136 ) = W( 136 ) + a*JVS( 727 )
  W( 137 ) = W( 137 ) + a*JVS( 728 )
  W( 139 ) = W( 139 ) + a*JVS( 729 )
  W( 140 ) = W( 140 ) + a*JVS( 730 )
  a = -W( 113 ) / JVS(          732  )
  W( 113 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 733 )
  W( 133 ) = W( 133 ) + a*JVS( 734 )
  W( 134 ) = W( 134 ) + a*JVS( 735 )
  W( 136 ) = W( 136 ) + a*JVS( 736 )
  a = -W( 114 ) / JVS(          748  )
  W( 114 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 749 )
  W( 128 ) = W( 128 ) + a*JVS( 750 )
  W( 133 ) = W( 133 ) + a*JVS( 751 )
  W( 134 ) = W( 134 ) + a*JVS( 752 )
  W( 135 ) = W( 135 ) + a*JVS( 753 )
  W( 136 ) = W( 136 ) + a*JVS( 754 )
  W( 139 ) = W( 139 ) + a*JVS( 755 )
  W( 140 ) = W( 140 ) + a*JVS( 756 )
  a = -W( 115 ) / JVS(          763  )
  W( 115 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 764 )
  W( 124 ) = W( 124 ) + a*JVS( 765 )
  W( 133 ) = W( 133 ) + a*JVS( 766 )
  W( 134 ) = W( 134 ) + a*JVS( 767 )
  W( 135 ) = W( 135 ) + a*JVS( 768 )
  W( 136 ) = W( 136 ) + a*JVS( 769 )
  W( 139 ) = W( 139 ) + a*JVS( 770 )
  W( 140 ) = W( 140 ) + a*JVS( 771 )
  a = -W( 116 ) / JVS(          774  )
  W( 116 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 775 )
  W( 134 ) = W( 134 ) + a*JVS( 776 )
  W( 135 ) = W( 135 ) + a*JVS( 777 )
  W( 136 ) = W( 136 ) + a*JVS( 778 )
  W( 139 ) = W( 139 ) + a*JVS( 779 )
  W( 140 ) = W( 140 ) + a*JVS( 780 )
  a = -W( 117 ) / JVS(          788  )
  W( 117 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 789 )
  W( 134 ) = W( 134 ) + a*JVS( 790 )
  W( 135 ) = W( 135 ) + a*JVS( 791 )
  W( 136 ) = W( 136 ) + a*JVS( 792 )
  W( 139 ) = W( 139 ) + a*JVS( 793 )
  a = -W( 118 ) / JVS(          795  )
  W( 118 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 796 )
  W( 133 ) = W( 133 ) + a*JVS( 797 )
  W( 134 ) = W( 134 ) + a*JVS( 798 )
  W( 135 ) = W( 135 ) + a*JVS( 799 )
  W( 139 ) = W( 139 ) + a*JVS( 800 )
  W( 140 ) = W( 140 ) + a*JVS( 801 )
  a = -W( 119 ) / JVS(          808  )
  W( 119 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 809 )
  W( 134 ) = W( 134 ) + a*JVS( 810 )
  W( 135 ) = W( 135 ) + a*JVS( 811 )
  W( 136 ) = W( 136 ) + a*JVS( 812 )
  W( 139 ) = W( 139 ) + a*JVS( 813 )
  W( 140 ) = W( 140 ) + a*JVS( 814 )
  a = -W( 120 ) / JVS(          819  )
  W( 120 ) = -a
  W( 122 ) = W( 122 ) + a*JVS( 820 )
  W( 124 ) = W( 124 ) + a*JVS( 821 )
  W( 129 ) = W( 129 ) + a*JVS( 822 )
  W( 131 ) = W( 131 ) + a*JVS( 823 )
  W( 133 ) = W( 133 ) + a*JVS( 824 )
  W( 134 ) = W( 134 ) + a*JVS( 825 )
  W( 135 ) = W( 135 ) + a*JVS( 826 )
  W( 136 ) = W( 136 ) + a*JVS( 827 )
  W( 139 ) = W( 139 ) + a*JVS( 828 )
  W( 140 ) = W( 140 ) + a*JVS( 829 )
  a = -W( 121 ) / JVS(          842  )
  W( 121 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 843 )
  W( 133 ) = W( 133 ) + a*JVS( 844 )
  W( 134 ) = W( 134 ) + a*JVS( 845 )
  W( 135 ) = W( 135 ) + a*JVS( 846 )
  W( 136 ) = W( 136 ) + a*JVS( 847 )
  W( 139 ) = W( 139 ) + a*JVS( 848 )
  W( 140 ) = W( 140 ) + a*JVS( 849 )
  a = -W( 122 ) / JVS(          853  )
  W( 122 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 854 )
  W( 133 ) = W( 133 ) + a*JVS( 855 )
  W( 134 ) = W( 134 ) + a*JVS( 856 )
  W( 135 ) = W( 135 ) + a*JVS( 857 )
  W( 136 ) = W( 136 ) + a*JVS( 858 )
  W( 139 ) = W( 139 ) + a*JVS( 859 )
  W( 140 ) = W( 140 ) + a*JVS( 860 )
  a = -W( 123 ) / JVS(          864  )
  W( 123 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 865 )
  W( 133 ) = W( 133 ) + a*JVS( 866 )
  W( 134 ) = W( 134 ) + a*JVS( 867 )
  W( 135 ) = W( 135 ) + a*JVS( 868 )
  W( 136 ) = W( 136 ) + a*JVS( 869 )
  W( 137 ) = W( 137 ) + a*JVS( 870 )
  W( 139 ) = W( 139 ) + a*JVS( 871 )
  W( 140 ) = W( 140 ) + a*JVS( 872 )
  a = -W( 124 ) / JVS(          881  )
  W( 124 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 882 )
  W( 134 ) = W( 134 ) + a*JVS( 883 )
  W( 135 ) = W( 135 ) + a*JVS( 884 )
  W( 136 ) = W( 136 ) + a*JVS( 885 )
  W( 139 ) = W( 139 ) + a*JVS( 886 )
  W( 140 ) = W( 140 ) + a*JVS( 887 )
  a = -W( 125 ) / JVS(          894  )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 895 )
  W( 128 ) = W( 128 ) + a*JVS( 896 )
  W( 129 ) = W( 129 ) + a*JVS( 897 )
  W( 131 ) = W( 131 ) + a*JVS( 898 )
  W( 132 ) = W( 132 ) + a*JVS( 899 )
  W( 133 ) = W( 133 ) + a*JVS( 900 )
  W( 134 ) = W( 134 ) + a*JVS( 901 )
  W( 135 ) = W( 135 ) + a*JVS( 902 )
  W( 136 ) = W( 136 ) + a*JVS( 903 )
  W( 137 ) = W( 137 ) + a*JVS( 904 )
  W( 138 ) = W( 138 ) + a*JVS( 905 )
  W( 139 ) = W( 139 ) + a*JVS( 906 )
  W( 140 ) = W( 140 ) + a*JVS( 907 )
  a = -W( 126 ) / JVS(          927  )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 928 )
  W( 128 ) = W( 128 ) + a*JVS( 929 )
  W( 129 ) = W( 129 ) + a*JVS( 930 )
  W( 131 ) = W( 131 ) + a*JVS( 931 )
  W( 132 ) = W( 132 ) + a*JVS( 932 )
  W( 133 ) = W( 133 ) + a*JVS( 933 )
  W( 134 ) = W( 134 ) + a*JVS( 934 )
  W( 135 ) = W( 135 ) + a*JVS( 935 )
  W( 136 ) = W( 136 ) + a*JVS( 936 )
  W( 137 ) = W( 137 ) + a*JVS( 937 )
  W( 138 ) = W( 138 ) + a*JVS( 938 )
  W( 139 ) = W( 139 ) + a*JVS( 939 )
  W( 140 ) = W( 140 ) + a*JVS( 940 )
  a = -W( 127 ) / JVS(          947  )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 948 )
  W( 131 ) = W( 131 ) + a*JVS( 949 )
  W( 133 ) = W( 133 ) + a*JVS( 950 )
  W( 134 ) = W( 134 ) + a*JVS( 951 )
  W( 135 ) = W( 135 ) + a*JVS( 952 )
  W( 136 ) = W( 136 ) + a*JVS( 953 )
  W( 137 ) = W( 137 ) + a*JVS( 954 )
  W( 139 ) = W( 139 ) + a*JVS( 955 )
  W( 140 ) = W( 140 ) + a*JVS( 956 )
  a = -W( 128 ) / JVS(          960  )
  W( 128 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 961 )
  W( 133 ) = W( 133 ) + a*JVS( 962 )
  W( 134 ) = W( 134 ) + a*JVS( 963 )
  W( 135 ) = W( 135 ) + a*JVS( 964 )
  W( 136 ) = W( 136 ) + a*JVS( 965 )
  W( 137 ) = W( 137 ) + a*JVS( 966 )
  W( 139 ) = W( 139 ) + a*JVS( 967 )
  W( 140 ) = W( 140 ) + a*JVS( 968 )
  a = -W( 129 ) / JVS(          972  )
  W( 129 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 973 )
  W( 133 ) = W( 133 ) + a*JVS( 974 )
  W( 134 ) = W( 134 ) + a*JVS( 975 )
  W( 135 ) = W( 135 ) + a*JVS( 976 )
  W( 136 ) = W( 136 ) + a*JVS( 977 )
  W( 137 ) = W( 137 ) + a*JVS( 978 )
  W( 139 ) = W( 139 ) + a*JVS( 979 )
  W( 140 ) = W( 140 ) + a*JVS( 980 )
  a = -W( 130 ) / JVS(         1028  )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 1029 )
  W( 132 ) = W( 132 ) + a*JVS( 1030 )
  W( 133 ) = W( 133 ) + a*JVS( 1031 )
  W( 134 ) = W( 134 ) + a*JVS( 1032 )
  W( 135 ) = W( 135 ) + a*JVS( 1033 )
  W( 136 ) = W( 136 ) + a*JVS( 1034 )
  W( 137 ) = W( 137 ) + a*JVS( 1035 )
  W( 138 ) = W( 138 ) + a*JVS( 1036 )
  W( 139 ) = W( 139 ) + a*JVS( 1037 )
  W( 140 ) = W( 140 ) + a*JVS( 1038 )
  a = -W( 131 ) / JVS(         1041  )
  W( 131 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1042 )
  W( 134 ) = W( 134 ) + a*JVS( 1043 )
  W( 135 ) = W( 135 ) + a*JVS( 1044 )
  W( 136 ) = W( 136 ) + a*JVS( 1045 )
  W( 137 ) = W( 137 ) + a*JVS( 1046 )
  W( 139 ) = W( 139 ) + a*JVS( 1047 )
  W( 140 ) = W( 140 ) + a*JVS( 1048 )
  a = -W( 132 ) / JVS(         1055  )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1056 )
  W( 134 ) = W( 134 ) + a*JVS( 1057 )
  W( 135 ) = W( 135 ) + a*JVS( 1058 )
  W( 136 ) = W( 136 ) + a*JVS( 1059 )
  W( 137 ) = W( 137 ) + a*JVS( 1060 )
  W( 138 ) = W( 138 ) + a*JVS( 1061 )
  W( 139 ) = W( 139 ) + a*JVS( 1062 )
  W( 140 ) = W( 140 ) + a*JVS( 1063 )
  a = -W( 133 ) / JVS(         1090  )
  W( 133 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 1091 )
  W( 135 ) = W( 135 ) + a*JVS( 1092 )
  W( 136 ) = W( 136 ) + a*JVS( 1093 )
  W( 137 ) = W( 137 ) + a*JVS( 1094 )
  W( 138 ) = W( 138 ) + a*JVS( 1095 )
  W( 139 ) = W( 139 ) + a*JVS( 1096 )
  W( 140 ) = W( 140 ) + a*JVS( 1097 )
  a = -W( 134 ) / JVS(         1203  )
  W( 134 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 1204 )
  W( 136 ) = W( 136 ) + a*JVS( 1205 )
  W( 137 ) = W( 137 ) + a*JVS( 1206 )
  W( 138 ) = W( 138 ) + a*JVS( 1207 )
  W( 139 ) = W( 139 ) + a*JVS( 1208 )
  W( 140 ) = W( 140 ) + a*JVS( 1209 )
  JVS( 1210) = W( 22 )
  JVS( 1211) = W( 23 )
  JVS( 1212) = W( 24 )
  JVS( 1213) = W( 25 )
  JVS( 1214) = W( 26 )
  JVS( 1215) = W( 27 )
  JVS( 1216) = W( 28 )
  JVS( 1217) = W( 30 )
  JVS( 1218) = W( 34 )
  JVS( 1219) = W( 35 )
  JVS( 1220) = W( 41 )
  JVS( 1221) = W( 42 )
  JVS( 1222) = W( 43 )
  JVS( 1223) = W( 44 )
  JVS( 1224) = W( 45 )
  JVS( 1225) = W( 46 )
  JVS( 1226) = W( 47 )
  JVS( 1227) = W( 48 )
  JVS( 1228) = W( 49 )
  JVS( 1229) = W( 51 )
  JVS( 1230) = W( 53 )
  JVS( 1231) = W( 54 )
  JVS( 1232) = W( 55 )
  JVS( 1233) = W( 56 )
  JVS( 1234) = W( 59 )
  JVS( 1235) = W( 60 )
  JVS( 1236) = W( 62 )
  JVS( 1237) = W( 64 )
  JVS( 1238) = W( 65 )
  JVS( 1239) = W( 69 )
  JVS( 1240) = W( 70 )
  JVS( 1241) = W( 72 )
  JVS( 1242) = W( 74 )
  JVS( 1243) = W( 76 )
  JVS( 1244) = W( 77 )
  JVS( 1245) = W( 78 )
  JVS( 1246) = W( 79 )
  JVS( 1247) = W( 80 )
  JVS( 1248) = W( 81 )
  JVS( 1249) = W( 82 )
  JVS( 1250) = W( 83 )
  JVS( 1251) = W( 84 )
  JVS( 1252) = W( 85 )
  JVS( 1253) = W( 86 )
  JVS( 1254) = W( 87 )
  JVS( 1255) = W( 88 )
  JVS( 1256) = W( 89 )
  JVS( 1257) = W( 90 )
  JVS( 1258) = W( 91 )
  JVS( 1259) = W( 92 )
  JVS( 1260) = W( 93 )
  JVS( 1261) = W( 94 )
  JVS( 1262) = W( 95 )
  JVS( 1263) = W( 96 )
  JVS( 1264) = W( 97 )
  JVS( 1265) = W( 98 )
  JVS( 1266) = W( 99 )
  JVS( 1267) = W( 100 )
  JVS( 1268) = W( 101 )
  JVS( 1269) = W( 102 )
  JVS( 1270) = W( 103 )
  JVS( 1271) = W( 104 )
  JVS( 1272) = W( 105 )
  JVS( 1273) = W( 106 )
  JVS( 1274) = W( 107 )
  JVS( 1275) = W( 108 )
  JVS( 1276) = W( 109 )
  JVS( 1277) = W( 111 )
  JVS( 1278) = W( 112 )
  JVS( 1279) = W( 113 )
  JVS( 1280) = W( 114 )
  JVS( 1281) = W( 115 )
  JVS( 1282) = W( 116 )
  JVS( 1283) = W( 117 )
  JVS( 1284) = W( 118 )
  JVS( 1285) = W( 119 )
  JVS( 1286) = W( 120 )
  JVS( 1287) = W( 121 )
  JVS( 1288) = W( 122 )
  JVS( 1289) = W( 123 )
  JVS( 1290) = W( 124 )
  JVS( 1291) = W( 125 )
  JVS( 1292) = W( 126 )
  JVS( 1293) = W( 127 )
  JVS( 1294) = W( 128 )
  JVS( 1295) = W( 129 )
  JVS( 1296) = W( 130 )
  JVS( 1297) = W( 131 )
  JVS( 1298) = W( 132 )
  JVS( 1299) = W( 133 )
  JVS( 1300) = W( 134 )
  JVS( 1301) = W( 135 )
  JVS( 1302) = W( 136 )
  JVS( 1303) = W( 137 )
  JVS( 1304) = W( 138 )
  JVS( 1305) = W( 139 )
  JVS( 1306) = W( 140 )
  IF ( ABS(  JVS( 1327 )) < TINY(a) ) THEN
         IER = 136                                     
         RETURN
  END IF
   W( 59 ) = JVS( 1307 )
   W( 63 ) = JVS( 1308 )
   W( 85 ) = JVS( 1309 )
   W( 86 ) = JVS( 1310 )
   W( 87 ) = JVS( 1311 )
   W( 92 ) = JVS( 1312 )
   W( 96 ) = JVS( 1313 )
   W( 104 ) = JVS( 1314 )
   W( 107 ) = JVS( 1315 )
   W( 108 ) = JVS( 1316 )
   W( 109 ) = JVS( 1317 )
   W( 113 ) = JVS( 1318 )
   W( 123 ) = JVS( 1319 )
   W( 128 ) = JVS( 1320 )
   W( 130 ) = JVS( 1321 )
   W( 131 ) = JVS( 1322 )
   W( 132 ) = JVS( 1323 )
   W( 133 ) = JVS( 1324 )
   W( 134 ) = JVS( 1325 )
   W( 135 ) = JVS( 1326 )
   W( 136 ) = JVS( 1327 )
   W( 137 ) = JVS( 1328 )
   W( 138 ) = JVS( 1329 )
   W( 139 ) = JVS( 1330 )
   W( 140 ) = JVS( 1331 )
  a = -W( 59 ) / JVS(          301  )
  W( 59 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 302 )
  W( 136 ) = W( 136 ) + a*JVS( 303 )
  a = -W( 63 ) / JVS(          323  )
  W( 63 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 324 )
  W( 134 ) = W( 134 ) + a*JVS( 325 )
  W( 136 ) = W( 136 ) + a*JVS( 326 )
  W( 138 ) = W( 138 ) + a*JVS( 327 )
  W( 139 ) = W( 139 ) + a*JVS( 328 )
  a = -W( 85 ) / JVS(          459  )
  W( 85 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 460 )
  W( 135 ) = W( 135 ) + a*JVS( 461 )
  W( 136 ) = W( 136 ) + a*JVS( 462 )
  W( 138 ) = W( 138 ) + a*JVS( 463 )
  W( 139 ) = W( 139 ) + a*JVS( 464 )
  a = -W( 86 ) / JVS(          465  )
  W( 86 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 466 )
  W( 136 ) = W( 136 ) + a*JVS( 467 )
  a = -W( 87 ) / JVS(          468  )
  W( 87 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 469 )
  W( 136 ) = W( 136 ) + a*JVS( 470 )
  a = -W( 92 ) / JVS(          498  )
  W( 92 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 499 )
  W( 134 ) = W( 134 ) + a*JVS( 500 )
  W( 136 ) = W( 136 ) + a*JVS( 501 )
  a = -W( 96 ) / JVS(          523  )
  W( 96 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 524 )
  W( 130 ) = W( 130 ) + a*JVS( 525 )
  W( 134 ) = W( 134 ) + a*JVS( 526 )
  W( 135 ) = W( 135 ) + a*JVS( 527 )
  W( 136 ) = W( 136 ) + a*JVS( 528 )
  W( 138 ) = W( 138 ) + a*JVS( 529 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  a = -W( 113 ) / JVS(          732  )
  W( 113 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 733 )
  W( 133 ) = W( 133 ) + a*JVS( 734 )
  W( 134 ) = W( 134 ) + a*JVS( 735 )
  W( 136 ) = W( 136 ) + a*JVS( 736 )
  a = -W( 123 ) / JVS(          864  )
  W( 123 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 865 )
  W( 133 ) = W( 133 ) + a*JVS( 866 )
  W( 134 ) = W( 134 ) + a*JVS( 867 )
  W( 135 ) = W( 135 ) + a*JVS( 868 )
  W( 136 ) = W( 136 ) + a*JVS( 869 )
  W( 137 ) = W( 137 ) + a*JVS( 870 )
  W( 139 ) = W( 139 ) + a*JVS( 871 )
  W( 140 ) = W( 140 ) + a*JVS( 872 )
  a = -W( 128 ) / JVS(          960  )
  W( 128 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 961 )
  W( 133 ) = W( 133 ) + a*JVS( 962 )
  W( 134 ) = W( 134 ) + a*JVS( 963 )
  W( 135 ) = W( 135 ) + a*JVS( 964 )
  W( 136 ) = W( 136 ) + a*JVS( 965 )
  W( 137 ) = W( 137 ) + a*JVS( 966 )
  W( 139 ) = W( 139 ) + a*JVS( 967 )
  W( 140 ) = W( 140 ) + a*JVS( 968 )
  a = -W( 130 ) / JVS(         1028  )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 1029 )
  W( 132 ) = W( 132 ) + a*JVS( 1030 )
  W( 133 ) = W( 133 ) + a*JVS( 1031 )
  W( 134 ) = W( 134 ) + a*JVS( 1032 )
  W( 135 ) = W( 135 ) + a*JVS( 1033 )
  W( 136 ) = W( 136 ) + a*JVS( 1034 )
  W( 137 ) = W( 137 ) + a*JVS( 1035 )
  W( 138 ) = W( 138 ) + a*JVS( 1036 )
  W( 139 ) = W( 139 ) + a*JVS( 1037 )
  W( 140 ) = W( 140 ) + a*JVS( 1038 )
  a = -W( 131 ) / JVS(         1041  )
  W( 131 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1042 )
  W( 134 ) = W( 134 ) + a*JVS( 1043 )
  W( 135 ) = W( 135 ) + a*JVS( 1044 )
  W( 136 ) = W( 136 ) + a*JVS( 1045 )
  W( 137 ) = W( 137 ) + a*JVS( 1046 )
  W( 139 ) = W( 139 ) + a*JVS( 1047 )
  W( 140 ) = W( 140 ) + a*JVS( 1048 )
  a = -W( 132 ) / JVS(         1055  )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1056 )
  W( 134 ) = W( 134 ) + a*JVS( 1057 )
  W( 135 ) = W( 135 ) + a*JVS( 1058 )
  W( 136 ) = W( 136 ) + a*JVS( 1059 )
  W( 137 ) = W( 137 ) + a*JVS( 1060 )
  W( 138 ) = W( 138 ) + a*JVS( 1061 )
  W( 139 ) = W( 139 ) + a*JVS( 1062 )
  W( 140 ) = W( 140 ) + a*JVS( 1063 )
  a = -W( 133 ) / JVS(         1090  )
  W( 133 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 1091 )
  W( 135 ) = W( 135 ) + a*JVS( 1092 )
  W( 136 ) = W( 136 ) + a*JVS( 1093 )
  W( 137 ) = W( 137 ) + a*JVS( 1094 )
  W( 138 ) = W( 138 ) + a*JVS( 1095 )
  W( 139 ) = W( 139 ) + a*JVS( 1096 )
  W( 140 ) = W( 140 ) + a*JVS( 1097 )
  a = -W( 134 ) / JVS(         1203  )
  W( 134 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 1204 )
  W( 136 ) = W( 136 ) + a*JVS( 1205 )
  W( 137 ) = W( 137 ) + a*JVS( 1206 )
  W( 138 ) = W( 138 ) + a*JVS( 1207 )
  W( 139 ) = W( 139 ) + a*JVS( 1208 )
  W( 140 ) = W( 140 ) + a*JVS( 1209 )
  a = -W( 135 ) / JVS(         1301  )
  W( 135 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 1302 )
  W( 137 ) = W( 137 ) + a*JVS( 1303 )
  W( 138 ) = W( 138 ) + a*JVS( 1304 )
  W( 139 ) = W( 139 ) + a*JVS( 1305 )
  W( 140 ) = W( 140 ) + a*JVS( 1306 )
  JVS( 1307) = W( 59 )
  JVS( 1308) = W( 63 )
  JVS( 1309) = W( 85 )
  JVS( 1310) = W( 86 )
  JVS( 1311) = W( 87 )
  JVS( 1312) = W( 92 )
  JVS( 1313) = W( 96 )
  JVS( 1314) = W( 104 )
  JVS( 1315) = W( 107 )
  JVS( 1316) = W( 108 )
  JVS( 1317) = W( 109 )
  JVS( 1318) = W( 113 )
  JVS( 1319) = W( 123 )
  JVS( 1320) = W( 128 )
  JVS( 1321) = W( 130 )
  JVS( 1322) = W( 131 )
  JVS( 1323) = W( 132 )
  JVS( 1324) = W( 133 )
  JVS( 1325) = W( 134 )
  JVS( 1326) = W( 135 )
  JVS( 1327) = W( 136 )
  JVS( 1328) = W( 137 )
  JVS( 1329) = W( 138 )
  JVS( 1330) = W( 139 )
  JVS( 1331) = W( 140 )
  IF ( ABS(  JVS( 1370 )) < TINY(a) ) THEN
         IER = 137                                     
         RETURN
  END IF
   W( 35 ) = JVS( 1332 )
   W( 41 ) = JVS( 1333 )
   W( 44 ) = JVS( 1334 )
   W( 50 ) = JVS( 1335 )
   W( 52 ) = JVS( 1336 )
   W( 66 ) = JVS( 1337 )
   W( 68 ) = JVS( 1338 )
   W( 73 ) = JVS( 1339 )
   W( 86 ) = JVS( 1340 )
   W( 87 ) = JVS( 1341 )
   W( 88 ) = JVS( 1342 )
   W( 89 ) = JVS( 1343 )
   W( 94 ) = JVS( 1344 )
   W( 101 ) = JVS( 1345 )
   W( 104 ) = JVS( 1346 )
   W( 107 ) = JVS( 1347 )
   W( 109 ) = JVS( 1348 )
   W( 110 ) = JVS( 1349 )
   W( 114 ) = JVS( 1350 )
   W( 115 ) = JVS( 1351 )
   W( 116 ) = JVS( 1352 )
   W( 117 ) = JVS( 1353 )
   W( 118 ) = JVS( 1354 )
   W( 119 ) = JVS( 1355 )
   W( 121 ) = JVS( 1356 )
   W( 123 ) = JVS( 1357 )
   W( 124 ) = JVS( 1358 )
   W( 125 ) = JVS( 1359 )
   W( 126 ) = JVS( 1360 )
   W( 127 ) = JVS( 1361 )
   W( 128 ) = JVS( 1362 )
   W( 129 ) = JVS( 1363 )
   W( 131 ) = JVS( 1364 )
   W( 132 ) = JVS( 1365 )
   W( 133 ) = JVS( 1366 )
   W( 134 ) = JVS( 1367 )
   W( 135 ) = JVS( 1368 )
   W( 136 ) = JVS( 1369 )
   W( 137 ) = JVS( 1370 )
   W( 138 ) = JVS( 1371 )
   W( 139 ) = JVS( 1372 )
   W( 140 ) = JVS( 1373 )
  a = -W( 35 ) / JVS(          200  )
  W( 35 ) = -a
  W( 41 ) = W( 41 ) + a*JVS( 201 )
  W( 134 ) = W( 134 ) + a*JVS( 202 )
  a = -W( 41 ) / JVS(          222  )
  W( 41 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 223 )
  a = -W( 44 ) / JVS(          230  )
  W( 44 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 231 )
  W( 139 ) = W( 139 ) + a*JVS( 232 )
  a = -W( 50 ) / JVS(          256  )
  W( 50 ) = -a
  W( 121 ) = W( 121 ) + a*JVS( 257 )
  W( 134 ) = W( 134 ) + a*JVS( 258 )
  W( 135 ) = W( 135 ) + a*JVS( 259 )
  a = -W( 52 ) / JVS(          264  )
  W( 52 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 265 )
  W( 134 ) = W( 134 ) + a*JVS( 266 )
  W( 135 ) = W( 135 ) + a*JVS( 267 )
  a = -W( 66 ) / JVS(          340  )
  W( 66 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 341 )
  W( 117 ) = W( 117 ) + a*JVS( 342 )
  W( 134 ) = W( 134 ) + a*JVS( 343 )
  W( 139 ) = W( 139 ) + a*JVS( 344 )
  a = -W( 68 ) / JVS(          354  )
  W( 68 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 355 )
  W( 134 ) = W( 134 ) + a*JVS( 356 )
  W( 135 ) = W( 135 ) + a*JVS( 357 )
  W( 137 ) = W( 137 ) + a*JVS( 358 )
  a = -W( 73 ) / JVS(          380  )
  W( 73 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 381 )
  W( 137 ) = W( 137 ) + a*JVS( 382 )
  W( 138 ) = W( 138 ) + a*JVS( 383 )
  a = -W( 86 ) / JVS(          465  )
  W( 86 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 466 )
  W( 136 ) = W( 136 ) + a*JVS( 467 )
  a = -W( 87 ) / JVS(          468  )
  W( 87 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 469 )
  W( 136 ) = W( 136 ) + a*JVS( 470 )
  a = -W( 88 ) / JVS(          471  )
  W( 88 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 472 )
  W( 134 ) = W( 134 ) + a*JVS( 473 )
  W( 135 ) = W( 135 ) + a*JVS( 474 )
  a = -W( 89 ) / JVS(          478  )
  W( 89 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 479 )
  W( 134 ) = W( 134 ) + a*JVS( 480 )
  W( 135 ) = W( 135 ) + a*JVS( 481 )
  W( 139 ) = W( 139 ) + a*JVS( 482 )
  a = -W( 94 ) / JVS(          510  )
  W( 94 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 511 )
  W( 135 ) = W( 135 ) + a*JVS( 512 )
  W( 139 ) = W( 139 ) + a*JVS( 513 )
  a = -W( 101 ) / JVS(          582  )
  W( 101 ) = -a
  W( 104 ) = W( 104 ) + a*JVS( 583 )
  W( 107 ) = W( 107 ) + a*JVS( 584 )
  W( 109 ) = W( 109 ) + a*JVS( 585 )
  W( 134 ) = W( 134 ) + a*JVS( 586 )
  W( 135 ) = W( 135 ) + a*JVS( 587 )
  W( 136 ) = W( 136 ) + a*JVS( 588 )
  W( 139 ) = W( 139 ) + a*JVS( 589 )
  a = -W( 104 ) / JVS(          605  )
  W( 104 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 606 )
  W( 134 ) = W( 134 ) + a*JVS( 607 )
  W( 136 ) = W( 136 ) + a*JVS( 608 )
  a = -W( 107 ) / JVS(          627  )
  W( 107 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 628 )
  W( 134 ) = W( 134 ) + a*JVS( 629 )
  W( 136 ) = W( 136 ) + a*JVS( 630 )
  a = -W( 109 ) / JVS(          635  )
  W( 109 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 636 )
  W( 134 ) = W( 134 ) + a*JVS( 637 )
  W( 136 ) = W( 136 ) + a*JVS( 638 )
  a = -W( 110 ) / JVS(          655  )
  W( 110 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 656 )
  W( 116 ) = W( 116 ) + a*JVS( 657 )
  W( 117 ) = W( 117 ) + a*JVS( 658 )
  W( 118 ) = W( 118 ) + a*JVS( 659 )
  W( 119 ) = W( 119 ) + a*JVS( 660 )
  W( 133 ) = W( 133 ) + a*JVS( 661 )
  W( 134 ) = W( 134 ) + a*JVS( 662 )
  W( 135 ) = W( 135 ) + a*JVS( 663 )
  W( 136 ) = W( 136 ) + a*JVS( 664 )
  W( 139 ) = W( 139 ) + a*JVS( 665 )
  W( 140 ) = W( 140 ) + a*JVS( 666 )
  a = -W( 114 ) / JVS(          748  )
  W( 114 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 749 )
  W( 128 ) = W( 128 ) + a*JVS( 750 )
  W( 133 ) = W( 133 ) + a*JVS( 751 )
  W( 134 ) = W( 134 ) + a*JVS( 752 )
  W( 135 ) = W( 135 ) + a*JVS( 753 )
  W( 136 ) = W( 136 ) + a*JVS( 754 )
  W( 139 ) = W( 139 ) + a*JVS( 755 )
  W( 140 ) = W( 140 ) + a*JVS( 756 )
  a = -W( 115 ) / JVS(          763  )
  W( 115 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 764 )
  W( 124 ) = W( 124 ) + a*JVS( 765 )
  W( 133 ) = W( 133 ) + a*JVS( 766 )
  W( 134 ) = W( 134 ) + a*JVS( 767 )
  W( 135 ) = W( 135 ) + a*JVS( 768 )
  W( 136 ) = W( 136 ) + a*JVS( 769 )
  W( 139 ) = W( 139 ) + a*JVS( 770 )
  W( 140 ) = W( 140 ) + a*JVS( 771 )
  a = -W( 116 ) / JVS(          774  )
  W( 116 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 775 )
  W( 134 ) = W( 134 ) + a*JVS( 776 )
  W( 135 ) = W( 135 ) + a*JVS( 777 )
  W( 136 ) = W( 136 ) + a*JVS( 778 )
  W( 139 ) = W( 139 ) + a*JVS( 779 )
  W( 140 ) = W( 140 ) + a*JVS( 780 )
  a = -W( 117 ) / JVS(          788  )
  W( 117 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 789 )
  W( 134 ) = W( 134 ) + a*JVS( 790 )
  W( 135 ) = W( 135 ) + a*JVS( 791 )
  W( 136 ) = W( 136 ) + a*JVS( 792 )
  W( 139 ) = W( 139 ) + a*JVS( 793 )
  a = -W( 118 ) / JVS(          795  )
  W( 118 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 796 )
  W( 133 ) = W( 133 ) + a*JVS( 797 )
  W( 134 ) = W( 134 ) + a*JVS( 798 )
  W( 135 ) = W( 135 ) + a*JVS( 799 )
  W( 139 ) = W( 139 ) + a*JVS( 800 )
  W( 140 ) = W( 140 ) + a*JVS( 801 )
  a = -W( 119 ) / JVS(          808  )
  W( 119 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 809 )
  W( 134 ) = W( 134 ) + a*JVS( 810 )
  W( 135 ) = W( 135 ) + a*JVS( 811 )
  W( 136 ) = W( 136 ) + a*JVS( 812 )
  W( 139 ) = W( 139 ) + a*JVS( 813 )
  W( 140 ) = W( 140 ) + a*JVS( 814 )
  a = -W( 121 ) / JVS(          842  )
  W( 121 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 843 )
  W( 133 ) = W( 133 ) + a*JVS( 844 )
  W( 134 ) = W( 134 ) + a*JVS( 845 )
  W( 135 ) = W( 135 ) + a*JVS( 846 )
  W( 136 ) = W( 136 ) + a*JVS( 847 )
  W( 139 ) = W( 139 ) + a*JVS( 848 )
  W( 140 ) = W( 140 ) + a*JVS( 849 )
  a = -W( 123 ) / JVS(          864  )
  W( 123 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 865 )
  W( 133 ) = W( 133 ) + a*JVS( 866 )
  W( 134 ) = W( 134 ) + a*JVS( 867 )
  W( 135 ) = W( 135 ) + a*JVS( 868 )
  W( 136 ) = W( 136 ) + a*JVS( 869 )
  W( 137 ) = W( 137 ) + a*JVS( 870 )
  W( 139 ) = W( 139 ) + a*JVS( 871 )
  W( 140 ) = W( 140 ) + a*JVS( 872 )
  a = -W( 124 ) / JVS(          881  )
  W( 124 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 882 )
  W( 134 ) = W( 134 ) + a*JVS( 883 )
  W( 135 ) = W( 135 ) + a*JVS( 884 )
  W( 136 ) = W( 136 ) + a*JVS( 885 )
  W( 139 ) = W( 139 ) + a*JVS( 886 )
  W( 140 ) = W( 140 ) + a*JVS( 887 )
  a = -W( 125 ) / JVS(          894  )
  W( 125 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 895 )
  W( 128 ) = W( 128 ) + a*JVS( 896 )
  W( 129 ) = W( 129 ) + a*JVS( 897 )
  W( 131 ) = W( 131 ) + a*JVS( 898 )
  W( 132 ) = W( 132 ) + a*JVS( 899 )
  W( 133 ) = W( 133 ) + a*JVS( 900 )
  W( 134 ) = W( 134 ) + a*JVS( 901 )
  W( 135 ) = W( 135 ) + a*JVS( 902 )
  W( 136 ) = W( 136 ) + a*JVS( 903 )
  W( 137 ) = W( 137 ) + a*JVS( 904 )
  W( 138 ) = W( 138 ) + a*JVS( 905 )
  W( 139 ) = W( 139 ) + a*JVS( 906 )
  W( 140 ) = W( 140 ) + a*JVS( 907 )
  a = -W( 126 ) / JVS(          927  )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 928 )
  W( 128 ) = W( 128 ) + a*JVS( 929 )
  W( 129 ) = W( 129 ) + a*JVS( 930 )
  W( 131 ) = W( 131 ) + a*JVS( 931 )
  W( 132 ) = W( 132 ) + a*JVS( 932 )
  W( 133 ) = W( 133 ) + a*JVS( 933 )
  W( 134 ) = W( 134 ) + a*JVS( 934 )
  W( 135 ) = W( 135 ) + a*JVS( 935 )
  W( 136 ) = W( 136 ) + a*JVS( 936 )
  W( 137 ) = W( 137 ) + a*JVS( 937 )
  W( 138 ) = W( 138 ) + a*JVS( 938 )
  W( 139 ) = W( 139 ) + a*JVS( 939 )
  W( 140 ) = W( 140 ) + a*JVS( 940 )
  a = -W( 127 ) / JVS(          947  )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 948 )
  W( 131 ) = W( 131 ) + a*JVS( 949 )
  W( 133 ) = W( 133 ) + a*JVS( 950 )
  W( 134 ) = W( 134 ) + a*JVS( 951 )
  W( 135 ) = W( 135 ) + a*JVS( 952 )
  W( 136 ) = W( 136 ) + a*JVS( 953 )
  W( 137 ) = W( 137 ) + a*JVS( 954 )
  W( 139 ) = W( 139 ) + a*JVS( 955 )
  W( 140 ) = W( 140 ) + a*JVS( 956 )
  a = -W( 128 ) / JVS(          960  )
  W( 128 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 961 )
  W( 133 ) = W( 133 ) + a*JVS( 962 )
  W( 134 ) = W( 134 ) + a*JVS( 963 )
  W( 135 ) = W( 135 ) + a*JVS( 964 )
  W( 136 ) = W( 136 ) + a*JVS( 965 )
  W( 137 ) = W( 137 ) + a*JVS( 966 )
  W( 139 ) = W( 139 ) + a*JVS( 967 )
  W( 140 ) = W( 140 ) + a*JVS( 968 )
  a = -W( 129 ) / JVS(          972  )
  W( 129 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 973 )
  W( 133 ) = W( 133 ) + a*JVS( 974 )
  W( 134 ) = W( 134 ) + a*JVS( 975 )
  W( 135 ) = W( 135 ) + a*JVS( 976 )
  W( 136 ) = W( 136 ) + a*JVS( 977 )
  W( 137 ) = W( 137 ) + a*JVS( 978 )
  W( 139 ) = W( 139 ) + a*JVS( 979 )
  W( 140 ) = W( 140 ) + a*JVS( 980 )
  a = -W( 131 ) / JVS(         1041  )
  W( 131 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1042 )
  W( 134 ) = W( 134 ) + a*JVS( 1043 )
  W( 135 ) = W( 135 ) + a*JVS( 1044 )
  W( 136 ) = W( 136 ) + a*JVS( 1045 )
  W( 137 ) = W( 137 ) + a*JVS( 1046 )
  W( 139 ) = W( 139 ) + a*JVS( 1047 )
  W( 140 ) = W( 140 ) + a*JVS( 1048 )
  a = -W( 132 ) / JVS(         1055  )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1056 )
  W( 134 ) = W( 134 ) + a*JVS( 1057 )
  W( 135 ) = W( 135 ) + a*JVS( 1058 )
  W( 136 ) = W( 136 ) + a*JVS( 1059 )
  W( 137 ) = W( 137 ) + a*JVS( 1060 )
  W( 138 ) = W( 138 ) + a*JVS( 1061 )
  W( 139 ) = W( 139 ) + a*JVS( 1062 )
  W( 140 ) = W( 140 ) + a*JVS( 1063 )
  a = -W( 133 ) / JVS(         1090  )
  W( 133 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 1091 )
  W( 135 ) = W( 135 ) + a*JVS( 1092 )
  W( 136 ) = W( 136 ) + a*JVS( 1093 )
  W( 137 ) = W( 137 ) + a*JVS( 1094 )
  W( 138 ) = W( 138 ) + a*JVS( 1095 )
  W( 139 ) = W( 139 ) + a*JVS( 1096 )
  W( 140 ) = W( 140 ) + a*JVS( 1097 )
  a = -W( 134 ) / JVS(         1203  )
  W( 134 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 1204 )
  W( 136 ) = W( 136 ) + a*JVS( 1205 )
  W( 137 ) = W( 137 ) + a*JVS( 1206 )
  W( 138 ) = W( 138 ) + a*JVS( 1207 )
  W( 139 ) = W( 139 ) + a*JVS( 1208 )
  W( 140 ) = W( 140 ) + a*JVS( 1209 )
  a = -W( 135 ) / JVS(         1301  )
  W( 135 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 1302 )
  W( 137 ) = W( 137 ) + a*JVS( 1303 )
  W( 138 ) = W( 138 ) + a*JVS( 1304 )
  W( 139 ) = W( 139 ) + a*JVS( 1305 )
  W( 140 ) = W( 140 ) + a*JVS( 1306 )
  a = -W( 136 ) / JVS(         1327  )
  W( 136 ) = -a
  W( 137 ) = W( 137 ) + a*JVS( 1328 )
  W( 138 ) = W( 138 ) + a*JVS( 1329 )
  W( 139 ) = W( 139 ) + a*JVS( 1330 )
  W( 140 ) = W( 140 ) + a*JVS( 1331 )
  JVS( 1332) = W( 35 )
  JVS( 1333) = W( 41 )
  JVS( 1334) = W( 44 )
  JVS( 1335) = W( 50 )
  JVS( 1336) = W( 52 )
  JVS( 1337) = W( 66 )
  JVS( 1338) = W( 68 )
  JVS( 1339) = W( 73 )
  JVS( 1340) = W( 86 )
  JVS( 1341) = W( 87 )
  JVS( 1342) = W( 88 )
  JVS( 1343) = W( 89 )
  JVS( 1344) = W( 94 )
  JVS( 1345) = W( 101 )
  JVS( 1346) = W( 104 )
  JVS( 1347) = W( 107 )
  JVS( 1348) = W( 109 )
  JVS( 1349) = W( 110 )
  JVS( 1350) = W( 114 )
  JVS( 1351) = W( 115 )
  JVS( 1352) = W( 116 )
  JVS( 1353) = W( 117 )
  JVS( 1354) = W( 118 )
  JVS( 1355) = W( 119 )
  JVS( 1356) = W( 121 )
  JVS( 1357) = W( 123 )
  JVS( 1358) = W( 124 )
  JVS( 1359) = W( 125 )
  JVS( 1360) = W( 126 )
  JVS( 1361) = W( 127 )
  JVS( 1362) = W( 128 )
  JVS( 1363) = W( 129 )
  JVS( 1364) = W( 131 )
  JVS( 1365) = W( 132 )
  JVS( 1366) = W( 133 )
  JVS( 1367) = W( 134 )
  JVS( 1368) = W( 135 )
  JVS( 1369) = W( 136 )
  JVS( 1370) = W( 137 )
  JVS( 1371) = W( 138 )
  JVS( 1372) = W( 139 )
  JVS( 1373) = W( 140 )
  IF ( ABS(  JVS( 1428 )) < TINY(a) ) THEN
         IER = 138                                     
         RETURN
  END IF
   W( 21 ) = JVS( 1374 )
   W( 40 ) = JVS( 1375 )
   W( 51 ) = JVS( 1376 )
   W( 56 ) = JVS( 1377 )
   W( 60 ) = JVS( 1378 )
   W( 62 ) = JVS( 1379 )
   W( 63 ) = JVS( 1380 )
   W( 64 ) = JVS( 1381 )
   W( 65 ) = JVS( 1382 )
   W( 70 ) = JVS( 1383 )
   W( 71 ) = JVS( 1384 )
   W( 73 ) = JVS( 1385 )
   W( 74 ) = JVS( 1386 )
   W( 75 ) = JVS( 1387 )
   W( 76 ) = JVS( 1388 )
   W( 78 ) = JVS( 1389 )
   W( 81 ) = JVS( 1390 )
   W( 85 ) = JVS( 1391 )
   W( 89 ) = JVS( 1392 )
   W( 90 ) = JVS( 1393 )
   W( 91 ) = JVS( 1394 )
   W( 93 ) = JVS( 1395 )
   W( 94 ) = JVS( 1396 )
   W( 95 ) = JVS( 1397 )
   W( 96 ) = JVS( 1398 )
   W( 97 ) = JVS( 1399 )
   W( 98 ) = JVS( 1400 )
   W( 102 ) = JVS( 1401 )
   W( 103 ) = JVS( 1402 )
   W( 105 ) = JVS( 1403 )
   W( 106 ) = JVS( 1404 )
   W( 108 ) = JVS( 1405 )
   W( 113 ) = JVS( 1406 )
   W( 114 ) = JVS( 1407 )
   W( 116 ) = JVS( 1408 )
   W( 117 ) = JVS( 1409 )
   W( 118 ) = JVS( 1410 )
   W( 119 ) = JVS( 1411 )
   W( 120 ) = JVS( 1412 )
   W( 121 ) = JVS( 1413 )
   W( 122 ) = JVS( 1414 )
   W( 124 ) = JVS( 1415 )
   W( 126 ) = JVS( 1416 )
   W( 127 ) = JVS( 1417 )
   W( 128 ) = JVS( 1418 )
   W( 129 ) = JVS( 1419 )
   W( 130 ) = JVS( 1420 )
   W( 131 ) = JVS( 1421 )
   W( 132 ) = JVS( 1422 )
   W( 133 ) = JVS( 1423 )
   W( 134 ) = JVS( 1424 )
   W( 135 ) = JVS( 1425 )
   W( 136 ) = JVS( 1426 )
   W( 137 ) = JVS( 1427 )
   W( 138 ) = JVS( 1428 )
   W( 139 ) = JVS( 1429 )
   W( 140 ) = JVS( 1430 )
  a = -W( 21 ) / JVS(          155  )
  W( 21 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 156 )
  W( 138 ) = W( 138 ) + a*JVS( 157 )
  a = -W( 40 ) / JVS(          219  )
  W( 40 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 220 )
  W( 138 ) = W( 138 ) + a*JVS( 221 )
  a = -W( 51 ) / JVS(          261  )
  W( 51 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 262 )
  W( 139 ) = W( 139 ) + a*JVS( 263 )
  a = -W( 56 ) / JVS(          281  )
  W( 56 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 282 )
  W( 135 ) = W( 135 ) + a*JVS( 283 )
  W( 138 ) = W( 138 ) + a*JVS( 284 )
  a = -W( 60 ) / JVS(          305  )
  W( 60 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 306 )
  W( 136 ) = W( 136 ) + a*JVS( 307 )
  W( 139 ) = W( 139 ) + a*JVS( 308 )
  a = -W( 62 ) / JVS(          316  )
  W( 62 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 317 )
  W( 135 ) = W( 135 ) + a*JVS( 318 )
  W( 139 ) = W( 139 ) + a*JVS( 319 )
  a = -W( 63 ) / JVS(          323  )
  W( 63 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 324 )
  W( 134 ) = W( 134 ) + a*JVS( 325 )
  W( 136 ) = W( 136 ) + a*JVS( 326 )
  W( 138 ) = W( 138 ) + a*JVS( 327 )
  W( 139 ) = W( 139 ) + a*JVS( 328 )
  a = -W( 64 ) / JVS(          332  )
  W( 64 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 333 )
  W( 135 ) = W( 135 ) + a*JVS( 334 )
  W( 139 ) = W( 139 ) + a*JVS( 335 )
  a = -W( 65 ) / JVS(          336  )
  W( 65 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 337 )
  W( 135 ) = W( 135 ) + a*JVS( 338 )
  W( 139 ) = W( 139 ) + a*JVS( 339 )
  a = -W( 70 ) / JVS(          365  )
  W( 70 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 366 )
  W( 135 ) = W( 135 ) + a*JVS( 367 )
  W( 139 ) = W( 139 ) + a*JVS( 368 )
  a = -W( 71 ) / JVS(          369  )
  W( 71 ) = -a
  W( 105 ) = W( 105 ) + a*JVS( 370 )
  W( 113 ) = W( 113 ) + a*JVS( 371 )
  W( 117 ) = W( 117 ) + a*JVS( 372 )
  W( 133 ) = W( 133 ) + a*JVS( 373 )
  W( 134 ) = W( 134 ) + a*JVS( 374 )
  W( 139 ) = W( 139 ) + a*JVS( 375 )
  a = -W( 73 ) / JVS(          380  )
  W( 73 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 381 )
  W( 137 ) = W( 137 ) + a*JVS( 382 )
  W( 138 ) = W( 138 ) + a*JVS( 383 )
  a = -W( 74 ) / JVS(          386  )
  W( 74 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 387 )
  W( 135 ) = W( 135 ) + a*JVS( 388 )
  W( 139 ) = W( 139 ) + a*JVS( 389 )
  a = -W( 75 ) / JVS(          392  )
  W( 75 ) = -a
  W( 114 ) = W( 114 ) + a*JVS( 393 )
  W( 126 ) = W( 126 ) + a*JVS( 394 )
  W( 130 ) = W( 130 ) + a*JVS( 395 )
  W( 133 ) = W( 133 ) + a*JVS( 396 )
  W( 134 ) = W( 134 ) + a*JVS( 397 )
  W( 138 ) = W( 138 ) + a*JVS( 398 )
  a = -W( 76 ) / JVS(          399  )
  W( 76 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 400 )
  W( 134 ) = W( 134 ) + a*JVS( 401 )
  W( 135 ) = W( 135 ) + a*JVS( 402 )
  W( 139 ) = W( 139 ) + a*JVS( 403 )
  a = -W( 78 ) / JVS(          412  )
  W( 78 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 413 )
  W( 135 ) = W( 135 ) + a*JVS( 414 )
  W( 138 ) = W( 138 ) + a*JVS( 415 )
  W( 139 ) = W( 139 ) + a*JVS( 416 )
  a = -W( 81 ) / JVS(          425  )
  W( 81 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 426 )
  W( 134 ) = W( 134 ) + a*JVS( 427 )
  W( 138 ) = W( 138 ) + a*JVS( 428 )
  a = -W( 85 ) / JVS(          459  )
  W( 85 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 460 )
  W( 135 ) = W( 135 ) + a*JVS( 461 )
  W( 136 ) = W( 136 ) + a*JVS( 462 )
  W( 138 ) = W( 138 ) + a*JVS( 463 )
  W( 139 ) = W( 139 ) + a*JVS( 464 )
  a = -W( 89 ) / JVS(          478  )
  W( 89 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 479 )
  W( 134 ) = W( 134 ) + a*JVS( 480 )
  W( 135 ) = W( 135 ) + a*JVS( 481 )
  W( 139 ) = W( 139 ) + a*JVS( 482 )
  a = -W( 90 ) / JVS(          484  )
  W( 90 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 485 )
  W( 134 ) = W( 134 ) + a*JVS( 486 )
  W( 135 ) = W( 135 ) + a*JVS( 487 )
  W( 139 ) = W( 139 ) + a*JVS( 488 )
  a = -W( 91 ) / JVS(          491  )
  W( 91 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 492 )
  W( 94 ) = W( 94 ) + a*JVS( 493 )
  W( 134 ) = W( 134 ) + a*JVS( 494 )
  W( 135 ) = W( 135 ) + a*JVS( 495 )
  W( 138 ) = W( 138 ) + a*JVS( 496 )
  W( 139 ) = W( 139 ) + a*JVS( 497 )
  a = -W( 93 ) / JVS(          504  )
  W( 93 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 505 )
  W( 135 ) = W( 135 ) + a*JVS( 506 )
  W( 139 ) = W( 139 ) + a*JVS( 507 )
  a = -W( 94 ) / JVS(          510  )
  W( 94 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 511 )
  W( 135 ) = W( 135 ) + a*JVS( 512 )
  W( 139 ) = W( 139 ) + a*JVS( 513 )
  a = -W( 95 ) / JVS(          517  )
  W( 95 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 518 )
  W( 135 ) = W( 135 ) + a*JVS( 519 )
  W( 138 ) = W( 138 ) + a*JVS( 520 )
  W( 139 ) = W( 139 ) + a*JVS( 521 )
  a = -W( 96 ) / JVS(          523  )
  W( 96 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 524 )
  W( 130 ) = W( 130 ) + a*JVS( 525 )
  W( 134 ) = W( 134 ) + a*JVS( 526 )
  W( 135 ) = W( 135 ) + a*JVS( 527 )
  W( 136 ) = W( 136 ) + a*JVS( 528 )
  W( 138 ) = W( 138 ) + a*JVS( 529 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 533 )
  W( 135 ) = W( 135 ) + a*JVS( 534 )
  W( 139 ) = W( 139 ) + a*JVS( 535 )
  W( 140 ) = W( 140 ) + a*JVS( 536 )
  a = -W( 98 ) / JVS(          540  )
  W( 98 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 541 )
  W( 135 ) = W( 135 ) + a*JVS( 542 )
  W( 138 ) = W( 138 ) + a*JVS( 543 )
  W( 139 ) = W( 139 ) + a*JVS( 544 )
  a = -W( 102 ) / JVS(          594  )
  W( 102 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 595 )
  W( 134 ) = W( 134 ) + a*JVS( 596 )
  W( 135 ) = W( 135 ) + a*JVS( 597 )
  W( 139 ) = W( 139 ) + a*JVS( 598 )
  W( 140 ) = W( 140 ) + a*JVS( 599 )
  a = -W( 103 ) / JVS(          600  )
  W( 103 ) = -a
  W( 108 ) = W( 108 ) + a*JVS( 601 )
  W( 133 ) = W( 133 ) + a*JVS( 602 )
  W( 135 ) = W( 135 ) + a*JVS( 603 )
  W( 139 ) = W( 139 ) + a*JVS( 604 )
  a = -W( 105 ) / JVS(          610  )
  W( 105 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 611 )
  W( 134 ) = W( 134 ) + a*JVS( 612 )
  W( 135 ) = W( 135 ) + a*JVS( 613 )
  W( 136 ) = W( 136 ) + a*JVS( 614 )
  W( 139 ) = W( 139 ) + a*JVS( 615 )
  a = -W( 106 ) / JVS(          619  )
  W( 106 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 620 )
  W( 133 ) = W( 133 ) + a*JVS( 621 )
  W( 134 ) = W( 134 ) + a*JVS( 622 )
  W( 135 ) = W( 135 ) + a*JVS( 623 )
  W( 136 ) = W( 136 ) + a*JVS( 624 )
  W( 139 ) = W( 139 ) + a*JVS( 625 )
  W( 140 ) = W( 140 ) + a*JVS( 626 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  a = -W( 113 ) / JVS(          732  )
  W( 113 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 733 )
  W( 133 ) = W( 133 ) + a*JVS( 734 )
  W( 134 ) = W( 134 ) + a*JVS( 735 )
  W( 136 ) = W( 136 ) + a*JVS( 736 )
  a = -W( 114 ) / JVS(          748  )
  W( 114 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 749 )
  W( 128 ) = W( 128 ) + a*JVS( 750 )
  W( 133 ) = W( 133 ) + a*JVS( 751 )
  W( 134 ) = W( 134 ) + a*JVS( 752 )
  W( 135 ) = W( 135 ) + a*JVS( 753 )
  W( 136 ) = W( 136 ) + a*JVS( 754 )
  W( 139 ) = W( 139 ) + a*JVS( 755 )
  W( 140 ) = W( 140 ) + a*JVS( 756 )
  a = -W( 116 ) / JVS(          774  )
  W( 116 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 775 )
  W( 134 ) = W( 134 ) + a*JVS( 776 )
  W( 135 ) = W( 135 ) + a*JVS( 777 )
  W( 136 ) = W( 136 ) + a*JVS( 778 )
  W( 139 ) = W( 139 ) + a*JVS( 779 )
  W( 140 ) = W( 140 ) + a*JVS( 780 )
  a = -W( 117 ) / JVS(          788  )
  W( 117 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 789 )
  W( 134 ) = W( 134 ) + a*JVS( 790 )
  W( 135 ) = W( 135 ) + a*JVS( 791 )
  W( 136 ) = W( 136 ) + a*JVS( 792 )
  W( 139 ) = W( 139 ) + a*JVS( 793 )
  a = -W( 118 ) / JVS(          795  )
  W( 118 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 796 )
  W( 133 ) = W( 133 ) + a*JVS( 797 )
  W( 134 ) = W( 134 ) + a*JVS( 798 )
  W( 135 ) = W( 135 ) + a*JVS( 799 )
  W( 139 ) = W( 139 ) + a*JVS( 800 )
  W( 140 ) = W( 140 ) + a*JVS( 801 )
  a = -W( 119 ) / JVS(          808  )
  W( 119 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 809 )
  W( 134 ) = W( 134 ) + a*JVS( 810 )
  W( 135 ) = W( 135 ) + a*JVS( 811 )
  W( 136 ) = W( 136 ) + a*JVS( 812 )
  W( 139 ) = W( 139 ) + a*JVS( 813 )
  W( 140 ) = W( 140 ) + a*JVS( 814 )
  a = -W( 120 ) / JVS(          819  )
  W( 120 ) = -a
  W( 122 ) = W( 122 ) + a*JVS( 820 )
  W( 124 ) = W( 124 ) + a*JVS( 821 )
  W( 129 ) = W( 129 ) + a*JVS( 822 )
  W( 131 ) = W( 131 ) + a*JVS( 823 )
  W( 133 ) = W( 133 ) + a*JVS( 824 )
  W( 134 ) = W( 134 ) + a*JVS( 825 )
  W( 135 ) = W( 135 ) + a*JVS( 826 )
  W( 136 ) = W( 136 ) + a*JVS( 827 )
  W( 139 ) = W( 139 ) + a*JVS( 828 )
  W( 140 ) = W( 140 ) + a*JVS( 829 )
  a = -W( 121 ) / JVS(          842  )
  W( 121 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 843 )
  W( 133 ) = W( 133 ) + a*JVS( 844 )
  W( 134 ) = W( 134 ) + a*JVS( 845 )
  W( 135 ) = W( 135 ) + a*JVS( 846 )
  W( 136 ) = W( 136 ) + a*JVS( 847 )
  W( 139 ) = W( 139 ) + a*JVS( 848 )
  W( 140 ) = W( 140 ) + a*JVS( 849 )
  a = -W( 122 ) / JVS(          853  )
  W( 122 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 854 )
  W( 133 ) = W( 133 ) + a*JVS( 855 )
  W( 134 ) = W( 134 ) + a*JVS( 856 )
  W( 135 ) = W( 135 ) + a*JVS( 857 )
  W( 136 ) = W( 136 ) + a*JVS( 858 )
  W( 139 ) = W( 139 ) + a*JVS( 859 )
  W( 140 ) = W( 140 ) + a*JVS( 860 )
  a = -W( 124 ) / JVS(          881  )
  W( 124 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 882 )
  W( 134 ) = W( 134 ) + a*JVS( 883 )
  W( 135 ) = W( 135 ) + a*JVS( 884 )
  W( 136 ) = W( 136 ) + a*JVS( 885 )
  W( 139 ) = W( 139 ) + a*JVS( 886 )
  W( 140 ) = W( 140 ) + a*JVS( 887 )
  a = -W( 126 ) / JVS(          927  )
  W( 126 ) = -a
  W( 127 ) = W( 127 ) + a*JVS( 928 )
  W( 128 ) = W( 128 ) + a*JVS( 929 )
  W( 129 ) = W( 129 ) + a*JVS( 930 )
  W( 131 ) = W( 131 ) + a*JVS( 931 )
  W( 132 ) = W( 132 ) + a*JVS( 932 )
  W( 133 ) = W( 133 ) + a*JVS( 933 )
  W( 134 ) = W( 134 ) + a*JVS( 934 )
  W( 135 ) = W( 135 ) + a*JVS( 935 )
  W( 136 ) = W( 136 ) + a*JVS( 936 )
  W( 137 ) = W( 137 ) + a*JVS( 937 )
  W( 138 ) = W( 138 ) + a*JVS( 938 )
  W( 139 ) = W( 139 ) + a*JVS( 939 )
  W( 140 ) = W( 140 ) + a*JVS( 940 )
  a = -W( 127 ) / JVS(          947  )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 948 )
  W( 131 ) = W( 131 ) + a*JVS( 949 )
  W( 133 ) = W( 133 ) + a*JVS( 950 )
  W( 134 ) = W( 134 ) + a*JVS( 951 )
  W( 135 ) = W( 135 ) + a*JVS( 952 )
  W( 136 ) = W( 136 ) + a*JVS( 953 )
  W( 137 ) = W( 137 ) + a*JVS( 954 )
  W( 139 ) = W( 139 ) + a*JVS( 955 )
  W( 140 ) = W( 140 ) + a*JVS( 956 )
  a = -W( 128 ) / JVS(          960  )
  W( 128 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 961 )
  W( 133 ) = W( 133 ) + a*JVS( 962 )
  W( 134 ) = W( 134 ) + a*JVS( 963 )
  W( 135 ) = W( 135 ) + a*JVS( 964 )
  W( 136 ) = W( 136 ) + a*JVS( 965 )
  W( 137 ) = W( 137 ) + a*JVS( 966 )
  W( 139 ) = W( 139 ) + a*JVS( 967 )
  W( 140 ) = W( 140 ) + a*JVS( 968 )
  a = -W( 129 ) / JVS(          972  )
  W( 129 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 973 )
  W( 133 ) = W( 133 ) + a*JVS( 974 )
  W( 134 ) = W( 134 ) + a*JVS( 975 )
  W( 135 ) = W( 135 ) + a*JVS( 976 )
  W( 136 ) = W( 136 ) + a*JVS( 977 )
  W( 137 ) = W( 137 ) + a*JVS( 978 )
  W( 139 ) = W( 139 ) + a*JVS( 979 )
  W( 140 ) = W( 140 ) + a*JVS( 980 )
  a = -W( 130 ) / JVS(         1028  )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 1029 )
  W( 132 ) = W( 132 ) + a*JVS( 1030 )
  W( 133 ) = W( 133 ) + a*JVS( 1031 )
  W( 134 ) = W( 134 ) + a*JVS( 1032 )
  W( 135 ) = W( 135 ) + a*JVS( 1033 )
  W( 136 ) = W( 136 ) + a*JVS( 1034 )
  W( 137 ) = W( 137 ) + a*JVS( 1035 )
  W( 138 ) = W( 138 ) + a*JVS( 1036 )
  W( 139 ) = W( 139 ) + a*JVS( 1037 )
  W( 140 ) = W( 140 ) + a*JVS( 1038 )
  a = -W( 131 ) / JVS(         1041  )
  W( 131 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1042 )
  W( 134 ) = W( 134 ) + a*JVS( 1043 )
  W( 135 ) = W( 135 ) + a*JVS( 1044 )
  W( 136 ) = W( 136 ) + a*JVS( 1045 )
  W( 137 ) = W( 137 ) + a*JVS( 1046 )
  W( 139 ) = W( 139 ) + a*JVS( 1047 )
  W( 140 ) = W( 140 ) + a*JVS( 1048 )
  a = -W( 132 ) / JVS(         1055  )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1056 )
  W( 134 ) = W( 134 ) + a*JVS( 1057 )
  W( 135 ) = W( 135 ) + a*JVS( 1058 )
  W( 136 ) = W( 136 ) + a*JVS( 1059 )
  W( 137 ) = W( 137 ) + a*JVS( 1060 )
  W( 138 ) = W( 138 ) + a*JVS( 1061 )
  W( 139 ) = W( 139 ) + a*JVS( 1062 )
  W( 140 ) = W( 140 ) + a*JVS( 1063 )
  a = -W( 133 ) / JVS(         1090  )
  W( 133 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 1091 )
  W( 135 ) = W( 135 ) + a*JVS( 1092 )
  W( 136 ) = W( 136 ) + a*JVS( 1093 )
  W( 137 ) = W( 137 ) + a*JVS( 1094 )
  W( 138 ) = W( 138 ) + a*JVS( 1095 )
  W( 139 ) = W( 139 ) + a*JVS( 1096 )
  W( 140 ) = W( 140 ) + a*JVS( 1097 )
  a = -W( 134 ) / JVS(         1203  )
  W( 134 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 1204 )
  W( 136 ) = W( 136 ) + a*JVS( 1205 )
  W( 137 ) = W( 137 ) + a*JVS( 1206 )
  W( 138 ) = W( 138 ) + a*JVS( 1207 )
  W( 139 ) = W( 139 ) + a*JVS( 1208 )
  W( 140 ) = W( 140 ) + a*JVS( 1209 )
  a = -W( 135 ) / JVS(         1301  )
  W( 135 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 1302 )
  W( 137 ) = W( 137 ) + a*JVS( 1303 )
  W( 138 ) = W( 138 ) + a*JVS( 1304 )
  W( 139 ) = W( 139 ) + a*JVS( 1305 )
  W( 140 ) = W( 140 ) + a*JVS( 1306 )
  a = -W( 136 ) / JVS(         1327  )
  W( 136 ) = -a
  W( 137 ) = W( 137 ) + a*JVS( 1328 )
  W( 138 ) = W( 138 ) + a*JVS( 1329 )
  W( 139 ) = W( 139 ) + a*JVS( 1330 )
  W( 140 ) = W( 140 ) + a*JVS( 1331 )
  a = -W( 137 ) / JVS(         1370  )
  W( 137 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 1371 )
  W( 139 ) = W( 139 ) + a*JVS( 1372 )
  W( 140 ) = W( 140 ) + a*JVS( 1373 )
  JVS( 1374) = W( 21 )
  JVS( 1375) = W( 40 )
  JVS( 1376) = W( 51 )
  JVS( 1377) = W( 56 )
  JVS( 1378) = W( 60 )
  JVS( 1379) = W( 62 )
  JVS( 1380) = W( 63 )
  JVS( 1381) = W( 64 )
  JVS( 1382) = W( 65 )
  JVS( 1383) = W( 70 )
  JVS( 1384) = W( 71 )
  JVS( 1385) = W( 73 )
  JVS( 1386) = W( 74 )
  JVS( 1387) = W( 75 )
  JVS( 1388) = W( 76 )
  JVS( 1389) = W( 78 )
  JVS( 1390) = W( 81 )
  JVS( 1391) = W( 85 )
  JVS( 1392) = W( 89 )
  JVS( 1393) = W( 90 )
  JVS( 1394) = W( 91 )
  JVS( 1395) = W( 93 )
  JVS( 1396) = W( 94 )
  JVS( 1397) = W( 95 )
  JVS( 1398) = W( 96 )
  JVS( 1399) = W( 97 )
  JVS( 1400) = W( 98 )
  JVS( 1401) = W( 102 )
  JVS( 1402) = W( 103 )
  JVS( 1403) = W( 105 )
  JVS( 1404) = W( 106 )
  JVS( 1405) = W( 108 )
  JVS( 1406) = W( 113 )
  JVS( 1407) = W( 114 )
  JVS( 1408) = W( 116 )
  JVS( 1409) = W( 117 )
  JVS( 1410) = W( 118 )
  JVS( 1411) = W( 119 )
  JVS( 1412) = W( 120 )
  JVS( 1413) = W( 121 )
  JVS( 1414) = W( 122 )
  JVS( 1415) = W( 124 )
  JVS( 1416) = W( 126 )
  JVS( 1417) = W( 127 )
  JVS( 1418) = W( 128 )
  JVS( 1419) = W( 129 )
  JVS( 1420) = W( 130 )
  JVS( 1421) = W( 131 )
  JVS( 1422) = W( 132 )
  JVS( 1423) = W( 133 )
  JVS( 1424) = W( 134 )
  JVS( 1425) = W( 135 )
  JVS( 1426) = W( 136 )
  JVS( 1427) = W( 137 )
  JVS( 1428) = W( 138 )
  JVS( 1429) = W( 139 )
  JVS( 1430) = W( 140 )
  IF ( ABS(  JVS( 1478 )) < TINY(a) ) THEN
         IER = 139                                     
         RETURN
  END IF
   W( 1 ) = JVS( 1431 )
   W( 19 ) = JVS( 1432 )
   W( 51 ) = JVS( 1433 )
   W( 60 ) = JVS( 1434 )
   W( 62 ) = JVS( 1435 )
   W( 64 ) = JVS( 1436 )
   W( 65 ) = JVS( 1437 )
   W( 70 ) = JVS( 1438 )
   W( 74 ) = JVS( 1439 )
   W( 76 ) = JVS( 1440 )
   W( 78 ) = JVS( 1441 )
   W( 84 ) = JVS( 1442 )
   W( 85 ) = JVS( 1443 )
   W( 89 ) = JVS( 1444 )
   W( 90 ) = JVS( 1445 )
   W( 91 ) = JVS( 1446 )
   W( 93 ) = JVS( 1447 )
   W( 94 ) = JVS( 1448 )
   W( 95 ) = JVS( 1449 )
   W( 96 ) = JVS( 1450 )
   W( 97 ) = JVS( 1451 )
   W( 98 ) = JVS( 1452 )
   W( 102 ) = JVS( 1453 )
   W( 103 ) = JVS( 1454 )
   W( 105 ) = JVS( 1455 )
   W( 106 ) = JVS( 1456 )
   W( 108 ) = JVS( 1457 )
   W( 113 ) = JVS( 1458 )
   W( 116 ) = JVS( 1459 )
   W( 117 ) = JVS( 1460 )
   W( 118 ) = JVS( 1461 )
   W( 119 ) = JVS( 1462 )
   W( 121 ) = JVS( 1463 )
   W( 122 ) = JVS( 1464 )
   W( 124 ) = JVS( 1465 )
   W( 127 ) = JVS( 1466 )
   W( 128 ) = JVS( 1467 )
   W( 129 ) = JVS( 1468 )
   W( 130 ) = JVS( 1469 )
   W( 131 ) = JVS( 1470 )
   W( 132 ) = JVS( 1471 )
   W( 133 ) = JVS( 1472 )
   W( 134 ) = JVS( 1473 )
   W( 135 ) = JVS( 1474 )
   W( 136 ) = JVS( 1475 )
   W( 137 ) = JVS( 1476 )
   W( 138 ) = JVS( 1477 )
   W( 139 ) = JVS( 1478 )
   W( 140 ) = JVS( 1479 )
  a = -W( 1 ) / JVS(            1  )
  W( 1 ) = -a
  a = -W( 19 ) / JVS(          149  )
  W( 19 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 150 )
  a = -W( 51 ) / JVS(          261  )
  W( 51 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 262 )
  W( 139 ) = W( 139 ) + a*JVS( 263 )
  a = -W( 60 ) / JVS(          305  )
  W( 60 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 306 )
  W( 136 ) = W( 136 ) + a*JVS( 307 )
  W( 139 ) = W( 139 ) + a*JVS( 308 )
  a = -W( 62 ) / JVS(          316  )
  W( 62 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 317 )
  W( 135 ) = W( 135 ) + a*JVS( 318 )
  W( 139 ) = W( 139 ) + a*JVS( 319 )
  a = -W( 64 ) / JVS(          332  )
  W( 64 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 333 )
  W( 135 ) = W( 135 ) + a*JVS( 334 )
  W( 139 ) = W( 139 ) + a*JVS( 335 )
  a = -W( 65 ) / JVS(          336  )
  W( 65 ) = -a
  W( 130 ) = W( 130 ) + a*JVS( 337 )
  W( 135 ) = W( 135 ) + a*JVS( 338 )
  W( 139 ) = W( 139 ) + a*JVS( 339 )
  a = -W( 70 ) / JVS(          365  )
  W( 70 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 366 )
  W( 135 ) = W( 135 ) + a*JVS( 367 )
  W( 139 ) = W( 139 ) + a*JVS( 368 )
  a = -W( 74 ) / JVS(          386  )
  W( 74 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 387 )
  W( 135 ) = W( 135 ) + a*JVS( 388 )
  W( 139 ) = W( 139 ) + a*JVS( 389 )
  a = -W( 76 ) / JVS(          399  )
  W( 76 ) = -a
  W( 106 ) = W( 106 ) + a*JVS( 400 )
  W( 134 ) = W( 134 ) + a*JVS( 401 )
  W( 135 ) = W( 135 ) + a*JVS( 402 )
  W( 139 ) = W( 139 ) + a*JVS( 403 )
  a = -W( 78 ) / JVS(          412  )
  W( 78 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 413 )
  W( 135 ) = W( 135 ) + a*JVS( 414 )
  W( 138 ) = W( 138 ) + a*JVS( 415 )
  W( 139 ) = W( 139 ) + a*JVS( 416 )
  a = -W( 84 ) / JVS(          451  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 452 )
  W( 130 ) = W( 130 ) + a*JVS( 453 )
  W( 134 ) = W( 134 ) + a*JVS( 454 )
  W( 136 ) = W( 136 ) + a*JVS( 455 )
  a = -W( 85 ) / JVS(          459  )
  W( 85 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 460 )
  W( 135 ) = W( 135 ) + a*JVS( 461 )
  W( 136 ) = W( 136 ) + a*JVS( 462 )
  W( 138 ) = W( 138 ) + a*JVS( 463 )
  W( 139 ) = W( 139 ) + a*JVS( 464 )
  a = -W( 89 ) / JVS(          478  )
  W( 89 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 479 )
  W( 134 ) = W( 134 ) + a*JVS( 480 )
  W( 135 ) = W( 135 ) + a*JVS( 481 )
  W( 139 ) = W( 139 ) + a*JVS( 482 )
  a = -W( 90 ) / JVS(          484  )
  W( 90 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 485 )
  W( 134 ) = W( 134 ) + a*JVS( 486 )
  W( 135 ) = W( 135 ) + a*JVS( 487 )
  W( 139 ) = W( 139 ) + a*JVS( 488 )
  a = -W( 91 ) / JVS(          491  )
  W( 91 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 492 )
  W( 94 ) = W( 94 ) + a*JVS( 493 )
  W( 134 ) = W( 134 ) + a*JVS( 494 )
  W( 135 ) = W( 135 ) + a*JVS( 495 )
  W( 138 ) = W( 138 ) + a*JVS( 496 )
  W( 139 ) = W( 139 ) + a*JVS( 497 )
  a = -W( 93 ) / JVS(          504  )
  W( 93 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 505 )
  W( 135 ) = W( 135 ) + a*JVS( 506 )
  W( 139 ) = W( 139 ) + a*JVS( 507 )
  a = -W( 94 ) / JVS(          510  )
  W( 94 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 511 )
  W( 135 ) = W( 135 ) + a*JVS( 512 )
  W( 139 ) = W( 139 ) + a*JVS( 513 )
  a = -W( 95 ) / JVS(          517  )
  W( 95 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 518 )
  W( 135 ) = W( 135 ) + a*JVS( 519 )
  W( 138 ) = W( 138 ) + a*JVS( 520 )
  W( 139 ) = W( 139 ) + a*JVS( 521 )
  a = -W( 96 ) / JVS(          523  )
  W( 96 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 524 )
  W( 130 ) = W( 130 ) + a*JVS( 525 )
  W( 134 ) = W( 134 ) + a*JVS( 526 )
  W( 135 ) = W( 135 ) + a*JVS( 527 )
  W( 136 ) = W( 136 ) + a*JVS( 528 )
  W( 138 ) = W( 138 ) + a*JVS( 529 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 533 )
  W( 135 ) = W( 135 ) + a*JVS( 534 )
  W( 139 ) = W( 139 ) + a*JVS( 535 )
  W( 140 ) = W( 140 ) + a*JVS( 536 )
  a = -W( 98 ) / JVS(          540  )
  W( 98 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 541 )
  W( 135 ) = W( 135 ) + a*JVS( 542 )
  W( 138 ) = W( 138 ) + a*JVS( 543 )
  W( 139 ) = W( 139 ) + a*JVS( 544 )
  a = -W( 102 ) / JVS(          594  )
  W( 102 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 595 )
  W( 134 ) = W( 134 ) + a*JVS( 596 )
  W( 135 ) = W( 135 ) + a*JVS( 597 )
  W( 139 ) = W( 139 ) + a*JVS( 598 )
  W( 140 ) = W( 140 ) + a*JVS( 599 )
  a = -W( 103 ) / JVS(          600  )
  W( 103 ) = -a
  W( 108 ) = W( 108 ) + a*JVS( 601 )
  W( 133 ) = W( 133 ) + a*JVS( 602 )
  W( 135 ) = W( 135 ) + a*JVS( 603 )
  W( 139 ) = W( 139 ) + a*JVS( 604 )
  a = -W( 105 ) / JVS(          610  )
  W( 105 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 611 )
  W( 134 ) = W( 134 ) + a*JVS( 612 )
  W( 135 ) = W( 135 ) + a*JVS( 613 )
  W( 136 ) = W( 136 ) + a*JVS( 614 )
  W( 139 ) = W( 139 ) + a*JVS( 615 )
  a = -W( 106 ) / JVS(          619  )
  W( 106 ) = -a
  W( 116 ) = W( 116 ) + a*JVS( 620 )
  W( 133 ) = W( 133 ) + a*JVS( 621 )
  W( 134 ) = W( 134 ) + a*JVS( 622 )
  W( 135 ) = W( 135 ) + a*JVS( 623 )
  W( 136 ) = W( 136 ) + a*JVS( 624 )
  W( 139 ) = W( 139 ) + a*JVS( 625 )
  W( 140 ) = W( 140 ) + a*JVS( 626 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  a = -W( 113 ) / JVS(          732  )
  W( 113 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 733 )
  W( 133 ) = W( 133 ) + a*JVS( 734 )
  W( 134 ) = W( 134 ) + a*JVS( 735 )
  W( 136 ) = W( 136 ) + a*JVS( 736 )
  a = -W( 116 ) / JVS(          774  )
  W( 116 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 775 )
  W( 134 ) = W( 134 ) + a*JVS( 776 )
  W( 135 ) = W( 135 ) + a*JVS( 777 )
  W( 136 ) = W( 136 ) + a*JVS( 778 )
  W( 139 ) = W( 139 ) + a*JVS( 779 )
  W( 140 ) = W( 140 ) + a*JVS( 780 )
  a = -W( 117 ) / JVS(          788  )
  W( 117 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 789 )
  W( 134 ) = W( 134 ) + a*JVS( 790 )
  W( 135 ) = W( 135 ) + a*JVS( 791 )
  W( 136 ) = W( 136 ) + a*JVS( 792 )
  W( 139 ) = W( 139 ) + a*JVS( 793 )
  a = -W( 118 ) / JVS(          795  )
  W( 118 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 796 )
  W( 133 ) = W( 133 ) + a*JVS( 797 )
  W( 134 ) = W( 134 ) + a*JVS( 798 )
  W( 135 ) = W( 135 ) + a*JVS( 799 )
  W( 139 ) = W( 139 ) + a*JVS( 800 )
  W( 140 ) = W( 140 ) + a*JVS( 801 )
  a = -W( 119 ) / JVS(          808  )
  W( 119 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 809 )
  W( 134 ) = W( 134 ) + a*JVS( 810 )
  W( 135 ) = W( 135 ) + a*JVS( 811 )
  W( 136 ) = W( 136 ) + a*JVS( 812 )
  W( 139 ) = W( 139 ) + a*JVS( 813 )
  W( 140 ) = W( 140 ) + a*JVS( 814 )
  a = -W( 121 ) / JVS(          842  )
  W( 121 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 843 )
  W( 133 ) = W( 133 ) + a*JVS( 844 )
  W( 134 ) = W( 134 ) + a*JVS( 845 )
  W( 135 ) = W( 135 ) + a*JVS( 846 )
  W( 136 ) = W( 136 ) + a*JVS( 847 )
  W( 139 ) = W( 139 ) + a*JVS( 848 )
  W( 140 ) = W( 140 ) + a*JVS( 849 )
  a = -W( 122 ) / JVS(          853  )
  W( 122 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 854 )
  W( 133 ) = W( 133 ) + a*JVS( 855 )
  W( 134 ) = W( 134 ) + a*JVS( 856 )
  W( 135 ) = W( 135 ) + a*JVS( 857 )
  W( 136 ) = W( 136 ) + a*JVS( 858 )
  W( 139 ) = W( 139 ) + a*JVS( 859 )
  W( 140 ) = W( 140 ) + a*JVS( 860 )
  a = -W( 124 ) / JVS(          881  )
  W( 124 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 882 )
  W( 134 ) = W( 134 ) + a*JVS( 883 )
  W( 135 ) = W( 135 ) + a*JVS( 884 )
  W( 136 ) = W( 136 ) + a*JVS( 885 )
  W( 139 ) = W( 139 ) + a*JVS( 886 )
  W( 140 ) = W( 140 ) + a*JVS( 887 )
  a = -W( 127 ) / JVS(          947  )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 948 )
  W( 131 ) = W( 131 ) + a*JVS( 949 )
  W( 133 ) = W( 133 ) + a*JVS( 950 )
  W( 134 ) = W( 134 ) + a*JVS( 951 )
  W( 135 ) = W( 135 ) + a*JVS( 952 )
  W( 136 ) = W( 136 ) + a*JVS( 953 )
  W( 137 ) = W( 137 ) + a*JVS( 954 )
  W( 139 ) = W( 139 ) + a*JVS( 955 )
  W( 140 ) = W( 140 ) + a*JVS( 956 )
  a = -W( 128 ) / JVS(          960  )
  W( 128 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 961 )
  W( 133 ) = W( 133 ) + a*JVS( 962 )
  W( 134 ) = W( 134 ) + a*JVS( 963 )
  W( 135 ) = W( 135 ) + a*JVS( 964 )
  W( 136 ) = W( 136 ) + a*JVS( 965 )
  W( 137 ) = W( 137 ) + a*JVS( 966 )
  W( 139 ) = W( 139 ) + a*JVS( 967 )
  W( 140 ) = W( 140 ) + a*JVS( 968 )
  a = -W( 129 ) / JVS(          972  )
  W( 129 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 973 )
  W( 133 ) = W( 133 ) + a*JVS( 974 )
  W( 134 ) = W( 134 ) + a*JVS( 975 )
  W( 135 ) = W( 135 ) + a*JVS( 976 )
  W( 136 ) = W( 136 ) + a*JVS( 977 )
  W( 137 ) = W( 137 ) + a*JVS( 978 )
  W( 139 ) = W( 139 ) + a*JVS( 979 )
  W( 140 ) = W( 140 ) + a*JVS( 980 )
  a = -W( 130 ) / JVS(         1028  )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 1029 )
  W( 132 ) = W( 132 ) + a*JVS( 1030 )
  W( 133 ) = W( 133 ) + a*JVS( 1031 )
  W( 134 ) = W( 134 ) + a*JVS( 1032 )
  W( 135 ) = W( 135 ) + a*JVS( 1033 )
  W( 136 ) = W( 136 ) + a*JVS( 1034 )
  W( 137 ) = W( 137 ) + a*JVS( 1035 )
  W( 138 ) = W( 138 ) + a*JVS( 1036 )
  W( 139 ) = W( 139 ) + a*JVS( 1037 )
  W( 140 ) = W( 140 ) + a*JVS( 1038 )
  a = -W( 131 ) / JVS(         1041  )
  W( 131 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1042 )
  W( 134 ) = W( 134 ) + a*JVS( 1043 )
  W( 135 ) = W( 135 ) + a*JVS( 1044 )
  W( 136 ) = W( 136 ) + a*JVS( 1045 )
  W( 137 ) = W( 137 ) + a*JVS( 1046 )
  W( 139 ) = W( 139 ) + a*JVS( 1047 )
  W( 140 ) = W( 140 ) + a*JVS( 1048 )
  a = -W( 132 ) / JVS(         1055  )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1056 )
  W( 134 ) = W( 134 ) + a*JVS( 1057 )
  W( 135 ) = W( 135 ) + a*JVS( 1058 )
  W( 136 ) = W( 136 ) + a*JVS( 1059 )
  W( 137 ) = W( 137 ) + a*JVS( 1060 )
  W( 138 ) = W( 138 ) + a*JVS( 1061 )
  W( 139 ) = W( 139 ) + a*JVS( 1062 )
  W( 140 ) = W( 140 ) + a*JVS( 1063 )
  a = -W( 133 ) / JVS(         1090  )
  W( 133 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 1091 )
  W( 135 ) = W( 135 ) + a*JVS( 1092 )
  W( 136 ) = W( 136 ) + a*JVS( 1093 )
  W( 137 ) = W( 137 ) + a*JVS( 1094 )
  W( 138 ) = W( 138 ) + a*JVS( 1095 )
  W( 139 ) = W( 139 ) + a*JVS( 1096 )
  W( 140 ) = W( 140 ) + a*JVS( 1097 )
  a = -W( 134 ) / JVS(         1203  )
  W( 134 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 1204 )
  W( 136 ) = W( 136 ) + a*JVS( 1205 )
  W( 137 ) = W( 137 ) + a*JVS( 1206 )
  W( 138 ) = W( 138 ) + a*JVS( 1207 )
  W( 139 ) = W( 139 ) + a*JVS( 1208 )
  W( 140 ) = W( 140 ) + a*JVS( 1209 )
  a = -W( 135 ) / JVS(         1301  )
  W( 135 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 1302 )
  W( 137 ) = W( 137 ) + a*JVS( 1303 )
  W( 138 ) = W( 138 ) + a*JVS( 1304 )
  W( 139 ) = W( 139 ) + a*JVS( 1305 )
  W( 140 ) = W( 140 ) + a*JVS( 1306 )
  a = -W( 136 ) / JVS(         1327  )
  W( 136 ) = -a
  W( 137 ) = W( 137 ) + a*JVS( 1328 )
  W( 138 ) = W( 138 ) + a*JVS( 1329 )
  W( 139 ) = W( 139 ) + a*JVS( 1330 )
  W( 140 ) = W( 140 ) + a*JVS( 1331 )
  a = -W( 137 ) / JVS(         1370  )
  W( 137 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 1371 )
  W( 139 ) = W( 139 ) + a*JVS( 1372 )
  W( 140 ) = W( 140 ) + a*JVS( 1373 )
  a = -W( 138 ) / JVS(         1428  )
  W( 138 ) = -a
  W( 139 ) = W( 139 ) + a*JVS( 1429 )
  W( 140 ) = W( 140 ) + a*JVS( 1430 )
  JVS( 1431) = W( 1 )
  JVS( 1432) = W( 19 )
  JVS( 1433) = W( 51 )
  JVS( 1434) = W( 60 )
  JVS( 1435) = W( 62 )
  JVS( 1436) = W( 64 )
  JVS( 1437) = W( 65 )
  JVS( 1438) = W( 70 )
  JVS( 1439) = W( 74 )
  JVS( 1440) = W( 76 )
  JVS( 1441) = W( 78 )
  JVS( 1442) = W( 84 )
  JVS( 1443) = W( 85 )
  JVS( 1444) = W( 89 )
  JVS( 1445) = W( 90 )
  JVS( 1446) = W( 91 )
  JVS( 1447) = W( 93 )
  JVS( 1448) = W( 94 )
  JVS( 1449) = W( 95 )
  JVS( 1450) = W( 96 )
  JVS( 1451) = W( 97 )
  JVS( 1452) = W( 98 )
  JVS( 1453) = W( 102 )
  JVS( 1454) = W( 103 )
  JVS( 1455) = W( 105 )
  JVS( 1456) = W( 106 )
  JVS( 1457) = W( 108 )
  JVS( 1458) = W( 113 )
  JVS( 1459) = W( 116 )
  JVS( 1460) = W( 117 )
  JVS( 1461) = W( 118 )
  JVS( 1462) = W( 119 )
  JVS( 1463) = W( 121 )
  JVS( 1464) = W( 122 )
  JVS( 1465) = W( 124 )
  JVS( 1466) = W( 127 )
  JVS( 1467) = W( 128 )
  JVS( 1468) = W( 129 )
  JVS( 1469) = W( 130 )
  JVS( 1470) = W( 131 )
  JVS( 1471) = W( 132 )
  JVS( 1472) = W( 133 )
  JVS( 1473) = W( 134 )
  JVS( 1474) = W( 135 )
  JVS( 1475) = W( 136 )
  JVS( 1476) = W( 137 )
  JVS( 1477) = W( 138 )
  JVS( 1478) = W( 139 )
  JVS( 1479) = W( 140 )
  IF ( ABS(  JVS( 1515 )) < TINY(a) ) THEN
         IER = 140                                     
         RETURN
  END IF
   W( 53 ) = JVS( 1480 )
   W( 67 ) = JVS( 1481 )
   W( 68 ) = JVS( 1482 )
   W( 73 ) = JVS( 1483 )
   W( 77 ) = JVS( 1484 )
   W( 84 ) = JVS( 1485 )
   W( 95 ) = JVS( 1486 )
   W( 97 ) = JVS( 1487 )
   W( 98 ) = JVS( 1488 )
   W( 102 ) = JVS( 1489 )
   W( 108 ) = JVS( 1490 )
   W( 110 ) = JVS( 1491 )
   W( 113 ) = JVS( 1492 )
   W( 114 ) = JVS( 1493 )
   W( 115 ) = JVS( 1494 )
   W( 116 ) = JVS( 1495 )
   W( 117 ) = JVS( 1496 )
   W( 118 ) = JVS( 1497 )
   W( 119 ) = JVS( 1498 )
   W( 121 ) = JVS( 1499 )
   W( 122 ) = JVS( 1500 )
   W( 124 ) = JVS( 1501 )
   W( 127 ) = JVS( 1502 )
   W( 128 ) = JVS( 1503 )
   W( 129 ) = JVS( 1504 )
   W( 130 ) = JVS( 1505 )
   W( 131 ) = JVS( 1506 )
   W( 132 ) = JVS( 1507 )
   W( 133 ) = JVS( 1508 )
   W( 134 ) = JVS( 1509 )
   W( 135 ) = JVS( 1510 )
   W( 136 ) = JVS( 1511 )
   W( 137 ) = JVS( 1512 )
   W( 138 ) = JVS( 1513 )
   W( 139 ) = JVS( 1514 )
   W( 140 ) = JVS( 1515 )
  a = -W( 53 ) / JVS(          268  )
  W( 53 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 269 )
  W( 135 ) = W( 135 ) + a*JVS( 270 )
  W( 140 ) = W( 140 ) + a*JVS( 271 )
  a = -W( 67 ) / JVS(          345  )
  W( 67 ) = -a
  W( 108 ) = W( 108 ) + a*JVS( 346 )
  W( 113 ) = W( 113 ) + a*JVS( 347 )
  W( 132 ) = W( 132 ) + a*JVS( 348 )
  W( 134 ) = W( 134 ) + a*JVS( 349 )
  W( 135 ) = W( 135 ) + a*JVS( 350 )
  W( 136 ) = W( 136 ) + a*JVS( 351 )
  W( 137 ) = W( 137 ) + a*JVS( 352 )
  W( 140 ) = W( 140 ) + a*JVS( 353 )
  a = -W( 68 ) / JVS(          354  )
  W( 68 ) = -a
  W( 132 ) = W( 132 ) + a*JVS( 355 )
  W( 134 ) = W( 134 ) + a*JVS( 356 )
  W( 135 ) = W( 135 ) + a*JVS( 357 )
  W( 137 ) = W( 137 ) + a*JVS( 358 )
  a = -W( 73 ) / JVS(          380  )
  W( 73 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 381 )
  W( 137 ) = W( 137 ) + a*JVS( 382 )
  W( 138 ) = W( 138 ) + a*JVS( 383 )
  a = -W( 77 ) / JVS(          404  )
  W( 77 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 405 )
  W( 113 ) = W( 113 ) + a*JVS( 406 )
  W( 134 ) = W( 134 ) + a*JVS( 407 )
  W( 136 ) = W( 136 ) + a*JVS( 408 )
  a = -W( 84 ) / JVS(          451  )
  W( 84 ) = -a
  W( 113 ) = W( 113 ) + a*JVS( 452 )
  W( 130 ) = W( 130 ) + a*JVS( 453 )
  W( 134 ) = W( 134 ) + a*JVS( 454 )
  W( 136 ) = W( 136 ) + a*JVS( 455 )
  a = -W( 95 ) / JVS(          517  )
  W( 95 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 518 )
  W( 135 ) = W( 135 ) + a*JVS( 519 )
  W( 138 ) = W( 138 ) + a*JVS( 520 )
  W( 139 ) = W( 139 ) + a*JVS( 521 )
  a = -W( 97 ) / JVS(          532  )
  W( 97 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 533 )
  W( 135 ) = W( 135 ) + a*JVS( 534 )
  W( 139 ) = W( 139 ) + a*JVS( 535 )
  W( 140 ) = W( 140 ) + a*JVS( 536 )
  a = -W( 98 ) / JVS(          540  )
  W( 98 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 541 )
  W( 135 ) = W( 135 ) + a*JVS( 542 )
  W( 138 ) = W( 138 ) + a*JVS( 543 )
  W( 139 ) = W( 139 ) + a*JVS( 544 )
  a = -W( 102 ) / JVS(          594  )
  W( 102 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 595 )
  W( 134 ) = W( 134 ) + a*JVS( 596 )
  W( 135 ) = W( 135 ) + a*JVS( 597 )
  W( 139 ) = W( 139 ) + a*JVS( 598 )
  W( 140 ) = W( 140 ) + a*JVS( 599 )
  a = -W( 108 ) / JVS(          631  )
  W( 108 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 632 )
  W( 134 ) = W( 134 ) + a*JVS( 633 )
  W( 136 ) = W( 136 ) + a*JVS( 634 )
  a = -W( 110 ) / JVS(          655  )
  W( 110 ) = -a
  W( 115 ) = W( 115 ) + a*JVS( 656 )
  W( 116 ) = W( 116 ) + a*JVS( 657 )
  W( 117 ) = W( 117 ) + a*JVS( 658 )
  W( 118 ) = W( 118 ) + a*JVS( 659 )
  W( 119 ) = W( 119 ) + a*JVS( 660 )
  W( 133 ) = W( 133 ) + a*JVS( 661 )
  W( 134 ) = W( 134 ) + a*JVS( 662 )
  W( 135 ) = W( 135 ) + a*JVS( 663 )
  W( 136 ) = W( 136 ) + a*JVS( 664 )
  W( 139 ) = W( 139 ) + a*JVS( 665 )
  W( 140 ) = W( 140 ) + a*JVS( 666 )
  a = -W( 113 ) / JVS(          732  )
  W( 113 ) = -a
  W( 128 ) = W( 128 ) + a*JVS( 733 )
  W( 133 ) = W( 133 ) + a*JVS( 734 )
  W( 134 ) = W( 134 ) + a*JVS( 735 )
  W( 136 ) = W( 136 ) + a*JVS( 736 )
  a = -W( 114 ) / JVS(          748  )
  W( 114 ) = -a
  W( 117 ) = W( 117 ) + a*JVS( 749 )
  W( 128 ) = W( 128 ) + a*JVS( 750 )
  W( 133 ) = W( 133 ) + a*JVS( 751 )
  W( 134 ) = W( 134 ) + a*JVS( 752 )
  W( 135 ) = W( 135 ) + a*JVS( 753 )
  W( 136 ) = W( 136 ) + a*JVS( 754 )
  W( 139 ) = W( 139 ) + a*JVS( 755 )
  W( 140 ) = W( 140 ) + a*JVS( 756 )
  a = -W( 115 ) / JVS(          763  )
  W( 115 ) = -a
  W( 118 ) = W( 118 ) + a*JVS( 764 )
  W( 124 ) = W( 124 ) + a*JVS( 765 )
  W( 133 ) = W( 133 ) + a*JVS( 766 )
  W( 134 ) = W( 134 ) + a*JVS( 767 )
  W( 135 ) = W( 135 ) + a*JVS( 768 )
  W( 136 ) = W( 136 ) + a*JVS( 769 )
  W( 139 ) = W( 139 ) + a*JVS( 770 )
  W( 140 ) = W( 140 ) + a*JVS( 771 )
  a = -W( 116 ) / JVS(          774  )
  W( 116 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 775 )
  W( 134 ) = W( 134 ) + a*JVS( 776 )
  W( 135 ) = W( 135 ) + a*JVS( 777 )
  W( 136 ) = W( 136 ) + a*JVS( 778 )
  W( 139 ) = W( 139 ) + a*JVS( 779 )
  W( 140 ) = W( 140 ) + a*JVS( 780 )
  a = -W( 117 ) / JVS(          788  )
  W( 117 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 789 )
  W( 134 ) = W( 134 ) + a*JVS( 790 )
  W( 135 ) = W( 135 ) + a*JVS( 791 )
  W( 136 ) = W( 136 ) + a*JVS( 792 )
  W( 139 ) = W( 139 ) + a*JVS( 793 )
  a = -W( 118 ) / JVS(          795  )
  W( 118 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 796 )
  W( 133 ) = W( 133 ) + a*JVS( 797 )
  W( 134 ) = W( 134 ) + a*JVS( 798 )
  W( 135 ) = W( 135 ) + a*JVS( 799 )
  W( 139 ) = W( 139 ) + a*JVS( 800 )
  W( 140 ) = W( 140 ) + a*JVS( 801 )
  a = -W( 119 ) / JVS(          808  )
  W( 119 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 809 )
  W( 134 ) = W( 134 ) + a*JVS( 810 )
  W( 135 ) = W( 135 ) + a*JVS( 811 )
  W( 136 ) = W( 136 ) + a*JVS( 812 )
  W( 139 ) = W( 139 ) + a*JVS( 813 )
  W( 140 ) = W( 140 ) + a*JVS( 814 )
  a = -W( 121 ) / JVS(          842  )
  W( 121 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 843 )
  W( 133 ) = W( 133 ) + a*JVS( 844 )
  W( 134 ) = W( 134 ) + a*JVS( 845 )
  W( 135 ) = W( 135 ) + a*JVS( 846 )
  W( 136 ) = W( 136 ) + a*JVS( 847 )
  W( 139 ) = W( 139 ) + a*JVS( 848 )
  W( 140 ) = W( 140 ) + a*JVS( 849 )
  a = -W( 122 ) / JVS(          853  )
  W( 122 ) = -a
  W( 124 ) = W( 124 ) + a*JVS( 854 )
  W( 133 ) = W( 133 ) + a*JVS( 855 )
  W( 134 ) = W( 134 ) + a*JVS( 856 )
  W( 135 ) = W( 135 ) + a*JVS( 857 )
  W( 136 ) = W( 136 ) + a*JVS( 858 )
  W( 139 ) = W( 139 ) + a*JVS( 859 )
  W( 140 ) = W( 140 ) + a*JVS( 860 )
  a = -W( 124 ) / JVS(          881  )
  W( 124 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 882 )
  W( 134 ) = W( 134 ) + a*JVS( 883 )
  W( 135 ) = W( 135 ) + a*JVS( 884 )
  W( 136 ) = W( 136 ) + a*JVS( 885 )
  W( 139 ) = W( 139 ) + a*JVS( 886 )
  W( 140 ) = W( 140 ) + a*JVS( 887 )
  a = -W( 127 ) / JVS(          947  )
  W( 127 ) = -a
  W( 129 ) = W( 129 ) + a*JVS( 948 )
  W( 131 ) = W( 131 ) + a*JVS( 949 )
  W( 133 ) = W( 133 ) + a*JVS( 950 )
  W( 134 ) = W( 134 ) + a*JVS( 951 )
  W( 135 ) = W( 135 ) + a*JVS( 952 )
  W( 136 ) = W( 136 ) + a*JVS( 953 )
  W( 137 ) = W( 137 ) + a*JVS( 954 )
  W( 139 ) = W( 139 ) + a*JVS( 955 )
  W( 140 ) = W( 140 ) + a*JVS( 956 )
  a = -W( 128 ) / JVS(          960  )
  W( 128 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 961 )
  W( 133 ) = W( 133 ) + a*JVS( 962 )
  W( 134 ) = W( 134 ) + a*JVS( 963 )
  W( 135 ) = W( 135 ) + a*JVS( 964 )
  W( 136 ) = W( 136 ) + a*JVS( 965 )
  W( 137 ) = W( 137 ) + a*JVS( 966 )
  W( 139 ) = W( 139 ) + a*JVS( 967 )
  W( 140 ) = W( 140 ) + a*JVS( 968 )
  a = -W( 129 ) / JVS(          972  )
  W( 129 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 973 )
  W( 133 ) = W( 133 ) + a*JVS( 974 )
  W( 134 ) = W( 134 ) + a*JVS( 975 )
  W( 135 ) = W( 135 ) + a*JVS( 976 )
  W( 136 ) = W( 136 ) + a*JVS( 977 )
  W( 137 ) = W( 137 ) + a*JVS( 978 )
  W( 139 ) = W( 139 ) + a*JVS( 979 )
  W( 140 ) = W( 140 ) + a*JVS( 980 )
  a = -W( 130 ) / JVS(         1028  )
  W( 130 ) = -a
  W( 131 ) = W( 131 ) + a*JVS( 1029 )
  W( 132 ) = W( 132 ) + a*JVS( 1030 )
  W( 133 ) = W( 133 ) + a*JVS( 1031 )
  W( 134 ) = W( 134 ) + a*JVS( 1032 )
  W( 135 ) = W( 135 ) + a*JVS( 1033 )
  W( 136 ) = W( 136 ) + a*JVS( 1034 )
  W( 137 ) = W( 137 ) + a*JVS( 1035 )
  W( 138 ) = W( 138 ) + a*JVS( 1036 )
  W( 139 ) = W( 139 ) + a*JVS( 1037 )
  W( 140 ) = W( 140 ) + a*JVS( 1038 )
  a = -W( 131 ) / JVS(         1041  )
  W( 131 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1042 )
  W( 134 ) = W( 134 ) + a*JVS( 1043 )
  W( 135 ) = W( 135 ) + a*JVS( 1044 )
  W( 136 ) = W( 136 ) + a*JVS( 1045 )
  W( 137 ) = W( 137 ) + a*JVS( 1046 )
  W( 139 ) = W( 139 ) + a*JVS( 1047 )
  W( 140 ) = W( 140 ) + a*JVS( 1048 )
  a = -W( 132 ) / JVS(         1055  )
  W( 132 ) = -a
  W( 133 ) = W( 133 ) + a*JVS( 1056 )
  W( 134 ) = W( 134 ) + a*JVS( 1057 )
  W( 135 ) = W( 135 ) + a*JVS( 1058 )
  W( 136 ) = W( 136 ) + a*JVS( 1059 )
  W( 137 ) = W( 137 ) + a*JVS( 1060 )
  W( 138 ) = W( 138 ) + a*JVS( 1061 )
  W( 139 ) = W( 139 ) + a*JVS( 1062 )
  W( 140 ) = W( 140 ) + a*JVS( 1063 )
  a = -W( 133 ) / JVS(         1090  )
  W( 133 ) = -a
  W( 134 ) = W( 134 ) + a*JVS( 1091 )
  W( 135 ) = W( 135 ) + a*JVS( 1092 )
  W( 136 ) = W( 136 ) + a*JVS( 1093 )
  W( 137 ) = W( 137 ) + a*JVS( 1094 )
  W( 138 ) = W( 138 ) + a*JVS( 1095 )
  W( 139 ) = W( 139 ) + a*JVS( 1096 )
  W( 140 ) = W( 140 ) + a*JVS( 1097 )
  a = -W( 134 ) / JVS(         1203  )
  W( 134 ) = -a
  W( 135 ) = W( 135 ) + a*JVS( 1204 )
  W( 136 ) = W( 136 ) + a*JVS( 1205 )
  W( 137 ) = W( 137 ) + a*JVS( 1206 )
  W( 138 ) = W( 138 ) + a*JVS( 1207 )
  W( 139 ) = W( 139 ) + a*JVS( 1208 )
  W( 140 ) = W( 140 ) + a*JVS( 1209 )
  a = -W( 135 ) / JVS(         1301  )
  W( 135 ) = -a
  W( 136 ) = W( 136 ) + a*JVS( 1302 )
  W( 137 ) = W( 137 ) + a*JVS( 1303 )
  W( 138 ) = W( 138 ) + a*JVS( 1304 )
  W( 139 ) = W( 139 ) + a*JVS( 1305 )
  W( 140 ) = W( 140 ) + a*JVS( 1306 )
  a = -W( 136 ) / JVS(         1327  )
  W( 136 ) = -a
  W( 137 ) = W( 137 ) + a*JVS( 1328 )
  W( 138 ) = W( 138 ) + a*JVS( 1329 )
  W( 139 ) = W( 139 ) + a*JVS( 1330 )
  W( 140 ) = W( 140 ) + a*JVS( 1331 )
  a = -W( 137 ) / JVS(         1370  )
  W( 137 ) = -a
  W( 138 ) = W( 138 ) + a*JVS( 1371 )
  W( 139 ) = W( 139 ) + a*JVS( 1372 )
  W( 140 ) = W( 140 ) + a*JVS( 1373 )
  a = -W( 138 ) / JVS(         1428  )
  W( 138 ) = -a
  W( 139 ) = W( 139 ) + a*JVS( 1429 )
  W( 140 ) = W( 140 ) + a*JVS( 1430 )
  a = -W( 139 ) / JVS(         1478  )
  W( 139 ) = -a
  W( 140 ) = W( 140 ) + a*JVS( 1479 )
  JVS( 1480) = W( 53 )
  JVS( 1481) = W( 67 )
  JVS( 1482) = W( 68 )
  JVS( 1483) = W( 73 )
  JVS( 1484) = W( 77 )
  JVS( 1485) = W( 84 )
  JVS( 1486) = W( 95 )
  JVS( 1487) = W( 97 )
  JVS( 1488) = W( 98 )
  JVS( 1489) = W( 102 )
  JVS( 1490) = W( 108 )
  JVS( 1491) = W( 110 )
  JVS( 1492) = W( 113 )
  JVS( 1493) = W( 114 )
  JVS( 1494) = W( 115 )
  JVS( 1495) = W( 116 )
  JVS( 1496) = W( 117 )
  JVS( 1497) = W( 118 )
  JVS( 1498) = W( 119 )
  JVS( 1499) = W( 121 )
  JVS( 1500) = W( 122 )
  JVS( 1501) = W( 124 )
  JVS( 1502) = W( 127 )
  JVS( 1503) = W( 128 )
  JVS( 1504) = W( 129 )
  JVS( 1505) = W( 130 )
  JVS( 1506) = W( 131 )
  JVS( 1507) = W( 132 )
  JVS( 1508) = W( 133 )
  JVS( 1509) = W( 134 )
  JVS( 1510) = W( 135 )
  JVS( 1511) = W( 136 )
  JVS( 1512) = W( 137 )
  JVS( 1513) = W( 138 )
  JVS( 1514) = W( 139 )
  JVS( 1515) = W( 140 )
   
   END SUBROUTINE decomp_mozart_mosaic_4bin_aq
 


END MODULE mozart_mosaic_4bin_aq_Integrator
