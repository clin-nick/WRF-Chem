
































MODULE cb05_sorg_aq_Integrator

 USE cb05_sorg_aq_Parameters
 USE cb05_sorg_aq_Precision
 USE cb05_sorg_aq_JacobianSP

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

SUBROUTINE  cb05_sorg_aq_INTEGRATE( TIN, TOUT, &
  FIX, VAR,  RCONST, ATOL, RTOL, IRR_WRK,  &
  ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U  )

   USE cb05_sorg_aq_Parameters

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

   CALL cb05_sorg_aq_Rosenbrock(VAR, FIX, RCONST, TIN,TOUT,   &
         ATOL,RTOL,               &
         RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)

   STEPMIN = RCNTRL(ihexit)
   
   IF (PRESENT(ISTATUS_U)) ISTATUS_U(:) = ISTATUS(:)
   IF (PRESENT(RSTATUS_U)) RSTATUS_U(:) = RSTATUS(:)
   IF (PRESENT(IERR_U))    IERR_U       = IERR

END SUBROUTINE  cb05_sorg_aq_INTEGRATE


SUBROUTINE  cb05_sorg_aq_Rosenbrock(Y, FIX, RCONST, Tstart,Tend, &
           AbsTol,RelTol,            &
           RCNTRL,ICNTRL,RSTATUS,ISTATUS,IRR_WRK,IERR)







































































































  USE cb05_sorg_aq_Parameters

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
      CALL cb05_sorg_aq_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
      RETURN
   END IF


   IF (ICNTRL(4) == 0) THEN
      Max_no_steps = 100000
   ELSEIF (ICNTRL(4) > 0) THEN
      Max_no_steps=ICNTRL(4)
   ELSE
      PRINT * ,'User-selected max no. of steps: ICNTRL(4)=',ICNTRL(4)
      CALL cb05_sorg_aq_ros_ErrorMsg(-1,Tstart,ZERO,IERR)
      RETURN
   END IF


   Roundoff = cb05_sorg_aq_WLAMCH('E')


   IF (RCNTRL(1) == ZERO) THEN
      Hmin = ZERO
   ELSEIF (RCNTRL(1) > ZERO) THEN
      Hmin = RCNTRL(1)
   ELSE
      PRINT * , 'User-selected Hmin: RCNTRL(1)=', RCNTRL(1)
      CALL cb05_sorg_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(2) == ZERO) THEN
      Hmax = ABS(Tend-Tstart)
   ELSEIF (RCNTRL(2) > ZERO) THEN
      Hmax = MIN(ABS(RCNTRL(2)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hmax: RCNTRL(2)=', RCNTRL(2)
      CALL cb05_sorg_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(3) == ZERO) THEN
      Hstart = MAX(Hmin,DeltaMin)
   ELSEIF (RCNTRL(3) > ZERO) THEN
      Hstart = MIN(ABS(RCNTRL(3)),ABS(Tend-Tstart))
   ELSE
      PRINT * , 'User-selected Hstart: RCNTRL(3)=', RCNTRL(3)
      CALL cb05_sorg_aq_ros_ErrorMsg(-3,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(4) == ZERO) THEN
      FacMin = 0.2_dp
   ELSEIF (RCNTRL(4) > ZERO) THEN
      FacMin = RCNTRL(4)
   ELSE
      PRINT * , 'User-selected FacMin: RCNTRL(4)=', RCNTRL(4)
      CALL cb05_sorg_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF
   IF (RCNTRL(5) == ZERO) THEN
      FacMax = 6.0_dp
   ELSEIF (RCNTRL(5) > ZERO) THEN
      FacMax = RCNTRL(5)
   ELSE
      PRINT * , 'User-selected FacMax: RCNTRL(5)=', RCNTRL(5)
      CALL cb05_sorg_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(6) == ZERO) THEN
      FacRej = 0.1_dp
   ELSEIF (RCNTRL(6) > ZERO) THEN
      FacRej = RCNTRL(6)
   ELSE
      PRINT * , 'User-selected FacRej: RCNTRL(6)=', RCNTRL(6)
      CALL cb05_sorg_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

   IF (RCNTRL(7) == ZERO) THEN
      FacSafe = 0.9_dp
   ELSEIF (RCNTRL(7) > ZERO) THEN
      FacSafe = RCNTRL(7)
   ELSE
      PRINT * , 'User-selected FacSafe: RCNTRL(7)=', RCNTRL(7)
      CALL cb05_sorg_aq_ros_ErrorMsg(-4,Tstart,ZERO,IERR)
      RETURN
   END IF

    DO i=1,UplimTol
      IF ( (AbsTol(i) <= ZERO) .OR. (RelTol(i) <= 10.0_dp*Roundoff) &
         .OR. (RelTol(i) >= 1.0_dp) ) THEN
        PRINT * , ' AbsTol(',i,') = ',AbsTol(i)
        PRINT * , ' RelTol(',i,') = ',RelTol(i)
        CALL cb05_sorg_aq_ros_ErrorMsg(-5,Tstart,ZERO,IERR)
        RETURN
      END IF
    END DO



   SELECT CASE (Method)
     CASE (1)
       CALL cb05_sorg_aq_Ros2(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (2)
       CALL cb05_sorg_aq_Ros3(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (3)
       CALL cb05_sorg_aq_Ros4(ros_S, ros_A, ros_C, ros_M, ros_E,   &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (4)
       CALL cb05_sorg_aq_Rodas3(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE (5)
       CALL cb05_sorg_aq_Rodas4(ros_S, ros_A, ros_C, ros_M, ros_E, &
          ros_Alpha, ros_Gamma, ros_NewF, ros_ELO, ros_Name)
     CASE DEFAULT
       PRINT * , 'Unknown Rosenbrock method: ICNTRL(4)=', Method
       CALL cb05_sorg_aq_ros_ErrorMsg(-2,Tstart,ZERO,IERR)
       RETURN
   END SELECT


   CALL cb05_sorg_aq_ros_Integrator(Y,Tstart,Tend,Texit,      &
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



 SUBROUTINE  cb05_sorg_aq_ros_ErrorMsg(Code,T,H,IERR)



   USE cb05_sorg_aq_Precision

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

 END SUBROUTINE  cb05_sorg_aq_ros_ErrorMsg


 SUBROUTINE  cb05_sorg_aq_ros_Integrator (Y, Tstart, Tend, T,     &
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
      CALL cb05_sorg_aq_ros_ErrorMsg(-6,T,H,IERR)
      RETURN
   END IF
   IF ( ((T+0.1_dp*H) == T).OR.(H <= Roundoff) ) THEN  
      CALL cb05_sorg_aq_ros_ErrorMsg(-7,T,H,IERR)
      RETURN
   END IF


   Hexit = H
   H = MIN(H,ABS(Tend-T))


   CALL cb05_sorg_aq_FunTemplate(T,Y,Fcn0, RCONST, FIX, Nfun)
   IF( T == Tstart ) THEN
     CALL cb05_sorg_aq_IRRFun( Y, FIX, RCONST, IRR_WRK )
   ENDIF


   IF (.NOT.Autonomous) THEN
      CALL cb05_sorg_aq_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )
   END IF


   CALL cb05_sorg_aq_JacTemplate(T,Y,Jac0, FIX, Njac, RCONST)


UntilAccepted: DO

   CALL cb05_sorg_aq_ros_PrepareMatrix(H,Direction,ros_Gamma(1), &
          Jac0,Ghimj,Pivot,Singular, Ndec,  Nsng )
   IF (Singular) THEN 
       CALL cb05_sorg_aq_ros_ErrorMsg(-8,T,H,IERR)
       RETURN
   END IF


Stage: DO istage = 1, ros_S

      
       ioffset = NVAR*(istage-1)

      
       IF ( istage == 1 ) THEN
         CALL cb05_sorg_aq_WCOPY(NVAR,Fcn0,1,Fcn,1)
      
       ELSEIF ( ros_NewF(istage) ) THEN
         CALL cb05_sorg_aq_WCOPY(NVAR,Y,1,Ynew,1)
         DO j = 1, istage-1
           CALL cb05_sorg_aq_WAXPY(NVAR,ros_A((istage-1)*(istage-2)/2+j), &
            K(NVAR*(j-1)+1),1,Ynew,1)
         END DO
         Tau = T + ros_Alpha(istage)*Direction*H
         CALL cb05_sorg_aq_FunTemplate(Tau,Ynew,Fcn, RCONST, FIX, Nfun)
       END IF 
       CALL cb05_sorg_aq_WCOPY(NVAR,Fcn,1,K(ioffset+1),1)
       DO j = 1, istage-1
         HC = ros_C((istage-1)*(istage-2)/2+j)/(Direction*H)
         CALL cb05_sorg_aq_WAXPY(NVAR,HC,K(NVAR*(j-1)+1),1,K(ioffset+1),1)
       END DO
       IF ((.NOT. Autonomous).AND.(ros_Gamma(istage).NE.ZERO)) THEN
         HG = Direction*H*ros_Gamma(istage)
         CALL cb05_sorg_aq_WAXPY(NVAR,HG,dFdT,1,K(ioffset+1),1)
       END IF
       CALL cb05_sorg_aq_ros_Solve(Ghimj, Pivot, K(ioffset+1), Nsol)

   END DO Stage



   CALL cb05_sorg_aq_WCOPY(NVAR,Y,1,Ynew,1)
   DO j=1,ros_S
         CALL cb05_sorg_aq_WAXPY(NVAR,ros_M(j),K(NVAR*(j-1)+1),1,Ynew,1)
   END DO


   CALL cb05_sorg_aq_WSCAL(NVAR,ZERO,Yerr,1)
   DO j=1,ros_S
        CALL cb05_sorg_aq_WAXPY(NVAR,ros_E(j),K(NVAR*(j-1)+1),1,Yerr,1)
   END DO
   Err = cb05_sorg_aq_ros_ErrorNorm ( Y, Ynew, Yerr, AbsTol, RelTol, VectorTol )


   Fac  = MIN(FacMax,MAX(FacMin,FacSafe/Err**(ONE/ros_ELO)))
   Hnew = H*Fac


   Nstp = Nstp+1
   IF ( (Err <= ONE).OR.(H <= Hmin) ) THEN  
      Nacc = Nacc+1
      CALL cb05_sorg_aq_WCOPY(NVAR,Ynew,1,Y,1)
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

  END SUBROUTINE  cb05_sorg_aq_ros_Integrator



  REAL(kind=dp) FUNCTION  cb05_sorg_aq_ros_ErrorNorm ( Y, Ynew, Yerr, &
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

    cb05_sorg_aq_ros_ErrorNorm = Err

  END FUNCTION  cb05_sorg_aq_ros_ErrorNorm



  SUBROUTINE cb05_sorg_aq_ros_FunTimeDeriv ( T, Roundoff, Y, &
                Fcn0, dFdT, RCONST, FIX, Nfun )



   IMPLICIT NONE


   REAL(kind=dp), INTENT(IN) :: T, Roundoff, Y(NVAR), Fcn0(NVAR)
   REAL(kind=dp), INTENT(IN) :: RCONST(NREACT), FIX(NFIX)

   REAL(kind=dp), INTENT(OUT) :: dFdT(NVAR)

   INTEGER, INTENT(INOUT) ::Nfun

   REAL(kind=dp) :: Delta
   REAL(kind=dp), PARAMETER :: ONE = 1.0_dp, DeltaMin = 1.0E-6_dp

   Delta = SQRT(Roundoff)*MAX(DeltaMin,ABS(T))
   CALL cb05_sorg_aq_FunTemplate(T+Delta,Y,dFdT, RCONST, FIX, Nfun)
   CALL cb05_sorg_aq_WAXPY(NVAR,(-ONE),Fcn0,1,dFdT,1)
   CALL cb05_sorg_aq_WSCAL(NVAR,(ONE/Delta),dFdT,1)

  END SUBROUTINE  cb05_sorg_aq_ros_FunTimeDeriv



  SUBROUTINE  cb05_sorg_aq_ros_PrepareMatrix ( H, Direction, gam, &
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


     CALL cb05_sorg_aq_WCOPY(LU_NONZERO,Jac0,1,Ghimj,1)
     CALL cb05_sorg_aq_WSCAL(LU_NONZERO,(-ONE),Ghimj,1)
     ghinv = ONE/(Direction*H*gam)
     DO i=1,NVAR
       Ghimj(LU_DIAG(i)) = Ghimj(LU_DIAG(i))+ghinv
     END DO

     CALL cb05_sorg_aq_ros_Decomp( Ghimj, Pivot, ising, Ndec )
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

  END SUBROUTINE  cb05_sorg_aq_ros_PrepareMatrix



  SUBROUTINE  cb05_sorg_aq_ros_Decomp( A, Pivot, ising, Ndec )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(INOUT) :: A(LU_NONZERO)

   INTEGER, INTENT(OUT) :: Pivot(NVAR), ising
   INTEGER, INTENT(INOUT) :: Ndec 



CALL decomp_cb05_sorg_aq ( A, ising )
   Pivot(1) = 1
   Ndec = Ndec + 1

  END SUBROUTINE  cb05_sorg_aq_ros_Decomp



  SUBROUTINE  cb05_sorg_aq_ros_Solve( A, Pivot, b, Nsol )



   IMPLICIT NONE

   REAL(kind=dp), INTENT(IN) :: A(LU_NONZERO)
   INTEGER, INTENT(IN) :: Pivot(NVAR)

   INTEGER, INTENT(INOUT) :: nsol 

   REAL(kind=dp), INTENT(INOUT) :: b(NVAR)


   CALL cb05_sorg_aq_KppSolve( A, b )

   Nsol = Nsol+1

  END SUBROUTINE  cb05_sorg_aq_ros_Solve




  SUBROUTINE  cb05_sorg_aq_Ros2 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

 END SUBROUTINE  cb05_sorg_aq_Ros2



  SUBROUTINE  cb05_sorg_aq_Ros3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cb05_sorg_aq_Ros3





  SUBROUTINE  cb05_sorg_aq_Ros4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cb05_sorg_aq_Ros4


  SUBROUTINE  cb05_sorg_aq_Rodas3 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cb05_sorg_aq_Rodas3


  SUBROUTINE  cb05_sorg_aq_Rodas4 (ros_S,ros_A,ros_C,ros_M,ros_E,ros_Alpha,&
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

  END SUBROUTINE  cb05_sorg_aq_Rodas4




END SUBROUTINE  cb05_sorg_aq_Rosenbrock




SUBROUTINE  cb05_sorg_aq_FunTemplate( T, Y, Ydot, RCONST, FIX, Nfun )




   USE cb05_sorg_aq_Parameters




   REAL(kind=dp) :: T, Y(NVAR)
   REAL(kind=dp) :: RCONST(NREACT)
   REAL(kind=dp) :: FIX(NFIX)

   REAL(kind=dp) :: Ydot(NVAR)
   INTEGER :: Nfun









   CALL cb05_sorg_aq_Fun( Y, FIX, RCONST, Ydot )


   Nfun = Nfun+1

END SUBROUTINE  cb05_sorg_aq_FunTemplate



SUBROUTINE  cb05_sorg_aq_JacTemplate( T, Y, Jcb, FIX, Njac, RCONST )




 USE cb05_sorg_aq_Parameters
 
 USE cb05_sorg_aq_Jacobian



    REAL(kind=dp) :: T, Y(NVAR)
    REAL(kind=dp) :: FIX(NFIX)
    REAL(kind=dp) :: RCONST(NREACT)

    INTEGER :: Njac


    REAL(kind=dp) :: Jcb(LU_NONZERO)

    REAL(kind=dp) :: Told





    CALL cb05_sorg_aq_Jac_SP( Y, FIX, RCONST, Jcb )


    Njac = Njac+1

END SUBROUTINE  cb05_sorg_aq_JacTemplate

















SUBROUTINE cb05_sorg_aq_Fun ( V, F, RCT, Vdot )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: Vdot(NVAR)




  REAL(kind=dp) :: A(NREACT)


  A(1) = RCT(1)*V(87)
  A(2) = RCT(2)*V(93)
  A(3) = RCT(3)*V(86)*V(95)
  A(4) = RCT(4)*V(87)*V(93)
  A(5) = RCT(5)*V(87)*V(93)
  A(6) = RCT(6)*V(93)*V(95)
  A(7) = RCT(7)*V(86)*V(87)
  A(8) = RCT(8)*V(86)
  A(9) = RCT(9)*V(86)
  A(10) = RCT(10)*V(45)*F(2)
  A(11) = RCT(11)*V(45)*F(1)
  A(12) = RCT(12)*V(86)*V(96)
  A(13) = RCT(13)*V(86)*V(90)
  A(14) = RCT(14)*V(97)
  A(15) = RCT(15)*V(97)
  A(16) = RCT(16)*V(95)*V(97)
  A(17) = RCT(17)*V(87)*V(97)
  A(18) = RCT(18)*V(87)*V(97)
  A(19) = RCT(19)*V(46)*F(1)
  A(20) = RCT(20)*V(46)*F(1)*F(1)
  A(21) = RCT(21)*V(46)
  A(22) = RCT(22)*V(95)*V(95)
  A(23) = RCT(23)*V(87)*V(95)*F(1)
  A(24) = RCT(24)*V(95)*V(96)
  A(25) = RCT(25)*V(50)
  A(26) = RCT(26)*V(50)*V(96)
  A(27) = RCT(27)*V(50)*V(50)
  A(28) = RCT(28)*V(87)*V(96)
  A(29) = RCT(29)*V(69)*V(96)
  A(30) = RCT(30)*V(90)*V(95)
  A(31) = RCT(31)*V(87)*V(90)
  A(32) = RCT(32)*V(53)
  A(33) = RCT(33)*V(53)*V(96)
  A(34) = RCT(34)*V(90)*V(90)
  A(35) = RCT(35)*V(90)*V(90)*F(1)
  A(36) = RCT(36)*V(60)
  A(37) = RCT(37)*V(60)*V(96)
  A(38) = RCT(38)*V(44)*V(45)
  A(39) = RCT(39)*V(44)*V(96)
  A(40) = RCT(40)*V(93)*V(96)
  A(41) = RCT(41)*V(96)*V(96)
  A(42) = RCT(42)*V(96)*V(96)
  A(43) = RCT(43)*V(90)*V(96)
  A(44) = RCT(44)*V(90)*V(93)
  A(45) = RCT(45)*V(60)*V(93)
  A(46) = RCT(46)*V(93)*V(97)
  A(47) = RCT(47)*V(96)*V(97)
  A(48) = RCT(48)*V(90)*V(97)
  A(49) = RCT(49)*V(86)*V(97)
  A(50) = RCT(50)*V(97)*V(97)
  A(51) = RCT(51)*V(53)
  A(52) = RCT(52)*V(69)
  A(53) = RCT(53)*V(46)
  A(54) = RCT(54)*V(89)*V(95)
  A(55) = RCT(55)*V(78)*V(95)
  A(56) = RCT(56)*V(89)*V(90)
  A(57) = RCT(57)*V(78)*V(90)
  A(58) = RCT(58)*V(89)*V(89)
  A(59) = RCT(59)*V(78)*V(78)
  A(60) = RCT(60)*V(78)*V(89)
  A(61) = RCT(61)*V(84)*V(96)
  A(62) = RCT(62)*V(84)
  A(63) = RCT(63)*V(65)*V(96)
  A(64) = RCT(64)*V(65)
  A(65) = RCT(65)*V(70)*V(96)
  A(66) = RCT(66)*V(48)*V(96)
  A(67) = RCT(67)*V(92)*V(95)
  A(68) = RCT(68)*V(90)*V(92)
  A(69) = RCT(69)*V(92)*V(92)
  A(70) = RCT(70)*V(68)*V(96)
  A(71) = RCT(71)*V(68)
  A(72) = RCT(72)*V(58)*V(96)
  A(73) = RCT(73)*V(85)*V(96)
  A(74) = RCT(74)*V(85)
  A(75) = RCT(75)*V(85)
  A(76) = RCT(76)*V(85)*V(93)
  A(77) = RCT(77)*V(85)*V(97)
  A(78) = RCT(78)*V(85)*V(90)
  A(79) = RCT(79)*V(62)
  A(80) = RCT(80)*V(62)*V(95)
  A(81) = RCT(81)*V(62)*V(90)
  A(82) = RCT(82)*V(51)*V(96)
  A(83) = RCT(83)*V(83)*V(93)
  A(84) = RCT(84)*V(83)*V(96)
  A(85) = RCT(85)*V(83)*V(97)
  A(86) = RCT(86)*V(83)
  A(87) = RCT(87)*V(91)*V(95)
  A(88) = RCT(88)*V(87)*V(91)
  A(89) = RCT(89)*V(41)
  A(90) = RCT(90)*V(41)
  A(91) = RCT(91)*V(90)*V(91)
  A(92) = RCT(92)*V(91)*V(92)
  A(93) = RCT(93)*V(89)*V(91)
  A(94) = RCT(94)*V(91)*V(91)
  A(95) = RCT(95)*V(54)*V(96)
  A(96) = RCT(96)*V(54)
  A(97) = RCT(97)*V(55)*V(96)
  A(98) = RCT(98)*V(82)*V(93)
  A(99) = RCT(99)*V(82)*V(96)
  A(100) = RCT(100)*V(82)*V(97)
  A(101) = RCT(101)*V(82)
  A(102) = RCT(102)*V(94)*V(95)
  A(103) = RCT(103)*V(87)*V(94)
  A(104) = RCT(104)*V(56)
  A(105) = RCT(105)*V(56)
  A(106) = RCT(106)*V(56)*V(96)
  A(107) = RCT(107)*V(90)*V(94)
  A(108) = RCT(108)*V(92)*V(94)
  A(109) = RCT(109)*V(89)*V(94)
  A(110) = RCT(110)*V(94)*V(94)
  A(111) = RCT(111)*V(91)*V(94)
  A(112) = RCT(112)*V(79)*V(96)
  A(113) = RCT(113)*V(72)
  A(114) = 1600*V(72)
  A(115) = RCT(115)*V(72)*V(87)
  A(116) = RCT(116)*V(77)*V(93)
  A(117) = RCT(117)*V(77)*V(96)
  A(118) = RCT(118)*V(77)*V(86)
  A(119) = RCT(119)*V(77)*V(97)
  A(120) = RCT(120)*V(74)*V(93)
  A(121) = RCT(121)*V(74)*V(96)
  A(122) = RCT(122)*V(74)*V(86)
  A(123) = RCT(123)*V(74)*V(97)
  A(124) = RCT(124)*V(76)*V(93)
  A(125) = RCT(125)*V(76)*V(96)
  A(126) = RCT(126)*V(76)*V(86)
  A(127) = RCT(127)*V(76)*V(97)
  A(128) = RCT(128)*V(43)*V(96)
  A(129) = RCT(129)*V(52)*V(95)
  A(130) = 4.2*V(52)
  A(131) = RCT(131)*V(73)*V(96)
  A(132) = RCT(132)*V(73)*V(97)
  A(133) = RCT(133)*V(63)*V(87)
  A(134) = RCT(134)*V(63)*V(90)
  A(135) = RCT(135)*V(71)
  A(136) = RCT(136)*V(71)*V(96)
  A(137) = RCT(137)*V(71)*V(86)
  A(138) = RCT(138)*V(47)*V(96)
  A(139) = RCT(139)*V(66)*V(96)
  A(140) = RCT(140)*V(66)
  A(141) = RCT(141)*V(80)*V(93)
  A(142) = RCT(142)*V(80)*V(96)
  A(143) = RCT(143)*V(80)*V(86)
  A(144) = RCT(144)*V(80)*V(97)
  A(145) = RCT(145)*V(81)*V(96)
  A(146) = RCT(146)*V(81)*V(86)
  A(147) = RCT(147)*V(81)*V(97)
  A(148) = RCT(148)*V(81)
  A(149) = RCT(149)*V(75)*V(93)
  A(150) = RCT(150)*V(75)*V(96)
  A(151) = RCT(151)*V(75)*V(86)
  A(152) = RCT(152)*V(75)*V(97)
  A(153) = RCT(153)*V(40)*V(96)
  A(154) = RCT(154)*V(59)*V(96)
  A(155) = RCT(155)*V(57)*V(96)
  A(156) = RCT(156)*V(80)*V(87)
  A(157) = RCT(157)*V(39)
  A(158) = RCT(158)*V(42)
  A(159) = RCT(159)*V(86)*V(88)
  A(160) = RCT(160)*V(64)*V(64)
  A(161) = RCT(161)*V(64)*V(95)
  A(162) = RCT(162)*V(64)*V(90)
  A(163) = RCT(163)*V(67)*V(96)
  A(164) = RCT(164)*V(67)
  A(165) = RCT(165)*V(48)*V(88)
  A(166) = RCT(166)*V(79)*V(88)
  A(167) = RCT(167)*V(57)*V(88)
  A(168) = RCT(168)*V(74)*V(88)
  A(169) = RCT(169)*V(77)*V(88)
  A(170) = RCT(170)*V(76)*V(88)
  A(171) = RCT(171)*V(80)*V(88)
  A(172) = RCT(172)*V(85)*V(88)
  A(173) = RCT(173)*V(83)*V(88)
  A(174) = RCT(174)*V(82)*V(88)
  A(175) = RCT(175)*V(58)*V(88)
  A(176) = RCT(176)*V(59)*V(88)
  A(177) = RCT(177)*V(61)*V(96)
  A(178) = RCT(178)*V(49)*V(86)
  A(179) = RCT(179)*V(49)*V(96)
  A(180) = RCT(180)*V(49)*V(60)
  A(181) = RCT(181)*V(12)*V(96)
  A(182) = RCT(182)*V(15)*V(96)
  A(183) = RCT(183)*V(18)*V(96)
  A(184) = RCT(184)*V(23)*V(96)
  A(185) = RCT(185)*V(23)*V(86)
  A(186) = RCT(186)*V(29)*V(96)
  A(187) = RCT(187)*V(29)*V(86)
  A(188) = RCT(188)*V(29)*V(97)
  A(189) = RCT(189)*V(32)*V(96)
  A(190) = RCT(190)*V(34)*V(96)
  A(191) = RCT(191)*V(37)*V(96)


  Vdot(1) = 0.071*A(128)
  Vdot(2) = 0.138*A(128)
  Vdot(3) = A(131)+A(132)
  Vdot(4) = 0.038*A(138)
  Vdot(5) = 0.167*A(138)
  Vdot(6) = 0.232*A(141)+0.232*A(142)+0.232*A(143)+0.232*A(144)
  Vdot(7) = 0.0228*A(141)+0.0288*A(142)+0.0288*A(143)+0.0288*A(144)
  Vdot(8) = A(153)
  Vdot(9) = A(153)
  Vdot(10) = A(149)+A(150)+A(151)+A(152)
  Vdot(11) = A(181)
  Vdot(12) = -A(181)
  Vdot(13) = 0.239*A(182)
  Vdot(14) = 0.363*A(182)
  Vdot(15) = -A(182)
  Vdot(16) = 0.045*A(183)
  Vdot(17) = 0.149*A(183)
  Vdot(18) = -A(183)
  Vdot(19) = 0.038*A(184)
  Vdot(20) = 0.326*A(184)
  Vdot(21) = 0.125*A(185)
  Vdot(22) = 0.102*A(185)
  Vdot(23) = -A(184)-A(185)
  Vdot(24) = 0.13*A(186)
  Vdot(25) = 0.0406*A(186)
  Vdot(26) = 0.026*A(187)
  Vdot(27) = 0.485*A(187)
  Vdot(28) = A(188)
  Vdot(29) = -A(186)-A(187)-A(188)
  Vdot(30) = 0.091*A(189)
  Vdot(31) = 0.367*A(189)
  Vdot(32) = -A(189)
  Vdot(33) = 1.173*A(190)
  Vdot(34) = -A(190)
  Vdot(35) = 0.156*A(191)
  Vdot(36) = 0.777*A(191)
  Vdot(37) = -A(191)
  Vdot(38) = A(178)+A(179)+A(180)
  Vdot(39) = -A(157)+0.3*A(160)
  Vdot(40) = -A(153)
  Vdot(41) = A(88)-A(89)-A(90)
  Vdot(42) = -A(158)+A(162)
  Vdot(43) = -A(128)
  Vdot(44) = -A(38)-A(39)
  Vdot(45) = A(9)-A(10)-A(11)-A(38)
  Vdot(46) = A(18)-A(19)-A(20)-A(21)-A(53)
  Vdot(47) = -A(138)
  Vdot(48) = -A(66)-A(165)
  Vdot(49) = -A(178)-A(179)-A(180)
  Vdot(50) = 2*A(23)+A(24)-A(25)-A(26)-2*A(27)
  Vdot(51) = A(80)-A(82)+0.37*A(122)
  Vdot(52) = 0.56*A(128)-A(129)-A(130)+0.3*A(138)
  Vdot(53) = A(31)-A(32)-A(33)-A(51)
  Vdot(54) = 0.8*A(91)-A(95)-A(96)+0.8*A(107)
  Vdot(55) = 0.2*A(91)+0.1*A(92)+0.1*A(93)-A(97)+0.2*A(107)+0.1*A(108)+0.1*A(109)
  Vdot(56) = A(103)-A(104)-A(105)-A(106)
  Vdot(57) = -A(155)-A(167)
  Vdot(58) = 0.63*A(69)-A(72)-A(175)
  Vdot(59) = -A(154)-A(176)
  Vdot(60) = A(34)+A(35)-A(36)-A(37)+A(42)-A(45)-A(180)
  Vdot(61) = A(165)+A(166)+A(167)+0.3*A(170)+0.15*A(171)+A(172)+A(173)+A(174)+A(175)+A(176)-A(177)
  Vdot(62) = A(78)-A(79)-A(80)-A(81)
  Vdot(63) = 0.4*A(131)+A(132)-A(133)-A(134)
  Vdot(64) = A(159)-2*A(160)-A(161)-A(162)
  Vdot(65) = A(56)+A(57)-A(63)-A(64)
  Vdot(66) = 0.2*A(137)+0.8*A(138)-A(139)-A(140)+0.168*A(145)+0.85*A(146)
  Vdot(67) = -A(163)-A(164)+A(168)+A(169)+0.7*A(170)+0.85*A(171)
  Vdot(68) = A(68)-A(70)-A(71)+A(81)
  Vdot(69) = 2*A(19)+2*A(20)+A(28)-A(29)+A(48)-A(52)+A(61)+A(77)+A(85)+A(100)+A(132)+0.15*A(147)
  Vdot(70) = -A(65)+A(73)+A(74)+A(75)+A(76)+A(77)+A(86)+A(101)+0.2*A(116)+0.33*A(118)+A(120)+0.63*A(122)+0.1*A(124)+0.25&
               &*A(126)+A(135)+2*A(136)+0.69*A(137)+A(140)+0.066*A(143)+0.334*A(145)+0.225*A(146)+0.643*A(147)+0.333*A(148)&
               &+0.001*A(151)+A(163)+A(164)+A(172)
  Vdot(71) = 0.9*A(129)+0.3*A(131)-A(135)-A(136)-A(137)
  Vdot(72) = 0.76*A(112)-0.98*A(113)-A(114)-A(115)+0.76*A(166)
  Vdot(73) = 0.36*A(128)+A(130)-A(131)-A(132)+A(134)+0.2*A(138)
  Vdot(74) = -A(120)-A(121)-A(122)-A(123)-A(168)
  Vdot(75) = -A(149)-A(150)-A(151)-A(152)
  Vdot(76) = -A(124)-A(125)-A(126)-A(127)-A(170)
  Vdot(77) = -A(116)-A(117)-A(118)-A(119)-A(169)+0.3*A(170)
  Vdot(78) = -A(55)-A(57)-2*A(59)-A(60)+0.13*A(112)+0.04*A(113)+0.01*A(116)+0.09*A(119)+0.088*A(142)+0.25*A(150)+0.18&
               &*A(151)+0.25*A(152)+0.009*A(155)+0.13*A(166)+0.009*A(167)
  Vdot(79) = -0.66*A(61)-0.66*A(62)-1.11*A(112)-2.1*A(113)+0.2*A(116)-0.7*A(117)-A(118)-A(119)+0.1*A(124)+1.1*A(138)&
               &+0.35*A(143)+1.565*A(145)+0.36*A(146)+1.282*A(147)+0.832*A(148)+5.12*A(149)+1.66*A(150)+7*A(151)+2.4*A(156)&
               &-1.11*A(166)-A(169)+0.3*A(170)
  Vdot(80) = -A(141)-A(142)-A(143)-A(144)-A(156)-A(171)
  Vdot(81) = 0.75*A(141)+0.912*A(142)+0.65*A(143)+0.2*A(144)-A(145)-A(146)-A(147)-A(148)+0.2*A(156)+A(171)
  Vdot(82) = 0.33*A(61)+0.33*A(62)+0.5*A(63)+0.5*A(64)-A(98)-A(99)-A(100)-A(101)+0.05*A(112)+0.5*A(113)+0.3*A(116)+0.62&
               &*A(117)+0.32*A(118)+0.56*A(119)+0.22*A(121)+0.66*A(124)+0.7*A(125)+0.35*A(126)+0.64*A(127)+0.03*A(137)+0.15&
               &*A(143)+0.8*A(144)+0.12*A(145)+0.357*A(147)+0.47*A(150)+0.21*A(151)+0.47*A(152)+0.05*A(154)+0.8*A(156)+0.05&
               &*A(166)+0.67*A(169)+0.55*A(170)-A(174)
  Vdot(83) = 0.33*A(61)+0.33*A(62)+0.5*A(63)+0.5*A(64)-A(83)-A(84)-A(85)-A(86)+A(102)+A(106)+0.9*A(108)+0.9*A(109)+2&
               &*A(110)+A(111)+0.06*A(112)+0.6*A(113)+0.2*A(116)+0.33*A(117)+0.18*A(118)+0.35*A(119)+1.24*A(124)+1.3*A(125)&
               &+0.65*A(126)+1.18*A(127)+0.252*A(145)+0.02*A(146)+0.067*A(148)+0.9*A(154)+0.991*A(155)+0.06*A(166)+0.991&
               &*A(167)+0.33*A(169)+0.45*A(170)-A(173)+A(176)
  Vdot(84) = A(55)-A(61)-A(62)+A(115)+0.1*A(129)+A(133)+0.8*A(144)+0.85*A(147)+0.53*A(152)+0.8*A(156)
  Vdot(85) = 0.33*A(61)+0.33*A(62)+A(67)+1.37*A(69)+A(71)+A(72)-A(73)-A(74)-A(75)-A(76)-A(77)-A(78)+A(79)+A(92)+0.1&
               &*A(108)+0.2*A(116)+0.8*A(117)+0.74*A(118)+A(119)+A(120)+1.56*A(121)+A(122)+2*A(123)+0.25*A(126)+A(136)+0.7&
               &*A(137)+0.5*A(141)+0.629*A(142)+0.6*A(143)+0.167*A(145)+0.15*A(146)+0.282*A(147)+0.9*A(148)+0.28*A(150)+0.24&
               &*A(151)+0.1*A(154)+A(168)-A(172)+A(175)
  Vdot(86) = A(2)-A(3)-A(7)-A(8)-A(9)-A(12)-A(13)-A(49)+0.2*A(91)+0.2*A(107)-A(118)-A(122)-A(126)-A(137)-A(143)-A(146)&
               &-A(151)-A(159)-A(178)
  Vdot(87) = -A(1)+A(3)-A(4)-A(5)+A(6)-A(7)+A(14)+2*A(16)-A(18)+A(21)+2*A(22)-A(23)+A(26)+A(27)-A(28)+A(30)-A(31)+A(32)&
               &+A(33)+A(46)+A(47)+A(49)+2*A(50)+0.61*A(51)+A(52)+A(53)+A(54)+A(62)+A(67)+A(80)+A(87)-A(88)+A(89)+A(90)&
               &+A(102)-A(103)+A(104)+A(105)+A(106)-A(115)+A(119)+A(123)+A(127)+0.9*A(129)-A(133)+0.2*A(144)+0.47*A(152)&
               &-A(156)+A(161)
  Vdot(88) = 2*A(157)+A(158)-A(159)+1.4*A(160)+A(161)+A(163)+A(164)-A(165)-A(166)-A(167)-A(168)-A(169)-A(170)-A(171)&
               &-A(172)-A(173)-A(174)-A(175)-A(176)+A(177)
  Vdot(89) = -A(54)-A(56)-2*A(58)-A(60)+A(63)+0.3*A(70)-A(93)+A(102)+0.9*A(108)-A(109)+2*A(110)+A(111)+0.87*A(112)+0.96&
               &*A(113)+0.2*A(116)+0.8*A(117)+0.22*A(118)+0.91*A(119)+0.7*A(120)+A(121)+A(123)+0.1*A(124)+A(125)+0.08*A(128)&
               &+0.6*A(131)+A(136)+0.03*A(137)+0.5*A(138)+A(139)+0.25*A(141)+0.991*A(142)+0.2*A(143)+A(144)+0.713*A(145)&
               &+0.064*A(146)+0.075*A(147)+0.7*A(148)+1.25*A(150)+0.76*A(151)+1.03*A(152)+0.1*A(154)+0.991*A(155)+A(156)&
               &+0.87*A(166)+0.991*A(167)+2*A(168)+2*A(169)+1.7*A(170)+A(171)
  Vdot(90) = A(12)-A(13)-A(30)-A(31)+A(32)-2*A(34)-2*A(35)+A(37)+A(38)+A(39)+A(40)-A(43)-A(44)+A(45)+A(47)-A(48)+0.61&
               &*A(51)-A(56)-A(57)+A(61)+A(62)+A(64)+A(65)+A(67)-A(68)+0.74*A(69)+0.3*A(70)+A(71)+A(72)+A(73)+2*A(74)+A(76)&
               &+A(77)-A(78)+A(79)+A(80)-A(81)+A(82)+A(86)-A(91)+0.9*A(92)+A(101)+A(102)-A(107)+A(108)+2*A(110)+A(111)+0.11&
               &*A(112)+0.94*A(113)+A(114)+0.3*A(116)+0.95*A(117)+0.44*A(118)+1.7*A(120)+A(121)+0.13*A(122)+0.1*A(124)&
               &+A(125)+0.5*A(126)+A(127)+0.44*A(128)+0.9*A(129)+A(130)+0.6*A(131)-A(134)+A(135)+2*A(136)+0.76*A(137)+0.7&
               &*A(138)+A(140)+0.25*A(141)+0.912*A(142)+0.066*A(143)+0.8*A(144)+0.503*A(145)+0.154*A(146)+0.925*A(147)+1.033&
               &*A(148)+0.75*A(150)+0.07*A(151)+0.28*A(152)+A(153)+A(154)+A(155)+0.8*A(156)-A(162)+A(164)+0.11*A(166)+A(167)&
               &+A(168)+A(169)+A(170)+A(171)+A(172)+A(175)+A(176)
  Vdot(91) = A(83)+A(84)+A(85)-A(87)-A(88)+A(89)+A(90)-A(91)-A(92)-A(93)-2*A(94)+A(95)-A(111)+A(135)+A(136)+0.62*A(137)&
               &+A(139)+A(140)+0.21*A(145)+0.114*A(146)+0.967*A(148)+A(173)
  Vdot(92) = A(66)-A(67)-A(68)-2*A(69)+0.7*A(70)+A(86)+A(87)-0.1*A(92)+0.9*A(93)+2*A(94)+A(96)+A(97)+A(101)-A(108)&
               &+A(111)+A(165)
  Vdot(93) = A(1)-A(2)-A(4)-A(5)-A(6)+A(8)+A(10)+A(14)-A(40)+A(41)-A(44)-A(45)-A(46)-A(76)-A(83)-A(98)-A(116)-A(120)&
               &-A(124)+0.5*A(126)-A(141)-A(149)
  Vdot(94) = A(98)+A(99)+A(100)-A(102)-A(103)+A(104)+A(105)-A(107)-A(108)-A(109)-2*A(110)-A(111)+0.25*A(141)+0.2*A(143)&
               &+0.25*A(145)+0.075*A(147)+0.39*A(151)+A(174)
  Vdot(95) = A(1)-A(3)+A(4)-A(6)+A(15)-A(16)+A(17)-2*A(22)-A(23)-A(24)+A(25)+A(27)-A(30)-A(54)-A(55)-A(67)-A(80)-A(87)&
               &-A(102)-A(129)+0.2*A(156)-A(161)
  Vdot(96) = 2*A(11)-A(12)+A(13)-A(24)+A(25)-A(26)-A(28)-A(29)+A(30)-A(33)+2*A(36)-A(37)+A(38)-A(39)-A(40)-2*A(41)-2&
               &*A(42)-A(43)+A(44)+A(45)-A(47)+0.39*A(51)+A(52)-A(61)-A(63)+A(64)-A(65)-A(66)-A(70)+A(71)-A(72)-A(73)+A(76)&
               &-A(82)+A(83)-A(84)-A(95)+A(96)-A(97)+A(98)-A(99)-A(106)-A(112)+0.1*A(116)-A(117)+0.1*A(118)+0.3*A(120)&
               &-A(121)+0.13*A(122)-A(125)+0.5*A(126)-A(128)-A(131)-A(136)+0.08*A(137)-A(138)-A(139)-A(142)+0.266*A(143)&
               &-A(145)+0.268*A(146)+0.15*A(149)-A(150)+0.57*A(151)-A(153)-A(154)-A(155)+A(158)-A(163)-A(177)-A(179)
  Vdot(97) = A(5)+A(7)-A(14)-A(15)-A(16)-A(17)-A(18)+A(21)+A(29)-A(46)-A(47)-A(48)-A(49)-2*A(50)+0.39*A(51)+A(53)-A(77)&
               &-A(85)-A(100)-A(119)-A(123)-A(127)-A(132)-A(144)-A(147)-A(152)
      
END SUBROUTINE cb05_sorg_aq_Fun
















SUBROUTINE cb05_sorg_aq_IRRFun ( V, F, RCT, IRR )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: IRR(NREACT)



  IRR(1) = RCT(1)*V(87)
  IRR(2) = RCT(2)*V(93)
  IRR(3) = RCT(3)*V(86)*V(95)
  IRR(4) = RCT(4)*V(87)*V(93)
  IRR(5) = RCT(5)*V(87)*V(93)
  IRR(6) = RCT(6)*V(93)*V(95)
  IRR(7) = RCT(7)*V(86)*V(87)
  IRR(8) = RCT(8)*V(86)
  IRR(9) = RCT(9)*V(86)
  IRR(10) = RCT(10)*V(45)*F(2)
  IRR(11) = RCT(11)*V(45)*F(1)
  IRR(12) = RCT(12)*V(86)*V(96)
  IRR(13) = RCT(13)*V(86)*V(90)
  IRR(14) = RCT(14)*V(97)
  IRR(15) = RCT(15)*V(97)
  IRR(16) = RCT(16)*V(95)*V(97)
  IRR(17) = RCT(17)*V(87)*V(97)
  IRR(18) = RCT(18)*V(87)*V(97)
  IRR(19) = RCT(19)*V(46)*F(1)
  IRR(20) = RCT(20)*V(46)*F(1)*F(1)
  IRR(21) = RCT(21)*V(46)
  IRR(22) = RCT(22)*V(95)*V(95)
  IRR(23) = RCT(23)*V(87)*V(95)*F(1)
  IRR(24) = RCT(24)*V(95)*V(96)
  IRR(25) = RCT(25)*V(50)
  IRR(26) = RCT(26)*V(50)*V(96)
  IRR(27) = RCT(27)*V(50)*V(50)
  IRR(28) = RCT(28)*V(87)*V(96)
  IRR(29) = RCT(29)*V(69)*V(96)
  IRR(30) = RCT(30)*V(90)*V(95)
  IRR(31) = RCT(31)*V(87)*V(90)
  IRR(32) = RCT(32)*V(53)
  IRR(33) = RCT(33)*V(53)*V(96)
  IRR(34) = RCT(34)*V(90)*V(90)
  IRR(35) = RCT(35)*V(90)*V(90)*F(1)
  IRR(36) = RCT(36)*V(60)
  IRR(37) = RCT(37)*V(60)*V(96)
  IRR(38) = RCT(38)*V(44)*V(45)
  IRR(39) = RCT(39)*V(44)*V(96)
  IRR(40) = RCT(40)*V(93)*V(96)
  IRR(41) = RCT(41)*V(96)*V(96)
  IRR(42) = RCT(42)*V(96)*V(96)
  IRR(43) = RCT(43)*V(90)*V(96)
  IRR(44) = RCT(44)*V(90)*V(93)
  IRR(45) = RCT(45)*V(60)*V(93)
  IRR(46) = RCT(46)*V(93)*V(97)
  IRR(47) = RCT(47)*V(96)*V(97)
  IRR(48) = RCT(48)*V(90)*V(97)
  IRR(49) = RCT(49)*V(86)*V(97)
  IRR(50) = RCT(50)*V(97)*V(97)
  IRR(51) = RCT(51)*V(53)
  IRR(52) = RCT(52)*V(69)
  IRR(53) = RCT(53)*V(46)
  IRR(54) = RCT(54)*V(89)*V(95)
  IRR(55) = RCT(55)*V(78)*V(95)
  IRR(56) = RCT(56)*V(89)*V(90)
  IRR(57) = RCT(57)*V(78)*V(90)
  IRR(58) = RCT(58)*V(89)*V(89)
  IRR(59) = RCT(59)*V(78)*V(78)
  IRR(60) = RCT(60)*V(78)*V(89)
  IRR(61) = RCT(61)*V(84)*V(96)
  IRR(62) = RCT(62)*V(84)
  IRR(63) = RCT(63)*V(65)*V(96)
  IRR(64) = RCT(64)*V(65)
  IRR(65) = RCT(65)*V(70)*V(96)
  IRR(66) = RCT(66)*V(48)*V(96)
  IRR(67) = RCT(67)*V(92)*V(95)
  IRR(68) = RCT(68)*V(90)*V(92)
  IRR(69) = RCT(69)*V(92)*V(92)
  IRR(70) = RCT(70)*V(68)*V(96)
  IRR(71) = RCT(71)*V(68)
  IRR(72) = RCT(72)*V(58)*V(96)
  IRR(73) = RCT(73)*V(85)*V(96)
  IRR(74) = RCT(74)*V(85)
  IRR(75) = RCT(75)*V(85)
  IRR(76) = RCT(76)*V(85)*V(93)
  IRR(77) = RCT(77)*V(85)*V(97)
  IRR(78) = RCT(78)*V(85)*V(90)
  IRR(79) = RCT(79)*V(62)
  IRR(80) = RCT(80)*V(62)*V(95)
  IRR(81) = RCT(81)*V(62)*V(90)
  IRR(82) = RCT(82)*V(51)*V(96)
  IRR(83) = RCT(83)*V(83)*V(93)
  IRR(84) = RCT(84)*V(83)*V(96)
  IRR(85) = RCT(85)*V(83)*V(97)
  IRR(86) = RCT(86)*V(83)
  IRR(87) = RCT(87)*V(91)*V(95)
  IRR(88) = RCT(88)*V(87)*V(91)
  IRR(89) = RCT(89)*V(41)
  IRR(90) = RCT(90)*V(41)
  IRR(91) = RCT(91)*V(90)*V(91)
  IRR(92) = RCT(92)*V(91)*V(92)
  IRR(93) = RCT(93)*V(89)*V(91)
  IRR(94) = RCT(94)*V(91)*V(91)
  IRR(95) = RCT(95)*V(54)*V(96)
  IRR(96) = RCT(96)*V(54)
  IRR(97) = RCT(97)*V(55)*V(96)
  IRR(98) = RCT(98)*V(82)*V(93)
  IRR(99) = RCT(99)*V(82)*V(96)
  IRR(100) = RCT(100)*V(82)*V(97)
  IRR(101) = RCT(101)*V(82)
  IRR(102) = RCT(102)*V(94)*V(95)
  IRR(103) = RCT(103)*V(87)*V(94)
  IRR(104) = RCT(104)*V(56)
  IRR(105) = RCT(105)*V(56)
  IRR(106) = RCT(106)*V(56)*V(96)
  IRR(107) = RCT(107)*V(90)*V(94)
  IRR(108) = RCT(108)*V(92)*V(94)
  IRR(109) = RCT(109)*V(89)*V(94)
  IRR(110) = RCT(110)*V(94)*V(94)
  IRR(111) = RCT(111)*V(91)*V(94)
  IRR(112) = RCT(112)*V(79)*V(96)
  IRR(113) = RCT(113)*V(72)
  IRR(114) = 1600*V(72)
  IRR(115) = RCT(115)*V(72)*V(87)
  IRR(116) = RCT(116)*V(77)*V(93)
  IRR(117) = RCT(117)*V(77)*V(96)
  IRR(118) = RCT(118)*V(77)*V(86)
  IRR(119) = RCT(119)*V(77)*V(97)
  IRR(120) = RCT(120)*V(74)*V(93)
  IRR(121) = RCT(121)*V(74)*V(96)
  IRR(122) = RCT(122)*V(74)*V(86)
  IRR(123) = RCT(123)*V(74)*V(97)
  IRR(124) = RCT(124)*V(76)*V(93)
  IRR(125) = RCT(125)*V(76)*V(96)
  IRR(126) = RCT(126)*V(76)*V(86)
  IRR(127) = RCT(127)*V(76)*V(97)
  IRR(128) = RCT(128)*V(43)*V(96)
  IRR(129) = RCT(129)*V(52)*V(95)
  IRR(130) = 4.2*V(52)
  IRR(131) = RCT(131)*V(73)*V(96)
  IRR(132) = RCT(132)*V(73)*V(97)
  IRR(133) = RCT(133)*V(63)*V(87)
  IRR(134) = RCT(134)*V(63)*V(90)
  IRR(135) = RCT(135)*V(71)
  IRR(136) = RCT(136)*V(71)*V(96)
  IRR(137) = RCT(137)*V(71)*V(86)
  IRR(138) = RCT(138)*V(47)*V(96)
  IRR(139) = RCT(139)*V(66)*V(96)
  IRR(140) = RCT(140)*V(66)
  IRR(141) = RCT(141)*V(80)*V(93)
  IRR(142) = RCT(142)*V(80)*V(96)
  IRR(143) = RCT(143)*V(80)*V(86)
  IRR(144) = RCT(144)*V(80)*V(97)
  IRR(145) = RCT(145)*V(81)*V(96)
  IRR(146) = RCT(146)*V(81)*V(86)
  IRR(147) = RCT(147)*V(81)*V(97)
  IRR(148) = RCT(148)*V(81)
  IRR(149) = RCT(149)*V(75)*V(93)
  IRR(150) = RCT(150)*V(75)*V(96)
  IRR(151) = RCT(151)*V(75)*V(86)
  IRR(152) = RCT(152)*V(75)*V(97)
  IRR(153) = RCT(153)*V(40)*V(96)
  IRR(154) = RCT(154)*V(59)*V(96)
  IRR(155) = RCT(155)*V(57)*V(96)
  IRR(156) = RCT(156)*V(80)*V(87)
  IRR(157) = RCT(157)*V(39)
  IRR(158) = RCT(158)*V(42)
  IRR(159) = RCT(159)*V(86)*V(88)
  IRR(160) = RCT(160)*V(64)*V(64)
  IRR(161) = RCT(161)*V(64)*V(95)
  IRR(162) = RCT(162)*V(64)*V(90)
  IRR(163) = RCT(163)*V(67)*V(96)
  IRR(164) = RCT(164)*V(67)
  IRR(165) = RCT(165)*V(48)*V(88)
  IRR(166) = RCT(166)*V(79)*V(88)
  IRR(167) = RCT(167)*V(57)*V(88)
  IRR(168) = RCT(168)*V(74)*V(88)
  IRR(169) = RCT(169)*V(77)*V(88)
  IRR(170) = RCT(170)*V(76)*V(88)
  IRR(171) = RCT(171)*V(80)*V(88)
  IRR(172) = RCT(172)*V(85)*V(88)
  IRR(173) = RCT(173)*V(83)*V(88)
  IRR(174) = RCT(174)*V(82)*V(88)
  IRR(175) = RCT(175)*V(58)*V(88)
  IRR(176) = RCT(176)*V(59)*V(88)
  IRR(177) = RCT(177)*V(61)*V(96)
  IRR(178) = RCT(178)*V(49)*V(86)
  IRR(179) = RCT(179)*V(49)*V(96)
  IRR(180) = RCT(180)*V(49)*V(60)
  IRR(181) = RCT(181)*V(12)*V(96)
  IRR(182) = RCT(182)*V(15)*V(96)
  IRR(183) = RCT(183)*V(18)*V(96)
  IRR(184) = RCT(184)*V(23)*V(96)
  IRR(185) = RCT(185)*V(23)*V(86)
  IRR(186) = RCT(186)*V(29)*V(96)
  IRR(187) = RCT(187)*V(29)*V(86)
  IRR(188) = RCT(188)*V(29)*V(97)
  IRR(189) = RCT(189)*V(32)*V(96)
  IRR(190) = RCT(190)*V(34)*V(96)
  IRR(191) = RCT(191)*V(37)*V(96)
      
END SUBROUTINE cb05_sorg_aq_IRRFun
















SUBROUTINE cb05_sorg_aq_Jac_SP ( V, F, RCT, JVS )


  REAL(kind=dp) :: V(NVAR)

  REAL(kind=dp) :: F(NFIX)

  REAL(kind=dp) :: RCT(NREACT)

  REAL(kind=dp) :: JVS(LU_NONZERO)




  REAL(kind=dp) :: B(336)


  B(1) = RCT(1)

  B(2) = RCT(2)

  B(3) = RCT(3)*V(95)

  B(4) = RCT(3)*V(86)

  B(5) = RCT(4)*V(93)

  B(6) = RCT(4)*V(87)

  B(7) = RCT(5)*V(93)

  B(8) = RCT(5)*V(87)

  B(9) = RCT(6)*V(95)

  B(10) = RCT(6)*V(93)

  B(11) = RCT(7)*V(87)

  B(12) = RCT(7)*V(86)

  B(13) = RCT(8)

  B(14) = RCT(9)

  B(15) = RCT(10)*F(2)

  B(17) = RCT(11)*F(1)

  B(19) = RCT(12)*V(96)

  B(20) = RCT(12)*V(86)

  B(21) = RCT(13)*V(90)

  B(22) = RCT(13)*V(86)

  B(23) = RCT(14)

  B(24) = RCT(15)

  B(25) = RCT(16)*V(97)

  B(26) = RCT(16)*V(95)

  B(27) = RCT(17)*V(97)

  B(28) = RCT(17)*V(87)

  B(29) = RCT(18)*V(97)

  B(30) = RCT(18)*V(87)

  B(31) = RCT(19)*F(1)

  B(33) = RCT(20)*F(1)*F(1)

  B(35) = RCT(21)

  B(36) = RCT(22)*2*V(95)

  B(37) = RCT(23)*V(95)*F(1)

  B(38) = RCT(23)*V(87)*F(1)

  B(40) = RCT(24)*V(96)

  B(41) = RCT(24)*V(95)

  B(42) = RCT(25)

  B(43) = RCT(26)*V(96)

  B(44) = RCT(26)*V(50)

  B(45) = RCT(27)*2*V(50)

  B(46) = RCT(28)*V(96)

  B(47) = RCT(28)*V(87)

  B(48) = RCT(29)*V(96)

  B(49) = RCT(29)*V(69)

  B(50) = RCT(30)*V(95)

  B(51) = RCT(30)*V(90)

  B(52) = RCT(31)*V(90)

  B(53) = RCT(31)*V(87)

  B(54) = RCT(32)

  B(55) = RCT(33)*V(96)

  B(56) = RCT(33)*V(53)

  B(57) = RCT(34)*2*V(90)

  B(58) = RCT(35)*2*V(90)*F(1)

  B(60) = RCT(36)

  B(61) = RCT(37)*V(96)

  B(62) = RCT(37)*V(60)

  B(63) = RCT(38)*V(45)

  B(64) = RCT(38)*V(44)

  B(65) = RCT(39)*V(96)

  B(66) = RCT(39)*V(44)

  B(67) = RCT(40)*V(96)

  B(68) = RCT(40)*V(93)

  B(69) = RCT(41)*2*V(96)

  B(70) = RCT(42)*2*V(96)

  B(71) = RCT(43)*V(96)

  B(72) = RCT(43)*V(90)

  B(73) = RCT(44)*V(93)

  B(74) = RCT(44)*V(90)

  B(75) = RCT(45)*V(93)

  B(76) = RCT(45)*V(60)

  B(77) = RCT(46)*V(97)

  B(78) = RCT(46)*V(93)

  B(79) = RCT(47)*V(97)

  B(80) = RCT(47)*V(96)

  B(81) = RCT(48)*V(97)

  B(82) = RCT(48)*V(90)

  B(83) = RCT(49)*V(97)

  B(84) = RCT(49)*V(86)

  B(85) = RCT(50)*2*V(97)

  B(86) = RCT(51)

  B(87) = RCT(52)

  B(88) = RCT(53)

  B(89) = RCT(54)*V(95)

  B(90) = RCT(54)*V(89)

  B(91) = RCT(55)*V(95)

  B(92) = RCT(55)*V(78)

  B(93) = RCT(56)*V(90)

  B(94) = RCT(56)*V(89)

  B(95) = RCT(57)*V(90)

  B(96) = RCT(57)*V(78)

  B(97) = RCT(58)*2*V(89)

  B(98) = RCT(59)*2*V(78)

  B(99) = RCT(60)*V(89)

  B(100) = RCT(60)*V(78)

  B(101) = RCT(61)*V(96)

  B(102) = RCT(61)*V(84)

  B(103) = RCT(62)

  B(104) = RCT(63)*V(96)

  B(105) = RCT(63)*V(65)

  B(106) = RCT(64)

  B(107) = RCT(65)*V(96)

  B(108) = RCT(65)*V(70)

  B(109) = RCT(66)*V(96)

  B(110) = RCT(66)*V(48)

  B(111) = RCT(67)*V(95)

  B(112) = RCT(67)*V(92)

  B(113) = RCT(68)*V(92)

  B(114) = RCT(68)*V(90)

  B(115) = RCT(69)*2*V(92)

  B(116) = RCT(70)*V(96)

  B(117) = RCT(70)*V(68)

  B(118) = RCT(71)

  B(119) = RCT(72)*V(96)

  B(120) = RCT(72)*V(58)

  B(121) = RCT(73)*V(96)

  B(122) = RCT(73)*V(85)

  B(123) = RCT(74)

  B(124) = RCT(75)

  B(125) = RCT(76)*V(93)

  B(126) = RCT(76)*V(85)

  B(127) = RCT(77)*V(97)

  B(128) = RCT(77)*V(85)

  B(129) = RCT(78)*V(90)

  B(130) = RCT(78)*V(85)

  B(131) = RCT(79)

  B(132) = RCT(80)*V(95)

  B(133) = RCT(80)*V(62)

  B(134) = RCT(81)*V(90)

  B(135) = RCT(81)*V(62)

  B(136) = RCT(82)*V(96)

  B(137) = RCT(82)*V(51)

  B(138) = RCT(83)*V(93)

  B(139) = RCT(83)*V(83)

  B(140) = RCT(84)*V(96)

  B(141) = RCT(84)*V(83)

  B(142) = RCT(85)*V(97)

  B(143) = RCT(85)*V(83)

  B(144) = RCT(86)

  B(145) = RCT(87)*V(95)

  B(146) = RCT(87)*V(91)

  B(147) = RCT(88)*V(91)

  B(148) = RCT(88)*V(87)

  B(149) = RCT(89)

  B(150) = RCT(90)

  B(151) = RCT(91)*V(91)

  B(152) = RCT(91)*V(90)

  B(153) = RCT(92)*V(92)

  B(154) = RCT(92)*V(91)

  B(155) = RCT(93)*V(91)

  B(156) = RCT(93)*V(89)

  B(157) = RCT(94)*2*V(91)

  B(158) = RCT(95)*V(96)

  B(159) = RCT(95)*V(54)

  B(160) = RCT(96)

  B(161) = RCT(97)*V(96)

  B(162) = RCT(97)*V(55)

  B(163) = RCT(98)*V(93)

  B(164) = RCT(98)*V(82)

  B(165) = RCT(99)*V(96)

  B(166) = RCT(99)*V(82)

  B(167) = RCT(100)*V(97)

  B(168) = RCT(100)*V(82)

  B(169) = RCT(101)

  B(170) = RCT(102)*V(95)

  B(171) = RCT(102)*V(94)

  B(172) = RCT(103)*V(94)

  B(173) = RCT(103)*V(87)

  B(174) = RCT(104)

  B(175) = RCT(105)

  B(176) = RCT(106)*V(96)

  B(177) = RCT(106)*V(56)

  B(178) = RCT(107)*V(94)

  B(179) = RCT(107)*V(90)

  B(180) = RCT(108)*V(94)

  B(181) = RCT(108)*V(92)

  B(182) = RCT(109)*V(94)

  B(183) = RCT(109)*V(89)

  B(184) = RCT(110)*2*V(94)

  B(185) = RCT(111)*V(94)

  B(186) = RCT(111)*V(91)

  B(187) = RCT(112)*V(96)

  B(188) = RCT(112)*V(79)

  B(189) = RCT(113)

  B(190) = 1600

  B(191) = RCT(115)*V(87)

  B(192) = RCT(115)*V(72)

  B(193) = RCT(116)*V(93)

  B(194) = RCT(116)*V(77)

  B(195) = RCT(117)*V(96)

  B(196) = RCT(117)*V(77)

  B(197) = RCT(118)*V(86)

  B(198) = RCT(118)*V(77)

  B(199) = RCT(119)*V(97)

  B(200) = RCT(119)*V(77)

  B(201) = RCT(120)*V(93)

  B(202) = RCT(120)*V(74)

  B(203) = RCT(121)*V(96)

  B(204) = RCT(121)*V(74)

  B(205) = RCT(122)*V(86)

  B(206) = RCT(122)*V(74)

  B(207) = RCT(123)*V(97)

  B(208) = RCT(123)*V(74)

  B(209) = RCT(124)*V(93)

  B(210) = RCT(124)*V(76)

  B(211) = RCT(125)*V(96)

  B(212) = RCT(125)*V(76)

  B(213) = RCT(126)*V(86)

  B(214) = RCT(126)*V(76)

  B(215) = RCT(127)*V(97)

  B(216) = RCT(127)*V(76)

  B(217) = RCT(128)*V(96)

  B(218) = RCT(128)*V(43)

  B(219) = RCT(129)*V(95)

  B(220) = RCT(129)*V(52)

  B(221) = 4.2

  B(222) = RCT(131)*V(96)

  B(223) = RCT(131)*V(73)

  B(224) = RCT(132)*V(97)

  B(225) = RCT(132)*V(73)

  B(226) = RCT(133)*V(87)

  B(227) = RCT(133)*V(63)

  B(228) = RCT(134)*V(90)

  B(229) = RCT(134)*V(63)

  B(230) = RCT(135)

  B(231) = RCT(136)*V(96)

  B(232) = RCT(136)*V(71)

  B(233) = RCT(137)*V(86)

  B(234) = RCT(137)*V(71)

  B(235) = RCT(138)*V(96)

  B(236) = RCT(138)*V(47)

  B(237) = RCT(139)*V(96)

  B(238) = RCT(139)*V(66)

  B(239) = RCT(140)

  B(240) = RCT(141)*V(93)

  B(241) = RCT(141)*V(80)

  B(242) = RCT(142)*V(96)

  B(243) = RCT(142)*V(80)

  B(244) = RCT(143)*V(86)

  B(245) = RCT(143)*V(80)

  B(246) = RCT(144)*V(97)

  B(247) = RCT(144)*V(80)

  B(248) = RCT(145)*V(96)

  B(249) = RCT(145)*V(81)

  B(250) = RCT(146)*V(86)

  B(251) = RCT(146)*V(81)

  B(252) = RCT(147)*V(97)

  B(253) = RCT(147)*V(81)

  B(254) = RCT(148)

  B(255) = RCT(149)*V(93)

  B(256) = RCT(149)*V(75)

  B(257) = RCT(150)*V(96)

  B(258) = RCT(150)*V(75)

  B(259) = RCT(151)*V(86)

  B(260) = RCT(151)*V(75)

  B(261) = RCT(152)*V(97)

  B(262) = RCT(152)*V(75)

  B(263) = RCT(153)*V(96)

  B(264) = RCT(153)*V(40)

  B(265) = RCT(154)*V(96)

  B(266) = RCT(154)*V(59)

  B(267) = RCT(155)*V(96)

  B(268) = RCT(155)*V(57)

  B(269) = RCT(156)*V(87)

  B(270) = RCT(156)*V(80)

  B(271) = RCT(157)

  B(272) = RCT(158)

  B(273) = RCT(159)*V(88)

  B(274) = RCT(159)*V(86)

  B(275) = RCT(160)*2*V(64)

  B(276) = RCT(161)*V(95)

  B(277) = RCT(161)*V(64)

  B(278) = RCT(162)*V(90)

  B(279) = RCT(162)*V(64)

  B(280) = RCT(163)*V(96)

  B(281) = RCT(163)*V(67)

  B(282) = RCT(164)

  B(283) = RCT(165)*V(88)

  B(284) = RCT(165)*V(48)

  B(285) = RCT(166)*V(88)

  B(286) = RCT(166)*V(79)

  B(287) = RCT(167)*V(88)

  B(288) = RCT(167)*V(57)

  B(289) = RCT(168)*V(88)

  B(290) = RCT(168)*V(74)

  B(291) = RCT(169)*V(88)

  B(292) = RCT(169)*V(77)

  B(293) = RCT(170)*V(88)

  B(294) = RCT(170)*V(76)

  B(295) = RCT(171)*V(88)

  B(296) = RCT(171)*V(80)

  B(297) = RCT(172)*V(88)

  B(298) = RCT(172)*V(85)

  B(299) = RCT(173)*V(88)

  B(300) = RCT(173)*V(83)

  B(301) = RCT(174)*V(88)

  B(302) = RCT(174)*V(82)

  B(303) = RCT(175)*V(88)

  B(304) = RCT(175)*V(58)

  B(305) = RCT(176)*V(88)

  B(306) = RCT(176)*V(59)

  B(307) = RCT(177)*V(96)

  B(308) = RCT(177)*V(61)

  B(309) = RCT(178)*V(86)

  B(310) = RCT(178)*V(49)

  B(311) = RCT(179)*V(96)

  B(312) = RCT(179)*V(49)

  B(313) = RCT(180)*V(60)

  B(314) = RCT(180)*V(49)

  B(315) = RCT(181)*V(96)

  B(316) = RCT(181)*V(12)

  B(317) = RCT(182)*V(96)

  B(318) = RCT(182)*V(15)

  B(319) = RCT(183)*V(96)

  B(320) = RCT(183)*V(18)

  B(321) = RCT(184)*V(96)

  B(322) = RCT(184)*V(23)

  B(323) = RCT(185)*V(86)

  B(324) = RCT(185)*V(23)

  B(325) = RCT(186)*V(96)

  B(326) = RCT(186)*V(29)

  B(327) = RCT(187)*V(86)

  B(328) = RCT(187)*V(29)

  B(329) = RCT(188)*V(97)

  B(330) = RCT(188)*V(29)

  B(331) = RCT(189)*V(96)

  B(332) = RCT(189)*V(32)

  B(333) = RCT(190)*V(96)

  B(334) = RCT(190)*V(34)

  B(335) = RCT(191)*V(96)

  B(336) = RCT(191)*V(37)



  JVS(1) = 0

  JVS(2) = 0.071*B(217)

  JVS(3) = 0.071*B(218)

  JVS(4) = 0

  JVS(5) = 0.138*B(217)

  JVS(6) = 0.138*B(218)

  JVS(7) = 0

  JVS(8) = B(222)+B(224)

  JVS(9) = B(223)

  JVS(10) = B(225)

  JVS(11) = 0

  JVS(12) = 0.038*B(235)

  JVS(13) = 0.038*B(236)

  JVS(14) = 0

  JVS(15) = 0.167*B(235)

  JVS(16) = 0.167*B(236)

  JVS(17) = 0

  JVS(18) = 0.232*B(240)+0.232*B(242)+0.232*B(244)+0.232*B(246)

  JVS(19) = 0.232*B(245)

  JVS(20) = 0.232*B(241)

  JVS(21) = 0.232*B(243)

  JVS(22) = 0.232*B(247)

  JVS(23) = 0

  JVS(24) = 0.0228*B(240)+0.0288*B(242)+0.0288*B(244)+0.0288*B(246)

  JVS(25) = 0.0288*B(245)

  JVS(26) = 0.0228*B(241)

  JVS(27) = 0.0288*B(243)

  JVS(28) = 0.0288*B(247)

  JVS(29) = 0

  JVS(30) = B(263)

  JVS(31) = B(264)

  JVS(32) = 0

  JVS(33) = B(263)

  JVS(34) = B(264)

  JVS(35) = 0

  JVS(36) = B(255)+B(257)+B(259)+B(261)

  JVS(37) = B(260)

  JVS(38) = B(256)

  JVS(39) = B(258)

  JVS(40) = B(262)

  JVS(41) = 0

  JVS(42) = B(315)

  JVS(43) = B(316)

  JVS(44) = -B(315)

  JVS(45) = -B(316)

  JVS(46) = 0

  JVS(47) = 0.239*B(317)

  JVS(48) = 0.239*B(318)

  JVS(49) = 0

  JVS(50) = 0.363*B(317)

  JVS(51) = 0.363*B(318)

  JVS(52) = -B(317)

  JVS(53) = -B(318)

  JVS(54) = 0

  JVS(55) = 0.045*B(319)

  JVS(56) = 0.045*B(320)

  JVS(57) = 0

  JVS(58) = 0.149*B(319)

  JVS(59) = 0.149*B(320)

  JVS(60) = -B(319)

  JVS(61) = -B(320)

  JVS(62) = 0

  JVS(63) = 0.038*B(321)

  JVS(64) = 0.038*B(322)

  JVS(65) = 0

  JVS(66) = 0.326*B(321)

  JVS(67) = 0.326*B(322)

  JVS(68) = 0

  JVS(69) = 0.125*B(323)

  JVS(70) = 0.125*B(324)

  JVS(71) = 0

  JVS(72) = 0.102*B(323)

  JVS(73) = 0.102*B(324)

  JVS(74) = -B(321)-B(323)

  JVS(75) = -B(324)

  JVS(76) = -B(322)

  JVS(77) = 0

  JVS(78) = 0.13*B(325)

  JVS(79) = 0.13*B(326)

  JVS(80) = 0

  JVS(81) = 0.0406*B(325)

  JVS(82) = 0.0406*B(326)

  JVS(83) = 0

  JVS(84) = 0.026*B(327)

  JVS(85) = 0.026*B(328)

  JVS(86) = 0

  JVS(87) = 0.485*B(327)

  JVS(88) = 0.485*B(328)

  JVS(89) = 0

  JVS(90) = B(329)

  JVS(91) = B(330)

  JVS(92) = -B(325)-B(327)-B(329)

  JVS(93) = -B(328)

  JVS(94) = -B(326)

  JVS(95) = -B(330)

  JVS(96) = 0

  JVS(97) = 0.091*B(331)

  JVS(98) = 0.091*B(332)

  JVS(99) = 0

  JVS(100) = 0.367*B(331)

  JVS(101) = 0.367*B(332)

  JVS(102) = -B(331)

  JVS(103) = -B(332)

  JVS(104) = 0

  JVS(105) = 1.173*B(333)

  JVS(106) = 1.173*B(334)

  JVS(107) = -B(333)

  JVS(108) = -B(334)

  JVS(109) = 0

  JVS(110) = 0.156*B(335)

  JVS(111) = 0.156*B(336)

  JVS(112) = 0

  JVS(113) = 0.777*B(335)

  JVS(114) = 0.777*B(336)

  JVS(115) = -B(335)

  JVS(116) = -B(336)

  JVS(117) = 0

  JVS(118) = B(309)+B(311)+B(313)

  JVS(119) = B(314)

  JVS(120) = B(310)

  JVS(121) = B(312)

  JVS(122) = -B(271)

  JVS(123) = 0.3*B(275)

  JVS(124) = -B(263)

  JVS(125) = -B(264)

  JVS(126) = -B(149)-B(150)

  JVS(127) = B(147)

  JVS(128) = B(148)

  JVS(129) = -B(272)

  JVS(130) = B(278)

  JVS(131) = B(279)

  JVS(132) = -B(217)

  JVS(133) = -B(218)

  JVS(134) = -B(63)-B(65)

  JVS(135) = -B(64)

  JVS(136) = -B(66)

  JVS(137) = -B(63)

  JVS(138) = -B(15)-B(17)-B(64)

  JVS(139) = B(14)

  JVS(140) = 0

  JVS(141) = -B(31)-B(33)-B(35)-B(88)

  JVS(142) = B(29)

  JVS(143) = B(30)

  JVS(144) = -B(235)

  JVS(145) = -B(236)

  JVS(146) = -B(109)-B(283)

  JVS(147) = -B(284)

  JVS(148) = -B(110)

  JVS(149) = -B(309)-B(311)-B(313)

  JVS(150) = -B(314)

  JVS(151) = -B(310)

  JVS(152) = -B(312)

  JVS(153) = -B(42)-B(43)-2*B(45)

  JVS(154) = 2*B(37)

  JVS(155) = 2*B(38)+B(40)

  JVS(156) = B(41)-B(44)

  JVS(157) = -B(136)

  JVS(158) = B(132)

  JVS(159) = 0.37*B(205)

  JVS(160) = 0.37*B(206)

  JVS(161) = B(133)

  JVS(162) = -B(137)

  JVS(163) = 0.56*B(217)

  JVS(164) = 0.3*B(235)

  JVS(165) = -B(219)-B(221)

  JVS(166) = -B(220)

  JVS(167) = 0.56*B(218)+0.3*B(236)

  JVS(168) = -B(54)-B(55)-B(86)

  JVS(169) = B(52)

  JVS(170) = B(53)

  JVS(171) = -B(56)

  JVS(172) = -B(158)-B(160)

  JVS(173) = 0.8*B(151)+0.8*B(178)

  JVS(174) = 0.8*B(152)

  JVS(175) = 0.8*B(179)

  JVS(176) = -B(159)

  JVS(177) = -B(161)

  JVS(178) = 0.1*B(155)+0.1*B(182)

  JVS(179) = 0.2*B(151)+0.2*B(178)

  JVS(180) = 0.2*B(152)+0.1*B(153)+0.1*B(156)

  JVS(181) = 0.1*B(154)+0.1*B(180)

  JVS(182) = 0.2*B(179)+0.1*B(181)+0.1*B(183)

  JVS(183) = -B(162)

  JVS(184) = -B(174)-B(175)-B(176)

  JVS(185) = B(172)

  JVS(186) = B(173)

  JVS(187) = -B(177)

  JVS(188) = -B(267)-B(287)

  JVS(189) = -B(288)

  JVS(190) = -B(268)

  JVS(191) = -B(119)-B(303)

  JVS(192) = -B(304)

  JVS(193) = 0.63*B(115)

  JVS(194) = -B(120)

  JVS(195) = -B(265)-B(305)

  JVS(196) = -B(306)

  JVS(197) = -B(266)

  JVS(198) = -B(313)

  JVS(199) = -B(60)-B(61)-B(75)-B(314)

  JVS(200) = 0

  JVS(201) = B(57)+B(58)

  JVS(202) = -B(76)

  JVS(203) = -B(62)+B(70)

  JVS(204) = B(283)

  JVS(205) = B(287)

  JVS(206) = B(303)

  JVS(207) = B(305)

  JVS(208) = -B(307)

  JVS(209) = 0.3*B(293)

  JVS(210) = B(285)

  JVS(211) = 0.15*B(295)

  JVS(212) = B(301)

  JVS(213) = B(299)

  JVS(214) = B(297)

  JVS(215) = B(284)+B(286)+B(288)+0.3*B(294)+0.15*B(296)+B(298)+B(300)+B(302)+B(304)+B(306)

  JVS(216) = 0

  JVS(217) = -B(308)

  JVS(218) = -B(131)-B(132)-B(134)

  JVS(219) = B(129)

  JVS(220) = B(130)-B(135)

  JVS(221) = -B(133)

  JVS(222) = -B(226)-B(228)

  JVS(223) = 0.4*B(222)+B(224)

  JVS(224) = -B(227)

  JVS(225) = -B(229)

  JVS(226) = 0.4*B(223)

  JVS(227) = B(225)

  JVS(228) = -2*B(275)-B(276)-B(278)

  JVS(229) = B(273)

  JVS(230) = B(274)

  JVS(231) = -B(279)

  JVS(232) = -B(277)

  JVS(233) = -B(104)-B(106)

  JVS(234) = B(95)

  JVS(235) = B(93)

  JVS(236) = B(94)+B(96)

  JVS(237) = -B(105)

  JVS(238) = 0.8*B(235)

  JVS(239) = -B(237)-B(239)

  JVS(240) = 0.2*B(233)

  JVS(241) = 0.168*B(248)+0.85*B(250)

  JVS(242) = 0.2*B(234)+0.85*B(251)

  JVS(243) = 0.8*B(236)-B(238)+0.168*B(249)

  JVS(244) = -B(280)-B(282)

  JVS(245) = B(289)

  JVS(246) = 0.7*B(293)

  JVS(247) = B(291)

  JVS(248) = 0.85*B(295)

  JVS(249) = B(290)+B(292)+0.7*B(294)+0.85*B(296)

  JVS(250) = -B(281)

  JVS(251) = B(134)

  JVS(252) = -B(116)-B(118)

  JVS(253) = 0

  JVS(254) = B(113)+B(135)

  JVS(255) = B(114)

  JVS(256) = 0

  JVS(257) = -B(117)

  JVS(258) = 2*B(31)+2*B(33)

  JVS(259) = -B(48)-B(87)

  JVS(260) = B(224)

  JVS(261) = 0.15*B(252)

  JVS(262) = B(167)

  JVS(263) = B(142)

  JVS(264) = B(101)

  JVS(265) = B(127)

  JVS(266) = B(46)

  JVS(267) = B(81)

  JVS(268) = B(47)-B(49)+B(102)

  JVS(269) = B(82)+B(128)+B(143)+B(168)+B(225)+0.15*B(253)

  JVS(270) = B(239)

  JVS(271) = B(280)+B(282)

  JVS(272) = -B(107)

  JVS(273) = B(230)+2*B(231)+0.69*B(233)

  JVS(274) = B(201)+0.63*B(205)

  JVS(275) = 0.001*B(259)

  JVS(276) = 0.1*B(209)+0.25*B(213)

  JVS(277) = 0.2*B(193)+0.33*B(197)

  JVS(278) = 0.066*B(244)

  JVS(279) = 0.334*B(248)+0.225*B(250)+0.643*B(252)+0.333*B(254)

  JVS(280) = B(169)

  JVS(281) = B(144)

  JVS(282) = B(121)+B(123)+B(124)+B(125)+B(127)+B(297)

  JVS(283) = 0.33*B(198)+0.63*B(206)+0.25*B(214)+0.69*B(234)+0.066*B(245)+0.225*B(251)+0.001*B(260)

  JVS(284) = B(298)

  JVS(285) = B(126)+0.2*B(194)+B(202)+0.1*B(210)

  JVS(286) = -B(108)+B(122)+2*B(232)+0.334*B(249)+B(281)

  JVS(287) = B(128)+0.643*B(253)

  JVS(288) = 0.9*B(219)

  JVS(289) = -B(230)-B(231)-B(233)

  JVS(290) = 0.3*B(222)

  JVS(291) = -B(234)

  JVS(292) = 0.9*B(220)

  JVS(293) = 0.3*B(223)-B(232)

  JVS(294) = -0.98*B(189)-B(190)-B(191)

  JVS(295) = 0.76*B(187)+0.76*B(285)

  JVS(296) = -B(192)

  JVS(297) = 0.76*B(286)

  JVS(298) = 0.76*B(188)

  JVS(299) = 0.36*B(217)

  JVS(300) = 0.2*B(235)

  JVS(301) = B(221)

  JVS(302) = B(228)

  JVS(303) = -B(222)-B(224)

  JVS(304) = 0

  JVS(305) = B(229)

  JVS(306) = 0

  JVS(307) = 0.36*B(218)-B(223)+0.2*B(236)

  JVS(308) = -B(225)

  JVS(309) = -B(201)-B(203)-B(205)-B(207)-B(289)

  JVS(310) = -B(206)

  JVS(311) = -B(290)

  JVS(312) = -B(202)

  JVS(313) = -B(204)

  JVS(314) = -B(208)

  JVS(315) = -B(255)-B(257)-B(259)-B(261)

  JVS(316) = -B(260)

  JVS(317) = -B(256)

  JVS(318) = -B(258)

  JVS(319) = -B(262)

  JVS(320) = -B(209)-B(211)-B(213)-B(215)-B(293)

  JVS(321) = -B(214)

  JVS(322) = -B(294)

  JVS(323) = -B(210)

  JVS(324) = -B(212)

  JVS(325) = -B(216)

  JVS(326) = 0.3*B(293)

  JVS(327) = -B(193)-B(195)-B(197)-B(199)-B(291)

  JVS(328) = -B(198)

  JVS(329) = -B(292)+0.3*B(294)

  JVS(330) = -B(194)

  JVS(331) = -B(196)

  JVS(332) = -B(200)

  JVS(333) = 0.009*B(267)+0.009*B(287)

  JVS(334) = 0.04*B(189)

  JVS(335) = 0.25*B(257)+0.18*B(259)+0.25*B(261)

  JVS(336) = 0.01*B(193)+0.09*B(199)

  JVS(337) = -B(91)-B(95)-2*B(98)-B(99)

  JVS(338) = 0.13*B(187)+0.13*B(285)

  JVS(339) = 0.088*B(242)

  JVS(340) = 0.18*B(260)

  JVS(341) = 0

  JVS(342) = 0.13*B(286)+0.009*B(288)

  JVS(343) = -B(100)

  JVS(344) = -B(96)

  JVS(345) = 0.01*B(194)

  JVS(346) = -B(92)

  JVS(347) = 0.13*B(188)+0.088*B(243)+0.25*B(258)+0.009*B(268)

  JVS(348) = 0.09*B(200)+0.25*B(262)

  JVS(349) = 1.1*B(235)

  JVS(350) = -2.1*B(189)

  JVS(351) = 5.12*B(255)+1.66*B(257)+7*B(259)

  JVS(352) = 0.1*B(209)+0.3*B(293)

  JVS(353) = 0.2*B(193)-0.7*B(195)-B(197)-B(199)-B(291)

  JVS(354) = -1.11*B(187)-1.11*B(285)

  JVS(355) = 0.35*B(244)+2.4*B(269)

  JVS(356) = 1.565*B(248)+0.36*B(250)+1.282*B(252)+0.832*B(254)

  JVS(357) = -0.66*B(101)-0.66*B(103)

  JVS(358) = -B(198)+0.35*B(245)+0.36*B(251)+7*B(260)

  JVS(359) = 2.4*B(270)

  JVS(360) = -1.11*B(286)-B(292)+0.3*B(294)

  JVS(361) = 0.2*B(194)+0.1*B(210)+5.12*B(256)

  JVS(362) = -0.66*B(102)-1.11*B(188)-0.7*B(196)+1.1*B(236)+1.565*B(249)+1.66*B(258)

  JVS(363) = -B(200)+1.282*B(253)

  JVS(364) = -B(240)-B(242)-B(244)-B(246)-B(269)-B(295)

  JVS(365) = -B(245)

  JVS(366) = -B(270)

  JVS(367) = -B(296)

  JVS(368) = -B(241)

  JVS(369) = -B(243)

  JVS(370) = -B(247)

  JVS(371) = 0.75*B(240)+0.912*B(242)+0.65*B(244)+0.2*B(246)+0.2*B(269)+B(295)

  JVS(372) = -B(248)-B(250)-B(252)-B(254)

  JVS(373) = 0.65*B(245)-B(251)

  JVS(374) = 0.2*B(270)

  JVS(375) = B(296)

  JVS(376) = 0.75*B(241)

  JVS(377) = 0.912*B(243)-B(249)

  JVS(378) = 0.2*B(247)-B(253)

  JVS(379) = 0.05*B(265)

  JVS(380) = 0.5*B(104)+0.5*B(106)

  JVS(381) = 0.03*B(233)

  JVS(382) = 0.5*B(189)

  JVS(383) = 0

  JVS(384) = 0.22*B(203)

  JVS(385) = 0.47*B(257)+0.21*B(259)+0.47*B(261)

  JVS(386) = 0.66*B(209)+0.7*B(211)+0.35*B(213)+0.64*B(215)+0.55*B(293)

  JVS(387) = 0.3*B(193)+0.62*B(195)+0.32*B(197)+0.56*B(199)+0.67*B(291)

  JVS(388) = 0

  JVS(389) = 0.05*B(187)+0.05*B(285)

  JVS(390) = 0.15*B(244)+0.8*B(246)+0.8*B(269)

  JVS(391) = 0.12*B(248)+0.357*B(252)

  JVS(392) = -B(163)-B(165)-B(167)-B(169)-B(301)

  JVS(393) = 0.33*B(101)+0.33*B(103)

  JVS(394) = 0.32*B(198)+0.35*B(214)+0.03*B(234)+0.15*B(245)+0.21*B(260)

  JVS(395) = 0.8*B(270)

  JVS(396) = 0.05*B(286)+0.67*B(292)+0.55*B(294)-B(302)

  JVS(397) = 0

  JVS(398) = 0

  JVS(399) = -B(164)+0.3*B(194)+0.66*B(210)

  JVS(400) = 0

  JVS(401) = 0.33*B(102)+0.5*B(105)-B(166)+0.05*B(188)+0.62*B(196)+0.22*B(204)+0.7*B(212)+0.12*B(249)+0.47*B(258)+0.05&
               &*B(266)

  JVS(402) = -B(168)+0.56*B(200)+0.64*B(216)+0.8*B(247)+0.357*B(253)+0.47*B(262)

  JVS(403) = B(176)

  JVS(404) = 0.991*B(267)+0.991*B(287)

  JVS(405) = 0.9*B(265)+B(305)

  JVS(406) = 0.5*B(104)+0.5*B(106)

  JVS(407) = 0.6*B(189)

  JVS(408) = 1.24*B(209)+1.3*B(211)+0.65*B(213)+1.18*B(215)+0.45*B(293)

  JVS(409) = 0.2*B(193)+0.33*B(195)+0.18*B(197)+0.35*B(199)+0.33*B(291)

  JVS(410) = 0

  JVS(411) = 0.06*B(187)+0.06*B(285)

  JVS(412) = 0

  JVS(413) = 0.252*B(248)+0.02*B(250)+0.067*B(254)

  JVS(414) = -B(138)-B(140)-B(142)-B(144)-B(299)

  JVS(415) = 0.33*B(101)+0.33*B(103)

  JVS(416) = 0.18*B(198)+0.65*B(214)+0.02*B(251)

  JVS(417) = 0

  JVS(418) = 0.06*B(286)+0.991*B(288)+0.33*B(292)+0.45*B(294)-B(300)+B(306)

  JVS(419) = 0.9*B(182)

  JVS(420) = 0

  JVS(421) = B(185)

  JVS(422) = 0.9*B(180)

  JVS(423) = -B(139)+0.2*B(194)+1.24*B(210)

  JVS(424) = B(170)+0.9*B(181)+0.9*B(183)+2*B(184)+B(186)

  JVS(425) = B(171)

  JVS(426) = 0.33*B(102)+0.5*B(105)-B(141)+B(177)+0.06*B(188)+0.33*B(196)+1.3*B(212)+0.252*B(249)+0.9*B(266)+0.991&
               &*B(268)

  JVS(427) = -B(143)+0.35*B(200)+1.18*B(216)

  JVS(428) = 0.1*B(219)

  JVS(429) = B(226)

  JVS(430) = B(191)

  JVS(431) = 0

  JVS(432) = 0.53*B(261)

  JVS(433) = B(91)

  JVS(434) = 0

  JVS(435) = 0.8*B(246)+0.8*B(269)

  JVS(436) = 0.85*B(252)

  JVS(437) = -B(101)-B(103)

  JVS(438) = 0

  JVS(439) = B(192)+B(227)+0.8*B(270)

  JVS(440) = 0

  JVS(441) = 0

  JVS(442) = 0

  JVS(443) = 0

  JVS(444) = B(92)+0.1*B(220)

  JVS(445) = -B(102)

  JVS(446) = 0.8*B(247)+0.85*B(253)+0.53*B(262)

  JVS(447) = B(119)+B(303)

  JVS(448) = 0.1*B(265)

  JVS(449) = B(131)

  JVS(450) = B(118)

  JVS(451) = B(231)+0.7*B(233)

  JVS(452) = 0

  JVS(453) = B(201)+1.56*B(203)+B(205)+2*B(207)+B(289)

  JVS(454) = 0.28*B(257)+0.24*B(259)

  JVS(455) = 0.25*B(213)

  JVS(456) = 0.2*B(193)+0.8*B(195)+0.74*B(197)+B(199)

  JVS(457) = 0.5*B(240)+0.629*B(242)+0.6*B(244)

  JVS(458) = 0.167*B(248)+0.15*B(250)+0.282*B(252)+0.9*B(254)

  JVS(459) = 0.33*B(101)+0.33*B(103)

  JVS(460) = -B(121)-B(123)-B(124)-B(125)-B(127)-B(129)-B(297)

  JVS(461) = 0.74*B(198)+B(206)+0.25*B(214)+0.7*B(234)+0.6*B(245)+0.15*B(251)+0.24*B(260)

  JVS(462) = 0

  JVS(463) = B(290)-B(298)+B(304)

  JVS(464) = 0

  JVS(465) = -B(130)

  JVS(466) = B(153)

  JVS(467) = B(111)+1.37*B(115)+B(154)+0.1*B(180)

  JVS(468) = -B(126)+0.2*B(194)+B(202)+0.5*B(241)

  JVS(469) = 0.1*B(181)

  JVS(470) = B(112)

  JVS(471) = 0.33*B(102)+B(120)-B(122)+0.8*B(196)+1.56*B(204)+B(232)+0.629*B(243)+0.167*B(249)+0.28*B(258)+0.1*B(266)

  JVS(472) = -B(128)+B(200)+2*B(208)+0.282*B(253)

  JVS(473) = -B(309)

  JVS(474) = 0

  JVS(475) = -B(233)

  JVS(476) = 0

  JVS(477) = -B(205)

  JVS(478) = -B(259)

  JVS(479) = -B(213)

  JVS(480) = -B(197)

  JVS(481) = -B(244)

  JVS(482) = -B(250)

  JVS(483) = -B(3)-B(11)-B(13)-B(14)-B(19)-B(21)-B(83)-B(198)-B(206)-B(214)-B(234)-B(245)-B(251)-B(260)-B(273)-B(310)

  JVS(484) = -B(12)

  JVS(485) = -B(274)

  JVS(486) = -B(22)+0.2*B(151)+0.2*B(178)

  JVS(487) = 0.2*B(152)

  JVS(488) = B(2)

  JVS(489) = 0.2*B(179)

  JVS(490) = -B(4)

  JVS(491) = -B(20)

  JVS(492) = -B(84)

  JVS(493) = B(149)+B(150)

  JVS(494) = B(35)+B(88)

  JVS(495) = B(43)+B(45)

  JVS(496) = 0.9*B(219)

  JVS(497) = B(54)+B(55)+0.61*B(86)

  JVS(498) = B(174)+B(175)+B(176)

  JVS(499) = B(132)

  JVS(500) = -B(226)

  JVS(501) = B(276)

  JVS(502) = B(87)

  JVS(503) = -B(191)

  JVS(504) = 0

  JVS(505) = B(207)

  JVS(506) = 0.47*B(261)

  JVS(507) = B(215)

  JVS(508) = B(199)

  JVS(509) = 0

  JVS(510) = 0.2*B(246)-B(269)

  JVS(511) = 0

  JVS(512) = 0

  JVS(513) = 0

  JVS(514) = B(103)

  JVS(515) = 0

  JVS(516) = B(3)-B(11)+B(83)

  JVS(517) = -B(1)-B(5)-B(7)-B(12)-B(29)-B(37)-B(46)-B(52)-B(147)-B(172)-B(192)-B(227)-B(270)

  JVS(518) = 0

  JVS(519) = B(89)

  JVS(520) = B(50)-B(53)

  JVS(521) = B(145)-B(148)

  JVS(522) = B(111)

  JVS(523) = -B(6)-B(8)+B(9)+B(77)

  JVS(524) = B(170)-B(173)

  JVS(525) = B(4)+B(10)+2*B(25)+2*B(36)-B(38)+B(51)+B(90)+B(112)+B(133)+B(146)+B(171)+0.9*B(220)+B(277)

  JVS(526) = B(44)-B(47)+B(56)+B(79)+B(177)

  JVS(527) = B(23)+2*B(26)-B(30)+B(78)+B(80)+B(84)+2*B(85)+B(200)+B(208)+B(216)+0.2*B(247)+0.47*B(262)

  JVS(528) = 2*B(271)

  JVS(529) = B(272)

  JVS(530) = -B(283)

  JVS(531) = -B(287)

  JVS(532) = -B(303)

  JVS(533) = -B(305)

  JVS(534) = B(307)

  JVS(535) = 1.4*B(275)+B(276)

  JVS(536) = B(280)+B(282)

  JVS(537) = -B(289)

  JVS(538) = -B(293)

  JVS(539) = -B(291)

  JVS(540) = -B(285)

  JVS(541) = -B(295)

  JVS(542) = 0

  JVS(543) = -B(301)

  JVS(544) = -B(299)

  JVS(545) = 0

  JVS(546) = -B(297)

  JVS(547) = -B(273)

  JVS(548) = 0

  JVS(549) = -B(274)-B(284)-B(286)-B(288)-B(290)-B(292)-B(294)-B(296)-B(298)-B(300)-B(302)-B(304)-B(306)

  JVS(550) = 0

  JVS(551) = 0

  JVS(552) = 0

  JVS(553) = 0

  JVS(554) = 0

  JVS(555) = 0

  JVS(556) = B(277)

  JVS(557) = B(281)+B(308)

  JVS(558) = 0

  JVS(559) = 0.08*B(217)

  JVS(560) = 0.5*B(235)

  JVS(561) = 0.991*B(267)+0.991*B(287)

  JVS(562) = 0.1*B(265)

  JVS(563) = B(104)

  JVS(564) = B(237)

  JVS(565) = 0.3*B(116)

  JVS(566) = B(231)+0.03*B(233)

  JVS(567) = 0.96*B(189)

  JVS(568) = 0.6*B(222)

  JVS(569) = 0.7*B(201)+B(203)+B(207)+2*B(289)

  JVS(570) = 1.25*B(257)+0.76*B(259)+1.03*B(261)

  JVS(571) = 0.1*B(209)+B(211)+1.7*B(293)

  JVS(572) = 0.2*B(193)+0.8*B(195)+0.22*B(197)+0.91*B(199)+2*B(291)

  JVS(573) = -B(99)

  JVS(574) = 0.87*B(187)+0.87*B(285)

  JVS(575) = 0.25*B(240)+0.991*B(242)+0.2*B(244)+B(246)+B(269)+B(295)

  JVS(576) = 0.713*B(248)+0.064*B(250)+0.075*B(252)+0.7*B(254)

  JVS(577) = 0

  JVS(578) = 0

  JVS(579) = 0.22*B(198)+0.03*B(234)+0.2*B(245)+0.064*B(251)+0.76*B(260)

  JVS(580) = B(270)

  JVS(581) = 0.87*B(286)+0.991*B(288)+2*B(290)+2*B(292)+1.7*B(294)+B(296)

  JVS(582) = -B(89)-B(93)-2*B(97)-B(100)-B(155)-B(182)

  JVS(583) = -B(94)

  JVS(584) = -B(156)+B(185)

  JVS(585) = 0.9*B(180)

  JVS(586) = 0.2*B(194)+0.7*B(202)+0.1*B(210)+0.25*B(241)

  JVS(587) = B(170)+0.9*B(181)-B(183)+2*B(184)+B(186)

  JVS(588) = -B(90)+B(171)

  JVS(589) = B(105)+0.3*B(117)+0.87*B(188)+0.8*B(196)+B(204)+B(212)+0.08*B(218)+0.6*B(223)+B(232)+0.5*B(236)+B(238)&
               &+0.991*B(243)+0.713*B(249)+1.25*B(258)+0.1*B(266)+0.991*B(268)

  JVS(590) = 0.91*B(200)+B(208)+B(247)+0.075*B(253)+1.03*B(262)

  JVS(591) = B(263)

  JVS(592) = 0.44*B(217)

  JVS(593) = B(63)+B(65)

  JVS(594) = B(64)

  JVS(595) = 0.7*B(235)

  JVS(596) = B(136)

  JVS(597) = 0.9*B(219)+B(221)

  JVS(598) = B(54)+0.61*B(86)

  JVS(599) = B(267)+B(287)

  JVS(600) = B(119)+B(303)

  JVS(601) = B(265)+B(305)

  JVS(602) = B(61)+B(75)

  JVS(603) = B(131)+B(132)-B(134)

  JVS(604) = -B(228)

  JVS(605) = -B(278)

  JVS(606) = B(106)

  JVS(607) = B(239)

  JVS(608) = B(282)

  JVS(609) = 0.3*B(116)+B(118)

  JVS(610) = B(107)

  JVS(611) = B(230)+2*B(231)+0.76*B(233)

  JVS(612) = 0.94*B(189)+B(190)

  JVS(613) = 0.6*B(222)

  JVS(614) = 1.7*B(201)+B(203)+0.13*B(205)+B(289)

  JVS(615) = 0.75*B(257)+0.07*B(259)+0.28*B(261)

  JVS(616) = 0.1*B(209)+B(211)+0.5*B(213)+B(215)+B(293)

  JVS(617) = 0.3*B(193)+0.95*B(195)+0.44*B(197)+B(291)

  JVS(618) = -B(95)

  JVS(619) = 0.11*B(187)+0.11*B(285)

  JVS(620) = 0.25*B(240)+0.912*B(242)+0.066*B(244)+0.8*B(246)+0.8*B(269)+B(295)

  JVS(621) = 0.503*B(248)+0.154*B(250)+0.925*B(252)+1.033*B(254)

  JVS(622) = B(169)

  JVS(623) = B(144)

  JVS(624) = B(101)+B(103)

  JVS(625) = B(121)+2*B(123)+B(125)+B(127)-B(129)+B(297)

  JVS(626) = B(19)-B(21)+0.44*B(198)+0.13*B(206)+0.5*B(214)+0.76*B(234)+0.066*B(245)+0.154*B(251)+0.07*B(260)

  JVS(627) = -B(52)+0.8*B(270)

  JVS(628) = 0.11*B(286)+B(288)+B(290)+B(292)+B(294)+B(296)+B(298)+B(304)+B(306)

  JVS(629) = -B(93)

  JVS(630) = -B(22)-B(50)-B(53)-2*B(57)-2*B(58)-B(71)-B(73)-B(81)-B(94)-B(96)-B(113)-B(130)-B(135)-B(151)-B(178)-B(229)&
               &-B(279)

  JVS(631) = -B(152)+0.9*B(153)+B(185)

  JVS(632) = B(111)-B(114)+0.74*B(115)+0.9*B(154)+B(180)

  JVS(633) = B(67)-B(74)+B(76)+B(126)+0.3*B(194)+1.7*B(202)+0.1*B(210)+0.25*B(241)

  JVS(634) = B(170)-B(179)+B(181)+2*B(184)+B(186)

  JVS(635) = -B(51)+B(112)+B(133)+B(171)+0.9*B(220)

  JVS(636) = B(20)+B(62)+B(66)+B(68)-B(72)+B(79)+B(102)+B(108)+0.3*B(117)+B(120)+B(122)+B(137)+0.11*B(188)+0.95*B(196)&
               &+B(204)+B(212)+0.44*B(218)+0.6*B(223)+2*B(232)+0.7*B(236)+0.912*B(243)+0.503*B(249)+0.75*B(258)+B(264)&
               &+B(266)+B(268)

  JVS(637) = B(80)-B(82)+B(128)+B(216)+0.8*B(247)+0.925*B(253)+0.28*B(262)

  JVS(638) = B(149)+B(150)

  JVS(639) = B(158)

  JVS(640) = B(237)+B(239)

  JVS(641) = B(230)+B(231)+0.62*B(233)

  JVS(642) = 0

  JVS(643) = 0.21*B(248)+0.114*B(250)+0.967*B(254)

  JVS(644) = B(138)+B(140)+B(142)+B(299)

  JVS(645) = 0

  JVS(646) = 0.62*B(234)+0.114*B(251)

  JVS(647) = -B(147)

  JVS(648) = B(300)

  JVS(649) = -B(155)

  JVS(650) = -B(151)

  JVS(651) = -B(145)-B(148)-B(152)-B(153)-B(156)-2*B(157)-B(185)

  JVS(652) = -B(154)

  JVS(653) = B(139)

  JVS(654) = -B(186)

  JVS(655) = -B(146)

  JVS(656) = B(141)+B(159)+B(232)+B(238)+0.21*B(249)

  JVS(657) = B(143)

  JVS(658) = B(109)+B(283)

  JVS(659) = B(160)

  JVS(660) = B(161)

  JVS(661) = 0.7*B(116)

  JVS(662) = B(169)

  JVS(663) = B(144)

  JVS(664) = 0

  JVS(665) = 0

  JVS(666) = 0

  JVS(667) = 0

  JVS(668) = B(284)

  JVS(669) = 0.9*B(155)

  JVS(670) = -B(113)

  JVS(671) = B(145)-0.1*B(153)+0.9*B(156)+2*B(157)+B(185)

  JVS(672) = -B(111)-B(114)-2*B(115)-0.1*B(154)-B(180)

  JVS(673) = 0

  JVS(674) = -B(181)+B(186)

  JVS(675) = -B(112)+B(146)

  JVS(676) = B(110)+0.7*B(117)+B(162)

  JVS(677) = 0

  JVS(678) = B(15)

  JVS(679) = -B(75)

  JVS(680) = -B(201)

  JVS(681) = -B(255)

  JVS(682) = -B(209)+0.5*B(213)

  JVS(683) = -B(193)

  JVS(684) = -B(240)

  JVS(685) = -B(163)

  JVS(686) = -B(138)

  JVS(687) = 0

  JVS(688) = -B(125)

  JVS(689) = B(13)+0.5*B(214)

  JVS(690) = B(1)-B(5)-B(7)

  JVS(691) = 0

  JVS(692) = 0

  JVS(693) = -B(73)

  JVS(694) = 0

  JVS(695) = 0

  JVS(696) = -B(2)-B(6)-B(8)-B(9)-B(67)-B(74)-B(76)-B(77)-B(126)-B(139)-B(164)-B(194)-B(202)-B(210)-B(241)-B(256)

  JVS(697) = 0

  JVS(698) = -B(10)

  JVS(699) = -B(68)+B(69)

  JVS(700) = B(23)-B(78)

  JVS(701) = B(174)+B(175)

  JVS(702) = 0.39*B(259)

  JVS(703) = 0.25*B(240)+0.2*B(244)

  JVS(704) = 0.25*B(248)+0.075*B(252)

  JVS(705) = B(163)+B(165)+B(167)+B(301)

  JVS(706) = 0

  JVS(707) = 0.2*B(245)+0.39*B(260)

  JVS(708) = -B(172)

  JVS(709) = B(302)

  JVS(710) = -B(182)

  JVS(711) = -B(178)

  JVS(712) = -B(185)

  JVS(713) = -B(180)

  JVS(714) = B(164)+0.25*B(241)

  JVS(715) = -B(170)-B(173)-B(179)-B(181)-B(183)-2*B(184)-B(186)

  JVS(716) = -B(171)

  JVS(717) = B(166)+0.25*B(249)

  JVS(718) = B(168)+0.075*B(253)

  JVS(719) = B(42)+B(45)

  JVS(720) = -B(219)

  JVS(721) = -B(132)

  JVS(722) = -B(276)

  JVS(723) = -B(91)

  JVS(724) = 0

  JVS(725) = 0.2*B(269)

  JVS(726) = 0

  JVS(727) = 0

  JVS(728) = 0

  JVS(729) = -B(3)

  JVS(730) = B(1)+B(5)+B(27)-B(37)+0.2*B(270)

  JVS(731) = 0

  JVS(732) = -B(89)

  JVS(733) = -B(50)

  JVS(734) = -B(145)

  JVS(735) = -B(111)

  JVS(736) = B(6)-B(9)

  JVS(737) = -B(170)

  JVS(738) = -B(4)-B(10)-B(25)-2*B(36)-B(38)-B(40)-B(51)-B(90)-B(92)-B(112)-B(133)-B(146)-B(171)-B(220)-B(277)

  JVS(739) = -B(41)

  JVS(740) = B(24)-B(26)+B(28)

  JVS(741) = -B(263)

  JVS(742) = B(272)

  JVS(743) = -B(217)

  JVS(744) = B(63)-B(65)

  JVS(745) = 2*B(17)+B(64)

  JVS(746) = -B(235)

  JVS(747) = -B(109)

  JVS(748) = -B(311)

  JVS(749) = B(42)-B(43)

  JVS(750) = -B(136)

  JVS(751) = -B(55)+0.39*B(86)

  JVS(752) = -B(158)+B(160)

  JVS(753) = -B(161)

  JVS(754) = -B(176)

  JVS(755) = -B(267)

  JVS(756) = -B(119)

  JVS(757) = -B(265)

  JVS(758) = 2*B(60)-B(61)+B(75)

  JVS(759) = -B(307)

  JVS(760) = 0

  JVS(761) = 0

  JVS(762) = -B(104)+B(106)

  JVS(763) = -B(237)

  JVS(764) = -B(280)

  JVS(765) = -B(116)+B(118)

  JVS(766) = -B(48)+B(87)

  JVS(767) = -B(107)

  JVS(768) = -B(231)+0.08*B(233)

  JVS(769) = -B(222)

  JVS(770) = 0.3*B(201)-B(203)+0.13*B(205)

  JVS(771) = 0.15*B(255)-B(257)+0.57*B(259)

  JVS(772) = -B(211)+0.5*B(213)

  JVS(773) = 0.1*B(193)-B(195)+0.1*B(197)

  JVS(774) = 0

  JVS(775) = -B(187)

  JVS(776) = -B(242)+0.266*B(244)

  JVS(777) = -B(248)+0.268*B(250)

  JVS(778) = B(163)-B(165)

  JVS(779) = B(138)-B(140)

  JVS(780) = -B(101)

  JVS(781) = -B(121)+B(125)

  JVS(782) = -B(19)+B(21)+0.1*B(198)+0.13*B(206)+0.5*B(214)+0.08*B(234)+0.266*B(245)+0.268*B(251)+0.57*B(260)

  JVS(783) = -B(46)

  JVS(784) = 0

  JVS(785) = 0

  JVS(786) = B(22)+B(50)-B(71)+B(73)

  JVS(787) = 0

  JVS(788) = 0

  JVS(789) = -B(67)+B(74)+B(76)+B(126)+B(139)+B(164)+0.1*B(194)+0.3*B(202)+0.15*B(256)

  JVS(790) = 0

  JVS(791) = -B(40)+B(51)

  JVS(792) = -B(20)-B(41)-B(44)-B(47)-B(49)-B(56)-B(62)-B(66)-B(68)-2*B(69)-2*B(70)-B(72)-B(79)-B(102)-B(105)-B(108)&
               &-B(110)-B(117)-B(120)-B(122)-B(137)-B(141)-B(159)-B(162)-B(166)-B(177)-B(188)-B(196)-B(204)-B(212)-B(218)&
               &-B(223)-B(232)-B(236)-B(238)-B(243)-B(249)-B(258)-B(264)-B(266)-B(268)-B(281)-B(308)-B(312)

  JVS(793) = -B(80)

  JVS(794) = B(35)+B(88)

  JVS(795) = 0.39*B(86)

  JVS(796) = B(48)

  JVS(797) = -B(224)

  JVS(798) = -B(207)

  JVS(799) = -B(261)

  JVS(800) = -B(215)

  JVS(801) = -B(199)

  JVS(802) = -B(246)

  JVS(803) = -B(252)

  JVS(804) = -B(167)

  JVS(805) = -B(142)

  JVS(806) = 0

  JVS(807) = -B(127)

  JVS(808) = B(11)-B(83)

  JVS(809) = B(7)+B(12)-B(27)-B(29)

  JVS(810) = 0

  JVS(811) = 0

  JVS(812) = -B(81)

  JVS(813) = 0

  JVS(814) = 0

  JVS(815) = B(8)-B(77)

  JVS(816) = 0

  JVS(817) = -B(25)

  JVS(818) = B(49)-B(79)

  JVS(819) = -B(23)-B(24)-B(26)-B(28)-B(30)-B(78)-B(80)-B(82)-B(84)-2*B(85)-B(128)-B(143)-B(168)-B(200)-B(208)-B(216)&
               &-B(225)-B(247)-B(253)-B(262)
      
END SUBROUTINE cb05_sorg_aq_Jac_SP














SUBROUTINE cb05_sorg_aq_KppDecomp( JVS, IER )







      INTEGER  :: IER
      REAL(kind=dp) :: JVS(819), W(97), a
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
      
END SUBROUTINE cb05_sorg_aq_KppDecomp



SUBROUTINE cb05_sorg_aq_KppDecompCmplx( JVS, IER )







      INTEGER  :: IER
      DOUBLE COMPLEX :: JVS(819), W(97), a
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
      
END SUBROUTINE cb05_sorg_aq_KppDecompCmplx


SUBROUTINE cb05_sorg_aq_KppSolveIndirect( JVS, X )







      INTEGER i, j
      REAL(kind=dp) JVS(819), X(97), sum

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
      
END SUBROUTINE cb05_sorg_aq_KppSolveIndirect


SUBROUTINE cb05_sorg_aq_KppSolveCmplx( JVS, X )







      INTEGER i, j
      DOUBLE COMPLEX JVS(819), X(97), sum

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
      
END SUBROUTINE cb05_sorg_aq_KppSolveCmplx













SUBROUTINE cb05_sorg_aq_KppSolve ( JVS, X )


  REAL(kind=dp) :: JVS(LU_NONZERO)

  REAL(kind=dp) :: X(NVAR)

  X(45) = X(45)-JVS(137)*X(44)
  X(52) = X(52)-JVS(163)*X(43)-JVS(164)*X(47)
  X(60) = X(60)-JVS(198)*X(49)
  X(61) = X(61)-JVS(204)*X(48)-JVS(205)*X(57)-JVS(206)*X(58)-JVS(207)*X(59)
  X(66) = X(66)-JVS(238)*X(47)
  X(68) = X(68)-JVS(251)*X(62)
  X(69) = X(69)-JVS(258)*X(46)
  X(70) = X(70)-JVS(270)*X(66)-JVS(271)*X(67)
  X(71) = X(71)-JVS(288)*X(52)
  X(73) = X(73)-JVS(299)*X(43)-JVS(300)*X(47)-JVS(301)*X(52)-JVS(302)*X(63)
  X(77) = X(77)-JVS(326)*X(76)
  X(78) = X(78)-JVS(333)*X(57)-JVS(334)*X(72)-JVS(335)*X(75)-JVS(336)*X(77)
  X(79) = X(79)-JVS(349)*X(47)-JVS(350)*X(72)-JVS(351)*X(75)-JVS(352)*X(76)-JVS(353)*X(77)
  X(81) = X(81)-JVS(371)*X(80)
  X(82) = X(82)-JVS(379)*X(59)-JVS(380)*X(65)-JVS(381)*X(71)-JVS(382)*X(72)-JVS(383)*X(73)-JVS(384)*X(74)-JVS(385)*X(75)&
            &-JVS(386)*X(76)-JVS(387)*X(77)-JVS(388)*X(78)-JVS(389)*X(79)-JVS(390)*X(80)-JVS(391)*X(81)
  X(83) = X(83)-JVS(403)*X(56)-JVS(404)*X(57)-JVS(405)*X(59)-JVS(406)*X(65)-JVS(407)*X(72)-JVS(408)*X(76)-JVS(409)*X(77)&
            &-JVS(410)*X(78)-JVS(411)*X(79)-JVS(412)*X(80)-JVS(413)*X(81)
  X(84) = X(84)-JVS(428)*X(52)-JVS(429)*X(63)-JVS(430)*X(72)-JVS(431)*X(73)-JVS(432)*X(75)-JVS(433)*X(78)-JVS(434)*X(79)&
            &-JVS(435)*X(80)-JVS(436)*X(81)
  X(85) = X(85)-JVS(447)*X(58)-JVS(448)*X(59)-JVS(449)*X(62)-JVS(450)*X(68)-JVS(451)*X(71)-JVS(452)*X(73)-JVS(453)*X(74)&
            &-JVS(454)*X(75)-JVS(455)*X(76)-JVS(456)*X(77)-JVS(457)*X(80)-JVS(458)*X(81)-JVS(459)*X(84)
  X(86) = X(86)-JVS(473)*X(49)-JVS(474)*X(60)-JVS(475)*X(71)-JVS(476)*X(73)-JVS(477)*X(74)-JVS(478)*X(75)-JVS(479)*X(76)&
            &-JVS(480)*X(77)-JVS(481)*X(80)-JVS(482)*X(81)
  X(87) = X(87)-JVS(493)*X(41)-JVS(494)*X(46)-JVS(495)*X(50)-JVS(496)*X(52)-JVS(497)*X(53)-JVS(498)*X(56)-JVS(499)*X(62)&
            &-JVS(500)*X(63)-JVS(501)*X(64)-JVS(502)*X(69)-JVS(503)*X(72)-JVS(504)*X(73)-JVS(505)*X(74)-JVS(506)*X(75)&
            &-JVS(507)*X(76)-JVS(508)*X(77)-JVS(509)*X(79)-JVS(510)*X(80)-JVS(511)*X(81)-JVS(512)*X(82)-JVS(513)*X(83)&
            &-JVS(514)*X(84)-JVS(515)*X(85)-JVS(516)*X(86)
  X(88) = X(88)-JVS(528)*X(39)-JVS(529)*X(42)-JVS(530)*X(48)-JVS(531)*X(57)-JVS(532)*X(58)-JVS(533)*X(59)-JVS(534)*X(61)&
            &-JVS(535)*X(64)-JVS(536)*X(67)-JVS(537)*X(74)-JVS(538)*X(76)-JVS(539)*X(77)-JVS(540)*X(79)-JVS(541)*X(80)&
            &-JVS(542)*X(81)-JVS(543)*X(82)-JVS(544)*X(83)-JVS(545)*X(84)-JVS(546)*X(85)-JVS(547)*X(86)-JVS(548)*X(87)
  X(89) = X(89)-JVS(559)*X(43)-JVS(560)*X(47)-JVS(561)*X(57)-JVS(562)*X(59)-JVS(563)*X(65)-JVS(564)*X(66)-JVS(565)*X(68)&
            &-JVS(566)*X(71)-JVS(567)*X(72)-JVS(568)*X(73)-JVS(569)*X(74)-JVS(570)*X(75)-JVS(571)*X(76)-JVS(572)*X(77)&
            &-JVS(573)*X(78)-JVS(574)*X(79)-JVS(575)*X(80)-JVS(576)*X(81)-JVS(577)*X(84)-JVS(578)*X(85)-JVS(579)*X(86)&
            &-JVS(580)*X(87)-JVS(581)*X(88)
  X(90) = X(90)-JVS(591)*X(40)-JVS(592)*X(43)-JVS(593)*X(44)-JVS(594)*X(45)-JVS(595)*X(47)-JVS(596)*X(51)-JVS(597)*X(52)&
            &-JVS(598)*X(53)-JVS(599)*X(57)-JVS(600)*X(58)-JVS(601)*X(59)-JVS(602)*X(60)-JVS(603)*X(62)-JVS(604)*X(63)&
            &-JVS(605)*X(64)-JVS(606)*X(65)-JVS(607)*X(66)-JVS(608)*X(67)-JVS(609)*X(68)-JVS(610)*X(70)-JVS(611)*X(71)&
            &-JVS(612)*X(72)-JVS(613)*X(73)-JVS(614)*X(74)-JVS(615)*X(75)-JVS(616)*X(76)-JVS(617)*X(77)-JVS(618)*X(78)&
            &-JVS(619)*X(79)-JVS(620)*X(80)-JVS(621)*X(81)-JVS(622)*X(82)-JVS(623)*X(83)-JVS(624)*X(84)-JVS(625)*X(85)&
            &-JVS(626)*X(86)-JVS(627)*X(87)-JVS(628)*X(88)-JVS(629)*X(89)
  X(91) = X(91)-JVS(638)*X(41)-JVS(639)*X(54)-JVS(640)*X(66)-JVS(641)*X(71)-JVS(642)*X(73)-JVS(643)*X(81)-JVS(644)*X(83)&
            &-JVS(645)*X(84)-JVS(646)*X(86)-JVS(647)*X(87)-JVS(648)*X(88)-JVS(649)*X(89)-JVS(650)*X(90)
  X(92) = X(92)-JVS(658)*X(48)-JVS(659)*X(54)-JVS(660)*X(55)-JVS(661)*X(68)-JVS(662)*X(82)-JVS(663)*X(83)-JVS(664)*X(84)&
            &-JVS(665)*X(85)-JVS(666)*X(86)-JVS(667)*X(87)-JVS(668)*X(88)-JVS(669)*X(89)-JVS(670)*X(90)-JVS(671)*X(91)
  X(93) = X(93)-JVS(678)*X(45)-JVS(679)*X(60)-JVS(680)*X(74)-JVS(681)*X(75)-JVS(682)*X(76)-JVS(683)*X(77)-JVS(684)*X(80)&
            &-JVS(685)*X(82)-JVS(686)*X(83)-JVS(687)*X(84)-JVS(688)*X(85)-JVS(689)*X(86)-JVS(690)*X(87)-JVS(691)*X(88)&
            &-JVS(692)*X(89)-JVS(693)*X(90)-JVS(694)*X(91)-JVS(695)*X(92)
  X(94) = X(94)-JVS(701)*X(56)-JVS(702)*X(75)-JVS(703)*X(80)-JVS(704)*X(81)-JVS(705)*X(82)-JVS(706)*X(84)-JVS(707)*X(86)&
            &-JVS(708)*X(87)-JVS(709)*X(88)-JVS(710)*X(89)-JVS(711)*X(90)-JVS(712)*X(91)-JVS(713)*X(92)-JVS(714)*X(93)
  X(95) = X(95)-JVS(719)*X(50)-JVS(720)*X(52)-JVS(721)*X(62)-JVS(722)*X(64)-JVS(723)*X(78)-JVS(724)*X(79)-JVS(725)*X(80)&
            &-JVS(726)*X(81)-JVS(727)*X(84)-JVS(728)*X(85)-JVS(729)*X(86)-JVS(730)*X(87)-JVS(731)*X(88)-JVS(732)*X(89)&
            &-JVS(733)*X(90)-JVS(734)*X(91)-JVS(735)*X(92)-JVS(736)*X(93)-JVS(737)*X(94)
  X(96) = X(96)-JVS(741)*X(40)-JVS(742)*X(42)-JVS(743)*X(43)-JVS(744)*X(44)-JVS(745)*X(45)-JVS(746)*X(47)-JVS(747)*X(48)&
            &-JVS(748)*X(49)-JVS(749)*X(50)-JVS(750)*X(51)-JVS(751)*X(53)-JVS(752)*X(54)-JVS(753)*X(55)-JVS(754)*X(56)&
            &-JVS(755)*X(57)-JVS(756)*X(58)-JVS(757)*X(59)-JVS(758)*X(60)-JVS(759)*X(61)-JVS(760)*X(62)-JVS(761)*X(64)&
            &-JVS(762)*X(65)-JVS(763)*X(66)-JVS(764)*X(67)-JVS(765)*X(68)-JVS(766)*X(69)-JVS(767)*X(70)-JVS(768)*X(71)&
            &-JVS(769)*X(73)-JVS(770)*X(74)-JVS(771)*X(75)-JVS(772)*X(76)-JVS(773)*X(77)-JVS(774)*X(78)-JVS(775)*X(79)&
            &-JVS(776)*X(80)-JVS(777)*X(81)-JVS(778)*X(82)-JVS(779)*X(83)-JVS(780)*X(84)-JVS(781)*X(85)-JVS(782)*X(86)&
            &-JVS(783)*X(87)-JVS(784)*X(88)-JVS(785)*X(89)-JVS(786)*X(90)-JVS(787)*X(91)-JVS(788)*X(92)-JVS(789)*X(93)&
            &-JVS(790)*X(94)-JVS(791)*X(95)
  X(97) = X(97)-JVS(794)*X(46)-JVS(795)*X(53)-JVS(796)*X(69)-JVS(797)*X(73)-JVS(798)*X(74)-JVS(799)*X(75)-JVS(800)*X(76)&
            &-JVS(801)*X(77)-JVS(802)*X(80)-JVS(803)*X(81)-JVS(804)*X(82)-JVS(805)*X(83)-JVS(806)*X(84)-JVS(807)*X(85)&
            &-JVS(808)*X(86)-JVS(809)*X(87)-JVS(810)*X(88)-JVS(811)*X(89)-JVS(812)*X(90)-JVS(813)*X(91)-JVS(814)*X(92)&
            &-JVS(815)*X(93)-JVS(816)*X(94)-JVS(817)*X(95)-JVS(818)*X(96)
  X(97) = X(97)/JVS(819)
  X(96) = (X(96)-JVS(793)*X(97))/(JVS(792))
  X(95) = (X(95)-JVS(739)*X(96)-JVS(740)*X(97))/(JVS(738))
  X(94) = (X(94)-JVS(716)*X(95)-JVS(717)*X(96)-JVS(718)*X(97))/(JVS(715))
  X(93) = (X(93)-JVS(697)*X(94)-JVS(698)*X(95)-JVS(699)*X(96)-JVS(700)*X(97))/(JVS(696))
  X(92) = (X(92)-JVS(673)*X(93)-JVS(674)*X(94)-JVS(675)*X(95)-JVS(676)*X(96)-JVS(677)*X(97))/(JVS(672))
  X(91) = (X(91)-JVS(652)*X(92)-JVS(653)*X(93)-JVS(654)*X(94)-JVS(655)*X(95)-JVS(656)*X(96)-JVS(657)*X(97))/(JVS(651))
  X(90) = (X(90)-JVS(631)*X(91)-JVS(632)*X(92)-JVS(633)*X(93)-JVS(634)*X(94)-JVS(635)*X(95)-JVS(636)*X(96)-JVS(637)&
            &*X(97))/(JVS(630))
  X(89) = (X(89)-JVS(583)*X(90)-JVS(584)*X(91)-JVS(585)*X(92)-JVS(586)*X(93)-JVS(587)*X(94)-JVS(588)*X(95)-JVS(589)&
            &*X(96)-JVS(590)*X(97))/(JVS(582))
  X(88) = (X(88)-JVS(550)*X(89)-JVS(551)*X(90)-JVS(552)*X(91)-JVS(553)*X(92)-JVS(554)*X(93)-JVS(555)*X(94)-JVS(556)&
            &*X(95)-JVS(557)*X(96)-JVS(558)*X(97))/(JVS(549))
  X(87) = (X(87)-JVS(518)*X(88)-JVS(519)*X(89)-JVS(520)*X(90)-JVS(521)*X(91)-JVS(522)*X(92)-JVS(523)*X(93)-JVS(524)&
            &*X(94)-JVS(525)*X(95)-JVS(526)*X(96)-JVS(527)*X(97))/(JVS(517))
  X(86) = (X(86)-JVS(484)*X(87)-JVS(485)*X(88)-JVS(486)*X(90)-JVS(487)*X(91)-JVS(488)*X(93)-JVS(489)*X(94)-JVS(490)&
            &*X(95)-JVS(491)*X(96)-JVS(492)*X(97))/(JVS(483))
  X(85) = (X(85)-JVS(461)*X(86)-JVS(462)*X(87)-JVS(463)*X(88)-JVS(464)*X(89)-JVS(465)*X(90)-JVS(466)*X(91)-JVS(467)&
            &*X(92)-JVS(468)*X(93)-JVS(469)*X(94)-JVS(470)*X(95)-JVS(471)*X(96)-JVS(472)*X(97))/(JVS(460))
  X(84) = (X(84)-JVS(438)*X(86)-JVS(439)*X(87)-JVS(440)*X(88)-JVS(441)*X(89)-JVS(442)*X(90)-JVS(443)*X(93)-JVS(444)&
            &*X(95)-JVS(445)*X(96)-JVS(446)*X(97))/(JVS(437))
  X(83) = (X(83)-JVS(415)*X(84)-JVS(416)*X(86)-JVS(417)*X(87)-JVS(418)*X(88)-JVS(419)*X(89)-JVS(420)*X(90)-JVS(421)&
            &*X(91)-JVS(422)*X(92)-JVS(423)*X(93)-JVS(424)*X(94)-JVS(425)*X(95)-JVS(426)*X(96)-JVS(427)*X(97))/(JVS(414))
  X(82) = (X(82)-JVS(393)*X(84)-JVS(394)*X(86)-JVS(395)*X(87)-JVS(396)*X(88)-JVS(397)*X(89)-JVS(398)*X(90)-JVS(399)&
            &*X(93)-JVS(400)*X(95)-JVS(401)*X(96)-JVS(402)*X(97))/(JVS(392))
  X(81) = (X(81)-JVS(373)*X(86)-JVS(374)*X(87)-JVS(375)*X(88)-JVS(376)*X(93)-JVS(377)*X(96)-JVS(378)*X(97))/(JVS(372))
  X(80) = (X(80)-JVS(365)*X(86)-JVS(366)*X(87)-JVS(367)*X(88)-JVS(368)*X(93)-JVS(369)*X(96)-JVS(370)*X(97))/(JVS(364))
  X(79) = (X(79)-JVS(355)*X(80)-JVS(356)*X(81)-JVS(357)*X(84)-JVS(358)*X(86)-JVS(359)*X(87)-JVS(360)*X(88)-JVS(361)&
            &*X(93)-JVS(362)*X(96)-JVS(363)*X(97))/(JVS(354))
  X(78) = (X(78)-JVS(338)*X(79)-JVS(339)*X(80)-JVS(340)*X(86)-JVS(341)*X(87)-JVS(342)*X(88)-JVS(343)*X(89)-JVS(344)&
            &*X(90)-JVS(345)*X(93)-JVS(346)*X(95)-JVS(347)*X(96)-JVS(348)*X(97))/(JVS(337))
  X(77) = (X(77)-JVS(328)*X(86)-JVS(329)*X(88)-JVS(330)*X(93)-JVS(331)*X(96)-JVS(332)*X(97))/(JVS(327))
  X(76) = (X(76)-JVS(321)*X(86)-JVS(322)*X(88)-JVS(323)*X(93)-JVS(324)*X(96)-JVS(325)*X(97))/(JVS(320))
  X(75) = (X(75)-JVS(316)*X(86)-JVS(317)*X(93)-JVS(318)*X(96)-JVS(319)*X(97))/(JVS(315))
  X(74) = (X(74)-JVS(310)*X(86)-JVS(311)*X(88)-JVS(312)*X(93)-JVS(313)*X(96)-JVS(314)*X(97))/(JVS(309))
  X(73) = (X(73)-JVS(304)*X(87)-JVS(305)*X(90)-JVS(306)*X(95)-JVS(307)*X(96)-JVS(308)*X(97))/(JVS(303))
  X(72) = (X(72)-JVS(295)*X(79)-JVS(296)*X(87)-JVS(297)*X(88)-JVS(298)*X(96))/(JVS(294))
  X(71) = (X(71)-JVS(290)*X(73)-JVS(291)*X(86)-JVS(292)*X(95)-JVS(293)*X(96))/(JVS(289))
  X(70) = (X(70)-JVS(273)*X(71)-JVS(274)*X(74)-JVS(275)*X(75)-JVS(276)*X(76)-JVS(277)*X(77)-JVS(278)*X(80)-JVS(279)&
            &*X(81)-JVS(280)*X(82)-JVS(281)*X(83)-JVS(282)*X(85)-JVS(283)*X(86)-JVS(284)*X(88)-JVS(285)*X(93)-JVS(286)*X(96)&
            &-JVS(287)*X(97))/(JVS(272))
  X(69) = (X(69)-JVS(260)*X(73)-JVS(261)*X(81)-JVS(262)*X(82)-JVS(263)*X(83)-JVS(264)*X(84)-JVS(265)*X(85)-JVS(266)&
            &*X(87)-JVS(267)*X(90)-JVS(268)*X(96)-JVS(269)*X(97))/(JVS(259))
  X(68) = (X(68)-JVS(253)*X(85)-JVS(254)*X(90)-JVS(255)*X(92)-JVS(256)*X(95)-JVS(257)*X(96))/(JVS(252))
  X(67) = (X(67)-JVS(245)*X(74)-JVS(246)*X(76)-JVS(247)*X(77)-JVS(248)*X(80)-JVS(249)*X(88)-JVS(250)*X(96))/(JVS(244))
  X(66) = (X(66)-JVS(240)*X(71)-JVS(241)*X(81)-JVS(242)*X(86)-JVS(243)*X(96))/(JVS(239))
  X(65) = (X(65)-JVS(234)*X(78)-JVS(235)*X(89)-JVS(236)*X(90)-JVS(237)*X(96))/(JVS(233))
  X(64) = (X(64)-JVS(229)*X(86)-JVS(230)*X(88)-JVS(231)*X(90)-JVS(232)*X(95))/(JVS(228))
  X(63) = (X(63)-JVS(223)*X(73)-JVS(224)*X(87)-JVS(225)*X(90)-JVS(226)*X(96)-JVS(227)*X(97))/(JVS(222))
  X(62) = (X(62)-JVS(219)*X(85)-JVS(220)*X(90)-JVS(221)*X(95))/(JVS(218))
  X(61) = (X(61)-JVS(209)*X(76)-JVS(210)*X(79)-JVS(211)*X(80)-JVS(212)*X(82)-JVS(213)*X(83)-JVS(214)*X(85)-JVS(215)&
            &*X(88)-JVS(216)*X(92)-JVS(217)*X(96))/(JVS(208))
  X(60) = (X(60)-JVS(200)*X(86)-JVS(201)*X(90)-JVS(202)*X(93)-JVS(203)*X(96))/(JVS(199))
  X(59) = (X(59)-JVS(196)*X(88)-JVS(197)*X(96))/(JVS(195))
  X(58) = (X(58)-JVS(192)*X(88)-JVS(193)*X(92)-JVS(194)*X(96))/(JVS(191))
  X(57) = (X(57)-JVS(189)*X(88)-JVS(190)*X(96))/(JVS(188))
  X(56) = (X(56)-JVS(185)*X(87)-JVS(186)*X(94)-JVS(187)*X(96))/(JVS(184))
  X(55) = (X(55)-JVS(178)*X(89)-JVS(179)*X(90)-JVS(180)*X(91)-JVS(181)*X(92)-JVS(182)*X(94)-JVS(183)*X(96))/(JVS(177))
  X(54) = (X(54)-JVS(173)*X(90)-JVS(174)*X(91)-JVS(175)*X(94)-JVS(176)*X(96))/(JVS(172))
  X(53) = (X(53)-JVS(169)*X(87)-JVS(170)*X(90)-JVS(171)*X(96))/(JVS(168))
  X(52) = (X(52)-JVS(166)*X(95)-JVS(167)*X(96))/(JVS(165))
  X(51) = (X(51)-JVS(158)*X(62)-JVS(159)*X(74)-JVS(160)*X(86)-JVS(161)*X(95)-JVS(162)*X(96))/(JVS(157))
  X(50) = (X(50)-JVS(154)*X(87)-JVS(155)*X(95)-JVS(156)*X(96))/(JVS(153))
  X(49) = (X(49)-JVS(150)*X(60)-JVS(151)*X(86)-JVS(152)*X(96))/(JVS(149))
  X(48) = (X(48)-JVS(147)*X(88)-JVS(148)*X(96))/(JVS(146))
  X(47) = (X(47)-JVS(145)*X(96))/(JVS(144))
  X(46) = (X(46)-JVS(142)*X(87)-JVS(143)*X(97))/(JVS(141))
  X(45) = (X(45)-JVS(139)*X(86)-JVS(140)*X(96))/(JVS(138))
  X(44) = (X(44)-JVS(135)*X(45)-JVS(136)*X(96))/(JVS(134))
  X(43) = (X(43)-JVS(133)*X(96))/(JVS(132))
  X(42) = (X(42)-JVS(130)*X(64)-JVS(131)*X(90))/(JVS(129))
  X(41) = (X(41)-JVS(127)*X(87)-JVS(128)*X(91))/(JVS(126))
  X(40) = (X(40)-JVS(125)*X(96))/(JVS(124))
  X(39) = (X(39)-JVS(123)*X(64))/(JVS(122))
  X(38) = (X(38)-JVS(118)*X(49)-JVS(119)*X(60)-JVS(120)*X(86)-JVS(121)*X(96))/(JVS(117))
  X(37) = (X(37)-JVS(116)*X(96))/(JVS(115))
  X(36) = (X(36)-JVS(113)*X(37)-JVS(114)*X(96))/(JVS(112))
  X(35) = (X(35)-JVS(110)*X(37)-JVS(111)*X(96))/(JVS(109))
  X(34) = (X(34)-JVS(108)*X(96))/(JVS(107))
  X(33) = (X(33)-JVS(105)*X(34)-JVS(106)*X(96))/(JVS(104))
  X(32) = (X(32)-JVS(103)*X(96))/(JVS(102))
  X(31) = (X(31)-JVS(100)*X(32)-JVS(101)*X(96))/(JVS(99))
  X(30) = (X(30)-JVS(97)*X(32)-JVS(98)*X(96))/(JVS(96))
  X(29) = (X(29)-JVS(93)*X(86)-JVS(94)*X(96)-JVS(95)*X(97))/(JVS(92))
  X(28) = (X(28)-JVS(90)*X(29)-JVS(91)*X(97))/(JVS(89))
  X(27) = (X(27)-JVS(87)*X(29)-JVS(88)*X(86))/(JVS(86))
  X(26) = (X(26)-JVS(84)*X(29)-JVS(85)*X(86))/(JVS(83))
  X(25) = (X(25)-JVS(81)*X(29)-JVS(82)*X(96))/(JVS(80))
  X(24) = (X(24)-JVS(78)*X(29)-JVS(79)*X(96))/(JVS(77))
  X(23) = (X(23)-JVS(75)*X(86)-JVS(76)*X(96))/(JVS(74))
  X(22) = (X(22)-JVS(72)*X(23)-JVS(73)*X(86))/(JVS(71))
  X(21) = (X(21)-JVS(69)*X(23)-JVS(70)*X(86))/(JVS(68))
  X(20) = (X(20)-JVS(66)*X(23)-JVS(67)*X(96))/(JVS(65))
  X(19) = (X(19)-JVS(63)*X(23)-JVS(64)*X(96))/(JVS(62))
  X(18) = (X(18)-JVS(61)*X(96))/(JVS(60))
  X(17) = (X(17)-JVS(58)*X(18)-JVS(59)*X(96))/(JVS(57))
  X(16) = (X(16)-JVS(55)*X(18)-JVS(56)*X(96))/(JVS(54))
  X(15) = (X(15)-JVS(53)*X(96))/(JVS(52))
  X(14) = (X(14)-JVS(50)*X(15)-JVS(51)*X(96))/(JVS(49))
  X(13) = (X(13)-JVS(47)*X(15)-JVS(48)*X(96))/(JVS(46))
  X(12) = (X(12)-JVS(45)*X(96))/(JVS(44))
  X(11) = (X(11)-JVS(42)*X(12)-JVS(43)*X(96))/(JVS(41))
  X(10) = (X(10)-JVS(36)*X(75)-JVS(37)*X(86)-JVS(38)*X(93)-JVS(39)*X(96)-JVS(40)*X(97))/(JVS(35))
  X(9) = (X(9)-JVS(33)*X(40)-JVS(34)*X(96))/(JVS(32))
  X(8) = (X(8)-JVS(30)*X(40)-JVS(31)*X(96))/(JVS(29))
  X(7) = (X(7)-JVS(24)*X(80)-JVS(25)*X(86)-JVS(26)*X(93)-JVS(27)*X(96)-JVS(28)*X(97))/(JVS(23))
  X(6) = (X(6)-JVS(18)*X(80)-JVS(19)*X(86)-JVS(20)*X(93)-JVS(21)*X(96)-JVS(22)*X(97))/(JVS(17))
  X(5) = (X(5)-JVS(15)*X(47)-JVS(16)*X(96))/(JVS(14))
  X(4) = (X(4)-JVS(12)*X(47)-JVS(13)*X(96))/(JVS(11))
  X(3) = (X(3)-JVS(8)*X(73)-JVS(9)*X(96)-JVS(10)*X(97))/(JVS(7))
  X(2) = (X(2)-JVS(5)*X(43)-JVS(6)*X(96))/(JVS(4))
  X(1) = (X(1)-JVS(2)*X(43)-JVS(3)*X(96))/(JVS(1))
      
END SUBROUTINE cb05_sorg_aq_KppSolve
























      SUBROUTINE cb05_sorg_aq_WCOPY(N,X,incX,Y,incY)








      
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

      END SUBROUTINE cb05_sorg_aq_WCOPY



      SUBROUTINE cb05_sorg_aq_WAXPY(N,Alpha,X,incX,Y,incY)









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
      
      END SUBROUTINE cb05_sorg_aq_WAXPY




      SUBROUTINE cb05_sorg_aq_WSCAL(N,Alpha,X,incX)









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

      END SUBROUTINE cb05_sorg_aq_WSCAL


      REAL(kind=dp) FUNCTION cb05_sorg_aq_WLAMCH( C )








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
          CALL cb05_sorg_aq_WLAMCH_ADD(ONE,Eps,Sum)
          IF (Sum.LE.ONE) GOTO 10
        END DO
        PRINT*,'ERROR IN WLAMCH. EPS < ',Eps
        RETURN
10      Eps = Eps*2
        i = i-1      
      END IF

      cb05_sorg_aq_WLAMCH = Eps

      END FUNCTION cb05_sorg_aq_WLAMCH
     
      SUBROUTINE cb05_sorg_aq_WLAMCH_ADD( A, B, Sum )

      
      REAL(kind=dp) A, B, Sum
      Sum = A + B

      END SUBROUTINE cb05_sorg_aq_WLAMCH_ADD




      SUBROUTINE cb05_sorg_aq_SET2ZERO(N,Y)




      
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

      END SUBROUTINE cb05_sorg_aq_SET2ZERO



      REAL(kind=dp) FUNCTION cb05_sorg_aq_WDOT (N, DX, incX, DY, incY) 









      IMPLICIT NONE
      INTEGER :: N, incX, incY
      REAL(kind=dp) :: DX(N), DY(N) 

      INTEGER :: i, IX, IY, M, MP1, NS
                                 
      cb05_sorg_aq_WDOT = 0.0D0 
      IF (N .LE. 0) RETURN 
      IF (incX .EQ. incY) IF (incX-1) 5,20,60 



    5 IX = 1 
      IY = 1 
      IF (incX .LT. 0) IX = (-N+1)*incX + 1 
      IF (incY .LT. 0) IY = (-N+1)*incY + 1 
      DO i = 1,N 
        cb05_sorg_aq_WDOT = cb05_sorg_aq_WDOT + DX(IX)*DY(IY) 
        IX = IX + incX 
        IY = IY + incY 
      END DO 
      RETURN 





   20 M = MOD(N,5) 
      IF (M .EQ. 0) GO TO 40 
      DO i = 1,M 
         cb05_sorg_aq_WDOT = cb05_sorg_aq_WDOT + DX(i)*DY(i) 
      END DO 
      IF (N .LT. 5) RETURN 
   40 MP1 = M + 1 
      DO i = MP1,N,5 
          cb05_sorg_aq_WDOT = cb05_sorg_aq_WDOT + DX(i)*DY(i) + DX(i+1)*DY(i+1) +&
                   DX(i+2)*DY(i+2) +  &
                   DX(i+3)*DY(i+3) + DX(i+4)*DY(i+4)                   
      END DO 
      RETURN 



   60 NS = N*incX 
      DO i = 1,NS,incX 
        cb05_sorg_aq_WDOT = cb05_sorg_aq_WDOT + DX(i)*DY(i) 
      END DO 

      END FUNCTION cb05_sorg_aq_WDOT                                          




   SUBROUTINE decomp_cb05_sorg_aq( JVS, IER )
   
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
   W( 43 ) = JVS( 2 )
   W( 96 ) = JVS( 3 )
  JVS( 1) = W( 1 )
  JVS( 2) = W( 43 )
  JVS( 3) = W( 96 )
  IF ( ABS(  JVS( 4 )) < TINY(a) ) THEN
         IER = 2                                       
         RETURN
  END IF
   W( 2 ) = JVS( 4 )
   W( 43 ) = JVS( 5 )
   W( 96 ) = JVS( 6 )
  JVS( 4) = W( 2 )
  JVS( 5) = W( 43 )
  JVS( 6) = W( 96 )
  IF ( ABS(  JVS( 7 )) < TINY(a) ) THEN
         IER = 3                                       
         RETURN
  END IF
   W( 3 ) = JVS( 7 )
   W( 73 ) = JVS( 8 )
   W( 96 ) = JVS( 9 )
   W( 97 ) = JVS( 10 )
  JVS( 7) = W( 3 )
  JVS( 8) = W( 73 )
  JVS( 9) = W( 96 )
  JVS( 10) = W( 97 )
  IF ( ABS(  JVS( 11 )) < TINY(a) ) THEN
         IER = 4                                       
         RETURN
  END IF
   W( 4 ) = JVS( 11 )
   W( 47 ) = JVS( 12 )
   W( 96 ) = JVS( 13 )
  JVS( 11) = W( 4 )
  JVS( 12) = W( 47 )
  JVS( 13) = W( 96 )
  IF ( ABS(  JVS( 14 )) < TINY(a) ) THEN
         IER = 5                                       
         RETURN
  END IF
   W( 5 ) = JVS( 14 )
   W( 47 ) = JVS( 15 )
   W( 96 ) = JVS( 16 )
  JVS( 14) = W( 5 )
  JVS( 15) = W( 47 )
  JVS( 16) = W( 96 )
  IF ( ABS(  JVS( 17 )) < TINY(a) ) THEN
         IER = 6                                       
         RETURN
  END IF
   W( 6 ) = JVS( 17 )
   W( 80 ) = JVS( 18 )
   W( 86 ) = JVS( 19 )
   W( 93 ) = JVS( 20 )
   W( 96 ) = JVS( 21 )
   W( 97 ) = JVS( 22 )
  JVS( 17) = W( 6 )
  JVS( 18) = W( 80 )
  JVS( 19) = W( 86 )
  JVS( 20) = W( 93 )
  JVS( 21) = W( 96 )
  JVS( 22) = W( 97 )
  IF ( ABS(  JVS( 23 )) < TINY(a) ) THEN
         IER = 7                                       
         RETURN
  END IF
   W( 7 ) = JVS( 23 )
   W( 80 ) = JVS( 24 )
   W( 86 ) = JVS( 25 )
   W( 93 ) = JVS( 26 )
   W( 96 ) = JVS( 27 )
   W( 97 ) = JVS( 28 )
  JVS( 23) = W( 7 )
  JVS( 24) = W( 80 )
  JVS( 25) = W( 86 )
  JVS( 26) = W( 93 )
  JVS( 27) = W( 96 )
  JVS( 28) = W( 97 )
  IF ( ABS(  JVS( 29 )) < TINY(a) ) THEN
         IER = 8                                       
         RETURN
  END IF
   W( 8 ) = JVS( 29 )
   W( 40 ) = JVS( 30 )
   W( 96 ) = JVS( 31 )
  JVS( 29) = W( 8 )
  JVS( 30) = W( 40 )
  JVS( 31) = W( 96 )
  IF ( ABS(  JVS( 32 )) < TINY(a) ) THEN
         IER = 9                                       
         RETURN
  END IF
   W( 9 ) = JVS( 32 )
   W( 40 ) = JVS( 33 )
   W( 96 ) = JVS( 34 )
  JVS( 32) = W( 9 )
  JVS( 33) = W( 40 )
  JVS( 34) = W( 96 )
  IF ( ABS(  JVS( 35 )) < TINY(a) ) THEN
         IER = 10                                      
         RETURN
  END IF
   W( 10 ) = JVS( 35 )
   W( 75 ) = JVS( 36 )
   W( 86 ) = JVS( 37 )
   W( 93 ) = JVS( 38 )
   W( 96 ) = JVS( 39 )
   W( 97 ) = JVS( 40 )
  JVS( 35) = W( 10 )
  JVS( 36) = W( 75 )
  JVS( 37) = W( 86 )
  JVS( 38) = W( 93 )
  JVS( 39) = W( 96 )
  JVS( 40) = W( 97 )
  IF ( ABS(  JVS( 41 )) < TINY(a) ) THEN
         IER = 11                                      
         RETURN
  END IF
   W( 11 ) = JVS( 41 )
   W( 12 ) = JVS( 42 )
   W( 96 ) = JVS( 43 )
  JVS( 41) = W( 11 )
  JVS( 42) = W( 12 )
  JVS( 43) = W( 96 )
  IF ( ABS(  JVS( 44 )) < TINY(a) ) THEN
         IER = 12                                      
         RETURN
  END IF
   W( 12 ) = JVS( 44 )
   W( 96 ) = JVS( 45 )
  JVS( 44) = W( 12 )
  JVS( 45) = W( 96 )
  IF ( ABS(  JVS( 46 )) < TINY(a) ) THEN
         IER = 13                                      
         RETURN
  END IF
   W( 13 ) = JVS( 46 )
   W( 15 ) = JVS( 47 )
   W( 96 ) = JVS( 48 )
  JVS( 46) = W( 13 )
  JVS( 47) = W( 15 )
  JVS( 48) = W( 96 )
  IF ( ABS(  JVS( 49 )) < TINY(a) ) THEN
         IER = 14                                      
         RETURN
  END IF
   W( 14 ) = JVS( 49 )
   W( 15 ) = JVS( 50 )
   W( 96 ) = JVS( 51 )
  JVS( 49) = W( 14 )
  JVS( 50) = W( 15 )
  JVS( 51) = W( 96 )
  IF ( ABS(  JVS( 52 )) < TINY(a) ) THEN
         IER = 15                                      
         RETURN
  END IF
   W( 15 ) = JVS( 52 )
   W( 96 ) = JVS( 53 )
  JVS( 52) = W( 15 )
  JVS( 53) = W( 96 )
  IF ( ABS(  JVS( 54 )) < TINY(a) ) THEN
         IER = 16                                      
         RETURN
  END IF
   W( 16 ) = JVS( 54 )
   W( 18 ) = JVS( 55 )
   W( 96 ) = JVS( 56 )
  JVS( 54) = W( 16 )
  JVS( 55) = W( 18 )
  JVS( 56) = W( 96 )
  IF ( ABS(  JVS( 57 )) < TINY(a) ) THEN
         IER = 17                                      
         RETURN
  END IF
   W( 17 ) = JVS( 57 )
   W( 18 ) = JVS( 58 )
   W( 96 ) = JVS( 59 )
  JVS( 57) = W( 17 )
  JVS( 58) = W( 18 )
  JVS( 59) = W( 96 )
  IF ( ABS(  JVS( 60 )) < TINY(a) ) THEN
         IER = 18                                      
         RETURN
  END IF
   W( 18 ) = JVS( 60 )
   W( 96 ) = JVS( 61 )
  JVS( 60) = W( 18 )
  JVS( 61) = W( 96 )
  IF ( ABS(  JVS( 62 )) < TINY(a) ) THEN
         IER = 19                                      
         RETURN
  END IF
   W( 19 ) = JVS( 62 )
   W( 23 ) = JVS( 63 )
   W( 96 ) = JVS( 64 )
  JVS( 62) = W( 19 )
  JVS( 63) = W( 23 )
  JVS( 64) = W( 96 )
  IF ( ABS(  JVS( 65 )) < TINY(a) ) THEN
         IER = 20                                      
         RETURN
  END IF
   W( 20 ) = JVS( 65 )
   W( 23 ) = JVS( 66 )
   W( 96 ) = JVS( 67 )
  JVS( 65) = W( 20 )
  JVS( 66) = W( 23 )
  JVS( 67) = W( 96 )
  IF ( ABS(  JVS( 68 )) < TINY(a) ) THEN
         IER = 21                                      
         RETURN
  END IF
   W( 21 ) = JVS( 68 )
   W( 23 ) = JVS( 69 )
   W( 86 ) = JVS( 70 )
  JVS( 68) = W( 21 )
  JVS( 69) = W( 23 )
  JVS( 70) = W( 86 )
  IF ( ABS(  JVS( 71 )) < TINY(a) ) THEN
         IER = 22                                      
         RETURN
  END IF
   W( 22 ) = JVS( 71 )
   W( 23 ) = JVS( 72 )
   W( 86 ) = JVS( 73 )
  JVS( 71) = W( 22 )
  JVS( 72) = W( 23 )
  JVS( 73) = W( 86 )
  IF ( ABS(  JVS( 74 )) < TINY(a) ) THEN
         IER = 23                                      
         RETURN
  END IF
   W( 23 ) = JVS( 74 )
   W( 86 ) = JVS( 75 )
   W( 96 ) = JVS( 76 )
  JVS( 74) = W( 23 )
  JVS( 75) = W( 86 )
  JVS( 76) = W( 96 )
  IF ( ABS(  JVS( 77 )) < TINY(a) ) THEN
         IER = 24                                      
         RETURN
  END IF
   W( 24 ) = JVS( 77 )
   W( 29 ) = JVS( 78 )
   W( 96 ) = JVS( 79 )
  JVS( 77) = W( 24 )
  JVS( 78) = W( 29 )
  JVS( 79) = W( 96 )
  IF ( ABS(  JVS( 80 )) < TINY(a) ) THEN
         IER = 25                                      
         RETURN
  END IF
   W( 25 ) = JVS( 80 )
   W( 29 ) = JVS( 81 )
   W( 96 ) = JVS( 82 )
  JVS( 80) = W( 25 )
  JVS( 81) = W( 29 )
  JVS( 82) = W( 96 )
  IF ( ABS(  JVS( 83 )) < TINY(a) ) THEN
         IER = 26                                      
         RETURN
  END IF
   W( 26 ) = JVS( 83 )
   W( 29 ) = JVS( 84 )
   W( 86 ) = JVS( 85 )
  JVS( 83) = W( 26 )
  JVS( 84) = W( 29 )
  JVS( 85) = W( 86 )
  IF ( ABS(  JVS( 86 )) < TINY(a) ) THEN
         IER = 27                                      
         RETURN
  END IF
   W( 27 ) = JVS( 86 )
   W( 29 ) = JVS( 87 )
   W( 86 ) = JVS( 88 )
  JVS( 86) = W( 27 )
  JVS( 87) = W( 29 )
  JVS( 88) = W( 86 )
  IF ( ABS(  JVS( 89 )) < TINY(a) ) THEN
         IER = 28                                      
         RETURN
  END IF
   W( 28 ) = JVS( 89 )
   W( 29 ) = JVS( 90 )
   W( 97 ) = JVS( 91 )
  JVS( 89) = W( 28 )
  JVS( 90) = W( 29 )
  JVS( 91) = W( 97 )
  IF ( ABS(  JVS( 92 )) < TINY(a) ) THEN
         IER = 29                                      
         RETURN
  END IF
   W( 29 ) = JVS( 92 )
   W( 86 ) = JVS( 93 )
   W( 96 ) = JVS( 94 )
   W( 97 ) = JVS( 95 )
  JVS( 92) = W( 29 )
  JVS( 93) = W( 86 )
  JVS( 94) = W( 96 )
  JVS( 95) = W( 97 )
  IF ( ABS(  JVS( 96 )) < TINY(a) ) THEN
         IER = 30                                      
         RETURN
  END IF
   W( 30 ) = JVS( 96 )
   W( 32 ) = JVS( 97 )
   W( 96 ) = JVS( 98 )
  JVS( 96) = W( 30 )
  JVS( 97) = W( 32 )
  JVS( 98) = W( 96 )
  IF ( ABS(  JVS( 99 )) < TINY(a) ) THEN
         IER = 31                                      
         RETURN
  END IF
   W( 31 ) = JVS( 99 )
   W( 32 ) = JVS( 100 )
   W( 96 ) = JVS( 101 )
  JVS( 99) = W( 31 )
  JVS( 100) = W( 32 )
  JVS( 101) = W( 96 )
  IF ( ABS(  JVS( 102 )) < TINY(a) ) THEN
         IER = 32                                      
         RETURN
  END IF
   W( 32 ) = JVS( 102 )
   W( 96 ) = JVS( 103 )
  JVS( 102) = W( 32 )
  JVS( 103) = W( 96 )
  IF ( ABS(  JVS( 104 )) < TINY(a) ) THEN
         IER = 33                                      
         RETURN
  END IF
   W( 33 ) = JVS( 104 )
   W( 34 ) = JVS( 105 )
   W( 96 ) = JVS( 106 )
  JVS( 104) = W( 33 )
  JVS( 105) = W( 34 )
  JVS( 106) = W( 96 )
  IF ( ABS(  JVS( 107 )) < TINY(a) ) THEN
         IER = 34                                      
         RETURN
  END IF
   W( 34 ) = JVS( 107 )
   W( 96 ) = JVS( 108 )
  JVS( 107) = W( 34 )
  JVS( 108) = W( 96 )
  IF ( ABS(  JVS( 109 )) < TINY(a) ) THEN
         IER = 35                                      
         RETURN
  END IF
   W( 35 ) = JVS( 109 )
   W( 37 ) = JVS( 110 )
   W( 96 ) = JVS( 111 )
  JVS( 109) = W( 35 )
  JVS( 110) = W( 37 )
  JVS( 111) = W( 96 )
  IF ( ABS(  JVS( 112 )) < TINY(a) ) THEN
         IER = 36                                      
         RETURN
  END IF
   W( 36 ) = JVS( 112 )
   W( 37 ) = JVS( 113 )
   W( 96 ) = JVS( 114 )
  JVS( 112) = W( 36 )
  JVS( 113) = W( 37 )
  JVS( 114) = W( 96 )
  IF ( ABS(  JVS( 115 )) < TINY(a) ) THEN
         IER = 37                                      
         RETURN
  END IF
   W( 37 ) = JVS( 115 )
   W( 96 ) = JVS( 116 )
  JVS( 115) = W( 37 )
  JVS( 116) = W( 96 )
  IF ( ABS(  JVS( 117 )) < TINY(a) ) THEN
         IER = 38                                      
         RETURN
  END IF
   W( 38 ) = JVS( 117 )
   W( 49 ) = JVS( 118 )
   W( 60 ) = JVS( 119 )
   W( 86 ) = JVS( 120 )
   W( 96 ) = JVS( 121 )
  JVS( 117) = W( 38 )
  JVS( 118) = W( 49 )
  JVS( 119) = W( 60 )
  JVS( 120) = W( 86 )
  JVS( 121) = W( 96 )
  IF ( ABS(  JVS( 122 )) < TINY(a) ) THEN
         IER = 39                                      
         RETURN
  END IF
   W( 39 ) = JVS( 122 )
   W( 64 ) = JVS( 123 )
  JVS( 122) = W( 39 )
  JVS( 123) = W( 64 )
  IF ( ABS(  JVS( 124 )) < TINY(a) ) THEN
         IER = 40                                      
         RETURN
  END IF
   W( 40 ) = JVS( 124 )
   W( 96 ) = JVS( 125 )
  JVS( 124) = W( 40 )
  JVS( 125) = W( 96 )
  IF ( ABS(  JVS( 126 )) < TINY(a) ) THEN
         IER = 41                                      
         RETURN
  END IF
   W( 41 ) = JVS( 126 )
   W( 87 ) = JVS( 127 )
   W( 91 ) = JVS( 128 )
  JVS( 126) = W( 41 )
  JVS( 127) = W( 87 )
  JVS( 128) = W( 91 )
  IF ( ABS(  JVS( 129 )) < TINY(a) ) THEN
         IER = 42                                      
         RETURN
  END IF
   W( 42 ) = JVS( 129 )
   W( 64 ) = JVS( 130 )
   W( 90 ) = JVS( 131 )
  JVS( 129) = W( 42 )
  JVS( 130) = W( 64 )
  JVS( 131) = W( 90 )
  IF ( ABS(  JVS( 132 )) < TINY(a) ) THEN
         IER = 43                                      
         RETURN
  END IF
   W( 43 ) = JVS( 132 )
   W( 96 ) = JVS( 133 )
  JVS( 132) = W( 43 )
  JVS( 133) = W( 96 )
  IF ( ABS(  JVS( 134 )) < TINY(a) ) THEN
         IER = 44                                      
         RETURN
  END IF
   W( 44 ) = JVS( 134 )
   W( 45 ) = JVS( 135 )
   W( 96 ) = JVS( 136 )
  JVS( 134) = W( 44 )
  JVS( 135) = W( 45 )
  JVS( 136) = W( 96 )
  IF ( ABS(  JVS( 138 )) < TINY(a) ) THEN
         IER = 45                                      
         RETURN
  END IF
   W( 44 ) = JVS( 137 )
   W( 45 ) = JVS( 138 )
   W( 86 ) = JVS( 139 )
   W( 96 ) = JVS( 140 )
  a = -W( 44 ) / JVS(          134  )
  W( 44 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 135 )
  W( 96 ) = W( 96 ) + a*JVS( 136 )
  JVS( 137) = W( 44 )
  JVS( 138) = W( 45 )
  JVS( 139) = W( 86 )
  JVS( 140) = W( 96 )
  IF ( ABS(  JVS( 141 )) < TINY(a) ) THEN
         IER = 46                                      
         RETURN
  END IF
   W( 46 ) = JVS( 141 )
   W( 87 ) = JVS( 142 )
   W( 97 ) = JVS( 143 )
  JVS( 141) = W( 46 )
  JVS( 142) = W( 87 )
  JVS( 143) = W( 97 )
  IF ( ABS(  JVS( 144 )) < TINY(a) ) THEN
         IER = 47                                      
         RETURN
  END IF
   W( 47 ) = JVS( 144 )
   W( 96 ) = JVS( 145 )
  JVS( 144) = W( 47 )
  JVS( 145) = W( 96 )
  IF ( ABS(  JVS( 146 )) < TINY(a) ) THEN
         IER = 48                                      
         RETURN
  END IF
   W( 48 ) = JVS( 146 )
   W( 88 ) = JVS( 147 )
   W( 96 ) = JVS( 148 )
  JVS( 146) = W( 48 )
  JVS( 147) = W( 88 )
  JVS( 148) = W( 96 )
  IF ( ABS(  JVS( 149 )) < TINY(a) ) THEN
         IER = 49                                      
         RETURN
  END IF
   W( 49 ) = JVS( 149 )
   W( 60 ) = JVS( 150 )
   W( 86 ) = JVS( 151 )
   W( 96 ) = JVS( 152 )
  JVS( 149) = W( 49 )
  JVS( 150) = W( 60 )
  JVS( 151) = W( 86 )
  JVS( 152) = W( 96 )
  IF ( ABS(  JVS( 153 )) < TINY(a) ) THEN
         IER = 50                                      
         RETURN
  END IF
   W( 50 ) = JVS( 153 )
   W( 87 ) = JVS( 154 )
   W( 95 ) = JVS( 155 )
   W( 96 ) = JVS( 156 )
  JVS( 153) = W( 50 )
  JVS( 154) = W( 87 )
  JVS( 155) = W( 95 )
  JVS( 156) = W( 96 )
  IF ( ABS(  JVS( 157 )) < TINY(a) ) THEN
         IER = 51                                      
         RETURN
  END IF
   W( 51 ) = JVS( 157 )
   W( 62 ) = JVS( 158 )
   W( 74 ) = JVS( 159 )
   W( 86 ) = JVS( 160 )
   W( 95 ) = JVS( 161 )
   W( 96 ) = JVS( 162 )
  JVS( 157) = W( 51 )
  JVS( 158) = W( 62 )
  JVS( 159) = W( 74 )
  JVS( 160) = W( 86 )
  JVS( 161) = W( 95 )
  JVS( 162) = W( 96 )
  IF ( ABS(  JVS( 165 )) < TINY(a) ) THEN
         IER = 52                                      
         RETURN
  END IF
   W( 43 ) = JVS( 163 )
   W( 47 ) = JVS( 164 )
   W( 52 ) = JVS( 165 )
   W( 95 ) = JVS( 166 )
   W( 96 ) = JVS( 167 )
  a = -W( 43 ) / JVS(          132  )
  W( 43 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 133 )
  a = -W( 47 ) / JVS(          144  )
  W( 47 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 145 )
  JVS( 163) = W( 43 )
  JVS( 164) = W( 47 )
  JVS( 165) = W( 52 )
  JVS( 166) = W( 95 )
  JVS( 167) = W( 96 )
  IF ( ABS(  JVS( 168 )) < TINY(a) ) THEN
         IER = 53                                      
         RETURN
  END IF
   W( 53 ) = JVS( 168 )
   W( 87 ) = JVS( 169 )
   W( 90 ) = JVS( 170 )
   W( 96 ) = JVS( 171 )
  JVS( 168) = W( 53 )
  JVS( 169) = W( 87 )
  JVS( 170) = W( 90 )
  JVS( 171) = W( 96 )
  IF ( ABS(  JVS( 172 )) < TINY(a) ) THEN
         IER = 54                                      
         RETURN
  END IF
   W( 54 ) = JVS( 172 )
   W( 90 ) = JVS( 173 )
   W( 91 ) = JVS( 174 )
   W( 94 ) = JVS( 175 )
   W( 96 ) = JVS( 176 )
  JVS( 172) = W( 54 )
  JVS( 173) = W( 90 )
  JVS( 174) = W( 91 )
  JVS( 175) = W( 94 )
  JVS( 176) = W( 96 )
  IF ( ABS(  JVS( 177 )) < TINY(a) ) THEN
         IER = 55                                      
         RETURN
  END IF
   W( 55 ) = JVS( 177 )
   W( 89 ) = JVS( 178 )
   W( 90 ) = JVS( 179 )
   W( 91 ) = JVS( 180 )
   W( 92 ) = JVS( 181 )
   W( 94 ) = JVS( 182 )
   W( 96 ) = JVS( 183 )
  JVS( 177) = W( 55 )
  JVS( 178) = W( 89 )
  JVS( 179) = W( 90 )
  JVS( 180) = W( 91 )
  JVS( 181) = W( 92 )
  JVS( 182) = W( 94 )
  JVS( 183) = W( 96 )
  IF ( ABS(  JVS( 184 )) < TINY(a) ) THEN
         IER = 56                                      
         RETURN
  END IF
   W( 56 ) = JVS( 184 )
   W( 87 ) = JVS( 185 )
   W( 94 ) = JVS( 186 )
   W( 96 ) = JVS( 187 )
  JVS( 184) = W( 56 )
  JVS( 185) = W( 87 )
  JVS( 186) = W( 94 )
  JVS( 187) = W( 96 )
  IF ( ABS(  JVS( 188 )) < TINY(a) ) THEN
         IER = 57                                      
         RETURN
  END IF
   W( 57 ) = JVS( 188 )
   W( 88 ) = JVS( 189 )
   W( 96 ) = JVS( 190 )
  JVS( 188) = W( 57 )
  JVS( 189) = W( 88 )
  JVS( 190) = W( 96 )
  IF ( ABS(  JVS( 191 )) < TINY(a) ) THEN
         IER = 58                                      
         RETURN
  END IF
   W( 58 ) = JVS( 191 )
   W( 88 ) = JVS( 192 )
   W( 92 ) = JVS( 193 )
   W( 96 ) = JVS( 194 )
  JVS( 191) = W( 58 )
  JVS( 192) = W( 88 )
  JVS( 193) = W( 92 )
  JVS( 194) = W( 96 )
  IF ( ABS(  JVS( 195 )) < TINY(a) ) THEN
         IER = 59                                      
         RETURN
  END IF
   W( 59 ) = JVS( 195 )
   W( 88 ) = JVS( 196 )
   W( 96 ) = JVS( 197 )
  JVS( 195) = W( 59 )
  JVS( 196) = W( 88 )
  JVS( 197) = W( 96 )
  IF ( ABS(  JVS( 199 )) < TINY(a) ) THEN
         IER = 60                                      
         RETURN
  END IF
   W( 49 ) = JVS( 198 )
   W( 60 ) = JVS( 199 )
   W( 86 ) = JVS( 200 )
   W( 90 ) = JVS( 201 )
   W( 93 ) = JVS( 202 )
   W( 96 ) = JVS( 203 )
  a = -W( 49 ) / JVS(          149  )
  W( 49 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 150 )
  W( 86 ) = W( 86 ) + a*JVS( 151 )
  W( 96 ) = W( 96 ) + a*JVS( 152 )
  JVS( 198) = W( 49 )
  JVS( 199) = W( 60 )
  JVS( 200) = W( 86 )
  JVS( 201) = W( 90 )
  JVS( 202) = W( 93 )
  JVS( 203) = W( 96 )
  IF ( ABS(  JVS( 208 )) < TINY(a) ) THEN
         IER = 61                                      
         RETURN
  END IF
   W( 48 ) = JVS( 204 )
   W( 57 ) = JVS( 205 )
   W( 58 ) = JVS( 206 )
   W( 59 ) = JVS( 207 )
   W( 61 ) = JVS( 208 )
   W( 76 ) = JVS( 209 )
   W( 79 ) = JVS( 210 )
   W( 80 ) = JVS( 211 )
   W( 82 ) = JVS( 212 )
   W( 83 ) = JVS( 213 )
   W( 85 ) = JVS( 214 )
   W( 88 ) = JVS( 215 )
   W( 92 ) = JVS( 216 )
   W( 96 ) = JVS( 217 )
  a = -W( 48 ) / JVS(          146  )
  W( 48 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 147 )
  W( 96 ) = W( 96 ) + a*JVS( 148 )
  a = -W( 57 ) / JVS(          188  )
  W( 57 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 189 )
  W( 96 ) = W( 96 ) + a*JVS( 190 )
  a = -W( 58 ) / JVS(          191  )
  W( 58 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 192 )
  W( 92 ) = W( 92 ) + a*JVS( 193 )
  W( 96 ) = W( 96 ) + a*JVS( 194 )
  a = -W( 59 ) / JVS(          195  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 196 )
  W( 96 ) = W( 96 ) + a*JVS( 197 )
  JVS( 204) = W( 48 )
  JVS( 205) = W( 57 )
  JVS( 206) = W( 58 )
  JVS( 207) = W( 59 )
  JVS( 208) = W( 61 )
  JVS( 209) = W( 76 )
  JVS( 210) = W( 79 )
  JVS( 211) = W( 80 )
  JVS( 212) = W( 82 )
  JVS( 213) = W( 83 )
  JVS( 214) = W( 85 )
  JVS( 215) = W( 88 )
  JVS( 216) = W( 92 )
  JVS( 217) = W( 96 )
  IF ( ABS(  JVS( 218 )) < TINY(a) ) THEN
         IER = 62                                      
         RETURN
  END IF
   W( 62 ) = JVS( 218 )
   W( 85 ) = JVS( 219 )
   W( 90 ) = JVS( 220 )
   W( 95 ) = JVS( 221 )
  JVS( 218) = W( 62 )
  JVS( 219) = W( 85 )
  JVS( 220) = W( 90 )
  JVS( 221) = W( 95 )
  IF ( ABS(  JVS( 222 )) < TINY(a) ) THEN
         IER = 63                                      
         RETURN
  END IF
   W( 63 ) = JVS( 222 )
   W( 73 ) = JVS( 223 )
   W( 87 ) = JVS( 224 )
   W( 90 ) = JVS( 225 )
   W( 96 ) = JVS( 226 )
   W( 97 ) = JVS( 227 )
  JVS( 222) = W( 63 )
  JVS( 223) = W( 73 )
  JVS( 224) = W( 87 )
  JVS( 225) = W( 90 )
  JVS( 226) = W( 96 )
  JVS( 227) = W( 97 )
  IF ( ABS(  JVS( 228 )) < TINY(a) ) THEN
         IER = 64                                      
         RETURN
  END IF
   W( 64 ) = JVS( 228 )
   W( 86 ) = JVS( 229 )
   W( 88 ) = JVS( 230 )
   W( 90 ) = JVS( 231 )
   W( 95 ) = JVS( 232 )
  JVS( 228) = W( 64 )
  JVS( 229) = W( 86 )
  JVS( 230) = W( 88 )
  JVS( 231) = W( 90 )
  JVS( 232) = W( 95 )
  IF ( ABS(  JVS( 233 )) < TINY(a) ) THEN
         IER = 65                                      
         RETURN
  END IF
   W( 65 ) = JVS( 233 )
   W( 78 ) = JVS( 234 )
   W( 89 ) = JVS( 235 )
   W( 90 ) = JVS( 236 )
   W( 96 ) = JVS( 237 )
  JVS( 233) = W( 65 )
  JVS( 234) = W( 78 )
  JVS( 235) = W( 89 )
  JVS( 236) = W( 90 )
  JVS( 237) = W( 96 )
  IF ( ABS(  JVS( 239 )) < TINY(a) ) THEN
         IER = 66                                      
         RETURN
  END IF
   W( 47 ) = JVS( 238 )
   W( 66 ) = JVS( 239 )
   W( 71 ) = JVS( 240 )
   W( 81 ) = JVS( 241 )
   W( 86 ) = JVS( 242 )
   W( 96 ) = JVS( 243 )
  a = -W( 47 ) / JVS(          144  )
  W( 47 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 145 )
  JVS( 238) = W( 47 )
  JVS( 239) = W( 66 )
  JVS( 240) = W( 71 )
  JVS( 241) = W( 81 )
  JVS( 242) = W( 86 )
  JVS( 243) = W( 96 )
  IF ( ABS(  JVS( 244 )) < TINY(a) ) THEN
         IER = 67                                      
         RETURN
  END IF
   W( 67 ) = JVS( 244 )
   W( 74 ) = JVS( 245 )
   W( 76 ) = JVS( 246 )
   W( 77 ) = JVS( 247 )
   W( 80 ) = JVS( 248 )
   W( 88 ) = JVS( 249 )
   W( 96 ) = JVS( 250 )
  JVS( 244) = W( 67 )
  JVS( 245) = W( 74 )
  JVS( 246) = W( 76 )
  JVS( 247) = W( 77 )
  JVS( 248) = W( 80 )
  JVS( 249) = W( 88 )
  JVS( 250) = W( 96 )
  IF ( ABS(  JVS( 252 )) < TINY(a) ) THEN
         IER = 68                                      
         RETURN
  END IF
   W( 62 ) = JVS( 251 )
   W( 68 ) = JVS( 252 )
   W( 85 ) = JVS( 253 )
   W( 90 ) = JVS( 254 )
   W( 92 ) = JVS( 255 )
   W( 95 ) = JVS( 256 )
   W( 96 ) = JVS( 257 )
  a = -W( 62 ) / JVS(          218  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 219 )
  W( 90 ) = W( 90 ) + a*JVS( 220 )
  W( 95 ) = W( 95 ) + a*JVS( 221 )
  JVS( 251) = W( 62 )
  JVS( 252) = W( 68 )
  JVS( 253) = W( 85 )
  JVS( 254) = W( 90 )
  JVS( 255) = W( 92 )
  JVS( 256) = W( 95 )
  JVS( 257) = W( 96 )
  IF ( ABS(  JVS( 259 )) < TINY(a) ) THEN
         IER = 69                                      
         RETURN
  END IF
   W( 46 ) = JVS( 258 )
   W( 69 ) = JVS( 259 )
   W( 73 ) = JVS( 260 )
   W( 81 ) = JVS( 261 )
   W( 82 ) = JVS( 262 )
   W( 83 ) = JVS( 263 )
   W( 84 ) = JVS( 264 )
   W( 85 ) = JVS( 265 )
   W( 87 ) = JVS( 266 )
   W( 90 ) = JVS( 267 )
   W( 96 ) = JVS( 268 )
   W( 97 ) = JVS( 269 )
  a = -W( 46 ) / JVS(          141  )
  W( 46 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 142 )
  W( 97 ) = W( 97 ) + a*JVS( 143 )
  JVS( 258) = W( 46 )
  JVS( 259) = W( 69 )
  JVS( 260) = W( 73 )
  JVS( 261) = W( 81 )
  JVS( 262) = W( 82 )
  JVS( 263) = W( 83 )
  JVS( 264) = W( 84 )
  JVS( 265) = W( 85 )
  JVS( 266) = W( 87 )
  JVS( 267) = W( 90 )
  JVS( 268) = W( 96 )
  JVS( 269) = W( 97 )
  IF ( ABS(  JVS( 272 )) < TINY(a) ) THEN
         IER = 70                                      
         RETURN
  END IF
   W( 66 ) = JVS( 270 )
   W( 67 ) = JVS( 271 )
   W( 70 ) = JVS( 272 )
   W( 71 ) = JVS( 273 )
   W( 74 ) = JVS( 274 )
   W( 75 ) = JVS( 275 )
   W( 76 ) = JVS( 276 )
   W( 77 ) = JVS( 277 )
   W( 80 ) = JVS( 278 )
   W( 81 ) = JVS( 279 )
   W( 82 ) = JVS( 280 )
   W( 83 ) = JVS( 281 )
   W( 85 ) = JVS( 282 )
   W( 86 ) = JVS( 283 )
   W( 88 ) = JVS( 284 )
   W( 93 ) = JVS( 285 )
   W( 96 ) = JVS( 286 )
   W( 97 ) = JVS( 287 )
  a = -W( 66 ) / JVS(          239  )
  W( 66 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 240 )
  W( 81 ) = W( 81 ) + a*JVS( 241 )
  W( 86 ) = W( 86 ) + a*JVS( 242 )
  W( 96 ) = W( 96 ) + a*JVS( 243 )
  a = -W( 67 ) / JVS(          244  )
  W( 67 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 245 )
  W( 76 ) = W( 76 ) + a*JVS( 246 )
  W( 77 ) = W( 77 ) + a*JVS( 247 )
  W( 80 ) = W( 80 ) + a*JVS( 248 )
  W( 88 ) = W( 88 ) + a*JVS( 249 )
  W( 96 ) = W( 96 ) + a*JVS( 250 )
  JVS( 270) = W( 66 )
  JVS( 271) = W( 67 )
  JVS( 272) = W( 70 )
  JVS( 273) = W( 71 )
  JVS( 274) = W( 74 )
  JVS( 275) = W( 75 )
  JVS( 276) = W( 76 )
  JVS( 277) = W( 77 )
  JVS( 278) = W( 80 )
  JVS( 279) = W( 81 )
  JVS( 280) = W( 82 )
  JVS( 281) = W( 83 )
  JVS( 282) = W( 85 )
  JVS( 283) = W( 86 )
  JVS( 284) = W( 88 )
  JVS( 285) = W( 93 )
  JVS( 286) = W( 96 )
  JVS( 287) = W( 97 )
  IF ( ABS(  JVS( 289 )) < TINY(a) ) THEN
         IER = 71                                      
         RETURN
  END IF
   W( 52 ) = JVS( 288 )
   W( 71 ) = JVS( 289 )
   W( 73 ) = JVS( 290 )
   W( 86 ) = JVS( 291 )
   W( 95 ) = JVS( 292 )
   W( 96 ) = JVS( 293 )
  a = -W( 52 ) / JVS(          165  )
  W( 52 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 166 )
  W( 96 ) = W( 96 ) + a*JVS( 167 )
  JVS( 288) = W( 52 )
  JVS( 289) = W( 71 )
  JVS( 290) = W( 73 )
  JVS( 291) = W( 86 )
  JVS( 292) = W( 95 )
  JVS( 293) = W( 96 )
  IF ( ABS(  JVS( 294 )) < TINY(a) ) THEN
         IER = 72                                      
         RETURN
  END IF
   W( 72 ) = JVS( 294 )
   W( 79 ) = JVS( 295 )
   W( 87 ) = JVS( 296 )
   W( 88 ) = JVS( 297 )
   W( 96 ) = JVS( 298 )
  JVS( 294) = W( 72 )
  JVS( 295) = W( 79 )
  JVS( 296) = W( 87 )
  JVS( 297) = W( 88 )
  JVS( 298) = W( 96 )
  IF ( ABS(  JVS( 303 )) < TINY(a) ) THEN
         IER = 73                                      
         RETURN
  END IF
   W( 43 ) = JVS( 299 )
   W( 47 ) = JVS( 300 )
   W( 52 ) = JVS( 301 )
   W( 63 ) = JVS( 302 )
   W( 73 ) = JVS( 303 )
   W( 87 ) = JVS( 304 )
   W( 90 ) = JVS( 305 )
   W( 95 ) = JVS( 306 )
   W( 96 ) = JVS( 307 )
   W( 97 ) = JVS( 308 )
  a = -W( 43 ) / JVS(          132  )
  W( 43 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 133 )
  a = -W( 47 ) / JVS(          144  )
  W( 47 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 145 )
  a = -W( 52 ) / JVS(          165  )
  W( 52 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 166 )
  W( 96 ) = W( 96 ) + a*JVS( 167 )
  a = -W( 63 ) / JVS(          222  )
  W( 63 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 223 )
  W( 87 ) = W( 87 ) + a*JVS( 224 )
  W( 90 ) = W( 90 ) + a*JVS( 225 )
  W( 96 ) = W( 96 ) + a*JVS( 226 )
  W( 97 ) = W( 97 ) + a*JVS( 227 )
  JVS( 299) = W( 43 )
  JVS( 300) = W( 47 )
  JVS( 301) = W( 52 )
  JVS( 302) = W( 63 )
  JVS( 303) = W( 73 )
  JVS( 304) = W( 87 )
  JVS( 305) = W( 90 )
  JVS( 306) = W( 95 )
  JVS( 307) = W( 96 )
  JVS( 308) = W( 97 )
  IF ( ABS(  JVS( 309 )) < TINY(a) ) THEN
         IER = 74                                      
         RETURN
  END IF
   W( 74 ) = JVS( 309 )
   W( 86 ) = JVS( 310 )
   W( 88 ) = JVS( 311 )
   W( 93 ) = JVS( 312 )
   W( 96 ) = JVS( 313 )
   W( 97 ) = JVS( 314 )
  JVS( 309) = W( 74 )
  JVS( 310) = W( 86 )
  JVS( 311) = W( 88 )
  JVS( 312) = W( 93 )
  JVS( 313) = W( 96 )
  JVS( 314) = W( 97 )
  IF ( ABS(  JVS( 315 )) < TINY(a) ) THEN
         IER = 75                                      
         RETURN
  END IF
   W( 75 ) = JVS( 315 )
   W( 86 ) = JVS( 316 )
   W( 93 ) = JVS( 317 )
   W( 96 ) = JVS( 318 )
   W( 97 ) = JVS( 319 )
  JVS( 315) = W( 75 )
  JVS( 316) = W( 86 )
  JVS( 317) = W( 93 )
  JVS( 318) = W( 96 )
  JVS( 319) = W( 97 )
  IF ( ABS(  JVS( 320 )) < TINY(a) ) THEN
         IER = 76                                      
         RETURN
  END IF
   W( 76 ) = JVS( 320 )
   W( 86 ) = JVS( 321 )
   W( 88 ) = JVS( 322 )
   W( 93 ) = JVS( 323 )
   W( 96 ) = JVS( 324 )
   W( 97 ) = JVS( 325 )
  JVS( 320) = W( 76 )
  JVS( 321) = W( 86 )
  JVS( 322) = W( 88 )
  JVS( 323) = W( 93 )
  JVS( 324) = W( 96 )
  JVS( 325) = W( 97 )
  IF ( ABS(  JVS( 327 )) < TINY(a) ) THEN
         IER = 77                                      
         RETURN
  END IF
   W( 76 ) = JVS( 326 )
   W( 77 ) = JVS( 327 )
   W( 86 ) = JVS( 328 )
   W( 88 ) = JVS( 329 )
   W( 93 ) = JVS( 330 )
   W( 96 ) = JVS( 331 )
   W( 97 ) = JVS( 332 )
  a = -W( 76 ) / JVS(          320  )
  W( 76 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 93 ) = W( 93 ) + a*JVS( 323 )
  W( 96 ) = W( 96 ) + a*JVS( 324 )
  W( 97 ) = W( 97 ) + a*JVS( 325 )
  JVS( 326) = W( 76 )
  JVS( 327) = W( 77 )
  JVS( 328) = W( 86 )
  JVS( 329) = W( 88 )
  JVS( 330) = W( 93 )
  JVS( 331) = W( 96 )
  JVS( 332) = W( 97 )
  IF ( ABS(  JVS( 337 )) < TINY(a) ) THEN
         IER = 78                                      
         RETURN
  END IF
   W( 57 ) = JVS( 333 )
   W( 72 ) = JVS( 334 )
   W( 75 ) = JVS( 335 )
   W( 77 ) = JVS( 336 )
   W( 78 ) = JVS( 337 )
   W( 79 ) = JVS( 338 )
   W( 80 ) = JVS( 339 )
   W( 86 ) = JVS( 340 )
   W( 87 ) = JVS( 341 )
   W( 88 ) = JVS( 342 )
   W( 89 ) = JVS( 343 )
   W( 90 ) = JVS( 344 )
   W( 93 ) = JVS( 345 )
   W( 95 ) = JVS( 346 )
   W( 96 ) = JVS( 347 )
   W( 97 ) = JVS( 348 )
  a = -W( 57 ) / JVS(          188  )
  W( 57 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 189 )
  W( 96 ) = W( 96 ) + a*JVS( 190 )
  a = -W( 72 ) / JVS(          294  )
  W( 72 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 295 )
  W( 87 ) = W( 87 ) + a*JVS( 296 )
  W( 88 ) = W( 88 ) + a*JVS( 297 )
  W( 96 ) = W( 96 ) + a*JVS( 298 )
  a = -W( 75 ) / JVS(          315  )
  W( 75 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 316 )
  W( 93 ) = W( 93 ) + a*JVS( 317 )
  W( 96 ) = W( 96 ) + a*JVS( 318 )
  W( 97 ) = W( 97 ) + a*JVS( 319 )
  a = -W( 77 ) / JVS(          327  )
  W( 77 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 328 )
  W( 88 ) = W( 88 ) + a*JVS( 329 )
  W( 93 ) = W( 93 ) + a*JVS( 330 )
  W( 96 ) = W( 96 ) + a*JVS( 331 )
  W( 97 ) = W( 97 ) + a*JVS( 332 )
  JVS( 333) = W( 57 )
  JVS( 334) = W( 72 )
  JVS( 335) = W( 75 )
  JVS( 336) = W( 77 )
  JVS( 337) = W( 78 )
  JVS( 338) = W( 79 )
  JVS( 339) = W( 80 )
  JVS( 340) = W( 86 )
  JVS( 341) = W( 87 )
  JVS( 342) = W( 88 )
  JVS( 343) = W( 89 )
  JVS( 344) = W( 90 )
  JVS( 345) = W( 93 )
  JVS( 346) = W( 95 )
  JVS( 347) = W( 96 )
  JVS( 348) = W( 97 )
  IF ( ABS(  JVS( 354 )) < TINY(a) ) THEN
         IER = 79                                      
         RETURN
  END IF
   W( 47 ) = JVS( 349 )
   W( 72 ) = JVS( 350 )
   W( 75 ) = JVS( 351 )
   W( 76 ) = JVS( 352 )
   W( 77 ) = JVS( 353 )
   W( 79 ) = JVS( 354 )
   W( 80 ) = JVS( 355 )
   W( 81 ) = JVS( 356 )
   W( 84 ) = JVS( 357 )
   W( 86 ) = JVS( 358 )
   W( 87 ) = JVS( 359 )
   W( 88 ) = JVS( 360 )
   W( 93 ) = JVS( 361 )
   W( 96 ) = JVS( 362 )
   W( 97 ) = JVS( 363 )
  a = -W( 47 ) / JVS(          144  )
  W( 47 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 145 )
  a = -W( 72 ) / JVS(          294  )
  W( 72 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 295 )
  W( 87 ) = W( 87 ) + a*JVS( 296 )
  W( 88 ) = W( 88 ) + a*JVS( 297 )
  W( 96 ) = W( 96 ) + a*JVS( 298 )
  a = -W( 75 ) / JVS(          315  )
  W( 75 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 316 )
  W( 93 ) = W( 93 ) + a*JVS( 317 )
  W( 96 ) = W( 96 ) + a*JVS( 318 )
  W( 97 ) = W( 97 ) + a*JVS( 319 )
  a = -W( 76 ) / JVS(          320  )
  W( 76 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 93 ) = W( 93 ) + a*JVS( 323 )
  W( 96 ) = W( 96 ) + a*JVS( 324 )
  W( 97 ) = W( 97 ) + a*JVS( 325 )
  a = -W( 77 ) / JVS(          327  )
  W( 77 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 328 )
  W( 88 ) = W( 88 ) + a*JVS( 329 )
  W( 93 ) = W( 93 ) + a*JVS( 330 )
  W( 96 ) = W( 96 ) + a*JVS( 331 )
  W( 97 ) = W( 97 ) + a*JVS( 332 )
  JVS( 349) = W( 47 )
  JVS( 350) = W( 72 )
  JVS( 351) = W( 75 )
  JVS( 352) = W( 76 )
  JVS( 353) = W( 77 )
  JVS( 354) = W( 79 )
  JVS( 355) = W( 80 )
  JVS( 356) = W( 81 )
  JVS( 357) = W( 84 )
  JVS( 358) = W( 86 )
  JVS( 359) = W( 87 )
  JVS( 360) = W( 88 )
  JVS( 361) = W( 93 )
  JVS( 362) = W( 96 )
  JVS( 363) = W( 97 )
  IF ( ABS(  JVS( 364 )) < TINY(a) ) THEN
         IER = 80                                      
         RETURN
  END IF
   W( 80 ) = JVS( 364 )
   W( 86 ) = JVS( 365 )
   W( 87 ) = JVS( 366 )
   W( 88 ) = JVS( 367 )
   W( 93 ) = JVS( 368 )
   W( 96 ) = JVS( 369 )
   W( 97 ) = JVS( 370 )
  JVS( 364) = W( 80 )
  JVS( 365) = W( 86 )
  JVS( 366) = W( 87 )
  JVS( 367) = W( 88 )
  JVS( 368) = W( 93 )
  JVS( 369) = W( 96 )
  JVS( 370) = W( 97 )
  IF ( ABS(  JVS( 372 )) < TINY(a) ) THEN
         IER = 81                                      
         RETURN
  END IF
   W( 80 ) = JVS( 371 )
   W( 81 ) = JVS( 372 )
   W( 86 ) = JVS( 373 )
   W( 87 ) = JVS( 374 )
   W( 88 ) = JVS( 375 )
   W( 93 ) = JVS( 376 )
   W( 96 ) = JVS( 377 )
   W( 97 ) = JVS( 378 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  JVS( 371) = W( 80 )
  JVS( 372) = W( 81 )
  JVS( 373) = W( 86 )
  JVS( 374) = W( 87 )
  JVS( 375) = W( 88 )
  JVS( 376) = W( 93 )
  JVS( 377) = W( 96 )
  JVS( 378) = W( 97 )
  IF ( ABS(  JVS( 392 )) < TINY(a) ) THEN
         IER = 82                                      
         RETURN
  END IF
   W( 59 ) = JVS( 379 )
   W( 65 ) = JVS( 380 )
   W( 71 ) = JVS( 381 )
   W( 72 ) = JVS( 382 )
   W( 73 ) = JVS( 383 )
   W( 74 ) = JVS( 384 )
   W( 75 ) = JVS( 385 )
   W( 76 ) = JVS( 386 )
   W( 77 ) = JVS( 387 )
   W( 78 ) = JVS( 388 )
   W( 79 ) = JVS( 389 )
   W( 80 ) = JVS( 390 )
   W( 81 ) = JVS( 391 )
   W( 82 ) = JVS( 392 )
   W( 84 ) = JVS( 393 )
   W( 86 ) = JVS( 394 )
   W( 87 ) = JVS( 395 )
   W( 88 ) = JVS( 396 )
   W( 89 ) = JVS( 397 )
   W( 90 ) = JVS( 398 )
   W( 93 ) = JVS( 399 )
   W( 95 ) = JVS( 400 )
   W( 96 ) = JVS( 401 )
   W( 97 ) = JVS( 402 )
  a = -W( 59 ) / JVS(          195  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 196 )
  W( 96 ) = W( 96 ) + a*JVS( 197 )
  a = -W( 65 ) / JVS(          233  )
  W( 65 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 234 )
  W( 89 ) = W( 89 ) + a*JVS( 235 )
  W( 90 ) = W( 90 ) + a*JVS( 236 )
  W( 96 ) = W( 96 ) + a*JVS( 237 )
  a = -W( 71 ) / JVS(          289  )
  W( 71 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 290 )
  W( 86 ) = W( 86 ) + a*JVS( 291 )
  W( 95 ) = W( 95 ) + a*JVS( 292 )
  W( 96 ) = W( 96 ) + a*JVS( 293 )
  a = -W( 72 ) / JVS(          294  )
  W( 72 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 295 )
  W( 87 ) = W( 87 ) + a*JVS( 296 )
  W( 88 ) = W( 88 ) + a*JVS( 297 )
  W( 96 ) = W( 96 ) + a*JVS( 298 )
  a = -W( 73 ) / JVS(          303  )
  W( 73 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 304 )
  W( 90 ) = W( 90 ) + a*JVS( 305 )
  W( 95 ) = W( 95 ) + a*JVS( 306 )
  W( 96 ) = W( 96 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  a = -W( 74 ) / JVS(          309  )
  W( 74 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 310 )
  W( 88 ) = W( 88 ) + a*JVS( 311 )
  W( 93 ) = W( 93 ) + a*JVS( 312 )
  W( 96 ) = W( 96 ) + a*JVS( 313 )
  W( 97 ) = W( 97 ) + a*JVS( 314 )
  a = -W( 75 ) / JVS(          315  )
  W( 75 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 316 )
  W( 93 ) = W( 93 ) + a*JVS( 317 )
  W( 96 ) = W( 96 ) + a*JVS( 318 )
  W( 97 ) = W( 97 ) + a*JVS( 319 )
  a = -W( 76 ) / JVS(          320  )
  W( 76 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 93 ) = W( 93 ) + a*JVS( 323 )
  W( 96 ) = W( 96 ) + a*JVS( 324 )
  W( 97 ) = W( 97 ) + a*JVS( 325 )
  a = -W( 77 ) / JVS(          327  )
  W( 77 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 328 )
  W( 88 ) = W( 88 ) + a*JVS( 329 )
  W( 93 ) = W( 93 ) + a*JVS( 330 )
  W( 96 ) = W( 96 ) + a*JVS( 331 )
  W( 97 ) = W( 97 ) + a*JVS( 332 )
  a = -W( 78 ) / JVS(          337  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 86 ) = W( 86 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 88 ) = W( 88 ) + a*JVS( 342 )
  W( 89 ) = W( 89 ) + a*JVS( 343 )
  W( 90 ) = W( 90 ) + a*JVS( 344 )
  W( 93 ) = W( 93 ) + a*JVS( 345 )
  W( 95 ) = W( 95 ) + a*JVS( 346 )
  W( 96 ) = W( 96 ) + a*JVS( 347 )
  W( 97 ) = W( 97 ) + a*JVS( 348 )
  a = -W( 79 ) / JVS(          354  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 355 )
  W( 81 ) = W( 81 ) + a*JVS( 356 )
  W( 84 ) = W( 84 ) + a*JVS( 357 )
  W( 86 ) = W( 86 ) + a*JVS( 358 )
  W( 87 ) = W( 87 ) + a*JVS( 359 )
  W( 88 ) = W( 88 ) + a*JVS( 360 )
  W( 93 ) = W( 93 ) + a*JVS( 361 )
  W( 96 ) = W( 96 ) + a*JVS( 362 )
  W( 97 ) = W( 97 ) + a*JVS( 363 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  a = -W( 81 ) / JVS(          372  )
  W( 81 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 373 )
  W( 87 ) = W( 87 ) + a*JVS( 374 )
  W( 88 ) = W( 88 ) + a*JVS( 375 )
  W( 93 ) = W( 93 ) + a*JVS( 376 )
  W( 96 ) = W( 96 ) + a*JVS( 377 )
  W( 97 ) = W( 97 ) + a*JVS( 378 )
  JVS( 379) = W( 59 )
  JVS( 380) = W( 65 )
  JVS( 381) = W( 71 )
  JVS( 382) = W( 72 )
  JVS( 383) = W( 73 )
  JVS( 384) = W( 74 )
  JVS( 385) = W( 75 )
  JVS( 386) = W( 76 )
  JVS( 387) = W( 77 )
  JVS( 388) = W( 78 )
  JVS( 389) = W( 79 )
  JVS( 390) = W( 80 )
  JVS( 391) = W( 81 )
  JVS( 392) = W( 82 )
  JVS( 393) = W( 84 )
  JVS( 394) = W( 86 )
  JVS( 395) = W( 87 )
  JVS( 396) = W( 88 )
  JVS( 397) = W( 89 )
  JVS( 398) = W( 90 )
  JVS( 399) = W( 93 )
  JVS( 400) = W( 95 )
  JVS( 401) = W( 96 )
  JVS( 402) = W( 97 )
  IF ( ABS(  JVS( 414 )) < TINY(a) ) THEN
         IER = 83                                      
         RETURN
  END IF
   W( 56 ) = JVS( 403 )
   W( 57 ) = JVS( 404 )
   W( 59 ) = JVS( 405 )
   W( 65 ) = JVS( 406 )
   W( 72 ) = JVS( 407 )
   W( 76 ) = JVS( 408 )
   W( 77 ) = JVS( 409 )
   W( 78 ) = JVS( 410 )
   W( 79 ) = JVS( 411 )
   W( 80 ) = JVS( 412 )
   W( 81 ) = JVS( 413 )
   W( 83 ) = JVS( 414 )
   W( 84 ) = JVS( 415 )
   W( 86 ) = JVS( 416 )
   W( 87 ) = JVS( 417 )
   W( 88 ) = JVS( 418 )
   W( 89 ) = JVS( 419 )
   W( 90 ) = JVS( 420 )
   W( 91 ) = JVS( 421 )
   W( 92 ) = JVS( 422 )
   W( 93 ) = JVS( 423 )
   W( 94 ) = JVS( 424 )
   W( 95 ) = JVS( 425 )
   W( 96 ) = JVS( 426 )
   W( 97 ) = JVS( 427 )
  a = -W( 56 ) / JVS(          184  )
  W( 56 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 185 )
  W( 94 ) = W( 94 ) + a*JVS( 186 )
  W( 96 ) = W( 96 ) + a*JVS( 187 )
  a = -W( 57 ) / JVS(          188  )
  W( 57 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 189 )
  W( 96 ) = W( 96 ) + a*JVS( 190 )
  a = -W( 59 ) / JVS(          195  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 196 )
  W( 96 ) = W( 96 ) + a*JVS( 197 )
  a = -W( 65 ) / JVS(          233  )
  W( 65 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 234 )
  W( 89 ) = W( 89 ) + a*JVS( 235 )
  W( 90 ) = W( 90 ) + a*JVS( 236 )
  W( 96 ) = W( 96 ) + a*JVS( 237 )
  a = -W( 72 ) / JVS(          294  )
  W( 72 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 295 )
  W( 87 ) = W( 87 ) + a*JVS( 296 )
  W( 88 ) = W( 88 ) + a*JVS( 297 )
  W( 96 ) = W( 96 ) + a*JVS( 298 )
  a = -W( 76 ) / JVS(          320  )
  W( 76 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 93 ) = W( 93 ) + a*JVS( 323 )
  W( 96 ) = W( 96 ) + a*JVS( 324 )
  W( 97 ) = W( 97 ) + a*JVS( 325 )
  a = -W( 77 ) / JVS(          327  )
  W( 77 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 328 )
  W( 88 ) = W( 88 ) + a*JVS( 329 )
  W( 93 ) = W( 93 ) + a*JVS( 330 )
  W( 96 ) = W( 96 ) + a*JVS( 331 )
  W( 97 ) = W( 97 ) + a*JVS( 332 )
  a = -W( 78 ) / JVS(          337  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 86 ) = W( 86 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 88 ) = W( 88 ) + a*JVS( 342 )
  W( 89 ) = W( 89 ) + a*JVS( 343 )
  W( 90 ) = W( 90 ) + a*JVS( 344 )
  W( 93 ) = W( 93 ) + a*JVS( 345 )
  W( 95 ) = W( 95 ) + a*JVS( 346 )
  W( 96 ) = W( 96 ) + a*JVS( 347 )
  W( 97 ) = W( 97 ) + a*JVS( 348 )
  a = -W( 79 ) / JVS(          354  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 355 )
  W( 81 ) = W( 81 ) + a*JVS( 356 )
  W( 84 ) = W( 84 ) + a*JVS( 357 )
  W( 86 ) = W( 86 ) + a*JVS( 358 )
  W( 87 ) = W( 87 ) + a*JVS( 359 )
  W( 88 ) = W( 88 ) + a*JVS( 360 )
  W( 93 ) = W( 93 ) + a*JVS( 361 )
  W( 96 ) = W( 96 ) + a*JVS( 362 )
  W( 97 ) = W( 97 ) + a*JVS( 363 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  a = -W( 81 ) / JVS(          372  )
  W( 81 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 373 )
  W( 87 ) = W( 87 ) + a*JVS( 374 )
  W( 88 ) = W( 88 ) + a*JVS( 375 )
  W( 93 ) = W( 93 ) + a*JVS( 376 )
  W( 96 ) = W( 96 ) + a*JVS( 377 )
  W( 97 ) = W( 97 ) + a*JVS( 378 )
  JVS( 403) = W( 56 )
  JVS( 404) = W( 57 )
  JVS( 405) = W( 59 )
  JVS( 406) = W( 65 )
  JVS( 407) = W( 72 )
  JVS( 408) = W( 76 )
  JVS( 409) = W( 77 )
  JVS( 410) = W( 78 )
  JVS( 411) = W( 79 )
  JVS( 412) = W( 80 )
  JVS( 413) = W( 81 )
  JVS( 414) = W( 83 )
  JVS( 415) = W( 84 )
  JVS( 416) = W( 86 )
  JVS( 417) = W( 87 )
  JVS( 418) = W( 88 )
  JVS( 419) = W( 89 )
  JVS( 420) = W( 90 )
  JVS( 421) = W( 91 )
  JVS( 422) = W( 92 )
  JVS( 423) = W( 93 )
  JVS( 424) = W( 94 )
  JVS( 425) = W( 95 )
  JVS( 426) = W( 96 )
  JVS( 427) = W( 97 )
  IF ( ABS(  JVS( 437 )) < TINY(a) ) THEN
         IER = 84                                      
         RETURN
  END IF
   W( 52 ) = JVS( 428 )
   W( 63 ) = JVS( 429 )
   W( 72 ) = JVS( 430 )
   W( 73 ) = JVS( 431 )
   W( 75 ) = JVS( 432 )
   W( 78 ) = JVS( 433 )
   W( 79 ) = JVS( 434 )
   W( 80 ) = JVS( 435 )
   W( 81 ) = JVS( 436 )
   W( 84 ) = JVS( 437 )
   W( 86 ) = JVS( 438 )
   W( 87 ) = JVS( 439 )
   W( 88 ) = JVS( 440 )
   W( 89 ) = JVS( 441 )
   W( 90 ) = JVS( 442 )
   W( 93 ) = JVS( 443 )
   W( 95 ) = JVS( 444 )
   W( 96 ) = JVS( 445 )
   W( 97 ) = JVS( 446 )
  a = -W( 52 ) / JVS(          165  )
  W( 52 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 166 )
  W( 96 ) = W( 96 ) + a*JVS( 167 )
  a = -W( 63 ) / JVS(          222  )
  W( 63 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 223 )
  W( 87 ) = W( 87 ) + a*JVS( 224 )
  W( 90 ) = W( 90 ) + a*JVS( 225 )
  W( 96 ) = W( 96 ) + a*JVS( 226 )
  W( 97 ) = W( 97 ) + a*JVS( 227 )
  a = -W( 72 ) / JVS(          294  )
  W( 72 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 295 )
  W( 87 ) = W( 87 ) + a*JVS( 296 )
  W( 88 ) = W( 88 ) + a*JVS( 297 )
  W( 96 ) = W( 96 ) + a*JVS( 298 )
  a = -W( 73 ) / JVS(          303  )
  W( 73 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 304 )
  W( 90 ) = W( 90 ) + a*JVS( 305 )
  W( 95 ) = W( 95 ) + a*JVS( 306 )
  W( 96 ) = W( 96 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  a = -W( 75 ) / JVS(          315  )
  W( 75 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 316 )
  W( 93 ) = W( 93 ) + a*JVS( 317 )
  W( 96 ) = W( 96 ) + a*JVS( 318 )
  W( 97 ) = W( 97 ) + a*JVS( 319 )
  a = -W( 78 ) / JVS(          337  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 86 ) = W( 86 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 88 ) = W( 88 ) + a*JVS( 342 )
  W( 89 ) = W( 89 ) + a*JVS( 343 )
  W( 90 ) = W( 90 ) + a*JVS( 344 )
  W( 93 ) = W( 93 ) + a*JVS( 345 )
  W( 95 ) = W( 95 ) + a*JVS( 346 )
  W( 96 ) = W( 96 ) + a*JVS( 347 )
  W( 97 ) = W( 97 ) + a*JVS( 348 )
  a = -W( 79 ) / JVS(          354  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 355 )
  W( 81 ) = W( 81 ) + a*JVS( 356 )
  W( 84 ) = W( 84 ) + a*JVS( 357 )
  W( 86 ) = W( 86 ) + a*JVS( 358 )
  W( 87 ) = W( 87 ) + a*JVS( 359 )
  W( 88 ) = W( 88 ) + a*JVS( 360 )
  W( 93 ) = W( 93 ) + a*JVS( 361 )
  W( 96 ) = W( 96 ) + a*JVS( 362 )
  W( 97 ) = W( 97 ) + a*JVS( 363 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  a = -W( 81 ) / JVS(          372  )
  W( 81 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 373 )
  W( 87 ) = W( 87 ) + a*JVS( 374 )
  W( 88 ) = W( 88 ) + a*JVS( 375 )
  W( 93 ) = W( 93 ) + a*JVS( 376 )
  W( 96 ) = W( 96 ) + a*JVS( 377 )
  W( 97 ) = W( 97 ) + a*JVS( 378 )
  JVS( 428) = W( 52 )
  JVS( 429) = W( 63 )
  JVS( 430) = W( 72 )
  JVS( 431) = W( 73 )
  JVS( 432) = W( 75 )
  JVS( 433) = W( 78 )
  JVS( 434) = W( 79 )
  JVS( 435) = W( 80 )
  JVS( 436) = W( 81 )
  JVS( 437) = W( 84 )
  JVS( 438) = W( 86 )
  JVS( 439) = W( 87 )
  JVS( 440) = W( 88 )
  JVS( 441) = W( 89 )
  JVS( 442) = W( 90 )
  JVS( 443) = W( 93 )
  JVS( 444) = W( 95 )
  JVS( 445) = W( 96 )
  JVS( 446) = W( 97 )
  IF ( ABS(  JVS( 460 )) < TINY(a) ) THEN
         IER = 85                                      
         RETURN
  END IF
   W( 58 ) = JVS( 447 )
   W( 59 ) = JVS( 448 )
   W( 62 ) = JVS( 449 )
   W( 68 ) = JVS( 450 )
   W( 71 ) = JVS( 451 )
   W( 73 ) = JVS( 452 )
   W( 74 ) = JVS( 453 )
   W( 75 ) = JVS( 454 )
   W( 76 ) = JVS( 455 )
   W( 77 ) = JVS( 456 )
   W( 80 ) = JVS( 457 )
   W( 81 ) = JVS( 458 )
   W( 84 ) = JVS( 459 )
   W( 85 ) = JVS( 460 )
   W( 86 ) = JVS( 461 )
   W( 87 ) = JVS( 462 )
   W( 88 ) = JVS( 463 )
   W( 89 ) = JVS( 464 )
   W( 90 ) = JVS( 465 )
   W( 91 ) = JVS( 466 )
   W( 92 ) = JVS( 467 )
   W( 93 ) = JVS( 468 )
   W( 94 ) = JVS( 469 )
   W( 95 ) = JVS( 470 )
   W( 96 ) = JVS( 471 )
   W( 97 ) = JVS( 472 )
  a = -W( 58 ) / JVS(          191  )
  W( 58 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 192 )
  W( 92 ) = W( 92 ) + a*JVS( 193 )
  W( 96 ) = W( 96 ) + a*JVS( 194 )
  a = -W( 59 ) / JVS(          195  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 196 )
  W( 96 ) = W( 96 ) + a*JVS( 197 )
  a = -W( 62 ) / JVS(          218  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 219 )
  W( 90 ) = W( 90 ) + a*JVS( 220 )
  W( 95 ) = W( 95 ) + a*JVS( 221 )
  a = -W( 68 ) / JVS(          252  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 253 )
  W( 90 ) = W( 90 ) + a*JVS( 254 )
  W( 92 ) = W( 92 ) + a*JVS( 255 )
  W( 95 ) = W( 95 ) + a*JVS( 256 )
  W( 96 ) = W( 96 ) + a*JVS( 257 )
  a = -W( 71 ) / JVS(          289  )
  W( 71 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 290 )
  W( 86 ) = W( 86 ) + a*JVS( 291 )
  W( 95 ) = W( 95 ) + a*JVS( 292 )
  W( 96 ) = W( 96 ) + a*JVS( 293 )
  a = -W( 73 ) / JVS(          303  )
  W( 73 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 304 )
  W( 90 ) = W( 90 ) + a*JVS( 305 )
  W( 95 ) = W( 95 ) + a*JVS( 306 )
  W( 96 ) = W( 96 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  a = -W( 74 ) / JVS(          309  )
  W( 74 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 310 )
  W( 88 ) = W( 88 ) + a*JVS( 311 )
  W( 93 ) = W( 93 ) + a*JVS( 312 )
  W( 96 ) = W( 96 ) + a*JVS( 313 )
  W( 97 ) = W( 97 ) + a*JVS( 314 )
  a = -W( 75 ) / JVS(          315  )
  W( 75 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 316 )
  W( 93 ) = W( 93 ) + a*JVS( 317 )
  W( 96 ) = W( 96 ) + a*JVS( 318 )
  W( 97 ) = W( 97 ) + a*JVS( 319 )
  a = -W( 76 ) / JVS(          320  )
  W( 76 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 93 ) = W( 93 ) + a*JVS( 323 )
  W( 96 ) = W( 96 ) + a*JVS( 324 )
  W( 97 ) = W( 97 ) + a*JVS( 325 )
  a = -W( 77 ) / JVS(          327  )
  W( 77 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 328 )
  W( 88 ) = W( 88 ) + a*JVS( 329 )
  W( 93 ) = W( 93 ) + a*JVS( 330 )
  W( 96 ) = W( 96 ) + a*JVS( 331 )
  W( 97 ) = W( 97 ) + a*JVS( 332 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  a = -W( 81 ) / JVS(          372  )
  W( 81 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 373 )
  W( 87 ) = W( 87 ) + a*JVS( 374 )
  W( 88 ) = W( 88 ) + a*JVS( 375 )
  W( 93 ) = W( 93 ) + a*JVS( 376 )
  W( 96 ) = W( 96 ) + a*JVS( 377 )
  W( 97 ) = W( 97 ) + a*JVS( 378 )
  a = -W( 84 ) / JVS(          437  )
  W( 84 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 88 ) = W( 88 ) + a*JVS( 440 )
  W( 89 ) = W( 89 ) + a*JVS( 441 )
  W( 90 ) = W( 90 ) + a*JVS( 442 )
  W( 93 ) = W( 93 ) + a*JVS( 443 )
  W( 95 ) = W( 95 ) + a*JVS( 444 )
  W( 96 ) = W( 96 ) + a*JVS( 445 )
  W( 97 ) = W( 97 ) + a*JVS( 446 )
  JVS( 447) = W( 58 )
  JVS( 448) = W( 59 )
  JVS( 449) = W( 62 )
  JVS( 450) = W( 68 )
  JVS( 451) = W( 71 )
  JVS( 452) = W( 73 )
  JVS( 453) = W( 74 )
  JVS( 454) = W( 75 )
  JVS( 455) = W( 76 )
  JVS( 456) = W( 77 )
  JVS( 457) = W( 80 )
  JVS( 458) = W( 81 )
  JVS( 459) = W( 84 )
  JVS( 460) = W( 85 )
  JVS( 461) = W( 86 )
  JVS( 462) = W( 87 )
  JVS( 463) = W( 88 )
  JVS( 464) = W( 89 )
  JVS( 465) = W( 90 )
  JVS( 466) = W( 91 )
  JVS( 467) = W( 92 )
  JVS( 468) = W( 93 )
  JVS( 469) = W( 94 )
  JVS( 470) = W( 95 )
  JVS( 471) = W( 96 )
  JVS( 472) = W( 97 )
  IF ( ABS(  JVS( 483 )) < TINY(a) ) THEN
         IER = 86                                      
         RETURN
  END IF
   W( 49 ) = JVS( 473 )
   W( 60 ) = JVS( 474 )
   W( 71 ) = JVS( 475 )
   W( 73 ) = JVS( 476 )
   W( 74 ) = JVS( 477 )
   W( 75 ) = JVS( 478 )
   W( 76 ) = JVS( 479 )
   W( 77 ) = JVS( 480 )
   W( 80 ) = JVS( 481 )
   W( 81 ) = JVS( 482 )
   W( 86 ) = JVS( 483 )
   W( 87 ) = JVS( 484 )
   W( 88 ) = JVS( 485 )
   W( 90 ) = JVS( 486 )
   W( 91 ) = JVS( 487 )
   W( 93 ) = JVS( 488 )
   W( 94 ) = JVS( 489 )
   W( 95 ) = JVS( 490 )
   W( 96 ) = JVS( 491 )
   W( 97 ) = JVS( 492 )
  a = -W( 49 ) / JVS(          149  )
  W( 49 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 150 )
  W( 86 ) = W( 86 ) + a*JVS( 151 )
  W( 96 ) = W( 96 ) + a*JVS( 152 )
  a = -W( 60 ) / JVS(          199  )
  W( 60 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 200 )
  W( 90 ) = W( 90 ) + a*JVS( 201 )
  W( 93 ) = W( 93 ) + a*JVS( 202 )
  W( 96 ) = W( 96 ) + a*JVS( 203 )
  a = -W( 71 ) / JVS(          289  )
  W( 71 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 290 )
  W( 86 ) = W( 86 ) + a*JVS( 291 )
  W( 95 ) = W( 95 ) + a*JVS( 292 )
  W( 96 ) = W( 96 ) + a*JVS( 293 )
  a = -W( 73 ) / JVS(          303  )
  W( 73 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 304 )
  W( 90 ) = W( 90 ) + a*JVS( 305 )
  W( 95 ) = W( 95 ) + a*JVS( 306 )
  W( 96 ) = W( 96 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  a = -W( 74 ) / JVS(          309  )
  W( 74 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 310 )
  W( 88 ) = W( 88 ) + a*JVS( 311 )
  W( 93 ) = W( 93 ) + a*JVS( 312 )
  W( 96 ) = W( 96 ) + a*JVS( 313 )
  W( 97 ) = W( 97 ) + a*JVS( 314 )
  a = -W( 75 ) / JVS(          315  )
  W( 75 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 316 )
  W( 93 ) = W( 93 ) + a*JVS( 317 )
  W( 96 ) = W( 96 ) + a*JVS( 318 )
  W( 97 ) = W( 97 ) + a*JVS( 319 )
  a = -W( 76 ) / JVS(          320  )
  W( 76 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 93 ) = W( 93 ) + a*JVS( 323 )
  W( 96 ) = W( 96 ) + a*JVS( 324 )
  W( 97 ) = W( 97 ) + a*JVS( 325 )
  a = -W( 77 ) / JVS(          327  )
  W( 77 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 328 )
  W( 88 ) = W( 88 ) + a*JVS( 329 )
  W( 93 ) = W( 93 ) + a*JVS( 330 )
  W( 96 ) = W( 96 ) + a*JVS( 331 )
  W( 97 ) = W( 97 ) + a*JVS( 332 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  a = -W( 81 ) / JVS(          372  )
  W( 81 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 373 )
  W( 87 ) = W( 87 ) + a*JVS( 374 )
  W( 88 ) = W( 88 ) + a*JVS( 375 )
  W( 93 ) = W( 93 ) + a*JVS( 376 )
  W( 96 ) = W( 96 ) + a*JVS( 377 )
  W( 97 ) = W( 97 ) + a*JVS( 378 )
  JVS( 473) = W( 49 )
  JVS( 474) = W( 60 )
  JVS( 475) = W( 71 )
  JVS( 476) = W( 73 )
  JVS( 477) = W( 74 )
  JVS( 478) = W( 75 )
  JVS( 479) = W( 76 )
  JVS( 480) = W( 77 )
  JVS( 481) = W( 80 )
  JVS( 482) = W( 81 )
  JVS( 483) = W( 86 )
  JVS( 484) = W( 87 )
  JVS( 485) = W( 88 )
  JVS( 486) = W( 90 )
  JVS( 487) = W( 91 )
  JVS( 488) = W( 93 )
  JVS( 489) = W( 94 )
  JVS( 490) = W( 95 )
  JVS( 491) = W( 96 )
  JVS( 492) = W( 97 )
  IF ( ABS(  JVS( 517 )) < TINY(a) ) THEN
         IER = 87                                      
         RETURN
  END IF
   W( 41 ) = JVS( 493 )
   W( 46 ) = JVS( 494 )
   W( 50 ) = JVS( 495 )
   W( 52 ) = JVS( 496 )
   W( 53 ) = JVS( 497 )
   W( 56 ) = JVS( 498 )
   W( 62 ) = JVS( 499 )
   W( 63 ) = JVS( 500 )
   W( 64 ) = JVS( 501 )
   W( 69 ) = JVS( 502 )
   W( 72 ) = JVS( 503 )
   W( 73 ) = JVS( 504 )
   W( 74 ) = JVS( 505 )
   W( 75 ) = JVS( 506 )
   W( 76 ) = JVS( 507 )
   W( 77 ) = JVS( 508 )
   W( 79 ) = JVS( 509 )
   W( 80 ) = JVS( 510 )
   W( 81 ) = JVS( 511 )
   W( 82 ) = JVS( 512 )
   W( 83 ) = JVS( 513 )
   W( 84 ) = JVS( 514 )
   W( 85 ) = JVS( 515 )
   W( 86 ) = JVS( 516 )
   W( 87 ) = JVS( 517 )
   W( 88 ) = JVS( 518 )
   W( 89 ) = JVS( 519 )
   W( 90 ) = JVS( 520 )
   W( 91 ) = JVS( 521 )
   W( 92 ) = JVS( 522 )
   W( 93 ) = JVS( 523 )
   W( 94 ) = JVS( 524 )
   W( 95 ) = JVS( 525 )
   W( 96 ) = JVS( 526 )
   W( 97 ) = JVS( 527 )
  a = -W( 41 ) / JVS(          126  )
  W( 41 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 127 )
  W( 91 ) = W( 91 ) + a*JVS( 128 )
  a = -W( 46 ) / JVS(          141  )
  W( 46 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 142 )
  W( 97 ) = W( 97 ) + a*JVS( 143 )
  a = -W( 50 ) / JVS(          153  )
  W( 50 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 154 )
  W( 95 ) = W( 95 ) + a*JVS( 155 )
  W( 96 ) = W( 96 ) + a*JVS( 156 )
  a = -W( 52 ) / JVS(          165  )
  W( 52 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 166 )
  W( 96 ) = W( 96 ) + a*JVS( 167 )
  a = -W( 53 ) / JVS(          168  )
  W( 53 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 169 )
  W( 90 ) = W( 90 ) + a*JVS( 170 )
  W( 96 ) = W( 96 ) + a*JVS( 171 )
  a = -W( 56 ) / JVS(          184  )
  W( 56 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 185 )
  W( 94 ) = W( 94 ) + a*JVS( 186 )
  W( 96 ) = W( 96 ) + a*JVS( 187 )
  a = -W( 62 ) / JVS(          218  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 219 )
  W( 90 ) = W( 90 ) + a*JVS( 220 )
  W( 95 ) = W( 95 ) + a*JVS( 221 )
  a = -W( 63 ) / JVS(          222  )
  W( 63 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 223 )
  W( 87 ) = W( 87 ) + a*JVS( 224 )
  W( 90 ) = W( 90 ) + a*JVS( 225 )
  W( 96 ) = W( 96 ) + a*JVS( 226 )
  W( 97 ) = W( 97 ) + a*JVS( 227 )
  a = -W( 64 ) / JVS(          228  )
  W( 64 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 229 )
  W( 88 ) = W( 88 ) + a*JVS( 230 )
  W( 90 ) = W( 90 ) + a*JVS( 231 )
  W( 95 ) = W( 95 ) + a*JVS( 232 )
  a = -W( 69 ) / JVS(          259  )
  W( 69 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 260 )
  W( 81 ) = W( 81 ) + a*JVS( 261 )
  W( 82 ) = W( 82 ) + a*JVS( 262 )
  W( 83 ) = W( 83 ) + a*JVS( 263 )
  W( 84 ) = W( 84 ) + a*JVS( 264 )
  W( 85 ) = W( 85 ) + a*JVS( 265 )
  W( 87 ) = W( 87 ) + a*JVS( 266 )
  W( 90 ) = W( 90 ) + a*JVS( 267 )
  W( 96 ) = W( 96 ) + a*JVS( 268 )
  W( 97 ) = W( 97 ) + a*JVS( 269 )
  a = -W( 72 ) / JVS(          294  )
  W( 72 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 295 )
  W( 87 ) = W( 87 ) + a*JVS( 296 )
  W( 88 ) = W( 88 ) + a*JVS( 297 )
  W( 96 ) = W( 96 ) + a*JVS( 298 )
  a = -W( 73 ) / JVS(          303  )
  W( 73 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 304 )
  W( 90 ) = W( 90 ) + a*JVS( 305 )
  W( 95 ) = W( 95 ) + a*JVS( 306 )
  W( 96 ) = W( 96 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  a = -W( 74 ) / JVS(          309  )
  W( 74 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 310 )
  W( 88 ) = W( 88 ) + a*JVS( 311 )
  W( 93 ) = W( 93 ) + a*JVS( 312 )
  W( 96 ) = W( 96 ) + a*JVS( 313 )
  W( 97 ) = W( 97 ) + a*JVS( 314 )
  a = -W( 75 ) / JVS(          315  )
  W( 75 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 316 )
  W( 93 ) = W( 93 ) + a*JVS( 317 )
  W( 96 ) = W( 96 ) + a*JVS( 318 )
  W( 97 ) = W( 97 ) + a*JVS( 319 )
  a = -W( 76 ) / JVS(          320  )
  W( 76 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 93 ) = W( 93 ) + a*JVS( 323 )
  W( 96 ) = W( 96 ) + a*JVS( 324 )
  W( 97 ) = W( 97 ) + a*JVS( 325 )
  a = -W( 77 ) / JVS(          327  )
  W( 77 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 328 )
  W( 88 ) = W( 88 ) + a*JVS( 329 )
  W( 93 ) = W( 93 ) + a*JVS( 330 )
  W( 96 ) = W( 96 ) + a*JVS( 331 )
  W( 97 ) = W( 97 ) + a*JVS( 332 )
  a = -W( 79 ) / JVS(          354  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 355 )
  W( 81 ) = W( 81 ) + a*JVS( 356 )
  W( 84 ) = W( 84 ) + a*JVS( 357 )
  W( 86 ) = W( 86 ) + a*JVS( 358 )
  W( 87 ) = W( 87 ) + a*JVS( 359 )
  W( 88 ) = W( 88 ) + a*JVS( 360 )
  W( 93 ) = W( 93 ) + a*JVS( 361 )
  W( 96 ) = W( 96 ) + a*JVS( 362 )
  W( 97 ) = W( 97 ) + a*JVS( 363 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  a = -W( 81 ) / JVS(          372  )
  W( 81 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 373 )
  W( 87 ) = W( 87 ) + a*JVS( 374 )
  W( 88 ) = W( 88 ) + a*JVS( 375 )
  W( 93 ) = W( 93 ) + a*JVS( 376 )
  W( 96 ) = W( 96 ) + a*JVS( 377 )
  W( 97 ) = W( 97 ) + a*JVS( 378 )
  a = -W( 82 ) / JVS(          392  )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 393 )
  W( 86 ) = W( 86 ) + a*JVS( 394 )
  W( 87 ) = W( 87 ) + a*JVS( 395 )
  W( 88 ) = W( 88 ) + a*JVS( 396 )
  W( 89 ) = W( 89 ) + a*JVS( 397 )
  W( 90 ) = W( 90 ) + a*JVS( 398 )
  W( 93 ) = W( 93 ) + a*JVS( 399 )
  W( 95 ) = W( 95 ) + a*JVS( 400 )
  W( 96 ) = W( 96 ) + a*JVS( 401 )
  W( 97 ) = W( 97 ) + a*JVS( 402 )
  a = -W( 83 ) / JVS(          414  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 415 )
  W( 86 ) = W( 86 ) + a*JVS( 416 )
  W( 87 ) = W( 87 ) + a*JVS( 417 )
  W( 88 ) = W( 88 ) + a*JVS( 418 )
  W( 89 ) = W( 89 ) + a*JVS( 419 )
  W( 90 ) = W( 90 ) + a*JVS( 420 )
  W( 91 ) = W( 91 ) + a*JVS( 421 )
  W( 92 ) = W( 92 ) + a*JVS( 422 )
  W( 93 ) = W( 93 ) + a*JVS( 423 )
  W( 94 ) = W( 94 ) + a*JVS( 424 )
  W( 95 ) = W( 95 ) + a*JVS( 425 )
  W( 96 ) = W( 96 ) + a*JVS( 426 )
  W( 97 ) = W( 97 ) + a*JVS( 427 )
  a = -W( 84 ) / JVS(          437  )
  W( 84 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 88 ) = W( 88 ) + a*JVS( 440 )
  W( 89 ) = W( 89 ) + a*JVS( 441 )
  W( 90 ) = W( 90 ) + a*JVS( 442 )
  W( 93 ) = W( 93 ) + a*JVS( 443 )
  W( 95 ) = W( 95 ) + a*JVS( 444 )
  W( 96 ) = W( 96 ) + a*JVS( 445 )
  W( 97 ) = W( 97 ) + a*JVS( 446 )
  a = -W( 85 ) / JVS(          460  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 461 )
  W( 87 ) = W( 87 ) + a*JVS( 462 )
  W( 88 ) = W( 88 ) + a*JVS( 463 )
  W( 89 ) = W( 89 ) + a*JVS( 464 )
  W( 90 ) = W( 90 ) + a*JVS( 465 )
  W( 91 ) = W( 91 ) + a*JVS( 466 )
  W( 92 ) = W( 92 ) + a*JVS( 467 )
  W( 93 ) = W( 93 ) + a*JVS( 468 )
  W( 94 ) = W( 94 ) + a*JVS( 469 )
  W( 95 ) = W( 95 ) + a*JVS( 470 )
  W( 96 ) = W( 96 ) + a*JVS( 471 )
  W( 97 ) = W( 97 ) + a*JVS( 472 )
  a = -W( 86 ) / JVS(          483  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 90 ) = W( 90 ) + a*JVS( 486 )
  W( 91 ) = W( 91 ) + a*JVS( 487 )
  W( 93 ) = W( 93 ) + a*JVS( 488 )
  W( 94 ) = W( 94 ) + a*JVS( 489 )
  W( 95 ) = W( 95 ) + a*JVS( 490 )
  W( 96 ) = W( 96 ) + a*JVS( 491 )
  W( 97 ) = W( 97 ) + a*JVS( 492 )
  JVS( 493) = W( 41 )
  JVS( 494) = W( 46 )
  JVS( 495) = W( 50 )
  JVS( 496) = W( 52 )
  JVS( 497) = W( 53 )
  JVS( 498) = W( 56 )
  JVS( 499) = W( 62 )
  JVS( 500) = W( 63 )
  JVS( 501) = W( 64 )
  JVS( 502) = W( 69 )
  JVS( 503) = W( 72 )
  JVS( 504) = W( 73 )
  JVS( 505) = W( 74 )
  JVS( 506) = W( 75 )
  JVS( 507) = W( 76 )
  JVS( 508) = W( 77 )
  JVS( 509) = W( 79 )
  JVS( 510) = W( 80 )
  JVS( 511) = W( 81 )
  JVS( 512) = W( 82 )
  JVS( 513) = W( 83 )
  JVS( 514) = W( 84 )
  JVS( 515) = W( 85 )
  JVS( 516) = W( 86 )
  JVS( 517) = W( 87 )
  JVS( 518) = W( 88 )
  JVS( 519) = W( 89 )
  JVS( 520) = W( 90 )
  JVS( 521) = W( 91 )
  JVS( 522) = W( 92 )
  JVS( 523) = W( 93 )
  JVS( 524) = W( 94 )
  JVS( 525) = W( 95 )
  JVS( 526) = W( 96 )
  JVS( 527) = W( 97 )
  IF ( ABS(  JVS( 549 )) < TINY(a) ) THEN
         IER = 88                                      
         RETURN
  END IF
   W( 39 ) = JVS( 528 )
   W( 42 ) = JVS( 529 )
   W( 48 ) = JVS( 530 )
   W( 57 ) = JVS( 531 )
   W( 58 ) = JVS( 532 )
   W( 59 ) = JVS( 533 )
   W( 61 ) = JVS( 534 )
   W( 64 ) = JVS( 535 )
   W( 67 ) = JVS( 536 )
   W( 74 ) = JVS( 537 )
   W( 76 ) = JVS( 538 )
   W( 77 ) = JVS( 539 )
   W( 79 ) = JVS( 540 )
   W( 80 ) = JVS( 541 )
   W( 81 ) = JVS( 542 )
   W( 82 ) = JVS( 543 )
   W( 83 ) = JVS( 544 )
   W( 84 ) = JVS( 545 )
   W( 85 ) = JVS( 546 )
   W( 86 ) = JVS( 547 )
   W( 87 ) = JVS( 548 )
   W( 88 ) = JVS( 549 )
   W( 89 ) = JVS( 550 )
   W( 90 ) = JVS( 551 )
   W( 91 ) = JVS( 552 )
   W( 92 ) = JVS( 553 )
   W( 93 ) = JVS( 554 )
   W( 94 ) = JVS( 555 )
   W( 95 ) = JVS( 556 )
   W( 96 ) = JVS( 557 )
   W( 97 ) = JVS( 558 )
  a = -W( 39 ) / JVS(          122  )
  W( 39 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 123 )
  a = -W( 42 ) / JVS(          129  )
  W( 42 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 130 )
  W( 90 ) = W( 90 ) + a*JVS( 131 )
  a = -W( 48 ) / JVS(          146  )
  W( 48 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 147 )
  W( 96 ) = W( 96 ) + a*JVS( 148 )
  a = -W( 57 ) / JVS(          188  )
  W( 57 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 189 )
  W( 96 ) = W( 96 ) + a*JVS( 190 )
  a = -W( 58 ) / JVS(          191  )
  W( 58 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 192 )
  W( 92 ) = W( 92 ) + a*JVS( 193 )
  W( 96 ) = W( 96 ) + a*JVS( 194 )
  a = -W( 59 ) / JVS(          195  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 196 )
  W( 96 ) = W( 96 ) + a*JVS( 197 )
  a = -W( 61 ) / JVS(          208  )
  W( 61 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 209 )
  W( 79 ) = W( 79 ) + a*JVS( 210 )
  W( 80 ) = W( 80 ) + a*JVS( 211 )
  W( 82 ) = W( 82 ) + a*JVS( 212 )
  W( 83 ) = W( 83 ) + a*JVS( 213 )
  W( 85 ) = W( 85 ) + a*JVS( 214 )
  W( 88 ) = W( 88 ) + a*JVS( 215 )
  W( 92 ) = W( 92 ) + a*JVS( 216 )
  W( 96 ) = W( 96 ) + a*JVS( 217 )
  a = -W( 64 ) / JVS(          228  )
  W( 64 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 229 )
  W( 88 ) = W( 88 ) + a*JVS( 230 )
  W( 90 ) = W( 90 ) + a*JVS( 231 )
  W( 95 ) = W( 95 ) + a*JVS( 232 )
  a = -W( 67 ) / JVS(          244  )
  W( 67 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 245 )
  W( 76 ) = W( 76 ) + a*JVS( 246 )
  W( 77 ) = W( 77 ) + a*JVS( 247 )
  W( 80 ) = W( 80 ) + a*JVS( 248 )
  W( 88 ) = W( 88 ) + a*JVS( 249 )
  W( 96 ) = W( 96 ) + a*JVS( 250 )
  a = -W( 74 ) / JVS(          309  )
  W( 74 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 310 )
  W( 88 ) = W( 88 ) + a*JVS( 311 )
  W( 93 ) = W( 93 ) + a*JVS( 312 )
  W( 96 ) = W( 96 ) + a*JVS( 313 )
  W( 97 ) = W( 97 ) + a*JVS( 314 )
  a = -W( 76 ) / JVS(          320  )
  W( 76 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 93 ) = W( 93 ) + a*JVS( 323 )
  W( 96 ) = W( 96 ) + a*JVS( 324 )
  W( 97 ) = W( 97 ) + a*JVS( 325 )
  a = -W( 77 ) / JVS(          327  )
  W( 77 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 328 )
  W( 88 ) = W( 88 ) + a*JVS( 329 )
  W( 93 ) = W( 93 ) + a*JVS( 330 )
  W( 96 ) = W( 96 ) + a*JVS( 331 )
  W( 97 ) = W( 97 ) + a*JVS( 332 )
  a = -W( 79 ) / JVS(          354  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 355 )
  W( 81 ) = W( 81 ) + a*JVS( 356 )
  W( 84 ) = W( 84 ) + a*JVS( 357 )
  W( 86 ) = W( 86 ) + a*JVS( 358 )
  W( 87 ) = W( 87 ) + a*JVS( 359 )
  W( 88 ) = W( 88 ) + a*JVS( 360 )
  W( 93 ) = W( 93 ) + a*JVS( 361 )
  W( 96 ) = W( 96 ) + a*JVS( 362 )
  W( 97 ) = W( 97 ) + a*JVS( 363 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  a = -W( 81 ) / JVS(          372  )
  W( 81 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 373 )
  W( 87 ) = W( 87 ) + a*JVS( 374 )
  W( 88 ) = W( 88 ) + a*JVS( 375 )
  W( 93 ) = W( 93 ) + a*JVS( 376 )
  W( 96 ) = W( 96 ) + a*JVS( 377 )
  W( 97 ) = W( 97 ) + a*JVS( 378 )
  a = -W( 82 ) / JVS(          392  )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 393 )
  W( 86 ) = W( 86 ) + a*JVS( 394 )
  W( 87 ) = W( 87 ) + a*JVS( 395 )
  W( 88 ) = W( 88 ) + a*JVS( 396 )
  W( 89 ) = W( 89 ) + a*JVS( 397 )
  W( 90 ) = W( 90 ) + a*JVS( 398 )
  W( 93 ) = W( 93 ) + a*JVS( 399 )
  W( 95 ) = W( 95 ) + a*JVS( 400 )
  W( 96 ) = W( 96 ) + a*JVS( 401 )
  W( 97 ) = W( 97 ) + a*JVS( 402 )
  a = -W( 83 ) / JVS(          414  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 415 )
  W( 86 ) = W( 86 ) + a*JVS( 416 )
  W( 87 ) = W( 87 ) + a*JVS( 417 )
  W( 88 ) = W( 88 ) + a*JVS( 418 )
  W( 89 ) = W( 89 ) + a*JVS( 419 )
  W( 90 ) = W( 90 ) + a*JVS( 420 )
  W( 91 ) = W( 91 ) + a*JVS( 421 )
  W( 92 ) = W( 92 ) + a*JVS( 422 )
  W( 93 ) = W( 93 ) + a*JVS( 423 )
  W( 94 ) = W( 94 ) + a*JVS( 424 )
  W( 95 ) = W( 95 ) + a*JVS( 425 )
  W( 96 ) = W( 96 ) + a*JVS( 426 )
  W( 97 ) = W( 97 ) + a*JVS( 427 )
  a = -W( 84 ) / JVS(          437  )
  W( 84 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 88 ) = W( 88 ) + a*JVS( 440 )
  W( 89 ) = W( 89 ) + a*JVS( 441 )
  W( 90 ) = W( 90 ) + a*JVS( 442 )
  W( 93 ) = W( 93 ) + a*JVS( 443 )
  W( 95 ) = W( 95 ) + a*JVS( 444 )
  W( 96 ) = W( 96 ) + a*JVS( 445 )
  W( 97 ) = W( 97 ) + a*JVS( 446 )
  a = -W( 85 ) / JVS(          460  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 461 )
  W( 87 ) = W( 87 ) + a*JVS( 462 )
  W( 88 ) = W( 88 ) + a*JVS( 463 )
  W( 89 ) = W( 89 ) + a*JVS( 464 )
  W( 90 ) = W( 90 ) + a*JVS( 465 )
  W( 91 ) = W( 91 ) + a*JVS( 466 )
  W( 92 ) = W( 92 ) + a*JVS( 467 )
  W( 93 ) = W( 93 ) + a*JVS( 468 )
  W( 94 ) = W( 94 ) + a*JVS( 469 )
  W( 95 ) = W( 95 ) + a*JVS( 470 )
  W( 96 ) = W( 96 ) + a*JVS( 471 )
  W( 97 ) = W( 97 ) + a*JVS( 472 )
  a = -W( 86 ) / JVS(          483  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 90 ) = W( 90 ) + a*JVS( 486 )
  W( 91 ) = W( 91 ) + a*JVS( 487 )
  W( 93 ) = W( 93 ) + a*JVS( 488 )
  W( 94 ) = W( 94 ) + a*JVS( 489 )
  W( 95 ) = W( 95 ) + a*JVS( 490 )
  W( 96 ) = W( 96 ) + a*JVS( 491 )
  W( 97 ) = W( 97 ) + a*JVS( 492 )
  a = -W( 87 ) / JVS(          517  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 518 )
  W( 89 ) = W( 89 ) + a*JVS( 519 )
  W( 90 ) = W( 90 ) + a*JVS( 520 )
  W( 91 ) = W( 91 ) + a*JVS( 521 )
  W( 92 ) = W( 92 ) + a*JVS( 522 )
  W( 93 ) = W( 93 ) + a*JVS( 523 )
  W( 94 ) = W( 94 ) + a*JVS( 524 )
  W( 95 ) = W( 95 ) + a*JVS( 525 )
  W( 96 ) = W( 96 ) + a*JVS( 526 )
  W( 97 ) = W( 97 ) + a*JVS( 527 )
  JVS( 528) = W( 39 )
  JVS( 529) = W( 42 )
  JVS( 530) = W( 48 )
  JVS( 531) = W( 57 )
  JVS( 532) = W( 58 )
  JVS( 533) = W( 59 )
  JVS( 534) = W( 61 )
  JVS( 535) = W( 64 )
  JVS( 536) = W( 67 )
  JVS( 537) = W( 74 )
  JVS( 538) = W( 76 )
  JVS( 539) = W( 77 )
  JVS( 540) = W( 79 )
  JVS( 541) = W( 80 )
  JVS( 542) = W( 81 )
  JVS( 543) = W( 82 )
  JVS( 544) = W( 83 )
  JVS( 545) = W( 84 )
  JVS( 546) = W( 85 )
  JVS( 547) = W( 86 )
  JVS( 548) = W( 87 )
  JVS( 549) = W( 88 )
  JVS( 550) = W( 89 )
  JVS( 551) = W( 90 )
  JVS( 552) = W( 91 )
  JVS( 553) = W( 92 )
  JVS( 554) = W( 93 )
  JVS( 555) = W( 94 )
  JVS( 556) = W( 95 )
  JVS( 557) = W( 96 )
  JVS( 558) = W( 97 )
  IF ( ABS(  JVS( 582 )) < TINY(a) ) THEN
         IER = 89                                      
         RETURN
  END IF
   W( 43 ) = JVS( 559 )
   W( 47 ) = JVS( 560 )
   W( 57 ) = JVS( 561 )
   W( 59 ) = JVS( 562 )
   W( 65 ) = JVS( 563 )
   W( 66 ) = JVS( 564 )
   W( 68 ) = JVS( 565 )
   W( 71 ) = JVS( 566 )
   W( 72 ) = JVS( 567 )
   W( 73 ) = JVS( 568 )
   W( 74 ) = JVS( 569 )
   W( 75 ) = JVS( 570 )
   W( 76 ) = JVS( 571 )
   W( 77 ) = JVS( 572 )
   W( 78 ) = JVS( 573 )
   W( 79 ) = JVS( 574 )
   W( 80 ) = JVS( 575 )
   W( 81 ) = JVS( 576 )
   W( 84 ) = JVS( 577 )
   W( 85 ) = JVS( 578 )
   W( 86 ) = JVS( 579 )
   W( 87 ) = JVS( 580 )
   W( 88 ) = JVS( 581 )
   W( 89 ) = JVS( 582 )
   W( 90 ) = JVS( 583 )
   W( 91 ) = JVS( 584 )
   W( 92 ) = JVS( 585 )
   W( 93 ) = JVS( 586 )
   W( 94 ) = JVS( 587 )
   W( 95 ) = JVS( 588 )
   W( 96 ) = JVS( 589 )
   W( 97 ) = JVS( 590 )
  a = -W( 43 ) / JVS(          132  )
  W( 43 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 133 )
  a = -W( 47 ) / JVS(          144  )
  W( 47 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 145 )
  a = -W( 57 ) / JVS(          188  )
  W( 57 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 189 )
  W( 96 ) = W( 96 ) + a*JVS( 190 )
  a = -W( 59 ) / JVS(          195  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 196 )
  W( 96 ) = W( 96 ) + a*JVS( 197 )
  a = -W( 65 ) / JVS(          233  )
  W( 65 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 234 )
  W( 89 ) = W( 89 ) + a*JVS( 235 )
  W( 90 ) = W( 90 ) + a*JVS( 236 )
  W( 96 ) = W( 96 ) + a*JVS( 237 )
  a = -W( 66 ) / JVS(          239  )
  W( 66 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 240 )
  W( 81 ) = W( 81 ) + a*JVS( 241 )
  W( 86 ) = W( 86 ) + a*JVS( 242 )
  W( 96 ) = W( 96 ) + a*JVS( 243 )
  a = -W( 68 ) / JVS(          252  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 253 )
  W( 90 ) = W( 90 ) + a*JVS( 254 )
  W( 92 ) = W( 92 ) + a*JVS( 255 )
  W( 95 ) = W( 95 ) + a*JVS( 256 )
  W( 96 ) = W( 96 ) + a*JVS( 257 )
  a = -W( 71 ) / JVS(          289  )
  W( 71 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 290 )
  W( 86 ) = W( 86 ) + a*JVS( 291 )
  W( 95 ) = W( 95 ) + a*JVS( 292 )
  W( 96 ) = W( 96 ) + a*JVS( 293 )
  a = -W( 72 ) / JVS(          294  )
  W( 72 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 295 )
  W( 87 ) = W( 87 ) + a*JVS( 296 )
  W( 88 ) = W( 88 ) + a*JVS( 297 )
  W( 96 ) = W( 96 ) + a*JVS( 298 )
  a = -W( 73 ) / JVS(          303  )
  W( 73 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 304 )
  W( 90 ) = W( 90 ) + a*JVS( 305 )
  W( 95 ) = W( 95 ) + a*JVS( 306 )
  W( 96 ) = W( 96 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  a = -W( 74 ) / JVS(          309  )
  W( 74 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 310 )
  W( 88 ) = W( 88 ) + a*JVS( 311 )
  W( 93 ) = W( 93 ) + a*JVS( 312 )
  W( 96 ) = W( 96 ) + a*JVS( 313 )
  W( 97 ) = W( 97 ) + a*JVS( 314 )
  a = -W( 75 ) / JVS(          315  )
  W( 75 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 316 )
  W( 93 ) = W( 93 ) + a*JVS( 317 )
  W( 96 ) = W( 96 ) + a*JVS( 318 )
  W( 97 ) = W( 97 ) + a*JVS( 319 )
  a = -W( 76 ) / JVS(          320  )
  W( 76 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 93 ) = W( 93 ) + a*JVS( 323 )
  W( 96 ) = W( 96 ) + a*JVS( 324 )
  W( 97 ) = W( 97 ) + a*JVS( 325 )
  a = -W( 77 ) / JVS(          327  )
  W( 77 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 328 )
  W( 88 ) = W( 88 ) + a*JVS( 329 )
  W( 93 ) = W( 93 ) + a*JVS( 330 )
  W( 96 ) = W( 96 ) + a*JVS( 331 )
  W( 97 ) = W( 97 ) + a*JVS( 332 )
  a = -W( 78 ) / JVS(          337  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 86 ) = W( 86 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 88 ) = W( 88 ) + a*JVS( 342 )
  W( 89 ) = W( 89 ) + a*JVS( 343 )
  W( 90 ) = W( 90 ) + a*JVS( 344 )
  W( 93 ) = W( 93 ) + a*JVS( 345 )
  W( 95 ) = W( 95 ) + a*JVS( 346 )
  W( 96 ) = W( 96 ) + a*JVS( 347 )
  W( 97 ) = W( 97 ) + a*JVS( 348 )
  a = -W( 79 ) / JVS(          354  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 355 )
  W( 81 ) = W( 81 ) + a*JVS( 356 )
  W( 84 ) = W( 84 ) + a*JVS( 357 )
  W( 86 ) = W( 86 ) + a*JVS( 358 )
  W( 87 ) = W( 87 ) + a*JVS( 359 )
  W( 88 ) = W( 88 ) + a*JVS( 360 )
  W( 93 ) = W( 93 ) + a*JVS( 361 )
  W( 96 ) = W( 96 ) + a*JVS( 362 )
  W( 97 ) = W( 97 ) + a*JVS( 363 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  a = -W( 81 ) / JVS(          372  )
  W( 81 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 373 )
  W( 87 ) = W( 87 ) + a*JVS( 374 )
  W( 88 ) = W( 88 ) + a*JVS( 375 )
  W( 93 ) = W( 93 ) + a*JVS( 376 )
  W( 96 ) = W( 96 ) + a*JVS( 377 )
  W( 97 ) = W( 97 ) + a*JVS( 378 )
  a = -W( 84 ) / JVS(          437  )
  W( 84 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 88 ) = W( 88 ) + a*JVS( 440 )
  W( 89 ) = W( 89 ) + a*JVS( 441 )
  W( 90 ) = W( 90 ) + a*JVS( 442 )
  W( 93 ) = W( 93 ) + a*JVS( 443 )
  W( 95 ) = W( 95 ) + a*JVS( 444 )
  W( 96 ) = W( 96 ) + a*JVS( 445 )
  W( 97 ) = W( 97 ) + a*JVS( 446 )
  a = -W( 85 ) / JVS(          460  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 461 )
  W( 87 ) = W( 87 ) + a*JVS( 462 )
  W( 88 ) = W( 88 ) + a*JVS( 463 )
  W( 89 ) = W( 89 ) + a*JVS( 464 )
  W( 90 ) = W( 90 ) + a*JVS( 465 )
  W( 91 ) = W( 91 ) + a*JVS( 466 )
  W( 92 ) = W( 92 ) + a*JVS( 467 )
  W( 93 ) = W( 93 ) + a*JVS( 468 )
  W( 94 ) = W( 94 ) + a*JVS( 469 )
  W( 95 ) = W( 95 ) + a*JVS( 470 )
  W( 96 ) = W( 96 ) + a*JVS( 471 )
  W( 97 ) = W( 97 ) + a*JVS( 472 )
  a = -W( 86 ) / JVS(          483  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 90 ) = W( 90 ) + a*JVS( 486 )
  W( 91 ) = W( 91 ) + a*JVS( 487 )
  W( 93 ) = W( 93 ) + a*JVS( 488 )
  W( 94 ) = W( 94 ) + a*JVS( 489 )
  W( 95 ) = W( 95 ) + a*JVS( 490 )
  W( 96 ) = W( 96 ) + a*JVS( 491 )
  W( 97 ) = W( 97 ) + a*JVS( 492 )
  a = -W( 87 ) / JVS(          517  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 518 )
  W( 89 ) = W( 89 ) + a*JVS( 519 )
  W( 90 ) = W( 90 ) + a*JVS( 520 )
  W( 91 ) = W( 91 ) + a*JVS( 521 )
  W( 92 ) = W( 92 ) + a*JVS( 522 )
  W( 93 ) = W( 93 ) + a*JVS( 523 )
  W( 94 ) = W( 94 ) + a*JVS( 524 )
  W( 95 ) = W( 95 ) + a*JVS( 525 )
  W( 96 ) = W( 96 ) + a*JVS( 526 )
  W( 97 ) = W( 97 ) + a*JVS( 527 )
  a = -W( 88 ) / JVS(          549  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 550 )
  W( 90 ) = W( 90 ) + a*JVS( 551 )
  W( 91 ) = W( 91 ) + a*JVS( 552 )
  W( 92 ) = W( 92 ) + a*JVS( 553 )
  W( 93 ) = W( 93 ) + a*JVS( 554 )
  W( 94 ) = W( 94 ) + a*JVS( 555 )
  W( 95 ) = W( 95 ) + a*JVS( 556 )
  W( 96 ) = W( 96 ) + a*JVS( 557 )
  W( 97 ) = W( 97 ) + a*JVS( 558 )
  JVS( 559) = W( 43 )
  JVS( 560) = W( 47 )
  JVS( 561) = W( 57 )
  JVS( 562) = W( 59 )
  JVS( 563) = W( 65 )
  JVS( 564) = W( 66 )
  JVS( 565) = W( 68 )
  JVS( 566) = W( 71 )
  JVS( 567) = W( 72 )
  JVS( 568) = W( 73 )
  JVS( 569) = W( 74 )
  JVS( 570) = W( 75 )
  JVS( 571) = W( 76 )
  JVS( 572) = W( 77 )
  JVS( 573) = W( 78 )
  JVS( 574) = W( 79 )
  JVS( 575) = W( 80 )
  JVS( 576) = W( 81 )
  JVS( 577) = W( 84 )
  JVS( 578) = W( 85 )
  JVS( 579) = W( 86 )
  JVS( 580) = W( 87 )
  JVS( 581) = W( 88 )
  JVS( 582) = W( 89 )
  JVS( 583) = W( 90 )
  JVS( 584) = W( 91 )
  JVS( 585) = W( 92 )
  JVS( 586) = W( 93 )
  JVS( 587) = W( 94 )
  JVS( 588) = W( 95 )
  JVS( 589) = W( 96 )
  JVS( 590) = W( 97 )
  IF ( ABS(  JVS( 630 )) < TINY(a) ) THEN
         IER = 90                                      
         RETURN
  END IF
   W( 40 ) = JVS( 591 )
   W( 43 ) = JVS( 592 )
   W( 44 ) = JVS( 593 )
   W( 45 ) = JVS( 594 )
   W( 47 ) = JVS( 595 )
   W( 51 ) = JVS( 596 )
   W( 52 ) = JVS( 597 )
   W( 53 ) = JVS( 598 )
   W( 57 ) = JVS( 599 )
   W( 58 ) = JVS( 600 )
   W( 59 ) = JVS( 601 )
   W( 60 ) = JVS( 602 )
   W( 62 ) = JVS( 603 )
   W( 63 ) = JVS( 604 )
   W( 64 ) = JVS( 605 )
   W( 65 ) = JVS( 606 )
   W( 66 ) = JVS( 607 )
   W( 67 ) = JVS( 608 )
   W( 68 ) = JVS( 609 )
   W( 70 ) = JVS( 610 )
   W( 71 ) = JVS( 611 )
   W( 72 ) = JVS( 612 )
   W( 73 ) = JVS( 613 )
   W( 74 ) = JVS( 614 )
   W( 75 ) = JVS( 615 )
   W( 76 ) = JVS( 616 )
   W( 77 ) = JVS( 617 )
   W( 78 ) = JVS( 618 )
   W( 79 ) = JVS( 619 )
   W( 80 ) = JVS( 620 )
   W( 81 ) = JVS( 621 )
   W( 82 ) = JVS( 622 )
   W( 83 ) = JVS( 623 )
   W( 84 ) = JVS( 624 )
   W( 85 ) = JVS( 625 )
   W( 86 ) = JVS( 626 )
   W( 87 ) = JVS( 627 )
   W( 88 ) = JVS( 628 )
   W( 89 ) = JVS( 629 )
   W( 90 ) = JVS( 630 )
   W( 91 ) = JVS( 631 )
   W( 92 ) = JVS( 632 )
   W( 93 ) = JVS( 633 )
   W( 94 ) = JVS( 634 )
   W( 95 ) = JVS( 635 )
   W( 96 ) = JVS( 636 )
   W( 97 ) = JVS( 637 )
  a = -W( 40 ) / JVS(          124  )
  W( 40 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 125 )
  a = -W( 43 ) / JVS(          132  )
  W( 43 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 133 )
  a = -W( 44 ) / JVS(          134  )
  W( 44 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 135 )
  W( 96 ) = W( 96 ) + a*JVS( 136 )
  a = -W( 45 ) / JVS(          138  )
  W( 45 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 139 )
  W( 96 ) = W( 96 ) + a*JVS( 140 )
  a = -W( 47 ) / JVS(          144  )
  W( 47 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 145 )
  a = -W( 51 ) / JVS(          157  )
  W( 51 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 158 )
  W( 74 ) = W( 74 ) + a*JVS( 159 )
  W( 86 ) = W( 86 ) + a*JVS( 160 )
  W( 95 ) = W( 95 ) + a*JVS( 161 )
  W( 96 ) = W( 96 ) + a*JVS( 162 )
  a = -W( 52 ) / JVS(          165  )
  W( 52 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 166 )
  W( 96 ) = W( 96 ) + a*JVS( 167 )
  a = -W( 53 ) / JVS(          168  )
  W( 53 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 169 )
  W( 90 ) = W( 90 ) + a*JVS( 170 )
  W( 96 ) = W( 96 ) + a*JVS( 171 )
  a = -W( 57 ) / JVS(          188  )
  W( 57 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 189 )
  W( 96 ) = W( 96 ) + a*JVS( 190 )
  a = -W( 58 ) / JVS(          191  )
  W( 58 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 192 )
  W( 92 ) = W( 92 ) + a*JVS( 193 )
  W( 96 ) = W( 96 ) + a*JVS( 194 )
  a = -W( 59 ) / JVS(          195  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 196 )
  W( 96 ) = W( 96 ) + a*JVS( 197 )
  a = -W( 60 ) / JVS(          199  )
  W( 60 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 200 )
  W( 90 ) = W( 90 ) + a*JVS( 201 )
  W( 93 ) = W( 93 ) + a*JVS( 202 )
  W( 96 ) = W( 96 ) + a*JVS( 203 )
  a = -W( 62 ) / JVS(          218  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 219 )
  W( 90 ) = W( 90 ) + a*JVS( 220 )
  W( 95 ) = W( 95 ) + a*JVS( 221 )
  a = -W( 63 ) / JVS(          222  )
  W( 63 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 223 )
  W( 87 ) = W( 87 ) + a*JVS( 224 )
  W( 90 ) = W( 90 ) + a*JVS( 225 )
  W( 96 ) = W( 96 ) + a*JVS( 226 )
  W( 97 ) = W( 97 ) + a*JVS( 227 )
  a = -W( 64 ) / JVS(          228  )
  W( 64 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 229 )
  W( 88 ) = W( 88 ) + a*JVS( 230 )
  W( 90 ) = W( 90 ) + a*JVS( 231 )
  W( 95 ) = W( 95 ) + a*JVS( 232 )
  a = -W( 65 ) / JVS(          233  )
  W( 65 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 234 )
  W( 89 ) = W( 89 ) + a*JVS( 235 )
  W( 90 ) = W( 90 ) + a*JVS( 236 )
  W( 96 ) = W( 96 ) + a*JVS( 237 )
  a = -W( 66 ) / JVS(          239  )
  W( 66 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 240 )
  W( 81 ) = W( 81 ) + a*JVS( 241 )
  W( 86 ) = W( 86 ) + a*JVS( 242 )
  W( 96 ) = W( 96 ) + a*JVS( 243 )
  a = -W( 67 ) / JVS(          244  )
  W( 67 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 245 )
  W( 76 ) = W( 76 ) + a*JVS( 246 )
  W( 77 ) = W( 77 ) + a*JVS( 247 )
  W( 80 ) = W( 80 ) + a*JVS( 248 )
  W( 88 ) = W( 88 ) + a*JVS( 249 )
  W( 96 ) = W( 96 ) + a*JVS( 250 )
  a = -W( 68 ) / JVS(          252  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 253 )
  W( 90 ) = W( 90 ) + a*JVS( 254 )
  W( 92 ) = W( 92 ) + a*JVS( 255 )
  W( 95 ) = W( 95 ) + a*JVS( 256 )
  W( 96 ) = W( 96 ) + a*JVS( 257 )
  a = -W( 70 ) / JVS(          272  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 273 )
  W( 74 ) = W( 74 ) + a*JVS( 274 )
  W( 75 ) = W( 75 ) + a*JVS( 275 )
  W( 76 ) = W( 76 ) + a*JVS( 276 )
  W( 77 ) = W( 77 ) + a*JVS( 277 )
  W( 80 ) = W( 80 ) + a*JVS( 278 )
  W( 81 ) = W( 81 ) + a*JVS( 279 )
  W( 82 ) = W( 82 ) + a*JVS( 280 )
  W( 83 ) = W( 83 ) + a*JVS( 281 )
  W( 85 ) = W( 85 ) + a*JVS( 282 )
  W( 86 ) = W( 86 ) + a*JVS( 283 )
  W( 88 ) = W( 88 ) + a*JVS( 284 )
  W( 93 ) = W( 93 ) + a*JVS( 285 )
  W( 96 ) = W( 96 ) + a*JVS( 286 )
  W( 97 ) = W( 97 ) + a*JVS( 287 )
  a = -W( 71 ) / JVS(          289  )
  W( 71 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 290 )
  W( 86 ) = W( 86 ) + a*JVS( 291 )
  W( 95 ) = W( 95 ) + a*JVS( 292 )
  W( 96 ) = W( 96 ) + a*JVS( 293 )
  a = -W( 72 ) / JVS(          294  )
  W( 72 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 295 )
  W( 87 ) = W( 87 ) + a*JVS( 296 )
  W( 88 ) = W( 88 ) + a*JVS( 297 )
  W( 96 ) = W( 96 ) + a*JVS( 298 )
  a = -W( 73 ) / JVS(          303  )
  W( 73 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 304 )
  W( 90 ) = W( 90 ) + a*JVS( 305 )
  W( 95 ) = W( 95 ) + a*JVS( 306 )
  W( 96 ) = W( 96 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  a = -W( 74 ) / JVS(          309  )
  W( 74 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 310 )
  W( 88 ) = W( 88 ) + a*JVS( 311 )
  W( 93 ) = W( 93 ) + a*JVS( 312 )
  W( 96 ) = W( 96 ) + a*JVS( 313 )
  W( 97 ) = W( 97 ) + a*JVS( 314 )
  a = -W( 75 ) / JVS(          315  )
  W( 75 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 316 )
  W( 93 ) = W( 93 ) + a*JVS( 317 )
  W( 96 ) = W( 96 ) + a*JVS( 318 )
  W( 97 ) = W( 97 ) + a*JVS( 319 )
  a = -W( 76 ) / JVS(          320  )
  W( 76 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 93 ) = W( 93 ) + a*JVS( 323 )
  W( 96 ) = W( 96 ) + a*JVS( 324 )
  W( 97 ) = W( 97 ) + a*JVS( 325 )
  a = -W( 77 ) / JVS(          327  )
  W( 77 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 328 )
  W( 88 ) = W( 88 ) + a*JVS( 329 )
  W( 93 ) = W( 93 ) + a*JVS( 330 )
  W( 96 ) = W( 96 ) + a*JVS( 331 )
  W( 97 ) = W( 97 ) + a*JVS( 332 )
  a = -W( 78 ) / JVS(          337  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 86 ) = W( 86 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 88 ) = W( 88 ) + a*JVS( 342 )
  W( 89 ) = W( 89 ) + a*JVS( 343 )
  W( 90 ) = W( 90 ) + a*JVS( 344 )
  W( 93 ) = W( 93 ) + a*JVS( 345 )
  W( 95 ) = W( 95 ) + a*JVS( 346 )
  W( 96 ) = W( 96 ) + a*JVS( 347 )
  W( 97 ) = W( 97 ) + a*JVS( 348 )
  a = -W( 79 ) / JVS(          354  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 355 )
  W( 81 ) = W( 81 ) + a*JVS( 356 )
  W( 84 ) = W( 84 ) + a*JVS( 357 )
  W( 86 ) = W( 86 ) + a*JVS( 358 )
  W( 87 ) = W( 87 ) + a*JVS( 359 )
  W( 88 ) = W( 88 ) + a*JVS( 360 )
  W( 93 ) = W( 93 ) + a*JVS( 361 )
  W( 96 ) = W( 96 ) + a*JVS( 362 )
  W( 97 ) = W( 97 ) + a*JVS( 363 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  a = -W( 81 ) / JVS(          372  )
  W( 81 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 373 )
  W( 87 ) = W( 87 ) + a*JVS( 374 )
  W( 88 ) = W( 88 ) + a*JVS( 375 )
  W( 93 ) = W( 93 ) + a*JVS( 376 )
  W( 96 ) = W( 96 ) + a*JVS( 377 )
  W( 97 ) = W( 97 ) + a*JVS( 378 )
  a = -W( 82 ) / JVS(          392  )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 393 )
  W( 86 ) = W( 86 ) + a*JVS( 394 )
  W( 87 ) = W( 87 ) + a*JVS( 395 )
  W( 88 ) = W( 88 ) + a*JVS( 396 )
  W( 89 ) = W( 89 ) + a*JVS( 397 )
  W( 90 ) = W( 90 ) + a*JVS( 398 )
  W( 93 ) = W( 93 ) + a*JVS( 399 )
  W( 95 ) = W( 95 ) + a*JVS( 400 )
  W( 96 ) = W( 96 ) + a*JVS( 401 )
  W( 97 ) = W( 97 ) + a*JVS( 402 )
  a = -W( 83 ) / JVS(          414  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 415 )
  W( 86 ) = W( 86 ) + a*JVS( 416 )
  W( 87 ) = W( 87 ) + a*JVS( 417 )
  W( 88 ) = W( 88 ) + a*JVS( 418 )
  W( 89 ) = W( 89 ) + a*JVS( 419 )
  W( 90 ) = W( 90 ) + a*JVS( 420 )
  W( 91 ) = W( 91 ) + a*JVS( 421 )
  W( 92 ) = W( 92 ) + a*JVS( 422 )
  W( 93 ) = W( 93 ) + a*JVS( 423 )
  W( 94 ) = W( 94 ) + a*JVS( 424 )
  W( 95 ) = W( 95 ) + a*JVS( 425 )
  W( 96 ) = W( 96 ) + a*JVS( 426 )
  W( 97 ) = W( 97 ) + a*JVS( 427 )
  a = -W( 84 ) / JVS(          437  )
  W( 84 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 88 ) = W( 88 ) + a*JVS( 440 )
  W( 89 ) = W( 89 ) + a*JVS( 441 )
  W( 90 ) = W( 90 ) + a*JVS( 442 )
  W( 93 ) = W( 93 ) + a*JVS( 443 )
  W( 95 ) = W( 95 ) + a*JVS( 444 )
  W( 96 ) = W( 96 ) + a*JVS( 445 )
  W( 97 ) = W( 97 ) + a*JVS( 446 )
  a = -W( 85 ) / JVS(          460  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 461 )
  W( 87 ) = W( 87 ) + a*JVS( 462 )
  W( 88 ) = W( 88 ) + a*JVS( 463 )
  W( 89 ) = W( 89 ) + a*JVS( 464 )
  W( 90 ) = W( 90 ) + a*JVS( 465 )
  W( 91 ) = W( 91 ) + a*JVS( 466 )
  W( 92 ) = W( 92 ) + a*JVS( 467 )
  W( 93 ) = W( 93 ) + a*JVS( 468 )
  W( 94 ) = W( 94 ) + a*JVS( 469 )
  W( 95 ) = W( 95 ) + a*JVS( 470 )
  W( 96 ) = W( 96 ) + a*JVS( 471 )
  W( 97 ) = W( 97 ) + a*JVS( 472 )
  a = -W( 86 ) / JVS(          483  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 90 ) = W( 90 ) + a*JVS( 486 )
  W( 91 ) = W( 91 ) + a*JVS( 487 )
  W( 93 ) = W( 93 ) + a*JVS( 488 )
  W( 94 ) = W( 94 ) + a*JVS( 489 )
  W( 95 ) = W( 95 ) + a*JVS( 490 )
  W( 96 ) = W( 96 ) + a*JVS( 491 )
  W( 97 ) = W( 97 ) + a*JVS( 492 )
  a = -W( 87 ) / JVS(          517  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 518 )
  W( 89 ) = W( 89 ) + a*JVS( 519 )
  W( 90 ) = W( 90 ) + a*JVS( 520 )
  W( 91 ) = W( 91 ) + a*JVS( 521 )
  W( 92 ) = W( 92 ) + a*JVS( 522 )
  W( 93 ) = W( 93 ) + a*JVS( 523 )
  W( 94 ) = W( 94 ) + a*JVS( 524 )
  W( 95 ) = W( 95 ) + a*JVS( 525 )
  W( 96 ) = W( 96 ) + a*JVS( 526 )
  W( 97 ) = W( 97 ) + a*JVS( 527 )
  a = -W( 88 ) / JVS(          549  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 550 )
  W( 90 ) = W( 90 ) + a*JVS( 551 )
  W( 91 ) = W( 91 ) + a*JVS( 552 )
  W( 92 ) = W( 92 ) + a*JVS( 553 )
  W( 93 ) = W( 93 ) + a*JVS( 554 )
  W( 94 ) = W( 94 ) + a*JVS( 555 )
  W( 95 ) = W( 95 ) + a*JVS( 556 )
  W( 96 ) = W( 96 ) + a*JVS( 557 )
  W( 97 ) = W( 97 ) + a*JVS( 558 )
  a = -W( 89 ) / JVS(          582  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 583 )
  W( 91 ) = W( 91 ) + a*JVS( 584 )
  W( 92 ) = W( 92 ) + a*JVS( 585 )
  W( 93 ) = W( 93 ) + a*JVS( 586 )
  W( 94 ) = W( 94 ) + a*JVS( 587 )
  W( 95 ) = W( 95 ) + a*JVS( 588 )
  W( 96 ) = W( 96 ) + a*JVS( 589 )
  W( 97 ) = W( 97 ) + a*JVS( 590 )
  JVS( 591) = W( 40 )
  JVS( 592) = W( 43 )
  JVS( 593) = W( 44 )
  JVS( 594) = W( 45 )
  JVS( 595) = W( 47 )
  JVS( 596) = W( 51 )
  JVS( 597) = W( 52 )
  JVS( 598) = W( 53 )
  JVS( 599) = W( 57 )
  JVS( 600) = W( 58 )
  JVS( 601) = W( 59 )
  JVS( 602) = W( 60 )
  JVS( 603) = W( 62 )
  JVS( 604) = W( 63 )
  JVS( 605) = W( 64 )
  JVS( 606) = W( 65 )
  JVS( 607) = W( 66 )
  JVS( 608) = W( 67 )
  JVS( 609) = W( 68 )
  JVS( 610) = W( 70 )
  JVS( 611) = W( 71 )
  JVS( 612) = W( 72 )
  JVS( 613) = W( 73 )
  JVS( 614) = W( 74 )
  JVS( 615) = W( 75 )
  JVS( 616) = W( 76 )
  JVS( 617) = W( 77 )
  JVS( 618) = W( 78 )
  JVS( 619) = W( 79 )
  JVS( 620) = W( 80 )
  JVS( 621) = W( 81 )
  JVS( 622) = W( 82 )
  JVS( 623) = W( 83 )
  JVS( 624) = W( 84 )
  JVS( 625) = W( 85 )
  JVS( 626) = W( 86 )
  JVS( 627) = W( 87 )
  JVS( 628) = W( 88 )
  JVS( 629) = W( 89 )
  JVS( 630) = W( 90 )
  JVS( 631) = W( 91 )
  JVS( 632) = W( 92 )
  JVS( 633) = W( 93 )
  JVS( 634) = W( 94 )
  JVS( 635) = W( 95 )
  JVS( 636) = W( 96 )
  JVS( 637) = W( 97 )
  IF ( ABS(  JVS( 651 )) < TINY(a) ) THEN
         IER = 91                                      
         RETURN
  END IF
   W( 41 ) = JVS( 638 )
   W( 54 ) = JVS( 639 )
   W( 66 ) = JVS( 640 )
   W( 71 ) = JVS( 641 )
   W( 73 ) = JVS( 642 )
   W( 81 ) = JVS( 643 )
   W( 83 ) = JVS( 644 )
   W( 84 ) = JVS( 645 )
   W( 86 ) = JVS( 646 )
   W( 87 ) = JVS( 647 )
   W( 88 ) = JVS( 648 )
   W( 89 ) = JVS( 649 )
   W( 90 ) = JVS( 650 )
   W( 91 ) = JVS( 651 )
   W( 92 ) = JVS( 652 )
   W( 93 ) = JVS( 653 )
   W( 94 ) = JVS( 654 )
   W( 95 ) = JVS( 655 )
   W( 96 ) = JVS( 656 )
   W( 97 ) = JVS( 657 )
  a = -W( 41 ) / JVS(          126  )
  W( 41 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 127 )
  W( 91 ) = W( 91 ) + a*JVS( 128 )
  a = -W( 54 ) / JVS(          172  )
  W( 54 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 173 )
  W( 91 ) = W( 91 ) + a*JVS( 174 )
  W( 94 ) = W( 94 ) + a*JVS( 175 )
  W( 96 ) = W( 96 ) + a*JVS( 176 )
  a = -W( 66 ) / JVS(          239  )
  W( 66 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 240 )
  W( 81 ) = W( 81 ) + a*JVS( 241 )
  W( 86 ) = W( 86 ) + a*JVS( 242 )
  W( 96 ) = W( 96 ) + a*JVS( 243 )
  a = -W( 71 ) / JVS(          289  )
  W( 71 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 290 )
  W( 86 ) = W( 86 ) + a*JVS( 291 )
  W( 95 ) = W( 95 ) + a*JVS( 292 )
  W( 96 ) = W( 96 ) + a*JVS( 293 )
  a = -W( 73 ) / JVS(          303  )
  W( 73 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 304 )
  W( 90 ) = W( 90 ) + a*JVS( 305 )
  W( 95 ) = W( 95 ) + a*JVS( 306 )
  W( 96 ) = W( 96 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  a = -W( 81 ) / JVS(          372  )
  W( 81 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 373 )
  W( 87 ) = W( 87 ) + a*JVS( 374 )
  W( 88 ) = W( 88 ) + a*JVS( 375 )
  W( 93 ) = W( 93 ) + a*JVS( 376 )
  W( 96 ) = W( 96 ) + a*JVS( 377 )
  W( 97 ) = W( 97 ) + a*JVS( 378 )
  a = -W( 83 ) / JVS(          414  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 415 )
  W( 86 ) = W( 86 ) + a*JVS( 416 )
  W( 87 ) = W( 87 ) + a*JVS( 417 )
  W( 88 ) = W( 88 ) + a*JVS( 418 )
  W( 89 ) = W( 89 ) + a*JVS( 419 )
  W( 90 ) = W( 90 ) + a*JVS( 420 )
  W( 91 ) = W( 91 ) + a*JVS( 421 )
  W( 92 ) = W( 92 ) + a*JVS( 422 )
  W( 93 ) = W( 93 ) + a*JVS( 423 )
  W( 94 ) = W( 94 ) + a*JVS( 424 )
  W( 95 ) = W( 95 ) + a*JVS( 425 )
  W( 96 ) = W( 96 ) + a*JVS( 426 )
  W( 97 ) = W( 97 ) + a*JVS( 427 )
  a = -W( 84 ) / JVS(          437  )
  W( 84 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 88 ) = W( 88 ) + a*JVS( 440 )
  W( 89 ) = W( 89 ) + a*JVS( 441 )
  W( 90 ) = W( 90 ) + a*JVS( 442 )
  W( 93 ) = W( 93 ) + a*JVS( 443 )
  W( 95 ) = W( 95 ) + a*JVS( 444 )
  W( 96 ) = W( 96 ) + a*JVS( 445 )
  W( 97 ) = W( 97 ) + a*JVS( 446 )
  a = -W( 86 ) / JVS(          483  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 90 ) = W( 90 ) + a*JVS( 486 )
  W( 91 ) = W( 91 ) + a*JVS( 487 )
  W( 93 ) = W( 93 ) + a*JVS( 488 )
  W( 94 ) = W( 94 ) + a*JVS( 489 )
  W( 95 ) = W( 95 ) + a*JVS( 490 )
  W( 96 ) = W( 96 ) + a*JVS( 491 )
  W( 97 ) = W( 97 ) + a*JVS( 492 )
  a = -W( 87 ) / JVS(          517  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 518 )
  W( 89 ) = W( 89 ) + a*JVS( 519 )
  W( 90 ) = W( 90 ) + a*JVS( 520 )
  W( 91 ) = W( 91 ) + a*JVS( 521 )
  W( 92 ) = W( 92 ) + a*JVS( 522 )
  W( 93 ) = W( 93 ) + a*JVS( 523 )
  W( 94 ) = W( 94 ) + a*JVS( 524 )
  W( 95 ) = W( 95 ) + a*JVS( 525 )
  W( 96 ) = W( 96 ) + a*JVS( 526 )
  W( 97 ) = W( 97 ) + a*JVS( 527 )
  a = -W( 88 ) / JVS(          549  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 550 )
  W( 90 ) = W( 90 ) + a*JVS( 551 )
  W( 91 ) = W( 91 ) + a*JVS( 552 )
  W( 92 ) = W( 92 ) + a*JVS( 553 )
  W( 93 ) = W( 93 ) + a*JVS( 554 )
  W( 94 ) = W( 94 ) + a*JVS( 555 )
  W( 95 ) = W( 95 ) + a*JVS( 556 )
  W( 96 ) = W( 96 ) + a*JVS( 557 )
  W( 97 ) = W( 97 ) + a*JVS( 558 )
  a = -W( 89 ) / JVS(          582  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 583 )
  W( 91 ) = W( 91 ) + a*JVS( 584 )
  W( 92 ) = W( 92 ) + a*JVS( 585 )
  W( 93 ) = W( 93 ) + a*JVS( 586 )
  W( 94 ) = W( 94 ) + a*JVS( 587 )
  W( 95 ) = W( 95 ) + a*JVS( 588 )
  W( 96 ) = W( 96 ) + a*JVS( 589 )
  W( 97 ) = W( 97 ) + a*JVS( 590 )
  a = -W( 90 ) / JVS(          630  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 631 )
  W( 92 ) = W( 92 ) + a*JVS( 632 )
  W( 93 ) = W( 93 ) + a*JVS( 633 )
  W( 94 ) = W( 94 ) + a*JVS( 634 )
  W( 95 ) = W( 95 ) + a*JVS( 635 )
  W( 96 ) = W( 96 ) + a*JVS( 636 )
  W( 97 ) = W( 97 ) + a*JVS( 637 )
  JVS( 638) = W( 41 )
  JVS( 639) = W( 54 )
  JVS( 640) = W( 66 )
  JVS( 641) = W( 71 )
  JVS( 642) = W( 73 )
  JVS( 643) = W( 81 )
  JVS( 644) = W( 83 )
  JVS( 645) = W( 84 )
  JVS( 646) = W( 86 )
  JVS( 647) = W( 87 )
  JVS( 648) = W( 88 )
  JVS( 649) = W( 89 )
  JVS( 650) = W( 90 )
  JVS( 651) = W( 91 )
  JVS( 652) = W( 92 )
  JVS( 653) = W( 93 )
  JVS( 654) = W( 94 )
  JVS( 655) = W( 95 )
  JVS( 656) = W( 96 )
  JVS( 657) = W( 97 )
  IF ( ABS(  JVS( 672 )) < TINY(a) ) THEN
         IER = 92                                      
         RETURN
  END IF
   W( 48 ) = JVS( 658 )
   W( 54 ) = JVS( 659 )
   W( 55 ) = JVS( 660 )
   W( 68 ) = JVS( 661 )
   W( 82 ) = JVS( 662 )
   W( 83 ) = JVS( 663 )
   W( 84 ) = JVS( 664 )
   W( 85 ) = JVS( 665 )
   W( 86 ) = JVS( 666 )
   W( 87 ) = JVS( 667 )
   W( 88 ) = JVS( 668 )
   W( 89 ) = JVS( 669 )
   W( 90 ) = JVS( 670 )
   W( 91 ) = JVS( 671 )
   W( 92 ) = JVS( 672 )
   W( 93 ) = JVS( 673 )
   W( 94 ) = JVS( 674 )
   W( 95 ) = JVS( 675 )
   W( 96 ) = JVS( 676 )
   W( 97 ) = JVS( 677 )
  a = -W( 48 ) / JVS(          146  )
  W( 48 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 147 )
  W( 96 ) = W( 96 ) + a*JVS( 148 )
  a = -W( 54 ) / JVS(          172  )
  W( 54 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 173 )
  W( 91 ) = W( 91 ) + a*JVS( 174 )
  W( 94 ) = W( 94 ) + a*JVS( 175 )
  W( 96 ) = W( 96 ) + a*JVS( 176 )
  a = -W( 55 ) / JVS(          177  )
  W( 55 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 178 )
  W( 90 ) = W( 90 ) + a*JVS( 179 )
  W( 91 ) = W( 91 ) + a*JVS( 180 )
  W( 92 ) = W( 92 ) + a*JVS( 181 )
  W( 94 ) = W( 94 ) + a*JVS( 182 )
  W( 96 ) = W( 96 ) + a*JVS( 183 )
  a = -W( 68 ) / JVS(          252  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 253 )
  W( 90 ) = W( 90 ) + a*JVS( 254 )
  W( 92 ) = W( 92 ) + a*JVS( 255 )
  W( 95 ) = W( 95 ) + a*JVS( 256 )
  W( 96 ) = W( 96 ) + a*JVS( 257 )
  a = -W( 82 ) / JVS(          392  )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 393 )
  W( 86 ) = W( 86 ) + a*JVS( 394 )
  W( 87 ) = W( 87 ) + a*JVS( 395 )
  W( 88 ) = W( 88 ) + a*JVS( 396 )
  W( 89 ) = W( 89 ) + a*JVS( 397 )
  W( 90 ) = W( 90 ) + a*JVS( 398 )
  W( 93 ) = W( 93 ) + a*JVS( 399 )
  W( 95 ) = W( 95 ) + a*JVS( 400 )
  W( 96 ) = W( 96 ) + a*JVS( 401 )
  W( 97 ) = W( 97 ) + a*JVS( 402 )
  a = -W( 83 ) / JVS(          414  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 415 )
  W( 86 ) = W( 86 ) + a*JVS( 416 )
  W( 87 ) = W( 87 ) + a*JVS( 417 )
  W( 88 ) = W( 88 ) + a*JVS( 418 )
  W( 89 ) = W( 89 ) + a*JVS( 419 )
  W( 90 ) = W( 90 ) + a*JVS( 420 )
  W( 91 ) = W( 91 ) + a*JVS( 421 )
  W( 92 ) = W( 92 ) + a*JVS( 422 )
  W( 93 ) = W( 93 ) + a*JVS( 423 )
  W( 94 ) = W( 94 ) + a*JVS( 424 )
  W( 95 ) = W( 95 ) + a*JVS( 425 )
  W( 96 ) = W( 96 ) + a*JVS( 426 )
  W( 97 ) = W( 97 ) + a*JVS( 427 )
  a = -W( 84 ) / JVS(          437  )
  W( 84 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 88 ) = W( 88 ) + a*JVS( 440 )
  W( 89 ) = W( 89 ) + a*JVS( 441 )
  W( 90 ) = W( 90 ) + a*JVS( 442 )
  W( 93 ) = W( 93 ) + a*JVS( 443 )
  W( 95 ) = W( 95 ) + a*JVS( 444 )
  W( 96 ) = W( 96 ) + a*JVS( 445 )
  W( 97 ) = W( 97 ) + a*JVS( 446 )
  a = -W( 85 ) / JVS(          460  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 461 )
  W( 87 ) = W( 87 ) + a*JVS( 462 )
  W( 88 ) = W( 88 ) + a*JVS( 463 )
  W( 89 ) = W( 89 ) + a*JVS( 464 )
  W( 90 ) = W( 90 ) + a*JVS( 465 )
  W( 91 ) = W( 91 ) + a*JVS( 466 )
  W( 92 ) = W( 92 ) + a*JVS( 467 )
  W( 93 ) = W( 93 ) + a*JVS( 468 )
  W( 94 ) = W( 94 ) + a*JVS( 469 )
  W( 95 ) = W( 95 ) + a*JVS( 470 )
  W( 96 ) = W( 96 ) + a*JVS( 471 )
  W( 97 ) = W( 97 ) + a*JVS( 472 )
  a = -W( 86 ) / JVS(          483  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 90 ) = W( 90 ) + a*JVS( 486 )
  W( 91 ) = W( 91 ) + a*JVS( 487 )
  W( 93 ) = W( 93 ) + a*JVS( 488 )
  W( 94 ) = W( 94 ) + a*JVS( 489 )
  W( 95 ) = W( 95 ) + a*JVS( 490 )
  W( 96 ) = W( 96 ) + a*JVS( 491 )
  W( 97 ) = W( 97 ) + a*JVS( 492 )
  a = -W( 87 ) / JVS(          517  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 518 )
  W( 89 ) = W( 89 ) + a*JVS( 519 )
  W( 90 ) = W( 90 ) + a*JVS( 520 )
  W( 91 ) = W( 91 ) + a*JVS( 521 )
  W( 92 ) = W( 92 ) + a*JVS( 522 )
  W( 93 ) = W( 93 ) + a*JVS( 523 )
  W( 94 ) = W( 94 ) + a*JVS( 524 )
  W( 95 ) = W( 95 ) + a*JVS( 525 )
  W( 96 ) = W( 96 ) + a*JVS( 526 )
  W( 97 ) = W( 97 ) + a*JVS( 527 )
  a = -W( 88 ) / JVS(          549  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 550 )
  W( 90 ) = W( 90 ) + a*JVS( 551 )
  W( 91 ) = W( 91 ) + a*JVS( 552 )
  W( 92 ) = W( 92 ) + a*JVS( 553 )
  W( 93 ) = W( 93 ) + a*JVS( 554 )
  W( 94 ) = W( 94 ) + a*JVS( 555 )
  W( 95 ) = W( 95 ) + a*JVS( 556 )
  W( 96 ) = W( 96 ) + a*JVS( 557 )
  W( 97 ) = W( 97 ) + a*JVS( 558 )
  a = -W( 89 ) / JVS(          582  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 583 )
  W( 91 ) = W( 91 ) + a*JVS( 584 )
  W( 92 ) = W( 92 ) + a*JVS( 585 )
  W( 93 ) = W( 93 ) + a*JVS( 586 )
  W( 94 ) = W( 94 ) + a*JVS( 587 )
  W( 95 ) = W( 95 ) + a*JVS( 588 )
  W( 96 ) = W( 96 ) + a*JVS( 589 )
  W( 97 ) = W( 97 ) + a*JVS( 590 )
  a = -W( 90 ) / JVS(          630  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 631 )
  W( 92 ) = W( 92 ) + a*JVS( 632 )
  W( 93 ) = W( 93 ) + a*JVS( 633 )
  W( 94 ) = W( 94 ) + a*JVS( 634 )
  W( 95 ) = W( 95 ) + a*JVS( 635 )
  W( 96 ) = W( 96 ) + a*JVS( 636 )
  W( 97 ) = W( 97 ) + a*JVS( 637 )
  a = -W( 91 ) / JVS(          651  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 652 )
  W( 93 ) = W( 93 ) + a*JVS( 653 )
  W( 94 ) = W( 94 ) + a*JVS( 654 )
  W( 95 ) = W( 95 ) + a*JVS( 655 )
  W( 96 ) = W( 96 ) + a*JVS( 656 )
  W( 97 ) = W( 97 ) + a*JVS( 657 )
  JVS( 658) = W( 48 )
  JVS( 659) = W( 54 )
  JVS( 660) = W( 55 )
  JVS( 661) = W( 68 )
  JVS( 662) = W( 82 )
  JVS( 663) = W( 83 )
  JVS( 664) = W( 84 )
  JVS( 665) = W( 85 )
  JVS( 666) = W( 86 )
  JVS( 667) = W( 87 )
  JVS( 668) = W( 88 )
  JVS( 669) = W( 89 )
  JVS( 670) = W( 90 )
  JVS( 671) = W( 91 )
  JVS( 672) = W( 92 )
  JVS( 673) = W( 93 )
  JVS( 674) = W( 94 )
  JVS( 675) = W( 95 )
  JVS( 676) = W( 96 )
  JVS( 677) = W( 97 )
  IF ( ABS(  JVS( 696 )) < TINY(a) ) THEN
         IER = 93                                      
         RETURN
  END IF
   W( 45 ) = JVS( 678 )
   W( 60 ) = JVS( 679 )
   W( 74 ) = JVS( 680 )
   W( 75 ) = JVS( 681 )
   W( 76 ) = JVS( 682 )
   W( 77 ) = JVS( 683 )
   W( 80 ) = JVS( 684 )
   W( 82 ) = JVS( 685 )
   W( 83 ) = JVS( 686 )
   W( 84 ) = JVS( 687 )
   W( 85 ) = JVS( 688 )
   W( 86 ) = JVS( 689 )
   W( 87 ) = JVS( 690 )
   W( 88 ) = JVS( 691 )
   W( 89 ) = JVS( 692 )
   W( 90 ) = JVS( 693 )
   W( 91 ) = JVS( 694 )
   W( 92 ) = JVS( 695 )
   W( 93 ) = JVS( 696 )
   W( 94 ) = JVS( 697 )
   W( 95 ) = JVS( 698 )
   W( 96 ) = JVS( 699 )
   W( 97 ) = JVS( 700 )
  a = -W( 45 ) / JVS(          138  )
  W( 45 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 139 )
  W( 96 ) = W( 96 ) + a*JVS( 140 )
  a = -W( 60 ) / JVS(          199  )
  W( 60 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 200 )
  W( 90 ) = W( 90 ) + a*JVS( 201 )
  W( 93 ) = W( 93 ) + a*JVS( 202 )
  W( 96 ) = W( 96 ) + a*JVS( 203 )
  a = -W( 74 ) / JVS(          309  )
  W( 74 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 310 )
  W( 88 ) = W( 88 ) + a*JVS( 311 )
  W( 93 ) = W( 93 ) + a*JVS( 312 )
  W( 96 ) = W( 96 ) + a*JVS( 313 )
  W( 97 ) = W( 97 ) + a*JVS( 314 )
  a = -W( 75 ) / JVS(          315  )
  W( 75 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 316 )
  W( 93 ) = W( 93 ) + a*JVS( 317 )
  W( 96 ) = W( 96 ) + a*JVS( 318 )
  W( 97 ) = W( 97 ) + a*JVS( 319 )
  a = -W( 76 ) / JVS(          320  )
  W( 76 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 93 ) = W( 93 ) + a*JVS( 323 )
  W( 96 ) = W( 96 ) + a*JVS( 324 )
  W( 97 ) = W( 97 ) + a*JVS( 325 )
  a = -W( 77 ) / JVS(          327  )
  W( 77 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 328 )
  W( 88 ) = W( 88 ) + a*JVS( 329 )
  W( 93 ) = W( 93 ) + a*JVS( 330 )
  W( 96 ) = W( 96 ) + a*JVS( 331 )
  W( 97 ) = W( 97 ) + a*JVS( 332 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  a = -W( 82 ) / JVS(          392  )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 393 )
  W( 86 ) = W( 86 ) + a*JVS( 394 )
  W( 87 ) = W( 87 ) + a*JVS( 395 )
  W( 88 ) = W( 88 ) + a*JVS( 396 )
  W( 89 ) = W( 89 ) + a*JVS( 397 )
  W( 90 ) = W( 90 ) + a*JVS( 398 )
  W( 93 ) = W( 93 ) + a*JVS( 399 )
  W( 95 ) = W( 95 ) + a*JVS( 400 )
  W( 96 ) = W( 96 ) + a*JVS( 401 )
  W( 97 ) = W( 97 ) + a*JVS( 402 )
  a = -W( 83 ) / JVS(          414  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 415 )
  W( 86 ) = W( 86 ) + a*JVS( 416 )
  W( 87 ) = W( 87 ) + a*JVS( 417 )
  W( 88 ) = W( 88 ) + a*JVS( 418 )
  W( 89 ) = W( 89 ) + a*JVS( 419 )
  W( 90 ) = W( 90 ) + a*JVS( 420 )
  W( 91 ) = W( 91 ) + a*JVS( 421 )
  W( 92 ) = W( 92 ) + a*JVS( 422 )
  W( 93 ) = W( 93 ) + a*JVS( 423 )
  W( 94 ) = W( 94 ) + a*JVS( 424 )
  W( 95 ) = W( 95 ) + a*JVS( 425 )
  W( 96 ) = W( 96 ) + a*JVS( 426 )
  W( 97 ) = W( 97 ) + a*JVS( 427 )
  a = -W( 84 ) / JVS(          437  )
  W( 84 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 88 ) = W( 88 ) + a*JVS( 440 )
  W( 89 ) = W( 89 ) + a*JVS( 441 )
  W( 90 ) = W( 90 ) + a*JVS( 442 )
  W( 93 ) = W( 93 ) + a*JVS( 443 )
  W( 95 ) = W( 95 ) + a*JVS( 444 )
  W( 96 ) = W( 96 ) + a*JVS( 445 )
  W( 97 ) = W( 97 ) + a*JVS( 446 )
  a = -W( 85 ) / JVS(          460  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 461 )
  W( 87 ) = W( 87 ) + a*JVS( 462 )
  W( 88 ) = W( 88 ) + a*JVS( 463 )
  W( 89 ) = W( 89 ) + a*JVS( 464 )
  W( 90 ) = W( 90 ) + a*JVS( 465 )
  W( 91 ) = W( 91 ) + a*JVS( 466 )
  W( 92 ) = W( 92 ) + a*JVS( 467 )
  W( 93 ) = W( 93 ) + a*JVS( 468 )
  W( 94 ) = W( 94 ) + a*JVS( 469 )
  W( 95 ) = W( 95 ) + a*JVS( 470 )
  W( 96 ) = W( 96 ) + a*JVS( 471 )
  W( 97 ) = W( 97 ) + a*JVS( 472 )
  a = -W( 86 ) / JVS(          483  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 90 ) = W( 90 ) + a*JVS( 486 )
  W( 91 ) = W( 91 ) + a*JVS( 487 )
  W( 93 ) = W( 93 ) + a*JVS( 488 )
  W( 94 ) = W( 94 ) + a*JVS( 489 )
  W( 95 ) = W( 95 ) + a*JVS( 490 )
  W( 96 ) = W( 96 ) + a*JVS( 491 )
  W( 97 ) = W( 97 ) + a*JVS( 492 )
  a = -W( 87 ) / JVS(          517  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 518 )
  W( 89 ) = W( 89 ) + a*JVS( 519 )
  W( 90 ) = W( 90 ) + a*JVS( 520 )
  W( 91 ) = W( 91 ) + a*JVS( 521 )
  W( 92 ) = W( 92 ) + a*JVS( 522 )
  W( 93 ) = W( 93 ) + a*JVS( 523 )
  W( 94 ) = W( 94 ) + a*JVS( 524 )
  W( 95 ) = W( 95 ) + a*JVS( 525 )
  W( 96 ) = W( 96 ) + a*JVS( 526 )
  W( 97 ) = W( 97 ) + a*JVS( 527 )
  a = -W( 88 ) / JVS(          549  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 550 )
  W( 90 ) = W( 90 ) + a*JVS( 551 )
  W( 91 ) = W( 91 ) + a*JVS( 552 )
  W( 92 ) = W( 92 ) + a*JVS( 553 )
  W( 93 ) = W( 93 ) + a*JVS( 554 )
  W( 94 ) = W( 94 ) + a*JVS( 555 )
  W( 95 ) = W( 95 ) + a*JVS( 556 )
  W( 96 ) = W( 96 ) + a*JVS( 557 )
  W( 97 ) = W( 97 ) + a*JVS( 558 )
  a = -W( 89 ) / JVS(          582  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 583 )
  W( 91 ) = W( 91 ) + a*JVS( 584 )
  W( 92 ) = W( 92 ) + a*JVS( 585 )
  W( 93 ) = W( 93 ) + a*JVS( 586 )
  W( 94 ) = W( 94 ) + a*JVS( 587 )
  W( 95 ) = W( 95 ) + a*JVS( 588 )
  W( 96 ) = W( 96 ) + a*JVS( 589 )
  W( 97 ) = W( 97 ) + a*JVS( 590 )
  a = -W( 90 ) / JVS(          630  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 631 )
  W( 92 ) = W( 92 ) + a*JVS( 632 )
  W( 93 ) = W( 93 ) + a*JVS( 633 )
  W( 94 ) = W( 94 ) + a*JVS( 634 )
  W( 95 ) = W( 95 ) + a*JVS( 635 )
  W( 96 ) = W( 96 ) + a*JVS( 636 )
  W( 97 ) = W( 97 ) + a*JVS( 637 )
  a = -W( 91 ) / JVS(          651  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 652 )
  W( 93 ) = W( 93 ) + a*JVS( 653 )
  W( 94 ) = W( 94 ) + a*JVS( 654 )
  W( 95 ) = W( 95 ) + a*JVS( 655 )
  W( 96 ) = W( 96 ) + a*JVS( 656 )
  W( 97 ) = W( 97 ) + a*JVS( 657 )
  a = -W( 92 ) / JVS(          672  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 673 )
  W( 94 ) = W( 94 ) + a*JVS( 674 )
  W( 95 ) = W( 95 ) + a*JVS( 675 )
  W( 96 ) = W( 96 ) + a*JVS( 676 )
  W( 97 ) = W( 97 ) + a*JVS( 677 )
  JVS( 678) = W( 45 )
  JVS( 679) = W( 60 )
  JVS( 680) = W( 74 )
  JVS( 681) = W( 75 )
  JVS( 682) = W( 76 )
  JVS( 683) = W( 77 )
  JVS( 684) = W( 80 )
  JVS( 685) = W( 82 )
  JVS( 686) = W( 83 )
  JVS( 687) = W( 84 )
  JVS( 688) = W( 85 )
  JVS( 689) = W( 86 )
  JVS( 690) = W( 87 )
  JVS( 691) = W( 88 )
  JVS( 692) = W( 89 )
  JVS( 693) = W( 90 )
  JVS( 694) = W( 91 )
  JVS( 695) = W( 92 )
  JVS( 696) = W( 93 )
  JVS( 697) = W( 94 )
  JVS( 698) = W( 95 )
  JVS( 699) = W( 96 )
  JVS( 700) = W( 97 )
  IF ( ABS(  JVS( 715 )) < TINY(a) ) THEN
         IER = 94                                      
         RETURN
  END IF
   W( 56 ) = JVS( 701 )
   W( 75 ) = JVS( 702 )
   W( 80 ) = JVS( 703 )
   W( 81 ) = JVS( 704 )
   W( 82 ) = JVS( 705 )
   W( 84 ) = JVS( 706 )
   W( 86 ) = JVS( 707 )
   W( 87 ) = JVS( 708 )
   W( 88 ) = JVS( 709 )
   W( 89 ) = JVS( 710 )
   W( 90 ) = JVS( 711 )
   W( 91 ) = JVS( 712 )
   W( 92 ) = JVS( 713 )
   W( 93 ) = JVS( 714 )
   W( 94 ) = JVS( 715 )
   W( 95 ) = JVS( 716 )
   W( 96 ) = JVS( 717 )
   W( 97 ) = JVS( 718 )
  a = -W( 56 ) / JVS(          184  )
  W( 56 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 185 )
  W( 94 ) = W( 94 ) + a*JVS( 186 )
  W( 96 ) = W( 96 ) + a*JVS( 187 )
  a = -W( 75 ) / JVS(          315  )
  W( 75 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 316 )
  W( 93 ) = W( 93 ) + a*JVS( 317 )
  W( 96 ) = W( 96 ) + a*JVS( 318 )
  W( 97 ) = W( 97 ) + a*JVS( 319 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  a = -W( 81 ) / JVS(          372  )
  W( 81 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 373 )
  W( 87 ) = W( 87 ) + a*JVS( 374 )
  W( 88 ) = W( 88 ) + a*JVS( 375 )
  W( 93 ) = W( 93 ) + a*JVS( 376 )
  W( 96 ) = W( 96 ) + a*JVS( 377 )
  W( 97 ) = W( 97 ) + a*JVS( 378 )
  a = -W( 82 ) / JVS(          392  )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 393 )
  W( 86 ) = W( 86 ) + a*JVS( 394 )
  W( 87 ) = W( 87 ) + a*JVS( 395 )
  W( 88 ) = W( 88 ) + a*JVS( 396 )
  W( 89 ) = W( 89 ) + a*JVS( 397 )
  W( 90 ) = W( 90 ) + a*JVS( 398 )
  W( 93 ) = W( 93 ) + a*JVS( 399 )
  W( 95 ) = W( 95 ) + a*JVS( 400 )
  W( 96 ) = W( 96 ) + a*JVS( 401 )
  W( 97 ) = W( 97 ) + a*JVS( 402 )
  a = -W( 84 ) / JVS(          437  )
  W( 84 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 88 ) = W( 88 ) + a*JVS( 440 )
  W( 89 ) = W( 89 ) + a*JVS( 441 )
  W( 90 ) = W( 90 ) + a*JVS( 442 )
  W( 93 ) = W( 93 ) + a*JVS( 443 )
  W( 95 ) = W( 95 ) + a*JVS( 444 )
  W( 96 ) = W( 96 ) + a*JVS( 445 )
  W( 97 ) = W( 97 ) + a*JVS( 446 )
  a = -W( 86 ) / JVS(          483  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 90 ) = W( 90 ) + a*JVS( 486 )
  W( 91 ) = W( 91 ) + a*JVS( 487 )
  W( 93 ) = W( 93 ) + a*JVS( 488 )
  W( 94 ) = W( 94 ) + a*JVS( 489 )
  W( 95 ) = W( 95 ) + a*JVS( 490 )
  W( 96 ) = W( 96 ) + a*JVS( 491 )
  W( 97 ) = W( 97 ) + a*JVS( 492 )
  a = -W( 87 ) / JVS(          517  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 518 )
  W( 89 ) = W( 89 ) + a*JVS( 519 )
  W( 90 ) = W( 90 ) + a*JVS( 520 )
  W( 91 ) = W( 91 ) + a*JVS( 521 )
  W( 92 ) = W( 92 ) + a*JVS( 522 )
  W( 93 ) = W( 93 ) + a*JVS( 523 )
  W( 94 ) = W( 94 ) + a*JVS( 524 )
  W( 95 ) = W( 95 ) + a*JVS( 525 )
  W( 96 ) = W( 96 ) + a*JVS( 526 )
  W( 97 ) = W( 97 ) + a*JVS( 527 )
  a = -W( 88 ) / JVS(          549  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 550 )
  W( 90 ) = W( 90 ) + a*JVS( 551 )
  W( 91 ) = W( 91 ) + a*JVS( 552 )
  W( 92 ) = W( 92 ) + a*JVS( 553 )
  W( 93 ) = W( 93 ) + a*JVS( 554 )
  W( 94 ) = W( 94 ) + a*JVS( 555 )
  W( 95 ) = W( 95 ) + a*JVS( 556 )
  W( 96 ) = W( 96 ) + a*JVS( 557 )
  W( 97 ) = W( 97 ) + a*JVS( 558 )
  a = -W( 89 ) / JVS(          582  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 583 )
  W( 91 ) = W( 91 ) + a*JVS( 584 )
  W( 92 ) = W( 92 ) + a*JVS( 585 )
  W( 93 ) = W( 93 ) + a*JVS( 586 )
  W( 94 ) = W( 94 ) + a*JVS( 587 )
  W( 95 ) = W( 95 ) + a*JVS( 588 )
  W( 96 ) = W( 96 ) + a*JVS( 589 )
  W( 97 ) = W( 97 ) + a*JVS( 590 )
  a = -W( 90 ) / JVS(          630  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 631 )
  W( 92 ) = W( 92 ) + a*JVS( 632 )
  W( 93 ) = W( 93 ) + a*JVS( 633 )
  W( 94 ) = W( 94 ) + a*JVS( 634 )
  W( 95 ) = W( 95 ) + a*JVS( 635 )
  W( 96 ) = W( 96 ) + a*JVS( 636 )
  W( 97 ) = W( 97 ) + a*JVS( 637 )
  a = -W( 91 ) / JVS(          651  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 652 )
  W( 93 ) = W( 93 ) + a*JVS( 653 )
  W( 94 ) = W( 94 ) + a*JVS( 654 )
  W( 95 ) = W( 95 ) + a*JVS( 655 )
  W( 96 ) = W( 96 ) + a*JVS( 656 )
  W( 97 ) = W( 97 ) + a*JVS( 657 )
  a = -W( 92 ) / JVS(          672  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 673 )
  W( 94 ) = W( 94 ) + a*JVS( 674 )
  W( 95 ) = W( 95 ) + a*JVS( 675 )
  W( 96 ) = W( 96 ) + a*JVS( 676 )
  W( 97 ) = W( 97 ) + a*JVS( 677 )
  a = -W( 93 ) / JVS(          696  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 697 )
  W( 95 ) = W( 95 ) + a*JVS( 698 )
  W( 96 ) = W( 96 ) + a*JVS( 699 )
  W( 97 ) = W( 97 ) + a*JVS( 700 )
  JVS( 701) = W( 56 )
  JVS( 702) = W( 75 )
  JVS( 703) = W( 80 )
  JVS( 704) = W( 81 )
  JVS( 705) = W( 82 )
  JVS( 706) = W( 84 )
  JVS( 707) = W( 86 )
  JVS( 708) = W( 87 )
  JVS( 709) = W( 88 )
  JVS( 710) = W( 89 )
  JVS( 711) = W( 90 )
  JVS( 712) = W( 91 )
  JVS( 713) = W( 92 )
  JVS( 714) = W( 93 )
  JVS( 715) = W( 94 )
  JVS( 716) = W( 95 )
  JVS( 717) = W( 96 )
  JVS( 718) = W( 97 )
  IF ( ABS(  JVS( 738 )) < TINY(a) ) THEN
         IER = 95                                      
         RETURN
  END IF
   W( 50 ) = JVS( 719 )
   W( 52 ) = JVS( 720 )
   W( 62 ) = JVS( 721 )
   W( 64 ) = JVS( 722 )
   W( 78 ) = JVS( 723 )
   W( 79 ) = JVS( 724 )
   W( 80 ) = JVS( 725 )
   W( 81 ) = JVS( 726 )
   W( 84 ) = JVS( 727 )
   W( 85 ) = JVS( 728 )
   W( 86 ) = JVS( 729 )
   W( 87 ) = JVS( 730 )
   W( 88 ) = JVS( 731 )
   W( 89 ) = JVS( 732 )
   W( 90 ) = JVS( 733 )
   W( 91 ) = JVS( 734 )
   W( 92 ) = JVS( 735 )
   W( 93 ) = JVS( 736 )
   W( 94 ) = JVS( 737 )
   W( 95 ) = JVS( 738 )
   W( 96 ) = JVS( 739 )
   W( 97 ) = JVS( 740 )
  a = -W( 50 ) / JVS(          153  )
  W( 50 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 154 )
  W( 95 ) = W( 95 ) + a*JVS( 155 )
  W( 96 ) = W( 96 ) + a*JVS( 156 )
  a = -W( 52 ) / JVS(          165  )
  W( 52 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 166 )
  W( 96 ) = W( 96 ) + a*JVS( 167 )
  a = -W( 62 ) / JVS(          218  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 219 )
  W( 90 ) = W( 90 ) + a*JVS( 220 )
  W( 95 ) = W( 95 ) + a*JVS( 221 )
  a = -W( 64 ) / JVS(          228  )
  W( 64 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 229 )
  W( 88 ) = W( 88 ) + a*JVS( 230 )
  W( 90 ) = W( 90 ) + a*JVS( 231 )
  W( 95 ) = W( 95 ) + a*JVS( 232 )
  a = -W( 78 ) / JVS(          337  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 86 ) = W( 86 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 88 ) = W( 88 ) + a*JVS( 342 )
  W( 89 ) = W( 89 ) + a*JVS( 343 )
  W( 90 ) = W( 90 ) + a*JVS( 344 )
  W( 93 ) = W( 93 ) + a*JVS( 345 )
  W( 95 ) = W( 95 ) + a*JVS( 346 )
  W( 96 ) = W( 96 ) + a*JVS( 347 )
  W( 97 ) = W( 97 ) + a*JVS( 348 )
  a = -W( 79 ) / JVS(          354  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 355 )
  W( 81 ) = W( 81 ) + a*JVS( 356 )
  W( 84 ) = W( 84 ) + a*JVS( 357 )
  W( 86 ) = W( 86 ) + a*JVS( 358 )
  W( 87 ) = W( 87 ) + a*JVS( 359 )
  W( 88 ) = W( 88 ) + a*JVS( 360 )
  W( 93 ) = W( 93 ) + a*JVS( 361 )
  W( 96 ) = W( 96 ) + a*JVS( 362 )
  W( 97 ) = W( 97 ) + a*JVS( 363 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  a = -W( 81 ) / JVS(          372  )
  W( 81 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 373 )
  W( 87 ) = W( 87 ) + a*JVS( 374 )
  W( 88 ) = W( 88 ) + a*JVS( 375 )
  W( 93 ) = W( 93 ) + a*JVS( 376 )
  W( 96 ) = W( 96 ) + a*JVS( 377 )
  W( 97 ) = W( 97 ) + a*JVS( 378 )
  a = -W( 84 ) / JVS(          437  )
  W( 84 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 88 ) = W( 88 ) + a*JVS( 440 )
  W( 89 ) = W( 89 ) + a*JVS( 441 )
  W( 90 ) = W( 90 ) + a*JVS( 442 )
  W( 93 ) = W( 93 ) + a*JVS( 443 )
  W( 95 ) = W( 95 ) + a*JVS( 444 )
  W( 96 ) = W( 96 ) + a*JVS( 445 )
  W( 97 ) = W( 97 ) + a*JVS( 446 )
  a = -W( 85 ) / JVS(          460  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 461 )
  W( 87 ) = W( 87 ) + a*JVS( 462 )
  W( 88 ) = W( 88 ) + a*JVS( 463 )
  W( 89 ) = W( 89 ) + a*JVS( 464 )
  W( 90 ) = W( 90 ) + a*JVS( 465 )
  W( 91 ) = W( 91 ) + a*JVS( 466 )
  W( 92 ) = W( 92 ) + a*JVS( 467 )
  W( 93 ) = W( 93 ) + a*JVS( 468 )
  W( 94 ) = W( 94 ) + a*JVS( 469 )
  W( 95 ) = W( 95 ) + a*JVS( 470 )
  W( 96 ) = W( 96 ) + a*JVS( 471 )
  W( 97 ) = W( 97 ) + a*JVS( 472 )
  a = -W( 86 ) / JVS(          483  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 90 ) = W( 90 ) + a*JVS( 486 )
  W( 91 ) = W( 91 ) + a*JVS( 487 )
  W( 93 ) = W( 93 ) + a*JVS( 488 )
  W( 94 ) = W( 94 ) + a*JVS( 489 )
  W( 95 ) = W( 95 ) + a*JVS( 490 )
  W( 96 ) = W( 96 ) + a*JVS( 491 )
  W( 97 ) = W( 97 ) + a*JVS( 492 )
  a = -W( 87 ) / JVS(          517  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 518 )
  W( 89 ) = W( 89 ) + a*JVS( 519 )
  W( 90 ) = W( 90 ) + a*JVS( 520 )
  W( 91 ) = W( 91 ) + a*JVS( 521 )
  W( 92 ) = W( 92 ) + a*JVS( 522 )
  W( 93 ) = W( 93 ) + a*JVS( 523 )
  W( 94 ) = W( 94 ) + a*JVS( 524 )
  W( 95 ) = W( 95 ) + a*JVS( 525 )
  W( 96 ) = W( 96 ) + a*JVS( 526 )
  W( 97 ) = W( 97 ) + a*JVS( 527 )
  a = -W( 88 ) / JVS(          549  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 550 )
  W( 90 ) = W( 90 ) + a*JVS( 551 )
  W( 91 ) = W( 91 ) + a*JVS( 552 )
  W( 92 ) = W( 92 ) + a*JVS( 553 )
  W( 93 ) = W( 93 ) + a*JVS( 554 )
  W( 94 ) = W( 94 ) + a*JVS( 555 )
  W( 95 ) = W( 95 ) + a*JVS( 556 )
  W( 96 ) = W( 96 ) + a*JVS( 557 )
  W( 97 ) = W( 97 ) + a*JVS( 558 )
  a = -W( 89 ) / JVS(          582  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 583 )
  W( 91 ) = W( 91 ) + a*JVS( 584 )
  W( 92 ) = W( 92 ) + a*JVS( 585 )
  W( 93 ) = W( 93 ) + a*JVS( 586 )
  W( 94 ) = W( 94 ) + a*JVS( 587 )
  W( 95 ) = W( 95 ) + a*JVS( 588 )
  W( 96 ) = W( 96 ) + a*JVS( 589 )
  W( 97 ) = W( 97 ) + a*JVS( 590 )
  a = -W( 90 ) / JVS(          630  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 631 )
  W( 92 ) = W( 92 ) + a*JVS( 632 )
  W( 93 ) = W( 93 ) + a*JVS( 633 )
  W( 94 ) = W( 94 ) + a*JVS( 634 )
  W( 95 ) = W( 95 ) + a*JVS( 635 )
  W( 96 ) = W( 96 ) + a*JVS( 636 )
  W( 97 ) = W( 97 ) + a*JVS( 637 )
  a = -W( 91 ) / JVS(          651  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 652 )
  W( 93 ) = W( 93 ) + a*JVS( 653 )
  W( 94 ) = W( 94 ) + a*JVS( 654 )
  W( 95 ) = W( 95 ) + a*JVS( 655 )
  W( 96 ) = W( 96 ) + a*JVS( 656 )
  W( 97 ) = W( 97 ) + a*JVS( 657 )
  a = -W( 92 ) / JVS(          672  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 673 )
  W( 94 ) = W( 94 ) + a*JVS( 674 )
  W( 95 ) = W( 95 ) + a*JVS( 675 )
  W( 96 ) = W( 96 ) + a*JVS( 676 )
  W( 97 ) = W( 97 ) + a*JVS( 677 )
  a = -W( 93 ) / JVS(          696  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 697 )
  W( 95 ) = W( 95 ) + a*JVS( 698 )
  W( 96 ) = W( 96 ) + a*JVS( 699 )
  W( 97 ) = W( 97 ) + a*JVS( 700 )
  a = -W( 94 ) / JVS(          715  )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  W( 97 ) = W( 97 ) + a*JVS( 718 )
  JVS( 719) = W( 50 )
  JVS( 720) = W( 52 )
  JVS( 721) = W( 62 )
  JVS( 722) = W( 64 )
  JVS( 723) = W( 78 )
  JVS( 724) = W( 79 )
  JVS( 725) = W( 80 )
  JVS( 726) = W( 81 )
  JVS( 727) = W( 84 )
  JVS( 728) = W( 85 )
  JVS( 729) = W( 86 )
  JVS( 730) = W( 87 )
  JVS( 731) = W( 88 )
  JVS( 732) = W( 89 )
  JVS( 733) = W( 90 )
  JVS( 734) = W( 91 )
  JVS( 735) = W( 92 )
  JVS( 736) = W( 93 )
  JVS( 737) = W( 94 )
  JVS( 738) = W( 95 )
  JVS( 739) = W( 96 )
  JVS( 740) = W( 97 )
  IF ( ABS(  JVS( 792 )) < TINY(a) ) THEN
         IER = 96                                      
         RETURN
  END IF
   W( 40 ) = JVS( 741 )
   W( 42 ) = JVS( 742 )
   W( 43 ) = JVS( 743 )
   W( 44 ) = JVS( 744 )
   W( 45 ) = JVS( 745 )
   W( 47 ) = JVS( 746 )
   W( 48 ) = JVS( 747 )
   W( 49 ) = JVS( 748 )
   W( 50 ) = JVS( 749 )
   W( 51 ) = JVS( 750 )
   W( 53 ) = JVS( 751 )
   W( 54 ) = JVS( 752 )
   W( 55 ) = JVS( 753 )
   W( 56 ) = JVS( 754 )
   W( 57 ) = JVS( 755 )
   W( 58 ) = JVS( 756 )
   W( 59 ) = JVS( 757 )
   W( 60 ) = JVS( 758 )
   W( 61 ) = JVS( 759 )
   W( 62 ) = JVS( 760 )
   W( 64 ) = JVS( 761 )
   W( 65 ) = JVS( 762 )
   W( 66 ) = JVS( 763 )
   W( 67 ) = JVS( 764 )
   W( 68 ) = JVS( 765 )
   W( 69 ) = JVS( 766 )
   W( 70 ) = JVS( 767 )
   W( 71 ) = JVS( 768 )
   W( 73 ) = JVS( 769 )
   W( 74 ) = JVS( 770 )
   W( 75 ) = JVS( 771 )
   W( 76 ) = JVS( 772 )
   W( 77 ) = JVS( 773 )
   W( 78 ) = JVS( 774 )
   W( 79 ) = JVS( 775 )
   W( 80 ) = JVS( 776 )
   W( 81 ) = JVS( 777 )
   W( 82 ) = JVS( 778 )
   W( 83 ) = JVS( 779 )
   W( 84 ) = JVS( 780 )
   W( 85 ) = JVS( 781 )
   W( 86 ) = JVS( 782 )
   W( 87 ) = JVS( 783 )
   W( 88 ) = JVS( 784 )
   W( 89 ) = JVS( 785 )
   W( 90 ) = JVS( 786 )
   W( 91 ) = JVS( 787 )
   W( 92 ) = JVS( 788 )
   W( 93 ) = JVS( 789 )
   W( 94 ) = JVS( 790 )
   W( 95 ) = JVS( 791 )
   W( 96 ) = JVS( 792 )
   W( 97 ) = JVS( 793 )
  a = -W( 40 ) / JVS(          124  )
  W( 40 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 125 )
  a = -W( 42 ) / JVS(          129  )
  W( 42 ) = -a
  W( 64 ) = W( 64 ) + a*JVS( 130 )
  W( 90 ) = W( 90 ) + a*JVS( 131 )
  a = -W( 43 ) / JVS(          132  )
  W( 43 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 133 )
  a = -W( 44 ) / JVS(          134  )
  W( 44 ) = -a
  W( 45 ) = W( 45 ) + a*JVS( 135 )
  W( 96 ) = W( 96 ) + a*JVS( 136 )
  a = -W( 45 ) / JVS(          138  )
  W( 45 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 139 )
  W( 96 ) = W( 96 ) + a*JVS( 140 )
  a = -W( 47 ) / JVS(          144  )
  W( 47 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 145 )
  a = -W( 48 ) / JVS(          146  )
  W( 48 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 147 )
  W( 96 ) = W( 96 ) + a*JVS( 148 )
  a = -W( 49 ) / JVS(          149  )
  W( 49 ) = -a
  W( 60 ) = W( 60 ) + a*JVS( 150 )
  W( 86 ) = W( 86 ) + a*JVS( 151 )
  W( 96 ) = W( 96 ) + a*JVS( 152 )
  a = -W( 50 ) / JVS(          153  )
  W( 50 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 154 )
  W( 95 ) = W( 95 ) + a*JVS( 155 )
  W( 96 ) = W( 96 ) + a*JVS( 156 )
  a = -W( 51 ) / JVS(          157  )
  W( 51 ) = -a
  W( 62 ) = W( 62 ) + a*JVS( 158 )
  W( 74 ) = W( 74 ) + a*JVS( 159 )
  W( 86 ) = W( 86 ) + a*JVS( 160 )
  W( 95 ) = W( 95 ) + a*JVS( 161 )
  W( 96 ) = W( 96 ) + a*JVS( 162 )
  a = -W( 53 ) / JVS(          168  )
  W( 53 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 169 )
  W( 90 ) = W( 90 ) + a*JVS( 170 )
  W( 96 ) = W( 96 ) + a*JVS( 171 )
  a = -W( 54 ) / JVS(          172  )
  W( 54 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 173 )
  W( 91 ) = W( 91 ) + a*JVS( 174 )
  W( 94 ) = W( 94 ) + a*JVS( 175 )
  W( 96 ) = W( 96 ) + a*JVS( 176 )
  a = -W( 55 ) / JVS(          177  )
  W( 55 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 178 )
  W( 90 ) = W( 90 ) + a*JVS( 179 )
  W( 91 ) = W( 91 ) + a*JVS( 180 )
  W( 92 ) = W( 92 ) + a*JVS( 181 )
  W( 94 ) = W( 94 ) + a*JVS( 182 )
  W( 96 ) = W( 96 ) + a*JVS( 183 )
  a = -W( 56 ) / JVS(          184  )
  W( 56 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 185 )
  W( 94 ) = W( 94 ) + a*JVS( 186 )
  W( 96 ) = W( 96 ) + a*JVS( 187 )
  a = -W( 57 ) / JVS(          188  )
  W( 57 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 189 )
  W( 96 ) = W( 96 ) + a*JVS( 190 )
  a = -W( 58 ) / JVS(          191  )
  W( 58 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 192 )
  W( 92 ) = W( 92 ) + a*JVS( 193 )
  W( 96 ) = W( 96 ) + a*JVS( 194 )
  a = -W( 59 ) / JVS(          195  )
  W( 59 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 196 )
  W( 96 ) = W( 96 ) + a*JVS( 197 )
  a = -W( 60 ) / JVS(          199  )
  W( 60 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 200 )
  W( 90 ) = W( 90 ) + a*JVS( 201 )
  W( 93 ) = W( 93 ) + a*JVS( 202 )
  W( 96 ) = W( 96 ) + a*JVS( 203 )
  a = -W( 61 ) / JVS(          208  )
  W( 61 ) = -a
  W( 76 ) = W( 76 ) + a*JVS( 209 )
  W( 79 ) = W( 79 ) + a*JVS( 210 )
  W( 80 ) = W( 80 ) + a*JVS( 211 )
  W( 82 ) = W( 82 ) + a*JVS( 212 )
  W( 83 ) = W( 83 ) + a*JVS( 213 )
  W( 85 ) = W( 85 ) + a*JVS( 214 )
  W( 88 ) = W( 88 ) + a*JVS( 215 )
  W( 92 ) = W( 92 ) + a*JVS( 216 )
  W( 96 ) = W( 96 ) + a*JVS( 217 )
  a = -W( 62 ) / JVS(          218  )
  W( 62 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 219 )
  W( 90 ) = W( 90 ) + a*JVS( 220 )
  W( 95 ) = W( 95 ) + a*JVS( 221 )
  a = -W( 64 ) / JVS(          228  )
  W( 64 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 229 )
  W( 88 ) = W( 88 ) + a*JVS( 230 )
  W( 90 ) = W( 90 ) + a*JVS( 231 )
  W( 95 ) = W( 95 ) + a*JVS( 232 )
  a = -W( 65 ) / JVS(          233  )
  W( 65 ) = -a
  W( 78 ) = W( 78 ) + a*JVS( 234 )
  W( 89 ) = W( 89 ) + a*JVS( 235 )
  W( 90 ) = W( 90 ) + a*JVS( 236 )
  W( 96 ) = W( 96 ) + a*JVS( 237 )
  a = -W( 66 ) / JVS(          239  )
  W( 66 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 240 )
  W( 81 ) = W( 81 ) + a*JVS( 241 )
  W( 86 ) = W( 86 ) + a*JVS( 242 )
  W( 96 ) = W( 96 ) + a*JVS( 243 )
  a = -W( 67 ) / JVS(          244  )
  W( 67 ) = -a
  W( 74 ) = W( 74 ) + a*JVS( 245 )
  W( 76 ) = W( 76 ) + a*JVS( 246 )
  W( 77 ) = W( 77 ) + a*JVS( 247 )
  W( 80 ) = W( 80 ) + a*JVS( 248 )
  W( 88 ) = W( 88 ) + a*JVS( 249 )
  W( 96 ) = W( 96 ) + a*JVS( 250 )
  a = -W( 68 ) / JVS(          252  )
  W( 68 ) = -a
  W( 85 ) = W( 85 ) + a*JVS( 253 )
  W( 90 ) = W( 90 ) + a*JVS( 254 )
  W( 92 ) = W( 92 ) + a*JVS( 255 )
  W( 95 ) = W( 95 ) + a*JVS( 256 )
  W( 96 ) = W( 96 ) + a*JVS( 257 )
  a = -W( 69 ) / JVS(          259  )
  W( 69 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 260 )
  W( 81 ) = W( 81 ) + a*JVS( 261 )
  W( 82 ) = W( 82 ) + a*JVS( 262 )
  W( 83 ) = W( 83 ) + a*JVS( 263 )
  W( 84 ) = W( 84 ) + a*JVS( 264 )
  W( 85 ) = W( 85 ) + a*JVS( 265 )
  W( 87 ) = W( 87 ) + a*JVS( 266 )
  W( 90 ) = W( 90 ) + a*JVS( 267 )
  W( 96 ) = W( 96 ) + a*JVS( 268 )
  W( 97 ) = W( 97 ) + a*JVS( 269 )
  a = -W( 70 ) / JVS(          272  )
  W( 70 ) = -a
  W( 71 ) = W( 71 ) + a*JVS( 273 )
  W( 74 ) = W( 74 ) + a*JVS( 274 )
  W( 75 ) = W( 75 ) + a*JVS( 275 )
  W( 76 ) = W( 76 ) + a*JVS( 276 )
  W( 77 ) = W( 77 ) + a*JVS( 277 )
  W( 80 ) = W( 80 ) + a*JVS( 278 )
  W( 81 ) = W( 81 ) + a*JVS( 279 )
  W( 82 ) = W( 82 ) + a*JVS( 280 )
  W( 83 ) = W( 83 ) + a*JVS( 281 )
  W( 85 ) = W( 85 ) + a*JVS( 282 )
  W( 86 ) = W( 86 ) + a*JVS( 283 )
  W( 88 ) = W( 88 ) + a*JVS( 284 )
  W( 93 ) = W( 93 ) + a*JVS( 285 )
  W( 96 ) = W( 96 ) + a*JVS( 286 )
  W( 97 ) = W( 97 ) + a*JVS( 287 )
  a = -W( 71 ) / JVS(          289  )
  W( 71 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 290 )
  W( 86 ) = W( 86 ) + a*JVS( 291 )
  W( 95 ) = W( 95 ) + a*JVS( 292 )
  W( 96 ) = W( 96 ) + a*JVS( 293 )
  a = -W( 73 ) / JVS(          303  )
  W( 73 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 304 )
  W( 90 ) = W( 90 ) + a*JVS( 305 )
  W( 95 ) = W( 95 ) + a*JVS( 306 )
  W( 96 ) = W( 96 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  a = -W( 74 ) / JVS(          309  )
  W( 74 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 310 )
  W( 88 ) = W( 88 ) + a*JVS( 311 )
  W( 93 ) = W( 93 ) + a*JVS( 312 )
  W( 96 ) = W( 96 ) + a*JVS( 313 )
  W( 97 ) = W( 97 ) + a*JVS( 314 )
  a = -W( 75 ) / JVS(          315  )
  W( 75 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 316 )
  W( 93 ) = W( 93 ) + a*JVS( 317 )
  W( 96 ) = W( 96 ) + a*JVS( 318 )
  W( 97 ) = W( 97 ) + a*JVS( 319 )
  a = -W( 76 ) / JVS(          320  )
  W( 76 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 93 ) = W( 93 ) + a*JVS( 323 )
  W( 96 ) = W( 96 ) + a*JVS( 324 )
  W( 97 ) = W( 97 ) + a*JVS( 325 )
  a = -W( 77 ) / JVS(          327  )
  W( 77 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 328 )
  W( 88 ) = W( 88 ) + a*JVS( 329 )
  W( 93 ) = W( 93 ) + a*JVS( 330 )
  W( 96 ) = W( 96 ) + a*JVS( 331 )
  W( 97 ) = W( 97 ) + a*JVS( 332 )
  a = -W( 78 ) / JVS(          337  )
  W( 78 ) = -a
  W( 79 ) = W( 79 ) + a*JVS( 338 )
  W( 80 ) = W( 80 ) + a*JVS( 339 )
  W( 86 ) = W( 86 ) + a*JVS( 340 )
  W( 87 ) = W( 87 ) + a*JVS( 341 )
  W( 88 ) = W( 88 ) + a*JVS( 342 )
  W( 89 ) = W( 89 ) + a*JVS( 343 )
  W( 90 ) = W( 90 ) + a*JVS( 344 )
  W( 93 ) = W( 93 ) + a*JVS( 345 )
  W( 95 ) = W( 95 ) + a*JVS( 346 )
  W( 96 ) = W( 96 ) + a*JVS( 347 )
  W( 97 ) = W( 97 ) + a*JVS( 348 )
  a = -W( 79 ) / JVS(          354  )
  W( 79 ) = -a
  W( 80 ) = W( 80 ) + a*JVS( 355 )
  W( 81 ) = W( 81 ) + a*JVS( 356 )
  W( 84 ) = W( 84 ) + a*JVS( 357 )
  W( 86 ) = W( 86 ) + a*JVS( 358 )
  W( 87 ) = W( 87 ) + a*JVS( 359 )
  W( 88 ) = W( 88 ) + a*JVS( 360 )
  W( 93 ) = W( 93 ) + a*JVS( 361 )
  W( 96 ) = W( 96 ) + a*JVS( 362 )
  W( 97 ) = W( 97 ) + a*JVS( 363 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  a = -W( 81 ) / JVS(          372  )
  W( 81 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 373 )
  W( 87 ) = W( 87 ) + a*JVS( 374 )
  W( 88 ) = W( 88 ) + a*JVS( 375 )
  W( 93 ) = W( 93 ) + a*JVS( 376 )
  W( 96 ) = W( 96 ) + a*JVS( 377 )
  W( 97 ) = W( 97 ) + a*JVS( 378 )
  a = -W( 82 ) / JVS(          392  )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 393 )
  W( 86 ) = W( 86 ) + a*JVS( 394 )
  W( 87 ) = W( 87 ) + a*JVS( 395 )
  W( 88 ) = W( 88 ) + a*JVS( 396 )
  W( 89 ) = W( 89 ) + a*JVS( 397 )
  W( 90 ) = W( 90 ) + a*JVS( 398 )
  W( 93 ) = W( 93 ) + a*JVS( 399 )
  W( 95 ) = W( 95 ) + a*JVS( 400 )
  W( 96 ) = W( 96 ) + a*JVS( 401 )
  W( 97 ) = W( 97 ) + a*JVS( 402 )
  a = -W( 83 ) / JVS(          414  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 415 )
  W( 86 ) = W( 86 ) + a*JVS( 416 )
  W( 87 ) = W( 87 ) + a*JVS( 417 )
  W( 88 ) = W( 88 ) + a*JVS( 418 )
  W( 89 ) = W( 89 ) + a*JVS( 419 )
  W( 90 ) = W( 90 ) + a*JVS( 420 )
  W( 91 ) = W( 91 ) + a*JVS( 421 )
  W( 92 ) = W( 92 ) + a*JVS( 422 )
  W( 93 ) = W( 93 ) + a*JVS( 423 )
  W( 94 ) = W( 94 ) + a*JVS( 424 )
  W( 95 ) = W( 95 ) + a*JVS( 425 )
  W( 96 ) = W( 96 ) + a*JVS( 426 )
  W( 97 ) = W( 97 ) + a*JVS( 427 )
  a = -W( 84 ) / JVS(          437  )
  W( 84 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 88 ) = W( 88 ) + a*JVS( 440 )
  W( 89 ) = W( 89 ) + a*JVS( 441 )
  W( 90 ) = W( 90 ) + a*JVS( 442 )
  W( 93 ) = W( 93 ) + a*JVS( 443 )
  W( 95 ) = W( 95 ) + a*JVS( 444 )
  W( 96 ) = W( 96 ) + a*JVS( 445 )
  W( 97 ) = W( 97 ) + a*JVS( 446 )
  a = -W( 85 ) / JVS(          460  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 461 )
  W( 87 ) = W( 87 ) + a*JVS( 462 )
  W( 88 ) = W( 88 ) + a*JVS( 463 )
  W( 89 ) = W( 89 ) + a*JVS( 464 )
  W( 90 ) = W( 90 ) + a*JVS( 465 )
  W( 91 ) = W( 91 ) + a*JVS( 466 )
  W( 92 ) = W( 92 ) + a*JVS( 467 )
  W( 93 ) = W( 93 ) + a*JVS( 468 )
  W( 94 ) = W( 94 ) + a*JVS( 469 )
  W( 95 ) = W( 95 ) + a*JVS( 470 )
  W( 96 ) = W( 96 ) + a*JVS( 471 )
  W( 97 ) = W( 97 ) + a*JVS( 472 )
  a = -W( 86 ) / JVS(          483  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 90 ) = W( 90 ) + a*JVS( 486 )
  W( 91 ) = W( 91 ) + a*JVS( 487 )
  W( 93 ) = W( 93 ) + a*JVS( 488 )
  W( 94 ) = W( 94 ) + a*JVS( 489 )
  W( 95 ) = W( 95 ) + a*JVS( 490 )
  W( 96 ) = W( 96 ) + a*JVS( 491 )
  W( 97 ) = W( 97 ) + a*JVS( 492 )
  a = -W( 87 ) / JVS(          517  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 518 )
  W( 89 ) = W( 89 ) + a*JVS( 519 )
  W( 90 ) = W( 90 ) + a*JVS( 520 )
  W( 91 ) = W( 91 ) + a*JVS( 521 )
  W( 92 ) = W( 92 ) + a*JVS( 522 )
  W( 93 ) = W( 93 ) + a*JVS( 523 )
  W( 94 ) = W( 94 ) + a*JVS( 524 )
  W( 95 ) = W( 95 ) + a*JVS( 525 )
  W( 96 ) = W( 96 ) + a*JVS( 526 )
  W( 97 ) = W( 97 ) + a*JVS( 527 )
  a = -W( 88 ) / JVS(          549  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 550 )
  W( 90 ) = W( 90 ) + a*JVS( 551 )
  W( 91 ) = W( 91 ) + a*JVS( 552 )
  W( 92 ) = W( 92 ) + a*JVS( 553 )
  W( 93 ) = W( 93 ) + a*JVS( 554 )
  W( 94 ) = W( 94 ) + a*JVS( 555 )
  W( 95 ) = W( 95 ) + a*JVS( 556 )
  W( 96 ) = W( 96 ) + a*JVS( 557 )
  W( 97 ) = W( 97 ) + a*JVS( 558 )
  a = -W( 89 ) / JVS(          582  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 583 )
  W( 91 ) = W( 91 ) + a*JVS( 584 )
  W( 92 ) = W( 92 ) + a*JVS( 585 )
  W( 93 ) = W( 93 ) + a*JVS( 586 )
  W( 94 ) = W( 94 ) + a*JVS( 587 )
  W( 95 ) = W( 95 ) + a*JVS( 588 )
  W( 96 ) = W( 96 ) + a*JVS( 589 )
  W( 97 ) = W( 97 ) + a*JVS( 590 )
  a = -W( 90 ) / JVS(          630  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 631 )
  W( 92 ) = W( 92 ) + a*JVS( 632 )
  W( 93 ) = W( 93 ) + a*JVS( 633 )
  W( 94 ) = W( 94 ) + a*JVS( 634 )
  W( 95 ) = W( 95 ) + a*JVS( 635 )
  W( 96 ) = W( 96 ) + a*JVS( 636 )
  W( 97 ) = W( 97 ) + a*JVS( 637 )
  a = -W( 91 ) / JVS(          651  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 652 )
  W( 93 ) = W( 93 ) + a*JVS( 653 )
  W( 94 ) = W( 94 ) + a*JVS( 654 )
  W( 95 ) = W( 95 ) + a*JVS( 655 )
  W( 96 ) = W( 96 ) + a*JVS( 656 )
  W( 97 ) = W( 97 ) + a*JVS( 657 )
  a = -W( 92 ) / JVS(          672  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 673 )
  W( 94 ) = W( 94 ) + a*JVS( 674 )
  W( 95 ) = W( 95 ) + a*JVS( 675 )
  W( 96 ) = W( 96 ) + a*JVS( 676 )
  W( 97 ) = W( 97 ) + a*JVS( 677 )
  a = -W( 93 ) / JVS(          696  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 697 )
  W( 95 ) = W( 95 ) + a*JVS( 698 )
  W( 96 ) = W( 96 ) + a*JVS( 699 )
  W( 97 ) = W( 97 ) + a*JVS( 700 )
  a = -W( 94 ) / JVS(          715  )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  W( 97 ) = W( 97 ) + a*JVS( 718 )
  a = -W( 95 ) / JVS(          738  )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 739 )
  W( 97 ) = W( 97 ) + a*JVS( 740 )
  JVS( 741) = W( 40 )
  JVS( 742) = W( 42 )
  JVS( 743) = W( 43 )
  JVS( 744) = W( 44 )
  JVS( 745) = W( 45 )
  JVS( 746) = W( 47 )
  JVS( 747) = W( 48 )
  JVS( 748) = W( 49 )
  JVS( 749) = W( 50 )
  JVS( 750) = W( 51 )
  JVS( 751) = W( 53 )
  JVS( 752) = W( 54 )
  JVS( 753) = W( 55 )
  JVS( 754) = W( 56 )
  JVS( 755) = W( 57 )
  JVS( 756) = W( 58 )
  JVS( 757) = W( 59 )
  JVS( 758) = W( 60 )
  JVS( 759) = W( 61 )
  JVS( 760) = W( 62 )
  JVS( 761) = W( 64 )
  JVS( 762) = W( 65 )
  JVS( 763) = W( 66 )
  JVS( 764) = W( 67 )
  JVS( 765) = W( 68 )
  JVS( 766) = W( 69 )
  JVS( 767) = W( 70 )
  JVS( 768) = W( 71 )
  JVS( 769) = W( 73 )
  JVS( 770) = W( 74 )
  JVS( 771) = W( 75 )
  JVS( 772) = W( 76 )
  JVS( 773) = W( 77 )
  JVS( 774) = W( 78 )
  JVS( 775) = W( 79 )
  JVS( 776) = W( 80 )
  JVS( 777) = W( 81 )
  JVS( 778) = W( 82 )
  JVS( 779) = W( 83 )
  JVS( 780) = W( 84 )
  JVS( 781) = W( 85 )
  JVS( 782) = W( 86 )
  JVS( 783) = W( 87 )
  JVS( 784) = W( 88 )
  JVS( 785) = W( 89 )
  JVS( 786) = W( 90 )
  JVS( 787) = W( 91 )
  JVS( 788) = W( 92 )
  JVS( 789) = W( 93 )
  JVS( 790) = W( 94 )
  JVS( 791) = W( 95 )
  JVS( 792) = W( 96 )
  JVS( 793) = W( 97 )
  IF ( ABS(  JVS( 819 )) < TINY(a) ) THEN
         IER = 97                                      
         RETURN
  END IF
   W( 46 ) = JVS( 794 )
   W( 53 ) = JVS( 795 )
   W( 69 ) = JVS( 796 )
   W( 73 ) = JVS( 797 )
   W( 74 ) = JVS( 798 )
   W( 75 ) = JVS( 799 )
   W( 76 ) = JVS( 800 )
   W( 77 ) = JVS( 801 )
   W( 80 ) = JVS( 802 )
   W( 81 ) = JVS( 803 )
   W( 82 ) = JVS( 804 )
   W( 83 ) = JVS( 805 )
   W( 84 ) = JVS( 806 )
   W( 85 ) = JVS( 807 )
   W( 86 ) = JVS( 808 )
   W( 87 ) = JVS( 809 )
   W( 88 ) = JVS( 810 )
   W( 89 ) = JVS( 811 )
   W( 90 ) = JVS( 812 )
   W( 91 ) = JVS( 813 )
   W( 92 ) = JVS( 814 )
   W( 93 ) = JVS( 815 )
   W( 94 ) = JVS( 816 )
   W( 95 ) = JVS( 817 )
   W( 96 ) = JVS( 818 )
   W( 97 ) = JVS( 819 )
  a = -W( 46 ) / JVS(          141  )
  W( 46 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 142 )
  W( 97 ) = W( 97 ) + a*JVS( 143 )
  a = -W( 53 ) / JVS(          168  )
  W( 53 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 169 )
  W( 90 ) = W( 90 ) + a*JVS( 170 )
  W( 96 ) = W( 96 ) + a*JVS( 171 )
  a = -W( 69 ) / JVS(          259  )
  W( 69 ) = -a
  W( 73 ) = W( 73 ) + a*JVS( 260 )
  W( 81 ) = W( 81 ) + a*JVS( 261 )
  W( 82 ) = W( 82 ) + a*JVS( 262 )
  W( 83 ) = W( 83 ) + a*JVS( 263 )
  W( 84 ) = W( 84 ) + a*JVS( 264 )
  W( 85 ) = W( 85 ) + a*JVS( 265 )
  W( 87 ) = W( 87 ) + a*JVS( 266 )
  W( 90 ) = W( 90 ) + a*JVS( 267 )
  W( 96 ) = W( 96 ) + a*JVS( 268 )
  W( 97 ) = W( 97 ) + a*JVS( 269 )
  a = -W( 73 ) / JVS(          303  )
  W( 73 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 304 )
  W( 90 ) = W( 90 ) + a*JVS( 305 )
  W( 95 ) = W( 95 ) + a*JVS( 306 )
  W( 96 ) = W( 96 ) + a*JVS( 307 )
  W( 97 ) = W( 97 ) + a*JVS( 308 )
  a = -W( 74 ) / JVS(          309  )
  W( 74 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 310 )
  W( 88 ) = W( 88 ) + a*JVS( 311 )
  W( 93 ) = W( 93 ) + a*JVS( 312 )
  W( 96 ) = W( 96 ) + a*JVS( 313 )
  W( 97 ) = W( 97 ) + a*JVS( 314 )
  a = -W( 75 ) / JVS(          315  )
  W( 75 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 316 )
  W( 93 ) = W( 93 ) + a*JVS( 317 )
  W( 96 ) = W( 96 ) + a*JVS( 318 )
  W( 97 ) = W( 97 ) + a*JVS( 319 )
  a = -W( 76 ) / JVS(          320  )
  W( 76 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 321 )
  W( 88 ) = W( 88 ) + a*JVS( 322 )
  W( 93 ) = W( 93 ) + a*JVS( 323 )
  W( 96 ) = W( 96 ) + a*JVS( 324 )
  W( 97 ) = W( 97 ) + a*JVS( 325 )
  a = -W( 77 ) / JVS(          327  )
  W( 77 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 328 )
  W( 88 ) = W( 88 ) + a*JVS( 329 )
  W( 93 ) = W( 93 ) + a*JVS( 330 )
  W( 96 ) = W( 96 ) + a*JVS( 331 )
  W( 97 ) = W( 97 ) + a*JVS( 332 )
  a = -W( 80 ) / JVS(          364  )
  W( 80 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 365 )
  W( 87 ) = W( 87 ) + a*JVS( 366 )
  W( 88 ) = W( 88 ) + a*JVS( 367 )
  W( 93 ) = W( 93 ) + a*JVS( 368 )
  W( 96 ) = W( 96 ) + a*JVS( 369 )
  W( 97 ) = W( 97 ) + a*JVS( 370 )
  a = -W( 81 ) / JVS(          372  )
  W( 81 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 373 )
  W( 87 ) = W( 87 ) + a*JVS( 374 )
  W( 88 ) = W( 88 ) + a*JVS( 375 )
  W( 93 ) = W( 93 ) + a*JVS( 376 )
  W( 96 ) = W( 96 ) + a*JVS( 377 )
  W( 97 ) = W( 97 ) + a*JVS( 378 )
  a = -W( 82 ) / JVS(          392  )
  W( 82 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 393 )
  W( 86 ) = W( 86 ) + a*JVS( 394 )
  W( 87 ) = W( 87 ) + a*JVS( 395 )
  W( 88 ) = W( 88 ) + a*JVS( 396 )
  W( 89 ) = W( 89 ) + a*JVS( 397 )
  W( 90 ) = W( 90 ) + a*JVS( 398 )
  W( 93 ) = W( 93 ) + a*JVS( 399 )
  W( 95 ) = W( 95 ) + a*JVS( 400 )
  W( 96 ) = W( 96 ) + a*JVS( 401 )
  W( 97 ) = W( 97 ) + a*JVS( 402 )
  a = -W( 83 ) / JVS(          414  )
  W( 83 ) = -a
  W( 84 ) = W( 84 ) + a*JVS( 415 )
  W( 86 ) = W( 86 ) + a*JVS( 416 )
  W( 87 ) = W( 87 ) + a*JVS( 417 )
  W( 88 ) = W( 88 ) + a*JVS( 418 )
  W( 89 ) = W( 89 ) + a*JVS( 419 )
  W( 90 ) = W( 90 ) + a*JVS( 420 )
  W( 91 ) = W( 91 ) + a*JVS( 421 )
  W( 92 ) = W( 92 ) + a*JVS( 422 )
  W( 93 ) = W( 93 ) + a*JVS( 423 )
  W( 94 ) = W( 94 ) + a*JVS( 424 )
  W( 95 ) = W( 95 ) + a*JVS( 425 )
  W( 96 ) = W( 96 ) + a*JVS( 426 )
  W( 97 ) = W( 97 ) + a*JVS( 427 )
  a = -W( 84 ) / JVS(          437  )
  W( 84 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 438 )
  W( 87 ) = W( 87 ) + a*JVS( 439 )
  W( 88 ) = W( 88 ) + a*JVS( 440 )
  W( 89 ) = W( 89 ) + a*JVS( 441 )
  W( 90 ) = W( 90 ) + a*JVS( 442 )
  W( 93 ) = W( 93 ) + a*JVS( 443 )
  W( 95 ) = W( 95 ) + a*JVS( 444 )
  W( 96 ) = W( 96 ) + a*JVS( 445 )
  W( 97 ) = W( 97 ) + a*JVS( 446 )
  a = -W( 85 ) / JVS(          460  )
  W( 85 ) = -a
  W( 86 ) = W( 86 ) + a*JVS( 461 )
  W( 87 ) = W( 87 ) + a*JVS( 462 )
  W( 88 ) = W( 88 ) + a*JVS( 463 )
  W( 89 ) = W( 89 ) + a*JVS( 464 )
  W( 90 ) = W( 90 ) + a*JVS( 465 )
  W( 91 ) = W( 91 ) + a*JVS( 466 )
  W( 92 ) = W( 92 ) + a*JVS( 467 )
  W( 93 ) = W( 93 ) + a*JVS( 468 )
  W( 94 ) = W( 94 ) + a*JVS( 469 )
  W( 95 ) = W( 95 ) + a*JVS( 470 )
  W( 96 ) = W( 96 ) + a*JVS( 471 )
  W( 97 ) = W( 97 ) + a*JVS( 472 )
  a = -W( 86 ) / JVS(          483  )
  W( 86 ) = -a
  W( 87 ) = W( 87 ) + a*JVS( 484 )
  W( 88 ) = W( 88 ) + a*JVS( 485 )
  W( 90 ) = W( 90 ) + a*JVS( 486 )
  W( 91 ) = W( 91 ) + a*JVS( 487 )
  W( 93 ) = W( 93 ) + a*JVS( 488 )
  W( 94 ) = W( 94 ) + a*JVS( 489 )
  W( 95 ) = W( 95 ) + a*JVS( 490 )
  W( 96 ) = W( 96 ) + a*JVS( 491 )
  W( 97 ) = W( 97 ) + a*JVS( 492 )
  a = -W( 87 ) / JVS(          517  )
  W( 87 ) = -a
  W( 88 ) = W( 88 ) + a*JVS( 518 )
  W( 89 ) = W( 89 ) + a*JVS( 519 )
  W( 90 ) = W( 90 ) + a*JVS( 520 )
  W( 91 ) = W( 91 ) + a*JVS( 521 )
  W( 92 ) = W( 92 ) + a*JVS( 522 )
  W( 93 ) = W( 93 ) + a*JVS( 523 )
  W( 94 ) = W( 94 ) + a*JVS( 524 )
  W( 95 ) = W( 95 ) + a*JVS( 525 )
  W( 96 ) = W( 96 ) + a*JVS( 526 )
  W( 97 ) = W( 97 ) + a*JVS( 527 )
  a = -W( 88 ) / JVS(          549  )
  W( 88 ) = -a
  W( 89 ) = W( 89 ) + a*JVS( 550 )
  W( 90 ) = W( 90 ) + a*JVS( 551 )
  W( 91 ) = W( 91 ) + a*JVS( 552 )
  W( 92 ) = W( 92 ) + a*JVS( 553 )
  W( 93 ) = W( 93 ) + a*JVS( 554 )
  W( 94 ) = W( 94 ) + a*JVS( 555 )
  W( 95 ) = W( 95 ) + a*JVS( 556 )
  W( 96 ) = W( 96 ) + a*JVS( 557 )
  W( 97 ) = W( 97 ) + a*JVS( 558 )
  a = -W( 89 ) / JVS(          582  )
  W( 89 ) = -a
  W( 90 ) = W( 90 ) + a*JVS( 583 )
  W( 91 ) = W( 91 ) + a*JVS( 584 )
  W( 92 ) = W( 92 ) + a*JVS( 585 )
  W( 93 ) = W( 93 ) + a*JVS( 586 )
  W( 94 ) = W( 94 ) + a*JVS( 587 )
  W( 95 ) = W( 95 ) + a*JVS( 588 )
  W( 96 ) = W( 96 ) + a*JVS( 589 )
  W( 97 ) = W( 97 ) + a*JVS( 590 )
  a = -W( 90 ) / JVS(          630  )
  W( 90 ) = -a
  W( 91 ) = W( 91 ) + a*JVS( 631 )
  W( 92 ) = W( 92 ) + a*JVS( 632 )
  W( 93 ) = W( 93 ) + a*JVS( 633 )
  W( 94 ) = W( 94 ) + a*JVS( 634 )
  W( 95 ) = W( 95 ) + a*JVS( 635 )
  W( 96 ) = W( 96 ) + a*JVS( 636 )
  W( 97 ) = W( 97 ) + a*JVS( 637 )
  a = -W( 91 ) / JVS(          651  )
  W( 91 ) = -a
  W( 92 ) = W( 92 ) + a*JVS( 652 )
  W( 93 ) = W( 93 ) + a*JVS( 653 )
  W( 94 ) = W( 94 ) + a*JVS( 654 )
  W( 95 ) = W( 95 ) + a*JVS( 655 )
  W( 96 ) = W( 96 ) + a*JVS( 656 )
  W( 97 ) = W( 97 ) + a*JVS( 657 )
  a = -W( 92 ) / JVS(          672  )
  W( 92 ) = -a
  W( 93 ) = W( 93 ) + a*JVS( 673 )
  W( 94 ) = W( 94 ) + a*JVS( 674 )
  W( 95 ) = W( 95 ) + a*JVS( 675 )
  W( 96 ) = W( 96 ) + a*JVS( 676 )
  W( 97 ) = W( 97 ) + a*JVS( 677 )
  a = -W( 93 ) / JVS(          696  )
  W( 93 ) = -a
  W( 94 ) = W( 94 ) + a*JVS( 697 )
  W( 95 ) = W( 95 ) + a*JVS( 698 )
  W( 96 ) = W( 96 ) + a*JVS( 699 )
  W( 97 ) = W( 97 ) + a*JVS( 700 )
  a = -W( 94 ) / JVS(          715  )
  W( 94 ) = -a
  W( 95 ) = W( 95 ) + a*JVS( 716 )
  W( 96 ) = W( 96 ) + a*JVS( 717 )
  W( 97 ) = W( 97 ) + a*JVS( 718 )
  a = -W( 95 ) / JVS(          738  )
  W( 95 ) = -a
  W( 96 ) = W( 96 ) + a*JVS( 739 )
  W( 97 ) = W( 97 ) + a*JVS( 740 )
  a = -W( 96 ) / JVS(          792  )
  W( 96 ) = -a
  W( 97 ) = W( 97 ) + a*JVS( 793 )
  JVS( 794) = W( 46 )
  JVS( 795) = W( 53 )
  JVS( 796) = W( 69 )
  JVS( 797) = W( 73 )
  JVS( 798) = W( 74 )
  JVS( 799) = W( 75 )
  JVS( 800) = W( 76 )
  JVS( 801) = W( 77 )
  JVS( 802) = W( 80 )
  JVS( 803) = W( 81 )
  JVS( 804) = W( 82 )
  JVS( 805) = W( 83 )
  JVS( 806) = W( 84 )
  JVS( 807) = W( 85 )
  JVS( 808) = W( 86 )
  JVS( 809) = W( 87 )
  JVS( 810) = W( 88 )
  JVS( 811) = W( 89 )
  JVS( 812) = W( 90 )
  JVS( 813) = W( 91 )
  JVS( 814) = W( 92 )
  JVS( 815) = W( 93 )
  JVS( 816) = W( 94 )
  JVS( 817) = W( 95 )
  JVS( 818) = W( 96 )
  JVS( 819) = W( 97 )
   
   END SUBROUTINE decomp_cb05_sorg_aq
 


END MODULE cb05_sorg_aq_Integrator
