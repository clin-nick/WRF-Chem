

































MODULE module_input_tracer
USE module_input_tracer_data
USE module_state_description, only:tracer_smoke,tracer_test1,tracer_test2,param_first_scalar,p_tr17_1,p_tr17_2,p_tr17_3,p_tr17_4,p_tr17_5,p_tr17_6,p_tr17_7,p_tr17_8
CONTAINS
   SUBROUTINE initialize_tracer (chem,chem_in_opt,         &
                                       tracer_opt,num_chem,&
                               ids,ide, jds,jde, kds,kde,  & 
                               ims,ime, jms,jme, kms,kme,  & 
                               ips,ipe, jps,jpe, kps,kpe,  & 
                               its,ite, jts,jte, kts,kte )
      INTEGER,      INTENT(IN   )    :: chem_in_opt,tracer_opt,num_chem
      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: ips,ipe, jps,jpe, kps,kpe
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte
      REAL,  DIMENSION(ims:ime,kms:kme,jms:jme,num_chem ), INTENT(INOUT) :: chem
      if(chem_in_opt == 1 )return
      if     (tracer_opt == TRACER_TEST1)then
       chem(:,:,:,:)=.0
      else if(tracer_opt == TRACER_TEST2)then
       chem(:,:,:,:)=.0
      else if(tracer_opt == TRACER_SMOKE)then
       chem(:,:,:,:)=.08
      endif
   END SUBROUTINE initialize_tracer
   SUBROUTINE flow_dep_bdy_tracer  (  chem,                                       &
                               chem_bxs,chem_btxs,                                  &
                               chem_bxe,chem_btxe,                                  &
                               chem_bys,chem_btys,                                  &
                               chem_bye,chem_btye,                                  &
                               dt,                                              &
                               spec_bdy_width,z,                                &
                               have_bcs_chem,                        & 
                               u, v, tracer_opt, alt, & 
                               t,pb,p,t0,p1000mb,rcp,ph,phb,g, &
                               spec_zone, ic,           &
                               ids,ide, jds,jde, kds,kde,  & 
                               ims,ime, jms,jme, kms,kme,  & 
                               ips,ipe, jps,jpe, kps,kpe,  & 
                               its,ite, jts,jte, kts,kte )







      IMPLICIT NONE

      INTEGER,      INTENT(IN   )    :: tracer_opt
      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: ips,ipe, jps,jpe, kps,kpe
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte
      INTEGER,      INTENT(IN   )    :: spec_zone,spec_bdy_width,ic
      REAL,         INTENT(IN   )    :: dt


      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(INOUT) :: chem
      REAL,  DIMENSION( jms:jme , kds:kde , spec_bdy_width), INTENT(IN   ) :: chem_bxs, chem_bxe, chem_btxs, chem_btxe
      REAL,  DIMENSION( ims:ime , kds:kde , spec_bdy_width), INTENT(IN   ) :: chem_bys, chem_bye, chem_btys, chem_btye
      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: z
      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: alt
      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: u
      REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ), INTENT(IN   ) :: v
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,         &
          INTENT(IN   ) ::                                           &
                               ph,phb,t,pb,p
   real, INTENT (IN) :: g,rcp,t0,p1000mb

      INTEGER    :: i, j, k, numgas
      INTEGER    :: ibs, ibe, jbs, jbe, itf, jtf, ktf
      INTEGER    :: i_inner, j_inner
      INTEGER    :: b_dist
      integer    :: i_bdy_method
      real tempfac,convfac
      logical, optional    :: have_bcs_chem

      ibs = ids
      ibe = ide-1
      itf = min(ite,ide-1)
      jbs = jds
      jbe = jde-1
      jtf = min(jte,jde-1)
      ktf = kde-1



      i_bdy_method = 0
        if (tracer_opt == TRACER_TEST1 ) then
          i_bdy_method = 2
        end if   
        if (tracer_opt == TRACER_TEST2 ) then
          i_bdy_method = 2
        end if   
        if (tracer_opt == TRACER_SMOKE ) then
          i_bdy_method = 1
        end if
      if (have_bcs_chem) i_bdy_method =6
      if (ic .lt. param_first_scalar) i_bdy_method = 0

      IF (jts - jbs .lt. spec_zone) THEN

        DO j = jts, min(jtf,jbs+spec_zone-1)
          b_dist = j - jbs
          DO k = kts, ktf
            DO i = max(its,b_dist+ibs), min(itf,ibe-b_dist)
              i_inner = max(i,ibs+spec_zone)
              i_inner = min(i_inner,ibe-spec_zone)
              IF(v(i,k,j) .lt. 0.)THEN
                chem(i,k,j) = chem(i_inner,k,jbs+spec_zone)
              ELSE
                if (i_bdy_method .eq. 0) then
                   chem(i,k,j) = tracer_bv_def
                else if (i_bdy_method .eq. 1) then
                   chem(i,k,j)=tr_smoke_value
                else if (i_bdy_method .eq. 2) then
                   if (ic .eq. p_tr17_1 .or. ic .eq. p_tr17_2) then
                      chem(i,k,j)= tracer_bv_one
                   else
                      chem(i,k,j)= tracer_bv_def
                   endif
                else if (i_bdy_method .eq. 6) then
                   CALL bdy_tracer_value ( chem(i,k,j),chem_bys(i,k,1),chem_btys(i,k,1),dt,ic)
                else
                   chem(i,k,j) = tracer_bv_def
                endif
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF 
      IF (jbe - jtf .lt. spec_zone) THEN 

        DO j = max(jts,jbe-spec_zone+1), jtf 
          b_dist = jbe - j 
          DO k = kts, ktf 
            DO i = max(its,b_dist+ibs), min(itf,ibe-b_dist)
              i_inner = max(i,ibs+spec_zone)
              i_inner = min(i_inner,ibe-spec_zone)
              IF(v(i,k,j+1) .gt. 0.)THEN
                chem(i,k,j) = chem(i_inner,k,jbe-spec_zone)
              ELSE
                if (i_bdy_method .eq. 0) then
                   chem(i,k,j) = tracer_bv_def
                else if (i_bdy_method .eq. 1) then
                   chem(i,k,j)=tr_smoke_value
                else if (i_bdy_method .eq. 2) then
                   if (ic .eq. p_tr17_1 .or. ic .eq. p_tr17_2) then
                      chem(i,k,j)= tracer_bv_one
                   else
                      chem(i,k,j)= tracer_bv_def
                   endif
                else if (i_bdy_method .eq. 6) then
                   CALL bdy_tracer_value ( chem(i,k,j),chem_bye(i,k,1),chem_btye(i,k,1),dt,ic)
                else
                   chem(i,k,j) = tracer_bv_def
                endif
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF 

      IF (its - ibs .lt. spec_zone) THEN

        DO i = its, min(itf,ibs+spec_zone-1)
          b_dist = i - ibs
          DO k = kts, ktf
            DO j = max(jts,b_dist+jbs+1), min(jtf,jbe-b_dist-1)
              j_inner = max(j,jbs+spec_zone)
              j_inner = min(j_inner,jbe-spec_zone)
              IF(u(i,k,j) .lt. 0.)THEN
                chem(i,k,j) = chem(ibs+spec_zone,k,j_inner)
              ELSE
                if (i_bdy_method .eq. 0) then
                   chem(i,k,j) = tracer_bv_def
                else if (i_bdy_method .eq. 1) then
                   chem(i,k,j)=tr_smoke_value
                else if (i_bdy_method .eq. 2) then
                   if (ic .eq. p_tr17_1 .or. ic .eq. p_tr17_2) then
                      chem(i,k,j)= tracer_bv_one
                   else
                      chem(i,k,j)= tracer_bv_def
                   endif
                else if (i_bdy_method .eq. 6) then
                   CALL bdy_tracer_value ( chem(i,k,j),chem_bxs(j,k,1),chem_btxs(j,k,1),dt,ic)   
                else
                   chem(i,k,j) = tracer_bv_def
                endif
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF 

      IF (ibe - itf .lt. spec_zone) THEN

        DO i = max(its,ibe-spec_zone+1), itf
          b_dist = ibe - i
          DO k = kts, ktf
            DO j = max(jts,b_dist+jbs+1), min(jtf,jbe-b_dist-1)
              j_inner = max(j,jbs+spec_zone)
              j_inner = min(j_inner,jbe-spec_zone)
              IF(u(i+1,k,j) .gt. 0.)THEN
                chem(i,k,j) = chem(ibe-spec_zone,k,j_inner)
              ELSE
                if (i_bdy_method .eq. 0) then
                   chem(i,k,j) = tracer_bv_def
                else if (i_bdy_method .eq. 1) then
                   chem(i,k,j)=tr_smoke_value
                else if (i_bdy_method .eq. 2) then
                   if (ic .eq. p_tr17_1 .or. ic .eq. p_tr17_2) then
                      chem(i,k,j)= tracer_bv_one
                   else
                      chem(i,k,j)= tracer_bv_def
                   endif
                else if (i_bdy_method .eq. 6) then
                   CALL bdy_tracer_value ( chem(i,k,j),chem_bxe(j,k,1),chem_btxe(j,k,1),dt,ic)
                else
                   chem(i,k,j) = tracer_bv_def
                endif
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF 

   END SUBROUTINE flow_dep_bdy_tracer
   SUBROUTINE set_tracer(dtstep,ktau,pbl_h,tracer,t,tracer_opt,num_tracer,&
                         z,ht,ids,ide, jds,jde, kds,kde,                  & 
                               ims,ime, jms,jme, kms,kme,                 & 
                               its,ite, jts,jte, kts,kte                  )
      INTEGER,      INTENT(IN   )    :: ktau,tracer_opt,num_tracer
      INTEGER,      INTENT(IN   )    :: ids,ide, jds,jde, kds,kde
      INTEGER,      INTENT(IN   )    :: ims,ime, jms,jme, kms,kme
      INTEGER,      INTENT(IN   )    :: its,ite, jts,jte, kts,kte
      REAL,  DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer ), INTENT(INOUT) :: tracer
      REAL,  DIMENSION(ims:ime,kms:kme,jms:jme ), INTENT(IN) :: t,z
      REAL,  DIMENSION(ims:ime,jms:jme ), INTENT(IN) :: PBL_H,HT
      REAL,  INTENT(IN) :: dtstep
      INTEGER:: count_trop,count_pbl



    factor_decay = 1./(86400./dtstep)



    tracer(its:ite,kts:kte,jts:jte,p_tr17_2) = &
       tracer(its:ite,kts:kte,jts:jte,p_tr17_2) * (1. - factor_decay)

    tracer(its:ite,kts:kte,jts:jte,p_tr17_4) = &
       tracer(its:ite,kts:kte,jts:jte,p_tr17_4) * (1. - factor_decay)

    tracer(its:ite,kts:kte,jts:jte,p_tr17_6) = &
       tracer(its:ite,kts:kte,jts:jte,p_tr17_6) * (1. - factor_decay)

    tracer(its:ite,kts:kte,jts:jte,p_tr17_8) = &
       tracer(its:ite,kts:kte,jts:jte,p_tr17_8) * (1. - factor_decay)
 IF (ktau .ge. 2) THEN
    



    if(tracer_opt == TRACER_TEST1)then   
       tracer(its:ite,kts,jts:jte,p_tr17_3)     = 1.0
       tracer(its:ite,kts,jts:jte,p_tr17_4)     = 1.0
    endif
       
    do i= its,ite
    do j= jts,jte
 




       count_trop = minloc(t(i,kts:kte,j),1)

       tracer(i,count_trop:kte,j,p_tr17_5) = 1.0
       tracer(i,count_trop:kte,j,p_tr17_6) = 1.0





       count_pbl = 0

       do k=kts,kte
          if ( (z(i,k,j)-ht(i,j)) .le. pbl_h(i,j) ) then
             count_pbl = count_pbl + 1
          endif
       end do

       if (count_pbl .ge. 1) then
          tracer(i,kts:count_pbl,j,p_tr17_7) = 1.0
          tracer(i,kts:count_pbl,j,p_tr17_8) = 1.0
       endif

    end do   
    end do   

 ENDIF   
   END SUBROUTINE set_tracer

  SUBROUTINE bdy_tracer_value ( trac, trac_b, trac_bt, dt,ic)
                                  
    IMPLICIT NONE

    REAL,    intent(OUT)  :: trac
    REAL,    intent(IN)   :: trac_b
    REAL,    intent(IN)   :: trac_bt
    REAL,    intent(IN)   :: dt
    INTEGER, intent(IN)   :: ic

    REAL                  :: epsilc = 1.E-12








     
      trac=max(epsilc,trac_b + trac_bt * dt)

      RETURN
  END SUBROUTINE bdy_tracer_value

END MODULE module_input_tracer
