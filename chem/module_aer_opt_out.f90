MODULE module_aer_opt_out

   REAL,    PARAMETER, PRIVATE ::   afwalowv1   = 3.  
   REAL,    PARAMETER, PRIVATE ::   afwahiwv1   = 5.  
   REAL,    PARAMETER, PRIVATE ::   afwalowv2   = 8.  
   REAL,    PARAMETER, PRIVATE ::   afwahiwv2   = 12. 
CONTAINS
   SUBROUTINE aer_opt_out(dz8w  &
                   ,ext_coeff,bscat_coeff,asym_par                &
                   ,tauaer300,tauaer400,tauaer600,tauaer999       & 
                   ,gaer300,gaer400,gaer600,gaer999               & 
                   ,waer300,waer400,waer600,waer999               & 
                   ,num_ext_coef,num_bscat_coef,num_asym_par    &
                   ,ids,ide, jds,jde, kds,kde                     &
                   ,ims,ime, jms,jme, kms,kme                     &
                   ,its,ite, jts,jte, kts,kte )
USE module_configure, only:p_extcof3,p_extcof55,p_extcof106,p_extcof3_5,p_extcof8_12,p_bscof3,p_bscof55, &
           p_bscof106,p_asympar3,p_asympar55,p_asympar106

   IMPLICIT NONE
   INTEGER,    INTENT(IN   ) ::        ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte, &
                                       num_ext_coef,num_bscat_coef,num_asym_par
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_ext_coef ), INTENT (OUT) :: ext_coeff
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_bscat_coef ), INTENT (OUT) :: bscat_coeff
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_asym_par ), INTENT (OUT) :: asym_par
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),        &
         INTENT(IN    ) :: tauaer300,tauaer400,tauaer600,tauaer999, &
                                 gaer300,gaer400,gaer600,gaer999, & 
                              waer300,waer400,waer600,waer999,dz8w
   real :: ang,slope,slopeg,slopessa,onemang    
   integer :: i,j,k










      do j = jts,jte
      do k = kts,kte
      do i = its,ite



      ext_coeff(i,k,j,p_extcof3)=tauaer300(i,k,j)*1.E3/dz8w(i,k,j)  
      bscat_coeff(i,k,j,p_bscof3)=tauaer300(i,k,j)*waer300(i,k,j)*1.E3/dz8w(i,k,j)  
      asym_par(i,k,j,p_asympar3)=gaer300(i,k,j)  

           ang=log(tauaer300(i,k,j)/tauaer999(i,k,j))/log(999./300.)
           slopessa=(waer600(i,k,j)-waer400(i,k,j))/.2
           slopeg=(gaer600(i,k,j)-gaer400(i,k,j))/.2
      ext_coeff(i,k,j,p_extcof55)=tauaer400(i,k,j)*1.E3*((0.4/0.55)**ang)/dz8w(i,k,j)  
      slope= slopessa*(0.55-.6)+waer600(i,k,j) 
      slope=AMIN1(1.0,AMAX1(0.4,slope))   
      bscat_coeff(i,k,j,p_bscof55)=ext_coeff(i,k,j,p_extcof55)*slope  
      asym_par(i,k,j,p_asympar55)=AMIN1(1.,AMAX1(0.5,slopeg*(.55-.6)+gaer600(i,k,j)))  

           slopessa=(waer999(i,k,j)-waer600(i,k,j))/.399
           slopeg=(gaer999(i,k,j)-gaer600(i,k,j))/.399
      ext_coeff(i,k,j,p_extcof106)=tauaer400(i,k,j)*1.E3*((0.4/1.06)**ang)/dz8w(i,k,j)  
      slope= slopessa*(1.06-.999)+waer999(i,k,j) 
      slope=AMIN1(1.0,AMAX1(0.4,slope))   
      bscat_coeff(i,k,j,p_bscof106)=ext_coeff(i,k,j,p_extcof106)*slope  
      asym_par(i,k,j,p_asympar106)=AMIN1(1.,AMAX1(0.5,slopeg*(1.06-.999)+gaer600(i,k,j)))  

      onemang=1.-ang
      if(abs(onemang).gt.1.E-3)then   
      slope = tauaer400(i,k,j)*(0.4/afwalowv1)**ang  
      slopeg = tauaer400(i,k,j)*(0.4/afwahiwv1)**ang  
      ext_coeff(i,k,j,p_extcof3_5) = (slopeg*afwahiwv1-slope*afwalowv1)/(afwahiwv1-afwalowv1)/onemang
      slope = tauaer400(i,k,j)*(0.4/afwalowv2)**ang  
      slopeg = tauaer400(i,k,j)*(0.4/afwahiwv2)**ang  
      ext_coeff(i,k,j,p_extcof8_12) = (slopeg*afwahiwv2-slope*afwalowv2)/(afwahiwv2-afwalowv2)/onemang
      else                           
      ext_coeff(i,k,j,p_extcof3_5) = tauaer400(i,k,j)*0.4*log(afwahiwv1/afwalowv1)/(afwahiwv1-afwalowv1)
      ext_coeff(i,k,j,p_extcof8_12) = tauaer400(i,k,j)*0.4*log(afwahiwv2/afwalowv2)/(afwahiwv2-afwalowv2)
      endif

      ext_coeff(i,k,j,p_extcof3_5) = ext_coeff(i,k,j,p_extcof3_5)*1.E3/dz8w(i,k,j)
      ext_coeff(i,k,j,p_extcof8_12) = ext_coeff(i,k,j,p_extcof8_12)*1.E3/dz8w(i,k,j)
      end do
      end do
      end do
  END SUBROUTINE AER_OPT_OUT 
END MODULE module_aer_opt_out

