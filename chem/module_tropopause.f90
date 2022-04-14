













 module module_tropopause

      implicit none

      private

      public  :: tropopause_init
      public  :: tropopause_driver

      save

      integer,  parameter :: r8    = selected_real_kind(12) 

      real(r8), parameter :: pi    =3.14159265358979323846_r8
      real(r8), parameter :: mb2Pa = 100._r8
      real(r8), parameter :: d2r   = pi/180._r8
      real(r8), parameter :: zero  = 0._r8
      real(r8), parameter :: twopi = pi*2._r8

      integer,  parameter :: NOTFOUND = -1
 
      integer :: iend 
      integer :: jend

      integer :: tropo_month_n      
      integer :: tropo_lon_n        
      integer :: tropo_lat_n        

      type tropo_type
          real(r8), pointer :: tropo_p_loc(:,:,:)  
                                                   
          logical           :: is_allocated 
      end type tropo_type 

      type(tropo_type), allocatable :: tropo_bc(:)

      
      
      
      

      real(r8),parameter :: CONST_BOLTZ  = 1.38065e-23_r8  
      real(r8),parameter :: CONST_AVOGAD = 6.02214e26_r8   
      real(r8),parameter :: CONST_RGAS   = CONST_AVOGAD*CONST_BOLTZ 
      real(r8),parameter :: CONST_MWDAIR = 28.966_r8       
      real(r8),parameter :: CONST_RDAIR  = CONST_RGAS/CONST_MWDAIR 
      real(r8),parameter :: CONST_CPDAIR  = 1.00464e3_r8    

      real(r8),parameter :: gravit       = 9.80616_r8
      real(r8),parameter :: rair         = CONST_RDAIR
      real(r8),parameter :: cappa        = CONST_RDAIR/CONST_CPDAIR   

      real(r8),parameter :: cnst_kap    = cappa
      real(r8),parameter :: cnst_faktor = -gravit/rair
      real(r8),parameter :: cnst_ka1    = cnst_kap - 1._r8

      contains

      subroutine tropopause_driver( id, dtstep, current_date_char,  &
                                    t_phy, p_phy, p8w, zmid, z8w,   &
                                    tropo_lev, tropo_p, tropo_z,    &
                                    ids, ide, jds, jde, kds, kde,   &
                                    ims, ime, jms, jme, kms, kme,   &
                                    its, ite, jts, jte, kts, kte    )

      implicit none




      integer, intent(in   ) :: id,                             &
                                ids,ide, jds,jde, kds,kde,      &
                                ims,ime, jms,jme, kms,kme,      &
                                its,ite, jts,jte, kts,kte

      real, intent(in  ) :: dtstep

      real, dimension( ims:ime, kms:kme, jms:jme ),             &
               intent(in   ) :: t_phy,  & 
                                p_phy,  & 
                                zmid,   & 
                                z8w,    & 
                                p8w       

      real, dimension( ims:ime, jms:jme ),             &
               intent(inout) :: tropo_p,  & 
                                tropo_z     
      integer, dimension( ims:ime, jms:jme ),          &
               intent(inout) :: tropo_lev   

      CHARACTER (LEN=256),intent(in) :: current_date_char




      integer :: i, j, k





      iend = min(ite,ide-1)
      jend = min(jte,jde-1)

      tropo_lev(its:iend,jts:jend) = NOTFOUND


       
      call tropopause_twmo (id, t_phy, p_phy, p8w, zmid, z8w, &
                            tropo_lev, tropo_p, tropo_z,      &
                            ids,ide, jds,jde, kds,kde,        &
                            ims,ime, jms,jme, kms,kme,        &
                            its,ite, jts,jte, kts,kte         )



      if ( any(tropo_lev(its:iend,jts:jend) == NOTFOUND) ) then

         call tropopause_climate (id, current_date_char,       &
                                  p_phy, p8w, zmid, z8w,       &
                                  tropo_lev, tropo_p, tropo_z, &
                                  ids,ide, jds,jde, kds,kde,   &
                                  ims,ime, jms,jme, kms,kme,   &
                                  its,ite, jts,jte, kts,kte    )
      end if

   end subroutine tropopause_driver



      subroutine tropopause_init (id, xlat, xlon, config_flags, &
                                  ids,ide, jds,jde, kds,kde,    &
                                  ims,ime, jms,jme, kms,kme,    &
                                  its,ite, jts,jte, kts,kte     )



      use module_interpolate,  only : lininterp_init, lininterp, interp_type, lininterp_finish
      use module_configure,    only : grid_config_rec_type

      implicit none




      integer, intent(in) :: id,                          &
                             ids,ide, jds,jde, kds,kde,   &
                             ims,ime, jms,jme, kms,kme,   &
                             its,ite, jts,jte, kts,kte

      real,    intent(in) :: xlat(ims:ime, jms:jme)
      real,    intent(in) :: xlon(ims:ime, jms:jme)
      type(grid_config_rec_type), intent(in) :: config_flags





      type(interp_type) :: lon_wgts, lat_wgts

      integer :: max_dom
      integer :: astat
      integer :: ncid
      integer :: varid
      integer :: dimid(3)
      integer :: start(3)
      integer :: count(3)
      integer :: dimid_lon
      integer :: dimid_lat
      integer :: dimid_month
      integer :: ndims  

      character(len=128) :: err_msg
      character(len=64)  :: filename
      character(len=3)   :: id_num

      character(len=80) :: attribute

      real(r8), allocatable :: tropo_p_in(:,:,:)  
      real(r8), allocatable :: tropo_lat(:)       
      real(r8), allocatable :: tropo_lon(:)       

      real(r8) :: wrf_lon(1)       
      real(r8) :: wrf_lat(1)       
      real(r8) :: tmp_out(1)       

      integer  :: i, j, m
      CHARACTER(LEN=132) :: message

      LOGICAL , EXTERNAL      :: wrf_dm_on_monitor

include 'netcdf.inc'



      iend = min(ite,ide-1)
      jend = min(jte,jde-1)




   if( id == 1 .and. .not. allocated(tropo_bc) ) then
       CALL nl_get_max_dom( 1,max_dom )

       allocate(tropo_bc(max_dom),stat=astat )
       if( astat /= 0 ) then
          CALL wrf_message( 'tropopause_init: failed to allocate tropo_bc' )
          CALL wrf_abort
       end if

       tropo_bc(:)%is_allocated = .false.
  endif

tropo_bc_allocated : &
  if( .not. tropo_bc(id)%is_allocated ) then





master_proc : &
   IF( wrf_dm_on_monitor() ) THEN
      write(id_num,'(i3)') 100+id
      write(message,*) 'tropopause_init: intializing domain ' // id_num(2:3)
      call wrf_message( trim(message) )
       




      filename = config_flags%trop_lev_inname
      if( filename == ' ' ) then
        call wrf_message( 'tropopause_init: input filename not specified in namelist' )
        call wrf_abort
      endif

      err_msg = 'tropopause_init: failed to open file ' // trim(filename)
      call handle_ncerr( nf_open( trim(filename), nf_noclobber, ncid ), trim(err_msg) )
      write(message,*) 'tropopause_init: open filename= ',filename
      call wrf_message( trim(message) )



      err_msg = 'tropopause_init: failed to get time id'
      call handle_ncerr( nf_inq_dimid( ncid, 'time', dimid_month ), trim(err_msg) )
      err_msg = 'tropopause_init: failed to get time'
      call handle_ncerr( nf_inq_dimlen( ncid, dimid_month, tropo_month_n ), trim(err_msg) )
      if( tropo_month_n /= 12 )then
       write(message,*) 'tropopause_init: number of months = ',tropo_month_n,'; expecting 12'
       call wrf_message( trim(message) )
       call wrf_abort
      end if

      err_msg = 'tropopause_init: failed to get lat id'
      call handle_ncerr( nf_inq_dimid( ncid, 'lat', dimid_lat ), trim(err_msg) )
      err_msg = 'tropopause_init: failed to get lat'
      call handle_ncerr( nf_inq_dimlen( ncid, dimid_lat, tropo_lat_n ), trim(err_msg) )

      err_msg = 'tropopause_init: failed to get lon id'
      call handle_ncerr( nf_inq_dimid( ncid, 'lon', dimid_lon ), trim(err_msg) )
      err_msg = 'tropopause_init: failed to get lon'
      call handle_ncerr( nf_inq_dimlen( ncid, dimid_lon, tropo_lon_n ), trim(err_msg) )

   END IF master_proc



      CALL wrf_dm_bcast_integer ( tropo_month_n , 1 )
      CALL wrf_dm_bcast_integer ( tropo_lat_n   , 1 )
      CALL wrf_dm_bcast_integer ( tropo_lon_n   , 1 )





      allocate( tropo_lat(tropo_lat_n), stat=astat )
      if( astat /= 0 ) then
         call wrf_message( 'tropopause_init: failed to allocate tropo_lat' )
         call wrf_abort
      end if

      allocate( tropo_lon(tropo_lon_n), stat=astat )
      if( astat /= 0 ) then
         call wrf_message( 'tropopause_init: failed to allocate tropo_lon' )
         call wrf_abort
      end if

      allocate( tropo_p_in(tropo_lon_n,tropo_lat_n,tropo_month_n), stat=astat )
      if( astat /= 0 ) then
         call wrf_message( 'tropopause_init: failed to allocate tropo_p_in' )
         call wrf_abort
      end if
 




      allocate( tropo_bc(id)%tropo_p_loc(its:iend,jts:jend,tropo_month_n), stat=astat )
      if( astat /= 0 ) then
         call wrf_message( 'tropopause_init: failed to allocate tropo_bc(id)%tropo_p_loc' )
         call wrf_abort
      end if  

      tropo_bc(id)%is_allocated = .true.




master_proc_a : &
   IF ( wrf_dm_on_monitor() ) THEN


      err_msg = 'tropopause_init: failed to get lat variable id'
      call handle_ncerr( nf_inq_varid( ncid, 'lat', varid ), trim(err_msg) )
      err_msg = 'tropopause_init: failed to read lat variable'
      call handle_ncerr( nf_get_var_double( ncid, varid, tropo_lat ), trim(err_msg) )

      
     tropo_lat(:) = tropo_lat(:) * d2r



      err_msg = 'tropopause_init: failed to get lon variable id'
      call handle_ncerr( nf_inq_varid( ncid, 'lon', varid ), trim(err_msg) )
      err_msg = 'tropopause_init: failed to read lon variable'
      call handle_ncerr( nf_get_var_double( ncid, varid, tropo_lon ), trim(err_msg) )

      
      tropo_lon(:) = tropo_lon(:) * d2r


      err_msg = 'tropopause_init: failed to get trop_p variable id'
      call handle_ncerr( nf_inq_varid( ncid, 'trop_p', varid ), trim(err_msg) )

      

      err_msg = 'tropopause_init: failed to get ndims of trop_p variable'
      call handle_ncerr( nf_inq_varndims( ncid, varid, ndims ), trim(err_msg) )

      if( ndims /= 3 ) then
         write(message,*) 'tropopause_init: error! variable trop_p has ndims = ',ndims,', expecting 3'
         call wrf_message( trim(message) )
         call wrf_abort
      end if

      err_msg = 'tropopause_init: failed to get dimid of vmr variable'
      call handle_ncerr( nf_inq_vardimid( ncid, varid, dimid ), trim(err_msg) )

      if( dimid(1) /= dimid_lon   .or. dimid(2) /= dimid_lat .or. &
          dimid(3) /= dimid_month )then
          write(message,*) 'tropopause_init: error! dimensions in wrong order for variable trop_p,'// &
               'expecting (lon,lat,month)'
          call wrf_message( trim(message) )
          call wrf_abort
      end if

      
      
      
      err_msg = 'tropopause_init: failed to read trop_p variable'
      call handle_ncerr( nf_get_var_double( ncid, varid, tropo_p_in ), trim(err_msg) )




      err_msg = 'tropopause_init: failed to close file ' // trim(filename)
      call handle_ncerr( nf_close( ncid ), trim(err_msg) )

   END IF master_proc_a




      CALL wrf_dm_bcast_double ( tropo_lat, size(tropo_lat) )
      CALL wrf_dm_bcast_double ( tropo_lon, size(tropo_lon) )
      CALL wrf_dm_bcast_double ( tropo_p_in, size(tropo_p_in) )




      
      
      
      
      
      
      

      
      
      
      do i = its,iend
        do j = jts,jend   
         
         
         
          wrf_lat(1) = xlat(i,j)*d2r

          wrf_lon(1) = xlon(i,j)     
          if( wrf_lon(1)  < 0.0_r8 ) wrf_lon(1) = wrf_lon(1) + 360.0_r8
          wrf_lon(1) = wrf_lon(1)*d2r
   
         
         
         
          call lininterp_init( tropo_lon, tropo_lon_n, wrf_lon, 1, 2, lon_wgts, zero, twopi )
          call lininterp_init( tropo_lat, tropo_lat_n, wrf_lat, 1, 1, lat_wgts )

         
         
         
          do m = 1,tropo_month_n
            call lininterp( tropo_p_in(:,:,m), tropo_lon_n, tropo_lat_n, &
                            tmp_out, 1 , lon_wgts, lat_wgts)
            tropo_bc(id)%tropo_p_loc(i,j,m) = tmp_out(1)    
          end do

        end do
      end do

      call lininterp_finish( lon_wgts )
      call lininterp_finish( lat_wgts )

      deallocate(tropo_lon)
      deallocate(tropo_lat)
      deallocate(tropo_p_in)
      
      call wrf_message( ' ' )
      write(message,*) 'tropopause_init: DONE intialized domain ',id
      call wrf_message( trim(message) )
      call wrf_message( ' ' )

endif tropo_bc_allocated

      end subroutine tropopause_init



      subroutine handle_ncerr( ret, mes )




      implicit none




      integer, intent(in) :: ret
      character(len=*), intent(in) :: mes

include 'netcdf.inc'

      if( ret /= nf_noerr ) then
         call wrf_message( trim(mes) )
         call wrf_message( trim(nf_strerror(ret)) )
         call wrf_abort
      end if

      end subroutine handle_ncerr



  
  
  

  subroutine tropopause_twmo (id, t_phy, p_phy, p8w, zmid, z8w, & 
                              tropo_lev, tropo_p,tropo_z,       &
                              ids,ide, jds,jde, kds,kde,        &
                              ims,ime, jms,jme, kms,kme,        &
                              its,ite, jts,jte, kts,kte         )
      implicit none




      integer, intent(in) :: id,                         &
                             ids,ide, jds,jde, kds,kde,  &
                             ims,ime, jms,jme, kms,kme,  &
                             its,ite, jts,jte, kts,kte

      real, dimension( ims:ime , kms:kme , jms:jme ),    &
            intent(in)    :: t_phy,  & 
                             p_phy,  & 
                             zmid,   & 
                             z8w,    & 
                             p8w       

      real, dimension( ims:ime, jms:jme ),             &
               intent(inout) :: tropo_p,  & 
                                tropo_z     
      integer, dimension( ims:ime, jms:jme ),          &
               intent(inout) :: tropo_lev   



    real, dimension( kts:kte) :: t_temp    
    real, dimension( kts:kte) :: p_temp    
 
    real(r8), parameter     :: gamma    = -0.002_r8    
    real(r8), parameter     :: plimu    = 45000._r8    
    real(r8), parameter     :: pliml    = 7500._r8     
     
    integer                 :: i
    integer                 :: j
    integer                 :: k
    integer                 :: kk
    integer                 :: pLev

    real(r8)                :: tP                       
    real(r8)                :: dZdlogP


    pLev = kte - kts + 1

    
    do i = its,iend
    do j = jts,jend

       
       

       do k = kts,kte
          kk = pLev - k + 1
          t_temp(kk) = t_phy(i,k,j)
          p_temp(kk) = p_phy(i,k,j)
       end do
 
       

       call twmo(pLev, t_temp, p_temp, plimu, pliml, gamma, tP)
     
       
       if (tP > 0) then        
          
          do k = kts,kte-1
             if (tP >= p8w(i,k,j)) then
                tropo_lev(i,j) = k - 1      
                tropo_p  (i,j) = tP
                exit
             end if
          end do
          
          
          

          k = k - 1     
    
          
          if (tropo_p(i,j) == p_phy(i,k,j)) then
             tropo_z(i,j)= zmid(i,k,j)
          else if (tropo_p(i,j) < p_phy(i,k,j)) then

             if ( k > kts ) then 
                dZdlogP = (zmid(i,k,j) - z8w(i,k-1,j)) / &
                          (log(p_phy(i,k,j)) - log(p8w(i,k-1,j)) )
                tropo_z(i,j)= zmid(i,k,j) + (log(tropo_p(i,j)) - log(p_phy(i,k,j))) * dZdlogP
             end if
          else

             if ( k < kte ) then
                dZdlogP = (zmid(i,k,j) - z8w(i,k,j)) / &
                          (log(p_phy(i,k,j)) - log(p8w(i,k,j)) )
                tropo_z(i,j)= zmid(i,k,j) + (log(tropo_p(i,j)) - log(p_phy(i,k,j))) * dZdlogP
             end if
          end if
          
       end if





    end do
    end do

  end subroutine tropopause_twmo



  
  
  
  
  

  
  
  
  
  
  
  
  
  
  

  subroutine twmo(level, t, p, plimu, pliml, gamma, trp)

    integer, intent(in)                     :: level


    real,     intent(in), dimension(level)  :: t, p

    real(r8), intent(in)                    :: plimu, pliml, gamma
    real(r8), intent(out)                   :: trp
    
    real(r8), parameter                     :: deltaz = 2000.0_r8

    real(r8)                                :: pmk, pm, a, b, tm, dtdp, dtdz
    real(r8)                                :: ag, bg, ptph
    real(r8)                                :: pm0, pmk0, dtdz0
    real(r8)                                :: p2km, asum, aquer
    real(r8)                                :: pmk2, pm2, a2, b2, tm2, dtdp2, dtdz2
    integer                                 :: icount, jj
    integer                                 :: j
    

    trp=-99.0_r8                           

lev_loop : &
    do j=level,2,-1
    
      
      pmk= .5_r8 * (p(j-1)**cnst_kap+p(j)**cnst_kap)
      pm = pmk**(1/cnst_kap)               
      a = (t(j-1)-t(j))/(p(j-1)**cnst_kap-p(j)**cnst_kap)
      b = t(j)-(a*p(j)**cnst_kap)
      tm = a * pmk + b               
      dtdp = a * cnst_kap * (pm**cnst_ka1)
      dtdz = cnst_faktor*dtdp*pm/tm
      



      if( j == level .or. dtdz <= gamma .or. pm > plimu ) then
        pm0   = pm
        pmk0  = pmk
        dtdz0 = dtdz
        cycle lev_loop
      endif
  
      
      if (dtdz0 < gamma) then
        ag = (dtdz-dtdz0) / (pmk-pmk0)     
        bg = dtdz0 - (ag * pmk0)          
        ptph = exp(log((gamma-bg)/ag)/cnst_kap)
      else
        ptph = pm
      endif
  


      if( ptph < pliml .or. ptph > plimu ) then
        pm0   = pm
        pmk0  = pmk
        dtdz0 = dtdz
        cycle lev_loop
      endif
  
      
      p2km = ptph + deltaz*(pm/tm)*cnst_faktor     
      asum = 0.0_r8                                
      icount = 0                                   
  
      
lev_loop_a : &
      do jj=j,2,-1
    
        pmk2 = .5_r8 * (p(jj-1)**cnst_kap+p(jj)**cnst_kap) 
        pm2 = pmk2**(1/cnst_kap)                           
        if( pm2 > ptph ) then                     
          cycle lev_loop_a
        endif
        if( pm2 < p2km ) then                     
          trp = ptph
          exit lev_loop
        endif

        a2 = (t(jj-1)-t(jj))                     
        a2 = a2/(p(jj-1)**cnst_kap-p(jj)**cnst_kap)
        b2 = t(jj)-(a2*p(jj)**cnst_kap)          
        tm2 = a2 * pmk2 + b2                     
        dtdp2 = a2 * cnst_kap * (pm2**(cnst_kap-1))  
        dtdz2 = cnst_faktor*dtdp2*pm2/tm2
        asum = asum+dtdz2
        icount = icount+1
        aquer = asum/float(icount)               
   
        
        if( aquer <= gamma ) then                
          pm0   = pm
          pmk0  = pmk
          dtdz0 = dtdz
          exit lev_loop_a
        endif
    
      enddo lev_loop_a
    
    enddo lev_loop

  end subroutine twmo



  
  
  
  
  

  subroutine tropopause_climate (id, current_date_char,        &
                                 p_phy, p8w, zmid, z8w,        &
                                 tropo_lev, tropo_p,  tropo_z, &
                                 ids,ide, jds,jde, kds,kde,    &
                                 ims,ime, jms,jme, kms,kme,    &
                                 its,ite, jts,jte, kts,kte     )
      implicit none




      integer, intent(in) :: id,                         &
                             ids,ide, jds,jde, kds,kde,  &
                             ims,ime, jms,jme, kms,kme,  &
                             its,ite, jts,jte, kts,kte

      real, dimension( ims:ime , kms:kme , jms:jme ),    &
            intent(in)    :: p_phy,  & 
                             zmid,   & 
                             z8w,    & 
                             p8w       

      real, dimension( ims:ime, jms:jme ),             &
               intent(inout) :: tropo_p,  & 
                                tropo_z     
      integer, dimension( ims:ime, jms:jme ),          &
               intent(inout) :: tropo_lev   

      CHARACTER (LEN=256),intent(in) :: current_date_char


    

    real      :: del_time
    integer   :: last
    integer   :: next

    real(r8)  :: tP                       
    real(r8)  :: dZdlogP

    integer   :: i
    integer   :: j
    integer   :: k
    CHARACTER(LEN=132) :: message

    
    if (any(tropo_lev(its:iend,jts:jend) == NOTFOUND)) then

       
       
       

       call get_time_interp_factors( current_date_char, last, next, del_time )

       

       do i = its,iend
       do j = jts,jend
       
          

          if (tropo_lev(i,j) == NOTFOUND) then
        
             
             
             
             
             tP =  tropo_bc(id)%tropo_p_loc(i,j,last) &
                + (tropo_bc(id)%tropo_p_loc(i,j,next) &
                -  tropo_bc(id)%tropo_p_loc(i,j,last))*del_time

             if (tP > 0) then        
                
                do k = kts,kte-1
                   if (tP >= p8w(i,k,j)) then
                      tropo_lev(i,j) = k - 1   
                      tropo_p  (i,j) = tP
                      exit
                   end if
                end do

                
                
                

                k = k - 1     
    
                
                if (tropo_p(i,j) == p_phy(i,k,j)) then
                   tropo_z(i,j)= zmid(i,k,j)
    
                else if (tropo_p(i,j) < p_phy(i,k,j)) then

                   if ( k > kts ) then 
                      dZdlogP = (zmid(i,k,j) - z8w(i,k-1,j)) / &
                                (log(p_phy(i,k,j)) - log(p8w(i,k-1,j)) )
                      tropo_z(i,j)= zmid(i,k,j) + (log(tropo_p(i,j)) - log(p_phy(i,k,j))) * dZdlogP
                   end if
                else

                   if ( k < kte ) then 
                      dZdlogP = (zmid(i,k,j) - z8w(i,k,j)) / &
                                (log(p_phy(i,k,j)) - log(p8w(i,k,j)) )
                      tropo_z(i,j)= zmid(i,k,j) + (log(tropo_p(i,j)) - log(p_phy(i,k,j))) * dZdlogP
                   end if
                end if
                
             end if
          end if
       end do
       end do
    end if

    if (any(tropo_lev(its:iend,jts:jend) == NOTFOUND)) then       
       write(message,*) 'tropopause_climate: Warning: some tropopause levels still NOTFOUND' 
    else        
       write(message,*) 'tropopause_climate: Warning: Done finding tropopause'
    end if
    call wrf_message( trim(message) )
  
  end subroutine tropopause_climate
 


      subroutine get_time_interp_factors(current_date_char , last, next, del_time )




      use module_date_time, only : get_julgmt

      implicit none




      CHARACTER (LEN=256),intent(in) :: current_date_char

      integer, intent(out) :: next, last
      real,    intent(out) :: del_time




      integer, parameter :: day_of_year(12) = (/ 16, 45, 75, 105, 136, 166, 197, &
                                                 228, 258, 289, 319, 350        /)
      integer :: julyr , julday 
      real    :: gmt
      real    :: calday

      integer  :: m


      call get_julgmt ( current_date_char, julyr, julday, gmt )

      calday = real(julday) + gmt

      if( calday < day_of_year(1) ) then
         next = 1
         last = 12
         del_time = (365. + calday - day_of_year(12)) &
                / (365. + day_of_year(1) - day_of_year(12))
      else if( calday >= day_of_year(12) ) then
         next = 1
         last = 12
         del_time = (calday - day_of_year(12)) &
                / (365. + day_of_year(1) - day_of_year(12))
      else
         do m = 11,1,-1
            if( calday >= day_of_year(m) ) then
               exit
            end if
         end do
         last = m
         next = m + 1
         del_time = (calday - day_of_year(m)) / (day_of_year(m+1) - day_of_year(m))
      end if

      del_time = max( min( 1.,del_time ),0. )

      end subroutine get_time_interp_factors


      end module module_tropopause

