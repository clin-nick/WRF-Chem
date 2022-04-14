
      MODULE module_ftuv_subs

      use module_wave_data, only : nw

      implicit none

      integer, parameter :: dp = selected_real_kind(14,300)



      integer, private, save :: nintervals
      integer, private, allocatable, save :: xi(:), xcnt(:)
      real(dp), private, allocatable, save :: xfrac(:,:)

      real(dp), private, parameter :: km2cm = 1.e5_dp




      integer, private, parameter :: nla = 2
      integer, private, parameter :: ngast = 17




      integer, private, parameter  :: tdim = 501

      real(dp), private, parameter  :: tdim_real = 501.0
      
      
      real(dp), private, parameter :: t_del = 5._dp/(tdim_real-1._dp)
      real(dp), private, parameter :: t_fac = (tdim_real-1._dp)/5._dp

      integer, private :: ii, jj
      real(dp), private, save :: d_table(0:tdim,ngast-1)
      real(dp), private, save :: x_table(0:tdim,ngast-1)
      real(dp), private, save :: o2_table(tdim)
      real(dp), private, dimension(12,ngast-1), save :: a_schu, b_schu






      data ((a_schu(jj,ii),jj=1,12),ii=1,ngast-1) /                             &

       1.13402e-01,1.00088e-20,3.48747e-01,2.76282e-20,3.47322e-01,1.01267e-19, &
       1.67351e-01,5.63588e-19,2.31433e-02,1.68267e-18,0.00000e+00,0.00000e+00, &

       2.55268e-03,1.64489e-21,1.85483e-01,2.03591e-21,2.60603e-01,4.62276e-21, &
       2.50337e-01,1.45106e-20,1.92340e-01,7.57381e-20,1.06363e-01,7.89634e-19, &

       4.21594e-03,8.46639e-22,8.91886e-02,1.12935e-21,2.21334e-01,1.67868e-21, &
       2.84446e-01,3.94782e-21,2.33442e-01,1.91554e-20,1.63433e-01,2.25346e-19, &

       3.93529e-03,6.79660e-22,4.46906e-02,9.00358e-22,1.33060e-01,1.55952e-21, &
       3.25506e-01,3.43763e-21,2.79405e-01,1.62086e-20,2.10316e-01,1.53883e-19, &

       2.60939e-03,2.33791e-22,2.08101e-02,3.21734e-22,1.67186e-01,5.77191e-22, &
       2.80694e-01,1.33362e-21,3.26867e-01,6.10533e-21,1.96539e-01,7.83142e-20, &

       9.33711e-03,1.32897e-22,3.63980e-02,1.78786e-22,1.46182e-01,3.38285e-22, &
       3.81762e-01,8.93773e-22,2.58549e-01,4.28115e-21,1.64773e-01,4.67537e-20, &

       9.51799e-03,1.00252e-22,3.26320e-02,1.33766e-22,1.45962e-01,2.64831e-22, &
       4.49823e-01,6.42879e-22,2.14207e-01,3.19594e-21,1.45616e-01,2.77182e-20, &

       7.87331e-03,3.38291e-23,6.91451e-02,4.77708e-23,1.29786e-01,8.30805e-23, &
       3.05103e-01,2.36167e-22,3.35007e-01,8.59109e-22,1.49766e-01,9.63516e-21, &

       6.92175e-02,1.56323e-23,1.44403e-01,3.03795e-23,2.94489e-01,1.13219e-22, &
       3.34773e-01,3.48121e-22,9.73632e-02,2.10693e-21,5.94308e-02,1.26195e-20, &

       1.47873e-01,8.62033e-24,3.15881e-01,3.51859e-23,4.08077e-01,1.90524e-22, &
       8.08029e-02,9.93062e-22,3.90399e-02,6.38738e-21,8.13330e-03,9.93644e-22, &

       1.50269e-01,1.02621e-23,2.39823e-01,3.48120e-23,3.56408e-01,1.69494e-22, &
       1.61277e-01,6.59294e-22,8.89713e-02,2.94571e-21,3.25063e-03,1.25548e-20, &

       2.55746e-01,8.49877e-24,2.94733e-01,2.06878e-23,2.86382e-01,9.30992e-23, &
       1.21011e-01,3.66239e-22,4.21105e-02,1.75700e-21,0.00000e+00,0.00000e+00, &

       5.40111e-01,7.36085e-24,2.93263e-01,2.46742e-23,1.63417e-01,1.37832e-22, &
       3.23781e-03,2.15052e-21,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00, &

       8.18514e-01,7.17937e-24,1.82262e-01,4.17496e-23,0.00000e+00,0.00000e+00, &
       0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00, &

       8.73680e-01,7.13444e-24,1.25583e-01,2.77819e-23,0.00000e+00,0.00000e+00, &
       0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00, &

       3.32476e-04,7.00362e-24,9.89000e-01,6.99600e-24,0.00000e+00,0.00000e+00, &
       0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00 /
      data ((b_schu(jj,ii),jj=1,12),ii=1,ngast-1) /                             &

       1.07382e-21,9.95029e-21,7.19430e-21,2.48960e-20,2.53735e-20,7.54467e-20, &
       4.48987e-20,2.79981e-19,9.72535e-20,9.29745e-19,2.30892e-20,4.08009e-17, &

       3.16903e-22,1.98251e-21,5.87326e-22,3.44057e-21,2.53094e-21,8.81484e-21, &
       8.82299e-21,4.17179e-20,2.64703e-20,2.43792e-19,8.73831e-20,1.46371e-18, &

       1.64421e-23,9.26011e-22,2.73137e-22,1.33640e-21,9.79188e-22,2.99706e-21, &
       3.37768e-21,1.39438e-20,1.47898e-20,1.04322e-19,4.08014e-20,6.31023e-19, &

       8.68729e-24,7.31056e-22,8.78313e-23,1.07173e-21,8.28170e-22,2.54986e-21, &
       2.57643e-21,9.42698e-21,9.92377e-21,5.21402e-20,3.34301e-20,2.91785e-19, &

       1.20679e-24,2.44092e-22,2.64326e-23,4.03998e-22,2.53514e-22,8.53166e-22, &
       1.29834e-21,3.74482e-21,5.12103e-21,2.65798e-20,2.10948e-20,2.35315e-19, &

       2.79656e-24,1.40820e-22,3.60824e-23,2.69510e-22,4.02850e-22,8.83735e-22, &
       1.77198e-21,6.60221e-21,9.60992e-21,8.13558e-20,4.95591e-21,1.22858e-17, &

       2.36959e-24,1.07535e-22,2.83333e-23,2.16789e-22,3.35242e-22,6.42753e-22, &
       1.26395e-21,5.43183e-21,4.88083e-21,5.42670e-20,3.27481e-21,1.58264e-17, &

       8.65018e-25,3.70310e-23,1.04351e-23,6.43574e-23,1.17431e-22,2.70904e-22, &
       4.88705e-22,1.65505e-21,2.19776e-21,2.71172e-20,2.65257e-21,2.13945e-17, &

       9.63263e-25,1.54249e-23,4.78065e-24,2.97642e-23,6.40637e-23,1.46464e-22, &
       1.82634e-22,7.12786e-22,1.64805e-21,2.37376e-17,9.33059e-22,1.13741e-20, &

       1.08414e-24,8.37560e-24,9.15550e-24,2.99295e-23,9.38405e-23,1.95845e-22, &
       2.84356e-22,3.39699e-21,1.94524e-22,2.72227e-19,1.18924e-21,3.20246e-17, &

       1.52817e-24,1.01885e-23,1.22946e-23,4.16517e-23,9.01287e-23,2.34869e-22, &
       1.93510e-22,1.44956e-21,1.81051e-22,5.17773e-21,9.82059e-22,6.22768e-17, &

       2.12813e-24,8.48035e-24,5.23338e-24,1.93052e-23,1.99464e-23,7.48997e-23, &
       4.96642e-22,6.15691e-17,4.47504e-23,2.76004e-22,8.26788e-23,1.65278e-21, &

       3.81336e-24,7.32307e-24,5.60549e-24,2.04651e-23,3.36883e-22,6.15708e-17, &
       2.09877e-23,1.07474e-22,9.13562e-24,8.41252e-22,0.00000e+00,0.00000e+00, &

       5.75373e-24,7.15986e-24,5.90031e-24,3.05375e-23,2.97196e-22,8.92000e-17, &
       8.55920e-24,1.66709e-17,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00, &

       6.21281e-24,7.13108e-24,3.30780e-24,2.61196e-23,1.30783e-22,9.42550e-17, &
       2.69241e-24,1.46500e-17,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00, &

       6.81118e-24,6.98767e-24,7.55667e-25,2.75124e-23,1.94044e-22,1.45019e-16, &
       1.92236e-24,3.73223e-17,0.00000e+00,0.00000e+00,0.00000e+00,0.00000e+00 /




      integer, public,  parameter :: nwint = 105

      real(dp), public, save  :: wlint(nwint)
      real(dp), private, save :: xso2int(nwint)
      real(dp), private, save :: wlla(nla)
      real(dp), private, save :: wlgast(ngast)




      data (wlint(ii),ii=1,105)/  &
            0.0000E+00, 0.1166E+03, 0.1167E+03, 0.1173E+03, 0.1179E+03,  &
            0.1187E+03, 0.1194E+03, 0.1202E+03, 0.1208E+03, 0.1214E+03,  &
            0.1219E+03, 0.1223E+03, 0.1231E+03, 0.1238E+03, 0.1246E+03,  &
            0.1254E+03, 0.1262E+03, 0.1270E+03, 0.1278E+03, 0.1286E+03,  &
            0.1294E+03, 0.1303E+03, 0.1311E+03, 0.1320E+03, 0.1329E+03,  &
            0.1338E+03, 0.1346E+03, 0.1356E+03, 0.1365E+03, 0.1374E+03,  &
            0.1384E+03, 0.1399E+03, 0.1418E+03, 0.1439E+03, 0.1459E+03,  &
            0.1481E+03, 0.1504E+03, 0.1526E+03, 0.1550E+03, 0.1574E+03,  &
            0.1600E+03, 0.1626E+03, 0.1653E+03, 0.1681E+03, 0.1709E+03,  &
            0.1731E+03, 0.1746E+03, 0.1754E+03, 0.1770E+03, 0.1786E+03,  &
            0.1802E+03, 0.1818E+03, 0.1835E+03, 0.1852E+03, 0.1869E+03,  &
            0.1887E+03, 0.1905E+03, 0.1923E+03, 0.1942E+03, 0.1961E+03,  &
            0.1980E+03, 0.2000E+03, 0.2020E+03, 0.2041E+03, 0.2050E+03,  &
            0.2060E+03, 0.2070E+03, 0.2080E+03, 0.2090E+03, 0.2100E+03,  &
            0.2110E+03, 0.2120E+03, 0.2130E+03, 0.2140E+03, 0.2150E+03,  &
            0.2160E+03, 0.2170E+03, 0.2180E+03, 0.2190E+03, 0.2200E+03,  &
            0.2210E+03, 0.2220E+03, 0.2230E+03, 0.2240E+03, 0.2250E+03,  &
            0.2260E+03, 0.2270E+03, 0.2280E+03, 0.2290E+03, 0.2300E+03,  &
            0.2310E+03, 0.2320E+03, 0.2330E+03, 0.2340E+03, 0.2350E+03,  &
            0.2360E+03, 0.2370E+03, 0.2380E+03, 0.2390E+03, 0.2400E+03,  &
            0.2410E+03, 0.2420E+03, 0.2430E+03, 0.2430E+03, 0.1000E+39 /




      data (xso2int(ii),ii=1,105)/  &
            0.0000E+00, 0.1000E-19, 0.6350E-18, 0.7525E-18, 0.1425E-18,  &
            0.2025E-18, 0.2413E-17, 0.6400E-17, 0.5251E-17, 0.7329E-18,  &
            0.3445E-18, 0.3425E-18, 0.3925E-18, 0.8918E-17, 0.9197E-17,  &
            0.6625E-18, 0.2700E-18, 0.1575E-18, 0.3240E-18, 0.4990E-18,  &
            0.4875E-18, 0.5525E-18, 0.1068E-17, 0.1850E-17, 0.2275E-17,  &
            0.3425E-17, 0.5890E-17, 0.8365E-17, 0.1090E-16, 0.1275E-16,  &
            0.1340E-16, 0.1380E-16, 0.1440E-16, 0.1445E-16, 0.1350E-16,  &
            0.1220E-16, 0.1071E-16, 0.9075E-17, 0.7410E-17, 0.5775E-17,  &
            0.4210E-17, 0.2765E-17, 0.1655E-17, 0.9760E-18, 0.5900E-18,  &
            0.3660E-18, 0.2061E-18, 0.3889E-19, 0.4092E-20, 0.2273E-20,  &
            0.1256E-20, 0.6902E-21, 0.3750E-21, 0.2097E-21, 0.1226E-21,  &
            0.6508E-22, 0.4579E-22, 0.3220E-22, 0.2507E-22, 0.1941E-22,  &
            0.1609E-22, 0.1319E-22, 0.9907E-23, 0.7419E-23, 0.7240E-23,  &
            0.7090E-23, 0.6955E-23, 0.6770E-23, 0.6595E-23, 0.6375E-23,  &
            0.6145E-23, 0.5970E-23, 0.5805E-23, 0.5655E-23, 0.5470E-23,  &
            0.5240E-23, 0.5005E-23, 0.4760E-23, 0.4550E-23, 0.4360E-23,  &
            0.4175E-23, 0.3990E-23, 0.3780E-23, 0.3560E-23, 0.3330E-23,  &
            0.3095E-23, 0.2875E-23, 0.2875E-23, 0.2875E-23, 0.2700E-23,  &
            0.2530E-23, 0.2340E-23, 0.2175E-23, 0.2020E-23, 0.1860E-23,  &
            0.1705E-23, 0.1555E-23, 0.1410E-23, 0.1280E-23, 0.1160E-23,  &
            0.1055E-23, 0.9550E-24, 0.4500E-24, 0.0000E+00, 0.0000E+00 /




      data (wlla(ii),ii=1,2)/ 121.4_dp, 121.9_dp /




      data (wlgast(ii),ii=1,17)/  &
            0.1754E+03_dp, 0.1770E+03_dp, 0.1786E+03_dp, 0.1802E+03_dp, 0.1818E+03_dp,  &
            0.1835E+03_dp, 0.1852E+03_dp, 0.1869E+03_dp, 0.1887E+03_dp, 0.1905E+03_dp,  &
            0.1923E+03_dp, 0.1942E+03_dp, 0.1961E+03_dp, 0.1980E+03_dp, 0.2000E+03_dp,  &
            0.2020E+03_dp, 0.2041E+03_dp /




      real(dp) :: op_cb1w(nw-1), om_cb1w(nw-1), g_cb1w(nw-1)
      real(dp) :: op_cb2w(nw-1), om_cb2w(nw-1), g_cb2w(nw-1)
      real(dp) :: op_oc1w(nw-1), g_oc1w(nw-1)
      real(dp) :: op_oc2w(nw-1), om_oc2w(nw-1), g_oc2w(nw-1)
      real(dp) :: op_antw(nw-1), g_antw(nw-1)
      real(dp) :: op_so4w(nw-1)
      real(dp) :: op_salw(nw-1)

      CONTAINS

      subroutine photoin( chem_opt, nz, zen, o3toms, esfact,     &
                          o3top, o2top, albedo, z, tlev,         & 
                          tlay, airlev, rh, xlwc, o3,            &
                          acb1, acb2, aoc1, aoc2, aant,          & 
                          aso4, asal, tauaer300, tauaer400, tauaer600, &
                          tauaer999, waer300, waer400, waer600,  &
                          waer999, gaer300, gaer400, gaer600,    &
                          gaer999, aer_ra_feedback, prate, adjcoe, radfld )

      use module_wave_data, only : nw, tuv_jmax, deltaw, sflx,   &
                                   c20, c40, c60, c80, sq



      implicit none



      integer, intent(in)    :: chem_opt      
      integer, intent(in)    :: nz      
      real(dp), intent(in)    :: zen     
      real(dp), intent(in)    :: o3toms  
      real(dp), intent(in)    :: esfact  
      real(dp), intent(in)    :: o3top   
      real(dp), intent(in)    :: o2top  
      real(dp), intent(in)    :: albedo(nw-1)
      real(dp), intent(in)    :: z(nz)
      real(dp), intent(in)    :: tlev(nz)
      real(dp), intent(in)    :: tlay(nz-1)
      real(dp), intent(in)    :: airlev(nz)
      real(dp), intent(in)    :: rh(nz)
      real(dp), intent(in)    :: xlwc(nz)
      real(dp), intent(in)    :: o3(nz)
      real(dp), intent(in)    :: acb1(nz)
      real(dp), intent(in)    :: acb2(nz)
      real(dp), intent(in)    :: aoc1(nz)
      real(dp), intent(in)    :: aoc2(nz)
      real(dp), intent(in)    :: aant(nz)
      real(dp), intent(in)    :: aso4(nz)
      real(dp), intent(in)    :: asal(nz)

      real(dp), intent(in)    :: tauaer300(nz-1)
      real(dp), intent(in)    :: tauaer400(nz-1)
      real(dp), intent(in)    :: tauaer600(nz-1)
      real(dp), intent(in)    :: tauaer999(nz-1)
      real(dp), intent(in)    :: waer300(nz-1)
      real(dp), intent(in)    :: waer400(nz-1)
      real(dp), intent(in)    :: waer600(nz-1)
      real(dp), intent(in)    :: waer999(nz-1)
      real(dp), intent(in)    :: gaer300(nz-1)
      real(dp), intent(in)    :: gaer400(nz-1)
      real(dp), intent(in)    :: gaer600(nz-1)
      real(dp), intent(in)    :: gaer999(nz-1)
      INTEGER, INTENT(IN)     :: aer_ra_feedback

      real(dp), intent(out)   :: prate(nz,tuv_jmax)
      real(dp), intent(out)   :: adjcoe(nz,tuv_jmax)
      real(dp), intent(out)   :: radfld(nz,nw-1)



      integer :: i, j, k, iw, n
      real(dp) :: factor, delzint
      character(len=8 ) :: cdate(4)
      character(len=10) :: ctime(4)



      integer  :: iz, izi
      real(dp) :: colinc(nz)
      real(dp) :: vcol(nz)
      real(dp) :: scol(nz)
      real(dp) :: to3(nz)



      real(dp), dimension(nz,tuv_jmax) :: adjcoe1, adjcoe2




      integer  :: nid(0:nz-1)
      real(dp) :: dsdh(0:nz-1,nz-1)



      real(dp), dimension(nw-1) :: etf, delw, xsec



      real(dp) :: xso2(nw-1,nz)



      real(dp), dimension(nz-1,nw-1) :: dtrl, dto3, dto2
      real(dp), dimension(nz-1,nw-1) :: dtcld, omcld, gcld
      real(dp), dimension(nz-1,nw-1) :: dtcb1, omcb1, gcb1
      real(dp), dimension(nz-1,nw-1) :: dtcb2, omcb2, gcb2
      real(dp), dimension(nz-1,nw-1) :: dtoc1, omoc1, goc1
      real(dp), dimension(nz-1,nw-1) :: dtoc2, omoc2, goc2
      real(dp), dimension(nz-1,nw-1) :: dtant, omant, gant
      real(dp), dimension(nz-1,nw-1) :: dtso4, omso4, gso4
      real(dp), dimension(nz-1,nw-1) :: dtsal, omsal, gsal

      real(dp), dimension(nz-1,nw-1) :: dtaer, omaer, gaer



      real(dp), dimension(nz,nw-1) :: radxx



      integer :: m



      integer :: iyear, imonth, iday
      real(dp) :: dtime, ut0



      etf(1:nw-1) = sflx(1:nw-1) * esfact 



      call setair( nz, z, airlev, dtrl, colinc, o2top )



      call setozo( nz, z, tlay, dto3, to3, o3, airlev, o3top, o3toms )



      call setcld( nz, z, xlwc, dtcld, omcld, gcld )




















      dtaer(:,:) = 0._dp
      omaer(:,:) = 0._dp
       gaer(:,:) = 0._dp
      if(aer_ra_feedback == 1) then 
       call aer_wrf2ftuv(nz, z, tauaer300, tauaer400,    &
                         tauaer600, tauaer999, waer300, &
                         waer400, waer600, waer999,     &
                         gaer300, gaer400, gaer600,     &
                         gaer999, dtaer, omaer, gaer    )
      endif 







      call pchem( chem_opt, nz, tlev, airlev )




       call spheres( nz, z, zen, dsdh, nid )

       call airmas( nz, z, zen, dsdh, nid, colinc, vcol, scol )




      if( zen <= 20.5_dp ) then
         call setz( nz, to3, tlev, c20, 20, adjcoe )
      else if( zen > 20.5_dp .and. zen < 40.5_dp ) then
         call setz( nz, to3, tlev, c20, 20, adjcoe1 )
         call setz( nz, to3, tlev, c40, 40, adjcoe2 )
         factor = (zen - 20.5_dp) / 20._dp
         adjcoe(:nz,:tuv_jmax) = adjcoe1(:nz,:tuv_jmax) &
                               + (adjcoe2(:nz,:tuv_jmax) - adjcoe1(:nz,:tuv_jmax)) * factor
      else if( zen >= 40.5_dp .and. zen < 60.5_dp ) then
         call setz( nz, to3, tlev, c40, 40, adjcoe1 )
         call setz( nz, to3, tlev, c60, 60, adjcoe2 )
         factor = (zen - 40.5_dp) / 20._dp
         adjcoe(:nz,:tuv_jmax) = adjcoe1(:nz,:tuv_jmax) &
                               + (adjcoe2(:nz,:tuv_jmax) - adjcoe1(:nz,:tuv_jmax)) * factor
      else if( zen >= 60.5_dp .and. zen < 80._dp ) then
         call setz( nz, to3, tlev, c60, 60, adjcoe1 )
         call setz( nz, to3, tlev, c80, 80, adjcoe2 )
         factor = (zen - 60.5_dp) / 19.5_dp
         adjcoe(:nz,:tuv_jmax) = adjcoe1(:nz,:tuv_jmax) &
                               + (adjcoe2(:nz,:tuv_jmax) - adjcoe1(:nz,:tuv_jmax)) * factor
      else if( zen >= 80. ) then
         call setz( nz, to3, tlev, c80, 80, adjcoe )
      end if




       call set_o2_xsect( nz, z, colinc, vcol, scol, dto2, xso2 )
       sq(1:nw-1,1:nz,1) = xso2(1:nw-1,1:nz)




      delw(1:nw-1) = deltaw(1:nw-1) * etf(1:nw-1)




      call rtlink( nz,                  &
                   nw,                  &
                   zen,                 &
                   albedo,              & 
                   z,                   &
                   nid,                 &
                   dsdh,                & 
                   dtrl,                & 
                   dto3,                & 
                   dto2,                & 
                   dtcld, omcld, gcld,  & 







                   dtaer, omaer, gaer,  &  
                   radfld )                           



         radfld(:,:) = max( radfld(:,:),0._dp )
         delzint = (z(nz-2) - z(nz-3))/(z(nz-1) - z(nz-2))
         do iw = 1,nw-1
            radfld(1,iw) = radfld(2,iw) + (radfld(2,iw)-radfld(3,iw))*delzint
            radfld(1,iw) = max(radfld(1,iw),radfld(2,iw))
         end do








      do m = 1,tuv_jmax
        do iz = 1,nz
          izi = nz - iz + 1
          xsec(:nw-1) = sq(:nw-1,iz,m) * delw(:nw-1)
          prate(iz,m) = dot_product( radfld(izi,:nw-1), xsec(:nw-1) )
        end do
        prate(1:nz,m) = prate(1:nz,m) * adjcoe(1:nz,m)
      end do

      end subroutine photoin






      subroutine aer_wrf2ftuv( nzlev, z, tauaer300, tauaer400,     &
                               tauaer600, tauaer999,               &
                               waer300, waer400, waer600, waer999, &
                               gaer300, gaer400, gaer600, gaer999, &
                               dtaer, omaer, gaer )
       use module_wave_data, only : nw, wc
















      implicit none




      integer, intent(in)   :: nzlev
      real(dp), intent(in)  :: z(nzlev)
      real(dp), intent(in)  :: tauaer300(nzlev-1), tauaer400(nzlev-1),    &
                               tauaer600(nzlev-1), tauaer999(nzlev-1)
      real(dp), intent(in)  :: waer300(nzlev-1), waer400(nzlev-1),        &
                               waer600(nzlev-1), waer999(nzlev-1)
      real(dp), intent(in)  :: gaer300(nzlev-1), gaer400(nzlev-1),        &
                               gaer600(nzlev-1), gaer999(nzlev-1)

      real(dp), intent(out) :: dtaer(nzlev-1,nw-1),                   &
                               omaer(nzlev-1,nw-1), gaer(nzlev-1,nw-1)


      integer     :: k, wn, nloop
      real(dp)    :: ang, slope
      real(dp), parameter :: thresh = 1.e-9_dp

      ang   = 0._dp
      slope = 0._dp


      do wn = 1,nw-1       
       do k = 1,nzlev-1    


        if(tauaer300(k) .gt. thresh .and. tauaer999(k) .gt. thresh) then
         ang = log(tauaer300(k)/tauaer999(k))/log(0.999_dp/0.3_dp)
         dtaer(k,wn) = tauaer400(k)*(0.4_dp/(wc(wn)*1.e-3_dp))**ang



         slope = (waer600(k)-waer400(k))/0.2_dp
         omaer(k,wn) = slope*((wc(wn)*1.e-3_dp)-0.6_dp)+waer600(k)
         if(omaer(k,wn) .lt. 0.4_dp) omaer(k,wn)=0.4_dp
         if(omaer(k,wn) .ge. 1.0_dp) omaer(k,wn)=1.0_dp


         slope = (gaer600(k)-gaer400(k))/0.2_dp
         gaer(k,wn) = slope*((wc(wn)*1.e-3_dp)-0.6_dp)+gaer600(k)
         if(gaer(k,wn) .lt. 0.5_dp) gaer(k,wn) = 0.5_dp
         if(gaer(k,wn) .ge. 1.0_dp) gaer(k,wn) = 1.0_dp
        endif
       enddo  
      enddo   

      end subroutine aer_wrf2ftuv



      subroutine setaer( nzlev, z, airden, rh, acb1,    &
                         acb2, aoc1, aoc2, aant, aso4,  &
                         asal,                          &
                         op_cb1, om_cb1, g_cb1,         &
                         op_cb2, om_cb2, g_cb2,         &
                         op_oc1, om_oc1, g_oc1,         &
                         op_oc2, om_oc2, g_oc2,         &
                         op_ant, om_ant, g_ant,         &
                         op_so4, om_so4, g_so4,         &
                         op_sal, om_sal, g_sal )

      use module_wave_data, only : nw, wc






















      implicit none




      integer, intent(in)   :: nzlev
      real(dp), intent(in)  :: z(nzlev)
      real(dp), intent(in)  :: airden(nzlev)
      real(dp), intent(in)  :: rh(nzlev)
      real(dp), intent(in)  :: acb1(nzlev)
      real(dp), intent(in)  :: acb2(nzlev)
      real(dp), intent(in)  :: aoc1(nzlev)
      real(dp), intent(in)  :: aoc2(nzlev)
      real(dp), intent(in)  :: aant(nzlev)
      real(dp), intent(in)  :: aso4(nzlev)
      real(dp), intent(in)  :: asal(nzlev)
      real(dp), intent(out) :: op_cb1(nzlev-1,nw-1), om_cb1(nzlev-1,nw-1), g_cb1(nzlev-1,nw-1)
      real(dp), intent(out) :: op_cb2(nzlev-1,nw-1), om_cb2(nzlev-1,nw-1), g_cb2(nzlev-1,nw-1)
      real(dp), intent(out) :: op_oc1(nzlev-1,nw-1), om_oc1(nzlev-1,nw-1), g_oc1(nzlev-1,nw-1)
      real(dp), intent(out) :: op_oc2(nzlev-1,nw-1), om_oc2(nzlev-1,nw-1), g_oc2(nzlev-1,nw-1)
      real(dp), intent(out) :: op_ant(nzlev-1,nw-1), om_ant(nzlev-1,nw-1), g_ant(nzlev-1,nw-1)
      real(dp), intent(out) :: op_so4(nzlev-1,nw-1), om_so4(nzlev-1,nw-1), g_so4(nzlev-1,nw-1)
      real(dp), intent(out) :: op_sal(nzlev-1,nw-1), om_sal(nzlev-1,nw-1), g_sal(nzlev-1,nw-1)




      integer     :: k, wn, nz
      real(dp)    :: wcen
      real(dp)    :: dz(nzlev)
      real(dp)    :: sig_oc2(nzlev-1), sig_cb2(nzlev-1)
      real(dp)    :: opt_oc1(nzlev-1), opt_oc2(nzlev-1), opt_cb1(nzlev-1), opt_cb2(nzlev-1)
      real(dp)    :: opt_ant(nzlev-1), opt_so4(nzlev-1), opt_sal(nzlev-1)  
      real(dp)    :: rw(nzlev-1), num(nzlev-1), wght(nzlev-1)





      real(dp), parameter :: rm_cb1 = 0.035_dp
      real(dp), parameter :: rm_cb2 = 0.035_dp
      real(dp), parameter :: rm_oc1 = 0.10_dp
      real(dp), parameter :: rm_oc2 = 0.10_dp
      real(dp), parameter :: rm_ant = 0.15_dp
      real(dp), parameter :: rm_so4 = 0.24_dp
      real(dp), parameter :: rm_sal = 1.30_dp

      real(dp), parameter :: de_cb1 = 1.0_dp
      real(dp), parameter :: de_cb2 = 1.0_dp
      real(dp), parameter :: de_oc1 = 1.5_dp 
      real(dp), parameter :: de_oc2 = 1.5_dp 
      real(dp), parameter :: de_ant = 1.7_dp  
      real(dp), parameter :: de_so4 = 1.7_dp
      real(dp), parameter :: de_sal = 2.2_dp

      real(dp), parameter :: ml_cb1 = 12._dp
      real(dp), parameter :: ml_cb2 = 12._dp
      real(dp), parameter :: ml_oc1 = 48._dp 
      real(dp), parameter :: ml_oc2 = 48._dp 
      real(dp), parameter :: ml_ant = 90._dp  
      real(dp), parameter :: ml_so4 = 96._dp
      real(dp), parameter :: ml_sal = 22._dp

      real(dp), parameter :: qx_cb1 = 2.138_dp
      real(dp), parameter :: qx_cb2 = 2.138_dp
      real(dp), parameter :: qw_cb1 = 0.510_dp
      real(dp), parameter :: qw_cb2 = 0.510_dp
      real(dp), parameter :: sx_cb2 = 0.120_dp
      real(dp), parameter :: sw_cb2 = 1.000_dp

      real(dp), parameter :: qx_oc1 = 0.675_dp
      real(dp), parameter :: qx_oc2 = 0.675_dp
      real(dp), parameter :: qw_oc1 = 1.118_dp
      real(dp), parameter :: qw_oc2 = 1.118_dp
      real(dp), parameter :: sx_oc2 = 0.980_dp
      real(dp), parameter :: sw_oc2 = 1.000_dp

      real(dp), parameter :: qx_ant = 0.726_dp
      real(dp), parameter :: qw_ant = 1.770_dp

      real(dp), parameter :: qx_so4 = 2.830_dp
      real(dp), parameter :: qw_so4 = 3.605_dp

      real(dp), parameter :: qx_sal = 2.600_dp
      real(dp), parameter :: qw_sal = 2.575_dp

      real(dp), parameter :: pi = 3.1416_dp




      real(dp), parameter :: a1_cb2 = 1.000_dp
      real(dp), parameter :: a2_cb2 =-0.829_dp
      real(dp), parameter :: a3_cb2 = 1.350_dp

      real(dp), parameter :: a1_oc2 = 1.005_dp
      real(dp), parameter :: a2_oc2 =-0.0666_dp
      real(dp), parameter :: a3_oc2 = 0.7608_dp

      real(dp), parameter :: a1_ant = 1.000_dp
      real(dp), parameter :: a2_ant = 0.493_dp
      real(dp), parameter :: a3_ant = 0.527_dp

      real(dp), parameter :: a1_so4 = 1.000_dp
      real(dp), parameter :: a2_so4 = 0.493_dp
      real(dp), parameter :: a3_so4 = 0.528_dp

      real(dp), parameter :: a1_sal = 1.011_dp
      real(dp), parameter :: a2_sal = 0.457_dp
      real(dp), parameter :: a3_sal = 1.102_dp

      nz = nzlev - 1



      dz(1:nz) = abs( z(2:nz+1) - z(1:nz) ) * km2cm * pi * 1.e-8_dp



         rw(:nz) = rm_cb1
         call cnum( acb1, airden, de_cb1, ml_cb1, rm_cb1, &
                    rw, num, wght, nzlev )
      
         opt_cb1(:nz) = qx_cb1*wght(:nz) + (1._dp - wght(:nz))*qw_cb1     
         opt_cb1(:nz) = opt_cb1(:nz)*num(:nz)*rw(:nz)*rw(:nz)*dz(:nz)



         rw(:nz) = rm_cb2*(a1_cb2 + rh(:nz)*(a2_cb2 + a3_cb2*rh(:nz)))
         call cnum( acb2, airden, de_cb2, ml_cb2, rm_cb2, &
                    rw, num, wght, nzlev )
      
         opt_cb2(:nz) = qx_cb2*wght(:nz) + (1._dp - wght(:nz))*qw_cb2     
         opt_cb2(:nz) = opt_cb2(:nz)*num(:nz)*rw(:nz)*rw(:nz)*dz(:nz)
         sig_cb2(:nz) = sx_cb2*wght(:nz) + (1._dp - wght(:nz))*sw_cb2   



         rw(:nz) = rm_oc1
         call cnum( aoc1, airden, de_oc1, ml_oc1, rm_oc1, &
                    rw, num, wght, nzlev )
      
         opt_oc1(:nz) = qx_oc1*wght(:nz) + (1._dp - wght(:nz))*qw_oc1     
         opt_oc1(:nz) = opt_oc1(:nz)*num(:nz)*rw(:nz)*rw(:nz)*dz(:nz)



         rw(:nz) = rm_oc2*(a1_oc2 + rh(:nz)*(a2_oc2 + a3_oc2*rh(:nz)))
         call cnum( aoc2, airden, de_oc2, ml_oc2, rm_oc2, &
                    rw, num, wght, nzlev )
      
         opt_oc2(:nz) = qx_oc2*wght(:nz) + (1._dp - wght(:nz))*qw_oc2     
         opt_oc2(:nz) = opt_oc2(:nz)*num(:nz)*rw(:nz)*rw(:nz)*dz(:nz)
         sig_oc2(:nz) = sx_oc2*wght(:nz) + (1._dp - wght(:nz))*sw_oc2   



         rw(:nz) = rm_ant*(a1_ant + rh(:nz)*(a2_ant + a3_ant*rh(:nz)))
         call cnum( aant, airden, de_ant, ml_ant, rm_ant, &
                    rw, num, wght, nzlev )
      
         opt_ant(:nz) = qx_ant*wght(:nz) + (1._dp - wght(:nz))*qw_ant     
         opt_ant(:nz) = opt_ant(:nz)*num(:nz)*rw(:nz)*rw(:nz)*dz(:nz)



         rw(:nz) = rm_so4*(a1_so4 + rh(:nz)*(a2_so4 + a3_so4*rh(:nz)))
         call cnum( aso4, airden, de_so4, ml_so4, rm_so4, &
                    rw, num, wght, nzlev )
      
         opt_so4(:nz) = qx_so4*wght(:nz) + (1._dp - wght(:nz))*qw_so4     
         opt_so4(:nz) = opt_so4(:nz)*num(:nz)*rw(:nz)*rw(:nz)*dz(:nz)



         rw(:nz) = rm_sal*(a1_sal + rh(:nz)*(a2_sal + a3_sal*rh(:nz)))
         call cnum( asal, airden, de_sal, ml_sal, rm_sal, &
                    rw, num, wght, nzlev )
      
         opt_sal(:nz) = qx_sal*wght(:nz) + (1._dp - wght(:nz))*qw_sal     
         opt_sal(:nz) = opt_sal(:nz)*num(:nz)*rw(:nz)*rw(:nz)*dz(:nz)

WAVE_LENGTH_LOOP : &
         do wn = 1,nw-1
            wcen = wc(wn)
            op_cb1(:nz,wn) = opt_cb1(:nz)*op_cb1w(wn)
            om_cb1(:nz,wn) = om_cb1w(wn)
             g_cb1(:nz,wn) = g_cb1w(wn)

            op_cb2(:nz,wn) = opt_cb2(:nz)*op_cb2w(wn)
            om_cb2(:nz,wn) = sig_cb2(:nz)*om_cb2w(wn)
             g_cb2(:nz,wn) = g_cb2w(wn)

            op_oc1(:nz,wn) = opt_oc1(:nz)*op_oc1w(wn)
            om_oc1(:nz,wn) = 0.98_dp
             g_oc1(:nz,wn) = g_oc1w(wn)

            op_oc2(:nz,wn) = opt_oc2(:nz)*op_oc2w(wn)
            om_oc2(:nz,wn) = sig_oc2(:nz)
             g_oc2(:nz,wn) = g_oc2w(wn)

            op_ant(:nz,wn) = opt_ant(:nz)*op_antw(wn)
            om_ant(:nz,wn) = 1.00_dp
             g_ant(:nz,wn) = g_antw(wn)

            op_so4(:nz,wn) = opt_so4(:nz)*op_so4w(wn)
            om_so4(:nz,wn) = 1.00_dp
             g_so4(:nz,wn) = 0.75_dp

            op_sal(:nz,wn) = opt_sal(:nz)*op_salw(wn)
            om_sal(:nz,wn) = 1.00_dp
             g_sal(:nz,wn) = 0.77_dp
         enddo WAVE_LENGTH_LOOP

      end subroutine setaer

      subroutine aer_init




      use module_wave_data, only : nw, wc

      implicit none




      real(dp), parameter :: d1_cb = 1.61892_dp
      real(dp), parameter :: d2_cb =-0.0006853_dp
      real(dp), parameter :: d3_cb =-1.07806e-6_dp

      real(dp), parameter :: o1_cb = 1.29185_dp
      real(dp), parameter :: o2_cb = 7.6863e-5_dp
      real(dp), parameter :: o3_cb =-1.28567e-6_dp

      real(dp), parameter :: g1_cb = 2.87326_dp
      real(dp), parameter :: g2_cb =-0.004482_dp
      real(dp), parameter :: g3_cb = 1.50361e-6_dp

      real(dp), parameter :: d1_oc = 10.5098_dp
      real(dp), parameter :: d2_oc =-0.0301455_dp
      real(dp), parameter :: d3_oc = 2.22837e-5_dp

      real(dp), parameter :: g1_oc = 1.67125_dp
      real(dp), parameter :: g2_oc =-0.0008341_dp
      real(dp), parameter :: g3_oc =-1.18801e-6_dp

      real(dp), parameter :: d1_ant = 9.44794_dp
      real(dp), parameter :: d2_ant =-0.0268029_dp
      real(dp), parameter :: d3_ant = 1.97106e-5_dp

      real(dp), parameter :: g1_ant = 1.34742_dp
      real(dp), parameter :: g2_ant =-1.85850e-5_dp
      real(dp), parameter :: g3_ant =-1.54297e-6_dp

      real(dp), parameter :: d1_so4 = 0.942351_dp
      real(dp), parameter :: d2_so4 = 0.00220168_dp
      real(dp), parameter :: d3_so4 =-3.9622e-6_dp

      real(dp), parameter :: d1_sal = 1.0_dp
      real(dp), parameter :: d2_sal = 0.0_dp
      real(dp), parameter :: d3_sal = 0.0_dp

      op_cb1w(:nw-1) = d1_cb + wc(:nw-1)*(d2_cb + d3_cb*wc(:nw-1))
      om_cb1w(:nw-1) = .12_dp*(o1_cb + wc(:nw-1)*(o2_cb + o3_cb*wc(:nw-1)))
      g_cb1w(:nw-1)  = .12_dp*(g1_cb + wc(:nw-1)*(g2_cb + g3_cb*wc(:nw-1)))

      op_cb2w(:nw-1) = op_cb1w(:nw-1)
      om_cb2w(:nw-1) = om_cb1w(:nw-1)
      g_cb2w(:nw-1)  = g_cb1w(:nw-1)

      op_oc1w(:nw-1) = d1_oc + wc(:nw-1)*(d2_oc + d3_oc*wc(:nw-1))
      g_oc1w(:nw-1)  = 0.50_dp*(g1_oc + wc(:nw-1)*(g2_oc + g3_oc*wc(:nw-1)))

      op_oc2w(:nw-1) = op_oc1w(:nw-1)

      g_oc2w(:nw-1)  = g_oc1w(:nw-1)

      op_antw(:nw-1) = d1_ant + wc(:nw-1)*(d2_ant + d3_ant*wc(:nw-1))
      g_antw(:nw-1)  = 0.63_dp*(g1_ant + wc(:nw-1)*(g2_ant + g3_ant*wc(:nw-1)))

      op_so4w(:nw-1) = d1_so4 + wc(:nw-1)*(d2_so4 + d3_so4*wc(:nw-1))

      op_salw(:nw-1) = d1_sal + wc(:nw-1)*(d2_sal + d3_sal*wc(:nw-1))

      end subroutine aer_init

      subroutine cnum( mix, denx, cden, molw, rm, &
                       rw, num, weight, nzlev )





      implicit none



      integer,  intent(in)  :: nzlev
      real(dp), intent(in)  :: molw, rm
      real(dp), intent(in)  :: cden
      real(dp), intent(in)  :: rw(nzlev-1)
      real(dp), intent(in)  :: mix(nzlev-1), denx(nzlev-1)
      real(dp), intent(out) :: num(nzlev-1), weight(nzlev-1)



      REAL(dp), PARAMETER :: PI = 3.14159265358979324
      REAL(dp), PARAMETER :: vol_fac = 4._dp * PI / 3._dp
      integer  :: k
      real(dp) :: rm_cubic, factor
      real(dp) :: mas(nzlev-1), masw(nzlev-1)

      rm_cubic = rm**3
      factor   = 1._dp/(vol_fac * cden * 1.e-12_dp * rm_cubic)



      where( mix(:nzlev-1) >= 1.e-15_dp )






         mas(:nzlev-1)    = mix(:nzlev-1) * 1.e-12_dp                              
         num(:nzlev-1)    = mas(:nzlev-1) * factor
      elsewhere
        num(:nzlev-1)    = 0._dp
        weight(:nzlev-1) = 0._dp
      endwhere

      do k = 1,nzlev-1
         if( abs(rw(k) - rm) >= 1.e-6_dp ) then
           masw(k) = vol_fac*(rw(k)**3 - rm_cubic)*1.e-12_dp * num(k)   
         else
           masw(k) = 0._dp
         end if
      end do

      where( mix(:nzlev-1) >= 1.e-15_dp )
         weight(:nzlev-1) = mas(:nzlev-1)/( mas(:nzlev-1) + masw(:nzlev-1) )        
         weight(:nzlev-1) = max( 0._dp, min( 1._dp,weight(:nzlev-1) ) )
      endwhere

      end subroutine cnum

      subroutine setair( nz, z, airlev, dtrl, &
                         cz, o2top )

      use module_wave_data, only : nw, wc






















      implicit none



      integer, intent(in) :: nz
      real(dp), intent(in) :: o2top
      real(dp), intent(in) :: z(nz)
      real(dp), intent(in) :: airlev(nz)




      real(dp), intent(out) :: dtrl(nz-1,nw-1)
      real(dp), intent(out) :: cz(nz)



      integer :: i, ip1, iw
      real(dp) :: hscale
      real(dp) :: srayl
      real(dp) :: deltaz
      real(dp) :: wmicrn, xx



      do i = 1,nz-1
         ip1 = i + 1
         deltaz = km2cm * (z(ip1) - z(i))
         cz(i) = (airlev(ip1) - airlev(i)) / log(airlev(ip1)/airlev(i)) * deltaz
      end do






      cz(nz) = o2top/.2095_dp



      do iw = 1,nw-1






         wmicrn = 1.e-3_dp*wc(iw)
         if( wmicrn <= .55_dp ) then
            xx = 3.6772_dp + 0.389_dp*wmicrn + 0.09426_dp/wmicrn
         else
            xx = 4.04_dp
         end if
         srayl = 4.02e-28_dp/(wmicrn)**xx
         dtrl(1:nz-2,iw) = cz(1:nz-2)*srayl
         dtrl(nz-1,iw) = (cz(nz-1) + cz(nz))*srayl
      end do

      end subroutine setair

      subroutine zenith( lat, long, idate, ut, zen )


















      implicit none



      integer, intent(in) :: idate
      real(dp), intent(in) :: lat,long
      real(dp), intent(in) :: ut
      real(dp), intent(out) :: zen



      real(dp), parameter :: pi = 3.1415926
      real(dp), parameter :: d2r = pi/180.0
      real(dp), parameter :: r2d = 1.0/d2r
      integer :: i
      integer :: iiyear, imth, iday, ijd
      integer :: imn(12) = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
      real(dp) :: lbut,lzut
      real(dp) :: rlt
      real(dp) :: d, tz, rdecl, eqr, eqh, zpt
      real(dp) :: csz, zr, caz, raz
      real(dp) :: sintz, costz, sin2tz, cos2tz, sin3tz, cos3tz



      rlt = lat*d2r



      iiyear = idate/10000
      imth = (idate - iiyear*10000)/100
      iday = idate - iiyear*10000 - imth*100



      if( mod(iiyear,4) == 0 ) then
         imn(2) = 29
      else
         imn(2) = 28
      end if



      ijd = 0
      do i = 1,imth-1
         ijd = ijd + imn(i)
      end do
      ijd = ijd + iday



      d = real(ijd-1,kind=dp) + ut/24._dp



      tz = 2._dp*pi*d/365._dp













      sintz = sin( tz )
      costz = cos( tz )
      sin2tz = 2._dp*sintz*costz
      cos2tz = costz*costz - sintz*sintz
      sin3tz = sintz*cos2tz + costz*sin2tz
      cos3tz = costz*cos2tz - sintz*sin2tz



      rdecl = 0.006918_dp - 0.399912_dp*costz + 0.070257_dp*sintz &
                          - 0.006758_dp*cos2tz + 0.000907_dp*sin2tz &
                          - 0.002697_dp*cos3tz + 0.001480_dp*sin3tz



      eqr = 0.000075_dp + 0.001868_dp*costz - 0.032077_dp*sintz &
                        - 0.014615_dp*cos2tz - 0.040849_dp*sin2tz



      eqh = eqr*24._dp/(2._dp*pi)



      lbut = 12._dp - eqh - long*24._dp/360._dp



      lzut = 15._dp*(ut - lbut)
      zpt = lzut*d2r



      csz = sin(rlt)*sin(rdecl) + cos(rlt)*cos(rdecl)*cos(zpt)
      zr = acos(csz)
      zen = zr*r2d

      end subroutine zenith

      subroutine calc_zenith( lat, long, ijd, gmt, azimuth, zenith )













        integer, intent(in)   :: ijd
        real(dp), intent(in)  :: gmt, lat, long
        real(dp), intent(out) :: azimuth, zenith




      real(dp), parameter :: d2r = 3.1415926_dp/180.0_dp
      real(dp), parameter :: r2d = 1.0_dp/d2r

      integer  :: jd
      real(dp) :: caz, csz, cw, d, decl, ec, epsi, eqt, eyt, feqt, feqt1,       &
          feqt2, feqt3, feqt4, feqt5, feqt6, feqt7, lbgmt, lzgmt, ml, pepsi,     &
          pi, ra, raz, rdecl, reqt, rlt, rml, rphi, rra, ssw, sw, tab, w, wr,    &
          yt, zpt, zr



        intrinsic acos, atan, cos, min, sin, tan




        rlt = lat*d2r
        rphi = long*d2r

        jd = ijd

        d = jd + gmt/24.0_dp



        ml = 279.2801988_dp + d*(.9856473354_dp + 2.267E-13_dp*d)
        rml = ml*d2r







        w = 282.4932328_dp + d*(4.70684E-5_dp + 3.39E-13_dp*d)
        wr = w*d2r
        ec   = 1.6720041E-2_dp - d*(1.1444E-9_dp + 9.4E-17_dp*d)
        epsi = 23.44266511_dp - d*(3.5626E-7_dp + 1.23E-15_dp*d)
        pepsi = epsi*d2r
        yt = (tan(pepsi/2.0_dp))**2
        cw = cos(wr)
        sw = sin(wr)
        ssw = sin(2.0_dp*wr)
        eyt = 2._dp*ec*yt
        feqt1 = -sin(rml)*cw*(eyt + 2._dp*ec)
        feqt2 = cos(rml)*sw*(2._dp*ec - eyt)
        feqt3 = sin(2._dp*rml)*(yt - (5._dp*ec**2/4._dp)*(cw**2 - sw**2))
        feqt4 = cos(2._dp*rml)*(5._dp*ec**2*ssw/4._dp)
        feqt5 = sin(3._dp*rml)*(eyt*cw)
        feqt6 = -cos(3._dp*rml)*(eyt*sw)
        feqt7 = -sin(4._dp*rml)*(.5_dp*yt**2)
        feqt = feqt1 + feqt2 + feqt3 + feqt4 + feqt5 + feqt6 + feqt7
        eqt = feqt*13751.0_dp




        reqt = eqt/240._dp



        ra = ml - reqt
        rra = ra*d2r



        tab = 0.43360_dp*sin(rra)
        rdecl = atan(tab)
        decl = rdecl/d2r



        lbgmt = 12.0_dp - eqt/3600._dp + long*24._dp/360._dp
        lzgmt = 15.0_dp*(gmt-lbgmt)
        zpt = lzgmt*d2r
        csz = sin(rlt)*sin(rdecl) + cos(rlt)*cos(rdecl)*cos(zpt)
        if(csz.gt.1)print *,'calczen,csz ',csz
        csz = min( 1._dp,csz )


        zr = acos(csz)
        zenith = zr/d2r



        caz = (sin(rdecl)-sin(rlt)*cos(zr))/(cos(rlt)*sin(zr))
        if(caz < -0.999999_dp)then
          azimuth=180._dp
        elseif(caz > 0.999999_dp)then
          azimuth=0._dp
        else
          raz = acos(caz)
          azimuth = raz/d2r
        endif






        if( lzgmt > 0._dp ) azimuth = azimuth + (2._dp*(180._dp - azimuth))

      end subroutine calc_zenith

      subroutine sundis( julday, esrm2 )











      implicit none



      integer, intent(in) :: julday
      real(dp), intent(out) :: esrm2



      real(dp), parameter :: pi = 3.1415926_dp
      real(dp) :: dayn, thet0
      real(dp) :: sinth, costh, sin2th, cos2th



      dayn = real(julday - 1,kind=dp) + .5_dp



      thet0 = 2._dp*pi*dayn/365._dp







      sinth  = sin( thet0 )
      costh  = cos( thet0 )
      sin2th = 2._dp*sinth*costh
      cos2th = costh*costh - sinth*sinth
      esrm2  = 1.000110_dp + .034221_dp*costh + .001280_dp*sinth + .000719_dp*cos2th + .000077_dp*sin2th

      end subroutine sundis

      subroutine spheres( nz, z, zen, dsdh, nid )



















      implicit none



      integer, intent(in)  :: nz
      integer, intent(out) :: nid(0:nz-1)
      real(dp), intent(in) :: zen
      real(dp), intent(in) :: z(nz)
      real(dp), intent(out) :: dsdh(0:nz-1,nz-1)



      real(dp), parameter :: radius = 6.371e+3_dp         
      real(dp), parameter :: d2r = 3.1415926_dp/180.0_dp
      integer :: i, j, k
      integer :: id
      real(dp) :: re, ze(nz)
      real(dp) :: zd(0:nz-1)
      real(dp) :: zenrad, rpsinz, rj, rjp1, dsj, dhj, ga, gb, sm

      zenrad = zen*d2r



      re = radius + z(1)



      ze(1:nz) = z(1:nz) - z(1)



      zd(0) = ze(nz)
      do k = 1,nz-1
        zd(k) = ze(nz - k)
      end do



      dsdh(0:nz-1,1:nz-1) = 0._dp
      nid(0:nz-1) = 0



      do i = 0,nz-1
        rpsinz = (re + zd(i)) * sin(zenrad)
        if( zen > 90._dp .and. rpsinz < re ) then
           nid(i) = -1
        else



           id = i
           if( zen > 90._dp ) then
              do j = 1,nz-1
                 if( (rpsinz < ( zd(j-1) + re ) ) .and. (rpsinz >= ( zd(j) + re )) ) then
                    id = j
                 end if
              end do
           end if
           do j = 1,id
             if( j == id .and. id == i .and. zen > 90._dp ) then
                sm = -1._dp
             else
                sm = 1._dp
             end if
             rj = re + zd(j-1)
             rjp1 = re + zd(j)
             dhj = zd(j-1) - zd(j)
             ga = max( 0._dp,rj*rj - rpsinz*rpsinz )
             gb = max( 0._dp,rjp1*rjp1 - rpsinz*rpsinz )
             if( id > i .and. j == id ) then
                dsj = sqrt( ga )
             else
                dsj = sqrt( ga ) - sm*sqrt( gb )
             end if
             dsdh(i,j) = dsj / dhj
           end do
           nid(i) = id
        end if
      end do

      end subroutine spheres

      subroutine setz( nz, cz, tlev, c, ndx, adjcoe )

      use module_wave_data, only : tuv_jmax



      implicit none



      integer, intent(in) :: nz
      integer, intent(in) :: ndx
      real(dp), intent(in) :: cz(nz)
      real(dp), intent(in) :: tlev(nz)
      real(dp), intent(in) :: c(5,tuv_jmax)
      real(dp), intent(inout) :: adjcoe(nz,tuv_jmax)



      integer :: k, m
      real(dp) :: tt, adjin
      real(dp) :: c0, c1, c2
      real(dp) :: xz(nz)
































      xz(1:nz) = cz(1:nz)*1.e-18_dp
      do k = 1,tuv_jmax
         adjcoe(1:nz,k) = 1._dp
      end do
      tt = tlev(1)/281._dp
      do m = 1,tuv_jmax
         adjin = 1._dp
         if( m == 2 ) then




            select case( ndx )
               case( 20 )
                  c0 = 4.52372_dp ; c1 = -5.94317_dp ; c2 = 2.63156_dp
               case( 40 )
                  c0 = 4.99378_dp ; c1 = -7.92752_dp ; c2 = 3.94715_dp
               case( 60 )
                  c0 = .969867_dp ; c1 = -.841035_dp ; c2 = .878835_dp
               case( 80 )
                  c0 = 1.07801_dp ; c1 = -2.39580_dp ; c2 = 2.32632_dp
            end select
            adjin = c0 + tt*(c1 + c2*tt)
         else if( m == 11 ) then




            select case( ndx )
               case( 20 )
                  c0 = 2.43360_dp ; c1 = -3.61363_dp ; c2 = 2.19018_dp
               case( 40 )
                  c0 = 3.98265_dp ; c1 = -6.90516_dp ; c2 = 3.93602_dp
               case( 60 )
                  c0 = 3.49843_dp ; c1 = -5.98839_dp ; c2 = 3.50262_dp
               case( 80 )
                  c0 = 3.06312_dp ; c1 = -5.26281_dp ; c2 = 3.20980_dp
            end select
            adjin = c0 + tt*(c1 + c2*tt)
         end if
         call calcoe( nz, c(1,m), xz, tt, adjin, adjcoe(1,m) )
      end do

      end subroutine setz

      subroutine setozo( nz, z, tlay, dto3, &
                         to3, o3, airlev, o3top, o3toms )

      use module_wave_data, only : nw, wl, xso3, s226, s263, s298































      implicit none




      integer, intent(in) :: nz
      real(dp), intent(in) :: o3top
      real(dp), intent(in) :: o3toms
      real(dp), intent(in) :: z(nz)
      real(dp), intent(in) :: tlay(nz-1)
      real(dp), intent(in) :: airlev(nz)
      real(dp), intent(out) :: dto3(nz-1,nw-1)
      real(dp), intent(out) :: to3(nz)
      real(dp), intent(in) :: o3(nz)



      integer  :: k, iw, nd, ilat
      real(dp) :: so3
      real(dp) :: scale
      real(dp) :: div1, div2
      real(dp) :: o3den(nz)
      real(dp) :: cz(nz)







      o3den(1:nz) = o3(1:nz)*airlev(1:nz)
      cz(1:nz-1)  = 0.5_dp*(o3den(2:nz) + o3den(1:nz-1))*km2cm*(z(2:nz) - z(1:nz-1))

      to3(nz)  = o3top
      do k = nz-1,1,-1
         to3(k) = to3(k+1) + cz(k)
      end do


      cz(nz-1)    = cz(nz-1) + o3top



      if( o3toms > 0.0_dp ) then
        scale = o3toms/(to3(1)/2.687e16_dp)
        cz(1:nz)  = cz(1:nz)*scale
        to3(1:nz) = to3(1:nz)*scale
      endif




      div1 = 1._dp / ( 263._dp - 226._dp)
      div2 = 1._dp / ( 298._dp - 263._dp)
      do iw = 1,nw-1
         so3 = xso3(iw)
         do k = 1,nz-1
            if( wl(iw) > 240.5_dp .and. wl(iw+1) < 350._dp ) then
               if( tlay(k) < 263._dp ) then
                  so3 = s226(iw) + (s263(iw) - s226(iw)) * div1 * (tlay(k) - 226._dp)
               else
                  so3 = s263(iw) + (s298(iw) - s263(iw)) * div2 * (tlay(k) - 263._dp)
               end if
            end if
            dto3(k,iw) = cz(k)*so3
         end do
      end do

      end subroutine setozo

      subroutine set_o2_xsect( nz, z, cz, vcol, scol, dto2, xso2 )

      use module_wave_data, only : nw, wl

































      implicit none



      integer, intent(in) :: nz
      real(dp), intent(in) :: cz(nz)
      real(dp), intent(in) :: z(nz)
      real(dp), intent(in) :: vcol(nz)
      real(dp), intent(in) :: scol(nz)
      real(dp), intent(out) :: dto2(nz-1,nw-1)
      real(dp), intent(out) :: xso2(nw-1,nz)



      integer :: iw, iz, igast
      real(dp) :: secchi(nz)



      real(dp) :: dto2k(nz-1,ngast-1)
      real(dp) :: xso2k(nz,ngast-1)



      real(dp) :: dto2la(nz-1,nla-1)
      real(dp) :: xso2la(nz,nla-1)





      real(dp), dimension(2*nwint) :: dttmp, xstmp    
      real(dp), dimension(nwint) :: dtuser, xsuser    
      real(dp) :: o2col(nz)
      real(dp) :: x, y
      real(dp) :: delo2




      dto2(:nz-1,:nw-1) = 0._dp
      xso2(:nw-1,:nz) = 0._dp
      if( wl(1) > 243._dp ) then
         return
      end if











      o2col(1:nz) = 0.2095_dp * scol(1:nz)




      secchi(1:nz-1) = scol(1:nz-1)/vcol(1:nz-1)
      where( secchi(1:nz-1) == 0._dp )
         secchi(1:nz-1) = 2._dp
      endwhere
      secchi(nz) = secchi(nz-1)





      if( wl(1) < wlgast(ngast) .and. wl(nw) > wlgast(1) ) then
           call schu( nz, o2col, secchi, dto2k, xso2k )
      else
         dto2k(:,:) = 0._dp
         xso2k(:,:) = 0._dp
      end if




      if( wl(1) <= wlla(nla) .and. wl(nw) >= wlla(1) ) then
         call lymana( nz, o2col, secchi, dto2la, xso2la )
      else
         dto2la(:,:) = 0._dp
         xso2la(:,:) = 0._dp
      end if



      do iz = 1,nz
         igast = 0



         do iw = 1,nwint-1






            if( wlint(iw+1) <= wlgast(1) .or. wlint(iw) >= wlgast(ngast) ) then
              if( wlint(iw+1) <= wlla(1) .or. wlint(iw) >= wlla(nla) ) then
                 xstmp(iw) = xso2int(iw)
              else
                 xstmp(iw) = xso2la(iz,1)
              end if
            else
               igast = igast + 1
               xstmp(iw) = xso2k(iz,igast)
            end if



            xstmp(iw) = xstmp(iw) * (wlint(iw+1) - wlint(iw))
         end do



         call inter3( nw, wl, xsuser, nwint, wlint, xstmp )
         xso2(1:nw-1,iz) = xsuser(1:nw-1)/(wl(2:nw) - wl(1:nw-1))
      end do
      do iz = 1,nz-1
         igast = 0
         delo2 = .2095_dp * cz(iz) 



         do iw = 1,nwint-1






            if( wlint(iw+1) <= wlgast(1) .or. wlint(iw) >= wlgast(ngast) ) then
              if( wlint(iw+1) <= wlla(1) .or. wlint(iw) >= wlla(nla) ) then
                 dttmp(iw) = xso2int(iw) * delo2
              else
                 dttmp(iw) = dto2la(iz,1)
              end if
            else
               igast = igast + 1
               dttmp(iw) = dto2k(iz,igast)
            end if



            dttmp(iw) = dttmp(iw) * (wlint(iw+1) - wlint(iw))
         end do



         call inter3( nw, wl, dtuser, nwint, wlint, dttmp )
         dto2(iz,1:nw-1) = dtuser(1:nw-1)/(wl(2:nw) - wl(1:nw-1))
      end do

      end subroutine set_o2_xsect

      subroutine setcld( nz, z, xlwc, dtcld, omcld, gcld )

      use module_wave_data, only : nw




















      implicit none



      integer, intent(in) :: nz
      real(dp), intent(in) :: z(nz)
      real(dp), intent(in) :: xlwc(nz)
      real(dp), intent(out) :: dtcld(nz-1,nw-1)
      real(dp), intent(out) :: omcld(nz-1,nw-1)
      real(dp), intent(out) :: gcld(nz-1,nw-1)



      real(dp), parameter :: wden = 1.0_dp * 1.e6_dp 
      real(dp), parameter :: re = 10.0_dp * 1.e-6_dp 

      integer  :: i                                  
      real(dp) :: dz(nz-1)



      do i=1,nz-1
          dz(i) = (z(i+1)-z(i)) * 1.e3_dp 
      end do




      do i = 1,nz-1
        dtcld(i,:) = 1.5_dp * xlwc(i)*dz(i)/ (wden * re)
        dtcld(i,:) = max( dtcld(i,:), 0._dp )
        omcld(i,:) = .9999_dp
        gcld (i,:) = .85_dp
      end do

      end subroutine setcld

      subroutine schu_inti



      implicit none



      integer :: i, iw, k, j1, jp1, j
      real(dp) :: col
      real(dp), dimension(6) :: a0, a1, b0, b1
      do iw = 1,ngast-1
         x_table(0,iw) = sum( a_schu(1:11:2,iw) )
         d_table(0,iw) = sum( b_schu(1:11:2,iw) )
         do k = 1,tdim
            col = 22._dp + t_del*real(k-1,kind=dp)
            o2_table(k) = col
            col = 10._dp**col
            a1(:) = a_schu(2:12:2,iw)*col
            b1(:) = b_schu(2:12:2,iw)*col
            where( a1(:) < 500._dp )
               a0(:) = exp( -a1(:) )
            elsewhere
               a0(:) = 0._dp
            endwhere
            where( b1(:) < 500._dp )
               b0(:) = exp( -b1(:) )
            elsewhere
               b0(:) = 0._dp
            endwhere
            x_table(k,iw) = dot_product( a_schu(1:11:2,iw),a0(:) )
            d_table(k,iw) = dot_product( b_schu(1:11:2,iw),b0(:) )
         end do
      end do

      end subroutine schu_inti

      subroutine schu( nz, o2col, secchi, dto2, xscho2 )



















      implicit none



      integer, intent(in) :: nz
      real(dp), intent(inout) :: dto2(nz-1,ngast-1)
      real(dp), intent(inout) :: xscho2(nz,ngast-1)
      real(dp), intent(in) :: o2col(nz)
      real(dp), intent(in) :: secchi(nz)



      integer :: i, iw, j, j1, jp1, k, ki, kp1
      integer :: index(nz)
      integer :: minind(1), maxind(1)
      real(dp) :: a0, a1, b0, b1
      real(dp), dimension(6) :: ac, bc
      real(dp), dimension(nz,6) :: aa, bb
      real(dp), dimension(nz) :: rjm, rjo2
      real(dp), dimension(nz) :: rjmi, rjo2i, lo2col, dels



      rjm(1:nz) = 0._dp
      rjo2(1:nz) = 0._dp



      where( o2col(:) /= 0 )
         lo2col(:) = log10( o2col(:) )
      endwhere
      do ki = 1,nz
         if( o2col(ki) /= 0._dp ) then
            if( lo2col(ki) <= o2_table(1) ) then
               dels(ki) = 0._dp
               index(ki) = 1
            else if( lo2col(ki) >= o2_table(tdim) ) then
               dels(ki) = 1._dp
               index(ki) = tdim-1
            else
               do k = 2,tdim
                  if( lo2col(ki) <= o2_table(k) ) then
                     dels(ki) = t_fac*(lo2col(ki) - o2_table(k-1))
                     index(ki) = k-1
                     exit
                  end if
               end do
            end if
         else
            index(ki) = 0
            dels(ki) = 0._dp
         end if
      end do



      do iw = 1,ngast-1
         do k = 1,nz
            ki = index(k)
            rjm(k) = x_table(ki,iw) + dels(k)*(x_table(ki+1,iw) - x_table(ki,iw))
            rjo2(k) = d_table(ki,iw) + dels(k)*(d_table(ki+1,iw) - d_table(ki,iw))
         end do
         do k = 1,nz-1
            if( rjm(k) > 1.e-100_dp ) then
               kp1 = k + 1
               if( rjm(kp1) > 0._dp ) then
                  dto2(k,iw) = log( rjm(kp1) ) / secchi(kp1) - log( rjm(k) ) * secchi(k)
               else
                  dto2(k,iw) = 1000._dp
               end if
            else
               dto2(k,iw) = 1000._dp
            end if
         end do
         do k = 1,nz
            if( rjm(k) > 1.e-100_dp ) then
               if( rjo2(k) > 1.e-100_dp ) then
                  xscho2(k,iw) = rjo2(k)/rjm(k)
               else
                  xscho2(k,iw) = 0._dp
               end if
            else
               xscho2(k,iw) = 0._dp
            end if
         end do
      end do

      end subroutine schu

      subroutine rtlink( nz,                  &
                         nw,                  &
                         zen,                 &
                         albedo,              & 
                         z,                   &
                         nid,                 &
                         dsdh,                & 
                         dtrl,                &
                         dto3,                &
                         dto2,                &
                         dtcld, omcld, gcld,  &








                         dtaer, omaer, gaer,  &  
                         radfld )

      implicit none



      integer, intent(in) :: nz
      integer, intent(in) :: nw
      real(dp), intent(in)    :: z(nz)
      real(dp), intent(in)    :: zen
      real(dp), intent(in)    :: albedo(nw-1)
      integer, intent(in)    :: nid(0:nz-1)
      real(dp), intent(in)    :: dsdh(0:nz-1,nz-1)
      real(dp), intent(in)    :: dtrl(nz-1,nw-1)
      real(dp), intent(in)    :: dto3(nz-1,nw-1)
      real(dp), intent(in)    :: dto2(nz-1,nw-1)
      real(dp), intent(in)    :: dtcld(nz-1,nw-1), omcld(nz-1,nw-1), gcld(nz-1,nw-1)








      real(dp), intent(in)    :: dtaer(nz-1,nw-1), omaer(nz-1,nw-1), gaer(nz-1,nw-1) 
      real(dp), intent(out)   :: radfld(nz,nw-1)           



      real(dp), parameter :: smallest = tiny( 1.0 )
      integer :: i, ii, iw
      real(dp) :: dtsct, dtabs
      real(dp) :: dscld, dacld








      real(dp) :: dsaer, daaer
      real(dp), dimension(nz-1,nw-1) :: dt, om, g



      do iw = 1,nw-1
         do i = 1,nz-1

            dscld = dtcld(i,iw)*omcld(i,iw)
            dacld = dtcld(i,iw)*(1._dp - omcld(i,iw))























            dsaer = dtaer(i,iw)*omaer(i,iw)
            daaer = dtaer(i,iw)*(1._dp - omaer(i,iw))
            dtsct = dtrl(i,iw) + dscld + dsaer 

            dtabs = dto3(i,iw) + dto2(i,iw) + dacld + daaer


            dtabs = max( dtabs,smallest )
            dtsct = max( dtsct,smallest )



            ii = nz - i

            dt(ii,iw) = dtsct + dtabs
            if( dtsct /= smallest ) then
               om(ii,iw) = dtsct/(dtsct + dtabs)
               g(ii,iw) = ( gcld(i,iw)*dscld + &







                            gaer(i,iw)*dsaer ) / dtsct
               g(ii,iw) = max( smallest, g(ii,iw) )
               g(ii,iw) = min( 1.d0, g(ii,iw) )
            else
               om(ii,iw) = smallest
               g(ii,iw)  = smallest
            end if

         end do
      end do




      call ps2str( nz, nw, zen, albedo, dt, om, &
                   g, dsdh, nid, radfld )

      end subroutine rtlink

      subroutine tridec( syscnt, order, lower, main, upper )





















































































      implicit none



      integer, intent(in) :: syscnt, order
      real(dp), intent(in), dimension(syscnt,order) :: lower
      real(dp), intent(inout), dimension(syscnt,order) :: main, upper



      integer :: i



      main(:,1) = 1._dp / main(:,1)
      upper(:,1) = upper(:,1)*main(:,1)
      do i = 2,order-1
         main(:,i) = 1._dp / (main(:,i) - lower(:,i)*upper(:,i-1))
         upper(:,i) = upper(:,i)*main(:,i)
      end do

      end subroutine tridec

      subroutine trislv( syscnt, order, lower, main, upper, x )

      implicit none



      integer, intent(in) :: syscnt, order
      real(dp), intent(in), dimension(syscnt,order) :: lower, &
                                                       main, &
                                                       upper
      real(dp), intent(inout), dimension(syscnt,order) :: x



      integer :: i, im1, j, n, nm1
      nm1 = order - 1
      n = order



      x(:,1) = x(:,1)*main(:,1)
      do i = 2,nm1
         x(:,i) = (x(:,i) - lower(:,i)*x(:,i-1))*main(:,i)
      end do
      x(:,n) = (x(:,n) - lower(:,n)*x(:,nm1))/(main(:,n) - lower(:,n)*upper(:,nm1))
      do i = nm1,1,-1
         x(:,i) = x(:,i) - upper(:,i)*x(:,i+1)
      end do

      end subroutine trislv

      subroutine ps2str( nz, nw, zen, rsfc, tauu, omu, &
                         gu, dsdh, nid, radfld )






































      implicit none



      integer, intent(in) :: nz, nw
      integer, intent(in) :: nid(0:nz-1)
      real(dp), intent(in) :: zen
      real(dp), intent(in) :: rsfc(nw-1)
      real(dp), dimension(nz-1,nw-1), intent(in) :: tauu, omu, gu
      real(dp), intent(in) :: dsdh(0:nz-1,nz-1)
      real(dp), intent(out) :: radfld(nz,nw-1)













      real(dp), parameter :: smallest = tiny( 1.0_dp )
      real(dp), parameter :: largest  = huge( 1.0_dp )
      real(dp), parameter :: d2r = 3.1415926_dp/180.0_dp
      integer :: mrows, nzm1, nzm2
      real(dp), parameter :: eps = 1.e-3_dp
      real(dp), parameter :: pifs = 1._dp
      real(dp), parameter :: fdn0 = 0._dp
      integer :: row
      integer :: lev
      integer :: i, ip1, iw
      integer :: j, jl, ju
      real(dp) :: precis, wrk
      real(dp) :: tempg
      real(dp) :: mu, suma
      real(dp) :: g, om
      real(dp) :: gam1, gam2, gam3, gam4
      real(dp), dimension(nz-1) :: f, gi, omi
      real(dp), dimension(0:nz) :: tauc, mu2
      real(dp), dimension(nz-1) :: lam, taun, bgam
      real(dp), dimension(nz-1) :: cdn
      real(dp), dimension(0:nz,nw-1) :: tausla
      real(dp), dimension(nz-1,nw-1) :: cup, cuptn, cdntn
      real(dp), dimension(nz-1,nw-1) :: e1, e2, e3, e4
      real(dp), dimension(2*(nz-1)) :: a, b, d, e
      real(dp), dimension(nw-1,2*(nz-1)) :: sub, main, super, y









      real(dp) :: expon, expon0, expon1, divisr, temp, up, dn
      real(dp) :: ssfc



      nzm1 = nz - 1
      nzm2 = nz - 2
      mrows = 2*nzm1
      precis = epsilon( precis )
      mu = cos( zen*d2r )
wave_loop : &
      do iw = 1,nw-1









         tauc(0:nz) = 0._dp
         tausla(0:nz,iw) = 0._dp
         mu2(0:nz) = sqrt( smallest )





         f(1:nzm1) = gu(:,iw)*gu(:,iw)
         gi(1:nzm1) = (gu(:,iw) - f(1:nzm1))/(1._dp - f(1:nzm1))
         omi(1:nzm1) = (1._dp - f(1:nzm1))*omu(1:nzm1,iw)/(1._dp - omu(1:nzm1,iw)*f(1:nzm1))
         taun(1:nzm1) = (1._dp - omu(1:nzm1,iw)*f(1:nzm1))*tauu(1:nzm1,iw)





         if( zen > 90._dp ) then
            if( nid(0) < 0 ) then
               tausla(0,iw) = largest
            else
               ju = nid(0)
               tausla(0,iw) = 2._dp*dot_product( taun(1:ju),dsdh(0,1:ju) )
            end if
         end if
level_loop : &
         do i = 1,nzm1
            g = gi(i)
            om = omi(i)
            tauc(i) = tauc(i-1) + taun(i)



            tempg = min( abs(g),1._dp - precis )
            g = sign( tempg,g )
            om = min( om,1._dp - precis )



            if( nid(i) < 0 ) then
               tausla(i,iw) = largest
            else
               ju = min( nid(i),i )
               suma = dot_product( taun(1:ju),dsdh(i,1:ju) )
               jl = min( nid(i),i ) + 1
               tausla(i,iw) = suma + 2._dp*dot_product( taun(jl:nid(i)),dsdh(i,jl:nid(i)) )
               if( abs(tausla(i,iw)-tausla(i-1,iw)) < 1.0e-30_dp ) then
                 mu2(i) = sqrt( largest )
               else
                 mu2(i) = (tauc(i) - tauc(i-1))/(tausla(i,iw) - tausla(i-1,iw))
                 mu2(i) = sign( max( abs(mu2(i)),sqrt(smallest) ),mu2(i) )
               end if
            end if




            gam1 = .25_dp*(7._dp - om*(4._dp + 3._dp*g))
            gam2 = -.25_dp*(1._dp - om*(4._dp - 3._dp*g))
            gam3 = .25_dp*(2._dp - 3._dp*g*mu)
            gam4 = 1._dp - gam3




            lam(i) = sqrt( gam1*gam1 - gam2*gam2 )
            bgam(i) = (gam1 - lam(i))/gam2
            wrk = lam(i)*taun(i)
            if( wrk < 500._dp ) then
               expon = exp( -wrk )
            else
               expon = 0._dp
            end if



            e1(i,iw) = 1._dp + bgam(i)*expon
            e2(i,iw) = 1._dp - bgam(i)*expon
            e3(i,iw) = bgam(i) + expon
            e4(i,iw) = bgam(i) - expon






            if( tausla(i-1,iw) < 500._dp ) then
               expon0 = exp( -tausla(i-1,iw) )
            else
               expon0 = 0._dp
            end if
            if( tausla(i,iw) < 500._dp ) then
               expon1 = exp( -tausla(i,iw) )
            else
               expon1 = 0._dp
            end if

            divisr = lam(i)*lam(i) - 1._dp/(mu2(i)*mu2(i))
            temp = max( eps,abs(divisr) )
            divisr = 1._dp/sign( temp,divisr )
            up = om*pifs*((gam1 - 1._dp/mu2(i))*gam3 + gam4*gam2)*divisr
            dn = om*pifs*((gam1 + 1._dp/mu2(i))*gam4 + gam2*gam3)*divisr




            cup(i,iw) = up*expon0
            cdn(i) = dn*expon0
            cuptn(i,iw) = up*expon1
            cdntn(i,iw) = dn*expon1

         end do level_loop




        if( tausla(nzm1,iw) < 500._dp ) then
           ssfc = rsfc(iw)*mu*exp( -tausla(nzm1,iw) )*pifs
        else
           ssfc = 0._dp
        end if




        a(1) = 0._dp
        b(1) = e1(1,iw)
        d(1) = -e2(1,iw)
        e(1) = fdn0 - cdn(1)



        a(3:mrows-1:2) = e2(1:nzm2,iw)*e3(1:nzm2,iw) - e4(1:nzm2,iw)*e1(1:nzm2,iw)
        b(3:mrows-1:2) = e1(1:nzm2,iw)*e1(2:nzm1,iw) - e3(1:nzm2,iw)*e3(2:nzm1,iw)
        d(3:mrows-1:2) = e3(1:nzm2,iw)*e4(2:nzm1,iw) - e1(1:nzm2,iw)*e2(2:nzm1,iw)
        e(3:mrows-1:2) = e3(1:nzm2,iw)*(cup(2:nzm1,iw) - cuptn(1:nzm2,iw)) + e1(1:nzm2,iw)*(cdntn(1:nzm2,iw) - cdn(2:nzm1))



        a(2:mrows-2:2) = e2(2:nzm1,iw)*e1(1:nzm2,iw) - e3(1:nzm2,iw)*e4(2:nzm1,iw)
        b(2:mrows-2:2) = e2(1:nzm2,iw)*e2(2:nzm1,iw) - e4(1:nzm2,iw)*e4(2:nzm1,iw)
        d(2:mrows-2:2) = e1(2:nzm1,iw)*e4(2:nzm1,iw) - e2(2:nzm1,iw)*e3(2:nzm1,iw)
        e(2:mrows-2:2) = (cup(2:nzm1,iw) - cuptn(1:nzm2,iw))*e2(2:nzm1,iw) - (cdn(2:nzm1) - cdntn(1:nzm2,iw))*e4(2:nzm1,iw)



        a(mrows) = e1(nzm1,iw) - rsfc(iw)*e3(nzm1,iw)
        b(mrows) = e2(nzm1,iw) - rsfc(iw)*e4(nzm1,iw)
        d(mrows) = 0._dp
        e(mrows) = ssfc - cuptn(nzm1,iw) + rsfc(iw)*cdntn(nzm1,iw)
        sub(iw,1:mrows) = a(1:mrows)
        main(iw,1:mrows) = b(1:mrows)
        super(iw,1:mrows) = d(1:mrows)
        y(iw,1:mrows) = e(1:mrows)
      end do wave_loop



      call tridec( nw-1, mrows, sub, main, super )
      call trislv( nw-1, mrows, sub, main, super, y )



      do iw = 1,nw-1



         e(:mrows) = y(iw,:mrows)
         if( tausla(0,iw) < 500._dp ) then
            radfld(1,iw) = 2._dp*(fdn0 + e(1)*e3(1,iw) - e(2)*e4(1,iw) + cup(1,iw)) + exp( -tausla(0,iw) )
         else
            radfld(1,iw) = 2._dp*(fdn0 + e(1)*e3(1,iw) - e(2)*e4(1,iw) + cup(1,iw))
         end if
         where( tausla(1:nzm1,iw) < 500._dp )
            cdn(1:nzm1) = exp( -tausla(1:nzm1,iw) )
         elsewhere
            cdn(1:nzm1) = 0._dp
         endwhere
         radfld(2:nz,iw) = 2._dp*(e(1:mrows-1:2)*(e3(1:nzm1,iw) + e1(1:nzm1,iw)) &
                            + e(2:mrows:2)*(e4(1:nzm1,iw) + e2(1:nzm1,iw)) &
                            + cdntn(1:nzm1,iw) + cuptn(1:nzm1,iw)) + cdn(1:nzm1)
      end do

      end subroutine ps2str


      subroutine pchem( chem_opt, nz, tlev, airlev )

      use module_state_description, only : MOZART_KPP, MOZCART_KPP, T1_MOZCART_KPP, &
                                           MOZART_MOSAIC_4BIN_KPP, MOZART_MOSAIC_4BIN_AQ_KPP
      use module_wave_data, only : r01, r04, r44, r08, r06, r10,  &
                                   r11, r14, r15, r17, r18,       &
                                   xs_mvk, xs_macr, xs_hyac, xs_glyald, nj



























      implicit none




      integer, intent(in)  :: chem_opt
      integer, intent(in)  :: nz
      real(dp), intent(in) :: tlev(nz)
      real(dp), intent(in) :: airlev(nz)




      integer :: j
      character(len=40) :: jlabel(nj)





      j = 1
      jlabel(j) = 'o2 + hv -> o + o '



      call r01( nz, tlev, airlev, j, jlabel )



      j = j + 1
      jlabel(j) = 'no2 + hv -> no + o(3p) '



      j = j + 1
      jlabel(j) = 'no3 + hv -> no2 + o(3p) '



      j = j + 1
      jlabel(j) = 'no3 + hv -> no + o2 '



      call r04( nz,tlev,airlev,j,jlabel )



      call r44( nz,tlev,airlev,j,jlabel )



      j = j + 1
      jlabel(j) = 'ho2 + hv -> oh + o '



      call r08( nz,tlev,airlev,j,jlabel )



      j = j + 1
      jlabel(j) = 'hno2 + hv -> oh + no '



      call r06( nz,tlev,airlev,j,jlabel )



      j = j + 1
      jlabel(j) = 'hno4 + hv -> ho2 + no2 '



      call r10( nz,tlev,airlev,j,jlabel )



      call r11( nz,tlev,airlev,j,jlabel )



      j = j + 1
      jlabel(j) = 'c2h5cho + hv -> c2h5 + hco '



      j = j + 1
      jlabel(j) = 'chocho + hv -> products '



      call r14( nz,tlev,airlev,j,jlabel )



      call r15( nz,tlev,airlev,j,jlabel )



      j = j + 1
      jlabel(j) = 'ch3ooh + hv -> ch3o + oh '



      call r17( nz,tlev,airlev,j,jlabel )



      call r18( nz,tlev,airlev,j,jlabel )






is_mozart : &
      if( chem_opt == MOZART_KPP .or. chem_opt == MOZCART_KPP .or. &
          chem_opt == T1_MOZCART_KPP .or. &
          chem_opt == MOZART_MOSAIC_4BIN_KPP .or. &
          chem_opt == MOZART_MOSAIC_4BIN_AQ_KPP) then
        call xs_mvk( nz, tlev, airlev )
        call xs_macr( nz, tlev, airlev )
        call xs_hyac( nz, tlev, airlev )
        call xs_glyald( nz, tlev, airlev )
      else is_mozart



      j = j + 1
      jlabel(j) = 'cloo + hv -> products '



      j = j + 1
      jlabel(j) = 'clono2 + hv -> cl + no3'
      j = j + 1
      jlabel(j) = 'clono2 + hv -> clo + no2'



      j = j + 1
      jlabel(j) = 'ch3cl + hv -> products '



      j = j + 1
      jlabel(j) = 'ccl2o + hv -> products'



      j = j + 1
      jlabel(j) = 'ccl4 + hv -> products '



      j = j + 1
      jlabel(j) = 'cclfo + hv -> products '



      j = j + 1
      jlabel(j) = 'ccf2o + hv -> products '



      j = j + 1
      jlabel(j) = 'cf2clcfcl2 (cfc-113) + hv -> products'



      j = j + 1
      jlabel(j) = 'cf2clcf2cl (cfc-114) + hv -> products '



      j = j + 1
      jlabel(j) = 'cf3cf2cl (cfc-115) + hv -> products '



      j = j + 1
      jlabel(j) = 'ccl3f (cfc-111) + hv -> products'



      j = j + 1
      jlabel(j) = 'ccl2f2 (cfc-112) + hv -> products'



      j = j + 1
      jlabel(j)='ch3ccl3 + hv -> products'



      j = j + 1
      jlabel(j)='cf3chcl2 (hcfc-123) + hv -> products'



      j = j + 1
      jlabel(j)='cf3chfcl (hcfc-124) + hv -> products'



      j = j + 1
      jlabel(j)='ch3cfcl2 (hcfc-141b) + hv -> products'



      j = j + 1
      jlabel(j)='ch3cf2cl (hcfc-142b) + hv -> products'



      j = j + 1
      jlabel(j)='cf3cf2chcl2 (hcfc-225ca) + hv -> products'



      j = j + 1
      jlabel(j)='cf2clcf2chfcl (hcfc-225cb) + hv -> products'



      j = j + 1
      jlabel(j)='chclf2 (hcfc-22) + hv -> products'



      j = j + 1
      jlabel(j)='brono2 + hv -> br + o + no2'



      j = j + 1
      jlabel(j)='brono2 + hv -> bro + no2'



      j = j + 1
      jlabel(j)=' ch3br + hv -> products'



      j = j + 1
      jlabel(j)='chbr3 + hv -> products'



      j = j + 1
      jlabel(j)='cf3br (halon-1301) + hv -> products'



      j = j + 1
      jlabel(j)='cf2brcf2br (halon-2402) + hv -> products'



      j = j + 1
      jlabel(j)='cf2br2 (halon-1202) + hv -> products'



      j = j + 1
      jlabel(j)='cf2brcl (halon-1211) + hv -> products '



      j = j + 1
      jlabel(j)='cl2 + hv -> cl + cl '

! (57) hocl + hv -> oh + cl

      j = j + 1
      jlabel(j)='hocl + hv -> oh + cl'



      j = j + 1
      jlabel(j)='fmcl + hv -> cl + co + ho2'

      end if is_mozart


      end subroutine pchem

      subroutine lymana( nz, o2col, secchi, dto2la, xso2la )

















      implicit none



      integer, intent(in) :: nz
      real(dp), intent(in) :: o2col(nz)
      real(dp), intent(in) :: secchi(nz)
      real(dp), intent(out) :: dto2la(nz-1,nla-1)
      real(dp), intent(out) :: xso2la(nz,nla-1)



      integer :: i, k, kp1
      real(dp), dimension(nz) :: rm, ro2
      real(dp), save :: b(3), c(3), d(3), e(3)
      data b / 6.8431e-01_dp, 2.29841e-01_dp, 8.65412e-02_dp/, &
           c /8.22114e-21_dp, 1.77556e-20_dp, 8.22112e-21_dp/, &
           d / 6.0073e-21_dp, 4.28569e-21_dp, 1.28059e-20_dp/, &
           e /8.21666e-21_dp, 1.63296e-20_dp, 4.85121e-17_dp/



      rm(:) = 0._dp
      ro2(:) = 0._dp
      do k = 1,nz
        do i = 1,3
          rm(k) = rm(k) + b(i) * exp( -c(i)*o2col(k) )
          ro2(k) = ro2(k) + d(i) * exp( -e(i)*o2col(k) )
        end do
      end do



      do k = 1,nz-1
         if( rm(k) > 1.e-100_dp ) then
            kp1 = k + 1
            if( rm(kp1) > 0._dp ) then
               dto2la(k,1) = log( rm(kp1) )/secchi(kp1) - log( rm(k) )/secchi(k)
            else
               dto2la(k,1) = 1000._dp
            end if
         else
            dto2la(k,1) = 1000._dp
         end if
      end do
      do k = 1,nz
         if( rm(k) > 1.e-100_dp ) then
            if( ro2(k) > 1.e-100_dp ) then
               xso2la(k,1) = ro2(k)/rm(k)
            else
               xso2la(k,1) = 0._dp
            end if
         else
            xso2la(k,1) = 0._dp
         end if
      end do

      end subroutine lymana


      subroutine inter2( ng, xg, yg, n, x, y, ierr )





























      implicit none



      integer, intent(in) :: ng, n
      integer, intent(out) :: ierr
      real(dp), intent(in) :: x(n)
      real(dp), intent(in) :: y(n)
      real(dp), intent(in) :: xg(ng)
      real(dp), intent(out) :: yg(ng)



      integer :: ngintv
      integer :: i, k, jstart
      real(dp) :: area, xgl, xgu
      real(dp) :: darea, slope
      real(dp) :: a1, a2, b1, b2
      ierr = 0



      do i = 2,n
         if( x(i) <= x(i-1) ) then
            ierr = 1
            write(*,*) 'inter2: x coord not monotonically increasing'
            return
         end if
      end do
      do i = 2,ng
        if( xg(i) <= xg(i-1) ) then
           ierr = 2
           write(*,*) 'inter2: xg coord not monotonically increasing'
           return
        end if
      end do



      if( x(1) > xg(1) .or. x(n) < xg(ng) ) then
          write(*,*) 'inter2: data does not span grid'
          write(*,*) '        use addpnt to expand data and re-run'

      end if





      jstart = 1
      ngintv = ng - 1
      do i = 1,ngintv



         area = 0._dp
         xgl = xg(i)
         xgu = xg(i+1)







         k = jstart
         if( k <= n-1 ) then



            do
               if( x(k+1) <= xgl ) then
                  jstart = k - 1
                  k = k+1
                  if( k <= n-1 ) then
                     cycle
                  else
                     exit
                  end if
               else
                  exit
               end if
            end do




            do
               if( k <= n-1 .and. x(k) < xgu ) then
                  jstart = k-1



                  a1 = max( x(k),xgl )
                  a2 = min( x(k+1),xgu )



                  if( x(k+1) == x(k) ) then
                     darea = 0._dp
                  else
                     slope = (y(k+1) - y(k))/(x(k+1) - x(k))
                     b1 = y(k) + slope*(a1 - x(k))
                     b2 = y(k) + slope*(a2 - x(k))
                     darea = .5*(a2 - a1)*(b2 + b1)
                  end if



                  area = area + darea
                  k = k+1
                  cycle
               else
                  exit
               end if
            end do
         end if



         yg(i) = area/(xgu - xgl)
      end do

      end subroutine inter2

      subroutine inter_inti( ng, xg, n, x )



      implicit none



      integer, intent(in) :: ng, n
      real(dp), intent(in) :: xg(ng)
      real(dp), intent(in) :: x(n)



      integer :: i, ii, iil, astat
      integer :: ndim(1)
      allocate( xi(ng), xcnt(ng-1), stat=astat )
      if( astat /= 0 ) then
         write(*,*) 'inter_inti: failed to allocate wrk arrays; error = ',astat
         stop
      else
         xi(:) = 0
         xcnt(:) = 0
      end if
      iil = 1
      do i = 1,ng
         do ii = iil,n-1
            if( xg(i) < x(ii) ) then
               xi(i) = ii - 1
               iil = ii
               exit
            end if
         end do
      end do
      nintervals = count( xi(:) /= 0 )
      if( nintervals == 0 ) then
         write(*,*) 'inter_inti: wavelength grids do not overlap'
         stop
      else
         nintervals = nintervals - 1
      end if
      xcnt(1:nintervals) = xi(2:nintervals+1) - xi(1:nintervals) + 1
      ndim(:) = maxval( xcnt(1:nintervals) )
      allocate( xfrac(ndim(1),nintervals),stat=astat )
      if( astat /= 0 ) then
         write(*,*) 'inter_inti: failed to allocate wrk array; error = ',astat
         stop
      else
         xfrac(:,:) = 1._dp
      end if
      do i = 1,nintervals
        iil = xi(i)
        xfrac(1,i) = (min( x(iil+1),xg(i+1) ) - xg(i))/(x(iil+1) - x(iil))
        if( xcnt(i) > 1 ) then
           iil = xi(i) + xcnt(i) - 1
           xfrac(xcnt(i),i) = (xg(i+1) - x(iil))/(x(iil+1) - x(iil))
        end if
      end do










      end subroutine inter_inti

      subroutine inter3( ng, xg, yg, n, x, y )

































      implicit none



      integer, intent(in) :: n, ng
      real(dp), intent(in) :: xg(ng)
      real(dp), intent(in) :: x(n)
      real(dp), intent(in) :: y(n)
      real(dp), intent(out) :: yg(ng)



      integer :: i, ii, iil



      yg(:) = 0.
      do i = 1,nintervals
         iil = xi(i)
         ii = xcnt(i)
         if( ii == 1 ) then
            yg(i) = xfrac(1,i)*y(iil)
         else
            yg(i) = dot_product( xfrac(1:ii,i),y(iil:iil+ii-1) )
         end if
      end do

      end subroutine inter3

      subroutine calcoe( nz, c, xz, tt, adjin, adjcoe )






      implicit none



      integer, intent(in) :: nz
      real(dp), intent(in) :: adjin
      real(dp), intent(in) :: tt
      real(dp), intent(in) :: c(5)
      real(dp), intent(in) :: xz(nz)
      real(dp), intent(inout) :: adjcoe(nz)



      integer :: k
      real(dp) :: x
      do k = 1,nz
         x = xz(k)
         adjcoe(k) = adjin * (1._dp + .01_dp*(c(1) + x*(c(2) + x*(c(3) + x*(c(4) + x*c(5))))))
      end do

      end subroutine calcoe


      subroutine airmas( nz, z, zen, dsdh, nid, cz, &
                         vcol, scol )




















      implicit none



      integer, intent(in) :: nz
      integer, intent(in) :: nid(0:nz-1)
      real(dp), intent(in) :: z(nz)
      real(dp), intent(in) :: zen
      real(dp), intent(in) :: dsdh(0:nz-1,nz-1)
      real(dp), intent(in) :: cz(nz)
      real(dp), intent(out) :: vcol(nz)
      real(dp), intent(out) :: scol(nz)



      real(dp), parameter :: largest = huge(1.0_dp)



      integer :: id, j
      real(dp) :: sum, ssum, vsum, ratio



      vsum = 0._dp
      ssum = 0._dp
      do id = 1,nz-1
         vsum = vsum + cz(nz-id)
         vcol(nz-id) = vsum
         sum = 0.
         if( nid(id) < 0 ) then
            sum = largest
         else



            do j = 1,min( nid(id),id )
               sum = sum + cz(nz-j)*dsdh(id,j)
            end do



            do j = min( nid(id),id )+1,nid(id)
               sum = sum + 2._dp*cz(nz-j)*dsdh(id,j)
            end do
         end if
         scol(nz - id) = sum
      end do



      if( scol(nz-2) /= 0._dp ) then
         ratio = scol(nz-1)/scol(nz-2)
         scol(nz) = ratio * scol(nz-1)
      else
         scol(nz) = 0._dp
      end if

      end subroutine airmas

      END MODULE module_ftuv_subs 
