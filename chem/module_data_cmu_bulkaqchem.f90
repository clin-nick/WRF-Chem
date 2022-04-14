








      module module_data_cmu_bulkaqchem


      implicit none
































      double precision, parameter :: pi  = 3.14159
      double precision, parameter :: pi6 = pi/6.0

      double precision, parameter :: rho = 1.4e12       



      integer, parameter :: nas =  1           
      integer, parameter :: nah =  2           
      integer, parameter :: naa =  3           
      integer, parameter :: nan =  4           
      integer, parameter :: nac =  5           
      integer, parameter :: na4 =  6           
      integer, parameter :: naw =  7           
      integer, parameter :: nae =  8           
      integer, parameter :: nao =  9           
      integer, parameter :: nar = 10           
      integer, parameter :: nahso5 = 11        
      integer, parameter :: nahmsa = 12        
      integer, parameter :: naspec = 12        



      integer, parameter :: ngca =  1          
      integer, parameter :: ngcn =  2          
      integer, parameter :: ngcc =  3          
      integer, parameter :: ngc4 =  4          
      integer, parameter :: ngco =  5          
      integer, parameter :: ngcspec = 5        





      integer, parameter :: nga =  1           
      integer, parameter :: ngn =  2           
      integer, parameter :: ngc =  3           
      integer, parameter :: ng4 =  4           
      integer, parameter :: ngo =  5           
      integer, parameter :: ngspec = 5         




      integer, parameter :: ngtotal = 26 		
      integer, parameter :: ngas=ngtotal
      integer, parameter :: naers=naspec





















      integer, parameter :: ksod = nas               
      integer, parameter :: khyd = nah               
      integer, parameter :: kamm = naa               
      integer, parameter :: knit = nan               
      integer, parameter :: kchl = nac               
      integer, parameter :: ksvi = na4               
      integer, parameter :: kwat = naw               
      integer, parameter :: kec  = nae               
      integer, parameter :: koc  = nao               
      integer, parameter :: kcru = nar               






      integer, parameter :: ngso2     = 11
      integer, parameter :: ngh2o2    = 12
      integer, parameter :: nghcho    = 13
      integer, parameter :: nghcooh   = 14
      integer, parameter :: nghno2    = 15
      integer, parameter :: ngno      = 16
      integer, parameter :: ngno2     = 17
      integer, parameter :: ngo3      = 18
      integer, parameter :: ngpan     = 19
      integer, parameter :: ngoh      = 20
      integer, parameter :: ngho2     = 21
      integer, parameter :: ngno3     = 22
      integer, parameter :: ngch3o2   = 23
      integer, parameter :: ngch3o2h  = 24
      integer, parameter :: ngch3oh   = 25
      integer, parameter :: ngch3co3h = 26



      integer, parameter :: meqn1max = 20




      double precision, parameter :: dactiv = 0.7e-6       





      double precision, parameter :: avdiam = 20.e-6








      integer, parameter :: kiron = 1            










      double precision, parameter :: chlorine = 0.0






                                                    






      double precision, parameter :: caratio = 0.001





      double precision, parameter :: frac1 = 0.8               
      double precision, parameter :: frac2 = 0.2               








      double precision, parameter :: firon = 0.00003       
      double precision, parameter :: fman  = 0.00001       












	double precision, save :: wso2, wh2o2, whcho, whcooh, wnh3, whno3, whcl, wh2so4
	double precision, save :: wmol(29), amol(3), gmol(22)












      integer, parameter :: itol = 4
      integer, parameter :: itask = 1

      integer, parameter :: iopt = 1
      integer, parameter :: mf = 22
      integer, parameter :: worki = 100000             

      integer, parameter :: lrw1 = 22+9*meqn1max+2*meqn1max**2
      integer, parameter :: liw1 = 30+meqn1max



      double precision, parameter :: tola = 1.e-8             

      double precision, parameter :: tolr = 1.e-5             

      double precision, parameter :: workr = 300.0



































       integer, parameter :: numfunc = 7
       integer, parameter :: maxfev = 300*(numfunc+1)
       integer, parameter :: ml = numfunc - 1, mu = numfunc -1
       integer, parameter :: nprint = 0
       integer, parameter :: lr = numfunc*(numfunc+1)/2, ldfjac = numfunc
       integer, parameter :: mode = 2


       double precision, parameter :: xtol = 1.0e-3
       double precision, parameter :: epsfcn = 0.0e0, factor = 100.























































	integer, save :: maqurxn_all = 1
	integer, save :: maqurxn_sulf1 = 0
	integer, save :: mopt_eqrt_cons = 0
	integer, save :: mequlib_h2o2_ho2m = 0
	integer, save :: mgasrxn = 0
	integer, save :: mdiag_fullequil = 1
	integer, save :: mdiag_hybrd = 1
	integer, save :: mdiag_negconc = 1
	integer, save :: mdiag_rsrate = 1
	integer, save :: mdiag_svode = 1




      double precision, parameter :: rideal = 0.082058e0


      integer, parameter :: kaqx_siv   = 1
      integer, parameter :: kaqx_svi   = 2
      integer, parameter :: kaqx_no3m  = 4
      integer, parameter :: kaqx_h2o2  = 6
      integer, parameter :: kaqx_clm   = 15
      integer, parameter :: kaqx_nh4p  = 19
      integer, parameter :: kaqx_hso5m = 26
      integer, parameter :: kaqx_hmsa  = 27



      end module module_data_cmu_bulkaqchem

