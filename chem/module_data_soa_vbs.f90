MODULE module_data_soa_vbs






      IMPLICIT NONE
      INTEGER NP                
      PARAMETER (NP = 8)



      INTEGER MAXITS            
      PARAMETER (MAXITS = 100)

      REAL TOLF                 
      PARAMETER (TOLF = 1.E-09)

      REAL TOLMIN                 
      PARAMETER (TOLMIN = 1.E-12) 

      REAL TOLX                 
      PARAMETER (TOLX = 1.E-10)

      REAL STPMX                
      PARAMETER (STPMX = 100.)

      REAL c303, c302
      PARAMETER (c303=19.83, c302=5417.4)

      INTEGER lcva, lcvb, lspcv, ldesn
      PARAMETER (lcva=4,lcvb=4, lspcv=lcva+lcvb)
      PARAMETER (ldesn=13)




      INTEGER laerdvc, lnonaerdvc, l1ae, laero, imodes, aspec

       PARAMETER (laerdvc=46,lnonaerdvc=17+lspcv)


      PARAMETER (l1ae=laerdvc+lnonaerdvc)
      PARAMETER (laero=4,imodes=4,aspec=1)










      INTEGER aemiss
      PARAMETER (aemiss=4)


 
      INTEGER, PARAMETER :: ldroga=6    
      INTEGER, PARAMETER :: ldrogb=3    
      INTEGER, PARAMETER :: ldrogr=1    
      INTEGER, PARAMETER :: ldrog_vbs=ldroga+ldrogb+ldrogr 




























                            
      INTEGER orgaer
                            


      PARAMETER (orgaer=2)



      INTEGER n_ae_vis_spc
      PARAMETER (n_ae_vis_spc=2)


      INTEGER idcvw
      PARAMETER (idcvw=1)

      INTEGER ibext
      PARAMETER (ibext=2)




      INTEGER vso4aj
      PARAMETER (vso4aj=1)


      INTEGER vso4ai
      PARAMETER (vso4ai=2)


      INTEGER vnh4aj
      PARAMETER (vnh4aj=3)


      INTEGER vnh4ai
      PARAMETER (vnh4ai=4)


      INTEGER vno3aj
      PARAMETER (vno3aj=5)


      INTEGER vno3ai
      PARAMETER (vno3ai=6)


      INTEGER vnaaj
      PARAMETER (vnaaj=7)


      INTEGER vnaai
      PARAMETER (vnaai=8)


      INTEGER vclaj
      PARAMETER (vclaj=9)


      INTEGER vclai
      PARAMETER (vclai=10)




                INTEGER vcaaj
                PARAMETER (vcaaj=11)

                INTEGER vcaai
                PARAMETER (vcaai=12)


                INTEGER vkaj
                PARAMETER (vkaj=13)

                INTEGER vkai
                PARAMETER (vkai=14)


                INTEGER vmgaj
                PARAMETER (vmgaj=15)

                INTEGER vmgai
                PARAMETER (vmgai=16)



      INTEGER, PARAMETER ::  vasoa1j=17
      INTEGER, PARAMETER ::  vasoa1i=18

      INTEGER, PARAMETER ::  vasoa2j=19
      INTEGER, PARAMETER ::  vasoa2i=20

      INTEGER, PARAMETER ::  vasoa3j=21
      INTEGER, PARAMETER ::  vasoa3i=22

      INTEGER, PARAMETER ::  vasoa4j=23
      INTEGER, PARAMETER ::  vasoa4i=24


      INTEGER, PARAMETER ::  vbsoa1j=25
      INTEGER, PARAMETER ::  vbsoa1i=26

      INTEGER, PARAMETER ::  vbsoa2j=27
      INTEGER, PARAMETER ::  vbsoa2i=28

      INTEGER, PARAMETER ::  vbsoa3j=29
      INTEGER, PARAMETER ::  vbsoa3i=30

      INTEGER, PARAMETER ::  vbsoa4j=31
      INTEGER, PARAMETER ::  vbsoa4i=32



      INTEGER vorgpaj
      PARAMETER (vorgpaj=33)


      INTEGER vorgpai
      PARAMETER (vorgpai=34)


      INTEGER vecj
      PARAMETER (vecj=35)


      INTEGER veci
      PARAMETER (veci=36)


      INTEGER vp25aj
      PARAMETER (vp25aj=37)


      INTEGER vp25ai
      PARAMETER (vp25ai=38)


      INTEGER vantha
      PARAMETER (vantha=39)


      INTEGER vseas
      PARAMETER (vseas=40)


      INTEGER vsoila
      PARAMETER (vsoila=41)


      INTEGER vnu0
      PARAMETER (vnu0=42)


      INTEGER vac0
      PARAMETER (vac0=43)


      INTEGER vcorn
      PARAMETER (vcorn=44)


      INTEGER vh2oaj
      PARAMETER (vh2oaj=45)


      INTEGER vh2oai
      PARAMETER (vh2oai=46)


      INTEGER vnu3
      PARAMETER (vnu3=47)

      INTEGER vac3
      PARAMETER (vac3=48)


      INTEGER vcor3
      PARAMETER (vcor3=49)


      INTEGER vsulf
      PARAMETER (vsulf=50)


      INTEGER vhno3
      PARAMETER (vhno3=51)


      INTEGER vnh3
      PARAMETER (vnh3=52)


      INTEGER vhcl
      PARAMETER (vhcl=53)



        INTEGER vn2o5
        PARAMETER (vn2o5=54)

        INTEGER vclno2
        PARAMETER (vclno2=55)


        INTEGER vgamn2o5
        PARAMETER (vgamn2o5=56)


        INTEGER vcn2o5
        PARAMETER (vcn2o5=57)


        INTEGER vkn2o5
        PARAMETER (vkn2o5=58)


        INTEGER vyclno2
        PARAMETER (vyclno2=59)


        INTEGER vsnu
        PARAMETER (vsnu=60)


        INTEGER vsac
        PARAMETER (vsac=61)


        INTEGER vsco
        PARAMETER (vsco=62)




        INTEGER valt_in
        PARAMETER (valt_in=63)

INTEGER, PARAMETER :: vcvasoa1=64
INTEGER, PARAMETER :: vcvasoa2=65
INTEGER, PARAMETER :: vcvasoa3=66
INTEGER, PARAMETER :: vcvasoa4=67
INTEGER, PARAMETER :: vcvbsoa1=68
INTEGER, PARAMETER :: vcvbsoa2=69
INTEGER, PARAMETER :: vcvbsoa3=70
INTEGER, PARAMETER :: vcvbsoa4=71







      INTEGER naspcssed
      PARAMETER (naspcssed=6)


      INTEGER vsnnuc
      PARAMETER (vsnnuc=1)


      INTEGER vsnacc
      PARAMETER (vsnacc=2)


      INTEGER vsncor
      PARAMETER (vsncor=3)


      INTEGER vsmnuc
      PARAMETER (vsmnuc=4)


      INTEGER vsmacc
      PARAMETER (vsmacc=5)


      INTEGER vsmcor
      PARAMETER (vsmcor=6)





      INTEGER naspcsdep
      PARAMETER (naspcsdep=7)


      INTEGER vdnnuc
      PARAMETER (vdnnuc=1)


      INTEGER vdnacc
      PARAMETER (vdnacc=2)


      INTEGER vdncor
      PARAMETER (vdncor=3)


      INTEGER vdmnuc
      PARAMETER (vdmnuc=4)


      INTEGER vdmacc
      PARAMETER (vdmacc=5)


      INTEGER vdmfine
      PARAMETER (vdmfine=6)


      INTEGER vdmcor
      PARAMETER (vdmcor=7)



INTEGER, PARAMETER :: palk4=1
INTEGER, PARAMETER :: palk5=2
INTEGER, PARAMETER :: pole1=3
INTEGER, PARAMETER :: pole2=4
INTEGER, PARAMETER :: paro1=5
INTEGER, PARAMETER :: paro2=6


INTEGER, PARAMETER :: pisop=7
INTEGER, PARAMETER :: pterp=8
INTEGER, PARAMETER :: psesq=9


INTEGER, PARAMETER :: pbrch=10

 
INTEGER, PARAMETER :: pasoa1=1
INTEGER, PARAMETER :: pasoa2=2
INTEGER, PARAMETER :: pasoa3=3
INTEGER, PARAMETER :: pasoa4=4
      
INTEGER, PARAMETER :: pbsoa1=5
INTEGER, PARAMETER :: pbsoa2=6
INTEGER, PARAMETER :: pbsoa3=7
INTEGER, PARAMETER :: pbsoa4=8





















































      REAL*8 & 
        pirs
      PARAMETER (pirs=3.14159265358979324)





      REAL avo
      PARAMETER (avo=6.0221367E23)


      REAL rgasuniv
      PARAMETER (rgasuniv=8.314510)


      REAL stdatmpa
      PARAMETER (stdatmpa=101325.0)


      REAL stdtemp
      PARAMETER (stdtemp=273.15)


      REAL stfblz
      PARAMETER (stfblz=5.67051E-8)



      REAL grav
      PARAMETER (grav=9.80622)



      REAL molvol
      PARAMETER (molvol=22.41410)





      REAL mwair
                        

      PARAMETER (mwair=28.9628)


      REAL rdgas
      PARAMETER (rdgas=1.0E3*rgasuniv/mwair)


      REAL threepi
      PARAMETER (threepi=3.0*pirs)


      REAL f6dpi
      PARAMETER (f6dpi=6.0/pirs)


      REAL f6dpi9
      PARAMETER (f6dpi9=1.0E9*f6dpi)


      REAL f6dpim9
      PARAMETER (f6dpim9=1.0E-9*f6dpi)


      REAL sqrtpi
      PARAMETER (sqrtpi=1.7724539)


      REAL sqrt2
      PARAMETER (sqrt2=1.4142135623731)


      REAL lgsqt2
      PARAMETER (lgsqt2=0.34657359027997)


      REAL dlgsqt2
      PARAMETER (dlgsqt2=1.0/lgsqt2)


      REAL one3
      PARAMETER (one3=1.0/3.0)


      REAL two3
      PARAMETER (two3=2.0/3.0)





      REAL boltz
      PARAMETER (boltz=rgasuniv/avo)






      REAL rhoso4
      PARAMETER (rhoso4=1.8E3)


      REAL rhonh4
      PARAMETER (rhonh4=1.8E3)


      REAL rhono3
      PARAMETER (rhono3=1.8E3)


      REAL rhoh2o
      PARAMETER (rhoh2o=1.0E3)


      REAL rhoorg
      PARAMETER (rhoorg=1.0E3)


      REAL rhosoil
      PARAMETER (rhosoil=2.6E3)


      REAL rhoseas
      PARAMETER (rhoseas=2.2E3)


      REAL rhoanth
      PARAMETER (rhoanth=2.2E3)


      REAL rhona
      PARAMETER (rhona=2.2E3)


      REAL rhocl
      PARAMETER (rhocl=2.2E3)



                REAL rhoca
                PARAMETER (rhoca=2.6E3)


                REAL rhok
                PARAMETER (rhok=2.6E3)


                REAL rhomg
                PARAMETER (rhomg=2.6E3)




      REAL so4fac
      PARAMETER (so4fac=f6dpim9/rhoso4)

      REAL nh4fac
      PARAMETER (nh4fac=f6dpim9/rhonh4)

      REAL h2ofac
      PARAMETER (h2ofac=f6dpim9/rhoh2o)

      REAL no3fac
      PARAMETER (no3fac=f6dpim9/rhono3)

      REAL orgfac
      PARAMETER (orgfac=f6dpim9/rhoorg)

      REAL soilfac
      PARAMETER (soilfac=f6dpim9/rhosoil)

      REAL seasfac
      PARAMETER (seasfac=f6dpim9/rhoseas)

      REAL anthfac
      PARAMETER (anthfac=f6dpim9/rhoanth)

      REAL nafac
      PARAMETER (nafac=f6dpim9/rhona)

      REAL clfac
      PARAMETER (clfac=f6dpim9/rhocl)


                REAL cafac
                PARAMETER (cafac=f6dpim9/rhoca)

                REAL kfac
                PARAMETER (kfac=f6dpim9/rhok)

                REAL mgfac
                PARAMETER (mgfac=f6dpim9/rhomg)




      REAL pss0
      PARAMETER (pss0=101325.0)


      REAL tss0
      PARAMETER (tss0=288.15)


      REAL sginin
      PARAMETER (sginin=1.70)


      REAL sginia
      PARAMETER (sginia=2.00)


      REAL sginic
      PARAMETER (sginic=2.5)


      REAL dginin
      PARAMETER (dginin=0.01E-6)


      REAL dginia
      PARAMETER (dginia=0.07E-6)


      REAL dginic
      PARAMETER (dginic=1.0E-6)









      REAL en1

      REAL ea1

      REAL ec1


      REAL esn04

      REAL esa04

      REAL esc04


      REAL esn05

      REAL esa05


      REAL esn08

      REAL esa08

      REAL esc08


      REAL esn09

      REAL esa09


      REAL esn12

      REAL esa12

      REAL esc12


      REAL esn16

      REAL esa16

      REAL esc16


      REAL esn20

      REAL esa20

      REAL esc20


      REAL esn25

      REAL esa25


      REAL esn24

      REAL esa24

      REAL esc24


      REAL esn28

      REAL esa28

      REAL esc28


      REAL esn32

      REAL esa32

      REAL esc32


      REAL esn36

      REAL esa36

      REAL esc36


      REAL esn49

      REAL esa49


      REAL esn52

      REAL esa52


      REAL esn64

      REAL esa64

      REAL esc64


      REAL esn100


      REAL esnm20

      REAL esam20

      REAL escm20


      REAL esnm32

      REAL esam32

      REAL escm32


      REAL xxlsgn

      REAL xxlsga

      REAL xxlsgc


      REAL l2sginin

      REAL l2sginia

      REAL l2sginic







                            
      INTEGER inucl
                            
                            

      PARAMETER (inucl=2)



      LOGICAL icoarse
      PARAMETER (icoarse=.FALSE.) 





      REAL dgvem_i
      PARAMETER (dgvem_i=0.03E-6) 
      REAL sgem_i
      PARAMETER (sgem_i=1.7)


      REAL dgvem_j
      PARAMETER (dgvem_j=0.3E-6) 
      REAL sgem_j
      PARAMETER (sgem_j=2.0)


      REAL dgvem_c
      PARAMETER (dgvem_c=6.0E-6) 
      REAL sgem_c
      PARAMETER (sgem_c=2.2)



      REAL factnumn

      REAL factnuma

      REAL factnumc

      REAL facatkn_min, facacc_min
      PARAMETER (facatkn_min=0.04,facacc_min=1.0-facatkn_min)
      REAL xxm3
      REAL, PARAMETER ::  conmin = 1.E-16
      REAL, PARAMETER ::  epsilc = 1.E-16

      REAL*8 & 
        nummin_i
      REAL*8 & 
        nummin_j
      REAL*8 & 
        nummin_c
















      REAL alphsulf
      PARAMETER (alphsulf=1.0) 


      REAL mwh2so4
      PARAMETER (mwh2so4=98.07354E-3) 


      REAL diffsulf
      PARAMETER (diffsulf=9.362223E-06) 


      REAL alphaorg
      PARAMETER (alphaorg=1.0)                                    


      REAL mworg
      PARAMETER (mworg=175.0E-03)






      REAL difforg
      PARAMETER (difforg=5.151174E-06)


      REAL cconc
      PARAMETER (cconc=2.0*pirs*diffsulf) 


      REAL cconc_org
      PARAMETER (cconc_org=2.0*pirs*difforg) 


      REAL ccofm_org






      REAL ccofm

      REAL aeroconcmin
      PARAMETER (aeroconcmin=0.0001) 



































































	integer, parameter :: maxd_atype = 2
	integer, parameter :: maxd_asize = 2
	integer, parameter :: maxd_acomp = 19
	integer, parameter :: maxd_aphase = 2
	integer, save :: ai_phase 
	integer, save :: cw_phase 
	integer, save :: ci_phase 
	integer, save :: cr_phase 
	integer, save :: cs_phase 
	integer, save :: cg_phase 

	integer, save :: ntype_aer = 0 
	integer, save :: ntot_mastercomp_aer = 0 
	integer, save :: nphase_aer = 0 

	integer, save ::   &
      	  msectional, maerosolincw,   &
      	  nsize_aer( maxd_atype ),   & 
      	  ncomp_aer( maxd_atype ),   & 
      	  ncomp_aer_nontracer( maxd_atype ),   &
          mastercompptr_aer(maxd_acomp, maxd_atype), &   
      	  massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase ), & 
      	  waterptr_aer( maxd_asize, maxd_atype ), & 
      	  hygroptr_aer( maxd_asize, maxd_atype ), & 
      	  numptr_aer( maxd_asize, maxd_atype, maxd_aphase ), & 
          mprognum_aer(maxd_asize,maxd_atype,maxd_aphase)

	real, save ::   &
          dens_aer( maxd_acomp, maxd_atype ),   &
          dens_mastercomp_aer( maxd_acomp ),   &
      	  mw_mastercomp_aer( maxd_acomp ), &
      	  mw_aer( maxd_acomp, maxd_atype ),  &
      	  hygro_mastercomp_aer( maxd_acomp ), &
      	  hygro_aer( maxd_acomp, maxd_atype )
	character*10, save ::   &
      	  name_mastercomp_aer( maxd_acomp ), &
      	  name_aer( maxd_acomp, maxd_atype )

	real, save ::   &
          volumcen_sect( maxd_asize, maxd_atype ),   &
          volumlo_sect( maxd_asize, maxd_atype ),   &
          volumhi_sect( maxd_asize, maxd_atype ),   &
          dcen_sect( maxd_asize, maxd_atype ),   &
          dlo_sect( maxd_asize, maxd_atype ),   &
          dhi_sect( maxd_asize, maxd_atype ),   &
	  sigmag_aer(maxd_asize, maxd_atype)

	integer, save ::                     &
      	  lptr_so4_aer(maxd_asize,maxd_atype,maxd_aphase),        &
      	  lptr_nh4_aer(maxd_asize,maxd_atype,maxd_aphase),        &
      	  lptr_no3_aer(maxd_asize,maxd_atype,maxd_aphase),        &

       	  lptr_asoa1_aer(maxd_asize,maxd_atype,maxd_aphase),    &
      	  lptr_asoa2_aer(maxd_asize,maxd_atype,maxd_aphase),    &
      	  lptr_asoa3_aer(maxd_asize,maxd_atype,maxd_aphase),     &
      	  lptr_asoa4_aer(maxd_asize,maxd_atype,maxd_aphase),     &
      	  lptr_bsoa1_aer(maxd_asize,maxd_atype,maxd_aphase),     &
      	  lptr_bsoa2_aer(maxd_asize,maxd_atype,maxd_aphase),     &
      	  lptr_bsoa3_aer(maxd_asize,maxd_atype,maxd_aphase),     &
      	  lptr_bsoa4_aer(maxd_asize,maxd_atype,maxd_aphase),     &










      	  lptr_orgpa_aer(maxd_asize,maxd_atype,maxd_aphase),      &
      	  lptr_ec_aer(maxd_asize,maxd_atype,maxd_aphase),         &
      	  lptr_p25_aer(maxd_asize,maxd_atype,maxd_aphase),        &
          lptr_anth_aer(maxd_asize,maxd_atype,maxd_aphase),       &
      	  lptr_cl_aer(maxd_asize,maxd_atype,maxd_aphase),         &
      	  lptr_na_aer(maxd_asize,maxd_atype,maxd_aphase),         &
      	  lptr_seas_aer(maxd_asize,maxd_atype,maxd_aphase),       &
      	  lptr_soil_aer(maxd_asize,maxd_atype,maxd_aphase)

	logical, save ::                     &
      	  do_cloudchem_aer(maxd_asize,maxd_atype)



	real, parameter :: mw_so4_aer   = 96.066
	real, parameter :: mw_no3_aer   = 62.007
	real, parameter :: mw_nh4_aer   = 18.042
	real, parameter :: mw_oc_aer    = 250.0
	real, parameter :: mw_ec_aer    = 1.0
	real, parameter :: mw_oin_aer   = 1.0
	real, parameter :: mw_dust_aer  = 100.087
	real, parameter :: mw_seas_aer  = 58.440
	real, parameter :: mw_cl_aer    = 35.450
	real, parameter :: mw_na_aer    = 22.990
	real, parameter :: mw_water_aer = 18.016


	real, parameter :: dens_so4_aer  = 1.80   
	real, parameter :: dens_no3_aer  = 1.80   
	real, parameter :: dens_nh4_aer  = 1.80   
	real, parameter :: dens_oc_aer   = 1.5    
	real, parameter :: dens_ec_aer   = 1.70
	real, parameter :: dens_dust_aer = 2.60  
	real, parameter :: dens_oin_aer  = 2.20  
	real, parameter :: dens_seas_aer = 2.20  
	real, parameter :: dens_cl_aer   = 2.20
	real, parameter :: dens_na_aer   = 2.20


	real, parameter :: dens_water_aer  = 1.0


	real, parameter :: hygro_so4_aer  = 0.5
	real, parameter :: hygro_no3_aer  = 0.5
	real, parameter :: hygro_nh4_aer  = 0.5
	real, parameter :: hygro_oc_aer   = 0.14
	real, parameter :: hygro_ec_aer   = 1.e-6
	real, parameter :: hygro_oin_aer  = 0.14
	real, parameter :: hygro_dust_aer = 0.1
	real, parameter :: hygro_seas_aer = 1.16
	real, parameter :: hygro_cl_aer   = 1.16
	real, parameter :: hygro_na_aer   = 1.16


	real dlndg_nimptblgrow
	integer nimptblgrow_mind, nimptblgrow_maxd
	parameter (nimptblgrow_mind=-14, nimptblgrow_maxd=24)
     	real scavimptblnum(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype), &
     	     scavimptblvol(4, nimptblgrow_mind:nimptblgrow_maxd, maxd_asize, maxd_atype)


      INTEGER NGAUSdv
      PARAMETER( NGAUSdv = 7 )  
      REAL Y_GQ(NGAUSdv), WGAUS(NGAUSdv)








END Module module_data_soa_vbs
