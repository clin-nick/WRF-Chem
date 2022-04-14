










	module module_data_mosaic_other











	integer, parameter :: imaxd=1, jmaxd=1, kmaxd=100

	integer, parameter :: lmaxd=1200, l2maxd=1200 

	integer, parameter :: nsubareamaxd = 1



	integer, parameter :: k_pegbegin = 1


	integer, save :: khno3   = -999888777
	integer, save :: kh2so4  = -999888777
	integer, save :: knh3    = -999888777
	integer, save :: khcl    = -999888777
	integer, save :: kn2o5   = -999888777
	integer, save :: kclno2  = -999888777
	integer, save :: ko3     = -999888777
	integer, save :: kh2o    = -999888777
	integer, save :: ktemp   = -999888777
        integer, save :: kpcg1_b_c  = -999888777
        integer, save :: kpcg2_b_c  = -999888777
        integer, save :: kpcg3_b_c  = -999888777
        integer, save :: kpcg4_b_c  = -999888777
        integer, save :: kpcg5_b_c  = -999888777
        integer, save :: kpcg6_b_c  = -999888777
        integer, save :: kpcg7_b_c  = -999888777
        integer, save :: kpcg8_b_c  = -999888777
        integer, save :: kpcg9_b_c  = -999888777
        integer, save :: kpcg1_b_o  = -999888777
        integer, save :: kpcg2_b_o  = -999888777
        integer, save :: kpcg3_b_o  = -999888777
        integer, save :: kpcg4_b_o  = -999888777
        integer, save :: kpcg5_b_o  = -999888777
        integer, save :: kpcg6_b_o  = -999888777
        integer, save :: kpcg7_b_o  = -999888777
        integer, save :: kpcg8_b_o  = -999888777
        integer, save :: kpcg9_b_o  = -999888777
        integer, save :: kopcg1_b_c  = -999888777
        integer, save :: kopcg2_b_c  = -999888777
        integer, save :: kopcg3_b_c  = -999888777
        integer, save :: kopcg4_b_c  = -999888777
        integer, save :: kopcg5_b_c  = -999888777
        integer, save :: kopcg6_b_c  = -999888777
        integer, save :: kopcg7_b_c  = -999888777
        integer, save :: kopcg8_b_c  = -999888777
        integer, save :: kopcg1_b_o  = -999888777
        integer, save :: kopcg2_b_o  = -999888777
        integer, save :: kopcg3_b_o  = -999888777
        integer, save :: kopcg4_b_o  = -999888777
        integer, save :: kopcg5_b_o  = -999888777
        integer, save :: kopcg6_b_o  = -999888777
        integer, save :: kopcg7_b_o  = -999888777
        integer, save :: kopcg8_b_o  = -999888777
        integer, save :: kpcg1_f_c  = -999888777
        integer, save :: kpcg2_f_c  = -999888777
        integer, save :: kpcg3_f_c  = -999888777
        integer, save :: kpcg4_f_c  = -999888777
        integer, save :: kpcg5_f_c  = -999888777
        integer, save :: kpcg6_f_c  = -999888777
        integer, save :: kpcg7_f_c  = -999888777
        integer, save :: kpcg8_f_c  = -999888777
        integer, save :: kpcg9_f_c  = -999888777
        integer, save :: kpcg1_f_o  = -999888777
        integer, save :: kpcg2_f_o  = -999888777
        integer, save :: kpcg3_f_o  = -999888777
        integer, save :: kpcg4_f_o  = -999888777
        integer, save :: kpcg5_f_o  = -999888777
        integer, save :: kpcg6_f_o  = -999888777
        integer, save :: kpcg7_f_o  = -999888777
        integer, save :: kpcg8_f_o  = -999888777
        integer, save :: kpcg9_f_o  = -999888777
        integer, save :: kopcg1_f_c  = -999888777
        integer, save :: kopcg2_f_c  = -999888777
        integer, save :: kopcg3_f_c  = -999888777
        integer, save :: kopcg4_f_c  = -999888777
        integer, save :: kopcg5_f_c  = -999888777
        integer, save :: kopcg6_f_c  = -999888777
        integer, save :: kopcg7_f_c  = -999888777
        integer, save :: kopcg8_f_c  = -999888777
        integer, save :: kopcg1_f_o  = -999888777
        integer, save :: kopcg2_f_o  = -999888777
        integer, save :: kopcg3_f_o  = -999888777
        integer, save :: kopcg4_f_o  = -999888777
        integer, save :: kopcg5_f_o  = -999888777
        integer, save :: kopcg6_f_o  = -999888777
        integer, save :: kopcg7_f_o  = -999888777
        integer, save :: kopcg8_f_o  = -999888777
        integer, save :: ksmpa  = -999888777
        integer, save :: ksmpbb  = -999888777
        integer, save :: kant1_c  = -999888777
        integer, save :: kant2_c  = -999888777
        integer, save :: kant3_c  = -999888777
        integer, save :: kant4_c  = -999888777
        integer, save :: kant1_o  = -999888777
        integer, save :: kant2_o  = -999888777
        integer, save :: kant3_o  = -999888777
        integer, save :: kant4_o  = -999888777
        integer, save :: kbiog1_c  = -999888777
        integer, save :: kbiog2_c  = -999888777
        integer, save :: kbiog3_c  = -999888777
        integer, save :: kbiog4_c  = -999888777
        integer, save :: kbiog1_o  = -999888777
        integer, save :: kbiog2_o  = -999888777
        integer, save :: kbiog3_o  = -999888777
        integer, save :: kbiog4_o  = -999888777
        integer, save :: kasoaX = -999888777
        integer, save :: kasoa1 = -999888777
        integer, save :: kasoa2 = -999888777
        integer, save :: kasoa3 = -999888777
        integer, save :: kasoa4 = -999888777
        integer, save :: kbsoaX = -999888777
        integer, save :: kbsoa1 = -999888777
        integer, save :: kbsoa2 = -999888777
        integer, save :: kbsoa3 = -999888777
        integer, save :: kbsoa4 = -999888777

        integer, save :: kgly        = -999888777








	integer, save :: kso2    = -999888777
	integer, save :: kh2o2   = -999888777
	integer, save :: khcho   = -999888777
	integer, save :: khcooh  = -999888777
	integer, save :: koh     = -999888777
	integer, save :: kho2    = -999888777
	integer, save :: kno3    = -999888777
	integer, save :: kno     = -999888777
	integer, save :: kno2    = -999888777
	integer, save :: khono   = -999888777
	integer, save :: kpan    = -999888777
	integer, save :: kch3o2  = -999888777
	integer, save :: kch3oh  = -999888777
	integer, save :: kch3ooh = -999888777


	integer, save :: lunerr=-1, lunout=-1

	integer, save :: ltot=+999888777, ltot2=+999888777

	integer, save :: itot, jtot, ktot
	integer, save :: isvode, jsvode, ksvode, msvode
	integer, save :: iymdcur, ihmscur
	integer, save :: ncorecnt
	integer, save :: nsubareas


	real, parameter :: pi = 3.14159265

	real, save :: afracsubarea(kmaxd,nsubareamaxd)
	real, save :: cairclm(kmaxd)
	real, save :: ptotclm(kmaxd)
	real, save :: rclm(kmaxd,l2maxd)
	real, save :: relhumclm(kmaxd)
	real, save :: rcldwtr_sub(kmaxd,nsubareamaxd)
	real, save :: rsub(l2maxd,kmaxd,nsubareamaxd)
	real, save :: t


	character(len=20), save :: name(l2maxd)




	integer, save :: aboxtest_testmode = 0
	integer, save :: aboxtest_units_convert = 1
	integer, save :: aboxtest_rh_method = 1
	integer, save :: aboxtest_map_method = 1
	integer, save :: aboxtest_gases_fixed = 0

	real, save :: aboxtest_min_temp = 233.0
	real, save :: aboxtest_min_relhum = 0.05
	real, save :: aboxtest_max_relhum = 0.98


	end module module_data_mosaic_other
