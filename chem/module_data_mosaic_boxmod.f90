      module module_data_mosaic_boxmod








      use module_data_mosaic_constants, only:  &
          avogad, deg2rad, pi, piover4, piover6, third

      implicit none


      integer, save :: idiag_sect_coag = 0

      integer, save :: idiag_sect_movesect = 0

      integer, save :: idiag_sect_newnuc = 0

      character(len=20), allocatable :: name_rbox(:)









      integer, save :: kh2so4     = -999888777
      integer, save :: khno3      = -999888777
      integer, save :: knh3       = -999888777
      integer, save :: khcl       = -999888777
      integer, save :: kmsa       = -999888777
      integer, save :: karo1      = -999888777
      integer, save :: karo2      = -999888777
      integer, save :: kalk1      = -999888777
      integer, save :: kole1      = -999888777
      integer, save :: kapi1      = -999888777
      integer, save :: kapi2      = -999888777
      integer, save :: klim1      = -999888777
      integer, save :: klim2      = -999888777
      integer, save :: kn2o5      = -999888777
      integer, save :: kclno2     = -999888777
      integer, save :: ko3        = -999888777
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
      integer, save :: kopcg1_b_c = -999888777
      integer, save :: kopcg2_b_c = -999888777
      integer, save :: kopcg3_b_c = -999888777
      integer, save :: kopcg4_b_c = -999888777
      integer, save :: kopcg5_b_c = -999888777
      integer, save :: kopcg6_b_c = -999888777
      integer, save :: kopcg7_b_c = -999888777
      integer, save :: kopcg8_b_c = -999888777
      integer, save :: kopcg1_b_o = -999888777
      integer, save :: kopcg2_b_o = -999888777
      integer, save :: kopcg3_b_o = -999888777
      integer, save :: kopcg4_b_o = -999888777
      integer, save :: kopcg5_b_o = -999888777
      integer, save :: kopcg6_b_o = -999888777
      integer, save :: kopcg7_b_o = -999888777
      integer, save :: kopcg8_b_o = -999888777
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
      integer, save :: kopcg1_f_c = -999888777
      integer, save :: kopcg2_f_c = -999888777
      integer, save :: kopcg3_f_c = -999888777
      integer, save :: kopcg4_f_c = -999888777
      integer, save :: kopcg5_f_c = -999888777
      integer, save :: kopcg6_f_c = -999888777
      integer, save :: kopcg7_f_c = -999888777
      integer, save :: kopcg8_f_c = -999888777
      integer, save :: kopcg1_f_o = -999888777
      integer, save :: kopcg2_f_o = -999888777
      integer, save :: kopcg3_f_o = -999888777
      integer, save :: kopcg4_f_o = -999888777
      integer, save :: kopcg5_f_o = -999888777
      integer, save :: kopcg6_f_o = -999888777
      integer, save :: kopcg7_f_o = -999888777
      integer, save :: kopcg8_f_o = -999888777
      integer, save :: ksmpa      = -999888777
      integer, save :: ksmpbb     = -999888777


      integer, save :: kant1_c    = -999888777
      integer, save :: kant2_c    = -999888777
      integer, save :: kant3_c    = -999888777
      integer, save :: kant4_c    = -999888777
      integer, save :: kant1_o    = -999888777
      integer, save :: kant2_o    = -999888777
      integer, save :: kant3_o    = -999888777
      integer, save :: kant4_o    = -999888777
      integer, save :: kbiog1_c   = -999888777
      integer, save :: kbiog2_c   = -999888777
      integer, save :: kbiog3_c   = -999888777
      integer, save :: kbiog4_c   = -999888777
      integer, save :: kbiog1_o   = -999888777
      integer, save :: kbiog2_o   = -999888777
      integer, save :: kbiog3_o   = -999888777
      integer, save :: kbiog4_o   = -999888777
      integer, save :: kso2       = -999888777
      integer, save :: kh2o2      = -999888777
      integer, save :: khcho      = -999888777
      integer, save :: khcooh     = -999888777
      integer, save :: koh        = -999888777
      integer, save :: kho2       = -999888777
      integer, save :: kno3       = -999888777
      integer, save :: kno        = -999888777
      integer, save :: kno2       = -999888777
      integer, save :: khono      = -999888777
      integer, save :: kpan       = -999888777
      integer, save :: kch3o2     = -999888777
      integer, save :: kch3oh     = -999888777
      integer, save :: kch3ooh    = -999888777








































      end module module_data_mosaic_boxmod
