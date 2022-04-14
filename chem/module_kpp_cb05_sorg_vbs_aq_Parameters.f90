























MODULE cb05_sorg_vbs_aq_Parameters

  USE cb05_sorg_vbs_aq_Precision
  PUBLIC
  SAVE



  INTEGER, PARAMETER :: NSPEC = 107 

  INTEGER, PARAMETER :: NVAR = 105 

  INTEGER, PARAMETER :: NVARACT = 73 

  INTEGER, PARAMETER :: NFIX = 2 

  INTEGER, PARAMETER :: NREACT = 198 

  INTEGER, PARAMETER :: NVARST = 1 

  INTEGER, PARAMETER :: NFIXST = 106 

  INTEGER, PARAMETER :: NONZERO = 759 

  INTEGER, PARAMETER :: LU_NONZERO = 840 

  INTEGER, PARAMETER :: CNVAR = 106 

  INTEGER, PARAMETER :: NLOOKAT = 0 

  INTEGER, PARAMETER :: NMONITOR = 0 

  INTEGER, PARAMETER :: NMASS = 1 

  REAL(kind=dp), PARAMETER :: PI = 3.14159265358979 




  INTEGER, PARAMETER :: ind_tolaer1 = 1 
  INTEGER, PARAMETER :: ind_tolaer2 = 2 
  INTEGER, PARAMETER :: ind_cslaer = 3 
  INTEGER, PARAMETER :: ind_xylaer1 = 4 
  INTEGER, PARAMETER :: ind_xylaer2 = 5 
  INTEGER, PARAMETER :: ind_isoaer1 = 6 
  INTEGER, PARAMETER :: ind_isoaer2 = 7 
  INTEGER, PARAMETER :: ind_sulf = 8 
  INTEGER, PARAMETER :: ind_sulaer = 9 
  INTEGER, PARAMETER :: ind_terpaer = 10 
  INTEGER, PARAMETER :: ind_humaer = 11 
  INTEGER, PARAMETER :: ind_hum = 12 
  INTEGER, PARAMETER :: ind_limaer1 = 13 
  INTEGER, PARAMETER :: ind_limaer2 = 14 
  INTEGER, PARAMETER :: ind_lim = 15 
  INTEGER, PARAMETER :: ind_ociaer1 = 16 
  INTEGER, PARAMETER :: ind_ociaer2 = 17 
  INTEGER, PARAMETER :: ind_oci = 18 
  INTEGER, PARAMETER :: ind_apinaer1 = 19 
  INTEGER, PARAMETER :: ind_apinaer2 = 20 
  INTEGER, PARAMETER :: ind_apinaer3 = 21 
  INTEGER, PARAMETER :: ind_apinaer4 = 22 
  INTEGER, PARAMETER :: ind_apin = 23 
  INTEGER, PARAMETER :: ind_bpinaer1 = 24 
  INTEGER, PARAMETER :: ind_bpinaer2 = 25 
  INTEGER, PARAMETER :: ind_bpinaer3 = 26 
  INTEGER, PARAMETER :: ind_bpinaer4 = 27 
  INTEGER, PARAMETER :: ind_bpinaer5 = 28 
  INTEGER, PARAMETER :: ind_bpin = 29 
  INTEGER, PARAMETER :: ind_teraer1 = 30 
  INTEGER, PARAMETER :: ind_teraer2 = 31 
  INTEGER, PARAMETER :: ind_ter = 32 
  INTEGER, PARAMETER :: ind_alkhaer1 = 33 
  INTEGER, PARAMETER :: ind_alkh = 34 
  INTEGER, PARAMETER :: ind_pahaer1 = 35 
  INTEGER, PARAMETER :: ind_pahaer2 = 36 
  INTEGER, PARAMETER :: ind_pah = 37 
  INTEGER, PARAMETER :: ind_hg2 = 38 
  INTEGER, PARAMETER :: ind_CVASOA1 = 39 
  INTEGER, PARAMETER :: ind_CVASOA2 = 40 
  INTEGER, PARAMETER :: ind_CVASOA3 = 41 
  INTEGER, PARAMETER :: ind_CVASOA4 = 42 
  INTEGER, PARAMETER :: ind_CVBSOA1 = 43 
  INTEGER, PARAMETER :: ind_CVBSOA2 = 44 
  INTEGER, PARAMETER :: ind_CVBSOA3 = 45 
  INTEGER, PARAMETER :: ind_CVBSOA4 = 46 
  INTEGER, PARAMETER :: ind_cl2 = 47 
  INTEGER, PARAMETER :: ind_so2 = 48 
  INTEGER, PARAMETER :: ind_pan = 49 
  INTEGER, PARAMETER :: ind_hocl = 50 
  INTEGER, PARAMETER :: ind_tol = 51 
  INTEGER, PARAMETER :: ind_h2 = 52 
  INTEGER, PARAMETER :: ind_o1d = 53 
  INTEGER, PARAMETER :: ind_n2o5 = 54 
  INTEGER, PARAMETER :: ind_xyl = 55 
  INTEGER, PARAMETER :: ind_ch4 = 56 
  INTEGER, PARAMETER :: ind_hg0 = 57 
  INTEGER, PARAMETER :: ind_hono = 58 
  INTEGER, PARAMETER :: ind_facd = 59 
  INTEGER, PARAMETER :: ind_to2 = 60 
  INTEGER, PARAMETER :: ind_pacd = 61 
  INTEGER, PARAMETER :: ind_aacd = 62 
  INTEGER, PARAMETER :: ind_panx = 63 
  INTEGER, PARAMETER :: ind_pna = 64 
  INTEGER, PARAMETER :: ind_etha = 65 
  INTEGER, PARAMETER :: ind_meoh = 66 
  INTEGER, PARAMETER :: ind_etoh = 67 
  INTEGER, PARAMETER :: ind_h2o2 = 68 
  INTEGER, PARAMETER :: ind_hco3 = 69 
  INTEGER, PARAMETER :: ind_hcl = 70 
  INTEGER, PARAMETER :: ind_cro = 71 
  INTEGER, PARAMETER :: ind_rooh = 72 
  INTEGER, PARAMETER :: ind_clo = 73 
  INTEGER, PARAMETER :: ind_mgly = 74 
  INTEGER, PARAMETER :: ind_fmcl = 75 
  INTEGER, PARAMETER :: ind_mepx = 76 
  INTEGER, PARAMETER :: ind_hno3 = 77 
  INTEGER, PARAMETER :: ind_co = 78 
  INTEGER, PARAMETER :: ind_open = 79 
  INTEGER, PARAMETER :: ind_ror = 80 
  INTEGER, PARAMETER :: ind_eth = 81 
  INTEGER, PARAMETER :: ind_cres = 82 
  INTEGER, PARAMETER :: ind_terp = 83 
  INTEGER, PARAMETER :: ind_iole = 84 
  INTEGER, PARAMETER :: ind_ole = 85 
  INTEGER, PARAMETER :: ind_xo2n = 86 
  INTEGER, PARAMETER :: ind_par = 87 
  INTEGER, PARAMETER :: ind_isop = 88 
  INTEGER, PARAMETER :: ind_ispd = 89 
  INTEGER, PARAMETER :: ind_aldx = 90 
  INTEGER, PARAMETER :: ind_ald2 = 91 
  INTEGER, PARAMETER :: ind_form = 92 
  INTEGER, PARAMETER :: ind_ntr = 93 
  INTEGER, PARAMETER :: ind_o3 = 94 
  INTEGER, PARAMETER :: ind_no3 = 95 
  INTEGER, PARAMETER :: ind_xo2 = 96 
  INTEGER, PARAMETER :: ind_c2o3 = 97 
  INTEGER, PARAMETER :: ind_oh = 98 
  INTEGER, PARAMETER :: ind_no = 99 
  INTEGER, PARAMETER :: ind_cl = 100 
  INTEGER, PARAMETER :: ind_o = 101 
  INTEGER, PARAMETER :: ind_ho2 = 102 
  INTEGER, PARAMETER :: ind_cxo3 = 103 
  INTEGER, PARAMETER :: ind_meo2 = 104 
  INTEGER, PARAMETER :: ind_no2 = 105 




  INTEGER, PARAMETER :: ind_H2O = 106 
  INTEGER, PARAMETER :: ind_M = 107 




  INTEGER, PARAMETER :: indf_H2O = 1 
  INTEGER, PARAMETER :: indf_M = 2 

END MODULE cb05_sorg_vbs_aq_Parameters

