MODULE module_ctrans_aqchem

CONTAINS
























      SUBROUTINE AQCHEM ( TEMP, PRES_PA, TAUCLD, PRCRATE, &
                          WCAVG, WTAVG, AIRM, ALFA0, ALFA2, ALFA3, GAS, &
                          AEROSOL, LIQUID, GASWDEP, AERWDEP, HPWDEP )
























































































































      IMPLICIT NONE






































      REAL        PI 
      PARAMETER ( PI = 3.14159265358979324 )
 
      REAL        PI180 
      PARAMETER ( PI180  = PI / 180.0 )


 
      REAL        REARTH 
                         
                         

      PARAMETER ( REARTH = 6370000.0 )  
 
      REAL        SIDAY 
                        
      PARAMETER ( SIDAY = 86164.09 )
 
      REAL        GRAV 
                       
                       
      PARAMETER ( GRAV = 9.80622 )

      REAL        DG2M 
      PARAMETER ( DG2M = REARTH * PI180 )


      REAL        SOLCNST 
      PARAMETER ( SOLCNST = 1373.0 )



      REAL        AVO 
      PARAMETER ( AVO = 6.0221367E23 )

      REAL        RGASUNIV 
      PARAMETER ( RGASUNIV = 8.314510 )

      REAL        STDATMPA 
      PARAMETER ( STDATMPA = 101325.0 )

      REAL        STDTEMP 
      PARAMETER ( STDTEMP = 273.15 )

      REAL        STFBLZ 
      PARAMETER ( STFBLZ = 5.67051E-8 ) 



      REAL        MOLVOL 
      PARAMETER ( MOLVOL = 22.41410 ) 



      REAL        MWAIR 
                        
                        
      PARAMETER ( MWAIR = 28.9628 )

      REAL        RDGAS  
      PARAMETER ( RDGAS = 1.0E3 * RGASUNIV / MWAIR ) 

      REAL        MWWAT 
      PARAMETER ( MWWAT = 18.0153 )

      REAL        RWVAP 
      PARAMETER ( RWVAP = 1.0E3 * RGASUNIV / MWWAT ) 





      REAL        CPD 
      PARAMETER ( CPD = 7.0 * RDGAS / 2.0 )          

      REAL        CVD 
      PARAMETER ( CVD = 5.0 * RDGAS / 2.0 )          

      REAL        CPWVAP 
      PARAMETER ( CPWVAP = 4.0 * RWVAP )             

      REAL        CVWVAP 
      PARAMETER ( CVWVAP = 3.0 * RWVAP )             

      REAL        VP0 
      PARAMETER ( VP0 = 611.29 )



      REAL        LV0 
      PARAMETER ( LV0 = 2.501E6 )

      REAL        DLVDT 
                        
      PARAMETER ( DLVDT = 2370.0 ) 

      REAL        LF0 
      PARAMETER ( LF0 = 3.34E5 )





      INTEGER, PARAMETER :: NGAS = 12  
      INTEGER, PARAMETER :: NAER = 36  
      INTEGER, PARAMETER :: NLIQS = 41 



      INTEGER, PARAMETER :: LSO2    =  1  
      INTEGER, PARAMETER :: LHNO3   =  2  
      INTEGER, PARAMETER :: LN2O5   =  3  
      INTEGER, PARAMETER :: LCO2    =  4  
      INTEGER, PARAMETER :: LNH3    =  5  
      INTEGER, PARAMETER :: LH2O2   =  6  
      INTEGER, PARAMETER :: LO3     =  7  
      INTEGER, PARAMETER :: LFOA    =  8  
      INTEGER, PARAMETER :: LMHP    =  9  
      INTEGER, PARAMETER :: LPAA    = 10  
      INTEGER, PARAMETER :: LH2SO4  = 11  
      INTEGER, PARAMETER :: LHCL    = 12  



      INTEGER, PARAMETER :: LSO4AKN  =  1  
      INTEGER, PARAMETER :: LSO4ACC  =  2  
      INTEGER, PARAMETER :: LSO4COR  =  3  
      INTEGER, PARAMETER :: LNH4AKN  =  4  
      INTEGER, PARAMETER :: LNH4ACC  =  5  
      INTEGER, PARAMETER :: LNO3AKN  =  6  
      INTEGER, PARAMETER :: LNO3ACC  =  7  
      INTEGER, PARAMETER :: LNO3COR  =  8  
      INTEGER, PARAMETER :: LORGAAKN =  9  
      INTEGER, PARAMETER :: LORGAACC = 10  
      INTEGER, PARAMETER :: LORGPAKN = 11  
      INTEGER, PARAMETER :: LORGPACC = 12  
      INTEGER, PARAMETER :: LORGBAKN = 13  
      INTEGER, PARAMETER :: LORGBACC = 14  
      INTEGER, PARAMETER :: LECAKN   = 15  
      INTEGER, PARAMETER :: LECACC   = 16  
      INTEGER, PARAMETER :: LPRIAKN  = 17  
      INTEGER, PARAMETER :: LPRIACC  = 18  
      INTEGER, PARAMETER :: LPRICOR  = 19  
      INTEGER, PARAMETER :: LNAAKN   = 20  
      INTEGER, PARAMETER :: LNAACC   = 21  
      INTEGER, PARAMETER :: LNACOR   = 22  
      INTEGER, PARAMETER :: LCLAKN   = 23  
      INTEGER, PARAMETER :: LCLACC   = 24  
      INTEGER, PARAMETER :: LCLCOR   = 25  
      INTEGER, PARAMETER :: LNUMAKN  = 26  
      INTEGER, PARAMETER :: LNUMACC  = 27  
      INTEGER, PARAMETER :: LNUMCOR  = 28  
      INTEGER, PARAMETER :: LSRFAKN  = 29  
      INTEGER, PARAMETER :: LSRFACC  = 30  
      INTEGER, PARAMETER :: LNACL    = 31  
      INTEGER, PARAMETER :: LCACO3   = 32  
      INTEGER, PARAMETER :: LMGCO3   = 33  
      INTEGER, PARAMETER :: LA3FE    = 34  
      INTEGER, PARAMETER :: LB2MN    = 35  
      INTEGER, PARAMETER :: LK       = 36  



      INTEGER, PARAMETER :: LACL        =  1  
      INTEGER, PARAMETER :: LNH4L       =  2  
      INTEGER, PARAMETER :: LCAL        =  3  
      INTEGER, PARAMETER :: LNAACCL     =  4  
      INTEGER, PARAMETER :: LOHL        =  5  
      INTEGER, PARAMETER :: LSO4ACCL    =  6  
      INTEGER, PARAMETER :: LHSO4ACCL   =  7  
      INTEGER, PARAMETER :: LSO3L       =  8  
      INTEGER, PARAMETER :: LHSO3L      =  9  
      INTEGER, PARAMETER :: LSO2L       = 10  
      INTEGER, PARAMETER :: LCO3L       = 11  
      INTEGER, PARAMETER :: LHCO3L      = 12  
      INTEGER, PARAMETER :: LCO2L       = 13  
      INTEGER, PARAMETER :: LNO3ACCL    = 14  
      INTEGER, PARAMETER :: LNH3L       = 15  
      INTEGER, PARAMETER :: LCLACCL     = 16  
      INTEGER, PARAMETER :: LH2O2L      = 17  
      INTEGER, PARAMETER :: LO3L        = 18  
      INTEGER, PARAMETER :: LFEL        = 19  
      INTEGER, PARAMETER :: LMNL        = 20  
      INTEGER, PARAMETER :: LAL         = 21  
      INTEGER, PARAMETER :: LFOAL       = 22  
      INTEGER, PARAMETER :: LHCO2L      = 23  
      INTEGER, PARAMETER :: LMHPL       = 24  
      INTEGER, PARAMETER :: LPAAL       = 25  
      INTEGER, PARAMETER :: LHCLL       = 26  
      INTEGER, PARAMETER :: LPRIML      = 27  
      INTEGER, PARAMETER :: LMGL        = 28  
      INTEGER, PARAMETER :: LKL         = 29  
      INTEGER, PARAMETER :: LBL         = 30  
      INTEGER, PARAMETER :: LHNO3L      = 31  
      INTEGER, PARAMETER :: LPRIMCORL   = 32  
      INTEGER, PARAMETER :: LNUMCORL    = 33  
      INTEGER, PARAMETER :: LTS6CORL    = 34  
      INTEGER, PARAMETER :: LNACORL     = 35  
      INTEGER, PARAMETER :: LCLCORL     = 36  
      INTEGER, PARAMETER :: LNO3CORL    = 37  
      INTEGER, PARAMETER :: LORGAL      = 38  
      INTEGER, PARAMETER :: LORGPL      = 39  
      INTEGER, PARAMETER :: LORGBL      = 40  
      INTEGER, PARAMETER :: LECL        = 41  




      CHARACTER*16, SAVE :: SGRGAS  ( NGAS ) 
      CHARACTER*16, SAVE :: BUNTSGAS( NGAS ) 

      REAL, SAVE :: BGNDGAS( NGAS ) 

      DATA SGRGAS( LSO2   ), BGNDGAS( LSO2   ) /'SO2             ',   0.0 /
      DATA SGRGAS( LHNO3  ), BGNDGAS( LHNO3  ) /'HNO3            ',   0.0 /
      DATA SGRGAS( LN2O5  ), BGNDGAS( LN2O5  ) /'N2O5            ',   0.0 /
      DATA SGRGAS( LCO2   ), BGNDGAS( LCO2   ) /'CO2             ', 340.0 /
      DATA SGRGAS( LNH3   ), BGNDGAS( LNH3   ) /'NH3             ',   0.0 /
      DATA SGRGAS( LH2O2  ), BGNDGAS( LH2O2  ) /'H2O2            ',   0.0 /
      DATA SGRGAS( LO3    ), BGNDGAS( LO3    ) /'O3              ',   0.0 /
      DATA SGRGAS( LFOA   ), BGNDGAS( LFOA   ) /'FOA             ',   0.0 /
      DATA SGRGAS( LMHP   ), BGNDGAS( LMHP   ) /'MHP             ',   0.0 /
      DATA SGRGAS( LPAA   ), BGNDGAS( LPAA   ) /'PAA             ',   0.0 /
      DATA SGRGAS( LH2SO4 ), BGNDGAS( LH2SO4 ) /'H2SO4           ',   0.0 /
      DATA SGRGAS( LHCL   ), BGNDGAS( LHCL   ) /'HCL             ',   0.0 /

      DATA BUNTSGAS( LSO2   ) / 'ppm' /
      DATA BUNTSGAS( LHNO3  ) / 'ppm' /
      DATA BUNTSGAS( LN2O5  ) / 'ppm' /
      DATA BUNTSGAS( LCO2   ) / 'ppm' /
      DATA BUNTSGAS( LNH3   ) / 'ppm' /
      DATA BUNTSGAS( LH2O2  ) / 'ppm' /
      DATA BUNTSGAS( LO3    ) / 'ppm' /
      DATA BUNTSGAS( LFOA   ) / 'ppm' /
      DATA BUNTSGAS( LMHP   ) / 'ppm' /
      DATA BUNTSGAS( LPAA   ) / 'ppm' /
      DATA BUNTSGAS( LH2SO4 ) / 'ppm' /
      DATA BUNTSGAS( LHCL   ) / 'ppm' /




      CHARACTER*16, SAVE :: SGRAER  ( NAER ) 
      CHARACTER*16, SAVE :: BUNTSAER( NAER ) 

      REAL, SAVE :: SGRAERMW( NAER ) 
      REAL, SAVE :: BGNDAER ( NAER ) 

      DATA SGRAER( LSO4AKN  ), SGRAERMW( LSO4AKN  ) / 'SO4_AITKEN      ' ,  96.0 /
      DATA SGRAER( LSO4ACC  ), SGRAERMW( LSO4ACC  ) / 'SO4_ACCUM       ' ,  96.0 /
      DATA SGRAER( LSO4COR  ), SGRAERMW( LSO4COR  ) / 'SO4_COARSE      ' ,  96.0 /
      DATA SGRAER( LNH4AKN  ), SGRAERMW( LNH4AKN  ) / 'NH4_AITKEN      ' ,  18.0 /
      DATA SGRAER( LNH4ACC  ), SGRAERMW( LNH4ACC  ) / 'NH4_ACCUM       ' ,  18.0 /
      DATA SGRAER( LNO3AKN  ), SGRAERMW( LNO3AKN  ) / 'NO3_AITKEN      ' ,  62.0 /
      DATA SGRAER( LNO3ACC  ), SGRAERMW( LNO3ACC  ) / 'NO3_ACCUM       ' ,  62.0 /
      DATA SGRAER( LNO3COR  ), SGRAERMW( LNO3COR  ) / 'NO3_COARSE      ' ,  62.0 /
      DATA SGRAER( LORGAAKN ), SGRAERMW( LORGAAKN ) / 'ORGA_AITKEN     ' , 150.0 /
      DATA SGRAER( LORGAACC ), SGRAERMW( LORGAACC ) / 'ORGA_ACCUM      ' , 150.0 /
      DATA SGRAER( LORGPAKN ), SGRAERMW( LORGPAKN ) / 'ORGP_AITKEN     ' , 220.0 /
      DATA SGRAER( LORGPACC ), SGRAERMW( LORGPACC ) / 'ORGP_ACCUM      ' , 220.0 /
      DATA SGRAER( LORGBAKN ), SGRAERMW( LORGBAKN ) / 'ORGB_AITKEN     ' , 177.0 /
      DATA SGRAER( LORGBACC ), SGRAERMW( LORGBACC ) / 'ORGB_ACCUM      ' , 177.0 /
      DATA SGRAER( LECAKN   ), SGRAERMW( LECAKN   ) / 'EC_AITKEN       ' ,  12.0 /
      DATA SGRAER( LECACC   ), SGRAERMW( LECACC   ) / 'EC_ACCUM        ' ,  12.0 /
      DATA SGRAER( LPRIAKN  ), SGRAERMW( LPRIAKN  ) / 'PRI_AITKEN      ' , 200.0 /
      DATA SGRAER( LPRIACC  ), SGRAERMW( LPRIACC  ) / 'PRI_ACCUM       ' , 200.0 /
      DATA SGRAER( LPRICOR  ), SGRAERMW( LPRICOR  ) / 'PRI_COARSE      ' , 100.0 /
      DATA SGRAER( LNAAKN   ), SGRAERMW( LNAAKN   ) / 'NA_AITKEN       ' ,  23.0 /
      DATA SGRAER( LNAACC   ), SGRAERMW( LNAACC   ) / 'NA_ACCUM        ' ,  23.0 /
      DATA SGRAER( LNACOR   ), SGRAERMW( LNACOR   ) / 'NA_COARSE       ' ,  23.0 /
      DATA SGRAER( LCLAKN   ), SGRAERMW( LCLAKN   ) / 'CL_AITKEN       ' ,  35.5 /
      DATA SGRAER( LCLACC   ), SGRAERMW( LCLACC   ) / 'CL_ACCUM        ' ,  35.5 /
      DATA SGRAER( LCLCOR   ), SGRAERMW( LCLCOR   ) / 'CL_COARSE       ' ,  35.5 /
      DATA SGRAER( LNUMAKN  ), SGRAERMW( LNUMAKN  ) / 'NUM_AITKEN      ' ,   1.0 /
      DATA SGRAER( LNUMACC  ), SGRAERMW( LNUMACC  ) / 'NUM_ACCUM       ' ,   1.0 /
      DATA SGRAER( LNUMCOR  ), SGRAERMW( LNUMCOR  ) / 'NUM_COARSE      ' ,   1.0 /
      DATA SGRAER( LSRFAKN  ), SGRAERMW( LSRFAKN  ) / 'SRF_AITKEN      ' ,   1.0 /
      DATA SGRAER( LSRFACC  ), SGRAERMW( LSRFACC  ) / 'SRF_ACCUM       ' ,   1.0 /
      DATA SGRAER( LNACL    ), SGRAERMW( LNACL    ) / 'NACL            ' ,  58.4 /  
      DATA SGRAER( LCACO3   ), SGRAERMW( LCACO3   ) / 'CACO3           ' , 100.1 /
      DATA SGRAER( LMGCO3   ), SGRAERMW( LMGCO3   ) / 'MGCO3           ' ,  84.3 /
      DATA SGRAER( LA3FE    ), SGRAERMW( LA3FE    ) / 'A3FE            ' ,  55.8 /
      DATA SGRAER( LB2MN    ), SGRAERMW( LB2MN    ) / 'B2MN            ' ,  54.9 /
      DATA SGRAER( LK       ), SGRAERMW( LK       ) / 'K               ' ,  39.1 /

      DATA BGNDAER( LSO4AKN  ), BUNTSAER( LSO4AKN  ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LSO4ACC  ), BUNTSAER( LSO4ACC  ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LSO4COR  ), BUNTSAER( LSO4COR  ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LNH4AKN  ), BUNTSAER( LNH4AKN  ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LNH4ACC  ), BUNTSAER( LNH4ACC  ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LNO3AKN  ), BUNTSAER( LNO3AKN  ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LNO3ACC  ), BUNTSAER( LNO3ACC  ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LNO3COR  ), BUNTSAER( LNO3COR  ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LORGAAKN ), BUNTSAER( LORGAAKN ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LORGAACC ), BUNTSAER( LORGAACC ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LORGPAKN ), BUNTSAER( LORGPAKN ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LORGPACC ), BUNTSAER( LORGPACC ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LORGBAKN ), BUNTSAER( LORGBAKN ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LORGBACC ), BUNTSAER( LORGBACC ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LECAKN   ), BUNTSAER( LECAKN   ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LECACC   ), BUNTSAER( LECACC   ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LPRIAKN  ), BUNTSAER( LPRIAKN  ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LPRIACC  ), BUNTSAER( LPRIACC  ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LPRICOR  ), BUNTSAER( LPRICOR  ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LNAAKN   ), BUNTSAER( LNAAKN   ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LNAACC   ), BUNTSAER( LNAACC   ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LNACOR   ), BUNTSAER( LNACOR   ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LCLAKN   ), BUNTSAER( LCLAKN   ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LCLACC   ), BUNTSAER( LCLACC   ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LCLCOR   ), BUNTSAER( LCLCOR   ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LNUMAKN  ), BUNTSAER( LNUMAKN  ) /  0.0,   ' #/m3' /
      DATA BGNDAER( LNUMACC  ), BUNTSAER( LNUMACC  ) /  0.0,   ' #/m3' /
      DATA BGNDAER( LNUMCOR  ), BUNTSAER( LNUMCOR  ) /  0.0,   ' #/m3' /
      DATA BGNDAER( LSRFAKN  ), BUNTSAER( LSRFAKN  ) /  0.0,   'm2/m3' /
      DATA BGNDAER( LSRFACC  ), BUNTSAER( LSRFACC  ) /  0.0,   'm2/m3' /
      DATA BGNDAER( LNACL    ), BUNTSAER( LNACL    ) /  0.0,   'ug/m3' /  
      DATA BGNDAER( LCACO3   ), BUNTSAER( LCACO3   ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LMGCO3   ), BUNTSAER( LMGCO3   ) /  0.0,   'ug/m3' /
      DATA BGNDAER( LA3FE    ), BUNTSAER( LA3FE    ) /  0.010, 'ug/m3' /
      DATA BGNDAER( LB2MN    ), BUNTSAER( LB2MN    ) /  0.005, 'ug/m3' /
      DATA BGNDAER( LK       ), BUNTSAER( LK       ) /  0.0,   'ug/m3' /

      CHARACTER*120 XMSG           
      DATA          XMSG / ' ' /



      INTEGER      NUMOX           
      PARAMETER  ( NUMOX =  5 )

      REAL         H2ODENS         
      PARAMETER  ( H2ODENS = 1000.0 )  

      REAL         ONETHIRD       
      PARAMETER  ( ONETHIRD = 1.0 / 3.0 )

      REAL         TWOTHIRDS       
      PARAMETER  ( TWOTHIRDS = 2.0 / 3.0 )

      REAL         CONCMIN         
      PARAMETER  ( CONCMIN = 1.0E-30 )

      REAL         SEC2HR          
      PARAMETER  ( SEC2HR = 1.0 / 3600.0 )






      REAL         AIRM                       
      REAL         ALFA0                      
      REAL         ALFA2                      
      REAL         ALFA3                      
      REAL         HPWDEP                     
      REAL         PRCRATE                    
      REAL         PRES_PA                    
      REAL         TAUCLD                     
      REAL         TEMP                       
      REAL         WCAVG                      
      REAL         WTAVG                      
      
      REAL,  INTENT(INOUT) :: GAS ( NGAS )    
      REAL,  INTENT(INOUT) :: AEROSOL( NAER ) 
      REAL,  INTENT(OUT)   :: LIQUID( NLIQS ) 
      REAL,  INTENT(OUT)   :: GASWDEP( NGAS ) 
      REAL,  INTENT(OUT)   :: AERWDEP( NAER ) 



      LOGICAL, SAVE :: FIRSTIME = .TRUE. 

      CHARACTER*6  PNAME           
      DATA         PNAME / 'AQCHEM' /
      SAVE         PNAME
      CHARACTER( 16 ), SAVE :: AE_VRSN 

      INTEGER      I20C            
      INTEGER      I30C            
      INTEGER      ITERAT          
      INTEGER      I7777C          
      INTEGER      ICNTAQ          
      INTEGER      LIQ             
      INTEGER      IOX             

      REAL         DEPSUM
      REAL         BETASO4
      REAL         A               
      REAL         AC              
      REAL         ACT1            
      REAL         ACT2            
      REAL         ACTB            
      REAL         AE              
      REAL         B               
      REAL         PRES_ATM        
      REAL         BB              
      REAL         CA              
      REAL         CAA             
      REAL         CL              
      REAL         CLACC           
      REAL         CLACCA          
      REAL         CLAKNA          
      REAL         CLCOR           
      REAL         CLCORA          
      REAL         CO2H            
      REAL         CO21            
      REAL         CO22            
      REAL         CO212           
      REAL         CO212H          
      REAL         CO21H           
      REAL         CO2L            
      REAL         CO3             
      REAL         CO3A            
      REAL         CTHK1           
      REAL         DTRMV           
      REAL         DTS6            
      REAL         EBETASO4T       
      REAL         EALFA0T         
      REAL         EALFA2T         
      REAL         EALFA3T         
      REAL         EC              
      REAL         ECACCA          
      REAL         ECAKNA          
      REAL         FA              
      REAL         FB              
      REAL         FE              
      REAL         FEA             
      REAL         FNH3            
      REAL         FNH4ACC         
      REAL         FHNO3           
      REAL         FNO3ACC         
      REAL         FRACLIQ         
      REAL         FOA1            
      REAL         FOAH            
      REAL         FOA1H           
      REAL         FOAL            
      REAL         FTST            
      REAL         GM              
      REAL         GM1             
      REAL         GM1LOG          
      REAL         GM2             
      REAL         GM2LOG          
      REAL         HA              
      REAL         HB              
      REAL         H2OW            
      REAL         H2O2H           
      REAL         H2O2L           
      REAL         HCLH            
      REAL         HCL1            
      REAL         HCL1H           
      REAL         HCLL            
      REAL         HCO2            
      REAL         HCO3            
      REAL         HNO3H           
      REAL         HNO31           
      REAL         HNO31H          
      REAL         HNO3L           
      REAL         HSO3            
      REAL         HSO4            
      REAL         HSO4ACC         
      REAL         HSO4COR         
      REAL         HTST            
      REAL         K               
      REAL         KA              
      REAL         LGTEMP          
      REAL         M3NEW           
      REAL         M3OLD           
      REAL         MG              
      REAL         MGA             
      REAL         MHPH            
      REAL         MHPL            
      REAL         MN              
      REAL         MNA             
      REAL         NA              
      REAL         NAACC           
      REAL         NAACCA          
      REAL         NAAKNA          
      REAL         NACOR           
      REAL         NACORA          
      REAL         NH31            
      REAL         NH3H            
      REAL         NH3DH20         
      REAL         NH31HDH         
      REAL         NH3L            
      REAL         NH4             
      REAL         NH4AKNA         
      REAL         NH4ACCA         
      REAL         NITAER          
      REAL         NO3             
      REAL         NO3ACC          
      REAL         NO3ACCA         
      REAL         NO3AKNA         
      REAL         NO3CORA         
      REAL         NO3COR          
      REAL         NUMCOR          
      REAL         NUMCORA         
      REAL         O3H             
      REAL         O3L             
      REAL         OH              
      REAL         ORGA            
      REAL         ORGAACCA        
      REAL         ORGAAKNA        
      REAL         ORGP            
      REAL         ORGPACCA        
      REAL         ORGPAKNA        
      REAL         ORGB            
      REAL         ORGBACCA        
      REAL         ORGBAKNA        
      REAL         PAAH            
      REAL         PAAL            
      REAL         PCO20           
      REAL         PCO2F           
      REAL         PFOA0           
      REAL         PFOAF           
      REAL         PH2O20          
      REAL         PH2O2F          
      REAL         PHCL0           
      REAL         PHCLF           
      REAL         PHNO30          
      REAL         PHNO3F          
      REAL         PMHP0           
      REAL         PMHPF           
      REAL         PNH30           
      REAL         PNH3F           
      REAL         PO30            
      REAL         PO3F            
      REAL         PPAA0           
      REAL         PPAAF           
      REAL         PRIM            
      REAL         PRIMCOR         
      REAL         PRIACCA         
      REAL         PRIAKNA         
      REAL         PRICORA         
      REAL         PSO20           
      REAL         PSO2F           
      REAL         RATE            
      REAL         RECIPA1         
      REAL         RECIPA2         
      REAL         RECIPAP1        
      REAL         RH2O2           
      REAL         RMHP            
      REAL         RPAA            
      REAL         RT              
      REAL         SCVEFF          
      DATA         SCVEFF / 100.0 / 
      SAVE         SCVEFF
      REAL         SIV             
      REAL         SK6             
      REAL         SK6TS6          
      REAL         SO21            
      REAL         SO22            
      REAL         SO2H            
      REAL         SO212           
      REAL         SO212H          
      REAL         SO21H           
      REAL         SO2L            
      REAL         SO3             
      REAL         SO4             
      REAL         SO4ACC          
      REAL         SO4COR          
      REAL         STION           
      REAL         TAC             
      REAL         TEMP1           
      REAL         TIMEW           
      REAL         TOTOX           
      REAL         TOTAMM          
      REAL         TOTNIT          
      REAL         TS6             
      REAL         TS6AKNA         
      REAL         TS6ACC          
      REAL         TS6ACCA         
      REAL         TS6COR          
      REAL         TS6CORA         
      REAL         TSIV            
      REAL         TST             
      REAL         TWASH           
      REAL         WETFAC          
      REAL         XC1             
      REAL         XC2             
      REAL         XL              
      REAL         ONE_OVER_XL     
      REAL         PRES_ATM_OVER_XL 
      REAL         XLCO2           
      REAL         XLH2O2          
      REAL         XLHCL           
      REAL         XLHNO3          
      REAL         XLMHP           
      REAL         XLNH3           
      REAL         XLO3            
      REAL         XLPAA           
      REAL         XLSO2           



      REAL         WETDEP( NLIQS ) 
      REAL         DSIVDT( 0:NUMOX ) 
      REAL         DS4   ( 0:NUMOX ) 
      REAL         DTW   ( 0:NUMOX ) 

      REAL         ONE_OVER_TEMP     









      ONE_OVER_TEMP = 1.0 / TEMP



      IF ( TEMP .LE. 0.0 .OR. AIRM .LE. 0.0 .OR. PRES_PA .LE. 0.0 ) THEN
         XMSG = 'MET DATA ERROR, EXITING ROUTINE.'

        write(0,*) ''
        write(0,*) PNAME,' : ',XMSG
        write(0,*) ''
        write(0,*) 'TEMP :'
        write(0,*) TEMP
        write(0,*) 'PRES_PA :'
        write(0,*) PRES_PA
        write(0,*) 'TAUCLD :'
        write(0,*) TAUCLD
        write(0,*) 'PRCRATE :'
        write(0,*) PRCRATE
        write(0,*) 'WCAVG :'
        write(0,*) WCAVG
        write(0,*) 'WTAVG :'
        write(0,*) WTAVG
        write(0,*) 'AIRM :'
        write(0,*) AIRM
        write(0,*) 'ALFA0 :'
        write(0,*) ALFA0
        write(0,*) 'ALFA2 :'
        write(0,*) ALFA2
        write(0,*) 'ALFA3 :'
        write(0,*) ALFA3
        write(0,*) 'GAS :'
        write(0,*) GAS
        write(0,*) 'AEROSOL :'
        write(0,*) AEROSOL
        write(0,*) 'GASWDEP :'
        write(0,*) GASWDEP
        write(0,*) 'AERWDEP :'
        write(0,*) AERWDEP
        write(0,*) 'HPWDEP :'
        write(0,*) HPWDEP
        write(0,*) ''
        return
      END IF



      ICNTAQ = 0
      ITERAT = 0
      RT = ( MOLVOL / STDTEMP ) * TEMP             
      PRES_ATM = PRES_PA /  STDATMPA               
      CTHK1 = AIRM * RT / ( PRES_ATM * 1000.0 )    
      XL   = WCAVG * RT / H2ODENS     
      ONE_OVER_XL = 1.0 / XL
      PRES_ATM_OVER_XL = PRES_ATM / XL
      TST  = 0.999
      GM   = SCVEFF / 100.0
      ACT1 = 1.0
      ACT2 = 1.0
      GM2  = 1.0
      TIMEW = 0.0
      RECIPAP1 = 1.0 / PRES_ATM
      XC1  = 1.0 / ( WCAVG * CTHK1 )
      XC2  = RT / ( 1000.0 * CTHK1 )
      FRACLIQ = WCAVG / WTAVG
      TWASH = WTAVG * 1000.0 * CTHK1 * 3600.0 &
            / ( H2ODENS * MAX( 1.0E-20, PRCRATE ) )




      SO2H  = HLCONST( 'SO2             ', TEMP, .FALSE., 0.0 )
      CO2H  = HLCONST( 'CO2             ', TEMP, .FALSE., 0.0 )
      NH3H  = HLCONST( 'NH3             ', TEMP, .FALSE., 0.0 )
      H2O2H = HLCONST( 'H2O2            ', TEMP, .FALSE., 0.0 )
      O3H   = HLCONST( 'O3              ', TEMP, .FALSE., 0.0 )
      HCLH  = HLCONST( 'HCL             ', TEMP, .FALSE., 0.0 )
      HNO3H = HLCONST( 'HNO3            ', TEMP, .FALSE., 0.0 )
      MHPH  = HLCONST( 'METHYLHYDROPEROX', TEMP, .FALSE., 0.0 )
      PAAH  = HLCONST( 'PEROXYACETIC_ACI', TEMP, .FALSE., 0.0 )
      FOAH  = HLCONST( 'FORMIC_ACID     ', TEMP, .FALSE., 0.0 )

      TEMP1 = ONE_OVER_TEMP - 1.0 / 298.0



      FOA1  = 1.80E-04 * EXP( -2.00E+01 * TEMP1 )      
      SK6   = 1.02E-02 * EXP(  2.72E+03 * TEMP1 )      
      SO21  = 1.30E-02 * EXP(  1.96E+03 * TEMP1 )      
      SO22  = 6.60E-08 * EXP(  1.50E+03 * TEMP1 )      
      CO21  = 4.30E-07 * EXP( -1.00E+03 * TEMP1 )      
      CO22  = 4.68E-11 * EXP( -1.76E+03 * TEMP1 )      
      H2OW  = 1.00E-14 * EXP( -6.71E+03 * TEMP1 )      
      NH31  = 1.70E-05 * EXP( -4.50E+02 * TEMP1 )      
      HCL1  = 1.74E+06 * EXP(  6.90E+03 * TEMP1 )      
      HNO31 = 1.54E+01 * EXP(  8.70E+03 * TEMP1 )      






       RH2O2 = 7.45E+07 * EXP( -15.96E0 * ( ( 298.0E0 / TEMP )  - 1.0E0 ) )







      RMHP = 1.90E+07 * EXP( -12.75E0 * ( ( 298.0E0 / TEMP )  - 1.0E0 ) )
      RPAA = 3.67E+07 * EXP( -13.42E0 * ( ( 298.0E0 / TEMP )  - 1.0E0 ) )



      DO LIQ = 1, NLIQS
        WETDEP( LIQ ) = 0.0
      END DO

      DO IOX = 0, NUMOX
        DSIVDT( IOX ) = 0.0
        DTW   ( IOX ) = 0.0
        DS4   ( IOX ) = 0.0
      END DO




      M3OLD = ( AEROSOL( LSO4ACC  ) * SGRAERMW( LSO4ACC  ) / 1.8e6 &
            +   AEROSOL( LNH4ACC  ) * SGRAERMW( LNH4ACC  ) / 1.8e6 &
            +   AEROSOL( LNO3ACC  ) * SGRAERMW( LNO3ACC  ) / 1.8e6 &
            +   AEROSOL( LORGPACC ) * SGRAERMW( LORGPACC ) / 2.0e6 &
            +   AEROSOL( LECACC   ) * SGRAERMW( LECACC   ) / 2.2e6 &
            +   AEROSOL( LPRIACC  ) * SGRAERMW( LPRIACC  ) / 2.2e6 &
            +   AEROSOL( LNAACC   ) * SGRAERMW( LNAACC   ) / 2.2e6 &
            +   AEROSOL( LCLACC   ) * SGRAERMW( LCLACC   ) / 2.2e6 )




      TOTNIT = GAS( LHNO3 ) + AEROSOL( LNO3ACC )
      IF ( TOTNIT .GT. 0.0 ) THEN
        FHNO3   = GAS( LHNO3 ) / TOTNIT
        FNO3ACC = AEROSOL( LNO3ACC ) / TOTNIT
      ELSE
        FHNO3   = 1.0
        FNO3ACC = 0.0
      END IF

      TOTAMM = GAS( LNH3 ) + AEROSOL( LNH4ACC )
      IF ( TOTAMM .GT. 0.0 ) THEN
        FNH3    = GAS( LNH3 ) / TOTAMM
        FNH4ACC = AEROSOL( LNH4ACC ) / TOTAMM
      ELSE
        FNH3    = 1.0
        FNH4ACC = 0.0
      END IF





      TS6ACCA = ( AEROSOL( LSO4ACC ) &
              +   GAS    ( LH2SO4  ) ) * PRES_ATM_OVER_XL
      NO3ACCA =   AEROSOL( LNO3ACC )   * PRES_ATM_OVER_XL
      NH4ACCA =   AEROSOL( LNH4ACC )   * PRES_ATM_OVER_XL
      ORGAACCA =  AEROSOL( LORGAACC )  * PRES_ATM_OVER_XL
      ORGPACCA =  AEROSOL( LORGPACC )  * PRES_ATM_OVER_XL
      ORGBACCA =  AEROSOL( LORGBACC )  * PRES_ATM_OVER_XL
      ECACCA  =   AEROSOL( LECACC  )   * PRES_ATM_OVER_XL
      PRIACCA =   AEROSOL( LPRIACC )   * PRES_ATM_OVER_XL
      NAACCA  =   AEROSOL( LNAACC  )   * PRES_ATM_OVER_XL
      CLACCA  =   AEROSOL( LCLACC  )   * PRES_ATM_OVER_XL





      TS6CORA =   AEROSOL( LSO4COR )   * PRES_ATM_OVER_XL
      NO3CORA =   AEROSOL( LNO3COR )   * PRES_ATM_OVER_XL

      IF ( AE_VRSN .EQ. 'AE3' ) THEN
        CLCORA  = AEROSOL( LNACL   )   * PRES_ATM_OVER_XL
        NACORA  = AEROSOL( LNACL   )   * PRES_ATM_OVER_XL
      ELSE
        CLCORA  = AEROSOL( LCLCOR  )   * PRES_ATM_OVER_XL
        NACORA  = AEROSOL( LNACOR  )   * PRES_ATM_OVER_XL
      END IF

      KA      =   AEROSOL( LK      )   * PRES_ATM_OVER_XL
      CAA     =   AEROSOL( LCACO3  )   * PRES_ATM_OVER_XL
      MGA     =   AEROSOL( LMGCO3  )   * PRES_ATM_OVER_XL
      FEA     =   AEROSOL( LA3FE   )   * PRES_ATM_OVER_XL
      MNA     =   AEROSOL( LB2MN   )   * PRES_ATM_OVER_XL
      CO3A    = ( AEROSOL( LCACO3  ) &
              +   AEROSOL( LMGCO3  ) ) * PRES_ATM_OVER_XL
      PRICORA =   AEROSOL( LPRICOR )   * PRES_ATM_OVER_XL
      NUMCORA =   AEROSOL( LNUMCOR )   * PRES_ATM_OVER_XL



      XLH2O2  = H2O2H * XL
      XLO3    = O3H   * XL
      XLMHP   = MHPH  * XL
      XLPAA   = PAAH  * XL
      XLSO2   = SO2H  * XL
      XLNH3   = NH3H  * XL
      XLHCL   = HCLH  * XL
      XLHNO3  = HNO3H * XL
      XLCO2   = CO2H  * XL

      SO212   = SO21  * SO22
      SO21H   = SO21  * SO2H
      SO212H  = SO212 * SO2H
      CO212   = CO21  * CO22
      CO21H   = CO21  * CO2H
      CO212H  = CO22  * CO21H
      NH3DH20 = NH31  / H2OW
      NH31HDH = NH3H  * NH3DH20
      FOA1H   = FOA1  * FOAH
      HCL1H   = HCL1  * HCLH
      HNO31H  = HNO31 * HNO3H



      I20C = 0
20    CONTINUE

      I20C = I20C + 1
      IF ( I20C .GE. 10000 ) THEN
        XMSG = 'EXCESSIVE LOOPING AT I20C, EXITING ROUTINE.'

        write(0,*) ''
        write(0,*) PNAME,' : ',XMSG
        write(0,*) ''
        write(0,*) 'TEMP :'
        write(0,*) TEMP
        write(0,*) 'PRES_PA :'
        write(0,*) PRES_PA
        write(0,*) 'TAUCLD :'
        write(0,*) TAUCLD
        write(0,*) 'PRCRATE :'
        write(0,*) PRCRATE
        write(0,*) 'WCAVG :'
        write(0,*) WCAVG
        write(0,*) 'WTAVG :'
        write(0,*) WTAVG
        write(0,*) 'AIRM :'
        write(0,*) AIRM
        write(0,*) 'ALFA0 :'
        write(0,*) ALFA0
        write(0,*) 'ALFA2 :'
        write(0,*) ALFA2
        write(0,*) 'ALFA3 :'
        write(0,*) ALFA3
        write(0,*) 'GAS :'
        write(0,*) GAS
        write(0,*) 'AEROSOL :'
        write(0,*) AEROSOL
        write(0,*) 'GASWDEP :'
        write(0,*) GASWDEP
        write(0,*) 'AERWDEP :'
        write(0,*) AERWDEP
        write(0,*) 'HPWDEP :'
        write(0,*) HPWDEP
        write(0,*) ''
        return
      END IF



      NO3AKNA = AEROSOL( LNO3AKN ) * PRES_ATM_OVER_XL &
              * ( 1.0 - EXP( -ALFA3 * TIMEW ) )
      NH4AKNA = AEROSOL( LNH4AKN ) * PRES_ATM_OVER_XL &
              * ( 1.0 - EXP( -ALFA3 * TIMEW ) )
      TS6AKNA = AEROSOL( LSO4AKN ) * PRES_ATM_OVER_XL &
              * ( 1.0 - EXP( -ALFA3 * TIMEW ) )
      ORGAAKNA = AEROSOL( LORGAAKN ) * PRES_ATM_OVER_XL &
               * ( 1.0 - EXP( -ALFA3 * TIMEW ) )
      ORGPAKNA = AEROSOL( LORGPAKN ) * PRES_ATM_OVER_XL &
               * ( 1.0 - EXP( -ALFA3 * TIMEW ) )
      ORGBAKNA = AEROSOL( LORGBAKN ) * PRES_ATM_OVER_XL &
               * ( 1.0 - EXP( -ALFA3 * TIMEW ) )
      ECAKNA  = AEROSOL( LECAKN  ) * PRES_ATM_OVER_XL &
              * ( 1.0 - EXP( -ALFA3 * TIMEW ) )
      PRIAKNA = AEROSOL( LPRIAKN ) * PRES_ATM_OVER_XL &
              * ( 1.0 - EXP( -ALFA3 * TIMEW ) )
      NAAKNA  = AEROSOL( LNAAKN  ) * PRES_ATM_OVER_XL &
              * ( 1.0 - EXP( -ALFA3 * TIMEW ) )
      CLAKNA  = AEROSOL( LCLAKN  ) * PRES_ATM_OVER_XL &
              * ( 1.0 - EXP( -ALFA3 * TIMEW ) )




      PSO20  = GAS( LSO2  ) * PRES_ATM &
             + DS4( 0 ) * XL &
             - ( WETDEP( LSO3L ) + WETDEP( LHSO3L ) + WETDEP( LSO2L ) ) * XC2
      PNH30  = GAS( LNH3  ) * PRES_ATM &
             + ( NH4ACCA + NH4AKNA ) * XL &
             - ( WETDEP( LNH4L ) + WETDEP( LNH3L ) ) * XC2
      PHNO30 = ( GAS( LHNO3 ) + 2.0 * GAS( LN2O5 ) ) * PRES_ATM &
             + ( NO3ACCA + NO3CORA + NO3AKNA ) * XL &
             - ( WETDEP( LNO3ACCL ) + WETDEP( LHNO3L ) + WETDEP( LNO3CORL ) ) * XC2
      PHCL0  = GAS(  LHCL ) * PRES_ATM &
             + ( CLACCA + CLCORA + CLAKNA ) * XL & 
             - ( WETDEP( LCLACCL ) + WETDEP( LHCLL ) + WETDEP( LCLCORL ) ) * XC2
      PH2O20 = GAS( LH2O2 ) * PRES_ATM - WETDEP( LH2O2L ) * XC2
      PO30   = GAS( LO3   ) * PRES_ATM - WETDEP( LO3L   ) * XC2
      PFOA0  = GAS( LFOA  ) * PRES_ATM &
             - ( WETDEP( LFOAL ) + WETDEP( LHCO2L ) ) * XC2
      PMHP0  = GAS( LMHP  ) * PRES_ATM - WETDEP( LMHPL  ) * XC2
      PPAA0  = GAS( LPAA  ) * PRES_ATM - WETDEP( LPAAL  ) * XC2
      PCO20  = GAS( LCO2  ) * PRES_ATM &
             + CO3A * XL &
             - ( WETDEP( LCO3L ) + WETDEP( LHCO3L ) + WETDEP( LCO2L ) ) * XC2



      PSO20  = MAX( PSO20,  0.0 )
      PNH30  = MAX( PNH30,  0.0 )
      PH2O20 = MAX( PH2O20, 0.0 )
      PO30   = MAX( PO30,   0.0 )
      PFOA0  = MAX( PFOA0,  0.0 )
      PMHP0  = MAX( PMHP0,  0.0 )
      PPAA0  = MAX( PPAA0,  0.0 )
      PCO20  = MAX( PCO20,  0.0 )
      PHCL0  = MAX( PHCL0,  0.0 )
      PHNO30 = MAX( PHNO30, 0.0 )




      TS6COR  = MAX( TS6CORA - WETDEP( LTS6CORL ) * XC1, 0.0 )
      NO3COR  = MAX( NO3CORA - WETDEP( LNO3CORL ) * XC1, 0.0 )
      NACOR   = MAX( NACORA  - WETDEP( LNACORL  ) * XC1, 0.0 )
      CLCOR   = MAX( CLCORA  - WETDEP( LCLCORL  ) * XC1, 0.0 )

      TS6     = TS6ACCA  + TS6AKNA + TS6COR &
              - ( WETDEP( LSO4ACCL ) + WETDEP( LHSO4ACCL ) ) * XC1 &
              - DS4( 0 )
      NA      = NAACCA   + NAAKNA  + NACOR &
              - WETDEP( LNAACCL ) * XC1
      CA      = CAA      -   WETDEP( LCAL ) * XC1
      MG      = MGA      -   WETDEP( LMGL ) * XC1
      K       = KA       -   WETDEP( LKL  ) * XC1
      FE      = FEA      -   WETDEP( LFEL ) * XC1
      MN      = MNA      -   WETDEP( LMNL ) * XC1
      ORGA    = ORGAACCA + ORGAAKNA - WETDEP( LORGAL ) * XC1
      ORGP    = ORGPACCA + ORGPAKNA - WETDEP( LORGPL ) * XC1
      ORGB    = ORGBACCA + ORGBAKNA - WETDEP( LORGBL ) * XC1
      EC      = ECACCA   + ECAKNA   - WETDEP( LECL   ) * XC1
      PRIM    = PRIACCA  + PRIAKNA  - WETDEP( LPRIML ) * XC1
      PRIMCOR = PRICORA - WETDEP( LPRIMCORL ) * XC1
      NUMCOR  = NUMCORA - WETDEP( LNUMCORL  ) * XC1
      A       = 3.0 * FE
      B       = 2.0 * MN



      TS6     = MAX( TS6,     0.0 )
      NA      = MAX( NA,      0.0 )
      CA      = MAX( CA,      0.0 )
      MG      = MAX( MG,      0.0 )
      K       = MAX( K,       0.0 )
      FE      = MAX( FE,      0.0 )
      MN      = MAX( MN,      0.0 )
      ORGA    = MAX( ORGA,    0.0 )
      ORGP    = MAX( ORGP,    0.0 )
      ORGB    = MAX( ORGB,    0.0 )
      EC      = MAX( EC,      0.0 )
      PRIM    = MAX( PRIM,    0.0 )
      PRIMCOR = MAX( PRIMCOR, 0.0 )
      NUMCOR  = MAX( NUMCOR,  0.0 )
      A       = MAX( A,       0.0 )
      B       = MAX( B,       0.0 )

      SK6TS6 = SK6 * TS6




      HA =  0.01
      HB = 10.0

      I7777C = 0
7777  CONTINUE

      I7777C = I7777C + 1
      IF ( I7777C .GE. 10000 ) THEN
        XMSG = 'EXCESSIVE LOOPING AT I7777C, EXITING ROUTINE.'

        write(0,*) ''
        write(0,*) PNAME,' : ',XMSG
        write(0,*) ''
        write(0,*) 'TEMP :'
        write(0,*) TEMP
        write(0,*) 'PRES_PA :'
        write(0,*) PRES_PA
        write(0,*) 'TAUCLD :'
        write(0,*) TAUCLD
        write(0,*) 'PRCRATE :'
        write(0,*) PRCRATE
        write(0,*) 'WCAVG :'
        write(0,*) WCAVG
        write(0,*) 'WTAVG :'
        write(0,*) WTAVG
        write(0,*) 'AIRM :'
        write(0,*) AIRM
        write(0,*) 'ALFA0 :'
        write(0,*) ALFA0
        write(0,*) 'ALFA2 :'
        write(0,*) ALFA2
        write(0,*) 'ALFA3 :'
        write(0,*) ALFA3
        write(0,*) 'GAS :'
        write(0,*) GAS
        write(0,*) 'AEROSOL :'
        write(0,*) AEROSOL
        write(0,*) 'GASWDEP :'
        write(0,*) GASWDEP
        write(0,*) 'AERWDEP :'
        write(0,*) AERWDEP
        write(0,*) 'HPWDEP :'
        write(0,*) HPWDEP
        write(0,*) ''
        return
      END IF

      HA = MAX( HA - 0.8, 0.1 )
      HB = MIN( HB + 0.8, 9.9 )
      AE = 10.0**( -HA )

      RECIPA1 = 1.0 / ( AE * ACT1 )
      RECIPA2 = 1.0 / ( AE * AE * ACT2 )




      PSO2F = PSO20 / ( 1.0 + XLSO2 * ( 1.0 + SO21 * RECIPA1 &
            + SO212 * RECIPA2 ) )

      PNH3F = PNH30 / ( 1.0 + XLNH3 * ( 1.0 + NH3DH20 * AE ) )

      PHCLF = PHCL0 / ( 1.0 + XLHCL *  ( 1.0 + HCL1 * RECIPA1 ) )

      PFOAF = PFOA0 / ( 1.0 + XL * ( FOAH + FOA1H * RECIPA1 ) )

      PHNO3F = PHNO30 / ( 1.0 + XLHNO3 * ( 1.0 + HNO31 * RECIPA1 ) )

      PCO2F = PCO20 / ( 1.0 + XLCO2 * ( 1.0 + CO21 * RECIPA1 &
            + CO212 * RECIPA2 ) )



      SO4  = SK6TS6 / ( AE * GM2 + SK6 )
      HSO4 = TS6 - SO4
      SO3  = SO212H  * PSO2F  * RECIPA2
      HSO3 = SO21H   * PSO2F  * RECIPA1
      CO3  = CO212H  * PCO2F  * RECIPA2
      HCO3 = CO21H   * PCO2F  * RECIPA1
      OH   = H2OW    * RECIPA1
      NH4  = NH31HDH * PNH3F  * AE
      HCO2 = FOA1H   * PFOAF  * RECIPA1
      NO3  = HNO31H  * PHNO3F * RECIPA1
      CL   = HCL1H   * PHCLF  * RECIPA1 



      FA = AE + NH4 + NA + 2.0 * ( CA + MG - CO3 - SO3 - SO4 ) &
         - OH - HCO3 - HSO3 - NO3 - HSO4 - HCO2 - CL



      I30C = 0
30    CONTINUE

      I30C = I30C + 1
      IF ( I30C .GE. 10000 ) THEN
        XMSG = 'EXCESSIVE LOOPING AT I30C, EXITING ROUTINE.'

        write(0,*) ''
        write(0,*) PNAME,' : ',XMSG
        write(0,*) ''
        write(0,*) 'TEMP :'
        write(0,*) TEMP
        write(0,*) 'PRES_PA :'
        write(0,*) PRES_PA
        write(0,*) 'TAUCLD :'
        write(0,*) TAUCLD
        write(0,*) 'PRCRATE :'
        write(0,*) PRCRATE
        write(0,*) 'WCAVG :'
        write(0,*) WCAVG
        write(0,*) 'WTAVG :'
        write(0,*) WTAVG
        write(0,*) 'AIRM :'
        write(0,*) AIRM
        write(0,*) 'ALFA0 :'
        write(0,*) ALFA0
        write(0,*) 'ALFA2 :'
        write(0,*) ALFA2
        write(0,*) 'ALFA3 :'
        write(0,*) ALFA3
        write(0,*) 'GAS :'
        write(0,*) GAS
        write(0,*) 'AEROSOL :'
        write(0,*) AEROSOL
        write(0,*) 'GASWDEP :'
        write(0,*) GASWDEP
        write(0,*) 'AERWDEP :'
        write(0,*) AERWDEP
        write(0,*) 'HPWDEP :'
        write(0,*) HPWDEP
        write(0,*) ''
        return
      END IF

      BB = ( HA + HB ) / 2.0
      AE = 10.0**( -BB )

      ICNTAQ = ICNTAQ + 1
      IF ( ICNTAQ .GE. 60000 ) THEN
        XMSG = 'Maximum AQCHEM total iterations exceeded, EXITING ROUTINE.'

        write(0,*) ''
        write(0,*) PNAME,' : ',XMSG
        write(0,*) ''
        write(0,*) 'TEMP :'
        write(0,*) TEMP
        write(0,*) 'PRES_PA :'
        write(0,*) PRES_PA
        write(0,*) 'TAUCLD :'
        write(0,*) TAUCLD
        write(0,*) 'PRCRATE :'
        write(0,*) PRCRATE
        write(0,*) 'WCAVG :'
        write(0,*) WCAVG
        write(0,*) 'WTAVG :'
        write(0,*) WTAVG
        write(0,*) 'AIRM :'
        write(0,*) AIRM
        write(0,*) 'ALFA0 :'
        write(0,*) ALFA0
        write(0,*) 'ALFA2 :'
        write(0,*) ALFA2
        write(0,*) 'ALFA3 :'
        write(0,*) ALFA3
        write(0,*) 'GAS :'
        write(0,*) GAS
        write(0,*) 'AEROSOL :'
        write(0,*) AEROSOL
        write(0,*) 'GASWDEP :'
        write(0,*) GASWDEP
        write(0,*) 'AERWDEP :'
        write(0,*) AERWDEP
        write(0,*) 'HPWDEP :'
        write(0,*) HPWDEP
        write(0,*) ''
        return
      END IF

      RECIPA1 = 1.0 / ( AE * ACT1 )
      RECIPA2 = 1.0 / ( AE * AE * ACT2 )




      PSO2F = PSO20 / ( 1.0 + XLSO2 &
      	    * ( 1.0 + SO21 * RECIPA1 + SO212 * RECIPA2 ) )

      PNH3F = PNH30 / ( 1.0 + XLNH3 * ( 1.0 + NH3DH20 * AE ) )

      PHCLF = PHCL0  / ( 1.0 + XLHCL *  ( 1.0 + HCL1 * RECIPA1 ) )

      PHNO3F = PHNO30 / ( 1.0 + XLHNO3 * ( 1.0 + HNO31 * RECIPA1 ) )

      PFOAF = PFOA0 / ( 1.0 + XL * ( FOAH + FOA1H * RECIPA1 ) )

      PCO2F = PCO20 / ( 1.0 + XLCO2 * ( 1.0 + CO21 * RECIPA1 &
            + CO212 * RECIPA2 ) )



      SO4  = SK6TS6 / ( AE * GM2 + SK6 )
      HSO4 = TS6 - SO4
      SO3  = SO212H  * PSO2F  * RECIPA2
      HSO3 = SO21H   * PSO2F  * RECIPA1
      CO3  = CO212H  * PCO2F  * RECIPA2
      HCO3 = CO21H   * PCO2F  * RECIPA1
      OH   = H2OW    * RECIPA1
      NH4  = NH31HDH * PNH3F  * AE
      HCO2 = FOA1H   * PFOAF  * RECIPA1
      NO3  = HNO31H  * PHNO3F * RECIPA1
      CL   = HCL1H   * PHCLF  * RECIPA1 



      FB = AE + NH4 + NA + 2.0 * ( CA + MG - CO3 - SO3 - SO4 ) &
           - OH - HCO3 - HSO3 - NO3 - HSO4 - HCO2 - CL



      FTST = FA * FB
      IF ( FTST .LE. 0.0 ) THEN
        HB = BB
      ELSE
        HA = BB
        FA = FB
      END IF



      HTST = HA / HB
      IF ( HTST .LE. TST ) GO TO 30





      STION = 0.5 * (AE + NH4 + OH + HCO3 + HSO3 &
            + 4.0 * (SO4 + CO3 + SO3 + CA + MG + MN) &
            + NO3 + HSO4 + 9.0 * FE + NA + K + CL + A + B + HCO2)
      GM1LOG = -0.509 * ( SQRT( STION ) &
             / ( 1.0 + SQRT( STION ) ) - 0.2 * STION )
      GM2LOG = GM1LOG * 4.0
      GM1  = 10.0**GM1LOG
      GM2  = MAX( 10.0**GM2LOG, 1.0E-30 )
      ACTB = ACT1
      ACT1 = MAX( GM1 * GM1, 1.0E-30 )
      ACT2 = MAX( GM1 * GM1 * GM2, 1.0E-30 )




      TAC = ABS( ACTB - ACT1 ) / ACTB
      IF ( TAC .GE. 1.0E-2 ) GO TO 7777




      IF ( ( HA .LT. 0.1 ) .OR. ( HA .GT. 9.9 ) ) THEN
        print *, ha
        XMSG = 'PH VALUE OUT OF RANGE, EXITING ROUTINE.'

        write(0,*) ''
        write(0,*) PNAME,' : ',XMSG
        write(0,*) ''
        write(0,*) 'TEMP :'
        write(0,*) TEMP
        write(0,*) 'PRES_PA :'
        write(0,*) PRES_PA
        write(0,*) 'TAUCLD :'
        write(0,*) TAUCLD
        write(0,*) 'PRCRATE :'
        write(0,*) PRCRATE
        write(0,*) 'WCAVG :'
        write(0,*) WCAVG
        write(0,*) 'WTAVG :'
        write(0,*) WTAVG
        write(0,*) 'AIRM :'
        write(0,*) AIRM
        write(0,*) 'ALFA0 :'
        write(0,*) ALFA0
        write(0,*) 'ALFA2 :'
        write(0,*) ALFA2
        write(0,*) 'ALFA3 :'
        write(0,*) ALFA3
        write(0,*) 'GAS :'
        write(0,*) GAS
        write(0,*) 'AEROSOL :'
        write(0,*) AEROSOL
        write(0,*) 'GASWDEP :'
        write(0,*) GASWDEP
        write(0,*) 'AERWDEP :'
        write(0,*) AERWDEP
        write(0,*) 'HPWDEP :'
        write(0,*) HPWDEP
        write(0,*) ''
        return
      END IF




      SO2L = SO2H * PSO2F
      AC = 10.0**( -BB )
      SIV = SO3 + HSO3 + SO2L



      PH2O2F = ( PH2O20 + XL * DS4( 1 ) ) / ( 1.0 + XLH2O2 )
      PO3F   = ( PO30   + XL * DS4( 2 ) ) / ( 1.0 + XLO3   )
      PMHPF  = ( PMHP0  + XL * DS4( 4 ) ) / ( 1.0 + XLMHP  )
      PPAAF  = ( PPAA0  + XL * DS4( 5 ) ) / ( 1.0 + XLPAA  )

      PH2O2F = MAX( PH2O2F, 0.0 )
      PO3F   = MAX( PO3F,   0.0 )
      PMHPF  = MAX( PMHPF,  0.0 )
      PPAAF  = MAX( PPAAF,  0.0 )



      H2O2L = PH2O2F * H2O2H
      O3L   = PO3F   * O3H
      MHPL  = PMHPF  * MHPH
      PAAL  = PPAAF  * PAAH
      FOAL  = PFOAF  * FOAH
      NH3L  = PNH3F  * NH3H
      CO2L  = PCO2F  * CO2H
      HCLL  = PHCLF  * HCLH
      HNO3L = PHNO3F * HNO3H



      SO4COR  = SK6 * TS6COR / ( AE * GM2 + SK6 )
      HSO4COR = MAX( TS6COR - SO4COR, 0.0 )

      TS6ACC  = MAX( TS6  - TS6COR,   0.0 )
      SO4ACC  = MAX( SO4  - SO4COR,   0.0 )
      HSO4ACC = MAX( HSO4 - HSO4COR,  0.0 )
      NO3ACC  = MAX( NO3  - NO3COR,   0.0 )
      NAACC   = MAX( NA   - NACOR,    0.0 )
      CLACC   = MAX( CL   - CLCOR,    0.0 )



      LIQUID( LACL      ) = AC
      LIQUID( LNH4L     ) = NH4
      LIQUID( LCAL      ) = CA
      LIQUID( LNAACCL   ) = NAACC
      LIQUID( LOHL      ) = OH
      LIQUID( LSO4ACCL  ) = SO4ACC
      LIQUID( LHSO4ACCL ) = HSO4ACC
      LIQUID( LSO3L     ) = SO3
      LIQUID( LHSO3L    ) = HSO3
      LIQUID( LSO2L     ) = SO2L
      LIQUID( LCO3L     ) = CO3
      LIQUID( LHCO3L    ) = HCO3
      LIQUID( LCO2L     ) = CO2L
      LIQUID( LNO3ACCL  ) = NO3ACC
      LIQUID( LNH3L     ) = NH3L
      LIQUID( LCLACCL   ) = CLACC
      LIQUID( LH2O2L    ) = H2O2L
      LIQUID( LO3L      ) = O3L
      LIQUID( LFEL      ) = FE
      LIQUID( LMNL      ) = MN
      LIQUID( LAL       ) = A
      LIQUID( LFOAL     ) = FOAL
      LIQUID( LHCO2L    ) = HCO2
      LIQUID( LMHPL     ) = MHPL
      LIQUID( LPAAL     ) = PAAL
      LIQUID( LHCLL     ) = HCLL
      LIQUID( LORGAL    ) = ORGA
      LIQUID( LPRIML    ) = PRIM
      LIQUID( LMGL      ) = MG
      LIQUID( LKL       ) = K
      LIQUID( LBL       ) = B
      LIQUID( LHNO3L    ) = HNO3L
      LIQUID( LPRIMCORL ) = PRIMCOR
      LIQUID( LNUMCORL  ) = NUMCOR
      LIQUID( LTS6CORL  ) = TS6COR
      LIQUID( LNACORL   ) = NACOR
      LIQUID( LCLCORL   ) = CLCOR
      LIQUID( LNO3CORL  ) = NO3COR
      LIQUID( LORGPL    ) = ORGP
      LIQUID( LORGBL    ) = ORGB
      LIQUID( LECL      ) = EC




      IF ( TIMEW .LT. TAUCLD ) THEN




        DTRMV = 300.0
        IF ( ( CTHK1 .GT. 1.0E-10 ) .AND. ( PRCRATE .GT. 1.0E-10 ) ) &
           DTRMV = 3.6 * WTAVG * 1000.0 * CTHK1 / PRCRATE  
        DTRMV = MIN( DTRMV, 300.0 )
        ITERAT = ITERAT + 1



        TSIV = PSO20 * ONE_OVER_XL





        DSIVDT( 1 ) = -RH2O2 * H2O2L * HSO3 * AC / ( 0.1 + 13.0 * AC )
       
        TOTOX = PH2O20 * ONE_OVER_XL
        IF ( ( DSIVDT( 1 ) .EQ. 0.0 ) .OR. &
             ( TSIV  .LE. CONCMIN ) .OR. &
             ( TOTOX .LE. CONCMIN ) ) THEN
          DTW( 1 ) = DTRMV
        ELSE
          DTW( 1 ) = -0.05 * MIN( TOTOX, TSIV ) / DSIVDT( 1 )
        END IF









        DSIVDT( 2 ) = -( 2.4E4 * SO2L                                          + &
                        3.7E5 * EXP( -18.56 * ( ( 298.0E0 / TEMP ) - 1.0E0 ) ) * HSO3 + &
                        1.5E9 * EXP( -17.72 * ( ( 298.0E0 / TEMP ) - 1.0E0 ) ) * SO3 ) * O3L
        TOTOX = PO30 * ONE_OVER_XL
        IF ( ( DSIVDT( 2 ) .EQ. 0.0 ) .OR. &
             ( TSIV  .LE. CONCMIN ) .OR. &
             ( TOTOX .LE. CONCMIN ) ) THEN
          DTW( 2 ) = DTRMV
        ELSE
          DTW( 2 ) = -0.01 * MIN( TOTOX, TSIV ) / DSIVDT( 2 )
        END IF





























        DSIVDT( 3 ) = - ( 750.0E0  * MN * SIV +            &     
                         2600.0E0 * FE * SIV +             &    
                         1.0E10   * MN * FE * SIV )         

        IF ( ( DSIVDT( 3 ) .EQ. 0.0 ) .OR. ( TSIV .LE. CONCMIN ) ) THEN
          DTW( 3 ) = DTRMV
        ELSE
          DTW( 3 ) = -0.1 * TSIV / DSIVDT( 3 )
        END IF



        DSIVDT( 4 ) = -RMHP * AC * MHPL * HSO3
        TOTOX = PMHP0 * ONE_OVER_XL
        IF ( ( DSIVDT( 4 ) .EQ. 0.0 ) .OR. &
             ( TSIV  .LE. CONCMIN ) .OR. &
             ( TOTOX .LE. CONCMIN ) ) THEN
          DTW( 4 ) = DTRMV
        ELSE
          DTW( 4 ) = -0.1 * MIN( TOTOX, TSIV ) / DSIVDT( 4 )
        END IF





        DSIVDT( 5 ) = -( RPAA * AC + 7.00E2 ) * HSO3 * PAAL
        TOTOX = PPAA0 * ONE_OVER_XL
        IF ( ( DSIVDT( 5 ) .EQ. 0.0 ) .OR. &
             ( TSIV  .LE. CONCMIN ) .OR. &
             ( TOTOX .LE. CONCMIN ) ) THEN
          DTW( 5 ) = DTRMV
        ELSE
          DTW( 5 ) = -0.1 * MIN( TOTOX, TSIV ) / DSIVDT( 5 )
        END IF



        DSIVDT( 0 ) = 0.0
        DO IOX = 1, NUMOX
          DSIVDT( 0 ) = DSIVDT( 0 ) + DSIVDT( IOX )
        END DO



        DTW( 0 ) = MIN( DTW( 1 ), DTW( 2 ), DTW( 3 ), &
                        DTW( 4 ), DTW( 5 ) )



        IF ( DTW( 0 ) .GT. 8.0E+37 ) THEN
          WRITE(6,1001) PRCRATE, DSIVDT(0), TS6, DTW(0), CTHK1, WTAVG
        ELSE



60        DTS6 = ABS( DTW( 0 ) * ( -DSIVDT( 0 ) - TS6 * PRCRATE &
               / ( 3600.0 * CTHK1 * WTAVG ) ) )




          IF ( DTW( 0 ) .LE. TAUCLD ) THEN
            IF ( DTS6 .LT. 0.05 * TS6 ) THEN
              DTW( 0 ) = DTW( 0 ) * 2.0
	      GO TO 60
            END IF
          END IF
        END IF
        DTW( 0 ) = MIN( DTW( 0 ), DTRMV )




        IF ( DSIVDT( 0 ) .LT. 0.0 ) THEN
          DTW( 0 ) = MIN( DTW( 0 ), -TSIV * 1.00001 / DSIVDT( 0 ) )
        END IF




        IF ( TIMEW + DTW( 0 ) .GT. TAUCLD ) DTW( 0 ) = TAUCLD - TIMEW

        IF ( ITERAT .GT. 100 ) DTW( 0 ) = TAUCLD - TIMEW



        DTW( 0 ) = MIN( DTW( 0 ), TWASH )




        DO IOX = 0, NUMOX
          DS4( IOX ) = DS4( IOX ) + DTW( 0 ) * DSIVDT( IOX )
        END DO



        WETFAC = PRCRATE * FRACLIQ * DTW( 0 ) * SEC2HR
        DO LIQ = 1, NLIQS
          WETDEP( LIQ ) = WETDEP( LIQ ) + LIQUID( LIQ ) * WETFAC
        END DO

        TIMEW = TIMEW + DTW( 0 )



        GO TO 20
      END IF





      DEPSUM = ( WETDEP( LSO4ACCL ) + WETDEP( LHSO4ACCL ) ) * XC1

      IF ( ( TS6ACCA + TS6AKNA - DS4( 0 ) ) .NE. 0.0 ) THEN
        BETASO4 = DEPSUM / ( ( TS6ACCA + TS6AKNA - DS4( 0 ) ) * TAUCLD )
      ELSE
        BETASO4 = 0.0
      END IF

      EBETASO4T = EXP( -BETASO4 * TAUCLD )
      EALFA0T   = EXP( -ALFA0 * TAUCLD )
      EALFA2T   = EXP( -ALFA2 * TAUCLD )
      EALFA3T   = EXP( -ALFA3 * TAUCLD )



      TOTAMM = ( PNH3F  + ( NH4 + NH3L  ) * XL ) * RECIPAP1
      TOTNIT = ( PHNO3F + ( NO3ACC + HNO3L ) * XL ) * RECIPAP1



      GASWDEP( LSO2   ) = WETDEP( LSO3L  ) + WETDEP( LHSO3L ) &
                        + WETDEP( LSO2L  )
      GASWDEP( LNH3   ) = WETDEP( LNH3L  )
      GASWDEP( LH2O2  ) = WETDEP( LH2O2L )
      GASWDEP( LO3    ) = WETDEP( LO3L   )
      GASWDEP( LCO2   ) = WETDEP( LCO3L  ) + WETDEP( LHCO3L ) &
                        + WETDEP( LCO2L  )
      GASWDEP( LFOA   ) = WETDEP( LFOAL  ) + WETDEP( LHCO2L )
      GASWDEP( LMHP   ) = WETDEP( LMHPL  )
      GASWDEP( LPAA   ) = WETDEP( LPAAL  )
      GASWDEP( LHCL   ) = WETDEP( LHCLL  )
      GASWDEP( LHNO3  ) = WETDEP( LHNO3L )
      GASWDEP( LN2O5  ) = 0.0
      GASWDEP( LH2SO4 ) = 0.0



      GAS( LSO2   ) = ( PSO2F   + XL *  SIV )   * RECIPAP1
      GAS( LH2O2  ) = ( PH2O2F  + XL *  H2O2L ) * RECIPAP1
      GAS( LO3    ) = ( PO3F    + XL *  O3L )   * RECIPAP1
      GAS( LCO2   ) = ( PCO2F   + XL *  CO2L )  * RECIPAP1
      GAS( LFOA   ) = ( PFOAF   + XL * ( FOAL + HCO2 ) ) * RECIPAP1
      GAS( LMHP   ) = ( PMHPF   + XL *  MHPL )  * RECIPAP1
      GAS( LPAA   ) = ( PPAAF   + XL *  PAAL )  * RECIPAP1
      GAS( LHCL   ) = ( PHCLF   + XL *  HCLL )  * RECIPAP1

      GAS( LNH3   ) = FNH3  * TOTAMM
      GAS( LHNO3  ) = FHNO3 * TOTNIT
      GAS( LN2O5  ) = 0.0 
      GAS( LH2SO4 ) = 0.0 





      AERWDEP( LSO4AKN ) = 0.0
      AERWDEP( LNH4AKN ) = 0.0
      AERWDEP( LNO3AKN ) = 0.0
      AERWDEP( LECAKN  ) = 0.0
      AERWDEP( LPRIAKN ) = 0.0

      AERWDEP( LORGAAKN ) = 0.0
      AERWDEP( LORGPAKN ) = 0.0
      AERWDEP( LORGBAKN ) = 0.0

      AERWDEP( LSO4ACC ) = WETDEP( LSO4ACCL ) + WETDEP( LHSO4ACCL )
      AERWDEP( LNH4ACC ) = WETDEP( LNH4L    )
      AERWDEP( LNO3ACC ) = WETDEP( LNO3ACCL )
      AERWDEP( LECACC  ) = WETDEP( LECL     )
      AERWDEP( LPRIACC ) = WETDEP( LPRIML   )

      AERWDEP( LORGAACC ) = WETDEP( LORGAL )
      AERWDEP( LORGPACC ) = WETDEP( LORGPL )
      AERWDEP( LORGBACC ) = WETDEP( LORGBL )

      AERWDEP( LSO4COR ) = WETDEP( LTS6CORL  )
      AERWDEP( LNO3COR ) = WETDEP( LNO3CORL  )
      AERWDEP( LPRICOR ) = WETDEP( LPRIMCORL )

      IF ( AE_VRSN .EQ. 'AE3' ) THEN
        AERWDEP( LNACL  ) = WETDEP( LNACORL )
      ELSE
        AERWDEP( LNAAKN ) = 0.0
        AERWDEP( LCLAKN ) = 0.0
        AERWDEP( LNAACC ) = WETDEP( LNAACCL )
        AERWDEP( LCLACC ) = WETDEP( LCLACCL )
        AERWDEP( LNACOR ) = WETDEP( LNACORL )
        AERWDEP( LCLCOR ) = WETDEP( LCLCORL )
      END IF

      AERWDEP( LK      ) = WETDEP( LKL  )
      AERWDEP( LA3FE   ) = WETDEP( LFEL )
      AERWDEP( LB2MN   ) = WETDEP( LMNL )
      AERWDEP( LCACO3  ) = WETDEP( LCAL )
      AERWDEP( LMGCO3  ) = WETDEP( LMGL )

      AERWDEP( LNUMAKN ) = 0.0
      AERWDEP( LNUMACC ) = 0.0
      AERWDEP( LNUMCOR ) = 0.0
      AERWDEP( LSRFAKN ) = 0.0
      AERWDEP( LSRFACC ) = 0.0



      AEROSOL( LSO4AKN ) = AEROSOL( LSO4AKN ) * EALFA3T
      AEROSOL( LNH4AKN ) = AEROSOL( LNH4AKN ) * EALFA3T
      AEROSOL( LNO3AKN ) = AEROSOL( LNO3AKN ) * EALFA3T
      AEROSOL( LECAKN  ) = AEROSOL( LECAKN  ) * EALFA3T
      AEROSOL( LPRIAKN ) = AEROSOL( LPRIAKN ) * EALFA3T

      AEROSOL( LORGAAKN ) = AEROSOL( LORGAAKN ) * EALFA3T
      AEROSOL( LORGPAKN ) = AEROSOL( LORGPAKN ) * EALFA3T
      AEROSOL( LORGBAKN ) = AEROSOL( LORGBAKN ) * EALFA3T

      AEROSOL( LSO4ACC ) = TS6ACC * XL * RECIPAP1
      AEROSOL( LECACC  ) = EC     * XL * RECIPAP1
      AEROSOL( LPRIACC ) = PRIM   * XL * RECIPAP1

      AEROSOL( LORGAACC ) = ORGA * XL * RECIPAP1
      AEROSOL( LORGPACC ) = ORGP * XL * RECIPAP1
      AEROSOL( LORGBACC ) = ORGB * XL * RECIPAP1

      AEROSOL( LNH4ACC ) = FNH4ACC * TOTAMM
      AEROSOL( LNO3ACC ) = FNO3ACC * TOTNIT

      AEROSOL( LSO4COR ) = TS6COR * XL * RECIPAP1
      AEROSOL( LNO3COR ) = NO3COR * XL * RECIPAP1
      AEROSOL( LPRICOR ) = PRIMCOR* XL * RECIPAP1
      AEROSOL( LK      ) = K      * XL * RECIPAP1
      AEROSOL( LA3FE   ) = FE     * XL * RECIPAP1
      AEROSOL( LB2MN   ) = MN     * XL * RECIPAP1
      AEROSOL( LCACO3  ) = CA     * XL * RECIPAP1
      AEROSOL( LMGCO3  ) = MG     * XL * RECIPAP1

      IF ( AE_VRSN .EQ. 'AE3' ) THEN
        AEROSOL( LNACL  ) = NACOR * XL * RECIPAP1
      ELSE
        AEROSOL( LNAAKN ) = AEROSOL( LNAAKN ) * EALFA3T
        AEROSOL( LCLAKN ) = AEROSOL( LCLAKN ) * EALFA3T
        AEROSOL( LNAACC ) = NAACC * XL * RECIPAP1
        AEROSOL( LCLACC ) = CLACC * XL * RECIPAP1
        AEROSOL( LNACOR ) = NACOR * XL * RECIPAP1
        AEROSOL( LCLCOR ) = CLCOR * XL * RECIPAP1
      END IF

      AEROSOL( LNUMAKN ) = AEROSOL( LNUMAKN ) * EALFA0T
      AEROSOL( LNUMACC ) = AEROSOL( LNUMACC ) * EBETASO4T
      AEROSOL( LNUMCOR ) = NUMCOR * XL * RECIPAP1



      M3NEW = ( AEROSOL( LSO4ACC  ) * SGRAERMW( LSO4ACC  ) / 1.8e6 &
            +   AEROSOL( LNH4ACC  ) * SGRAERMW( LNH4ACC  ) / 1.8e6 &
            +   AEROSOL( LNO3ACC  ) * SGRAERMW( LNO3ACC  ) / 1.8e6 &
            +   AEROSOL( LORGPACC ) * SGRAERMW( LORGPACC ) / 2.0e6 &
            +   AEROSOL( LECACC   ) * SGRAERMW( LECACC   ) / 2.2e6 &
            +   AEROSOL( LPRIACC  ) * SGRAERMW( LPRIACC  ) / 2.2e6 &
            +   AEROSOL( LNAACC   ) * SGRAERMW( LNAACC   ) / 2.2e6 &
            +   AEROSOL( LCLACC   ) * SGRAERMW( LCLACC   ) / 2.2e6 )


      AEROSOL( LSRFAKN ) = AEROSOL( LSRFAKN ) * EALFA2T
      AEROSOL( LSRFACC ) = AEROSOL( LSRFACC ) &
                         * ( EXP( -BETASO4 * TAUCLD * ONETHIRD ) ) &
                         * ( M3NEW / MAX( M3OLD, CONCMIN) ) ** TWOTHIRDS



      HPWDEP = WETDEP( LACL )

      RETURN



1001  FORMAT (1X,'STORM RATE=', F6.3, 'DSIVDT(0) =', F10.5, &
             'TS6=', F10.5, 'DTW(0)=', F10.5, 'CTHK1=', F10.5, &
             'WTAVG=', F10.5)

      END SUBROUTINE AQCHEM

        INTEGER  FUNCTION  TRIMLEN ( STRING )

















      IMPLICIT NONE




        CHARACTER*(*)	STRING




        INTEGER  	L, K





        L = LEN( STRING )
        DO  11  K = L, 1, -1
            IF ( STRING( K:K ) .NE. ' ' ) THEN
                GO TO  12
            END IF
11      CONTINUE

        K = 1

12      CONTINUE

        TRIMLEN = K



END FUNCTION TRIMLEN


























      REAL FUNCTION HLCONST ( CNAME, TEMP, EFFECTIVE, HPLUS )


































      IMPLICIT NONE








      INTEGER, PARAMETER :: MXSPCS = 116     
      INTEGER, PARAMETER :: MXDSPCS = 17     



      INTEGER, PARAMETER :: LSO2       =  1  
      INTEGER, PARAMETER :: LHSO3      =  2  
      INTEGER, PARAMETER :: LHNO2      =  3  
      INTEGER, PARAMETER :: LHNO3      =  4  
      INTEGER, PARAMETER :: LCO2       =  5  
      INTEGER, PARAMETER :: LHCO3      =  6  
      INTEGER, PARAMETER :: LH2O2      =  7  
      INTEGER, PARAMETER :: LHCHO      =  8  
      INTEGER, PARAMETER :: LHCOOH     =  9  
      INTEGER, PARAMETER :: LHO2       = 10  
      INTEGER, PARAMETER :: LNH4OH     = 11  
      INTEGER, PARAMETER :: LH2O       = 12  
      INTEGER, PARAMETER :: LATRA      = 13  
      INTEGER, PARAMETER :: LCL2       = 14  
      INTEGER, PARAMETER :: LHOCL      = 15  ! HOCL
      INTEGER, PARAMETER :: LHCL       = 16  
      INTEGER, PARAMETER :: LHYDRAZINE = 17  



      CHARACTER*(*) CNAME               
      REAL          TEMP                
      LOGICAL       EFFECTIVE           
      REAL          HPLUS               



      CHARACTER(  7 ), SAVE :: PNAME = 'HLCONST'  
      CHARACTER( 16 ), SAVE :: SUBNAME( MXSPCS )  

      CHARACTER( 120 ) :: XMSG = ' '    

      INTEGER       SPC                 

      REAL          HPLUSI              
      REAL          HPLUS2I             
      REAL          CLMINUS             
      REAL          CLMINUSI            
      REAL          TFAC                
      REAL          AKEQ1               
      REAL          AKEQ2               
      REAL          OHION               
      REAL          KH                  





      REAL, SAVE :: A( MXSPCS )         
      REAL, SAVE :: E( MXSPCS )         




      REAL, SAVE :: B( MXDSPCS )        
      REAL, SAVE :: D( MXDSPCS )        

      DATA SUBNAME(  1), A(  1), E(  1) / 'O3              ', 1.14E-02, 2.3E+03 / 
      DATA SUBNAME(  2), A(  2), E(  2) / 'HO2             ', 4.0E+03, 5.9E+03 /  
      DATA SUBNAME(  3), A(  3), E(  3) / 'H2O2            ', 8.3E+04, 7.4E+03 /  
      DATA SUBNAME(  4), A(  4), E(  4) / 'NH3             ', 6.1E+01, 4.2E+03 /  
      DATA SUBNAME(  5), A(  5), E(  5) / 'NO              ', 1.9E-03, 1.4E+03 /  
      DATA SUBNAME(  6), A(  6), E(  6) / 'NO2             ', 1.2E-02, 2.5E+03 /  
      DATA SUBNAME(  7), A(  7), E(  7) / 'NO3             ', 0.6E+00, 0.0E+00 /  
      DATA SUBNAME(  8), A(  8), E(  8) / 'N2O5            ', 1.0E+30, 0.0E+00 /  
      DATA SUBNAME(  9), A(  9), E(  9) / 'HNO2            ', 5.0E+01, 4.9E+03 /  
      DATA SUBNAME( 10), A( 10), E( 10) / 'HNO3            ', 2.1E+05, 8.7E+03 /  
      DATA SUBNAME( 11), A( 11), E( 11) / 'HNO4            ', 1.2E+04, 6.9E+03 /  
      DATA SUBNAME( 12), A( 12), E( 12) / 'SO2             ', 1.4E+00, 2.9E+03 /  
      DATA SUBNAME( 13), A( 13), E( 13) / 'H2SO4           ', 1.0E+30, 0.0E+00 /  
      DATA SUBNAME( 14), A( 14), E( 14) / 'METHANE         ', 1.4E-03, 1.6E+03 /  
      DATA SUBNAME( 15), A( 15), E( 15) / 'ETHANE          ', 1.9E-03, 2.3E+03 /  
      DATA SUBNAME( 16), A( 16), E( 16) / 'PROPANE         ', 1.5E-03, 2.7E+03 /  
      DATA SUBNAME( 17), A( 17), E( 17) / 'BUTANE          ', 1.1E-03, 0.0E+00 /  
      DATA SUBNAME( 18), A( 18), E( 18) / 'PENTANE         ', 8.1E-04, 0.0E+00 /  
      DATA SUBNAME( 19), A( 19), E( 19) / 'HEXANE          ', 0.1E-03, 7.5E+03 /  
      DATA SUBNAME( 20), A( 20), E( 20) / 'OCTANE          ', 2.9E-03, 7.8E+03 /  
      DATA SUBNAME( 21), A( 21), E( 21) / 'NONANE          ', 2.4E-03, 2.1E+02 /  
      DATA SUBNAME( 22), A( 22), E( 22) / 'DECANE          ', 1.4E-04, 0.0E+00 /  
      DATA SUBNAME( 23), A( 23), E( 23) / 'ETHENE          ', 4.7E-03, 0.0E+00 /  
      DATA SUBNAME( 24), A( 24), E( 24) / 'PROPENE         ', 4.8E-03, 0.0E+00 /  
      DATA SUBNAME( 25), A( 25), E( 25) / 'ISOPRENE        ', 2.8E-02, 0.0E+00 /  
      DATA SUBNAME( 26), A( 26), E( 26) / 'ACETYLENE       ', 4.1E-02, 1.8E+03 /  
      DATA SUBNAME( 27), A( 27), E( 27) / 'BENZENE         ', 1.6E-01, 4.1E+03 /  
      DATA SUBNAME( 28), A( 28), E( 28) / 'TOLUENE         ', 1.5E-01, 4.0E+03 /  
      DATA SUBNAME( 29), A( 29), E( 29) / 'O-XYLENE        ', 1.9E-01, 4.0E+03 /  
      DATA SUBNAME( 30), A( 30), E( 30) / 'METHANOL        ', 2.2E+02, 5.2E+03 /  
      DATA SUBNAME( 31), A( 31), E( 31) / 'ETHANOL         ', 1.9E+02, 6.6E+03 /  
      DATA SUBNAME( 32), A( 32), E( 32) / '2-CRESOL        ', 8.2E+02, 0.0E+00 /  
      DATA SUBNAME( 33), A( 33), E( 33) / '4-CRESOL        ', 1.3E+02, 0.0E+00 /  
      DATA SUBNAME( 34), A( 34), E( 34) / 'METHYLHYDROPEROX', 3.1E+02, 5.2E+03 /  
      DATA SUBNAME( 35), A( 35), E( 35) / 'FORMALDEHYDE    ', 3.2E+03, 6.8E+03 /  
      DATA SUBNAME( 36), A( 36), E( 36) / 'ACETALDEHYDE    ', 1.4E+01, 5.6E+03 /  
      DATA SUBNAME( 37), A( 37), E( 37) / 'GENERIC_ALDEHYDE', 4.2E+03, 0.0E+00 /  
      DATA SUBNAME( 38), A( 38), E( 38) / 'GLYOXAL         ', 3.6E+05, 0.0E+00 /  
      DATA SUBNAME( 39), A( 39), E( 39) / 'ACETONE         ', 3.0E+01, 4.6E+03 /  
      DATA SUBNAME( 40), A( 40), E( 40) / 'FORMIC_ACID     ', 8.9E+03, 6.1E+03 /  
      DATA SUBNAME( 41), A( 41), E( 41) / 'ACETIC_ACID     ', 4.1E+03, 6.3E+03 /  
      DATA SUBNAME( 42), A( 42), E( 42) / 'METHYL_GLYOXAL  ', 3.2E+04, 0.0E+00 /  
      DATA SUBNAME( 43), A( 43), E( 43) / 'CO              ', 9.9E-04, 1.3E+03 /  
      DATA SUBNAME( 44), A( 44), E( 44) / 'CO2             ', 3.6E-02, 2.2E+03 /  
      DATA SUBNAME( 45), A( 45), E( 45) / 'PAN             ', 2.8E+00, 6.5E+03 /  
      DATA SUBNAME( 46), A( 46), E( 46) / 'MPAN            ', 1.7E+00, 0.0E+00 /  
      DATA SUBNAME( 47), A( 47), E( 47) / 'OH              ', 3.0E+01, 4.5E+03 /  
      DATA SUBNAME( 48), A( 48), E( 48) / 'METHYLPEROXY_RAD', 2.0E+03, 6.6E+03 /  
      DATA SUBNAME( 49), A( 49), E( 49) / 'PEROXYACETIC_ACI', 8.4E+02, 5.3E+03 /  
      DATA SUBNAME( 50), A( 50), E( 50) / 'PROPANOIC_ACID  ', 5.7E+03, 0.0E+00 /  
      DATA SUBNAME( 51), A( 51), E( 51) / '2-NITROPHENOL   ', 7.0E+01, 4.6E+03 /  
      DATA SUBNAME( 52), A( 52), E( 52) / 'PHENOL          ', 1.9E+03, 7.3E+03 /  
      DATA SUBNAME( 53), A( 53), E( 53) / 'BIACETYL        ', 7.4E+01, 5.7E+03 /  
      DATA SUBNAME( 54), A( 54), E( 54) / 'BENZALDEHYDE    ', 3.9E+01, 4.8E+03 /  
      DATA SUBNAME( 55), A( 55), E( 55) / 'PINENE          ', 4.9E-02, 0.0E+00 /  
      DATA SUBNAME( 56), A( 56), E( 56) / 'ATRA            ', 4.1E+05, 6.0E+03 /  
      DATA SUBNAME( 57), A( 57), E( 57) / 'DATRA           ', 4.1E+05, 6.0E+03 /  
      DATA SUBNAME( 58), A( 58), E( 58) / 'ADIPIC_ACID     ', 2.0E+08, 0.0E+00 /  
      DATA SUBNAME( 59), A( 59), E( 59) / 'ACROLEIN        ', 8.2E+00, 0.0E+00 /  
      DATA SUBNAME( 60), A( 60), E( 60) / '1,3-BUTADIENE   ', 1.4E-02, 0.0E+00 /  
      DATA SUBNAME( 61), A( 61), E( 61) / 'ACRYLONITRILE   ', 7.3E+00, 0.0E+00 /  
      DATA SUBNAME( 62), A( 62), E( 62) / 'CARBONTETRACHLOR', 3.4E-02, 4.2E+03 /  
      DATA SUBNAME( 63), A( 63), E( 63) / 'PROPYLENE_DICHLO', 3.4E-01, 4.3E+03 /  
      DATA SUBNAME( 64), A( 64), E( 64) / '1,3DICHLORPROPEN', 6.5E-01, 4.2E+03 /  
      DATA SUBNAME( 65), A( 65), E( 65) / '1,1,2,2-CL4ETHAN', 2.4E+00, 3.2E+03 /  
      DATA SUBNAME( 66), A( 66), E( 66) / 'CHLOROFORM      ', 2.5E-01, 4.5E+03 /  
      DATA SUBNAME( 67), A( 67), E( 67) / '1,2DIBROMOETHANE', 1.5E+00, 3.9E+03 /  
      DATA SUBNAME( 68), A( 68), E( 68) / '1,2DICHLOROETHAN', 7.3E-01, 4.2E+03 /  
      DATA SUBNAME( 69), A( 69), E( 69) / 'METHYLENE_CHLORI', 3.6E-01, 4.1E+03 /  
      DATA SUBNAME( 70), A( 70), E( 70) / 'PERCHLOROETHYLEN', 5.9E-02, 4.8E+03 /  
      DATA SUBNAME( 71), A( 71), E( 71) / 'TRICHLOROETHENE ', 1.0E-01, 4.6E+03 /  
      DATA SUBNAME( 72), A( 72), E( 72) / 'VINYL_CHLORIDE  ', 3.9E-02, 3.1E+03 /  
      DATA SUBNAME( 73), A( 73), E( 73) / 'ETHYLENE_OXIDE  ', 8.4E+00, 0.0E+00 /  
      DATA SUBNAME( 74), A( 74), E( 74) / 'PPN             ', 2.9E+00, 0.0E+00 /  
      DATA SUBNAME( 75), A( 75), E( 75) / 'NAPHTHALENE     ', 2.0E+00, 3.6E+03 /  
      DATA SUBNAME( 76), A( 76), E( 76) / 'QUINOLINE       ', 3.7E+03, 5.4E+03 /  
      DATA SUBNAME( 77), A( 77), E( 77) / 'MEK             ', 2.0E+01, 5.0E+03 /  
      DATA SUBNAME( 78), A( 78), E( 78) / 'MVK             ', 4.1E+01, 0.0E+00 /  
      DATA SUBNAME( 79), A( 79), E( 79) / 'METHACROLEIN    ', 6.5E+00, 0.0E+00 /  
      DATA SUBNAME( 80), A( 80), E( 80) / 'CL2             ', 8.6E-02, 2.0E+03 /  
      DATA SUBNAME( 81), A( 81), E( 81) / 'HOCL            ', 6.6E+02, 5.9E+03 /  
      DATA SUBNAME( 82), A( 82), E( 82) / 'HCL             ', 1.9E+01, 6.0E+02 /  
      DATA SUBNAME( 83), A( 83), E( 83) / 'FMCL            ', 1.1E+00, 0.0E+00 /  
      DATA SUBNAME( 84), A( 84), E( 84) / 'ICL1            ', 6.9E+01, 0.0E+00 /  
      DATA SUBNAME( 85), A( 85), E( 85) / 'ICL2            ', 6.9E+01, 0.0E+00 /  
      DATA SUBNAME( 86), A( 86), E( 86) / 'HG              ', 1.11E-01, 4.97E+03 /
      DATA SUBNAME( 87), A( 87), E( 87) / 'HGIIGAS         ', 1.41E+06, 5.26E+03 /
      DATA SUBNAME( 88), A( 88), E( 88) / 'TECDD_2378      ', 5.1E+00, 3.6E+03 /  
      DATA SUBNAME( 89), A( 89), E( 89) / 'PECDD_12378     ', 4.6E+00, 3.2E+03 /  
      DATA SUBNAME( 90), A( 90), E( 90) / 'HXCDD_123478    ', 8.1E+00, 2.9E+03 /  
      DATA SUBNAME( 91), A( 91), E( 91) / 'HXCDD_123678    ', 2.9E+00, 2.8E+03 /  
      DATA SUBNAME( 92), A( 92), E( 92) / 'HXCDD_123789    ', 6.5E+00, 2.7E+03 /  
      DATA SUBNAME( 93), A( 93), E( 93) / 'HPCDD_1234678   ', 1.2E+01, 2.4E+03 /  
      DATA SUBNAME( 94), A( 94), E( 94) / 'OTCDD           ', 9.8E+00, 2.3E+03 /  
      DATA SUBNAME( 95), A( 95), E( 95) / 'TECDF_2378      ', 8.5E+01, 3.7E+03 /  
      DATA SUBNAME( 96), A( 96), E( 96) / 'PECDF_12378     ', 5.2E+01, 2.9E+03 /  
      DATA SUBNAME( 97), A( 97), E( 97) / 'PECDF_23478     ', 1.8E+02, 3.0E+03 /  
      DATA SUBNAME( 98), A( 98), E( 98) / 'HXCDF_123478    ', 3.8E+01, 2.4E+03 /  
      DATA SUBNAME( 99), A( 99), E( 99) / 'HXCDF_123678    ', 9.0E+01, 2.9E+03 /  
      DATA SUBNAME(100), A(100), E(100) / 'HXCDF_234678    ', 1.0E+02, 2.6E+03 /  
      DATA SUBNAME(101), A(101), E(101) / 'HXCDF_123789    ', 5.6E+01, 2.6E+03 /  
      DATA SUBNAME(102), A(102), E(102) / 'HPCDF_1234678   ', 2.8E+01, 1.6E+03 /  
      DATA SUBNAME(103), A(103), E(103) / 'HPCDF_1234789   ', 8.0E+01, 2.1E+03 /  
      DATA SUBNAME(104), A(104), E(104) / 'OTCDF           ', 7.6E+01, 2.4E+03 /  
      DATA SUBNAME(105), A(105), E(105) / 'NAPHTHOL        ', 3.60E+03, 0.0E+00 / 
      DATA SUBNAME(106), A(106), E(106) / '1NITRONAPHTHALEN', 5.68E+02, 0.0E+00 / 
      DATA SUBNAME(107), A(107), E(107) / '2NITRONAPHTHALEN', 6.42E+02, 0.0E+00 / 
      DATA SUBNAME(108), A(108), E(108) / '14NAPHTHOQUINONE', 5.08E+05, 0.0E+00 / 
      DATA SUBNAME(109), A(109), E(109) / '2,4-TOLUENE_DIIS', 7.25E+00, 0.0E+00 / 
      DATA SUBNAME(110), A(110), E(110) / 'HEXAMETHYLE_DIIS', 2.08E+01, 0.0E+00 / 
      DATA SUBNAME(111), A(111), E(111) / 'HYDRAZINE       ', 1.14E+03, 0.0E+00 / 
      DATA SUBNAME(112), A(112), E(112) / 'MALEIC_ANHYDRIDE', 2.54E+02, 0.0E+00 / 
      DATA SUBNAME(113), A(113), E(113) / 'TRIETHYLAMINE   ', 6.71E+00, 0.0E+00 / 
      DATA SUBNAME(114), A(114), E(114) / 'P_DICHLOROBENZEN', 2.38E+00, 0.0E+00 / 
      DATA SUBNAME(115), A(115), E(115) / 'M-XYLENE        ', 1.43E-01, 3.9E+03 / 
      DATA SUBNAME(116), A(116), E(116) / 'P-XYLENE        ', 1.35E-01, 3.7E+03 / 

      DATA B( LSO2   ), D( LSO2   ) / 1.30E-02,  1.96E+03 /  
      DATA B( LHSO3  ), D( LHSO3  ) / 6.60E-08,  1.50E+03 /  
      DATA B( LHNO2  ), D( LHNO2  ) / 5.10E-04, -1.26E+03 /  
      DATA B( LHNO3  ), D( LHNO3  ) / 1.54E+01,  8.70E+03 /  
      DATA B( LCO2   ), D( LCO2   ) / 4.30E-07, -1.00E+03 /  
      DATA B( LHCO3  ), D( LHCO3  ) / 4.68E-11, -1.76E+03 /  
      DATA B( LH2O2  ), D( LH2O2  ) / 2.20E-12, -3.73E+03 /  
      DATA B( LHCHO  ), D( LHCHO  ) / 2.53E+03,  4.02E+03 /  
      DATA B( LHCOOH ), D( LHCOOH ) / 1.80E-04, -2.00E+01 /  
      DATA B( LHO2   ), D( LHO2   ) / 3.50E-05,  0.00E+00 /  
      DATA B( LNH4OH ), D( LNH4OH ) / 1.70E-05, -4.50E+02 /  
      DATA B( LH2O   ), D( LH2O   ) / 1.00E-14, -6.71E+03 /  
      DATA B( LATRA  ), D( LATRA  ) / 2.09E-02,  0.00E+00 /  
      DATA B( LCL2   ), D( LCL2   ) / 5.01E-04,  0.00E+00 /  ! CL2*H2O <=> HOCL + H + CL : LIN AND PEHKONEN, JGR, 103, D21, 28093-28102, NOVEMBER 20, 1998. ALSO SEE NOTE BELOW
      DATA B( LHOCL  ), D( LHOCL  ) / 3.16E-08,  0.00E+00 /  ! HOCL <=>H + OCL      : LIN AND PEHKONEN, JGR, 103, D21, 28093-28102, NOVEMBER 20, 1998
      DATA B( LHCL   ), D( LHCL   ) / 1.74E+06,  6.90E+03 /  
      DATA B( LHYDRAZINE), D( LHYDRAZINE) / 1.11E-08,  0.00E+00 /  

! Note for dissociation constant for equation 14: CL2*H2O <=> HOCL + H + CL












      SPC = INDEX1( CNAME, MXSPCS, SUBNAME )



      IF ( SPC .LE. 0 ) THEN
        XMSG = CNAME( 1:TRIMLEN( CNAME ) ) // ' not found in Henry''s Law Constant table, aborting.'

        write(0,*) ''
        write(0,*) PNAME,' : ',XMSG
        stop
      END IF



      TFAC = ( 298.0 - TEMP) / ( 298.0 * TEMP )
      KH = A( SPC ) * EXP( E( SPC ) * TFAC )
      HLCONST = KH



      IF ( EFFECTIVE ) THEN

        IF ( HPLUS .LE. 0.0 ) THEN
          XMSG = 'Negative or Zero [H+] concentration specified, aborting.'

          write(0,*) ''
          write(0,*) PNAME,' : ',XMSG
          stop
        END IF

        HPLUSI = 1.0 / HPLUS
        HPLUS2I = HPLUSI * HPLUSI



        CLMINUS   = 2.0E-03                
        CLMINUSI  = 1.0 / CLMINUS          

        CHECK_NAME: SELECT CASE ( CNAME( 1:TRIMLEN( CNAME ) ) )

        CASE ('SO2')            
                                

          AKEQ1 = B( LSO2 )  * EXP( D( LSO2 )  * TFAC )
          AKEQ2 = B( LHSO3 ) * EXP( D( LHSO3 ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI + AKEQ1 * AKEQ2 * HPLUS2I )

        CASE ('HNO2')           

          AKEQ1 = B( LHNO2 ) * EXP( D( LHNO2 ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        CASE ('HNO3')           

          AKEQ1 = B( LHNO3 ) * EXP( D( LHNO3 ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        CASE ('CO2')            
                                

          AKEQ1 = B( LCO2 )  * EXP( D( LCO2 )  * TFAC )
          AKEQ2 = B( LHCO3 ) * EXP( D( LHCO3 ) * TFAC )
          HLCONST = KH &
                  * ( 1.0 + AKEQ1 * HPLUSI + AKEQ1 * AKEQ2 * HPLUS2I )

        CASE ('H2O2')           

          AKEQ1 = B( LH2O2 ) * EXP( D( LH2O2 ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        CASE ('FORMALDEHYDE')   

          AKEQ1 = B( LHCHO ) * EXP( D( LHCHO ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 )

        CASE ('FORMIC_ACID')    

          AKEQ1 = B( LHCOOH ) * EXP( D( LHCOOH ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        CASE ('HO2')            

          AKEQ1 = B( LHO2 ) * EXP( D( LHO2 ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        CASE ('NH3')            

          AKEQ1 = B( LNH4OH ) * EXP( D( LNH4OH ) * TFAC )
          AKEQ2 = B( LH2O ) * EXP( D( LH2O ) * TFAC )
          OHION = AKEQ2 * HPLUSI
          HLCONST = KH * ( 1.0 + AKEQ1 / OHION )

        CASE ('HYDRAZINE')      

          AKEQ1 = B( LHYDRAZINE ) * EXP( D( LHYDRAZINE ) * TFAC )
          AKEQ2 = B( LH2O ) * EXP( D( LH2O ) * TFAC )
          OHION = AKEQ2 * HPLUSI
          HLCONST = KH * ( 1.0 + AKEQ1 / OHION )

        CASE ('ATRA', 'DATRA')  
                                

          AKEQ1   = B( LATRA ) * EXP( D( LATRA ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        CASE ( 'CL2' )          ! CL2*H2O <=> HOCL + H + CL
                                ! HOCL <=>H + OCL

          AKEQ1   = B( LCL2 )  * EXP( D( LCL2 ) * TFAC )
          AKEQ2   = B( LHOCL ) * EXP( D( LHOCL ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI * CLMINUSI &
                  + AKEQ1 * AKEQ2 * HPLUS2I * CLMINUSI )

        CASE ( 'HCL' )          

          AKEQ1   = B( LHCL ) * EXP( D( LHCL ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        CASE ( 'HOCL' )         ! HOCL <=> H+ + OCL-

          AKEQ1   = B( LHOCL ) * EXP( D( LHOCL ) * TFAC )
          HLCONST = KH * ( 1.0 + AKEQ1 * HPLUSI )

        END SELECT CHECK_NAME

      END IF

      RETURN
      END FUNCTION HLCONST








      INTEGER FUNCTION INDEX1 (NAME, N, NLIST)





















      IMPLICIT NONE



      CHARACTER*(*) NAME        
      INTEGER       N           
      CHARACTER*(*) NLIST(*)    



      INTEGER       I   




      DO 100 I = 1, N

          IF ( NAME .EQ. NLIST( I ) ) THEN    
              INDEX1 = I
              RETURN
          ENDIF

100   CONTINUE

      INDEX1 = 0        
      RETURN

END FUNCTION INDEX1

END MODULE module_ctrans_aqchem
