
  ñ&  O   k820309    a          16.0        ¥iXb                                                                                                           
       module_mosaic_soa_vbs.f90 MODULE_MOSAIC_SOA_VBS                                                     
       R8                                                                                                       #         @                                                      #MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE    #DTCHEM    #P_ATM    #T_K    #SWDOWN_CELL    #JAEROSOLSTATE    #AER     #GAS !   #WATER_A "   #AREA_WET_A #   #DP_WET_A $   #KG %   #SAT_SOA &   #TOTAL_SPECIES '   #MA (   #MC )   #MOSAIC_VARS_AA *                                                                                                                                                                                                                                                                                                          @               @                'X                   #IT_HOST    #IT_MOSAIC    #HOSTGRIDINFO    #IDIAGBB_HOST    #F_MOS_FAIL 	   #ISTEPS_ASTEM 
   #ISTEPS_ASTEM_MAX    #JASTEM_CALL    #JASTEM_FAIL    #JMESA_CALL    #JMESA_FAIL    #NITER_MESA_MAX    #NMAX_ASTEM    #NMAX_MESA    #FIX_ASTEM_NEGATIVE    #FLAG_ITR_KEL    #ZERO_WATER_FLAG    #CUMUL_STEPS_ASTEM    #NITER_MESA    #SWDOWN    #XNERR_ASTEM_NEGATIVE    #ITER_MESA                                                                                                                                                                                                            p          p            p                                                                                                                                     	     $                                                         
     (                                                              ,                                                              0                                                              4       	                                                       8       
                                                       <                                                              @                                                              D                                                              H                                                              L                                                              P                                                              T                                                             X          
                                                   `          
                                                   h          
                                                          p                 
  p          p          p            p          p                                                                                                            &                                                                                          
                                                      
                                                      
                                                      
                 
                                                         p          p            p                                    
D @                                    
             
     p Á        p p         p          p            p p         p          p                                    
D @                              !     k              
     p          p k           p k                                   
                                "                   
     p          p            p                                    
                                #                   
     p          p            p                                    
                                $                   
     p          p            p                                    
D @                              %     X             
     p l         p k         p            p k         p                                    
D @                              &     g              
     p          p g           p g                                   
D @                              '     g              
 	    p          p g           p g                                   
                                (     (              
     p          p          p            p          p                                    
                                )                    
 
    p          p          p            p          p                                    
                                 *     X              #MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE    #         @                                  +                   #EQUILIBRIUM%NBIN_A ,   #START_IND -   #N .   #AER /   #GAS 0   #KG 1   #SAT_SOA 2   #TOTAL_SPECIES 3                                                                                                                                                                                     ,                      
                                  -                     
                                  .                     
D                                /     
             
     p Á        p p         p          p            p p         p          p                                    
D                                0     k              
     p          p k           p k                                   
                                1     X             
     p l         p k         p            p k         p                                    
                                2     g              
     p          p g           p g                                   
D                                3     g              
     p          p g           p g                         #         @                                  4                    #N 5   #CTOT 6   #CSAT 7   #CA 8   #CGAS 9   #CPX :             D @                     @         5                     D @                              6                    
 (    p          5  p        r 5       5  p        r 5                              D @                              7                    
 )    p          5  p        r 5       5  p        r 5                              D @                              8                    
 *    p          5  p        r 5       5  p        r 5                              D                                9                    
 +    p          5  p        r 5       5  p        r 5                               D @                              :     
       #         @                                  ;                    #N <   #CTOT =   #CSAT >   #CA ?   #CPX @   #TOM A   #FVAL B                                    @         <                                                     =                    
 %    p          5  p        r <       5  p        r <                                                              >                    
 &    p          5  p        r <       5  p        r <                              D                                ?                    
 '    p          5  p        r <       5  p        r <                                                               @     
                                                 A     
                 D                                B     
       #         @                                   C                    #IPASS D                                                    
                                  D           #         @                                   E                    #T_K F   #PO_SOA G   #SAT_SOA H                                                                                                                                       
  @                              F     
                
D                                G     g              
 ,    p          p g           p g                                   
D                                H     g              
 -    p          p g           p g                         %         @                               I                    
       #PO_298 J   #DH K   #T L             
                                 J     
                
                                 K     
                
                                 L     
      %         @                               M                    
       #LOG10_CSAT_298 N             
                                 N     
             8      fn#fn (   Ø   C   J  MODULE_DATA_MOSAIC_KIND +     p       R8+MODULE_DATA_MOSAIC_KIND $     W      MOSAIC_SOA_VBS_INTR P   â  Ø     MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE+MODULE_DATA_MOSAIC_AERO X   º  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%IT_HOST+MODULE_DATA_MOSAIC_AERO Z     H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%IT_MOSAIC+MODULE_DATA_MOSAIC_AERO ]   J     a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%HOSTGRIDINFO+MODULE_DATA_MOSAIC_AERO ]   æ  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%IDIAGBB_HOST+MODULE_DATA_MOSAIC_AERO [   .  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%F_MOS_FAIL+MODULE_DATA_MOSAIC_AERO ]   v  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%ISTEPS_ASTEM+MODULE_DATA_MOSAIC_AERO a   ¾  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%ISTEPS_ASTEM_MAX+MODULE_DATA_MOSAIC_AERO \     H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%JASTEM_CALL+MODULE_DATA_MOSAIC_AERO \   N  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%JASTEM_FAIL+MODULE_DATA_MOSAIC_AERO [     H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%JMESA_CALL+MODULE_DATA_MOSAIC_AERO [   Þ  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%JMESA_FAIL+MODULE_DATA_MOSAIC_AERO _   &	  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%NITER_MESA_MAX+MODULE_DATA_MOSAIC_AERO [   n	  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%NMAX_ASTEM+MODULE_DATA_MOSAIC_AERO Z   ¶	  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%NMAX_MESA+MODULE_DATA_MOSAIC_AERO c   þ	  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%FIX_ASTEM_NEGATIVE+MODULE_DATA_MOSAIC_AERO ]   F
  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%FLAG_ITR_KEL+MODULE_DATA_MOSAIC_AERO `   
  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%ZERO_WATER_FLAG+MODULE_DATA_MOSAIC_AERO b   Ö
  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%CUMUL_STEPS_ASTEM+MODULE_DATA_MOSAIC_AERO [     H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%NITER_MESA+MODULE_DATA_MOSAIC_AERO W   f  H   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%SWDOWN+MODULE_DATA_MOSAIC_AERO e   ®  ¼   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%XNERR_ASTEM_NEGATIVE+MODULE_DATA_MOSAIC_AERO Z   j     a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA_TYPE%ITER_MESA+MODULE_DATA_MOSAIC_AERO +   þ  @   a   MOSAIC_SOA_VBS_INTR%DTCHEM *   >  @   a   MOSAIC_SOA_VBS_INTR%P_ATM (   ~  @   a   MOSAIC_SOA_VBS_INTR%T_K 0   ¾  @   a   MOSAIC_SOA_VBS_INTR%SWDOWN_CELL 2   þ     a   MOSAIC_SOA_VBS_INTR%JAEROSOLSTATE (     Ô   a   MOSAIC_SOA_VBS_INTR%AER (   f     a   MOSAIC_SOA_VBS_INTR%GAS ,   ú     a   MOSAIC_SOA_VBS_INTR%WATER_A /        a   MOSAIC_SOA_VBS_INTR%AREA_WET_A -   "     a   MOSAIC_SOA_VBS_INTR%DP_WET_A '   ¶  ´   a   MOSAIC_SOA_VBS_INTR%KG ,   j     a   MOSAIC_SOA_VBS_INTR%SAT_SOA 2   þ     a   MOSAIC_SOA_VBS_INTR%TOTAL_SPECIES '     ´   a   MOSAIC_SOA_VBS_INTR%MA '   F  ´   a   MOSAIC_SOA_VBS_INTR%MC 3   ú  u   a   MOSAIC_SOA_VBS_INTR%MOSAIC_VARS_AA    o  8      EQUILIBRIUM ;   §  @     EQUILIBRIUM%NBIN_A+MODULE_DATA_MOSAIC_AERO &   ç  @   a   EQUILIBRIUM%START_IND    '  @   a   EQUILIBRIUM%N     g  Ô   a   EQUILIBRIUM%AER     ;     a   EQUILIBRIUM%GAS    Ï  ´   a   EQUILIBRIUM%KG $        a   EQUILIBRIUM%SAT_SOA *        a   EQUILIBRIUM%TOTAL_SPECIES    «  ~       SOAP    )  @   a   SOAP%N    i  ´   a   SOAP%CTOT      ´   a   SOAP%CSAT    Ñ  ´   a   SOAP%CA      ´   a   SOAP%CGAS    9  @   a   SOAP%CPX    y         SPFCN       @   a   SPFCN%N    @  ´   a   SPFCN%CTOT    ô  ´   a   SPFCN%CSAT    ¨   ´   a   SPFCN%CA    \!  @   a   SPFCN%CPX    !  @   a   SPFCN%TOM    Ü!  @   a   SPFCN%FVAL $   "  z       SOA_VBS_LOAD_PARAMS *   "  @   a   SOA_VBS_LOAD_PARAMS%IPASS )   Ö"  ä       SOA_VBS_UPDATE_THERMCONS -   º#  @   a   SOA_VBS_UPDATE_THERMCONS%T_K 0   ú#     a   SOA_VBS_UPDATE_THERMCONS%PO_SOA 1   $     a   SOA_VBS_UPDATE_THERMCONS%SAT_SOA    "%  k       FN_PO    %  @   a   FN_PO%PO_298    Í%  @   a   FN_PO%DH    &  @   a   FN_PO%T    M&  d       DHR_APPROX *   ±&  @   a   DHR_APPROX%LOG10_CSAT_298 