
  ªF     k820309    a          16.0        åiXb                                                                                                           
       module_cmu_bulkaqchem.f90 MODULE_CMU_BULKAQCHEM                                                     
       PEG_ERROR_FATAL PEG_MESSAGE #         @                                                      #LUN    #STR              
                                                       
                                                     1 #         @                                                       #LUN    #STR              
                                                       
                                                     1 #         @                                                       #ISTAT_FATAL 	   #ISTAT_WARN 
   #TBEG_SEC    #TEND_SEC    #GAS    #AEROSOL    #GAS_AQFRAC    #TEMP    #P    #LWC    #RH    #CO2_MIXRAT_IN    #PHOTO_IN    #XPRESCRIBE_PH    #IRADICAL_IN    #IDECOMP_HMSA_HSO5    #IAQ    #YAQ_BEG    #YAQ_END    #PH_CMUAQ_CUR                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                D                                 	                      D                                 
                                                            
                                                       
                 D                                                    
     p          p            p                                    D                                                    
     p          p            p                                    D                                                    
     p          p            p                                    D @                                    
                                                       
                                                       
                                                       
                                                       
                                                       
                                                       
                                                                                                                               D                                                       D                                                    
     p          p            p                                    D                                                    
     p          p            p                                    D @                                    
       #         @                                                       #         @                                                      #TEMP    #AKEQ     #AKHEN !   #AKRE "                                                   
                 D                                                     
 W    p          p            p                                    D                                 !                   
 X    p          p            p                                    D                                 "     x              
 Y    p          p x           p x                         #         @                                  #                    #MEQN1 $   #Y %   #STIME &   #STOUT '   #RPAR (   #IPAR )                                                                                                                                                                                                                                              D @                               $                      D @                               %                   
     p          p            p                                    D @                               &     
                 D @                               '     
              @  D @                               (                    
     p          1     1                          @  D @                               )                         p          1     1                   #         @    @                             *                    #MEQN1 +   #T ,   #Y -   #F .   #RPAR /   #IPAR 0                                D @                               +                      D @                               ,     
                 D @                               -                   
     p          p            p                                    D @                               .                   
     p          p            p                                 @  D @                               /                    
     p          1     1                          @  D @                               0                         p          1     1                   #         @    @                             1                 	   #MEQN1 2   #T 3   #Y 4   #ML 5   #MU 6   #PD 7   #NROWPD 8   #RPAR 9   #IPAR :                                              2                                                       3     
                                                 4                    
     p          5  p        r 2       5  p        r 2                                                                5                                                       6                                                      7                    
       p        5  p        r 8   p          5  p        r 8     5  p        r 2       5  p        r 8     5  p        r 2                                                                8                   @                                   9                    
     p          1     1                          @                                   :                         p          1     1                   #         @                                  ;                    #MEQN1 <   #T =   #YAQ >   #AQPROD ?   #AQDEST @   #YAQPRIME A   #RPAR B   #IPAR C                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           <                                                       =     
                                                  >                   
     p          p            p                                    D                                 ?                   
 !    p          p            p                                    D                                 @                   
 "    p          p            p                                    D                                 A                   
      p          p            p                                 @  D                                 B                    
 #    p          1     1                          @  D @                               C                         p          1     1                   #         @                                  D                    #TEMP E   #QSAT F                                              E     
                 D                                 F     
       #         @                                  G                 
   #ACON H   #ASPRES I   #ACMET J   #AAKEQ K   #AAKHEN L   #AWV M   #ATEMP N   #AXSOL O   #XPRESCRIBE_PH P   #NERR_FULLEQUIL Q                                              H                   
 G    p          p            p                                                                     I                   
 H    p          p            p                                                                     J                   
 I    p          p            p                                                                     K                   
 J    p          p            p                                                                     L                   
 K    p          p            p                                                                     M     
                                                  N     
                 D                                 O     
                 D @                               P     
                 D                                 Q            #         @                                  R                    #X S   #CON T   #CMET U   #AKEQ V   #CC W                                              S     
                                                  T                   
 `    p          p            p                                                                     U                   
 a    p          p            p                                                                     V                   
 b    p          p            p                                    D                                 W     .              
 c    p          p .           p .                         #         @                                  X                    #C Y   #CMET Z   #CON [   #PHOTO \   #ZKRE ]   #RR ^   #ARYTM _                                                  D                                 Y     .              
 d    p          p .           p .                                   D                                 Z                   
 e    p          p            p                                                                     [                   
 f    p          p            p                                                                     \     
                                                  ]     x              
 g    p          p x           p x                                   D                                 ^     x              
 h    p          p x           p x                                   D                                 _     
       #         @                                  `                    #RR a   #ARYTM b   #RP c   #RL d                                              a     x              
 {    p          p x           p x                                                                    b     
                 D                                 c                   
 |    p          p            p                                    D                                 d                   
 }    p          p            p                          #         @                                  e                    #WL f   #RADIUS g   #TEMP h   #GCON i   #CON j   #C k   #AKEQ l   #AKHEN m   #FGL n   #FLG o   #GFGL p   #GFLG q                                              f     
                                                  g     
                                                  h     
                                                  i                   
 i    p          p            p                                                                     j                   
 j    p          p            p                                                                     k     .              
 k    p          p .           p .                                                                    l                   
 l    p          p            p                                                                     m                   
 m    p          p            p                                    D                                 n                   
 n    p          p            p                                    D                                 o                   
 o    p          p            p                                    D                                 p                   
 p    p          p            p                                    D                                 q                   
 q    p          p            p                          #         @                                  r                    #RP s   #RL t   #FGL u   #FLG v   #GFGL w   #GFLG x   #DP y   #DL z                                              s                   
 s    p          p            p                                                                     t                   
 t    p          p            p                                                                     u                   
 u    p          p            p                                                                     v                   
 v    p          p            p                                                                     w                   
 w    p          p            p                                                                     x                   
 x    p          p            p                                    D                                 y     1              
 y    p          p 1           p 1                                   D                                 z     1              
 z    p          p 1           p 1                         #         @                                   {                    #RADIUS |   #TEMP }   #C ~   #CON    #GCON    #AKEQ    #AKHEN    #AKRE                                               |     
                                                  }     
                 D                                 ~     .              
 @    p          p .           p .                                   D                                                    
 E    p          p            p                                                                                        
 A    p          p            p                                                                                        
 B    p          p            p                                                                                        
 C    p          p            p                                                                          x              
 D    p          p x           p x                         #         @                                                   
   #X    #F    #CON    #SPRES    #CMET    #ZKEQ    #ZKHEN    #WV    #TEMP    #XPRESCRIBE_PH                                                                        
                 D                                      
                                                                     
 Q    p          p            p                                                                                        
 R    p          p            p                                                                                        
 S    p          p            p                                                                                        
 T    p          p            p                                                                                        
 U    p          p            p                                                                          
                                                       
                                                       
              8      fn#fn     Ø   \   J  MODULE_PEG_UTIL 0   4  Z       PEG_ERROR_FATAL+MODULE_PEG_UTIL 4     @   a   PEG_ERROR_FATAL%LUN+MODULE_PEG_UTIL 4   Î  L   a   PEG_ERROR_FATAL%STR+MODULE_PEG_UTIL ,     Z       PEG_MESSAGE+MODULE_PEG_UTIL 0   t  @   a   PEG_MESSAGE%LUN+MODULE_PEG_UTIL 0   ´  L   a   PEG_MESSAGE%STR+MODULE_PEG_UTIL       p      AQOPERATOR1 (   p  @   a   AQOPERATOR1%ISTAT_FATAL '   °  @   a   AQOPERATOR1%ISTAT_WARN %   ð  @   a   AQOPERATOR1%TBEG_SEC %   0  @   a   AQOPERATOR1%TEND_SEC     p     a   AQOPERATOR1%GAS $        a   AQOPERATOR1%AEROSOL '        a   AQOPERATOR1%GAS_AQFRAC !   ,	  @   a   AQOPERATOR1%TEMP    l	  @   a   AQOPERATOR1%P     ¬	  @   a   AQOPERATOR1%LWC    ì	  @   a   AQOPERATOR1%RH *   ,
  @   a   AQOPERATOR1%CO2_MIXRAT_IN %   l
  @   a   AQOPERATOR1%PHOTO_IN *   ¬
  @   a   AQOPERATOR1%XPRESCRIBE_PH (   ì
  @   a   AQOPERATOR1%IRADICAL_IN .   ,  @   a   AQOPERATOR1%IDECOMP_HMSA_HSO5     l  @   a   AQOPERATOR1%IAQ $   ¬     a   AQOPERATOR1%YAQ_BEG $   @     a   AQOPERATOR1%YAQ_END )   Ô  @   a   AQOPERATOR1%PH_CMUAQ_CUR      H       DROPINIT    \  q       CONSTANTS    Í  @   a   CONSTANTS%TEMP         a   CONSTANTS%AKEQ     ¡     a   CONSTANTS%AKHEN    5     a   CONSTANTS%AKRE    É  e      AQINTEGR1     .  @   a   AQINTEGR1%MEQN1    n     a   AQINTEGR1%Y       @   a   AQINTEGR1%STIME     B  @   a   AQINTEGR1%STOUT         a   AQINTEGR1%RPAR         a   AQINTEGR1%IPAR             FEX1      @   a   FEX1%MEQN1    Y  @   a   FEX1%T         a   FEX1%Y    -     a   FEX1%F    Á     a   FEX1%RPAR    E     a   FEX1%IPAR    É         JEX    b  @   a   JEX%MEQN1    ¢  @   a   JEX%T    â  ´   a   JEX%Y      @   a   JEX%ML    Ö  @   a   JEX%MU      $  a   JEX%PD    :  @   a   JEX%NROWPD    z     a   JEX%RPAR    þ     a   JEX%IPAR      º      AQRATESA    <  @   a   AQRATESA%MEQN1    |  @   a   AQRATESA%T    ¼     a   AQRATESA%YAQ     P     a   AQRATESA%AQPROD     ä     a   AQRATESA%AQDEST "   x      a   AQRATESA%YAQPRIME    !     a   AQRATESA%RPAR    !     a   AQRATESA%IPAR    "  \       QSATURATION !   p"  @   a   QSATURATION%TEMP !   °"  @   a   QSATURATION%QSAT    ð"  Æ       FULLEQUIL    ¶#     a   FULLEQUIL%ACON !   J$     a   FULLEQUIL%ASPRES     Þ$     a   FULLEQUIL%ACMET     r%     a   FULLEQUIL%AAKEQ !   &     a   FULLEQUIL%AAKHEN    &  @   a   FULLEQUIL%AWV     Ú&  @   a   FULLEQUIL%ATEMP     '  @   a   FULLEQUIL%AXSOL (   Z'  @   a   FULLEQUIL%XPRESCRIBE_PH )   '  @   a   FULLEQUIL%NERR_FULLEQUIL    Ú'  t       VALUES    N(  @   a   VALUES%X    (     a   VALUES%CON    ")     a   VALUES%CMET    ¶)     a   VALUES%AKEQ    J*     a   VALUES%CC    Þ*  ¯       REACT    +     a   REACT%C    !,     a   REACT%CMET    µ,     a   REACT%CON    I-  @   a   REACT%PHOTO    -     a   REACT%ZKRE    .     a   REACT%RR    ±.  @   a   REACT%ARYTM    ñ.  k       ADDIT    \/     a   ADDIT%RR    ð/  @   a   ADDIT%ARYTM    00     a   ADDIT%RP    Ä0     a   ADDIT%RL    X1  »       MASS    2  @   a   MASS%WL    S2  @   a   MASS%RADIUS    2  @   a   MASS%TEMP    Ó2     a   MASS%GCON    g3     a   MASS%CON    û3     a   MASS%C    4     a   MASS%AKEQ    #5     a   MASS%AKHEN    ·5     a   MASS%FGL    K6     a   MASS%FLG    ß6     a   MASS%GFGL    s7     a   MASS%GFLG    8         DIFFER    8     a   DIFFER%RP    )9     a   DIFFER%RL    ½9     a   DIFFER%FGL    Q:     a   DIFFER%FLG    å:     a   DIFFER%GFGL    y;     a   DIFFER%GFLG    <     a   DIFFER%DP    ¡<     a   DIFFER%DL    5=         STEADY    Ì=  @   a   STEADY%RADIUS    >  @   a   STEADY%TEMP    L>     a   STEADY%C    à>     a   STEADY%CON    t?     a   STEADY%GCON    @     a   STEADY%AKEQ    @     a   STEADY%AKHEN    0A     a   STEADY%AKRE    ÄA  Â       ELECTRO    B  @   a   ELECTRO%X    ÆB  @   a   ELECTRO%F    C     a   ELECTRO%CON    C     a   ELECTRO%SPRES    .D     a   ELECTRO%CMET    ÂD     a   ELECTRO%ZKEQ    VE     a   ELECTRO%ZKHEN    êE  @   a   ELECTRO%WV    *F  @   a   ELECTRO%TEMP &   jF  @   a   ELECTRO%XPRESCRIBE_PH 