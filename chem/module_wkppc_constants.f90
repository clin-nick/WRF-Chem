MODULE module_wkppc_constants




     REAL, PARAMETER  ::  navgdro = 6.022e23   

     REAL, PARAMETER  ::   &   
           mwh = 1.0079,  mwo = 15.9994,  mwair = 28.97


     REAL, PARAMETER :: mwh2o = 2*mwh + mwo



      REAL, PARAMETER ::          dens2con_a = 1.e-3     &
                                   * (1./mwair)          &
                                   * navgdro              
   



      REAL, PARAMETER ::          dens2con_w = 1.e-3     &
                                    * (1./mwh2o)         &
                                    * navgdro             











 
     REAL, PARAMETER ::  rtols=1.E-3 
     REAL, PARAMETER ::  atols=1.





END MODULE module_wkppc_constants
