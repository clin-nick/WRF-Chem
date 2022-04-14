      MODULE ISRPIA





















































      USE module_data_isrpia_data

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
       PARAMETER (NCOMP=5, NIONS=7, NGASAQ=3, NSLDS=9, NPAIR=13,   &
                  NZSR=100, NERRMX=25)


      DOUBLE PRECISION  W(NCOMP), WAER(NCOMP), TEMP, RH, IPROB
      INTEGER METSTBL, NADJ

      DOUBLE PRECISION  DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL,   DRNANO3,   &
                        DRNA2SO4, DRNH4HS4, DRLC,     DRNH4NO3, DRNH4CL

      DOUBLE PRECISION :: DRMLCAB,DRMLCAS,DRMASAN,DRMG1,  DRMG2,     &
                    DRMG3,    DRMH1,    DRMH2,    DRMI1,    DRMI2,   &
                    DRMI3,    DRMQ1,    DRMR1,    DRMR2,    DRMR3,   &
                    DRMR4,    DRMR5,    DRMR6,    DRMR7,    DRMR8,   &
                    DRMR9,    DRMR10,   DRMR11,   DRMR12,   DRMR13
      INTEGER WFTYP










      CHARACTER SCASE*15
      DOUBLE PRECISION SULRATW, SULRAT, SODRAT




      CHARACTER*40 ERRMSG(NERRMX)
      INTEGER   ERRSTK(NERRMX), NOFER
      LOGICAL   STKOFL
      DATA ERRSTK/NERRMX*0/,   ERRMSG/NERRMX*' '/,  NOFER/0/,   &
           STKOFL/.FALSE./





















































































      END MODULE ISRPIA

