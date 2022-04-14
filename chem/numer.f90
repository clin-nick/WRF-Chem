










      SUBROUTINE inter1(ng,xg,yg, n,x,y)
























      IMPLICIT NONE


      INTEGER n, ng
      REAL xg(ng)
      REAL x(n), y(n)


      REAL yg(ng)


      REAL slope
      INTEGER jsave, i, j

      jsave = 1
      DO i = 1, ng
         yg(i) = 0.
         j = jsave
   10    CONTINUE
         IF ((x(j) .GT. xg(i)) .OR. (xg(i) .GE. x(j+1))) THEN
            j = j+1
            IF (j .LE. n-1) GOTO 10

         ELSE
            slope = (y(j+1)-y(j)) / (x(j+1)-x(j))
            yg(i) = y(j) + slope * (xg(i) - x(j))
            jsave = j
         ENDIF
      ENDDO

      END SUBROUTINE inter1



      SUBROUTINE inter2(ng,xg,yg,n,x,y,ierr)































      IMPLICIT NONE


      INTEGER ng, n
      REAL x(n), y(n), xg(ng)


      REAL yg(ng)


      REAL area, xgl, xgu
      REAL darea, slope
      REAL a1, a2, b1, b2
      INTEGER ngintv
      INTEGER i, k, jstart
      INTEGER ierr

      ierr = 0



      DO i = 2, n
        IF (x(i) .LE. x(i-1)) THEN
          ierr = 1
          call wrf_debug( 0,'inter2: ERROR <<< x-grid not sorted' )
          RETURN
        ENDIF
      ENDDO

      DO i = 2, ng
        IF (xg(i) .LE. xg(i-1)) THEN
          ierr = 2
          call wrf_debug( 0,'inter2: ERROR <<< xg-grid not sorted!' )
          RETURN
        ENDIF
      ENDDO



      IF ( (x(1) .GT. xg(1)) .OR. (x(n) .LT. xg(ng)) ) THEN
        call wrf_error_fatal3("<stdin>",143,&
'inter2: <<<  Data do not span grid; Use ADDPNT to expand data and re-run.' )
      ENDIF





      jstart = 1
      ngintv = ng - 1
      DO i = 1,ngintv


            area = 0.0
            xgl = xg(i)
            xgu = xg(i+1)







            k = jstart
            IF (k .LE. n-1) THEN


   30         CONTINUE
              IF (x(k+1) .LE. xgl) THEN
                jstart = k - 1
                k = k+1
                IF (k .LE. n-1) GO TO 30
              ENDIF



   40         CONTINUE
              IF ((k .LE. n-1) .AND. (x(k) .LT. xgu)) THEN          
                jstart = k-1

                a1 = MAX(x(k),xgl)
                a2 = MIN(x(k+1),xgu)

                IF (x(k+1).EQ.x(k)) THEN
                  darea = 0.e0
                ELSE
                  slope = (y(k+1) - y(k))/(x(k+1) - x(k))
                  b1 = y(k) + slope*(a1 - x(k))
                  b2 = y(k) + slope*(a2 - x(k))
                  darea = (a2 - a1)*(b2 + b1)/2.
                ENDIF


                area = area + darea

                k = k+1
                GO TO 40
              ENDIF
            ENDIF


            yg(i) = area/(xgu - xgl)

      ENDDO

      END SUBROUTINE inter2



      SUBROUTINE inter3(ng,xg,yg, n,x,y, FoldIn)






































      IMPLICIT NONE
      

      INTEGER n, ng
      REAL xg(ng)
      REAL x(n), y(n)

      INTEGER FoldIn


      REAL yg(ng)


      REAL a1, a2, sum
      REAL tail
      INTEGER jstart, i, j, k


      IF ((FoldIn .NE. 0) .AND. (FoldIn .NE. 1)) THEN
         call wrf_error_fatal3("<stdin>",270,&
'inter3: ERROR <<<  Value for FOLDIN invalid. Must be 0 or 1' )
      ENDIF



      jstart = 1
      DO i = 1, ng - 1
         yg(i) = 0.
         sum = 0.
         j = jstart
         IF (j .LE. n-1) THEN
   20      CONTINUE

           IF (x(j+1) .LT. xg(i)) THEN
              jstart = j
              j = j+1
              IF (j .LE. n-1) GO TO 20
           ENDIF               

   25      CONTINUE

           IF ((x(j) .LE. xg(i+1)) .AND. (j .LE. n-1)) THEN
              a1 = MAX(x(j),xg(i))
              a2 = MIN(x(j+1),xg(i+1))
              sum = sum + y(j) * (a2-a1)/(x(j+1)-x(j))
              j = j+1
              GO TO 25
           ENDIF
           yg(i) = sum 
         ENDIF
      ENDDO
      



      IF (FoldIn .EQ. 1) THEN
         j = j-1
         a1 = xg(ng)            
         a2 = x(j+1)            

         IF ((a2 .GT. a1) .OR. (j+1 .LT. n)) THEN
            tail = y(j) * (a2-a1)/(x(j+1)-x(j))
            DO k = j+1, n-1
               tail = tail + y(k) * (x(k+1)-x(k))
            ENDDO
            yg(ng-1) = yg(ng-1) + tail
         ENDIF
      ENDIF

      END SUBROUTINE inter3



      SUBROUTINE inter4(ng,xg,yg, n,x,y, FoldIn)






































      IMPLICIT NONE
      

      INTEGER n, ng
      REAL xg(ng)
      REAL x(n), y(n)

      INTEGER FoldIn


      REAL yg(ng)


      REAL a1, a2, sum
      REAL tail
      INTEGER jstart, i, j, k


      IF ((FoldIn .NE. 0) .AND. (FoldIn .NE. 1)) THEN
         call wrf_error_fatal3("<stdin>",382,&
'inter3: ERROR <<<  Value for FOLDIN invalid. Must be 0 or 1' )
      ENDIF



      jstart = 1
      DO i = 1, ng - 1
         yg(i) = 0.
         sum = 0.
         j = jstart
         IF (j .LE. n-1) THEN
   20      CONTINUE
             IF (x(j+1) .LT. xg(i)) THEN
                jstart = j
                j = j+1
                IF (j .LE. n-1) GO TO 20
             ENDIF               
   25      CONTINUE
           IF ((x(j) .LE. xg(i+1)) .AND. (j .LE. n-1)) THEN
              a1 = MAX(x(j),xg(i))
              a2 = MIN(x(j+1),xg(i+1))
              sum = sum + y(j) * (a2-a1)
              j = j+1
              GO TO 25
           ENDIF
           yg(i) = sum /(xg(i+1)-xg(i))
        ENDIF
      ENDDO




      IF (FoldIn .EQ. 1) THEN
         j = j-1
         a1 = xg(ng)     
         a2 = x(j+1)     

         IF ((a2 .GT. a1) .OR. (j+1 .LT. n)) THEN
           tail = y(j) * (a2-a1)/(x(j+1)-x(j))
           DO k = j+1, n-1
              tail = tail + y(k) * (x(k+1)-x(k))
           ENDDO
           yg(ng-1) = yg(ng-1) + tail
         ENDIF
      ENDIF

      END SUBROUTINE inter4



      SUBROUTINE addpnt ( x, y, ld, n, xnew, ynew )
















      IMPLICIT NONE



      INTEGER, intent(in) :: ld
      INTEGER, intent(inout) :: n
      REAL, intent(inout)    :: x(ld), y(ld)
      REAL, intent(in) :: xnew, ynew



      INTEGER insert
      INTEGER i
      CHARACTER(len=256) :: emsg


      IF (n .GE. ld) THEN
         call wrf_error_fatal3("<stdin>",467,&
'addpnt: ERROR <<<  Cannot expand array All elements used.' )
      ENDIF

      insert = 1
      i = 2





 10   CONTINUE
      IF (i .LT. n) THEN
        IF (x(i) .LT. x(i-1)) THEN
           call wrf_error_fatal3("<stdin>",481,&
'addpnt: ERROR <<<  x-data must be in ascending order!' )
        ELSE
           IF (xnew .GT. x(i-1)) insert = i 
        ENDIF
        i = i+1
        GOTO 10
      ENDIF




      IF ( xnew .GT. x(n) ) THEN
         x(n+1) = xnew
         y(n+1) = ynew
      ELSE

         DO i = n, insert, -1
           x(i+1) = x(i)
           y(i+1) = y(i)
         ENDDO

         x(insert) = xnew
         y(insert) = ynew
      ENDIF



      n = n+1

      END SUBROUTINE addpnt
