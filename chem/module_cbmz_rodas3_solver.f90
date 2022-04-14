







      module module_cbmz_rodas3_solver

      contains




      subroutine rodas3_ff_x2( n, t, tnext, hmin, hmax, hstart,   &
            y, abstol, reltol, yposlimit, yneglimit,   &
            yfixed, rconst,   &
            lu_nonzero_v, lu_crow_v, lu_diag_v, lu_icol_v,   &
            info, iok, lunerr,   &
            dydtsubr, yjacsubr, decompsubr, solvesubr )





































































      use module_peg_util, only:  peg_message, peg_error_fatal

      implicit none

      integer nvar_maxd, lu_nonzero_v_maxd
      parameter (nvar_maxd=99)
      parameter (lu_nonzero_v_maxd=999)




      integer    n, info(6), iok, lunerr, lu_nonzero_v
      integer    lu_crow_v(n), lu_diag_v(n), lu_icol_v(lu_nonzero_v)
      real       t, tnext, hmin, hmax, hstart
      real       y(n), abstol(n), reltol(n), yposlimit(n), yneglimit(n)
      real       yfixed(*), rconst(*)
      external   dydtsubr, yjacsubr, decompsubr, solvesubr


      logical    isreject, autonom
      integer    nfcn, njac, naccept, nreject, nnocnvg, i, j
      integer    ier

      real       k1(nvar_maxd), k2(nvar_maxd)
      real       k3(nvar_maxd), k4(nvar_maxd)
      real       f1(nvar_maxd), ynew(nvar_maxd)
      real       jac(lu_nonzero_v_maxd)
      real       ghinv, uround
      real       tin, tplus, h, hold, hlowest
      real       err, factor, facmax
      real       dround, c43, tau, x1, x2, ytol
      real       ylimit_solvesubr

      character*80 errmsg

      ylimit_solvesubr = 1.0e18


      if (n .gt. nvar_maxd) then
          call peg_message( lunerr, '*** rodas3 dimensioning problem' )
          write(errmsg,9050) 'n, nvar_maxd = ', n, nvar_maxd
          call peg_message( lunerr, errmsg )
          call peg_error_fatal( lunerr, '*** rodas3 fatal error' )
      else if (lu_nonzero_v .gt. lu_nonzero_v_maxd) then
          call peg_message( lunerr, '*** rodas3 dimensioning problem' )
          write(errmsg,9050) 'lu_nonvero_v, lu_nonzero_v_maxd = ',   &
                lu_nonzero_v, lu_nonzero_v_maxd
          call peg_message( lunerr, errmsg )
          call peg_error_fatal( lunerr, '*** rodas3 fatal error' )
      end if
9050  format( a, 2(1x,i6) )


      uround = 1.e-7
      dround = sqrt(uround)


      hlowest = dround
      if (info(1) .eq. 1) hlowest = 1.0e-7
      if (hmin .lt. hlowest) then
          call peg_message( lunerr, '*** rodas3 -- hmin is too small' )
          write(errmsg,9060) 'hmin and minimum allowed value = ',   &
                hmin, hlowest
          call peg_message( lunerr, errmsg )
          call peg_error_fatal( lunerr, '*** rodas3 fatal error' )
      else if (hmin .ge. hmax) then
          call peg_message( lunerr, '*** rodas3 -- hmin >= hmax' )
          write(errmsg,9060) 'hmin, hmax = ', hmin, hmax
          call peg_message( lunerr, errmsg )
          call peg_error_fatal( lunerr, '*** rodas3 fatal error' )
      end if
9060  format( a, 1p, 2e14.4 )


      autonom = info(1) .eq. 1

      c43 = - 8.e0/3.e0

      h = max( hstart, hmin )
      tplus = t
      tin = t
      isreject = .false.
      naccept  = 0
      nreject  = 0
      nnocnvg  = 0
      nfcn     = 0
      njac     = 0



 10    continue
       tplus = t + h
       if ( tplus .gt. tnext ) then
          h = tnext - t
          tplus = tnext
       end if

       call yjacsubr( n, t, y, jac, yfixed, rconst )
       njac = njac+1
       ghinv = -2.0e0/h
       do 20 j=1,n
         jac(lu_diag_v(j)) = jac(lu_diag_v(j)) + ghinv
 20    continue
       call decompsubr( n, jac, ier, lu_crow_v, lu_diag_v, lu_icol_v )

       if (ier.ne.0) then
         if ( h.gt.hmin) then
            h = 5.0e-1*h
            go to 10
         else


            iok = -1
            goto 200
         end if
       end if

       call dydtsubr( n, t, y, f1, yfixed, rconst )


       if (.not. autonom) then

         tau = dround*max( 1.0e-6, abs(t) )
         call dydtsubr( n, t+tau, y, k2, yfixed, rconst )
         nfcn=nfcn+1
         do 30 j = 1,n
           k3(j) = ( k2(j)-f1(j) )/tau
 30      continue


         x1 = 0.5*h
         do 40 j = 1,n
           k1(j) =  f1(j) + x1*k3(j)
           if (abs(k1(j)) .gt. ylimit_solvesubr) then
               iok = -(5000+j)
               goto 135
           end if
 40      continue
         call solvesubr( jac, k1 )


         x1 = 4.e0/h
         x2 = 1.5e0*h
         do 50 j = 1,n
           k2(j) = f1(j) - x1*k1(j) + x2*k3(j)
           if (abs(k2(j)) .gt. ylimit_solvesubr) then
               iok = -(5000+j)
               goto 135
           end if
 50      continue
         call solvesubr( jac, k2 )


       else

         do 60 j = 1,n
           k1(j) =  f1(j)
           if (abs(k1(j)) .gt. ylimit_solvesubr) then
               iok = -(5000+j)
               goto 135
           end if
 60      continue
         call solvesubr( jac, k1 )


         x1 = 4.e0/h
         do 70 j = 1,n
           k2(j) = f1(j) - x1*k1(j)
           if (abs(k2(j)) .gt. ylimit_solvesubr) then
               iok = -(5000+j)
               goto 135
           end if
 70      continue
         call solvesubr( jac, k2 )
       end if


       do 80 j = 1,n
         ynew(j) = y(j) - 2.0e0*k1(j)
 80    continue
       call dydtsubr( n, t+h, ynew, f1, yfixed, rconst )
       nfcn=nfcn+1
       do 90 j = 1,n
         k3(j) = f1(j) + ( -k1(j) + k2(j) )/h
         if (abs(k3(j)) .gt. ylimit_solvesubr) then
             iok = -(5000+j)
             goto 135
         end if
 90    continue
       call solvesubr( jac, k3 )


       do 100 j = 1,n
         ynew(j) = y(j) - 2.0e0*k1(j) - k3(j)
 100   continue
       call dydtsubr( n, t+h, ynew, f1, yfixed, rconst )
       nfcn=nfcn+1
       do 110 j = 1,n
         k4(j) = f1(j) + ( -k1(j) + k2(j) - c43*k3(j)  )/h
         if (abs(k4(j)) .gt. ylimit_solvesubr) then
             iok = -(5000+j)
             goto 135
         end if
 110   continue
       call solvesubr( jac, k4 )



       do 120 j = 1,n
         ynew(j) = y(j) - 2.0e0*k1(j) - k3(j) - k4(j)
 120   continue




        err=0.e0
        do 130 i=1,n
           ytol = abstol(i) + reltol(i)*abs(ynew(i))
           err = err + ( k4(i)/ytol )**2
 130    continue
        err = max( uround, sqrt( err/n ) )



        factor = 0.9/err**(1.e0/3.e0)
        if (isreject) then
            facmax=1.0
        else
            facmax=10.0
        end if
        factor = max( 1.0e-1, min(factor,facmax) )
        hold = h
        h = min( hmax, max(hmin,factor*h) )






	iok = 1
	do i = 1, n
	    if (y(i) .ne. y(i)) then
		iok = -(1000+i)
		goto 135
	    else if (y(i) .lt. yneglimit(i)) then
		iok = -(2000+i)
		goto 135
	    else if (y(i) .gt. yposlimit(i)) then
		iok = -(3000+i)
		goto 135
	    end if
	end do
135	if (iok .lt. 0) then
	    if (hold .le. hmin) then
		goto 200
	    else
		isreject = .true.
		nreject  = nreject+1
                h = max(hmin, 0.5*hold)
		if (t.eq.tin) h = max(hmin, 1.0e-1*hold)
		goto 10
	    end if
	end if
		



        if ( (err.gt.1).and.(hold.gt.hmin) ) then
          isreject = .true.
          nreject  = nreject+1
          if (t.eq.tin) h = max(hmin, 1.0e-1*hold)
        else
          isreject = .false.
          do 140 i=1,n
             y(i)  = ynew(i)
 140      continue
          t = tplus
          if (err.gt.1) then
            nnocnvg = nnocnvg+1
          else
            naccept = naccept+1
          end if
        end if


      if ( t .lt. tnext ) go to 10

      iok = 1
      if (nnocnvg .gt. 0) iok = 2



200   info(2) = nfcn
      info(3) = njac
      info(4) = naccept
      info(5) = nreject
      info(6) = nnocnvg

      return
      end subroutine rodas3_ff_x2                                


      end module module_cbmz_rodas3_solver
