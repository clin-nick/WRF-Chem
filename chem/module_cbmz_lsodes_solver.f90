








      module module_cbmz_lsodes_solver





































      contains





      subroutine lsodes_solver (   &
                  f, neq, y, t, tout, itol, rtol, atol, itask,   &
                  istate, iopt, rwork, lrw, iwork, liw, jac, mf,   &
                  ruserpar, nruserpar, iuserpar, niuserpar )
      external f, jac
      integer neq, itol, itask, istate, iopt, lrw, iwork, liw, mf
      integer nruserpar, iuserpar, niuserpar
      real y, t, tout, rtol, atol, rwork
      real ruserpar

      dimension neq(*), y(*), rtol(*), atol(*), rwork(lrw), iwork(liw)
      dimension ruserpar(nruserpar), iuserpar(niuserpar)




































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,   &
         mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,   &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,   &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,   &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, i1, i2, iflag, imax, imul, imxer, ipflag, ipgo, irem,   &
         j, kgo, lenrat, lenyht, leniw, lenrw, lf0, lia, lja,   &
         lrtem, lwtem, lyhd, lyhn, mf1, mord, mxhnl0, mxstp0, ncolm
      real rowns,   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real con0, conmin, ccmxj, psmall, rbig, seth



      real atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli,   &
         tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0
      dimension mord(2)
      logical ihit
















      common /ls0001/ rowns(209),   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,   &
         illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,   &
         mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),   &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu

      common /lss001/ con0, conmin, ccmxj, psmall, rbig, seth,   &
         iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,   &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,   &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,   &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu

      integer iok_vnorm
      common / lsodes_cmn_iok_vnorm / iok_vnorm

      data mord(1),mord(2)/12,5/, mxstp0/500/, mxhnl0/10/







      data lenrat/1/








      iok_vnorm = 1

      if (istate .lt. 1 .or. istate .gt. 3) go to 601
      if (itask .lt. 1 .or. itask .gt. 5) go to 602
      if (istate .eq. 1) go to 10
      if (init .eq. 0) go to 603
      if (istate .eq. 2) go to 200
      go to 20
 10   init = 0
      if (tout .eq. t) go to 430
 20   ntrep = 0











      if (neq(1) .le. 0) go to 604
      if (istate .eq. 1) go to 25
      if (neq(1) .gt. n) go to 605
 25   n = neq(1)
      if (itol .lt. 1 .or. itol .gt. 4) go to 606
      if (iopt .lt. 0 .or. iopt .gt. 1) go to 607
      moss = mf/100
      mf1 = mf - 100*moss
      meth = mf1/10
      miter = mf1 - 10*meth
      if (moss .lt. 0 .or. moss .gt. 2) go to 608
      if (meth .lt. 1 .or. meth .gt. 2) go to 608
      if (miter .lt. 0 .or. miter .gt. 3) go to 608
      if (miter .eq. 0 .or. miter .eq. 3) moss = 0

      if (iopt .eq. 1) go to 40
      maxord = mord(meth)
      mxstep = mxstp0
      mxhnil = mxhnl0
      if (istate .eq. 1) h0 = 0.0e0
      hmxi = 0.0e0
      hmin = 0.0e0
      seth = 0.0e0
      go to 60
 40   maxord = iwork(5)
      if (maxord .lt. 0) go to 611
      if (maxord .eq. 0) maxord = 100
      maxord = min0(maxord,mord(meth))
      mxstep = iwork(6)
      if (mxstep .lt. 0) go to 612
      if (mxstep .eq. 0) mxstep = mxstp0
      mxhnil = iwork(7)
      if (mxhnil .lt. 0) go to 613
      if (mxhnil .eq. 0) mxhnil = mxhnl0
      if (istate .ne. 1) go to 50
      h0 = rwork(5)
      if ((tout - t)*h0 .lt. 0.0e0) go to 614
 50   hmax = rwork(6)
      if (hmax .lt. 0.0e0) go to 615
      hmxi = 0.0e0
      if (hmax .gt. 0.0e0) hmxi = 1.0e0/hmax
      hmin = rwork(7)
      if (hmin .lt. 0.0e0) go to 616
      seth = rwork(8)
      if (seth .lt. 0.0e0) go to 609

 60   rtoli = rtol(1)
      atoli = atol(1)
      do 65 i = 1,n
        if (itol .ge. 3) rtoli = rtol(i)
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        if (rtoli .lt. 0.0e0) go to 619
        if (atoli .lt. 0.0e0) go to 620
 65     continue
















      lrat = lenrat
      if (istate .eq. 1) nyh = n
      lwmin = 0
      if (miter .eq. 1) lwmin = 4*n + 10*n/lrat
      if (miter .eq. 2) lwmin = 4*n + 11*n/lrat
      if (miter .eq. 3) lwmin = n + 2
      lenyh = (maxord+1)*nyh
      lrest = lenyh + 3*n
      lenrw = 20 + lwmin + lrest
      iwork(17) = lenrw
      leniw = 30
      if (moss .eq. 0 .and. miter .ne. 0 .and. miter .ne. 3)   &
         leniw = leniw + n + 1
      iwork(18) = leniw
      if (lenrw .gt. lrw) go to 617
      if (leniw .gt. liw) go to 618
      lia = 31
      if (moss .eq. 0 .and. miter .ne. 0 .and. miter .ne. 3)   &
         leniw = leniw + iwork(lia+n) - 1
      iwork(18) = leniw
      if (leniw .gt. liw) go to 618
      lja = lia + n + 1
      lia = min0(lia,liw)
      lja = min0(lja,liw)
      lwm = 21
      if (istate .eq. 1) nq = 1
      ncolm = min0(nq+1,maxord+2)
      lenyhm = ncolm*nyh
      lenyht = lenyh
      if (miter .eq. 1 .or. miter .eq. 2) lenyht = lenyhm
      imul = 2
      if (istate .eq. 3) imul = moss
      if (moss .eq. 2) imul = 3
      lrtem = lenyht + imul*n
      lwtem = lwmin
      if (miter .eq. 1 .or. miter .eq. 2) lwtem = lrw - 20 - lrtem
      lenwk = lwtem
      lyhn = lwm + lwtem
      lsavf = lyhn + lenyht
      lewt = lsavf + n
      lacor = lewt + n
      istatc = istate
      if (istate .eq. 1) go to 100









      lyhd = lyh - lyhn
      imax = lyhn - 1 + lenyhm

      if (lyhd) 70,80,74
 70   do 72 i = lyhn,imax
        j = imax + lyhn - i
 72     rwork(j) = rwork(j+lyhd)
      go to 80
 74   do 76 i = lyhn,imax
 76     rwork(i) = rwork(i+lyhd)
 80   lyh = lyhn
      iwork(22) = lyh
      if (miter .eq. 0 .or. miter .eq. 3) go to 92
      if (moss .ne. 2) go to 85

      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 82 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0e0) go to 621
 82     rwork(i+lewt-1) = 1.0e0/rwork(i+lewt-1)
 85   continue

      lsavf = min0(lsavf,lrw)
      lewt = min0(lewt,lrw)
      lacor = min0(lacor,lrw)
      call iprep (neq, y, rwork, iwork(lia), iwork(lja), ipflag, f, jac,   &
                  ruserpar, nruserpar, iuserpar, niuserpar)
      lenrw = lwm - 1 + lenwk + lrest
      iwork(17) = lenrw
      if (ipflag .ne. -1) iwork(23) = ipian
      if (ipflag .ne. -1) iwork(24) = ipjan
      ipgo = -ipflag + 1
      go to (90, 628, 629, 630, 631, 632, 633), ipgo
 90   iwork(22) = lyh
      if (lenrw .gt. lrw) go to 617

 92   jstart = -1
      if (n .eq. nyh) go to 200

      i1 = lyh + l*nyh
      i2 = lyh + (maxord + 1)*nyh - 1
      if (i1 .gt. i2) go to 200
      do 95 i = i1,i2
 95     rwork(i) = 0.0e0
      go to 200








 100  continue
      lyh = lyhn
      iwork(22) = lyh
      tn = t
      nst = 0
      h = 1.0e0
      nnz = 0
      ngp = 0
      nzl = 0
      nzu = 0

      do 105 i = 1,n
 105    rwork(i+lyh-1) = y(i)

      lf0 = lyh + nyh
      call f (neq, t, y, rwork(lf0),   &
          ruserpar, nruserpar, iuserpar, niuserpar)
      nfe = 1

      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 110 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0e0) go to 621
 110    rwork(i+lewt-1) = 1.0e0/rwork(i+lewt-1)
      if (miter .eq. 0 .or. miter .eq. 3) go to 120

      lacor = min0(lacor,lrw)
      call iprep (neq, y, rwork, iwork(lia), iwork(lja), ipflag, f, jac,   &
                  ruserpar, nruserpar, iuserpar, niuserpar)
      lenrw = lwm - 1 + lenwk + lrest
      iwork(17) = lenrw
      if (ipflag .ne. -1) iwork(23) = ipian
      if (ipflag .ne. -1) iwork(24) = ipjan
      ipgo = -ipflag + 1
      go to (115, 628, 629, 630, 631, 632, 633), ipgo
 115  iwork(22) = lyh
      if (lenrw .gt. lrw) go to 617

 120  continue
      if (itask .ne. 4 .and. itask .ne. 5) go to 125
      tcrit = rwork(1)
      if ((tcrit - tout)*(tout - t) .lt. 0.0e0) go to 625
      if (h0 .ne. 0.0e0 .and. (t + h0 - tcrit)*h0 .gt. 0.0e0)   &
         h0 = tcrit - t

 125  uround = r1mach(4)
      jstart = 0
      if (miter .ne. 0) rwork(lwm) = sqrt(uround)
      msbj = 50
      nslj = 0
      ccmxj = 0.2e0
      psmall = 1000.0e0*uround
      rbig = 0.01e0/psmall
      nhnil = 0
      nje = 0
      nlu = 0
      nslast = 0
      hu = 0.0e0
      nqu = 0
      ccmax = 0.3e0
      maxcor = 3
      msbp = 20
      mxncf = 10
















      lf0 = lyh + nyh
      if (h0 .ne. 0.0e0) go to 180
      tdist = abs(tout - t)
      w0 = amax1(abs(t),abs(tout))
      if (tdist .lt. 2.0e0*uround*w0) go to 622
      tol = rtol(1)
      if (itol .le. 2) go to 140
      do 130 i = 1,n
 130    tol = amax1(tol,rtol(i))
 140  if (tol .gt. 0.0e0) go to 160
      atoli = atol(1)
      do 150 i = 1,n
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        ayi = abs(y(i))
        if (ayi .ne. 0.0e0) tol = amax1(tol,atoli/ayi)
 150    continue
 160  tol = amax1(tol,100.0e0*uround)
      tol = amin1(tol,0.001e0)
      sum = vnorm (n, rwork(lf0), rwork(lewt))
      if (iok_vnorm .lt. 0) then
          istate = -901
          return
      end if
      sum = 1.0e0/(tol*w0*w0) + tol*sum**2
      h0 = 1.0e0/sqrt(sum)
      h0 = amin1(h0,tdist)
      h0 = sign(h0,tout-t)

 180  rh = abs(h0)*hmxi
      if (rh .gt. 1.0e0) h0 = h0/rh

      h = h0
      do 190 i = 1,n
 190    rwork(i+lf0-1) = h0*rwork(i+lf0-1)
      go to 270





 200  nslast = nst
      go to (210, 250, 220, 230, 240), itask
 210  if ((tn - tout)*h .lt. 0.0e0) go to 250
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 220  tp = tn - hu*(1.0e0 + 100.0e0*uround)
      if ((tp - tout)*h .gt. 0.0e0) go to 623
      if ((tn - tout)*h .lt. 0.0e0) go to 250
      go to 400
 230  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0e0) go to 624
      if ((tcrit - tout)*h .lt. 0.0e0) go to 625
      if ((tn - tout)*h .lt. 0.0e0) go to 245
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 240  tcrit = rwork(1)
      if ((tn - tcrit)*h .gt. 0.0e0) go to 624
 245  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0e0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0e0 + 4.0e0*uround)
      if ((tnext - tcrit)*h .le. 0.0e0) go to 250
      h = (tcrit - tn)*(1.0e0 - 4.0e0*uround)
      if (istate .eq. 2) jstart = -2











 250  continue
      if ((nst-nslast) .ge. mxstep) go to 500
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do 260 i = 1,n
        if (rwork(i+lewt-1) .le. 0.0e0) go to 510
 260    rwork(i+lewt-1) = 1.0e0/rwork(i+lewt-1)
 270  tolsf = uround*vnorm (n, rwork(lyh), rwork(lewt))
      if (tolsf .le. 1.0e0) go to 280

      tolsf = tolsf*2.0e0
      if (nst .eq. 0) go to 626
      go to 520
 280  if ((tn + h) .ne. tn) go to 290
      nhnil = nhnil + 1
      if (nhnil .gt. mxhnil) go to 290
      call xerrwv('lsodes-- warning..internal t (=r1) and h (=r2) are',   &
         50, 101, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(   &
        '      such that in the machine, t + h = t on the next step  ',   &
         60, 101, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv('      (h = step size). solver will continue anyway',   &
         50, 101, 0, 0, 0, 0, 2, tn, h)
      if (nhnil .lt. mxhnil) go to 290
      call xerrwv('lsodes-- above warning has been issued i1 times.  ',   &
         50, 102, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv('      it will not be issued again for this problem',   &
         50, 102, 0, 1, mxhnil, 0, 0, 0.0e0, 0.0e0)
 290  continue



      call stode_lsodes (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt),   &
         rwork(lsavf), rwork(lacor), rwork(lwm), rwork(lwm),   &
         f, jac, prjs, slss,   &
         ruserpar, nruserpar, iuserpar, niuserpar )
      kgo = 1 - kflag
      go to (300, 530, 540, 550), kgo





 300  init = 1
      go to (310, 400, 330, 340, 350), itask

 310  if ((tn - tout)*h .lt. 0.0e0) go to 250
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420

 330  if ((tn - tout)*h .ge. 0.0e0) go to 400
      go to 250

 340  if ((tn - tout)*h .lt. 0.0e0) go to 345
      call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
      t = tout
      go to 420
 345  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0e0*uround*hmx
      if (ihit) go to 400
      tnext = tn + h*(1.0e0 + 4.0e0*uround)
      if ((tnext - tcrit)*h .le. 0.0e0) go to 250
      h = (tcrit - tn)*(1.0e0 - 4.0e0*uround)
      jstart = -2
      go to 250

 350  hmx = abs(tn) + abs(h)
      ihit = abs(tn - tcrit) .le. 100.0e0*uround*hmx









 400  do 410 i = 1,n
 410    y(i) = rwork(i+lyh-1)
      t = tn
      if (itask .ne. 4 .and. itask .ne. 5) go to 420
      if (ihit) t = tcrit
 420  istate = 2
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      iwork(19) = nnz
      iwork(20) = ngp
      iwork(21) = nlu
      iwork(25) = nzl
      iwork(26) = nzu
      if (iok_vnorm .lt. 0) istate = -912
      return

 430  ntrep = ntrep + 1

      if (ntrep .lt. 5) then
          if (iok_vnorm .lt. 0) istate = -913
          return
      end if
      call xerrwv(   &
        'lsodes-- repeated calls with istate = 1 and tout = t (=r1)  ',   &
         60, 301, 0, 0, 0, 0, 1, t, 0.0e0)
      go to 800










 500  call xerrwv('lsodes-- at current t (=r1), mxstep (=i1) steps   ',   &
         50, 201, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv('      taken on this call before reaching tout     ',   &
         50, 201, 0, 1, mxstep, 0, 1, tn, 0.0e0)
      istate = -1
      go to 580

 510  ewti = rwork(lewt+i-1)
      call xerrwv('lsodes-- at t (=r1), ewt(i1) has become r2 .le. 0.',   &
         50, 202, 0, 1, i, 0, 2, tn, ewti)
      istate = -6
      go to 580

 520  call xerrwv('lsodes-- at t (=r1), too much accuracy requested  ',   &
         50, 203, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv('      for precision of machine..  see tolsf (=r2) ',   &
         50, 203, 0, 0, 0, 0, 2, tn, tolsf)
      rwork(14) = tolsf
      istate = -2
      go to 580

 530  call xerrwv('lsodes-- at t(=r1) and step size h(=r2), the error',   &
         50, 204, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv('      test failed repeatedly or with abs(h) = hmin',   &
         50, 204, 0, 0, 0, 0, 2, tn, h)
      istate = -4
      go to 560

 540  call xerrwv('lsodes-- at t (=r1) and step size h (=r2), the    ',   &
         50, 205, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv('      corrector convergence failed repeatedly     ',   &
         50, 205, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv('      or with abs(h) = hmin   ',   &
         30, 205, 0, 0, 0, 0, 2, tn, h)
      istate = -5
      go to 560

 550  call xerrwv('lsodes-- at t (=r1) and step size h (=r2), a fatal',   &
         50, 207, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv('      error flag was returned by cdrv (by way of  ',   &
         50, 207, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv('      subroutine prjs or slss)',   &
         30, 207, 0, 0, 0, 0, 2, tn, h)
      istate = -7
      go to 580

 560  big = 0.0e0
      imxer = 1
      do 570 i = 1,n
        size = abs(rwork(i+lacor-1)*rwork(i+lewt-1))
        if (big .ge. size) go to 570
        big = size
        imxer = i
 570    continue
      iwork(16) = imxer

 580  do 590 i = 1,n
 590    y(i) = rwork(i+lyh-1)
      t = tn
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      iwork(19) = nnz
      iwork(20) = ngp
      iwork(21) = nlu
      iwork(25) = nzl
      iwork(26) = nzu
      if (iok_vnorm .lt. 0) istate = -914
      return








 601  call xerrwv('lsodes-- istate (=i1) illegal ',   &
         30, 1, 0, 1, istate, 0, 0, 0.0e0, 0.0e0)
      go to 700
 602  call xerrwv('lsodes-- itask (=i1) illegal  ',   &
         30, 2, 0, 1, itask, 0, 0, 0.0e0, 0.0e0)
      go to 700
 603  call xerrwv('lsodes-- istate .gt. 1 but lsodes not initialized ',   &
         50, 3, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      go to 700
 604  call xerrwv('lsodes-- neq (=i1) .lt. 1     ',   &
         30, 4, 0, 1, neq(1), 0, 0, 0.0e0, 0.0e0)
      go to 700
 605  call xerrwv('lsodes-- istate = 3 and neq increased (i1 to i2)  ',   &
         50, 5, 0, 2, n, neq(1), 0, 0.0e0, 0.0e0)
      go to 700
 606  call xerrwv('lsodes-- itol (=i1) illegal   ',   &
         30, 6, 0, 1, itol, 0, 0, 0.0e0, 0.0e0)
      go to 700
 607  call xerrwv('lsodes-- iopt (=i1) illegal   ',   &
         30, 7, 0, 1, iopt, 0, 0, 0.0e0, 0.0e0)
      go to 700
 608  call xerrwv('lsodes-- mf (=i1) illegal     ',   &
         30, 8, 0, 1, mf, 0, 0, 0.0e0, 0.0e0)
      go to 700
 609  call xerrwv('lsodes-- seth (=r1) .lt. 0.0  ',   &
         30, 9, 0, 0, 0, 0, 1, seth, 0.0e0)
      go to 700
 611  call xerrwv('lsodes-- maxord (=i1) .lt. 0  ',   &
         30, 11, 0, 1, maxord, 0, 0, 0.0e0, 0.0e0)
      go to 700
 612  call xerrwv('lsodes-- mxstep (=i1) .lt. 0  ',   &
         30, 12, 0, 1, mxstep, 0, 0, 0.0e0, 0.0e0)
      go to 700
 613  call xerrwv('lsodes-- mxhnil (=i1) .lt. 0  ',   &
         30, 13, 0, 1, mxhnil, 0, 0, 0.0e0, 0.0e0)
      go to 700
 614  call xerrwv('lsodes-- tout (=r1) behind t (=r2)      ',   &
         40, 14, 0, 0, 0, 0, 2, tout, t)
      call xerrwv('      integration direction is given by h0 (=r1)  ',   &
         50, 14, 0, 0, 0, 0, 1, h0, 0.0e0)
      go to 700
 615  call xerrwv('lsodes-- hmax (=r1) .lt. 0.0  ',   &
         30, 15, 0, 0, 0, 0, 1, hmax, 0.0e0)
      go to 700
 616  call xerrwv('lsodes-- hmin (=r1) .lt. 0.0  ',   &
         30, 16, 0, 0, 0, 0, 1, hmin, 0.0e0)
      go to 700
 617  call xerrwv('lsodes-- rwork length is insufficient to proceed. ',   &
         50, 17, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(   &
        '        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)',   &
         60, 17, 0, 2, lenrw, lrw, 0, 0.0e0, 0.0e0)
      go to 700
 618  call xerrwv('lsodes-- iwork length is insufficient to proceed. ',   &
         50, 18, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(   &
        '        length needed is .ge. leniw (=i1), exceeds liw (=i2)',   &
         60, 18, 0, 2, leniw, liw, 0, 0.0e0, 0.0e0)
      go to 700
 619  call xerrwv('lsodes-- rtol(i1) is r1 .lt. 0.0        ',   &
         40, 19, 0, 1, i, 0, 1, rtoli, 0.0e0)
      go to 700
 620  call xerrwv('lsodes-- atol(i1) is r1 .lt. 0.0        ',   &
         40, 20, 0, 1, i, 0, 1, atoli, 0.0e0)
      go to 700
 621  ewti = rwork(lewt+i-1)
      call xerrwv('lsodes-- ewt(i1) is r1 .le. 0.0         ',   &
         40, 21, 0, 1, i, 0, 1, ewti, 0.0e0)
      go to 700
 622  call xerrwv(   &
        'lsodes-- tout (=r1) too close to t(=r2) to start integration',   &
         60, 22, 0, 0, 0, 0, 2, tout, t)
      go to 700
 623  call xerrwv(   &
        'lsodes-- itask = i1 and tout (=r1) behind tcur - hu (= r2)  ',   &
         60, 23, 0, 1, itask, 0, 2, tout, tp)
      go to 700
 624  call xerrwv(   &
        'lsodes-- itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   ',   &
         60, 24, 0, 0, 0, 0, 2, tcrit, tn)
      go to 700
 625  call xerrwv(   &
        'lsodes-- itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   ',   &
         60, 25, 0, 0, 0, 0, 2, tcrit, tout)
      go to 700
 626  call xerrwv('lsodes-- at start of problem, too much accuracy   ',   &
         50, 26, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(   &
        '      requested for precision of machine..  see tolsf (=r1) ',   &
         60, 26, 0, 0, 0, 0, 1, tolsf, 0.0e0)
      rwork(14) = tolsf
      go to 700
 627  call xerrwv('lsodes-- trouble from intdy. itask = i1, tout = r1',   &
         50, 27, 0, 1, itask, 0, 1, tout, 0.0e0)
      go to 700
 628  call xerrwv(   &
        'lsodes-- rwork length insufficient (for subroutine prep).   ',   &
         60, 28, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(   &
        '        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)',   &
         60, 28, 0, 2, lenrw, lrw, 0, 0.0e0, 0.0e0)
      go to 700
 629  call xerrwv(   &
        'lsodes-- rwork length insufficient (for subroutine jgroup). ',   &
         60, 29, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(   &
        '        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)',   &
         60, 29, 0, 2, lenrw, lrw, 0, 0.0e0, 0.0e0)
      go to 700
 630  call xerrwv(   &
        'lsodes-- rwork length insufficient (for subroutine odrv).   ',   &
         60, 30, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(   &
        '        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)',   &
         60, 30, 0, 2, lenrw, lrw, 0, 0.0e0, 0.0e0)
      go to 700
 631  call xerrwv(   &
        'lsodes-- error from odrv in yale sparse matrix package      ',   &
         60, 31, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      imul = (iys - 1)/n
      irem = iys - imul*n
      call xerrwv(   &
        '      at t (=r1), odrv returned error flag = i1*neq + i2.   ',   &
         60, 31, 0, 2, imul, irem, 1, tn, 0.0e0)
      go to 700
 632  call xerrwv(   &
        'lsodes-- rwork length insufficient (for subroutine cdrv).   ',   &
         60, 32, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      call xerrwv(   &
        '        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)',   &
         60, 32, 0, 2, lenrw, lrw, 0, 0.0e0, 0.0e0)
      go to 700
 633  call xerrwv(   &
        'lsodes-- error from cdrv in yale sparse matrix package      ',   &
         60, 33, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      imul = (iys - 1)/n
      irem = iys - imul*n
      call xerrwv(   &
        '      at t (=r1), cdrv returned error flag = i1*neq + i2.   ',   &
         60, 33, 0, 2, imul, irem, 1, tn, 0.0e0)
      if (imul .eq. 2) call xerrwv(   &
        '        duplicate entry in sparsity structure descriptors   ',   &
         60, 33, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)
      if (imul .eq. 3 .or. imul .eq. 6) call xerrwv(   &
        '        insufficient storage for nsfc (called by cdrv)      ',   &
         60, 33, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)

 700  if (illin .eq. 5) go to 710
      illin = illin + 1
      istate = -3
      if (iok_vnorm .lt. 0) istate = -915
      return
 710  call xerrwv('lsodes-- repeated occurrences of illegal input    ',   &
         50, 302, 0, 0, 0, 0, 0, 0.0e0, 0.0e0)

 800  call xerrwv('lsodes-- run aborted.. apparent infinite loop     ',   &
         50, 303, 2, 0, 0, 0, 0, 0.0e0, 0.0e0)
      if (iok_vnorm .lt. 0) istate = -916
      return

      end subroutine lsodes_solver
      subroutine adjlr (n, isp, ldif)
      integer n, isp, ldif

      dimension isp(*)







      integer ip, jlmax, jumax, lnfc, lsfc, nzlu

      ip = 2*n + 1

      jlmax = isp(ip)
      jumax = isp(ip+ip)

      nzlu = isp(n+1) - isp(1) + isp(ip+n+1) - isp(ip+1)
      lsfc = 12*n + 3 + 2*max0(jlmax,jumax)
      lnfc = 9*n + 2 + jlmax + jumax + nzlu
      ldif = max0(0, lsfc - lnfc)
      return

      end subroutine adjlr               
      subroutine cdrv   &
           (n, r,c,ic, ia,ja,a, b, z, nsp,isp,rsp,esp, path, flag)





































































































































































      integer  r(*), c(*), ic(*),  ia(*), ja(*),  isp(*), esp,  path,   &
         flag,  d, u, q, row, tmp, ar,  umax
      real  a(*), b(*), z(*), rsp(*)






      data lratio/1/

      if (path.lt.1 .or. 5.lt.path)  go to 111

      il   = 1
      ijl  = il  + (n+1)
      iu   = ijl +   n
      iju  = iu  + (n+1)
      irl  = iju +   n
      jrl  = irl +   n
      jl   = jrl +   n


      if ((path-1) * (path-5) .ne. 0)  go to 5
        max = (lratio*nsp + 1 - jl) - (n+1) - 5*n
        jlmax = max/2
        q     = jl   + jlmax
        ira   = q    + (n+1)
        jra   = ira  +   n
        irac  = jra  +   n
        iru   = irac +   n
        jru   = iru  +   n
        jutmp = jru  +   n
        jumax = lratio*nsp  + 1 - jutmp
        esp = max/lratio
        if (jlmax.le.0 .or. jumax.le.0)  go to 110

        do 1 i=1,n
          if (c(i).ne.i)  go to 2
   1      continue
        go to 3
   2    ar = nsp + 1 - n
        call  nroc   &
           (n, ic, ia,ja,a, isp(il), rsp(ar), isp(iu), flag)
        if (flag.ne.0)  go to 100

   3    call  nsfc   &
           (n, r, ic, ia,ja,   &
            jlmax, isp(il), isp(jl), isp(ijl),   &
            jumax, isp(iu), isp(jutmp), isp(iju),   &
            isp(q), isp(ira), isp(jra), isp(irac),   &
            isp(irl), isp(jrl), isp(iru), isp(jru),  flag)
        if(flag .ne. 0)  go to 100

        jlmax = isp(ijl+n-1)
        ju    = jl + jlmax
        jumax = isp(iju+n-1)
        if (jumax.le.0)  go to 5
        do 4 j=1,jumax
   4      isp(ju+j-1) = isp(jutmp+j-1)


   5  jlmax = isp(ijl+n-1)
      ju    = jl  + jlmax
      jumax = isp(iju+n-1)
      l     = (ju + jumax - 2 + lratio)  /  lratio    +    1
      lmax  = isp(il+n) - 1
      d     = l   + lmax
      u     = d   + n
      row   = nsp + 1 - n
      tmp   = row - n
      umax  = tmp - u
      esp   = umax - (isp(iu+n) - 1)

      if ((path-1) * (path-2) .ne. 0)  go to 6
        if (umax.lt.0)  go to 110
        call nnfc   &
           (n,  r, c, ic,  ia, ja, a, z, b,   &
            lmax, isp(il), isp(jl), isp(ijl), rsp(l),  rsp(d),   &
            umax, isp(iu), isp(ju), isp(iju), rsp(u),   &
            rsp(row), rsp(tmp),  isp(irl), isp(jrl),  flag)
        if(flag .ne. 0)  go to 100

   6  if ((path-3) .ne. 0)  go to 7
        call nnsc   &
           (n,  r, c,  isp(il), isp(jl), isp(ijl), rsp(l),   &
            rsp(d),    isp(iu), isp(ju), isp(iju), rsp(u),   &
            z, b,  rsp(tmp))

   7  if ((path-4) .ne. 0)  go to 8
        call nntc   &
           (n,  r, c,  isp(il), isp(jl), isp(ijl), rsp(l),   &
            rsp(d),    isp(iu), isp(ju), isp(iju), rsp(u),   &
            z, b,  rsp(tmp))
   8  return


 100  return

 110  flag = 10*n + 1
      return

 111  flag = 11*n + 1
      return
      end subroutine cdrv
      subroutine cfode (meth, elco, tesco)

      integer meth
      integer i, ib, nq, nqm1, nqp1
      real elco, tesco
      real agamq, fnq, fnqm1, pc, pint, ragq,   &
         rqfac, rq1fac, tsign, xpin
      dimension elco(13,12), tesco(3,12)


























      dimension pc(12)

      go to (100, 200), meth

 100  elco(1,1) = 1.0e0
      elco(2,1) = 1.0e0
      tesco(1,1) = 0.0e0
      tesco(2,1) = 2.0e0
      tesco(1,2) = 1.0e0
      tesco(3,12) = 0.0e0
      pc(1) = 1.0e0
      rqfac = 1.0e0
      do 140 nq = 2,12





        rq1fac = rqfac
        rqfac = rqfac/float(nq)
        nqm1 = nq - 1
        fnqm1 = float(nqm1)
        nqp1 = nq + 1

        pc(nq) = 0.0e0
        do 110 ib = 1,nqm1
          i = nqp1 - ib
 110      pc(i) = pc(i-1) + fnqm1*pc(i)
        pc(1) = fnqm1*pc(1)

        pint = pc(1)
        xpin = pc(1)/2.0e0
        tsign = 1.0e0
        do 120 i = 2,nq
          tsign = -tsign
          pint = pint + tsign*pc(i)/float(i)
 120      xpin = xpin + tsign*pc(i)/float(i+1)

        elco(1,nq) = pint*rq1fac
        elco(2,nq) = 1.0e0
        do 130 i = 2,nq
 130      elco(i+1,nq) = rq1fac*pc(i)/float(i)
        agamq = rqfac*xpin
        ragq = 1.0e0/agamq
        tesco(2,nq) = ragq
        if (nq .lt. 12) tesco(1,nqp1) = ragq*rqfac/float(nqp1)
        tesco(3,nqm1) = ragq
 140    continue
      return

 200  pc(1) = 1.0e0
      rq1fac = 1.0e0
      do 230 nq = 1,5





        fnq = float(nq)
        nqp1 = nq + 1

        pc(nqp1) = 0.0e0
        do 210 ib = 1,nq
          i = nq + 2 - ib
 210      pc(i) = pc(i-1) + fnq*pc(i)
        pc(1) = fnq*pc(1)

        do 220 i = 1,nqp1
 220      elco(i,nq) = pc(i)/pc(2)
        elco(2,nq) = 1.0e0
        tesco(1,nq) = rq1fac
        tesco(2,nq) = float(nqp1)/elco(1,nq)
        tesco(3,nq) = float(nq+2)/elco(1,nq)
        rq1fac = rq1fac/fnq
 230    continue
      return

      end subroutine cfode                    
      subroutine cntnzu (n, ia, ja, nzsut)
      integer n, ia, ja, nzsut

      dimension ia(*), ja(*)







      integer ii, jj, j, jmin, jmax, k, kmin, kmax, num

      num = 0
      do 50 ii = 1,n
        jmin = ia(ii)
        jmax = ia(ii+1) - 1
        if (jmin .gt. jmax) go to 50
        do 40 j = jmin,jmax
          if (ja(j) - ii) 10, 40, 30
 10       jj =ja(j)
          kmin = ia(jj)
          kmax = ia(jj+1) - 1
          if (kmin .gt. kmax) go to 30
          do 20 k = kmin,kmax
            if (ja(k) .eq. ii) go to 40
 20         continue
 30       num = num + 1
 40       continue
 50     continue
      nzsut = num
      return

      end subroutine cntnzu                   
      subroutine ewset (n, itol, rtol, atol, ycur, ewt)







      integer n, itol
      integer i
      real rtol, atol, ycur, ewt

      dimension rtol(*), atol(*), ycur(n), ewt(n)

      go to (10, 20, 30, 40), itol
 10   continue
      do 15 i = 1,n
 15     ewt(i) = rtol(1)*abs(ycur(i)) + atol(1)
      return
 20   continue
      do 25 i = 1,n
 25     ewt(i) = rtol(1)*abs(ycur(i)) + atol(i)
      return
 30   continue
      do 35 i = 1,n
 35     ewt(i) = rtol(i)*abs(ycur(i)) + atol(1)
      return
 40   continue
      do 45 i = 1,n
 45     ewt(i) = rtol(i)*abs(ycur(i)) + atol(i)
      return

      end subroutine ewset                                 
      subroutine intdy (t, k, yh, nyh, dky, iflag)

      integer k, nyh, iflag
      integer iownd, iowns,   &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, ic, j, jb, jb2, jj, jj1, jp1
      real t, yh, dky
      real rowns,   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real c, r, s, tp

      dimension yh(nyh,*), dky(*)
      common /ls0001/ rowns(209),   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,   &
         iownd(14), iowns(6),   &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu




















      iflag = 0
      if (k .lt. 0 .or. k .gt. nq) go to 80
      tp = tn - hu -  100.0e0*uround*(tn + hu)
      if ((t-tp)*(t-tn) .gt. 0.0e0) go to 90

      s = (t - tn)/h
      ic = 1
      if (k .eq. 0) go to 15
      jj1 = l - k
      do 10 jj = jj1,nq
 10     ic = ic*jj
 15   c = float(ic)
      do 20 i = 1,n
 20     dky(i) = c*yh(i,l)
      if (k .eq. nq) go to 55
      jb2 = nq - k
      do 50 jb = 1,jb2
        j = nq - jb
        jp1 = j + 1
        ic = 1
        if (k .eq. 0) go to 35
        jj1 = jp1 - k
        do 30 jj = jj1,j
 30       ic = ic*jj
 35     c = float(ic)
        do 40 i = 1,n
 40       dky(i) = c*yh(i,jp1) + s*dky(i)
 50     continue
      if (k .eq. 0) return
 55   r = h**(-k)
      do 60 i = 1,n
 60     dky(i) = r*dky(i)
      return

 80   call xerrwv('intdy--  k (=i1) illegal      ',   &
         30, 51, 0, 1, k, 0, 0, 0.0e0, 0.0e0)
      iflag = -1
      return
 90   call xerrwv('intdy--  t (=r1) illegal      ',   &
         30, 52, 0, 0, 0, 0, 1, t, 0.0e0)
      call xerrwv(   &
        '      t not in interval tcur - hu (= r1) to tcur (=r2)      ',   &
         60, 52, 0, 0, 0, 0, 2, tp, tn)
      iflag = -2
      return

      end subroutine intdy                            
      subroutine iprep (neq, y, rwork, ia, ja, ipflag, f, jac,   &
                        ruserpar, nruserpar, iuserpar, niuserpar )

      external f, jac
      integer neq, ia, ja, ipflag
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,   &
         mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,   &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,   &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,   &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, imax, lewtn, lyhd, lyhn
      integer nruserpar, iuserpar, niuserpar
      real y, rwork
      real rowns,   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real rlss
      real ruserpar

      dimension neq(*), y(*), rwork(*), ia(*), ja(*)
      dimension ruserpar(nruserpar), iuserpar(niuserpar)
      common /ls0001/ rowns(209),   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,   &
         illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,   &
         mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),   &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /lss001/ rlss(6),   &
         iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,   &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,   &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,   &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu













      ipflag = 0

      call prep_lsodes (neq, y, rwork(lyh), rwork(lsavf), rwork(lewt),   &
         rwork(lacor), ia, ja, rwork(lwm), rwork(lwm), ipflag, f, jac,   &
         ruserpar, nruserpar, iuserpar, niuserpar )
      lenwk = max0(lreq,lwmin)
      if (ipflag .lt. 0) return

      lyhn = lwm + lenwk
      if (lyhn .gt. lyh) return
      lyhd = lyh - lyhn
      if (lyhd .eq. 0) go to 20
      imax = lyhn - 1 + lenyhm
      do 10 i = lyhn,imax
 10     rwork(i) = rwork(i+lyhd)
      lyh = lyhn

 20   lsavf = lyh + lenyh
      lewtn = lsavf + n
      lacor = lewtn + n
      if (istatc .eq. 3) go to 40

      if (lewtn .gt. lewt) return
      do 30 i = 1,n
 30     rwork(i+lewtn-1) = rwork(i+lewt-1)
 40   lewt = lewtn
      return

      end subroutine iprep                                        
      subroutine jgroup (n,ia,ja,maxg,ngrp,igp,jgp,incl,jdone,ier)

      integer n, ia, ja, maxg, ngrp, igp, jgp, incl, jdone, ier

      dimension ia(*), ja(*), igp(*), jgp(n), incl(n), jdone(n)




















      integer i, j, k, kmin, kmax, ncol, ng

      ier = 0
      do 10 j = 1,n
 10     jdone(j) = 0
      ncol = 1
      do 60 ng = 1,maxg
        igp(ng) = ncol
        do 20 i = 1,n
 20       incl(i) = 0
        do 50 j = 1,n

          if (jdone(j) .eq. 1) go to 50
          kmin = ia(j)
          kmax = ia(j+1) - 1
          do 30 k = kmin,kmax

            i = ja(k)
            if (incl(i) .eq. 1) go to 50
 30         continue

          jgp(ncol) = j
          ncol = ncol + 1
          jdone(j) = 1
          do 40 k = kmin,kmax
            i = ja(k)
 40         incl(i) = 1
 50       continue

        if (ncol .eq. igp(ng)) go to 70
 60     continue

      if (ncol .le. n) go to 80
      ng = maxg
 70   ngrp = ng - 1
      return
 80   ier = 1
      return

      end subroutine jgroup                                           
      subroutine md   &
           (n, ia,ja, max, v,l, head,last,next, mark, flag)






















































































      integer  ia(*), ja(*),  v(*), l(*),  head(*), last(*), next(*),   &
         mark(*),  flag,  tag, dmin, vk,ek, tail
      equivalence  (vk,ek)


      tag = 0
      call  mdi   &
         (n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
      if (flag.ne.0)  return

      k = 0
      dmin = 1


   1  if (k.ge.n)  go to 4


   2    if (head(dmin).gt.0)  go to 3
          dmin = dmin + 1
          go to 2


   3    vk = head(dmin)
        head(dmin) = next(vk)
        if (head(dmin).gt.0)  last(head(dmin)) = -dmin


        k = k+1
        next(vk) = -k
        last(ek) = dmin - 1
        tag = tag + last(ek)
        mark(vk) = tag


        call  mdm   &
           (vk,tail, v,l, last,next, mark)


        call  mdp   &
           (k,ek,tail, v,l, head,last,next, mark)


        call  mdu   &
           (ek,dmin, v,l, head,last,next, mark)

        go to 1


   4  do 5 k=1,n
        next(k) = -next(k)
   5    last(next(k)) = k

      return
      end subroutine md
      subroutine mdi   &
           (n, ia,ja, max,v,l, head,last,next, mark,tag, flag)






      integer  ia(*), ja(*),  v(*), l(*),  head(*), last(*), next(*),   &
         mark(*), tag,  flag,  sfs, vi,dvi, vj


      do 1 vi=1,n
        mark(vi) = 1
        l(vi) = 0
   1    head(vi) = 0
      sfs = n+1



      do 6 vi=1,n
        jmin = ia(vi)
        jmax = ia(vi+1) - 1
        if (jmin.gt.jmax)  go to 6
        do 5 j=jmin,jmax
          vj = ja(j)
          if (vj-vi) 2, 5, 4



   2      lvk = vi
          kmax = mark(vi) - 1
          if (kmax .eq. 0) go to 4
          do 3 k=1,kmax
            lvk = l(lvk)
            if (v(lvk).eq.vj) go to 5
   3        continue

   4        if (sfs.ge.max)  go to 101


            mark(vi) = mark(vi) + 1
            v(sfs) = vj
            l(sfs) = l(vi)
            l(vi) = sfs
            sfs = sfs+1


            mark(vj) = mark(vj) + 1
            v(sfs) = vi
            l(sfs) = l(vj)
            l(vj) = sfs
            sfs = sfs+1
   5      continue
   6    continue


      do 7 vi=1,n
        dvi = mark(vi)
        next(vi) = head(dvi)
        head(dvi) = vi
        last(vi) = -dvi
        nextvi = next(vi)
        if (nextvi.gt.0)  last(nextvi) = vi
   7    mark(vi) = tag

      return


 101  flag = 9*n + vi
      return
      end subroutine mdi
      subroutine mdm   &
           (vk,tail, v,l, last,next, mark)






      integer  vk, tail,  v(*), l(*),   last(*), next(*),   mark(*),   &
         tag, s,ls,vs,es, b,lb,vb, blp,blpmax
      equivalence  (vs, es)


      tag = mark(vk)
      tail = vk


      ls = l(vk)
   1  s = ls
      if (s.eq.0)  go to 5
        ls = l(s)
        vs = v(s)
        if (next(vs).lt.0)  go to 2



          mark(vs) = tag
          l(tail) = s
          tail = s
          go to 4



   2      lb = l(es)
          blpmax = last(es)
          do 3 blp=1,blpmax
            b = lb
            lb = l(b)
            vb = v(b)



            if (mark(vb).ge.tag)  go to 3
              mark(vb) = tag
              l(tail) = b
              tail = b
   3        continue


          mark(es) = tag

   4    go to 1


   5  l(tail) = 0

      return
      end subroutine mdm
      subroutine mdp   &
           (k,ek,tail, v,l, head,last,next, mark)






      integer  ek, tail,  v(*), l(*),  head(*), last(*), next(*),   &
         mark(*),  tag, free, li,vi,lvi,evi, s,ls,es, ilp,ilpmax


      tag = mark(ek)


      li = ek
      ilpmax = last(ek)
      if (ilpmax.le.0)  go to 12
      do 11 ilp=1,ilpmax
        i = li
        li = l(i)
        vi = v(li)


        if (last(vi).eq.0)  go to 3
          if (last(vi).gt.0)  go to 1
            head(-last(vi)) = next(vi)
            go to 2
   1        next(last(vi)) = next(vi)
   2      if (next(vi).gt.0)  last(next(vi)) = last(vi)


   3    ls = vi
   4    s = ls
        ls = l(s)
        if (ls.eq.0)  go to 6
          es = v(ls)
          if (mark(es).lt.tag)  go to 5
            free = ls
            l(s) = l(ls)
            ls = s
   5      go to 4


   6    lvi = l(vi)
        if (lvi.ne.0)  go to 7
          l(i) = l(li)
          li = i

          k = k+1
          next(vi) = -k
          last(ek) = last(ek) - 1
          go to 11



   7      if (l(lvi).ne.0)  go to 9
            evi = v(lvi)
            if (next(evi).ge.0)  go to 9
              if (mark(evi).lt.0)  go to 8




                last(vi) = evi
                mark(evi) = -1
                l(tail) = li
                tail = li
                l(i) = l(li)
                li = i
                go to 10



   8            last(vi) = 0
                mark(evi) = mark(evi) - 1
                go to 10


   9            last(vi) = -ek


  10      v(free) = ek
          l(free) = l(vi)
          l(vi) = free
  11    continue


  12  l(tail) = 0

      return
      end subroutine mdp
      subroutine mdu   &
           (ek,dmin, v,l, head,last,next, mark)







      integer  ek, dmin,  v(*), l(*),  head(*), last(*), next(*),   &
         mark(*),  tag, vi,evi,dvi, s,vs,es, b,vb, ilp,ilpmax,   &
         blp,blpmax
      equivalence  (vs, es)


      tag = mark(ek) - last(ek)


      i = ek
      ilpmax = last(ek)
      if (ilpmax.le.0)  go to 11
      do 10 ilp=1,ilpmax
        i = l(i)
        vi = v(i)
        if (last(vi))  1, 10, 8



   1      tag = tag + 1
          dvi = last(ek)


          s = l(vi)
   2      s = l(s)
          if (s.eq.0)  go to 9
            vs = v(s)
            if (next(vs).lt.0)  go to 3


              mark(vs) = tag
              dvi = dvi + 1
              go to 5



   3          if (mark(es).lt.0)  go to 6


              b = es
              blpmax = last(es)
              do 4 blp=1,blpmax
                b = l(b)
                vb = v(b)


                if (mark(vb).ge.tag)  go to 4
                  mark(vb) = tag
                  dvi = dvi + 1
   4            continue

   5        go to 2



   6      last(vi) = 0
          mark(es) = mark(es) - 1
   7      s = l(s)
          if (s.eq.0)  go to 10
            es = v(s)
            if (mark(es).lt.0)  mark(es) = mark(es) - 1
            go to 7



   8      evi = last(vi)
          dvi = last(ek) + last(evi) + mark(evi)
          mark(evi) = 0


   9    next(vi) = head(dvi)
        head(dvi) = vi
        last(vi) = -dvi
        if (next(vi).gt.0)  last(next(vi)) = vi
        if (dvi.lt.dmin)  dmin = dvi

  10    continue

  11  return
      end subroutine mdu
      subroutine nnfc   &
           (n, r,c,ic, ia,ja,a, z, b,   &
            lmax,il,jl,ijl,l, d, umax,iu,ju,iju,u,   &
            row, tmp, irl,jrl, flag)































      integer rk,umax




      integer  r(*), c(*), ic(*), ia(*), ja(*), il(*), jl(*), ijl(*)
      integer  iu(*), ju(*), iju(*), irl(*), jrl(*), flag
      real  a(*), l(*), d(*), u(*), z(*), b(*), row(*)
      real tmp(*), lki, sum, dk




      if(il(n+1)-1 .gt. lmax) go to 104
      if(iu(n+1)-1 .gt. umax) go to 107
      do 1 k=1,n
        irl(k) = il(k)
        jrl(k) = 0
   1    continue


      do 19 k=1,n

        row(k) = 0
        i1 = 0
        if (jrl(k) .eq. 0) go to 3
        i = jrl(k)
   2    i2 = jrl(i)
        jrl(i) = i1
        i1 = i
        row(i) = 0
        i = i2
        if (i .ne. 0) go to 2

   3    jmin = iju(k)
        jmax = jmin + iu(k+1) - iu(k) - 1
        if (jmin .gt. jmax) go to 5
        do 4 j=jmin,jmax
   4      row(ju(j)) = 0

   5    rk = r(k)
        jmin = ia(rk)
        jmax = ia(rk+1) - 1
        do 6 j=jmin,jmax
          row(ic(ja(j))) = a(j)
   6      continue

        sum = b(rk)
        i = i1
        if (i .eq. 0) go to 10

   7      lki = -row(i)

          l(irl(i)) = -lki
          sum = sum + lki * tmp(i)
          jmin = iu(i)
          jmax = iu(i+1) - 1
          if (jmin .gt. jmax) go to 9
          mu = iju(i) - jmin
          do 8 j=jmin,jmax
   8        row(ju(mu+j)) = row(ju(mu+j)) + lki * u(j)
   9      i = jrl(i)
          if (i .ne. 0) go to 7


  10    if (row(k) .eq. 0.0e0) go to 108
        dk = 1.0e0 / row(k)
        d(k) = dk
        tmp(k) = sum * dk
        if (k .eq. n) go to 19
        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin .gt. jmax)  go to 12
        mu = iju(k) - jmin
        do 11 j=jmin,jmax
  11      u(j) = row(ju(mu+j)) * dk
  12    continue


        i = i1
        if (i .eq. 0) go to 18
  14    irl(i) = irl(i) + 1
        i1 = jrl(i)
        if (irl(i) .ge. il(i+1)) go to 17
        ijlb = irl(i) - il(i) + ijl(i)
        j = jl(ijlb)
  15    if (i .gt. jrl(j)) go to 16
          j = jrl(j)
          go to 15
  16    jrl(i) = jrl(j)
        jrl(j) = i
  17    i = i1
        if (i .ne. 0) go to 14
  18    if (irl(k) .ge. il(k+1)) go to 19
        j = jl(ijl(k))
        jrl(k) = jrl(j)
        jrl(j) = k
  19    continue


      k = n
      do 22 i=1,n
        sum =  tmp(k)
        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin .gt. jmax)  go to 21
        mu = iju(k) - jmin
        do 20 j=jmin,jmax
  20      sum = sum - u(j) * tmp(ju(mu+j))
  21    tmp(k) =  sum
        z(c(k)) =  sum
  22    k = k-1
      flag = 0
      return


 104  flag = 4*n + 1
      return

 107  flag = 7*n + 1
      return

 108  flag = 8*n + k
      return
      end subroutine nnfc
      subroutine nnsc   &
           (n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)



















      integer r(*), c(*), il(*), jl(*), ijl(*), iu(*), ju(*), iju(*)
      real l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk, sum



      do 1 k=1,n
   1    tmp(k) = b(r(k))

      do 3 k=1,n
        jmin = il(k)
        jmax = il(k+1) - 1
        tmpk = -d(k) * tmp(k)
        tmp(k) = -tmpk
        if (jmin .gt. jmax) go to 3
        ml = ijl(k) - jmin
        do 2 j=jmin,jmax
   2      tmp(jl(ml+j)) = tmp(jl(ml+j)) + tmpk * l(j)
   3    continue

      k = n
      do 6 i=1,n
        sum = -tmp(k)
        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin .gt. jmax) go to 5
        mu = iju(k) - jmin
        do 4 j=jmin,jmax
   4      sum = sum + u(j) * tmp(ju(mu+j))
   5    tmp(k) = -sum
        z(c(k)) = -sum
        k = k - 1
   6    continue
      return
      end subroutine nnsc
      subroutine nntc   &
           (n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)




















      integer r(*), c(*), il(*), jl(*), ijl(*), iu(*), ju(*), iju(*)
      real l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum



      do 1 k=1,n
   1    tmp(k) = b(c(k))

      do 3 k=1,n
        jmin = iu(k)
        jmax = iu(k+1) - 1
        tmpk = -tmp(k)
        if (jmin .gt. jmax) go to 3
        mu = iju(k) - jmin
        do 2 j=jmin,jmax
   2      tmp(ju(mu+j)) = tmp(ju(mu+j)) + tmpk * u(j)
   3    continue

      k = n
      do 6 i=1,n
        sum = -tmp(k)
        jmin = il(k)
        jmax = il(k+1) - 1
        if (jmin .gt. jmax) go to 5
        ml = ijl(k) - jmin
        do 4 j=jmin,jmax
   4      sum = sum + l(j) * tmp(jl(ml+j))
   5    tmp(k) = -sum * d(k)
        z(r(k)) = tmp(k)
        k = k - 1
   6    continue
      return
      end subroutine nntc
      subroutine nroc (n, ic, ia, ja, a, jar, ar, p, flag)







































































































































































































      integer  ic(*), ia(*), ja(*), jar(*), p(*), flag
      real  a(*), ar(*)



      do 5 k=1,n
        jmin = ia(k)
        jmax = ia(k+1) - 1
        if(jmin .gt. jmax) go to 5
        p(n+1) = n + 1

        do 3 j=jmin,jmax
          newj = ic(ja(j))
          i = n + 1
   1      if(p(i) .ge. newj) go to 2
            i = p(i)
            go to 1
   2      if(p(i) .eq. newj) go to 102
          p(newj) = p(i)
          p(i) = newj
          jar(newj) = ja(j)
          ar(newj) = a(j)
   3      continue

        i = n + 1
        do 4 j=jmin,jmax
          i = p(i)
          ja(j) = jar(i)
   4      a(j) = ar(i)
   5    continue
      flag = 0
      return


 102  flag = n + k
      return
      end subroutine nroc                                     
      subroutine nsfc   &
            (n, r, ic, ia,ja, jlmax,il,jl,ijl, jumax,iu,ju,iju,   &
             q, ira,jra, irac, irl,jrl, iru,jru, flag)





















































      integer cend, qm, rend, rk, vj



      integer ia(*), ja(*), ira(*), jra(*), il(*), jl(*), ijl(*)
      integer iu(*), ju(*), iju(*), irl(*), jrl(*), iru(*), jru(*)
      integer r(*), ic(*), q(*), irac(*), flag


      np1 = n + 1
      jlmin = 1
      jlptr = 0
      il(1) = 1
      jumin = 1
      juptr = 0
      iu(1) = 1
      do 1 k=1,n
        irac(k) = 0
        jra(k) = 0
        jrl(k) = 0
   1    jru(k) = 0

      do 2 k=1,n
        rk = r(k)
        iak = ia(rk)
        if (iak .ge. ia(rk+1))  go to 101
        jaiak = ic(ja(iak))
        if (jaiak .gt. k)  go to 105
        jra(k) = irac(jaiak)
        irac(jaiak) = k
   2    ira(k) = iak


      do 41 k=1,n


        q(np1) = np1
        luk = -1

        vj = irac(k)
        if (vj .eq. 0)  go to 5
   3      qm = np1
   4      m = qm
          qm =  q(m)
          if (qm .lt. vj)  go to 4
          if (qm .eq. vj)  go to 102
            luk = luk + 1
            q(m) = vj
            q(vj) = qm
            vj = jra(vj)
            if (vj .ne. 0)  go to 3

   5    lastid = 0
        lasti = 0
        ijl(k) = jlptr
        i = k
   6      i = jru(i)
          if (i .eq. 0)  go to 10
          qm = np1
          jmin = irl(i)
          jmax = ijl(i) + il(i+1) - il(i) - 1
          long = jmax - jmin
          if (long .lt. 0)  go to 6
          jtmp = jl(jmin)
          if (jtmp .ne. k)  long = long + 1
          if (jtmp .eq. k)  r(i) = -r(i)
          if (lastid .ge. long)  go to 7
            lasti = i
            lastid = long

   7      do 9 j=jmin,jmax
            vj = jl(j)
   8        m = qm
            qm = q(m)
            if (qm .lt. vj)  go to 8
            if (qm .eq. vj)  go to 9
              luk = luk + 1
              q(m) = vj
              q(vj) = qm
              qm = vj
   9        continue
            go to 6


  10    qm = q(np1)
        if (qm .ne. k)  go to 105
        if (luk .eq. 0)  go to 17
        if (lastid .ne. luk)  go to 11

        irll = irl(lasti)
        ijl(k) = irll + 1
        if (jl(irll) .ne. k)  ijl(k) = ijl(k) - 1
        go to 17

  11    if (jlmin .gt. jlptr)  go to 15
        qm = q(qm)
        do 12 j=jlmin,jlptr
          if (jl(j) - qm)  12, 13, 15
  12      continue
        go to 15
  13    ijl(k) = j
        do 14 i=j,jlptr
          if (jl(i) .ne. qm)  go to 15
          qm = q(qm)
          if (qm .gt. n)  go to 17
  14      continue
        jlptr = j - 1

  15    jlmin = jlptr + 1
        ijl(k) = jlmin
        if (luk .eq. 0)  go to 17
        jlptr = jlptr + luk
        if (jlptr .gt. jlmax)  go to 103
          qm = q(np1)
          do 16 j=jlmin,jlptr
            qm = q(qm)
  16        jl(j) = qm
  17    irl(k) = ijl(k)
        il(k+1) = il(k) + luk


        q(np1) = np1
        luk = -1

        rk = r(k)
        jmin = ira(k)
        jmax = ia(rk+1) - 1
        if (jmin .gt. jmax)  go to 20
        do 19 j=jmin,jmax
          vj = ic(ja(j))
          qm = np1
  18      m = qm
          qm = q(m)
          if (qm .lt. vj)  go to 18
          if (qm .eq. vj)  go to 102
            luk = luk + 1
            q(m) = vj
            q(vj) = qm
  19      continue

  20    lastid = 0
        lasti = 0
        iju(k) = juptr
        i = k
        i1 = jrl(k)
  21      i = i1
          if (i .eq. 0)  go to 26
          i1 = jrl(i)
          qm = np1
          jmin = iru(i)
          jmax = iju(i) + iu(i+1) - iu(i) - 1
          long = jmax - jmin
          if (long .lt. 0)  go to 21
          jtmp = ju(jmin)
          if (jtmp .eq. k)  go to 22

            long = long + 1
            cend = ijl(i) + il(i+1) - il(i)
            irl(i) = irl(i) + 1
            if (irl(i) .ge. cend)  go to 22
              j = jl(irl(i))
              jrl(i) = jrl(j)
              jrl(j) = i
  22      if (lastid .ge. long)  go to 23
            lasti = i
            lastid = long

  23      do 25 j=jmin,jmax
            vj = ju(j)
  24        m = qm
            qm = q(m)
            if (qm .lt. vj)  go to 24
            if (qm .eq. vj)  go to 25
              luk = luk + 1
              q(m) = vj
              q(vj) = qm
              qm = vj
  25        continue
          go to 21

  26    if (il(k+1) .le. il(k))  go to 27
          j = jl(irl(k))
          jrl(k) = jrl(j)
          jrl(j) = k


  27    qm = q(np1)
        if (qm .ne. k)  go to 105
        if (luk .eq. 0)  go to 34
        if (lastid .ne. luk)  go to 28

        irul = iru(lasti)
        iju(k) = irul + 1
        if (ju(irul) .ne. k)  iju(k) = iju(k) - 1
        go to 34

  28    if (jumin .gt. juptr)  go to 32
        qm = q(qm)
        do 29 j=jumin,juptr
          if (ju(j) - qm)  29, 30, 32
  29      continue
        go to 32
  30    iju(k) = j
        do 31 i=j,juptr
          if (ju(i) .ne. qm)  go to 32
          qm = q(qm)
          if (qm .gt. n)  go to 34
  31      continue
        juptr = j - 1

  32    jumin = juptr + 1
        iju(k) = jumin
        if (luk .eq. 0)  go to 34
        juptr = juptr + luk
        if (juptr .gt. jumax)  go to 106
          qm = q(np1)
          do 33 j=jumin,juptr
            qm = q(qm)
  33        ju(j) = qm
  34    iru(k) = iju(k)
        iu(k+1) = iu(k) + luk


        i = k
  35      i1 = jru(i)
          if (r(i) .lt. 0)  go to 36
          rend = iju(i) + iu(i+1) - iu(i)
          if (iru(i) .ge. rend)  go to 37
            j = ju(iru(i))
            jru(i) = jru(j)
            jru(j) = i
            go to 37
  36      r(i) = -r(i)
  37      i = i1
          if (i .eq. 0)  go to 38
          iru(i) = iru(i) + 1
          go to 35


  38    i = irac(k)
        if (i .eq. 0)  go to 41
  39      i1 = jra(i)
          ira(i) = ira(i) + 1
          if (ira(i) .ge. ia(r(i)+1))  go to 40
          irai = ira(i)
          jairai = ic(ja(irai))
          if (jairai .gt. i)  go to 40
          jra(i) = irac(jairai)
          irac(jairai) = i
  40      i = i1
          if (i .ne. 0)  go to 39
  41    continue

      ijl(n) = jlptr
      iju(n) = juptr
      flag = 0
      return


 101  flag = n + rk
      return

 102  flag = 2*n + rk
      return

 103  flag = 3*n + k
      return

 105  flag = 5*n + k
      return

 106  flag = 6*n + k
      return
      end subroutine nsfc
      subroutine odrv   &
           (n, ia,ja,a, p,ip, nsp,isp, path, flag)







































































































































      integer  ia(*), ja(*),  p(*), ip(*),  isp(*),  path,  flag,   &
         v, l, head,  tmp, q
      real  a(*)

      logical  dflag


      flag = 0
      if (path.lt.1 .or. 5.lt.path)  go to 111


      if ((path-1) * (path-2) * (path-4) .ne. 0)  go to 1
        max = (nsp-n)/2
        v    = 1
        l    = v     +  max
        head = l     +  max
        next = head  +  n
        if (max.lt.n)  go to 110

        call  md   &
           (n, ia,ja, max,isp(v),isp(l), isp(head),p,ip, isp(v), flag)
        if (flag.ne.0)  go to 100


   1  if ((path-2) * (path-3) * (path-4) * (path-5) .ne. 0)  go to 2
        tmp = (nsp+1) -      n
        q   = tmp     - (ia(n+1)-1)
        if (q.lt.1)  go to 110

        dflag = path.eq.4 .or. path.eq.5
        call sro   &
           (n,  ip,  ia, ja, a,  isp(tmp),  isp(q),  dflag)

   2  return


 100  return

 110  flag = 10*n + 1
      return

 111  flag = 11*n + 1
      return
      end subroutine odrv



      subroutine prjs (neq,y,yh,nyh,ewt,ftem,savf,wk,iwk,f,jac,   &
                       ruserpar, nruserpar, iuserpar, niuserpar )

      external f,jac
      integer neq, nyh, iwk
      integer iownd, iowns,   &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,   &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,   &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,   &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, imul, j, jj, jok, jmax, jmin, k, kmax, kmin, ng
      integer nruserpar, iuserpar, niuserpar
      real y, yh, ewt, ftem, savf, wk
      real rowns,   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real con0, conmin, ccmxj, psmall, rbig, seth


      real con, di, fac, hl0, pij, r, r0, rcon, rcont,   &
         srur
      real ruserpar


      dimension neq(*), y(*), yh(nyh,*), ewt(*), ftem(*), savf(*),   &
         wk(*), iwk(*)
      dimension ruserpar(nruserpar), iuserpar(niuserpar)
      common /ls0001/ rowns(209),   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,   &
         iownd(14), iowns(6),   &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /lss001/ con0, conmin, ccmxj, psmall, rbig, seth,   &
         iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,   &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,   &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,   &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu






































      hl0 = h*el0
      con = -hl0
      if (miter .eq. 3) go to 300

      jok = 1
      if (nst .eq. 0 .or. nst .ge. nslj+msbj) jok = 0
      if (icf .eq. 1 .and. abs(rc - 1.0e0) .lt. ccmxj) jok = 0
      if (icf .eq. 2) jok = 0
      if (jok .eq. 1) go to 250


 20   jcur = 1
      nje = nje + 1
      nslj = nst
      iplost = 0
      conmin = abs(con)
      go to (100, 200), miter


 100  continue
      kmin = iwk(ipian)
      do 130 j = 1, n
        kmax = iwk(ipian+j) - 1
        do 110 i = 1,n
 110      ftem(i) = 0.0e0
        call jac (neq, tn, y, j, iwk(ipian), iwk(ipjan), ftem,   &
            ruserpar, nruserpar, iuserpar, niuserpar)
        do 120 k = kmin, kmax
          i = iwk(ibjan+k)
          wk(iba+k) = ftem(i)*con
          if (i .eq. j) wk(iba+k) = wk(iba+k) + 1.0e0
 120      continue
        kmin = kmax + 1
 130    continue
      go to 290


 200  continue
      fac = vnorm(n, savf, ewt)
      r0 = 1000.0e0 * abs(h) * uround * float(n) * fac
      if (r0 .eq. 0.0e0) r0 = 1.0e0
      srur = wk(1)
      jmin = iwk(ipigp)
      do 240 ng = 1,ngp
        jmax = iwk(ipigp+ng) - 1
        do 210 j = jmin,jmax
          jj = iwk(ibjgp+j)
          r = amax1(srur*abs(y(jj)),r0/ewt(jj))
 210      y(jj) = y(jj) + r
        call f (neq, tn, y, ftem,   &
            ruserpar, nruserpar, iuserpar, niuserpar)
        do 230 j = jmin,jmax
          jj = iwk(ibjgp+j)
          y(jj) = yh(jj,1)
          r = amax1(srur*abs(y(jj)),r0/ewt(jj))
          fac = -hl0/r
          kmin =iwk(ibian+jj)
          kmax =iwk(ibian+jj+1) - 1
          do 220 k = kmin,kmax
            i = iwk(ibjan+k)
            wk(iba+k) = (ftem(i) - savf(i))*fac
            if (i .eq. jj) wk(iba+k) = wk(iba+k) + 1.0e0
 220        continue
 230      continue
        jmin = jmax + 1
 240    continue
      nfe = nfe + ngp
      go to 290


 250  jcur = 0
      rcon = con/con0
      rcont = abs(con)/conmin
      if (rcont .gt. rbig .and. iplost .eq. 1) go to 20
      kmin = iwk(ipian)
      do 275 j = 1,n
        kmax = iwk(ipian+j) - 1
        do 270 k = kmin,kmax
          i = iwk(ibjan+k)
          pij = wk(iba+k)
          if (i .ne. j) go to 260
          pij = pij - 1.0e0
          if (abs(pij) .ge. psmall) go to 260
            iplost = 1
            conmin = amin1(abs(con0),conmin)
 260      pij = pij*rcon
          if (i .eq. j) pij = pij + 1.0e0
          wk(iba+k) = pij
 270      continue
        kmin = kmax + 1
 275    continue


 290  nlu = nlu + 1
      con0 = con
      ierpj = 0
      do 295 i = 1,n
 295    ftem(i) = 0.0e0
      call cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),   &
         wk(ipa),ftem,ftem,nsp,iwk(ipisp),wk(iprsp),iesp,2,iys)
      if (iys .eq. 0) return
      imul = (iys - 1)/n
      ierpj = -2
      if (imul .eq. 8) ierpj = 1
      if (imul .eq. 10) ierpj = -1
      return


 300  continue
      jcur = 1
      nje = nje + 1
      wk(2) = hl0
      ierpj = 0
      r = el0*0.1e0
      do 310 i = 1,n
 310    y(i) = y(i) + r*(h*savf(i) - yh(i,2))
      call f (neq, tn, y, wk(3),   &
          ruserpar, nruserpar, iuserpar, niuserpar)
      nfe = nfe + 1
      do 320 i = 1,n
        r0 = h*savf(i) - yh(i,2)
        di = 0.1e0*r0 - h*(wk(i+2) - savf(i))
        wk(i+2) = 1.0e0
        if (abs(r0) .lt. uround/ewt(i)) go to 320
        if (abs(di) .eq. 0.0e0) go to 330
        wk(i+2) = 0.1e0*r0/di
 320    continue
      return
 330  ierpj = 2
      return

      end subroutine prjs                                          
      subroutine slss (wk, iwk, x, tem)

      integer iwk
      integer iownd, iowns,   &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,   &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,   &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,   &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i
      real wk, x, tem
      real rowns,   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real rlss
      real di, hl0, phl0, r

      dimension wk(*), iwk(*), x(*), tem(*)

      common /ls0001/ rowns(209),   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,   &
         iownd(14), iowns(6),   &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /lss001/ rlss(6),   &
         iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,   &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,   &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,   &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu


























      iersl = 0
      go to (100, 100, 300), miter
 100  call cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),   &
         wk(ipa),x,x,nsp,iwk(ipisp),wk(iprsp),iesp,4,iersl)
      if (iersl .ne. 0) iersl = -1
      return

 300  phl0 = wk(2)
      hl0 = h*el0
      wk(2) = hl0
      if (hl0 .eq. phl0) go to 330
      r = hl0/phl0
      do 320 i = 1,n
        di = 1.0e0 - r*(1.0e0 - 1.0e0/wk(i+2))
        if (abs(di) .eq. 0.0e0) go to 390
 320    wk(i+2) = 1.0e0/di
 330  do 340 i = 1,n
 340    x(i) = wk(i+2)*x(i)
      return
 390  iersl = 1
      return


      end subroutine slss                  
      subroutine sro   &
           (n, ip, ia,ja,a, q, r, dflag)

































      integer  ip(*),  ia(*), ja(*),  q(*), r(*)
      real  a(*),  ak

      logical  dflag




      do 1 i=1,n
  1     q(i) = 0


      do 3 i=1,n
        jmin = ia(i)
        jmax = ia(i+1) - 1
        if (jmin.gt.jmax)  go to 3
        do 2 j=jmin,jmax


          k = ja(j)
          if (ip(k).lt.ip(i))  ja(j) = i
          if (ip(k).ge.ip(i))  k = i
          r(j) = k


  2       q(k) = q(k) + 1
  3     continue




      do 4 i=1,n
        ia(i+1) = ia(i) + q(i)
  4     q(i) = ia(i+1)



      ilast = 0
      jmin = ia(1)
      jmax = ia(n+1) - 1
      j = jmax
      do 6 jdummy=jmin,jmax
        i = r(j)
        if (.not.dflag .or. ja(j).ne.i .or. i.eq.ilast)  go to 5


          r(j) = ia(i)
          ilast = i
          go to 6


  5       q(i) = q(i) - 1
          r(j) = q(i)

  6     j = j-1



      do 8 j=jmin,jmax
  7     if (r(j).eq.j)  go to 8
          k = r(j)
          r(j) = r(k)
          r(k) = k
          jak = ja(k)
          ja(k) = ja(j)
          ja(j) = jak
          ak = a(k)
          a(k) = a(j)
          a(j) = ak
          go to 7
  8     continue

      return
      end subroutine sro



      real function vnorm (n, v, w)







      integer n,   i
      real v, w,   sum
      dimension v(n), w(n)
      integer iok_vnorm
      common / lsodes_cmn_iok_vnorm / iok_vnorm
      sum = 0.0e0
      do 10 i = 1,n
        if (abs(v(i)*w(i)) .ge. 1.0e18) then
            vnorm = 1.0e18
            iok_vnorm = -1
            return
        end if
 10     sum = sum + (v(i)*w(i))**2
      vnorm = sqrt(sum/float(n))
      return

      end function vnorm          
      subroutine xerrwv (msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)
      use module_peg_util, only:  peg_message, peg_error_fatal

      integer      nmes, nerr, level, ni, i1, i2, nr,   &
         i, lun, lunit, mesflg, ncpw, nch, nwds
      real r1, r2
      character(*) msg
      character*80 errmsg






















































      common /eh0001/ mesflg, lunit


	mesflg = 1
	lunit = 6












      data ncpw/4/




      if (mesflg .eq. 0) go to 100

      lun = lunit

      nch = min0(nmes,60)
      nwds = nch/ncpw
      if (nch .ne. nwds*ncpw) nwds = nwds + 1



      call peg_message( lun, msg )














  10  format(1x,a)



      errmsg = ' '

      if (ni .eq. 1) write (errmsg, 20) i1
 20   format(6x,23hin above message,  i1 =,i10)


      if (ni .eq. 2) write (errmsg, 30) i1,i2
 30   format(6x,23hin above message,  i1 =,i10,3x,4hi2 =,i10)


      if (nr .eq. 1) write (errmsg, 40) r1
 40   format(6x,23hin above message,  r1 =,e21.13)


      if (nr .eq. 2) write (errmsg, 50) r1,r2
 50   format(6x,15hin above,  r1 =,e21.13,3x,4hr2 =,e21.13)

      if (errmsg .ne. ' ') call peg_message( lun, errmsg )


 100  if (level .ne. 2) return
      call peg_error_fatal( lun, '*** subr xerrwv fatal error' )


      end subroutine xerrwv                                                 

      real function r1mach(i)
      use module_peg_util, only:  peg_error_fatal





























      integer mach_small(2)
      integer mach_large(2)
      integer mach_right(2)
      integer mach_diver(2)
      integer mach_log10(2)
      integer sc

      character*80 errmsg

      real rmach(5)

      equivalence (rmach(1), mach_small(1))
      equivalence (rmach(2), mach_large(1))
      equivalence (rmach(3), mach_right(1))
      equivalence (rmach(4), mach_diver(1))
      equivalence (rmach(5), mach_log10(1))
















       data rmach(1) / 1.1754944000E-38 /
       data rmach(2) / 3.4028235000E+38 /
       data rmach(3) / 5.9604645000E-08 /
       data rmach(4) / 1.1920929000E-07 /
       data rmach(5) / 3.0103001000E-01 /
       data sc / 987 /





























































































































































































      real dum




      if (sc .ne. 987) then
          call peg_error_fatal( -1,   &
          '*** func r1mach fatal error -- all data statements inactive' )
      end if

      if (i .lt. 1  .or.  i .gt. 5) goto 999

      r1mach = rmach(i)












































      return



  999 write(errmsg,1999) i
 1999 format('*** func r1mach fatal error -- i out of bounds',i10)
      call peg_error_fatal( -1, errmsg )
      end function r1mach   



      subroutine xsetf (mflag)



      integer mflag, mesflg, lunit
      common /eh0001/ mesflg, lunit

      if (mflag .eq. 0 .or. mflag .eq. 1) mesflg = mflag
      return

      end subroutine xsetf        



      subroutine set_lsodes_common_vars()



      common /eh0001/ mesflg, lunit
      common /ls0001/ rowns(209),   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,   &
         illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,   &
         mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6),   &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu


      illin = 0
      ntrep = 0
      mesflg = 1
      lunit = 6

      return

      end subroutine set_lsodes_common_vars


      end module module_cbmz_lsodes_solver











      subroutine stode_lsodes (neq, y, yh, nyh, yh1, ewt, savf, acor,   &
         wm, iwm, f, jac, pjac, slvs,   &
         ruserpar, nruserpar, iuserpar, niuserpar )
      use module_cbmz_lsodes_solver, only:  cfode, prjs, slss, r1mach, vnorm

      external f, jac, pjac, slvs
      integer neq, nyh, iwm
      integer iownd, ialth, ipup, lmax, meo, nqnyh, nslp,   &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer i, i1, iredo, iret, j, jb, m, ncf, newq
      integer nruserpar, iuserpar, niuserpar
      real y, yh, yh1, ewt, savf, acor, wm
      real conit, crate, el, elco, hold, rmax, tesco,   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround


      real dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup,   &
         r, rh, rhdn, rhsm, rhup, told
      real ruserpar


      dimension neq(*), y(*), yh(nyh,*), yh1(*), ewt(*), savf(*),   &
         acor(*), wm(*), iwm(*)
      dimension ruserpar(nruserpar), iuserpar(niuserpar)
      common /ls0001/ conit, crate, el(13), elco(13,12),   &
         hold, rmax, tesco(3,12),   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, iownd(14),   &
         ialth, ipup, lmax, meo, nqnyh, nslp,   &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu







































































      kflag = 0
      told = tn
      ncf = 0
      ierpj = 0
      iersl = 0
      jcur = 0
      icf = 0
      delp = 0.0e0
      if (jstart .gt. 0) go to 200
      if (jstart .eq. -1) go to 100
      if (jstart .eq. -2) go to 160








      lmax = maxord + 1
      nq = 1
      l = 2
      ialth = 2
      rmax = 10000.0e0
      rc = 0.0e0
      el0 = 1.0e0
      crate = 0.7e0
      hold = h
      meo = meth
      nslp = 0
      ipup = miter
      iret = 3
      go to 140













 100  ipup = miter
      lmax = maxord + 1
      if (ialth .eq. 1) ialth = 2
      if (meth .eq. meo) go to 110
      call cfode (meth, elco, tesco)
      meo = meth
      if (nq .gt. maxord) go to 120
      ialth = l
      iret = 1
      go to 150
 110  if (nq .le. maxord) go to 160
 120  nq = maxord
      l = lmax
      do 125 i = 1,l
 125    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5e0/float(nq+2)
      ddn = vnorm (n, savf, ewt)/tesco(1,l)
      exdn = 1.0e0/float(l)
      rhdn = 1.0e0/(1.3e0*ddn**exdn + 0.0000013e0)
      rh = amin1(rhdn,1.0e0)
      iredo = 3
      if (h .eq. hold) go to 170
      rh = amin1(rh,abs(h/hold))
      h = hold
      go to 175





 140  call cfode (meth, elco, tesco)
 150  do 155 i = 1,l
 155    el(i) = elco(i,nq)
      nqnyh = nq*nyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5e0/float(nq+2)
      go to (160, 170, 200), iret






 160  if (h .eq. hold) go to 200
      rh = h/hold
      h = hold
      iredo = 3
      go to 175
 170  rh = amax1(rh,hmin/abs(h))
 175  rh = amin1(rh,rmax)
      rh = rh/amax1(1.0e0,abs(h)*hmxi*rh)
      r = 1.0e0
      do 180 j = 2,l
        r = r*rh
        do 180 i = 1,n
 180      yh(i,j) = yh(i,j)*r
      h = h*rh
      rc = rc*rh
      ialth = l
      if (iredo .eq. 0) go to 690








 200  if (abs(rc-1.0e0) .gt. ccmax) ipup = miter
      if (nst .ge. nslp+msbp) ipup = miter
      tn = tn + h
      i1 = nqnyh + 1
      do 215 jb = 1,nq
        i1 = i1 - nyh
!dir$ ivdep
        do 210 i = i1,nqnyh
 210      yh1(i) = yh1(i) + yh1(i+nyh)
 215    continue






 220  m = 0
      do 230 i = 1,n
 230    y(i) = yh(i,1)
      call f (neq, tn, y, savf,   &
          ruserpar, nruserpar, iuserpar, niuserpar)
      nfe = nfe + 1
      if (ipup .le. 0) go to 250





      call pjac (neq, y, yh, nyh, ewt, acor, savf, wm, iwm, f, jac,   &
                 ruserpar, nruserpar, iuserpar, niuserpar )
      ipup = 0
      rc = 1.0e0
      nslp = nst
      crate = 0.7e0
      if (ierpj .ne. 0) go to 430
 250  do 260 i = 1,n
 260    acor(i) = 0.0e0
 270  if (miter .ne. 0) go to 350




      do 290 i = 1,n
        savf(i) = h*savf(i) - yh(i,2)
 290    y(i) = savf(i) - acor(i)
      del = vnorm (n, y, ewt)
      do 300 i = 1,n
        y(i) = yh(i,1) + el(1)*savf(i)
 300    acor(i) = savf(i)
      go to 400





 350  do 360 i = 1,n
 360    y(i) = h*savf(i) - (yh(i,2) + acor(i))
      call slvs (wm, iwm, y, savf)
      if (iersl .lt. 0) go to 430
      if (iersl .gt. 0) go to 410
      del = vnorm (n, y, ewt)
      do 380 i = 1,n
        acor(i) = acor(i) + y(i)
 380    y(i) = yh(i,1) + el(1)*acor(i)




 400  if (m .ne. 0) crate = amax1(0.2e0*crate,del/delp)
      dcon = del*amin1(1.0e0,1.5e0*crate)/(tesco(2,nq)*conit)
      if (dcon .le. 1.0e0) go to 450
      m = m + 1
      if (m .eq. maxcor) go to 410
      if (m .ge. 2 .and. del .gt. 2.0e0*delp) go to 410
      delp = del
      call f (neq, tn, y, savf,   &
          ruserpar, nruserpar, iuserpar, niuserpar)
      nfe = nfe + 1
      go to 270







 410  if (miter .eq. 0 .or. jcur .eq. 1) go to 430
      icf = 1
      ipup = miter
      go to 220
 430  icf = 2
      ncf = ncf + 1
      rmax = 2.0e0
      tn = told
      i1 = nqnyh + 1
      do 445 jb = 1,nq
        i1 = i1 - nyh
!dir$ ivdep
        do 440 i = i1,nqnyh
 440      yh1(i) = yh1(i) - yh1(i+nyh)
 445    continue
      if (ierpj .lt. 0 .or. iersl .lt. 0) go to 680
      if (abs(h) .le. hmin*1.00001e0) go to 670
      if (ncf .eq. mxncf) go to 670
      rh = 0.25e0
      ipup = miter
      iredo = 1
      go to 170






 450  jcur = 0
      if (m .eq. 0) dsm = del/tesco(2,nq)
      if (m .gt. 0) dsm = vnorm (n, acor, ewt)/tesco(2,nq)
      if (dsm .gt. 1.0e0) go to 500










      kflag = 0
      iredo = 0
      nst = nst + 1
      hu = h
      nqu = nq
      do 470 j = 1,l
        do 470 i = 1,n
 470      yh(i,j) = yh(i,j) + el(j)*acor(i)
      ialth = ialth - 1
      if (ialth .eq. 0) go to 520
      if (ialth .gt. 1) go to 700
      if (l .eq. lmax) go to 700
      do 490 i = 1,n
 490    yh(i,lmax) = acor(i)
      go to 700







 500  kflag = kflag - 1
      tn = told
      i1 = nqnyh + 1
      do 515 jb = 1,nq
        i1 = i1 - nyh
!dir$ ivdep
        do 510 i = i1,nqnyh
 510      yh1(i) = yh1(i) - yh1(i+nyh)
 515    continue
      rmax = 2.0e0
      if (abs(h) .le. hmin*1.00001e0) go to 660
      if (kflag .le. -3) go to 640
      iredo = 2
      rhup = 0.0e0
      go to 540









 520  rhup = 0.0e0
      if (l .eq. lmax) go to 540
      do 530 i = 1,n
 530    savf(i) = acor(i) - yh(i,lmax)
      dup = vnorm (n, savf, ewt)/tesco(3,nq)
      exup = 1.0e0/float(l+1)
      rhup = 1.0e0/(1.4e0*dup**exup + 0.0000014e0)
 540  exsm = 1.0e0/float(l)
      rhsm = 1.0e0/(1.2e0*dsm**exsm + 0.0000012e0)
      rhdn = 0.0e0
      if (nq .eq. 1) go to 560
      ddn = vnorm (n, yh(1,l), ewt)/tesco(1,nq)
      exdn = 1.0e0/float(nq)
      rhdn = 1.0e0/(1.3e0*ddn**exdn + 0.0000013e0)
 560  if (rhsm .ge. rhup) go to 570
      if (rhup .gt. rhdn) go to 590
      go to 580
 570  if (rhsm .lt. rhdn) go to 580
      newq = nq
      rh = rhsm
      go to 620
 580  newq = nq - 1
      rh = rhdn
      if (kflag .lt. 0 .and. rh .gt. 1.0e0) rh = 1.0e0
      go to 620
 590  newq = l
      rh = rhup
      if (rh .lt. 1.1e0) go to 610
      r = el(l)/float(l)
      do 600 i = 1,n
 600    yh(i,newq+1) = acor(i)*r
      go to 630
 610  ialth = 3
      go to 700
 620  if ((kflag .eq. 0) .and. (rh .lt. 1.1e0)) go to 610
      if (kflag .le. -2) rh = amin1(rh,0.2e0)





      if (newq .eq. nq) go to 170
 630  nq = newq
      l = nq + 1
      iret = 2
      go to 150









 640  if (kflag .eq. -10) go to 660
      rh = 0.1e0
      rh = amax1(hmin/abs(h),rh)
      h = h*rh
      do 645 i = 1,n
 645    y(i) = yh(i,1)
      call f (neq, tn, y, savf,   &
          ruserpar, nruserpar, iuserpar, niuserpar)
      nfe = nfe + 1
      do 650 i = 1,n
 650    yh(i,2) = h*savf(i)
      ipup = miter
      ialth = 5
      if (nq .eq. 1) go to 200
      nq = 1
      l = 2
      iret = 3
      go to 150




 660  kflag = -1
      go to 720
 670  kflag = -2
      go to 720
 680  kflag = -3
      go to 720
 690  rmax = 10.0e0
 700  r = 1.0e0/tesco(2,nqu)
      do 710 i = 1,n
 710    acor(i) = acor(i)*r
 720  hold = h
      jstart = 1
      return

      end subroutine stode_lsodes 



      subroutine prep_lsodes (neq, y, yh, savf, ewt, ftem, ia, ja,   &
                           wk, iwk, ipper, f, jac,   &
                           ruserpar, nruserpar, iuserpar, niuserpar )
      use module_cbmz_lsodes_solver, only:  adjlr, cdrv, cntnzu, jgroup,   &
                                       odrv

      external f,jac
      integer neq, ia, ja, iwk, ipper
      integer iownd, iowns,   &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,   &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,   &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,   &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, ibr, ier, ipil, ipiu, iptt1, iptt2, j, jfound, k,   &
         knew, kmax, kmin, ldif, lenigp, liwk, maxg, np1, nzsut
      integer nruserpar, iuserpar, niuserpar
      real y, yh, savf, ewt, ftem, wk
      real rowns,   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      real con0, conmin, ccmxj, psmall, rbig, seth
      real dq, dyj, erwt, fac, yj
      real ruserpar


      dimension neq(*), y(*), yh(*), savf(*), ewt(*), ftem(*),   &
         ia(*), ja(*), wk(*), iwk(*)
      dimension ruserpar(nruserpar), iuserpar(niuserpar)
      common /ls0001/ rowns(209),   &
         ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,   &
         iownd(14), iowns(6),   &
         icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,   &
         maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu
      common /lss001/ con0, conmin, ccmxj, psmall, rbig, seth,   &
         iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,   &
         ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,   &
         lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,   &
         nslj, ngp, nlu, nnz, nsp, nzl, nzu



































      ibian = lrat*2
      ipian = ibian + 1
      np1 = n + 1
      ipjan = ipian + np1
      ibjan = ipjan - 1
      liwk = lenwk*lrat
      if (ipjan+n-1 .gt. liwk) go to 210
      if (moss .eq. 0) go to 30

      if (istatc .eq. 3) go to 20

      do 10 i = 1,n
        erwt = 1.0e0/ewt(i)
        fac = 1.0e0 + 1.0e0/(float(i)+1.0e0)
        y(i) = y(i) + fac*sign(erwt,y(i))
 10     continue
      go to (70, 100), moss

 20   continue

      do 25 i = 1,n
 25     y(i) = yh(i)
      go to (70, 100), moss


 30   knew = ipjan
      kmin = ia(1)
      iwk(ipian) = 1
      do 60 j = 1,n
        jfound = 0
        kmax = ia(j+1) - 1
        if (kmin .gt. kmax) go to 45
        do 40 k = kmin,kmax
          i = ja(k)
          if (i .eq. j) jfound = 1
          if (knew .gt. liwk) go to 210
          iwk(knew) = i
          knew = knew + 1
 40       continue
        if (jfound .eq. 1) go to 50
 45     if (knew .gt. liwk) go to 210
        iwk(knew) = j
        knew = knew + 1
 50     iwk(ipian+j) = knew + 1 - ipjan
        kmin = kmax + 1
 60     continue
      go to 140


 70   continue

      call f (neq, tn, y, savf,   &
          ruserpar, nruserpar, iuserpar, niuserpar)
      k = ipjan
      iwk(ipian) = 1
      do 90 j = 1,n
        if (k .gt. liwk) go to 210
        iwk(k) = j
        k = k + 1
        do 75 i = 1,n
 75       savf(i) = 0.0e0
        call jac (neq, tn, y, j, iwk(ipian), iwk(ipjan), savf,   &
            ruserpar, nruserpar, iuserpar, niuserpar)
        do 80 i = 1,n
          if (abs(savf(i)) .le. seth) go to 80
          if (i .eq. j) go to 80
          if (k .gt. liwk) go to 210
          iwk(k) = i
          k = k + 1
 80       continue
        iwk(ipian+j) = k + 1 - ipjan
 90     continue
      go to 140


 100  k = ipjan
      iwk(ipian) = 1
      call f (neq, tn, y, savf,   &
          ruserpar, nruserpar, iuserpar, niuserpar)
      do 120 j = 1,n
        if (k .gt. liwk) go to 210
        iwk(k) = j
        k = k + 1
        yj = y(j)
        erwt = 1.0e0/ewt(j)
        dyj = sign(erwt,yj)
        y(j) = yj + dyj
        call f (neq, tn, y, ftem,   &
            ruserpar, nruserpar, iuserpar, niuserpar)
        y(j) = yj
        do 110 i = 1,n
          dq = (ftem(i) - savf(i))/dyj
          if (abs(dq) .le. seth) go to 110
          if (i .eq. j) go to 110
          if (k .gt. liwk) go to 210
          iwk(k) = i
          k = k + 1
 110      continue
        iwk(ipian+j) = k + 1 - ipjan
 120    continue

 140  continue
      if (moss .eq. 0 .or. istatc .ne. 1) go to 150

      do 145 i = 1,n
 145    y(i) = yh(i)
 150  nnz = iwk(ipian+n) - 1
      lenigp = 0
      ipigp = ipjan + nnz
      if (miter .ne. 2) go to 160


      maxg = np1
      ipjgp = ipjan + nnz
      ibjgp = ipjgp - 1
      ipigp = ipjgp + n
      iptt1 = ipigp + np1
      iptt2 = iptt1 + n
      lreq = iptt2 + n - 1
      if (lreq .gt. liwk) go to 220
      call jgroup (n, iwk(ipian), iwk(ipjan), maxg, ngp, iwk(ipigp),   &
         iwk(ipjgp), iwk(iptt1), iwk(iptt2), ier)
      if (ier .ne. 0) go to 220
      lenigp = ngp + 1


 160  ipr = ipigp + lenigp
      ipc = ipr
      ipic = ipc + n
      ipisp = ipic + n
      iprsp = (ipisp - 2)/lrat + 2
      iesp = lenwk + 1 - iprsp
      if (iesp .lt. 0) go to 230
      ibr = ipr - 1
      do 170 i = 1,n
 170    iwk(ibr+i) = i
      nsp = liwk + 1 - ipisp
      call odrv (n, iwk(ipian), iwk(ipjan), wk, iwk(ipr), iwk(ipic),   &
         nsp, iwk(ipisp), 1, iys)
      if (iys .eq. 11*n+1) go to 240
      if (iys .ne. 0) go to 230


      ipa = lenwk + 1 - nnz
      nsp = ipa - iprsp
      lreq = max0(12*n/lrat, 6*n/lrat+2*n+nnz) + 3
      lreq = lreq + iprsp - 1 + nnz
      if (lreq .gt. lenwk) go to 250
      iba = ipa - 1
      do 180 i = 1,nnz
 180    wk(iba+i) = 0.0e0
      ipisp = lrat*(iprsp - 1) + 1
      call cdrv (n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),   &
         wk(ipa),wk(ipa),wk(ipa),nsp,iwk(ipisp),wk(iprsp),iesp,5,iys)
      lreq = lenwk - iesp
      if (iys .eq. 10*n+1) go to 250
      if (iys .ne. 0) go to 260
      ipil = ipisp
      ipiu = ipil + 2*n + 1
      nzu = iwk(ipil+n) - iwk(ipil)
      nzl = iwk(ipiu+n) - iwk(ipiu)
      if (lrat .gt. 1) go to 190
      call adjlr (n, iwk(ipisp), ldif)
      lreq = lreq + ldif
 190  continue
      if (lrat .eq. 2 .and. nnz .eq. n) lreq = lreq + 1
      nsp = nsp + lreq - lenwk
      ipa = lreq + 1 - nnz
      iba = ipa - 1
      ipper = 0
      return

 210  ipper = -1
      lreq = 2 + (2*n + 1)/lrat
      lreq = max0(lenwk+1,lreq)
      return

 220  ipper = -2
      lreq = (lreq - 1)/lrat + 1
      return

 230  ipper = -3
      call cntnzu (n, iwk(ipian), iwk(ipjan), nzsut)
      lreq = lenwk - iesp + (3*n + 4*nzsut - 1)/lrat + 1
      return

 240  ipper = -4
      return

 250  ipper = -5
      return

 260  ipper = -6
      lreq = lenwk
      return

      end subroutine prep_lsodes                                      
