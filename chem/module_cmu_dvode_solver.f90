































      module module_cmu_dvode_solver



      type dvode_cmn_vars
         integer ::   &
            icf, init, ipup, jcur, jstart, jsv, kflag, kuth,   &
            l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm,   &
            locjs, maxord, meth, miter, msbj, mxhnil, mxstep,   &
            n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj,   &
            nslp, nyh,   &
            ncfn, netf, nfe, nje, nlu, nni, nqu, nst
         double precision ::   &
            acnrm, ccmxj, conp, crate, drc, el(13),   &
            eta, etamax, h, hmin, hmxi, hnew, hscal, prl1,   &
            rc, rl1, tau(13), tq(5), tn, uround,   &
            hu,   &
            etaq, etaqm1









      end type dvode_cmn_vars



      contains








      subroutine dvode (f, neq, y, t, tout, itol, rtol, atol, itask,   &
                  istate, iopt, rwork, lrw, iwork, liw, jac, mf,   &
                  rpar, ipar)
      implicit none
      external f, jac
      double precision y, t, tout, rtol, atol, rwork, rpar
      integer neq, itol, itask, istate, iopt, lrw, iwork, liw,   &
              mf, ipar
      dimension y(*), rtol(*), atol(*), rwork(lrw), iwork(liw),   &
                rpar(*), ipar(*)









































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































      logical ihit
      double precision atoli, big, ewti, four, h0, hmax, hmx, hun, one,   &
         pt2, rh, rtoli, size, tcrit, tnext, tolsf, tp, two, zero
      integer i, ier, iflag, imxer, jco, kgo, leniw, lenj, lenp, lenrw,   &
         lenwm, lf0, mband, mfa, ml, mord, mu, mxhnl0, mxstp0, niter,   &
         nslast
      character*80 msg






      dimension mord(2)




      save mord, mxhnl0, mxstp0
      save zero, one, two, four, pt2, hun











































































































      type( dvode_cmn_vars ) :: cmn12
      double precision h
      integer n


      data  mord(1) /12/, mord(2) /5/, mxstp0 /500/, mxhnl0 /10/
      data zero /0.0d0/, one /1.0d0/, two /2.0d0/, four /4.0d0/,   &
           pt2 /0.2d0/, hun /100.0d0/










      if (istate .ne. 1) then
         msg = '*** module_cmu_dvode_solver/dvode -- ' // &
               'fatal error, istate must be 1'
         call wrf_error_fatal3("<stdin>",1264,&
msg )
      end if

      if (istate .lt. 1 .or. istate .gt. 3) go to 601
      if (itask .lt. 1 .or. itask .gt. 5) go to 602
      if (istate .eq. 1) go to 10
      if (cmn12%init .ne. 1) go to 603
      if (istate .eq. 2) go to 200
      go to 20
 10   cmn12%init = 0
      if (tout .eq. t) return









 20   if (neq .le. 0) go to 604
      if (istate .eq. 1) go to 25
      if (neq .gt. n) go to 605
 25   n = neq
      if (itol .lt. 1 .or. itol .gt. 4) go to 606
      if (iopt .lt. 0 .or. iopt .gt. 1) go to 607
      cmn12%jsv = sign(1,mf)
      mfa = abs(mf)
      cmn12%meth = mfa/10
      cmn12%miter = mfa - 10*cmn12%meth
      if (cmn12%meth .lt. 1 .or. cmn12%meth .gt. 2) go to 608
      if (cmn12%miter .lt. 0 .or. cmn12%miter .gt. 5) go to 608
      if (cmn12%miter .le. 3) go to 30
      ml = iwork(1)
      mu = iwork(2)
      if (ml .lt. 0 .or. ml .ge. n) go to 609
      if (mu .lt. 0 .or. mu .ge. n) go to 610
 30   continue

      if (iopt .eq. 1) go to 40
      cmn12%maxord = mord(cmn12%meth)
      cmn12%mxstep = mxstp0
      cmn12%mxhnil = mxhnl0
      if (istate .eq. 1) h0 = zero
      cmn12%hmxi = zero
      cmn12%hmin = zero
      go to 60
 40   cmn12%maxord = iwork(5)
      if (cmn12%maxord .lt. 0) go to 611
      if (cmn12%maxord .eq. 0) cmn12%maxord = 100
      cmn12%maxord = min(cmn12%maxord,mord(cmn12%meth))
      cmn12%mxstep = iwork(6)
      if (cmn12%mxstep .lt. 0) go to 612
      if (cmn12%mxstep .eq. 0) cmn12%mxstep = mxstp0
      cmn12%mxhnil = iwork(7)
      if (cmn12%mxhnil .lt. 0) go to 613
      if (cmn12%mxhnil .eq. 0) cmn12%mxhnil = mxhnl0
      if (istate .ne. 1) go to 50
      h0 = rwork(5)
      if ((tout - t)*h0 .lt. zero) go to 614
 50   hmax = rwork(6)
      if (hmax .lt. zero) go to 615
      cmn12%hmxi = zero
      if (hmax .gt. zero) cmn12%hmxi = one/hmax
      cmn12%hmin = rwork(7)
      if (cmn12%hmin .lt. zero) go to 616







 60   cmn12%lyh = 21
      if (istate .eq. 1) cmn12%nyh = n
      cmn12%lwm = cmn12%lyh + (cmn12%maxord + 1)*cmn12%nyh
      jco = max(0,cmn12%jsv)
      if (cmn12%miter .eq. 0) lenwm = 0
      if (cmn12%miter .eq. 1 .or. cmn12%miter .eq. 2) then
        lenwm = 2 + (1 + jco)*n*n
        cmn12%locjs = n*n + 3
      endif
      if (cmn12%miter .eq. 3) lenwm = 2 + n
      if (cmn12%miter .eq. 4 .or. cmn12%miter .eq. 5) then
        mband = ml + mu + 1
        lenp = (mband + ml)*n
        lenj = mband*n
        lenwm = 2 + lenp + jco*lenj
        cmn12%locjs = lenp + 3
        endif
      cmn12%lewt = cmn12%lwm + lenwm
      cmn12%lsavf = cmn12%lewt + n
      cmn12%lacor = cmn12%lsavf + n
      lenrw = cmn12%lacor + n - 1
      iwork(17) = lenrw
      cmn12%liwm = 1
      leniw = 30 + n
      if (cmn12%miter .eq. 0 .or. cmn12%miter .eq. 3) leniw = 30
      iwork(18) = leniw
      if (lenrw .gt. lrw) go to 617
      if (leniw .gt. liw) go to 618

      rtoli = rtol(1)
      atoli = atol(1)
      do 70 i = 1,n
        if (itol .ge. 3) rtoli = rtol(i)
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        if (rtoli .lt. zero) go to 619
        if (atoli .lt. zero) go to 620
 70     continue
      if (istate .eq. 1) go to 100

      cmn12%jstart = -1
      if (cmn12%nq .le. cmn12%maxord) go to 90

      call dcopy (n, rwork(cmn12%lwm), 1, rwork(cmn12%lsavf), 1)

 90   if (cmn12%miter .gt. 0) rwork(cmn12%lwm) = sqrt(cmn12%uround)
      go to 200







 100  cmn12%uround = dumach()
      cmn12%tn = t
      if (itask .ne. 4 .and. itask .ne. 5) go to 110
      tcrit = rwork(1)
      if ((tcrit - tout)*(tout - t) .lt. zero) go to 625
      if (h0 .ne. zero .and. (t + h0 - tcrit)*h0 .gt. zero)   &
         h0 = tcrit - t
 110  cmn12%jstart = 0
      if (cmn12%miter .gt. 0) rwork(cmn12%lwm) = sqrt(cmn12%uround)
      cmn12%ccmxj = pt2
      cmn12%msbj = 50
      cmn12%nhnil = 0
      cmn12%nst = 0
      cmn12%nje = 0
      cmn12%nni = 0
      cmn12%ncfn = 0
      cmn12%netf = 0
      cmn12%nlu = 0
      cmn12%nslj = 0
      nslast = 0
      cmn12%hu = zero
      cmn12%nqu = 0

      lf0 = cmn12%lyh + cmn12%nyh
      call f (n, t, y, rwork(lf0), rpar, ipar)
      cmn12%nfe = 1

      call dcopy (n, y, 1, rwork(cmn12%lyh), 1)

      cmn12%nq = 1
      h = one
      call dewset (n, itol, rtol, atol, rwork(cmn12%lyh), rwork(cmn12%lewt))
      do 120 i = 1,n
        if (rwork(i+cmn12%lewt-1) .le. zero) go to 621
 120    rwork(i+cmn12%lewt-1) = one/rwork(i+cmn12%lewt-1)
      if (h0 .ne. zero) go to 180

      call dvhin (n, t, rwork(cmn12%lyh), rwork(lf0), f, rpar, ipar, tout,   &
         cmn12%uround, rwork(cmn12%lewt), itol, atol, y, rwork(cmn12%lacor), h0,   &
         niter, ier)
      cmn12%nfe = cmn12%nfe + niter
      if (ier .ne. 0) go to 622

 180  rh = abs(h0)*cmn12%hmxi
      if (rh .gt. one) h0 = h0/rh

      h = h0
      call dscal (n, h0, rwork(lf0), 1)
      go to 270





 200  nslast = cmn12%nst
      cmn12%kuth = 0
      go to (210, 250, 220, 230, 240), itask
 210  if ((cmn12%tn - tout)*h .lt. zero) go to 250
      cmn12%h = h ; cmn12%n = n   
      call dvindy (tout, 0, rwork(cmn12%lyh), cmn12%nyh, y, iflag, cmn12)
      h = cmn12%h ; n = cmn12%n   
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 220  tp = cmn12%tn - cmn12%hu*(one + hun*cmn12%uround)
      if ((tp - tout)*h .gt. zero) go to 623
      if ((cmn12%tn - tout)*h .lt. zero) go to 250
      go to 400
 230  tcrit = rwork(1)
      if ((cmn12%tn - tcrit)*h .gt. zero) go to 624
      if ((tcrit - tout)*h .lt. zero) go to 625
      if ((cmn12%tn - tout)*h .lt. zero) go to 245
      cmn12%h = h ; cmn12%n = n   
      call dvindy (tout, 0, rwork(cmn12%lyh), cmn12%nyh, y, iflag, cmn12)
      h = cmn12%h ; n = cmn12%n   
      if (iflag .ne. 0) go to 627
      t = tout
      go to 420
 240  tcrit = rwork(1)
      if ((cmn12%tn - tcrit)*h .gt. zero) go to 624
 245  hmx = abs(cmn12%tn) + abs(h)
      ihit = abs(cmn12%tn - tcrit) .le. hun*cmn12%uround*hmx
      if (ihit) go to 400
      tnext = cmn12%tn + cmn12%hnew*(one + four*cmn12%uround)
      if ((tnext - tcrit)*h .le. zero) go to 250
      h = (tcrit - cmn12%tn)*(one - four*cmn12%uround)
      cmn12%kuth = 1











 250  continue
      if ((cmn12%nst-nslast) .ge. cmn12%mxstep) go to 500
      call dewset (n, itol, rtol, atol, rwork(cmn12%lyh), rwork(cmn12%lewt))
      do 260 i = 1,n
        if (rwork(i+cmn12%lewt-1) .le. zero) go to 510
 260    rwork(i+cmn12%lewt-1) = one/rwork(i+cmn12%lewt-1)
 270  tolsf = cmn12%uround*dvnorm (n, rwork(cmn12%lyh), rwork(cmn12%lewt))
      if (tolsf .le. one) go to 280
      tolsf = tolsf*two
      if (cmn12%nst .eq. 0) go to 626
      go to 520
 280  if ((cmn12%tn + h) .ne. cmn12%tn) go to 290
      cmn12%nhnil = cmn12%nhnil + 1
      if (cmn12%nhnil .gt. cmn12%mxhnil) go to 290
      msg = 'dvode--  warning: internal t (=r1) and h (=r2) are'
      call xerrwd (msg, 50, 101, 1, 0, 0, 0, 0, zero, zero)
      msg='      such that in the machine, t + h = t on the next step  '
      call xerrwd (msg, 60, 101, 1, 0, 0, 0, 0, zero, zero)
      msg = '      (h = step size). solver will continue anyway'
      call xerrwd (msg, 50, 101, 1, 0, 0, 0, 2, cmn12%tn, h)
      if (cmn12%nhnil .lt. cmn12%mxhnil) go to 290
      msg = 'dvode--  above warning has been issued i1 times.  '
      call xerrwd (msg, 50, 102, 1, 0, 0, 0, 0, zero, zero)
      msg = '      it will not be issued again for this problem'
      call xerrwd (msg, 50, 102, 1, 1, cmn12%mxhnil, 0, 0, zero, zero)
 290  continue




      cmn12%h = h ; cmn12%n = n   
      call dvstep (y, rwork(cmn12%lyh), cmn12%nyh, rwork(cmn12%lyh), rwork(cmn12%lewt),   &
         rwork(cmn12%lsavf), y, rwork(cmn12%lacor), rwork(cmn12%lwm), iwork(cmn12%liwm),   &
         f, jac, f, dvnlsd, rpar, ipar, cmn12)
      h = cmn12%h ; n = cmn12%n   
      kgo = 1 - cmn12%kflag


      go to (300, 530, 540), kgo





 300  cmn12%init = 1
      cmn12%kuth = 0
      go to (310, 400, 330, 340, 350), itask

 310  if ((cmn12%tn - tout)*h .lt. zero) go to 250
      cmn12%h = h ; cmn12%n = n   
      call dvindy (tout, 0, rwork(cmn12%lyh), cmn12%nyh, y, iflag, cmn12)
      h = cmn12%h ; n = cmn12%n   
      t = tout
      go to 420

 330  if ((cmn12%tn - tout)*h .ge. zero) go to 400
      go to 250

 340  if ((cmn12%tn - tout)*h .lt. zero) go to 345
      cmn12%h = h ; cmn12%n = n   
      call dvindy (tout, 0, rwork(cmn12%lyh), cmn12%nyh, y, iflag, cmn12)
      h = cmn12%h ; n = cmn12%n   
      t = tout
      go to 420
 345  hmx = abs(cmn12%tn) + abs(h)
      ihit = abs(cmn12%tn - tcrit) .le. hun*cmn12%uround*hmx
      if (ihit) go to 400
      tnext = cmn12%tn + cmn12%hnew*(one + four*cmn12%uround)
      if ((tnext - tcrit)*h .le. zero) go to 250
      h = (tcrit - cmn12%tn)*(one - four*cmn12%uround)
      cmn12%kuth = 1
      go to 250

 350  hmx = abs(cmn12%tn) + abs(h)
      ihit = abs(cmn12%tn - tcrit) .le. hun*cmn12%uround*hmx







 400  continue
      call dcopy (n, rwork(cmn12%lyh), 1, y, 1)
      t = cmn12%tn
      if (itask .ne. 4 .and. itask .ne. 5) go to 420
      if (ihit) t = tcrit
 420  istate = 2
      rwork(11) = cmn12%hu
      rwork(12) = cmn12%hnew
      rwork(13) = cmn12%tn
      iwork(11) = cmn12%nst
      iwork(12) = cmn12%nfe
      iwork(13) = cmn12%nje
      iwork(14) = cmn12%nqu
      iwork(15) = cmn12%newq
      iwork(19) = cmn12%nlu
      iwork(20) = cmn12%nni
      iwork(21) = cmn12%ncfn
      iwork(22) = cmn12%netf
      return









 500  msg = 'dvode--  at current t (=r1), mxstep (=i1) steps   '
      call xerrwd (msg, 50, 201, 1, 0, 0, 0, 0, zero, zero)
      msg = '      taken on this call before reaching tout     '
      call xerrwd (msg, 50, 201, 1, 1, cmn12%mxstep, 0, 1, cmn12%tn, zero)
      istate = -1
      go to 580

 510  ewti = rwork(cmn12%lewt+i-1)
      msg = 'dvode--  at t (=r1), ewt(i1) has become r2 .le. 0.'
      call xerrwd (msg, 50, 202, 1, 1, i, 0, 2, cmn12%tn, ewti)
      istate = -6
      go to 580

 520  msg = 'dvode--  at t (=r1), too much accuracy requested  '
      call xerrwd (msg, 50, 203, 1, 0, 0, 0, 0, zero, zero)
      msg = '      for precision of machine:   see tolsf (=r2) '
      call xerrwd (msg, 50, 203, 1, 0, 0, 0, 2, cmn12%tn, tolsf)
      rwork(14) = tolsf
      istate = -2
      go to 580

 530  msg = 'dvode--  at t(=r1) and step size h(=r2), the error'
      call xerrwd (msg, 50, 204, 1, 0, 0, 0, 0, zero, zero)
      msg = '      test failed repeatedly or with abs(h) = hmin'
      call xerrwd (msg, 50, 204, 1, 0, 0, 0, 2, cmn12%tn, h)
      istate = -4
      go to 560

 540  msg = 'dvode--  at t (=r1) and step size h (=r2), the    '
      call xerrwd (msg, 50, 205, 1, 0, 0, 0, 0, zero, zero)
      msg = '      corrector convergence failed repeatedly     '
      call xerrwd (msg, 50, 205, 1, 0, 0, 0, 0, zero, zero)
      msg = '      or with abs(h) = hmin   '
      call xerrwd (msg, 30, 205, 1, 0, 0, 0, 2, cmn12%tn, h)
      istate = -5

 560  big = zero
      imxer = 1
      do 570 i = 1,n
        size = abs(rwork(i+cmn12%lacor-1)*rwork(i+cmn12%lewt-1))
        if (big .ge. size) go to 570
        big = size
        imxer = i
 570    continue
      iwork(16) = imxer

 580  continue
      call dcopy (n, rwork(cmn12%lyh), 1, y, 1)
      t = cmn12%tn
      rwork(11) = cmn12%hu
      rwork(12) = h
      rwork(13) = cmn12%tn
      iwork(11) = cmn12%nst
      iwork(12) = cmn12%nfe
      iwork(13) = cmn12%nje
      iwork(14) = cmn12%nqu
      iwork(15) = cmn12%nq
      iwork(19) = cmn12%nlu
      iwork(20) = cmn12%nni
      iwork(21) = cmn12%ncfn
      iwork(22) = cmn12%netf
      return







 601  msg = 'dvode--  istate (=i1) illegal '
      call xerrwd (msg, 30, 1, 1, 1, istate, 0, 0, zero, zero)
      if (istate .lt. 0) go to 800
      go to 700
 602  msg = 'dvode--  itask (=i1) illegal  '
      call xerrwd (msg, 30, 2, 1, 1, itask, 0, 0, zero, zero)
      go to 700
 603  msg='dvode--  istate (=i1) .gt. 1 but dvode not initialized      '
      call xerrwd (msg, 60, 3, 1, 1, istate, 0, 0, zero, zero)
      go to 700
 604  msg = 'dvode--  neq (=i1) .lt. 1     '
      call xerrwd (msg, 30, 4, 1, 1, neq, 0, 0, zero, zero)
      go to 700
 605  msg = 'dvode--  istate = 3 and neq increased (i1 to i2)  '
      call xerrwd (msg, 50, 5, 1, 2, n, neq, 0, zero, zero)
      go to 700
 606  msg = 'dvode--  itol (=i1) illegal   '
      call xerrwd (msg, 30, 6, 1, 1, itol, 0, 0, zero, zero)
      go to 700
 607  msg = 'dvode--  iopt (=i1) illegal   '
      call xerrwd (msg, 30, 7, 1, 1, iopt, 0, 0, zero, zero)
      go to 700
 608  msg = 'dvode--  mf (=i1) illegal     '
      call xerrwd (msg, 30, 8, 1, 1, mf, 0, 0, zero, zero)
      go to 700
 609  msg = 'dvode--  ml (=i1) illegal:  .lt.0 or .ge.neq (=i2)'
      call xerrwd (msg, 50, 9, 1, 2, ml, neq, 0, zero, zero)
      go to 700
 610  msg = 'dvode--  mu (=i1) illegal:  .lt.0 or .ge.neq (=i2)'
      call xerrwd (msg, 50, 10, 1, 2, mu, neq, 0, zero, zero)
      go to 700
 611  msg = 'dvode--  maxord (=i1) .lt. 0  '
      call xerrwd (msg, 30, 11, 1, 1, cmn12%maxord, 0, 0, zero, zero)
      go to 700
 612  msg = 'dvode--  mxstep (=i1) .lt. 0  '
      call xerrwd (msg, 30, 12, 1, 1, cmn12%mxstep, 0, 0, zero, zero)
      go to 700
 613  msg = 'dvode--  mxhnil (=i1) .lt. 0  '
      call xerrwd (msg, 30, 13, 1, 1, cmn12%mxhnil, 0, 0, zero, zero)
      go to 700
 614  msg = 'dvode--  tout (=r1) behind t (=r2)      '
      call xerrwd (msg, 40, 14, 1, 0, 0, 0, 2, tout, t)
      msg = '      integration direction is given by h0 (=r1)  '
      call xerrwd (msg, 50, 14, 1, 0, 0, 0, 1, h0, zero)
      go to 700
 615  msg = 'dvode--  hmax (=r1) .lt. 0.0  '
      call xerrwd (msg, 30, 15, 1, 0, 0, 0, 1, hmax, zero)
      go to 700
 616  msg = 'dvode--  hmin (=r1) .lt. 0.0  '
      call xerrwd (msg, 30, 16, 1, 0, 0, 0, 1, cmn12%hmin, zero)
      go to 700
 617  continue
      msg='dvode--  rwork length needed, lenrw (=i1), exceeds lrw (=i2)'
      call xerrwd (msg, 60, 17, 1, 2, lenrw, lrw, 0, zero, zero)
      go to 700
 618  continue
      msg='dvode--  iwork length needed, leniw (=i1), exceeds liw (=i2)'
      call xerrwd (msg, 60, 18, 1, 2, leniw, liw, 0, zero, zero)
      go to 700
 619  msg = 'dvode--  rtol(i1) is r1 .lt. 0.0        '
      call xerrwd (msg, 40, 19, 1, 1, i, 0, 1, rtoli, zero)
      go to 700
 620  msg = 'dvode--  atol(i1) is r1 .lt. 0.0        '
      call xerrwd (msg, 40, 20, 1, 1, i, 0, 1, atoli, zero)
      go to 700
 621  ewti = rwork(cmn12%lewt+i-1)
      msg = 'dvode--  ewt(i1) is r1 .le. 0.0         '
      call xerrwd (msg, 40, 21, 1, 1, i, 0, 1, ewti, zero)
      go to 700
 622  continue
      msg='dvode--  tout (=r1) too close to t(=r2) to start integration'
      call xerrwd (msg, 60, 22, 1, 0, 0, 0, 2, tout, t)
      go to 700
 623  continue
      msg='dvode--  itask = i1 and tout (=r1) behind tcur - hu (= r2)  '
      call xerrwd (msg, 60, 23, 1, 1, itask, 0, 2, tout, tp)
      go to 700
 624  continue
      msg='dvode--  itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   '
      call xerrwd (msg, 60, 24, 1, 0, 0, 0, 2, tcrit, cmn12%tn)
      go to 700
 625  continue
      msg='dvode--  itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   '
      call xerrwd (msg, 60, 25, 1, 0, 0, 0, 2, tcrit, tout)
      go to 700
 626  msg = 'dvode--  at start of problem, too much accuracy   '
      call xerrwd (msg, 50, 26, 1, 0, 0, 0, 0, zero, zero)
      msg='      requested for precision of machine:   see tolsf (=r1) '
      call xerrwd (msg, 60, 26, 1, 0, 0, 0, 1, tolsf, zero)
      rwork(14) = tolsf
      go to 700
 627  msg='dvode--  trouble from dvindy.  itask = i1, tout = r1.       '
      call xerrwd (msg, 60, 27, 1, 1, itask, 0, 1, tout, zero)

 700  continue
      istate = -3
      return

 800  msg = 'dvode--  run aborted:  apparent infinite loop     '
      call xerrwd (msg, 50, 303, 2, 0, 0, 0, 0, zero, zero)
      return

      end subroutine dvode

      subroutine dvhin (n, t0, y0, ydot, f, rpar, ipar, tout, uround,   &
         ewt, itol, atol, y, temp, h0, niter, ier)
      external f
      double precision t0, y0, ydot, rpar, tout, uround, ewt, atol, y,   &
         temp, h0
      integer n, ipar, itol, niter, ier
      dimension y0(*), ydot(*), ewt(*), atol(*), y(*),   &
         temp(*), rpar(*), ipar(*)









































      double precision afi, atoli, delyi, h, half, hg, hlb, hnew, hrat,   &
           hub, hun, pt1, t1, tdist, tround, two, yddnrm
      integer i, iter









      save half, hun, pt1, two
      data half /0.5d0/, hun /100.0d0/, pt1 /0.1d0/, two /2.0d0/

      niter = 0
      tdist = abs(tout - t0)
      tround = uround*max(abs(t0),abs(tout))
      if (tdist .lt. two*tround) go to 100


      hlb = hun*tround

      hub = pt1*tdist
      atoli = atol(1)
      do 10 i = 1, n
        if (itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
        delyi = pt1*abs(y0(i)) + atoli
        afi = abs(ydot(i))
        if (afi*hub .gt. delyi) hub = delyi/afi
 10     continue


      iter = 0
      hg = sqrt(hlb*hub)

      if (hub .lt. hlb) then
        h0 = hg
        go to 90
      endif


 50   continue

      h = sign (hg, tout - t0)
      t1 = t0 + h
      do 60 i = 1, n
 60     y(i) = y0(i) + h*ydot(i)
      call f (n, t1, y, temp, rpar, ipar)
      do 70 i = 1, n
 70     temp(i) = (temp(i) - ydot(i))/h
      yddnrm = dvnorm (n, temp, ewt)

      if (yddnrm*hub*hub .gt. two) then
        hnew = sqrt(two/yddnrm)
      else
        hnew = sqrt(hg*hub)
      endif
      iter = iter + 1







      if (iter .ge. 4) go to 80
      hrat = hnew/hg
      if ( (hrat .gt. half) .and. (hrat .lt. two) ) go to 80
      if ( (iter .ge. 2) .and. (hnew .gt. two*hg) ) then
        hnew = hg
        go to 80
      endif
      hg = hnew
      go to 50


 80   h0 = hnew*half
      if (h0 .lt. hlb) h0 = hlb
      if (h0 .gt. hub) h0 = hub
 90   h0 = sign(h0, tout - t0)
      niter = iter
      ier = 0
      return

 100  ier = -1
      return

      end subroutine dvhin

      subroutine dvindy (t, k, yh, ldyh, dky, iflag, cmn12)
      implicit none
      double precision t, yh, dky
      integer k, ldyh, iflag
      dimension yh(ldyh,*), dky(*)
      type( dvode_cmn_vars ) :: cmn12


















































      double precision c, hun, r, s, tfuzz, tn1, tp, zero
      integer i, ic, j, jb, jb2, jj, jj1, jp1
      character*80 msg




      save hun, zero











      data hun /100.0d0/, zero /0.0d0/

      iflag = 0
      if (k .lt. 0 .or. k .gt. cmn12%nq) go to 80
      tfuzz = hun*cmn12%uround*(cmn12%tn + cmn12%hu)
      tp = cmn12%tn - cmn12%hu - tfuzz
      tn1 = cmn12%tn + tfuzz
      if ((t-tp)*(t-tn1) .gt. zero) go to 90

      s = (t - cmn12%tn)/cmn12%h
      ic = 1
      if (k .eq. 0) go to 15
      jj1 = cmn12%l - k
      do 10 jj = jj1, cmn12%nq
 10     ic = ic*jj
 15   c = real(ic)
      do 20 i = 1, cmn12%n
 20     dky(i) = c*yh(i,cmn12%l)
      if (k .eq. cmn12%nq) go to 55
      jb2 = cmn12%nq - k
      do 50 jb = 1, jb2
        j = cmn12%nq - jb
        jp1 = j + 1
        ic = 1
        if (k .eq. 0) go to 35
        jj1 = jp1 - k
        do 30 jj = jj1, j
 30       ic = ic*jj
 35     c = real(ic)
        do 40 i = 1, cmn12%n
 40       dky(i) = c*yh(i,jp1) + s*dky(i)
 50     continue
      if (k .eq. 0) return
 55   r = cmn12%h**(-k)
      call dscal (cmn12%n, r, dky, 1)
      return

 80   msg = 'dvindy-- k (=i1) illegal      '
      call xerrwd (msg, 30, 51, 1, 1, k, 0, 0, zero, zero)
      iflag = -1
      return
 90   msg = 'dvindy-- t (=r1) illegal      '
      call xerrwd (msg, 30, 52, 1, 0, 0, 0, 1, t, zero)
      msg='      t not in interval tcur - hu (= r1) to tcur (=r2)      '
      call xerrwd (msg, 60, 52, 1, 0, 0, 0, 2, tp, cmn12%tn)
      iflag = -2
      return

      end subroutine dvindy

      subroutine dvstep (y, yh, ldyh, yh1, ewt, savf, vsav, acor,   &
                        wm, iwm, f, jac, psol, vnls, rpar, ipar, cmn12)
      implicit none
      external f, jac, psol, vnls
      double precision y, yh, yh1, ewt, savf, vsav, acor, wm, rpar
      integer ldyh, iwm, ipar
      dimension y(*), yh(ldyh,*), yh1(*), ewt(*), savf(*), vsav(*),   &
         acor(*), wm(*), iwm(*), rpar(*), ipar(*)
      type( dvode_cmn_vars ) :: cmn12


















































































      double precision addon, bias1,bias2,bias3, cnquot, ddn, dsm, dup,   &
           etacf, etamin, etamx1, etamx2, etamx3, etamxf,   &
           etaqp1, flotl, one, onepsm,   &
           r, thresh, told, zero
      integer i, i1, i2, iback, j, jb, kfc, kfh, mxncf, ncf, nflag














      save addon, bias1, bias2, bias3,   &
           etacf, etamin, etamx1, etamx2, etamx3, etamxf,   &
           kfc, kfh, mxncf, onepsm, thresh, one, zero











      data kfc/-3/, kfh/-7/, mxncf/10/
      data addon  /1.0d-6/,    bias1  /6.0d0/,     bias2  /6.0d0/,   &
           bias3  /10.0d0/,    etacf  /0.25d0/,    etamin /0.1d0/,   &
           etamxf /0.2d0/,     etamx1 /1.0d4/,     etamx2 /10.0d0/,   &
           etamx3 /10.0d0/,    onepsm /1.00001d0/, thresh /1.5d0/
      data one/1.0d0/, zero/0.0d0/

      cmn12%kflag = 0
      told = cmn12%tn
      ncf = 0
      cmn12%jcur = 0
      nflag = 0
      if (cmn12%jstart .gt. 0) go to 20
      if (cmn12%jstart .eq. -1) go to 100








      cmn12%lmax = cmn12%maxord + 1
      cmn12%nq = 1
      cmn12%l = 2
      cmn12%nqnyh = cmn12%nq*ldyh
      cmn12%tau(1) = cmn12%h
      cmn12%prl1 = one
      cmn12%rc = zero
      cmn12%etamax = etamx1
      cmn12%nqwait = 2
      cmn12%hscal = cmn12%h
      go to 200









 20   continue
      if (cmn12%kuth .eq. 1) then
        cmn12%eta = min(cmn12%eta,cmn12%h/cmn12%hscal)
        cmn12%newh = 1
        endif
 50   if (cmn12%newh .eq. 0) go to 200
      if (cmn12%newq .eq. cmn12%nq) go to 150
      if (cmn12%newq .lt. cmn12%nq) then
        call dvjust (yh, ldyh, -1, cmn12)
        cmn12%nq = cmn12%newq
        cmn12%l = cmn12%nq + 1
        cmn12%nqwait = cmn12%l
        go to 150
        endif
      if (cmn12%newq .gt. cmn12%nq) then
        call dvjust (yh, ldyh, 1, cmn12)
        cmn12%nq = cmn12%newq
        cmn12%l = cmn12%nq + 1
        cmn12%nqwait = cmn12%l
        go to 150
      endif












 100  continue
      cmn12%lmax = cmn12%maxord + 1
      if (cmn12%n .eq. ldyh) go to 120
      i1 = 1 + (cmn12%newq + 1)*ldyh
      i2 = (cmn12%maxord + 1)*ldyh
      if (i1 .gt. i2) go to 120
      do 110 i = i1, i2
 110    yh1(i) = zero
 120  if (cmn12%newq .le. cmn12%maxord) go to 140
      flotl = real(cmn12%lmax)
      if (cmn12%maxord .lt. cmn12%nq-1) then
        ddn = dvnorm (cmn12%n, savf, ewt)/cmn12%tq(1)
        cmn12%eta = one/((bias1*ddn)**(one/flotl) + addon)
        endif
      if (cmn12%maxord .eq. cmn12%nq .and. cmn12%newq .eq. cmn12%nq+1) cmn12%eta = cmn12%etaq
      if (cmn12%maxord .eq. cmn12%nq-1 .and. cmn12%newq .eq. cmn12%nq+1) then
        cmn12%eta = cmn12%etaqm1
        call dvjust (yh, ldyh, -1, cmn12)
        endif
      if (cmn12%maxord .eq. cmn12%nq-1 .and. cmn12%newq .eq. cmn12%nq) then
        ddn = dvnorm (cmn12%n, savf, ewt)/cmn12%tq(1)
        cmn12%eta = one/((bias1*ddn)**(one/flotl) + addon)
        call dvjust (yh, ldyh, -1, cmn12)
        endif
      cmn12%eta = min(cmn12%eta,one)
      cmn12%nq = cmn12%maxord
      cmn12%l = cmn12%lmax
 140  if (cmn12%kuth .eq. 1) cmn12%eta = min(cmn12%eta,abs(cmn12%h/cmn12%hscal))
      if (cmn12%kuth .eq. 0) cmn12%eta = max(cmn12%eta,cmn12%hmin/abs(cmn12%hscal))
      cmn12%eta = cmn12%eta/max(one,abs(cmn12%hscal)*cmn12%hmxi*cmn12%eta)
      cmn12%newh = 1
      cmn12%nqwait = cmn12%l
      if (cmn12%newq .le. cmn12%maxord) go to 50

 150  r = one
      do 180 j = 2, cmn12%l
        r = r*cmn12%eta
        call dscal (cmn12%n, r, yh(1,j), 1 )
 180    continue
      cmn12%h = cmn12%hscal*cmn12%eta
      cmn12%hscal = cmn12%h
      cmn12%rc = cmn12%rc*cmn12%eta
      cmn12%nqnyh = cmn12%nq*ldyh






 200  cmn12%tn = cmn12%tn + cmn12%h
      i1 = cmn12%nqnyh + 1
      do 220 jb = 1, cmn12%nq
        i1 = i1 - ldyh
        do 210 i = i1, cmn12%nqnyh
 210      yh1(i) = yh1(i) + yh1(i+ldyh)
 220  continue
      call dvset( cmn12 )
      cmn12%rl1 = one/cmn12%el(2)
      cmn12%rc = cmn12%rc*(cmn12%rl1/cmn12%prl1)
      cmn12%prl1 = cmn12%rl1



      call vnls (y, yh, ldyh, vsav, savf, ewt, acor, iwm, wm,   &
                 f, jac, psol, nflag, rpar, ipar, cmn12)

      if (nflag .eq. 0) go to 450






        ncf = ncf + 1
        cmn12%ncfn = cmn12%ncfn + 1
        cmn12%etamax = one
        cmn12%tn = told
        i1 = cmn12%nqnyh + 1
        do 430 jb = 1, cmn12%nq
          i1 = i1 - ldyh
          do 420 i = i1, cmn12%nqnyh
 420        yh1(i) = yh1(i) - yh1(i+ldyh)
 430      continue
        if (nflag .lt. -1) go to 680
        if (abs(cmn12%h) .le. cmn12%hmin*onepsm) go to 670
        if (ncf .eq. mxncf) go to 670
        cmn12%eta = etacf
        cmn12%eta = max(cmn12%eta,cmn12%hmin/abs(cmn12%h))
        nflag = -1
        go to 150




 450  continue
      dsm = cmn12%acnrm/cmn12%tq(2)
      if (dsm .gt. one) go to 500






      cmn12%kflag = 0
      cmn12%nst = cmn12%nst + 1
      cmn12%hu = cmn12%h
      cmn12%nqu = cmn12%nq
      do 470 iback = 1, cmn12%nq
        i = cmn12%l - iback
 470    cmn12%tau(i+1) = cmn12%tau(i)
      cmn12%tau(1) = cmn12%h
      do 480 j = 1, cmn12%l
        call daxpy (cmn12%n, cmn12%el(j), acor, 1, yh(1,j), 1 )
 480    continue
      cmn12%nqwait = cmn12%nqwait - 1
      if ((cmn12%l .eq. cmn12%lmax) .or. (cmn12%nqwait .ne. 1)) go to 490
      call dcopy (cmn12%n, acor, 1, yh(1,cmn12%lmax), 1 )
      cmn12%conp = cmn12%tq(5)
 490  if (cmn12%etamax .ne. one) go to 560
      if (cmn12%nqwait .lt. 2) cmn12%nqwait = 2
      cmn12%newq = cmn12%nq
      cmn12%newh = 0
      cmn12%eta = one
      cmn12%hnew = cmn12%h
      go to 690







 500  cmn12%kflag = cmn12%kflag - 1
      cmn12%netf = cmn12%netf + 1
      nflag = -2
      cmn12%tn = told
      i1 = cmn12%nqnyh + 1
      do 520 jb = 1, cmn12%nq
        i1 = i1 - ldyh
        do 510 i = i1, cmn12%nqnyh
 510      yh1(i) = yh1(i) - yh1(i+ldyh)
 520  continue
      if (abs(cmn12%h) .le. cmn12%hmin*onepsm) go to 660
      cmn12%etamax = one
      if (cmn12%kflag .le. kfc) go to 530

      flotl = real(cmn12%l)
      cmn12%eta = one/((bias2*dsm)**(one/flotl) + addon)
      cmn12%eta = max(cmn12%eta,cmn12%hmin/abs(cmn12%h),etamin)
      if ((cmn12%kflag .le. -2) .and. (cmn12%eta .gt. etamxf)) cmn12%eta = etamxf
      go to 150








 530  if (cmn12%kflag .eq. kfh) go to 660
      if (cmn12%nq .eq. 1) go to 540
      cmn12%eta = max(etamin,cmn12%hmin/abs(cmn12%h))
      call dvjust (yh, ldyh, -1, cmn12)
      cmn12%l = cmn12%nq
      cmn12%nq = cmn12%nq - 1
      cmn12%nqwait = cmn12%l
      go to 150
 540  cmn12%eta = max(etamin,cmn12%hmin/abs(cmn12%h))
      cmn12%h = cmn12%h*cmn12%eta
      cmn12%hscal = cmn12%h
      cmn12%tau(1) = cmn12%h
      call f (cmn12%n, cmn12%tn, y, savf, rpar, ipar)
      cmn12%nfe = cmn12%nfe + 1
      do 550 i = 1, cmn12%n
 550    yh(i,2) = cmn12%h*savf(i)
      cmn12%nqwait = 10
      go to 200











 560  flotl = real(cmn12%l)
      cmn12%etaq = one/((bias2*dsm)**(one/flotl) + addon)
      if (cmn12%nqwait .ne. 0) go to 600
      cmn12%nqwait = 2
      cmn12%etaqm1 = zero
      if (cmn12%nq .eq. 1) go to 570

      ddn = dvnorm (cmn12%n, yh(1,cmn12%l), ewt)/cmn12%tq(1)
      cmn12%etaqm1 = one/((bias1*ddn)**(one/(flotl - one)) + addon)
 570  etaqp1 = zero
      if (cmn12%l .eq. cmn12%lmax) go to 580

      cnquot = (cmn12%tq(5)/cmn12%conp)*(cmn12%h/cmn12%tau(2))**cmn12%l
      do 575 i = 1, cmn12%n
 575    savf(i) = acor(i) - cnquot*yh(i,cmn12%lmax)
      dup = dvnorm (cmn12%n, savf, ewt)/cmn12%tq(3)
      etaqp1 = one/((bias3*dup)**(one/(flotl + one)) + addon)
 580  if (cmn12%etaq .ge. etaqp1) go to 590
      if (etaqp1 .gt. cmn12%etaqm1) go to 620
      go to 610
 590  if (cmn12%etaq .lt. cmn12%etaqm1) go to 610
 600  cmn12%eta = cmn12%etaq
      cmn12%newq = cmn12%nq
      go to 630
 610  cmn12%eta = cmn12%etaqm1
      cmn12%newq = cmn12%nq - 1
      go to 630
 620  cmn12%eta = etaqp1
      cmn12%newq = cmn12%nq + 1
      call dcopy (cmn12%n, acor, 1, yh(1,cmn12%lmax), 1)

 630  if (cmn12%eta .lt. thresh .or. cmn12%etamax .eq. one) go to 640
      cmn12%eta = min(cmn12%eta,cmn12%etamax)
      cmn12%eta = cmn12%eta/max(one,abs(cmn12%h)*cmn12%hmxi*cmn12%eta)
      cmn12%newh = 1
      cmn12%hnew = cmn12%h*cmn12%eta
      go to 690
 640  cmn12%newq = cmn12%nq
      cmn12%newh = 0
      cmn12%eta = one
      cmn12%hnew = cmn12%h
      go to 690




 660  cmn12%kflag = -1
      go to 720
 670  cmn12%kflag = -2
      go to 720
 680  if (nflag .eq. -2) cmn12%kflag = -3
      if (nflag .eq. -3) cmn12%kflag = -4
      go to 720
 690  cmn12%etamax = etamx3
      if (cmn12%nst .le. 10) cmn12%etamax = etamx2
 700  r = one/cmn12%tq(2)
      call dscal (cmn12%n, r, acor, 1)
 720  cmn12%jstart = 1
      return

      end subroutine dvstep

      subroutine dvset( cmn12 )
      implicit none
      type( dvode_cmn_vars ) :: cmn12



























































      double precision ahatn0, alph0, cnqm1, cortes, csum, elp, em,   &
           em0, floti, flotl, flotnq, hsum, one, rxi, rxis, s, six,   &
           t1, t2, t3, t4, t5, t6, two, xi, zero
      integer i, iback, j, jp1, nqm1, nqm2

      dimension em(13)




      save cortes, one, six, two, zero










      data cortes /0.1d0/
      data one  /1.0d0/, six /6.0d0/, two /2.0d0/, zero /0.0d0/

      flotl = real(cmn12%l)
      nqm1 = cmn12%nq - 1
      nqm2 = cmn12%nq - 2
      go to (100, 200), cmn12%meth


 100  if (cmn12%nq .ne. 1) go to 110
      cmn12%el(1) = one
      cmn12%el(2) = one
      cmn12%tq(1) = one
      cmn12%tq(2) = two
      cmn12%tq(3) = six*cmn12%tq(2)
      cmn12%tq(5) = one
      go to 300
 110  hsum = cmn12%h
      em(1) = one
      flotnq = flotl - one
      do 115 i = 2, cmn12%l
 115    em(i) = zero
      do 150 j = 1, nqm1
        if ((j .ne. nqm1) .or. (cmn12%nqwait .ne. 1)) go to 130
        s = one
        csum = zero
        do 120 i = 1, nqm1
          csum = csum + s*em(i)/real(i+1)
 120      s = -s
        cmn12%tq(1) = em(nqm1)/(flotnq*csum)
 130    rxi = cmn12%h/hsum
        do 140 iback = 1, j
          i = (j + 2) - iback
 140      em(i) = em(i) + em(i-1)*rxi
        hsum = hsum + cmn12%tau(j)
 150    continue

      s = one
      em0 = zero
      csum = zero
      do 160 i = 1, cmn12%nq
        floti = real(i)
        em0 = em0 + s*em(i)/floti
        csum = csum + s*em(i)/(floti+one)
 160    s = -s

      s = one/em0
      cmn12%el(1) = one
      do 170 i = 1, cmn12%nq
 170    cmn12%el(i+1) = s*em(i)/real(i)
      xi = hsum/cmn12%h
      cmn12%tq(2) = xi*em0/csum
      cmn12%tq(5) = xi/cmn12%el(cmn12%l)
      if (cmn12%nqwait .ne. 1) go to 300

      rxi = one/xi
      do 180 iback = 1, cmn12%nq
        i = (cmn12%l + 1) - iback
 180    em(i) = em(i) + em(i-1)*rxi

      s = one
      csum = zero
      do 190 i = 1, cmn12%l
        csum = csum + s*em(i)/real(i+1)
 190    s = -s
      cmn12%tq(3) = flotl*em0/csum
      go to 300


 200  do 210 i = 3, cmn12%l
 210    cmn12%el(i) = zero
      cmn12%el(1) = one
      cmn12%el(2) = one
      alph0 = -one
      ahatn0 = -one
      hsum = cmn12%h
      rxi = one
      rxis = one
      if (cmn12%nq .eq. 1) go to 240
      do 230 j = 1, nqm2

        hsum = hsum + cmn12%tau(j)
        rxi = cmn12%h/hsum
        jp1 = j + 1
        alph0 = alph0 - one/real(jp1)
        do 220 iback = 1, jp1
          i = (j + 3) - iback
 220      cmn12%el(i) = cmn12%el(i) + cmn12%el(i-1)*rxi
 230    continue
      alph0 = alph0 - one/real(cmn12%nq)
      rxis = -cmn12%el(2) - alph0
      hsum = hsum + cmn12%tau(nqm1)
      rxi = cmn12%h/hsum
      ahatn0 = -cmn12%el(2) - rxi
      do 235 iback = 1, cmn12%nq
        i = (cmn12%nq + 2) - iback
 235    cmn12%el(i) = cmn12%el(i) + cmn12%el(i-1)*rxis
 240  t1 = one - ahatn0 + alph0
      t2 = one + real(cmn12%nq)*t1
      cmn12%tq(2) = abs(alph0*t2/t1)
      cmn12%tq(5) = abs(t2/(cmn12%el(cmn12%l)*rxi/rxis))
      if (cmn12%nqwait .ne. 1) go to 300
      cnqm1 = rxis/cmn12%el(cmn12%l)
      t3 = alph0 + one/real(cmn12%nq)
      t4 = ahatn0 + rxi
      elp = t3/(one - t4 + t3)
      cmn12%tq(1) = abs(elp/cnqm1)
      hsum = hsum + cmn12%tau(cmn12%nq)
      rxi = cmn12%h/hsum
      t5 = alph0 - one/real(cmn12%nq+1)
      t6 = ahatn0 - rxi
      elp = t2/(one - t6 + t5)
      cmn12%tq(3) = abs(elp*rxi*(flotl + one)*t5)
 300  cmn12%tq(4) = cortes*cmn12%tq(2)
      return

      end subroutine dvset

      subroutine dvjust (yh, ldyh, iord, cmn12)
      implicit none
      double precision yh
      integer ldyh, iord
      dimension yh(ldyh,*)
      type( dvode_cmn_vars ) :: cmn12

































      double precision alph0, alph1, hsum, one, prod, t1, xi,xiold, zero
      integer i, iback, j, jp1, lp1, nqm1, nqm2, nqp1




      save one, zero










      data one /1.0d0/, zero /0.0d0/

      if ((cmn12%nq .eq. 2) .and. (iord .ne. 1)) return
      nqm1 = cmn12%nq - 1
      nqm2 = cmn12%nq - 2
      go to (100, 200), cmn12%meth




 100  continue
      if (iord .eq. 1) go to 180

      do 110 j = 1, cmn12%lmax
 110    cmn12%el(j) = zero
      cmn12%el(2) = one
      hsum = zero
      do 130 j = 1, nqm2

        hsum = hsum + cmn12%tau(j)
        xi = hsum/cmn12%hscal
        jp1 = j + 1
        do 120 iback = 1, jp1
          i = (j + 3) - iback
 120      cmn12%el(i) = cmn12%el(i)*xi + cmn12%el(i-1)
 130    continue

      do 140 j = 2, nqm1
 140    cmn12%el(j+1) = real(cmn12%nq)*cmn12%el(j)/real(j)

      do 170 j = 3, cmn12%nq
        do 160 i = 1, cmn12%n
 160      yh(i,j) = yh(i,j) - yh(i,cmn12%l)*cmn12%el(j)
 170    continue
      return


 180  continue
      lp1 = cmn12%l + 1
      do 190 i = 1, cmn12%n
 190    yh(i,lp1) = zero
      return




 200  continue
      if (iord .eq. 1) go to 300

      do 210 j = 1, cmn12%lmax
 210    cmn12%el(j) = zero
      cmn12%el(3) = one
      hsum = zero
      do 230 j = 1,nqm2

        hsum = hsum + cmn12%tau(j)
        xi = hsum/cmn12%hscal
        jp1 = j + 1
        do 220 iback = 1, jp1
          i = (j + 4) - iback
 220      cmn12%el(i) = cmn12%el(i)*xi + cmn12%el(i-1)
 230    continue

      do 250 j = 3,cmn12%nq
        do 240 i = 1, cmn12%n
 240      yh(i,j) = yh(i,j) - yh(i,cmn12%l)*cmn12%el(j)
 250    continue
      return

 300  do 310 j = 1, cmn12%lmax
 310    cmn12%el(j) = zero
      cmn12%el(3) = one
      alph0 = -one
      alph1 = one
      prod = one
      xiold = one
      hsum = cmn12%hscal
      if (cmn12%nq .eq. 1) go to 340
      do 330 j = 1, nqm1

        jp1 = j + 1
        hsum = hsum + cmn12%tau(jp1)
        xi = hsum/cmn12%hscal
        prod = prod*xi
        alph0 = alph0 - one/real(jp1)
        alph1 = alph1 + one/xi
        do 320 iback = 1, jp1
          i = (j + 4) - iback
 320      cmn12%el(i) = cmn12%el(i)*xiold + cmn12%el(i-1)
        xiold = xi
 330    continue
 340  continue
      t1 = (-alph0 - alph1)/prod

      lp1 = cmn12%l + 1
      do 350 i = 1, cmn12%n
 350    yh(i,lp1) = t1*yh(i,cmn12%lmax)

      nqp1 = cmn12%nq + 1
      do 370 j = 3, nqp1
        call daxpy (cmn12%n, cmn12%el(j), yh(1,lp1), 1, yh(1,j), 1 )
 370  continue
      return

      end subroutine dvjust

      subroutine dvnlsd (y, yh, ldyh, vsav, savf, ewt, acor, iwm, wm,   &
                       f, jac, pdum, nflag, rpar, ipar, cmn12)
      implicit none
      external f, jac, pdum
      double precision y, yh, vsav, savf, ewt, acor, wm, rpar
      integer ldyh, iwm, nflag, ipar
      dimension y(*), yh(ldyh,*), vsav(*), savf(*), ewt(*), acor(*),   &
                iwm(*), wm(*), rpar(*), ipar(*)
      type( dvode_cmn_vars ) :: cmn12












































































      double precision ccmax, crdown, cscale, dcon, del, delp, one,   &
           rdiv, two, zero
      integer i, ierpj, iersl, m, maxcor, msbp









      save ccmax, crdown, maxcor, msbp, rdiv, one, two, zero











      data ccmax /0.3d0/, crdown /0.3d0/, maxcor /3/, msbp /20/,   &
           rdiv  /2.0d0/
      data one /1.0d0/, two /2.0d0/, zero /0.0d0/





      if (cmn12%jstart .eq. 0) cmn12%nslp = 0
      if (nflag .eq. 0) cmn12%icf = 0
      if (nflag .eq. -2) cmn12%ipup = cmn12%miter
      if ( (cmn12%jstart .eq. 0) .or. (cmn12%jstart .eq. -1) ) cmn12%ipup = cmn12%miter

      if (cmn12%miter .eq. 0) then
        cmn12%crate = one
        go to 220
      endif






      cmn12%drc = abs(cmn12%rc-one)
      if (cmn12%drc .gt. ccmax .or. cmn12%nst .ge. cmn12%nslp+msbp) cmn12%ipup = cmn12%miter






 220  m = 0
      delp = zero
      call dcopy (cmn12%n, yh(1,1), 1, y, 1 )
      call f (cmn12%n, cmn12%tn, y, savf, rpar, ipar)
      cmn12%nfe = cmn12%nfe + 1
      if (cmn12%ipup .le. 0) go to 250





      call dvjac (y, yh, ldyh, ewt, acor, savf, wm, iwm, f, jac, ierpj,   &
                 rpar, ipar, cmn12)
      cmn12%ipup = 0
      cmn12%rc = one
      cmn12%drc = zero
      cmn12%crate = one
      cmn12%nslp = cmn12%nst

      if (ierpj .ne. 0) go to 430
 250  do 260 i = 1,cmn12%n
 260    acor(i) = zero

 270  if (cmn12%miter .ne. 0) go to 350




      do 280 i = 1,cmn12%n
 280    savf(i) = cmn12%rl1*(cmn12%h*savf(i) - yh(i,2))
      do 290 i = 1,cmn12%n
 290    y(i) = savf(i) - acor(i)
      del = dvnorm (cmn12%n, y, ewt)
      do 300 i = 1,cmn12%n
 300    y(i) = yh(i,1) + savf(i)
      call dcopy (cmn12%n, savf, 1, acor, 1)
      go to 400






 350  do 360 i = 1,cmn12%n
 360    y(i) = (cmn12%rl1*cmn12%h)*savf(i) - (cmn12%rl1*yh(i,2) + acor(i))
      call dvsol (wm, iwm, y, iersl, cmn12)
      cmn12%nni = cmn12%nni + 1
      if (iersl .gt. 0) go to 410
      if (cmn12%meth .eq. 2 .and. cmn12%rc .ne. one) then
        cscale = two/(one + cmn12%rc)
        call dscal (cmn12%n, cscale, y, 1)
      endif
      del = dvnorm (cmn12%n, y, ewt)
      call daxpy (cmn12%n, one, y, 1, acor, 1)
      do 380 i = 1,cmn12%n
 380    y(i) = yh(i,1) + acor(i)




 400  if (m .ne. 0) cmn12%crate = max(crdown*cmn12%crate,del/delp)
      dcon = del*min(one,cmn12%crate)/cmn12%tq(4)
      if (dcon .le. one) go to 450
      m = m + 1
      if (m .eq. maxcor) go to 410
      if (m .ge. 2 .and. del .gt. rdiv*delp) go to 410
      delp = del
      call f (cmn12%n, cmn12%tn, y, savf, rpar, ipar)
      cmn12%nfe = cmn12%nfe + 1
      go to 270

 410  if (cmn12%miter .eq. 0 .or. cmn12%jcur .eq. 1) go to 430
      cmn12%icf = 1
      cmn12%ipup = cmn12%miter
      go to 220

 430  continue
      nflag = -1
      cmn12%icf = 2
      cmn12%ipup = cmn12%miter
      return


 450  nflag = 0
      cmn12%jcur = 0
      cmn12%icf = 0
      if (m .eq. 0) cmn12%acnrm = del
      if (m .gt. 0) cmn12%acnrm = dvnorm (cmn12%n, acor, ewt)
      return

      end subroutine dvnlsd

      subroutine dvjac (y, yh, ldyh, ewt, ftem, savf, wm, iwm, f, jac,   &
                       ierpj, rpar, ipar, cmn12)
      implicit none
      external f, jac
      double precision y, yh, ewt, ftem, savf, wm, rpar
      integer ldyh, iwm, ierpj, ipar
      dimension y(*), yh(ldyh,*), ewt(*), ftem(*), savf(*),   &
         wm(*), iwm(*), rpar(*), ipar(*)
      type( dvode_cmn_vars ) :: cmn12











































































      double precision con, di, fac, hrl1, one, pt1, r, r0, srur, thou,   &
           yi, yj, yjj, zero
      integer i, i1, i2, ier, ii, j, j1, jj, jok, lenp, mba, mband,   &
              meb1, meband, ml, ml3, mu, np1









      save one, pt1, thou, zero











      data one /1.0d0/, thou /1000.0d0/, zero /0.0d0/, pt1 /0.1d0/

      ierpj = 0
      hrl1 = cmn12%h*cmn12%rl1

      jok = cmn12%jsv
      if (cmn12%jsv .eq. 1) then
        if (cmn12%nst .eq. 0 .or. cmn12%nst .gt. cmn12%nslj+cmn12%msbj) jok = -1
        if (cmn12%icf .eq. 1 .and. cmn12%drc .lt. cmn12%ccmxj) jok = -1
        if (cmn12%icf .eq. 2) jok = -1
      endif


      if (jok .eq. -1 .and. cmn12%miter .eq. 1) then

      cmn12%nje = cmn12%nje + 1
      cmn12%nslj = cmn12%nst
      cmn12%jcur = 1
      lenp = cmn12%n*cmn12%n
      do 110 i = 1,lenp
 110    wm(i+2) = zero
      call jac (cmn12%n, cmn12%tn, y, 0, 0, wm(3), cmn12%n, rpar, ipar)
      if (cmn12%jsv .eq. 1) call dcopy (lenp, wm(3), 1, wm(cmn12%locjs), 1)
      endif

      if (jok .eq. -1 .and. cmn12%miter .eq. 2) then

      cmn12%nje = cmn12%nje + 1
      cmn12%nslj = cmn12%nst
      cmn12%jcur = 1
      fac = dvnorm (cmn12%n, savf, ewt)
      r0 = thou*abs(cmn12%h)*cmn12%uround*real(cmn12%n)*fac
      if (r0 .eq. zero) r0 = one
      srur = wm(1)
      j1 = 2
      do 230 j = 1,cmn12%n
        yj = y(j)
        r = max(srur*abs(yj),r0/ewt(j))
        y(j) = y(j) + r
        fac = one/r
        call f (cmn12%n, cmn12%tn, y, ftem, rpar, ipar)
        do 220 i = 1,cmn12%n
 220      wm(i+j1) = (ftem(i) - savf(i))*fac
        y(j) = yj
        j1 = j1 + cmn12%n
 230    continue
      cmn12%nfe = cmn12%nfe + cmn12%n
      lenp = cmn12%n*cmn12%n
      if (cmn12%jsv .eq. 1) call dcopy (lenp, wm(3), 1, wm(cmn12%locjs), 1)
      endif

      if (jok .eq. 1 .and. (cmn12%miter .eq. 1 .or. cmn12%miter .eq. 2)) then
      cmn12%jcur = 0
      lenp = cmn12%n*cmn12%n
      call dcopy (lenp, wm(cmn12%locjs), 1, wm(3), 1)
      endif

      if (cmn12%miter .eq. 1 .or. cmn12%miter .eq. 2) then

      con = -hrl1
      call dscal (lenp, con, wm(3), 1)
      j = 3
      np1 = cmn12%n + 1
      do 250 i = 1,cmn12%n
        wm(j) = wm(j) + one
 250    j = j + np1
      cmn12%nlu = cmn12%nlu + 1
      call dgefa (wm(3), cmn12%n, cmn12%n, iwm(31), ier)
      if (ier .ne. 0) ierpj = 1
      return
      endif


      if (cmn12%miter .eq. 3) then

      cmn12%nje = cmn12%nje + 1
      cmn12%jcur = 1
      wm(2) = hrl1
      r = cmn12%rl1*pt1
      do 310 i = 1,cmn12%n
 310    y(i) = y(i) + r*(cmn12%h*savf(i) - yh(i,2))
      call f (cmn12%n, cmn12%tn, y, wm(3), rpar, ipar)
      cmn12%nfe = cmn12%nfe + 1
      do 320 i = 1,cmn12%n
        r0 = cmn12%h*savf(i) - yh(i,2)
        di = pt1*r0 - cmn12%h*(wm(i+2) - savf(i))
        wm(i+2) = one
        if (abs(r0) .lt. cmn12%uround/ewt(i)) go to 320
        if (abs(di) .eq. zero) go to 330
        wm(i+2) = pt1*r0/di
 320    continue
      return
 330  ierpj = 1
      return
      endif



      ml = iwm(1)
      mu = iwm(2)
      ml3 = ml + 3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*cmn12%n

      if (jok .eq. -1 .and. cmn12%miter .eq. 4) then

      cmn12%nje = cmn12%nje + 1
      cmn12%nslj = cmn12%nst
      cmn12%jcur = 1
      do 410 i = 1,lenp
 410    wm(i+2) = zero
      call jac (cmn12%n, cmn12%tn, y, ml, mu, wm(ml3), meband, rpar, ipar)
      if (cmn12%jsv .eq. 1)   &
         call dacopy (mband, cmn12%n, wm(ml3), meband, wm(cmn12%locjs), mband)
      endif

      if (jok .eq. -1 .and. cmn12%miter .eq. 5) then

      cmn12%nje = cmn12%nje + 1
      cmn12%nslj = cmn12%nst
      cmn12%jcur = 1
      mba = min(mband,cmn12%n)
      meb1 = meband - 1
      srur = wm(1)
      fac = dvnorm (cmn12%n, savf, ewt)
      r0 = thou*abs(cmn12%h)*cmn12%uround*real(cmn12%n)*fac
      if (r0 .eq. zero) r0 = one
      do 560 j = 1,mba
        do 530 i = j,cmn12%n,mband
          yi = y(i)
          r = max(srur*abs(yi),r0/ewt(i))
 530      y(i) = y(i) + r
        call f (cmn12%n, cmn12%tn, y, ftem, rpar, ipar)
        do 550 jj = j,cmn12%n,mband
          y(jj) = yh(jj,1)
          yjj = y(jj)
          r = max(srur*abs(yjj),r0/ewt(jj))
          fac = one/r
          i1 = max(jj-mu,1)
          i2 = min(jj+ml,cmn12%n)
          ii = jj*meb1 - ml + 2
          do 540 i = i1,i2
 540        wm(ii+i) = (ftem(i) - savf(i))*fac
 550      continue
 560    continue
      cmn12%nfe = cmn12%nfe + mba
      if (cmn12%jsv .eq. 1)   &
         call dacopy (mband, cmn12%n, wm(ml3), meband, wm(cmn12%locjs), mband)
      endif

      if (jok .eq. 1) then
      cmn12%jcur = 0
      call dacopy (mband, cmn12%n, wm(cmn12%locjs), mband, wm(ml3), meband)
      endif


      con = -hrl1
      call dscal (lenp, con, wm(3), 1 )
      ii = mband + 2
      do 580 i = 1,cmn12%n
        wm(ii) = wm(ii) + one
 580    ii = ii + meband
      cmn12%nlu = cmn12%nlu + 1
      call dgbfa (wm(3), meband, cmn12%n, ml, mu, iwm(31), ier)
      if (ier .ne. 0) ierpj = 1
      return



      end subroutine dvjac

      subroutine dacopy (nrow, ncol, a, nrowa, b, nrowb)
      double precision a, b
      integer nrow, ncol, nrowa, nrowb
      dimension a(nrowa,ncol), b(nrowb,ncol)












      integer ic

      do 20 ic = 1,ncol
        call dcopy (nrow, a(1,ic), 1, b(1,ic), 1)
 20     continue

      return

      end subroutine dacopy

      subroutine dvsol (wm, iwm, x, iersl, cmn12)
      implicit none
      double precision wm, x
      integer iwm, iersl
      dimension wm(*), iwm(*), x(*)
      type( dvode_cmn_vars ) :: cmn12












































      integer i, meband, ml, mu
      double precision di, hrl1, one, phrl1, r, zero




      save one, zero










      data one /1.0d0/, zero /0.0d0/

      iersl = 0
      go to (100, 100, 300, 400, 400), cmn12%miter
 100  call dgesl (wm(3), cmn12%n, cmn12%n, iwm(31), x, 0)
      return

 300  phrl1 = wm(2)
      hrl1 = cmn12%h*cmn12%rl1
      wm(2) = hrl1
      if (hrl1 .eq. phrl1) go to 330
      r = hrl1/phrl1
      do 320 i = 1,cmn12%n
        di = one - r*(one - one/wm(i+2))
        if (abs(di) .eq. zero) go to 390
 320    wm(i+2) = one/di

 330  do 340 i = 1,cmn12%n
 340    x(i) = wm(i+2)*x(i)
      return
 390  iersl = 1
      return

 400  ml = iwm(1)
      mu = iwm(2)
      meband = 2*ml + mu + 1
      call dgbsl (wm(3), meband, cmn12%n, ml, mu, iwm(31), x, 0)
      return

      end subroutine dvsol































































      subroutine dewset (n, itol, rtol, atol, ycur, ewt)





















      integer n, itol
      integer i
      double precision rtol, atol, ycur, ewt
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

      end subroutine dewset

      double precision function dvnorm (n, v, w)





















      integer n,   i
      double precision v, w,   sum
      dimension v(n), w(n)


      sum = 0.0d0
      do 10 i = 1,n
 10     sum = sum + (v(i)*w(i))**2
      dvnorm = sqrt(sum/n)
      return

      end function dvnorm

      subroutine xerrwd (msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)



























































      double precision r1, r2
      integer nmes, nerr, level, ni, i1, i2, nr
      character*(*) msg





      integer lunit,        mesflg




      lunit = ixsav (1, 0, .false.)
      mesflg = ixsav (2, 0, .false.)
      if (mesflg .eq. 0) go to 100



      write (lunit,10)  msg
 10   format(1x,a)
      if (ni .eq. 1) write (lunit, 20) i1
 20   format(6x,'in above message,  i1 =',i10)
      if (ni .eq. 2) write (lunit, 30) i1,i2
 30   format(6x,'in above message,  i1 =',i10,3x,'i2 =',i10)
      if (nr .eq. 1) write (lunit, 40) r1
 40   format(6x,'in above message,  r1 =',d21.13)
      if (nr .eq. 2) write (lunit, 50) r1,r2
 50   format(6x,'in above,  r1 =',d21.13,3x,'r2 =',d21.13)



 100  if (level .ne. 2) return
         call wrf_error_fatal3("<stdin>",3723,&
" ")


      end subroutine xerrwd

      subroutine xsetf (mflag)




























      integer mflag, junk


      if (mflag .eq. 0 .or. mflag .eq. 1) junk = ixsav (2,mflag,.true.)
      return

      end subroutine xsetf

      subroutine xsetun (lun)


























      integer lun, junk


      if (lun .gt. 0) junk = ixsav (1,lun,.true.)
      return

      end subroutine xsetun

      integer function ixsav (ipar, ivalue, iset)













































      logical iset
      integer ipar, ivalue



      integer         lunit, mesflg




      save lunit, mesflg
      data lunit/-1/, mesflg/1/


      if (ipar .eq. 1) then
        if (lunit .eq. -1) lunit = iumach()
        ixsav = lunit
        if (iset) lunit = ivalue
        endif

      if (ipar .eq. 2) then
        ixsav = mesflg
        if (iset) mesflg = ivalue
        endif

      return

      end function ixsav

      integer function iumach()


























      iumach = 6

      return

      end function iumach

      double precision function dumach ()



























      double precision u, comp

      u = 1.0d0
 10   u = u*0.5d0
      call dumsum(1.0d0, u, comp)
      if (comp .ne. 1.0d0) go to 10
      dumach = u*2.0d0
      return

      end function dumach
      subroutine dumsum(a,b,c)

      double precision a, b, c
      c = a + b
      return
      end subroutine dumsum












      subroutine daxpy(n,da,dx,incx,dy,incy)






      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n

      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20




      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return






   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end subroutine daxpy



      subroutine  dcopy(n,dx,incx,dy,incy)






      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20




      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return






   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end subroutine  dcopy


      double precision function ddot(n,dx,incx,dy,incy)






      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n

      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20




      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return






   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +   &
         dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end function ddot



      double precision function dnrm2 ( n, x, incx )

      integer                           incx, n

      double precision                  x( * )















      double precision      one         , zero
      parameter           ( one = 1.0d+0, zero = 0.0d+0 )

      integer               ix
      double precision      absxi, norm, scale, ssq

      intrinsic             abs, sqrt


      if( n.lt.1 .or. incx.lt.1 )then
         norm  = zero
      else if( n.eq.1 )then
         norm  = abs( x( 1 ) )
      else
         scale = zero
         ssq   = one




         do 10, ix = 1, 1 + ( n - 1 )*incx, incx
            if( x( ix ).ne.zero )then
               absxi = abs( x( ix ) )
               if( scale.lt.absxi )then
                  ssq   = one   + ssq*( scale/absxi )**2
                  scale = absxi
               else
                  ssq   = ssq   +     ( absxi/scale )**2
               end if
            end if
   10    continue
         norm  = scale * sqrt( ssq )
      end if

      dnrm2 = norm
      return



      end function dnrm2



      subroutine  dscal(n,da,dx,incx)







      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20



      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return






   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end subroutine  dscal



      integer function idamax(n,dx,incx)






      double precision dx(*),dmax
      integer i,incx,ix,n

      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20



      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return



   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end function idamax



      subroutine dgbfa(abd,lda,n,ml,mu,ipvt,info)
      integer lda,n,ml,mu,ipvt(1),info
      double precision abd(lda,1)



















































































      double precision t


      integer i,       i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1


      m = ml + mu + 1
      info = 0



      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0d0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0



      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1



         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0d0
   40       continue
   50    continue



         lm = min0(ml,n-k)
         l = idamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m



         if (abd(l,k) .eq. 0.0d0) go to 100



            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue



            t = -1.0d0/abd(m,k)
            call dscal(lm,t,abd(m+1,k),1)



            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0d0) info = n
      return
      end subroutine dgbfa



      subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job)
      integer lda,n,ml,mu,ipvt(1),job
      double precision abd(lda,1),b(1)


































































      double precision      t
      integer k,kb,l,la,lb,lm,m,nm1

      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50




         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue



         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call daxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue




         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = ddot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   60    continue



         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end subroutine dgbsl



      subroutine dgefa(a,lda,n,ipvt,info)


      integer lda,n,ipvt(n),info

      double precision a(lda,n)














































      double precision t


      integer        j,k,kp1,l,nm1




      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1



         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l



         if (a(l,k) .eq. 0.0d0) go to 40



            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue



            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)



            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end subroutine dgefa



      subroutine dgesl(a,lda,n,ipvt,b,job)


      integer lda,n,ipvt(n),job

      double precision a(lda,n),b(n)



























































      double precision      t
      integer k,kb,l,nm1

      nm1 = n - 1
      if (job .ne. 0) go to 50




         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue



         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue




         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue



         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end subroutine dgesl



      end module module_cmu_dvode_solver

