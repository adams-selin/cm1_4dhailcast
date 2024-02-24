  MODULE hailcast_module

!-----------------------------------------------------------------------------
!
!  HAILCAST hail trajectory model
!  https://github.com/adams-selin/cm1_4dhailcast
!
!  (c)2019 - Atmospheric and Environmental Research
!
!-----------------------------------------------------------------------------
!
!   Created by Becky Adams-Selin as an embedding of the Adams-Selin and Ziegler
!       (2016, MWR) HAILCAST hail trajectory model within CM1. This code
!       reuses the CM1 parcel code to instead move around hailstones.
!
!   A few differences with the ASZ16 hail trajectory model: this subroutine
!    allows for melting at every timestep.  If inside cloud, the melting is
!    determined via the RH87 heat balance equation. If outside or below cloud,
!    melting is determined using Eq. 3 of Goyer et al. (1969), assuming a 
!    spherical hailstone melting in dry air. This replaces the mean below-
!    cloud melting calculations performed in ASZ16.
!
!   The CM1 parcel code was repurposed here for advecting the hailstone.
!    Hailstones are first advected using a two-step Runge-Kutta process.
!    Hailstone terminal velocity is calculated and incorporated at both steps. 
!    Following the advection, physical variables are interpolated at the new
!    hailstone location and hailstone growth/melting is determined via the
!    hailstone_trajectory subroutine.  Once the hailstone reaches one of 
!    the edges of the domain, hailstone_trajectory is no longer called and
!    the hailstone remains there for the rest of the simulation.
!
!   Due to the nature of hailstones, some of the largest trajectory speeds
!    occur near the surface.  I had to add a check for below ground heights
!    within the 1st rk loop, not just after the 2nd as in the parcel code.
!
!   One possibility I haven't accounted for: I'm still allowing the hailstone
!    to fall and be advected even if it has reached the edge of the domain.
!    That could have an impact if it shot off the top, as it might eventually
!    fall back down into the simulation and start growing again.  I'll see
!    how much of a problem stratospheric hailstones are.
!
!   Finally, the implementation of the code for sigma vertical coordinates
!    is not complete and should not be used.
!
!-----------------------------------------------------------------------------

  implicit none

  private
  public :: hailcast_driver,setup_hailcast_vars,hailcast_write,gethailtrajzs

  CONTAINS
  
      subroutine hailcast_driver(dt,xh,uh,ruh,xf,yh,vh,rvh,yf,zh,mh,rmh,zf,mf,zs,    &
                               sigma,sigmaf,znt,                                     &
                               temp,rho,prs,qv,qia,qwa,ua,va,wa,haildata,            &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,        &
                               nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,reqs_p, &
                               tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2)
                               
      use input
      use constants
      use bc_module
      use comm_module
#ifdef MPI
      use mpi
#endif
      implicit none

!-----------------------------------------------------------------------
!  This subroutine updates the hailstone locations
!-----------------------------------------------------------------------

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(ib:ie,jb:je) :: znt
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: temp,rho,prs,qv,qia,qwa
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(nhailtrajs,nhvals) :: haildata
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      integer, intent(inout), dimension(rmp) :: reqs_p
      real, intent(inout), dimension(cmp,jmp,kmt) :: tkw1,tkw2,tke1,tke2
      real, intent(inout), dimension(imp,cmp,kmt) :: tks1,tks2,tkn1,tkn2

      integer :: n,np,i,j,k,iflag,jflag,kflag
      real :: uval,vval,wval,rx,ry,rz,w1,w2,w3,w4,w5,w6,w7,w8,wsum
      real :: tempval,rhoval,prsval,qvval,qiaval,qwaval
      real :: rxu,ryv,rzw,rxs,rys,rzs
      real :: x3d,y3d,z3d
      real :: hdense3d !hailstone density
      real*8 :: hd3d !hailstone diameter
      integer :: nrkp
      real :: dt2,uu1,vv1,ww1
      real :: z0,rznt,var
      real :: sigdot,sig1,zsp,sig3d

      logical, parameter :: debug = .true.

      real, parameter :: a_shedd = 1.3978
      !Shedd et al. 2021 constant (ignoring the y intercept in Fig. 9)
      !For an oblate spheroid, Dmax = a_shedd * Deq, where Deq
      ! is the equivalent spherical radius.
      !Used in TERMINL, HEATBUD, and MASSAGR depending on namelist
      ! options.


!----------------------------------------------------------------------
!  apply bottom/top boundary conditions:
!  [Note:  for u,v the array index (i,j,0) means the surface, ie z=0]
!     (for the parcel subroutines only!)


!$omp parallel do default(shared)  &
!$omp private(i,j)
  DO j=jb,je+1

    IF(bbc.eq.1)THEN
      ! free slip ... extrapolate:
      IF(j.le.je)THEN
          do i=ib,ie+1
            ua(i,j,0) = cgs1*ua(i,j,1)+cgs2*ua(i,j,2)+cgs3*ua(i,j,3)
          enddo
      ENDIF
      do i=ib,ie
        va(i,j,0) = cgs1*va(i,j,1)+cgs2*va(i,j,2)+cgs3*va(i,j,3)
      enddo
    ELSEIF(bbc.eq.2)THEN
      ! no slip:
      if( imove.eq.1 )then
        IF(j.le.je)THEN
            do i=ib,ie+1
              ua(i,j,0) = 0.0 - umove
            enddo
        ENDIF
        do i=ib,ie
          va(i,j,0) = 0.0 - vmove
        enddo
      else
        IF(j.le.je)THEN
            do i=ib,ie+1
              ua(i,j,0) = 0.0
            enddo
        ENDIF
        do i=ib,ie
          va(i,j,0) = 0.0
        enddo
      endif
    ELSEIF(bbc.eq.3)THEN
      ! u,v near sfc are determined below using log-layer equations
    ENDIF

!----------

    IF(tbc.eq.1)THEN
      ! free slip ... extrapolate:
      IF(j.le.je)THEN
          do i=ib,ie+1
            ua(i,j,nk+1) = cgt1*ua(i,j,nk)+cgt2*ua(i,j,nk-1)+cgt3*ua(i,j,nk-2)
          enddo
      ENDIF
      do i=ib,ie
        va(i,j,nk+1) = cgt1*va(i,j,nk)+cgt2*va(i,j,nk-1)+cgt3*va(i,j,nk-2)
      enddo
    ELSEIF(tbc.eq.2)THEN
      ! no slip:
      IF(j.le.je)THEN
          do i=ib,ie+1
            ua(i,j,nk+1) = 0.0
          enddo
      ENDIF
      do i=ib,ie
        va(i,j,nk+1) = 0.0
      enddo
    ENDIF

!----------

      IF(j.le.je)THEN
      do i=ib,ie
        wa(i,j,nk+1) = 0.0
      enddo
      ENDIF

  ENDDO

!----------------------------------------------------------------------
!  Determine state variables at the top/bottom boundaries


  call prepcorners(temp,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1), &
                             se1(1),se2(1),pw1,pw2,pe1,pe2,ps1,ps2,&
                             pn1,pn2,reqs_p,1)
  call prepcorners(rho,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
  call prepcorners(prs,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
  call prepcorners(qv,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
  call prepcorners(qia,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)
  call prepcorners(qwa,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1), &
                             pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,reqs_p,1)

!!$omp parallel do default(shared)  &
!!$omp private(i,j)
!  DO j=jb,je
!      DO i=ib,ie
!        temp(i,j,0) = cgs1*temp(i,j,1)+cgs2*temp(i,j,2)+cgs3*temp(i,j,3)
!        rho(i,j,0) = cgs1*rho(i,j,1)+cgs2*rho(i,j,2)+cgs3*rho(i,j,3)
!        prs(i,j,0) = cgs1*prs(i,j,1)+cgs2*prs(i,j,2)+cgs3*prs(i,j,3)
!        qv(i,j,0)  = cgs1*qv(i,j,1)+ cgs2*qv(i,j,2)+ cgs3*qv(i,j,3)
!        qia(i,j,0) = cgs1*qia(i,j,1)+cgs2*qia(i,j,2)+cgs3*qia(i,j,3)
!        qwa(i,j,0) = cgs1*qwa(i,j,1)+cgs2*qwa(i,j,2)+cgs3*qwa(i,j,3)
!      ENDDO
!  ENDDO

!  DO j=jb,je
!      DO i=ib,ie
!        temp(i,j,nk+1) = cgt1*temp(i,j,nk)+cgt2*temp(i,j,nk-1)+cgt3*temp(i,j,nk-2)
!        rho(i,j,nk+1) = cgt1*rho(i,j,nk)+cgt2*rho(i,j,nk-1)+cgt3*rho(i,j,nk-2)
!        prs(i,j,nk+1) = cgt1*prs(i,j,nk)+cgt2*prs(i,j,nk-1)+cgt3*prs(i,j,nk-2)
!        qv(i,j,nk+1) = cgt1*qv(i,j,nk)+cgt2*qv(i,j,nk-1)+cgt3*qv(i,j,nk-2)
!        qia(i,j,nk+1) = cgt1*qia(i,j,nk)+cgt2*qia(i,j,nk-1)+cgt3*qia(i,j,nk-2)
!        qwa(i,j,nk+1) = cgt1*qwa(i,j,nk)+cgt2*qwa(i,j,nk-1)+cgt3*qwa(i,j,nk-2)
!      ENDDO
!  ENDDO


!----------------------------------------------------------------------
!  Loop through all hailstones:  if you have it, update its location:

    dt2 = dt/2.0

    nploop:  &
    DO np=1,nhailtrajs

      x3d = haildata(np,hailx)
      y3d = haildata(np,haily)
      if( .not. terrain_flag )then
        z3d = haildata(np,hailz)
      else
        sig3d = haildata(np,hailsig)
      endif
      hd3d = haildata(np,haild)
      hdense3d = haildata(np,haildense)

      iflag = -100
      jflag = -100
      kflag = 0



      ! cm1r19:  skip if we already know this processor doesnt have this parcel
      haveit1:  &
      IF( x3d.ge.xf(1) .and. x3d.le.xf(ni+1) .and.  &
          y3d.ge.yf(1) .and. y3d.le.yf(nj+1) )THEN

      !if (np.eq.1282828) print *, 'hail.F: x3d,xf(1),xf(ni+1): ', x3d,xf(1),xf(ni+1)

        IF(nx.eq.1)THEN
          iflag = 1
        ELSE
          ! cm1r19:
          i = ni+1
          do while( iflag.lt.0 .and. i.gt.1 )
            i = i-1
            if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
              iflag = i
            endif
          enddo
        ENDIF

        IF(axisymm.eq.1.or.ny.eq.1)THEN
          jflag = 1
        ELSE
          ! cm1r19:
          j = nj+1
          do while( jflag.lt.0 .and. j.gt.1 )
            j = j-1
            if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
              jflag = j
            endif
          enddo
        ENDIF

      ENDIF  haveit1

#ifdef MPI
      ! check for conflict:
      IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
          (jflag.ge.1.and.jflag.le.nj) )THEN
          
!            if (np.eq.3089694) print *, 'conflict check: i,iflag,haildatax,xf: ', i,iflag,&
!             haildata(np,hailx),xf(iflag+1)
 
          IF( iflag.eq.ni .and. haildata(np,hailx).eq.xf(iflag+1) .and. nodex.gt.1 .and.  myi.ne.nodex ) iflag = -1
          IF( jflag.eq.nj .and. haildata(np,haily).eq.yf(jflag+1) .and. nodey.gt.1 .and.  myj.ne.nodey ) jflag = -1
      ENDIF
#endif


      !myparcel:  IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
      !               (jflag.ge.1.and.jflag.le.nj) )THEN
      
      !Added an extra check here for hail trajectories requiring
      !z3d to be above 0.  If not, we can skip the rest of the 
      !code.
      myparcel:  IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
                     (jflag.ge.1.and.jflag.le.nj) .and.   &
                     (z3d.gt.0) )THEN

      !Added an check here for mtime greater than hail_start, but 
      ! then decided against it
      !myparcel:  IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
                     !(jflag.ge.1.and.jflag.le.nj) .and.   &
                     !(z3d.gt.0) )THEN

          rkloop:  DO nrkp = 1,2

          !if (np.eq.313793) write(outfile,*) ' '
          !if (np.eq.313793) write(outfile,'(a28,i7,3f9.1)'), 'start rkl nrkp,x3d,y3d,z3d: ', nrkp, x3d, y3d, z3d

            IF( nrkp.eq.1 )THEN
                i=iflag
                j=jflag
            ELSE
                iflag = -100
                jflag = -100
                IF(nx.eq.1)THEN
                  iflag = 1
                ELSE
                  ! cm1r19:
                  i = ni+2
                  do while( iflag.lt.0 .and. i.gt.0 )
                    i = i-1
                    if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
                      iflag = i
                    endif
                  enddo
                ENDIF
                IF(axisymm.eq.1.or.ny.eq.1)THEN
                  jflag = 1
                ELSE
                  do j=0,nj+1
                    if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) ) jflag=j
                  enddo
                  ! cm1r19:
                  j = nj+2
                  do while( jflag.lt.0 .and. j.gt.0 )
                    j = j-1
                    if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
                      jflag = j
                    endif
                  enddo
                ENDIF
                i=iflag
                j=jflag
            ENDIF


            IF(debug)THEN
                if( i.lt.0 .or. i.gt.(ni+1) .or. j.lt.0 .or. j.gt.(nj+1) )then
                  print *,'  myid,i,j = ',myid,i,j
                  print *,'  x,x1     = ',x3d,haildata(np,hailx)
                  print *,'  y,y1     = ',y3d,haildata(np,haily)
                  do i=0,ni+1
                    print *,i,abs(xh(i)-x3d),0.5*dx*ruh(i)
                  enddo
                  do j=0,nj+1
                    print *,j,abs(yh(j)-y3d),0.5*dy*rvh(j)
                  enddo
                  call stopcm1
                endif
            ENDIF

            kflag = 1
            if( .not. terrain_flag )then
              do while( z3d.ge.zf(iflag,jflag,kflag+1) )
                kflag = kflag+1
              enddo
            else
              do while( sig3d.ge.sigmaf(kflag+1) )
                kflag = kflag+1
              enddo
            endif

            IF(debug)THEN
            if( kflag.le.0 .or. kflag.ge.(nk+1) )then
              print *,myid,nrkp,np
              print *,iflag,jflag,kflag
              print *,haildata(np,hailx),haildata(np,haily),haildata(np,hailz)
              print *,x3d,y3d,z3d
              print *,uval,vval,wval
              print *,zf(iflag,jflag,kflag),z3d,zf(iflag,jflag,kflag+1)
              print *,'  16667 '
              call stopcm1
            endif
            ENDIF
            
            
            !if (np.eq.313793) WRITE(outfile,'(a35,i6,3i4)') ' after flagset: iflag,jflag,kflag: ',np,iflag,jflag,kflag

            !----------------------------------------------------------------------
            !  Data on u points

            i=iflag
            j=jflag
            k=kflag

            if( y3d.lt.yh(j) )then
              j=j-1
            endif
            if( .not. terrain_flag )then
              if( z3d.lt.zh(iflag,jflag,k) )then
                k=k-1
              endif
              rz = ( z3d-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )
            else
              if( sig3d.lt.sigma(k) )then
                k=k-1
              endif
              rz = ( sig3d-sigma(k) )/( sigma(k+1)-sigma(k) )
            endif

            rx = ( x3d-xf(i) )/( xf(i+1)-xf(i) )
            ry = ( y3d-yh(j) )/( yh(j+1)-yh(j) )

            ! saveit:
            rxu = rx
            rys = ry
            rzs = rz

            w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
            w2=rx*(1.0-ry)*(1.0-rz)
            w3=(1.0-rx)*ry*(1.0-rz)
            w4=(1.0-rx)*(1.0-ry)*rz
            w5=rx*(1.0-ry)*rz
            w6=(1.0-rx)*ry*rz
            w7=rx*ry*(1.0-rz)
            w8=rx*ry*rz

            IF(debug)THEN
                wsum = w1+w2+w3+w4+w5+w6+w7+w8
                if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
                    ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
                    rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
                    wsum.le.0.99999 .or.                  &
                    wsum.ge.1.00001 .or.                  &
                    i.lt.0 .or. i.gt.(ni+1)   .or.        &
                    j.lt.-1 .or. j.gt.(nj+1)   .or.       &
                    k.lt.0 .or. k.gt.nk                   )then
                  print *
                  print *,'  13333a: '
                  print *,'  np          = ',np
                  print *,'  myid,i,j,k  = ',myid,i,j,k
                  print *,'  rx,ry,rz    = ',rx,ry,rz
                  print *,'  wsum        = ',wsum
                  print *,'  xf1,x3d,xf2 = ',xf(i),x3d,xf(i+1)
                  print *,'  yh1,y3d,yh2 = ',yh(j),y3d,yh(j+1)
                  print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),z3d,zh(iflag,jflag,k+1)
                  print *
                  call stopcm1
                endif
            ENDIF

            !if (np.eq.3089694) print *, ' hail.F, i, j, k: ', i, j, k
            !print *, ' hail.F, ni, nj, nk: ', ni, nj, nk
            !print *, ' hail.F, ua(i,j,k): ', ua(i,j,k), ua(i+i,j,k), ua(i,j+1,k)
            call tri_interp(ni+1,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,ua,uval)
             !if (np.eq.1282828) print *, ' hail.F, tri_interp uval:', uval

            !----------------------------------------------------------------------
            !  Data on v points

            i=iflag
            j=jflag
            k=kflag

            if( x3d.lt.xh(i) )then
              i=i-1
            endif
            if( .not. terrain_flag )then
              if( z3d.lt.zh(iflag,jflag,k) )then
                k=k-1
              endif
            else
              if( sig3d.lt.sigma(k) )then
                k=k-1
              endif
            endif

            rx = ( x3d-xh(i) )/( xh(i+1)-xh(i) )
            ry = ( y3d-yf(j) )/( yf(j+1)-yf(j) )
            rz = rzs

            ! saveit:
            rxs = rx
            ryv = ry

            w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
            w2=rx*(1.0-ry)*(1.0-rz)
            w3=(1.0-rx)*ry*(1.0-rz)
            w4=(1.0-rx)*(1.0-ry)*rz
            w5=rx*(1.0-ry)*rz
            w6=(1.0-rx)*ry*rz
            w7=rx*ry*(1.0-rz)
            w8=rx*ry*rz

            IF(debug)THEN
                wsum = w1+w2+w3+w4+w5+w6+w7+w8
                if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
                    ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
                    rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
                    wsum.le.0.99999 .or.                  &
                    wsum.ge.1.00001 .or.                  &
                    i.lt.-1 .or. i.gt.(ni+1)   .or.       &
                    j.lt.0 .or. j.gt.(nj+1)   .or.        &
                    k.lt.0 .or. k.gt.nk                   )then
                  print *
                  print *,'  23333b: '
                  print *,'  np          = ',np
                  print *,'  myid,i,j,k  = ',myid,i,j,k
                  print *,'  rx,ry,rz    = ',rx,ry,rz
                  print *,'  wsum        = ',wsum
                  print *,'  xh1,x3d,xh2 = ',xh(i),x3d,xh(i+1)
                  print *,'  yf1,y3d,yh2 = ',yf(j),y3d,yf(j+1)
                  print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),z3d,zh(iflag,jflag,k+1)
                  print *
                  call stopcm1
                endif
            ENDIF

            call tri_interp(ni,nj+1,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,va,vval)
            !print *, ' hail.F, tri_interp vval:', vval

            !----------------------------------------------------------------------
            !  Data on w points

            i=iflag
            j=jflag
            k=kflag

            if( x3d.lt.xh(i) )then
              i=i-1
            endif
            if( y3d.lt.yh(j) )then
              j=j-1
            endif

            rx = rxs
            ry = rys
            if( .not. terrain_flag )then
              rz = ( z3d-zf(iflag,jflag,k) )/( zf(iflag,jflag,k+1)-zf(iflag,jflag,k) )
            else
              rz = ( sig3d-sigmaf(k) )/( sigmaf(k+1)-sigmaf(k) )
            endif

            ! saveit:
            rzw = rz

            w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
            w2=rx*(1.0-ry)*(1.0-rz)
            w3=(1.0-rx)*ry*(1.0-rz)
            w4=(1.0-rx)*(1.0-ry)*rz
            w5=rx*(1.0-ry)*rz
            w6=(1.0-rx)*ry*rz
            w7=rx*ry*(1.0-rz)
            w8=rx*ry*rz

            IF(debug)THEN
                wsum = w1+w2+w3+w4+w5+w6+w7+w8
                if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
                    ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
                    rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
                    wsum.le.0.99999 .or.                  &
                    wsum.ge.1.00001 .or.                  &
                    i.lt.-1 .or. i.gt.(ni+1)   .or.       &
                    j.lt.-1 .or. j.gt.(nj+1)   .or.       &
                    k.lt.1 .or. k.gt.nk                   )then
                  print *
                  print *,'  43333a: '
                  print *,'  np          = ',np
                  print *,'  myid,i,j,k  = ',myid,i,j,k
                  print *,'  rx,ry,rz    = ',rx,ry,rz
                  print *,'  wsum        = ',wsum
                  print *,'  xh1,x3d,xh2 = ',xh(i),x3d,xh(i+1)
                  print *,'  yh1,y3d,yh2 = ',yh(j),y3d,yh(j+1)
                  print *,'  zh1,z3d,zh2 = ',zf(iflag,jflag,k),z3d,zf(iflag,jflag,k+1)
                  print *
                  call stopcm1
                endif
            ENDIF

            call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,wa,wval)
            !print *, ' HAILCAST tri_interp np, w: ', np, wval
            !print *, ' HAILCAST np, wa values: ', np, wa(iflag,jflag,k)
            !if (np.eq.1282828) print *, ' hail.F, i, j, k: ', i, j, k
            !if (np.eq.1282828) print *, ' hail.F, xh(i), yh(j): ', xh(i), yh(j)

            if( terrain_flag )then
                  call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,wa,sigdot)
                  call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,zs,zsp)
                  z3d = zsp + sig3d*((zt-zsp)*rzt)
            endif

            !----------------------------------------------------------------------
            !  Data on scalar points

            i=iflag
            j=jflag
            k=kflag

            if(x3d.lt.xh(i) )then
              i=i-1
            endif
            if(y3d.lt.yh(j) )then
              j=j-1
            endif
            if( .not. terrain_flag )then
              if(z3d.lt.zh(iflag,jflag,k) )then
                k=k-1
              endif
            else
              if( sig3d.lt.sigma(k) )then
                k=k-1
              endif
            endif

            !k can be 0 here, because rz ignores it.
            rx = rxs
            ry = rys
            rz = rzs

            w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
            w2=rx*(1.0-ry)*(1.0-rz)
            w3=(1.0-rx)*ry*(1.0-rz)
            w4=(1.0-rx)*(1.0-ry)*rz
            w5=rx*(1.0-ry)*rz
            w6=(1.0-rx)*ry*rz
            w7=rx*ry*(1.0-rz)
            w8=rx*ry*rz

            IF(debug)THEN
            wsum = w1+w2+w3+w4+w5+w6+w7+w8
            if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
                ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
                rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
                wsum.le.0.99999 .or.                  &
                wsum.ge.1.00001 .or.                  &
                i.lt.-1 .or. i.gt.ni   .or.           &
                j.lt.-1 .or. j.gt.nj   .or.           &
                k.lt.0 .or. k.gt.nk                   )then
              print *
              print *,'  15558: '
              print *,'  np          = ',np
              print *,'  myid,i,j,k  = ',myid,i,j,k
              print *,'  rx,ry,rz    = ',rx,ry,rz
              print *,'  wsum        = ',wsum
              print *
              call stopcm1
            endif
            ENDIF

            call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,temp,tempval)
            call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,rho,rhoval)



            !----------------------------------------------------------------------
            !  uv for parcels below lowest model level:
            !   (note - non wind parcels were already set above 
            !     via extrapolation from the free-slip conditions)

            IF( bbc.eq.3 )THEN
                ! semi-slip lower boundary condition:
                if( z3d.lt.zh(1,1,1) )then
                  ! re-calculate velocities if parcel is below lowest model level:
                  !------
                  ! u at lowest model level:
                  i=iflag
                  j=jflag
                  if( y3d.lt.yh(j) )then
                    j=j-1
                  endif
                  call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 1, 0, 1, 0,ua(ib,jb,1),uval)
                  !------
                  ! v at lowest model level:
                  i=iflag
                  j=jflag
                  if( x3d.lt.xh(i) )then
                    i=i-1
                  endif
                  call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 1, 0, 1,va(ib,jb,1),vval)
                  !------
                  ! z0:
                  i=iflag
                  j=jflag
                  if( x3d.lt.xh(i) )then
                    i=i-1
                  endif
                  if( y3d.lt.yh(j) )then
                    j=j-1
                  endif
                  call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,znt,z0)
                  !------
                  ! get u,v from (neutral) log-layer equation:
                  rznt = 1.0/z0
                  var = alog((z3d+z0)*rznt)/alog((zh(1,1,1)+z0)*rznt)
                  if( imove.eq.1 )then
                    uval = (uval+umove)*var - umove
                    vval = (vval+vmove)*var - vmove
                  else
                    uval = uval*var
                    vval = vval*var
                  endif
                endif
            ENDIF


            !-----------------------------------------------------
            !  Account for terminal velocity of hailstone in 
            !  vertical motion calculation
            !  Note TERMINL returns a positive value for downward
            !-----------------------------------------------------
            ! if (np.eq.313793) THEN
            !     print *, ' HAILCAST before TERMINL, np,wval: ',np,wval
            ! endif
            CALL TERMINL(rhoval, hdense3d, hd3d, haildata(np,hailvt), tempval,&
                         oblate_vt, a_shedd)
            if( terrain_flag )then
                !Convert m/s to sigma units....not sure about this
                sigdot = sigdot - haildata(np,hailvt) / ((zt-zsp)*rzt)
            else
                wval = wval - haildata(np,hailvt)
            endif
            ! if (np.eq.313793) THEN
            !     print *, ' HAILCAST after TERMINL, np,wval,vt: ',np,wval,haildata(np,hailvt)
            ! endif


            !-----------------------------------------------------
            !  Update parcel positions:
            !-----------------------------------------------------



            ! RK2 scheme:
            IF(nrkp.eq.1)THEN
                IF(nx.eq.1)THEN
                  x3d=0.0
                ELSE
                  x3d=haildata(np,hailx)+dt*uval
                ENDIF
                IF(axisymm.eq.1.or.ny.eq.1)THEN
                  y3d=0.0
                ELSE
                  y3d=haildata(np,haily)+dt*vval
                ENDIF
                if( terrain_flag )then
                  !again, sig code just copied without comment from parcel code
                  sig3d = haildata(np,hailsig) + dt*sigdot
                  sig1 = sigdot
                  IF( sig3d.lt.0.0 )THEN
                    sig3d=1.0e-6
                  ENDIF
                else
                  z3d = haildata(np,hailz)+dt*wval
                  ! 200129 RAS - needed to add check for negative z values
                  !  during first rk loop as well.  Hailstones move the fastest
                  !  near the ground.
                  IF( z3d.lt.0.0 )THEN
                    z3d = 0.0
                  ENDIF
                  ww1=wval
                endif
                uu1=uval
                vv1=vval
            ELSE
                IF(nx.eq.1)THEN
                  x3d=0.0
                ELSE
                  x3d=haildata(np,hailx)+dt2*(uu1+uval)
                ENDIF
                IF(axisymm.eq.1.or.ny.eq.1)THEN
                  y3d=0.0
                ELSE
                  y3d=haildata(np,haily)+dt2*(vv1+vval)
                ENDIF
                if( terrain_flag )then
                  sig3d = haildata(np,hailsig) + dt2*(sig1+sigdot)
                  !This sig3d check code below has not been tested or confirmed useful,
                  ! and was merely copied from the parcel code.
                  IF( sig3d.lt.0.0 )THEN
                    print *,'  hailstone is below surface:  np,x3d,y3d,sig3d = ',np,x3d,y3d,sig3d
                    sig3d=1.0e-6
                    !sig3d = 0.0
                  ENDIF
                  sig3d=min(sig3d,maxz)
                else
                  z3d = haildata(np,hailz)+dt2*(ww1+wval)
                  IF( z3d.lt.0.0 )THEN
                    !print *,'  hailstone is below surface:  np,x3d,y3d,z3d = ',np,x3d,y3d,z3d
                    !z3d=1.0e-6
                    z3d = 0.0
                  ENDIF
                  z3d=min(z3d,maxz)
                endif
            ENDIF

          !if (np.eq.313793) write(*,'(a28,i7,6f8.1)'), ' HAILCAST end rkloop: ', &
          ! nrkp, x3d, y3d, z3d, uval, vval, wval

          ENDDO  rkloop

          !-----------------------------------------------------
          !  Account for boundary conditions (if necessary)
          !-----------------------------------------------------

          ! From parcel code: New for cm1r17:  if parcel exits domain,
          ! just assume periodic lateral boundary conditions
          ! (no matter what actual settings are for wbc,ebc,sbc,nbc)
          
          ! For HAILCAST, we don't want hailstones periodically 
          ! circling the domain.Follow the convention used with Z
          ! above instead.

          if(x3d.lt.minx)then
              !x3d=x3d+(maxx-minx)
              x3d = minx
          endif
          if(x3d.gt.maxx)then
              !x3d=x3d-(maxx-minx)
              x3d=maxx
          endif

          if( (y3d.gt.maxy).and.(axisymm.ne.1).and.(ny.ne.1) )then
              !y3d=y3d-(maxy-miny)
              y3d = maxy
          endif
          if( (y3d.lt.miny).and.(axisymm.ne.1).and.(ny.ne.1) )then
              !y3d=y3d+(maxy-miny)
              y3d = miny
          endif

          haildata(np,hailx)=x3d
          haildata(np,haily)=y3d
          if( .not. terrain_flag )then
              haildata(np,hailz)=z3d
          else
              haildata(np,hailsig)=sig3d
          endif

          !if (np.eq.1282828) write(*,'(a28,i7,4f8.1)'), ' HAILCAST after checkk: ', &
          !  np, x3d, y3d, z3d, dt2
          !-----------------------------------------------------
          !  Now that we're done moving the hailstone around,
          !  let's let it interact with the environment.
          !-----------------------------------------------------

          ! First, store our u, v, w values
          haildata(np,hailu) = 0.5*(uu1+uval)
          haildata(np,hailv) = 0.5*(vv1+vval)
          !In our output, we want to know the w the hailstone is experiencing,
          ! which does not include its terminal velocity
          if( terrain_flag ) then
              haildata(np,hailw) = 0.5*(sig1+sigdot) + haildata(np,hailvt)
          else
              haildata(np,hailw) = 0.5*(ww1+wval) + haildata(np,hailvt)
          endif
          
          !if((myid.eq.0).and.(np.ge.20000).and.(np.lt.20100))THEN
          !  print *, ' HAILCAST, start interact, np,hd w: ',np,haildata(np,hailw)      
          !endif
          !----------------------------------------------------------------------
          !  Data on scalar points

          i=iflag
          j=jflag
          k=kflag

          if(x3d.lt.xh(i) )then
            i=i-1
          endif
          if(y3d.lt.yh(j) )then
            j=j-1
          endif
          if( .not. terrain_flag )then
            if(z3d.lt.zh(iflag,jflag,k) )then
              k=k-1
            endif
          else
            if( sig3d.lt.sigma(k) )then
              k=k-1
            endif
          endif

          rx = rxs
          ry = rys
          rz = rzs

          w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
          w2=rx*(1.0-ry)*(1.0-rz)
          w3=(1.0-rx)*ry*(1.0-rz)
          w4=(1.0-rx)*(1.0-ry)*rz
          w5=rx*(1.0-ry)*rz
          w6=(1.0-rx)*ry*rz
          w7=rx*ry*(1.0-rz)
          w8=rx*ry*rz

          IF(debug)THEN
          wsum = w1+w2+w3+w4+w5+w6+w7+w8
          if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
              ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
              rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
              wsum.le.0.99999 .or.                  &
              wsum.ge.1.00001 .or.                  &
              i.lt.-1 .or. i.gt.ni   .or.           &
              j.lt.-1 .or. j.gt.nj   .or.           &
              k.lt.0 .or. k.gt.nk                   )then
            print *
            print *,'  15558: '
            print *,'  np          = ',np
            print *,'  myid,i,j,k  = ',myid,i,j,k
            print *,'  rx,ry,rz    = ',rx,ry,rz
            print *,'  wsum        = ',wsum
            print *
            call stopcm1
          endif
          ENDIF

          call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,temp,tempval)
          call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,rho,rhoval)
          call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,prs,prsval)
          call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qv,qvval)
          call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qia,qiaval)
          call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qwa,qwaval)

!           if (np.eq.51187)  write(*,'(a46,i7,6f12.4)'),&
!               'HAILCAST before hailtraj,np,t,rho,p,qv,qi,qw: ',&
!               np,tempval,rhoval,prsval,qvval,qiaval,qwaval

          !-----------------------------------------------------
          !  Call hailstone driver and store the results.
          !  Note - still need to update even if there isn't water so
          !   the falling hailstone can melt.
          !  However, we don't want to call this if the hailstone 
          !     is at any of the edges of the model domain, so check first.
          !-----------------------------------------------------
          
          IF ( terrain_flag ) THEN
              IF ( (x3d.gt.minx) .AND. (x3d.lt.maxx) .AND. &
                   (y3d.gt.miny) .AND. (y3d.lt.maxy) .AND. &
                   (sig3d.gt.0.0) .AND. (sig3d.lt.maxz) ) THEN
 
                   CALL hailstone_trajectory(tempval,prsval,rhoval,   &
                       qvval,qiaval,qwaval,dt,haildata(np,:),np,      &
                       a_shedd)                  
              ENDIF
          ELSE
              IF ( (x3d.gt.minx) .AND. (x3d.lt.maxx) .AND. &
                   (y3d.gt.miny) .AND. (y3d.lt.maxy) .AND. &
                   (z3d.gt.0.0)  .AND. (z3d.lt.maxz) ) THEN

                    ! if (np.eq.313793) write(*,*),' '
                    ! if (np.eq.313793) write(*,'(a33,3(1x,f12.4))'),&
                    !     'HAILCAST before traj d,fw,dns:',&
                    !     haildata(np,haild), &
                    !     haildata(np,hailfw),haildata(np,haildense)

                   CALL hailstone_trajectory(tempval,prsval,rhoval,   &
                       qvval,qiaval,qwaval,dt,haildata(np,:),np,      &
                       a_shedd)                  

                    ! if (np.eq.313793) write(*,'(a33,3(1x,f12.4))'),&
                    !     'HAILCAST after traj d,fw,dns:',&
                    !     haildata(np,haild), &
                    !     haildata(np,hailfw),haildata(np,haildense)
                    ! if (np.eq.313793) write(*,*),' '

                      !if (np.eq.313793) write(*,'(a36,i7,4f12.4)'),'HAILCAST after traj, np,d,ts,vt,de: ',&
!                       np,haildata(np,haild),haildata(np,hailts),haildata(np,hailvt),&
!                        haildata(np,haildense)


              ENDIF
          ENDIF

#ifdef MPI

      ELSE

            ! set to really small number (so we can use the allreduce command below)
            haildata(np,hailx) = -1.0e30
            haildata(np,haily) = -1.0e30
            if( .not. terrain_flag )then
              haildata(np,hailz) = -1.0e30
            else
              haildata(np,hailsig) = -1.0e30
            endif
            haildata(np,hailu) = -1.0e30
            haildata(np,hailv) = -1.0e30
            haildata(np,hailw) = -1.0e30
            haildata(np,haild) = -1.0e30
            haildata(np,haildense) = -1.0e30
            haildata(np,haildice) = -1.0e30
            haildata(np,hailts) = -1.0e30
            haildata(np,hailvt) = -1.0e30
            haildata(np,hailitype) = -1.0e30
            haildata(np,hailqice) = -1.0e30
            haildata(np,hailqliq) = -1.0e30
            haildata(np,hailtc) = -1.0e30
#endif

      ENDIF  myparcel

    ENDDO  nploop

!----------------------------------------------------------------------
!  communicate data  (for MPI runs)

#ifdef MPI
    if( .not. terrain_flag )then
        call MPI_ALLREDUCE(MPI_IN_PLACE,haildata(1,1),16*nhailtrajs,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
    else
        call MPI_ALLREDUCE(MPI_IN_PLACE,haildata(1,1),11*nhailtrajs,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE,haildata(1,hailsig),nhailtrajs,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
    endif
#endif

!----------------------------------------------------------------------
!  get height ASL:

    if( terrain_flag )then
        call gethailtrajzs(xh,uh,ruh,xf,yh,vh,rvh,yf,zs,haildata)
        DO np=1,nhailtrajs
          ! get z from sigma:
          ! (see Section 3 of "The governing equations for CM1", 
          !  http://www2.mmm.ucar.edu/people/bryan/cm1/cm1_equations.pdf)

          ! Think this is right...
          haildata(np,hailz) = haildata(np,hailzs) + &
            haildata(np,hailsig)*((zt-haildata(np,hailzs))*rzt)
        ENDDO
    endif

!----------------------------------------------------------------------

      !WRITE (*,*) 'HAILCAST, MAXVAL(w), w(314): ', MAXVAL(haildata(:,hailw)), haildata(314,hailw)

      end subroutine hailcast_driver


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine tri_interp(iz,jz,kz,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,s,haildata)
      use input
      implicit none

      integer :: iz,jz,kz,i,j,k
      real :: w1,w2,w3,w4,w5,w6,w7,w8
      real, dimension(1-ngxy:iz+ngxy,1-ngxy:jz+ngxy,1-ngz:kz+ngz) :: s
      real :: haildata

      haildata=s(i  ,j  ,k  )*w1    &
           +s(i+1,j  ,k  )*w2    &
           +s(i  ,j+1,k  )*w3    &
           +s(i  ,j  ,k+1)*w4    &
           +s(i+1,j  ,k+1)*w5    &
           +s(i  ,j+1,k+1)*w6    &
           +s(i+1,j+1,k  )*w7    &
           +s(i+1,j+1,k+1)*w8

      !  WRITE(*,'(a14,8f12.4)') 'tri_interp: ', s(i  ,j  ,k  ), &
      !      s(i+1,j  ,k  ),    &
      !      s(i  ,j+1,k  ),    &
      !      s(i  ,j  ,k+1),    &
      !      s(i+1,j  ,k+1),    &
      !      s(i  ,j+1,k+1),   &
      !      s(i+1,j+1,k  ),    &
      !      s(i+1,j+1,k+1)


      end subroutine tri_interp


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine get2d(i,j,x3d,y3d,xh,xf,yh,yf,xs,ys,is,js,s,sval)
    use input
    implicit none

    integer, intent(in) :: i,j
    real, intent(in) :: x3d,y3d
    real, intent(in), dimension(ib:ie) :: xh
    real, intent(in), dimension(ib:ie+1) :: xf
    real, intent(in), dimension(jb:je) :: yh
    real, intent(in), dimension(jb:je+1) :: yf

    ! 0 = scalar point
    ! 1 = velocity point
    integer, intent(in) :: xs,ys
    integer, intent(in) :: is,js

    real, intent(in), dimension(ib:ie+is,jb:je+js) :: s
    real, intent(out) :: sval

    real :: wg1,wg2,wg3,wg4
    real :: x13,x23,x33,x43
    real :: w1,w2,w3,w7,rx,ry,rz

    logical, parameter :: debug = .false.

!-----------------------------------------------------------------------
      ! tri-linear interp:

      IF(xs.eq.1)THEN
        rx = ( x3d-xf(i) )/( xf(i+1)-xf(i) )
      ELSE
        rx = ( x3d-xh(i) )/( xh(i+1)-xh(i) )
      ENDIF

      IF(ys.eq.1)THEN
        ry = ( y3d-yf(j) )/( yf(j+1)-yf(j) )
      ELSE
        ry = ( y3d-yh(j) )/( yh(j+1)-yh(j) )
      ENDIF

        w1=(1.0-rx)*(1.0-ry)
        w2=rx*(1.0-ry)
        w3=(1.0-rx)*ry
        w7=rx*ry

      IF( debug )THEN
        if( rx.lt.-0.000001 .or. rx.gt.1.000001 .or.        &
            ry.lt.-0.000001 .or. ry.gt.1.000001 .or.        &
            (w1+w2+w3+w7).lt.0.999999 .or.  &
            (w1+w2+w3+w7).gt.1.000001       &
          )then
          print *,'  x3d,y3d     = ',x3d,y3d
          print *,'  i,j         = ',i,j
          print *,'  rx,ry       = ',rx,ry
          print *,'  w1,w2,w3,w7 = ',w1,w2,w3,w7
          print *,'  w1+w2+w3+w7 = ',w1+w2+w3+w7
          print *,' 22346 '
          call stopcm1
        endif
      ENDIF

      sval =s(i  ,j  )*w1    &
           +s(i+1,j  )*w2    &
           +s(i  ,j+1)*w3    &
           +s(i+1,j+1)*w7

!-----------------------------------------------------------------------

    end subroutine get2d


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine gethailtrajzs(xh,uh,ruh,xf,yh,vh,rvh,yf,zs,haildata)
      use input
#ifdef MPI
      use mpi
#endif
      implicit none

      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(ib:ie,jb:je) :: zs
      real, intent(inout), dimension(nhailtrajs,nhvals) :: haildata

      integer :: i,j,iflag,jflag,np
      real :: x3d,y3d

    zsnploop:  &
    DO np=1,nhailtrajs

      x3d = haildata(np,hailx)
      y3d = haildata(np,haily)

      iflag = -100
      jflag = -100

  ! cm1r19:  skip if we already know this processor doesnt have this parcel
  zshaveit1:  &
  IF( x3d.ge.xf(1) .and. x3d.le.xf(ni+1) .and.  &
      y3d.ge.yf(1) .and. y3d.le.yf(nj+1) )THEN

    IF(nx.eq.1)THEN
      iflag = 1
    ELSE
      ! cm1r19:
      i = ni+1
      do while( iflag.lt.0 .and. i.gt.1 )
        i = i-1
        if( x3d.ge.xf(i) .and. x3d.le.xf(i+1) )then
          iflag = i
        endif
      enddo
    ENDIF

    IF(axisymm.eq.1.or.ny.eq.1)THEN
      jflag = 1
    ELSE
      ! cm1r19:
      j = nj+1
      do while( jflag.lt.0 .and. j.gt.1 )
        j = j-1
        if( y3d.ge.yf(j) .and. y3d.le.yf(j+1) )then
          jflag = j
        endif
      enddo
    ENDIF

  ENDIF  zshaveit1

#ifdef MPI
      ! check for conflict:
    IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
        (jflag.ge.1.and.jflag.le.nj) )THEN
      IF( iflag.eq.ni .and. haildata(np,hailx).eq.xf(iflag+1) .and. nodex.gt.1 .and.  myi.ne.nodex ) iflag = -1
      IF( jflag.eq.nj .and. haildata(np,haily).eq.yf(jflag+1) .and. nodey.gt.1 .and.  myj.ne.nodey ) jflag = -1
    ENDIF
#endif

      zsmyparcel:  IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
                       (jflag.ge.1.and.jflag.le.nj) )THEN

        i=iflag
        j=jflag

        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( y3d.lt.yh(j) )then
          j=j-1
        endif

        call get2d(i,j,x3d,y3d,xh,xf,yh,yf, 0, 0, 0, 0,zs,haildata(np,hailzs))

      ELSE  zsmyparcel

        haildata(np,hailzs) = -1.0e30

      ENDIF  zsmyparcel

    ENDDO  zsnploop

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE,haildata(1,hailzs) ,nhailtrajs,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
#endif

      end subroutine gethailtrajzs


      subroutine setup_hailcast_vars(name_hail,desc_hail,unit_hail,qname,tdef)
      use input
      implicit none

      character(len=40), intent(inout), dimension(maxvars) :: name_hail,desc_hail,unit_hail
      character(len=3), intent(in), dimension(maxq) :: qname
      character(len=15), intent(inout) :: tdef

      integer :: n,n2
      character(len=8) :: text1
      character(len=30) :: text2

      hail_out = 0

      hail_out = hail_out+1
      name_hail(hail_out) = 'x'
      desc_hail(hail_out) = 'x position'
      unit_hail(hail_out) = 'm'

      hail_out = hail_out+1
      name_hail(hail_out) = 'y'
      desc_hail(hail_out) = 'y position'
      unit_hail(hail_out) = 'm'

      hail_out = hail_out+1
      name_hail(hail_out) = 'z'
      desc_hail(hail_out) = 'z position (above sea level)'
      unit_hail(hail_out) = 'm'

      hail_out = hail_out+1
      name_hail(hail_out) = 'u'
      desc_hail(hail_out) = 'u velocity'
      unit_hail(hail_out) = 'm/s'

      hail_out = hail_out+1
      name_hail(hail_out) = 'v'
      desc_hail(hail_out) = 'v velocity'
      unit_hail(hail_out) = 'm/s'

      hail_out = hail_out+1
      name_hail(hail_out) = 'w'
      desc_hail(hail_out) = 'w velocity'
      unit_hail(hail_out) = 'm/s'

      hail_out = hail_out+1
      name_hail(hail_out) = 'd'
      desc_hail(hail_out) = 'hailstone diameter'
      unit_hail(hail_out) = 'm'
      
      hail_out = hail_out+1
      name_hail(hail_out) = 'dense'
      desc_hail(hail_out) = 'hailstone density'
      unit_hail(hail_out) = 'kg/m3'

      hail_out = hail_out+1
      name_hail(hail_out) = 'tv'
      desc_hail(hail_out) = 'hailstone terminal velocity'
      unit_hail(hail_out) = 'm/s'

      hail_out = hail_out+1
      name_hail(hail_out) = 'ts'
      desc_hail(hail_out) = 'hailstone temperature'
      unit_hail(hail_out) = 'K'

      hail_out = hail_out+1
      name_hail(hail_out) = 'fw'
      desc_hail(hail_out) = 'hailstone water fraction'
      unit_hail(hail_out) = 'unitless'

      hail_out = hail_out+1
      name_hail(hail_out) = 'itype'
      desc_hail(hail_out) = 'wet or dry growth mode'
      unit_hail(hail_out) = 'unitless'

      if( haildice.ge.1 )then
        hail_out = hail_out+1
        name_hail(hail_out) = 'dice'
        desc_hail(hail_out) = 'hailstone ice-only diameter'
        unit_hail(hail_out) = 'm'
      endif

      if( hailqice.ge.1 )then
        hail_out = hail_out+1
        name_hail(hail_out) = 'qice'
        desc_hail(hail_out) = 'in-cloud ice (qi+qs) at hailstone loc.'
        unit_hail(hail_out) = 'kg/kg'
      endif

      if( hailqliq.ge.1 )then
        hail_out = hail_out+1
        name_hail(hail_out) = 'qliq'
        desc_hail(hail_out) = 'in-cloud liquid (qc) at hailstone loc.'
        unit_hail(hail_out) = 'mm'
      endif

      if( hailtc.ge.1 )then
        hail_out = hail_out+1
        name_hail(hail_out) = 'tc'
        desc_hail(hail_out) = 'in-cloud temp. at hailstone loc.'
        unit_hail(hail_out) = 'K'
      endif

      if( hailzs.ge.1 )then
        hail_out = hail_out+1
        name_hail(hail_out) = 'zs'
        desc_hail(hail_out) = 'terrain height'
        unit_hail(hail_out) = 'm'
      endif

      if( hailsig.ge.1 )then
        hail_out = hail_out+1
        name_hail(hail_out) = 'sigma'
        desc_hail(hail_out) = 'sigma (nondimensional height)'
        unit_hail(hail_out) = 'nondimensional'
      endif

      if (myid.eq.0) print *, 'setup hail_out: ', hail_out

!-----------------------------------------------------------------------

      if( hail_out.gt.0 .and. output_format.eq.1 )then
        ! write GrADS descriptor file:
        call write_hailctl(name_hail,desc_hail,unit_hail,tdef)
      endif

!-----------------------------------------------------------------------

      end subroutine setup_hailcast_vars

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write_hailctl(name_hail,desc_hail,unit_hail,tdef)
      use input
      use constants , only : grads_undef
      implicit none

      character(len=40), intent(inout), dimension(maxvars) :: name_hail,desc_hail,unit_hail
      character(len=15), intent(inout) :: tdef

      integer :: n,nn
      character(len=16) :: a16

      !---------------------------------------------------------------
      ! This subroutine writes the GrADS descriptor file for hailstones
      !---------------------------------------------------------------

    idcheck:  &
    IF( myid.eq.0 )THEN

      string(totlen+1:totlen+22) = '_haildata.ctl            '
      if(dowr) write(outfile,*) string
      open(unit=50,file=string,status='unknown')

      sstring(baselen+1:baselen+1+12) = '_haildata.dat  '

      write(50,401) sstring
      write(50,402)
      write(50,403) grads_undef
      write(50,404) nhailtrajs
      write(50,405)
      write(50,406)
      if( hailfrq.gt.0 )then
        write(50,407) 1+int(timax/hailfrq),tdef,max(1,int(hailfrq/60.0))
      else
        write(50,407) 1000000000,tdef,max(1,int(hailfrq/60.0))
      endif

      write(50,408) hail_out

      DO n = 1 , hail_out
        a16 = '                '
        nn = len(trim(unit_hail(n)))
        write(a16(2:15),214) unit_hail(n)
        write(a16(1:1),201 )       '('
        write(a16(nn+2:nn+2),201 ) ')'
        write(50,409) name_hail(n),desc_hail(n),a16
      ENDDO

      write(50,410)

      close(unit=50)

    ENDIF  idcheck

201   format(a1)
214   format(a14)

401   format('dset ^',a70)
402   format('title cm1r19 output, hail data')
403   format('undef ',f10.1)
404   format('xdef ',i10,' linear 1 1')
405   format('ydef          1 linear 1 1')
406   format('zdef          1 linear 1 1')
407   format('tdef ',i10,' linear ',a15,' ',i5,'MN')
408   format('vars ',i6)
409   format(a12,' 1 99 ',a40,1x,a16)
410   format('endvars')

      end subroutine write_hailctl


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine hailcast_write(hrec,rtime,qname,name_hail,desc_hail,unit_hail,haildata,hailloc)
      use input
#ifdef NETCDF
      use writeout_nc_module, only : writehaildata_nc
#endif
      implicit none

      integer, intent(inout) :: hrec
      real, intent(in) :: rtime
      character(len=3), intent(in), dimension(maxq) :: qname
      character(len=40), intent(in), dimension(maxvars) :: name_hail,desc_hail,unit_hail
      real, intent(in), dimension(nhailtrajs,nhvals) :: haildata
      real, intent(inout), dimension(nhailtrajs,3) :: hailloc

      integer :: n,np

!----------------------------------------------------------------------
!  write out data

    IF(myid.eq.0)THEN

      IF(output_format.eq.1)THEN
        ! GrADS format:

        string(totlen+1:totlen+22) = '_haildata.dat         '
        if(dowr) write(outfile,*) string
        open(unit=61,file=string,form='unformatted',access='direct',   &
             recl=4*nhvals*nhailtrajs,status='unknown')

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  hdata hrec = ',hrec

        write(61,rec=hrec) ((haildata(np,n),np=1,nhailtrajs),n=1,nhvals)

        close(unit=61)

#ifdef NETCDF
      ELSEIF(output_format.eq.2)THEN

        call     writehaildata_nc(hrec,rtime,qname,name_hail,desc_hail,unit_hail,&
                    haildata,hailloc(1,1))

#endif
      ENDIF
      if(dowr) write(outfile,*)

    ENDIF   ! endif for myid=0

      return
      end subroutine hailcast_write




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  SUBROUTINE hailstone_trajectory ( TCA, PA, rho, RA, RIA, RWA, dt,  &
                                haildata1d, np, a_shedd)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!! Hailstone driver, adapted from hailstone subroutine in HAILCAST
!!!!  Inputs:
!!!!     TCA          temperature (K) 
!!!!     PA           total pressure (Pa)
!!!!     rho          density (kg/m3)
!!!!     RA           vapor mixing ratio (kg/kg)
!!!!     RWA          cloud water mixing ratio (kg/kg)
!!!!     RIA          cloud ice and snow mixing ratio (kg/kg)
!!!!
!!!!  Input/Output:
!!!!     haildata1d:  array of hailstone characteristics. 
!!!!      including:  hail diameter in ***m***
!!!!
!!!!  13 Aug 2013 .................................Becky Adams-Selin AER
!!!!     adapted from hailstone subroutine in SPC's HAILCAST
!!!!  18 Mar 2014 .................................Becky Adams-Selin AER
!!!!     added variable rime layer density, per Ziegler et al. (1983)
!!!!  4 Jun 2014 ..................................Becky Adams-Selin AER
!!!!     removed initial embryo size dependency on microphysic scheme
!!!!  5 Jun 2014 ..................................Becky Adams-Selin AER
!!!!     used smaller initial embryo sizes
!!!!  25 Jun 2015..................................Becky Adams-Selin AER
!!!!     Significant revamping.  Fixed units bug in HEATBUD that caused
!!!!     hailstone temperature instabilities.  Similar issue fixed in BREAKUP
!!!!     subroutine.  Removed graupel from ice content.  Changed initial
!!!!     embryo size and location to better match literature.  Added
!!!!     enhanced melting when hailstone collides with liquid water
!!!!     in regions above freezing.  Final diameter returned is ice diameter
!!!!     only. Added hailstone temperature averaging over previous timesteps 
!!!!     to decrease initial temperature instability at small embyro diameters.  
!!!!  3 Sep 2015...................................Becky Adams-Selin AER
!!!!    Insert embryos at -13C; interpret pressure and other variables to
!!!!    that exact temperature level.
!!!! 16 Nov 2015...................................Becky Adams-Selin AER
!!!!     Hailstone travels horizontally through updraft instead of being
!!!!     locked in the center.
!!!!  9 Jul 2016...................................Becky Adams-Selin AER
!!!!     Uses an adiabatic liquid cloud water profile generated from
!!!!     the saturation vapor pressure using the model temperature
!!!!     profile.
!!!! 15 Apr 2017...................................Becky Adams-Selin AER
!!!!     Converted to a trajectory model for use in 3D steady-state simulations.
!!!!     Uses cloud water directly from netcdf file, no adiabatic profile.
!!!! 15 Mar 2018...................................Becky Adams-Selin AER
!!!!     Converted to allow multiple trajectories to be run at once
!!!! 2 Oct 2019....................................Becky Adams-Selin AER
!!!!     Added time-dependent melting outside of cloud.
!!!! 9 Oct 2019....................................Becky Adams-Selin AER
!!!!     A few melting bug fixes.
!!!! 31 Oct 2019...................................Becky Adams-Selin AER
!!!!     Incorporated into CM1 as a time-dependent trajectory model.
!!!!     Removed advection and interpolation code as that happens in 
!!!!     hailcast_driver now. Subroutine only advances hailstone 
!!!!     interactions for 1 timestep (dt), so removed temperature
!!!!     differencing.
!!!! 18 Mar 2022...................................Becky Adams-Selin AER
!!!!     Updates to fix units in the hail melt calculation during wet
!!!!     growth with air temperatures > 0C. Fix dynamic viscosity calculation 
!!!!     exponent error, units error in melting subroutine, and SI conversion 
!!!!     inside HEATBUD.
!!!! 20 Apr 2022...................................Becky Adams-Selin AER
!!!!     Bug fix with sign of AFACTOR from HP85.
!!!! 3 Jan 2024...................................Becky Adams-Selin AER
!!!!     Add ability to model oblate hailstones!
!!!!
!!!! There are several options set in namelist &param17 that determine
!!!! what physics is used herein. Consult the readme.md file for an 
!!!! explanation.
!!!!
!!!! See Adams-Selin and Ziegler 2016, MWR; Adams-Selin 2023 MWR; and 
!!!!  Pounds et al. 2024a, b JAS for further documentation.
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use input

    IMPLICIT NONE

    REAL, INTENT(IN   ) ::                                 TCA  & ! temperature (K)
                                              ,            rho  & ! density (kg/m3)
                                              ,             PA  & ! pressure (Pa)
                                              ,             RA  & ! vapor mixing ratio (kg/kg)
                                              ,            RIA  & ! kg/kg
                                              ,            RWA    ! kg/kg
    real, intent(in) :: dt  !time step 
    
    !OUTPUT:  hailstone information for *this* hailstone
    REAL, DIMENSION (nhvals), INTENT(INOUT) :: haildata1d
    !see param.F for which variables are defined to be where

    integer, intent(in) :: np  !hailstone trajectory ID, for debugging

    !Local variables
    REAL g                     ! gravity (m/s)
    REAL r_d                   ! constant
    REAL pi                    ! yup
    !hailstone parameters
    REAL*8 D, D_ICE            ! hail diameter (m)
    REAL VT                    ! terminal velocity (m/s)
    REAL V                     ! actual stone velocity (m/s)
    REAL TS                    ! hailstone temperature (K)
    REAL X, Y, Z               ! hail location (m)
    !RAS15.10 HAILSTONE temperature differencing
    !REAL TSm1, TSm2            ! hailstone temperature at previous 3 timesteps
    REAL FW                    ! fraction of stone that is liquid
    REAL WATER                 ! mass of stone that is liquid
    REAL CRIT                  ! critical water fraction before shedding
    REAL DENSE                 ! hailstone density (kg/m3)
    INTEGER ITYPE              ! wet (2) or dry (1) growth regime
    !in-cloud updraft parameters at location of hailstone
    REAL P                     ! in-cloud pressure (Pa)
    REAL RS                    ! in-cloud saturation mixing ratio 
    REAL RI, RW                ! ice, liquid water mix. ratio (kg/kg)
    REAL RT                    ! total hydromeoteor mix. ratio (kg/kg)
    REAL XI, XW                ! ice, liquid water content (kg/m3 air)
    REAL PC                    ! in-cloud fraction of frozen water
    REAL TC                    ! in-cloud temperature (K)
    REAL VU                    ! in-cloud updraft speed (m/s)
    REAL UC                    ! in-cloud u wind speed (m/s)
    REAL VC                    ! in-cloud v wind speed (m/s)
    REAL DENSA                 ! in-cloud updraft density (kg/m3)
    REAL DELRW                 ! diff in sat vap. dens. between hail and air (kg/m3)
                                     !used as initial hail embryo size
    REAL a_shedd               !Constant multiplied by equivalent spherical diameter 
                               ! (Deq) to get maximum diameter (Dmax) for 
                               ! non-spherical hailstones. Retrieved from Shedd
                               ! et al. 2021 Fig. 9. from 3D-scanned stones.
    !internal function variables
    REAL GM,GM1,GMW,GMI,DGM,DGMW,DGMI,DGMV,DI,ANU,RE,AE
    REAL dum, DIFF, icefactor, adiabatic_factor
    REAL ellipticity, sa_factor
      
    REAL sec, secdel           ! time step, increment in seconds
    INTEGER tcounter           ! counter for timestep
    INTEGER i, j, k, IFOUT, ind(1), inject, t
    CHARACTER*256 :: message

    !set some constants
    !secdel = 1.0
    secdel = dt
    g=9.81
    r_d = 287.
    pi = 3.14159265359
    ! CRTICAL MASS CAPABLE OF BEING "SUPPORTED" ON THE STONE'S SURFACE 
    CRIT = 2.0E-4

    !Now, calculate the difference between the surface area of a sphere
    ! and an oblate spheroid with equivalent mass. 
    ! SA_spheroid = sa_factor * SA_eq.sphere
    if ( (oblate_heat .gt. 0) .or. (oblate_massagg .gt. 0) .or. &
         (oblate_vt .gt. 0) ) then 
        ellipticity = (1 - phi_min**2.)**0.5
        sa_factor = (a_shedd**2. / 2.) * &
                    ( 1 + (phi_min**2./(2*ellipticity)) * &
                          LOG((1+ellipticity)/(1-ellipticity)) )
    else
        sa_factor = 1.0
    endif

    !Transfer input data to "interpreted" variables.  These use to be arrays,
    ! hence the "A" at the end.
    P = PA
    RS = RA
    DENSA = rho
    TC = TCA
    RI = RIA
    RW = RWA

    !Set initial hailstone parameters based on haildata
    DENSE = haildata1d(haildense)  !kg/m3
    IF (haildata1d(hailts).lt.-99) THEN   !never changed from initiation
        TS = TC
        IF (TS.GE.273.155) THEN
            ITYPE = 2  !1 dry, 2 wet
        ELSE
            ITYPE = 1  !1 dry, 2 wet
        ENDIF
    ELSE
        TS = haildata1d(hailts) !K
        ITYPE = haildata1d(hailitype)  !1 dry, 2 wet
    ENDIF
    !TSm1 = TS  !temperature differencing to account for instabilities
    !TSm2 = TS
    D = haildata1d(haild)  !m
    FW = haildata1d(hailfw)  !dimensionless
    VT = haildata1d(hailvt)  !terminal velocity calculated in hailcast_driver
    !Calculate our hailstone's mass given these variables
    GM = pi/6.*(D**3.)*DENSE
    

    !!!!!!!!!!!!!!!!!!!!!!!!! START INTERACTING !!!!!!!!!!!!!!!!!!!!!!!
    ! if (np.eq.313793) write(*,'(a30,2f12.9)') &
    !     'Entering hailtraj D, RI+RW: ', D, RI+RW
    ! if (np.eq.313793) write(*,'(a30,f5.1,1x,f12.9,1x,f9.6)') &
    !     'Entering hailtraj TC, GM, FW: ', TC, GM, FW
    

    !Has the hailstone has melted entirely? If so (and we 
    ! are also below the melting level), we can skip the rest.
    IF ((FW .EQ. 1.0) .AND. (TC .GT. 273.15)) THEN
        GOTO 300
    ENDIF

    !Are we inside the cloud?
    IF ( (RI+RW) .LT. 1.E-12) THEN
        !No water anywhere?  Let the hailstone fall, but first check to see 
        ! if it should be melting.  Right now just using a flat 273.15
        ! temperature check to avoid an iterative wet bulb calculation
        ! every time.
        IF (TC .GT. 273.15) THEN
            CALL MELT(D,TC,P,RS,secdel,VT,TS,DENSE,GM,FW,np,oblate_heat,sa_factor)
            !!!  TEST IF RETAINED WATER ON STONE IS GREATER THAN CRITICAL LIMIT, IF SO  
            !!!  BREAK UP 
            WATER=FW*GM  !KG
            IF (WATER.GT.CRIT)THEN
               CALL BREAKUP(DENSE,D,GM,FW,CRIT,np)
            ENDIF
            !calculate ice-only diameter size
            D_ICE = ( (6*GM*(1.-FW)) / (3.141592654*DENSE) )**0.33333333 
        ENDIF
        !We're done
        GOTO 400
    ENDIF
        
    !Calculate ice/liquid water content
    XI = RI * DENSA
    XW = RW * DENSA 
    IF( (XW+XI).GT.0) THEN
    PC = XI / (XW+XI)
    ELSE
    PC = 1.
    ENDIF
        
    ! SATURATION VAPOUR DENSITY DIFFERENCE BETWTEEN STONE AND CLOUD
    CALL VAPORCLOSE(DELRW,PC,TS,TC,ITYPE,np)
 
    !if (np.eq.313793) write(*,'(a30,3f9.3)') 'Before MASSAGR, FW,TS,DENSE: ', FW,TS,DENSE

    
    !!!!!!!!!!!!!!!!!!!!   STONE'S MASS GROWTH !!!!!!!!!!!!!!!!!!!!
    CALL MASSAGR(D,GM,GM1,GMW,GMI,DGM,DGMW,DGMI,DGMV,DI,ANU,RE,AE,&
             TC,TS,P,DENSE,DENSA,FW,VT,XW,XI,secdel,ITYPE,DELRW,np,&
             oblate_massagg, chi, a_shedd, sa_factor) 
             
    !if (np.eq.313793) write(*,'(a30,3f12.9)') 'After MASSAGR, D,GM,DENSE: ', D,GM,DENSE
    !if (np.eq.313793) write(*,'(a30,2f12.9)') 'After MASSAGR, FW,TS: ', FW,TS


    !!!!!!!!!!!!!!!!!!!!   HEAT BUDGET OF HAILSTONE !!!!!!!!!!!!!!!
    !CALL HEATBUD(TS,TSm1,TSm2,FW,TC,VT,DELRW,D,DENSA,GM1,GM,DGM,DGMW,  & 
    !             DGMV,DGMI,GMW,GMI,DI,ANU,RE,AE,secdel,ITYPE,P)
    CALL HEATBUD(TS,FW,TC,VT,DELRW,D,DENSA,GM1,GM,DGM,DGMW,  & 
                 DGMV,DGMI,GMW,GMI,DI,ANU,RE,AE,secdel,ITYPE,P,np, &
                 oblate_heat, chi, sa_factor)

    ! if (np.eq.313793) write(*,'(a30,3f12.9)') 'After HEATBUD, D,GM,DENSE: ', D,GM,DENSE
    ! if (np.eq.313793) write(*,'(a30,2f12.9)') 'After HEATBUD, FW,TS: ', FW,TS

 
    !!!!!!!! TEST DIAMETER OF STONE AND HEIGHT ABOVE GROUND !!!!!!!
    !!!  TEST IF DIAMETER OF STONE IS GREATER THAN CRITICAL LIMIT, IF SO  
    !!!  BREAK UP 
    WATER=FW*GM  !KG
    IF (WATER.GT.CRIT)THEN
       CALL BREAKUP(DENSE,D,GM,FW,CRIT,np)
    ENDIF
 
                
    !calculate ice-only diameter size
    D_ICE = ( (6*GM*(1.-FW)) / (3.141592654*DENSE) )**0.33333333 


    !Has the stone entirely melted and it's below the freezing level?  
    IF ((D_ICE .LT. 1.E-8) .AND. (TC.GT.273.155)) GOTO 300

    !move values to previous timestep value
    !TSm1 = TS
    !TSm2 = TSm1


!100   CONTINUE !outside model domain (second in-loop check)
!200   CONTINUE !stone reached below cloud base
300   CONTINUE !stone has entirely melted and is below freezing level
400   CONTINUE !hailstone not inside cloud


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !store interesting data in output arrays
      haildata1d(haild) = D
      haildata1d(haildense) = DENSE
      haildata1d(hailts) = TS
      haildata1d(hailfw) = FW
      haildata1d(hailvt) = VT
      haildata1d(hailitype) = ITYPE
      
      IF (haildice.ge.1) haildata1d(haildice) = D_ICE
      IF (hailqice.ge.1) haildata1d(hailqice) = RI
      IF (hailqliq.ge.1) haildata1d(hailqliq) = RW
      IF (hailtc.ge.1)   haildata1d(hailtc) = TC

  END SUBROUTINE hailstone_trajectory


  
  SUBROUTINE TERMINL(DENSA,DENSE,D,VT,TC,oblate_vt,a_shedd)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!
  !!!! INTERP: Calculate terminal velocity of the hailstone
  !!!!
  !!!! INPUT: DENSA  density of updraft air (kg/m3)
  !!!!        DENSE  density of hailstone
  !!!!        D      diameter of hailstone (m)
  !!!!        TC     updraft temperature (K)
  !!!!        oblate_vt   namelist option indicating what hailstone
  !!!!                    vt, shape relationships to use
  !!!!        a_shedd     Factor multiplied by equivalent spherical diameter 
  !!!!                    (Deq) to get maximum diameter (Dmax) for 
  !!!!                    non-spherical hailstones.
  !!!! OUTPUT:VT     hailstone terminal velocity (m/s)
  !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
      
      REAL*8 D
      REAL DENSA, DENSE, TC, VT
      REAL GMASS, GX, RE, WDUM, Y
      REAL, PARAMETER :: PI = 3.141592654, G = 9.78956
      REAL ANU, a_shedd
      REAL*8 DMAX  !DMAX = a_shedd * D (where D = Deq)
      REAL*8 DCMAX !DMAX in cm

      INTEGER oblate_vt

      if (oblate_vt .gt. 0) then
        !Use Heymsfield et al. 2020 Corregendum Vt-Dmax relationship.
        !This relationship equates roughly to the RH87 spherical 
        !relationship (see their Fig. 3), but assuming no critical
        !transition in the Cd-Re relationship (due to roughness)
        !and a lower Cd at higher Re. See Pounds et al. 2024b
        !discussion for more info.

        DMAX = D*a_shedd
        DCMAX = DMAX*1.E2  !m to cm

        IF (DCMAX .ge. 1.5) THEN
            VT = 8.4*DCMAX**0.67
        ELSE  ! if .lt. 1.5
            VT = 7.6*DCMAX**0.89
        ENDIF

    
      else  !oblate_vt == 0
        !Use RH87 Appendix B Best-Reynolds number relationship
        ! for spherical hailstones. 

        !Mass of stone in kg
        GMASS = (DENSE * PI * (D**3.)) / 6.
        
        !Dynamic viscosity
        ANU = (0.00001718)*(273.155+120.)/(TC+120.)*(TC/273.155)**(1.5)
        
        !CALC THE BEST NUMBER, X AND REYNOLDS NUMBER, RE 
        GX=(8.0*GMASS*G*DENSA)/(PI*(ANU*ANU))
        RE=(GX/0.6)**0.5

        !SELECT APPROPRIATE EQUATIONS FOR TERMINAL VELOCITY DEPENDING ON 
        !THE BEST NUMBER
        IF (GX.LT.550) THEN
            WDUM=LOG10(GX)
            Y= -1.7095 + 1.33438*WDUM - 0.11591*(WDUM**2.0)      
            RE=10**Y
            VT=ANU*RE/(D*DENSA)
        ELSE IF (GX.GE.550.AND.GX.LT.1800) THEN
            WDUM=LOG10(GX)
            Y= -1.81391 + 1.34671*WDUM - 0.12427*(WDUM**2.0) + 0.0063*(WDUM**3.0)
            RE=10**Y
            VT=ANU*RE/(D*DENSA)
        ELSE IF (GX.GE.1800.AND.GX.LT.3.45E08) THEN
            RE=0.4487*(GX**0.5536)
            VT=ANU*RE/(D*DENSA)
        ELSE 
            RE=(GX/0.6)**0.5
            VT=ANU*RE/(D*DENSA)
        ENDIF
      endif

  END SUBROUTINE TERMINL   

   
   
  SUBROUTINE VAPORCLOSE(DELRW,PC,TS,TC,ITYPE,np)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  VAPORCLOSE: CALC THE DIFFERENCE IN SATURATION VAPOUR DENSITY 
  !!!  BETWEEN THAT OVER THE HAILSTONE'S SURFACE AND THE IN-CLOUD 
  !!!  AIR, DEPENDS ON THE WATER/ICE RATIO OF THE UPDRAFT, 
  !!!  AND IF THE STONE IS IN WET OR DRY GROWTH REGIME
  !!!
  !!!  INPUT:  PC    fraction of updraft water that is frozen
  !!!          TS    temperature of hailstone (K)
  !!!          TC    temperature of updraft air (K)
  !!!          ITYPE wet (2) or dry (1) growth regime
  !!!  OUTPUT: DELRW difference in sat vap. dens. between hail and air
  !!!          (kg/m3)
  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE
      REAL DELRW, PC, TS, TC
      INTEGER ITYPE

      integer np

      !local variables
      REAL RV, ALV, ALS, RATIO
      DATA RV/461.48/,ALV/2500000./,ALS/2836050./ 
      REAL ESAT, RHOKOR, ESATW, RHOOMGW, ESATI, RHOOMGI, RHOOMG

      !!!  FOR HAILSTONE:  FIRST TEST IF STONE IS IN WET OR DRY GROWTH
      RATIO = 1./273.155
      IF(ITYPE.EQ.2) THEN !!WET GROWTH
        ESAT=611.*EXP(ALV/RV*(RATIO-1./TS)) !From (2.11) of RY89. Units Pa
      ELSE  !!DRY GROWTH
        ESAT=611.*EXP(ALS/RV*(RATIO-1./TS))
      ENDIF
      RHOKOR=ESAT/(RV*TS) !kg /m3
      
      !!!  NOW FOR THE AMBIENT/IN-CLOUD CONDITIONS 
      ESATW=611.*EXP(ALV/RV*(RATIO-1./TC))
      RHOOMGW=ESATW/(RV*TC)
      ESATI=611.*EXP(ALS/RV*(RATIO-1./TC))
      RHOOMGI=ESATI/(RV*TC)
      !RHOOMG=PC*(RHOOMGI-RHOOMGW)+RHOOMGW
      RHOOMG = RHOOMGI  !done as in hailtraj.f

      !!!  CALC THE DIFFERENCE(KG/M3): <0 FOR CONDENSATION, 
      !!!  >0 FOR EVAPORATION
      DELRW=(RHOKOR-RHOOMG) !hailstone - environment. Units kg/m3.

  END SUBROUTINE VAPORCLOSE
     
      

  SUBROUTINE MASSAGR(D,GM,GM1,GMW,GMI,DGM,DGMW,DGMI,DGMV,DI,ANU,RE,AE,& 
                 TC,TS,P,DENSE,DENSA,FW,VT,XW,XI,SEKDEL,ITYPE,DELRW,np,&
                 oblate_massagg, chi, a_shedd, sa_factor) 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! CALC THE STONE'S INCREASE IN MASS 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
      IMPLICIT NONE
      REAL*8 D
      REAL GM,GM1,GMW,GMI,DGM,DGMW,DGMI,DI,DI_SI,ANU,RE,AE,  &
                 TC,TS,P,DENSE,DENSA,FW,VT,XW,XI,SEKDEL,DELRW, &
                 chi, a_shedd, sa_factor
      INTEGER ITYPE, oblate_massagg

      integer np

      !local variables
      REAL PI, D0, GMW2, GMI2, EW, EI,DGMV
      REAL DENSEL, DENSELI, DENSELW 
      REAL DC !MEAN CLOUD DROPLET DIAMETER (MICRONS, 1E-6M)
      REAL VOLL, VOLT !VOLUME OF NEW LAYER, TOTAL (M3)
      REAL VOL1, DGMW_NOSOAK, SOAK, SOAKM
      REAL DENSAC, E, AFACTOR, NS, TSCELSIUS, VIMP, W, ANUC
      REAL*8 DMAX !max diameter of oblate spheroid (m)
      PI=3.141592654

      !!!  CALCULATE THE DIFFUSIVITY DI 
      !!D0=0.226*1.E-4  ! change to m2/s, not cm2/s
      !!DI=D0*(TC/273.155)**1.81*(100000./P)
      !Returning to RH87 Appendix A calculation for DI (their Dv).
      ! They want temperature in K, pressure in hPa, and produces
      ! results in cm2 s-1.
      DI = 0.211*(TC/273.15)**1.94 * (1013.25/(P*1E-2))  !cm2 s-1
      DI_SI = DI * (1E-2)**2.  !m2 s-1, b/c multiplied by (1m/100cm)^2


      !!!  COLLECTION EFFICIENCY FOR WATER AND ICE 
      EW=1.0      
      !!!  Linear function for ice accretion efficiency
      IF (TC .GE. 273.155) THEN
         EI=1.00
      ELSE IF (TC.GE.233.155) THEN
         EI= 1.0 - ( (273.155 - TS) / 40. )
      ELSE  !cooler than -40C
         EI = 0.0
      ENDIF

      !Something to consider - Pruppacher and Klett 1996 (2nd ed.) have an
      ! equation to parameterize EW (14-34) based on laboratory experiments
      ! from Cober and List (1993), that reduces efficiency by ~25% due to
      ! surface roughness.

      !!!  CALCULATE THE VENTILATION COEFFICIENT - NEEDED FOR GROWTH FROM VAPOR
      !The coefficients in the ventilation coefficient equations have been 
      !experimentally derived, and are expecting cal-C-g units.  However,
      !the Reynolds number is a ratio, so as long as all the units cancel
      !we should be good.
      !
      !Dynamic viscosity 
      !ANU=1.717E-4*(393.0/(TC+120.0))*(TC/273.155)**1.5
      !RAS exponent fix to agree with Roger and Yau ps. 102
      ANU=1.717E-5*(393.0/(TC+120.0))*(TC/273.155)**1.5   !kg m-1 s-1     
      !!!  CALCULATE THE REYNOLDS NUMBER - unitless, and all in SI units
      RE=D*VT*DENSA/ANU   
      E=(0.60)**(0.333333333)*(RE**0.50) !ventilation coefficient vapor (fv)
      !!!   SELECT APPROPRIATE VALUES OF AE ACCORDING TO RE
      IF(RE.LT.6000.0)THEN
         AE=0.78+0.308*E
      ELSE IF (RE.GE.6000.0.AND.RE.LT.20000.0)THEN
        if (oblate_massagg .gt. 0) then
            AE=0.5*chi*E
        else 
            AE=0.5*0.76*E !RH87 Table 1 eq. 4 shifts from factor of 4 to factor of 2
        endif
      ELSE IF (RE.GE.20000.0) THEN
        if (oblate_massagg .gt. 0) then 
            AE = 0.5*( (chi - 0.19) +9.0E-6*RE)*E !!RH87 Table 1 eq. 4 shifts from factor of 4 to factor of 2
        else
            AE = 0.5*(0.57+9.0E-6*RE)*E !!RH87 Table 1 eq. 4 shifts from factor of 4 to factor of 2
        endif
      ENDIF

      !!!  CALC HAILSTONE'S MASS (GM), MASS OF WATER (GMW) AND THE  
      !!!  MASS OF ICE IN THE STONE (GMI)
      GM=PI/6.*(D**3.)*DENSE
      GMW=FW*GM
      GMI=GM-GMW
  
      !!!  STORE THE MASS
      GM1=GM
      
      !!! NEW MASS GROWTH CALCULATIONS WITH VARIABLE RIME 
      !!! LAYER DENSITY BASED ON ZIEGLER ET AL. (1983)

      !!! If oblate sperhoid, determine new maximum diameter.
      if ( oblate_massagg .gt. 0 ) then
          DMAX = a_shedd * D !m
      else
          DMAX = D
      endif
      
      !!! CALCULATE INCREASE IN MASS DUE INTERCEPTED CLD WATER, USE
      !!! ORIGINAL DIAMETER
      GMW2=GMW+SEKDEL*(PI/4.*DMAX**2.*VT*XW*EW) !DMAX = D if oblate_massagg off
      DGMW=GMW2-GMW 
      GMW=GMW2

      !!!  CALCULATE THE INCREASE IN MASS DUE INTERCEPTED CLOUD ICE
      GMI2=GMI+SEKDEL*(PI/4.*DMAX**2.*VT*XI*EI) !DMAX = D if oblate_massagg off
      DGMI=GMI2-GMI 
      GMI=GMI2
  
      !!! CALCULATE INCREASE IN MASS DUE TO SUBLIMATION/CONDENSATION OF 
      !!! WATER VAPOR
      if ( oblate_massagg .gt. 0) then
          DGMV = -sa_factor*SEKDEL*PI*D*AE*DI_SI*DELRW
      else
          DGMV = -SEKDEL*PI*D*AE*DI_SI*DELRW
      endif
      IF (DGMV .LT. 0) DGMV=0

      !!!  CALCULATE THE TOTAL MASS CHANGE 
      DGM=DGMW+DGMI+DGMV

      !if (np.eq.313793) write(*,'(a30,i3,2(1x,f12.9))') 'Inside MASSAGR, ITYPE,DGM,GM: ', ITYPE, DGM, GM

    
      !!! CALCULATE DENSITY OF NEW LAYER, DEPENDS ON FW AND ITYPE
      IF (ITYPE.EQ.1) THEN !DRY GROWTH
          !If hailstone encountered supercooled water, calculate new layer density 
          ! using Macklin form
          IF ((DGMW.GT.0).OR.(DGMV.GT.0)) THEN
             !MEAN CLOUD DROPLET RADIUS (M), ASSUME CLOUD DROPLET CONC OF 3E8 M-3 (300 CM-3)
             DC = (0.74*XW / (3.14159*1000.*3.E8))**0.33333333 * 1.E6 !MICRONS
             !!! FIND THE STOKES NUMBER  (rasmussen heymsfield 1985)
             ! NS = 2 * VT * r_droplet * rho_droplet / (9 * dyn_viscosity * r_hail) 
             ! VT (m/s) * 100 = vt (cm/s)
             ! r_droplet = DC (micron) * 10**-4 = DC (cm)
             ! rho_droplet = 1 g/cm3
             ! D (m) * 50 = R (cm), where D (cm) = 2 * R (cm)
             ! dyn_viscosity (g cm-1 s-1), OR 
             !    ANU (kg m-1 s-1) * (1000 g/1kg) * (1m/100cm) --> 10 * ANU g cm-1 s-1
             NS = 2*VT*100.*(DC*1.E-4)**2. / (90*ANU*D*50)  !need hail radius in cm
             !!! FIND IMPACT VELOCITY (rasmussen heymsfield 1985)
             IF (NS.LT.0.1)THEN
                W=-1.
             ELSE
                W = LOG10(NS)
             ENDIF
             IF (RE.GT.200) THEN
                IF (NS.LE.0.1) THEN
                   VIMP = 0.0
                ELSEIF ((NS.GT.0.1).AND.(NS.LE.10)) THEN
                   VIMP = (0.356 + 0.4738*W - 0.1233*W**2. &
                           -0.1618*W**3. + 0.0807*W**4.)*VT
                ELSEIF (NS.GT.10) THEN
                   VIMP = 0.63*VT
                ENDIF
             ELSEIF ((RE.GT.65).AND.(RE.LE.200)) THEN
                IF (NS.LE.0.1) THEN
                   VIMP = 0.0
                ELSEIF ((NS.GT.0.1).AND.(NS.LE.10)) THEN
                   VIMP = (0.3272 + 0.4907*W - 0.09452*W**2. &
                           -0.1906*W**3. + 0.07105*W**4.)*VT
                ELSEIF (NS.GT.10) THEN
                   VIMP = 0.61*VT
                ENDIF
             ELSEIF ((RE.GT.20).AND.(RE.LE.65)) THEN
                IF (NS.LE.0.1) THEN
                   VIMP = 0.0
                ELSEIF ((NS.GT.0.1).AND.(NS.LE.10)) THEN
                   VIMP = (0.2927 + 0.5085*W - 0.03453*W**2. &
                           -0.2184*W**3. + 0.03595*W**4.)*VT
                ELSEIF (NS.GT.10) THEN
                   VIMP = 0.59*VT
                ENDIF
             ELSEIF (RE.LE.20) THEN
                IF (NS.LE.0.4) THEN
                   VIMP = 0.0
                ELSEIF ((NS.GT.0.4).AND.(NS.LE.10)) THEN
                   VIMP = (0.1701 + 0.7246*W + 0.2257*W**2. &
                           -1.13*W**3. + 0.5756*W**4.)*VT
                ELSEIF (NS.GT.10) THEN
                   VIMP = 0.57*VT
                ENDIF
             ENDIF
 
              
             !RIME LAYER DENSITY, HEYMSFIELD AND PFLAUM 1985 FORM
             TSCELSIUS = TS - 273.16 ! The constant here, 273.16, is purposefully chosen
                                     ! to be slightly larger than 273.155. That way a TS of
                                     ! exactly 273.155 doesn't result in an undefined
                                     ! AFACTOR. When in transition, for a few timesteps 
                                     ! hailstones can have a TS of 273.155 but no retained 
                                     ! liquid water on their surface. Is such a hailstone 
                                     ! in officially wet or dry growth? Not clear. This
                                     ! kludge is the result.
             AFACTOR = -DC*VIMP/TSCELSIUS
             !IF ((TSCELSIUS.LE.-5.).OR.(AFACTOR.GE.-1.60)) THEN
             !error in AFACTOR sign in HF85
             IF ((TSCELSIUS.LE.-5.).AND.(AFACTOR.GE.1.60)) THEN  
                 DENSELW = 0.30*(AFACTOR)**0.44
             !ELSEIF (TSCELSIUS.GT.-5.) THEN
             ELSE
                 DENSELW = EXP(-0.03115 - 1.7030*AFACTOR + &
                               0.9116*AFACTOR**2. - 0.1224*AFACTOR**3.)
             ENDIF

             DENSELW = DENSELW * 1000. !KG M-3

             !BOUND POSSIBLE DENSITIES
             IF (DENSELW.LT.500) DENSELW=500  !match with Kumjian and Lombardo 2020
             IF (DENSELW.GT.900) DENSELW=900
             !WRITE(12,*) 'MASSAGR, PFLAUM, MACKLIN: ', DENSELW, &
             !             MACKLIN_DENSELW
          ENDIF
          IF (DGMI.GT.0) THEN
             !Ice collection main source of growth, so set new density layer
             DENSELI = 700.
          ENDIF
          
          !All liquid water contributes to growth, none is soaked into center.
          DGMW_NOSOAK = DGMW  !All liquid water contributes to growth,
                              ! none of it is soaked into center.

      ELSE !WET GROWTH
          !Collected liquid water can soak into the stone before freezing,
          ! increasing mass and density but leaving volume constant.
          !Volume of current drop, before growth 
          VOL1 = GM/DENSE
          !Difference b/w mass of stone if density is 900 kg/m3, and
          ! current mass
          SOAK = 900*VOL1 - GM
          !Liquid mass available
          SOAKM = DGMW
          !Soak up as much liquid as we can, up to a density of 900 kg/m3
          IF (SOAKM.GT.SOAK) SOAKM=SOAK
          GM = GM+SOAKM  !Mass of current drop, plus soaking
          !New density of current drop, including soaking but before growth
          DENSE = GM/VOL1 
          !Mass increment of liquid water growth that doesn't
          ! include the liquid water we just soaked into the stone.
          DGMW_NOSOAK = DGMW - SOAKM
          
          !Whatever growth does occur has high density
          DENSELW = 900.  !KG M-3
          DENSELI = 900.
         
      ENDIF

      !!!VOLUME OF NEW LAYER
      IF (DGM.GT.0) THEN
          IF (DGMI.LE.0) THEN
             VOLL = (DGMW_NOSOAK+DGMV) / DENSELW
          ELSE IF (DGMW.LE.0) THEN
             VOLL = (DGMI) / DENSELI
          ELSE
             VOLL = (DGMI) / DENSELI + (DGMW_NOSOAK+DGMV) / DENSELW
          ENDIF
      ELSE
          VOLL = 0.
      ENDIF

      !!!NEW TOTAL VOLUME, DENSITY, DIAMETER
      VOLT = VOLL + GM/DENSE
      DENSE = (GM+DGMI+DGMV+DGMW_NOSOAK) / VOLT
      GM = GM+DGMI+DGMW_NOSOAK+DGMV
      D = ( (6*GM) / (PI*DENSE) )**0.33333333 


  END SUBROUTINE MASSAGR



  SUBROUTINE HEATBUD(TS,FW,TC,VT,DELRW,D,DENSA,GM1,GM,DGM,DGMW,     &
                     DGMV,DGMI,GMW,GMI,DI,ANU,RE,AE,SEKDEL,ITYPE,P, np, &
                     oblate_heat, chi, sa_factor)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! CALCULATE HAILSTONE'S HEAT BUDGET 
  !!! See Rasmussen and Heymsfield 1987; JAS
  !!! Original Hailcast's variable units
  !!! TS - Celsius
  !!! FW - unitless, between 0 and 1
  !!! TC - Celsius
  !!! VT - m/s
  !!! D  - m
  !!! DELRW - g/cm3 (per comment)
  !!! DENSA - g/cm3 (per comment)
  !!! GM1, DMG, DGMW, DGMV, DGMI, GMW, GMI - should all be kg
  !!! DI - cm2 / sec
  !!! P  - hPa
  !!! Original HAILCAST HEATBUD subroutine uses c-g-s units, so do some conversions
  !!! 
  !!! oblate_heat: namelist option deteriming if should use a spherical (0)
  !!!              or oblate (1) hailstone. Oblate hailstone assumed to 
  !!!              fall with its maximum diameter horizontal.
  !!! chi: heat transfer rate, as in Macklin 1963.
  !!! sa_factor: scaling factor between surface area of oblate spheroid and
  !!!            surface area of equivalent sphere of same mass.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
      REAL*8 D
      REAL TS,TSm1,TSm2,FW,TC,VT,DELRW,DENSA,GM1,GM,DGM,DGMW,DGMV,  &
                    DGMI,GMW,GMI,DI,ANU,RE,AE,SEKDEL,P, chi, &
                    sa_factor
      INTEGER ITYPE, oblate_heat

      integer np
      
      REAL RV, RD, G, PI, ALF, ALV, ALS, CI, CW, AK
      REAL H, AH, TCC, TSC, DELRWC, DENSAC, TDIFF
      REAL DMLT
      REAL TSCm1, TSCm2
      REAL DCM, GMC, GM1C, DGMWC, DGMIC, DGMC
      DATA RV/461.48/,RD/287.04/,G/9.78956/
      DATA PI/3.141592654/,ALF/79.7/,ALV/597.3/
      DATA ALS/677.0/,CI/0.5/,CW/1./
      
      !Convert values to non-SI units here
      TSC = TS - 273.155
      !TSCm1 = TSm1 - 273.155
      !TSCm2 = TSm2 - 273.155
      TCC = TC - 273.155
      DELRWC = DELRW * (1.E3) * (1.E-6)  !kg/m3 --> g/cm3
      DENSAC = DENSA * (1.E3) * (1.E-6)
      !DI still in cm2/sec
      !RAS 20220318
      DCM = D * 100.  !m to cm
      GMC = GM * 1000.  !kg to g
      GM1C = GM1 * 1000.  !kg to g
      DGMWC = DGMW * 1000.  !kg to g
      DGMIC = DGMI * 1000.  !kg to g
      DGMC = DGM * 1000. !kg to g



      !!!  CALCULATE THE CONSTANTS 
      AK=(5.8+0.0184*TCC)*1.E-5  !thermal conductivity - cal/(cm*sec*K)
      !dynamic viscosity - calculated in MASSAGR
      !ANU=1.717E-5*(393.0/(TC+120.0))*(TC/273.155)**1.5 kg m-1 s-1

      !!!  CALCULATE THE REYNOLDS NUMBER - unitless. 
      !RE=D*VT*DENSA/ANU  - calculated in MASSAGR, where all inputs are in SI units.
      
      H=(0.71)**(0.333333333)*(RE**0.50) !ventilation coefficient heat (fh)
      !E=(0.60)**(0.333333333)*(RE**0.50) !ventilation coefficient vapor (fv)

      !!!   SELECT APPROPRIATE VALUES OF AH AND AE ACCORDING TO RE
      IF(RE.LT.6000.0)THEN
         AH=0.78+0.308*H !RH87 Table 1 eq. 3
      ELSE IF(RE.GE.6000.0.AND.RE.LT.20000.0)THEN
         if ( oblate_heat .gt. 0 ) then
            AH=0.5 * chi*H
         else            
            AH=0.5 * 0.76*H  !RH87 Table 1 eq. 4 shifts from factor of 4 to factor of 2
         end if
      ELSE IF(RE.GE.20000.0) THEN
        if ( oblate_heat .gt. 0 ) then
           AH = 0.5 * ( (chi - 0.19) + 9.0E-6*RE)*H
        else
           AH = 0.5 * (0.57 + 9.0E-6*RE)*H !RH87 Table 1 eq. 5 shifts from factor of 4 to factor of 2
        endif
      ENDIF

      !!!  FOR DRY GROWTH FW=0, CALCULATE NEW TS, ITIPE=1 
      !!!  FOR WET GROWTH TS=0, CALCULATE NEW FW, ITIPE=2

      IF(ITYPE.EQ.1) THEN
      !!!  DRY GROWTH; CALC NEW TEMP OF THE STONE 
         !Original Hailcast algorithm (no time differencing)
         !TSC=TSC-TSC*DGM/GM1+SEKDEL/(GM1*CI)*                &
         !   (2.*PI*D*(AH*AK*(TCC-TSC)-AE*ALS*DI*DELRWC)+     &
         !   DGMW/SEKDEL*(ALF+CW*TCC)+DGMI/SEKDEL*CI*TCC)
         !Time differencing option
         !TSC=0.6*(TSC-TSC*DGM/GM1+SEKDEL/(GM1*CI)*                &
         !   (2.*PI*D*(AH*AK*(TCC-TSC)-AE*ALS*DI*DELRWC)+     &
         !   DGMW/SEKDEL*(ALF+CW*TCC)+DGMI/SEKDEL*CI*TCC)) + &
         !   0.2*TSCm1 + 0.2*TSCm2
        !Note about signs, RAS - !The second and third lines
        !  in the coded equation below are equivalent to
        ! -dq/dt in RH87 Table 1 Eqs. 3-5. dq/dt is positive
        ! when heat is fluxing AWAY from the hailstone, which 
        ! would result in a decrease in TSC. Using a -dq/dt
        ! term ensures it is appropriately related to TSC.
        !The TCC term below is T_env - T_hailstone. 
        ! DELRWC, conversely, is rho_hailstone - rho_env. Because
        ! these are different, a different sign is used before the 
        ! term with DELRWC.
         if ( oblate_heat .gt. 0) then
            TSC=TSC-TSC*DGMC/GM1C+SEKDEL/(GM1C*CI)*                &
                (sa_factor*2.*PI*DCM*(AH*AK*(TCC-TSC)-AE*ALS*DI*DELRWC)+     &
                DGMWC/SEKDEL*(ALF+CW*TCC)+DGMIC/SEKDEL*CI*TCC)
         else
            TSC=TSC-TSC*DGMC/GM1C+SEKDEL/(GM1C*CI)*                &
            (2.*PI*DCM*(AH*AK*(TCC-TSC)-AE*ALS*DI*DELRWC)+     &
            DGMWC/SEKDEL*(ALF+CW*TCC)+DGMIC/SEKDEL*CI*TCC)           
         endif
         TS = TSC+273.155
         IF (TS.GE.273.155) THEN 
            TS=273.155
         ENDIF
         !TDIFF = ABS(TS-273.155)         
         !IF (TDIFF.LE.1.E-6) ITYPE=2  !NOW IN WET GROWTH
         IF ((TS.GE.273.155).OR.(TC.GT.273.155)) ITYPE=2  !NOW IN WET GROWTH
     
      ELSE IF (ITYPE.EQ.2) THEN
      !!!  WET GROWTH; CALC NEW FW          
         
         IF (TCC.LT.0.) THEN
            !Original Hailcast algorithm
            !FW=FW-FW*DGM/GM1+SEKDEL/(GM1*ALF)*               &
            !    (2.*PI*D*(AH*AK*TCC-AE*ALV*DI*DELRWC)+          &
            !    DGMW/SEKDEL*(ALF+CW*TCC)+DGMI/SEKDEL*CI*TCC)
            !Note about signs, RAS - !The second and third lines
            !  in the coded equation below are equivalent to
            ! -dq/dt in RH87 Table 1 Eqs. 3-5. dq/dt is positive
            ! when heat is fluxing AWAY from the hailstone, which 
            ! would result in a decrease in FW. Using a -dq/dt
            ! term ensures it is appropriately related to FW.
            !The TCC term below is T_env - T_hailstone. 
            ! DELRWC, conversely, is rho_hailstone - rho_env. Because
            ! these are different, a different sign is used before the 
            ! term with DELRWC.
            if (oblate_heat .gt. 0) then
                FW = FW-FW*DGMC/GMC + &
                    SEKDEL/(GMC*ALF)*&
                    (sa_factor*2.*PI*DCM*(AH*AK*TCC-AE*ALV*DI*DELRWC)+ &
                    DGMWC/SEKDEL*(ALF+CW*TCC)+DGMIC/SEKDEL*CI*TCC)
            else
                FW = FW-FW*DGMC/GMC + &
                    SEKDEL/(GMC*ALF)*&
                    (2.*PI*DCM*(AH*AK*TCC-AE*ALV*DI*DELRWC)+ &
                    DGMWC/SEKDEL*(ALF+CW*TCC)+DGMIC/SEKDEL*CI*TCC)
            endif                
         ELSE
            !Calculate decrease in ice mass due to melting
            !Note about signs, RAS - when TCC is positive it means an
            ! increase in melting. When DELRWC is positive, it means
            ! a flux of vapor away from the hailstone, but that's a decrease
            ! in temperature. Hence, second sign remains negative.
            if (oblate_heat .gt. 0) then
                DMLT = SEKDEL * (sa_factor*2.*PI*DCM*AH*AK*TCC - &
                                 sa_factor*2.*PI*DCM*AE*ALV*DI*DELRWC + &
                                 DGMWC/SEKDEL*CW*TCC) / ALF
            else
                DMLT = SEKDEL * (2.*PI*DCM*AH*AK*TCC - &
                                 2.*PI*DCM*AE*ALV*DI*DELRWC + &
                                 DGMWC/SEKDEL*CW*TCC) / ALF
            endif                
            FW = (FW*GM1C + DMLT + DGMWC) / GMC
         ENDIF
         
         IF(FW.GT.1.)FW=1.
         IF(FW.LT.0.)FW=0.

         !IF ALL OUR ACCRETED WATER WAS FROZEN, WE ARE BACK IN DRY GROWTH
         IF(FW.LE.1.E-6) THEN
            ITYPE=1  
         ENDIF
         
      ENDIF

  END SUBROUTINE HEATBUD


  
  SUBROUTINE BREAKUP(DENSE,D,GM,FW,CRIT,np)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  TEST IF AMOUNT OF WATER ON SURFACE EXCEEDS CRTICAL LIMIT- 
  !!!  IF SO INVOKE SHEDDING SCHEME 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE
      REAL*8 D
      REAL DENSE, GM, FW, CRIT
      integer np

      !local variables
      REAL WATER, GMI, WAT, PI
      DATA PI/3.141592654/

      WATER=FW*GM  !KG
      !GMI=(GM-WATER) !KG
      !REMAIN = CRIT*GM

      ! CALC CRTICAL MASS CAPABLE OF BEING "SUPPORTED" ON THE STONE'S 
      ! SURFACE 
      !CRIT=0.268+0.1389*GMI 
      !CRIT=0.268*1.E-3 + 0.1389*1.E-3*GMI  !mass now in kg instead of g
      !CRIT = 1.0E-10
      !CRIT - now passed from main subroutine

      WAT=WATER-CRIT
      GM=GM-WAT
      FW=(CRIT)/GM
    
      IF(FW.GT.1.0) FW=1.0
      IF(FW.LT.0.0) FW=0.0


      ! RECALCULATE DIAMETER AFTER SHEDDING 
      ! Assume density remains the same
      D=(6.*GM/(PI*DENSE))**(0.333333333)


  END SUBROUTINE BREAKUP
  
  
  
  SUBROUTINE MELT(D,TLAYER,PLAYER,RLAYER,SEC,VT,TS,DENSE,GM,FW,np,&
                  oblate_heat, sa_factor)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!  This is a spherical hail melting estimate based on the Goyer 
  !!!  et al. (1969) eqn (3).  The depth of the warm layer, estimated 
  !!!  terminal velocity, and mean temperature of the warm layer are 
  !!!  used.  DRB.  11/17/2003.
  !!!
  !!!  Incorporated some known variables into the subroutine
  !!!  to use instead of constants (VT, TS, DENSE).  RAS 10/2/2019
  !!!
  !!!  Increases water fraction, then determines in main subroutine
  !!!   what should be shed.
  !!!
  !!!  Oblate hailstone note, RAS 17 Jan 2024: Goyer found a chi
  !!!  of 1 to better match with his observations of spherical hail, 
  !!!  rather than the chi of 0.76 found for spheric hail by Macklin.
  !!!  That was different enough I didn't want to muck with changing
  !!!  around the Goyer assumptions about chi for non-spherical hail,
  !!!  as I really don't have any idea what they could be.
  !!!  So, the only change made was to add the change in surface area.
  !!! 
  !!!  INPUT:  TLAYER   mean sub-cloud layer temperature (K)
  !!!          PLAYER   mean sub-cloud layer pressure (Pa)
  !!!          RLAYER   mean sub-cloud layer mixing ratio (kg/kg)
  !!!          VT       terminal velocity of stone (m/s)
  !!!          TS       temperature of the hailstone (K)
  !!!          DENSE    density of the hailstone (kg/m3)
  !!!          D        diameter (m)
  !!!          GM       mass of the hailstone (kg)
  !!!          SEC      time interval over which hailstone is falling (s)
  !!!        sa_factor: scaling factor between surface area of oblate spheroid and
  !!!                   surface area of equivalent sphere of same mass.
  !!!      oblate_heat: namelist option, set > 0 if want to use oblate hailstones
  !!!  OUTPUT: FW       water fraction (b/w 0 and 1)
  !!!          
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE

      REAL*8 D
      REAL TLAYER, PLAYER, RLAYER, SEC, VT, TS, DENSE, GM, FW
      REAL eenv, delta, ewet, de, der, wetold, wetbulb, wetbulbk
      REAL tdclayer, tclayer, eps, b, hplayer
      REAL*8 a
      REAL sd, lt, ka, lf, lv, t0, dv, pi, rv, rhoice, &
           tres, re, delt, esenv, rhosenv, essfc, rhosfc, dsig, &
           dmdt, mass, massorg, newmass, gamma, r, rho, &
           sa_factor
      INTEGER wcnt, oblate_heat

      integer np
      
      !Convert temp to Celsius, calculate dewpoint in celsius
      tclayer = TLAYER - 273.155
      a = 2.53E11
      b = 5.42E3
      eps = 0.622
      tdclayer = b / LOG(a*eps / (rlayer*player))
      hplayer = player / 100.
      
      !Calculate partial vapor pressure
      eenv = (player*rlayer) / (rlayer+eps)
      eenv = eenv / 100.  !convert to mb
      
      !Estimate wet bulb temperature (C)
      gamma = 6.6E-4*player
      delta = (4098.0*eenv)/((tdclayer+237.7)*(tdclayer+237.7))
      wetbulb = ((gamma*tclayer)+(delta*tdclayer))/(gamma+delta)
      
      !Iterate to get exact wet bulb
      wcnt = 0
      DO WHILE (wcnt .lt. 11)
        ewet = 6.108*(exp((17.27*wetbulb)/(237.3 + wetbulb))) 
        de = (0.0006355*hplayer*(tclayer-wetbulb))-(ewet-eenv)
        der= (ewet*(.0091379024 - (6106.396/(273.155+wetbulb)**2))) &
             - (0.0006355*hplayer)
        wetold = wetbulb
        wetbulb = wetbulb - de/der
        wcnt = wcnt + 1
        IF ((abs(wetbulb-wetold)/wetbulb.gt.0.0001)) THEN
           EXIT
        ENDIF
      ENDDO
      
      wetbulbk = wetbulb + 273.155  !convert to K
      ka = .02 ! thermal conductivity of air
      lf = 3.34e5 ! latent heat of melting/fusion
      lv = 2.5e6  ! latent heat of vaporization
      t0 = 273.155 ! temp of ice/water melting interface
      dv = 0.25e-4 ! diffusivity of water vapor (m2/s)
      pi = 3.1415927
      rv = 1004. - 287. ! gas constant for water vapor
      !rhoice = 917.0 ! density of ice (kg/m**3)
      rhoice = DENSE  ! we know the stone density, let's use it
      r = D/2. ! radius of stone (m)
      
      !Compute residence time in warm layer
      !tres = LDEPTH / VT
      tres = SEC  !we know the residence time directly, use it 
        
      !Calculate dmdt based on eqn (3) of Goyer et al. (1969)
      !Reynolds number...from pg 317 of Atmo Physics (Salby 1996)
      !Just use the density of air at 850 mb...close enough.
      !rho = 85000./(287.*TLAYER)
      rho = player/(287.*TLAYER)  !we have the layer pressure, just use it
      !re = rho*r*VT*.01/1.7e-5
      !units fix - r is now in meters, not mm. Plus we need D, not r. RAS 20220308
      re = rho*r*2*VT/1.7e-5
      
      !Temperature difference between environment and hailstone surface
      !delt = wetbulb !- 0.0 !assume stone surface is at 0C
      !                      !wetbulb is in Celsius
      !We know the stone surface temperature, let's use it
      delt = wetbulb - (TS - 273.155) !again, wet bulb is in C


      !Difference in vapor density of air stream and equil vapor
      !density at the sfc of the hailstone
      esenv = 610.8*(exp((17.27*wetbulb)/  &
               (237.3 + wetbulb))) ! es environment in Pa
      rhosenv = esenv/(rv*wetbulbk)
      essfc = 610.8*(exp((17.27*(t0-273.155))/  &
               (237.3 + (t0-273.155)))) ! es environment in Pa
      rhosfc = essfc/(rv*t0)
      dsig = rhosenv - rhosfc  !kg /m3

      !Calculate new mass growth
      if (oblate_heat .gt. 0 ) then
         dmdt = (-sa_factor*1.7*pi*r*(re**0.5)/lf)*((ka*delt)+((lv-lf)*dv*dsig))
      else
        dmdt = (- 1.7*pi*r*(re**0.5)/lf)*((ka*delt)+((lv-lf)*dv*dsig))
      endif
      IF (dmdt.gt.0.) dmdt = 0
      mass = - dmdt*tres  !flip from <0; now represents frozen mass lost due to melt

      !!Find the new hailstone diameter
      !!massorg = 1.33333333*pi*r*r*r*rhoice   !reminder now rhoice = DENSE
      !!newmass = massorg + mass
      !!if (newmass.lt.0.0) newmass = 0.0
      !!D = 2.*(0.75*newmass/(pi*rhoice))**0.333333333
      
      !instead of finding a new hailstone diameter, increase the fraction
      ! of water being retained on the hailstone surface.
      FW = (FW*GM + mass) / GM

      IF(FW.GT.1.0) FW=1.0
      IF(FW.LT.0.0) FW=0.0

  END SUBROUTINE MELT






  END MODULE hailcast_module
