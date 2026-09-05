!--------------
  module march
!--------------
  use eig
  use scalar3
  use legops
  use layout
  
  implicit none
  integer, parameter :: p8=selected_real_kind(p=11)
  real(p8),parameter:: pi=3.141592653589793238462643_p8
!-------------------------------------------------------------
  private
  public:: scalar
  public:: timedata,      tim
  public:: imposed_field, adv
  public:: removal,       rmv
  public:: viscosity,     visc
  public:: iofiles,       files
  public:: diagdata,      diagv
  public:: lineardata,    lin

  public:: adamsb, rich, euler, rich2
  public:: diagnost
  public:: visc1,visc2,hyperv

  public:: remove
  public:: freeadj

  public:: prodctm, prodctk
     ! compute integration of two function for each azimuthal/axial 
     ! wave number

  public:: enemon
     ! compute energy for each azimuthal/axial wavenumber

  public:: hypadj

  public:: postProcessKit, postProcess
  public:: monitorEigenValuesKit, monitorData 
!-------------------------------------------------------------
  type timedata
    real(p8):: dt,t,limit
    integer:: n
    integer:: scheme
 endtype timedata
  type(timedata):: tim

  type imposed_field
    integer:: sw, int  
    real(p8):: ux,uy,uz
    real(p8):: x,y
    real(p8):: t0,x0,y0,t1,x1,y1
 endtype imposed_field
  type(imposed_field):: adv

  type removal
    integer:: sw,int
 endtype removal
  type(removal):: rmv

  type viscosity
    integer:: sw,p,adjsw,adjint
    real(p8):: nu,nup
 endtype viscosity
  type(viscosity):: visc

  type lineardata
    integer:: mm,kk
 endtype lineardata
  type(lineardata):: lin

  type iofiles
    character(len=72):: psi0,chi0
    character(len=72):: psii,chii
    real(p8),dimension(:),pointer:: t
    character(len=72),dimension(:),pointer:: psi,chi
    integer:: n,ne
 endtype iofiles
  type(iofiles):: files

  type diagdata
    integer:: it1,it2,iniflag
    real(p8):: vori,zmomi,enei,angi,heli,xcen,ycen
 endtype diagdata
  type(diagdata):: diagv

  type postProcessKit
     integer :: start, finish, job, sliceInt
  endtype postProcessKit
  type(postProcessKit):: postProcess

  type monitorEigenValuesKit
     integer :: start, interval
  endtype monitorEigenValuesKit
  type(monitorEigenValuesKit):: monitorData


!-------------------------------------------------------------
  contains
! --------

   subroutine adamsb(psi,chi,psino,chino)
!  update psi,chi by one time step by adams-bashforth method
!  psi: toroidal component
!  chi: poloidal component
!  psino,chino: nonlinear components of one time step before
!-----------------------------------------------------------
   type(scalar):: psi,chi,psino,chino
   type(scalar):: psin,chin
   type(scalar):: psi0,chi0
   real(p8):: dt

   dt = tim%dt
   tim%t = tim%t + dt
   tim%n = tim%n + 1
   adv%x = adv%x + adv%ux*dt
   adv%y = adv%y + adv%uy*dt

   call allocate( psin )
   call allocate( chin )
   call nonlin(psi,chi,psin,chin)

   call allocate( psi0 )
   call allocate( chi0 )
   psi0 = psi
   chi0 = chi

   psi%e = psi%e +dt*(1.5_p8*psin%e -0.5_p8*psino%e)
   chi%e = chi%e +dt*(1.5_p8*chin%e -0.5_p8*chino%e)
   psino=psin
   chino=chin
   !there you go, cool eh?

   call deallocate( psin )
   call deallocate( chin )

   call visc2(psi,chi,psi0,chi0)
   call hyperv(psi,chi,dt)

   call deallocate( psi0 )
   call deallocate( chi0 )

   return
   end subroutine


   subroutine rich(psi,chi,psin,chin)
! ----------------------------------------------------
!  update psi,chi by one time step. 
!  by richardson extrapolation
!  > input
!  psi: toroidal component
!  chi: poloidal component
!  ischem: designate what psio,chio are. see below
!  > output
!  psi,chi: updated fields.
!  psio,chio: 
!    ischem=2: psin,chin are copied to psio and chio.(adamsb)
!  all the arrays in f/f space.
! ----------------------------------------------------
   type(scalar):: psi,chi,psin,chin
   type(scalar):: psi2,chi2
   type(scalar):: pn,cn
   real(p8):: dt, hdt

   dt = tim%dt
   hdt = dt*0.5_p8
   tim%t = tim%t + dt
   tim%n = tim%n + 1
   adv%x = adv%x + adv%ux*dt
   adv%y = adv%y + adv%uy*dt

   call allocate(psi2)
   call allocate(chi2)
   call allocate( pn )
   call allocate( cn )

   ! first half-step
   call nonlin(psi,chi,psin,chin)
   psi2%e =psi%e + (hdt*psin%e)
   chi2%e =chi%e + (hdt*chin%e)
   call visc1(psi2,chi2,hdt)
   call hyperv(psi2,chi2,hdt)

   ! second half-step
   call nonlin(psi2,chi2,pn,cn)
   psi2%e =psi2%e + (hdt*pn%e)
   chi2%e =chi2%e + (hdt*cn%e)
   call visc1(psi2,chi2,hdt)
   call hyperv(psi2,chi2,hdt)

   ! full step
   psi%e =psi%e +(dt)*psin%e
   chi%e =chi%e +(dt)*chin%e
   call visc1(psi,chi,dt)
   call hyperv(psi,chi,dt)

   psi%e =(2.0_p8*psi2%e)-psi%e
   chi%e =(2.0_p8*chi2%e)-chi%e

   call deallocate( psi2 )
   call deallocate( chi2 )
   call deallocate( pn )
   call deallocate( cn )

   return
 end subroutine rich

   !--------------------------------------------------------
   !Added by Sujit for testing. Nov 2002
   subroutine rich2(psi,chi,psin,chin)
     ! ----------------------------------------------------
     !  update psi,chi by one time step. 
     !  by forward Eulers 
     !  > input
     !  psi: toroidal component
     !  chi: poloidal component
     !  ischem: designate what psio,chio are. see below
     !  > output
     !  psi,chi: updated fields.
     !  psio,chio: 
     !    ischem=2: psin,chin are copied to psio and chio.(adamsb)
     !  all the arrays in f/f space.
     ! ----------------------------------------------------
     type(scalar):: psi,chi,psin,chin
     type(scalar):: psi2,chi2
     type(scalar):: pn,cn
     real(p8):: dt, hdt
     complex(p8),dimension(:),allocatable:: temp
     complex(p8),parameter:: ii=(0.0_p8,1.0_p8)

     allocate(temp(16))
     
     dt = tim%dt
     hdt = dt*0.5_p8
     tim%t = tim%t + dt
     tim%n = tim%n + 1
     adv%x = adv%x + adv%ux*dt
     adv%y = adv%y + adv%uy*dt
     
     call allocate(psi2)
     call allocate(chi2)
     call allocate( pn )
     call allocate( cn )
     
     call nonlin(psi,chi,psin,chin)

     temp(1:16) = ii*chin%e(1:16,2,81)/chi%e(1:16,2,81)
     print*, temp(:)
     print*,''
!     call msave(temp, 'm1')

     temp(1:16) = ii*chin%e(1:16,1,2)/chi%e(1:16,1,2)	

! Special:
!temp(1) = 0.0_p8	

     print*, temp(:)
     print*,''
!     call msave(temp, 'm2')

     call msave(psin%e(1:16,2,3),'non1')
     
stop
     temp(1:16) = ii*chin%e(1:16,2,3)/chi%e(1:16,2,3)
     print*, temp(:)
     print*,''
     call msave(temp, 'mBar')

     temp(1:16) = ii*chin%e(1:16,2,82)/chi%e(1:16,2,82)
     print*, temp(:)
     print*,''
     call msave(temp, 'I3mBar')

     call deallocate( pn )
     call deallocate( cn )
     call deallocate(psi2)
     call deallocate(chi2)

     stop 

 end subroutine rich2
!-----------------------------------------------------


   subroutine euler(psi,chi,psin,chin)
   type(scalar):: psi,chi,psin,chin
   real(p8):: dt

   dt = tim%dt
   tim%t = tim%t + dt
   tim%n = tim%n + 1
   adv%x = adv%x + adv%ux*dt
   adv%y = adv%y + adv%uy*dt

   call nonlin(psi,chi,psin,chin)
   psi%e =psi%e + (dt*psin%e)
   chi%e =chi%e + (dt*chin%e)
   call visc1(psi,chi,dt)
   call hyperv(psi,chi,dt)

   return
   end subroutine


   subroutine visc2(psi,chi,psio,chio)
! ----------------------------------------------------
!  viscous step by crank-nicolson method
!
!    psi(1)-psi(1/2)               psi(1)+psi(0)
!    ------------- = viscnu(del2)(--------------)
!          dt                            2
!
!  call in f/f space
!  in:  psi(0), chi(0) comes as psio,chio
!       psi(1/2), chi(1/2) comes as psi,chi
!  out: psi(1), chi(1) returns in psi,chi
!       psio, chio returns intact
!  p167
! ----------------------------------------------------
   type(scalar):: psi,chi,psio,chio,wk1
   real(p8):: alp

   if(visc%sw.ne.1) return
   if(visc%nu.eq.0.0) then
     print *,'visc: visc%nu can''t be zero'
     stop
   endif

   call allocate( wk1 )

   alp=2/visc%nu/tim%dt

   call chopset(2)
   call del2(psio,wk1)
   call chopset(-2)
   call chopdo(wk1)
   wk1%e = -wk1%e-alp*psi%e
   wk1%ln= -alp*psi%ln
   call ihelm(wk1,psi,alp)

   call chopset(2)
   call del2(chio,wk1)
   call chopset(-2)
   call chopdo(wk1)
   wk1%e = -wk1%e-alp*chi%e
   wk1%ln= -alp*chi%ln
   call ihelm(wk1,chi,alp)

   call deallocate( wk1 )

   return
   end subroutine


   subroutine visc1(psi,chi,dt)
!     viscous step by implicit euler
!
!       psi(1)-psi(0)
!       ------------- = viscnu(del2)(psi(1))
!             dt
!
!     call in f/f space
!     in:  psi(0), chi(0) comes as psi,chi
!          delt: time step passed explicitly.
!     out: psi(1), chi(1) returns in psi,chi
! ----------------------------------------------------
   type(scalar):: psi,chi,wk1
   real(p8):: dt,alp
      
   if(visc%sw.ne.1) return
   if(visc%nu.eq.0.0) then
     print *,'visc: visc%nu can''t be zero'
     stop
   endif

   call allocate( wk1 )
   alp=1/(visc%nu*dt)

   wk1%ln=-alp*psi%ln
   wk1%e =-alp*psi%e
   call ihelm(wk1,psi,alp)

   wk1%ln=-alp*chi%ln
   wk1%e =-alp*chi%e
   call ihelm(wk1,chi,alp)

   call deallocate( wk1 )

   return
   end subroutine


   subroutine hyperv(psi,chi,dt)
!  do hyperviscosity dumping of the high coefficients
! ----------------------------------------------------
   type(scalar):: psi,chi,w
   real(p8):: dt,alp,bet
   integer:: ph

   if(visc%sw.eq.0) then
     print *,'hyperv: warning: no dissipation at all.'
   else if(visc%sw.eq.2) then
     call allocate( w )

     ph = visc%p/2
     alp = 1/dt/visc%nup *(-1)**ph
     bet = -(-1)**ph*visc%nu/visc%nup

     w%e = alp*psi%e
     w%ln = alp*psi%ln
     call ihelmp(visc%p,w,psi,-alp,bet)
     !
     w%e = alp*chi%e
     w%ln = alp*chi%ln
     call ihelmp(visc%p,w,chi,-alp,bet)

     call deallocate( w )

   endif

   return
   end subroutine

   
   subroutine diagnost(psi,chi)
! ----------------------------------------------------
!  diagnostics
!  psi : toroidal component
!  chi : poloidal component
!  call in f/f space
!  call this routine when tim%n=0 to initialize 
! ----------------------------------------------------
  type(scalar):: psi,chi
  type(scalar):: rur,rup,uzh
  type(scalar):: ror,rop,ozh
  complex(p8),dimension(nxchopdim):: psiat1,chiat1
  real(p8):: ene1,ene2,valvor,valzmo,valang,valene,hel1,valhel
  real(p8):: pat1max,pat1min,cat1max,cat1min
  integer:: nn,mm,kk


! ------------
! INITIALIZED?
! ------------
  !> this subroutine must be called once while tim%n = 0

  ! tim is of "timedata" data-type, where n is an integer.
  ! diagv is of "diagdata" data-type; iniflag is an int.

  if(tim%n.ne.0 .and. diagv%iniflag.ne.100) then
    print *,'diag: WARNING:initialize incomplete !!'
  endif

! ----------
! NORM CHECK
! ----------

  ! He is doing the following:
  ! He takes the m=0, k=0 modes and sees what their imaginary part is.
  ! The arrays start at 1 (not 0: ndimr-1), hence *,1,1 corresponds 
  ! to the m=0 k=0 data.
  ! However the first piece of the legendre polynomial series is a constant
  ! hence he checks to see if the psi and chi are real for the 
  ! second piece in the r-series.. hence (2,1,1)

  if(abs(aimag(psi%e(2,1,1))+aimag(chi%e(2,1,1))).gt.1e-18) then
    write(6,*) '### diag: norm incomplete !!!'
    write(6,*) '### psi%e(2,1,1)=',psi%e(2,1,1)
    write(6,*) '### chi%e(2,1,1)=',chi%e(2,1,1)
    write(6,*) '###'
  endif
! ----------------
! DEALIASING CHECK
! ----------------
  if(mod(tim%n,1000).eq.20) then
    if(maxval(real(psi%e(nzchops(1)+1:,1,:))).ne.0.0 .or.  &
       maxval(real(psi%e(nzchops(2)+1:,2,:))).ne.0.0 ) then
       print *,'diag:dealiasing incomplete'
       call spy(psi%e(:,1,:))
       stop
    endif
  endif
! ------------------------
! CALCULATION BLOWING UP ?
! ------------------------
  do kk=1,min(5,nxchopdim)
  do nn=1,min(5,nzchop)
  do mm=2,min(6,ntchop)
    if(abs(psi%e(nn,mm,kk)).gt.100.) then
       write(6,*) '### CALCULATION DIVERGING !!!'
       write(6,60) tim%t,tim%n
       write(6,*) '(nn,mm,kk,psi)=',nn,mm,kk,psi%e(nn,mm,kk)
       !call msave(psi,files%psi(files%n))
       !call msave(chi,files%chi(files%n))
       stop
     endif
  enddo
  enddo
  enddo
! -----------
! DIAGNOSTICS 
! -----------
  if(mod(tim%n,diagv%it1).eq.0) then
    call allocate( rur )
    call allocate( rup )
    call allocate( uzh )
    call allocate( ror )
    call allocate( rop )
    call allocate( ozh )
    call chopset(3)
    call pc2vel(psi,chi,rur,rup,uzh,'h')
    call pc2vor(psi,chi,ror,rop,ozh,'h')

    !> vorticity 
    valvor= integ(ozh)
     
    !> z-direction momentum 
    valzmo= integ(uzh)
     
    !> angular momentum 
    rur=psi
    rur%ln=0
    call xxdx(rur,rup)
    valang = -integ(rup)

    !> energy p157
    rup=psi
    rup%ln=2*rup%ln
    call delsqh(rup,rur)
    call rtran(rur,1)
    rup=psi

    rup%ln=0
    call rtran(rup,1)
    ene1= prodct(rur,rup)
    !
    call del2(chi,rup)
    call rtran(rup,1)
    call delsqh(chi,rur)
    call rtran(rur,1)
    ene2= prodct(rur,rup)
    valene = -ene1+ene2

    !> helicity p157
    call delsqh(psi,rur)
    call rtran(rur,1)
    hel1 =prodct(rur,rup)
    valhel = 2*hel1
    call deallocate( rur )
    call deallocate( rup )
    call deallocate( uzh )

    !> store the initial values
    if(tim%n .eq.0) then
	diagv%iniflag = 100
        diagv%vori = valvor 
	diagv%zmomi = valzmo
        diagv%enei = valene
	diagv%heli = valhel
	diagv%angi = valang
    endif
 
    print *,'<<< diagnostics >>>'
    write(6,60) tim%t,tim%n
    60 format('tim%t=',e13.4,' tim%n=',i8)
    write(6,62) adv%x,adv%y
    62 format('adv%x=',e13.5,'  adv%y=',e13.5)
    70 format(a12,'=',e23.15, '  fractional change=',e23.15)

    psiat1= calcat1(psi)
    chiat1= calcat1(chi)
    pat1max = maxval(abs(psiat1))
    pat1min = minval(abs(psiat1))
    cat1max = maxval(abs(chiat1))
    cat1min = minval(abs(chiat1))
    if(pat1max-pat1min.gt.1.0e-10 .or. cat1max-cat1min.gt.1.0e-10) then
      write(6,81) 'psiat1(kk=1)=',abs(psiat1(1)),pat1max-pat1min
      write(6,81) 'chiat1(kk=1)=',abs(chiat1(1)),cat1max-cat1min
      81 format(a,e18.9, ' maxmin=',e18.9)
    endif

    if(diagv%vori.ne.0.0) then
      write(6,70) 'vorticity',valvor,(valvor-diagv%vori)/diagv%vori
    else
      write(6,70) 'vorticity',valvor
    endif

    if(diagv%zmomi.ne.0.0) then
      write(6,70) 'z-momentum',valzmo,(valzmo-diagv%zmomi)/diagv%zmomi
    else
      write(6,70) 'z-momentum',valzmo
    endif

    write(6,70) 'energy',valene,-(diagv%enei-valene)/diagv%enei
    if(diagv%heli.ne.0.0) then
      write(6,70) 'helicity',valhel,-(diagv%heli-valhel)/diagv%heli
    else
      write(6,70) 'helicity',valhel
    endif

    if(diagv%angi.ne.0.0) then
      write(6,70) 'integ_angmom ',valang,-(diagv%angi-valang)/diagv%angi
    else
      write(6,70) 'integ_angmom ',valang
    endif

    call allocate(rur)
    call mulxm(ozh,rur)
    call mulxm(rur,ozh)
    call deallocate(rur)
    call calcxy(ror,rop,ozh)
    call chopset(-3)
    call deallocate( ror )
    call deallocate( rop )
    call deallocate( ozh )
 end if

 return
end subroutine diagnost



  subroutine calcxy(ror,rop,oz)
! ----------------------------------------------------
! compute the location of the vortex (<x>,<y>) by weighting x and y
! with enstrophy. This routine is for monitoring purpose and does not
! affect the nonlinear computation. 
! Frequency of the monitoring is controled by variables in /diagnos/. 
! This routine is inserted here to save the time for transforms that
! would happen if this routine were called from subroutine diag.
! or,op,ozh in physical space
! or,op,ozh returns intact
! ----------------------------------------------------
  type(scalar):: ror,rop,oz
  real(p8):: valens,valx,valy,valr

  open(unit=14,file='coresize.dat',position='APPEND')

  call tofp(ror)
  call tofp(rop)
  call tofp(oz)
  call absvec2(ror,rop,oz,oz)

  ror=oz
  call toff(ror)
  valens=integ(ror)

  ror=oz
  call mulrcosth(ror)
  call toff(ror)
  valx=integ(ror)/valens

  ror=oz
  call mulrsinth(ror)
  call toff(ror)
  valy=integ(ror)/valens

  ror=oz
  call mulr2(ror)
  call toff(ror)
  valr=sqrt(integ(ror)/valens)

  write(6,60) tim%t,valx,valy,valr
  write(14,'(4e14.6)') tim%t,valx,valy,valr
  close(14)
  60 format('t= ',e10.4,' <x>= ',e10.4,' <y>= ',e10.4,' <r>= ',e10.4)
  
  if(tim%n.eq.0) then
    adv%t0=tim%t-tim%dt
    adv%x0=valx
    adv%y0=valy
    adv%t1=tim%t
    adv%x1=adv%x0
    adv%y1=adv%y0
  else
    adv%t0=adv%t1
    adv%x0=adv%x1
    adv%y0=adv%y1
    adv%t1=tim%t
    adv%x1=valx
    adv%y1=valy
  endif
  return
  end subroutine
!
!
!
    real(p8) function rcosth(r,th,z)
! ----------------------------------------------------
    real(p8):: r,th,z,x,y,r2
    rcosth = r*cos(th)
    return
    end function
!
!
!
   subroutine mulrcosth(a)
! ----------------------------------------------------
   type(scalar):: a
   real(p8):: r,z,thr,thi,aaa,bbb
   integer:: kk,nn,mm
   call tofp(a)
   do kk=1,nx
   do nn=1,nz
   do mm=1,nth
     z  = tfm%z(kk)
     r  = tfm%r(nn)
     thr = tfm%thr(mm)
     thi = tfm%thi(mm)
     aaa = r*cos(thr)*real(a%e(nn,mm,kk))
     bbb = r*cos(thi)*aimag(a%e(nn,mm,kk))
     a%e(nn,mm,kk)= cmplx(aaa,bbb,p8)
   enddo
   enddo
   enddo
   return
   end subroutine
!
!
!
   subroutine mulrsinth(a)
! ----------------------------------------------------
   type(scalar):: a
   real(p8):: r,z,thr,thi,aaa,bbb
   integer:: kk,nn,mm
   call tofp(a)
   do kk=1,nx
   do nn=1,nz
   do mm=1,nth
     z  = tfm%z(kk)
     r  = tfm%r(nn)
     thr = tfm%thr(mm)
     thi = tfm%thi(mm)
     aaa = r*sin(thr)*real(a%e(nn,mm,kk))
     bbb = r*sin(thi)*aimag(a%e(nn,mm,kk))
     a%e(nn,mm,kk)= cmplx(aaa,bbb,p8)
   enddo
   enddo
   enddo
   return
   end subroutine
!
!
!
   subroutine mulr2(a)
! ----------------------------------------------------
   type(scalar):: a
   real(p8):: r,z,thr,thi,aaa,bbb
   integer:: kk,nn,mm
   call tofp(a)
   do kk=1,nx
   do nn=1,nz
   do mm=1,nth
     z  = tfm%z(kk)
     r  = tfm%r(nn)
     thr = tfm%thr(mm)
     thi = tfm%thi(mm)
     aaa = r**2*real(a%e(nn,mm,kk))
     bbb = r**2*aimag(a%e(nn,mm,kk))
     a%e(nn,mm,kk)= cmplx(aaa,bbb,p8)
   enddo
   enddo
   enddo
   return
   end subroutine
!
!
!
    subroutine absvec2(u,v,w,a)
!----------------------------------------------------
!   this is custom routine for calcxy.
!   call in physical space
!   u,v,w in the form ( r*u_r, u*u_phi, u_z )
!   a can be one of u,v,w
!   (1-x)**2 division made
!----------------------------------------------------
    type(scalar):: u,v,w,a
    integer:: i,mm,kk
    real(p8):: aaa,bbb
    real(p8):: rur,rvr,wr
    real(p8):: rui,rvi,wi,r,f

    do mm=1,nth
    do kk=1,nx
    do i=1,nz
      r  = tfm%r(i)
      f = (1-tfm%x(i))**2
      rur=real(u%e(i,mm,kk))
      rui=aimag(u%e(i,mm,kk))
      rvr=real(v%e(i,mm,kk))
      rvi=aimag(v%e(i,mm,kk))
       wr=real(w%e(i,mm,kk))
       wi=aimag(w%e(i,mm,kk))
      aaa = ( (rur**2+rvr**2)/r**2 + wr**2)
      bbb = ( (rui**2+rvi**2)/r**2 + wi**2)
      a%e(i,mm,kk) = cmplx(aaa,bbb,p8)/f
    enddo
    enddo
    enddo

    return
    end subroutine
!
!
!
  subroutine freeadj
! ----------------------------------------------------
!  adjust freestream velocity such that the center of
!  vorticity distribution is at the center of the 
!  computational domain
!  in: ozh : (axial vorticity)/(1-x)**2
! ----------------------------------------------------
  real(p8):: freedt,vx,vy
  freedt = adv%int*tim%dt
  vx = (adv%x1-adv%x0)/(adv%t1-adv%t0)
  vy = (adv%y1-adv%y0)/(adv%t1-adv%t0)
  adv%ux = adv%ux - adv%x1/freedt - vx
  adv%uy = adv%uy - adv%y1/freedt - vy
  !
  nadd%u = sqrt(adv%ux**2+adv%uy**2)
  if(nadd%u.eq.0.0) then
    nadd%ang = 0.0
  else
     !got rid of to compile
         if(adv%uy >= 0) then
        nadd%ang = acos(adv%ux/nadd%u)
     else 
        nadd%ang = -1.0*acos(adv%ux/nadd%u)
     endif


     !nadd%ang = sign(1.0,adv%uy)*acos(adv%ux/nadd%u)
  endif

  write(6,60) tim%n,tim%t
  60 format('freeadj: tim%n=',i5,' tim%t=',f10.6)
  write(6,61) adv%ux,adv%uy,adv%uz
  61 format('freeadj: adv%ux=',e18.10,'  uy=',e18.10,' uz=',e18.10)
  write(6,62) nadd%u,nadd%ang
  62 format('freeadj: nadd%u=',e18.10,' ang=',e18.10)

  return
  end subroutine
!
!
!
  subroutine remove(psi,chi)
! kill the vortices which moves away from the computational domain
! ----------------------------------------------------
  use fd
  type(scalar):: psi,chi,wk
  integer:: ns,ns0,mm,kk,i
  real(p8):: ln
  if(rmv%sw.eq.1) then
    call allocate( wk )
    ns = nz*2/3
    ns0 = min(ns+4,nz)

    call delsqh(psi,wk)
    call rtran(wk,1)
    wk%e(ns0:nz,2:ntchop,1:nx)=0
    if (nx>1) wk%e(ns0:nz,1,2:nx)=0
    do mm=1,ntchop
    do kk=1,nx
      if(kk.gt.nxchop .and. kk.lt.nxchopu) cycle
      wk%e(ns:nz,mm,kk)= wk%e(ns:nz,mm,kk)*tfm%w(ns:nz)
      do i=1,5
        wk%e(ns:nz,mm,kk)=smooth(wk%e(ns:nz,mm,kk))
      enddo
      wk%e(ns:nz,mm,kk)= wk%e(ns:nz,mm,kk)/tfm%w(ns:nz)
    enddo
    enddo
    call rtran(wk,-1)
    call idelsqh(wk,psi)

    call chopset(2)
    call del2(chi,wk)
    call rtran(wk,1)
    wk%e(ns0:nz,2:ntchop,1:nx)=0
    if (nx>1) wk%e(ns0:nz,1,2:nx)=0
    do mm=1,ntchop
    do kk=1,nx
      if(kk.gt.nxchop .and. kk.lt.nxchopu) cycle
      wk%e(ns:nz,mm,kk)= wk%e(ns:nz,mm,kk)*tfm%w(ns:nz) &
		      /(1-tfm%x(ns:nz))**2
      do i=1,5
        wk%e(ns:nz,mm,kk)=smooth(wk%e(ns:nz,mm,kk))
      enddo
      wk%e(ns:nz,mm,kk)= wk%e(ns:nz,mm,kk)/tfm%w(ns:nz) &
		      *(1-tfm%x(ns:nz))**2
    enddo
    enddo
    call chopset(-2)
    call rtran(wk,-1)
    ln=chi%ln
    call idel2(wk,chi,ln)

    write(6,60) tim%n,tfm%r(ns)
    60 format('remove: tim%n=',i5,'  r >=',f6.2)

    call deallocate( wk )
  endif

  return
  end subroutine
!! !
!! !
!! !
!!   subroutine remove(psi,chi)
!! ! kill the vortices which moves away from the computational domain
!! ! ----------------------------------------------------
!!   type(scalar):: psi,chi,wk
!!   integer:: ns
!!   real(p8):: ln
!!   if(rmv%sw.eq.1) then
!!     call allocate( wk )
!! !   ns = nz*3/4
!!     ns = nz*2/3
!! 
!!     call delsqh(psi,wk)
!!     call tofp(wk)
!!     wk%e(ns:nz,1:nth,1:nx) = 0
!!     call toff(wk)
!!     call idelsqh(wk,psi)
!! 
!!     call chopset(2)
!!     call del2(chi,wk)
!!     call tofp(wk)
!!     wk%e(ns:nz,1:nth,1:nx) = 0
!!     call chopset(-2)
!!     call toff(wk)
!!     call idel2ln(wk,chi)
!! 
!!     write(6,60) tim%n,tfm%r(ns)
!!     60 format('remove: tim%n=',i5,'  r >=',f6.2)
!! 
!!     call deallocate( wk )
!!   endif
!! 
!!   return
!!   end subroutine
!
!
!

   function prodctm(a,b)
!----------------------------------------------------
!  calculate the product and integrate over the domain
!  for each azimuthal wavenumber.
!  call in r-physical / phi,z-fourier space
!  what's integrated is f=a*b*(1-mu)^2
!----------------------------------------------------
   type(scalar):: a,b
   complex(p8):: prod(nz,ntchop)
   real(p8),dimension(ntchop):: prodctm
   integer:: mm,kk
       
   if(a%space.ne.pff_space .or. b%space.ne.pff_space) then
     print *,'prodct:not in pff_space'
     print *,'a%space,b%space=',a%space,b%space
     stop
   endif
   if(a%ln.ne.0.0 .or. b%ln.ne.0.0) then
     print *,'prodct:logterm not zero'
   endif

   prod=0
   do kk = 1,nx
     if(kk.gt.nxchop .and. kk.lt.nxchopu) cycle
     prod(:,1) = prod(:,1)+a%e(:nz,1,kk)*conjg(b%e(:nz,1,kk))
   enddo

   do kk = 1,nx
     if(kk.gt.nxchop .and. kk.lt.nxchopu) cycle
     do mm = 2,ntchop
        prod(:,mm) = prod(:,mm)+2*(real(a%e(:nz,mm,kk))*real(b%e(:nz,mm,kk)) &
                      + aimag(a%e(:nz,mm,kk))*aimag(b%e(:nz,mm,kk)))
     enddo
   enddo

   do mm=1,ntchop
     prodctm(mm) = sum((prod(:,mm)*tfm%w)*tfm%pf(1,1,1))
     prodctm(mm) = 4*pi*zlen*eel2*prodctm(mm)*tfm%norm(1,1)
   enddo

   return
   end function
!
!
!
   function prodctk(a,b)
!----------------------------------------------------
!  calculate the product and integrate over the domain
!  for each axial wavenumber.
!  call in r-physical / phi,z-fourier space
!  what's integrated is f=a*b*(1-mu)^2
!----------------------------------------------------
   type(scalar):: a,b
   complex(p8):: prod(nz,nxchop)
   real(p8),dimension(nxchop):: prodctk
   integer:: mm,kk,kh,ck,cm
       
   if(a%space.ne.pff_space .or. b%space.ne.pff_space) then
     print *,'prodct:not in pff_space'
     print *,'a%space,b%space=',a%space,b%space
     stop
   endif
   if(a%ln.ne.0.0 .or. b%ln.ne.0.0) then
     print *,'prodct:logterm not zero'
   endif

   prod=0
   do kk = 1,nxchop
     if(kk.eq.1) then
       ck=2
       kh=1
     else
       ck=1
       kh=nx-kk+2
     endif
     do mm = 1,ntchop
       cm=1
       if(mm.eq.1) cm=2
       prod(:,kk)=prod(:,kk)+1.0_p8/ck/cm* & 
	       2*real(a%e(:nz,mm,kk)*conjg(b%e(:nz,mm,kk)) &
	       +b%e(:nz,mm,kh)*conjg(a%e(:nz,mm,kh)))
     enddo
   enddo

   do kk=1,nxchop
     prodctk(kk) = sum((prod(:,kk)*tfm%w)*tfm%pf(1,1,1))
     prodctk(kk) = 4*pi*zlen*eel2*prodctk(kk)*tfm%norm(1,1)
   enddo

   return
   end function


   subroutine enemon(psi,chi,em,ek)
!  compute energy for each azimuthal and axial wave number
!  -------------------------------------------------------
   type(scalar):: psi,chi, w1,w2
   real(p8):: em(ntchop),ek(nxchop)

   call allocate(w1)
   call allocate(w2)

    !> energy p157
    w1=psi
    w1%ln=2*psi%ln
    call delsqh(w1,w2)
    call rtran(w2,1)
    w1=psi
    w1%ln=0
    call rtran(w1,1)
    em= prodctm(w2,w1)
    ek= prodctk(w2,w1)
    !
    call del2(chi,w1)
    call rtran(w1,1)
    call delsqh(chi,w2)
    call rtran(w2,1)
    em= -em+prodctm(w2,w1)
    ek= -ek+prodctk(w2,w1)

   call deallocate(w1)
   call deallocate(w2)

   end subroutine


   subroutine hypadj(psi,chi,nup)
   ! adjust hyperviscosity from m-energy spectrum
   type(scalar):: psi,chi
   real(p8),dimension(:),allocatable:: em,ek
   real(p8):: nup,summ,sumk,fac,dx,x,sum,sumo
   real(p8),dimension(:),allocatable:: spm,spk
   integer:: is,ie,nm,nk,i
   real(p8),save:: summo,sumko,nupb
   integer,save:: first=0

   allocate( em(ntchop), ek(nxchop) )
   call chopset(2)
   call enemon(psi,chi,em,ek)
   call chopset(-2)

   is = size(em)/3
   ie = size(em)-2
   nm = ie-is+1
   allocate(spm(nm))
   
   spm = log10(em(is:ie))
   summ = 0
   is = size(ek)/3
   ie = size(ek)
   nk = ie-is+1
   allocate(spk(nk))
   spk = log10(ek(is:ie))
   sumk = 0

   !> dot product with legendre polynomial P_2(x)
   dx = 2.0_p8/(nm-1)
   do i=1,nm
     x  = -1+2.0_p8*(i-1)/(nm-1)
     fac = 1
     if(i.eq.1 .or. i.eq.nm) fac=0.5_p8
     summ=summ+ fac*spm(i)*dx*sqrt(5.0_p8/8)*(3*x**2-1)
   enddo
   dx = 2.0_p8/(nk-1)
   do i=1,nk
     x  = -1+2.0_p8*(i-1)/(nk-1)
     fac = 1
     if(i.eq.1 .or. i.eq.nk) fac=0.5_p8
     sumk=sumk+ fac*spk(i)*dx*sqrt(5.0_p8/8)*(3*x**2-1)
   enddo

   if(first.eq.0) then
     first=1
     summo=summ
     sumko=sumk
     nupb = 0.2/(1.5_p8*nzchop**2/eel**2)**(visc%p/2)
   endif

   if(summ.gt.sumk) then
     sumo = (3*summo+sumko)/4
     sum  = (3*summ+sumk)/4
   else
     sumo = (summo+3*sumko)/4
     sum  = (summ+3*sumk)/4
   endif
   !sum = 2*(sum-sumo) +sum + 0.3
   sum = 2*(sum-sumo) +sum + 0.2
   sum = min(3.0_p8,max(-0.5_p8,1.5*sum))
   nup = (1+sum)*nup
   if(nup.lt.nupb) nup=nupb
   write(6,60) tim%t,nup,summ,sumk
   60 format('hypadj:t= ',1pe10.4,' nup= ',1pe10.4,' summ= ',1pe10.3, & 
             ' sumk= ',1pe10.3)
   summo=summ
   sumko=sumk

   deallocate( em, ek )

 end subroutine hypadj

  !--------------
end module march
 !--------------
 
