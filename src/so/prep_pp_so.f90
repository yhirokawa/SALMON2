!
!  Copyright 2019 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module prep_pp_so_sub

  use spin_orbit_global, only: SPIN_ORBIT_ON

  implicit none

  private
  public :: calc_uv_so
  public :: SPIN_ORBIT_ON

  integer,allocatable :: lma_tbl_so(:,:)
  integer :: Nlma_so

contains

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

  subroutine set_nlma_so( pp, ppg )
    use salmon_global,only : natom,kion
    use structures,only : s_pp_info,s_pp_grid
    implicit none 
    type(s_pp_info) :: pp
    type(s_pp_grid) :: ppg
    integer :: lma, lm, lm_max
    integer :: a,ik,m,l,ll,l0

write(*,*) "----------- set_nlma_so"

    lm_max=0
    lma=0
    do a=1,natom
      ik=kion(a)
      lm=0
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        if ( pp%inorm(l,ik) == 0 .and. pp%inorm_so(l,ik) == 0 ) cycle
        do m=-ll-1,ll
          lm=lm+1
          lma=lma+1
        end do
      end do
      l0=l
      end do
      lm_max=max(lm_max,lm)
    end do

    Nlma_so = lma

    write(*,*) "Nlma_so=", Nlma_so, ppg%nlma

    allocate( lma_tbl_so(lm_max,natom) ); lma_tbl_so=0
    allocate( ppg%ia_tbl_so(Nlma_so) ); ppg%ia_tbl_so=0

    call set_lma_tbl( lma_tbl_so, ppg%ia_tbl_so, pp )

    allocate( ppg%uv_so(ppg%nps,Nlma_so,2,2) ); ppg%uv_so=(0.0d0,0.0d0)
    allocate( ppg%duv_so(ppg%nps,Nlma_so,3,2,2) ); ppg%duv_so=(0.0d0,0.0d0)

  end subroutine set_nlma_so

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

  subroutine set_lma_tbl( lma_tbl_so, ia_tbl_so, pp )
    use salmon_global,only : natom, kion
    use structures,only : s_pp_info, s_pp_grid
    implicit none
    integer,intent(out) :: lma_tbl_so(:,:)
    integer,intent(out) :: ia_tbl_so(:)
    type(s_pp_info) :: pp
    integer :: lm,lma
    integer :: a,ik,m,l,l0,ll
write(*,*) "----------- set_lma_tbl"
    lma_tbl_so=0
    ia_tbl_so=0
    lma=0
    do a=1,natom
      ik=kion(a)
      lm=0
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        if ( pp%inorm(l,ik) == 0 .and. pp%inorm_so(l,ik) == 0 ) cycle
        do m=-ll-1,ll
          lm=lm+1
          lma=lma+1
          lma_tbl_so(lm,a)=lma
          ia_tbl_so(lma)=a
        end do
      end do
      l0=l
      end do
    end do !a
  end subroutine set_lma_tbl

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
  subroutine calc_uv_so(pp,ppg,lx,ly,lz,nl,hx,hy,hz,flag_use_grad_wf_on_force,property,hvol0)
  use salmon_global,  only : natom, kion, iperiodic, yn_domain_parallel
  use math_constants, only : pi
  use structures,     only : s_pp_info, s_pp_grid
  implicit none
  type(s_pp_info) :: pp
  type(s_pp_grid) :: ppg
  integer,intent(in) :: nl
  integer,intent(in) :: lx(nl),ly(nl),lz(nl)
  real(8),intent(in) :: hx,hy,hz
  logical,intent(in) :: flag_use_grad_wf_on_force
  character(17),intent(in) :: property
  real(8),intent(in),optional :: hvol0
  real(8) :: save_udvtbl_a(pp%nrmax,0:2*pp%lmax+1,natom)
  real(8) :: save_udvtbl_b(pp%nrmax,0:2*pp%lmax+1,natom)
  real(8) :: save_udvtbl_c(pp%nrmax,0:2*pp%lmax+1,natom)
  real(8) :: save_udvtbl_d(pp%nrmax,0:2*pp%lmax+1,natom)
  integer :: a,ik,j,l,lm,m,ll,l0
  integer :: ilma,intr,ir,j_angular_momentum
  real(8),allocatable :: xn(:),yn(:),an(:),bn(:),cn(:),dn(:)  
  real(8) :: dudvtbl_a(pp%nrmax,0:2*pp%lmax+1), dudvtbl_b(pp%nrmax,0:2*pp%lmax+1)
  real(8) :: dudvtbl_c(pp%nrmax,0:2*pp%lmax+1), dudvtbl_d(pp%nrmax,0:2*pp%lmax+1)
  real(8) :: uvr(0:2*pp%lmax+1),duvr(0:2*pp%lmax+1)
  real(8) :: r,x,y,z
  real(8) :: xx
  real(8) :: rshift(3)
  real(8) :: hvol, coef

  hvol=hx*hy*hz
  if ( present(hvol0) ) hvol = hvol0

  if ( iperiodic == 0 ) then
    if ( mod(lx(nl)-lx(1)+1,2) == 1 ) then
      rshift(1)=0.0d0
    else
      rshift(1)=-0.5d0*Hx
    end if
    if ( mod(ly(nl)-ly(1)+1,2) == 1 ) then
      rshift(2)=0.0d0
    else
      rshift(2)=-0.5d0*Hy
    end if
    if ( mod(lz(nl)-lz(1)+1,2) == 1 ) then
      rshift(3)=0.0d0
    else
      rshift(3)=-0.5d0*Hz
    end if
  else if ( iperiodic == 3 ) then 
    if ( yn_domain_parallel == 'y' ) then
      rshift(1)=-Hx
      rshift(2)=-Hy
      rshift(3)=-Hz
    else
      rshift(:)=0.0d0
    end if
  end if

  call set_nlma_so( pp, ppg )

  do j_angular_momentum = 1, 2 ![ 1:j=l+1/2, 2:j=l-1/2 ]

  if ( property /= 'update_wo_realloc' ) then
    do a=1,natom
      ik=kion(a)
      allocate( xn(0:pp%nrps(ik)-1),yn(0:pp%nrps(ik)-1),an(0:pp%nrps(ik)-2) &
               ,bn(0:pp%nrps(ik)-2),cn(0:pp%nrps(ik)-2),dn(0:pp%nrps(ik)-2) )
      xn(0:pp%nrps(ik)-1) = pp%radnl(1:pp%nrps(ik),ik)
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        select case( j_angular_momentum )
        case(1); yn(0:pp%nrps(ik)-1) = pp%udvtbl(1:pp%nrps(ik),l,ik)
        case(2); yn(0:pp%nrps(ik)-1) = pp%udvtbl_so(1:pp%nrps(ik),l,ik)
        end select
        call spline(pp%nrps(ik),xn,yn,an,bn,cn,dn)
        save_udvtbl_a(1:pp%nrps(ik)-1,l,a) = an(0:pp%nrps(ik)-2)
        save_udvtbl_b(1:pp%nrps(ik)-1,l,a) = bn(0:pp%nrps(ik)-2)
        save_udvtbl_c(1:pp%nrps(ik)-1,l,a) = cn(0:pp%nrps(ik)-2)
        save_udvtbl_d(1:pp%nrps(ik)-1,l,a) = dn(0:pp%nrps(ik)-2)
      end do
      l0=l
      end do
      deallocate(xn,yn,an,bn,cn,dn)
    end do !a
  end if
  
  do a=1,natom

    ik=kion(a)
    if ( .not.flag_use_grad_wf_on_force ) then !legacy for ion-force
      allocate( xn(0:pp%nrps(ik)-1),yn(0:pp%nrps(ik)-1),an(0:pp%nrps(ik)-2) &
               ,bn(0:pp%nrps(ik)-2),cn(0:pp%nrps(ik)-2),dn(0:pp%nrps(ik)-2) )
      xn(0:pp%nrps(ik)-1) = pp%radnl(1:pp%nrps(ik),ik)
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        select case( j_angular_momentum )
        case(1); yn(0:pp%nrps(ik)-1) = pp%dudvtbl(1:pp%nrps(ik),l,ik)
        case(2); yn(0:pp%nrps(ik)-1) = pp%dudvtbl_so(1:pp%nrps(ik),l,ik)
        end select
        call spline(pp%nrps(ik),xn,yn,an,bn,cn,dn)
        dudvtbl_a(1:pp%nrps(ik)-1,l) = an(0:pp%nrps(ik)-2)
        dudvtbl_b(1:pp%nrps(ik)-1,l) = bn(0:pp%nrps(ik)-2)
        dudvtbl_c(1:pp%nrps(ik)-1,l) = cn(0:pp%nrps(ik)-2)
        dudvtbl_d(1:pp%nrps(ik)-1,l) = dn(0:pp%nrps(ik)-2)
      end do
      l0=l
      end do
      deallocate(xn,yn,an,bn,cn,dn)
    end if
 
    do j=1,ppg%mps(a)

      x=ppg%rxyz(1,j,a)
      y=ppg%rxyz(2,j,a)
      z=ppg%rxyz(3,j,a)
      r=sqrt(x*x+y*y+z*z)+1d-50
      do ir=1,pp%nrps(ik)
        if ( pp%radnl(ir,ik) > r ) exit
      end do
      intr=ir-1
      if (intr < 0 .or. intr >= pp%nrps(ik) ) stop 'bad intr at prep_ps'
      xx = r - pp%radnl(intr,ik)

      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1
        uvr(l) = save_udvtbl_a(intr,l,a)*xx**3 + save_udvtbl_b(intr,l,a)*xx**2 &
               + save_udvtbl_c(intr,l,a)*xx    + save_udvtbl_d(intr,l,a)
        if ( .not.flag_use_grad_wf_on_force ) then !legacy for ion-force
          duvr(l) = dudvtbl_a(intr,l)*xx**3 + dudvtbl_b(intr,l)*xx**2 &
                  + dudvtbl_c(intr,l)*xx    + dudvtbl_d(intr,l)
        end if
      end do
      l0=l
      end do
 
      lm=0
      l0=0
      do ll=0,pp%mlps(ik)
      do l=l0,l0+pp%nproj(ll,ik)-1

        if ( pp%inorm(l,ik) == 0 .and. pp%inorm_so(l,ik) == 0 ) then
           write(*,*) "skil(l,ll)",l,ll
           cycle
        end if

        select case( j_angular_momentum )
        case( 1 )

          do m=-ll-1,ll
            lm=lm+1
            ilma=lma_tbl_so(lm,a)
!
! j=l+1/2, alpha spin
!
            coef=sqrt( dble(ll+m+1)/dble(2*ll+1) )
if ( a == 1 .and. j == 11 ) write(*,'("(1,1)",4i4,2f15.8)') l,ll,m,m,coef,abs(zylm(x,y,z,ll,m))
            ppg%uv_so(j,ilma,1,1)=coef*uvr(l)*zylm(x,y,z,ll,m)
            if ( .not.flag_use_grad_wf_on_force ) then !legacy for ion-force
              if ( r > 1.0d-6 ) then
                ppg%duv_so(j,ilma,1,1,1) = coef*(duvr(l)*(x/r)*zylm(x,y,z,ll,m)+uvr(l)*dzylm(x,y,z,ll,m,1))
                ppg%duv_so(j,ilma,2,1,1) = coef*(duvr(l)*(y/r)*zylm(x,y,z,ll,m)+uvr(l)*dzylm(x,y,z,ll,m,2))
                ppg%duv_so(j,ilma,3,1,1) = coef*(duvr(l)*(z/r)*zylm(x,y,z,ll,m)+uvr(l)*dzylm(x,y,z,ll,m,3))
              else
                ppg%duv_so(j,ilma,1,1,1) = coef*uvr(l)*dzylm(x,y,z,ll,m,1)
                ppg%duv_so(j,ilma,2,1,1) = coef*uvr(l)*dzylm(x,y,z,ll,m,2)
                ppg%duv_so(j,ilma,3,1,1) = coef*uvr(l)*dzylm(x,y,z,ll,m,3)
              end if
            end if
!
! j=l+1/2, beta spin
!
            coef=sqrt( dble(ll-m)/dble(2*ll+1) )
if ( a == 1 .and. j == 11 ) write(*,'("(1,2)",4i4,2f15.8)') l,ll,m,m+1,coef,abs(zylm(x,y,z,ll,m+1))
            ppg%uv_so(j,ilma,1,2)=coef*uvr(l)*zylm(x,y,z,ll,m+1)
            if ( .not.flag_use_grad_wf_on_force ) then !legacy for ion-force
              if ( r > 1.0d-6 ) then
                ppg%duv_so(j,ilma,1,1,2) = coef*(duvr(l)*(x/r)*zylm(x,y,z,ll,m+1)+uvr(l)*dzylm(x,y,z,ll,m+1,1))
                ppg%duv_so(j,ilma,2,1,2) = coef*(duvr(l)*(y/r)*zylm(x,y,z,ll,m+1)+uvr(l)*dzylm(x,y,z,ll,m+1,2))
                ppg%duv_so(j,ilma,3,1,2) = coef*(duvr(l)*(z/r)*zylm(x,y,z,ll,m+1)+uvr(l)*dzylm(x,y,z,ll,m+1,3))
              else
                ppg%duv_so(j,ilma,1,1,2) = coef*uvr(l)*dzylm(x,y,z,ll,m+1,1)
                ppg%duv_so(j,ilma,2,1,2) = coef*uvr(l)*dzylm(x,y,z,ll,m+1,2)
                ppg%duv_so(j,ilma,3,1,2) = coef*uvr(l)*dzylm(x,y,z,ll,m+1,3)
              end if
            end if

          end do !m

        case( 2 )

          do m=-ll,ll
            lm=lm+1
            ilma=lma_tbl_so(lm,a)
!
! j=l-1/2, alpha spin
!
            coef=sqrt( dble(ll-m+1)/dble(2*ll+1) )
if ( a == 1 .and. j == 11 ) write(*,'("(2,1)",4i4,2f15.8)') l,ll,m,m-1,coef,abs(zylm(x,y,z,ll,m-1))
            ppg%uv_so(j,ilma,2,1)=coef*uvr(l)*zylm(x,y,z,ll,m-1)
            if ( .not.flag_use_grad_wf_on_force ) then !legacy for ion-force
              if ( r > 1.0d-6 ) then
                ppg%duv_so(j,ilma,1,2,1) = coef*(duvr(l)*(x/r)*zylm(x,y,z,ll,m-1)+uvr(l)*dzylm(x,y,z,ll,m-1,1))
                ppg%duv_so(j,ilma,2,2,1) = coef*(duvr(l)*(y/r)*zylm(x,y,z,ll,m-1)+uvr(l)*dzylm(x,y,z,ll,m-1,2))
                ppg%duv_so(j,ilma,3,2,1) = coef*(duvr(l)*(z/r)*zylm(x,y,z,ll,m-1)+uvr(l)*dzylm(x,y,z,ll,m-1,3))
              else
                ppg%duv_so(j,ilma,1,2,1) = coef*uvr(l)*dzylm(x,y,z,ll,m-1,1)
                ppg%duv_so(j,ilma,2,2,1) = coef*uvr(l)*dzylm(x,y,z,ll,m-1,2)
                ppg%duv_so(j,ilma,3,2,1) = coef*uvr(l)*dzylm(x,y,z,ll,m-1,3)
              end if
            end if
!
! j=l-1/2, beta spin
!
            coef=-sqrt( dble(ll+m)/dble(2*ll+1) )
if ( a == 1 .and. j == 11 ) write(*,'("(2,2)",4i4,2f15.8)') l,ll,m,m,coef,abs(zylm(x,y,z,ll,m))
            ppg%uv_so(j,ilma,2,2)=coef*uvr(l)*zylm(x,y,z,ll,m)
            if ( .not.flag_use_grad_wf_on_force ) then !legacy for ion-force
              if ( r > 1.0d-6 ) then
                ppg%duv_so(j,ilma,1,2,2) = coef*(duvr(l)*(x/r)*zylm(x,y,z,ll,m)+uvr(l)*dzylm(x,y,z,ll,m,1))
                ppg%duv_so(j,ilma,2,2,2) = coef*(duvr(l)*(y/r)*zylm(x,y,z,ll,m)+uvr(l)*dzylm(x,y,z,ll,m,2))
                ppg%duv_so(j,ilma,3,2,2) = coef*(duvr(l)*(z/r)*zylm(x,y,z,ll,m)+uvr(l)*dzylm(x,y,z,ll,m,3))
              else
                ppg%duv_so(j,ilma,1,2,2) = coef*uvr(l)*dzylm(x,y,z,ll,m,1)
                ppg%duv_so(j,ilma,2,2,2) = coef*uvr(l)*dzylm(x,y,z,ll,m,2)
                ppg%duv_so(j,ilma,3,2,2) = coef*uvr(l)*dzylm(x,y,z,ll,m,3)
              end if
            end if

          end do !m

        end select  
      end do
      l0=l
      end do

    end do !j

  end do !a

  !lma=0
  !do a=1,natom
    !ik=kion(a)
    !if ( property /= 'update_wo_realloc' ) then
    !  l0=0
    !  do ll=0,pp%mlps(ik)
    !  do l=l0,l0+pp%nproj(ll,ik)-1
    !    if ( pp%inorm(l,ik) == 0 ) cycle
    !    do m=-ll,ll
    !      lma=lma+1
    !      ppg%rinv_uvu(lma)=dble(pp%inorm(l,ik))*hvol
    !    end do
    !  end do
    !  l0=l
    !  end do
    !end if
  !end do !a

  end do !j_angular_momentum

  end subroutine calc_uv_so


subroutine spline(Np,xn,yn,an,bn,cn,dn)
  integer,intent(in) :: Np
  real(8),intent(in) :: xn(0:Np-1),yn(0:Np-1)
  real(8),intent(out) :: an(0:Np-2),bn(0:Np-2),cn(0:Np-2),dn(0:Np-2)
  integer :: i,Npm2,info
  real(8) :: dxn(0:Np-1),dyn(0:Np-1),u(1:Np-2),v(1:Np-2),Amat(1:Np-2,1:Np-2)
  real(8) :: Amat_t(1:Np-2,1:Np-2)
! for lapack
  integer :: LWORK
  integer, allocatable :: IPIV(:) ! dimension N
  real(8), allocatable :: WORK(:) ! dimension LWORK
! for check inverse matrix problem
!  integer :: j,k
!  real(8) :: Amat_chk(1:Np-2,1:Np-2)
!  real(8) :: ss

  Npm2 = Np-2
  LWORK = Npm2*Npm2*6
  allocate(IPIV(Npm2),WORK(LWORK))


  do i = 0,Np-2
    dxn(i) = xn(i+1) - xn(i)
    dyn(i) = yn(i+1) - yn(i)
  end do

  do i = 1,Npm2
    v(i) = 6d0*(dyn(i)/dxn(i) - dyn(i-1)/dxn(i-1))
  end do

  Amat = 0d0
  Amat(1,1) = 2d0*(dxn(1) + dxn(0))
  Amat(1,2) = dxn(1)
  do i = 2,Npm2-1
    Amat(i,i+1) = dxn(i)
    Amat(i,i  ) = 2d0*(dxn(i)+dxn(i-1))
    Amat(i,i-1) = dxn(i-1)
  end do
  Amat(Npm2,Npm2  ) = 2d0*(dxn(Npm2)+dxn(Npm2-1))
  Amat(Npm2,Npm2-1) = dxn(Npm2-1)

! inverse matrix problem
  Amat_t = Amat


  call DGETRF(Npm2, Npm2, Amat_t, Npm2, IPIV, info)  ! factorize
  call DGETRI(Npm2, Amat_t, Npm2, IPIV, WORK, LWORK, info)  ! inverse

!  check inverse matrix problem
!  do i = 1,Npm2
!    do j = 1,Npm2
!      ss = 0d0
!      do k = 1,Npm2
!        ss = ss + Amat(i,k)*Amat_t(k,j)
!      end do
!      Amat_chk(i,j) = ss
!    end do
!  end do
!
!  do i = 1,Npm2
!    write(*,'(999e16.6e3)')(Amat_chk(i,j),j=1,Npm2)
!  end do
!
!  stop


  do i = 1,Npm2
    u(i) = sum(Amat_t(i,:)*v(:))
  end do

! for b
  bn(0) = 0d0
  bn(1:Np-2) = 0.5d0*u(1:Np-2)
! for a
  do i = 0,Npm2-1
    an(i) = (u(i+1) -2d0*bn(i))/(6d0*dxn(i))
  end do
  an(Npm2) = (0d0 -2d0*bn(Npm2))/(6d0*dxn(Npm2))
! for d
  dn(0:Npm2) = yn(0:Npm2)
! for c
  i=0
  cn(i) = dyn(i)/dxn(i) - dxn(i)*(u(i+1)+2d0*0.d0)/6d0
  do i = 1,Npm2-1
     cn(i) = dyn(i)/dxn(i) - dxn(i)*(u(i+1)+2d0*u(i))/6d0
  end do
  cn(Npm2) = dyn(Npm2)/dxn(Npm2) - dxn(Npm2)*(0d0+2d0*u(Npm2))/6d0

  return
end subroutine spline

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine bisection(xx,inode,iak,nr,rad_psl)
  use salmon_global,only : nelem
  implicit none
  integer,intent(out) :: inode
  integer,intent(in)  :: iak
  integer,intent(in)  :: nr
  real(8),intent(in)  :: rad_psl(nr,nelem)
  real(8),intent(in)  :: xx
  integer :: imin,imax
  
  imin=1
  imax=nr
  do while (imax-imin>1)
    inode=(imin+imax)/2
    if(xx>rad_psl(inode,iak))then
      imin=inode
    else
      imax=inode
    end if
  end do
  inode=imin

end subroutine bisection

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine init_uvpsi_summation(ppg,icomm_r)
  use structures,    only: s_pp_grid
  use salmon_global, only: natom
  use communication, only: comm_get_groupinfo &
                                 ,comm_allgather &
                                 ,comm_create_group_byid
  implicit none
  type(s_pp_grid),intent(inout) :: ppg
  integer,intent(in) :: icomm_r

  integer :: ilma,ia
  integer :: irank_r,isize_r,n,i
  integer :: nlma
  logical,allocatable :: ireferred_atom_comm_r(:,:)
  integer,allocatable :: iranklist(:)
stop "xxxx"
  call comm_get_groupinfo(icomm_r, irank_r, isize_r)

  nlma = ppg%Nlma

  allocate(iranklist(isize_r))
  allocate(ireferred_atom_comm_r(natom,isize_r))
  allocate(ppg%ireferred_atom(natom))
  allocate(ppg%icomm_atom(natom))
  allocate(ppg%irange_atom(2,natom))

  ppg%ireferred_atom = .false.
  do ilma=1,nlma
    ia = ppg%ia_tbl(ilma)
    ppg%ireferred_atom(ia) = ppg%ireferred_atom(ia) .or. (ppg%mps(ia) > 0)
  end do
  call comm_allgather(ppg%ireferred_atom, ireferred_atom_comm_r, icomm_r)

  ppg%irange_atom(1,:) = 1
  ppg%irange_atom(2,:) = 0
  do ia=1,natom
    ! forward search
    do ilma=1,nlma
      if (ppg%ia_tbl(ilma) == ia) then
        ppg%irange_atom(1,ia) = ilma
        exit
      end if
    end do

    ! backward search
    do ilma=nlma,1,-1
      if (ppg%ia_tbl(ilma) == ia) then
        ppg%irange_atom(2,ia) = ilma
        exit
      end if
    end do
  end do

  do ia=1,natom
    n = 0
    do i=1,isize_r
      if (ireferred_atom_comm_r(ia,i)) then
        n = n + 1
        iranklist(n) = i - 1
      end if
    end do

    ppg%icomm_atom(ia) = comm_create_group_byid(icomm_r, iranklist(1:n))
  end do
end subroutine init_uvpsi_summation

subroutine finalize_uvpsi_summation(ppg)
  use structures,    only: s_pp_grid
  use communication, only: comm_free_group
  implicit none
  type(s_pp_grid),intent(inout) :: ppg
  integer :: ia

  if (allocated(ppg%irange_atom))    deallocate(ppg%irange_atom)
  if (allocated(ppg%ireferred_atom)) deallocate(ppg%ireferred_atom)
  if (allocated(ppg%icomm_atom)) then
    do ia=1,size(ppg%icomm_atom)
      call comm_free_group(ppg%icomm_atom(ia))
    end do
    deallocate(ppg%icomm_atom)
  end if
end subroutine

  function zylm( x,y,z,l,m )
    use salmon_math, only: Ylm
    implicit none
    real(8),intent(in) :: x,y,z
    integer,intent(in) :: l,m
    complex(8) :: zylm
    real(8),parameter :: factor=0.70710678118654752d0 ! =sqrt(0.5)
    zylm=(0.0d0,0.0d0)
    if ( l < 0 .or. abs(m) > l ) return
    if ( m < 0 ) then
      zylm = factor*dcmplx( -Ylm(x,y,z,l,m), Ylm(x,y,z,l,-m) )
    else if ( m == 0 ) then
      zylm = Ylm(x,y,z,l,m)
    else if ( m > 0 ) then
      zylm = factor*dcmplx(  Ylm(x,y,z,l,m), Ylm(x,y,z,l,-m) )
    end if
  end function zylm

  function dzylm( x,y,z,l,m,idir )
    use salmon_math, only: dYlm
    implicit none
    real(8),intent(in) :: x,y,z
    integer,intent(in) :: l,m,idir
    complex(8) :: dzylm
    real(8),parameter :: factor=0.70710678118654752d0 ! =sqrt(0.5)
    dzylm=(0.0d0,0.0d0)
    if ( l < 0 .or. abs(m) > l ) return
    if ( m < 0 ) then
      dzylm = factor*dcmplx(-dYlm(x,y,z,l,m,idir), dYlm(x,y,z,l,-m,idir) )
    else if ( m == 0 ) then
      dzylm = dYlm(x,y,z,l,m,idir)
    else if ( m > 0 ) then
      dzylm = factor*dcmplx( dYlm(x,y,z,l,m,idir), dYlm(x,y,z,l,-m,idir) )
    end if
  end function dzylm

end module prep_pp_so_sub
