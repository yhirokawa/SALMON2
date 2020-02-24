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

! for single-scale Maxwell-TDDFT method
module fdtd_coulomb_gauge
  implicit none

contains

subroutine fdtd_singlescale(itt,lg,mg,ng,system,info,info_field,rho,Vh,j_e,fg,poisson,srg_ng,Ac,div_Ac,fw)
  use structures
  use math_constants,only : zi,pi
  use phys_constants, only: cspeed_au
  use salmon_global, only: dt
  use sendrecv_grid, only: update_overlap_real8
  use stencil_sub, only: calc_gradient_field
  use communication, only: comm_is_root, comm_summation
  use inputoutput, only: t_unit_time
  use timer
  implicit none
  integer                 ,intent(in) :: itt
  type(s_rgrid)           ,intent(in) :: lg,mg,ng
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_field_parallel)  ,intent(in) :: info_field
  type(s_scalar)          ,intent(in) :: rho,Vh ! electron number density & Hartree potential
  type(s_vector)          ,intent(in) :: j_e    ! electron number current density (without rho*A/c)
  type(s_reciprocal_grid) ,intent(in) :: fg
  type(s_poisson)         ,intent(in) :: poisson
  type(s_sendrecv_grid)               :: srg_ng
  type(s_vector)                      :: Ac     ! A/c, A: vector potential, c: speed of light
  type(s_scalar)                      :: div_Ac ! div(A/c)
  type(s_singlescale)                 :: fw     ! FDTD working arrays, etc.
  !
  integer,parameter :: mstep=100
  integer :: ix,iy,iz,i1,ii,krd(3,3),lcs(3,3,3),dr(3)
  real(8) :: Hvol,hgs(3),dt_m,tm,coef,lap_A,Energy_em,diff_A,coef2 &
  & ,e_em,e_em_wrk,e_joule,e_joule_wrk,e_poynting(2),e_poynting_wrk(2),rho_t
  real(8),dimension(3) :: out_curr,out_Aext,out_Ab1,out_Ab2,wrk,wrk2,wrk3,wrk4,vec_je,Aext0,Aext1,Aext0_old,Aext1_old
  real(8) :: e_poy1,e_poy2,rtmp1(6),rtmp2(6)

  call timer_begin(LOG_SS_FDTD_CALC)

  krd = 0
  krd(1,1) = 1; krd(2,2) = 1; krd(3,3) = 1

  lcs = 0
  lcs(1,2,3) = 1; lcs(3,1,2) = 1; lcs(2,3,1) = 1
  lcs(1,3,2) = -1; lcs(2,1,3) = -1; lcs(3,2,1) = -1

  hgs = system%hgs
  Hvol = system%Hvol
  dt_m = dt / dble(mstep)

!-----------------------------------------------------------------------------------------------------------------------------------

  fw%box = 0d0
  wrk = 0d0
!$OMP parallel do collapse(2) private(ix,iy,iz,vec_je,rho_t) reduction(+:wrk)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    vec_je = ( j_e%v(1:3,ix,iy,iz) + fw%vec_je_old(1:3,ix,iy,iz) )*0.5d0 ! j_e(t) = ( j_e(t+dt/2) + j_e(t-dt/2) )/2
    rho_t  = ( rho%f(ix,iy,iz)     + fw%rho_old(ix,iy,iz)        )*0.5d0 ! rho(t) = ( rho(t+dt/2) + rho(t-dt/2) )/2
    fw%curr(ix,iy,iz,1:3) = vec_je + rho_t * fw%vec_Ac_m(1,ix,iy,iz,1:3) ! curr(t): electron number current density
    wrk = wrk + vec_je ! definition of out_curr?

    fw%box(ix,iy,iz) = Vh%f(ix,iy,iz) ! Vh(t+dt/2)
  end do
  end do
  end do
  wrk = wrk/dble(lg%num(1)*lg%num(2)*lg%num(3))
  call timer_end(LOG_SS_FDTD_CALC)

  call timer_begin(LOG_SS_FDTD_COMM_COLL)
  call comm_summation(wrk,out_curr,3,info_field%icomm_all)
  call timer_end(LOG_SS_FDTD_COMM_COLL)

! gradient of d(Vh)/dt (Vh: Hartree potential)
  call timer_begin(LOG_SS_FDTD_COMM)
  call update_overlap_real8(srg_ng, ng, fw%box)
  call timer_end(LOG_SS_FDTD_COMM)

  call timer_begin(LOG_SS_FDTD_CALC)
  call calc_gradient_field(ng,fw%coef_nab,fw%box,fw%grad_Vh) ! grad[Vh(t+dt/2)]

  !$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    fw%current4pi(ix,iy,iz,1:3) = ( fw%grad_Vh(1:3,ix,iy,iz) - fw%grad_Vh_old(1:3,ix,iy,iz) ) /dt & ! d(grad(Vh))/dt
                                & - 4d0*pi * fw%curr(ix,iy,iz,1:3) ! 4*pi* (electric current density)
  end do
  end do
  end do

!-----------------------------------------------------------------------------------------------------------------------------------

  ! FDTD loop: A(t) --> A(t+dt)

  fw%rotation_A = 0d0
  fw%divergence_A = 0d0
  do ii=1,mstep

  !$OMP parallel do collapse(2) private(ix,iy,iz)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)

      fw%vec_Ac_m(-1,ix,iy,iz,1:3) = fw%vec_Ac_m(0,ix,iy,iz,1:3)
      fw%vec_Ac_m(0 ,ix,iy,iz,1:3) = fw%vec_Ac_m(1,ix,iy,iz,1:3)

    end do
    end do
    end do

    do i1=1,3

  !$OMP parallel do collapse(2) private(ix,iy,iz)
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        fw%box(ix,iy,iz) = fw%vec_Ac_m(0,ix,iy,iz,i1)
      end do
      end do
      end do
      call timer_end(LOG_SS_FDTD_CALC)

      call timer_begin(LOG_SS_FDTD_COMM)
      call update_overlap_real8(fw%srg_eg, ng, fw%box)
      call timer_end(LOG_SS_FDTD_COMM)

      call timer_begin(LOG_SS_FDTD_CALC)
      if(ng%is(3)==lg%is(3))then
  !$OMP parallel do collapse(2) private(ix,iy,iz)
        do iy=ng%is(2),ng%ie(2)
        do ix=ng%is(1),ng%ie(1)
          fw%box(ix,iy,lg%is(3)-1) = fw%vec_Ac_boundary_bottom(ix,iy,i1)
        end do
        end do
      end if
      if(ng%ie(3)==lg%ie(3))then
  !$OMP parallel do collapse(2) private(ix,iy,iz)
        do iy=ng%is(2),ng%ie(2)
        do ix=ng%is(1),ng%ie(1)
          fw%box(ix,iy,lg%ie(3)+1) = fw%vec_Ac_boundary_top(ix,iy,i1)
        end do
        end do
      end if

  !$OMP parallel do collapse(2) private(ix,iy,iz,lap_A)
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        lap_A = ( - 2d0* fw%box(ix,iy,iz) + fw%box(ix-1,iy,iz) + fw%box(ix+1,iy,iz) ) / Hgs(1)**2 &
              + ( - 2d0* fw%box(ix,iy,iz) + fw%box(ix,iy-1,iz) + fw%box(ix,iy+1,iz) ) / Hgs(2)**2 &
              + ( - 2d0* fw%box(ix,iy,iz) + fw%box(ix,iy,iz-1) + fw%box(ix,iy,iz+1) ) / Hgs(3)**2
        fw%vec_Ac_m(1,ix,iy,iz,i1) = ( cspeed_au * dt_m )**2 * lap_A &
                                  + 2.d0* fw%box(ix,iy,iz) - fw%vec_Ac_m(-1,ix,iy,iz,i1) &
                                  + dt_m**2 * fw%current4pi(ix,iy,iz,i1)
      end do
      end do
      end do

  !   rotation & divergence of A
      if(ii==mstep/2) then
  !$OMP parallel do collapse(2) private(ix,iy,iz,wrk,dr,diff_A)
        do iz=ng%is(3),ng%ie(3)
        do iy=ng%is(2),ng%ie(2)
        do ix=ng%is(1),ng%ie(1)
        ! rot(A)
          wrk(1) = ( fw%box(ix+1,iy,iz) - fw%box(ix-1,iy,iz) ) / ( 2d0* Hgs(1) )
          wrk(2) = ( fw%box(ix,iy+1,iz) - fw%box(ix,iy-1,iz) ) / ( 2d0* Hgs(2) )
          wrk(3) = ( fw%box(ix,iy,iz+1) - fw%box(ix,iy,iz-1) ) / ( 2d0* Hgs(3) )
          fw%rotation_A(:,ix,iy,iz) = fw%rotation_A(:,ix,iy,iz) + lcs(:,1,i1) * wrk(1) + lcs(:,2,i1) * wrk(2) + lcs(:,3,i1) * wrk(3)

        ! div(A)
          dr = krd(:,i1)
          diff_A = ( fw%box(ix+dr(1),iy+dr(2),iz+dr(3)) - fw%box(ix-dr(1),iy-dr(2),iz-dr(3)) ) / ( 2d0* Hgs(i1) )
          fw%divergence_A(ix,iy,iz) = fw%divergence_A(ix,iy,iz) + diff_A
        end do
        end do
        end do
      end if

    end do ! i1 (spacial )

  ! external field
    tm = ( dble(itt-1) + dble(ii-1)/dble(mstep) ) *dt
    call pulse(tm,0d0,   Aext0_old)
    call pulse(tm,Hgs(3),Aext1_old)
    call pulse(tm+dt_m,0d0,   Aext0)
    call pulse(tm+dt_m,Hgs(3),Aext1)
  !$OMP parallel do collapse(2) private(ix,iy)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      fw%vec_Ac_ext_old(ix,iy,0,1:3) = Aext0_old
      fw%vec_Ac_ext_old(ix,iy,1,1:3) = Aext1_old
      fw%vec_Ac_ext(ix,iy,0,1:3) = Aext0
      fw%vec_Ac_ext(ix,iy,1,1:3) = Aext1
    end do
    end do
    out_Aext = Aext1

  ! z axis: Mur absorbing boundary condition
    coef = ( cspeed_au * dt_m - Hgs(3) ) / ( cspeed_au * dt_m + Hgs(3) )
    if(ng%is(3)==lg%is(3))then
  !$OMP parallel do collapse(2) private(ix,iy,iz)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
      ! absorbing boundary condition with the incident field vec_Ac_ext
        fw%vec_Ac_boundary_bottom(ix,iy,1:3) = fw%vec_Ac_ext(ix,iy,0,1:3) &
                                        + ( fw%vec_Ac_m(0,ix,iy,lg%is(3),1:3) - fw%vec_Ac_ext_old(ix,iy,1,1:3) )  &
                                        + coef* ( ( fw%vec_Ac_m(1,ix,iy,lg%is(3),1:3) - fw%vec_Ac_ext(ix,iy,1,1:3) ) &
                                                - ( fw%vec_Ac_boundary_bottom_old(ix,iy,1:3) - fw%vec_Ac_ext_old(ix,iy,0,1:3) ) )
      end do
      end do
    end if
    if(ng%ie(3)==lg%ie(3))then
  !$OMP parallel do collapse(2) private(ix,iy,iz)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        fw%vec_Ac_boundary_top(ix,iy,1:3) = fw%vec_Ac_m(0,ix,iy,lg%ie(3),1:3)   &
                                    + coef* ( fw%vec_Ac_m(1,ix,iy,lg%ie(3),1:3) - fw%vec_Ac_boundary_top_old(ix,iy,1:3) )
      end do
      end do
    end if

  !$OMP parallel do collapse(2) private(ix,iy,iz)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      fw%vec_Ac_boundary_bottom_old(ix,iy,1:3) = fw%vec_Ac_boundary_bottom(ix,iy,1:3)
      fw%vec_Ac_boundary_top_old   (ix,iy,1:3) = fw%vec_Ac_boundary_top   (ix,iy,1:3)
    end do
    end do

  end do ! ii=1,mstep

!-----------------------------------------------------------------------------------------------------------------------------------

  coef = cspeed_au / (4d0*pi)
  coef2 = Hvol / (8d0*pi)
  e_em_wrk = 0d0 ! for Electro-Magnetic energy
  e_joule_wrk = 0d0 ! for Joule dissipated power

  fw%lgbox1 = 0d0
  fw%vbox = 0d0

!$OMP parallel do collapse(2) private(ix,iy,iz,wrk,wrk2,wrk3,wrk4) reduction(+:e_em_wrk,e_joule_wrk)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)

    fw%vbox(1:3,ix,iy,iz) = ( fw%vec_Ac_m(1,ix,iy,iz,1:3) + fw%vec_Ac_old(1:3,ix,iy,iz) ) * 0.5d0 ! ( A(t+dt) + A(t) )/2
    fw%lgbox1(ix,iy,iz)   = fw%divergence_A(ix,iy,iz)

    wrk4 = ( fw%vec_Ac_m(1,ix,iy,iz,:) - fw%vec_Ac_old(:,ix,iy,iz) ) / dt ! (A(t+dt)-A(t))/dt
    wrk  = - (-1d0)*fw%grad_Vh(:,ix,iy,iz) - wrk4 ! E
    wrk2 = cspeed_au * fw%rotation_A(:,ix,iy,iz)  ! B
    wrk3 = - fw%curr(ix,iy,iz,1:3)                ! j
    fw%poynting_vector(:,ix,iy,iz) = coef * ( lcs(:,1,2) * wrk(1) * wrk2(2) + lcs(:,1,3) * wrk(1) * wrk2(3) &
                                            + lcs(:,2,1) * wrk(2) * wrk2(1) + lcs(:,2,3) * wrk(2) * wrk2(3) &
                                            + lcs(:,3,1) * wrk(3) * wrk2(1) + lcs(:,3,2) * wrk(3) * wrk2(2) ) ! E x B
    e_em_wrk = e_em_wrk + coef2 * ( sum(wrk**2) + sum(wrk2**2) ) ! ( E^2 + B^2 )/(8*pi)
    e_joule_wrk = e_joule_wrk + sum(wrk3*wrk) * Hvol ! j*E

  end do
  end do
  end do
  call timer_end(LOG_SS_FDTD_CALC)

  call timer_begin(LOG_SS_FDTD_COMM_COLL)
  call comm_summation(fw%vbox,  fw%vec_Ac,3*lg%num(1)*lg%num(2)*lg%num(3),info_field%icomm_all)
  call comm_summation(fw%lgbox1,fw%lgbox2,  lg%num(1)*lg%num(2)*lg%num(3),info_field%icomm_all)
  call comm_summation(e_em_wrk,e_em,info_field%icomm_all)
  call comm_summation(e_joule_wrk,e_joule,info_field%icomm_all)
  call timer_end(LOG_SS_FDTD_COMM_COLL)

  call timer_begin(LOG_SS_FDTD_CALC)
!$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
    Ac%v(:,ix,iy,iz)   = fw%vec_Ac(:,ix,iy,iz) ! Ac(t+dt/2)
    div_Ac%f(ix,iy,iz) = fw%lgbox2(ix,iy,iz)
  end do
  end do
  end do

! stock old Ac, grad_Vh, j_e, & rho

!$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    fw%vec_Ac_old(:,ix,iy,iz)    = fw%vec_Ac_m(1,ix,iy,iz,1:3) ! Ac(t+dt) --> Ac(t) of next step
    fw%grad_Vh_old(1:3,ix,iy,iz) = fw%grad_Vh(1:3,ix,iy,iz) ! grad[Vh(t-dt/2)] of next step
    fw%vec_je_old(1:3,ix,iy,iz)  = j_e%v(1:3,ix,iy,iz) ! j_e(t-dt/2) of next step
    fw%rho_old(ix,iy,iz)         = rho%f(ix,iy,iz)     ! rho(t-dt/2) of next step
  end do
  end do
  end do

!-----------------------------------------------------------------------------------------------------------------------------------

  ! integral(A) @ z = 0 (bottom boundary)
  wrk = 0d0
  if (ng%is(3) == lg%is(3)) then
!$omp parallel do private(ix,iy,iz) reduction(+:wrk)
    do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        wrk(1) = wrk(1) + fw%vec_Ac(1,ix,iy,ng%is(3))
        wrk(2) = wrk(2) + fw%vec_Ac(2,ix,iy,ng%is(3))
        wrk(3) = wrk(2) + fw%vec_Ac(3,ix,iy,ng%is(3))
      end do
    end do
  end if

  ! integral(A) @ z = az (top boundary)
  wrk3 = 0d0
  if (ng%ie(3) == lg%ie(3)) then
!$omp parallel do private(ix,iy,iz) reduction(+:wrk3)
    do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        wrk3(1) = wrk3(1) + fw%vec_Ac(1,ix,iy,ng%ie(3))
        wrk3(2) = wrk3(2) + fw%vec_Ac(2,ix,iy,ng%ie(3))
        wrk3(3) = wrk3(2) + fw%vec_Ac(3,ix,iy,ng%ie(3))
      end do
    end do
  end if
  call timer_end(LOG_SS_FDTD_CALC)

  call timer_begin(LOG_SS_FDTD_COMM_COLL)
  rtmp1 = [wrk, wrk3]
  call comm_summation(rtmp1,rtmp2,6,info_field%icomm_xy)
  wrk (1:3) = rtmp2(1:3)
  wrk3(1:3) = rtmp2(4:6)
  call timer_end(LOG_SS_FDTD_COMM_COLL)

  call timer_begin(LOG_SS_FDTD_CALC)
  out_Ab1 = wrk  / dble(lg%num(1)*lg%num(2))
  out_Ab2 = wrk3 / dble(lg%num(1)*lg%num(2))

! Surface integral of the poynting vector S
  e_poy1 = 0d0
  e_poy2 = 0d0
  coef = Hgs(1)*Hgs(2)
  if(ng%is(3)==lg%is(3)) then ! integral(S) @ z = 0 (bottom boundary)
!$omp parallel do private(ix,iy,iz) reduction(+:e_poy1)
    do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        e_poy1 = e_poy1 + fw%poynting_vector(3,ix,iy,lg%ie(3)) * coef
      end do
    end do
  end if
  if(ng%ie(3)==lg%ie(3)) then ! integral(S) @ z = az (top boundary)
!$omp parallel do private(ix,iy,iz) reduction(+:e_poy2)
    do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        e_poy2 = e_poy2 + fw%poynting_vector(3,ix,iy,lg%ie(3)) * coef
      end do
    end do
  end if
  call timer_end(LOG_SS_FDTD_CALC)

  call timer_begin(LOG_SS_FDTD_COMM_COLL)
  e_poynting_wrk = [e_poy1, e_poy2]
  call comm_summation(e_poynting_wrk,e_poynting,2,info_field%icomm_xy)
  call timer_end(LOG_SS_FDTD_COMM_COLL)

  call timer_begin(LOG_SS_FDTD_CALC)
  Energy_em = e_em
  fw%Energy_joule = fw%Energy_joule + dt*e_joule
  fw%Energy_poynting = fw%Energy_poynting + dt*e_poynting

  if(comm_is_root(info%id_rko)) write(fw%fh_rt_micro,'(99(1X,E23.15E3))') &
    dble(itt)*dt*t_unit_time%conv,out_Ab1,out_Ab2,out_Aext,out_curr,fw%E_electron,fw%Energy_poynting,Energy_em,fw%Energy_joule

! FIXME: is following result unused?
! for spatial distribution of excitation energy
!  coef = Hgs(1)*Hgs(2)
!  fw%integral_poynting_tmp = 0d0
!!$omp parallel do collapse(2) private(iz,iy,ix)
!  do iy=ng%is(2),ng%ie(2)
!  do ix=ng%is(1),ng%ie(1)
!    do iz=ng%is(3),ng%ie(3)
!      fw%integral_poynting_tmp(iz) = fw%integral_poynting_tmp(iz) + fw%poynting_vector(3,ix,iy,iz) * coef
!    end do
!  end do
!  end do
!  call timer_end(LOG_SS_FDTD_CALC)
!
!  call timer_begin(LOG_SS_FDTD_COMM_COLL)
!  call comm_summation(fw%integral_poynting_tmp,fw%integral_poynting_tmp2,lg_num(3),info_field%icomm_all)
!  call timer_end(LOG_SS_FDTD_COMM_COLL)
!
!  call timer_begin(LOG_SS_FDTD_CALC)
!  fw%integral_poynting = fw%integral_poynting + dt * fw%integral_poynting_tmp2

  fw%Ac_zt_t = 0d0
! for the vector potential Ax(z,t)
  do iz=ng%is(3),ng%ie(3)
!$omp parallel do collapse(2) private(iy,ix)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
        fw%Ac_zt_t(1,iz) = fw%Ac_zt_t(1,iz) + fw%vec_Ac(1,ix,iy,iz) / (lg%num(1)*lg%num(2))
        fw%Ac_zt_t(2,iz) = fw%Ac_zt_t(2,iz) + fw%vec_Ac(2,ix,iy,iz) / (lg%num(1)*lg%num(2))
        fw%Ac_zt_t(3,iz) = fw%Ac_zt_t(3,iz) + fw%vec_Ac(3,ix,iy,iz) / (lg%num(1)*lg%num(2))
    end do
    end do
  end do
  call timer_end(LOG_SS_FDTD_CALC)

  call timer_begin(LOG_SS_FDTD_COMM_COLL)
  call comm_summation(fw%Ac_zt_t,fw%Ac_zt,size(fw%Ac_zt),info_field%icomm_all)
  call timer_end(LOG_SS_FDTD_COMM_COLL)

  if(comm_is_root(info%id_rko)) then
    do iz=lg%is(3),lg%ie(3)
      write(fw%fh_Ac_zt,fmt='(99(1X,E23.15E3))',advance='no') dble(iz)*hgs(3),fw%Ac_zt(1,iz),fw%Ac_zt(2,iz),fw%Ac_zt(3,iz)
    end do
    write(fw%fh_Ac_zt,'()')
  end if

!-----------------------------------------------------------------------------------------------------------------------------------

  return

contains

  subroutine pulse(t,r,A_ext)
    use em_field, only: calc_Ac_ext
    implicit none
    real(8),intent(in)  :: t,r
    real(8),intent(out) :: A_ext(3)
    !
    real(8) :: tt

    tt = t - r/cspeed_au
    call calc_Ac_ext(tt,A_ext)

    return
  end subroutine pulse

end subroutine fdtd_singlescale

!===================================================================================================================================

subroutine init_singlescale(comm,ng,mg,lg,hgs,rho,Vh,srg_ng,fw)
  use structures
  use sendrecv_grid, only: update_overlap_real8
  use stencil_sub, only: calc_gradient_field
  use salmon_global, only: sysname,base_directory,yn_restart
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root
  use initialization_sub, only: set_bn
  use filesystem, only: open_filehandle
  use inputoutput, only: t_unit_time
  use checkpoint_restart_sub, only: restart_singlescale
  implicit none
  integer       ,intent(in) :: comm
  type(s_rgrid) ,intent(in) :: lg,mg,ng
  real(8)       ,intent(in) :: hgs(3)
  type(s_scalar),intent(in) :: rho,Vh ! electron number density & Hartree potential
  type(s_sendrecv_grid)     :: srg_ng
  type(s_singlescale)       :: fw
  !
  character(100) :: filename
  integer :: ii,jj,ix,iy,iz
  real(8) :: bnmat(4,4)

  fw%Energy_poynting = 0d0
  fw%Energy_joule = 0d0

  call set_bn(bnmat)
  do jj=1,3
    do ii=1,4
      fw%coef_nab(ii,jj) = bnmat(ii,4)/hgs(jj)
    end do
  end do

  allocate( fw%vec_Ac     (3,lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)) )
  allocate( fw%vec_Ac_old (3,ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)) )
  allocate( fw%curr         (ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3),3) )
  allocate( fw%vec_je_old (3,ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)) )
  allocate( fw%rho_old      (ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)) )
  allocate( fw%current4pi   (ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3),3) )
  allocate( fw%grad_Vh    (3,ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)) )
  allocate( fw%grad_Vh_old(3,ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)) )

!1st element: time step (-1->m-1, 0->m, 1->m+1)
!5th element: components of A vector (1->Ax, 2->Ay, 3->Az)
  allocate( fw%vec_Ac_m(-1:1,ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3),1:3) )

  allocate( fw%vec_Ac_boundary_bottom(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),1:3) &
           ,fw%vec_Ac_boundary_bottom_old(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),1:3) &
           ,fw%vec_Ac_boundary_top(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),1:3) &
           ,fw%vec_Ac_boundary_top_old(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),1:3) )

  allocate(fw%integral_poynting(lg%is(3):lg%ie(3)))

  allocate(fw%rotation_A(3,ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)) &
   & ,fw%poynting_vector(3,ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)) &
        & ,fw%divergence_A(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)) &
        & ,fw%vbox      (3,lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)) &
        & ,fw%lgbox1      (lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)) &
        & ,fw%lgbox2      (lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)) &
        & ,fw%integral_poynting_tmp(lg%num(3)),fw%integral_poynting_tmp2(lg%num(3)) &
        & ,fw%vec_Ac_ext    (ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),0:1,1:3) &
        & ,fw%vec_Ac_ext_old(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),0:1,1:3) &
        & ,fw%Ac_zt(3,lg%is(3):lg%ie(3)) &
        & ,fw%Ac_zt_t(3,lg%is(3):lg%ie(3)))

  fw%vec_Ac = 0d0
  fw%vec_Ac_old = 0d0
  fw%vec_Ac_m = 0d0
  fw%vec_Ac_boundary_bottom = 0d0
  fw%vec_Ac_boundary_bottom_old = 0d0
  fw%vec_Ac_boundary_top = 0d0
  fw%vec_Ac_boundary_top_old = 0d0
  fw%integral_poynting = 0d0
  fw%curr = 0d0
  fw%vec_je_old = 0d0
  fw%Ac_zt = 0d0
  fw%Ac_zt_t = 0d0

  if(comm_is_root(nproc_id_global)) then
    write(filename,"(2A,'_rt_micro.data')") trim(base_directory),trim(SYSname)
    fw%fh_rt_micro = open_filehandle(filename)
    write(fw%fh_rt_micro, '("#",99(1X,I0,":",A,"[",A,"]"))') &
      & 1, "time", trim(t_unit_time%name), &
      & 2, "Ac_tot_x(z=0)", "a.u.", &
      & 3, "Ac_tot_y(z=0)", "a.u.", &
      & 4, "Ac_tot_z(z=0)", "a.u.", &
      & 5, "Ac_tot_x(z=L)", "a.u.", &
      & 6, "Ac_tot_y(z=L)", "a.u.", &
      & 7, "Ac_tot_z(z=L)", "a.u.", &
      & 8, "Ac_ext_x",      "a.u.", &
      & 9, "Ac_ext_y",      "a.u.", &
     & 10, "Ac_ext_z",      "a.u.", &
     & 11, "J_x(w/o Ac)",   "a.u.", &
     & 12, "J_y(w/o Ac)",   "a.u.", &
     & 13, "J_z(w/o Ac)",   "a.u.", &
     & 14, "E_electron",    "a.u.", &
     & 15, "E_poynting(z=0)", "a.u.", &
     & 16, "E_poynting(z=L)", "a.u.", &
     & 17, "E_em",            "a.u.", &
     & 18, "E_joule",         "a.u."

  ! for spatial distribution of excitation energy
    write(filename,"(2A,'_excitation.data')") trim(base_directory),trim(SYSname)
    fw%fh_excitation = open_filehandle(filename)

  ! for the vector potential Ac(z,t)
    write(filename,"(2A,'_Ac_zt.data')") trim(base_directory),trim(SYSname)
    fw%fh_Ac_zt = open_filehandle(filename)
  end if
  
  allocate(fw%box(ng%is_array(1):ng%ie_array(1), &
  & ng%is_array(2):ng%ie_array(2), &
  & ng%is_array(3):ng%ie_array(3)))
  fw%box = 0d0

!$OMP parallel do collapse(2) private(ix,iy,iz)
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    fw%box(ix,iy,iz) = Vh%f(ix,iy,iz)
    fw%rho_old(ix,iy,iz) = rho%f(ix,iy,iz)
  end do
  end do
  end do
  call update_overlap_real8(srg_ng, ng, fw%box)
  call calc_gradient_field(ng,fw%coef_nab,fw%box,fw%grad_Vh_old)
  
  if(yn_restart=='y') then
    call restart_singlescale(comm,lg,ng,fw)
  end if

  return
end subroutine init_singlescale

!===================================================================================================================================

! line integral \int_{{\bf r}_0 \rightarrow {\bf r}_1} {\bf A}({\bf r}) \cdot d{\bf x}
! path: r0 --> r1 = (ix1*Hx,iy1*Hy,iz1*Hz)
subroutine line_integral(integral,r0,A,nx,ny,nz,ix1,iy1,iz1,Hx,Hy,Hz &
                        ,A_lerp,line,wrk,n_max,index)
  implicit none
  integer,intent(in)  :: nx,ny,nz,ix1,iy1,iz1,n_max
  real(8),intent(in)  :: r0(3),A(3,0:nx-1,0:ny-1,0:nz-1),Hx,Hy,Hz
  real(8),intent(out) :: integral
  integer             :: index(n_max)
  real(8)             :: A_lerp(3,n_max),line(3,n_max),wrk(n_max)
  !
  integer :: ix,iy,iz,ix0,iy0,iz0,i,j,k,n,ixp,iyp,izp
  real(8) :: r1(3),p(3),q(3),r(3),h(3),l(3),t

# define IMOD(ix,iy,iz) modulo(ix,nx),modulo(iy,ny),modulo(iz,nz)

  r1 = (/ ix1*Hx, iy1*Hy, iz1*Hz /)

  ix0 = floor( r0(1)/Hx )
  iy0 = floor( r0(2)/Hy )
  iz0 = floor( r0(3)/Hz )

  h = (/ Hx, Hy, Hz /)

! trilinear interpolation @ r0
  ix = ix0
  iy = iy0
  iz = iz0
  ixp = ix + 1
  iyp = iy + 1
  izp = iz + 1
  p = (/ ix*Hx, iy*Hy, iz*Hz /)
  r = ( r0 - p )/h
  A_lerp(:,1) = A(:,IMOD(ix,iy,iz)) * ( 1d0 - r(1) ) * ( 1d0 - r(2) ) * ( 1d0 - r(3) ) &
              + A(:,IMOD(ixp,iy,iz)) * r(1) * ( 1d0 - r(2) ) * ( 1d0 - r(3) ) &
              + A(:,IMOD(ix,iyp,iz)) * ( 1d0 - r(1) ) * r(2) * ( 1d0 - r(3) ) &
              + A(:,IMOD(ix,iy,izp)) * ( 1d0 - r(1) ) * ( 1d0 - r(2) ) * r(3) &
              + A(:,IMOD(ix,iyp,izp)) * ( 1d0 - r(1) ) * r(2) * r(3) &
              + A(:,IMOD(ixp,iy,izp)) * r(1) * ( 1d0 - r(2) ) * r(3) &
              + A(:,IMOD(ixp,iyp,iz)) * r(1) * r(2) * ( 1d0 - r(3) ) &
              + A(:,IMOD(ixp,iyp,izp)) * r(1) * r(2) * r(3)
  line(:,1) = 0d0
  wrk(1) = 0d0

! bilinear interpolation @ intersections
  i = 1
  do ix=min(ix0,ix1),max(ix0,ix1)
    if(ix0==ix1) cycle
    t = ( ix*Hx - r0(1) ) / ( r1(1) - r0(1) )
    if(t<0d0 .or. 1d0<t) cycle
    l = r0 + t * ( r1 - r0 )
    iy = floor( l(2)/Hy )
    iz = floor( l(3)/Hz )
    iyp = iy + 1
    izp = iz + 1
    q = (/ ix*Hx, iy*Hy, iz*Hz /)
    r = ( l - q )/h
    i = i + 1
    if(i>n_max) stop "n_max is too small"
    A_lerp(:,i) = A(:,IMOD(ix,iy,iz)) * ( 1d0 - r(2) ) * ( 1d0 - r(3) ) &
                + A(:,IMOD(ix,iyp,iz)) * r(2) * ( 1d0 - r(3) ) &
                + A(:,IMOD(ix,iy,izp)) * ( 1d0 - r(2) ) * r(3) &
                + A(:,IMOD(ix,iyp,izp)) * r(2) * r(3)
    line(:,i) = l - r0
    wrk(i) = line(1,i)**2 + line(2,i)**2 + line(3,i)**2
  end do
  do iy=min(iy0,iy1),max(iy0,iy1)
    if(iy0==iy1) cycle
    t = ( iy*Hy - r0(2) ) / ( r1(2) - r0(2) )
    if(t<0d0 .or. 1d0<t) cycle
    l = r0 + t * ( r1 - r0 )
    ix = floor( l(1)/Hx )
    iz = floor( l(3)/Hz )
    ixp = ix + 1
    izp = iz + 1
    q = (/ ix*Hx, iy*Hy, iz*Hz /)
    r = ( l - q )/h
    i = i + 1
    if(i>n_max) stop "n_max is too small"
    A_lerp(:,i) = A(:,IMOD(ix,iy,iz)) * ( 1d0 - r(1) ) * ( 1d0 - r(3) ) &
                + A(:,IMOD(ixp,iy,iz)) * r(1) * ( 1d0 - r(3) ) &
                + A(:,IMOD(ix,iy,izp)) * ( 1d0 - r(1) ) * r(3) &
                + A(:,IMOD(ixp,iy,izp)) * r(1) * r(3)
    line(:,i) = l - r0
    wrk(i) = line(1,i)**2 + line(2,i)**2 + line(3,i)**2
  end do
  do iz=min(iz0,iz1),max(iz0,iz1)
    if(iz0==iz1) cycle
    t = ( iz*Hz - r0(3) ) / ( r1(3) - r0(3) )
    if(t<0d0 .or. 1d0<t) cycle
    l = r0 + t * ( r1 - r0 )
    ix = floor( l(1)/Hx )
    iy = floor( l(2)/Hy )
    ixp = ix + 1
    iyp = iy + 1
    q = (/ ix*Hx, iy*Hy, iz*Hz /)
    r = ( l - q )/h
    i = i + 1
    if(i>n_max) stop "n_max is too small"
    A_lerp(:,i) = A(:,IMOD(ix,iy,iz)) * ( 1d0 - r(1) ) * ( 1d0 - r(2) ) &
                + A(:,IMOD(ixp,iy,iz)) * r(1) * ( 1d0 - r(2) ) &
                + A(:,IMOD(ix,iyp,iz)) * ( 1d0 - r(1) ) * r(2) &
                + A(:,IMOD(ixp,iyp,iz)) * r(1) * r(2)
    line(:,i) = l - r0
    wrk(i) = line(1,i)**2 + line(2,i)**2 + line(3,i)**2
  end do
  if(i==1) then
    i = i + 1
    if(i>n_max) stop "n_max is too small"
    A_lerp(:,i) = A(:,IMOD(ix0,iy0,iz0))
    line(:,i) = p - r0
    wrk(i) = line(1,i)**2 + line(2,i)**2 + line(3,i)**2
  end if
  n = i

  call heapsort(n,wrk,index)

  integral = 0d0
  do i=2,n
    j = index(i)
    k = index(i-1)
    integral = integral + 0.5d0* ( A_lerp(1,j) + A_lerp(1,k) ) &
                        * ( line(1,j) - line(1,k) ) &
                        + 0.5d0* ( A_lerp(2,j) + A_lerp(2,k) ) &
                        * ( line(2,j) - line(2,k) ) &
                        + 0.5d0* ( A_lerp(3,j) + A_lerp(3,k) ) &
                        * ( line(3,j) - line(3,k) )
  end do

  return
contains

  subroutine heapsort(n,array,index)
    implicit none
    integer,intent(in)    :: n
    real(8),intent(inout) :: array(n)
    integer,intent(out)   :: index(n)
    !
    integer :: i,k,j,l,m
    real(8) :: t

    if(n < 2) return

    do i=1,n
      index(i) = i
    end do

    l = n/2 + 1
    k = n
    do while(k /= 1)
      if(l > 1)then
        l = l - 1
        t = array(l)
        m = index(l)
      else
        t = array(k)
        m = index(k)
        array(k) = array(1)
        index(k) = index(1)
        k = k - 1
        if(k == 1) then
          array(1) = t
          index(1) = m
          exit
        end if
      end if
      i = l
      j = l + l
      do while(j <= k)
        if(j < k) then
          if(array(j) < array(j+1)) j = j + 1
        endif
        if (t < array(j))then
          array(i) = array(j)
          index(i) = index(j)
          i = j
          j = j + j
        else
          j = k + 1
        end if
      end do
      array(i) = t
      index(i) = m
    end do

    return
  end subroutine heapsort
end subroutine line_integral

end module fdtd_coulomb_gauge
