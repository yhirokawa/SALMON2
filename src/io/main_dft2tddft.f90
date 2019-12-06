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

#include "config.h"

subroutine main_dft2tddft
use structures
use salmon_global, only: ispin,cval,xc,xname,cname,directory_read_data,calc_mode, &
                         target_nproc_k,target_nproc_ob,target_nproc_domain_orbital, &
                         target_nproc_domain_general
use salmon_parallel, only: nproc_group_global
use salmon_communication, only: comm_get_globalinfo,comm_bcast,comm_is_root,comm_sync_all
use salmon_xc
use timer
use initialization_sub
use init_communicator
use checkpoint_restart_sub
use set_numcpu
use filesystem, only: create_directory
implicit none
integer :: Miter
character(100) :: file_atoms_coo

type(s_rgrid) :: lg_scf,lg_rt
type(s_rgrid) :: mg_scf,mg_rt
type(s_rgrid) :: ng_scf,ng_rt
type(s_process_info) :: pinfo_scf,pinfo_rt
type(s_orbital_parallel) :: info_scf,info_rt
type(s_field_parallel) :: info_field_scf,info_field_rt
type(s_sendrecv_grid) :: srg, srg_ng
type(s_orbital) :: spsi,shpsi,sttpsi
type(s_dft_system) :: system_scf,system_rt
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_scalar) :: srho,sVh,sVpsl
type(s_scalar),allocatable :: V_local(:),srho_s(:),sVxc(:)
type(s_reciprocal_grid) :: fg
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
type(s_ofile)  :: ofl

integer :: icomm,irank,nprocs
logical :: if_stop
character(256) :: dir_file_out,gdir,pdir
integer,parameter :: fh = 41

call comm_get_globalinfo(icomm,irank,nprocs)

if(comm_is_root(irank))then
  print *, '==============================================='
  print *, 'DFT2TDDFT'
  print *, '  DFT calulated data convert to calculate TDDFT'
  print *, '==============================================='
end if

call init_xc(xc_func, ispin, cval, xcname=xc, xname=xname, cname=cname)

call timer_begin(LOG_TOTAL)

call timer_begin(LOG_INIT_GS)

call convert_input_scf(file_atoms_coo)

! please move folloings into initialization_dft
call init_dft(nproc_group_global,pinfo_scf,info_scf,info_field_scf,lg_scf,mg_scf,ng_scf,system_scf,stencil,fg,poisson,srg,srg_ng,ofl)
allocate( srho_s(system_scf%nspin),V_local(system_scf%nspin),sVxc(system_scf%nspin) )

call initialization1_dft( system_scf, energy, stencil, fg, poisson,  &
                          lg_scf, mg_scf, ng_scf,  &
                          pinfo_scf, info_scf, info_field_scf,  &
                          srg, srg_ng,  &
                          srho, srho_s, sVh, V_local, sVpsl, sVxc,  &
                          spsi, shpsi, sttpsi,  &
                          pp, ppg, ppn,  &
                          ofl )

call generate_restart_directory_name(directory_read_data,gdir,pdir)
call read_bin(pdir,lg_scf,mg_scf,ng_scf,system_scf,info_scf,spsi,Miter)

call timer_end(LOG_INIT_GS)


! Redistributed write to use TDDFT calculation.
! ---------------------------------------------------------
call timer_begin(LOG_WRITE_GS_DATA)

calc_mode = 'RT' ! FIXME
call init_dft_system(lg_rt,system_rt,stencil)

pinfo_rt%npk                 = target_nproc_k
pinfo_rt%nporbital           = target_nproc_ob
pinfo_rt%npdomain_orbital    = target_nproc_domain_orbital
pinfo_rt%npdomain_general    = target_nproc_domain_general
if((target_nproc_ob + sum(target_nproc_domain_orbital) + sum(target_nproc_domain_general)) == 0) then
  if (system_rt%ngrid > 16**3) then
    call set_numcpu_general(iprefer_domain_distribution,system_rt%nk,system_rt%no,icomm,pinfo_rt)
  else
    call set_numcpu_general(iprefer_orbital_distribution,system_rt%nk,system_rt%no,icomm,pinfo_rt)
  end if
end if

if (comm_is_root(irank)) then
  if_stop = .not. check_numcpu(icomm, pinfo_rt)
end if
call comm_bcast(if_stop, icomm)
if (if_stop) stop 'fail: check_numcpu'

if (ispin==1) then
  pinfo_rt%nporbital_spin(1)=(pinfo_rt%nporbital+1)/2
  pinfo_rt%nporbital_spin(2)= pinfo_rt%nporbital   /2
else
  pinfo_rt%nporbital_spin = 0
end if
pinfo_rt%npdomain_general_dm(1:3)=pinfo_rt%npdomain_general(1:3)/pinfo_rt%npdomain_orbital(1:3)
pinfo_rt%npdomain_general_dm(1:3)=pinfo_rt%npdomain_general(1:3)/pinfo_rt%npdomain_orbital(1:3)

call init_communicator_dft(icomm,pinfo_rt,info_rt,info_field_rt)
call init_orbital_parallel_singlecell(system_rt,info_rt,pinfo_rt)
call init_grid_parallel(irank,nprocs,pinfo_rt,lg_rt,mg_rt,ng_rt)

call deallocate_orbital(shpsi)
call allocate_orbital_complex(system_rt%nspin,mg_rt,info_rt,shpsi)

call convert_wave_function

! TODO: move to checkpoint_restart.f90
call generate_restart_directory_name(ofl%dir_out_restart,gdir,pdir)
call create_directory(pdir)

if(comm_is_root(irank))then
!information
  dir_file_out = trim(pdir)//"info.bin"
  open(fh,file=dir_file_out,form='unformatted')
  write(fh) system_scf%nk
  write(fh) system_scf%no
  write(fh) Miter
  write(fh) nprocs
  close(fh)

!occupation
  dir_file_out = trim(pdir)//"occupation.bin"
  open(fh,file=dir_file_out,form='unformatted')
  write(fh) system_scf%rocc(1:system_scf%no,1:system_scf%nk,1:system_scf%nspin)
  close(fh)
end if

call write_wavefunction(pdir,lg_rt,mg_rt,system_rt,info_rt,shpsi,.false.)
! TODO: move to checkpoint_restart.f90

call timer_end(LOG_WRITE_GS_DATA)


call finalize_xc(xc_func)

call timer_end(LOG_TOTAL)

contains

subroutine convert_wave_function
use salmon_communication, only: comm_summation, comm_create_group_byid, &
                                comm_proc_null, comm_send, comm_recv, comm_bcast
implicit none
integer,allocatable    :: is_ref(:,:,:),is_ref_srank(:,:,:),is_ref_rrank(:,:,:)
integer,allocatable    :: irank_src(:,:),irank_dst(:,:)
real(8),allocatable    :: dbuf1(:,:,:,:),dbuf2(:,:,:,:)
complex(8),allocatable :: zbuf1(:,:,:,:),zbuf2(:,:,:,:)
integer :: ik,io,jrank,krank

! get referrer
allocate(is_ref      (system_scf%no,system_scf%nk,0:nprocs-1))
allocate(is_ref_srank(system_scf%no,system_scf%nk,0:nprocs-1)) ! SCF
allocate(is_ref_rrank(system_scf%no,system_scf%nk,0:nprocs-1)) ! RT

is_ref = 0
is_ref(info_scf%io_s:info_scf%io_e, &
       info_scf%ik_s:info_scf%ik_e, irank) = 1
call comm_summation(is_ref,is_ref_srank,size(is_ref),icomm)

is_ref = 0
is_ref(info_rt%io_s:info_rt%io_e, &
       info_rt%ik_s:info_rt%ik_e, irank) = 1
call comm_summation(is_ref,is_ref_rrank,size(is_ref),icomm)


! find src rank and dst rank
allocate(irank_src(system_scf%no,system_scf%nk))
allocate(irank_dst(system_scf%no,system_scf%nk))

irank_src = -1
irank_dst = -1

do ik=1,system_scf%nk
do io=1,system_scf%no
  ! src (SCF) rank
  do jrank=0,nprocs-1
    if (is_ref_srank(io,ik,jrank) >= 1) then
      irank_src(io,ik) = jrank
      exit
    end if
  end do

  ! dst (RT) rank
  do jrank=0,nprocs-1
    if (is_ref_rrank(io,ik,jrank) >= 1) then
      irank_dst(io,ik) = jrank
      exit
    end if
  end do
end do
end do


! conversion
if (allocated(spsi%rwf)) then
  allocate(dbuf1(lg_scf%is(1):lg_scf%ie(1), &
                 lg_scf%is(2):lg_scf%ie(2), &
                 lg_scf%is(3):lg_scf%ie(3),system_scf%nspin))
  allocate(dbuf2(lg_scf%is(1):lg_scf%ie(1), &
                 lg_scf%is(2):lg_scf%ie(2), &
                 lg_scf%is(3):lg_scf%ie(3),system_scf%nspin))
end if
if (allocated(spsi%zwf)) then
  allocate(zbuf1(lg_scf%is(1):lg_scf%ie(1), &
                 lg_scf%is(2):lg_scf%ie(2), &
                 lg_scf%is(3):lg_scf%ie(3),system_scf%nspin))
  allocate(zbuf2(lg_scf%is(1):lg_scf%ie(1), &
                 lg_scf%is(2):lg_scf%ie(2), &
                 lg_scf%is(3):lg_scf%ie(3),system_scf%nspin))
end if

do ik=1,system_scf%nk
do io=1,system_scf%no
  ! gather rgrid data
  if (info_scf%ik_s <= ik .and. ik <= info_scf%ik_e .and. &
      info_scf%io_s <= io .and. io <= info_scf%io_e) then
    if (allocated(dbuf1)) then
      dbuf1 = 0d0
      dbuf1(mg_scf%is(1):mg_scf%ie(1),mg_scf%is(2):mg_scf%ie(2),mg_scf%is(3):mg_scf%ie(3),:) &
        = spsi%rwf(mg_scf%is(1):mg_scf%ie(1),mg_scf%is(2):mg_scf%ie(2),mg_scf%is(3):mg_scf%ie(3),:,io,ik,1)
      call comm_summation(dbuf1,dbuf2,size(dbuf1),info_scf%icomm_r)
    end if
    if (allocated(zbuf1)) then
      zbuf1 = 0d0
      zbuf1(mg_scf%is(1):mg_scf%ie(1),mg_scf%is(2):mg_scf%ie(2),mg_scf%is(3):mg_scf%ie(3),:) &
        = spsi%zwf(mg_scf%is(1):mg_scf%ie(1),mg_scf%is(2):mg_scf%ie(2),mg_scf%is(3):mg_scf%ie(3),:,io,ik,1)
      call comm_summation(zbuf1,zbuf2,size(zbuf1),info_scf%icomm_r)
    end if
  end if

  ! send representive SCF rank to RT rank
  jrank = irank_src(io,ik)
  krank = irank_dst(io,ik)
  if (jrank == irank) then
    if (allocated(dbuf2)) call comm_send(dbuf2, krank, io*system_scf%nk+ik, icomm)
    if (allocated(zbuf2)) call comm_send(zbuf2, krank, io*system_scf%nk+ik, icomm)
  else if (krank == irank) then
    if (allocated(dbuf2)) call comm_recv(dbuf2, jrank, io*system_scf%nk+ik, icomm)
    if (allocated(zbuf2)) call comm_recv(zbuf2, jrank, io*system_scf%nk+ik, icomm)
  end if

  ! scatter rgrid data
  if (info_rt%ik_s <= ik .and. ik <= info_rt%ik_e .and. &
      info_rt%io_s <= io .and. io <= info_rt%io_e) then
    if (allocated(dbuf2)) then
      call comm_bcast(dbuf2, info_rt%icomm_r)
      shpsi%zwf(mg_rt%is(1):mg_rt%ie(1),mg_rt%is(2):mg_rt%ie(2),mg_rt%is(3):mg_rt%ie(3),:,io,ik,1) &
        = cmplx(dbuf2(mg_rt%is(1):mg_rt%ie(1),mg_rt%is(2):mg_rt%ie(2),mg_rt%is(3):mg_rt%ie(3),:))
    end if
    if (allocated(zbuf2)) then
      call comm_bcast(zbuf2, info_rt%icomm_r)
      shpsi%zwf(mg_rt%is(1):mg_rt%ie(1),mg_rt%is(2):mg_rt%ie(2),mg_rt%is(3):mg_rt%ie(3),:,io,ik,1) &
        = zbuf2(mg_rt%is(1):mg_rt%ie(1),mg_rt%is(2):mg_rt%ie(2),mg_rt%is(3):mg_rt%ie(3),:)
    end if
  end if
end do
end do

end subroutine convert_wave_function

end subroutine main_dft2tddft
