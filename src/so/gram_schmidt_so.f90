!
!  Copyright 2019-2020 SALMON developers
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
!=======================================================================

module gram_schmidt_so_sub

  use structures, only: s_dft_system, s_rgrid, s_parallel_info, s_orbital
  use pack_unpack, only: copy_data
  use timer
  use communication, only: comm_bcast, comm_summation
  use spin_orbit_global, only: SPIN_ORBIT_ON

  implicit none

  private
  public :: gram_schmidt_so
  public :: SPIN_ORBIT_ON

  complex(8),allocatable :: wf_io1(:,:,:,:)
  complex(8),allocatable :: wf_exc(:,:,:,:)
  complex(8),allocatable :: wf_exc_tmp(:,:,:,:)
  complex(8),allocatable :: coeff_tmp(:), coeff(:)

  complex(8),parameter :: zero=(0.0d0,0.0d0)
  complex(8),parameter :: one =(1.0d0,0.0d0)

contains

  subroutine gram_schmidt_so( sys, rg, wfi, wf )
    implicit none
    type(s_dft_system),       intent(in)    :: sys
    type(s_rgrid),            intent(in)    :: rg
    type(s_parallel_info)   , intent(in)    :: wfi
    type(s_orbital),          intent(inout) :: wf
    integer :: im, ik, is, io1, io2
    real(8) :: norm2, norm2_tmp

    call timer_begin(LOG_CALC_GRAM_SCHMIDT)

    if ( .not.allocated(wf_io1) ) then
      allocate( wf_io1(size(wf%zwf,1),size(wf%zwf,2),size(wf%zwf,3),size(wf%zwf,4)) )
      allocate( wf_exc(size(wf%zwf,1),size(wf%zwf,2),size(wf%zwf,3),size(wf%zwf,4)) )
      allocate( wf_exc_tmp(size(wf%zwf,1),size(wf%zwf,2),size(wf%zwf,3),size(wf%zwf,4)) )
      wf_io1=zero
      wf_exc=zero
      wf_exc_tmp=zero
      allocate( coeff_tmp(sys%no-1) ); coeff_tmp=zero
      allocate( coeff(sys%no-1) ); coeff=zero
    end if

    do im = wfi%im_s, wfi%im_e
    do ik = wfi%ik_s, wfi%ik_e

      do io1 = 1, sys%no

        ! Retrieve orbit #io1 into "wf_io1":
        if ( wfi%io_s <= io1 .and. io1 <= wfi%io_e ) then
          call copy_data( wf%zwf(:,:,:,:,io1,ik,im), wf_io1 )
        end if
        if ( wfi%if_divide_orbit ) then
          ! comm_bcast for double_complex_4d has not been implemented yet
          do is=1,size(wf_io1,4)
            call comm_bcast( wf_io1(:,:,:,is), wfi%icomm_o, wfi%irank_io(io1) )
          end do
        end if

        ! Calculate overlap coefficients: 

        coeff_tmp=zero

        do io2 = 1, io1 - 1
          if ( wfi%io_s <= io2 .and. io2 <= wfi%io_e ) then
            coeff_tmp(io2) = dot_wf( wf%zwf(:,:,:,:,io2,ik,im), wf_io1 )
          end if
        end do 

        if ( wfi%if_divide_rspace .or. wfi%if_divide_orbit ) then
          call comm_summation( coeff_tmp, coeff, size(coeff), wfi%icomm_ro )
        else
          coeff = coeff_tmp
        end if

        ! Calculate exclusion term "wf_exc":
        wf_exc_tmp = 0.0d0
        do io2 = 1, io1 - 1
          if ( wfi%io_s <= io2 .and. io2 <= wfi%io_e ) then
            call axpy_wf( coeff(io2), wf%zwf(:,:,:,:,io2,ik,im), wf_exc_tmp )
          end if
        end do
        if ( wfi%if_divide_orbit ) then
          call comm_summation( wf_exc_tmp, wf_exc, size(wf_exc), wfi%icomm_o )
        else
          call copy_data( wf_exc_tmp, wf_exc )
        end if

        if ( wfi%io_s <= io1 .and. io1 <= wfi%io_e ) then
          ! Exclude non-orthonormal component:
          call axpy_wf( -one, wf_exc, wf_io1 )
          ! Normalization:
          norm2_tmp = real( dot_wf(wf_io1,wf_io1) )
          if ( wfi%if_divide_rspace ) then
            call comm_summation( norm2_tmp, norm2, wfi%icomm_r )
          else
            norm2 = norm2_tmp
          end if
          call scal_wf( one/sqrt(norm2), wf_io1 )
          ! Write back to "zwf":
          call copy_data( wf_io1, wf%zwf(:,:,:,:,io1,ik,im) ) 
        end if
      end do !io1

    end do !ik
    end do !im
    
    call timer_end(LOG_CALC_GRAM_SCHMIDT)

    return

  contains

  ! Dot product of two wavefunctions:
  complex(8) function dot_wf(x, y) result(p)
    implicit none
    complex(8), intent(in) :: x(:,:,:,:)
    complex(8), intent(in) :: y(:,:,:,:)
    integer :: i1, i2, i3, i4
    p=zero
    !$omp parallel do collapse(2) default(shared) private(i1,i2,i3) reduction(+:p)
    do i4 = 1, size(x,4)
      do i3 = 1, size(x,3)
        do i2 = 1, size(x,2)
          do i1 = 1, size(x,1)
            p = p + conjg( x(i1,i2,i3,i4) )*y(i1,i2,i3,i4)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    p = p * sys%hvol
    return
  end function dot_wf

  ! Constant times a wavefunction plus a wavefunction:
  subroutine axpy_wf(a, x, y)
    implicit none
    complex(8), intent(in) :: a
    complex(8), intent(in) :: x(:,:,:,:)
    complex(8), intent(inout) :: y(:,:,:,:)
    integer :: i1, i2, i3, i4
    !$omp parallel do collapse(2) default(shared) private(i1,i2,i3)
    do i4 = 1, size(x,4)
      do i3 = 1, size(x,3)
        do i2 = 1, size(x,2)
          do i1 = 1, size(x,1)
            y(i1, i2, i3, i4) = a * x(i1, i2, i3, i4) + y(i1, i2, i3, i4)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    return
  end subroutine axpy_wf


  ! Scales a wavefunction by a constant 
  subroutine scal_wf(a, x)
    implicit none
    complex(8), intent(in) :: a
    complex(8), intent(inout) :: x(:,:,:,:)
    integer :: i1, i2, i3, i4
    !$omp parallel do collapse(2) default(shared) private(i1,i2,i3)
    do i4 = 1, size(x,4)
      do i3 = 1, size(x,3)
        do i2 = 1, size(x,2)
          do i1 = 1, size(x,1)
            x(i1, i2, i3, i4) = a * x(i1, i2, i3, i4)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    return
  end subroutine scal_wf

  end subroutine gram_schmidt_so

end module gram_schmidt_so_sub
