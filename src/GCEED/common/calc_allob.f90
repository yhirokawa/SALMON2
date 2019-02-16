!
!  Copyright 2017 SALMON developers
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
subroutine calc_allob(iob,iob_allob)
  use salmon_parallel, only: nproc_id_kgrid
  use scf_data
  use new_world_sub
  implicit none
  integer,intent(in) :: iob
  integer,intent(out) ::  iob_allob
  integer :: iob_tmp
  
  if(ilsda==0)then
    if(iparaway_ob==1)then
      iob_allob=nproc_id_kgrid*itotMST/nproc_ob+iob
    else if(iparaway_ob==2)then
      iob_allob=(iob-1)*nproc_ob+nproc_id_kgrid+1
    end if
  else
    if(iparaway_ob==1)then
      if(iob<=iobnum/2)then
        iob_allob=nproc_id_kgrid*MST(1)/nproc_ob+iob
      else
        iob_tmp=iob-iobnum/2
        iob_allob=nproc_id_kgrid*MST(2)/nproc_ob+iob_tmp+MST(1)
      end if
    else if(iparaway_ob==2)then
      if(iob<=iobnum/2)then
        iob_allob=(iob-1)*nproc_ob+nproc_id_kgrid+1
      else
        iob_tmp=iob-iobnum/2
        iob_allob=(iob_tmp-1)*nproc_ob+nproc_id_kgrid+1+MST(1)
      end if
    end if
  end if
  
end subroutine calc_allob
