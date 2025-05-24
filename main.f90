! Copyright (C) 2015-2024 L. C. Ducati
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

program periodic
use param
implicit none
integer :: z, mm
!character (len=80) :: pinp
! conditional compilation block for OpenMP                                          
!$      integer :: nthreads, tid
!$      integer omp_get_num_threads,omp_get_thread_num                              
!$      external omp_get_num_threads,omp_get_thread_num                             
!$omp parallel private(tid,nthreads)                                                
!$      nthreads = omp_get_num_threads()                                            
!$      tid = omp_get_thread_num()                                                  
!$      if (tid == 1) then                                                          
!$         print '("PBC - Periodic Boundary Condition program")'
!$         print '(" Parallel execution with number of threads: ",i3)', nthreads              
!$      end if                                                                      
!$omp end parallel

if (iargc() < 1)then
  stop 'Please provide at least one cartesian coordinate input file'
end if
print '(i3," parsed input file(s)")', iargc()

do z=1, iargc()
  call getarg(z,arg)
  if(z >  1) print *,''
  print '(" Opening input file ",i3," -> ",a20)', z, arg

  open(unit=1,file='pbc.input',iostat=mm,status='old')
  if(mm /= 0)then
    linp=.True.
  else
    print *, 'External pbc.input file was found'
    call read_inp
    rewind(1)
    close(1,status='keep')
    linp=.False.
  end if  

  open(unit=1,file=arg,status='old')
  !
  !
  ! SUBROUTINE READ INP
  !
  if(linp) call read_inp
  !
  !
  if(pcart)then
    pinp2=arg(1:index(arg,'.')-1)//'_geompar.xyz'
    print *,'  New info file... ', pinp2
    open(unit=2,file=pinp2,status='replace')
  end if
  !
  !
  if(debug)then                                                                  
    pinp2=arg(1:index(arg,'.')-1)//'_debug.xyz'                                  
    print *,'  New info file... ', pinp2                                         
    open(unit=4,file=pinp2,status='replace')                                     
  end if  
  !
  ! SUBROUTINE READ CART
  !
  call read_cart()
  !
  !
  ! SUBROUTINE INCELL
  ! 
  if(in_cell) call incell()
  !
  !
  ! SUBROUTINE SCELL
  !
  if(full_cell) call scell()
  !
  !
  ! SUBROUTINE PRINT SCELL
  !
  if(pcell) call print_scell()
  !
  !
  ! SUBROUTINE PRINT CART
  !
  if(pcart) call print_cart()
  !
  !
  ! SUBROUTINE BONDS
  !
  if(cbonds) call bonds()
  !
  !
  !
  ! SUBROUTINE TOPOLOGY
  !
  if(topos) call topology()
  !
  !
  ! SUBROUTINE PRINT BONDS
  !
  if(pbonds) call print_bonds()
  !
  !
  ! SUBROUTINE PRINT ANGLE
  !
  if(pangs)call print_angle()
  !
  !
  ! SUBROUTINE PRINT DIHEDRAL
  !
    if(pdih) call print_dihedral()
  !
  !
  ! SUBROUTINE PRINT WITHIN
  !
    if(within) call print_within()
  !
  !
  ! SUBROUTINE MOLECULES
  !
  if(molecule) call molecules()
  !
  !
  ! SUBROUTINE CENTERMASS
  !
  if(cmass) call centermass()
  !
  !
  ! SUBROUTINE INCELL CENTER OF MASS
  !
  if(in_cellcmass) then
    ! DETOUR TO WRITE THE CMASS INFO AND RECLASSIFY THE MOLECULES IN THE NEW POSITION
    call incellcmass()
    ! REWIND 2 FOR NEW WRITING 
    rewind(2)
    ! PRINTING MOLECULES INSIDE THE CELL BASED ON CMASS CRITERIA 
    call print_cart()
    !
    if(allocated(bond)) deallocate(bond,stat=m)
    if(allocated(lbond)) deallocate(lbond,stat=m)
    call bonds()
    !
    if(allocated(molec)) deallocate(molec,stat=m)
    if(allocated(imolec)) deallocate(imolec,stat=m)
    if(allocated(molecmass)) deallocate(molecmass,stat=m)
    call molecules()
    !
    if(allocated(cmassx)) deallocate(cmassx,stat=m)
    if(allocated(cmassy)) deallocate(cmassy,stat=m)
    if(allocated(cmassz)) deallocate(cmassz,stat=m)
    call centermass()
  end if
  !
  !
  ! SUBROUTINE PRINT MOLECULES
  !
  if(pmolecule) call print_molecules()
  !
  !
  ! SUBROUTINE PRINT SOLUTE
  !
  if(psolute) call print_solute()
  !
  !
  ! SUBROUTINE NEIGHBOR
  !
  if(neighb) call neighbor()
  ! 
  !
  ! SUBROUTINE PRINT NEIGHBOR
  !
  if(pnns) call print_neighbor()
  !
  !
  ! SUBROUTINE PRINT TRAJECTORY
  !
  if(ptrajec) call print_trajectory()
  !
  !
  ! SUBROUTINE PRINT AXSF
  if(paxsf) call print_axsf()
  !
  !
  ! Closing units
  close(1,status='keep')
  if(pcart) close(2,status='keep')
!
!
!
  ! dealloccating dynamic memory
  if(allocated(rsymb)) deallocate(rsymb,stat=m)
  if(m /= 0) print *, "allocating rsymb failed"
  if(allocated(rxcoord)) deallocate(rxcoord,stat=m)
  if(m /= 0) print *, "allocating rxcoord failed"
  if(allocated(rycoord)) deallocate(rycoord,stat=m)
  if(m /= 0) print *, "allocating rycoord failed"
  if(allocated(rzcoord)) deallocate(rzcoord,stat=m)
  if(m /= 0) print *, "allocating rzcoord failed"
  if(allocated(rmap)) deallocate(rmap,stat=m)
  if(m /= 0) print *, "allocating rmap failed"
  if(allocated(symb)) deallocate(symb,stat=m)
  if(m /= 0) print *, "allocating symb failed"
  if(allocated(xcoord)) deallocate(xcoord,stat=m)
  if(m /= 0) print *, "allocating xcoord failed"
  if(allocated(ycoord)) deallocate(ycoord,stat=m)
  if(m /= 0) print *, "allocating ycoord failed"
  if(allocated(zcoord)) deallocate(zcoord,stat=m)
  if(m /= 0) print *, "allocating zcoord failed"
  if(allocated(amap)) deallocate(amap,stat=m)
  if(m /= 0) print *, "allocating amap failed"
  if(allocated(acharge)) deallocate(acharge,stat=m)
  if(m /= 0) print *, "allocating acharge failed"
  if(allocated(moleccharge)) deallocate(moleccharge,stat=m)
  if(m /= 0) print *, "allocating moleccharge failed"
  if(allocated(molec)) deallocate(molec,stat=m)
  if(m /= 0) print *, "allocating molec failed"
  if(allocated(imolec)) deallocate(imolec,stat=m)
  if(m /= 0) print *, "allocating imolec failed"
  if(allocated(bond)) deallocate(bond,stat=m)
  if(m /= 0) print *, "allocating bond failed"
  if(allocated(lbond)) deallocate(lbond,stat=m)
  if(m /= 0) print *, "allocating lbond failed"
  if(allocated(molec_bond)) deallocate(molec_bond,stat=m)
  if(m /= 0) print *, "allocating molec_bond failed"
  if(allocated(molec_map)) deallocate(molec_map,stat=m)
  if(m /= 0) print *, "allocating molec_map failed"
  if(allocated(molec_mask)) deallocate(molec_mask,stat=m)
  if(m /= 0) print *, "allocating molec_mask failed"
  if(allocated(sbond)) deallocate(sbond,stat=m)
  if(m /= 0) print *, "allocating sbond failed"
  if(allocated(sang)) deallocate(sang,stat=m)
  if(m /= 0) print *, "allocating sang failed"
  if(allocated(sdih)) deallocate(sdih,stat=m)
  if(m /= 0) print *, "allocating sdih failed"
  if(allocated(iamap)) deallocate(iamap,stat=m)
  if(m /= 0) print *, "allocating iamap failed"
  if(allocated(samap)) deallocate(samap,stat=m)
  if(m /= 0) print *, "allocating samap failed"
  if(allocated(iacharge)) deallocate(iacharge,stat=m)
  if(m /= 0) print *, "allocating iacharge failed"
  if(allocated(sacharge)) deallocate(sacharge,stat=m)
  if(m /= 0) print *, "allocating sacharge failed"
  if(allocated(molecmass)) deallocate(molecmass,stat=m)
  if(m /= 0) print *, "allocating molecmass failed"
  if(allocated(cmassx)) deallocate(cmassx,stat=m)
  if(m /= 0) print *, "allocating cmassx failed"
  if(allocated(cmassy)) deallocate(cmassy,stat=m)
  if(m /= 0) print *, "allocating cmassy failed"
  if(allocated(cmassz)) deallocate(cmassz,stat=m)
  if(m /= 0) print *, "allocating cmassz failed"
  if(allocated(itbody)) deallocate(itbody,stat=m)
  if(m /= 0) print *, "allocating itbody failed"
  if(allocated(ltbody)) deallocate(ltbody,stat=m)
  if(m /= 0) print *, "allocating ltbody failed"
end do
print '("Normal termination of PBC")'
end program
