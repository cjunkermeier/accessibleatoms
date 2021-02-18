! ==================================================================
! AccessibleAtoms.f90
! Computes the number of atoms on which each molecular orbital resides.
! ==================================================================
! copyright info:
! @Copyright 2019 Chad Junkermeier
!
! GPL-3.0-or-later licensing:
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
! ==================================================================
! Code written by:
! Chad Junkermeier, PhD
! ===========================================================================
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
! ===========================================================================

! Program Declaration
! ===========================================================================

        program AccessibleAtoms


! How to call the program:
! ./AccessibleAtoms.x 95 156 48 C 4 N 4 H 1

! What each number means
! ./AccessibleAtoms.x  nkpoint neigenvec natom C nvalenceOrbitals N nvalenceOrbitals H nvalenceOrbitals


! Type Declaration
! ===========================================================================
       type T_atom
          character (len = 2) :: species             ! atomic species
          integer :: norbitals                       ! number of orbitals of species
       end type T_atom









! Variable Declaration and Description
! ==================================================================
!        implicit none
        integer :: a,i,j,k,l
        integer :: narguments
        integer :: kcheck
        integer :: eigveccheck

        character (len = 60) :: bah
        character (len = 60) :: natomchar
        character (len = 60) :: nkpointschar
        character (len = 60) :: neigenvecchar
        character (len = 2)  :: speciestempchar
        character (len = 1)  :: orbitalstempchar
        character (len = 2)  :: speciescheck
        integer :: natom
        integer :: nkpoints
        integer :: neigenvec
        integer :: nspecies
        integer :: orbitalstemp

        type(T_atom), dimension(:), allocatable :: atom_species

        real :: entropy, entropytemp, entropytemp2

        CHARACTER(LEN=18) :: keiFormatcheck
        CHARACTER(LEN=21) :: atomformat
        CHARACTER(LEN=10) :: orbitalformat
        CHARACTER(LEN=20) :: outputformat


! Parse command-line arguements
! ==================================================================
      IF(COMMAND_ARGUMENT_COUNT().lt.5)THEN
        WRITE(*,*)
        WRITE(*,*)'ERROR, TOO FEW ARGUEMENTS, STOPPING'
        WRITE(*,*)
        STOP
      ELSE IF (Mod( COMMAND_ARGUMENT_COUNT(), 2 ) .eq. 0) THEN
        WRITE(*,*)
        WRITE(*,*)'ERROR, THERE SHOULD BE AN ODD NUMBER OF ARGUEMENTS, STOPPING'
        WRITE(*,*)
        STOP
      END IF

      narguments = command_argument_count()

      !parse first three command-line values
      CALL GET_COMMAND_ARGUMENT(1,nkpointschar)
      CALL GET_COMMAND_ARGUMENT(2,neigenvecchar)
      CALL GET_COMMAND_ARGUMENT(3,natomchar)
      !convert first three command-line values to integers
      READ(nkpointschar,*)nkpoints
      READ(neigenvecchar,*)neigenvec
      READ(natomchar,*)natom

      nspecies = (narguments - 3)/2

      allocate(atom_species(nspecies))

      !loop through the rest of the command-line values
      j = 1
      do i=4, narguments, 2
          CALL GET_COMMAND_ARGUMENT(i, speciestempchar)
          CALL GET_COMMAND_ARGUMENT(i + 1, orbitalstempchar)

          READ(orbitalstempchar,*)orbitalstemp

          atom_species(j)%species = speciestempchar
          atom_species(j)%norbitals = orbitalstemp
          j = j + 1
     end do !i


! Input Output Format statements
! ===========================================================================

      keiFormatcheck = "(8X,I5.1,16X,I5.1)"
      atomformat = "(I5.1,1X,A2,40X,F8.1)"
      orbitalformat = "(48X,F8.1)"
      outputformat = "(I5.1,',',I5.1,',',F8.1)"


! Procedure
! ===========================================================================


! opening of eigenvec.out file and get to the first line of data.
      open(unit=13,file='eigenvec.out',status='old')
      read(13,*)
      read(13,*)
      j = 2 !data checking starts on the third line of eigenvec.out



! opening the output file accessible_atoms.csv
      open(unit=23,file='accessible_atoms.csv',status='unknown')
      write(23,*)"kpoint,eigenvec,W,IPR"




      kpointloop: DO k = 1, nkpoints
      eigenvecloop: DO i = 1, neigenvec

        read(13,keiFormatcheck) kcheck, eigveccheck             !Format check for kpoint and eigenvec
        read(13,*)                                              !Format check for kpoint and eigenvec
        j = j+2   !increase eigenvec.out line count

        if ((kcheck .ne. k) .and. (i .ne. eigveccheck)) THEN    !Format check for kpoint and eigenvec
          WRITE(*,*)
          WRITE(*,*)"ERROR! kcheck, eigveccheck don't match."   !Format check for kpoint and eigenvec
          WRITE(*,*)"On line:", j                               !Format check for kpoint and eigenvec
          WRITE(*,*)"STOPPING"
          WRITE(*,*)
          STOP
        END IF



      entropy = 0.0      !initialize the entropy-like value of accessible atoms for the kth kpoint and ith eigenvector.
      squaredsum = 0.0    !initialize the squared sum value of IRP for the kth kpoint and ith eigenvector.

      atomloop: DO a = 1, natom

        read(13,atomformat) natomcheck, speciescheck, entropytemp !parse atom number (as a check), species (to determine the number of orbitals on that atom), and the mulliken charge

        speciescheck = trim(speciescheck)

        if (entropytemp .ge. 1.e-6) THEN
           entropytemp = abs(entropytemp)
        end if




        if (a .ne. natomcheck) THEN
          WRITE(*,*)
          WRITE(*,*)"ERROR! natomcheck and atomloop don't match."
          WRITE(*,*)"On line:", j
          WRITE(*,*)"STOPPING"
          WRITE(*,*)
          STOP
        END IF


          speciesloop: DO l = 1, nspecies                            !tells the program how many orbitals the atom has
             if (speciescheck .eq. atom_species(l)%species) THEN     !tells the program how many orbitals the atom has
                 orbitalstemp = atom_species(l)%norbitals            !tells the program how many orbitals the atom has
             end if                                                  !tells the program how many orbitals the atom has
          end do speciesloop                                         !tells the program how many orbitals the atom has


          orbitalloop: Do l =2, orbitalstemp                         !loops over remaining orbitals for atoms with more than one orbital

            j = j + 1                                                !increment eigenvec.out line count
            read(13,orbitalformat)  entropytemp2                      !parse mulliken charge



        if (entropytemp2 .ge. 1.0e-6) THEN                                         !if nonzero compute the following.
           entropytemp = entropytemp + entropytemp2

        end if

          end do orbitalloop

          j = j + 2    !increment eigenvec.out line count
          read(13,*)   !parse a blank line



       if (entropytemp .ge. 1.0e-6) THEN                                         !if nonzero compute the following.
           entropy = entropy - abs(entropytemp) * log(abs(entropytemp)) !W
           squaredsum = squaredsum + abs(entropytemp) * abs(entropytemp) !IPR
        end if


      END DO atomloop

        write(23,'(I0,",",I0,",",F0.5,",",F0.5)') k,  i, exp(entropy), 1/squaredsum !write results to accessible_atoms.csv

      END DO eigenvecloop
      END DO kpointloop











! Closing files
      close(13)
      close(23)

! Deallocate
      deallocate(atom_species)


! end program
! ===========================================================================
      stop
      end program AccessibleAtoms




