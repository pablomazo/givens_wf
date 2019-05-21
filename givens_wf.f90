program givens
!program that calculates the WF for a given one-dimensional potential using the givens method

use diagonalization

implicit none

!ham: "hamiltonian" matrix
!wf: eigenvectors of the hamiltonian matrix
!eigval: eigenvalues of the hamiltonian matrix (not the energy!!!)
!npoints: number of points in which the wf is calculated
!intv: interval in which the wf is calculated
!step: space between two consecutive points
!rmass:reduced mass
!mass1,mass2: mass of the oscillators
!v: value of the potential in a given point of space
!norm: normalization factors for every functions

real(8), dimension(:,:), allocatable :: ham, wf
real(8), dimension(:), allocatable ::  eigval, norm
real(8),dimension(2) :: intv
real(8):: step, rmass, mass1, mass2, v
integer :: npoints,i ,j
!this vlaues will not be use
integer :: it_num, rot_num

!constants
real(8), parameter :: hbar=1.0545718d-34 !J*s

!---------------------------------------------------------
!definition of the number of points and distances
npoints=500
intv(1)=0  !m
intv(2)=3d-10 !m
step=(intv(2)-intv(1))/(npoints-1)

!allocation
allocate (ham(npoints,npoints))
allocate (wf(npoints,npoints))
allocate (eigval(npoints))
allocate (norm(npoints))

!all the values of matrices and vectors are set to 0
ham=0
wf=0
eigval=0

!calculation of the reduced mass
mass1= 1.66d-27    !kg   hidrógeno
mass2= 5.89d-26    !kg   cloro
rmass=(mass1*mass2)/(mass1+mass2) !kg

!the hamiltonian matrix is filled

do i=1,npoints-1
 ham(i,i+1)=1
 ham(i+1,i)=1
enddo

do i=1, npoints
! call v(intv(1)+step*(i-1))
 ham(i,i)=-2-2*(rmass*v(intv(1)+step*(i-1))*step**2)/hbar**2
enddo

!the ham matrix is diagonalize
call jacobi_eigenvalue ( npoints, ham, 30, wf, eigval, it_num, rot_num )

!once we have obtained the wf it must be normalise
norm=0

do i=1,npoints
 do j=1,npoints
  norm(i)=norm(i)+wf(j,i)**2
 enddo
 
 norm(i)=sqrt(norm(i)*step)

enddo

!the results are printed in a external file
open(10,file='wf.dat',status='replace')

write(10,*) '#  Position      Potential         Energy          Functions'

do i=1,npoints
!the crazy thing done in eigval and wf is just to order the orbital from the one with lowest energy to the one with highest
! Energy; Potential; Energy; Eigenvalues
 write(10,*) intv(1)+step*(i-1),v(intv(1)+step*(i-1)) , &
-(eigval(npoints-i+1)*hbar**2)/(2*rmass*step**2), (wf(i,npoints-j+1)/norm(npoints-j+1), j=1,npoints)
enddo
endprogram 
!---------------------------------------------------------
