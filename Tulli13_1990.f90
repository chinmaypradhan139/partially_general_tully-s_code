module mod_global
implicit none
  integer :: nsize,INFO,total_dimensions,lambda,Hi 
  real*8 :: mass
  real*8 :: dt,time,rnd
  complex*16,dimension(:), allocatable :: c
  real*8,dimension(:,:), allocatable :: H,Energy_hamil,acc1,g
  real*8,dimension(:,:,:), allocatable :: Gradient_hamil,acw
  real*8,dimension(:), allocatable :: Energy
  real*8,dimension(:), allocatable :: pos,v,old_acc2
  
contains
!.......................................................................
subroutine setup_initial_values1
implicit none
Hi=2 !Hamiltonian dimensions
total_dimensions=1
mass=2000   
  allocate(pos(total_dimensions))
  allocate(v(total_dimensions))
  allocate(acc1(Hi,total_dimensions))
  allocate(old_acc2(total_dimensions))
  allocate(acw(Hi,Hi,total_dimensions))
  allocate(H(Hi,Hi))
  allocate(Energy_hamil(Hi,Hi)) 
  allocate(Energy(Hi))
  allocate(Gradient_hamil(Hi,Hi,total_dimensions))
  allocate(c(Hi))
  allocate(g(Hi,Hi))
end subroutine
!.........................................................

subroutine setup_initial_values2
dt=1
v(1)=(1/mass)  
pos(1)=-4.0
time=0.00
c(1)=1.0
c(2)=0.0
lambda=1
end subroutine 
!........................................................................
subroutine diag_wrapper(matrix,nsize,eigen_values,eigen_vectors)
real*8, intent(inout) :: matrix(nsize,nsize)
real*8, dimension(nsize,nsize) :: mat
integer LWORK,nsize
real*8, allocatable :: WORK(:)
real*8, intent(out) :: eigen_vectors(nsize,nsize),eigen_values(nsize)
mat=matrix
LWORK=3*nsize-1
allocate(WORK(LWORK))
call dsyev('V','U',nsize,mat,nsize,eigen_values,WORK,LWORK,INFO)
eigen_vectors=mat

end subroutine

!..........................................................................
subroutine potential
implicit none
real*8 :: V11,V22,V12,gV11,gV12,gV22
real*8 ::A1=0.01,B1=1.6,C1=0.005,D1=1.0


 if (pos(1)>0) then
   V11=A1*(1-exp(-B1*pos(1)))
 else 
   V11=-A1*(1-exp(B1*pos(1)))
 endif
 V22=-V11

 V12=C1*exp(-D1*(pos(1)**2))



nsize=2
H(1,1)=V11
H(1,2)=V12
H(2,1)=V12
H(2,2)=V22
call diag_wrapper(H,nsize,Energy,Energy_hamil)

if (pos(1)>0) then
  gV11=A1*B1*exp(-B1*pos(1))
 else
  gV11=A1*B1*exp(B1*pos(1))
 endif
 gV12=-C1*D1*exp(-D1*pos(1)*pos(1))*2*pos(1)
 gV22=-gV11
Gradient_hamil(1,1,1)=gV11
Gradient_hamil(1,2,1)=gV12
Gradient_hamil(2,1,1)=gV12
Gradient_hamil(2,2,1)=gV22


end subroutine

!...........................................................................
subroutine nonadiabaticvector
implicit none
integer :: acwi,i,j
do acwi=1,total_dimensions
do i=1,Hi
do j=1,Hi
if (i.ne.j) then
acw(i,j,acwi)=(sum(Energy_hamil(:,i)*matmul(Gradient_hamil(:,:,acwi),Energy_hamil(:,j))))/(Energy(j)-Energy(i))
else
acw(i,j,acwi)=0
endif
end do
end do
end do
end subroutine
!..........................................................................
subroutine force
implicit none
integer :: acc1i,p
do acc1i=1,total_dimensions
do p=1,Hi
acc1(p,acc1i)=(-sum(Energy_hamil(:,p)*matmul(Gradient_hamil(:,:,acc1i),Energy_hamil(:,p))))/mass
end do
end do
end subroutine
!..........................................................................
subroutine Rungekutta
implicit none
integer :: i,j,p,q
complex*16,dimension(Hi,Hi) :: Vad
complex*16,dimension(Hi) ::k1,k2,k3,k4
complex*16,dimension(Hi,Hi) ::M

do i=1,Hi
do j=1,Hi 
if (i.eq.j) then
Vad(i,j)=Energy(i)/(cmplx(0,1))
else
Vad(i,j)=0
end if
end do
end do

do p=1,Hi
do q=1,Hi
M(p,q)=Vad(p,q)-sum(v*acw(p,q,:))
enddo
enddo

k1=dt*matmul(M,c)
k2=dt*matmul(M,(c+k1/2))
k3=dt*matmul(M,(c+k2/2))
k4=dt*matmul(M,(c+k3))
c=c+(k1+2*k2+2*k3+k4)/6
end subroutine
!.........................................................................
subroutine gs
implicit none
integer :: i,j,l,m
real*8, dimension(Hi,Hi) :: b
do i=1,Hi
do j=1,Hi
if (i.ne.j) then
b(i,j)=-2*real(conjg(c(i))*c(j)*sum(v*acw(i,j,:)))
else
b(i,j)=0
end if
end do
end do


do l=1,Hi
do m=1,Hi
g(l,m)=(dt*b(m,l))/(conjg(c(l))*(c(l)))
end do
end do
end subroutine
!...........................................................................
subroutine evolve(ip,res,r1,ser,r0,per,l1,rep,l0,ks)
implicit none
integer :: r,s,t
integer, intent(inout) :: ip,res,r1,ser,r0,per,l1,rep,l0,ks
real*8 :: para_v          !para_v= component of v parallel to acw
real*8, dimension(:),allocatable :: perp_v 
real*8, dimension(Hi,Hi) :: hp
allocate(perp_v(total_dimensions))
ip=ks                     
v=v*ip
 
 call potential
 call nonadiabaticvector
 call Rungekutta
 call random_number(rnd)
do while(abs(pos(1))<6.00)
!time=time+dt
call force
call random_number(rnd)
call gs

do r=1,Hi
hp(lambda,r)=g(lambda,r)
!hp(lambda,r)=sum(g(lambda,1:r))
s=r-1
if (s>0) then
hp(lambda,r)=hp(lambda,r)+g(lambda,s)
s=s-1
end if
enddo


!...........................................................................e.

do t=1,Hi
if (t.ne.lambda) then
if (hp(lambda,t)>rnd.and.((sum(v*acw(lambda,t,:))/norm2(acw(lambda,t,:)))**2+(2*Energy(lambda)/mass)-(2*Energy(t)/mass))>0) then
    para_v=sum(v*acw(lambda,t,:))/norm2(acw(lambda,t,:))
    perp_v=v-para_v*acw(lambda,t,:)/norm2(acw(lambda,t,:))
    para_v=((para_v)/abs(para_v))*sqrt(para_v**2+(2*Energy(lambda)/mass)-(2*Energy(t)/mass))
    v=(para_v)*acw(lambda,t,:)/norm2(acw(lambda,t,:))+perp_v
    lambda=t
end if
end if
end do
!................................................................................

pos=pos+v*dt+0.5*acc1(lambda,:)*dt*dt
old_acc2=acc1(lambda,:)
call potential
call force
v=v+0.5*(acc1(lambda,:)+old_acc2)*dt
call nonadiabaticvector
call Rungekutta
call gs
!if (lambda.eq.1) then
!write(4,*) pos(1),'pos',0.5*mass*v**2+Energy(1),'TE1',acw(1,2,1),'acw',lambda, & 
!'lambda',c(1),'c1',c(2),'c2',conjg(c(1))*c(1)+conjg(c(2))*c(2),'c(1)^2+c(2)^2',hp(1,2),'hp12'
!else
!write(4,*)pos(1),'pos',0.5*mass*v**2+Energy(2),'TE2',acw(1,2,1),'acw',lambda, & 
!'lambda',c(1),'c(1)',c(2),'c(2)',conjg(c(1))*c(1)+conjg(c(2))*c(2),'c1^2+c2^2',hp(1,2),'hp12'
!end if
enddo
call traj(res,r1,ser,r0,per,l1,rep,l0,ip,ks)
!endif
!ip=ip+1


end subroutine
!................................................................................
subroutine traj(res,r1,ser,r0,per,l1,rep,l0,ip,ks)
implicit none
integer, intent(inout) :: res,r1,ser,r0,per,l1,rep,l0,ip,ks
if (pos(1)>0.and.(lambda.eq.2)) then
res=res+1
endif
 do while(res>0)
    res=res-1
    r1=r1+1
 end do
if (pos(1)>0.and.(lambda.eq.1)) then
ser=ser+1
endif
 do while(ser>0)
    ser=ser-1
    r0=r0+1
enddo
if (pos(1)<0.and.(lambda.eq.2)) then
per=per+1
endif 
 do while(per>0)
    per=per-1
    l1=l1+1
 enddo
if (pos(1)<0.and.(lambda.eq.1)) then
rep=rep+1
endif
 do while(rep>0)
    rep=rep-1
    l0=l0+1
 enddo
end subroutine
end module
!....................................................................................

