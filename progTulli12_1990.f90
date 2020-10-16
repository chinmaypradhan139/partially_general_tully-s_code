program tulli12
use mod_global
implicit none
integer :: i
integer :: res,r1,ser,r0,per,l1,rep,l0,ip,ks,mp,jk
res=0
r1=0
ser=0
r0=0
per=0
l1=0
rep=0
l0=0
i=1
mp=1
call setup_initial_values1
do mp=1,30
ks=mp
r1=0;r0=0;l0=0;l1=0
do i=1,1000
call setup_initial_values2
call evolve(ip,res,r1,ser,r0,per,l1,rep,l0,ks)
!write(*,*) i,pos,lambda
end do
write(*,*) ip,l0,'l0',r0,'r0',l1,'l1',r1,'r1'
write(1,*) mp,'mp',l0/5,'l0',r0/5,'r0',l1/5,'l1',r1/5,'r1'

enddo

end program
!..............................................................................


