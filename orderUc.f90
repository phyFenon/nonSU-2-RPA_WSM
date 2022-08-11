program FILEREADER
integer,parameter::n=48**3,m=4  !48 is the q-mesh in BZ.
double precision,dimension(n,m)::x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13
double precision,dimension(n,m)::x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27
double precision::t1,t2

open (unit=5+0, file='g0.txt', status='old', action='read')
open (unit=5+28, file='g0.5.txt', status='old', action='read')
open (unit=5+2, file='g0.6.txt', status='old', action='read')
open (unit=5+3, file='g1.0.txt', status='old', action='read')
open (unit=5+4, file='g1.2.txt', status='old', action='read')

! open (unit=5+5, file='g1.03.txt', status='old', action='read')
! open (unit=5+6, file='g1.04.txt', status='old', action='read')
! open (unit=5+7, file='g1.05.txt', status='old', action='read')
! g1.10-1.65  NO:12
! open (unit=5+8, file='g1.10.txt', status='old', action='read')
! open (unit=5+9, file='g1.15.txt', status='old', action='read')
! open (unit=5+10, file='g1.20.txt', status='old', action='read')
! open (unit=5+11, file='g1.25.txt', status='old', action='read')
! open (unit=5+12, file='g1.30.txt', status='old', action='read')
! open (unit=5+13, file='g1.35.txt', status='old', action='read')
! open (unit=5+14, file='g1.40.txt', status='old', action='read')
! open (unit=5+15, file='g1.45.txt', status='old', action='read')
! open (unit=5+16, file='g1.50.txt', status='old', action='read')
! open (unit=5+17, file='g1.55.txt', status='old', action='read')
! open (unit=5+18, file='g1.60.txt', status='old', action='read')
! open (unit=5+19, file='g1.65.txt', status='old', action='read')
! open (unit=5+20, file='g1.70.txt', status='old', action='read')
! open (unit=5+21, file='g1.80.txt', status='old', action='read')
! 
! open (unit=5+22, file='g2.0.txt', status='old', action='read')
! open (unit=5+23, file='g2.2.txt', status='old', action='read')
! open (unit=5+24, file='g2.4.txt', status='old', action='read')
! open (unit=5+25, file='g2.6.txt', status='old', action='read')
! open (unit=5+26, file='g2.8.txt', status='old', action='read')
! open (unit=5+27, file='g3.0.txt', status='old', action='read')


open (unit=35+0, file='g0_Uc.txt')
open (unit=35+28, file='g0.5_Uc.txt')
open (unit=35+2, file='g0.6_Uc.txt')

open (unit=35+3, file='g1.0_Uc.txt')
open (unit=35+4, file='g1.2_Uc.txt')

! open (unit=35+5, file='g1.03_Uc.txt')
! open (unit=35+6, file='g1.04_Uc.txt')
! open (unit=35+7, file='g1.05_Uc.txt')
! g1.10-1.65  NO:12
! open (unit=35+8, file='g1.10_Uc.txt')
! open (unit=35+9, file='g1.15_Uc.txt')
! open (unit=35+10, file='g1.20_Uc.txt')
! open (unit=35+11, file='g1.25_Uc.txt')
! open (unit=35+12, file='g1.30_Uc.txt')
! open (unit=35+13, file='g1.35_Uc.txt')
! open (unit=35+14, file='g1.40_Uc.txt')
! open (unit=35+15, file='g1.45_Uc.txt')
! open (unit=35+16, file='g1.50_Uc.txt')
! open (unit=35+17, file='g1.55_Uc.txt')
! open (unit=35+18, file='g1.60_Uc.txt')
! open (unit=35+19, file='g1.65_Uc.txt')
! 
! 
! open (unit=35+20, file='g1.70_Uc.txt')
! open (unit=35+21, file='g1.80_Uc.txt')
! open (unit=35+22, file='g2.0_Uc.txt')
! open (unit=35+23, file='g2.2_Uc.txt')
! open (unit=35+24, file='g2.4_Uc.txt')
! open (unit=35+25, file='g2.6_Uc.txt')
! open (unit=35+26, file='g2.8_Uc.txt')
! open (unit=35+27, file='g3.0_Uc.txt')

call cpu_time(t1)
do I=1,n,1
read(5+0,*) x0(I,:)
read(5+28,*) x1(I,:)
read(5+2,*) x2(I,:)
read(5+3,*) x3(I,:)
read(5+4,*) x4(I,:)
! read(5+5,*) x5(I,:)
! read(5+6,*) x6(I,:)
! read(5+7,*) x7(I,:)
! read(5+8,*) x8(I,:)
! read(5+9,*) x9(I,:)
! read(5+10,*) x10(I,:)
! read(5+11,*) x11(I,:)
! read(5+12,*) x12(I,:)
! read(5+13,*) x13(I,:)
! read(5+14,*) x14(I,:)
! read(5+15,*) x15(I,:)
! read(5+16,*) x16(I,:)
! read(5+17,*) x17(I,:)
! read(5+18,*) x18(I,:)
! read(5+19,*) x19(I,:)
! read(5+20,*) x20(I,:)
! read(5+21,*) x21(I,:)
! read(5+22,*) x22(I,:)
! read(5+23,*) x23(I,:)
! read(5+24,*) x24(I,:)
! read(5+25,*) x25(I,:)
! read(5+26,*) x26(I,:)
! read(5+27,*) x27(I,:)
enddo

!x0,x1 x2 x3 x4 x5 x6 x7 x8 x9 x10
call order_mintomax(x0,n,m)
call order_mintomax(x1,n,m)
call order_mintomax(x2,n,m)
call order_mintomax(x3,n,m)
call order_mintomax(x4,n,m)
! call order_mintomax(x5,n,m)
! call order_mintomax(x6,n,m)
! call order_mintomax(x7,n,m)
! call order_mintomax(x8,n,m)
! call order_mintomax(x9,n,m)
! call order_mintomax(x10,n,m)
! call order_mintomax(x11,n,m)
! call order_mintomax(x12,n,m)
! call order_mintomax(x13,n,m)
! call order_mintomax(x14,n,m)
! call order_mintomax(x15,n,m)
! call order_mintomax(x16,n,m)
! call order_mintomax(x17,n,m)
! call order_mintomax(x18,n,m)
! call order_mintomax(x19,n,m)
! call order_mintomax(x20,n,m)
! call order_mintomax(x21,n,m)
! call order_mintomax(x22,n,m)
! call order_mintomax(x23,n,m)
! call order_mintomax(x24,n,m)
! call order_mintomax(x25,n,m)
! call order_mintomax(x26,n,m)
! call order_mintomax(x27,n,m)



do I=1,n,1
write(35+0,*) x0(I,:)
write(35+28,*) x1(I,:)
write(35+2,*) x2(I,:)
write(35+3,*) x3(I,:)
write(35+4,*) x4(I,:)
! write(35+5,*) x5(I,:)
! write(35+6,*) x6(I,:)
! write(35+7,*) x7(I,:)
! write(35+8,*) x8(I,:)
! write(35+9,*) x9(I,:)
! write(35+10,*) x10(I,:)
! write(35+11,*) x11(I,:)
! write(35+12,*) x12(I,:)
! write(35+13,*) x13(I,:)
! write(35+14,*) x14(I,:)
! write(35+15,*) x15(I,:)
! write(35+16,*) x16(I,:)
! write(35+17,*) x17(I,:)
! write(35+18,*) x18(I,:)
! write(35+19,*) x19(I,:)
! write(35+20,*) x20(I,:)
! write(35+21,*) x21(I,:)
! write(35+22,*) x22(I,:)
! write(35+23,*) x23(I,:)
! write(35+24,*) x24(I,:)
! write(35+25,*) x25(I,:)
! write(35+26,*) x26(I,:)
! write(35+27,*) x27(I,:)

enddo

call cpu_time(t2)

print*,(t2-t1)/60.0,"min"
end

!!!----EIGENVALUE AND EIGENVECTOR IN a decending order(bigger to smaller)----!!!!
subroutine order_mintomax(a,row,col)
implicit none
integer::row,col,i,j
double precision,DIMENSION(row,col)::a
double precision,dimension(1:1,1:col)::temp
do i=1,row-1
do j=i+1,row
if (a(j,col) .lt. a(i,col))then
temp(1:1,1:col) = a(i:i,1:col)
a(i:i,1:col) = a(j:j,1:col)
a(j:j,1:col) = temp(1:1,1:col)

endif
enddo
enddo
end
