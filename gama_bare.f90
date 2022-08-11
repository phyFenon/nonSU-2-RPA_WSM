!!!this is for weyl semimetals system and calculate its cdw and sdw
Program wssusceptability
use omp_lib
Implicit None
Double Precision, Parameter :: pi = atan(1.d0)*4
integer,parameter::steps=24
Integer::i,j,k,N,qnx,qny,qnz
Double precision::kx,ky,kz,hubbardu
Double precision::beta,wn !1/Temperature,the frequency
Complex *16,Dimension(2,2)::uq,u ! the eigenvectors correspond to k+q and k
Double precision,Dimension(2)::eigenq,eigen   !the eigenvalues of k+q and k
double precision:: delta,gama,eps,t1,t2
double precision,dimension(steps)::qx,qy,qz
Complex *16,Dimension(4,4)::kai,kai0,kaivr,kairpa! suscepatability
Complex *16,Dimension(4)::kaiw_temp! for the single k
Complex *16,Dimension(steps,steps,steps,4)::kaiw! the eigenvalue
double precision,dimension(4,4)::umatrix,imatrix
Complex*16::E1,E2,E3,E4,kaidd,kaiss,kaipm,kaimp




 gama=0.1


 open(10,file='g0.1.txt')


 !---the routine of Q--------!
 !    ky     *
 !       *   *
 !    *======* kx
 !---------------------------!
 wn=0.0d0 !the transvese frequency
 delta=0.01d0 ! to avoid the singularity in numerator wn+i delta
 beta=100.d0 !temperature



imatrix=0.d0
do i=1,4
imatrix(i,i)=1.d0
enddo




umatrix=0.d0
umatrix(1,4)=1.d0
umatrix(2,3)=-1.d0
umatrix(3,2)=umatrix(2,3)
umatrix(4,1)=umatrix(1,4)

t1=omp_get_wtime()

!$omp parallel do private(kx,ky,kz,E1,E2,E3,E4,u,uq,eigen,eigenq,kai,kai0,kairpa,kaiw_temp)
do qnz=1,steps,1
qz(qnz) = -pi+(qnz*2.d0/steps)*pi
print*,qnz

do qnx=1,steps,1
qx(qnx) = -pi+(qnx*2.d0/steps)*pi


DO qny=1,steps,1

qy(qny) = -pi+(qny*2.d0/steps)*pi

!envalue the kai!
 kai=0.d0

 DO i=1,steps,1
 kx = (- 1.d0 + i * 2.d0 / steps)*pi
 Do j = 1, steps, 1
 ky= (- 1.d0 + j * 2.d0 / steps)*pi
 Do k = 1, steps, 1
 kz= (- 1.d0 + k * 2.d0 / steps)*pi


call eigensystems(kx,ky,kz,gama,u,eigen)
call eigensystems(kx+qx(qnx),ky+qy(qny),kz+qz(qnz),gama,uq,eigenq)

 E1=1.d0/(exp(beta*eigen(1))+1.d0)-1.d0/(exp(beta*eigenq(1))+1.d0)
 E1=-E1/cmplx(eigen(1)-eigenq(1)+wn,delta)


 E2=1.d0/(exp(beta*eigen(1))+1.d0)-1.d0/(exp(beta*eigenq(2))+1.d0)
 E2=-E2/cmplx(eigen(1)-eigenq(2)+wn,delta)

 E3=1.d0/(exp(beta*eigen(2))+1.d0)-1.d0/(exp(beta*eigenq(1))+1.d0)
 E3=-E3/cmplx(eigen(2)-eigenq(1)+wn,delta)

 E4=1.d0/(exp(beta*eigen(2))+1.d0)-1.d0/(exp(beta*eigenq(2))+1.d0)
 E4=-E4/cmplx(eigen(2)-eigenq(2)+wn,delta)


!uuuu
kai(1,1)=kai(1,1)+E1*u(1,1)*conjg(u(1,1))*uq(1,1)*conjg(uq(1,1))&
                &+E2*u(1,1)*conjg(u(1,1))*uq(1,2)*conjg(uq(1,2))&
                &+E3*u(1,2)*conjg(u(1,2))*uq(1,1)*conjg(uq(1,1))&
                &+E4*u(1,2)*conjg(u(1,2))*uq(1,2)*conjg(uq(1,2))
!uuud!
kai(1,2)=kai(1,2)+E1*u(1,1)*conjg(u(1,1))*uq(1,1)*conjg(uq(2,1))&
                &+E2*u(1,1)*conjg(u(1,1))*uq(1,2)*conjg(uq(2,2))&
                &+E3*u(1,2)*conjg(u(1,2))*uq(1,1)*conjg(uq(2,1))&
                &+E4*u(1,2)*conjg(u(1,2))*uq(1,2)*conjg(uq(2,2))
!duuu
kai(1,3)=kai(1,3)+E1*u(2,1)*conjg(u(1,1))*uq(1,1)*conjg(uq(1,1))&
                &+E2*u(2,1)*conjg(u(1,1))*uq(1,2)*conjg(uq(1,2))&
                &+E3*u(2,2)*conjg(u(1,2))*uq(1,1)*conjg(uq(1,1))&
                &+E4*u(2,2)*conjg(u(1,2))*uq(1,2)*conjg(uq(1,2))

!duud
kai(1,4)=kai(1,4)+E1*u(2,1)*conjg(u(1,1))*uq(1,1)*conjg(uq(2,1))&
                &+E2*u(2,1)*conjg(u(1,1))*uq(1,2)*conjg(uq(2,2))&
                &+E3*u(2,2)*conjg(u(1,2))*uq(1,1)*conjg(uq(2,1))&
                &+E4*u(2,2)*conjg(u(1,2))*uq(1,2)*conjg(uq(2,2))
!uduu
kai(2,1)=kai(2,1)+E1*u(1,1)*conjg(u(2,1))*uq(1,1)*conjg(uq(1,1))&
                &+E2*u(1,1)*conjg(u(2,1))*uq(1,2)*conjg(uq(1,2))&
                &+E3*u(1,2)*conjg(u(2,2))*uq(1,1)*conjg(uq(1,1))&
                &+E4*u(1,2)*conjg(u(2,2))*uq(1,2)*conjg(uq(1,2))
!udud
kai(2,2)=kai(2,2)+E1*u(1,1)*conjg(u(2,1))*uq(1,1)*conjg(uq(2,1))&
                &+E2*u(1,1)*conjg(u(2,1))*uq(1,2)*conjg(uq(2,2))&
                &+E3*u(1,2)*conjg(u(2,2))*uq(1,1)*conjg(uq(2,1))&
                &+E4*u(1,2)*conjg(u(2,2))*uq(1,2)*conjg(uq(2,2))
!dduu
kai(2,3)=kai(2,3)+E1*u(2,1)*conjg(u(2,1))*uq(1,1)*conjg(uq(1,1))&
                &+E2*u(2,1)*conjg(u(2,1))*uq(1,2)*conjg(uq(1,2))&
                &+E3*u(2,2)*conjg(u(2,2))*uq(1,1)*conjg(uq(1,1))&
                &+E4*u(2,2)*conjg(u(2,2))*uq(1,2)*conjg(uq(1,2))
!ddud
kai(2,4)=kai(2,4)+E1*u(2,1)*conjg(u(2,1))*uq(1,1)*conjg(uq(2,1))&
                &+E2*u(2,1)*conjg(u(2,1))*uq(1,2)*conjg(uq(2,2))&
                &+E3*u(2,2)*conjg(u(2,2))*uq(1,1)*conjg(uq(2,1))&
                &+E4*u(2,2)*conjg(u(2,2))*uq(1,2)*conjg(uq(2,2))
!uudu
kai(3,1)=kai(3,1)+E1*u(1,1)*conjg(u(1,1))*uq(2,1)*conjg(uq(1,1))&
                &+E2*u(1,1)*conjg(u(1,1))*uq(2,2)*conjg(uq(1,2))&
                &+E3*u(1,2)*conjg(u(1,2))*uq(2,1)*conjg(uq(1,1))&
                &+E4*u(1,2)*conjg(u(1,2))*uq(2,2)*conjg(uq(1,2))
!uudd
kai(3,2)=kai(3,2)+E1*u(1,1)*conjg(u(1,1))*uq(2,1)*conjg(uq(2,1))&
                &+E2*u(1,1)*conjg(u(1,1))*uq(2,2)*conjg(uq(2,2))&
                &+E3*u(1,2)*conjg(u(1,2))*uq(2,1)*conjg(uq(2,1))&
                &+E4*u(1,2)*conjg(u(1,2))*uq(2,2)*conjg(uq(2,2))
!dudu
kai(3,3)=kai(3,3)+E1*u(2,1)*conjg(u(1,1))*uq(2,1)*conjg(uq(1,1))&
                &+E2*u(2,1)*conjg(u(1,1))*uq(2,2)*conjg(uq(1,2))&
                &+E3*u(2,2)*conjg(u(1,2))*uq(2,1)*conjg(uq(1,1))&
                &+E4*u(2,2)*conjg(u(1,2))*uq(2,2)*conjg(uq(1,2))
!dudd
kai(3,4)=kai(3,4)+E1*u(2,1)*conjg(u(1,1))*uq(2,1)*conjg(uq(2,1))&
                &+E2*u(2,1)*conjg(u(1,1))*uq(2,2)*conjg(uq(2,2))&
                &+E3*u(2,2)*conjg(u(1,2))*uq(2,1)*conjg(uq(2,1))&
                &+E4*u(2,2)*conjg(u(1,2))*uq(2,2)*conjg(uq(2,2))
!uddu
kai(4,1)=kai(4,1)+E1*u(1,1)*conjg(u(2,1))*uq(2,1)*conjg(uq(1,1))&
                &+E2*u(1,1)*conjg(u(2,1))*uq(2,2)*conjg(uq(1,2))&
                &+E3*u(1,2)*conjg(u(2,2))*uq(2,1)*conjg(uq(1,1))&
                &+E4*u(1,2)*conjg(u(2,2))*uq(2,2)*conjg(uq(1,2))
!uddd
kai(4,2)=kai(4,2)+E1*u(1,1)*conjg(u(2,1))*uq(2,1)*conjg(uq(2,1))&
                &+E2*u(1,1)*conjg(u(2,1))*uq(2,2)*conjg(uq(2,2))&
                &+E3*u(1,2)*conjg(u(2,2))*uq(2,1)*conjg(uq(2,1))&
                &+E4*u(1,2)*conjg(u(2,2))*uq(2,2)*conjg(uq(2,2))
!dddu
kai(4,3)=kai(4,3)+E1*u(2,1)*conjg(u(2,1))*uq(2,1)*conjg(uq(1,1))&
                &+E2*u(2,1)*conjg(u(2,1))*uq(2,2)*conjg(uq(1,2))&
                &+E3*u(2,2)*conjg(u(2,2))*uq(2,1)*conjg(uq(1,1))&
                &+E4*u(2,2)*conjg(u(2,2))*uq(2,2)*conjg(uq(1,2))
!dddd
kai(4,4)=kai(4,4)+E1*u(2,1)*conjg(u(2,1))*uq(2,1)*conjg(uq(2,1))&
                &+E2*u(2,1)*conjg(u(2,1))*uq(2,2)*conjg(uq(2,2))&
                &+E3*u(2,2)*conjg(u(2,2))*uq(2,1)*conjg(uq(2,1))&
                &+E4*u(2,2)*conjg(u(2,2))*uq(2,2)*conjg(uq(2,2))
enddo
enddo
enddo

kai=kai/(steps**3)

!Diagnose the bare kai
kai0=kai
kai0=matmul(kai0,umatrix)
kairpa=kai

call T_diagnose(kai0,4,kaiw_temp,kaivr)


kaiw(qnz,qnx,qny,1:4)=kaiw_temp


enddo
enddo
enddo
!$omp end parallel do

do i=1,steps
do j=1,steps
do k=1,steps

write(10,*)qz(i),qx(j),qy(k),-1.0/real(kaiw(i,j,k,1)),real(kaiw(i,j,k,1:4))

enddo
enddo
enddo

t2=omp_get_wtime()
print*,'total time:',(t2-t1)/60.0,"min"




End Program

!!!!get the eigenvalues and eigenvectors of Hamiltonian.
Subroutine eigensystems(kx,ky,kz,gama,H,W)
      Implicit None
      Double Precision, Parameter :: pi = atan(1.d0)*4
      Integer, Parameter :: LWORK = 5000 ,n=2!LWORK>=MAX(2*2*N-1)
      Integer :: INFO
      Complex * 16, Dimension (n, n) :: H
      Double Precision, Dimension (n) :: W
      Complex * 16, Dimension (LWORK) :: WORK
      Double Precision, Dimension (3*100-2) :: RWORK
      Complex * 16, Dimension(2,2) :: sigmax,sigmay,sigmaz,sigma0
      double precision::kx0,kx,ky,kz,m,tx,gama,gamax ! unit of hopping t

      m=2.d0
      gamax=0.5d0
      tx=1.d0
      kx0=pi/2.d0

      sigma0=0.d0
      sigmax=0.d0
      sigmay=0.d0
      sigmaz=0.d0 
      
      sigma0(1,1)=1.d0
      sigma0(2,2)=1.d0

      sigmax(1,2)=1.d0
      sigmax(2,1)=1.d0

      sigmaz(1,1)=1.d0
      sigmaz(2,2)=-1.d0
	 
      sigmay(1,2)=cmplx(0.d0,-1.d0)
      sigmay(2,1)=cmplx(0.d0,1.d0)

      H = 0.d0
!---------the whole H get valued-----!
      H=gama*(cos(kx)-cos(kx0))*sigma0-gamax*(cos(3*kx)-cos(3*kx0))*sigmax-m*(2.d0-cos(ky)-cos(kz))*sigmax&
       &-2.d0*(cos(kx)-cos(kx0))*sigmax-2.d0*sin(ky)*sigmay-2.d0*sin(kz)*sigmaz
!---------the end---------------------------!

Call zheev ('V', 'U', n, H, n, W, WORK, LWORK, RWORK, INFO)

End


!!-----TO USE THE LAPCAK DIAGNOSE THE NONSYMMETRIC MATRIX!!!!
Subroutine T_diagnose (TN, N, WN,VRN)
Implicit None
Complex * 16, Dimension (N, N) :: TN, VL,VRN
Complex * 16, Dimension (N) :: WN
Integer :: LDA, LDVL, LDVR
Integer, Parameter :: LWORK = 5000
Integer :: INFO, N
Complex * 16, Dimension (LWORK) :: WORK
Double Precision, Dimension (2*N) :: RWORK
LDA = N
LDVL = N
LDVR = N

Call zgeev ('N', 'V', N, TN, LDA, WN, VL, LDVL, VRN, LDVR, WORK, &
& LWORK, RWORK, INFO)

call  order_WS_AND_WN(VRN,WN,N)
End

!=====================================================!
! a more general way to get the inversion of a complex matrix
! use the zgtri in lapack library.
!=====================================================!
Subroutine inver(A,N)
Implicit None
INTEGER::N
Integer,PARAMETER :: LWORK=100000
Complex * 16, Dimension (N, N) :: A! INPUT AND OUTPUT
Complex * 16, Dimension (LWORK) :: WORK
integer,dimension(N)::IPIV
integer::LDA,INFO
LDA=N
call  zgetrf (N,N,A,LDA,IPIV,INFO)
call  zgetri (N,A,LDA,IPIV,WORK,LWORK,INFO)
end
!!!!!!!!==========the end===========================!!!!!

!!!----EIGENVALUE AND EIGENVECTOR IN a decending order(bigger to smaller)----!!!!
subroutine order_WS_AND_WN(b,a,N)
implicit none
integer::N,I,J,K
COMPLEX*16,DIMENSION(N,N)::b
COMPLEX*16,DIMENSION(N)::a
COMPLEX*16::temp
COMPLEX*16,DIMENSION(N,1)::temp_vector
do i=1,N-1
do j=i+1,N
if (real(a(j)) .lt. real(a(i)))then
temp = a(i)
a(i) = a(j)
a(j) = temp
temp_vector(:,1)=b(:,i)
b(:,i)=b(:,j)
b(:,j)=temp_vector(:,1)
endif
enddo
enddo
end
!!---------THE END-----------------------------!!!!
