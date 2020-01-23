      program reader
      
      implicit none

      !Read in the number, type, and positions of atoms.

      integer:: num,i
      character (len=2)::atm(2)
      double precision::xcoord(2),ycoord(2),zcoord(2)
      double precision:: x,y,z

      open(unit = 2, file ="h2.XYZ")

      read(2,*)num

      do i=1,num
      read(2,*)atm(i),x,y,z
      xcoord(i)=x*1.8897259886        
      ycoord(i)=y*1.8897259886
      zcoord(i)=z*1.8897259886
      end do

      close(2)

             
      write(*,2000)"Number of atoms =",num
     
2000  format(a17,i2)

      do i=1,num
      write(*,1000)atm(i),xcoord(i),ycoord(i),zcoord(i)
      end do
1000  format(a2,"=",1x,3(f11.9,2x))
      
      !Basis Set Setup


      !Number of Gaussians used to form a contracted Gaussian orbital.

      integer::STOnG

      STOnG=3
      
      !Max quantum number

      integer::H

      H=1

      !Gaussian contraction coefficients

      !D(1,*) represents 1s, D(2,*) represents 2s.

      double precision::D(2,3)

      D(1,1)=0.444635
      D(1,2)=0.535328
      D(1,3)=0.154329
      D(2,1)=0.700115
      D(2,2)=0.399513
      D(2,3)=-0.0999672

      !Gaussian orbital exponents.

      !alpha(1,*) represents 1s, alpha(2,*) represents 2s.

      double precision:: alpha(2,3)

      alpha(1,1)=0.109818
      alpha(1,2)=0.405771
      alpha(1,3)=2.22766
      alpha(2,1)=0.0751386
      alpha(2,2)=0.231031
      alpha(2,3)=0.994203
      

      !Basis set size

      integer::B

      B=H

      !Number of electrons

      integer::N

      N=2

      !Atom charge

      integer::charge

      charge=1

      end program reader

      !Integrals between two Gaussian functions

      !Subroutine to calculate the product of two Gaussians.

      subroutine gauss_product(gauss_A,gauss_B,xi,yi,zi,p,diff,K,Rp)

      implicit none
      
      double precision,intent(in)::gauss_A,gauss_B,xi(2),yi(2),zi(2)
      double precision::a,b,Ra,Rb,N
      double precision,intent(out)::p,diff,K,Rp
      double precision,parameter::Pi=3.1415927
      
      gauss_A=a
      gauss_B=b
      p=a+b
      diff=(dsqrt((xi(2)-xi(1))**2+(yi(2)-yi(1))**2+(zi(2)-zi(1))**2)
     +**2)
      N=(4*a*b/(Pi**2))**0.75   
      K=N*exp(-a*b/p*diff)  
      Rp=(a*Ra+b*Rb)/p

      end subroutine gauss_product

      !Subroutine to calculate the overlap integral.

      subroutine overlap(gauss_A,gauss_B,xi,yi,zi,pre-k)
      double precision,intent(in)::gauss_A,gauss_B,xi(2),yi(2),zi(2)
      double precision::p,diff,K,Rp,prefactor
      double precision,intent(out)::pre-k
      double precision,parameter::Pi=3.1415927      

      call gauss_product(gauss_A,gauss_B,xi,yi,zi,p,diff,K,Rp)
      prefactor=(Pi/p)**1.5
      pre-k=prefactor*K

      end subroutine overlap

      !Subroutine to calculate the kinetic integral.

      subroutine kinetic (xi,yi,zi)
      double precision,intent(in)::xi(2),yi(2),zi(2)
      
      call gauss_product(
      
      

