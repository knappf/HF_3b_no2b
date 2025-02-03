      subroutine interaction_binary

       USE technical
       USE geom
       use Tcm_2body
       USE v3body_no2b
       use NN_bin
       use NNN_no2b_bin
!      Here for compression of 2b MEs       
       use v2body

       include 'define.inc'
       include 'parameters.inc'
       include 'commons.inc'

       double precision :: unitar_matrix(id,id)

       integer(kind=1) :: ija,ita
       integer(kind=1) :: j2_abFF,j2_deFF,J2_3FF,it_abFF,it_deFF,IT2_3FF
       integer(kind=2) :: i1a,i2a,i3a,i4a,i5a,i6a
       real(kind=4) :: vint1,vint4
!      Here for compression of 2b MEs       
       integer :: id_val,jmax_val
       
       character*2 :: ch

       id2=2*id

       allocate(tran_p(id,id),tran_n(id,id))
       tran_p=0.d0
       tran_n=0.d0


       allocate(Vpp(id,id,id,id,0:jmax))
       allocate(Vnn(id,id,id,id,0:jmax))
       allocate(Vpn(id,id,id,id,0:jmax))
       Vpp=0.0
       Vnn=0.0
       Vpn=0.0
!      Here for compression of 2b MEs       
       id_val=id
       jmax_val=jmax
       
!       allocate(Vpp_2b(id,id,id,id,0:jmax))
!       allocate(Vnn_2b(id,id,id,id,0:jmax))
!       allocate(Vpn_2b(id,id,id,id,0:jmax))
!       Vpp_2b=0.d0
!       Vnn_2b=0.d0
!       Vpn_2b=0.d0

       ACM = dble(AZ+AN)
 
       write(*,*)'File space:',noscmax,noscmax12,noscmax123
       write(*,*)'Model space:', dmax,dmax12,dmax123

!       if   (1==1) then
!***********************************************************************
!      Reading of 2-body NN elements from binary file generated by NuHamil code
       call read_2b_bin(noscmax,noscmax12,hbarom,ACM,dmax,dmax12)

!    2-body part of CM Hamiltonian
!    be careful with the factor of dsqrt(2.d0) 
!    

       icm_2b=1

       if (icm_2b /= 0) then

       write(*,*)'Calculation of 2-body part of CM Hamiltonian'
        factab=1.d0
        factcd=1.d0

!!$OMP PARALLEL DEFAULT(SHARED) 
!!$OMP& PRIVATE(i,j,k,l,JJ,factab,factcd,jmax_ij_min,jmax_ij_max)
!!$OMP DO
        do i=1,id
        write(*,*)'i=',i,id
         do j=1,id
          factab=1.d0
          if(i.eq.j) factab=dsqrt(2.d0)
          jmax_ij_min=iabs((levn(i)%j2-levn(j)%j2)/2)
          jmax_ij_max=(levn(i)%j2+levn(j)%j2)/2
          do k=1,id
           do l=1,id
            
            if ((100000*i+j) >= (100000*k+l)) then
             
            factcd=1.d0
            if(k.eq.l) factcd=dsqrt(2.d0)

!             do JJ=0,jmax
             do  JJ=jmax_ij_min,jmax_ij_max

!       2-body part of CM Hamiltonian
!         if (dabs(Vpn(i,j,k,l,JJ)) > 1.d-10) then 
          Vpn(i,j,k,l,JJ)=Vpn(i,j,k,l,JJ)
     &   +hbarom/ACM*T_nn(i,j,k,l,JJ)
        
          Vpn(k,l,i,j,JJ) = Vpn(i,j,k,l,JJ)

!         endif 
!         if (dabs(Vpp(i,j,k,l,JJ)) > 1.d-10) then 
          Vpp(i,j,k,l,JJ)=Vpp(i,j,k,l,JJ)
     &   +hbarom/ACM*T_nn_asym(i,j,k,l,JJ)
     &   *factab*factcd  
!         endif

          Vpp(k,l,i,j,JJ) = Vpp(i,j,k,l,JJ)

!         if (dabs(Vnn(i,j,k,l,JJ)) > 1.d-10) then
          Vnn(i,j,k,l,JJ)=Vnn(i,j,k,l,JJ)
     &    +hbarom/ACM*T_nn_asym(i,j,k,l,JJ)  
     &    *factab*factcd

           Vnn(k,l,i,j,JJ) = Vnn(i,j,k,l,JJ)
!          endif

            enddo

            endif 
           enddo
          enddo
         enddo
        enddo 
!!$OMP END DO
!!$OMP END PARALLEL

      endif

!***********************************************************************
!      Compressing of 2-body  elements
       call Compress_2b(noscmax,noscmax12,dmax,dmax12,id_val,jmax_val)     
!***********************************************************************
!      Reading of 3-body NO2B elements from binary file generated by NuHamil code
       call read_3b_no2b_bin(noscmax,noscmax12,noscmax123,
     & dmax,dmax12,dmax123)
***********************************************************************
!    Here the DD interaction elements were calculated                 *
!    obsolete in the new version of the code                          *
!    shoule be put to another place


       unitar_matrix=0.d0
       do i1=1,id
        unitar_matrix(i1,i1)=1.d0
       enddo

       igrid=20      ! number of the node points
       igrid2=250     ! number of the grid points

!       sizebox=10.d0    ! the interval in which the radial integral will be numerically summed
       sizebox=max(6.d0,2.0d0*dble(AZ+AN)**0.33333333d0)

       bos1=dsqrt(0.5d0*(zmp+zmn)*hbarom/(hbarc**2.d0))
       bos2=dsqrt(zmY*hbarom/(hbarc**2.d0))

       dx=sizebox/dble(igrid2)

       allocate(zcross(igrid))
       zcross(1)=0.070539889692d0 
       zcross(2)=0.372126818002d0
       zcross(3)=0.916582102483d0
       zcross(4)=1.70730653103d0
       zcross(5)=2.74919925531d0
       zcross(6)=4.04892531384d0
       zcross(7)=5.61517497087d0
       zcross(8)=7.45901745389d0
       zcross(9)=9.59439286749d0
       zcross(10)=12.0388025566d0
       zcross(11)=14.8142934155d0
       zcross(12)=17.9488955686d0
       zcross(13)=21.4787881904d0
       zcross(14)=25.4517028094d0
       zcross(15)=29.9325546634d0
       zcross(16)=35.0134341868d0
       zcross(17)=40.8330570974d0
       zcross(18)=47.6199940299d0
       zcross(19)=55.8107957541d0
       zcross(20)=66.5244165252d0
      
       do i=1,igrid
        zcross(i)=dsqrt(zcross(i))/bos1
       enddo

       allocate(weight(igrid))
       weight(1)=0.181080062419d0
       weight(2)=0.422556767879d0
       weight(3)=0.666909546702d0
       weight(4)=0.9153523727d0
       weight(5)=1.1695397071d0
       weight(6)=1.43135498624d0
       weight(7)=1.7029811359d0
       weight(8)=1.98701589585d0
       weight(9)=2.28663576323d0
       weight(10)=2.60583465152d0
       weight(11)=2.94978381794d0
       weight(12)=3.32539569477d0
       weight(13)=3.74225636246d0
       weight(14)=4.21424053477d0
       weight(15)=4.76252016007d0
       weight(16)=5.42172779036d0
       weight(17)=6.25401146407d0
       weight(18)=7.38731523837d0
       weight(19)=9.15132879607d0
       weight(20)=12.8933886244d0
  
       do i=1,igrid
        weight(i)=weight(i)/(2.d0*zcross(i)*bos1**2.d0)
       enddo

       allocate(rad_den1(igrid))
       rad_den1=0.d0

       do i=1,igrid
        val=0.d0
        radi=zcross(i)   !dble(i)*sizebox/dble(igrid)
        do j=1,id
         val_p=0.d0
         val_n=0.d0
         do k=1,id
          val_p=val_p+unitar_matrix(k,j)
     &                *R_val(levp(k)%nn,levp(k)%l,levp(k)%j2,bos1,radi)
          val_n=val_n+unitar_matrix(k,j)
     &                *R_val(levn(k)%nn,levn(k)%l,levn(k)%j2,bos1,radi)
         enddo
         val=val+lhfp(j)%vi**2.d0*val_p**2.d0*dble(lhfp(j)%j2+1)
         val=val+lhfn(j)%vi**2.d0*val_n**2.d0*dble(lhfn(j)%j2+1)
        enddo
        rad_den1(i)=val/(4.d0*pi)
       enddo


       return
      end
