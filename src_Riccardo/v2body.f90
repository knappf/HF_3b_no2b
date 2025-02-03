module v2body

use technical

contains

  
subroutine Compress_2b(e1max,e2max,d1max,d2max,id_val,jmax_val)

integer :: e1max, e2max
integer :: d1max, d2max  
integer :: i,j,k,l,m
integer(kind=8) :: c2,c3,c4,c5
integer(kind=8) :: c_pp, c_nn, c_pn, n, o
integer(kind=8) :: tmp_lenght
real(kind=8), allocatable :: Vpp_tmp(:), Vnn_tmp(:), Vpn_tmp(:) 
integer(kind=8), allocatable :: PP_tmp(:), NN_tmp(:), PN_tmp(:)
integer :: id_val,jmax_val

tmp_lenght = id_val*id_val*id_val*id_val*(jmax_val+1)

allocate(Vpp_tmp(tmp_lenght))
allocate(PP_tmp(tmp_lenght))
allocate(Vnn_tmp(tmp_lenght))
allocate(NN_tmp(tmp_lenght))
allocate(Vpn_tmp(tmp_lenght))
allocate(PN_tmp(tmp_lenght))

Vpp_tmp=0.0
PP_tmp=0
Vnn_tmp=0.0
NN_tmp=0
Vpn_tmp=0.0
PN_tmp=0

c2=0
c3=0
c4=0
c5=0

c_pp=0
c_nn=0
c_pn=0

do i = 1, id_val
   do j = 1, id_val
      do k = 1, id_val
         do l = 1, id_val
            do m = 0, jmax_val
               c2=Cantor_pair(int(i,8),int(j,8))
               c3=Cantor_pair(c2,int(k,8))
               c4=Cantor_pair(c3,int(l,8))
               c5=Cantor_pair(c4,int(m,8))
               if (abs(Vpp(i,j,k,l,m)) > 1.d-9) then
                  c_pp=c_pp+1
                  Vpp_tmp(c_pp)=Vpp(i,j,k,l,m)
                  PP_tmp(c_pp)=c5
               endif   
               if (abs(Vnn(i,j,k,l,m)) > 1.d-9) then
                  c_nn=c_nn+1
                  Vnn_tmp(c_nn)=Vnn(i,j,k,l,m)
                  NN_tmp(c_nn)=c5
               endif   
               if (abs(Vpn(i,j,k,l,m)) > 1.d-9) then
                  c_pn=c_pn+1
                  Vpn_tmp(c_pn)=Vpn(i,j,k,l,m)
                  PN_tmp(c_pn)=c5
               endif
            enddo
         enddo  
      enddo    
   enddo
enddo

deallocate(Vpp,Vnn,Vpn)

allocate(Vpp_ar(c_pp))
allocate(Vpp_pair(c_pp))
Vpp_ar=0.0
Vpp_pair=0
do n = 1, c_pp
   Vpp_ar(n)=Vpp_tmp(n)
   Vpp_pair(n)=PP_tmp(n)
enddo
deallocate(Vpp_tmp)
deallocate(PP_tmp)

allocate(Vnn_ar(c_nn))
allocate(Vnn_pair(c_nn))
Vnn_ar=0.0
Vnn_pair=0
do n = 1, c_nn
   Vnn_ar(n)=Vnn_tmp(n)
   Vnn_pair(n)=NN_tmp(n)
enddo
deallocate(Vnn_tmp)
deallocate(NN_tmp)

allocate(Vpn_ar(c_pn))
allocate(Vpn_pair(c_pn))
Vpn_ar=0.0
Vpn_pair=0
do n = 1, c_pn
   Vpn_ar(n)=Vpn_tmp(n)
   Vpn_pair(n)=PN_tmp(n)
enddo
deallocate(Vpn_tmp)
deallocate(PN_tmp)

call sort2b(c_pp,Vpp_pair,Vpp_ar)
call sort2b(c_nn,Vnn_pair,Vnn_ar)
call sort2b(c_pn,Vpn_pair,Vpn_ar)

write(*,*)'Compression of V2B'
write(*,*)'from', tmp_lenght, tmp_lenght, tmp_lenght
write(*,*)'to', c_pp, c_nn, c_pn

!do i = 1, id_val
!   do j = 1, id_val
!      do k = 1, id_val
!         do l = 1, id_val
!            do m = 0, jmax_val
!               if (abs(Vpp(i,j,k,l,m)-Vpp_me(i,j,k,l,m)) > 1.d-9) then
!                  write(*,*)'Error at V2Bpp',i,j,k,l,m,Vpp(i,j,k,l,m),Vpp_me(i,j,k,l,m)
!               endif   
!               if (abs(Vnn(i,j,k,l,m)-Vnn_me(i,j,k,l,m)) > 1.d-9) then
!                  write(*,*)'Error at V2Bnn',i,j,k,l,m,Vnn(i,j,k,l,m),Vnn_me(i,j,k,l,m)
!               endif   
!               if (abs(Vpn(i,j,k,l,m)-Vpn_me(i,j,k,l,m)) > 1.d-9) then
!                  write(*,*)'Error at V2Bpn',i,j,k,l,m,Vpn(i,j,k,l,m),Vpn_me(i,j,k,l,m)
!               endif   
!            enddo
!         enddo  
!      enddo    
!   enddo
!enddo

end subroutine Compress_2b


subroutine Unpack_2b(id_val,jmax_val)

integer :: i,j,k,l,m
integer(kind=8) :: c2,c3,c4,c5
integer(kind=8) :: c_pp, c_nn, c_pn, n, o
integer :: id_val,jmax_val

allocate(Vpp(id_val,id_val,id_val,id_val,0:jmax_val))
allocate(Vnn(id_val,id_val,id_val,id_val,0:jmax_val))
allocate(Vpn(id_val,id_val,id_val,id_val,0:jmax_val))

Vpp=0.0
Vnn=0.0
Vpn=0.0

do i = 1, id_val
   do j = 1, id_val
      do k = 1, id_val
         do l = 1, id_val
            do m = 0, jmax_val
               Vpp(i,j,k,l,m)=Vpp_me(i,j,k,l,m)
               Vnn(i,j,k,l,m)=Vnn_me(i,j,k,l,m)
               Vpn(i,j,k,l,m)=Vpn_me(i,j,k,l,m)
            enddo
         enddo  
      enddo    
   enddo
enddo

deallocate(Vpp_ar,Vnn_ar,Vpn_ar)
deallocate(Vpp_pair,Vnn_pair,Vpn_pair)

end subroutine Unpack_2b



function Vpp_no(i) result(val)
implicit none
integer(kind=8) :: arrl
integer(kind=8) :: low,high,middle,rang
integer(kind=8) :: i
real(kind=8) :: val

val=0.0

arrl=size(Vpp_ar,DIM=1,KIND=8)
low=1
high=arrl
rang=high-low
middle=(low+high)/2

   do while((Vpp_pair(middle).ne.i).and.(rang.gt.0))
        if(i>Vpp_pair(middle)) then
           low=middle+1
        else
           high=middle-1
        end if
        
        rang=high-low
        middle=(high+low)/2
   end do
   if(Vpp_pair(middle)==i) then
       val=Vpp_ar(middle)
   end if

 end function Vpp_no
 
function Vpp_me(i,j,k,l,m) result(res)

integer :: i,j,k,l,m
integer(kind=8) :: c2,c3,c4,c5
real(kind=8) :: res

res=0.0

c2=0
c3=0
c4=0
c5=0

c2=Cantor_pair(int(i,8),int(j,8))
c3=Cantor_pair(c2,int(k,8))
c4=Cantor_pair(c3,int(l,8))
c5=Cantor_pair(c4,int(m,8))

res=Vpp_no(c5)

end function Vpp_me



function Vnn_no(i) result(val)
implicit none
integer(kind=8) :: arrl
integer(kind=8) :: low,high,middle,rang
integer(kind=8) :: i
real(kind=8) :: val

val=0.0

arrl=size(Vnn_ar,DIM=1,KIND=8)
low=1
high=arrl
rang=high-low
middle=(low+high)/2


   do while((Vnn_pair(middle).ne.i).and.(rang.gt.0))
        if(i>Vnn_pair(middle)) then
           low=middle+1
        else
           high=middle-1
        end if
        
        rang=high-low
        middle=(high+low)/2
   end do
   if(Vnn_pair(middle)==i) then
       val=Vnn_ar(middle)
   end if

end function Vnn_no

function Vnn_me(i,j,k,l,m) result(res)

integer :: i,j,k,l,m
integer(kind=8) :: c2,c3,c4,c5
real(kind=8) :: res

res=0.0

c2=0
c3=0
c4=0
c5=0

c2=Cantor_pair(int(i,8),int(j,8))
c3=Cantor_pair(c2,int(k,8))
c4=Cantor_pair(c3,int(l,8))
c5=Cantor_pair(c4,int(m,8))

res=Vnn_no(c5)

end function Vnn_me



function Vpn_no(i) result(val)
implicit none
integer(kind=8) :: arrl
integer(kind=8) :: low,high,middle,rang
integer(kind=8) :: i
real(kind=8) :: val

val=0.0

arrl=size(Vpn_ar,DIM=1,KIND=8)
low=1
high=arrl
rang=high-low
middle=(low+high)/2


   do while((Vpn_pair(middle).ne.i).and.(rang.gt.0))
        if(i>Vpn_pair(middle)) then
           low=middle+1
        else
           high=middle-1
        end if
        
        rang=high-low
        middle=(high+low)/2
   end do
   if(Vpn_pair(middle)==i) then
       val=Vpn_ar(middle)
   end if

end function Vpn_no

function Vpn_me(i,j,k,l,m) result(res)

integer :: i,j,k,l,m
integer(kind=8) :: c2,c3,c4,c5
real(kind=8) :: res

res=0.0

c2=0
c3=0
c4=0
c5=0

c2=Cantor_pair(int(i,8),int(j,8))
c3=Cantor_pair(c2,int(k,8))
c4=Cantor_pair(c3,int(l,8))
c5=Cantor_pair(c4,int(m,8))

res=Vpn_no(c5)

end function Vpn_me



function Cantor_pair(ia,ja) result(ra)

integer(kind=8) :: ia,ja,ra,rr

ra=0
rr=0.0
rr=(ia+ja)*(ia+ja+1)
ra=rr/2+ja

end function Cantor_pair


     SUBROUTINE sort2b(n,arr,brr)
     INTEGER (kind=8) M,NSTACK
     PARAMETER (M=7,NSTACK=500)
     INTEGER (kind=8) n,i,ir,j,jstack,k,l,istack(NSTACK)
     INTEGER (kind=8) a, tempa
     REAL (kind=8) b, tempb
     INTEGER (kind=8) arr(n)
     REAL (kind=8) brr(n)
!Sorts an array arr(1:n) into ascending order using Quicksort, while making the corresponding rearrangement of the array brr(1:n).
     jstack=0
     l=1
     ir=n

1    if(ir-l.lt.M) then
        do j=l+1,ir
           a=arr(j)
           b=brr(j)
           do i=j-1,l,-1
              if(arr(i).le.a) goto 2
              arr(i+1)=arr(i)
              brr(i+1)=brr(i)
           enddo
           i=l-1
2          arr(i+1)=a
           brr(i+1)=b
        enddo
        if(jstack.eq.0) return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
     else
        k=(l+ir)/2
        tempa=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=tempa
        tempb=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=tempb
        if(arr(l).gt.arr(ir)) then
           tempa=arr(l)
           arr(l)=arr(ir)
           arr(ir)=tempa
           tempb=brr(l)
           brr(l)=brr(ir)
           brr(ir)=tempb
        endif
        if (arr(l+1).gt.arr(ir)) then
           tempa=arr(l+1)
           arr(l+1)=arr(ir)
           arr(ir)=tempa
           tempb=brr(l+1)
           brr(l+1)=brr(ir)
           brr(ir)=tempb
        endif
        if(arr(l).gt.arr(l+1)) then
           tempa=arr(l)
           arr(l)=arr(l+1)
           arr(l+1)=tempa
           tempb=brr(l)
           brr(l)=brr(l+1)
           brr(l+1)=tempb
        endif
        i=l+1
        j=ir
        a=arr(l+1)
        b=brr(l+1)
3       continue
           i=i+1
        if(arr(i).lt.a) goto 3
4       continue
           j=j-1
        if(arr(j).gt.a) goto 4
        if(j.lt.i) goto 5
        tempa=arr(i)
        arr(i)=arr(j)
        arr(j)=tempa
        tempb=brr(i)
        brr(i)=brr(j)
        brr(j)=tempb
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        brr(l+1)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK) pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l) then
           istack(jstack)=ir
           istack(jstack-1)=i
           ir=j-1
        else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
        endif
     endif   
     goto 1
     END SUBROUTINE sort2b


end module v2body
