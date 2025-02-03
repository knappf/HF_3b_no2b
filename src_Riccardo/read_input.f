      subroutine read_input

       include 'define.inc'
       include 'parameters.inc' 

       open(unit=1, file='input.dat', status= 'old', form= 'formatted')
        read(1,*) i,j
        read(1,*) k,k12,k123
        read(1,*) kk,kk12,kk123
        read(1,*) if_2b
        read(1,*) k1,k2,k3,k4
        read(1,*) if_q
        read(1,*) d
        read(1,*) e
        read(1,*) qp,qn,qpn,F0
        read(1,*) v1
        read(1,*) if_o
       close(1)

       AZ = j
       AN = i-j
       noscmax = k
       noscmax12 = k12
       noscmax123 = k123
       dmax = kk
       dmax12 = kk12
       dmax123 = kk123
       min_p = k1
       max_p = k2
       min_n = k3
       max_n = k4
       if_QTDA = if_q
       precis = d
       hbarom = e
       quenp = qp
       quenn = qn
       quenpn = qpn
       V3b_NNN = v1
       F_0=F0
       if_ort = if_o

       write(*,*) 'HF input parameters:'
       write(*,*) 'A=',AZ+AN,'Z=',AZ
       write(*,*) 
       write(*,*) 'Model space:'
       write(*,*) 'noscmax1=',noscmax
       write(*,*) 'noscmax2=',noscmax12
       write(*,*) 'noscmax3=',noscmax123       
       write(*,*) 'dmax1=',dmax
       write(*,*) 'dmax2=',dmax12
       write(*,*) 'dmax3=',dmax123       
       write(*,*)
!       write(*,*) 'nosc_TDA=',nosc_TDA
       write(*,*) 'precis=',precis
       write(*,*) 'hbarom=',hbarom,' MeV'
       write(*,*) 'CM_ort=',if_ort
       write(*,*) 
!       write(*,*) 'quenching factors =',quenp,quenn,quenpn,F_0

       if(mod(AZ,2).eq.1) then
        write(*,*) 'Error! Wrong input AZ - just even nuclei allowed.'
        stop
       endif

       if(mod(AN,2).eq.1) then
        write(*,*) 'Error! Wrong input AN - just even nuclei allowed.'
        stop
       endif

       if(dmax123.gt.noscmax123) then
        write(*,*) 'Error! Mismatch of model space and file.'
        stop
       endif
       
       if(AZ.eq.2.or.(AZ.eq.8.or.(AZ.eq.20.or.(AZ.eq.28.or.
     &                   (AZ.eq.50.or.(AZ.eq.82.or.AZ.eq.126)))))) then
        ifp_hfb = .false.
       else
        ifp_hfb = .true.
       endif

       if(AN.eq.2.or.(AN.eq.8.or.(AN.eq.20.or.(AN.eq.28.or.
     &                   (AN.eq.50.or.(AN.eq.82.or.AN.eq.126)))))) then
        ifn_hfb = .false.
       else
        ifn_hfb = .true.
       endif

       if(ifp_hfb.and.ifn_hfb) then
        write(*,*) 
        write(*,*) 'This option (open-shell nucleus in both protons and 
     &neutrons) '
        write(*,*) 'was not tested and may be it will not run properly.'
        write(*,*) 
!        stop
       endif

       if(if_QTDA.ne.0.and.if_QTDA.ne.1) then
        write(*,*) 'The if_QTDA must be =0 or =1!!!'
        stop
       endif

       if(if_QTDA.eq.0.and.(ifp_hfb.or.ifn_hfb)) then
        write(*,*) 'The particle-hole TDA cannot be calculated if 
     &non-zero pairing is active!!!'
        stop
       endif

       if(.not.(if_2b.eq.0.or.if_2b.eq.1)) then
        write(*,*) 'if_2b must be = 0 or = 1 !! '
        stop
       endif

      return
      end
