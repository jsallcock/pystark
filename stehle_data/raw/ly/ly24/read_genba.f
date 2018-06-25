C Decembre 1998:
! how to read inside the tables of H-lines
! and to calculate the line shape
! for a given density DENS and temperature T
! (cm^-3 and K)
! id_max is the number of input densities
! or the number of input files profil1.dat...

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)

         parameter(id_maxi=30) ! maximum number of densities
      
	dimension tempe(10),jtot(id_maxi,10)
        dimension din(10,100),sprof(10,100),sprofs(10,100)
        dimension dl12(10),dl12s(10)
 
        dimension dom(id_maxi,10,100)
        dimension oline(id_maxi,10,100),olines(id_maxi,10,100)
        dimension dom0(10000)
        dimension domm(100),tprof(id_maxi,10,100),tprofs(id_maxi,10,100)
        dimension dense(id_maxi),f00(id_maxi)
        dimension pr0(id_maxi,10), uprof(2,100),uprofs(2,100)

        open(15,FILE='nraie.dat',status='old')
        read(15,1600) id_max
        close(15)

 1604	FORMAT(5(1X,D13.6))
 1600	FORMAT(I3,1X,I3)
 1601	FORMAT(2I3,F9.5,1X,D12.5)
 1602	FORMAT(2(D13.6))
 5044   FORMAT(1X)
 5045	FORMAT(1X,1PE10.3,1X,1PE10.3,1X,'(',1PE10.3,')',
     S 3X,1PE10.3,1X,'(',1PE10.3,')',
     S 3X,1PE10.3,1X,'(',1PE10.3,')')
 5054	format(24x,1pe10.3)
 7100	format(8x, 2(1x,i3))
 1801	FORMAT(21X,I2,11X,I2,26x,F8.2)
 1802	format(49x,e10.3)
 3000	format(1X,E10.3,1x,e10.3,2x,e10.3,3x,e10.3,2x,e10.3,
     s    3x,e10.3,2x,e10.3,2x,e10.3,2x,e10.3)
 6000	format(20X,1pe10.3,2(15x,1pe10.3))
 6002	format(20x,1pe10.3,2(15x,1pe10.3))
 5300	FORMAT(1X,1PE10.3,3(1x,1PE10.3,1x,'(',1pe10.3,' )'))
 8000	  format(1X,1PE10.3,3(1x,1pe10.3,1x,'(',1pe10.3,' )'))
 8003	  format(1X,1PE10.3,2(1x,1pe10.3,1x,'(',1pe10.3,' )'))
 8001  	  format(1x,1Pe10.3,25x,2(1x,1pe10.3,1x,'(',1pe10.3,' )'))
 8002	  format(1x,1Pe10.3,1x,1pe10.3,1x,'(',1pe10.3,' )',
     s     25x,2(1x,1pe10.3,1x,'(',1pe10.3,' )'))
 8005	  format(1x,1Pe10.3,24x,2(1x,1pe10.3,1x,'(',1pe10.3,' )'))
 8004  	  format(1x,1Pe10.3,50x,2(1x,1pe10.3,1x,'(',1pe10.3,' )'))
 8006	  format(1x,1Pe10.3,1x,1pe10.3,1x,'(',1pe10.3,' )')

	OPI=3.1415926536D00
        cspeed=2.9979d18 ! velocity of light in Ansgtroms/s
        cspeed_pi=2.*opi*cspeed

        read(5,*) TEMP ! temperature in 
        read(5,*) DENS ! electronic density in cm-3
        PR0_exp= 0.0898*(DENS**(1./6.))/dsqrt(TEMP) !=(r0/debye)
        F00_exp=1.25d-9*(DENS**(2./3.)) ! normal field value in ues
        if(PR0_exp.gt.1.) then
           write(91,*) 'PR0_exp',PR0_exp
           stop
        endif

        tprof(:,:,:)=0.
        tprofs(:,:,:)=0.

! each input file contains tables for a given density
! but 10 temperatures.
! the parameter pr0 = r0/debye= 0.0898 Ne^{1/6} T^{-1/2}, 
! (Ne in cm^-3, T in K) . For pr0 >1, we do not give a tabulation
! yet now 
          jtot(:,:)=0
          dom(:,:,:)=0.
          oline(:,:,:)=0.
          olines(:,:,:)=0
          dom0(:)=0.
          pr0(:,:)=0.

        do 444 id=1,id_max ! loop over the density
          write(90,*) id

          tempe(:)=0.
          dl12(:)=0.
          dl12s(:)=0
          din(:,:)=0.
          sprof(:,:)=0.
          sprofs(:,:)=0.
! input files 
! profil1.dat for N1
! profil2.dat for N2
! profil3.dat for N3 ..., with N1<N2<N3

        if(id.eq.1) then
  	  OPEN(16,FILE='profil1.dat',STATUS='old')
          open(17,file='index1.dat',status='old')
        endif
        if(id.eq.2) then
	  OPEN(16,FILE='profil2.dat',STATUS='old')
          open(17,file='index2.dat',status='old')
        endif
        if(id.eq.3) then
	  OPEN(16,FILE='profil3.dat',STATUS='old')
          open(17,file='index3.dat',status='old')
        endif
        if(id.eq.4) then
	  OPEN(16,FILE='profil4.dat',STATUS='old')
          open(17,file='index4.dat',status='old')
        endif
        if(id.eq.5) then
	  OPEN(16,FILE='profil5.dat',STATUS='old')
          open(17,file='index5.dat',status='old')
        endif
        if(id.eq.6) then
	  OPEN(16,FILE='profil6.dat',STATUS='old')
          open(17,file='index6.dat',status='old')
        endif
        if(id.eq.7) then
	  OPEN(16,FILE='profil7.dat',STATUS='old')
          open(17,file='index7.dat',status='old')
        endif
        if(id.eq.8) then
	  OPEN(16,FILE='profil8.dat',STATUS='old')
          open(17,file='index8.dat',status='old')
        endif
        if(id.eq.9) then
	  OPEN(16,FILE='profil9.dat',STATUS='old')
          open(17,file='index9.dat',status='old')
        endif
        if(id.eq.10) then
	  OPEN(16,FILE='profil10.dat',STATUS='old')
          open(17,file='index10.dat',status='old')
        endif
        if(id.eq.11) then
	  OPEN(16,FILE='profil11.dat',STATUS='old')
          open(17,file='index11.dat',status='old')
        endif
        if(id.eq.12) then
	  OPEN(16,FILE='profil12.dat',STATUS='old')
          open(17,file='index12.dat',status='old')
        endif
        if(id.eq.13) then
	  OPEN(16,FILE='profil13.dat',STATUS='old')
          open(17,file='index13.dat',status='old')
        endif
        if(id.eq.14) then
	  OPEN(16,FILE='profil14.dat',STATUS='old')
          open(17,file='index14.dat',status='old')
        endif
        if(id.eq.15) then
	  OPEN(16,FILE='profil15.dat',STATUS='old')
          open(17,file='index15.dat',status='old')
        endif
        if(id.eq.16) then
	  OPEN(16,FILE='profil16.dat',STATUS='old')
          open(17,file='index16.dat',status='old')
        endif
        if(id.eq.17) then
	  OPEN(16,FILE='profil17.dat',STATUS='old')
          open(17,file='index17.dat',status='old')
        endif
        if(id.eq.18) then
	  OPEN(16,FILE='profil18.dat',STATUS='old')
          open(17,file='index18.dat',status='old')
        endif
        if(id.eq.19) then
	  OPEN(16,FILE='profil19.dat',STATUS='old')
          open(17,file='index19.dat',status='old')
        endif
        if(id.eq.20) then
	  OPEN(16,FILE='profil20.dat',STATUS='old')
          open(17,file='index20.dat',status='old')
        endif


! input tables Lyman/Balmer, and indexes for reading inside
! the tabulations which are given in dalfa units
! dalfa=dlambda(Angstroms)/F00(ues)
! corespond to delta alpha >0 (ie delta omega=domega <0)
! Note that I(-domega) =I(domega
! F00 = 1.25 10^{-9} Ne^{2/3}  (cm-3,ues) = normal field strength
! fainu = wings factor in alfa units
!      ie (Idalfa)_wings= fainu/(dalfa**2.5)
!      ( but not to far in the wings)

	read(16,1801) n,np,olam0
        read(16,1802) dense(id) ! electronic density in cm^{-3}
        read(16,1802) F00(id)   ! = 1.25 10^{-9} Ne^{2/3}
        read(16,1802) fainu

	ON=N
	ONP=NP
	DNU=1.D0/ON**2 - 1.D0/ONP**2
	AMBDA=911.7633455*(ON*ONP)**2/((ONP-ON)*(ONP+ON))

        omega=cspeed_pi/ambda
        read(17,5044 )
        read(17,5044 )
        read(17,5044 )
        read(17,5044 )
        read(17,5044 )

        if(id.gt.17) go to 1717
	read(16,5044)
        read(16,6000) tempe(1),tempe(2),tempe(3)
        read(16,6002) pr0(id,1),pr0(id,2),pr0(id,3)
        read(16,6002) dl12s(1),dl12s(2),dl12s(3)
        read(16,6002) dl12(1),dl12(2),dl12(3)
        read(16,5044 )

        read(17,7100) ifm0,ifm1

         itot=ifm1-ifm0+1
         do i=1,itot
            read(16,3000) din(1,i),sprof(1,i),sprofs(1,i),
     s      sprof(2,i),sprofs(2,i),sprof(3,i),sprofs(3,i)
         end do

         din(2,:)=din(1,:)
         din(3,:)=din(1,:)

c=======> jtot(id,it)= number of wave lengths for the couple (T,Ne)
        do j=1,3
          do i=1,itot
            jtot(id,j)=i
            if(sprof(j,i).eq.0.) then
              jtot(id,j)=i-1
              exit
            endif
          end do
        end do

	read(16,5044)
        read(16,6000) tempe(4),tempe(5),tempe(6)
        read(16,6002) pr0(id,4),pr0(id,5),pr0(id,6)
        read(16,6002) dl12s(4),dl12s(5),dl12s(6)
        read(16,6002) dl12(4),dl12(5),dl12(6)
        read(16,5044 )

        read(17,7100) ifm0,ifm1
        itot=ifm1-ifm0+1
        do i=1,itot
           read(16,3000) din(4,i),sprof(4,i),sprofs(4,i),
     s     sprof(5,i),sprofs(5,i),sprof(6,i),sprofs(6,i)
         end do

        din(5,:)=din(4,:)
        din(6,:)=din(4,:)
c=======>
        do j=4,6
          do i=1,itot
            jtot(id,j)=i
            if(sprof(j,i).eq.0.) then
              jtot(id,j)=i-1
              exit
            endif
          end do
        end do
 1717	continue
	read(16,5044)
        read(16,6000) tempe(7),tempe(8),tempe(9)
        read(16,6002) pr0(id,7),pr0(id,8),pr0(id,9)
        read(16,6002) dl12s(7),dl12s(8),dl12s(9)
        read(16,6002) dl12(7),dl12(8),dl12(9)
        read(16,5044 )

        read(17,7100) ifm0,ifm1

        itot=ifm1-ifm0+1
        do i=1,itot
            read(16,3000) din(7,i),sprof(7,i),sprofs(7,i),
     s      sprof(8,i),sprofs(8,i),sprof(9,i),sprofs(9,i)
         end do

        din(8,:)=din(7,:)
        din(9,:)=din(7,:)
c=======>
        do j=7,9
          do i=1,itot
            jtot(id,j)=i
            if(sprof(j,i).eq.0.) then
              jtot(id,j)=i-1
              exit
            endif
          end do
        end do

	read(16,5044)
        read(16,6000) tempe(10)
        read(16,6002) pr0(id,10)
        read(16,6002) dl12s(10)
        read(16,6002) dl12(10)
        read(16,5044 )

        read(17,7100) ifm0,ifm1
        itot=ifm1-ifm0+1
        do i=1,itot
          read(16,3000) din(10,i),sprof(10,i),sprofs(10,i)
        end do
c=======>
        do j=10,10
          do i=1,itot
            jtot(id,j)=i
            if(sprof(j,i).eq.0.) then
              jtot(id,j)=i-1
              exit
            endif
          end do
        end do

!===
	Tempe(1) =   2500.
        Tempe(2) =   5000.
	Tempe(3) =  10000.
	Tempe(4) =  19950.
	Tempe(5) =  39810.
        Tempe(6) =  79430.
        Tempe(7) = 158500.
        Tempe(8) = 316200.
        Tempe(9) = 631000.
        Tempe(10)=1259600.

        do it=1,10
          if(pr0(id,it).eq.0.) then ! non tabulated case
            pr0(id,it)=0.0898*(dense(id)**(1./6.))/dsqrt(Tempe(it))
          endif
	end do

        PR00=0.0898*(dense(id)**(1./6.))/dsqrt(TEMP)


! conversion from alfa=dlambda(Angstroms)/F00 units
! to normalized domega units (domega/F00) in rd/(s,ues),
! normalized intensities are in units of (ues*s)/rd
c=======>
        do j=1,10
          do i=1,jtot(id,j)
! detunings for 1-np transition (alfa, omega, lambda units)
             dalfa=din(j,i)
             dlambda=f00(id)*dalfa   ! angstroms
             domega=-cspeed_pi*dlambda/((ambda+dlambda)*ambda) !rd/s
             dom(id,j,i)=domega/F00(id) ! (rd/(s*ues)
             dom(id,j,i)=dabs(dom(id,j,i))
             otrans=-cspeed_pi/(ambda*ambda)

             olines(id,j,i)=sprofs(j,i)/dabs(otrans)
             oline(id,j,i) =sprof(j,i)/dabs(otrans)
             write(73,1941) Tempe(j),dense(id),
     s          dom(id,j,i),olines(id,j,i),oline(id,j,i)
           end do
         end do
 1941	 format(5(1x,e10.3))
         write(73,*)
!  asymptotic wing factor in normalized omega units
! in the line wings : I(domega)=fainom/(domega**2.5)
! In these units (rd/(s*ues)), the asymptotic constant
! fainom is independnt of Ne and T
! in the wings, one has I=fainom/(dom**2.5)
         fainom=fainu* ((dabs(otrans))**1.5)
! in true angular frequency units one has
! in the wings, one has I=fainom_exp/(domega**2.5)
         fainom_exp=fainom*(F00_exp**1.5)
! in frequency units (1/s) nu=omega/(2*pi)
! in the wings, one has I=fainum_exp/(dnu**2.5)
         fainum_exp=fainom_exp/( (opi*2.)**1.5) ! pour Lyman et Balmer

	close(17)
     	close(16)
 444	continue

        
! define an unique detunings grid - domm -  for the tabulated
! profiles ( various temperatures , same density)
         inc=1
         domm(inc)=0.
         do id=1,id_max
            do j=1,10
               do i=2,jtot(id,j)
                 inc=inc+1
                 dom0(inc)=dom(id,j,i)
               end do
            end do
         end do
         npik=inc
         write(71,*)'npik',npik

         nut=10000
         call piksrt(npik,nut,dom0)

         inc=1
         domm(1)=0.
         do i=2,nut
           dif=(dom0(i)-dom0(i-1))
           if(dif.le.1.e-6) cycle
	   if((dif/dabs(dom0(i))).le.0.1) cycle
           inc=inc+1
           domm(inc)=dom0(i)
         end do
         jdom=inc
         write(71,*) 'inc',inc

        do 110 id=1,id_max 
          do 111 j=1,10
             if(pr0(id,j).gt.1.) go to 111
             write(71,*)
             tprof(id,j,1) =oline(id,j,1)
             tprofs(id,j,1)=olines(id,j,1)
             write(71,9100) F00(id),domm(1),oline(id,j,1),olines(id,j,1)
             do 112 i=2,jdom
                 domeg=domm(i)
                 ij_max=jtot(id,j)
                 do ij=2,ij_max-1
                    test=(domeg-dom(id,j,ij))*(domeg-dom(id,j,ij-1))
                    if(test.le.0.) then
                      x1=dom(id,j,ij-1)
                      x2=dom(id,j,ij)
                      x3=dom(id,j,ij+1)
                      y1=oline(id,j,ij-1)
                      y2=oline(id,j,ij)
                      y3=oline(id,j,ij+1)
                      tprof(id,j,i)= FINTRP(X1,X2,X3,Y1,Y2,Y3,domeg)
                      y1=olines(id,j,ij-1)
                      y2=olines(id,j,ij)
                      y3=olines(id,j,ij+1)
                      tprofs(id,j,i)= FINTRP(X1,X2,X3,Y1,Y2,Y3,domeg)
                      write(71,9100)
     s                  F00(id),domeg,tprof(id,j,i),tprofs(id,j,i)
                      go to 112
                    endif
		 end do

                 test=
     s            (domeg-dom(id,j,ij_max-1))*(domeg-dom(id,j,ij_max))
                 if(test.le.0.) then
                      x1=dom(id,j,ij_max-2)
                      x2=dom(id,j,ij_max-1)
                      x3=dom(id,j,ij_max)
                      y1=oline(id,j,ij_max-2)
                      y2=oline(id,j,ij_max-1)
                      y3=oline(id,j,ij_max)
                      tprof(id,j,i)= FINTRP(X1,X2,X3,Y1,Y2,Y3,domeg)
                      y1=olines(id,j,ij_max-2)
                      y2=olines(id,j,ij_max-1)
                      y3=olines(id,j,ij_max)
                      tprofs(id,j,i)= FINTRP(X1,X2,X3,Y1,Y2,Y3,domeg)
                      write(71,9100)
     s                  F00(id),domeg,tprof(id,j,i),tprofs(id,j,i)
                      go to 112
                    endif

                 if(domeg.gt.dom(id,j,ij_max)) then
                      tprof(id,j,i) =fainom/(domeg**2.5)
                      tprofs(id,j,i)=tprof(id,j,i)
                      write(71,9100)
     s                   F00(id),domeg,tprof(id,j,i),tprofs(id,j,i)
                      go to 112
                 endif
 112	 end do
 111	end do
 110	end do

        if(DENS.ge.10.*dense(id_max)) stop
        if(DENS.le.dense(1))      stop
        if(TEMP.ge.Tempe(10))     stop
        if(TEMP.le.Tempe(1))      stop

        do id=1,id_max-1
          otest_dens=(DENS-dense(id))*(DENS-dense(id+1))
          if (otest_dens.le.0.) then
            dense1=dense(id)
            dense2=dense(id+1)
            id1=id
            id2=id+1
            exit
          endif
        end do

         if(DENS.gt.dense(id_max)) then
           dense1=dense(id_max-1)
           dense2=dense(id_max)
           id1=id_max-1
           id2=id_max
        endif

        do it=1,9         
          otest=(TEMP-tempe(it))*(TEMP-tempe(it+1))
          if( otest.le.0. ) then
            it1=it
            it2=it+1
            pr01=pr0(id2,it1) ! max value of pr0 for T1,T2,dense1,dense2
            tempe1=tempe(it)
	    tempe2=tempe(it+1)
            exit
          endif
         end do

         if(pr01.gt.1.) stop

         write(81,*) ' temperature, density and r0/debye F00'
         write(81,9100) TEMP,DENS,PR0_exp,F00_exp
         write(81,*) ' interpolation gives T1 and T2, dense1 and dense2'
         write(81,9100) Tempe1,Tempe2,dense1,dense2
         write(81,*) id1,id2,it1,it2

 


! interpolation in temperature       
        do id=id1,id2
!            write(72,*) 'dense',dense(id),'T',TEMP
            do i=1,jdom
               uprof(id,i)=tprof(id,it1,i)+(TEMP-tempe1)*
     s            (tprof(id,it2,i)-tprof(id,it1,i))/(tempe2-tempe1)
               uprofs(id,i)=tprofs(id,it1,i)+(TEMP-tempe1)*
     s            (tprofs(id,it2,i)-tprofs(id,it1,i))/(tempe2-tempe1)
!               write(72,9100)
!     s             F00(id),domm(i),uprof(id,i),uprofs(id,i),fainom
            end do
        end do

! si on utilise Lyman pour Balmer on_new=2
! pour Paschen ce serait on_new=3
        on_nw=2.
        ambda_nw=911.7633455*(on_nw*ONP)**2/((ONP-on_nw)*(ONP+on_nw))
        otrans=-cspeed_pi/(ambda_nw*ambda_nw)



            write(72,*) 'DENS, T, fainum_exp'
            write(72,9100) DENS,TEMP,fainum_exp
            write(72,*) 'lAMBDA, n, np'
            write(72,9100) AMBDA_nw, ON_nw, ONP
            write(72,*)
 1666	    format(i3)
            write(74,*) 'DENS, T, fainum_exp'
            write(74,9100) DENS,TEMP,fainum_exp
            write(74,*) 'lAMBDA, n, np'
            write(74,9100) AMBDA_nw, ON_nw, ONP
            write(74,1666) jdom
            write(6,*) 'DENS, T, fainum_exp'
            write(6,9100) DENS,TEMP,fainum_exp
            write(6,*) 'lAMBDA, n, np'
            write(6,9100) AMBDA_nw, ON_nw, ONP
            write(6,*)

! interpolation in density
            do i=1,jdom
               wprof=uprof(id1,i)+(DENS-dense1)*
     s            (uprof(id2,i)-uprof(id1,i))/(dense2-dense1)
               wprofs=uprofs(id1,i)+(DENS-dense1)*
     s            (uprofs(id2,i)-uprofs(id1,i))/(dense2-dense1)
               write(72,9100)
     s             domm(i),wprof,wprofs,fainom
! 
! further conversions
! domega(rd.s-1) = F00_exp * domm
!        I(domega)(s/rd) =wprof/F00_exp
! dnu (s-1)      = domega/(2.*pi)
!        I(nu) (s)= I(domega) *(2*pi)
! dlambda(Angstroms) = - ambda * domega(rd s^_1) / (omega+domega)
!       different for domega and -domega
!       omega= cspeed_pi/ambda
            omega_nw=cspeed_pi/ambda_nw
             delta_omega=domm(i)*F00_exp
             delta_nu=delta_omega/(2*opi)
             delta_lambda= AMBDA_nw*delta_omega/(omega_nw + delta_omega)
             wprof_nu=(wprof/F00_exp)*(2.*opi)
             wprofs_nu=(wprofs/F00_exp)*(2.*opi)
             write(74,9100) delta_lambda,delta_nu,wprof_nu,wprofs_nu     
             write(6,9100) delta_lambda,delta_nu,wprof_nu,wprofs_nu
            end do

 9100	format(5(1x,e10.3))

	STOP
	END
C**********************************************************************
        DOUBLE PRECISION FUNCTION FINTRP(X1,X2,X3,Y1,Y2,Y3,X)
C
 	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
        IF(X.EQ.X2) GO TO 10
C
C L'INTERPOLATION EST :
C SI LA FONCTION EST MONOTONE (CROISSANTE OU DECROISSANTE)
C    HYPERBOLIQUE : Y=A+B/(X-C)
C    OU LINEAIRE  : Y=A*X+B
C SI LA FONCTION N'EST PAS MONOTONE
C    PARABOLIQUE  : Y=A*X*X+B*X+C
C
        A12=X1-X2
        A22=X1-X3
        V1=Y1-Y2
        V2=Y1-Y3
        IF(Y1.LT.Y2.AND.Y2.LT.Y3) GO TO 20
        IF(Y1.GT.Y2.AND.Y2.GT.Y3) GO TO 20
        X1C=X1*X1
        A11=X1C-X2*X2
        A21=X1C-X3*X3
        DETER=A11*A22-A12*A21
        IF(DABS(DETER).LT.1.D-40) GO TO 30
        A=(A22*V1-A12*V2)/DETER
        B=(-A21*V1+A11*V2)/DETER
        C=Y1-A*X1C-B*X1
        FINTRP=(A*X+B)*X+C
        RETURN
 20     DETER=V1*A22-V2*A12
        IF(DABS(DETER).LT.1.D-40) GO TO 40
        A21=X1*Y1
        A11=A21-X2*Y2
        A21=A21-X3*Y3
        C=(A22*A11-A12*A21)/DETER
        A=(-V2*A11+V1*A21)/DETER
        B=(Y1-A)*(X1-C)
        FINTRP=A+B/(X-C)
        RETURN
 10     FINTRP=Y2
        RETURN
 30     WRITE(5,100)
        STOP
 100    FORMAT(' ERREUR DANS FINTRP : DEUX DES XI SONT IDENTIQUES')
 40     FINTRP=Y1+(X-X1)*(Y3-Y1)/(X3-X1)
        RETURN
        END
C/////////////////////////////////////////////
	SUBROUTINE PIKSRT(N,nut,ARR)
C
C SORTING BY STRAIGHT INSERTION ; CF NUMERICAL RECIPES P 226
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION ARR(Nut)
	DO 12 J=2,N
	  A=ARR(J)
	  DO 11 I=J-1,1,-1
	  IF(ARR(I).LE.A) GO TO 10
	  ARR(I+1)=ARR(I)
  11	CONTINUE
	I=0
  10	ARR(I+1)=A
  12	CONTINUE
	RETURN
	END
C//////////////////////
