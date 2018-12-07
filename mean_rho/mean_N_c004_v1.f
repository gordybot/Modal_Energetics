      PROGRAM mean_rho
       
      IMPLICIT NONE
      
c---USED TO TEST THE EXTRACTION OF A SINGLE FIELD: THKNSS---------      
c mean_N_c004_v1.f, MCB, USM, 2017-11-13
c   for 1/25 grid
c mean_N_c008_v1.f, USM, 2017-9-11
c   load T and S tiles and compute N
c   at local potential density
c mean_rho_c008_v1.f, USM, 2017-04-27
c global extraction of z gridded density
c loads a global *.a file, adds halo, and saves chunks
       
      integer idm,jdm,kdm,blki,blkj,halocd,halo
      integer jblks,jblke !1,35
      integer iblks,iblke !1,60
      real flag
      real maxvar 
      real rho
      dimension rho(1:2)

c CHANGE THIS!!
      real, parameter :: dz=25.0

c glbc008
c      parameter(idm=4500,jdm=3528,kdm=240) 
c      parameter (blki=150,blkj=196,halo=3)
c glbc004
      parameter(idm=9000,jdm=7055,kdm=240) !150 x 60 = 9000 and 200 x 35 = 7000
      parameter (blki=150,blkj=200,halo=3)

      parameter(flag=2.0**100)
      real, parameter :: grav = 9.81
      character flnm*240,outfil*240,flnm2*240,flnm3*240,flnm4*240
      character flnm5*240
      character blank*5,runid*3,intid*2
      character*2 numstr(90)
      integer i,j,k,n2drec,nrecl,numfls,iblk,jblk,n,ic,flip_mcb,cc,kk
      integer l,m,dd,ee

      INTEGER I0,ID
      dimension I0(1:kdm+1)
      dimension ID(1:kdm+1)


c 2D variables
      integer blk_coast
      dimension blk_coast(1-halo:blki+halo,1-halo:blkj+halo)

      real blk_depth
      dimension blk_depth(1-halo:blki+halo,1-halo:blkj+halo)

c 3D variables
      real blk_Savg,blk_Tavg,N2
      dimension blk_Savg(1-halo:blki+halo,1-halo:blkj+halo,1:kdm)
      dimension blk_Tavg(1-halo:blki+halo,1-halo:blkj+halo,1:kdm)
      dimension N2(1-halo:blki+halo,1-halo:blkj+halo,1:kdm+1)

      character infil10*240,infil11*240
      character outfile13*60



c == EQ of STATE == EQ of STATE == EQ of STATE == EQ of STATE == EQ of STATE ==

      real salin,temp,c1,c2,c3,c4,c5,c6,c7,c8,c9,rc6
      integer, parameter ::
     &  sigver=6  !17-term sigma-2
      real    sig_n,sig_d,sig_q, dsigdt_n,dsigdt_d, dsigds_n,dsigds_d
      real    kappaf1
      real    sigloc_n,sigloc_d,sigloc_q,
     &        dsiglocdt_n,dsiglocdt_d, dsiglocds_n,dsiglocds_d
      real    r,pdb,prs
      integer kkf
      real, parameter ::
     &   aone =1.0,
     &   ahalf=1.0/2.0,
     &   a3rd =1.0/3.0, athird =a3rd,
     &   a4th =1.0/4.0, afourth=a4th
      real, parameter ::
     &   c001= 9.9984085444849347d+02,     !num. constant    coefficent
     &   c002= 7.3471625860981584d+00,     !num.    T        coefficent
     &   c003=-5.3211231792841769d-02,     !num.    T^2      coefficent
     &   c004= 3.6492439109814549d-04,     !num.    T^3      coefficent
     &   c005= 2.5880571023991390d+00,     !num.       S     coefficent
     &   c006= 6.7168282786692355d-03,     !num.    T  S     coefficent
     &   c007= 1.9203202055760151d-03,     !num.       S^2   coefficent
     &   c008= 1.0000000000000000d+00,     !den. constant    coefficent
     &   c009= 7.2815210113327091d-03,     !den.    T        coefficent
     &   c010=-4.4787265461983921d-05,     !den.    T^2      coefficent
     &   c011= 3.3851002965802430d-07,     !den.    T^3      coefficent
     &   c012= 1.3651202389758572d-10,     !den.    T^4      coefficent
     &   c013= 1.7632126669040377d-03,     !den.       S     coefficent
     &   c014= 8.8066583251206474d-06,     !den.    T  S     coefficent
     &   c015= 1.8832689434804897d-10,     !den.    T^3S     coefficent
     &   c016= 5.7463776745432097d-06,     !den.    T  S^1.5 coefficent
     &   c017= 1.4716275472242334d-09      !den.    T^3S^1.5 coefficent
      real, parameter ::
     &   c018= 1.1798263740430364d-02,     !num. P           coefficent
     &   c019= 9.8920219266399117d-08,     !num. P  T^2      coefficent
     &   c020= 4.6996642771754730d-06,     !num. P     S     coefficent
     &   c021= 2.5862187075154352d-08,     !num. P^2         coefficent
     &   c022= 3.2921414007960662d-12,     !num. P^2T^2      coefficent
     &   c023= 6.7103246285651894d-06,     !den. P           coefficent
     &   c024= 2.4461698007024582d-17,     !den. P^2T^3      coefficent
     &   c025= 9.1534417604289062d-18      !den. P^3T        coefficent
c --- additional coefficients for dsiglocdt().
      real, parameter ::
     &   c031= 7.3471625860981580d+00,     !num. constant    coefficent
     &   c032=-1.0642246358568354d-01,     !num.    T        coefficent
     &   c033= 1.0947731732944364d-03,     !num.    T^2      coefficent
     &   c034= 6.7168282786692355d-03,     !num.       S     coefficent
     &   c035= 7.2815210113327090d-03,     !den. constant    coefficent
     &   c036=-8.9574530923967840d-05,     !den.    T        coefficent
     &   c037= 1.0155300889740728d-06,     !den.    T^2      coefficent
     &   c038= 5.4604809559034290d-10,     !den.    T^3      coefficent
     &   c039=-8.8066583251206470d-06,     !den.       S     coefficent
     &   c040= 5.6498068304414700d-10,     !den.    T^2S     coefficent
     &   c041= 2.9432550944484670d-09,     !den.    T  S^1.5 coefficent
     &   c042= 1.9784043853279823d-07,     !num. P  T        coefficent
     &   c043= 6.5842828015921320d-12,     !num. P^2T        coefficent
     &   c044= 7.3385094021073750d-17,     !den. P^2T^2      coefficent
     &   c045= 9.1534417604289060d-18      !den. P^3         coefficent
c --- additional coefficients for dsiglocds().
      real, parameter ::
     &   c051= 2.5880571023991390d+00,     !num. constant    coefficent
     &   c052= 6.7168282786692355d-03,     !num.    T        coefficent
     &   c053= 3.8406404111520300d-03,     !num.       S     coefficent
     &   c054= 1.7632126669040377d-03,     !den. constant    coefficent
     &   c055=-8.8066583251206470d-06,     !den.    T        coefficent
     &   c056= 1.8832689434804897d-10,     !den.    T^3      coefficent
     &   c057= 8.6195665118148150d-06,     !den.       S^0.5 coefficent
     &   c058= 2.2074413208363504d-09,     !den.    T^2S^0.5 coefficent
     &   c059= 4.6996642771754730d-06      !num. P           coefficent
c
      real, parameter :: sqrmin=0.d0       !sqrt arg can't be negative
c --- reference pressure.
      real, parameter :: prs2pdb=1.d-4     !Pascals to dbar
      real pref                            !ref. pressure in Pascals, sigma2
      real rpdb                            !ref. pressure in dbar
c --- coefficients for 17-term rational function sig() at rpdb.
      real c101,c103,c105,c108,c109,c111
c --- additional coefficients for dsigdt().
      real c132,c135,c137           
c --- additional coefficients for dsigds().
      real c151            
c
c --- coefficients for kappa^(theta)
c --- new values (w.r.t. t-toff,s-soff,prs) from Shan Sun, Sep.2004
c --- 1=Arctic/Antarctic; 2=Atlantic; 3=Mediterranean
      real, parameter ::
     &   rhoref=1.d3  !rhoref=qthref kg/m^3
      real, parameter ::
     &   sclkap=1.e-11
      real, parameter, dimension(3) ::
     &  toff = (/  0.0,             3.0,            13.0 /)
     & ,soff = (/ 34.5,            35.0,            38.5 /)
     & ,qttt = (/ -3.03869354E-05, -3.03869352E-05, -3.03869353E-05 /)
     & ,qtt  = (/  4.56625601E-03,  4.29277358E-03,  3.38116552E-03 /)
     & ,qt   = (/ -2.88801209E-01, -2.61828868E-01, -1.81335007E-01 /)
     & ,qs   = (/ -1.08670290E-01, -1.05131061E-01, -9.33336309E-02 /)
     & ,qst  = (/  7.90503772E-04,  7.71096940E-04,  1.07270585E-03 /)
     & ,qpt  = (/  1.07813750E-09,  1.00638435E-09,  7.57239852E-10 /)
     & ,qpst = (/  1.41541548E-11,  1.48598578E-11,  3.89226107E-12 /)
     & ,qptt = (/ -1.31383708E-11, -1.31383707E-11, -1.31383708E-11 /)

c == EQ of STATE == EQ of STATE == EQ of STATE == EQ of STATE == EQ of STATE ==

      data blank/'     '/
      data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/


! read begin and end of blocks
      read(*,*) jblks
      read(*,*) jblke
      read(*,*) iblks
      read(*,*) iblke
      read(*,*) runid
      read(*,*) intid

      write (6,'(a,a3)') ' runid: ',runid
      call flush(6)

      write (6,'(a,a3)') ' intid: ',intid
      call flush(6)

      write(6,'(a)') blank 
      call flush(6)
      write (6,'(a,4i3)') ' jblks jblke iblks iblke: ',jblks,jblke,
     & iblks,iblke
      call flush(6)

!! loop over blocks ----------------------------------------------------
      do jblk=jblks,jblke
        do iblk=iblks,iblke
!! loop over blocks ----------------------------------------------------
          write(6,*) 'blk_',jblk,iblk
          call flush(6)

          write(6,*) 'opening files for reading and writing ...'
          call flush(6)

!!--------------------------------------------------------------------------
!! OPEN 2D GRIDDATA FILES            
!!--------------------------------------------------------------------------

! open the time invariant data files and read data these files 
! are written as 3-dimesional arrays
! (lon,lat,depth) - (i,j,k)

! coast
         infil10='griddata/coast_'//runid//
     &'_blk_'
     &                //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

          write(6,*) 'read',infil10
          open(unit=10,file=infil10,status='old',
     &                           form='unformatted')

          read(10) blk_coast

          close(10)

! depth
         infil11='griddata/'
     &  //'depth_'//runid//
     &'_blk_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF' 

          write(6,*) 'read',infil11     
          open(unit=11,file=infil11,status='old',
     &                           form='unformatted')

          read(11) blk_depth

          write(6,*) 'max depth is ',maxval(maxval(blk_depth(:,:),1))
c          write(6,*) 'depth is ',blk_depth(1,:)

          close(11)

!!--------------------------------------------------------------------------
!! OPEN 3D mean T and S FILES            
!!--------------------------------------------------------------------------

         infil10='zSavg/Savg_'//intid//'_'//runid//
     &'_blk_'
     &                //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

          write(6,*) 'read',infil10
          open(unit=10,file=infil10,status='old',
     &                           form='unformatted')
          do k=1,kdm
            read(10) blk_Savg(:,:,k)
c            write(6,*) k,blk_Savg(100,100,k)
          enddo   

          close(10)

         infil10='zTavg/Tavg_'//intid//'_'//runid//
     &'_blk_'
     &                //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

          write(6,*) 'read',infil10
          open(unit=10,file=infil10,status='old',
     &                           form='unformatted')
          do k=1,kdm
            read(10) blk_Tavg(:,:,k)
!            write(6,*) k,blk_Tavg(33,-2,k),blk_Savg(33,-2,k)
            write(6,*) k,blk_Tavg(45,-2,k),blk_Savg(45,-2,k)
          enddo    

          close(10)

!!--------------------------------------------------------------------------
!! OPEN 3D output files         
!!--------------------------------------------------------------------------
! store surface down at interfaces, starting at k+1

       outfile13='zN2avg/N2avg_'//intid//'_'//runid//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        open(unit=113,file=outfile13,status='replace',
     &     form='unformatted',access='sequential')

c        open(unit=129,file=outfile29,status='new',
c     &     form='unformatted',access='sequential')


!!--------------------------------------------------------------------------
!! LOOP within tile and compute N         
!!--------------------------------------------------------------------------
      
c compute local buoyancy frequency ================================
c first compute the local potential density of subsequent pairs
c loop over i,j,k, run EOS twice, for each pair to compute gradient 
c check if land
c N2 computed at faces, surface at k=1 is N2=0 

       do j=1-halo,blkj+halo
         do i=1-halo,blki+halo

           N2(i,j,:) = 99.9
           if (blk_coast(i,j).eq.1) then

           N2(i,j,1) = 0.0    ! surface value
           I0(:)=0

           pref = 0.0
           cc = count(abs(blk_Tavg(i,j,:)).LT.99999.9)

c           if (i.eq.33.and.j.eq.-2) then
c             write(*,*)cc,blk_Tavg(i,j,:),blk_Savg(i,j,:),
c           endif

           l=0
           m=0
           do k=1,cc-1

             ! at first interface
             pref = dz*rhoref*grav + pref
!             pref = 2000.d4 !no change in pressure!!!

             do kk=1,2
               temp  = blk_Tavg(i,j,k+kk-1)
               salin = blk_Savg(i,j,k+kk-1)

      rpdb=pref*prs2pdb                !ref. pressure in dbar
      c101=c001+(c018-c021*rpdb)*rpdb  !num. constant    coefficent
      c103=c003+(c019-c022*rpdb)*rpdb  !num.    T^2      coefficent
      c105=c005+c020*rpdb              !num.       S     coefficent
      c108=c008+c023*rpdb              !den. constant    coefficent
      c109=c009-c025*rpdb**3           !den.    T        coefficent
      c111=c011-c024*rpdb**2           !den.    T^3      coefficent
      c132=c032+(c042-c043*rpdb)*rpdb  !num.    T        coefficent
      c135=c035-c045*rpdb**3           !den. constant    coefficent
      c137=c037-c044*rpdb**2           !den.    T^2      coefficent
      c151=c051+c059*rpdb

      sig_n=c101+temp*(c002+temp*(c103+temp*c004))+
     &                    salin*(c105-temp*c006+salin*c007)
      sig_d=c108+temp*(c109+temp*(c010+temp*
     &                   (c111+temp*c012))) +
     &                    salin*(c013-temp*(c014+temp
     &                   *temp*c015) +
     &  sqrt(max(sqrmin,salin))*(c016+temp*temp*c017))
      sig_q = aone/sig_d
      rho(kk) = sig_n*sig_q
             enddo !kk

             ! compute N
             N2(i,j,k+1) = grav/rhoref*(rho(2)-rho(1))/dz

             ! where there is an instability
             ! set N2 to small value
c             if (N2(i,j,k+1).LT.0.0) N2(i,j,k+1) = 1.0E-10

c             if(N2(i,j,k+1).EQ.0.0)then !default
             if(N2(i,j,k+1).LE.0.0)then
               l=l+1
               I0(l) = k+1
             endif

             if (i.eq.45.and.j.eq.-2) then
               write(*,*)k+1,cc,N2(i,j,k+1),salin,temp
             endif

           enddo !k


! replace bottom and zero N2 values with nearest good value
         dd = count(I0.GT.0)
         if (dd.GT.0) then
           ! loop over dd
           ! first find index differences
           ID(:) = 0
           do k=1,dd-1
             ID(k) = I0(k+1)-I0(k)
           enddo 

           if (i.eq.45.and.j.eq.-2) then
             write(*,*) I0
             write(*,*) ID
           endif

           ! replace singles
           ee = count(ID.EQ.0)
           k=1
           if(ee.EQ.kdm+1)then      !1 I0 value
             if (I0(k).LT.cc)then !first single bad AND ONLY value
               N2(i,j,I0(k)) = (N2(i,j,I0(k)-1)+N2(i,j,I0(k)+1))*0.5
             elseif (I0(k).EQ.cc)then !last single bad AND ONLY value
               N2(i,j,I0(k)) = N2(i,j,I0(k)-1)
             endif
           else                     !more than 1 I0 value
             if (ID(k).GT.1.AND.I0(k).LT.cc)then !first single bad value
                N2(i,j,I0(k)) = (N2(i,j,I0(k)-1)+N2(i,j,I0(k)+1))*0.5
             endif

             ! more than 2 singles 
             do k=2,dd-1
               if (ID(k).GT.1.AND.ID(k-1).GT.1)then  !singles in the middle
                 N2(i,j,I0(k)) = (N2(i,j,I0(k)-1)+N2(i,j,I0(k)+1))*0.5
               endif
             enddo

             k=dd
             if (ID(k-1).GT.1.AND.I0(k).LT.cc)then !last single bad value
               N2(i,j,I0(k)) = (N2(i,j,I0(k)-1)+N2(i,j,I0(k)+1))*0.5
             endif

             k=dd
        if (ID(k-1).GT.1.AND.I0(k).EQ.cc.AND.N2(i,j,I0(k)).LE.0.0)then !single at 1 above the bottom
               N2(i,j,I0(k)) = N2(i,j,I0(k)-1)
             endif
           endif !ee

           ! then take care of the 1s
           k=0
           do while (k.LT.dd)
             k=k+1 
             if (ID(k).EQ.1)then !find first non 1
               l=1
               do while (ID(k+l).EQ.1) 
                 l=l+1
               enddo    
               write(*,*)k,l,I0(k),I0(k+l)
               write(*,*)N2(i,j,I0(k)-1),N2(i,j,I0(k+l)+1)
               ! in case the last 1 is not at the bottom
               ! I0(k) = cc
               if (I0(k+l).EQ.cc) then ! boundary values of 1s
!                 N2(i,j,I0(k):I0(k+l)) = N2(i,j,I0(k)-1) 

!                 or take mean value; NOTE, this only works if N2(i,j,I0(k)-2 has been filled  
!                 this to reduce effect of small N values for mid-lattitude files
!         N2(i,j,I0(k):I0(k+l)) = 0.5*(N2(i,j,I0(k)-1)+N2(i,j,I0(k)-2)) 

!                 This skips one positive bad value that is the result of the interpolation
                  N2(i,j,I0(k)-1:I0(k+l)) = N2(i,j,I0(k)-2) 
                  write(*,*)'g1',I0(k),I0(k+l)
               else 
                 N2(i,j,I0(k):I0(k+l)) = 
     &                 0.5*(N2(i,j,I0(k)-1)+N2(i,j,I0(k+l)+1)) 
         write(*,*)'g2',I0(k),I0(k+l)
               endif  
               k = k+l-1 !skip ahead ...   
             endif 
           enddo

         endif !dd

         ! last face value at cc+1
         N2(i,j,cc+1) = N2(i,j,cc) 

         if (i.eq.45.and.j.eq.-2) then
           do k=1,kdm+1       
               write(*,*)k,cc,N2(i,j,k),blk_depth(i,j)
           enddo
         endif

         dd = count(isnan(N2(i,j,:)))
         if (dd.GT.0) write(*,*)'NaN in ',i,j,blk_depth(i,j)
 
         dd = count(N2(i,j,:).LT.0.0)
         if (dd.GT.0) write(*,*)'negative values in ',i,j,blk_depth(i,j)

         endif   !coast
         enddo   !i
       enddo     !j

! write the data
       write(*,*) 'store the data +.+.+'

       do k=1,kdm+1
          write(113) N2(:,:,k)
       enddo
       close(113)

!! loop over blocks ----------------------------------------------------
      enddo ! iblk
      enddo ! jblk
!! loop over blocks ----------------------------------------------------

      write(*,*) 'main program finished ~~~~~~~~~~~~~~~~~~~~~~~~' 

      stop
      end
