      PROGRAM NEIG_c008 
! NEIG_c004_v1.f, MCB, USM, 2017-11-13
!   eigenmodes for 1/25, depth cutoff is 200m
! NEIG_c008_v3.f, MCB, USM, 2017-9-22
!   import N computed using local reference pressure with mean_N_c008_v1.f 
! NEIG_c008_v2.f, MCB, USM, 2017-5-16
!   compute eigenvalues and vectors, independent of latitude
!   wave numbers and phasespeeds are only affected by that
! NEIG_c008_v1.f, MCB, USM, 2017-4-29
!   compute buoyancy frequency and extract eigen values and modes


      IMPLICIT NONE
      
      integer idm,jdm,kdm,blki,blkj ! block sizes
      integer maxobs                ! tend-tstart+1
      integer maxobstot             ! number of observations
      integer tstart                ! start time
      integer tend                  ! end time
      integer mbdy                  ! size of halo
      real grav                     ! gravitational acceleration
      real RHO_0                    ! reference density
      real flag                     ! missing value flag
      integer jblks,jblke !1,5
      integer iblks,iblke !1,10
      character runid*3,intid*2
      character timeid*1
      integer MEIG,MEIGout,IU,IL,INFO
      real VL, VU, ABSTOL
      real dz,omf,omm2,pi,orth1,orth2
     
! expt_05.1
c      parameter(idm=4500,jdm=3528,kdm=240)
c      parameter(blki=150,blkj=196,mbdy=3) 
c glbc004
      parameter(idm=9000,jdm=7055,kdm=240) !150 x 60 = 9000 and 200 x 35 = 7000
      parameter (blki=150,blkj=200,mbdy=3)
      
      parameter(grav=9.806)
      parameter(RHO_0=1034.0)
      parameter(flag=2.0**100)
      parameter(dz=25.0)
      parameter(omf=7.292115900231276/100000.0)
      parameter(omm2=1.405189025081018/10000.0)
      parameter(pi=3.141592653589793)
      parameter(IL=1)    ! lower and upper mode number
      parameter(IU=5)


! filter parameters
      character infil10*240,infil12*240
      character outfile12*60,outfile13*60,outfile14*60
      character outfile15*60,outfile16*60
      character*2 numstr(90)
      character*2 constr(8)
      character*3 cnt
      character(len=8) :: fmt

      integer i,j,k,l,m,iblk,jblk,joe,joe2,JJ,ll,numt,numloop,cc,dd
      integer ii
      real fcor 

c 1D variables
      real zc, zf, N2,WVAL,WORK
      dimension zc(1:kdm),zf(1:kdm-1)      
      dimension N2(1:kdm-1),WVAL(1:kdm) 
      dimension WORK(1:35*kdm) 
      
      INTEGER I0, ID, IWORK,IFAIL
      dimension I0(1:kdm-1), ID(1:kdm-1),IWORK(1:5*kdm) 
      dimension IFAIL(1:kdm-1)


c 2D variables
      real A,B
      dimension A(1:kdm,1:kdm),B(1:kdm,1:kdm)

      integer blk_coast,blk_boti
      dimension blk_coast(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_boti(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:2)  !location of bottom matrix

      real blk_plat,blk_depth,blk_ortho
      dimension blk_plat(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_depth(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)
      dimension blk_ortho(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy)

c 3D variables
      real  blk_mean_N2
      dimension blk_mean_N2(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,kdm+1) 

c     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH

      real, allocatable :: Zvec(:,:),blk_Evec(:,:,:,:)
      real, allocatable :: blk_EVAL(:,:,:),LW(:)

      data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/



      MEIG = IU-IL+1 

      allocate(Zvec(1:kdm,1:MEIG))
      allocate(blk_Evec(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:kdm,1:MEIG))
      allocate(blk_EVAL(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:MEIG))
      allocate(LW(1:MEIG))

c z grid, not used in this code .....
      zc(1)=dz/2.0
      do k=2,kdm
        zc(k) = zc(k-1) + dz
        zf(k-1) = (k-1)*dz
c        write(*,*)k,zc(k),zf(k-1)
      enddo

c ------------------------------------------------------------
! read begin and end of blocks
      read(*,*) jblks
      read(*,*) jblke
      read(*,*) iblks
      read(*,*) iblke
      read(*,*) runid
      read(*,*) intid
!DO_FILT, maxobs, FC1, FC2, ftnm, tidcon; make sure numcon

      write (6,'(a,a3)') ' runid: ',runid
      call flush(6)

      write (6,'(a,a3)') ' intid: ',intid
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

      ! initialize field
      blk_EVAL(:,:,:) = 0.0  
      blk_Evec(:,:,:,:) = 0.0  
      blk_ortho(:,:) = -99.9  
      blk_boti(:,:,:) = -99  
          
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

! lat and lon
         infil10='griddata/plat_'//runid//
     &'_blk_'
     &                //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

          write(6,*) 'read',infil10
          open(unit=10,file=infil10,status='old',
     &                           form='unformatted')

          read(10) blk_plat

          close(10)

! depth
         infil10='griddata/depth_'//runid//
     &'_blk_'
     &                //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

          write(6,*) 'read',infil10
          open(unit=10,file=infil10,status='old',
     &                           form='unformatted')

          read(10) blk_depth

          close(10)


!!--------------------------------------------------------------------------
!! OPEN 3D mean density          
!!--------------------------------------------------------------------------
          infil12='zN2avg/'
     &  //'N2avg_'//intid//'_'//runid//
     &'_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        write(6,*) 'read',infil12
        open(unit=12,file=infil12,status='old',form='unformatted')


!!--------------------------------------------------------------------------
!! OPEN output files         
!!--------------------------------------------------------------------------

       outfile13='modes/'
     &  //'EVAL_'//intid//'_'//runid//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        open(unit=113,file=outfile13,status='replace',
     &     form='unformatted',access='sequential')


       outfile14='modes/'
     &  //'EVEC_'//intid//'_'//runid//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        open(unit=114,file=outfile14,status='replace',
     &     form='unformatted',access='sequential')



       outfile15='modes/'
     &  //'ORTHO_'//intid//'_'//runid//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        open(unit=115,file=outfile15,status='replace',
     &     form='unformatted',access='sequential')


       outfile16='modes/'
     &  //'BOTI_'//intid//'_'//runid//'_'
     &  //numstr(jblk)//'_'//numstr(iblk)//'.BinF'

        open(unit=116,file=outfile16,status='replace',
     &     form='unformatted',access='sequential')



!!--------------------------------------------------------------------------
!! Load variables 
!!--------------------------------------------------------------------------

      ! initialize values
      blk_mean_N2(:,:,:) = 0.0

      ! read some files
      do k=1,kdm+1 !depth 1

        read(12) blk_mean_N2(:,:,k)	     
  
      enddo !depth 1

!!--------------------------------------------------------------------------
! loop over x,y, compute N profile and eigen modes
! for the non-land values
!!--------------------------------------------------------------------------

      do j=1-mbdy,blkj+mbdy
        do i=1-mbdy,blki+mbdy

c          if (j.eq.38.and.i.eq.29)then !test

          if(blk_coast(i,j).EQ.1.AND.blk_depth(i,j).GT.200) then !ocean

            ! detect number of cell centers above blk_depth(i,j)
            dd = count(blk_depth(i,j)-zc.GT.0.0)

c            do k=1,kdm 
c              write(*,*)k,blk_depth(i,j),zc(k),dd
c            enddo   
 
            ! find first good cell above bottom (cell center)
c            cc = kdm - count(blk_mean_N2(i,j,:).EQ.flag)
            cc = count(abs(blk_mean_N2(i,j,:)).LT.99.9)-1

            if (i.eq.33.and.j.eq.-2) then
            do k=1,kdm+1
              write(*,*)k,blk_mean_N2(i,j,k),(k-1)*dz,blk_depth(i,j)
            enddo
            endif


            ! store number of layers, cc is the most relevant
            blk_boti(i,j,1) = cc
            blk_boti(i,j,2) = dd

            ! check if dd differs with cc by more than 1
            ! if more than 1, skip! 
            ! else do eigenvalue calculation   
            if(abs(dd-cc).GT.1) then 
              write(*,*)'error depth',i,j,dd,cc,blk_depth(i,j) 
            else
              
c            ! find first bad cell
c            cc = minloc(abs(blk_mean_N2(i,j,:)), 1, 
c     &                mask=abs(blk_mean_N2(i,j,:)).EQ.flag)

            ! print cc is num layers cc+1 is num faces
            ! cc+2 should always have 99.9 
       write(*,*)'N2 vals',blk_mean_N2(i,j,cc+2),blk_mean_N2(i,j,cc+1),
     &             blk_mean_N2(i,j,cc)

! -------------------------------------------------------
! compute eigenvalues and eigenvectors
! A*phi - lambda*B*phi = 0 
! phi    = eigenvector blk_Evec
! lambda = eigenvalue blk_EVAL
! k = sqrt(blk_EVAL*(om^2-f^2)) 

            ! prepare A and B matrices
            ! A and B are layed out Matlab style
            ! vertical,horizontal
            ! at cc-1 faces (excluding top and bottom)
            ! recall, cc is then number of layers

            ! N2 is at faces excluding the surface and bottom
            ! NOTE N2 should be bottom up
            N2(1:cc-1)   = blk_mean_N2(i,j,cc:2:-1) !bottom up
            N2(cc:kdm-1) = 99.9          !filler

            A(:,:)=0.0
            B(:,:)=0.0
            ii=1
            A(ii,ii)   = -2.*1./dz**2
            A(ii,ii+1) =  1.*1./dz**2
            do ii=2,cc-2 !loop over number of layers
              A(ii,-1+ii) =  1.*1./dz**2
              A(ii, 0+ii) = -2.*1./dz**2
              A(ii, 1+ii) =  1.*1./dz**2    
            enddo
            ! index is already advanced to ii+i!! 
            A(ii,-1+ii) =  1.*1./dz**2
            A(ii, 0+ii) = -2.*1./dz**2
 
            ! B = diag(-N2(2:end-1))/(om^2-f^2);  % N to power 2
            ! over faces, excluding bottom and surface
            do ii=1,cc-1
              B(ii,ii) = -N2(ii)
            enddo

            ! solve for phase speed and wave length
            ! http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssygvx.f
            ! SUBROUTINE SSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
            !                    VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
            !                    LWORK, IWORK, IFAIL, INFO )
            WVAL(:)   = 1.0 
            Zvec(:,:) = 0.0
            IWORK(:)  = 0
            IFAIL(:)  = 0
            ABSTOL = 2*SLAMCH('S')
c            write(*,*) ABSTOL,VL, VU, IL, IU,MEIG,cc
!            ABSTOL = 1e-9

            ! note A and B are multiplied by -1 as B has to be positive
            call SSYGVX( 1, 'V', 'I', 'U',cc-1,
     &      -A(1:cc-1,1:cc-1), cc-1,
     &      -B(1:cc-1,1:cc-1), cc-1,                        
     &      VL, VU, IL, IU, ABSTOL, MEIGout, WVAL(1:cc-1),
     &      Zvec(1:cc-1,1:MEIG), cc-1, WORK(1:35*(cc-1)),35*(cc-1), 
     &      IWORK(1:5*(cc-1)), IFAIL(1:cc-1), INFO )    
        
c            if(j.eq.100.and.i.eq.100)then
c              do k=1,kdm-1
c           write(*,*)Zvec(k,1),Zvec(k,2),Zvec(k,3),Zvec(k,4),Zvec(k,5)
c              enddo
c              stop
c            endif  

            ! check orthogonality
            fcor = 2.0*omf*sin( blk_plat(i,j)*pi/180.0 )
            orth1 = 0.0 
            orth2 = 0.0 
            do ii=1,cc-1
c          orth1 = orth1 + N2(ii)/(omm2**2-fcor**2)*Zvec(ii,1)*Zvec(ii,2)
c          orth2 = orth2 + N2(ii)/(omm2**2-fcor**2)*Zvec(ii,1)*Zvec(ii,1)
          orth1 = orth1 + N2(ii)*Zvec(ii,1)*Zvec(ii,2)
          orth2 = orth2 + N2(ii)*Zvec(ii,1)*Zvec(ii,1)
            enddo
            write(*,*)'orthogo: ',orth1,orth2 
            blk_ortho(i,j) = orth2 

            ! store the eigenvectors
            ! surface down, which is the HYCOM default
            do l=1,MEIG 
              blk_Evec(i,j,1:cc-1,l) = Zvec(cc-1:1:-1,l) 
            enddo 

            ! store eigenvalues blk_EVAL = k^2/(om^2-f^2)
            blk_EVAL(i,j,:) = WVAL(1:MEIG)

            ! compute wave lengths; not stored  
            LW   = 2.0*pi/sqrt( WVAL(1:MEIG)*(omm2**2.0-fcor**2.0) )

            write(*,*)'L=',LW/1000.
c            write(*,*)'work=',WORK(1:10)
            write(*,*)'info=',INFO,WORK(1),WORK(1)/(cc-1),MEIGout
c            write(*,*)'IFAIL=',IFAIL

            write(*,*)'    '

          endif !dd-cc=0
          endif !ocean
c          endif !test

        enddo !i
      enddo   !j

! -----------------------------------------------
! store output
! write each mode into one record
! -----------------------------------------------
      do l=1,MEIG

        write(113) blk_EVAL(:,:,l)
        write(114) blk_Evec(:,:,:,l)

      enddo 

      write(115) blk_ortho(:,:)
      write(116) blk_boti(:,:,:) 

      ! close the files
      close(113)
      close(114)
      close(115)
      close(116)

!! loop over blocks ----------------------------------------------------
      enddo ! iblk
      enddo ! jblk
!! loop over blocks ----------------------------------------------------

      write(6,*) 'Finished...'

      stop

      end program NEIG_c008 

