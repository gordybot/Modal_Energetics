      PROGRAM PSURF
       
      !!! Compute the surface pressure by loading the eigenvectors W,
      !!! renormalizing dW/dz, and etracting the surface value of U.
      !!! Load the pressure mode amplitudes and multiply by the surface
      !!! value to get Psurf. Then divide by (rho0 g) to get surface
      !!!  displacement.

      IMPLICIT NONE

      ! Declare variables.
      ! Declare and type variables.
      integer i,j,k,l,q,cc, mbdy
      integer idm, jdm, kdm, blki, blkj
      integer iblk, iblks, iblke
      integer jblk, jblks, jblke
      integer IL, IU, MEIG
      character*3 runid
      character infileE*240, infileB*240
      character infileP*240, outfile*240
      character*2 numstr(90)
      integer signUvec, offset_i, offset_j
      real dz, grav, sum1
      real mean_dWdz, var_dWdz, meanA
      real mean_Pcell, H, tempSum1, tempSum2

      ! Declare array variables.
      integer, allocatable, dimension(:,:,:)  :: blk_boti
      real, allocatable, dimension(:) :: W
      real, allocatable, dimension(:) :: dWdz
      real, allocatable, dimension(:,:,:,:) :: blk_Evec
      real, allocatable, dimension(:) :: Uvec
      real, allocatable, dimension(:,:,:) :: esurf
      real, allocatable, dimension(:,:,:,:) :: blk_modeP

      real, allocatable, dimension(:,:) :: A
      real, allocatable, dimension(:,:,:,:) :: glb_modeP
      real, allocatable, dimension(:,:,:,:) :: glb_psurf

      ! Assign parameters.
       parameter(idm=9000,jdm=7055,kdm=240)
       parameter(blki=150,blkj=200,mbdy=3)
       parameter(IL=1,IU=5)
       parameter(dz=25.0)
       parameter(grav=9.81)

       ! Define some basic variables.
       MEIG=IU-IL+1
       data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/

      ! Allocate variables.
      allocate( blk_boti(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy, 1:2) )
      allocate( W(0:kdm) )
      allocate( dWdz(1:kdm) )
      allocate( Uvec(1:kdm) )
      allocate( A(1:kdm, 1:MEIG) )
      allocate( blk_Evec(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy, 1:kdm,
     &           1:MEIG) )
      allocate( esurf(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:MEIG) )
      allocate( glb_psurf( 1:idm, 1:jdm, 1:2, 1:MEIG) )
      !allocate( glb_modeP( 1:idm, 1:jdm, 1:2, 1:MEIG) )
      allocate( blk_modeP( 1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy, 1:2, 
     &           1:MEIG) )

      ! Read from input file.
       read(*,*) jblks
       read(*,*) jblke
       read(*,*) iblks
       read(*,*) iblke
       read(*,*) runid

       write(6,'(a,a3)'),'runid: ',runid
       call flush(6)

       write(6,'(a,4i3)') 'jblks jblke iblks iblke: ',jblks, jblke,
     & iblks, iblke
       call flush(6)

       ! Initialize global fields to 0.0
!       glb_modeP(:,:,:,:) = 0.0
       glb_psurf(:,:,:,:) = 0.0

       ! Read global pressure mode amplitudes.
!       infileP = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
!     &     //'modes_global/glb_phfvm_may.BinF'
!       write(6,*) 'Reading global pressure modes from ',infileP,'.'
!       call flush(6)
!       open(unit=10,file=infileP,status='old',form='unformatted')
!       do l=1,MEIG
!         read(10) glb_modeP(:,:,:,l)
!       enddo !loop over l
!       close(10)
       

      ! Loop over iblk, jblk
       do jblk=jblks,jblke
       do iblk=iblks,iblke
         write(6,*) 'iblk, jblk = ',iblk,jblk
         call flush(6)

         offset_i = blki*(iblk-1)
         offset_j = blkj*(jblk-1)       
 
       ! Initialize blk variables:
         blk_modeP(:,:,:,:)=0.0
         !Initialize blk_Evec
         blk_Evec(:,:,:,:)=0.0
         ! Initialize blk_boti
         blk_boti(:,:,:)=0 ! integer.
         ! Reset ebot to 0.0
         esurf(:,:,:) = 0.0
         !modeP, pbot

      ! READ a bunch of FILES.
      ! Load blk_boti, blk_evec
         ! Read in the eigenvectors.       
         infileE='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/modes/'
     &   //'EVEC_1z_'//runid//'_'//numstr(jblk)//'_'//numstr(iblk)//
     &   '.BinF'
         write(6,*) 'reading ',infileE
         call flush(6)
         open(unit=11,file=infileE,status='old',form='unformatted')
         !Read blk_Evec in same way it was written out.
         do l=1,MEIG
           read(11) blk_Evec(:,:,:,l)
         enddo !loop over l
         close(11)

         ! Read in the bottom bin locations.
         infileB='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/modes/'
     &   //'BOTI_1z_'//runid//'_'//numstr(jblk)//'_'//numstr(iblk)//
     &   '.BinF'
         write(6,*) 'reading ',infileB
         call flush(6)
         open(unit=12,file=infileB,status='old',form='unformatted')
         ! Read blk_boti
         read(12) blk_boti(:,:,:)
         close(12)
         write(6,*) 'Finished reading ',infileB
         call flush(6)

         infileP='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &     //'modes/phfvm_221_blk_'//numstr(jblk)//'_'//numstr(iblk)
     &     //'_apr.BinF'
         write(6,*) 'Reading blk pressure modes from ',infileP,'.'
         call flush(6)
         open(unit=13,file=infileP,status='old',form='unformatted')
         do l=1,MEIG
           read(13) blk_modeP(:,:,:,l)
         enddo !loop over l
         close(13)

         ! Loop over i,j.
         do i=1-mbdy,blki+mbdy
         do j=1-mbdy,blkj+mbdy
           ! Reinitialize A to 0
           A(:,:) = 0.0

           cc=blk_boti(i,j,1) ! number of deepest good bin.
           ! Check to see that there ARE any good cells:
           if (cc.GT.1) then

           ! Normalize the eigenvectors.
           do l=1,MEIG
            ! Initialize W, dWdz.
            W(:)=0.0
            dWdz(:) = 0.0
            ! Define dW/dz ~ U,V,P eigenvector.
            do k=1,cc
              W(k) = blk_Evec(i,j,k,l)  ! face-centered
              dWdz(k) = ( W(k)-W(k-1) ) / dz ! cell-centered.
            enddo
            ! We will normalize so that the integral of Uvec**2 dz == H.
            sum1 = sum( dWdz(1:cc)**2 ) * 1.0/cc

            ! Define  Uvec, to be the normalized eigenvector for U
            Uvec(:) = 0.0
            ! Possibly prevent some Infinity values.
            if (sum1.GT.0.0) then
              do k=1,cc
                Uvec(k) = dWdz(k) / (sum1**0.5)
              enddo !k
            endif

            ! Define so vector is positive at the seabed.
            signUvec = 1
            if (Uvec(cc).LT.0.0) then
              signUvec = -1
            endif

            A(:,l) = signUvec * Uvec(:)
            ! Make sure mean(A(:,l))==0.0
            meanA = sum(A(:,l))*1.0/cc
            do k=1,cc
              A(k,l) = A(k,l) - meanA
            enddo !k

            ! Assign top value of eigenvector to esurf.
            esurf(i,j,l) = A(1,l)
          enddo !l loop over modes.
 
          ! Continue check on any good cells.
          else if (cc.LT.2) then
            esurf(i,j,:) = 0.0
          endif !  Check cc.GT.1, at least one good cell.

          ! Store surface pressure in global field.
          do l=1,MEIG
            glb_psurf(offset_i + i, offset_j + j, 1, l) = 
     &        blk_modeP( i, j, 1, l) * esurf(i,j,l)
            glb_psurf(offset_i + i, offset_j + j, 2, l) =
     &        blk_modeP( i, j, 2, l) * esurf(i,j,l)
          enddo !l
         enddo !j
         enddo !i
       enddo !blkj
       enddo !blki

        ! Write out surface pressure mode-by-mode.
        outfile = '/p/work1/grstephe/hycom/GLBc0.03/expt_22.1/'
     &      //'modes_global/glb_psurf_221_mode.BinF'
        write(6,*) 'Writing psurf to ',outfile,'.'
        call flush(6)
        open(unit=97,file=outfile,status='replace',
     &         form='unformatted',access='sequential')
        do l=1,MEIG
          write(97) glb_modeP(:,:,:,l)
        enddo !l
        close(97)

        write(6,*) 'Finished...'
        call flush(6)

        STOP
        END PROGRAM PSURF
