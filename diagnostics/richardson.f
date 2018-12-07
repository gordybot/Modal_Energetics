       PROGRAM GRADRI

       ! Program to compute the gradient Richardson number 
       ! associated with particular wave modes.
       ! I want the maximum within the water column at a given i,j
       ! for each mode -- so the maximum over time and depth.
       ! Obviously, this won't necessarily predict the shear
       ! and mixing perfectly at any one point, 
       ! but if we assume that it reflects propagating wave modes
       ! I think it could be a useful measure.
       ! Started at Palo Alto Library, 28-June-2018.

       IMPLICIT NONE

       ! Declare variables.
       character*2 numstr(90)
       character*3 runid
       character infileE*240, infileB*240, infileU*240
       character infileV*240, infileN*240, outfile*240
       integer iblk, jblk, iblks, iblke, jblks, jblke, blki, blkj
       integer i,j,k,l, MEIG, cc, offset_i, offset_j, idm, jdm, kdm
       integer mbdy, signUvec
       real dz, meanA, sum1

       integer, allocatable, dimension(:,:,:) :: blk_boti
       real, allocatable, dimension(:,:,:,:) :: blk_EVEC
       real, allocatable, dimension(:,:,:) :: blk_Nsqr
       real, allocatable, dimension(:,:,:,:) :: blk_modeU
       real, allocatable, dimension(:,:,:,:) :: blk_modeV
       real, allocatable, dimension(:,:,:) :: glb_min_gradri
       real, allocatable, dimension(:,:,:) :: min_gradri
       real, allocatable, dimension(:,:) :: velo
       real, allocatable, dimension(:,:) :: prof_gradri
       real, allocatable, dimension(:,:) :: shear
       real, allocatable, dimension(:) :: Uvec
       real, allocatable, dimension(:,:) :: A
       real, allocatable, dimension(:) :: W
       real, allocatable, dimension(:) :: dWdz

       ! Define parameters.
       parameter(dz=25.0, MEIG=5)
       parameter(blki=150, blkj=200, kdm=240, mbdy=3)
!       parameter(jblks=10, jblke=10, iblks=10, iblke=10)
!       parameter(runid='221')
       parameter(idm=9000,jdm=7055)

       data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/ 

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

       ! Allocate variables.
       allocate( blk_boti(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,1:2) )
       allocate( blk_EVEC(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,
     &           1:kdm,1:MEIG) )
       allocate( blk_Nsqr(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 0:kdm ) )
       allocate( blk_modeU(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 
     &            1:2, 1:MEIG) )
       allocate( blk_modeV(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,
     &            1:2, 1:MEIG) )
       allocate( glb_min_gradri(1:idm, 1:jdm, 1:MEIG) )
       allocate( min_gradri( 1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,
     &            1:MEIG) )
       allocate( velo( 1:kdm, 1:4) )
       allocate( prof_gradri( 1:kdm, 1:2) )
       allocate( shear( 1:kdm, 1:2) )
       allocate( Uvec(1:kdm) )
       allocate( A(1:kdm, 1:MEIG) )
       allocate( W(0:kdm) )
       allocate( dWdz(1:kdm) )

       ! Initialize global variables.
       glb_min_gradri(:,:,:) = 0.0

       ! Loop over tiles.
       do jblk = jblks, jblke
       do iblk = iblks, iblke
         ! Define offset for transfer from tile to global field.
         offset_i = blki*(iblk-1)
         offset_j = blkj*(jblk-1)

         ! Initialize block variables to 0.
         blk_EVEC(:,:,:,:) = 0.0
         blk_Nsqr(:,:,:) = 0.0
         blk_modeU(:,:,:,:) = 0.0
         blk_modeV(:,:,:,:) = 0.0
         blk_boti(:,:,:) = 0

         ! Load Eigenvectors, Nsquared, Velocity-amplitudes, blk_boti
         ! Read in the eigenvectors.       
         infileE='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/modes/'
     &   //'EVEC_1z_'//runid//'_'//numstr(jblk)//'_'//numstr(iblk)//
     &   '.BinF'
         write(6,*) 'reading ',infileE
         call flush(6)
         open(unit=11,file=infileE,status='old',form='unformatted')
         !Initialize blk_Evec
         blk_EVEC(:,:,:,:)=0.0
         !Read blk_Evec in same way it was written out.
         do l=1,MEIG
           read(11) blk_EVEC(:,:,:,l)
         enddo !loop over l
         close(11)
         write(6,*) 'Finished reading ',infileE
         call flush(6)

         ! Read in the bottom bin locations.
         infileB='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/modes/'
     &   //'BOTI_1z_'//runid//'_'//numstr(jblk)//'_'//numstr(iblk)//
     &   '.BinF'
         write(6,*) 'reading ',infileB
         call flush(6)
         open(unit=12,file=infileB,status='old',form='unformatted')
         ! Initialize blk_boti
         blk_boti(:,:,:)=0
         ! Read blk_boti
         read(12) blk_boti(:,:,:)
         close(12)
         write(6,*) 'Finished reading ',infileB
         call flush(6)

       infileN = '/p/work1/mbui/hycom/GLBc0.04/expt_22.1/zN2avg/'
     &      //'N2avg_1z_221_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
         write(6,*) 'reading ',infileN
         call flush(6)
         open(unit=12,file=infileN,status='old',form='unformatted')
         ! Initialize blk_boti
         blk_Nsqr(:,:,:)=0
         ! Read blk_boti
         read(12) blk_Nsqr(:,:,:)
         close(12)
         write(6,*) 'Finished reading ',infileN
         call flush(6)

         infileU='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &  'modes/uhfvm_'//runid//'_blk_'//numstr(jblk)//
     &  '_'//numstr(iblk)//'_feb8.BinF'
         write(6,*) 'reading ',infileU
         call flush(6)
         open(unit=12,file=infileU,status='old',form='unformatted')
         blk_modeU(:,:,:,:)=0
         read(12) blk_modeU(:,:,:,:)
         close(12)
         write(6,*) 'Finished reading ',infileU
         call flush(6)

         infileV='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &  'modes/vhfvm_'//runid//'_blk_'//numstr(jblk)//
     &  '_'//numstr(iblk)//'_feb8.BinF'
         write(6,*) 'reading ',infileV
         call flush(6)
         open(unit=12,file=infileV,status='old',form='unformatted')
         blk_modeV(:,:,:,:)=0
         read(12) blk_modeV(:,:,:,:)
         close(12)
         write(6,*) 'Finished reading ',infileV
         call flush(6)

         ! For each point in the tile, compute the shear.
         do i=1,blki
         do j=1,blkj
           cc = blk_boti(i,j,1)
           ! Convert EVEC to horizontal-velocity eigenvector, U
           if (cc.GT.1) then
             ! Normalize the eigenvectors.
             do l=1,MEIG
               ! Initialize W, dWdz.
               W(:)=0.0
               dWdz(:) = 0.0
               ! Define dW/dz ~ U,V,P eigenvector.
               do k=1,cc
                 W(k) = blk_EVEC(i,j,k,l)  ! face-centered
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
             enddo ! loop on modes.

             do l=1,MEIG
               ! Reinitialize velocity profile to 0.
               velo(:,:) = 0.0

               ! Compute velocity profile.
               velo(1:cc, 1) = blk_modeU(i,j,1,l) * A(1:cc,l)
               velo(1:cc, 2) = blk_modeU(i,j,2,l) * A(1:cc,l)
               velo(1:cc, 3) = blk_modeV(i,j,1,l) * A(1:cc,l)
               velo(1:cc, 4) = blk_modeV(i,j,2,l) * A(1:cc,l)

               ! Compute shear in an odd way. 1,3 are real components
               ! 2,4 are imaginary components.
               prof_gradri(:,:) = 0.0
               shear(:,:) = 0.0
               do k=2,cc
                 shear(k-1,1) = ( ( velo(k,1) - velo(k-1,1) )**2 + 
     &                         ( velo(k,3) - velo(k-1,3) )**2 ) **1.0
                 shear(k-1,2) = ( ( velo(k,2) - velo(k-1,2) )**2 + 
     &                         ( velo(k,4) - velo(k-1,4) )**2 ) **1.0
                 ! TODO: Interpolate Nsqr in vertical?
                 prof_gradri(k-1,1) = blk_Nsqr(i,j,k-1) / shear(k-1,1)
                 prof_gradri(k-1,2) = blk_Nsqr(i,j,k-1) / shear(k-1,2)
               enddo ! k

              ! Store min value of the gradient Richardson number
               min_gradri(i,j,l) = MINVAL( prof_gradri(1:k-1, 1:2) )

              ! Store in a global field.
            glb_min_gradri(offset_i+i, offset_j+j, l)= min_gradri(i,j,l)
             enddo ! l
           endif ! cc.GT.1, there are some good cells.
         enddo ! j
         enddo ! i
       enddo ! iblk
       enddo ! jblk

       ! Write min / max gradient Richardson number to file.
       outfile = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &         //'modes_global/glb_minRi_221.BinF'
       write(6,*) 'Writing gradRi to ',outfile,'.'
       call flush(6)
       open(unit=99,file=outfile,status='replace',
     &      form='unformatted',access='sequential')
       do l=1,MEIG
         write(99) glb_min_gradri(:,:,l)
       enddo !l

       write(6,*) 'Finished...'
       call flush(6)
     
       STOP
       END PROGRAM GRADRI

