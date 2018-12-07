       PROGRAM COMPUTE_CMN

       !!! Compute the mode-mode conversion. 
       !!! 1/2 * integral( u_m* . del P_n - u_n* . del P_m) dz
       !!! 3.15.18.GS
       IMPLICIT NONE
       
       !!! Declare variables.
       character*3 runid
       character*2 numstr(90)
       character outfile*240, infileB*240, infileX*240, infileY*240
       character infileE*240, infileU*240, infileV*240, infileP*240
       integer i, j, k, l, m, n, mn, cc, cc1, cc2, MEIG, MNdm
       integer blki, blkj, mbdy, kdm, offset_i, offset_j, idm, jdm
       integer iblk, iblks, iblke, jblk, jblks, jblke
       real dz
       complex umk_conj, unk_conj, vmk_conj, vnk_conj
       complex dPmk_dx, dPmk_dy, dPnk_dx, dPnk_dy
       real meanA, sum1, signUvec

       !!! Temporary holding variables.
       complex, allocatable, dimension(:) :: Cmn

       !!! Block variables to load.
       real, allocatable, dimension(:,:) :: pscx
       real, allocatable, dimension(:,:) :: pscy
       integer, allocatable, dimension(:,:,:) :: blk_boti
       real, allocatable, dimension(:,:,:,:) :: blk_Evec
       real, allocatable, dimension(:,:,:,:) :: modeU
       real, allocatable, dimension(:,:,:,:) :: modeV
       real, allocatable, dimension(:,:,:,:) :: modeP
       real, allocatable, dimension(:,:,:,:) :: Aij
       real, allocatable, dimension(:,:) :: A
       real, allocatable, dimension(:) :: Uvec
       real, allocatable, dimension(:) :: W
       real, allocatable, dimension(:) :: dWdz
       real, allocatable, dimension(:,:,:,:) :: blk_ReCmn
       real, allocatable, dimension(:,:,:,:) :: glb_ReCmn

       !!! Output variables.
       complex, allocatable, dimension(:,:,:) :: blk_Cmn       

       !!! Define some parameters.
       parameter(dz=25.0, MEIG=5, MNdm=10)
       parameter(idm=9000, jdm=7055)
       parameter(blki=150, blkj=200, kdm=240, mbdy=3)
       parameter(jblks=1, jblke=35, iblks=1, iblke=60)
       parameter(runid='221')
       data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/

       !!! Allocate variables.
       allocate( Cmn( 1:MNdm ) )
       allocate( pscx( 1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy ) )
       allocate( pscy( 1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy ) )
       allocate( blk_boti(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:2) )
       allocate( blk_Evec(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,
     &                    1:kdm, 1:MEIG ) )
       allocate( modeU(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,1:MEIG,1:2) )
       allocate( modeV(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,1:MEIG,1:2) )
       allocate( modeP(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,1:MEIG,1:2) )
       allocate( Aij(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,1:kdm,1:MEIG) )
       allocate( A(1:kdm, 1:MEIG) )
       allocate( W(0:kdm) )
       allocate( dWdz(1:kdm) )
       allocate( Uvec(1:kdm) )
       allocate( blk_Cmn(1:blki, 1:blkj, 1:MEIG) )
       allocate( blk_ReCmn(1:blki, 1:blkj, 1:MEIG, 1:2) )
       allocate( glb_ReCmn(1:idm, 1:jdm, 1:MNdm, 1:2) )
       !!! Initialize global variables.
       glb_ReCmn(:,:,:,:) = 0.0

       do jblk=jblks,jblke
         do iblk=iblks,iblke
           !!! 
           offset_i = blki*(iblk-1)
           offset_j = blkj*(jblk-1)

           !!! Initialize block variables.
           blk_Cmn(:,:,:) = 0.0
           blk_ReCmn(:,:,:,:) = 0.0
           blk_boti(:,:,:) = 0
           blk_Evec(:,:,:,:) = 0.0
           pscx(:,:) = 0.0
           pscy(:,:) = 0.0
           modeU(:,:,:,:) = 0.0
           modeV(:,:,:,:) = 0.0
           modeP(:,:,:,:) = 0.0

           !!! Load eigenvector tile -- modeP, modeU, modeV.
           !!!  Also pscx, pscy, blk_boti.
           infileB = '/p/work1/mbui/hycom/GLBc0.04/expt_22.1/modes/'
     &            //'BOTI_1z_221_'//numstr(jblk)//'_'//numstr(iblk)
     &            //'.BinF'
           write(6,*) 'reading ',infileB
           call flush(6)
           open(unit=10,file=infileB,status='old',form='unformatted')
           ! Read blk_boti
           read(10) blk_boti(:,:,:)
           close(10)
           write(6,*) 'Finished reading ',infileB
           call flush(6)

           infileE = '/p/work1/mbui/hycom/GLBc0.04/expt_22.1/modes/'
     &            //'EVEC_1z_221_'//numstr(jblk)//'_'//numstr(iblk)
     &            //'.BinF'
           write(6,*) 'reading ',infileE
           call flush(6)
           open(unit=11,file=infileE,status='old',form='unformatted')
           !Read blk_Evec in same way it was written out.
           do l=1,MEIG
             read(11) blk_Evec(:,:,:,l)
           enddo !loop over l
           close(11)

           infileX= '/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &            //'griddata/pscx_221_blk_'//numstr(jblk)//'_'//
     &             numstr(iblk)//'.BinF'
           open(unit=12, file=infileX, form='unformatted',
     &             status='OLD', access='sequential',action='READ')
           write(6,*) 'Reading ',infileX
           call flush(6)
           read(12) pscx(:,:)
           close(12)

           infileY= '/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &             //'griddata/pscy_221_blk_'//numstr(jblk)//'_'//
     &             numstr(iblk)//'.BinF'
           open(unit=13, file=infileY, form='unformatted',
     &              status='OLD', access='sequential',action='READ')
           write(6,*) 'Reading ',infileY
           call flush(6)
           read(13) pscy(:,:)
           close(13)

           infileU = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &            //'modes/uhfvm_221_blk_'//numstr(jblk)
     &            //'_'//numstr(iblk)//'_may.BinF'
           ! Changed from _feb8.BinF to _may.BinF // May 7.

           write(6,*) 'reading ',infileU
           call flush(6)
           modeU(:,:,:,:) = 0.0
           open(unit=14, file=infileU,status='old',form='unformatted')
           read(14) modeU(:,:,:,:)
           close(14)

           infileV = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &            //'modes/vhfvm_221_blk_'//numstr(jblk)
     &            //'_'//numstr(iblk)//'_may.BinF'
           write(6,*) 'reading ',infileV
           call flush(6)
           ! Initialize blk_mHFv. Indices are (i,j,l,[mean amp pha])
           modeV(:,:,:,:) = 0.0
           open(unit=15, file=infileV,status='old',form='unformatted')
           read(15) modeV(:,:,:,:)
           close(15)

           infileP = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &            //'modes/phfvm_221_blk_'//numstr(jblk)
     &            //'_'//numstr(iblk)//'_apr.BinF'
           write(6,*) 'reading ',infileP
           call flush(6)
           ! Initialize blk_mHFp. Indices are (i,j,l,[mean amp pha])
           modeP(:,:,:,:) = 0.0
           open(unit=16, file=infileP,status='old',form='unformatted')
           read(16) modeP(:,:,:,:)
           close(16)

           !!! Compute normalized dW/dz for all i,j into A. 
           !!! Store in Aij so I can take horizontal derivatives.
           Aij(:,:,:,:) = 0.0
           do i=1-mbdy,blki+mbdy
             do j=1-mbdy,blkj+mbdy
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
                 ! We will normalize so that the 
                 ! integral of Uvec**2 dz == H.
                 sum1 = sum( dWdz(1:cc)**2 ) * 1.0/cc

                 ! Define  Uvec, to be the normalized eigenvector for U
                 Uvec(:) = 0.0
                 ! Possibly prevent some Infinity values.
                 if (sum1.GT.0.0) then
                   do k=1,cc
                     Uvec(k) = dWdz(k) / (sum1**0.5)
                   enddo !k
                 endif

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
               enddo !l
               Aij(i,j,:,:) = A(:,:)
               endif ! cc.GT.1
             enddo !i
           enddo !j

           do i=1,blki
             do j=1,blkj
               cc1 = min( blk_boti(i-1, j, 1), blk_boti(i, j-1, 1) )
               cc2 = min( blk_boti(i+1, j, 1), blk_boti(i, j+1, 1) )
               cc = min( cc1, cc2 )
               if (cc.GT.1) then

               do k=1,cc
                 mn = 0
                 Cmn( 1:10 ) = 0.0
                 do m=1,4
                   umk_conj = cmplx( modeU(i,j,m,1), -modeU(i,j,m,2) )
     &                        * Aij(i,j,k,m)
                   vmk_conj = cmplx( modeU(i,j,m,1), -modeU(i,j,m,2) )
     &                        * Aij(i,j,k,m)
                  dPmk_dx=( cmplx( modeP(i+1,j,m,1), modeP(i+1,j,m,2) )
     &                       * Aij(i-1,j,k,m)
     &                    - cmplx( modeP(i-1,j,m,1), modeP(i-1,j,m,2) )
     &                       * Aij(i-1,j,k,m) )
     &                     * 1.0 / ( 2.0*pscx(i,j) )
                  dPmk_dy=( cmplx( modeP(i,j+1,m,1), modeP(i,j+1,m,2) )
     &                       * Aij(i,j+1,k,m)
     &                    - cmplx( modeP(i,j-1,m,1), modeP(i,j-1,m,2) )
     &                       * Aij(i,j-1,k,m) )
     &                     * 1.0 / ( 2.0*pscy(i,j) )

                   do n=m+1,5
                     unk_conj= cmplx( modeU(i,j,n,1), -modeU(i,j,n,2) )
     &                                *Aij(i,j,k,n)
                     vnk_conj= cmplx( modeV(i,j,n,1), -modeV(i,j,n,2) )
     &                                *Aij(i,j,k,n) 
                  dPnk_dx=( cmplx( modeP(i+1,j,n,1), modeP(i+1,j,n,2) )
     &                       * Aij(i-1,j,k,n)
     &                    - cmplx( modeP(i-1,j,n,1), modeP(i-1,j,n,2) )
     &                       * Aij(i-1,j,k,n) )
     &                     * 1.0 / ( 2.0*pscx(i,j) )
                  dPnk_dy=( cmplx( modeP(i,j+1,n,1), modeP(i,j+1,n,2) )
     &                       * Aij(i,j+1,k,n) 
     &                    - cmplx( modeP(i,j-1,n,1), modeP(i,j-1,n,2) )
     &                       * Aij(i,j-1,k,n) )
     &                     * 1.0 / ( 2.0*pscy(i,j) )
                  !!! Compute mode-mode conversion - store in Cmn(1:10)
                     mn = mn + 1
                     Cmn( mn ) = Cmn( mn ) + 0.5 * dz * 
     &                       ( umk_conj * dPnk_dx + vmk_conj * dPnk_dy
     &                       - unk_conj * dPmk_dx - vnk_conj * dPmk_dy)
                   enddo !n 
                 enddo !m 
               enddo !k
               !!! Store Cmn( 1:10 ) into blk_Cmn.
               blk_Cmn( i, j, 1:MNdm ) = Cmn
               !!! Store block values in global array.           
               glb_ReCmn(offset_i + i, offset_j + j, 1:10, 1) = 
     &                  real(real( Cmn ) )
               glb_ReCmn(offset_i+i, offset_j+j, 1:10, 2) = 
     &                  real( aimag( Cmn ) )
               endif !cc.GT.
             enddo !j
           enddo !i
           
           blk_ReCmn(:,:,:,1) = real( real( blk_Cmn(:,:,:) ) )
           blk_ReCmn(:,:,:,2) = real( aimag( blk_Cmn(:,:,:) ) )

           !!! Save a tile.
           outfile ='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/modes/'
     &       //'blk_Cmn_'//numstr(jblk)//'_'//numstr(iblk)//'_may.BinF'
           write(6,*) 'Writing blk_ReCmn (i,j,m,1:2) to ',outfile
           call flush(6)
           open(unit=98,file=outfile,status='replace',
     &          form='unformatted',access='sequential')
           write(98) blk_ReCmn(:,:,:,:)
           close(98)

         enddo !iblk
       enddo !jblk

       !!! Write a global file.
       do l=1,MNdm

       outfile = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &         //'modes_global/glb_221_Cmn_MN'//numstr(l)//'_may.BinF'
       write(6,*) 'Writing glb_ReCmn to ', outfile
       call flush(6)
       open( unit=99, file=outfile,status='replace',
     &       form='unformatted',access='sequential')
       write(99) glb_ReCmn(:,:,l,:)
       close(99)
       enddo ! l

       STOP
       END PROGRAM COMPUTE_CMN

