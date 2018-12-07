       program ebot

       ! Do a modal decomposition of u, v velocities and 
       ! (densities converted to pressure perturbations)
       ! using Maarten's eigenvectors.
       ! Normalize U mode structure, as this changes the answer 
       ! by about a factor H (1000 ish)
       ! 8-22-2017 GRS. Updated to include pressure, disinclude
       ! harmonic fits, on 10-03-2017.
       ! 10-25-2017 GRS. Update to not include pressure.
       ! Takes too long and I suspect an error in p-modes.
       !
       ! Read INPUTS:
       ! W-Eigenvectors
       !'/p/work1/mbui/hycom/GLBc0.04/expt_05.1/modes/'
       ! ... EVEC_051_JJ_II.BinF 
       ! New modes have a "1z" in the filename!
       ! 
       ! BOTI
       ! /p/work1/mbui/hycom/GLBc0.04/expt_05.1/modes/'
       ! ... BOIT_051_JJ_II.BinF
       ! Look for "1z", just in case.
       !
       ! U, V, density
       ! /p/work1/grstephe/hycom/GLBc0.04/expt_05.1/'
       !  ...  u_iso/u_051_blk_JJ_II.BinF
       !  ...  v_iso/v_051_blk_JJ_II.BinF
       !  ...  sig/sig_051_blk_JJ_II.BinF
       !
       ! Write OUTPUTS:
       ! '/p/work1/grstephe/hycom/GLBc0.04/expt_05.1/modes/'
       !  ... modeU_051_blk_JJ_II.BinF
       !  ... modeV_051_blk_JJ_II.BinF
       !  ... modeP_051_blk_JJ_II.BinF
       !  These hold modal decomposition amplitudes as timeseries
       !  in an array indexed by (i, j, t, l=mode number), sized
       !  (1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:maxobstot, 1:MEIG+1)
       !  Written as e.g., write() mode(:,:,:,:) (in one write statement)

       implicit none

       !Declare variables' types here:
       character*3 runid
       character infileE*240, infileB*240
       character outfilE*240
       logical isThere, OVERWRITE
       !character infileP*240, outfilP*240
       character*2 numstr(90)
       integer blki, blkj, mbdy
       integer idm, jdm, kdm
       integer IL,IU,MEIG
       integer i,j,k,l
       integer iblk, iblks, iblke
       integer jblk, jblks, jblke
       integer cc
       integer offset_i, offset_j

       real, allocatable, dimension(:,:,:,:)  :: blk_Evec
       integer, allocatable, dimension(:,:,:) :: blk_boti

       ! Declare work variables used in SGELS.
       real eBottom
       real dz
       real sumAu, tempA, signAu1
       real, allocatable, dimension(:,:) :: Au

       ! Declare output variables.
       real, allocatable, dimension(:,:,:) :: glb_ebot
       
       !Pass in parameters.
       parameter(OVERWRITE=.TRUE.)
       parameter(idm=9000,jdm=7055,kdm=240)
       parameter(blki=150,blkj=200,mbdy=3)
       parameter(IL=1,IU=5)
       parameter(MEIG=5)
       parameter(dz=25.0)

       ! TODO: Temporary parameters jblks, jblke, iblks, iblke
       parameter(runid='221')
       parameter(jblks=1,jblke=35,iblks=1,iblke=60)

       ! Define some variables.
       data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/

       ! Allocate space for 2D, 3D, 4D variables.
       allocate( blk_Evec(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:kdm,
     & 1:MEIG))
       allocate( glb_ebot(1:idm, 1:jdm, 1:MEIG) )
       allocate( blk_boti(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:2) )

       ! Allocate space for subroutine variables.
       allocate( Au(1:kdm,1:MEIG) )

       write(6,'(a,a3)'),'runid: ',runid
       call flush(6)

       write(6,'(a,4i3)') 'jblks jblke iblks iblke: ',jblks, jblke,
     & iblks, iblke
       call flush(6)

       ! Loop over blocks from input file.
       do jblk=jblks,jblke
       do iblk=iblks,iblke
         offset_i = blki*(iblk-1)
         offset_j = blkj*(jblk-1)

         write(6,*) 'iblk, jblk = ',iblk,jblk
         call flush(6)

         ! Read in the eigenvectors.       
         infileE='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/modes/'
     &   //'EVEC_1z_'//runid//'_'//numstr(jblk)//'_'//numstr(iblk)//
     &   '.BinF'

         write(6,*) 'reading ',infileE
         call flush(6)

         open(unit=11,file=infileE,status='old',form='unformatted')

         !Initialize blk_Evec
         blk_Evec(:,:,:,:)=0.0
         !Read blk_Evec in same way it was written out.
         do l=1,MEIG
           read(11) blk_Evec(:,:,:,l)
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

         !Loop over i,j, compute dW/dz, modal decomposition.
         do j=1-mbdy,blkj+mbdy
           do i=1-mbdy,blki+mbdy

             ! Find first good cell above the bottom.
             cc=blk_boti(i,j,1)
!             write(6,*) 'i,j, cc = ',i,' ',j,' ',cc
!             call flush(6)
             ! Do a check on boti -- any good cells?
             if (cc.GT.1) then   
         
              ! Initialize the A,BU, BV matrices to 0.
              ! Split A into Au, Av because SGELS alters A input.
              Au(:,:)=0.0

              ! Make columns of A = dW/dz.
              do k=1,cc-1
                do l=1,MEIG
                  Au(k,l)=1.0/dz*(blk_Evec(i,j,k+1,l)-blk_Evec(i,j,k,l))
!                  write(6,*) 'A: ', Au(k,l)
!                  call flush(6)
                enddo !l
              enddo !k

              ! Normalize Au, Av.
              ! I should load Depth, and normalize using
              ! (sum(Ueigl.^2.*DZ,1)./H).^(1/2)
              ! But in this case, DZ=25 is constant, and 
              ! approximating H as sum(DZ where Ueigl!=0) feels ok.
              ! 8-23-17: Forgot something! H=DZ*cc, so I should divide
              ! by cc inside the square root -- otherwise my answers are 
              ! O(10^5) instead of O(10^4) J. GS.
              do l=1,MEIG
                sumAu = 0.0
!                write(6,*) 'real cc: ', real(cc-1)
!                call flush(6)
                sumAu = ( sum(Au(1:cc-1,l)**2)/real(cc-1))**0.5
!                write(6,*) 'sumA: ',sum(Au(1:cc-1,l),1)
!                call flush(6)
!                write(6,*) 'sumA2: ', sum( Au(1:cc-1,l)**2,1)
!                call flush(6)
!                write(6,*) 'sumA: ', sumAu
!                call flush(6)

                ! Flip around so signs match.
                signAu1=1
                if (Au(1,l).lt.0) then
                  signAu1=-1
                endif

                ! Normalize -- divide by sum of squares.
                eBottom = 0.0
                if (sumAu.eq.0.0) then
                  Au(k,l)=0
                else
                  do k=1,(cc-1)
                    tempA = Au(k,l)
                    Au(k,l) = signAu1*tempA / sumAu
                    
                    ! I want to get the last non-zero value of the
                    ! normalized eigenvector.
                    if ( Au(k,l) .NE. 0.0) then
                      eBottom = Au(k,l)
                    endif ! check if Au(k,l) is 0.
                  enddo !k
                endif ! Check sumAu==0.
!                write(6,*) 'Au: ',Au(cc-1:cc+1,l)
!               call flush(6)
                glb_ebot( offset_i+i,offset_j+j,l) = eBottom
       
               ! write(6,*)  Au(cc-1,l),' ', Au(cc,l),' ',eBottom
                !call flush(6)
              enddo !l
              ! Finished normalization.
             endif ! Test whether cc.GT.0.
           enddo !i
         enddo !j
       enddo !iblk
       enddo !jblk

         outfilE='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &   //'modes_global/'
     &   //'glb_ebot_'//runid//'.BinF'

         write(6,*) 'Writing ebot to ',outfilE
         call flush(6)
         open(unit=98,file=outfilE,status='replace',
     &       form='unformatted',access='sequential')

         write(98) glb_ebot(:,:,:)
          
         close(98)
       write(6,*) 'Finished ...'
       call flush(6)

       stop
       end program ebot

