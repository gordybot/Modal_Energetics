       program vert_mode_UV

       implicit none

       !Declare variables' types here:
       character*3 runid
       character infileE*240, infileB*240
       character infileU*240, infileV*240
       character outfilU*240, outfilV*240
       character*2 numstr(90)
       integer blki, blkj, mbdy
       integer idm, jdm, kdm
       integer maxobstot
       integer IL,IU,MEIG
       integer i,j,k,t,l,t1,t2, q
       integer iblk, iblks, iblke
       integer jblk, jblks, jblke
       integer cc
       real omm2    ! M2 frequency in per hours.
       real grav, H, meanA    ! accel. due to gravity.
       real tempSumU1, tempSumU2, tempSumV1, tempSumV2

       integer blk_boti !location of the bottom
       real, allocatable, dimension(:,:,:,:) :: u_iso !i,j,k,"t=1:3"
       real, allocatable, dimension(:,:,:,:) :: v_iso !i,j,k
       real, allocatable, dimension(:,:,:,:) :: blk_Evec   
       real, allocatable, dimension(:,:) :: x
       real, allocatable, dimension(:,:) :: y

       ! Declare work variables used in SGELS.
       real dz
       real sumA, sumAp, tempA
       integer signA1, signAv1, signAp1
       real, allocatable, dimension(:,:) :: A
       real, allocatable, dimension(:,:) :: BU
       real, allocatable, dimension(:,:) :: BV
       real, allocatable, dimension(:) :: WORK
       integer LWORK, INFO1, INFO2, INFO3

       ! Declare output variables.
       ! modeU, modeV [i,j,t,l=mode number] are the
       ! amplitudes of the modal fits.
       real, allocatable, dimension(:,:,:,:) :: modeU
       real, allocatable, dimension(:,:,:,:) :: modeV

       !Pass in parameters.
       parameter(idm=9000,jdm=7055,kdm=240)
       parameter(blki=150,blkj=200,mbdy=3)
       parameter(IL=1,IU=5)
       parameter(dz=25.0)
       parameter(grav=9.81)
       parameter(maxobstot=3)! corresponding to 1, cos, sin terms.
       !parameter(maxobstot=679) !TODO: Should be 744 eventually?
       ! This number is based off files in input3z directory.

       !Specify dimensions.
       dimension blk_boti(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:2)

       ! Define some variables.
       omm2= 2*3.14159265/(12.42)
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

       ! Allocate space for 2D, 3D, 4D variables.
       allocate( u_iso(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:kdm,
     & 1:3) )
       allocate( v_iso(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:kdm,
     & 1:3) )
       allocate( blk_Evec(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:kdm,
     & 1:MEIG))

       ! Allocate space for subroutine variables.
       LWORK=4*kdm*MEIG
       INFO1=-10
       INFO2=-10
       allocate( A(1:kdm,1:(MEIG)) )
       allocate( BU(1:kdm, 1:maxobstot) )
       allocate( BV(1:kdm, 1:maxobstot) )
       allocate( WORK(1:LWORK) )
 
       ! Allocate. 
       allocate( x(1:MEIG, 1:2) )
       allocate( y(1:MEIG, 1:2) )

       ! Allocate space for output variables.
       allocate( modeU(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:2,
     & 1:MEIG) )
       allocate( modeV(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:2,
     & 1:MEIG) )

       ! TODO: Temporary parameters jblks, jblke, iblks, iblke
       !parameter(runid='051')
       !parameter(jblks=16,jblke=18,iblks=1,iblke=30)
       ! TODO: Read begin/end of blocks from com-generated file.
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

       ! Loop over blocks from input file.
       do jblk=jblks,jblke
       do iblk=iblks,iblke
         ! Initialize modeU, modeV
         modeU(:,:,:,:) = 0.0
         modeV(:,:,:,:) = 0.0

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

         ! Open the file for  u-velocity.
!         infileU='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/u_iso/'
         infileU='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &   'HarmonicFit/chf_u_iso_'//runid//'_blk_'//
     &   numstr(jblk)//'_'//numstr(iblk)//'_apr.BinF'
!   WAS  numstr(jblk)//'_'//nustr(iblk)//'_v18.BinF' before May 7, 2018.
         write(6,*) 'reading ',infileU
         call flush(6)
         open(unit=13,file=infileU,status='old',form='unformatted')

          ! e.g., chf_v_iso_061_blk_01_04.BinF
         ! Open the file for v-velocity.
!         infileV='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/v_iso/'
         infileV='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &   'HarmonicFit/chf_v_iso_'//runid//'_blk_'//
     &   numstr(jblk)//'_'//numstr(iblk)//'_apr.BinF'
!       CHANGED infileV from *_v18.BinF   on May 7 2018.
         write(6,*) 'reading ',infileV
         call flush(6)
         open(unit=14,file=infileV,status='old',form='unformatted')

         ! Initialize values.
         u_iso(:,:,:,:)=0.0
         v_iso(:,:,:,:)=0.0

         ! Read file all at once. (HF_out(:,:,:,:) )
         read(13) u_iso(:,:,:,:)
         read(14) v_iso(:,:,:,:)
         ! Close u, v, rho files.
         close(13)
         close(14)

         !Loop over i,j, compute dW/dz, modal decomposition.
         do j=1-mbdy,blkj+mbdy
           do i=1-mbdy,blki+mbdy
             ! Find first good cell above the bottom.
             ! NEIG says use the first, not the second.
             cc=blk_boti(i,j,1)
             ! Do a check on boti -- any good cells?
             if(cc.GT.1) then   
         
              ! Initialize the A,BU, BV matrices to 0.
              ! Split A into A, Av because SGELS alters A input.
              A(:,:)=0.0
              BU(:,:)=0.0
              BV(:,:)=0.0

              ! Make columns of A = dW/dz.
              ! Use 0 as first value of blk_Evec
              do l=1,MEIG
                A(1,l)=1.0/dz*(blk_Evec(i,j,1,l))
              enddo !l
              BU(1,1:2)=u_iso(i,j,1,2:3)
              BV(1,1:2)=v_iso(i,j,1,2:3)
              do k=2,cc
                do l=1,MEIG
                  A(k,l)=1.0/dz*(blk_Evec(i,j,k,l)-blk_Evec(i,j,k-1,l))
                enddo !l
                ! Interpolate u,v to cell edge to match dW/dz depth.
                do q=2,3
                  BU(k,q-1)=u_iso(i,j,k,q)
                  BV(k,q-1)=v_iso(i,j,k,q)
                enddo !t
              enddo !k

              ! Normalize A, Av.
              ! I should load Depth, and normalize using
              ! (sum(Ueigl.^2.*DZ,1)./H).^(1/2)
              ! But in this case, DZ=25 is constant, and 
              ! approximating H as sum(DZ where Ueigl!=0) feels ok.
              ! 8-23-17: Forgot something! H=DZ*cc, so I should divide
              ! by cc inside the square root -- otherwise my answers  
              ! O(10^5) instead of O(10^4) J. GS.
              do l=1,MEIG
                sumA = sum((A(:,l)**2)/cc)**0.5

                ! Flip around so signs match at seabed.
                signA1=1
                if (A(cc,l).lt.0) then
                  signA1=-1
                endif

                ! Normalize -- divide by sum of squares.
                do k=1,cc
                  tempA = A(k,l)
                  if (sumA.eq.0) then
                    A(k,l)=0
                  else
                    A(k,l) = signA1*tempA / sumA
                  endif
                
                !  if (i.EQ.50 .AND. j.EQ.50) then
                !    write(6,*) 'sum(A^2): ', sum(A(:,1)**2)
                !  endif
                enddo !k

                ! Make sure mean(A(:,l))==0.0
                meanA = sum(A(:,l))*1.0/cc
                do k=1,cc
                  A(k,l) = A(k,l) - meanA
                enddo !k
              enddo !l
              ! Finished normalization.

              ! Thoughts -- I should probably divide by cc rather than H
              ! and I should remove a vertical mean of BU(:,2) and BU(:,3)
              ! before fitting the vertical modes.
              ! Project normal modes onto B.
              x(:,:) = 0.0 
              y(:,:) = 0.0
              do l=1,MEIG
                ! Remove vertical mean of BU, BV
                !tempSumU1 = 1.0/cc*sum(BU(:,2))
                !tempSumU2 = 1.0/cc*sum(BU(:,3))
                !tempSumV1 = 1.0/cc*sum(BV(:,2))
                !tempSumV2 = 1.0/cc*sum(BV(:,3))
                !do k=1,cc
                !  BU(k,2) = BU(k,2) - tempSumU1
                !  BU(k,3) = BU(k,3) - tempSumU2
                !  BV(k,2) = BV(k,2) - tempSumV1
                !  BV(k,3) = BV(k,3) - tempSumV2
                !enddo !k
                tempSumU1 = 0.0
                tempSumU2 = 0.0
                tempSumV1 = 0.0
                tempSumV2 = 0.0
                do k=1,cc
                  tempSumU1 = tempSumU1 + 1.0/cc * A(k,l) * BU(k,1)
                  tempSumU2 = tempSumU2 + 1.0/cc * A(k,l) * BU(k,2)
                  ! 
                  tempSumV1 = tempSumV1 + 1.0/cc * A(k,l) * BV(k,1)
                  tempSumV2 = tempSumV2 + 1.0/cc * A(k,l) * BV(k,2)
                enddo !k
                x(l,1) = tempSumU1
                x(l,2) = tempSumU2
                y(l,1) = tempSumV1
                y(l,2) = tempSumV2
                ! Remove last mode from the data.
                do k=1,kdm
                  BU(k,1) = BU(k,1) - tempSumU1*A(k,l)
                  BU(k,2) = BU(k,2) - tempSumU2*A(k,l)
                  ! 
                  BV(k,1) = BV(k,1) - tempSumV1*A(k,l)
                  BV(k,2) = BV(k,2) - tempSumV2*A(k,l)
                enddo !k
              enddo !l

              ! Store the first MEIG rows of BU, BV -- these are the
              ! solutions to A x = B.
              do l=1,MEIG
                do t=1,2
                  modeU(i,j,t,l)=x(l,t)
                  modeV(i,j,t,l)=y(l,t)
                enddo !t
              enddo !l
             endif ! Test whether cc.GT.0.
           enddo !i
         enddo !j

         outfilU='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/modes/'
     &   //'uhfvm_'//runid//'_blk_'//numstr(jblk)//'_'//numstr(iblk)
     &   //'_may.BinF'
         write(6,*) 'Writing modeU to ',outfilU
         call flush(6)
         open(unit=98,file=outfilU,status='replace',
     &       form='unformatted',access='sequential')

         outfilV='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/modes/'
     &   //'vhfvm_'//runid//'_blk_'//numstr(jblk)//'_'//numstr(iblk)
     &   //'_may.BinF'
         write(6,*) 'Writing modeV to ',outfilV
         call flush(6)
         open(unit=99,file=outfilV,status='replace',
     &       form='unformatted',access='sequential')

         write(98) modeU(:,:,:,:)
         write(99) modeV(:,:,:,:)
         close(98)
         close(99)

         write(6,*) 'iblk / jblk / MIN U / MAX U: '
         call flush(6)
         write(6,*) iblk,' ',jblk,' ',MINVAL(modeU),' ', MAXVAL(modeU)
         call flush(6)

         write(6,*) 'iblk / jblk / MIN V / MAX V: '
         call flush(6)
         write(6,*) iblk,' ',jblk,' ',MINVAL(modeV),' ', MAXVAL(modeV)
         call flush(6)

       enddo !iblk
       enddo !jblk

       write(6,*) 'Finished ...'
       call flush(6)

       stop
       end program vert_mode_UV
