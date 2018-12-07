       program charmfit2

       IMPLICIT NONE
       ! TODO: Found a floating MEIG bug in charmfit_ubaro61.f?
       ! What does that impact?
       ! Nevermind, I specified it to be 1.

       ! Do a harmonic fit with complex output on THING.
       ! For now, assume THING has size blki, blkj, maxobstot

       ! Modify infile_THING at line 90.
       ! CHANGE THINGstr --- or risk overwriting another variable.
       ! CHANGE infile_THING -- to read correct variable.
       ! DON'T mix up 'ubaro' and 'vbaro'!

       ! Declare variables' types.
       character*3 runid
       character infileE*240, infileB*240
       character infile_THING*240, outfileHF*240
       character*2 numstr(90)
       character*5 THINGstr
       integer blki, blkj, mbdy
       integer idm, jdm, kdm
       integer maxobstot
       integer IL,IU,MEIG
       integer i,j,k,t,l,t1,t2, cc
       integer iblk, iblks, iblke
       integer jblk, jblks, jblke
       real omm2    ! M2 frequency in per hours.
       real meanHF
       integer maxTime

       ! Declare output variables.
       real, allocatable, dimension(:,:,:,:) :: HF_out
       real, allocatable, dimension(:,:,:,:) :: THING
       ! Declare work variables used in SGELS.
       real, allocatable, dimension(:,:) :: Ahf
       ! SGELS modifies Ahf, so I introduce a variable to be 
       ! reinitialized at every call of SGELS.
       real, allocatable, dimension(:,:) :: Awork
       real, allocatable, dimension(:) :: Bhf
       real, allocatable, dimension(:) :: WORK
       integer LWORK, INFO
       integer, allocatable, dimension(:,:,:) :: blk_boti

       ! Pass in parameters.
       parameter(idm=9000,jdm=7055, kdm=240)
       parameter(blki=150,blkj=200,mbdy=3)
       parameter(IL=1,IU=5)
       parameter(maxTime=679)
       parameter(maxobstot=383) ! 167 for expt_05.1
       parameter(THINGstr='v_iso')
       ! This number is based off files in input3z directory.

       ! Read in some parameters.
       read(*,*) jblks
       read(*,*) jblke
       read(*,*) iblks
       read(*,*) iblke
       read(*,*) runid
       !read(*,*) THINGstr

       MEIG = 1

       ! Allocate space for variables.
       allocate( THING(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:kdm,
     &                    1:maxobstot))
       allocate( HF_out(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,1:kdm,
     &           1:3) )
       allocate( blk_boti(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,1:2) )

       ! Allocate space for subroutine variables.
       LWORK=4*maxobstot*MEIG
       INFO=-10
       allocate( Ahf(1:maxobstot,3) )
       allocate( Awork(1:maxobstot,3) )
       allocate( Bhf( 1:maxobstot) )
       allocate( WORK(1:LWORK) )

       write(6,'(a,a3)'),'runid: ',runid
       call flush(6)

       write(6,'(a,4i3)') 'jblks jblke iblks iblke: ',jblks, jblke,
     & iblks, iblke
       call flush(6)

       ! Define some variables.
       omm2= 2*3.14159265/(12.42)
       data numstr/'01','02','03','04','05','06','07','08','09','10',
     &             '11','12','13','14','15','16','17','18','19','20',
     &             '21','22','23','24','25','26','27','28','29','30',
     &             '31','32','33','34','35','36','37','38','39','40',
     &             '41','42','43','44','45','46','47','48','49','50',
     &             '51','52','53','54','55','56','57','58','59','60',
     &             '61','62','63','64','65','66','67','68','69','70',
     &             '71','72','73','74','75','76','77','78','79','80',
     &             '81','82','83','84','85','86','87','88','89','90'/

       ! Loop over blocks from input file.
       do jblk=jblks,jblke
       do iblk=iblks,iblke
         write(6,*) 'iblk, jblk = ',iblk,jblk
         call flush(6)

         ! Load BOTI
         ! Read in the bottom bin locations.
         infileB='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/modes/'
     &      //'BOTI_1z_'//runid//'_'//numstr(jblk)//'_'//numstr(iblk)//
     &      '.BinF'
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

         ! Load THING.
!         infile_THING='/p/work1/grstephe/hycom/GLBc0.04/expt_05.1/'
!     &   //THINGstr//'/'
!     &   //THINGstr//'_'//runid//'_blk_'//numstr(jblk)//'_'//
!     &   numstr(iblk)
!     &   //'.BinF'
         infile_THING = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &      //THINGstr//'/v_iso_221_blk_'//numstr(jblk)//'_'//
     &        numstr(iblk)//'.BinF'
         write(6,*) 'Reading THING from ',infile_THING
         call flush(6)
         open(unit=13,file=infile_THING,status='old',
     &     form='unformatted',action='read')

         ! Initialize values.
           THING(:,:,:,:)=0.0

         ! Read values from file.
           do t=1,maxobstot
             do k=1,kdm
               read(13) THING(:,:,k,t)
             enddo !k
           enddo !t
         close(13)

       ! Initialize storage variables for mean HFit properties.
       ! mHFu(:,:,:,:)=0.0
       ! mHFv(:,:,:,:)=0.0

       !Loop over i,j, compute dW/dz, modal decomposition.
       do j=1-mbdy,blkj+mbdy
         do i=1-mbdy,blki+mbdy
           cc = blk_boti(i,j,1)

           Ahf(:,:)=0.0
           Bhf(:)=0.0
         do t2=1,maxobstot
           ! I had the next line = t2 instead of 1.0
           Ahf(t2,1)=1.0
           Ahf(t2,2)=COS(omm2*t2)
           Ahf(t2,3)=SIN(omm2*t2)
         enddo !t2
         
         if (cc.GT.1) then
         ! Loop over depths.
         do k=1,cc
           ! Re-initialize Awork to Ahf -- 
           ! This matrix changes everytime SGELS runs.
           Awork(:,:) = Ahf(:,:)

           ! Fill in Bhf
           !Bhf = modeU(i,j,(6*t+1-6):(6*t+1+6),:) 
           !Bhf = modeU(i,j,1:maxobstot,:)
           Bhf = THING(i,j,k,1:maxobstot)
!          call SGELS('N',kdm,MEIG+1,maxobstot,A,kdm,BV,kdm,
           !call SGELS('N',13,3,MEIG,Ahf,13,Bhf,13,
           call SGELS('N',maxobstot,3,MEIG,Awork,maxobstot,Bhf,
     &            maxobstot,WORK(1:LWORK), LWORK, INFO)

            ! Store mean, amp, phase of Bhf to HFu_out.
           if(INFO.EQ.0) then
             meanHF = Bhf(1)
             !ampHF =  (Bhf(2,l)**2+Bhf(3,l)**2)**0.5
             !phaseHF = ATAN2(Bhf(2,l),Bhf(3,l))
             HF_out(i,j,k,1)=meanHF
             HF_out(i,j,k,2)=Bhf(2) !cosine term
             HF_out(i,j,k,3)=Bhf(3) !sine term
            ! write(6,*) 'HF MAP (line 317): ',meanHF, ampHF, phaseHF
            ! call flush(6)
           end if !INFO
         enddo !k
         endif ! cc<1.
       enddo !i
       enddo !j

       ! Write harmonic fit output...
       outfileHF='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     & //'HarmonicFit/'
     & //'chf_'//THINGstr//'_'//runid//'_blk_'//numstr(jblk)//
     &   '_'//numstr(iblk)
     & //'_apr.BinF'
       write(6,*) 'Writing HF_out to', outfileHF
       call flush(6)
       open(unit=88,file=outfileHF,status='replace',
     &       form='unformatted',access='sequential')
       write(88) HF_out(:,:,:,:)
       close(88)

       enddo ! iblk
       enddo ! jblk

       write(6,*) 'Finished ...'
       call flush(6)

       stop
       end program charmfit2
