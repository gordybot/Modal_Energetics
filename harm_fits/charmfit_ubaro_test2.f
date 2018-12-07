       program charmfit2

       IMPLICIT NONE

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
       integer idm, jdm
       integer maxobstot
       integer IL,IU,MEIG
       integer i,j,k,t,l,t1,t2
       integer iblk, iblks, iblke
       integer jblk, jblks, jblke
       real omm2    ! M2 frequency in per hours.
       real meanHF, meanB, costerm, sinterm
       integer maxTime

       ! Declare output variables.
       real, allocatable, dimension(:,:,:) :: HF_out
       real, allocatable, dimension(:,:,:) :: THING
       ! Declare work variables used in SGELS.
       real, allocatable, dimension(:,:) :: Ahf
       real, allocatable, dimension(:) :: Bhf
       real, allocatable, dimension(:) :: WORK
       integer LWORK, INFO

       ! Pass in parameters.
       parameter(idm=9000,jdm=7055)
       parameter(blki=150,blkj=200,mbdy=3)
       parameter(IL=1,IU=5)
       parameter(maxTime=679)
       parameter(maxobstot=383) ! 167 for expt_05.1
       parameter(THINGstr='ubaro')
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
       allocate( THING(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:maxTime))
       allocate( HF_out(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,
     &            1:3) )
       ! Allocate space for subroutine variables.
       LWORK=4*maxobstot*MEIG
       INFO=-10
       allocate( Ahf(1:maxobstot,3) )
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

         ! Load THING.
!         infile_THING='/p/work1/grstephe/hycom/GLBc0.04/expt_05.1/'
!     &   //THINGstr//'/'
!     &   //THINGstr//'_'//runid//'_blk_'//numstr(jblk)//'_'//
!     &   numstr(iblk)
!     &   //'.BinF'
         infile_THING = '/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &      //'ubaro/ubaro_221_blk_'//numstr(jblk)//'_'//
     &        numstr(iblk)//'.BinF'
         write(6,*) 'Reading THING from ',infile_THING
         call flush(6)
         open(unit=13,file=infile_THING,status='old',
     &     form='unformatted',action='read')

         ! Initialize values.
           THING(:,:,:)=0.0

         ! Read values from file.
           do t=1,maxobstot
             read(13) THING(:,:,t)
           enddo !t
         close(13)

       ! Initialize storage variables for mean HFit properties.
       ! mHFu(:,:,:,:)=0.0
       ! mHFv(:,:,:,:)=0.0

       !Loop over i,j, compute dW/dz, modal decomposition.
       do j=1-mbdy,blkj+mbdy
         do i=1-mbdy,blki+mbdy

         Ahf(:,:)=0.0
         Bhf(:)=0.0
         do k=1,maxobstot
           t2=k
           ! I had the next line = t2 instead of 1.0
           Ahf(k,1)=1.0
           Ahf(k,2)=COS(omm2*t2)
           Ahf(k,3)=SIN(omm2*t2)
         enddo !k
         ! Fill in Bhf
         !Bhf = modeU(i,j,(6*t+1-6):(6*t+1+6),:) 
         !Bhf = modeU(i,j,1:maxobstot,:)
         Bhf = THING(i,j,1:maxobstot)
!        call SGELS('N',kdm,MEIG+1,maxobstot,A,kdm,BV,kdm,
         !call SGELS('N',13,3,MEIG,Ahf,13,Bhf,13,
!         call SGELS('N',maxobstot,3,MEIG,Ahf,maxobstot,Bhf,
!     &            maxobstot,WORK(1:LWORK), LWORK, INFO)
         meanB = 1.0/maxobstot * sum(Bhf)
         Bhf(:) = Bhf(:) - meanB
         costerm = 2.0/maxobstot * sum( Bhf(:)*Ahf(:,2) )
         sinterm = 2.0/maxobstot * sum( Bhf(:)*Ahf(:,3) )
         
         INFO=0
         ! Store mean, amp, phase of Bhf to HFu_out.
         if(INFO.EQ.0) then
!           do l=1,MEIG
             meanHF = Bhf(1)
             !ampHF =  (Bhf(2,l)**2+Bhf(3,l)**2)**0.5
             !phaseHF = ATAN2(Bhf(2,l),Bhf(3,l))
             HF_out(i,j,1)= meanB !meanHF
             HF_out(i,j,2)= costerm ! Bhf(2) ! COSINE term.
             HF_out(i,j,3)= sinterm !Bhf(3) ! SINE term.
            ! write(6,*) 'HF MAP (line 317): ',meanHF, ampHF, phaseHF
            ! call flush(6)
!           enddo !l
         ! Take average
         end if !INFO
       enddo !i
       enddo !j

       ! Write harmonic fit output...
       outfileHF='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     & //'HarmonicFit/'
     & //'chf_'//THINGstr//'_'//runid//'_blk_'//numstr(jblk)//
     &   '_'//numstr(iblk)
     & //'_test2.BinF'
       write(6,*) 'Writing HF_out to', outfileHF
       call flush(6)
       open(unit=88,file=outfileHF,status='replace',
     &       form='unformatted',access='sequential')
       write(88) HF_out(:,:,:)
       close(88)

       enddo ! iblk
       enddo ! jblk

       write(6,*) 'Finished ...'
       call flush(6)

       stop
       end program charmfit2
