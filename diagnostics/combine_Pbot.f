       program combine_Ptiles
       !!! Program to combine the blocks I just made into
       !!! global fields. 
       !!! Input files are the harmonic fits
       !!! to the vertical mode decompositions 
       !!! from evec_hfit.f
       !!! Written in Cedar Falls, Iowa. 6/30/2017 GRS
       !!! Modified to combine pressure modes only, 
       !!! check if file exists before trying to open.
       !!! Written in Oakland, 09/01/2017 GRS.
       !!!  - Compute APE using DEPTH, EVAL, and P-mode amplitude.
       !!! APE = rho_0^-1 * p_amp^2 * H/4 * EVAL
       !!! where EVAL = k_n^2/(omm2^2-f^2), H is depth,
       !!! p_amp is harmonic-fit amplitude of nth pressure-mode
       !!! Updated, 09/26/17 GRS.

       implicit none
       ! Declare variables and types.
       character runid*3
       character infilePbot*240, infilePsurf*240
       character outfilePbot*240, outfilePsurf*240
       character*2 numstr(90)
       integer MEIG, IU, IL
       integer idm, jdm, blki, blkj, mbdy, iblk, jblk, l
       integer i,j
       logical isThere
       real PI, omm2, rho0, f0, f, p_n, eigval, zbot
       real, allocatable, dimension(:,:,:,:) :: glb_pbot
       real, allocatable, dimension(:,:,:,:) :: glb_psurf
       real, allocatable, dimension(:,:,:,:) :: blk_pbot
       real, allocatable, dimension(:,:,:,:) :: blk_psurf

       ! Define some useful numbers.
       !Pass in parameters.
       parameter(runid='221')
       parameter(idm=9000,jdm=7055)
       parameter(blki=150,blkj=200,mbdy=3)
       parameter(IL=1,IU=5)
       MEIG = IU-IL+1
       data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/

       ! Allocate space for blk fields.
       allocate( blk_pbot(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,
     & 1:2,1:MEIG) )
       allocate( blk_psurf(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,
     & 1:2,1:MEIG) )

       ! Allocate space for global fields.
       allocate( glb_pbot(1:idm,1:jdm,1:2,1:MEIG) )
       allocate( glb_psurf(1:idm,1:jdm,1:2,1:MEIG) )
 
       ! Initialize global fields to 0.0
       glb_psurf(:,:,:,:) = 0.0
       glb_pbot(:,:,:,:) = 0.0
  
       do jblk=1,35
         do iblk=1,60
           ! Initialize block variables.
           blk_pbot(:,:,:,:)=0.0
           blk_psurf(:,:,:,:) = 0.0

           infilePbot='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &     'modes/pbot_'//runid//'_blk_'//numstr(jblk)//
     &     '_'//numstr(iblk)//'_apr.BinF'
           
           write(6,*) 'i,j = ',iblk,jblk, infilePbot
           call flush(6)

           ! Check if file exists.
           inquire(file=infilePbot, EXIST=isThere )
           if (isThere) then        
             ! Read one block.
             write(6,*) 'reading ',infilePbot
             call flush(6)
          open(unit=11, file=infilePbot,status='old',form='unformatted')

             blk_pbot(:,:,:,:) = 0.0
             read(11) blk_pbot(:,:,:,:)
             close(11)

             ! Move blk data into global fields.
             glb_pbot( (blki*(iblk-1)+1):(blki*iblk),
     &         (blkj*(jblk-1)+1):(blkj*jblk),1,1:MEIG ) =
     &         blk_pbot(1:iblk,1:jblk,1,1:MEIG)

             glb_pbot( (blki*(iblk-1)+1):(blki*iblk),
     &         (blkj*(jblk-1)+1):(blkj*jblk),2,1:MEIG ) =
     &         blk_pbot(1:iblk,1:jblk,2,1:MEIG)

           else ! File does not exist (yet).
             write(6,*) infilePbot,' does not exist (yet).'
             call flush(6)
           end if

           ! Initialize block variables.
           blk_psurf(:,:,:,:)=0.0

           infilePsurf='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &     'modes/psurf_hfvm_'//runid//'_blk_'//numstr(jblk)//
     &     '_'//numstr(iblk)//'_may.BinF'

           write(6,*) 'i,j = ',iblk,jblk, infilePsurf
           call flush(6)

           ! Check if file exists.
           inquire(file=infilePsurf, EXIST=isThere )
           if (isThere) then
             ! Read one block.
             write(6,*) 'reading ',infilePsurf
             call flush(6)
         open(unit=11, file=infilePsurf,status='old',form='unformatted')

             blk_psurf(:,:,:,:) = 0.0
             read(11) blk_psurf(:,:,:,:)
             close(11)

             ! Move blk data into global fields.
             glb_psurf( (blki*(iblk-1)+1):(blki*iblk),
     &         (blkj*(jblk-1)+1):(blkj*jblk),1,1:MEIG ) =
     &         blk_psurf(1:iblk,1:jblk,1,1:MEIG)

             glb_psurf( (blki*(iblk-1)+1):(blki*iblk),
     &         (blkj*(jblk-1)+1):(blkj*jblk),2,1:MEIG ) =
     &         blk_psurf(1:iblk,1:jblk,2,1:MEIG)

           else ! File does not exist (yet).
             write(6,*) infilePsurf,' does not exist (yet).'
             call flush(6)
           end if


         enddo !iblk
       enddo ! jblk

       ! Write out bottom pressure fields,  one mode at a time.
       outfilePbot='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     & 'modes_global/glb_pbot_may.BinF'

       write(6,*) 'Writing global fields to ',outfilePbot
       call flush(6)
       open(unit=98,file=outfilePbot,status='replace',
     &   form='unformatted',access='sequential')

       do l=1,MEIG
         write(98) glb_pbot(:,:,:,l)
       enddo !l 
       close(98)

       ! Write out surface pressure fields,  one mode at a time.
       outfilePsurf='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     & 'modes_global/glb_psurf_may.BinF'

       write(6,*) 'Writing global fields to ',outfilePsurf
       call flush(6)
       open(unit=97,file=outfilePsurf,status='replace',
     &   form='unformatted',access='sequential')

       do l=1,MEIG
         write(97) glb_psurf(:,:,:,l)
       enddo !l 
       close(97)

       write(6,*) 'Finished ...'
       call flush(6)

       stop
       end program combine_Ptiles
