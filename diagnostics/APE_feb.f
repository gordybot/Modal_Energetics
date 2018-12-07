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
       character infileP*240, infileK*240
       character infile_plat*240, infile_depth*240
       character outfileP*240
       character outfileAPE*240
       character*2 numstr(90)
       integer MEIG, IU, IL
       integer idm, jdm, blki, blkj, mbdy, iblk, jblk, l
       integer i,j
       logical isThere
       real PI, omm2, rho0, f0, f, p_n, eigval, zbot
       real, allocatable, dimension(:,:,:,:) :: glb_pbot
       real, allocatable, dimension(:,:,:) :: glb_mHFp_cos
       real, allocatable, dimension(:,:,:) :: glb_mHFp_sin

       real, allocatable, dimension(:,:,:,:) :: blk_mHFp

       ! Variables used to compute APE.
       ! Need global plat, depth
       ! blk_ EVAL and blk_Pmode
       real, allocatable, dimension(:,:) :: plat
       real, allocatable, dimension(:,:) :: depth
       real, allocatable, dimension(:,:,:) :: eval
       real, allocatable, dimension(:,:,:) :: glb_APE
       
       ! Define some useful numbers.
       !Pass in parameters.
       parameter(runid='221')
       parameter(idm=9000,jdm=7055)
       parameter(blki=150,blkj=200,mbdy=3)
       parameter(IL=1,IU=5)
       MEIG = IU-IL+1
       PI = 3.1415926535897932384  ! Approximately.
       omm2 = 2*PI / (12.42*3600)  ! M2 frequency
       rho0 = 1034.0               ! Close enough.
       f0 = 2*(2*PI/(24*3600))     ! f = f0*sin(latitude)
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
       allocate( blk_mHFp(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,
     & 1:2,1:MEIG) )
       allocate( eval(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:MEIG) )
       allocate( plat(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy) )
       allocate( depth(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy) )

       ! Allocate space for global fields.
       allocate( glb_pbot(1:idm,1:jdm,1:2,1:MEIG) )
       allocate( glb_mHFp_cos(1:idm,1:jdm,1:MEIG) )
       allocate( glb_mHFp_sin(1:idm,1:jdm,1:MEIG) )
       allocate( glb_APE(1:idm, 1:jdm, 1:MEIG) )
 
       ! Initialize global fields to 0.0
!       glb_mHFp_mean(:,:,:) = 0.0
       glb_mHFp_cos(:,:,:) = 0.0
       glb_mHFp_sin(:,:,:) = 0.0
       glb_APE(:,:,:) = 0.0
  
       do jblk=1,35
         do iblk=1,60
           ! Initialize block variables.
           blk_mHFp(:,:,:,:)=0.0
           plat(:,:) = 0.0
           depth(:,:) =0.0
           eval(:,:,:) = 0.0
!           infileP='/p/work1/grstephe/hycom/GLBc0.04/expt_05.1/'//
!     & 'modes/mean_vmhfP_'//runid//'_blk_'//numstr(jblk)//
!     & '_'//numstr(iblk)//'_v918.BinF'

           infileP='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
!     &     'HarmonicFit/hfvmP_'//runid//'_blk_'//numstr(jblk)//
!     &     'modes/P_vmhf_'//runid//'_blk_'//numstr(jblk)//
     &     'modes/pbot_hfvm_'//runid//'_blk_'//numstr(jblk)//
     &     '_'//numstr(iblk)//'_apr.BinF'
           
           write(6,*) 'i,j = ',iblk,jblk, infileP
           call flush(6)

           ! Check if file exists.
           inquire(file=infileP, EXIST=isThere )
           if (isThere) then        
             ! Read one block.
             write(6,*) 'reading ',infileP
             call flush(6)
             open(unit=11, file=infileP,status='old',form='unformatted')

             blk_mHFp(:,:,:,:) = 0.0
             read(11) blk_mHFp(:,:,:,:)
             close(11)

             ! Move blk data into global fields.
             glb_pbot( (blki*(iblk-1)+1):(blki*iblk),
     &         (blkj*(jblk-1)+1):(blkj*jblk),1,1:MEIG ) =
     &         blk_mHFp(1:iblk,1:jblk,1,1:MEIG)

             glb_pbot( (blki*(iblk-1)+1):(blki*iblk),
     &         (blkj*(jblk-1)+1):(blkj*jblk),2,1:MEIG ) =
     &         blk_mHFp(1:iblk,1:jblk,2,1:MEIG)

!             glb_mHFp_cos( (blki*(iblk-1)+1):(blki*iblk),
!     &         (blkj*(jblk-1)+1):(blkj*jblk),1:MEIG ) =
!     &         blk_mHFp(1:iblk,1:jblk,1:MEIG,1)

!             glb_mHFp_sin( (blki*(iblk-1)+1):(blki*iblk),
!     &         (blkj*(jblk-1)+1):(blkj*jblk),1:MEIG ) =
!     &          blk_mHFp(1:iblk,1:jblk,1:MEIG,2)

            ! Load Maarten's eigenvalues 
            ! From mbui/hycom/GLBc0.04/eigenval/NEIG_c008_v3.f 
            ! EVAL = k_n^2/(omm2^2-f^2)
            infileK = '/p/work1/mbui/hycom/GLBc0.04/expt_22.1/modes/'
     &        //'EVAL_1z_'//runid//'_'//numstr(jblk)//'_'//
     &          numstr(iblk)//'.BinF'
             open(unit=12,file=infileK,status='old',form='unformatted')
             write(6,*) 'Reading ',infileK
             call flush(6)
             eval(:,:,:) = 0.0
             do l=1,MEIG
               read(12) eval(:,:,l)
             enddo !l. 
             close(12)

            ! Load plat to compute Coriolis frequency.
            ! f= 2*(2*PI/(24*3600))*SIN(lat_i*PI/180)
            ! use f0 = 2*(2*PI/(24*3600))
            infile_plat='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &       //'griddata/plat_221_blk_'//numstr(jblk)//'_'//
     &         numstr(iblk)//'.BinF'
             write(6,*) 'Reading ',infile_plat
             call flush(6)
             open(unit=13, file=infile_plat, form='unformatted',
     &            status='OLD', access='sequential',action='READ')
            plat(:,:) = 0.0
            read(13) plat(:,:)
            close(13)

            ! Load depth.
            infile_depth='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &       //'griddata/depth_221_blk_'//numstr(jblk)//'_'//
     &         numstr(iblk)//'.BinF'
             write(6,*) 'Reading ',infile_depth
             call flush(6)
             open(unit=14, file=infile_depth, form='unformatted',
     &            status='OLD', access='sequential',action='READ')
            depth(:,:) = 0.0
            read(14) depth(:,:)
            close(14)

            ! Calculate APE.
            do i = 1,blki
            do j = 1,blkj
              f = f0 * SIN( plat(i,j)*PI/180 )
              zbot = depth(i,j)

              
            do l = 1,MEIG
              ! Change third index from [1,2] to [2,3]. EDIT 2.2.18.GS
            p_n = ( blk_mHFp(i,j,1,l)**2 + blk_mHFp(i,j,2,l)**2 )**0.5
              eigval = eval(i,j,l)

              if (l.EQ.1 .AND. i.EQ.50 .AND. j.EQ.50) then
               write(6,*) 'jblk, iblk, eval: ',jblk,',',iblk,', ',eigval
               call flush(6)
              endif

              glb_APE( blki*(iblk-1)+i, blkj*(jblk-1)+j,l ) = 1.0/rho0 
     &           * p_n**2 * zbot/4 * eigval
            enddo !l
            enddo !j
            enddo !i

           else ! File does not exist (yet).
             write(6,*) infileP,' does not exist (yet).'
             call flush(6)
           end if
         enddo !iblk
       enddo ! jblk

       ! Write out global fields one mode at a time.
       outfileP='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     & 'modes_global/glb_pbot_hfvm_may.BinF'

       write(6,*) 'Writing global fields to ',outfileP
       call flush(6)
       
       open(unit=98,file=outfileP,status='replace',
     &   form='unformatted',access='sequential')

       do l=1,MEIG
         write(98) glb_pbot(:,:,:,l)
!         write(98) glb_mHFp_cos(:,:,l)
!         write(98) glb_mHFp_sin(:,:,l)
       enddo !l
       
       close(98)

       ! Write the Available Potential Energy to its own file.
       ! One file for each mode.
       do l=1,MEIG 
         outfileAPE='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     & 'modes_global/glb_APE_mode'//numstr(l)//'_may.BinF'

         write(6,*) 'Writing global APE to ',outfileAPE
         call flush(6)
       
         open(unit=99,file=outfileAPE,status='replace',
     &   form='unformatted',access='sequential')

         write(99) glb_APE(:,:,l)
       close(99)
       enddo !l
       
       write(6,*) 'Finished ...'
       call flush(6)

       stop
       end program combine_Ptiles
