       program combine_UVtiles
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
       character infileU*240, infileV*240, infileK*240
       character infile_plat*240, infile_depth*240
       character outfileP*240
       character outfileKE*240
       character*2 numstr(90)
       integer MEIG, IU, IL
       integer idm, jdm, blki, blkj, mbdy, iblk, jblk, l
       integer i,j
       logical isThere
       real PI, omm2, rho0, f0, f, p_n, eigval, zbot

       real, allocatable, dimension(:,:,:,:) :: blk_modeU
       real, allocatable, dimension(:,:,:,:) :: blk_modeV

       ! Variables used to compute APE.
       ! Need global plat, depth
       ! blk_ EVAL and blk_Pmode
       real, allocatable, dimension(:,:) :: plat
       real, allocatable, dimension(:,:) :: depth
       real, allocatable, dimension(:,:,:) :: glb_KE
       
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
       allocate( blk_modeU(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,
     & 1:MEIG, 1:3) )
       allocate( blk_modeV(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,
     & 1:MEIG, 1:3) )
       allocate( plat(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy) )
       allocate( depth(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy) )

       ! Allocate space for global fields.
       allocate( glb_KE(1:idm, 1:jdm, 1:MEIG) )
 
       ! Initialize global fields to 0.0
       glb_KE(:,:,:) = 0.0
  
       do jblk=1,35
         do iblk=1,60
           ! Initialize block variables.
           blk_modeU(:,:,:,:)=0.0
           blk_modeV(:,:,:,:)=0.0
           plat(:,:) = 0.0
           depth(:,:) =0.0
!           infileP='/p/work1/grstephe/hycom/GLBc0.08/expt_05.1/'//
!     & 'modes/mean_vmhfP_'//runid//'_blk_'//numstr(jblk)//
!     & '_'//numstr(iblk)//'_v918.BinF'

           infileU='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     & 'HarmonicFit/chf_modeU_'//runid//'_blk_'//numstr(jblk)//
     &  '_'//numstr(iblk)//'.BinF'
           
           infileV='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     & 'HarmonicFit/chf_modeV_'//runid//'_blk_'//numstr(jblk)//
     &  '_'//numstr(iblk)//'.BinF'

           write(6,*) 'i,j = ',iblk,jblk, infileU
           call flush(6)

           ! Check if file exists.
           inquire(file=infileU, EXIST=isThere )
           if (isThere) then        
             ! Read one block.
             write(6,*) 'reading ',infileU
             call flush(6)
             open(unit=11, file=infileU,status='old',form='unformatted')
             read(11) blk_modeU(:,:,:,:)
             close(11)

             write(6,*) 'reading ',infileV
             call flush(6)
             open(unit=12, file=infileV,status='old',form='unformatted')
             read(12) blk_modeV(:,:,:,:)
             close(12)

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

             ! Move blk data into global fields.
             do l=1,MEIG
               glb_KE( (blki*(iblk-1)+1):(blki*iblk),
     &         (blkj*(jblk-1)+1):(blkj*jblk), l ) =
     &           ( blk_modeU(1:blki,1:blkj, l ,2)**2 +
     &             blk_modeU(1:blki,1:blkj, l ,3)**2 +
     &             blk_modeV(1:blki,1:blkj, l ,2)**2 +
     &             blk_modeV(1:blki,1:blkj, l ,3)**2 )
     &            * rho0 * depth(1:blki, 1:blkj) / 4
             enddo !l

           else ! File does not exist (yet).
             write(6,*) infileU,' does not exist (yet).'
             call flush(6)
           end if
         enddo !iblk
       enddo ! jblk

       ! Write out global fields one mode at a time.

       ! Write the Kinetic Energy to its own file.
       ! One file for each mode.
       do l=1,MEIG 
         outfileKE='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     & 'modes_global/glb_KE_mode'//numstr(l)//'.BinF'

         write(6,*) 'Writing global KE to ',outfileKE
         call flush(6)
       
         open(unit=99,file=outfileKE,status='replace',
     &   form='unformatted',access='sequential')

         write(99) glb_KE(:,:,l)
         close(99)
       enddo !l
       
       write(6,*) 'Finished ...'
       call flush(6)

       stop
       end program combine_UVtiles
