       program save_pgrid
       !! Program to combine the plon, plat tiles into a global field
       !! of lat, lon, depth, that I can read in a binary format.
       implicit none

       ! Declare variables and types.
       character runid*3
       character infileU*240, infileV*240, infileK*240
       character infile_lat*240, infile_depth*240
       character infile_lon*240
       character outfile*240
       character*2 numstr(90)
       integer MEIG, IU, IL
       integer idm, jdm, blki, blkj, mbdy, iblk, jblk, l
       integer i,j, offset_i, offset_j
       logical isThere
       real PI, omm2, rho0, f0, f, p_n, eigval, zbot

       ! Variables used to compute APE.
       ! Need global plat, depth
       ! blk_ EVAL and blk_Pmode
       real, allocatable, dimension(:,:) :: plat
       real, allocatable, dimension(:,:) :: plon
       real, allocatable, dimension(:,:) :: depth
       real, allocatable, dimension(:,:) :: glb_lon
       real, allocatable, dimension(:,:) :: glb_lat
       real, allocatable, dimension(:,:) :: glb_depth
       
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
       ! Corrected allocation size for blk_modeU / blk_modeV.
       allocate( plat(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy) )
       allocate( plon(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy) )
       allocate( depth(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy) )

       ! Allocate space for global fields.
       allocate( glb_lat(1:idm, 1:jdm) )
       allocate( glb_lon(1:idm, 1:jdm) )
       allocate( glb_depth(1:idm, 1:jdm) )
 
       ! Initialize global fields to 0.0
       glb_lat(:,:) = 0.0
       glb_lon(:,:) = 0.0
       glb_depth(:,:) = 0.0
  
       do jblk=1,35
         do iblk=1,60
           ! Initialize block variables.
           plat(:,:) = 0.0
           plon(:,:) = 0.0
           depth(:,:) =0.0

           write(6,*) 'i,j = ',iblk,jblk, infileU
           call flush(6)
            
            infile_lon='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &       //'griddata/plon_221_blk_'//numstr(jblk)//'_'//
     &         numstr(iblk)//'.BinF'
            write(6,*) 'Reading ',infile_lon
            call flush(6)
            open(unit=11, file=infile_lon, form='unformatted',
     &           status='OLD',access='sequential',action='READ')
            plon(:,:) = 0.0
            read(11) plon(:,:)
            close(11)

            infile_lat='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &       //'griddata/plat_221_blk_'//numstr(jblk)//'_'//
     &         numstr(iblk)//'.BinF'
            write(6,*) 'Reading ',infile_lat
            call flush(6)
            open(unit=12, file=infile_lat, form='unformatted',
     %           status='OLD',access='sequential',action='READ')
            plon(:,:) = 0.0
            read(12) plon(:,:)
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
             ! Changed the third index from [2,3] to [1,2].
             !   EDIT / TODO check. 2.2.18.GS
             offset_i = blki*(iblk-1)
             offset_j = blkj*(jblk-1)
             glb_lon( offset_i + 1:offset_i + blki, 
     &                offset_j + 1:offset_j + blkj) = plon
             glb_lat( offset_i + 1:offset_i + blki,
     &                offset_j + 1:offset_j + blkj) = plat
             glb_depth( offset_i + 1:offset_i + blki,
     &                  offset_j + 1:offset_j + blkj) = depth

         enddo !iblk
       enddo ! jblk

       ! Write out global fields one at a time.
         outfile='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     & 'modes_global/glb_lon_221.BinF'
         write(6,*) 'Writing global longitude to ',outfile
         call flush(6)
         open(unit=99,file=outfile,status='replace',
     &        form='unformatted',access='sequential')
         write(99) glb_lon(:,:)
         close(99)
      
         outfile='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &  'modes_global/glb_lat_221.BinF'
         write(6,*) 'Writing global latitude to ',outfile
         call flush(6)
         open(unit=99,file=outfile,status='replace',
     &        form='unformatted',access='sequential')
         write(99) glb_lat(:,:)
         close(99)

         outfile='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &      'modes_global/glb_depth_221.BinF'
         write(6,*) 'Writing global depth to ',outfile
         call flush(6)
         open(unit=99,file=outfile,status='replace',
     &       form='unformatted',access='sequential')
         write(99) glb_depth(:,:)
         close(99)

       write(6,*) 'Finished ...'
       call flush(6)

       stop
       end program save_pgrid
