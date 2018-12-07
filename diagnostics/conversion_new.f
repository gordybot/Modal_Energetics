      program conversion_modes
      !!! Compute conversion using harmonic fits to ubaro, vbaro
      !!! and pbot.
      !!! Harmonic fit output assumed to be [mean, cos, sin].
      !!! C_0n = -1/2*del(H) dot u_0* p_n(at z=-H)
      !!! conversion( ) = 0.5* (udhdx + vdhdy) * pbot
      !!! 11/17/17 - GS.
      
      implicit none

      ! Declare variables.
      integer i,j,l,iblk,jblk,blki,blkj, idm, jdm, kdm
      integer mbdy,MEIG,offset_i,offset_j
      character infile_ubaro*240, infile_vbaro*240, infile_depth*240
      character infilePSCX*240, infilePSCY*240, infile_pbot*240
      character outfile_conv*240
      character*3 runid 
      character*2 numstr(90)
      logical isThere
 
      real, allocatable, dimension(:,:) :: depth
      real, allocatable, dimension(:,:) :: blk_pscx
      real, allocatable, dimension(:,:) :: blk_pscy

      real, allocatable, dimension(:,:,:) :: ubaro
      real, allocatable, dimension(:,:,:) :: vbaro

      real, allocatable, dimension(:,:,:,:) :: pbot
      real, allocatable, dimension(:,:,:) :: glb_conv

      complex, allocatable, dimension(:,:) :: udhdx
      complex, allocatable, dimension(:,:) :: vdhdy

      ! Assign size parameters and define some useful variables.
      parameter(idm=9000, jdm=7055, kdm=240, mbdy=3)
      parameter(blki=150 ,blkj= 200, MEIG = 5)
      parameter(runid='221') 

      ! Define numstr.
       data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/

      ! Allocate memory for variables.
      allocate( depth(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy)  )
      allocate( blk_pscx( 1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy) )
      allocate( blk_pscy( 1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy) )

      allocate( ubaro( 1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:3) )
      allocate( vbaro( 1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:3) )

      allocate( pbot(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:3, 1:MEIG) )
      allocate( glb_conv(1:idm, 1:jdm, 1:MEIG) )

      allocate( udhdx(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy) )
      allocate( vdhdy(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy) )

      ! Initialize global fields.
      glb_conv(:,:,:) = 0.0

      do jblk=1,35
      do iblk=1,60
        ! i,j offsets for entry to global fields.
        offset_i = blki*(iblk-1)
        offset_j = blkj*(jblk-1)

        ! Initialize computed block fields.
        udhdx(:,:) = 0.0
        vdhdy(:,:) = 0.0
 
        !!! Read the depth.
        ! Load depth.
        infile_depth='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &       //'griddata/depth_221_blk_'//numstr(jblk)//'_'//
     &         numstr(iblk)//'.BinF'
        write(6,*) 'Reading ',infile_depth
        call flush(6)
        ! Initialize depth to 0.
        depth(:,:) = 0.0
        open(unit=14, file=infile_depth, form='unformatted',
     &            status='OLD', access='sequential',action='READ')
        read(14) depth(:,:)
        close(14)

        ! Read pscx, pscy from gridfiles.
        infilePSCX='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &      //'griddata/pscx_'//runid//'_blk_'//numstr(jblk)//'_'//
     &       numstr(iblk)//'.BinF'
        open(unit=15,file=infilePSCX,status='old',form='unformatted')
        blk_pscx(:,:) = 1.0
        read(15) blk_pscx(:,:)
        close(15)

        infilePSCY='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &       //'griddata/pscy_'//runid//'_blk_'//numstr(jblk)//'_'//
     &       numstr(iblk)//'.BinF'
        open(unit=16,file=infilePSCY,status='old',form='unformatted')
        blk_pscy(:,:) = 1.0
        read(16) blk_pscy(:,:)
        close(16)

        !load complex harmonic fits to ubaro, vbaro, pbot
        infile_ubaro='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &       'HarmonicFit/'//'chf_ubaro_'//runid//
     &       '_blk_'//numstr(jblk)//'_'//
     &          numstr(iblk)//'.BinF'
        write(6,*) 'Reading ',infile_ubaro
        call flush(6)
        ! Initialize ubaro to 0.
        ubaro(:,:,:) = 0.0
        open(unit=17, file=infile_ubaro, form='unformatted',
     &            status='OLD', access='sequential',action='READ')
        read(17) ubaro(:,:,:)
        close(17)

        !load complex harmonic fits to vbaro
        infile_vbaro='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &       'HarmonicFit/'//'chf_vbaro_'//runid//
     &       '_blk_'//numstr(jblk)//'_'//
     &          numstr(iblk)//'.BinF'
        write(6,*) 'Reading ',infile_vbaro
        call flush(6)
        ! Initialize vbaro to 0.
        vbaro(:,:,:) = 0.0
        open(unit=18, file=infile_vbaro, form='unformatted',
     &            status='OLD', access='sequential',action='READ')
        read(18) vbaro(:,:,:)
        close(18)

        !load complex harmonic fits to pbot
        infile_pbot='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &       'HarmonicFit/'//'chf_pbot_'//runid//
     &       '_blk_'//numstr(jblk)//'_'//
     &          numstr(iblk)//'_v2.BinF'
        write(6,*) 'Reading ',infile_pbot
        call flush(6)
        ! Initialize pbot to 0.
        pbot(:,:,:,:) = 0.0

      inquire( file=infile_pbot, EXIST = isThere)
      if (isThere) then
        open(unit=19, file=infile_pbot, form='unformatted',
     &            status='OLD', access='sequential',action='READ')
        read(19) pbot(:,:,:,:)
        close(19)

        !conversion( ) = 0.5* (udhdx + vdhdy) * pbot

        do i=1,blki+1
        do j=1,blkj+1
          ! Compute u_conjugate * dh/dx and v_conjugate * dh/dy
          udhdx(i,j) =  cmplx( ubaro(i,j,2), -ubaro(i,j,3)) * 
     &                   (depth(i+1,j)-depth(i,j))/blk_pscx(i,j)
          vdhdy(i,j) =  cmplx( vbaro(i,j,2), -vbaro(i,j,3)) * 
     &                   (depth(i,j+1)-depth(i,j))/blk_pscy(i,j)
        enddo !j
        enddo !i

        do i=1,blki
        do j=1,blkj
        do l=1,MEIG
          glb_conv( offset_i + i, offset_j + j, l) = real( 0.5*
     &               ( 0.5* (udhdx(i,j) + udhdx(i+1,j))
     &               + 0.5* (vdhdy(i,j) + vdhdy(i,j+1)) )
     &              *  cmplx( pbot(i,j,l,2), pbot(i,j,l,3)) )
        enddo !l 
        enddo !j
        enddo !i

       else
           write(6,*) 'File does not exist (yet).'
       endif ! infile_pbot exists.
      enddo !iblk
      enddo !jblk

      !!!!! Write the conversion to output file.
      ! Write the glb_conv variable by mode.
       do l=1,MEIG
         outfile_conv = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &        'modes_global/'//'glb_conv_1z_'//runid//
     &        '_mode'//numstr(l)//'.BinF'
       !
         write(6,*) 'Writing global bt conversion to ',outfile_conv
         call flush(6)
         open(unit=99,file=outfile_conv,status='replace',
     &     form='unformatted',access='sequential')
           write(99) glb_conv(:,:,l)
         close(99)
       enddo !l

       write(6,*) 'Finished ...'
       call flush(6)

       stop
       end program conversion_modes
