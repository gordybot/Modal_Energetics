      PROGRAM CONVERSION

      IMPLICIT NONE
      ! Inputs:
      ! blk_ubaro, blk_vbaro, z, blk_pbot
      ! These are the results of complex_harmonic_fit
      ! aka charmfit_pbot.f

      ! Grid cell has       N
      !                  W  C  E
      !                     S
      ! U is defined at W, V at S, P and Z at C.
      ! E_i = W_i+1.

      ! Declare: z, pscx, pscy, dhdx_W, dhdhx_E, dhdy_N, dhdy_S
      ! Define z, pscx, pscy
      integer i,j,l,MEIG, offset_i, offset_j
      integer iblk,jblk,iblks,iblke,jblks,jblke
      integer idm,jdm, halo, blki, blkj
      character*2 numstr(90)
      character*3 runid
      character outfil*240, infilePbot*240,infileUbt*240,infileVbt*240
      character infileZ*240,infileDX*240,infileDY*240
      complex u_conj_C, v_conj_C, p_botn_C
      real dhdx_W, dhdx_E, dhdy_N, dhdy_S

      real, allocatable, dimension(:,:,:,:) :: blk_pbot
      real, allocatable, dimension(:,:,:) :: blk_ubaro
      real, allocatable, dimension(:,:,:) :: blk_vbaro
      real, allocatable, dimension(:,:,:) :: blk_conv
      real, allocatable, dimension(:,:,:) :: glb_conv
      real, allocatable, dimension(:,:) :: pscx
      real, allocatable, dimension(:,:) :: pscy
      real, allocatable, dimension(:,:) :: z

      ! Define some parameters and useful values.
      parameter(idm=9000,jdm=7055,MEIG=5)
      parameter(blki=150, blkj=200, halo=3)
      parameter(runid='221')
      parameter(jblks=1,jblke=35,iblks=1,iblke=60)

      ! Allocate 2D and higher variables.
      allocate( z( 1-halo:blki+halo, 1-halo:blkj+halo) )
      allocate( pscx( 1-halo:blki+halo, 1-halo:blkj+halo) )
      allocate( pscy( 1-halo:blki+halo, 1-halo:blkj+halo) )
      allocate( blk_conv( 1-halo:blki+halo, 1-halo:blkj+halo, 1:MEIG) )
      allocate( blk_ubaro( 1-halo:blki+halo, 1-halo:blkj+halo, 1:3) )
      allocate( blk_vbaro( 1-halo:blki+halo, 1-halo:blkj+halo, 1:3) )
      allocate( blk_pbot( 1-halo:blki+halo, 1-halo:blkj+halo, 
     &                    1:MEIG, 1:3) )
      allocate( glb_conv( 1:idm, 1:jdm, 1:MEIG ) )

      data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/
      
      ! Initialize global conversion field to 0.
      glb_conv(:,:,:) = 0.0

      do jblk=jblks,jblke
        do iblk=iblks,iblke
          ! Initialize blk/tile variables.
          pscx(:,:) = 0.0
          pscy(:,:) = 0.0
          z(:,:) = 0.0
          
          blk_ubaro(:,:,:) = 0.0
          blk_vbaro(:,:,:) = 0.0
          blk_pbot(:,:,:,:) = 0.0
          blk_conv(:,:,:) = 0.0

          ! Read input files / assign values to 2D variables.
          infileZ = '/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &       //'griddata/depth_221_blk_'//numstr(jblk)//'_'//
     &         numstr(iblk)//'.BinF'
          open(unit=11, file=infileZ, form='unformatted',
     &            status='OLD', access='sequential',action='READ')
          write(6,*) 'Reading ',infileZ
          call flush(6)
          read(11) z(:,:)
          close(11)

          infileDX= '/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &       //'griddata/pscx_221_blk_'//numstr(jblk)//'_'//
     &         numstr(iblk)//'.BinF'
          open(unit=12, file=infileDX, form='unformatted',
     &            status='OLD', access='sequential',action='READ')
          write(6,*) 'Reading ',infileDX
          call flush(6)
          read(12) pscx(:,:)
          close(12)

          infileDY= '/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &       //'griddata/pscy_221_blk_'//numstr(jblk)//'_'//
     &         numstr(iblk)//'.BinF'
          open(unit=13, file=infileDY, form='unformatted',
     &            status='OLD', access='sequential',action='READ')
          write(6,*) 'Reading ',infileDY
          call flush(6)
          read(13) pscy(:,:)
          close(13)
          
          infileUbt='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &       'HarmonicFit/'//'chf_ubaro_'//runid//
     &       '_blk_'//numstr(jblk)//'_'//
     &          numstr(iblk)//'.BinF'
           write(6,*) 'Reading ',infileUbt
           call flush(6)
           open(unit=14, file=infileUbt, form='unformatted',
     &            status='OLD', access='sequential',action='READ')
          write(6,*) 'Reading ',infileUbt
          call flush(6)
           read(14) blk_ubaro(:,:,:)
           close(14)

          infileVbt='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &       'HarmonicFit/'//'chf_vbaro_'//runid//
     &       '_blk_'//numstr(jblk)//'_'//
     &          numstr(iblk)//'.BinF'
           write(6,*) 'Reading ',infileVbt
           call flush(6)
           open(unit=15, file=infileVbt, form='unformatted',
     &            status='OLD', access='sequential',action='READ')
          write(6,*) 'Reading ',infileVbt
          call flush(6)
           read(15) blk_vbaro(:,:,:)
           close(15)

          infilePbot='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &       'HarmonicFit/'//'chf_pbot_'//runid//
     &       '_blk_'//numstr(jblk)//'_'//
     &          numstr(iblk)//'.BinF'
           write(6,*) 'Reading ',infilePbot
           call flush(6)
           open(unit=16, file=infilePbot, form='unformatted',
     &            status='OLD', access='sequential',action='READ')
          write(6,*) 'Reading ',infilePbot
          call flush(6)
           read(16) blk_pbot(:,:,:,:)
           close(16)

          ! Loop over values within tile/blk.
          ! Interpolate u, v to point C.
          ! Interpolate
          do i=1,blki
          do j=1,blkj
          do l=1,MEIG
 
            dhdx_W = ( z(i,j) - z(i-1,j) ) / 
     &                (0.5*(pscx(i,j) + pscx(i-1,j)))
            dhdx_E = ( z(i+1,j)-z(i,j) ) /  
     &                (0.5*(pscx(i+1,j) + pscx(i,j) ))

            dhdy_S = ( z(i,j) - z(i,j-1) ) / 
     &                (0.5*(pscy(i,j) + pscy(i,j-1)))
            dhdy_N = ( z(i,j+1)-z(i,j) )  / 
     &                (0.5*(pscy(i,j+1)+pscy(i,j)))

            u_conj_C = cmplx( 0.5*(blk_ubaro(i,j,2)+blk_ubaro(i+1,j,2)),
     &                      -0.5*(blk_ubaro(i,j,3)+blk_ubaro(i+1,j,3)) )
            v_conj_C = cmplx( 0.5*(blk_vbaro(i,j,2)+blk_vbaro(i,j+1,2)),
     &                      -0.5*(blk_vbaro(i,j,3)+blk_vbaro(i,j+1,3)) )

            p_botn_C = cmplx( blk_pbot(i,j,l,2), blk_pbot(i,j,l,3) )

            blk_conv(i,j,l) = real(
     &        0.5*(dhdx_W + dhdx_E)* 0.5*u_conj_C*p_botn_C
     &      + 0.5*(dhdy_S + dhdy_N)* 0.5*v_conj_C*p_botn_C )

            ! Store blk_conv into glb_conv.
            offset_i = blki*(iblk-1)
            offset_j = blkj*(jblk-1)
            glb_conv(offset_i+i, offset_j+j, l) = 
     &             blk_conv(i,j,l)
          enddo !l
          enddo !j
          enddo !l
        enddo !iblk
      enddo !jblk

      ! Write out the global conversion variable mode-by-mode.
      do l=1,MEIG
        outfil= '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &        'modes_global/'//'glb_conv_1z_'//runid//
     &        '_mode'//numstr(l)//'.BinF'
         write(6,*) 'Writing global bt conversion to ',outfil
         call flush(6)
         open(unit=99,file=outfil,status='replace',
     &   form='unformatted',access='sequential')
         write(99) glb_conv(:,:,l)
         close(99)
      enddo !l

      write(6,*) 'Finished ...'
      call flush(6)

      stop
      END PROGRAM CONVERSION
