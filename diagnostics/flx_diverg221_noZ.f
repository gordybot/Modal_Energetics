      PROGRAM FLUX_DIVERG

      IMPLICIT NONE
      ! Load fluxes and pscy, pscx,
      ! Then compute the flux divergence.

      ! Declare and type variables.
      character infilF*240, infilePSCX*240, infilePSCY*240
      character infileZ*240, outfil*240
      integer MEIG, i,j,l, mbdy, blki, blkj, iblk, jblk
      integer iblks, iblke, jblks, jblke, idm, jdm
      integer offset_i, offset_j
      character*3 runid
      character*2 numstr(90)
      real, allocatable, dimension(:,:) :: glb_flxu
      real, allocatable, dimension(:,:) :: glb_flxv
      real, allocatable, dimension(:,:) :: glb_flx_div
      real, allocatable, dimension(:,:) :: blk_pscx
      real, allocatable, dimension(:,:) :: blk_pscy
      !real, allocatable, dimension(:,:) :: blk_depth

      ! Parameters
      parameter(runid='221')
      parameter(idm=9000,jdm=7055,blki=150,blkj=200,mbdy=3)
      parameter(iblks=1,iblke=60,jblks=1,jblke=35)
      parameter(MEIG=5)

      data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/

      allocate( blk_pscx(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy) )
      allocate( blk_pscy(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy) )
      allocate( glb_flxu(1-mbdy:idm+mbdy, 1-mbdy:jdm+mbdy) )
      allocate( glb_flxv(1-mbdy:idm+mbdy, 1-mbdy:jdm+mbdy) )
      allocate( glb_flx_div(1:idm,1:jdm) )
      !allocate( blk_depth(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy) )

      ! For each mode
      do l=1,MEIG
        ! Initialize global flux field.
        glb_flxu(:,:) = 0.0
        glb_flxv(:,:) = 0.0
        glb_flx_div(:,:) = 0.0

        ! Load global flux file for mode l.
        infilF = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &         'modes_global/bcflux_'//runid//'_mode'//numstr(l)//
     &         '_mar.BinF'

        open(unit=17,file=infilF,status='old',
     &              form='unformatted')
        read(17) glb_flxu(1:idm,1:jdm)
        read(17) glb_flxv(1:idm,1:jdm)
        close(17)

        ! Wraparound glb_flx
        glb_flxu(1-mbdy:0,1:jdm) = glb_flxu(idm+1-mbdy:idm,1:jdm)
        glb_flxu(idm+1-mbdy:idm,1:jdm) = glb_flxu(1:mbdy,1:jdm)

        glb_flxv(1-mbdy:0,1:jdm) = glb_flxv(idm+1-mbdy:idm,1:jdm)
        glb_flxv(idm+1-mbdy:idm,1:jdm) = glb_flxv(1:mbdy,1:jdm)

        do jblk=jblks,jblke
        do iblk=iblks,iblke
          ! Initialize blk fields.
          blk_pscy(:,:) = 0.0
          blk_pscx(:,:) = 0.0
          
          offset_i = blki*(iblk-1)
          offset_j = blkj*(jblk-1)

          write(6,*) 'iblk, jblk = ',iblk,' ',jblk
          call flush(6)

          ! Read in the pscx and pscy files.
          ! Reinitialize blk_pscx, blk_pscy
          blk_pscx(:,:)=0.0
          blk_pscy(:,:)=0.0
          !blk_depth(:,:)=0.0

          ! Read pscx, pscy from gridfiles.
          infilePSCX='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &   //'griddata/pscx_'//runid//'_blk_'//numstr(jblk)//'_'//
     &    numstr(iblk)//'.BinF'
         open(unit=15,file=infilePSCX,status='old',form='unformatted')
          read(15) blk_pscx(:,:)
          close(15)

          infilePSCY='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &  //'griddata/pscy_'//runid//'_blk_'//numstr(jblk)//'_'//
     &    numstr(iblk)//'.BinF'
         open(unit=16,file=infilePSCY,status='old',form='unformatted')
          read(16) blk_pscy(:,:)
          close(16)

!          infileZ='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
!     &   //'griddata/depth_'//runid//'_blk_'//numstr(jblk)//'_'//
!     &     numstr(iblk)//'.BinF'
!          open(unit=17,file=infileZ,status='old',form='unformatted')
!          read(17) blk_depth(:,:)
!          close(17)
         
          ! Compute flux divergence.
          do i=1,blki
          do j=1,blkj
            glb_flx_div( offset_i+i, offset_j + j) = 
     &        ( glb_flxu( offset_i+i+1, offset_j + j)
     &        - glb_flxu( offset_i+i-1, offset_j + j) )
     &        / ( blk_pscx( i, j) + 
     &            0.5*( blk_pscx( i-1, j)
     &                 +blk_pscx( i+1, j) ) ) ! * blk_depth(i,j)
     &      + ( glb_flxv( offset_i+i, offset_j + j+1)
     &        - glb_flxv( offset_i+i, offset_j + j-1) )
     &        / ( blk_pscy( i,j) +
     &            0.5*( blk_pscy( i, j+1)
     &                 +blk_pscy( i, j-1) ) ) ! * blk_depth(i,j)
          enddo !j
          enddo !i
        enddo !jblk
        enddo !iblk

        ! Write out global flux divergence for mode l.
        outfil = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &         //'modes_global/glb_flx_diverg_'//runid//'_mode'//
     &           numstr(l)//'_mar.BinF'
        open(unit=99,file=outfil,status='replace',
     &          form='unformatted', access='sequential')
        write(99) glb_flx_div(:,:) 
        close(99)

        write(6,*) 'Writing glb flux div. mod',numstr(l),' ',outfil
        call flush(6)
      enddo !l

      write(6,*) 'Finished...'
      call flush(6)

      STOP
      END PROGRAM FLUX_DIVERG
