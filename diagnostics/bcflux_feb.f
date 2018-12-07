       program bcflux
       !!! input modeX has format (i j t=1:3, l=1:MEIG)
       !!! Use the complex harmonic fits to compute baroclinic fluxes.
       implicit none

       ! Declare and type variables.
       character*3 runid
       character*2 numstr(90)
       character infileP*240, infileU*240, infileV*240
       character outfileF*240, infile_depth*240
       logical isThere
       integer i, j, l, iblk, jblk, mbdy
       integer idm, jdm, blki, blkj
       integer MEIG
       integer offset_i, offset_j
       real z
       real, allocatable, dimension(:,:) :: depth
       
       real, allocatable, dimension(:,:,:,:) :: blk_modeU
       real, allocatable, dimension(:,:,:,:) :: blk_modeV
       real, allocatable, dimension(:,:,:,:) :: blk_modeP

       complex, allocatable, dimension(:,:,:) :: cblk_fluxU
       complex, allocatable, dimension(:,:,:) :: cblk_fluxV
       complex, allocatable, dimension(:,:,:,:) :: fluxUV

       ! Assign parameters, define useful constants.
       parameter(idm=9000, jdm=7055, mbdy=3)
       parameter(blki=150 ,blkj= 200, MEIG = 5)
       parameter(runid='221')
       data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/


       ! Allocate variables.
       allocate( blk_modeP(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,
     &                    1:2, 1:MEIG) )
       allocate( blk_modeU(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,
     &                    1:2, 1:MEIG) )
       allocate( blk_modeV(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,
     &                    1:2, 1:MEIG) )
       allocate( cblk_fluxU(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,
     &                    1:MEIG ) )
       allocate( cblk_fluxV(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy,
     &                    1:MEIG ) )
       allocate( fluxUV(1:idm, 1:jdm, 1:2, 1:MEIG) )
       allocate( depth(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy) )

       ! Initialize global variables.
       fluxUV(:,:,:,:) = 0.0

       ! Loop over blocks.
       do jblk=1,35
         do iblk=1,60

       ! Initialize blk variables.
         offset_i = blki * (iblk-1)
         offset_j = blkj * (jblk-1)
         cblk_fluxU(:,:,:) = 0.0
         cblk_fluxV(:,:,:) = 0.0
         depth(:,:) = 0.0

       ! Read in blk_modeU, V, P.
          !!! Read the pressure mode amplitudes.
           infileP='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
!     &        'HarmonicFit/chf_modeP_'
!     &        'HarmonicFit/hfvmP_'
     &        'modes/phfvm_'
     &        //runid//'_blk_'//numstr(jblk)//
     &        '_'//numstr(iblk)//'_apr.BinF'
           ! Read one block.
           write(6,*) 'reading ',infileP
           call flush(6)
           ! Initialize blk_mHFp. Indices are (i,j,l,[mean amp pha])
           blk_modeP(:,:,:,:) = 0.0
         inquire(file=infileP, EXIST=isThere)
           open(unit=11, file=infileP,status='old',form='unformatted')
           read(11) blk_modeP(:,:,:,:)
           close(11)

           ! Read u_mode
          !!! Read the U mode amplitudes.
           infileU='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
!     &        'HarmonicFit/chf_modeU_'
     &        'modes/uhfvm_'
     &        //runid//'_blk_'//numstr(jblk)//
     &        '_'//numstr(iblk)//'_may.BinF'
           ! Read one block.
           write(6,*) 'reading ',infileU
           call flush(6)
           ! Initialize blk_mHFu. Indices are (i,j,l,[mean amp pha])
           blk_modeU(:,:,:,:) = 0.0
         if (isThere) then
             inquire(file=infileU, EXIST=isThere)
         endif
           open(unit=12, file=infileU,status='old',form='unformatted')
           read(12) blk_modeU(:,:,:,:)
           close(12)

           ! Read v_mode
          !!! Read the V mode amplitudes.
           infileV='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
!     &        'HarmonicFit/chf_modeV_'
     &        'modes/vhfvm_'
     &        //runid//'_blk_'//numstr(jblk)//
     &        '_'//numstr(iblk)//'_may.BinF'
           ! Read one block.
           write(6,*) 'reading ',infileV
           call flush(6)
           ! Initialize blk_mHFv. Indices are (i,j,l,[mean amp pha])
           blk_modeV(:,:,:,:) = 0.0
         if (isThere) then
             inquire(file=infileV, EXIST=isThere)
         endif
           open(unit=13, file=infileV,status='old',form='unformatted')
           read(13) blk_modeV(:,:,:,:)
           close(13)
 
          !!! Read depth
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

       ! If U, V, and P modes exist.
         if (isThere) then

       ! Compute fluxes = complex amplitudes of harmonic fits for 
       ! modeU * modeP * depth.
           do i=1,blki
           do j=1,blkj
             z = depth(i,j)
           do l=1,MEIG
             if (z.EQ.0) then
               cblk_fluxU(i,j,l) = 0.0
             else
        ! Interpolate U* and V* to Center of cell.
        ! Change third index from [2,3] to [1,2].
        ! EDIT / TODO: 2.2.18.GS
        cblk_fluxU(i,j,l) = 0.5* real(
     &     0.5* ( cmplx(blk_modeU(i,j,1,l), -1*blk_modeU(i,j,2,l))
     &           +cmplx(blk_modeU(i+1,j,1,l),-1*blk_modeU(i+1,j,2,l)))
     &       * cmplx( blk_modeP(i,j,1,l), blk_modeP(i,j,2,l) ) * z )
        cblk_fluxV(i,j,l) = 0.5* real(
     &     0.5*( cmplx(blk_modeV(i,j,1,l), -1*blk_modeV(i,j,2,l))
     &          +cmplx(blk_modeV(i,j+1,1,l),-1*blk_modeV(i,j+1,2,l)))
     &       * cmplx( blk_modeP(i,j,1,l), blk_modeP(i,j,2,l) ) * z )
             endif

         ! Assign to storage variable.
       fluxUV(offset_i+i, offset_j+j,1,l) = cblk_fluxU(i,j,l)
       fluxUV(offset_i+i, offset_j+j,2,l) = cblk_fluxV(i,j,l)

         ! End loop over i,j within block.
           enddo !l
           enddo !j
           enddo !i

         endif ! If isThere

       ! End loop over blocks.
         enddo !iblk
       enddo !jblk


       do l=1,MEIG
         ! Write out storage variable.
         outfileF = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'//
     &        'modes_global/'//'bcflux_'//runid//'_mode'
     &      //numstr(l)//'_may.BinF'
         write(6,*) 'Writing global baroclinic fluxes to ',outfileF
         call flush(6)

         open(unit=99,file=outfileF,status='replace',
     &   form='unformatted',access='sequential')
         write(99) real(fluxUV(:,:,1,l))
         write(99) real(fluxUV(:,:,2,l))
         close(99)
       enddo !l

       !!! Done.       
       write(6,*) 'Finished ...'
       call flush(6)

       end program bcflux
