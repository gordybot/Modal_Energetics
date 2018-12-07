      PROGRAM PSURF

      ! Compute the surface pressure and displacement from the modal pressure amplitudes.
      ! P_surf = rho_0 eta_0 integral( Nsqr(z) Phi(z) dz )      (Z.Zhao Appendix A).
      ! rho_0 is a reference density. eta_0 is the modal displacement amplitude.
      ! Integral from -H to 0. Nsqr is buoyancy frequency squared. Phi(z) is W, the w eigenvector.
      ! From Gerkema's book, I think eta_0 = 1/( rho_* (omega^2 - f^2) n*pi / (k_n^2 H) )
      ! Using EVAL = k_n**2/(omm2**2-f**2), we get
      !    eta0 = p0 / rho0 * H / (n PI) * EVAL.
      ! Surface displacement a = P_surf / (rho_0 g).
      IMPLICIT NONE 

      ! Variables I will need:
      character*3 runid
      character*2 numstr(90)
      character outfile*240, infileNsqr*240, infileEvec*240
      character infileP*240, infileB*240, infileEval*240
      real f0, omm2, PI, g, rho0, dz
      integer idm, jdm, kdm, mbdy, MEIG, blki, blkj, iblk, jblk
      integer iblks, iblke, jblks, jblke
      integer offset_i, offset_j, i,j,k,l, cc

      real, allocatable, dimension(:) :: eta0
      real, allocatable, dimension(:) :: P_surf
      real, allocatable, dimension(:,:,:) :: blk_Nsqr
      real, allocatable, dimension(:,:,:,:) :: blk_Evec
      real, allocatable, dimension(:,:,:) :: blk_Eval
      real, allocatable, dimension(:,:,:,:) :: modeP
      integer, allocatable, dimension(:,:,:) :: blk_boti
      real, allocatable, dimension(:,:,:,:) :: glb_psurf ! i,j,1:2, l

      ! Define some parameters.
      parameter(iblks=1,iblke=60,jblks=1,jblke=35)
      parameter(runid='221')
      parameter(dz=25.0, g=9.81, rho0=1034.0)
      parameter(idm=9000, jdm=7055, kdm=240, mbdy=3, MEIG=5)
      parameter( blki=150, blkj=200 )
      PI = 3.1415926535897932384
      omm2 = 2*PI / (12.42*3600)
      f0 = 2*(2*PI/(24*3600) )

      data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/

      ! Allocate memory for arrays.
      allocate( eta0(1:2) )
      allocate( P_surf(1:2) )

      allocate( blk_Nsqr( 1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 0:kdm) )
      allocate( blk_Evec( 1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy,
     &              1:kdm, 1:MEIG))
      allocate( blk_Eval(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:MEIG) )
      allocate( blk_boti( 1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:2) )
      allocate( modeP(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:2, 1:MEIG))
      allocate( glb_psurf(idm, jdm, 1:2, 1:MEIG) )

      ! Initialize global variables.
      glb_psurf(:,:,:,:) = 0.0

      do jblk=jblks,jblke
      do iblk=iblks,iblke
        offset_i = blki*(iblk-1)
        offset_j = blkj*(jblk-1)

        blk_boti(:,:,:) = 0
        blk_Nsqr(:,:,:) = 0.0
        blk_Evec(:,:,:,:) = 0.0
        blk_Eval(:,:,:) = 0.0
        modeP(:,:,:,:) = 0.0

        infileB = '/p/work1/mbui/hycom/GLBc0.04/expt_22.1/modes/'
     &         //'BOTI_1z_221_'//numstr(jblk)//'_'//numstr(iblk)
     &         //'.BinF'
        write(6,*) 'reading ',infileB
        call flush(6)
        open(unit=10,file=infileB,status='old',form='unformatted')
        ! Read blk_boti
        read(10) blk_boti(:,:,:)   
        close(10)
        write(6,*) 'Finished reading ',infileB
        call flush(6)

        infileNsqr = '/p/work1/mbui/hycom/GLBc0.04/expt_22.1/zN2avg/'
     &      //'N2avg_1z_221_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
        write(6,*) 'reading ',infileNsqr
        call flush(6)
        open(unit=11,file=infileNsqr,status='old',form='unformatted')
        ! Read blk_Nsqr
        read(11) blk_Nsqr(:,:,:)
        write(6,*) 'Finished reading ',infileNsqr
        call flush(6) 

        infileEvec = '/p/work1/mbui/hycom/GLBc0.04/expt_22.1/modes/'
     &     //'EVEC_1z_221_'//numstr(jblk)//'_'//numstr(iblk)//'.BinF'
        write(6,*) 'reading ',infileEvec
        call flush(6)
        open(unit=12,file=infileEvec,status='old',form='unformatted')
        ! Read blk_Evec
        do l=1,MEIG
          read(12) blk_Evec(:,:,:,l)
        enddo !loop over l
        close(12)
 
        infileP = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &        //'modes/phfvm_221_blk_'//numstr(jblk)//'_'//numstr(iblk)
     &        //'_apr.BinF'
        write(6,*) 'reading ',infileP
        call flush(6)
        open(unit=13,file=infilep,status='old',form='unformatted')
        ! Read modeP
        read(13) modeP(:,:,:,:)
        close(13)

        ! Load plat to compute Coriolis frequency.
        ! f= 2*(2*PI/(24*3600))*SIN(lat_i*PI/180)
        ! use f0 = 2*(2*PI/(24*3600))
        infileEval='/p/work1/mbui/hycom/GLBc0.04/expt_22.1/'
     &     //'modes/EVAL_1z_221_'//numstr(jblk)//'_'//
     &       numstr(iblk)//'.BinF'
        write(6,*) 'Reading ',infileEval
        call flush(6)
        open(unit=14, file=infileEval, form='unformatted',
     &            status='OLD', access='sequential',action='READ')
        do l=1,MEIG
          read(14) blk_Eval(:,:,l)
        enddo !l.
        close(14)

        ! Initialize P_surf 
        ! Compute the integral
        !P_surf = rho_0 eta_0 integral( Nsqr(z) Phi(z) dz )      (Z.Zhao Appendix A).
        do i=1,blki
        do j=1,blkj
          cc = blk_boti(i,j,1) ! Number of deepest good bin.
          P_surf(:) = 0.0

          if (cc.GT.1) then
          do l=1,MEIG
            eta0(:) = 0.0
            P_surf(:) = 0.0
           !eta_0 = 1/( rho_* (omega^2 - f^2) n*pi / (k_n^2 H) )
           ! TODO: Does rho* = rho0? Ask Gerkema's book...
            ! TODO: Check indices on modeP...
            eta0(1) = modeP(i,j,1,l) / rho0 
     &            * (cc*dz) / (l * PI) * blk_Eval(i,j,l)
            eta0(2) = modeP(i,j,2,l) / rho0
     &            * (cc*dz) / (l * PI) * blk_Eval(i,j,l)

            ! Integral over depth of Nsqr(z) Evec(z) dz.
            P_surf =  rho0 * eta0 
     &        * sum( blk_Evec(i,j,1:cc,l) * blk_Nsqr(i,j,1:cc) ) * dz 

            ! Store in global array.
            glb_psurf( offset_i + i, offset_j + j, 1:2, l) = P_surf
          enddo !l             
          endif ! Any good cells / cc.GT.1.

        enddo !j
        enddo !i

      enddo !iblk
      enddo !jblk             
   

       outfile ='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &           //'modes_global/'
     &           //'glb_Psurf.BinF'
       write(6,*) 'Writing glb_psurf to ',outfile
       call flush(6)
       open(unit=98,file=outfile,status='replace',
     &          form='unformatted',access='sequential')
       do l=1,MEIG
         write(98) glb_psurf(:,:,:,l)
       enddo !l
       close(98)

       outfile ='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &           //'modes_global/'
     &           //'glb_surfH.BinF'
       write(6,*) 'Writing glb_surfH to ',outfile
       call flush(6)
       open(unit=98,file=outfile,status='replace',
     &          form='unformatted',access='sequential')
       do l=1,MEIG
         write(98) glb_psurf(:,:,:,l)/(rho0*g)
       enddo !l
       close(98)


       STOP
       END PROGRAM PSURF
