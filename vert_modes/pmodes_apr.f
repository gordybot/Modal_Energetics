      PROGRAM PRESSURE_MODES
      ! Fixed ERROR on line 202 -- was not normalizing sign
      ! to match at bottom.

      ! New attempt at making the pressure modes work.
      ! Hawk and Pony, 1/17/18, GS.

      ! Remove each mode after it is fit.
      ! 1/23/18, GS.

      IMPLICIT NONE

      ! Declare and type variables.
      integer i,j,k,l,q,cc, mbdy
      integer idm, jdm, kdm, blki, blkj
      integer iblk, iblks, iblke
      integer jblk, jblks, jblke
      integer IL, IU, MEIG
      character*3 runid
      character infileE*240, infileB*240
      character infileP*240, outfile*240
      character outfilePbot*240, outfilePsurf*240
      character*2 numstr(90)
      integer signUvec
      real dz, grav, sum1
      real mean_dWdz, var_dWdz, meanA
      real mean_Pcell, H, tempSum1, tempSum2

      ! Declare array variables.
      integer, allocatable, dimension(:,:,:)  :: blk_boti
      real, allocatable, dimension(:) :: W
      real, allocatable, dimension(:) :: dWdz
      real, allocatable, dimension(:,:,:,:) :: blk_Evec
      real, allocatable, dimension(:,:,:,:) :: blk_sig
      real, allocatable, dimension(:) :: Uvec
      real, allocatable, dimension(:,:,:) :: ebot
      real, allocatable, dimension(:,:,:) :: etop
      real, allocatable, dimension(:) :: Pface
      real, allocatable, dimension(:) :: Pcell

      real, allocatable, dimension(:,:) :: x
      real, allocatable, dimension(:,:) :: A
      real, allocatable, dimension(:,:) :: B
      real, allocatable, dimension(:) :: WORK
      integer LWORK, INFO

      real, allocatable, dimension(:,:,:,:) :: modeP
      real, allocatable, dimension(:,:,:,:) :: pbot
      real, allocatable, dimension(:,:,:,:) :: psurf

      real, allocatable, dimension(:) :: S
      integer rank

      ! Assign parameters.
       parameter(idm=9000,jdm=7055,kdm=240)
       parameter(blki=150,blkj=200,mbdy=3)
       parameter(IL=1,IU=5)
       parameter(dz=25.0)
       parameter(grav=9.81)

       ! Define some basic variables.
       MEIG=IU-IL+1
       data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/

       LWORK=4*4*kdm*MEIG
 
      ! Allocate variables.
      allocate( blk_sig(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:kdm,
     & 1:3) )
      allocate( blk_boti(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy, 1:2) )
      allocate( W(0:kdm) )
      allocate( Pface(0:kdm) )
      allocate( Pcell(1:kdm) )
      allocate( dWdz(1:kdm) )
      allocate( Uvec(1:kdm) )
      allocate( A(1:kdm, 1:MEIG) )
      allocate( B(1:kdm, 1:2) )
      allocate( modeP(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:2, 
     &          1:MEIG) )
      allocate( pbot(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:2, 1:MEIG) )
      allocate( psurf(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy, 1:2, 1:MEIG) )
      allocate( blk_Evec(1-mbdy:blki+mbdy,1-mbdy:blkj+mbdy, 1:kdm,
     &           1:MEIG) )
      allocate( ebot(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:MEIG) )
      allocate( etop(1-mbdy:blki+mbdy, 1-mbdy:blkj+mbdy, 1:MEIG) )
      allocate( WORK(1:LWORK) )

      ! NEW! TODO!
      ! allocate( S(1:kdm) )
      allocate( x(1:MEIG, 1:2) )

      ! Read from input file.
       read(*,*) jblks
       read(*,*) jblke
       read(*,*) iblks
       read(*,*) iblke
       read(*,*) runid

       write(6,'(a,a3)'),'runid: ',runid
       call flush(6)

       write(6,'(a,4i3)') 'jblks jblke iblks iblke: ',jblks, jblke,
     & iblks, iblke
       call flush(6)

      ! Loop over iblk, jblk
       do jblk=jblks,jblke
       do iblk=iblks,iblke
         write(6,*) 'iblk, jblk = ',iblk,jblk
         call flush(6)
        
         INFO = -10
       ! Initialize blk variables:
         !Initialize blk_Evec
         blk_Evec(:,:,:,:)=0.0
         ! Initialize blk_boti
         blk_boti(:,:,:)=0 ! integer.
         ! Initialize sigma.
         blk_sig(:,:,:,:) = 0.0 
         ! Reset ebot to 0.0
         ebot(:,:,:) = 0.0
         etop(:,:,:) = 0.0
         !modeP, pbot
         pbot(:,:,:,:) = 0.0
         psurf(:,:,:,:) = 0.0
         modeP(:,:,:,:) = 0.0 

      ! READ a bunch of FILES.
      ! Load blk_boti, blk_evec, blk_chf_sigma
         ! Read in the eigenvectors.       
         infileE='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/modes/'
     &   //'EVEC_1z_'//runid//'_'//numstr(jblk)//'_'//numstr(iblk)//
     &   '.BinF'
         write(6,*) 'reading ',infileE
         call flush(6)
         open(unit=11,file=infileE,status='old',form='unformatted')
         !Read blk_Evec in same way it was written out.
         do l=1,MEIG
           read(11) blk_Evec(:,:,:,l)
         enddo !loop over l
         close(11)

         ! Read in the bottom bin locations.
         infileB='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/modes/'
     &   //'BOTI_1z_'//runid//'_'//numstr(jblk)//'_'//numstr(iblk)//
     &   '.BinF'
         write(6,*) 'reading ',infileB
         call flush(6)
         open(unit=12,file=infileB,status='old',form='unformatted')
         ! Read blk_boti
         read(12) blk_boti(:,:,:)
         close(12)
         write(6,*) 'Finished reading ',infileB
         call flush(6)

         ! Open the file for density.
!         infileP='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/sig/'
         infileP='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &       //'HarmonicFit/chf_'
     &   //'sig_'//runid//'_blk_'//numstr(jblk)//'_'//numstr(iblk)//
     &   '_apr.BinF'
         write(6,*) 'reading ',infileP
         call flush(6)
         open(unit=15,file=infileP,status='old',form='unformatted')
         read(15) blk_sig(:,:,:,:)
         close(15)

      ! Loop over i,j.
      do i=1-mbdy,blki+mbdy
        do j=1-mbdy,blkj+mbdy
          ! Reinitialize A to 0
          A(:,:) = 0.0

          cc=blk_boti(i,j,1) ! number of deepest good bin.
          ! Check to see that there ARE any good cells:
          if (cc.GT.1) then

           ! Normalize the eigenvectors.
           do l=1,MEIG
            ! Initialize W, dWdz.
            W(:)=0.0
            dWdz(:) = 0.0
            ! Define dW/dz ~ U,V,P eigenvector.
            do k=1,cc
              W(k) = blk_Evec(i,j,k,l)  ! face-centered
              dWdz(k) = ( W(k)-W(k-1) ) / dz ! cell-centered.
            enddo
            ! We will normalize so that the integral of Uvec**2 dz == H.
            sum1 = sum( dWdz(1:cc)**2 ) * 1.0/cc

            ! Define  Uvec, to be the normalized eigenvector for U
            Uvec(:) = 0.0
            ! Possibly prevent some Infinity values.
            if (sum1.GT.0.0) then
              do k=1,cc
                !Uvec(k) = ( dWdz(k) - mean_dWdz ) / (var_dWdz**0.5)
                Uvec(k) = dWdz(k) / (sum1**0.5)
              enddo !k
            endif

            signUvec = 1
            if (Uvec(cc).LT.0.0) then
!            if (A(cc,l).LT.0) then  ! ERROR -- A=0.0
              signUvec = -1
            endif

             ! Moved this below make sure meanA=0.
!            A(:,l) = signUvec * Uvec(:) 
!            A(:,l) = Uvec(:) 
!            ebot(i,j,l) = A(cc,l)
                       
            A(:,l) = signUvec * Uvec(:)
            ! Make sure mean(A(:,l))==0.0
            meanA = sum(A(:,l))*1.0/cc
            do k=1,cc
              A(k,l) = A(k,l) - meanA
            enddo !k

             ! Moved from 8 lines up. 2.13.18.GS
!            A(:,l) = Uvec(:)
!             write(6,*) 'A(cc,l), i,j,l ', A(cc,l),' ',i,' ',j,' ',l
!             call flush(6)
             ebot(i,j,l) = A(cc,l)
             etop(i,j,l) = A( 1,l)
          enddo !l loop over modes.

          ! Compute perturbation pressure and prepare B matrix
          ! A*x = B -- B has the data in it.
          ! q loops over cosine and sine terms in the harmonic fit
          ! input in chf_..sig_blk_JJ_II
          B(:,:) = 0.0
          do q=2,3
            Pface(:) = 0.0 ! P has indices (0:kdm)
            ! Convert density perturbation to pressure perturbation.
            do k=1,cc
              Pface(k) = Pface(k-1) + blk_sig(i,j,k,q) * grav * dz
            enddo !k 
 
            Pcell(:) = 0.0 ! Pcell [1:kdm] is cell-centered
            ! Interpolate face-centered P to cell-centered.
            do k=1,cc
              Pcell(k) = 0.5 * ( Pface(k) + Pface(k-1) )
            enddo !k

            ! Remove depth-mean of pressure and assign to B matrix.
            mean_Pcell = 1.0/cc * sum( Pcell(1:cc) )
            do k=1,cc
              B(k,q-1) = Pcell(k) - mean_Pcell
            enddo !k
          enddo !q

          ! Call subroutine SGELS to solve A*x=B.
              ! Call subroutine SGELS('N',cc,MEIG,maxobstot,
              ! A, 240, BU, 240, WORK(1:LWORK), LWORK, INFO)
              ! where LWORK >= MN+max(MN,NRHS)*NB
          ! * SGELS modifies A, B, WORK, and INFO.
!          call SGELS('N',kdm, MEIG, 2, A, kdm, B, kdm,
!     &               WORK(1:LWORK), LWORK, INFO)

          ! Attempt a different subroutine to see what happens.
          ! args: M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
          !       WORK, LWORK, INFO
          ! A is M by N. LDA = max(1,M)
          ! B is M by NRHS. LDB = max(1,max(M,N))
          ! S (out) is real, (min(M,N)) kdm
          ! RCOND (in) set to -1 ==> use machine precision.
          ! rank is an integer (out)
!          S(:) = 0.0
!          rank = -10
!          WORK(:) = 0.0
!          call SGELSS(kdm, MEIG, 2, A, kdm, B, kdm, S, -1, rank,
!     &                WORK(1:LWORK), LWORK, INFO)
          ! INFO == 0 on a successful exit.
          ! If something went wrong, make a note of it.
          !if (INFO.NE.0) then
          !  write(6,*) 'INFO, i,j: ', INFO, i, j
          !  call flush(6)
          !endif         

          ! Project normal modes onto B.
          ! x = 1/H * A dot B.
          ! ^^ Sort of ^^. Technically, I normalized so that 
          ! x = 1/cc * A dot B, I believe... TODO: 1.31.18.GS
          x(:,:) = 0.0
          ! Moved the do l=1,MEIG down below the vertical mean.
!          do l=1,MEIG
            ! Remove vertical mean of B(:,1:2)
           !! Except that I already removed the vertical mean -- 
           !! Either remove it again after every modal fit, or 
           !! don't worry about it!   2.9.18.GS.
!            tempSum1 = 1.0/cc * sum( B(1:cc,1) )
!            tempSum2 = 1.0/cc * sum( B(1:cc,2) )
!            tempVal = 0.0
!            do k=1,cc
!              tempVal = B(k,1) - tempSum1
!              B(k,1) = tempVal
!              tempVal = B(k,2) - tempSum2
!              B(k,2) = tempVal
!            enddo !k
          do l=1,MEIG
            tempSum1 = 0.0
            tempSum2 = 0.0
            ! H = cc*dz
            ! Compute the dot product of eigenvector with data.
            do k=1,kdm
              tempSum1 = tempSum1 + 1.0/cc * A(k,l) * B(k,1)
              tempSum2 = tempSum2 + 1.0/cc * A(k,l) * B(k,2)
            enddo !k
            x(l,1) = tempSum1
            x(l,2) = tempSum2
            ! Remove last mode from the data.
            do k=1,kdm
              B(k,1) = B(k,1) - tempSum1*A(k,l)
              B(k,2) = B(k,2) - tempSum2*A(k,l)
            enddo !k
          enddo !l
          
          ! Store the first MEIG rows of B -- these are solutions.
          do l=1,MEIG
            do q=1,2
              modeP(i,j,q,l) = x(l,q)
              pbot(i,j,q,l)  = x(l,q)*ebot(i,j,l)
              psurf(i,j,q,l) = x(l,q)*etop(i,j,l)
            enddo !q: cos, sin
          enddo !l: mode number.

         ! Continue check on any good cells.
         else if (cc.LT.2) then
           modeP(i,j,:,:) = 0.0
           pbot(i,j,:,:) = 0.0
           psurf(i,j,:,:) = 0.0
         endif !  
       enddo !j
      enddo !i

      ! Write pressure mode amplitudes to outfile.
      outfile='/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/modes/'
     &   //'phfvm_'//runid//'_blk_'//numstr(jblk)//'_'//numstr(iblk)
     &   //'_apr.BinF'
      write(6,*) 'Writing modeU to ',outfile
      call flush(6)
      open(unit=98,file=outfile,status='replace',
     &       form='unformatted',access='sequential')
      write(98) modeP(:,:,:,:)
      close(98)

      outfilePbot = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/modes/'
     &    //'pbot_hfvm_'//runid//'_blk_'//numstr(jblk)//'_'
     &    //numstr(iblk)//'_apr.BinF'
      write(6,*) 'Writing pbot to ',outfilePbot
      call flush(6)
      open(unit=97,file=outfilePbot,status='replace',
     &     form='unformatted',access='sequential')
      write(97) pbot(:,:,:,:)
      close(97)

      outfilePsurf = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/modes/'
     &    //'psurf_hfvm_'//runid//'_blk_'//numstr(jblk)//'_'
     &    //numstr(iblk)//'_may.BinF'
      write(6,*) 'Writing psurf to ',outfilePsurf
      call flush(6)
      open(unit=96,file=outfilePsurf,status='replace',
     &     form='unformatted',access='sequential')
      write(96) psurf(:,:,:,:)
      close(96)      

         write(6,*) 'iblk / jblk / MIN P / MAX P: '
         call flush(6)
         write(6,*) iblk,' ',jblk,' ',MINVAL(modeP),' ', MAXVAL(modeP)
         call flush(6)

         write(6,*) 'iblk / jblk / MIN pbot / MAX pbot: '
         call flush(6)
         write(6,*) iblk,' ',jblk,' ',MINVAL(pbot),' ', MAXVAL(pbot)
         call flush(6)

         write(6,*) 'MIN Psurf / MAX Psurf:'
         call flush(6)
         write(6,*) MINVAL(psurf),' ',MAXVAL(psurf)
         call flush(6)
      enddo !iblk
      enddo !jblk

      write(6,*) 'Finished ...'
      call flush(6)

      STOP
      END PROGRAM PRESSURE_MODES

