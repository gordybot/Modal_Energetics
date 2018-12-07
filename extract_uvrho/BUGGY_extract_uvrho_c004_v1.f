	PROGRAM EXTRACT_UVRHO

	IMPLICIT NONE

!c---Used to extract baroclinic u, v, and rho
!c into tiles. Based on 
!c extract_vars_c008_v12.f AND mean_rho_c008_v1.f
!c GRS - 05/23/2017
!c
!c load global *.a files interpolated to 240 layers
!c /p/work1/mbui/hycom/GLBc0.08/expt_05.1/input3z/

	integer idm,jdm,kdm,blki,blkj,halocd,halo
	integer jblks,jblke !1,35
	integer iblks,iblke !1,60
	real flag
	real maxvar
!c glbc008
	parameter(idm=9000,jdm=7055,kdm=240)
	parameter(blki=150,blkj=200,halo=3)

	parameter(flag=2.0**100)

	character flnm*240,outfil*240,flnm2*240,flnm3*240,flnm4*240
	character flnm5*240
	character blank*5,runid*3,intid*3
	character*2 numstr(90)
	integer i,j,k,n2drec,nrecl,numfls,iblk,jblk,n,flip_mcb

	real, allocatable, dimension(:) :: w
	real, allocatable, dimension(:,:,:) :: sig
	real, allocatable, dimension(:,:,:) :: u
	real, allocatable, dimension(:,:,:) :: v

	data blank/'     '/

        data numstr/'01','02','03','04','05','06','07','08','09','10',
     &            '11','12','13','14','15','16','17','18','19','20',
     &            '21','22','23','24','25','26','27','28','29','30',
     &            '31','32','33','34','35','36','37','38','39','40',
     &            '41','42','43','44','45','46','47','48','49','50',
     &            '51','52','53','54','55','56','57','58','59','60',
     &            '61','62','63','64','65','66','67','68','69','70',
     &            '71','72','73','74','75','76','77','78','79','80',
     &            '81','82','83','84','85','86','87','88','89','90'/

	!n2drec = ((idm*jdm+4095)/4096)*4096
        n2drec = (((blki+6)*(blkj+6)+4095)/4096)*4096
!c
!c	initialize I/O buffer
!c
	!allocate( w(n2drec) )
	!allocate( sig( -halo:idm+halo, -halo:jdm+halo,kdm) )
	!allocate( u( -halo:idm+halo, -halo:jdm+halo,kdm) )
	!allocate( v( -halo:idm+halo, -halo:jdm+halo,kdm) )
	
        allocate( w(n2drec) )
        allocate( sig(1-halo:blki+halo, 1-halo:blkj+halo,kdm) )
        allocate( u( 1-halo:blki+halo,  1-halo:blkj+halo,kdm) )
        allocate( v( 1-halo:blki+halo,  1-halo:blkj+halo, kdm) )

	sig(:,:,:) = 0.0
!c
	inquire(iolength=nrecl) w
!c
! read beginning and end of blocks
	read(*,*) jblks
	read(*,*) jblke
	read(*,*) iblks
	read(*,*) iblke
	read(*,*) runid
	read(*,*) intid
	
!c	print out some of the info
	write (6,'(a,a3)' ) 'runid: ',runid
	call flush(6)

	write (6,'(a,a3)') 'intid: ',intid
	call flush(6)

	write (6,'(a)') blank
	call flush(6)
	write (6,'(a,4i3)' ) 'jblks jblke iblks iblke: ',jblks,jblke,
     &iblks,iblke
	call flush(6)
	read(*,*) numfls
	write (6,'(a)') blank
	call flush(6)
	write (6,'(a,i5)') 'number of time steps: ',numfls
	call flush(6)

! skip over the first files //depth, grid, tidal, cb
! GS-5/23/17 -- modified .com file to not include these lines.
!	read (*,'(a)') flnm2
!	read (*,'(a)') flnm3
!	read (*,'(a)') flnm4
!	read (*,'(a)') flnm5
	write(6,'(a)') blank

!c begin loop over files listed in *[blkj_blki].in file.
!c in order, they are listed *_r.a, *_u.a, *_v.a, so track the
!c modulus of n to know which variable to load.
 
	do n=1,numfls
	   
	   read (*,'(a)') flnm
	   write(6,'(a)') blank
	   call flush(6)
	   write (6,'(i5,a,a)') n,') input file: ',trim(flnm)
	   call flush(6)
	   write (6,'(a)') blank
	   call flush(6)

	   open (unit=9, file=flnm,form='unformatted', status='OLD',
     &access='DIRECT', recl=nrecl, action='READ')

	 if (mod(n,3).EQ.1) then
!  Initialize density.
           sig(:,:,:) = 0.0

!  read density ------------------------------------------
	   do k=1,kdm !loop over depth
	      read(9,rec=k) w
	      
	      do j=1,blkj

		 do i=1,blki

		sig(i,j,k)=w(i+(j-1)*blki)

		 enddo 	! loop over i
	      enddo 	! loop over j
	   enddo 	! loop over depth (k=1,kdm)
	 close(9)

!c	   add in the 'halo'
!	   do j=1,halo
!	     do i=1,blki
!c Take a value that is 'reflected at the northern boundary'
!c       Oh wait, that's dumb. Don't reflect those values N/S.
!	      !sig(i,blkj+j,:)=sig(i,blkj-j-1,:)
!	     enddo !i
!  	   enddo !j

!c do the x-s on both sides
!	   sig(-2:0,:,:) = sig(blki-halo+1:blki,:,:)
!	   sig(blkj+1:idm+halo,:,:) = sig(1:3,:,:)

!c pad the southern boundary with zeros
!	   sig(:,1:3,:)=0.0

!c save sig density ------------------------------------------
       do jblk=jblks,jblke
        do iblk=iblks,iblke
		outfil='sig/sig_'//runid//'_blk_'//numstr(jblk)//'_'//
     &numstr(iblk)//'.BinF'
	 
	   open (unit=99,file=outfil,status='unknown',form='unformatted',
     &access='APPEND')

	   do k=1,kdm
!	     write(99) sig(1+blki*(iblk-1)-halo:blki*iblk+halo,
!     &1+blkj*(jblk-1)-halo:blkj*jblk+halo,k)
             write(99) sig(1-halo:blki+halo,1-halo:blkj+halo,k)
	   enddo !k
	   close(99)

	   write(6,'(2a)') 'saved to: ', outfil
	   write(6,'(a)') blank
	   call flush(6)
        enddo !iblk
       enddo !jblk
	 else if (mod(n,3).EQ.2) then


! read u-velocity ------------------------------------------

	   do k=1,kdm !loop over depth
	      read(9,rec=k) w
	      
	      do j=1,blkj

		 do i=1,blki

		u(i,j,k)=w(i+(j-1)*blki)

		 enddo 	! loop over i
	      enddo 	! loop over j
	   enddo 	! loop over depth (k=1,kdm)
	 close(9)

!c	   add in the 'halo'
!	   do j=1,halo
!	     do i=1,blki

!c Take a value that is 'reflected at the northern boundary'
!	      u(i,jdm+j,:)=u(i,jdm-j-1,:)
!	     enddo !i
!  	   enddo !j

!c do the x-s on both sides
!	   u(-2:0,:,:) = u(idm-halo+1:idm,:,:)
!	   u(idm+1:idm+halo,:,:) = u(1:3,:,:)

!c pad the southern boundary with zeros
!	   u(:,1:3,:)=0.0

!c save u velocity ------------------------------------------
      do jblk=jblks,jblke
       do iblk=iblks,iblke
		outfil='u_iso/u_'//runid//'_blk_'//numstr(jblk)//'_'//
     &numstr(iblk)//'.BinF'
	 
	   open (unit=99,file=outfil,status='unknown',form='unformatted',
     &access='APPEND')
	   do k=1,kdm
!	     write(99) u(1+blki*(iblk-1)-halo:blki*iblk+halo,
!     &1+blkj*(jblk-1)-halo:blkj*jblk+halo,k)
             write(99) u(1-halo:blki+halo, 1-halo:blkj+halo, k)
	   enddo !k
	   close(99)

	   write(6,'(2a)') 'saved to: ', outfil
	   write(6,'(a)') blank
	   call flush(6)
       enddo !iblk
      enddo !jblk

	 else if (mod(n,3).EQ.0) then

! read v-velocity ------------------------------------------

	   do k=1,kdm !loop over depth
	      read(9,rec=k) w
	      
	      do j=1,blkj

		 do i=1,blki

		v(i,j,k)=w(i+(j-1)*blki)

		 enddo 	! loop over i
	      enddo 	! loop over j
	   enddo 	! loop over depth (k=1,kdm)
	 close(9)

!c	   add in the 'halo'
!	   do j=1,halo
!	     do i=1,blki
!c Take a value that is 'reflected at the northern boundary'
!	      v(i,jdm+j,:)=v(i,jdm-j-1,:)
!	     enddo !i
!  	   enddo !j

!c do the x-s on both sides
!	   v(-2:0,:,:) = v(idm-halo+1:idm,:,:)
!	   v(idm+1:idm+halo,:,:) = v(1:3,:,:)

!c pad the southern boundary with zeros
!	   v(:,1:3,:)=0.0

!c save v velocity ------------------------------------------
       do jblk=jblks,jblke
        do iblk=iblks,iblke
		outfil='v_iso/v_'//runid//'_blk_'//numstr(jblk)//'_'//
     &numstr(iblk)//'.BinF'
	 
	   open (unit=99,file=outfil,status='unknown',form='unformatted',
     &access='APPEND')
	   do k=1,kdm
!	     write(99) v(1+blki*(iblk-1)-halo:blki*iblk+halo,
!     &1+blkj*(jblk-1)-halo:blkj*jblk+halo,k)
             write(99) v(1-halo:blki+halo, 1-halo:blkj+halo, k)
	   enddo !k
	   close(99)

	   write(6,'(2a)') 'saved to: ', outfil
	   write(6,'(a)') blank
	   call flush(6)
        enddo !iblk
       enddo !jblk
	 endif ! check on which variable to read.
! 
	enddo ! n=1,numfls

	stop
	end
