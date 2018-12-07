       PROGRAM EXTRACT_U

       IMPLICIT NONE
       ! Redo the extraction routine for experiment 22.1.
       ! At Stennis Space Center, 12/5/2017. 11:11 AM.

       ! One issue is, can't load the whole archv file into
       ! into memory (57 Gb). 
       ! Read one time step = one archv file
       !   read one depth =    k=1,240
       !      for each tile, append relevant indices.
       ! Tiled data will look like:
       !     1:blki, 1:blkj, t=1, k=1
       !     1:blki, 1:blkj, t=1, k=2
       !       ..     ...
       !     1:blki, 1:blkj, t=1, k=240
       !     1:blki, 1:blkj, t=2, k=1
       !               ...        ...
       !     1:blki, 1:blkj, t=2, k=240

       ! Read global, single-depth into w.

       ! Declare and type variables.
       integer idm,jdm,kdm,blki,blkj,halo
       integer i,j,k,n, iblk, jblk, numfls
       integer iblks, iblke, jblks, jblke
       integer offset_i, offset_j
       integer nrecl, n2drec
       real flag

       character*2 numstr(90)
       character blank*5, runid*3, intid*3
       character flnm*240, outfil*240

       real, allocatable, dimension(:) :: w
       real, allocatable, dimension(:,:) :: glb_u
       real, allocatable, dimension(:,:) :: blk_u

       ! Assign size parameters and define some useful variables.
       parameter(idm=9000, jdm=7055, kdm=240, halo=3)
       parameter(blki=150, blkj=200)
       parameter(flag=2.0**100)
       
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

       n2drec = ((idm*jdm+4095)/4096)*4096

       allocate( w(1:n2drec) )
       allocate( glb_u(1-halo:idm+halo, 1-halo:jdm+halo) )
       allocate( blk_u(1-halo:blki+halo, 1-halo:blkj+halo) )

       inquire(iolength=nrecl) w

       ! Read beginning and end of blocks.
       read(*,*) jblks
       read(*,*) jblke
       read(*,*) iblks
       read(*,*) iblke
       read(*,*) runid
       read(*,*) intid

!c     print out some of the info
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
       write(6,'(a)') blank
       call flush(6)

       ! Read the 
       do n=1,numfls
         read(*,'(a)') flnm
         write(6,*) blank
         call flush(6)
         write (6,'(i5,a,a)') n,') input file: ',trim(flnm)
         call flush(6)
         write (6,'(a)') blank
         call flush(6)

         open (unit=9, file=flnm, form='unformatted',status='OLD',
     &     access='DIRECT', recl=nrecl, action='READ')

         glb_u(:,:) = 0.0
         
         ! Read density, one depth at a time.
         do k=1,kdm !loop over depth.
           read(9,rec=k) w
           do j=1,jdm
             do i=1,idm
               glb_u(i,j) = w( i + (j-1) * idm )
             enddo !i
           enddo !j

           ! Wrap-around to give glb_u a buffer/boundary/halo.
           do j=1,jdm
             do i=1-halo,0
               glb_u(i,j) = glb_u( idm+i ,j)
               glb_u(idm + i, j) = glb_u(i,j)
             enddo !i in halo
           enddo !j

           ! Global field at depth k has been assigned.
           ! Tile data.
           do jblk=jblks,jblke
             do iblk=iblks,iblke
               outfil='u_iso/u_iso_'//runid//'_blk_'//numstr(jblk)//
     &                '_'//numstr(iblk)//'.BinF'
               open (unit=99, file=outfil, status='unknown',
     &              form='unformatted',access='APPEND')
               
               offset_i = blki*(iblk-1)
               offset_j = blkj*(jblk-1)
               blk_u(1-halo:blki+halo,1-halo:blkj+halo) = 
     &               glb_u( (offset_i+1-halo):(offset_i+blki+halo),
     &                      (offset_j+1-halo):(offset_j+blkj+halo) )
               write(99) blk_u(1-halo:blki+halo, 1-halo:blkj+halo)
               close(99)
             enddo !iblk
           enddo !jblk
           write(6,*) 'Finished depth level ',k
           call flush(6)
         enddo !k /  each depth level.
         write(6,*) 'Finished time-step ',n
         call flush(6)
       enddo !numfils / each time step. 

       write(6,*) 'Finished...'
       call flush(6)

       stop
       END PROGRAM EXTRACT_U
