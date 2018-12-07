      PROGRAM mean_rho
      use omp_lib 

      IMPLICIT NONE
      
c---USED TO TEST THE EXTRACTION OF A SINGLE FIELD: THKNSS---------      
c mean_T_S_c004_v1.f, USM, 2017-11-06
c   extract values for c004 grid
c mean_T_S_c008_v1.f, USM, 2017-9-11
c   extract T and S and compute time-mean values
c mean_rho_c008_v1.f, USM, 2017-04-27
c global extraction of z gridded density
c loads a global *.a file, adds halo, and saves chunks

      integer idm,jdm,kdm,blki,blkj,halocd,halo
      integer jblks,jblke !1,35
      integer iblks,iblke !1,90
      real flag
      real maxvar 

c CHANGE THIS!!

c glbc008
c      parameter(idm=4500,jdm=3528,kdm=240) 
c      parameter (blki=150,blkj=196,halo=3)
c glbc004
      parameter(idm=9000,jdm=7055,kdm=240) !100 x 90 = 9000 and 200 x 35 = 7000
      parameter (blki=150,blkj=200,halo=3)

      parameter(flag=2.0**100)
      character flnm*240,outfil*240,flnm2*240,flnm3*240,flnm4*240
      character flnm5*240
      character blank*5,runid*3,intid*2
      character*2 numstr(90)
      integer i,j,k,n2drec,nrecl,numfls,iblk,jblk,n,ic,flip_mcb,cc,kk

      real, allocatable, dimension(:) :: w
      real, allocatable, dimension(:) :: ww
      real, allocatable, dimension(:,:,:) :: Ta
      real, allocatable, dimension(:,:,:) :: Sa

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
c
c      initialize I/O buffer
c
      allocate( w(n2drec) )
      allocate( ww(n2drec) )
      allocate( Ta(1-halo:idm+halo,1-halo:jdm+halo,kdm) )
      allocate( Sa(1-halo:idm+halo,1-halo:jdm+halo,kdm) )

      Ta(:,:,:) = 0.0
      Sa(:,:,:) = 0.0
c
      inquire(iolength=nrecl) w
c
! read begin and end of blocks
      read(*,*) jblks
      read(*,*) jblke
      read(*,*) iblks
      read(*,*) iblke
      read(*,*) runid
      read(*,*) intid

      write (6,'(a,a3)') ' runid: ',runid
      call flush(6)

      write (6,'(a,a3)') ' intid: ',intid
      call flush(6)

      write(6,'(a)') blank 
      call flush(6)
      write (6,'(a,4i3)') ' jblks jblke iblks iblke: ',jblks,jblke,
     & iblks,iblke
      call flush(6)

      read(*,*) numfls
c      numfls = numfls - 4 !exclude grid, depth, and 2 drag files

      write(6,'(a)') blank 
      call flush(6)
      write (6,'(a,i5)') ' number of time steps: ',numfls
      call flush(6)

! skip over the first files
      read (*,'(a)') flnm2
      read (*,'(a)') flnm3
      read (*,'(a)') flnm4
      read (*,'(a)') flnm5
      write(6,'(a)') blank

      do n=1,numfls

        read (*,'(a)') flnm
	write(6,'(a)') blank 
        call flush(6)
        write (6,'(i5,a,a)') n,') input file: ',trim(flnm)
        call flush(6)
	write(6,'(a)') blank
	call flush(6)

        open(unit=9, file=flnm,form='unformatted', status='OLD',
     &     access='direct', recl=nrecl, action='READ')

        read (*,'(a)') flnm
        write(6,'(a)') blank
        call flush(6)
        write (6,'(i5,a,a)') n,') input file: ',trim(flnm)
        call flush(6)
        write(6,'(a)') blank
        call flush(6)

        open(unit=10, file=flnm,form='unformatted', status='OLD',
     &     access='direct', recl=nrecl, action='READ')


! read T and S  

        do k=1,kdm ! loop over depth

          read(9,rec=k) w
          read(10,rec=k) ww

          do j=1,jdm

            do i=1,idm

      Ta(i,j,k) = Ta(i,j,k) + w(i+(j-1)*idm)/float(numfls)
      Sa(i,j,k) = Sa(i,j,k) + ww(i+(j-1)*idm)/float(numfls)

            if (k.eq.1.and.i.eq.1450.AND.j.eq.1864) then
                write(6,*) w(i+(j-1)*idm),ww(i+(j-1)*idm)
                call flush(6) 
            end if

	    enddo

	  enddo

        enddo ! loop over depth

        close(9)
        close(10)
c        stop



      enddo ! n=1,numfls

      do j=1,halo
        do i=1,idm

c padd the northern boundary
c does not matter with what
              Ta(i,jdm+j,:) =Ta(i,jdm-j-1,:)
              Sa(i,jdm+j,:) =Sa(i,jdm-j-1,:)
        enddo !j
      enddo !i

c padd the east west boundaries
      Ta(-2:0,:,:)   =Ta(idm-halo+1:idm,:,:)
      Sa(-2:0,:,:)   =Sa(idm-halo+1:idm,:,:)
      Ta(idm+1:idm+halo,:,:)   =Ta(1:3,:,:)
      Sa(idm+1:idm+halo,:,:)   =Sa(1:3,:,:)

c padd the south boundary with zeros
      Ta(:,1:3,:)=0.0
      Sa(:,1:3,:)=0.0

c      stop


c save data =======================================================
      do jblk=jblks,jblke !1,35

       do iblk=iblks,iblke !1,90

c save time-mean S and T

          outfil='zTavg/Tavg_'//intid//'_'//runid//
     &'_blk_'//numstr(jblk)//'_'
     &//numstr(iblk)//'.BinF'

         open(unit=99,file=outfil,status='replace',form='unformatted',
     &             access='append')


          outfil='zSavg/Savg_'//intid//'_'//runid//
     &'_blk_'//numstr(jblk)//'_'
     &//numstr(iblk)//'.BinF'

         open(unit=100,file=outfil,status='replace',form='unformatted',
     &             access='append')

         do k=1,kdm 
 
           write(99) Ta(1+blki*(iblk-1)-halo:blki*iblk+halo,
     &              1+blkj*(jblk-1)-halo:blkj*jblk+halo,k)

           write(100) Sa(1+blki*(iblk-1)-halo:blki*iblk+halo,
     &              1+blkj*(jblk-1)-halo:blkj*jblk+halo,k)

         enddo
         close(99) 
         close(100) 

         write(6,'(2a)') 'saved to: ',  outfil
         write(6,'(a)') blank 
         call flush(6)

        enddo

      enddo
 

      do k=1,kdm
        write(6,*) k,Ta(1450,1864,k),Sa(1450,1864,k)
        call flush(6)
      enddo

      write(6,*) 'calculating mean T and S is finished'
      stop
      end
