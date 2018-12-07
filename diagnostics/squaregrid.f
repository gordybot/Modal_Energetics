      PROGRAM SQUAREPGRID
      IMPLICIT NONE

      ! A program to move data from the tripolar p-grid to a more regular square grid.
      ! Let's assume we're dealing with a 1/25 degree resolution field.

      ! Input variable: plon, plat, pVar.
      ! Output variables: glon, glat, gVar.

      ! Declare variables. 
      ! INPUTS
      real, allocatable, dimension(:,:) :: plon
      real, allocatable, dimension(:,:) :: plat
      real, allocatable, dimension(:,:) :: pVar
      ! OUTPUTS
      real, allocatable, dimension(:,:) :: glon    ! 
      real, allocatable, dimension(:,:) :: glat
      real, allocatable, dimension(:,:) :: gVar    ! Gridded variable

      ! Helper variables:
      character infileLon*240, infileLat*240, infileVar*240
      character outfilLon*240, outfilLat*240, outfileVar*240
      real minLon, maxLon, minLat, maxLat
      real p_resln, g_resln
      integer idm, jdm, gnx, gny, i, j, pind, pjnd, di, dj, mbdy
      real lon, lat, lonp, latp
      logical inbounds

      ! Load plon, plat. 
      ! Longitude ranges from 74 to 434. 
      ! Gridded latitude ranges from -80 to 80, I guess.
      mbdy = 3
      minLon = 74.14
      maxLon = 434.14
      minLat = -80.0
      maxLat = 80.0
      p_resln = 1.0/25
      g_resln = 1.0/25

      idm = 9000   ! size of pgrid x-dim / lon
      jdm = 7055   ! size of pgrid y-dim / lat

      ! TODO cast these to integers.
      gnx = NINT( (maxLon - minLon) / g_resln )
      gny = NINT( (maxLat - minLat) / g_resln )
      write(6,*) 'Lon: ',minLon,' ',maxLon
      call flush(6)

      write(6,*) 'Lat: ',minLat,' ',maxLat
      call flush(6)

      write(6,*) 'Grid: ',gnx,' by ',gny
      call flush(6)
       
      write(6,*) 'p_res: ',p_resln,' deg.| g_res: ',g_resln,' deg.'
      call flush(6)

      ! Size of gridded array.
      allocate( plon(1-mbdy:idm+mbdy, 1:jdm) )
      allocate( plat(1-mbdy:idm+mbdy, 1:jdm) )
      allocate( pVar(1-mbdy:idm+mbdy, 1:jdm) )
      allocate( glon(1:gnx, 1:gny) )
      allocate( glat(1:gnx, 1:gny) )
      allocate( gVar(1:gnx, 1:gny) )

      infileVar = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &      //'modes_global/glb_KE_mode01_feb8.BinF'

      ! Read in the global pgrid lon, lat, Var.
       infileLon = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &  //'modes_global/glb_lon_221.BinF'
!     &   //'_'//runid//'_blk_'//numstr(jblk)//'_'//numstr(iblk)
!     &   //'.BinF'
       open(unit=10,file=infileLon,status='old',form='unformatted')
       write(6,*) 'Reading ',infileLon
       call flush(6)
       plon(:,:) = 0.0
       read(10) plon(1:idm,1:jdm)
       close(10)
       plon(1-mbdy:0,:) = plon(idm+1-mbdy:idm,:)
       plon(idm+1:idm+mbdy,:) = plon(1:mbdy,:)

       infileLat = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &    //'modes_global/glb_lat_221.BinF'
!     &   //'_'//runid//'_blk_'//numstr(jblk)//'_'//numstr(iblk)
!     &   //'.BinF'
       open(unit=11,file=infileLat,status='old',form='unformatted')
       write(6,*) 'Reading ',infileLat
       call flush(6)
       plat(:,:) = 0.0
       read(11) plat(1:idm,1:jdm)
       close(11)
       plat(1-mbdy:0,:) = plat(idm+1-mbdy:idm,:)
       plat(idm+1:idm+mbdy,:) = plat(1:mbdy,:)

       open(unit=12,file=infileVar,status='old',form='unformatted')
       write(6,*) 'Reading ',infileVar
       call flush(6)
       pVar(:,:) = 0.0
       read(12) pVar(1:idm,1:jdm)
       close(12) 
       pVar(1-mbdy:0,:) = pVar(idm+1-mbdy:idm,:)
       pVar(idm+1:idm+mbdy,:) = pVar(1:mbdy,:)

      glon(:,:) = 0.0
      glat(:,:) = 0.0
      gVar(:,:) = 0.0

      ! Fill in values of regular gridded glon, glat.
      do i = 1,gnx
        write(6,*) (1.0*i/gnx)
        call flush(6)
        do j = 1,gny
          lon = minLon + i*g_resln
          lat = minLat + j*g_resln

          ! Fill in values of regular gridded glon, glat.
          glon(i,j) = lon
          glat(i,j) = lat
          
          ! For each point in ggrid, do the hop and wiggle to find 
          ! the closest point on the pgrid.
          ! Start close to where the point would be if pgrid were rectangular.
          pind = NINT( 1.0 * i * idm/gnx )
          pjnd = NINT( 1.0 * j * jdm/gny )

          lonp = plon(pind, pjnd)
          latp = plat(pind, pjnd)

          ! Which way to hop? Hop using the nominal resolution of the pgrid to estimate
          ! number of indices to jump.
          di = - NINT( (lon - lonp / p_resln ) )
          dj = - NINT( (lat - latp / p_resln ) )
          ! Example: lon = 101.2. -> i = 272. pind = NINT( 272*2.5) = 680. Suppose plon(pind,pjnd) = 105.4
          ! Then need to shift west on the pgrid to find a closer point. As a first estimate,
          ! we need to shift westward by 42 boxes -- but this will be refined in the next step.

          pind = pind + di
          pjnd = pjnd + dj
          lonp = plon( pind, pjnd)
          latp = plat( pind, pjnd)
          inbounds = .TRUE.
          ! Next wiggle -- take steps one by one until the points are close.
          do while ( ( abs( lon - lonp).GT.p_resln .OR.
     &                abs( lat - latp).GT.p_resln) .AND. inbounds )
            pind = pind + SIGN(1.0, lon - lonp )
            pjnd = pjnd + SIGN(1.0, lat - latp )
            if ( (pind.LT.1) .OR. (pjnd.LT.1) .OR. 
     &            (pind.GT.idm) .OR. (pjnd.GT.jdm ) ) then
              inbounds = .FALSE.
            else
              lonp = plon( pind, pjnd )
              latp = plat( pind, pjnd )
            endif
          enddo ! while loop.
          ! Use pind, pjnd to assign value to gVar. 
          ! Could change this to interpolate, but I'll just try assigning a "close enough" point
          ! at this stage.
          if (inbounds) then
            gVar(i,j) = pVar( pind, pjnd ) 
          else
            gVar(i,j) = -9999.9
          endif
        enddo !j
      enddo !i

      ! Write out the regridded file.
       outfileVar = '/p/work1/grstephe/hycom/GLBc0.04/expt_22.1/'
     &         //'regridded/glb_KE_221_sqrgrid.BinF' 
       write(6,*) 'Writing regridded variable to ',outfileVar
       call flush(6)
       open(unit=99,file=outfileVar,status='replace',
     &       form='unformatted',access='sequential')
       write(99) gVar(:,:)
       close(99)

      STOP
      END PROGRAM SQUAREPGRID 
