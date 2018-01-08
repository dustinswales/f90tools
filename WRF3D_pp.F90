! #############################################################################
!
! The purpose of this program is to post-process WRF 3D data. Non-standard WRF
! outputs (e.g IVT and freezing-level height) are computed and output to a
! netCDF file.
!
! Dustin Swales (4/2017) - Initial version (dustin.swales@noaa.gov)
! Dustin Swales (1/2018) - Added computation of freezing-level height.
!
! To compile:
! ifort -I/usr/local/ifort/include -L/usr/local/ifort/lib -lnetcdff WRF3D_pp.F90 -o WRF3D_pp.o
!
! To run:
! ./WRF3D_pp fileIN fileOUT
!
! Inputs
! fileIN     - WRF input file
!
! Outputs
! fileOUT    - Output file
! #############################################################################
program WRF3D_pp
  use netcdf
  implicit none

  ! Parameters
  logical, parameter :: &
       verbose=.true.  ! On/off
  real,parameter :: &
       z_soil = 0.1     ! Depth of WRF soil moisture to output (m)

  ! Arguments
  character(len=256) :: &
       fileIN, & ! Input file
       fileOUT   ! Output file

  ! WRF fields
  integer :: &
       nTime,     & ! WRF file dimension: Number of times
       nLon,      & ! WRF file dimension: Number of longitudes
       nLat,      & ! WRF file dimension: Number of latitudes
       nLev,      & ! WRF file dimension: Number of vertical levels
       nLon_stag, & ! WRF file dimension: Number of longitudes (staggerd grid)
       nLat_stag, & ! WRF file dimension: Number of latitude (staggerd grid)
       nLev_stag, & ! WRF file dimension: Number of vertical levels (staggerd grid)
       nSoil_stag   ! WRF file dimension: Number of soil layers (staggered grid) 
  real,dimension(:,:,:,:),allocatable :: &!                              WRF NAME   UNITS             
       q,         & ! WRF input field: Water vapor mixing ratio          (QVAPOR)  (kg/kg)
       u,         & ! WRF input field: Zonal component of the wind       (U)       (m/s)
       v,         & ! WRF input field: Meridional component of the wind  (V)       (m/s)
       pp,        & ! WRF input field: Pertubation pressure              (P)       (Pa)
       pb,        & ! WRF input field: Base state pressure               (PB)      (Pa)
       ta,        & ! WRF input field: Pertubation potential temp.       (T)       (K)
       ph,        & ! WRF input field: Pertubation geopotential          (PH)      (m2/s2)
       phb,       & ! WRF input field: Base-state geopotential           (PHB)     (m2/s2)
       smois        ! WRF input field: Soil mosisture                    (SMOIS)   (m3/m3)
  real,dimension(:,:,:),allocatable :: &
       lon,       & ! WRF input field: Longitude                         (XLONG)   (degree_east)
       lat,       & ! WRF input field: Latitude                          (XLAT)    (degree_north)
       lon_u,     & ! WRF input field: Longitude (staggerd u-grid)       (XLONG_U) (degree_east)
       lat_u,     & ! WRF input field: Latitude (staggerd u-grid)        (XLAT_U)  (degree_north)
       lon_v,     & ! WRF input field: Longitude (staggerd v-grid)       (XLONG_V) (degree_east)
       lat_v,     & ! WRF input field: Latitude (staggerd v-grid)        (XLAT_V)  (degree_north)
       psfc,      & ! WRF input field: Surface pressure                  (PSFC)    (Pa)
       terrainZ     ! WRF input field: Terrain height                    (HGT)     (m)
  real,dimension(:),allocatable :: &
       t00,       & ! WRF input field: Base state temperature            (T00)     (K)
       zs           ! WRF input field: Soil levels                       (ZS)      (m)

  ! Local fields
  real,dimension(:,:,:,:),allocatable :: &
       p,         & ! WRF pressue (computed from pp and pb)  (pa)
       hgt          ! WRF heights (computed from ph and phb) (m)
  real,dimension(:,:,:),allocatable :: &
       ivtU,      & ! WRF IVT (u-component)                  (kg/m/s)
       ivtV,      & ! WRF IVT (v-component)                  (kg/m/s)
       z0k,       & ! WRF Freezing level height              (m)
       soilMoisture ! WRF soil mositure @ z_soil              (m3/m3)
  
  integer,dimension(8) :: dimID,varIDout
  integer :: status,fileID,varID,ii,ij,ik,il
  real,dimension(:),allocatable :: ui,vi,temp,phm,a,b
  real :: wt1,dp,wt2
  integer,dimension(1) :: xi,xf

  ! Get filenames.
  call getarg(1,fileIN)
  call getarg(2,fileOUT)
  
  ! #############################################################################
  ! Part A: Data ingest
  ! #############################################################################
  if (verbose) print*,'Reading in data ....'
  
  ! 0) Open file
  status = nf90_open(trim(fileIN),nf90_NoWrite,fileID)
  if (status /= nf90_NoErr) then
     print*,'ERROR: Cannot open input file ',trim(fileIN)
     print*,'EXITING!!!'
     goto 101
  end if

  ! 1) Get dimensions
  ! "Time" 
  status = nf90_inq_dimid(fileID,"Time",dimID(1))
  if (status /= nf90_NoErr) print*,'ERROR: Dimension not present in file, Time'
  if (status == nf90_NoErr) then
     status = nf90_inquire_dimension(fileID,dimID(1),len=nTime)
  endif
  ! "west_east"
  status = nf90_inq_dimid(fileID,"west_east",dimID(2))
  if (status /= nf90_NoErr) print*,'ERROR: Dimension not present in file, west_east'
  if (status == nf90_NoErr) then
     status = nf90_inquire_dimension(fileID,dimID(2),len=nLon)
  endif
  ! "south_north"
  status = nf90_inq_dimid(fileID,"south_north",dimID(3))
  if (status /= nf90_NoErr) print*,'ERROR: Dimension not present in file, south_north'
  if (status == nf90_NoErr) then
     status = nf90_inquire_dimension(fileID,dimID(3),len=nLat)
  endif
  ! "bottom_top"
  status = nf90_inq_dimid(fileID,"bottom_top",dimID(4))
  if (status /= nf90_NoErr) print*,'ERROR: Dimension not present in file, bottom_top'
  if (status == nf90_NoErr) then
     status = nf90_inquire_dimension(fileID,dimID(4),len=nLev)
  endif
  ! "west_east_stag"
  status = nf90_inq_dimid(fileID,"west_east_stag",dimID(5))
  if (status /= nf90_NoErr) print*,'ERROR: Dimension not present in file, west_east_stag'
  if (status == nf90_NoErr) then
     status = nf90_inquire_dimension(fileID,dimID(5),len=nLon_stag)
  endif
  ! "south_north_stag"
  status = nf90_inq_dimid(fileID,"south_north_stag",dimID(6))
  if (status /= nf90_NoErr) print*,'ERROR: Dimension not present in file, south_north_stag'
  if (status == nf90_NoErr) then
     status = nf90_inquire_dimension(fileID,dimID(6),len=nLat_stag)
  endif
  ! "bottom_top_stag"
  status = nf90_inq_dimid(fileID,"bottom_top_stag",dimID(7))
  if (status /= nf90_NoErr) print*,'ERROR: Dimension not present in file, bottom_top_stag'
  if (status == nf90_NoErr) then
     status = nf90_inquire_dimension(fileID,dimID(7),len=nLev_stag)
  endif
  ! "soil_layers_stag"
  status = nf90_inq_dimid(fileID,"soil_layers_stag",dimID(8))
  if (status /= nf90_NoErr) print*,'ERROR: Dimension not present in file, soil_layers_stag'
  if (status == nf90_NoErr) then
     status = nf90_inquire_dimension(fileID,dimID(8),len=nSoil_stag)
  endif
    
  ! 2) Allocate space for input fields
  allocate(q(nLon,nLat,nLev,nTime), u(nLon_stag,nLat,nLev,nTime), v(nLon,nLat_stag,nLev,nTime), &
          lat(nLon,nLat,nTime), lon(nLon,nLat,nTime), lon_u(nLon_stag,nLat,nTime),              & 
          lat_u(nLon_stag,nLat,nTime), lon_v(nLon,nLat_stag,nTime), lat_v(nLon,nLat_stag,nTime),&
          pp(nLon,nLat,nLev,nTime), pb(nLon,nLat,nLev,nTime) ,p(nLon,nLat,nLev,nTime),          &
          psfc(nLon,nLat,nTime), ta(nLon,nLat,nLev,nTime), ph(nLon,nLat,nLev_stag,nTime),       &
          phb(nLon,nLat,nLev_stag,nTime), hgt(nLon,nLat,nLev_stag,nTime), t00(nTime),           &
          terrainZ(nLon,nLat,nTime),zs(nSoil_stag),smois(nLon,nLat,nSoil_stag,nTime))

  ! 3) Read in fields
  ! Longitude
  status = nf90_inq_varid(fileID,"XLONG",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, XLONG'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,lon)
  endif
  
  ! Latitude
  status = nf90_inq_varid(fileID,"XLAT",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, XLAT'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,lat)
  endif
  
  ! Longitude (staggered in zonal)
  status = nf90_inq_varid(fileID,"XLONG_U",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, XLONG_U'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,lon_u)
  endif
  
  ! Latitude (staggered in zonal)
  status = nf90_inq_varid(fileID,"XLAT_U",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, XLAT_U'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,lat_u)
  endif
  
  ! Longitude (staggered in meridional)
  status = nf90_inq_varid(fileID,"XLONG_V",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, XLONG_V'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,lon_v)
  endif
  
  ! Latitude (staggered in meridional)
  status = nf90_inq_varid(fileID,"XLAT_V",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, XLAT_V'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,lat_v)
  endif

  ! Soil depth (staggered in vertical)
  status = nf90_inq_varid(fileID,"ZS",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, ZS'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,zs,count=(/nSoil_stag,1/))
  endif
  
  ! Terrain height
  status = nf90_inq_varid(fileID,"HGT",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, HGT'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,terrainZ)
  endif
  
  ! Base temperature
  status = nf90_inq_varid(fileID,"T00",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, T00'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,t00)
  endif
  
  ! Temperature
  status = nf90_inq_varid(fileID,"T",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, T'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,ta,count=(/nLon,nLat,nLev,nTime/))
  endif

  ! Specific humidity
  status = nf90_inq_varid(fileID,"QVAPOR",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, QVAPOR'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,q,count=(/nLon,nLat,nLev,nTime/))
  endif
  
  ! Eastward component of wind
  status = nf90_inq_varid(fileID,"U",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, U'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,u,count=(/nLon_stag,nLat,nLev,nTime/))
  endif
  
  ! Northward component of wind
  status = nf90_inq_varid(fileID,"V",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, V'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,v,count=(/nLon,nLat_stag,nLev,nTime/))
  endif

  ! Base state (reference) height
  status = nf90_inq_varid(fileID,"PHB",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, PHB'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,phb,count=(/nLon,nLat,nLev_stag,nTime/))
  endif
  
  ! Pertubation height
  status = nf90_inq_varid(fileID,"PH",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, PH'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,ph,count=(/nLon,nLat,nLev_stag,nTime/))
  endif
  
  ! Base state (reference) pressure
  status = nf90_inq_varid(fileID,"PB",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, PB'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,pb,count=(/nLon,nLat,nLev,nTime/))
  endif
  
  ! Pertubation pressure
  status = nf90_inq_varid(fileID,"P",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, P'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,pp,count=(/nLon,nLat,nLev,nTime/))
  endif
  
  ! Surface pressure
  status = nf90_inq_varid(fileID,"PSFC",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, PSFC'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,psfc,count=(/nLon,nLat,nTime/))
  endif

  ! Soil moisture
  status = nf90_inq_varid(fileID,"SMOIS",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, SMOIS'
  if (status == nf90_NoErr) then
     status = nf90_get_var(fileID,varID,smois)
  endif
  
  ! 4) Close file
  status=nf90_close(fileID)
  if (status /= nf90_NoErr) print*,'ERROR: Could not close file, ',trim(fileiN)
  if (verbose) print*,'Finished reading in data'

  ! #############################################################################
  ! Part B: Computations
  ! Interpolate staggered velocity grid to mass-centered grid points, compute
  ! the integrated vapor transport (IVT), and compute freezing-level height (Z0K)
  ! #############################################################################

  ! Compute pressure and height.
  p   = pb+pp
  hgt = (ph+phb)/9.8
  deallocate(pb,pp,ph,phb)

  ! Allocate
  allocate(ui(nLev),vi(nLev), ivtU(nLon,nLat,nTime), ivtV(nLon,nLat,nTime), &
       z0k(nLon,nLat,nTime),temp(nLev),phm(nLev),a(nLev),b(nSoil_stag),soilMoisture(nLon,nLat,nTime))

  ! Initialize
  ivtU(:,:,:) = 0.
  ivtV(:,:,:) = 0.
  z0k(:,:,:)  = 0.

  ! Loop over all points/times and compute IVT and freezing level heights.
  if (verbose) print*,'Regridding velocity fields, computing IVT and freezing-level height: '
  do ii=1,nTime
     if (verbose) write(*,"(a12,i2,a4,i2)"),'   @Timestep ',ii,' of ',nTime
     do ij=1,nLon
        do ik=1,nLat
           ! ###################################################################
           ! First, put winds on mass centered grid points, in WRF the velocity
           ! components are on staggered grids.
           ! ###################################################################
           ! Eastward component of wind (on mass-centered point)
           wt1 = (lon_u(ij+1,ik,ii)-lon(ij,ik,ii))/(lon_u(ij+1,ik,ii)-lon_u(ij,ik,ii))
           ui = wt1*u(ij,ik,:,ii)+(1-wt1)*u(ij+1,ik,:,ii)
           ! Northward component of wind (on mass-centered point)
           wt1 = (lat_v(ij,ik+1,ii)-lat(ij,ik,ii)) / (lat_v(ij,ik+1,ii)-lat_v(ij,ik,ii))
           vi = wt1*v(ij,ik,:,ii)+(1-wt1)*v(ij,ik+1,:,ii)

           ! ###################################################################
           ! Compute integrated vapor transport (IVT)
           ! ###################################################################
           do il=1,nLev
              ! First, compute pressure change across layer.
              if (il .eq. 1)                    dp = psfc(ij,ij,ii)-sum(p(ij,ik,il:il+1,ii))/2.              ! Bottommost level
              if (il .ne. 1 .and. il .ne. nLev) dp = sum(p(ij,ik,il-1:il,ii))/2.-sum(p(ij,ik,il:il+1,ii))/2. ! Middle levels
              if (il .eq. nLev)                 dp = 0.                                                      ! Top level
              ! Integrate vertically.
              ivtU(ij,ik,ii) = ivtU(ij,ik,ii) + ui(il)*q(ij,ik,il,ii)*dp/9.8
              ivtV(ij,ik,ii) = ivtV(ij,ik,ii) + vi(il)*q(ij,ik,il,ii)*dp/9.8
           enddo

           ! ###################################################################
           ! Compute freezing-level height.
           ! ###################################################################
           ! Compute temperature from pertubation potential temperature.
           temp = (ta(ij,ik,:,ii)+t00(ii))*(101325.0/p(ij,ik,:,ii))**(-0.286)

           ! Geopotential heights are on a staggerd vertical grid, so compute mass-
           ! centered geopotential height.
           phm  = (hgt(ij,ik,1:nlev_stag-1,ii)+hgt(ij,ik,2:nlev_stag,ii))*0.5

           ! Only compute if freezing level is present.
           if (any(temp .gt. 273.16)) then
              ! What is the highest level that is above freezing?
              a  = temp-273.15
              xi = minloc(a,temp .gt. 273.16)
              xf = xi+1
              ! Interpolate to find freezing level height.
              wt1 = (273.15-temp(xf(1)))/(temp(xi(1))-temp(xf(1)))
              z0k(ij,ik,ii) = phm(xi(1))*wt1 + phm(xf(1))*(1-wt1)
           else
              ! If it's freezing at all layers, set freezing level to terrain height.
              z0k(ij,ik,ii) = terrainZ(ij,ik,ii)
           endif

           ! ###################################################################
           ! Compute soil moisture at z_soil.
           ! ###################################################################
           b  = zs-z_soil
           xf = minloc(b,b .gt. 0)
           xi = xf-1
           if (xi(1) .lt. 1) then
              write(*,"(a37,f5.2,a17)"),'ERROR: Requested soil moisture depth, ',z_soil,' is out of bounds.'
           endif
           wt2 = (zs(xf(1))-z_soil)/(zs(xf(1))-zs(xi(1)))
           soilMoisture(ij,ik,ii) = smois(ij,ik,xi(1),ii)*wt2+smois(ij,ik,xf(1),ii)*(1-wt2)

        enddo ! Latitude
     enddo    ! Longitude
  enddo       ! Time

  ! Cleanup a tad.
  deallocate(u,v,ui,vi,a,b)
  if (verbose) print*,'Finished computations'

  ! #############################################################################
  ! Part C: Output
  ! #############################################################################
  ! Create output file
  status = nf90_create(trim(fileOUT),cmode=0,ncid=fileID)
  if (status /= nf90_NoErr) print*,'ERROR: Failure making output file'
  if (verbose) print*,'Writing to output file: ',trim(fileOUT)

  ! Define dimensions
  status = nf90_def_dim(fileID,"Time",nTime,dimID(1))
  if (status /= nf90_NoErr) print*,'ERROR: Failure making dimension ID, Time'
  status = nf90_def_dim(fileID,"west_east",nLon,dimID(2))
  if (status /= nf90_NoErr) print*,'ERROR: Failure making dimension ID, west_east'
  status = nf90_def_dim(fileID,"south_north",nLat,dimID(3))
  if (status /= nf90_NoErr) print*,'ERROR: Failure making dimension ID, south_north'
  status = nf90_def_dim(fileID,"bottom_top",nLev,dimID(4))
  if (status /= nf90_NoErr) print*,'ERROR: Failure making dimension ID, bottom_top'

  ! Define variables
  status = nf90_def_var(fileID,'XLAT',nf90_float, (/dimID(2),dimID(3)/),varIDout(1))
  if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, XLAT'
  status = nf90_def_var(fileID,'XLONG',nf90_float,(/dimID(2),dimID(3)/),varIDout(2))
  if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, XLONG'
  status = nf90_def_var(fileID,'IVTU',nf90_float, (/dimID(2),dimID(3),dimID(1)/),varIDout(3))
  if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, IVTU'
  status = nf90_def_var(fileID,'IVTV',nf90_float, (/dimID(2),dimID(3),dimID(1)/),varIDout(4))
  if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, IVTV'
  status = nf90_def_var(fileID,'Z0K',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(5))
  if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Z0K'
  status = nf90_def_var(fileID,'SMOIS',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(6))
  if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, SMOIS'  
  
  ! Add attributes to variables
  status = nf90_put_att(fileID,varIDout(1),"FieldType",104)
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable XLAT'
  status = nf90_put_att(fileID,varIDout(1),"MemoryOrder","XY")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable XLAT'
  status = nf90_put_att(fileID,varIDout(1),"description","LATITUDE, SOUTH IS NEGATIVE")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable XLAT'
  status = nf90_put_att(fileID,varIDout(1),"units","degree_north")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable XLAT'
  status = nf90_put_att(fileID,varIDout(1),"stagger","")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable XLAT'  
  status = nf90_put_att(fileID,varIDout(2),"FieldType",104)
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable XLONG'
  status = nf90_put_att(fileID,varIDout(2),"MemoryOrder","XY")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable XLONG'
  status = nf90_put_att(fileID,varIDout(2),"description","LONGITUDE, WEST IS NEGATIVE")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable XLONG'
  status = nf90_put_att(fileID,varIDout(2),"units","degree_east")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable XLONG'
  status = nf90_put_att(fileID,varIDout(2),"stagger","")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable XLONG'
  status = nf90_put_att(fileID,varIDout(3),"FieldType",104)
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable IVTU'
  status = nf90_put_att(fileID,varIDout(3),"MemoryOrder","XY")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable IVTU'
  status = nf90_put_att(fileID,varIDout(3),"description","ZONAL COMPONENT OF THE INTEGRATED VAPOR TRANSPORT")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable IVTU'
  status = nf90_put_att(fileID,varIDout(3),"units","kg m-2 s-1")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable IVTU'
  status = nf90_put_att(fileID,varIDout(3),"stagger","")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable IVTU'
  status = nf90_put_att(fileID,varIDout(4),"FieldType",104)
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable IVTV'
  status = nf90_put_att(fileID,varIDout(4),"MemoryOrder","XY")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable IVTV'
  status = nf90_put_att(fileID,varIDout(4),"description","MERIDIONAL COMPONENT OF THE INTEGRATED VAPOR TRANSPORT")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable IVTV'
  status = nf90_put_att(fileID,varIDout(4),"units","kg m-2 s-1")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable IVTV'
  status = nf90_put_att(fileID,varIDout(4),"stagger","")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable IVTV'
  status = nf90_put_att(fileID,varIDout(5),"FieldType",104)
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Z0K'
  status = nf90_put_att(fileID,varIDout(5),"MemoryOrder","XY")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Z0K'
  status = nf90_put_att(fileID,varIDout(5),"description","FREEZING LEVEL HEIGHT")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Z0K'
  status = nf90_put_att(fileID,varIDout(5),"units","m")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Z0K'
  status = nf90_put_att(fileID,varIDout(5),"stagger","")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Z0K'

  status = nf90_put_att(fileID,varIDout(6),"FieldType",104)
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable SMOIS'
  status = nf90_put_att(fileID,varIDout(6),"MemoryOrder","XY")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable SMOIS'
  status = nf90_put_att(fileID,varIDout(6),"description","SOIL MOISTURE @ Z_soil")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable SMOIS'
  status = nf90_put_att(fileID,varIDout(6),"units","m")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable SMOIS'
  status = nf90_put_att(fileID,varIDout(6),"stagger","")
  if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable SMOIS'

  ! Exit define mode
  status = nf90_enddef(fileID)
  if (status /= nf90_NoErr) print*,'ERROR: Failure exiting define mode'

  ! Populate file
  status = nf90_put_var(fileID,varIDout(1),lat(:,:,1))
  if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, XLAT'
  status = nf90_put_var(fileID,varIDout(2),lon(:,:,1))
  if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, XLONG'
  status = nf90_put_var(fileID,varIDout(3),ivtU)
  if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, IVTU'
  status = nf90_put_var(fileID,varIDout(4),ivtV)
  if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, IVTV'
  status = nf90_put_var(fileID,varIDout(5),z0k)
  if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z0K'
  status = nf90_put_var(fileID,varIDout(6),soilMoisture)
  if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, SMOIS'
  
  ! Close output file
  status = nf90_close(fileID)
  if (status /= nf90_NoErr) print*,'ERROR: Failure closing file'

  ! Cleanup
  deallocate(q,lon,lat,lon_u,lat_u,lon_v,lat_v,p,psfc,ivtV,ivtU)
  
101 continue
  ! #############################################################################
  ! END PROGRAM
  ! #############################################################################
end program WRF3D_pp
