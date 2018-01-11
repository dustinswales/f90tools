! #############################################################################
!
! The purpose of this program is to post-process WRF 3D data. Non-standard WRF
! outputs (e.g IVT and freezing-level height) are computed and output to a
! netCDF file. Inputs/Outputs controlled by input_namelist, wrf3d_pp_nl.txt.
!
! Dustin Swales (4/2017) - Initial version (dustin.swales@noaa.gov)
!
! To compile:
! ifort -I/usr/local/ifort/include -L/usr/local/ifort/lib -lnetcdff WRF3D_pp.F90 -o WRF3D_pp.o
!
! #############################################################################
program WRF3D_pp
  use netcdf
  implicit none

  ! Parameters
  character(len=128),parameter :: &
       input_namelist = 'wrf3d_pp_nl.txt'
  
  ! Namelist
  character(len=256) :: &
       fileIN, & ! Input file
       fileOUT   ! Output file
  logical :: &
       verbose, & ! Turn on/off info print to screen
       l_ivt,   & ! Compute IVT?
       l_z0k,   & ! Compute freezing level height?
       l_smois, & ! Compute soil moisture?
       l_z1000, & ! Output 1000hPa geopotential?
       l_z950,  & ! Output 950hPa geopotential?
       l_z900,  & ! Output 900hPa geopotential?
       l_z850,  & ! Output 850hPa geopotential?
       l_z800,  & ! Output 800hPa geopotential?
       l_z750,  & ! Output 750hPa geopotential?
       l_z700,  & ! Output 700hPa geopotential?
       l_z600,  & ! Output 600hPa geopotential?
       l_z500,  & ! Output 500hPa geopotential?
       l_z250,  & ! Output 250hPa geopotential?
       l_q1000, & ! Output 1000hPa specific-humidity?
       l_q950,  & ! Output 950hPa specific-humidity?
       l_q900,  & ! Output 900hPa specific-humidity?
       l_q850,  & ! Output 850hPa specific-humidity?
       l_q800,  & ! Output 800hPa specific-humidity?
       l_q750,  & ! Output 750hPa specific-humidity?
       l_q700,  & ! Output 700hPa specific-humidity?
       l_q600,  & ! Output 600hPa specific-humidity?
       l_q500,  & ! Output 500hPa specific-humidity?
       l_q250,  & ! Output 250hPa specific-humidity?
       l_u1000, & ! Output 1000hPa U-component of wind??
       l_u950,  & ! Output 950hPa U-component of wind??
       l_u900,  & ! Output 900hPa U-component of wind??
       l_u850,  & ! Output 850hPa U-component of wind??
       l_u800,  & ! Output 800hPa U-component of wind??
       l_u750,  & ! Output 750hPa U-component of wind??
       l_u700,  & ! Output 700hPa U-component of wind??
       l_u600,  & ! Output 600hPa U-component of wind??
       l_u500,  & ! Output 500hPa U-component of wind??
       l_u250,  & ! Output 250hPa U-component of wind??
       l_v1000, & ! Output 1000hPa V-component of wind??
       l_v950,  & ! Output 950hPa V-component of wind??
       l_v900,  & ! Output 900hPa V-component of wind??
       l_v850,  & ! Output 850hPa V-component of wind??
       l_v800,  & ! Output 800hPa V-component of wind??
       l_v750,  & ! Output 750hPa V-component of wind??
       l_v700,  & ! Output 700hPa V-component of wind??
       l_v600,  & ! Output 600hPa V-component of wind??
       l_v500,  & ! Output 500hPa V-component of wind??
       l_v250     ! Output 250hPa V-component of wind??
  real :: &
       z_soil
  namelist/nmlist/fileIN,fileOUT,z_soil,verbose,l_ivt,l_z0k,l_smois,&
       l_z1000,l_z950,l_z900,l_z850,l_z800,l_z750,l_z700,l_z600,l_z500,l_z250,&
       l_q1000,l_q950,l_q900,l_q850,l_q800,l_q750,l_q700,l_q600,l_q500,l_q250,&
       l_u1000,l_u950,l_u900,l_u850,l_u800,l_u750,l_u700,l_u600,l_u500,l_u250,&
       l_v1000,l_v950,l_v900,l_v850,l_v800,l_v750,l_v700,l_v600,l_v500,l_v250

  ! WRF fields
  integer :: &
       nTime,     & ! WRF file dimension: Number of times
       nLon,      & ! WRF file dimension: Number of longitudes
       nLat,      & ! WRF file dimension: Number of latitudes
       nLev,      & ! WRF file dimension: Number of vertical levels
       nLon_stag, & ! WRF file dimension: Number of longitudes (staggerd grid)
       nLat_stag, & ! WRF file dimension: Number of latitude (staggerd grid)
       nLev_stag, & ! WRF file dimension: Number of vertical levels (staggerd grid)
       DateStrLen,& ! WRF file dimension: String length for date.
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
  character(len=19),dimension(:),allocatable :: &
       times        ! WRF file time.
  integer,dimension(:),allocatable :: &
       year,      & ! WRF data year
       month,     & ! WRF data month
       day,       & ! WRF data day
       hour         ! WRF data hour
  
  ! Local fields
  real,dimension(:,:,:,:),allocatable :: &
       p,         & ! WRF pressue (computed from pp and pb)  (pa)
       hgt          ! WRF heights (computed from ph and phb) (m)
  real,dimension(:,:,:),allocatable :: &
       ivtU,      & ! WRF IVT (u-component)                  (kg/m/s)
       ivtV,      & ! WRF IVT (v-component)                  (kg/m/s)
       z0k,       & ! WRF Freezing level height              (m)
       z1000,     & ! WRF geopotential height @ 1000hPa      (m)
       z950,      & ! WRF geopotential height @ 950hPa       (m)
       z900,      & ! WRF geopotential height @ 900hPa       (m)
       z850,      & ! WRF geopotential height @ 850hPa       (m)
       z800,      & ! WRF geopotential height @ 800hPa       (m)
       z750,      & ! WRF geopotential height @ 750hPa       (m)
       z700,      & ! WRF geopotential height @ 700hPa       (m)
       z600,      & ! WRF geopotential height @ 600hPa       (m)
       z500,      & ! WRF geopotential height @ 500hPa       (m)
       z250,      & ! WRF geopotential height @ 2500hPa      (m)
       q1000,     & ! WRF specific-humidity @ 1000hPa        (kg/kg)
       q950,      & ! WRF specific-humidity @ 950hPa         (kg/kg)
       q900,      & ! WRF specific-humidity @ 900hPa         (kg/kg)
       q850,      & ! WRF specific-humidity @ 850hPa         (kg/kg)
       q800,      & ! WRF specific-humidity @ 800hPa         (kg/kg)
       q750,      & ! WRF specific-humidity @ 750hPa         (kg/kg)
       q700,      & ! WRF specific-humidity @ 700hPa         (kg/kg)
       q600,      & ! WRF specific-humidity @ 600hPa         (kg/kg)
       q500,      & ! WRF specific-humidity @ 500hPa         (kg/kg)
       q250,      & ! WRF specific-humidity @ 250hPa         (kg/kg)
       u1000,     & ! WRF u-wind @ 1000hPa                   (m/s)
       u950,      & ! WRF u-wind @ 950hPa                    (m/s)
       u900,      & ! WRF u-wind @ 900hPa                    (m/s)
       u850,      & ! WRF u-wind @ 850hPa                    (m/s)
       u800,      & ! WRF u-wind @ 800hPa                    (m/s)
       u750,      & ! WRF u-wind @ 750hPa                    (m/s)
       u700,      & ! WRF u-wind @ 700hPa                    (m/s)
       u600,      & ! WRF u-wind @ 600hPa                    (m/s)
       u500,      & ! WRF u-wind @ 500hPa                    (m/s)
       u250,      & ! WRF u-wind @ 250hPa                    (m/s)     
       v1000,     & ! WRF v-wind @ 1000hPa                   (m/s)
       v950,      & ! WRF v-wind @ 950hPa                    (m/s)
       v900,      & ! WRF v-wind @ 900hPa                    (m/s)
       v850,      & ! WRF v-wind @ 850hPa                    (m/s)
       v800,      & ! WRF v-wind @ 800hPa                    (m/s)
       v750,      & ! WRF v-wind @ 750hPa                    (m/s)
       v700,      & ! WRF v-wind @ 700hPa                    (m/s)
       v600,      & ! WRF v-wind @ 600hPa                    (m/s)
       v500,      & ! WRF v-wind @ 500hPa                    (m/s)
       v250,      & ! WRF v-wind @ 250hPa                    (m/s)
       soilMoisture ! WRF soil mositure @ z_soil             (m3/m3)
  
  integer,dimension(50) :: dimID,varIDout
  character(len=19) :: tempTime
  integer :: status,fileID,varID,ii,ij,ik,il
  real,dimension(:),allocatable :: ui,vi,temp,phm,a,b
  real :: wt1,dp,wt2
  integer,dimension(1) :: xi,xf
  logical :: &
       lread_QVAPOR = .false., &
       lread_U      = .false., &
       lread_V      = .false., &
       lread_PSFC   = .false., &
       lread_PB     = .false., &
       lread_P      = .false., &
       lread_PH     = .false., &
       lread_PHB    = .false., &
       lread_T      = .false., &
       lread_T00    = .false., &
       lread_HGT    = .false., &
       lread_ZS     = .false., &
       lread_SMOIS  = .false.

  ! #############################################################################
  ! A) Read in namelist
  ! #############################################################################
  open(10,file=trim(input_namelist),status='unknown')
  read(10,nml=nmlist)
  close(10)

  ! #############################################################################
  ! B) What fields need to be read in?
  ! #############################################################################
  if (l_ivt) then
     lread_QVAPOR = .true.
     lread_U      = .true.
     lread_V      = .true.
     lread_PSFC   = .true.
     lread_PB     = .true.
     lread_P      = .true.
  endif
  if (l_z0k) then
     lread_T      = .true.
     lread_T00    = .true.
     lread_HGT    = .true.
     lread_PB     = .true.
     lread_P      = .true.
     lread_PH     = .true.
     lread_PHB    = .true.
  endif
  if (l_smois) then
     lread_ZS     = .true.
     lread_SMOIS  = .true.
  endif
  if (any([l_z1000,l_z950,l_z900,l_z850,l_z800,l_z750,l_z700,l_z600,l_z500,l_z250])) then
     lread_PB     = .true.
     lread_P      = .true.
     lread_PH     = .true.
     lread_PHB    = .true.
  endif
  if (any([l_q1000,l_q950,l_q900,l_q850,l_q800,l_q750,l_q700,l_q600,l_q500,l_q250])) then
     lread_QVAPOR = .true.
     lread_PSFC   = .true.
     lread_PB     = .true.
     lread_P      = .true.
  endif
  if (any([l_u1000,l_u950,l_u900,l_u850,l_u800,l_u750,l_u700,l_u600,l_u500,l_u250])) then
     lread_U      = .true.
  endif
  if (any([l_v1000,l_v950,l_v900,l_v850,l_v800,l_v750,l_v700,l_v600,l_v500,l_v250])) then
     lread_V      = .true.
  endif
  ! #############################################################################
  ! C) Data ingest
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
  ! "DateStrLen"
  status = nf90_inq_dimid(fileID,"soil_layers_stag",dimID(9))
  if (status /= nf90_NoErr) print*,'ERROR: Dimension not present in file, DateStrLen'
  if (status == nf90_NoErr) then
     status = nf90_inquire_dimension(fileID,dimID(9),len=DateStrLen)
  endif
    
  ! 2) Read in fields
  ! Times
  status = nf90_inq_varid(fileID,"Times",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, Times'
  if (status == nf90_NoErr) then
     allocate(times(nTime),year(nTime),month(nTime),day(nTime),hour(nTime))
     status = nf90_get_var(fileID,varID,times)
     
     ! Pull out year, month, day and hour from timestamp
     do ij=1,nTime
        tempTime = times(ij)
        call str2int(tempTime(1:4),  year(ij), stat=status)
        call str2int(tempTime(6:7),  month(ij),stat=status)
        call str2int(tempTime(9:10), day(ij),  stat=status)
        call str2int(tempTime(12:13),hour(ij), stat=status)
     end do
  endif

  ! Longitude
  status = nf90_inq_varid(fileID,"XLONG",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, XLONG'
  if (status == nf90_NoErr) then
     allocate(lon(nLon,nLat,nTime))
     status = nf90_get_var(fileID,varID,lon)
  endif
  
  ! Latitude
  status = nf90_inq_varid(fileID,"XLAT",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, XLAT'
  if (status == nf90_NoErr) then
     allocate(lat(nLon,nLat,nTime))
     status = nf90_get_var(fileID,varID,lat)
  endif
  
  ! Longitude (staggered in zonal)
  status = nf90_inq_varid(fileID,"XLONG_U",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, XLONG_U'
  if (status == nf90_NoErr) then
     allocate(lon_u(nLon_stag,nLat,nTime))
     status = nf90_get_var(fileID,varID,lon_u)
  endif
  
  ! Latitude (staggered in zonal)
  status = nf90_inq_varid(fileID,"XLAT_U",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, XLAT_U'
  if (status == nf90_NoErr) then
     allocate(lat_u(nLon_stag,nLat,nTime))
     status = nf90_get_var(fileID,varID,lat_u)
  endif
  
  ! Longitude (staggered in meridional)
  status = nf90_inq_varid(fileID,"XLONG_V",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, XLONG_V'
  if (status == nf90_NoErr) then
     allocate(lon_v(nLon,nLat_stag,nTime))
     status = nf90_get_var(fileID,varID,lon_v)
  endif
  
  ! Latitude (staggered in meridional)
  status = nf90_inq_varid(fileID,"XLAT_V",varID)
  if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, XLAT_V'
  if (status == nf90_NoErr) then
     allocate(lat_v(nLon,nLat_stag,nTime))
     status = nf90_get_var(fileID,varID,lat_v)
  endif

  ! Only read in fields which are required for requested computations.
  if (lread_QVAPOR) then
     ! Water vapor mixing-ratio.
     status = nf90_inq_varid(fileID,"QVAPOR",varID)
     if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, QVAPOR'
     if (status == nf90_NoErr) then
        allocate(q(nLon,nLat,nLev,nTime))
        status = nf90_get_var(fileID,varID,q,count=(/nLon,nLat,nLev,nTime/))
     endif
  endif
  if (lread_U) then
     ! Eastward component of wind
     status = nf90_inq_varid(fileID,"U",varID)
     if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, U'
     if (status == nf90_NoErr) then
        allocate(u(nLon,nLat,nLev,nTime))
        status = nf90_get_var(fileID,varID,u,count=(/nLon_stag,nLat,nLev,nTime/))
     endif
  endif
  if (lread_V) then
     ! Northward component of wind
     status = nf90_inq_varid(fileID,"V",varID)
     if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, V'
     if (status == nf90_NoErr) then
        allocate(v(nLon,nLat,nLev,nTime))
        status = nf90_get_var(fileID,varID,v,count=(/nLon,nLat_stag,nLev,nTime/))
     endif
  endif
  if (lread_PSFC) then
     ! Surface pressure
     status = nf90_inq_varid(fileID,"PSFC",varID)
     if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, PSFC'
     if (status == nf90_NoErr) then
        allocate(psfc(nLon,nLat,nTime))
        status = nf90_get_var(fileID,varID,psfc,count=(/nLon,nLat,nTime/))
     endif
  endif
  if (lread_T00) then
     ! Base temperature
     status = nf90_inq_varid(fileID,"T00",varID)
     if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, T00'
     if (status == nf90_NoErr) then
        allocate(t00(nTime))
        status = nf90_get_var(fileID,varID,t00)
     endif
  endif
  if (lread_T) then
     ! Temperature
     status = nf90_inq_varid(fileID,"T",varID)
     if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, T'
     if (status == nf90_NoErr) then
        allocate(ta(nLon,nLat,nLev,nTime))
        status = nf90_get_var(fileID,varID,ta,count=(/nLon,nLat,nLev,nTime/))
     endif
  endif
  if (lread_HGT) then
     ! Terrain height
     status = nf90_inq_varid(fileID,"HGT",varID)
     if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, HGT'
     if (status == nf90_NoErr) then
        allocate(terrainZ(nLon,nLat,nTime))
        status = nf90_get_var(fileID,varID,terrainZ)
     endif
  endif
  if (lread_PB) then
     ! Base state (reference) pressure
     status = nf90_inq_varid(fileID,"PB",varID)
     if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, PB'
     if (status == nf90_NoErr) then
        allocate(pb(nLon,nLat,nLev,nTime))
        status = nf90_get_var(fileID,varID,pb,count=(/nLon,nLat,nLev,nTime/))
     endif
  endif
  if (lread_P) then
     ! Pertubation pressure
     status = nf90_inq_varid(fileID,"P",varID)
     if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, P'
     if (status == nf90_NoErr) then
        allocate(pp(nLon,nLat,nLev,nTime))
        status = nf90_get_var(fileID,varID,pp,count=(/nLon,nLat,nLev,nTime/))
     endif
  endif
  if (lread_PHB) then
     ! Base state (reference) height
     status = nf90_inq_varid(fileID,"PHB",varID)
     if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, PHB'
     if (status == nf90_NoErr) then
        allocate(phb(nLon,nLat,nLev_stag,nTime))
        status = nf90_get_var(fileID,varID,phb,count=(/nLon,nLat,nLev_stag,nTime/))
     endif
  endif
  if (lread_PH) then
     ! Pertubation height
     status = nf90_inq_varid(fileID,"PH",varID)
     if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, PH'
     if (status == nf90_NoErr) then
        allocate(ph(nLon,nLat,nLev_stag,nTime))
        status = nf90_get_var(fileID,varID,ph,count=(/nLon,nLat,nLev_stag,nTime/))
     endif
  endif
  if (lread_ZS) then
     ! Soil depth (staggered in vertical)
     status = nf90_inq_varid(fileID,"ZS",varID)
     if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, ZS'
     if (status == nf90_NoErr) then
        allocate(zs(nSoil_stag))
        status = nf90_get_var(fileID,varID,zs,count=(/nSoil_stag,1/))
     endif
  endif
  if (lread_SMOIS) then
     ! Soil moisture
     status = nf90_inq_varid(fileID,"SMOIS",varID)
     if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, SMOIS'
     if (status == nf90_NoErr) then
        allocate(smois(nLon,nLat,nSoil_stag,nTime))
        status = nf90_get_var(fileID,varID,smois)
     endif
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
  !if (l_ivt .or. l_z0k .or. l_z500 .or. l_q700) then
  if (allocated(pb) .and. allocated(pp)) then
     allocate(p(nLon,nLat,nLev,nTime))
     p   = pb+pp
  endif
  !if (l_z0k .or. l_z500) then
  if (allocated(ph) .and. allocated(phb)) then
     allocate(hgt(nLon,nLat,nLev_stag,nTime))
     hgt = (ph+phb)/9.8
  endif

  ! Allocate and initialize
  if (l_ivt) then
     allocate(ui(nLev),vi(nLev), ivtU(nLon,nLat,nTime), ivtV(nLon,nLat,nTime),a(nLev))
     ivtU(:,:,:) = 0.
     ivtV(:,:,:) = 0.
  endif
  if (l_z0k) then
     allocate(z0k(nLon,nLat,nTime),temp(nLev),phm(nLev))
     z0k(:,:,:)  = 0.
  endif
  if (l_smois) then
     allocate(b(nSoil_stag),soilMoisture(nLon,nLat,nTime))
     soilMoisture(:,:,:) = 0.
  endif
  if (l_z1000) then
     allocate(z1000(nLon,nLat,nTime))
     z1000(:,:,:) = 0.
  endif
  if (l_z950) then
     allocate(z950(nLon,nLat,nTime))
     z950(:,:,:) = 0.
  endif
  if (l_z900) then
     allocate(z900(nLon,nLat,nTime))
     z900(:,:,:) = 0.
  endif
  if (l_z850) then
     allocate(z850(nLon,nLat,nTime))
     z850(:,:,:) = 0.
  endif
  if (l_z800) then
     allocate(z800(nLon,nLat,nTime))
     z800(:,:,:) = 0.
  endif
  if (l_z750) then
     allocate(z750(nLon,nLat,nTime))
     z750(:,:,:) = 0.
  endif
  if (l_z700) then
     allocate(z700(nLon,nLat,nTime))
     z700(:,:,:) = 0.
  endif
  if (l_z600) then
     allocate(z600(nLon,nLat,nTime))
     z600(:,:,:) = 0.
  endif
  if (l_z500) then
     allocate(z500(nLon,nLat,nTime))
     z500(:,:,:) = 0.
  endif
  if (l_z250) then
     allocate(z250(nLon,nLat,nTime))
     z250(:,:,:) = 0.
  endif
  if (l_q1000) then
     allocate(q1000(nLon,nLat,nTime))
     q1000(:,:,:) = 0.
  endif
  if (l_q950) then
     allocate(q950(nLon,nLat,nTime))
     q950(:,:,:) = 0.
  endif
  if (l_q900) then
     allocate(q900(nLon,nLat,nTime))
     q900(:,:,:) = 0.
  endif
  if (l_q850) then
     allocate(q850(nLon,nLat,nTime))
     q850(:,:,:) = 0.
  endif
  if (l_q800) then
     allocate(q800(nLon,nLat,nTime))
     q800(:,:,:) = 0.
  endif
  if (l_q750) then
     allocate(q750(nLon,nLat,nTime))
     q750(:,:,:) = 0.
  endif
  if (l_q700) then
     allocate(q700(nLon,nLat,nTime))
     q700(:,:,:) = 0.
  endif
  if (l_q600) then
     allocate(q600(nLon,nLat,nTime))
     q600(:,:,:) = 0.
  endif
  if (l_q500) then
     allocate(q500(nLon,nLat,nTime))
     q500(:,:,:) = 0.
  endif
  if (l_q250) then
     allocate(q250(nLon,nLat,nTime))
     q250(:,:,:) = 0.
  endif
  if (l_u1000) then
     allocate(u1000(nLon,nLat,nTime))
     u1000(:,:,:) = 0.
  endif
  if (l_u950) then
     allocate(u950(nLon,nLat,nTime))
     u950(:,:,:) = 0.
  endif
  if (l_u900) then
     allocate(u900(nLon,nLat,nTime))
     u900(:,:,:) = 0.
  endif
  if (l_u850) then
     allocate(u850(nLon,nLat,nTime))
     u850(:,:,:) = 0.
  endif
  if (l_u800) then
     allocate(u800(nLon,nLat,nTime))
     u800(:,:,:) = 0.
  endif
  if (l_u750) then
     allocate(u750(nLon,nLat,nTime))
     u750(:,:,:) = 0.
  endif
  if (l_u700) then
     allocate(u700(nLon,nLat,nTime))
     u700(:,:,:) = 0.
  endif
  if (l_u600) then
     allocate(u600(nLon,nLat,nTime))
     u600(:,:,:) = 0.
  endif
  if (l_u500) then
     allocate(u500(nLon,nLat,nTime))
     u500(:,:,:) = 0.
  endif
  if (l_u250) then
     allocate(u250(nLon,nLat,nTime))
     u250(:,:,:) = 0.
  endif
  if (l_v1000) then
     allocate(v1000(nLon,nLat,nTime))
     v1000(:,:,:) = 0.
  endif
  if (l_v950) then
     allocate(v950(nLon,nLat,nTime))
     v950(:,:,:) = 0.
  endif
  if (l_v900) then
     allocate(v900(nLon,nLat,nTime))
     v900(:,:,:) = 0.
  endif
  if (l_v850) then
     allocate(v850(nLon,nLat,nTime))
     v850(:,:,:) = 0.
  endif
  if (l_v800) then
     allocate(v800(nLon,nLat,nTime))
     v800(:,:,:) = 0.
  endif
  if (l_v750) then
     allocate(v750(nLon,nLat,nTime))
     v750(:,:,:) = 0.
  endif
  if (l_v700) then
     allocate(v700(nLon,nLat,nTime))
     v700(:,:,:) = 0.
  endif
  if (l_v600) then
     allocate(v600(nLon,nLat,nTime))
     v600(:,:,:) = 0.
  endif
  if (l_v500) then
     allocate(v500(nLon,nLat,nTime))
     v500(:,:,:) = 0.
  endif
  if (l_v250) then
     allocate(v250(nLon,nLat,nTime))
     v250(:,:,:) = 0.
  endif
  
  ! Loop over all points/times and compute IVT and freezing level heights.
  if (verbose) print*,'Regridding velocity fields, computing IVT and freezing-level height: '
  do ii=1,nTime
     if (verbose) write(*,"(a12,i2,a4,i2)"),'   @Timestep ',ii,' of ',nTime
     do ij=1,nLon
        do ik=1,nLat
           if (any([l_u1000,l_u950,l_u900,l_u850,l_u800,l_u750,l_u700,l_u600,l_u500,l_u250]) .or. &
               any([l_v1000,l_v950,l_v900,l_v850,l_v800,l_v750,l_v700,l_v600,l_v500,l_v250]) .or. &
               l_ivt) then
              ! First, put winds on mass centered grid points, in WRF the velocity
              ! components are on staggered grids.
              ! Eastward component of wind (on mass-centered point)
              wt1 = (lon_u(ij+1,ik,ii)-lon(ij,ik,ii))/(lon_u(ij+1,ik,ii)-lon_u(ij,ik,ii))
              ui = wt1*u(ij,ik,:,ii)+(1-wt1)*u(ij+1,ik,:,ii)
              ! Northward component of wind (on mass-centered point)
              wt1 = (lat_v(ij,ik+1,ii)-lat(ij,ik,ii)) / (lat_v(ij,ik+1,ii)-lat_v(ij,ik,ii))
              vi = wt1*v(ij,ik,:,ii)+(1-wt1)*v(ij,ik+1,:,ii)
           endif
           
           ! ######################################################################
           ! Compute integrated vapor transport (IVT)
           ! ######################################################################
           if (l_ivt) then
              do il=1,nLev
                 ! Compute pressure change across layer.
                 if (il .eq. 1)                    dp = psfc(ij,ij,ii)-sum(p(ij,ik,il:il+1,ii))/2.              ! Bottommost level
                 if (il .ne. 1 .and. il .ne. nLev) dp = sum(p(ij,ik,il-1:il,ii))/2.-sum(p(ij,ik,il:il+1,ii))/2. ! Middle levels
                 if (il .eq. nLev)                 dp = 0.                                                      ! Top level
                 ! Integrate vertically.
                 ivtU(ij,ik,ii) = ivtU(ij,ik,ii) + ui(il)*q(ij,ik,il,ii)*dp/9.8
                 ivtV(ij,ik,ii) = ivtV(ij,ik,ii) + vi(il)*q(ij,ik,il,ii)*dp/9.8
              enddo
           endif

           ! ######################################################################
           ! Compute freezing-level height.
           ! ######################################################################
           if (l_z0k) then
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
           endif
           
           ! ###################################################################
           ! Compute soil moisture at z_soil.
           ! ###################################################################
           if (l_smois) then
              b  = zs-z_soil
              xf = minloc(b,b .gt. 0)
              xi = xf-1
              if (xi(1) .lt. 1) then
                 write(*,"(a37,f5.2,a17)"),'ERROR: Requested soil moisture depth, ',z_soil,' is out of bounds.'
              endif
              wt2 = (zs(xf(1))-z_soil)/(zs(xf(1))-zs(xi(1)))
              soilMoisture(ij,ik,ii) = smois(ij,ik,xi(1),ii)*wt2+smois(ij,ik,xf(1),ii)*(1-wt2)
           endif

           ! ###################################################################
           ! Compute u-winds @ various levels
           ! ###################################################################
           if (l_u1000) then
              xi = minloc(p(ij,ik,:,ii)-100000.,p(ij,ik,:,ii)-100000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-100000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              u1000(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
           endif
           if (l_u950) then
              xi = minloc(p(ij,ik,:,ii)-95000.,p(ij,ik,:,ii)-95000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-95000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              u950(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
           endif
           if (l_u900) then
              xi = minloc(p(ij,ik,:,ii)-90000.,p(ij,ik,:,ii)-90000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-90000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              u900(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
           endif
           if (l_u850) then
              xi = minloc(p(ij,ik,:,ii)-85000.,p(ij,ik,:,ii)-85000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-85000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              u850(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
           endif
           if (l_u800) then
              xi = minloc(p(ij,ik,:,ii)-80000.,p(ij,ik,:,ii)-80000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-80000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              u800(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
           endif
           if (l_u750) then
              xi = minloc(p(ij,ik,:,ii)-75000.,p(ij,ik,:,ii)-75000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-75000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              u750(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
           endif
           if (l_u700) then
              xi = minloc(p(ij,ik,:,ii)-70000.,p(ij,ik,:,ii)-70000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-70000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              u700(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
           endif
           if (l_u600) then
              xi = minloc(p(ij,ik,:,ii)-60000.,p(ij,ik,:,ii)-60000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-60000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              u600(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
           endif
           if (l_u500) then
              xi = minloc(p(ij,ik,:,ii)-50000.,p(ij,ik,:,ii)-50000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-50000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              u500(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
           endif
           if (l_u250) then
              xi = minloc(p(ij,ik,:,ii)-25000.,p(ij,ik,:,ii)-25000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-25000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              u250(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
           endif

           ! ###################################################################
           ! Compute v-winds @ various levels
           ! ###################################################################
           if (l_v1000) then
              xi = minloc(p(ij,ik,:,ii)-100000.,p(ij,ik,:,ii)-100000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-100000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              v1000(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
           endif
           if (l_v950) then
              xi = minloc(p(ij,ik,:,ii)-95000.,p(ij,ik,:,ii)-95000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-95000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              v950(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
           endif
           if (l_v900) then
              xi = minloc(p(ij,ik,:,ii)-90000.,p(ij,ik,:,ii)-90000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-90000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              v900(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
           endif
           if (l_v850) then
              xi = minloc(p(ij,ik,:,ii)-85000.,p(ij,ik,:,ii)-85000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-85000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              v850(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
           endif
           if (l_v800) then
              xi = minloc(p(ij,ik,:,ii)-80000.,p(ij,ik,:,ii)-80000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-80000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              v800(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
           endif
           if (l_v750) then
              xi = minloc(p(ij,ik,:,ii)-75000.,p(ij,ik,:,ii)-75000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-75000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              v750(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
           endif
           if (l_v700) then
              xi = minloc(p(ij,ik,:,ii)-70000.,p(ij,ik,:,ii)-70000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-70000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              v700(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
           endif
           if (l_v600) then
              xi = minloc(p(ij,ik,:,ii)-60000.,p(ij,ik,:,ii)-60000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-60000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              v600(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
           endif
           if (l_v500) then
              xi = minloc(p(ij,ik,:,ii)-50000.,p(ij,ik,:,ii)-50000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-50000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              v500(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
           endif
           if (l_v250) then
              xi = minloc(p(ij,ik,:,ii)-25000.,p(ij,ik,:,ii)-25000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-25000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              v250(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
           endif
           
           ! ###################################################################
           ! Compute geopotential heights @ various levels
           ! ###################################################################
           if (l_z1000) then
              xi = minloc(p(ij,ik,:,ii)-100000.,p(ij,ik,:,ii)-100000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-100000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              z1000(ij,ik,ii) = hgt(ij,ik,xi(1),ii)*wt2 + hgt(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_z950) then
              xi = minloc(p(ij,ik,:,ii)-95000.,p(ij,ik,:,ii)-95000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-95000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              z950(ij,ik,ii) = hgt(ij,ik,xi(1),ii)*wt2 + hgt(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_z900) then
              xi = minloc(p(ij,ik,:,ii)-90000.,p(ij,ik,:,ii)-90000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-90000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              z900(ij,ik,ii) = hgt(ij,ik,xi(1),ii)*wt2 + hgt(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_z850) then
              xi = minloc(p(ij,ik,:,ii)-85000.,p(ij,ik,:,ii)-85000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-85000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              z850(ij,ik,ii) = hgt(ij,ik,xi(1),ii)*wt2 + hgt(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_z800) then
              xi = minloc(p(ij,ik,:,ii)-80000.,p(ij,ik,:,ii)-80000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-80000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              z800(ij,ik,ii) = hgt(ij,ik,xi(1),ii)*wt2 + hgt(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_z750) then
              xi = minloc(p(ij,ik,:,ii)-75000.,p(ij,ik,:,ii)-75000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-75000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              z750(ij,ik,ii) = hgt(ij,ik,xi(1),ii)*wt2 + hgt(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_z700) then
              xi = minloc(p(ij,ik,:,ii)-70000.,p(ij,ik,:,ii)-70000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-70000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              z700(ij,ik,ii) = hgt(ij,ik,xi(1),ii)*wt2 + hgt(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_z600) then
              xi = minloc(p(ij,ik,:,ii)-60000.,p(ij,ik,:,ii)-60000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-60000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              z600(ij,ik,ii) = hgt(ij,ik,xi(1),ii)*wt2 + hgt(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_z500) then
              xi = minloc(p(ij,ik,:,ii)-50000.,p(ij,ik,:,ii)-50000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-50000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              z500(ij,ik,ii) = hgt(ij,ik,xi(1),ii)*wt2 + hgt(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_z250) then
              xi = minloc(p(ij,ik,:,ii)-25000.,p(ij,ik,:,ii)-25000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-25000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              z250(ij,ik,ii) = hgt(ij,ik,xi(1),ii)*wt2 + hgt(ij,ik,xf(1),ii)*(1-wt2)
           endif

           ! ###################################################################
           ! Compute specific-humidity @ various levels
           ! ###################################################################         
           if (l_q1000) then
              xi = minloc(p(ij,ik,:,ii)-100000.,p(ij,ik,:,ii)-100000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-100000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              q1000(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_q950) then
              xi = minloc(p(ij,ik,:,ii)-95000.,p(ij,ik,:,ii)-95000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-95000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              q950(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_q900) then
              xi = minloc(p(ij,ik,:,ii)-90000.,p(ij,ik,:,ii)-90000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-90000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              q900(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_q850) then
              xi = minloc(p(ij,ik,:,ii)-85000.,p(ij,ik,:,ii)-85000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-85000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              q850(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_q800) then
              xi = minloc(p(ij,ik,:,ii)-80000.,p(ij,ik,:,ii)-80000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-80000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              q800(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_q750) then
              xi = minloc(p(ij,ik,:,ii)-75000.,p(ij,ik,:,ii)-75000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-75000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              q750(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_q700) then
              xi = minloc(p(ij,ik,:,ii)-70000.,p(ij,ik,:,ii)-70000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-70000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              q700(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_q600) then
              xi = minloc(p(ij,ik,:,ii)-60000.,p(ij,ik,:,ii)-60000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-60000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              q600(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_q500) then
              xi = minloc(p(ij,ik,:,ii)-50000.,p(ij,ik,:,ii)-50000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-50000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              q500(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
           endif
           if (l_q250) then
              xi = minloc(p(ij,ik,:,ii)-25000.,p(ij,ik,:,ii)-25000. .gt. 0)
              xf = xi + 1
              wt2 = (p(ij,ik,xf(1),ii)-25000.)/(p(ij,ik,xf(1),ii)-p(ij,ik,xi(1),ii))
              q250(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
           endif

        enddo ! Latitude
     enddo    ! Longitude
  enddo       ! Time

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

  ! Define variables and add attributes.
  ! Time
  status = nf90_def_var(fileID,'Year',nf90_int, (/dimID(1)/),varIDout(10))
  if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Year'
  status = nf90_def_var(fileID,'Month',nf90_int, (/dimID(1)/),varIDout(11))
  if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Month'
  status = nf90_def_var(fileID,'Day',nf90_int, (/dimID(1)/),varIDout(12))
  if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Day'
  status = nf90_def_var(fileID,'Hour',nf90_int, (/dimID(1)/),varIDout(13))
  if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Hour'
  
  ! Latitude
  status = nf90_def_var(fileID,'XLAT',nf90_float, (/dimID(2),dimID(3)/),varIDout(1))
  if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, XLAT'
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
  ! Longitude
  status = nf90_def_var(fileID,'XLONG',nf90_float,(/dimID(2),dimID(3)/),varIDout(2))
  if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, XLONG'
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
  ! IVT
  if (l_ivt) then
     ! x-component
     status = nf90_def_var(fileID,'IVTU',nf90_float, (/dimID(2),dimID(3),dimID(1)/),varIDout(3))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, IVTU'
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
     ! y-component
     status = nf90_def_var(fileID,'IVTV',nf90_float, (/dimID(2),dimID(3),dimID(1)/),varIDout(4))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, IVTV'
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
  endif
  ! Freezing-level height
  if (l_z0k) then
     status = nf90_def_var(fileID,'Z0K',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(5))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Z0K'
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
  endif
  ! Soil mositure @ z_soil
  if (l_smois) then
     status = nf90_def_var(fileID,'SMOIS',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(6))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, SMOIS'
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
  endif

  ! 1000hPa geopotential height
  if (l_z1000) then
     status = nf90_def_var(fileID,'Z1000',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(14))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Z1000'
     status = nf90_put_att(fileID,varIDout(14),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Z1000'
     status = nf90_put_att(fileID,varIDout(14),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Z1000'
     status = nf90_put_att(fileID,varIDout(14),"description","GEOPOTENTIAL HEIGHT @ 1000hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Z1000'
     status = nf90_put_att(fileID,varIDout(14),"units","m")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Z1000'
     status = nf90_put_att(fileID,varIDout(14),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Z1000'
  endif
  
  ! 950hPa geopotential height
  if (l_z950) then
     status = nf90_def_var(fileID,'Z950',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(15))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Z950'
     status = nf90_put_att(fileID,varIDout(15),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Z950'
     status = nf90_put_att(fileID,varIDout(15),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Z950'
     status = nf90_put_att(fileID,varIDout(15),"description","GEOPOTENTIAL HEIGHT @ 950hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Z950'
     status = nf90_put_att(fileID,varIDout(15),"units","m")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Z950'
     status = nf90_put_att(fileID,varIDout(15),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Z950'
  endif
  
  ! 900hPa geopotential height
  if (l_z900) then
     status = nf90_def_var(fileID,'Z900',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(16))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Z900'
     status = nf90_put_att(fileID,varIDout(16),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Z900'
     status = nf90_put_att(fileID,varIDout(16),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Z900'
     status = nf90_put_att(fileID,varIDout(16),"description","GEOPOTENTIAL HEIGHT @ 900hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Z900'
     status = nf90_put_att(fileID,varIDout(16),"units","m")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Z900'
     status = nf90_put_att(fileID,varIDout(16),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Z900'
  endif
  
  ! 850hPa geopotential height
  if (l_z850) then
     status = nf90_def_var(fileID,'Z850',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(17))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Z850'
     status = nf90_put_att(fileID,varIDout(17),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Z850'
     status = nf90_put_att(fileID,varIDout(17),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Z850'
     status = nf90_put_att(fileID,varIDout(17),"description","GEOPOTENTIAL HEIGHT @ 850hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Z850'
     status = nf90_put_att(fileID,varIDout(17),"units","m")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Z850'
     status = nf90_put_att(fileID,varIDout(17),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Z850'
  endif
  
  ! 800hPa geopotential height
  if (l_z800) then
     status = nf90_def_var(fileID,'Z800',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(18))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Z800'
     status = nf90_put_att(fileID,varIDout(18),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Z800'
     status = nf90_put_att(fileID,varIDout(18),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Z800'
     status = nf90_put_att(fileID,varIDout(18),"description","GEOPOTENTIAL HEIGHT @ 800hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Z800'
     status = nf90_put_att(fileID,varIDout(18),"units","m")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Z800'
     status = nf90_put_att(fileID,varIDout(18),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Z800'
  endif
  
  ! 750hPa geopotential height
  if (l_z750) then
     status = nf90_def_var(fileID,'Z750',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(19))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Z750'
     status = nf90_put_att(fileID,varIDout(19),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Z750'
     status = nf90_put_att(fileID,varIDout(19),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Z750'
     status = nf90_put_att(fileID,varIDout(19),"description","GEOPOTENTIAL HEIGHT @ 750hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Z750'
     status = nf90_put_att(fileID,varIDout(19),"units","m")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Z750'
     status = nf90_put_att(fileID,varIDout(19),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Z750'
  endif
  
  ! 700hPa geopotential height
  if (l_z700) then
     status = nf90_def_var(fileID,'Z700',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(20))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Z700'
     status = nf90_put_att(fileID,varIDout(20),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Z700'
     status = nf90_put_att(fileID,varIDout(20),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Z700'
     status = nf90_put_att(fileID,varIDout(20),"description","GEOPOTENTIAL HEIGHT @ 700hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Z700'
     status = nf90_put_att(fileID,varIDout(20),"units","m")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Z700'
     status = nf90_put_att(fileID,varIDout(20),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Z700'
  endif
  
  ! 600hPa geopotential height
  if (l_z600) then
     status = nf90_def_var(fileID,'Z600',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(21))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Z600'
     status = nf90_put_att(fileID,varIDout(21),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Z600'
     status = nf90_put_att(fileID,varIDout(21),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Z600'
     status = nf90_put_att(fileID,varIDout(21),"description","GEOPOTENTIAL HEIGHT @ 600hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Z600'
     status = nf90_put_att(fileID,varIDout(21),"units","m")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Z600'
     status = nf90_put_att(fileID,varIDout(21),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Z600'
  endif
  
  ! 500hPa geopotential height
  if (l_z500) then
     status = nf90_def_var(fileID,'Z500',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(22))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Z500'
     status = nf90_put_att(fileID,varIDout(22),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Z500'
     status = nf90_put_att(fileID,varIDout(22),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Z500'
     status = nf90_put_att(fileID,varIDout(22),"description","GEOPOTENTIAL HEIGHT @ 500hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Z500'
     status = nf90_put_att(fileID,varIDout(22),"units","m")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Z500'
     status = nf90_put_att(fileID,varIDout(22),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Z500'
  endif
  
  ! 250hPa geopotential height
  if (l_z250) then
     status = nf90_def_var(fileID,'Z250',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(23))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Z250'
     status = nf90_put_att(fileID,varIDout(23),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Z250'
     status = nf90_put_att(fileID,varIDout(23),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Z250'
     status = nf90_put_att(fileID,varIDout(23),"description","GEOPOTENTIAL HEIGHT @ 250hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Z250'
     status = nf90_put_att(fileID,varIDout(23),"units","m")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Z250'
     status = nf90_put_att(fileID,varIDout(23),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Z250'
  endif

  ! 1000hPa specific-humidity
  if (l_q1000) then
     status = nf90_def_var(fileID,'Q1000',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(24))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Q1000'
     status = nf90_put_att(fileID,varIDout(24),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Q1000'
     status = nf90_put_att(fileID,varIDout(24),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Q1000'
     status = nf90_put_att(fileID,varIDout(24),"description","SPECIFIC HUMIDITY @ 1000hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Q1000'
     status = nf90_put_att(fileID,varIDout(24),"units","kg/kg")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Q1000'
     status = nf90_put_att(fileID,varIDout(24),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Q1000'
  endif
  
  ! 950hPa specific-humidity
  if (l_q950) then
     status = nf90_def_var(fileID,'Q950',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(25))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Q950'
     status = nf90_put_att(fileID,varIDout(25),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Q950'
     status = nf90_put_att(fileID,varIDout(25),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Q950'
     status = nf90_put_att(fileID,varIDout(25),"description","SPECIFIC HUMIDITY @ 950hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Q950'
     status = nf90_put_att(fileID,varIDout(25),"units","kg/kg")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Q950'
     status = nf90_put_att(fileID,varIDout(25),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Q950'
  endif
  
  ! 900hPa specific-humidity
  if (l_q900) then
     status = nf90_def_var(fileID,'Q900',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(26))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Q900'
     status = nf90_put_att(fileID,varIDout(26),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Q900'
     status = nf90_put_att(fileID,varIDout(26),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Q900'
     status = nf90_put_att(fileID,varIDout(26),"description","SPECIFIC HUMIDITY @ 900hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Q900'
     status = nf90_put_att(fileID,varIDout(26),"units","kg/kg")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Q900'
     status = nf90_put_att(fileID,varIDout(26),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Q900'
  endif
  
  ! 850hPa specific-humidity
  if (l_q850) then
     status = nf90_def_var(fileID,'Q850',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(27))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Q850'
     status = nf90_put_att(fileID,varIDout(27),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Q850'
     status = nf90_put_att(fileID,varIDout(27),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Q850'
     status = nf90_put_att(fileID,varIDout(27),"description","SPECIFIC HUMIDITY @ 850hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Q850'
     status = nf90_put_att(fileID,varIDout(27),"units","kg/kg")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Q850'
     status = nf90_put_att(fileID,varIDout(27),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Q850'
  endif
  
  ! 800hPa specific-humidity
  if (l_q800) then
     status = nf90_def_var(fileID,'Q800',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(28))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Q800'
     status = nf90_put_att(fileID,varIDout(28),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Q800'
     status = nf90_put_att(fileID,varIDout(28),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Q800'
     status = nf90_put_att(fileID,varIDout(28),"description","SPECIFIC HUMIDITY @ 800hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Q800'
     status = nf90_put_att(fileID,varIDout(28),"units","kg/kg")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Q800'
     status = nf90_put_att(fileID,varIDout(28),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Q800'
  endif
  
  ! 750hPa specific-humidity
  if (l_q750) then
     status = nf90_def_var(fileID,'Q750',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(29))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Q750'
     status = nf90_put_att(fileID,varIDout(29),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Q750'
     status = nf90_put_att(fileID,varIDout(29),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Q750'
     status = nf90_put_att(fileID,varIDout(29),"description","SPECIFIC HUMIDITY @ 750hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Q750'
     status = nf90_put_att(fileID,varIDout(29),"units","kg/kg")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Q750'
     status = nf90_put_att(fileID,varIDout(29),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Q750'
  endif
  
  ! 700hPa specific-humidity
  if (l_q700) then
     status = nf90_def_var(fileID,'Q700',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(30))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Q700'
     status = nf90_put_att(fileID,varIDout(30),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Q700'
     status = nf90_put_att(fileID,varIDout(30),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Q700'
     status = nf90_put_att(fileID,varIDout(30),"description","SPECIFIC HUMIDITY @ 700hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Q700'
     status = nf90_put_att(fileID,varIDout(30),"units","kg/kg")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Q700'
     status = nf90_put_att(fileID,varIDout(30),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Q700'
  endif
  
  ! 600hPa specific-humidity
  if (l_q600) then
     status = nf90_def_var(fileID,'Q600',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(31))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Q600'
     status = nf90_put_att(fileID,varIDout(31),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Q600'
     status = nf90_put_att(fileID,varIDout(31),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Q600'
     status = nf90_put_att(fileID,varIDout(31),"description","SPECIFIC HUMIDITY @ 600hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Q600'
     status = nf90_put_att(fileID,varIDout(31),"units","kg/kg")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Q600'
     status = nf90_put_att(fileID,varIDout(31),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Q600'
  endif
  
  ! 500hPa specific-humidity
  if (l_q500) then
     status = nf90_def_var(fileID,'Q500',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(32))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Q500'
     status = nf90_put_att(fileID,varIDout(32),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Q500'
     status = nf90_put_att(fileID,varIDout(32),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Q500'
     status = nf90_put_att(fileID,varIDout(32),"description","SPECIFIC HUMIDITY @ 500hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Q500'
     status = nf90_put_att(fileID,varIDout(32),"units","kg/kg")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Q500'
     status = nf90_put_att(fileID,varIDout(32),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Q500'
  endif
  
  ! 250hPa specific-humidity
  if (l_q250) then
     status = nf90_def_var(fileID,'Q250',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(33))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, Q250'
     status = nf90_put_att(fileID,varIDout(33),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable Q250'
     status = nf90_put_att(fileID,varIDout(33),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable Q250'
     status = nf90_put_att(fileID,varIDout(33),"description","SPECIFIC HUMIDITY @ 250hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable Q250'
     status = nf90_put_att(fileID,varIDout(33),"units","kg/kg")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable Q250'
     status = nf90_put_att(fileID,varIDout(33),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable Q250'
  endif

  ! 1000hPa u-winds
  if (l_u1000) then
     status = nf90_def_var(fileID,'U1000',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(34))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, U1000'
     status = nf90_put_att(fileID,varIDout(34),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable U1000'
     status = nf90_put_att(fileID,varIDout(34),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable U1000'
     status = nf90_put_att(fileID,varIDout(34),"description","ZONAL WIND @ 1000hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable U1000'
     status = nf90_put_att(fileID,varIDout(34),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable U1000'
     status = nf90_put_att(fileID,varIDout(34),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable U1000'
  endif
  
  ! 950hPa u-winds
  if (l_u950) then
     status = nf90_def_var(fileID,'U950',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(35))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, U950'
     status = nf90_put_att(fileID,varIDout(35),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable U950'
     status = nf90_put_att(fileID,varIDout(35),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable U950'
     status = nf90_put_att(fileID,varIDout(35),"description","ZONAL WIND @ 950hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable U950'
     status = nf90_put_att(fileID,varIDout(35),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable U950'
     status = nf90_put_att(fileID,varIDout(35),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable U950'
  endif
  
  ! 900hPa u-winds
  if (l_u900) then
     status = nf90_def_var(fileID,'U900',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(36))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, U900'
     status = nf90_put_att(fileID,varIDout(36),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable U900'
     status = nf90_put_att(fileID,varIDout(36),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable U900'
     status = nf90_put_att(fileID,varIDout(36),"description","ZONAL WIND @ 900hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable U900'
     status = nf90_put_att(fileID,varIDout(36),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable U900'
     status = nf90_put_att(fileID,varIDout(36),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable U900'
  endif
  
  ! 850hPa u-winds
  if (l_u850) then
     status = nf90_def_var(fileID,'U850',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(37))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, U850'
     status = nf90_put_att(fileID,varIDout(37),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable U850'
     status = nf90_put_att(fileID,varIDout(37),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable U850'
     status = nf90_put_att(fileID,varIDout(37),"description","ZONAL WIND @ 850hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable U850'
     status = nf90_put_att(fileID,varIDout(37),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable U850'
     status = nf90_put_att(fileID,varIDout(37),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable U850'
  endif
  
  ! 800hPa u-winds
  if (l_u800) then
     status = nf90_def_var(fileID,'U800',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(38))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, U800'
     status = nf90_put_att(fileID,varIDout(38),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable U800'
     status = nf90_put_att(fileID,varIDout(38),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable U800'
     status = nf90_put_att(fileID,varIDout(38),"description","ZONAL WIND @ 800hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable U800'
     status = nf90_put_att(fileID,varIDout(38),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable U800'
     status = nf90_put_att(fileID,varIDout(38),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable U800'
  endif
  
  ! 750hPa u-winds
  if (l_u750) then
     status = nf90_def_var(fileID,'U750',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(39))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, U750'
     status = nf90_put_att(fileID,varIDout(39),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable U750'
     status = nf90_put_att(fileID,varIDout(39),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable U750'
     status = nf90_put_att(fileID,varIDout(39),"description","ZONAL WIND @ 750hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable U750'
     status = nf90_put_att(fileID,varIDout(39),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable U750'
     status = nf90_put_att(fileID,varIDout(39),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable U750'
  endif
  
  ! 700hPa u-winds
  if (l_u700) then
     status = nf90_def_var(fileID,'U700',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(40))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, U700'
     status = nf90_put_att(fileID,varIDout(40),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable U700'
     status = nf90_put_att(fileID,varIDout(40),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable U700'
     status = nf90_put_att(fileID,varIDout(40),"description","ZONAL WIND @ 700hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable U700'
     status = nf90_put_att(fileID,varIDout(40),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable U700'
     status = nf90_put_att(fileID,varIDout(40),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable U700'
  endif
  
  ! 600hPa u-winds
  if (l_u600) then
     status = nf90_def_var(fileID,'U600',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(41))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, U600'
     status = nf90_put_att(fileID,varIDout(41),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable U600'
     status = nf90_put_att(fileID,varIDout(41),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable U600'
     status = nf90_put_att(fileID,varIDout(41),"description","ZONAL WIND @ 600hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable U600'
     status = nf90_put_att(fileID,varIDout(41),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable U600'
     status = nf90_put_att(fileID,varIDout(41),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable U600'
  endif
  
  ! 500hPa u-winds
  if (l_u500) then
     status = nf90_def_var(fileID,'U500',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(42))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, U500'
     status = nf90_put_att(fileID,varIDout(42),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable U500'
     status = nf90_put_att(fileID,varIDout(42),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable U500'
     status = nf90_put_att(fileID,varIDout(42),"description","ZONAL WIND @ 500hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable U500'
     status = nf90_put_att(fileID,varIDout(42),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable U500'
     status = nf90_put_att(fileID,varIDout(42),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable U500'
  endif
  
  ! 250hPa u-winds
  if (l_u250) then
     status = nf90_def_var(fileID,'U250',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(43))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, U250'
     status = nf90_put_att(fileID,varIDout(43),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable U250'
     status = nf90_put_att(fileID,varIDout(43),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable U250'
     status = nf90_put_att(fileID,varIDout(43),"description","ZONAL WIND @ 250hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable U250'
     status = nf90_put_att(fileID,varIDout(43),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable U250'
     status = nf90_put_att(fileID,varIDout(43),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable U250'
  endif

















  



  ! 1000hPa v-winds
  if (l_v1000) then
     status = nf90_def_var(fileID,'V1000',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(44))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, V1000'
     status = nf90_put_att(fileID,varIDout(44),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable V1000'
     status = nf90_put_att(fileID,varIDout(44),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable V1000'
     status = nf90_put_att(fileID,varIDout(44),"description","MERIDIONAL WIND @ 1000hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable V1000'
     status = nf90_put_att(fileID,varIDout(44),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable V1000'
     status = nf90_put_att(fileID,varIDout(44),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable V1000'
  endif
  
  ! 950hPa v-winds
  if (l_v950) then
     status = nf90_def_var(fileID,'V950',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(45))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, V950'
     status = nf90_put_att(fileID,varIDout(45),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable V950'
     status = nf90_put_att(fileID,varIDout(45),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable V950'
     status = nf90_put_att(fileID,varIDout(45),"description","MERIDIONAL WIND @ 950hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable V950'
     status = nf90_put_att(fileID,varIDout(45),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable V950'
     status = nf90_put_att(fileID,varIDout(45),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable V950'
  endif
  
  ! 900hPa v-winds
  if (l_v900) then
     status = nf90_def_var(fileID,'V900',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(46))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, V900'
     status = nf90_put_att(fileID,varIDout(46),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable V900'
     status = nf90_put_att(fileID,varIDout(46),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable V900'
     status = nf90_put_att(fileID,varIDout(46),"description","MERIDIONAL WIND @ 900hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable V900'
     status = nf90_put_att(fileID,varIDout(46),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable V900'
     status = nf90_put_att(fileID,varIDout(46),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable V900'
  endif
  
  ! 850hPa v-winds
  if (l_v850) then
     status = nf90_def_var(fileID,'V850',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(47))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, V850'
     status = nf90_put_att(fileID,varIDout(47),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable V850'
     status = nf90_put_att(fileID,varIDout(47),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable V850'
     status = nf90_put_att(fileID,varIDout(47),"description","MERIDIONAL WIND @ 850hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable V850'
     status = nf90_put_att(fileID,varIDout(47),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable V850'
     status = nf90_put_att(fileID,varIDout(47),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable V850'
  endif
  
  ! 800hPa v-winds
  if (l_v800) then
     status = nf90_def_var(fileID,'V800',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(48))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, V800'
     status = nf90_put_att(fileID,varIDout(48),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable V800'
     status = nf90_put_att(fileID,varIDout(48),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable V800'
     status = nf90_put_att(fileID,varIDout(48),"description","MERIDIONAL WIND @ 800hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable V800'
     status = nf90_put_att(fileID,varIDout(48),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable V800'
     status = nf90_put_att(fileID,varIDout(48),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable V800'
  endif
  
  ! 750hPa v-winds
  if (l_v750) then
     status = nf90_def_var(fileID,'V750',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(49))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, V750'
     status = nf90_put_att(fileID,varIDout(49),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable V750'
     status = nf90_put_att(fileID,varIDout(49),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable V750'
     status = nf90_put_att(fileID,varIDout(49),"description","MERIDIONAL WIND @ 750hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable V750'
     status = nf90_put_att(fileID,varIDout(49),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable V750'
     status = nf90_put_att(fileID,varIDout(49),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable V750'
  endif
  
  ! 700hPa v-winds
  if (l_v700) then
     status = nf90_def_var(fileID,'V700',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(50))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, U700'
     status = nf90_put_att(fileID,varIDout(50),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable V700'
     status = nf90_put_att(fileID,varIDout(50),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable V700'
     status = nf90_put_att(fileID,varIDout(50),"description","MERIDIONAL WIND @ 700hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable V700'
     status = nf90_put_att(fileID,varIDout(50),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable V700'
     status = nf90_put_att(fileID,varIDout(50),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable V700'
  endif
  
  ! 600hPa v-winds
  if (l_v600) then
     status = nf90_def_var(fileID,'V600',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(51))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, V600'
     status = nf90_put_att(fileID,varIDout(51),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable V600'
     status = nf90_put_att(fileID,varIDout(51),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable V600'
     status = nf90_put_att(fileID,varIDout(51),"description","MERIDIONAL WIND @ 600hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable V600'
     status = nf90_put_att(fileID,varIDout(51),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable V600'
     status = nf90_put_att(fileID,varIDout(51),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable V600'
  endif
  
  ! 500hPa v-winds
  if (l_v500) then
     status = nf90_def_var(fileID,'V500',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(52))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, V500'
     status = nf90_put_att(fileID,varIDout(52),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable V500'
     status = nf90_put_att(fileID,varIDout(52),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable V500'
     status = nf90_put_att(fileID,varIDout(52),"description","MERIDIONAL WIND @ 500hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable V500'
     status = nf90_put_att(fileID,varIDout(52),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable V500'
     status = nf90_put_att(fileID,varIDout(52),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable V500'
  endif
  
  ! 250hPa v-winds
  if (l_v250) then
     status = nf90_def_var(fileID,'V250',nf90_float,  (/dimID(2),dimID(3),dimID(1)/),varIDout(53))
     if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable, V250'
     status = nf90_put_att(fileID,varIDout(53),"FieldType",104)
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable V250'
     status = nf90_put_att(fileID,varIDout(53),"MemoryOrder","XY")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable V250'
     status = nf90_put_att(fileID,varIDout(53),"description","MERIDIONAL WIND @ 250hPa")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable V250'
     status = nf90_put_att(fileID,varIDout(53),"units","m/s")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable V250'
     status = nf90_put_att(fileID,varIDout(53),"stagger","")
     if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable V250'
  endif
  













  
  
  ! Exit define mode
  status = nf90_enddef(fileID)
  if (status /= nf90_NoErr) print*,'ERROR: Failure exiting define mode'

  ! Populate file
  status = nf90_put_var(fileID,varIDout(1),lat(:,:,1))
  if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, XLAT'
  status = nf90_put_var(fileID,varIDout(2),lon(:,:,1))
  if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, XLONG'
  status = nf90_put_var(fileID,varIDout(10),year)
  if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Year'
  status = nf90_put_var(fileID,varIDout(11),month)
  if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Month'
  status = nf90_put_var(fileID,varIDout(12),day)
  if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Day'
  status = nf90_put_var(fileID,varIDout(13),hour)
  if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Hour'
  if (l_ivt) then
     status = nf90_put_var(fileID,varIDout(3),ivtU)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, IVTU'
     status = nf90_put_var(fileID,varIDout(4),ivtV)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, IVTV'
  endif
  if (l_z0k) then
     status = nf90_put_var(fileID,varIDout(5),z0k)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z0K'
  endif
  if (l_smois) then
     status = nf90_put_var(fileID,varIDout(6),soilMoisture)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, SMOIS'
  endif
  ! Geopotential height
  if (l_z1000) then
     status = nf90_put_var(fileID,varIDout(14),z1000)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z1000'
  endif
  if (l_z950) then
     status = nf90_put_var(fileID,varIDout(15),z950)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z950'
  endif
    if (l_z900) then
     status = nf90_put_var(fileID,varIDout(16),z900)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z900'
  endif
  if (l_z850) then
     status = nf90_put_var(fileID,varIDout(17),z850)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z850'
  endif
  if (l_z800) then
     status = nf90_put_var(fileID,varIDout(18),z800)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z800'
  endif
  if (l_z750) then
     status = nf90_put_var(fileID,varIDout(19),z750)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z750'
  endif
  if (l_z700) then
     status = nf90_put_var(fileID,varIDout(20),z700)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z700'
  endif
  if (l_z600) then
     status = nf90_put_var(fileID,varIDout(21),z600)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z600'
  endif
  if (l_z500) then
     status = nf90_put_var(fileID,varIDout(22),z500)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z500'
  endif
  if (l_z250) then
     status = nf90_put_var(fileID,varIDout(23),z250)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z250'
  endif
  ! Specific humidity
  if (l_q1000) then
     status = nf90_put_var(fileID,varIDout(24),q1000)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q1000'
  endif
  if (l_q950) then
     status = nf90_put_var(fileID,varIDout(25),q950)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q950'
  endif
  if (l_q900) then
     status = nf90_put_var(fileID,varIDout(26),q900)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q900'
  endif
  if (l_q850) then
     status = nf90_put_var(fileID,varIDout(27),q850)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q850'
  endif
  if (l_q800) then
     status = nf90_put_var(fileID,varIDout(28),q800)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q800'
  endif
  if (l_q750) then
     status = nf90_put_var(fileID,varIDout(29),q750)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q750'
  endif
  if (l_q700) then
     status = nf90_put_var(fileID,varIDout(30),q700)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q700'
  endif
  if (l_q600) then
     status = nf90_put_var(fileID,varIDout(31),q600)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q600'
  endif
  if (l_q500) then
     status = nf90_put_var(fileID,varIDout(32),q500)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q500'
  endif
  if (l_q250) then
     status = nf90_put_var(fileID,varIDout(33),q250)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q250'
  endif
  ! Zonal wind
  if (l_u1000) then
     status = nf90_put_var(fileID,varIDout(34),u1000)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U1000'
  endif
  if (l_u950) then
     status = nf90_put_var(fileID,varIDout(35),u950)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U950'
  endif
  if (l_u900) then
     status = nf90_put_var(fileID,varIDout(36),u900)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U900'
  endif
  if (l_u850) then
     status = nf90_put_var(fileID,varIDout(37),u850)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U850'
  endif
  if (l_u800) then
     status = nf90_put_var(fileID,varIDout(38),u800)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U800'
  endif
  if (l_u750) then
     status = nf90_put_var(fileID,varIDout(39),u750)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U750'
  endif
  if (l_u700) then
     status = nf90_put_var(fileID,varIDout(40),u700)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U700'
  endif
  if (l_u600) then
     status = nf90_put_var(fileID,varIDout(41),u600)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U600'
  endif
  if (l_u500) then
     status = nf90_put_var(fileID,varIDout(42),u500)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U500'
  endif
  if (l_u250) then
     status = nf90_put_var(fileID,varIDout(43),u250)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U250'
  endif

  ! Meridional wind
  if (l_v1000) then
     status = nf90_put_var(fileID,varIDout(44),v1000)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V1000'
  endif
  if (l_v950) then
     status = nf90_put_var(fileID,varIDout(45),v950)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V950'
  endif
  if (l_v900) then
     status = nf90_put_var(fileID,varIDout(46),v900)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V900'
  endif
  if (l_v850) then
     status = nf90_put_var(fileID,varIDout(47),v850)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V850'
  endif
  if (l_v800) then
     status = nf90_put_var(fileID,varIDout(48),v800)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V800'
  endif
  if (l_v750) then
     status = nf90_put_var(fileID,varIDout(49),v750)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V750'
  endif
  if (l_v700) then
     status = nf90_put_var(fileID,varIDout(50),v700)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V700'
  endif
  if (l_v600) then
     status = nf90_put_var(fileID,varIDout(51),v600)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V600'
  endif
  if (l_v500) then
     status = nf90_put_var(fileID,varIDout(52),v500)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V500'
  endif
  if (l_v250) then
     status = nf90_put_var(fileID,varIDout(53),v250)
     if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V250'
  endif
  
  ! Close output file
  status = nf90_close(fileID)
  if (status /= nf90_NoErr) print*,'ERROR: Failure closing file'

101 continue
  ! #############################################################################
  ! END PROGRAM
  ! #############################################################################
contains
  elemental subroutine str2int(str,int,stat)
    implicit none
    ! Arguments
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int
    integer,intent(out)         :: stat
    read(str,*,iostat=stat)  int
  end subroutine str2int
  
end program WRF3D_pp
