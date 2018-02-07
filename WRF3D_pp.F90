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
    real,parameter :: &
         output_fill_value = -999.
    
    ! Namelist
    character(len=256) :: &
         fileIN, & ! Input file
         fileOUT   ! Output file
    logical :: &
         verbose, & ! Turn on/off info print to screen
         toa2sfc, & ! If vertical ordering is from TOA-2_SFC, set to .true.
         l_ivt,   & ! Compute IVT?
         l_z0k,   & ! Compute freezing level height?
         l_smois, & ! Output soil moisture?
         l_tslb,  & ! Output soil temperature?
         l_rainnc,& ! Extract accumulated grid-scale precipitations?
         l_rainc, & ! Extract accumulated grid-scale cummulus precipitation?
         l_psfc,  & ! Output surface pressure?
         l_tskin, & ! Output skin temperature?
         l_sst,   & ! Output SST?
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
         l_q2m,   & ! Output 2m specific-humidity?
         l_u1000, & ! Output 1000hPa U-component of wind?
         l_u950,  & ! Output 950hPa U-component of wind?
         l_u900,  & ! Output 900hPa U-component of wind?
         l_u850,  & ! Output 850hPa U-component of wind?
         l_u800,  & ! Output 800hPa U-component of wind?
         l_u750,  & ! Output 750hPa U-component of wind?
         l_u700,  & ! Output 700hPa U-component of wind?
         l_u600,  & ! Output 600hPa U-component of wind?
         l_u500,  & ! Output 500hPa U-component of wind?
         l_u250,  & ! Output 250hPa U-component of wind?
         l_u10m,  & ! Output 10m U-component of wind?
         l_v1000, & ! Output 1000hPa V-component of wind?
         l_v950,  & ! Output 950hPa V-component of wind?
         l_v900,  & ! Output 900hPa V-component of wind?
         l_v850,  & ! Output 850hPa V-component of wind?
         l_v800,  & ! Output 800hPa V-component of wind?
         l_v750,  & ! Output 750hPa V-component of wind?
         l_v700,  & ! Output 700hPa V-component of wind?
         l_v600,  & ! Output 600hPa V-component of wind?
         l_v500,  & ! Output 500hPa V-component of wind?
         l_v250,  & ! Output 250hPa V-component of wind?
         l_v10m,  & ! Output 10m V-component of wind?
         l_t1000, & ! Output 1000hPa temperature?
         l_t950,  & ! Output 950hPa temperature?
         l_t900,  & ! Output 900hPa temperature?
         l_t850,  & ! Output 850hPa temperature?
         l_t800,  & ! Output 800hPa temperature?
         l_t750,  & ! Output 750hPa temperature?
         l_t700,  & ! Output 700hPa temperature?
         l_t600,  & ! Output 600hPa temperature?
         l_t500,  & ! Output 500hPa temperature?
         l_t250,  & ! Output 250hPa temperature?
         l_t2m      ! Output 2m temperature?
    namelist/nmlist/fileIN, fileOUT, verbose, toa2sfc, l_ivt, l_z0k, l_smois, l_tslb, l_rainnc,   &
         l_rainc, l_psfc, l_tskin, l_sst,                                                         &
         l_z1000, l_z950, l_z900, l_z850, l_z800, l_z750, l_z700, l_z600, l_z500, l_z250,         &
         l_q1000, l_q950, l_q900, l_q850, l_q800, l_q750, l_q700, l_q600, l_q500, l_q250, l_q2m,  &
         l_u1000, l_u950, l_u900, l_u850, l_u800, l_u750, l_u700, l_u600, l_u500, l_u250, l_u10m, &
         l_v1000, l_v950, l_v900, l_v850, l_v800, l_v750, l_v700, l_v600, l_v500, l_v250, l_v10m, &
         l_t1000, l_t950, l_t900, l_t850, l_t800, l_t750, l_t700, l_t600, l_t500, l_t250, l_t2m
    
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
         smois,     & ! WRF input field: Soil mosisture                    (SMOIS)   (m3/m3)
         tslb         ! WRF input field: Soil temperature                  (TSLB)    (K)
    real,dimension(:,:,:),allocatable :: &
         lon,       & ! WRF input field: Longitude                         (XLONG)   (degree_east)
         lat,       & ! WRF input field: Latitude                          (XLAT)    (degree_north)
         lon_u,     & ! WRF input field: Longitude (staggerd u-grid)       (XLONG_U) (degree_east)
         lat_u,     & ! WRF input field: Latitude (staggerd u-grid)        (XLAT_U)  (degree_north)
         lon_v,     & ! WRF input field: Longitude (staggerd v-grid)       (XLONG_V) (degree_east)
         lat_v,     & ! WRF input field: Latitude (staggerd v-grid)        (XLAT_V)  (degree_north)
         psfc,      & ! WRF input field: Surface pressure                  (PSFC)    (Pa)
         t2m,       & ! WRF input field: Temperature @ 2m                  (T2)      (K)
         q2m,       & ! WRF input field: Specific humidity @ 2m            (Q2)      (kg/kg)
         tsk,       & ! WRF input field: Skin temperature                  (TSK)     (K)
         sst,       & ! WRF input field: SST                               (SST)     (K)
         v10m,      & ! WRF input field: v-wind @ 10m                      (U10)     (m/s)     
         u10m,      & ! WRF input field: v-wind @ 10m                      (V10)     (m/s)     
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
         rainnc,    & ! WRF total precipitaion                 (mm)
         rainc,     & ! WRF total cummulus precipitaion        (mm)
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
         t1000,     & ! WRF temperature @ 1000hPa              (K)
         t950,      & ! WRF temperature @ 950hPa               (K)
         t900,      & ! WRF temperature @ 900hPa               (K)
         t850,      & ! WRF temperature @ 850hPa               (K)
         t800,      & ! WRF temperature @ 800hPa               (K)
         t750,      & ! WRF temperature @ 750hPa               (K)
         t700,      & ! WRF temperature @ 700hPa               (K)
         t600,      & ! WRF temperature @ 600hPa               (K)
         t500,      & ! WRF temperature @ 500hPa               (K)
         t250         ! WRF temperature @ 250hPa               (K)
    integer,dimension(:,:,:) ,allocatable:: &
         i_rainnc,  & ! WRF total precipitaion bucket          (1)
         i_rainc      ! WRF total cummulus precipitaion bucket (1)
    
    integer,dimension(80) :: dimID,varIDout
    character(len=19) :: tempTime
    integer :: status,fileID,varID,ii,ij,ik,il
    real,dimension(:),allocatable :: ui,vi,temp,phm,a,b
    real :: wt1,dp,wt2,Bucket_size
    integer,dimension(1) :: xi,xf

    ! By default, turn all output to off.
    logical :: &
         lread_QVAPOR = .false., &
         lread_Q2M    = .false., &
         lread_U      = .false., &
         lread_U10    = .false., &
         lread_V10    = .false., &
         lread_V      = .false., &
         lread_PSFC   = .false., &
         lread_PB     = .false., &
         lread_P      = .false., &
         lread_PH     = .false., &
         lread_PHB    = .false., &
         lread_T      = .false., &
         lread_T00    = .false., &
         lread_T2M    = .false., &
         lread_TSK    = .false., &
         lread_HGT    = .false., &
         lread_ZS     = .false., &
         lread_SMOIS  = .false., &
         lread_TSLB   = .false., &
         lread_RAINNC = .false., &
         lread_RAINC  = .false., &
         lread_SST    = .false.
    
    ! #############################################################################
    ! A) Read in namelist
    ! #############################################################################
    open(10,file=trim(input_namelist),status='unknown')
    read(10,nml=nmlist)
    close(10)
    
    ! #############################################################################
    ! B) What WRF fields need to be read in?
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
    if (l_tslb) then
       lread_ZS     = .true.
       lread_TSLB   = .true.
    endif
    if (l_u10m) then
       lread_U10    = .true.
    endif
    if (l_v10m) then
       lread_V10    = .true.
    endif
    if (l_q2m) then
       lread_Q2M    = .true.
    endif
    if (l_t2m) then
       lread_T2M    = .true.
    endif
    if (l_rainnc) then
       Lread_RAINNC = .true.
    endif
    if (l_rainc) then
       Lread_RAINC = .true.
    endif
    if (l_tskin) then
       Lread_TSK   = .true.
    endif
    if (l_sst) then
       Lread_SST   = .true.
    endif
    if (any([l_z1000,l_z950,l_z900,l_z850,l_z800,l_z750,l_z700,l_z600,l_z500,l_z250])) then
       lread_PB     = .true.
       lread_P      = .true.
       lread_PH     = .true.
       lread_PHB    = .true.
       lread_HGT    = .true.
       lread_PSFC   = .true.
    endif
    if (any([l_q1000,l_q950,l_q900,l_q850,l_q800,l_q750,l_q700,l_q600,l_q500,l_q250])) then
       lread_QVAPOR = .true.
       lread_PB     = .true.
       lread_P      = .true.
       lread_PSFC   = .true.
       lread_HGT    = .true.
    endif
    if (any([l_u1000,l_u950,l_u900,l_u850,l_u800,l_u750,l_u700,l_u600,l_u500,l_u250])) then
       lread_U      = .true.
       lread_PSFC   = .true.
       lread_HGT    = .true.
    endif
    if (any([l_v1000,l_v950,l_v900,l_v850,l_v800,l_v750,l_v700,l_v600,l_v500,l_v250])) then
       lread_V      = .true.
       lread_PSFC   = .true.
       lread_HGT    = .true.
    endif
    if (any([l_t1000,l_t950,l_t900,l_t850,l_t800,l_t750,l_t700,l_t600,l_t500,l_t250])) then
       lread_T      = .true.
       lread_T00    = .true.
       lread_PB     = .true.
       lread_P      = .true.
       lread_PSFC   = .true.
       lread_HGT    = .true.
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
          allocate(u(nLon_stag,nLat,nLev,nTime))
          status = nf90_get_var(fileID,varID,u,count=(/nLon_stag,nLat,nLev,nTime/))
       endif
    endif
    if (lread_V) then
       ! Northward component of wind
       status = nf90_inq_varid(fileID,"V",varID)
       if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, V'
       if (status == nf90_NoErr) then
          allocate(v(nLon,nLat_stag,nLev,nTime))
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
    
    if (lread_TSLB) then
       ! Soil temperature
       status = nf90_inq_varid(fileID,"TSLB",varID)
       if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, TSLB'
       if (status == nf90_NoErr) then
          allocate(tslb(nLon,nLat,nSoil_stag,nTime))
          status = nf90_get_var(fileID,varID,tslb)
       endif
    endif
    
    if (lread_RAINNC) then
       ! Accumulated total grid scale precipitaiton.
       status = nf90_inq_varid(fileID,"RAINNC",varID)
       if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, RAINNC'
       if (status == nf90_NoErr) then
          allocate(rainnc(nLon,nLat,nTime))
          status = nf90_get_var(fileID,varID,rainnc)
       endif
       ! Bucket count
       status = nf90_inq_varid(fileID,"I_RAINNC",varID)
       if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, I_RAINNC'
       if (status == nf90_NoErr) then
          allocate(i_rainnc(nLon,nLat,nTime))
          status = nf90_get_var(fileID,varID,i_rainnc)
       endif
       ! Read in global attribute desribing "bucket size"
       status = nf90_get_att(fileID, NF90_GLOBAL, "BUCKET_MM", bucket_size)
       if (status /= nf90_NoErr) then
          print*,'ERROR: Requested global attribute not in file, BUCKET_MM, setting to default vaule of 100. mm'
          bucket_size = 100.
       endif
       ! Compute precipitation.
       where(i_rainnc .gt. 0) rainnc = rainnc+i_rainnc*bucket_size
    endif

    if (lread_RAINC) then
       ! Accumulated total grid scale cummulus precipitaiton.
       status = nf90_inq_varid(fileID,"RAINC",varID)
       if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, RAINC'
       if (status == nf90_NoErr) then
          allocate(rainc(nLon,nLat,nTime))
          status = nf90_get_var(fileID,varID,rainc)
       endif
       ! Bucket count
       status = nf90_inq_varid(fileID,"I_RAINC",varID)
       if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, I_RAINC'
       if (status == nf90_NoErr) then
          allocate(i_rainc(nLon,nLat,nTime))
          status = nf90_get_var(fileID,varID,i_rainc)
       endif
       ! Read in global attribute desribing "bucket size"
       status = nf90_get_att(fileID, NF90_GLOBAL, "BUCKET_MM", bucket_size)
       if (status /= nf90_NoErr) then
          print*,'ERROR: Requested global attribute not in file, BUCKET_MM, setting to default vaule of 100. mm'
          bucket_size = 100.
       endif
       ! Compute precipitation.
       where(i_rainc .gt. 0) rainc = rainc+i_rainc*bucket_size
    endif

    if (lread_T2M) then
       ! Temperature @ 2 meters
       status = nf90_inq_varid(fileID,"T2",varID)
       if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, T2'
       if (status == nf90_NoErr) then
          allocate(t2m(nLon,nLat,nTime))
          status = nf90_get_var(fileID,varID,t2m)
       endif
    endif
    
    if (lread_Q2M) then
       ! Temperature @ 2 meters
       status = nf90_inq_varid(fileID,"Q2",varID)
       if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, Q2'
       if (status == nf90_NoErr) then
          allocate(q2m(nLon,nLat,nTime))
          status = nf90_get_var(fileID,varID,q2m)
       endif
    endif

    if (lread_TSK) then
       ! Skin temperature
       status = nf90_inq_varid(fileID,"TSK",varID)
       if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, TSK'
       if (status == nf90_NoErr) then
          allocate(tsk(nLon,nLat,nTime))
          status = nf90_get_var(fileID,varID,tsk)
       endif
    endif

    if (lread_SST) then
       ! SST
       status = nf90_inq_varid(fileID,"SST",varID)
       if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, SST'
       if (status == nf90_NoErr) then
          allocate(sst(nLon,nLat,nTime))
          status = nf90_get_var(fileID,varID,sst)
       endif
    endif

    if (lread_U10) then
       ! Zonal wind @ 10 meters
       status = nf90_inq_varid(fileID,"U10",varID)
       if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, U10'
       if (status == nf90_NoErr) then
          allocate(u10m(nLon,nLat,nTime))
          status = nf90_get_var(fileID,varID,u10m)
       endif
    endif
    
    if (lread_V10) then
       ! Meridional wind @ 10 meters
       status = nf90_inq_varid(fileID,"V10",varID)
       if (status /= nf90_NoErr) print*,'ERROR: Requested variable not in file, V10'
       if (status == nf90_NoErr) then
          allocate(v10m(nLon,nLat,nTime))
          status = nf90_get_var(fileID,varID,v10m)
       endif
    endif
    
    ! 4) Close file
    status=nf90_close(fileID)
    if (status /= nf90_NoErr) print*,'ERROR: Could not close file, ',trim(fileiN)
    if (verbose) print*,'Finished reading in data'
    
    ! #############################################################################
    ! Part B: Computations.
    ! #############################################################################
    
    ! Compute pressure and height.
    if (allocated(pb) .and. allocated(pp)) then
       allocate(p(nLon,nLat,nLev,nTime))
       p   = pb+pp
    endif
    if (allocated(ph) .and. allocated(phb)) then
       allocate(hgt(nLon,nLat,nLev_stag,nTime))
       hgt = (ph+phb)/9.8
    endif
    
    ! Allocate space.
    allocate(ui(nLev),vi(nLev),temp(nLev),phm(nLev))
    if (l_ivt)   allocate(ivtU(nLon,  nLat, nTime))
    if (l_ivt)   allocate(ivtV(nLon,  nLat, nTime))
    if (l_z0k)   allocate(z0k( nLon,  nLat, nTime))
    if (l_z1000) then
       allocate(z1000(nlon, nLat, nTime))
       z1000(:,:,:) = output_fill_value
    endif
    if (l_z950)  then
       allocate(z950(nLon,  nLat, nTime))
       z950(:,:,:) = output_fill_value
    endif
    if (l_z900)  then
       allocate(z900(nLon,  nLat, nTime))
       z900(:,:,:) = output_fill_value
    endif
    if (l_z850)  then
       allocate(z850(nLon,  nLat, nTime))
       z850(:,:,:) = output_fill_value
    endif
    if (l_z800) then
       allocate(z800(nLon,  nLat, nTime))
       z800(:,:,:) = output_fill_value
    endif
    if (l_z750)  then
       allocate(z750(nLon,  nLat, nTime))
       z750(:,:,:) = output_fill_value
    endif
    if (l_z700)  then
       allocate(z700(nLon,  nLat, nTime))
       z700(:,:,:) = output_fill_value
    endif
    if (l_z600) then
       allocate(z600(nLon,  nLat, nTime))
       z600(:,:,:) = output_fill_value
    endif
    if (l_z500) then
       allocate(z500(nLon,  nLat, nTime))
       z500(:,:,:) = output_fill_value
    endif
    if (l_z250)  then
       allocate(z250(nLon,  nLat, nTime))
       z250(:,:,:) = output_fill_value
    endif
    if (l_q1000) then
       allocate(q1000(nlon, nLat, nTime))
       q1000(:,:,:) = output_fill_value
    endif
    if (l_q950)  then
       allocate(q950(nLon,  nLat, nTime))
       q950(:,:,:) = output_fill_value
    endif
    if (l_q900)  then
       allocate(q900(nLon,  nLat, nTime))
       q900(:,:,:) = output_fill_value
    endif
    if (l_q850)  then
       allocate(q850(nLon,  nLat, nTime))
       q850(:,:,:) = output_fill_value
    endif
    if (l_q800) then
       allocate(q800(nLon,  nLat, nTime))
       q800(:,:,:) = output_fill_value
    endif
    if (l_q750)  then
       allocate(q750(nLon,  nLat, nTime))
       q750(:,:,:) = output_fill_value
    endif
    if (l_q700)  then
       allocate(q700(nLon,  nLat, nTime))
       q700(:,:,:) = output_fill_value
    endif
    if (l_q600) then
       allocate(q600(nLon,  nLat, nTime))
       q600(:,:,:) = output_fill_value
    endif
    if (l_q500) then
       allocate(q500(nLon,  nLat, nTime))
       q500(:,:,:) = output_fill_value
    endif
    if (l_q250)  then
       allocate(q250(nLon,  nLat, nTime))
       q250(:,:,:) = output_fill_value
    endif
    if (l_u1000) then
       allocate(u1000(nlon, nLat, nTime))
       u1000(:,:,:) = output_fill_value
    endif
    if (l_u950)  then
       allocate(u950(nLon,  nLat, nTime))
       u950(:,:,:) = output_fill_value
    endif
    if (l_u900)  then
       allocate(u900(nLon,  nLat, nTime))
       u900(:,:,:) = output_fill_value
    endif
    if (l_u850)  then
       allocate(u850(nLon,  nLat, nTime))
       u850(:,:,:) = output_fill_value
    endif
    if (l_u800) then
       allocate(u800(nLon,  nLat, nTime))
       u800(:,:,:) = output_fill_value
    endif
    if (l_u750)  then
       allocate(u750(nLon,  nLat, nTime))
       u750(:,:,:) = output_fill_value
    endif
    if (l_u700)  then
       allocate(u700(nLon,  nLat, nTime))
       u700(:,:,:) = output_fill_value
    endif
    if (l_u600) then
       allocate(u600(nLon,  nLat, nTime))
       u600(:,:,:) = output_fill_value
    endif
    if (l_u500) then
       allocate(u500(nLon,  nLat, nTime))
       u500(:,:,:) = output_fill_value
    endif
    if (l_u250)  then
       allocate(u250(nLon,  nLat, nTime))
       u250(:,:,:) = output_fill_value
    endif
    if (l_v1000) then
       allocate(v1000(nlon, nLat, nTime))
       v1000(:,:,:) = output_fill_value
    endif
    if (l_v950)  then
       allocate(v950(nLon,  nLat, nTime))
       v950(:,:,:) = output_fill_value
    endif
    if (l_v900)  then
       allocate(v900(nLon,  nLat, nTime))
       v900(:,:,:) = output_fill_value
    endif
    if (l_v850)  then
       allocate(v850(nLon,  nLat, nTime))
       v850(:,:,:) = output_fill_value
    endif
    if (l_v800) then
       allocate(v800(nLon,  nLat, nTime))
       v800(:,:,:) = output_fill_value
    endif
    if (l_v750)  then
       allocate(v750(nLon,  nLat, nTime))
       v750(:,:,:) = output_fill_value
    endif
    if (l_v700)  then
       allocate(v700(nLon,  nLat, nTime))
       v700(:,:,:) = output_fill_value
    endif
    if (l_v600) then
       allocate(v600(nLon,  nLat, nTime))
       v600(:,:,:) = output_fill_value
    endif
    if (l_v500) then
       allocate(v500(nLon,  nLat, nTime))
       v500(:,:,:) = output_fill_value
    endif
    if (l_v250)  then
       allocate(v250(nLon,  nLat, nTime))
       v250(:,:,:) = output_fill_value
    endif
    if (l_t1000) then
       allocate(t1000(nlon, nLat, nTime))
       t1000(:,:,:) = output_fill_value
    endif
    if (l_t950)  then
       allocate(t950(nLon,  nLat, nTime))
       t950(:,:,:) = output_fill_value
    endif
    if (l_t900)  then
       allocate(t900(nLon,  nLat, nTime))
       t900(:,:,:) = output_fill_value
    endif
    if (l_t850)  then
       allocate(t850(nLon,  nLat, nTime))
       t850(:,:,:) = output_fill_value
    endif
    if (l_t800) then
       allocate(t800(nLon,  nLat, nTime))
       t800(:,:,:) = output_fill_value
    endif
    if (l_t750)  then
       allocate(t750(nLon,  nLat, nTime))
       t750(:,:,:) = output_fill_value
    endif
    if (l_t700)  then
       allocate(t700(nLon,  nLat, nTime))
       t700(:,:,:) = output_fill_value
    endif
    if (l_t600) then
       allocate(t600(nLon,  nLat, nTime))
       t600(:,:,:) = output_fill_value
    endif
    if (l_t500) then
       allocate(t500(nLon,  nLat, nTime))
       t500(:,:,:) = output_fill_value
    endif
    if (l_t250)  then
       allocate(t250(nLon,  nLat, nTime))
       t250(:,:,:) = output_fill_value
    endif
    
    ! Loop over all points/times.
    if (verbose) print*,'Begin computations... '
    do ii=1,nTime
       if (verbose) write(*,"(a12,i2,a4,i2)"),'   @Timestep ',ii,' of ',nTime
       do ij=1,nLon
          do ik=1,nLat
             ! If needed, put winds on mass centered grid points, in WRF the velocity components are on staggered grids.
             ! Eastward component of wind (on mass-centered point)
             if (any([l_u1000,l_u950,l_u900,l_u850,l_u800,l_u750,l_u700,l_u600,l_u500,l_u250]) .or. l_ivt) then
                wt1 = (lon_u(ij+1,ik,ii)-lon(ij,ik,ii))/(lon_u(ij+1,ik,ii)-lon_u(ij,ik,ii))
                ui = wt1*u(ij,ik,:,ii)+(1-wt1)*u(ij+1,ik,:,ii)
             endif
             ! Northward component of wind (on mass-centered point)
             if (any([l_v1000,l_v950,l_v900,l_v850,l_v800,l_v750,l_v700,l_v600,l_v500,l_v250]) .or. l_ivt) then
                wt1 = (lat_v(ij,ik+1,ii)-lat(ij,ik,ii)) / (lat_v(ij,ik+1,ii)-lat_v(ij,ik,ii))
                vi = wt1*v(ij,ik,:,ii)+(1-wt1)*v(ij,ik+1,:,ii)
             endif

             ! If needed, compute temperature from pertubation potential temperature.
             if (any([l_t1000,l_t950,l_t900,l_t850,l_t800,l_t750,l_t700,l_t600,l_t500,l_t250]) .or. l_z0k) then
                temp = (ta(ij,ik,:,ii)+t00(ii))*(101325.0/p(ij,ik,:,ii))**(-0.286)
             endif

             ! Geopotential heights are on a staggerd vertical grid, so compute mass-
             ! centered geopotential height.
             if (any([l_u1000,l_u950,l_u900,l_u850,l_u800,l_u750,l_u700,l_u600,l_u500,l_u250]) .or.  &
                 any([l_v1000,l_v950,l_v900,l_v850,l_v800,l_v750,l_v700,l_v600,l_v500,l_v250]) .or.  &
                 any([l_t1000,l_t950,l_t900,l_t850,l_t800,l_t750,l_t700,l_t600,l_t500,l_t250]) .or.  &
                 any([l_q1000,l_q950,l_q900,l_q850,l_q800,l_q750,l_q700,l_q600,l_q500,l_q250]) .or.  &
                 any([l_z1000,l_z950,l_z900,l_z850,l_z800,l_z750,l_z700,l_z600,l_z500,l_z250]) .or.  &
                 l_z0k) then
                phm  = (hgt(ij,ik,1:nlev_stag-1,ii)+hgt(ij,ik,2:nlev_stag,ii))*0.5
             endif
             
             ! ######################################################################
             ! Compute integrated vapor transport (IVT)
             ! ######################################################################
             if (l_ivt) then
                do il=1,nLev
                   ! Compute pressure change across layer.
                   if (.not. toa2sfc) then
                      if (il .eq. 1)                    dp = psfc(ij,ik,ii)-sum(p(ij,ik,il:il+1,ii))/2.              ! Bottom level
                      if (il .ne. 1 .and. il .ne. nLev) dp = sum(p(ij,ik,il-1:il,ii))/2.-sum(p(ij,ik,il:il+1,ii))/2. ! Middle levels
                      if (il .eq. nLev)                 dp = 0.                                                      ! Top level
                   else
                      if (il .eq. 1)                    dp = 0.                                                      ! Top level
                      if (il .ne. 1 .and. il .ne. nLev) dp = sum(p(ij,ik,il-1:il,ii))/2.-sum(p(ij,ik,il:il+1,ii))/2. ! Middle levels
                      if (il .eq. nLev)                 dp = psfc(ij,ik,ii)-sum(p(ij,ik,il:il+1,ii))/2.              ! Bottom level
                   endif
                   ! Integrate vertically.
                   ivtU(ij,ik,ii) = ivtU(ij,ik,ii) + ui(il)*q(ij,ik,il,ii)*dp/9.8
                   ivtV(ij,ik,ii) = ivtV(ij,ik,ii) + vi(il)*q(ij,ik,il,ii)*dp/9.8
                enddo
             endif

             ! ######################################################################
             ! Compute freezing-level height.
             ! ######################################################################
             if (l_z0k) then
                ! Only compute if freezing level is present.
                if (any(temp .gt. 273.16)) then
                   ! What is the highest level that is above freezing?
                   xi = minloc(temp-273.15,temp .gt. 273.16)
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
             ! Compute fields at synoptic levels.
             ! ###################################################################
             ! 1000hPa
             if (any([l_u1000, l_v1000, l_q1000, l_t1000, l_z1000])) then
                if (psfc(ij,ik,ii) .ge. 100000.) then
                   call compute_interpolation_wts(p(ij,ik,:,ii),100000.,xi,xf,wt2)
                   if (l_u1000) u1000(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
                   if (l_v1000) v1000(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
                   if (l_q1000) q1000(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
                   if (l_t1000) t1000(ij,ik,ii) = temp(xi(1))*wt2 + temp(xf(1))*(1-wt2)
                   if (l_z1000) z1000(ij,ik,ii) = phm(xi(1))*wt2 + phm(xf(1))*(1-wt2)
                endif
             endif
             ! 950hPa
             if (any([l_u950, l_v950, l_q950, l_t950, l_z950])) then
                if (psfc(ij,ik,ii) .ge. 95000.) then
                   call compute_interpolation_wts(p(ij,ik,:,ii),95000.,xi,xf,wt2)
                   if (l_u950) u950(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
                   if (l_v950) v950(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
                   if (l_q950) q950(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
                   if (l_t950) t950(ij,ik,ii) = temp(xi(1))*wt2 + temp(xf(1))*(1-wt2)
                   if (l_z950) z950(ij,ik,ii) = phm(xi(1))*wt2 + phm(xf(1))*(1-wt2)
                endif
             endif
             ! 900hPa
             if (any([l_u900, l_v900, l_q900, l_t900, l_z900])) then
                if (psfc(ij,ik,ii) .ge. 90000.) then
                   call compute_interpolation_wts(p(ij,ik,:,ii),90000.,xi,xf,wt2)
                   if (l_u900) u900(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
                   if (l_v900) v900(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
                   if (l_q900) q900(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
                   if (l_t900) t900(ij,ik,ii) = temp(xi(1))*wt2 + temp(xf(1))*(1-wt2)
                   if (l_z900) z900(ij,ik,ii) = phm(xi(1))*wt2 + phm(xf(1))*(1-wt2)
                endif
             endif
             ! 850hPa
             if (any([l_u850, l_v850, l_q850, l_t850, l_z850])) then
                if (psfc(ij,ik,ii) .ge. 85000.) then
                   call compute_interpolation_wts(p(ij,ik,:,ii),85000.,xi,xf,wt2)
                   if (l_u850) u850(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
                   if (l_v850) v850(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
                   if (l_q850) q850(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
                   if (l_t850) t850(ij,ik,ii) = temp(xi(1))*wt2 + temp(xf(1))*(1-wt2)
                   if (l_z850) z850(ij,ik,ii) = phm(xi(1))*wt2 + phm(xf(1))*(1-wt2)
                endif
             endif
             ! 800hPa
             if (any([l_u800, l_v800, l_q800, l_t800, l_z800])) then
                if (psfc(ij,ik,ii) .ge. 80000.) then
                   call compute_interpolation_wts(p(ij,ik,:,ii),80000.,xi,xf,wt2)
                   if (l_u800) u800(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
                   if (l_v800) v800(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
                   if (l_q800) q800(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
                   if (l_t800) t800(ij,ik,ii) = temp(xi(1))*wt2 + temp(xf(1))*(1-wt2)
                   if (l_z800) z800(ij,ik,ii) = phm(xi(1))*wt2 + phm(xf(1))*(1-wt2)
                endif
             endif
             ! 750hPa
             if (any([l_u750, l_v750, l_q750, l_t750, l_z750])) then
                if (psfc(ij,ik,ii) .ge. 75000.) then
                   call compute_interpolation_wts(p(ij,ik,:,ii),75000.,xi,xf,wt2)
                   if (l_u750) u750(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
                   if (l_v750) v750(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
                   if (l_q750) q750(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
                   if (l_t750) t750(ij,ik,ii) = temp(xi(1))*wt2 + temp(xf(1))*(1-wt2)
                   if (l_z750) z750(ij,ik,ii) = phm(xi(1))*wt2 + phm(xf(1))*(1-wt2)
                endif
             endif
             ! 700hPa
             if (any([l_u700, l_v700, l_q700, l_t700, l_z700])) then
                if (psfc(ij,ik,ii) .ge. 70000.) then
                   call compute_interpolation_wts(p(ij,ik,:,ii),70000.,xi,xf,wt2)
                   if (l_u700) u700(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
                   if (l_v700) v700(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
                   if (l_q700) q700(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
                   if (l_t700) t700(ij,ik,ii) = temp(xi(1))*wt2 + temp(xf(1))*(1-wt2)
                   if (l_z700) z700(ij,ik,ii) = phm(xi(1))*wt2 + phm(xf(1))*(1-wt2)
                endif
             endif
             ! 600hPa
             if (any([l_u600, l_v600, l_q600, l_t600, l_z600])) then
                if (psfc(ij,ik,ii) .ge. 60000.) then
                   call compute_interpolation_wts(p(ij,ik,:,ii),60000.,xi,xf,wt2)
                   if (l_u600) u600(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
                   if (l_v600) v600(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
                   if (l_q600) q600(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
                   if (l_t600) t600(ij,ik,ii) = temp(xi(1))*wt2 + temp(xf(1))*(1-wt2)
                   if (l_z600) z600(ij,ik,ii) = phm(xi(1))*wt2 + phm(xf(1))*(1-wt2)
                endif
             endif
             ! 500hPa
             if (any([l_u500, l_v500, l_q500, l_t500, l_z500])) then
                if (psfc(ij,ik,ii) .ge. 50000.) then
                   call compute_interpolation_wts(p(ij,ik,:,ii),50000.,xi,xf,wt2)
                   if (l_u500) u500(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
                   if (l_v500) v500(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
                   if (l_q500) q500(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
                   if (l_t500) t500(ij,ik,ii) = temp(xi(1))*wt2 + temp(xf(1))*(1-wt2)
                   if (l_z500) z500(ij,ik,ii) = phm(xi(1))*wt2 + phm(xf(1))*(1-wt2)
                endif
             endif
             ! 250hPa
             if (any([l_u250, l_v250, l_q250, l_t250, l_z250])) then
                if (psfc(ij,ik,ii) .ge. 25000.) then
                   call compute_interpolation_wts(p(ij,ik,:,ii),25000.,xi,xf,wt2)
                   if (l_u250) u250(ij,ik,ii) = ui(xi(1))*wt2 + ui(xf(1))*(1-wt2)
                   if (l_v250) v250(ij,ik,ii) = vi(xi(1))*wt2 + vi(xf(1))*(1-wt2)
                   if (l_q250) q250(ij,ik,ii) = q(ij,ik,xi(1),ii)*wt2 + q(ij,ik,xf(1),ii)*(1-wt2)
                   if (l_t250) t250(ij,ik,ii) = temp(xi(1))*wt2 + temp(xf(1))*(1-wt2)
                   if (l_z250) z250(ij,ik,ii) = phm(xi(1))*wt2 + phm(xf(1))*(1-wt2)
                endif
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
    status = nf90_def_dim(fileID,"soil_layers_stag",nSoil_stag,dimID(5))
    if (status /= nf90_NoErr) print*,'ERROR: Failure making dimension ID, soil_layers_stag'
    
    ! Define variables and add attributes.
    ! 1D
    call init_ncdfOutVar1Dint(fileID,varIDout(10),dimID(1),"Year"," "," ")
    call init_ncdfOutVar1Dint(fileID,varIDout(11),dimID(1),"Month"," "," ")
    call init_ncdfOutVar1Dint(fileID,varIDout(12),dimID(1),"Day"," "," ")
    call init_ncdfOutVar1Dint(fileID,varIDout(13),dimID(1),"Hour"," "," ")
    if (l_smois .or. l_tslb) then
       call init_ncdfOutVar1Dflt(fileID,varIDout(7),dimID(5),"ZS","DEPTHS OF CENTERS OF SOIL LAYERS","m")
    endif
    
    ! 2D
    call init_ncdfOutVar2D(fileID,varIDout(1),dimID(2),dimID(3),"XLAT", "LATITUDE, SOUTH IS NEGATIVE","degree_north")
    call init_ncdfOutVar2D(fileID,varIDout(2),dimID(2),dimID(3),"XLONG","LONGITUDE, WEST IS NEGATIVE","degree_east")
    
    ! 3D
    if (l_ivt) then
       call init_ncdfOutVar3D(fileID,varIDout(3),dimID(2),dimID(3),dimID(1),"IVTU","ZONAL COMPONENT OF THE INTEGRATED VAPOR TRANSPORT","kg m-2 s-1")
       call init_ncdfOutVar3D(fileID,varIDout(4),dimID(2),dimID(3),dimID(1),"IVTV","MERIDIONAL COMPONENT OF THE INTEGRATED VAPOR TRANSPORT","kg m-2 s-1")
    endif
    if (l_z0k)     call init_ncdfOutVar3D(fileID,varIDout(5), dimID(2),dimID(3),dimID(1),"Z0K",   "FREEZING LEVEL HEIGHT",         "m")
    if (l_rainnc)  call init_ncdfOutVar3D(fileID,varIDout(68),dimID(2),dimID(3),dimID(1),"RAINNC","TOTAL GRID SCALE PRECIPITATION","mm")
    if (l_rainc)   call init_ncdfOutVar3D(fileID,varIDout(69),dimID(2),dimID(3),dimID(1),"RAINC", "TOTAL CUMMULUS PRECIPITATION",  "mm")
    if (l_psfc)    call init_ncdfOutVar3D(fileID,varIDout(70),dimID(2),dimID(3),dimID(1),"PSFC",  "SFC PRESSURE",                  "Pa")
    if (l_tskin)   call init_ncdfOutVar3D(fileID,varIDout(71),dimID(2),dimID(3),dimID(1),"TSK",   "SKIN TEMPERATURE",              "K")
    if (l_sst)     call init_ncdfOutVar3D(fileID,varIDout(72),dimID(2),dimID(3),dimID(1),"SST",   "SEA SURFACE TEMPERATURE",       "K")
    if (l_z1000)   call init_ncdfOutVar3D(fileID,varIDout(14),dimID(2),dimID(3),dimID(1),"Z1000", "GEOPOTENTIAL HEIGHT @ 1000hPa", "m")
    if (l_z950)    call init_ncdfOutVar3D(fileID,varIDout(15),dimID(2),dimID(3),dimID(1),"Z950",  "GEOPOTENTIAL HEIGHT @ 950hPa","  m")
    if (l_z900)    call init_ncdfOutVar3D(fileID,varIDout(16),dimID(2),dimID(3),dimID(1),"Z900",  "GEOPOTENTIAL HEIGHT @ 900hPa",  "m")
    if (l_z850)    call init_ncdfOutVar3D(fileID,varIDout(17),dimID(2),dimID(3),dimID(1),"Z850",  "GEOPOTENTIAL HEIGHT @ 850hPa",  "m")
    if (l_z800)    call init_ncdfOutVar3D(fileID,varIDout(18),dimID(2),dimID(3),dimID(1),"Z800",  "GEOPOTENTIAL HEIGHT @ 800hPa",  "m")
    if (l_z750)    call init_ncdfOutVar3D(fileID,varIDout(19),dimID(2),dimID(3),dimID(1),"Z750",  "GEOPOTENTIAL HEIGHT @ 750hPa",  "m")
    if (l_z700)    call init_ncdfOutVar3D(fileID,varIDout(20),dimID(2),dimID(3),dimID(1),"Z700",  "GEOPOTENTIAL HEIGHT @ 700hPa",  "m")
    if (l_z600)    call init_ncdfOutVar3D(fileID,varIDout(21),dimID(2),dimID(3),dimID(1),"Z600",  "GEOPOTENTIAL HEIGHT @ 600hPa",  "m")
    if (l_z500)    call init_ncdfOutVar3D(fileID,varIDout(22),dimID(2),dimID(3),dimID(1),"Z500",  "GEOPOTENTIAL HEIGHT @ 500hPa",  "m")
    if (l_z250)    call init_ncdfOutVar3D(fileID,varIDout(23),dimID(2),dimID(3),dimID(1),"Z250",  "GEOPOTENTIAL HEIGHT @ 250hPa",  "m")
    if (l_q1000)   call init_ncdfOutVar3D(fileID,varIDout(24),dimID(2),dimID(3),dimID(1),"Q1000", "SPECIFIC HUMIDITY @ 1000hPa",   "kg/kg")
    if (l_q950)    call init_ncdfOutVar3D(fileID,varIDout(25),dimID(2),dimID(3),dimID(1),"Q950",  "SPECIFIC HUMIDITY @ 950hPa",    "kg/kg")
    if (l_q900)    call init_ncdfOutVar3D(fileID,varIDout(26),dimID(2),dimID(3),dimID(1),"Q900",  "SPECIFIC HUMIDITY @ 900hPa",    "kg/kg")
    if (l_q850)    call init_ncdfOutVar3D(fileID,varIDout(27),dimID(2),dimID(3),dimID(1),"Q850",  "SPECIFIC HUMIDITY @ 850hPa",    "kg/kg")
    if (l_q800)    call init_ncdfOutVar3D(fileID,varIDout(28),dimID(2),dimID(3),dimID(1),"Q800",  "SPECIFIC HUMIDITY @ 800hPa",    "kg/kg")
    if (l_q750)    call init_ncdfOutVar3D(fileID,varIDout(29),dimID(2),dimID(3),dimID(1),"Q750",  "SPECIFIC HUMIDITY @ 750hPa",    "kg/kg")
    if (l_q700)    call init_ncdfOutVar3D(fileID,varIDout(30),dimID(2),dimID(3),dimID(1),"Q700",  "SPECIFIC HUMIDITY @ 700hPa",    "kg/kg")
    if (l_q600)    call init_ncdfOutVar3D(fileID,varIDout(31),dimID(2),dimID(3),dimID(1),"Q600",  "SPECIFIC HUMIDITY @ 600hPa",    "kg/kg")
    if (l_q500)    call init_ncdfOutVar3D(fileID,varIDout(32),dimID(2),dimID(3),dimID(1),"Q500",  "SPECIFIC HUMIDITY @ 500hPa",    "kg/kg")
    if (l_q250)    call init_ncdfOutVar3D(fileID,varIDout(33),dimID(2),dimID(3),dimID(1),"Q250",  "SPECIFIC HUMIDITY @ 250hPa",    "kg/kg")
    if (l_q2m)     call init_ncdfOutVar3D(fileID,varIDout(67),dimID(2),dimID(3),dimID(1),"Q2m",   "SPECIFIC HUMIDITY @ 2m",        "kg/kg")
    if (l_u1000)   call init_ncdfOutVar3D(fileID,varIDout(34),dimID(2),dimID(3),dimID(1),"U1000", "ZONAL WIND @ 1000hPa",          "m/s")
    if (l_u950)    call init_ncdfOutVar3D(fileID,varIDout(35),dimID(2),dimID(3),dimID(1),"U950",  "ZONAL WIND @ 950hPa",           "m/s")
    if (l_u900)    call init_ncdfOutVar3D(fileID,varIDout(36),dimID(2),dimID(3),dimID(1),"U900",  "ZONAL WIND @ 900hPa",           "m/s")
    if (l_u850)    call init_ncdfOutVar3D(fileID,varIDout(37),dimID(2),dimID(3),dimID(1),"U850",  "ZONAL WIND @ 850hPa",           "m/s")
    if (l_u800)    call init_ncdfOutVar3D(fileID,varIDout(38),dimID(2),dimID(3),dimID(1),"U800",  "ZONAL WIND @ 800hPa",           "m/s")
    if (l_u750)    call init_ncdfOutVar3D(fileID,varIDout(39),dimID(2),dimID(3),dimID(1),"U750",  "ZONAL WIND @ 750hPa",           "m/s")
    if (l_u700)    call init_ncdfOutVar3D(fileID,varIDout(40),dimID(2),dimID(3),dimID(1),"U700",  "ZONAL WIND @ 700hPa",           "m/s")
    if (l_u600)    call init_ncdfOutVar3D(fileID,varIDout(41),dimID(2),dimID(3),dimID(1),"U600",  "ZONAL WIND @ 600hPa",           "m/s")
    if (l_u500)    call init_ncdfOutVar3D(fileID,varIDout(42),dimID(2),dimID(3),dimID(1),"U500",  "ZONAL WIND @ 500hPa",           "m/s")
    if (l_u250)    call init_ncdfOutVar3D(fileID,varIDout(43),dimID(2),dimID(3),dimID(1),"U250",  "ZONAL WIND @ 250hPa",           "m/s")
    if (l_u10m)    call init_ncdfOutVar3D(fileID,varIDout(54),dimID(2),dimID(3),dimID(1),"U10m",  "ZONAL WIND @ 10m",              "m/s")
    if (l_v1000)   call init_ncdfOutVar3D(fileID,varIDout(44),dimID(2),dimID(3),dimID(1),"V1000", "MERIDIONAL WIND @ 1000hPa",     "m/s")
    if (l_v950)    call init_ncdfOutVar3D(fileID,varIDout(45),dimID(2),dimID(3),dimID(1),"V950",  "MERIDIONAL WIND @ 950hPa",      "m/s")
    if (l_v900)    call init_ncdfOutVar3D(fileID,varIDout(46),dimID(2),dimID(3),dimID(1),"V900",  "MERIDIONAL WIND @ 900hPa",      "m/s")
    if (l_v850)    call init_ncdfOutVar3D(fileID,varIDout(47),dimID(2),dimID(3),dimID(1),"V850",  "MERIDIONAL WIND @ 850hPa",      "m/s")
    if (l_v800)    call init_ncdfOutVar3D(fileID,varIDout(48),dimID(2),dimID(3),dimID(1),"V800",  "MERIDIONAL WIND @ 800hPa",      "m/s")
    if (l_v750)    call init_ncdfOutVar3D(fileID,varIDout(49),dimID(2),dimID(3),dimID(1),"V750",  "MERIDIONAL WIND @ 750hPa",      "m/s")
    if (l_v700)    call init_ncdfOutVar3D(fileID,varIDout(50),dimID(2),dimID(3),dimID(1),"V700",  "MERIDIONAL WIND @ 700hPa",      "m/s")
    if (l_v600)    call init_ncdfOutVar3D(fileID,varIDout(51),dimID(2),dimID(3),dimID(1),"V600",  "MERIDIONAL WIND @ 600hPa",      "m/s")
    if (l_v500)    call init_ncdfOutVar3D(fileID,varIDout(52),dimID(2),dimID(3),dimID(1),"V500",  "MERIDIONAL WIND @ 500hPa",      "m/s")
    if (l_v250)    call init_ncdfOutVar3D(fileID,varIDout(53),dimID(2),dimID(3),dimID(1),"V250",  "MERIDIONAL WIND @ 250hPa",      "m/s")
    if (l_v10m)    call init_ncdfOutVar3D(fileID,varIDout(65),dimID(2),dimID(3),dimID(1),"V10m",  "MERIDIONAL WIND @ 10m",         "m/s")
    if (l_t1000)   call init_ncdfOutVar3D(fileID,varIDout(64),dimID(2),dimID(3),dimID(1),"T1000", "TEMPERATURE @ 1000hPa",         "K")
    if (l_t950)    call init_ncdfOutVar3D(fileID,varIDout(55),dimID(2),dimID(3),dimID(1),"T950",  "TEMPERATURE @ 950hPa",          "K")
    if (l_t900)    call init_ncdfOutVar3D(fileID,varIDout(56),dimID(2),dimID(3),dimID(1),"T900",  "TEMPERATURE @ 900hPa",          "K")
    if (l_t850)    call init_ncdfOutVar3D(fileID,varIDout(57),dimID(2),dimID(3),dimID(1),"T850",  "TEMPERATURE @ 850hPa",          "K")
    if (l_t800)    call init_ncdfOutVar3D(fileID,varIDout(58),dimID(2),dimID(3),dimID(1),"T800",  "TEMPERATURE @ 800hPa",          "K")
    if (l_t750)    call init_ncdfOutVar3D(fileID,varIDout(59),dimID(2),dimID(3),dimID(1),"T750",  "TEMPERATURE @ 750hPa",          "K")
    if (l_t700)    call init_ncdfOutVar3D(fileID,varIDout(60),dimID(2),dimID(3),dimID(1),"T700",  "TEMPERATURE @ 700hPa",          "K")
    if (l_t600)    call init_ncdfOutVar3D(fileID,varIDout(61),dimID(2),dimID(3),dimID(1),"T600",  "TEMPERATURE @ 600hPa",          "K")
    if (l_t500)    call init_ncdfOutVar3D(fileID,varIDout(62),dimID(2),dimID(3),dimID(1),"T500",  "TEMPERATURE @ 500hPa",          "K")
    if (l_t250)    call init_ncdfOutVar3D(fileID,varIDout(63),dimID(2),dimID(3),dimID(1),"T250",  "TEMPERATURE @ 250hPa",          "K")
    if (l_t2m)     call init_ncdfOutVar3D(fileID,varIDout(66),dimID(2),dimID(3),dimID(1),"T2m",   "TEMPERATURE @ 2m",              "K")
    
    ! 4D Outputs
    if (l_smois)   call init_ncdfOutVar4D(fileID,varIDout(6), dimID(2),dimID(3),dimID(5),dimID(1),"SMOIS", "SOIL MOISTURE","m3/m3")
    if (l_tslb)    call init_ncdfOutVar4D(fileID,varIDout(8), dimID(2),dimID(3),dimID(5),dimID(1),"TSLB",  "SOIL TEMPERATURE","K")
    
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
    
    if (l_smois .or. l_tslb) status = nf90_put_var(fileID,varIDout(7),zs)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, ZS'
    if (l_ivt) then
       status = nf90_put_var(fileID,varIDout(3),ivtU)
       if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, IVTU'
       status = nf90_put_var(fileID,varIDout(4),ivtV)
       if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, IVTV'
    endif
    if (l_z0k) status = nf90_put_var(fileID,varIDout(5),z0k)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z0K'
    if (l_smois) status = nf90_put_var(fileID,varIDout(6),smois)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, SMOIS'
    if (l_tslb) status = nf90_put_var(fileID,varIDout(8),tslb)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, TSLB'
    if (l_rainnc) status = nf90_put_var(fileID,varIDout(68),rainnc)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, RAINNC'
    if (l_rainc) status = nf90_put_var(fileID,varIDout(69),rainc)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, RAINC'
    if (l_t2m) status = nf90_put_var(fileID,varIDout(70),psfc)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, PSFC'
    if (l_t2m) status = nf90_put_var(fileID,varIDout(71),tsk)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, TSK'
    if (l_t2m) status = nf90_put_var(fileID,varIDout(72),sst)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, SST'
    if (l_z1000) status = nf90_put_var(fileID,varIDout(14),z1000)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z1000'
    if (l_z950) status = nf90_put_var(fileID,varIDout(15),z950)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z950'
    if (l_z900) status = nf90_put_var(fileID,varIDout(16),z900)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z900'
    if (l_z850) status = nf90_put_var(fileID,varIDout(17),z850)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z850'
    if (l_z800) status = nf90_put_var(fileID,varIDout(18),z800)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z800'
    if (l_z750) status = nf90_put_var(fileID,varIDout(19),z750)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z750'
    if (l_z700)status = nf90_put_var(fileID,varIDout(20),z700)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z700'
    if (l_z600) status = nf90_put_var(fileID,varIDout(21),z600)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z600'
    if (l_z500) status = nf90_put_var(fileID,varIDout(22),z500)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z500'
    if (l_z250) status = nf90_put_var(fileID,varIDout(23),z250)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Z250'
    if (l_q1000) status = nf90_put_var(fileID,varIDout(24),q1000)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q1000'
    if (l_q950) status = nf90_put_var(fileID,varIDout(25),q950)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q950'
    if (l_q900) status = nf90_put_var(fileID,varIDout(26),q900)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q900'
    if (l_q850) status = nf90_put_var(fileID,varIDout(27),q850)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q850'
    if (l_q800)  status = nf90_put_var(fileID,varIDout(28),q800)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q800'
    if (l_q750) status = nf90_put_var(fileID,varIDout(29),q750)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q750'
    if (l_q700) status = nf90_put_var(fileID,varIDout(30),q700)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q700'
    if (l_q600) status = nf90_put_var(fileID,varIDout(31),q600)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q600'
    if (l_q500) status = nf90_put_var(fileID,varIDout(32),q500)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q500'
    if (l_q250) status = nf90_put_var(fileID,varIDout(33),q250)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q250'
    if (l_q2m)  status = nf90_put_var(fileID,varIDout(67),q2m)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, Q2m'    
    if (l_u1000) status = nf90_put_var(fileID,varIDout(34),u1000)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U1000'
    if (l_u950) status = nf90_put_var(fileID,varIDout(35),u950)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U950'
    if (l_u900) status = nf90_put_var(fileID,varIDout(36),u900)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U900'
    if (l_u850) status = nf90_put_var(fileID,varIDout(37),u850)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U850'
    if (l_u800) status = nf90_put_var(fileID,varIDout(38),u800)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U800'
    if (l_u750) status = nf90_put_var(fileID,varIDout(39),u750)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U750'
    if (l_u700) status = nf90_put_var(fileID,varIDout(40),u700)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U700'
    if (l_u600) status = nf90_put_var(fileID,varIDout(41),u600)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U600'
    if (l_u500) status = nf90_put_var(fileID,varIDout(42),u500)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U500'
    if (l_u250) status = nf90_put_var(fileID,varIDout(43),u250)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U250'
    if (l_u10m) status = nf90_put_var(fileID,varIDout(54),u10m)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, U10m'
    if (l_v1000) status = nf90_put_var(fileID,varIDout(44),v1000)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V1000'
    if (l_v950)  status = nf90_put_var(fileID,varIDout(45),v950)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V950'
    if (l_v900)  status = nf90_put_var(fileID,varIDout(46),v900)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V900'
    if (l_v850)  status = nf90_put_var(fileID,varIDout(47),v850)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V850'
    if (l_v800) status = nf90_put_var(fileID,varIDout(48),v800)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V800'
    if (l_v750) status = nf90_put_var(fileID,varIDout(49),v750)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V750'
    if (l_v700) status = nf90_put_var(fileID,varIDout(50),v700)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V700'
    if (l_v600) status = nf90_put_var(fileID,varIDout(51),v600)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V600'
    if (l_v500) status = nf90_put_var(fileID,varIDout(52),v500)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V500'
    if (l_v250) status = nf90_put_var(fileID,varIDout(53),v250)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V250'
    if (l_v10m)  status = nf90_put_var(fileID,varIDout(65),v10m)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, V10m'    
    if (l_t1000) status = nf90_put_var(fileID,varIDout(64),t1000)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, T1000'
    if (l_t950) status = nf90_put_var(fileID,varIDout(55),t950)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, T950'
    if (l_t900) status = nf90_put_var(fileID,varIDout(56),t900)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, T900'
    if (l_t850) status = nf90_put_var(fileID,varIDout(57),t850)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, T850'
    if (l_t800) status = nf90_put_var(fileID,varIDout(58),t800)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, T800'
    if (l_t750) status = nf90_put_var(fileID,varIDout(59),t750)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, T750'
    if (l_t700) status = nf90_put_var(fileID,varIDout(60),t700)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, T700'
    if (l_t600) status = nf90_put_var(fileID,varIDout(61),t600)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, T600'
    if (l_t500) status = nf90_put_var(fileID,varIDout(62),t500)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, T500'
    if (l_t250) status = nf90_put_var(fileID,varIDout(63),t250)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, T250'
    if (l_t2m) status = nf90_put_var(fileID,varIDout(66),t2m)
    if (status /= nf90_NoErr) print*,'ERROR: Failure populating output field, T2m'
    
    ! Close output file
    status = nf90_close(fileID)
    if (status /= nf90_NoErr) print*,'ERROR: Failure closing file'
    
101 continue
    ! #############################################################################
    ! END PROGRAM
    ! #############################################################################
  contains

    ! #############################################################################
    ! Subroutine to convert a character into integer.
    ! #############################################################################
    elemental subroutine str2int(str,int,stat)
      implicit none
      ! Arguments
      character(len=*),intent(in) :: str
      integer,intent(out)         :: int
      integer,intent(out)         :: stat
      read(str,*,iostat=stat)  int
    end subroutine str2int
    ! #############################################################################
    ! Subroutine to compute interpolation weight for linear interpolation.
    ! #############################################################################
    subroutine compute_interpolation_wts(p,pi,xi1,xi2,wt)
      real,dimension(:),intent(in) :: p
      real,intent(in) :: pi
      integer,dimension(1),intent(out) :: xi1,xi2
      real,intent(out) :: wt
      xi1 = minloc(p-pi,p-pi.gt. 0)
      xi2 = xi1 + 1
      wt = (p(xi2(1))-pi)/(p(xi2(1))-p(xi1(1)))
    end subroutine compute_interpolation_wts
    
    ! #############################################################################
    ! Subroutines to initialize netCDF output fields
    ! #############################################################################
    subroutine init_ncdfOutVar1Dflt(fileID,varID,dimID1,varName,description,units)
      character(len=*),intent(in) :: &
           varName,     & ! Name for output variable
           description, & ! Metadata information
           units          ! Units
      integer,intent(in) :: fileID,dimID1
      integer,intent(inout) :: varID
      status = nf90_def_var(fileID,varName,nf90_float, (/dimID1/),varID)
      if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable,',varName
      status = nf90_put_att(fileID,varID,"FieldType",104)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable, ',varName
      status = nf90_put_att(fileID,varID,"MemoryOrder","XY")
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable, ',varName
      status = nf90_put_att(fileID,varID,"description",description)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable, ',varName
      status = nf90_put_att(fileID,varID,"units",units)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable, ',varName
      status = nf90_put_att(fileID,varID,"stagger","")
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable, ',varName
    end subroutine init_ncdfOutVar1Dflt
    subroutine init_ncdfOutVar1Dint(fileID,varID,dimID1,varName,description,units)
      character(len=*),intent(in) :: &
           varName,     & ! Name for output variable
           description, & ! Metadata information
           units          ! Units
      integer,intent(in) :: fileID,dimID1
      integer,intent(inout) :: varID
      status = nf90_def_var(fileID,varName,nf90_int, (/dimID1/),varID)
      if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable,',varName
      status = nf90_put_att(fileID,varID,"FieldType",104)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable, ',varName
      status = nf90_put_att(fileID,varID,"MemoryOrder","XY")
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable, ',varName
      status = nf90_put_att(fileID,varID,"description",description)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable, ',varName
      status = nf90_put_att(fileID,varID,"units",units)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable, ',varName
      status = nf90_put_att(fileID,varID,"stagger","")
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable, ',varName
    end subroutine init_ncdfOutVar1Dint
    subroutine init_ncdfOutVar2D(fileID,varID,dimID1,dimID2,varName,description,units)
      character(len=*),intent(in) :: &
           varName,     & ! Name for output variable
           description, & ! Metadata information
           units          ! Units
      integer,intent(in) :: fileID,dimID1,dimID2
      integer,intent(inout) :: varID
      status = nf90_def_var(fileID,varName,nf90_float, (/dimID1,dimID2/),varID)
      if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable,',varName
      status = nf90_put_att(fileID,varID,"FieldType",104)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable, ',varName
      status = nf90_put_att(fileID,varID,"MemoryOrder","XY")
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable, ',varName
      status = nf90_put_att(fileID,varID,"description",description)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable, ',varName
      status = nf90_put_att(fileID,varID,"units",units)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable, ',varName
      status = nf90_put_att(fileID,varID,"stagger","")
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable, ',varName
      status = nf90_put_att(fileID,varID,"FillValue",output_fill_value)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable, ',varName
    end subroutine init_ncdfOutVar2D
    subroutine init_ncdfOutVar3D(fileID,varID,dimID1,dimID2,dimID3,varName,description,units)
      character(len=*),intent(in) :: &
           varName,     & ! Name for output variable
           description, & ! Metadata information
           units          ! Units
      integer,intent(in) :: fileID,dimID1,dimID2,dimID3
      integer,intent(inout) :: varID
      status = nf90_def_var(fileID,varName,nf90_float, (/dimID1,dimID2,dimID3/),varID)
      if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable,',varName
      status = nf90_put_att(fileID,varID,"FieldType",104)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable, ',varName
      status = nf90_put_att(fileID,varID,"MemoryOrder","XY")
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable, ',varName
      status = nf90_put_att(fileID,varID,"description",description)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable, ',varName
      status = nf90_put_att(fileID,varID,"units",units)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable, ',varName
      status = nf90_put_att(fileID,varID,"stagger","")
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable, ',varName
      status = nf90_put_att(fileID,varID,"FillValue",output_fill_value)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable, ',varName
    end subroutine init_ncdfOutVar3D
    subroutine init_ncdfOutVar4D(fileID,varID,dimID1,dimID2,dimID3,dimID4,varName,description,units)
      character(len=*),intent(in) :: &
           varName,     & ! Name for output variable
           description, & ! Metadata information
           units          ! Units
      integer,intent(in) :: fileID,dimID1,dimID2,dimID3,dimID4
      integer,intent(inout) :: varID
      status = nf90_def_var(fileID,varName,nf90_float, (/dimID1,dimID2,dimID3,dimID4/),varID)
      if (status /= nf90_NoErr) print*,'ERROR: Failure defining output variable,',varName
      status = nf90_put_att(fileID,varID,"FieldType",104)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, FieldType, to variable, ',varName
      status = nf90_put_att(fileID,varID,"MemoryOrder","XY")
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, MemoryOrder, to variable, ',varName
      status = nf90_put_att(fileID,varID,"description",description)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, description, to variable, ',varName
      status = nf90_put_att(fileID,varID,"units",units)
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, units, to variable, ',varName
      status = nf90_put_att(fileID,varID,"stagger","")
      if (status /= nf90_NoErr) print*,'ERROR: Failure adding attribute, stagger, to variable, ',varName
    end subroutine init_ncdfOutVar4D
  end program WRF3D_pp
  
