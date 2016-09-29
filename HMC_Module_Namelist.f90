!--------------------------------------------------------------------------------  
! File:   HMC_Module_Namelist.f90
! Author(s): Fabio Delogu, Francesco Silvestro, Simone Gabellani
! Created on May, 20 2014, 9:57 AM
!
! Module to allocate and read namelist
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
! Module Header
module HMC_Module_Namelist
    
    !--------------------------------------------------------------------------------
    ! External module(s) and implicit none
    use HMC_Module_Tools_Debug  ! to import global variable(s) declaration
    
    implicit none
    
    include "HMC_Type_Namelist.inc"
    integer, parameter :: iMaxDomain=1
    type(HMC_Type_Namelist), dimension(iMaxDomain) :: oHMC_Namelist
    save oHMC_Namelist
    !--------------------------------------------------------------------------------
    
contains 
    
    !--------------------------------------------------------------------------------
    ! Subroutine to read namelist
    subroutine HMC_Namelist_Read(dUc, dUh, dCt, dCf, dCPI, dWTableHbr, dKSatRatio, dSlopeMax, &
                                 sDomainName, & 
                                 oHMC_Namelist_Init) 
        
        !--------------------------------------------------------------------------------
        ! External module(s) and implicit none
        implicit none

        type(HMC_Type_Namelist) oHMC_Namelist_Init
        !--------------------------------------------------------------------------------
        
        !--------------------------------------------------------------------------------
        ! Variable(s) declaration    
        integer(kind = 4)       :: iErr
        character(len = 256)    :: sFileName
        logical                 :: bFileExist
        
        integer(kind = 4)       :: iFlagDebugSet, iFlagDebugLevel
        
        integer(kind = 4)       :: iFlagOs
    
        integer(kind = 4)       :: iFlagTypeData_Static
        integer(kind = 4)       :: iFlagTypeData_Forcing_Gridded, iFlagTypeData_Forcing_Point
        integer(kind = 4)       :: iFlagTypeData_Output_Gridded, iFlagTypeData_Output_Point
        integer(kind = 4)       :: iFlagTypeData_State_Gridded, iFlagTypeData_State_Point
        integer(kind = 4)       :: iFlagTypeData_Restart_Gridded, iFlagTypeData_Restart_Point
        integer(kind = 4)       :: iFlagRestart, iFlagFlowDeep
        integer(kind = 4)       :: iFlagVarDtPhysConv, iFlagVarUc
        integer(kind = 4)       :: iFlagLAI, iFlagAlbedo, iFlagCH
        integer(kind = 4)       :: iFlagSnow, iFlagSnowAssim
        integer(kind = 4)       :: iFlagGrid
        
        logical                 :: bGridCheck

        integer(kind = 4)       :: iSimLength, iDtModel, iDtPhysConv
        integer(kind = 4)       :: iDtData_Forcing
        integer(kind = 4)       :: iDtData_Output_Gridded, iDtData_Output_Point
        integer(kind = 4)       :: iDtData_State_Gridded, iDtData_State_Point
        integer(kind = 4)       :: iDtPhysPrev
        integer(kind = 4)       :: iScaleFactor, iTcMax, iTc

        integer(kind = 4)       :: iRowsL, iColsL
        real(kind = 4)          :: dXLLCornerL, dYLLCornerL, dXCellSizeL, dYCellSizeL, dNoDataL
        
        integer(kind = 4)       :: iRowsF, iColsF
        real(kind = 4)          :: dXLLCornerF, dYLLCornerF, dXCellSizeF, dYCellSizeF, dNoDataF
        
        integer(kind = 4)       :: iNSection
        integer(kind = 4)       :: iNLake, iNDam, iNPlant, iNJoint, iNCatch, iNRelease
        integer(kind = 4)       :: iDaySteps, iTMarkedSteps
        
        integer(kind = 4)       :: iTdeepShift, iNTime, iETime
        integer(kind = 4)       :: iNData
        real(kind = 4)          :: dWTableHMin, dWTableHUSoil, dWTableHUChannel, dWTableSlopeBM, dWTableHOBedRock
        real(kind = 4)          :: dRateMin, dBc
        real(kind = 4)          :: dTRef, dEpsS, dSigma, dBFMin, dBFMax
        real(kind = 4)          :: dZRef, dG, dCp, dRd, dRhoS, dRhoW, dCpS, dCpW
        real(kind = 4)          :: dKq, dKw, dKo, dPorS, dFqS 
        real(kind = 4)          :: dTV, dDamSpillH
        
        integer(kind = 4)       :: iGlacierValue
        real(kind = 4)          :: dRhoSnowMax, dSnowQualityThr
        
        integer(kind = 4), target, dimension(3) :: a1iFlagS
        integer(kind = 4), target, dimension(3) :: a1iFlagO
        
        integer(kind = 4), target, dimension(2)     :: a1iDimsForcing
        real(kind = 4), target, dimension(2)        :: a1dGeoForcing
        real(kind = 4), target, dimension(2)        :: a1dResForcing
        
        real(kind = 4), target, dimension(4)        :: a1dArctUp
        real(kind = 4), target, dimension(4)        :: a1dExpRhoLow
        real(kind = 4), target, dimension(4)        :: a1dExpRhoHigh
        real(kind = 4), target, dimension(4)        :: a1dAltRange
                                        
        real(kind = 4), target, dimension(12)       :: a1dAlbedoMonthly
        real(kind = 4), target, dimension(12)       :: a1dLAIMonthly
        real(kind = 4), target, dimension(12)       :: a1dCHMonthly
        
        character(len = 12)     :: sTimeStart
        character(len = 12)     :: sTimeRestart
        character(len = 12)     :: sTimeStatus
        
        character(len = 19)     :: sTimeStartLong
        character(len = 19)     :: sTimeRestartLong
        character(len = 19)     :: sTimeStatusLong
        
        character(len = 256)    :: sPathData_Static_Gridded
        character(len = 256)    :: sPathData_Static_Point
        character(len = 256)    :: sPathData_Forcing_Gridded
        character(len = 256)    :: sPathData_Forcing_Point
        character(len = 256)    :: sPathData_Output_Gridded
        character(len = 256)    :: sPathData_Output_Point
        character(len = 256)    :: sPathData_State_Gridded
        character(len = 256)    :: sPathData_State_Point
        character(len = 256)    :: sPathData_Restart_Gridded
        character(len = 256)    :: sPathData_Restart_Point
        
        character(len = 1)      :: sPathBar
        character(len = 700)    :: sCommandUnzipFile
        character(len = 700)    :: sCommandZipFile
        character(len = 700)    :: sCommandRemoveFile
        character(len = 700)    :: sCommandCreateFolder
        
        character(len = 10)     :: sReleaseDate
        character(len = 700)    :: sAuthorNames
        character(len = 5)      :: sReleaseVersion
        
        real(kind = 4)          :: dUc, dUh, dCt, dCf, dCPI, dWTableHbr, dKSatRatio, dSlopeMax
        character(len = 256)    :: sDomainName
        
        character(len = 256)    :: sStrCf, sStrCt, sStrUh, sStrUc
        !--------------------------------------------------------------------------------

        !--------------------------------------------------------------------------------
        ! Read namelist(s)
        namelist /HMC_Namelist/         iFlagTypeData_Static, &
                                        iFlagTypeData_Forcing_Gridded, iFlagTypeData_Forcing_Point, &
                                        iFlagTypeData_Output_Gridded, iFlagTypeData_Output_Point, &
                                        iFlagTypeData_State_Gridded, iFlagTypeData_State_Point, &
                                        iFlagTypeData_Restart_Gridded, iFlagTypeData_Restart_Point, &
                                        iFlagDebugSet, iFlagDebugLevel, &
                                        iFlagOs, iFlagFlowDeep, iFlagRestart, &
                                        iFlagVarDtPhysConv, iFlagVarUc, &
                                        a1iFlagS, a1iFlagO, &
                                        iFlagLAI, iFlagAlbedo, iFlagCH, &
                                        iFlagSnow, iFlagSnowAssim, &
                                        a1dGeoForcing, a1dResForcing, a1iDimsForcing, &
                                        iScaleFactor, iTcMax, &
                                        iSimLength, iDtModel, iDtPhysConv, &
                                        iDtData_Forcing, &
                                        iDtData_Output_Gridded, iDtData_Output_Point, &
                                        iDtData_State_Gridded, iDtData_State_Point, &
                                        sTimeStart, sTimeStatus, sTimeRestart, &
                                        sPathData_Static_Gridded, sPathData_Static_Point, &
                                        sPathData_Forcing_Gridded, sPathData_Forcing_Point, &
                                        sPathData_Output_Gridded, sPathData_Output_Point, &
                                        sPathData_State_Gridded, sPathData_State_Point, &
                                        sPathData_Restart_Gridded, sPathData_Restart_Point
                                        
        namelist /HMC_Snow/             a1dArctUp, a1dExpRhoLow, a1dExpRhoHigh, a1dAltRange, &
                                        iGlacierValue, dRhoSnowMax, dSnowQualityThr
                                        
        namelist /HMC_Constants/        a1dAlbedoMonthly, a1dLAIMonthly, a1dCHMonthly, &
                                        dWTableHMin, dWTableHUSoil, dWTableHUChannel, dWTableSlopeBM, dWTableHOBedRock, &
                                        dRateMin, dBc, &
                                        dTRef, iTdeepShift, dEpsS, dSigma, dBFMin, dBFMax, &
                                        dZRef, dG, dCp, dRd, dRhoS, dRhoW, dCpS, dCpW, & 
                                        dKq, dKw, dKo, dPorS, dFqS, &
                                        dTV, dDamSpillH
                                        
        namelist /HMC_Command/          sCommandZipFile, &
                                        sCommandUnzipFile, &
                                        sCommandRemoveFile, &
                                        sCommandCreateFolder
        
        namelist /HMC_Info/             sReleaseDate, &
                                        sAuthorNames, &
                                        sReleaseVersion
        !--------------------------------------------------------------------------------
        
        !--------------------------------------------------------------------------------
        ! Set defaults for namelist variables
        iFlagTypeData_Static = -9999; 
        iFlagTypeData_Forcing_Gridded = -9999; iFlagTypeData_Forcing_Point = -9999; 
        iFlagTypeData_Output_Gridded = -9999; iFlagTypeData_Output_Point = -9999; 
        iFlagTypeData_State_Gridded = -9999; iFlagTypeData_State_Point = -9999; 
        iFlagTypeData_Restart_Gridded = -9999; iFlagTypeData_Restart_Point = -9999; 
        iFlagDebugSet = -9999; iFlagDebugLevel = -9999;
        iFlagOs = -9999; iFlagFlowDeep = -9999; iFlagRestart = -9999;
        iFlagVarDtPhysConv = -9999; iFlagVarUc = -9999; a1iFlagS = -9999; a1iFlagO = -9999;
        iFlagLAI = -9999; iFlagAlbedo = -9999; iFlagCH = -9999; 
        iFlagSnow = -9999; iFlagSnowAssim = -9999;
        a1dGeoForcing = -9999.0; a1dResForcing = -9999.0; a1iDimsForcing = -9999; 
        iScaleFactor = -9999; iTcMax = -9999; iTc = -9999
        iSimLength = -9999; iDtModel = -9999; iDtPhysConv = -9999; 
        iDtData_Forcing = -9999; 
        iDtData_Output_Gridded = -9999; iDtData_Output_Point = -9999; 
        iDtData_State_Gridded = -9999; iDtData_State_Point = -9999; 
        sTimeStart = ""; sTimeStatus = ""; sTimeRestart = ""; 
        sPathData_Static_Gridded = ""; sPathData_Static_Point = ""; 
        sPathData_Forcing_Gridded = ""; sPathData_Forcing_Point = ""; 
        sPathData_Output_Gridded = ""; sPathData_Output_Gridded = ""; 
        sPathData_State_Gridded = ""; sPathData_State_Gridded = ""; 
        sPathData_Restart_Gridded = ""; sPathData_Restart_Point = ""; 
        
        iNTime = 0; iNData = 0; iETime = 0;
        
        a1dArctUp = -9999.0; a1dExpRhoLow = -9999.0; a1dExpRhoHigh = -9999.0; a1dAltRange = -9999.0;
        iGlacierValue = -9999; dRhoSnowMax = -9999.0; dSnowQualityThr = -9999.0;
        
        a1dAlbedoMonthly = -9999.0; a1dLAIMonthly = -9999.0; a1dCHMonthly = -9999.0
        dWTableHMin = -9999.0; dWTableHUSoil = -9999.0; dWTableHUChannel = -9999.0; 
        dWTableSlopeBM = -9999.0; dWTableHOBedRock = -9999.0;
        dRateMin = -9999.0; dBc = -9999.0;
        dTRef = -9999.0; iTdeepShift = -9999; dEpsS = -9999.0; 
        dSigma = -9999.0; dBFMin = -9999.0; dBFMax = -9999.0
        dZRef = -9999.0; dG = -9999.0; dCp = -9999.0; dRd = -9999.0; dRhoS = -9999.0; 
        dRhoW = -9999.0; dCpS = -9999.0; dCpW = -9999.0; 
        dKq = -9999.0; dKw = -9999.0; dKo = -9999.0; dPorS = -9999.0; dFqS = -9999.0;
        dTV = -9999.0; dDamSpillH = -9999.0
        
        sCommandZipFile = ""; sCommandUnzipFile = ""; sCommandRemoveFile = ""; sCommandCreateFolder = ""
        
        sReleaseDate = ""; sAuthorNames = ""; sReleaseVersion = "";
        !--------------------------------------------------------------------------------  
                                        
        !--------------------------------------------------------------------------------           
        ! Define namelist file
        sFileName = trim(sDomainName)//'.info.txt' 

        ! Checking file input availability
        inquire (file = trim(sFileName), exist = bFileExist)
        if ( .not. bFileExist ) then
            write(6, *) "No file info found ", trim(sFileName)
            stop "Stopped"
        else
            open(20,file = trim(sFileName))
            read(20, HMC_Namelist, iostat=iErr)
            read(20, HMC_Snow, iostat=iErr)
            read(20, HMC_Constants, iostat=iErr)
            read(20, HMC_Command, iostat=iErr)
            read(20, HMC_Info, iostat=iErr)
            close(20) 
        endif
        !--------------------------------------------------------------------------------
        
        !--------------------------------------------------------------------------------
        ! Set debug level
        call HMC_Tools_Debug_SetLevel(iFlagDebugSet, iFlagDebugLevel)
        ! Set file unit debug
        call HMC_Tools_Debug_SetUnit(iDebugUnit)
        !--------------------------------------------------------------------------------

        !-----------------------------------------------------------------------------------
        ! HMC Version
        call mprintf(.true., iINFO_Basic, '************************************************************************')
        call mprintf(.true., iINFO_Basic, ' Hydrological Model Continuum (Version: '//trim(sReleaseVersion)//')')
        call mprintf(.true., iINFO_Basic, ' Author(s): '//trim(sAuthorNames))
        call mprintf(.true., iINFO_Basic, ' Release Date: '//trim(sReleaseDate))
        call mprintf(.true., iINFO_Basic, '************************************************************************')
        !-----------------------------------------------------------------------------------        
        
        !--------------------------------------------------------------------------------
        ! Time long definition
        sTimeStartLong = sTimeStart(1:4)//"-"//sTimeStart(5:6)//"-"//sTimeStart(7:8)//"_" &
                         //sTimeStart(9:10)//":"//sTimeStart(11:12)//":"//"00"
        
        sTimeRestartLong = sTimeRestart(1:4)//"-"//sTimeRestart(5:6)//"-"//sTimeRestart(7:8)//"_" &
                           //sTimeRestart(9:10)//":"//sTimeRestart(11:12)//":"//"00"
        
        sTimeStatusLong = sTimeStatus(1:4)//"-"//sTimeStatus(5:6)//"-"//sTimeStatus(7:8)//"_" &
                           //sTimeStatus(9:10)//":"//sTimeStatus(11:12)//":"//"00"
        !--------------------------------------------------------------------------------
                                   
        !--------------------------------------------------------------------------------
        ! Define OS features
        if (iFlagOS.eq.10) then
            call mprintf(.true., iINFO_Main, ' Operative System: GNU/Linux')
            sPathBar = '/'
        elseif (iFlagOS.eq.0) then
            call mprintf(.true., iINFO_Main, ' Operative System: WinOS')
            sPathBar = '\'
        else
            call mprintf(.true., iERROR, ' Incorrect OS definition (Allowed value: GNU/Linux = 10, WinOS = 0)')
        endif
        !--------------------------------------------------------------------------------               
        
        !------------------------------------------------------------------------------------------
        ! Compute simulation time length
        iNTime = (iSimLength + 1)*3600./nint(real(iDtModel))
        ! Compute data time steps
        iNData = int((iNTime*3600))/int(iDtData_Forcing)
        !------------------------------------------------------------------------------------------
        
        !--------------------------------------------------------------------------------  
        ! Check snow physics flag(s)
        if (iFlagSnow.eq.0) iFlagSnowAssim = 0;
        !--------------------------------------------------------------------------------  
        
        !--------------------------------------------------------------------------------
        ! Namelist definition
        call mprintf(.true., iINFO_Main, ' Read Namelist ... ')
        
        ! Model dims (initialization)
        oHMC_Namelist_Init%iRowsL = 0
        oHMC_Namelist_Init%iColsL = 0
        oHMC_Namelist_Init%iRowsF = 0
        oHMC_Namelist_Init%iColsF = 0
        
        oHMC_Namelist_Init%iNSection = 0
        
        oHMC_Namelist_Init%iNLake = 0
        oHMC_Namelist_Init%iNDam = 0
        oHMC_Namelist_Init%iNPlant = 0
        oHMC_Namelist_Init%iNJoint = 0
        oHMC_Namelist_Init%iNCatch = 0
        oHMC_Namelist_Init%iNRelease = 0
        oHMC_Namelist_Init%iDaySteps = 0
        oHMC_Namelist_Init%iTMarkedSteps = 0
        
        ! Flag(s) info
        oHMC_Namelist_Init%iFlagTypeData_Static = iFlagTypeData_Static
        oHMC_Namelist_Init%iFlagTypeData_Forcing_Gridded = iFlagTypeData_Forcing_Gridded
        oHMC_Namelist_Init%iFlagTypeData_Forcing_Point = iFlagTypeData_Forcing_Point
        oHMC_Namelist_Init%iFlagTypeData_Output_Gridded = iFlagTypeData_Output_Gridded
        oHMC_Namelist_Init%iFlagTypeData_Output_Point = iFlagTypeData_Output_Point
        oHMC_Namelist_Init%iFlagTypeData_State_Gridded = iFlagTypeData_State_Gridded
        oHMC_Namelist_Init%iFlagTypeData_State_Point = iFlagTypeData_State_Point
        oHMC_Namelist_Init%iFlagTypeData_Restart_Gridded = iFlagTypeData_Restart_Gridded
        oHMC_Namelist_Init%iFlagTypeData_Restart_Point = iFlagTypeData_Restart_Point
        
        oHMC_Namelist_Init%iFlagDebugLevel = iFlagDebugLevel
        oHMC_Namelist_Init%iFlagOs = iFlagOs
        oHMC_Namelist_Init%iFlagFlowDeep = iFlagFlowDeep
        oHMC_Namelist_Init%iFlagRestart = iFlagRestart
        oHMC_Namelist_Init%a1iFlagS = a1iFlagS
        oHMC_Namelist_Init%a1iFlagO = a1iFlagO
        oHMC_Namelist_Init%iFlagVarDtPhysConv = iFlagVarDtPhysConv
        oHMC_Namelist_Init%iFlagVarUc = iFlagVarUc
        oHMC_Namelist_Init%iFlagLAI = iFlagLAI
        oHMC_Namelist_Init%iFlagAlbedo = iFlagAlbedo
        oHMC_Namelist_Init%iFlagCH = iFlagCH
        oHMC_Namelist_Init%iFlagSnow = iFlagSnow
        oHMC_Namelist_Init%iFlagSnowAssim = iFlagSnowAssim
        oHMC_Namelist_Init%iFlagGrid = -9999

        ! Geographical land and forcing info
        oHMC_Namelist_Init%bGridCheck = .false.
        
        oHMC_Namelist_Init%dXLLCornerL = 0.0
        oHMC_Namelist_Init%dYLLCornerL = 0.0
        oHMC_Namelist_Init%dXCellSizeL = 0.0
        oHMC_Namelist_Init%dYCellSizeL = 0.0
        oHMC_Namelist_Init%dNoDataL = 0.0
        oHMC_Namelist_Init%dXLLCornerF = 0.0
        oHMC_Namelist_Init%dYLLCornerF = 0.0
        oHMC_Namelist_Init%dXCellSizeF = 0.0
        oHMC_Namelist_Init%dYCellSizeF = 0.0
        oHMC_Namelist_Init%dNoDataF = 0.0
        
        oHMC_Namelist_Init%a1dGeoForcing = a1dGeoForcing
        oHMC_Namelist_Init%a1dResForcing = a1dResForcing
        oHMC_Namelist_Init%a1iDimsForcing = a1iDimsForcing
        
        ! S3M parameter(s) and constant(s)
        oHMC_Namelist_Init%a1dArctUp = a1dArctUp
        oHMC_Namelist_Init%a1dExpRhoLow = a1dExpRhoLow
        oHMC_Namelist_Init%a1dExpRhoHigh = a1dExpRhoHigh
        oHMC_Namelist_Init%a1dAltRange = a1dAltRange
        oHMC_Namelist_Init%iGlacierValue = iGlacierValue
        oHMC_Namelist_Init%dRhoSnowMax = dRhoSnowMax
        oHMC_Namelist_Init%dSnowQualityThr = dSnowQualityThr

        ! Time, dt and step(s) info
        oHMC_Namelist_Init%iNTime = iNTime
        oHMC_Namelist_Init%iETime = 0
        oHMC_Namelist_Init%iTcMax = iTcMax
        oHMC_Namelist_Init%iTc = 0
        oHMC_Namelist_Init%iSimLength = iSimLength
        
        oHMC_Namelist_Init%iDtModel = iDtModel
        oHMC_Namelist_Init%iDtPhysConv = iDtPhysConv
        oHMC_Namelist_Init%iDtPhysConvPrevious = iDtPhysConv
        
        oHMC_Namelist_Init%iDtData_Forcing = iDtData_Forcing
        oHMC_Namelist_Init%iDtData_Output_Gridded = iDtData_Output_Gridded 
        oHMC_Namelist_Init%iDtData_Output_Point = iDtData_Output_Point                            
        oHMC_Namelist_Init%iDtData_State_Gridded = iDtData_State_Gridded
        oHMC_Namelist_Init%iDtData_State_Point = iDtData_State_Point
        
        ! Time reference info
        oHMC_Namelist_Init%sTimeStart = sTimeStartLong
        oHMC_Namelist_Init%sTimeStatus = sTimeStatusLong
        oHMC_Namelist_Init%sTimeRestart = sTimeRestartLong
        
        ! Data info
        oHMC_Namelist_Init%iNData = iNData
        oHMC_Namelist_Init%iScaleFactor = iScaleFactor
        
        ! Path(s) info
        oHMC_Namelist_Init%sPathData_Static_Gridded = sPathData_Static_Gridded
        oHMC_Namelist_Init%sPathData_Static_Point = sPathData_Static_Point
        oHMC_Namelist_Init%sPathData_Forcing_Gridded = sPathData_Forcing_Gridded
        oHMC_Namelist_Init%sPathData_Forcing_Point = sPathData_Forcing_Point
        oHMC_Namelist_Init%sPathData_Output_Gridded = sPathData_Output_Gridded
        oHMC_Namelist_Init%sPathData_Output_Point = sPathData_Output_Point
        oHMC_Namelist_Init%sPathData_State_Gridded = sPathData_State_Gridded
        oHMC_Namelist_Init%sPathData_State_Point = sPathData_State_Point
        oHMC_Namelist_Init%sPathData_Restart_Gridded = sPathData_Restart_Gridded
        oHMC_Namelist_Init%sPathData_Restart_Point = sPathData_Restart_Point
        
        ! Monthly value(s)
        oHMC_Namelist_Init%a1dAlbedoMonthly = a1dAlbedoMonthly
        oHMC_Namelist_Init%a1dLAIMonthly = a1dLAIMonthly
        oHMC_Namelist_Init%a1dCHMonthly = a1dCHMonthly
        
        ! Command line
        oHMC_Namelist_Init%sPathBar = sPathBar
        oHMC_Namelist_Init%sCommandZipFile = sCommandZipFile
        oHMC_Namelist_Init%sCommandUnzipFile = sCommandUnzipFile
        oHMC_Namelist_Init%sCommandRemoveFile = sCommandRemoveFile
        oHMC_Namelist_Init%sCommandCreateFolder = sCommandCreateFolder
        
        ! Model parameter(s)
        oHMC_Namelist_Init%dUc = dUc
        oHMC_Namelist_Init%dUh = dUh
        oHMC_Namelist_Init%dCt = dCt
        oHMC_Namelist_Init%dCf = dCf
        oHMC_Namelist_Init%dCPI = dCPI
        oHMC_Namelist_Init%dWTableHbr = dWTableHbr
        oHMC_Namelist_Init%dKSatRatio = dKSatRatio
        oHMC_Namelist_Init%dSlopeMax = dSlopeMax
        oHMC_Namelist_Init%sDomainName = sDomainName
        
        ! Water-table constant(s)
        oHMC_Namelist_Init%dWTableHMin = dWTableHMin
        oHMC_Namelist_Init%dWTableHUSoil = dWTableHUSoil
        oHMC_Namelist_Init%dWTableHUChannel = dWTableHUChannel
        oHMC_Namelist_Init%dWTableSlopeBM = dWTableSlopeBM
        oHMC_Namelist_Init%dWTableHOBedRock = dWTableHOBedRock
        
        ! Convolution constant(s)
        oHMC_Namelist_Init%dRateMin = dRateMin
        oHMC_Namelist_Init%dBc = dBc
        
        ! LSM constant(s)
        oHMC_Namelist_Init%dTRef = dTRef
        oHMC_Namelist_Init%iTdeepShift = iTdeepShift
        oHMC_Namelist_Init%dEpsS = dEpsS
        oHMC_Namelist_Init%dSigma = dSigma
        oHMC_Namelist_Init%dBFMin = dBFMin
        oHMC_Namelist_Init%dBFMax = dBFMax
        oHMC_Namelist_Init%dZRef = dZRef
        oHMC_Namelist_Init%dG = dG
        oHMC_Namelist_Init%dCp = dCp
        oHMC_Namelist_Init%dRd = dRd
        oHMC_Namelist_Init%dRhoS = dRhoS
        oHMC_Namelist_Init%dRhoW = dRhoW
        oHMC_Namelist_Init%dCpS = dCpS
        oHMC_Namelist_Init%dCpW = dCpW
        oHMC_Namelist_Init%dKq = dKq
        oHMC_Namelist_Init%dKw = dKw
        oHMC_Namelist_Init%dKo = dKo
        oHMC_Namelist_Init%dPorS = dPorS
        oHMC_Namelist_Init%dFqS = dFqS
        oHMC_Namelist_Init%dTV = dTV
        oHMC_Namelist_Init%dDamSpillH = dDamSpillH
        
        ! Model info
        oHMC_Namelist_Init%sReleaseDate = sReleaseDate
        oHMC_Namelist_Init%sAuthorNames = sAuthorNames
        oHMC_Namelist_Init%sReleaseVersion = sReleaseVersion
        
        ! Info model
        write(sStrUc, *) dUc; write(sStrUh, *) dUh; write(sStrCt, *) dCt; write(sStrCf, *) dCf;
        call mprintf(.true., iINFO_Basic, ' PARAMETER(S) INFO --- dUc: '//trim(sStrUc)//' - dUh: '//trim(sStrUh)// &
                    ' dCt: '//trim(sStrCt)//' dCf: '//trim(sStrCf) )           
        ! Info
        call mprintf(.true., iINFO_Main, ' Read Namelist ... OK')
        !--------------------------------------------------------------------------------
        
    end subroutine HMC_Namelist_Read
    !--------------------------------------------------------------------------------

end module HMC_Module_Namelist
!--------------------------------------------------------------------------------