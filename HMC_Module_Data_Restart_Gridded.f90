!------------------------------------------------------------------------------------------    
! File:   HMC_Module_Data_Restart_Gridded.f90
! Author(s): Fabio Delogu, Francesco Silvestro, Simone Gabellani
!
! Created on May 7, 2015, 1:27 PM
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
! Module Header
module HMC_Module_Data_Restart_Gridded
    
    !------------------------------------------------------------------------------------------
    ! External module(s) for all subroutine in this module
#ifdef LIB_NC
    use netcdf
#endif
    
    use HMC_Module_Namelist,        only:   oHMC_Namelist
    use HMC_Module_Vars_Loader,     only:   oHMC_Vars
    
    use HMC_Module_Tools_Debug
    
#ifdef LIB_NC
    use HMC_Module_Tools_IO,        only:   HMC_Tools_IO_Get2d_Binary, &
                                            HMC_Tools_IO_Get3d_Binary, &
                                            HMC_Tools_IO_Get2d_NC, &
                                            HMC_Tools_IO_Get3d_NC, &
                                            check
#else
    use HMC_Module_Tools_IO,        only:   HMC_Tools_IO_Get2d_Binary, &
                                            HMC_Tools_IO_Get3d_Binary                                        
#endif                                   
                                            
    
    use HMC_Module_Tools_Generic,   only:   HMC_Tools_Generic_ReplaceText, & 
                                            HMC_Tools_Generic_CreateFolder, &
                                            HMC_Tools_Generic_ZipFile, &
                                            HMC_Tools_Generic_UnzipFile, &
                                            HMC_Tools_Generic_RemoveFile, &
                                            transpose3Dvar, &
                                            checkdomainvar
    
    ! Implicit none for all subroutines in this module
    implicit none
    !------------------------------------------------------------------------------------------
    
contains
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to manage restart gridded data
    subroutine HMC_Data_Restart_Gridded_Cpl( iID, sTime, &
                                             iRowsStart, iRowsEnd, &
                                             iColsStart, iColsEnd, &
                                             iDaySteps, iTMarkedSteps)

        !------------------------------------------------------------------------------------------
        ! Variable(s)                                    
        integer(kind = 4)           :: iID
        integer(kind = 4)           :: iFlagRestart, iFlagSnow
        integer(kind = 4)           :: iRows, iCols
        integer(kind = 4)           :: iRowsStart, iColsStart, iRowsEnd, iColsEnd
        integer(kind = 4)           :: iDaySteps, iTMarkedSteps

        integer(kind = 4)           :: iFlagTypeData_Restart 
        integer(kind = 4)           :: iScaleFactor
        
        character(len = 19)         :: sTime
        character(len = 700)        :: sPathData_Restart
        character(len = 700)        :: sFileNameData_Restart, sFileNameData_Restart_Zip
        character(len = 700)        :: sCommandUnzipFile
        
        real(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1) :: a2dVarDEM
        
        real(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1)    :: a2dVarVTot, a2dVarVRet, &
                                                                                              a2dVarHydro, a2dVarRouting, &
                                                                                              a2dVarFlowDeep, &
                                                                                              a2dVarWTable, a2dVarLST, &
                                                                                              a2dVarLat, a2dVarLon
                                                                                           
        integer(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1) :: a2iVarAgeS
        real(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1)    :: a2dVarSWE, a2dVarAlbedoS, a2dVarRhoS
                                                                                           
        real(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1) :: a2dVarWTableUpd                         

        real(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1, iTMarkedSteps) :: a3dVarTaKMarked
        real(kind = 4), dimension(iRowsEnd - iRowsStart + 1, iColsEnd - iColsStart + 1, iDaySteps) :: a3dVarTaK24   
        
        logical                     :: bFileExist, bCheckRestart, bCheckRestartS
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarVTot = 0.0; a2dVarVRet = 0.0; a2dVarHydro = 0.0; a2dVarRouting = 0.0; a2dVarFlowDeep = 0.0; 
        a2dVarLST = 0.0; a2dVarWTable = 0.0; a3dVarTaKMarked = 0.0; a3dVarTaK24  = 0.0;
        a2dVarLat = 0.0; a2dVarLon = 0.0
        
        a2iVarAgeS = 0; a2dVarSWE = 0.0; a2dVarAlbedoS = 0.0; a2dVarRhoS = 0.0;
        
        a2dVarWTableUpd = 0.0;
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Defining iRows and iCols (output data)
        iRows = iRowsEnd - iRowsStart + 1
        iCols = iColsEnd - iColsStart + 1
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        iFlagRestart = oHMC_Namelist(iID)%iFlagRestart
        iFlagSnow = oHMC_Namelist(iID)%iFlagSnow
        sPathData_Restart = oHMC_Namelist(iID)%sPathData_Restart_Gridded
        iFlagTypeData_Restart = oHMC_Namelist(iID)%iFlagTypeData_Restart_Gridded
        iScaleFactor = oHMC_Namelist(iID)%iScaleFactor
        sCommandUnzipFile = oHMC_Namelist(iID)%sCommandUnzipFile
        ! Get glabal variable(s)
        a2dVarDem = oHMC_Vars(iID)%a2dDem
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Debug
        if (iDEBUG.gt.0) then
            call mprintf(.true., iINFO_Extra, ' ========= RESTART GRIDDED START =========== ')
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dVTot, oHMC_Vars(iID)%a2iMask, 'VTOT START') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dVRet, oHMC_Vars(iID)%a2iMask, 'VRET START') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dHydro, oHMC_Vars(iID)%a2iMask, 'HYDRO START') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dRouting, oHMC_Vars(iID)%a2iMask, 'ROUTING START') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dWTable, oHMC_Vars(iID)%a2iMask, 'WTABLE START') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dLST, oHMC_Vars(iID)%a2iMask, 'LST START') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a3dTaKMarked(:,:,1), oHMC_Vars(iID)%a2iMask, 'TAMk START') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a3dTaK24(:,:,1), oHMC_Vars(iID)%a2iMask, 'TA24 END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dFlowDeep, oHMC_Vars(iID)%a2iMask, 'FLOWDEEP START') )
            call mprintf(.true., iINFO_Extra, checkvar(real(oHMC_Vars(iID)%a2iAge), oHMC_Vars(iID)%a2iMask, 'AGES START') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dSWE, oHMC_Vars(iID)%a2iMask, 'SWE START') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dRhoS, oHMC_Vars(iID)%a2iMask, 'RHOS START') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dAlbedo_Snow, oHMC_Vars(iID)%a2iMask, 'ALBEDOS START') )
            call mprintf(.true., iINFO_Extra, '')
        endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Check restart flag value
        if (iFlagRestart == 1) then
        
            !------------------------------------------------------------------------------------------
            ! Replace general path with specific time feature(s)
            call HMC_Tools_Generic_ReplaceText(sPathData_Restart, '$yyyy', sTime(1:4))
            call HMC_Tools_Generic_ReplaceText(sPathData_Restart, '$mm', sTime(6:7))
            call HMC_Tools_Generic_ReplaceText(sPathData_Restart, '$dd', sTime(9:10))
            call HMC_Tools_Generic_ReplaceText(sPathData_Restart, '$HH', sTime(12:13))
            call HMC_Tools_Generic_ReplaceText(sPathData_Restart, '$MM', sTime(15:16))
            !------------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------------
            ! Subroutine for reading netCDF restart data 
            if (iFlagTypeData_Restart == 2) then

                !------------------------------------------------------------------------------------------
                ! Call subroutine to read data in netCDF format
#ifdef LIB_NC
                call HMC_Data_Restart_Gridded_NC(iID, &
                                        sPathData_Restart, &
                                        iRows, iCols, &
                                        iDaySteps, iTMarkedSteps, &
                                        sTime, iFlagSnow, &
                                        a2dVarVTot, a2dVarVRet, &
                                        a2dVarHydro, a2dVarRouting, &
                                        a2dVarFlowDeep, &
                                        a2dVarWTable, &
                                        a2dVarLST, a3dVarTaKMarked, a3dVarTaK24, &
                                        a2iVarAgeS, a2dVarSWE, a2dVarAlbedoS, a2dVarRhoS, &
                                        a2dVarLat, a2dVarLon, &
                                        bCheckRestart, bCheckRestartS)
#else
                ! Redefinition of forcing data flag (if netCDF library is not linked)
                iFlagTypeData_Restart = 1   
                call mprintf(.true., iWARN, ' ATTENTION: '// &
                                            'restart gridded data type selected was netCDF but library is not linked! '// &
                                            'Will be used data in binary format!')
#endif
                !------------------------------------------------------------------------------------------
                                            
            endif
            !------------------------------------------------------------------------------------------                                    

            !------------------------------------------------------------------------------------------
            ! Subroutine for reading restart data in binary format
            if (iFlagTypeData_Restart == 1) then

                !------------------------------------------------------------------------------------------
                ! Call subroutine to read data in binary format
                call HMC_Data_Restart_Gridded_Binary(iID, &
                                        sPathData_Restart, &
                                        iRows, iCols, &
                                        iDaySteps, iTMarkedSteps, &
                                        sTime, iFlagSnow, &
                                        a2dVarVTot, a2dVarVRet, &
                                        a2dVarHydro, a2dVarRouting, &
                                        a2dVarFlowDeep, &
                                        a2dVarWTable, &
                                        a2dVarLST, a3dVarTaKMarked, a3dVarTaK24, &
                                        a2iVarAgeS, a2dVarSWE, a2dVarAlbedoS, a2dVarRhoS, &
                                        bCheckRestart, bCheckRestartS)
                !------------------------------------------------------------------------------------------

            endif
            !------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------
            ! Check restart flag on data availability
            if (bCheckRestart .eqv. .true.) then

                !------------------------------------------------------------------------------------------
                ! Variable(s) conversion (Watertable)
                where(a2dVarDem.gt.0.0)
                    a2dVarWTableUpd = a2dVarDem - a2dVarWTable/1000
                endwhere
                
                ! Check limit(s)
                where(a2dVarWTableUpd.gt.a2dVarDem)
                    a2dVarWTableUpd = a2dVarDem
                endwhere
                where(a2dVarWTableUpd.lt.oHMC_Vars(iID)%a2dWTableMax)
                    a2dVarWTableUpd = oHMC_Vars(iID)%a2dWTableMax
                endwhere
                
                ! Check hydro initialization
                where (a2dVarHydro .lt. 0.0000001)
                    a2dVarHydro = 0.0000001
                endwhere
                where (a2dVarHydro .gt. 100000.0)
                    a2dVarHydro = 0.0000001
                endwhere
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! Check variable(s) domain
                a2dVarVTot = checkdomainvar(a2dVarVTot, oHMC_Vars(iID)%a2iMask,             -9999.0 )
                a2dVarVRet = checkdomainvar(a2dVarVRet, oHMC_Vars(iID)%a2iMask,             0.001 )
                a2dVarHydro = checkdomainvar(a2dVarHydro, oHMC_Vars(iID)%a2iMask,           0.0 )
                a2dVarRouting = checkdomainvar(a2dVarRouting, oHMC_Vars(iID)%a2iMask,       0.0 )
                a2dVarWTableUpd = checkdomainvar(a2dVarWTableUpd, oHMC_Vars(iID)%a2iMask,   -9999.0 )
                a2dVarLST = checkdomainvar(a2dVarLST, oHMC_Vars(iID)%a2iMask,               -9999.0 )
                a2dVarFlowDeep = checkdomainvar(a2dVarFlowDeep, oHMC_Vars(iID)%a2iMask,     0.0 )
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! Pass variable(s) to global workspace
                oHMC_Vars(iID)%a2dVTot = a2dVarVTot;
                oHMC_Vars(iID)%a2dVRet = a2dVarVRet;
                oHMC_Vars(iID)%a2dHydro = a2dVarHydro;
                oHMC_Vars(iID)%a2dRouting = a2dVarRouting;
                oHMC_Vars(iID)%a2dWTable = a2dVarWTableUpd;
                oHMC_Vars(iID)%a2dLST = a2dVarLST;
                oHMC_Vars(iID)%a2dFlowDeep = a2dVarFlowDeep;
                oHMC_Vars(iID)%a3dTaKMarked = a3dVarTaKMarked;
                oHMC_Vars(iID)%a3dTaK24 = a3dVarTaK24

                ! Info end
                call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded ... OK' )
                !------------------------------------------------------------------------------------------
                
            else
                !------------------------------------------------------------------------------------------
                ! Exit message for not using restart data
                call mprintf(.true., iINFO_Verbose, ' Restart flag selected but data are N/A (gridded data)')
                ! Info end
                call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded ... SKIPPED ' )
                !------------------------------------------------------------------------------------------
                
            endif
            !------------------------------------------------------------------------------------------  
            
            !------------------------------------------------------------------------------------------
            ! Check restart flag on snow data availability
            if (bCheckRestartS .eqv. .true.) then
            
                !------------------------------------------------------------------------------------------
                ! Pass variable(s) to global workspace
                oHMC_Vars(iID)%a2iAge = a2iVarAgeS
                oHMC_Vars(iID)%a2dSWE = a2dVarSWE
                oHMC_Vars(iID)%a2dRhoS = a2dVarRhoS
                oHMC_Vars(iID)%a2dAlbedo_Snow = a2dVarAlbedoS

                ! Info end
                call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded ... OK' )
                !------------------------------------------------------------------------------------------
                
            else
                !------------------------------------------------------------------------------------------
                ! Exit message for not using restart data
                call mprintf(.true., iINFO_Verbose, ' Restart flag selected but snow data are N/A (gridded data)')
                ! Info end
                call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded for snow data ... SKIPPED ' )
                !------------------------------------------------------------------------------------------
                
            endif
            !------------------------------------------------------------------------------------------  
            
        else
            !------------------------------------------------------------------------------------------
            ! Exit message for not using restart data
            call mprintf(.true., iINFO_Verbose, ' No restart run selected (gridded data)')
            ! Info end
            call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded ... SKIPPED ' )
            !------------------------------------------------------------------------------------------
        endif
        !------------------------------------------------------------------------------------------
       
        !------------------------------------------------------------------------------------------
        ! Debug
        if (iDEBUG.gt.0) then
            call mprintf(.true., iINFO_Extra, '')
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dVTot, oHMC_Vars(iID)%a2iMask, 'VTOT END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dVRet, oHMC_Vars(iID)%a2iMask, 'VRET END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dHydro, oHMC_Vars(iID)%a2iMask, 'HYDRO END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dRouting, oHMC_Vars(iID)%a2iMask, 'ROUTING END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dWTable, oHMC_Vars(iID)%a2iMask, 'WTABLE END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dLST, oHMC_Vars(iID)%a2iMask, 'LST END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a3dTaKMarked(:,:,1), oHMC_Vars(iID)%a2iMask, 'TAMk END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a3dTaK24(:,:,1), oHMC_Vars(iID)%a2iMask, 'TA24 END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dFlowDeep, oHMC_Vars(iID)%a2iMask, 'FLOWDEEP END') )
            call mprintf(.true., iINFO_Extra, checkvar(real(oHMC_Vars(iID)%a2iAge), oHMC_Vars(iID)%a2iMask, 'AGES END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dSWE, oHMC_Vars(iID)%a2iMask, 'SWE END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dRhoS, oHMC_Vars(iID)%a2iMask, 'RHOS END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dAlbedo_Snow, oHMC_Vars(iID)%a2iMask, 'ALBEDOS END') )
            call mprintf(.true., iINFO_Extra, ' ========= RESTART GRIDDED END =========== ')
        endif
        !------------------------------------------------------------------------------------------
        
    end subroutine HMC_Data_Restart_Gridded_Cpl
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to read netCDF data restart
#ifdef LIB_NC
    subroutine HMC_Data_Restart_Gridded_NC(iID, &
                                           sPathData_Restart, &
                                           iRows, iCols, &
                                           iDaySteps, iTMarkedSteps, &
                                           sTime, iFlagSnow, &
                                           a2dVarVTot, a2dVarVRet, &
                                           a2dVarHydro, a2dVarRouting, &
                                           a2dVarFlowDeep, &
                                           a2dVarWTable, &
                                           a2dVarLST, a3dVarTaKMarked, a3dVarTaK24, &
                                           a2iVarAgeS, a2dVarSWE, a2dVarAlbedoS, a2dVarRhoS, &
                                           a2dVarLat, a2dVarLon, &
                                           bCheckRestart, bCheckRestartS)
                                      
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)                       :: iID                  
                                  
        character(len = 256), intent(in)        :: sPathData_Restart
        character(len = 700)                    :: sFileNameData_Restart, sFileNameData_Restart_Zip
        character(len = 700)                    :: sCommandUnzipFile
        character(len = 256)                    :: sVarName
        integer(kind = 4), intent(in)           :: iRows, iCols
        integer(kind = 4), intent(in)           :: iDaySteps, iTMarkedSteps
        integer(kind = 4), intent(in)           :: iFlagSnow

        character(len = 19), intent(in)         :: sTime

        real(kind = 4), dimension(iCols, iRows)                                :: a2dVar
        real(kind = 4), dimension(iCols, iRows, iTMarkedSteps)                 :: a3dVar1
        real(kind = 4), dimension(iCols, iRows, iDaySteps)                     :: a3dVar2
        
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarVTot
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarVRet
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarHydro
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarRouting      
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarFlowDeep       
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarWTable
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarLST
        real(kind = 4), dimension(iRows, iCols, iTMarkedSteps), intent(out)    :: a3dVarTaKMarked
        real(kind = 4), dimension(iRows, iCols, iDaySteps),     intent(out)    :: a3dVarTaK24
        integer(kind = 4), dimension(iRows, iCols),             intent(out)    :: a2iVarAgeS
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarSWE
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarAlbedoS
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarRhoS
       
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarLat
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarLon
        
        logical, dimension(9)   :: a1bVarCheck
        logical, dimension(4)   :: a1bVarCheckS

        character(len = 256)    :: sVarUnits
        integer(kind = 4)       :: iErr
        integer(kind = 4)       :: iFileID
        
        logical                 :: bFileExist
        
        logical                 :: bCheckRestart, bCheckRestartS
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarVTot = -9999.0; a2dVarVRet = -9999.0; 
        a2dVarHydro = -9999.0; a2dVarRouting = -9999.0; 
        a2dVarFlowDeep = -9999.0; a2dVarWTable = -9999.0; 
        a2dVarLST = -9999.0; a3dVarTaKMarked = -9999.0; a3dVarTaK24 = -9999.0; 
        a2iVarAgeS = -9999; a2dVarSWE = -9999.0; a2dVarAlbedoS = -9999.0; a2dVarRhoS = -9999.0;
        a2dVarLat = -9999.0; a2dVarLon = -9999.0;
        
        a1bVarCheck = .false.; a1bVarCheckS = .false.; bCheckRestart = .false.; 
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sCommandUnzipFile = oHMC_Namelist(iID)%sCommandUnzipFile
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded :: NetCDF ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Filename restart (example: hmc.state.201404300000.nc)
        sFileNameData_Restart = trim(sPathData_Restart)//"hmc.state-grid."// &
        sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
        sTime(12:13)//sTime(15:16)// &
        ".nc"

        ! Info netCDF filename
        call mprintf(.true., iINFO_Basic, ' Get filename (restart gridded): '//trim(sFileNameData_Restart)//' ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Checking file input availability
        sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
        inquire (file = trim(sFileNameData_Restart)//'.gz', exist = bFileExist)
        if ( .not. bFileExist ) then
            !------------------------------------------------------------------------------------------
            ! Warning message
            call mprintf(.true., iWARN, ' No compressed restart netCDF data found: '//trim(sFileNameData_Restart_Zip) )
            call mprintf(.true., iINFO_Verbose, &
                         ' Get filename (restart gridded): '//trim(sFileNameData_Restart)//' ... FAILED' )
            a1bVarCheck = .false.
            !------------------------------------------------------------------------------------------
        else
            !------------------------------------------------------------------------------------------
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Restart, .true.)
            !------------------------------------------------------------------------------------------
        
            !------------------------------------------------------------------------------------------
            ! Opening netCDF file
            iErr = nf90_open(trim(sFileNameData_Restart), NF90_NOWRITE, iFileID)
            if (iErr /= 0) then
                
                !------------------------------------------------------------------------------------------
                ! Condition for no file restart found
                call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed netCDF file: '// &
                             trim(sFileNameData_Restart)//' --> Undefined restart data values' ) 
                call mprintf(.true., iINFO_Verbose, &
                            ' Get filename (restart gridded): '//trim(sFileNameData_Restart)//' ... FAILED' )
                
                ! Flag check restart 
                a1bVarCheck = .false.
                !------------------------------------------------------------------------------------------
                            
            else
                
                !------------------------------------------------------------------------------------------
                ! Condition for file restart found
                ! VTot
                sVarName = 'VTot';
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .true., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get volume total data failed! ')
                    a2dVarVTot = -9999.0;
                    a1bVarCheck(1) = .false.
                else
                    a2dVarVTot = transpose(a2dVar)
                    a1bVarCheck(1) = .true.
                endif

                ! VRet
                sVarName = 'VRet';
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .true., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get volume retention data failed! ')
                    a2dVarVRet = -9999.0;
                    a1bVarCheck(2) = .false.
                else
                    a2dVarVRet = transpose(a2dVar)
                    a1bVarCheck(2) = .true.
                endif

                ! HydroLevel
                sVarName = 'HydroLevel';
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .true., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get hydro level data failed! ')
                    a2dVarHydro = -9999.0;
                    a1bVarCheck(3) = .false.
                else
                    a2dVarHydro = transpose(a2dVar)
                    a1bVarCheck(3) = .true.
                endif

                ! Routing
                sVarName = 'Routing';
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .true., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get routing data failed! ')
                    a2dVarRouting = -9999.0;
                    a1bVarCheck(4) = .false.
                else
                    a2dVarRouting = transpose(a2dVar)
                    a1bVarCheck(4) = .true.
                endif

                ! DFE
                sVarName = 'DFE';
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .true., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get flow deep data failed! ')
                    a2dVarFlowDeep = -9999.0;
                    a1bVarCheck(5) = .false.
                else
                    a2dVarFlowDeep = transpose(a2dVar)
                    a1bVarCheck(5) = .true.
                endif

                ! WTLevel
                sVarName = 'WTLevel';
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .true., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get water-table level data failed! ')
                    a2dVarWTable = -9999.0;
                    a1bVarCheck(6) = .false.
                else
                    a2dVarWTable = transpose(a2dVar)
                    a1bVarCheck(6) = .true.
                endif

                ! LST
                sVarName = 'LST';
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .true., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get land surface temperature data failed! ')
                    a2dVarLST = -9999.0;
                    a1bVarCheck(7) = .false.
                else
                    a2dVarLST = transpose(a2dVar)
                    a1bVarCheck(7) = .true.
                endif

                ! Tmk
                sVarName = 'Tmk';
                call HMC_Tools_IO_Get3d_NC((sVarName), iFileID, a3dVar1, sVarUnits, iTMarkedSteps, iCols, iRows, .true., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get temperature marked data failed! ')
                    a3dVarTaKMarked = -9999.0;
                    a1bVarCheck(8) = .false.
                else
                    a3dVarTaKMarked = transpose3Dvar(a3dVar1)
                    a1bVarCheck(8) = .true.
                endif

                ! T24
                sVarName = 'T24';
                call HMC_Tools_IO_Get3d_NC((sVarName), iFileID, a3dVar2, sVarUnits, iDaySteps, iCols, iRows, .true., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get temperature last 24 hours data failed! ')
                    a3dVarTaK24 = -9999.0;
                    a1bVarCheck(9) = .false.
                else
                    a3dVarTaK24 = transpose3Dvar(a3dVar2)
                    a1bVarCheck(9) = .true.
                endif

                ! Snow variable(s)                
                if (iFlagSnow.eq.1) then
                    ! SWE
                    sVarName = 'SWE'
                    call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                    if(iErr /= 0) then
                        call mprintf(.true., iWARN, ' ATTENTION: get snow water equivalent data failed! ')
                        a2dVarSWE = -9999.0;
                        a1bVarCheckS(1) = .false.
                    else
                        a2dVarSWE = transpose(a2dVar)
                        a1bVarCheckS(1) = .true.
                    endif
                   
                    ! Snow density
                    sVarName = 'RhoS';
                    call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                    if(iErr /= 0) then
                        call mprintf(.true., iWARN, ' ATTENTION: get snow density data failed! ')
                        a2dVarRhoS = -9999.0;
                        a1bVarCheckS(2) = .false.
                    else
                        a2dVarRhoS = transpose(a2dVar)
                        a1bVarCheckS(2) = .true.
                    endif
                    
                    ! Snow albedo
                    sVarName = 'AlbedoS';
                    call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                    if(iErr /= 0) then
                        call mprintf(.true., iWARN, ' ATTENTION: get snow albedo data failed! ')
                        a2dVarAlbedoS = -9999.0;
                        a1bVarCheckS(3) = .false.
                    else
                        a2dVarAlbedoS = transpose(a2dVar)
                        a1bVarCheckS(3) = .true.
                    endif
                    
                    ! Snow age
                    sVarName = 'AgeS';
                    call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                    if(iErr /= 0) then
                        call mprintf(.true., iWARN, ' ATTENTION: get snow age data failed! ')
                        a2iVarAgeS = -9999;
                        a1bVarCheckS(4) = .false.
                    else
                        a2iVarAgeS = int(transpose(a2dVar))
                        a1bVarCheckS(4) = .true.
                    endif
                
                else
                    ! Condition snow not activated
                    a2dVarSWE = -9999.0; a1bVarCheckS(1) = .true.
                    a2dVarRhoS = -9999.0; a1bVarCheckS(2) = .true.
                    a2dVarAlbedoS = -9999.0; a1bVarCheckS(3) = .true.
                    a2iVarAgeS = -9999; a1bVarCheckS(4) = .true.
                endif

                ! Closing netcdf file (drops db)
                iErr = nf90_close(iFileID)
                !------------------------------------------------------------------------------------------
                    
                !------------------------------------------------------------------------------------------
                ! Info filename
                call mprintf(.true., iINFO_Basic, ' Get filename (restart gridded): '//trim(sFileNameData_Restart)//' ... OK' )
                ! Info end
                call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded :: NetCDF ... OK' )
                !------------------------------------------------------------------------------------------

            endif
            !------------------------------------------------------------------------------------------
            
        endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Check restart
        if (all(a1bVarCheck .eqv. .true.) ) then
            call mprintf(.true., iINFO_Basic, ' Data :: Restart gridded :: NetCDF :: All variable(s) are loaded! ' )
            bCheckRestart = .true.
        else
            call mprintf(.true., iINFO_Basic, ' Data :: Restart gridded :: NetCDF :: Some/All variable(s) are N/A! ' )
            call mprintf(.true., iWARN, ' ATTENTION: restart flag activated but some data restart are not available! ')
            call mprintf(.true., iWARN, ' ATTENTION: restart conditions are null! ')
            bCheckRestart = .false.
        endif
        
        ! Check restart snow
        if (all(a1bVarCheckS .eqv. .true.) ) then
            call mprintf(.true., iINFO_Verbose, ' Data :: Restart gridded :: NetCDF :: All snow variable(s) are loaded! ' )
            bCheckRestartS = .true.
        else
            call mprintf(.true., iINFO_Verbose, ' Data :: Restart gridded :: NetCDF :: Some/All snow variable(s) are N/A! ' )
            call mprintf(.true., iWARN, ' ATTENTION: restart flag activated but some data snow restart are not available! ')
            call mprintf(.true., iWARN, ' ATTENTION: restart snow conditions are null! ')
            bCheckRestartS = .false.
        endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Debug
        if (iDEBUG.gt.0) then
            call mprintf(.true., iINFO_Extra, ' ========= CHECK FORCING GRIDDED NC =========== ')
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarVTot, int(oHMC_Vars(iID)%a2dDEM), 'VTOT NC') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarVRet, oHMC_Vars(iID)%a2iMask, 'VRET NC') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarHydro, oHMC_Vars(iID)%a2iMask, 'HYDRO NC') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRouting, oHMC_Vars(iID)%a2iMask, 'ROUTING NC') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarWTable, oHMC_Vars(iID)%a2iMask, 'WTABLE NC') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarLST, oHMC_Vars(iID)%a2iMask, 'LST NC') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarFlowDeep, oHMC_Vars(iID)%a2iMask, 'FLOWDEEP NC') )
            call mprintf(.true., iINFO_Extra, checkvar(real(a2iVarAgeS), oHMC_Vars(iID)%a2iMask, 'AGES NC') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWE, oHMC_Vars(iID)%a2iMask, 'SWE NC') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRhoS, oHMC_Vars(iID)%a2iMask, 'RHOS NC') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarAlbedoS, oHMC_Vars(iID)%a2iMask, 'ALBEDOS NC') )
            call mprintf(.true., iINFO_Extra, ' ========= CHECK FORCING GRIDDED NC =========== ')
        endif
        !------------------------------------------------------------------------------------------
        
    end subroutine HMC_Data_Restart_Gridded_NC
#endif
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to read binary data restart
    subroutine HMC_Data_Restart_Gridded_Binary(iID, &
                                               sPathData_Restart, &
                                               iRows, iCols, &
                                               iDaySteps, iTMarkedSteps, &
                                               sTime, iFlagSnow, &
                                               a2dVarVTot, a2dVarVRet, &
                                               a2dVarHydro, a2dVarRouting, &
                                               a2dVarFlowDeep, &
                                               a2dVarWTable, &
                                               a2dVarLST, a3dVarTaKMarked, a3dVarTaK24, &
                                               a2iVarAgeS, a2dVarSWE, a2dVarAlbedoS, a2dVarRhoS, &
                                               bCheckRestart, bCheckRestartS)
    
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)                       :: iID                  
                                  
        character(len = 700), intent(in)        :: sPathData_Restart
        character(len = 700)                    :: sFileNameData_Restart, sFileNameData_Restart_Zip
        character(len = 700)                    :: sCommandUnzipFile
        character(len = 256)                    :: sVarName
        integer(kind = 4), intent(in)           :: iRows, iCols
        integer(kind = 4), intent(in)           :: iDaySteps, iTMarkedSteps
        integer(kind = 4), intent(in)           :: iFlagSnow
        
        integer(kind = 4)                       :: iScaleFactor
        character(len = 19), intent(in)         :: sTime

        real(kind = 4), dimension(iRows, iCols)                                :: a2dVar
        real(kind = 4), dimension(iRows, iCols, iTMarkedSteps)                 :: a3dVar1
        real(kind = 4), dimension(iRows, iCols, iDaySteps)                     :: a3dVar2
        
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarVTot
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarVRet
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarHydro
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarRouting      
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarFlowDeep       
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarWTable
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarLST
        real(kind = 4), dimension(iRows, iCols, iTMarkedSteps), intent(out)    :: a3dVarTaKMarked
        real(kind = 4), dimension(iRows, iCols, iDaySteps),     intent(out)    :: a3dVarTaK24
        integer(kind = 4), dimension(iRows, iCols),             intent(out)    :: a2iVarAgeS
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarSWE
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarAlbedoS
        real(kind = 4), dimension(iRows, iCols),                intent(out)    :: a2dVarRhoS
        
        logical, dimension(9)  :: a1bVarCheck
        logical, dimension(4)   :: a1bVarCheckS
        
        character(len = 256)    :: sVarUnits
        integer(kind = 4)       :: iErr
        integer(kind = 4)       :: iFileID
        
        logical                 :: bFileExist, bCheckRestart, bCheckRestartS
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarVTot = -9999.0; a2dVarVRet = -9999.0; a2dVarHydro = -9999.0; 
        a2dVarRouting = -9999.0; a2dVarFlowDeep = -9999.0; a2dVarWTable = -9999.0; 
        a2dVarLST = -9999.0; a3dVarTaKMarked = -9999.0; a3dVarTaK24 = -9999.0; 
        a2iVarAgeS = -9999; a2dVarSWE = -9999.0; a2dVarAlbedoS = -9999.0; a2dVarRhoS = -9999.0;
        
        sFileNameData_Restart = ""; sCommandUnzipFile = "";
        
        a1bVarCheck = .false.
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sCommandUnzipFile = oHMC_Namelist(iID)%sCommandUnZipFile
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded :: Binary ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info filename(s) at each time step
        call mprintf(.true., iINFO_Basic, ' Get (restart gridded) at time '//trim(sTime)//' ... ')
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! VTot  (example: V_201405010000.bin.gz)
        iScaleFactor = 10000
        sFileNameData_Restart = trim(sPathData_Restart)//"V_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"  
        call mprintf(.true., iINFO_Extra, ' Get filename: '//trim(sFileNameData_Restart) )

        ! Checking file input availability
        sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
        inquire (file = sFileNameData_Restart_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Restart_Zip)//' --> Undefined restart data values' )
            a2dVar = -9999.0
            a1bVarCheck(1) = .false.
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Restart, .true.)
            ! Read binary data
            a2dVar = -9999.0
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Restart, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            a1bVarCheck(1) = .true.
        endif
        a2dVarVTot = a2dVar
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! VRet  (example: Ret_201405010000.bin.gz)
        iScaleFactor = 10000
        sFileNameData_Restart = trim(sPathData_Restart)//"Ret_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"  
        call mprintf(.true., iINFO_Extra, ' Get filename: '//trim(sFileNameData_Restart) )

        ! Checking file input availability
        sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
        inquire (file = sFileNameData_Restart_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Restart_Zip)//' --> Undefined restart data values' )
            a2dVar = -9999.0
            a1bVarCheck(2) = .false.
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Restart, .true.)
            ! Read binary data
            a2dVar = -9999.0
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Restart, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            a1bVarCheck(2) = .true.
        endif
        a2dVarVRet = a2dVar
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! HydroLevel  (example: Wl_201405010000.bin.gz) 
        iScaleFactor = 100000
        sFileNameData_Restart = trim(sPathData_Restart)//"Wl_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"  
        call mprintf(.true., iINFO_Extra, ' Get filename: '//trim(sFileNameData_Restart) )

        ! Checking file input availability
        sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
        inquire (file = sFileNameData_Restart_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Restart_Zip)//' --> Undefined restart data values' )
            a2dVar = -9999.0
            a1bVarCheck(3) = .false.
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Restart, .true.)
            ! Read binary data
            a2dVar = -9999.0
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Restart, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            a1bVarCheck(3) = .true.
        endif
        a2dVarHydro = a2dVar
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Routing  (example: Rou_201405010000.bin.gz)
        iScaleFactor = 100000
        sFileNameData_Restart = trim(sPathData_Restart)//"Rou_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"  
        call mprintf(.true., iINFO_Extra, ' Get filename: '//trim(sFileNameData_Restart) )

        ! Checking file input availability
        sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
        inquire (file = sFileNameData_Restart_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Restart_Zip)//' --> Undefined restart data values' )
            a2dVar = -9999.0
            a1bVarCheck(4) = .false.
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Restart, .true.)
            ! Read binary data
            a2dVar = -9999.0
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Restart, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            a1bVarCheck(4) = .true.
        endif
        a2dVarRouting = a2dVar
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! DFE  (example: DFE_201405010000.bin.gz)
        iScaleFactor = 10000
        sFileNameData_Restart = trim(sPathData_Restart)//"DFE_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"  
        call mprintf(.true., iINFO_Extra, ' Get filename: '//trim(sFileNameData_Restart) )

        ! Checking file input availability
        sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
        inquire (file = sFileNameData_Restart_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Restart_Zip)//' --> Undefined restart data values' )
            a2dVar = -9999.0
            a1bVarCheck(5) = .false.
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Restart, .true.)
            ! Read binary data
            a2dVar = -9999.0
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Restart, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            a1bVarCheck(5) = .true.
        endif
        a2dVarFlowDeep = a2dVar
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! WTLevel  (example: Vw_201405010000.bin.gz)
        iScaleFactor = 10000
        sFileNameData_Restart = trim(sPathData_Restart)//"Vw_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"  
        call mprintf(.true., iINFO_Extra, ' Get filename: '//trim(sFileNameData_Restart) )

        ! Checking file input availability
        sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
        inquire (file = sFileNameData_Restart_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Restart_Zip)//' --> Undefined restart data values' )
            a2dVar = -9999.0
            a1bVarCheck(6) = .false.
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Restart, .true.)
            ! Read binary data
            a2dVar = -9999.0
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Restart, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            a1bVarCheck(6) = .true.
        endif
        a2dVarWTable = a2dVar
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! LST  (example: Ts_201405010000.bin.gz)
        iScaleFactor = 10000
        sFileNameData_Restart = trim(sPathData_Restart)//"Ts_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"  
        call mprintf(.true., iINFO_Extra, ' Get filename: '//trim(sFileNameData_Restart) )

        ! Checking file input availability
        sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
        inquire (file = sFileNameData_Restart_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Restart_Zip)//' --> Undefined restart data values' )
            a2dVar = -9999.0
            a1bVarCheck(7) = .false.
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Restart, .true.)
            ! Read binary data
            a2dVar = -9999.0
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Restart, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            a1bVarCheck(7) = .true.
        endif
        a2dVarLST = a2dVar
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! TMarked  (example: Tmk_201405010000.bin.gz)
        iScaleFactor = 10000
        sFileNameData_Restart = trim(sPathData_Restart)//"Tmk_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"  
        call mprintf(.true., iINFO_Extra, ' Get filename: '//trim(sFileNameData_Restart) )

        ! Checking file input availability
        sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
        inquire (file = sFileNameData_Restart_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Restart_Zip)//' --> Undefined restart data values' )
            a3dVar1 = -9999.0
            a1bVarCheck(8) = .false.
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Restart, .true.)
            ! Read binary data
            a3dVar1 = -9999.0
            call HMC_Tools_IO_Get3d_Binary(sFileNameData_Restart, a3dVar1, &
                                           iRows, iCols, iTMarkedSteps, iScaleFactor, .true., iErr) 
            a1bVarCheck(8) = .true.
        endif
        a3dVarTaKMarked = a3dVar1
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! T24  (example: T24_201405010000.bin.gz)
        iScaleFactor = 10000
        sFileNameData_Restart = trim(sPathData_Restart)//"T24_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"  
        call mprintf(.true., iINFO_Extra, ' Get filename: '//trim(sFileNameData_Restart) )

        ! Checking file input availability
        sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
        inquire (file = sFileNameData_Restart_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Restart_Zip)//' --> Undefined restart data values' )
            a3dVar2 = -9999.0
            a1bVarCheck(9) = .false.
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Restart, .true.)
            ! Read binary data
            a3dVar2 = -9999.0
            call HMC_Tools_IO_Get3d_Binary(sFileNameData_Restart, a3dVar2, &
                                           iRows, iCols, iDaySteps, iScaleFactor, .true., iErr) 
            a1bVarCheck(9) = .true.
        endif
        a3dVarTaK24 = a3dVar2
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Snow variable(s)                
        if (iFlagSnow.eq.1) then
            
            !------------------------------------------------------------------------------------------
            ! SWE  (example: SWE_201405010000.bin.gz)
            iScaleFactor = 1
            sFileNameData_Restart = trim(sPathData_Restart)//"SWE_"// &
                sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
                sTime(12:13)//sTime(15:16)// &
                ".bin"  
            call mprintf(.true., iINFO_Extra, ' Get filename: '//trim(sFileNameData_Restart) )

            ! Checking file input availability
            sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
            inquire (file = sFileNameData_Restart_Zip, exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                             trim(sFileNameData_Restart_Zip)//' --> Undefined restart data values' )
                a2dVar = -9999.0
                a1bVarCheckS(1) = .false.
            else
                ! Unzip file
                call HMC_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Restart, .true.)
                ! Read binary data
                a2dVar = -9999.0
                call HMC_Tools_IO_Get2d_Binary(sFileNameData_Restart, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
                a1bVarCheckS(1) = .true.
            endif
            a2dVarSWE = a2dVar
            !------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------
            ! Snow density  (example: Density_201405010000.bin.gz)
            iScaleFactor = 1
            sFileNameData_Restart = trim(sPathData_Restart)//"Density_"// &
                sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
                sTime(12:13)//sTime(15:16)// &
                ".bin"  
            call mprintf(.true., iINFO_Extra, ' Get filename: '//trim(sFileNameData_Restart) )

            ! Checking file input availability
            sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
            inquire (file = sFileNameData_Restart_Zip, exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                             trim(sFileNameData_Restart_Zip)//' --> Undefined restart data values' )
                a2dVar = -9999.0
                a1bVarCheckS(2) = .false.
            else
                ! Unzip file
                call HMC_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Restart, .true.)
                ! Read binary data
                a2dVar = -9999.0
                call HMC_Tools_IO_Get2d_Binary(sFileNameData_Restart, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
                a1bVarCheckS(2) = .true.
            endif
            a2dVarRhoS = a2dVar
            !------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------
            ! Snow albedo  (example: Density_201405010000.bin.gz)
            iScaleFactor = 1
            sFileNameData_Restart = trim(sPathData_Restart)//"AlbedoS_"// &
                sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
                sTime(12:13)//sTime(15:16)// &
                ".bin"  
            call mprintf(.true., iINFO_Extra, ' Get filename: '//trim(sFileNameData_Restart) )

            ! Checking file input availability
            sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
            inquire (file = sFileNameData_Restart_Zip, exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                             trim(sFileNameData_Restart_Zip)//' --> Undefined restart data values' )
                a2dVar = -9999.0
                a1bVarCheckS(3) = .false.
            else
                ! Unzip file
                call HMC_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Restart, .true.)
                ! Read binary data
                a2dVar = -9999.0
                call HMC_Tools_IO_Get2d_Binary(sFileNameData_Restart, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
                a1bVarCheckS(3) = .true.
            endif
            a2dVarAlbedoS = a2dVar
            !------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------
            ! Snow age  (example: Age_201405010000.bin.gz)
            iScaleFactor = 1
            sFileNameData_Restart = trim(sPathData_Restart)//"Age_"// &
                sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
                sTime(12:13)//sTime(15:16)// &
                ".bin"  
            call mprintf(.true., iINFO_Extra, ' Get filename: '//trim(sFileNameData_Restart) )

            ! Checking file input availability
            sFileNameData_Restart_Zip = trim(sFileNameData_Restart)//'.gz'
            inquire (file = sFileNameData_Restart_Zip, exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                             trim(sFileNameData_Restart_Zip)//' --> Undefined restart data values' )
                a2dVar = -9999.0
                a1bVarCheckS(4) = .false.
            else
                ! Unzip file
                call HMC_Tools_Generic_UnzipFile(sCommandUnzipFile, &
                                             sFileNameData_Restart_Zip, &
                                             sFileNameData_Restart, .true.)
                ! Read binary data
                a2dVar = -9999.0
                call HMC_Tools_IO_Get2d_Binary(sFileNameData_Restart, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
                a1bVarCheckS(4) = .true.
            endif
            a2iVarAgeS = int(a2dVar)
            !------------------------------------------------------------------------------------------
            
        else
            
            !------------------------------------------------------------------------------------------
            ! Condition(s) if snow not activated
            a2dVarSWE = -9999.0; a1bVarCheckS(1) = .true.
            a2dVarRhoS = -9999.0; a1bVarCheckS(2) = .true.
            a2dVarAlbedoS = -9999.0; a1bVarCheckS(3) = .true.
            a2iVarAgeS = -9999; a1bVarCheckS(4) = .true.
            !------------------------------------------------------------------------------------------
            
        endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info filename(s) at each time step
        call mprintf(.true., iINFO_Basic, ' Get (restart gridded) at time '//trim(sTime)//' ... OK')
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Restart gridded :: Binary ... OK' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Check restart
        if (all(a1bVarCheck .eqv. .true.) ) then
            call mprintf(.true., iINFO_Basic, ' Data :: Restart gridded :: Binary :: All variable(s) are loaded! ' )
            bCheckRestart = .true.
        else
            call mprintf(.true., iINFO_Basic, ' Data :: Restart gridded :: Binary :: Some/All variable(s) are N/A! ' )
            call mprintf(.true., iWARN, ' ATTENTION: restart flag activated but data restart are not available! ')
            call mprintf(.true., iWARN, ' ATTENTION: restart conditions are null! ')
            bCheckRestart = .false.
        endif
        ! Check restart snow
        if (all(a1bVarCheckS .eqv. .true.) ) then
            call mprintf(.true., iINFO_Verbose, ' Data :: Restart gridded :: Binary :: All snow variable(s) are loaded! ' )
            bCheckRestartS = .true.
        else
            call mprintf(.true., iINFO_Verbose, ' Data :: Restart gridded :: Binary :: Some/All snow variable(s) are N/A! ' )
            call mprintf(.true., iWARN, ' ATTENTION: restart flag activated but snow data restart are not available! ')
            call mprintf(.true., iWARN, ' ATTENTION: restart snow conditions are null! ')
            bCheckRestartS = .false.
        endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Debug
        if (iDEBUG.gt.0) then
            call mprintf(.true., iINFO_Extra, ' ========= CHECK FORCING GRIDDED BINARY =========== ')
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarVTot, int(oHMC_Vars(iID)%a2dDEM), 'VTOT BIN') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarVRet, oHMC_Vars(iID)%a2iMask, 'VRET BIN') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarHydro, oHMC_Vars(iID)%a2iMask, 'HYDRO BIN') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRouting, oHMC_Vars(iID)%a2iMask, 'ROUTING BIN') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarWTable, oHMC_Vars(iID)%a2iMask, 'WTABLE BIN') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarLST, oHMC_Vars(iID)%a2iMask, 'LST BIN') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarFlowDeep, oHMC_Vars(iID)%a2iMask, 'FLOWDEEP BIN') )
            call mprintf(.true., iINFO_Extra, checkvar(real(a2iVarAgeS), oHMC_Vars(iID)%a2iMask, 'AGES BIN') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWE, oHMC_Vars(iID)%a2iMask, 'SWE BIN') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarRhoS, oHMC_Vars(iID)%a2iMask, 'RHOS BIN') )
            call mprintf(.true., iINFO_Extra, checkvar(a2dVarAlbedoS, oHMC_Vars(iID)%a2iMask, 'ALBEDOS BIN') )
            call mprintf(.true., iINFO_Extra, ' ========= CHECK FORCING GRIDDED BINARY =========== ')
        endif
        !------------------------------------------------------------------------------------------
        
    end subroutine HMC_Data_Restart_Gridded_Binary
    !------------------------------------------------------------------------------------------
    
end module HMC_Module_Data_Restart_Gridded
!------------------------------------------------------------------------------------------