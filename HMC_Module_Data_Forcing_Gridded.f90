!------------------------------------------------------------------------------------------    
! File:   HMC_Module_Data_Forcing_Gridded.f90
! Author(s): Fabio Delogu, Francesco Silvestro, Simone Gabellani
!
! Created on April 22, 2015, 5:19 PM
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
! Module Header
module HMC_Module_Data_Forcing_Gridded
    
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
                                            HMC_Tools_IO_Get2d_NC, &
                                            check
#else
    use HMC_Module_Tools_IO,        only:   HMC_Tools_IO_Get2d_Binary                                      
#endif                                    
                                                                                  
    use HMC_Module_Tools_Generic,   only:   HMC_Tools_Generic_ReplaceText, &
                                            HMC_Tools_Generic_SwitchGrid, &
                                            HMC_Tools_Generic_UnzipFile, &
                                            HMC_Tools_Generic_RemoveFile, &
                                            check2Dvar
                                            
    use HMC_Module_Tools_Time,      only:   HMC_Tools_Time_MonthVal
                             
    ! Implicit none for all subroutines in this module
    implicit none
    !------------------------------------------------------------------------------------------
    
contains
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to manage forcing gridded data
    subroutine HMC_Data_Forcing_Gridded_Cpl( iID, sTime, &
                                     iRowsStartL, iRowsEndL, iColsStartL, iColsEndL, &
                                     iRowsStartF, iRowsEndF, iColsStartF, iColsEndF)
        
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)           :: iID
                                    
        integer(kind = 4)           :: iRowsStartL, iRowsEndL, iColsStartL, iColsEndL
        integer(kind = 4)           :: iRowsStartF, iRowsEndF, iColsStartF, iColsEndF
        integer(kind = 4)           :: iRowsL, iColsL, iRowsF, iColsF
        integer(kind = 4)           :: iFlagTypeData_Forcing
        integer(kind = 4)           :: iScaleFactor
        
        character(len = 19)         :: sTime
        character(len = 256)        :: sPathData_Forcing
        
        real(kind = 4), dimension(iRowsEndL - iRowsStartL + 1, &
                                  iColsEndL - iColsStartL + 1) ::   a2dVarRainL, a2dVarTaL, &
                                                                    a2dVarIncRadL, a2dVarWindL, & 
                                                                    a2dVarRelHumL, a2dVarPaL, &
                                                                    a2dVarAlbedoL, a2dVarLAIL, &
                                                                    a2dVarSnowHeightL, a2dVarSnowKernelL, &
                                                                    a2dVarSnowCAL, a2dVarSnowQAL
        real(kind = 4), dimension(iRowsEndL - iRowsStartL + 1, &
                                  iColsEndL - iColsStartL + 1) ::   a2dVarRainF, a2dVarTaF, &
                                                                    a2dVarIncRadF, a2dVarWindF, & 
                                                                    a2dVarRelHumF, a2dVarPaF, &
                                                                    a2dVarAlbedoF, a2dVarLAIF, &
                                                                    a2dVarSnowHeightF, a2dVarSnowKernelF, &
                                                                    a2dVarSnowCAF, a2dVarSnowQAF                          
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Forcing data (mandatory):
        !   a2dVarRainL         : rain [mm]
        !   a2dVarPaL           : atmospheric pressure [kPa]
        !   a2dVarTaL           : air temperature [C]
        !   a2dVarIncRadL       : incoming radiation [W/m^2]
        !   a2dVarWindL         : wind speed [m/s]
        !   a2dVarRelHumL       : relative humidity [%]
        
        ! Forcing data (optional):
        !   a2dVarAlbedoF       : albedo [0,1]
        !   a2dVarLAIF          : leaf area index [0,8]
        !   a2dVarSnowHeightF   : snow height [cm]                                                    
        !   a2dVarSnowKernelF   : snow kernel [0,1]   
        !   a2dVarSnowCAF       : snow cover area [-2,2] 
        !   a2dVarSnowQAF       : snow cover quality [0,1] 
        !------------------------------------------------------------------------------------------
                                                                                                                        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarRainF = -9999.0; a2dVarTaF = -9999.0; a2dVarIncRadF = -9999.0;  
        a2dVarWindF = -9999.0; a2dVarRelHumF = -9999.0;  a2dVarPaF = -9999.0; 
        a2dVarAlbedoF = -9999.0; a2dVarLAIF = -9999.0; 
        a2dVarSnowHeightF = -9999.9; a2dVarSnowKernelF = -9999.0; 
        a2dVarSnowCAF = -9999.0; a2dVarSnowQAF = -9999.0;
        
        a2dVarRainL = -9999.0; a2dVarTaL = -9999.0; a2dVarIncRadL = -9999.0;
        a2dVarWindL = -9999.0; a2dVarRelHumL = -9999.0; a2dVarPaL = -9999.0;
        a2dVarAlbedoL = -9999.0; a2dVarLAIL = -9999.0;
        a2dVarSnowHeightL = -9999.9; a2dVarSnowKernelL = -9999.0; 
        a2dVarSnowCAL = -9999.0; a2dVarSnowQAL = -9999.0;
        !------------------------------------------------------------------------------------------
                                                                                                
        !------------------------------------------------------------------------------------------
        ! Defining iRows and iCols (Land data)
        iRowsL = iRowsEndL - iRowsStartL + 1
        iColsL = iColsEndL - iColsStartL + 1
        ! Defining iRows and iCols (forcing data)
        iRowsF = iRowsEndF - iRowsStartF + 1
        iColsF = iColsEndF - iColsStartF + 1
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sPathData_Forcing = oHMC_Namelist(iID)%sPathData_Forcing_Gridded
        iFlagTypeData_Forcing = oHMC_Namelist(iID)%iFlagTypeData_Forcing_Gridded
        iScaleFactor = oHMC_Namelist(iID)%iScaleFactor
        
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Forcing Gridded ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Replace general path with specific time feature(s)
        call HMC_Tools_Generic_ReplaceText(sPathData_Forcing, '$yyyy', sTime(1:4))
        call HMC_Tools_Generic_ReplaceText(sPathData_Forcing, '$mm', sTime(6:7))
        call HMC_Tools_Generic_ReplaceText(sPathData_Forcing, '$dd', sTime(9:10))
        call HMC_Tools_Generic_ReplaceText(sPathData_Forcing, '$HH', sTime(12:13))
        call HMC_Tools_Generic_ReplaceText(sPathData_Forcing, '$MM', sTime(15:16))
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Check time step (iT)
        if (oHMC_Vars(iID)%iTime .lt. oHMC_Namelist(iID)%iNTime) then
 
            !------------------------------------------------------------------------------------------
            ! Subroutine for reading sequential netCDF forcing data 
            if (iFlagTypeData_Forcing == 2) then

                !------------------------------------------------------------------------------------------
                ! Call subroutine to get forcing data in netCDF format
#ifdef LIB_NC
                call HMC_Data_Forcing_Gridded_NC(iID, &
                                        sPathData_Forcing, &
                                        iRowsF, iColsF, &
                                        sTime, &
                                        a2dVarRainF, a2dVarTaF, a2dVarIncRadF, &
                                        a2dVarWindF, a2dVarRelHumF, a2dVarPaF, &
                                        a2dVarAlbedoF, a2dVarLAIF, &
                                        a2dVarSnowHeightF, a2dVarSnowKernelF, &
                                        a2dVarSnowCAF, a2dVarSnowQAF)
#else   
                ! Redefinition of forcing data flag (if netCDF library is not linked)
                iFlagTypeData_Forcing = 1 
                call mprintf(.true., iWARN, ' ATTENTION: '// &
                                            'forcing gridded data type selected was netCDF but library is not linked! '// &
                                            'Will be used data in binary format!')
#endif
                !------------------------------------------------------------------------------------------
                           
            endif
            !------------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------------
            ! Subroutine for reading sequential binary forcing data
            if (iFlagTypeData_Forcing == 1) then

                !------------------------------------------------------------------------------------------
                ! Calling subroutine to read data retrieved from drops db
                call HMC_Data_Forcing_Gridded_Binary(iID, &
                                            sPathData_Forcing, &
                                            iRowsF, iColsF, &
                                            sTime, &
                                            a2dVarRainF, a2dVarTaF, a2dVarIncRadF, &
                                            a2dVarWindF, a2dVarRelHumF, a2dVarPaF, &
                                            a2dVarAlbedoF, a2dVarLAIF, &
                                            a2dVarSnowHeightF, a2dVarSnowKernelF, &
                                            a2dVarSnowCAF, a2dVarSnowQAF, &
                                            iScaleFactor)
                !------------------------------------------------------------------------------------------

            endif
            !------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------
            ! Debug
            if (iDEBUG.gt.0) then
                call mprintf(.true., iINFO_Extra, ' ========= FORCING GRIDDED START =========== ')
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarRainF, oHMC_Vars(iID)%a2iMask, 'RAIN START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarTaF, oHMC_Vars(iID)%a2iMask, 'TA START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarIncRadF, oHMC_Vars(iID)%a2iMask, 'INCRAD START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarWindF, oHMC_Vars(iID)%a2iMask, 'WIND START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarRelHumF, oHMC_Vars(iID)%a2iMask, 'RELHUM START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarPaF, oHMC_Vars(iID)%a2iMask, 'PA START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarAlbedoF, oHMC_Vars(iID)%a2iMask, 'ALBEDO START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarLAIF, oHMC_Vars(iID)%a2iMask, 'LAI START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowHeightF, oHMC_Vars(iID)%a2iMask, 'SNOWHEIGHT START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowKernelF, oHMC_Vars(iID)%a2iMask, 'SNOWKERNEL START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowCAF, oHMC_Vars(iID)%a2iMask, 'SNOWCA START') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowQAF, oHMC_Vars(iID)%a2iMask, 'SNOWQA START') )
                call mprintf(.true., iINFO_Extra, '')
            endif
            !------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------
            ! Grid switcher land-forcing
            call HMC_Tools_Generic_SwitchGrid(oHMC_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarRainL, &
                                              iRowsF, iColsF, a2dVarRainF, &
                                              oHMC_Vars(iID)%a2dDem, &
                                              oHMC_Vars(iID)%a2iXIndex, oHMC_Vars(iID)%a2iYIndex)
            call HMC_Tools_Generic_SwitchGrid(oHMC_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarTaL, &
                                              iRowsF, iColsF, a2dVarTaF, &
                                              oHMC_Vars(iID)%a2dDem, &
                                              oHMC_Vars(iID)%a2iXIndex, oHMC_Vars(iID)%a2iYIndex)
            call HMC_Tools_Generic_SwitchGrid(oHMC_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarIncRadL, &
                                              iRowsF, iColsF, a2dVarIncRadF, &
                                              oHMC_Vars(iID)%a2dDem, &
                                              oHMC_Vars(iID)%a2iXIndex, oHMC_Vars(iID)%a2iYIndex)
            call HMC_Tools_Generic_SwitchGrid(oHMC_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarWindL, &
                                              iRowsF, iColsF, a2dVarWindF, &
                                              oHMC_Vars(iID)%a2dDem, &
                                              oHMC_Vars(iID)%a2iXIndex, oHMC_Vars(iID)%a2iYIndex)
            call HMC_Tools_Generic_SwitchGrid(oHMC_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarRelHumL, &
                                              iRowsF, iColsF, a2dVarRelHumF, &
                                              oHMC_Vars(iID)%a2dDem, &
                                              oHMC_Vars(iID)%a2iXIndex, oHMC_Vars(iID)%a2iYIndex)    
            call HMC_Tools_Generic_SwitchGrid(oHMC_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarPaL, &
                                              iRowsF, iColsF, a2dVarPaF, &
                                              oHMC_Vars(iID)%a2dDem, &
                                              oHMC_Vars(iID)%a2iXIndex, oHMC_Vars(iID)%a2iYIndex)                            
            call HMC_Tools_Generic_SwitchGrid(oHMC_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarAlbedoL, &
                                              iRowsF, iColsF, a2dVarAlbedoF, &
                                              oHMC_Vars(iID)%a2dDem, &
                                              oHMC_Vars(iID)%a2iXIndex, oHMC_Vars(iID)%a2iYIndex) 
            call HMC_Tools_Generic_SwitchGrid(oHMC_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarLAIL, &
                                              iRowsF, iColsF, a2dVarLAIF, &
                                              oHMC_Vars(iID)%a2dDem, &
                                              oHMC_Vars(iID)%a2iXIndex, oHMC_Vars(iID)%a2iYIndex)                     
            call HMC_Tools_Generic_SwitchGrid(oHMC_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarSnowHeightL, &
                                              iRowsF, iColsF, a2dVarSnowHeightF, &
                                              oHMC_Vars(iID)%a2dDem, &
                                              oHMC_Vars(iID)%a2iXIndex, oHMC_Vars(iID)%a2iYIndex)   
            call HMC_Tools_Generic_SwitchGrid(oHMC_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarSnowKernelL, &
                                              iRowsF, iColsF, a2dVarSnowKernelF, &
                                              oHMC_Vars(iID)%a2dDem, &
                                              oHMC_Vars(iID)%a2iXIndex, oHMC_Vars(iID)%a2iYIndex) 
            call HMC_Tools_Generic_SwitchGrid(oHMC_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarSnowCAL, &
                                              iRowsF, iColsF, a2dVarSnowCAF, &
                                              oHMC_Vars(iID)%a2dDem, &
                                              oHMC_Vars(iID)%a2iXIndex, oHMC_Vars(iID)%a2iYIndex)  
            call HMC_Tools_Generic_SwitchGrid(oHMC_Namelist(iID)%iFlagGrid, &
                                              iRowsL, iColsL, a2dVarSnowQAL, &
                                              iRowsF, iColsF, a2dVarSnowQAF, &
                                              oHMC_Vars(iID)%a2dDem, &
                                              oHMC_Vars(iID)%a2iXIndex, oHMC_Vars(iID)%a2iYIndex)  
            !------------------------------------------------------------------------------------------
                                            
            !------------------------------------------------------------------------------------------
            ! Check variable(s) limits and domain
            a2dVarRainL = check2Dvar(a2dVarRainL,               oHMC_Vars(iID)%a2iMask,     0.0,    850.0,  0.0)
            a2dVarTaL = check2Dvar(a2dVarTaL,                   oHMC_Vars(iID)%a2iMask,     -70.0,  60.0,   0.0 )    
            a2dVarIncRadL = check2Dvar(a2dVarIncRadL,           oHMC_Vars(iID)%a2iMask,     0.0,    1412.0, 0.0 ) 
            a2dVarWindL = check2Dvar(a2dVarWindL,               oHMC_Vars(iID)%a2iMask,     0.0,    80.0,   0.0 ) 
            a2dVarRelHumL = check2Dvar(a2dVarRelHumL,           oHMC_Vars(iID)%a2iMask,     0.0,    100.0,  0.0 )
            a2dVarPaL = check2Dvar(a2dVarPaL,                   oHMC_Vars(iID)%a2iMask,     50.0,   101.3,  101.3 )
            a2dVarSnowHeightL = check2Dvar(a2dVarSnowHeightL,   oHMC_Vars(iID)%a2iMask,     0.0,    1000.0, -9999.0 ) 
            a2dVarSnowKernelL = check2Dvar(a2dVarSnowKernelL,   oHMC_Vars(iID)%a2iMask,     0.0,    1.0,    -9999.0 )
            a2dVarSnowCAL = check2Dvar(a2dVarSnowCAL,           oHMC_Vars(iID)%a2iMask,     -1.0,   2.0,    -9999.0 )  
            a2dVarSnowQAL = check2Dvar(a2dVarSnowQAL,           oHMC_Vars(iID)%a2iMask,     0.0,    1.0 ,   -9999.0 ) 
            !------------------------------------------------------------------------------------------

        else
            !------------------------------------------------------------------------------------------
            ! Extra steps condition
            a2dVarRainL = 0.0;
            a2dVarTaL = oHMC_Vars(iID)%a2dTa; 
            a2dVarIncRadL = oHMC_Vars(iID)%a2dK;  
            a2dVarWindL = oHMC_Vars(iID)%a2dW; 
            a2dVarRelHumL = oHMC_Vars(iID)%a2dRHum
            a2dVarPaL = oHMC_Vars(iID)%a2dPres
            a2dVarAlbedoL = oHMC_Vars(iID)%a2dAlbedo; 
            a2dVarLAIL = oHMC_Vars(iID)%a2dLAI; 
            a2dVarSnowHeightL = oHMC_Vars(iID)%a2dSHeight; 
            a2dVarSnowKernelL = oHMC_Vars(iID)%a2dSKernel; 
            a2dVarSnowCAL = -9999.0; 
            a2dVarSnowQAL = -9999.0;
            
            ! Info message for extra time step(s)
            call mprintf(.true., iINFO_Extra, ' Extra time step ---> Forcing data are set constant to last real value')
            call mprintf(.true., iINFO_Extra, ' Extra time step ---> Rain data are set to 0.0')
            call mprintf(.true., iINFO_Extra, ' Extra time step ---> Snow data are set to -9999.0')
            !------------------------------------------------------------------------------------------
        endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Check update and save forcing data to local variable(s) to global workspace
        ! Rain
        if ( .not. all(a2dVarRainL.eq.-9999.0) ) then
            oHMC_Vars(iID)%a2dRain = a2dVarRainL
        else
            call mprintf(.true., iWARN, ' ATTENTION: all rain values are undefined! Check forcing data!' )
        endif
        ! Air temperature
        if ( .not. all(a2dVarTaL.eq.-9999.0) ) then
            oHMC_Vars(iID)%a2dTa = a2dVarTaL
        else
            call mprintf(.true., iWARN, ' ATTENTION: all air temperature values are undefined! Check forcing data!' )
        endif
        ! Incoming radiation
        if ( .not. all(a2dVarIncRadL.eq.-9999.0) ) then
            oHMC_Vars(iID)%a2dK = a2dVarIncRadL
        else
            call mprintf(.true., iWARN, ' ATTENTION: all incoming radiation values are undefined! Check forcing data!' )
        endif
        ! Wind
        if ( .not. all(a2dVarWindL.eq.-9999.0) ) then
            oHMC_Vars(iID)%a2dW = a2dVarWindL
        else
            call mprintf(.true., iWARN, ' ATTENTION: all wind values are undefined! Check forcing data!' )
        endif
        ! Relative humidity
        if ( .not. all(a2dVarRelHumL.eq.-9999.0) ) then
            oHMC_Vars(iID)%a2dRHum = a2dVarRelHumL
        else
            call mprintf(.true., iWARN, ' ATTENTION: all relative humidity values are undefined! Check forcing data!' )
        endif
        ! Air pressure
        if ( .not. all(a2dVarPaL.eq.101.3) ) then
            oHMC_Vars(iID)%a2dPres = a2dVarPaL
        else
            if ( all(a2dVarPaL.eq.101.3) ) then
                call mprintf(.true., iWARN, ' ATTENTION: all air pressure values are equal to 101.3 [kPa]! '// &
                                                         'AirPressure initializes using DEM information!')
                where (oHMC_Vars(iID)%a2dDem.gt.0.0)
                    a2dVarPaL = 101.3*((293-0.0065*oHMC_Vars(iID)%a2dDem)/293)**5.26 ![kPa]
                elsewhere
                    a2dVarPaL = 0.0
                endwhere
                oHMC_Vars(iID)%a2dPres = a2dVarPaL
            else
                call mprintf(.true., iWARN, ' ATTENTION: all air pressure values are undefined! Check forcing data!' )
            endif
        endif 
        ! LAI
        if ( .not. all(a2dVarLAIL.eq.-9999.0) ) then
            oHMC_Vars(iID)%a2dLAI = a2dVarLAIL
        else
            call mprintf(.true., iWARN, ' ATTENTION: all LAI values are undefined! Check forcing data!' )
        endif 
        ! Albedo
        if ( .not. all(a2dVarAlbedoL.eq.-9999.0) ) then
            oHMC_Vars(iID)%a2dAlbedo = a2dVarAlbedoL
        else
            call mprintf(.true., iWARN, ' ATTENTION: all albedo values are undefined! Check forcing data!' )
        endif 
        ! Snow Height
        if ( .not. all(a2dVarSnowHeightL.eq.-9999.0) ) then
            oHMC_Vars(iID)%a2dSHeight = a2dVarSnowHeightL
        else
            call mprintf(.true., iWARN, ' ATTENTION: all snow height values are undefined! Check forcing data!' )
        endif 
        ! Snow Kernel
        if ( .not. all(a2dVarSnowKernelL.eq.-9999.0) ) then
            oHMC_Vars(iID)%a2dSKernel = a2dVarSnowKernelL
        else
            call mprintf(.true., iWARN, ' ATTENTION: all snow kernel values are undefined! Check forcing data!' )
        endif 
        ! Snow Cover Area
        if ( .not. all(a2dVarSnowCAL.eq.-9999.0) ) then
            oHMC_Vars(iID)%a2dSCA = a2dVarSnowCAL
        else
            call mprintf(.true., iWARN, ' ATTENTION: all snow CA values are undefined! Check forcing data!' )
            oHMC_Vars(iID)%a2dSCA = -9999.0
        endif 
        ! Snow Quality Assessment
        if ( .not. all(a2dVarSnowQAL.eq.-9999.0) ) then
            oHMC_Vars(iID)%a2dSQA = a2dVarSnowQAL
        else
            call mprintf(.true., iWARN, ' ATTENTION: all snow QA values are undefined! Check forcing data!' )
            oHMC_Vars(iID)%a2dSQA = -9999.0
        endif 
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Debug
        if (iDEBUG.gt.0) then
            call mprintf(.true., iINFO_Extra, '')
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dRain, oHMC_Vars(iID)%a2iMask, 'RAIN END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dTa, oHMC_Vars(iID)%a2iMask, 'TA END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dK, oHMC_Vars(iID)%a2iMask, 'INCRAD END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dW, oHMC_Vars(iID)%a2iMask, 'WIND END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dRHum, oHMC_Vars(iID)%a2iMask, 'RELHUM END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dPres, oHMC_Vars(iID)%a2iMask, 'PA END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dAlbedo, oHMC_Vars(iID)%a2iMask, 'ALBEDO END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dLAI, oHMC_Vars(iID)%a2iMask, 'LAI END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dSHeight, oHMC_Vars(iID)%a2iMask, 'SNOWHEIGHT END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dSKernel, oHMC_Vars(iID)%a2iMask, 'SNOWKERNEL END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dSCA, oHMC_Vars(iID)%a2iMask, 'SNOWCA END') )
            call mprintf(.true., iINFO_Extra, checkvar(oHMC_Vars(iID)%a2dSQA, oHMC_Vars(iID)%a2iMask, 'SNOWQA END') )
            call mprintf(.true., iINFO_Extra, ' ========= FORCING GRIDDED END =========== ')
        endif
        
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Data :: Forcing Gridded ... OK' )
        !------------------------------------------------------------------------------------------
        
    end subroutine HMC_Data_Forcing_Gridded_Cpl
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to read NC data forcing
#ifdef LIB_NC
    subroutine HMC_Data_Forcing_Gridded_NC(iID,  &
                                  sPathData_Forcing, &
                                  iRows, iCols, sTime, &
                                  a2dVarRain, a2dVarTa, a2dVarIncRad, &
                                  a2dVarWind, a2dVarRelHum, a2dVarPa, &
                                  a2dVarAlbedo, a2dVarLAI, &
                                  a2dVarSnowHeight, a2dVarSnowKernel, &
                                  a2dVarSnowCA, a2dVarSnowQA)
                                  
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)                       :: iID                  
        
        character(len = 256), intent(in)        :: sPathData_Forcing
        character(len = 700)                    :: sFileNameData_Forcing, sFileNameData_Forcing_Zip
        character(len = 700)                    :: sCommandUnzipFile, sCommandRemoveFile
        character(len = 256)                    :: sVarName
        integer(kind = 4), intent(in)           :: iRows, iCols

        character(len = 19), intent(in)         :: sTime
        character(len = 12)                      :: sTimeMonth
        
        integer(kind = 4)                       :: iFlagSnow, iFlagSnowAssim
        real(kind = 4)                          :: dVarLAI, dVarAlbedo
        
        real(kind = 4), dimension(iCols, iRows)                 :: a2dVar
        
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarRain
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarTa
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarIncRad
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarWind      
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarRelHum       
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarPa 
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarAlbedo
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarLAI
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowHeight
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowKernel
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowCA
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowQA
       
        character(len = 256):: sVarUnits
        integer(kind = 4)   :: iErr
        integer(kind = 4)   :: iFileID
        
        logical             :: bFileExist
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarRain = -9999.0; a2dVarTa = -9999.0; a2dVarIncRad = -9999.0; a2dVarWind = -9999.0
        a2dVarRelHum = -9999.0; a2dVarPa = -9999.0; a2dVarAlbedo = -9999.0; a2dVarLAI = -9999.0
        a2dVarSnowHeight = -9999.0; a2dVarSnowKernel = -9999.0
        sFileNameData_Forcing = ''; sFileNameData_Forcing_Zip = ''; sTimeMonth = ''
        
        ! Checking date
        write(sTimeMonth,'(A,A,A)') sTime(1:4), sTime(6:7), sTime(9:10)
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sCommandUnzipFile = oHMC_Namelist(iID)%sCommandUnzipFile
        sCommandRemoveFile = oHMC_Namelist(iID)%sCommandRemoveFile
        iFlagSnow = oHMC_Namelist(iID)%iFlagSnow
        iFlagSnowAssim = oHMC_Namelist(iID)%iFlagSnowAssim
        
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Forcing Gridded :: NetCDF ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Filename forcing (example: hmc.dynamicdata.201404300000.nc.gz)
        sFileNameData_Forcing = trim(sPathData_Forcing)//"hmc.forcing-grid."// &
        sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
        sTime(12:13)//sTime(15:16)// &
        ".nc"

        ! Info netCDF filename
        call mprintf(.true., iINFO_Verbose, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing)//' ... ' )
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Checking file input availability
        sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
        inquire (file = trim(sFileNameData_Forcing)//'.gz', exist = bFileExist)
        if ( .not. bFileExist ) then
            !------------------------------------------------------------------------------------------
            ! Warning message
            call mprintf(.true., iWARN, ' No compressed forcing netCDF data found: '//trim(sFileNameData_Forcing_Zip) )
            ! Info netCDF filename
            call mprintf(.true., iINFO_Verbose, &
                         ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing)//' ... FAILED' )
            ! Info end
            call mprintf(.true., iINFO_Extra, ' Data :: Forcing Gridded :: NetCDF ... SKIPPED!' )
            !------------------------------------------------------------------------------------------
        else
            
            !------------------------------------------------------------------------------------------
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(oHMC_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Forcing_Zip, &
                                             sFileNameData_Forcing, .true.)
            !------------------------------------------------------------------------------------------
        
            !------------------------------------------------------------------------------------------
            ! Open netCDF file
            iErr = nf90_open(trim(sFileNameData_Forcing), NF90_NOWRITE, iFileID)
            if (iErr /= 0) then
                call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed netCDF file: '// &
                             trim(sFileNameData_Forcing)//' --> Undefined forcing data values' )
                call mprintf(.true., iINFO_Verbose, &
                             ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing)//' ... FAILED' )
            else
                
                !------------------------------------------------------------------------------------------
                ! RAIN
                sVarName = 'Rain'
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get rain data failed! ')
                    a2dVarRain = -9999.0;
                else
                    a2dVarRain = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! AIR TEMPERATURE
                sVarName = 'AirTemperature'
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get air temperature data failed! ')
                    a2dVarTa = -9999.0;
                else
                    a2dVarTa = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! INCOMING RADIATION
                sVarName = 'IncRadiation'
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get incoming radiation data failed! ')
                    a2dVarIncRad = -9999.0;
                else
                    a2dVarIncRad = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! RELATIVE HUMIDITY
                sVarName = 'RelHumidity'
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get relative humidity data failed! ')
                    a2dVarRelHum = -9999.0;
                else
                    a2dVarRelHum = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! WIND 
                sVarName = 'Wind'
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get wind data failed! ')
                    a2dVarWind = -9999.0;
                else
                    a2dVarWind = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! AIR PRESSURE
                sVarName = 'AirPressure'
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get air pressure data failed! Initialize with altitude function!')
                    a2dVarPa = -9999.0;
                    !a2dVarPa = 101.3*((293-0.0065*oHMC_Vars(iID)%a2dDem)/293)**5.26	![kPa]
                else
                    a2dVarPa = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! ALBEDO
                sVarName = 'Albedo'
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get albedo data failed! Initialize with monthly values!')
                    a2dVarAlbedo = -9999.0; dVarAlbedo = -9999.0
                    call HMC_Tools_Time_MonthVal(oHMC_Namelist(iID)%a1dAlbedoMonthly, sTimeMonth, dVarAlbedo)
                    a2dVarAlbedo = dVarAlbedo;
                else
                    a2dVarAlbedo = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! LAI
                sVarName = 'LAI'
                call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                if(iErr /= 0) then
                    call mprintf(.true., iWARN, ' ATTENTION: get LAI data failed! Initialize with monthly values!')
                    a2dVarLAI = -9999.0; dVarLAI = -9999.0
                    call HMC_Tools_Time_MonthVal(oHMC_Namelist(iID)%a1dLAIMonthly, sTimeMonth, dVarLAI)
                    a2dVarLAI = dVarLAI;
                else
                    a2dVarLAI = transpose(a2dVar)
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! SNOW VARIABLE(S)
                if ( (iFlagSnow.eq.1) .and. (iFlagSnowAssim.eq.1) ) then
                    
                    !------------------------------------------------------------------------------------------
                    ! SNOW HEIGHT
                    sVarName = 'SnowHeight'
                    call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                    if(iErr /= 0) then
                        call mprintf(.true., iWARN, ' ATTENTION: get snow height data failed! ')
                        a2dVarSnowHeight = -9999.0;
                    else
                        a2dVarSnowHeight = transpose(a2dVar)
                    endif
                    !------------------------------------------------------------------------------------------
                    
                    !------------------------------------------------------------------------------------------
                    ! SNOW KERNEL
                    sVarName = 'Kernel'
                    call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                    if(iErr /= 0) then
                        call mprintf(.true., iWARN, ' ATTENTION: get snow kernel data failed! ')
                        a2dVarSnowKernel = -9999.0;
                    else
                        a2dVarSnowKernel = transpose(a2dVar)
                    endif
                    !------------------------------------------------------------------------------------------
                    
                    !------------------------------------------------------------------------------------------
                    ! SNOW COVER AREA
                    sVarName = 'SCA'
                    call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                    if(iErr /= 0) then
                        call mprintf(.true., iWARN, ' ATTENTION: get snow CA data failed! ')
                        a2dVarSnowCA = -9999.0;
                    else
                        a2dVarSnowCA = transpose(a2dVar)
                    endif
                    !------------------------------------------------------------------------------------------

                    !------------------------------------------------------------------------------------------
                    ! SNOW QUALITY ASSESSMENT
                    sVarName = 'SQA'
                    call HMC_Tools_IO_Get2d_NC((sVarName), iFileID, a2dVar, sVarUnits, iCols, iRows, .false., iErr)
                    if(iErr /= 0) then
                        call mprintf(.true., iWARN, ' ATTENTION: get snow QA data failed! ')
                        a2dVarSnowQA = -9999.0;
                    else
                        a2dVarSnowQA = transpose(a2dVar)
                    endif
                    !------------------------------------------------------------------------------------------
                    
                elseif ( (iFlagSnow.eq.1) .and. (iFlagSnowAssim.eq.0) ) then
                    
                    !------------------------------------------------------------------------------------------
                    ! Exit without snow assimilation
                    a2dVarSnowKernel = -9999.0
                    a2dVarSnowHeight = -9999.0
                    a2dVarSnowCA = -9999.0
                    a2dVarSnowQA = -9999.0
                    !------------------------------------------------------------------------------------------
                    
                else
                    
                    !------------------------------------------------------------------------------------------
                    ! Exit without snow physics
                    a2dVarSnowKernel = -9999.0
                    a2dVarSnowHeight = -9999.0
                    a2dVarSnowCA = -9999.0
                    a2dVarSnowQA = -9999.0
                    !------------------------------------------------------------------------------------------
                    
                endif
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! Closing netCDF file
                iErr = nf90_close(iFileID)
                ! Remove uncompressed file (to save space on disk)
                call HMC_Tools_Generic_RemoveFile(oHMC_Namelist(iID)%sCommandRemoveFile, &
                                                  sFileNameData_Forcing, .false.)
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! Info netCDF filename
                call mprintf(.true., iINFO_Verbose, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing)//' ... OK' )
                ! Info end
                call mprintf(.true., iINFO_Extra, ' Data :: Forcing Gridded :: NetCDF ... OK' )
                !------------------------------------------------------------------------------------------
                
            endif
            !------------------------------------------------------------------------------------------
            
        endif
        !------------------------------------------------------------------------------------------
        
    end subroutine HMC_Data_Forcing_Gridded_NC
#endif
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to read binary forcing data
    subroutine HMC_Data_Forcing_Gridded_Binary(iID, &
                                      sPathData_Forcing, &
                                      iRows, iCols, sTime, &
                                      a2dVarRain, a2dVarTa, a2dVarIncRad, &
                                      a2dVarWind, a2dVarRelHum, a2dVarPa, &
                                      a2dVarAlbedo, a2dVarLAI, &
                                      a2dVarSnowHeight, a2dVarSnowKernel, &
                                      a2dVarSnowCA, a2dVarSnowQA, &
                                      iScaleFactor)
    
        !------------------------------------------------------------------------------------------
        ! Variable(s)
        integer(kind = 4)                   :: iID
                                      
        character(len = 256), intent(in)    :: sPathData_Forcing
        character(len = 700)                :: sFileNameData_Forcing, sFileNameData_Forcing_Zip
        character(len = 700)                :: sCommandUnzipFile
        character(len = 256)                :: sVarName
        integer(kind = 4), intent(in)       :: iRows, iCols
        real(kind = 4)                      :: dVar
        
        integer(kind = 4)                   :: iFlagSnow, iFlagSnowAssim
        
        character(len = 19), intent(in)     :: sTime
        character(len = 12)                 :: sTimeMonth
        
        real(kind = 4), dimension(iRows, iCols)                 :: a2dVar

        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarRain
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarTa
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarIncRad
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarWind      
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarRelHum       
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarPa 
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarAlbedo
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarLAI
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowHeight
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowKernel
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowCA
        real(kind = 4), dimension(iRows, iCols), intent(out)    :: a2dVarSnowQA
       
        character(len = 256):: sVarUnits
        integer(kind = 4)   :: iErr
        integer(kind = 4)   :: iFileID, iScaleFactor
        
        logical             :: bFileExist
        !------------------------------------------------------------------------------------------
	
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2dVarRain = -9999.0; a2dVarTa = -9999.0; a2dVarIncRad = -9999.0; a2dVarWind = -9999.0
        a2dVarRelHum = -9999.0; a2dVarPa = -9999.0; a2dVarAlbedo = -9999.0; a2dVarLAI = -9999.0
        a2dVarSnowKernel = -9999.0; a2dVarSnowHeight = -9999.0; 
        a2dVarSnowCA = -9999.0; a2dVarSnowCA = -9999.0 
        sFileNameData_Forcing = ''; sFileNameData_Forcing_Zip = ''; sTimeMonth = ''
        
        ! Checking date
        write(sTimeMonth,'(A,A,A)') sTime(1:4), sTime(6:7), sTime(9:10)
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Get global information
        sCommandUnzipFile = oHMC_Namelist(iID)%sCommandUnzipFile
        iFlagSnow = oHMC_Namelist(iID)%iFlagSnow
        iFlagSnowAssim = oHMC_Namelist(iID)%iFlagSnowAssim
        
        ! Info start
        call mprintf(.true., iINFO_Extra, ' Data :: Forcing Gridded :: Binary ... ' )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info binary file(s) time step
        call mprintf(.true., iINFO_Verbose, ' Get (forcing gridded) at time '//trim(sTime)//' ... ')
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Rain  (example: rain_201405010000.bin.gz)
        sFileNameData_Forcing = trim(sPathData_Forcing)//"Rain_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"  
        call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

        ! Checking file input availability
        sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
        inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values' )
            a2dVar = -9999.0
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(oHMC_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Forcing_Zip, &
                                             sFileNameData_Forcing, .true.)
            ! Read binary data
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Forcing, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            ! Remove uncompressed file (to save space on disk)
            call HMC_Tools_Generic_RemoveFile(oHMC_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Forcing, .false.)
            
        endif
        a2dVarRain = a2dVar
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Temperature  (example: temperature_201405010000.bin.gz)
        sFileNameData_Forcing = trim(sPathData_Forcing)//"Temperature_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

        ! Checking file input availability
        sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
        inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values' )
            a2dVar = -9999.0
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(oHMC_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Forcing_Zip, &
                                             sFileNameData_Forcing, .true.)
            ! Read binary data
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Forcing, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            ! Remove uncompressed file (to save space on disk)
            call HMC_Tools_Generic_RemoveFile(oHMC_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Forcing, .false.)
        endif
        a2dVarTa = a2dVar
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Radiation  (example: radiation_201405010000.bin.gz)
        sFileNameData_Forcing = trim(sPathData_Forcing)//"Radiation_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

        ! Checking file input availability
        sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
        inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values' )
            a2dVar = -9999.0
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(oHMC_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Forcing_Zip, &
                                             sFileNameData_Forcing, .true.)
            ! Read binary data
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Forcing, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            ! Remove uncompressed file (to save space on disk)
            call HMC_Tools_Generic_RemoveFile(oHMC_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Forcing, .false.)
        endif
        a2dVarIncRad = a2dVar
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Wind  (example: wind_201405010000.bin.gz)
        sFileNameData_Forcing = trim(sPathData_Forcing)//"Wind_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

        ! Checking file input availability
        sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
        inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values' )
            a2dVar = -9999.0
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(oHMC_Namelist(iID)%sCommandUnzipFile, & 
                                             sFileNameData_Forcing_Zip, &
                                             sFileNameData_Forcing, .true.)
            ! Read binary data
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Forcing, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            ! Remove uncompressed file (to save space on disk)
            call HMC_Tools_Generic_RemoveFile(oHMC_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Forcing, .false.)
        endif
        a2dVarWind = a2dVar
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! RelHum  (example: wind_201405010000.bin.gz)
        sFileNameData_Forcing = trim(sPathData_Forcing)//"RelUmid_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

        ! Checking file input availability
        sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
        inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values' )
            a2dVar = -9999.0
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(oHMC_Namelist(iID)%sCommandUnzipFile, & 
                                             sFileNameData_Forcing_Zip, & 
                                             sFileNameData_Forcing, .true.)
            ! Read binary data
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Forcing, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            ! Remove uncompressed file (to save space on disk)
            call HMC_Tools_Generic_RemoveFile(oHMC_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Forcing, .false.)
        endif
        a2dVarRelHum = a2dVar
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! AirPressure (example: pressure_201405010000.bin.gz)
        sFileNameData_Forcing = trim(sPathData_Forcing)//"Pressure_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

        ! Checking file input availability
        sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
        inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values' )
            a2dVar = -9999.0
            !a2dVar = 101.3*((293-0.0065*oHMC_Vars(iID)%a2dDem)/293)**5.26	![kPa]
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(oHMC_Namelist(iID)%sCommandUnzipFile, & 
                                             sFileNameData_Forcing_Zip, &
                                             sFileNameData_Forcing, .true.)
            ! Read binary data
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Forcing, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            ! Remove uncompressed file (to save space on disk)
            call HMC_Tools_Generic_RemoveFile(oHMC_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Forcing, .false.)
        endif
        a2dVarPa = a2dVar
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Albedo (example: albedo_201405010000.bin.gz)
        sFileNameData_Forcing = trim(sPathData_Forcing)//"Albedo_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

        ! Checking file input availability
        sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
        inquire (file = trim(sFileNameData_Forcing)//'.gz', exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values' )
            a2dVar = -9999.0; dVar = -9999.0;
            call HMC_Tools_Time_MonthVal(oHMC_Namelist(iID)%a1dAlbedoMonthly, sTimeMonth, dVar)
            a2dVar = dVar;
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(oHMC_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Forcing_Zip, &
                                             sFileNameData_Forcing, .true.)
            ! Read binary data
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Forcing, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            ! Remove uncompressed file (to save space on disk)
            call HMC_Tools_Generic_RemoveFile(oHMC_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Forcing, .false.)
        endif
        a2dVarAlbedo = a2dVar
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Albedo (example: lai_201405010000.bin.gz)
        sFileNameData_Forcing = trim(sPathData_Forcing)//"LAI_"// &
            sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
            sTime(12:13)//sTime(15:16)// &
            ".bin"
        call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

        ! Checking file input availability
        sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
        inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
        if ( .not. bFileExist ) then
            call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                         trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values' )
            a2dVar = -9999.0; dVar = -9999.0
            call HMC_Tools_Time_MonthVal(oHMC_Namelist(iID)%a1dLAIMonthly, sTimeMonth, dVar)
            a2dVar = dVar;
        else
            ! Unzip file
            call HMC_Tools_Generic_UnzipFile(oHMC_Namelist(iID)%sCommandUnzipFile, &
                                             sFileNameData_Forcing_Zip, &
                                             sFileNameData_Forcing, .true.)
            ! Read binary data
            call HMC_Tools_IO_Get2d_Binary(sFileNameData_Forcing, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
            ! Remove uncompressed file (to save space on disk)
            call HMC_Tools_Generic_RemoveFile(oHMC_Namelist(iID)%sCommandRemoveFile, &
                                              sFileNameData_Forcing, .false.)
        endif
        a2dVarLAI = a2dVar
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Snow variable(s)
        if ( (iFlagSnow.eq.1) .and. (iFlagSnowAssim.eq.1) ) then 
            
            !------------------------------------------------------------------------------------------
            ! SnowHeight (example: snowheight_201405010000.bin.gz)
            sFileNameData_Forcing = trim(sPathData_Forcing)//"SnowHeight_"// &
                sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
                sTime(12:13)//sTime(15:16)// &
                ".bin"
            call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

            ! Checking file input availability
            sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
            inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                             trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values' )
                a2dVar = -9999.0; 
            else
                ! Unzip file
                call HMC_Tools_Generic_UnzipFile(oHMC_Namelist(iID)%sCommandUnzipFile, &
                                                 sFileNameData_Forcing_Zip, &
                                                 sFileNameData_Forcing, .true.)
                ! Read binary data
                call HMC_Tools_IO_Get2d_Binary(sFileNameData_Forcing, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
                ! Remove uncompressed file (to save space on disk)
                call HMC_Tools_Generic_RemoveFile(oHMC_Namelist(iID)%sCommandRemoveFile, &
                                                  sFileNameData_Forcing, .false.)
            endif
            a2dVarSnowHeight = a2dVar
            !------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------
            ! Snow Kernel (example: kernel_201405010000.bin.gz)
            sFileNameData_Forcing = trim(sPathData_Forcing)//"Kernel_"// &
                sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
                sTime(12:13)//sTime(15:16)// &
                ".bin"
            call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

            ! Checking file input availability
            sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
            inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                             trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values' )
                a2dVar = -9999.0; 
            else
                ! Unzip file
                call HMC_Tools_Generic_UnzipFile(oHMC_Namelist(iID)%sCommandUnzipFile, &
                                                 sFileNameData_Forcing_Zip, &
                                                 sFileNameData_Forcing, .true.)
                ! Read binary data
                call HMC_Tools_IO_Get2d_Binary(sFileNameData_Forcing, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
                ! Remove uncompressed file (to save space on disk)
                call HMC_Tools_Generic_RemoveFile(oHMC_Namelist(iID)%sCommandRemoveFile, &
                                                  sFileNameData_Forcing, .false.)
            endif
            a2dVarSnowKernel = a2dVar
            !------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------
            ! Snow Cover Area (example: SCA_201405010000.bin.gz)
            sFileNameData_Forcing = trim(sPathData_Forcing)//"SCA_"// &
                sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
                sTime(12:13)//sTime(15:16)// &
                ".bin"
            call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

            ! Checking file input availability
            sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
            inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                             trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values' )
                a2dVar = -9999.0; 
            else
                ! Unzip file
                call HMC_Tools_Generic_UnzipFile(oHMC_Namelist(iID)%sCommandUnzipFile, &
                                                 sFileNameData_Forcing_Zip, &
                                                 sFileNameData_Forcing, .true.)
                ! Read binary data
                call HMC_Tools_IO_Get2d_Binary(sFileNameData_Forcing, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
                ! Remove uncompressed file (to save space on disk)
                call HMC_Tools_Generic_RemoveFile(oHMC_Namelist(iID)%sCommandRemoveFile, &
                                                  sFileNameData_Forcing, .false.)
            endif
            a2dVarSnowCA = a2dVar
            !------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------
            ! Snow Quality Assessment (example: SQA_201405010000.bin.gz)
            sFileNameData_Forcing = trim(sPathData_Forcing)//"SQA_"// &
                sTime(1:4)//sTime(6:7)//sTime(9:10)// & 
                sTime(12:13)//sTime(15:16)// &
                ".bin"
            call mprintf(.true., iINFO_Extra, ' Get filename (forcing gridded): '//trim(sFileNameData_Forcing) )

            ! Checking file input availability
            sFileNameData_Forcing_Zip = trim(sFileNameData_Forcing)//'.gz'
            inquire (file = sFileNameData_Forcing_Zip, exist = bFileExist)
            if ( .not. bFileExist ) then
                call mprintf(.true., iWARN, ' ATTENTION:  Problem opening uncompressed binary file: '// &
                             trim(sFileNameData_Forcing_Zip)//' --> Undefined forcing data values' )
                a2dVar = -9999.0; 
            else
                ! Unzip file
                call HMC_Tools_Generic_UnzipFile(oHMC_Namelist(iID)%sCommandUnzipFile, &
                                                 sFileNameData_Forcing_Zip, &
                                                 sFileNameData_Forcing, .true.)
                ! Read binary data
                call HMC_Tools_IO_Get2d_Binary(sFileNameData_Forcing, a2dVar, iRows, iCols, iScaleFactor, .true., iErr) 
                ! Remove uncompressed file (to save space on disk)
                call HMC_Tools_Generic_RemoveFile(oHMC_Namelist(iID)%sCommandRemoveFile, &
                                                  sFileNameData_Forcing, .false.)
            endif
            a2dVarSnowQA = a2dVar
            !------------------------------------------------------------------------------------------
            
        elseif ( (iFlagSnow.eq.1) .and. (iFlagSnowAssim.eq.0) ) then
            
            !------------------------------------------------------------------------------------------
            ! Exit without snow assimilation
            a2dVarSnowKernel = -9999.0
            a2dVarSnowHeight = -9999.0
            a2dVarSnowCA = -9999.0
            a2dVarSnowQA = -9999.0
            !------------------------------------------------------------------------------------------
            
        else
            
            !------------------------------------------------------------------------------------------
            ! Exit without snow physics
            a2dVarSnowKernel = -9999.0
            a2dVarSnowHeight = -9999.0
            a2dVarSnowCA = -9999.0
            a2dVarSnowQA = -9999.0
            !------------------------------------------------------------------------------------------
            
        endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info binary file(s) time step
        call mprintf(.true., iINFO_Verbose, ' Get (forcing gridded) at time '//trim(sTime)//' ... OK')
        ! Info end
        call mprintf(.true., iINFO_Extra, ' Data :: Forcing Gridded :: Binary ... OK' )
        !------------------------------------------------------------------------------------------
        
        
    end subroutine HMC_Data_Forcing_Gridded_Binary
    !------------------------------------------------------------------------------------------
    
end module HMC_Module_Data_Forcing_Gridded
!-----------------------------------------------------------------------------------------