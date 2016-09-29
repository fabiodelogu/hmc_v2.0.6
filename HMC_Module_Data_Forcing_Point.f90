!------------------------------------------------------------------------------------------    
! File:   HMC_Module_Data_Forcing_Point.f90
! Author(s): Fabio Delogu, Francesco Silvestro, Simone Gabellani
!
! Created on April 22, 2015, 5:19 PM
!------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------
! Module Header
module HMC_Module_Data_Forcing_Point
    
    !------------------------------------------------------------------------------------------
    ! External module(s) for all subroutine in this module
    use HMC_Module_Namelist,                    only:   oHMC_Namelist
    use HMC_Module_Vars_Loader,                 only:   oHMC_Vars
    
    use HMC_Module_Tools_Debug
    use HMC_Module_Tools_Generic,               only:   HMC_Tools_Generic_SmoothTimeSeries
    
    use HMC_Module_Tools_IO,                    only:   HMC_Tools_IO_Get2d_ASCII
                             
    ! Implicit none for all subroutines in this module
    implicit none
    !------------------------------------------------------------------------------------------
    
contains
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to manage forcing point data
    subroutine HMC_Data_Forcing_Point_Cpl(iID, sTime, &
                                          iNData, iETime, &  
                                          iNLake, iNDam, iNPlant, iNJoint, iNCatch, iNRelease)
        
        !------------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4)           :: iID
        character(len = 19)         :: sTime
        
        integer(kind = 4)           :: iNData, iETime
        integer(kind = 4)           :: iNLake, iNDam, iNPlant, iNJoint, iNCatch, iNRelease
        
        integer(kind = 4)           :: iFlagTypeData_Forcing
        character(len = 256)        :: sPathData_Forcing
        !------------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------------
        ! Check logical variable to control forcing point section
        if (.not. oHMC_Vars(iID)%bLogForcingPoint) then
        
            !------------------------------------------------------------------------------------------
            ! Get global information
            sPathData_Forcing = oHMC_Namelist(iID)%sPathData_Forcing_Point
            iFlagTypeData_Forcing = oHMC_Namelist(iID)%iFlagTypeData_Forcing_Point
            
            ! Info start
            call mprintf(.true., iINFO_Extra, ' Data :: Forcing Point ... ' )
            !------------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------------
            ! Subroutine for reading sequential ASCII forcing data point
            if (iFlagTypeData_Forcing == 1) then

                !------------------------------------------------------------------------------------------
                ! Choosing data type
                call mprintf(.true., iINFO_Extra, ' Data :: Forcing Point :: ASCII ... ' )
                !------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------
                ! Get point time-series for plant(s)
                call HMC_Data_Forcing_Point_Plant(iID, iNPlant, iNData, iETime)
                !------------------------------------------------------------------------------------------
                        
                !------------------------------------------------------------------------------------------
                ! Get point time-series for intake(s)
                call HMC_Data_Forcing_Point_Intake(iID, iNCatch, iNRelease, iNData, iETime)
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! Info end
                call mprintf(.true., iINFO_Extra, ' Data :: Forcing Point :: ASCII ... OK' )
                !------------------------------------------------------------------------------------------
                
            endif
            !------------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------------
            ! Subroutine for reading sequential unknown forcing data point
            if (iFlagTypeData_Forcing == 2) then

                !------------------------------------------------------------------------------------------
                ! Choosing data type
                call mprintf(.true., iERROR, ' Using UNKNOWN data point forcing. Check settings file!')
                !------------------------------------------------------------------------------------------

            endif
            !------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------
            ! Update logical variable to control forcing point section
            oHMC_Vars(iID)%bLogForcingPoint = .true.
            ! Info end
            call mprintf(.true., iINFO_Extra, ' Data :: Forcing Point ... OK' )
            !------------------------------------------------------------------------------------------
            
        else
            !------------------------------------------------------------------------------------------
            ! Info forcing point data
            call mprintf(.true., iINFO_Extra, ' Forcing point data loaded previously! Skipping this step!')
            ! Info end
            call mprintf(.true., iINFO_Extra, ' Data :: Forcing Point ... SKIPPED!' )
            !------------------------------------------------------------------------------------------
        endif
        !------------------------------------------------------------------------------------------

    end subroutine HMC_Data_Forcing_Point_Cpl
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to get time series of plant(s)
    subroutine HMC_Data_Forcing_Point_Plant(iID, iNPlant, iNData, iETime)

        !------------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4)           :: iID
        integer(kind = 4)           :: iP, iI, iJ, iErr
        integer(kind = 4)           :: iNPlant, iNData, iDtDataForcing, iETime
        real(kind = 4)              :: dAreaCell
        
        character(len = 256)        :: sPathData
        character(len = 256)        :: sNamePlant
        character(len = 256)        :: sFilePlant
        
        character(len = 256), dimension(iETime)                     :: a1sVar 
        real(kind = 4), dimension(iNData)                           :: a1dVar, a1dVarHydro
        real(kind = 4), dimension(iNPlant, iETime)                  :: a2dVarHydro

        logical                     :: bFileExist
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a1sVar = ''; a1dVar = 0.0; a1dVarHydro = 0.0
        a2dVarHydro = 0.0; 
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info
        call mprintf(.true., iINFO_Verbose, ' Get plant forcing filename')
        
        ! Debug
        call mprintf(.true., iINFO_Extra, checkarray(a2dVarHydro(:,2), 'HYDRO PLANT START') )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Check plant(s) availability
        if (iNPlant .gt. 0) then
        
            !------------------------------------------------------------------------------------------
            ! Get global information
            sPathData = oHMC_Namelist(iID)%sPathData_Forcing_Point
            iDtDataForcing = oHMC_Namelist(iID)%iDtData_Forcing
            !------------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------------
            ! Cycle on plant(s)
            do iP = 1, iNPlant

                !------------------------------------------------------------------------------------------
                ! Initialize step variable(s)
                a1dVar = 0.0; iI = 0; iJ = 0; dAreaCell = 0.0
                !------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------
                ! Get global information
                iI = oHMC_Vars(iID)%a2iXYPlant(iP,2); iJ = oHMC_Vars(iID)%a2iXYPlant(iP,1);
                dAreaCell = oHMC_Vars(iID)%a2dAreaCell(iI, iJ)
                !------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------
                ! Plant name and filename data
                sNamePlant = oHMC_Vars(iID)%a1sNamePlant(iP)
                sFilePlant = trim(sPathData)//'hmc.forcing-point.plant_'//trim(sNamePlant)//'.txt'
                !------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------
                ! Read ASCII file
                call HMC_Tools_IO_Get2d_ASCII(sFilePlant, a1sVar, a1dVar, iNData, .false. , iErr)
                !------------------------------------------------------------------------------------------
                
                !------------------------------------------------------------------------------------------
                ! Check data
                if (.not. all(a1dVar.eq.-9999.0) ) then

                    !------------------------------------------------------------------------------------------
                    ! Control data 
                    where (a1dVar .gt. 0.0)
                        a1dVarHydro = a1dVar*1000*real(iDtDataForcing)/dAreaCell
                    elsewhere
                        a1dVarHydro = 0.0
                    endwhere

                    ! Smooth data
                    call HMC_Tools_Generic_SmoothTimeSeries(a1dVarHydro, iNData, 2)
                    !------------------------------------------------------------------------------------------

                    !------------------------------------------------------------------------------------------
                    ! Save all information in one variable
                    a2dVarHydro(iP, :) = a1dVarHydro
                    !------------------------------------------------------------------------------------------
                else
                    !------------------------------------------------------------------------------------------
                    ! No data available
                    call mprintf(.true., iWARN, ' ATTENTION: all plant(s) values are undefined! Check forcing data!' )
                    a2dVarHydro(iP, :) = -9999.0 ! NoData dam value == -9999.0 NOT CHANGE
                    !------------------------------------------------------------------------------------------
                endif
                !------------------------------------------------------------------------------------------

            enddo 
            !------------------------------------------------------------------------------------------
        
        else
            !------------------------------------------------------------------------------------------
            ! No plant(s)
            a2dVarHydro = -9999.0 ! NoData dam value == -9999.0 NOT CHANGE
            !------------------------------------------------------------------------------------------
        endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Debug
        call mprintf(.true., iINFO_Extra, checkarray(a2dVarHydro(:,2), 'HYDRO PLANT END') )
        
        ! Pass local variable(s) to global workspace
        oHMC_Vars(iID)%a2dHydroPlant = a2dVarHydro
        !------------------------------------------------------------------------------------------
        
    end subroutine HMC_Data_Forcing_Point_Plant
    !------------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to get time series of intake(s)
    subroutine HMC_Data_Forcing_Point_Intake(iID, iNCatch, iNRelease, iNData, iETime)
        
        !------------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4)           :: iID
        integer(kind = 4)           :: iR, iP, iStep, iI, iJ, iErr, iTimeShift, iETime
        integer(kind = 4)           :: iNCatch, iNRelease, iNData, iDtDataForcing
        real(kind = 4)              :: dAreaCell, dTCorrCatch, dWeigthCatch
        
        character(len = 256)        :: sPathData
        character(len = 256)        :: sNameRelease, sNameCatch
        character(len = 256)        :: sFileRelease
        
        character(len = 256), dimension(iNData)                     :: a1sVar 
        real(kind = 4), dimension(iNData)                           :: a1dVar
        real(kind = 4), dimension(iNData)                           :: a1dVarHydroO
        real(kind = 4), dimension(iNData)                           :: a1dVarHydroR
        real(kind = 4), dimension(iNData)                           :: a1dVarHydroC
        
        real(kind = 4), dimension(iNRelease, iETime)                :: a2dVarHydroR
        real(kind = 4), dimension(iNCatch, iETime)                  :: a2dVarHydroC
   
        logical                     :: bFileExist
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a1sVar = ''; a1dVar = 0.0; 
        a1dVarHydroO = 0.0; a1dVarHydroR = 0.0; a1dVarHydroC = 0.0;
        a2dVarHydroR = 0.0; a2dVarHydroC = 0.0; 
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Info
        call mprintf(.true., iINFO_Verbose, ' Get intake (releases and catches) forcing filename')
        ! Debug
        call mprintf(.true., iINFO_Extra, checkarray(a2dVarHydroC(:,2), 'HYDRO CATCH START') )
        call mprintf(.true., iINFO_Extra, checkarray(a2dVarHydroR(:,2), 'HYDRO RELEASE START') )
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Check release(s) availability
        if (iNRelease .gt. 0) then
        
            !------------------------------------------------------------------------------------------
            ! Get global information
            sPathData = oHMC_Namelist(iID)%sPathData_Forcing_Point
            iDtDataForcing = oHMC_Namelist(iID)%iDtData_Forcing
            !------------------------------------------------------------------------------------------

            !------------------------------------------------------------------------------------------
            ! Cycle on release(s)
            do iR = 1, iNRelease

                !------------------------------------------------------------------------------------------
                ! Initialize step variable(s)
                a1dVar = -9999.0; iI = 0; iJ = 0; dAreaCell = 0.0
                !------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------
                ! Get global information
                iI = oHMC_Vars(iID)%a2iXYRelease(iR,2); iJ = oHMC_Vars(iID)%a2iXYRelease(iR,1);
                dAreaCell = oHMC_Vars(iID)%a2dAreaCell(iI, iJ)
                !------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------
                ! Plant name and filename data
                sNameRelease = oHMC_Vars(iID)%a1sNameRelease(iR)
                sFileRelease = trim(sPathData)//'hmc.forcing-point.plant_'//trim(sNameRelease)//'.txt'
                !------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------
                ! Read ASCII file
                call HMC_Tools_IO_Get2d_ASCII(sFileRelease, a1sVar, a1dVar, iNData, .false. , iErr)
                !------------------------------------------------------------------------------------------

                !------------------------------------------------------------------------------------------
                ! Check data
                if (.not. all(a1dVar.eq.-9999.)) then

                    !------------------------------------------------------------------------------------------
                    ! Control data 
                    where (a1dVar .gt. 0.0)
                        a1dVarHydroO = a1dVar*1.0
                        a1dVarHydroR = a1dVar*1000*real(iDtDataForcing)/dAreaCell
                    elsewhere
                        a1dVarHydroO = 0.0
                        a1dVarHydroR = 0.0
                    endwhere

                    ! Smooth data
                    call HMC_Tools_Generic_SmoothTimeSeries(a1dVarHydroO, iNData, 2)
                    call HMC_Tools_Generic_SmoothTimeSeries(a1dVarHydroR, iNData, 2)
                    !------------------------------------------------------------------------------------------

                    !------------------------------------------------------------------------------------------
                    ! Cycle on catch(es)
                    do iP = 1, iNCatch

                        !------------------------------------------------------------------------------------------
                        ! Get release and catch name
                        sNameRelease = oHMC_Vars(iID)%a1sNameRelease(iR)
                        sNameCatch = oHMC_Vars(iID)%a1sNameCatch(iP)
                        !------------------------------------------------------------------------------------------

                        !------------------------------------------------------------------------------------------
                        ! Check release and catch name(s9
                        if ( trim(sNameRelease) .eq. trim(sNameCatch) ) then

                            !------------------------------------------------------------------------------------------
                            ! Get info catch
                            dTCorrCatch = oHMC_Vars(iID)%a1dTCorrCatch(iP)
                            dWeigthCatch = oHMC_Vars(iID)%a1dWeigthCatch(iP)

                            iTimeShift = int(dTCorrCatch*60/iDtDataForcing)
                            !------------------------------------------------------------------------------------------

                            !------------------------------------------------------------------------------------------
                            ! Compute hydro arrays for catch(es)
                            do iStep = 1, iNData - 1
                                if ( (iStep - iTimeShift).gt.1 ) then
                                    a1dVarHydroC(iStep - iTimeShift) = a1dVarHydroO(iStep)*dWeigthCatch ! [m^3/s]
                                endif
                            enddo
                            do iStep = iNData - 1 - iTimeShift, iNData -1
                                a1dVarHydroC(iStep) = a1dVarHydroO(iNData-1)*dWeigthCatch
                            enddo

                            a2dVarHydroC(iP, :) = a1dVarHydroC
                            !------------------------------------------------------------------------------------------
                        else
                            !------------------------------------------------------------------------------------------
                            ! No data available
                            call mprintf(.true., iWARN, ' ATTENTION: release and catch name(s) link not available!' )
                            a2dVarHydroC(iP, :) = 0.0
                            !------------------------------------------------------------------------------------------
                        endif
                        !------------------------------------------------------------------------------------------

                    enddo
                    !------------------------------------------------------------------------------------------

                    !------------------------------------------------------------------------------------------
                    ! Save all information in one variable
                    a2dVarHydroR(iR, :) = a1dVarHydroR
                    !------------------------------------------------------------------------------------------
                else
                    !------------------------------------------------------------------------------------------
                    ! No data available
                    call mprintf(.true., iWARN, ' ATTENTION: all catch(es) values are undefined! Check forcing data!' )
                    a2dVarHydroR(iR, :) = 0.0
                    !------------------------------------------------------------------------------------------
                endif
                !------------------------------------------------------------------------------------------

            enddo 
            !------------------------------------------------------------------------------------------
        
        else
            !------------------------------------------------------------------------------------------
            ! No release(s) and no catch(es)
            a2dVarHydroR = 0.0
            a2dVarHydroC = 0.0
            !------------------------------------------------------------------------------------------
        endif
        !------------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------------
        ! Debug
        !call mprintf(.true., iINFO_Extra, checkarray(a2dVarHydroC(:,2), 'HYDRO CATCH END') )
        !call mprintf(.true., iINFO_Extra, checkarray(a2dVarHydroR(:,2), 'HYDRO RELEASE END') )
        
        ! Pass local variable(s) to global workspace
        oHMC_Vars(iID)%a2dHydroRelease = a2dVarHydroR
        oHMC_Vars(iID)%a2dHydroCatch = a2dVarHydroC
        !------------------------------------------------------------------------------------------
        
    end subroutine HMC_Data_Forcing_Point_Intake
    !------------------------------------------------------------------------------------------
    
end module HMC_Module_Data_Forcing_Point
!------------------------------------------------------------------------------------------