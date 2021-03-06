 !------------------------------------------------------------------------------------
! File:   HMC_Module_Phys_Snow.f90
!
! Author:   Simone Gabellani, Fabio Delogu
! Date:     20150714
!
! Physics subroutine(s) for snow model (S3M)
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Module Header
module HMC_Module_Phys_Snow

    !------------------------------------------------------------------------------------
    ! External module(s) 
    use HMC_Module_Namelist,            only:   oHMC_Namelist
    use HMC_Module_Vars_Loader,         only:   oHMC_Vars
    
    use HMC_Module_Tools_Debug
    
    use HMC_Module_Phys_Snow_Apps,       only:  HMC_Phys_Snow_Apps_TMean, &
                                                HMC_Phys_Snow_Apps_MeltingOL, &
                                                HMC_Phys_Snow_Apps_SWEAssim, &
                                                HMC_Phys_Snow_Apps_Age, &
                                                HMC_Phys_Snow_Apps_Albedo, &
                                                HMC_Phys_Snow_Apps_Rho
 
    ! Implicit none for all subroutines in this module
    implicit none
    !------------------------------------------------------------------------------------------

contains

    !------------------------------------------------------------------------------------------
    ! Subroutine to calculate snow updating 
    subroutine HMC_Phys_Snow_Cpl(iID, iRows, iCols)
        
        !-------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4)   :: iID, iRows, iCols
        integer(kind = 4)   :: iDaySteps, iDaySteps1Days, iDaySteps5Days
        integer(kind = 4)   :: iTime, iDtDataForcing, iDt
        integer(kind = 4)   :: iGlacierValue
        
        integer(kind = 4)   :: iFlagSnow, iFlagSnowAssim
        real(kind = 4)      :: dVarRhoW, dVarRhoSMax, dVarTRif
        real(kind = 4)      :: dDt
        
        character(len = 19) :: sTime
        
        integer(kind = 4), dimension(iRows, iCols)      :: a2iVarMask, a2iVarNature
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarDEM
        
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarRain, a2dVarSnowFall
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarTa, a2dVarRelHum, a2dVarIncRad 
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSnowHeight, a2dVarSnowKernel
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSnowCA, a2dVarSnowQA
        
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarArctUp
        
        integer(kind = 4), dimension(iRows, iCols)      :: a2iVarAgeS
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSWE, a2dVarRhoS, a2dVarAlbedoS
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarMeltingS, a2dVarMeltingSDayCum
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarMeltingSc
        
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarTaC120
        real(kind = 4), dimension(iRows, iCols)         :: a2dVarSepCoeff, a2dVarCloudFactor
        !-------------------------------------------------------------------------------------
        
        !-------------------------------------------------------------------------------------
        ! Snow cover area code(s)
        ! 'missing_data'    : { 'satellite' : 0, 	's3m' : -1 	},
        ! 'no_decision'     : { 'satellite' : 1, 	's3m' : -1 	},
        ! 'night'           : { 'satellite' : 11, 	's3m' : -1 	},
        ! 'no_snow'         : { 'satellite' : 25, 	's3m' : 0 	},
        ! 'lake'		    : { 'satellite' : 37, 	's3m' : -1 	},
        ! 'ocean'	    : { 'satellite' : 39, 	's3m' : -1 	},
        ! 'cloud'	    : { 'satellite' : 50, 	's3m' : 2 	},
        ! 'lake_ice'	    : { 'satellite' : 100, 	's3m' : -1 	},
        ! 'snow'		    : { 'satellite' : 200, 	's3m' : 1 	},
        ! 'detector_sat'    : { 'satellite' : 254, 	's3m' : -1 	},
        ! 'fill'		    : { 'satellite' : 255, 	's3m' : -1 	},
        !-------------------------------------------------------------------------------------
        
        !-------------------------------------------------------------------------------------
        ! Local variable(s) initialization
        a2iVarNature = -9999; a2dVarDEM = -9999.0;
        a2dVarArctUp = -9999.0;
        a2dVarRain = -9999.0; a2dVarTa = -9999.0; a2dVarRelHum = -9999.0; a2dVarIncRad = -9999.0; 
        a2dVarSnowHeight = -9999.0; a2dVarSnowKernel = -9999.0; a2dVarSnowCA = -9999.0; a2dVarSnowQA = -9999.0
        
        a2dVarSnowFall = 0.0;
        
        a2dVarSWE = -9999.0; a2dVarRhoS = -9999.0; a2dVarAlbedoS = -9999.0; a2iVarAgeS = -9999.0;
        a2dVarMeltingS = -9999.0; a2dVarMeltingSDayCum = -9999.0;
        
        a2dVarSepCoeff = -9999.0; a2dVarCloudFactor = -9999.0; a2dVarTaC120 = -9999.0;
        a2dVarMeltingSc = -9999.0;
        !-------------------------------------------------------------------------------------                         
                                     
        !-------------------------------------------------------------------------------------                      
        ! Snow flag value
        iFlagSnow = oHMC_Namelist(iID)%iFlagSnow
        ! Snow flag assimilation value
        iFlagSnowAssim = oHMC_Namelist(iID)%iFlagSnowAssim
        ! Water density [kg m^-3]
        dVarRhoW = oHMC_Namelist(iID)%dRhoW
        ! Define days steps [-]
        iDaySteps1Days = oHMC_Namelist(iID)%iDaySteps
        iDaySteps5Days = oHMC_Namelist(iID)%iDaySteps*5
        ! Glacier value 
        iGlacierValue = oHMC_Namelist(iID)%iGlacierValue
        ! Dt data forcing
        iDtDataForcing = oHMC_Namelist(iID)%iDtData_Forcing
        
        ! Threshold for snow melting [C]
        dVarTRif = 1.0
        
        ! Info start
        call mprintf(.true., iINFO_Verbose, ' Phys :: Snow model ... ' )
        !-------------------------------------------------------------------------------------  
        
        !-------------------------------------------------------------------------------------  
        ! Flag to activate snow physics
        if (iFlagSnow.eq.1) then
            
            ! Static variable(s)
            a2dVarDEM = oHMC_Vars(iID)%a2dDem
            a2iVarMask = oHMC_Vars(iID)%a2iMask
            a2iVarNature = oHMC_Vars(iID)%a2iNature
            a2dVarArctUp = oHMC_Vars(iID)%a2dArctUp

            ! Time information
            sTime = oHMC_Vars(iID)%sTimeStep
            iTime = oHMC_Vars(iID)%iTime

            ! Extracting dynamic forcing variable(s)
            a2dVarRain = oHMC_Vars(iID)%a2dRain
            a2dVarTa = oHMC_Vars(iID)%a2dTa
            a2dVarRelHum = oHMC_Vars(iID)%a2dRHum
            a2dVarIncRad = oHMC_Vars(iID)%a2dK
            a2dVarSnowHeight = oHMC_Vars(iID)%a2dSHeight
            a2dVarSnowKernel = oHMC_Vars(iID)%a2dSKernel
            a2dVarSnowCA = oHMC_Vars(iID)%a2dSCA
            a2dVarSnowQA = oHMC_Vars(iID)%a2dSQA

            ! Extracting state variable(s)
            a2dVarSWE = oHMC_Vars(iID)%a2dSWE
            a2dVarRhoS = oHMC_Vars(iID)%a2dRhoS
            a2dVarAlbedoS = oHMC_Vars(iID)%a2dAlbedo_Snow
            a2iVarAgeS = oHMC_Vars(iID)%a2iAge
            a2dVarMeltingS = oHMC_Vars(iID)%a2dMelting
            a2dVarMeltingSDayCum = oHMC_Vars(iID)%a2dMeltingDayCum
            
            ! Dt snow model step
            iDt = iDaySteps1Days; iDt = iDt*3600/86400
            !-------------------------------------------------------------------------------------
            
            !-----------------------------------------------------------------------------------------
            ! Debug
            if (iDEBUG.gt.0) then
                call mprintf(.true., iINFO_Extra, ' ========= SNOW START =========== ')  
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarRain, a2iVarMask, 'RAIN START ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarTa, a2iVarMask, 'TA START ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarRelHum, a2iVarMask, 'RELHUM START ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarIncRad, a2iVarMask, 'INCRAD START ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowHeight, a2iVarMask, 'SNOWH START ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowKernel, a2iVarMask, 'SNOWK START ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowCA, a2iVarMask, 'SNOWCA START ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowQA, a2iVarMask, 'SNOWQA START ') )
                
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowFall, a2iVarMask, 'SNOWFALL START ') )

                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWE, a2iVarMask, 'SWE START ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarRhoS, a2iVarMask, 'RHOS START ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarAlbedoS, a2iVarMask, 'ALBEDOS START ') )
                call mprintf(.true., iINFO_Extra, checkvar(real(a2iVarAgeS), a2iVarMask, 'AGES START ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingS, a2iVarMask, 'MELTINGS START ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingSDayCum, a2iVarMask, 'MELTINGS DAYCUM START ') )
                call mprintf(.true., iINFO_Extra, ' ') 
            endif
            !-----------------------------------------------------------------------------------------
            
            !-------------------------------------------------------------------------------------
            ! Check variable(s)
            where(a2dVarSWE.lt.0.0) a2dVarSWE = 0.0 ! just in case
            !-------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------
            ! Define cloud factor
            where( (a2dVarDEM.ge.0.0) .and. (a2dVarRain.gt.0.2) ) a2dVarCloudFactor = 1.1
            where( (a2dVarDEM.ge.0.0) .and. (a2dVarRain.ge.1.0) ) a2dVarCloudFactor = 1.2
            where( (a2dVarDEM.ge.0.0) .and. (a2dVarRain.le.0.2) ) a2dVarCloudFactor = 0.9
            !-------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------
            ! Compute SepCoeff
            where(a2dVarDEM.ge.0.0)
                a2dVarSepCoeff = 1/(1 + exp(22 - 2.7*a2dVarTa - 0.2*a2dVarRelHum))
            endwhere    
            ! Compute SnowFall
            where(a2dVarDEM.ge.0.0 .and. a2dVarRain.gt.0.0)
                a2dVarSnowFall = (1 - a2dVarSepCoeff)*a2dVarRain 
            endwhere
            ! Compute SWE and update Rain
            where( (a2dVarDEM.ge.0.0) .and. (a2dVarRain.gt.0.0) )
                a2dVarSWE = a2dVarSWE + a2dVarSnowFall	
                a2dVarRain = a2dVarSepCoeff*a2dVarRain
            endwhere
            ! Convert SnowHeigth units (from cm to mm)
            where(a2dVarDEM.ge.0.0)
                a2dVarSnowHeight = a2dVarSnowHeight*10;
            endwhere
            !-------------------------------------------------------------------------------------
            
            !-------------------------------------------------------------------------------------
            ! Update daily cumulated snow melting 
            if ( (sTime(12:13).eq.'00') .and. (iTime.gt.1) ) then
                ! Update cumulated daily melting
                a2dVarMeltingSDayCum = 0.0
            endif
            !-------------------------------------------------------------------------------------
            
            !-------------------------------------------------------------------------------------
            ! Call subroutine to compute average temperature over n days
            call HMC_Phys_Snow_Apps_TMean(iID, iRows, iCols, & 
                                            iDaySteps5Days, &
                                            a2iVarMask, &
                                            a2dVarTa, a2dVarTaC120)
            ! Debug
            ! call mprintf(.true., iINFO_Extra, checkvar(a2dVarTaC120, a2iVarMask, 'TA120') ) 
            !-------------------------------------------------------------------------------------
                                     
            !-------------------------------------------------------------------------------------
            ! Call subroutine to compute snow age
            call HMC_Phys_Snow_Apps_Age(iID, iRows, iCols, &
                                        sTime, iTime, &
                                        a2dVarDem, &
                                        a2dVarSnowFall, a2dVarSWE, &
                                        a2iVarAgeS) 
            ! Debug
            ! call mprintf(.true., iINFO_Extra, checkvar(float(a2iVarAgeS), a2iVarMask, 'AGES') ) 
            !-------------------------------------------------------------------------------------     
                                      
            !-------------------------------------------------------------------------------------
            ! Call subroutine to compute snow albedo
            call HMC_Phys_Snow_Apps_Albedo(iID, iRows, iCols, &
                                           sTime, iTime, iGlacierValue, &
                                           a2dVarDem, &
                                           a2dVarTaC120, a2iVarAgeS, &
                                           a2dVarAlbedoS)
            ! Debug                            
            ! call mprintf(.true., iINFO_Extra, checkvar(a2dVarAlbedoS, a2iVarMask, 'ALBEDOS') )                              
            !-------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------
            ! Call subroutine to compute snow density
            call HMC_Phys_Snow_Apps_Rho(iID, iRows, iCols, &
                                            sTime, iTime, &
                                            iDt, &
                                            a2dVarDem, &
                                            a2dVarTa, a2dVarSnowFall, a2dVarSWE, a2dVarMeltingSDayCum, &
                                            a2dVarRhoS) 

            ! Nullify snowfall
            !a2dVarSnowFall = 0.0
            
            ! Debug
            !call mprintf(.true., iINFO_Extra, checkvar(a2dVarRhoS, a2iVarMask, 'RHOS ') )
            !-------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------
            ! Call subroutine to compute open-loop melting
            call HMC_Phys_Snow_Apps_MeltingOL(iID, iRows, iCols, iDtDataForcing, iDaySteps1Days, &
                                              sTime, iTime, &
                                              iGlacierValue, dVarRhoW, dVarTRif, &
                                              a2dVarDem, a2iVarNature, a2dVarArctUp, &
                                              a2dVarTa, a2dVarIncRad, &
                                              a2dVarTaC120, a2dVarCloudFactor, &
                                              a2dVarAlbedoS, a2dVarSWE, a2dVarMeltingS)
            !-------------------------------------------------------------------------------------
                                          
            !-------------------------------------------------------------------------------------                                   
            ! Compute SWE, melting and rain
            where( (a2dVarDem.gt.0.0) .and. (a2dVarSWE.gt.0.0) .and. &
                   (a2dVarSWE.le.a2dVarMeltingS) .and. (a2iVarNature.ne.iGlacierValue) ) 
                   
                ! All melted snow condition
                a2dVarMeltingS = a2dVarSWE
                a2dVarSWE = 0.0
                a2dVarRain = a2dVarRain + a2dVarMeltingS

            elsewhere( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.gt.0.0) .and. (a2iVarNature.ne.iGlacierValue))

                a2dVarSWE = a2dVarSWE - a2dVarMeltingS
                a2dVarRain = a2dVarRain + a2dVarMeltingS
            endwhere

            where( (a2dVarDem.ge.0.0) .and. (a2iVarNature.ne.iGlacierValue) )
                
                ! Glacier(s) condition
                a2dVarSWE = a2dVarSWE - a2dVarMeltingS
                a2dVarRain = a2dVarRain + a2dVarMeltingS

            endwhere
            
            ! Check SWE values
            where( (a2dVarDem.ge.0.0) .and. (a2dVarSWE.lt.0.0) ) a2dVarSWE = 0.0
            
            ! Compute daily cumulated melting
            where( (a2dVarDem.ge.0.0) .and. (a2dVarMeltingS.GT.0.0) )
                a2dVarMeltingSDayCum = a2dVarMeltingSDayCum + a2dVarMeltingS
            endwhere
            !-------------------------------------------------------------------------------------

            !-------------------------------------------------------------------------------------
            ! Call subroutine to compute SWE assimilation
            call mprintf(.true., iINFO_Verbose, ' Phys :: Snow :: Assimilation ... ' )
            if (iFlagSnowAssim.eq.1) then 
                
                !-------------------------------------------------------------------------------------
                ! Check forcing(s) to use assimilation method
                if ( any(a2dVarSnowHeight.ne.-9999.0) .and. any(a2dVarSnowKernel.ne.-9999.0) .and. &
                     any(a2dVarSnowCA.ne.-9999.0) .and. any(a2dVarSnowQA.ne.-9999.0) ) then
                    
                    !-------------------------------------------------------------------------------------
                    ! Call subroutine to compute SWE assimilation
                    call HMC_Phys_Snow_Apps_SWEAssim(iID, iRows, iCols, &
                                                        dVarRhoW, &
                                                        a2iVarMask, &
                                                        a2dVarSnowHeight, a2dVarSnowKernel, &
                                                        a2dVarSnowCA, a2dVarSnowQA, &
                                                        a2dVarSWE, a2iVarAgeS, a2dVarAlbedoS, a2dVarRhoS)

                    ! Info start assimilation
                    call mprintf(.true., iINFO_Verbose, ' Phys :: Snow :: Assimilation ... OK' )
                    !-------------------------------------------------------------------------------------
                                
                else
              
                    !-------------------------------------------------------------------------------------
                    ! Info assimilation no data available
                    call mprintf(.true., iINFO_Verbose, ' Phys :: Snow : Assimilation ... DATA NO AVAILABLE ' )
                    !-------------------------------------------------------------------------------------
                    
                endif
                !-------------------------------------------------------------------------------------

            else
                
                !-------------------------------------------------------------------------------------
                ! Info assimilation no data available
                call mprintf(.true., iINFO_Verbose, ' Phys :: Snow : Assimilation ... NOT ACTIVATED ' )
                !-------------------------------------------------------------------------------------
                
            endif
            !-------------------------------------------------------------------------------------
            
            !-------------------------------------------------------------------------------------
            ! Check variable(s) space domain and physical boundaries 
            where(a2dVarDem.lt.0.0)
                a2dVarRain = -9999.0;  a2dVarSWE = -9999.0; a2dVarSnowFall = -9999.0; a2dVarRhos = -9999.0;
                a2dVarTaC120 = -9999.0; a2iVarAgeS = -9999; a2dVarMeltingS = -9999.0; a2dVarMeltingSDayCum = -9999.0
            endwhere
            
            where( (a2dVarDem.gt.0.0) .and. (a2dVarSWE.eq.0.0) ) 
                a2iVarAgeS = 0; a2dVarAlbedoS = 0.0; a2dVarRhos = 0.0; a2dVarMeltingS = 0.0;
            endwhere
            !-------------------------------------------------------------------------------------
            
            !-----------------------------------------------------------------------------------------
            ! Debug
            if (iDEBUG.gt.0) then
                call mprintf(.true., iINFO_Extra, ' ')
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarRain, a2iVarMask, 'RAIN END ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarTaC120, a2iVarMask, 'TA MEAN END ') )
                
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSnowFall, a2iVarMask, 'SNOWFALL END ') )
                
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarSWE, a2iVarMask, 'SWE END ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarRhoS, a2iVarMask, 'RHOS END ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarAlbedoS, a2iVarMask, 'ALBEDOS END ') )
                call mprintf(.true., iINFO_Extra, checkvar(real(a2iVarAgeS), a2iVarMask, 'AGES END ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingS, a2iVarMask, 'MELTINGS END ') )
                call mprintf(.true., iINFO_Extra, checkvar(a2dVarMeltingSDayCum, a2iVarMask, 'MELTINGS DAYCUM END ') )
                call mprintf(.true., iINFO_Extra, ' ========= SNOW END =========== ')  
            endif
            !-----------------------------------------------------------------------------------------
            
            !-----------------------------------------------------------------------------------------
            ! Update forcing variable(s)
            oHMC_Vars(iID)%a2dRain = a2dVarRain

            ! Update state variable(s)
            oHMC_Vars(iID)%a2dSWE = a2dVarSWE
            oHMC_Vars(iID)%a2dRhoS = a2dVarRhoS
            oHMC_Vars(iID)%a2dAlbedo_Snow = a2dVarAlbedoS
            oHMC_Vars(iID)%a2iAge = a2iVarAgeS
            oHMC_Vars(iID)%a2dMelting = a2dVarMeltingS
            oHMC_Vars(iID)%a2dMeltingSc = a2dVarMeltingSc
            oHMC_Vars(iID)%a2dMeltingDayCum = a2dVarMeltingSDayCum

            ! Info end
            call mprintf(.true., iINFO_Verbose, ' Phys :: Snow model ... OK' )
            !-----------------------------------------------------------------------------------------
        
        else
            
            !-----------------------------------------------------------------------------------------
            ! Update state variable(s)
            oHMC_Vars(iID)%a2dSWE = -9999.0
            oHMC_Vars(iID)%a2dRhoS = -9999.0
            oHMC_Vars(iID)%a2dAlbedo_Snow = -9999.0
            oHMC_Vars(iID)%a2iAge = -9999.0
            
            oHMC_Vars(iID)%a2dMelting = -9999.0
            oHMC_Vars(iID)%a2dMeltingSc = -9999.0
            oHMC_Vars(iID)%a2dMeltingDayCum  = -9999.0
            
            ! Info end
            call mprintf(.true., iINFO_Verbose, ' Phys :: Snow model ... NOT ACTIVATED' )
            !-----------------------------------------------------------------------------------------
            
        endif
        !-----------------------------------------------------------------------------------------
        
        !-----------------------------------------------------------------------------------------
        ! Debug
        if ( (sTime(12:13).eq.'00') .and. (iTime.gt.1) ) then
            write(*,*)
        endif
        !-----------------------------------------------------------------------------------------
        
    end subroutine HMC_Phys_Snow_Cpl
    !------------------------------------------------------------------------------------------

end module HMC_Module_Phys_Snow
!------------------------------------------------------------------------------------