!------------------------------------------------------------------------------------
! File:         HMC_Module_Function_IO.f90
! Author(s):    Fabio Delogu, Francesco Silvestro, Simone Gabellani
! Created on May 11, 2015, 10:27 AM
!
! Module to perform IO actions
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Module Header
module HMC_Module_Tools_IO
    
    !------------------------------------------------------------------------------------
    ! External module(s) and implicit none
#ifdef LIB_NC
    use netcdf
#endif

    use HMC_Module_Tools_Debug
    
    implicit none
    
    integer(kind = 4), parameter    :: iDeflateLevelNC = 9;
    !------------------------------------------------------------------------------------

contains 

    !------------------------------------------------------------------------------------
    ! Subroutine for getting 1d variable NC
#ifdef LIB_NC
    subroutine HMC_Tools_IO_Get1d_NC(sVarName, iFileID, a1dVar, sVarUnits, iX, &
                               bFatalError, iErr) 

        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4) :: iFileID, iVarID
        integer(kind = 4) :: iX

        character(len = 256) :: sVarName
        character(len = 256), intent(out)               :: sVarUnits
        real(kind = 4), dimension(iX), intent(out)      :: a1dVar 

        logical, intent(in) :: bFatalError
        integer(kind = 4)   :: iErr, iRet
        !------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------
        ! Initializing units attribute
        sVarUnits = "";
        iVarID = 0;
        a1dVar = -9999.0;

        ! Finding 2d variable ID
        iRet = nf90_inq_varid(iFileID,  trim(sVarName), iVarID)
        if (iRet /= 0) then
            if (bFatalError) then
                call mprintf(.true., iERROR, ' INQUIRE FAILED :: Get1d_NC ---> nf90_inq_varid (varname: '//trim(sVarName)//')')
            else
                call mprintf(.true., iWARN, ' INQUIRE FAILED :: Get1d_NC ---> nf90_inq_varid (varname: '//trim(sVarName)//')')
                iErr = iRet
                return
            endif
        else
            call mprintf(.true., iINFO_Extra, ' INQUIRE OK :: Get1d_NC ---> nf90_inq_varid (varname: '//trim(sVarName)//')')
        endif

        ! Extracting 2d variable
        iRet = nf90_get_var(iFileID, iVarID, a1dVar) 
        if (iRet /= 0) then
            if (bFatalError) then
                call mprintf(.true., iERROR, ' READ FAILED :: Get1d_NC ---> nf90_get_var (varname: '//trim(sVarName)//')')
            else
                call mprintf(.true., iWARN, ' READ FAILED :: Get1d_NC ---> nf90_get_var (varname: '//trim(sVarName)//')')
            endif  
        else
            call mprintf(.true., iINFO_Extra, ' READ OK :: Get1d_NC ---> nf90_get_var (varname: '//trim(sVarName)//')')
        endif
 
        ! Algorithm conclusion
        iErr = 0;
        !------------------------------------------------------------------------------------

    end subroutine HMC_Tools_IO_Get1d_NC
#endif
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine for getting 2d variable NC
#ifdef LIB_NC
    subroutine HMC_Tools_IO_Get2d_NC(sVarName, iFileID, a2dVar, sVarUnits, iX, iY, &
                               bFatalError, iErr) 

        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4) :: iFileID, iVarID
        integer(kind = 4) :: iX, iY

        character(len = 256) :: sVarName
        character(len = 256), intent(out)                   :: sVarUnits
        real(kind = 4), dimension(iX, iY), intent(out)      :: a2dVar 

        logical, intent(in) :: bFatalError
        integer(kind = 4)   :: iErr, iRet
        !------------------------------------------------------------------------------------

        !------------------------------------------------------------------------------------
        ! Initializing units attribute
        sVarUnits = "";
        iVarID = 0;
        a2dVar = -9999.0;

        ! Finding 2d variable ID
        iRet = nf90_inq_varid(iFileID,  trim(sVarName), iVarID)
        if (iRet /= 0) then
            if (bFatalError) then
                call mprintf(.true., iERROR, ' INQUIRE FAILED :: Get2d_NC ---> nf90_inq_varid (varname: '//trim(sVarName)//')')
            else
                call mprintf(.true., iWARN, ' INQUIRE FAILED :: Get2d_NC ---> nf90_inq_varid (varname: '//trim(sVarName)//')')
                iErr = iRet
                return
            endif
        else
            call mprintf(.true., iINFO_Extra, ' INQUIRE OK :: Get2d_NC ---> nf90_inq_varid (varname: '//trim(sVarName)//')')
        endif

        ! Extracting 2d variable
        iRet = nf90_get_var(iFileID, iVarID, a2dVar) 
        if (iRet /= 0) then
            if (bFatalError) then
                call mprintf(.true., iERROR, ' READ FAILED ::  Get2d_NC ---> nf90_get_var (varname: '//trim(sVarName)//')')
            else
                call mprintf(.true., iWARN, ' READ FAILED ::  Get2d_NC ---> nf90_get_var (varname: '//trim(sVarName)//')')
            endif
        else
            call mprintf(.true., iINFO_Extra, ' READ OK :: Get2d_NC ---> nf90_get_var (varname: '//trim(sVarName)//')')
        endif
 
        ! Algorithm conclusion
        iErr = 0;
        !------------------------------------------------------------------------------------

    end subroutine HMC_Tools_IO_Get2d_NC
#endif
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to get 3d variable
#ifdef LIB_NC
    subroutine HMC_Tools_IO_Get3d_NC(sVarName, iFileID, a3dVar, sVarUnits, iT, iX, iY, &
                                     bFatalError, iErr) 

        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4) :: iFileID, iVarID
        integer(kind = 4) :: iX, iY, iT

        character(len = 256) :: sVarName
        character(len = 256), intent(out)                   :: sVarUnits
        real(kind = 4), dimension(iX, iY, iT), intent(out)  :: a3dVar 

        logical, intent(in) :: bFatalError
        integer(kind = 4)   :: iErr, iRet
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Initializing units attribute
        sVarUnits = "";
        iVarID = 0;
        a3dVar = -9999.0;

        ! Finding 3d variable ID
        iRet = nf90_inq_varid(iFileID,  trim(sVarName), iVarID)
        if (iRet /= 0) then
            if (bFatalError) then
                call mprintf(.true., iERROR, ' INQUIRE FAILED :: Get3d_NC ---> nf90_inq_varid (varname: '//trim(sVarName)//')')
            else
                call mprintf(.true., iWARN, ' INQUIRE FAILED :: Get3d_NC ---> nf90_inq_varid (varname: '//trim(sVarName)//')')
                iErr = iRet
                return
            endif
        else
            call mprintf(.true., iINFO_Extra, ' INQUIRE OK :: Get3d_NC ---> nf90_inq_varid (varname: '//trim(sVarName)//')')
        endif

        ! Extracting 3d variable
        iRet = nf90_get_var(iFileID, iVarID, a3dVar) 
        if (iRet /= 0) then
            if (bFatalError) then
                call mprintf(.true., iERROR, ' READ FAILED ::  Get3d_NC ---> nf90_get_var (varname: '//trim(sVarName)//')')
            else
                call mprintf(.true., iWARN, ' READ FAILED ::  Get3d_NC ---> nf90_get_var (varname: '//trim(sVarName)//')')
            endif
        else
            call mprintf(.true., iINFO_Extra, ' READ OK :: Get3d_NC ---> nf90_get_var (varname: '//trim(sVarName)//')')
        endif
 
        ! Algorithm conclusion
        iErr = 0;
        !------------------------------------------------------------------------------------
        
    end subroutine HMC_Tools_IO_Get3d_NC
#endif
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to put 1d variable
#ifdef LIB_NC
    subroutine HMC_Tools_IO_Put1d_NC(iFileID, iID_DimX, &
                               sVarName, sVarNameLong, sVarDescription, &
                               sVarUnits, &
                               sVarCoords, sVarGridMap, &
                               dVarMissingValue, &
                               iVarX, a1dVar) 
            
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4) :: iFileID
        integer(kind = 4) :: iID_DimX
        
        character(len = 256) :: sVarName, sVarNameLong, sVarDescription
        character(len = 256) :: sVarUnits
        character(len = 256) :: sVarCoords, sVarGridMap
        
        real(kind = 4)       :: dVarMissingValue
        
        integer(kind = 4) :: iVarID, iVarX
        real(kind = 4), dimension (iVarX)   :: a1dVar
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Initializing variable(s)
        iVarID = 0;
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Definition mode ON - Data mode OFF
        call check( nf90_redef(iFileID) )
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Allocating variable 
        call check( nf90_def_var(iFileID, trim(sVarName), nf90_float, (/iID_DimX/), iVarID, &
                                 deflate_level = iDeflateLevelNC) )
        
        !  Writing variable attribute(s)
        call check( nf90_put_att(iFileID, iVarID, 'long_name'     , trim(sVarNameLong)) )
        call check( nf90_put_att(iFileID, iVarID, 'description'   , trim(sVarDescription)) )
        call check( nf90_put_att(iFileID, iVarID, 'units'         , trim(sVarUnits)) )
        call check( nf90_put_att(iFileID, iVarID, 'coordinates'   , trim(sVarCoords)) )
        call check( nf90_put_att(iFileID, iVarID, 'grid_mapping'  , trim(sVarGridMap)) )
        call check( nf90_put_att(iFileID, iVarID, 'missing_value' , real(dVarMissingValue)) )
        
        ! Writing variable value(s)
        call check( nf90_inq_varid(iFileID, trim(sVarName), iVarID) )
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Definition mode OFF - Data mode ON
        call check( nf90_enddef(iFileID))
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Writing variable(s)
        call check( nf90_put_var(iFileID, iVarID, a1dVar, &
                                 start = (/ 1 /), count = (/ iVarX /)) )
        !------------------------------------------------------------------------------------

    end subroutine HMC_Tools_IO_Put1d_NC
#endif
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to put 2d variable
#ifdef LIB_NC
    subroutine HMC_Tools_IO_Put2d_NC(iFileID, iID_DimX, iID_DimY, &
                               sVarName, sVarNameLong, sVarDescription, &
                               sVarUnits, &
                               sVarCoords, sVarGridMap, &
                               dVarMissingValue, &
                               iVarX, iVarY, a2dVar) 
            
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4) :: iFileID
        integer(kind = 4) :: iID_DimX, iID_DimY
        
        character(len = 256) :: sVarName, sVarNameLong, sVarDescription
        character(len = 256) :: sVarUnits
        character(len = 256) :: sVarCoords, sVarGridMap
        
        real(kind = 4)       :: dVarMissingValue
        
        integer(kind = 4) :: iVarID, iVarX, iVarY
        real(kind = 4), dimension (iVarX, iVarY)   :: a2dVar
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Subroutine input example
        ! sVarName = 'SoilM'
        ! sVarNameLong = 'soil_moisture'
        ! sVarDescription = 'moisture content'
        ! sVarUnits = 'm^3/m^3'
        ! sVarCoords = 'x y z'
        ! sVarGridMap = 'lambert_conformal_conic'
        ! dVarMissingValue = -9E15
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Initializing variable(s)
        iVarID = 0;
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Definition mode ON - Data mode OFF
        call check( nf90_redef(iFileID) )
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Allocating variable 
        call check( nf90_def_var(iFileID, trim(sVarName), nf90_float, (/iID_DimX, iID_DimY/), iVarID, &
                                 deflate_level = iDeflateLevelNC) )
        
        !  Writing variable attribute(s)
        call check( nf90_put_att(iFileID, iVarID, 'long_name'     , trim(sVarNameLong)) )
        call check( nf90_put_att(iFileID, iVarID, 'description'   , trim(sVarDescription)) )
        call check( nf90_put_att(iFileID, iVarID, 'units'         , trim(sVarUnits)) )
        call check( nf90_put_att(iFileID, iVarID, 'coordinates'   , trim(sVarCoords)) )
        call check( nf90_put_att(iFileID, iVarID, 'grid_mapping'  , trim(sVarGridMap)) )
        call check( nf90_put_att(iFileID, iVarID, 'missing_value' , real(dVarMissingValue)) )
        
        ! Writing variable value(s)
        call check( nf90_inq_varid(iFileID, trim(sVarName), iVarID) )
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Definition mode OFF - Data mode ON
        call check( nf90_enddef(iFileID))
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Writing variable(s)
        call check( nf90_put_var(iFileID, iVarID, a2dVar, &
                                 start = (/ 1, 1/), count = (/ iVarX, iVarY/)) )
        !------------------------------------------------------------------------------------

    end subroutine HMC_Tools_IO_Put2d_NC
#endif
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to put 3d variable
#ifdef LIB_NC
    subroutine HMC_Tools_IO_Put3d_NC(iFileID, iID_DimX, iID_DimY, iID_DimT, &
                                     sVarName, sVarNameLong, sVarDescription, &
                                     sVarUnits, &
                                     sVarCoords, sVarGridMap, &
                                     dVarMissingValue, &
                                     iVarX, iVarY, iVarT, a3dVar) 
                        
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4)       :: iFileID
        integer(kind = 4)       :: iID_DimX, iID_DimY, iID_DimT
        
        character(len = 256)    :: sVarName, sVarNameLong, sVarDescription
        character(len = 256)    :: sVarUnits
        character(len = 256)    :: sVarCoords, sVarGridMap
        
        real(kind = 4)          :: dVarMissingValue
        
        integer(kind = 4)       :: iT, iVarID, iVarX, iVarY, iVarT
        real(kind = 4), dimension (iVarX, iVarY, iVarT)   :: a3dVar
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Initializing variable(s)
        iVarID = 0;
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Definition mode ON - Data mode OFF
        call check( nf90_redef(iFileID) )
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Allocating variable 
        call check( nf90_def_var(iFileID, trim(sVarName), nf90_float, &
                                 (/iID_DimX, iID_DimY, iID_DimT/), iVarID, &
                                 deflate_level = iDeflateLevelNC) )
        
        !  Writing variable attribute(s)
        call check( nf90_put_att(iFileID, iVarID, 'long_name'     , trim(sVarNameLong)) )
        call check( nf90_put_att(iFileID, iVarID, 'description'   , trim(sVarDescription)) )
        call check( nf90_put_att(iFileID, iVarID, 'units'         , trim(sVarUnits)) )
        call check( nf90_put_att(iFileID, iVarID, 'coordinates'   , trim(sVarCoords)) )
        call check( nf90_put_att(iFileID, iVarID, 'grid_mapping'  , trim(sVarGridMap)) )
        call check( nf90_put_att(iFileID, iVarID, 'missing_value' , real(dVarMissingValue)) )
        
        ! Writing variable value(s)
        call check( nf90_inq_varid(iFileID, trim(sVarName), iVarID) )
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Definition mode OFF - Data mode ON
        call check( nf90_enddef(iFileID))
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Writing variable(s)
        call check( nf90_put_var(iFileID, iVarID, a3dVar, &
                                 start = (/ 1, 1, 1/), count = (/ iVarX, iVarY, iVarT/)) )
        !------------------------------------------------------------------------------------
                                                  
    end subroutine HMC_Tools_IO_Put3d_NC
#endif
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to get 1d variable ASCII (1 columns)
    subroutine HMC_Tools_IO_Get1d_ASCII(sFileName, a1dVar, iX, &
                               bFatalError, iErr)
        
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4) :: iX, iI, iJ

        character(len = 256) :: sFileName
        real(kind = 4), dimension(iX), intent(out)      :: a1dVar 

        logical, intent(in) :: bFatalError
        logical             :: bFileExist
        integer(kind = 4)   :: iErr, iRet
        !------------------------------------------------------------------------------------
                               
        !------------------------------------------------------------------------------------
        ! Open and read ascii file
        inquire (file = trim(sFileName), exist = bFileExist, iostat = iRet)
        if (.not. bFileExist ) then
            
            if (bFatalError) then
                call mprintf(.true., iERROR, &
                            ' INQUIRE FAILED :: Get1d_ASCII ---> file not found (filename: '//trim(sFileName)//')')
            else
                call mprintf(.true., iWARN, &
                            ' INQUIRE FAILED :: Get1d_ASCII ---> file not found (filename: '//trim(sFileName)//')')
                iErr = -1
                return
            endif
            
        elseif ( bFileExist ) then

            open(unit=1, file=trim(sFileName), status='old', iostat=iRet)
            if (iRet == 0) then 

                call mprintf(.true., iINFO_Extra, &
                            ' READ OK :: Get1d_ASCII ---> file readable (filename: '//trim(sFileName)//')')
                do iI = 1,iX
                    read(1,*) a1dVar(iI) 
                enddo
                close(1)
                
            elseif (iRet /= 0) then
                
                if (bFatalError) then
                    call mprintf(.true., iERROR, &
                                ' READ FAILED :: Get1d_ASCII ---> file not readable (filename: '//trim(sFileName)//')')
                else
                    call mprintf(.true., iWARN, &
                                ' READ FAILED :: Get1d_ASCII ---> file not readable (filename: '//trim(sFileName)//')')
                    iErr = iRet
                    return
                endif
                   
            endif

        endif
        !------------------------------------------------------------------------------------
        
    end subroutine HMC_Tools_IO_Get1d_ASCII
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    ! Subroutine to get 2d variable ASCII (2 columns)
    subroutine HMC_Tools_IO_Get2d_ASCII(sFileName, a1sVar, a1dVar, iX,  bFatalError, iErr)
        
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4) :: iX, iI, iJ

        character(len = 256) :: sFileName
        real(kind = 4), dimension(iX), intent(out)      :: a1dVar 
        character(len = 256), dimension(iX)             :: a1sVar 

        logical, intent(in) :: bFatalError
        logical             :: bFileExist
        integer(kind = 4)   :: iErr, iRet
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a1sVar = ''; a1dVar = -9999.0
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Open and read ascii file
        inquire (file = trim(sFileName), exist = bFileExist, iostat = iRet)
        if (.not. bFileExist ) then
            
            if (bFatalError) then
                call mprintf(.true., iERROR, &
                            ' INQUIRE FAILED :: Get2d_ASCII ---> file not found (filename: '//trim(sFileName)//')')
            else
                call mprintf(.true., iWARN, &
                            ' INQUIRE FAILED :: Get1d_ASCII ---> file not found (filename: '//trim(sFileName)//')')
                iErr = -1
                return
            endif
            
        elseif ( bFileExist ) then
            
            open(unit=1, file=trim(sFileName), status='old', iostat=iRet)
            
            call mprintf(.true., iINFO_Extra, &
                            ' READ OK :: Get2d_ASCII ---> file readable (filename: '//trim(sFileName)//')')
            if (iRet == 0) then 
                
                do iI = 1, iX
                    read(1,*) a1sVar(iI), a1dVar(iI)
                enddo
                close(1)
                
            elseif (iRet /= 0) then
                
                if (bFatalError) then
                    call mprintf(.true., iERROR, &
                            ' READ FAILED :: Get2d_ASCII ---> file not readable (filename: '//trim(sFileName)//')')
                else
                    call mprintf(.true., iWARN, &
                            ' READ FAILED :: Get2d_ASCII ---> file not readable (filename: '//trim(sFileName)//')')
                    iErr = iRet
                    return
                endif
                   
            endif

        endif
        !------------------------------------------------------------------------------------

    end subroutine HMC_Tools_IO_Get2d_ASCII
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to get ArcGrid variable ASCII
    subroutine HMC_Tools_IO_GetArcGrid_ASCII(sFileName, a2dVar, iX, iY, &
                               bFatalError, iErr)
        
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4) :: iX, iY, iI, iJ

        character(len = 256) :: sFileName
        !real(kind = 4), dimension(iX, iY), intent(out)      :: a2dVar 
        real(kind = 4), dimension(iY, iX), intent(out)      :: a2dVar 
        
        logical, intent(in) :: bFatalError
        logical             :: bFileExist
        integer(kind = 4)   :: iErr, iRet
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Init variable(s)
        iErr = 0; a2dVar = 0.0
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Open and read ascii file
        inquire (file = trim(sFileName), exist = bFileExist, iostat = iRet)
        if (.not. bFileExist ) then
            
            if (bFatalError) then
                call mprintf(.true., iERROR, &
                            ' INQUIRE FAILED :: GetArcGrid ---> file not found (filename: '//trim(sFileName)//')')
            else
                call mprintf(.true., iWARN, &
                            ' INQUIRE FAILED :: GetArcGrid ---> file not found (filename: '//trim(sFileName)//')')
                iErr = -1
                return
            endif
            
        elseif ( bFileExist ) then
            
            open(unit=1, file=trim(sFileName), status='old', iostat=iRet)
            
            if (iRet == 0) then 
                
                call mprintf(.true., iINFO_Extra, &
                            ' READ OK :: GetArcGrid ---> file readable (filename: '//trim(sFileName)//')')
                
                do iI = 1,6; read(1,*); enddo
                !do iI = 1,iX,1
                !        read(1,*) (a2dVar(iX - iI + 1,iJ), iJ = 1,iY)
                !enddo
                
                do iI = 1,iY,1
                        read(1,*) (a2dVar(iY - iI + 1,iJ), iJ = 1,iX)
                enddo
                
                close(1)
                iErr = 0
                
            elseif (iRet /= 0) then
                
                if (bFatalError) then
                    call mprintf(.true., iERROR, &
                            ' READ FAILED :: GetArcGrid ---> file not readable (filename: '//trim(sFileName)//')')
                else
                    call mprintf(.true., iWARN, &
                            ' READ FAILED :: GetArcGrid ---> file not readable (filename: '//trim(sFileName)//')')
                    iErr = iRet
                    return
                endif
                   
            endif

        endif
        !------------------------------------------------------------------------------------
        
    end subroutine HMC_Tools_IO_GetArcGrid_ASCII
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to put 1d variable ASCII
    subroutine HMC_Tools_IO_Put1d_ASCII(sFileName, a1dVar, iX, bFatalError, iErr, sFMT)
        
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4) :: iX, iI

        character(len = 256) :: sFileName
        real(kind = 4), dimension(iX)   :: a1dVar 

        logical, intent(in) :: bFatalError
        logical             :: bFileExist
        integer(kind = 4)   :: iErr, iRet
        
        character(len = 20) :: sFMT
        !------------------------------------------------------------------------------------
                               
        !------------------------------------------------------------------------------------
        ! Open and write ASCII file
        open(unit = 20, file = trim(sFileName), status = 'new', iostat = iRet)
        inquire (file = trim(sFileName), exist = bFileExist, iostat = iRet)
        if (.not. bFileExist ) then
            
            if (bFatalError) then
                call mprintf(.true., iERROR, &
                            ' INQUIRE FAILED :: Put1d_ASCII ---> file not readable (filename: '//trim(sFileName)//')')
            else
                call mprintf(.true., iWARN, &
                            ' INQUIRE FAILED :: Put1d_ASCII ---> file not readable (filename: '//trim(sFileName)//')')
                iErr = 0
                return
            endif
            
        elseif ( bFileExist ) then
            
            if (iRet == 0) then 
                
                call mprintf(.true., iINFO_Extra, &
                            ' WRITE OK :: Put1d_ASCII ---> file writable (filename: '//trim(sFileName)//')')
                do iI = 1,iX
                    write(20 , sFMT) a1dVar(iI) 
                enddo
                close(20)
                
            elseif (iRet /= 0) then
                
                if (bFatalError) then
                    call mprintf(.true., iERROR, &
                            ' WRITE FAILED :: Put1d_ASCII ---> file not writable (filename: '//trim(sFileName)//')')
                else
                    call mprintf(.true., iWARN, &
                            ' WRITE FAILED :: Put1d_ASCII ---> file not writable (filename: '//trim(sFileName)//')')
                    iErr = iRet
                    return
                endif
                   
            endif

        endif
        !------------------------------------------------------------------------------------
        
    end subroutine HMC_Tools_IO_Put1d_ASCII
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to put 2d variable ASCII
    subroutine HMC_Tools_IO_Put2d_ASCII(sFileName, a2dVar, iX, iY, &
                               bFatalError, iErr, sFMT)
        
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4) :: iX, iY, iI, iJ

        character(len = 256) :: sFileName
        real(kind = 4), dimension(iX, iY) :: a2dVar 

        logical, intent(in) :: bFatalError
        logical             :: bFileExist
        integer(kind = 4)   :: iErr, iRet
        
        character(len = 40), optional :: sFMT
        !------------------------------------------------------------------------------------
                               
        !------------------------------------------------------------------------------------
        ! Open and write ASCII file
        open(unit = 20, file = trim(sFileName), status = 'new', iostat = iRet)
        inquire (file = trim(sFileName), exist = bFileExist, iostat = iRet)
        if (.not. bFileExist ) then
            
            if (bFatalError) then
                call mprintf(.true., iERROR, &
                            ' INQUIRE FAILED :: Put2d_ASCII ---> file not readable (filename: '//trim(sFileName)//')')
            else
                call mprintf(.true., iWARN, &
                            ' INQUIRE FAILED :: Put2d_ASCII ---> file not readable (filename: '//trim(sFileName)//')')
                iErr = iRet
                return
            endif
            
        elseif ( bFileExist ) then
            
            if (iRet == 0) then 
                
                call mprintf(.true., iINFO_Extra, &
                            ' WRITE OK :: Put2d_ASCII ---> file writable (filename: '//trim(sFileName)//')')
                
                !do iJ = 1,iY
                    do iI = 1,iX
                        if (present(sFMT)) then
                            write(20 , sFMT) a2dVar(iI, :) 
                        else
                            write(20 , *) a2dVar(iI, :) 
                        endif
                    enddo
                !enddo
                close(20)
                
            elseif (iRet /= 0) then
                
                if (bFatalError) then
                    call mprintf(.true., iERROR, &
                            ' WRITE FAILED :: Put2d_ASCII ---> file not writable (filename: '//trim(sFileName)//')')
                else
                    call mprintf(.true., iWARN, &
                            ' WRITE FAILED :: Put2d_ASCII ---> file not writable (filename: '//trim(sFileName)//')')
                    iErr = iRet
                    return
                endif
                   
            endif

        endif
        !------------------------------------------------------------------------------------
        
    end subroutine HMC_Tools_IO_Put2d_ASCII
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to get 2d variable binary
    subroutine HMC_Tools_IO_Get2d_Binary(sFileName, a2dVar, iX, iY, iVarScale, bFatalError, iErr)
        
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4) :: iX, iY, iI, iJ

        character(len = 256) :: sFileName
        
        real(kind = 4), dimension(iX, iY), intent(out)          :: a2dVar 
        
        integer(kind = 4), dimension(iX, iY)                    :: a2iVar 
        real(kind = 8), dimension(iX, iY)                       :: a2dVarPass 

        logical, intent(in) :: bFatalError
        logical             :: bFileExist
        integer(kind = 4)   :: iErr, iRet, iVarScale
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Initialized variable(s)
        a2iVar = 0; a2dVar = 0.0; a2dVarPass = 0.0
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Open and read binary file
        inquire (file = trim(sFileName), exist = bFileExist, iostat = iRet)
        if (.not. bFileExist ) then
            
            if (bFatalError) then
                call mprintf(.true., iERROR, &
                            ' INQUIRE FAILED :: Get2d_Binary ---> file not found (filename: '//trim(sFileName)//')')
            else
                call mprintf(.true., iWARN, &
                            ' INQUIRE FAILED :: Get2d_Binary ---> file not found (filename: '//trim(sFileName)//')')
                iErr = iRet
                return
            endif
            
        elseif ( bFileExist ) then
            
            open(unit=1, file=trim(sFileName), form='unformatted', access='direct', &
                 recl=iX*4, iostat=iRet)
                 
            if (iRet == 0) then 
                
                call mprintf(.true., iINFO_Extra, &
                            ' READ OK :: Get2d_Binary ---> file readable (filename: '//trim(sFileName)//')')
                
                a2iVar = 0
                do iJ = 1,iY
                   read(1,rec=iJ)(a2iVar(iX - iI + 1, iJ), iI = 1,iX)
                enddo
                close(1)
                
                a2dVarPass = real(a2iVar)/real(iVarScale)
                
            elseif (iRet /= 0) then
                
                if (bFatalError) then
                    call mprintf(.true., iERROR, &
                            ' READ FAILED :: Get2d_Binary ---> file not readable (filename: '//trim(sFileName)//')')
                else
                    call mprintf(.true., iWARN, &
                            ' READ FAILED :: Get2d_Binary ---> file not readable (filename: '//trim(sFileName)//')')
                    iErr = iRet
                    return
                endif
                   
            endif

        endif
        
        ! Convert from real*8 to real*4
        a2dVar = sngl(a2dVarPass)
        !------------------------------------------------------------------------------------
        
    end subroutine HMC_Tools_IO_Get2d_Binary
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to get 3d variable binary
    subroutine HMC_Tools_IO_Get3d_Binary(sFileName, a3dVar, iX, iY, iT, iVarScale, bFatalError, iErr)
        
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4) :: iX, iY, iT, iI, iJ, iK

        character(len = 256) :: sFileName
        real(kind = 4), dimension(iX, iY, iT), intent(out)          :: a3dVar 
        integer(kind = 4), dimension(iX, iY, iT)                    :: a3iVar 

        logical, intent(in) :: bFatalError
        logical             :: bFileExist
        integer(kind = 4)   :: iErr, iRet, iVarScale
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a3iVar = 0; a3dVar = 0.0
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Open and read binary file
        inquire (file = trim(sFileName), exist = bFileExist, iostat = iRet)
        if (.not. bFileExist ) then
            
            if (bFatalError) then
                call mprintf(.true., iERROR, &
                            ' INQUIRE FAILED :: Get3d_Binary ---> file not found (filename: '//trim(sFileName)//')')
            else
                call mprintf(.true., iWARN, &
                            ' INQUIRE FAILED :: Get3d_Binary ---> file not found (filename: '//trim(sFileName)//')')
                iErr = iRet
                return
            endif
            
        elseif ( bFileExist ) then
            
            open(unit=1, file=trim(sFileName), form='unformatted', access='direct', &
                 recl=iX*4, iostat=iRet)
            
            if (iRet == 0) then 
                
                call mprintf(.true., iINFO_Extra, &
                            ' READ OK :: Get3d_Binary ---> file readable (filename: '//trim(sFileName)//')')
                           
                do iK = 1,iT
                    do iJ = 1,iY
                       read(1,rec=iJ)(a3iVar(iX - iI + 1, iJ, iK),iI = 1,iX)
                    enddo
                enddo
                close(1)

                a3dVar = a3iVar/real(iVarScale)
                
            elseif (iRet /= 0) then
                
                if (bFatalError) then
                    call mprintf(.true., iERROR, &
                            ' READ FAILED :: Get3d_Binary ---> file not readable (filename: '//trim(sFileName)//')')
                else
                    call mprintf(.true., iWARN, &
                            ' READ FAILED :: Get3d_Binary ---> file not readable (filename: '//trim(sFileName)//')')
                    iErr = iRet
                    return
                endif
                   
            endif

        endif
        !------------------------------------------------------------------------------------

    end subroutine
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to put 2d variable binary
    subroutine HMC_Tools_IO_Put2d_Binary(sFileName, a2dVar, iX, iY, iVarScale, bFatalError, iErr)
    
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4) :: iX, iY, iI, iJ

        character(len = 256) :: sFileName
        real(kind = 4), dimension(iX, iY)                       :: a2dVar 
        integer(kind = 4), dimension(iX, iY)                    :: a2iVar 

        logical             :: bFileExist, bFatalError
        integer(kind = 4)   :: iErr, iRet, iVarScale
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a2iVar = 0;
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Open and write binary file
        open(20, file=trim(sFileName), form='unformatted', access='direct', recl=iX*4, iostat=iRet)
        
        ! Check file availability
        inquire (file = trim(sFileName), exist = bFileExist, iostat = iRet)
        if (.not. bFileExist ) then
            
            if (bFatalError) then
                call mprintf(.true., iERROR, &
                            ' INQUIRE FAILED :: Put2d_Binary ---> file not found (filename: '//trim(sFileName)//')')
            else
                call mprintf(.true., iWARN, &
                            ' INQUIRE FAILED :: Put2d_Binary ---> file not found (filename: '//trim(sFileName)//')')
                iErr = iRet
                return
            endif
            
        elseif ( bFileExist ) then
            
            if (iRet == 0) then 
                
                call mprintf(.true., iINFO_Extra, &
                            ' WRITE OK :: Put2d_Binary ---> file writable (filename: '//trim(sFileName)//')')
                
                ! Compute integer array variable
                a2iVar = int(a2dVar*real(iVarScale))
                ! Save data
                do iJ = 1,iY
                    write(20,rec=iJ)(a2iVar(iX - iI + 1,iJ),iI = 1,iX)
                enddo
                ! Close file
                close(20)
                
            elseif (iRet /= 0) then
                
                if (bFatalError) then
                    call mprintf(.true., iERROR, &
                            ' WRITE FAILED :: Put2d_Binary ---> file not writable (filename: '//trim(sFileName)//')')
                else
                    call mprintf(.true., iWARN, &
                            ' WRITE FAILED :: Put2d_Binary ---> file not writable (filename: '//trim(sFileName)//')')
                    iErr = iRet
                    return
                endif
                
            endif
            
        endif
        !------------------------------------------------------------------------------------
        
    end subroutine HMC_Tools_IO_Put2d_Binary
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to put 3d variable binary
    subroutine HMC_Tools_IO_Put3d_Binary(sFileName, a3dVar, iX, iY, iT, iVarScale, bFatalError, iErr)
    
        !------------------------------------------------------------------------------------
        ! Variable(s) declaration
        integer(kind = 4) :: iX, iY, iT, iI, iJ, iK

        character(len = 256) :: sFileName
        real(kind = 4), dimension(iX, iY, iT), intent(out)          :: a3dVar 
        integer(kind = 4), dimension(iX, iY, iT)                    :: a3iVar 

        logical             :: bFileExist, bFatalError
        integer(kind = 4)   :: iErr, iRet, iVarScale
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Initialize variable(s)
        a3iVar = 0;
        !------------------------------------------------------------------------------------
        
        !------------------------------------------------------------------------------------
        ! Open file
        open(20, file=trim(sFileName), form='unformatted', access='direct', recl=iX*4, iostat=iRet)
        
        ! Check file availability
        inquire (file = trim(sFileName), exist = bFileExist, iostat = iErr)
        if (.not. bFileExist ) then
            
            if (bFatalError) then
                call mprintf(.true., iERROR, &
                            ' INQUIRE FAILED :: Put3d_Binary ---> file not found (filename: '//trim(sFileName)//')')
            else
                call mprintf(.true., iWARN, &
                            ' INQUIRE FAILED :: Put3d_Binary ---> file not found (filename: '//trim(sFileName)//')')
                iErr = iRet
                return
            endif
            
        elseif ( bFileExist ) then
            
            if (iRet == 0) then 
                
                call mprintf(.true., iINFO_Extra, &
                            ' WRITE OK :: Put3d_Binary ---> file writable (filename: '//trim(sFileName)//')')
                
                ! Compute integer array variable
                a3iVar = int(a3dVar*real(iVarScale))
                ! Save data
                do iK = 1,iT
                    do iJ = 1,iY
                        write(20,rec=iJ)(a3iVar(iX - iI + 1, iJ, iK),iI = 1,iX)
                    enddo
                enddo
                ! Close file
                close(20)
                
            elseif (iRet /= 0) then
                
                if (bFatalError) then
                    call mprintf(.true., iERROR, &
                            ' WRITE FAILED :: Put3d_Binary ---> file not writable (filename: '//trim(sFileName)//')')
                else
                    call mprintf(.true., iWARN, &
                            ' WRITE FAILED :: Put3d_Binary ---> file not writable (filename: '//trim(sFileName)//')')
                    iErr = iRet
                    return
                endif
                
            endif
            
        endif
        !------------------------------------------------------------------------------------
        
    end subroutine HMC_Tools_IO_Put3d_Binary
    !------------------------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------
    ! Subroutine to check netCDF file error(s)
#ifdef LIB_NC
    subroutine check(status)

        integer, intent (in) :: status

        if(status /= nf90_noerr) then 
            call mprintf(.true., iERROR, ' netCDF error found: '//trim(nf90_strerror(status)) )
        end if

    end subroutine check 
#endif
    !------------------------------------------------------------------------------------
    
end module HMC_Module_Tools_IO
!------------------------------------------------------------------------------------