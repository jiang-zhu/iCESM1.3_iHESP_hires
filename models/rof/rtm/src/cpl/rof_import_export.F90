module rof_import_export

  use shr_kind_mod     , only : r8 => shr_kind_r8
  use shr_file_mod     , only : shr_file_setLogUnit, shr_file_setLogLevel, &
                                shr_file_getLogUnit, shr_file_getLogLevel, &
                                shr_file_getUnit, shr_file_setIO
  use shr_sys_mod      , only : shr_sys_abort
  use shr_const_mod    , only : SHR_CONST_REARTH
  use RunoffMod        , only : runoff
  use RtmVar           , only : ice_runoff, iulog, wiso_runoff
  use RtmSpmd          , only : masterproc
  use perf_mod         , only : t_startf, t_stopf, t_barrierf
  use rtm_cpl_indices  , only : nt_rtm, rtm_tracers 
  use rtm_cpl_indices  , only : index_x2r_Flrl_rofl,index_x2r_Flrl_rofi 
  use rtm_cpl_indices  , only : index_r2x_Forr_rofl, index_r2x_Forr_rofi
  use rtm_cpl_indices  , only : index_r2x_Flrr_flood, index_r2x_Flrr_volr
  use rtm_cpl_indices  , only : index_x2r_Flrl_rofl_16O, index_x2r_Flrl_rofi_16O, &
                                index_x2r_Flrl_rofl_18O, index_x2r_Flrl_rofi_18O, &
                                index_x2r_Flrl_rofl_HDO, index_x2r_Flrl_rofi_HDO, &
                                index_r2x_Forr_rofl_16O, index_r2x_Forr_rofi_16O, &
                                index_r2x_Forr_rofl_18O, index_r2x_Forr_rofi_18O, &
                                index_r2x_Forr_rofl_HDO, index_r2x_Forr_rofi_HDO, &
                                index_r2x_Flrr_flood_16O,index_r2x_Flrr_volr_16O, &
                                index_r2x_Flrr_flood_18O,index_r2x_Flrr_volr_18O, &
                                index_r2x_Flrr_flood_HDO,index_r2x_Flrr_volr_HDO

  implicit none
  public

contains

  subroutine rof_import( x2r, totrunin)

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Obtain the runoff input from the coupler
    !
    ! ARGUMENTS:
    real(r8), intent(inout) :: x2r(:,:)         
    real(r8), intent(out)   :: totrunin( runoff%begr: ,: ) 
    !
    ! LOCAL VARIABLES
    integer :: n2, n, nt
    integer :: begr, endr
    integer :: nliq, nfrz
    character(len=32), parameter :: sub = 'rof_import'
    !---------------------------------------------------------------------------
    
    ! Note that totrunin is a flux

    nliq = 0
    nfrz = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') then
          nliq = nt
       endif
       if (trim(rtm_tracers(nt)) == 'ICE') then
          nfrz = nt
       endif
    enddo
    if (nliq == 0 .or. nfrz == 0) then
       write(iulog,*) trim(sub),': ERROR in rtm_tracers LIQ ICE ',nliq,nfrz,rtm_tracers
       call shr_sys_abort()
    endif

    begr = runoff%begr
    endr = runoff%endr
    do n = begr,endr
       n2 = n - begr + 1
       totrunin(n,nliq) = x2r(index_x2r_Flrl_rofl,n2)
       totrunin(n,nfrz) = x2r(index_x2r_Flrl_rofi,n2)

       ! adding water isotopes
       if ( wiso_runoff ) then
          totrunin(n,3) = x2r(index_x2r_Flrl_rofl_16O,n2)
          totrunin(n,4) = x2r(index_x2r_Flrl_rofi_16O,n2)
          totrunin(n,5) = x2r(index_x2r_Flrl_rofl_18O,n2)
          totrunin(n,6) = x2r(index_x2r_Flrl_rofi_18O,n2)
          totrunin(n,7) = x2r(index_x2r_Flrl_rofl_HDO,n2)
          totrunin(n,8) = x2r(index_x2r_Flrl_rofi_HDO,n2)
       end if

    enddo

  end subroutine rof_import

  !====================================================================================

  subroutine rof_export(r2x)

    !---------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Send the runoff model export state to the coupler
    !
    ! ARGUMENTS:
    real(r8), intent(inout) :: r2x(:,:)  ! Runoff to coupler export state
    !
    ! LOCAL VARIABLES
    integer :: ni, n, nt
    integer :: nliq, nfrz
    logical :: first_time = .true.
    character(len=32), parameter :: sub = 'rof_export'
    !---------------------------------------------------------------------------
    
    nliq = 0
    nfrz = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') then
          nliq = nt
       endif
       if (trim(rtm_tracers(nt)) == 'ICE') then
          nfrz = nt
       endif
    enddo
    if (nliq == 0 .or. nfrz == 0) then
       write(iulog,*) trim(sub),': ERROR in rtm_tracers LIQ ICE ',nliq,nfrz,rtm_tracers
       call shr_sys_abort()
    endif

    r2x(:,:) = 0._r8

    if (first_time) then
       if (masterproc) then
          if ( ice_runoff )then
             write(iulog,*)'Snow capping will flow out in frozen river runoff'
          else
             write(iulog,*)'Snow capping will flow out in liquid river runoff'
          endif
       endif
       first_time = .false.
    end if

    ni = 0
    if ( ice_runoff )then
       do n = runoff%begr,runoff%endr
          ni = ni + 1
          if (runoff%mask(n) == 2) then
             ! liquid and ice runoff are treated separately - this is what goes to the ocean
             r2x(index_r2x_Forr_rofl,ni) = &
                  runoff%runoff(n,nliq)/(runoff%area(n)*1.0e-6_r8*1000._r8)
             r2x(index_r2x_Forr_rofi,ni) = &
                  runoff%runoff(n,nfrz)/(runoff%area(n)*1.0e-6_r8*1000._r8)

             ! adding water isotopes
             if ( wiso_runoff ) then
                r2x(index_r2x_Forr_rofl_16O,ni) = &
                   runoff%runoff(n,3)/(runoff%area(n)*1.0e-6_r8*1000._r8)
                r2x(index_r2x_Forr_rofi_16O,ni) = &
                   runoff%runoff(n,4)/(runoff%area(n)*1.0e-6_r8*1000._r8)

                r2x(index_r2x_Forr_rofl_18O,ni) = &
                   runoff%runoff(n,5)/(runoff%area(n)*1.0e-6_r8*1000._r8)
                r2x(index_r2x_Forr_rofi_18O,ni) = &
                   runoff%runoff(n,6)/(runoff%area(n)*1.0e-6_r8*1000._r8)

                r2x(index_r2x_Forr_rofl_HDO,ni) = &
                   runoff%runoff(n,7)/(runoff%area(n)*1.0e-6_r8*1000._r8)
                r2x(index_r2x_Forr_rofi_HDO,ni) = &
                   runoff%runoff(n,8)/(runoff%area(n)*1.0e-6_r8*1000._r8)
             end if

             if (ni > runoff%lnumr) then
                write(iulog,*) sub, ' : ERROR runoff count',n,ni
                call shr_sys_abort( sub//' : ERROR runoff > expected' )
             endif
          endif
       end do
    else
       do n = runoff%begr,runoff%endr
          ni = ni + 1
          if (runoff%mask(n) == 2) then
             ! liquid and ice runoff are bundled together to liquid runoff
             ! and then ice runoff set to zero
             r2x(index_r2x_Forr_rofl,ni) =   &
                  (runoff%runoff(n,nfrz)+runoff%runoff(n,nliq))&
                 /(runoff%area(n)*1.0e-6_r8*1000._r8)
             r2x(index_r2x_Forr_rofi,ni) = 0._r8

             ! adding water isotopes
             if ( wiso_runoff ) then
                r2x(index_r2x_Forr_rofl_16O,ni) =   &
                   (runoff%runoff(n,3)+runoff%runoff(n,4))/(runoff%area(n)*1.0e-6_r8*1000._r8)
                r2x(index_r2x_Forr_rofi_16O,ni) = 0._r8

                r2x(index_r2x_Forr_rofl_18O,ni) =   &
                   (runoff%runoff(n,5)+runoff%runoff(n,6))/(runoff%area(n)*1.0e-6_r8*1000._r8)
                r2x(index_r2x_Forr_rofi_18O,ni) = 0._r8

                r2x(index_r2x_Forr_rofl_HDO,ni) =   &
                   (runoff%runoff(n,7)+runoff%runoff(n,8))/(runoff%area(n)*1.0e-6_r8*1000._r8)
                r2x(index_r2x_Forr_rofi_HDO,ni) = 0._r8
             end if

             if (ni > runoff%lnumr) then
                write(iulog,*) sub, ' : ERROR runoff count',n,ni
                call shr_sys_abort( sub//' : ERROR runoff > expected' )
             endif
          endif
       end do
    end if

    ! Flooding back to land, sign convention is positive in land->rof direction
    ! so if water is sent from rof to land, the flux must be negative.
    ni = 0
    do n = runoff%begr, runoff%endr
       ni = ni + 1
       r2x(index_r2x_Flrr_flood,ni) = -runoff%flood(n)

       ! adding water isotopes
       ! if ( wiso_runoff ) then
       !    r2x(index_r2x_Flrr_flood,ni) = -runoff%flood(n,1)
       !    r2x(index_r2x_Flrr_flood_16O,ni) = -runoff%flood(n,2)
       !    r2x(index_r2x_Flrr_flood_18O,ni) = -runoff%flood(n,3)
       !    r2x(index_r2x_Flrr_flood_HDO,ni) = -runoff%flood(n,4)
       ! end if

    end do

    ! Want volr on land side to do a correct water balance
    ni = 0
    do n = runoff%begr, runoff%endr
      ni = ni + 1
         r2x(index_r2x_Flrr_volr,ni) = runoff%volr_nt1(n) / (runoff%area(n))

       ! adding water isotopes
       if ( wiso_runoff ) then
          r2x(index_r2x_Flrr_volr_16O,ni) = runoff%volr_nt3(n) / (runoff%area(n))
          r2x(index_r2x_Flrr_volr_18O,ni) = runoff%volr_nt5(n) / (runoff%area(n))
          r2x(index_r2x_Flrr_volr_HDO,ni) = runoff%volr_nt7(n) / (runoff%area(n))
       end if

    end do

  end subroutine rof_export

end module rof_import_export
