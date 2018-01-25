!==========================================================================================!
!     This module contains the Photosynthesis model.  The references are:                  !
!                                                                                          !
!       Xu et al. (2016) Diversity in plant hydraulic traits explains seasonal and         !
!           inter-annual variations in vegetation dynamics in seasonally dry tropical      !
!           forests.  New Phytologist                                                      !
!                                                                                          !
!       Manzoni, S., G. Vico, et al. (2011). "Optimizing stomatal conductance for maximum  !
!           carbon gain under water stress: a meta-analysis across plant functional types  !
!           and climates." Functional Ecology 25(3): 456-467.                              ! 
!                                                                                          !
!       Vico, G., S. Manzoni, et al. (2013). "A perspective on optimal leaf                !
!           stomatal conductance under CO2 and light co-limitations." Agricultural         !
!           and Forest Meteorology 182–183(0): 191-199.                                    !
!                                                                                          !
!       Katul, G., S. Manzoni, et al. (2010). "A stomatal optimization theory to describe  !
!           the effects of atmospheric CO2 on leaf photosynthesis and transpiration."      !
!           Annals of Botany 105(3): 431-442.                                              !
!                                                                                          !
!------------------------------------------------------------------------------------------!
Module farq_katul

Contains

  subroutine katul_lphys(can_prss,can_rhos,can_shv,can_co2,ipft,leaf_par,leaf_temp        &
                        ,lint_shv,green_leaf_factor,leaf_aging_factor,llspan,vm_bar,vm0in &
                        ,leaf_gbw,psi_leaf,last_gV,last_gJ,A_open,A_closed,gsw_open       &
                        ,gsw_closed,lsfc_shv_open,lsfc_shv_closed,lsfc_co2_open           &
                        ,lsfc_co2_closed,lint_co2_open,lint_co2_closed,leaf_resp          &
                        ,vmout,comppout,limit_flag)
      
    use rk4_coms       , only : tiny_offset              & ! intent(in)
                              , effarea_transp           ! ! intent(in)

    use pft_coms       , only : photosyn_pathway         & ! intent(in)
                              , phenology                & ! intent(in)
                              , Vm0                      & ! intent(in)
                              , Rd0                      & ! intent(in)
                              , cuticular_cond           & ! intent(in)
                              , dark_respiration_factor  & ! intent(in)
                              , quantum_efficiency       & ! intent(in)
                              , TLP                      & ! intent(in)
                              , stoma_lambda             & ! intent(in)
                              , stoma_beta               & ! intent(in)
                              , leaf_psi_min             &
                              , photosyn_pathway         ! ! intent(in)
    use consts_coms    , only : t00                      & ! intent(in)
                              , mmdry1000                & ! intent(in)
                              , mmh2o                    & ! intent(in)
                              , mmdry                    ! ! intent(in)
    use therm_lib      , only : eslf                     ! ! function
    use ed_misc_coms   , only : current_time             ! ! intent(in)
    use physiology_coms, only : gbw_2_gbc                & ! intent(in)
                              , gsw_2_gsc                & ! intent(in)
                              , track_plant_hydro        & ! intent(in)
                              , o2_ref                   ! ! intent(in)
    use phenology_coms , only : vm0_tran                 & ! intent(in)
                              , vm0_slope                & ! intent(in)
                              , vm0_amp                  & ! intent(in)
                              , vm0_min                  ! ! intent(in)

    implicit none

      !------ Arguments. ------------------------------------------------------------------!
      real(kind=4), intent(in)    :: can_prss          ! Canopy air pressure    [       Pa]
      real(kind=4), intent(in)    :: can_rhos          ! Canopy air density     [    kg/m³]
      real(kind=4), intent(in)    :: can_shv           ! Canopy air sp. hum.    [    kg/kg]
      real(kind=4), intent(in)    :: can_co2           ! Canopy air CO2         [ µmol/mol]
      integer     , intent(in)    :: ipft              ! Plant functional type  [      ---]
      real(kind=4), intent(in)    :: leaf_par          ! Absorbed PAR           [     W/m²]
      real(kind=4), intent(in)    :: leaf_temp         ! Leaf temperature       [        K]
      real(kind=4), intent(in)    :: lint_shv          ! Leaf interc. sp. hum.  [    kg/kg]
      real(kind=4), intent(in)    :: green_leaf_factor ! Frac. of on-allom. gr. [      ---]
      real(kind=4), intent(in)    :: leaf_aging_factor ! Ageing parameter       [      ---]
      real(kind=4), intent(in)    :: llspan            ! Leaf life span         [     mnth]
      real(kind=4), intent(in)    :: vm_bar            ! Average Vm function    [µmol/m²/s]
      real(kind=4), intent(in)    :: vm0in             ! Vm0 of the leaf        [µmol/m²/s]
      real(kind=4), intent(in)    :: leaf_gbw          ! B.lyr. cnd. of H2O     [  kg/m²/s]
      real(kind=4), intent(in)    :: psi_leaf          ! leaf water potential   [        m]
      real(kind=4), intent(inout) :: last_gV           ! gs for last timestep   [  kg/m²/s]
      real(kind=4), intent(inout) :: last_gJ           ! gs for last timestep   [  kg/m²/s]
      real(kind=4), intent(out)   :: A_open            ! Photosyn. rate (op.)   [µmol/m²/s]
      real(kind=4), intent(out)   :: A_closed          ! Photosyn. rate (cl.)   [µmol/m²/s]
      real(kind=4), intent(out)   :: gsw_open          ! St. cnd. of H2O  (op.) [  kg/m²/s]
      real(kind=4), intent(out)   :: gsw_closed        ! St. cnd. of H2O  (cl.) [  kg/m²/s]
      real(kind=4), intent(out)   :: lsfc_shv_open     ! Leaf sfc. sp.hum.(op.) [    kg/kg] 
      real(kind=4), intent(out)   :: lsfc_shv_closed   ! Leaf sfc. sp.hum.(cl.) [    kg/kg]
      real(kind=4), intent(out)   :: lsfc_co2_open     ! Leaf sfc. CO2    (op.) [ µmol/mol]
      real(kind=4), intent(out)   :: lsfc_co2_closed   ! Leaf sfc. CO2    (cl.) [ µmol/mol]
      real(kind=4), intent(out)   :: lint_co2_open     ! Intercell. CO2   (op.) [ µmol/mol]
      real(kind=4), intent(out)   :: lint_co2_closed   ! Intercell. CO2   (cl.) [ µmol/mol]
      real(kind=4), intent(out)   :: leaf_resp         ! Leaf respiration rate  [µmol/m²/s]
      real(kind=4), intent(out)   :: vmout             ! Max. Rubisco capacity  [µmol/m²/s]
      real(kind=4), intent(out)   :: comppout          ! GPP compensation point [ µmol/mol]
      integer     , intent(out)   :: limit_flag        ! Photosyn. limit. flag  [      ---]
      !----- External function. -----------------------------------------------------------!
      real(kind=4)    , external      :: sngloff     ! Safe double -> single precision
      !----- Local Variables    -----------------------------------------------------------!
      real(kind=4)                :: leaf_o2            ! O2 concentration  mmol/mol
      real(kind=4)                :: par                ! PAR          micromol/m2/s
      real(kind=4)                :: leaf_temp_degC     ! leaf_temperature in degree C
      real(kind=4)                :: can_vpr_prss       ! canopy vapor pressure in kPa
      real(kind=4)                :: temp_coef          ! temperautre coefficient
      real(kind=4)                :: Vcmax              ! current Vcmax  umol/m2/s
      real(kind=4)                :: Vcmax25            ! current Vcmax at 25 degC umol/m2/s
      real(kind=4)                :: Vcmax15            ! current Vcmax at 15 degC umol/m2/s
      real(kind=4)                :: Jmax               ! current Jmax  umol/m2/s
      real(kind=4)                :: Jmax25             ! current Jmax at 25 degC umol/m2/s
      real(kind=4)                :: Jrate              ! current Jrate umol/m2/s
      real(kind=4)                :: Rdark              ! current dark respiration rate umol/m2/s
      real(kind=4)                :: cuticular_gsc      ! current cuticular_conductance for CO2 umol/m2/s
      real(kind=4)                :: lambda             ! current lambda factor    numo/mol/kPa
      real(kind=4)                :: down_factor        ! photosynthetic down-regulation factor
      real(kind=4)                :: aero_resistance    ! aerodynamic resistance
      real(kind=4)                :: accepted_gsc       ! gsc solved from optimization scheme umol/m2/s
      real(kind=4)                :: accepted_fc        ! CO2 flux solved from optimization scheme umol/m2/s
      real(kind=4)                :: accepted_ci        ! Internal CO2 concentration ppm
      real(kind=4)                :: leaf_vpr_prss      ! leaf internal vapor pressure  in kPa
      real(kind=4)                :: cp                 ! constants in photosynthesis  umol/mol
      real(kind=4)                :: kc                 ! constants in photosynthesis  umol/mol
      real(kind=4)                :: ko                 ! constants in photosynthesis  umol/mol
      real(kind=4)                :: dfcdg              ! variables used in optimization scheme
      real(kind=4)                :: dfedg              ! variables used in optimization scheme
      real(kind=4)                :: d2fcdg2            ! variables used in optimization scheme
      real(kind=4)                :: d2fedg2            ! variables used in optimization scheme
      real(kind=4)                :: delta_g            ! variables used in optimization scheme
      real(kind=4)                :: a1gk,a2gk          ! variables used in optimization scheme
      real(kind=4)                :: k1ci,k2ci,k3ci,k4ci! variables used in optimization scheme
      real(kind=4)                :: resid,test_gsc     ! variables used in optimization scheme
      real(kind=4)                :: test_gV            ! variables used in optimization scheme
      real(kind=4)                :: test_fcV           ! variables used in optimization scheme
      real(kind=4)                :: test_ciV           ! variables used in optimization scheme
      real(kind=4)                :: test_gJ            ! variables used in optimization scheme
      real(kind=4)                :: test_fcJ           ! variables used in optimization scheme
      real(kind=4)                :: test_ciJ           ! variables used in optimization scheme
      integer                     :: iter               ! variables used in optimization scheme
      real(kind=4)                :: testfc1,testfc2    ! variables used in optimization scheme
      real(kind=4)                :: testci             ! variables used in optimization scheme
      integer                     :: k
      logical                     :: is_resolvable
      logical, parameter          :: quality_check = .false.
    
      !------------------     Make sure leaf psi is computed
      if (track_plant_hydro  == 0 ) then   
            write (unit=*,fmt='(80a)')         ('=',k=1,80)
            write (unit=*,fmt='(a)')           'Leaf water potential is tracked. Turn TRACK_PLANT_HYDRO to 1      '
            call fatal_error('Leaf water potential is not tracked          ','katul_lphys'               &
                            &,'farq_katul.f90')
      endif

      !------------------     Make sure only C3 plants go into this scheme
!      if (photosyn_pathway(ipft) == 4) then
!            write (unit=*,fmt='(80a)')         ('=',k=1,80)
!            write (unit=*,fmt='(a)')           'The optimization stomatal scheme for C4 plants are under development.....'
!            call fatal_error('Optimizing stomatal conductance for C4 plants','katul_lphys'               &
!                            &,'farq_katul.f90')
!      endif

      !-------------------    Define some constants....
      leaf_o2           = o2_ref * 1000.                             ! convert to mmol/mol
      par               = leaf_par * 4.6                             ! convert to micromol/m2/s
      leaf_temp_degC    = leaf_temp - t00                            ! convert to degC
      leaf_vpr_prss     = eslf(leaf_temp) * 0.001                    ! in kPa
      can_vpr_prss      = can_shv * can_prss * 0.001                 ! in kPa
  
    
      !------------------------------------------------------------------------------------!
      ! correcting for light-phenology
      select case(phenology(ipft))
      case (3)
         !------ Light-controlled phenology. ----------------------------------------------!
         Vcmax15    = vm0_amp / (1.0 + (llspan/vm0_tran)**vm0_slope) + vm0_min
      case default
         !------ Other phenologies, no distinction on Vm0. --------------------------------!
         Vcmax15    = vm0in
      end select
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      ! correcting for water stress
      Vcmax15       = Vcmax15

      cuticular_gsc = cuticular_cond(ipft) * 1.0e-6 * &
                    max(0.0,min(1.0, (psi_leaf - leaf_psi_min(ipft)) &
                    / (TLP(ipft) - leaf_psi_min(ipft))))
      ! when leaf water potential is too low, there can't be any water loss...
      ! mimicking the cease of soil evaporation near sh

      lambda =  stoma_lambda(ipft) * can_co2 / 400. * exp(stoma_beta(ipft) * psi_leaf)
       
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! Set temperature dependence and other derived variables
    
      ! Vcmax, Jmax Temperature dependence according to Farquhar

      ! CLM parameters for Rd
      ! Bonan et al. 2011

      call temp_fun(15.+273.15,298.15,&
                            58.52,  & ! Hv
                            0.710,   & ! Sv
                            220.0,  & ! Hd
                            temp_coef)
      Vcmax25 = Vcmax15 / temp_coef
      Jmax25 = Vcmax25 * 1.97

      ! calculate the Vcmax Jmax and Rd at current T
      call temp_fun(leaf_temp,298.15,&
                            58.52,  & ! Hv
                            0.710,   & ! Sv
                            220.0,  & ! Hd
                            temp_coef)
      Vcmax = Vcmax25 * temp_coef

      ! Jmax
      call temp_fun(leaf_temp,298.15,&
                            37.10,  & ! Hv
                            0.710,   & ! Sv
                            220.,  & ! Hd
                            temp_coef)
      Jmax = Jmax25 * temp_coef

      ! Rd
      call temp_fun(leaf_temp,298.15,&
                            64.5,  & ! Hv
                            0.71,  &!0.0,   & ! Sv
                            220.,  & ! Hd
                            temp_coef)
      Rdark = Vcmax25 * dark_respiration_factor(ipft) * temp_coef
    
      cp = exp(19.02-37830./(8.314* leaf_temp ))
      kc = exp(38.05-79430./(8.314* leaf_temp ))
      ko = exp(20.30-36380./(8.314* leaf_temp ))

      !------------------------------------------------------------------------------------!

  	  down_factor   = max(1e-6,min(1.0,1. / (1. + (psi_leaf / TLP(ipft)) ** 6.0)))
      Vcmax = Vcmax * down_factor
      Jmax = Jmax * down_factor

      ! Solve the quadratic function for light-limited photosynthesis
      if (photosyn_pathway(ipft) == 3) then
          Jrate = ((Jmax+0.385*par) -   &
                  sqrt((Jmax+0.385*par)**2-4.*0.7*0.385*Jmax*par))/1.4
      elseif (photosyn_pathway(ipft) == 4) then
          ! this is still under test
          Jrate = quantum_efficiency(ipft) * par  
      endif

      ! calcualte aerodynamic resistance
      if (leaf_gbw > 0.) then
          aero_resistance = mmdry / (leaf_gbw / gbw_2_gbc)
      else
          aero_resistance = 1e10
      endif

      is_resolvable = (Jmax /= 0.) .and. (Vcmax /= 0.) .and. &
            (aero_resistance < 1e8) .and.(cuticular_gsc > 1e-8)

      if (is_resolvable) then
        ! Rubisco-limited photosynthesis
        a1gk = Vcmax
        a2gk = kc * (1. + leaf_o2 / ko)
        k1ci = a1gk / can_co2 - Rdark / can_co2
        k2ci = a1gk * aero_resistance / can_co2 - Rdark * aero_resistance / can_co2 - 1. + a2gk / can_co2
        k3ci = -a1gk*cp/can_co2/can_co2-Rdark*a2gk/can_co2/can_co2
        k4ci = -a1gk*cp/can_co2*aero_resistance/can_co2-Rdark*a2gk/can_co2*aero_resistance/can_co2-a2gk/can_co2
       

        if(leaf_vpr_prss > can_vpr_prss) then
            ! Use newton's method to find the zero point of 
            ! start with gsw from last time
            if ((.not. isnan(last_gV)) .and. last_gV > 1e-10) then
               test_gsc = last_gV
            else
               test_gsc = cuticular_gsc
            endif
           
            do iter = 1, 500
                ! calculate dfcdg - dfedg
                call fluxsolver(test_gsc, aero_resistance, can_co2, k1ci, k2ci, k3ci, k4ci, testci, testfc1)
                call fcderiv(test_gsc, aero_resistance, k1ci, k2ci, k3ci, k4ci, can_co2, testci, dfcdg)
            
                dfedg = lambda / gsw_2_gsc * (leaf_vpr_prss - can_vpr_prss) / &
                       (1. + gbw_2_gbc/gsw_2_gsc * test_gsc * aero_resistance) ** 2

                call dfcdgderiv(test_gsc, aero_resistance, k1ci, k2ci, k3ci, k4ci,can_co2,testci,dfcdg,d2fcdg2)
                d2fedg2 = - 2. * gbw_2_gbc/gsw_2_gsc * aero_resistance / &
                          (gbw_2_gbc/gsw_2_gsc * test_gsc * aero_resistance + 1) * dfedg

                ! calculate the derivative of dfcdg - dfedg
                if (d2fcdg2 - d2fedg2 == 0) then
                    delta_g = 0.
                else
                    delta_g = - (dfcdg - dfedg) / (d2fcdg2 - d2fedg2)
                endif

                ! control exit

                if (abs(dfcdg - dfedg) < 1e-4 .or.                              &  ! converge
                    (test_gsc < cuticular_gsc .and. dfcdg-dfedg < 0.) .or.      &  ! close stomatal
                    (test_gsc < cuticular_gsc .and. isnan(dfcdg-dfedg)) .or.    &  ! close stomatal
                    (test_gsc + delta_g <= 0.)                 .or.             &  ! unrealistic values
                    (isnan(delta_g) < 0.)                      .or.             &  ! unrealistic values
                    (delta_g == 0.)                            .or.             &  ! trapped
                    (test_gsc > 1.0 .and. dfcdg-dfedg > 0.)                     &  ! fully open stomatal
                   ) then
                    exit
                endif

                test_gsc = test_gsc + delta_g

            enddo

            ! check the case that a negative or no optimal value is found
            if (test_gsc < cuticular_gsc .or. isnan(test_gsc)) then 
                test_gsc = cuticular_gsc
            endif
                
            if (test_gsc > 1.0) then 
               if (par > 0. .and. dfcdg-dfedg > 0.) then ! light
                    test_gsc = 1.0
               else ! dark
                    test_gsc = cuticular_gsc 
               endif
                
            endif

            call fluxsolver(test_gsc, aero_resistance, can_co2, k1ci, k2ci, k3ci, k4ci, testci, testfc1)

            last_gV = test_gsc
            
            test_gV = test_gsc
            test_fcV = testfc1
            test_ciV = testci
          
        else
          ! canopy air is saturated...
            test_ciV = can_co2
            test_gV = 1000. / aero_resistance
            call fluxsolver(test_gV, aero_resistance, can_co2, k1ci, k2ci, k3ci, k4ci, testci,test_fcV)
            
        endif


        ! Repeat calculation for light-limited case
        Jrate = Jrate * 0.25
        a1gk = Jrate
        a2gk = 2. * cp
        k1ci = a1gk / can_co2 - Rdark / can_co2
        k2ci = a1gk * aero_resistance / can_co2 - Rdark * aero_resistance / can_co2 - 1. + a2gk / can_co2
        k3ci = -a1gk*cp/can_co2/can_co2-Rdark*a2gk/can_co2/can_co2
        k4ci = -a1gk*cp/can_co2*aero_resistance/can_co2-Rdark*a2gk/can_co2*aero_resistance/can_co2-a2gk/can_co2
    
        if(leaf_vpr_prss > can_vpr_prss)then
        ! Use newton's method to find the zero point of 
        ! start with gsw from last time
            if ((.not. isnan(last_gJ)) .and. last_gJ > 1e-10) then
               test_gsc = last_gJ
            else
               test_gsc = cuticular_gsc 
            endif

            do iter = 1, 500
                ! calculate dfcdg - dfedg
                call fluxsolver(test_gsc, aero_resistance, can_co2, k1ci, k2ci, k3ci, k4ci, testci, testfc1)
                call fcderiv(test_gsc, aero_resistance, k1ci, k2ci, k3ci, k4ci, can_co2, testci, dfcdg)
            
                dfedg = lambda / gsw_2_gsc * (leaf_vpr_prss - can_vpr_prss) / &
                        (1. + gbw_2_gbc/gsw_2_gsc * test_gsc * aero_resistance) ** 2
                

                call dfcdgderiv(test_gsc, aero_resistance, k1ci, k2ci, k3ci, k4ci,can_co2,testci,dfcdg,d2fcdg2)
                d2fedg2 = -2. * gbw_2_gbc/gsw_2_gsc * aero_resistance / &
                          (gbw_2_gbc/gsw_2_gsc * test_gsc * aero_resistance + 1) * dfedg
    
                ! calculate the derivative of dfcdg - dfedg
                if (d2fcdg2 - d2fedg2 == 0.) then
                    delta_g = 0.
                else
                    delta_g = - (dfcdg - dfedg) / (d2fcdg2 - d2fedg2)
                endif

                ! control exit
                if (abs(dfcdg - dfedg) < 1e-4 .or.                              &   ! converge
                    (test_gsc < cuticular_gsc .and. dfcdg-dfedg < 0.) .or.      &   ! close stomatal
                    (test_gsc < cuticular_gsc .and. isnan(dfcdg-dfedg)) .or.    &   ! close stomatal
                    (test_gsc + delta_g <= 0.)                 .or.             &   ! unrealistic values
                    (isnan(delta_g) < 0.)                   .or.                &   ! unrealistic values
                    (delta_g == 0.)                 .or.                        &   ! trapped
                    (test_gsc > 1.0 .and. dfcdg-dfedg > 0.)                     &   ! fullyopen stomatal
                   ) then
                            exit
                endif

                test_gsc = test_gsc + delta_g

            enddo

            if (test_gsc < cuticular_gsc .or. isnan(test_gsc)) then
                test_gsc = cuticular_gsc
            endif
            
            if (test_gsc > 1.0) then 
                if (par > 0. .and. dfcdg-dfedg>0.) then ! light
                    test_gsc = 1.0
                else ! dark
                    test_gsc = cuticular_gsc 
                endif
            endif

            call fluxsolver(test_gsc, aero_resistance, can_co2, k1ci, k2ci, k3ci, k4ci, testci, testfc1)
            last_gJ = test_gsc
            test_gJ = test_gsc
            test_fcJ = testfc1
            test_ciJ = testci

        else
            test_ciJ = can_co2
            test_gJ = 1000. / aero_resistance
            call fluxsolver(test_gJ, aero_resistance, can_co2, k1ci, k2ci, k3ci, k4ci, testci,test_fcJ)
        endif


        if(test_fcV < test_fcJ)then
            accepted_fc  = test_fcV
            accepted_gsc = test_gV
            accepted_ci  = test_ciV
        else
           accepted_fc   = test_fcJ
           accepted_gsc  = test_gJ
           accepted_ci   = test_ciJ
        endif

    else  ! not resolvable
        accepted_gsc     = cuticular_gsc
        accepted_fc      = -Rdark
        accepted_ci      = can_co2

    endif

    if (quality_check) then
       write (unit=*,fmt='(80a)')         ('=',k=1,80)
       write (unit=*,fmt='(a)')           'Katul Stomatal Scheme Quality Check:'
       write (unit=*,fmt='(a,1x,i9)')   ' + HOUR:                ',current_time%hour
       write (unit=*,fmt='(a,1x,i9)')   ' + PFT:                 ',ipft
       write (unit=*,fmt='(a,1x,l9)')   ' + RESOLVABLE:          ',is_resolvable
       write (unit=*,fmt='(a,1x,es12.4)')   ' + PSI_LEAF:            ',psi_leaf
       write (unit=*,fmt='(a,1x,es12.4)')   ' + PAR:                 ',par
       write (unit=*,fmt='(a,1x,es12.4)')   ' + Vcmax25:             ',Vcmax25
       write (unit=*,fmt='(a,1x,es12.4)')   ' + Jmax25:              ',Jmax25
       write (unit=*,fmt='(a,1x,es12.4)')   ' + test_gV:             ',test_gV
       write (unit=*,fmt='(a,1x,es12.4)')   ' + test_gJ:             ',test_gJ
       write (unit=*,fmt='(a,1x,es12.4)')   ' + test_fcV:            ',test_fcV
       write (unit=*,fmt='(a,1x,es12.4)')   ' + test_fcJ:            ',test_fcJ
       write (unit=*,fmt='(a,1x,es12.4)')   ' + aero_resistance      ',aero_resistance
       write (unit=*,fmt='(a,1x,es12.4)')   ' + cuticular_gsc        ',cuticular_gsc
    endif
!-------------------- copy the solution to output-----------------------!
    A_closed        = -Rdark                ! umol/m2/s
    A_open          = accepted_fc           ! umol/m2/s
    leaf_resp       = Rdark

    gsw_closed      = cuticular_gsc / gsw_2_gsc  * mmdry / effarea_transp(ipft)  ! convert to kg/m2/s
    gsw_open        = accepted_gsc / gsw_2_gsc  * mmdry / effarea_transp(ipft)  ! convert to kg/m2/s

    !------------------- these variables are not tracked....
    lsfc_shv_closed = can_shv
    lsfc_shv_open   = can_shv

    lsfc_co2_open   = can_co2
    lsfc_co2_closed = can_co2

    lint_co2_closed = accepted_ci
    lint_co2_open   = accepted_ci

    vmout           = Vcmax
    comppout        = 0.
    limit_flag      = 0.

  end subroutine katul_lphys

!==========================================================

  
  subroutine fluxsolver(g, ra, ca, k1, k2, k3, k4, ci,fc)
    implicit none
    
    real, intent(in) :: g, k1,k2, k3,k4, ra, ca
    real, intent(out) :: ci, fc
    real :: cip, cim, rad
    
    if (g == 0) print*,'g is 0 in fluxsolver'
    rad = sqrt((k1/g+k2)**2 - 4. * (k3/g + k4))
    
    cip = ca * (-(k1/g+k2) + rad)/2.
    cim = ca * (-(k1/g+k2) - rad)/2.
    
    ci = cip
    
    fc = (ca - ci) / (1./g + ra)
    
    return
  end subroutine fluxsolver
  
  subroutine fcderiv(g, ra, k1, k2, k3, k4, ca, ci, dfcdg)
    implicit none
    
    real, intent(in) :: g, ra, k1, k2, k3, k4, ca, ci
    real, intent(out) :: dfcdg
    real :: dcidg, myroot
    
    myroot = sqrt((k1/g+k2)**2-4.*(k3/g+k4))
    
    dcidg = ca * (0.5*k1/g**2 +   &
         0.25/myroot*(-2.*k1**2/g**3-2.*k1*k2/g**2+4.*k3/g**2))
    
    dfcdg = ((1./g+ra)*(-dcidg) - (ca-ci)*(-1./g**2))/(1./g+ra)**2
    
    return
  end subroutine fcderiv

  subroutine dfcdgderiv(g, ra, k1, k2, k3, k4,ca,ci,dfcdg,d2fcdg2)
    implicit none
    real, intent(in) :: g,ra,k1,k2,k3,k4,ca,ci,dfcdg
    real, intent(out) :: d2fcdg2
    real :: t1, t2

    t1 = 4. * k3 / g - 2. * k1 * (k2 + k1/g) / g
    t2 = sqrt((k1/g + k2) ** 2 - 4. * (k3/g + k4))

    d2fcdg2 = 1. / ((ra + 1./g) * g ** 2) * ( &
                ca / 2. * (t1 ** 2 / (4. * t2 ** 3) - &
                (k1 ** 2 / g ** 2 - t1) / t2 + 2. * k1 * g) + &
                2. * dfcdg - 2. * (ca - ci) / (g * (ra + 1./g)))

  end subroutine dfcdgderiv
  !======================================================

  subroutine temp_fun(Tleaf,Tref,Hv,Sv,Hd,T_coef) ! based on Harley et al. 1991
  ! does not consider low temperature cut-off
  implicit none
  real, intent(in) :: Tleaf  ! K
  real, intent(in) :: Tref   ! K
  real, intent(in) :: Hv     ! kJ/mol
  real, intent(in) :: Sv     ! kJ/mol/K
  real, intent(in) :: Hd     ! kJ/mol
  real, intent(out) :: T_coef  ! unitless

  real,  parameter :: R = 8.314e-3 !kJ/mol/K
  ! local
  real             :: temp_T

  temp_T = Tleaf
  if (temp_T == 0.) temp_T = 0.1

  T_coef = exp(Hv/(R * Tref) * (1 - Tref/temp_T)) / &
           (1 + exp((Sv * temp_T - Hd)/ (R * temp_T)))

  end subroutine temp_fun

end Module farq_katul
