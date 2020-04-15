!==========================================================================================!
!==========================================================================================!
module mortality
   !=======================================================================================!
   !=======================================================================================!


   contains



   !=======================================================================================!
   !=======================================================================================!
   !    This subroutine computes the total PFT-dependent mortality rate:                   !
   !---------------------------------------------------------------------------------------!
   subroutine mortality_rates(cpatch,ico,avg_daily_temp, patch_age)
      use ed_state_vars , only : patchtype                  ! ! Structure
      use pft_coms      , only : mort0                      & ! intent(in)
                               , mort1                      & ! intent(in)
                               , mort2                      & ! intent(in)
                               , mort3                      & ! intent(in)
                               , mort_alpha                 & ! intent(in)
                               , mort_beta                  & ! intent(in)
                               , mort_plc_max               & ! intent(in)
                               , mort_plc_th                & ! intent(in)
                               , plant_min_temp             & ! intent(in)
                               , frost_mort                 ! ! intent(in)
      use disturb_coms  , only : treefall_disturbance_rate  & ! intent(in)
                               , treefall_hite_threshold    & ! intent(in)
                               , time2canopy                ! ! intent(in)
      use ed_max_dims   , only : n_pft                      & ! intent(in)
                               , n_mort                     ! ! intent(in)
      use consts_coms   , only : lnexp_min                  & ! intent(in)
                               , lnexp_max                  ! ! intent(in)
      use physiology_coms,only : imort_scheme               ! ! intent(in)
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(patchtype), target     :: cpatch          ! Current patch
      integer        , intent(in) :: ico             ! Current cohort ID
      real           , intent(in) :: avg_daily_temp  ! Mean temperature yesterday
      real           , intent(in) :: patch_age       ! Patch age
      !----- Local variables --------------------------------------------------------------!
      integer                     :: ipft            ! PFT 
      integer                     :: imonth          ! Month
      real                        :: temp_dep        ! Temp. function  (frost mortality)
      real                        :: expmort         ! Carbon-balance term
      real                        :: plc_avg         ! average percentage loss of conductance
      real                        :: plc_mort        ! plc mortality at monthly scale
      real                        :: ddbh_avg        ! average DBH growth rates              
      !------------------------------------------------------------------------------------!


      !----- Assume happy end, all plants survive... --------------------------------------!
      cpatch%mort_rate(1:4,ico) = 0.0
      ! don't include disturbance mortality
      cpatch%mort_rate(6:7,ico) = 0.0
      ipft = cpatch%pft(ico)

      !------------------------------------------------------------------------------------!
      ! 1.  Ageing, PFT-dependent but otherwise constant.                                  !
      !------------------------------------------------------------------------------------!
      cpatch%mort_rate(1,ico) = mort3(ipft)
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! 2.  Mortality rates due to negative carbon balance.                                !
      !------------------------------------------------------------------------------------!
      expmort = max( lnexp_min, min( lnexp_max                                             &
                                   , mort2(ipft) * ( cpatch%cbr_bar(ico) - mort0(ipft) ) ) )
      cpatch%mort_rate(2,ico) = mort1(ipft) / (1. + exp(expmort))
      !------------------------------------------------------------------------------------!


      !------------------------------------------------------------------------------------!
      ! 3.  Mortality due to treefall.                                                     !
      !------------------------------------------------------------------------------------!
      if (cpatch%hite(ico) <= treefall_hite_threshold .and. patch_age > time2canopy) then
         cpatch%mort_rate(3,ico) = treefall_disturbance_rate
      else
         cpatch%mort_rate(3,ico) = 0.
      end if
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! 4.   Mortality due to cold, after:                                                 !
      !      Albani, M.; D. Medvigy; G. C. Hurtt; P. R. Moorcroft, 2006: The contributions !
      !           of land-use change, CO2 fertilization, and climate variability to the    !
      !           Eastern US carbon sink.  Glob. Change Biol., 12, 2370-2390,              !
      !           doi: 10.1111/j.1365-2486.2006.01254.x                                    !
      !------------------------------------------------------------------------------------!
      temp_dep = max( 0.0, min( 1.0, 1.0 - (avg_daily_temp - plant_min_temp(ipft)) / 5.0) )
      cpatch%mort_rate(4,ico) = frost_mort(ipft) * temp_dep
      !------------------------------------------------------------------------------------!



      !------------------------------------------------------------------------------------!
      ! 5. Disturbance rate mortality.  This is not used by the cohort dynamics, instead   !
      !    this is just to account for the lost density due to the patch creation.  This   !
      !    mortality will be determined by the disturbance_mortality subroutine, not here. !
      !------------------------------------------------------------------------------------!
      !cpatch%mort_rate(5,ico) = TBD
      !------------------------------------------------------------------------------------!

      select case (imort_scheme)
      case (1)
      !------------------------------------------------------------------------------------!
      ! 6. Hydraulic failure mortality.                                                    !
      !------------------------------------------------------------------------------------!
      plc_mort = 0.
      do imonth = 1,12
        plc_mort = plc_mort + hydro_mort_rate(cpatch%plc_monthly(imonth,ico),ipft) / 12.
      enddo
      cpatch%mort_rate(6,ico) = plc_mort
      !------------------------------------------------------------------------------------!
      case (2)
      !------------------------------------------------------------------------------------!
      ! 7. Growth-based mortality from Camac et al. BioRxiv.
      ! Note we need to zero the negative carbon mortality in this 
      ! case to avoid double counting                                                      !
      !------------------------------------------------------------------------------------!
      ! calculate ddbh_avg
      ddbh_avg = sum(cpatch%ddbh_monthly(1:12,ico)) / 12.
        
      cpatch%mort_rate(2,ico) = 0. ! reset the original negative carbon mortality
      expmort = max( lnexp_min, min( lnexp_max                                             &
                                   , mort_beta(ipft) * ddbh_avg  ) )
      cpatch%mort_rate(7,ico) = mort_alpha(ipft) * exp(expmort)
      !------------------------------------------------------------------------------------!

      case (3)
      ! Both mortality is considered

      plc_mort = 0.
      do imonth = 1,12
        plc_mort = plc_mort + hydro_mort_rate(cpatch%plc_monthly(imonth,ico),ipft) / 12.
      enddo
      cpatch%mort_rate(6,ico) = plc_mort

      cpatch%mort_rate(2,ico) = 0. ! reset the original negative carbon mortality
      ! calculate ddbh_avg
      ddbh_avg = sum(cpatch%ddbh_monthly(1:12,ico)) / 12.
        
      expmort = max( lnexp_min, min( lnexp_max                                             &
                                   , mort_beta(ipft) * ddbh_avg  ) )
      cpatch%mort_rate(7,ico) = mort_alpha(ipft) * exp(expmort)

      end select

      return
   end subroutine mortality_rates
   !=======================================================================================!
   !=======================================================================================!






   !=======================================================================================!
   !=======================================================================================!
   !     This subroutine determines the mortality rates associated with the current        !
   ! disturbance.                                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine disturbance_mortality(csite,ipa,disturbance_rate,new_lu,dist_path            &
                                   ,mindbh_harvest)
      use ed_state_vars, only : sitetype  & ! structure
                              , patchtype ! ! structure
      use ed_max_dims  , only : n_pft     ! ! intent(in)
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(sitetype)                   , target     :: csite
      integer                          , intent(in) :: ipa
      real                             , intent(in) :: disturbance_rate
      integer                          , intent(in) :: new_lu
      integer                          , intent(in) :: dist_path
      real           , dimension(n_pft), intent(in) :: mindbh_harvest
      !----- Local variables. -------------------------------------------------------------!
      type(patchtype)                  , pointer    :: cpatch
      integer                                       :: ico
      real                                          :: f_survival
      !------------------------------------------------------------------------------------!

      cpatch => csite%patch(ipa)
      do ico=1,cpatch%ncohorts
         f_survival = survivorship(new_lu,dist_path,mindbh_harvest,cpatch,ico)
         cpatch%mort_rate(5,ico) = cpatch%mort_rate(5,ico)                                 &
                                 - log( f_survival                                         &
                                      + (1.0 - f_survival) * exp(- disturbance_rate) )
      end do
      return
   end subroutine disturbance_mortality
   !=======================================================================================!
   !=======================================================================================!





   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the survivorship rate associated with a disturbance.       !
   !  Input variables:                                                                     !
   !  -- new_lu: the disturbance/land use type after disturbance:                          !
   !     1. Clear cut (cropland and pasture).                                              !
   !     2. Forest plantation.                                                             !
   !     3. Tree fall.                                                                     !
   !     4. Fire.                                                                          !
   !     5. Forest regrowth.                                                               !
   !     6. Logged forest.                                                                 !
   !  -- dist_path: the pathway for the disturbance.  The flags depend on new_lu.  See     !
   !        comments at the select case (new_lu) block for additional details.             !
   !  -- mindbh_harvest: minimum DBH for harvesting (selective logging and forest          !
   !        plantantions).  If the tree DBH is greater than mindbh_harvest, the tree may   !
   !        be harvested, otherwise it may be damaged by logging but not harvested.        !
   !  -- cpatch: current patch.                                                            !
   !  -- ico: index for current cohort.                                                    !
   !---------------------------------------------------------------------------------------!
   real function survivorship(new_lu,dist_path,mindbh_harvest,cpatch,ico)
      use ed_state_vars, only : patchtype                ! ! structure
      use disturb_coms , only : treefall_hite_threshold  & ! intent(in)
                              , fire_hite_threshold      ! ! intent(in)
      use pft_coms     , only : treefall_s_ltht          & ! intent(in)
                              , treefall_s_gtht          & ! intent(in)
                              , fire_s_ltht              & ! intent(in)
                              , fire_s_gtht              ! ! intent(in)
      use ed_max_dims  , only : n_pft                    ! ! intent(in)
      
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      type(patchtype)                 , target     :: cpatch
      real          , dimension(n_pft), intent(in) :: mindbh_harvest
      integer                         , intent(in) :: ico
      integer                         , intent(in) :: new_lu
      integer                         , intent(in) :: dist_path
      !----- Local variables. -------------------------------------------------------------!
      integer                                      :: ipft
      !------------------------------------------------------------------------------------!


      !----- Alias for current PFT. -------------------------------------------------------!
      ipft = cpatch%pft(ico)
      !------------------------------------------------------------------------------------!

      !------------------------------------------------------------------------------------!
      !      Survivorship depends on the disturbance type (new_lu), and within each        !
      ! disturbance type, it may depend on the disturbance pathway (dist_path), the PFT,   !
      ! and size.                                                                          !
      !------------------------------------------------------------------------------------!
      select case(new_lu)
      case (1)
         !---------------------------------------------------------------------------------!
         !     Clear cut (cropland/pasture).  For now, nothing survives.                   !
         !---------------------------------------------------------------------------------!
         survivorship = 0.0
         !---------------------------------------------------------------------------------!

      case (2)
         !---------------------------------------------------------------------------------!
         !     Forest plantation.  Two types of mortality may exist: treefall disturbance  !
         ! rates (which will maintain a plantation set as a plantation as it is a managed  !
         ! land), and harvesting.                                                          !
         !---------------------------------------------------------------------------------!
         select case (dist_path)
         case (20)
            !----- Harvesting, assumes that nothing survives. -----------------------------!
            survivorship = 0.0
            !------------------------------------------------------------------------------!

         case (21)
            !----- Tree fall, assumes typical tree fall mortality. ------------------------!
            if (cpatch%hite(ico) < treefall_hite_threshold) then
               survivorship = treefall_s_ltht(ipft)
            else
               survivorship = treefall_s_gtht(ipft)
            end if
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!

      case (3)
         !---------------------------------------------------------------------------------!
         !     Tree fall.  Mortality depends on the cohort height and PFT.                 !
         !---------------------------------------------------------------------------------!
         if (cpatch%hite(ico) < treefall_hite_threshold) then
            survivorship = treefall_s_ltht(ipft)
         else
            survivorship = treefall_s_gtht(ipft)
         end if
         !---------------------------------------------------------------------------------!

      case (4)
         !---------------------------------------------------------------------------------!
         !     Fire.  For now fire mortality depends on the cohort height and PFT.         !
         !---------------------------------------------------------------------------------!
         if (cpatch%hite(ico) < fire_hite_threshold) then
            survivorship = fire_s_ltht(ipft)
         else
            survivorship = fire_s_gtht(ipft)
         end if
         !---------------------------------------------------------------------------------!

       case (5)
         !---------------------------------------------------------------------------------!
         !     Abandonment (secondary regrowth).  Two paths are possible: abandonment      !
         ! occurs after one last harvest, or the field/plantation is left as is.           !
         !---------------------------------------------------------------------------------!
         select case (dist_path)
         case (50)
            !----- Agriculture field is left as is. ---------------------------------------!
            survivorship = 1.0
            !------------------------------------------------------------------------------!
         case (51)
            !----- Forest plantation abandoned following fire. ----------------------------!
            if (cpatch%hite(ico) < fire_hite_threshold) then
               survivorship = fire_s_ltht(ipft)
            else
               survivorship = fire_s_gtht(ipft)
            end if
            !------------------------------------------------------------------------------!
         case (52)
            !----- Harvest precedes abandonment.  Nothing survives. -----------------------!
            survivorship = 0.0
            !------------------------------------------------------------------------------!
         end select
         !---------------------------------------------------------------------------------!

       case (6)
         !---------------------------------------------------------------------------------!
         !     Logging.  At this point a single pathway exists: cohorts above threshold    !
         ! are completely removed, and small cohorts have the same survivorship as small   !
         ! cohorts at a treefall site.  Both could be re-visited in the future, and        !
         ! different pathways for conventional and reduced-impact logging could be         !
         ! applied.                                                                        !
         !---------------------------------------------------------------------------------!
         if (cpatch%dbh(ico) >= mindbh_harvest(ipft)) then
            survivorship = 0.0
         else
            survivorship = treefall_s_ltht(ipft)
         end if
         !---------------------------------------------------------------------------------!
      end select
      !------------------------------------------------------------------------------------!

      return
   end function survivorship
   !=======================================================================================!
   !=======================================================================================!

   !=======================================================================================!
   !=======================================================================================!
   !     This function computes the hydraulic failure mortality                    .       !
   !---------------------------------------------------------------------------------------!
   real function hydro_mort_rate(plc,ipft)
       use pft_coms      , only : mort_plc_max               & ! intent(in)
                               , mort_plc_th                 & ! intent(in)
                               , is_grass                    ! ! intent(in)
     
      implicit none
      !----- Arguments. -------------------------------------------------------------------!
      real                         , intent(in) :: plc
      integer                      , intent(in) :: ipft
      !------------------------------------------------------------------------------------!
      ! local variables
      real                                      :: slope, intercept

      ! Here we assume the mort rate when plc is 100% (all conductance is lost) is 24 (die in half
      ! month)

      ! the mort rate when plc == mort_plc_th is mort_plc_max (die in one year by default)

      if (is_grass(ipft)) then
          ! is grass
          ! no need to calculate hydraulic failure moratlity since they have no stems
          hydro_mort_rate = 0.
          return
      endif

      ! TODO: need to move slope/intercept calculations into ed_params.f90

      ! we then extrapolate it linearly in log space
      slope = (log(24.) - log(mort_plc_max(ipft))) / (log(1.) - log(mort_plc_th(ipft)))
      intercept = log(mort_plc_max(ipft)) - slope * log(mort_plc_th(ipft))

      if (plc <= 0.) then
          hydro_mort_rate = 0.
      else
          hydro_mort_rate = exp(slope * log(plc) + intercept)
      endif

      return
   end function hydro_mort_rate

end module mortality
!==========================================================================================!
!==========================================================================================!
