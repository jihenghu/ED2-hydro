module stem_resp_driv
  contains

!==========================================================================================!
!==========================================================================================!
!     This subroutine will control the stem respiration scheme                             !
! This is called every step, but not every sub-step.                                       !
!------------------------------------------------------------------------------------------!
subroutine stem_respiration(csite,ipa)
   use ed_state_vars  , only : sitetype           & ! structure
                             , patchtype          ! ! structure
   use ed_max_dims    , only : n_pft              ! ! intent(in)
   use pft_coms       , only : stem_respiration_factor  ! ! intent(in)
   use consts_coms    , only : pi1                & ! intent(in)
                             , umols_2_kgCyr      & ! intent(in)
                             , yr_day             ! ! intent(in)
   use ed_misc_coms   , only : dtlsm              & ! intent(in)
                             , frqsum             ! ! intent(in)
   use physiology_coms, only : istem_respiration_scheme ! ! intent(in)
   implicit none
   !----- Arguments -----------------------------------------------------------------------!
   type(sitetype)            , target      :: csite             ! Current site
   integer                   , intent(in)  :: ipa               ! Current patch #
   !----- Local variables -----------------------------------------------------------------!
   type(patchtype)           , pointer     :: cpatch             ! Current site
   integer                                 :: ico                ! Current cohort #
   integer                                 :: ipft
   real                                    :: stem_area        ! m2
   !----- Locally saved variables. --------------------------------------------------------!
   real                          , save    :: dtlsm_o_frqsum
   logical                       , save    :: first_time = .true.
   !---------------------------------------------------------------------------------------!


   !----- Assign the constant scaling factor. ---------------------------------------------!
   if (first_time) then
      first_time     = .false.
      dtlsm_o_frqsum = dtlsm / frqsum
   end if
   !---------------------------------------------------------------------------------------!


   !----- Point to the cohort structures --------------------------------------------------!
   cpatch => csite%patch(ipa)
   !---------------------------------------------------------------------------------------!

   !---------------------------------------------------------------------------------------!
   !    Loop over all cohorts, from tallest to shortest.                                   !
   !---------------------------------------------------------------------------------------!
   cohortloop: do ico = 1,cpatch%ncohorts

        !----- Alias for PFT and root layer. ------------------------------------------!
        ipft  = cpatch%pft(ico)

        ! assume stem is a cylinder
        stem_area = pi1 * cpatch%dbh(ico) / 100.            & ! meter
                  * cpatch%hite(ico)                        & ! meter squared
                  * cpatch%nplant(ico)                      ! ! m2 stem/m2 ground


        select case (istem_respiration_scheme)
        case (0)
            cpatch%stem_respiration(ico) = 0.
        case (1)
            cpatch%stem_respiration(ico) = stem_resp_norm(ipft,cpatch%dbh(ico),cpatch%wood_temp(ico)) &
                                         * stem_area
                                         ! umol/m2 ground/s
        end select

        cpatch%today_stem_resp(ico) = cpatch%today_stem_resp(ico)                      &
                                    + cpatch%stem_respiration(ico)
                                    ! umol/m2 ground/s
        
        !----- The output variable must be in [kgC/plant/yr]. -------------------------!
        cpatch%fmean_stem_resp(ico)  = cpatch%fmean_stem_resp (ico)                    &
                                        + cpatch%stem_respiration(ico)                    &
                                        * dtlsm_o_frqsum * umols_2_kgCyr                  &
                                        / cpatch%nplant          (ico)

   end do cohortloop
   !---------------------------------------------------------------------------------------!


   return
end subroutine stem_respiration
!==========================================================================================!
!==========================================================================================!

!==========================================================================================!
!==========================================================================================!
!     This function determines the normalised stem respiration (umol/m2 stem surface/s)    !
!------------------------------------------------------------------------------------------!
real function stem_resp_norm(ipft,dbh,wood_temp)
   use pft_coms       , only : stem_respiration_factor  & ! intent(in)
                             , stem_resp_size_factor    & ! intent(in)
                             , rrf_low_temp             & ! intent(in)
                             , rrf_high_temp            & ! intent(in)
                             , rrf_decay_e              & ! intent(in)
                             , rrf_hor                  & ! intent(in)
                             , rrf_q10                  ! ! intent(in)
   use farq_leuning   , only : arrhenius                & ! function
                             , collatz                  ! ! function
   use rk4_coms       , only : tiny_offset              ! ! intent(in)
   use physiology_coms, only : iphysiol                 ! ! intent(in)
   use consts_coms    , only : lnexp_min8               & ! intent(in)
                             , lnexp_max8               & ! intent(in)
                             , t008                     ! ! intent(in)
   implicit none
   !----- Arguments. ----------------------------------------------------------------------!
   integer     , intent(in) :: ipft
   real(kind=4), intent(in) :: dbh
   real(kind=4), intent(in) :: wood_temp
   !----- Local variables. ----------------------------------------------------------------!
   real(kind=8)             :: wood_temp8
   real(kind=8)             :: srf08
   real(kind=8)             :: rrf_low_temp8
   real(kind=8)             :: rrf_high_temp8
   real(kind=8)             :: rrf_decay_e8
   real(kind=8)             :: rrf_hor8
   real(kind=8)             :: rrf_q108
   real(kind=8)             :: lnexplow
   real(kind=8)             :: lnexphigh
   real(kind=8)             :: tlow_fun
   real(kind=8)             :: thigh_fun
   real(kind=8)             :: srf8
   !----- External functions. -------------------------------------------------------------!
   real(kind=4)             :: sngloff
   !---------------------------------------------------------------------------------------!

   !----- Copy some variables to double precision temporaries. ----------------------------!
   wood_temp8      = dble(wood_temp                    )
   srf08           = dble(stem_respiration_factor(ipft))                                   &
                   * dble(exp(stem_resp_size_factor(ipft) * dbh))
   rrf_low_temp8   = dble(rrf_low_temp           (ipft)) + t008
   rrf_high_temp8  = dble(rrf_high_temp          (ipft)) + t008
   rrf_decay_e8    = dble(rrf_decay_e            (ipft))
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Compute the functions that will control the Rrf function for low and high temper-  !
   ! ature.  In order to avoid floating point exceptions, we check whether the temperature !
   ! will make the exponential too large or too small.                                     !
   !---------------------------------------------------------------------------------------!
   !----- Low temperature. ----------------------------------------------------------------!
   lnexplow  = rrf_decay_e8 * (rrf_low_temp8  - wood_temp8)
   lnexplow  = max(lnexp_min8,min(lnexp_max8,lnexplow))
   tlow_fun  = 1.d0 +  exp(lnexplow)
   !----- High temperature. ---------------------------------------------------------------!
   lnexphigh = rrf_decay_e8 * (wood_temp8 - rrf_high_temp8)
   lnexphigh = max(lnexp_min8,min(lnexp_max8,lnexphigh))
   thigh_fun = 1.d0 + exp(lnexphigh)
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !    Decide which functional form to use based on the physiology.  This is just to make !
   ! it look similar to the leaf respiration respiration.                                  !
   !---------------------------------------------------------------------------------------!
   select case (iphysiol)
   case (0,1,4)
      rrf_hor8 = dble(rrf_hor(ipft))
      srf8     = arrhenius(wood_temp8,srf08,rrf_hor8) / (tlow_fun * thigh_fun)
   case (2,3)
      rrf_q108 = dble(rrf_q10(ipft))
      srf8     = collatz(wood_temp8,srf08,rrf_q108)   / (tlow_fun * thigh_fun)
   end select
   !---------------------------------------------------------------------------------------!



   !---------------------------------------------------------------------------------------!
   !     Convert result to single precision.                                               !
   !---------------------------------------------------------------------------------------!
   stem_resp_norm = sngloff(srf8,tiny_offset)
   !---------------------------------------------------------------------------------------!


   return
end function stem_resp_norm
!==========================================================================================!
!==========================================================================================!





end module stem_resp_driv
