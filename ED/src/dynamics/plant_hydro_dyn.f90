!==========================================================================================!
!==========================================================================================!
!     This module contains the wrappers for the Runge-Kutta integration scheme.            !
!==========================================================================================!
!==========================================================================================!
module plant_hydro_dyn

   contains
   !=======================================================================================!
   !=======================================================================================!
   !      Driver of plant hydrodynamics. This subroutine tracks changes of leaf/stem water !
   ! water potential as well as water flows from root to stem and from stem to canopy.     !
   ! This subroutine currently works at DTLSM level, should be integrated into
   ! RK4 scheme in the future.   [XXT]
   !---------------------------------------------------------------------------------------!
   subroutine update_plant_hydrodynamics(csite,initp,ipa,hdid) 
      use ed_state_vars   , only : sitetype             & ! structure
                                 , patchtype            ! ! structure
      use ed_misc_coms    , only : dtlsm                & ! intent(in)
                                 , frqsum               & ! intent(in)        
                                 , current_time         ! ! intent(in) 
      use soil_coms       , only : dslzi8               &
                                 , dslz8                &
                                 , nzg                  &
                                 , soil8                &
                                 , slz                  &
                                 , dslz                 &
                                 , soil                    ! ! intent(in)
      use consts_coms     , only : wdns8                &
                                 , wdnsi8               &
                                 , wdns                 &
                                 , pi1                  
      use rk4_coms        , only : rk4patchtype         &
                                 , rk4site              &
                                 , tiny_offset
      use therm_lib8     , only  : tl2uint8             & ! function
                                 , uextcm2tl8           ! ! function
      use pft_coms       , only  : Ks_stem              &
                                 , Ks_stem_b            &
                                 , psi50                &
                                 , Cap_leaf             &
                                 , Cap_stem             &
                                 , SRA                  &
                                 , rho                  &
                                 , root_beta            &
                                 , xylem_fraction       &
                                 , vessel_curl_factor   &
                                 , is_grass             
      implicit none
      !----- Arguments --------------------------------------------------------------------!
      type(sitetype)        , target      :: csite
      type(rk4patchtype)    , target      :: initp
      integer               , intent(in)  :: ipa
      real(kind=8)          , intent(in)  :: hdid                   !past timestep

      !----- Local Vars  --------------------------------------------------------------------!
      type(patchtype)       , pointer       :: cpatch
      !----- temporary state variables for water flow within plant
      real                                  :: transp               !transpiration
      real                                  :: J_sr                 !soil-root flow
      real                                  :: J_rl                 !root-leaf flow
      real                                  :: org_psi_leaf         !original leaf water potential
      real                                  :: org_psi_stem         !original stem water potential
      !----- variables necessary to calcualte water flow
      real                                  :: ap
      real                                  :: bp
      real                                  :: stem_cond            !stem conductance
      real                                  :: c_leaf               !leaf water capacitance
      real                                  :: c_stem               !stem water capacitance
      real                                  :: RAI                  !root area index
      real                                  :: cohort_crown_area
      real                                  :: wgpfrac              !relative soil water
      real                                  :: soil_water_cond
      real   , dimension(nzg)               :: layer_psi
      real                                  :: weighted_soil_psi
      real                                  :: weighted_soil_cond
      real                                  :: above_layer_depth
      real                                  :: current_layer_depth
      real                                  :: total_water_supply
      real(kind=8)      ,dimension(nzg)     :: layer_water_supply
      !----- variables to deal with energy budget of plant water absorption
      real(kind=8)                          :: dsw                  !change of soil water
      real(kind=8)                          :: org_soil_tempk 
      !----- variables for loops
      integer                               :: ico
      integer                               :: k
      integer                               :: nsoil
      integer                               :: ipft
      !--------------- other
      logical,parameter                     :: quality_check = .false.
      logical                               :: update_psi
      !----- External function ------------------------------------------------------------!
      real        , external          :: sngloff
      !----- Locally saved variables. -----------------------------------------------------!
      real        , save              :: dtlsm_o_frqsum
      logical     , save              :: first_time = .true.
      !------------------------------------------------------------------------------------!

      !----- Assign the constant scaling factor. ------------------------------------------!
      if (first_time) then
         first_time     = .false.
         dtlsm_o_frqsum = dtlsm / frqsum
      end if
      !------------------------------------------------------------------------------------!
      !------------ Check whether it is time to updated psi-----------------------------------------
      ! We update the dmin/dmax psi, one timestep after new_day.
      ! Because we need to conserve the value for phenology_driv and output
      !------------------------------------------------------------
      update_psi       = current_time%time == dtlsm

      cpatch => csite%patch(ipa)

      !----------------------------------------------------------------------
      ! Preparation for plant hydrodynamics
      !-----------------------------------------------------------------------
      
      ! Calculate layer_psi, soil water potential of each soil layer
      do k = 1,nzg
        nsoil = rk4site%ntext_soil(k)
        
        !get relative soil moisture
        wgpfrac = min(1.0,real(initp%soil_water(k) * initp%soil_fracliq(k) / soil8(nsoil)%slmsts))

        !To ensure numerical stability, we assign a very large negative number
        !to soil layer water potential if soil moisture is very low.

        if (wgpfrac < 1e-3) then
            layer_psi(k) = -1e6
        else
            layer_psi(k) = soil(nsoil)%slpots / wgpfrac ** soil(nsoil)%slbs
        endif

      enddo
      

      !----------------------------------------------------------------------
      ! Update plant hydrodynamics
      !-----------------------------------------------------------------------

      !Loop over all the cohorts
      cohortloop:do ico = 1,cpatch%ncohorts
        ipft = cpatch%pft(ico)
      
        !1. get transpiration rate
        if (initp%leaf_resolvable(ico)) then
            transp = sngloff(initp%fs_open(ico) * initp%psi_open(ico)              &
                    + (1.0d0 - initp%fs_open(ico)) *                               &
                    initp%psi_closed(ico),tiny_offset) * cpatch%lai(ico) / cpatch%nplant(ico)   &
                    / sngl(hdid)
        else
            transp = 0.
        endif

        ! the unit of transp is kg/sec
          

        !2. update leaf, stem psi and water flows

        !2.1 update leaf psi while assuming stem psi is constant....
        !    only calcualte the leaf psi if it is not grass. For grass, we treat
        !    leaf_psi the same as stem_psi
        

        if (.not. is_grass(ipft)) then

            ! Calculate the change of psi_leaf, assuming stem_psi as constant
            c_leaf = Cap_leaf(ipft) * cpatch%lai(ico) / cpatch%nplant(ico)

            stem_cond =  Ks_stem(ipft)   &                                                    
                         *  1 / (1 + (cpatch%psi_stem(ico) / psi50(ipft)) ** Ks_stem_b(ipft))  &  !cavitation effect
                         * xylem_fraction(ipft) * (pi1 * (cpatch%dbh(ico) / 200.) ** 2)         &  !conducting area     m2
                         / (cpatch%hite(ico) * vessel_curl_factor(ipft))                          !conducting length   m
            
            if (c_leaf > 0. .and. stem_cond > 0.) then
                ! if there are leaves and the tree has at least some level of
                ! conductivity
                ap = - stem_cond / c_leaf
                ! the unit of ap is s-1

                bp = ((cpatch%psi_stem(ico) - cpatch%hite(ico)) * stem_cond -   &                       
                         transp)                                                &
                         / c_leaf                                                           
                ! the unit of bp is m s-1

                ! calculate new psi_leaf
                org_psi_leaf = cpatch%psi_leaf(ico)
                cpatch%psi_leaf(ico) = ((ap * org_psi_leaf + bp) *  &!  m s-1
                                         exp(ap * dtlsm) - bp) / ap

                ! calculate the average sapflow rate from stem to leaf within the
                ! time step
                J_rl = (cpatch%psi_leaf(ico) - org_psi_leaf) * c_leaf / dtlsm + &
                        transp ! kgH2O s-1
            elseif (c_leaf > 0.) then
                ! If there is no conductance in the trunk but the cohort still
                ! has some leaves...

                ! changes of psi_leaf is solely due to transpiration
                    
                org_psi_leaf = cpatch%psi_leaf(ico)
                cpatch%psi_leaf(ico) = org_psi_leaf -                           &
                             transp * dtlsm  / c_leaf
                J_rl = 0.

            else
                ! If there is no leaves...
                ! psi_leaf is set to be equal to psi_stem minus the gravitational
                ! effect
                org_psi_leaf = cpatch%psi_leaf(ico)
                ! only include gravitational effect
                cpatch%psi_leaf(ico) = cpatch%psi_stem(ico) - cpatch%hite(ico) 
                J_rl = 0.

            endif

        else   ! If it is grass....

            org_psi_leaf = cpatch%psi_leaf(ico)
            J_rl = transp
            ! For the convenience of calculation later, we attribute J_rl as
            ! transp here...

            ! psi_leaf will be updated later

        endif

        ! For debugging purpose
        if(isnan(J_rl)) then
            write (unit=*,fmt='(80a)')         ('=',k=1,80)
            write (unit=*,fmt='(a)')           'Sapflow is NaN!!'
            write (unit=*,fmt='(a,1x,i9)')   ' + ico:                 ',ico
            write (unit=*,fmt='(a,1x,es12.4)')   ' + LEAF CAPACITANCE:    ',c_leaf
            write (unit=*,fmt='(a,1x,es12.4)')   ' + LAI:                 ',cpatch%lai(ico)
            write (unit=*,fmt='(a,1x,es12.4)')   ' + NPLANT:              ',cpatch%nplant(ico)
            write (unit=*,fmt='(a,1x,es12.4)')   ' + TRANSPIRATION:       ',transp
            write (unit=*,fmt='(a,1x,es12.4)')   ' + STEM CONDUCTANCE:    ',stem_cond
            write (unit=*,fmt='(a,1x,es12.4)')   ' + ORG LEAF PSI:        ',org_psi_leaf
            write (unit=*,fmt='(a,1x,es12.4)')   ' + ORG STEM PSI:        ',cpatch%psi_stem(ico)
           
            call fatal_error('Sapflow is wrong','update_plant_hydrodynamics'               &
                            &,'plant_hydrodynamics.f90')
        endif

        
        !2.2 update stem psi while assuming soil psi is constant....

        weighted_soil_psi  = 0.
        weighted_soil_cond = 0.
        layer_water_supply = 0.
        total_water_supply = 0.

        ! loop over all soil layers to get the aggregated water conductance
        do k = cpatch%krdepth(ico),nzg
            current_layer_depth = -slz(k)
            if (k+1 .le. nzg) then
                above_layer_depth = -slz(k+1)
            else
                above_layer_depth = 0.0
            endif
                                       
            !  Calculate RAI in each layer
            !  Based on Katul et al. 2003 PCE
            cohort_crown_area = cpatch%crown_area(ico) / cpatch%nplant(ico)

            if (cohort_crown_area == 0.) then
                RAI = 0.
            else
                RAI = cpatch%broot(ico) * & !kgC
                    (root_beta(ipft) ** (above_layer_depth / -slz(cpatch%krdepth(ico))) - &
                    root_beta(ipft) ** (current_layer_depth / -slz(cpatch%krdepth(ico))) ) * &
                    SRA(ipft) / &  !    m2/kgC
                    (2.0 * cohort_crown_area)     ! m2
            endif

            ! The unit of RAI is m2/m2

            !  Calculate soil water conductance
            nsoil = rk4site%ntext_soil(k)
            wgpfrac = min(1.0,initp%soil_water(k) * initp%soil_fracliq(k) / soil(nsoil)%slmsts)
            
            soil_water_cond = soil(nsoil)%slcons * wgpfrac ** (2.0 * soil(nsoil)%slbs + 3.0) * 1.e3 &   ! kgH2O m-2 s-1
                                  * sqrt(RAI) / (pi1 * dslz(k))  &  ! m-1
                                  * 2.0 * cohort_crown_area   ! m2


            ! disable hydraulic redistribution
            if (layer_psi(k) <= cpatch%psi_stem(ico)) then
                soil_water_cond = 0.0
            endif
            ! The unit of soil_water_cond is kgH2O m-1 s-1

            
            ! Calculate weighted conductance, weighted psi, and
            ! water_supply_layer_frac
            weighted_soil_cond = weighted_soil_cond + soil_water_cond
            weighted_soil_psi  = weighted_soil_psi + soil_water_cond * layer_psi(k) ! kgH2O s-1


            layer_water_supply(k) = &
                    dble(soil_water_cond * (layer_psi(k) - cpatch%psi_stem(ico)))!kgH2O s-1 
        enddo

        ! Update psi_stem
        org_psi_stem = cpatch%psi_stem(ico)
        
        if (is_grass(ipft)) then
            ! for grass, lump leaf and stem water capacity together
             c_stem = Cap_stem(ipft) * & ! kg H2O m-3 m-1
                    (cpatch%broot(ico) + cpatch%bdead(ico) * xylem_fraction(ipft)) / (rho(ipft) / 2.) + &
                    Cap_leaf(ipft) * cpatch%lai(ico) / cpatch%nplant(ico)
        else
            ! for trees, only use stem capacitance
             c_stem = Cap_stem(ipft) * & ! kg H2O m-3 m-1
                    (cpatch%broot(ico) + cpatch%bdead(ico) * xylem_fraction(ipft)) / (rho(ipft) / 2.)
        endif


        if (c_stem > 0. .and. weighted_soil_cond > 0.) then
            ! If the tree stem has some kind of capacity

            ap = - weighted_soil_cond  & !kgH2O m-1 s-1
                / c_stem  ! kgH2O m-1
            !the unit of ap is s-1

            bp = (weighted_soil_psi - J_rl) & ! kgH2O s-1
                / c_stem                    ! kgH2O m-1
            !the unit of bp is m s-1
            cpatch%psi_stem(ico) = ((ap * org_psi_stem + bp) * exp(ap * dtlsm)- bp) &
                                / ap
            J_sr = (cpatch%psi_stem(ico) - org_psi_stem) * c_stem  / dtlsm + &! kgH2O s-1
                    J_rl
        elseif (c_stem > 0.) then
                    
            ! the plants cannot take up water
            ! change of psi_stem is solely due to J_rl
            cpatch%psi_stem(ico) = org_psi_stem - J_rl * dtlsm /c_stem
            J_sr = 0.0
                
        else
            ! This is when the tree stem has no capacity
            ! This means the all the sapwood and fine roots are dead
            ! psi_stem cannot change in this case
            cpatch%psi_stem(ico) = org_psi_stem
            J_sr = 0.0

        endif

        ! For debugging purpose
        if(isnan(J_sr)) then
            write (unit=*,fmt='(80a)')         ('=',k=1,80)
            write (unit=*,fmt='(a)')           'Soil to root water flow is NaN!!'
            write (unit=*,fmt='(a,1x,i9)')   ' + ico:                 ',ico
            write (unit=*,fmt='(a,1x,es12.4)')   ' + STEM CAPACITANCE:    ',c_stem
            write (unit=*,fmt='(a,1x,es12.4)')   ' + RAI:                 ',RAI
            write (unit=*,fmt='(a,1x,es12.4)')   ' + NPLANT:              ',cpatch%nplant(ico)
            write (unit=*,fmt='(a,1x,es12.4)')   ' + SAPFLOW:             ',J_rl
            write (unit=*,fmt='(a,1x,es12.4)')   ' + SOIL CONDUCTANCE:    ',weighted_soil_cond
            write (unit=*,fmt='(a,1x,es12.4)')   ' + ORG LEAF PSI:        ',org_psi_leaf
            write (unit=*,fmt='(a,1x,es12.4)')   ' + ORG STEM PSI:        ',cpatch%psi_stem(ico)
           
            call fatal_error('soil to root water flow is wrong','update_plant_hydrodynamics'               &
                            &,'plant_hydrodynamics.f90')
        endif


        !now record the water potential and water flow values

        if (is_grass(ipft)) then
            ! for grass, we assume leaf and stem are the same
            cpatch%psi_leaf(ico) = cpatch%psi_stem(ico)
            J_rl = 0.
        endif

        cpatch%water_flux_rl(ico) = J_rl
        cpatch%water_flux_sr(ico) = J_sr

        !----------------------------------------------------------------------
        ! integrate fast-analysis and dmax/dmin variables
        !----------------------------------------------------------------------
        cpatch%fmean_psi_leaf(ico)      = cpatch%fmean_psi_leaf(ico)        &
                                        + cpatch%psi_leaf(ico) * dtlsm_o_frqsum
        cpatch%fmean_psi_stem(ico)      = cpatch%fmean_psi_stem(ico)        &
                                        + cpatch%psi_stem(ico) * dtlsm_o_frqsum
        cpatch%fmean_water_flux_rl(ico) = cpatch%fmean_water_flux_rl(ico)   &
                                        + cpatch%water_flux_rl(ico) * dtlsm_o_frqsum
        cpatch%fmean_water_flux_sr(ico) = cpatch%fmean_water_flux_sr(ico)   &
                                        + cpatch%water_flux_sr(ico) * dtlsm_o_frqsum

        if (cpatch%dmax_psi_leaf(ico) == 0. .or. update_psi) then
            cpatch%dmax_psi_leaf(ico) = cpatch%psi_leaf(ico)
        else
            cpatch%dmax_psi_leaf(ico) = max(cpatch%dmax_psi_leaf(ico),      &
                                            cpatch%psi_leaf(ico))
        endif

        if (cpatch%dmin_psi_leaf(ico) == 0. .or. update_psi) then
            cpatch%dmin_psi_leaf(ico) = cpatch%psi_leaf(ico)
        else
            cpatch%dmin_psi_leaf(ico) = min(cpatch%dmin_psi_leaf(ico),      &
                                            cpatch%psi_leaf(ico))
        endif

        if (cpatch%dmax_psi_stem(ico) == 0. .or. update_psi) then
            cpatch%dmax_psi_stem(ico) = cpatch%psi_stem(ico)
        else
            cpatch%dmax_psi_stem(ico) = max(cpatch%dmax_psi_stem(ico),      &
                                            cpatch%psi_stem(ico))
        endif

        if (cpatch%dmin_psi_stem(ico) == 0. .or. update_psi) then
            cpatch%dmin_psi_stem(ico) = cpatch%psi_stem(ico)
        else
            cpatch%dmin_psi_stem(ico) = min(cpatch%dmin_psi_stem(ico),      &
                                            cpatch%psi_stem(ico))
        endif

        !----------------------------------------------------------------------
        ! Extract water from soil layers
        !-----------------------------------------------------------------------
  
        ! update layer_water supply according to the final J_sr value
        ! the unit should be kgH2O m-2
  
        total_water_supply = J_sr * cpatch%nplant(ico)  !kgH2O m-2 s-1
  
        if (sum(layer_water_supply) /= 0.) then
          layer_water_supply = layer_water_supply / sum(layer_water_supply) * &
              dble(total_water_supply * dtlsm)
        else
          layer_water_supply = dble(0.)
        endif
  
  
        do k = nzg, 1,-1
           nsoil = rk4site%ntext_soil(k)
  
           ! changes of volumetric soil water content
           dsw = max(min(initp%soil_water(k) - soil8(nsoil)%soilcp, &
                      layer_water_supply(k) * wdnsi8 * dslzi8(k)),&
                      initp%soil_water(k) - soil8(nsoil)%slmsts)
           
           ! if this layer is not enough, extracts from the next layer
           if (k - 1 > 0) then
              layer_water_supply(k-1) = &
                  layer_water_supply(k-1) +           &
                  layer_water_supply(k) - dsw * &
                  dslz8(k) * wdns8
           endif
  
           layer_water_supply(k) = dsw * wdns8 * dslz8(k)
  
           initp%soil_water(k) = initp%soil_water(k) - dsw
           initp%soil_energy(k) = initp%soil_energy(k) - &
                                   dsw * tl2uint8(initp%soil_tempk(k),1.d0)
    

        enddo
  
  
        if (quality_check) then
           write (unit=*,fmt='(80a)')         ('=',k=1,80)
           write (unit=*,fmt='(a)')           'Plant Hydrodynamics Quality Check:'
           write (unit=*,fmt='(a,1x,i9)')   ' + HOUR:                ',current_time%hour
           write (unit=*,fmt='(a,1x,i9)')   ' + ICO:                 ',ico
           write (unit=*,fmt='(a,1x,i9)')   ' + PFT:                 ',cpatch%pft(ico)
           write (unit=*,fmt='(a,1x,f9.4)')   ' + DBH:                 ',cpatch%dbh(ico)
           write (unit=*,fmt='(a,1x,i9)')   ' + KRDEPTH:             ',cpatch%krdepth(ico)
           write (unit=*,fmt='(a,1x,es12.4)')   ' + PSI_STEM:            ',cpatch%psi_stem(ico)
           write (unit=*,fmt='(a,1x,es12.4)')   ' + PSI_LEAF:            ',cpatch%psi_leaf(ico)
           write (unit=*,fmt='(a,1x,es12.4)')   ' + TRANSP:              ',transp
           write (unit=*,fmt='(a,1x,es12.4)')   ' + WATER_FLUX_RL:       ',cpatch%water_flux_rl(ico)
           write (unit=*,fmt='(a,1x,es12.4)')   ' + WATER_FLUX_SR:       ',cpatch%water_flux_sr(ico)
           write (unit=*,fmt='(a,1x,3f9.4)')  ' + SOIL_WATER (top 3 layer):       ',real(initp%soil_water(nzg-2:nzg),4)
           write (unit=*,fmt='(a,1x,3f9.4)')  ' + SOIL_TEMPK (top 3 layer):       ',real(initp%soil_tempk(nzg-2:nzg),4)
        endif

    enddo cohortloop
                
    !Finally, update soil temp and fracliq in initp although soil
    !temperature and fracliq should not change in this case
    do k = nzg, 1,-1
       nsoil = rk4site%ntext_soil(k)
       call uextcm2tl8(initp%soil_energy(k),initp%soil_water(k) * wdns8,&
                       soil8(nsoil)%slcpd,initp%soil_tempk(k),initp%soil_fracliq(k))
    enddo


    
   end subroutine update_plant_hydrodynamics

end module plant_hydro_dyn

!==========================================================================================!
!==========================================================================================!
