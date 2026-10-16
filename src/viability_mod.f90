module viability_mod

   !=============================================================================
   ! Global Imports:
   ! These are available to every subroutine and function below 'contains'
   !=============================================================================
   use par_mod         ! Parameter definitions (numpf, numwfmem, dp)
   use com_mod         ! Global meteorological fields (ssr, tt, memtime, etc.)
   use particle_mod    ! Particle derived type (part)
   use date_mod        ! Calendar and timing variables/functions
   use qvsat_mod       ! f_qvsat, ew (Saturation vapor pressure functions)
   use windfields_mod
   use interpol_mod
   ! If you plan to calculate zval inside this module instead of passing it in:
! #ifdef ETA
!    use coord_ecmwf_mod ! zeta_to_z (Only if doing ETA transforms internally)
! #endif

   implicit none

   !=============================================================================
   ! Module-Level Variables (Optional)
   ! You can define your global viability parameters here so they are
   ! easily tunable or readable from a settings file later.
   !=============================================================================
   ! real :: uv_half_life = 14.0       ! MJ/m2 standard survival dose
   ! real :: temp_optimal_min = 15.0   ! Celsius
   ! real :: temp_optimal_max = 25.0   ! Celsius
   real :: k_sr = 0.5
   real :: k_isard = 1.0/14.0
   real :: sr_half = 12.0

contains

   subroutine calculate_decay_sr(jpart, sr_interp, sr_decay)

      ! Fits a decaying logistic function to calculate viable fraction due to solar radiation.
      ! returns the rate of decay at the current time instance
      implicit none
      integer, intent(in) :: jpart
      real, intent(in) :: sr_interp
      real, intent(out) :: sr_decay
      real :: V_sr

      V_sr = 1.0/(1.0 + exp(k_sr*(part(jpart)%sr_acc) - sr_half))
      sr_decay = k_sr*(1.0 - V_sr)*(sr_interp*1.0E-6)

   end subroutine calculate_decay_sr

   subroutine calculate_decay_sr_rh(jpart, sr_interp, rh_interp, sr_decay)

      ! Fits a decaying logistic function to calculate viable fraction due to solar radiation.
      ! returns the rate of decay at the current time instance
      implicit none
      integer, intent(in) :: jpart
      real, intent(in) :: sr_interp, rh_interp
      real, intent(out) :: sr_decay
      real :: V_sr, s
      if (rh_interp >= 0.9) then
         s = 2.0
      else
         s = 1.0
      end if

      V_sr = 1.0/(1.0 + exp(k_sr*(part(jpart)%sr_acc) - sr_half))
      sr_decay = k_sr*(1.0 - V_sr)*(s*sr_interp*1.0E-6)

   end subroutine calculate_decay_sr_rh

   subroutine calculate_decay_isard(sr_interp, sr_decay)

      ! Fits a decaying logistic function to calculate viable fraction due to solar radiation
      ! based on Isard 2005 paper. Returns the rate of decay at the current time instance.
      implicit none
      real, intent(in) :: sr_interp
      real, intent(out) :: sr_decay

      sr_decay = k_isard*(sr_interp*1.0E-6)

   end subroutine calculate_decay_isard

   subroutine calculate_decay_uv_rh_agpads(jpart, sr_interp, rh_interp, decay_rate)

      ! Fits a decaying logistic function to calculate viable fraction due to solar radiation.
      ! returns the rate of decay at the current time instance
      implicit none
      integer, intent(in) :: jpart
      real, intent(in) :: sr_interp, rh_interp
      real, intent(out) :: decay_rate
      real :: s
      if (rh_interp >= 0.9) then
         s = 2.0
      else
         s = 1.0
      end if

      decay_rate = s*sr_interp*1.0E-6

   end subroutine calculate_decay_uv_rh_agpads

   ! subroutine calculate_decay_temp(jpart, temp_interp, temp_decay)
   !    ! Fits a decaying logistic function to calculate viable fraction due to temperature.
   !    ! returns the rate of decay at the current time instance
   !
   !    implicit none
   !    integer, intent(in) :: jpart
   !    real, intent(in) :: temp_interp
   !    real, intent(out) :: temp_decay
   !    real :: V_sr
   !
   !    V_sr = 1.0/(1 + exp(k_sr*(part(jpart)%sr_acc) - sr_half))
   !    sr_decay = k_sr*(1 - V_sr)*(sr_interp/1.0E6)
   !
   ! end subroutine calculate_decay_sr

   ! subroutine calculate_decay_rh(jpart, rh_interp, rh_decay)
   !    ! Fits a decaying logistic function to calculate viable fraction due to relative humidity.
   !    ! returns the rate of decay at the current time instance
   !
   !    implicit none
   !    integer, intent(in) :: jpart
   !    real, intent(in) :: rh_interp
   !    real, intent(out) :: rh_decay
   !    real :: V_sr
   !
   !    V_sr = 1.0/(1 + exp(k_sr*(part(jpart)%sr_acc) - sr_half))
   !    sr_decay = k_sr*(1 - V_sr)*(sr_interp/1.0E6)
   !
   ! end subroutine calculate_decay_sr

   subroutine calculate_viability(itime, jpart, dt_particle, decay_rate)

      implicit none
      integer, intent(in) :: itime, jpart
      real, intent(in) :: dt_particle
      real intent(out) :: decay_rate
      real :: interp_sr, interp_rh, interp_temp

      ! This is where you will eventually call interpol_viability,
      ! compute your lambdas, and do: part(jpart)%mass = part(jpart)%mass * decay

   end subroutine calculate_viability

end module viability_mod
