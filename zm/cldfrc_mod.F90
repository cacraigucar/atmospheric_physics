module cldfrc_mod

  use ccpp_kinds, only:  kind_phys

  implicit none
  private
  save

  public cldfrc_fice

contains

  subroutine cldfrc_fice(ncol, pver, tmelt, t, fice, fsnow)
!
! Compute the fraction of the total cloud water which is in ice phase.
! The fraction depends on temperature only.
! This is the form that was used for radiation, the code came from cldefr originally
!
! Author: B. A. Boville Sept 10, 2002
!  modified: PJR 3/13/03 (added fsnow to ascribe snow production for convection )
!-----------------------------------------------------------------------

!CACNOTE - This use statement is temporary until the entire cloud_fraction module is converted to CCPP
    use cloud_fraction, only: top_lev

! Arguments
    integer,  intent(in)  :: ncol                 ! number of active columns
    integer,  intent(in)  :: pver                 ! number of vertical levels
    real(kind_phys),intent(in) :: tmelt           ! Freezing point of water (K)
    real(kind_phys), intent(in)  :: t(:,:)        ! temperature

    real(kind_phys), intent(out) :: fice(:,:)     ! Fractional ice content within cloud
    real(kind_phys), intent(out) :: fsnow(:,:)    ! Fractional snow content for convection

! Local variables
    real(kind_phys) :: tmax_fice                         ! max temperature for cloud ice formation
    real(kind_phys) :: tmin_fice                         ! min temperature for cloud ice formation
    real(kind_phys) :: tmax_fsnow                        ! max temperature for transition to convective snow
    real(kind_phys) :: tmin_fsnow                        ! min temperature for transition to convective snow

    integer :: i,k                                ! loop indexes

!-----------------------------------------------------------------------

    tmax_fice = tmelt - 10._kind_phys        ! max temperature for cloud ice formation
    tmin_fice = tmax_fice - 30._kind_phys    ! min temperature for cloud ice formation
    tmax_fsnow = tmelt                ! max temperature for transition to convective snow
    tmin_fsnow = tmelt - 5._kind_phys        ! min temperature for transition to convective snow

    fice(:,:top_lev-1) = 0._kind_phys
    fsnow(:,:top_lev-1) = 0._kind_phys

! Define fractional amount of cloud that is ice
    do k=top_lev,pver
       do i=1,ncol

! If warmer than tmax then water phase
          if (t(i,k) > tmax_fice) then
             fice(i,k) = 0.0_kind_phys

! If colder than tmin then ice phase
          else if (t(i,k) < tmin_fice) then
             fice(i,k) = 1.0_kind_phys

! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
          else
             fice(i,k) =(tmax_fice - t(i,k)) / (tmax_fice - tmin_fice)
          end if

! snow fraction partitioning

! If warmer than tmax then water phase
          if (t(i,k) > tmax_fsnow) then
             fsnow(i,k) = 0.0_kind_phys

! If colder than tmin then ice phase
          else if (t(i,k) < tmin_fsnow) then
             fsnow(i,k) = 1.0_kind_phys

! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
          else
             fsnow(i,k) =(tmax_fsnow - t(i,k)) / (tmax_fsnow - tmin_fsnow)
          end if

       end do
    end do

  end subroutine cldfrc_fice
end module cldfrc_mod
