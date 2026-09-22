#include "symbol.inc"
!> Provide routines to analyze the band gap of a VASP calculation
!>
!> The general idea is to analyze the bandgap at every k point individually such that
!> the system behaves continuous across a metal-insulator transition. If one occupied
!> the same number of states at all k points, one can define a "valence" band maximum
!> and "conduction" band minimum even for metallic systems. Note that such a definition
!> results in a negative bandgap for metals.
module bandgap_tools

    use bandgap_struct
    use prec, only: q
    use wave, only: wavespin

    implicit none

    type band_edge
        type(orbital_property) valence_band, conduction_band
    end type band_edge

    integer, parameter :: max_ = 1, min_ = 2
    integer, parameter :: vbm = max_, cbm = min_
    integer, parameter :: finished = 10, update_lower = 11, update_upper = 12
    character(len=*), parameter :: undefined = repeat(" ", 11) // "undefined"

    private
    public analyze_and_write_bandgap
    ! for unit tests
    public find_bandgap, band_indices_from_local_fermi_energy, band_indices_from_fermi_weights
    public write_to_unit_compact, write_to_unit_verbose

contains


    !> Write the bandgap information to the given unit and to the HDF5 file
    !> If a negative unit is passed no output is produced.
    subroutine analyze_and_write_bandgap(unit_, w, num_electron, format_)
        use string, only: str
        use tutor, only: vtutor
        use wave, only: wavespin
#ifdef VASP_HDF5
        use vhdf5, only: vh5_write_bandgap, ih5outfileid, grp_intermediate
#endif
        integer, intent(in) :: unit_
        type(wavespin), intent(in) :: w
        real(q), intent(in) :: num_electron
        character(len=*), intent(in) :: format_
        type(bandgap_info), allocatable :: bandgap_weight(:), bandgap_kpoint(:)
        logical is_open
        !
        if (unit_ < 0) return
        !
        PROFILING_START('bandgap')
        inquire(unit=unit_, opened=is_open)
        if (.not.is_open) &
            call vtutor%bug("The unit=" // str(unit_) // " is not associated with an open file.", __FILE__, __LINE__)
        !
        if (.not.is_sorted(real(w%celtot))) then
            call vtutor%warning("Could not determine fundamental gap because eigenvalues are not sorted.")
            PROFILING_STOP('bandgap')
            return
        end if
        !
        call find_bandgap(w, num_electron, use_fermi_weights, bandgap_weight)
        call find_bandgap(w, num_electron, use_local_fermi_energy, bandgap_kpoint)
        call write_to_unit(unit_, format_, bandgap_weight, bandgap_kpoint)
#ifdef VASP_HDF5
        call vh5_write_bandgap(ih5outfileid, grp_intermediate, use_fermi_weights, bandgap_weight)
        call vh5_write_bandgap(ih5outfileid, grp_intermediate, use_local_fermi_energy, bandgap_kpoint)
#endif
        PROFILING_STOP('bandgap')
        !
    end subroutine analyze_and_write_bandgap

    logical function is_sorted(band_structure)
        real(q), intent(in) :: band_structure(:,:,:)
        is_sorted = all(band_structure(2:,:,:) - band_structure(:size(band_structure, 1) - 1,:,:) >= 0)
    end function is_sorted


    subroutine write_to_unit(unit_, format_, bandgap_weight, bandgap_kpoint)
        integer, intent(in) :: unit_
        character(len=*), intent(in) :: format_
        type(bandgap_info), intent(in) :: bandgap_weight(:), bandgap_kpoint(:)
        if (format_ == "COMPACT") then
            call write_to_unit_compact(unit_, bandgap_weight)
        else if (format_ == "WEIGHT") then
            call write_to_unit_verbose(unit_, bandgap_weight)
        else if (format_ == "KPOINT") then
            call write_to_unit_verbose(unit_, bandgap_kpoint)
        end if
    end subroutine write_to_unit


    subroutine write_to_unit_compact(unit_, bandgap)
        integer, intent(in) :: unit_
        type(bandgap_info), intent(in) :: bandgap(:)
        call write_orbital_compact(unit_, "val. band max:", bandgap(1)%fundamental_valence)
        call write_orbital_compact(unit_, "cond. band min:", bandgap(1)%fundamental_conduction)
        call write_gap_compact(unit_, "fundamental gap:", bandgap(1)%fundamental_gap)
    end subroutine write_to_unit_compact

    subroutine write_orbital_compact(unit_, label, orbital)
        integer, intent(in) :: unit_
        character(len=*), intent(in) :: label
        type(orbital_property), intent(in) :: orbital
        character(len=*), parameter :: valid_format = "(1x,a,t19,f17.10,3x,a,3(1x,f8.4))"
        character(len=*), parameter :: invalid_format = "(1x,a,t19,6x,a,t39,a,a)"
        !
        if (orbital%band /= band_not_found) then
            write(unit_, valid_format) label, orbital%eigenvalue, "@ k =",  orbital%kvector
        else
            write(unit_, invalid_format) label, adjustl(undefined), "@ k =", undefined
        end if
    end subroutine write_orbital_compact

    subroutine write_gap_compact(unit_, label, gap)
        integer, intent(in) :: unit_
        character(len=*), intent(in) :: label
        real(q), intent(in) :: gap
        character(len=*), parameter :: valid_format = "(1x,a,t19,f17.10)"
        character(len=*), parameter :: invalid_format = "(1x,a,t19,6x,a)"
        !
        if (gap /= huge(1.0_q)) then
            write(unit_, valid_format) label, gap
        else
            write(unit_, invalid_format) label, trim(adjustl(undefined))
        end if
    end subroutine write_gap_compact


    subroutine write_to_unit_verbose(unit_, bandgap)
        integer, intent(in) :: unit_
        type(bandgap_info), intent(in) :: bandgap(:)
        character(len=*), parameter :: header = "(1x,a)"
        logical spin_polarized
        !
        spin_polarized = size(bandgap) == 3
        write(unit_, header) "Band structure"
        write(unit_, header) "--------------"
        if (spin_polarized) then
            write(unit_, header) "                       spin independent             spin component 1             spin component 2"
        else
            write(unit_, header) "                       spin independent"
        end if
        call write_band_verbose(unit_, "val. band max:", bandgap(:)%fundamental_valence)
        call write_band_verbose(unit_, "cond. band min:", bandgap(:)%fundamental_conduction)
        call write_gap_verbose(unit_, "fundamental gap:", bandgap(:)%fundamental_gap)
        call write_kpoint_verbose(unit_, "VBM @ kpoint:", bandgap(:)%fundamental_valence)
        call write_kpoint_verbose(unit_, "CBM @ kpoint:", bandgap(:)%fundamental_conduction)
        write(unit_, header)
        call write_band_verbose(unit_, "lower band:", bandgap(:)%direct_valence)
        call write_band_verbose(unit_, "upper band:", bandgap(:)%direct_conduction)
        call write_gap_verbose(unit_, "direct gap:", bandgap(:)%direct_gap)
        call write_kpoint_verbose(unit_, "@ kpoint:", bandgap(:)%direct_valence)
        write(unit_, header)
    end subroutine write_to_unit_verbose

    subroutine write_band_verbose(unit_, label, orbital)
        use string, only: str
        integer, intent(in) :: unit_
        character(len=*), intent(in) :: label
        type(orbital_property), intent(in) :: orbital(:)
        character(len=*), parameter :: band = "(1x,a,t19,3a29))"
        character(len=29) eigenvalue(size(orbital))
        integer component
        !
        do component = 1, size(orbital)
            if (orbital(component)%band /= band_not_found) then
                eigenvalue(component) = str(orbital(component)%eigenvalue, "(f20.6)")
            else
                eigenvalue(component) = undefined
            end if
        end do
        write(unit_, band) label, eigenvalue
    end subroutine write_band_verbose

    subroutine write_gap_verbose(unit_, label, gap)
        use string, only: str
        integer, intent(in) :: unit_
        character(len=*), intent(in) :: label
        real(q), intent(in) :: gap(:)
        character(len=*), parameter :: format_ = "(1x,a,t19,3a29))"
        character(len=29) gap_str(size(gap))
        integer component
        !
        do component = 1, size(gap)
            if (gap(component) /= huge(1.0_q)) then
                gap_str(component) = str(gap(component), "(f20.6)")
            else
                gap_str(component) = undefined
            end if
        end do
        write(unit_, format_) label, gap_str
    end subroutine write_gap_verbose

    subroutine write_kpoint_verbose(unit_, label, orbital)
        integer, intent(in) :: unit_
        character(len=*), intent(in) :: label
        type(orbital_property), intent(in) :: orbital(:)
        character(len=*), parameter :: format_ = "(1x,a,t19,3a29))"
        character(len=29) kvector(size(orbital))
        integer component
        !
        do component = 1, size(orbital)
            if (orbital(component)%band /= band_not_found) then
                write(kvector(component), "(1x,3(f8.4,1x))") orbital(component)%kvector
            else
                kvector(component) = undefined
            end if
        end do
        write(unit_, format_) label, kvector
    end subroutine write_kpoint_verbose


    !> Determine the bandgap of the system
    !>
    !> Creates a bandgap_info object that contains the bandgap of the system as
    !> well as the corresponding k-point coordinates and the Fermi energy for
    !> reference.
    subroutine find_bandgap(w, num_electron, algorithm, bandgap)
        type(wavespin), intent(in) :: w
        real(q), intent(in) :: num_electron
        integer, intent(in) :: algorithm
        type(bandgap_info), allocatable, intent(out) :: bandgap(:)
        integer, allocatable :: band_indices(:,:,:)
        integer, allocatable :: indices_vbm(:), indices_cbm(:), indices_direct(:)
        integer spin
        !
        call setup_bandgap(size(w%celtot, 3), bandgap)
        if (algorithm == use_fermi_weights) then
            call band_indices_from_fermi_weights(w%fertot, band_indices)
        else
            call band_indices_from_local_fermi_energy(w, num_electron, band_indices)
        end if
        call find_band_edges(w, band_indices, bandgap)
        bandgap(1)%fermi_energy = w%efermi(1)
        if (size(bandgap) == 3) &
            bandgap(2:3)%fermi_energy = w%efermi
    end subroutine find_bandgap


    subroutine setup_bandgap(num_spin, bandgap)
        integer, intent(in) :: num_spin
        type(bandgap_info), allocatable, intent(out) :: bandgap(:)
        integer component
        !
        if (num_spin == 1) then
            allocate(bandgap(1))
        else
            allocate(bandgap(3))
        end if
        !
        bandgap%fundamental_valence%eigenvalue = -huge(1.0_q)
        bandgap%fundamental_conduction%eigenvalue = huge(1.0_q)
        bandgap%fundamental_gap = huge(1.0_q)
        bandgap%direct_gap = huge(1.0_q)
    end subroutine setup_bandgap


    subroutine band_indices_from_fermi_weights(fermi_weights, band_indices)
        real(q), intent(in) :: fermi_weights(:,:,:)
        integer, allocatable, intent(out) :: band_indices(:,:,:)
        real(q), parameter :: half_filled = 0.5_q
        integer kpoint, num_kpoint, spin, num_spin
        !
        num_kpoint = size(fermi_weights, 2)
        num_spin = size(fermi_weights, 3)
        allocate(band_indices(2, num_kpoint, num_spin))
        do spin = 1, num_spin
            do kpoint = 1, num_kpoint
                band_indices(vbm, kpoint, spin) = count(fermi_weights(:, kpoint, spin) > half_filled)
                band_indices(cbm, kpoint, spin) = band_indices(vbm, kpoint, spin) + 1
            end do
        end do
    end subroutine band_indices_from_fermi_weights


    subroutine band_indices_from_local_fermi_energy(w, num_electron, band_indices)
        type(wavespin), intent(in) :: w
        real(q), intent(in) :: num_electron
        integer, allocatable, intent(out) :: band_indices(:,:,:)
        integer num_band, num_kpoint, num_spin, num_elec, spin_degen
        integer kpoint, lower, upper, middle
        !
        num_band = size(w%celtot, 1)
        num_kpoint = size(w%celtot, 2)
        num_spin = size(w%celtot, 3)
        num_elec = nint(num_electron)
        spin_degen = nint(w%wdes%rspin)
        allocate(band_indices(2, num_kpoint, num_spin))
        !
        if (num_spin == 1) then
            band_indices(vbm,:,:) = (num_elec + spin_degen - 1) / spin_degen
            band_indices(cbm,:,:) = num_elec / spin_degen + 1
            return
        end if
        !
        do kpoint = 1, num_kpoint
            lower = max(0, num_elec - num_band)
            upper = min(num_elec, num_band)
            do
                middle = (lower + upper) / 2
                if (upper <= lower) &
                    exit
                select case (bisect_step(w, num_elec, middle, kpoint))
                case (finished)
                    exit
                case (update_lower)
                    lower = middle + 1
                case (update_upper)
                    upper = middle - 1
                end select
            end do
            band_indices(vbm, kpoint, 1) = middle
            band_indices(cbm, kpoint, 1) = middle + 1
            band_indices(vbm, kpoint, 2) = num_elec - middle
            band_indices(cbm, kpoint, 2) = num_elec - middle + 1
        end do
        !
    end subroutine band_indices_from_local_fermi_energy

    integer function bisect_step(w, num_elec, band_up, kpoint) result (status_)
        use tutor, only: vtutor
        type(wavespin), intent(in) :: w
        integer, intent(in) :: num_elec, band_up, kpoint
        integer band_down, num_band
        logical lower_correct, upper_correct
        !
        num_band = size(w%celtot, 1)
        band_down = num_elec - band_up
        lower_correct = band_down == 0 .or. band_up == num_band
        if (.not.lower_correct) &
            lower_correct = real(w%celtot(band_down, kpoint, 2)) <= real(w%celtot(band_up + 1, kpoint, 1))
        upper_correct = band_up == 0 .or. band_down == num_band
        if (.not.upper_correct) &
            upper_correct = real(w%celtot(band_up, kpoint, 1)) <= real(w%celtot(band_down + 1, kpoint, 2))
        !
        if (lower_correct .and. upper_correct) then
            status_ = finished
        else if (lower_correct) then
            status_ = update_upper
        else if (upper_correct) then
            status_ = update_lower
        else
            call vtutor%bug("Could not determine new bound in interval bisection. Most " // &
                "likely this means that the eigenvalues are not sorted", __FILE__, __LINE__)
        end if
    end function bisect_step


    subroutine find_band_edges(w, band_indices, bandgap)
        type(wavespin), intent(in) :: w
        integer, intent(in) :: band_indices(:,:,:)
        type(bandgap_info), intent(inout) :: bandgap(:)
        type(band_edge), allocatable :: edge(:)
        integer kpoint, component
        real(q) not_found
        !
        do kpoint = 1, size(w%celtot, 2)
            call find_band_edges_current_kpoint(w, band_indices, kpoint, edge)
            call update_valence_band_maximum(edge, bandgap) !valence_band_maximum)
            call update_conduction_band_minimum(edge, bandgap) !conduction_band_minimum)
            call update_direct_gap(edge, bandgap) !direct_gap)
        end do
        !
        where (bandgap%fundamental_valence%band /= band_not_found .and. &
               bandgap%fundamental_conduction%band /= band_not_found)
            bandgap%fundamental_gap = bandgap%fundamental_conduction%eigenvalue - bandgap%fundamental_valence%eigenvalue
        end where
        !
    end subroutine find_band_edges


    subroutine find_band_edges_current_kpoint(w, band_indices, kpoint, edge)
        type(wavespin), intent(in) :: w
        integer, intent(in) :: band_indices(:,:,:), kpoint
        type(band_edge), allocatable, intent(out) :: edge(:)
        integer spin, num_spin, bands(2), num_band, component, num_component, optimal_band
        logical valid_bands(2)
        !
        num_band = size(w%celtot, 1)
        num_spin = size(w%celtot, 3)
        num_component = 2 * num_spin - 1
        allocate(edge(num_component))
        edge%valence_band%eigenvalue = -huge(1.0_q)
        edge%conduction_band%eigenvalue = huge(1.0_q)
        !
        do spin = 1, num_spin
            component = spin + (num_component - 1) / 2  ! equivalent to component = spin or spin + 1 for num_spin = 1 or 2
            bands = band_indices(:, kpoint, spin)
            valid_bands = 0 < bands .and. bands <= num_band
            if (valid_bands(vbm)) then
                edge(component)%valence_band = orbital_property( &
                    eigenvalue = real(w%celtot(bands(vbm), kpoint, spin)), &
                    kvector = w%wdes%vkpt(:,kpoint), &
                    band = bands(vbm), kpoint = kpoint, spin = spin &
                )
            end if
            if (valid_bands(cbm)) then
                edge(component)%conduction_band = orbital_property( &
                    eigenvalue = real(w%celtot(bands(cbm), kpoint, spin)), &
                    kvector = w%wdes%vkpt(:,kpoint), &
                    band = bands(cbm), kpoint = kpoint, spin = spin &
                )
            end if
        end do
        !
        if (num_component == 3) then
            optimal_band = maxloc(edge(:)%valence_band%eigenvalue, dim=1)
            edge(1)%valence_band = edge(optimal_band)%valence_band
            optimal_band = minloc(edge(:)%conduction_band%eigenvalue, dim=1)
            edge(1)%conduction_band = edge(optimal_band)%conduction_band
        end if
        !
    end subroutine find_band_edges_current_kpoint


    elemental subroutine update_valence_band_maximum(edge, bandgap)
        type(band_edge), intent(in) :: edge
        type(bandgap_info), intent(inout) :: bandgap
        if (edge%valence_band%eigenvalue > bandgap%fundamental_valence%eigenvalue) &
            bandgap%fundamental_valence = edge%valence_band
    end subroutine update_valence_band_maximum

    elemental subroutine update_conduction_band_minimum(edge, bandgap)
        type(band_edge), intent(in) :: edge
        type(bandgap_info), intent(inout) :: bandgap
        if (edge%conduction_band%eigenvalue < bandgap%fundamental_conduction%eigenvalue) &
            bandgap%fundamental_conduction = edge%conduction_band
    end subroutine update_conduction_band_minimum

    elemental subroutine update_direct_gap(edge, bandgap)
        type(band_edge), intent(in) :: edge
        type(bandgap_info), intent(inout) :: bandgap
        !
        if (direct_gap(edge) < bandgap%direct_gap) then
            bandgap%direct_valence = edge%valence_band
            bandgap%direct_conduction = edge%conduction_band
            bandgap%direct_gap = direct_gap(edge)
        end if
    end subroutine update_direct_gap

    elemental real(q) function direct_gap(edge)
        type(band_edge), intent(in) :: edge
        if (edge%valence_band%band == band_not_found .or. edge%conduction_band%band == band_not_found) then
            direct_gap = huge(1.0_q)
        else
            direct_gap = edge%conduction_band%eigenvalue - edge%valence_band%eigenvalue
        end if
    end function direct_gap

end module bandgap_tools

