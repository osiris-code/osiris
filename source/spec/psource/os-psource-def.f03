
! This file is meant to be included through an #include statement

!-----------------------------------------------------------------------------------------
! Particle source abstract class definitions
!-----------------------------------------------------------------------------------------
type, abstract :: t_psource

contains
    procedure(if_inject_int), deferred  :: if_inject
    procedure(read_input_int), deferred :: read_input
    procedure(cleanup_int), deferred    :: cleanup
    procedure(inject_int), deferred     :: inject
end type

public :: t_psource


abstract interface

function if_inject_int(this)
    import t_psource
    logical :: if_inject_int
    class(t_psource), intent(in) :: this
end function if_inject_int

end interface


abstract interface

subroutine read_input_int(this, input_file, coordinates )
    import t_psource, t_input_file
    class(t_psource), intent(inout) :: this
    class(t_input_file), intent(inout) :: input_file
    integer, intent(in) :: coordinates
end subroutine read_input_int

end interface

abstract interface

subroutine cleanup_int( this )
    import t_psource
    class(t_psource), intent(inout) :: this
end subroutine cleanup_int

end interface


abstract interface

function inject_int( this, species, ig_xbnd_inj, jay, no_co, bnd_cross, node_cross, &
                             send_msg, recv_msg )
    import t_psource, t_species, t_current, t_node_conf, t_part_idx, t_spec_msg
    class(t_psource), intent(inout) :: this
    class(t_species), intent(inout), target :: species
    integer, dimension(:, :), intent(in) :: ig_xbnd_inj
    class( t_current ), intent(inout)   :: jay
    class( t_node_conf ), intent(in)     :: no_co
      type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
      type( t_part_idx ), intent(inout) :: node_cross
      type( t_spec_msg ), dimension(2), intent(inout) :: send_msg, recv_msg

    integer :: inject_int
end function inject_int

end interface

