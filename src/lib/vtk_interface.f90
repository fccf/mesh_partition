module vtk_interface
  use txt_interface, only: txt_file
  implicit none

  public :: vtk_file

  private

  type, extends(txt_file) ::  vtk_file
  contains
    procedure :: write => vtk_write
    procedure :: write_header => vtk_write_header
    procedure :: write_coord  => vtk_write_coord
    procedure :: write_cell   => vtk_write_cell
    procedure :: write_cell_value => vtk_write_cell_value
    procedure :: write_node_value => vtk_write_node_value
  end type vtk_file

contains
  !=============================================================================
  subroutine vtk_write(this,file,title,coord,en,data)
    class(vtk_file), intent(inout) :: this
    character(*), intent(in) :: file,title
    real, intent(in) :: coord(:,:),data(:)
    integer, intent(in) :: en(:,:)

    integer :: nn,ne,nnv
    nn = size(coord,2)
    ne = size(en,2)
    nnv = size(data)

    call this%open(file)
    call this%write_header(title)

    call this%write_coord(coord)

    call this%write_cell(en)

    if(nnv==nn) then
      call this%write_node_value(data)
    elseif(nnv==ne) then
      call this%write_cell_value(data)
    else
      error stop 'error : dimension error'
    endif

    call this%close()

  end subroutine vtk_write
  !=============================================================================
  subroutine vtk_write_header(this,title)
    class(vtk_file), intent(in) :: this
    character(*), intent(in) :: title

    write (this%unit, '(a)' ) '# vtk DataFile Version 3.0'
    write (this%unit, '(a)' ) title
    write (this%unit, '(a)' ) 'ASCII'
    write (this%unit, '(a)' ) 'DATASET UNSTRUCTURED_GRID'

  end subroutine vtk_write_header
  !=============================================================================
  subroutine vtk_write_coord(this,coord)
    class(vtk_file), intent(in) :: this
    real, intent(in) :: coord(:,:)

    integer :: nn,nd,i

    nn = size(coord,2)
    nd = size(coord,1)

    write (this%unit, '(a)' ) 'POINTS '//str(nn)//' float'

    if(nd == 3) then
      do i = 1,nn
        write(this%unit,*) coord(1:3,i)
      enddo
    elseif(nd == 2) then
      do i = 1,nn
        write(this%unit,*) coord(1:2,i),0.0d0
      enddo
    elseif(nd== 1) then
      do i = 1,nn
        write(this%unit,*) coord(1:1,i),0.0d0,0.0d0
      enddo
    else

      error stop 'error : dimension error'
    endif

  end subroutine vtk_write_coord
  !=============================================================================
  ! VTK cell type https://www.vtk.org/doc/nightly/html/vtkCellType_8h.html
  subroutine vtk_write_cell(this,en)
    class(vtk_file), intent(in) :: this
    integer, intent(in) :: en(:,:)

    integer :: ne,enn,i

    ne = size(en,2)
    enn= size(en,1)

    !>Note that the element node indices must be converted from 1-based to 0-based
    !  before being written out.

    write (this%unit, '(a)' ) ' '
    write (this%unit, '(a)' ) 'CELLS '//str(ne)//' '//str((enn+1)*ne)
    do i = 1,ne
      write (this%unit, *) enn, en(1:enn,i) - 1
    end do

    write (this%unit, '(a)' ) ' '
    write (this%unit, '(a)' ) 'CELL_TYPES '//str(ne)

    if(enn == 2) then
      do i = 1,ne
        write (this%unit, *) 3
      end do
    elseif(enn==3) then
      do i = 1,ne
        write (this%unit, *) 5
      end do
    elseif(enn==4) then
      do i = 1,ne
        write (this%unit, *) 10
      end do
    else
      error stop 'error : Unsupported cell types'
    endif

  end subroutine vtk_write_cell
  !=============================================================================
  subroutine vtk_write_cell_value(this,ev)
    class(vtk_file), intent(in) :: this
    real, intent(in) :: ev(:)

    integer :: ne, i

    ne = size(ev)

    write (this%unit, '(a)' ) ' '
    write (this%unit, '(a)' ) 'CELL_DATA '//str(ne)
    write (this%unit, '(a)' ) 'SCALARS flux float'
    write (this%unit, '(a)' ) 'LOOKUP_TABLE default'
    do i = 1, ne
      write (this%unit,*) ev(i)
    end do

  end subroutine vtk_write_cell_value
  !=============================================================================
  subroutine vtk_write_node_value(this,nv)
    class(vtk_file), intent(in) :: this
    real, intent(in) :: nv(:)

    integer :: nn, i

    nn = size(nv)

    write (this%unit, '(a)' ) ' '
    write (this%unit, '(a)' ) 'POINT_DATA '//str(nn)
    write (this%unit, '(a)' ) 'SCALARS flux float'
    write (this%unit, '(a)' ) 'LOOKUP_TABLE default'
    do i = 1, nn
      write (this%unit,*) nv(i)
    end do

  end subroutine vtk_write_node_value

  !=============================================================================
  !> convert integer to string
  function str(i4) result(s)
    integer, intent(in)          :: i4
    character(len=:),allocatable :: s
    character(len=20) :: str_temp

    write(str_temp,*) i4
    s = trim(adjustl(str_temp))

  end function str

end module vtk_interface
