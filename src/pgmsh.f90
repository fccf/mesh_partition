program pgmsh
  use gmsh_partitioner
  use system
  use hdf5_interface
  implicit none

  character(:),allocatable :: file, h5fname
  type(hdf5_file) :: h5f
  integer :: np, opt

  file = '../data/gmsh.msh'
  h5fname = 'fish.h5'

  np= 4
  opt = 0

  call get_option('-np', np)
  call get_option('-opt', opt)
  call gmsh_partition(file, np, opt)

  call write_part(6)

  call h5f%initialize()
  call h5f%open_file(h5fname,status='new',action='rw')
  call write_part_h5(h5f)
  call h5f%close_file()

  call h5f%finalize()

end program pgmsh
