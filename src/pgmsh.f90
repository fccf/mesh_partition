program pgmsh
  use gmsh_partitioner
  use system
  implicit none

  character(:),allocatable :: file
  integer :: np, opt

  file = '../data/gmsh.msh'
  np= 4
  opt = 0

  call get_option('-np', np)
  call get_option('-opt', opt)
  call gmsh_partition(file, np, opt)

  call write_part(6)

end program pgmsh
