program mesh2graph
  use gmsh_interface
  use metis_interface
  use string
  use iso_c_binding
  implicit none

  character(:), allocatable :: file
  integer(c_int) :: options(0:39)
  type(msh_file) :: gmshf
  integer :: nd, nn, ne, enn, ie

  integer, allocatable :: eptr(:), eind(:)
  type(c_ptr) :: xadj, adjncy
  integer(c_int), pointer :: fxadj(:) => null(), fadjncy(:) => null()
  integer :: info

  file = '../data/gmsh.msh'


  call gmshf%init(file)
  nn = gmshf%get_nn()
  ne = gmshf%get_ne()
  enn = gmshf%get_enn()

  call gmshf%get_eni1(eptr, eind)

  info = METIS_SetDefaultOptions(options)
  options(17) = 1

  info = METIS_MeshToNodal(ne,nn,eptr,eind,0,xadj,adjncy)

  call c_f_pointer(xadj,fxadj,shape=[nn+1])
  call c_f_pointer(adjncy,fadjncy,shape=[fxadj(nn+1)])

  print*, 'nn = '//to_str(nn)
  print*, 'ne = '//to_str(ne)
  print*, 'xadj = '//to_str(fxadj)
  print*, 'adjncy = '//to_str(fadjncy)


end program mesh2graph
