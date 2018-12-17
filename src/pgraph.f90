program pmesh
  use gmsh_interface
  use metis_interface
  use string
  use iso_c_binding
  implicit none

  character(:), allocatable :: file
  integer(c_int) :: options(0:39)
  type(msh_file) :: gmshf
  integer, allocatable :: eni(:,:)
  integer :: nd, nn, ne, enn
  integer :: ie

  integer :: info

  integer, allocatable :: eptr(:), eind(:), epart(:), npart(:)
  type(c_ptr) :: xadj, adjncy
  integer(c_int), pointer :: fxadj(:) => null(), fadjncy(:) => null()
  integer :: ist,ied, objval

  file = '../data/gmsh.msh'
  ! type(vtk_file) :: vtk

  call gmshf%init(file)
  nn = gmshf%get_nn()
  ne = gmshf%get_ne()
  enn = gmshf%get_enn()

  call gmshf%get_eni1(eptr, eind)
  allocate(epart(ne), npart(nn))



  info = METIS_SetDefaultOptions(options)
  options(17) = 1

  info = METIS_MeshToNodal(ne,nn,eptr,eind,0,xadj,adjncy)

  call c_f_pointer(xadj,fxadj,shape=[nn+1])
  call c_f_pointer(adjncy,fadjncy,shape=[fxadj(nn+1)])

  info = METIS_PartGraphRecursive(nn,ncon=1,xadj=fxadj,adjncy=fadjncy,nparts=4,objval=objval,part=npart)

  print*, 'nn = '//to_str(nn)
  print*, 'ne = '//to_str(ne)
  print*, 'npart = '//to_str(npart)


end program pmesh
