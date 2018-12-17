program pmesh1
  use gmsh_interface
  use metis_interface
  use vtk_interface
  use string
  use iso_c_binding
  use system
  implicit none

  character(:), allocatable :: m_file, v_file
  type(msh_file) :: msh
  type(vtk_file) :: vtk

  integer, allocatable :: eptr(:), eind(:), epart(:), npart(:)
  integer, allocatable :: eni(:,:)
  real, allocatable    :: coo(:,:)
  integer(c_int) :: options(0:39)
  integer :: nd, nn, ne, enn
  integer :: nparts,objval
  integer :: info

  m_file = '../data/gmsh.msh'

  nparts = 4
  if(have_option('-np')) call get_option('-np',nparts)

  v_file = 'METIS_PartMeshDual_'//to_str(nparts)//'.vtk'

  call msh%init(m_file)
  nn = msh%get_nn()
  ne = msh%get_ne()
  enn = msh%get_enn()
  eni = msh%get_eni2()
  coo = msh%get_coord()
  call msh%get_eni1(eptr, eind)

  allocate(epart(ne), npart(nn))

  info = METIS_SetDefaultOptions(options)
  options(17) = 1

  info = METIS_PartMeshDual(ne,nn,eptr,eind,ncommon = 1,nparts=nparts,options=options,&
      objval=objval,epart=epart,npart=npart)

  call vtk%write(v_file,'mesh_part',coo,eni,real(epart))

end program pmesh1
