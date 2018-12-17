program main
  use psb_base_mod
  use psb_prec_mod
  use psb_krylov_mod
  use psb_util_mod

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

  block
    ! parallel environment
    integer(psb_ipk_) :: ictxt, iam, np, info
    ! sparse matrix
    type(psb_dspmat_type) :: a
    ! descriptor
    type(psb_desc_type)   :: desc_a

    integer :: x(3), y(3)

    call psb_init(ictxt)
    call psb_info(ictxt,iam,np)

    if (iam == psb_root_) then
      write(*,*) '==> Test Partition to PSBLAS'
      print*, to_str(npart-1)
    end if

    call psb_cdall(ictxt,desc_a,info,vg=npart-1)
    call psb_cdasb(desc_a, info)

    x = [197, 369, 293]
    print*, '# process_'//to_str(iam)
    print*, to_str(x)
    call psb_glob_to_loc(x, desc_a, info)
    print*, to_str(x)

    call psb_exit(ictxt)

  end block

end program main
