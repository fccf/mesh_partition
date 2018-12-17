module gmsh_partitioner
  use gmsh_interface
  use metis_interface
  use string
  use search
  use hdf5_interface
  use iso_c_binding
  use psb_util_mod
  implicit none


  ! private
  type part_t
    integer :: na = 0
    integer :: nn = 0
    integer :: nh = 0
    integer :: ne = 0
    integer :: nf = 0

    integer, allocatable :: all_idx(:)
    integer, allocatable :: node_idx(:)
    integer, allocatable :: halo_idx(:)
    integer, allocatable :: elem_idx(:)
    integer, allocatable :: face_idx(:)
  end type part_t

  type(part_t), allocatable :: part(:)
  type(msh_file) :: msh
  logical :: partitioned = .FALSE.

contains
  !=============================================================================
  subroutine gmsh_partition(file, np, method)
    character(*), intent(in) :: file
    integer, intent(in) :: np
    integer, intent(in), optional :: method

    integer, allocatable :: eptr(:), eind(:), fni(:,:)
    integer, allocatable :: epart(:), npart(:), bpart(:)
    integer :: ne, nn, nf, fnn, enn
    integer :: options(0:39)
    integer :: info, objval
    integer :: ie, in, ip
    integer :: method_

    method_ = 0
    if(present(method)) then
      method_ = method
    endif

    call msh%init(file)

    ne = msh%get_ne()
    nn = msh%get_nn()
    enn = msh%get_enn()
    nf  = msh%get_nf()
    fni = msh%get_fni2()
    call msh%get_eni1(eptr,eind)

    info = METIS_SetDefaultOptions(options)
    options(17) = 1

    allocate(epart(ne))
    allocate(npart(nn))
    allocate(bpart(nf))

    if(method_ == 0) then
      info = METIS_PartMeshDual(ne,nn,eptr,eind,ncommon = 0,nparts=np,options=options,&
      &    objval=objval,epart=epart,npart=npart)
    elseif(method_ == 1) then
      info = METIS_PartMeshNodal(ne,nn,eptr,eind,nparts=np,options=options,&
      &    objval=objval,epart=epart,npart=npart)
    else
      error stop "Unsupported partition method!"
    endif

    !> allocate
    allocate(part(np))
    do ip = 1, np
      part(ip)%all_idx  = [integer :: ]
      part(ip)%node_idx = [integer :: ]
      part(ip)%halo_idx = [integer :: ]
      part(ip)%elem_idx = [integer :: ]
      part(ip)%face_idx = [integer :: ]
    enddo

    !> make part(ip)%node_idx
    do in = 1, nn
      ip = npart(in)
      part(ip)%node_idx = [part(ip)%node_idx, in]
    enddo

    !> make part(ip)%elem_idx
    do ie = 1, ne
      ip = epart(ie)
      part(ip)%elem_idx = [part(ip)%elem_idx, ie]
    enddo

    !> make part(ip)%face_idx
    do ie = 1, nf
      in = fni(1,ie)
      bpart(ie) = npart(in)
    enddo
    do ie = 1, nf
      ip = bpart(ie)
      part(ip)%face_idx = [part(ip)%face_idx, ie]
    enddo

    do ip = 1, np
      part(ip)%nn = size(part(ip)%node_idx)
      part(ip)%ne = size(part(ip)%elem_idx)
      part(ip)%nf = size(part(ip)%face_idx)
    enddo

    !> make part(ip)%all_idx
    block
      integer :: ne_p, ele, ist, ied, temp, i
      integer, allocatable :: eind_p(:), eni(:,:) !< eind in process ip

      eni = msh%get_eni2()
      do ip = 1, np
        ne_p = part(ip)%ne
        allocate(eind_p(enn*ne_p))

        do ie = 1, ne_p
          ied = ie*enn
          ist = ied - enn + 1
          ele = part(ip)%elem_idx(ie)
          eind_p(ist:ied) =  eni(:,ele)  !< get eind_p
        enddo

        !> sort
        call psb_msort(eind_p)

        temp = eind_p(1)
        part(ip)%all_idx = [part(ip)%all_idx, temp]
        do i = 2, size(eind_p)
          if(eind_p(i) /= temp) then
            temp = eind_p(i)
            part(ip)%all_idx = [part(ip)%all_idx, temp]
          endif
        enddo

        part(ip)%na = size(part(ip)%all_idx)

        deallocate(eind_p)
      enddo
    end block

    !> part(ip)%halo_idx
    block
      integer :: ii, jj

      do ip = 1, np
        ii = 1
        jj = 1
        do
          if(part(ip)%all_idx(jj) == part(ip)%node_idx(ii)) then
            ii = ii + 1
            jj = jj + 1
          else
            part(ip)%halo_idx = [part(ip)%halo_idx, part(ip)%all_idx(jj)]
            jj = jj + 1
          endif

          if(jj == size(part(ip)%all_idx)) exit
        enddo

        part(ip)%nh = size(part(ip)%halo_idx)
      enddo
    end block

    partitioned = .TRUE.
    !> make face index

  end subroutine gmsh_partition
  !=============================================================================
  subroutine write_part_h5(h5f, prefix)
    type(hdf5_file),intent(in) :: h5f
    character(*), intent(in), optional :: prefix

    integer :: ip, np, nn_lc, ne_lc
    real, allocatable :: coo(:,:), xyz_lc(:)
    integer, allocatable :: eni(:,:), er(:), eind_lc(:), er_lc(:)

    character(:),allocatable :: pre_, pre
    integer :: nd, enn, i, k

    if(.not. h5f%is_open) error stop "hdf5 file has not been opened."
    if(.not. partitioned) error stop "the mesh has not been partitioned."

    pre_ = '/'
    if(present(prefix)) pre_ = trim(adjustl(prefix))

    coo = msh%get_coord()
    eni = msh%get_eni2()
    er = msh%get_er()
    nd = msh%get_nd()
    enn = msh%get_enn()
    np = size(part)

    print*, np

    do ip = 0, np-1
      nn_lc = part(ip+1)%na
      ne_lc = part(ip+1)%ne

      allocate(xyz_lc(nn_lc*nd))
      allocate(eind_lc(enn*ne_lc))
      allocate(er_lc(ne_lc))

      er_lc = er(part(ip+1)%elem_idx)
      ! xyz_lc = reshape(source=coo(:,part(ip+1)%node_idx), shape=(/nn_lc*nd/))
      xyz_lc = reshape(source=coo(:,part(ip+1)%all_idx), shape=(/nn_lc*nd/))
      eind_lc = reshape(source = eni(:,part(ip+1)%elem_idx), shape=(/ne_lc*enn/))

      print*, to_str(part(ip+1)%all_idx)
      print*, to_str(eind_lc)
      print*, size(eind_lc)
      do i=1, size(eind_lc)
        k = eind_lc(i)
        eind_lc(i) = binary_search(k, part(ip+1)%all_idx)
        print*, i, k, eind_lc(i), size(eind_lc)
      enddo

      print*, 'xcd'
      pre = pre_//'process_'//to_str(ip)

      print*, pre

      call h5f%add(pre//'/node_idx',part(ip+1)%node_idx)
      call h5f%add(pre//'/all_idx',part(ip+1)%all_idx)
      call h5f%add(pre//'/elem_idx',part(ip+1)%elem_idx)
      call h5f%add(pre//'/halo_idx',part(ip+1)%halo_idx)
      call h5f%add(pre//'/face_idx',part(ip+1)%face_idx)

      call h5f%add(pre//'/xyz',xyz_lc)
      call h5f%add(pre//'/eind',eind_lc)
      call h5f%add(pre//'/er',er_lc)

      deallocate(xyz_lc)
      deallocate(eind_lc)
      deallocate(er_lc)
    enddo

    ! print*, 'hdf5 '

  end subroutine write_part_h5
  !=============================================================================
  subroutine write_part(unit)
    integer, intent(in) :: unit

    integer :: i, nn, ne, nf

    write(unit,*) 'nn = '//to_str(msh%get_nn())
    write(unit,*) 'ne = '//to_str(msh%get_ne())
    write(unit,*) 'nf = '//to_str(msh%get_nf())

    nn = 0
    ne = 0
    nf = 0

    do i = 1, size(part)
      write(unit,*) '# part '//to_str(i)
      write(unit,*) '  nn = '//to_str(part(i)%nn)
      write(unit,*) '  ne = '//to_str(part(i)%ne)
      write(unit,*) '  na = '//to_str(part(i)%na)
      write(unit,*) '  nh = '//to_str(part(i)%nh)
      write(unit,*) '  nf = '//to_str(part(i)%nf)
      write(unit,*)
      write(unit,*) '  node_idx = '//to_str(part(i)%node_idx)
      write(unit,*) '  all_idx  = '//to_str(part(i)%all_idx)
      write(unit,*) '  halo_idx = '//to_str(part(i)%halo_idx)
      write(unit,*) '  elem_idx = '//to_str(part(i)%elem_idx)
      write(unit,*) '  face_idx = '//to_str(part(i)%face_idx)

      nn = nn + part(i)%nn
      ne = ne + part(i)%ne
      nf = nf + part(i)%nf
    enddo

    write(unit,*) 'nn = '//to_str(nn)
    write(unit,*) 'ne = '//to_str(ne)
    write(unit,*) 'nf = '//to_str(nf)

  end subroutine write_part

end module gmsh_partitioner
