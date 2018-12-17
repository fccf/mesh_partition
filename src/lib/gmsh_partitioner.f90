module gmsh_partitioner
  use gmsh_interface
  use metis_interface
  use string
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

    !> make face index

  end subroutine gmsh_partition
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
