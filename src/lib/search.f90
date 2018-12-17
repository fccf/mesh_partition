module search

  implicit none

contains
  !=============================================================================
  pure function binary_search(v, array) result(idx)
    integer, intent(in) :: v
    integer, intent(in) :: array(:)
    integer :: idx
    integer :: low, high, mid

    low = lbound(array, dim=1)
    high = ubound(array, dim=1)

    idx = 0

    do while ( low <= high )
      mid = (low + high)/2
      if( array(mid) == v) then
        idx = mid
        exit
      elseif(array(mid)>v ) then
        high = mid - 1
      else
        low = mid +1
      endif
    end do

  end function binary_search


end module search
