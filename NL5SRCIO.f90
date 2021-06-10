!/ ------------------------------------------------------------------- /
! Compile this module for python by:
! $ f2py -c --fcompiler=gnu95 --f90flags="-fconvert=big-endian" -m NL5SRCIO NL5SRCIO.f90
! $ f2py -c --fcompiler=gnu95 -m NL5SRCIO NL5SRCIO.f90
!
! →  Somehow, the f90flags option does not work on f2py ! But, anyway,
!    `convert='big_endian'` option in the `open` command works for both
!    gnu and intel fortran compilers.
!/ ------------------------------------------------------------------- /
! Purpose: Fortran subroutines for reading the binary file containing GKE
!          nonlinear source term Snl and Kurtosis
!
! Org: Qingxiang Liu (31-Mar-2019) [binary only]
! Mod: Qingxiang Liu (28-Apr-2019) [ascii option]
! Mod: Qingxiang Liu (10-Jun-2021) [upload to Github]
!
! Email: liuqingxiang@ouc.edu.cn; qingxiang.liu@unimelb.edu.au
!/ ------------------------------------------------------------------- /
subroutine HasNaN(nk, nth, arr2D, flag)
!/
!/  IsNaN check (added on 24-Apr-2019)
!/
    implicit none
!
    integer, intent(in)  :: nk, nth ! # of freq. & dirc.
    real, intent(in)     :: arr2D(nk, nth)
    logical, intent(out) :: flag
!
! Intialization
    flag = .true.
!
    if ( all(arr2D .ge. -huge(arr2D(1, 1))) .and.  &
         all(arr2D .le.  huge(arr2D(1, 1))) ) then
        flag = .false.
    end if
!/
end subroutine HasNaN

!/ ------------------------------------------------------------------- /
subroutine CountRecords(fnm, nt, nk, nth, FBIN)
!/
    implicit none
!
! Parameter list
    character(len=*), intent(in) :: fnm
    integer, intent(out)         :: nt,  & ! # of time
                                    nk,  & ! # of frequency
                                    nth    ! # of directions
    integer, intent(in)          :: FBIN   ! flag for binary & ascii
                                           ! 0: binary
                                           ! 1: ascii
!
! local vars.
    integer                      :: temp(2), rnum
    real                         :: plon, plat, pdpt
    real, allocatable            :: freq(:), th(:), spec(:, :), tsnl(:, :)
    real                         :: tkurt
    logical                      :: flagE, flagS
!
    if (FBIN .eq. 0) then
! Open binary file (convert option is supported by gnu & intel compilers)
        open(12, file=trim(fnm), form='unformatted', access='stream',  &
             status='old', action='read', convert='big_endian')
! Location of the point
        read(12) plon, plat, pdpt
! f, θ
        read(12) nk, nth
!
    else if (FBIN .eq. 1) then
! ascii
        open(12, file=trim(fnm), form='formatted', status='old',       &
             action='read')
        read(12, *) plon, plat
        read(12, *) pdpt
        read(12, *) nk, nth
    end if
!
    if (allocated(freq)) deallocate(freq); allocate(freq(nk))
    if (allocated(th  )) deallocate(th  ); allocate(th(nth))
    if (allocated(spec)) deallocate(spec); allocate(spec(nk, nth))
    if (allocated(tsnl)) deallocate(tsnl); allocate(tsnl(nk, nth))
! f, θ
    if (FBIN .eq. 0) then
        read(12) freq, th
    else if (FBIN .eq. 1) then
        read(12, *) freq
        read(12, *) th
    end if
!
! Get the number of records
    nt = 0
    do
        if (FBIN .eq. 0) then             ! binary
            read(12, iostat=rnum) temp
            if (rnum == -1) exit
            read(12) tkurt
            read(12) spec
            read(12) tsnl
        else if (FBIN .eq. 1) then
            read(12, *, iostat=rnum) temp
            if (rnum == -1) exit
            read(12, *) tkurt
            read(12, *) spec
            read(12, *) tsnl
        end if
        nt = nt + 1
! NaN check
        call HasNaN(nk, nth, spec, flagE)
        call HasNaN(nk, nth, tsnl, flagS)
! Screen output
        if (flagE .or. flagS) then
            write(*, *)
            write(*, '(A)')                         '|----------------------------------------->'
            write(*, '(A, I5, 4X, I8.8, 4X, I6.6)') '|    time: ', nt, temp
            write(*, '(A)', advance='no')           '|    spec: '
            write(*, '(8E12.4)')                                   spec
            write(*, *)
            write(*, '(A)', advance='no')           '|    snl : '
            write(*, '(8E12.4)')                                   tsnl
            write(*, *)
        end if
    end do
!
! close file
    close(12)
!
! Screen output
    write(*, *)
    write(*, '(A, I5, 3A)') "⊚ → Find a total of ", nt, " records in |", trim(fnm), "|"
    write(*, '(A, I5)')     "⊚ → # of frequency :", nk
    write(*, '(A, I5)')     "⊚ → # of direction :", nth
    write(*, *)
!/
end subroutine CountRecords

!/ ------------------------------------------------------------------- /
subroutine ReadNL5SRC(fnm, nt, nk, nth, FBIN, plon, plat, pdpt, &
                      freq, th, time, kurt, efth, snl)
!/
    implicit none
!
! Parameter list
    character(len=*), intent(in) :: fnm
    integer, intent(in)          :: nt, nk, nth
    integer, intent(in)          :: FBIN
    real, intent(out)            :: plon
    real, intent(out)            :: plat
    real, intent(out)            :: pdpt
    real, intent(out)            :: freq(nk)
    real, intent(out)            :: th(nth)
    integer, intent(out)         :: time(2, nt)
    real, intent(out)            :: kurt(nt)
    real, intent(out)            :: efth(nt, nk, nth)
    real, intent(out)            :: snl(nt, nk, nth)
!
! local vars.
    integer                      :: it, tnk, tnth
!
    if (FBIN .eq. 0) then
! Open binary file
        open(12, file=trim(fnm), form='unformatted', access='stream',  &
             status='old', action='read', convert='big_endian')
! Location of the point
        read(12) plon, plat, pdpt
! σ, θ
        read(12) tnk, tnth
        read(12) freq, th
! Snl(it, θ, σ)
        do it = 1, nt
            read(12) time(:, it)
            read(12) kurt(it)
            read(12) efth(it, :, :)
            read(12) snl (it, :, :)
        end do
!
    else if (FBIN .eq. 1) then
        open(12, file=trim(fnm), form='formatted',  status='old',      &
             action='read')
! Location of the point
        read(12, *) plon, plat
        read(12, *) pdpt
! σ, θ
        read(12, *) tnk, tnth
        read(12, *) freq
        read(12, *) th
! Snl(it, θ, σ)
        do it = 1, nt
            read(12, *) time(:, it)
            read(12, *) kurt(it)
            read(12, *) efth(it, :, :)
            read(12, *) snl (it, :, :)
        end do
    end if
!
    close(12)
!/
end subroutine ReadNL5SRC
!/ ------------------------------------------------------------------- /
