program matrix
    implicit none
    integer :: dim = 10 ! dimension of the input matrix
    complex*16, allocatable :: mat(:,:), mat2(:,:), mat3(:,:)
    double precision, allocatable :: eigenvals(:)
    complex*16, allocatable :: WORK(:)
    double precision, allocatable :: WORK2(:)
    integer :: ret
    integer :: i
    allocate(mat(dim,dim))
    allocate(mat2(dim,dim))
    allocate(mat3(dim,dim))
    allocate(eigenvals(dim))
    ! allocate workspaces for LAPACK subroutines
    allocate(WORK(dim * 2))
    allocate(WORK2(dim * 3 - 2))
    open(1, file = "input.txt", status = "old")
    ! read by columns
    read(1, *) mat
    close(1)
    ! transpose to restore the matrix
    mat = transpose(mat)
    ! save the original matrix
    mat2(:,:) = mat(:,:)
    do i = 1 , dim
        eigenvals(i) = 0 ! initialise the array
    end do
    call ZHEEV ('V', 'L', dim, mat, dim, eigenvals, WORK, dim * 2, WORK2, ret)
    ! print *, ret
    if (ret == 0) then
        open(2, file = "output.txt", status = "replace")
        write (2, *) "Eigenvalues: "
        do i = 1, dim
            write (2, *) "Lambda ", i, " ", eigenvals(i)
        end do
        write (2, *) " "
        write (2, *) "Orthonormal eigenvectors: "
        do i = 1, dim
            write (2, *) "Lambda ", i, " ", mat(i, :)
        end do
        mat3 = matmul(transpose(mat), mat2)
        mat = matmul(mat3, mat)
        write (2, *) " "
        write (2, *) "Diagonalised: "
        do i = 1, dim
            write (2, *) mat(i, :)
        end do
        close(2)
    else
        print *, "Computation Failed !"
    end if
end program matrix