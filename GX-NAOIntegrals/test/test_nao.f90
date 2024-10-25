program test 

    use wigner, only: threej_table_init, threej_lookup
    use gaunt, only: r_gaunt

    implicit none 

    real(kind=8) :: a 

    !call threej_table_init(2, 40)
    !a =  threej_lookup(4, 4, 4, 0, 0, 0)

    a = r_gaunt(10, 10, 10, 0, 0, 0)

    print *, a 
end program test 