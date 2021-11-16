program  main

    use m_pade

    !write(*,*) pade(m, z, n, iw, ih, h)

    complex(8),parameter :: zz = (1.0, 1.0) 
    !OPEN(10, file='output.dat', status='new')
    !WRITE(10,*) test()
    !WRITE(unit=10, fmt=*) test()
!    integer n 
!    integer i
!    complex x,y
!    n = 10
!    do i = 1,n
!        x = i* (0.0,1.0)
!        y = cos(x)
!        write(*,*) y
!    enddo
    
    
    write(*,*) 'pade =', test() , 'original_function=', cos(zz) 


end program


    

