program  main

    use m_simple_pade

    
    !write(*,*) pade(n,z,f,zz)

    complex(dpc),parameter :: zz = (1.0,1.0) 


    write(*,*) 'pade =', test(), 'original_function=', cos(zz) 


end program


    

