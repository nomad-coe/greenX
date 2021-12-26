! COPYRIGHT
!  Copyright (C) 2020-2021 Green-X library (MA, MG, XG)
!  This file is distributed under the terms of the
!  APACHE2 License.

program  main

    use pade_approximant
    use test
    
    complex(dpc),parameter :: xx = (1.0,1.0) 
    complex(dpc), parameter :: w = (2.0,2.0)

    write(*,*) 'pade1 =', test1(), 'original_function1=', -1./(xx-w) 
    write(*,*) 'pade2 =', test2(), 'original_function2=', -1./sqrt(xx-w) 
end program


    

