! Module provides a consistent API for unit testing, allowing the backend
! test framework to be switched with no or minimal rewriting of the existing
! unit tests.
module unit_test_base_class
    implicit none
    private

    ! Abstract class
    type, abstract :: base_test_type
        contains
            procedure(base_test_type_assert), deferred :: assert
    end type

    ! Explicit interface for procedures (methods) to be defined by the child class.
    interface
        !TODO(Alex) How does one declare overloads for the child class?
        subroutine base_test_type_assert(this, condition)
            import :: base_test_type
            class(base_test_type), intent (in) :: this
            logical, intent(in) :: condition
        end subroutine

    end interface

end module unit_test_base_class
