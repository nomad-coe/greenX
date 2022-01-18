module zofu_wrapper
    ! External test framework, Zofu. https://github.com/acroucher/zofu
    use zofu, only: unit_test_type

    use unit_test_base_class, only: base_test_type
    use constants, only: dp

    implicit none
    private

    !> Unit-testing Zofu wrapper for serial code
    type, public, extends(base_test_type) :: zofu_wrapper_type
        !> Test object of Zofu's `unit_test_type` type
        type(unit_test_type) :: test
        contains
            !> Assert two scalars are the same, or a logical is true
            procedure :: assert => assert_logical, assert_integer, assert_real, assert_complex
    end type zofu_wrapper_type


contains
    ! Procedures for our Zofu wrapper class, which are themselves wrappers
    ! of the procedures associated with Zofu's test type

    ! Zofu unit test type does not need initialising
    !subroutine init(this)
    !    class(zofu_wrapper_type), intent(in) :: this
    !end subroutine

    !> Wrapper for logical assert
    subroutine assert_logical(this, condition, msg)
        class(zofu_wrapper_type), intent(in) :: this
        !> Logical condition
        logical, intent(in) :: condition
        !> Optional message
        character(len=*), optional, intent(in) :: msg
        if (present(msg)) then
            call this%test%assert(condition, name=msg)
        else
            call this%test%assert(condition)
        end if
    end subroutine

    !> Wrapper for integers-equal assertion
    subroutine assert_integer(this, x, y, msg)
        class(zofu_wrapper_type), intent(in) :: this
        !> Integers to compare
        integer, intent(in) :: x, y
        !> Optional message
        character(len=*), optional, intent(in) :: msg
        if (present(msg)) then
            call this%test%assert(x, y, name=msg)
        else
            call this%test%assert(x, y)
        end if
    end subroutine

    !> Wrapper for real(dp)-equal assertion
    subroutine assert_real(this, x, y, tol, msg)
        class(zofu_wrapper_type), intent(in) :: this
        !> Reals to compare
        real(dp), intent(in) :: x, y
        !> Tolerance below which two reals are considered the same
        real(dp), intent(in) :: tol
        !> Optional message
        character(len=*), optional, intent(in) :: msg
        if (present(msg)) then
            call this%test%assert(x, y, tol, name=msg)
        else
            call this%test%assert(x, y, tol)
        endif
    end subroutine

    !> Wrapper for complex(dp)-equal assertion
    subroutine assert_complex(this, x, y, tol, msg)
        class(zofu_wrapper_type), intent(in) :: this
        logical, intent(in) :: condition
        !> Complex values to compare
        complex(dp), intent(in) :: x, y
        !> Tolerance below which two complex values are considered the same
        complex(dp), intent(in) :: tol
        !> Optional message
        character(len=*), optional, intent(in) :: msg
        if (present(msg)) then
            call this%test%assert(x, y, tol)
        else
            call this%test%assert(x, y, tol, name=msg)
        end if
    end subroutine

end module zofu_wrapper
