MODULE kinds

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: sp, dp, default_string_length
   PUBLIC :: angs, pi

   INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6, 30)
   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 200)
   INTEGER, PARAMETER :: default_string_length = 100

   REAL(KIND=dp), PARAMETER ::  pi = 3.14159265358979323846264338_dp
   REAL(KIND=dp), PARAMETER ::  angs = 0.52917720859_dp 

END MODULE kinds
