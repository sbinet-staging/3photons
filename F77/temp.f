      program temp
      implicit none
c     Purpose: to compute and display Hex value of NaN
      double precision not_a_number

      not_a_number = sqrt (-1d0)
      write (*, *) not_a_number
      write (*, '(Z16)') not_a_number

      stop
      end
