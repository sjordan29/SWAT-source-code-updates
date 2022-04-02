      real*8 function interp_lin (x,y,x_interp) result (y_interp)
        real*8, dimension(:), allocatable, intent(in) :: x,y
        real*8, intent(in):: x_interp
        real*8 :: yyout
        !!real*8 :: y_interp

        real*8:: mult

        integer :: i_val, j, k, dim 
        real*8 :: delta, min_d, a, b
        real*8 :: xxmn, xxmn2, xxmx, xxmx2
        real*8 :: yymn, yymn2, yymx, yymx2
        !! print*, "x array", x
        !! print*, "y array", y 

        dim = size(x)
        !! print*, "dim", dim
        delta = 0.
        min_d = 9999999999999999.
        j = -99. 
        y_interp = 0.


        xxmn = minval(x, dim=1)
        xxmn2 = minval(x(2:dim-1))
        xxmx = maxval(x, dim=1)
        xxmx2 = maxval(x(2:dim-1), dim=1)
        yymn = minval(y,dim=1)
        yymn2 = minval(y(2:dim), dim=1)
        yymx = maxval(y, dim=1)
        yymx2 = maxval(y(2:dim-1), dim=1)

        !! print*, "Starting interp function"
        !! print*, "xinterp", x_interp
        !! print*, "x min/max", xxmn, xxmx
        !! print*, "x 1 in", xxmn2, xxmx2
        !! print*, "y min/max", yymn, yymx
        !! print*, "y 1 in", yymn2, yymx2

        !! extreme cases extrapolation
        if (x_interp <= xxmn) then
            yyout = ((yymn2 - yymn) / (xxmn2 - xxmn)) * (x_interp - xxmn) + yymn

        elseif (x_interp >= xxmx) then
            mult = ((yymx - yymx2 ) / (xxmx - xxmx2)) * (x_interp - xxmx)
            yyout = yymx + mult 

        else
          do i_val = 1,dim
            if(x(i_val) == x_interp) then 
                yyout = y(i_val)
                exit 
            else 

              delta = ABS(x(i_val) - x_interp)

              if(delta < min_d) then
                min_d = delta 
                j = i_val
              end if 

            end if
            
          end do


              k = 0
    
              if(x(j) < x_interp) then
                  k = j
              else
                  k = j-1
              end if

    
              a = (y(k+1) - y(k)) / (x(k+1) - x(k))
              b = y(k) - a*x(k)
              yyout= a*x_interp + b

        end if
        y_interp = yyout 
      return
      end