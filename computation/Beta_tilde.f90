double precision function Beta_tilde(k,l,m,aratio)

     implicit none
     integer             :: k,l,m
     double precision    :: aratio
     
     double precision    :: twooverpi, q, exp_q
     
     parameter(twooverpi = 0.63661977236758134d0)
     
     q = -3.141592653589793238d0 / aratio
     exp_q = dexp(q)
     
     
     if ( (k .eq. m-l) .or. (k .eq. l-m) .or. (k .eq. l+m) ) then
          Beta_tilde = 0.0d0
     else
          Beta_tilde = -2.0d0/q * ( &
                              (exp_q**m - exp_q**(k+l))/(k+l-m) &
                              + (exp_q**l - exp_q**(k+m))/(k-l+m) &
                              + (exp_q**(l+m) - exp_q**k)/(k-l-m) &
                              + (1.0d0 - exp_q**(k+l+m))/(k+l+m) &
                              )
     end if

end function Beta_tilde
