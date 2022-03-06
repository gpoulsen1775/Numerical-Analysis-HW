 function sol = NMN(F,J,x0,tol,N) 
 
      k = 1;
      sol = x0;
      
      while k <= N
          
          y = J(sol)\(-F(sol));
          sol = sol + y;
          
          if(norm(y) < tol)
              
              break;
          end
          
          k = k + 1;
      end
 end