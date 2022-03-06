function [eig,eVec] = invPowMethod (A,x,iterations)

  %Give Initial Guess @ EigenVector
  eVec = x;
  
  %Iterate 
  for i = 1: iterations
    
    %Determing Components of new eVec
    v = A\eVec;
    normV = norm(v);
    
    %Set New eVec
    eVec = v/normV;
  end
  
  %Determine eigenvalue from norm of eVec
  eig = 1/normV;

end
