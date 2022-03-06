function sol = SteepestDescent(g,x,tol,N,J,F)
%SteepestDescent Finds Solution to System Using Steepest Descent Method
% g -> Error Function
% x -> Initial Approximation
% tol -> Tolerance For Convergence
% N -> Number of Allowed Iterations
% J -> The Jacobian
% F -> The Cost Function

%Iter 1
k = 1;

%While Below the Max Iterations
while(k<=N)
    
    %Starting Point
    g1 = g(x);
    z = 2*transpose(J(x))*F(x); %Actually Gradient Just Noticed Analytically
    z0 = norm(z,2);
    
    if z0 == 0
        
        sol = x;
        disp("Zero Gradient");
        return;
    end
    
    z = z/z0;
    alpha1 = 0;
    alpha3 = 1;
    g3 = g(x-alpha3*z);
    
    while g3 >= g1
        
        alpha3 = alpha3/2;
        g3 = g(x-alpha3*z);
        
        if alpha3 < tol/2
           
            sol = x;
            disp("No Likely Improvement");
            return;
        end
    end
    
    alpha2 = alpha3/2;
    g2 = g(x-alpha2*z);
    
    h1 = (g2-g1)/alpha2;
    h2 = (g3-g2)/(alpha3-alpha2);
    h3 = (h2-h1)/alpha3;
    
    alpha0 = .5*(alpha2 - h1/h3);
    g0 = g(x-alpha0*z);
    
    alphas = zeros(4,2);
    alphas(:,1) = [alpha0; alpha1; alpha2; alpha3];
    
    for i = 1: 4
        
        alphas(i,2) = g(x-alphas(i,1)*z);
    end
    
    m = min(alphas(:,2));
    alpha = 0;
    
    for i = 1: 4
        
        if alphas(i,2) == m
            
            alpha = alphas(i,1);
            break;
        end
    end
    
       
    gf = g(x-alpha*z);
    x = x - alpha*z;
    
    if abs(gf-g1) < tol
        
        sol = x;
        disp("It worked");
        return;
    end
    
    k = k+1;
end
   
   sol = x;
   disp("Max Iterations Exceeded");
   
end

