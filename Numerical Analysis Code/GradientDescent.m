function params = GradientDescent(grad, observation, parameter0, N, Tol)

%grad: The User must pass in the gradient of the function they are using.
%The gradient should have one vector input containing parameters and 
%another for 
%observation: This is the observed data in a 2xn matrix. The first column
%are the observed x factors and the second column is the observed y values
%(y hat and x hat). 
%parameter0: This is our initial guess at paramters by which we update each
%iteration.
%N: Max number of iterations allowed. 
%Tol: Tolerance for when we can quit iterating. 


k = 1;

alpha = .0000000001;

while(k<=N)
    
    %Calculate MSE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Initialize MSE
    mse = 0;
    
    %Sum MSE 
    for i = 1:length(observation)
        [a, b, c] = grad(parameter0,observation(i,:));
        r = [a; b; c];
        mse = mse + r;
    end
    
    %Normalize by num of points
    mse = mse * 1/length(observation);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Check For Convergence%%%%%%%%%%%%%%%%%%%%%%
    if(norm(mse)<=Tol)
        
        params = parameter0;
        return
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Update Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Otherwise keep updating parameters
    parameter0 = parameter0 + alpha * -mse; 
    
    k = k +1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


params = parameter0;

return

end

