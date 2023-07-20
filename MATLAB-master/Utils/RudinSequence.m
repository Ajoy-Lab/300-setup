function retval = RudinSequence (m)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    numOfSteps = int16(round(log(m) / log(2)));
    
    if 2^numOfSteps < m
        numOfSteps = numOfSteps + 1;
    end
    
    retval =  zeros(1, 2^numOfSteps);
    retval(1:2) = 1;            
    numOfSteps = numOfSteps - 1;   
    
    for n=1:numOfSteps
        retval((2^n + 1) : (2^n + 2^(n-1))) = retval(1 : 2^(n-1));
        retval((2^n + 2^(n-1)+1) : (2^(n+1))) = -retval(2^(n-1)+1 : 2^n);
    end 
    
    retval = retval(1:m);
    
    retval = -0.5 * pi .* (retval - 1);
end

