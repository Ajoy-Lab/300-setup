%% Author: Joan <Joan@ARBITRARYLAPTOP>
%% Created: 2020-01-22

function retval = IqIdealInterpolation (myArray, xFactor)
  myArray = myArray.';
  %expansion by zero-padding
  retval = zeros(1, xFactor * length(myArray));
  retval([1:xFactor:end]) = myArray;
  % "Ideal" Interpolation filter
  lenSinc = 40;  
  mySinc = -lenSinc * xFactor : 1 : lenSinc * xFactor ;  
  mySinc = sinc(mySinc / xFactor);
  myWindow = blackman(length(mySinc));
  myWindow = myWindow.'; 
  mySinc = mySinc .* myWindow;
  %convolution
  retval = cconv(retval, mySinc, length(retval));
  
  %retval = real(retval);
  retval = retval.';
end
