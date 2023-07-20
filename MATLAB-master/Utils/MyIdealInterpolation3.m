%% Author: Joan <Joan@ARBITRARYLAPTOP>
%% Created: 2020-01-22

function retval = MyIdealInterpolation2 (myArray, xFactor, quality, bwFraction)
  
  %expansion by zero-padding
  retval = zeros(1, xFactor * length(myArray));
  retval([1:xFactor:end]) = myArray;
  % "Ideal" Interpolation filter
  lenSinc = quality; 
  
  mySinc = -lenSinc * xFactor : 1 : lenSinc * xFactor ;  
  mySinc = sinc(bwFraction .* mySinc / xFactor);
  %The FIR taps must be normalized so the DC response is 0dB
  normFactor = sum(mySinc(1:xFactor:length(mySinc)));
  normFactor = 1/ normFactor;
  mySinc = normFactor .* mySinc;
  myWindow = blackman(length(mySinc));
  myWindow = myWindow.'; 
  mySinc = mySinc .* myWindow;
  %convolution
  retval = cconv(retval, mySinc, length(retval));
  %retval = real(retval);
  
end
