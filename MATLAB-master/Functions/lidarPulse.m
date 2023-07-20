function modPulse = gaussianPulse(sclk, cycles, segLen, amp, bits, sigma)
  
  verticalScale = ((2^bits))-1
  time = -(segLen-1)/2:(segLen-1)/2;
  omega = 2 * pi() * cycles;
  rawSine = amp* sin(omega*time/segLen); 
  %rawSine = amp* cos(omega*time/segLen); 
  format shortEng;
  (sclk * cycles) / segLen;
  
  variance=sigma^2; % is the pulse half-duration
  rawPulse=(exp(-(time/segLen).^2/(2*variance)));
  
  modRawPulse = rawSine .* rawPulse;
  %modRawPulse = rawSine;
  
  dcShift = modRawPulse + amp;
  modPulseCol = fix(dcShift * verticalScale/2);
  
  modPulse = modPulseCol.'; % col to row for Octave
  %modPulse = modPulseCol; % for MatLab
  
  plot(modPulse);
  
  fileName = 'files\gaussianPulse.csv';
  csvwrite (fileName, modPulse);
  
end