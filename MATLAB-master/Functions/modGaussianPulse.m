function dacSignal = modGaussianPulse(sclk, cycles, segLen, amp, bits, sigma, offset)
  
  time = -(segLen-1)/2:(segLen-1)/2;
  omega = 2 * pi() * cycles;
  rawSine = amp* sin(omega*time/segLen); 
  %rawSine = amp* cos(omega*time/segLen); 
  format shortEng;
  freq = (sclk * cycles) / segLen
  
  variance=sigma^2; % is the pulse half-duration
  rawPulse=(exp(-(time/segLen).^2/(2*variance)));
  
  rawSignal = amp * (rawSine .* rawPulse);
  
  offsetDacSignal = rawSignal(1:offset);
  rawSignal = rawSignal(offset+1:end);
  rawSignal = [rawSignal, offsetDacSignal];
  
  dacSignal = ampScale(bits, rawSignal);
  
  dacSignal = dacSignal.'; % col to row for Octave
  %dacSignal = dacSignal; % for MatLabcsv';
  %csvwrite (fileName, dacSignal);
  
end