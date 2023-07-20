function dacSignal = gaussianPulse(sclk, segLen, amp, bits, sigma, offset)
  
  
  format shortEng;
  time = -(segLen-1)/2:(segLen-1)/2;
  variance=sigma^2; % is the pulse half-duration
  rawPulse=(exp(-(time/segLen).^2/(2*variance)));
  
  rawSignal = amp *(rawPulse);
  
  offsetDacSignal = rawSignal(1:offset);
  rawSignal = rawSignal(offset+1:end);
  rawSignal = [rawSignal, offsetDacSignal];
  
  dacSignal = ampScale(bits, rawSignal);
  
  dacSignal = dacSignal.'; % col to row for Octave
  %dacSignal = dacSignal; % for MatLab

end