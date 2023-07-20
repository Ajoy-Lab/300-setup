function dacSignal = pulse(sclk, segLen, amp, bits, onTime)

  time = -(segLen-1)/2:(segLen-1)/2;
  
  onTime = int32(onTime*segLen);          

  rawSignal = rectpuls(time,onTime);
  
  dacSignal = ampScale(bits, rawSignal);
  
  dacSignal = dacSignal.'; % for MatLab
  %dacSignal = dacSignal; % for Octave
  
end
