%% Author: Joan <Joan@ARBITRARYLAPTOP>
%% Created: 2020-01-21

%% It returns an IQ modulated carrier from the complex baseband signal

function retval = IqModulator2 (iqWfm, cFreq, cPhase, sRate)
  
  % The carrier complex waveform is calculated
  carrier = 0:(length(iqWfm) - 1);
  carrier = carrier / sRate;
  %
  %carrier = cos(2*pi*cFreq * carrier) + 1i * sin(2*pi*cFreq * carrier);
  carrier = exp(1i * 2 * pi * cFreq * carrier +  cPhase);
  % IQ modulation is performed by complex multiplication of the complex baseband
  % and the complex carrier
  retval = iqWfm .* carrier;
  % The actual modulated waveform is just the real part of the product
  retval = real(retval);

end
