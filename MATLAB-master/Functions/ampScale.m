
function dacSignal = ampScale(bits, rawSignal)
 
  maxSig = max(rawSignal);
  verticalScale = ((2^bits)/2)-1;

  vertScaled = (rawSignal / maxSig) * verticalScale;
  dacSignal = uint8(vertScaled + verticalScale);
  %plot(dacSignal);

  if bits > 8
      dacSignal16 = [];
      sigLen = length(dacSignal);
      k=1;
      for j = 1:2:sigLen*2;
        dacSignal16(j) = bitand(dacSignal(k), 255);
        dacSignal16(j+1) = bitshift(dacSignal(k),-8);
        k = k + 1;
      end
      dacSignal = dacSignal16;
  end
  
  dacSignal = typecast(dacSignal, 'uint8');
end
