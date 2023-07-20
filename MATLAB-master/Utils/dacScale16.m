
maxSig = max(waveformReSamp);
minSig = min(waveformReSamp);
avgSig = mean(waveformReSamp);
signalRangePkPk = (minSig*-1)+maxSig;

dacBits = 16;
verticalScale = ((2^dacBits)/2)-1;

vertScaled = (waveformReSamp / maxSig) * verticalScale;
dacSignal = int32(vertScaled + verticalScale);
disp(fileName);
%Save
csvwrite (fileName, dacSignal);