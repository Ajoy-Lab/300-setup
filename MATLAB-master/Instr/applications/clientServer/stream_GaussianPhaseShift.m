addpath '..\..\..\Utils';
addpath '..\..\drivers';
addpath '..\..\..\Functions';
ip = "192.168.0.197";
%ip = "76.102.69.4";

format shortg
sampleRate = 2.5e9
bits = 16
segLen = 4096

instInit(ip, 1, 1, sampleRate) %ip, openPXI, doReset, setSample
%instInit(ip, 0, 0, 0) %ip, openPXI, doReset, setSample
segDef(ip, 1, segLen);
useSeg(ip, 1)
setCh(ip, 1)
output(ip, 1)
markers
%setCW(ip, 400E6)

for c = 0:1:100
    %tic
    dacSignal_1 = gaussianPulse(sampleRate, segLen/8, 1, bits, 0.1, 0);
    dacSignal_2 = gaussianPulse(sampleRate, segLen/8, 1, bits, 0.1, c);
    dacSignal_3 = gaussianPulse(sampleRate, segLen/8, 1, bits, 0.1, 0);
    dacSignal_4 = gaussianPulse(sampleRate, segLen/8, 1, bits, 0.1, 0);
    dacSignal_5 = gaussianPulse(sampleRate, segLen/8, 1, bits, 0.1, 0);
    dacSignal_6 = gaussianPulse(sampleRate, segLen/8, 1, bits, 0.1, 0);
    dacSignal_7 = gaussianPulse(sampleRate, segLen/8, 1, bits, 0.1, 0);
    dacSignal_8 = gaussianPulse(sampleRate, segLen/8, 1, bits, 0.1, 0);
    dacSignal = [dacSignal_1; dacSignal_2];
    dacSignal = [dacSignal; dacSignal_3];
    dacSignal = [dacSignal; dacSignal_4];
    dacSignal = [dacSignal; dacSignal_5];
    dacSignal = [dacSignal; dacSignal_6];
    dacSignal = [dacSignal; dacSignal_7];
    dacSignal = [dacSignal; dacSignal_8];
    plot(dacSignal);
    dataWrite (ip, dacSignal);
    
    
    %toc
    
end

instClose