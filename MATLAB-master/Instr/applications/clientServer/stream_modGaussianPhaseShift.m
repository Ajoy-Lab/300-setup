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

setCh(ip, 1)
segDef(ip, 1, segLen);
useSeg(ip, 1)
output(ip, 1)
markers
%setCW(ip, 400E6)

for c = 0:1:100
    %tic
    dacSignal_1 = modGaussianPulse(sampleRate, 25, segLen/2, 1, bits, 0.1, 0);
    dacSignal_2 = modGaussianPulse(sampleRate, 25, segLen/2, 1, bits, 0.1, c);
    dacSignal = [dacSignal_1; dacSignal_2];
    plot(dacSignal);
    dataWrite (ip, dacSignal);
    %toc
    
end

instClose