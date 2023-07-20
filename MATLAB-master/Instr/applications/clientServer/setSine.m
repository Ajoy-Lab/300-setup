addpath '..\..\..\Utils';
addpath '..\..\drivers';
addpath '..\..\..\Functions';

ip = "192.168.0.197";
%ip = "76.102.69.4";
sampleRate = 1e9
bits = 16
segLen = 4096
trace = 1;

instInit(ip, 1, 1, sampleRate) %ip, openPXI, doReset, setSample
%instInit(ip, 0, 0, 0) %ip, openPXI, doReset, setSample
format shortg
segDef(ip, trace, segLen);
setCh(ip, 1)
useSeg(ip, trace)
output(ip, 1)

markers
%setCW(ip, 900E6)
%Frequency = Cycles x Sample Rate / Length  
scpiWrite (ip, ":TRACe:SEL " + trace);
dacSignal = sine(sampleRate, 1, 0, segLen, bits); % sclk, cycles, phase, segLen, bits
dataWrite (ip, dacSignal);   
plot(dacSignal);

dataWrite (ip, dacSignal);

%instClose