addpath '..\..\..\Utils';
addpath '..\..\drivers';
addpath '..\..\..\Functions';
ip = "192.168.0.197";
%ip = "76.102.69.4";

sampleRate = 2.5e9
bits = 16
segLen = 4096

%instInit(ip, 1, 1, sampleRate) %ip, openPXI, doReset, setSample
instInit(ip, 0, 0, 0) %ip, openPXI, doReset, setSample
format shortg
segDef(ip, 1, segLen); %ip, trace, segmentLength
segDef(ip, 2, segLen); %ip, trace, segmentLength
%dacSignal = pulse(sampleRate, segLen, 1, bits, 0.5); % SCLK, Length, Amplitude, bits, Duty Cycle (%)
%dataWrite (ip, dacSignal);
%plot(dacSignal);
setCh(ip, 1)
useSeg(ip, 1)
output(ip, 1)

setCh(ip, 2)
useSeg(ip, 2)
output(ip, 1)

markers % CH 2 Marker

for duty = 0:0.01:0.9
    duty
    trace=1
    scpiWrite (ip, ":TRACe:SEL " + trace);
    dacSignal = pulse(sampleRate, segLen, 1, bits, duty); % SCLK, Length, Amplitude, bits, Duty Cycle (%)
    dataWrite (ip, dacSignal);
    trace=2
    scpiWrite (ip, ":TRACe:SEL " + trace);
    dacSignal = pulse(sampleRate, segLen, 1, bits, 0.5); % SCLK, Length, Amplitude, bits, Duty Cycle (%)
    dataWrite (ip, dacSignal);
    scpiWrite (ip, ":MARK:SEL 1")
    scpiWrite (ip, ":MARK:STAT ON")
end

%instClose