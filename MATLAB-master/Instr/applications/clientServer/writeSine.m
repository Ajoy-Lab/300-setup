%% Proteus P2584M Example
% Connect to P2584M instrument, and download waveform.

addpath '..\..\..\Utils';
addpath '..\..\drivers';
addpath '..\..\..\Functions';

sampleRate = 2.5E9;
segLen = 1024;
bits = 16;

openInst;
    
res = inst.SendScpi('*IDN?');
assert(res.ErrCode == 0);
fprintf(1, '\nConnected to ''%s''\n', convertCharsToStrings(char(res.RespStr)));

res = inst.SendScpi('*CLS');
assert(res.ErrCode == 0);

res = inst.SendScpi('*RST');
assert(res.ErrCode == 0);

res = inst.SendScpi(':TRAC:DEL 1');
assert(res.ErrCode == 0);

res = inst.SendScpi(':TRAC:DEL 2');
assert(res.ErrCode == 0);

res = inst.SendScpi(':FREQ:RAST 2.5E9');
assert(res.ErrCode == 0);     

% Set function mode = 'AWG'
res = inst.SendScpi(':TRAC:DEF:TYPE NORM');
assert(res.ErrCode == 0);

res = inst.SendScpi(':INST:CHAN 1');
assert(res.ErrCode == 0);

res = inst.SendScpi(':INIT:CONT ON"');
assert(res.ErrCode == 0);
res = inst.SendScpi(':SOUR:FUNC:MODE ARB');
assert(res.ErrCode == 0);

% Define segment 1 of 1024 points
res = inst.SendScpi(':TRAC:DEF 1, 1024');
assert(res.ErrCode == 0);

dacSignal = sine(sampleRate, 5, 0, segLen, bits); % sclk, cycles, phase, segLen, bits

% select segmen 1 as the the programmable segment
res = inst.SendScpi(':TRAC:SEL 1');
assert(res.ErrCode == 0); 

% Download the binary data to segment 1
res = inst.WriteBinaryData(':TRAC:DATA 0,#', dacSignal);
assert(res.ErrCode == 0);

res = inst.SendScpi(':SOUR:FUNC:SEG 1');
assert(res.ErrCode == 0); 

res = inst.SendScpi(':OUTP ON');
assert(res.ErrCode == 0); 

closeInst;
    