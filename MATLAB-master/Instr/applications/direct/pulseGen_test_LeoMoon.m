%% Clear everything
clear;
close;

%% Set default variables
% global sampleRateDAC
% sampleRateDAC = 1e9;
global inst
global bits
bits = 16;

remoteAddr = '192.168.10.5'; % new computer
remotePort = 2020;
localPort = 9090;

dll_path = 'C:\\Windows\\System32\\TEPAdmin.dll';
asm = NET.addAssembly(dll_path);
import TaborElec.Proteus.CLI.*
import TaborElec.Proteus.CLI.Admin.*
import System.*

admin = CProteusAdmin(@OnLoggerEvent);
rc = admin.Open();
assert(rc == 0);


try
    slotIds = admin.GetSlotIds();
    numSlots = length(size(slotIds));
    assert(numSlots > 0);

    % If there are multiple slots, let the user select one ..
    sId = slotIds(1);
    if numSlots > 1
        fprintf('\n%d slots were found\n', numSlots);
        for n = 1:numSlots
            sId = slotIds(n);
            slotInfo = admin.GetSlotInfo(sId);
            if ~slotInfo.IsSlotInUse
                modelName = slotInfo.ModelName;
                if slotInfo.IsDummySlot
                    fprintf(' * Slot Number:%d Model %s [Dummy Slot].\n', sId, modelName);
                else
                    fprintf(' * Slot Number:%d Model %s.\n', sId, modelName);
                end
            end
        end
        pause(0.1);
        choice = 8%input('Enter SlotId ');
        fprintf('\n');
        sId = uint32(choice);
    end

    % Connect to the selected instrument ..
    should_reset = true;
    inst = admin.OpenInstrument(sId, should_reset);
    instId = inst.InstrId;
        
    catch ME
        admin.Close();
        rethrow(ME) 
end
%%Until above is good

% run_square_pulse(inst);
% run_NCO(inst)
set_marker(inst);

admin.CloseInstrument(inst.InstrId);
admin.Close();



function run_NCO(inst)
global sampleRateDAC
ch = 1;
inst.SendScpi(sprintf(':INST:CHAN %d',ch));
inst.SendScpi([':FREQ:RAST ' num2str(2.5E9)]);
inst.SendScpi(':SOUR:INT X8');
inst.SendScpi(':MODE NCO');
% inst.SendScpi(':IQM ONE');
% sampleRateDAC = sampleRateInterp;
inst.SendScpi(sprintf(':FREQ:RAST %d', 5400e6));
inst.SendScpi('SOUR:NCO:SIXD1 ON');
inst.SendScpi(sprintf(':SOUR:NCO:CFR1 %d',75.38E4));
inst.SendScpi(sprintf(':SOUR:NCO:PHAS1 %d',90));
resp = inst.SendScpi(':OUTP ON');

end

function set_marker(inst)
    max_sampleRateDAC = 9e9;
%     sampleRateDAC = 1.125e9; % 9e9 / 8
    dac_res = 16;
    granularity = 32;
    half_dac = floor((2^16 - 1) / 2);
    %% Create a marker waveform that has 
%     marker_on_pts = granularity*round(sampleRateDAC * marker_on/granularity);
%     marker_off_pts = granularity*round(sampleRateDAC * marker_off/granularity);
    marker_on_pts = 6400;
    marker_off_pts = 6400;
    signal_on_pts = marker_on_pts * 2;
    signal_off_pts = marker_off_pts * 2;
    %% turn on marker
    temp_mrkr_on = uint8(zeros(1, marker_on_pts) + 1)*half_dac;
    temp_mrkr_off = uint8(zeros(1, marker_off_pts));
    temp_on = uint8(zeros(1, signal_on_pts));
    temp_off = uint8(zeros(1, signal_off_pts));
    myWfm = [temp_on temp_off];
    myMkr = [temp_mrkr_on temp_mrkr_off];
    inst.SendScpi(sprintf(':FREQ:RAST %d', 2.5e9));
    myMkr = myMkr(1:2:length(myMkr)) + 16 * myMkr(2:2:length(myMkr));
    
    %% Download to channel 2 data
    inst.SendScpi(sprintf(':INST:CHAN %d',2));
    inst.SendScpi(':TRAC:FORM U16');
    inst.SendScpi(sprintf(':TRAC:DEF %d, %d', 1, length(myWfm)));
    inst.SendScpi(sprintf(':TRAC:SEL %d',1));
%     prefix = ':TRAC:DATA 0,';
%     myWfm = uint16(myWfm);
%     myWfm = typecast(myWfm, 'uint8');
%     res = inst.WriteBinaryData(prefix, myWfm);
    
    %% Setting into a segment in channel 2
    inst.SendScpi(sprintf(':TRAC:SEL %d', 1));
    res = inst.WriteBinaryData(':MARK:DATA 0,',myMkr);
    inst.SendScpi(sprintf(':MARK:SEL %d', 1));
    inst.SendScpi(':MARK:VOLT:PTOP 0.5');
    inst.SendScpi(':MARK:VOLT:OFFS 0.0');
    inst.SendScpi(':MARK:STAT ON');
    assert(res.ErrCode == 0); 
    
    %% creating task table for infinite repetition
    inst.SendScpi(sprintf(':INST:CHAN %d', 2));
    inst.SendScpi('TASK:ZERO:ALL');
    inst.SendScpi(sprintf(':TASK:COMP:LENG %d',1));

    inst.SendScpi(sprintf(':TASK:COMP:SEL %d',1));
    inst.SendScpi(':TASK:COMP:TYPE SING');
    inst.SendScpi(sprintf(':TASK:COMP:SEGM %d',1));
    inst.SendScpi(':TASK:COMP:ENAB NONE');
    inst.SendScpi(sprintf(':TASK:COMP:LOOP %d',10^6));
    inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',1));

    inst.SendScpi('TASK:COMP:WRITE');
    inst.SendScpi('SOUR:FUNC:MODE TASK');
    inst.SendScpi(':OUTP ON');
end

function run_square_pulse(inst)
    % maximum sampling rate for IQ One mode is 2.5Gsa/sec for the baseband
    % waveform(not NCO!) -- because IQ One mode has I and Q interleaved,
    % the complex waveform has sample rate of 1.25GSa/sec
    max_sampleRateDAC = 9e9;
    sampleRateDAC = 1.125e9;
    % resolution for DAC
    dac_res = 16;
    % granularity of a waveform (Waveform length must be 32)
    granularity = 32;
    res = inst.SendScpi('*IDN?');
    assert(res.ErrCode == 0);

    fprintf(1, '\nConnected to ''%s''\n', netStrToStr(res.RespStr));

    pause(0.01);

    res = inst.SendScpi('*CLS'); % clear
    assert(res.ErrCode == 0);

    res = inst.SendScpi('*RST'); % reset
    assert(res.ErrCode == 0);
    fprintf('Reset complete\n');

    %% Initialization of the AWG
    % Sets the voltage value to maximum -- not sure why we need this line though
    inst.SendScpi(':SOUR:VOLT MAX');
    % initializes to the continous mode where no triggers are needed
    inst.SendScpi(':INIT:CONT ON');
    % delete all segments stored previously
    res = inst.SendScpi(':TRAC:DEL:ALL');
    assert(res.ErrCode == 0);

    max_dac = 2^dac_res - 1;
    half_dac = floor(max_dac / 2);
    
    pulse_on_len = 60e-6;
    pulse_off_len = 40e-6;

    pulse_on_pts = granularity*round(sampleRateDAC * pulse_on_len/granularity);
    pulse_off_pts = granularity*round(sampleRateDAC * pulse_off_len/granularity);
    dacWaveI_on = (zeros(1, pulse_on_pts) + 1) * max_dac;
    fprintf("%d, %d \n", pulse_off_pts, uint32(pulse_off_pts));
    dacWaveI_off = (zeros(1, pulse_off_pts)+1)* half_dac;
    dacWaveQ_on = (zeros(1, pulse_on_pts)+ 1) * max_dac;
    dacWaveQ_off = (zeros(1, pulse_off_pts)+1)* half_dac;
    dacWaveI = [dacWaveI_on dacWaveI_off];
    dacWaveQ = [dacWaveQ_on dacWaveQ_off];
    dacWaveIQ = [dacWaveI ; dacWaveQ];
    dacWaveIQ = dacWaveIQ(:)';
    ch = 1;
    inst.SendScpi(sprintf(':INST:CHAN %d',ch));
    inst.SendScpi(':TRAC:FORM U16');
    inst.SendScpi(sprintf(':TRAC:DEF %d, %d', 1, length(dacWaveIQ)));
    inst.SendScpi(sprintf(':TRAC:SEL %d',1));
    res = inst.SendScpi("TRAC:DEF?");
    prefix = ':TRAC:DATA 0,';
    myWfm = uint16(dacWaveIQ);
    myWfm = typecast(myWfm, 'uint8');
    res = inst.WriteBinaryData(prefix, myWfm);
    
    inst.SendScpi(sprintf(':INST:CHAN %d',ch));
    inst.SendScpi('TASK:ZERO:ALL');
    inst.SendScpi(sprintf(':TASK:COMP:LENG %d',1));

    inst.SendScpi(sprintf(':TASK:COMP:SEL %d',1));
    inst.SendScpi(':TASK:COMP:TYPE SING');
    inst.SendScpi(sprintf(':TASK:COMP:SEGM %d',1));
    inst.SendScpi(':TASK:COMP:ENAB NONE');
    inst.SendScpi(sprintf(':TASK:COMP:LOOP %d',10^6));
    inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',1));

    inst.SendScpi('TASK:COMP:WRITE');
    inst.SendScpi(':INST:CHAN 1');
    inst.SendScpi(sprintf(':FREQ:RAST %d', sampleRateDAC));
    inst.SendScpi(':SOUR:INT X8');
    inst.SendScpi(sprintf(':FREQ:RAST %d', max_sampleRateDAC));
    inst.SendScpi(':MODE DUC');
    inst.SendScpi(':IQM ONE');
    inst.SendScpi(sprintf(':SOUR:NCO:CFR1 %d',50E4));
    inst.SendScpi(sprintf(':SOUR:NCO:PHAS1 %d',90));
    inst.SendScpi('SOUR:NCO:SIXD1 ON')
    inst.SendScpi('SOUR:FUNC:MODE TASK');
    inst.SendScpi(':OUTP ON');
    
end

function str = netStrToStr(netStr)
    try
        str = convertCharsToStrings(char(netStr));
    catch
        str = '';
    end
end

function dacRes = getDacResolution(inst)
%Proteus unit internal data format. The unit used in 300 setup is 16 bit
    dacRes = inst.SendScpi(':TRAC:FORM?');
    dacRes = strtrim(netStrToStr(dacRes.RespStr));
    if contains(dacRes, 'U8')
        dacRes = 8;
    else
        dacRes = 16;
    end
end

%granularity of the Proteus model is 32

function maxSr = getMaxSamplingRate(inst)
%This function returns 9e-9 as expected
    maxSr = inst.SendScpi(':FREQ:RAST MAX?');
    maxSr = strtrim(netStrToStr(maxSr.RespStr));
    maxSr = str2double(maxSr);
end