

%% ADC Example
% Capture Data with ADC

clear all;
close all;
clc;

%% Load TEPAdmin.dll which is a .Net Assembly
% The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32 

% Currently the tested DLL is installed in the following path
asm = NET.addAssembly('C:\Windows\System32\TEPAdmin.dll');

import TaborElec.Proteus.CLI.*
import TaborElec.Proteus.CLI.Admin.*

admin = CProteusAdmin(@OnLoggerEvent);
rc = admin.Open();
assert(rc == 0);

%% Use the administrator, and close it at the end
try
    slotIds = admin.GetSlotIds();
    numSlots = length(slotIds);
    assert(numSlots > 0);
    
    % If there are multiple slots, let the user select one ..
    sId = slotIds(1);
    if numSlots > 1
        fprintf(1, '\n%d slots were found\n', numSlots);
        for n = 1:numSlots
            sId = slotIds(n);
            slotInfo = admin.GetSlotInfo(sId);
            if ~slotInfo.IsSlotInUse
                modelName = slotInfo.ModelName;
                if slotInfo.IsDummySlot
                    fprintf(1, ' * Slot Number: Model %s [Dummy Slot].\n', sId, ModelName);
                else
                    fprintf(1, ' * Slot Number: Model %s.\n', sId, ModelName);
                end
            end
        end
        choice = input('Enter SlotId ');
        fprintf(1, '\n');
        sId = uint32(choice);
    end
    
    % Connect to the selected instrument ..
    should_reset = true;
    inst = admin.OpenInstrument(sId, should_reset);
    instId = inst.InstrId;
    %% 
    
    % ---------------------------------------------------------------------
    % Send SCPI commands
    % ---------------------------------------------------------------------
    
    res = inst.SendScpi('*IDN?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nConnected to ''%s''\n', netStrToStr(res.RespStr));
    
    res = inst.SendScpi('*CLS'); % clear
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi('*RST'); % reset
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi(':FREQ:RAST 1E9'); % set sample clock
    assert(res.ErrCode == 0);    
    
    fprintf('Reset complete\n');

    % ---------------------------------------------------------------------
    % Download sine waveform of 1024 points to segment 1 of channel 1
    % ---------------------------------------------------------------------

    res = inst.SendScpi(':INST:CHAN 1'); % select channel 1
    assert(res.ErrCode == 0);
    
    segLen = 40960;
    segLenDC = 50048;
    bits = 8;
    cycles = 800;
    
    dacSignalMod = sine(cycles, 0, segLen, bits); % sclk, cycles, phase, segLen, amp, bits, onTime(%)
    
%     dacSignalDC = [];
%     for i = 1 : segLenDC
%       dacSignalDC(i) = 127;
%     end
    
%     dacSignal = [dacSignalDC dacSignalMod];

    dacSignal = dacSignalMod;

    % Define segment 1 
    res = inst.SendScpi(strcat(':TRAC:DEF 1,',num2str(length(dacSignal)))); % define memory location 1 of length dacSignal
    assert(res.ErrCode == 0);

    % select segmen 1 as the the programmable segment
    res = inst.SendScpi(':TRAC:SEL 1');
    assert(res.ErrCode == 0); 

    % Download the binary data to segment 1
    res = inst.WriteBinaryData(':TRAC:DATA 0,#', dacSignal);
    assert(res.ErrCode == 0);
    
    % ---------------------------------------------------------------------
    % Play segment 1 in channel 1
    % ---------------------------------------------------------------------

    res = inst.SendScpi(':INST:CHAN 1');
    assert(res.ErrCode == 0);    

    res = inst.SendScpi(':SOUR:FUNC:MODE:SEGM 1');
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':SOUR:VOLT 0.2');
    assert(res.ErrCode == 0);
    
    %Enable output
    res = inst.SendScpi(':OUTP ON');
    assert(res.ErrCode == 0);

    fprintf('Waveform generated and playing\n');
    
    onTime = 4096;
    offTime = int16(length(dacSignal)/8)-onTime;
    markerWave = [];
    for i = 1 : onTime
        markerWave(i) = 1;
    end
    for i = onTime : offTime;
        markerWave(i) = 0;
    end
    
    % Download the binary data to segment 1
    res = inst.WriteBinaryData(':MARK:DATA 0,#', uint8(markerWave));
    assert(res.ErrCode == 0);

    res = inst.SendScpi(":MARK:VOLT:PTOP 1");
    assert(res.ErrCode == 0);
    res = inst.SendScpi(":MARK:SEL 1");
    assert(res.ErrCode == 0);
    res = inst.SendScpi(":MARK ON");
    assert(res.ErrCode == 0);

    % ---------------------------------------------------------------------
    % Operate the ADC with direct functions 
    % ---------------------------------------------------------------------

    sampleRateADC = 1e9;

    readLen = 16896; % integer div of 96
    readLenSeconds = readLen/1e9;
    numOfFramesToRead=1;
    chanIndex = 0; %CH 1 of Digitizer
    off = 0;
    on = 1;
    numberOfFrames = 749500;
    frameNumber = 10;
    netArray = NET.createArray('System.UInt16', readLen*numOfFramesToRead);

    % Turn on ADC dual-channels mode (state = 1)
    rc = inst.SetAdcDualChanMode(on);
    assert(rc == 0);

    rc = inst.SetAdcSamplingRate(sampleRateADC);
    assert(rc == 0);
    
    rc = inst.SetAdcDualChanMode(on);
    assert(rc == 0);
    
    rc = inst.SetAdcFullScaleMilliVolts(chanIndex,1000); % mv - 500, 800, 1000
    assert(rc == 0);
    
    rc = inst.SetAdcAcquisitionEn(1,1);
    assert(rc == 0);
    
    rc = inst.SetAdcFramesLayout(numberOfFrames, readLen);
    assert(rc == 0);
    
      
    fprintf('ReadLen %d Number of frames %d\n', readLen, numberOfFrames);
    fprintf('Capture size set %d \n', numberOfFrames*readLen);
    
    fprintf('Capture size set\n');
 
%     fprintf('Trigger set - ');
%     trigSource = 1; % 2 = Self, 1 = External
    
    trigSource = 1;%input('Select trigger source - enter 1 for External or 2 for Self :');
    
    rc = inst.SetAdcCaptureTrigSource(chanIndex, trigSource); % 1 means external trigger, 2 means self trigger
    assert(rc == 0); 

    if trigSource == 1
        rc = inst.SetAdcExternTrigPattern(0);
        assert(rc==0);
        %Enable Trig gate mode
        rc = inst.SetAdcExternTrigGateModeEn(1);
        assert(rc==0);
        %Set trigger polarity 0-positive 1-negative
        rc = inst.SetAdcExternTrigPolarity(1);
        assert(rc==0);
        %Set Trig level to 2V
        rc = inst.SetAdcExternTrigThresh(0,1);
        assert(rc==0);
        prompt= 'Set the required trigger delay ';
        %Set trigger delay
        rc = inst.SetAdcExternTrigDelay(0.0000172);
        assert(rc==0);       
    end
   
    fprintf('Aquisition Enabled, starting meas......\n');
    
    rc = inst.SetAdcCaptureEnable(on);
    assert(rc == 0);
    
    % Wait until the capture completes
    status = inst.ReadAdcCompleteFramesCount();
    
    while status ~= numberOfFrames
        pause(0.01);
        status = inst.ReadAdcCompleteFramesCount();
    end
 
    prompt = 'Wait for shuttle to complete ';
    x = input(prompt);    
    fprintf('Reading Frames....\n'); 
    % <param name="adcChanInd">The zero-based index ADC channel.</param>
	% <param name="frame1st">The frame-number of the 1st frame.</param>
	% <param name="numFrames">The number of frames (-1: till the last).</param>
	% <param name="outBuf">Array for the data.</param>
    peakSamp = [];
    
    pulseAmp = [];
    power_of_2 = floor(log2(readLen)); % this is normally 15 for 32us captures
    padded_len= 2^(power_of_2) ;%2^15;
    dF = sampleRateADC/padded_len; %set the discretization of freq in terms of sampleRate
    f = -sampleRateADC/2:dF:sampleRateADC/2-dF; %sampleRate sets the 'bandwidth'

    [v,b1]=min(abs(f-20e6)); %picks out the 20MHz component
    [v,b2]=min(abs(f+20e6));

    
    for i = 1 : numberOfFrames
        rc = inst.ReadMultipleAdcFrames(chanIndex, i, numOfFramesToRead , netArray);
        assert(rc == 0);
        samples = double(netArray);
        peakSamp(i) = max(samples);
        if mod(i,50000)==0
            figure(2);clf;
            plot(samples);
            hold on;
            yline(2048);
        end
        pulse = samples;
        pulseMean = mean(pulse);
        pulseDC = pulse - pulseMean; % remove DC
        X = fftshift(fft(pulseDC,padded_len)); % perform FFT and zero-pad to the padded_len
        linFFT = (abs(X)/readLen);
%                         [amp, loc] = max(linFFT);
        amp=(linFFT(b1) + linFFT(b2));
%             linFFT_vec(:,i)=linFFT;
        pulseAmp(i) = amp;
    end
    
    ivec=1:numberOfFrames;
    time_cycle=40*1e-6+(4+2+2+32)*1e-6;
    time_axis=time_cycle.*ivec;
    %drop first point
    time_axis(1)=[];pulseAmp(1)=[];
    
    figure(1);clf;
    plot(time_axis,(smooth(pulseAmp,10)),'r-');
    set(gca,'ylim',[0 max(pulseAmp)*1.1]);  
    

%     plot(samples);
%     axis([0 length(samples) 0 2^12]);
%     yline(2047);

    pause(0.1);

    rc = inst.SetAdcCaptureEnable(off);
    assert(rc == 0);
        
    %Wipe frames
    rc = inst.WipeMultipleAdcFrames(0,1,-1,0);
    assert(rc == 0);  
    
    % Free the memory space that was allocated for ADC capture
    rc = inst.FreeAdcReservedSpace();
    assert(rc == 0);
   
    fprintf('ADC mem cleared\n');

% ---------------------------------------------------------------------
% End of the example
% ---------------------------------------------------------------------
    
    res = inst.SendScpi(':SYST:ERR?');
    fprintf(1, '\nEnd of Example - %s\n', netStrToStr(res.RespStr));
    % It is recommended to disconnect from instrumet at the end
    rc = admin.CloseInstrument(instId);    
    
    % Close the administrator at the end ..
    admin.Close();
catch ME
    admin.Close();
    rethrow(ME)
end

% Function netStrToStr
function str = netStrToStr(netStr)
    try
        str = convertCharsToStrings(char(netStr));
    catch        
        str = '';
    end
end

%% Function - On Logger-Event
function OnLoggerEvent(~, e)
    try       
        % print only:
        % TaborElec.Proteus.CLI.LogLevel.Critical (1)
        % TaborElec.Proteus.CLI.LogLevel.Error    (2)
        % TaborElec.Proteus.CLI.LogLevel.Warning  (3)
        if int32(e.Level) >= 1 &&  int32(e.Level) <= 3
            msg = netStrToStr(e.Message);
            if strlength(msg) > 0
                fprintf(2, '\n ~ %s\n', msg(1));
            end
        end
        System.Diagnostics.Trace.WriteLine(e.Message);
    catch ME
       rethrow(ME) 
    end
end
%% Signall Processing Functions

%Create a sine wave
function sine = sine(cycles, phase, segLen, bits)
   
  verticalScale = ((2^bits))-1; % 2584 16 bits , 9082 8 bits
  time = -(segLen-1)/2:(segLen-1)/2;
  omega = 2 * pi() * cycles;
  rawSignal = sin(omega*time/segLen); 
  %rawSine = amp* cos(omega*time/segLen); 
 
  sine = ampScale(bits, rawSignal);
  
  %plot(sine);
  
end

% Sacle to FSD

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
end
