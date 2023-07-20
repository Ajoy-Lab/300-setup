%% Clear everything
clear;
close;

%% Set defaults Vars
savealldata=false;
savesinglechunk=false;
savemultwind=false;

sampleRate = 1000e6;
sampleRateDAC = 9e9;
adcDualChanMode = 1;
fullScaleMilliVolts =1000;
trigSource = 1; % 1 = external-trigger
dacChanInd = 1; % 1 = chan 1 and 2 = chan 2
adcChanInd = 0; % ADC Channel 1
measurementTimeSeconds = 7; %Integer
%delay = 0.0000178; % dead time
%delay = 0.0000348; % dead time
%delay=0.00000148;
%delay = 0.0000108; % dead time
%delay = 0.0000028; % dead time
delay = 0.0000038; % dead time


% remoteAddr = '192.168.1.2'; % old computer
remoteAddr = '192.168.10.5'; % new computer
remotePort = 2020;
localPort = 9090;

off = 0;
on = 1;


%% Load TEPAdmin.dll which is a .Net Assembly
% The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32

% Currently the tested DLL is installed in the following path
asm = NET.addAssembly('C:\Windows\System32\TEPAdmin.dll');

import TaborElec.Proteus.CLI.*
import TaborElec.Proteus.CLI.Admin.*
%% Create Administrator
% Create instance of the CProteusAdmin class, and open it.
% Note that only a single CProteusAdmin instance can be open,
% and it must be kept open till the end of the session.

admin = CProteusAdmin(@OnLoggerEvent);

rc = admin.Open();
assert(rc == 0);

%% Use the administrator, and close it at the end
try
    slotIds = admin.GetSlotIds();
    numSlots = length(slotIds);
    assert(numSlots > 0);
    
    connect = 1;
    % If there are multiple slots, let the user select one ..
    sId = slotIds(2);
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
    
     sId = 8;
    
    % Connect to the selected instrument ..
    should_reset = false;
    inst = admin.OpenInstrument(sId);
    instId = inst.InstrId;
    
%     In other code...
%     % Connect to the selected instrument ..
%     should_reset = true;
%     inst = admin.OpenInstrument(sId, should_reset);
%     instId = inst.InstrId;
    
    % ---------------------------------------------------------------------
    % Setup instrument
    % ---------------------------------------------------------------------
    
    res = inst.SendScpi('*IDN?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nConnected to ''%s''\n', netStrToStr(res.RespStr));
    
    pause(0.01);
    
    res = inst.SendScpi('*CLS'); % clear
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi('*RST'); % reset
    assert(res.ErrCode == 0);
    
%     sampleRateDAC = 9e9;
%     sampleRateDAC_str = [':FREQ:RAST ' sprintf('%0.2e', sampleRateDAC)];
%     res = inst.SendScpi(sampleRateDAC_str); % set sample clock
%     assert(res.ErrCode == 0);
    
    fprintf('Reset complete\n');
    
    % ---------------------------------------------------------------------
    % ADC Config
    % ---------------------------------------------------------------------
    
    rc = inst.SetAdcDualChanMode(1); %set to dual frame granularity = 48
    assert(rc == 0);
    
    rc = inst.SetAdcSamplingRate(sampleRate);
    assert(rc == 0);
    
    rc = inst.SetAdcFullScaleMilliVolts(0,fullScaleMilliVolts);
    assert(rc == 0);
    
    fprintf('ADC Configured\n');
    
    %% Measurement Loop
    u2 = udp(remoteAddr, 'RemotePort', remotePort, 'LocalPort', localPort);
    try
        fopen(u2);
        connect=on;
        fprintf('Server connected and started\n');
    catch
        fclose(u2);
        connect=off;
        fprintf('Server not connected\n');
    end
    
    % Loop to repeatedly wait for messages and send replies
    % Break or Ctrl+C to get out of loop
    while ( connect==on )
        % Wait for message
        while(u2.BytesAvailable == 0)
            % If no bytes in u2 buffer, wait 10ms then check again
            pause(0.01);
        end
        %         cmdBytes = fread(u2);
        readBytes = fscanf(u2);
        dataBytes=1;counter=1;
        while ~isempty(readBytes)
            [tempBytes,readBytes]=strtok(readBytes,',');
            cmdBytes(counter)=str2num(tempBytes);
            counter=counter+1;
        end
        %dataBytes=strtok(dataBytes,',');
        cmdByte=cmdBytes(1);
        
        % 1 - Init
        % 2 - Aquire on Trig
        % 3 - Measure
        % 4 - Cleanup - prepare for next cycle
        % 5 - Disconnect Intrument
        % 6 - Program MW chirp waveform
        % 7 - Play MW chirp waveform
        
        %measure_time=str2double(dataBytes); %in sec
        switch cmdByte
            case 1 % Init
                
                rc = inst.SetAdcCaptureTrigSource(adcChanInd, trigSource);
                assert(rc == 0);
                
                rc = inst.SetAdcExternTrigPattern(0);
                assert(rc==0)
                
                rc = inst.SetAdcExternTrigGateModeEn(on);
                assert(rc==0);
                
                rc = inst.SetAdcExternTrigPolarity(1);
                assert(rc==0);
                
                rc = inst.SetAdcExternTrigThresh(0,0.3);
                assert(rc==0);
                
                % rc = inst.SetAdcExternTrigDelay(0.0000025);
                rc = inst.SetAdcExternTrigDelay(delay+0.000001);
                assert(rc==0);
                
                rc = inst.SetAdcAcquisitionEn(on,off);
                assert(rc == 0);
                
                fprintf('Instr setup complete and ready to aquire\n');
                
            case 2 % Aquire
 
                fprintf('Calculate and set data structures...\n');

%                  NB: In FID, there will only be one pulse-acquire cycle
                pw= cmdBytes(2)*1e-6; % pulse width in seconds, set by make_pw_nmr "pw" in microseconds
%                 pw= 40*1e-6;
                
%                 numberOfPulses_total= cmdBytes(3); % total number of
%                 pulses (usually around ~200000)
%                 Tmax= cmdBytes(4); % total acquisition time, in terms of number of 12-pulse cycles; irrelevant for FID
                
%                 length of acuisition in microseconds, set by make_pw_nmr "pw" in microseconds
                tacq=cmdBytes(5);
                Tmax=tacq;
                
                readLen = round2((tacq)*1e-6*sampleRate,96)-96; %constraint: has to be multiple of 96, add 4 of deadtime
                
                windows_power_of_2 = floor(log2(readLen/1e4));
                numberOfWindows = 2^(windows_power_of_2); %how many points do you want on the final FID curve?
                sampleWindow = round2(readLen/numberOfWindows,96)-96;
                
                numberOfPulses = 1;% floor(numberOfPulses_total/Tmax); %in one second  %will be 1 for FID
                loops = 1; % Tmax; % loops are now equivalent to 1-second
%                 'chunks' % irrelevant for FID --> loops=1
                
%                 Can I changed sampleRate to 6e9? -- check if it messes up the code elsewhere...
                
                offLen = 0;

                netArray = NET.createArray('System.UInt16', sampleWindow*numberOfWindows); %total array in one dimension -- all memory needed
                
                rc = inst.SetAdcAcquisitionEn(on,off);
                assert(rc == 0);
                
                rc = inst.SetAdcFramesLayout(numberOfPulses*loops, sampleWindow*numberOfWindows); %set memory of the AWG %for FID, only 1 total number of pulses, so only 1 frame needed
                assert(rc == 0);

                fprintf('Waiting... Listen for Shuttle\n');
                rc = inst.SetAdcCaptureEnable(on);
                assert(rc == 0);
                
            case 3
                
                status = inst.ReadAdcCaptureStatus();
                for i = 1 : 250
                    if status ~= 0
                        break;
                    end
                    pause(0.01);
                    status = inst.ReadAdcCaptureStatus();
                end
                
                rc = inst.SetAdcCaptureEnable(on);
                assert(rc == 0);

                rc = inst.ReadAdcCaptureStatus();

                readSize = uint64(sampleWindow*numberOfWindows);
                readOffset = uint64(offLen);

                chanIndex = 0;
                pulseAmp = [];
                
%                 Automatically calculate padded_len based on readLen,
%                 which is just based on tacq
%                 padded_len=2^20;
                power_of_2 = floor(log2(sampleWindow));
                padded_len = 2^(power_of_2);
                
                dF = sampleRate/padded_len; %set the discretization of freq in terms of sampleRate
                f = -sampleRate/2:dF:sampleRate/2-dF; %sampleRate sets the 'bandwidth'

                [v,b1]=min(abs(f-20e6)); %picks out the 20MHz component
                [v,b2]=min(abs(f+20e6));
                
%                 What does this line do:
                eff_read=100:sampleWindow*numberOfWindows-100;
                cyclesPoints = 50;
                
                fprintf('Shuttle complete\n')
                fprintf('Transfering aquired data to computer....\n')
                
%                 This is where the 1s-chunk 'loops' would have started:
                n=1;
                    fprintf('Start Read %d .... ', n);
                    firstIndex = ((n-1)*numberOfPulses)+1;
                    tic
%                     Perhaps we can take advantage of the buffering
%                     options? See the different ReadMultipleAdcFrames
%                     functions listed in the .h file open in Notepad
                    rc = inst.ReadMultipleAdcFrames(chanIndex, firstIndex, numberOfPulses, netArray); %this is where the device reads
                    assert(rc == 0);
                    samples = double(netArray); %get the data (1s chunk)
                    fprintf('Read %d Complete\n', n);
                    toc
                    
                    tic
                    %delete mem
%                     Do we need this section in here?
                    fprintf('Clear mem %d .... ', n);
                    rc =inst.WipeMultipleAdcFrames(chanIndex, ((n-1)*numberOfPulses)+1, numberOfPulses, 0);
                    assert(rc == 0);
                    fprintf('Clear mem %d Complete\n', n);
                    toc

                    tic 
                    fprintf('Starting iteration %d data processing....', n);
                    pulses = reshape(samples, [], numberOfWindows); % reshape samples into a more useful array (2 dimensions)

                    clear samples;
                    
%                     This is where the second for loop would have started:
                    for i = 1:numberOfWindows
                        pulse = pulses(:, i);
                        
                        if mod(i,96) == 0
                            figure(4);clf;
                            plot(pulse);
                            hold on;
                            yline(2048);
                        end
                        
                        if i == 1
                            figure(1);clf;
                            plot(pulse);
                            hold on;
                            yline(2048);
                        end
                        
                        if i == 2
                            figure(2);clf;
                            plot(pulse);
                            hold on;
                            yline(2048);
                        end
                        
                        if i == 3
                            figure(3);clf;
                            plot(pulse);
                            hold on;
                            yline(2048);
                        end
                        
%                         What were these two commented-out lines for:
%                         pulse = pulse(1024:readLen);
%                         readLen=length(pulse);
                        
                        pulseMean = mean(pulse);
                        pulseDC = pulse - pulseMean; % remove DC
                        X = fftshift(fft(pulseDC,padded_len)); % perform FFT and zero-pad to the padded_len
                        linFFT = (abs(X)/(sampleWindow*numberOfWindows));
                        
%                         Would  location be useful information?
%                         [amp, loc] = max(linFFT);

                        amp=(linFFT(b1) + linFFT(b2));
            %             linFFT_vec(:,i)=linFFT;
                        pulseAmp(i+(1*(n-1))) = amp;
                    end
                    clear pulses;
                    fprintf('Data processing iteration %d complete!\n', n);
                    toc
%                 end
                
                ivec=1:numberOfWindows*loops;
%                 time_cycle=pw+(4+2+2+tacq)*1e-6;
%                 time_axis=time_cycle.*ivec;
                time_window = sampleWindow/sampleRate;
                time_axis = time_window.*ivec;
                
%                 %drop first point
%                 time_axis(1)=[];pulseAmp(1)=[];
                try
                figure(5);clf;
%                plot(time_axis,(pulseAmp),'b-');hold on;
                plot(time_axis,(smooth(pulseAmp,10)),'r-');
                set(gca,'ylim',[0 max(pulseAmp)*1.1]);  
                catch
                    disp('Plot error occured');
                end
                
                %fn=dataBytes; %filename
                a = datestr(now,'yyyy-mm-dd-HHMMSS');
                fn = sprintf([a,'_Proteus']);
                % Save data
                fprintf('Writing data to Z:.....\n');
                save(['Z:\' fn],'pulseAmp','time_axis');
                
            case 4 % Cleanup, save and prepare for next experiment
                
                rc = inst.SetAdcCaptureEnable(off);
                assert(rc == 0);
                % Free the memory space that was allocated for ADC capture
                rc = inst.FreeAdcReservedSpace();
                assert(rc == 0);
                
                %                 clear pulseAmp time_axis;
                fprintf('Complete ready for next aquisition!\n');
                
            case 5 % Shutdown disconnect instr
                
                connect = 0;
            
            case 6 % Program MW chirp waveform
                
                % ---------------------------------------------------------------------
                % Download chirp waveform of 1024 points to segment 1 of channel 1
                % ---------------------------------------------------------------------
                
                awg_center_freq = cmdBytes(2);
                awg_bw_freq = cmdBytes(3);
                awg_amp = cmdBytes(4);
                sweep_freq = cmdBytes(5);
                sweep_sigma = cmdBytes(6);
                symm = cmdBytes(7);
                srs_freq = cmdBytes(8); % should be 3.4e9
                srs_amp = cmdBytes(9);
                
                res = inst.SendScpi(['INST:CHAN ' num2str(dacChanInd)]); % select channel 2
                assert(res.ErrCode == 0);
                
                bits = 8;
                
                rampTime = 1/sweep_freq;
                fCenter = awg_center_freq - srs_freq;
                fStart = fCenter - 0.5*awg_bw_freq;
                fStop = fCenter + 0.5*awg_bw_freq;
                dt = 1/sampleRateDAC;
                dacSignal = makeChirp(sampleRateDAC, rampTime, dt, fStart, fStop, bits);   
                fprintf('waveform length - ');
                fprintf(num2str(length(dacSignal)));
                fprintf('\n') ;
                
                % Define segment 1 
                res = inst.SendScpi(strcat(':TRAC:DEF 1,',num2str(length(dacSignal)))); % define memory location 1 of length dacSignal
                assert(res.ErrCode == 0);
                
                % select segmen 1 as the the programmable segment
                res = inst.SendScpi(':TRAC:SEL 1');
                assert(res.ErrCode == 0); 
                
                % Download the binary data to segment 1
                res = inst.WriteBinaryData(':TRAC:DATA 0,#', dacSignal);
                assert(res.ErrCode == 0);
                
                srs_freq_str = [':SOUR:CFR ' sprintf('%0.2e', srs_freq)];
                res = inst.SendScpi(srs_freq_str);
                assert(res.ErrCode == 0);
                
                res = inst.SendScpi(':SOUR:MODE IQM'); % IQ MODULATOR --
%                 THINK OF A MIXER, BUT DIGITAL
                assert(res.ErrCode == 0);
                
                try
                sampleRateDAC_str = [':FREQ:RAST ' sprintf('%0.2e', sampleRateDAC)];
                res = inst.SendScpi(sampleRateDAC_str); % set sample clock
                assert(res.ErrCode == 0);
                catch
                sampleRateDAC_str = [':FREQ:RAST ' sprintf('%0.2e', sampleRateDAC)];
                res = inst.SendScpi(sampleRateDAC_str); % set sample clock
                assert(res.ErrCode == 0);
                end
                
                segLen = 40960;
                segLenDC = 50048;
                cycles = 800;
                
                dacSignalMod = sine(cycles, 0, segLen, bits); % sclk, cycles, phase, segLen, amp, bits, onTime(%)
                
                dacSignalDC = [];
                for i = 1 : segLenDC
                  dacSignalDC(i) = 127;
                end
                
                dacSignal = [dacSignalMod dacSignalDC];
                
            case 7 % Play MW chirp waveform
                
                % ---------------------------------------------------------------------
                % Play segment 1 in channel 1
                % ---------------------------------------------------------------------
                
                res = inst.SendScpi(['INST:CHAN ' num2str(dacChanInd)]); % select channel 2
                assert(res.ErrCode == 0);
                
                res = inst.SendScpi(':SOUR:FUNC:MODE:SEG 1');
                assert(res.ErrCode == 0);
                
                amp_str = [':SOUR:VOLT ' sprintf('%0.2f', awg_amp)];
                res = inst.SendScpi(amp_str);
                assert(res.ErrCode == 0);
                
                res = inst.SendScpi(':OUTP ON');
                assert(res.ErrCode == 0);
                
                fprintf('Waveform generated and playing\n');
                
            case 8
                
                % Disable MW chirp output
                res = inst.SendScpi(':OUTP OFF');
                assert(res.ErrCode == 0);
                
                fprintf('MW Chirp Waveform stopped playing (on purpose)');
                
        end % end switch
        
    end % end while
    
    res = inst.SendScpi(':SYST:ERR?');
    fprintf(1, '\nEnd - server stopped!! \nInstrument Error Status = %s\n', netStrToStr(res.RespStr));
    
    % It is recommended to disconnect from instrumet at the end
    rc = admin.CloseInstrument(instId);
    
    % Close the administrator at the end ...
    admin.Close();
catch ME
    admin.Close();
    rethrow(ME)
end

fclose(u2);

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

%% Signal Processing Functions

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

function sine = addSinePulse(segment, starttime, dt, pulseDuration, freq, phase, bits)            
    rawSignal = segment;
    npoints = length(segment);
    starttimeIdx = roundCheck(starttime, dt)+1;
    durationPts = roundCheck(pulseDuration, dt);

    if starttimeIdx > 0 && starttimeIdx + durationPts <= npoints
        rawSignal(starttimeIdx:starttimeIdx+durationPts) = sin(2*pi*freq*dt*(0:durationPts) + 2*pi*freq*starttime + phase);
        sine = ampScale(bits, rawSignal);
    else
        error('Pulse out of bounds');
    end
end

function dacWav = makeChirp(sampleRateDAC, rampTime, dt, fStart, fStop, bits)            

    t = 0:1/sampleRateDAC:rampTime;
    dacWave = chirp(t,fStart,rampTime,fStop);
    seglenTrunk = (floor(length(dacWave)/ 64))*64;
    dacWave = dacWave(1:seglenTrunk);
    dacWav = ampScale(bits, dacWave);

end

% Scale to FSD
function dacSignal = ampScale(bits, rawSignal)
 
  maxSig = max(rawSignal);
  verticalScale = ((2^bits)/2)-1;

  vertScaled = (rawSignal / maxSig) * verticalScale;
  dacSignal = uint8(vertScaled + verticalScale);
  %plot(dacSignal);

%   if bits > 8
%       dacSignal16 = [];
%       sigLen = length(dacSignal);
%       k=1;
%       for j = 1:2:sigLen*2;
%         dacSignal16(j) = bitand(dacSignal(k), 255);
%         dacSignal16(j+1) = bitshift(dacSignal(k),-8);
%         k = k + 1;
%       end
%       dacSignal = dacSignal16;
%   end
end

function rNb = roundCheck(nb, units)
    global debug
    rNb = round(nb/units);
    if ~isempty(debug) && debug > 1 && abs(rNb*units - nb) > units/100
       disp([inputname(1) ' = ' num2str(nb/units) ' rounded to ' num2str(rNb)]);
    end
end
