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
                %                 rc = inst.SetAdcExternTrigDelay(delay+0.000001+33.8e-6+33.8e-6 - 12e-6);
                rc = inst.SetAdcExternTrigDelay(4.1e-6);
                assert(rc==0);
                
                rc = inst.SetAdcAcquisitionEn(on,off);
                assert(rc == 0);
                
                fprintf('Instr setup complete and ready to aquire\n');
                
            case 2 % Aquire on trig
                
                fprintf('Calculate and set data structures...\n');
                
                pw = cmdBytes(2)*1e-6
                %                 pw=110*1e-6;
                %                 pw= 30.0*1e-6
                %pw=32.5e-6
                %pw=42.0e-6;
                %
                
                
                pw3=cmdBytes(3)*1e-6
                %                 nc=135135;
                %                 nc2=117187;
                nc=cmdBytes(4);
                nc2=cmdBytes(5);
                nc3=cmdBytes(9);
                %nc2=11160*12;
                numberOfPulses_total= nc; %number of pulses for the first rigol time
                numberOfPulses_total2= nc2; %number of pulses for the second rigol time
                numberOfPulses_total3= nc3; %number of pulses for the remaining time in the experiment
                
                
                %Tmax=cmdBytes(4);
                %                 Tmax=15;
                %                 Tmax2=15;
                T1=cmdBytes(6); %first rigol time
                T2=cmdBytes(7); %second rigol time
                
                
                tacq=cmdBytes(8);
                
                rep=cmdBytes(10); %number of times rigol1 and 2 are repeated
                Ttot=cmdBytes(11); %total experiment time
                Trem=Ttot-rep*(T1+T2);
                
                numberOfPulses = floor(numberOfPulses_total/T1);
%                 numberOfPulses2 =floor(numberOfPulses_total2/T2); %in 1 second %will be 1 for FID
%                 numberOfPulses3 = floor(numberOfPulses_total3/Trem);
                %loops=T1;
                %loops2=T2;
                %remloops=Ttot-rep*(T1+T2);
                %looptot=Ttot;
                
                %trying to do it without loops determied by time
                totpulses=rep*(nc+nc2)+nc3;
                noofframes=10000;
                loops=floor(totpulses/noofframes);
                remframes=mod(totpulses,noofframes);
                
                readLen = round2((tacq+2)*1e-6/1e-9,96)-96;
                %readLen = round2((tacq)*1e-6/1e-9,96)-96;
                %                 readLen = 33888; % round2((tacq+2)*1e-6/1e-9,96)-96 %constraint: has to be multiple of 96, add 4 of deadtime %number of points in a frame
                %readLen = 16896; % for tacq=16
                %                 readLen = 65856; % for tacq=64
                %                 readLen = 129888; % for tacq=128
                
                %readLen = 66912; % actual tacq=64
                %                 readLen = 131040;
                %                 readLen = 129600;
                %                 readLen = 126912; % actual tacq=128
                %                 readLen = 97920; % for tacq=96
                
                offLen = 0;
                
     
                netArray = NET.createArray('System.UInt16', readLen*noofframes); %total array -- all memory needed
                netArrayRem = NET.createArray('System.UInt16', readLen*remframes); %total array -- all memory needed
                
                
                rc = inst.SetAdcAcquisitionEn(on,off);
                assert(rc == 0);
                
%                 rc = inst.SetAdcFramesLayout((numberOfPulses*loops+numberOfPulses2*loops2)*rep+numberOfPulses3, readLen); %set memory of the AWG
%                 assert(rc == 0);
                
                rc = inst.SetAdcFramesLayout(totpulses, readLen); %set memory of the AWG
                assert(rc == 0);
                
                fprintf('Waiting... Listen for Shuttle\n');
                rc = inst.SetAdcCaptureEnable(on);
                assert(rc == 0);
                
                
            case 3 % Measure
                
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
                
                
                chanIndex = 0;
                pulseAmp = [];
                relPhase = [];
                
                power_of_2 = floor(log2(readLen)); % this is normally 15 for 32us captures
                padded_len= 2^(power_of_2) ;%2^15;
                dF = sampleRate/padded_len; %set the discretization of freq in terms of sampleRate
                f = -sampleRate/2:dF:sampleRate/2-dF; %sampleRate sets the 'bandwidth'
                
                [v,b1]=min(abs(f-20.00706e6)); %picks out the 20MHz component
                [v,b2]=min(abs(f+20.00706e6));
                
                eff_read=100:readLen-100;
                cyclesPoints = 50;
                fprintf('Shuttle complete\n')
                fprintf('Transfering aquired data to computer....\n')
                %fprintf('Total number of Loops: %s\n', num2str(loops));
                    for n = 1:loops
                        fprintf('Start Read %d/%d .... ', n,loops);
                        firstIndex=(n-1)*noofframes+1;
                        tic
                        rc = inst.ReadMultipleAdcFrames(chanIndex, firstIndex, noofframes, netArray); %this is where the device reads
                        assert(rc == 0);
                        samples = double(netArray); %get the data (1s chunk)
                        fprintf('Read %d/%d Complete\n', n,loops);
                        toc
                        
                        tic
                        %delete mem
                        fprintf('Clear mem %d/%d .... ', n,loops);
                        rc =inst.WipeMultipleAdcFrames(chanIndex, firstIndex, noofframes, 0);
                        assert(rc == 0);
                        fprintf('Clear mem %d/%d Complete\n', n,loops);
                        toc
                        
                        tic
                        fprintf('Starting iteration %d/%d data processing....', n,loops);
                        pulses = reshape(samples, [], noofframes); % reshape samples into a more useful array (2 dimensions)
                        
                        
                        clear samples;
                        for i = 1:noofframes
                            pulse = pulses(:, i);
                         
                            pulseMean = mean(pulse);
                            pulseDC = pulse - pulseMean; % remove DC
                            X = fftshift(fft(pulseDC,padded_len)); % perform FFT and zero-pad to the padded_len
                            linFFT = (abs(X)/readLen);
                            %                         [amp, loc] = max(linFFT);
                            amp=(linFFT(b1) + linFFT(b2));
                            phase1=angle(X(b1));
                            phase2=angle(X(b2));
                            phase=phase1-phase2;
                            %             linFFT_vec(:,i)=linFFT;
                            pulseAmp(i+(noofframes*(n-1))) = amp;
                            relPhase(i+(noofframes*(n-1))) = phase1;
                        end
                 
                        clear pulses;
                        fprintf('Data processing iteration %d/%d complete!\n', n,loops);
                        toc
                    end
                    
                    
                      fprintf('Start Read Remaining .... ');
                        firstIndex=loops*noofframes+1;
                        tic
                        rc = inst.ReadMultipleAdcFrames(chanIndex, firstIndex, remframes, netArrayRem); %this is where the device reads
                        assert(rc == 0);
                        samples = double(netArrayRem); %get the data (1s chunk)
                        fprintf('Read Remaining Complete\n');
                        toc
                        
                        tic
                        %delete mem
                        fprintf('Clear mem remaining .... ');
                        rc =inst.WipeMultipleAdcFrames(chanIndex, firstIndex, remframes, 0);
                        assert(rc == 0);
                        fprintf('Clear mem Remaining Complete\n', n);
                        toc
                        
                        tic
                        fprintf('Starting remaining data processing....');
                        pulses = reshape(samples, [], remframes); % reshape samples into a more useful array (2 dimensions)
                        
                        
                        clear samples;
                        for i = 1:remframes
                            pulse = pulses(:, i);
                         
                            pulseMean = mean(pulse);
                            pulseDC = pulse - pulseMean; % remove DC
                            X = fftshift(fft(pulseDC,padded_len)); % perform FFT and zero-pad to the padded_len
                            linFFT = (abs(X)/readLen);
                            %                         [amp, loc] = max(linFFT);
                            amp=(linFFT(b1) + linFFT(b2));
                            phase1=angle(X(b1));
                            phase2=angle(X(b2));
                            phase=phase1-phase2;
                            %             linFFT_vec(:,i)=linFFT;
                            pulseAmp(i+(noofframes*loops)) = amp;
                            relPhase(i+(noofframes*loops)) = phase1;
                        end
                 
                        clear pulses;
                        fprintf('Data processing the remaining frames complete!\n');
                        toc
                    
                  
                
                %ivec=1:numberOfPulses*loops;
                ivec=1:length(pulseAmp);
                delay2 = 0.000003; % dead time the unknown one, this is actually rof3 -Ozgur
                
                %time_cycle=pw+96+(tacq+2+4+2+delay2)*1e-6;
                time_cycle=pw+delay2+(tacq+2+4+2)*1e-6 %for pentagon
                %time_cycle=pw+extraDelay+(4+2+2+tacq+17)*1e-6;
                time_cycle2=pw3+delay2+(tacq+2+4+2)*1e-6 %for square
                
                %time_axis=[time_cycle:time_cycle:nc*time_cycle nc*time_cycle+time_cycle2:time_cycle2:nc*time_cycle+time_cycle2+nc2*time_cycle2];
                time_axis=[];
                for ind=1:rep
                    time_axis=[time_axis, (ind-1)*(T1+T2)+[time_cycle*(1:nc) time_cycle*nc+time_cycle2*(1:nc2)]];
                end
                time_axis= [time_axis rep*(T1+T2)+time_cycle2*(1:nc3)];
                %drop first point
                time_axis(1)=[];pulseAmp(1)=[];relPhase(1)=[];
                try
                    figure(1);clf;
                    %                plot(time_axis,(pulseAmp),'b-');hold on;
                    %time_axis=time_axis(1:length(pulseAmp));
                    plot(time_axis,pulseAmp,'r-');
                    %plot(1:1:length(pulseAmp),pulseAmp,'r-');
                    set(gca,'ylim',[0 max(pulseAmp)*1.1]);
                    
                    figure(12);clf;
                    plot(time_axis,relPhase);
                    %plot(1:1:length(relPhase),relPhase);
                    
                    %                 figure(3);clf;
                    %                 plot(pulse);
                    %                 hold on;
                    %                 yline(2048);
                catch
                    disp('Plot error occured');
                end
                
                %fn=dataBytes; %filename
                a = datestr(now,'yyyy-mm-dd-HHMMSS');
                fn = sprintf([a,'_Proteus']);
                % Save data
                fprintf('Writing data to Z:.....\n');
                save(['Z:\' fn],'pulseAmp','time_axis','relPhase');
                fprintf('Save complete\n');
                
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
                
                % Define segment 2
                res = inst.SendScpi(strcat(':TRAC:DEF 2,',num2str(length(dacSignal)))); % define memory location 1 of length dacSignal
                assert(res.ErrCode == 0);
                
                % select segmen 2 as the the programmable segment
                res = inst.SendScpi(':TRAC:SEL 2');
                assert(res.ErrCode == 0);
                
                % Download the binary data to segment 2
                res = inst.WriteBinaryData(':TRAC:DATA 0,#', fliplr(dacSignal));
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
                
                %                 dacSignalMod = sine(cycles, 0, segLen, bits); % sclk, cycles, phase, segLen, amp, bits, onTime(%)
                %
                %                 dacSignalDC = [];
                %                 for i = 1 : segLenDC
                %                   dacSignalDC(i) = 127;
                %                 end
                %
                %                 dacSignal = [dacSignalMod dacSignalDC];
                
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
                
                fprintf('MW Chirp Waveform stopped playing (on purpose)\n');
                
            case 9
                
                % ---------------------------------------------------------------------
                % Play segment 2 in channel 1
                % ---------------------------------------------------------------------
                
                %                 res = inst.SendScpi(['INST:CHAN ' num2str(dacChanInd)]); % select channel 2
                %                 assert(res.ErrCode == 0);
                
                res = inst.SendScpi(':SOUR:FUNC:MODE:SEG 2');
                assert(res.ErrCode == 0);
                
                %                 amp_str = [':SOUR:VOLT ' sprintf('%0.2f', awg_amp)];
                %                 res = inst.SendScpi(amp_str);
                %                 assert(res.ErrCode == 0);
                
                %                 res = inst.SendScpi(':OUTP ON');
                %                 assert(res.ErrCode == 0);
                
                fprintf('Waveform in Segment 2 playing\n');
                
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
