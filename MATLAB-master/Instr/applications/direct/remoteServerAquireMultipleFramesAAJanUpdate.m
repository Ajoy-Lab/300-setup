%% Clear everything
clear;
close;

%% Set defaults Vars
sampleRate = 1000e6;
adcDualChanMode = 1;
fullScaleMilliVolts =1000;
trigSource = 1; % trigger for ADC channel-1
adcChanInd = 0; % ADC Channel 1
measurementTimeSeconds = 7; %Integer

% start of the pulse from trigger is 16.1us (delay)
trigToPulseStart = 0.0000161;
measStart = 0;

remoteAddr = '192.168.1.2';
remotePort = 2020;
localPort = 9090;


%% Load TEPAdmin.dll which is a .Net Assembly
% The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32 

% Currently the tested DLL is installed in the following path
asm = NET.addAssembly('C:\Users\Simon\Documents\Tabor Electronics\ProteusAwg-10204-Customer\TEPAdmin.dll');

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
    should_reset = false;
    inst = admin.OpenInstrument(sId);
    instId = inst.InstrId;

    % ---------------------------------------------------------------------
    % Setup instrument
    % ---------------------------------------------------------------------

    res = inst.SendScpi('*IDN?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nConnected to ''%s''\n', netStrToStr(res.RespStr));
   
    pause(0.01);
    
    
    res = inst.SendScpi('*CLS');
    assert(res.ErrCode == 0);

    res = inst.SendScpi('*RST');
    assert(res.ErrCode == 0);

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
    fprintf('ADC Configured\n');
    u2 = udp(remoteAddr, 'RemotePort', remotePort, 'LocalPort', localPort);
    try
        fopen(u2);
        connect=1 
        fprintf('Server Started\n');
    catch
        fclose(u2);
        connect=0;
        fprintf('Server Not Started!\n');
    end
    
    % Loop to repeatedly wait for messages and send replies
    % Break or Ctrl+C to get out of loop
    while ( connect==1 )
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
        %cmdByte
        % 1 - Init
        % 2 - Aquire on Trig
        % 3 - Measure
        % 4 - Cleanup - prepare for next cycle
        % 5 - Disconnect Intrument
        
        %measure_time=str2double(dataBytes); %in sec
        
        switch cmdByte
            case 1 % Setup
%                 numberOfPulses = 2717*4; % 1 second of pulses
%                 numberOfPulses = floor(11111*12/12); % 1 second of pulses
%                 loops = measurementTimeSeconds;
% %                 readLen = 1056*34; %Granuality 48 2CH, 96 1CH
% %                 readLen = round2(36e-6/1e-9,96)-96; %Granuality 48 2CH, 96 1CH
%                 offLen = 0;

%                 rc = inst.SetAdcFullScaleMilliVolts(0,fullScaleMilliVolts);
%                 assert(rc == 0);

                fprintf('Setting triggers\n');

                rc = inst.SetAdcCaptureTrigSource(adcChanInd, trigSource);
                assert(rc == 0);
                
                % trigSource
                % - 0: internal-trigger (generated by software).
                % - 1: external-trigger.
                % - 2: self-trigger from ADC channel-1 (what we originally used)
                % - 3: self-trigger from ADC channel-2.
                
                % start of the pulse from trigger is 16.1us (delay)
%                 rc = inst.SetAdcExternTrigParams(0,1,0,0,0,trigToPulseStart+measStart);
%                 assert(rc == 0);
                
                rc = inst.SetAdcExternTrigParams(0,1,0,1,0,17e-6);
                assert(rc == 0);

                rc = inst.SetAdcExternTrigThresh(0,0.5); % voltage threshold
                assert(rc==0);
                
                
%                 rc = inst.SetAdcFramesLayout(numberOfPulses*loops, readLen);
%                 assert(rc == 0);       
%                 
%                 rc = inst.SetAdcMultiFrameModeEn(1);
%                 assert(rc == 0);

                fprintf('Instr Setup Complete\n');
                
            case 2 % Aquire
 
                fprintf('Waiting for meas data...\n');

                pw= cmdBytes(2)*1e-6;
                numberOfPulses_total= cmdBytes(3);
                Tmax=cmdBytes(4);
                tacq=cmdBytes(5);
                numberOfPulses= floor(numberOfPulses_total/Tmax); %in one second
                loops=Tmax;
                
                readLen = round2((tacq+2)*1e-6/1e-9,96)-96; %constaint: has to be multiple of 96, add 4 of deadtime
                offLen = 0;
                                
                fprintf('Measurement started\n');
                rc = inst.SetAdcAcquisitionEn(1,0);
                assert(rc == 0);
                
                rc = inst.SetAdcFramesLayout(numberOfPulses*loops, readLen);
                assert(rc == 0);

                
%                 rc = inst.SetAdcMultiFrameModeEn(1);
%                 assert(rc == 0);
                
                 rc = inst.SetAdcCaptureEnable(1);
                assert(rc == 0);

            case 3
                
                status = inst.ReadAdcCaptureStatus();
                for i = 1 : 250
                    if status ~= 0
                        break;
                    end
                    pause(0.01);
                    status = inst.ReadAdcCaptureStatus()
                end

                %readSize = uint64(100000);
                %readOffset = uint64(100000);

                readSize = uint64(readLen);
                readOffset = uint64(offLen);


                chanIndex = 0;
                %netArray = NET.createArray('System.UInt16', 100000);

                netArray = NET.createArray('System.UInt16', readLen*numberOfPulses);

                %pulseAmp = [numberOfPulses*seconds];
                pulseAmp = [];
                padded_len=2^15;
                dF = sampleRate/padded_len;
                f = -sampleRate/2:dF:sampleRate/2-dF;

                [v,b1]=min(abs(f-20e6));
                [v,b2]=min(abs(f+20e6));

                eff_read=100:readLen-100;
                cyclesPoints = 50;

                for n = 1:loops
                    n
                    fprintf('Start Read\n');
                    tic
                    rc = inst.ReadMultipleAdcFrames(chanIndex, ((n-1)*numberOfPulses)+1, numberOfPulses, netArray);
                    assert(rc == 0);
                    samples = double(netArray);
                    fprintf('Read Complete\n');
                    toc
                    
                    tic
                    %delete mem
                    rc =inst.WipeMultipleAdcFrames(chanIndex, ((n-1)*numberOfPulses)+1, numberOfPulses, 0);
                     assert(rc == 0);
                     toc

                    tic 
                    fprintf('Start Processing\n');
                    pulses = reshape(samples, [], numberOfPulses);
                    clear samples;
                    for i = 1:numberOfPulses
                        pulse = pulses(:, i);
                        pulseMean = mean(pulse);
                        pulseDC = pulse - pulseMean; % remove DC
                         X = fftshift(fft(pulseDC,padded_len));
                        linFFT = (abs(X)/readLen);
%                         [amp, loc] = max(linFFT);
                        amp=(linFFT(b1) + linFFT(b2));
            %             linFFT_vec(:,i)=linFFT;
                        pulseAmp(i+(numberOfPulses*(n-1))) = amp;
                        if n==1 && i==1
                            firstPulse = pulse;
                        end
                    end
                    clear pulses;
                    fprintf('End Processing\n');
                    toc
                end

                ivec=1:numberOfPulses*loops;
                time_cycle=pw+(4+2+2+tacq)*1e-6;
                time_axis=time_cycle.*ivec;
                %drop first point
                time_axis(1)=[];pulseAmp(1)=[];
                try
                figure(2);clf;
                plot(pulse);
                figure(3);clf;
                plot(firstPulse);
                
                figure(1);clf;
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
            
            case 4 % Cleanup and save

                rc = inst.SetAdcCaptureEnable(0);
                assert(rc == 0);
                % Free the memory space that was allocated for ADC capture
                rc = inst.FreeAdcReservedSpace();
                assert(rc == 0);
                
%                 clear pulseAmp time_axis;
                fprintf('Complete ready for next aquisition!\n');
            
             case 5 % Shutdown disconnect instr
                 
                 connect = 0;
                 
        end % end switch
        
    end % end while
    
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
