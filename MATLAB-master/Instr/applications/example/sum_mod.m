%% make signal

segLen = 1024;
cycles = 1;
cycles2 =100;
phase = 0;
sclk=1e9;

tic
time = -(segLen-1)/2:(segLen-1)/2;
omega = 2 * pi() * cycles;
omega2 = 2 * pi() * cycles2;

rawSignal_1 = sin(((omega*time)/segLen)+(phase)); 
rawSignal_2 = sawtooth(((omega2*time)/segLen)+(phase),1/2);
rawSignal = rawSignal_1+rawSignal_2;
rawSignal = rawSignal_1;

dacBits = 8;
maxSig = max(rawSignal);
verticalScale = ((2^dacBits)/2)-1;
vertScaled = (rawSignal / maxSig) * verticalScale;
dacSignal = uint8(vertScaled + verticalScale);
dacSignal = typecast(dacSignal, 'uint8');
toc
plot(dacSignal);

%% Send Signal to instrument

asm = NET.addAssembly('C:\Users\Simon\Documents\Tabor Electronics\ProteusAwg-10204-Customer\TEPAdmin.dll');
import TaborElec.Proteus.CLI.*
import TaborElec.Proteus.CLI.Admin.*

admin = CProteusAdmin(@OnLoggerEvent);
rc = admin.Open();
assert(rc == 0);

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
    
    % ---------------------------------------------------------------------
    % Send SCPI commands
    % ---------------------------------------------------------------------
    
    res = inst.SendScpi('*IDN?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nConnected to ''%s''\n', netStrToStr(res.RespStr));
    
    res = inst.SendScpi('*CLS');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi('*RST');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi(':FREQ:RAST 1E9');
    assert(res.ErrCode == 0);    
    
    fprintf('Reset complete\n');

    % ---------------------------------------------------------------------
    % Download sine waveform of 1024 points to segment 1 of channel 1
    % ---------------------------------------------------------------------

    res = inst.SendScpi('INST:CHAN 1');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi('TRAC:DEL 1');
    assert(res.ErrCode == 0);
    
    sampleRateDAC = 1E9;
    segLen = 2048;
    bits = 8;
    amplitude = 1;
    cycles = 1;

    % Define segment 1 
    res = inst.SendScpi(':TRAC:DEF 1,1024');
    assert(res.ErrCode == 0);
    
    % select segmen 1 as the the programmable segment
    res = inst.SendScpi(':TRAC:SEL 1');
    assert(res.ErrCode == 0); 
    tic
    % Download the binary data to segment 1
    res = inst.WriteBinaryData(':TRAC:DATA 0,#', dacSignal);
    assert(res.ErrCode == 0);

    toc

    % ---------------------------------------------------------------------
    % Play segment 1 in channel 1
    % ---------------------------------------------------------------------

    res = inst.SendScpi('INST:CHAN 1');
    assert(res.ErrCode == 0);  

    res = inst.SendScpi(':SOUR:FUNC:MODE:SEG 1');
    assert(res.ErrCode == 0);

%     res = inst.SendScpi(':SOUR:VOLT 0.3');
%     assert(res.ErrCode == 0);

    res = inst.SendScpi(':OUTP ON');
    assert(res.ErrCode == 0);

    fprintf('Waveform generated and playing\n');
    
    %end
    
    res = inst.SendScpi(':SYST:ERR?');
    fprintf(1, '\nEnd of Example - %s\n', netStrToStr(res.RespStr));
    %close all % Close all figures
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