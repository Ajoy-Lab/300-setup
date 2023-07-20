%% Proteus P2584M Example
% Connect to P2584M instrument

%% Load TEPAdmin.dll which is a .Net Assembly
% The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32 
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
    catch ME
       rethrow(ME) 
    end
end