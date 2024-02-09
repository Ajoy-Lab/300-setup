import os
import sys
import inspect
import clr  
import numpy as np
import time

datapath = os.path.dirname(sys.argv[0])
maxScpiResponse = 65535

#########################  - Set driver path - #########################
########################################################################
winpath = R'C:\Windows\System32'
winpath = os.path.join(winpath, 'TEPAdmin.dll')
########################################################################

if (datapath):
    datapath = datapath + "\\"
print(datapath)


def getSlotId(admin):
    try:
        if admin.IsOpen():
            rc = admin.Close()
            Validate(rc, __name__, inspect.currentframe().f_back.f_lineno)

        rc = admin.Open()
        Validate(rc, __name__, inspect.currentframe().f_back.f_lineno)

        slotIds = admin.GetSlotIds()
        n = 0
        for i in range(0, slotIds.Length, 1):
            slotId = slotIds[i]
            slotInfo = admin.GetSlotInfo(slotId)
            if slotInfo:
                if not slotInfo.IsDummySlot:
                    n = n + 1
                    print(("{0}. Slot-ID {1} [chassis {2}, slot {3}], "
                           "IsDummy={4}, IsInUse={5}, IDN=\'{6}\'").
                          format(i + 1,
                                 slotId,
                                 slotInfo.ChassisIndex,
                                 slotInfo.SlotNumber,
                                 'Yes' if slotInfo.IsDummySlot != 0 else 'No',
                                 'Yes' if slotInfo.IsSlotInUse != 0 else 'No',
                                 slotInfo.GetIdnStr()))
                else:
                    dummy = 1
                    # print("{0}. Slot-ID {1} - Failed to acquire Slot Info".
                    # .format(i + 1,slotId))

        if n == 1:
            sel = slotIds[0]
        else:
            sel = input("Please select slot-Id:")
        slotId = np.uint32(sel)
    except Exception as e:  # pylint: disable=broad-except
        print(e)
    return slotId


def OnLoggerEvent(sender, e):
    del sender
    print(e.Message.Trim())
    if (e.Level <= LogLevel.Warning):  # @UndefinedVariable
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print(e.Message.Trim())
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~")


def Validate(rc, condExpr, funcName="", lineNumber=0):
    _ = condExpr

    # cond = (rc == 0)

    if rc != 0:
        errMsg = "Assertion \"{0}\" Failed at line {1} of {2}."
        errMsg = errMsg.format(rc, lineNumber, funcName)
        raise Exception(errMsg)


def loadDLL():
  
    clr.AddReference(winpath)

    # pylint: disable=import-error
    from TaborElec.Proteus.CLI.Admin import CProteusAdmin  # @UnresolvedImport

    from TaborElec.Proteus.CLI.Admin import IProteusInstrument  # @UnusedImport @UnresolvedImport @IgnorePep8

    return CProteusAdmin(OnLoggerEvent)


def SendBinScpi(inst, prefix, path, query_err=False):

    err_code, resp_str = -1, ''
    try:
        print(prefix)
        inBinDat = bytearray(path, "utf8")
        inBinDatSz = np.uint64(len(inBinDat))
        DummyParam = np.uint64(0)
        res = inst.WriteBinaryData(prefix, inBinDat, inBinDatSz, DummyParam)

        err_code = int(res.ErrCode)
        resp_str = str(res.RespStr).strip()

        if 0 != err_code:
            print("Error {0} ({1})".format(err_code, resp_str))
        elif len(resp_str) > 0:
            print("{0}".format(resp_str))
        if query_err:
            err = inst.SendScpi(':SYST:ERR?')
            if not err.RespStr.startswith('0'):
                print(err.RespStr)
                err = inst.SendScpi('*CLS')

    except Exception as e:  # pylint: disable=broad-except
        print(e)

    return err_code, resp_str


def SendScpi(inst, line, query_err=False, print_line=True):

    err_code, resp_str = -1, ''
    try:
        if print_line:
            print(line)
        line = line + "\n"
        res = inst.SendScpi(str(line))
        err_code = int(res.ErrCode)
        resp_str = str(res.RespStr).strip()

        if 0 != err_code:
            if not print_line:
                print(line.strip())
            print("Error {0} - ({1})".format(err_code, resp_str))
        elif len(resp_str) > 0 and print_line:
            print("{0}".format(resp_str))
        if query_err:
            err = inst.SendScpi(':SYST:ERR?')
            if not err.RespStr.startswith('0'):
                print(err.RespStr)
                err = inst.SendScpi('*CLS')
    except Exception as e:  # pylint: disable=broad-except
        print(e)

    return err_code, resp_str
    
def scaleWaveform(rawSignal, model):

    if model == "P9082M":  # 9GS/s
        bits = 8
        wpt_type = np.uint8
    else:  # 2GS/s or 1.25GS/s models
        bits = 16
        wpt_type = np.uint16

    maxSig = max(rawSignal)
    verticalScale = ((pow(2, bits))/2)-1
    vertScaled = (rawSignal/maxSig) * verticalScale
    dacSignal = (vertScaled + verticalScale)
    dacSignal = dacSignal.astype(wpt_type)
    
    if max(dacSignal) > 256:
        dacSignal16 = []
        for i in range(0,len(dacSignal)*2):
            dacSignal16.append(0)
            
        j=0
        for i in range(0,len(dacSignal)):
            dacSignal16[j] = dacSignal[i] & 0x0f
            dacSignal16[j+1] = dacSignal[i] >> 8
            j=j+2
        dacSignal = dacSignal16
    
    return(dacSignal);
    
def sineWave(segmentLength, cycles):
    
    time = np.linspace(0, segmentLength-1, segmentLength)
    omega = 2 * np.pi * cycles
    rawSignal = np.sin(omega*time/segmentLength)
    
    return(rawSignal);

def main():

    admin = loadDLL()
    try:
        slotId = getSlotId(admin)

        if not slotId:
            print("Invalid choice!")
        else:
            inst = admin.OpenInstrument(slotId)
            if not inst:
                print("Failed to Open instrument with slot-Id {0}".format(slotId))  # @IgnorePep8
                print("\n")
            else:
                instId = inst.InstrId
                ###########  CALL THE SCPI ########### @IgnorePep8
                instrumentCalls(inst)  #     <<<-------   # @IgnorePep8
                ######################################
                rc = admin.CloseInstrument(instId)
                Validate(rc, __name__, inspect.currentframe().f_back.f_lineno)
    finally:
        if admin is not None:
            rc = admin.Close()
            Validate(rc, __name__, inspect.currentframe().f_back.f_lineno)


def instrumentCalls(inst):

    query_syst_err = True  
    
    res = SendScpi(inst, ":SYST:INF:MODel?")
    
    model = res[1]

    amp = 1    
    segmentLength = (int(5E6/64))*64
    cycles = 5096

    dacSignal = scaleWaveform(sineWave(segmentLength, cycles), model)
     
    sampleRateDAC = 2E9  # override max sample rate above 

    # Module 1
    SendScpi(inst, ":INST:ACTive 1", query_syst_err)

    # get hw option
    SendScpi(inst, "*IDN?")
    # reset - must!
    SendScpi(inst, "*CLS", query_syst_err)
    # reset - must!
    SendScpi(inst, "*RST", query_syst_err)

    # set sampling DAC freq.
    SendScpi(inst, ":FREQ:RAST {0}".format(sampleRateDAC), query_syst_err)  # takes a long time - not required until power cycle @IgnorePep8
    
    # ---------------------------------------------------------------------
    # DAC functions CH 1 
    # ---------------------------------------------------------------------

    # select channel
    SendScpi(inst, ":INST:CHAN 1", query_syst_err)


    start_time = time.time()   

    # load I waveform into instrument
    segNum = 1
    SendScpi(inst, ":TRACe:DEF {0},{1}".format(segNum, len(dacSignal)), query_syst_err)
    SendScpi(inst, ":TRACe:SEL {0}".format(segNum), query_syst_err)

    prefix = ':TRAC:DATA 0,#'
    print(prefix, end=' .. ')

    res = inst.WriteBinaryData(prefix, dacSignal.tobytes())
    Validate(res.ErrCode, __name__, inspect.currentframe().f_back.f_lineno)
    
    print("\n--- download time: %s seconds ---" % (time.time() - start_time))

    res = SendScpi(inst, ":SYST:ERR?", False, False)
    print(res[1])

    # Vpp for output
    SendScpi(inst, ":SOUR:VOLT 0.2")

    # sel segment 1 - play I
    SendScpi(inst, ":SOUR:FUNC:MODE:SEGM {0}".format(segNum), query_syst_err) # Nadav (2021-01-21): 'SEGM' rather than 'SEG'
    # connect ouput
    SendScpi(inst, ":OUTP ON", query_syst_err)
    
    # enble marker 1 CH 1
    # onTime = 4096
    # offTime = int(len(dacSignal_I)/8)-4096
    # markerWave = []
    # for i in range(0,onTime):
        # markerWave.append(1)
    # for i in range(0,offTime):
        # markerWave.append(0)  
    # markerWave = np.uint8(markerWave)
    
    # prefix = ':MARK:DATA 0,#'
    # print(prefix, end=' .. ')
    # res = inst.WriteBinaryData(prefix, markerWave.tobytes())
    # Validate(res.ErrCode, __name__, inspect.currentframe().f_back.f_lineno)

    # SendScpi(inst, ":MARK:VOLT:PTOP 1")
    # SendScpi(inst, ":MARK:SEL 1")
    # SendScpi(inst, ":MARK ON")
    


if __name__ == '__main__':
    print('Process Id = {0}'.format(os.getpid()))
    main()
    print('End of example')
