
function output = output(ip, openPXI, doReset, setSample)

    if openPXI == 1
        scpiWrite (ip, "OPENPXI")
        scpiWrite (ip, ":INST:ACTive 1")
    end
    scpiWrite (ip, ":*OPT?")
    if doReset == 1
        scpiWrite (ip, ":*RST")
    end
    if setSample > 0
        scpiWrite (ip, ":FREQ:RAST " + setSample)
    end
    scpiWrite (ip, ":SOUR:FUNC:MODE ARB")
    scpiWrite (ip, ":TRAC:DEF:TYPE NORM")
end
    



