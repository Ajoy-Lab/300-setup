function segSel = segSel (ip, trace)
    scpiWrite (ip, ":SOUR:FUNC:SEG " + trace)
end