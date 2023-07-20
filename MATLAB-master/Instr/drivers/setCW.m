    
function setCW = setCW(ip, cwFreq)
    scpiWrite (ip, ":SOUR:CFR " + cwFreq)
    scpiWrite (ip, ":SOUR:MODE IQM")
end
