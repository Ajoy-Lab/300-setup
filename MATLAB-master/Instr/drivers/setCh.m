    
function setCh = setCh(ip, channel)
    scpiWrite (ip, ":INST:CHAN " + channel)
    scpiWrite (ip, ":INIT:CONT ON")
    scpiWrite (ip, ":SOUR:VOLT 1.0")
end




    