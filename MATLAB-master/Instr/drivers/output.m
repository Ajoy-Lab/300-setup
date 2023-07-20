
function output = output(ip, state)
    if state == 1
        scpiWrite (ip, ":OUTP ON")
    else
        scpiWrite (ip, ":OUTP OFF")
    end