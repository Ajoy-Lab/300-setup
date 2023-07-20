
function segDef = segDef (ip, trace, segmentLength)
  
    scpiWrite (ip, ":TRAC:DEL " + trace);
    scpiWrite (ip, ":TRACe:DEF " + trace + ", " + segmentLength);
    
end