function [xValues, yValues] = CCDF (iqWfm, numOfBins)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    pwrWfm = abs(iqWfm) .* abs(iqWfm);
    meanPwr = 10 * log10(mean(pwrWfm));
    numOfSamples = length(iqWfm);
    papr = 10 * log10(max(pwrWfm)) - meanPwr;
    
    xValues = 0 : (papr / (numOfBins - 1)) : papr;
    yValues = zeros(1, length(xValues));
    
    for i=1:numOfSamples
        if pwrWfm(i) > 0
            pwrDB = 10 * log10(pwrWfm(i)) - meanPwr;
        else
            pwrDB = -100;
        end
        
        for j=1:length(xValues)            
            if pwrDB >= xValues(j)
                yValues(j) = yValues(j) + 1;
            end                      
        end
    end
    
    yValues = yValues * 100 / numOfSamples;    
end

