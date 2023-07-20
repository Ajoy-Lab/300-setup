function retval = myQuantization (myArray, dacRes)
  
  minLevel = 0;
  maxLevel = 2 ^ dacRes - 1;  
  numOfLevels = maxLevel - minLevel + 1;
  
  retval = round((numOfLevels .* (myArray + 1) - 1) ./ 2);
  retval = retval + minLevel;
  
  retval(retval > maxLevel) = maxLevel;
  retval(retval < minLevel) = minLevel;

end
