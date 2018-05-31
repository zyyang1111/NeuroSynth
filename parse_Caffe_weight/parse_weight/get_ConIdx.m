function [idex] = get_ConIdx(line)

digit_idx = isstrprop(line, 'digit');
idex = str2num(line(digit_idx == 1)); 
idex = idex + 1;