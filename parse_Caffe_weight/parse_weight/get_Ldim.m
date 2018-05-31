function [Ldim] = get_Ldim(line)
% get the digit of the line 

digit_idx = isstrprop(line, 'digit');
count = 0;
flag = 0;
for i = 1 : length(digit_idx)-1
    if ~digit_idx(i) && digit_idx(i+1)
        s = i+1;
        flag = 1; 
    end
    if digit_idx(i) && ~digit_idx(i+1)
        if flag
            e = i;
            count = count + 1;
            Ldim(count) = str2num(line(s:e));
            flag = 0;
        end
    end
end
