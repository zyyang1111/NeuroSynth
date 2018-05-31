function [NMer, NM] = assgin_delay(D, d_i, N)

% Zhiyuan Yang   May 27 2017
% D: the delay we want to meet 
% d_i: the delay of the current multiplier
% N: the number of multiplications 
% NMer: the list of #-multipliers 
% NM: the list of loop for each type of multiplier

if fix(D/d_i) >= N
    NMer = 1; 
    NM = N;
else
    t = ceil(N / fix(D/d_i)); 
    n = fix(N / t); 
    res = mod(N, t); 
    
    NMer = t; 
    NM = n; 
    if res 
        NMer = cat(1, NMer, 1);
        NM = cat(1, NM, res);
    end
end