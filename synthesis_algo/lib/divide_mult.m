function [NMer, NM] = divide_mult(P, N)

% Zhiyuan Yang   May 27 2017
% P: the number of multipliers that can be provided 
% N: the number of multiplications 
% NMer: the list of #-multipliers 
% NM: the list of loop for each type of multiplier

if N < P
    NMer = N; 
    NM = 1; 
else
    NMer = []; NM = []; 
    res_N = N; 
    res_P = P; 
    while res_N
        if res_N < res_P 
            NMer = cat(1, NMer, res_N); 
            NM = cat(1, NM, 1); 
            break;
        end
        t = ceil(res_N/res_P); 
        n_mer = fix(res_N/t); 
        NMer = cat(1, NMer, n_mer); 
        NM = cat(1, NM, t); 
        
        res_N = res_N - n_mer * t; 
        res_P = res_P - n_mer; 
    end
end

