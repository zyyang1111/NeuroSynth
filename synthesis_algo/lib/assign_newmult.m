function [NMer, NM] = assign_newmult(iter, N)

NMer = fix(N / iter); 
NM = iter; 

if NMer * NM < N
    rest = N - NMer * NM; 
    NMer = cat(1, NMer, 1);
    NM = cat(1, NM, rest);
end

if NMer(1) == 0
    NMer(1) = [];
    NM(1) = [];
end
