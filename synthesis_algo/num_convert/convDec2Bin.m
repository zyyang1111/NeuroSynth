function [binA] = convDec2Bin(decN, bitS, Q)

% Zhiyuan Yang  May 25 2017
% decN is the decimal number, either positive or negative 
% bitS is the length of the binary array 
% Q is the number of decimal bits 

a = round(decN * 2^Q); 
if a < 0 
    sign_flag = 1; 
else
    sign_flag = 0;
end

if abs(a-decN*2^Q) ~= 0
    fprintf('decN = %.6f;\ta = %.6f\n', decN, a/2^Q);
end

% check if there is a overflow 
if a > 2^(bitS-1)-1
    a = 2^(bitS-1)-1;
end
if a < -2^(bitS-1)-1
    a = -2^(bitS-1)-1;
end

if sign_flag 
    a = 2^(bitS) + a; 
end

binA = de2bi(abs(a), bitS, 'left-msb');
 
