function [B] = get_mult_use(A)

% yzy : 8/25/2017
% givent a set of weights, calculated the use of multipliers 

B = unique(A); 
for i = 1 : size(B,1)
    id = B(i,1);
    B(i,2) = length(find(A == id));
end