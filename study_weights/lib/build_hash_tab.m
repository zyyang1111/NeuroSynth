function [hash_tab] = build_hash_tab(w_u)

for i = 1 : size(w_u,1)
    hash_tab(i).done = 0;
    hash_tab(i).tot_D = 0;
    hash_tab(i).tot_A = 0;
    hash_tab(i).wc = []; 
end