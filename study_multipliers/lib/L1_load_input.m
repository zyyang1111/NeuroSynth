function [L1_store_in_size, L1_load_in_size, L1_in_tag] = L1_load_input(rr, cc, ti, trr, tcc, tn, S, K)

L1_store_in_size = trr * tcc * tn * K^2; 
% convert L1_store_in_size to bytes 
L1_store_in_size = L1_store_in_size * 2; 

row_e = rr+trr-1; 
row_in_e = (row_e-1)*S+K; 
row_in_s = (rr-1)*S+1; 
col_e = cc+tcc-1;
col_in_e = (col_e-1)*S+K;
col_in_s = (cc-1)*S+1; 
L1_load_in_size = (row_in_e - row_in_s + 1) * (col_in_e - col_in_s + 1) * tn;
% convert L1_load_in_size to bytes 
L1_load_in_size = L1_load_in_size * 2;

L1_in_tag = [row_in_s, row_in_e;
             col_in_s, col_in_e;
             ti, ti + tn - 1]; 