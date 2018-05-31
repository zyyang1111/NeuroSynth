function [L2_in_size, L2_in_tag] = L2_load_input(row, col, ti, tr, tc, tn, S, K)

row_e = row+tr-1; 
row_in_e = (row_e-1)*S+K; 
row_in_s = (row-1)*S+1; 
col_e = col+tc-1;
col_in_e = (col_e-1)*S+K;
col_in_s = (col-1)*S+1; 
L2_in_size = (row_in_e - row_in_s + 1) * (col_in_e - col_in_s + 1) * tn;
% convert L2_in_size in bytes 
L2_in_size = L2_in_size * 2;
L2_in_tag = [row_in_s, row_in_e;
             col_in_s, col_in_e;
             ti, ti + tn - 1]; 