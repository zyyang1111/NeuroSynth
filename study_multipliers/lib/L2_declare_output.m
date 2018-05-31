function [L2_out_size, L2_out_tag] = L2_declare_output(row, col, to, tr, tc, tm)

L2_out_size = tr * tc * tm; 
% convert L2_out_size to bytes 
L2_out_size = L2_out_size * 2; 
L2_out_tag = [row, row+tr-1;
              col, col+tc-1;
              to, to+tm-1]; 