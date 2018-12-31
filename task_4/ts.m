function [R_z_Tsmooth ] = ts(R_mat, Nc, Q); 
%**************************************************** 
% Function to do temporal smoothing on a Matrix  (Equation 18) 
% Written by Farrukh Rashid 
% INPUTS 
% R_mat = Covariance Matrix which is to be smoothed  
% Nc = Code Size (31) 
% Q = Number of submatrices ( >= Ki) 
 
% OUTPUTS 
% R_z_Tsmooth = Covariance Matrix 
%**************************************************** 
d = 2*Nc - 1 - Q;    % Length of the overlapping submatrices  
N = size(R_mat,1)/(2*Nc); % (N = 5) 
 
Num_segment = N*N; % Total number of segmented matrices  
 
for i = 0 : N - 1 , % Parse the Rows of R_mat (5 rows) 
    for j = 0: N - 1, % Parse the columns of R_mat (5 columns) 
         
        %% Extract the segmented matrix of size 2*Nc 
     R_seg_mat(:,:)=R_mat((1+i*2*Nc):((i+1)*2*Nc),(1+j*2*Nc):((j+1)*2*Nc)); 
         
        Rz_acc_mat = zeros(d,d);  
        for k = 1:Q, % For each segment parse Average a set of Q submarices 
             
            %% Extract the Submatrices (N_1 x d) 
            Rz_sub_mat(:,:) = R_seg_mat( k:(k+d-1), k:(k+d-1));  
 
            %% Accumulate the Submatrices 
            Rz_acc_mat = Rz_acc_mat + Rz_sub_mat; 
        end 
 
        % Divide by the number of submatrices 
        Rz_acc_mat = Rz_acc_mat/Q;  
         
        R_z_Tsmooth((1+i*d):((i+1)*d),(1+j*d):((j+1)*d))=Rz_acc_mat(:,:);  
    end   
end 
flag = 1;
