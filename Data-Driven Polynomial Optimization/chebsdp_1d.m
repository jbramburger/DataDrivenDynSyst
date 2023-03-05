function At = chebsdp_1d(d)

% Create SDP data matrices for an SOS constraint in the chebyshev basis
% Precisely, returns coefficient matrices such that
%
% [T_0]*[T_0 T_1 ... T_d] = A_0 T_0 + A_1 T_1 + ... + A_2d T_2d
% [T_1]
% [...]
% [T_d]
%
% The matrices A0,...A_2d are stored as columns of a sparse matrix At.

% Written by Giovanni Fantuzzi

[i,j] = meshgrid(1:d+1,1:d+1);
s = sub2ind([d+1,d+1],i,j);
i = i-1; j = j-1;
ipj = i+j;
imj = abs(i-j);
row = [];
col = [];
val = [];
for k = 0:2*d+1
    B = zeros(d+1,d+1);
    idx= [s(ipj==k); s(imj==k)];
    row = [row; idx];
    col = [col; k*ones(size(idx))+1];
    val = [val; 0.5*ones(size(idx))];
end
At = sparse(row,col,val,(d+1)^2,2*d+1);
