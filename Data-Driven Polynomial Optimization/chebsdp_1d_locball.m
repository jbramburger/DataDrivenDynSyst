function Bt = chebsdp_1d_locball(d)

% Linear operator for S-procedure on the interval [-1,1] for the Chebyshev
% basis, where we represent [-1,1] as {x: T_0(x)-T_2(x)>=0}

% Written by Giovanni Fantuzzi

At = chebsdp_1d(d-1);
Bt = At;
Bt(:,end+1:end+2) = 0;
T2 = chebpoly(2);
for i = 1:size(Bt,2)-2
    Ti = chebpoly(i-1);
    TiT2 = chebcoeffs(T2.*Ti);
    TiT2 = TiT2(:).';
    TiT2(abs(TiT2)<1e-14) = 0;
    tmp = At(:,i).*TiT2;
    Bt(:,1:i+2) = Bt(:,1:i+2) - tmp;
end