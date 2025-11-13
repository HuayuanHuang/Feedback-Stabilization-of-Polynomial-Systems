%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verification

vP11 = 7.0383;
vP12 = 0;
vP13 = -0.9984;
vP14 = -0.8188;
vP22 = 0.4356;
vP23 = 0.2750;
vP24 = [0.5281; 0.1023];
vP24 = vP24'*P24v;
vP33 = 0.6237;
vP34 = [0.4214; -0.0375];
vP34 = vP34'*P34v;
vP44 = [2.9631; 0.2017];
vP44 = vP44'*P44v;

vPx = [vP11, vP12, vP13, vP14;
       vP12, vP22, vP23, vP24;
       vP13, vP23, vP33, vP34;
       vP14, vP24, vP34, vP44];

vLx1 = -0.4531;
vLx2 = 1.8808;
vLx3 = -0.0707;
vLx4 = [-1.6343; -0.7992; -0.8037; -1.1876];
vLx4 = vLx4'*Lx4v;

vLx = [vLx1,vLx2,vLx3,vLx4];

vR12 = -2*kron(y,Hx*vPx*y);
vR3 = -2*y(4)*vLx*y;
vR4 = kron((jacobian(y'*vPx*y,x1))',Fx);
vRx = [vR12+[vR4;zeros(4,1)];
    vR3];
vRxbar = vRx([2,6,7],:);
vRx([1,2,3,4,5,6,7,8,11,14],:) = [];

vSxMatrix = [-vRx'*(N22\N21) + vRxbar'*[1;1;-1], sqrtSchurN*vRx', epsilon2*y'*vPx; 
    sqrtSchurN*vRx, (vRx'*(N22\N21) - vRxbar'*[1;1;-1])*N22, zeros(7,4);
    epsilon2*vPx*y,zeros(4,7), epsilon2*epsilon3*I_n];

vF = [sos(vPx - epsilon1*I_n), sos(vSxMatrix)];
info2 = optimize(vF,[],ops,[])

