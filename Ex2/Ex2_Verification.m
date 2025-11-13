% Verification

% valueP11 = value(P11c)'*P11v;
% valueP12 = value(P12c)'*P12v;
% valueP22 = value(P22c)'*P22v;

valueP11 = 0.4547*P11v;
valueP12 = [-0.6397, -0.2666]*P12v;
valueP22 = [2.7295, 0.4405]*P22v;

valuePx = [valueP11,valueP12;
           valueP12,valueP22];

% valueLx1 = value(Lx1c)'*Lx1v;
% valueLx2 = value(Lx2c)'*Lx2v;

valueLx1 = -2.3608*Lx1v;
valueLx2 = [-1.6172, -1.4640, -1.8819]*Lx2v;
valueLx = [valueLx1,valueLx2];

vR11 = jacobian(y'*valuePx*y,x(1))*Fx;
vR12 = -2*y(1)*Hx*valuePx*y;
vR2 = -2*y(2)*Hx*valuePx*y;
vR3 = -2*y(2)*valueLx*y;
vR4 = jacobian(y'*valuePx*y,x(1))*Fx;
vRxy = [vR11+vR12;
        vR2;
        vR3];

vSxMatrix = [-vRxy'*(N22\N21), sqrtSchurN*vRxy', epsilon2*y'*valuePx; 
    sqrtSchurN*vRxy, vRxy'*(N22\N21)*N22, zeros(7,2);
    epsilon2*valuePx*y,zeros(2,7), epsilon2*epsilon3*I_n];

info2 = optimize([sos(valuePx - epsilon1*I_n),sos(vSxMatrix)],[],ops,[])








