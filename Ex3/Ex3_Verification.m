%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verification

% valueP11 = value(P11c)'*P11v;
% valueP12 = value(P12c)'*P12v;
% valueP22 = value(P22c)'*P22v;

valuePx = [0.2965, -0.2330-0.1433*x(1)^2;
    -0.2330-0.1433*x(1)^2, 1.2106+0.2346*x(1)^4];

% valueLx1 = value(Lx1c)'*Lx1v;
% valueLx2 = value(Lx2c)'*Lx2v;

valueLx1 = -0.3329;
valueLx2 = -2.8145-0.8602*x(2)^2-1.4494*x(1)^6;
valueLx = [valueLx1,valueLx2];

vR1 = -2*y(1)*Hx*valuePx*y;
vR2 = -2*y(2)*Hx*valuePx*y;
vR3 = -2*y(2)*valueLx*y;
vR4 = jacobian(y'*valuePx*y,x(1))*Fx;
vRx = [vR1+vR4;
       vR2;
       vR3];

vSxMatrix = [-vRx'*(N22\N21), sqrtSchurN*vRx', epsilon2*y'*valuePx; 
    sqrtSchurN*vRx, vRx'*(N22\N21)*N22, zeros(7,2);
    epsilon2*valuePx*y,zeros(2,7), epsilon2*epsilon3*I_n];

info2 = optimize([sos(vSxMatrix), sos(valuePx - epsilon1*I_n)],[],ops,[])



