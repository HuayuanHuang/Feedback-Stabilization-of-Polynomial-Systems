clearvars -except out;

time = out.data.Time(:,1)';
Data_U = out.data.Data(:,1)'; 
Data_DX = [out.data.Data(:,2)';
    out.data.Data(:,3)';
    out.data.Data(:,4)';
    out.data.Data(:,5)'];  
Data_F = [out.data.Data(:,7)';
    out.data.Data(:,8)';
    out.data.Data(:,9)';
    out.data.Data(:,10)']; 
Data_X = [out.data.Data(:,6)';
    out.data.Data(:,7)';
    out.data.Data(:,8)';
    out.data.Data(:,9)'];  

T = 5;
Omega = 0.2;

% i = randperm(501, T);
% i = sort(i);
i = [82,165,235,271,473];

dataU = Data_U(:,i);
dataDX = Data_DX(:,i);
dataF = Data_F(:,i);
dataX = Data_X(:,i);


n = 4;
x1 = sdpvar(3,1);
x2 = sdpvar(1,1);
x = [x1;x2];
Fx = [x(2);x(3);x(4);x(2)^3];
Hx = [zeros(3,1),eye(3);0,x(2)^2,0,0];
I_n = eye(n);

y = sdpvar(4,1);

% P(x)
[P11,P11c,P11v] = polynomial(x1,0);
[P12,P12c,P12v] = polynomial(x1,0);
[P13,P13c,P13v] = polynomial(x1,0);
[P14,P14c,P14v] = polynomial(x1,0);

[P22,P22c,P22v] = polynomial(x1,0);
[P23,P23c,P23v] = polynomial(x1,0);
% [P24,P24c,P24v] = polynomial(x1,2);
P24v = [1;x(2)^2];
P24c = sdpvar(size(P24v,1),1);
P24 = P24c'*P24v;
[P33,P33c,P33v] = polynomial(x1,0);
% [P34,P34c,P34v] = polynomial(x1,2);
P34v = [1;x(2)^2];
P34c = sdpvar(size(P34v,1),1);
P34 = P34c'*P34v;
% [P44,P44c,P44v] = polynomial(x1,4);
P44v = [1;x(2)^4];
P44c = sdpvar(size(P44v,1),1);
P44 = P44c'*P44v;

Px = [P11, P12, P13, P14;
    P12, P22, P23, P24;
    P13, P23, P33, P34;
    P14, P24, P34, P44];

%L(x)
[Lx1,Lx1c,Lx1v] = polynomial(x,0);
[Lx2,Lx2c,Lx2v] = polynomial(x,0);
[Lx3,Lx3c,Lx3v] = polynomial(x,0);
% [Lx4,Lx4c,Lx4v] = polynomial(x,6);
Lx4v = [1;x(2)^2*x(3)^2;x(2)^2*x(4)^2;x(2)^6];
Lx4c = sdpvar(size(Lx4v,1),1);
Lx4 = Lx4c'*Lx4v;
Lx = [Lx1,Lx2,Lx3,Lx4];


X = [dataDX(1,:)';dataDX(2,:)';dataDX(3,:)';dataDX(4,:)'];

D = [kron(I_n,dataF);
    zeros(1,3*T),dataU];
Dbar = D([2,6,7],:);
D([1,2,3,4,5,6,7,8,11,14],:) = [];

Xbar = X-Dbar'*[1;1;-1];

Phi11 = Omega^2*T;

N11 = Phi11 - Xbar' * Xbar;
N12 = Xbar' * D';
N21 = N12';
N22 = -D * D';
SchurN = N11 - N12 /N22 * N21;
sqrtSchurN = sqrt(SchurN);

R12 = -2*kron(y,Hx*Px*y);
R3 = -2*y(4)*Lx*y;
R4 = kron((jacobian(y'*Px*y,x1))',Fx);
Rx = [R12+[R4;zeros(4,1)];
    R3];
Rxbar = Rx([2,6,7],:);
Rx([1,2,3,4,5,6,7,8,11,14],:) = [];

%epsilonx
epsilon1 = 0.1;
epsilon2 = 0.01; 
epsilon3 = 2*(x(1)^2 + x(2)^2 + x(3)^2 + 1);

SxMatrix = [-Rx'*(N22\N21) + Rxbar'*[1;1;-1], sqrtSchurN*Rx', epsilon2*y'*Px; 
    sqrtSchurN*Rx, (Rx'*(N22\N21) - Rxbar'*[1;1;-1])*N22, zeros(7,4);
    epsilon2*Px*y,zeros(4,7), 0.5*epsilon2*epsilon3*I_n];


F = [sos(Px - epsilon1*I_n), sos(SxMatrix)];

ops = sdpsettings('solver','mosek');
info = optimize(F,[],ops,[P11c;P12c;P13c;P14c;P22c;P23c;P24c;P33c;P34c;P44c;Lx1c;Lx2c;Lx3c;Lx4c])


value(P11c)
value(P12c)
value(P13c)
value(P14c)
value(P22c)
value(P23c)
value(P24c)
value(P33c)
value(P34c)
value(P44c)
value(Lx1c)
value(Lx2c)
value(Lx3c)
value(Lx4c)











