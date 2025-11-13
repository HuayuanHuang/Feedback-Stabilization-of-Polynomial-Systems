clearvars -except out;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = out.data.Time(:,1)';
Data_U = out.data.Data(:,1)'; 
Data_DX = [out.data.Data(:,2)';
    out.data.Data(:,3)'];  
Data_F = [out.data.Data(:,4)';
    out.data.Data(:,5)';
    out.data.Data(:,6)'];  

T=4; 
Omega = 0.2; 

% i = randperm(501, T);
% i = sort(i);
i = [40,108,202,465];

dataU = Data_U(:,i);
dataDX = Data_DX(:,i);
dataF = Data_F(:,i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=2;
x = sdpvar(n,1);
Fx = [x(1);x(2);x(1)^2];
I_n = eye(n);
Hx = [I_n;x(1),0];

y = sdpvar(2,1);

% P(x)
% [P11,P11c,P11v] = polynomial(x(1),0);
% [P12,P12c,P12v] = polynomial(x(1),2);
% P21 = P12;
% [P22,P22c,P22v] = polynomial(x(1),4);

[P11,P11c,P11v] = polynomial(x(1),0);

P12v = [1;x(1)^2];
P12c = sdpvar(size(P12v,1),1);
P12 = P12c'*P12v;
P21 = P12;

P22v = [1; x(1)^4];
P22c = sdpvar(size(P22v,1),1);
P22 = P22c'*P22v;

Px = [P11, P12;
    P21, P22];

% L(x)
% [Lx1,Lx1c,Lx1v] = polynomial(x,3); 
% [Lx2,Lx2c,Lx2v] = polynomial(x,6);

Lx1v = 1;
Lx1c = sdpvar(size(Lx1v,1),1);
Lx1 = Lx1c'*Lx1v;

Lx2v = [1; x(2)^2; x(1)^6];
Lx2c = sdpvar(size(Lx2v,1),1);
Lx2 = Lx2c'*Lx2v;

Lx = [Lx1,Lx2];

% R(x,y)
R11 = jacobian(y'*Px*y,x(1))*Fx;
R12 = -2*y(1)*Hx*Px*y;
R2 = -2*y(2)*Hx*Px*y;
R3 = -2*y(2)*Lx*y;
Rxy = [R11+R12;
    R2;
    R3];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = [dataDX(1,:)';dataDX(2,:)'];
D = [dataF, zeros(3,T);
    zeros(4,T),[dataF;dataU]];
Phi11 = Omega^2*T;

%N matrix
N11 = Phi11 - X' * X;
N12 = X' * D';
N21 = N12';
N22 = -D * D';
SchurN = N11 - N12 /N22 * N21;
sqrtSchurN = sqrt(SchurN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon1 = 0.1;
epsilon2 = 0.01; 
epsilon3 = 2*(x(1)^2 + 1);

SxMatrix = [-Rxy'*(N22\N21), sqrtSchurN*Rxy', epsilon2*y'*Px; 
    sqrtSchurN*Rxy, Rxy'*(N22\N21)*N22, zeros(7,2);
    epsilon2*Px*y,zeros(2,7), 0.5*epsilon2*epsilon3*I_n];

F = [sos(Px - epsilon1*I_n), sos(SxMatrix)];

ops = sdpsettings('solver','mosek');
info = optimize(F,[],ops,[P11c;P12c;P22c;Lx1c;Lx2c])

value(P11c)
value(P12c)
value(P22c)
value(Lx1c)
value(Lx2c)




