
Kt = .498192/85;
R = 12/85/2;
Kv = (17800/60*2*3.141592)/(12 - 12/85/2*1.4);
J = .0032;
G = 11/34;
%^971's values

%time step
dt = .1;

%not affected by anything but input and inertia
A = [0      1
     0 -Kt/(Kv*R*J*G^2)];
%can directly control acceleration
B = [   0
     Kt/(R*J*G)];
%can observe position
C = [1 0];

%put A and B in discrete time space
A_d = expm(A*dt);
B_d = pinv(A)*(A_d - eye(size(A_d)))*B;

P_K = [.6 .981];
P_L = [.45-.07i,.45+.07i];
K = place(A_d,B_d,P_K);
L = place(A_d.',C.',P_L).';

%initialize
x = [0;0];
y = C*x;
x_hat = [0;0];

t=0;
while t<100;
    Rs = [x_hat(1);1]; 

    u = K*(Rs-x_hat);
    u(u>12) = 12;
    u(u<-12) = -12;

    x_hat = A_d*x_hat +B_d*u + L*(y - C*x_hat);
    x = A_d*x + B_d*u;
    y = C*x;

    t = t+dt;
end;
