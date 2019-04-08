%% Assignment 4: Part 1
% 
% This assignment serves as an introduction to circuit simulation and
% compact models. The simulation is done using Modified Nodal Analysis (MNA) on
% circuits. This analysis is done for transient and frequency domain
% simulations, the two simulations run in SPICE applications. The MNA
% technique is applied using common node voltages to create a matrix of
% circuit equations. The nodel equations are created using Kirchoffs
% current law and voltage law. The Nodal matrix is craeted from the
% equation $ G_n*V = I_ns$ where V and I are the voltage and currents of
% the branches. Each device is given a stamp, whereby it's current is
% defined in terms of an applied voltage. The first part of the assignemnt
% provides a circuit with a set of components that will be used for nodal
% analysis. Stamps can be found below: 
close all;
clear; 

R1 =1;
R2 = 2;
R3 = 10;
R4 = 0.1;
R0 = 1000;
cap = 0.25;
L = 0.2;
alpha = 100;

G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
G0 = 1/R0;

%% 
% The Matrices are used to satisfy the equation $C(dV/dt) + GV = F$ where C 
% is a capacitive matrix, looking at energy storing devices and G is the conduction
% matrix, taking into account linear relationships. The G matrix was then created using the nodal equations generated from
% the stamps above. The G matrix was created based on the number of
% equations as well as number of unkowns. Initialization of the G matrix
% can be found below. 

G = zeros(8,8); 
C = zeros(8, 8);

G(1, 1) = -G1;
G(1, 2) =  G1;
G(2, 1) =  G1;
G(1, 3) =  G1;
G(2, 3) = -G1-G2;
G(3, 4) = -G3;
G(2, 7) = -1;
G(3, 7) = 1;
G(4, 3) = 1;
G(4, 4) = -1;
G(5, 6) = G4;
G(5, 7) = -alpha*G4;
G(5, 8) = 1;
G(6, 6) = -G4-G0;
G(6, 7) = alpha*G4;
G(7, 1) = 1;
G(8, 5) = 1; 
G(8, 7) = -alpha;

C(1, 1) = -cap;
C(1, 3) = cap;
C(2, 1) = cap;
C(2, 3) = -cap;
C(4, 7) = -L;
 
%%
% The nodes were then analized in the frequency domain to create similar
% equations, however taking into account the frequency aspects of the
% circuit which includes capacitors and inductors.  initiliazation of the C matrix
% can be found below. The F vector was also initilized which would serve as the resulting
% vector for the simulation. 


F = zeros(8,1); 

% Solve F matrix
Vin = linspace(-10,10,100);
V3 = zeros(length(Vin),1);
V0 = zeros(length(Vin),1);

for i = 1:length(Vin)
    F(7,1) = Vin(i);
    V = G\F;
    V3(i) = V(4);
    V0(i) = V(6);
end

figure(1)
plot(Vin,V3);
xlabel('Vin')
ylabel('V3')
title('V3 vs Vin sweep')

figure(2)
plot(Vin,V0);
xlabel('Vin')
ylabel('V0')
title('V0 vs Vin sweep')

%% AC Gain Plot
% 
% A plot of Vin vs Vo was then plotted for the AC case of the circuit. This
% involved taking into account the frequency components in the circuit. 

 w = 2*pi*linspace(0,80,100);
 V0 = zeros(length(w),1);
 gain = zeros(length(w),1);

for i = 1:length(w)
    s = 1i*w(i);
    M = inv((G +((s).*C)))*F; 
    V0(i) = abs(M(5));
    gain(i) = 20*log10(abs(V0(i))/abs(M(1)));
end

figure(3)
plot(w,V0);
xlabel('w (rads/sec)')
ylabel('V0')
title('AC plot for V0')


figure(4)
semilogx(w,gain);
xlabel('w (rads/sec)')
ylabel('V0/V1 (dB)')
title('Gain');

V0 = zeros(length(w),1);
gain = zeros(length(w),1);

for i = 1:length(gain)
    pert = randn()*0.05;
    C(1,1) = cap*pert;
    C(1,3) = -cap*pert;
    C(2,1) = -cap*pert;
    C(2,3) = cap*pert;

    s = 1i*2*pi*pi;
    M = inv((G +((s).*C)))*F; 
    V0(i) = abs(M(5));
    gain(i) = 20*log10((V0(i))/abs(M(1)));
end

figure(5);
hist(gain,100);
xlabel('Gain')
ylabel('Counts')
title('Hist C')