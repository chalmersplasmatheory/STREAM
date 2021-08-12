clc
% DYON
M = 2.49e-6;
L_MK2 = 9.1e-6;
L_p = 5.4e-6;
R_MK2 = 7.5e-4;
% STREAM (todo)
% ...

R_p = @(t) 1e-5;
V = @(t) 11;

[t,I] = ode45(@(t,I) diffsys(t,I,M,L_MK2,L_p,R_MK2,R_p,V),[0 0.5],[1000; 0]);

figure(1)
clf
plot(t,I(:,1),t,I(:,2))
xlabel('Time t');
ylabel('Solution I');
legend('I_p','I_MK2')
axis([0 0.5 0 7.5e5])

M = 4.17e-6;
L_MK2 = 4.17e-6;
L_p = 8.3e-6;
R_MK2 = 3.4e-4;

[t,I] = ode45(@(t,I) diffsys(t,I,M,L_MK2,L_p,R_MK2,R_p,V),[0 0.5],[1000; 0]);

figure(2)
clf
plot(t,I(:,1),t,I(:,2))
xlabel('Time t');
ylabel('Solution I');
legend('I_p','I_MK2')
axis([0 0.5 0 7.5e5])
% LEK
M = 4.17e-6;
L_MK2 = 4.17e-6;
L_p = 8.3e-6;
R_MK2 = 2*3.4e-4;

[t,I] = ode45(@(t,I) diffsys(t,I,M,L_MK2,L_p,R_MK2,R_p,V),[0 0.5],[1000; 0]);

figure(3)
clf
plot(t,I(:,1),t,I(:,2))
xlabel('Time t');
ylabel('Solution I');
legend('I_p','I_MK2')
axis([0 0.5 0 7.5e5])

function svar = diffsys(t,I,M,L_MK2,L_p,R_MK2,R_p,V)

%ekv1=1/(M - L_MK2 * L_p / M)*(V(t) - R_p(t)   * I(1) - L_p   / M * ( V(t) - R_MK2 * I(2)));
%ekv2=1/(M - L_MK2 * L_p / M)*(V(t) - R_MK2 * I(2) - L_MK2 / M * ( V(t) - R_p(t)   * I(1)));
ekv1=1/(M - L_MK2 * L_p / M)*(V(t) - R_MK2   * I(2) - L_MK2 / M * ( V(t) - R_p(t) * I(1)));
ekv2=1/(M - L_MK2 * L_p / M)*(V(t) - R_p(t)  * I(1) - L_p   / M * ( V(t) - R_MK2  * I(2)));

svar = [ekv1; ekv2];
end