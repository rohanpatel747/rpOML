%% rpOML : Verify PropKepElip Function
%  C: 07APR22

clear; close all; clc; format long g; rpOMLstart();
pcd = constants();
mu = pcd.Earth.mu;


%% Problem Setup
r1 =  7500;
r2 = 10000;


x1_ = [r1; 0; 0; 0; sqrt(mu/r1); 0];
if true
    % Elliptical Solution
    x2_ = [0; -r2; 0; sqrt(mu/r2);  0; 0];
    dt = 86400/16;
else
    % Hyperbolic Solution
    x2_ = [0; r2; 0; -sqrt(mu/r2);  0; 0];
    dt = 86400/100;
end



%% Lambert Solution
l1 = lambert0rev(x1_,x2_, dt, mu);
xl1_ = [x1_(1:3); l1.vi];


elementsData = conv_state2ele(xl1_,mu,true);


xi_ = xl1_;


%% Testing Function
j=2; xf_(1,1:6) = xi_.';
for i=1:1:dt
    xf_(j,1:6) = propKep(xi_, i, mu).';
    j=j+1;
end








%% Integration Checks
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8) ;
[~,s1] = ode45(@(t,Y) eom2BP(t,Y,mu), [0,  2*pi*sqrt((r1^3)/mu)], x1_ , options);
[~,s2] = ode45(@(t,Y) eom2BP(t,Y,mu), [0,  2*pi*sqrt((r2^3)/mu)], x2_ , options);
[~,sL] = ode45(@(t,Y) eom2BP(t,Y,mu), [0,                    dt], xl1_, options);


% These Must Match
disp('Difference Between Lambert/Integration and PropKep()')
for i=1:6
    disp(num2str(sL(end,i)-xf_(end,i), 15));
end



%% Plotting
figure()
hold on
scatter3(0,0,0,'filled');
scatter3(x1_(1),x1_(2),x1_(3),'filled');
scatter3(x2_(1),x2_(2),x2_(3),'filled');
plot3(s1(:,1),s1(:,2),s1(:,3),'k')
plot3(s2(:,1),s2(:,2),s2(:,3),'k')
plot3(sL(:,1),sL(:,2),sL(:,3),'r')
plot3(xf_(:,1),xf_(:,2),xf_(:,3),'b')
e__ = 50000*elementsData.e_;
quiver3(0,0,0,e__(1),e__(2),e__(3),'autoscale','off');
hold off
grid on; box on; axis equal;





















