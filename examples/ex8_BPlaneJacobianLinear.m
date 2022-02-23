%% rpOML Examples - B-Plane Jacobian Linear Perturbations
clear; close all; clc; format long g; rpOMLstart();


%% Overview
% This script is intended to demonstrate the acceptable values of the
% perturbation magnitude (to construct a B-Plane Jacobian Matrix) that
% falls in the linear region. Any arbitrary value for the perturbation DV
% magnitude can be used, but some large (>=1 km/s) or small (<1e-14 km/s)
% values will yield unusable numerical partial values.




%% Problem Setup
% Constants
pcd = constants();
mu = pcd.Earth.mu;

% Suppose we have an Earth gravity assist (EGA). The spacecraft is at the
% Sphere of Influence (SOI) with the following state:
r_= [546507.344255845;  -527978.380486028;  531109.066836708];
v_= [-4.9220589268733;   5.36316523097915; -5.22166308425181];
x_= [r_;v_];




%% Compute Nominal B-Plane Values
% The Nominal B-Plane values are shown in the BN structure for the flyby.
% In the command window, type "help bplanefromRV" for a full description of
% the outputs coming from this function. The expression is unsuppressed in
% MATLAB to show the values in the command window.
help bplanefromRV
BN = bplanefromRV(mu,x_)




%% Determining Acceptable Perturbation Vector Magnitude
% We will sample 100 perturbation DV magnitudes from 1e-15 to 1 km/s.
dV = logspace(-15,1,100);

% Find the 2x2 Jacobian Matrix (varying Vx and Vy to find new B*R and B*T)
for i=1:length(dV)
jBjVitr(:,:,i) = bplane_getJacobian(mu, x_, dV(i),'npS',2);
end

% Plot Results
% It is evident that values from 1e-15 to around 1e-11 km/s are unstable
% and yield very different partial numerical values. On the upper end, the
% values begin to diverge after around 1 km/s. The linear region (where the
% least changes in the numerical value versus perturbation magnitude) occur
% between these bounds. This means when computing the Jacobian, these
% values are useful and will yield consistent numerical values.
figure()
hold on
plot(dV,squeeze(jBjVitr(1,1,:)),'linewidth',1.5);
plot(dV,squeeze(jBjVitr(1,2,:)),'linewidth',1.5);
plot(dV,squeeze(jBjVitr(2,1,:)),'linewidth',1.5);
plot(dV,squeeze(jBjVitr(2,2,:)),'linewidth',1.5);
xline(1e-11);
xline(1e-1);
hold off
legend({'\partialB_T/\partialV_\infty_x', ...
    '\partialB_T/\partialV_\infty_y', ...
    '\partialB_R/\partialV_\infty_x', ...
    '\partialB_R/\partialV_\infty_y', ...
    'Linear Region'},'fontsize',12);
set(gca, 'XScale', 'log');
xlim([1e-15 1e1]);
xlabel('Perturbation Magnitude (km/s)','fontsize',12);
ylabel('Numerical Partial Value','fontsize',12);
title('Jacobian Components versus Perturbation Magnitude','fontsize',16); 
grid on; box on; set(gcf, 'color','w');
