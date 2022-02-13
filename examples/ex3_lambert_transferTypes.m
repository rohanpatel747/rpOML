%% rpOML Examples - Studying Transfer Types Between Earth and Venus
clear; clc; format long g; rpOMLstart();

%% Overview
% This script is intended to show how the transfer time of flight varies
% between two bodies as a function of the universal variable. 


%% Problem Setup
% Get Planetary Constants Data
pcd = constants();
mu  = pcd.Sun.mu;

% Specify Body 3 (Earth) and Body 2 (Mars)
bdy1 = 3;
bdy2 = 2;

% Get state vector of the bodies at some arbitrary times
sv1 = getStatePlanet(bdy1,2454085.5,'meeus').x;
sv2 = getStatePlanet(bdy2,2454085.5+450,'meeus').x;

% For this problem, we will see how the tranfer time of flight varies for
% different orbit conics and revolutions. The universal variable (psi) is
% the driving variable in the computation of the flight time. Below is how
% the value of psi defines the conic type and revolutions
%{

***NOTE***
There is no clear distinction between Type I/II transfers in this
demonstration. This is due to the nature of the universal variable.

Bounds of Psi Based on Conic:

    Hyperbolic :     -4*pi   < psi < 0.00
    Parabolic  :               psi = 0.00
    Elliptical : 
        Type I/II   : 0.00   < psi < 4*pi^2
        Type III/IV : 4*pi^2 < psi < 4*(2^2)*(pi^2)

%}
% The values of psi are setup below according to the values above.
[psiL, psiU] = lambert_getPsiNRev(1);
psi_hyp = (-4*pi):0.05:0.000;
psi_par = 0.000;
psi_t12 = 0.00:0.05:psiL;
psi_t34 = (psiL+0.001):0.05:(psiU);
psi = [psi_hyp, psi_par, psi_t12, psi_t34];



%% Computing TOF from Psi Values
% Find TOF givne Psi Value
for i=1:length(psi)
    tof(i) = (lambert_getTOFfromPsi(psi(i), sv1, sv2, mu))/86400;
end

% Plot Figure
figure()

idx1 = length(psi_hyp);
idx2 = idx1+1;
idx3 = idx2+length(psi_t12);
idx4 = idx3+length(psi_t34);

hold on
plot(psi_hyp, tof(1:idx1),'b');                 % Hyperbolic
plot(psi_t12, tof(idx2+1:idx3),'g');            % Type I/II
plot(psi_t34, tof(idx3+1:idx4),'r');            % Type III/IV
scatter(psi_par, tof(idx2),'r','filled');       % Parabolic
xline(0,'--')                                   % Hyp./El. Dividing Line
xline(psiL,'--')                                % Type I/II & III/IV
xline(psiU,'--')                                % Type III/IV & V/VI
hold off
grid on; box on; set(gcf,'color','w');
ylim([0 1200]);
ylabel('TOF (days)','fontsize',12,'interpreter','latex')
xlabel('$\Psi$ (radians)','fontsize',12,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex');
title(['Transfer $\Psi$ versus TOF Between Bodies: ',num2str(bdy1), ...
    ' and ', num2str(bdy2)],'fontsize',14, 'interpreter','Latex');
legend({'Hyperbolic','Type I/II','Type III/IV','Parabolic'}, ...
    'fontsize',12,'interpreter','latex','location','southeast');



















