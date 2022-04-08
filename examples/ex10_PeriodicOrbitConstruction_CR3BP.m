%% rpOML Examples - Constructing Periodic Orbits with a Visualization
clear; clc; close all; format long g; 

%% Overview
% This script is intended to show a single shooting method can converge on
% a periodic orbit solution in the CR3BP. We will start with some 'good'
% initial guesses for the state at the y-axis crossing and then use the STM
% to correct our orbit to achieve periodicity. This is defined as when the
% next y-axis crossing has a dx and dz components of 0. Because we are
% using a numerical method that relies on linearization (use of the STM), 
% we cannot reach 0 exactly, though it is possible to converge on a
% solution within some tolerance parameter (using 1e-10).

% The problem is based loosly off the James Webb Space Telescope's Halo
% Orbit. The initial conditions and single shooting method is provided by
% professor Kate Davis at CU Boulder in ASEN6008.





%% Problem Setup
% We define some tolerance criteria to determine if the orbit is periodic.
% Also, we want to evaluate the state at the subsequent y-axis crossing.
% Therefore, in our integrator's options, we will call "Events" and specify
% to terminate integration at the y-axis using the library's termination
% function.

tolVal = 1e-10;
options = odeset('Events', @cr3bp_yaxstopcond);





%% Initial Conditions
% For the following Problem we are in the Sun-Earth/Moon 3-Body system.
% We can query the constantscr3bp() function in the library to get our
% system parameters:
pcd3 = constantscr3bp();

% CR3BP System (Sun-Earth/Moon) Constants
sem = pcd3.SunEarthMoon;
mu  = sem.mu;


% Our initial conditions are as follows:
x0 = 1.0076245;
z0 = -0.0029;
dy0= 0.013;





%% Initial State Vector (42x1)
% The full state vector contains the position [3x1], velocity [3x1], and
% the initial time state transition matrix [6x6] (which is identity at t=0.
% This matrix is to be turned into a [36x1] to turn it into a vector for
% the initial conditions.
X0 = [x0; 0; z0; 0; dy0; 0];
J0 = eye(6);
J0_= reshape(J0, 36,1);
X0_ = [X0; J0_];





%% Running Correction
% With the STM we are able to correct either the x0 and dy0 initial states
% or we are able to correct z0 and dy0. The following boolean will tell the
% algorithm to correct either of these to compute our solution.
correct_z0 = false;         % <--- false=correct x0 and dy0.


% Visualize Corrections
hold on
scatter3(sem.b1(1), sem.b1(2), sem.b1(3), 'r', 'filled')
scatter3(sem.b2(1), sem.b2(2), sem.b2(3), 'k', 'filled')

% Compute and Plot Each Iteration
ittr = 0;
uwm(1:2,1) = 1;
while (norm(uwm(1))>tolVal) && (norm(uwm(2))>tolVal)

% Integrate w/ ODE113
[~,  x,~,~,~] = ode113(@(t,Y) eomCR3BPwSTM(t,Y,mu), [0 8*pi], X0_,options);


% State and STM at 1st Y-Axis Crossing
YT2_ = x(end,:);
stm  = cr3bp_getSTMfromY(YT2_, 1);  % STM(T/2, 0)


% Get Accels at Time T/2
dYdT2 = eomCR3BP(0, YT2_, mu);


% Unwanted Motion
ddx = X0(4) - x(end,4);
ddz = X0(6) - x(end,6);
uwm = [ddx; ddz];


% Correct Unwanted Motion
%correct_z0 = true;

if correct_z0
    % Hold x0 fixed, vary z0 and dy0
    Y0dy0_1 = [stm(4,3) stm(4,5);
               stm(6,3) stm(6,5)];
    Y0dy0_2 = [stm(2,3) stm(2,5)];
else
    % Hold z0 fixed, vary x0 and dy0
    Y0dy0_1 = [stm(4,1) stm(4,5);
               stm(6,1) stm(6,5)];
    Y0dy0_2 = [stm(2,1) stm(2,5)];
end
uwmc = Y0dy0_1 - (1/dYdT2(2))*[dYdT2(4); dYdT2(6)]*Y0dy0_2;


% Compute Correction
uwmCor = inv(uwmc)*uwm;
if correct_z0
    X0_ = [X0_(1);           X0_(2); X0_(3)+uwmCor(1);
           X0_(4); X0_(5)+uwmCor(2);           X0_(6)];
else
    X0_ = [X0_(1)+uwmCor(1);           X0_(2); X0_(3);
                     X0_(4); X0_(5)+uwmCor(2); X0_(6)];
end


% Update I.C.
X0_(7:42) = reshape(eye(6), 36,1);

ittr = ittr+1;



if ittr==1
    plot3(x(:,1),x(:,2),x(:,3),'r','linewidth',1.5);
else
    plot3(x(:,1),x(:,2),x(:,3),'Color',[0.6 0.6 0.6],'linewidth',0.5)
end

end





%% Resulting Periodic Orbit (and subsequent Initial Conditions Vector X0_)
% Final Trajectory (Converged)
disp(['Number of Iterations: ', num2str(ittr)]);
disp(' ');
disp(['Unwanted Motion Components:']);
disp(['    dx0(T/2): ',num2str(uwm(1),15)]);
disp(['    dz0(T/2): ',num2str(uwm(2),15)]);
disp(' ');
[~,  x,te,~,~] = ode113(@(t,Y) eomCR3BPwSTM(t,Y,mu), [0 8*pi], X0_,options);
plot3(x(:,1),x(:,2),x(:,3),'g','linewidth',1.5)

hold off
grid on; box on; axis equal; set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex');
xlabel('$X_{ND}$','fontsize',12,'interpreter','latex');
ylabel('$Y_{ND}$','fontsize',12,'interpreter','latex');
zlabel('$Z_{ND}$','fontsize',12,'interpreter','latex');
title('Periodic Orbit about Sun-Earth/Moon L2', ...
    'fontsize',16,'interpreter','latex');
disp('Plot is 3-D! Rotate it in the Figure Window');
xlim([0.999 1.017]);















