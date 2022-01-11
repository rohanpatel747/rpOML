%% Launch Vehicle Performance Curves
% C: 21JAN20

clear; close all; clc; format long g;
disp(' ');
disp('- - - Caution - - -');
disp('MATLAB CURVE FIT TOOLBOX REQUIRED');
disp('- - - - - - - - - - ');
disp(' ');

%% Inputs
info = 1;       % Plot mass/c3          (0=off, 1=on)
plt  = 1;       % Mass delivered info   (0=off, 1=on)
c3dep = 7.8^2; % Input departure C3 (km^2/s^2)

%% Launch Vehicle Data
% Courtesy of Dr. Gregory Lantoine (JPL) and Joshua Fofrich (Cal Poly Pomona)
% SLS Block-2 performance (From NASA/JPL NEO Deflection App): 
% https://cneos.jpl.nasa.gov/nda/)
C3_slsb2 = [15.040; 19.397; 30.569; 50.662; 76.424; 99.743; 130.376] ;
mass_slsb2 = [39335.88; 36856.18; 30977.73; 22111.90; 13701.39; 8493.88; 4200.52] ;
%SLS Block 1B
%https://sites.nationalacademies.org/cs/groups/ssbsite/documents/webpage/ssb_183081.pdf
C3_sls1b = [-8; 12; 38; 74; 76; 100; 140] ;
mass_sls1b = [40000; 30000; 20000; 10172; 9655; 5000; 0] ;
%Falcon Heavy Expendable
%KSC Site
C3_fhe = [0; 20; 30; 40 ; 60; 90; 100] ;
mass_fhe = [15010; 10115; 8230; 6640 ; 4100; 1425; 755] ;
%Falcon Heavy Recoverable
%KSC Site
C3_fhr = [0; 20; 30; 40; 50 ;55; 60] ;
mass_fhr = [6690; 3845; 2740; 1805; 1005; 660; 320] ;
%Delta Quatro Heavy
%KSC Site
C3_d4h = [-10; 20; 30; 40; 60; 90; 100] ;
mass_d4h = [12225; 6995; 5755; 4700; 3000; 1180; 705] ;
%Atlas Cinco 551
%KSC Site
C3_a5551 = [-10; 20; 30; 40; 50; 55; 60] ;
mass_a5551 = [7275; 4140; 3345; 2670; 2120; 1910; 1695] ;

%% Fit the data using a 5th-order polynomial
c3 = 0:1:140;
f1 = fit(C3_slsb2,mass_slsb2,'poly4');
f2 = fit(C3_sls1b,mass_sls1b,'poly4');
f3 = fit(C3_fhe,mass_fhe,'poly4');
f4 = fit(C3_fhr,mass_fhr,'poly4');
f5 = fit(C3_d4h,mass_d4h,'poly4');
f6 = fit(C3_a5551,mass_a5551,'poly4');

%% Plot/Outputs 
y = [0,5e4];
x = [c3dep,c3dep];
if plt == 1
    figure
    hold on
    f0 = line(x,y);
    set(f0,'color','k','linewidth',1.5,'HandleVisibility','off');
    f1_ = plot(f1(c3));
    set(f1_,'LineWidth',2,'color',[0.4 0.4 0.4]);
    f2_ = plot(f2(c3));
    set(f2_,'LineWidth',2,'color','b');
    f3_ = plot(f3(c3(1:101)));
    set(f3_,'LineWidth',2,'color','c');
    f4_ = plot(f4(c3(1:61)));
    set(f4_,'LineWidth',2,'color','r');
    f5_ = plot(f5(c3(1:101)));
    set(f5_,'LineWidth',2,'color','g');
    f6_ = plot(f6(c3(1:61)));
    set(f6_,'LineWidth',2,'color','m');
    hold off
    title('Launch Performance Curves','fontsize',20)
    legend({'SLS Block 2','SLS Block 1-B','Falcon Heavy Expendable',...
        'Falcon Heavy Recoverable','Delta IV Heavy','Atlas V 551'}, ...
        'fontsize',14)
    xlabel('Launch C3 (km^2/s2)','fontsize',14)
    ylabel('Launch mass (kg)','fontsize',14)
    ylim([0,5e4])
    set(gcf,'color','w');
    grid on
end
if info == 1
    disp('------------------------------------------------------')
    disp(['Launch Vehicle Performance for ',num2str(c3dep),' km^2/s^2'])
    disp(' ')    
    disp(['SLS Block 2: ', num2str(f1(c3dep)),' kg'])
    disp(['SLS Block 1B: ', num2str(f2(c3dep)),' kg'])
    disp(['Falcon Heavy Expendable: ', num2str(f3(c3dep)),' kg'])
    disp(['Falcon Heavy Recovery: ', num2str(f4(c3dep)),' kg'])
    disp(['Delta IV Heavy: ', num2str(f5(c3dep)),' kg'])
    disp(['Atlas V 551: ', num2str(f6(c3dep)),' kg'])
    disp('------------------------------------------------------')
end