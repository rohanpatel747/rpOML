function out = lvperformance(c3dep,showPlot,dispData)
%LVPERFORMANCE Compute Mass to Orbit of Various LVs given a launch C3.
%
%   Assumptions/Warnings:
%   	1. MATLAB Curve Fit Toolbox REQUIRED
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. C3dep    [1x1]     Required Departure C3 (km2/s2)
%       2. showPlot [t/f]     Show Plot of C3 versus Mass for each LV
%       3. dispData [t/f]     Show Data of Mass delivered to console
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. vi  [3x1]        Departure Velocity Vector (km/s)
%       2. vf  [3x1]        Arrival   Velocity Vector (km/s)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Credit/References:
%       1. Dr. Gregory Lantoine (JPL) and Joshua Fofrich (Cal Poly Pomona)
%       2. SLS Block-2 performance (From NASA/JPL NEO Deflection App):
%               https://cneos.jpl.nasa.gov/nda/
%       3. SLS Block 1B:
%               https://sites.nationalacademies.org/cs/groups/ssbsite/documents/webpage/ssb_183081.pdf
%       4. KSC Website for: 
%               Falcon Heavy (expendable & recoverable)
%               Delta IV Heavy
%               Atlas V 551
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    % LV Curve Fit Data
    C3_slsb2   = [15.040; 19.397; 30.569; 50.662; 76.424; 99.743; 130.376] ;
    mass_slsb2 = [39335.88; 36856.18; 30977.73; 22111.90; 13701.39; 8493.88; 4200.52] ;
    C3_sls1b   = [-8; 12; 38; 74; 76; 100; 140] ;
    mass_sls1b = [40000; 30000; 20000; 10172; 9655; 5000; 0] ;
    C3_fhe     = [0; 20; 30; 40 ; 60; 90; 100] ;
    mass_fhe   = [15010; 10115; 8230; 6640 ; 4100; 1425; 755] ;
    C3_fhr     = [0; 20; 30; 40; 50 ;55; 60] ;
    mass_fhr   = [6690; 3845; 2740; 1805; 1005; 660; 320] ;
    C3_d4h     = [-10; 20; 30; 40; 60; 90; 100] ;
    mass_d4h   = [12225; 6995; 5755; 4700; 3000; 1180; 705] ;
    C3_a5551   = [-10; 20; 30; 40; 50; 55; 60] ;
    mass_a5551 = [7275; 4140; 3345; 2670; 2120; 1910; 1695] ;

    % 5th-order Polynomial
    c3 = 0:1:140;
    f1 = fit(C3_slsb2,mass_slsb2,'poly4');
    f2 = fit(C3_sls1b,mass_sls1b,'poly4');
    f3 = fit(C3_fhe,mass_fhe,'poly4');
    f4 = fit(C3_fhr,mass_fhr,'poly4');
    f5 = fit(C3_d4h,mass_d4h,'poly4');
    f6 = fit(C3_a5551,mass_a5551,'poly4');

    % Plotting
    if showPlot
        figure
        hold on
        f0 = xline(c3dep);
        set(f0,'color','k','linewidth',1.5,'HandleVisibility','off');
        f1_ = plot(f1(c3));
        set(f1_,'LineWidth',1.5,'color',[0.4 0.4 0.4]);
        f2_ = plot(f2(c3));
        set(f2_,'LineWidth',1.5,'color','b');
        f3_ = plot(f3(c3(1:101)));
        set(f3_,'LineWidth',1.5,'color','c');
        f4_ = plot(f4(c3(1:61)));
        set(f4_,'LineWidth',1.5,'color','r');
        f5_ = plot(f5(c3(1:101)));
        set(f5_,'LineWidth',1.5,'color','g');
        f6_ = plot(f6(c3(1:61)));
        set(f6_,'LineWidth',1.5,'color','m');
        hold off
        title('Launch Performance Curves','fontsize',14,'interpreter','latex')
        legend({'SLS Block 2','SLS Block 1-B','Falcon Heavy Expendable',...
            'Falcon Heavy Recoverable','Delta IV Heavy','Atlas V 551'}, ...
            'fontsize',12,'interpreter','latex')
        xlabel('Departure $C_3$ ($km^2/s2$)','fontsize',12,'interpreter','latex')
        ylabel('Mass Delivered ($kg$)','fontsize',12,'interpreter','latex')
        ylim([0,5e4])
        set(gcf,'color','w');
        set(gca,'TickLabelInterpreter','latex');
        grid on; box on;
    end
    if dispData
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

    % Output Data
    out                = struct;
    out.slsb2          = f1(c3dep);
    out.slsb1          = f2(c3dep);
    out.falconheavyexp = f3(c3dep);
    out.falconheavyrec = f4(c3dep);
    out.deltaivheavy   = f5(c3dep);
    out.atlasv551      = f6(c3dep);

end

