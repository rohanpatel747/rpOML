function f = pkchp_plt(in)
%PKCHP_PLT Summary of this function goes here
%   Detailed explanation goes here

    c3dep   = in.c3dep;
    vinfarr = in.vinfarr;
    vinfdep = in.vinfdep;
    RLA     = in.RLAdep;
    DLA     = in.DLAdep;
    tof     = in.tof;
    ntof    = in.ntof;
    dp1     = in.dp1;
    dp2     = in.dp2;
    depDate = in.depDate;
    arrDate = in.arrDate;
    depBody = in.depBody;
    arrBody = in.arrBody;
    
    
    % Julian or Calendar Date
    if isnumeric(depDate)
       depDate = datestr(datetime(depDate,'convertfrom','juliandate')); 
    end
    
    if isnumeric(arrDate)
        arrDate = datestr(datetime(arrDate,'convertfrom','juliandate'));
    end
    
    
    
    % Plotting Constraints
    if isfield(in,'limC3')
        c3lim = in.limC3;
    else
        c3lim = [5:10:100];
    end
    
    if isfield(in,'limVinfDep')
        vinfdeplim = in.limVinfDep;
    else
        vinfdeplim = [3:1:12];
    end
    
    if isfield(in,'limVinfArr')
        vinfarrlim = in.limVinfArr;
    else
        vinfarrlim = [5:1:15];
    end
    
    if isfield(in,'limRLA')
        rlalim = in.limRLA;
    else
        rlalim = [-180:20:180];
    end
    
    if isfield(in,'limDLA')
        dlalim = in.limDLA;
    else
        dlalim = [-35:.5:30];
    end
    
    

    % Figure Plotting
    f = figure();
    set(gcf,'color','w'); box on;
    hold on
    
    
    % Time of Flight
    [C0, H0] = contour(dp1, dp2, tof    , ntof, 'k','ShowText','on');
    clabel(C0,H0,'FontSize',10,'Color','k','FontName','Courier')
    lgnd = {'TOF (days)'}; i=1;
    
    
    % Plot C3
    if contains("C3",in.plts)
        [C1, H1] = contour(dp1, dp2, c3dep  , c3lim, 'r','ShowText','on');
        clabel(C1,H1,'FontSize',10,'Color','r','FontName','Courier');
        lgnd{1,i+1} = '$C_3$ ($km^2/s^2$)'; i=i+1;
    end
    
    % Plot Departure Vinfinity
    if contains("VinfDep",in.plts)
        [C5, H5] = contour(dp1, dp2, vinfdep  , vinfdeplim, 'r','ShowText','on');
        clabel(C5,H5,'FontSize',10,'Color','r','FontName','Courier');
        lgnd{1,i+1} = 'Departure $V_\infty$ ($km/s$)'; i=i+1;
    end
    
    % Plot Arrival Vinfinity
    if contains("VinfArr",in.plts)
        [C2, H2] = contour(dp1, dp2, vinfarr, vinfarrlim, 'b','ShowText','on');
        clabel(C2,H2,'FontSize',10,'Color','b','FontName','Courier')
        lgnd{1,i+1} = 'Arrival $V_\infty$ ($km/s$)'; i=i+1;
    end
    
    % Plot RLA
    if contains("RLA",in.plts)
        [C3, H3] = contour(dp1, dp2, RLA, rlalim, 'Color',[0.3010 0.7450 0.9330],'ShowText','on');
        clabel(C3,H3,'FontSize',10,'Color',[0.3010 0.7450 0.9330],'FontName','Courier')
        lgnd{1,i+1} = 'Departure RLA ($deg.$)'; i=i+1;
    end
    
    % Plot DLA
    if contains("DLA",in.plts)
        [C4, H4] = contour(dp1, dp2, DLA, dlalim, 'Color',[0.4660 0.6740 0.1880],'ShowText','on');
        clabel(C4,H4,'FontSize',10,'Color',[0.4660 0.6740 0.1880],'FontName','Courier')
        lgnd{1,i+1} = 'Departure DLA ($deg.$)'; i=i+1;
    end
    
    
        
    
    
    hold off
    
    
    
    % Plot Labels/Data
    legend(lgnd,'fontsize',12,'interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    xlabel(['Days Past: ', depDate],'fontsize',12, 'interpreter','Latex');
    ylabel(['Days Past: ', arrDate],'fontsize',12, 'interpreter','Latex');
    title(['Type I/II Transfer Between Bodies: ',num2str(depBody), ...
        ' and ', num2str(arrBody)],'fontsize',14, 'interpreter','Latex');

end

