function f = pkchp_plt(in)
%PKCHP_PLT Summary of this function goes here
%   Detailed explanation goes here

    c3dep   = in.c3dep;
    vinfarr = in.vinfarr;
    tof     = in.tof;
    ctr1    = in.ctr1;
    ctr2    = in.ctr2;
    ntof    = in.ntof;
    dp1     = in.dp1;
    dp2     = in.dp2;
    depDate = in.depDate;
    arrDate = in.arrDate;
    depBody = in.depBody;
    arrBody = in.arrBody;


    f = figure();
    set(gcf,'color','w'); box on;
    % Plot Contours of C3, Vinf, and TOF
    hold on
    [C1, H1] = contour(dp1, dp2, c3dep  , ctr1, 'r','ShowText','on');
    [C2, H2] = contour(dp1, dp2, vinfarr, ctr2, 'b','ShowText','on');
    [C3, H3] = contour(dp1, dp2, tof    , ntof, 'k','ShowText','on');
    hold off
    % Contour Line Text Size/Color
    clabel(C1,H1,'FontSize',10,'Color','r','FontName','Courier')
    clabel(C2,H2,'FontSize',10,'Color','b','FontName','Courier')
    clabel(C3,H3,'FontSize',10,'Color','k','FontName','Courier')
    % Plot Labels/Data
    legend({'$C_3$ ($km^2/s^2$)','Arrival $V_\infty$ ($km/s$)', 'TOF (days)'}, ...
        'fontsize',12,'interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    xlabel(['Days Past: ', depDate],'fontsize',12, 'interpreter','Latex');
    ylabel(['Days Past: ', arrDate],'fontsize',12, 'interpreter','Latex');
    title(['Type I/II Transfer Between Bodies: ',num2str(depBody), ...
        ' and ', num2str(arrBody)],'fontsize',14, 'interpreter','Latex');

end

