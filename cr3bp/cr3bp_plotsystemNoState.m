function cr3bp_plotsystemNoState(varargin)
%CR3BP_PLOTSYSTEMNOSTATE Plot N.D. System Only (without trajectory)
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. sys      [struct]    Structure Containing System L,V,T,mu vals.
%       4. optional
%               J               [1x1] Jacobi Constant Value for ZVC
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: (none)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    if length(varargin)==1
        sys = varargin{1};
        plotZVC = false;
    elseif length(varargin)==2
        sys = varargin{1};
        J   = varargin{2};
        plotZVC = true;
    end
    

    tiledFigTotal = 1;


    figure()
    ft = tiledlayout(tiledFigTotal,1,'Tilespacing','compact','Padding','compact');
    title(ft, [sys.name,' System'], ...
        'fontsize',14,'interpreter','latex');
    set(gcf,'color','w');

    nexttile
    lgndNum = 0;
    hold on
    

    % Plots Primary and Secondary Bodies
    scatter3(sys.b1(1),sys.b1(2),sys.b1(3),100,[0.4660 0.6740 0.1880],'filled');
    scatter3(sys.b2(1),sys.b2(2),sys.b2(3),'k','filled');
    lgnd{1,lgndNum+1} = 'Primary Body';
    lgnd{1,lgndNum+2} = 'Secondary Body';
    lgndNum = lgndNum+2;


    % Plots Zero Velocity Curves for Initial State Jacobi Constant
    if plotZVC
        title(ft, [sys.name,' System  |  Jacobi Constant (C) = ', num2str(J,8)], ...
        'fontsize',14,'interpreter','latex');
        zvs = cr3bp_jacobiZVC(sys,J,'onlyplotcurve',true,'simplePlot',true);
        lgnd{1,lgndNum+1} = 'ZVC';
        lgndNum = lgndNum + 1;
    end


    % Plots Lagrange Points of the System
    if true
        sys = cr3bp_computeLagrangePoints(sys);
        Lpts  = [sys.L1.';
            sys.L2.';
            sys.L3.';
            sys.L4.';
            sys.L5.'];
        disp('hi')
        hold on
        scatter(Lpts(:,1),Lpts(:,2),'d','r','filled');
        hold off
        lgnd{1,lgndNum+1} = 'Lagrange Points';
        lgndNum = lgndNum + 1;
    end
    
    hold off
    grid on; box on; axis equal;
    legend(lgnd,'fontsize',12,'interpreter','latex');
    xlabel('X [ND]','fontsize',12,'interpreter','Latex');
    ylabel('Y [ND]','fontsize',12,'interpreter','Latex');
    zlabel('Z [ND]','fontsize',12,'interpreter','Latex');
    
end