function figureDefaults(axlbl, opts)
%FIGUREDEFAULTS
%
%   'figType' array with 3 elements:
%       1 = X/Y/Z Rot. (ND)             figureDefaults(ax1, "figType",  1);
%      -1 = X/Y/Z Inertial (km)         figureDefaults(ax1, "figType", -1);
%
%
%
%
%
%
%
%
%
%
%
%
%

    
    
    arguments
        axlbl;
        opts.figType       = 0;
        opts.axFontSize    = 16;
        opts.titleFontSize = 20;
        opts.legendFontSize= 16;
        opts.length        = 800;
        opts.height        = 800;
        opts.figColor      = 'normal';
    end
    
    
    set(gcf,'color','w','Position',[100, 100, 100+opts.length, 100+opts.height]);
    
    grid on; box on;
    axlbl.TickLabelInterpreter = 'latex';
    axlbl.FontSize = opts.axFontSize;
    
    
    axlbl.XLabel.Interpreter = 'latex';
    axlbl.YLabel.Interpreter = 'latex';
    axlbl.ZLabel.Interpreter = 'latex';
    axlbl.XLabel.FontSize = opts.axFontSize;
    axlbl.YLabel.FontSize = opts.axFontSize;
    axlbl.ZLabel.FontSize = opts.axFontSize;

    axlbl.Title.Interpreter = 'latex';

    % Plot 3-D Axes
    if     opts.figType ==  1       
        axlbl.XLabel.String = 'X (ND)';
        axlbl.YLabel.String = 'Y (ND)';
        axlbl.ZLabel.String = 'Z (ND)';
    elseif opts.figType == -1
        axlbl.XLabel.String = 'X Inertial (km)';
        axlbl.YLabel.String = 'Y Inertial (km)';
        axlbl.ZLabel.String = 'Z Inertial (km)';
    elseif opts.figType == -2
        set(gca, 'YDir','reverse')
        axlbl.XLabel.String = 'B$\cdot$T (km)';
        axlbl.YLabel.String = 'B$\cdot$R (km)';
    end
    axis equal;
    

    if contains(opts.figColor, 'invert')
        set(gcf(),'Color','k');
        set(axlbl.XLabel,'Color','w');
        set(axlbl.YLabel,'Color','w');
        set(axlbl.ZLabel,'Color','w');
        set(axlbl.Title, 'Color','w');
        set(axlbl,'Color','w','GridColor','w','MinorGridColor','w',...
            'XColor','w','YColor','w','ZColor','w');
        set(axlbl,'Color', [0 0 0])
    end
    
end