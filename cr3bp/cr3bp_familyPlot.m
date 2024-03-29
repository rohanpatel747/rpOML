function out = cr3bp_familyPlot(c3sys, A, sortbyIndex, colorType, pltFig, initialPlotNum, finalPlotNum, stepPlotNum, opts)
%CR3BP_FAMILYPLOT Plot Periodic Orbit Families from: @CR3BP_FAMILYNPCSTATEX 
%
%   Assumptions/Warnings: 
%       *** Currently only tested for L1/L2 Lyaponov Orbits! ***
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. c3sys      [struct  ]  Structure Containing at least Sys. mu val
%       2. A          [nx9 Dbl.]  Array of all PO families merged row-wise
%       3. sortbyIndex[1x1 Dbl.]  Row to Color Data By:
%                                   (9)    Iterations Required
%                                   (8)    Jacobi Constant
%                                   (7)    Orbit Period
%                                   (6->1) State Vars. [dZ dY dX Z Y X]
%       4. colorType  [str.]      MATLAB Colormap Type:
%                                   'parula','cool','jet', etc.
%       5. pltFig     [T/F]       Plot Figure (or just compute data)?
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. A          [nx10 Dbl.] Array of all PO families merged row-wise
%                                 Including the PO color index (row 1)
%       2. colordata  [nx1  Dbl.] RGB color pair w.r.t A(:,1)
%       3. colorname  [str.]      String of MATLAB Colormap used.
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
    
    arguments
        c3sys;
        A;
        sortbyIndex;
        colorType;
        pltFig;
        initialPlotNum;
        finalPlotNum;
        stepPlotNum;
        opts.plotDim = false;
    end

    options = odeset('reltol',1e-13,'abstol',1e-13);

    
    mu = c3sys.mu;


    if opts.plotDim
        if sortbyIndex == 7
            A(:,7) = A(:,7) .* (c3sys.T / (2*pi));
            A(:,7) = A(:,7) ./ 86400;
        end
        if sortbyIndex == 0
            A(:,end+1) = (c3sys.b2(1) - A(:,1)).*c3sys.L;
        end
    end

    if opts.plotDim
        if sortbyIndex == 0
            minval = min(A(:,end));
            maxval = max(A(:,end));
        else
            minval = min(A(:,sortbyIndex));
            maxval = max(A(:,sortbyIndex));
        end
    else
        minval = min(A(:,sortbyIndex));
        maxval = max(A(:,sortbyIndex));
    end
    

    
    % Normalize Sort Index Values to 0-1 Scale
    if opts.plotDim
        if sortbyIndex == 0
            B = A(:,end);
        else
            B = A(:,7);
        end
    else
        B = A(:,sortbyIndex);
    end

    B = B-minval;
    B = (1/max(B)).*B;
    
    % Initialize Colormap and Normalize Colormap RGB to 0-1 Scale
    cmap     = colormap(colorType);
    cmap     = [((1:length(cmap)).'./length(cmap)) , cmap];
    ncolors  = length(cmap);
    clrIdx   = ((1/ncolors).*(1:ncolors)).';
    
    % Match ClrIdx values to B matrix Values
    d        = abs(bsxfun(@minus, B(:).', clrIdx(:))); 
    [~, ind] = min(d);
    clrIdxA  = clrIdx(ind);
    
    % Join B with A and Find Associated RGB Triplet (stored in colorIndex) 
    A = [clrIdxA, A];
    for i=1:length(cmap)
    
        ci = cmap(i);
        ai = find(A(:,1)==ci);
        for j=1:length(ai)
            ind = ai(j);
            %colorIndex(ind,1:3) = cmap(ncolors+1-i,2:4);\
            colorIndex(ind,1:3) = cmap(i,2:4);
        end
    
    end
    
    % Color Bar Legend Variable Names
    if     sortbyIndex == 8
        cbarName = 'Jacobi Constant (C)';
    elseif sortbyIndex == 9
        cbarName = 'Number of Iterations (itr)';
    elseif sortbyIndex == 7
        cbarName = 'Orbit Period (T) (ND)';
    elseif sortbyIndex == 6
        cbarName = 'Initial State (dZ) (ND)';
    elseif sortbyIndex == 5
        cbarName = 'Initial State (dY) (ND)';
    elseif sortbyIndex == 4
        cbarName = 'Initial State (dX) (ND)';
    elseif sortbyIndex == 3
        cbarName = 'Initial State (Z) (ND)';
    elseif sortbyIndex == 2
        cbarName = 'Initial State (Y) (ND)';
    elseif sortbyIndex == 1
        cbarName = 'Initial State (X) (ND)';
    elseif sortbyIndex == 0
        cbarName = 'Amplitude w.r.t Moon (km)';
    end
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Create Figure
    if pltFig
        if ~opts.plotDim
            ft = tiledlayout(1,1,'Tilespacing','compact','Padding','compact');
            %title(ft, 'Lyapunov Periodic Orbit Family | Earth-Moon $L_1$', ...
            %    'fontsize',14,'interpreter','latex');
            
            hold on
            scatter3(c3sys.b1(1),c3sys.b1(2),c3sys.b1(3),'b','filled');
            scatter3(c3sys.b2(1),c3sys.b2(2),c3sys.b2(3),'k','filled');
            scatter3(c3sys.L1(1),c3sys.L1(2),c3sys.L1(3),'d','r','filled');
            scatter3(c3sys.L2(1),c3sys.L2(2),c3sys.L2(3),'d','r','filled');
            scatter3(c3sys.L3(1),c3sys.L3(2),c3sys.L3(3),'d','r','filled');
            scatter3(c3sys.L4(1),c3sys.L4(2),c3sys.L4(3),'d','r','filled');
            scatter3(c3sys.L5(1),c3sys.L5(2),c3sys.L5(3),'d','r','filled');
            
            
            for i=initialPlotNum:stepPlotNum:finalPlotNum
                [~ ,xx] = ode113(@(t,Y) eomCR3BP(t,Y,mu), [0 A(i,8)], A(i,2:7).', options);
                plot3( xx(:,1), xx(:,2), xx(:,3),'linewidth',1.5,'color',colorIndex(i,1:3));
            end
            hold off 
            
            cb = colorbar;
            caxis([minval maxval]);
            ylabel(cb,cbarName,'fontsize',16,'interpreter','Latex'); 
            
            grid on; box on; axis equal; set(gcf,'color','w');
            xlabel('X (ND)','fontsize',16,'interpreter','Latex');
            ylabel('Y (ND)','fontsize',16,'interpreter','Latex');
            zlabel('Z (ND)','fontsize',16,'interpreter','Latex');
        else
            ft = tiledlayout(1,1,'Tilespacing','compact','Padding','compact');
            b1 = c3sys.b1(1:3).*c3sys.L;
            b2 = c3sys.b2(1:3).*c3sys.L;
            L1 = c3sys.L1(1:3).*c3sys.L;
            L2 = c3sys.L2(1:3).*c3sys.L;
            L3 = c3sys.L3(1:3).*c3sys.L;
            L4 = c3sys.L4(1:3).*c3sys.L;
            L5 = c3sys.L5(1:3).*c3sys.L;

            hold on;
            scatter3(b1(1),b1(2),b1(3),'b','filled');
            scatter3(b2(1),b2(2),b2(3),'k','filled');
            scatter3(L1(1),L1(2),L1(3),'d','r','filled');
            scatter3(L2(1),L2(2),L2(3),'d','r','filled');
            scatter3(L3(1),L3(2),L3(3),'d','r','filled');
            scatter3(L4(1),L4(2),L4(3),'d','r','filled');
            scatter3(L5(1),L5(2),L5(3),'d','r','filled');
            
            for i=initialPlotNum:stepPlotNum:finalPlotNum
                [~ ,xx] = ode113(@(t,Y) eomCR3BP(t,Y,mu), [0 A(i,8)], A(i,2:7).', options);
                plot3( xx(:,1).*c3sys.L, xx(:,2).*c3sys.L, xx(:,3).*c3sys.L,'linewidth',1.5,'color',colorIndex(i,1:3));
            end
            hold off;

            cb = colorbar;
            caxis([minval maxval]);
            ylabel(cb,cbarName,'fontsize',16,'interpreter','Latex'); 
            
            grid on; box on; axis equal; set(gcf,'color','w');
            xlabel('X (km)','fontsize',16,'interpreter','Latex');
            ylabel('Y (km)','fontsize',16,'interpreter','Latex');
            zlabel('Z (km)','fontsize',16,'interpreter','Latex');

        end
    end
    
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Outputs
    out           = struct;
    out.data      = A;
    out.colordata = colorIndex;
    out.colorname = colorType;

end
