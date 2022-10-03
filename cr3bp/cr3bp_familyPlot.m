function out = cr3bp_familyPlot(c3sys, A, sortbyIndex, colorType, pltFig, initialPlotNum, finalPlotNum, stepPlotNum)
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
    
    mu = c3sys.mu;
    
    minval = min(A(:,sortbyIndex));
    maxval = max(A(:,sortbyIndex));
    
    % Normalize Sort Index Values to 0-1 Scale
    B = A(:,sortbyIndex);
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
        cbarName = 'Orbit Period (T) (ND Rot.)';
    elseif sortbyIndex == 6
        cbarName = 'Initial State (dZ) (ND Rot.)';
    elseif sortbyIndex == 5
        cbarName = 'Initial State (dY) (ND Rot.)';
    elseif sortbyIndex == 4
        cbarName = 'Initial State (dX) (ND Rot.)';
    elseif sortbyIndex == 3
        cbarName = 'Initial State (Z) (ND Rot.)';
    elseif sortbyIndex == 2
        cbarName = 'Initial State (Y) (ND Rot.)';
    elseif sortbyIndex == 1
        cbarName = 'Initial State (X) (ND Rot.)';
    end
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Create Figure
    if pltFig
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
        
        options = odeset('reltol',1e-13,'abstol',1e-13);
        for i=initialPlotNum:stepPlotNum:finalPlotNum
            [~ ,xx] = ode113(@(t,Y) eomCR3BP(t,Y,mu), [0 A(i,8)], A(i,2:7).', options);
            plot3( xx(:,1), xx(:,2), xx(:,3),'linewidth',1.5,'color',colorIndex(i,1:3));
        end
        hold off 
        
        cb = colorbar;
        caxis([minval maxval]);
        ylabel(cb,cbarName,'fontsize',12,'interpreter','Latex'); 
        
        grid on; box on; axis equal; set(gcf,'color','w');
        xlabel('X Rot. (ND)','fontsize',12,'interpreter','Latex');
        ylabel('Y Rot. (ND)','fontsize',12,'interpreter','Latex');
        zlabel('Z Rot. (ND)','fontsize',12,'interpreter','Latex');
        legend({'Primary','Secondary','Lagrange Pts.'}, ...
         'fontsize',12,'interpreter','Latex','location','southwest');
    end
    
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Outputs
    out           = struct;
    out.data      = A;
    out.colordata = colorIndex;
    out.colorname = colorType;

end
