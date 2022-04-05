function broadsearch_plot(out)
%BROADSEARCH_PLOT Plots Results from broadsearch()
%
%   Assumptions/Warnings:
%   	[NOTE: This function is not meant to be standalone]
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. out [struct] Results from broadsearch()
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Dependencies:
%       1. broadsearch_plotTraj()
%       2. broadsearch()
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   TODO:
%       1. Single trajectory plotting
%

    if out.error == false
        gaseqname = broadsearch_sequencename(out.seq);
        Jmax = max(out.costsum);


        figure()
        hold on
        for i=1:height(out.encs)
            broadsearch_plotTraj(out.seq,out.encs(i,:),out.mu,out.pcd,'r');
        end
        hold off
        title([gaseqname, ' Trajectory'],'fontsize',12,'interpreter','latex');
        xlabel('X EC','fontsize',12,'interpreter','latex');
        ylabel('Y EC','fontsize',12,'interpreter','latex');
        zlabel('Z EC','fontsize',12,'interpreter','latex');
        set(gca,'TickLabelInterpreter','latex');
        grid on; box on;
        set(gcf,'color','w'); axis equal;


        if ~isnan(out.cost)
            figure()
            hold on
            scatter([1:height(out.cost)],out.costsum,'r','filled');
            lgnd = {'Total'};

            for i=1:width(out.cost)
                scatter([1:height(out.cost)],out.cost(:,i),'x');
                lgnd{1,i+1} = [gaseqname(i+1),'GA']; %[getPlanetName(out.seq(i)),' GA'];
            end
            hold off
            if max(out.costsum)+1 < Jmax
                ylim([-0.0 max(out.costsum)+1]);
            else
                ylim([-0.0 Jmax]);
            end
            xticks([1:1:height(out.cost)]);
            title([gaseqname, ' Sequence Cost'],'fontsize',12,'interpreter','latex');
            xlabel('Trajectory Number','fontsize',12,'interpreter','latex');
            ylabel('Velocity Discontinuity at Flyby (Cost) (km/s)','fontsize',12,'interpreter','latex');
            legend(lgnd,'fontsize',12,'interpreter','latex');
            set(gca,'TickLabelInterpreter','latex');
            grid on; box on;
            set(gcf,'color','w');
        end
    elseif out.forceplot == true
        % No Complete Trajectory, Plotting Initial and Final leg Points
        
        lps = fieldnames(out.dataTrajUF);
        
        
        figure()
        hold on
        for i=1:height(out.seq)-1
            leg = out.dataTrajUF.(lps{i}).legdata;
            ti = leg(1,:);
            tf = leg(end,:);
            
            broadsearch_plotTraj([out.seq(i);out.seq(i+1)],[ti(2); ti(3)],out.mu,out.pcd,'r')
            broadsearch_plotTraj([out.seq(i);out.seq(i+1)],[tf(2); tf(3)],out.mu,out.pcd,'b')
     
        end
        hold off
        grid on; box on; axis equal;
        
        
        
    else
        disp('No Plots Available. Broad Search Error')
    end

end

























