function [dt, stateX3] = plotLambertTransfer(mu, x1_, x2_, dv1_, dt)

    x1_t = x1_;
    x1_t(4) = x1_(4) + dv1_(1);
    x1_t(5) = x1_(5) + dv1_(2);
    x1_t(6) = x1_(6) + dv1_(3);
    

    fX1 = create_state(x1_, mu, 'rv', false, false, true);          % Initial
    fX2 = create_state(x2_, mu, 'rv', false, false, true);          % Final
    fX3 = create_state(x1_t, mu, 'rv', false, false, true);         % Transfer
    [~,stateX1] = plotTraj(fX1, [0 dt], 1000, 'none');              % Integ. Initial  [0 dt]
    [~,stateX2] = plotTraj(fX2, [dt 0], 1000, 'none');              % Integ. Final    [dt 0]
    [~,stateX11] = plotTraj(fX1, [0 fX1.T], 1000, 'none');          % Integ. Initial  [0 T]
    [~,stateX22] = plotTraj(fX2, [fX2.T 0], 1000, 'none');          % Integ. Final    [0 T]
    [~,stateX3] = plotTraj(fX3, [0 dt], 1000, 'none');              % Integ. Transfer [0 dt]
    [~,stateX33] = plotTraj(fX3, [0 fX3.T], 1000, 'none');          % Integ. Transfer [0 T]

    
    %figure
    hold on
    % Original/Final Orbit
    plot3(stateX11(:,1),stateX11(:,2),stateX11(:,3),'-.k')
    plot3(stateX22(:,1),stateX22(:,2),stateX22(:,3),'-.k')
    plot3(stateX33(:,1),stateX33(:,2),stateX33(:,3),'-.k')

    % Original/Final Traveled Orbit
    plot3(stateX1(:,1),stateX1(:,2),stateX1(:,3),'b','linewidth',1.5)
    plot3(stateX2(:,1),stateX2(:,2),stateX2(:,3),'r','linewidth',1.5)

    % r1_ and r2_
    quiver3(0,0,0,x1_(1),x1_(2),x1_(3),'autoscale','off','color','b')
    quiver3(0,0,0,x2_(1),x2_(2),x2_(3),'autoscale','off','color','r')
    
    %quiver3(0,0,0,fX3.e_(1)*1e9,fX3.e_(2)*1e9,fX3.e_(3)*1e9,'autoscale','off','color','g')
    %quiver3(0,0,0,-fX3.e_(1)*1e9,-fX3.e_(2)*1e9,-fX3.e_(3)*1e9,'autoscale','off','color','g')

    % Transfer Orbit
    plot3(stateX3(:,1),stateX3(:,2),stateX3(:,3),'g','linewidth',1.5)

    hold off
    grid on; box on; axis equal;
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    set(gcf,'color','w');

end
