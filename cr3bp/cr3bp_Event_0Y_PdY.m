function [value,isterminal,direction] = cr3bp_Event_0Y_PdY(~,x)
%CR3BP_EVENT_0Y_PdY Terminate ODE113 at first y=0 and dy>0.
    value = [x(2), x(5)<0];
    isterminal = [1, 0];
    direction = [1 1];
end