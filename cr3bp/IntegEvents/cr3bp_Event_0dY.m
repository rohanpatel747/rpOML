function [value,isterminal,direction] = cr3bp_Event_0dY(~,x)
%CR3BP_EVENT_0dY Terminate ODE113 at first dy=0.
    value = x(5);
    isterminal = 1;
    direction = 0;
end

