function [value,isterminal,direction] = cr3bp_Event_xdesired(~,x,xdesired)
%CR3BP_EVENT_0Y Terminate ODE113 when x-state = desired x-state value.
    value = x(1) - xdesired;
    isterminal = 1;
    direction = 0;
end

