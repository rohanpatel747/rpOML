function [value,isterminal,direction] = cr3bp_Event_0Y(~,x)
%CR3BP_EVENT_0Y Terminate ODE113 at first y=0 occurance.
    value = x(2);
    isterminal = 1;
    direction = 0;
end

