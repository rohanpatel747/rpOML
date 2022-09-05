function [value,isterminal,direction] = cr3bp_yaxstopcond(~,x)
%CR3BP_YAXSTOPCOND Terminate ODE113 at first y=0 (OLD VERSION - NO USE).
    value = x(2);
    isterminal = 1;
    direction = 0;

end