function [value,isterminal,direction] = cr3bp_yaxstopcond(~,x)
    value = x(2);
    isterminal = 1;
    direction = 0;   
end