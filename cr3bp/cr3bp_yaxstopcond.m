function [value,isterminal,direction] = cr3bp_yaxstopcond(~,x)

    value = x(5);
    isterminal = 1;
    direction = 0;

    value = [x(2), x(5)<0];
    isterminal = [1, 0];
    direction = [1 1];   
end