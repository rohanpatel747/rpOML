function dYdt = eomNBP(t,Y, mu, pList, muList, eti, cbname, backOffset)
%EOMNBP N-Body EOM [6x1] State Space of the Rates for Integration
%
%   *** WORK IN PROGRESS, CODE NOT UNIT TESTED YET ***
%
%
%   Assumptions/Warnings:
%   	1. First perturbing body in pList is the departure body and so only
%   	it is affected by the "backOffset"
%
%       2. BackOffset is used to remove a singularity at dt=0 when the
%       departure body and the the satellite have a cooresponding position
%       vector.
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. t          [1x1 Int.]  dt time step from ODE45
%       2. Y          [6x1 Dbl.]  State to pass into integ. instant
%       3. mu         [1x1 Dbl.]  Central Body Gravitational Param. (km3/s2)
%       4. pList      [nx1 Int.]  Perturbing Bodies List
%       5. muList     [nx1 Dbl.]  Perturbing Bodies Grav. Param. (km3/s2)
%       6. eti        [1x1 Dbl.]  Ephemeris Time of Initial Epoch (s)
%       7. cbname     [    Str.]  Central Body SPICE Name
%       8. backOffset [1x1 Dbl.]  Offset of Initial Body State by #Days (s)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. dydt       [6x1 Dbl.]  State Space Vector for Integration
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Dependencies:
%       1. mice  - mice_spkezr()


    % Spacecraft State
    r_ = Y(1:3);   r  = norm(r_);
    v_ = Y(4:6);

    % Gravitational Force from Additional Bodies
    for i=1:length(pList)
        
        % Body Position Vector
        if i==1
            %b_ = mice_spkezr(num2str(pList(i)), eti+t-backOffset, 'J2000', ...
            %'NONE', cbname).state(1:3);
            b_ = getStatePlanet(pList(i),(eti + (t/86400))).x.';
        else
           %b_ = mice_spkezr(num2str(pList(i)), eti+t, 'J2000', ...
           % 'NONE', cbname).state(1:3);
           b_ = getStatePlanet(pList(i),(eti + (t/86400))).x.';
        end
        b  = b_/(norm(b_)^3);
        
        % Relative Position of S/C to Body
        bb_= b_ - r_;
        bb = bb_/(norm(bb_)^3);
        
        % Force Contribution
        fb(i,1:3) = (muList(i)*(b - bb)).';
        
    end
    
    % DE EOM
    Dx   = v_;
    D2x  = (-mu/r^3)*r_;
    for i=1:length(pList)
       D2x = D2x - fb(i,1:3).';
    end

    dYdt = [Dx; D2x];
end