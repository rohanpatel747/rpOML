function out = getStatePlanet(ID, JDE, method)
%GETSTATEPLANET Returns [6x1] Cartesian State for given Julian Date(s) Vector
%
%   Assumptions/Warnings:
%   	1. If method is not specified, algorithm will use Meeus Method
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. ID  [1x1 Int.]  Planet ID Number (3=Earth, 1-9 inc. Pluto)
%       2. JDE [nx1 Dbl.]  Julian Date Row or Col. Vector
%       3. method   [str]  (Optional) Computation Method
%                   ***DEFAULT USES MEEUS METHOD***
%                   'meeus' = Meeus Method w/ Planetary Coeffs.
%                   'de'    = Query DE Ephemerides File (DE430.bsp)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. t [nx1] Julian Date
%       2. x [nx6] Inertial State Vector [x,y,z,vx,vy,vz] km and km/s
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Dependencies:
%       1. rpOML - getMeeusData()
%       2. mice  - cspcie_str2et()
%       3. mice  - mice_spkezr()
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Example:
%       out = getPlanetState(3, 2455450)
%       out = getPlanetState(3, [2455450; 2455610], 'de')
%

arguments
   ID
   JDE
   method = 'meeus'; % Use Meeus Method by Default
end

    % NAIF ID PLANETARY DATA FOR PLANETARY STATES
    if contains('de',method)
        
        % Condition Julian Dates for MICE
        % JDE(:).'
        for j=1:length(JDE)
           jd = JDE(j);
           jds(j,:) = ['JD',num2str(jd)];
        end
        et = cspice_str2et(jds);
        
        % Get State from DE Ephemeris Files
        JDE = JDE(:).';
        states = mice_spkezr(num2str(ID), et, 'J2000', 'NONE', 'SOLAR SYSTEM BARYCENTER');
        
        % Output Data
        out = struct;
        for j=1:length(JDE)
           out.t(j,1) = JDE(j);
           out.x(j,:) = states(j).state.';
        end


    % MEEUS METHOD FOR PLANETARY STATE
    elseif contains('meeus',method)
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Constants

        % Mu Sun
        mu = 1.32712440018e11;  % Taken from Constants planetaryconstants()

        % AUKM
        aukm = 1.49597870700e8; % Taken from Constants planetaryconstants()

        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Loop Through Julian Dates Given
        out = struct;
        
        for j=1:length(JDE)
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Time
            T = (JDE(j) - 2451545.0)/36525;

            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Get Planet Data
            oe = getMeeusData(ID);
            E  = zeros(6,1);
            for i=1:6
                E(i,:) = oe(i,1) + oe(i,2)*T + oe(i,3)*(T^2) + oe(i,4)*(T^3);
            end

            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Variables
            L = E(1);
            a = E(2);
            e = E(3);
            i = E(4);
            O = E(5);
            II= E(6);

            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Mean Anomaly
            M = L - II;
            M = M*(pi/180);

            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % True Anomaly
            Ccen = ((2*e)         - ((e^3)/4) + ((5/96)*(e^5)))  * sin(  M)  + ...
                   ((5/4)*(e^2)   - (11/24)*(e^4))               * sin(2*M)  + ...
                   ((13/12)*(e^3) - (43/64)*(e^5))               * sin(3*M)  + ...
                   ((103/96)*(e^4))                              * sin(4*M)  + ...
                   ((1097/960)*(e^5))                            * sin(5*M);
            ta = M+Ccen;
            ta = ta*(180/pi);
            
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Argument of Perihelion
            w = II-O;
            
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Output
            st2ele     = conv_ele2state([a*aukm,e,i,O,w,ta], mu, true, false);
            out.x(j,:) = st2ele.xi;
            out.t(j,1) = T;
        end
    end
    
end