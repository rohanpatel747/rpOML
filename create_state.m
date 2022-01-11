function out = create_state(varargin)
%CREATE_STATE Create a Full State Structure given a partial state
%   Partial State = Either the R&V vectors, OR the a,e,i,o,w,Ta elements
%                   (in degrees)
%
%   _______________________________________________________________________
%   Required Inputs:
%       xi   = [r;v] or [a;e;i;o;w;ta] (ta [0 360] degrees) or full state Structure;
%       mu   = Gravitational parameter of central body (km3/s2)
%       type = 'rv' or 'kep' or 'fullX'
%   Optional Inputs:
%       displayOut = Print Output (Default: FALSE), [true/false]
%       isDeg      = i,o,w in Degrees (Default: TRUE), [true/false]
%       addDefs    = Add Variable Definitions (Default: FALSE), [true/false]
%
%
%   Outputs: 
%       out = Structure Containing Full State
%
%   _______________________________________________________________________
%   Use Cases:
%   (1) Find full state given a state vector
%       The x (6x1) vector is known [r;v] and the full state is desired.
%       We also want the to know what each variable means and their units.
%       
%       r1 = [-720000;  670000; 310000];
%       v1 = [ 2.160 ; -3.360 ; 0.620 ];
%       x1 = [r1;v1];
%       displayOut = false;
%       isDeg = false;
%       addDefs = true;
%       fX1 = create_state(x1, mu_s, 'rv', displayOut, isDeg, addDefs);
%
%
%
%   (2) Find full state given Keplerian Orbital Elements
%       A, e, i, o, w, and ta are known and angles are in degrees.
%       
%       aeiowta = [a; e; i; o; w; ta];
%       fX1 = create_state(aeiowta, mu_s, 'kep') 
%
%
%
%   (3) Modify/Recalculate State w/ Changed Parameters:
%       A produced full state structure can be modified, and refed into
%       this function to produce a new full state structure. For example:
%
%       fX1 = create_state(x1, mu_s, 'rv');      % State at some initial x1
%       fX2 = fX1;                               % Make a copy to modify
%       fX2.ta = -13.1687;                       % Change True Anomaly
%       fX2.changed = 'aeiouta';                 % Tell Function what you have changed (either 'aeiowta' or 'rv')
%       fX2 = create_state(fX2, mu_s, 'fullX');  % Recalculate full state
%
%   _______________________________________________________________________
%   Definitions and Units
%        xi        : Inertial State Vector                  (km, km/s)
%        xp        : Perifocal State Vector                 (km, km/s)
%        R313      : Perifocal to Inertial Rotation DCM     (rad)
%        mu        : Central Body Gravitaitonal Parameter   (km3/s2)
%        r_        : Inertial Position Vector               (km)
%        v_        : Inertial Velocity Vector               (km/s)
%        h_        : Angular Momentum Vector                (km2/s)
%        h         : Angular Momentum Magnitude             (km2/s)
%        a         : Semi-major axis                        (km)
%        e_        : Eccentricity Vector                    (n/a)
%        e         : Eccentricity                           (n/a)
%        i         : Inclination                            (deg)
%        o         : RAAN                                   (deg)
%        w         : Argument of Periapsis                  (deg)
%        ta        : True Anomaly                           (deg)
%        ma        : Mean Anomaly                           (deg)
%        ea        : Eccentric Anomaly                      (deg)
%        p         : Semi-latus Rectum                      (km)
%        rp        : Radius at Periapsis                    (km)
%        ra        : Radius at Apoapsis                     (km)
%        T         : Period                                 (s)
%        dt        : Time since Periapsis                   (s)
%        nodesVec  : Node Vector                            (km)
%        E         : energy                                 (km2/s2)
%        orbitType : Orbit Direction                        (pro/retro)
%
%

    % _____________________________________________________________________
    % Inputs Filtering 
    switch nargin
        case 3
            if class(varargin{1,1}) == 'struct'
               fX = varargin{1,1};
            else
                xi = cell2mat(varargin(1));
            end
            mu = cell2mat(varargin(2));
            type = cell2mat(varargin(3));
            displayOut = false;
            isDeg = true;
            addDefs = false;
        case 4
            if class(varargin{1,1}) == 'struct'
               fX = varargin{1,1}; 
            else
                xi = cell2mat(varargin(1));
            end
            mu = cell2mat(varargin(2));
            type = cell2mat(varargin(3));
            displayOut = cell2mat(varargin(4));
            isDeg = true;
            addDefs = false;
        case 5
            if class(varargin{1,1}) == 'struct'
               fX = varargin{1,1}; 
            else
                xi = cell2mat(varargin(1));
            end
            mu = cell2mat(varargin(2));
            type = cell2mat(varargin(3));
            displayOut = cell2mat(varargin(4));
            isDeg = cell2mat(varargin(5));
            addDefs = false;
        case 6
            if class(varargin{1,1}) == 'struct'
               fX = varargin{1,1}; 
            else
                xi = cell2mat(varargin(1));
            end
            mu = cell2mat(varargin(2));
            type = cell2mat(varargin(3));
            displayOut = cell2mat(varargin(4));
            isDeg = cell2mat(varargin(5));
            addDefs = cell2mat(varargin(6));
        otherwise
            error('Unexpected number of inputs. See function help')
    end

    % _____________________________________________________________________
    % Get State
    if strcmp(type, 'rv')
        % Only r and v vector given
        rvele = conv_state2ele(xi, mu, false);
        kpele= conv_ele2state([
            rvele.a;
            rvele.e;
            rvele.i;
            rvele.o;
            rvele.w;
            rvele.ta], mu, isDeg, false);   
    elseif strcmp(type, 'kep')
        % Only Keplerian elements (a,e,i,o,w,ta) given
        kpele = conv_ele2state(xi, mu, isDeg, false);
        rvele = conv_state2ele(kpele.xi, mu, false);
    elseif strcmp(type, 'fullX')
        % Full State Structure Given
        if strcmp(fX.changed, 'aeiowta')
            % Compute New Full State Structure w changes to Keplerian Ele.
            aeiouwta = [
              fX.a;
              fX.e;
              fX.i;
              fX.o;
              fX.w;
              fX.ta;
            ];
            kpele = conv_ele2state(aeiouwta, mu, isDeg, false);
            rvele = conv_state2ele(kpele.xi, mu, false);
        elseif strcmp(fX.changed, 'rv')
            % Compute New Full State Structure w changes to xi vector.
            rvele = conv_state2ele(fX.xi, mu, false);
            kpele= conv_ele2state([
                rvele.a;
                rvele.e;
                rvele.i;
                rvele.o;
                rvele.w;
                rvele.ta], mu, isDeg, false);
        end
    end

    anele = getAnomalyandDt(rvele,false);
    
    % _____________________________________________________________________
    % Output Structure
    out = struct;
    out.xi  = rvele.xi;
    out.xp  = kpele.xp;
    out.R313= kpele.R313;   
    out.mu  = mu;
    out.r_  = rvele.r_;
    out.v_  = rvele.v_;
    out.h_  = rvele.h_;
    out.h   = rvele.h;
    out.a   = rvele.a;
    out.e_  = rvele.e_;
    out.e   = rvele.e;
    out.i   = rvele.i;
    out.o   = rvele.o;
    out.w   = rvele.w;
    out.ta  = rvele.ta;
    out.rp  = rvele.rp;
    out.ra  = rvele.ra;
    out.T   = rvele.T;
    out.orbitType = rvele.orbitType;
    out.N_  = rvele.nodesVec;
    out.E   = rvele.E;
    out.p   = rvele.p;
    out.ea  = anele.ea;
    out.ma  = anele.ma;
    out.dt  = anele.dt;
    disp('All angles in DEGREES')
    
    if addDefs
        out.unitNames = [
            "xi : Inertial State Vector (km, km/s)";
            "xp : Perifocal State Vector (km, km/s)";
            "R313 : Perifocal to Inertial Rotation DCM";
            "mu : Central Body Gravitaitonal Parameter (km3/s2)";
            "r_ : Inertial Position Vector (km)";
            "v_ : Inertial Velocity Vector (km/s)";
            "h_ : Angular Momentum Vector (km2/s)";
            "h : Angular Momentum Magnitude (km2/s)";
            "a : Semi-major axis (km)";
            "e_ : Eccentricity Vector (n/a)";
            "e : Eccentricity";
            "i : Inclination (deg)";
            "o : RAAN (deg)";
            "w : Argument of Periapsis (deg)";
            "ta : True Anomaly (deg)";
            "ma : Mean Anomaly (rad)";
            "ea : Eccentric Anomaly (rad)";
            "p : Semi-latus Rectum (km)";
            "rp : Radius at Periapsis (km)";
            "ra : Radius at Apoapsis (km)";
            "T : Period (s)";
            "dt : Time since Periapsis (s)";
            "nodesVec : Node Vector (km)";
            "E : energy (km2/s2)";
            "orbitType : Orbit Direction (pro/retro)"];
    end 
    
    
    
    % _____________________________________________________________________
    % Display Out
    if round(out.e, 5) == 1;     conicType = 'Parabolic';
    elseif round(out.e, 5) < 1;  conicType = 'Elliptical';
    elseif round(out.e, 5) == 0; conicType = 'Circular';
    else;                        conicType = 'Hyperbolic';
    end
    
    if displayOut
        disp('-------------------------------------------')
        disp('Inertial State Vector: [x;y;z;vx;vy;vz]')
        disp(' ')
        disp(out.xi)
        disp('-------------------------------------------')
        disp('Perifocal State Vector: [x;y;z;vx;vy;vz]')
        disp(' ')
        disp(out.xp)
        disp('-------------------------------------------')
        disp('Keplerian Elements:')
        disp(' ')
        disp(['(a)  Semimajor Axis   : ',num2str(out.a),' km'])
        disp(['(e)  Eccentricity     : ',num2str(out.e),' (n/a)'])
        disp(['(i)  Inclination      : ',num2str(out.i),' deg'])
        disp(['(o)  RAAN             : ',num2str(out.o),' deg'])
        disp(['(w)  Arg. Periapsis   : ',num2str(out.w),' deg'])
        disp(['(ta) True Anomaly     : ',num2str(out.ta),' deg'])
        disp(['(ma) Mean Anomaly     : ',num2str(out.ma),' deg'])
        disp(['(ma) Eccentric Anomaly: ',num2str(out.ma),' deg'])
        disp(' ')
        disp('-------------------------------------------')
        disp('Orbit Geometry:')
        disp(' ')
        disp(['( )  Conic            : ',conicType])
        disp(['( )  Direction        : ',out.orbitType])
        disp(['(rp) Periapsis Radius : ',num2str(out.rp),' km'])
        disp(['(ra) Apoapsis Radius  : ',num2str(out.ra),' km'])
        disp(['(p)  Semilatus Rectum : ',num2str(out.p),' km'])
        disp(['(N)  Nodes Vector (km):'])
        disp(out.N_)
        disp('-------------------------------------------')
        disp('Time:')
        disp(' ')
        disp(['(T)  Period           : ',num2str(out.T),' sec'])
        disp(['(dt) dt from Periapsis: ',num2str(out.dt),' sec'])
        disp('-------------------------------------------')
        disp(' ')
    end

end

