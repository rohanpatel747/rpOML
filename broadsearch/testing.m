
clear; close all; clc; format long g; rpOMLstart();

ei = getJulianDate('01-08-2037');
ji = getJulianDate('01-02-2039');
si = getJulianDate('01-09-2040');
c3 = 125;

% _________________________________________________________________________
% Inputs

in = struct;
in.verbose = true;
%in.spacing = [10; 10];
in.sequence = [
    3, ei, ei+60;
    5, ji, ji+200;
    6, si, si+400;
];
in.constrainVinf   = [sqrt(c3); 16];
in.constrainRflyby = 0;

% Constrain Vinf Setup:
% EVEEJ Example:
%   1. in.constrainVinf = [sqrt(c3); 9];    % <--- constrain E&J dep/arr
%   2. in.constrainVinf = [sqrt(c3); 5; 9; 9; 9] % <--- cstrn all legs
%   3. ________________ = % <--- do not specify any constraints








% _________________________________________________________________________
% Function
global pcd
pcd = constants();

Jmax = 0.5; 
sequence = in.sequence;
verbose  = in.verbose;

mu  = pcd.Sun.mu;

dispLnum = 10;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Setup Date Spacing

defaultSpacing = [
    10; % Mercury
    20; % Venus
    20; % Earth
    20; % Mars
    50; % Jupiter
    80; % Saturn
   120; % Uranus
   150; % Neptune
   200; % Pluto
];

if isfield(in,'spacing')
    spacing = in.spacing;
else
    for i=1:height(sequence)
        spacing(i,:) = defaultSpacing(sequence(i));
    end
end

times = {};
for i=1:height(sequence)
    times{i} = [sequence(i,2):spacing(i):sequence(i,3), sequence(i,3)];
    
    if times{i}(end) == times{i}(end-1)
        times{i} = times{i}(1:end-1);
    end
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Setup V-Infinity Constraints
cvinf = 1000*ones(height(sequence),1);
if isfield(in,'constrainVinf')
    if numel(in.constrainVinf) == 2         % if only the first/last vinf are entered
        cvinf(1)   = in.constrainVinf(1);
        cvinf(end) = in.constrainVinf(2);
    else
        cvinf = in.constrainVinf;
    end
end



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Lambert and Data Filtering
nLambertEst = 0;
for i=1:length(times)-1
    nLambertEst = nLambertEst + (length(times{i})*length(times{i+1}));
end
if verbose
    disp(' ');
    disp(['    Number of Lambert Arcs To Compute: ',num2str(nLambertEst)]);
    disp(' ')
    disp('    Starting Lambert Calculations:');
end


dataTrajUF = struct;
nLamberts = 0;
for k=1:height(sequence)-1
    
    % Departure/Arrival Body
    db = sequence(k  ,1);
    ab = sequence(k+1,1);
    
    dataTrajUF.(['l',num2str(k)]).db = db;
    dataTrajUF.(['l',num2str(k)]).ab = ab;
    
    c=1;
    
    
    for i=1:length(times{k})
        for j=1:length(times{k+1})

            % Lambert Arc
            dt = times{k}(i);
            at = times{k+1}(j);
            tof= (at-dt)*86400;
            xi = getStatePlanet(db,dt,'meeus').x;
            xf = getStatePlanet(ab,at,'meeus').x;
            l  = lambert0rev(xi,xf,tof,mu);
            
            
            % Post Processing Lambert Sol.
            if isnan(l.vi)
                if verbose
                   disp(['    [Failed Lambert]  Leg: ', num2str(k), ...
                       '      |     ', ...
                       '[',num2str(db),' - ',num2str(ab),'] : [', ...
                       num2str(dt,10),' - ',num2str(at,10),']']);
                end
 
            else
                
                % V-Infinity Computation
                vi1_ = l.vi - xi(4:6);    vi1 = norm(vi1_);
                vi2_ = l.vf - xf(4:6);    vi2 = norm(vi2_);
                
                % Valid Solution
                if (vi1<cvinf(k)) && (vi2<cvinf(k+1))
                    legdata(c,1) = c;
                    legdata(c,2) = dt;
                    legdata(c,3) = at;
                    legdata(c,4) = tof;
                    legdata(c,5) = vi1;
                    legdata(c,6) = vi2;
                    legdata(c,7:9) = vi1_.';
                    legdata(c,10:12) = vi2_.';
                    c=c+1;
                else
                   %disp('nope') 
                end
            end

            % Update #of Lamberts Calculated
            nLamberts=nLamberts+1;
            if rem(nLamberts,dispLnum) == 0
                disp(['        #L: ',num2str(nLamberts)]);
            end
        end
    end
    
    dataTrajUF.(['l',num2str(k)]).legdata = legdata;
    
    clear legdata;
    
end



dataTrajUF = broadsearch_patchLegs(dataTrajUF, 'l1', 'l2', 'p1', Jmax); 





function dataTrajUF = broadsearch_patchLegs(dataTrajUF, leg1, leg2, patchName,Jmax)

global pcd

leg1 = dataTrajUF.(leg1);
leg2 = dataTrajUF.(leg2);

dataTrajUF.([patchName]).db = leg1.db;
dataTrajUF.([patchName]).ab = leg2.ab;
dataTrajUF.([patchName]).fb = leg1.ab;
c=1;
for i=1:height(leg1.legdata)
    
    for j=1:height(leg2.legdata)
        
        % Match Dates
        date1 = leg1.legdata(i,3);
        date2 = leg2.legdata(j,2);
        
        if date1==date2
            % Get Cost
            J1 = leg1.legdata(i,6);
            J2 = leg2.legdata(j,5);
            J  = abs(J2-J1);
            
            if round(J,5)<=Jmax
                
                % Compute Flyby Parameters
                fbb  = leg1.ab;
                vi1_ = leg1.legdata(i,7:9).';
                vi2_ = leg2.legdata(j,10:12).';   
                mufb = pcd.(getPlanetName(fbb)).mu;
                bpdata = bplanefromVi1Vi2(mufb, vi1_, vi2_);
                
                % Close Approach Radius and Turning Angle
                rp = bpdata.rp;
                d  = bpdata.d;
                
                if rp > pcd.(getPlanetName(fbb)).r
                    cost(c,1)  = J;
                    patchdata(c,1) = c;
                    patchdata(c,2) = leg1.legdata(i,2);
                    patchdata(c,3) = leg2.legdata(j,3);
                    patchdata(c,4) = patchdata(c,3) - patchdata(c,2);
                    patchdata(c,5) = NaN; %vi1;
                    patchdata(c,6) = NaN; %vi2;
                    patchdata(c,7:9) = NaN; %vi1_.';
                    patchdata(c,10:12) = NaN; %vi2_.';
                    refleg(c,1) = i;
                    refleg(c,2) = j;
                    flyby(c,1) = rp;
                    flyby(c,2) = rp-pcd.(getPlanetName(fbb)).r;
                    flyby(c,3) = d;
                    c=c+1;
                end
            end
        end
    end
end
dataTrajUF.([patchName]).patchdata = patchdata;
dataTrajUF.([patchName]).flyby     = flyby;
dataTrajUF.([patchName]).refleg    = refleg;
dataTrajUF.([patchName]).cost      = cost;

end







