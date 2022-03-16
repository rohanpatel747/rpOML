
clear; close all; clc; format long g; rpOMLstart();



% _________________________________________________________________________
% Inputs
if false
ei = getJulianDate('01-08-2037');
ji = getJulianDate('01-02-2039');
si = getJulianDate('01-09-2040');
c3 = 125;
    
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


end





e1 = getJulianDate('01-03-2023');
e2 = getJulianDate('01-09-2023');
e3 = getJulianDate('01-08-2024');
e4 = getJulianDate('01-11-2027');
e5 = getJulianDate('01-12-2030');

c3 = 25;
    
in = struct;
in.verbose = true;
%in.spacing = [10; 10];
in.sequence = [
    3, e1, e1+100;
    2, e2, e2+100;
    3, e3, e3+100;
    3, e4, e4+150;
    5, e5, e5+300;
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
% Lambert Arcs
nLambertEst = 0;
for i=1:length(times)-1
    nLambertEst = nLambertEst + (length(times{i})*length(times{i+1}));
end
if verbose
    disp(' ');
    disp(' ');
    disp(['Lambert Arcs Process']);
    disp(['    Number of Arcs To Compute: ',num2str(nLambertEst)]);
    disp('    Starting Lambert Calculations:');
end


dataTrajUF = struct;
nLamberts = 0;
ndiscarded= 0;
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
            
            if (db==ab) && (tof>pcd.(getPlanetName(db)).t)
                l = lambertNrev(xi,xf,tof,mu,3);
            else
                l = lambert0rev(xi,xf,tof,mu);
            end
            

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
                   ndiscarded=ndiscarded+1;
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
if verbose
disp(['    Number of Solutions Discarded: ',num2str(ndiscarded)]);
disp(' ');
end








% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Patching Arcs
if verbose
    disp(['Patching Lambert Arcs']);
    disp(['    Patching Legs:']);
end
noPatches = 1;
for i=1:height(sequence)-2
    
    if verbose
       disp(['        ', num2str(i), '-', num2str(i+1)]); 
    end
    
    
    if i==1
        legi = ['l',num2str(i)];
        legf = ['l',num2str(i+1)];
    else
        legi = ['p',num2str(i-1)];
        legf = ['l',num2str(i+1)];
    end
    pi = ['p',num2str(noPatches)];

    %disp(legi); disp(legf); disp(pi); disp(' ');
    dataTrajUF = broadsearch_patchLegs(dataTrajUF, legi, legf, pi, Jmax);
    
    if dataTrajUF.([pi]).error
        disp('        No Solution. Trajectory Incomplete');
        break;
    end
    
    noPatches = noPatches+1;
end
noPatches  = noPatches-1;



% Output Data
fp = dataTrajUF.(['p',num2str(noPatches)]);

j=1;
for i=1:height(fp.cost)
    if ~isnan(fp.cost(i,:))

        cost(j,:) = fp.cost(i,:);
        costsum(j,1) = sum(cost(j,:));
        
        encs = [fp.legdata(i,2)];
        encsd= [datetime(encs,'convertfrom','juliandate')];
        fbalt= [NaN];
        fbta = [NaN];
        for k=1:4:width(fp.flyby)
            encs = [encs , fp.flyby(i,k)];
            encsd= [encsd; datetime(fp.flyby(i,k),'convertfrom','juliandate')];
            fbalt= [fbalt, fp.flyby(i,k+2)];
            fbta = [fbta , fp.flyby(i,k+3)];
        end
        encs  = [encs, fp.legdata(i,3)];
        fbalt = [fbalt, NaN];
        fbta  = [fbta, NaN];
        encsd= [encsd; datetime(encs(end),'convertfrom','juliandate')];

        j=j+1;
    end

end

out            = struct;
out.seq        = sequence(:,1);
out.encs       = encs;
out.encsd      = encsd;
out.cost       = cost;
out.costsum    = costsum;
out.fbalt      = fbalt;
out.fbta       = fbta;
out.dataTrajUF = dataTrajUF;
out.fp         = fp;









function dataTrajUF = broadsearch_patchLegs(dataTrajUF, leg1, leg2, patchName,Jmax)

global pcd

leg1 = dataTrajUF.(leg1);
leg2 = dataTrajUF.(leg2);
pnum = str2double(patchName(2:end));

dataTrajUF.([patchName]).db = leg1.db;
dataTrajUF.([patchName]).ab = leg2.ab;
dataTrajUF.([patchName]).fb = leg1.ab;

if pnum==1
    cost(1:height(leg1.legdata),1)     = NaN;
    refleg(1:height(leg1.legdata),1:2) = NaN;
    flyby(1:height(leg1.legdata),1:4)  = NaN;
elseif pnum>1
    
    cost   = dataTrajUF.(['p',num2str(pnum-1)]).cost;
    cost(:,((pnum-1)*1)+1) = NaN;
    
    refleg = dataTrajUF.(['p',num2str(pnum-1)]).refleg;
    refleg(:,((pnum-1)*2)+1:((pnum-1)*2)+2) = NaN;
    
    flyby  = dataTrajUF.(['p',num2str(pnum-1)]).flyby;
    flyby(:,((pnum-1)*4)+1:((pnum-1)*4)+4) = NaN;
end


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
                    %c=i;
                    legdata(i,1)     = i;
                    legdata(i,2)     = leg1.legdata(i,2);
                    legdata(i,3)     = leg2.legdata(j,3);
                    legdata(i,4)     = (legdata(i,3) - legdata(i,2))*86400;
                    legdata(i,5)     = leg1.legdata(i,5);
                    legdata(i,6)     = leg2.legdata(j,6);
                    legdata(i,7:9)   = leg1.legdata(i,7:9);
                    legdata(i,10:12) = leg2.legdata(j,10:12);
                    %if true
                        cost(i,((pnum-1)*1)+1)  = J;
                        refleg(i,((pnum-1)*2)+1) = i;
                        refleg(i,((pnum-1)*2)+1+1) = j;
                        flyby(i,((pnum-1)*4)+1) = date1;
                        flyby(i,((pnum-1)*4)+2) = rp;
                        flyby(i,((pnum-1)*4)+3) = rp-pcd.(getPlanetName(fbb)).r;
                        flyby(i,((pnum-1)*4)+4) = d;
                    %{
                    else
                        cost(c,1)  = J;
                        refleg(c,1) = i;
                        refleg(c,2) = j;
                        flyby(c,1)  = rp;
                        flyby(c,2)  = rp-pcd.(getPlanetName(fbb)).r;
                        flyby(c,3)  = d;
                    end
                    c=c+1;
                    %}
                else
                    legdata(i,1:12) = NaN;
                    cost(i,((pnum-1)*1)+1)  = NaN;
                    refleg(i,((pnum-1)*2)+1) = NaN;
                    refleg(i,((pnum-1)*2)+1+1) = NaN;
                    flyby(i,((pnum-1)*3)+1) = NaN;
                    flyby(i,((pnum-1)*3)+1+1) = NaN;
                    flyby(i,((pnum-1)*3)+1+2) = NaN;
                end
            end
        end
    end
end

if isnan(cost(:,((pnum-1)*1)+1))
    dataTrajUF.([patchName]).error = true;
else
    dataTrajUF.([patchName]).error   = false;
    dataTrajUF.([patchName]).legdata = legdata;
    dataTrajUF.([patchName]).flyby   = flyby;
    dataTrajUF.([patchName]).refleg  = refleg;
    dataTrajUF.([patchName]).cost    = cost;
end

end







