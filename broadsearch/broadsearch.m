function out = broadsearch(in)
%BROADSEARCH Computes Multi-flyby Sequences Search of Solar System Planets
%
%   Assumptions/Warnings:
%   	1. This method does not support resonant orbits yet!
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%   Required Inputs: 'in' structure containing fields:
%       1. sequence     [nx3]   List of Flyby Planets and Dates Format:
%                                   in.sequence = [
%                                       Planet 1, MinDate1, MaxDate2;
%                                       Planet 2, MinDate1, MaxDate2;
%                                       Planet 3, MinDate3, MaxDate2;
%                                   ];
%                               (dates are Julian, planets are numeric 1-9)
%   Optional Inputs (fields in same structure as above) 
%       1. spacing          [nx1]   Grid spacing between dates (days)
%       2. constrainVinf    [nx1] [2x1]  Encounter V-Infinity Upper Limit
%                                   [nx1] : Constain every encounter
%                                   [2x1] : Constain Launch and Arrival
%       3. verbose          [T/F]   Print Verbose Statements to Cmd Window
%       4. Ephemeis Type    ['meeus','ephem'] Switch between Meeus Alg. or
%                                             loaded NAIF SPICE SPK for
%                                             body states.
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%   Output: out structure containing fields:
%       1. seq        [nx1]    Sequence List (Set of Numbers)
%       2. encs       [pxn]    Julian Date Array of Encounters (inc. Flybys)
%       3. encsd      [pxn]    MATLAB Calendar Date Array of Encounters (inc. FB)
%       4. cost       [pxn]    V-Infinity Mismatch at Each Flyby
%       5. costsum    [px1]    V-Infinity Mismatch Sum for Each Trajectory
%       6. fbalt      [pxn]    Flyby Altitude w.r.t FB Body (km)
%       7. fbta       [pxn]    Flyby Turning Angle (radians)
%       8. Jmax       [1x1]    Cost Maximum Per Flyby (km/s)
%
%   [RAW Output Data]
%
%       9. dataTrajUF [struct] [RAW] Output Data Structure (for debugging)
%            l#       [strcut]      Leg # Connecting Encoutners i to i+1 (format is identical to 'fp' below)
%            p#       [struct]      Patch Legs i and i+1 (format same as ^^)
%      10. fp         [strcut] [RAW] Patched Structure With Fields:
%            db       [1x1]         Departure Body Number
%            ab       [1x1]         Arrival   Body Number
%            fb       [1x1]         Last Flyby Body Number Patched
%            error    [T/F]         Error In Flyby Patching? i.e. No results found
%            legdata  [px12]        Storing Leg Data (This is explained in later documentation)
%            flyby    [qx4]         Flyby Data       (This is explained in later documentation)
%            refleg   [qx2]         Sol. From w/ Traj in Legs Joined
%            cost     [qx1]         Patch Cost 
%
%      (p = #of Valid Trajectories)     (q = #of Sequences Attempted)
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%   Dependencies:
%       1. broadsearch_patchLegs()
%       2. getPlanetName()
%       3. lambert0rev()
%       4. lambertNrev()
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%   TODO:
%       1. Inclusion of Resonant Orbits in Sequence
%       2. Inclusion of DVEGA or general orbit leveraging
%       2. Documentation On Outputs and Use Examples
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
%   Example Input:
%
%       EVEEJ Trajectory
%         e1 = juliandate('01-feb-2023','dd-mmm-yyyy');
%         e2 = juliandate('01-sep-2023','dd-mmm-yyyy');
%         e3 = juliandate('01-aug-2024','dd-mmm-yyyy');
%         e4 = juliandate('01-nov-2027','dd-mmm-yyyy');
%         e5 = juliandate('01-dec-2030','dd-mmm-yyyy') - 300;
% 
%         c3 = 30;
% 
%         in = struct;
%         in.verbose = true;                    % <--- Print to Cmd. Wndw.
%         in.spacing = [10;10;10;10;10];        % <--- Grid Search 10 days
%         in.sequence = [                       % <--- Search Input
%             3, e1, e1+150;
%             2, e2, e2+150;
%             3, e3, e3+150;
%             3, e4, e4+150;
%             5, e5, e5+600;
%         ];
%         in.constrainVinf   = [sqrt(c3); 16];  % <--- Dep. C3 and Arr Vi.
%         in.Jmax = 0.05;                       % <--- Constrain V-Inf.
%                                                      Disconnect by 50m/s
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

    %global pcd
    pcd = constants();

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Setup Date Spacing and Inputs

    sequence = in.sequence;
    mu  = pcd.Sun.mu;
    dispLnum = 50;

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

    if isfield(in,'Jmax')
        Jmax = in.Jmax;
    else
        Jmax = 1;
    end

    if isfield(in,'verbose')
        verbose = in.verbose;
    else
        verbose = false;
    end

    if isfield(in,'spacing')
        for i=1:height(sequence)
            spacing(i,:) = in.spacing;
        end
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

    if isfield(in,'ephemType')
        ephemType = in.ephemType;
    else
        ephemType = 'meeus';
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
    %if verbose
        disp(' ');
        disp(' ');
        disp(['Lambert Arcs Process']);
        disp(['    Number of Arcs To Compute: ',num2str(nLambertEst)]);
        disp('    Starting Lambert Calculations:');
    %end


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
                xi = getStatePlanet(db,dt,ephemType).x;
                xf = getStatePlanet(ab,at,ephemType).x;

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

        if exist('legdata','var')
            dataTrajUF.(['l',num2str(k)]).legdata = legdata;
        else
            disp('***No Solutions in bounds found. Broad Search will not complete.***')
            disp(' ');
            %dataTrajUF.l1.legdata

        end

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
        dataTrajUF = broadsearch_patchLegs(dataTrajUF, legi, legf, pi, Jmax, pcd);

        noPatches = noPatches+1;
        
        if dataTrajUF.([pi]).error
            disp('        No Solution. Trajectory Incomplete');
            break;
        end

    end
    noPatches  = noPatches-1;



    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Output Data
    out = struct;
    
    if height(sequence)-2 == 0 
        dataTrajUF.p1 = dataTrajUF.l1;

        for i=1:height(dataTrajUF.p1.legdata)
            out.encs(i,1)  = dataTrajUF.p1.legdata(i,2);
            out.encs(i,2)  = dataTrajUF.p1.legdata(i,3);
            out.encsd(i,1) = datetime(out.encs(i,1),'convertfrom','juliandate');
            out.encsd(i,2) = datetime(out.encs(i,2),'convertfrom','juliandate');
        end
        
        out.cost       = NaN;
        out.costsum    = NaN;
        out.fbalt      = NaN;
        out.fbta       = NaN;
        out.fp         = NaN;
        out.Jmax       = NaN;
        out.error      = false;
        
    elseif dataTrajUF.(['p',num2str(noPatches)]).error == 0
        
        fp = dataTrajUF.(['p',num2str(noPatches)]);

        j=1; seqL = height(sequence);
        for i=1:height(fp.cost)
            if ~isnan(fp.cost(i,:))

                cost(j,:)    = fp.cost(i,:);
                costsum(j,1) = sum(cost(j,:));
                encs(j,1)    = fp.legdata(i,2);
                encsd(j,1)   = datetime(encs(1),'convertfrom','juliandate');
                fbalt(j,1 )  = NaN;
                fbta(j,1)    = NaN;
                kk=2;
                for k=1:4:width(fp.flyby)
                    encs(j,kk)  = fp.flyby(i,k);
                    encsd(j,kk) = datetime(fp.flyby(i,k),'convertfrom','juliandate');
                    fbta(j,kk)  = fp.flyby(i,k+3);
                    fbalt(j,kk) = fp.flyby(i,k+2);
                    kk=kk+1;
                end
                encs(j,seqL)  = fp.legdata(i,3);
                fbalt(j,seqL) = NaN;
                fbta(j,seqL)  = NaN;
                encsd(j,seqL) = datetime(encs(end),'convertfrom','juliandate');

                j=j+1;
            end

        end
        
        out.encs       = encs;
        out.encsd      = encsd;
        out.cost       = cost;
        out.costsum    = costsum;
        out.fbalt      = fbalt;
        out.fbta       = fbta;
        out.fp         = fp;
        out.Jmax       = Jmax;
        out.error      = false;
    else
        disp('Unable to patch trajectory and output results. Error flag thrown.')
        disp('Returning Individual Leg Data.')
        out.error      = true;
        out.forceplot  = true;
    end

    
    out.pcd        = pcd;
    out.mu         = mu;
    out.seq        = sequence(:,1);
    out.dataTrajUF = dataTrajUF;
    out.ephemType  = ephemType;


end

