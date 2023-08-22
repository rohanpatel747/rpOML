function in = broadsearch_optTraj(in, opts)
%BROADSEARCH_OPTTRAJ Minimize V-Inf. Discontinuities Between Encounters
%
%   *** Heliocentric and Planetary ONLY! ***
%
%   Inputs: 
%   1. "in" structure containing fields:
%       in.sequence = [
%           body1ID, JD lower, JD guess, JD upper;
%           body2ID, "   "   , "   "   , "   "   ;
%               ...
%           bodyNID, "   "   , "   "   , "   "   ];
%
%       in.optType = 'minTOF' or 'minC3';
%       in.c3max   = max launch energy (km3/s2)
%       in.rminFB  = minimum flyby altitude (for all bodies)
%       in.ephemType = 'meeus' or 'ephem'
%
%   2. useParallel = Parallel Processing? Auto-detects #Cores
%   3. cspiceSPK   = {{body1.bsp}, {body2.bsp}} etc... 
%
%

    arguments
        in;
        opts.useParallel = false;
        opts.cspiceSPK   = {};
    end


    if true
        disp(' '); disp(' ');
        disp('! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !');
        disp(' ');
        disp('***WARNING!***')
        disp('    If you run into errors in Lambert solutions not being computed,')
        disp('    it is likely the trajectory is very sensitive and the current')
        disp('    Lambert algorithm cannot find a solution.')
        disp(' ')
        disp('    You can try adding this Lambert algorithm:')
        disp('        https://www.mathworks.com/matlabcentral/fileexchange/26348-robust-solver-for-lambert-s-orbital-boundary-value-problem')
        disp('    ');
        disp('    to your pwd. The opt function is equiped to handle it ^^^');
        disp('    The Lancaster-Blanchard method is more robust but slower');
        disp(' ');
        disp('! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !');
        disp(' ');
        disp(' ');
    end


    
    pcd = constants();
    seq = in.sequence(:,1);
    SPK = opts.cspiceSPK;

    
    % Opt Problem Adtl. Variables
    noEnc= height(in.sequence);
    noFb = noEnc-2;
    encs = optimvar("encs",noEnc,1, ...
        "LowerBound",in.sequence(:,2), ...
        "UpperBound",in.sequence(:,4));
    
    % Cost and Constraints
    c1        = @(encs) jC3(encs, seq, pcd.Sun.mu);
    c2        = @(encs) jTOF(encs);
    nlcstr1   = @(encs) cstr2bp(encs, seq, pcd, noFb, in.ephemType);
    [vM, rM, c3, de]  = fcn2optimexpr(nlcstr1, encs, ...
        'OutputSize',{[noFb,1],[noFb,1],[1,1], [noEnc-1, 1]});
    
    % Problem Setup
    prob                = optimproblem;
    prob.ObjectiveSense = 'minimize';
    if contains(in.optType, 'minTOF')
        prob.Objective.o1   = fcn2optimexpr(c2,encs); % Min. Trip Time
    elseif contains(in.optType, 'minC3')
        prob.Objective.o1   = c3;
    else
        disp(' ');
        disp('Cost Function Not Specified. Using Minimum Launch C3.');
        disp(' ');
        prob.Objective.o1   = c3;
    end
    
    prob.Constraints.c1 = vM == 0;
    prob.Constraints.c2 = rM >= in.rminFb;
    prob.Constraints.c3 = c3 <= in.c3max;
    prob.Constraints.c4 = de >= 50;
    
    optopts = optimoptions("fmincon");
    optopts.ConstraintTolerance    = 1e-06;
    optopts.MaxIterations          = 1e+04;
    optopts.MaxFunctionEvaluations = 1e+05;
    optopts.StepTolerance          = 1e-12;
    optopts.PlotFcn                = 'optimplotfval';
    optopts.Display                = 'iter';
    optopts.Algorithm              = 'sqp';
    

    % Parallel Processing
    optopts.UseParallel            = opts.useParallel;
    if ~isempty(opts.cspiceSPK)
        core_info = evalc('feature(''numcores'')');
        corecount = str2num(core_info(end-3:end-2));
        parfor i=1:corecount
            rpOMLstart('cspice',true);
            for j=1:length(SPK)
                cspice_furnsh(SPK(j));
            end
        end
    else
        rpOMLstart('cspice',true);
        for j=1:length(SPK)
            cspice_furnsh(SPK(j));
        end

    end

    
    % Solve Problem
    tic
    x0.encs = in.sequence(:,3);
    [sol,fval,exitflag,output,lambda] = solve(prob,x0, ...
        'Options',optopts, ...
        'ObjectiveDerivative' ,"auto", ...
        'ConstraintDerivative',"auto");
    eSol = sol.encs;
    in.solution = eSol;
    toc
    
    %[vM, rM, c3, de] = trajOpt_c1(in.sequence(:,2), seq, pcd, noFb);
    
    in.seq  = in.sequence(:,1);
    in.disp = true;
    in.plt  = true;
    in.fval = fval;
    in.exitflag = exitflag;
    in.output = output;
    in.lambda = lambda;


    % ---------------------------------------------------------------------
    % Optimization Functions


    % [COST] MinTOF
    function dtTotal = jTOF(encs)
        dtTotal = (encs(end)-encs(1));
    end
    
    % [COST] MinC3
    function c3 = jC3(encs, seq, mu) 
        sv1 = getStatePlanet(seq(1), encs(1)).x.';
        sv2 = getStatePlanet(seq(2), encs(2)).x.';
        dt = (encs(2)-encs(1))*86400;
        l1 = lambert0rev(sv1, sv2, dt, mu);
        vi = l1.vi - sv1(4:6);
        c3 = norm(vi)^2; 
    end
    
    % [CONSTRAINT] Heliocentric 2-Body Dynamics
    function [vM, rM, c3, de] = cstr2bp(encs, seq, pcd, noFb, ephemType)
        mu = pcd.Sun.mu;
        
        vM = zeros(noFb,1);
        rM = zeros(noFb,1);
        
        for i=1:noFb
            
            b1 = seq(i);
            b2 = seq(i+1);
            b3 = seq(i+2);
            
            e1 = encs(i);
            e2 = encs(i+1);
            e3 = encs(i+2);
            
            dt1= (e2-e1)*86400;
            dt2= (e3-e2)*86400;
            
            % Query Ephemeris Data (DE421.bsp) for Planet & KBO States
            stateType = ephemType;
            sv1= getStatePlanet(b1, e1, stateType).x.';
            sv2= getStatePlanet(b2, e2, stateType).x.';
            if b3 == 2020000
                stateType = ephemType;
                sv3 = getStatePlanet(b3, e3, stateType).x.';
            else
                sv3 = getStatePlanet(b3, e3, stateType).x.';
            end
        
            
        
            % Computing Lambert Arcs
            if (b1==b2) && (dt1>pcd.(getPlanetName(b2)).t)
                l1 = lambertNrev(sv1, sv2, dt1, mu, 3);
        
                if isnan(l1.vi(1))
                    l1 = lambertLancasterBlanchard(sv1, sv2, dt1, mu, 3);
                end
        
            else
                l1 = lambert0rev(sv1, sv2, dt1, mu);
        
                if isnan(l1.vi(1))
                    l1 = lambertLancasterBlanchard(sv1, sv2, dt1, mu, 0);
                end
        
            end
        
            if (b2==b3) && (dt2>pcd.(getPlanetName(b3)).t)
                l2 = lambertNrev(sv2, sv3, dt2, mu, 3);
        
                if isnan(l2.vi(1))
                    l2 = lambertLancasterBlanchard(sv2, sv3, dt2, mu, 3);
                end
        
            else
                l2 = lambert0rev(sv2, sv3, dt2, mu);
        
                if isnan(l2.vi(1))
                    l2 = lambertLancasterBlanchard(sv2, sv3, dt2, mu, 0);
                end
        
            end
        
            % V-Infinity    
            
            if (isnan(l1.vi(1))) || (isnan(l2.vi(1)))
                
                vM(i,1) = 1e10;
                rM(i,1) = -1e10;
                
                if i==1; c3 = 1e10; end
        
            else
                
                vinf1 = l1.vi - sv1(4:6);
                vinf2n= l1.vf - sv2(4:6);
                vinf2p= l2.vi - sv2(4:6);
                vinf3 = l2.vf - sv3(4:6);
        
                if i==1; c3 = norm(vinf1)^2; end
        
                muplanet= pcd.(getPlanetName(b2)).mu;
                rplanet = pcd.(getPlanetName(b2)).r;
        
                bp2 = bplanefromVi1Vi2(muplanet, vinf2n, vinf2p);
                rp  = bp2.rp;
                
                vM(i,1) = norm(vinf2p) - norm(vinf2n);
                rM(i,1) = rp - rplanet;
        
            end
            
        end
        
        for i=1:length(encs)-1
           de(i,1) = encs(i+1)-encs(i);
        end
        

            function out = lambertLancasterBlanchard(sv1, sv2, dt1, mu, rev)
                
                % lambert revs
                if rev==3; rNum = 1; else; rNum = 0; end
        
                r1vec = sv1(1:3).';
                r2vec = sv2(1:3).';
                muC   = mu;
                tf    = dt1;
        
                [V1, V2, ~, ~] = lambert_LancasterBlanchard( ...
                    r1vec, r2vec, tf, rNum, muC);
        
                out = struct;
                out.vi = V1.';
                out.vf = V2.';
        
            end   
        
        end

end
