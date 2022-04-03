function out = pkchp_porkchop(in)

    mu = in.mu;
    depBody = in.depBody;
    arrBody = in.arrBody;
    depDate = in.depDate;
    arrDate = in.arrDate;
    dayStep = in.dayStep;
    daysPastDep = in.daysPastDep;
    daysPastArr = in.daysPastArr;
    ntof = in.ntof;

    dp1 = 0:dayStep:daysPastDep;
    dp2 = 0:dayStep:daysPastArr;

    % Accept Either Julian Date or Calendar Date
    if ischar(depDate)
        jddep = juliandate(depDate,'dd-mmm-yyyy');
    elseif isnumeric(depDate)
        jddep = depDate;
    end
    if ischar(arrDate)
        jdarr = juliandate(arrDate,'dd-mmm-yyyy');
    elseif isnumeric(arrDate)
        jdarr = arrDate;
    end  
    
    jdd = jddep:dayStep:(jddep+daysPastDep);
    jda = jdarr:dayStep:(jdarr+daysPastArr);


    for i=1:length(jdd)
        for j=1:length(jda)
            
            % Lambert's Problem Setup
            xd = getStatePlanet(depBody, jdd(i), 'meeus').x;
            xa = getStatePlanet(arrBody, jda(j), 'meeus').x;
            
            rd = xd(1:3);
            vd = xd(4:6);
            ra = xa(1:3);
            va = xa(4:6);

            tof(j,i) = (jda(j) - jdd(i));
            dt = tof(j,i)*86400;

            % Lambert Arc
            l = lambert0rev(rd, ra, dt, mu);

            % Outgoing/Incomig Vinfinity 
            vio_ = l.vi - vd;
            vif_ = l.vf - va;
            
            viox = vio_(1);
            vioy = vio_(2);
            vioz = vio_(3);
            
            % Storing Data
            vinfdep(j,i) = norm(vio_);
            vinfarr(j,i) = norm(vif_);
            c3dep(j,i)   = vinfdep(j,i)^2;
            rla(j,i)     = atan2(vioy,viox) * (180/pi);
            dla(j,i)     = asin(vioz/vinfdep(j,i)) * (180/pi);

        end

    end

    %out         = struct;
    out         = in;
    out.c3dep   = c3dep;
    out.vinfdep = vinfdep;
    out.vinfarr = vinfarr;
    out.tof     = tof;
    out.ntof    = ntof;
    out.dp1     = dp1;
    out.dp2     = dp2;
    out.depDate = depDate;
    out.arrDate = arrDate;
    out.depBody = depBody;
    out.arrBody = arrBody;
    out.RLAdep  = rla;
    out.DLAdep  = dla;
    %out.plts    = in.plts;
end
