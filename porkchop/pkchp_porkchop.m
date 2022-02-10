function out = pkchp_porkchop(in)

    mu = in.mu;
    depBody = in.depBody;
    arrBody = in.arrBody;
    depDate = in.depDate;
    arrDate = in.arrDate;
    dayStep = in.dayStep;
    daysPastDep = in.daysPastDep;
    daysPastArr = in.daysPastArr;
    ctr1 = in.ctr1;
    ctr2 = in.ctr2;
    ntof = in.ntof;

    dp1 = 0:dayStep:daysPastDep;
    dp2 = 0:dayStep:daysPastArr;

    jddep = juliandate(depDate,'dd-mmm-yyyy');
    jdarr = juliandate(arrDate,'dd-mmm-yyyy');
    jdd = jddep:dayStep:(jddep+daysPastDep);
    jda = jdarr:dayStep:(jdarr+daysPastArr);


    for i=1:length(jdd)
        for j=1:length(jda)

            xd = getStatePlanet(depBody, jdd(i), 'meeus').x;
            xa = getStatePlanet(arrBody, jda(j), 'meeus').x;

            rd = xd(1:3);
            vd = xd(4:6);

            ra = xa(1:3);
            va = xa(4:6);


            tof(j,i) = (jda(j) - jdd(i));
            dt = tof(j,i)*86400;
            
            
            %disp(jdd(i))
            %disp(jda(j))
            %disp(tof(j,i))
            %disp(' ')
            %disp(' ')
            %disp(' ')
            
            
            out = lambert0rev(rd, ra, dt, mu);

            vinfdep(j,i) = norm(out.vi - vd);
            vinfarr(j,i) = norm(out.vf - va);
            c3dep(j,i)   = vinfdep(j,i)^2;        

        end

    end


    out         = struct;
    out.c3dep   = c3dep;
    out.vinfarr = vinfarr;
    out.tof     = tof;
    out.ctr1    = ctr1;
    out.ctr2    = ctr2;
    out.ntof    = ntof;
    out.dp1     = dp1;
    out.dp2     = dp2;
    out.depDate = depDate;
    out.arrDate = arrDate;
    out.depBody = depBody;
    out.arrBody = arrBody;
end
