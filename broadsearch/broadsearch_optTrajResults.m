function out = broadsearch_optTrajResults(in)    

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

    pcd  = constants(); 

    s = in.seq;
    e = in.solution;
    ed= datetime(e,'convertfrom','juliandate');
    ed.Format = 'MMM dd, yyyy HH:mm:ss.SSS';

    mu= pcd.Sun.mu;

    noFb = height(e)-2;
    for i=1:noFb

        b1 = s(i);   e1 = e(i);
        b2 = s(i+1); e2 = e(i+1);
        b3 = s(i+2); e3 = e(i+2);

        dt1 = (e2-e1)*86400;
        dt2 = (e3-e2)*86400;

        % Query Ephemeris Data for Celestial Body States
        stateType = in.ephemType;
        sv1= getStatePlanet(b1, e1, stateType).x.';
        sv2= getStatePlanet(b2, e2, stateType).x.';
        sv3 = getStatePlanet(b3, e3, stateType).x.';

        
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
        vi1_  = l1.vi - sv1(4:6);     vi1  = norm(vi1_);
        vi2n_ = l1.vf - sv2(4:6);     vi2n = norm(vi2n_);
        vi2p_ = l2.vi - sv2(4:6);     vi2p = norm(vi2p_);
        vi3_  = l2.vf - sv3(4:6);     vi3  = norm(vi3_);

        dvi2 = vi2p-vi2n;

        if i==1
            c3   = vi1^2;
            viox = vi1_(1); vioy = vi1_(2); vioz = vi1_(3);
            vinf1_ = vi1_;
        end
        if i==noFb
            vinff_ = vi3_;
            vinff  = vi3;
            dtArr  = e3-e2;
        end

        mup = pcd.(getPlanetName(b2)).mu;
        rpl = pcd.(getPlanetName(b2)).r;
        bp  = bplanefromVi1Vi2(mup, vi2n_, vi2p_);
        rp  = bp.rp;


        % Outputs
        dtVal(i,1) = e2-e1;
        dtAdd(i,1) = e2-e(1);
        dvi(i,1) = dvi2;  
        rpca(i,1)= rp;
        ralt(i,1)= rp - pcd.(getPlanetName(b2)).r;
        br(i,1) = bp.BR;
        bt(i,1) = bp.BT;
        b_(i,:) = bp.b_.';
        vinf2n_(i,1:3) = vi2n_.';
        vinf2p_(i,1:3) = vi2p_.';
        vinfmag(i,1) = (norm(vinf2n_(i,1:3))+norm(vinf2p_(i,1:3)))/2;
        xstates(i,1:6)  = [sv1(1:3); l1.vi];
        xstates(i+1,1:6)= [sv2(1:3); l2.vi];

    end


    % Outputs
    out = struct;
    out.totaltof= e(end)-e(1);
    out.depdla  = asin(vioz/sqrt(c3))   * (180/pi);
    out.deprla  = atan2(vioy,viox) * (180/pi);
    out.depc3   = c3;
    out.xstates = xstates;
    out.br = br;
    out.bt = bt;
    out.b_ = b_;
    out.ralt = ralt;
    out.rpca = rpca;
    


    if in.plt
        options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8) ;

        figure()
        ft = tiledlayout(1,1,'Tilespacing','compact','Padding','compact');
        hold on
        for i=1:height(xstates)

            sv = getStatePlanet(s(i), e(i)).x.';
            dt = (e(i+1)-e(i))*86400;
            x_ = xstates(i,:).';

            [~,s1] = ode45(@(t,Y) eom2BP(t,Y,mu), [pcd.(getPlanetName(s(i))).t 0], sv , options);
            [~,s2] = ode45(@(t,Y) eom2BP(t,Y,mu), [0 dt], x_ , options);
            plot3(s2(:,1)./pcd.aukm,s2(:,2)./pcd.aukm,s2(:,3)./pcd.aukm)
            plot3(s1(:,1)./pcd.aukm,s1(:,2)./pcd.aukm,s1(:,3)./pcd.aukm, 'k')

        end
        %sv = getStatePlanet(s(end), e(end)).x.';
        %[~,s1] = ode45(@(t,Y) eom2BP(t,Y,mu), [pcd.(getPlanetName(s(end))).t 0], sv , options);
        %plot3(s1(:,1)./pcd.aukm,s1(:,2)./pcd.aukm,s1(:,3)./pcd.aukm, 'k')

        hold off
        grid on; box on; axis equal;
        xlabel('$X_{Inertial}$ (AU)','fontsize',12,'interpreter','Latex');
        ylabel('$Y_{Inertial}$ (AU)','fontsize',12,'interpreter','Latex');
        zlabel('$Z_{Inertial}$ (AU)','fontsize',12,'interpreter','Latex');
        set(gcf,'color','w');
        set(gca,'TickLabelInterpreter','latex');

        %print(gcf,'asen6008_labFN_EV300EJUopt.png','-dpng','-r100')



    end

    if in.disp
        dec=15;
        disp(' ');
        disp(' ');
        disp('Optimization Output');
        disp('______________________________________________________________________');
        disp(' ');
        disp(['Planetary Sequence and Encounters:']);
        for i=1:length(in.seq)
            disp(['    ',num2str(in.seq(i)),' - ',myOCD(in.seq(i)), ...
                '  -  ', char(ed(i)), '   -   ', num2str(e(i),dec), ...
                '  -  ', num2str(getGMATmodJD(e(i)),dec)]);
        end
        disp(' ');
        disp('______________________________________________________________________');
        disp(' ');
        disp('Launch Parameters:')
        disp(['    C3         (km2/s2) : ', num2str(out.depc3       ,dec)]);
        disp(['    RLA        (deg.)   : ', num2str(out.deprla      ,dec)]);
        disp(['    DLA        (deg.)   : ', num2str(out.depdla      ,dec)]);
        disp(['    |Vinf.|    (km/s)   : ', num2str(sqrt(out.depc3) ,dec)]);
        disp(['    Out. Vinf. (km/s)   : ', num2str(vinf1_.'        ,dec)]);

        disp(' ');    
        disp('______________________________________________________________________');
        disp(' ');
        disp('Flyby Parameters:');
        for i=1:noFb
            disp(['    ',num2str(in.seq(i+1)),' - ',getPlanetName(in.seq(i+1))]);
            disp(['        ','T+Days     (days)  : ',num2str(dtAdd(i)      ,dec)]);
            disp(['        ','deltaDays  (days)  : ',num2str(dtVal(i)      ,dec)]);
            disp(['        ','rp         (km)    : ',num2str(rpca(i)      ,dec)]);
            disp(['        ','ralt       (km)    : ',num2str(ralt(i)      ,dec)]);
            disp(['        ','B*R        (km)    : ',num2str(br(i)        ,dec)]);
            disp(['        ','B*T        (km)    : ',num2str(bt(i)        ,dec)]);
            disp(['        ','C3         (km2/s2): ',num2str(vinfmag(i)^2 ,dec)]);
            disp(['        ','|Vinf.|    (km/s)  : ',num2str(vinfmag(i)   ,dec)]);
            disp(['        ','Inc. Vinf. (km/s)  : ',num2str(vinf2n_(i,:) ,dec)]);
            disp(['        ','Out. Vinf. (km/s)  : ',num2str(vinf2p_(i,:) ,dec)]);
            disp('    ');
        end
        disp(' ');
        disp('______________________________________________________________________');
        disp(' ');
        disp('Arrival Parameters:')
        disp(['    T+Days     (days)  : ', num2str(out.totaltof,dec)]);
        disp(['    deltaDays  (days)  : ', num2str(dtArr       ,dec)]);
        disp(['    C3         (km2/s2): ', num2str(vinff^2     ,dec)]);
        disp(['    |Vinf.|    (km/s)  : ', num2str(vinff       ,dec)]);
        disp(['    Out. Vinf. (km/s)  : ', num2str(vinff_.'    ,dec)]);
        disp(' ');    
        disp('______________________________________________________________________');
        disp(' ');
        disp('Optimization Statistics:');
        disp('    V-Infinity Discontinuity Magnitudes (km/s)');
        for i=1:length(dvi)
            disp(['        ',num2str(dvi(i),dec)]);
        end
        disp(' ');
        disp('______________________________________________________________________');

    end

    function out = myOCD(ID)   
        p = {'Mercury';
             'Venus  ';
             'Earth  ';
             'Mars   ';
             'Jupiter';
             'Saturn ';
             'Uranus ';
             'Neptune';
             'Pluto  '};
        if ID==2020000
            out = 'Varuna';
        else
            out = p{ID};
        end
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    