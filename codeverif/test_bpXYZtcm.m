%% rpOML : Testing B-Plane Functions with TCA
%  C: 09MAR22
clear; close all; clc; format long g; rpOMLstart();

%% Problem Setup
% Constants
pcd = constants();
mu = pcd.Earth.mu;

% Initial State Vector
x_ = [    546507.344255845
         -527978.380486028
          531109.066836708
          -4.9220589268733
          5.36316523097915
         -5.22166308425181 ];

% Desired B-Plane State     
BRt = 5022.26511510685;
BTt = 13135.7982982557;
TCAt= 102405.381807303;         % Holding TCA Constant as TCA from x_
BTarg = [BRt; BTt; TCAt];


%% B-Plane and Associated Functions
% Hyperbolic Anomaly Parameters
s2e = conv_state2ele(x_,mu,false);      % State to Elements
an  = getAnomalyandDt(s2e, false);      % Mean, Eccentric Anomalies, and Time to Periapsis


% Full B-Plane State (BR, BT, TCA)
[BRn, BTn] = bplaneBRBTfromRV(mu,x_);
TCAn       = bplaneTCA(mu,x_);


% Full 3x3 Jacobian
jBjV = bplane_getJacobian(mu, x_, 0.001, 'npS', 3);


% TCM Target VX VY and VZ 
out = bplane_computeXYZTCM(mu,x_,BTarg);


% Check Values:
% Finalized B-Plane Values (Post-Maneuver)
TCMDVMag   = out.DV
[BRf, BTf] = bplaneBRBTfromRV(mu,out.xf_)
TCAf       = bplaneTCA(mu,out.xf_)



%% GMAT State Data
%
%   DefaultSC     : Initial State w/o TCM Applied
%   DefaultSC_wTCM: Initial State w/  TCM Applied
%
%   Mission Sequence Summary:
%       1. [Both s/cs] Propagate 1 second
%       2. [Both s/cs] Propagate to Computed TCA (TCAn = TCAt)
%
%
%
%
%   Command Summary after 1 second propagation:
%{
******  Changes made to the mission will not be reflected ******
******  in the data displayed until the mission is rerun  ******

        Propagate Command: Propagate1
        Spacecraft       : DefaultSC
        Coordinate System: EarthMJ2000Ec

        Time System   Gregorian                     Modified Julian  
        ----------------------------------------------------------------------    
        UTC Epoch:    01 Jan 2000 11:59:29.000      21544.9996412037
        TAI Epoch:    01 Jan 2000 12:00:01.000      21545.0000115741
        TT  Epoch:    01 Jan 2000 12:00:33.184      21545.0003840741
        TDB Epoch:    01 Jan 2000 12:00:33.184      21545.0003840732

        Cartesian State                       Keplerian State
        ---------------------------           -------------------------------- 
        X  =   546502.42219678 km             SMA  =  -5020.4044856722 km
        Y  =  -527973.01732067 km             ECC  =   9.4352521824899
        Z  =   531103.84517349 km             INC  =   37.657087231646 deg
        VX =  -4.9220592002510 km/sec         RAAN =   200.90708855255 deg
        VY =   5.3631654950880 km/sec         AOP  =   203.51066523767 deg
        VZ =  -5.2216633499268 km/sec         TA   =   266.82044230827 deg
                                              MA   =  -10413.626494605 deg
                                              HA   =  -210.38932058248 deg

        Spherical State                       Other Orbit Data
        ---------------------------           --------------------------------
        RMAG =   927087.75142578 km           Mean Motion        =   1.774846606e-03 deg/sec
        RA   =  -44.012028164513 deg          Orbit Energy       =   39.698040530158 km^2/s^2
        DEC  =   34.950835662072 deg          C3                 =   79.396081060316 km^2/s^2
        VMAG =   8.9585701447522 km/s         Semilatus Rectum   =   441916.00285108 km   
        AZI  =   105.00915365790 deg          Angular Momentum   =   419699.79013856 km^2/s
        VFPA =   177.10341115343 deg          Beta Angle         =  -59.778187109961 deg  
        RAV  =   132.54423929382 deg          Periapsis Altitude =   35970.241594748 km   
        DECV =  -35.652512307961 deg          VelPeriapsis       =   9.9106461924392 km/s

        Planetodetic Properties               Hyperbolic Parameters
        ---------------------------           --------------------------------
        LST       =   308.15315466661 deg     BdotT          =   45892.323795435 km   
        MHA       =   280.33258104394 deg     BdotR          =   10606.210428744 km   
        Latitude  =   17.402759898789 deg     B Vector Angle =   13.013194891089 deg  
        Longitude =   27.820573622675 deg     B Vector Mag   =   47101.985977280 km   
        Altitude  =   920711.52498354 km      DLA            =  -32.088775292655 deg  
                                              RLA            =   146.56286750075 deg  

   ========================================================================

        Spacecraft       : DefaultSC_wTCM
        Coordinate System: EarthMJ2000Ec

        Time System   Gregorian                     Modified Julian  
        ----------------------------------------------------------------------    
        UTC Epoch:    01 Jan 2000 11:59:29.000      21544.9996412037
        TAI Epoch:    01 Jan 2000 12:00:01.000      21545.0000115741
        TT  Epoch:    01 Jan 2000 12:00:33.184      21545.0003840741
        TDB Epoch:    01 Jan 2000 12:00:33.184      21545.0003840732

        Cartesian State                       Keplerian State
        ---------------------------           -------------------------------- 
        X  =   546502.19928783 km             SMA  =  -5074.9503606625 km
        Y  =  -527973.23602741 km             ECC  =   2.9460061409372
        Z  =   531103.92227933 km             INC  =   40.301208276968 deg
        VX =  -5.1449681514028 km/sec         RAAN =   191.48674522828 deg
        VY =   5.1444587506178 km/sec         AOP  =   226.63882398979 deg
        VZ =  -5.1445575108288 km/sec         TA   =   251.02387841654 deg
                                              MA   =  -10246.188990192 deg
                                              HA   =  -276.49929796768 deg

        Spherical State                       Other Orbit Data
        ---------------------------           --------------------------------
        RMAG =   927087.78874936 km           Mean Motion        =   1.746309392e-03 deg/sec
        RA   =  -44.012051702549 deg          Orbit Energy       =   39.271363577236 km^2/s^2
        DEC  =   34.950839863724 deg          C3                 =   78.542727154471 km^2/s^2
        VMAG =   8.9108150636169 km/s         Semilatus Rectum   =   38970.301147782 km   
        AZI  =   111.49167703177 deg          Angular Momentum   =   124633.78050470 km^2/s
        VFPA =   179.13555655304 deg          Beta Angle         =  -63.334032850288 deg  
        RAV  =   135.00283655400 deg          Periapsis Altitude =   3497.7482668006 km   
        DECV =  -35.263570963194 deg          VelPeriapsis       =   12.620011874549 km/s

        Planetodetic Properties               Hyperbolic Parameters
        ---------------------------           --------------------------------
        LST       =   308.15313405827 deg     BdotT          =   13135.798299791 km   
        MHA       =   280.33258104394 deg     BdotR          =   5022.2651156885 km   
        Latitude  =   17.402758123844 deg     B Vector Angle =   20.923565597896 deg  
        Longitude =   27.820553014326 deg     B Vector Mag   =   14063.155544367 km   
        Altitude  =   920711.56230675 km      DLA            =  -14.957485399549 deg  
                                              RLA            =   173.12570305467 deg  

%}
%
%
%
%
%   Command Summary at TCA propagation:
%{
******  Changes made to the mission will not be reflected ******
******  in the data displayed until the mission is rerun  ******

        Propagate Command: Propagate2
        Spacecraft       : DefaultSC
        Coordinate System: EarthMJ2000Ec

        Time System   Gregorian                     Modified Julian  
        ----------------------------------------------------------------------    
        UTC Epoch:    02 Jan 2000 16:26:14.382      21546.1848886783
        TAI Epoch:    02 Jan 2000 16:26:46.382      21546.1852590487
        TT  Epoch:    02 Jan 2000 16:27:18.566      21546.1856315487
        TDB Epoch:    02 Jan 2000 16:27:18.566      21546.1856315483

        Cartesian State                       Keplerian State
        ---------------------------           -------------------------------- 
        X  =   31497.116035577 km             SMA  =  -5020.3901069198 km
        Y  =   26356.724791587 km             ECC  =   9.4352675344847
        Z  =  -10326.480186618 km             INC  =   37.657149307226 deg
        VX =  -6.2609083211972 km/sec         RAAN =   200.90704797390 deg
        VY =   5.3100493936274 km/sec         AOP  =   203.51083393524 deg
        VZ =  -5.5520768160946 km/sec         TA   =   0.0133418981541 deg
                                              MA   =   0.1011845317601 deg
                                              HA   =   0.0119954145520 deg

        Spherical State                       Other Orbit Data
        ---------------------------           --------------------------------
        RMAG =   42348.334717468 km           Mean Motion        =   1.774854230e-03 deg/sec
        RA   =   39.922526711505 deg          Orbit Energy       =   39.698154228153 km^2/s^2
        DEC  =  -14.113657560191 deg          C3                 =   79.396308456307 km^2/s^2
        VMAG =   9.9106586330584 km/s         Semilatus Rectum   =   441916.19158364 km   
        AZI  =   125.28168800217 deg          Angular Momentum   =   419699.87976081 km^2/s
        VFPA =   89.987936640962 deg          Beta Angle         =  -59.908808478220 deg  
        RAV  =   139.69780162231 deg          Periapsis Altitude =   35970.197379349 km   
        DECV =  -34.070508532061 deg          VelPeriapsis       =   9.9106586563398 km/s

        Planetodetic Properties               Hyperbolic Parameters
        ---------------------------           --------------------------------
        LST       =   41.928967866394 deg     BdotT          =   45892.239551170 km   
        MHA       =   348.18990399521 deg     BdotR          =   10606.320063631 km   
        Latitude  =   1.3655240858234 deg     B Vector Angle =   13.013347903493 deg  
        Longitude =   53.739063871183 deg     B Vector Mag   =   47101.928583808 km   
        Altitude  =   35970.210528824 km      DLA            =  -32.088769192549 deg  
                                              RLA            =   146.56302471478 deg  

   ========================================================================

        Spacecraft       : DefaultSC_wTCM
        Coordinate System: EarthMJ2000Ec

        Time System   Gregorian                     Modified Julian  
        ----------------------------------------------------------------------    
        UTC Epoch:    02 Jan 2000 15:44:06.868      21546.1556350437
        TAI Epoch:    02 Jan 2000 15:44:38.868      21546.1560054141
        TT  Epoch:    02 Jan 2000 15:45:11.052      21546.1563779141
        TDB Epoch:    02 Jan 2000 15:45:11.052      21546.1563779136

        Cartesian State                       Keplerian State
        ---------------------------           -------------------------------- 
        X  =   24868.764816363 km             SMA  =  -5074.9077949105 km
        Y  =  -7790.2111812607 km             ECC  =   2.9460199517975
        Z  =   10674.606126115 km             INC  =   40.301494348296 deg
        VX =  -6.2522622979530 km/sec         RAAN =   191.48629833200 deg
        VY =   5.7806925450098 km/sec         AOP  =   226.63912577313 deg
        VZ =  -5.8603772600018 km/sec         TA   =   277.48563472265 deg
                                              MA   =  -252.79520110383 deg
                                              HA   =  -82.340751936949 deg

        Spherical State                       Other Orbit Data
        ---------------------------           --------------------------------
        RMAG =   28161.854869449 km           Mean Motion        =   1.746331363e-03 deg/sec
        RA   =  -17.393313118202 deg          Orbit Energy       =   39.271692965510 km^2/s^2
        DEC  =   22.274531171518 deg          C3                 =   78.543385931020 km^2/s^2
        VMAG =   10.336885980404 km/s         Semilatus Rectum   =   38970.387252699 km   
        AZI  =   124.49748918030 deg          Angular Momentum   =   124633.91819385 km^2/s
        VFPA =   154.65045597811 deg          Beta Angle         =  -63.224499371536 deg  
        RAV  =   137.24426480174 deg          Periapsis Altitude =   3497.7355224287 km   
        DECV =  -34.537003545529 deg          VelPeriapsis       =   12.620042102086 km/s

        Planetodetic Properties               Hyperbolic Parameters
        ---------------------------           --------------------------------
        LST       =   335.38589368584 deg     BdotT          =   13135.717683337 km   
        MHA       =   337.62976193068 deg     BdotR          =   5022.3543267135 km   
        Latitude  =   13.772120260868 deg     B Vector Angle =   20.924022388525 deg  
        Longitude =  -2.2438682448431 deg     B Vector Mag   =   14063.112103635 km   
        Altitude  =   21784.926815485 km      DLA            =  -14.957449757046 deg  
                                              RLA            =   173.12549607547 deg  

%}
