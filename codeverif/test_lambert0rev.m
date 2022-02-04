%% rpOML + ASEN6008 Planetary States Check File 
%  C: 03FEB22

clear; clc; format long g; rpOMLstart();

pcd = constants(); 
mu  = pcd.Sun.mu;

use6008testData = true;


%% Function:
%       ./lambert/lambert0rev()
        help lambert0rev


%% Test Data from ASEN6008
if use6008testData
testCases = test_planetdata();

    j=1;
    for i=1:2:5
        sv1 = testCases(i  , 3:end);                     % Departure Plaent SV
        sv2 = testCases(i+1, 3:end);                     % Arrival   Planet SV
        tof = (testCases(i+1,2) - testCases(i,2))*86400; % Julian Date Diff * dayToS
        
        lam = lambert0rev(sv1,sv2,tof,mu);
        
        resultCases(j,1:3) = lam.vi.';
        resultCases(j,4:6) = lam.vf.';
        j=j+1;
    end


    % Test Data from Professor Kate Davis Canvas: "Lambert-Checks.xlsx" 
    testCaseResult = [4.651443497	 26.08241441	-1.393060432 ...
                      16.79262045	-33.35167484	 1.523021504;
                      13.74077736	 28.83099312	 0.691285008 ...
                     -0.883933069	 -7.983627014	-0.240770598;
                      11.18326152	 -8.90233011 	 0.420697886 ...
                      7.522127215	  4.928368894	-0.474069569];


    mm = round(testCaseResult, 5);
    tm = round(resultCases   , 5);

    disp('Test Cases Check')
    disp('   Pretty much all the values from col. 3:end should be 1')
    disp(mm == tm)
    disp(' ')
    disp(' ')
    
else
%% Test with CPSICE Compare with David Eagle's Lambert Alg.
 
    % Lambert Sample Problem
    depBdy = '3';
    arrBdy = '4';
    etp    = cspice_str2et( {'Jul 30, 2020', 'Feb 18, 2021'} );

    [ctr_bdy] = mice_bodc2n(0);
    tofLc = etp(2) - etp(1);
    sv1 = mice_spkezr(depBdy, etp(1), 'J2000', 'NONE', ctr_bdy.name ).state;
    sv2 = mice_spkezr(arrBdy, etp(2), 'J2000', 'NONE', ctr_bdy.name ).state;
    [vito, vfto] = glambert(pcd.Sun.mu, sv1, sv2, tofLc, 0);

    out = lambert0rev(sv1,sv2,tofLc,mu);

    disp(' ')
    disp(' ')
    disp(' ')
    disp(' ')
    disp('David Eagle on top, lambert0rev() on bottom')
    disp(' ')
    disp(' ')
    vito
    out.vi
    disp(' ')
    vfto
    out.vf
    disp(' ')
    disp(' ')
    disp(' ')
    disp(' ')
    
end


