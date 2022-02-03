%% rpOML + ASEN6008 Planetary States Check File 
%  C: 03FEB22

clear; clc; format long g; rpOMLstart();


%% Function:
%       ./getPlanetaryState()
        help getStatePlanet



%% Meeus Method Test Case
% Data from Professor Kate Davis (Canvas: "LambertChecks.xlsx")
%
% Format:
%   Planet ID   Julian Date     Inertial State: X Y Z Vx Vy Vz in km and km/s
testCases = [
    2 2455610   -88002509.16	-62680223.13	4220331.525	20.0705936	-28.68982987	-1.551291815;
    3 2455450   147084764.9	-32521189.65	467.1900914	5.94623924	28.97464121	-0.000715915;
    4 2456300   170145121.3	-117637192.8	-6642044.272	14.70149986	22.00292904	0.100109562;
    5 2457500   -803451694.7	121525767.1	17465211.78	-2.110465959	-12.31199244	0.098198408;
    6 2455940   -1334047119	-571391392.8	63087187.14	3.265660977	-8.899950822	0.02505182;
    8 2461940   4446562425	484989501.5	-111833872.5	-0.627466452	5.427326309	-0.097899482;
    3 2460545   130423562.1	-76679031.85	3624.816561	14.61294123	25.56747613	-0.001503446;
    2 2460919   19195371.67	106029328.4	348953.802	-34.57913611	6.064190776	2.078550651];

% Meeus Method Check Algorithm
for i=1:length(testCases)
    st = getStatePlanet(testCases(i,1), testCases(i,2));
    resultCases(i,1)   = testCases(i,1);
    resultCases(i,2)   = testCases(i,2);
    resultCases(i,3:8) = st.x;
end

mm = round(resultCases, 2);
tm = round(testCases  , 2);

disp('Test Cases Check')
disp('   Pretty much all the values from col. 3:end should be 1')
disp(mm == tm)
disp(' ')
disp(' ')