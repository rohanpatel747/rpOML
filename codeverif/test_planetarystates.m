%% rpOML + ASEN6008 Planetary States Check File 
%  C: 03FEB22

clear; clc; format long g; rpOMLstart();


%% Function:
%       ./getPlanetaryState()
        help getStatePlanet



%% Meeus Method Test Case
% Data from Professor Kate Davis (Canvas: "LambertChecks.xlsx")
testCases = test_planetdata();

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