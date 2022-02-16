function D = getDatePastDate(d,t)
%GETDATEPASTDATE Get Julian Date using built-in MATLAB function
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. d [string] Date Character format: 'dd-mm-yyyy'
%       2. t [1x1]    Number of Days Past the Date
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. D  [1x1]  Julian Date from MATLAB
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
    jd = juliandate(d,'dd-mm-yyyy');
    jdn= jd+t;
    D = datetime(jdn,'convertfrom','juliandate','Format','dd-MM-yyyy');
end

