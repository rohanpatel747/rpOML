function jd = getJulianDate(d)
%GETJULIANDATE Get Julian Date using built-in MATLAB function
%
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. d [string] Date Character format: 'dd-Mmm-yyyy'
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. jd  [1x1]  Julian Date from MATLAB
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
jd = juliandate(d,'dd-mm-yyyy');
end

