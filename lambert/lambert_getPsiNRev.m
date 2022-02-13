function [psiL, psiU] = lambert_getPsiNRev(n)
%LAMBERTPSI Compute Bounds on Psi given 'n>0' Number of Revolutions
%
%   Assumptions/Warnings:
%   	1. n>0 (n=1 or more)
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Inputs:
%       1. n     [1x1]        Number of Revolutions
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Output: out structure containing fields:
%       1. psiL  [1x1]        Lower Bound of Psi
%       2. psiU  [1x1]        Upper Bound of Psi
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   References:
%       1. ASEN6008 "Lambert's Problem" Professor Kate Davis
%       "Lambert-Handout-3.pdf"
%       2. Bate, Muller, and White "Fundamentals of Astrodynamics"
%       3. Vallado "Fundamentals of Astrodynamics and Applications"
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
    psiL = 4*(n^2)*(pi^2);
    psiU = 4*((n+1)^2)*(pi^2);
end