function b = broadsearch_sequencename(seq)
% BROADSEARCH_SEQUENCENAME Generates String of Planets for Naming
%
%   Assumptions/Warnings:
%   	[NOTE: This function is not meant to be standalone]
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   Dependencies:
%       1. getPlanetName()
%   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
    b = [];
    for i=1:length(seq)

        if seq(i)<9
            a = getPlanetName(seq(i));
            a = a(1);
        else
            a = '(other)';
        end

        b = [b,a];

    end
end