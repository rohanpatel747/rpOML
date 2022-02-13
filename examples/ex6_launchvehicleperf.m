%% rpOML Examples - Launch Vehicle Performance
clear; clc; format long g; rpOMLstart();
tic
%% Overview
% This script is intended to show how to use the launch vehicle performance
% function.

% Credit to Dr. Gregory Lantoine and Joshua Fofrich for the compile data
% for the launch vehicle C3 and mass delivered parameters


%% Problem Setup
departureC3 = 12.5;     % Required C3 in (km2/s2)

showPlot = true;        % Show various LV cruves for mass delivered v. C3
dispData = true;        % Display in Command Window Mass Delivered at C3

lvp = lvperformance(departureC3,showPlot,dispData);

disp(['The Atlas V 551 can deliver: ', num2str(lvp.atlasv551), ' kg'])


