%% SacOverShootRemover.m
% Use P1 trace as the genuine saccade trace, and we aim to create a map between the "eye trace" computed from P4-P1 to the P1 trace, 
% so that for any given saccadic "eye trace", the map will generate the correct saccadic trace

%% prepare data set
dataFolder = 'F:/SmoothPursuit/DDPI-MK2';
