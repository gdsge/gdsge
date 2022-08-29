% Needs older version of the Symbolic Math toolbox which MATLAB 2017b+ does not
% support. For license issues, cannot distribute MATLAB's library in m files.
% The temporary solution is to compile the sym class into p files and distribute the p
% files instead of the source codes, and list Symbolic Math toolbox as one of the dependencies.
% Licensed MATLAB users can acquire a copy of the older version (<=2017b) of the Symbolic
% Math toolbox and put it in the folder for explicit uses.
pcode('@sym/private');
pcode('@sym');