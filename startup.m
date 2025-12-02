%% MatLAB startup script
%
% This script defines paths and default settings for MINDI. Edit as
% required.
%
% This script will be run automatically if, under windows, the MatLAB icon
% properties 'Start in' field is set to the scripts directory where this
% file is located.  Under UNIX using ~/matlab as your scripts directory will
% achieve the same purpose.  Alternatively you may change to your scripts
% directory from the MatLAB prompt and type startup.


%% EDIT THIS PATH FOR YOUR SETUP
path('/Users/chartat1/itsi/nebula/utils/', path)
path('./utils/', path)
path('./ais/', path)
path('./dx_analysis/', path)



%% Set plotting defaults (edit as required)
S = get(0,'ScreenSize');
set(0,'DefaultFigurePosition',[8,82,S(3:4)/2]);    % Figure position and size
set(0,'DefaultFigureColor',[1,1,1]);               % Figure background color
set(0,'DefaultTextColor',[0,0,0]);                 % Color for text and titles
set(0,'DefaultAxesColor',[1,1,1]);              % Color fill for plot
set(0,'DefaultAxesYColor',[0,0,0]);                % Color of y axis text and ticks
set(0,'DefaultAxesXColor',[0,0,0]);                % Color of x axis text and ticks
set(0,'DefaultAxesZColor',[0,0,0]);                % Color of z axis text and ticks
set(0,'DefaultTextFontSize',14)                    % Font size of titles
set(0,'DefaultAxesFontSize',14)                    % Font size of axes
set(0,'DefaultAxesFontName','arial');              % Axis Font
set(0,'DefaultTextFontName','arial');              % Axis Font



%% envprop toolkit
path('/Users/chartat1/envprop-matlab/', path)
path('/Users/chartat1/envprop-matlab/mfiles/trh5read/', path)
path('/Users/chartat1/envprop-matlab/mfiles/cmapfun/', path)
path('/Users/chartat1/envprop-matlab/mfiles/geofun/', path)

%% pharlap stuff
path('/Users/chartat1/pharlap/',path);
path('/Users/chartat1/pharlap/src/matlab/',path);
path('/Users/chartat1/pharlap/mex/',path);
setenv('DIR_MODELS_REF_DAT', '/Users/chartat1/pharlap/dat/')


%% Startup output (please leave)
V = ver('mindi');
fprintf('_______________________________________________________________________\n\n');
fprintf(' %s, version %s\n',V.Name, V.Version);
fprintf(' %s %s\n',V.Release, V.Date);
fprintf('_______________________________________________________________________\n\n');
fprintf('   Running startup from   : %s\n',which('startup.m'));
fprintf('   Working directory      : %s\n',pwd);
fprintf('   Computer type and host : %s %s\n\n',computer,getenv('HOSTNAME'));
fprintf('_______________________________________________________________________\n');

