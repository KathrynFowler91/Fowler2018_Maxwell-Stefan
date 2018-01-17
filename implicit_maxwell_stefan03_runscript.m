% script to run the implicit_maxwell_stefan03 model

%--------------------------------------------------------------------------

fickian = false;            % true==fickian / false==MS
diffusion = 'darken';       % 'darken' / 2 component 'vignes'

n_components = 2;
dt = 10.^(-9);             % time-step
ip = 1e2;                   % number of grid points
ntm = 1e3;                  % number of time steps

DSelf = 1e-18;
aerosol = {'Sucrose'}; % see MS variable switch (only for 2 components)
CNum = 6;

r0 = 0;                     % lowest radius edge
rN = 1.2e-7;                % upper radius edge
R = 1e-7;                   % initial radius of aerosol
T = 293.15;                 % kelvin (20degreeC)

% initial particle conditions - to change initial concentration profile see below
Xw_init = 0.1;              % initial water mole fraction
Xw_shell = zeros(ntm,1);
Xw_shell(:) = 0.8;

%--------------------------------------------------------------------------

[u_save, RN] = implicit_maxwell_stefan03_updated(...
    fickian,...
    diffusion,...
    n_components,...
    dt,...
    ip,...
    ntm,...
    rN,...
    R,...
    aerosol{1},...
    CNum,...
    DSelf,...
    Xw_init,...
    Xw_shell,...
    T);

%--------------------------------------------------------------------------

% user input

% if nargin > 0 fickian=varargin{1}; end
% if nargin > 1 diffusion=varargin{2}; end
% if nargin > 2 n_components=varargin{3}; end
% if nargin > 3 dt=varargin{4}; end
% if nargin > 4 ip=varargin{5}; end
% if nargin > 5 ntm=varargin{6}; end
% if nargin > 6 rN=varargin{7}; end
% if nargin > 7 R=varargin{8}; end
% if nargin > 8 aerosol=varargin{9}; end
% if nargin > 9 CNum=varargin{10}; end
% if nargin > 10 DSelf=varargin{11}; end
% if nargin > 11 Xw_init=varargin{12}; end
% if nargin > 12 Xw_shell=varargin{13}; end
% if nargin > 13 T=varargin{14}; end

%--------------------------------------------------------------------------