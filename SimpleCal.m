%% Interface for calculation of SAW velocity
%
% This interface is used to calculate SAW velocity of layered structure
% with any number of layers.
% 
% Usage:
%
% First step: define materials using the following grammer
%
%   m1 = Material(name_of_material, elastic_constants, density, symmetry)
%     
%     name_of_material: string as the name of the material
%     density: unit g/cm^3
%     elastic_constants: unit GPa
%     
%     For isotropic material: 
%                   elastic_constants has format: [c11, c44]
%                   symmetry: 'iso'
%     
%     For cubic material:
%                   elastic_constants has format: [c11, c12, c44]
%                   symmetry: 'cubic' (001 top plane)
%                             'texture111', 
%                             'texture110'
%     
%     For hexagonal material: 
%                   elastic_constants has format: [c11, c12, c13, c33, c44]
%                   symmetry: 'hex' (perpendicular c-axis)
% 
% 
% Second step: define the layered structure using following grammer
%   
%   sample = FilmStructure();
%   sample = sample.AddLayer(m1, thickness, in_plane_orientation);
%   
%       m1: the material defined before.
%       thickness: thickness of the current layer, unit um.
%       in_plane_orientation: in-plane Eular angle. If it's in-plane
%                   isotropic, [0,0] is used. Unit degree.
% 
%   Add the layers in the order of from top to substrate. The last layer
%   will be assume as substrate so the thickness doesn't matter.
%%
clear all

%% Define the layered structure.

Si = Material('Si', [167.4, 65.2, 79.6], 2.33, 'cubic');
Al = Material('Al', [107, 60.8, 28.3], 2.7, 'texture111');
Polymer = Material('Polymer', [5.8, 1.3], 1.03, 'iso');
Cu = Material('Cu', [169, 75.3, 122], 8.96, 'texture111');
Al2O3 = Material('Al2O3', [497, 163, 116, 501, 147], 3.95, 'hex');
V = Material('V', [230, 120, 43.1], 6.0, 'texture110');

% Add the substrate last
sample = FilmStructure();
sample = sample.AddLayer(V, 0.5, [0,0]);
% sample = sample.AddLayer(Polymer, 0.173, [0,0]);
sample = sample.AddLayer(Si, 1000, [0,0]);
% sample = sample.AddLayer(Cu, 0.2, [0,0]);
% sample = sample.AddLayer(Al2O3, 1000, [0,0]);


%% Change the following if needed

theta = 0 * pi/180; % propagation direction of SAW on surface

% Range of frequency and wave vector for Green Function plot
k_resolu = 128; 
kMax = 10; % Upper limit of the k in calculation

omega_resolu = 256;
omega_Max = 80; % Upper limit of the omega in calculation

omega = linspace(1, omega_Max, omega_resolu);
k = linspace(0, kMax, k_resolu);

waveLength = 0.7; % Wavelength of the SAW that is measured. Unit: um.
kIn = 2 * pi / waveLength; % k vector that is measured. Unit: 1/um.

%% Calculation flow. Don't need to change unless really necessary.

% Calculate And Plot Green Function
dispersionCode = input('Do you want to calculate dispersion curve?(Y/N)\n','s');
if dispersionCode=='Y'|| dispersionCode=='y'
    display('Draw dispersion curve')
    [GreenFun, omegaCen] = SawUtility.SawDispersion(sample, omega, k, theta, kIn, 1, [0, 0.1]);
end
% Calculate velocity accurate
velocityCode = input('Do you want to calculate velocity?(Y/N)\n','s');
if velocityCode == 'Y' || velocityCode == 'y'
    [velocity, angularFreq] = SawUtility.AccurateSAW(sample, kIn, theta);
    disp(['The velocity is: ', num2str(velocity), ...
        'nm/ps, at bangular frequency of: ', num2str(angularFreq), 's^-1']);
end