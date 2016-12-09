clear all

Si = Material('Si', [167.4, 65.2, 79.6], 2.33, 'cubic');
Al = Material('Al', [107, 60.8, 28.3], 2.7, 'texture111');
% MoSe2 = Material('MoSe2', [196.1, 42.3, 9.8, 32.9, 3], 6.96, 'hex');

% Ferrite = Material('NiZnAl', [335.5, 134.2, 49], 5.12, 'cubic');
% MgAlO = Material('MgAl2O4', [286.3, 157.2, 153.4], 3.64, 'cubic');

% Guess elastic constant of (SnSe)(MoSe2)
c11 = 109.4;
c66 = 40.5;
c12 = c11 - 2 * c66;
c13 = 4.5;
SnSeMoSe2 = Material('SnSeMoSe2', [c11, c12, c13, 38.8 * 1.05, 1], 6.57, 'hex');

% Usage: layers.AddLayer(material, thickness, Eular Angle-[alpha, beta])
% Add the substrate last

sample = FilmStructure();
% sample = sample.AddLayer(Si, 0.1, [0,0]);
sample = sample.AddLayer(Al, 0.147, [0,0]);
sample = sample.AddLayer(SnSeMoSe2, 0.0596, [0,0]);
% sample = sample.AddLayer(MoSe2, 0.053, [0,0]);
sample = sample.AddLayer(Si, 1000, [0,0]);
% sample = sample.AddLayer(Ferrite, 0.085, [0,0]);
% sample = sample.AddLayer(MgAlO, 1000, [0,0]);


% Range of frequency and wave vector for Green Function plot
N = 100; 
M = 256;
kMax = 10; 
k = linspace(0, kMax, N);
omega = linspace(1, 60, M);
theta = 45 * pi/180; % propagation direction of SAW on surface
waveLength = 0.7;
kIn = 2 * pi / waveLength;

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
        'nm/ps, at angular frequency of: ', num2str(angularFreq), 's^-1']);
end

% Calculate displacement
dispCode = input('Do you want to calculate displacement?(Y/N)\n','s');
if dispCode == 'Y' || dispCode == 'y'
    pointNum = 1000;
    depth = (0:1/(pointNum-1):1) * waveLength * 3;
    U = SawUtility.Displacement(sample, kIn, theta, depth);
    SawUtility.PlotDisplacement(U, depth, waveLength);
end

% Calculate strain
strainCode = input('Do you want to calculate strain?(Y/N)\n','s');
if strainCode == 'Y' || strainCode == 'y'
    pointNum = 1000;
    depth = (0:1/(pointNum-1):1) * waveLength * 3;
    strain = SawUtility.Strain(sample, kIn, theta, depth);
    SawUtility.PlotStrain(strain, depth, waveLength);
end

% Calculate density of strain energy
energyCode = input('Do you want to calculate strain energy?(Y/N)\n','s');
% Lamb viscosity
% visco = [1.505, -0.532, 0.553];
% King's viscosity
visco = [5.9, 5.16, 0.62];
if energyCode == 'Y' || energyCode == 'y'
    pointNum = 1000;
    depth = (0:1/(pointNum-1):1) * waveLength * 1.1;
    [energy, loss, strain_power] = SawUtility.CubicEnergy(sample, kIn, theta, depth, visco);
    % Time average of all losses
    % loss = loss / 2;
    % No 0.5 for elastic energy. This elastic energy is not time averaging
    % but amplitude
    [compress, tetra, shear] = SawUtility.PlotEnergy(energy, depth, waveLength);
    [compress_loss, tetra_loss, shear_loss] = SawUtility.PlotLoss(loss, depth, waveLength);
    
    strain_comp = sum(squeeze(strain_power(1,:)));
    strain_tetra = sum(squeeze(strain_power(2,:)));
    strain_shear = sum(squeeze(strain_power(3,:)));
end

kineticCode = input('Do you want to calculate Kinetic energy?(Y/N)\n','s');
if kineticCode == 'Y' || kineticCode == 'y'
    pointNum = 1000;
    depth = (0:1/(pointNum-1):1) * waveLength * 1.1;
    KE_depth = SawUtility.KineticEnergy(sample, kIn, theta, depth);
    kinetic_Energy = sum(KE_depth);
end


