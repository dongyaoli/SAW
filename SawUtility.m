classdef SawUtility
    
    methods
        function obj = SawUtility()
        end
    end
    methods (Static)
        
        function [GreenFun, omegaCenZ] = SawDispersion(structure, omegaRange, kRange, theta, kIn, method, clims)
            k1Range = kRange .* cos(theta); 
            k2Range = kRange .* sin(theta);
            M = length(omegaRange);
            N = length(kRange);
            GreenFun = zeros(3,M,N);
            parfor i = 1:N
                for j = 1:M
                    wave = SurfAcoustWave(structure, omegaRange(j), k1Range(i), k2Range(i), method);
                    GreenFun(:,j,i) = wave.SurfGreenFunc();
                end
            end
            verticalGreen = abs(squeeze(GreenFun(3,:,:)));           
            SawUtility.PlotGreen(verticalGreen, kRange, omegaRange, clims, 'Vertical Green Function');
            omegaCenZ = SawUtility.PlotPeak(verticalGreen, kRange, omegaRange, kIn);
            display(strcat('The interested frequency in z direction is : ', num2str(omegaCenZ)));
            
%             yGreen = abs(squeeze(GreenFun(2,:,:)));
%             SawUtility.PlotGreen(yGreen, kRange, omegaRange, clims, 'y Green Function');
%             omegaCenY = SawUtility.PlotPeak(yGreen, kRange, omegaRange, kIn);
%             display(strcat('The interested frequency in y direction is : ', num2str(omegaCenY)));
        end
        
        function PlotGreen(twoDMap, kRange, omegaRange, clims, ti)
            figure
            imagesc(kRange*1e6, omegaRange*1e9, twoDMap, clims);
            set(gca, 'YDir', 'normal'); 
            xlabel('Wave Vector (1/m)'); 
            ylabel('Angular Frequency (1/s)');
            title(ti);
        end
        
        function omegaCen = PlotPeak(twoDMap, kRange, omegaRange, kIn)
            % Find the peak in a 2D color map
            M = length(omegaRange);
            N = length(kRange);
            map = diff(twoDMap);
            frequency = zeros(N);
            for i = 1:N
                m = 1;
                for j = 2:M-3
                    if( (map(j,i)*map(j+1,i)<0) && ((map(j,i)-map(j+1,i))>0.01*max(map(:,i))) )
                        frequency(m,i) = real(omegaRange(j));
                        m = m + 1;
                    end
                end
            end
            frequency = frequency';
            figure
            plot(abs(kRange*1e6), frequency*1e9, 'o')
            % Frequency matrix contains information on dispersion relationship
            % the column corresponds to K, the data is the frequency where oscillation
            % occurs    
            mode = 1;
            interestK = floor(kIn/max(kRange)*N);
            omegaCen = frequency(interestK, mode);           
        end
        
        function [velocity, frequency] = AccurateSAW(filmStructure, kIn, theta)
            DoAgain = 1;
            while DoAgain == 1
                centerFrequency = input('Please enter center frequency\n');
                kInX = kIn * cos(theta); kInY = kIn * sin(theta);
                pointNumber = 1500;
                omegaIn = linspace(0.975*centerFrequency, 1.025*centerFrequency, pointNumber);
                GreenFun = zeros(3,pointNumber);
                parfor i = 1:pointNumber
                    wave = SurfAcoustWave(filmStructure, omegaIn(i), kInX, kInY, 2);
                    GreenFun(:,i) = wave.SurfGreenFunc();
                end
                x = abs(omegaIn)';
                y = abs(squeeze(GreenFun(3,:)));
                figure
                plot(x, y, 'o')
                VsawCode = input('Do you want to do it again?(Y/N)\n','s');
                if VsawCode=='N' || VsawCode=='n'
                    DoAgain = 0;
                end
            end
            fo_ = fitoptions('method','NonlinearLeastSquares','Lower',[-Inf -Inf 0]);
            st_ = [1e-3 centerFrequency 0.5];
            set(fo_,'Startpoint',st_);
            ft_ = fittype('gauss1');
            cf_ = fit(abs(omegaIn)',abs(squeeze(GreenFun(3,:)))',ft_,fo_);
            velocity = cf_.b1./kIn;   
            frequency = (kIn * 10^6) * velocity;
        end
        
        function disp = Displacement(structure, kIn, theta, depth)
            k1 = kIn * cos(theta); k2 = kIn * sin(theta);
            pointNum = length(depth);
            omega = input('Please enter center frequency\n');           
            dispCarte = zeros(3, pointNum);
            parfor i = 1:pointNum
                wave = SurfAcoustWave(structure, omega, k1, k2, 2);
                dispCarte(:,i) = wave.GreenFuncDepth(depth(i));               
            end
            % The coordinates change to: propagation direction becomes x
            disp = zeros(3, pointNum);
            for j = 1:pointNum
                disp(1,j) = dispCarte(1,j) * cos(theta) + dispCarte(2,j) * sin(theta);
                disp(2,j) = dispCarte(1,j) * sin(theta) - dispCarte(2,j) * cos(theta);
                disp(3,j) = dispCarte(3,j);
            end
        end
        
        function PlotDisplacement(U, depth, wavelength)
            x = U(1,:);
            y = U(2,:);
            z = U(3,:);
            ref = max ( max(abs(real(z))), max(abs(imag(z))));
            
            figure
            p1 = plot(depth/wavelength, real(x)/ref, 'Color', 'k', 'LineStyle', '--');
            hold on
            p2 = plot(depth/wavelength, imag(x)/ref, 'Color', 'r', 'LineStyle', ':');
            legend([p1 p2], {'real(// disp)', 'imag(// disp)'})
            xlabel('Depth(Wavelength)')
            ylabel('Relative Amplitude')
            hold off
            
            figure
            p3 = plot(depth/wavelength, real(y)/ref, 'Color', 'k', 'LineStyle', '--');
            hold on
            p4 = plot(depth/wavelength, imag(y)/ref, 'Color', 'r', 'LineStyle', ':');
            legend([p3 p4], {'real(\perp disp)', 'imag(\perp disp)'})
            xlabel('Depth(Wavelength)')
            ylabel('Relative Amplitude')
            hold off
            
            figure
            p5 = plot(depth/wavelength, real(z)/ref, 'Color', 'k', 'LineStyle', '--');
            hold on
            p6 = plot(depth/wavelength, imag(z)/ref, 'Color', 'r', 'LineStyle', ':');
            legend([p5 p6], {'real(z disp)', 'imag(z disp)'})
            xlabel('Depth(Wavelength)')
            ylabel('Relative Amplitude')
            hold off
            
        end
        
        function strain = Strain(structure, kIn, theta, depth)
            k1 = kIn * cos(theta); k2 = kIn * sin(theta);
            pointNum = length(depth);
            omega = input('Please enter center frequency\n');
            strainCarte = zeros(3, 3, pointNum);
            parfor i = 1:pointNum
                wave = SurfAcoustWave(structure, omega, k1, k2, 2);
                strainCarte(:,:,i) = wave.StrainDepth(depth(i));
            end
            % The coordinates change to: propagation direction becomes x
            strain = zeros(3, 3, pointNum);
            rotationMat = [cos(theta) sin(theta) 0; ...
                            -sin(theta) cos(theta) 0; ...
                            0 0 1];
            for j = 1:pointNum
                strain(:,:,j) = FilmStructure.transform(strainCarte(:,:,j), rotationMat);
            end
        end
        
        function PlotStrain(strain, depth, wavelength)
            epsilon11 = squeeze(strain(1,1,:));
            epsilon22 = squeeze(strain(2,2,:));
            epsilon33 = squeeze(strain(3,3,:));
            epsilon12 = squeeze(strain(1,2,:));
            epsilon13 = squeeze(strain(1,3,:));
            epsilon23 = squeeze(strain(2,3,:));
            
            figure
            p1 = plot(depth/wavelength, real(epsilon11), 'Color', 'k', 'LineStyle', '-');
            hold on
            p2 = plot(depth/wavelength, imag(epsilon11), 'Color', 'r', 'LineStyle', '-');
            p3 = plot(depth/wavelength, real(epsilon22), 'Color', 'k', 'LineStyle', '--');
            p4 = plot(depth/wavelength, imag(epsilon22), 'Color', 'r', 'LineStyle', '--');
            p5 = plot(depth/wavelength, real(epsilon33), 'Color', 'k', 'LineStyle', ':');
            p6 = plot(depth/wavelength, imag(epsilon33), 'Color', 'r', 'LineStyle', ':');
            legend([p1 p2 p3 p4 p5 p6], {'real(\epsilon_{11})', 'imag(\epsilon_{11})', 'real(\epsilon_{22})', 'imag(\epsilon_{22})', 'real(\epsilon_{33})', 'imag(\epsilon_{33})'})
            xlabel('Depth(Wavelength)')
            ylabel('Compressive Strain')
            hold off
            
            figure
            p1 = plot(depth/wavelength, real(epsilon12), 'Color', 'k', 'LineStyle', '-');
            hold on
            p2 = plot(depth/wavelength, imag(epsilon12), 'Color', 'r', 'LineStyle', '-');
            p3 = plot(depth/wavelength, real(epsilon13), 'Color', 'k', 'LineStyle', '--');
            p4 = plot(depth/wavelength, imag(epsilon13), 'Color', 'r', 'LineStyle', '--');
            p5 = plot(depth/wavelength, real(epsilon23), 'Color', 'k', 'LineStyle', ':');
            p6 = plot(depth/wavelength, imag(epsilon23), 'Color', 'r', 'LineStyle', ':');
            
            legend([p1 p2 p3 p4 p5 p6], {'real(\epsilon_{12})', 'imag(\epsilon_{12})', 'real(\epsilon_{13})', 'imag(\epsilon_{13})', 'real(\epsilon_{23})', 'imag(\epsilon_{23})'})
            xlabel('Depth(Wavelength)')
            ylabel('Shear Strain')
            hold off
            
        end
        
        function Kinetic_E = KineticEnergy(structure, kIn, theta, depth)
            k1 = kIn * cos(theta); k2 = kIn * sin(theta);
            pointNum = length(depth);
            omega = input('Please enter center frequency\n');
            Kinetic_E = zeros(1, pointNum);
            parfor i = 1:pointNum
                wave = SurfAcoustWave(structure, omega, k1, k2, 2);
                Kinetic_E(i) = wave.KineticDepth(depth(i));
            end
        end
        
        function [energy, loss, strain] = CubicEnergy(structure, kIn, theta, depth, visco)
            k1 = kIn * cos(theta); k2 = kIn * sin(theta);
            pointNum = length(depth);
            omega = input('Please enter center frequency\n');
            energy = zeros(3, pointNum);
            loss = zeros(3, pointNum);
            strain = zeros(3, pointNum);
            parfor i = 1:pointNum
                wave = SurfAcoustWave(structure, omega, k1, k2, 2);
                [energy(:,i), loss(:,i), strain(:,i)] = wave.CubicEnergy(depth(i), visco);
            end
        end
            
        function [total_compress, total_tetra, total_shear] = PlotEnergy(energy, depth, wavelength)
            compress = squeeze(energy(1,:));
            tetraShear = squeeze(energy(2,:));
            normalShear = squeeze(energy(3,:));
            figure
            p1 = plot(depth/wavelength, compress, 'Color', 'k', 'LineStyle', '-');
            hold on
            p2 = plot(depth/wavelength, tetraShear, 'Color', 'r', 'LineStyle', '-');
            p3 = plot(depth/wavelength, normalShear, 'Color', 'b', 'LineStyle', '-');
            p4 = plot(depth/wavelength, tetraShear + normalShear, 'Color', 'r', 'LineStyle', '--');
            legend([p1 p2 p3 p4], {'Compressive', 'Tetra Shear', 'Normal Shear', 'Total Shear'})
            xlabel('Depth(Wavelength)')
            ylabel('Energy')
            hold off
            total_compress = sum(compress);
            total_tetra = sum(tetraShear);
            total_shear = sum(normalShear);
        end
        
        function [loss_compress, loss_tetra, loss_shear] = PlotLoss(loss, depth, wavelength)
            compress = squeeze(loss(1,:));
            tetraShear = squeeze(loss(2,:));
            normalShear = squeeze(loss(3,:));
            figure
            p1 = plot(depth/wavelength, compress, 'Color', 'k', 'LineStyle', '-');
            hold on
            p2 = plot(depth/wavelength, tetraShear, 'Color', 'r', 'LineStyle', '-');
            p3 = plot(depth/wavelength, normalShear, 'Color', 'b', 'LineStyle', '-');
            p4 = plot(depth/wavelength, tetraShear + normalShear, 'Color', 'r', 'LineStyle', '--');
            legend([p1 p2 p3 p4], {'Compressive', 'Tetra Shear', 'Normal Shear', 'Total Shear'})
            xlabel('Depth(Wavelength)')
            ylabel('loss')
            hold off
            loss_compress = sum(compress);
            loss_tetra = sum(tetraShear);
            loss_shear = sum(normalShear);
        end
        
    end
    
end