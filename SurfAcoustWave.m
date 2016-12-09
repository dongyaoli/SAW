classdef SurfAcoustWave
    % SurfaceAcoustWave Class summary:
    % This is the main class to do the calculation of SAW
    % 
    % SurfAcoustWave Unit System:
    %   Angular frequency omega: 10^9 s^-1
    %   density: g/cm^3
    %   wave vector k: 1/um = 10^6 m^-1
    %   elastic constants: GPa
    %   Distance: um = 10^-6 m
    
    
    properties (GetAccess = private)
        structure
        layerNumber
        k1
        k2
        omega
        method
        eigenModeList
        k3List
        coeffiList
        scaleConst
        compareTolerance = 1e-8;
    end
    
    methods
        
        function obj = SurfAcoustWave(structure, omega, k1, k2, method)
            obj.structure = structure;
            obj.omega = omega;
            obj.k1 = k1;
            obj.k2 = k2;
            obj.method = method;            
            obj.layerNumber = structure.GetNumLayers();

            obj.eigenModeList = cell(1, obj.layerNumber);
            obj.k3List = cell(1, obj.layerNumber);
            obj.coeffiList = cell(1, obj.layerNumber);
            
            deltaForce = [0;0;1];
            
            listThickness = structure.GetThickList();
            depth = 0;
            listLayerMat = cell(1, obj.layerNumber);
            listLayerMatUp = cell(1, obj.layerNumber);
            listLayerMatDown = cell(1, obj.layerNumber);
            layers = obj.structure.GetLayers();
            % Loop through each layer
            for i = 1:obj.layerNumber
                thisLayer = layers{i};
                elastMat = thisLayer.GetElastMat();
                density = thisLayer.GetDensity;
                thickness = listThickness{i};
                [KK, UU] = SurfAcoustWave.Cristoffel(obj.k1, obj.k2, obj.omega, ...
                    density, elastMat, obj.method);
                if i == obj.layerNumber % The substrate
                    [K, U] = SubstrateFilter(obj, KK, UU);
                else
                    U = UU;
                    K = KK;
                end
                [~, modeNumber] = size(U);
                stressDown = zeros(3,modeNumber); 
                stressUp = zeros(3,modeNumber);
                phaseUp = exp(1i*K* (depth + thickness)).'; 
                phaseDown = exp(1i * K * depth).';
                dispUp = U .* [phaseUp; phaseUp; phaseUp];
                dispDown = U .* [phaseDown; phaseDown; phaseDown];
                
                for m = 1:3
                    for n = 1:modeNumber
                        stressDown(m,n) = sum(sum(elastMat(:,:,m,3) .* (dispDown(:,n) * [1i*obj.k1 1i*obj.k2 1i*K(n)])));
                        stressUp(m,n) = sum(sum(elastMat(:,:,m,3) .* (dispUp(:,n) * [1i*obj.k1 1i*obj.k2 1i*K(n)])));
                    end
                end
                             
                if i == 1                       % surface layer, no displacement considered at surface                   
                    layerMatDown = stressDown;                    
                    if obj.layerNumber == 1         % Also substrate                       
                        layerMat = layerMatDown;
                        listLayerMatDown{1} = layerMatDown;
                    else
                        obj.scaleConst = abs(mean(mean(stressUp))/mean(mean(dispUp)));
                        layerMatUp = [stressUp; obj.scaleConst .* dispUp];
                        reverseUp = eye(6)/layerMatUp;
                        layerMat = layerMatDown * reverseUp;
                        listLayerMatDown{1} = layerMatDown;
                        listLayerMatUp{1} = layerMatUp;
                    end
                else
                    layerMatDown = [stressDown; obj.scaleConst .* dispDown];
                    if i == obj.layerNumber        % Substrate, no up matrix
                        layerMat = layerMatDown;
                        listLayerMatDown{i} = layerMatDown;
                    else
                        layerMatUp = [stressUp; obj.scaleConst .* dispUp];
                        reverse = eye(6)/layerMatUp;
                        layerMat = layerMatDown * reverse;
                        listLayerMatUp{i} = layerMatUp;
                        listLayerMatDown{i} = layerMatDown;
                    end
                end    
                obj.eigenModeList{i} = U;
                obj.k3List{i} = K;
                listLayerMat{i} = layerMat;
                depth = depth + thickness;
            end
            
            % Calculate coefficient of each layer
            S = 1;
            for i = 1:obj.layerNumber
                S = S * listLayerMat{i};
            end
            % delta force
            A_substrate = S\deltaForce;
            if obj.layerNumber == 1
                obj.coeffiList{1} = A_substrate;
            else
                obj.coeffiList{obj.layerNumber} = A_substrate;
                for i = fliplr(1:(obj.layerNumber-1))
                    reverse = eye(6)/listLayerMatUp{i};
                    SS = reverse * listLayerMatDown{i+1};
                    obj.coeffiList{i} = SS * obj.coeffiList{i+1};
                end
            end
        end
        
        function [K, U] = SubstrateFilter(obj, KK, UU)
            j = 1;
            K = zeros(3,1);
            U = zeros(3);
            for i = 1:6
                % The wave attenuates as it goes deepper, or it travels
                % into substrate
                if(abs(imag(KK(i))) > obj.compareTolerance && ...
                    imag(KK(i)) > 0) || (abs(imag(KK(i))) < obj.compareTolerance && real(KK(i)) > 0) 
                   K(j) = KK(i);
                   U(:,j) = UU(:,i);
                   j = j+1;
                end
            end
            [dimen, ~] = size(K);
            if dimen ~= 3
                display('Wrong k3 Choosen!!!')
                return
            end
            if (obj.k1 == 0 && obj.k2 == 0)
                U = eye(3);
            end 
        end
        
        function G = SurfGreenFunc(obj)
            G = obj.eigenModeList{1} * obj.coeffiList{1};            
        end
        
        function G = GreenFuncDepth(obj, depth)
            level = obj.LayerPosition(depth);
            k = obj.k3List{level};
            U = obj.eigenModeList{level};
            A = obj.coeffiList{level};
            G = U * (A .* exp(1i * k * depth));
        end
        
        function strain = StrainDepth(obj, depth)
            level = obj.LayerPosition(depth);
            k3 = obj.k3List{level};
            U = obj.eigenModeList{level};
            A = obj.coeffiList{level};
            % Here the coordinates is not changed to propagation direction.
            strain = zeros(3, 3);
            strain(1,1) = 1i * obj.k1 * (U(1,:) * (A .* exp(1i * k3 * depth)));
            strain(2,2) = 1i * obj.k2 * (U(2,:) * (A .* exp(1i * k3 * depth)));
            strain(3,3) = U(3,:) * (1i * A .* exp(1i * k3 * depth) .* k3);
            strain(1,2) = 0.5 * (1i * obj.k1 * (U(2,:) * (A .* exp(1i * k3 * depth))) ...
                + 1i * obj.k2 * (U(1,:) * (A .* exp(1i * k3 * depth))));
            strain(1,3) = 0.5 * (1i * obj.k1 * (U(3,:) * (A .* exp(1i * k3 * depth))) ...
                + 1i * (U(1,:) * (A .* exp(1i * k3 * depth) .* k3)));
            strain(2,3) = 0.5 * (1i * obj.k2 * (U(3,:) * (A .* exp(1i * k3 * depth))) ...
                + 1i * (U(2,:) * (A .* exp(1i * k3 * depth) .* k3)));
            strain(2,1) = strain(1,2);
            strain(3,1) = strain(1,3);
            strain(3,2) = strain(2,3);
        end
        
        function Kinetic_E = KineticDepth(obj, depth)
            level = obj.LayerPosition(depth);
            k3 = obj.k3List{level};
            U = obj.eigenModeList{level};
            A = obj.coeffiList{level};
            % Here the coordinates is not changed to propagation direction
            u1 = U(1,:) * (A .* exp(1i * k3 * depth));
            u2 = U(2,:) * (A .* exp(1i * k3 * depth));
            u3 = U(3,:) * (A .* exp(1i * k3 * depth));
            layers = obj.structure.GetLayers();
            thisLayer = layers{level};
            density = thisLayer.GetDensity;
            Kinetic_E = 0.5 * obj.omega^2 * density * (abs(u1)^2 + abs(u2)^2 + abs(u3)^2);
        end
        
        function [energy, loss, mode_strain] = CubicEnergy(obj, depth, visco)
            level = obj.LayerPosition(depth);
            layers = obj.structure.GetLayers();
            thisLayer = layers{level};
            elastMat = thisLayer.GetElastMat();
            % The strain is not changed to propagation coordinates here.
            strain = obj.StrainDepth(depth);
            energy = zeros(3,1);
            loss = zeros(3,1);
            mode_strain = zeros(3,1);
            % Bulk modulus
            bulkModulus = (elastMat(1,1,1,1) + 2 * elastMat(1,1,2,2))/3;
            bulkVisco = (visco(1) + 2 * visco(2))/3;
            % Tetragonal shear modulus
            shearTetra = (elastMat(1,1,1,1) - elastMat(1,1,2,2))/2;
            shearTetraVisco = (visco(1) - visco(2))/2;
            % Shear modulus
            shearModu = elastMat(2,3,2,3);
            shearVisco = visco(3);
            
            % Compressive strain energy
            energy(1) = 0.5 * bulkModulus * abs((strain(1,1) + strain(2,2) + strain(3,3)))^2;
            loss(1) = 0.5 * bulkVisco * abs((strain(1,1) + strain(2,2) + strain(3,3)))^2 * obj.omega^2;
            mode_strain(1) = 0.5 * abs((strain(1,1) + strain(2,2) + strain(3,3)))^2;
            
            % Shear strain energy
            energy(2) = 0.5 * shearTetra * (abs(strain(1,1) - strain(2,2))^2 + ...
                1/3 * abs(2 * strain(3,3) - strain(1,1) - strain(2,2))^2);
            energy(3) = 2 * shearModu * (abs(strain(1,2))^2 + abs(strain(2,3))^2 + abs(strain(1,3))^2);
       
            loss(2) = 0.5 * shearTetraVisco * (abs(strain(1,1) - strain(2,2))^2 + ...
                1/3 * abs(2 * strain(3,3) - strain(1,1) - strain(2,2))^2) * obj.omega^2;
            loss(3) = 2 * shearVisco * (abs(strain(1,2))^2 + abs(strain(2,3))^2 + abs(strain(1,3))^2) * obj.omega^2;
            
            mode_strain(2) = 0.5 * (abs(strain(1,1) - strain(2,2))^2 + 1/3 * abs(2 * strain(3,3) - strain(1,1) - strain(2,2))^2);
            mode_strain(3) = 2 * (abs(strain(1,2))^2 + abs(strain(2,3))^2 + abs(strain(1,3))^2);
        end
        
        function level = LayerPosition(obj, depth)
            thickList = obj.structure.GetThickList();
            interface = zeros(1, obj.layerNumber); % includes surface
            for i = 1:obj.layerNumber-1
                interface(i+1) = interface(i) + thickList{i};
            end
            if depth >= interface(obj.layerNumber)
                level = obj.layerNumber;
            else
                level = find(interface>depth, 1) - 1;
            end
        end
    end
    
    methods (Static)
        
        function [K, U] = Cristoffel(kx, ky, omega, density, C, flag)
            if flag == 1
                % Use Peng's method to calculate k3
                a=zeros(3);b=zeros(3);c=zeros(3); d=zeros(3,3,3);
                for m=1:3
                    for n=1:3
                        a(m,n)=-C(m,1,n,1)*kx*kx-(C(m,1,n,2)+C(m,2,n,1))*kx*ky-C(m,2,n,2)*ky*ky;
                        b(m,n)=-1i*(C(m,1,n,3)+C(m,3,n,1))*kx-1i*(C(m,2,n,3)+C(m,3,n,2))*ky; %the minus sign is from the p
                        c(m,n)=C(m,3,n,3);
                        d(:,m,n)=[c(m,n) b(m,n) a(m,n)+density.*(m==n)*omega.^2];
                    end
                end
                Det=conv(d(:,3,3),conv(d(:,1,1),d(:,2,2)))+conv(d(:,3,1),conv(d(:,1,2),d(:,2,3)))+conv(d(:,2,1),conv(d(:,1,3),d(:,3,2)))-conv(d(:,3,1),conv(d(:,1,3),d(:,2,2)))-conv(d(:,3,3),conv(d(:,1,2),d(:,2,1)))-conv(d(:,1,1),conv(d(:,2,3),d(:,3,2)));
                P=roots(Det);
                K=1i*P; % because the definition of P here is -P=1i*k3
            else
                % Symbolic calculation by DYL
                if flag == 2 
                    syms k3
                    % k1 = kx; 
                    % k2 = ky;
                    % k = [k1 k2 k3];
                    k = [kx ky k3];
                    Gama = sym([0 0 0; 0 0 0; 0 0 0]);
                    for m = 1:1:3
                        for n = 1:1:3
                            a = 0;
                            for mm = 1:1:3
                                for nn = 1:1:3
                                    a = a + C(m,mm,n,nn) * k(mm) * k(nn);
                                end
                            end
                            Gama(m,n) = a - density * (m == n) * omega^2;
                        end
                    end
                    K=double(solve (det(Gama),'k3'));
                else
                    display('Please Choose Correct Calculation Algorithm')
                    return
                end
            end

            % Define what is ZERO
            if flag==1
                ZeroTolerance=1e-4;
                OmegaTolerance=5e-3;
            else
                ZeroTolerance=1e-10;
                OmegaTolerance=1e-8;
            end

            % Determine if degenerent happend
            IndenticalRecord = zeros(1,6);
            if mean(abs(real(K))) > ZeroTolerance
                [~,position] = sort(real(K));
                gap = diff(K(position));
                counts = 1;
                for ii = 1:1:5
                    if abs(gap(ii)) < ZeroTolerance
                        IndenticalRecord(counts) = position(ii);
                        IndenticalRecord(counts + 1) = position(ii + 1);
                        counts = counts + 2;
                    end
                end
            else 
                [~,position] = sort(imag(K));
                gap = diff(K(position));
                counts = 1;
                for ii = 1:1:5
                    if abs(gap(ii)) < ZeroTolerance
                        IndenticalRecord(counts) = position(ii);
                        IndenticalRecord(counts + 1) = position(ii+1);
                        counts = counts + 2;
                    end
                end
            end
            if IndenticalRecord(5) ~= 0
                display('Wrong Identity Number Recogonition')
            end

            % Calculate all displacement vector corresponds to all 6 k3
            U = zeros(3,6);  %Store 6 eigenvector for 6 k3 value, with density1*Omega^2 as eigenvalue
            for i = 1:1:6
                % rebuild Cristoffel matrix
                k3 = K(i);
                k = [kx ky k3];  
                Gama = zeros(3);
                for m = 1:1:3
                    for n = 1:1:3
                        for mm = 1:1:3
                            for nn = 1:1:3
                                Gama(m,n) = Gama(m,n) + C(m,mm,n,nn) * k(mm) * k(nn);
                            end
                        end
                    end
                end
                if ~ismember(i,IndenticalRecord)
                    % degenerate doesn't happen for this k3, Use svd to calculate null
                    % vector
                    Temp = Gama - density * omega^2 * eye(3);
                    [~,~,VV] = svd(Temp);
                    U(:,i) = VV(:,end);
                else
                    % degenerate happend for this k3. Use eig to calculate two
                    % orthogonal eigenvector for two identical eigenvalue
                   pair1position = find(IndenticalRecord == i); 
                   if mod(pair1position,2) == 1
                       pair2 = IndenticalRecord(pair1position + 1);
                   else
                       pair2 = IndenticalRecord(pair1position - 1);
                   end
                   % Get rid of Round off error and make Gama really symmetric
                   % And remember use .', which is matrix transpose rather than complex
                   % transpose
                   % Gama2=(Gama.'+Gama)/2;   
                   [VV,DD] = eig(Gama);
                   Freq = sqrt(DD/density);
                   DegenLocation = zeros(1,3);
                   count = 1;
                   for j = 1:1:3
                       if abs(Freq(j,j) - omega) < OmegaTolerance  
                           DegenLocation(count) = j;
                           count = count + 1;
                       end
                   end
                   if DegenLocation(1)==0 || DegenLocation(2)==0 || DegenLocation(3)~=0
                       display('Cannot find correct eigenvalue!!!! Serious mistake!!!!!!!!!!!')
                       break
                   end
                   if abs(U(:,pair2)'*VV(:,DegenLocation(1)))<ZeroTolerance || rank([U(:,pair2);VV(:,DegenLocation(1))])==2
                       U(:,i)=VV(:,DegenLocation(1));
                   else
                       if abs(U(:,pair2)'*VV(:,DegenLocation(2)))<ZeroTolerance || rank([U(:,pair2);VV(:,DegenLocation(2))])==2
                           U(:,i)=VV(:,DegenLocation(2));
                       else
                           %  display('Non-orthogonal eigenvector found')
                           if U(:,pair2)~=VV(:,DegenLocation(1))
                               U(:,i)=VV(:,DegenLocation(1));
                           else
                               if U(:,pair2)~=VV(:,DegenLocation(2))
                                   U(:,i)=VV(:,DegenLocation(2));
                               else
                                   display('Cannot find orthogonal eigenvector!!!!')
                               end
                           end
                       end
                   end
                end
            end
        end
                
    end
    
 
end