classdef FilmStructure
    % FilmStructure Class summary:
    % The class to form the layered structure of SAW calculation. Used to
    % define elastic constants, thickness, density, and position in the
    % structure.
    
    properties (Access = private)
        layers                      % layers in the structure
        thickness                   % thickness of each layer
        names                       % name of each layer
        layerNumber                 % number of layers
    end
    
    methods (Access = public)
        
        function obj = FilmStructure()
            % Initialize all properties with empty
            
            obj.layers = {};        % cell array to store each layer as a Material Object with order
            obj.thickness = {};     % cell array to store thickness of each layer with order
            obj.names = {};         % name of each layer
            obj.layerNumber = 0;    % total number of layers
        end
        
        function obj = AddLayer(obj, material, thick, orientation)
            % Add a layer intothe structure, with specified material,
            % thickness, and orientation. Material is defined by the
            % Material Class.
            %
            % Orientation: two Eular angle to specify the surface
            % orientation of the layer
            
            alpha = orientation(1);
            beta = orientation(2);
            rotatMat = [[cos(alpha)*cos(beta)  cos(alpha)*sin(beta)  -sin(alpha)];
                        [-sin(beta)            cos(beta)             0];
                        [sin(alpha)*cos(beta)  sin(alpha)*sin(beta)  cos(alpha)]];
            rotatElastMat = FilmStructure.transform(material.GetElastMat(), rotatMat);
            thisLayer = Material();
            thisLayer = thisLayer.SetDensity(material.GetDensity());
            thisLayer = thisLayer.SetElastMat(rotatElastMat);
            thisLayer = thisLayer.SetName(material.GetName());
            obj.layers = [obj.layers, {thisLayer}];
            obj.thickness = [obj.thickness, {thick}];
            obj.names = [obj.names, {material.GetName()}];
            obj.layerNumber = obj.layerNumber + 1;
        end
        
        function layerNameList = GetNames(obj)
            layerNameList = obj.names;
        end
        
        function layers = GetLayers(obj)
            layers = obj.layers;
        end
        
        function num = GetNumLayers(obj)
            num = obj.layerNumber;
        end
        
        function thick = GetThickList(obj)
            thick = obj.thickness;
        end
    end
    
    methods (Static)
        
        function otr = transform(itr,tmx)
            % FUNCTION
            % otr = transform(itr,tmx)
            %
            % DESCRIPTION
            % transform 3D-tensor (Euclidean or Cartesion tensor) of any order (>0) to another coordinate system
            %
            % PARAMETERS
            % otr = output tensor, after transformation; has the same dimensions as the input tensor
            % itr = input tensor, before transformation; should be a 3-element vector, a 3x3 matrix, or a 3x3x3x... multidimensional array, each dimension containing 3 elements
            % tmx = transformation matrix, 3x3 matrix that contains the direction cosines between the old and the new coordinate system
            %
            ne = numel(itr);                % number of tensor elements
            nd = ndims(itr);                % number of tensor dimensions, i.e. order of tensor
            if (ne==3), nd = 1; end         % order of tensor is 1 in case of a 3x1 or 1x3 vector

            otr = itr;                      % create output tensor
            otr(:) = 0;                     % fill output tensor with zeros; this way a symbolic tensor remains symbolic

            iie = zeros(nd,1);              % initialise vector with indices of input tensor element
            ioe = zeros(nd,1);              % initialise vector with indices of output tensor element
            cne = cumprod(3*ones(nd,1))/3;  % vector with cumulative number of elements for each dimension (divided by three)

            for oe = 1:ne,                                 % loop over all output elements
               ioe = mod(floor((oe-1)./cne),3)+1;          % calculate indices of current output tensor element
               for ie = 1:ne,                              % loop over all input elements
                  pmx = 1;                                 % initialise product of transformation matrices
                  iie = mod(floor((ie-1)./cne),3)+1;       % calculate indices of current input tensor element
                  for id = 1:nd,                           % loop over all dimensions
                     pmx = pmx * tmx( ioe(id), iie(id) );  % create product of transformation matrices
                  end
                  otr(oe) = otr(oe) + pmx * itr(ie);       % add product of transformation matrices and input tensor element to output tensor element
               end
            end
        end
        
    end
end