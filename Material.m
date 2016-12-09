classdef Material
    properties (GetAccess = private)
        density
        elastMat
        symmetry
        name
    end
    
    methods
        function obj = Material(varargin)
            if nargin == 4
                obj.name = varargin{1};
                elastCons = varargin{2};
                obj.density = varargin{3};
                obj.symmetry = varargin{4};
                if strcmp(obj.symmetry, 'iso')
                    obj.elastMat = Material.Iso(elastCons);
                elseif strcmp(obj.symmetry, 'cubic')
                    obj.elastMat = Material.Cubic(elastCons);    
                elseif strcmp(obj.symmetry, 'hex')
                    obj.elastMat = Material.Hex(elastCons);
                elseif strcmp(obj.symmetry, 'texture111')
                    obj.elastMat = Material.Texture111(elastCons);
                elseif strcmp(obj.symmetry, 'texture110')
                    obj.elastMat = Material.Texture110(elastCons);
                else
                    error('Wrong Symmetry. (Supportted symmetry: iso, cubic, hex, texture111)') 
                end
            elseif nargin == 0
                obj.density = NaN;
                obj.elastMat = NaN;
                obj.symmetry = NaN;
            else
                error('Wrong number of argument for constructor')
            end
        end
                
        function dens = GetDensity(obj)
            dens = obj.density;
        end
        
        function ec = GetElastMat(obj)
            ec = obj.elastMat;
        end
        
        function sym = GetSymmetry(obj)
            sym = obj.symmetry;
        end
        
        function name = GetName(obj)
            name = obj.name;
        end
        
        function obj = SetElastMat(obj, elastMat)
            obj.elastMat = elastMat;
        end
        
        function obj = SetDensity(obj, density)
            obj.density = density;
        end
        
        function obj = SetSymmetry(obj, sym)
            if strcmp(sym, 'iso') ||  strcmp(sym, 'cubic') || strcmp(sym, 'hex') || strcmp(sym, 'texture111')
                 obj.symmetry = sym;
            else    
                 error('Symmetry not supportted. (Supportted symmetry: iso, cubic, hex, texture111)') 
            end
           
        end
        
        function obj = SetName(obj, name)
            obj.name = name;
        end
    end
    
    methods (Static)
        function ec = Iso(elastCons)
            if length(elastCons) ~= 2
                print 'Wrong number of elastic constants for isotropic crystal'
                return
            end
            c11 = elastCons(1);
            c44 = elastCons(2);           
            c12 = c11 - 2 * c44;
            ec = zeros(3,3,3,3);    
            ec(1,1,1,1) = c11; ec(2,2,2,2) = c11; ec(3,3,3,3) = c11;
            ec(1,1,2,2) = c12; ec(1,1,3,3) = c12; ec(2,2,3,3) = c12;
            ec(2,2,1,1) = c12; ec(3,3,1,1) = c12; ec(3,3,2,2) = c12;

            ec(2,3,2,3) = c44; ec(1,3,1,3) = c44; ec(1,2,1,2) = c44;
            ec(3,2,2,3) = c44; ec(3,1,1,3) = c44; ec(2,1,1,2) = c44;
            ec(2,3,3,2) = c44; ec(1,3,3,1) = c44; ec(1,2,2,1) = c44;
            ec(3,2,3,2) = c44; ec(3,1,3,1) = c44; ec(2,1,2,1) = c44;          
        end
        
        function ec = Cubic(elastCons)
            if length(elastCons) ~= 3
                print 'Wrong number of elastic constants for cubic crystal'
                return
            end
            c11 = elastCons(1);
            c12 = elastCons(2);
            c44 = elastCons(3);
            ec = zeros(3,3,3,3);    
            ec(1,1,1,1) = c11; ec(2,2,2,2) = c11; ec(3,3,3,3) = c11;
            ec(1,1,2,2) = c12; ec(1,1,3,3) = c12; ec(2,2,3,3) = c12;
            ec(2,2,1,1) = c12; ec(3,3,1,1) = c12; ec(3,3,2,2) = c12;
            
            ec(2,3,2,3) = c44; ec(1,3,1,3) = c44; ec(1,2,1,2) = c44;
            ec(3,2,2,3) = c44; ec(3,1,1,3) = c44; ec(2,1,1,2) = c44;
            ec(2,3,3,2) = c44; ec(1,3,3,1) = c44; ec(1,2,2,1) = c44;
            ec(3,2,3,2) = c44; ec(3,1,3,1) = c44; ec(2,1,2,1) = c44;           
        end
        
        function ec = Hex(elastCons)
            if length(elastCons) ~= 5
                print 'Wrong number of elastic constants for hexagonal crystal'
                return
            end
            c11 = elastCons(1);
            c12 = elastCons(2);
            c13 = elastCons(3);
            c33 = elastCons(4);
            c44 = elastCons(5);
            
            c66 = 0.5 * (c11 - c12);
            ec = zeros(3,3,3,3);    
            ec(1,1,1,1) = c11; ec(2,2,2,2) = c11; ec(3,3,3,3) = c33;
            ec(1,1,2,2) = c12; ec(1,1,3,3) = c13; ec(2,2,3,3) = c13;
            ec(2,2,1,1) = c12; ec(3,3,1,1) = c13; ec(3,3,2,2) = c13;    

            ec(2,3,2,3) = c44; ec(1,3,1,3) = c44; ec(1,2,1,2) = c66;
            ec(3,2,2,3) = c44; ec(3,1,1,3) = c44; ec(2,1,1,2) = c66;
            ec(2,3,3,2) = c44; ec(1,3,3,1) = c44; ec(1,2,2,1) = c66;
            ec(3,2,3,2) = c44; ec(3,1,3,1) = c44; ec(2,1,2,1) = c66;           
        end
        
        function ec = Texture111(elastCons)
            if length(elastCons) ~= 3
                print 'Wrong number of elastic constants for cubic crystal'
                return
            end
            c11 = elastCons(1);
            c12 = elastCons(2);
            c44 = elastCons(3);
            
            a11 = (c11+c12+2*c44)/2;
            a33 = (c11+2*c12+4*c44)/3;
            a12 = (c11+5*c12-2*c44)/6;
            a13 = (c11+2*c12-2*c44)/3;
            a44 = (c11-c12+c44)/3;
            ec = Material.Hex([a11, a12, a13, a33, a44]);
        end
        
        function ec = Texture110(elastCons)
            if length(elastCons) ~= 3
                print 'Wrong number of elastic constants for cubic crystal'
                return
            end
            c11 = elastCons(1);
            c12 = elastCons(2);
            c44 = elastCons(3);
            
            a11 = (9 * c11 + 7 * c12 + 14 * c44)/16;
            a33 = (c11 + c12 + 2 * c44)/2;
            a12 = (3 * c11 + 13 * c12 - 6 * c44)/16;
            a13 = (c11 + 3 * c12 - 2 * c44)/4;
            a44 = (c11 - c12 + 2 * c44)/4;
            ec = Material.Hex([a11, a12, a13, a33, a44]);
        end
               
    end
    
end
