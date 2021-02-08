classdef Dual
    properties
        value
        deriv
    end
    methods
        function obj = Dual(v, d)
            obj.value = v;
            obj.deriv= d;
        end
        
        
        function z = sin(x)
            z = Dual(sin(x.value), cos(x.value)*x.deriv);
        end
        
        function z = cos(x)
            z = Dual( , );
        end
        
        function z = exp(x)
            z = Dual( , );
        end
        
        function z = sqrt(x)
            z = Dual( , );
        end
        
        function z = mpower(x, p)
           z = Dual( , ); 
        end
        
        function z = mrdivide(x, y)
           z = Dual( , );
        end
        
        function z = plus(x, y)
           z = Dual(x.value + y.value, x.deriv + y.deriv);
        end
    end
    
end
