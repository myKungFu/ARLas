function [epsilon] = findOptimalBeta(obj)
    obj.calculate;
    if obj.A<= 1 || obj.B<=1
        %warning('Alpha and Beta cannot be <=1 in the current setup!')
        epsilon = 1e10;
    end
    epsilon = obj.epsilon;
end