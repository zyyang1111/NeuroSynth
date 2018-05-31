function [E] = getL2Energy(L, type)

if type == 'r'
    E = 0.05537 * L + 7.069;
elseif type == 'w'
    E = 0.1138 * L + 7.939; 
else
    fprintf('Error: unknown activity of L2\n');
    E = 0;
end
    