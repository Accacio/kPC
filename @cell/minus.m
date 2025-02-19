function C = minus(A,B)
% Define .* operator to cells
%   Detailed explanation goes here
C = cellfun(@minus,A,B,'UniformOutput',false);
end
