function Value = returnNameValue(Name,VARARGIN)
%RETURNNAMEVALUE find next element of cell array with name Name and return
% it as corresponding value.
% Input:
%       Name:[(string)] Name
%       VARARGIN: [vector of cell array] variable arguments.
% Output:
%       Value: value after specific name.

validateattributes(VARARGIN,{'cell'},{'vector','nonempty'});
validateattributes(Name,{'char'},{'vector','nonempty'});

for i=1:length(VARARGIN)
    if strcmp(Name,VARARGIN{i})
        Value = VARARGIN{i+1};
        return
    end
end

% there is no Name-Value pair
Value = [];

end

