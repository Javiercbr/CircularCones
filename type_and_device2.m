function varargout = type_and_device2(varargin)

n = nargin;
varargout = cell(1,n);

for k = 1:n

    tmp = varargin{k};

    if ~isa(tmp,'single')
        tmp = single(tmp);
    end

    if ~isa(tmp,'gpuArray')
        tmp = gpuArray(tmp);
    end

    varargout{k} = tmp;

end

end