function varargout = type_and_device2(varargin)

    config = varargin{end};
    varargin(end) = [];

    n = nargin - 1;
    varargout = cell(1,n);

    for k = 1:n

        tmp = single(varargin{k});

        if config.gpu.useGPU
            tmp = gpuArray(tmp);
        end

        varargout{k} = tmp;

    end

end