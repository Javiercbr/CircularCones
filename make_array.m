function t = make_array(t_m, t_M, Nt,   config, extra)

    % Default value for the optional argument
    if nargin < 5
        extra = false;
    end

    [t_m, t_M, Nt] = type_and_device2(t_m, t_M, Nt, config);

    dt = (t_M - t_m) / (Nt - 1);

    if extra
        t = t_m : dt : (t_M + dt);
    else
        t = t_m : dt : t_M;
    end

end