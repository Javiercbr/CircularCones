
config = load_config('config/config_base.json');
R = 1;
N = 32;

        xM = R/(2*sqrt(3)) ;
        x = make_array(-xM, xM, N, config);
        y = make_array( 0, 2*xM, N, config);
        z = y;