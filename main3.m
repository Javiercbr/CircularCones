clear all, clc, close all, reset(gpuDevice)


% config file

% fake recosntruction
% rec_fake = true
% data_stored = false
% 
% % true recosntruction with stored projections
% rec_fake = false
% data_stored = true
% 
% % true recosntructions and calculates new data
% rec_fake = false
% data_stored = false

config = load_config('config/config_base.json');
if config.gpu.useGPU
    gpuDevice(config.gpu.device);
end

N = config.sampling.N;
show_legend = config.enable_legend;
% filename = fullfile(config.paths.input, ['phSL3D' num2str(N) '.mat']);
% load(filename);
 
p = phantom3d(N);
p = p.*single(p>0);

M = config.sampling.M;
R = config.geometry.R;
Ktheta = config.sampling.Ktheta;
Ntheta  = Ktheta*M +1;



theta = make_array(0, 2*pi, Ntheta, config, true);





    if config.data_stored == false
%         xM = R/(2*sqrt(2) + 0.01);
%         x = make_array(-xM, xM, N, config);
%         y = make_array( 0.02, 2*xM, N, config);
%         z = y;

        yM = R/(sqrt(3)) ;
        ym = R/50;
        L = yM - ym;
        x = make_array(-L/2, L/2, N, config);
        y = make_array( ym, yM, N, config);
        z = y;

        [rMax] = get_max_r(theta, x, y, z, R);

        Nr =  config.sampling.Nr;
        Ks =  config.sampling.Ks;
        Kpsi00 =  config.sampling.Kpsi00;
        Kvarphi =  config.sampling.Kvarphi;

        if config.rec_fake == false

        J = config.sampling.J; 
        Nvarphi =Kvarphi*M +1;
        Ns =Ks*M +1;
        Npsi = Kpsi00*M +1;
        [p, theta, R, Nr, x, y, z, rMax, Ns, Nvarphi, Npsi] = type_and_device2(p, theta, R, Nr, x, y, z, rMax, Ns, Nvarphi, Npsi, config);



        outputFolder = string(config.paths.input);
        [gz0, varphi, s] = CircularCones_z0_v5(p,theta, Nvarphi, Ns, Nr, Npsi, R, x, y, z, rMax, J, show_legend );
        outputFile = fullfile(outputFolder, "gz0.mat");
        save(outputFile, 'gz0', 'theta', 'varphi', 's', 'R', 'x', 'y', 'z');

        [gy0, varphi, s] = CircularCones_y0_v5(p,theta, Nvarphi, Ns, Nr, Npsi, R, x, y, z, rMax, J, show_legend );
        outputFile = fullfile(outputFolder, "gy0.mat");
        save(outputFile, 'gy0', 'theta', 'varphi', 's', 'R', 'x', 'y', 'z');
    end


end



config = load_config('config/config_base.json');
N = config.sampling.N;

filename = fullfile(config.paths.input, ['phSL3D' num2str(N) '.mat']);
load(filename);
 

M = config.sampling.M;
R = config.geometry.R;

Kphi = config.sampling.Kphi;
Kpsi = config.sampling.Kpsi;


% variables para el W
Nphi = Kphi*M-1; % 2* M+1;
Npsi = Kpsi*M+2; %  azimuth integration variable


phi = make_array(   0, 2*pi, Nphi, config, true);
psi = make_array(pi/4, pi/2, Npsi, config, true);


[p, theta,  phi, psi, R ] = type_and_device2(p, theta,  phi, psi, R, config );


if config.rec_fake == false
filename = fullfile(config.paths.input, ['gz0.mat']);
load(filename);
[ gf_z0 ] = filtrar4( varphi, theta, s, gz0 );

%gf_z0(:,:,1:512) =0;


[ varphi, theta, s, gf_z0 ] = type_and_device2(varphi, theta, s, gf_z0, config );
[Pfz0] =  CalcularP_z0_radon_clasica2(  theta,  varphi, s, gf_z0,  psi,  phi, show_legend);
Pfz0 = gather(Pfz0);
Pfz0(Pfz0<0)=0;
else
[Pfz0] = CalcularP_z0(p, theta,  psi, phi, R, Nr, x, y, z, rMax);

end

clear gz0;




%%
clc
% nuevas dimensiones
Kalpha  = config.sampling.Kalpha;
Nalpha = Kalpha*M+1;


Npsi0  = config.sampling.Npsi0; %67;
Nphi0  = config.sampling.Nphi0; %101;

 


phi0 = make_array(   pi,   2*pi, Nphi0, config);
psi0 = make_array( pi/4, 3*pi/4, Npsi0, config);

[Fz1]  = full_plane_integration_polar_z0( x, y, z, psi, phi,  Pfz0, theta,  Nalpha, R,  1,  psi0, phi0,  show_legend);
Fz1 = gather(Fz1);


phi0 = make_array(   0, pi, Nphi0, config);
[Fz2]  = full_plane_integration_polar_z0( x, y, z, psi, phi,  Pfz0, theta,  Nalpha, R, -1,  psi0, phi0,  show_legend);
Fz2 = gather(Fz2);


clear Pfz0
%%
clc


phi = make_array(   0, pi, Nphi, config, true);
psi = make_array(0, pi/4, Npsi,  config, true);

if config.rec_fake == false
filename = fullfile(config.paths.input, ['gy0.mat']);
load(filename);
[ gf_y0 ] = filtrar4( varphi, theta, s, gy0 );
%gf_y0(:,:,1:512) =0;

[ varphi, theta, s, gf_y0 ] = type_and_device2(varphi, theta, s, gf_y0, config );
[Pfy0] =  CalcularP_y0_radon_clasica2( theta,  varphi, s, gf_y0,  psi,  phi, show_legend);
Pfy0 = gather(Pfy0);
Pfy0 (Pfy0<0)=0;
else


[Pfy0] = CalcularP_y0(p, theta,  psi, phi, R, Nr, x, y, z, rMax);
end


phi0 = make_array(   0, 2*pi, Nphi0, config);
psi0 = make_array(   0, pi/4, Npsi0, config);
[Fy1 ] = full_plane_integration_polar_y0( x, y, z, psi, phi,  Pfy0, theta,  Nalpha, R,  1 ,  psi0, phi0,  show_legend);
Fy1=gather(Fy1);


psi0 = make_array(3*pi/4, pi, Npsi0, config);
[Fy2 ] = full_plane_integration_polar_y0( x, y, z, psi, phi,  Pfy0, theta,  Nalpha, R,  1 ,  psi0, phi0,  show_legend);
Fy2=gather(Fy2);
clear Pfy0


%%
clc

F =  Fz1 + Fz2  + Fy1 + Fy2 ; %+ Fy3 + Fy4;

dx = x(2) - x(1);
dy = y(2) - y(1);
dz = z(2) - z(1);


f_ = -6*del2( F  , dx, dy, dz)/(8*pi^2);

%%

outputFolder = string(config.paths.output);

outputFile = fullfile(outputFolder, "f_.mat");
save(outputFile, 'f_');

mostrar_muchos('gray','3',f_, p);

