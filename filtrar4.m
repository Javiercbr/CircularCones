 function [gf] = filtrar4( varphi, theta, s, g )


Ns = length(s);
Nfft = 2^nextpow2(2*Ns);

 
 
ds = s(2)-s(1);
%f_sample = 1/ds;
Nvarphi = length(varphi);
Ntheta = length(theta);




v = (-Nfft/2:Nfft/2-1)/Nfft;

v=single(v);

%ds = (s(end)-s(1))/Nfft;
ds = s(2)-s(1);
f_sample = 1/ds;
 
rampa = abs(f_sample*v);             % filtro rampa

H = 0.5 + 0.5*cos( 2*pi*v); % filtro de von Hann
filtro = rampa.*H;          % filtro final
 
size(filtro)

filtro = repmat( filtro, Nvarphi, 1 );     % repite el array con el filtro
filtro = squeeze(filtro);

a = zeros(  Nvarphi, Nfft);

gf = zeros(  Ntheta, Nvarphi, Ns);
 
for i = 1 : Ntheta


    waitbar(gather(i/Ntheta) );


    a(: ,1:Ns) = g(i,:,:);
    ga=a; 


 
G=fftshift(fft(ga,[],2),2);
% transformada de Fourier 1D de las filas
% algoritmo FFT

Gf = ifftshift(filtro.*G,2);     % filtra y hace swap
gaf = real(ifft(Gf,[],2));        % antritransforma las filas


% Obs: se puede comentar esta linea para ver el efecto de la
% retroproyeccion sin filtro

% mostrar_muchos('jet','1',rampa,filtro);

gf(i,:,:) = gaf(:,1:Ns);


end