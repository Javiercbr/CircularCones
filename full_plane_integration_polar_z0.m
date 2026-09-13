function [F] = full_plane_integration_polar_z0( x, y, z, psi, phi,  Pz0, theta,  Nalpha , R, signo,  psi0, phi0, show_legend )

 

dalpha = 2*pi/Nalpha;
alpha = 0 : dalpha : 2*pi - dalpha;

dpsi0 = psi0(2) - psi0(1);
dphi0 = phi0(2) - phi0(1);


% ATENCION
% varphi es la vsariable de la transformada (angulo del eje del cono)
% alpha es el angulo de la integracion polar sobre la esfera en los n
% ortogonales
% psi y phi parametrizan el vector n
% beta es el angulo de la transformada de Radon clasica

%%
 
[z_,  psi_, phi_, alpha_ ]= ndgrid(z,   psi0, phi0, alpha);

Nx = length(x);
Ny = length(y);
Nz = length(z);
Npsi = length(psi0);
Nphi = length(phi0);


F = zeros(Nx,Ny, Nz, Npsi, Nphi);



barra = waitbar(0,'Planar integration');




for i = 1 : Nx

    waitbar(gather(i/Nx), barra,sprintf('Planar integration z'));
        if show_legend
        fprintf('Planar integration z %d%%\n', ...
             round(100*i/Nx));
    end

    for k = 1 : Ny

        t_ = sin(psi_).*cos(phi_).*x(i) + sin(psi_).*sin(phi_).*y(k) + cos(psi_).*z_;
        

                THETAz0 = phi_   +signo* acos(t_./(  complex(  R*abs( sin(psi_)  ))));
                
 
                THETAz0(imag(THETAz0) ~= 0) = NaN;
                THETAz0 = real(THETAz0);
                THETAz0 = mod(THETAz0, 2*pi);
                %THETAz0 = mod(THETAz0, phi(end));
                % o es lo mismo 
                % THETAz0 = mod(THETAz0, 2*pi- dphi); % no es dphi0
                
               
 

       

        

        % GUARDSA CON LAS VARIABLES


        alpha1 = 1*(  cos(psi_).*cos(alpha_).*cos(phi_) - sin(alpha_).*sin(phi_) ) ;
        alpha2 = 1*(  cos(psi_).*cos(alpha_).*sin(phi_) + sin(alpha_).*cos(phi_) ) ;
        alpha3 = 1*(- sin(psi_).*cos(alpha_) );


         polar = acos(-sin(psi_).*cos(alpha_));

         %polar = polar.*(polar<=pi/2) + (pi-polar).*(polar>pi/2);


         
         az =  atan2(alpha2,alpha1); 
      
         az =  mod( az, 2*pi);
         %az =  mod( az, phi(end));
         

    
         
         
         G0 = interpn(theta, psi, phi, Pz0, THETAz0, polar  , az   );

    

         
 
 
        

          G0((isnan(G0)))=0;

          G0 =G0.*sin(psi_);
          G0 = squeeze(sum(G0, 4))*dalpha;  
         
        
         F(i, k, :, :, : ) = reshape(G0, 1 ,1, size(G0,1), size(G0,2), size(G0,3));

    end

   
end

F = squeeze(sum(sum(F,5),4))*dpsi0*dphi0; 

 close(barra)


end