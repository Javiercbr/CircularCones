function [P] = CalcularP_y0_radon_clasica2(theta, varphi, s, gf_y0, psi, phi, show_legend)

 


dvarphi = varphi(2)- varphi(1);
dvarphi0 = 8*dvarphi;
varphi0  = -pi/2:  dvarphi0  :  pi/2;
 

barra = waitbar(0, 'Processing...'); 

Ntheta = length(theta);

P = zeros(Ntheta, length(psi), length(phi));
 
for i = 1 : Ntheta

    waitbar(gather(i/Ntheta), barra,sprintf('Processing Pfy theta %d/%d', i, Ntheta));
        
    
    if show_legend
        fprintf('Processing Pfy %d%%\n', ...
             round(100*i/Ntheta));
    end

[  psi_, phi_, varphi_  ]= ndgrid(   psi, phi, varphi0 );




t_ = -(sin(psi_).*cos(phi_).*cos(theta(i) - varphi_) + cos(psi_).*sin(theta(i) - varphi_ ));

 
Py0_prima = interpn(  varphi, s, squeeze( gf_y0(i,:,:)),  varphi_, t_);
Py0_prima((isnan(Py0_prima)))=0;
Py0_prima = sum(Py0_prima.*abs(sin(psi_).*sin(phi_)),3);
Py0_prima = squeeze(Py0_prima*dvarphi0);



P(i,:,:)= Py0_prima;

end

close(barra);


end
