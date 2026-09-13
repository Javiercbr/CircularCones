function [P] = CalcularP_z0_radon_clasica2(theta, varphi, s, gf_z0, psi, phi)


dvarphi = varphi(2)- varphi(1);
dvarphi0 = 8*dvarphi;
varphi0  = -pi/2:  dvarphi0  :  pi/2;
 
Ntheta = length(theta);

P = zeros(Ntheta, length(psi), length(phi));
 
for i = 1 : Ntheta

    waitbar(gather(i/Ntheta) );
[  psi_, phi_, varphi_  ]= ndgrid(   psi, phi, varphi0 );

 
t_ = -sin(psi_).*cos(phi_ - varphi_ - theta(i));
 

Pz0_prima = interpn(  varphi, s, squeeze(gf_z0(i,:,:)),  varphi_, t_);
Pz0_prima((isnan(Pz0_prima)))=0;

%M = 0.02;
%a = M*(abs(cos(psi_))<M) +  (abs(cos(psi_))>=M).*(abs(cos(psi_)));

%Pz0_prima = sum(Pz0_prima.*a,3);



Pz0_prima = sum(Pz0_prima.*abs(cos(psi_)),3);

Pz0_prima = squeeze(Pz0_prima*dvarphi0);

P(i,:,:)= Pz0_prima;

end



end