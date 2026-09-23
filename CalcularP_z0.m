function [Pf_z0] = CalcularP_z0(f, theta,  psi, phi, R, Nr, x, y, z, rMax)


[ theta_, psi_, phi_  ]= ndgrid( theta, psi, phi);

alpha1 = sin(psi_).*cos(phi_); 
alpha2 = sin(psi_).*sin(phi_);
alpha3 = cos(psi_);


dr = rMax/(Nr-1);
r  = 0 : dr : rMax ; %+ dr; % full array

 

barra = waitbar(0,'Calculation of the projection');


P = 0;


for m=1:Nr                
    m
 waitbar(gather(m/Nr));

 
    

                              
                                   

                         
                                    xi =   R*cos(theta_) + r(m)*alpha1 ;
                                    yi =   R*sin(theta_) + r(m)*alpha2 ; 
                                    zi =                   r(m)*alpha3; 
    
 


                                        
                                    a = interpn(x,y,z,f,xi,yi,zi);  % interpolation
                                    a( isnan(a))=0;                 % for points outside the limits(returned as NaN)                                    
                                    
                                    %sum(sum(sum(a)))

                                    %rm = (r(m+1) - r(m))/2 + r(m);
                                     P= P + r(m)*a;                          % acumulation
                                    % P= P + rm*a;                          % acumulation
                                    

end

  
Pf_z0=squeeze(P*dr);                                      % tu adjust dimension 

close(barra)

end


