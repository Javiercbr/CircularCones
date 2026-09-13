function [ ] =mostrar_muchos(mapa,vista,f,g)


Mg = max(max(max(g)));
Mf = max(max(max(f)));

mg = min(min(min(g)));
mf = min(min(min(f)));


M = max(Mf,Mg);
m = max(mf,mg);

clims = [m M];
 



% mapas: jet, gray, hot, bones

close all




    
    switch vista
    
    
    case '1'
        
        N = size(f,1);
 
    
    case '2'

        N = size(f,2);

        otherwise    

            N = size(f,3);
    end
    
    for i=1:1:N
    
    
    figure,
    
    
    switch vista
    
    
    case '1'
        
    subplot(121),imagesc(squeeze(f(i,:,:)), clims),colormap(mapa), colorbar%, axis equal
    subplot(122),imagesc(squeeze(g(i,:,:)), clims),colormap(mapa), colorbar%, axis equal
    
    
    case '2'
       
    
    subplot(121),imagesc(squeeze(f(:,i,:)),clims),colormap(mapa), colorbar, axis equal
    subplot(122),imagesc(squeeze(g(:,i,:)),clims),colormap(mapa), colorbar, axis equal   
        
        
    otherwise
        
    subplot(121),imagesc(squeeze(f(:,:,i)),clims),colormap(mapa), colorbar, axis equal
    subplot(122),imagesc(squeeze(g(:,:,i)),clims),colormap(mapa), colorbar, axis equal   
        
end

    
    
    end


end