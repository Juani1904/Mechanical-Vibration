function ParcialA
    clc, clear
    k = 40*pi^2;
    P0 = 40;
    m = 1;
    Tp = 0.3;
    w = 2*pi/Tp;
    wn = sqrt(k/m);
    beta = w/wn;
    n = 4;
    zitta = 0;

    dt = 0.025;
    t = 0:dt:Tp;



    A=zeros(length(t),n);
    B=zeros(length(t),n);
    an=zeros(n,1);
    bn=zeros(n,1);
    a0=4*P0/pi;

    % Calcula los coeficientes de la serie de Fourier
    for j=1:n
      if mod(j,2) == 0
        an(j)= -4*P0/pi * sum(1/(j^2-1) * cos(w*j*t));
       endif
     end


    fourier=zeros(length(t),n+1);
    for(i=1:length(t))
        fourier(i,1)=a0/2;
        for j=1:(n)
             tgtitam = (2*zitta*j*w/wn)/ (1-(j*w/wn)^2);
            titam = atan(tgtitam);
            fourier(i,j+1)=bn(j)*sin(t(i)*w*j - titam)+an(j)*cos(t(i)*w*j - titam);
        end
    end

     suma_fourier = a0/2*ones(size(t));
    for j = 1:n
          suma_fourier = suma_fourier + bn(j)*sin(w*j*t) + an(j)*cos(w*j*t);
    end

     %Calculo de la respuesta del sistema
    x=zeros(length(t),1);
    for i=1:length(t)
        x(i)=fourier(i,1)/k;
        for j=1:(n)
        x(i)=x(i)+fourier(i,j+1)*(1/sqrt((1-(j*w/wn)^2)^2+(2*zitta*j*w/wn)^2))*1/k;
        end
    end

      x

    figure(1)
    plot(t, suma_fourier, 'b');
    figure(2)
    plot(t,x,'r')
endfunction

