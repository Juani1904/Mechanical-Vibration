function fourierJuancho
clc, clear
    k = 1000;
    m = 100;
    zitta = 0.1;
    Tp = 3*pi;
    wn = sqrt(k/m);
    w = 3*pi/Tp;
    beta = w/wn;
    n = 10; %Armonicos
    N = 1000; %Numero de puntos
    dt=(Tp-0)/(N-1); %Paso del tiempo
    t = 0:dt:Tp;
    tg=0:dt:10;   %Tiempo Grafico

    A=zeros(length(t),n);
    B=zeros(length(t),n);
    an=zeros(n,1);
    bn=zeros(n,1);
    a0=2/Tp*trapz(dt,P(t));

% Calcula los coeficientes de la serie de Fourier
    for j=1:n
       A(:,j)= P(t).*cos(w*j*t);
       B(:,j)= P(t).*sin(w*j*t);
       an(j)=2/Tp*trapz(dt,A(:,j));
       bn(j)=2/Tp*trapz(dt,B(:,j));
     end
     disp(an)

    % Suma los términos de la serie de Fourier
    fourier=zeros(length(tg),n+1);
    for(i=1:length(tg))
        fourier(i,1)=a0/2;
        for j=1:(n)
        tgtitam = (2*zitta*j*w/wn)/ (1-(j*w/wn)^2);
        titam = atan(tgtitam);
        fourier(i,j+1)=bn(j)*sin(tg(i)*w*j - titam)+an(j)*cos(tg(i)*w*j - titam);
        end
    end

    % Suma los términos de la serie de Fourier
    suma_fourier = a0/2*ones(size(tg));
    for j = 1:n
        suma_fourier = suma_fourier + bn(j)*sin(w*j*tg) + an(j)*cos(w*j*tg);
    end

    %Calculo de la respuesta del sistema
    x=zeros(length(tg),1);
    for i=1:length(tg)
        x(i)=fourier(i,1)/k;
        for j=1:(n)
        x(i)=x(i)+fourier(i,j+1)*(1/sqrt((1-(j*w/wn)^2)^2+(2*zitta*j*w/wn)^2))*1/k;
        end
    end

    % Graficar la señal original y la aproximación de Fourier
    figure (1)
    plot(tg,P(tg),'b');
    hold on; % Indica que se van a graficar varias señales en la misma figura
    plot(tg,suma_fourier,'r');
    title('Señal original y aproximación de Fourier');
    xlabel('Tiempo (s)');
    ylabel('Amplitud');

    figure (2)
    % Graficar la respuesta del sistema x(t)
    plot(tg,x);
    title('Respuesta del sistema x(t)');
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
end

function y = P(t)
        Tp = 3*pi;
        t_mod = mod(t, Tp);
        y = 10*sin(3*pi*t_mod/Tp).*(t_mod>=0 & t_mod<=2*pi) + 0.*(t_mod<0 | t_mod>2*pi);
        if 0
      % Aquí colocas el bloque de código que deseas comentar

          Para el primer tramo, desde t=0 hasta t=Tp/2, se puede usar la ecuación de una recta con pendiente m1 que pasa por el punto (0,0) de la forma:
         y = m1 * t
          Para encontrar la pendiente m1, se puede notar que la señal triangular aumenta desde y=0 hasta y=pi en el primer tramo. Esto corresponde a una pendiente positiva de pi/(Tp/2) = 2*pi/Tp. Por lo tanto, la ecuación de la recta para este tramo sería:
        y = (2*pi/Tp) * t

        Para el segundo tramo, desde t=Tp/2 hasta t=Tp, se puede usar la ecuación de una recta con pendiente m2 que pasa por el punto (Tp/2, pi) de la forma:
        y = m2 * (t - Tp/2) + pi
        Para encontrar la pendiente m2, se puede notar que la señal triangular disminuye desde y=pi hasta y=0 en el segundo tramo. Esto corresponde a una pendiente negativa de pi/(Tp/2) = -2*pi/Tp. Por lo tanto, la ecuación de la recta para este tramo sería:
        y = (-2*pi/Tp) * (t - Tp/2) + pi
        Juntando estas dos ecuaciones de rectas, se obtiene la expresión final para la señal triangular con las mismas pendientes que la original:
   endif
      end
