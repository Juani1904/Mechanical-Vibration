function ParcialA3
    clc, clear

  po = 100000;
  w = zeros(1,100);
  w1 = [0,0.5,1,0.5,0]*po;
  w (1:length(w1)) = w1;
  nsteps = length(w);

  t1 = 0.4;
  dt = t1/4;
  t = 0:dt:(nsteps-1)*dt;

 m = 500000;
 k = 2*10^11 * 4*10^-5 / 1;
 zitta = 0.02;
ome = sqrt(k/m);
f = ome/(2*pi);
T = 1/f;
c = zitta*2*m*ome;
omed = ome*(1-zitta^2)^0.5;

y0 = 0;
yd0 = 0;

cw = cos(omed*t);
sw = sin(omed*t);

yc = w.*cw;
ys = w.*sw;

A = zeros(1,length(w));
B = zeros(1,length(w));

for i = 2:length(w)
  A(i) = A(i-1)*e^(-zitta*ome*dt)+ dt/(2*ome*m) * (yc(i-1)*e^(-zitta*ome*dt) +yc(i));
  B(i) = B(i-1)*e^(-zitta*ome*dt)+ dt/(2*ome*m) * (ys(i-1)*e^(-zitta*ome*dt) +ys(i));
endfor
X = A.*sin(omed*t) - B.*cos(omed*t);

P = interp1(t, w, t);

figure(1)
subplot(2,1,1)
plot(t, P, 'b')
subplot(2,1,2)
plot(t,X,'r')

% Encuentra el máximo solo en los valores positivos de X(1,:)
%indices_positivos = X(1,t <= 3) > 0;  % Crea una máscara lógica para los valores positivos
%valores_positivos = X(1, indices_positivos);  % Filtra los valores positivos de X(1,:)
% Encuentra el máximo y su índice
%[maximo, indice_maximo] = max(valores_positivos);
% Obtiene el tiempo correspondiente al máximo
%tiempo_maximo = t(indices_positivos);
%tiempo_maximo = tiempo_maximo(indice_maximo);
%maximo
%tiempo_maximo


  [maximo, indice_maximo] = max(X( t <= 3));
      tiempo_maximo = t(indice_maximo);

      maximo
      tiempo_maximo

%Para carga impulsiva, es decir, t1 menor a 1/4 T

if(t1<0.3*T)
imp = t1*po/2;
xaprox = imp/(omed*m) * sin(omed*(t - t1)) .* exp(-zitta*ome*(t - t1));

figure(2)
subplot(2,1,1)
plot(t, P, 'b')
subplot(2,1,2)
plot(t,xaprox,'r')

  [maximo, indice_maximo] = max(xaprox(t <= 3));
   tiempo_maximo = t(indice_maximo);

      maximo
      tiempo_maximo
endif

endfunction
