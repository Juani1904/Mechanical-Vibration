function cargasGenericas
  clc
  clear
  close all
  pkg load signal
  %-------VARIABLES A MODIFICAR---------------------%
  ti=0.01;  #Tiempo inicial del que partiremos.Este seria el Tau
  tf=3;  #Ingrese hasta que tiempo vamos a medir desplazamiento (en segundos)
  dt=0.01;  #Ingrese el paso de muestreo de la fuerza
  tiempo=ti:dt:tf;
  Ftau=Fcarga(tiempo); #Vector de fuerza tau que ira en yc e ys
  #Ftau=cat(2,Ftau,zeros(1,length(tiempo)-length(Ftau)));
  #La frecuencia natural puede ser dato o no.
  m=500000; #Masa
  E=2e11; #Modulo de young
  A=4e-5;
  L=1;
  k=(E*A)/L; #Ctte de ridigez
  wn=sqrt(k/m);
  #c=0; #Ctte de amortiguamiento c. En caso de no haber amortiguamiento dejar en 0
  Ramort=0.02; #Si c no es dato simplemente colocar el valor de Ramort(Zitta)

  %----------------------------------------------------%
  wd=sqrt(1-(Ramort)^2)*wn; #Frecuencia amortiguada, si el amortiguamiento es nulo sera igual a wn
  #Ahora debemos calcular el An y Bn
  #Primero definimos los vectores An y Bn
  An=zeros(1,length(Ftau));
  Bn=zeros(1,length(Ftau));
  #Tambien definimos los yc e ys
  yc=Ftau.*cos(wd.*tiempo);
  ys=Ftau.*sin(wd.*tiempo);
  #Ahora calculamos los coeficientes
  for i = 2:length(Ftau)
    An(i) = An(i-1)*exp(-Ramort*wn*dt)+ dt/(2*wd*m) * (yc(i-1)*e^(-Ramort*wn*dt) +yc(i));
    Bn(i) = Bn(i-1)*exp(-Ramort*wn*dt)+ dt/(2*wd*m) * (ys(i-1)*e^(-Ramort*wn*dt) +ys(i));
  endfor
  #Finalmente calculamos el vector de respuesta
  X = An.*sin(wd*tiempo) - Bn.*cos(wd*tiempo);

  #Printeamos los valores de la carga generica, el desplazamiento y la fuerza elastica
  #disp("La carga generica discretizada es:")
  #disp(Ftau);
  disp(['El desplazamiento en el tiempo tf=',num2str(tf),' es:']);
  disp(X(end));
  disp(['La fuerza elastica en el tiempo tf=',num2str(tf),' es:']);
  disp(X(end)*k);

  #Primero filtramos los valores negativos y nos quedamos solo con los positivos
  index1=1;
  index2=1;
  for t=tiempo
    if X(index2) > 0
      respuestaXPositiva(index1)=X(index2);
      tiempoPositivo(index1)=t;
      index1++;
    endif
    index2++;
  endfor
  [picos,index] = findpeaks(respuestaXPositiva);
  disp(['Los 2 valores maximos en la fase POSITVA de la respuesta entre ',num2str(ti),' y ',num2str(tf),' son:']);
  disp("Pico 1");
  disp(['Primer pico: ',num2str(picos(1)),' en t=',num2str(tiempoPositivo(index(1)))]);
  disp(['Segundo pico: ',num2str(picos(2)),' en t=',num2str(tiempoPositivo(index(2)))]);

  #Finalmente graficamos la carga y la respuesta con respecto al tiempo y entregamos las respuestas
  figure(1)
  subplot(2,1,1);
  plot(tiempo,Ftau(1:length(tiempo)),'b'); #Con esa implementacion de 1:length(tiempo) logramos graficar la fuerza para t< al t final
  xlabel('Tiempo');
  ylabel('Carga');
  title('Carga generica');
  subplot(2,1,2)
  plot(tiempo,X,'r');
  xlabel('Tiempo');
  ylabel('Desplazamiento');
  title('Respuesta');
  grid on
endfunction


function y=Fcarga(tiempo)

#En este caso no nos dan la carga discretizada por lo que tenemos que hacer una discretizacion nosotros
Fo=100000;
to=0.4;
i=1;
for t=tiempo
  if t<=to/2
    y(i)=((Fo-0)/((to/2)-0))*(t-0)+0;
  elseif t>to/2 && t<=to
    y(i)=((0-Fo)/(to-(to/2)))*(t-(to/2))+Fo;
  else
    y(i)=0;
  endif
  i++;
endfor
endfunction
