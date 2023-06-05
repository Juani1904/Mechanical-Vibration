function cargaFourier=Pperiodica(tf,dt)
  clc
  clear
  close all

  %---------Parametros a modificar------------------%

  k = 1000; #Rigidez del sistema
  m = 100; #Masa del sistema
  Ramort = 0.1; #Relacion de amortiguamiento zitta
  Fo=1;
  m=1/(2*pi)^2;
  k=1;
  wn=sqrt(k/m);
  wc=(3/2)*wn;
  Tp=(2*pi)/wc
  #Tp = 3*pi; #Periodo de la carga
  tf=10; #Tiempo final hasta el que queremos printear la carga y respuesta
  numArmonicos = 10; %Numero de Armonicos
  N = 1000; %Numero de puntos en los que discretizamos
  frecWf = 2*pi/Tp; #Frecuencia forzada.Si nos dan Tp calcularla asi, si no colocar el valor que nos dan

  frecW = sqrt(k/m); #Frecuencia natural del sistema. La calculamos mediante k y m
  beta = frecWf/frecW; #Relacion de frecuencias
  dt=Tp/N; %Paso del tiempo. Viene de Tp=N*dt
  tiempo =0:dt:Tp; #Este es el paso del tiempo para calcular la fuerza hasta su periodo. La usamos para calculas los coeficientes
  tiempoGrafica=0:dt:tf; #Este es el paso del tiempo para graficar la fuerza y la respuesta mas veces para cualquier tiempo

  %-------------------COEFICIENTES DE FOURIER---------------------%

  #Ahora calculamos los coeficientes de fourier
  #Primero calculamos el coeficiente a0
  contador=1;
    for t=tiempo
      vectorP(contador)=P(t);
      contador++;
    endfor
    a0=(2/Tp)*trapz(dt,vectorP);

  for i=1:numArmonicos
    contador=1;
    for t=tiempo
      #Ahora seteamos la matriz A y B, cuyas filas corresponden al numero de armonicos y columnas la discretizacion
      #en el intervalo determinado por el vector t. Estas son las funciones que luego integramos con el metodo del trapecio
      A(i,contador)=P(t)*cos(frecWf*i*t); #Esta es un producto de funciones en forma discreta
      B(i,contador)=P(t)*sin(frecWf*i*t); #Esto es un producto de funciones en forma discreta
      contador++;
    endfor
    #Ahora calculamos los coeficientes am y bm para cada armonico
    am(i)=(2/Tp)*trapz(dt,A(i,:)); #Para cada modo esto es un escalar
    bm(i)=(2/Tp)*trapz(dt,B(i,:)); #Para cada modo esto es un escalar
    #La funcion trapz calcula la integral numerica de A y B por el metodo del trapecio con un paso de tiempo dt
  endfor

  %-------------------------APROXIMACION DE LA CARGA-----------------------------------%
  #Una vez calculado los numArmonicos armonicos, podemos aplicar la serie trigonometrica de Fourier para aproximar nuestra funcion
  #Aca si usamos el tiempoGrafica, para printear la carga y respuesta hasta el tiempo que queramos
  contador=1;
  for t=tiempoGrafica
    sumaAuxiliarFourier=a0/2;
    for i=1:numArmonicos
      sumaAuxiliarFourier=sumaAuxiliarFourier+am(i)*cos(frecWf*i*t)+bm(i)*sin(frecWf*i*t);
    endfor
    cargaFourier(contador)=sumaAuxiliarFourier;
    contador++;
  endfor

  %-------------------------RESPUESTA DEL SISTEMA-------------------------------%
  #Ahora vamos a calcular la respuesta del sistema
  contador=1;
  for t=tiempoGrafica
    sumaAuxiliarRespuestaFourier=a0/(2*k);
    for i=1:numArmonicos
      #Primero calculamos los parametros Dm y titam
      tanTita=(2*Ramort*beta*i)/(1-(beta*i)^2);
      tita=atan(tanTita);
      Dm=1/sqrt(((1-(beta*i)^2)^2)+(2*Ramort*beta*i)^2);
      #Ahora la serie de fourier
      sumaAuxiliarRespuestaFourier=sumaAuxiliarRespuestaFourier+(1/k)*(am(i)*Dm*cos(frecWf*i*t-tita)+bm(i)*Dm*sin(frecWf*i*t-tita));
    endfor
    respuestaFourier(contador)=sumaAuxiliarRespuestaFourier;
    contador++;
  endfor

  #Printeamos el output de respuesta
  disp(['El valor de la carga en el instante tf= ',num2str(tf),'s es:'])
  disp(cargaFourier(end));
  disp(['El valor de la respuesta en el instante tf= ',num2str(tf),'s es:'])
  disp(respuestaFourier(end));

  %--------------------------PLOTEO------------------------------------%
  figure(1)
  subplot(2,1,1);
  hold on
  grid on
  plot(tiempo,vectorP,'b');
  plot(tiempoGrafica,cargaFourier,'r');
  title('Carga real y aprox con Fourier');
  xlabel('Tiempo');
  ylabel('Desplazamiento');
  hold off
  subplot(2,1,2);
  grid on
  plot(tiempoGrafica,respuestaFourier,'r');
  title('Respuesta aprox con Fourier');
  xlabel('Tiempo');
  ylabel('Desplazamiento');
endfunction


function y = P(t)
  %------------Parametros a modificar--------------%
  Fo=1;
  m=1/(2*pi)^2;
  k=1;
  wn=sqrt(k/m);
  wc=(3/2)*wn;
  T=(2*pi)/wc
  if(t>=0) && (t<=T)
    y=(Fo/T)*t;
  else
    y=0;
  endif
  %------------------------------------------------------%

endfunction

