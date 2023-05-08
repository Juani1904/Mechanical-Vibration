function cargaPeriodica
  clc
  clear
  close all

  %---------Parametros a modificar------------------%

  k = 40*pi^2; #Rigidez del sistema
  m = 1; #Masa del sistema
  Ramort = 0; #Relacion de amortiguamiento zitta
  Tp = 0.3; #Periodo de la carga
  tf=Tp; #Tiempo final hasta el que queremos printear la carga y respuesta
  numArmonicos = 8; %Numero de Armonicos
  dt=0.025;
  Fo=40;
  frecWf = (2*pi)/Tp; #Frecuencia forzada.Si nos dan Tp calcularla asi, si no colocar el valor que nos dan

  frecW = sqrt(k/m); #Frecuencia natural del sistema. La calculamos mediante k y m
  beta = frecWf/frecW; #Relacion de frecuencias
  N = Tp/dt; %Numero de puntos en los que discretizamos
  tiempo =0:dt:Tp; #Este es el paso del tiempo para calcular la fuerza hasta su periodo. La usamos para calculas los coeficientes
  tiempoGrafica=0:dt:tf; #Este es el paso del tiempo para graficar la fuerza y la respuesta mas veces para cualquier tiempo

  %-------------------COEFICIENTES DE FOURIER---------------------%

  #Ahora calculamos los coeficientes de fourier
  #Primero calculamos el coeficiente a0
  #contador=1;
    #for t=tiempo
      #vectorP(contador)=P(t);
      #contador++;
    #endfor
    a0=(4*Fo)/pi; #OJO,mirar bien el ejercicio con este valor porque confunde
  contador=1;
  for i=2:2:numArmonicos
    #for t=tiempo
      #Ahora seteamos la matriz A y B, cuyas filas corresponden al numero de armonicos y columnas la discretizacion
      #en el intervalo determinado por el vector t. Estas son las funciones que luego integramos con el metodo del trapecio
      #A(i,contador)=P(t)*cos(frecWf*i*t); #Esta es un producto de funciones en forma discreta
      #B(i,contador)=P(t)*sin(frecWf*i*t); #Esto es un producto de funciones en forma discreta
      #contador++;
    #endfor
    #Ahora calculamos los coeficientes am y bm para cada armonico
    am(contador)=-((4*Fo)/(pi*(i^2-1))); #Para cada modo esto es un escalar
    bm(contador)=0; #Para cada modo esto es un escalar
    contador++;
    #La funcion trapz calcula la integral numerica de A y B por el metodo del trapecio con un paso de tiempo dt
  endfor
  %-------------------------APROXIMACION DE LA CARGA-----------------------------------%
  #Una vez calculado los numArmonicos armonicos, podemos aplicar la serie trigonometrica de Fourier para aproximar nuestra funcion
  #Aca si usamos el tiempoGrafica, para printear la carga y respuesta hasta el tiempo que queramos
  contador=1;
  for t=tiempoGrafica
    sumaAuxiliarFourier=a0/2;
    contador2=1;
    for i=2:2:numArmonicos
      sumaAuxiliarFourier=sumaAuxiliarFourier+am(contador2)*cos(frecWf*i*t)+bm(contador2)*sin(frecWf*i*t);
      contador2++;
    endfor
    cargaFourier(contador)=sumaAuxiliarFourier;
    contador++;
  endfor

  %-------------------------RESPUESTA DEL SISTEMA-------------------------------%
  #Ahora vamos a calcular la respuesta del sistema
  contador=1;
  for t=tiempoGrafica
    sumaAuxiliarRespuestaFourier=a0/(2*k);
    contador2=1;
    for i=2:2:numArmonicos
      #Primero calculamos los parametros Dm y titam
      tanTita=(2*Ramort*beta*i)/(1-(beta*i)^2);
      tita=atan(tanTita);
      Dm=1/sqrt(((1-(beta*i)^2)^2)+(2*Ramort*beta*i)^2);
      #Ahora la serie de fourier
      sumaAuxiliarRespuestaFourier=sumaAuxiliarRespuestaFourier+(1/k)*(am(contador2)*Dm*cos(frecWf*i*t-tita)+bm(contador2)*Dm*sin(frecWf*i*t-tita));
      contador2++;
    endfor
    respuestaFourier(contador)=sumaAuxiliarRespuestaFourier;
    contador++;
  endfor

  #Printeamos el output de respuesta
  disp(['El valor de la carga en el instante tf= ',num2str(tf),'s es:'])
  disp(cargaFourier(end));
  disp(['El valor de la respuesta en el instante tf= ',num2str(tf),'s es:'])
  disp(respuestaFourier(end));
  disp(['El vector de carga aproximado hasta tf= ',num2str(tf),'s es:'])
  disp(cargaFourier);
  disp(['El vector de respuesta aproximado hasta tf= ',num2str(tf),'s es:'])
  disp(respuestaFourier);

  %--------------------------PLOTEO------------------------------------%
  figure(1)
  subplot(2,1,1);
  hold on
  grid on
  #plot(tiempo,vectorP,'b');
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
  Tp=3*pi;
  Po=10;
  if(t>=0) && (t<=2*pi)
    y=Po*sin((3*pi/Tp)*t);
  else
    y=0;
  endif
  %------------------------------------------------------%

endfunction

