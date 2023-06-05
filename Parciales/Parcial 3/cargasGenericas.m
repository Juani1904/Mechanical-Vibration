function cargasGenericas
  clc
  clear
  close all
  pkg load signal
  %-------VARIABLES A MODIFICAR---------------------%
  ti=0.0;  #Tiempo inicial del que partiremos.Este seria el Tau
  tf=10;  #Ingrese hasta que tiempo vamos a medir desplazamiento (en segundos)
  dt=0.0556;  #Ingrese el paso de muestreo de la fuerza
  tiempo=ti:dt:tf;
  Fperiodica=Pperiodica();
  #Ftau=cat(2,Ftau,zeros(1,length(tiempo)-length(Ftau)));
  #La frecuencia natural puede ser dato o no.
  m=1/(2*pi)^2; #Masa
  k=1; #Ctte de ridigez
  wn=sqrt(k/m);
  wc=(3/2)*wn;
  Tp=(2*pi)/wc;
  N=1000;
  dtFou=0.6667/N; #0.667 es el Tp
  contador=1;
  for t=tiempo
    multi=1;
    while(t> multi*dtFou)
      multi++;
    endwhile
    Ftau(contador)=Fperiodica(multi); #Vector de fuerza tau que ira en yc e ys
    contador++;
  endfor

  #c=0; #Ctte de amortiguamiento c. En caso de no haber amortiguamiento dejar en 0
  Ramort=0.1; #Si c no es dato simplemente colocar el valor de Ramort(Zitta)

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
  hold on
  #subplot(2,1,1);
  plot(tiempo,Ftau(1:length(tiempo)),'b'); #Con esa implementacion de 1:length(tiempo) logramos graficar la fuerza para t< al t final
  xlabel('Tiempo');
  ylabel('Carga');
  title('Carga generica');
  #subplot(2,1,2)
  plot(tiempo,X,'r');
  xlabel('Tiempo');
  ylabel('Desplazamiento');
  title('Respuesta');
  hold off
  grid on
endfunction


function cargaFourier=Pperiodica()
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
  tf=10;
  #Tp = 3*pi; #Periodo de la carga
  #tf=20; #Tiempo final hasta el que queremos printear la carga y respuesta
  numArmonicos = 15; %Numero de Armonicos
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
  #figure(1)
  #subplot(2,1,1);
  #hold on
  #grid on
  #plot(tiempo,vectorP,'b');
  #plot(tiempoGrafica,cargaFourier,'r');
  #title('Carga real y aprox con Fourier');
  #xlabel('Tiempo');
 # ylabel('Desplazamiento');
  #hold off
  #subplot(2,1,2);
 # grid on
  #plot(tiempoGrafica,respuestaFourier,'r');
  #title('Respuesta aprox con Fourier');
  #xlabel('Tiempo');
  #ylabel('Desplazamiento');
endfunction


function y = P(t)
  %------------Parametros a modificar--------------%
  Fo=1;
  m=1/(2*pi)^2;
  k=1;
  wn=sqrt(k/m);
  wc=(3/2)*wn;
  T=(2*pi)/wc;
  if(t>=0) && (t<=T)
    y=(Fo/T)*t;
  else
    y=0;
  endif
  %------------------------------------------------------%

endfunction

