function cargasGenericas
  clc
  clear
  %-------VARIABLES A MODIFICAR---------------------%
  ti=0;  #Tiempo inicial del que partiremos.Este seria el Tau
  tf=2;  #Ingrese hasta que tiempo vamos a medir desplazamiento (en segundos)
  dt=0.12;  #Ingrese el paso de muestreo de la fuerza
  tiempo=ti:dt:tf;
  Ftau=[0,1,4,9,9,6,0]; #Vector de fuerza tau que ira en yc e ys
  Ftau=cat(2,Ftau,zeros(1,length(tiempo)-length(Ftau)));
  #La frecuencia natural puede ser dato o no.
  m=0.2; #Masa
  k=8 ; #Ctte de ridigez
  wn=sqrt(k/m);
  c=0.4; #Ctte de amortiguamiento c. En caso de no haber amortiguamiento dejar en 0
  Ramort=c/(2*m*wn); #Si c no es dato simplemente colocar el valor de Ramort(Zitta)
  %----------------------------------------------------%

  #Debemos calcular el An y Bn
  Aviejo=0; #Inicializamos el valor del acumulador A (Como en la teoria, admitimos A0=0)
  Bviejo=0; #Inicializamos el valor del acumulador B (Como en la teoria, adminitos A0=0)
  t=0; #contador de tiempo
  if (ti==0)
    ycviejo=0;
    ysviejo=0;
  else
    Fvieja=input(['Ingrese el valor de la fuerza en el instante ',num2str(ti-dt),' segundos']); #Esto es solo si el tiempo inicial no es 0
    ycviejo=Fvieja*cos(wn*(ti-dt));
    ysviejo=Fvieja*sin(wn*(ti-dt));
  endif
  #Para resolver el problema utilizamos regla del trapecio
  comp=1;
  for (t=tiempo)
    ycnuevo=Ftau(comp)*cos(wn*t);
    ysnuevo=Ftau(comp)*sin(wn*t);
    Anuevo=Aviejo+(dt/(2*m*wn))*(ycviejo+ycnuevo);
    Bnuevo=Bviejo+(dt/(2*m*wn))*(ysviejo+ysnuevo);
    #Tambien debemos calcular las componentes del vector X(tn), el me entrega el desplazamiento en forma discreta
    X(comp)=Anuevo*sin(wn*t)-Bnuevo*cos(wn*t);
    #Ahora debemos intercambias las variables viejas por las nuevas
    Aviejo=Anuevo;
    ycviejo=ycnuevo;
    ysviejo=ysnuevo;
    comp++; #Incrementamos en 1 el contador de componentes
  endfor

  #Printeamos los valores de la carga generica, el desplazamiento y la fuerza elastica
  disp("La carga generica discretizada es:")
  disp(Ftau);
  disp(['El desplazamiento en el tiempo tf=',num2str(tf),' es:']);
  disp(X(end));
  disp(['La fuerza elastica en el tiempo tf=',num2str(tf),' es:']);
  disp(X(end)*k);
  #Finalmente graficamos la carga y la respuesta con respecto al tiempo y entregamos las respuestas
  figure(1)
  subplot(1,2,1);
  plot(tiempo,Ftau(1:length(tiempo)),'b'); #Con esa implementacion de 1:length(tiempo) logramos graficar la fuerza para t< al t final
  xlabel('Tiempo');
  ylabel('Carga');
  title('Carga generica');
  subplot(1,2,2);
  plot(tiempo,X,'r');
  xlabel('Tiempo');
  ylabel('Desplazamiento');
  title('Respuesta');
  grid on
endfunction
