function ejercicio1
  clc
  clear all
  close all

  %----------------------------SISTEMA----------------------------------------------

  #Lo primero que hacemos es tomar los datos que nos proporciona el ejercicio
  m=2;
  k1=600;
  k2=1200;
  k3=2400;
  Ramort=[0,0,0]; #Asignar los valores de zitta (relacion de amortiguamiento) si los tuviera
  #En este caso se trata de un oscilador que consta de varias masas
  #Lo modelamos con un carrito con 3 resortes de distinta rigidez k y la misma masas
  #Primero debemos hallar la ecuacion de movimiento en papel mediante algun metodo


  #Armamos las matrices M y K
  M=[m 0 0;
     0 m 0;
     0 0 m];

  K=[(k2+k3) -k2 0;
     -k2 (k1+k2) -k1
      0 -k1 k1];
  #Asignamos una variable que sera el orden de la matriz cuadrada
  numGL=size(M)(1);

  #Una vez teniendo las matrices M y K debemos resolver el problema de vectores propios
  [Xnorm,lambda]=eig(K,M); #Esto me entrega el autovector normalizado con respecto a la matriz de masa modal y los autovalores

  #Hallamos los autovalores, que son la raiz cuadrada de los valores propios. Armamos un vector columna tambien
  frecW=diag(lambda.^0.5);
  #Ahora normalizamos nuevamente los autovectores, pero esta vez respecto a su primer elemento
  for i = 1:length(Xnorm(:,1))
    XnormPrima(:,i) = Xnorm(:,i)/Xnorm(1,i);
  endfor

  #Creamos un vector de valores y que represente la altura de la torre y el paso por las masas
  ejeY=0:numGL;
  hold on
  figure(1);
  for i=1:numGL
    subplot(1,numGL,i);
    plot([0,Xnorm(:,i)'],ejeY);
    title(['Modo',num2str(i)]);
    xlabel('Desplazamiento');
    ylabel('Altura');
  endfor

  #Ahora debemos Calcular la matriz de masa modal,la matriz de ridigez modal y la matriz de amortiguamiento modal
  #MATRIZ MASA MODAL
  Mm=Xnorm'*M*Xnorm; #Siempre y cuando usemos como matriz de vectores respuesta, la matriz de vectores normalizados respecto a la masa modal
                     #(la que entrega la funcion eig) esto nos va a dar la matriz identidad

  #MATRIZ RIGIDEZ MODAL
  Km=Xnorm'*K*Xnorm; #Lo mismo que arriba, pero va a dar una matriz diagonal de frecuencias modales al cuadrado

  #MATRIZ MODAL DE AMORTIGUAMIENTO
  for i=1:numGL
    auxCm(i)=2*Ramort(i)*frecW(i);
  endfor
  Cm=diag(auxCm);
  frecWd=frecW.*(1-Ramort'.^2).^0.5; #La frecuencia de amortiguamiento

  %-----------------------------------CARGAS--------------------------------------
  for i=1:numGL
    disp(['Ingrese tipo de carga para GL',num2str(i),'1 Libre / 2 Armonica / 3 Periodica / 4 Generica ']);
    tipoDeCarga=input("Ingrese 1, 2, 3 o 4: ");

    #PARAMETROS GENERALES A COMPLETAR (SE DEBE MODIFICARRR)

    tf=2*pi/frecW(1); #Tiempo final (introducir el valor al que se quiere obtener)
    dt=0.01; #Paso del tiempo
    #Entonces vamos a tener un paso del tiempo comun a todos los demas de:
    tiempo=0:dt:tf; #Esta variable no se modifica

    #Setear valores iniciales con tantos elementos como GL tenga el sistema y en su orden correspondiente
    #Solamente se utiliza en sistema libre, cuando no sea ese el caso simplemente colocar 0 en ese GL
    x0=[0.3;-0.8;0.3];
    x0prima=[0;0;0];


    #LIBRE
    if tipoDeCarga==1
      #Pasamos a coordenadas modales
      y0=Xnorm'*x0;
      y0prima=Xnorm'*x0prima;
      contador=1;
      for t=tiempo
        y(i,contador)=(y0(i)*cos(frecWd(i)*t)+((y0prima(i)+Ramort(i)*frecW(i)*y0(i))/frecWd(i))*sin(frecWd(i)*t))*exp(-Ramort(i)*frecW(i)*t);
        contador++;
      endfor
    endif

    #CARGA ARMONICA

    #CARGA PERIODICA

    #CARGA GENERICA
  endfor

  #Ahora pasamos la respuesta y a coordenadas geometricas nuevamente
  respuestaX=Xnorm*y;
  #Printeamos las cargas a las que se ve sometido el sistema
  figure(2);
  for i=1:numGL
    subplot(numGL,1,i)
    plot(tiempo,diag(zeros(length(respuestaX))));
    title(['Carga ',num2str(i)])
    xlabel("Tiempo")
    ylabel("Desplazamiento")
  endfor
  #Printeamos la respuesta del sistema (los 3 G.L), utilizando como eje de abscisas el tiempo y como ordenadas la respuesta en coord geometricas
  figure(3);
  for i=1:3
    subplot(3,1,i)
    plot(tiempo,respuestaX(i,:));
    title(['Grado de libertad ',num2str(i)])
    xlabel("Tiempo")
    ylabel("Desplazamiento")
  endfor
  hold off
  #Entregamos el valor del ultimo punto (en tf)
  disp(['El valor de todos los grados de libertad en el tiempo t= ',num2str(tf),' es: '])
  disp(respuestaX(:,end));

  endfunction
