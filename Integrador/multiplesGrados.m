function multiplesGrados
  clc
  clear all
  close all

  %----------------------------SISTEMA----------------------------------------------

  #Lo primero que hacemos es tomar los datos que nos proporciona el ejercicio
  ms=38500;
  J=2.446e6;
  m1=4330;
  m2=m1;
  a1=4.2;
  a2=a1;
  ks1=2.535e6;
  ks2=ks1;
  cs1=1.97e5;
  cs2=cs1;
  kt1=4.28e6;
  kt2=kt1;
  ct1=9.8e4;
  ct2=ct1;
  Ramort=[0.1,0.1,0.02,0.02]; #Asignar los valores de zitta (relacion de amortiguamiento) si los tuviera
  #En este caso se trata de un oscilador que consta de varias masas
  #Lo modelamos con un carrito con 3 resortes de distinta rigidez k y la misma masas
  #Primero debemos hallar la ecuacion de movimiento en papel mediante algun metodo


  #Armamos las matrices M y K
  M=[J 0 0 0;
     0 ms 0 0;
     0 0 m1 0;
     0 0 0 m2];

  K=[((a1^2)*ks1+(a2^2)*ks2) a1*ks1-a2*ks2 -a1*ks1 a2*ks2;
    a1*ks1-a2*ks2 ks1+ks2 -ks1 -ks2;
   -a1*ks1 -ks1 ks1+kt1 0;
  a2*ks2 -ks2 0 ks2+kt2];
  #Asignamos una variable que sera el orden de la matriz cuadrada
  numGL=size(M)(1);

  #Una vez teniendo las matrices M y K debemos resolver el problema de vectores propios
  [Xnorm,lambda]=eig(K,M); #Esto me entrega el autovector normalizado con respecto a la matriz de masa modal y los autovalores

  #Hallamos los autovalores, que son la raiz cuadrada de los valores propios. Armamos un vector columna tambien
  frecW=diag(lambda.^0.5);
  #Printeamos los autovectores y frecuencias
  disp("Los autovalores son: ");
  disp(diag(lambda));
  disp("Los modos o frecuencias naturales de cada nodo son: ");
  disp(frecW);
  disp("Los autovectores normalizados segun la matriz de masa son: ");
  disp("\n")
  for i=1:numGL
    disp(Xnorm(:,i));
    disp("\n")
  endfor
  #Ahora normalizamos nuevamente los autovectores, pero esta vez respecto a su primer elemento
  for i = 1:length(Xnorm(:,1))
    XnormPrima(:,i) = Xnorm(:,i)/Xnorm(1,i);
  endfor
  #Printeamos los autovectores normalizados segun el primer elemento
  disp("Los autovectores normalizados segun su primer elemento son: ");
  disp("\n")
  for i=1:numGL
    disp(XnormPrima(:,i));
    disp("\n")
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
    #Segun el tipo de carga que nos encontremos en el vector de cargas el metodo que utilizaremos para resolver
    tipoDeCarga=2;
    #Si el vector de cargas es nulo, el movimiento del sistema en todos sus GL sera libre. Colocar 1
    #Si el vector de carga tiene una o mas cargas armonicas, o combinacion entre 0 y cargas armonicas Coloca 2
    #Si el vector de carga tiene una o mas cargas periodicas, o combinacion entre cargas armonicas, periodicas y 0 Coloca 3
    #Si el vector de carga tiene cargas genericas o mezcla entre genericas, armonicas, periodicas y 0. Coloca 4

    #PARAMETROS GENERALES A COMPLETAR (SE DEBE MODIFICARRR)

    tf=2; #Tiempo final (introducir el valor al que se quiere obtener)
    dt=0.001; #Paso del tiempo
    #Entonces vamos a tener un paso del tiempo comun a todos los demas de:
    tiempo=0:dt:tf; #Esta variable no se modifica

    #Solamente se utiliza en sistema libre, cuando no sea ese el caso simplemente colocar 0 en ese GL
    x0=[0;0;0;0];
    x0prima=[0;0;0;0];

    #Vector de frecuencias forzadas y relacion beta. Colocar la fuerza en el modo que corresponda, en otro caso dejar en 0)
    longOnda=0.1; #Longitud de onda en metros
    velMovil=40*1000/3600; #Velocidad movil en metros por segundo

    frecWf=[0,0,2*pi*velMovil/longOnda,2*pi*velMovil/longOnda];
    beta=frecWf./frecW;
    contador=1;
    for t=tiempo
      F(:,contador)=[0,0,(t<= 1.501)*0.03*sin(frecWf(3)*t),(t>= 0.729 && 2.2292>t)*0.03*sin(frecWf(4)*t)]; #Esta solo sirve para printear la carga
      phi3=atan((2*Ramort(3)*beta(3))/(1-beta(3)^2));
      phi4=atan((2*Ramort(4)*beta(4))/(1-beta(4)^2));
      Fauxresp(:,contador)=[0,0,(t<= 1.501)*0.03*sin(frecWf(3)*t-phi3),(t>= 0.729 && 2.2292>t)*0.03*sin(frecWf(4)*t-phi4)]; #Esta es la que pasamos porque tiene el desfase
      contador++;
    endfor

    #LIBRE (SOLO CUANDO NO HAY CARGA EN NINGUN GL)
    if tipoDeCarga==1
      #Pasamos a coordenadas modales
      y0=Xnorm'*M*x0;
      y0prima=Xnorm'*M*x0prima;

      for i=1:numGL
        contador=1;
        for t=tiempo
          y(i,contador)=(y0(i)*cos(frecWd(i)*t)+((y0prima(i)+Ramort(i)*frecW(i)*y0(i))/frecWd(i))*sin(frecWd(i)*t))*exp(-Ramort(i)*frecW(i)*t);
          contador++;
        endfor
      endfor
    endif

    #CARGA ARMONICA
    if tipoDeCarga==2

      #Pasamos la carga a coordenadas modales
      Fmodal=Xnorm'*Fauxresp;
      for i=1:numGL
        contador=1;
        for t=tiempo
          #Calculamos los parametros para posteriormente calcular y
          D=1/(sqrt((1-beta(i)^2)^2+(2*Ramort(i)*beta(i))^2));
          #Calculamos la respuesta en coordenadas modales y con la formula correspondiente a respuesta a cargas armonicas
          y(i,contador)=(Fmodal(i,contador)/Km(i,i))*D; #Esta formula se hace teniendo en cuenta que la division de una matriz por otra en realidad es el producto de una por la inversa de la otra
          contador++;
        endfor
        #Adicionalmente creamos el vector de fases de cada modo
        phi(i)=atan((2*Ramort(i)*beta(i))/(1-beta(i)^2));
      endfor

    endif

    #CARGA PERIODICA

    #CARGA GENERICA

  #Ahora pasamos la respuesta y a coordenadas geometricas nuevamente
  respuestaX=Xnorm*y;
  #Printeamos las cargas a las que se ve sometido el sistema
  figure(2);
  for i=1:numGL
    subplot(numGL,1,i)
    plot(tiempo,F(i,:));
    title(['Carga ',num2str(i)])
    xlabel("Tiempo")
    ylabel("Desplazamiento")
  endfor
  #Printeamos la respuesta del sistema (los 3 G.L), utilizando como eje de abscisas el tiempo y como ordenadas la respuesta en coord geometricas
  figure(3);
  for i=1:numGL
    subplot(numGL,1,i)
    plot(tiempo,y(i,:));
    title(['Grado de libertad ',num2str(i)])
    xlabel("Tiempo")
    ylabel("Desplazamiento")
  endfor
  hold off
  #Entregamos el valor del ultimo punto (en tf)
  disp(['El valor de desplazamiento de todos los grados de libertad en el tiempo t= ',num2str(tf),' es: '])
  disp(respuestaX(:,end));
  if tipoDeCarga==2
    disp("La fase de cada modo en radianes es: ");
    disp(phi');
    disp("La fase de cada modo en grados es: ");
    disp(rad2deg(phi'));
  endif

  endfunction
