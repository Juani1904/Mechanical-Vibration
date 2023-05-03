function ejercicio5
  clc
  t=[0:2*pi/50:10*pi];
  x=2*cos(t);
  y=5*sin(t);

  #Primero ploteamos x e y en funcion de t
  subplot(1,2,1)
  plot(t,x)
  title("Grafica x-t")
  subplot(1,2,2)
  plot(t,y)
  title("Grafica y-t")

  #Si quisieramos graficar una elipse en forma parametrica lo hacemos con lo siguiente
  plot(x,y)
  axis equal
  #Graficamos la curva que nos piden en forma parametrica, para ello hacemos uso de una subrutina
  index=1;
  for i=t
    coordenadasX(index)=funcionX(i);
    coordenadasY(index)=funcionY(i);
    index++;
  endfor
  #Ploteamos ambas funciones
  subplot(2,1,1)
  plot(t,coordenadasX)
  subplot(2,1,2)
  plot(t,coordenadasY)

  #Pegamos el codigo que nos dejaron en el tp para realizar la animación
  figure, axis([-3 3 -4 4]), hold on
  for i=1:length(t)-1
  plot(x(i:i+1),y(i:i+1));
  pause(1/10);
  end
  hold off
endfunction

function[valorSalida]=funcionX(t)
  valorSalida=cos(t)*(exp(cos(t))-2*cos(4*t)-(sin(t/12))^2);
end

function[valorSalida]=funcionY(t)
  valorSalida=sin(t)*(exp(cos(t))-2*cos(4*t)-(sin(t/12))^2);
end
