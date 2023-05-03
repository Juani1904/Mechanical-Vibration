function ejercicio2
  clc
  clear all
  contador=1;
  vectorIndex=1;
  vector=[];
  do
    if (  mod(contador,2)==0 && mod(contador,7)==0 && mod(contador,13)==0 )
      vector(vectorIndex)=contador;
      vectorIndex++;
      disp(contador);
    endif
    contador++;
  until (contador==2000)

  disp("El vector es: "), disp(vector);

  endfunction
