function ejercicio3
  clc
  a=[0.1 0.25 1];
  b=1;
  c=1;

  for i=1:1:length(a)
    #Escribimos el determinante
    determinante=b^2-4*a(i)*b;
    if determinante > 0
      disp(['La ecuacion ',num2str(a(i)),'x^2+',num2str(b),'x+',num2str(c),' tiene raices reales distintas']);
    elseif determinante==0
      disp(['La ecuacion ',num2str(a(i)),'x^2+',num2str(b),'x+',num2str(c),' tiene raices reales iguales']);
    elseif determinante < 0
      disp(['La ecuacion ',num2str(a(i)),'x^2+',num2str(b),'x+',num2str(c),' tiene raices complejas conjugadas']);
    endif
  endfor


endfunction
