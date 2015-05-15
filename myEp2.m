1;

#=======================================================================
# As seguintes suposicoes foram feitas para a elaboracao do EP:
# - O problema possui pelo menos uma solucao viavel basica
# - Nao ha solucoes viaveis basicas degeneradas
#=======================================================================


#=======================================================================
# Teste com uma solucao viavel basica otima: (extraida do livro p.88)
# A = [1,1,1,1 ; 2,0,3,4]  x = [1;1;0;0]   b = [2;2]
# c = [2,0,0,0]  m = 2  n = 4
#
#
# Teste p.55
# A = [1,1,2,1,0,0,0; 0,1,6,0,1,0,0; 1,0,0,0,0,1,0; 0,1,0,0,0,0,1]
# x = [0,0,0,8,12,4,6]   otima
# c = [2,1,1,0, 0,0,0]  m = 4   n = 7
#=======================================================================

#=======================================================================
# Essa funcao recebe a matriz A, os vetores b e c, o numero de restricoes
# m, o numero de variaveis do problema n, e uma solucao viavel basica x.
#
# Devolve na variavel 'ind' -1 se o problema for ilimitado e em v a direçao
# que leva o custo a -Inf. Devolve 0 no 'ind' se o problema tiver uma 
# solucao otima, armazenada em v. (Exraida da descricao do EP2).
#=======================================================================
function [ind, v] = simplex(A, b, c, m, n, x)
   ind = 2;   #dummy
   count = 0;
   v = zeros(n,1);        

  [basic, Nbasic] = indBasico(x, n);
  [B,cb] = matrizB(A, basic, c, m, n, x);

  invB = inv(B);

  while(ind == 2)
    printIter(count);
    printCost(c, x);

    cost = custo(A, invB, Nbasic, c, cb, m, n, x);
 #   j  = checkNegativo(cost, n) ;

    if (j != 0)                    #Existe uma direcao que diminui o custo
   [ind, v, x] = direction(A, invB, basic, Nbasic,cost, m, n, x);

    else                           #A solucao e otima
      ind = 0;    fim = 1;
      v = x;
    endif
    count++;
  end
  printf("\n");
endfunction



#=======================================================================
# Esta funcao recebe um vetor x de tamanho n e devolve um vetor basic 
# com indices das variaveis basicas e Nbasic com indices das variaves 
# nao basicas de x.
#=======================================================================
function [basic, Nbasic] = indBasico(x, n)
  basic = zeros(n, 1);    k = 1;
  Nbasic = zeros(n, 1);   s = 1;   

  printf("Variaveis basicas:\n");
  for (i = 1 : n)
    if (x(i) == 0)
      Nbasic(k++) = i;

    else
      basic(s++) = i;
      printf("%d:  %f\n", i, x(i));
    endif
  end  
endfunction



#=======================================================================
# Esta funcao recebe uma matriz A de tamanho m x n, vetores c e x de  
# tamanho n e o vetor basic com indices das variaveis basicas de x.
# Devolve uma matriz B com as colunas basicas de A, associado a solucao 
# x e cb, o vetor de custo dos indices basicos.
#=======================================================================
function [B, cb] = matrizB(A, basic, c, m, n, x) 
  k = 1;  index = basic(k);

  B = A(:, index);
  cb = c(index);

  index = basic(++k);
  while(index != 0)
    B = horzcat(B, A(:, index));
    cb = vertcat(cb, c(index));
    index = basic(++k);
  end
endfunction



#=======================================================================
# Esta funcao calcula o vetor de custos reduzidos para cada direcao nao
# basica, onde: 
# p' = cb' * inv(B)   e   (cj* = cj - p' * Aj)
#=======================================================================
function [cost] = custo(A, B, Nbasic, c, cb, m, n, x)

  p = zeros(m,1);
  for (i = 1 : m)
    p(i) = cb' * B(:, i);
  end


  printf("\nCustos reduzidos:\n");
  cost = zeros(n,1);
  k = 1;  index = Nbasic(k);  

  while (index != 0)
    cost(index) =  c(index) - p' * A(:, index);
    printf("%d: %f\n", index, cost(index));
    index = Nbasic(++k);
  end  
endfunction



#=======================================================================
# Esta funcao retorna o indice do primeiro elemento do vetor v tal 
# que v[k] < 0, k = 1...n. Caso nao encontre nenhum elemento, retorna 0.
#=======================================================================
function [fim] = checkNegativo(v, n)
  fim = 0;
  for (k = 1 : n)
    if (v(k) < 0)
      fim = k;
      break;
    endif
  end  
endfunction



#=======================================================================
# Esta funcao determina se a direcao viavel basica que diminui o valor da 
# funcao de custo ou se leva a funcao de custo ao -Inf
#=======================================================================
function [ind, v, x] = direction(A, invB, basic, Nbasic, cost, m, n, x)
  j = checkNegativo(cost, n);  

  if (j == 0)     #A solucao e otima
    ind = 0;      #nao ha direcao basica que diminua o custo
    v = x;
    return;
  endif

  u = invB * A(:, j);
  s = checkNegativo(-u, n);

  printf("\nEntra na base: %d\n", j);
  printDirecao(basic, u);

  if (s != 0)                                   
    [theta, index] = thetaMax(u, basic, m, x);
    [B] = swap(A, B, basic, index, j, m, n);     #Coluna j de A entra da base
                                                 #Coluna index de B sai na base
    x = newX(theta, u, basic, index, j, x, m, n);
    ind = 2;  #dummy
    v = x;    fim = 0;
      
  else                                           #custo infinito 
    v = u;
    ind = -1;
    return;
  endif
endfunction



#=======================================================================
# A funcao calcula thetaMax (o maximo que se pode andar na direcao 
# viavel basica respeitando todas as restricoes) e o indice j tal que 
# thetaMax = min [xB(j) / u(j)], j = 1...m, u(j) > 0
#=======================================================================

####Arrumar [xb]i / ui < theta

function [theta, j] = thetaMax(u, basic, m, x)
  t = zeros(m,1);
  theta = -Inf;

  k = 1;   index = basic(k);
  while (index != 0)
#**    if (u(k) > 0)
#**       t(index) = (x(index) / u(index));
#*    endif
  index = basic(++k);
  end
t
  [theta, j] = minimo(-t, m);
  printf("\nTheta*: %f\n\n", theta);
  printf("Sai da base: %d\n\n", j);

endfunction



#=======================================================================
# A funcao retorna o indice e o elemento de menor valor do vetor 
# v de tamanho n. 
#=======================================================================
function [min, index] = minimo(v, n)
  min = v(1);
  index = 1;
  
  for (i = 2 : n) 
    if (min > v(i))
      min = v(i);
      index = i;
    endif
  end
endfunction



#=======================================================================
########PASSIVEL DE ERROS(?)
# A funcao troca a coluna B[sai] pela coluna A[entra]
#=======================================================================
function [B] = swap(A, B, basic, sai, entra, m, n)
  k = 1;  index = basic(k); 

  while (index != 0)
    if (index == sai)
      B(:, k) = A(:, entra:entra);        
      break;
    endif

    index = basic(++k);
  end
endfunction



#=======================================================================
# A funcao atualiza a antiga solucao basica pela nova que diminui a 
# funcao de custo.
#=======================================================================
function [x] = newX(theta, db, basic, sai, entra, x, m, n)
  k = 1;  index = basic(k);
  
  while(index != 0)
    x(index) = x(index) + theta * db(index);
    index = basic(++k);
  end

  x(entra) = - theta;
endfunction



#=======================================================================
# Funcoes auxiliares de impressao
#=======================================================================
function printIter(count)
    printf("#====================================\n");
    printf("           Iteracao: %d\n", count);
    printf("#====================================\n");
endfunction


function printCost(c, x)
    printf("Valor funcao objetivo: %f\n", rot90(c) * x);
endfunction


function printDirecao(basic, u)
    k = 1; index = basic(k);
    printf("\nDirecao\n");
    while (index != 0)
      printf("%d: %f\n", index, u(index));
      index = basic(++k);
    end
endfunction

