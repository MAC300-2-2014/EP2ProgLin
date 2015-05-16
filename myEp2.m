1;

#=======================================================================
# As seguintes suposicoes foram feitas para a elaboracao do EP:
# - O problema possui pelo menos uma solucao viavel basica
# - Nao ha solucoes viaveis basicas degeneradas
#=======================================================================


#=======================================================================
# Teste com uma solucao viavel basica otima: (extraida do livro p.88)
# A = [1,1,1,1 ; 2,0,3,4]  x = [1;1;0;0]   b = [2;2]
# c = [2;0;0;0]  m = 2  n = 4
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

  [basic, Nbasic] = indBasico(x, n);
  [B, InvB] = matrizB(A, basic, m, n, x);

  while(ind == 2)
    printIter(count);
    printBasic(basic, x);
    printCost(c, x);
    
    cb = findCb(basic, c);
    cost = custo(A, InvB, Nbasic, c, cb, m, n, x);
    [ind, v, j] = direction(A, InvB, cost, m, n, x);

    if (j != 0)
      [B, InvB, basic, Nbasic, x] = newB(A, B, InvB, basic, Nbasic, v, j, x, m, n);
    endif


    count++;
  end

  printSolucao(x, v, basic, c, ind, m, n);
endfunction



#=======================================================================
# Esta funcao recebe um vetor x de tamanho n e devolve um vetor basic 
# com indices das variaveis basicas e Nbasic com indices das variaves 
# nao basicas de x.
#=======================================================================
function [basic, Nbasic] = indBasico(x, n)
  basic = zeros(n, 1);    k = 1;
  Nbasic = zeros(n, 1);   s = 1;   

  for (i = 1 : n)
    if (x(i) == 0)
      Nbasic(k++) = i;

    else
      basic(s++) = i;
    endif
  end  
endfunction



#=======================================================================
# Esta funcao recebe uma matriz A de tamanho m x n, vetor x de  
# tamanho n e o vetor basic com indices das variaveis basicas de x.
# Devolve uma matriz B com as colunas basicas de A associado a solucao 
# x e a matriz inversa de B 
#=======================================================================
function [B, InvB] = matrizB(A, basic, m, n, x) 
  k = 1;  index = basic(k);

  B = A(:, index);
  index = basic(++k);
  while(index != 0)
    B = horzcat(B, A(:, index));
    index = basic(++k);
  end
 
  InvB = inv(B);
endfunction


#=======================================================================
# Esta funcao determina cb, o vetor de custo dos indices basicos.
#=======================================================================
function [cb] = findCb(basic, c) 
  k = 1;  index = basic(k);

  cb = c(index);
  index = basic(++k);
  while(index != 0)
    cb = vertcat(cb, c(index));
    index = basic(++k);
  end
 
endfunction



#=======================================================================
# Esta funcao calcula o vetor de custos reduzidos para cada direcao nao
# basica, onde:  p' = cb' * inv(B)   e   (cj* = cj - p' * Aj)
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
# Esta funcao devolve em v a solucao otima caso x seja otimo. 
# Se x nao for otimo, devolve o indice (k) da variavel nao basica que
# diminua o valor da funcao de custo. (k = 0 se esta indice nao existir)
#=======================================================================
function [ind, v, k] = direction(A, invB, cost, m, n, x)
  ind = 2;  #dummy
  k = checkNegativo(cost, n);  

  if (k == 0)     #A solucao e otima
    ind = 0;  v = x;
    j = 0;
    return;
  endif
                  
  v = invB * A(:, k);            
  j = checkNegativo(-v, n);      

  if (j == 0)
    printf("O valor da funcao objetivo vai para -Inf\n");
    ind = -1;
    return;
  endif


endfunction


#=======================================================================
# Esta funcao gera a matriz B atualizada e associada com a nova solucao 
# x. Atualiza tambem os indices das novas variaveis basicas e nao basicas.
#=======================================================================
function [B, InvB, basic, Nbasic, x] = newB (A, B, InvB, basic, Nbasic, u, j, x, m, n)
  printf("\nEntra na base: %d\n", j);
  printDirecao(basic, u, m);

  [theta, index] = thetaMax(u, basic, m, x);
  x = newX(theta, basic, u, j, x, m, n);
  [B, basic, Nbasic] = swap(A, B, basic, Nbasic, index, j, m, n);   
                                                  #Coluna index de B sai na base
                                                  #Coluna j de A entra da base
#   InvB = newInvB(invB, v);

  InvB = inv(B);

endfunction



#=======================================================================
# A funcao calcula thetaMax (o maximo que se pode andar na direcao 
# viavel basica respeitando todas as restricoes) e o indice j tal que 
# thetaMax = min [xB(j) / u(j)], j = 1...m, u(j) > 0
#=======================================================================
function [theta, j] = thetaMax(u, basic, m, x)
  theta = Inf;
 
  j = 0;
  k = 1;   index = basic(k);
   
  while (index != 0)
    if (u(k) > 0)
      t = (x(index) / u(k));
      
      if (t < theta )
	theta = t;
        j = index;
      endif
    endif
  index = basic(++k);
  end

  printf("\nTheta*: %f\n\n", theta);
  printf("Sai da base: %d\n\n", j);
endfunction


#=======================================================================
# A funcao troca a coluna B[sai] pela coluna A[entra] e atualiza os
# indices basicos e nao basicos
#=======================================================================
function [B, basic, Nbasic] = swap(A, B, basic, Nbasic, sai, entra, m, n)
  k = 1;  index = basic(k); 

  while (index != sai)
    index = basic(++k);
  end

  B(:, k) = A(:, entra);        
  basic(k) = entra;
  k = 1; index = Nbasic(k);
  while (index != entra)
   index = Nbasic(++k);
  end
  Nbasic(k) = sai;


endfunction



#=======================================================================
# A funcao atualiza a antiga solucao basica pela nova que diminui a 
# funcao de custo.
#=======================================================================
function [x] = newX(theta, basic, u, entra, x, m, n)

  k = 1; index = basic(k);  
  while (index != 0)
    x(index) -= theta * u(k);
    index = basic(++k);
  end

  x(entra) = theta;
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
# Funcoes auxiliares de impressao
#=======================================================================
function printBasic(basic, x)
  printf("Variaveis basicas:\n");
  k = 1; index = basic(k);
  while (index != 0)
    printf("%d:  %f\n", index, x(index));
    index = basic(++k);
  end
    printf("\n");
endfunction


function printIter(count)
    printf("#====================================\n");
    printf("           Iteracao: %d\n", count);
    printf("#====================================\n");
endfunction



function printDirecao(basic, u, m)
  printf("\nDirecao\n");
  for (i = 1 : m)
    printf("%d: %f\n", basic(i), u(i));
  end
endfunction



function printSolucao(x, v, basic, c, ind, m,  n)
  if (ind == 0)
    printf("\nSolucao otima com custo %f:\n", c' * x);
    for (i = 1 : n)
       printf("%d  %f\n",i ,x(i));
    end

  else
    printf("\nA direcao que leva o custo a -Inf:\n");
    printDirecao(basic, u, m);
  endif

endfunction



function printCost(c, x)
    printf("Valor funcao objetivo: %f\n", c' * x);
endfunction



