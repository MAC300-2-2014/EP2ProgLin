1;
#=======================================================================
# Essa funcao recebe a matriz A, os vetores b e c, o numero de restricoes
# m, o numero de variaveis do problema n.

# Itera sobre o metodo simplex revisado de duas fases, devolvendo em 
# ind -1 se o problema for ilimitado, 0 se o problema tiver soluçao otima 
# e 1 se o problema for inviavel. 

# Se o problema apresentar uma solucao viavel otima, devolve em x. Caso
# seja ilimitado devolve a direcao em d. (Exraida da descricao do EP3).
#=======================================================================
function [ind, x, d] = simplex (A, b, c, m, n)  
  #INICIO DA FASE 1
  phase = 1;

  # Organiza o problema original
  for (i = 1 : m)
    if (b(i) < 0) 
      A(i,:) = -1 * A(i,:);
      b(i) = -1 * b(i);
    endif
  end
  
  # Acrescenta as váriaveis artificiais
  newA = [A, eye(m, m)];
  newx = [linspace(0, 0, n), b']';
  newc = [zeros(n, 1)', ones(m, 1)']';
  
  # Calcula a solução ótima do problema auxiliar
  printPhase(phase);
  [ind, u] = subsimplex(newA, b, newc, m, n + m, newx, phase);

  # Verifica a viabilidade
  if (newc' * u > 0)
    ind = 1;
    u = 0;    x = 0;
    printf("O problema é inviável \n");
    return;
  endif
  
  [A, x, m, n] = resize(newA, u, ind, n);


  #INICIO DA FASE 2
  phase = 2;
  printPhase(phase);
  [ind, u] = subsimplex(A, b, c, m, n, x, phase);

  [x, d] = solucao(ind, u);
endfunction
  
  

#=======================================================================
# Essa funcao recebe a matriz A, os vetores b e c, o numero de restricoes
# m, o numero de variaveis do problema n, uma solucao viavel basica x
# nao degenerada e a fase que se encontra quando foi chamado.

# Devolve na variavel 'ind' -1 se o problema for ilimitado e em v a direçao
# que leva o custo a -Inf. Devolve 0 no 'ind' se o problema tiver uma 
# solucao otima, armazenada em v. (Exraida da descricao do EP2).
#=======================================================================
function [ind, v] = subsimplex(A, b, c, m, n, x, phase)
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
      
    [ind, v, j] = direction(A, InvB, cost, basic, m, n, x);

    if (j != 0)
      [B, InvB, basic, Nbasic, x] = newB(A, B, InvB, basic, Nbasic, v, j, x, m, n);
    endif

    count++;
  end
  
  if (phase == 1)
    [ind, x] = artificialOut(A, B, InvB, b, basic, Nbasic, x, m, n);
  endif
endfunction



#=======================================================================
# Esta funcao recebe um vetor x de tamanho n e devolve um vetor basic 
# com indices das variaveis basicas e Nbasic com indices das variaves 
# nao basicas de x.
#=======================================================================
function [basic, Nbasic] = indBasico(x, n)
  basic = zeros(n + 1, 1);    k = 1;
  Nbasic = zeros(n + 1, 1);   s = 1;   

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
  p = zeros(m, 1);
  for (i = 1 : m)
    p(i) = cb' * B(:, i);
  end

  printf("\nCustos reduzidos:\n");
  cost = zeros(n,1);
  k = 1;  index = Nbasic(k);  

  while (index != 0)
    cost(index) =  c(index) - p' * A(:, index);
    printf("%d: %f\n", index, cost(index));
    if (cost(index) < 0)
        return;
    endif
    index = Nbasic(++k);
  end  
endfunction



#=======================================================================
# Esta funcao devolve ind = 0, k = 0 e em v a solucao otima caso x seja 
# otimo. Se x nao for otimo, devolve o indice (k) da variavel nao basica 
# que diminua o valor da funcao de custo. Se este indice nao existir, 
# significa que o o valor da funcao objetivo sera -Inf (a funcao devolve 
# ind = -1 e k = 0)
#=======================================================================
function [ind, v, k] = direction(A, invB, cost, basic, m, n, x)
  ind = 2;  #dummy
  k = checkNegativo(cost, n);  

  if (k == 0)     #A solucao e otima
    ind = 0;  v = x;
    k = 0;
    return;
  endif
                  
  v = invB * A(:, k);            
  j = checkNegativo(-v, m);      

  if (j == 0)
    printf("\nO valor da funcao objetivo vai para -Inf\n");
    ind = -1;
    k = 0;
    v = setV(basic, v, n);
    return;
  endif
endfunction



#=======================================================================
# Esta funcao gera a matriz B atualizada e associada com a nova solucao 
# x e InvB da nova matriz (modo revisado). Atualiza tambem os indices 
# das novas variaveis basicas e nao basicas em basic e Nbasic.
#=======================================================================
function [B, InvB, basic, Nbasic, x] = newB (A, B, InvB, basic, Nbasic, u, j, x, m, n)
  printf("\nEntra na base: %d\n", j);
  printDirecao(basic, u, m);

  [theta, index] = thetaMax(u, basic, m, x);
  x = newX(theta, basic, u, j, x, m, n);

  [B, basic, Nbasic] = swap(A, B, basic, Nbasic, index, j, m, n);   
                                                  #Coluna index de B sai na base
                                                  #Coluna j de A entra da base

  for (i = 1 : m)
    if (basic(i) == j)
      k = i;
      break;
    endif
  end


  #Calculo do inverso do novo B
  D = horzcat(InvB, u);

  D(k,:) = (D(k,: ) / D(k,m+1));
  for (i = 1 : m)
    if (i != k)
     D(i,:) =  D(i,:) - (D(k,:) * D(i,m+1));
    endif
  end

  D(:, m+1) = [];
  InvB = D;
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
    if (abs(u(k)) > 1e-4 && u(k) > 0)
      t = (x(index) / u(k));
      
      if (t < theta)
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
# Esta funcao retorna o indice primeiro elemento do vetor v tal 
# que v[k] < 0, k = 1...n. Caso nao encontre nenhum elemento, retorna 0.
#=======================================================================
function [fim] = checkNegativo(v, m)
  fim = 0;
  value = Inf;
  
  for (k = 1 : m)
    if (v(k) < 0  && abs(v(k)) > 1e-4)
      fim = k;
      return
    endif
  end  
endfunction



#=======================================================================
# Esta funcao recupera a dimencao original de u 
#=======================================================================
function [v] = setV(basic, v, n);
  s = zeros(n, 1);
  
  for (i = 1 : n)
    index = basic(i);

    if (index != 0)
      s(index) = v(i);
    else
      v(i) = 0; 
    endif
  end
  v = s;
endfunction



#=======================================================================
# Esta funcao retira as variaveis artificiais basicas caso a i-esima linha
# da variavel xi nao seja nula. Devolve os indices das variaveis artificiais 
# basicas a serem eliminadas de A em ind
#=======================================================================
function [ind, x] = artificialOut(A, B, InvB, b, basic, Nbasic, x, m, n)
  artificial = checkArt(basic, m, n);

  if (artificial == 0)                #Nao ha variavel artificial na base
    ind = 0;
    x = x;
    return;
  endif

  while (artificial != 0) 
    k = 1; id = basic(k);
    while (id != artificial)
      id = basic(++k);
    end

    index = 0;  j = 1;
    #Verifica se a linha da variavel artificial basica xi e nula
    while (index == 0 && j <= n)

      v = InvB * A(:, j);
      if (v(k) > 1e-4 && j != artificial)
	index = j;
      endif
      j++;
    end

                                #index entra na base no lugar de artificial
    if (index == 0)             #linha redundante e pode ser eliminado
       B(k, :) = [];
       basic(k) = [];
       A(k, :) = [];				  
       A(:, artificial) = [];				  
       n = n - 1;
       m = m - 1;
       ind(end + 1) = k;   
    endif

  artificial = checkArt(basic, m, n);
  end
endfunction



#=======================================================================
# Checa se a solucao possui variavel artificial na base
# Retorna o menor indice referente a variavel artificial na matriz base. 
# Caso nao encontre, retorna zero.
#=======================================================================
function [artificial] = checkArt(basic, m, n)
  artificial = n + m;

  k = 1; index = basic(k);

  while (index != 0)
    if (index > n - m  && index < artificial)
	artificial = index;
    endif
  index = basic(++k);
  end
   
  if (artificial == n + m)
    artificial = 0;
  endif
endfunction



#=======================================================================
# Esta funcao elimina as linhas redundantes de A (se existir) indicadas 
# no vetor ind e retorna a dimensao correta da matriz A e do vetor x 
# para a segunda fase do simplex
#=======================================================================
function [A, x, m, n] = resize(newA, u, ind, n)
  A = newA(:, 1);
  x = u(1);  

  for (i = 2 : n)
    A = horzcat(A, newA(:, i));
    x = vertcat(x, u(i));
  end

  n = columns(ind);
  for (i = 1 : n)
    index = ind(i);
    if (index > 0)
       A(index, :) = [];				  
    endif
  end

  m = rows(A);
  n = columns(A);
endfunction



#=======================================================================
# Esta funcao indica se a solucao u encontrada e otima ou devolve a
# direcao para o qual o custo e -Inf 
#=======================================================================
function [x, d] = solucao(ind, u);
  if(ind == 0) 
    printf("\nSolucao e otima \n");
    x = u;
    d = 0;

  else
    printf("\nSolucao possui custo -Inf\n");
    x = 0;
    d = u;
  endif
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


function printCost(c, x)
    printf("Valor funcao objetivo: %f\n", c' * x);
endfunction


function printPhase(phase)
    printf("#====================================\n");
    printf("          Simplex Fase %d\n", phase);

endfunction
