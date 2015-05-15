%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Passos do algoritmo: 
% 1 - Recebe os dados e encontra a matriz B (basicamente é ver
% quais elementos não são zero e pegar a coluna de A
% correspondente);
% 2 - Calcula o vetor de custos reduzidos para toda variável não
% básica. Se nenhuma componente é negativa, então a svb atual é
% ótima, acabou. Caso contrário, escolha algum j, onde c_j <0 (de preferência
% o menor?);
% 3 - Calcule u = B^(-1)A_j. Se nenhum componente de u é positivo,
% temos que o custo ótimo é - infinito, acabou;
% 4 - Caso contrário, tomamos theta* = min{x_(B(i))/u_i}
% 5 - Seja l o indice onde o minimo foi encontrado. Forma uma nova
% base substituindo A_(B(l)) por A_j. Se y é a nova svb, os valores
% das novas variáveis básicas são y_j = theta*, y_B(i) - theta*u_i, i!=l. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[1 2 2 1 0 0 ; 2 1 2 0 1 0; 2 2 1 0 0 1];
b = [20; 20; 20];
c = [-10; -12; -12; 0; 0; 0];
m = 3;
n = 6;
x = [0; 0; 0; 20; 20; 20];

%----------------------------------------------------------------
% Devolve -1 se o custo for ilimitado e a direção em que vai para
% infinito é devolvida em v. 
% Devolve 0 se houver solução ótima em ind e em v devolve o a
% solução ótima encontrada.
%----------------------------------------------------------------
function [ind, v] = simplex(A,b,c,m,n,x)
    var = tipoVar(x, n);
    B = encontraBase(A, m, n, var);
    invB = inv(B);
    red = custosReduzidos(A, b, c, m, n, x, invB);
    
endfunction
    

%----------------------------------------------------------------
% Recebe uma solução viável básica e seu tamanho e retorna um 
% vetor que contém 1 nos índices das variáveis básicas e 0 nos 
% índices das variáveis não básicas.
%----------------------------------------------------------------
function var = tipoVar(x, n)
    for i = 1:n
        if (x(i) != 0)
            var(i) = 1;
        else 
	    var(i) = 0;
        endif 
    endfor
endfunction



%-----------------------------------------------------------------
% Recebe a matriz A, um vetor v descrevendo os tipos das váriáveis
% e o seu tamanho n, e encontra a matriz base da solução viável
% básica relacionada.
%-----------------------------------------------------------------
function B = encontraBase(A, m, n, v)
    j = 1
    for i = 1:n
        if (v(i) != 0)
            B(:,j) = A(:,i);
            j++;            
        endif;
     endfor;
 endfunction

%----------------------------------------------------------------
	% Recebe uma matriz A, 
% Calcula o vetor de custos reduzidos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function red = custosReduzidos(A, b, c, m, n, x, invB, var)
    k = 1;
    dmin = 1;
 
    % Calcula cb
    for i = 1:n
        if (var(i) != 0)
            cb(k) = c(i);
            k++;
        endif
    endfor
            
    
    % Calcula custos reduzidos
    for j = 1:n
        if (var(j) == 0)
            red(j) = c(j) - cb'*invB*A(:,j); %'
            if (red(j) < red(dmin))
	        dmin = j;
            endif 
        endif
    endfor

    for i = 1:n
	if (red(dmin) >= 0)
	    disp('x é solução ótima!');
        else 
            [u, ind] = calculaU(invB, A(dmin), m);
	    if (ind == -1)
	        disp('O custo ótimo é -infinito!');
            else
	        y = encontraSolucaoMelhor(x, u, m, n, var, dmin);
                simplex(A, b, c, m, n, y);
            endif 
        endif
    endfor
endfunction


%--------------------------------------------------------------------
% Recebe a inversa de B e a coluna Aj correspondente a direção viável
% que reduz o custo e calcula o vetor u = invB*Aj. Devolve -1 caso u 
% não possua algum componente positivo, caso contrário, devolve 0.
%--------------------------------------------------------------------
function [u, ind] = calculaU(invB, Aj, m)
    u = invB*Aj;
    
    for i = 1:m
        if (u(i) > 0)
            ind = 0;
            return;
        endif
    endfor
    
    ind = -1;
    return;
endfunction


%----------------------------------------------------------------------
% Recebe uma 
%----------------------------------------------------------------------
function y = encontraSolucaoMelhor(x, u, m, n, var, entra)
    for i = 1:n
	if var(i) != 0;
            xb(i) = x(i);
        endif
    endfor

    theta = xb(1)\u(1);
    for i = 2:m
        new = xb(i)\u(i)
        if (new < theta)
	    theta = new;
            l = i;
        endif
    endfor
    
    for i = 1:m
        xb(i) = xb(i)-theta*u(i)
    endfor

    j = 0;
    for i = 1:n
        if (var(i) != 0)
	  y(i) = xb(j);
          j++;
        endif
    endif

    y(entra) = theta;

endfunction
        
[ind, v] = simplex(A,b,c,m,n,x)
clear all
