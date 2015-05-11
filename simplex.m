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

1;

% Testes
A=[1 2 2 1 0 0 ; 2 1 2 0 1 0; 2 2 1 0 0 1];
b = [20; 20; 20];
c = [-10 -12 -12 0 0 0];
m = 3;
n = 6;
x = [0 0 0 20 20 20];


% Devolve -1 se o custo for ilimitado e a direção em que vai para
% infinito é devolvida em v. 
% Devolve 0 se houver solução ótima em ind e em v devolve o a
% solução ótima encontrada.
function [ind v] = simplex(A,b,c,m,n,x)
    B = encontraBase(A, m, n, x)
    % disp(B)
    red = custosReduzidos(A, b, c, m, n, x, inv(B))
    endfunction
    

% Encontra a matriz base da svb x dada.
function B = encontraBase(A, m, n, x)
    j = 1
    for i = 1:n
        if (x(i) != 0)
            B(:,j) = A(:,i);
            j++;            
            endif;
        endfor;
        endfunction

        
% Calcula o vetor de custos reduzidos
function red = custosReduzidos(A, b, c, m, n, x, B1)
    k = 1
    dmin = 0
    for i = 1:n
        if (x(i) != 0)
            cb(k) = c(i);
            k++;
            endif
        endfor
            
    for j = 1:n
        if (x(j) == 0)
            red(j) = c(j) - cb*B1*A(:,j);
            if (red(j) < red(dmin))
                dmin = j
                endif 
        else 
            red(j) = 0;
            endif

        endfor

    for i = 1:n
	if (red(dmin) == 0)
	    disp('x é solução ótima!');
        else 
	    if (calculaU(B1, A(dmin), m) == -1)
	        disp('O custo ótimo é -infinito!')
            else
                y = encontraSolucaoMelhor
	        simplex(A, b, c, m, n, y)
                endif 
            endif
        endfor
 
    endfunction

function [u ind] = calculaU(B1, Aj, m)
    for i = 1:m
        u = B1*Aj;
        if (u(i) <= 0)
            ind = -1;
            return;
            endif
        endfor
    ind = 0;
    return;

    endfunction

function y = encontraSolucaoMelhor(xb, u, m)
    theta = xb(1)\u(1);
    for i = 2:m
        new = xb(i)\u(i)
            if (new < theta)
	        theta = new;
                l = i;
                endif
        endfor
    for i = 1:n
        if (x(i) > 0)
	    yb(i) = xb(i)-theta*u(i)
            endif
        endfor
    endfunction
        
[ind v] = simplex(A,b,c,m,n,x)
clear all
