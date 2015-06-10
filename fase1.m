%---------------------------------------------
% Simplex - Fase 1 e 2
%
% Recebe um problema de programação linear
% encontra uma solução viável básica inicial 
% caso o problema seja viável, c.c. indica que
% o problema é inviável.
% Devolve a solução ótima caso exista, c.c.
% devolve a direção para a qual o custo ótimo
% é infinito.
%---------------------------------------------

%---------------------------------------------
% Só é necessária uma função que calcula a 
% solução ótima.
% 
%---------------------------------------------
1;
%---------------------------------------------
% Função inicial
%---------------------------------------------
function [ind x d] = simplex(A, b, c, m, n)
    % INICIO DA FASE 1

    % Organiza o problema original
    for i = 1: m
        if (b(i) < 0) 
            A(i,:) = -1 * A(i,:);
            b(i) = -1 * b(i);
        endif
    endfor

    % Acrescenta as váriaveis artificiais
    newA = [A eye(m)]
    newx = [linspace(0, 0, n), b']'
    newcb =  [linspace(1, 1, m)]'
    newc = [linspace(0, 0, n), newcb']'

    % Calcula a solução ótima do problema auxiliar
    [tgt, u] = simplex(newA, b, newc, m, n+m, newx)
    
    % Verifica a viabilidade
    if (tgt == 0)
        if (newc'*u > 0)
            ind = 1
            return
        endif
    endif
    
    % O problema é viável, retira restrições extras, se existirem
    % Remove as variáveis artificiais
    for i = n: n+m
        if (u(i) == 0)
            newA = [newA(1:(i-1)-n, :); newA((i+1)-n:m, :)]
            newx = [newx(1:(i-1)-n), newx((i+1)-n: end)]
        endif
 
    newA = [newA(,:1:i-1), newA(,: i+1:end)]
    newx = [newx(1: i-1), newx(i+1: end)]
    
    [m n] = size(newA)
    % FIM DA FASE 1
    
    % INICIO DA FASE 2
    [ind, u] = simplex(newA, b, c, m, n, newx)
    %FIM DA FASE 2
endfunction
