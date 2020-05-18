fprintf("\n\n------------Resultados del problema DIXMAANA J------------\n");

% f = @(x) dixmaana(x);
% dimensions = [200; 1000]; 
% 
% for index = 1:length(dimensions)
%     n = dimensions(index);
%     
%     x0 = [];
%     for aux = 1:n
%         x0 = [x0; 2];
%     end
%     
%     fprintf("\n\n n = %d \n\n", length(x0)); 
%     
%     fprintf ("\nMétodo lineBGFS (n = %d): \n", n)
%     
%     %Método BGFS:
%     tic; 
%     [xk, iter] = lineBGFS(f, x0, tol, maxIter); 
%     time = toc; 
%     fprintf("\t iter = %d, tiempo = %f \n", iter, time);
%     
%     %Método mRCSR1:
%     fprintf("\nMétodo mRCSR1 (n = %d): \n", n); 
%     tic; 
%     [xk, iter] = mRCSR1(f, x0, tol, maxIter, deltaMax); 
%     time = toc;
%     fprintf("\t iter = %d, tiempo = %f \n", iter, time);
%     
%     %Método lineBGFSLM con memoria cíclica:
%     
%     vecm = [1,3,5,17,29];
%     
%     for auxm = 1:length(vecm)
%         m = vecm(auxm);
%         fprintf("\nMétodo lineBGFSLMCyclic(n = %d, m = %d) \n", n, m);
%         tic 
%         [xk, iter] = lineBGFSLMCyclic(f, x0, tol, maxIter, m); 
%         time = toc; 
%         fprintf("\t iter = %d, tiempo = %f \n", iter, time);
%     end
% end
%     
%     
%     