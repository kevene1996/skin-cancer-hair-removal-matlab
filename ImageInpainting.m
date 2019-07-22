% Image inpainting is a technique that replaces a pixel with the average
% weighted value of the neighboor pixels, in this case the hair
% representet by the mask is replaced by the neighboor pixels or colors not
% included in the mask. Image inpainting will be implemented using the 
% mumford shah algorithm that is used to segment images.
% Image inpainting works by adding the orginal image with a mask
% representing the objects that want to be removed.
function ImageInpainting = ImageInpainting(Image,Imagemask)
ImagetoInpaint = Imagemask.*Image; % multiply original image with mask
imshow(ImagetoInpaint)
maxiter       = 5;     % number of iterations for inpainting 
tol           = 1e-14;
param.lambda  = 10^9;   % weight on data fidelity (fixed).
param.alpha   = 1;      % regularisation parameters \alpha.(fixed)
param.gamma   = 0.5;    % regularisation parameters \gamma.(fixed)
param.epsilon = 0.05;   % accuracy of approximation of the edge set.(usually 0.05)
% call for mumford shah algorithm
ImageInpainting = inpainting_mumford_shah(ImagetoInpaint,Imagemask,maxiter,tol,param);
end

function [u, chi] = inpainting_mumford_shah(input,mask,maxiter,tol,param)

% import input and mask, M and N are the size parameters of the image
% C is the number of color channels which is equal to 3
[M,N,C] = size(input);
% GRID INTERVAL FOR AXIS ij
h1 = 1; h2 = 1;

% Create forward and backward diagonal matrixes with (-1,1)
% creates an M-by-N sparse matrix from the columns of B and places 
% them along the diagonals specified by d. --> spdiags(B,d,M,N)
d1i_forward  = spdiags([-ones(M,1),ones(M,1)],[0,1],M,M)/h1; 
d1j_forward  = spdiags([-ones(N,1),ones(N,1)],[0,1],N,N)/h2;
d1i_backward = spdiags([-ones(M,1),ones(M,1)],[-1,0],M,M)/h1;
d1j_backward = spdiags([-ones(N,1),ones(N,1)],[-1,0],N,N)/h2;
% set boundary conditions
d1i_forward(end,[1 end]) = [1 -1]/h1;
d1j_forward(end,[1 end]) = [1 -1]/h2;
d1i_backward(1,[1 end]) = [-1 1]/h1;
d1j_backward(1,[1 end]) = [-1 1]/h2;
% get the kroneker tensor product for every matrix
matrices.Dif  = kron(speye(N),d1i_forward);
matrices.Dib  = kron(speye(N),d1i_backward);
matrices.Djf  = kron(d1j_forward,speye(M));
matrices.Djb  = kron(d1j_backward,speye(M));
% define center of the matrix
matrices.Dic = (matrices.Dif+matrices.Dib)/2;
matrices.Djc = (matrices.Djf+matrices.Djb)/2;
% obtain the toeplitz matrix
d2i = toeplitz(sparse([1,1],[1,2],[-2,1],1,M))/h1^2;
d2j = toeplitz(sparse([1,1],[1,2],[-2,1],1,N))/h2^2;
% set periodic boundaries
d2i(1,2)       = 2/h1^2;
d2i(end,end-1) = 2/h1^2;
d2j(1,2)       = 2/h2^2;
d2j(end,end-1) = 2/h2^2;
% apply laplace transform to the matrixes
matrices.LAP = (kron(speye(N),d2i)+kron(d2j,speye(M)));
% FREE MEMORY
clear d1i_forward dji_forward d1i_backward dji_backward
clear d2i d2j

% mumford shah ALGORITHM
u   = reshape(input,[],C); % make sure that matrix is M-by-N
chi = reshape(mask,[],C); % make sure that matrix is M-by-N
param.lambda = param.lambda*chi;
rhsL = (param.lambda/param.gamma).*reshape(u,[],C);
rhsM = ones(M*N,C); 
% for each color channel
for c = 1:C
    % iteration on every pixel where every pixel will be changed by the
    % weighted sum of the neighboor pixels
    for iter = 1:maxiter
        MM       = matrixM(param,M*N,matrices,u(:,c));
        chinew   = MM\rhsM(:,c);
        diff_chi = norm(chinew-chi(:,c))/norm(chinew);
        chi(:,c) = chinew;
        clear chinew
        LL     = matrixL(param,M*N,matrices,chi(:,c));
        unew   = LL\rhsL(:,c);
        diff_u = norm(unew-u(:,c))/norm(unew);
        u(:,c) = unew;
        clear unew
        if diff_u<tol
            break
        end    
        
    end    
end
u   = reshape(u,M,N,C);
chi = reshape(chi,M,N,C);
return
end

% Auxiliary functions for the algorithm
function M = matrixM(param,N,matrices,u)
% Definition of (\nabla u)^2:
nablau2 = (matrices.Dic*u).^2 + (matrices.Djc*u).^2;
M = speye(N)...
    + 2 * param.epsilon * param.gamma/param.alpha * spdiags(nablau2,0,N,N)...
    - 4*param.epsilon^2*matrices.LAP;
return
end
function L = matrixL(param,N,matrices,chi)
% Definition of the nonlinear diffusion weighted by \chi^2:
z  = chi.^2 + param.epsilon^2; % coefficient of nonlinear diffusion
zx = matrices.Dic*z;
zy = matrices.Djc*z;
Z  = spdiags(z, 0,N,N);
Zx = spdiags(zx,0,N,N);
Zy = spdiags(zy,0,N,N);
NonlinearDelta = Z*matrices.LAP + Zx*matrices.Dic + Zy*matrices.Djc;
L =  -NonlinearDelta + spdiags(param.lambda/param.gamma,0,N,N);
return
end
