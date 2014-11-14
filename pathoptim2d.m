clf
clear all
theta = pi/2;
Rz90 = [cos(theta) -sin(theta);
sin(theta) cos(theta)];

M_ws = 14;
t = linspace(0,1,M_ws);
F_ws = [t.^0;t.^1;t.^2;t.^3;t.^4;t.^5;t.^6;t.^7;t.^8;t.^9;t.^10];
dF_ws = [0*t; 1*t.^0; 2*t.^1; 3*t.^2; 4*t.^3; 5*t.^4; 6*t.^5; 7*t.^6; 8*t.^7;9*t.^8;10*t.^9];
K_ws = size(F_ws,1);

Mk = 2;
t = linspace(0,1,Mk);
Fk= [t.^0;t.^1;];
Kk= size(Fk,1);

Loffset=-200;
maxDistKnee=0.6
cvx_begin 
        variable X(2,M_ws)
        variable W_ws(K_ws,2)

        variable Z(2,3,M_ws)
        variable W_wk(Kk,2,M_ws)
        variable W_wk90(Kk,2,M_ws)

        variable lambda(M_ws+1)

        minimize sum(lambda)
        
        subject to
                %%% objective value epigraph minimizations
                for i=1:M_ws
                        norm(Z(:,2,i)-Z(:,3,i))+norm(Z(:,2,i)-Z(:,1,i)) <= lambda(i)
                end
                norm(W_ws) <= lambda(M_ws+1)

                %%% orthogonal distance constraint
                for i=1:M_ws
                        j=i+1;
                        distLambda=-Loffset;
                        if i==M_ws
                                j=i-1;
                                distLambda=Loffset;
                        end
                        X(:,i) == W_wk(:,:,i)'*Fk(:,1)
                        X(:,j) == W_wk(:,:,i)'*Fk(:,2)
                        W_wk90(:,:,i) == [W_wk(1,1,i) W_wk(1,2,i); W_wk(2,2,i) -W_wk(2,1,i)]

                        Z(:,3,i) == W_wk90(:,:,i)'*[Fk(1,2);distLambda]
                        norm(Z(:,2,i) - Z(:,1,i)) <= maxDistKnee
                        Z(:,1,i) == W_wk90(:,:,i)'*[Fk(1,2);0]
                end

                %%% Constraints on ground points X
                for i=1:M_ws-1
                        norm(X(:,i)-X(:,i+1)) <= 0.8
                end

                X(:,:) == W_ws'*F_ws
                X(:,1) == [1;1]
                X(:,5) == [3;1]
                X(:,10) == [3;3]
                X(:,M_ws) == [5;4]
cvx_end
%% evaluate the orthogonality of the knee vectors
dsum = 0
for i=1:M_ws-1
        if norm(X(:,i)-Z(:,1,i))<0.01
                z = Z(:,2,i)-Z(:,1,i);
                x = X(:,i+1)-Z(:,1,i);
                d = z'*x;
                dsum = dsum+d;
                if d > 0.001
                        disp(sprintf('Warning: non-orthogonality for i=%d', i))
                end
        else
                disp(sprintf('Warning: Z and X values not aligned for i=%d', i))
        end
end
disp(sprintf('Orthogonality SCORE=%f (should be 0)', dsum))

co=cvx_optval;
if co < +inf
        figure(1)
        plot(X(1,:),X(2,:),'*-r')
        hold on;
        for i=1:size(Z,3)
                plot(Z(1,1:2,i),Z(2,1:2,i),'*-b')
                hold on;
        end
        axis equal
end
