% for n=8:1:100
%     for alpha=1:20
%         for beta=1:20
%             abfoa2(@testSphereFunc, 10, 0, 1,n,alpha,beta);
%         end
%     end
% end

abfoa2(@testSphereFunc, 10, 0, 1,100,10,10);
