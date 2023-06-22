function ListLLBnew = updateLLB3(ListLLB,z,tolerance)
%updateLLB3 Updates the list of local lower bounds
% This function is an implementation of Algorithm 3 from:
% K. Klamroth et al., On the representation of the search region in
% multi-objective optimization

%% Initialization
% size_LLB = size(ListLLB,2);
p = size(z,1);
if nargin < 4
    tolerance = 1e-6;
end

% A is the set off all search zones {l}+K that contain z
% A_indexlist = all(repmat(z,1,size_LLB)>ListLLB);
A_indexlist = all(bsxfun(@gt,z,ListLLB));
A = ListLLB(:,A_indexlist);
sizeA = size(A,2);

B = cell(1,p);
P = cell(1,p);

%% Update procedure
for j=1:p
    % B(j) contains all search zones whose boundary contains z with respect
    % to the j-th component (u_j = z_j)
%     B{j} = ListLLB(:,all([repmat(z(1:j-1),1,size_LLB)>ListLLB(1:j-1,:);...
%                 repmat(z(j),1,size_LLB)==ListLLB(j,:);...
%                 repmat(z(j+1:p),1,size_LLB)>ListLLB(j+1:p,:)]));
    B{j} = ListLLB(:,all([bsxfun(@gt,z(1:j-1),ListLLB(1:j-1,:));...
                bsxfun(@eq,z(j),ListLLB(j,:));...
                bsxfun(@gt,z(j+1:p),ListLLB(j+1:p,:))]));       
    % P(j) contains the projections of z on all local upper bounds in A
    % along the j-th dimension / component
    P{j} = [A(1:j-1,:);z(j).*ones(1,sizeA);A(j+1:p,:)];
end

% The following loop filters all redundant points out of P(j) for every
% j=1:p where redundant bascially means dominated
P_indexlist = false(p,sizeA);
for j=1:p
    Pj = P{j};
    PB = [Pj B{j}];
%     sizePB = size(PB,2);
    for i=1:sizeA
%         P_indexlist(j,i) = ~any(all([all(repmat(Pj(:,i),1,sizePB) >= PB);any(repmat(Pj(:,i),1,sizePB) > PB)]));
        P_indexlist(j,i) = ~any(all([all(bsxfun(@ge,Pj(:,i),PB));any(bsxfun(@gt,Pj(:,i),PB))]));
    end
    P{j} = Pj(:,P_indexlist(j,:));
end

% P is transformed from a cell array to a matrix
P = [P{1:p}];

% The new list of local upper bounds is computed and made unique
ListLLBnew = unique([ListLLB(:,~A_indexlist),P]','rows')';
% ListLLBnew = [ListLLB(:,~A_indexlist),P];
% sizeLLBnew = size(ListLLBnew,2);
% if sizeLLBnew > 2
%     LLB_indexlist = true(1,sizeLLBnew);
%     for i=1:sizeLLBnew
%         llb = ListLLBnew(:,i);
% %         LLB_indexlist(i) = ~any(all((repmat(llb,1,sizeLLBnew-1)>=ListLLBnew(:,[1:i-1,i+1:sizeLLBnew])-tol/2)));
%         LLB_indexlist(i) = false;
%         LLB_indexlist(i) = (sum(max(abs(ListLLBnew(:,LLB_indexlist)-llb))<tolerance)<p);
%     end
%     ListLLBnew = ListLLBnew(:,LLB_indexlist);
% end
end