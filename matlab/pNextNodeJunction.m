% Run through nodes until a hairpin is reached, then return the node number
% Returns the index of the last hairpin after N
% When a junction is encountered, call recursively and record the next node values

function [Node,n] = pNextNodeJunction(Node,N)

    n = N;
    while n <= length(Node)
        if strcmp(Node(n).type,'Junction')
            fprintf('pNextNodeJunction: Found a junction at node %d\n',n);
            Node(n).nextnode(1) = n+1;
            k = n;
            fprintf('pNextNodeJunction: Found a junction child at node %d\n',k+1);
            for j = 2:length(Node(n).nextnode)
                [Node, k] = pNextNodeJunction(Node,k+1);
                Node(n).nextnode(j) = k+1;
                fprintf('pNextNodeJunction: Found a junction child at node %d\n',k+1);
            end
            [Node, k] = pNextNodeJunction(Node,k+1);
            n = k+1;
        elseif strcmp(Node(n).type,'Hairpin')
            fprintf('pNextNodeJunction: Found a hairpin at node %d\n',n);
            return
        else
            n = n + 1;
        end
    end
    n = min(n,length(Node));
    fprintf('pNextNodeJunction: Finished at node %d out of %d, type %s\n',n, length(Node), Node(n).type);
