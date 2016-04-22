function xc = msgmInherit(x, inds)

    if (isempty(x))
        
        xc = [];
    else
        
        xc = x(inds(:));
    end
end