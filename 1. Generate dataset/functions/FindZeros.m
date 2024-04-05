function z = FindZeros(d, mat, init)
% Find the force zero point after getting the force curve.
%
% In more detail, we are finding the closest point with value that has the
% opposite sign from the init value, and find a more accurate zero point
% with linear interpolation.

    n = size(mat,3);
    z = zeros(size(mat,[1,2]));
    i_init = find(abs(d-init)==min(abs(d-init)),1);
    
    for i = 1:size(mat,1)
        for j = 1:size(mat,2)
            l = reshape(mat(i,j,:),[],1);
            l2 = l(1:end-1) .* l(2:end) < 0;
            ids = find(l2);
            
            if ~isempty(ids)
                dists = abs(ids-i_init);
                iids = find(dists == min(dists),1);
                ids = ids(iids);
                
                % refine the force zero point, linear interpolation
                assert(l(ids)*l(ids+1)<0)
                
                slope = (l(ids) - l(ids+1)) / (d(ids) - d(ids+1));
                Fz = d(ids) - l(ids)/ slope;
                z(i,j) = Fz;
                
            else
                disp(['Zero point error: ',num2str(i), ' ', num2str(j)])
            end
            
        end
    end
    
%     [~,z] = min(abs(mat),[],3);
%     z = d(z);
end