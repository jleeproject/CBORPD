function pred_result = putin_array_4dim(ind_pred_result, pred_result, idx1, idx2,idx3, idx4 , dim1, dim2, dim3, dim4)
%     putin_array(ind_pred_result, pred_result, 1, 1);
    
    fnames = fieldnames(ind_pred_result);
    dim = max(size(fnames));
    for i= 1:dim
        fname = fnames{i};
        
%         if(strcmp(fname,'nEvaluations'))
%             disp('');
%         end
        
        if(~isfield(pred_result, fname) && nargin==8)
            str_eval = sprintf('pred_result.%s = cell(%d,%d,%d,%d);',fname, dim1, dim2,dim3,dim4);
            eval(str_eval);
        end
        
        
        field = getfield(ind_pred_result,fname);
%         if(prod(size(field))>1)
            str_eval = sprintf('pred_result.%s{%d,%d,%d,%d}=getfield(ind_pred_result,fname);',fname, idx1, idx2,idx3,idx4);
            eval(str_eval);
%         else
%             str_eval = sprintf('pred_result.%s(%d,%d,%d)=getfield(ind_pred_result,fname);',fname, idx1, idx2,idx3);
%             eval(str_eval);
%         end
    end
    pred_result.indRes(idx1,idx2,idx3,idx4) = ind_pred_result;
end