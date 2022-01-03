function throwUndefinedTypeError(this, type)
    if(isnumeric(type))
        str =num2str;
    elseif(isenum(type))
        str1 =class(type);
        str2 =type.char;
        str = sprintf('%s.%s',str1,str2);
    else
        str = 'Unknown Type';
    end
    fprintf('[Error] Undefined Type ''%s'' in: %s\n', str, class(this));
    error(sprintf('Undefined Type ''%s'' in: %s', str, class(this)));
end