function out = askYesNo(str_prompt)
    YesOrNo = input(sprintf('%s [y/n]',str_prompt),'s');
    if strcmpi(YesOrNo, 'y')
        out = 1;
    else
        out = 0;
    end
end