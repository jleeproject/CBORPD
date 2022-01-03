function index = fnGetIdxOfSubplotWithRowCol(subplotSize, row, col)
%     [spRow, spCol] = size(subplotSize);
    spRow = subplotSize(1);
    spCol = subplotSize(2);
    if(row>spRow || col >spCol)
        disp('ERROR: row or column is larger than the size of subplot');
        fprintf('index:[%d, %d], plot:[%d, %d]\n', row, col, spRow, spCol);
        index = [];
        return;
    end
    index = (row-1)*spCol + col;
end