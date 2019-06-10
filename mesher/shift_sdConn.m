function [sdConn] = shift_sdConn(sht, sdConn)
sdConnLen = cellfun(@length,sdConn);
sdConnUpdate = cell2mat(sdConn)+sht;
sdConn = mat2cell(sdConnUpdate,sdConnLen,2);