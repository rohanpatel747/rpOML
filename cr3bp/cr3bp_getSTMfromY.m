function stm = cr3bp_getSTMfromY(Y, i)
%CR3BP_GETSTMFROMY Returns STM [6x6] given integration output and index (i)

    stm = reshape(Y(i,7:end), 6, 6);

end