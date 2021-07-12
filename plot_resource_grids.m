function hDownlinkEstimationEqualizationResults(rxGrid, eqGrid)
    
    % Plot received grid error on logarithmic scale
    figure;
    dims = size(rxGrid);
    surf(20*log10(abs(rxGrid)));
    title('Received resource grid');
    ylabel('Subcarrier');
    xlabel('Symbol');
    zlabel('absolute value (dB)');
    axis([1 dims(2) 1 dims(1) -40 10]);

    % Plot equalized grid error on logarithmic scale
    figure
    surf(20*log10(abs(eqGrid)));
    title('Equalized resource grid');
    ylabel('Subcarrier');
    xlabel('Symbol');
    zlabel('absolute value (dB)');
    axis([1 dims(2) 1 dims(1) -40 10]);

end
