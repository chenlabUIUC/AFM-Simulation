function zmap_rescale = ZmapAdjust(zmap, adatom_offset, vacancy_offset,adatom_sigma, vacancy_sigma)
% Rescale Zmap for more reasonable intensities.
%
% We assume the original atomic layer gives a Gaussian intensity
% distribution, and adatom and vacancies just deviate slightly from that,
% as observed in real AFM data. We first determine the Gaussian
% distribution of the original atomic layer as 'Gaussian_base', and the
% adatom and vacancy distribution is controled by
% 'adatom_offset', 'adatom_sigma', 'vacancy_offset','vacancy_sigma' as the
% input. For those values larger than mu+3*sigma of Gaussian_base, we
% rescale them as Gaussian_adatom. For those value smaller than mu-3*sigma
% of Gaussian_base, we rescale them as Gaussian_vacancy.


    I = zmap;
    I_squeeze = reshape(I,1,[]);

    [h,c] = hist(I_squeeze,[-3:0.1:10]);
    [xData, yData] = prepareCurveData( c, h );
    % Set up fittype and options.
    ft = fittype( 'gauss1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0];
    opts.StartPoint = [2257 2.7 0.649649520227623];
    [fitresult, gof] = fit( xData, yData, ft, opts );
%     figure()
%     hold on;
%     plot(c,h)
%     plot(c, fitresult(c))

    % Prepare 3 Gaussian distribution
    Gaussian_base = [fitresult.b1, fitresult.c1];
    Gaussian_adatom = Gaussian_base + [adatom_offset, adatom_sigma*Gaussian_base(2)];
    Gaussian_vacancy = Gaussian_base + [vacancy_offset, vacancy_sigma*Gaussian_base(2)];

    % Rescaling for adatoms
    select_adatom = I_squeeze > (Gaussian_base(1)+3*Gaussian_base(2));
    selectData = I_squeeze(select_adatom);
    [sortedData,ind] = sort(selectData);
        while 1
            % Make sure no duplicate values
            select = sortedData(2:end) - sortedData(1:end-1);
            select = find(select==0)+1;
            if isempty(select); break; end
            sortedData(select) = sortedData(select)+eps;
            inverseIndex = zeros(1,length(ind));
            inverseIndex(ind) = 1:length(ind);
            selectData = sortedData(inverseIndex);
        end
    n = sum(select_adatom);
    quantiles = (1:n) / (n+1);
    standardQuantiles = norminv(quantiles, Gaussian_adatom(1), Gaussian_adatom(2));
    rescaledData = interp1(sortedData, standardQuantiles, selectData, 'linear', 'extrap');
    I_squeeze(select_adatom) = rescaledData;

    % Rescaling for vacancy
    select_vacancy = I_squeeze < (Gaussian_base(1)-3*Gaussian_base(2));
    selectData = I_squeeze(select_vacancy);
    [sortedData,ind] = sort(selectData);
        while 1
            % Make sure no duplicate values
            select = sortedData(2:end) - sortedData(1:end-1);
            select = find(select==0)+1;
            if isempty(select); break; end
            sortedData(select) = sortedData(select)+eps;
            inverseIndex = zeros(1,length(ind));
            inverseIndex(ind) = 1:length(ind);
            selectData = sortedData(inverseIndex);
        end
    n = sum(select_vacancy);
    quantiles = (1:n) / (n+1);
    standardQuantiles = norminv(quantiles, Gaussian_vacancy(1), Gaussian_vacancy(2));
    rescaledData = interp1(sortedData, standardQuantiles, selectData, 'linear', 'extrap');
    I_squeeze(select_vacancy) = rescaledData;


%     [h,c] = hist(I_squeeze,[-3:0.1:10]);
%     [xData, yData] = prepareCurveData( c, h );
%     % Set up fittype and options.
%     ft = fittype( 'gauss1' );
%     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%     opts.Display = 'Off';
%     opts.Lower = [-Inf -Inf 0];
%     opts.StartPoint = [2257 2.7 0.649649520227623];
%     [fitresult, gof] = fit( xData, yData, ft, opts );
%     figure()
%     hold on;
%     plot(c,h)
%     plot(c, fitresult(c))

    zmap_rescale = reshape(I_squeeze,size(zmap,1),[]);
end