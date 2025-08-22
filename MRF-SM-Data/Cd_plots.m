T = readtable('Cd.csv','PreserveVariableNames',true);

%% --- Prepare columns and data (flexible) ---
vars = T.Properties.VariableNames;
nv = numel(vars);

% Preferred columns by name
haveType = any(strcmpi(vars,'Type'));
haveMesh = any(strcmpi(vars,'Mesh'));

% If Type/Mesh do not exist, use columns 2 and 3 if available
useCol2 = 2; useCol3 = 3;

% Main columns (Results, Data, Comparative, RelativeError)
lowerVars = lower(vars);
findLike = @(pat) find(contains(lowerVars,lower(pat)));

iRes = findLike('result'); 
assert(~isempty(iRes),'Column "Results" not found.');

iVal = findLike('data'); 
if isempty(iVal), iVal = findLike('value'); end
assert(~isempty(iVal),'Column with data values not found (e.g., "Data").');

iComp = findLike('comparative'); 
assert(~isempty(iComp),'Column "Comparative" not found.');

iErr = findLike('relativeerror'); 
hasErr = ~isempty(iErr);

results = string(T{:, iRes(1)});
values  = double(T{:, iVal(1)});
comparative_raw = T{:, iComp(1)};
if hasErr
    relErrRaw = double(T{:, iErr(1)});
else
    relErrRaw = nan(height(T),1);
end

% Normalize Comparative column to numeric (0,1,2)
comp = nan(height(T),1);
if isnumeric(comparative_raw)
    comp = double(comparative_raw);
else
    s = string(comparative_raw);
    comp(s=="0" | s=="0.0") = 0;
    comp(s=="1" | s=="1.0") = 1;
    comp(s=="2" | s=="2.0") = 2;
end

% Experimental values (Cd and Cl)
expCdIdx = comp==0 & (contains(results,'Cd','IgnoreCase',true) | contains(results,'C_d','IgnoreCase',true));
expClIdx = comp==0 & (contains(results,'Cl','IgnoreCase',true) | contains(results,'C_l','IgnoreCase',true));

expCd = NaN; expCl = NaN;
if any(expCdIdx), expCd = values(find(expCdIdx,1,'first')); end
if any(expClIdx), expCl = values(find(expClIdx,1,'first')); end

%% --- Visual parameters ---
maxShortLabelLen = 18;
expColor = [0.7 0.7 0.7];
simColorMap = @lines;
fontSz = 12;

% Create separate figures:
createSeparateFigure(1, 'C_d', results, values, comp, relErrRaw, T, haveType, haveMesh, useCol2, useCol3, expCd, expColor, simColorMap, maxShortLabelLen, fontSz);
createSeparateFigure(1, 'C_l', results, values, comp, relErrRaw, T, haveType, haveMesh, useCol2, useCol3, expCl, expColor, simColorMap, maxShortLabelLen, fontSz);
createSeparateFigure(2, 'C_d', results, values, comp, relErrRaw, T, haveType, haveMesh, useCol2, useCol3, expCd, expColor, simColorMap, maxShortLabelLen, fontSz);
createSeparateFigure(2, 'C_l', results, values, comp, relErrRaw, T, haveType, haveMesh, useCol2, useCol3, expCl, expColor, simColorMap, maxShortLabelLen, fontSz);

disp('Done: 4 separate figures generated (if data available).');

%% -------------------- FUNCTIONS --------------------
function createSeparateFigure(cmpVal, metricKey, results, values, comp, relErrRaw, T, haveType, haveMesh, useCol2, useCol3, expVal, expColor, simColorMap, maxShortLabelLen, fontSz)
    rows = find(comp==cmpVal & (contains(results, metricKey,'IgnoreCase',true) | contains(lower(results),lower(metricKey))));
    figName = sprintf('Comparative %d - %s', cmpVal, metricKey);
    if isempty(rows)
        fprintf('No data for %s (Comparative %d). Skipping.\n', metricKey, cmpVal);
        return
    end
    fh = figure('Name',figName,'NumberTitle','off','Units','normalized','Position',[0.2 0.2 0.5 0.5]);
    plotMetricSeparate(rows, metricKey, expVal, results, values, relErrRaw, T, haveType, haveMesh, useCol2, useCol3, expColor, simColorMap, maxShortLabelLen, fontSz, figName);
end

function plotMetricSeparate(rows, metricName, expVal, results, values, relErrRaw, T, haveType, haveMesh, useCol2, useCol3, expColor, simColorMap, maxShortLabelLen, fontSz, figName)
    N = numel(rows);
    simVals = values(rows);
    simRel = relErrRaw(rows);
    mat = [repmat(expVal,N,1), simVals];
    
    % Labels
    fullNames = cell(N,1);
    shortNames = cell(N,1);
    for k=1:N
        fullNames{k} = labelFromCols(rows(k), T, haveType, haveMesh, useCol2, useCol3, results, maxShortLabelLen);
        shortNames{k} = shortenLabel(fullNames{k}, maxShortLabelLen);
    end
    
    cmap = simColorMap(N);
    hb = bar(mat,'grouped','BarWidth',0.85); hold on;
    ax = gca; ax.FontSize = fontSz; ax.TickLabelInterpreter = 'latex';
    if numel(hb) >= 2
        hb(1).FaceColor = expColor; hb(1).EdgeColor='k'; hb(1).LineWidth=1.0;
        hb(2).FaceColor = 'flat'; hb(2).CData = cmap; hb(2).EdgeColor='k'; hb(2).LineWidth=0.9;
    end
    set(gca,'XTickLabel',shortNames,'XTickLabelRotation',25,'TickLabelInterpreter','latex');
    ylabel(sprintf('$%s$', metricName),'Interpreter','latex'); 
    title(figName,'Interpreter','latex');
    grid on; box on;
    
    if numel(hb) >= 2
        xExp = hb(1).XEndPoints; yExp = hb(1).YEndPoints;
        xSim = hb(2).XEndPoints; ySim = hb(2).YEndPoints;
    else
        xExp = (1:N); yExp = mat(:,1);
        xSim = (1:N); ySim = mat(:,2);
    end
    
    % Absolute errors
    absErr = zeros(N,1);
    for k=1:N
        r = simRel(k);
        if isnan(r), absErr(k)=0;
        else
            if abs(r)>1, pct = r/100; else pct = r; end
            absErr(k) = abs(expVal * pct);
        end
    end
    
    % Error bars
    valid = absErr>0;
    if any(valid)
        hErr = errorbar(xSim(valid), ySim(valid), absErr(valid), absErr(valid), 'k','LineStyle','none','CapSize',8,'LineWidth',1.2);
        uistack(hErr,'top');
    end
    
    % Simulation value + error text
    for k=1:N
        if absErr(k)>0, offset = absErr(k)*1.05 + 0.01*max(mat(:));
        else offset = 0.02*max(mat(:)); end
        r = simRel(k);
        if ~isnan(r)
            if abs(r)>1
                pctVal = r;
            else
                pctVal = 100*r;
            end
            pctStr = sprintf('%.1f\\%%', pctVal);
        else
            pctStr = '';
        end
        if isempty(pctStr)
            labeltxt = sprintf('%.4g', simVals(k));
        else
            labeltxt = sprintf('%.4g\n%s', simVals(k), pctStr);
        end
        text(xSim(k), ySim(k) + offset, labeltxt, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize',10, 'Interpreter','latex');
    end
    
    % Legend
    hLeg = gobjects(1+N,1); labels = cell(1+N,1);
    if ~isnan(expVal), labels{1} = sprintf('Experimental: %.4g', expVal); else labels{1}='Experimental: n/a'; end
    hLeg(1) = plot(nan,nan,'s','MarkerFaceColor',expColor,'MarkerEdgeColor','k','MarkerSize',10);
    for k=1:N
        hLeg(1+k) = plot(nan,nan,'s','MarkerFaceColor',cmap(k,:),'MarkerEdgeColor','k','MarkerSize',10);
        labels{1+k} = fullNames{k};
    end
    lg = legend(hLeg, labels, 'Location','eastoutside','Interpreter','latex');
    lg.Box = 'on';
    
    % Adjust Y limits
    ymax_candidate = max([mat(:,1); mat(:,2) + absErr]);
    ylim([0 max(0.001, ymax_candidate*1.22)]);
end

function s = labelFromCols(rowIdx, T, haveType, haveMesh, useCol2, useCol3, results, maxShortLabelLen)
    if haveType
        v1 = string(T{rowIdx, find(strcmpi(T.Properties.VariableNames,'Type'),1)});
    elseif size(T,2) >= useCol2
        v1 = string(T{rowIdx, useCol2});
    else
        v1 = "";
    end
    if haveMesh
        v2 = string(T{rowIdx, find(strcmpi(T.Properties.VariableNames,'Mesh'),1)});
    elseif size(T,2) >= useCol3
        v2 = string(T{rowIdx, useCol3});
    else
        v2 = "";
    end
    v1 = strtrim(char(v1)); v2 = strtrim(char(v2));
    if isempty(v1) && isempty(v2)
        s = shortenLabel(results(rowIdx), maxShortLabelLen);
    else
        if isempty(v1), s = v2; elseif isempty(v2), s = v1; else s = sprintf('%s | %s', v1, v2); end
    end
end

function sShort = shortenLabel(sFull, maxLen)
    s = char(sFull); s = strtrim(s);
    if numel(s) <= maxLen, sShort = s;
    else sShort = [s(1:maxLen-3) '...']; end
end