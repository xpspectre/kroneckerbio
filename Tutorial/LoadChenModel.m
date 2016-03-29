%% Chen model
opts.Verbose = 1;
opts.Order = 2;

%% Load SBML
mSimbio = sbmlimport('Chen2009_ErbB_A431.xml');

% Purge notes from parameters
for i = 1:numel(mSimbio.Parameters)
    mSimbio.Parameters(i).Notes = '';
end

clear i

%% Purge XML in notes
for i = 1:numel(mSimbio.Species)
    notes = mSimbio.Species(i).Notes;
    notes = regexp(notes, '(?<=<notes>\s*)[^<>\s]+', 'match', 'once');
    mSimbio.Species(i).Notes = strtrim(notes);
end

clear i notes

%% Specify outputs
% Outputs: 
Ynames = {
            'ErbB1';
            'ErbB1#P';
            'ErbB2';
            'ErbB2#P';
            'ErbB3';
            'ErbB3#P';
            'ErbB4';
            'ErbB4#P';
            'GAP';
            'Grb2';
            'Gab1';
            'Gab1#P';
            'Gab1#P#P;'
            'Sos';
            'Sos#P';
            'Ras:GTP';
            'Ras:GDP';
            'Shc';
            'Shc#P';
            'Raf';
            'Raf#P';
            'MEK';
            'MEK#P';
            'MEK#P#P';
            'ERK';
            'ERK#P';
            'ERK#P#P';
            'Pase1';
            'Pase2';
            'Pase3';
            'PI3K';
            'AKT';
            'AKT#P';
            'AKT#P#P';
            'PDK1';
            'Pase4';
            'RTK_Pase';
            'Shp';
            'Shp2';
            'PTEN';
            'Ser';
            'Pase9t';
            'cPP';
            'PIP2';
            'PIP3';
            'ATP';
         };
     
Ymembers = {
            cat(1, ... % ErbB1
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^ErbB1$'),'Name'), ... % unique
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^ErbB1:'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^EGF:ErbB1:'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^\(EGF:ErbB1:ErbB.\)$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^\(EGF:ErbB1:ErbB.\)[^#]'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '::EGF:ErbB1:Inh'),'Name'), ... % unique
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '::EGF:ErbB1:Inh'),'Name'), ... % repeat
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^[^2]*ErbB1\)$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^[^2]*ErbB1\)[^#]'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '2\(EGF:ErbB1:'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '2\(EGF:ErbB1:'),'Name'), ... % repeat
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '2\(EGF:ErbB1\)$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '2\(EGF:ErbB1\)$'),'Name'), ... % repeat
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '2\(EGF:ErbB1\)[^#]'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '2\(EGF:ErbB1\)[^#]'),'Name'), ... % repeat
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^ErbB1_h$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^ErbB1_h:'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^EGF:ErbB1_h:'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '::EGF:ErbB1_h'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '::EGF:ErbB1_h'),'Name'), ... % repeat
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '2\(EGF:ErbB1_h:'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '2\(EGF:ErbB1_h:'),'Name') ... % repeat
            );
            cat(1, ... % ErbB1#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^EGF:ErbB1#P$'),'Name'), ... % unique
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^.*ErbB1.*\)#P'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^2.*ErbB1.*\)#P'),'Name') ... % repeat when there is a 2
            );
            cat(1, ... % ErbB2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB2$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB2\)$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB2\)[^#]'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB2:ATP'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB2:Inh'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB2:ErbB2:Inh'),'Name'), ... % repeat for this species
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^2.*ErbB2$'),'Name'), ... % repeat with 2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^2.*ErbB2\)$'),'Name'), ... % repeat with 2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^2.*ErbB2\)[^#]'),'Name'), ... % repeat with 2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB2:ErbB.$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB2:ErbB.\)$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB2:ErbB.\)[^#]'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB2:ErbB2#P'),'Name') ... % has one clean one
            );
            cat(1, ... % ErbB2#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^ErbB2#P$'),'Name'), ... % unique
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^.*ErbB2.*\)#P'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^2.*ErbB2.*\)#P'),'Name'), ... % repeat when there is a 2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB2:ErbB2#P'),'Name') ... % has one phospho one
            );
            cat(1, ... % ErbB3
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB3$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB3[^\)]*\)$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB3[^\)]*\)[^#]'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^2.*ErbB3$'),'Name'), ... % repeat with 2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^2.*ErbB3\)$'),'Name'), ... % repeat with 2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^2.*ErbB3\)[^#]'),'Name'), ... % repeat with 2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB3:ATP'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB3:Inh'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB3:ErbB.:Inh'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB3:ErbB.$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '\(ErbB3:ErbB.\)$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '\(ErbB3:ErbB.\)[^#]'),'Name') ...
            );
            cat(1, ... % ErbB2#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^ErbB3#P$'),'Name'), ... % unique
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^.*ErbB3.*\)#P'),'Name') ...
            );
            cat(1, ... % ErbB4
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB4$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB4[^\)]*\)$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB4[^\)]*\)[^#]'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^2.*ErbB4$'),'Name'), ... % repeat with 2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^2.*ErbB4\)$'),'Name'), ... % repeat with 2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^2.*ErbB4\)[^#]'),'Name'), ... % repeat with 2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB4:ATP'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB4:Inh'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB4:ErbB.:Inh'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ErbB4:ErbB.$'),'Name') ...
            );
            cat(1, ... % ErbB4#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^ErbB4#P$'),'Name'), ... % unique
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', '^.*ErbB4.*\)#P'),'Name') ...
            );
            cat(1, ... % GAP
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'GAP'),'Name') ...
            );
            cat(1, ... % Grb2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Grb2'),'Name') ...
            );
            cat(1, ... % Gab1
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Gab1$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Gab1[^#]'),'Name') ...
            );
            cat(1, ... % Gab1#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Gab1#P$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Gab1#P[^#]'),'Name') ...
            );
            cat(1, ... % Gab1#P#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Gab1#P#P'),'Name') ...
            );
            cat(1, ... % Sos
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Sos$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Sos[^#]'),'Name') ...
            );
            cat(1, ... % Sos#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Sos#P'),'Name') ...
            );
            cat(1, ... % Ras:GTP
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Ras:GTP'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Ras_activated:GTP'),'Name') ...
            );
            cat(1, ... % Ras:GDP
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Ras:GDP'),'Name') ...
            );
            cat(1, ... % Shc
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Shc$'),'Name') ...
            );
            cat(1, ... % Shc#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Shc#P'),'Name') ...
            );
            cat(1, ... % Raf
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Raf$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Raf[^#]'),'Name') ...
            );
            cat(1, ... % Raf#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Raf#P'),'Name') ...
            );
            cat(1, ... % MEK
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'MEK$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'MEK[^#]'),'Name') ...
            );
            cat(1, ... % MEK#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'MEK#P$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'MEK#P[^#]'),'Name') ...
            );
            cat(1, ... % MEK#P#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'MEK#P#P'),'Name') ...
            );
            cat(1, ... % ERK
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ERK$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ERK[^#]'),'Name') ...
            );
            cat(1, ... % ERK#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ERK#P$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ERK#P[^#]'),'Name') ...
            );
            cat(1, ... % ERK#P#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ERK#P#P'),'Name') ...
            );
            cat(1, ... % Pase1
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Pase1'),'Name') ...
            );
            cat(1, ... % Pase2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Pase2'),'Name') ...
            );
            cat(1, ... % Pase3
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Pase3'),'Name') ...
            );
            cat(1, ... % PI3K
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'PI3K'),'Name') ...
            );
            cat(1, ... % AKT
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'AKT$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'AKT:PDK1'),'Name') ... % unique
            );
            cat(1, ... % AKT#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'AKT#P$'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'AKT#P[^#]'),'Name') ...
            );
            cat(1, ... % AKT#P#P
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'AKT#P#P'),'Name') ...
            );
            cat(1, ... % PDK1
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'PDK1'),'Name') ...
            );
            cat(1, ... % Pase4
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Pase4'),'Name') ...
            );
            cat(1, ... % RTK_Pase
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'RTK_Pase'),'Name') ...
            );
            cat(1, ... % Shp
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Shp$'),'Name') ...
            );
            cat(1, ... % Shp2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Shp2'),'Name') ...
            );
            cat(1, ... % PTEN
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'PTEN'),'Name') ...
            );
            cat(1, ... % Ser
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Ser'),'Name') ...
            );
            cat(1, ... % Pase9t
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'Pase9t'),'Name') ...
            );
            cat(1, ... % cPP
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'cPP'),'Name') ...
            );
            cat(1, ... % PIP2
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'PIP2'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'PIP2)[23456]'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'PIP2)[3456]'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'PIP2)[456]'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'PIP2)[56]'),'Name'), ...
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'PIP2)[6]'),'Name') ... % copied a correct number of times
            );
            cat(1, ... % PIP3
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'PIP3'),'Name') ...
            );
            cat(1, ... % ATP
                dealfield2cell(mSimbio.sbioselect('WHERE', 'Notes', 'regexp', 'ATP'),'Name') ...
            );
           };
       
%% Create Kronecker Bio mass action model
m = LoadModelSbmlMassAction(mSimbio, opts);

%% Clean up the model
% Add Outputs
for iy = 1:numel(Ymembers)
    [exprs, values] = count_unique(Ymembers{iy});
    exprs = strcat('\.', exprs, '$'); % Exact matches on species
    m = AddOutput(m, Ynames{iy}, exprs, values);
end

% Remove unecessary parameters
nullParameters = false(m.nk,1);
for ik = 1:m.nk
    if m.Parameters(ik).Value == 0
        nullParameters(ik) = true;
        rInd = strcmp(m.Parameters(ik).Name, {m.Reactions.Parameter});
        m.Reactions(rInd) = [];
    end
end

% Unused paramters
unusedParameters = false(m.nk,1);
for ik = 1:m.nk
    found = find(strcmp({m.Reactions.Parameter}, m.Parameters(ik).Name), 1);
    if isempty(found)
        unusedParameters(ik) = true;
    end
end

m.Parameters(nullParameters | unusedParameters) = [];

% Rearrange inputs
m.Inputs = m.Inputs([1;4;3;2]); % EGF, HRG, Inh, endosomal HRG

m = FinalizeModel(m);

clear iy exprs values ik nullParameters unusedParameters rInd found ik temp
