function mrprot = parse_xprot(rawarr, parse_arrays)
% parse XProtocol structure stored in a char array
%  E. Auerbach, CMRR, 2013
%    some changes by Sebastian Stenmark, 2012 (SS)
%    some fixes by Steven Baete, NYU LMC CBI, 2013 (SB)

mrprot = [];
tagName = cell(1,10);
nest = -1;

% first, make sure it is a char array
if ~ischar(rawarr), error('parse_xprot: input is not a char array!'); end

tstart = tic;

% find <XProtocol>{
xstart = strfind(rawarr,'<XProtocol>');

if (xstart)
    workarr = extractBraceString(rawarr(xstart(1):end)); % <XProtocol> {}
    if (parse_arrays), nest = 0; end
    mrprot = parse_loop(mrprot, workarr, tagName, 0, nest, []);
else
    error('parse_xprot() failed to find <XProtocol>!');
end

telapsed = toc(tstart);
fprintf('   ** parse_xprot completed in %.3f s\n', telapsed);

%--------------------------------------------------------------------------

function mrprot = parse_loop(mrprot, workarr, tagName, level, nest_level, array_data)
% this function is called recursively to parse the xprotocol

if (level > 10), fprintf('exceeded preallocated tag depth; increase size of tagName!\n'); end

% eja debug
%if (nest_level > 0)
%    fprintf('array level %d: process this:\n', nest_level);
%    disp(workarr);
%    fprintf('with these received values:\n');
%    disp(array_data);
%end

[tagStr, tagType, workarr] = findNextTag(workarr);

while (~isempty(workarr))
    % for parammap, add the name as another level of the current array name
    % and spawn off another copy of this function to deal with them
    if (strncmpi(tagType, 'ParamMap', 8))
        if (~isempty(tagStr))
            level = level + 1;
            tagName{level} = make_safe_fieldname(tagStr);
        end
        
        [stubArr, pos] = extractBraceString(workarr);
        workarr = workarr(pos+1:end);
        
        array_data_stub = [];
        if (nest_level > 0)
            [array_data_stub, pos] = extractBraceString(array_data);
            array_data = array_data(pos+1:end);
        end
        
        mrprot = parse_loop(mrprot, stubArr, tagName, level, nest_level, array_data_stub);
        
        if (~isempty(tagStr))
            level = level - 1;
        end

    % for paramarray, special processing is required
    % can just skip them since the most important ones are still in the text mrprot
    % but sometimes we don't have the text mrport and need to process them
    elseif (strncmpi(tagType, 'ParamArray', 10))
        [localarr, pos] = extractBraceString(workarr);
        workarr = workarr(pos+1:end);
        
        %debug = false;
        
        % only process if option specified
        if (nest_level >= 0)
            nest_level = nest_level + 1;
        
            %fprintf('found ParamArray %s\n', tagStr);
            if (isempty(tagStr))
                if (nest_level > 1)
                    localtagName = tagName;
                else
                    error('found empty ParamArray at top level');
                end
            else
                localtagName = cell(1,10);
                localtagName{1} = make_safe_fieldname(tagStr);
            end
            
            % skip the datatype block and get the values
            [~, pos] = extractBraceString(localarr);
            datadefarr = localarr(1:pos);
            
            if (nest_level == 1)
                % this is the top-level array which contains all of the values
                %fprintf('found array %s\n',localtagName{1});
                valuearr = localarr(pos+1:end);
            else
                % for nested array, get values passed to us by parent(s)
                %fprintf('found nested array %s (level %d)\n',localtagName{1},nest_level);
                %if (strcmp(localtagName{1},'aFFT_SCALE')), debug = true; end
                [valuearr, pos] = extractBraceString(array_data);
                array_data = array_data(pos+1:end);
                
                % in principle this shouldn't be necessary, but it solves some error I
                % have with ParamMap/ParamArray nesting order
                if (isempty(valuearr))
                    localarr = localarr(pos+1:end);
                    [~, pos] = extractBraceString(localarr);
                    valuearr = localarr(pos+1:end);
                end
                
                %if (debug), disp(array_data); end
            end
            
            %if (debug)
            %    disp(localarr);
            %    fprintf('process these data types:\n');
            %    disp(datadefarr);
            %    fprintf('with these values:\n');
            %    disp(valuearr);
            %end

            dim = 0;
            local_mrprot = cell(128,1);
            while (~isempty(valuearr))
                [stubArr, pos] = extractBraceString(valuearr);
                if (isempty(stubArr)), break; end
                valuearr = valuearr(pos+1:end);
                
                % in case there are e.g. <MinSize> tags at the start of the arrayed values, strip everything before the opening brace
                spos = strfind(stubArr,'{');
                if (spos), stubArr = stubArr(spos(1):end); end
            
                dim = dim + 1;
                local_mrprot{dim} = parse_loop([], datadefarr, localtagName, 0, nest_level, ['{ ' stubArr ' }']);
            end
            
            if (dim)
                if (dim > 128), fprintf('exceeded preallocated array depth (%d > 128); increase size of local_mrprot!\n', dim); end
                
                % truncate array if not all fields are filled (should always indicate "empty" values)
                if (dim > 1)
                    expected_fields = size(fieldnames(local_mrprot{1}),1);
                    for m=2:dim
                        if (isempty(local_mrprot{m}) || (size(fieldnames(local_mrprot{m}),1) < expected_fields))
                            dim = m - 1;
                            break;
                        end
                    end
                end
                
                localtagName = [tagName(1:level) localtagName(1)];
                mrprot = setfield(mrprot, localtagName{:}, cell2mat(local_mrprot(1:dim)));
            end
            
            %if (strcmp(localtagName{1},'iGlobalRFactor')), error('stop'); end
            %if (debug), error('stop'); end
        end

    % these are values that we know how to process, so do so
    elseif ((strncmpi(tagType, 'ParamBool', 9))     || ...
            (strncmpi(tagType, 'ParamLong', 9))     || ...
            (strncmpi(tagType, 'ParamDouble', 11))  || ...
            (strncmpi(tagType, 'ParamString', 11)))
        tagName{level+1} = make_safe_fieldname(tagStr);
        
        % special case for arrayed single parameters (not within parammap)
        if ((nest_level > 0) && isempty(tagStr))
            tagName{level+1} = 'placeholder';
        end

        [stubStr, pos] = extractBraceString(workarr);
        workarr = workarr(pos+1:end);
        
        %fprintf('%s (%s) = %s\n',tagName{level+1},tagType,stubStr);
        
        value = get_xprot_value(stubStr, tagType);
        if (nest_level > 0)
            [arrstubStr, pos] = extractBraceString(array_data);
            array_data = array_data(pos+1:end);
            
            value_from_arr = get_xprot_value(arrstubStr, tagType);
            
            % use value from array unless it is empty and we have a default
            if ((isempty(value)) || (~isempty(value_from_arr)))
                value = value_from_arr;
            end
        end
        
        fields = {tagName{1:level+1}, value};
        mrprot = setfield(mrprot, fields{:});
        
    % if we are processing an array, we may come across some
    % <Default> tags preceding the data type that we can ignore
    % known: ParamMap ParamString ParamArray; throw error for others
    elseif ((nest_level > 0)                            && ...
            (strncmpi(tagStr, 'Default', 7)))
        [~, tmpTagType, ~] = findNextTag(workarr);
        if ((~strncmpi(tmpTagType, 'ParamMap', 8))      && ...
            (~strncmpi(tmpTagType, 'ParamString', 11))  && ...
            (~strncmpi(tmpTagType, 'ParamLong', 9))     && ...
            (~strncmpi(tmpTagType, 'ParamArray', 10)))
            error('found <Default> tag for arrayed %s, unexpected!', tmpTagType);
        end
        
    % ignore these also
    % skip to end of line for these
    elseif ((strncmpi(tagStr, 'Name', 4))             || ...
            (strncmpi(tagStr, 'ID', 2))               || ...
            (strncmpi(tagStr, 'Comment', 7))          || ...
            (strncmpi(tagStr, 'Label', 5))            || ...
            (strncmpi(tagStr, 'Visible', 7))          || ...
            (strncmpi(tagStr, 'Userversion', 11)))
        lend = strfind(workarr, char(10));
        workarr = workarr(lend(1)+1:end);

    % if we are processing an array skip to end of line for these
    elseif ((nest_level > 0)                            && ...
            (strncmpi(tagStr, 'MinSize', 7))            || ...
            (strncmpi(tagStr, 'MaxSize', 7)))
        lend = strfind(workarr, char(10));
        workarr = workarr(lend(1)+1:end);

    % if we are processing an array, only the above parameters are
    % supported; stop with error if we see anything else
    elseif (nest_level > 0)
        error('unsupported arrayed tag type %s', tagStr); 

    % we don't care about the things below, but acknowledge that we know
    % about them and skip over their data
    elseif ((strncmpi(tagType, 'ParamChoice', 11))      || ...
            (strncmpi(tagType, 'Pipe', 4))              || ...
            (strncmpi(tagType, 'Connection', 10))       || ...
            (strncmpi(tagType, 'Event', 5))             || ...
            (strncmpi(tagType, 'ParamFunctor', 12))     || ...
            (strncmpi(tagType, 'Method', 6))            || ...
            (strncmpi(tagStr,  'EVAStringTable', 14))   || ...
            (strncmpi(tagType, 'ProtocolComposer', 16)) || ...
            (strncmpi(tagType, 'Dependency', 10))       || ...
            (strncmpi(tagType, 'ParamCardLayout', 15))  || ...
            (strncmpi(tagType, 'EVACardLayout', 13)))
        [~, pos] = extractBraceString(workarr);
        workarr = workarr(pos+1:end);
        %if (nest_level > 0)
        %    [~, pos] = extractBraceString(array_data);
        %    array_data = array_data(pos+1:end);
        %end

    % something unknown has happened if we reach this point, so throw a warning
    else
        fprintf('parse_xprot(): WARNING: found unknown tag %s (%s)\n',tagStr,tagType);
    end

    [tagStr, tagType, workarr] = findNextTag(workarr);
end

%--------------------------------------------------------------------------

function value = get_xprot_value(stubStr, tagType)
% stubStr contains the values of interest, but also might include
% modifiers. except for <Default>, we will just ignore the modifiers,
% which seem to always be terminated by CR.

[tagStr, newtagType, remStr] = findNextTag(stubStr);
while (~isempty(remStr)) % found a tag
    skip = false;
    
    if (~isempty(newtagType)) % not a modifier tag???
        error('parse_xprot()::get_xprot_value(): ERROR: found unknown tag!')
    else
        % these are the modifier tags we know about
        if ((strncmpi(tagStr, 'Precision', 9))      || ...
            (strncmpi(tagStr, 'LimitRange', 10))    || ...
            (strncmpi(tagStr, 'MinSize', 7))        || ...
            (strncmpi(tagStr, 'MaxSize', 7))        || ...
            (strncmpi(tagStr, 'Limit', 5))          || ...
            (strncmpi(tagStr, 'Default', 7))        || ...
            (strncmpi(tagStr, 'InFile', 6))         || ...
            (strncmpi(tagStr, 'Context', 7))        || ...
            (strncmpi(tagStr, 'Dll', 3))            || ...
            (strncmpi(tagStr, 'Class', 5))          || ...
            (strncmpi(tagStr, 'Comment', 7))        || ...
            (strncmpi(tagStr, 'Label', 5))          || ...
            (strncmpi(tagStr, 'Tooltip', 7))        || ...
            (strncmpi(tagStr, 'Visible', 7))        || ...
            (strncmpi(tagStr, 'Unit', 4)))
            % acknowledge that we know about these tags
        else
            fprintf('parse_xprot(): WARNING: found unknown modifier ''%s''\n',tagStr);
        end
        
        if (strncmpi(tagStr, 'Default', 7))
            % for <Default>, check whether an actual value is present.
            % if so, remove the line, otherwise just skip the tag
            %fprintf('found default tag\n%s',stubStr);
            
            testStr = strip_tag_line(stubStr, tagStr);
            if (isempty(deblank(testStr)))
                % no value found, use default
                stpos = strfind(stubStr, ['<' tagStr '>']);
                stubStr = stubStr(stpos(1)+length(tagStr)+2:end);
                %fprintf('using default\n%s\n', stubStr);
            else
                % found value, use it
                stubStr = testStr;
                %fprintf('using value\n%s\n', stubStr);
            end
        else
            % otherwise, remove the entire line containing the tag
            stubStr = strip_tag_line(stubStr, tagStr);
        end
    end

    if (~skip), [tagStr, newtagType, remStr] = findNextTag(stubStr); end
end

if (isempty(stubStr))
    value = [];
    ok = true;
else
    if (strncmpi(tagType, 'ParamString', 11))
        value = getQuotString(stubStr);
        ok = true;
    else
        stubStr = strrep(stubStr, char(10), ' '); % remove newlines
        if (strncmpi(tagType, 'ParamBool', 9))
            stubStr = strrep(stubStr, '"true"', '1');
            [value,ok] = str2num(stubStr); %#ok<ST2NM>
            value = (value ~= 0);
        elseif (strncmpi(tagType, 'ParamLong', 9))
            [value,ok] = str2num(stubStr); %#ok<ST2NM>
        elseif (strncmpi(tagType, 'ParamDouble', 11))
            [value,ok] = str2num(stubStr); %#ok<ST2NM>
        end
    end
end

if (~ok)
    fprintf('WARNING: get_xprot_value failed: %s\n', stubStr);
    disp(value);
    %error('stop')
end

%--------------------------------------------------------------------------

function stvar = make_safe_fieldname(tagStr)
% this function checks potential fieldnames and makes sure they are valid
% for MATLAB syntax, e.g. 2DInterpolation -> x2DInterpolation (must begin
% with a letter)

tagStr = strtrim(tagStr);

if isempty(tagStr)
    stvar = '';
else
    if (isletter(tagStr(1)))
        stvar = tagStr;
    else
        stvar = strcat('x', tagStr);
    end
    
    if strfind(stvar, ';'), stvar = strrep(stvar, ';', '_'); end
    if strfind(stvar, '@'), stvar = strrep(stvar, '@', '_'); end % VD13 (SS)
    if strfind(stvar, '-'), stvar = strrep(stvar, '-', '_'); end % (SB)
end

%--------------------------------------------------------------------------

function [stvar, lastpos] = extractBraceString(text)
% extracts string from within curly braces, ignoring nested braces

tlen = length(text);
tstart = strfind(text,'{');
lastpos = 0;

stvar = [];
if (tstart)
    tnests = 1;
    for x=tstart(1)+1:tlen
        if (text(x) == '{')
            tnests = tnests + 1;
        elseif (text(x) == '}')
            tnests = tnests - 1;
            lastpos = x;
        end
        if (tnests == 0)
            stvar = text(tstart(1)+1:x-1);
            break;
        end
    end
end

%--------------------------------------------------------------------------

function [tagStr, tagType, remStr] = findNextTag(inStr)
% returns <tag> name and the remainder of the string following the tag.
% for e.g. <Tag>, returns tagStr='Tag', tagType=''
% for e.g. <ParamLong."Tag">, returns tagStr='Tag', tagType='ParamLong'
% if no tag is found, returns null strings

tagStr = [];
tagType = [];
remStr = [];

startPos = strfind(inStr,'<'); % look for start of tag
if (startPos)
    endPos = strfind(inStr,'>'); % look for end of tag
    if (endPos)
        % found complete tag, but check to make sure the < and > are not
        % quoted/escaped (exception for e.g. ICEOut_DataRoleSeriesMap)
        if (startPos > 1)
            if ('''' == inStr(startPos-1)), return; end
        end
        
        % extract tag
        fullTag = inStr(startPos(1)+1:endPos(1)-1);
        
        % now check for name/type
        dotPos = strfind(fullTag,'."');
        if (dotPos)
            tagStr = getQuotString(fullTag(dotPos(1)+1:end));
            tagType = fullTag(1:dotPos(1)-1);
        else
            tagStr = fullTag;
        end
        
        % return remainder
        remStr = inStr(endPos(1)+1:end);
    end
end

%--------------------------------------------------------------------------

function stvar = getQuotString(text)
% extracts string between double quotes, e.g. "string"
%  also works with double-double quotes, e.g. ""string""

idx = strfind(text,'"');

if ( (length(idx) == 4) && (idx(1)+1 == idx(2)) && (idx(3)+1 == idx(4)) ) % double-double quotes
    stvar = text(idx(2)+1:idx(3)-1);
elseif (length(idx) >= 2) % double quotes, or ??? just extract between first and last quotes
    stvar = text(idx(1)+1:idx(end)-1);
else % malformed?
    stvar = text;
end

%--------------------------------------------------------------------------

% function stvar = delQuotString(text)
% % deletes string between double quotes, e.g. "string"
% %  also works with double-double quotes, e.g. ""string""
% 
% idx = strfind(text,'"');
% 
% if ( (length(idx) == 4) && (idx(1)+1 == idx(2)) && (idx(3)+1 == idx(4)) ) % double-double quotes
%     stvar = [text(1:idx(2)) text(idx(3):end)];
% elseif (length(idx) >= 2) % double quotes, or ??? just extract between first and last quotes
%     stvar = [text(1:idx(1)) text(idx(end):end)];
% else % malformed?
%     stvar = text;
% end

%--------------------------------------------------------------------------
function outlin = strip_tag_line(inlin, tag)
% find a <tag> and remove everything from the tag to eol

lstart = strfind(inlin, ['<' tag '>']);
lend = strfind(inlin(lstart(1)+1:end), char(10));
if (isempty(lend))
    outlin = '';
else
    if (lstart > 1)
        outlin = [inlin(1:lstart(1)-1) inlin(lend(1)+lstart(1)+1:end)];
    else
        outlin = inlin(lend(1)+1:end);
    end
end
