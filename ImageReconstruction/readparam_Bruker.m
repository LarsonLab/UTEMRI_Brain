function out=readparam(filename)

    fid = 0;
    fid = fopen(filename, 'r'); 
    if(fid < 0)

        error('couldn''t open file!');
        return;
    end
    
    rem = fread(fid, inf, 'char'); 
    rem=char(rem);
    out = struct;
    
    while 1
        [str, rem] = strtokb(rem, '##');
        if isempty(str), break; end;
        
        str=char(str');
        
        %fprintf(1, '{%s}\n', strtrim(char(str(1:end-1))));
        
        out = addparameter(out, str);
        
        
    end
    
    fclose(fid);

    
    
function out = addparameter(infields, str)

    out = infields;
    [fieldname, strrem] = strtok(str, '$=');
    
    specialfields={'RecoStageEdges','RecoStageNode','RecoStagePasses'};
    for i=1:numel(specialfields)
        if contains(str,specialfields{i})
           
           try 
               val=splitbrack(strrem);
                out=setfield(infields,fieldname,val); 
           catch
           end
            return; 
        end
    end

    if isempty(fieldname) return; end;
   


    strrem = strtrim(strtok(strrem, '=')); % in case of leading = 
    
    if strfind(strrem,'$$')
        strrem=strtrim(strrem(1:strfind(strrem,'$$')-1));
    end
    
    if isempty(strrem) return; end;
    
    hascomp=regexp(strrem,'\@[0-9]+\*(.*)');
    if ~isempty(hascomp)
        strrem=uncomp(strrem);
    end


    % parse strrem to decide what is going on 
    if(strrem(1) == '(' && strrem(end) == ')')
        % it is a structure thing, put each thing separated by commas in a
        % different cell

        %unless, it's an expandable array - compression repeats many things
        %with @X*(Y). this has no commas and @s in it 
        if and(sum(strrem=='@')>0, sum(strrem==',')==0)
            % guess compressed list. first thing in brackets is length of
            % list

            [thinginbrackets,strrem] = strtok(strrem, '()');
            contents=zeros([str2double(thinginbrackets) 1]);cum=0; 
            mystart=find(strrem=='@',1);
            if ~isempty(mystart)
                strrem=strrem(mystart+1:end);

                while numel(strrem)>1
                    [thisthing,strrem]=strtok(strrem,'@');

                    if any(thisthing=='*')
                        [a,b]=strtok(thisthing,'*');
                        a=str2double(a); 
                        b=str2double(b(~ismember(b,'*() ')));
                        
                        contents(cum+(1:a))=b;
                        cum=cum+a;
                    end
                end
            end
           
        else

            
            %remove the brackets
            strrem = strrem(2:end-1);
            
            cell = {};
            i = 1;
            strrem2 = strtrim(strrem);
            while 1 
                [item, strrem2] = strtok(strrem2);
                if isempty(item)
                    if isempty(strrem2) break; end
                    continue;
                end;
                
                
                
                if ~or(isempty(str2double(item)), isnan(str2double(item)))
                    cell{i} = str2double(item);
                else
                    if strcmp(item,'On')
                        item = 1;
                    elseif strcmp(item,'Off')
                        item = 0;
                    end
                    
                    v = str2double(item); 
                    if(~or(isempty(v),isnan(v)))
                        item = str2double(item);
                    end
                    
                    if(item(end) == ',') item=item(1:end-1); end;
                    cell{i} = item;
                end
                i=i+1;
            end
            
            contents = cell;
        end
        
    elseif(strrem(1) == '(')
        % this number in brackets will define the space of what is to come
        
        % get the thing in brackets
        [thinginbrackets,strrem] = strtok(strrem, '()');
        
        % this line may cause trouble?
        strrem = strrem(3:end); % delete cr + )
        
        thinginbrackets = strtrim(thinginbrackets);
        
        % are there any commas in it?
        if strfind(thinginbrackets,',')
            % yes, there are multiple dims
            i = 1;
            while 1 
                [dd, thinginbrackets] = strtok(thinginbrackets, ',');
                if(isempty(dd)) break; end;
                dim(i) = str2double(dd);
                i = i+1; 
                
            end
        else
            % no, just a number
            dim = [str2double(thinginbrackets) 1];
            
        end
        
        dim = [fliplr(dim) 1 1 1];
        contents = zeros(dim); 
        i=1; 
        if numel(strrem)<50
        while 1
            [it,strrem] = strtok(strrem);
            
            if isempty(it)
                break;
            end;
            
            if isempty(str2num(it))
                contents = strcat(it,strrem);
                break;
            end
            
            contents(i) = str2double(it);
            i=i+1;
        end
        else
            
            %hope it's floats? 
            contents=textscan(strrem,'%f');
            contents=contents{1}; 
            
        end
        
        contents = permute(contents, [2 1 3:ndims(contents)]);
    else
        contents = strrem;
        
        if strcmp(strrem,'On')
            contents = 1;

        end 
        if strcmp(strrem,'Off')
            contents = 0;

        end;
        
        v = str2double(strrem); 
        if isnan(v) v = []; end;
        if ~isempty(v)
            contents = v;
        end

    end
        
            
    
    out = setfield(infields, fieldname, contents);
    
    
    
    
    function [out,rem] = strtokb(in, del)
        
        if ~ischar(in) 
            in=char(in);
        end
        
        if(size(in,1)==1) in==in'; end;
        
        in = strtrim(in); 
        n  = strfind(in',del); 
        
        if numel(n) == 0
            out = in;
            rem = []; 
            return;
        end
        
        if n(1)==1
            [out,rem]=strtokb(in((numel(del)+1):end),del);
            return;
        end
        
        out=in(1:(n(1)-1));
        rem=in((n(1)+numel(del)):end);
        
            
function out=splitbrack(in)
    %splits a string into cells of self-contained brackets
    opens=find(in=='(');
    closes=find(in==')');

    %find the nesting level of brackets 
    nest=zeros(size(opens));
    for i=1:numel(nest)
        nest(i)=sum(opens<opens(i)) - sum(closes<opens(i));
    end

    %kill the next close after each nested bracket
    todo=find(nest>0);kill=zeros(size(todo));
    for i=1:numel(todo)
        kill(i)=closes(find(closes>opens(todo(i)),1));
        closes=closes(closes~=kill(i));
    end

    %now create a cell array for each of the bracketed portions
    opens=opens(nest==0); 
    for i=1:numel(opens)
        myitem=in(opens(i)+1:closes(i)-1);
        myitem=strtrim(myitem);
        myitem=myitem(myitem~=10); %kill new lines
        out{i}=myitem;
    end



function out=uncomp(in)
    cs=regexp(in,'\@[0-9]+\*(.*)');

    atend=cs+find(in(cs:end)==')',1)-1;

    cinst=in(cs+1:atend);
    num=str2double(cinst(1:find(cinst=='*',1)-1));
    what=[cinst(1+find(cinst=='(',1):find(cinst==')',1',"last")-1) ' '];
    rep=repmat(what,[1 num]);
    rep=rep(1:end-1); %kill last space

    out=[in(1:cs-1) rep in(atend+1:end)];

