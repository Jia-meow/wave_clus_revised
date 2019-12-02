function fileName = loader(varargin)

directory = varargin{1};
if nargin > 1
    extracommand = varargin{2};
end

if ~iscell(directory)

d = pwd;

ind = strfind(d,'Analysis');
base = d(1:ind+7);
flag = 0;
if isempty(ind)
    ind = strfind(d,'UIC Data');
    base = d(1:ind+7);
    flag = 1;
end

if isempty(ind)
    ind = strfind(d,'uic data');
    base = d(1:ind+7);
    flag = 1;
end


if isempty(ind)
    ind = strfind(d,'CloudStation');
    base = d(1:ind+11);
    flag = 1;
end

[ret, name] = system('hostname');
if strcmp(name(1:end-1),'SAGGERELAB20605')
    base = 'E:\CloudStation';
    flag = 1;
    if strncmp(directory,'Labv',4)|strncmp(directory,'Code',4)
        base = 'C:\Users\Retina\CloudStation';
    end
end

if strncmp(name,'DESKTOP-VG3P',12)
    base = 'T:\UIC Data';
    flag = 1;
end



if nargin >1
   if varargin{2} == 1
        fileName = sprintf('%s/%s',base,directory);
   elseif varargin{2} == 2
        d = pwd;
        ind = strfind(d,'Dropbox');
      
        base = d(1:ind+6);
        fileName = sprintf('%s/%s',base,directory);
   end
       
else
    if flag == 0
        fileName = sprintf('%s/Data/%s',base,directory);
    else
        fileName = sprintf('%s/MEA/%s',base,directory);
    end
              
end

else
    
    for ii = 1:length(directory)
        if nargin == 1
            fileName(ii) = cellstr(loader(char(directory(ii))));
        else
            fileName(ii) = cellstr(loader(char(directory(ii)),extracommand));
        end
    end
end



