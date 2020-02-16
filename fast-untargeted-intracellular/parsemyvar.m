function [ composite ] = parsemyvar( myvar )
%PARSEMYVAR parsing the myvar output for the yeast conditions
%   Detailed explanation goes here
    labels = myvar.datasetLabels;

    myhand = @(x)strsplit(x,'/');

    output = cellfun(myhand, labels, 'UniformOutput', 0);

    %loop through and extract info. is there a better way to do this?

    composite(1).condition = output{1}{1};
    composite(1).data = [];
    composite(1).ODs = [];

    for i=1:length(output)
        current_cond = output{i}{1};
        result = find(strcmp({composite.condition},current_cond));

        if(length(result))
            composite(result).data = [composite(result).data; myvar.data(i,:)];
            composite(result).ODs = [composite(result).ODs; str2num(output{i}{2})];
        else
            composite(end+1).condition = current_cond;
            composite(end).data = myvar.data(i,:);
            composite(end).ODs = [str2num(output{i}{2})];
        end
    end

    [dontcare indices] = sort({composite.condition});

    %put the 30s condition in the front
    indices = [indices(end) indices(1:end-1)];

    composite = composite(indices);

end

