function struct3 = appendStruct(struct1,struct2,flagger)
% Function for combining the fields from struct1 and struct2 into a single
% structure struct3

    if nargin == 2
        flagger = false;
    end
    struct3 = [];
    if ~isempty(struct1)
        fields = fieldnames(struct1);
        if flagger
            for ii = 1:length(fields)
                [m1(ii),n1(ii)] = size(getfield(struct1,char(fields(ii))));
                [m2(ii),n2(ii)] = size(getfield(struct2,char(fields(ii))));
            end
            mM1 = mode(m1);
            mN1 = mode(n1);
            if mM1<mN1
                for ii = 1:length(fields)
                    if mN1==n1(ii)
                        temp = getfield(struct1,char(fields(ii)));
                        struct1 = setfield(struct1,char(fields(ii)),temp');
                    end
                end
            end

            mM2 = mode(m2);
            mN2 = mode(n2);
            if mM2<mN2
                for ii = 1:length(fields)
                    if mN2==n2(ii)
                        temp = getfield(struct2,char(fields(ii)));
                        struct2 = setfield(struct2,char(fields(ii)),temp');
                    end
                end
            end
        end
        
        for ii = 1:length(fields)
            FOI1 = getfield(struct1,char(fields(ii)));
            FOI2 = getfield(struct2,char(fields(ii)));
            if flagger
                struct3 = setfield(struct3,char(fields(ii)),[FOI1; FOI2]);
            else
                struct3 = setfield(struct3,char(fields(ii)),[FOI1(:); FOI2(:)]);
            end
        end
    else
        struct3 = struct2;
    end

end
