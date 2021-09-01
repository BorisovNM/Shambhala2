% This file presents example code lines on how to use the algorithm Shambhala2.
% It is therefore a guided example on how to use Shambhala2, not a script meant to
%   be ran.
%

    inputFileName = 'InputP0.txt'
    inData=readExpressionData(inputFileName,'log2');   
    Exp = inData.Samples;
    SYMBOL = inData.GeneList;
    SN = inData.SamplesName;
    NS = length(SN);
    NP = 39;
    % NP is the number of samples in the P-dataset
    NH = 222;
    % NH is the number of samples to be harmonized
    
    for ( i = 1:NH ) 
        
        disp(i);
        
        i0 = i;
        
        for ( jjj = 1:NP )             
            i0 = [i0 (jjj+NH)];
        end   
        
        EXP = Exp(:,i0);
        
        EXP = quantilenorm(EXP);
        
        dataN = CuBlock(EXP,[],4);

        log2e = log2(exp(1));

        DataN = dataN/log2e;

        EXPN = exp(DataN)-1;
        
        vecN = EXPN(:,1);
        
        if ( i == 1 ) 
            OUT = vecN;
        else    
            OUT = [OUT vecN];
        end
    end        
    
    outFileName = 'OutputP0.txt'; 
    outFile=fopen(outFileName,'w');
    nG=size(OUT,1);
    nS=size(OUT,2)
    for i=1:nG
        fprintf(outFile,'%s',SYMBOL{i,1});        
        for j=1:nS
            fprintf(outFile,' %f',OUT(i,j));        
        end
        fprintf(outFile, '\n');
    end
    fclose(outFile);
    
    
