function [bestMethod, index] = detBestMethod(finalTimeStep, finalSpaceStep)

%Function to conduct a simple analysis to determine the best method

%Input Arguments:
%finalTimeStep - stores the final time step found from the
%investMethodsTime function
%finalSpaceStep - stores the final space step found from the
%investMethodsSpace function

%Return Arguments:
%bestMethod - stores the best method of the four numerical methods
%index - stores the index value of the four methods used throughout the
%main code

    %Traverses through each of the four methods determining the difference
    %between the time and space steps, as time step is to be maximised and
    %space step is minimised ideally
    for i = 1:4
        finalMethod(i) = finalTimeStep(i) - finalSpaceStep(i);
    end

    %Determines which method has the highest final value 
    [value, index] = max(finalMethod);

    switch index
        case 1
            bestMethod = 'forward';
        case 2
            bestMethod = 'backward';
        case 3  
            bestMethod = 'dufort-frankel';
        case 4
            bestMethod = 'crank-nicholson';
    end
end