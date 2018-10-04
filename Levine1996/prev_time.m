% This function calcultates the number of time steps of state "value" in the history,
% starting from the last state from the history, and count back to when an alternative stae appears
function num=prev_time(history,value)
num=0;
history(history<0)=[]; % delete points without cells
for i=length(history):-1:1
    if history(i)==value
        num=num+1;
    else
        break;
    end
end

end
