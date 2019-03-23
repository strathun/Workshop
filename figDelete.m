% figDelete
% This script allows the user to choose to delete any selected open figure

keepgoing = 1;

while keepgoing == 1
    prompt = 'Delete currently selected figure? Type Y or N or finished.';
    deleteYN = input(prompt,'s');
    if strcmp(deleteYN,'Y')
        h1=gca;
        titre = h1.Title.String
        prompt = 'Correct title? Y or N.';
        deleteYN = input(prompt,'s');
            if strcmp(deleteYN,'Y')
                delete(titre)
            end
    elseif strcmp(deleteYN,'finished')
        keepgoing = 0;
    end
end
    