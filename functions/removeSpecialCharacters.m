function cleanedText = removeSpecialCharacters(inputText)
    % This function removes special characters from the input string
    % and keeps only alphanumeric characters and spaces.

    % Define a pattern to match special characters
    pattern = '[^a-zA-Z0-9 ]';

    % Replace special characters with an empty string
    cleanedText = regexprep(inputText, pattern, '');
end