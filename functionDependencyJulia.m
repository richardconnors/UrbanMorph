
% list of julia files
FF = dir('Q:\Dropbox\_LuxWork\MEVERST\code\Flexbus\*.jl');

allFunctions = strings();
for i=1:length(FF)
  buffer  = fileread(['Q:\Dropbox\_LuxWork\MEVERST\code\Flexbus\',FF(i).name]);
  pattern_tok = 'function\s(\w+)' ;
  tokens = regexp(buffer, pattern_tok, 'tokens');
  FF(i).subfunctions = string(tokens');
  allFunctions = [allFunctions; string(tokens')];


  pattern_tok = '\s*include\(\"(\w+)' ;
  tokens = regexp(buffer, pattern_tok, 'tokens');
  FF(i).include = string(tokens');
end



