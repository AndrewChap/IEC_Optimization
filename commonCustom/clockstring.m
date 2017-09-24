function [ out ] = clockstring(~)
% Drew Chap 9-12-2014
% gives a date and time string good for uniquely naming figures and videos

aclock = clock;
out = [num2str(aclock(1)),'-',num2str(aclock(2),'%02d'),'-',num2str(aclock(3),'%02d'),'_',...
       num2str(aclock(4),'%02d'),'-',num2str(aclock(5),'%02d'),'-',num2str(round(aclock(6)),'%02d')];

end

