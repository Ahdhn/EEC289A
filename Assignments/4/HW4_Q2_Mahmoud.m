%%%%%%%ECE 289A - An Introduction to Reinforcement Learning
%%%HW 4 - Q1 by Ahmed H. Mahmoud 

%%Regenerating figure 7.2

clc;
clear;
close all;
global  START_STATE END_STATE GAMMA currentStateValues
N_STATES =19; %all states
GAMMA =1; %discount
stateValues = zeros(1,N_STATES+2); %initial state values
states = 1:N_STATES; %all states but terminal states
START_STATE =11;%start from the middle 
END_STATE = [1 N_STATES+2]; 

realStateValues = (1/20.0).*[-20:2:21]; %true state value from bellman equation
realStateValues(1) = 0;
realStateValues(end) = 0;
truncateValue = 0.55; %truncate value for better display
steps = 2.^[0:9]; %all possible steps 
alphas = 0:0.1:1.0; 
episodes =10; %each run has 10 episodes 
runs =100;%performe 100 indepenent runs 
errors = zeros(length(steps), length(alphas)); %track the errors for each (step, alpha) combination

for run = 1:runs 
    for stepInd = 1:length(steps)
        step = steps(stepInd);
        for alphaInd = 1:length(alphas)
            alpha = alphas(alphaInd);
            str = ['run: ', num2str(run), ' step: ', num2str(step), ' alpha: ', num2str(alpha)];
            disp(str);
            currentStateValues = stateValues;
            for ep = 1:episodes
                temporalDifference(step, alpha);
                errors(stepInd, alphaInd) = errors(stepInd, alphaInd)+...
                    sqrt((1/N_STATES)*sum([currentStateValues-realStateValues].^2));                
            end
        end
    end 
end

errors = (1/(episodes*runs)).*errors ; %averaging 
errors(errors>truncateValue) = truncateValue;

plot(alphas,errors(1,:));
hold on
[err blah] = size(errors);
for n=2:err
    plot(alphas, errors(n,:));
end

legend(['n= ' num2str(steps(1))], ['n= ' num2str(steps(2))], ['n= ' num2str(steps(3))],...
       ['n= ' num2str(steps(4))],['n= ' num2str(steps(5))],['n= ' num2str(steps(6))],...
       ['n= ' num2str(steps(7))],['n= ' num2str(steps(8))],['n= ' num2str(steps(9))],...
       ['n= ' num2str(steps(10))]);
xlabel('\alpha','FontSize', 13,'FontWeight','bold');
ylabel('Avergae RMS error over 19 states and first 10 eps','FontSize', 13,'FontWeight','bold');
