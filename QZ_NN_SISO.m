% This script assumes these variables are defined:
%
%   TrainDataSet_SISO - input time series.
%   TrainDataSet_SISO - feedback time series.

data = TrainDataSet_SISO_Envelope;
data = [data{:}]
figure
plot(data)
xlabel("Symbols")
ylabel("Envelop")
% title("Time series channel envelop changes")

%partition dataset for training and test
numTimeStepsTrain = floor(0.9*numel(data));
dataTrain = data(1:numTimeStepsTrain+1);
dataTest = data(numTimeStepsTrain+1:end);
dataTrainStandardized = dataTrain;

%prepare predictors and responses
XTrain = dataTrainStandardized(1:end-1);
YTrain = dataTrainStandardized(2:end);

%define LSTM Network Architecture
numFeatures = 1;%1
numResponses = 1;%1
numHiddenUnits= 100;%200


%Create an LSTM regression network. Specify the LSTM layer to have 200 hidden units.
layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits)
    fullyConnectedLayer(numResponses)
    regressionLayer];
%Specify the training options. Set the solver to 'adam' and train for 250 
%epochs. To prevent the gradients from exploding, set the gradient threshold 
%to 1. Specify the initial learn rate 0.005, and drop the learn rate after 
%125 epochs by multiplying by a factor of 0.2.
options = trainingOptions('adam', ...
    'MaxEpochs',500, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.01, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',200, ...
    'LearnRateDropFactor',0.5, ...
    'Verbose',0, ...
    'Plots','training-progress');

%Train the LSTM network with the specified training options by using trainNetwork.
net = trainNetwork(XTrain,YTrain,layers,options);

%To forecast the values of multiple time steps in the future, use the 
%predictAndUpdateState function to predict time steps one at a time and 
%update the network state at each prediction. For each prediction, use the
%previous prediction as input to the function.
%Standardize the test data using the same parameters as the training data.

XTest = dataTest(1:end-1);
%To initialize the network state, first predict on the training data XTrain. 
%Next, make the first prediction using the last time step of the training 
%response YTrain(end). Loop over the remaining predictions and input the 
%previous prediction to predictAndUpdateState.
net = predictAndUpdateState(net,XTrain);
[net,YPred] = predictAndUpdateState(net,YTrain(end));

numTimeStepsTest = numel(XTest);
for i = 2:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net,YPred(:,i-1),'ExecutionEnvironment','cpu');
end

%The training progress plot reports the root-mean-square error (RMSE)
%calculated from the standardized data. Calculate the RMSE from the 
%unstandardized predictions.
YTest = dataTest(2:end);
%rmse = sqrt(mean((YPred-YTest).^2))

%Plot the training time series with the forecasted values.

figure
plot(dataTrain(1:end-1))
hold on
idx = numTimeStepsTrain:(numTimeStepsTrain+numTimeStepsTest);
plot(idx,[data(numTimeStepsTrain) YPred],'.-')
hold off
xlabel("Samples")
ylabel("Envelop")
%title("Prediction")
legend(["Actural" "Predict"])

%Compare the forecasted values with the test data.
figure
%subplot(2,1,1)
plot(YTest)
hold on
plot(YPred,'.-')
hold off
legend(["Actual" "Predict"])
ylabel("Envelop")
title("Prediction")

% subplot(2,1,2)
% stem(YPred - YTest)
% xlabel("Symbols")
% ylabel("Envelop")
% title("RMSE = " + rmse)

%Update Network State with Observed Values
%First, initialize the network state. To make predictions on a new sequence, 
%reset the network state using resetState. Resetting the network state 
%prevents previous predictions from affecting the predictions on the new data. 
%Reset the network state, and then initialize the network state by predicting 
%on the training data
net = resetState(net);
net = predictAndUpdateState(net,XTrain);
%Predict on each time step. For each prediction, predict the next time step 
%using the observed value of the previous time step. Set the 'ExecutionEnvironment' 
%option of predictAndUpdateState to 'cpu'.
YPred = [];
numTimeStepsTest = numel(XTest);
for i = 1:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net,XTest(:,i),'ExecutionEnvironment','cpu');
end

rmse = sqrt(mean((YPred-YTest).^2))


figure
%subplot(2,1,1)
plot(YTest)
hold on
plot(YPred,'.-')
hold off
legend(["Actual" "Predict"])
ylabel("Envelop")
xlabel("Samples")
% subplot(2,1,2)
% stem(YPred - YTest)
% xlabel("Symbols")
% ylabel("Error")
% title("RMSE = " + rmse)


