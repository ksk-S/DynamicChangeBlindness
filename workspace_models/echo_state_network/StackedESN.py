"""
    Stacked ESN model.
    Replace "feedback of output" by "feedback of output from upper network"
"""

import numpy as np
import numpy.random as rnd
from FeedbackBaseESN import BaseESN
from FeedbackPredictionESN import PredictionESN

from easyesn import backend as B

from sklearn.linear_model import Ridge
from sklearn.svm import SVR
from sklearn.linear_model import LogisticRegression
import progressbar

from easyesn.optimizers import GradientOptimizer
from easyesn.optimizers import GridSearchOptimizer


class StackedESN(BaseESN):
    #n_output1: for 1st layer
    #n_output2: for second layer
    def __init__(self, n_input, n_output1, n_output2, n_reservoir,
                 spectralRadius=1.0, noiseLevel=0.0, inputScaling=None,
                 leakingRate=1.0, feedbackScaling = 1.0, reservoirDensity=0.2, randomSeed=None,
                 out_activation=lambda x: x, out_inverse_activation=lambda x: x,
                 weightGeneration='naive', bias=1.0, outputBias=1.0,
                 outputInputScaling=1.0, inputDensity=1.0, solver='pinv', regressionParameters={}, activation = B.tanh, activationDerivation=lambda x: 1.0/B.cosh(x)**2):

        #probably not needed
        super(StackedESN, self).__init__(n_input=n_input, n_reservoir=n_reservoir, n_output=n_output1, spectralRadius=spectralRadius,
                                  noiseLevel=noiseLevel, inputScaling=inputScaling, leakingRate=leakingRate, feedbackScaling = feedbackScaling, reservoirDensity=reservoirDensity,
                                  randomSeed=randomSeed, feedback = False, out_activation=out_activation, out_inverse_activation=out_inverse_activation,
                                  weightGeneration=weightGeneration, bias=bias, outputBias=outputBias, outputInputScaling=outputInputScaling,
                                  inputDensity=inputDensity, activation=activation, activationDerivation=activationDerivation)

        #layers 
        self._esn1 = PredictionESN(n_input= (n_input + n_output2), n_reservoir=n_reservoir, n_output= n_output1, spectralRadius=spectralRadius,
                                  noiseLevel=noiseLevel, inputScaling=inputScaling, leakingRate=leakingRate, feedbackScaling = feedbackScaling, reservoirDensity=reservoirDensity,
                                  randomSeed=randomSeed, feedback = False, out_activation=out_activation, out_inverse_activation=out_inverse_activation,
                                  weightGeneration=weightGeneration, bias=bias, outputBias=outputBias, outputInputScaling=outputInputScaling,
                                  inputDensity=inputDensity, activation=activation, activationDerivation=activationDerivation)

        #esn_2 takes coordinates in (estimated and actual)
        self._esn2 = PredictionESN(n_input=n_output1*2, n_reservoir=n_reservoir, n_output=n_output2, spectralRadius=spectralRadius,
                                  noiseLevel=noiseLevel, inputScaling=inputScaling, leakingRate=leakingRate, feedbackScaling = feedbackScaling, reservoirDensity=reservoirDensity,
                                  randomSeed=randomSeed, feedback = False, out_activation=out_activation, out_inverse_activation=out_inverse_activation,
                                  weightGeneration=weightGeneration, bias=bias, outputBias=outputBias, outputInputScaling=outputInputScaling,
                                  inputDensity=inputDensity, activation=activation, activationDerivation=activationDerivation)

        self._solver = solver
        self._regressionParameters = regressionParameters
        self._transientTime = 0

        self._n_output1 = n_output1
        self._n_output2 = n_output2
        self._training_res = B.empty((2,2))


        """
            allowed values for the solver:
                pinv
                lsqr (will only be used in the thesis)

                sklearn_auto
                sklearn_svd
                sklearn_cholesky
                sklearn_lsqr
                sklearn_sag
        """
    """
        Fits the ESN so that by applying the inputData the outputData will be produced.
    """
    def fit(self, inputData, outputData1, outputData2, transientTime="AutoReduce", transientTimeCalculationEpsilon = 1e-3, transientTimeCalculationLength = 20, verbose=0):
        #check the input data
        if self.n_input != 0:
            if len(inputData.shape) == 3 and len(outputData1.shape) > 1:
                #multiple time series are used with a shape (timeseries, time, dimension) -> (timeseries, time, dimension)
                if inputData.shape[0] != outputData1.shape[0]:
                    raise ValueError("Amount of input and output datasets is not equal - {0} != {1}".format(inputData.shape[0], outputData1.shape[0]))
                if inputData.shape[1] != outputData1.shape[1]:
                    raise ValueError("Amount of input and output time steps is not equal - {0} != {1}".format(inputData.shape[1], outputData1.shape[1]))
            else:
                if inputData.shape[0] != outputData1.shape[0]:
                    raise ValueError("Amount of input and output time steps is not equal - {0} != {1}".format(inputData.shape[0], outputData1.shape[0]))
        else:
            if inputData is not None:
                raise ValueError("n_input has been set to zero. Therefore, the given inputData will not be used.")

        inputData = B.array(inputData)
        outputData1 = B.array(outputData1)      
        outputData2 = B.array(outputData2)  
        self._transientTime = transientTime
        
        self._esn1.resetState()
        self._esn2.resetState()

        total_length = inputData.shape[0]
        print("total_length ", total_length)

        aggregated_y1 = B.empty((outputData1.shape))
        aggregated_y2 = B.empty((outputData2.shape))

        #put pixels and switch data together
        inputDataWithFeedback = B.zeros((total_length, self.n_input + self._n_output2))
        inputDataWithFeedback[:, :self.n_input] = inputData
        inputDataWithFeedback[:, self.n_input:] = outputData2

        X1, _ = self._esn1.propagate(inputData=inputDataWithFeedback, outputData=outputData1, transientTime=transientTime, verbose=verbose-1)
        self._esn1._X = X1

        Y_target = self.out_inverse_activation(outputData1).T
        self._esn1._WOut = B.dot(Y_target, B.pinv(X1))
        aggregated_y1 = self.out_activation((B.dot(self._esn1._WOut, X1)).T)
        training_error1 = B.sqrt(B.mean((aggregated_y1.reshape((total_length, self._n_output1))- outputData1)**2))
        training_error2 = 0

        training_error1 = B.sqrt(B.mean((aggregated_y1 - outputData1)**2))
        training_error2 = B.sqrt(B.mean((aggregated_y2 - outputData2)**2))

        return training_error1, training_error2


    """
        Use the ESN in the predictive mode to predict the output signal by using an input signal.
        Used for performance estimation after training
    """
    def predict(self, inputData, outputData1, outputData2, continuation=True, initialData=None, update_processor=lambda x:x, verbose=0):
        inputData = B.array(inputData)
        
        #let some input run through the ESN to initialize its states from a new starting value
        if (not continuation):
            self._esn1._x = B.zeros(self._esn1._x.shape)
            self._esn2._x = B.zeros(self._esn2._x.shape)

        total_length = inputData.shape[0]
        aggregated_y1 = B.empty((total_length, self._n_output1))
        aggregated_y2 = B.empty((total_length, self._n_output2))

        #put pixels and switch data together
        inputDataWithFeedback = B.zeros((total_length, self.n_input + self._n_output2))
        inputDataWithFeedback[:, :self.n_input] = inputData
        inputDataWithFeedback[:, self.n_input:] = outputData2

        X1, _ = self._esn1.propagate(inputData=inputDataWithFeedback, outputData=None, transientTime=self._transientTime, verbose=verbose-1)
        aggregated_y1 = B.dot(self._esn1._WOut, X1)
        aggregated_y1 = update_processor(self.out_activation(aggregated_y1)).T

        training_error1 = B.sqrt(B.mean((aggregated_y1 - outputData1)**2))

        print("predict errors")
        print(training_error1)

        return aggregated_y1, aggregated_y2


    """
        Fits the ESN so that by applying the inputData the outputData will be produced.
    """
    def fit_loop(self, inputData, outputData1, outputData2, transientTime="AutoReduce", transientTimeCalculationEpsilon = 1e-3, transientTimeCalculationLength = 20, verbose=0):
        #check the input data
        if self.n_input != 0:
            if len(inputData.shape) == 3 and len(outputData1.shape) > 1:
                #multiple time series are used with a shape (timeseries, time, dimension) -> (timeseries, time, dimension)
                if inputData.shape[0] != outputData1.shape[0]:
                    raise ValueError("Amount of input and output datasets is not equal - {0} != {1}".format(inputData.shape[0], outputData1.shape[0]))
                if inputData.shape[1] != outputData1.shape[1]:
                    raise ValueError("Amount of input and output time steps is not equal - {0} != {1}".format(inputData.shape[1], outputData1.shape[1]))
            else:
                if inputData.shape[0] != outputData1.shape[0]:
                    raise ValueError("Amount of input and output time steps is not equal - {0} != {1}".format(inputData.shape[0], outputData1.shape[0]))
        else:
            if inputData is not None:
                raise ValueError("n_input has been set to zero. Therefore, the given inputData will not be used.")

        inputData = B.array(inputData)
        outputData1 = B.array(outputData1)      
        outputData2 = B.array(outputData2)  
        self._transientTime = transientTime
        
        self._esn1.resetState()
        self._esn2.resetState()

        total_length = inputData.shape[0]
        print("total_length ", total_length)


        if (verbose > 0):
            bar = progressbar.ProgressBar(max_value=total_length, redirect_stdout=True, poll_interval=0.0001)
            bar.update(0)

        #should be named aggregated_y
        aggregated_y1 = B.empty((outputData1.shape))
        aggregated_y2 = B.empty((outputData2.shape))


        # X1, _ = self._esn1.propagate(inputData=inputData, outputData=outputData1, transientTime=transientTime, verbose=verbose-1)
        # self._esn1._X = X1

        # Y_target = self.out_inverse_activation(outputData1).T
        # self._esn1._WOut = B.dot(Y_target, B.pinv(X1))
        # y1 = self.out_activation((B.dot(self._esn1._WOut, X1)).T)
        # training_error1 = B.sqrt(B.mean((y1.reshape((total_length, self._n_output1))- outputData1)**2))
        # training_error2 = 0

        y2 = self.out_activation(B.dot(self._esn2._WOut, self._esn2._X).T)

        for i in range((total_length)):

            #input: input + output from layer 2
            inputDataWithFeedback = B.zeros((self.n_input + self._n_output2))
            inputDataWithFeedback[:self.n_input] = inputData[i,:]
            inputDataWithFeedback[self.n_input:] = y2

            #update models
            X1, untrained_y1 = self._esn1.propagateOneStep(inputData=inputDataWithFeedback, outputData=outputData1[i], step = i, transientTime=transientTime, verbose=verbose-1, learn = True)
            #self._esn1._X = X1
            #y after model update (ideally we would use y before model update)
            #y1 = self.out_activation((B.dot(self._esn1._WOut, X1)).T)
            aggregated_y1[i,:]  = untrained_y1.reshape(self._n_output1)

            #output from 1st layer and correct output
            # input2 = B.vstack((y1.reshape(self._n_output1), outputData1[i]))
            # x2, y2 = self._esn2.propagateOneStep(inputData=input2, outputData=outputData2[i], step=i, transientTime=transientTime, verbose=verbose-1, learn=True)
            # Y_target2[i,:]  = y2.reshape(self._n_output2)

            if (verbose > 0):
                bar.update(i)
        if (verbose > 0):
            bar.finish()

        self._training_res = aggregated_y1
        training_error1 = B.sqrt(B.mean((aggregated_y1 - outputData1)**2))
        training_error2 = B.sqrt(B.mean((aggregated_y2 - outputData2)**2))

        print("training errors")
        print(training_error1)
        print(training_error2)

        return training_error1, training_error2


    """
        Use the ESN in the predictive mode to predict the output signal by using an input signal.
        Used for performance estimation after training
    """
    def predict_loop(self, inputData, outputData1, outputData2, continuation=True, initialData=None, update_processor=lambda x:x, verbose=0):
        inputData = B.array(inputData)
        
        #let some input run through the ESN to initialize its states from a new starting value
        if (not continuation):
            self._esn1._x = B.zeros(self._esn1._x.shape)
            self._esn2._x = B.zeros(self._esn2._x.shape)


        total_length = inputData.shape[0]
        aggregated_y1 = B.empty((total_length, self._n_output1))
        aggregated_y2 = B.empty((total_length, self._n_output2))

        # X1, y1 = self._esn1.propagate(inputData=inputData, outputData=None, transientTime=self._transientTime, verbose=verbose-1)
        # Y1 = B.dot(self._esn1._WOut, X1)
        # Y1 = update_processor(self.out_activation(Y1)).T

        y2 = self.out_activation(B.dot(self._esn2._WOut, self._esn2._X).T)

        for i in range(total_length):

            inputDataWithFeedback = B.zeros((self.n_input + self._n_output2))
            inputDataWithFeedback[:self.n_input] = inputData[i,:]
            inputDataWithFeedback[self.n_input:] = y2

            #update models
            X1, _ = self._esn1.propagateOneStep(inputData=inputDataWithFeedback, outputData=None, step = i, transientTime=self._transientTime, verbose=verbose-1, learn=False)
            y1 = B.dot(self._esn1._WOut, X1)
            aggregated_y1[i,:]  = update_processor(self.out_activation(y1)).T

            #output from 1st layer and correct output
            #input2 = outputData1[i] - y1.reshape(self._n_output1)
            # input2 = B.vstack((y1.reshape(self._n_output1), outputData1[i]))
            # x2, y2 = self._esn2.propagateOneStep(inputData=input2, outputData=None, step=i, transientTime=self._transientTime, verbose=verbose-1, learn=False)
            # Y_target2[i,:]  = y2.reshape(self._n_output2)

        training_error1 = B.sqrt(B.mean((aggregated_y1 - outputData1)**2))

        print("training errors")
        print(training_error1)

        return aggregated_y1, aggregated_y2


    """
        Use the ESN in the generative mode to generate a signal autonomously.
        NO remove
    """
    def generate(self, n, inputData=None, initialOutputData=None, continuation=True, initialData=None, update_processor=lambda x:x, verbose=0):
        #initialOutputData is the output of the last step BEFORE the generation shall start, e.g. the last step of the training's output

        #check the input data
        #if (self.n_input != self.n_output):
        #    raise ValueError("n_input does not equal n_output. The generation mode uses its own output as its input. Therefore, n_input has to be equal to n_output - please adjust these numbers!")

        if inputData is not None:
            inputData = B.array(inputData)

        if initialOutputData is not None:
            initialOutputData = B.array(initialOutputData)

        if initialData is not None:
            initialData = B.array(initialData)

        if initialOutputData is None and initialData is None:
            raise ValueError("Either intitialOutputData or initialData must be different from None, as the network needs an initial output value")

        if initialOutputData is None and initialData is not None:
            initialOutputData = initialData[1][-1]

        if inputData is not None:
            inputData = B.array(inputData)
        if initialData is not None:
            initialData = B.array(initialData)

        #let some input run through the ESN to initialize its states from a new starting value
        if not continuation:
            self._x = B.zeros(self._x.shape)

            if initialData is not None:
                if type(initialData) is tuple:
                    initialDataInput, initialDataOutput = initialData 
                    if initialDataInput is not None and len(initialDataInput) != len(initialDataOutput):
                        raise ValueError("Length of the inputData and the outputData of the initialData tuple do not match.")
                else:
                    raise ValueError("initialData has to be a tuple consisting out of the input and the output data.")

                for t in range(initialDataInput.shape[0]):
                    super(PredictionESN, self).update(initialDataInput[t], initialDataOutput[t])

        if self.n_input != 0:
            if inputData is None:
                raise ValueError("inputData must not be None.")
            elif len(inputData) < n:
                raise ValueError("Length of inputData has to be >= n.")

        _, Y = self.propagate(inputData, None, verbose=verbose, steps=n, previousOutputData=initialOutputData)
        Y = update_processor(Y)
        
        #return the result
        return Y.T

    def optimize(self, trainingInput, trainingOutput, validationInput, validationOutput, verbose):
        gridSearch = GridSearchOptimizer()
        gradientOptimizer = GradientOptimizer()
        pipe = Pipeline(gridSearch, gradientOptimizer)

        pipe.fit(trainingInput, trainingOutput, validationInput, validationOutput, verbose)