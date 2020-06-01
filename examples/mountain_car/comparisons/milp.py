'''

Authors: Radoslav Ivanov, Taylor J Carpenter, James Weimer, Rajeev Alur, George J. Pappa, Insup Lee

This script computes the output range of a neural network given input
constraints. It takes a neural network (encoded in the YAML format) as
a command line input. The input constraints can be modified in the
main() function below. If one wishes to use this script on a different
example, one would have to modify the inputBounds list accordingly.

INSTALLATION NOTE: this script uses the gurobipy library from the
mixed integer linear optimization tool Gurobi. Gurobi is free to use
for academic purposes but one still needs to obtain a licence from the
Gurobi website:
https://www.gurobi.com/academia/academic-program-and-licenses/

Example usage: python milp.py ../sig16x16.yml

'''


from gurobipy import *
import yaml
import numpy as np
import time

def computeReachableOutput(weights, offsets, activations, inputBounds, eps, maxStepSize, maximize):
    
    try:

        numLayers = len(offsets)

        # Create a new model
        m = Model("nn_" + str(len(offsets[1])) + 'x' + str(len(offsets[1])))

        #add variables-------------------------------------------------------------------------------------------------
        varDic = {}

        #add input variables
        for i in range(len(inputBounds)):
            varName = 'x' + str(i+1)
            varDic[varName] = m.addVar(lb=-GRB.INFINITY,vtype=GRB.CONTINUOUS, name=varName)

        #add neuron variables
        for layer in offsets:
            for neuron in range(len(offsets[layer])):
                contVarName = 'f' + str(layer) + '_' + str(neuron + 1)
                varDic[contVarName] = m.addVar(lb=-GRB.INFINITY,vtype=GRB.CONTINUOUS, name=contVarName)

                # if activations[layer] == 'Sigmoid' or activations[layer] == 'Tanh':
                #     for sigPiece in range(numLinPieces):
                #         binVarName = 'b' + str(layer) + '_' + str(neuron + 1) + '_' + str(sigPiece + 1)
                #         varDic[binVarName] = m.addVar(vtype=GRB.BINARY, name=binVarName)


        #set objective-------------------------------------------------------------------------------------------------
        #For now, assume we are interested in the first output variable
        outVarName = 'f' + str(numLayers) + '_1'

        if maximize:
            m.setObjective(varDic[outVarName], GRB.MAXIMIZE)
        else:
            m.setObjective(varDic[outVarName], GRB.MINIMIZE)

        #input constraints---------------------------------------------------------------------------------------------
        count = 1
        for var in inputBounds:
            varName = 'x' + str(count)
            m.addConstr(varDic[varName] >= var[0], 'i' + str(count) + 'l')
            m.addConstr(varDic[varName] <= var[1], 'i' + str(count) + 'u')

            count += 1

        #neuron constraints--------------------------------------------------------------------------------------------

        curBounds = inputBounds

        for layer in offsets:

            #print(curBounds)

            newBounds = []
            
            for neuron in range(len(offsets[layer])):
                curWeights = weights[layer][neuron]
                curOffset = offsets[layer][neuron]
                curActivation = activations[layer]

                Mx = getInputMx(curBounds, curWeights, curOffset)
                #print(Mx)

                appr = getEpsLinAppr(eps, maxStepSize, curActivation, Mx)
                linPieces = appr[0]
                apprPoints = appr[1]

                My = [0,0]
                if activations[layer] == 'Sigmoid':
                    My[0] = 1/(1 + np.exp(-Mx[0]))
                    My[1] = 1/(1 + np.exp(-Mx[1]))
                elif activations[layer] == 'Tanh':
                    My = np.tanh(Mx)
                else:
                    My = Mx

                addNeuronConstraints(m, curWeights, curOffset, curActivation, varDic,\
                                     layer, neuron + 1, linPieces, apprPoints, eps, Mx, My)

                newBounds.append(My)

            curBounds = newBounds

        m.write('temp.lp')
        
        m.optimize()

    except GurobiError as e:
        print('Error code ' + str(e.errno) + ": " + str(e))

    except AttributeError as e2:
        print('Encountered an attribute error: ' + str(e2))

def addNeuronConstraints(gurobiModel, weights, offset, activation, varDic, layerInd,
                         neurInd, linPieces, apprPoints, eps, Mx, My):
    
    if activation == 'Sigmoid' or activation == 'Tanh':

        #binary constraint
        binConstr = LinExpr()
        
        for piece in range(len(apprPoints)):

            #binary
            binVarName = 'b' + str(layerInd) + '_' + str(neurInd) + '_' + str(piece + 1)
            varDic[binVarName] = gurobiModel.addVar(vtype=GRB.BINARY, name=binVarName)
            binConstr += varDic[binVarName]
            
            #input constraints
            inConstrL = LinExpr()
            inConstrU = LinExpr()

            #output constraints
            outConstrL = LinExpr()
            outConstrU = LinExpr()            

            neurVarName = 'f' + str(layerInd) + '_' + str(neurInd)

            outConstrL += varDic[neurVarName]
            outConstrU += varDic[neurVarName]

            for weight in range(len(weights)):
                if layerInd == 1:
                    varName = 'x' + str(weight + 1)
                else:
                    varName = 'f' + str(layerInd - 1) + '_' + str(weight + 1)

                #input
                inConstrL += weights[weight] * varDic[varName]
                inConstrU += weights[weight] * varDic[varName]

                #output
                outConstrL -= linPieces[piece][0] * weights[weight] * varDic[varName]
                outConstrU -= linPieces[piece][0] * weights[weight] * varDic[varName]

            #input
            inConstrL += offset
            inConstrU += offset

            #output
            outConstrL -= linPieces[piece][0] * offset
            outConstrU -= linPieces[piece][0] * offset
            outConstrL -= linPieces[piece][1]
            outConstrU -= linPieces[piece][1]

            #input
            inConstrL -= (1 - varDic[binVarName]) * (Mx[0] - Mx[1])
            inConstrU -= (1 - varDic[binVarName]) * (Mx[1] - Mx[0])

            #output
            outConstrL -= (1 - varDic[binVarName]) * (My[0] - My[1])
            outConstrU -= (1 - varDic[binVarName]) * (My[1] - My[0])

            #input
            gurobiModel.addConstr(inConstrL, GRB.GREATER_EQUAL, apprPoints[piece][0], 'f' + str(layerInd) + '_' + str(neurInd) + 'in' + '_' + str(piece + 1) + 'l')
            gurobiModel.addConstr(inConstrU, GRB.LESS_EQUAL, apprPoints[piece][1], 'f' + str(layerInd) + '_' + str(neurInd) + 'in' + '_' + str(piece + 1) + 'u')

            #output
            gurobiModel.addConstr(outConstrL, GRB.GREATER_EQUAL, -eps, 'f' + str(layerInd) + '_' + str(neurInd) + 'out' + '_' + str(piece + 1) + 'l')
            gurobiModel.addConstr(outConstrU, GRB.LESS_EQUAL, eps, 'f' + str(layerInd) + '_' + str(neurInd) + 'out' + '_' + str(piece + 1) + 'u')

        gurobiModel.addConstr(binConstr, GRB.EQUAL, 1, 'b' + str(layerInd) + '_' + str(neurInd))

    else:
        constr = LinExpr()

        neurVarName = 'f' + str(layerInd) + '_' + str(neurInd)

        constr += varDic[neurVarName]

        for weight in range(len(weights)):
            if layerInd == 1:
                varName = 'x' + str(weight + 1)
            else:
                varName = 'f' + str(layerInd - 1) + '_' + str(weight + 1)

            #input
            constr -= weights[weight] * varDic[varName]

        constr -= offset

        gurobiModel.addConstr(constr, GRB.EQUAL, 0.0, 'f' + str(layerInd) + '_' + str(neurInd) + 'out' + '_l')

def getEpsLinAppr(eps, maxStepSize, activation, Mx):
    if activation == 'Sigmoid':
        appr = getSigAppr(eps, maxStepSize, activation, Mx)
    elif activation == 'Tanh':
        appr = getTanhAppr(eps, maxStepSize, activation, Mx)
    else:
        return (0, 0)

    linearPieces = appr[0]
    intervals = appr[1]

    outPieces = []
    outIntervals = []

    count = 0
    for interval in intervals:
        count += 1
        
        if interval[1] < Mx[0]:
            continue

        if interval[0] > Mx[1]:
            break

        if interval[1] >= Mx[1] and interval[0] <= Mx[1] and interval[0] >= Mx[0]:
            outIntervals.append((interval[0], Mx[1]))
            outPieces.append(linearPieces[count-1])

        elif Mx[1] >= interval[1] and Mx[0] <= interval[1] and Mx[0] >= interval[0]:
            outIntervals.append((Mx[0], interval[1]))
            outPieces.append(linearPieces[count-1])

        elif Mx[1] <= interval[1] and Mx[0] >= interval[0]:
            outIntervals.append((Mx[0], Mx[1]))
            outPieces.append(linearPieces[count-1])

        else:
            outIntervals.append(interval)
            outPieces.append(linearPieces[count-1])

    return (outPieces, outIntervals)

def getSigAppr(eps, maxStepSize, activation, Mx):

    apprMargin = 1e-2 * eps

    startPoint = -max(abs(Mx[0] - 1), abs(Mx[1] + 1))
    endPoint = 0

    argmax2ndDer = -1.317
    max2ndDer = 0.1

    allPoints = []
    allSlopes = []

    curPoint = startPoint

    while (curPoint <= endPoint):
        sig = 1/(1 + np.exp(-curPoint))
        slope = sig * (1 - sig)

        allPoints.append(curPoint)
        allSlopes.append(slope)
   
        if curPoint <= argmax2ndDer:
            if curPoint + maxStepSize <= argmax2ndDer:
                y = 1/(1 + np.exp(-curPoint - maxStepSize))
                M = y * (1 - y) * (1 - 2*y)
            else:
                M = max2ndDer
        else:
            y = 1/(1 + np.exp(-curPoint))
            M = y * (1 - y) * (1 - 2*y)   
   
        sqrtTerm = np.sqrt(M * (curPoint * curPoint * M - 2 * (M * curPoint * curPoint / 2 - eps)))
   
        maxX = (curPoint * M + sqrtTerm) / M
   
        if (maxX - curPoint > maxStepSize):
            curPoint = curPoint + maxStepSize - apprMargin
        else:
            curPoint = maxX - apprMargin

    allPoints.extend(np.flipud(allPoints) * -1)
    allSlopes.extend(np.flipud(allSlopes[0: len(allSlopes) - 1]))

    allApprPoints = []
    allPieces = []

    for i in range(len(allSlopes)):
        y = 1/(1 + np.exp(-allPoints[i]))
        inter = y - allSlopes[i] * allPoints[i]
        allPieces.append((allSlopes[i], inter))
        allApprPoints.append((allPoints[i], allPoints[i+1]))

    return (allPieces, allApprPoints)

def getTanhAppr(eps, maxStepSize, activation, Mx):

    apprMargin = 1e-2 * eps

    startPoint = -max(abs(Mx[0]), abs(Mx[1]))
    endPoint = 0

    argmax2ndDer = -0.6585
    max2ndDer = 0.8

    allPoints = []
    allSlopes = []

    curPoint = startPoint

    while (curPoint <= endPoint):

        slope = (1 - np.tanh(curPoint) * np.tanh(curPoint))

        allPoints.append(curPoint)
        allSlopes.append(slope)
   
        if curPoint <= argmax2ndDer:
            if curPoint + maxStepSize <= argmax2ndDer:
                y = np.tanh(curPoint + maxStepSize)
                M = -2 * y * (1 - y*y)
            else:
                M = max2ndDer
        else:
            y = np.tanh(curPoint)
            M = -2 * y * (1 - y*y)
   
        sqrtTerm = np.sqrt(M * (curPoint * curPoint * M - 2 * (M * curPoint * curPoint / 2 - eps)))
   
        maxX = (curPoint * M + sqrtTerm) / M
   
        if (maxX - curPoint > maxStepSize):
            curPoint = curPoint + maxStepSize - apprMargin
        else:
            curPoint = maxX - apprMargin

    allPoints.extend(np.flipud(allPoints) * -1)
    allSlopes.extend(np.flipud(allSlopes[0: len(allSlopes) - 1]))

    allApprPoints = []
    allPieces = []

    for i in range(len(allSlopes)):
        y = np.tanh(allPoints[i])
        inter = y - allSlopes[i] * allPoints[i]
        allPieces.append((allSlopes[i], inter))
        allApprPoints.append((allPoints[i], allPoints[i+1]))

    return (allPieces, allApprPoints)

def getInputMx(inBounds, weights, offset):
    lbSum = 0
    ubSum = 0

    varIndex = 0
    for inVar in inBounds:
        weight = weights[varIndex]
        if weight >= 0:
            lbSum += weight * inVar[0]
            ubSum += weight * inVar[1]
        else:
            lbSum += weight * inVar[1]
            ubSum += weight * inVar[0]

        varIndex += 1

    lb = lbSum + offset
    ub = ubSum + offset

    bounds = [lb, ub]

    return bounds


def main(argv):

    dnnYaml = argv[0]    

    # this eps controls the precision of the piecewise linear approximation
    eps = 0.0001

    # maxStepSize is a hyperparameter used in the sigmoid/tanh approximation
    maxStepSize = 1

    # change this to True to maximize the output range
    maximize = False

    # inputBouds stores the constraints on the DNN inputs as a list of tuples
    inputBounds = []
    inputBounds.append((-0.52, -0.5))
    inputBounds.append((0, 0))

    with open(dnnYaml, 'rb') as f:

        dnn = yaml.load(f)

        start = time.time()
        
        computeReachableOutput(dnn['weights'], dnn['offsets'], dnn['activations'], inputBounds, eps, maxStepSize, maximize)
        
        end = time.time()
        print(str(end - start) + ' seconds')

if __name__ == '__main__':
    main(sys.argv[1:])
