from scipy import *

def findClosest(arr, num):
    closeind = None
    mindist = num
    for elem in enumerate(arr):
        if abs(elem[1]-num) < mindist:
            closeind = elem[0]
            mindist = abs(elem[1]-num)
    return (closeind, arr[closeind])

def coarseArr(arr, ref):
    newarr = []
    newarrInd = []
    for elem in ref:
        ind, num = findClosest(arr, elem)
        if ind != None:
            newarr.append(num)
            newarrInd.append(ind)
    return (newarr, newarrInd)
