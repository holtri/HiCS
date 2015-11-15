import binarySubspace
import sets
import math
import sequtils

proc generateRandomPopulation(N: int, totalDim: int): BinaryPopulation =
  result = initBinaryPopulation(N)
  while(len(result) < N):
    result.incl(initBinarySolution(randomBinarySubspace(totalDim,0.1)))
  return result

when isMainModule:
  var a = generateRandomPopulation(2, 4)
  var test = a.mapIt($it & "\n")
  echo a
