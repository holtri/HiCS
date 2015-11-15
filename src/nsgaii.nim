import binarySubspace
import sets
import math
import sequtils
import tables
import unittest
import utils

proc generateRandomPopulation(N: int, totalDim: int): BinaryPopulation =
  result = initBinaryPopulation(N)
  while(len(result) < N):
    result.incl(initBinarySolution(randomBinarySubspace(totalDim,0.1)))
  return result

proc fastNonDominatedSort(population: BinaryPopulation): Table[int, BinaryPopulation] =
  var pm = population
  var F = initTable[int, BinaryPopulation]()

  var i = 1
  F[i] = initBinaryPopulation(population.len)

  for p in population:

    var Sp: BinaryPopulation = initBinaryPopulation(population.len)
    for q in population:
      if pm[p].dominates(pm[q]):
        Sp.incl(pm[q])
        #echo ifmt("$p dominates $q")
      elif q.dominates(pm[p]):
        pm[p].dominatedCount += 1
      pm[p].dominationSet = Sp
      pm[p].rank = 1

    if pm[p].dominatedCount == 0:
      F[i].incl(pm[p])

  while F[i].len > 0:
    var Q = initBinaryPopulation(population.len)
    for p in F[i]:
      for q in pm[p].dominationSet:
        pm[q].dominatedCount -= 1
        if pm[q].dominatedCount == 0:
          pm[q].rank = i + 1
          Q.incl(pm[q])
    inc(i)
    if Q.len > 0:
      F[i] = Q
    else: break
  return F

discard """ when isMainModule:

  var a: BinarySolution = (BinarySubspace(@[1,0,0,1,1]),@[0.5,0,0,0.5,0.5])
  echo a
 """
suite "nsgaii testing":

  test "fastNonDominatedSort":
    var a: BinarySolution = initBinarySolution(@[1,0,0,1,1],@[0.5,0.0,0.0,0.5,0.5]) #not dominated
    var b: BinarySolution = initBinarySolution(@[0,1,1,0,0],@[0.0,0.5,0.5,0.0,0.0]) #not dominated
    var c: BinarySolution = initBinarySolution(@[1,1,1,0,0],@[0.4,0.4,0.4,0.0,0.4]) #dominated by a and b
    var d: BinarySolution = initBinarySolution(@[1,1,1,1,1],@[0.3,0.3,0.4,0.6,0.6]) #dominated by b and c
    var e: BinarySolution = initBinarySolution(@[0,1,0,0,1],@[0.0,0.1,0.0,0.0,0.3]) #dominated by a, b, c, d

    var pop: BinaryPopulation = initBinaryPopulation(5)

    pop.incl(a)
    pop.incl(b)
    pop.incl(c)
    pop.incl(d)
    pop.incl(e)

    var actual = fastNonDominatedSort(pop)
    check(actual[1].len == 2)
    check(actual[2].len == 1)
    check(actual[3].len == 1)
    check(actual[4].len == 1)

    #check invaraiance on permutation
    pop = initBinaryPopulation(5)
    pop.incl(c)
    pop.incl(e)
    pop.incl(a)
    pop.incl(d)
    pop.incl(b)

    actual = fastNonDominatedSort(pop)
    check(actual[1].contains(a) and actual[1].contains(b))
    check(actual[2].contains(c))
    check(actual[3].contains(d))
    check(actual[4].contains(e))
    
