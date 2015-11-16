import binarySubspace
import sets
import math
import sequtils
import tables
import unittest
import utils
import algorithm
import preprocessing
import binaryHeap
import dataset
import hics
import stattest
import subspace

proc generateRandomPopulation(N: int, totalDim: int): BinaryPopulation =
  result = initBinaryPopulation(N)
  while(len(result) < N):
    result.incl(initBinarySolution(randomBinarySubspace(totalDim,0.1)))
  return result

proc calculateDeviations(b: BinarySolution, ds: Dataset, preproData: PreproData, params: Parameters, statTest: KSTest): seq[float] =
  result = b.deviations
  for i in 0..b.deviations.high:
    if b.binarySubspace[i] == 1:
      let cmpAttr:int = i
      let subspace = b.binarySubspace.asSubspace
      result[i] = computeAverageDeviation(subspace, ds, preproData, params, statTest, cmpAttr)

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

proc crowdingDistance(population: BinaryPopulation, totalDim: int): BinaryPopulation =
  var pm = population.asSeq
  for i in 0..pm.high:
    pm[i].crowdingDistance = 0

  for i in 0..totalDim-1:
    pm = pm.sorted(proc (x,y: BinarySolution): int=
      result = - cmp(x.deviations[i], y.deviations[i]))

    var range = abs(pm[0].deviations[i] - pm[pm.high].deviations[i])
    #only increase distance if there is any positive deviation
    if range < 0.001:
      range = 1
    else:
      pm[0].crowdingDistance += 999
      pm[pm.high].crowdingDistance += 999

    for j in 1..pm.high-1:
      pm[j].crowdingDistance += abs(pm[j-1].deviations[i] - pm[j+1].deviations[i]) / range

  return pm.toSet

proc crowdingSelection(a:BinarySolution, b:BinarySolution): BinarySolution =
  if a.rank < b.rank:
    return a
  elif a.rank > b.rank:
    return b
  else:
    if a.crowdingDistance > b.crowdingDistance:
      return a
    else:
      return b

proc calculateDeviations(population: BinaryPopulation, ds: Dataset, preproData: PreproData, params: Parameters, statTest: KSTest): BinaryPopulation =
  var pop = population
  for p in pop:
    #let a = calculateDeviations(p, ds, preproData, params, statTest)
    pop[p].deviations = calculateDeviations(p, ds, preproData, params, statTest)
  return pop

proc runNsga*(N: int, ds: Dataset, preproData: PreproData, params: Parameters, statTest: KSTest): SubspaceSet =
  let totalDim = ds.ncols
  var pop = initRandomBinaryPopulation(N, totalDim, 0.1)
  debug pop
  pop = calculateDeviations(pop,ds, preproData, params, statTest)
  debug pop
  return result



discard """ when isMainModule:

  var a: BinarySolution = (BinarySubspace(@[1,0,0,1,1]),@[0.5,0,0,0.5,0.5])
  echo a
 """
suite "nsgaii testing":

  test "crowdingSelection":
    var a: BinarySolution = initBinarySolution(@[1,1,1,0])
    var b: BinarySolution = initBinarySolution(@[1,1,1,1])
    var c: BinarySolution = initBinarySolution(@[0,1,1,0])

    a.rank = 1
    b.rank = 2
    c.rank = 2

    a.crowdingDistance = 4.3
    b.crowdingDistance = 3.3
    c.crowdingDistance = 2.3
    check(crowdingSelection(a,b) == a)
    check(crowdingSelection(b,c) == b)

  test "crowdingDistance":
    var a: BinarySolution = initBinarySolution(@[1,1,1,0],@[0.8,0.6,0.2,0.0]) #always ranked in the middle
    var b: BinarySolution = initBinarySolution(@[1,1,1,1],@[0.9,0.5,0.3,0.0]) #ranked 2x top, 1xbottom
    var c: BinarySolution = initBinarySolution(@[0,1,1,0],@[0.0,0.8,0.1,0.0]) #ranked 1x top, 2xbottom

    var pop: BinaryPopulation = initBinaryPopulation(5)

    pop.incl(a)
    pop.incl(b)
    pop.incl(c)

    var crowdingPop = crowdingDistance(pop,4)

    check(crowdingPop[a].crowdingDistance ==3)
    check(crowdingPop[b].crowdingDistance ==2997)
    check(crowdingPop[c].crowdingDistance ==2997)

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

    
