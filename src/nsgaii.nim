import binarySubspace as bs
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

proc generateRandomPopulation*(N: int, totalDim: int, prop: float): BinaryPopulation =
  result = initBinaryPopulation(N)
  while(len(result) < N):
    result.incl(initBinarySolution(randomBinarySubspace(totalDim,prop)))
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
  debug pm.len
  var i = 1
  F[i] = initBinaryPopulation(population.len)
  var minDominatedCount = 999
  for p in population:
    var Sp: BinaryPopulation = initBinaryPopulation(population.len)
    for q in population:

      if pm[p].dominates(pm[q]):
        Sp.incl(q)
      elif q.dominates(pm[p]):
        pm[p].dominatedCount += 1
    pm[p].dominationSet = Sp
    minDominatedCount = min(minDominatedCount, pm[p].dominatedCount)
    #debug pm[p].dominatedCount
  #debug minDominatedCount
  for p in population:
    if pm[p].dominatedCount == minDominatedCount:
      F[i].incl(pm[p])
      pm[p].rank = 1
      #debug pm[p].dominationSet.len

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
  var sum:int = 0
  for i in F.keys:
    sum += F[i].len
  debug sum
  assert sum == population.len
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
    pop[p].deviations = calculateDeviations(p, ds, preproData, params, statTest)
  return pop

proc selectMatingPool(population: BinaryPopulation): seq[BinarySolution] =
  var pop = population.asSeq
  var matingPool:seq[BinarySolution] = @[]#initBinaryPopulation(population.len)

  while matingPool.len < population.len:
    matingPool.add(crowdingSelection(pop[random(population.len)],pop[random(population.len)]))
  return matingPool

proc performMating(pop: seq[BinarySolution]): BinaryPopulation =

  var offsprings = initBinaryPopulation(pop.len)
  var i:int = 0
  let totalDim = pop[0].binarySubspace.len
  let mutationProb:float = 1.0 / totalDim

  while i < pop.high-1:
    let cross = onePointCrossover(pop[i].binarySubspace, pop[i+1].binarySubspace, random(totalDim))
    let os1 = cross[0].bitStringMutation(mutationProb)
    let os2 = cross[0].bitStringMutation(mutationProb)

    if bs.isValid(os1):
      offsprings.incl(initBinarySolution(os1))
    if bs.isValid(os2):
      offsprings.incl(initBinarySolution(os2))
    i+=2
  for os in offsprings:
    assert bs.isValid(os.binarySubspace)
  return offsprings

proc sortByCrowdingDistance(pop: BinaryPopulation): seq[BinarySolution] =
  var popm = pop.asSeq
  popm = popm.sorted(proc (x,y: BinarySolution): int=
      result = - cmp(x.crowdingDistance, y.crowdingDistance))
  return popm

proc mergeFronts(paretoFronts: Table[int, BinaryPopulation], N: int): BinaryPopulation =

  var pop = initBinaryPopulation(N)
  var i = 1
  while paretoFronts[i].len + pop.len <= N:
    pop = pop.union(paretoFronts[i])
    if i == paretoFronts.len:
      break
    else:
      inc(i)

  let sortedFront = sortByCrowdingDistance(paretoFronts[i])
  var j = 0
  while pop.len <= N:
    pop.incl(sortedFront[j])
    inc(j)
  return pop

proc resetPopulation(pop: BinaryPopulation): BinaryPopulation =
  var pm = pop
  for p in pop:
    pm[p].dominationSet = initSet[BinarySolution]()
    pm[p].dominatedCount = 0
    pm[p].rank = 0
    pm[p].crowdingDistance = 0.0
  return pm

proc runNsga*(N: int, ds: Dataset, preproData: PreproData, params: Parameters, statTest: KSTest): SubspaceSet =
  let totalDim = ds.ncols

  var parents = initRandomBinaryPopulation(N, totalDim, 1/N)
  parents = calculateDeviations(parents,ds, preproData, params, statTest)

  let maxIteration = 100

  for i in 1..maxIteration:
    echo ifmt("iteration $i \n----------------")
    debug parents.len
    let matingPool = selectMatingPool(parents)
    var offsprings = performMating(matingPool)
    offsprings = calculateDeviations(offsprings,ds, preproData, params, statTest)
    let unionGeneration = parents.union(offsprings)
    debug unionGeneration.len
    let paretoFronts = fastNonDominatedSort(unionGeneration)
    parents = mergeFronts(paretoFronts,N)
    parents = resetPopulation(parents)
    if i %% 10 == 0:
      for p in parents:
        debug p.toReal
  return result

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

    
