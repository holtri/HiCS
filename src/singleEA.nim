import binarySubspace as bs
import nsgaii
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
import strutils

type
  SingleBinarySolution = tuple
    subspace: BinarySubspace
    deviation: float
    referenceDim: int
  SingleBinaryPopulation = HashSet[SingleBinarySolution]

proc randomSingleBinarySubspace(totalDim: int, proportion: float, referenceDim): BinarySubspace =
  var subspace = newSeqWith(totalDim, flip(0, proportion))
  subspace[referencedim] = 1
  while not subspace.isValid:
    subspace = newSeqWith(totalDim, flip(0, proportion))
    subspace[referencedim] = 1
  return subspace

proc calculateDeviations(population: seq[SingleBinarySolution], ds: Dataset, preproData: PreproData, params: Parameters, statTest: KSTest, referenceDim): seq[SingleBinarySolution] =
  var pop = population
  for i in 0..pop.high:
    pop[i].deviation = computeAverageDeviation(pop[i].subspace.asSubspace, ds, preproData, params, statTest, referenceDim)
  return pop

proc generateRandomPopulation(N: int, totalDim: int, referenceDim: int): SingleBinaryPopulation =
  result = initSet[SingleBinarySolution]()
  while(len(result) < N):
    let subspace = randomSingleBinarySubspace(totalDim,1.0/totalDim, referenceDim)
    result.incl((subspace,random(1.0),referenceDim))
  return result

proc dominates(a: SingleBinarySolution, b: SingleBinarySolution): bool =
  assert a.subspace.len == b.subspace.len
  if a.deviation > b.deviation:
    return true
  else:
    return false

proc asSeq(pop: SingleBinaryPopulation): seq[SingleBinarySolution] =
  result = @[]
  for p in pop:
    result.add(p)

proc asSet(pop: seq[SingleBinarySolution]): SingleBinaryPopulation =
  result = initSet[SingleBinarySolution]()
  for p in pop:
    result.incl(p)

proc sorted(pop: SingleBinaryPopulation): seq[SingleBinarySolution] =
  var sortedPop = pop.asSeq

  sortedPop = sortedPop.sorted(proc (x,y: SingleBinarySolution): int=
      result = -cmp(x.deviation, y.deviation))
  return sortedPop

proc sorted(pop: seq[SingleBinarySolution]): seq[SingleBinarySolution] =
  var sortedPop = pop
  sortedPop = sortedPop.sorted(proc (x,y: SingleBinarySolution): int=
      result = -cmp(x.deviation, y.deviation))
  return sortedPop

proc toReal*(bs: SingleBinarySolution): string =
  result = "["
  for index,value in bs.subspace:
    if value == 1:
      result &= $index & ", "
  result.removeSuffix(", ")
  result &= "] : " & $bs.deviation & "\n"


proc `$`(x: SingleBinaryPopulation): string =
  result = "\n"
  for p in x:
    result &= $p.subspace & " " & $p.deviation & "\n"

proc `$`(x: seq[SingleBinarySolution]): string =
  result = "\n"
  for p in x:
    result &= $p.subspace & " " & $p.deviation & "\n"

proc bitStringMutation*(p: BinarySubspace, prob: float, referenceDim: int): BinarySubspace =
  result = p
  for i in 0..p.high:
    if not i==referenceDim:
      result[i] = flip(result[i], prob)

proc tournamentSelection(pop: seq[SingleBinarySolution], number: int): SingleBinarySolution =
  var candidates:seq[SingleBinarySolution] = @[]
  while candidates.len <= number:
    candidates.add(pop[random(pop.len)])
  return candidates.sorted[0]

proc stochasticUniversalSampling(pop: seq[SingleBinarySolution], sizeToKeep: int): seq[SingleBinarySolution] =
  var totalFitness = 0.0
  var keep:seq[SingleBinarySolution] = @[]

  for i in 0..pop.high:
    totalFitness += pop[i].deviation
  let v = random(totalFitness / sizeToKeep)

  var popm = pop.sorted

  var k = 1
  var sum = popm[0].deviation
  var indexPointer = 0
  # debug totalFitness, sizeToKeep, v, pop.len
  while keep.len <= sizeToKeep:
    while k*v < sum:
      keep.add(popm[indexPointer])
      k += 1
    indexPointer += 1
    sum += popm[indexPointer].deviation

  return keep

proc selectMatingPool(population: seq[SingleBinarySolution]): seq[SingleBinarySolution] =
  var pop = population
  var matingPool:seq[SingleBinarySolution] = @[]

  while matingPool.len < population.len:
    matingPool.add(tournamentSelection(pop,5))

  return matingPool

proc performMating(pop: seq[SingleBinarySolution]): seq[SingleBinarySolution] =

  var offsprings:seq[SingleBinarySolution] = @[]

  let totalDim = pop[0].subspace.len
  let referenceDim = pop[0].referenceDim
  let mutationProb:float = 0.01#5.0 / totalDim

  while offsprings.len < pop.len:
    let cross = onePointCrossover(pop[random(pop.len)].subspace, pop[random(pop.len)].subspace, random(totalDim))
    let os1 = cross[0].bitStringMutation(mutationProb, referenceDim)
    let os2 = cross[0].bitStringMutation(mutationProb, referenceDim)

    if bs.isValid(os1):
      offsprings.add((os1,0.0, referenceDim))
    if bs.isValid(os2):
      offsprings.add((os2,0.0, referenceDim))
  for os in offsprings:
    assert isValid(os.subspace)
  return offsprings


proc singleDimensionOptimization*(N: int, referenceDim: int, maxIteration: int, ds: Dataset, preproData: PreproData, params: Parameters, statTest: KSTest): SubspaceSet =
  debug referenceDim, N
  let totalDim = ds.ncols
  #generation p(0)
  var parents = generateRandomPopulation(N,totalDim,referenceDim).asSeq
  parents = calculateDeviations(parents,ds, preproData, params, statTest, referenceDim)

  var i = 0
  var stop = false
  var counter = 0
  var bestSolution:SingleBinarySolution = parents.sorted[0]
  while ((not stop) or (i > maxIteration)):
    #echo ifmt("iteration $i \n---------------")

    let matingPool = selectMatingPool(parents)
    var offsprings = performMating(matingPool)
    offsprings = calculateDeviations(offsprings,ds, preproData, params, statTest, referenceDim)

    #(mu + lambda)-Reproduction
    let unionGeneration = sorted(parents.concat(offsprings))

    if abs(bestSolution.deviation - unionGeneration[0].deviation) < 0.001:
      counter += 1
    else:
      counter = 0
    if counter > 5:
      stop = true

    bestSolution = unionGeneration[0]
    #debug bestSolution.toReal
    parents = unionGeneration[0..N-1]
    #parents = @[]
    #while parents.len <= N:
    #  parents.add(tournamentSelection(unionGeneration, 2))
    #parents = offsprings
    #bestSolution = parents[0]
    #parents = stochasticUniversalSampling(unionGeneration, N)
    #output
    discard """ if i %% 10 == 0:
      for p in parents.sorted:
        debug p.toReal
    inc(i) """
  debug bestSolution.toReal

when isMainModule:
  var parents = generateRandomPopulation(10,5,0).asSeq
  let a = stochasticUniversalSampling(parents)
  debug a


  
