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

proc calculateDeviations(population: SingleBinaryPopulation, ds: Dataset, preproData: PreproData, params: Parameters, statTest: KSTest, referenceDim): SingleBinaryPopulation =
  var pop = population
  for p in pop:
    pop[p].deviation = computeAverageDeviation(p.subspace.asSubspace, ds, preproData, params, statTest, referenceDim)
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

proc selectMatingPool(population: SingleBinaryPopulation): seq[SingleBinarySolution] =
  var pop = population.asSeq
  var matingPool:seq[SingleBinarySolution] = @[]

  while matingPool.len < population.len:
    let a = pop[random(population.len)]
    let b = pop[random(population.len)]
    if a.dominates(b):
      matingPool.add(a)
    else:
      matingPool.add(b)
  return matingPool

proc performMating(pop: seq[SingleBinarySolution]): SingleBinaryPopulation =

  var offsprings = initSet[SingleBinarySolution]()
  var i:int = 0
  let totalDim = pop[0].subspace.len
  let referenceDim = pop[0].referenceDim
  let mutationProb:float = 1.0 / totalDim

  while i < pop.high-1:
    let cross = onePointCrossover(pop[i].subspace, pop[i+1].subspace, random(totalDim))
    let os1 = cross[0].bitStringMutation(mutationProb, referenceDim)
    let os2 = cross[0].bitStringMutation(mutationProb, referenceDim)

    if bs.isValid(os1):
      offsprings.incl((os1,0.0, referenceDim))
    if bs.isValid(os2):
      offsprings.incl((os2,0.0, referenceDim))
    i+=2
  for os in offsprings:
    assert isValid(os.subspace)
  return offsprings


proc singleDimensionOptimization*(N: int, referenceDim: int, maxIteration: int, ds: Dataset, preproData: PreproData, params: Parameters, statTest: KSTest): SubspaceSet =
  debug referenceDim
  let totalDim = ds.ncols
  var parents = generateRandomPopulation(N,totalDim,referenceDim)
  parents = calculateDeviations(parents,ds, preproData, params, statTest, referenceDim)
  let sorted = parents.sorted

  for i in 1..maxIteration:
    echo ifmt("iteration $i \n---------------")

    let matingPool = selectMatingPool(parents)
    var offsprings = performMating(matingPool)
    offsprings = calculateDeviations(offsprings,ds, preproData, params, statTest, referenceDim)

    #(mu + lambda)-Reproduction
    let unionGeneration = sorted(parents.union(offsprings))
    parents = unionGeneration[0..N-1].asSet

    #output
    if i %% 10 == 0:
      for p in parents.sorted:
        debug p.toReal

when isMainModule:
  var a = @[1,2,3,4,5]
  var b = a[0..3]
  debug a, b


  
