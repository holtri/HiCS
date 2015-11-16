import sequtils
import future
import strutils
import unittest
import math
import sets
import hashes
import subspace
import hics
import utils

type
  BinarySubspace* = seq[int]
  BinarySolution* = tuple
    binarySubspace: BinarySubspace
    deviations: seq[float]
    dominationSet: HashSet[BinarySolution]
    dominatedCount: int
    rank: int
    crowdingDistance: float
  BinaryPopulation* = HashSet[BinarySolution]

proc hash*(x: BinarySolution): Hash =
  var h: Hash = 0
  h = h !& hash(x.binarySubspace)
  result = !$h

proc `==`*(x: BinarySolution, y:BinarySolution): bool =
  return x.hash == y.hash

template isBinary(b: int): bool =
  b == 0 or b == 1

proc flip(bit: int): int =
  assert bit.isBinary
  if bit == 0: 1 else: 0

proc flip(bit: int, prob: float): int =
  assert prob <= 1.0 and prob >= 0.0
  assert bit.isBinary
  if random(1.0) < prob:
    return flip(bit)
  else:
    return bit

proc asSeq*(pop: BinaryPopulation): seq[BinarySolution] =
  result = @[]
  for p in pop:
    result.add(p)

proc isValid*(s: BinarySubspace): bool =
  return foldl(s, a + b) > 1

proc randomBinarySubspace*(totalDim: int, proportion: float): BinarySubspace =
  return newSeqWith(totalDim, flip(0, proportion))

proc initBinarySolution*(binarySubspace: BinarySubspace, deviations: seq[float]): BinarySolution =
  return (binarySubspace, deviations, initSet[BinarySolution](),0,0,0.0)

proc initBinarySolution*(binarySubspace: BinarySubspace): BinarySolution =
  return initBinarySolution(binarySubspace, newSeqWith(len(binarySubspace), 0.0))

proc initBinaryPopulation*(populationSize: int): BinaryPopulation =
  return initSet[BinarySolution](nextPowerOfTwo(populationSize))

proc initRandomBinaryPopulation*(populationSize: int, totalDim: int, proportion: float): BinaryPopulation =
  result = initBinaryPopulation(populationSize)
  for i in 1..populationSize:
    var binarySubspace = randomBinarySubspace(totalDim, proportion)
    while not binarySubspace.isValid:
      binarySubspace = randomBinarySubspace(totalDim, proportion)
    let binarySolution = initBinarySolution(binarySubspace)
    result.incl(binarySolution)

proc dominates*(a: BinarySolution, b: BinarySolution): bool =
  assert a.binarySubspace.len == b.binarySubspace.len
  result = false
  for index in 0..high(a.binarySubspace):
    if a.binarySubspace[index] == 1:
      if a.deviations[index] < b.deviations[index]:
        result = false
        break
      elif a.deviations[index] > b.deviations[index]:
        result = true

proc onePointCrossover*(p1: BinarySubspace, p2:BinarySubspace, crossIndex: int): (BinarySubspace, BinarySubspace) =
  assert len(p1) == len(p2)
  assert crossIndex >= low(p1) and crossIndex <= high(p1)
  if crossIndex == 0:
    result = (p2,p1)
  elif crossIndex == high(p1):
    result = (p1,p2)
  else:
    var c1 = p1[0..crossIndex].concat(p2[(crossIndex + 1)..high(p2)])
    var c2 = p2[0..crossIndex].concat(p1[(crossIndex + 1)..high(p2)])
    result = (c1,c2)

proc bitStringMutation*(p: BinarySubspace, prob: float): BinarySubspace =
  result = p
  result.applyIt(flip(it, prob))

proc toReal(s: BinarySubspace): string =
  result = "["
  for index,value in s:
    if value == 1:
      addSep(result, startLen=len("["))
      add(result, intToStr(index))
  add(result, "]")

proc asSubspace*(s: BinarySubspace): Subspace =
  result.init()
  for index, value in s:
    if value == 1:
      result.incl(index)

proc `$`*(bs: BinarySolution): string =
  result ="["
  for index in 0..high(bs.binarySubspace):
    result &= $index & ":(" & $bs.binarySubspace[index] & "," & $bs.deviations[index] & ") "
  result.removeSuffix(" ")
  result &= "]"

proc toReal*(bs: BinarySolution): string =
  result = "["
  for index,value in bs.binarySubspace:
    if value == 1:
      result &= $index & ": " & $bs.deviations[index] & ", "
  result.removeSuffix(", ")
  result &= "]"

proc `$`*(bp: BinaryPopulation): string =
  result = ""
  for item in bp:
    result &= $item & ",\n"
  result.removeSuffix(",\n")

proc toReal*(bp: BinaryPopulation): string =
  result = ""
  for item in bp:
    result &= item.toReal & ",\n"
  result.removeSuffix(",\n")

#when isMainModule:

suite "Binary subspace testing":
  setup:
    var a = BinarySubspace(@[1,0,0,1,1])
    var b = BinarySubspace(@[0,1,1,0,0])

  test "onePointCrossovverTest1":
    let actual = onePointCrossover(a,b,0)
    check(actual[0] == b)
    check(actual[1] == a)

  test "onePointCrossovverTest2":
    let actual = onePointCrossover(a,b,high(a))
    check(actual[0] == a)
    check(actual[1] == b)

  test "onePointCrossovverTest3":
    let actual = onePointCrossover(a,b,2)
    check(actual[0] == @[1,0,0,0,0])
    check(actual[1] == @[0,1,1,1,1])

  test "flipBitTest1":
    let tmp = 1
    check(tmp.flip == 0)
    check(tmp.flip.flip == 1)

  test "flipBitTest2":
    let tmp = 0
    check(tmp.flip(1.0) == 1)
    check(tmp.flip(0) == 0)

  test "hashBinarySubspaceTest":
    check(hash(a) == hash(@[1,0,0,1,1]))
    check(hash(b) == hash(@[0,1,1,0,0]))

  test "randomBinarySubspaceTest":
    var actual = randomBinarySubspace(5, 0.3)
    check(len(actual) == 5)
    actual = @[1,1,1]
    check(actual == randomBinarySubspace(3,1.0))

  test "dominatesTest":
    var c: BinarySolution = initBinarySolution(@[1,0,0,1,1],@[0.6,0,0,0.5,0.5])
    var d: BinarySolution = initBinarySolution(@[1,0,0,1,1],@[0.5,0,0,0.5,0.5])
    check(c.dominates(d))
    check(not d.dominates(c))

    c = initBinarySolution(@[1,0,0,1,1],@[0.6,0,0,0.4,0.5])
    check(not c.dominates(d))
    check(not d.dominates(c))

    c = initBinarySolution(@[1,0,0,0,1],@[0.6,0,0,0,0.7])

    check(c.dominates(d))
    check(not d.dominates(c))

    c = initBinarySolution(@[1,0,0,0,1],@[0.5,0,0,0.5,0.5])
    check(not c.dominates(d))
    check(not d.dominates(c))

  test "asSeqTest":
    var e: BinarySolution = initBinarySolution(@[1,1,0],@[0.5,0.5,0.0])
    var f: BinarySolution = initBinarySolution(@[0,1,1],@[0.0,0.5,0.5])
    var g: BinarySolution = initBinarySolution(@[1,0,1],@[0.5,0.0,0.5])

    var pop: BinaryPopulation = initBinaryPopulation(5)
    pop.incl(e)
    pop.incl(f)
    pop.incl(g)

    var seqPop = pop.asSeq
    check(seqPop.len == pop.len)

  test "toSubspace":
    let actual = a.asSubspace
    check(actual.len == 3)
    check(actual.contains(0))
    check(actual.contains(3))
    check(actual.contains(4))

  test "isValid":
    var h = BinarySubspace(@[0,0,1,0,0])
    check(a.isValid)
    check(not h.isValid)

    
