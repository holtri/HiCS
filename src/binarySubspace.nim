import sequtils
import future
import strutils
import unittest
import math

type
  BinarySubspace* = seq[int]

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

proc toReal(s: BinarySubspace): string =
  result = "["
  for index,value in s:
    if value == 1:
      addSep(result, startLen=len("["))
      add(result, intToStr(index))
  add(result, "]")

proc onePointCrossover(p1: BinarySubspace, p2:BinarySubspace, crossIndex: int): (BinarySubspace, BinarySubspace) =
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

proc onePointMutation(p: BinarySubspace, prob: float): BinarySubspace =
  result = p
  result.applyIt(flip(it, prob))

when isMainModule:
  var nums = BinarySubspace(@[1,0,0,1,1])
  randomize()
  echo nums.onePointMutation(0.5)


suite "Binary subspace testing":
  setup:
    var a = BinarySubspace(@[1,0,0,1,1])
    var b = BinarySubspace(@[0,0,1,0,0])

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
    check(actual[1] == @[0,0,1,1,1])

  test "flipBit1":
    let tmp = 1
    check(tmp.flip == 0)
    check(tmp.flip.flip == 1)

  test "flipBit2":
    let tmp = 0
    check(tmp.flip(1.0) == 1)
    check(tmp.flip(0) == 0)
