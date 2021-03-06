
import sets
import tables
import hashes
import sequtils
import utils
import math
import algorithm

# to allow s.incl(i)
export sets.incl

type
  Subspace* = HashSet[int]
  SubspaceSet = HashSet[Subspace]

when NimVersion == "0.11.2":
  type
    HashAlias = THash
else:
  type
    HashAlias = Hash

proc hash*(s: Subspace): HashAlias =
  #echo "chalculating hash for ", s
  for it in items(s):
    #echo "   hash item is ", it
    #result = result !& hash(it)
    result = result +% hash(it)
  result = !$result

proc toSubspace*(s: openarray[int]): Subspace =
  s.toSet

proc toSubspace*(s: Slice[int]): Subspace =
  result.init()
  for dim in s.a .. s.b:
    result.incl(dim)

proc dimensionality*(s: Subspace): int =
  s.len

proc asSeq*(s: Subspace): seq[int] =
  ## ensures dimensions are ordered
  result = sequtils.toSeq(items(s))
  result.sort(system.cmp)

proc `$`*(s: Subspace): string =
  let sortedStr = $s.asSeq
  sortedStr[1..^1] # drop the @

proc randomDim*(s: Subspace): int {.inline.} =
  let randomI = random(s.len)
  var i = 0
  for x in s:
    if i == randomI:
      return x
    i.inc
  assert(false, "randomDim should always return")



iterator lowDimProjections*(s: Subspace): Subspace =
  for dim in s:
    var space = s
    space.excl(dim)
    yield space


proc newSubspaceSet*(): SubspaceSet =
  result = initSet[Subspace]()

proc generate2DSubspaces*(D: int): SubspaceSet =
  result = initSet[Subspace]()
  ijForLoop(D):
    let subspace = [i,j].toSubspace
    result.incl(subspace)




proc aprioriMerge*(subspaces: SubspaceSet): SubspaceSet =
  var prefixSuffixMap = initTable[Subspace, seq[int]]()

  # generate a map of prefixes and suffixes
  for subspace in subspaces:
    let subspaceSeq = subspace.asSeq
    if subspaceSeq.len == 0:
      continue
    let prefix = subspaceSeq[0..^2].toSubspace
    let suffix = subspaceSeq[^1]

    try:
      prefixSuffixMap.mget(prefix).add(suffix)
    except KeyError:
      prefixSuffixMap[prefix] = @[suffix]
    #debug subspace, prefix, suffix, prefixSuffixMap[prefix]

  #debug prefixSuffixMap

  # generate candidates
  # sorting of suffices like in original Apriori is actually not necessary,
  # since we only have to ensure to iterate over unique suffix pairs.
  # This is accomplished here by using a special for loop, which ensures i < j.
  var candidates: seq[Subspace] = @[]
  for prefix, suffices in prefixSuffixMap:
    #debug prefix, suffices
    ijForLoop(suffices.len):
      #debug suffices[i], suffices[j]
      var space: Subspace = prefix
      space.incl(suffices[i])
      space.incl(suffices[j])
      candidates.add(space)
  #debug candidates

  # check all lower-dim projections:
  var validCandidates: seq[Subspace] = @[]
  for candidate in candidates:
    var allTrue = true
    for lowDimProjection in candidate.lowDimProjections:
      #debug lowDimProjection, subspaces.contains(lowDimProjection)
      if not subspaces.contains(lowDimProjection):
        allTrue = false
        break
    if allTrue:
      validCandidates.add(candidate)
  #debug validCandidates
  
  return validCandidates.toSet



UnitTests("aprioriMerge"):

  test "order independence of hash function":
    var s1 = [].toSubspace
    var s2 = [].toSubspace
    s1.incl(0); s1.incl(1); s1.incl(2); s1.incl(3); s1.incl(4)
    s2.incl(0); s2.incl(1); s2.incl(2); s2.incl(4); s2.incl(3)
    check s1.hash == s2.hash
    block:
      for iter in 1 .. 100:
        let orig = toSeq(1..100)
        let perm = shuffle(orig)
        var s1 = [].toSubspace
        var s2 = [].toSubspace
        for x in orig: s1.incl(x)
        for x in perm: s2.incl(x)
        check s1.hash == s2.hash

  test "generate2DSubspaces":
    check generate2DSubspaces(1).len == 0
    check generate2DSubspaces(2).len == 1
    check generate2DSubspaces(3).len == 3
    check generate2DSubspaces(4).len == 6
    check generate2DSubspaces(5).len == 10
    check generate2DSubspaces(6).len == 15

  test "4-dim test":
    let s1 = [1,2,3].toSubspace
    let s2 = [1,2,4].toSubspace
    let s3 = [2,3,4].toSubspace
    let s4 = [1,3,4].toSubspace
    let res = aprioriMerge([s1,s2,s3,s4].toSet)
    check res.len == 1
    check ([1,2,3,4].toSubspace in res)

  test "5-dim test":
    let s1 = [0, 1, 3, 4].toSubspace
    let s2 = [0, 2, 3, 4].toSubspace
    let s3 = [0, 1, 2, 3].toSubspace
    let s4 = [1, 2, 3, 4].toSubspace
    let s5 = [0, 1, 2, 4].toSubspace 
    let res = aprioriMerge([s1,s2,s3,s4,s5].toSet)
    check res.len == 1
    check contains(res, ([0,1,2,3,4].toSubspace))

  test "created from lower dims (one step)":
    for N in 2..5:
      let base = (1..N).toSubspace
      var subspaces = initSet[Subspace]()
      for s in base.lowDimProjections:
        subspaces.incl(s)
      #debug subspaces
      let merged = aprioriMerge(subspaces)
      check merged.len == 1
    
  test "created from lower dims (multiple steps)":
    let N = 5
    let base = (1..N).toSubspace
    var subspaces = initSet[Subspace]()
    for s1 in base.lowDimProjections: # 4 dim
      for s2 in s1.lowDimProjections: # 3 dim
        for s3 in s2.lowDimProjections: # 2 dim
          subspaces.incl(s3)
    #debug subspaces, subspaces.len
    let merged = aprioriMerge(subspaces)
    assert(merged.len == 10) # a 5 dim subspace has 10 3-dim subspaces, created from the 10 2-dim subspaces

