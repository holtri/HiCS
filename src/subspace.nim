
import sets
import tables
import hashes
import sequtils
import utils

type
  Subspace* = HashSet[int]
  SubspaceSet = HashSet[Subspace]


proc hash*(s: Subspace): THash =
  var h: THash = 0
  for xAtom in s:
    h = h !& xAtom
  result = !$h


proc toSubspace*(s: openarray[int]): Subspace =
  s.toSet

proc toSubspace*(s: Slice[int]): Subspace =
  result.init()
  for dim in s.a .. s.b:
    result.incl(dim)

proc asSeq*(s: Subspace): seq[int] =
  result = sequtils.toSeq(items(s))

iterator lowDimProjections(s: Subspace): Subspace =
  for dim in s:
    var space = s
    space.excl(dim)
    yield space



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


runUnitTest("aprioriMerge"):
  let s1 = [1,2,3].toSubspace
  let s2 = [1,2,4].toSubspace
  let s3 = [2,3,4].toSubspace
  let s4 = [1,3,4].toSubspace
  discard aprioriMerge([s1,s2,s3,s4].toSet)

  for N in 2..5:
    let base = (1..N).toSubspace
    var subspaces = initSet[Subspace]()
    for s in base.lowDimProjections:
      subspaces.incl(s)
    #debug subspaces
    let merged = aprioriMerge(subspaces)
    assert(merged.len == 1)
    
  block:
    let N = 5
    let base = (1..5).toSubspace
    var subspaces = initSet[Subspace]()
    for s1 in base.lowDimProjections: # 4 dim
      for s2 in s1.lowDimProjections: # 3 dim
        for s3 in s2.lowDimProjections: # 2 dim
          subspaces.incl(s3)
    #debug subspaces, subspaces.len
    let merged = aprioriMerge(subspaces)
    assert(merged.len == 10) # a 5 dim subspace has 10 3-dim subspaces, created from the 10 2-dim subspaces

