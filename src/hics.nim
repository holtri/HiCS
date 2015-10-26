
import dataset
import subspace
import slicing
import stattest
import preprocessing
import utils
import math
import future
import sets
import topk
import tables
import sequtils
import algorithm

type
  IndexMap = seq[int] # not nil

  Parameters* = object
    numIterations*: int
    alpha*: float
    numCandidates*: int
    maxOutputSpaces*: int


proc initParameters*(
    numIterations = 100,
    alpha = 0.1,
    numCandidates = 500,
    maxOutputSpaces = 1000,
  ): Parameters =
  Parameters(
    numIterations: numIterations,
    alpha: alpha,
    numCandidates: numCandidates,
    maxOutputSpaces: maxOutputSpaces
  )

proc computeAverageDeviation*[T](subspace: Subspace, ds: Dataset, preproData: PreproData, params: Parameters, statTest: T, cmpAttr: int): float =
  let D = subspace.len
  let N = ds.nrows
  let M = (pow(params.alpha, 1 / (D-1)) * N.toFloat).toInt

  var iselAll = newIndexSelection(N)
  var iselCur = newIndexSelection(N)

  var totalDeviation: float = 0
  var numIterations: int = (params.numIterations/subspace.len).int

  for iter in 0 .. < numIterations:
    iselAll.reset(true)

    for j in subspace.asSeq:
      if j != cmpAttr:
        iselCur.selectRandomBlock(M)

        for indRank, used in iselCur:
          if not used:
            let indObject = preproData.indexMaps[j][indRank]
            iselAll[indObject] = false

    let deviation: float = statTest.computeDeviation(ds, cmpAttr, iselAll)

    totalDeviation += deviation
  return totalDeviation / numIterations.toFloat

proc computeAverageDeviation*[T](subspace: Subspace, ds: Dataset, preproData: PreproData, params: Parameters, statTest: T): Table[int,float] =

  let D = subspace.len
  let N = ds.nrows
  let M = (pow(params.alpha, 1 / (D-1)) * N.toFloat).toInt

  var iselAll = newIndexSelection(N)
  var iselCur = newIndexSelection(N)

  var deviations = initTable[int,float]()

  for s in subspace:

    var cmpAttr = s
    var totalDeviation: float = 0
    var numIterations: int = params.numIterations #(params.numIterations/subspace.len).int
    for iter in 0 .. < numIterations:
      iselAll.reset(true)

      for j in subspace.asSeq:
        if j != cmpAttr:
          iselCur.selectRandomBlock(M)

          for indRank, used in iselCur:
            if not used:
              let indObject = preproData.indexMaps[j][indRank]
              iselAll[indObject] = false

      let deviation: float = statTest.computeDeviation(ds, cmpAttr, iselAll)

      totalDeviation += deviation

    deviations[cmpAttr] = totalDeviation / numIterations.toFloat
  return deviations

proc greedyDeviation*[T](ds: Dataset, params: Parameters, focusDim: int, statTest: T, preproData: PreproData): StoreTopK[(float, Subspace)] =

  let D = ds.ncols
  var outputSpaces = newTupleStoreTopK[float,Subspace](params.maxOutputSpaces, keepLarge=true)

  type startElement = tuple[deviation: float, referenceDim: int]
  var start2Dim: seq[startElement] = @[]

  for i in 1..D-1:
    if i == focusDim:
      continue
    let deviation = computeAverageDeviation([i,focusDim].toSubspace, ds, preproData, params, statTest, focusDim)
    start2Dim.add((deviation, i))
  start2Dim = start2Dim.sortedByIt(it.deviation).reversed
  #debug start2Dim

  var maxDeviation: float = -Inf
  var maxSubspace: seq[int] = @[focusDim]

  for s in start2Dim:
    var tmp = maxSubspace
    tmp.add(s.referenceDim)
    let deviation = computeAverageDeviation(tmp.toSubspace, ds, preproData, params, statTest, focusDim)
    #debug tmp, deviation
    if deviation > maxDeviation:
      maxDeviation = deviation
      maxSubspace = tmp
  outputSpaces.add((maxDeviation, maxSubspace.toSubspace))
  maxSubspace.sort(system.cmp)
  debug maxSubspace , maxDeviation
  result = outputSpaces

proc deviation2DimSpaces*(ds: Dataset, params: Parameters, focusDim: int, verbose = false): StoreTopK[(float, Subspace)] =
  let N = ds.nrows
  let D = ds.ncols

  let preproData = ds.generatePreprocessingData()
  let statTest = initKSTest(ds, preproData, (params.alpha * N).toInt, verbose)

  var outputSpaces = newTupleStoreTopK[float,Subspace](params.maxOutputSpaces, keepLarge=true)
  var spaces = generate2DSubspaces(D)
  for s in spaces:
    var t = s
    if not s.contains(focusDim):
        continue

    let deviation = computeAverageDeviation(s, ds, preproData, params, statTest, focusDim)
    t.excl(focusDim)
    outputSpaces.add((deviation, t))

  result = outputSpaces

proc hicsFramework*(ds: Dataset, params: Parameters, focusDim: int, verbose = false): StoreTopK[(float, Subspace)] =

  let N = ds.nrows
  let D = ds.ncols

  let preproData = ds.generatePreprocessingData()
  let statTest = initKSTest(ds, preproData, (params.alpha * N).toInt, verbose)

  var outputSpaces = newTupleStoreTopK[float,Subspace](params.maxOutputSpaces, keepLarge=true)

  var d = 2

  var spaces = generate2DSubspaces(D)
  #var focusDim = 2
  #var spaces = generate2DSubspaces(D,focusDim)

  while spaces.len > 0:
    if verbose: echo ifmt" * processing subspaces of dim $d [number of spaces: ${spaces.len}]"

    var spacesForAprioriMerge = newTupleStoreTopK[float,Subspace](params.numCandidates)
    #debug spaces
    for s in spaces:
      #for 2 dim spaces
      var t = s
      if not s.contains(focusDim):
        continue

      let deviation = computeAverageDeviation(s, ds, preproData, params, statTest, focusDim)
      if d == 2:

        t.excl(focusDim)
        debug focusDim, t, deviation

      if verbose: debug s, deviation
      spacesForAprioriMerge.add((deviation, s))
      outputSpaces.add((deviation, s))

    var candidateSet = newSubspaceSet()
    for deviation, subspace in spacesForAprioriMerge.items:
      candidateSet.incl(subspace)
    assert candidateSet.len == spacesForAprioriMerge.size

    let oldNumSpaces = spaces.len
    spaces = candidateSet.aprioriMerge(focusDim)
    for s in spaces:
      assert s.dimensionality == d+1


    if verbose:
      let info = if oldNumSpaces <= params.numCandidates:
        ifmt"using all spaces; number of $d-dim spaces: $oldNumSpaces <= numCandidates: ${params.numCandidates}"
        else:
          ifmt"number $d-dim spaces [$oldNumSpaces] exceeds numCandidates [${params.numCandidates}] => limiting to top ${params.numCandidates}]"
      echo ifmt" => number of merged 3-dim spaces: ${spaces.len} ($info)"
    inc d

  #for contrast, subspace in outputSpaces.sortedItems:
  #  echo contrast, subspace

  result = outputSpaces

