import os
import parseutils
import strutils
import sequtils
import sets
import utils
import dataset
import preprocessing
import testdata
import hics
import slicing
import subspace
import stattest as module_stattest
import typetraits
import topk
import threadpool
import locks
import math
#import optionals


when false:
  #let ds = loadDataset("/tmp/data.csv")
  let ds = testData2DLinear(1000)
  let preproData = ds.generatePreprocessingData()
  let params = initParameters(numIterations = 100, alpha = 0.05)
  let statTest = initKSTest(ds, preproData, (params.alpha * ds.nrows).toInt)
  let contrast = computeContrast([0,1].toSubspace, ds, preproData, params, statTest)
  debug contrast

const
  version = "0.9.0"

  usage = """
Usage:
  HiCS [options] --csvIn <FILE>

Options:
  --csvIn <FILE>                  File name of CSV input data (requires delimiter: ';')
  --csvOut <FILE>                 File name of subspace mining output data
  --hasHeader                     Specifies if input CSV has a header row
  --numRuns <INT>                 Number of Monte Carlo iterations to perform (default: 100)
  --numCandidates <INT>           Number of candidates in subspace beam search (default: 500)
  --maxOutputSpaces <INT>         Maximum number of subspaces in the final output (default: 1000)
  --alpha <DOUBLE>                Size of the test statistic specified as a fraction (default 0.1)
  --onlySubspace <INT,INT,...>]   Using this command will not perform a subspace search.
                                  Instead, HiCS will only compute the contrast of the specified
                                  subspace. Note: The list of INTs must not contain whitespace and
                                  indexing starts at zero. For example: 0,3,12,51
  -h, --help                      Print this help message.
  -v, --version                   Print version information.
  -s, --silent                    Disables debug output on stdout.
"""

# ----------------------------------------
# Misc helpers
# ----------------------------------------

proc storeResults(filename: string, results: StoreTopK[(float, Subspace)]) =
  var o = open(filename, fmWrite)
  for contrast, subspace in results.sortedItems:
    let subspaceSorted = subspace.asSeq.mapIt(string, $it).join(";")
    o.writeln(ifmt"${contrast}%20.12f; ${subspaceSorted}")
  o.close


# ----------------------------------------
# Command line parsing
# ----------------------------------------

proc parse*[T](s: string): tuple[okay: bool, value: T] =
  ## semi-generic parsing, TODO: use Option[T]
  when T is int:
    var x: int
    if parseInt(s, x) > 0:
      return (true, x)
    else:
      return (false, 0)
  elif T is float:
    var x: float
    if parseFloat(s, x) > 0:
      return (true, x)
    else:
      return (false, 0.0)
  elif T is string:
    return (true, s)
  elif T is Subspace:
    let fields = s.split(",")
    debug fields
    var s = [].toSubspace
    for f in fields:
      let (b,i) = parse[int](f)
      if b: # should we check if dimensions are ok? 0 <= i < dimOfFile
        s.incl(i)
      else:
        return (false, s)
    return (true, s)


proc echoUsage(andQuit: bool = true) =
  echo usage
  if andQuit:
    quit 0


# must be defined here, the following procs are closures...
let args = commandLineParams()

proc hasArg(checkEither: varargs[string]): bool =
  for s in checkEither:
    if s in args:
      return true
  return false

proc pairExtractor[T](name: string): T =
  for i in 0 ..< args.len:
    if args[i] == name:
      if i == args.len-1:
        quit ifmt"Error: The parameter '$name' must be followed by a value of type ${T.name}."
      let parse = parse[T](args[i+1])
      if parse.okay:
        return parse.value
      else:
        quit ifmt"Error: Value for parameter '$name' cannot be cast to ${T.name} (value given was: '${args[i+1]}')."
  quit ifmt"Error: The parameter '$name' was not specified but is mandatory."

proc pairExtractor[T](name: string, default: T): T =
  for i in 0 ..< args.len:
    if args[i] == name:
      if i == args.len-1:
        quit ifmt"Error: The parameter '$name' must be followed by a value of type ${T.name}."
      let parse = parse[T](args[i+1])
      if parse.okay:
        return parse.value
      else:
        quit ifmt"Error: Value for parameter '$name' cannot be cast to ${T.name} (value given was: '${args[i+1]}')."
  return default


# "Usage java -jar hics.jar [--numRuns <INT>] [--numCandidates <INT>] [--alpha <DOUBLE>] [--hasHeader] [--onlySubspace <INT,INT,...>] [--csvOut <FILE>] --csvIn <FILE>")
# "Note: The list of INTs must not contain whitespace and indexing starts at zero, i.e.: 0,3,12,51")



if args.len == 0:
  echoUsage()

if hasArg("-h", "--help"):
  echoUsage()

if hasArg("-v", "--version"):
  echo "HiCS version ", version
  quit 0

let silent = hasArg("-s", "--silent")
let hasHeader = hasArg("--hasHeader")

let numRuns         = pairExtractor("--numRuns", 100)
let numCandidates   = pairExtractor("--numCandidates", 500)
let maxOutputSpaces = pairExtractor("--numCandidates", 1000)
let alpha           = pairExtractor[float]("--alpha", 0.1)

let fileI           = pairExtractor[string]("--csvIn")
let fileO           = pairExtractor("--csvOut", "out.csv")

let onlySubspace    = pairExtractor("--onlySubspace", [].toSubspace)
let allSubspaces    = pairExtractor("--allSubspaces", 1)
let startDimension  = pairExtractor("--startDimension", 0)
let maxSubSpaceDim  = pairExtractor("--maxSubspaceDim", 0)

let randomSeed = pairExtractor("--randomSeed", -1)

if not silent:
  echo ifmt"Running with parameters:"
  echo ifmt"  csvIn         = $fileI"
  echo ifmt"  csvOut        = $fileO"
  echo ifmt"  hasHeader     = $hasHeader"
  echo ifmt"  numRuns       = $numRuns"
  echo ifmt"  numCandidates = $numCandidates"
  echo ifmt"  alpha         = $alpha"
  echo ifmt"  randomSeed    = $randomSeed"


let params = initParameters(
  numIterations = numRuns,
  alpha = alpha,
  numCandidates = numCandidates,
  maxOutputSpaces = maxOutputSpaces
)

#initialize random numbers
if randomSeed > 0:
  randomize(randomSeed)
else:
  randomize()

let ds = loadDataset(fileI, hasHeader, silent)
if not silent:
  echo ifmt"Data dimensions: ${ds.nrows} x ${$ds.ncols}"


proc storeBruteForceResults(filename: string, output: string, index: int) =
  var outFile = ifmt("$filename_$index.csv")
  echo outFile
  var o = open(outFile, fmWrite)
  o.writeln(output)
  o.close

proc expandSubspace (currentSeq: seq[int], maxSubspaceDim: int, maxFullSpaceDim: int, preproData: PreproData, statTest: KSTest, ds: Dataset): string =
  var resultString: string = ""

  if(len(currentSeq) < maxSubspaceDim):
    for i in currentSeq[high(currentSeq)]+1..maxFullSpaceDim:
      let tmp = expandSubspace(currentSeq & (i), maxSubspaceDim, maxFullSpaceDim, preproData, statTest, ds)
      resultString = resultString & tmp

  if(len(currentSeq) > 1):
    let contrast = computeContrast(currentSeq.toSubspace, ds, preproData, params, statTest)
    let subspaceString = foldr(currentSeq.mapIt(string,";" & $it & " "), a & b)
    resultString = ifmt("\n$contrast $subspacestring") & resultString

  return resultString

if onlySubspace.dimensionality == 0:

  if allSubspaces > 1:
    let preproData = ds.generatePreprocessingData()
    let statTest = initKSTest(ds, preproData, (params.alpha * ds.nrows).toInt, verbose=not silent)
    let maxFullSpaceDim = ds.nCols - 1

    debug maxFullSpaceDim, maxSubspaceDim, startDimension
    let initialSeq = @[startDimension]
    var res = expandSubspace(initialSeq, maxSubspaceDim, maxFullSpaceDim, preproData, statTest, ds)
    storeBruteForceResults(fileO,res, initialSeq[0])
    sync()

  else:
    # regular mode:
    let results = hicsFramework(ds, params, verbose = not silent)

    storeResults(fileO, results)

else:

  let preproData = ds.generatePreprocessingData()
  let statTest = initKSTest(ds, preproData, (params.alpha * ds.nrows).toInt, verbose=not silent)
  let contrast = computeContrast(onlySubspace, ds, preproData, params, statTest)

  #echo ifmt"Runtime with preprocessing: ${(t3-t1).toDouble / 1000}%.3f"
  #echo ifmt"Runtime contrast calculation only: ${(t3-t2).toDouble / 1000}%.3f"
  echo ifmt"Contrast of subspace $onlySubspace"
  echo contrast






