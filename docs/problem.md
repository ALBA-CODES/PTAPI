# The problem class

Everything that depends on the problem is in two types written by the user:
the solution struct and a class derived from `Problem<S>`.

```cpp
struct mySolution: public solution{ /* the encoding */ };

class MyProblem: public Problem<mySolution>{
public:
    mySolution construction();                 // required
    mySolution neighbor(mySolution sol);       // required
    double evaluate(const mySolution& sol);    // one of three forms
    // optional: neighborEval, constructionReplica, moves with probabilities
};
```

The functions of this class are called from several threads at the same time
(one replica per thread), each call with its own solution. Members that are
only read (the instance data) need no protection; anything the functions write
outside the solution they receive does.

## The solution

`solution` has one field that matters: `double evalSol`, the cost. The API
minimizes it; to maximize, store the negated value. `Nup` and `Ndown` exist
only so that code written for version 1 compiles.

The struct is copied when a replica records a new best solution and moved or
swapped elsewhere, so it must be copyable. Fields that speed up the evaluation
(partial sums, positions) can live in it: see `diff` in the quick start.

## `construction`

```cpp
virtual S construction() = 0;
```

Returns an initial solution **already evaluated** (`evalSol` filled in). It is
called once for each replica.

```cpp
virtual S constructionReplica(int index, int count);   // optional; default: construction()
```

Override it to give each replica its own constructive heuristic: `index` goes
from 0 (first temperature) to `count - 1` (last). A usual choice is a greedy
solution for the cold replicas and random ones for the hot replicas. With
`PT::setParallelConstruction(true)` the replicas are built by several threads
at the same time, which pays off when the heuristic is expensive. Example:
[`exampleInitialSolutions.cpp`](../example/functions/exampleInitialSolutions.cpp).

## `neighbor`

```cpp
virtual S neighbor(S sol) = 0;
```

Receives a copy of the current solution of a replica and returns a neighbor.
It does not need to evaluate it. Apply the move to the argument and return it,
without a second copy.

## `evaluate`

Declare it in any of these forms; the API detects which one at compile time.

```cpp
S      evaluate(S sol);          // may also change the solution (e.g. fill auxiliary fields)
double evaluate(S sol);          // form of the examples published with version 1
double evaluate(const S& sol);   // no copy of the solution: the fastest of the three
```

`evaluate` is not needed at all if the class overrides `neighborEval` and no
initial solutions are given to `PT`. Example:
[`exampleEvaluate.cpp`](../example/functions/exampleEvaluate.cpp).

## `neighborEval`: a step without copies and with delta evaluation

```cpp
virtual void neighborEval(const S& current, S& candidate);   // optional
```

One step of the Markov chain: write into `candidate` an **evaluated** neighbor
of `current`. The default does `candidate = evaluate(neighbor(current))`.
Override it for speed:

- reuse the memory of `candidate` (`candidate.sol = current.sol` copies the
  values without allocating, because `candidate` already holds a solution of
  the same size; its content on entry is unspecified);
- compute `candidate.evalSol` from `current.evalSol` and the change the move
  causes, instead of evaluating the whole solution.

```cpp
void neighborEval(const solTSP& cur, solTSP& cand){
    int a, b;                                    // 2-opt: reverse the tour between a and b
    do{ a = ptapi::randInt(0, n-1); b = ptapi::randInt(0, n-1); }while(a == b);
    if(a > b) std::swap(a, b);
    cand.sol = cur.sol;
    std::reverse(cand.sol.begin()+a, cand.sol.begin()+b+1);
    if(a == 0 && b == n-1){ cand.evalSol = cur.evalSol; return; }
    int pa = cur.sol[(a+n-1)%n], ca = cur.sol[a], cb = cur.sol[b], nb = cur.sol[(b+1)%n];
    cand.evalSol = cur.evalSol - dist(pa,ca) - dist(cb,nb) + dist(pa,cb) + dist(ca,nb);
}
```

If both `neighbor`/`evaluate` and `neighborEval` exist, the chain uses
`neighborEval`; keep the two consistent (the tests of the API compare them).

## Random numbers

```cpp
ptapi::rng()          // random engine (std::mt19937_64) of the replica running on this thread
ptapi::randInt(a, b)  // uniform integer in [a, b]
ptapi::rand01()       // uniform real in [0, 1)
```

Use them inside `construction`, `neighbor`, `evaluate` and `neighborEval`
instead of creating a `std::random_device` and an engine on every call, which
is slow and makes the run impossible to repeat. A run in which the problem
only uses these functions and that stops by PTL is reproducible with
`PT::setSeed()`, whatever the number of threads.

## Several moves with probabilities

```cpp
MyProblem(){ setMoveProbabilities({0.6, 0.3, 0.1}); }   // weights; they are normalized

mySolution neighbor(mySolution s){
    switch(chooseMove()){
        case 0: /* swap */    break;
        case 1: /* insert */  break;
        case 2: /* 2-opt */   break;
    }
    return s;
}
```

`chooseMove()` draws the index of a move. After the run,
`PT::getMoveStatistics()` gives, for each move, how many times it was proposed
and accepted, which shows the moves that are worth their probability.
`getMoveCount()` and `getMoveProbabilities()` return what was registered.
Example: [`exampleMoves.cpp`](../example/functions/exampleMoves.cpp).

## Errors

An exception thrown by the problem code inside a worker thread stops the run
and is thrown again by `start()` in the calling thread
(`getStopReason()` is then `PT_STOP_ERROR`).
