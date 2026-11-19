# QHAT Developer's Guide

This guide is intended for developers who are contributing to the QHAT software.  This is mostly members of a LANL-internal team, but external merge requests will also be considered.

## Table of Contents

- [Development Policies](#development-policies)
- [Environment](#environment)
- [Qualtran Bloqs](#qualtran-bloqs)
- [Hamiltonian Terms](#hamiltonian-terms)

## Development Policies

This section covers some general policies on the development workflow for QHAT.

An important point to bear in mind up front: QHAT is research software that is intended to move quickly.  It will not have the same robustness and careful process that you would expect in a long-lived, high-reliability software project.  This may change in the future as the demands on QHAT change.

### Merge Requests

All code changes should be done through merge requests on GitHub.  This allows other developers to see what is in flight, helps us prevent breaking the `main` branch, and minimizes how much parallel efforts will step on each others’ toes.  A merge request does not need to be a complete effort; an intermediate state of a larger effort can be merged to `main` so long as it leaves `main` in a working state.

### CI Status

All merge requests must pass the full CI test suite before being merged to `main`.  We want to ensure that the `main` branch is always in a working state.

### Code Review

Merge requests contributed from outside of the core developer team and/or non-LANL developers must be reviewed by a core developer before merging.

Merge requests from core developers are not required to undergo code review.  This is left to the judgement of each core developer.  You are always welcome to ask for a review if you would like something to be checked over, but if you are confident in your work and it passes the CI pipeline you may merge without code review.

### Testing

All code should have tests where possible, and those tests should be added to the CI pipeline.  That lets us know when something is broken and needs to be fixed before being pushed to `main`.  This is a good use-case for agentic LLM coding tools (see [AI Usage](#ai-usage))

### AI Usage

QHAT makes heavy use of LLMs to support code development.  That spans the range from using LLMs to accelerate research up to agentic coding to implement well-understood features.  This goes back to the point above that QHAT is a fast-moving “research-grade” code rather than a high-reliability code.  However, ultimately the human who signs their name to the commits is responsible for the code changes.

## Environment

QHAT uses the pyLIQTR tools, which have very restrictive dependencies (either directly or transitively through its own dependencies).  Development is typically done in a Python 3.11 environment, but it appears that Python 3.14 may also work.  However QHAT does not guarantee that it will work on the latest and greatest version of Python because of the restrictions that propagate through the dependencies.  We hope to address this in the future.

## Qualtran Bloqs

The core tool that QHAT is built on is Qualtran and it’s “Bloq” framework.  A Bloq represents a piece of a quantum computing algorithm that can be composed together into a larger algorithm or used alone, and can be analyzed in a variety of ways.  Using the Bloq framework, QHAT builds an abstract representation of the algorithm by composing Bloqs into a tree of components that can be analyzed.  The basic functionality and interface for a Bloq is reasonably well-documented by Qualtran, and experience tells us that LLM tools are good at understanding and explaining the framework to assist QHAT developers.

### Protocols

QHAT leverages several analysis protocols supported by Bloqs:

- `t_complexity`: Walks the tree and aggregates counts of qubits, Clifford gates, and T gates to provide an estimate of the “cost” of running an algorithm on a quantum computer.
- `to_matrix`: Converts an algorithm to a unitary matrix, which enables numerical simulation, error analysis, etc.
  - We hope to add matrix-free operations as well, because dense-matrix operations do not scale well enough to analyze large algorithms.

Qualtran provides default mechanisms to do a variety of analyses, mostly through the `build_composite_bloq` method that recrusively decomposes a Bloq into its component sub-Bloqs until it reaches Bloqs that it knows how to analyze.  In order for a Bloq to be analyzed, it typically needs this method to be implemented.  In general this tends to be robust but memory-intensive and slow.  If you know something about the structure of your Bloq, there are a number of customization points that you can take advantage of to optimize analyses.  Sometimes even very simple optimizations can yield significant improvements in runtime and/or memory usage.

The most “complete” Bloq in QHAT is probably the `Trotterization` Bloq in `common/trotter_flattened.py`, which demonstrates a number of the optimizations available.

Some optimization points to be aware of

- `build_call_graph`: Instead of decomposing the tree, this method accumulates which Bloqs are called by the current Bloq, and how many times each is called, and returns Bloqs and counts.  For something like resource estimation, this is more efficient than `build_composite_bloq`.  But for something like computing the unitary matrix, this is not enough information to execute the analysis.
- `_t_complexity_`: If you know how to efficiently do resource estimation in terms of qubit count, Clifford gate count, and T gate count, implement this function and it will override the default decompose-and-analyze behavior.
- `tensor_contract`: Returns the unitary matrix equivalent to the algorithm represented by the Bloq, so this is the main entry point to optimize the `to_matrix` analysis protocol.
- `add_my_tensors`: Also used to combine the tensors of the sub-Bloqs, particularly when dealing with meta-operations (see below).  This is part of optimizing the `to_matrix` analysis protocol.

If you find that analysis of your block is very slow or consumes a significant amount of memory, finding the right hook to optimize a protocol is often the first step.

### Meta-Operations

There are a variety of what we call “meta-operations”: operations on operators.  This includes

- adjoint: $U^\dagger$
- controlled: $C[U]$
- power: $U^k$
- sign-controlled: a variant on controlled that may be more efficient in some situations

The default handling for these options tends to be robust but inefficient.

For adjoint, Qualtran wraps a Bloq in the adjoint meta-Bloq that reverses the order of the sub-Bloqs and applies the adjoint to each of them.  In some cases there may be a more efficient approach.  Consider, for example, the time-evolution operator of a Pauli string.  The adjoint can be computed simply by changing the sign of the rotation angle.

For control, Qualtran wraps a Bloq in the control meta-Bloq that applies a control to every sub-Bloq.  In some cases this may be very inefficient.  For example

- the time-evolution operator of a Pauli string: only the rotation needs to be controlled, not the basis change or parity operations;
- a Trotterization can share control information across all of the terms, which can be especially efficient if it is multi-controlled.

For both of these meta-Bloqs, not only are there the obvious optimizations, but applying a meta-Bloq forces the wrapped Bloq to decompose instead of using optimizations already in place.  For example, Trotterization typically involves many repetitions of the same pattern of sub-Bloqs so the analysis can be sped up by analyzing the repeating subset once and then applying that multiple times.  To do this for resource estimation, you would override the `_t_complexity_` protocol to understand and leverage that structure.  But if you need to analyze a controlled Trotterization or the adjoint of a Trotterization, the meta-Bloqs go around the `_t_complexity_` protocol and decompose the Trotterization so that the control or adjoint can be applied to each sub-Bloq.  In this case, you need to override the control or adjoint behavior to preserve the optimization.  You can overwrite the `adjoint` method to customize the adjoint behavior, and you can implement the `get_ctrl_system` method to customize the controlled behavior.  In the case of the `Trotterization` Bloq in QHAT, it doesn’t wrap the `Trotterization` in a meta-Bloq but instead returns a Bloq that represents controlled-Trotterization or the adjoint of Trotterization (both of which are simply the `Trotterization` block with different parameters).

For integers powers, Qualtran can simply apply the Bloq multiple times.  But you may be able to optimize this by implementing the `__pow__` method on your Bloq.  For example, the `Trotterization` Bloq simply takes more Trotter steps and runs for a longer evolution time.  And, thanks to the repetitive nature of Trotter steps, analyzing one `Trotterization` Bloq with more steps and a longer evolution time is less expensive than analyzing multiple `Trotterization` Bloqs.

TODO: sign-control (which may change in the code, so I’m leaving this for later)

## Hamiltonian Terms

QHAT represents Hamiltonians as an ordered collection of terms that are assumed to be summed together.  For product formulas, each term must be Hermitian and its time evolution operator must be exactly implementable.  For block-encoding, the requirements vary.  QHAT represents the terms using classes that all derive from the `BaseTerm` abstract base class and conform to the required interface:

- `num_qubits()` returns the number of qubits that the term operates on.  As a simplifying assumption, QHAT requires that each term present an interface as if it operates on all qubits in the Hamiltonian even if it only actually operates on a subset of those qubits.
- `generate_evolution_operator(time, hbar)` returns a Bloq that implements the time-evolution operator for that term.
- `__str__()` and `__repr__()` have their usual Python meaning and need to be implemented

Currently `BaseTerm` requires a `to_pauli_strings()` method, but this is expected to be a temporary patch until the term framework is completed.

TODO: sign-control (which may change in the code, so I’m laving this for later)

We expect that as the term framework is expanded and finalized, terms may also need to implement the following methods

- a method that returns the matrix equivalent to the term (not to the time-evolution operator of the term, but the term itself — used to compute the exact Hamiltonian in matrix form in certain analyses)
- some sort of “summability” operation (as a method or free function) that checks if two terms can be added together; for example, adding two Pauli string terms that have the same Pauli string can be achieved by simply combining them to a single Pauli string term and the new coefficient is the sum of the coefficients of the original Pauli string terms
- some sort of “commutation” operations, likely as free functions
  - one operation to compute the commutator of two terms
  - one operation to quickly check if two terms commute (default to checking if the commutator equals zero, but allows for more efficient implementations)
- a lexicographic sorting key (used in lexicographic ordering of Hamiltonian terms)
- a method that computes the spectral norm of the term (used in magnitude ordering of Hamiltonian terms, computing eigenvalue bounds, etc)