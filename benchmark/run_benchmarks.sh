#!/bin/bash
# =====================================================================
# HENDOS against Starlab's kira: the whole benchmark in one go
#
#   cd HENDOS/benchmark
#   ./run_benchmarks.sh            the standard set (an hour or so)
#   ./run_benchmarks.sh quick      one small test, to try things out (minutes)
#   ./run_benchmarks.sh big        the standard set plus 16,000 stars (many hours)
#   ./run_benchmarks.sh check      only build, make the stars and run the check
#
# What it does:
#   1. builds hermite and makeics (make in ../run)
#   2. makes the initial conditions of every test (files in benchmark/ics)
#   3. runs  bench_kira.py --check  on the first of them: kira and HENDOS
#      are run for one output interval and the results compared. If that
#      fails the script stops, as nothing after it would mean anything.
#   4. runs bench_kira.py for every test, one directory each in
#      benchmark/results/<test>_s<seed>, then puts all the tables in
#      benchmark/results/summary.txt
#
# Settings, given in front of the command, e.g.
#   KIRA=/opt/starlab/bin/kira SEEDS="1 2 3" ./run_benchmarks.sh
#   KIRA      the kira executable                      [kira, on your path]
#   SEEDS     random number seeds: every test is run once for each
#             (single runs vary, so use several for anything you quote)
#                                                      [12345]
#   THREADS   OpenMP threads for HENDOS. 1 is the fair comparison with
#             kira, which is serial                    [1]
#   MAKE      0 = use the executables already in ../run, do not run make
#                                                      [1]
#   FORCE     1 = carry on even if the check fails     [0]
#
# THE TESTS (the table below; times in Myr). Change them here.
#   name             tend  HENDOS modes             what it tests
#   plummer_eq_1k     10   direct neighbours        the basic integrators
#   plummer_imf_1k    10   direct neighbours        the same with a mass spectrum
#   fractal_cool_1k    5   direct nofew nobin       binaries and multiples forming
#   plummer_eq_4k      4   direct neighbours        where the neighbour scheme pays
#   plummer_eq_16k     1   direct neighbours tree   (big only) large N, and the tree
# Each runs to a few crossing times; the energy errors are measured in the
# same way for both codes from the snapshots (see README.md).
# =====================================================================
set -u
cd "$(dirname "$0")"
BENCH=$(pwd)
RUN=$(cd ../run && pwd)
KIRA=${KIRA:-kira}
SEEDS=${SEEDS:-12345}
THREADS=${THREADS:-1}
MAKE=${MAKE:-1}
FORCE=${FORCE:-0}
SET=${1:-standard}

# >>> CHANGE HERE: the tests - name, tend (Myr), HENDOS modes
QUICK=(  "plummer_eq_1k   2  direct neighbours" )
STANDARD=( "plummer_eq_1k  10  direct neighbours"
           "plummer_imf_1k 10  direct neighbours"
           "fractal_cool_1k 5  direct nofew nobin"
           "plummer_eq_4k   4  direct neighbours" )
BIG=(      "plummer_eq_16k  1  direct neighbours tree" )

case "$SET" in
    quick)    TESTS=("${QUICK[@]}") ;;
    standard) TESTS=("${STANDARD[@]}") ;;
    big)      TESTS=("${STANDARD[@]}" "${BIG[@]}") ;;
    check)    TESTS=("${QUICK[@]}") ;;
    *) echo "usage: $0 [quick|standard|big|check]"; exit 1 ;;
esac

echo "=============================================================================="
echo "  HENDOS against kira: '$SET' set, seeds $SEEDS, HENDOS on $THREADS thread(s)"
echo "=============================================================================="

# ---- 1. build
if [ "$MAKE" = 1 ]; then
    echo "building hermite and makeics ..."
    if ! make -C "$RUN" hermite makeics > "$BENCH/build.log" 2>&1; then
        echo "ERROR: the build failed - see benchmark/build.log"
        echo "       (hermite needs a C compiler with OpenMP, makeics gfortran)"
        exit 1
    fi
fi
for exe in hermite makeics; do
    [ -x "$RUN/$exe" ] || { echo "ERROR: no $RUN/$exe"; exit 1; }
done
if ! command -v "$KIRA" > /dev/null 2>&1 && [ ! -x "$KIRA" ]; then
    echo "ERROR: cannot find kira ('$KIRA'). Put it on your path or give KIRA=/path/to/kira"
    exit 1
fi

# ---- 2. initial conditions: one file per test and seed
mkdir -p "$RUN/data" "$BENCH/results"
icname() { echo "data/bench_$1_s$2.dat"; }
FAILED=()
for t in "${TESTS[@]}"; do
    read -r name tend modes <<< "$t"
    for seed in $SEEDS; do
        ic=$(icname "$name" "$seed")
        [ -f "$RUN/$ic" ] && continue
        sed -e "s|^\( *seed *=\).*|\1 $seed|" -e "s|^\( *outfile *=\).*|\1 '$ic'|" \
            "$BENCH/ics/$name.nml" > "$BENCH/results/ics_${name}_s$seed.nml"
        if ! (cd "$RUN" && ./makeics "$BENCH/results/ics_${name}_s$seed.nml" \
                > "$BENCH/results/ics_${name}_s$seed.log" 2>&1) || [ ! -f "$RUN/$ic" ]; then
            # (a fractal can fail for some seeds - see FRACTAL_README.md)
            echo "  ** could not make $ic - see results/ics_${name}_s$seed.log; test skipped"
            FAILED+=("${name}_s$seed")
        else
            echo "  made $ic"
        fi
    done
done

# ---- 3. the check
read -r name tend modes <<< "${TESTS[0]}"
seed=${SEEDS%% *}
echo
echo "checking kira and HENDOS on $(icname "$name" "$seed") ..."
(cd "$RUN" && python3 "$BENCH/bench_kira.py" "$(icname "$name" "$seed")" --tend "$tend" \
    --kira "$KIRA" --threads "$THREADS" --out "$BENCH/results/check" --check)
ok=$?
if [ $ok != 0 ] && [ "$FORCE" != 1 ]; then
    echo "The check failed: stopping. Look at the logs in benchmark/results/check"
    echo "(FORCE=1 carries on regardless)."
    exit 1
fi
[ "$SET" = check ] && exit $ok

# ---- 4. the benchmarks
for t in "${TESTS[@]}"; do
    read -r name tend modes <<< "$t"
    for seed in $SEEDS; do
        [[ " ${FAILED[*]-} " == *" ${name}_s$seed "* ]] && continue
        out="$BENCH/results/${name}_s$seed"
        echo
        echo "------------------------------------------------------------------------------"
        echo "  $name, seed $seed, to $tend Myr, HENDOS modes: $modes"
        echo "------------------------------------------------------------------------------"
        (cd "$RUN" && python3 "$BENCH/bench_kira.py" "$(icname "$name" "$seed")" \
            --tend "$tend" --modes $modes --kira "$KIRA" --threads "$THREADS" --out "$out")
    done
done

# ---- summary: every table in one file
summary="$BENCH/results/summary.txt"
{
    echo "HENDOS against kira - '$SET' set, seeds $SEEDS, HENDOS on $THREADS thread(s)"
    echo "kira: $(command -v "$KIRA" || echo "$KIRA")    date: $(date '+%Y-%m-%d %H:%M')"
    for d in "$BENCH"/results/*_s*/; do
        [ -f "$d/results.txt" ] && { echo; cat "$d/results.txt"; }
    done
} > "$summary"
echo
echo "Done. All the tables are in benchmark/results/summary.txt, and each test's"
echo "plots (bench_time_error.png, bench_halfmass.png) are in its own directory there."
