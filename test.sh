#!/bin/bash

# Check the required files exist
testFiles=("brock200_1.clq" "brock200_2.clq" "brock200_3.clq" "brock200_4.clq")
for f in ${testFiles[*]}; do
    if [[ ! -f ${f} ]]; then
        echo "Could not find required test file ${f}"
        exit -1
    fi
done

fail=0

# -- Maximising Tests -- #
res=$(./maxclique -f brock200_1.clq --hpx:threads=1 --hpx:queuing=local)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 21 ]]; then
    echo "FAILED Maxclique brock200_1"
    fail=$((fail+1))
fi

res=$(./maxclique -f brock200_2.clq --hpx:threads=1 --hpx:queuing=local)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 12 ]]; then
    echo "FAILED Maxclique brock200_2"
    fail=$((fail+1))
fi

res=$(./maxclique -f brock200_3.clq --hpx:threads=1 --hpx:queuing=local)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 15 ]]; then
    echo "FAILED Maxclique brock200_3"
    fail=$((fail+1))
fi

res=$(./maxclique -f brock200_4.clq --hpx:threads=1 --hpx:queuing=local)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 17 ]]; then
    echo "FAILED Maxclique brock200_4"
    fail=$((fail+1))
fi

res=$(./maxclique -f brock200_1.clq --hpx:threads=all --hpx:queuing=local)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 21 ]]; then
    echo "FAILED Maxclique brock200_1"
    fail=$((fail+1))
fi

res=$(./maxclique -f brock200_2.clq --hpx:threads=all --hpx:queuing=local)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 12 ]]; then
    echo "FAILED Maxclique brock200_2"
    fail=$((fail+1))
fi

res=$(./maxclique -f brock200_3.clq --hpx:threads=all --hpx:queuing=local)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 15 ]]; then
    echo "FAILED Maxclique brock200_3"
    fail=$((fail+1))
fi

res=$(./maxclique -f brock200_4.clq --hpx:threads=all --hpx:queuing=local)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 17 ]]; then
    echo "FAILED Maxclique brock200_4"
    fail=$((fail+1))
fi

# -- Find Tests -- #
res=$(./maxclique --find-clique -s 22 -f brock200_1.clq --hpx:threads=1 --hpx:queuing=local --hpx:ini=hpx.stacks.small_size=0x40000)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 0 ]]; then
    echo "FAILED Find No Exist Maxclique brock200_1"
    fail=$((fail+1))
fi

res=$(./maxclique --find-clique -s 13 -f brock200_2.clq --hpx:threads=1 --hpx:queuing=local --hpx:ini=hpx.stacks.small_size=0x40000)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 0 ]]; then
    echo "FAILED Find No Exist Maxclique brock200_2"
    fail=$((fail+1))
fi

res=$(./maxclique --find-clique -s 26 -f brock200_3.clq --hpx:threads=1 --hpx:queuing=local --hpx:ini=hpx.stacks.small_size=0x40000)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 0 ]]; then
    echo "FAILED Find No Exist Maxclique brock200_3"
    fail=$((fail+1))
fi

res=$(./maxclique --find-clique -s 18 -f brock200_4.clq --hpx:threads=1 --hpx:queuing=local --hpx:ini=hpx.stacks.small_size=0x40000)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 0 ]]; then
    echo "FAILED Find No Exist Maxclique brock200_4"
    fail=$((fail+1))
fi

res=$(./maxclique --find-clique -s 18 -f brock200_1.clq --hpx:threads=all --hpx:queuing=local --hpx:ini=hpx.stacks.small_size=0x40000)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 18 ]]; then
    echo "FAILED Find Exist Maxclique brock200_1"
    fail=$((fail+1))
fi

res=$(./maxclique --find-clique -s 11 -f brock200_2.clq --hpx:threads=all --hpx:queuing=local --hpx:ini=hpx.stacks.small_size=0x40000)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 11 ]]; then
    echo "FAILED Find Exist Maxclique brock200_2"
    fail=$((fail+1))
fi

res=$(./maxclique --find-clique -s 12 -f brock200_3.clq --hpx:threads=all --hpx:queuing=local --hpx:ini=hpx.stacks.small_size=0x40000)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 12 ]]; then
    echo "FAILED Find Exist Maxclique brock200_3"
    fail=$((fail+1))
fi

res=$(./maxclique --find-clique -s 16 -f brock200_4.clq --hpx:threads=all --hpx:queuing=local --hpx:ini=hpx.stacks.small_size=0x40000)
size=$(echo "${res}" | awk '/Size:/ {print $2}')
if [[ $size != 16 ]]; then
    echo "FAILED Find Exist Maxclique brock200_4"
    fail=$((fail+1))
fi

if [[ fail -gt 0 ]]; then
    echo "${fail} tests failed!"
else
    echo "Tests Succeeded"
fi
