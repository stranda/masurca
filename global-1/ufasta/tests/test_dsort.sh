# Skip for now. sort and C++ std::sort do not agree on the order!
exit 77
NB=$(grep -c '^>' $TEST1)
for i in $(seq 0 $((NB - 1))); do
    echo $i
    ufasta dsort $TEST1 | ufasta extract -n "read$i" | tail -n +2 | sort -C
    EXPECT_SUCCESS "Sorted entry $i"
done
