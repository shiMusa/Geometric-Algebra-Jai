println :: inline (msg: string, args: ..Any) {
    print(msg, ..args);
    print("\n");
}

println :: inline (arg: Any) {
    print("%\n", arg);
}

NewResizableArray :: (count: s64, $T: Type, $initialized := true, allocator: Allocator = .{}, alignment: s32 = -1) -> [..]T {
    arr := NewArray(count, T, initialized, allocator, alignment); // ! alloc
    res := resizable(arr);
    res.count = 0;
    return res;
}

powi :: (base: int, exponent: int) -> int {
    res := 1;
    for 0..exponent-1
        res *= base;
    return res;
}

hamming_weight :: (num: $T) -> int {
    count := 0;
    while num {
        num &= num-1;
        count += 1;
    }
    return count;
}